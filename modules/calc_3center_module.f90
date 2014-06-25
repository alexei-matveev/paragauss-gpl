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
  module calc_3center_module
  ! Purpose : the calculation of 3-center integrals (nuclear,
  !           coulomb) for [na,la,nb,lb] qadrupel using factorized
  !           expressions.
  !           (2-center also must be here !)
  ! The angular factors for nuclear and coulomb integrals are
  ! calculated SEPARATELY and should be joined later to use
  ! particular case of Ang(Ma,Mb,Mc=1) for nuclear integrals !!!
  ! The order could be:
  ! ==================
  ! - all Ang(Ma,Mb,Mc)%J
  ! - R_nuc%J
  ! - R_nuc%J * Ang(Ma,Mb,1)%J -> nuclear attraction
  ! - symmetrize Ang over (c,Mc)
  ! - R_fit%J
  ! - R_fit%J * Ang(Ma,Mb,t)
!================================================================
! End of public interface of module -- dont believe that, it`s a bad joke!
!================================================================
#include "def.h"
    use type_module
    use datatype, only: arrmat5
    use unique_atom_module
    use solid_harmonics_module, only : solid_harmonics_scalar
    use gamma_module
    use int_data_2cob3c_module
    use options_module, only: options_integral_expmax
    use integralpar_module
    use iounitadmin_module
    use potential_module,N_points_pot=>N_points
!    use fitcontract_module
#ifdef WITH_EPE
    use ewaldpc_module
#endif
    use timer_module
    use time_module
    use solv_cavity_module, only: to_calc_grads,with_pc,fixed_pc !!!!!!!!!!!
    use gradient_data_module, only: gradient_index, gradient_data_n_gradients, &
                                  gradient_data_n_spin_gradients
    use ll_calculate_grads_module, only: grad_mat &
                                       , grad_mat_p &
                                       , dervs_mat &
                                       , nuc_grad_gr &
                                       , ima &
                                       , imb &
                                       , imc &
                                       , moving_a &
                                       , moving_b &
                                       , moving_c &
                                       , do_rotation &
                                       , rotmat &
                                       , k_gr_0 &
                                       , k_gr_1 &
                                       , solv_grad_gr &
                                       , ca_dervs_mat &
                                       , coul_int_grad &
                                       , coul_int_grad_totsym &
                                       , coul_int_dervs &
                                       , coul_int_dervs_totsym &
                                       , coul_int_ca_dervs &
                                       , help_arr_gr1 &
                                       , pointer_prim_int &
                                       , add_nuc_or_pvsp



    use calc3c_switches
    use pointcharge_module !!!!!!!!!!!!!!
    use elec_static_field_module, only: totsym_field_length,surf_points_grad_index !!!!!
    use error_module, only:MyID

    USE_MEMLOG

    implicit none
    save
    private

    !
    ! FIXME: these are re-exports form other modules:
    !
    public :: prim_int_3cob_grad
    public :: prim_int_coul_dervs
    public :: prim_int_3cob_coul_grad
    public :: prim_int_2cob_nuc

    public :: unique_atoms
    public :: unique_atom_grad_info
    public :: n_unique_atoms

    public :: ewpcdervs
#ifdef WITH_EPE
    public :: ewpc_n
#endif

    public :: i4_kind, r8_kind

    public :: write_to_output_units

    public :: start_timer
    public :: stop_timer

    public :: integralpar_dervs
    public :: integralpar_gradients
    public :: integralpar_3c_co
    public :: integralpar_3cob_grad
    public :: integralpar_cpksdervs
    public :: integralpar_pot_for_secderiv

    public :: options_integral_expmax

#ifdef no_cpks_coul_grads
    real(kind=r8_kind),pointer :: cpks_grad_mat(:,:)
    real(kind=r8_kind),pointer :: cpks_grad_mat_p(:)
#else
    real(kind=r8_kind),pointer :: cpks_grad_mat(:,:,:,:,:,:)
    real(kind=r8_kind),pointer :: cpks_grad_mat_p(:,:,:,:,:)
#endif

  real(kind=r8_kind),     allocatable  :: rsh(:) !to be calculated once in CSH_scalar
  real(kind=r8_kind),     allocatable  :: rshg(:,:) !to be calculated once in CSHG_scalar
  real(kind=r8_kind),     allocatable  :: rshgg(:,:,:) !to be calculated once in CSHG_scalar

  integer(kind=i4_kind),public:: equ_c,uni_c

  logical, allocatable:: do_rotation_eq(:)
  logical::check_ab,check_bc
  logical:: new_nuc=.false.

  real(kind=r8_kind), allocatable :: two_fact2_vec(:,:),two_aexp_arr(:),two_bexp_arr(:)

!     temporary for Ang_radial procedures
      real(kind=r8_kind), allocatable :: grad_totsymM(:,:,:,:),   gradM(:,:,:,:)
      real(kind=r8_kind), allocatable, dimension(:,:,:,:,:) ::&
                               dervs_totsymM_temp,dervs_totsymM,dervsM,ca_dervsM

      real(kind=r8_kind), allocatable :: Gt(:,:,:),Ga(:,:),Gb(:,:),GtFvec(:,:,:), &
                                         GGt(:,:,:),GGa(:,:,:),GGb(:,:,:),GaGb(:,:,:), &
                                         GtGa(:,:,:,:),GtGb(:,:,:,:) ! sol harm Gt derivative

      real(kind=r8_kind), allocatable :: radang_temp(:,:,:),radangF(:), &
                                         radang_temp_ga(:,:,:,:),radang_temp_gb(:,:,:,:), &
                                         radang_temp_gx(:,:,:,:)
      real(kind=r8_kind), allocatable,dimension(:,:,:,:,:) :: &
                          radang_temp_gga,radang_temp_ggb,radang_temp_gagb
      real(kind=r8_kind), allocatable :: c_exp_arg2(:),g_shift_fac(:)


  real(kind=r8_kind)::dval
  real(kind=r8_kind),  allocatable  :: temp_overlap(:),gamma_arg2_vec(:,:)
  real(kind=r8_kind),  allocatable  :: temp_overlap_nuc(:)
  type(arrmat5), allocatable, target :: coul_int(:) ! (-1:lmax_ch)
  real(kind=r8_kind),pointer,dimension(:,:,:,:,:) :: pointer_coul, &
                     pointer_coul_a,pointer_coul_b,pointer_coul_c

    complex(kind=c16_kind), allocatable   :: A_factor(:,:,:), CSH(:,:,:)
    complex(kind=c16_kind), allocatable, dimension(:,:,:,:) :: A_factorGa,CSHG,&
                                                               A_factorGb,CSHGb
    complex(kind=c16_kind), allocatable, dimension(:,:,:,:,:) :: &
                     CSHGGa,CSHGGb,CSHGaGb,A_factorGGa,A_factorGGb,A_factorGaGb

    complex(kind=c16_kind), parameter     :: czero = (0.0_r8_kind,0.0_r8_kind)
    real(kind=r8_kind), public            :: xa(3), xb(3), xc(3), xd(3)
    real(kind=r8_kind), pointer, public   :: aexps(:), bexps(:), cexps(:)
    real(kind=r8_kind),     allocatable   :: fact0(:),fact1(:),fact2(:),fact3(:),fact4(:)
    real(kind=r8_kind),allocatable,public :: fact0_arr(:,:), fact1_arr(:,:), fact2_arr(:,:)
    real(kind=r8_kind),     allocatable   :: aexp_arr(:), bexp_arr(:)
    real(kind=r8_kind),     allocatable   :: gamma_arg(:,:), gamma_arg2(:), &
                                             gamma_help(:,:) ! calc_radial_nucP other not 3c procs
    real(kind=r8_kind),allocatable        :: gamma_help_co(:,:,:),gamma_arg_xc(:,:,:), &
                                             gamma_argXC(:,:)

    real(kind=r8_kind),     pointer   :: nuc(:,:,:),rotmat_eq(:,:,:)
    real(kind=r8_kind),     pointer   :: nuc_pc_timps(:,:,:)
    real(kind=r8_kind),     allocatable   :: nucEWGa(:,:,:,:),nucEW(:,:,:),nucEWGc(:,:,:,:)
    real(kind=r8_kind),pointer  :: potential(:,:,:), field(:,:,:,:), intermed_3c(:)
    real(kind=r8_kind),allocatable  :: potentialGc(:,:,:,:),potentialGa(:,:,:,:), &
                                       potential_sum(:,:,:)
    real(kind=r8_kind),parameter,public   :: pi         = 3.14159265358979324_r8_kind  &
                                           , very_big   = 1.0e100_r8_kind              &
                                           , very_small = 1.0e-100_r8_kind
    real(kind=r8_kind),parameter,private  :: zero       = 0.0_r8_kind,                 &
                                             one        = 1.0_r8_kind,                 &
                                             two        = 2.0_r8_kind,                 &
                                             four       = 4.0_r8_kind,                 &
                                             six        = 6.0_r8_kind
    real(kind=r8_kind), public            :: z,arg
    real(kind=r8_kind),    pointer        :: prim_int_ewpc(:,:,:,:)

    integer(kind=i4_kind)                 :: lc ! angular momentum of unique atom c
    integer(kind=i4_kind), public         :: naexps,nbexps,ncexps
    integer(kind=i4_kind)                 :: num,l_max
    integer(kind=i4_kind), public         :: alloc_stat(80) ! 36 empty now
    integer(kind=i4_kind)                 :: a,b,c
    integer(kind=i4_kind)                 :: i,j,lm,n_equals,j1,j2,j3,N_max,N,ix,nc,N_points
    integer(kind=i4_kind)                 :: ma,mb,mc,i1,i2,m_up,   ll,mm
    integer(kind=i4_kind)                 :: lmax_ch,k,n_jj,llam_m
    integer(kind=i4_kind)                 :: max_triple_terms

    logical, allocatable, public          :: cutoff(:,:)

    type, public :: csh_map_type
       integer (i4_kind), allocatable :: l(:) ! (lmax + 1)**2
       integer (i4_kind), allocatable :: m(:) ! (lmax + 1)**2
    end type csh_map_type

    type J_index_Complex
       complex(kind=c16_kind), pointer :: J(:,:,:,:,:,:)
    end type J_index_Complex
    type J_index_Real
       real(kind=r8_kind),     pointer :: J(:,:,:,:,:,:)
       real(kind=r8_kind),     pointer :: j_packed(:)
       integer(kind=i4_kind)           :: n_terms
       real(kind=r8_kind),     pointer :: m_c(:,:)
       real(kind=r8_kind),     pointer :: coeff(:)
       type(J_index_Complex), pointer  :: C(:)
    end type J_index_Real
    type, public :: jj_coeff_type
       integer(kind=i4_kind)           :: Num_M_terms
       integer(kind=i4_kind),  pointer :: m1(:)
       integer(kind=i4_kind),  pointer :: m2(:)
       integer(kind=i4_kind),  pointer :: m3(:)
       real(kind=r8_kind),     pointer :: value(:)
    end type jj_coeff_type
    type cr_table_type
       integer(kind=i4_kind)           :: N_terms
       integer(kind=i4_kind),  pointer :: sig(:,:)
       complex(kind=c16_kind), pointer :: Coeff(:)
    end type cr_table_type

    type (csh_map_type), public, protected :: csh_map ! -> solhrules_module
    type(unique_atom_type),   pointer            :: ua_pointer
    type(J_index_Complex), allocatable,target    :: Ang_C(:,:,:)
    type(J_index_Complex), allocatable,target    :: Ang_Cgc(:,:,:)
    type(J_index_Complex), allocatable,target    :: Ang_C3c(:,:,:)
    type(J_index_Complex), pointer               :: Ang_C_p(:,:,:)

    type(J_index_Real),       pointer            :: Ang_R(:,:,:)
    type(J_index_Real),       pointer            :: Ang_Rgc(:,:,:)
    type(J_index_Real),       pointer            :: Ang_R3c(:,:,:)
    type(J_index_Real),       pointer            :: Ang_R_p(:,:,:)

    complex(kind=c16_kind), pointer ::                       AngCMat(:,:,:,:)
    complex(kind=c16_kind), pointer,dimension(:,:,:,:,:)::   AngCMatGa, AngCMatGb
    complex(kind=c16_kind), pointer,dimension(:,:,:,:,:,:):: &
                                                AngCMatGGa,AngCMatGGb,AngCMatGaGb

    real(kind=r8_kind), pointer :: AngRMat(:,:,:,:)
    real(kind=r8_kind), pointer :: AngRMatGa(:,:,:,:,:),AngRMatGb(:,:,:,:,:)
    real(kind=r8_kind), pointer, dimension(:,:,:,:,:,:) :: &
                                 AngRMatGGa,AngRMatGGb,AngRMatGaGb

    real(kind=r8_kind), pointer :: AngR3cMat(:,:,:,:,:)
    real(kind=r8_kind),dimension(:,:,:,:,:,:), &
                         pointer::AngR3cMatGa,AngR3cMatGb ,AngR3cMatGx
    real(kind=r8_kind),dimension(:,:,:,:,:,:,:), &
       pointer::AngR3cMatGGa,AngR3cMatGGb,AngR3cMatGaGb

    real(kind=r8_kind),allocatable    :: exp_arg(:,:)

    type(jj_coeff_type),      allocatable,target :: jj_coeff(:,:,:)
    type(cr_table_type),      target             :: triple(-1:1,-1:1,-1:1)

    type(J_index_Real), allocatable,target       :: Radial(:)
    type(J_index_Real), allocatable,target       :: Radial3c(:)
    type(J_index_Real), allocatable,target       :: RadialGc(:)
    type(J_index_Real),       pointer            :: Radial_p(:)

    real(kind=r8_kind), pointer        :: radial_mat(:,:)
    real(kind=r8_kind), pointer        :: radialGc_mat(:,:)
    real(kind=r8_kind), pointer        :: radial3cmat(:,:,:), &
                                          radial3cmatG(:,:,:),radial3cmatGG(:,:,:)
    real(kind=r8_kind), pointer        :: radial3cmat_g(:,:,:) !gamma deriv
    real(kind=r8_kind), pointer        :: radialNuc_mat(:,:)
    real(kind=r8_kind), pointer        :: radialNuc_matGa(:,:)
    real(kind=r8_kind), pointer        :: radialNuc_matGGa(:,:)
    real(kind=r8_kind), allocatable :: radialN(:,:)
    real(kind=r8_kind),allocatable:: jj_fac(:,:,:)
    real(kind=r8_kind),allocatable:: jj_facGa(:,:,:)

    integer(kind=i4_kind) :: nc_fit_only
    integer(kind=i4_kind) :: alph_pw_min,alph_pw_max
    integer(kind=i4_kind) :: beta_pw_min,beta_pw_max
    integer(kind=i4_kind) :: fact2_pw_min,fact2_pw_max
    logical,allocatable ::   even_triangle(:,:,:)

    integer (i4_kind), private, parameter :: MAXFACN = 20
    integer (i8_kind), private :: fact(-1:MAXFACN) ! 13! doenst fit into 32 bits
    integer (i8_kind), public :: df(-1:MAXFACN)    ! 21! doenst fit into 64 bits

    integer(kind=i4_kind),allocatable :: radial_m_max(:)
    real(kind=r8_kind),allocatable:: alph_pw(:,:),beta_pw(:,:),fact2_pw(:,:)
    real(kind=r8_kind),allocatable:: i3_func_fac(:),i3_fac(:,:,:,:),nucR_fac(:,:,:,:), &
                                     nucR_facGa(:,:,:,:,:)
    real (r8_kind), allocatable :: bin_fac(:, :)
    real(kind=r8_kind),allocatable::nuc_n(:,:,:)
      logical,allocatable              ::m_c_required(:,:,:,:)
      logical,allocatable              ::mc_required(:,:,:,:)
      integer(kind=i4_kind), pointer   :: magn(:),eq_atom(:)
      integer(kind=i4_kind)            :: n_independent_fcts,i_ind,n_contributing_fcts,i_cnt
      integer(kind=i4_kind),allocatable :: csh_lm_of(:,:)
      integer(kind=i4_kind),allocatable:: mup_llam(:)
      real(kind=r8_kind), allocatable :: afac_sig(:,:,:),valjj(:,:,:,:)
      integer(kind=i4_kind),allocatable,dimension(:,:,:,:) :: m1jj,m2jj,m3jj
      integer(kind=i4_kind)::          maxNum_M_terms
!     integer(kind=i4_kind), private :: i_grad
      logical, public :: c_grad=.true.,l_print=.false.,l_time=.true.,a_grad=.true.
      integer(kind=i4_kind):: grad_dim,i_ma,index
      real(kind=r8_kind),allocatable    :: gamma_help_fac(:),r2_jjj_fac(:),gamma_help_fac_g(:)
      real(kind=r8_kind),allocatable    :: ghelp_shift_fac(:)
      real(kind=r8_kind),allocatable    :: jj_fac_SA(:,:,:,:)

      real(kind=r8_kind),allocatable,dimension(:,:,:) :: &
                               gamma_help_g,gamma_help_gG,gamma_help_gGG  ! used in radial_r2
      real(kind=r8_kind),allocatable:: jj_facG(:,:,:,:),jj_facGG(:,:,:,:) ! used in radial_r2
      real(kind=r8_kind),allocatable,dimension(:,:,:) :: gamma_help_r2    ! used in radial_r2

      public :: calc_3c_init!()

    !
    ! FIXME: which of these are re-exports from other modules?
    !

    ! these are used already in calc_3center.f90:
    public :: a
    public :: aexp_arr
    public :: afac_sig
    public :: A_factor
    public :: A_factorGa
    public :: A_factorGaGb
    public :: A_factorGb
    public :: A_factorGGa
    public :: A_factorGGb
    public :: alph_pw
    public :: alph_pw_max
    public :: AngCMat
    public :: AngCMatGa
    public :: AngCMatGaGb
    public :: AngCMatGb
    public :: AngCMatGGa
    public :: AngCMatGGb
    public :: AngR3cMat
    public :: AngR3cMatGa
    public :: AngR3cMatGaGb
    public :: AngR3cMatGb
    public :: AngR3cMatGGa
    public :: AngR3cMatGGb
    public :: AngRMat
    public :: AngRMatGa
    public :: AngRMatGaGb
    public :: AngRMatGb
    public :: AngRMatGGa
    public :: AngRMatGGb
    public :: b
    public :: beta_pw
    public :: beta_pw_max
    public :: bexp_arr
    public :: bin_fac
    public :: c
    public :: ca_dervsM
    public :: c_exp_arg2
    public :: coul_int
    public :: cpks_grad_mat
    public :: CSH
    public :: CSHG
    public :: CSHGaGb
    public :: CSHGb
    public :: CSHGGa
    public :: CSHGGb
    public :: cshgg_scalar
    public :: cshg_scalar
    public :: csh_lm_of
    public :: csh_scalar
    public :: czero
    public :: dervsM
    public :: dervs_mat
    public :: dervs_totsymM
    public :: dervs_totsymM_temp
    public :: do_rotation_eq
    public :: even_triangle
    public :: exp_arg
    public :: fact0
    public :: fact1
    public :: fact2
    public :: fact2_pw
    public :: fact2_pw_max
    public :: f_even_triangle
    public :: Ga
    public :: GaGb
    public :: gamma_arg
    public :: gamma_arg2
    public :: gamma_arg2_vec
    public :: gamma_arg_xc
    public :: gamma_argXC
    public :: gamma_help
    public :: gamma_help_co
    public :: gamma_help_fac
    public :: gamma_help_fac_g
    public :: gamma_help_g
    public :: gamma_help_gG
    public :: gamma_help_ggg
    public :: gamma_help_r2
    public :: gamma_var_fixed
    public :: Gb
    public :: GGa
    public :: GGb
    public :: GGt
    public :: ghelp_shift_fac
    public :: grad_dim
    public :: gradient_data_n_gradients
    public :: gradient_data_n_spin_gradients
    public :: gradient_index
    public :: gradM
    public :: grad_mat
    public :: grad_totsymM
    public :: g_shift_fac
    public :: gt
    public :: GtFvec
    public :: GtGa
    public :: GtGb
    public :: harmonic_var_fixed
    public :: i
    public :: i1
    public :: i2
    public :: i3_fac
    public :: i3_func_fac
    public :: ix
    public :: j
    public :: j1
    public :: j2
    public :: j3
    public :: jj_coeff
    public :: jj_fac
    public :: jj_facg
    public :: jj_facGa
    public :: jj_facgg
    public :: jj_fac_SA
    public :: k
    public :: lc
    public :: llam_m
    public :: llam_min
    public :: lm
    public :: l_max
    public :: lmax_ch
    public :: m1jj
    public :: m2jj
    public :: m3jj
    public :: ma
    public :: mb
    public :: mc
    public :: m_c_required
    public :: mc_required
    public :: m_up
    public :: mup_llam
    public :: nc
    public :: nc_fit_only
    public :: n_equals
    public :: new_3c_co
    public :: new_ewpc
    public :: new_nuc
    public :: n_independent_fcts
    public :: n_jj
    public :: n_jj_terms
    public :: n_max
    public :: nuc
    public :: nucR_fac
    public :: num
    public :: overlap_var_fixed
    public :: pointer_coul
    public :: pointer_prim_int
    public :: r2_jjj_fac
    public :: radangF
    public :: radang_temp
    public :: radang_temp_ga
    public :: radang_temp_gagb
    public :: radang_temp_gb
    public :: radang_temp_gga
    public :: radang_temp_ggb
    public :: radial3cmat
    public :: radial3cmat_g
    public :: radial3cmatG
    public :: radial3cmatGG
    public :: radialNuc_mat
    public :: radialNuc_matGa
    public :: radialNuc_matGGa
    public :: rotmat_eq
    public :: rsh
    public :: rshg
    public :: rshgg
    public :: temp_overlap
    public :: temp_overlap_nuc
    public :: timer_ang3c
    public :: timer_prod_nested
    public :: timer_rad3c
    public :: two_aexp_arr
    public :: two_bexp_arr
    public :: two_fact2_vec
    public :: ua_pointer
    public :: valjj

    ! these are used in calc_radial_r2.f90
    public :: gamma

    ! these are used in calc_3c_colc_setup.f90
    public :: coul_int_grad
    public :: coul_int_dervs
    public :: coul_int_grad_totsym
    public :: coul_int_ca_dervs
    public :: coul_int_dervs_totsym
    public :: ca_dervs_mat

    ! these are used in calc_3c_fitcontract.f90
    public :: index
    public :: cpks_grad_mat_p
    public :: grad_mat_p
    public :: ima, imb, imc
    public :: k_gr_0, k_gr_1

    ! subroutines:
    public :: Ang_Cmplx_to_RealGa
    public :: Ang_Cmplx_to_RealSA
#ifdef WITH_EPE
    public :: calc_3c_ewpc
#endif
    public :: calc_and_pack_3j
    public :: calc_angular_cmplx3cGa
    public :: calc_angular_cmplx3cSA
    public :: calc_csh_lm_of
    public :: calc_mc_required
    public :: calculate_factor
    public :: calculate_factor_Ga
    public :: calculate_factor_GaGb
    public :: calculate_factor_Gb
    public :: calculate_factor_GGa
    public :: calculate_factor_GGb
    public :: csh_lm_map
    public :: Radial_Ang_3cS
    public :: Radial_Ang_3cSA
    public :: Radial_Ang_nuc
    public :: set_triplets
    public :: shutdown_triplets

  contains

   subroutine calc_3c_init()
     !
     ! Initialize global lookup arrays for factorial
     ! and double factorial
     !
     implicit none
     ! *** end of interface ***

     integer (i4_kind) :: k

     fact(-1:0) = 1
     df(-1:0) = 1
     do k = 1, ubound (fact, 1)
       fact(k) = k * fact(k - 1)
       df(k) = k * df(k - 2)
     enddo
   end subroutine calc_3c_init

#ifdef WITH_EPE
   subroutine calc_3c_ewpc(la,lb,equalb)
     implicit none
     integer(kind=i4_kind), intent(in)::la,lb,equalb
     ! *** en dof interface ***

     integer (i4_kind) :: k_gr
     integer (i4_kind) :: k2dr
     integer (i4_kind) :: lp


           bin_fac = 0.0_r8_kind
           n_jj=n_jj_terms(la,lb,0)

           allocate(alph_pw(num,0:alph_pw_max),beta_pw(num,0:beta_pw_max) &
                   ,fact2_pw(num,0:fact2_pw_max),stat=alloc_stat(17))
           ASSERT(alloc_stat(17).eq.0)
           MEMLOG(size(alph_pw)+size(beta_pw)+size(fact2_pw))
           alloc_stat(17)=1

           do j1=0,alph_pw_max
            alph_pw(:,j1)=(two*aexp_arr(:))**j1
           enddo
           do j2=0,beta_pw_max
            beta_pw(:,j2)=(two*bexp_arr(:))**j2
           enddo
           do j3=0,fact2_pw_max
            fact2_pw(:,j3)=(two*fact2(:))**j3
           enddo

           allocate(mup_llam(0:La+Lb), &
            i3_func_fac(num),nucR_fac(num,0:la+lb,0:La,0:Lb),stat=alloc_stat(25))
            MEMLOG(size(mup_llam)+size(i3_func_fac)+size(nucR_fac))
            ASSERT(alloc_stat(25).eq.0)
            alloc_stat(25)=1

            llam_m=llam_min(la,lb,0)
            do m_up=0,la+lb
             i3_func_fac(:)=temp_overlap_nuc(:)*fact2_pw(:,m_up)
             do i1=0,La
              do i2=0,Lb
               nucR_fac(:,m_up,i1,i2)=i3_func_fac(:) * alph_pw(:,i1) * beta_pw(:,i2)
               enddo
             enddo
            enddo

    if(integralpar_3cob_grad) then
     grad_mat=0.0_r8_kind

     if(integralpar_dervs.or.ewpcdervs) then
      allocate( dervsM(num,2*Lb+1,2*La+1,6,6), stat=alloc_stat(15))
      ASSERT(alloc_stat(15).eq.0)
      MEMLOG(size(dervsM))
               dervsM=0.0_r8_kind
      N_max = la + lb +2
     else
      N_max = la + lb +1
     endif

     allocate(radang_temp(num,(2*La+1),(2*Lb+1)), &
              nucEWGa(num,2*la+1,2*lb+1,3), &
              nucEW(num,2*la+1,2*lb+1), stat=alloc_stat(47))
     MEMLOG(size(radang_temp)+size(nucEWGa)+size(nucEW))
     ASSERT(alloc_stat(47).eq.0)
     alloc_stat(47)=1
     nucEWGa=0.0_r8_kind

     A_factorGa=czero
    else
     N_max = la + lb
     allocate(radang_temp(num,(2*La+1),(2*Lb+1)), &
              nucEW(num,2*la+1,2*lb+1), stat=alloc_stat(47))
     MEMLOG(size(radang_temp)+size(nucEW))
     ASSERT(alloc_stat(47).eq.0)
    endif
     alloc_stat(47)=1

     nucEW=0.0_r8_kind

    nc = 0
    ewpc_unique_points : do i=1,ewpc_n
       z= ewpc_array(i)%z
       n_equals = ewpc_array(i)%n_equal_charges

       allocate ( gamma_help(num,1+N_max), gamma_arg2(num), stat=alloc_stat(18))
       MEMLOG(size(gamma_help)+size(gamma_arg2))
        ASSERT(alloc_stat(18).eq.0)
        alloc_stat(18)=1

       ewpc_equal_points  : do j=1,n_equals

         CSH(1:2,3,:)=czero
         CSH(3,1:2,:)=czero
         CSH(3,3,:)=czero

         A_factor = czero

         xc=ewpc_array(i)%position(:,j) ; nc=nc+1

          gamma_argXC(:,1)=gamma_arg(:,1)-xc(1)
          gamma_argXC(:,2)=gamma_arg(:,2)-xc(2)
          gamma_argXC(:,3)=gamma_arg(:,3)-xc(3)

          CSH(b,c,:) = CSH_scalar(l_max,xb-xc)
          CSH(c,a,:) = CSH_scalar(l_max,xc-xa)


          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(c,b,lm) = (-1)**lp * CSH(b,c,lm)
             CSH(a,c,lm) = (-1)**lp * CSH(c,a,lm)
          end do

          call calculate_factor(1,b,c,a,la)    ! A
          call calculate_factor(2,c,a,b,lb)    ! B
          call calculate_factor(3,a,b,c, 0)    ! C

           allocate(AngCMat(n_jj,-la:la, -lb:lb, 0:0), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("calc_3c AngCMat allocation failed")
           endif
           alloc_stat(22)=1

        if(integralpar_3cob_grad) then

         CSHG(b,c,:,:) = czero
         CSHG(c,a,:,:) = -CSHG_scalar(l_max,xc-xa)


          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSHG(c,b,lm,:) = (-1)**lp * CSHG(b,c,lm,:)
             CSHG(a,c,lm,:) = (-1)**lp * CSHG(c,a,lm,:)
          end do
          call calculate_factor_Ga(1,b,c,a,la)
          call calculate_factor_Ga(2,c,a,b,lb)    ! B
          call calculate_factor_Ga(3,a,b,c, 0)    ! C

         if(integralpar_dervs.or.ewpcdervs) then

           CSHGb(c,a,:,:) = czero  ! c,a do not depend on b
           CSHGb(a,c,:,:) = czero

!           CSHGb(b,c,:,:) = CSHG_scalar(l_max,xb-xc)

           CSHGGa(b,c,:,:,:) = czero ! (b,c) element do not depend on a
           CSHGGa(c,b,:,:,:) = czero

           CSHGGb(a,c,:,:,:) = czero ! (a,c) element do not depend on b
           CSHGGb(c,a,:,:,:) = czero

           CSHGaGb(a,c,:,:,:) = czero ! (a,c) element do not depend on b
           CSHGaGb(c,a,:,:,:) = czero

           CSHGaGb(c,b,:,:,:) = czero ! (c,b) element do not depend on a
           CSHGaGb(b,c,:,:,:) = czero

           CSHGGa(c,a,:,:,:) =  CSHGG_scalar(l_max,xc-xa)
           CSHGGb(c,b,:,:,:) =  CSHGG_scalar(l_max,xc-xb)

           do lm=1,(l_max+1)**2
!            CSHGb(c,b,lm,:) = (-1)**csh_map%l(lm) * CSHGb(b,c,lm,:)
            CSHGGa(a,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGa(c,a,lm,:,:)
            CSHGGb(b,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGb(c,b,lm,:,:)
           enddo


          call calculate_factor_Gb(1,b,c,a,la)    ! A
          call calculate_factor_Gb(2,c,a,b,lb)    ! B
          call calculate_factor_Gb(3,a,b,c,0)    ! C
          !   A_factorGb calculated

           call calculate_factor_GGa(1,b,c,a,la)    ! A
           call calculate_factor_GGa(2,c,a,b,lb)    ! B
           call calculate_factor_GGa(3,a,b,c,0)    ! C
           !  A_factorGGa calculated

           call calculate_factor_GGb(1,b,c,a,la)    ! A
           call calculate_factor_GGb(2,c,a,b,lb)    ! B
           call calculate_factor_GGb(3,a,b,c,0)    ! C
           !  A_factorGGb calculated

           call calculate_factor_GaGb(1,b,c,a,la)    ! A
           call calculate_factor_GaGb(2,c,a,b,lb)    ! B
           call calculate_factor_GaGb(3,a,b,c,0)    ! C
           !  A_factorGaGb calculated

         endif

           allocate(AngCMatGa(n_jj,-la:la, -lb:lb, 0:0,3), stat=alloc_stat(48))
           MEMLOG(size(AngCMatGa))

           if(alloc_stat(48).eq.0) then
           alloc_stat(48)=1
            AngCMatGa = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("48 AngCMat allocation failed")
           endif

           if(integralpar_dervs.or.ewpcdervs) then
           allocate(AngCMatGGa(n_jj,-la:la,-lb:lb,0:0,3,3), &
                    AngCMatGGb(n_jj,-la:la,-lb:lb,0:0,3,3), &
                    AngCMatGaGb(n_jj,-la:la,-lb:lb,0:0,3,3), &
                    AngCMatGb(n_jj,-la:la,-lb:lb,0:0,3), &
                                         stat=alloc_stat(48))
           MEMLOG(size(AngCMatGGa)*3+size(AngCMatGb))
                    ASSERT(alloc_stat(48).eq.0)
                           alloc_stat(48)=1

           AngCMatGGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGGb = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGaGb = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGb = (0.0_r8_kind,0.0_r8_kind)

             allocate(AngRMatGa(n_jj,2*La+1,2*Lb+1,1,3),     &
                      AngRMatGb(n_jj,2*La+1,2*Lb+1,1,3),     &
                      AngRMatGGa(n_jj,2*La+1,2*Lb+1,3,3,1),  &
                      AngRMatGGb(n_jj,2*La+1,2*Lb+1,3,3,1),  &
                      AngRMatGaGb(n_jj,2*La+1,2*Lb+1,3,3,1), &
                      stat=alloc_stat(49))
              ASSERT(alloc_stat(49).eq.0)
              MEMLOG(size(AngRMatGa)*2+size(AngRMatGGa)*3)
              AngRMatGGa=0.0_r8_kind
              AngRMatGGb=0.0_r8_kind
              AngRMatGaGb=0.0_r8_kind
              AngRMatGb=0.0_r8_kind
           else
             allocate(AngRMatGa(n_jj,2*La+1,2*Lb+1,1,3), stat=alloc_stat(49))
             MEMLOG(size(AngRMatGa))
              ASSERT(alloc_stat(49).eq.0)
           endif
              alloc_stat(49)=1
              AngRMatGa=0.0_r8_kind

        endif

           call calc_angular_cmplx3cGa(la,lb,0,n_jj,3, &    ! calc_3c_ewpc
                integralpar_gradients.and.ewpcdervs.or.integralpar_dervs, &
                                       AngCMatGb=AngCMatGb)
           ! AngCMatGGa calculated

           allocate(AngRMat(n_jj,2*La+1,2*Lb+1,1), stat=alloc_stat(27))
           ASSERT(alloc_stat(27).eq.0)
           MEMLOG(size(AngRMat))
           alloc_stat(27)=1
           AngRMat=0.0_r8_kind

           call Ang_Cmplx_to_RealGa(la,lb,0, &
                integralpar_gradients.and.ewpcdervs.or.integralpar_dervs, &
                                    AngCMatGb=AngCMatGb)
                                   ! AngRMatGGa calculated
           MEMLOG(-size(AngCMat))
           deallocate(AngCMat, stat=alloc_stat(22))
           ASSERT(alloc_stat(22).eq.0)

       if(integralpar_3cob_grad) then
           MEMLOG(-size(AngCMatGa))
           deallocate(AngCMatGa, stat=alloc_stat(48))
!          if(i.eq.1) print*, 'ewpc grad angle part calculated'
!
          if(integralpar_dervs.or.ewpcdervs) then
           MEMLOG(-size(AngCMatGGa)*3-size(AngCMatGb))
           deallocate(AngCMatGGa,AngCMatGGb,AngCMatGaGb,AngCMatGb, &
                      stat=alloc_stat(48))
          endif
           ASSERT(alloc_stat(48).eq.0)

          if(integralpar_dervs.or.ewpcdervs) then
           allocate(radialNuc_matGGa(num,n_jj), &
                    radialNuc_matGa(num,n_jj),  &
                    jj_facGa(num,0:la+lb,llam_m:la+lb), & ! calc_3c_ewpc
                    stat=alloc_stat(50))
           MEMLOG(size(radialNuc_matGGa)*2+size(jj_facGa))
           ASSERT(alloc_stat(50).eq.0)
          else
           allocate(radialNuc_matGa(num,n_jj), &
                    jj_facGa(num,0:la+lb,llam_m:la+lb), & ! calc_3c_ewpc
                    stat=alloc_stat(50))
           ASSERT(alloc_stat(50).eq.0)
           MEMLOG(size(radialNuc_matGa)+size(jj_facGa))
          endif
           alloc_stat(50)=1
        endif

            allocate(radialNuc_mat(num,n_jj), &                   ! 1
                     jj_fac(num,0:la+lb,llam_m:la+lb), &          ! 2
                     stat=alloc_stat(23))
               ASSERT(alloc_stat(23).eq.0)
               MEMLOG(size(radialNuc_mat)+size(jj_fac))
                      alloc_stat(23)=1 ! radialNuc_mat
                      alloc_stat(43)=1 ! jj_fac

            call calc_radial_nucP(la,lb)
            if(integralpar_3cob_grad) then
             call calc_radial_nucGa(la,lb) ! ewpc 1
             if(integralpar_dervs.or.ewpcdervs) call calc_radial_nucGGa(la,lb) !result: radialNuc_matGGa
            endif

             MEMLOG(-size(jj_fac))
             deallocate(jj_fac, stat=alloc_stat(43))  ! 2 3
              ASSERT(alloc_stat(43).eq.0)

            call dgemm('n','n',num,(2*La+1)*(2*Lb+1),n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMat(:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp,num)

            if(ewpcints) nucEW=nucEW+radang_temp

           if(integralpar_3cob_grad) then

             call Radial_Ang_EWGa(la,lb)  ! result:  nucEWGa without overlap contrib given below

             if(integralpar_dervs.or.ewpcdervs) then
             MEMLOG(-size(AngRMatGa)*2-size(AngRMatGGa)*3)
             deallocate(radialNuc_matGGa,radialNuc_matGa,jj_facGa,AngRMatGa, &
                        AngRMatGGa,AngRMatGGb,AngRMatGaGb,AngRMatGb,         &
                        stat=alloc_stat(50)) ! calc_3c_ewpc
             else
             MEMLOG(-size(AngRMatGa)-size(radialNuc_matGa)-size(jj_facGa))
             deallocate(radialNuc_matGa,jj_facGa,AngRMatGa, stat=alloc_stat(50)) ! calc_3c_ewpc
             endif
             ASSERT(alloc_stat(50).eq.0)
                alloc_stat(49) = 0 ! AngRMatGa,AngRMatGGa
            endif

            MEMLOG(-size(radialNuc_mat)-size(AngRMat))
            deallocate(radialNuc_mat,AngRMat, stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
            alloc_stat(27)=0 ! AngRMat

    enddo ewpc_equal_points
       MEMLOG(-size(gamma_help)-size(gamma_arg2))
       deallocate ( gamma_help, gamma_arg2, stat=alloc_stat(18))
       ASSERT(alloc_stat(18).eq.0)
   enddo ewpc_unique_points
!  DPRINT 'all ewpc_unique_points treated'

    if(integralpar_3cob_grad) then
         do ma=1,2*La+1
         do mb=1,2*Lb+1

         if(ewpcgrads) then
          do k_gr=1,3
          nucEWGa(:,ma,mb,k_gr)=nucEWGa(:,ma,mb,k_gr) &
                - two_fact2_vec(:,k_gr)*nucEW(:,ma,mb)
          enddo
         endif

         abmap_nuc_dervs: if(integralpar_dervs.and.ewpcdervs) then
           do k_gr=1,6
            do k2dr=1,6
            dervs_mat(:,:,mb,ma,k_gr,k2dr)=dervs_mat(:,:,mb,ma,k_gr,k2dr) &
                -unpack(dervsM(:,mb,ma,k_gr,k2dr),cutoff,zero)
            enddo
           enddo
         endif abmap_nuc_dervs
         end do

      enddo
!    print*,'nucEWGa nucEWGa',sum(abs(nucEWGa(:,:,:,1))),sum(abs(nucEWGa(:,:,:,2))),sum(abs(nucEWGa(:,:,:,3)))
    else
    endif

       MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw))
       deallocate(alph_pw,beta_pw,fact2_pw,stat=alloc_stat(17))
       ASSERT(alloc_stat(17).eq.0)

       MEMLOG(-size(nucR_fac)-size(i3_func_fac)-size(mup_llam))
       deallocate(nucR_fac,i3_func_fac,mup_llam,stat=alloc_stat(25))
       ASSERT(alloc_stat(25).eq.0)

!  print*,' start EWgc calculations'
   EWgc:if(integralpar_3cob_grad) then
    if(ewpcgrads) then
       allocate(nucEWGc(num,2*la+1,2*lb+1,3), stat=alloc_stat(51))
       MEMLOG(size(nucEWGc))
       ASSERT(alloc_stat(51).eq.0)
       alloc_stat(51)=1
       nucEWGc=0.0_r8_kind

           bin_fac = 0.0_r8_kind
           n_jj=n_jj_terms(la,lb,1)
           allocate(alph_pw(num,0:alph_pw_max),beta_pw(num,0:beta_pw_max) &
                   ,fact2_pw(num,0:fact2_pw_max),stat=alloc_stat(52))
           MEMLOG(size(alph_pw)+size(beta_pw)+size(fact2_pw))
           ASSERT(alloc_stat(52).eq.0)
           alloc_stat(52)=1
           do j1=0,alph_pw_max
           alph_pw(:,j1)=(two*aexp_arr(:))**j1
           enddo
           do j2=0,beta_pw_max
           beta_pw(:,j2)=(two*bexp_arr(:))**j2
           enddo
           do j3=0,fact2_pw_max
            fact2_pw(:,j3)=(two*fact2(:))**j3
           enddo
           allocate(mup_llam(0:La+Lb+1), i3_func_fac(num), &
                    nucR_fac(num,0:la+lb,0:La+1,0:Lb+1),stat=alloc_stat(53))
           ASSERT(alloc_stat(53).eq.0)
           MEMLOG(size(mup_llam)+size(i3_func_fac)+size(nucR_fac))
           alloc_stat(53)=1
           llam_m=llam_min(la,lb,1)
                do m_up=0,la+lb
                   i3_func_fac(:)=temp_overlap_nuc(:)*fact2_pw(:,m_up)
                   do i1=0,La+1
                   do i2=0,Lb+1
                   nucR_fac(:,m_up,i1,i2)=i3_func_fac(:) *alph_pw(:,i1)*beta_pw(:,i2)
                  enddo
                  enddo
                  enddo
    nc = 0
    EWgc_unique_points : do i=1,ewpc_n
       n_equals = ewpc_array(i)%n_equal_charges
       z= ewpc_array(i)%z

        N_max = la + lb  + 1
        allocate ( gamma_help(num,1+N_max), gamma_arg2(num), stat=alloc_stat(54))
        MEMLOG(size(gamma_help)+size(gamma_arg2))
        ASSERT(alloc_stat(54).eq.0)
        alloc_stat(54)=1

       EWGc_equal_points  : do j=1,n_equals

          CSH      = czero
          A_factor = czero
          xc=ewpc_array(i)%position(:,j)
          nc=nc+1
!          if(i.eq.1) xc(3)=xc(3)+0.001_r8_kind

          CSH(a,b,:) = CSH_scalar(l_max,xa-xb)
          CSH(b,c,:) = CSH_scalar(l_max,xb-xc)
          CSH(c,a,:) = CSH_scalar(l_max,xc-xa)
          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(b,a,lm) = (-1)**lp * CSH(a,b,lm)
             CSH(c,b,lm) = (-1)**lp * CSH(b,c,lm)
             CSH(a,c,lm) = (-1)**lp * CSH(c,a,lm)
          end do
          call calculate_factor(1,b,c,a,la)    ! A
          call calculate_factor(2,c,a,b,lb)    ! B
          call calculate_factor(3,a,b,c, 1)    ! C
          !**************************************!

           allocate(AngCMat(n_jj,-la:la, -lb:lb, -1:1), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("calc_3c AngCMat allocation failed")
           endif
            alloc_stat(22)=1
           call calc_angular_cmplx3cPP(la,lb,1,AngCMat,n_jj,3) !!!! 2

           allocate(AngRMat(n_jj,2*La+1,2*Lb+1,3), stat=alloc_stat(27))
           ASSERT(alloc_stat(27).eq.0)
           MEMLOG(size(AngRMat))
           alloc_stat(27)=1
           AngRMat=0.0_r8_kind

           call Ang_Cmplx_to_RealPP(la,lb,1)
           MEMLOG(-size(AngCMat))
           deallocate(AngCMat, stat=alloc_stat(22))
           ASSERT(alloc_stat(22).eq.0)

            allocate(radialNuc_mat(num,n_jj), &
                     jj_fac(num,0:la+lb,llam_m:la+lb+1), &
                                      stat=alloc_stat(23))
            MEMLOG(size(radialNuc_mat)+size(jj_fac))
               ASSERT(alloc_stat(23).eq.0)
               alloc_stat(23)=1
            call calc_radial_potG(la,lb,1) ! 1 jj_fac alloc
            call Radial_Ang_EWGG(la,lb,1)

            MEMLOG(-size(radialNuc_mat)-size(jj_fac)-size(AngRMat))
            deallocate(radialNuc_mat,jj_fac,AngRMat, stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
                   alloc_stat(27)=0 ! AngRMat

       enddo EWGc_equal_points
       MEMLOG(-size(gamma_help)-size(gamma_arg2))
       deallocate ( gamma_help, gamma_arg2, stat=alloc_stat(54))
       ASSERT(alloc_stat(54).eq.0)
    enddo EWgc_unique_points

   allocate(nuc_grad_gr(num,2*lb+1,2*la+1,6),stat=alloc_stat(45))
   ASSERT(alloc_stat(45).eq.0)
   MEMLOG(size(nuc_grad_gr))
    alloc_stat(45)=1
     do ma=1,2*la+1
      do mb=1,2*lb+1
       if(moving_a) then
        nuc_grad_gr(:,mb,ma,1)=nucEWGa(:,ma,mb,1)
        nuc_grad_gr(:,mb,ma,2)=nucEWGa(:,ma,mb,2)
        nuc_grad_gr(:,mb,ma,3)=nucEWGa(:,ma,mb,3)
       endif
       if (moving_b) then
        nuc_grad_gr(:,mb,ma,4)=-nucEWGa(:,ma,mb,1)-nucEWGc(:,ma,mb,2)
        nuc_grad_gr(:,mb,ma,5)=-nucEWGa(:,ma,mb,2)-nucEWGc(:,ma,mb,3)
        nuc_grad_gr(:,mb,ma,6)=-nucEWGa(:,ma,mb,3)-nucEWGc(:,ma,mb,1)
        endif
      enddo
     enddo
!print*,'nuc_grad_gr b',sum(abs(nuc_grad_gr(:,:,:,4))),sum(abs(nuc_grad_gr(:,:,:,5))),sum(abs(nuc_grad_gr(:,:,:,6)))

     ! what is k_gr_0 and k_gr_1?
     do k_gr=k_gr_0,k_gr_1
        do ma=1,2*la+1
           do mb=1,2*lb+1
              grad_mat(:,:,mb,ma,k_gr)=grad_mat(:,:,mb,ma,k_gr)-&
                   unpack(nuc_grad_gr(:,mb,ma,k_gr),cutoff,zero)
           enddo
        enddo
     end do
    MEMLOG(-size(nuc_grad_gr))
    deallocate(nuc_grad_gr,stat=alloc_stat(45))
    ASSERT(alloc_stat(45).eq.0)

    if(integralpar_relativistic.and.pseudopot_present) then
!!     pointer_prim_int=>prim_int_2cob_pseudo_grad
    else
!!     pointer_prim_int=>prim_int_3cob_grad
    endif

!!     if(new_ewpc) call add_nuc_or_pvsp(equalb,grad_mat)

      MEMLOG(-size(mup_llam)-size(i3_func_fac)-size(nucR_fac))
      deallocate(mup_llam, i3_func_fac, nucR_fac,stat=alloc_stat(53))
      ASSERT(alloc_stat(53).eq.0)
      MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw))
      deallocate(alph_pw,beta_pw,fact2_pw,stat=alloc_stat(52))
      ASSERT(alloc_stat(52).eq.0)
       MEMLOG(-size(nucEWGc))
       deallocate(nucEWGc, stat=alloc_stat(51))
       ASSERT(alloc_stat(51).eq.0)

     MEMLOG(-size(nucEWGa)-size(nucEW)-size(radang_temp))
     deallocate(nucEWGa,nucEW,radang_temp, stat=alloc_stat(47))
       ASSERT(alloc_stat(47).eq.0)

     if(integralpar_dervs.or.ewpcdervs) then
      MEMLOG(-size(dervsM))
      deallocate( dervsM,stat=alloc_stat(15))
      ASSERT(alloc_stat(15).eq.0)
     endif

    endif

   else EWgc
    if(ewpcints) then

     if(pseudopot_present.and.integralpar_relativistic) then
!     prim_int_ewpc=>prim_int_2cob_nuc_pseudo
     else
!     prim_int_ewpc=>prim_int_2cob_nuc
     endif
!      do mb=1,2*lb+1
!         do ma=1,2*la+1
!         prim_int_ewpc(:,:,mb,ma)=prim_int_ewpc(:,:,mb,ma) &
!                                    +unpack(nucEW(:,ma,mb),cutoff,zero)
!         enddo
!      end do
     endif
     MEMLOG(-size(nucEW)-size(radang_temp))
     deallocate(nucEW,radang_temp, stat=alloc_stat(47))
     ASSERT(alloc_stat(47).eq.0)
   endif EWgc

   end subroutine calc_3c_ewpc
#endif

    subroutine calc_3c_solv(la,lb,equalb)
      implicit none
      integer(kind=i4_kind),intent(in):: la,lb,equalb
      ! *** end of interface ***

      integer (i4_kind) :: k_gr, lp

!    potent_calc: if(integralpar_2cob_potential.or.integralpar_solv_grad) then

!      allocate(temp_overlap_nuc(num),gamma_argXC(num,3), stat=alloc_stat(12))
!      ASSERT(alloc_stat(12).eq.0)
!      alloc_stat(12)=1

!      temp_overlap_nuc = exp(-fact2*arg) / &
!                    (fact0*sqrt(aexp_arr**La*df(2*La-1)*bexp_arr**Lb*df(2*Lb-1) )) &
!                                    * 2*pi * (4*fact1(:)/pi2)**0.75_r8_kind

           bin_fac = 0.0_r8_kind
           n_jj=n_jj_terms(la,lb,0)
           allocate(alph_pw(num,0:alph_pw_max),beta_pw(num,0:beta_pw_max) &
                   ,fact2_pw(num,0:fact2_pw_max),stat=alloc_stat(17))
           MEMLOG(size(alph_pw)+size(beta_pw)+size(fact2_pw))
           ASSERT(alloc_stat(17).eq.0)
           alloc_stat(17)=1

           do j1=0,alph_pw_max
            alph_pw(:,j1)=(two*aexp_arr(:))**j1
           enddo
           do j2=0,beta_pw_max
            beta_pw(:,j2)=(two*bexp_arr(:))**j2
           enddo
           do j3=0,fact2_pw_max
            fact2_pw(:,j3)=(two*fact2(:))**j3
           enddo
           allocate(mup_llam(0:La+Lb), i3_func_fac(num), &
                    nucR_fac(num,0:la+lb,0:La,0:Lb),stat=alloc_stat(25))
           MEMLOG(size(mup_llam)+size(i3_func_fac)+size(nucR_fac))
           ASSERT(alloc_stat(25).eq.0)
           alloc_stat(25)=1

           llam_m=llam_min(la,lb,0)
                do m_up=0,la+lb
                   i3_func_fac(:)=fact2_pw(:,m_up)*temp_overlap_nuc(:)
                   do i1=0,La
                   do i2=0,Lb
                   nucR_fac(:,m_up,i1,i2)=i3_func_fac(:) *    alph_pw(:,i1) * beta_pw(:,i2)
                  enddo
                  enddo
                  enddo

    nc = 0
       if(integralpar_solv_grad) then
         allocate(potentialGa(num,2*la+1,2*lb+1,3), &
                  potential_sum(num,2*la+1,2*lb+1),stat=alloc_stat(40))
         if (alloc_stat(40).ne.0) call error_handler  ("40 allocation failed")
             alloc_stat(40)=1
         potentialGa=0.0_r8_kind
         potential_sum=0.0_r8_kind
         N_max = la + lb +1
         N_points=to_calc_grads%n_points

        else
         N_points=N_points_pot
         N_max = la + lb
        endif

        allocate(potential(num,2*la+1,2*lb+1),intermed_3c(num), stat=alloc_stat(28))
        if (alloc_stat(28).ne.0) call error_handler  ("28: potential allocation failed")
        alloc_stat(28)=1

          CSH(1:2,1:2,:)=czero
          CSH(a,b,:) = CSH_scalar(l_max,xa-xb)

          if(integralpar_solv_grad) then
           CSHG      =czero
           A_factorGa=czero
           CSHG(a,b,:,:) = CSHG_scalar(l_max,xa-xb,i_p=1)
          endif

          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(b,a,lm) = (-1)**lp * CSH(a,b,lm)
             if(integralpar_solv_grad) CSHG(b,a,lm,:) = (-1)**lp * CSHG(a,b,lm,:)
          end do

    unique_points : do i=1,N_points
    potential=0.0_r8_kind ! to be calculated for each i point
       if(integralpar_solv_grad) then
        n_equals=to_calc_grads%n_equal(i)
        z=to_calc_grads%Q(i)
       else
        n_equals = point_in_space(i)%n_equal_points
        z= one
       endif
        allocate ( gamma_help(num,1+N_max), gamma_arg2(num), stat=alloc_stat(18))
        MEMLOG(size(gamma_help)+size(gamma_arg2))
        ASSERT(alloc_stat(18).eq.0)
        alloc_stat(18)=1
       equal_points  : do j=1,n_equals

    CSH(1:2,3,:)=czero
    CSH(3,1:2,:)=czero
    CSH(3,3,:)=czero

          A_factor = czero
       if(integralpar_solv_grad) then
          xc =to_calc_grads%xyz(i,j,:)
       else
          xc=point_in_space(i)%position(:,j)
       endif
          nc=nc+1
!         if(i.eq.1) xc(3)=xc(3)+0.001_r8_kind

          gamma_argXC(:,1)=gamma_arg(:,1)-xc(1)
          gamma_argXC(:,2)=gamma_arg(:,2)-xc(2)
          gamma_argXC(:,3)=gamma_arg(:,3)-xc(3)


          CSH(b,c,:) = CSH_scalar(l_max,xb-xc)
          CSH(c,a,:) = CSH_scalar(l_max,xc-xa)

          if(integralpar_solv_grad) CSHG(b,c,:,:) = czero
          if(integralpar_solv_grad) CSHG(c,a,:,:) = -CSHG_scalar(l_max,xc-xa)

          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(c,b,lm) = (-1)**lp * CSH(b,c,lm)
             CSH(a,c,lm) = (-1)**lp * CSH(c,a,lm)
          end do
          call calculate_factor(1,b,c,a,la)    ! A
          call calculate_factor(2,c,a,b,lb)    ! B
          call calculate_factor(3,a,b,c, 0)    ! C

           allocate(AngCMat(n_jj,-la:la, -lb:lb, 0:0), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("calc_3c AngCMat allocation failed")
           endif
           alloc_stat(22)=1


        if(integralpar_solv_grad) then
          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSHG(c,b,lm,:) = (-1)**lp * CSHG(b,c,lm,:)
             CSHG(a,c,lm,:) = (-1)**lp * CSHG(c,a,lm,:)
          end do
          call calculate_factor_Ga(1,b,c,a,la)
          call calculate_factor_Ga(2,c,a,b,lb)    ! B
          call calculate_factor_Ga(3,a,b,c, 0)    ! C
           allocate(AngCMatGa(n_jj,-la:la, -lb:lb, 0:0,3), stat=alloc_stat(37))
           MEMLOG(size(AngCMatGa))
           if(alloc_stat(37).eq.0) then
            AngCMatGa = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("37 AngCMat allocation failed")
           endif
           alloc_stat(37)=1
        endif

            call calc_angular_cmplx3cGa(la,lb,0,n_jj,3,integralpar_dervs) ! in solvation part

           allocate(AngRMat(n_jj,2*La+1,2*Lb+1,1), stat=alloc_stat(27))
           ASSERT(alloc_stat(27).eq.0)
           MEMLOG(size(AngRMat))
           alloc_stat(27)=1
           AngRMat=0.0_r8_kind

           if(integralpar_solv_grad) then
              allocate(AngRMatGa(n_jj,2*La+1,2*Lb+1,1,3), stat=alloc_stat(38))
              ASSERT(alloc_stat(38).eq.0)
              MEMLOG(size(AngRMatGa))
              alloc_stat(38)=1
              AngRMatGa=0.0_r8_kind
           endif
           call Ang_Cmplx_to_RealGa(la,lb,0,integralpar_dervs)  ! solv
           MEMLOG(-size(AngCMat))
           deallocate(AngCMat, stat=alloc_stat(22))
           ASSERT(alloc_stat(22).eq.0)

            allocate(radialNuc_mat(num,n_jj), &
                     jj_fac(num,0:la+lb,llam_m:la+lb), &
                                   stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
            MEMLOG(size(radialNuc_mat)+size(jj_fac))
               alloc_stat(23)=1 ! radialNuc_mat
               alloc_stat(43)=1 ! jj_fac

       if(integralpar_solv_grad) then
!          if(i.eq.1.and.j.eq.1) print*,AngRMat(1,1,1,1),AngRMatGa(1,1,1,1,1),' a da'
           MEMLOG(-size(AngCMatGa))
           deallocate(AngCMatGa, stat=alloc_stat(37))
           ASSERT(alloc_stat(37).eq.0)
            allocate(radialNuc_matGa(num,n_jj), &
                     jj_facGa(num,0:la+lb,llam_m:la+lb), & ! temp for calc_radial_nucGa calc_3c_solv
                     stat=alloc_stat(39))
            ASSERT(alloc_stat(39).eq.0)
            MEMLOG(size(radialNuc_matGa)+size(jj_facGa))
            alloc_stat(39)=1
       endif

            call calc_radial_nucP(la,lb)
            if(integralpar_solv_grad) call calc_radial_nucGa(la,lb)  ! 2
            call Radial_Ang_PotPP(la,lb,i)

        if(integralpar_solv_grad) then
            call Radial_Ang_PotGa(la,lb,i)

            MEMLOG(-size(radialNuc_matGa)-size(jj_facGa)-size(AngRMatGa))
            deallocate(radialNuc_matGa,jj_facGa,AngRMatGa, stat=alloc_stat(39)) !calc_3c_solv
            ASSERT(alloc_stat(39).eq.0)
               alloc_stat(38)=0 ! AngRMatGa
        endif

            MEMLOG(-size(jj_fac)-size(radialNuc_mat)-size(AngRMat))
            deallocate(jj_fac,radialNuc_mat,AngRMat,stat=alloc_stat(43))
            ASSERT(alloc_stat(43).eq.0)
               alloc_stat(23)=0 ! radialNuc_mat
               alloc_stat(27)=0 ! AngRMat
    end do equal_points

       MEMLOG(-size(gamma_help)-size(gamma_arg2))
       deallocate ( gamma_help, gamma_arg2, stat=alloc_stat(18))
       ASSERT(alloc_stat(18).eq.0)

      if(integralpar_2cob_potential) then
       do mb=1,2*lb+1
          do ma=1,2*la+1
             intermed_3c(:)=potential(:,ma,mb)
             prim_int_2cob_poten(:,:,i,mb,ma)=unpack(intermed_3c(:),cutoff,zero)
          enddo
       end do
      endif


       if(integralpar_solv_grad) then
       potential_sum=potential_sum+potential
         do ma=1,2*La+1
         do mb=1,2*Lb+1
          do k_gr=1,3
          potentialGa(:,ma,mb,k_gr)=potentialGa(:,ma,mb,k_gr) &
                   - potential(:,ma,mb)*fact2*two*(xa(k_gr)-xb(k_gr))
          enddo
         end do
      end do
      endif
   end do unique_points
!!!     if(a_grad.and.integralpar_solv_grad)  &
!!!                               print*, potential_sum(1,1,1),potentialGa(1,1,1,1),' p dp'
           MEMLOG(-size(mup_llam)-size(nucR_fac)-size(i3_func_fac))
           deallocate(mup_llam,nucR_fac,i3_func_fac,stat=alloc_stat(25))
           ASSERT(alloc_stat(25).eq.0)
           MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw))
           deallocate(alph_pw,beta_pw,fact2_pw,stat=alloc_stat(17))
           ASSERT(alloc_stat(17).eq.0)
    g_pot: if(integralpar_solv_grad) then
       allocate(potentialGc(num,2*la+1,2*lb+1,3), stat=alloc_stat(30))
       ASSERT(alloc_stat(30).eq.0)
       MEMLOG(size(potentialGc))
        alloc_stat(30)=1
        potentialGc=0.0_r8_kind
           bin_fac = 0.0_r8_kind
           n_jj=n_jj_terms(la,lb,1)
           allocate(alph_pw(num,0:alph_pw_max),beta_pw(num,0:beta_pw_max) &
                   ,fact2_pw(num,0:fact2_pw_max),stat=alloc_stat(17))
           ASSERT(alloc_stat(17).eq.0)
           MEMLOG(+size(alph_pw)+size(beta_pw)+size(fact2_pw))
           alloc_stat(17)=1
           do j1=0,alph_pw_max
            alph_pw(:,j1)=(two*aexp_arr(:))**j1
           enddo
           do j2=0,beta_pw_max
            beta_pw(:,j2)=(two*bexp_arr(:))**j2
           enddo
           do j3=0,fact2_pw_max
            fact2_pw(:,j3)=(two*fact2(:))**j3
           enddo
           allocate(mup_llam(0:La+Lb+1), i3_func_fac(num), &
                    nucR_fac(num,0:la+lb,0:La+1,0:Lb+1),stat=alloc_stat(25))
           ASSERT(alloc_stat(25).eq.0)
           MEMLOG(size(mup_llam)+size(i3_func_fac)+size(nucR_fac))
           alloc_stat(25)=1
           llam_m=llam_min(la,lb,1)
                do m_up=0,la+lb
                   i3_func_fac(:)=temp_overlap_nuc(:)*fact2_pw(:,m_up)
                   do i1=0,La+1
                   do i2=0,Lb+1
                   nucR_fac(:,m_up,i1,i2)=i3_func_fac(:) * alph_pw(:,i1) * beta_pw(:,i2)
                  enddo
                  enddo
                  enddo
    nc = 0
    N_points=to_calc_grads%n_points
    gpot_unique_points : do i=1,N_points
        n_equals=to_calc_grads%n_equal(i)
        z=to_calc_grads%Q(i)

        N_max = la + lb  + 1
        allocate ( gamma_help(num,1+N_max), gamma_arg2(num), stat=alloc_stat(18))
        MEMLOG(size(gamma_help)+size(gamma_arg2))
        ASSERT(alloc_stat(18).eq.0)
        alloc_stat(18)=1

       gpot_equal_points  : do j=1,n_equals

          CSH      = czero
          A_factor = czero
          xc =to_calc_grads%xyz(i,j,:)
          nc=nc+1
!          if(i.eq.1) xc(3)=xc(3)+0.001_r8_kind

          !*************** A B C ****************!
          CSH(a,b,:) = CSH_scalar(l_max,xa-xb)
          CSH(b,c,:) = CSH_scalar(l_max,xb-xc)
          CSH(c,a,:) = CSH_scalar(l_max,xc-xa)
          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(b,a,lm) = (-1)**lp * CSH(a,b,lm)
             CSH(c,b,lm) = (-1)**lp * CSH(b,c,lm)
             CSH(a,c,lm) = (-1)**lp * CSH(c,a,lm)
          end do
          call calculate_factor(1,b,c,a,la)    ! A
          call calculate_factor(2,c,a,b,lb)    ! B
          call calculate_factor(3,a,b,c, 1)    ! C
          !**************************************!

           allocate(AngCMat(n_jj,-la:la, -lb:lb, -1:1), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("calc_3c AngCMat allocation failed")
           endif
            alloc_stat(22)=1
           call calc_angular_cmplx3cPP(la,lb,1,AngCMat,n_jj,3) !!! 3

           allocate(AngRMat(n_jj,2*La+1,2*Lb+1,3), stat=alloc_stat(27))
           ASSERT(alloc_stat(27).eq.0)
           MEMLOG(size(AngRMat))
           alloc_stat(27)=1
           AngRMat=0.0_r8_kind

           call Ang_Cmplx_to_RealPP(la,lb,1)
           MEMLOG(-size(AngCMat))
           deallocate(AngCMat, stat=alloc_stat(22))
           ASSERT(alloc_stat(22).eq.0)

            allocate(radialNuc_mat(num,n_jj), &
                     jj_fac(num,0:la+lb,llam_m:la+lb+1), &
                     stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
            MEMLOG(size(radialNuc_mat)+size(jj_fac))
               alloc_stat(23)=1
            call calc_radial_potG(la,lb,1) ! 2 jj_fac alloc
            call Radial_Ang_PotGG(la,lb,1)
            MEMLOG(-size(radialNuc_mat)-size(jj_fac)-size(AngRMat))
            deallocate(radialNuc_mat,jj_fac,AngRMat, stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
            alloc_stat(27)=0
       enddo gpot_equal_points
       MEMLOG(-size(gamma_help)-size(gamma_arg2))
       deallocate ( gamma_help, gamma_arg2, stat=alloc_stat(18))
       ASSERT(alloc_stat(18).eq.0)
    enddo gpot_unique_points

  if(integralpar_solv_grad) then
   allocate(solv_grad_gr(num,2*lb+1,2*la+1,6),stat=alloc_stat(45))
   ASSERT(alloc_stat(45).eq.0)
   MEMLOG(size(solv_grad_gr))
    alloc_stat(45)=1
     do ma=1,2*la+1
      do mb=1,2*lb+1
       if(moving_a) then
        solv_grad_gr(:,mb,ma,1)=potentialGa(:,ma,mb,1)
        solv_grad_gr(:,mb,ma,2)=potentialGa(:,ma,mb,2)
        solv_grad_gr(:,mb,ma,3)=potentialGa(:,ma,mb,3)
       endif
       if (moving_b) then
        solv_grad_gr(:,mb,ma,4)=-potentialGa(:,ma,mb,1)-potentialGc(:,ma,mb,2)
        solv_grad_gr(:,mb,ma,5)=-potentialGa(:,ma,mb,2)-potentialGc(:,ma,mb,3)
        solv_grad_gr(:,mb,ma,6)=-potentialGa(:,ma,mb,3)-potentialGc(:,ma,mb,1)
        endif
      enddo
     enddo
     do k_gr=k_gr_0,k_gr_1
!       print*,solv_grad_gr(1,1,1,k_gr), ' solv gr ab'
        do ma=1,2*la+1
           do mb=1,2*lb+1
              grad_mat(:,:,mb,ma,k_gr)=grad_mat(:,:,mb,ma,k_gr)-&
                   unpack(solv_grad_gr(:,mb,ma,k_gr),cutoff,zero)
           enddo
        enddo
     end do
    MEMLOG(-size(solv_grad_gr))
    deallocate(solv_grad_gr,stat=alloc_stat(45))
    if(alloc_stat(45).ne.0) call error_handler("45: deallocate solv_grad_gr failed")

     pointer_prim_int=>prim_int_3cob_solv_grad
     if(new_solv_grad) call add_nuc_or_pvsp(equalb,grad_mat) !!!
    endif


           MEMLOG(-size(mup_llam)-size(i3_func_fac)-size(nucR_fac))
           deallocate(mup_llam, i3_func_fac, nucR_fac,stat=alloc_stat(25))
           ASSERT(alloc_stat(25).eq.0)

           MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw))
           deallocate(alph_pw,beta_pw,fact2_pw,stat=alloc_stat(17))
           if(alloc_stat(17) /= 0) call error_handler("deallocate g_pw failed")
       if(l_print) print*,'potentialGc ',potentialGc(1,1,1,:)
       MEMLOG(-size(potentialGc))
       deallocate(potentialGc, stat=alloc_stat(30))
       ASSERT(alloc_stat(30).eq.0)
   endif g_pot

!       deallocate(temp_overlap_nuc,gamma_argXC, stat=alloc_stat(12))
!       if(alloc_stat(12) /= 0) call error_handler("pot: temp_overlap_nuc deallocation failed")

       deallocate(potential,intermed_3c,stat=alloc_stat(28))
       if(alloc_stat(28)/=0) call error_handler("28 potential deallocation  failed")

       if(integralpar_solv_grad) then
         deallocate(potentialGa,potential_sum,stat=alloc_stat(40))
         if(alloc_stat(40)/=0) call error_handler("40 potentialGa deallocation  failed")
       endif

!    endif potent_calc
   end subroutine calc_3c_solv

    subroutine calculate_factor(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm
      ! A_factor is a tenzor of rank L => M=-L,L
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factor(n,LM,j)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factor(n,LM,j) = A_factor(n,LM,j)   +        &
                       CSH(a,c,csh_lm_of(L-j,M-mm) )    *        &
                       CSH(b,c,csh_lm_of(j,mm) ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do

    end subroutine calculate_factor

    subroutine calculate_factor_Ga(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm
      ! A_factor is a tenzor of rank L => M=-L,L
      if(a.eq.1) then
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGa(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGa(n,LM,j,:) = A_factorGa(n,LM,j,:)   +        &
                       CSHG(a,c,csh_lm_of(L-j,M-mm),: )    *        &
                       CSH(b,c,csh_lm_of(j,mm) ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
        elseif(b.eq.1) then
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGa(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGa(n,LM,j,:) = A_factorGa(n,LM,j,:)   +        &
                       CSH(a,c,csh_lm_of(L-j,M-mm) )    *        &
                       CSHG(b,c,csh_lm_of(j,mm),: ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
        else
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGa(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGa(n,LM,j,:) = A_factorGa(n,LM,j,:)   +        &
                      (CSH(a,c,csh_lm_of(L-j,M-mm) )*CSHG(b,c,csh_lm_of(j,mm),:) &
                      +CSHG(a,c,csh_lm_of(L-j,M-mm),:)*CSH(b,c,csh_lm_of(j,mm))  )  /    &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
     endif

    end subroutine calculate_factor_Ga

    subroutine calculate_factor_GGa(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm,k_gr,k2dr
      ! A_factor is a tenzor of rank L => M=-L,L
     if(a.eq.1) then
      A_factorGGa(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGa(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGa(n,LM,j,:,:) = A_factorGGa(n,LM,j,:,:) +  &
                       CSHGGa(a,c,csh_lm_of(L-j,M-mm),:,:) * &
                       CSH(b,c,csh_lm_of(j,mm) ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do

!     if(abs(A_factorGGa(n,LM,j,1,1)).gt.1.0e-10) then
!      print*,'A_facGGA n',n,A_factor(n,LM,j),A_factorGa(n,LM,j,1),A_factorGGa(n,LM,j,1,1)
!     endif
            enddo
         end do
      end do
!     print*,'sum A_facGGA n',n,sum(A_factor(n,:,:)),sum(A_factorGa(n,:,:,1)),sum(A_factorGGa(n,:,:,1,1))
     elseif(b.eq.1) then
      A_factorGGa(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGa(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGa(n,LM,j,:,:) = A_factorGGa(n,LM,j,:,:) + &
                       CSH(a,c,csh_lm_of(L-j,M-mm) )    *        &
                       CSHGGa(b,c,csh_lm_of(j,mm),:,:) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
!     if(abs(A_factorGGa(n,LM,j,1,1)).gt.1.0e-10) then
!     print*,'A_facGGA n',n,A_factor(n,LM,j),A_factorGa(n,LM,j,1),A_factorGGa(n,LM,j,1,1)
!     endif
            end do
         end do
      end do
!     print*,'sum A_facGGA n',n,sum(A_factor(n,:,:)),sum(A_factorGa(n,:,:,1)),sum(A_factorGGa(n,:,:,1,1))
     else
      A_factorGGa(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGa(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGa(n,LM,j,:,:) = A_factorGGa(n,LM,j,:,:) + &
                      (CSH(a,c,csh_lm_of(L-j,M-mm) )*CSHGGa(b,c,csh_lm_of(j,mm),:,:) &
                      +CSHGGa(a,c,csh_lm_of(L-j,M-mm),:,:)*CSH(b,c,csh_lm_of(j,mm))  )  /  &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
              do k2dr=1,3
               do k_gr=1,3
                A_factorGGa(n,LM,j,k_gr,k2dr)=A_factorGGa(n,LM,j,k_gr,k2dr) + &
                 ( CSHG(a,c,csh_lm_of(L-j,M-mm),k2dr)*CSHG(b,c,csh_lm_of(j,mm),k_gr) &
                  +CSHG(a,c,csh_lm_of(L-j,M-mm),k_gr)*CSHG(b,c,csh_lm_of(j,mm),k2dr) ) / &
                   sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm)*fact(j+mm)*fact(j-mm), r8_kind))
               enddo
              enddo
             end do
!     if(abs(A_factorGGa(n,LM,j,1,1)).gt.1.0e-10) then
!      print*,'A_facGGA n',n,A_factor(n,LM,j),A_factorGa(n,LM,j,1),A_factorGGa(n,LM,j,1,1)
!     endif
            end do
         end do
      end do
!     print*,'sum A_facGGA n',n,sum(A_factor(n,:,:)),sum(A_factorGa(n,:,:,1)),sum(A_factorGGa(n,:,:,1,1))
     endif

    end subroutine calculate_factor_GGa

    subroutine calculate_factor_GaGb(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm,k_gr,k2dr
      ! A_factor is a tenzor of rank L => M=-L,L

      ! six possible cases
      !  1 2 3
      !  2 1 3
      !  1 3 2
      !  2 3 1
      !  3 1 2
      !  3 2 1

     if(a.eq.1.and.b.eq.2) then         ! 1 2 3   - (1) of 6
      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) +  &
                       !   ^ ^                              ^ ^
                       CSHG(a,c,csh_lm_of(L-j,M-mm),k_gr) * &
                       !  ^ 1
                       CSHGb(b,c,csh_lm_of(j,mm),k2dr) /               &
                       !   ^ 2
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo
               end do
            enddo
         end do
      end do

     elseif(b.eq.1.and.a.eq.2) then   ! 2 1 3    - (2) of 6
      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) + &
                       !   ^ ^                              ^ ^
                       CSHGb(a,c,csh_lm_of(L-j,M-mm),k2dr)    *        &
                       !   ^ 2
                       CSHG(b,c,csh_lm_of(j,mm),k_gr) /               &
                       !  ^ 1
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo
               end do
            end do
         end do
      end do

     elseif(a.eq.1.and.c.eq.2) then      ! case 1 3 2  - (3) of 6
      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle

                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) + &
                      !    ^ ^        ^    ^                ^ ^        ^    ^
                      (CSHGaGb(a,c,csh_lm_of(L-j,M-mm),k_gr,k2dr )*CSH(b,c,csh_lm_of(j,mm)) &
                      !    ^ ^ 1 2
                      +CSHG(a,c,csh_lm_of(L-j,M-mm),k_gr)*CSHGb(b,c,csh_lm_of(j,mm),k2dr)   &
                      !   ^ 1                                 ^   2
                      )/  sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo

             end do
            end do
         end do
      end do
     elseif(c.eq.1.and.a.eq.2) then     !  case 2 3 1 = (4) of 6
      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle

                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) + &
                      !    ^ ^        ^    ^                ^ ^        ^    ^
                      (CSHGaGb(a,c,csh_lm_of(L-j,M-mm),k_gr,k2dr)*CSH(b,c,csh_lm_of(j,mm)) &
                      !    ^ ^ 2 1
                      +CSHGb(a,c,csh_lm_of(L-j,M-mm),k2dr)*CSHG(b,c,csh_lm_of(j,mm),k_gr)   &
                      !    ^ 2                                ^   1
                      )/  sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo

             end do
            end do
         end do
      end do
      elseif(b.eq.1.and.c.eq.2) then     ! case 3 1 2 - (5) - of 6
      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle

                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) + &
                      !    ^ ^        ^    ^                ^ ^        ^    ^
                      (CSHGb(a,c,csh_lm_of(L-j,M-mm),k2dr)*CSHG(b,c,csh_lm_of(j,mm),k_gr) &
                      !    ^   ^                              ^ 1
                      +CSH(a,c,csh_lm_of(L-j,M-mm))*CSHGaGb(b,c,csh_lm_of(j,mm),k_gr,k2dr) &
                      !                                 ^   1 2
                      )/  sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo

             end do
            end do
         end do
      end do
     elseif(c.eq.1.and.b.eq.2) then   ! case 3 2 1 - (6) of 6

      A_factorGaGb(n,:,:,:,:)=czero
      !        ^ ^
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGaGb(n,LM,j,:,:)=czero
               !        ^ ^
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle

                  do k_gr=1,3
                  do k2dr=1,3
                  A_factorGaGb(n,LM,j,k_gr,k2dr) = A_factorGaGb(n,LM,j,k_gr,k2dr) + &
                      !    ^ ^        ^    ^                ^ ^        ^    ^
                      (CSHG(a,c,csh_lm_of(L-j,M-mm),k_gr)*CSHGb(b,c,csh_lm_of(j,mm),k2dr) &
                      !   ^   ^                     ^         ^ ^                   ^
                      +CSH(a,c,csh_lm_of(L-j,M-mm))*CSHGaGb(b,c,csh_lm_of(j,mm),k_gr,k2dr) &
                      !                                 ^ ^ 2 1                ^
                      )/  sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
                  enddo
                  enddo

             end do
            end do
         end do
      end do

     endif

    end subroutine calculate_factor_GaGb

    subroutine calculate_factor_Gb(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm
      ! A_factor is a tenzor of rank L => M=-L,L
      if(a.eq.2) then
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGb(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGb(n,LM,j,:) = A_factorGb(n,LM,j,:)   +        &
                       CSHGb(a,c,csh_lm_of(L-j,M-mm),: )    *        &
                       CSH(b,c,csh_lm_of(j,mm) ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
        elseif(b.eq.2) then
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGb(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGb(n,LM,j,:) = A_factorGb(n,LM,j,:)   +        &
                       CSH(a,c,csh_lm_of(L-j,M-mm) )    *        &
                       CSHGb(b,c,csh_lm_of(j,mm),: ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
        else
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGb(n,LM,j,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGb(n,LM,j,:) = A_factorGb(n,LM,j,:)   +        &
                      (CSH(a,c,csh_lm_of(L-j,M-mm) )*CSHGb(b,c,csh_lm_of(j,mm),:) &
                      +CSHGb(a,c,csh_lm_of(L-j,M-mm),:)*CSH(b,c,csh_lm_of(j,mm))  )  /    &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
            end do
         end do
      end do
        endif

    end subroutine calculate_factor_Gb

    subroutine calculate_factor_GGb(n,a,b,c,l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: n,a,b,c,l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,j,M,mm,k_gr,k2dr
      ! A_factor is a tenzor of rank L => M=-L,L
     if(a.eq.2) then
      A_factorGGb(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGb(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGb(n,LM,j,:,:) = A_factorGGb(n,LM,j,:,:) +  &
                       CSHGGb(a,c,csh_lm_of(L-j,M-mm),:,:) * &
                       CSH(b,c,csh_lm_of(j,mm) ) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do

!     if(abs(A_factorGGb(n,LM,j,1,1)).gt.1.0e-10) then
!      print*,'A_facGGb n',n,A_factor(n,LM,j),A_factorGb(n,LM,j,1),A_factorGGb(n,LM,j,1,1)
!     endif
            enddo
         end do
      end do
!     print*,'sum A_facGGA n',n,sum(A_factor(n,:,:)),sum(A_factorGa(n,:,:,1)),sum(A_factorGGa(n,:,:,1,1))
     elseif(b.eq.2) then
      A_factorGGb(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGb(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGb(n,LM,j,:,:) = A_factorGGb(n,LM,j,:,:) + &
                       CSH(a,c,csh_lm_of(L-j,M-mm) )    *        &
                       CSHGGb(b,c,csh_lm_of(j,mm),:,:) /               &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
               end do
!     if(abs(A_factorGGa(n,LM,j,1,1)).gt.1.0e-10) then
!     print*,'A_facGGA n',n,A_factor(n,LM,j),A_factorGa(n,LM,j,1),A_factorGGa(n,LM,j,1,1)
!     endif
            end do
         end do
      end do
!     print*,'sum A_facGGA n',n,sum(A_factor(n,:,:)),sum(A_factorGa(n,:,:,1)),sum(A_factorGGa(n,:,:,1,1))

     else
      A_factorGGb(n,:,:,:,:)=czero
      do L=0,l_max
         do M=-L,L
            do j=0,L
               LM  = csh_lm_of(L,M)
               A_factorGGb(n,LM,j,:,:)=czero
               do mm=-j,j
                  if(abs(M-mm) > L-j) cycle
                  A_factorGGb(n,LM,j,:,:) = A_factorGGb(n,LM,j,:,:) + &
                      (CSH(a,c,csh_lm_of(L-j,M-mm) )*CSHGGb(b,c,csh_lm_of(j,mm),:,:) &
                      +CSHGGb(a,c,csh_lm_of(L-j,M-mm),:,:)*CSH(b,c,csh_lm_of(j,mm))  )  /  &
                       sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm) * &
                       fact(j+mm)*fact(j-mm), r8_kind) )
              do k2dr=1,3
               do k_gr=1,3
                A_factorGGb(n,LM,j,k_gr,k2dr)=A_factorGGb(n,LM,j,k_gr,k2dr) + &
                 ( CSHGb(a,c,csh_lm_of(L-j,M-mm),k2dr)*CSHGb(b,c,csh_lm_of(j,mm),k_gr) &
                  +CSHGb(a,c,csh_lm_of(L-j,M-mm),k_gr)*CSHGb(b,c,csh_lm_of(j,mm),k2dr) ) / &
                   sqrt(real(fact(L-j+M-mm)*fact(L-j-M+mm)*fact(j+mm)*fact(j-mm), r8_kind))
               enddo
              enddo
             end do
!     if(abs(A_factorGGb(n,LM,j,1,1)).gt.1.0e-10) then
!      print*,'A_facGGb n',n,A_factor(n,LM,j),A_factorGb(n,LM,j,1),A_factorGGb(n,LM,j,1,1)
!     endif
            end do
         end do
      end do
!     print*,'sum A_facGGb n',n,sum(A_factor(n,:,:)),sum(A_factorGb(n,:,:,1)),sum(A_factorGGb(n,:,:,1,1))
     endif

    end subroutine calculate_factor_GGb

    subroutine calc_csh_lm_of(l_max)
      ! Purpose : the calculation of A-factors (Eq.15)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: l_max
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: L,M
      do L=0,l_max
         do M=-L,L
               csh_lm_of(L,M)  = fcsh_lm_of(L,M)
         end do
      end do

    end subroutine calc_csh_lm_of

    subroutine CSH_lm_map (l_max, key)
      ! Purpose : lm mapping for CSH
      ! (to be moved e.g. to "solhrules_module" later)
      implicit none
      integer(kind=i4_kind), intent(in) :: l_max
      character(len=*)                  :: key
      integer(kind=i4_kind)             :: l,m,lm

      if(trim(key) == "setup") then
         allocate (csh_map % l((l_max + 1)**2), csh_map % m((l_max + 1)**2), stat=alloc_stat(33))
         if(alloc_stat(33) /= 0) &
              call error_handler("CSH_lm_map : chs_map allocation failed")
         alloc_stat(33)=1

         !
         ! Correspondence between tuples and indices:
         !
         ! (l, m) <=> fcsh_lm_of (l, m) == 1 + l**2 + l + m
         !
         lm = 0
         do l = 0, l_max
            do m = -l, l
               lm = lm + 1
               csh_map % l(lm) = l
               csh_map % m(lm) = m
               ASSERT(lm==fcsh_lm_of(l,m))
            end do
         end do
      else if(trim(key) == "close") then
         deallocate (csh_map % l, csh_map % m, stat=alloc_stat(33))
         if(alloc_stat(33) /= 0) &
              call error_handler("CSH_lm_map : chs_map deallocation failed")
      else
         call error_handler("CSH_lm_map : unknown key "//trim(key))
      end if

    end subroutine CSH_lm_map

    pure function fcsh_lm_of (l, m)
      !
      ! Purpose: Inverse mapping, returns "lm" of CSH for (l, m) pair.
      !
      implicit none
      integer (i4_kind), intent (in) :: l, m
      integer (i4_kind) :: fcsh_lm_of
      !** End of interface *****************************************

      fcsh_lm_of = 1 + l**2 + l + m
    end function fcsh_lm_of

    function CSH_scalar(lmax, arg) result(CSH)
      ! Purpose : calculation of Complex Spherical Harmonics,
      !           just a complex variant of "solid_harmonics_scalar"
      !           (to be moved to e.g. "solhrules_module")
      !           CSH(l,0) = RSH(l,1)
      !           CSH(l,m) = (-1)**m * ( RSH(l,2m) + i*RSH(l,2m+1) )/sqrt(2)
      !           m=1...l
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: lmax
      real(kind=r8_kind), intent(in) :: arg(:) ! (3)
      !** End of interface *****************************************

      complex(kind=c16_kind)               :: CSH((lmax+1)**2)
      complex(kind=c16_kind), parameter    :: i=(0.0_r8_kind,1.0_r8_kind)
      !      real(kind=r8_kind),     allocatable  :: rsh(:)
      real(kind=r8_kind)                   :: rsh((lmax+1)**2)
      real(kind=r8_kind)                   :: sqrt_two
      integer(kind=i4_kind)                :: l,m,lm,lm_csh

      ASSERT(size(arg)==3)

      sqrt_two = sqrt(2.0_r8_kind)
      rsh = solid_harmonics_scalar(lmax,arg)
      do l=0,lmax
         lm_csh = l*(l+1)+1
         CSH(lm_csh) = RSH(l**2+1)
         lm = l**2 + 2
         do m=1,l
            lm_csh = lm_csh + 1
            CSH( lm_csh ) = (-1)**m * ( RSH(lm) + i*RSH(lm+1) ) / sqrt_two
            CSH( lm_csh - 2*m ) = (-1)**m * conjg( CSH(lm_csh) )
            lm = lm + 2
         end do
      end do

    end function CSH_scalar

    function CSHG_scalar(lmax, arg, i_p) result(CSHG)
      ! Purpose : calculation of Complex Spherical Harmonics,
      !           just a complex variant of "solid_harmonics_scalar"
      !           (to be moved to e.g. "solhrules_module")
      !           CSH(l,0) = RSH(l,1)
      !           CSH(l,m) = (-1)**m * ( RSH(l,2m) + i*RSH(l,2m+1) )/sqrt(2)
      !           m=1...l
      use solhrules_module
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: lmax
      real(kind=r8_kind), intent(in) :: arg(:) ! (3)
      integer(kind=i4_kind), intent(in), optional :: i_p
      !** End of interface *****************************************

      complex(kind=c16_kind)               :: CSHG((lmax+1)**2,3)
      complex(kind=c16_kind), parameter    :: i=(0.0_r8_kind,1.0_r8_kind)
      !      real(kind=r8_kind),     allocatable  :: rsh(:) ! to be calculated once in CSH_scalar
      real(kind=r8_kind)                   :: rsh((lmax+1)**2)
      real(kind=r8_kind)                   :: rshg((lmax+1)**2,3)
      real(kind=r8_kind)                   :: sqrt_two
      integer(kind=i4_kind)                :: l,m,lm,lm_csh,k_gr
      integer(kind=i4_kind)                :: xyz_map(3)=(/3,4,2/)

      ASSERT(size(arg)==3)

      sqrt_two = sqrt(2.0_r8_kind)
      rsh = solid_harmonics_scalar(lmax,arg) ! to be precalculated in CSH_scalar
      rshg = 0
      do l=1,(lmax+1)**2
        do k_gr=1,3
          do k=1,solhrules_differential(xyz_map(k_gr),l)%n_summands
           rshg(l,k_gr)=rshg(l,k_gr)+solhrules_differential(xyz_map(k_gr),l)%coef(k) &
                       *rsh(solhrules_differential(xyz_map(k_gr),l)%lm_sh(k))
          enddo
        enddo
      enddo

   do k_gr=1,3
      do l=0,lmax
         lm_csh = l*(l+1)+1
         CSHg(lm_csh,k_gr) = RSHg(l**2+1,k_gr)
         lm = l**2 + 2
         do m=1,l
           lm_csh = lm_csh + 1
            CSHg( lm_csh,k_gr ) = (-1)**m * ( RSHg(lm,k_gr) + i*RSHg(lm+1,k_gr) ) / sqrt_two
            CSHg( lm_csh - 2*m, k_gr ) = (-1)**m * conjg( CSHg(lm_csh, k_gr) )
           lm = lm + 2
         end do
      end do
    enddo

    end function CSHG_scalar

   function CSHGG_scalar(lmax, arg, i_p) result(CSHGG)
      ! Purpose : calculation of Complex Spherical Harmonics,
      !           just a complex variant of "solid_harmonics_scalar"
      !           (to be moved to e.g. "solhrules_module")
      !           CSH(l,0) = RSH(l,1)
      !           CSH(l,m) = (-1)**m * ( RSH(l,2m) + i*RSH(l,2m+1) )/sqrt(2)
      !           m=1...l
      use solhrules_module
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: lmax
      real(kind=r8_kind), intent(in) :: arg(:) ! (3)
      integer(kind=i4_kind), intent(in), optional :: i_p
      !** End of interface *****************************************

      complex(kind=c16_kind)               :: CSHGG((lmax+1)**2,3,3)
      complex(kind=c16_kind), parameter    :: i=(0.0_r8_kind,1.0_r8_kind)
      real(kind=r8_kind)                   :: sqrt_two
      integer(kind=i4_kind)                :: l,m,lm,lm_csh,k_gr,k2dr
      integer(kind=i4_kind)                :: xyz_map(3)=(/3,4,2/)

      ASSERT(size(arg)==3)

      sqrt_two = sqrt(2.0_r8_kind)
      rsh = solid_harmonics_scalar(lmax,arg) ! to be precalculated in CSH_scalar
      rshg = 0.0_r8_kind
      do l=1,(lmax+1)**2
        do k_gr=1,3
         do k=1,solhrules_differential(xyz_map(k_gr),l)%n_summands
          rshg(l,k_gr)=rshg(l,k_gr)+solhrules_differential(xyz_map(k_gr),l)%coef(k) &
               *rsh(solhrules_differential(xyz_map(k_gr),l)%lm_sh(k))
          enddo
        enddo
      enddo

      rshgg = 0.0_r8_kind
      do l=1,(lmax+1)**2
       do k2dr=1,3
        do k_gr=1,3
         do k=1,solhrules_differential(xyz_map(k2dr),l)%n_summands
          rshgg(l,k_gr,k2dr)=rshgg(l,k_gr,k2dr) &
                            +solhrules_differential(xyz_map(k2dr),l)%coef(k) &
               *rshg(solhrules_differential(xyz_map(k2dr),l)%lm_sh(k),k_gr)
          enddo
        enddo
       enddo
      enddo

       if(present(i_p)) then
       print*,i_p,uni_c,sum(rshg(:,1)),sum(rshgg(:,1,1))
       if(i_p.eq.uni_c) then
        do l=1,(lmax+1)**2
         print*,rshg(l,1),rshgg(l,1,1), '1 rshg rshgg',l
!         print*,rshg(l,2),rshgg(l,2,1), '2 rshg rshgg',l
!         print*,rshg(l,3),rshgg(l,3,1), '3 rshg rshgg',l
        enddo
       endif
       endif



  do k2dr=1,3
   do k_gr=1,3
      do l=0,lmax
         lm_csh = l*(l+1)+1
         CSHgg(lm_csh,k_gr,k2dr) = RSHgg(l**2+1,k_gr,k2dr)
         lm = l**2 + 2
         do m=1,l
           lm_csh = lm_csh + 1
            CSHgg(lm_csh,k_gr,k2dr)=(-1)**m * (RSHgg(lm,k_gr,k2dr)+i*RSHgg(lm+1,k_gr,k2dr))/sqrt_two
            CSHgg(lm_csh - 2*m, k_gr, k2dr) = (-1)**m * conjg( CSHgg(lm_csh, k_gr, k2dr) )
           lm = lm + 2
         end do
      end do
    enddo
  enddo
  end function CSHGG_scalar

    function bin (n1, n2)
      ! Binomial coefficient | n1 |
      !                      | n2 |
      implicit none
      integer (i4_kind), intent (in) :: n1, n2
      integer (i8_kind) :: bin
      !** End of interface *****************************************

      if (n1 >= n2 .and. n2 >=0) then
         bin = fact(n1) / (fact(n2) * fact(n1 - n2))
      else
         bin = 0
      end if
    end function bin

    subroutine calc_and_pack_3j()
      ! Puspose : Calculate and pack 3j-symbols products
      ! (to be moved outside "calc_3center" to setup only once)
      ! three(j1,m1,j2,m2,j3,j3) => three(jm1,jm2,jm3) !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      !** End of interface *****************************************
      real(kind=r8_kind), allocatable :: three(:,:,:,:,:,:)
      integer(kind=i4_kind)  :: j1,j2,j3,m1,m2,m3,N_nonzero,i_term

      allocate(three(0:l_max,-l_max:l_max, &
                     0:l_max,-l_max:l_max, &
                     0:l_max,-l_max:l_max),&
                     stat=alloc_stat(31))
      ASSERT(alloc_stat(31).eq.0)
      alloc_stat(31)=1
      MEMLOG(size(three))
      three = zero
      do j1 = 0, l_max
         do j2 = 0, l_max
            do j3 = 0, l_max
               do m1 = -j1, j1
                  do m2 = -j2, j2
                     do m3 = -j3, j3
                        three(j1,m1,j2,m2,j3,m3) = wigner_3j(j1,j2,j3, m1, m2, m3) * &
                                                   wigner_3j(j1,j2,j3, 0 , 0 , 0 )
                     end do
                  end do
               end do
            end do
         end do
      end do
      !**************************! Pack 3j-symbols !**************************!
      allocate( jj_coeff(0:l_max,0:l_max,0:l_max), stat=alloc_stat(31))
      MEMLOG(size(jj_coeff))
      ASSERT(alloc_stat(31).eq.0)
      maxNum_M_terms=0
      do j1=0,l_max
         do j2=0,l_max
            do j3=0,l_max
            N_nonzero = count(abs(three(j1,:,j2,:,j3,:)) > 1.0e-9_r8_kind)
            maxNum_M_terms=max(maxNum_M_terms,N_nonzero)
            enddo
         enddo
      enddo
        allocate(m1jj(0:l_max,0:l_max,0:l_max,maxNum_M_terms), &
                 m2jj(0:l_max,0:l_max,0:l_max,maxNum_M_terms), &
                 m3jj(0:l_max,0:l_max,0:l_max,maxNum_M_terms), &
                 valjj(0:l_max,0:l_max,0:l_max,maxNum_M_terms),stat=alloc_stat(31))
        ASSERT(alloc_stat(31).eq.0)
        alloc_stat(31)=1
        MEMLOG(size(m1jj)*4)

      do j1=0,l_max
         do j2=0,l_max
            do j3=0,l_max
               N_nonzero = count(abs(three(j1,:,j2,:,j3,:)) > 1.0e-9_r8_kind)
               jj_coeff(j1,j2,j3)%Num_M_terms = N_nonzero
               if(N_nonzero > 0) then
                  allocate(jj_coeff(j1,j2,j3)%value(N_nonzero), stat=alloc_stat(31))
                  ASSERT(alloc_stat(31).eq.0)
                  alloc_stat(31)=1
                  MEMLOG(size(jj_coeff(j1,j2,j3)%value)*4)
                  allocate( jj_coeff(j1,j2,j3)%m1(N_nonzero), &
                            jj_coeff(j1,j2,j3)%m2(N_nonzero), &
                            jj_coeff(j1,j2,j3)%m3(N_nonzero), stat=alloc_stat(31))
                  ASSERT(alloc_stat(31).eq.0)
                  i_term = 0
                  do m1 = -j1, j1
                     do m2 = -j2, j2
                        do m3 = -j3, j3
                           if( abs(three(j1,m1,j2,m2,j3,m3)) > 1.0e-9_r8_kind) then
                              i_term = i_term + 1
                              if(i_term > jj_coeff(j1,j2,j3)%Num_M_terms) &
                                   call error_handler("calc_and_pack_3j : Num_M_terms exceeded ?")
                              jj_coeff(j1,j2,j3)%m1(i_term)    = m1
                              jj_coeff(j1,j2,j3)%m2(i_term)    = m2
                              jj_coeff(j1,j2,j3)%m3(i_term)    = m3
                              jj_coeff(j1,j2,j3)%value(i_term) = three(j1,m1,j2,m2,j3,m3)

                              m1jj(j1,j2,j3,i_term)=m1
                              m2jj(j1,j2,j3,i_term)=m2
                              m3jj(j1,j2,j3,i_term)=m3
                              valjj(j1,j2,j3,i_term)=three(j1,m1,j2,m2,j3,m3)
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      MEMLOG(-size(three))
      deallocate(three, stat=alloc_stat(31))
      ASSERT(alloc_stat(31).eq.0)

    end subroutine calc_and_pack_3j

    subroutine calc_angular_cmplx(L_a,L_b,L_c,tp)
      ! Purpose : calculation of complex angular factor
      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      logical, optional,intent(in) :: tp
      !** End of interface *****************************************
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind), pointer    :: mp1,mp2,mp3

      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               call J_Complex("allocate",Ang_C(Ma,Mb,Mc)%J,L_a,L_b,L_c,"Ang_C%J")
               sqrt_LM = sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
                if(.true.) then
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                        if(.NOT. even_triangle(j1,j2,j3) ) cycle
                        d_fact = df(j1 + j2 + j3 + 1)
                        do i1 = 0,L_a-j1
                           do i2 = 0,L_b-j2
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                 do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                                    mp1 => jj_coeff(j1,j2,j3)%m1(i_term)
                                    mp2 => jj_coeff(j1,j2,j3)%m2(i_term)
                                    mp3 => jj_coeff(j1,j2,j3)%m3(i_term)
                                    if( abs(Ma-mp1) > L_a-j1 .or. &
                                        abs(Mb-mp2) > L_b-j2 .or. &
                                        abs(Mc-mp3) > L_c-j3      ) cycle
                                    !write(*,*) "A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) ", A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1)
                                    !write(*,*) "A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) ", A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                                    !write(*,*) "A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3) ", A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                    ! write(*,*) " "
                                    Ang_C(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3) =                    &
                                         Ang_C(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3) +               &
                                             d_fact *                                       &
                                             jj_coeff(j1,j2,j3)%value(i_term) *             &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *      &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                             A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3) *      &
                                             sqrt_LM                      /                 &
                                             sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                        fact(j2+mp2)*fact(j2-mp2) *         &
                                                        fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                                    !write(*,"('Ang ',5x,4i4,4x,2(d12.6,d12.6))") j1,j2,i1,i2, Ang_C(ma,mb,0)%J(j1,j2,0,i1,i2,0)
                                    !write(*,*) "Ang ", Ang_C(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3)
                                    !write(*,*) " "
                                 end do
                                 !**********************************************************!
                if(present(tp)) then
                             if(tp)    write(*,"('Ang ',3i2,5x,4i4,4x,2(d12.6,d12.6))") &
                                      ma,mb,mc,j1,j2,i1,i2, Ang_C(ma,mb,0)%J(j1,j2,0,i1,i2,0)
                endif
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
               if(Ma>0) Ang_C(Ma,Mb,Mc)%J =  Ang_C(Ma,Mb,Mc)%J*(-1)**Ma
               if(Mb>0) Ang_C(Ma,Mb,Mc)%J =  Ang_C(Ma,Mb,Mc)%J*(-1)**Mb
               if(Mc>0) Ang_C(Ma,Mb,Mc)%J =  Ang_C(Ma,Mb,Mc)%J*(-1)**Mc
                endif
            end do
         end do
      end do
    end subroutine calc_angular_cmplx

    subroutine calc_angular_cmplxP(L_a,L_b,L_c)
      ! Purpose : calculation of complex angular factor
      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind), pointer    :: mp1,mp2,mp3
      integer(kind=i4_kind)             :: jj

      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                        if(.NOT. even_triangle(j1,j2,j3) ) cycle
                        d_fact = df(j1 + j2 + j3 + 1)
                        do i1 = 0,L_a-j1
                           do i2 = 0,L_b-j2
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                 jj=jj+1
                                 do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                                    mp1 => jj_coeff(j1,j2,j3)%m1(i_term)
                                    mp2 => jj_coeff(j1,j2,j3)%m2(i_term)
                                    mp3 => jj_coeff(j1,j2,j3)%m3(i_term)
                                    if( abs(Ma-mp1) > L_a-j1 .or. &
                                        abs(Mb-mp2) > L_b-j2 .or. &
                                        abs(Mc-mp3) > L_c-j3      ) cycle
                                    AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) +             &
                                             d_fact *                                       &
                                             jj_coeff(j1,j2,j3)%value(i_term) *             &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *      &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                             A_factor(5,csh_lm_of(L_c-j3,Mc-mp3),i3) *      &
                                             sqrt_LM                      /                 &
                                             sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                        fact(j2+mp2)*fact(j2-mp2) *         &
                                                        fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
               if(Ma>0) AngCMat(:,Ma,Mb,Mc) =  AngCMat(:,Ma,Mb,Mc)*(-1)**Ma
               if(Mb>0) AngCMat(:,Ma,Mb,Mc) =  AngCMat(:,Ma,Mb,Mc)*(-1)**Mb
               if(Mc>0) AngCMat(:,Ma,Mb,Mc) =  AngCMat(:,Ma,Mb,Mc)*(-1)**Mc
            end do
         end do
      end do
    end subroutine calc_angular_cmplxP

    subroutine calc_angular_cmplx3cPP(L_a,L_b,L_c,AngCMat,n_jj,n)
      ! Purpose : calculation of complex angular factor
      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      integer(kind=i4_kind), intent(in) :: n,n_jj
      complex(kind=c16_kind),intent(inout) :: AngCMat(n_jj,-L_a:L_a,-L_b:L_b,-L_c:L_c)
      !** End of interface *****************************************
      complex(kind=c16_kind)            :: z1_fact,z2_fact
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind)             :: mp1,mp2,mp3
      integer(kind=i4_kind)             :: jj,jj_off

    if(.false.) then
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                  AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
    else
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
     endif
    end subroutine calc_angular_cmplx3cPP

    subroutine calc_angular_cmplx3cGa(L_a,L_b,L_c,n_jj,n,dervs,AngCMatGb)
      ! Purpose : calculation of complex angular factor
      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!

      ! hystory:  now extednded with GaGb



      implicit none
      !------------ Declaration of formal parameters ---------------
      complex(kind=c16_kind),optional,intent(inout) :: &
                              AngCMatGb(n_jj,-L_a:L_a,-L_b:L_b,-L_c:L_c,3)
      logical, intent(in) :: dervs
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      integer(kind=i4_kind), intent(in) :: n,n_jj
      !** End of interface *****************************************
      complex(kind=c16_kind)            :: z1_fact,z2_fact,A23
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind)             :: mp1,mp2,mp3
      integer(kind=i4_kind)             :: jj,jj_off
      integer(kind=i4_kind)             ::k_gr,k2dr

   dr2dr1els: if(dervs) then
    if(present(AngCMatGb)) then ! present in calculations of coul integrals
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                       valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                       jj=jj_off
                        do i1 = 0,L_a-j1
                        z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                                      !^ now only 3

                                  A23=A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)

                                   AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      & !A23
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        & !A23
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )

                                   ! first A_factorGa differintiated

                                   AngCMatGGa(jj,Ma,Mb,Mc,:,:)=AngCMatGGa(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (A23*A_factorGGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   )

                                  do k2dr=1,3
                                   do k_gr=1,3
                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
                                             A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr) *      &
                                               A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        + &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                               A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   )

                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
                                             A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )
                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
                                             A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )
                                   enddo
                                  enddo


                                  AngCMatGb(jj,Ma,Mb,Mc,:)=AngCMatGb(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:)    &
                                             *A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        &
                                             *A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   &
                                                      )

                                   ! GGb derivatives first A_factorGa differintiated 3->n

                                   AngCMatGGb(jj,Ma,Mb,Mc,:,:)=AngCMatGGb(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (  A23*A_factorGGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   &
                                             )
                                k2drvs_b:  do k2dr=1,3
                                   do k_gr=1,3
                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
                                             A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)       &
                                              *A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        + &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                               A_factorGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   &
                                                                                              )

                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )

                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )

                                   enddo
                                  enddo k2drvs_b

                                   ! GaGb derivatives first A_factorGa differintiated

                                   AngCMatGaGb(jj,Ma,Mb,Mc,:,:)=AngCMatGaGb(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (A23*A_factorGaGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                            +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGaGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGaGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   )

                                gagb_k2drvs:  do k2dr=1,3
                                   do k_gr=1,3

                                    AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
          A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)*A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)+ &
          A_factor  (2,csh_lm_of(L_b-j2,Mb-mp2),i2) *    A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   )
                                             !-------------------------------------------------

                                    AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
          A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
          A_factor  (1,csh_lm_of(L_a-j1,Ma-mp1),i1) *      A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )
                                               !-----------------------------------

                                   AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
          A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) *A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
          A_factor  (1,csh_lm_of(L_a-j1,Ma-mp1),i1) *    A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )
                                               !--------------------------------------------------

                                   enddo
                                  enddo gagb_k2drvs

                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do

    else ! i.e. not. present Gb

      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                              jj=jj+1
                              AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                  AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                              end do
                           end do
                        end do
                       endif
                      enddo
                       jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
     endif

   elseif(integralpar_gradients) then
    if(present(AngCMatGb)) then
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                  AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                                  AngCMatGb(jj,Ma,Mb,Mc,:)=AngCMatGb(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGb(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
    else
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                              jj=jj+1
                              AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                  AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                              d_fact * (                                     &
                                              A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) *  &
                                              A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(n,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                              end do
                           end do
                        end do
                       endif
                      enddo
                       jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
     endif

    else dr2dr1els
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(n,csh_lm_of(L_c-j3,Mc-mp3),i3)
                              end do
                           end do
                        end do
                       endif
                      enddo
                    jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
            end do
         end do
      end do
    endif dr2dr1els
    end subroutine calc_angular_cmplx3cGa

    subroutine calc_angular_cmplx3cSA(L_a,L_b,L_c,AngCMat,n_jj)

      ! Purpose : calculation of complex angular factor if l_c gt 0

      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      integer(kind=i4_kind), intent(in) :: n_jj
      complex(kind=c16_kind),intent(inout) :: AngCMat(n_jj,-L_a:L_a,-L_b:L_b,-L_c:L_c)
      !** End of interface *****************************************
      complex(kind=c16_kind)            :: z1_fact,z2_fact,A23
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind)             :: mp1,mp2,mp3
      integer(kind=i4_kind)             :: jj,jj_off,k_gr,k2dr

     dr2dr1els: if(integralpar_dervs) then
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
            if(m_c_required(Ma,Mb,Mc,j)) then
               sqrt_LM = afac_sig(ma,mb,mc)* &
                         sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                       mp1 =m1jj(j1,j2,j3,i_term)
                       mp2 =m2jj(j1,j2,j3,i_term)
                       mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                        z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)

                                  A23=A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)

                                   AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                   d_fact * (A23*A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )

                                   ! first A_factorGa differintiated

                                   AngCMatGGa(jj,Ma,Mb,Mc,:,:)=AngCMatGGa(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (A23*A_factorGGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   )
                                  do k2dr=1,3
                                   do k_gr=1,3
                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &

                                    d_fact * A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
                                             A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr) *      &
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        + &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                               A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   )

                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
                                             A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )
                                    AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGa(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
                                             A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )
                                   enddo
                                  enddo

                                   AngCMatGb(jj,Ma,Mb,Mc,:)=AngCMatGb(jj,Ma,Mb,Mc,:) +        &
                                   d_fact * (A23*A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) &

                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &

                                              A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )

                                   ! GGb derivatives first A_factorGb differintiated

                                   AngCMatGGb(jj,Ma,Mb,Mc,:,:)=AngCMatGGb(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (A23*A_factorGGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   )
                                k2drvs:  do k2dr=1,3
                                   do k_gr=1,3
                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
                                             A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr) *      &
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        + &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                               A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   )

                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )
                                    AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    d_fact * A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )

                                   enddo
                                  enddo k2drvs

                                   ! GaGb derivatives first A_factorGa differintiated

                                   AngCMatGaGb(jj,Ma,Mb,Mc,:,:)=AngCMatGaGb(jj,Ma,Mb,Mc,:,:) + &
                                   d_fact * (A23*A_factorGaGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:,:) &
                                            !             ^ ^
                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &
                                              A_factorGaGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:,:) *  &
                                            !          ^ ^
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGaGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:,:))   )
                                            !          ^ ^
                                gagb_k2drvs:  do k2dr=1,3
                                   do k_gr=1,3

                                    AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    !       ^ ^                                ^ ^
                                    d_fact * A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k_gr)* ( &
                                             !        ^
                                             A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr) *      &
                                             !        ^
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        + &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                               A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr)   )
                                             !          ^
                                             !-------------------------------------------------

                                    AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    !       ^ ^                                ^ ^
                                    d_fact * A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k_gr) * ( &
                                             !        ^
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               !      ^
                                               A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)     + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k2dr) )
                                               !        ^
                                               !-----------------------------------

                                    AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr)=AngCMatGaGb(jj,Ma,Mb,Mc,k_gr,k2dr) + &
                                    !       ^ ^                                ^ ^
                                    d_fact * A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,k_gr) * ( &
                                             !        ^
                                             A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,k2dr) * &
                                               !      ^
                                               A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)        + &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) * &
                                               A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,k2dr)   )
                                               !        ^
                                               !--------------------------------------------------

                                   enddo
                                  enddo gagb_k2drvs

                              end do
                           end do
                        end do
                       endif
                      enddo
                     jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
             endif
            end do
         end do
      end do
!      print*,'sum abs AngCMatGGb',sum(abs(AngCMatGGb))

     elseif(integralpar_gradients) then
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
            if(m_c_required(Ma,Mb,Mc,j)) then
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then

                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * &
                                     valjj(j1,j2,j3,i_term) /&
                                     sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact*A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *      &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)

                                  A23=A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                  AngCMatGa(jj,Ma,Mb,Mc,:)=AngCMatGa(jj,Ma,Mb,Mc,:) +        &
                                   d_fact * (A23*A_factorGa(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) &

                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &

                                              A_factorGa(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGa(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )

                                  AngCMatGb(jj,Ma,Mb,Mc,:)=AngCMatGb(jj,Ma,Mb,Mc,:) +        &
                                   d_fact * (A23*A_factorGb(1,csh_lm_of(L_a-j1,Ma-mp1),i1,:) &

                                             +A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *  (   &

                                              A_factorGb(2,csh_lm_of(L_b-j2,Mb-mp2),i2,:) *  &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)        &
                                             +A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                              A_factorGb(3,csh_lm_of(L_c-j3,Mc-mp3),i3,:))   )
                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do
             endif
            end do
         end do
      end do

     else dr2dr1els
      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
             if(m_c_required(Ma,Mb,Mc,j)) then
               sqrt_LM = afac_sig(ma,mb,mc)*sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
               jj_off=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      if(.NOT. even_triangle(j1,j2,j3) ) cycle
                      do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                        mp1 =m1jj(j1,j2,j3,i_term)
                        mp2 =m2jj(j1,j2,j3,i_term)
                        mp3 =m3jj(j1,j2,j3,i_term)
                       if(.not.(abs(Ma-mp1) > L_a-j1 .or. &
                           abs(Mb-mp2) > L_b-j2 .or. &
                           abs(Mc-mp3) > L_c-j3)     ) then
                        d_fact = df(j1 + j2 + j3 + 1) * sqrt_LM * valjj(j1, j2, j3, i_term) &
                               / sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                         fact(j2+mp2)*fact(j2-mp2) *         &
                                                         fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                        jj=jj_off
                        do i1 = 0,L_a-j1
                         z1_fact = d_fact * A_factor(1, csh_lm_of(L_a - j1, Ma - mp1), i1)
                           do i2 = 0,L_b-j2
                            z2_fact=z1_fact *  A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                  jj=jj+1
                                  AngCMat(jj,Ma,Mb,Mc)=AngCMat(jj,Ma,Mb,Mc) + z2_fact *       &
                                              A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)
                              end do
                           end do
                        end do
                       endif
                      enddo
                      jj_off=jj_off+(L_c-j3+1)*(L_b-j2+1)*(L_a-j1+1)

                     end do
                  end do
               end do

            endif
            end do
         end do
      end do
    endif dr2dr1els
    end subroutine calc_angular_cmplx3cSA

    function n_jj_terms(L_a,L_b,L_c)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind) :: n_jj_terms
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: j1,j2,j3,i1,i2,i3
      integer(kind=i4_kind)             :: jj,Lam,m

    alph_pw_max=-99
    beta_pw_max=-99
    fact2_pw_max=-99

               jj=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                      Lam = (j1 + j2 + j3)/2
                        if(.NOT. even_triangle(j1,j2,j3) ) cycle
                        do i1 = 0,L_a-j1
                           do i2 = 0,L_b-j2
                              do i3 = 0,L_c-j3
                                jj=jj+1
                               do m=0,L_a-i1+i2+j2-Lam
                                bin_fac(L_a - i1 + i2 + j2 - Lam, m) = &
                                     bin (L_a - i1 + i2 + j2 - Lam, m) * &
                                     (-1)**(L_a - i1 + i2 + j2 - Lam - m)
                               enddo
                         alph_pw_max=max(L_c-i3+j1+i1-Lam,alph_pw_max)
                         beta_pw_max=max(L_b-i2+j3+i3-Lam,beta_pw_max)
                         fact2_pw_max=max(L_a-i1+i2+j2-Lam,fact2_pw_max)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
        n_jj_terms=jj
    end function n_jj_terms

    subroutine calc_angular_cmplxGc(L_a,L_b,L_c,tp)
      ! Purpose : calculation of complex angular factor
      ! On definition : Ang(LaMa, LbMb, LcMc)%J(j1,j2,j3,j1',j2',j3')
      ! (LaMa, LbMb) are fixed
      ! (M1 M2 M3) -> (-M1 -M2 -M3) symmetry IS NOT YET USED !!!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      logical,optional,intent(in):: tp
      !** End of interface *****************************************
      real(kind=r8_kind)                :: sqrt_LM
      integer(kind=i4_kind)             :: Ma,Mb,Mc,j1,j2,j3,i1,i2,i3
      integer (i4_kind) :: i_term
      real (r8_kind) :: d_fact
      integer(kind=i4_kind), pointer    :: mp1,mp2,mp3

      do Ma=-L_a,L_a
         do Mb=-L_b,L_b
            do Mc  = -L_c,L_c
               call J_Complex("allocate",Ang_Cgc(Ma,Mb,Mc)%J,L_a,L_b,L_c,"Ang_Cgc%J")
               sqrt_LM = sqrt(real( fact(L_a+Ma)*fact(L_a-Ma) * &
                                    fact(L_b+Mb)*fact(L_b-Mb) * &
                                    fact(L_c+Mc)*fact(L_c-Mc),r8_kind) )
!               jj=0
               do j1 = 0,L_a
                  do j2 = 0,L_b
                     do j3 = 0,L_c
                        if(.NOT. even_triangle(j1,j2,j3) ) cycle
                        d_fact = df(j1 + j2 + j3 + 1)
                        do i1 = 0,L_a-j1
                           do i2 = 0,L_b-j2
                              do i3 = 0,L_c-j3
                                 !*********************************************************!
                                 do i_term = 1,jj_coeff(j1,j2,j3)%Num_M_terms
                                    mp1 => jj_coeff(j1,j2,j3)%m1(i_term)
                                    mp2 => jj_coeff(j1,j2,j3)%m2(i_term)
                                    mp3 => jj_coeff(j1,j2,j3)%m3(i_term)
                                    if( abs(Ma-mp1) > L_a-j1 .or. &
                                        abs(Mb-mp2) > L_b-j2 .or. &
                                        abs(Mc-mp3) > L_c-j3      ) cycle
                                    !write(*,*) "A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) ", A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1)
                                    !write(*,*) "A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) ", A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2)
                                    !write(*,*) "A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3) ", A_factor(3,csh_lm_of(L_c-j3,Mc-mp3),i3)
                                    ! write(*,*) " "
                                    Ang_Cgc(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3) =                    &
                                         Ang_Cgc(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3) +               &
                                             d_fact *                                       &
                                             jj_coeff(j1,j2,j3)%value(i_term) *             &
                                             A_factor(1,csh_lm_of(L_a-j1,Ma-mp1),i1) *      &
                                             A_factor(2,csh_lm_of(L_b-j2,Mb-mp2),i2) *      &
                                             A_factor(4,csh_lm_of(L_c-j3,Mc-mp3),i3) *      &
                                             sqrt_LM                      /                 &
                                             sqrt(real( fact(j1+mp1)*fact(j1-mp1) *         &
                                                        fact(j2+mp2)*fact(j2-mp2) *         &
                                                        fact(j3+mp3)*fact(j3-mp3), r8_kind ))
                                    !write(*,"('Ang ',5x,4i4,4x,2(d12.6,d12.6))") j1,j2,i1,i2, Ang_C(ma,mb,0)%J(j1,j2,0,i1,i2,0)
                                    !write(*,*) "Ang ", Ang_C(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3)
                                    !write(*,*) " "
                                 end do
!                               jj=jj+1
!                               AngCMat(jj,Ma,Mb,Mc)=Ang_Cgc(Ma,Mb,Mc)%J(j1,j2,j3,i1,i2,i3)
                                 !**********************************************************!
                                 !write(*,"('Ang ',3i2,5x,4i4,4x,2(d12.6,d12.6))") &
                                 !     ma,mb,mc,j1,j2,i1,i2, Ang_C(ma,mb,0)%J(j1,j2,0,i1,i2,0)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
               if(Ma>0) Ang_Cgc(Ma,Mb,Mc)%J =  Ang_Cgc(Ma,Mb,Mc)%J*(-1)**Ma
               if(Mb>0) Ang_Cgc(Ma,Mb,Mc)%J =  Ang_Cgc(Ma,Mb,Mc)%J*(-1)**Mb
               if(Mc>0) Ang_Cgc(Ma,Mb,Mc)%J =  Ang_Cgc(Ma,Mb,Mc)%J*(-1)**Mc
!               if(Ma>0) AngCMat(jj,Ma,Mb,Mc) =  AngCMat(jj,Ma,Mb,Mc)*(-1)**Ma
!               if(Mb>0) AngCMat(jj,Ma,Mb,Mc) =  AngCMat(jj,Ma,Mb,Mc)*(-1)**Mb
!               if(Mc>0) AngCMat(jj,Ma,Mb,Mc) =  AngCMat(jj,Ma,Mb,Mc)*(-1)**Mc
            end do
         end do
      end do
!       if(present(tp)) then
!        if(tp) print*,sum(AngCMat(:,0,0,1)),'AngCMat m=1',jj,n_jj
!       endif
    end subroutine calc_angular_cmplxGc
    !**************************************************************

    !**************************************************************
    function f_even_triangle(j1,j2,j3)
      ! Purpose : cehck if (j1,j2,j3) are even triangle
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: j1,j2,j3
      !** End of interface *****************************************
      logical                           :: f_even_triangle

      f_even_triangle = j1<=j2+j3 .and. j2<=j1+j3 .and. j3<=j1+j2 .and. &
                                           (j1+j2+j3)/2*2 ==j1+j2+j3
    end function f_even_triangle

    !**************************************************************
    function wigner_3j(j1,j2,j3,m1,m2,m3)
      ! Pupose : Wigner 3j-symbol  | j1 j2 j3 |
      !                            | m1 m2 m3 |
      ! (integer "j" and "m" only)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: j1,j2,j3,m1,m2,m3
      !** End of interface *****************************************
      real(kind=r8_kind)                :: wigner_3j
      integer(kind=i4_kind)             :: z,z_min,z_max
      logical                           :: triangle_j

      ASSERT(abs(m1)<=j1)
      ASSERT(abs(m2)<=j2)
      ASSERT(abs(m3)<=j3)

      wigner_3j = 0.0_r8_kind
      triangle_j = j1+j2 >= j3 .AND. j2+j3>=j1 .AND. j3+j1>=j2
      if( (m1+m2+m3) /= 0 .or. .NOT.triangle_j ) return
      z_max = minval( (/j1+j2-j3, j1-m1, j2+m2/) )
      z_min = abs( minval( (/0, j3-j2+m1, j3-j1-m2/) ) )
      if(z_min > z_max) return

      ! make sure not to get out of array bounds:
#     define OK(x) ((x)<=MAXFACN.and.(x)>=-1)

ASSERT(OK(z_min))
ASSERT(OK(z_max))

      ! make sure to use floating-point division and avoid integer overflows:
#     define FACT(n) real(fact(n),r8_kind)

      do z = z_min,z_max

ASSERT(OK(j1+j2-j3-z))
ASSERT(OK(j1-m1-z))
ASSERT(OK(j2+m2-z))
ASSERT(OK(j3-j2+m1+z))
ASSERT(OK(j3-j1-m2+z))

         ! make sure to use floating-point division and avoid integer overflows:
         wigner_3j = wigner_3j + (-1)**z                                     &
                     / ( FACT(z)       * FACT(j1+j2-j3-z) * FACT(j1-m1-z)*   &
                         FACT(j2+m2-z) * FACT(j3-j2+m1+z) * FACT(j3-j1-m2+z) )
      end do

ASSERT(OK(j1+j2-j3))
ASSERT(OK(j1-j2+j3))
ASSERT(OK(-j1+j2+j3))
ASSERT(OK(j1+j2+j3+1))
ASSERT(OK(j1-m1))
ASSERT(OK(j1+m1))
ASSERT(OK(j2-m2))
ASSERT(OK(j2+m2))
ASSERT(OK(j3-m3))
ASSERT(OK(j3+m3))

      ! make sure to use floating-point division and avoid integer overflows:
      wigner_3j = wigner_3j *                                              &
                 (-1)**(j1-j2+m3) *                                        &
                 sqrt(   FACT(j1+j2-j3) * FACT(j1-j2+j3) * FACT(-j1+j2+j3) &
                       / FACT(j1+j2+j3+1)                                  &
                       * FACT(j1-m1) * FACT(j1+m1)                         &
                       * FACT(j2-m2) * FACT(j2+m2)                         &
                       * FACT(j3-m3) * FACT(j3+m3)                         )
#     undef FACT
#     undef OK
    end function wigner_3j
    !****************************************************************!

    !****************************************************************!
    subroutine set_triplets()
      ! Purpose : calculation of CSH -> RSH conversion coefficients
      implicit none
      !** End of interface *****************************************
      integer(kind=i4_kind)               :: n_cr(-1:1)
      complex(kind=c16_kind)              :: cf(-1:1,2)
      complex(kind=c16_kind), allocatable :: Coeff(:)
      complex(kind=c16_kind), parameter   :: ic=(0.0_8,1.0_8)
      integer(kind=i4_kind),  allocatable :: sig(:,:)
      integer(kind=i4_kind)               :: ii,jj,kk,n_i,n_j,n_k,ns_i,ns_j
      integer(kind=i4_kind)               :: ns_k,N_terms

      n_cr(-1) = 2; n_cr(0) = 1; n_cr(1) = 2
      cf( 0,1) =  1.0_8;              cf( 0,2) =  0.0_8
      cf( 1,1) =  1.0_8/sqrt(2.0_8);  cf( 1,2) =  cf( 1,1)
      cf(-1,1) = -ic*cf(1,1);         cf(-1,2) = -cf(-1,1)
      allocate(Coeff(8), stat=alloc_stat(32))
      ASSERT(alloc_stat(32).eq.0)
      MEMLOG(size(Coeff))
      allocate(sig(3,8), stat=alloc_stat(32))
      ASSERT(alloc_stat(32).eq.0)
      MEMLOG(size(sig))
      do ii=-1,1
         do jj=-1,1
            do kk=-1,1
               N_terms = 0
               ns_i=-1
               do n_i=1,n_cr(ii)
                  ns_i=-ns_i
                  ns_j=-1
                  do n_j=1,n_cr(jj)
                     ns_j=-ns_j
                     ns_k=-1
                     do n_k=1,n_cr(kk)
                        ns_k=-ns_k
                        N_terms = N_terms + 1
                        Coeff(N_terms) = cf(ii,n_i)*cf(jj,n_j)*cf(kk,n_k)
                        sig(1,N_terms) = ns_i
                        sig(2,N_terms) = ns_j
                        sig(3,N_terms) = ns_k
                     end do
                  end do
               end do
               allocate( triple(ii,jj,kk)%Coeff(N_terms), stat=alloc_stat(32))
               ASSERT(alloc_stat(32).eq.0)
               MEMLOG(size( triple(ii,jj,kk)%Coeff))
               allocate( triple(ii,jj,kk)%sig(3,N_terms), stat=alloc_stat(32))
               ASSERT(alloc_stat(32).eq.0)
               MEMLOG(size(triple(ii,jj,kk)%sig))
               triple(ii,jj,kk)%N_terms = N_terms
               triple(ii,jj,kk)%Coeff   = Coeff(1:N_terms)
               triple(ii,jj,kk)%sig     = sig(:,1:N_terms)
            end do
         end do
      end do

      MEMLOG(-size(Coeff))
      deallocate(Coeff, stat=alloc_stat(32))
      ASSERT(alloc_stat(32).eq.0)
      MEMLOG(-size(sig))
      deallocate(sig, stat=alloc_stat(32))
      ASSERT(alloc_stat(32).eq.0)

    alloc_stat(32)=1

    end subroutine set_triplets

    subroutine shutdown_triplets
      integer(kind=i4_kind)               :: ii,jj,kk
      do ii=-1,1
         do jj=-1,1
            do kk=-1,1
               MEMLOG(-size(triple(ii,jj,kk)%Coeff))
               deallocate( triple(ii,jj,kk)%Coeff, stat=alloc_stat(32))
               ASSERT(alloc_stat(32).eq.0)
               MEMLOG(-size(triple(ii,jj,kk)%sig))
               deallocate( triple(ii,jj,kk)%sig, stat=alloc_stat(32))
               ASSERT(alloc_stat(32).eq.0)
            end do
         end do
      end do
    end subroutine shutdown_triplets

    function par_index(m)
      ! Purpose : calculation of Index "parity" :
      ! par_index(0)=0  par_index(2m)=1  par_index(2m+1)=-1 *
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: m
      !** End of interface *****************************************
      integer(kind=i4_kind) :: par_index

      if(m == 1) then
         par_index = 0
      else if( (m/2)*2 == m) then
         par_index= 1
      else
         par_index = -1
      end if

    end function par_index
    !****************************************************************!

    !****************************************************************!
    function par_indexv(m) result(p)
      ! Purpose : calculation of Index "parity" :
      ! par_index(0)=0  par_index(2m)=1  par_index(2m+1)=-1
      ! TO SUBSTITUTE par_index LATER!
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: m(3)
      !** End of interface *****************************************
      integer(kind=i4_kind) :: p(3)

      where( m == 1 )                 p = 0
      where( m>1 .and. (m/2)*2 == m ) p = 1
      where( m>1 .and. (m/2)*2 /= m ) p =-1

    end function par_indexv


    subroutine Ang_Cmplx_to_RealPP(L_a,L_b,L_c)
      ! Puspose : Complex -> Real conversion of Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      complex(kind=c16_kind), pointer  :: Coeff(:)
      integer(kind=i4_kind)            :: N_terms
      integer(kind=i4_kind)            :: m_r(3),m_c(3),p(3)
      integer(kind=i4_kind)            :: i_term,ma,mb,mc
      integer(kind=i4_kind)            :: kk
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               ! p = par_indexv(m_r) !!!
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                  AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) + Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3))
                enddo
               end do
            end do
         end do
      end do
    end subroutine Ang_Cmplx_to_RealPP

    subroutine Ang_Cmplx_to_RealGa(L_a,L_b,L_c,dervs,AngCMatGb)

      ! Puspose : Complex -> Real conversion of Angular factors
      ! called for lc.eq.0 case

      ! history: now extended with GGa, GGb and GaGb dervs

      implicit none

      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c

      !** End of interface *****************************************

      logical, intent(in) :: dervs
      complex(kind=c16_kind),optional,intent(inout) :: &
                              AngCMatGb(n_jj,-L_a:L_a,-L_b:L_b,-L_c:L_c,3)
      complex(kind=c16_kind), pointer  :: Coeff(:)
      integer(kind=i4_kind)            :: N_terms
      integer(kind=i4_kind)            :: m_r(3),m_c(3),p(3)
      integer(kind=i4_kind)            :: i_term,ma,mb,mc
      integer(kind=i4_kind)            :: kk,k2dr,k_gr

   alldervs: if(dervs) then
     present_gb: if(present(AngCMatGb)) then ! AngCMatGb now present only for coulomb part
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                 AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
               + real(Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3)),kind=r8_kind)

                 AngRMatGa(kk,ma,mb,1,mc) = AngRMatGa(kk,ma,mb,1,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,2,mc) = AngRMatGa(kk,ma,mb,2,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,3,mc) = AngRMatGa(kk,ma,mb,3,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)

                 AngRMatGb(kk,ma,mb,1,mc) = AngRMatGb(kk,ma,mb,1,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,2,mc) = AngRMatGb(kk,ma,mb,2,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,3,mc) = AngRMatGb(kk,ma,mb,3,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)

               k2dr_ab: do k2dr=1,3
                 do k_gr=1,3
                 AngRMatGGa(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGGa(kk,ma,mb,k_gr,k2dr,mc)  &
                  + real(Coeff(i_term)*AngCMatGGa(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)

                 AngRMatGGb(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGGb(kk,ma,mb,k_gr,k2dr,mc)  &
                  + real(Coeff(i_term)*AngCMatGGb(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)

                 AngRMatGaGb(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGaGb(kk,ma,mb,k_gr,k2dr,mc)  &
                  + real(Coeff(i_term)*AngCMatGaGb(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)
                enddo
               enddo k2dr_ab

              enddo
               end do
            end do
         end do
      end do
     endif present_gb
   elseif(integralpar_gradients) then
     if(present(AngCMatGb)) then
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                 AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
                    + real(Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3)),kind=r8_kind)

                 AngRMatGa(kk,ma,mb,1,mc) = AngRMatGa(kk,ma,mb,1,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,2,mc) = AngRMatGa(kk,ma,mb,2,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,3,mc) = AngRMatGa(kk,ma,mb,3,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)

                 AngRMatGb(kk,ma,mb,1,mc) = AngRMatGb(kk,ma,mb,1,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,2,mc) = AngRMatGb(kk,ma,mb,2,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,3,mc) = AngRMatGb(kk,ma,mb,3,mc)  &
                                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)
                enddo
               end do
            end do
         end do
      end do
     else
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                 AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
                    + real(Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3)),kind=r8_kind)
!                 AngRMatGa(kk,ma,mb,mc,1) = AngRMatGa(kk,ma,mb,mc,1)  &
                 AngRMatGa(kk,ma,mb,1,mc) = AngRMatGa(kk,ma,mb,1,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,2,mc) = AngRMatGa(kk,ma,mb,2,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,3,mc) = AngRMatGa(kk,ma,mb,3,mc)  &
                                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)
                enddo
               end do
            end do
         end do
      end do
     endif
    else alldervs
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               ! p = par_indexv(m_r) !!!
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                  AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) + Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3))
                enddo
               end do
            end do
         end do
      end do
    endif alldervs
    end subroutine Ang_Cmplx_to_RealGa


    subroutine Ang_Cmplx_to_RealSA(L_a,L_b,L_c)

      ! Puspose : Complex -> Real conversion of Angular factors
      ! called only for coulomb lc.ne.0 case

      ! history: GGb terms added
      !          GaGb terms added

      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      complex(kind=c16_kind), pointer  :: Coeff(:)
      integer(kind=i4_kind)            :: N_terms
      integer(kind=i4_kind)            :: m_r(3),m_c(3),p(3)
      integer(kind=i4_kind)            :: i_term,ma,mb,mc
      integer(kind=i4_kind)            :: kk,k2dr,k_gr
     if(integralpar_dervs) then
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
             if(mc_required(ma,mb,mc,j)) then
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               terms: do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
             jjs:   do kk=1,n_jj
                 AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
                    + real(Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3)),kind=r8_kind)

                 AngRMatGa(kk,ma,mb,1,mc) = AngRMatGa(kk,ma,mb,1,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,2,mc) = AngRMatGa(kk,ma,mb,2,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,3,mc) = AngRMatGa(kk,ma,mb,3,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)


                 AngRMatGb(kk,ma,mb,1,mc) = AngRMatGb(kk,ma,mb,1,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,2,mc) = AngRMatGb(kk,ma,mb,2,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,3,mc) = AngRMatGb(kk,ma,mb,3,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)

                  ! *** b-atom  2nd dervs
               k2drv: do k2dr=1,3
                 do k_gr=1,3
                  ! *** aa  2nd dervs
                 AngRMatGGa(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGGa(kk,ma,mb,k_gr,k2dr,mc)  &
                    + real(Coeff(i_term)*AngCMatGGa(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)

                  ! *** bb 2nd dervs
                 AngRMatGGb(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGGb(kk,ma,mb,k_gr,k2dr,mc)  &
                 + real(Coeff(i_term)*AngCMatGGb(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)

                  ! *** ab 2nd dervs
                 AngRMatGaGb(kk,ma,mb,k_gr,k2dr,mc) = AngRMatGaGb(kk,ma,mb,k_gr,k2dr,mc)  &
                 + real(Coeff(i_term)*AngCMatGaGb(kk,m_c(1),m_c(2),m_c(3),k_gr,k2dr),kind=r8_kind)

                enddo
               enddo k2drv

                enddo jjs
               enddo terms
              endif
            end do
         end do
      end do
     elseif(integralpar_gradients) then
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
            if(mc_required(ma,mb,mc,j)) then
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                 AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
                    + real(Coeff(i_term)*AngCMat(kk,m_c(1),m_c(2),m_c(3)),kind=r8_kind)

                 AngRMatGa(kk,ma,mb,1,mc) = AngRMatGa(kk,ma,mb,1,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,2,mc) = AngRMatGa(kk,ma,mb,2,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGa(kk,ma,mb,3,mc) = AngRMatGa(kk,ma,mb,3,mc)  &
                    + real(Coeff(i_term)*AngCMatGa(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)

                 AngRMatGb(kk,ma,mb,1,mc) = AngRMatGb(kk,ma,mb,1,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),1),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,2,mc) = AngRMatGb(kk,ma,mb,2,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),2),kind=r8_kind)
                 AngRMatGb(kk,ma,mb,3,mc) = AngRMatGb(kk,ma,mb,3,mc)  &
                    + real(Coeff(i_term)*AngCMatGb(kk,m_c(1),m_c(2),m_c(3),3),kind=r8_kind)
                enddo
               end do
              endif
            end do
         end do
      end do
     else
      do mc=1,2*L_c+1
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
            if(mc_required(ma,mb,mc,j)) then
               m_r      = (/ma,mb,mc/)
               p(1) = par_index(m_r(1))
               p(2) = par_index(m_r(2))
               p(3) = par_index(m_r(3))
               N_terms  =  triple( p(1), p(2), p(3) )%N_terms
               Coeff    => triple(p(1),p(2),p(3)) % Coeff(1:N_terms)
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                do kk=1,n_jj
                  AngRMat(kk,ma,mb,mc) = AngRMat(kk,ma,mb,mc) &
                                       + AngCMat(kk,m_c(1),m_c(2),m_c(3))*Coeff(i_term)
                enddo
               end do
            endif
            end do
         end do
      end do
     endif
    end subroutine Ang_Cmplx_to_RealSA

    subroutine calc_mc_required(L_a,L_b,L_c)
      ! Puspose : Complex -> Real conversion of Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind),  pointer  :: N_terms
      integer(kind=i4_kind)            :: m_r(3),m_c(3),p(3)
      integer(kind=i4_kind)            :: i_term,ma,mb,mc
      n_independent_fcts  = &
                      ua_pointer%symadapt_partner(1,L_c)%n_independent_fcts
      independents: do i_ind=1,n_independent_fcts
       n_contributing_fcts = &
            unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%n_fcts
       magn => unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%m
       eq_atom => unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%I_equal_atom
       contributing: do i_cnt=1,n_contributing_fcts
         mc=magn(i_cnt)
         do ma=1,2*L_a+1
            do mb=1,2*L_b+1
               mc_required(ma,mb,mc,eq_atom(i_cnt))=.true.
               m_r      = (/ma,mb,mc/)
               p(1)     =  par_index(m_r(1)); p(2) = par_index(m_r(2)); p(3) = par_index(m_r(3))
               N_terms  => triple( p(1), p(2), p(3) )%N_terms
               do i_term=1,N_terms
                  m_c(:)   =  triple(p(1),p(2),p(3)) % sig(:,i_term) * m_r(:)/2
                  m_c_required(m_c(1),m_c(2),m_c(3),eq_atom(i_cnt))=.true.
               end do
            end do
         end do
       enddo contributing
      enddo independents
    end subroutine calc_mc_required


    function llam_min(L_a,L_b,L_c)
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind)             :: Lam,j1,j2,j3,i1,i2,m_up,LLam,LLam_min

      mup_llam=-1
      llam_min=L_a+L_b+L_c
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                  m_up=L_a-i1+i2+j2-Lam
                  mup_llam(LLam)=max(mup_llam(LLam),m_up)
        llam_min=min(LLam-m_up,llam_min)
        enddo
        enddo
        enddo
        enddo
        enddo
    end function llam_min

    subroutine Radial_Ang_PotPP(L_a,L_b,i_p)
      ! Purpose : combines Radial and Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b,i_p
      !** End of interface *****************************************

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMat(:,:,:,:),n_jj, &
                        1.0_r8_kind,potential(:,:,:),num)

    end subroutine Radial_Ang_PotPP

    subroutine Radial_Ang_nuc(L_a,L_b)
      ! Purpose : combines Radial and Angular factors
      ! output: ca_dervs_mat
      !         dervs_mat
      !         prim_int_coul_dervs
      use shgi_cntrl, only: is_on, INUSD
      implicit none

      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b
      !** End of interface *****************************************

      integer(kind=i4_kind)              :: ma,mb,k_gr,k2dr,i_grad,ind_c
!     real(kind=r8_kind) :: cart_dervs(9,9)

      grads: if(integralpar_3cob_grad) then

       init_dervsM_equal: if(integralpar_dervs) then

        ! dervs matrices are filled in for each equal atom (c)

        dervsM=0.0_r8_kind
        dervs_totsymM=0.0_r8_kind
        dervs_totsymM_temp=0.0_r8_kind
        ca_dervsM=0.0_r8_kind

       endif init_dervsM_equal

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMat(:,:,:,:),n_jj, &
                        0.0_r8_kind,radang_temp(:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGa(:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_ga(:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMatGb(:,:,:,:,:),n_jj, &
                        0.0_r8_kind,radang_temp_gb(:,:,:,:),num)

          dervs: if(integralpar_dervs) then

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGGa(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_gga(:,:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGGb(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_ggb(:,:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGaGb(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_gagb(:,:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_matGGa(:,:),num,&
                       AngRMat(:,:,:,:),n_jj, &
                       0.0_r8_kind,GGt(:,:,:),num)  ! second gamma derivative

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                       AngRMatGa(:,:,:,:,:),n_jj, &
                       0.0_r8_kind,GtGa(:,:,:,:),num)
                                    !^ ga spherical harmonic gamma Gt derivative

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                        AngRMatGb(:,:,:,:,:),n_jj, &
                        0.0_r8_kind,GtGb(:,:,:,:),num)
                                    !^ gb spherical harmonic gamma Gt derivative

          endif dervs


           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                       AngRMat(:,:,:,:),n_jj, &
                       0.0_r8_kind,Gt(:,:,:),num)

        do k_gr=1,3
        ! gamma argument a and b derivatives

    ! gamma_arg2=gamma_arg_xc**2*(a+b)
    ! gamma_arg_xc = (a*vec_a + b*vec_b)/(a + b)-c
    ! d/da gamma_arg = a / (a + b) : aexp_arr/fact0
    ! d/da gamma_arg2 = - 2*a*gamma_arg_xc

          GtFvec(:,1,k_gr)= two_aexp_arr(:)*gamma_arg_xc(:,k_gr,j)
                                        !   ^                                    ! coordinate dependent factor to be
                                        !                                        ! deffirenciated for 2dervs
                                        !      d/da GtFvec = two_aexp_arr(:)*aexp_arr(:)/fact0(:)

          GtFvec(:,2,k_gr)= two_bexp_arr(:)*gamma_arg_xc(:,k_gr,j)
                                        !   ^                                    ! coordinate dependent factor to be
                                        !                                        ! deffirenciated for 2dervs
                                        !      d/da GtFvec = two_aexp_arr(:)*aexp_arr(:)/fact0(:)
        enddo ! k_gr

       moving: if(moving_a.and.moving_b.and.moving_c) then

               !--------------------------
               !  this is true in standart case
               ! and only this part is eloborated in details
               ! second derivatives are done only fot this part
               !------------------------
!      print*,'Ga,radang_temp_ga,-radangF,-GtFvec,Gt',shape(Ga),shape(radang_temp_ga),shape(radangF),&
!                shape(GtFvec),shape(Gt)

            ma_loop: do ma=1,2*L_a+1
              do mb=1,2*L_b+1

               do k_gr=1,3
                 radangF(:)=radang_temp(:,ma,mb)*two_fact2_vec(:,k_gr)
                         !  coordinate depended via harmonic gamma and overlap

!--------------------------------------------------------------------
                         !  coordinate depended via harmonic gamma and overlap
                 Ga(:,k_gr)= (radang_temp_ga(:,ma,mb,k_gr) -radangF(:) &
                ! ^                        ^               ^
                                      -GtFvec(:,1,k_gr)*Gt(:,ma,mb) ) ! ga gamma contribs -> 2 contribs to 2dervs
                                       !        ^                     ! radang_temp_ga ->3-contribs:
                                                                      ! radang_temp_gga,overlap,gamma


                 Gb(:,k_gr)= ( radang_temp_gb(:,ma,mb,k_gr)+radangF(:) &
                ! ^                         ^              ^          ! (b) overl = -(a) overl
                                      -GtFvec(:,2,k_gr)*Gt(:,ma,mb) ) ! gb gamma contribs
                                       !        ^
!----------------------------------------------------------------------

!                 gradM(:,mb,ma,k_gr)  =gradM(:,mb,ma,k_gr)  +Ga(:,k_gr)
!                 gradM(:,mb,ma,k_gr+3)=gradM(:,mb,ma,k_gr+3)+Gb(:,k_gr)
               enddo ! k_gr


   sd2: if(integralpar_dervs) then
          !   we have one overlap contrib for all terms but
          !   separate harmonic and gamma contribs
          !            + prefacror contribs
             do k2dr=1,3
             do k_gr=1,3
              GGa(:,k_gr,k2dr)= ( radang_temp_gga(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,1,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                              -   GtFvec(:,1,k_gr)*GtGa(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr) & ! 3rd term Gt gamma contib + second term below
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr) & ! gamma radangF
                         - radang_temp_ga(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                                   - Ga(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs

              GGb(:,k_gr,k2dr)= ( radang_temp_ggb(:,ma,mb,k_gr,k2dr)  &         ! 1st term (radang_temp_ga) harm contrib
                                              ! ^                                 (b) derivative
                                -GtGb(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  &         ! 1st term (radang_temp_ga) gamma contrib
                                 !  ^                     !  ^                    (b) derivative
                              -   GtFvec(:,2,k_gr)*GtGb(:,ma,mb,k2dr) &         ! 3rd term Gt hamonic contrib
                                         ! ^          ^                           (b) derivative
                     + GGt(:,ma,mb)*GtFvec(:,2,k_gr)*GtFvec(:,2,k2dr) &         ! 3rd term Gt gamma contib + second term below
                                          !  ^                ^                   (b) derivative
                         - Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & ! gamma radangF
                     !   ^                                            ^         ! 2nd term with different signs in Ga and Gb
                         + radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! 2nd sol-harm radangF
                     !   ^              ^                                       ! sign of (b) deriv as in Gb
                                   + Gb(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs
                                 ! ^  ^                                 ! sign changed because (b) overl derv
                                                                        ! should be equal to -(a) derv

              !----------------------------------------------------------
              GaGb(:,k_gr,k2dr)= ( radang_temp_gagb(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                 !              ^ ^
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                                 !  ^                        ^
                              -   GtFvec(:,1,k_gr)*GtGb(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                                 !         ^a         ^
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,2,k2dr) & ! 3rd term Gt gamma contib + second term below
                                 !           ^a               ^b
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & !2nd term gamma radangF
                     !   ^a-conrib sign                               ^b
                         - radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                     !   ^a-contr sign  ^
                                   + Ga(:,k_gr)*two_fact2_vec(:,k2dr) &! all GA overlap contribs
                                !  ^b ^    sign of b-cotrib
                                  )
                     !-----------------------------------------
             enddo
             enddo

           do k2dr=1,3

            !    2nd term  (radangF) prefactor contribs
            !    ( derivatives by two_fact2_vec fact2*2* (a-b) )
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib,
                                                                                 ! sign of radangF contrib
            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !                                ^  (sign of radangF contrib)*(sign (-) of b in (a-b) argument)

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)+radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !^ ^               ^ ^             ^ ((-) sing of a-contrib)x((-) sign of (a-b) argument)

            !    comment to check Gt contribution whan harmonic and overlap fixed

            !   g4 - (d/da GtFvec )xGt
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr) &
                            -Gt(:,ma,mb)*two_aexp_arr(:)*aexp_arr(:)/fact0(:)
            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr) &
                            -Gt(:,ma,mb)*two_bexp_arr(:)*bexp_arr(:)/fact0(:)

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr) &
                             -Gt(:,ma,mb)*two_aexp_arr(:)*bexp_arr(:)/fact0(:)
           enddo


!         if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i) then

!           if(k.eq.1.and.i.eq.1.and.abs(sum(GGa(:,1,1))).gt.1.0e-10) &
!             print*,'i,lc,ma,mb, Ga GGa',i,lc,mc,ma,mb, sum(Ga(:,1)),sum(GGa(:,1,1)),' nuc'

!!            if(k.eq.1.and.i.eq.1.and.abs(sum(GGb(:,1,1))).gt.1.0e-10) &
!!              print*,'i,lc,ma,mb, Gb GGb',i,lc,mc,ma,mb, sum(Gb(:,1)),sum(GGb(:,1,1)),' nuc'

!!            if(k.eq.1.and.i.eq.1.and.abs(sum(GaGb(:,1,1))).gt.1.0e-10) &
!!              print*,'i,lc,ma,mb, Ga GaGb',i,lc,mc,ma,mb, sum(Ga(:,1)),sum(GaGb(:,1,1)),' nuc'
!          endif


             ! probable summing up here for dervsM is not actually required

              dervsM(:,mb,ma,1:3,1:3)=dervsM(:,mb,ma,1:3,1:3)+GGa(:,:,:)
              dervsM(:,mb,ma,4:6,4:6)=dervsM(:,mb,ma,4:6,4:6)+GGb(:,:,:)
              dervsM(:,mb,ma,1:3,4:6)=dervsM(:,mb,ma,1:3,4:6)+GaGb(:,:,:)

              do k_gr=1,3    ! for (4:6,1:3) one has transposed  GaGb
              do k2dr=1,3
               dervsM(:,mb,ma,k2dr+3,k_gr)=dervsM(:,mb,ma,k2dr+3,k_gr)+GaGb(:,k_gr,k2dr)
              enddo
              enddo
!           print*,'ab GGa GGb GaGb', sum(GGa(:,3,3)),sum(GGb(:,3,3)),sum(GaGb(:,3,3))


          !----------------------------------------
          ! FILL IN ca_dervsM FOR EACH EQUAL (c) AND THEN REMAP
          ! AND SUM IT UP IN  ca_dervs_mat FOR EACH UNIQUE (c)

      call fill_ca_dervsM(ma,mb)


          ! DONE FILL IN ca_dervsM

          !-------------------------------------------
          ! *** symetrize (CC) dervs by i_grad index  and
          !     store result in dervs_totsymM_temp for each equal (c)

      rot_dervs1: if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim
        do k2dr=1,3                 !(4)
         dervs_totsymM_temp(:,mb,ma,i_grad,k2dr)=dervs_totsymM_temp(:,mb,ma,i_grad,k2dr) &
              +rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GGb(:,1,k2dr)+GaGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
              +rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GGb(:,2,k2dr)+GaGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
              +rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GGb(:,3,k2dr)+GaGb(:,3,k2dr)+GaGb(:,k2dr,3))
        enddo
       enddo
!      print*,'dervs cc',j,sum(GGa(:,3,3)+GGb(:,3,3)+GaGb(:,3,3)+GaGb(:,3,3))

      else if(grad_dim.gt.0) then
        do k_gr=1,3   !(5)
        do k2dr=1,3
         dervs_totsymM_temp(:,mb,ma,k_gr,k2dr)=dervs_totsymM_temp(:,mb,ma,k_gr,k2dr) &
        +(GGa(:,k_gr,k2dr)+GGb(:,k_gr,k2dr)+GaGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
        enddo
        enddo
      endif rot_dervs1

          ! symetrize (cc) dervs in dervs_totsymM_temp by k_gr index and
          ! remap it to grad_totsymM again for each equal (c)

       rot_dervs2: if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim
         dervs_totsymM(:,mb,ma,:,i_grad)=dervs_totsymM(:,mb,ma,:,i_grad) &
            +rotmat_eq(i_grad,1,j)*dervs_totsymM_temp(:,mb,ma,:,1) &
            +rotmat_eq(i_grad,2,j)*dervs_totsymM_temp(:,mb,ma,:,2) &
            +rotmat_eq(i_grad,3,j)*dervs_totsymM_temp(:,mb,ma,:,3)
       enddo

       else if(grad_dim.gt.0) then
          dervs_totsymM(:,mb,ma,:,:)=dervs_totsymM(:,mb,ma,:,:) &
              +dervs_totsymM_temp(:,mb,ma,:,:)
       endif rot_dervs2

      endif sd2

    enddo
   enddo ma_loop


   ! NUCL: NUCLEAR ATTRACTION DERIVATIVES HERE:
   remap_dervs: if(integralpar_dervs .and. is_on(INUSD) ) then ! if false then no nuclear attraction dervs

         if( .true. .or. &
            (quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i)) then !!!!

         ! put .false. to see only real 3center contribs , to be just commented after all checks

         ma_lp: do ma=1,2*L_a+1
               do mb=1,2*L_b+1

         abmap_nuc_dervs: if(.true.) then !(3) - 3rd  of 3 dervs contribs
           do k_gr=1,6
            do k2dr=1,6
            dervs_mat(:,:,mb,ma,k_gr,k2dr)=dervs_mat(:,:,mb,ma,k_gr,k2dr) &
                -unpack(dervsM(:,mb,ma,k_gr,k2dr),cutoff,zero)
            enddo
           enddo
!           print*,'ab GGa GGb GaGb', sum(GGa(:,3,3)),sum(GGb(:,3,3)),sum(GaGb(:,3,3))

         endif abmap_nuc_dervs

!          now dervs_mat is left to be remaped to prim_coul_dervs in due turn



        fst_grad_dim_ind: do i_grad=1,grad_dim ! only if moving_c

           camap_nuc_dervs: if(.true.) then   ! second of 3 dervs contins
             do k2dr=1,6
              ca_dervs_mat(:,:,mb,ma,i_grad,k2dr)=ca_dervs_mat(:,:,mb,ma,i_grad,k2dr) &
                  -unpack(ca_dervsM(:,mb,ma,i_grad,k2dr),cutoff,zero)
             enddo
!          print*,'ca_dervsM a',sum(ca_dervsM(:,:,:,1,3)), &
!             sum(GGa(:,3,3)),sum(GaGb(:,3,3)),rotmat_eq(1,3,j),j, &
!             sum(ca_dervs_mat(:,:,:,:,1,3))
                   !  now after leaving program ca_dervs_mat filled in for each unique (c)
                   !  will be remaped to prim_coul_dervs
                   !  in calc_3c_fitcont route together with coul and other contribs
          endif camap_nuc_dervs

         ccmap_nuc_dervs: if(.true.) then  ! first of 3 dervs contribs
             ind_c=gradient_index(imc)-1
             do k2dr=1,grad_dim
               prim_int_coul_dervs(ind_c+i_grad,ind_c+k2dr)%m(:,:,mb,ma)= &    ! (1) - first of 3 contribs
                 prim_int_coul_dervs(ind_c+i_grad,ind_c+k2dr)%m(:,:,mb,ma) &
                   -unpack(dervs_totsymM(:,mb,ma,i_grad,k2dr),cutoff,zero)
             enddo
         endif ccmap_nuc_dervs

            enddo fst_grad_dim_ind

          enddo
       enddo ma_lp

!     print*,'cartesian standard order dervs '
!     ! if uncommented can be used for check in calculation
!     ! of SS rype integrals
!
!     cart_dervs(1:3,1:3)=GGa(1,:,:)
!     cart_dervs(4:6,4:6)=GGb(1,:,:)
!     cart_dervs(1:3,4:6)=GaGb(1,:,:)
!     cart_dervs(4:6,1:3)=transpose(GaGb(1,:,:))
!     cart_dervs(1:3,7:9)=-GGa(1,:,:)-GaGb(1,:,:)
!     cart_dervs(4:6,7:9)=-GGb(1,:,:)-transpose(GaGb(1,:,:))
!     cart_dervs(7:9,1:3)=-GGa(1,:,:)-transpose(GaGb(1,:,:))
!     cart_dervs(7:9,4:6)=-GGb(1,:,:)-GaGb(1,:,:)
!     cart_dervs(7:9,7:9)=GGa(1,:,:)+GGb(1,:,:)+GaGb(1,:,:)+transpose(GaGb(1,:,:))
!
!     do k_gr=1,9
!      write(*,'(10f8.4)') cart_dervs(k_gr,:),sum(cart_dervs(k_gr,:))
!     enddo
!      write(*,'(10f8.4)') sum(cart_dervs(:,:),1)

    endif ! real 3 center only   !!!!!!!!!! to be commented

   endif remap_dervs

  endif  moving

 endif grads

 end subroutine Radial_Ang_nuc

      subroutine fill_ca_dervsM(ma,mb)
      integer(kind=i4_kind), intent(in):: ma,mb
      integer(kind=i4_kind)              :: k_gr,k2dr,i_grad
      cx_rot_dervs: if(do_rotation_eq(j)) then
      if(moving_a) then
       do i_grad=1,grad_dim
!        ca_dervsM(:,mb,ma,i_grad,1:3)=ca_dervsM(:,mb,ma,i_grad,1:3) &
!              -rotmat_eq(i_grad,1,j)*(GGa(:,1,:)+GaGb(:,1,:)) &
!              -rotmat_eq(i_grad,2,j)*(GGa(:,2,:)+GaGb(:,2,:)) &
!              -rotmat_eq(i_grad,3,j)*(GGa(:,3,:)+GaGb(:,3,:))
         do k2dr=1,3   ! transposed GaGb, intended to (7:9,1:3) block
          ca_dervsM(:,mb,ma,i_grad,k2dr)=ca_dervsM(:,mb,ma,i_grad,k2dr) &
              -rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GaGb(:,k2dr,1)) &
              -rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GaGb(:,k2dr,2)) &
              -rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GaGb(:,k2dr,3))
         enddo

       enddo
      endif

      if(moving_b) then              !(2)
       do i_grad=1,grad_dim
!         do k2dr=1,3
!          ca_dervsM(:,mb,ma,i_grad,3+k2dr)=ca_dervsM(:,mb,ma,i_grad,3+k2dr) &
!              -rotmat_eq(i_grad,1,j)*(GGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
!              -rotmat_eq(i_grad,2,j)*(GGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
!              -rotmat_eq(i_grad,3,j)*(GGb(:,3,k2dr)+GaGb(:,k2dr,3))
!         enddo  ! not transposed GaGb, intended for (7:9,4:6) block
        ca_dervsM(:,mb,ma,i_grad,4:6)=ca_dervsM(:,mb,ma,i_grad,4:6) &
              -rotmat_eq(i_grad,1,j)*(GGb(:,1,:)+GaGb(:,1,:)) &
              -rotmat_eq(i_grad,2,j)*(GGb(:,2,:)+GaGb(:,2,:)) &
              -rotmat_eq(i_grad,3,j)*(GGb(:,3,:)+GaGb(:,3,:))

       enddo
      endif

      else if(grad_dim.gt.0) then
        if(moving_a) &
         ca_dervsM(:,mb,ma,:,1:3)=ca_dervsM(:,mb,ma,:,1:3)-(GGa(:,:,:)+GaGb(:,:,:))

        if(moving_b) then           !(3)
        do k_gr=1,3
         do k2dr=1,3
          ca_dervsM(:,mb,ma,k_gr,3+k2dr)=ca_dervsM(:,mb,ma,k_gr,3+k2dr) &
                                 -(GGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
         enddo
        enddo
        endif

      endif cx_rot_dervs
      end subroutine fill_ca_dervsM

    subroutine Radial_Ang_PotGa(L_a,L_b,i_p)
      ! Purpose : combines Radial and Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b,i_p
      !** End of interface *****************************************
      integer(kind=i4_kind)              :: ma,mb,k_gr
      real(kind=r8_kind),allocatable :: potentialGt(:,:,:)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMatGa(:,:,:,:,:),n_jj, &
                        1.0_r8_kind,potentialGa(:,:,:,:),num)

           allocate(potentialGt(num,(2*L_a+1),(2*L_b+1)),stat=alloc_stat(79))
           ASSERT(alloc_stat(79).eq.0)
           MEMLOG(size(potentialGt))
           potentialGt=0.0_r8_kind

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                        AngRMat(:,:,:,:),n_jj, &
                        1.0_r8_kind,potentialGt(:,:,:),num)

      do ma=1,2*L_a+1
         do mb=1,2*L_b+1
          do k_gr=1,3
          potentialGa(:,ma,mb,k_gr)=potentialGa(:,ma,mb,k_gr)  &
                -two*aexp_arr(:)*potentialGt(:,ma,mb)*(gamma_arg(:,k_gr)-xc(k_gr))

          enddo
         end do
      end do
         MEMLOG(-size(potentialGt))
         deallocate(potentialGt,stat=alloc_stat(79))

    end subroutine Radial_Ang_PotGa

    subroutine Radial_Ang_EWGa(L_a,L_b)
      ! Purpose : combines Radial and Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b
      !** End of interface *****************************************
      integer(kind=i4_kind)              :: ma,mb,k_gr,k2dr
      real(kind=r8_kind),allocatable :: Gt(:,:,:)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMatGa(:,:,:,:,:),n_jj, &
                        0.0_r8_kind,radang_temp_ga,num)
                      ! nucEWGa is a sum over all PCs while
                      ! radang_temp_ga is a point contrib

           if(integralpar_dervs.or.ewpcdervs) then
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGa(:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_gb,num)
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGGa(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_gga(:,:,:,:,:),num)
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGGb(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_ggb(:,:,:,:,:),num)
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMatGaGb(:,:,:,:,:,:),n_jj, &
                       0.0_r8_kind,radang_temp_gagb(:,:,:,:,:),num)
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_matGGa(:,:),num,&
                       AngRMat(:,:,:,:),n_jj, &
                       0.0_r8_kind,GGt(:,:,:),num)  ! second gamma derivative
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                       AngRMatGa(:,:,:,:,:),n_jj, &
                       0.0_r8_kind,GtGa(:,:,:,:),num)
                                    !^ ga spherical harmonic gamma Gt derivative

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                        AngRMatGb(:,:,:,:,:),n_jj, &
                        0.0_r8_kind,GtGb(:,:,:,:),num)
                                    !^ gb spherical harmonic gamma Gt derivative


            GtFvec(:,2,:)=spread(two_bexp_arr(:),2,3)*gamma_argXC

           endif
            GtFvec(:,1,:)=spread(two_aexp_arr(:),2,3)*gamma_argXC

           allocate(Gt(num,(2*L_a+1),(2*L_b+1)))
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       1.0_r8_kind, radialNuc_matGa(:,:),num,&
                       !                         ^
                       AngRMat(:,:,:,:),n_jj, &
                       0.0_r8_kind,Gt(:,:,:),num)


      do ma=1,2*L_a+1
         do mb=1,2*L_b+1
          Ga(:,:)=radang_temp_ga(:,ma,mb,:)-GtFvec(:,1,:)*spread(Gt(:,ma,mb),2,3)
          if(ewpcgrads) nucEWGa(:,ma,mb,:)=nucEWGa(:,ma,mb,:) + Ga

      if(ewpcdervs) then
       Ga=Ga-two_fact2_vec*spread(radang_temp(:,ma,mb),2,3)
       Gb=radang_temp_gb(:,ma,mb,:)-GtFvec(:,2,:)*spread(Gt(:,ma,mb),2,3) &
         +two_fact2_vec*spread(radang_temp(:,ma,mb),2,3)
       GGa=radang_temp_gga(:,ma,mb,:,:)
       GGb=radang_temp_ggb(:,ma,mb,:,:)
       GaGb=radang_temp_gagb(:,ma,mb,:,:)
       do k2dr=1,3
       do k_gr=1,3
        GGa(:,k_gr,k2dr)=GGa(:,k_gr,k2dr)-GtGa(:,ma,mb,k_gr)*GtFvec(:,1,k2dr) &
                                         -GtGa(:,ma,mb,k2dr)*GtFvec(:,1,k_gr) &
                        +GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr)       &
                        + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr)  &
                        -radang_temp_ga(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr)   &
                        -Ga(:,k_gr)*two_fact2_vec(:,k2dr)
        GGb(:,k_gr,k2dr)=GGb(:,k_gr,k2dr)-GtGb(:,ma,mb,k_gr)*GtFvec(:,2,k2dr) &
                                         -GtGb(:,ma,mb,k2dr)*GtFvec(:,2,k_gr) &
                        +GGt(:,ma,mb)*GtFvec(:,2,k_gr)*GtFvec(:,2,k2dr)       &
                        - Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr)  &
                        +radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr)   &
                        +Gb(:,k_gr)*two_fact2_vec(:,k2dr)
        GaGb(:,k_gr,k2dr)=GaGb(:,k_gr,k2dr)-GtGa(:,ma,mb,k_gr)*GtFvec(:,2,k2dr) &
                                           -GtGb(:,ma,mb,k2dr)*GtFvec(:,1,k_gr) &
                        +GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,2,k2dr)         &
                        + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr)    &
                        -radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr)     &
                        + Ga(:,k_gr)*two_fact2_vec(:,k2dr)
       enddo
        GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)   &
                        -Gt(:,ma,mb)*two_aexp_arr(:)*aexp_arr(:)/fact0(:)
        GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)   &
                        -Gt(:,ma,mb)*two_bexp_arr(:)*bexp_arr(:)/fact0(:)
        GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)+radang_temp(:,ma,mb)*two*fact2(:) &
                         -Gt(:,ma,mb)*two_aexp_arr(:)*bexp_arr(:)/fact0(:)
       enddo
             if(ewpcdervs) then
              dervsM(:,mb,ma,1:3,1:3)=dervsM(:,mb,ma,1:3,1:3)+GGa(:,:,:)
              dervsM(:,mb,ma,4:6,4:6)=dervsM(:,mb,ma,4:6,4:6)+GGb(:,:,:)
              dervsM(:,mb,ma,1:3,4:6)=dervsM(:,mb,ma,1:3,4:6)+GaGb(:,:,:)
             endif
      endif
         end do
      enddo

         deallocate(Gt)
    end subroutine Radial_Ang_EWGa

    subroutine Radial_Ang_PotGG(L_a,L_b,L_c)
      ! Purpose : combines Radial and Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind) :: ma,mb,ism,i_ma,index,grad_dim,nn,i_grad,N_pc,km
      real(kind=r8_kind), allocatable:: potentialGt(:,:,:,:)
      integer(kind=i4_kind) :: alloc_stat

      N_pc=0
      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

        allocate(potentialGt(num,2*L_a+1,2*L_b+1,3), &
             help_arr_gr1(num,3,2*L_b+1,2*L_a+1),stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
        potentialGt=0
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMat(:,:,:,:),n_jj, &
                        1.0_r8_kind,potentialGt(:,:,:,:),num)
         potentialGc=potentialGc+potentialGt
!!$         if(i.eq.1.and.j.eq.1) print*,sum(abs(potentialGt(:,:,:,:))),' Gc'

           ism=to_calc_grads%i_symm_sort(i,j)
           do i_ma=1,N_moving_unique_atoms+N_pc
              if(i_ma <= N_moving_unique_atoms) then
                 index = gradient_index(i_ma) - 1
                 grad_dim=gradient_index(i_ma+1)-gradient_index(i_ma)
              else
                 km=i_ma-N_moving_unique_atoms
                 index = surf_points_grad_index(km) - 1
                 grad_dim=surf_points_grad_index(km+1)-surf_points_grad_index(km)
              end if

            help_arr_gr1=0.0_r8_kind
            do nn=1,grad_dim
             do ma=1,2*l_a+1
              do mb=1,2*l_b+1
               help_arr_gr1(:,nn,mb,ma)=help_arr_gr1(:,nn,mb,ma) - &
               potentialGt(:,ma,mb,2)*to_calc_grads%dxyz_totsyms(nn,i_ma)%m(1,ism)
               help_arr_gr1(:,nn,mb,ma)=help_arr_gr1(:,nn,mb,ma) - &
               potentialGt(:,ma,mb,3)*to_calc_grads%dxyz_totsyms(nn,i_ma)%m(2,ism)
               help_arr_gr1(:,nn,mb,ma)=help_arr_gr1(:,nn,mb,ma) - &
               potentialGt(:,ma,mb,1)*to_calc_grads%dxyz_totsyms(nn,i_ma)%m(3,ism)
              enddo
             enddo
            enddo

        if(new_solv_grad) then
           if(i_ma <= N_moving_unique_atoms) then
              do i_grad=1,grad_dim
                 grad_mat_p=>prim_int_3cob_solv_grad(index+i_grad)%m
                 do ma=1,2*l_a+1
                    do mb=1,2*l_b+1
                       grad_mat_p(:,:,mb,ma)=grad_mat_p(:,:,mb,ma)&
                            +unpack(help_arr_gr1(:,i_grad,mb,ma),cutoff,zero) !!!!
                    end do
                 end do
!!$             if(i.eq.1.and.j.eq.1) print*, sum(abs(grad_mat_p(1,1,:,:))), 'gr mat'
!!$             if(i.eq.2) print*,sum(abs(help_arr_gr1)),'3c sum abs help_arr_gr1'
              end do
           else
              do i_grad=1,grad_dim
                 grad_mat_p=>prim_int_3cob_solv_grad_pc(index+i_grad)%m
                 do ma=1,2*l_a+1
                    do mb=1,2*l_b+1
                       grad_mat_p(:,:,mb,ma)=grad_mat_p(:,:,mb,ma)&
                            -unpack(help_arr_gr1(:,i_grad,mb,ma),cutoff,zero) !!!!
                    end do
                 end do
              end do
           end if
        endif

           enddo

        deallocate(potentialGt,help_arr_gr1,stat=alloc_stat)
        if(alloc_stat.ne.0) call error_handler("deallocate potentialGt failed")

    end subroutine Radial_Ang_PotGG

    subroutine Radial_Ang_EWGG(L_a,L_b,L_c)
      ! Purpose : combines Radial and Angular factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      real(kind=r8_kind), allocatable:: Gt(:,:,:,:)
      integer(kind=i4_kind) :: alloc_stat

        allocate(Gt(num,2*L_a+1,2*L_b+1,3), stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                        AngRMat(:,:,:,:),n_jj, &
                        0.0_r8_kind,Gt(:,:,:,:),num)
         nucEWGc=nucEWGc+Gt

        deallocate(Gt,stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)

    end subroutine Radial_Ang_EWGG

    subroutine Radial_Ang_3cSA(L_a,L_b,L_c)

      ! Purpose : combines Radial and Angular factors

      ! history: 2nd GGb derivatives
      !          now extended with GaGb derivatives
      !          A x B derivatives complited

      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind),  intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind)              :: ma,mb,mc,k_gr,k2dr,i_grad
      real(kind=r8_kind), pointer      :: coeff(:)

         n_independent_fcts  = &
                      ua_pointer%symadapt_partner(1,L_c)%n_independent_fcts
         independents: do i_ind=1,n_independent_fcts

         n_contributing_fcts = &
            unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%n_fcts

         magn => unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%m
         eq_atom => unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%I_equal_atom
         coeff => unique_atoms(i)%symadapt_partner(1,L_c)%symadapt(i_ind,1)%c

 if(integralpar_gradients) then
   grad_totsymM=0.0_r8_kind
   gradM=0.0_r8_kind
   if(integralpar_dervs) then
   dervsM=0.0_r8_kind
   dervs_totsymM=0.0_r8_kind
   dervs_totsymM_temp=0.0_r8_kind
   ca_dervsM=0.0_r8_kind
   endif
 endif

         contributing: do i_cnt=1,n_contributing_fcts
                                 ! depend on coordinates via overlap, gamma and
                                 ! spherical harmonics parts
           mc=magn(i_cnt)
           j=eq_atom(i_cnt)

          !                        energy contribution
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMat(:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,radang_temp(:,:,:),num)

 gr_els:if(integralpar_gradients) then
   gr: if(radang_lco_grad.and.new_3c_co) then

        !                    1st  ga harmonic derivatives
           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMatGa(:,:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,radang_temp_ga(:,:,:,:),num)


           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMatGb(:,:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,radang_temp_gb(:,:,:,:),num)
        !                           gb harmonic derivatives

          if(integralpar_dervs) then

             ! 2nd deriv radang_temp_ga harmonic contrib (1)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMatGGa(:,:,:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,radang_temp_gga(:,:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMatGGb(:,:,:,:,:,mc,j),n_jj, &
                                !  ^
                        0.0_r8_kind,radang_temp_ggb(:,:,:,:,:),num)
                                   !              ^

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*9,n_jj, &
                       coeff(i_cnt), radial3cmat(:,:,j),num,&
                        AngR3cMatGaGb(:,:,:,:,:,mc,j),n_jj, &
                                !  ^
                        0.0_r8_kind,radang_temp_gagb(:,:,:,:,:),num)

             ! to check radang_temp_gga  gamma and overlap contribs
             ! have to be fixed


!--------------------------------------------------
!         if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i) then
!          if(k.eq.1.and.j.eq.1) print*, 'Ang lc mc',l_c,mc, sum(AngR3cMat(:,:,:,mc,j)), &
!                                                        sum(AngR3cMatGa(:,:,:,1,mc,j)), &
!                                                     sum(AngR3cMatGGa(:,:,:,1,1,mc,j))

!          if(k.eq.1.and.j.eq.1) print*, 'radang_temp ga gga i l_c', &
!                          i,l_c, sum(radang_temp),sum(radang_temp_ga(:,:,:,1)), &
!                                                 sum(radang_temp_gga(:,:,:,1,1))
!         endif
!-----------------------------------------------------


         ! now this SH contrib checked

!!----------------------------------------------------------------
!!
!!         this fragment commented after initial check of GGa contrib
!!
!          ! radang_temp_ga gamma contrib
!!           radang_temp_gadg=0.0_r8_kind
!!           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
!!                       coeff(i_cnt), radial3cmatG(:,:,j),num,&
!!                        AngR3cMatGa(:,:,:,:,mc,j),n_jj, &
!!                        1.0_r8_kind,radang_temp_gadg(:,:,:,:),num)
!!---------------------------------------------------

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       coeff(i_cnt), radial3cmatGG(:,:,j),num,&
                        AngR3cMat(:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,GGt(:,:,:),num)

           ! next find what is derivative by argument
           !                        gamma derivatives

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       coeff(i_cnt), radial3cmatG(:,:,j),num,&
                        AngR3cMatGa(:,:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,GtGa(:,:,:,:),num)

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1)*3,n_jj, &
                       coeff(i_cnt), radial3cmatG(:,:,j),num,&
                        AngR3cMatGb(:,:,:,:,mc,j),n_jj, &
                              !   ^
                        0.0_r8_kind,GtGb(:,:,:,:),num)
                                   !   ^

          endif ! integralpar_2dervs

           call dgemm('n','n',num,(2*L_a+1)*(2*L_b+1),n_jj, &
                       coeff(i_cnt), radial3cmatG(:,:,j),num,&
                        AngR3cMat(:,:,:,mc,j),n_jj, &
                        0.0_r8_kind,Gt(:,:,:),num)

         !----------------------------------------------
         !              first  gamma derivatives
         !    to 2dervs will contrubute
         !    g1(d/da AngR3cMat)+g2(d/da gamma Gt)+g3(d/da overlap Gt)
         !    for numerical  check fix SH and overlap variation
         !-----------------------------------------------

        do k_gr=1,3
        ! gamma argument a and b derivatives
        ! also to be applied for GGt case

          GtFvec(:,1,k_gr)= two_aexp_arr(:)*gamma_arg_xc(:,k_gr,j)*exp_arg(:,1)
                                        !   ^                                    ! coordinate dependent factor to be
                                        !                                        ! deffirenciated for 2dervs
                                        !      d/da GtFvec = two_aexp_arr(:)*exp_arg(:,1)*aexp_arr(:)/fact0(:)

          GtFvec(:,2,k_gr)= two_bexp_arr(:)*gamma_arg_xc(:,k_gr,j)*exp_arg(:,1)
                                        !     d/db  GtFvec = two_bexp_arr(:)*exp_arg(:,1)*bexp_arr(:)/fact0(:)
        enddo

!--------------------------------------------
        !here we can check gamma derivatives of Gt while harmonic and overlap
        ! variations fixed

!         if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i) then
!          if(k.eq.1.and.j.eq.1) then
!          do ma=1,2*L_a+1
!           do mb=1,2*L_b+1
!           if(abs(sum(GGt(:,ma,mb)*GtFvec(:,1,1))).gt.1.0e-10) then
!            print*,'Gt gamma var',ma,mb,l_c,mc,i, sum(Gt(:,ma,mb)),sum(GGt(:,ma,mb)*GtFvec(:,1,1))
!           endif
!
!            if(abs(sum(radang_temp_gga(:,ma,mb,1,1)-GtFvec(:,1,1)*radang_temp_gadg(:,ma,mb,1))).gt.1.0e-10) then
!             print*, 'radangGa', sum(radang_temp_ga(:,ma,mb,1)), &
!                                sum(radang_temp_gga(:,ma,mb,1,1)-GtFvec(:,1,1)*radang_temp_gadg(:,ma,mb,1))
!            endif
!          enddo
!          enddo
!
!           if(abs(sum(Ga(:,1))).gt.1.0e-10) then
!           print*,'sum abs AngR3cMatGa dadg gga', sum(abs(AngR3cMatGa)), &
!                                              sum(abs(radang_temp_gadg(:,ma,mb,1))), &
!                                              sum(abs(radang_temp_gga))
!           print*, 'radangGa ', sum(radang_temp_ga(:,:,:,1)),sum(Ga(:,1))
!           endif
!
!          endif
!         endif
!---------------------------------------------
         !now this contribution checked

       moving: if(moving_a.and.moving_b.and.moving_c) then

               !--------------------------
               !  this is true in standart case
               ! and only this part is eloborated in details
               ! second derivatives are done only fot this part
               !------------------------
            do ma=1,2*L_a+1
              do mb=1,2*L_b+1

               do k_gr=1,3
                 radangF(:)=radang_temp(:,ma,mb)*two_fact2_vec(:,k_gr)
                         !  coordinate depended via harmonic gamma and overlap

!--------------------------------------------------------------------
                         !  coordinate depended via harmonic gamma and overlap
                 Ga(:,k_gr)= (radang_temp_ga(:,ma,mb,k_gr) -radangF(:) &
                ! ^                        ^               ^
                                      -GtFvec(:,1,k_gr)*Gt(:,ma,mb) ) ! ga gamma contribs -> 2 contribs to 2dervs
                  ! radang_temp_ga ->3-contribs: radang_temp_gga,overlap,gamma
!----------------------------------------------------------------------


                 Gb(:,k_gr)= ( radang_temp_gb(:,ma,mb,k_gr)+radangF(:) &
                ! ^                         ^              ^          ! (b) overl = -(a) overl
                                      -GtFvec(:,2,k_gr)*Gt(:,ma,mb) ) ! gb gamma contribs
                 gradM(:,mb,ma,k_gr)  =gradM(:,mb,ma,k_gr)  +Ga(:,k_gr)
                 gradM(:,mb,ma,k_gr+3)=gradM(:,mb,ma,k_gr+3)+Gb(:,k_gr)
               enddo

!---------------------------------------------------------------
!          if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i) then
!            if(k.eq.1.and.i.eq.1.and.abs(sum(Ga(:,1))).gt.1.0e-10) then
!            print*,'radangF ', sum(radang_temp(:,ma,mb)*two_fact2_vec(:,1)), &
!                               sum(Ga(:,1)*two_fact2_vec(:,1)+radang_temp(:,ma,mb)*two*fact2(:) )
!            print*,'Ga overlap', sum(Ga(:,1)),-sum(Ga(:,1)*two_fact2_vec(:,1))
!            print*,'ovl radang_temp_ga',sum(radang_temp_ga(:,ma,mb,1)), &
!                                       -sum(radang_temp_ga(:,ma,mb,1)*two_fact2_vec(:,1))
!            print*,'ovl Gt',sum(Gt(:,ma,mb)),-sum(Gt(:,ma,mb)*two_fact2_vec(:,1))
!
!            print*,'two_fact2_vec',sum(two_fact2_vec(:,1)),sum(two*fact2(:))
!            endif
!          endif
!-----------------------------------------------------

         sd2: if(integralpar_dervs) then
          !   we have one overlap contrib for all terms but
          !   separate harmonic and gamma contribs
          !            + prefacror contribs
             do k2dr=1,3
             do k_gr=1,3
              GGa(:,k_gr,k2dr)= ( radang_temp_gga(:,ma,mb,k_gr,k2dr)  &         ! 1st term (radang_temp_ga) harm contrib
          !     ^                               ^
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,1,k2dr)  &         ! 1st term (radang_temp_ga) gamma contrib
                              -   GtFvec(:,1,k_gr)*GtGa(:,ma,mb,k2dr) &         ! 3rd term Gt hamonic contrib
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr) &         ! 3rd term Gt gamma contib + second term below
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr) & ! gamma radangF
          !              ^                                            ^         ! sign of (a) deriv
                         - radang_temp_ga(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! 2nd solharm radangF
                                   - Ga(:,k_gr)*two_fact2_vec(:,k2dr) )         ! all GA overlap contribs

!------------------------------------------------------------------------------
!!                    -radang_temp_gadg(:,ma,mb,k_gr)*GtFvec(:,1,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
!!                                   - two_fact2_vec(:,k_gr)*Ga(:,k2dr) & ! 2nd term radang_temp contribs,sign of radangF cont
!!              + radang_temp(:,ma,mb)*two_fact2_vec(:,k_gr)*two_fact2_vec(:,k2dr) &  ! ovl radangF
!!               GGa(:,k_gr,k2dr)= radang_temp_gga(:,ma,mb,k_gr,k2dr) + &
!!                                 GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr)*GGt(:,ma,mb)
!-----------------------------------------------------------------------------

              GGb(:,k_gr,k2dr)= ( radang_temp_ggb(:,ma,mb,k_gr,k2dr)  &         ! 1st term (radang_temp_ga) harm contrib
                                              ! ^                                 (b) derivative
                                -GtGb(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  &         ! 1st term (radang_temp_ga) gamma contrib
                                 !  ^                     !  ^                    (b) derivative
                              -   GtFvec(:,2,k_gr)*GtGb(:,ma,mb,k2dr) &         ! 3rd term Gt hamonic contrib
                                         ! ^          ^                           (b) derivative
                     + GGt(:,ma,mb)*GtFvec(:,2,k_gr)*GtFvec(:,2,k2dr) &         ! 3rd term Gt gamma contib + second term below
                                          !  ^                ^                   (b) derivative
                         - Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & ! gamma radangF
                     !   ^                                            ^         ! 2nd term with different signs in Ga and Gb
                         + radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! 2nd sol-harm radangF
                     !   ^              ^                                       ! sign of (b) deriv as in Gb
                                   + Gb(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs
                                 ! ^  ^                                 ! sign changed because (b) overl derv
                                                                        ! should be equal to -(a) derv

              !----------------------------------------------------------
              GaGb(:,k_gr,k2dr)= ( radang_temp_gagb(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                 !              ^ ^
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                                 !  ^                        ^
                              -   GtFvec(:,1,k_gr)*GtGb(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                                 !         ^a         ^
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,2,k2dr) & ! 3rd term Gt gamma contib + second term below
                                 !           ^a               ^b
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & !2nd term gamma radangF
                     !   ^a-conrib sign                               ^b
                         - radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                     !   ^a-contr sign  ^
                                   + Ga(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs
                                !  ^b ^    sign of b-cotrib
                     !-----------------------------------------

             enddo
             enddo

!----------------------------------------------------------------
!
!           exist if grad and dervs indexes coincide

            ! now go two_fact2_vec and GtFvec prefactor contribs

           do k2dr=1,3

            !    2nd term  (radangF) prefactor contribs
            !    ( derivatives by two_fact2_vec fact2*2* (a-b) )
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib, sign of radangF contrib
            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !                                ^  (sign of radangF contrib)*(sign (-) of b in (a-b) argument)

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)+radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !^ ^               ^ ^             ^ ((-) sing of a-contrib)x((-) sign of (a-b) argument)

            !    comment to check Gt contribution whan harmonic and overlap fixed

            !   g4 - (d/da GtFvec )xGt
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1)* &
                                                         aexp_arr(:)/fact0(:)
            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-Gt(:,ma,mb)*two_bexp_arr(:)*exp_arg(:,1)  &
                                                        !     ^
                                                        *bexp_arr(:)/fact0(:)
                                                        !^
                                                        !-----------------------

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1) &
                                                        !       ^
                                                        *bexp_arr(:)/fact0(:)
                                                        !^
                                                        !----------------------
           enddo

              dervsM(:,mb,ma,1:3,1:3)=dervsM(:,mb,ma,1:3,1:3)+GGa(:,:,:)
              dervsM(:,mb,ma,4:6,4:6)=dervsM(:,mb,ma,4:6,4:6)+GGb(:,:,:)
              dervsM(:,mb,ma,1:3,4:6)=dervsM(:,mb,ma,1:3,4:6)+GaGb(:,:,:)
              do k_gr=1,3 !(6)
              do k2dr=1,3
               dervsM(:,mb,ma,k2dr+3,k_gr)=dervsM(:,mb,ma,k2dr+3,k_gr)+GaGb(:,k_gr,k2dr)
              enddo
              enddo

      call fill_ca_dervsM(ma,mb)
!      cx_rot_dervs: if(do_rotation_eq(j)) then
!      if(moving_a) then
!       do i_grad=1,grad_dim
!!         ca_dervsM(:,mb,ma,i_grad,1:3)=ca_dervsM(:,mb,ma,i_grad,1:3) &
!!              -rotmat_eq(i_grad,1,j)*(GGa(:,1,:)+GaGb(:,1,:)) &
!!              -rotmat_eq(i_grad,2,j)*(GGa(:,2,:)+GaGb(:,2,:)) &
!!              -rotmat_eq(i_grad,3,j)*(GGa(:,3,:)+GaGb(:,3,:))
!         do k2dr=1,3   ! transposed GaGb, intended to (7:9,1:3) block
!          ca_dervsM(:,mb,ma,i_grad,k2dr)=ca_dervsM(:,mb,ma,i_grad,k2dr) &
!              -rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GaGb(:,k2dr,1)) &
!              -rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GaGb(:,k2dr,2)) &
!              -rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GaGb(:,k2dr,3))
!         enddo
!       enddo
!      endif
!
!      if(moving_b) then
!       do i_grad=1,grad_dim
!!        do k2dr=1,3   !(7)
!!        ca_dervsM(:,mb,ma,i_grad,3+k2dr)=ca_dervsM(:,mb,ma,i_grad,3+k2dr) &
!!              -rotmat_eq(i_grad,1,j)*(GGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
!!              -rotmat_eq(i_grad,2,j)*(GGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
!!              -rotmat_eq(i_grad,3,j)*(GGb(:,3,k2dr)+GaGb(:,k2dr,3))
!!         enddo
!        ca_dervsM(:,mb,ma,i_grad,4:6)=ca_dervsM(:,mb,ma,i_grad,4:6) &
!              -rotmat_eq(i_grad,1,j)*(GGb(:,1,:)+GaGb(:,1,:)) &
!              -rotmat_eq(i_grad,2,j)*(GGb(:,2,:)+GaGb(:,2,:)) &
!              -rotmat_eq(i_grad,3,j)*(GGb(:,3,:)+GaGb(:,3,:))
!       enddo
!      endif
!
!      else if(grad_dim.gt.0) then
!        if(moving_a) &
!        ca_dervsM(:,mb,ma,:,1:3)=ca_dervsM(:,mb,ma,:,1:3)-(GGa(:,:,:)+GaGb(:,:,:))
!        if(moving_b) then !(8)
!         do k_gr=1,3
!          do k2dr=1,3
!           ca_dervsM(:,mb,ma,k_gr,3+k2dr)= &
!             ca_dervsM(:,mb,ma,k_gr,3+k2dr)-(GGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
!          enddo
!         enddo
!        endif
!      endif cx_rot_dervs

          ! *** symetrize by i_grad index and store result in temp

      rot_dervs1: if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim              ! ! no sum by j to avoid double couning, see rot_dervs2
         do k2dr=1,3    !(9)
         dervs_totsymM_temp(:,mb,ma,i_grad,k2dr)= & !!! dervs_totsymM_temp(:,mb,ma,i_grad,k2dr) &
              +rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GGb(:,1,k2dr)+GaGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
              +rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GGb(:,2,k2dr)+GaGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
              +rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GGb(:,3,k2dr)+GaGb(:,3,k2dr)+GaGb(:,k2dr,3))
         enddo
       enddo

      else if(grad_dim.gt.0) then
       do k_gr=1,3   !(10)
        do k2dr=1,3
        dervs_totsymM_temp(:,mb,ma,k_gr,k2dr)=dervs_totsymM_temp(:,mb,ma,k_gr,k2dr) &
        +(GGa(:,k_gr,k2dr)+GGb(:,k_gr,k2dr)+GaGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
        enddo
       enddo
      endif rot_dervs1

          ! symetrize by k_gr index using temp of previous step and remap to grad_totsymM
       rot_dervs2: if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim
         dervs_totsymM(:,mb,ma,:,i_grad)=dervs_totsymM(:,mb,ma,:,i_grad) &
            +rotmat_eq(i_grad,1,j)*dervs_totsymM_temp(:,mb,ma,:,1) &
            +rotmat_eq(i_grad,2,j)*dervs_totsymM_temp(:,mb,ma,:,2) &
            +rotmat_eq(i_grad,3,j)*dervs_totsymM_temp(:,mb,ma,:,3)
       enddo

       else if(grad_dim.gt.0) then
          dervs_totsymM(:,mb,ma,:,:)=dervs_totsymM(:,mb,ma,:,:) &
              +dervs_totsymM_temp(:,mb,ma,:,:)
       endif rot_dervs2

!!----------------------------------------------------------------------------------------------
!!            !    2nd term  (radangF) prefactor contribs
!!            !    ( derivatives by two_fact2_vec fact2*2* (a-b) )
!!            GGa(:,1,1)=GGa(:,1,1)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
!!            ! ^          ^       ^                                   ! for (a) sign of radangF contrib
!!            GGa(:,2,2)=GGa(:,2,2)-radang_temp(:,ma,mb)*two*fact2(:)
!!            GGa(:,3,3)=GGa(:,3,3)-radang_temp(:,ma,mb)*two*fact2(:)
!
!!    commented to check Gt contribution whan harmonic and overlap fixed
!
!!            !   g4 - d/da GtFvec
!!            ! here the signs of contributions changed for consistency check
!!            GGa(:,1,1)=GGa(:,1,1)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1)* &
!!                                                  aexp_arr(:)/fact0(:)
!!            GGa(:,2,2)=GGa(:,2,2)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1)* &
!!                                                  aexp_arr(:)/fact0(:)
!!            GGa(:,3,3)=GGa(:,3,3)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1)* &
!!                                                  aexp_arr(:)/fact0(:)
!!---------------------------------------------
!
!!            !    2nd term  (radangF) prefactor contribs
!!            !    ( derivatives by two_fact2_vec fact2*2* (a-b) )
!!            GGb(:,1,1)=GGb(:,1,1)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
!!            ! ^          ^       ^                                   ! for (b) sign of 2nd contrib in Gb changed  (+ -> -)
!!            GGb(:,2,2)=GGb(:,2,2)-radang_temp(:,ma,mb)*two*fact2(:)
!!            GGb(:,3,3)=GGb(:,3,3)-radang_temp(:,ma,mb)*two*fact2(:)
!!!             ^          ^   !   ^   *** probable for (b) derivative sign should be oposite to (a) contrib
!!
!!
!!            ! ***    g4 - d/db GtFvec
!
!!            GGb(:,1,1)=GGb(:,1,1)-Gt(:,ma,mb)*two_bexp_arr(:)*exp_arg(:,1)* &
!!!                                                 ^
!!                                                         bexp_arr(:)/fact0(:)
!!                                                       ! ^
!!            GGb(:,2,2)=GGb(:,2,2)-Gt(:,ma,mb)*two_bexp_arr(:)*exp_arg(:,1)* &
!                                                         bexp_arr(:)/fact0(:)
!!            GGb(:,3,3)=GGb(:,3,3)-Gt(:,ma,mb)*two_bexp_arr(:)*exp_arg(:,1)* &
!!                                                         bexp_arr(:)/fact0(:)
!!!---------------------------------------------------------------
!!---------------------------------------------------------------------------------------


          gen_quad: if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i) then

!----------------------------------------------------------
  ! testing final GGa and GGb contributions

!            testing_a: if(k.eq.1.and.i.eq.1.and.abs(sum(GGa(:,1,1))).gt.1.0e-10) then
!              print*,'i,lc,ma,mb, Ga GGa',i,l_c,mc,ma,mb, sum(Ga(:,1)),sum(GGa(:,1,1))


!!----------------------------------
!   testing partial contribs
!!
!!              print*,'radangF all', -sum(radang_temp(:,ma,mb)*two_fact2_vec(:,1)), &
!!                                    sum(Gt(:,ma,mb)*two_fact2_vec(:,1)*GtFvec(:,1,1) & ! gamma radangF
!!                         - radang_temp_ga(:,ma,mb,1)*two_fact2_vec(:,1) & ! harm radangF
!!                         +radang_temp(:,ma,mb)*two_fact2_vec(:,1)*two_fact2_vec(:,1) &
!!                         -radang_temp(:,ma,mb)*two*fact2(:))
!
!!              print*,'GtFvec ',sum(GtFvec(:,1,1)),sum(two_aexp_arr(:)*exp_arg(:,1)* &
!!                                                         aexp_arr(:)/fact0(:))

!  endif testing_a

!            if(k.eq.1.and.i.eq.1.and.abs(sum(GGb(:,1,1))).gt.1.0e-10) &
!              print*,'i,lc,ma,mb, Gb GGb',i,l_c,mc,ma,mb, sum(Gb(:,1)),sum(GGb(:,1,1))

!            if(k.eq.1.and.i.eq.1.and.abs(sum(GaGb(:,1,1))).gt.1.0e-10) &
!              print*,'i,lc,ma,mb, Ga GaGb',i,l_c,mc,ma,mb, sum(Ga(:,1)),sum(GaGb(:,1,1))

!  done
!--------------------------------------

           endif gen_quad
          endif sd2

    if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim ! only if moving_c

          grad_totsymM(:,mb,ma,i_grad)=grad_totsymM(:,mb,ma,i_grad) &
              -rotmat_eq(i_grad,1,j)*(ga(:,1)+gb(:,1)) &
              -rotmat_eq(i_grad,2,j)*(ga(:,2)+gb(:,2)) &
              -rotmat_eq(i_grad,3,j)*(ga(:,3)+gb(:,3))
       enddo

    elseif(grad_dim.gt.0) then
         grad_totsymM(:,mb,ma,1)=grad_totsymM(:,mb,ma,1)-(ga(:,1)+gb(:,1))
         grad_totsymM(:,mb,ma,2)=grad_totsymM(:,mb,ma,2)-(ga(:,2)+gb(:,2))
         grad_totsymM(:,mb,ma,3)=grad_totsymM(:,mb,ma,3)-(ga(:,3)+gb(:,3))
    end if
              end do
            end do

      else moving

            do ma=1,2*L_a+1
              do mb=1,2*L_b+1
               do k_gr=1,3

              Ga(:,k_gr)= ( radang_temp_ga(:,ma,mb,k_gr) &
                -radang_temp(:,ma,mb)*fact2*two*(xa(k_gr)-xb(k_gr)) &
                -two*aexp_arr(:)*Gt(:,ma,mb) &
                 *(gamma_arg(:,k_gr)-unique_atoms(i)%position(k_gr,j)) &
                                            *exp_arg(:,1) )/c_exp_arg2(:)
              Gb(:,k_gr)= ( radang_temp_gb(:,ma,mb,k_gr) &
                +radang_temp(:,ma,mb)*fact2*two*(xa(k_gr)-xb(k_gr)) &
                -two*bexp_arr(:)*Gt(:,ma,mb) &
                 *(gamma_arg(:,k_gr)-unique_atoms(i)%position(k_gr,j)) &
                                            *exp_arg(:,1) )/c_exp_arg2(:)
               if(moving_a) &
                coul_int_grad(k_gr)%l(L_c)%m(:,k,i_ind,mb,ma)=  &
                        coul_int_grad(k_gr)%l(L_c)%m(:,k,i_ind,mb,ma)+Ga(:,k_gr)
               if(moving_b) &
                coul_int_grad(k_gr+3)%l(L_c)%m(:,k,i_ind,mb,ma)= &
                      coul_int_grad(k_gr+3)%l(L_c)%m(:,k,i_ind,mb,ma)+Gb(:,k_gr)
             enddo
    if(moving_c) then
    if(do_rotation) then
       ! make c gradient totalsymmetric before adding
       do i_grad=1,grad_dim ! only if moving_c
          coul_int_grad_totsym(i_grad)%l(l_c)%m(:,k,i_ind,mb,ma) = &
               coul_int_grad_totsym(i_grad)%l(l_c)%m(:,k,i_ind,mb,ma) - &
               rotmat(i_grad,1)*(ga(:,1)+gb(:,1)) - &
               rotmat(i_grad,2)*(ga(:,2)+gb(:,2)) - &
               rotmat(i_grad,3)*(ga(:,3)+gb(:,3))
       enddo
    elseif(grad_dim.gt.0) then
       coul_int_grad_totsym(1)%l(l_c)%m(:,k,i_ind,mb,ma)= &
                     coul_int_grad_totsym(1)%l(l_c)%m(:,k,i_ind,mb,ma)-(ga(:,1)+gb(:,1))
       coul_int_grad_totsym(2)%l(l_c)%m(:,k,i_ind,mb,ma)= &
                     coul_int_grad_totsym(2)%l(l_c)%m(:,k,i_ind,mb,ma)-(ga(:,2)+gb(:,2))
       coul_int_grad_totsym(3)%l(l_c)%m(:,k,i_ind,mb,ma)= &
                     coul_int_grad_totsym(3)%l(l_c)%m(:,k,i_ind,mb,ma)-(ga(:,3)+gb(:,3))
    end if !do_rotation
    end if ! moving_c
              end do
            end do
   endif moving
  endif gr

 else gr_els ! 1
            do ma=1,2*L_a+1
              do mb=1,2*L_b+1
        pointer_coul(:,k,i_ind,mb,ma)= pointer_coul(:,k,i_ind,mb,ma)    &
                                     + radang_temp(:,ma,mb)/c_exp_arg2(:)
              end do
            end do
  endif gr_els

         enddo contributing

        gr2cograd: if(integralpar_gradients) then
            do ma=1,2*L_a+1
              do mb=1,2*L_b+1

               do k_gr=1,6
               coul_int_grad(  k_gr)%l(L_c)%m(:,k,i_ind,mb,ma)=gradM(:,mb,ma,k_gr)/c_exp_arg2(:)
!              coul_int_grad(k_gr+3)%l(L_c)%m(:,k,i_ind,mb,ma)=gradM(:,mb,ma,k_gr+3)/c_exp_arg2(:)
               enddo

         fmap_ab_dervs: if(integralpar_dervs) then ! if false no ab lc derivative 1st of 3 contribs
           do k_gr=1,6  ! A x B block of derivatives
           do k2dr=1,6
            coul_int_dervs(k_gr,k2dr)%l(Lc)%m(:,k,i_ind,mb,ma)= &
              coul_int_dervs(k_gr,k2dr)%l(Lc)%m(:,k,i_ind,mb,ma) &
                              +dervsM(:,mb,ma,k_gr,k2dr)/c_exp_arg2(:)
           enddo
           enddo
         endif fmap_ab_dervs

       do i_grad=1,grad_dim            ! only if moving_c
          coul_int_grad_totsym(i_grad)%l(l_c)%m(:,k,i_ind,mb,ma) =  &
                 grad_totsymM(:,mb,ma,i_grad)/c_exp_arg2(:)

            fmap_cc_and_ca_dervs: if(.true..and.integralpar_dervs) then
             ! if false no cc and ca lccase contins

             do k2dr=1,grad_dim
              coul_int_dervs_totsym(i_grad,k2dr)%l(L_c)%m(:,k,i_ind,mb,ma)= &
                              dervs_totsymM(:,mb,ma,i_grad,k2dr)/c_exp_arg2(:)
             enddo

             do k2dr=1,6
              coul_int_ca_dervs(i_grad,k2dr)%l(L_c)%m(:,k,i_ind,mb,ma)= &
                                  ca_dervsM(:,mb,ma,i_grad,k2dr)/c_exp_arg2(:)
             enddo

            endif fmap_cc_and_ca_dervs

       enddo

       enddo
     enddo
    endif gr2cograd

   enddo independents

  end subroutine Radial_Ang_3cSA



    subroutine Radial_Ang_3cS(radial3cmat,la,lb,lc)
    ! used in nolc loop


    ! RESULT:
    !            coul_int_grad            - a and b coul grads
    !            coul_int_grad_totsym     - symmetrized c gradients

    ! local temps:
    !           (1)  gradM - (a) and (b) coul gradients
    !           (2)  grad_totsymM - symmetrized (c) gradients
    !           (3)  dervsM - (ab) block coul dervs

    ! history:
    ! now extended with GaGb derivs
    ! dervsM added




          integer(kind=i4_kind),intent(in):: lc,la,lb
           real(kind=r8_kind), intent(in)::radial3cmat(num,n_jj,n_equals)
!         real(kind=r8_kind), intent(inout)::radial3cmat(num,n_jj,n_equals) !!!temporary to fix this contrib

          integer(kind=i4_kind):: k_gr,k2dr,i_grad
!!!          real(kind=r8_kind):: dervs_cc(num,2*la+1,2*lb+1,3,3)

      if(integralpar_gradients) then
        grad_totsymM=0.0_r8_kind
        gradM=0.0_r8_kind

          if(integralpar_dervs) then
           dervsM=0.0_r8_kind
           dervs_totsymM=0.0_r8_kind
           dervs_totsymM_temp=0.0_r8_kind
           ca_dervsM=0.0_r8_kind
          endif

      endif

     equals: do j=1,n_equals

           radang_temp=0.0_r8_kind

            call dgemm('n','n',num,(2*La+1)*(2*Lb+1),n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                        AngR3cMat(:,:,:,:,j),n_jj, &
                        1.0_r8_kind,radang_temp(:,:,:),num)

   grsco: if(integralpar_3cob_grad) then
    gr:      if(radang_sco_grad) then

!----------------------------------------------------------------------------------
!             do ma=1,2*La+1  !!!!!
!              do mb=1,2*Lb+1
!               pointer_coul(:,k,1,mb,ma)= pointer_coul(:,k,1,mb,ma)    &
!                                             + radang_temp(:,ma,mb)/c_exp_arg2(:)
!              end do
!             end do   !!!!!

!           if(present(equala)) then
!           check_bc=(quadrupel%ua2.eq.i).and.(equalb.eq.j)
!           if(check_ab.and.check_bc) cycle
!           endif
!--------------------------------------------------------------


!            radial3cmat(:,:,j)=1.0_r8_kind  !!!! now it is fixed to check radang_temp_ggb

            call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*3,n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                        AngR3cMatGa(:,:,:,:,1,j),n_jj, &
                        0.0_r8_kind,radang_temp_ga(:,:,:,:),num)

            call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*3,n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                        AngR3cMatGb(:,:,:,:,1,j),n_jj, &
                        0.0_r8_kind,radang_temp_gb(:,:,:,:),num)

          dervs: if(integralpar_dervs) then

             ! 2nd deriv radang_temp_ga harmonic contrib (1)

           call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*9,n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                       ! ^ at diff to SA code above
                        AngR3cMatGGa(:,:,:,:,:,1,j),n_jj, &
                                              !^ at diff to SA code above
                        0.0_r8_kind,radang_temp_gga(:,:,:,:,:),num)


           call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*9,n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                       ! ^ at diff to SA code above
                        AngR3cMatGGb(:,:,:,:,:,1,j),n_jj, &
                       !           ^           ^
                        0.0_r8_kind,radang_temp_ggb(:,:,:,:,:),num)
! if(k.eq.1.and.i.eq.3) print*,sum(AngR3cMatGb(:,:,:,1,1,j)),sum(AngR3cMatGGb(:,:,:,1,1,1,j)),'AngR3cMatGb-AngR3cMatGGb'

           call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*9,n_jj, &
                       1.0_r8_kind, radial3cmat(:,:,j),num,&
                       ! ^ at diff to SA code above
                        AngR3cMatGaGb(:,:,:,:,:,1,j),n_jj, &
                       !          ^ ^           ^
                        0.0_r8_kind,radang_temp_gagb(:,:,:,:,:),num)
                       !                         ^ ^


             ! to check radang_temp_gga  gamma and overlap contribs
             ! have to be fixed

           call dgemm('n','n',num,(2*La+1)*(2*Lb+1),n_jj, &
                       1.0_r8_kind, radial3cmatGG(:,:,j),num,&
                        AngR3cMat(:,:,:,1,j),n_jj, &
                        0.0_r8_kind,GGt(:,:,:),num)  ! second gamma derivative

           ! next find what is derivative by argument
           !                        gamma derivatives


           call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*3,n_jj, &
                       1.0_r8_kind, radial3cmatG(:,:,j),num,&
                        AngR3cMatGa(:,:,:,:,1,j),n_jj, &
                        0.0_r8_kind,GtGa(:,:,:,:),num)
                                    !^ ga spherical harmonic gamma Gt derivative

           call dgemm('n','n',num,(2*La+1)*(2*Lb+1)*3,n_jj, &
                       1.0_r8_kind, radial3cmatG(:,:,j),num,&
                        AngR3cMatGb(:,:,:,:,1,j),n_jj, &
                        0.0_r8_kind,GtGb(:,:,:,:),num)
                                    !^ gb spherical harmonic gamma Gt derivative
          endif dervs

            Gt=0.0_r8_kind
            call dgemm('n','n',num,(2*La+1)*(2*Lb+1),n_jj, &
                       1.0_r8_kind, radial3cmatG(:,:,j),num,&
                        AngR3cMat(:,:,:,1,j),n_jj, &
                        1.0_r8_kind,Gt(:,:,:),num)

         do k_gr=1,3
          GtFvec(:,1,k_gr)= two_aexp_arr(:)*gamma_arg_xc(:,k_gr,j)*exp_arg(:,1)
                                           !^ coordinate dependent part
          GtFvec(:,2,k_gr)= two_bexp_arr(:)*gamma_arg_xc(:,k_gr,j)*exp_arg(:,1)
        enddo

           moving_cases: if(moving_a.and.moving_b.and.moving_c) then
             ma_loop: do ma=1,2*La+1
               do mb=1,2*Lb+1

             Gx_grads: do k_gr=1,3

             radangF(:)=radang_temp(:,ma,mb)*two_fact2_vec(:,k_gr)

               Ga(:,k_gr)= ( radang_temp_ga(:,ma,mb,k_gr)-radangF(:) &
                 -GtFvec(:,1,k_gr)*Gt(:,ma,mb) )!!!/c_exp_arg2(:)

               Gb(:,k_gr)= ( radang_temp_gb(:,ma,mb,k_gr)+radangF(:) &
                 -GtFvec(:,2,k_gr)*Gt(:,ma,mb) )!!!/c_exp_arg2(:)

                ! a & b gradients
                gradM(:,mb,ma,k_gr)=  gradM(:,mb,ma,k_gr)  +Ga(:,k_gr)
                gradM(:,mb,ma,k_gr+3)=gradM(:,mb,ma,k_gr+3)+Gb(:,k_gr)

              enddo Gx_grads

   sd2: if(integralpar_dervs) then
          !   we have one overlap contrib for all terms but
          !   separate harmonic and gamma contribs
          !            + prefacror contribs
             do k2dr=1,3
             do k_gr=1,3
              GGa(:,k_gr,k2dr)= ( radang_temp_gga(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,1,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                              -   GtFvec(:,1,k_gr)*GtGa(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr) & ! 3rd term Gt gamma contib + second term below
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr) & ! gamma radangF
                         - radang_temp_ga(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                                   - Ga(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs



              GGb(:,k_gr,k2dr)= ( radang_temp_ggb(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                -GtGb(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                              -   GtFvec(:,2,k_gr)*GtGb(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                     + GGt(:,ma,mb)*GtFvec(:,2,k_gr)*GtFvec(:,2,k2dr) & ! 3rd term Gt gamma contib + second term below
                         - Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & ! gamma radangF
                         + radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                                   + Gb(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs

!           to check radang_temp_ggb one fixes radial3c contib to one




              GaGb(:,k_gr,k2dr)= ( radang_temp_gagb(:,ma,mb,k_gr,k2dr)  & ! 1st term (radang_temp_ga) harm contrib
                                 !              ^ ^
                                -GtGa(:,ma,mb,k_gr)*GtFvec(:,2,k2dr)  & ! 1st term (radang_temp_ga) gamma contrib
                                 !  ^                        ^
                              -   GtFvec(:,1,k_gr)*GtGb(:,ma,mb,k2dr) & ! 3rd term Gt hamonic contrib
                                 !         ^a         ^
                     + GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,2,k2dr) & ! 3rd term Gt gamma contib + second term below
                                 !           ^a               ^b
                         + Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,2,k2dr) & !2nd term gamma radangF
                     !   ^a-conrib sign                               ^b
                         - radang_temp_gb(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr) & ! harm radangF
                     !   ^a-contr sign  ^
                                   + Ga(:,k_gr)*two_fact2_vec(:,k2dr) ) ! all GA overlap contribs
                                !  ^b ^    sign of b-cotrib
                     !-----------------------------------------

            !    comment to check Gt contribution whan harmonic and overlap fixed
             enddo ! k_gr

            !   g4 - (d/da GtFvec )xGt
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1)* &
                                                         aexp_arr(:)/fact0(:)

            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-Gt(:,ma,mb)*two_bexp_arr(:)*exp_arg(:,1)  &
                                                        *bexp_arr(:)/fact0(:)

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)-Gt(:,ma,mb)*two_aexp_arr(:)*exp_arg(:,1) &
                                                        !       ^
                                                        *bexp_arr(:)/fact0(:)
                                                        !^
                                                        !----------------------

!--------------------------------------------------------------------
            ! now go two_fact2_vec and GtFvec prefactor contribs


            !    2nd term  (radangF) prefactor contribs
            !    ( derivatives by two_fact2_vec fact2*2* (a-b) )
            GGa(:,k2dr,k2dr)=GGa(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib, sign of radangF contrib


            GGb(:,k2dr,k2dr)=GGb(:,k2dr,k2dr)-radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !                                ^  (sign of radangF contrib)*(sign (-) of b in (a-b) argument)

            GaGb(:,k2dr,k2dr)=GaGb(:,k2dr,k2dr)+radang_temp(:,ma,mb)*two*fact2(:)  !last 4th radanfF contrib
            !^ ^               ^ ^             ^ ((-) sing of a-contrib)x((-) sign of (a-b) argument)

           enddo ! k2dr

          if(quadrupel%ua1.ne.quadrupel%ua2.and.quadrupel%ua2.ne.i.and.quadrupel%ua1.ne.i) then
          ! all 3 indexes are different and it is easy to check with finite difference  calc

!            if(k.eq.1.and.i.eq.1.and.abs(sum(GGa(:,1,1))).gt.1.0e-10) &
!              if(k.eq.1) &
!              print*,'i,lc,ma,mb, Ga GGa',i,lc,ma,mb, sum(Ga(:,1)),sum(GGa(:,1,1)),' nolc'
              ! here 1 is x for cenrer (a)

!              if(k.eq.1) &
!              print*,'i,lc,ma,mb, radang_temp_ga radang_temp_gga(:,ma,mb,1,1)', &
!                   sum(radang_temp_ga(:,ma,mb,1)), sum(radang_temp_gga(:,ma,mb,1,1))

!            if(k.eq.1.and.i.eq.1.and.abs(sum(GGb(:,1,1))).gt.1.0e-10) &
!              if(k.eq.1) &
!              print*,'i,lc,ma,mb, Gb GGb',i,lc,ma,mb, sum(Gb(:,1)),sum(GGb(:,1,1)),' nolc'
              ! here 1 is x for cenrer (b)

!              if(k.eq.1) &
!              print*,' radang_temp_gb radang_temp_ggb(:,ma,mb,1,1)', &
!                   sum(radang_temp_gb(:,ma,mb,1)), sum(radang_temp_ggb(:,ma,mb,1,1))


!            if(k.eq.1.and.i.eq.1.and.abs(sum(GaGb(:,1,1))).gt.1.0e-10) &
!              if(k.eq.1) &
!              print*,'i,lc,ma,mb, Ga GaGb',i,lc,mc,ma,mb, sum(Ga(:,1)),sum(GaGb(:,1,1)),' nolc'
           endif

!        if(k.eq.1) then
!          print*,'sum components',&! sum(radang_temp_gga(:,ma,mb,:,:)), &
!                  sum( -GtGa(:,ma,mb,k_gr)),& !*GtFvec(:,1,k2dr)), &
!                  sum( - GtFvec(:,1,k_gr)*GtGa(:,ma,mb,k2dr)), &
!                  sum( GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr)), &
!                  sum(  Gt(:,ma,mb)),& !*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr)), &
!                  sum(  - radang_temp_ga(:,ma,mb,k2dr)),& !*two_fact2_vec(:,k_gr)), &
!                  sum(            - Ga(:,k_gr)) !*two_fact2_vec(:,k2dr) )
!        endif

!         do k_gr=1,3
!         do k2dr=1,3
!          if(abs(sum(GGa(:,k_gr,k2dr))-sum(GGa(:,k2dr,k_gr))).gt.1.e-6) then
!           print*,k_gr,k2dr,'k_gr,k2dr GGa 12 21 ',sum(GGa(:,k_gr,k2dr)),sum(GGa(:,k2dr,k_gr)),ma,mb
!           print*,sum(radang_temp_gga(:,ma,mb,k_gr,k2dr)), &
!                  sum( -GtGa(:,ma,mb,k_gr)*GtFvec(:,1,k2dr)), &
!                  sum( - GtFvec(:,1,k_gr)*GtGa(:,ma,mb,k2dr)), &
!                  sum( GGt(:,ma,mb)*GtFvec(:,1,k_gr)*GtFvec(:,1,k2dr)), &
!                  sum(  Gt(:,ma,mb)*two_fact2_vec(:,k_gr)*GtFvec(:,1,k2dr)), &
!                  sum(  - radang_temp_ga(:,ma,mb,k2dr)*two_fact2_vec(:,k_gr)), &
!                  sum(            - Ga(:,k_gr)*two_fact2_vec(:,k2dr) )
!          print*,sum(GGt(:,ma,mb)),sum(GtFvec(:,1,k_gr)),sum(GtFvec(:,1,k2dr)),' term 4', &
!              lc,sum(radial3cmatGG(:,:,j))
!
!          endif
!         enddo
!         enddo



           ! GGa GGb GaGb derivatives calculated, other can found combining this 3 matrices
           !  Gc = -Ga-Gb -> (d2/dcda) co= - GGa -GaGb
           !  GcGa=GaGc=-GGa-GaGb
           !  GcGb=GbGc=-GGb-GaGb
           !  GcGc=-GaGc-GbGc=GGa+GGb+2GaGb


              dervsM(:,mb,ma,1:3,1:3)=dervsM(:,mb,ma,1:3,1:3)+GGa(:,:,:)
              dervsM(:,mb,ma,4:6,4:6)=dervsM(:,mb,ma,4:6,4:6)+GGb(:,:,:)
              dervsM(:,mb,ma,1:3,4:6)=dervsM(:,mb,ma,1:3,4:6)+GaGb(:,:,:)
              do k_gr=1,3 !(11)
              do k2dr=1,3
               dervsM(:,mb,ma,k2dr+3,k_gr)=dervsM(:,mb,ma,k2dr+3,k_gr)+GaGb(:,k_gr,k2dr)
              enddo
              enddo


      call fill_ca_dervsM(ma,mb)

!      cx_rot_dervs: if(do_rotation_eq(j)) then
!      if(moving_a) then
!       do i_grad=1,grad_dim
!!       ca_dervsM(:,mb,ma,i_grad,1:3)=ca_dervsM(:,mb,ma,i_grad,1:3) &
!!              -rotmat_eq(i_grad,1,j)*(GGa(:,1,:)+GaGb(:,1,:)) &
!!              -rotmat_eq(i_grad,2,j)*(GGa(:,2,:)+GaGb(:,2,:)) &
!!              -rotmat_eq(i_grad,3,j)*(GGa(:,3,:)+GaGb(:,3,:))
!         do k2dr=1,3   ! transposed GaGb, intended to (7:9,1:3) block
!          ca_dervsM(:,mb,ma,i_grad,k2dr)=ca_dervsM(:,mb,ma,i_grad,k2dr) &
!              -rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GaGb(:,k2dr,1)) &
!              -rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GaGb(:,k2dr,2)) &
!              -rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GaGb(:,k2dr,3))
!         enddo
!       enddo
!      endif
!
!      if(moving_b) then
!       do i_grad=1,grad_dim
!
!!        do k2dr=1,3 !(12)
!!         ca_dervsM(:,mb,ma,i_grad,3+k2dr)=ca_dervsM(:,mb,ma,i_grad,3+k2dr) &
!!              -rotmat_eq(i_grad,1,j)*(GGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
!!              -rotmat_eq(i_grad,2,j)*(GGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
!!              -rotmat_eq(i_grad,3,j)*(GGb(:,3,k2dr)+GaGb(:,k2dr,3))
!!         enddo
!        ca_dervsM(:,mb,ma,i_grad,4:6)=ca_dervsM(:,mb,ma,i_grad,4:6) &
!              -rotmat_eq(i_grad,1,j)*(GGb(:,1,:)+GaGb(:,1,:)) &
!              -rotmat_eq(i_grad,2,j)*(GGb(:,2,:)+GaGb(:,2,:)) &
!              -rotmat_eq(i_grad,3,j)*(GGb(:,3,:)+GaGb(:,3,:))
!       enddo
!      endif
!
!      else if(grad_dim.gt.0) then
!        if(moving_a) &
!        ca_dervsM(:,mb,ma,:,1:3)=ca_dervsM(:,mb,ma,:,1:3)-(GGa(:,:,:)+GaGb(:,:,:))
!
!        if(moving_b) then   !(13)
!         do k_gr=1,3
!          do k2dr=1,3
!           ca_dervsM(:,mb,ma,k_gr,3+k2dr)= &
!             ca_dervsM(:,mb,ma,k_gr,3+k2dr)-(GGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
!          enddo
!         enddo
!        endif
!
!      endif cx_rot_dervs

          ! *** symetrize by i_grad index and store result in temp

      rot_dervs1: if(do_rotation_eq(j)) then
       do k2dr=1,3 !(14)                            ! no sum by j here it is individual contrib
        do i_grad=1,grad_dim                        ! sum by j is in rot_dervs2 block
         dervs_totsymM_temp(:,mb,ma,i_grad,k2dr)= & !!dervs_totsymM_temp(:,mb,ma,i_grad,k2dr) &
              +rotmat_eq(i_grad,1,j)*(GGa(:,1,k2dr)+GGb(:,1,k2dr)+GaGb(:,1,k2dr)+GaGb(:,k2dr,1)) &
              +rotmat_eq(i_grad,2,j)*(GGa(:,2,k2dr)+GGb(:,2,k2dr)+GaGb(:,2,k2dr)+GaGb(:,k2dr,2)) &
              +rotmat_eq(i_grad,3,j)*(GGa(:,3,k2dr)+GGb(:,3,k2dr)+GaGb(:,3,k2dr)+GaGb(:,k2dr,3))
        enddo
       enddo

      else if(grad_dim.gt.0) then
       do k_gr=1,3 !(15)
        do k2dr=1,3
         dervs_totsymM_temp(:,mb,ma,k_gr,k2dr)= & !!!dervs_totsymM_temp(:,mb,ma,k_gr,k2dr) &
        +(GGa(:,k_gr,k2dr)+GGb(:,k_gr,k2dr)+GaGb(:,k_gr,k2dr)+GaGb(:,k2dr,k_gr))
        enddo
       enddo
      endif rot_dervs1

          ! symetrize by k_gr index using temp of previous step and remap to grad_totsymM
       rot_dervs2: if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim
       ! here sum over all j-equal_c contribs, no separate contrib for each j
         dervs_totsymM(:,mb,ma,:,i_grad)= dervs_totsymM(:,mb,ma,:,i_grad) &
            +rotmat_eq(i_grad,1,j)*dervs_totsymM_temp(:,mb,ma,:,1) &
            +rotmat_eq(i_grad,2,j)*dervs_totsymM_temp(:,mb,ma,:,2) &
            +rotmat_eq(i_grad,3,j)*dervs_totsymM_temp(:,mb,ma,:,3)
       enddo

       else if(grad_dim.gt.0) then
          dervs_totsymM(:,mb,ma,:,:)=dervs_totsymM(:,mb,ma,:,:) &
              +dervs_totsymM_temp(:,mb,ma,:,:)
       endif rot_dervs2



     endif sd2

    if(do_rotation_eq(j)) then
       do i_grad=1,grad_dim ! only if moving_c
           grad_totsymM(:,mb,ma,i_grad)=grad_totsymM(:,mb,ma,i_grad) &
              -rotmat_eq(i_grad,1,j)*(ga(:,1)+gb(:,1)) &
              -rotmat_eq(i_grad,2,j)*(ga(:,2)+gb(:,2)) &
              -rotmat_eq(i_grad,3,j)*(ga(:,3)+gb(:,3))
       enddo

    else if(grad_dim.gt.0) then! do_rotation
         grad_totsymM(:,mb,ma,1)=grad_totsymM(:,mb,ma,1)-(ga(:,1)+gb(:,1))
         grad_totsymM(:,mb,ma,2)=grad_totsymM(:,mb,ma,2)-(ga(:,2)+gb(:,2))
         grad_totsymM(:,mb,ma,3)=grad_totsymM(:,mb,ma,3)-(ga(:,3)+gb(:,3))
    endif ! do_rotation

  end do
 end do ma_loop

!     if(k.eq.1.and.integralpar_dervs) then
!       print*,MyID//'dervs_totsymM_temp just calculated unique_c equal_c',i,j,lc
!       do k2dr=1,size(dervs_totsymM,4)
!        print*,sum(dervs_totsymM_temp(:,:,:,k2dr,1)), &
!               sum(dervs_totsymM_temp(:,:,:,k2dr,2)), &
!               sum(dervs_totsymM_temp(:,:,:,k2dr,3))
!       enddo
!      endif

   else moving_cases
             do ma=1,2*La+1
               do mb=1,2*Lb+1
               do k_gr=1,3

               Ga(:,k_gr)= ( radang_temp_ga(:,ma,mb,k_gr) &
                 -radang_temp(:,ma,mb)*fact2*two*(xa(k_gr)-xb(k_gr)) &
                 -two*aexp_arr(:)*Gt(:,ma,mb) &
                  *(gamma_arg(:,k_gr)-unique_atoms(i)%position(k_gr,j)) &
                                            *exp_arg(:,1) )/c_exp_arg2(:)
               Gb(:,k_gr)= ( radang_temp_gb(:,ma,mb,k_gr) &
                 +radang_temp(:,ma,mb)*fact2*two*(xa(k_gr)-xb(k_gr)) &
                 -two*bexp_arr(:)*Gt(:,ma,mb) &
                  *(gamma_arg(:,k_gr)-unique_atoms(i)%position(k_gr,j)) &
                                            *exp_arg(:,1) )/c_exp_arg2(:)
                if(moving_a) &
                 coul_int_grad(k_gr)%l(Lc)%m(:,k,1,mb,ma)= &
                                coul_int_grad(k_gr)%l(Lc)%m(:,k,1,mb,ma)+Ga(:,k_gr)
                if(moving_b) &
                 coul_int_grad(k_gr+3)%l(Lc)%m(:,k,1,mb,ma)= &
                                coul_int_grad(k_gr+3)%l(Lc)%m(:,k,1,mb,ma)+Gb(:,k_gr)
              enddo
    if(moving_c) then
       if(do_rotation_eq(j)) then
       ! make gradient totalsymmetric before adding
       do i_grad=1,grad_dim ! only if moving_c
          coul_int_grad_totsym(i_grad)%l(lc)%m(:,k,1,mb,ma) = &
               coul_int_grad_totsym(i_grad)%l(lc)%m(:,k,1,mb,ma) - &
               rotmat(i_grad,1)*(ga(:,1)+gb(:,1)) - &
               rotmat(i_grad,2)*(ga(:,2)+gb(:,2)) - &
               rotmat(i_grad,3)*(ga(:,3)+gb(:,3))
       enddo
    elseif(grad_dim.gt.0) then ! else do_rotation
       coul_int_grad_totsym(1)%l(lc)%m(:,k,1,mb,ma)= &
                        coul_int_grad_totsym(1)%l(lc)%m(:,k,1,mb,ma)-(ga(:,1)+gb(:,1))
       coul_int_grad_totsym(2)%l(lc)%m(:,k,1,mb,ma)= &
                        coul_int_grad_totsym(2)%l(lc)%m(:,k,1,mb,ma)-(ga(:,2)+gb(:,2))
       coul_int_grad_totsym(3)%l(lc)%m(:,k,1,mb,ma)= &
                        coul_int_grad_totsym(3)%l(lc)%m(:,k,1,mb,ma)-(ga(:,3)+gb(:,3))
     endif ! do_rotation

    endif
              end do
            end do

!       if(k.eq.1.and.i_ind.eq.1) &
!          print*,sum(coul_int_grad(3)%l(Lc)%m(:,1,1,:,:)),' a ',sum(pointer_coul(:,1,1,:,:))
!       if(k.eq.1.and.i_ind.eq.1) &
!          print*,coul_int_grad(6)%l(Lc)%m(1,1,1,1,1),' b ',pointer_coul(1,1,1,1,1)
!       if(k.eq.1.and.i_ind.eq.1) &
!          print*,coul_int_grad_totsym(3)%l(Lc)%m(1,1,1,1,1),' c ',pointer_coul(1,1,1,1,1)

        endif moving_cases

     endif gr

 else grsco
            do ma=1,2*La+1
              do mb=1,2*Lb+1
                pointer_coul(:,k,1,mb,ma)=pointer_coul(:,k,1,mb,ma)+ &
                                          radang_temp(:,ma,mb)/c_exp_arg2(:)
              end do
            end do
 endif grsco

 enddo equals



   fm_grads: if(integralpar_gradients) then

         do ma=1,2*La+1
               do mb=1,2*Lb+1

           do k_gr=1,6
                 coul_int_grad(k_gr)%l(Lc)%m(:,k,1,mb,ma)= &
                  coul_int_grad(k_gr)%l(Lc)%m(:,k,1,mb,ma)+gradM(:,mb,ma,k_gr)/c_exp_arg2(:)

           enddo

         fmap_ab_dervs: if(.true..and.integralpar_dervs) then
           do k_gr=1,6
           do k2dr=1,6
            coul_int_dervs(k_gr,k2dr)%l(Lc)%m(:,k,1,mb,ma)= &
              coul_int_dervs(k_gr,k2dr)%l(Lc)%m(:,k,1,mb,ma) &
                          +dervsM(:,mb,ma,k_gr,k2dr)/c_exp_arg2(:)

           enddo
           enddo
         endif fmap_ab_dervs

            do i_grad=1,grad_dim ! only if moving_c
              coul_int_grad_totsym(i_grad)%l(lc)%m(:,k,1,mb,ma) = &
                                       grad_totsymM(:,mb,ma,i_grad)/c_exp_arg2(:)

            fmap_ca_and_cc_dervs: if(.true..and.integralpar_dervs) then
            !if false no lc=0 fit ca and cc derivatives

             do k2dr=1,grad_dim
              coul_int_dervs_totsym(i_grad,k2dr)%l(Lc)%m(:,k,1,mb,ma)= &
                                  dervs_totsymM(:,mb,ma,i_grad,k2dr)/c_exp_arg2(:)
             enddo

             do k2dr=1,6
              coul_int_ca_dervs(i_grad,k2dr)%l(Lc)%m(:,k,1,mb,ma)= &
                                  ca_dervsM(:,mb,ma,i_grad,k2dr)/c_exp_arg2(:)
             enddo

            endif fmap_ca_and_cc_dervs

            enddo

          enddo
       enddo
   endif fm_grads

 end subroutine Radial_Ang_3cS


    subroutine J_Complex(todo,J_array,L1,L2,L3,arr_name)
      ! Purpose : (De)Allocation of complex J_index pointer
      !           (Complex Angular factor only)
      implicit none
      !------------ Declaration of formal parameters ---------------
      complex(kind=c16_kind), pointer    :: J_array(:,:,:,:,:,:)
      integer(kind=i4_kind),  intent(in) :: L1,L2,L3
      character(len=*),       intent(in) :: todo, arr_name
      !** End of interface *****************************************

      select case(trim(todo))
      case("allocate") !   j1   j2   j3   i1   i2   i3
         allocate(J_array(0:L1,0:L2,0:L3,0:L1,0:L2,0:L3), stat=alloc_stat(80))
         MEMLOG(size(J_array))
         ASSERT(alloc_stat(80).eq.0)
         alloc_stat(80)=1
         J_array = (0.0_r8_kind,0.0_r8_kind)
      case("deallocate")
         deallocate(J_array, stat=alloc_stat(80))
         ASSERT(alloc_stat(80).eq.0)
         MEMLOG(-size(J_array))
      end select
    end subroutine J_Complex

    subroutine J_Real(todo,J_array,L1,L2,L3,arr_name)
      ! Purpose : (De)Allocation of real J_index pointer (common to
      !           Real Angular and Exponential factors
      implicit none
      !------------ Declaration of formal parameters ---------------
      real(kind=r8_kind),     pointer    :: J_array(:,:,:,:,:,:)
      integer(kind=i4_kind),  intent(in) :: L1,L2,L3
      character(len=*),       intent(in) :: todo, arr_name
      !** End of interface *****************************************

      select case(trim(todo))
      case("allocate") !   j1   j2   j3   i1   i2   i3
         allocate(J_array(0:L1,0:L2,0:L3,0:L1,0:L2,0:L3),stat=alloc_stat(80))
         MEMLOG(size(J_array))
         ASSERT(alloc_stat(80).eq.0)
         alloc_stat(80)=1
         J_array = 0.0_r8_kind
      case("deallocate")
         MEMLOG(-size(J_array))
         deallocate(J_array, stat=alloc_stat(80))
         ASSERT(alloc_stat(80).eq.0)
      end select

    end subroutine J_Real

  end module calc_3center_module



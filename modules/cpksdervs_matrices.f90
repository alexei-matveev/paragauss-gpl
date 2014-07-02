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
 module cpksdervs_matrices
 use type_module
 use occupied_levels_module

  implicit none
  private         ! by default, all names are private
  save

 integer(kind=i4_kind), public:: cpks_coul_grad_size=0, max_cpks_coul_grad_size=0
 integer(kind=i4_kind), public:: ol_dervs_size=0, max_ol_dervs_size=0
 integer(kind=i4_kind), public:: size_coul_dervs=0,max_size_coul_dervs=0
 integer(kind=i4_kind), public:: co_ai_size=0
 integer(kind=i4_kind), public:: cpks_size=0,cpks_fitmat_size=0
 integer(kind=i4_kind), public:: cpksalloc(170)=1

  ! 158 prim_int_2cob_ol_dervs(i_grad,k2dr)%m
  ! 55 symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int
  ! 45 prim_int_2cob_ol_dervs+symadapt_int_2cob_ol_dervs

  ! 43 ls_cal_grad overlap_dervs
  ! 42 dervs_matt ll_calc_grads
  ! 36 prim_int_coul_dervs,symadapt_int_coul_dervs
  ! 38 symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int
  ! 35 prim_int_coul_dervs(i_grad,k2dr)%m
  ! 25 cpks_grad_fit
  ! 23 symadapt_int_cpks_coul_grad(i_grad,i_ir)%int
  ! 22 simadd
  ! 21 cont
  ! 19 cpks_inv_chmat
  ! 18 fmat

  ! 17 prim_int_cpks_coul_grad,symadapt_int_cpks_coul_grad
  ! 15 prim_int_cpks_coul_grad
  ! 14 cpks_grad_xc

  ! 3 - cpks_gvec
  ! 4 - cpks_gmat
  ! 5 - linsis dual
  ! 6 - coul_int
  ! 7 - cpks3c
  ! 8 - cpks3c%
  ! 9 - dvdrho
  ! 10 - rho etc
  ! 11 nuc_help_arr
  ! 12 nuc_help
  ! 13 cpks_gvec_temp

! real(kind=r8_kind),public:: avi(2) a test quontity

 logical, public :: comm_dervs_xc=.false.

 integer(kind=i4_kind),public,pointer :: cpks_gradient_index(:)
 integer(kind=i4_kind),public:: i_cpks
 integer(kind=i4_kind),public, parameter:: n_cpks=9
 integer(kind=i4_kind),public, parameter:: n_cpks_split=4

 real(kind=r8_kind),allocatable,public:: cpks_gmat(:,:)
 real(kind=r8_kind),allocatable,public:: dprime_coeff_charge(:)

 logical, public :: totsym_w_ol_grad=.true.

! USE integralpar_cpks_contribs INSTEAD:
! logical, public :: dervs_cpks_contribs=.false.
 logical, public :: cpks_noxc=.false.
 logical, public :: cpks_noco=.false.
 logical, public :: cpks_gmat_calculated=.false.
 logical, public :: cpks_Qxc_calculated=.false.
! logical, public :: calculate_cpks_gvec=.false.
 logical, public :: cpks_energies=.false.

 real(kind=r8_kind), public :: h0_xc
 real(kind=r8_kind), public :: cpksQxc
 real(kind=r8_kind), public :: cpks_ham(2,2)

#ifndef WITH_EXPERIMENTAL
type(arrmat2),pointer,public :: cpks_grad_xc(:)
#endif

                    ! G_kj to calc G_pg_munu integrals
 real(kind=r8_kind),allocatable,public:: cpks_gvec(:,:)
                    ! intermediate to calculate G_aikl*S_kl contrib
 real(kind=r8_kind),allocatable,public:: cpks_fitcoeff_grads(:,:)

 real(kind=r8_kind),public:: dervs_totalsym_fit(9,9)=0.0_r8_kind
 real(kind=r8_kind),public:: hai=0.0_r8_kind
! is used to see dervs_totalsym_fit contib only,
! now commened to be used for larger molecules

 type cpks3c_type
#define no_co_ai
#ifndef no_co_ai
  real(kind=r8_kind),pointer:: co_ai(:,:,:)
#endif
  real(kind=r8_kind),pointer:: eigvec_vir_and_holes(:,:)
  real(kind=r8_kind),pointer:: eigval_vir_and_holes(:)
  real(kind=r8_kind),pointer:: parocc_vir_and_holes(:)
  real(kind=r8_kind),pointer:: parocc_occ(:)
  integer(kind=i4_kind):: n_holes,n_ai,occ_dim,vir_dim
  real(kind=r8_kind),pointer:: V_ai(:,:,:) ! occup x unoccup components of elec.pot
  real(kind=r8_kind),pointer:: V_mn(:,:,:) ! occup x occup components of elec.pot
#if 0
  real(kind=r8_kind),pointer:: solv_ham(:,:)
#endif
 end type cpks3c_type


 type(cpks3c_type),public, allocatable:: cpks3c(:,:)

 type cpks_matrices
  ! Volodya, please correct the comments:
  integer(kind=i4_kind):: occ_dim
  real(kind=r8_kind)::Qxc
  real(kind=r8_kind),pointer:: s1(:,:)     ! (ni,ni), 1st dervs of overlap matrix, occ-occ block
                                           ! are used in calc_cpks_gvec for coul part
                                           ! and to calc nuc_impgrarho (Qai) and nuc_imp_dervs
                                           ! in XC part

  real(kind=r8_kind),pointer:: h1(:,:)     ! (ni,ni), 1st dervs of the Fock matrix, occ-occ block
                                           ! this is an EXPLICIT derivative of the Fock matrix?!?
!  real(kind=r8_kind),pointer:: h1c(:,:)    !?????????????????AS

  real(kind=r8_kind),pointer:: h1ai(:,:)   ! dervs of one particle hamiltonian+coulomb func derivatives
                                           ! (ni,na), 1st dervis of the Fock matrix, occ-vir block
                                           ! this is an EXPLICIT derivative of the Fock matrix?!?

  real(kind=r8_kind),pointer:: s1ai(:,:)   ! (ni,ni), 1st dervs of overlap matrix, occ-vir block

  real(kind=r8_kind),pointer:: Qai(:,:)    ! virtual occupied Q factors
  real(kind=r8_kind),pointer:: B(:,:,:)
  real(kind=r8_kind),pointer:: HBH(:,:)    ! hole occupation modified B matrix
  real(kind=r8_kind),pointer:: ABi(:,:)
  real(kind=r8_kind),pointer:: P1(:,:)  !dervs of densmat
  integer(kind=i4_kind)     :: n_fin_cyc = 0 ! PLEASE comment what that is!

  ! storage for CPKS-specific linear equation matrix:
  real(kind=r8_kind)        :: m(0:n_cpks,0:n_cpks)
  ! for calculations see sum_up_gradient in  integral_calc_quad_2cob3c

  ![[=============================================================
  ! An ugly solution to save 1st derivatives of
  ! the relativisticly transformed Hamiltonian
  !              Trel+Vrel
  ! till after CPKS:
  real(kind=r8_kind),pointer:: hr1(:,:) => NULL()
                             ! hr1(n,n), 1st dervs of rel. ham. H=T+V
                             ! full matrix in the symmetrized basis
  ! NOTE: Since these matrices are spin-independent, only the first
  ! s=1 component of the cpks(igr,irr,s) array will be used.
  !]]==============================================================
 end type cpks_matrices

 public cpks_matrices,cpksalloc_stat,cpks3c_type,print_cpksalloc
 type (cpks_matrices), public, allocatable :: cpks(:,:,:) ! (n_gr,n_irr,n_spin)
                       ! given for each grad proj and irr
  !== Interrupt end of public interface of module ====================
 contains

   subroutine cpksalloc_stat()
     ! presently  called from gradient_data_module  to be  called from
     ! main-slave as there are allocations on slaves also
     use error_module, only: MyID
     implicit none

     integer(i4_kind) :: i

     do i = 1, size(cpksalloc)
        if(cpksalloc(i).eq.0) print *, MyID // 'WARNING: cpksalloc(', i,') == 0'
     enddo
   end subroutine cpksalloc_stat

   subroutine print_cpksalloc()
     !
     ! FIXME: does not work with zero-sized arrays.
     !
#ifdef FPP_DEBUG
     use error_module, only: MyID
#endif
     implicit none

#ifdef FPP_DEBUG
     print *, MyID // 'cpksalloc:', associated(cpks(1,1,1)%s1),associated(cpks(1,1,1)%h1), &
          associated(cpks(1,1,1)%s1ai),associated(cpks(1,1,1)%h1ai), &
          associated(cpks(1,1,1)%Qai),associated(cpks(1,1,1)%B), &
          associated(cpks(1,1,1)%HBH),associated(cpks(1,1,1)%ABi), &
          allocated(cpks_gmat)
#endif
   end subroutine print_cpksalloc

 end module cpksdervs_matrices

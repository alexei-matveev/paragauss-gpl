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
subroutine ss_calculate_grads(i_ea1,i_ea2,imode)
  !-------------------------------------------------------------------
  !
  !  Purpose: calculate integrals for gradients for 
  !           orbitals with both angular momenta = 0.
  !           Symmetry adaption for fitfunctions
  !           is included.
  !           This routines calculates gradients for
  !           a given pair of uniques AND a given
  !           pair of atomic centers
  !
  !  Subroutine called by: ...
  !
  !  References: ...
  !
  !  Author: FN
  !  Date: ...
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Author: AH
  ! Date:   4/99
  ! Description: gradients for pseudopotentials has 
  !              been added
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !------------ Modules used -----------------------------------------
#include "def.h"
  use type_module ! type specification parameters
  use unique_atom_module, noname=>pseudopot_present
  use gamma_module
  use operations_module, only: operations_core_density
  use type_module
  use datatype
  use solid_harmonics_module, only : solid_harmonics_calc
  use solhrules_module, only: solhrules_differential,&
       diff_rule
  use fitcontract_module, only: fitcontract
  use integralpar_module
  use iounitadmin_module
  use gradient_data_module
  use integralpar_module
  use options_module, only: options_integral_expmax, options_split_gradients, &
       options_xcmode, xcmode_model_density, &
       options_spin_restricted
  use int_data_2cob3c_module, only: &
       ua1,ua2,&    ! indices of unique atoms of quadrupel
       ua1_basis,ua2_basis, & ! => uai%l_ob(0)
       center1,center2,& ! coordinates of centers
       quadrupel,&
       prim_int_2cob_ol_grad,&
       prim_int_3cob_grad,&
       prim_int_3cob_nuc_grad,&
       prim_int_3cob_epe,&
       prim_int_3cob_coul_grad,&
       prim_int_2cob_ks_grad,&
       prim_int_2cob_kin_grad,&
       prim_int_2cob_nuc_grad,&
       prim_int_2cob_pvsp_grad, &
       prim_int_3cob_solv_grad, &
       prim_int_3cob_solv_grad_pc, &
       prim_int_2cob_pseudo_grad, &
       prim_int_coul_dervs, &
       prim_int_2cob_ol_dervs 
  use pointcharge_module
#ifdef WITH_EPE
  use ewaldpc_module
#endif

  use solv_cavity_module, only: to_calc_grads, with_pc, fixed_pc
  use elec_static_field_module, only: totsym_field_length,surf_points_grad_index
  use calc3c_switches, only: old_3c_co,old_3c_fc,old_solv_grad,integralpar_dervs,new_3c_co_grad
  use shgi_cntrl, only: IPSEU


  implicit none


  !== Interrupt end of public interface of module ====================

  !------------ Declaration of formal parameters ---------------------
  integer(kind=i4_kind),intent(in)  :: i_ea1,i_ea2
  integer(kind=i8_kind),intent(in)  :: imode ! for control

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !..............................................................................
  ! << OUTPUT ARRAYS >>
  ! ===================
  ! prim_int_2cob_ol_grad  (  1:N) : dsym/dRa < xi_i | 1      | xi_j >
  ! prim_int_2cob_ol_grad  (N+1:M) : dsym/dRb < xi_i | 1      | xi_j >
  ! << relativistic calculation >>
  ! prim_int_2cob_kin_grad (  1:N) : dsym/dRa < xi_i | T      | xi_j >
  ! prim_int_2cob_kin_grad (N+1:M) : dsym/dRb < xi_i | T      | xi_j >
  ! prim_int_2cob_nuc_grad ( ia: ) : dsym/dRa < xi_i | V_nuc  | xi_j >
  ! prim_int_2cob_nuc_grad ( ib: ) : dsym/dRb < xi_i | V_nuc  | xi_j >
  ! prim_int_2cob_nuc_grad ( ic: ) : dsym/dRc < xi_i | V_nuc  | xi_j >
  ! prim_int_2cob_pvsp_grad( ia: ) : dsym/dRa < xi_i | V_pvsp | xi_j >
  ! prim_int_2cob_pvsp_grad( ib: ) : dsym/dRb < xi_i | V_pvsp | xi_j >
  ! prim_int_2cob_pvsp_grad( ic: ) : dsym/dRc < xi_i | V_pvsp | xi_j >
  ! prim_int_3cob_grad     ( ia: ) : dsym/dRa < xi_i | V_H    | xi_j >
  ! prim_int_3cob_grad     ( ib: ) : dsym/dRb < xi_i | V_H    | xi_j >
  ! prim_int_3cob_grad     ( ic: ) : dsym/dRc < xi_i | V_H    | xi_j >
  ! << non-relativistic calculation with total gradients >>
  ! prim_int_3cob_grad     ( ia: ) : dsym/dRa < xi_i | T + V_nuc + V_H | xi_j >
  ! prim_int_3cob_grad     ( ib: ) : dsym/dRb < xi_i | T + V_nuc + V_H | xi_j >
  ! prim_int_3cob_grad     ( ic: ) : dsym/dRc < xi_i |     V_nuc + V_H | xi_j >
  ! << non-relativistic calculation with split gradients >>
  ! prim_int_2cob_ks_grad  (  1:N) : dsym/dRa < xi_i | T + V_nuc + V_H | xi_j >
  ! prim_int_2cob_ks_grad  (N+1:M) : dsym/dRb < xi_i | T + V_nuc + V_H | xi_j >
  ! prim_int_3cob_nuc_grad ( ic: ) : dsym/dRc < xi_i |     V_nuc       | xi_j >
  ! prim_int_3cob_grad     ( ic: ) : dsym/dRc < xi_i |             V_H | xi_j >
  !
  ! Model_Density_Approach
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! << relativistic calculation >>
  ! prim_int_3cob_grad     ( toff+ia: ) : dsym/dRa <xi_i|        V_H+V_X,t|xi_j>
  ! prim_int_3cob_grad     ( toff+ib: ) : dsym/dRb <xi_i|        V_H+V_X,t|xi_j>
  ! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|        V_H+V_X,t|xi_j>
  ! << non-relativistic calculation with total gradients >>
  ! prim_int_3cob_grad     ( toff+ia: ) : dsym/dRa <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
  ! prim_int_3cob_grad     ( toff+ib: ) : dsym/dRb <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
  ! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|  V_nuc+V_H+V_X,t|xi_j>
  ! << non-relativistic calculation with split gradients >>
  ! prim_int_2cob_ks_grad  (toff+  1:N) : dsym/dRa <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
  ! prim_int_2cob_ks_grad  (toff+N+1:M) : dsym/dRb <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
  ! prim_int_3cob_nuc_grad (      ic: ) : dsym/dRc <xi_i|  V_nuc          |xi_j>
  ! prim_int_3cob_coul_grad(      ic: ) : dsym/dRc <xi_i|        V_H      |xi_j>
  ! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|            V_X,t|xi_j>
  !
  ! << WORKING ARRAYS >>
  ! ====================
  ! grad_overlap_a(:,1:3       ) :    d/dRa < xi_i | 1      | xi_j >
  ! grad_kinetic  (:,1:3       ) :    d/dRa < xi_i | T      | xi_j >
  ! grad_nuclear_a(:,1:3       ) :    d/dRa < xi_i | V_nuc  | xi_j >
  ! grad_nuclear_b(:,1:3       ) :    d/dRb < xi_i | V_nuc  | xi_j >
  ! grad_nuclear_c(:,1:3       ) :    d/dRc < xi_i | V_nuc  | xi_j >
  ! nuc_grad      (:,1:grad_dim) : dsym/dRc < xi_i | V_nuc  | xi_j >
  ! grad_pvsp_a   (:,1:3       ) :    d/dRa < xi_i | V_pvsp | xi_j >
  ! grad_pvsp_b   (:,1:3       ) :    d/dRb < xi_i | V_pvsp | xi_j >
  ! help_arr_c    (:,1:3       ) :    d/dRc < xi_i | V_pvsp | xi_j >
  ! pvsp_grad     (:,1:grad_dim) : dsym/dRc < xi_i | V_pvsp | xi_j >
  ! coul_int_a    (  1:3       ) :    d/dRa [ xi_i | f_k    | xi_j ]
  ! coul_int_b    (  1:3       ) :    d/dRb [ xi_i | f_k    | xi_j ]
  ! help_arr      (:,1:3       ) :    d/dRb [ xi_i | f_k    | xi_j ]
  ! coul_int_c    (  1:grad_dim) : dsym/dRb [ xi_i | f_k    | xi_j ]
  !
  ! grad_mat_[xyz]a  (:,:,:,:)   :    d/dRa < xi_i | ...    | xi_j > 
  ! grad_mat_[xyz]b  (:,:,:,:)   :    d/dRb < xi_i | ...    | xi_j > 
  ! grad_spmat_[xyz]a(:,:,:,:)   :    d/dRa < xi_i | ...    | xi_j > 
  ! grad_spmat_[xyz]b(:,:,:,:)   :    d/dRb < xi_i | ...    | xi_j > 
  !..............................................................................

  !------------ Declaration of subroutines ------------------------
  external error_handler
  intrinsic max,maxval
  !------------ Declaration of local constants --------------------
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: one=1.0_r8_kind,&
       two=2.0_r8_kind,&
       three=3.0_r8_kind, &
       four=4.0_r8_kind
! real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
  real(kind=r8_kind),dimension(3,3) :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  !------------ Declaration of local variables --------------------
  integer(kind=i4_kind)          :: naexps,nbexps,ncexps,grad_dim
  real(kind=r8_kind),pointer     :: aexps(:),bexps(:)
  real(kind=r8_kind),pointer     :: cexps(:)
  ! mapping of exponents to one dimension and cutoff of small integrals
  logical,allocatable   :: cutoff(:,:) ! (naexps,nbexps)
  integer(kind=i4_kind) :: num ! metaindex for (naexps,nbexps) > cutoff
  ! help factors
  integer(kind=i4_kind) :: na,nb, &! number of unique atoms
       i_grad,index,counter
  integer(kind=i4_kind)          :: ima, imb, imc, spin_index
  logical                        :: moving_a, moving_b, moving_c, &
       split_gradients, model_density, &
       spin_polarized
  logical                        ::  iam_ppot
  real(kind=r8_kind),allocatable,dimension(:,:):: fact0_arr, &
       fact1_arr,fact2_arr ! (naexps,nbexps)
  real(kind=r8_kind),allocatable,dimension(:)  :: fact0,fact1, &
       fact2,fact3,fac,aexp_arr,bexp_arr ! (num) metaindex for (naexps,nbexps) > cutoff
  real(kind=r8_kind),dimension(3)              :: dist
  ! xa - xb
  ! help factors for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:):: gamma_arg, &
       gamma_help
  real(kind=r8_kind),allocatable,dimension(:)  :: gamma_arg2
  ! help arrays for solid harmonics
  real(kind=r8_kind),allocatable,dimension(:,:)   :: yl_arg
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: yl_arr
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: yl_arr_grad
  ! help arrays for symmetry adaption
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: sym_coef1
  real(kind=r8_kind),allocatable,dimension(:,:,:,:,:) :: sym_coef2
  ! integrals
  real(kind=r8_kind),allocatable,dimension(:) :: overlap1

  type(unique_atom_type), pointer  :: ua_pointer

  real(kind=r8_kind),allocatable,dimension(:,:) :: &
       grad_overlap_a, &
       grad_kinetic, &
       grad_nuclear_a, &
       grad_nuclear_b, &
       grad_nuclear_c, &
       grad_nuc_pseudo_a, &
       grad_nuc_pseudo_b, &
       grad_nuc_pseudo_c, &
       grad_solv_a, &  
       grad_solv_b, &  
       grad_solv_c, &  
       grad_solv_c_help
  real(kind=r8_kind),allocatable,dimension(:,:) :: &
       grad_pseudo_a, &
       grad_pseudo_b, &
       grad_pseudo_c

  logical                         :: pseudopot_present ! same name as in UA module

  type(arrmat4),pointer :: pointer_prim_int(:) ! help pointer to point on primitive integrals

  !  ! 3-center integrals
  !  type(three_center_l_grad)  :: grad_coul_1, grad_coul_2
  type(three_center_l),allocatable  :: coul_int_a(:),coul_int_b(:),coul_int_c(:)

  ! help pointers to unique atoms module data
  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coef(:),rotmat(:,:)
  real(kind=r8_kind),allocatable    :: nuc_grad(:,:),help_arr(:,:)

  real(kind=r8_kind),pointer :: &
       grad_mat_xa(:,:,:,:),grad_mat_ya(:,:,:,:),grad_mat_za(:,:,:,:),&
       grad_mat_xb(:,:,:,:),grad_mat_yb(:,:,:,:),grad_mat_zb(:,:,:,:)
  real(kind=r8_kind),pointer :: &
       grad_spmat_xa(:,:,:,:),grad_spmat_ya(:,:,:,:),grad_spmat_za(:,:,:,:),&
       grad_spmat_xb(:,:,:,:),grad_spmat_yb(:,:,:,:),grad_spmat_zb(:,:,:,:)

  integer(kind=i4_kind) :: alloc_stat,i,ua3,max_order,&
       max_gamma,lmax_ch,k,m1,m2, &
       n_equal_c,lm,n_indep_max,n_indep,i_l,i_sum, &
       n_indep_fcts,n_contributing_fcts,i_ind,i_cont,i_ea3

  integer(kind=i4_kind) :: allocstat(7),i_alo
  integer(kind=i4_kind) :: k2dr
  real(kind=r8_kind),dimension(3)   :: xa,xb,xc
  real(kind=r8_kind)    :: arg,zc
  real(kind=r8_kind) :: zcc ! core charge
  logical :: check_ab,check_bc,check_ac,do_rotation
  !-------------------------------------------------------------------
  !------------ Executable code -----------------------------------

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'ss_grads: PP=',pseudopot_present,' imode=',imode

  allocstat=0
  na=quadrupel%ua1
  nb=quadrupel%ua2
  ima = unique_atoms(na)%moving_atom
  imb = unique_atoms(nb)%moving_atom
  moving_a = ima > 0
  moving_b = imb > 0
  split_gradients = options_split_gradients()
  model_density = options_xcmode() == xcmode_model_density
  spin_polarized = .not. options_spin_restricted()

  ! first get the exponent data ------------------------
  !write(*,'(A22,4I5)') 'ss_calculate_grads:',na,nb,1,1
  naexps = ua1_basis%n_exponents
  nbexps = ua2_basis%n_exponents

  allocate(fact0_arr(nbexps,naexps),fact1_arr(nbexps,naexps), &
           fact2_arr(nbexps,naexps),cutoff(nbexps,naexps), &
                                STAT=allocstat(1))
  ASSERT(allocstat(1).eq.0)
  allocstat(1)=1 ! fact0_arr fact1_arr fact2_arr 
  allocstat(2)=1 ! cutoff

  xa = center1
  xb = center2

  aexps => ua1_basis%exponents(:)
  bexps => ua2_basis%exponents(:)
  arg=sum((xa-xb)**2)

  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  where(fact2_arr*arg>=options_integral_expmax()) ! cutoff: where almost no overlap
     cutoff=.false.            ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)

  all_zero: if(num==0) then ! all integrals are equal zero

     if (integralpar_2cob_ol_grad) then
        do i_grad=1,size(prim_int_2cob_ol_grad)

           prim_int_2cob_ol_grad(i_grad)%m = 0.0_r8_kind

         if(integralpar_dervs) then
          ! FIXME: init prim_int_2cob_ol_dervs elsewhere:
          do k2dr=1,size(prim_int_2cob_ol_grad)
           prim_int_2cob_ol_dervs(i_grad,k2dr)%m = 0.0_r8_kind
          enddo
         endif

        end do
     end if

     if (integralpar_solv_grad) then                        !!!!!!!!!!!
        do i_grad=1,gradient_data_n_spin_gradients          !!!!!!!!!!!
           prim_int_3cob_solv_grad(i_grad)%m = 0.0_r8_kind  !!!!!!!!!!!
        enddo                                               !!!!!!!!!!!
        if (with_pc .and. .not.fixed_pc) then
           do i_grad=1,totsym_field_length          !!!!!!!!!!!
              prim_int_3cob_solv_grad_pc(i_grad)%m = 0.0_r8_kind  !!!!!!!!!!!
           enddo                                               !!!!!!!!!!!
        end if
     endif

     if (integralpar_3cob_grad) then
        do i_grad=1,gradient_data_n_spin_gradients
           prim_int_3cob_grad(i_grad)%m = 0.0_r8_kind
         if(integralpar_dervs) then
           ! FIXME: init prim_int_coul_dervs elsewhere:
           do k2dr=1,gradient_data_n_spin_gradients
            prim_int_coul_dervs(i_grad,k2dr)%m = 0.0_r8_kind
           enddo
         endif
        enddo
     end if

     if (split_gradients) then
        do i_grad=1,size(prim_int_2cob_ks_grad)
           prim_int_2cob_ks_grad(i_grad)%m = 0.0_r8_kind
        end do
        do i_grad=1,gradient_data_n_gradients
           prim_int_3cob_nuc_grad(i_grad)%m = 0.0_r8_kind
        end do
        if (model_density) then
           do i_grad=1,gradient_data_n_gradients
              prim_int_3cob_coul_grad(i_grad)%m = 0.0_r8_kind
           end do
        endif
     endif

     deallocate(fact0_arr,fact1_arr,fact2_arr,cutoff,stat=allocstat(1))
     if (allocstat(1).ne.0) call error_handler &
          ("SS_CALCULATE: deallocation failed")
     allocstat(2)=0 ! cutoff
     return

  end if all_zero

  allocate (fact0(num),fact1(num),fact2(num),&
            fact3(num),fac(num),help_arr(num,3),&
            aexp_arr(num),bexp_arr(num), STAT=allocstat(3))
     ASSERT(allocstat(3).eq.0)
     allocstat(3)=1

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)
  aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  bexp_arr=pack(spread(bexps,2,naexps),cutoff)
  help_arr=zero
  deallocate(fact0_arr,fact1_arr,fact2_arr,STAT=allocstat(1))
  ASSERT(allocstat(1).eq.0)

  allocate(grad_mat_xa(nbexps,naexps,1,1),&
       grad_mat_ya(nbexps,naexps,1,1),&
       grad_mat_za(nbexps,naexps,1,1),&
       grad_mat_xb(nbexps,naexps,1,1),&
       grad_mat_yb(nbexps,naexps,1,1),&
       grad_mat_zb(nbexps,naexps,1,1),stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

  grad_mat_xa=0.0_r8_kind
  grad_mat_ya=0.0_r8_kind
  grad_mat_za=0.0_r8_kind
  grad_mat_xb=0.0_r8_kind
  grad_mat_yb=0.0_r8_kind
  grad_mat_zb=0.0_r8_kind
  
  if (model_density .and. spin_polarized) then
     allocate(grad_spmat_xa(nbexps,naexps,1,1),&
          grad_spmat_ya(nbexps,naexps,1,1),&
          grad_spmat_za(nbexps,naexps,1,1),&
          grad_spmat_xb(nbexps,naexps,1,1),&
          grad_spmat_yb(nbexps,naexps,1,1),&
          grad_spmat_zb(nbexps,naexps,1,1),stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation  failed grad_spmat")
     grad_spmat_xa=0.0_r8_kind
     grad_spmat_ya=0.0_r8_kind
     grad_spmat_za=0.0_r8_kind
     grad_spmat_xb=0.0_r8_kind
     grad_spmat_yb=0.0_r8_kind
     grad_spmat_zb=0.0_r8_kind
  endif

  ! now prepare arguments for incomplete Gamma-Fct.
  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  ! allocation: num*max_order + num*3 + num
  lmax_ch = maxval(unique_atoms(:)%lmax_ch)
  max_order = max(lmax_ch+2,3)
  max_gamma=2
  if(integralpar_relativistic) then
     max_gamma=4
     max_order=max(max_order,4)
  end if
  allocate (gamma_help(num,max_order), gamma_arg(num,3), &
            gamma_arg2(num), &
                                           STAT=allocstat(4))
  if(allocstat(4).ne.0) call error_handler &
       ("SS_CALCULATE_GRAD : allocation gamma_help failed ")
     allocstat(4)=1

  do i=1,3
     gamma_arg(:,i) = pack( spread(aexps,1,nbexps)*xa(i) + &
          spread(bexps,2,naexps)*xb(i),cutoff )/fact0
  enddo

  ! ---------------------------------------------------

  m1=1
  m2=1

  ! first calculating 2-center integrals----------------
  ! PROCESS <xi_i|xi_j> and <xi_i|T|xi_j> 
  call calculate_ol_and_kin

  notsolv0: if( .not. integralpar_solv_grad) then    
  allocate( grad_nuclear_a(num,3),STAT=alloc_stat )
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE_GRAD: allocation grad_nuclear_a failed")
  grad_nuclear_a = zero
  allocate( grad_nuclear_b(num,3),STAT=alloc_stat )
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE_GRAD: allocation grad_nuclear_b failed")
  grad_nuclear_b = zero

  allocate( grad_nuclear_c(num,3),STAT=allocstat(6) )
  if (allocstat(6).ne.0) call error_handler &
       ("SS_CALCULATE_GRAD: allocation grad_nuclear_c failed")
      allocstat(6)=1
  grad_nuclear_c = zero

  if (pseudopot_present) then
     allocate( grad_pseudo_a(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_a failed")
     grad_pseudo_a = zero
     allocate( grad_pseudo_b(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_b failed")
     grad_pseudo_b = zero
     allocate( grad_pseudo_c(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_c failed")
     grad_pseudo_c = zero
  endif
  if (pseudopot_present.and.integralpar_relativistic) then
     allocate( grad_nuc_pseudo_a(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_a failed")
     grad_nuc_pseudo_a = zero
     allocate( grad_nuc_pseudo_b(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_b failed")
     grad_nuc_pseudo_b = zero
     allocate( grad_nuc_pseudo_c(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation grad_nuc_pseudo_c failed")
     grad_nuc_pseudo_c = zero
  endif
  endif  notsolv0
  fac = sqrt(fact0/pi)
!  prim_int_3cob_grad(2)%m(:,:,1,1)=unpack(overlap1,cutoff,zero) !no nuc

  notsolv1 : if( .not. integralpar_solv_grad) then      
!  prim_int_3cob_grad(2)%m(:,:,1,1)=unpack(overlap1,cutoff,zero) !no nuc
  ! loop over atomic centers
  ! PROCESS <xi_i|V_nuc|xi_j>, <xi_i|V_pvsp|xi_j>, and <xi_i|V_H|xi_j>
!  prim_int_3cob_grad(2)%m(:,:,1,1)=unpack(overlap1,cutoff,zero) !no nuc
  unique_3: do ua3 = 1,N_unique_atoms+n_timps
     if(ua3<=n_unique_atoms) then
        ua_pointer=>unique_atoms(ua3)
        zc = ua_pointer%Z
        zcc= ua_pointer%ZC
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'ss_grads: ua=',ua3,', zero its charge!'
        zcc = zero
        zc  = zero
        
        n_equal_c = unique_atoms(ua3)%N_equal_atoms
        imc = unique_atoms(ua3)%moving_atom
        moving_c = imc > 0
        if (moving_c) then
           grad_dim = gradient_index(imc+1) - gradient_index(imc)
        else
           grad_dim = 0
        endif

! ** allocate calculate and deallocate nuc_grad intermediate
        allocate(nuc_grad(num,grad_dim),stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)
        nuc_grad=zero

        
        ! PROCESS <xi_i|Zc/|r-Rc||xi_j> and <xi_i|V_pvsp[c]|xi_j>
        iam_ppot=zcc/=0.0_r8_kind &
           .and.(.not.operations_core_density).and. &
               pseudopot_present
        call calculate_nuc() ! requires iam_ppot
     else ! it is TIMP
        ua_pointer=> unique_timps(ua3-n_unique_atoms)
        zcc= ua_pointer%ZC 
!!$        moving_c=.false. !!!! not false for EPE calculations
        imc = ua_pointer%moving_atom !!!+N_moving_unique_atoms
        moving_c = imc > 0
        if (moving_c) then
           grad_dim = gradient_index(N_moving_unique_atoms+imc+1) &
                - gradient_index(N_moving_unique_atoms+imc)
           allocate(nuc_grad(num,grad_dim),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("SS_CACLULATE : allocation failed nuc_grad")
           nuc_grad=zero
        else
           grad_dim = 0
        endif
        n_equal_c=ua_pointer%n_equal_atoms
     end if !atom/timp
     ! PROCESS  pseudopotentials
     iam_ppot=zcc/=0.0_r8_kind.and. &
            (.not.operations_core_density).and. &
          pseudopot_present
     if (iam_ppot) then
        ABORT('not supported')
     endif
        ! result is in grad_nuc_pseudo_a, grad_nuc_pseudo_b and nuc_grad for c
   
     if(ua3<=n_unique_atoms) then
        ! ADD dsym/dRc <xi_i|V_nuc[c]|xi_j>, dsym/dRc <xi_i|V_pvsp[c]|xi_j>
        if (moving_c) index = gradient_index(imc) - 1 ! may be shifted up ???
        do i_grad=1,grad_dim ! only if moving_c
           if(.not.integralpar_relativistic) then
              if (split_gradients) then
                 prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,1,1) = &
                      prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,1,1) - &
                      unpack(nuc_grad(:,i_grad),cutoff,zero)
              else
!              if(.not.cpks_energies) then
                 prim_int_3cob_grad(index+i_grad)%m(:,:,1,1)=&        ! no nuc grad if commented
                      prim_int_3cob_grad(index+i_grad)%m(:,:,1,1) &
                            -unpack(nuc_grad(:,i_grad),cutoff,zero)
!               endif
              endif

           else ! integralpar_relativistic and regular PP
              ppot: if(iam_ppot) then

                 prim_int_2cob_pseudo_grad(index+i_grad)%m(:,:,1,1)=&
                   prim_int_2cob_pseudo_grad(index+i_grad)%m(:,:,1,1) &
                           -unpack(nuc_grad(:,i_grad),cutoff,zero)
              else

                 prim_int_2cob_nuc_grad(index+i_grad)%m(:,:,1,1)=&
                      prim_int_2cob_nuc_grad(index+i_grad)%m&
                      (:,:,1,1)-unpack(nuc_grad(:,i_grad),cutoff,zero)
              end if ppot
           end if
        end do

        deallocate(nuc_grad,stat=alloc_stat)        
        if (alloc_stat/=0) call error_handler &
             ("SS_CACLULATE : deallocation failed nuc_grad")

        ! now fitfunctions

        ! Orbital gradients:
        lmax_ch = int(unique_atoms(ua3)%lmax_ch,kind=i4_kind)
        max_order = max(3,lmax_ch+2)

      if(.not.new_3c_co_grad) then

        allocate(coul_int_a(3),coul_int_b(3),coul_int_c(grad_dim),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_a")

        allocate(coul_int_a(1)%l(-1:lmax_ch),coul_int_a(2)%l(-1:lmax_ch),&
             coul_int_a(3)%l(-1:lmax_ch),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_a (1) ")

        allocate(coul_int_b(1)%l(-1:lmax_ch),&
             coul_int_b(2)%l(-1:lmax_ch),&
             coul_int_b(3)%l(-1:lmax_ch),&
             STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_b (1) ")



        ncexps = unique_atoms(ua3)%r2_ch%n_exponents
        allocate(coul_int_a(1)%l(-1)%m(num,ncexps,1,1,1),&
             coul_int_a(2)%l(-1)%m(num,ncexps,1,1,1),&
             coul_int_a(3)%l(-1)%m(num,ncexps,1,1,1),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_a ")
        coul_int_a(1)%l(-1)%m = zero
        coul_int_a(2)%l(-1)%m = zero
        coul_int_a(3)%l(-1)%m = zero
        allocate(coul_int_b(1)%l(-1)%m(num,ncexps,1,1,1),&
             coul_int_b(2)%l(-1)%m(num,ncexps,1,1,1),&
             coul_int_b(3)%l(-1)%m(num,ncexps,1,1,1),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_b ")
        coul_int_b(1)%l(-1)%m = zero
        coul_int_b(2)%l(-1)%m = zero
        coul_int_b(3)%l(-1)%m = zero

        ncexps = unique_atoms(ua3)%l_ch(0)%n_exponents
        allocate(coul_int_a(1)%l(0)%m(num,ncexps,1,1,1),&
             coul_int_a(2)%l(0)%m(num,ncexps,1,1,1),&
             coul_int_a(3)%l(0)%m(num,ncexps,1,1,1),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_a ")     
        coul_int_a(1)%l(0)%m = zero
        coul_int_a(2)%l(0)%m = zero
        coul_int_a(3)%l(0)%m = zero
        allocate(coul_int_b(1)%l(0)%m(num,ncexps,1,1,1),&
             coul_int_b(2)%l(0)%m(num,ncexps,1,1,1),&
             coul_int_b(3)%l(0)%m(num,ncexps,1,1,1),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation coul_int_b ")     

        coul_int_b(1)%l(0)%m = zero
        coul_int_b(2)%l(0)%m = zero
        coul_int_b(3)%l(0)%m = zero

        do i_grad=1,grad_dim ! only if moving_c
           allocate(coul_int_c(i_grad)%l(-1:lmax_ch),&
                stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("SS_CACLULATE : allocation failed &
                &coul_int_c(i_grads)%l(-1:lmax_ch")
           ncexps = unique_atoms(ua3)%r2_ch%n_exponents
           allocate(coul_int_c(i_grad)%l(-1)%m&
                (num,ncexps,1,1,1),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("SS_CACLULATE : allocation failed &
                &coul_int_c(i_grad)%l(-1)%m")
           ncexps = unique_atoms(ua3)%l_ch(0)%n_exponents
           allocate(coul_int_c(i_grad)%l(0)%m&
                (num,ncexps,1,1,1),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("SS_CACLULATE : allocation failed &
                &coul_int_c(i_grad)%l(0)%m")
           coul_int_c(i_grad)%l(0)%m=0.0_r8_kind
           coul_int_c(i_grad)%l(-1)%m=0.0_r8_kind
        end do
      endif

        if (lmax_ch.gt.0) then
           allocate(yl_arr(num,(lmax_ch+1)**2,n_equal_c),STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: allocation yl_arr failed")
           yl_arr = zero
           allocate(yl_arr_grad(num,n_equal_c,(lmax_ch+1)**2,3)&
                ,STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: allocation yl_arr_grad failed")
           yl_arr_grad = zero
           n_indep_max = maxval(unique_atoms(ua3)%symadapt_partner&
                (1,:)%n_independent_fcts)
           allocate(sym_coef1(num,n_equal_c,n_indep_max,lmax_ch),&
                STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: allocation sym_coef1 failed")
           sym_coef1 = zero
           allocate(sym_coef2(num,n_equal_c,n_indep_max,lmax_ch,3),&
                STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: allocation sym_coef2 failed")
           sym_coef2 = zero

          if(.not.new_3c_co_grad) then
           do i_l = 1,lmax_ch
              ncexps = unique_atoms(ua3)%l_ch(i_l)%n_exponents
              n_indep =  &
                   unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
              allocate( coul_int_a(1)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   coul_int_a(2)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   coul_int_a(3)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   STAT=alloc_stat)
              if (alloc_stat.ne.0)  call error_handler &
                   ("SS_CALCULATE_GRAD: allocation coul_int_a, ang. mom. failed")
              coul_int_a(1)%l(i_l)%m = zero
              coul_int_a(2)%l(i_l)%m = zero
              coul_int_a(3)%l(i_l)%m = zero
              allocate (coul_int_b(1)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   coul_int_b(2)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   coul_int_b(3)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                   STAT=alloc_stat)
              if (alloc_stat.ne.0)  call error_handler &
                   ("SS_CALCULATE_GRAD: allocation coul_int_b, ang. mom. failed")
              coul_int_b(1)%l(i_l)%m = zero
              coul_int_b(2)%l(i_l)%m = zero
              coul_int_b(3)%l(i_l)%m = zero
              do i_grad=1,grad_dim ! only if moving_c
                 allocate (coul_int_c(i_grad)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                      STAT=alloc_stat)
                 if (alloc_stat.ne.0)  call error_handler &
                      ("SS_CALCULATE_GRAD: allocation coul_int_c, ang. mom. failed")
                 coul_int_c(i_grad)%l(i_l)%m = zero
              end do
           enddo
          endif
        endif

        equal_3_coul: do i_ea3 = 1,n_equal_c
           if (moving_c) then
              if(grad_dim==3) then
                 if(sum((unique_atom_grad_info(imc)%m(:,:,i_ea3)-unity_matrix)**2)<&
                      1.0e-7_r8_kind) then
                    do_rotation=.false.
                 else
                    do_rotation=.true.
                    rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
                 endif
              else
                 do_rotation=.true.
                 rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
              end if
           endif
           xc = unique_atoms(ua3)%position(:,i_ea3)
           check_ab =( (quadrupel%ua1.eq.quadrupel%ua2).and. &
                (i_ea1.eq.i_ea2))
           check_ac =( (quadrupel%ua1.eq.ua3).and. &
                (i_ea1.eq.i_ea3))
           check_bc =( (quadrupel%ua2.eq.ua3).and. &
                (i_ea2.eq.i_ea3))


           ! now do a precalculation of solid harmonics
           if (lmax_ch.gt.0) then

              allocate(yl_arg(num,3),STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_grad : allocation yl_arg failed")

              yl_arg(:,1) = gamma_arg(:,1) - xc(1)
              yl_arg(:,2) = gamma_arg(:,2) - xc(2)
              yl_arg(:,3) = gamma_arg(:,3) - xc(3)

              yl_arr(:,:,i_ea3) = solid_harmonics_calc(lmax_ch,yl_arg)

              deallocate(yl_arg,STAT=alloc_stat) ! this is not needed anymore
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_grad : deallocation yl_arg failed")

              do lm=1,(lmax_ch+1)**2
                 do i_sum=1,solhrules_differential(3,lm)%n_summands
                    yl_arr_grad(:,i_ea3,lm,1) = &
                         yl_arr_grad(:,i_ea3,lm,1) + &
                         solhrules_differential(3,lm)%coef(i_sum)* &
                         yl_arr(:,solhrules_differential(3,lm)%lm_sh(i_sum),i_ea3)
                 enddo
                 do i_sum=1,solhrules_differential(4,lm)%n_summands
                    yl_arr_grad(:,i_ea3,lm,2) = &
                         yl_arr_grad(:,i_ea3,lm,2) + &
                         solhrules_differential(4,lm)%coef(i_sum)* &
                         yl_arr(:,solhrules_differential(4,lm)%lm_sh(i_sum),i_ea3)
                 enddo
                 do i_sum=1,solhrules_differential(2,lm)%n_summands
                    yl_arr_grad(:,i_ea3,lm,3) = &
                         yl_arr_grad(:,i_ea3,lm,3) + &
                         solhrules_differential(2,lm)%coef(i_sum)* &
                         yl_arr(:,solhrules_differential(2,lm)%lm_sh(i_sum),i_ea3)
                 enddo
              enddo

           endif


           ! if all three centers are the same, the integrals are zero
           if ( check_ab.and.check_bc ) then
              cycle equal_3_coul
           endif

           ! PROCESS [xi_i|f_k|xi_j]
           ! first evaluate s-type coulomb integrals
           if(old_3c_co) call s_coulomb()
           ! now treating r2-type coloumb integrals
           if(old_3c_co) call r2_coulomb()

        enddo equal_3_coul

        if (lmax_ch.gt.0) then

           call calc_sym_coef()

           ! finally calculate l-type coulomb integrals
           if(old_3c_co) call l_coulomb()

        endif

        ! now let`s call fitcontract

        if(integralpar_3cob_grad) then
           ! coulombic contributions of derivatives with respect
           ! to center a and b
           ! CONTRACT d/dRa [xi_i|f_k|xi_j] and d/dRb [xi_i|f_k|xi_j]
           if (model_density .and. spin_polarized) then
             if(old_3c_fc.and.old_3c_co) then
              if (moving_a) then
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(1), &
                      grad_mat_xa,grad_spmat_xa)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(2), &
                      grad_mat_ya,grad_spmat_ya)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(3), &
                      grad_mat_za,grad_spmat_za)
              endif
              if (moving_b) then
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(1), &
                      grad_mat_xb,grad_spmat_xb)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(2), &
                      grad_mat_yb,grad_spmat_yb)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(3), &
                      grad_mat_zb,grad_spmat_zb)
              endif
             endif
           else ! standard SCF
	if(old_3c_co) then
              if (moving_a) then
!	 	do  i_l=-1,lmax_ch
!	       if(old_3c_co.and.associated(coul_int_a(3)%l(i_l)%m)) &
!		 print*,sum(coul_int_a(3)%l(i_l)%m),i_l,' ga'
!		enddo
		if(old_3c_fc) then
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(1),&
                      grad_mat_xa)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(2),&
                      grad_mat_ya)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_a(3),&
                      grad_mat_za)
	        endif
              endif
              if (moving_b) then
		if(old_3c_fc) then
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(1),&
                      grad_mat_xb)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(2),&
                      grad_mat_yb)
                 call fitcontract('grad',num,ua3,cutoff,coul_int_b(3),&
                      grad_mat_zb)
		endif
              endif
           endif
	endif

           ! CONTRACT dsym/dRc [xi_i|f_k|xi_j]
           do i_grad=1,grad_dim ! only if noving_c
              index = gradient_index(imc) + i_grad - 1
              if (model_density) then
                 if (spin_polarized) then
                    spin_index = index + gradient_data_n_gradients
                    if (split_gradients) then
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,ua3,cutoff, &
                            coul_int_c(i_grad), &
                            prim_int_3cob_coul_grad(index)%m, & ! V_H
                            prim_int_3cob_grad(spin_index)%m, & ! V_X,spin
                            prim_int_3cob_grad(index)%m )       ! V_X,tot
                    else ! total gradients only
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,ua3,cutoff, &
                            coul_int_c(i_grad), &
                            prim_int_3cob_grad(index)%m, &      ! V_H+V_x,tot
                            prim_int_3cob_grad(spin_index)%m )  ! V_X,spin
                    endif
                 else ! spin_restricted
                    if (split_gradients) then
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,ua3,cutoff, &
                            coul_int_c(i_grad), &
                            prim_int_3cob_coul_grad(index)%m, & ! V_H
                            mda_xcpot_gradients = &
                            prim_int_3cob_grad(index)%m )       ! V_X,tot
                    else ! total gradients only
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,ua3,cutoff, &
                            coul_int_c(i_grad), &
                            prim_int_3cob_grad(index)%m )       ! V_H+V_X,tot
                    endif
                 end if
              else ! standard SCF
	if(old_3c_fc.and.old_3c_co) then
!		do i_l=-1,lmax_ch
!		 if(old_3c_co.and.associated(coul_int_c(i_grad)%l(i_l)%m)) &
!                       print*,'c ',i_l,i_grad,sum(coul_int_c(i_grad)%l(i_l)%m)
!	        enddo
                 call fitcontract('grad',num,ua3,cutoff,coul_int_c(i_grad),&
                      prim_int_3cob_grad(index)%m)
	endif
              endif

      if(.not.new_3c_co_grad) then
              do i_l = -1, lmax_ch
                 deallocate(coul_int_c(i_grad)%l(i_l)%m,&
                      STAT=alloc_stat)
                 if(alloc_stat/=0) call error_handler &
                      ("SS_CALCULATE : deallocation coul_int_c&
                      &%l%m failed")
              enddo

              deallocate (coul_int_c(i_grad)%l,STAT=alloc_stat)

              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation coul_int_c&
                   &(i_grad)%l failed")
       endif
           end do

      if(.not.new_3c_co_grad) then
           deallocate (coul_int_c,STAT=alloc_stat)

           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : deallocation coul_int_c failed")
       endif

      if(.not.new_3c_co_grad) then
           do i_l=-1,lmax_ch
              deallocate(coul_int_a(1)%l(i_l)%m,&
                   coul_int_a(2)%l(i_l)%m,&
                   coul_int_a(3)%l(i_l)%m,&
                   STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation coul_int_a (1) failed")
              deallocate(coul_int_b(1)%l(i_l)%m,&
                   coul_int_b(2)%l(i_l)%m,&
                   coul_int_b(3)%l(i_l)%m,&
                   STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation coul_int_b (1) failed")

           enddo

           deallocate( coul_int_a(1)%l,&
                coul_int_a(2)%l,&
                coul_int_a(3)%l,&
                STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: deallocation coul_int_a (2) failed")
           deallocate( coul_int_b(1)%l,&
                coul_int_b(2)%l,&
                coul_int_b(3)%l,&
                STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: deallocation coul_int_b (2) failed")

           deallocate( coul_int_a,&
                coul_int_b,&
                STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: deallocation coul_int_a failed")
       endif

           if (lmax_ch.gt.0) then
              deallocate(sym_coef2,STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation symcoef2 failed")
              deallocate(sym_coef1,STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation symcoef1 failed")
              deallocate(yl_arr_grad,STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation yl_arr_grad failed")
              deallocate(yl_arr,STAT=alloc_stat)
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE_GRAD: deallocation yl_arr failed")
           endif
        end if

     else !treat tipms
        if(moving_c) then
           index = gradient_index(imc+n_moving_unique_atoms) - 1
           do i_grad=1,grad_dim ! only if moving_c
              prim_int_3cob_grad(index+i_grad)%m(:,:,1,1)=&
                   prim_int_3cob_grad(index+i_grad)%m&
                   (:,:,1,1)-unpack(nuc_grad(:,i_grad),cutoff,zero)
           end do

           deallocate(nuc_grad,STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE_GRAD: deallocation nuc_grad for timps failed")
        end if
     
     end if
  enddo unique_3
  endif notsolv1


  ! Add contribution of point charges to gradients of nuc and pvsp
  ! PROCESS <xi_i|Zp/|r-Rp||xi_j> and <xi_i|V_pvsp[p]|xi_j>

  if (pointcharge_N+n_timps.gt. 0) call calculate_nuc_pc

  if (integralpar_solv_grad.and.old_solv_grad ) call calculate_solv()  !!!!!!!!!!!!!

  map: if (integralpar_3cob_grad) then
     if(.not.integralpar_relativistic) then 
        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
        if (moving_a) then
           grad_mat_xa(:,:,1,1)=grad_mat_xa(:,:,1,1)-&
                unpack(grad_nuclear_a(:,1),cutoff,zero)
           grad_mat_ya(:,:,1,1)=grad_mat_ya(:,:,1,1)-&
                unpack(grad_nuclear_a(:,2),cutoff,zero)
           grad_mat_za(:,:,1,1)=grad_mat_za(:,:,1,1)-&
                unpack(grad_nuclear_a(:,3),cutoff,zero)
        endif
        if (moving_b) then
           grad_mat_xb(:,:,1,1)=grad_mat_xb(:,:,1,1)-&
                unpack(grad_nuclear_b(:,1),cutoff,zero)
           grad_mat_yb(:,:,1,1)=grad_mat_yb(:,:,1,1)-&
                unpack(grad_nuclear_b(:,2),cutoff,zero)
           grad_mat_zb(:,:,1,1)=grad_mat_zb(:,:,1,1)-&
                unpack(grad_nuclear_b(:,3),cutoff,zero)
        endif

        deallocate(grad_nuclear_a,grad_nuclear_b,stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)

        if (pseudopot_present)then ! for .not.integralpar_relativistic
           if (moving_a) then
              grad_mat_xa(:,:,1,1)=grad_mat_xa(:,:,1,1)+&
                   unpack(grad_pseudo_a(:,1),cutoff,zero)
              grad_mat_ya(:,:,1,1)=grad_mat_ya(:,:,1,1)+&
                   unpack(grad_pseudo_a(:,2),cutoff,zero)
              grad_mat_za(:,:,1,1)=grad_mat_za(:,:,1,1)+&
                   unpack(grad_pseudo_a(:,3),cutoff,zero)
           endif
           if (moving_b) then
              grad_mat_xb(:,:,1,1)=grad_mat_xb(:,:,1,1)+&
                   unpack(grad_pseudo_b(:,1),cutoff,zero)
              grad_mat_yb(:,:,1,1)=grad_mat_yb(:,:,1,1)+&
                   unpack(grad_pseudo_b(:,2),cutoff,zero)
              grad_mat_zb(:,:,1,1)=grad_mat_zb(:,:,1,1)+&
                   unpack(grad_pseudo_b(:,3),cutoff,zero)
           endif
        endif

     end if ! .not. integralpar_relativistic

     ! LOAD prim_int_3cob_grad(...,ia:) and prim_int_3cob_grad(...,ib: ) or
     !      prim_int_2cob_ks_grad( 1:N) and prim_int_2cob_ks_grad(N+1:M)
     if (split_gradients) then
        pointer_prim_int=>prim_int_2cob_ks_grad
        call add_grads(compact_storage=.TRUE.)
        if (model_density .and. spin_polarized) then
           spin_index = size(prim_int_2cob_ks_grad)/2 + 1
           pointer_prim_int=>prim_int_2cob_ks_grad(spin_index:)
           call grad_spmat_to_grad_mat()
           call add_grads(compact_storage=.TRUE.)
        end if

     else

        pointer_prim_int=>prim_int_3cob_grad
        call add_grads() 

        if (model_density .and. spin_polarized) then
           spin_index = gradient_data_n_gradients + 1
           pointer_prim_int=>prim_int_3cob_grad(spin_index:)
           call grad_spmat_to_grad_mat()
           call add_grads()
        endif
     endif

     relmap: if(integralpar_relativistic) then
        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
       if (moving_a) then

           grad_mat_xa(:,:,1,1)=-unpack(grad_nuclear_a(:,1),cutoff,zero)
           grad_mat_ya(:,:,1,1)=-unpack(grad_nuclear_a(:,2),cutoff,zero)
           grad_mat_za(:,:,1,1)=-unpack(grad_nuclear_a(:,3),cutoff,zero)

       endif
       if (moving_b) then

           grad_mat_xb(:,:,1,1)=-unpack(grad_nuclear_b(:,1),cutoff,zero)
           grad_mat_yb(:,:,1,1)=-unpack(grad_nuclear_b(:,2),cutoff,zero)
           grad_mat_zb(:,:,1,1)=-unpack(grad_nuclear_b(:,3),cutoff,zero)

      endif
        deallocate(grad_nuclear_a,grad_nuclear_b,stat=alloc_stat)
        if(alloc_stat/=0) call error_handler &
             ("SS_CALCULATE: deallocation grad_nuclear_a failed")
        ! LOAD prim_int_2cob_nuc_grad(...,ia:), 
        !      prim_int_2cob_nuc_grad(...,ib:)
        pointer_prim_int=>prim_int_2cob_nuc_grad
               call add_grads()

        if(pseudopot_present) then
          if (moving_a) then  
           grad_mat_xa(:,:,1,1)=unpack(grad_pseudo_a(:,1)& 
                                  -grad_nuc_pseudo_a(:,1),cutoff,zero)
	   grad_mat_ya(:,:,1,1)=unpack(grad_pseudo_a(:,2) &
                                  -grad_nuc_pseudo_a(:,2),cutoff,zero)
	   grad_mat_za(:,:,1,1)=unpack(grad_pseudo_a(:,3)&
                                  -grad_nuc_pseudo_a(:,3),cutoff,zero)  
          end if
       if (moving_b) then
           grad_mat_xb(:,:,1,1)=unpack(grad_pseudo_b(:,1)&
                                   -grad_nuc_pseudo_b(:,1),cutoff,zero)
	   grad_mat_yb(:,:,1,1)=unpack(grad_pseudo_b(:,2)&
                                   -grad_nuc_pseudo_b(:,2),cutoff,zero)
	   grad_mat_zb(:,:,1,1)=unpack(grad_pseudo_b(:,3)&
                                   -grad_nuc_pseudo_b(:,3),cutoff,zero)
       end if
      
          pointer_prim_int=>prim_int_2cob_pseudo_grad
       

          call add_grads()
       end if ! pseudopot_present
       
     end if relmap


     deallocate(&
          grad_mat_xa,grad_mat_xb,&
          grad_mat_ya,grad_mat_yb,&
          grad_mat_za,grad_mat_zb,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

  endif map

  if(.not.integralpar_solv_grad) &
  deallocate ( grad_nuclear_c,STAT=allocstat(6)) ! overlap1 grad_overlap_a
       if (allocstat(6).ne.0) call error_handler('ss grad deallocate 6 failed')

!!$  notsolv2: if(.not.(integralpar_solv_grad.and.old_solv_grad)) then
  notsolv2: if(.not.integralpar_solv_grad) then
  deallocate (aexp_arr,bexp_arr,fact0,fact1,fact2,fact3,fac,help_arr, &
                                                     STAT=allocstat(3))
  if (allocstat(3).ne.0) call error_handler &       
       ("SS_GRAD: deallocate (3) failed") 

  deallocate (cutoff,stat=allocstat(2))
  if (allocstat(2).ne.0) call error_handler('ss_grad deallocate 2 failed')
 
  deallocate (gamma_help,gamma_arg,gamma_arg2, &
                              STAT=allocstat(4)) ! gamma_help gamma_arg gamma_arg2
       if (allocstat(4).ne.0) call error_handler('ss grad deallocate 4 failed')
  deallocate ( overlap1, grad_overlap_a, STAT=allocstat(5)) ! overlap1 grad_overlap_a
       if (allocstat(5).ne.0) call error_handler('ss grad deallocate 5 failed')

  if (pseudopot_present) then
     deallocate(grad_pseudo_c,&
          grad_pseudo_a,grad_pseudo_b, stat=alloc_stat)
!       print*,'nuclear ss grad',sum(prim_int_3cob_grad(2)%m(:,:,1,1))
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: deallocation of grad_pseudo_c failed")

     if(integralpar_relativistic) then
     deallocate(grad_nuc_pseudo_a,grad_nuc_pseudo_b,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ("SS_CALCULATE_GRAD: deallocation grad_pseudo_ab failed")
     deallocate(grad_nuc_pseudo_c, stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE_GRAD: allocation of grad_nuc_pseudo_c failed")
     endif
  end if
  endif notsolv2

   do i_alo=1,size(allocstat)
    if(allocstat(i_alo) .ne.0) print*,i_alo ,'ss grad allocstat ne 0'
   enddo

contains

  subroutine grad_spmat_to_grad_mat()
    ! Purpose: deallocate grad_mat[xyz][ab] and link
    !          grad_spmat[xyz][ab] to grad_mat[xyz][ab]

!!! MF Bug:
!   deallocate( grad_mat_xa, grad_mat_ya, grad_mat_za, &
!        grad_mat_xa, grad_mat_ya, grad_mat_za, STAT=alloc_stat )
    deallocate( grad_mat_xa, grad_mat_ya, grad_mat_za, &
         grad_mat_xb, grad_mat_yb, grad_mat_zb, STAT=alloc_stat )
    if (alloc_stat /= 0) call error_handler("SS_CALCULATE&
         &[grad_spmat_to_grad_mat]: deallocation grad_mat_xa failed")

    grad_mat_xa => grad_spmat_xa
    grad_mat_ya => grad_spmat_ya
    grad_mat_za => grad_spmat_za
    grad_mat_xb => grad_spmat_xb
    grad_mat_yb => grad_spmat_yb
    grad_mat_zb => grad_spmat_zb

  end subroutine grad_spmat_to_grad_mat

  subroutine add_grads(compact_storage)
    logical, optional :: compact_storage
    ! Purpose: add gradients to the appropriate data structures

    ! LOAD dsym/dRa < xi_i | ... | xi_j >
    ! Input: grad_mat_[xyz]a(...) , Output: pointer_prim_int(ia:) or
    !                                       pointer_prim_int(1:N)
    if (moving_a) then
       if (present(compact_storage)) then
          index = 1
       else
          index = gradient_index(ima)
       endif
       grad_dim = gradient_index(ima+1) - gradient_index(ima)
       if(grad_dim==3) then
          if(sum((unique_atom_grad_info(ima)%m(:,:,i_ea1)-unity_matrix)**2)<&
               1.0e-7_r8_kind) then
             do_rotation=.false.
          else
             do_rotation=.true.
             rotmat=>unique_atom_grad_info(ima)%m(:,:,i_ea1)
          endif
       else
          do_rotation=.true.
          rotmat=>unique_atom_grad_info(ima)%m(:,:,i_ea1)
       end if

       if(do_rotation) then
          do i_grad=1,grad_dim ! only if moving_a
             pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
                  rotmat(i_grad,1)*grad_mat_xa+&
                  rotmat(i_grad,2)*grad_mat_ya+&
                  rotmat(i_grad,3)*grad_mat_za
             index=index+1
          end do
       else
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               grad_mat_xa
          pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
               grad_mat_ya
          pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
               grad_mat_za
       endif
    else
       grad_dim = 0 ! required for compact storage mode
    endif

    ! LOAD dsym/dRb < xi_i | ... | xi_j >
    ! Input: grad_mat_gr(...,4:6) , Output: pointer_prim_int(ib:) or
    !                                       pointer_prim_int(N+1:M)
    if (moving_b) then
       if (present(compact_storage)) then
          index = grad_dim + 1
       else
          index = gradient_index(imb)
       endif
       grad_dim = gradient_index(imb+1) - gradient_index(imb)
       if(grad_dim==3) then
          if(sum((unique_atom_grad_info(imb)%m(:,:,i_ea2)-unity_matrix)**2)<&
               1.0e-7_r8_kind) then
             do_rotation=.false.
          else
             do_rotation=.true.
             rotmat=>unique_atom_grad_info(imb)%m(:,:,i_ea2)
          endif
       else
          do_rotation=.true.
          rotmat=>unique_atom_grad_info(imb)%m(:,:,i_ea2)
       end if

       if(do_rotation) then
          do i_grad=1,grad_dim ! only if moving_b
             pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
                  rotmat(i_grad,1)*grad_mat_xb+&
                  rotmat(i_grad,2)*grad_mat_yb+&
                  rotmat(i_grad,3)*grad_mat_zb
             index=index+1
          end do
       else
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               grad_mat_xb
          pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
               grad_mat_yb
          pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
               grad_mat_zb
       endif
    endif
  end subroutine add_grads

  subroutine calc_sym_coef()
    ! Purpose: piece of the code
    !          calculate symmetry coefficients for the l-type 3 center
    !          fitintegrals
    ang_momentum_symadapt: do i_l=1,lmax_ch
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       do i_ind=1,n_indep_fcts
          n_contributing_fcts = &
               unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%N_fcts
          magn => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          coef =>  unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          eq_atom => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%I_equal_atom

          do i_cont=1,n_contributing_fcts
             lm = i_l**2 + magn(i_cont)
             sym_coef1(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coef1(:,eq_atom(i_cont),i_ind,i_l) + &
                  yl_arr(:,lm,eq_atom(i_cont))*&
                  coef(i_cont)
             do i = 1,3
                sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) = &
                     sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) + &
                     yl_arr_grad(:,eq_atom(i_cont),lm,i)*coef(i_cont)
             enddo
          enddo
       enddo
    enddo ang_momentum_symadapt
  end subroutine calc_sym_coef

  subroutine l_coulomb()
    ! Purpose: claculate gradients of l-type fit integrals
    real(kind=r8_kind), dimension(num) :: help_vec1, help_vec2, &
         help_vec3, help_vec4, help_vec5, help_vec6, help_vec7
    ang_momentum: do i_l=1,lmax_ch
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms(ua3)%l_ch(i_l)%N_exponents
       cexps => unique_atoms(ua3)%l_ch(i_l)%exponents(:)
       do i_ind=1,n_indep_fcts
          equal_3_l: do i_ea3=1,n_equal_c
             if (moving_c) then
                if(grad_dim==3) then
                   if(sum((unique_atom_grad_info(imc)%m(:,:,i_ea3)-unity_matrix)**2)<&
                        1.0e-7_r8_kind) then
                      do_rotation=.false.
                   else
                      do_rotation=.true.
                      rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
                   endif
                else
                   do_rotation=.true.
                   rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
                end if
             endif
             xc =  unique_atoms(ua3)%position(:,i_ea3)
             check_ab =( (quadrupel%ua1.eq.quadrupel%ua2).and. &
                  (i_ea1.eq.i_ea2))
             check_ac =( (quadrupel%ua1.eq.ua3).and. &
                  (i_ea1.eq.i_ea3))
             check_bc =( (quadrupel%ua2.eq.ua3).and. &
                  (i_ea2.eq.i_ea3))

             ! if all three centers are the same, the integrals
             ! are zero
             if ( check_ab.and.check_bc ) then
                cycle equal_3_l
             endif


             do k = 1,ncexps
                gamma_arg2 = cexps(k)*fact0/(fact0 + cexps(k))* &
                     (  (gamma_arg(:,1) - xc(1))**2 + &
                     (gamma_arg(:,2) - xc(2))**2 + &
                     (gamma_arg(:,3) - xc(3))**2) 
                gamma_help(:,1:lmax_ch+2) = gamma(lmax_ch+2,gamma_arg2)
                help_vec1=two*pi/cexps(k)*&
                     sqrt( fact0/(fact0+cexps(k)) ) * &
                     (two*cexps(k)*fact0/(fact0+cexps(k)) )**i_l
                help_vec2=overlap1*gamma_help(:,i_l+1)/fact0
                help_vec3=overlap1*two*cexps(k)/(fact0+cexps(k))* &
                     gamma_help(:,i_l+2)*sym_coef1(:,i_ea3,i_ind,i_l)
                help_vec5=sym_coef1(:,i_ea3,i_ind,i_l)* gamma_help(:,i_l+1)
                do i = 1,3
                   help_vec4=help_vec3*(gamma_arg(:,i)-xc(i))
                   help_vec6=help_vec5*grad_overlap_a(:,i)
                   help_vec7=help_vec2*sym_coef2(:,i_ea3,i_ind,i_l,i)
                   if (moving_a) &
                        coul_int_a(i)%l(i_l)%m(:,k,i_ind,m1,m2) = &
                        coul_int_a(i)%l(i_l)%m(:,k,i_ind,m1,m2) + &
                        help_vec1* ( help_vec6 + &
                        help_vec7*aexp_arr - help_vec4*aexp_arr)
                   if (moving_b) &
                        coul_int_b(i)%l(i_l)%m(:,k,i_ind,m1,m2) = &
                        coul_int_b(i)%l(i_l)%m(:,k,i_ind,m1,m2)  + &
                        help_vec1* ( -help_vec6 + &
                        help_vec7*bexp_arr-help_vec4*bexp_arr)
                   if (moving_c) &
                        help_arr(:,i) = &
                        help_vec1* ( -help_vec7*fact0 + &
                        fact0*help_vec4)
                enddo
                if (moving_c) then
                   if(do_rotation) then
                      ! make gradient totalsymmetric before adding
                      do i_grad=1,grad_dim ! only if moving_c
                         coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)=&
                              coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)+&
                              rotmat(i_grad,1)*help_arr(:,1)+&
                              rotmat(i_grad,2)*help_arr(:,2)+&
                              rotmat(i_grad,3)*help_arr(:,3)
                      enddo
                   else
                      coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)=&
                           coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)+&
                           help_arr(:,1)
                      coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)=&
                           coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)+&
                           help_arr(:,2)
                      coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)=&
                           coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)+&
                           help_arr(:,3)
                   end if
                endif
             enddo! loop over k
          enddo equal_3_l
       enddo
    enddo ang_momentum
  end subroutine l_coulomb
  !**************************************************************

  !**************************************************************
  subroutine r2_coulomb()
    ! Purpose: calculate gradients of r2 coulomb integrals
    !---** r2-type Fitfct. **---
    real(kind=r8_kind), dimension(num) :: help_vec1, help_vec2, help_vec3, &
         help_vec4, help_vec5
    ncexps = unique_atoms(ua3)%r2_ch%n_exponents
    cexps => unique_atoms(ua3)%r2_ch%exponents(:)

    do k=1,ncexps

       gamma_arg2 = cexps(k)*fact0/(fact0 + cexps(k))* &
            (  (gamma_arg(:,1) - xc(1))**2 + &
            (gamma_arg(:,2) - xc(2))**2 + &
            (gamma_arg(:,3) - xc(3))**2)
       gamma_help(:,1:3) = gamma(3,gamma_arg2)

       fact3 = (fact0+three/two*cexps(k))/fact0   ! (a+b+3/2c)/(a+b)
       help_vec1=two*pi/(cexps(k)**2)*(fact0/(fact0+cexps(k)))* &
            sqrt( fact0/(fact0+cexps(k)))
       help_vec3= two*(((1-fact3)*gamma_help(:,2)-gamma_arg2*gamma_help(:,3))*&
            cexps(k)/(fact0+cexps(k)))*overlap1
       help_vec4=(fact3*gamma_help(:,1)+gamma_arg2*gamma_help(:,2))
       if(check_ab) then
          help_vec1=help_vec1*help_vec3
          do i = 1,3
             help_vec2=help_vec1*(gamma_arg(:,i)-xc(i))
             if (moving_a) &
                  coul_int_a(i)%l(-1)%m(:,k,1,m1,m2) = &
                  coul_int_a(i)%l(-1)%m(:,k,1,m1,m2) + & 
                  help_vec2*aexp_arr
             if (moving_b) &
                  coul_int_b(i)%l(-1)%m(:,k,1,m1,m2) = &
                  coul_int_b(i)%l(-1)%m(:,k,1,m1,m2) + &
                  help_vec2*bexp_arr
             if (moving_c) &
                  help_arr(:,i) = -help_vec2*fact0
          enddo
       else
          do i = 1,3
             help_vec2=help_vec3*(gamma_arg(:,i)-xc(i))
             help_vec5=help_vec4*grad_overlap_a(:,i)
             if (moving_a) &
                  coul_int_a(i)%l(-1)%m(:,k,1,m1,m2) = &
                  coul_int_a(i)%l(-1)%m(:,k,1,m1,m2) + & 
                  help_vec1*(help_vec2*aexp_arr+help_vec5)
             if (moving_b) &
                  coul_int_b(i)%l(-1)%m(:,k,1,m1,m2) = &
                  coul_int_b(i)%l(-1)%m(:,k,1,m1,m2) + &
                  help_vec1*(help_vec2*bexp_arr-help_vec5)
             if (moving_c) &
                  help_arr(:,i) = -help_vec1*help_vec2*fact0
          end do
       end if

       if (moving_c) then
          if(do_rotation) then
             ! make gradient totalsymmetric before adding
             do i_grad=1,grad_dim ! only if moving_c
                coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)=&
                     coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)+&
                     rotmat(i_grad,1)*help_arr(:,1)+&
                     rotmat(i_grad,2)*help_arr(:,2)+&
                     rotmat(i_grad,3)*help_arr(:,3)
             enddo
          else
             coul_int_c(1)%l(-1)%m(:,k,1,1,1)=&
                  coul_int_c(1)%l(-1)%m(:,k,1,1,1)+&
                  help_arr(:,1)
             coul_int_c(2)%l(-1)%m(:,k,1,1,1)=&
                  coul_int_c(2)%l(-1)%m(:,k,1,1,1)+&
                  help_arr(:,2)
             coul_int_c(3)%l(-1)%m(:,k,1,1,1)=&
                  coul_int_c(3)%l(-1)%m(:,k,1,1,1)+&
                  help_arr(:,3)
          end if
       endif
    enddo! r2-exponents, third center
  end subroutine r2_coulomb
  !**************************************************************

  !**************************************************************
  subroutine s_coulomb
    ! Purpose: calculate gradients of s type coulomb fitintegrals
    ! loop over exponents of third center ---** s-type Fitfct. **---
    real(kind=r8_kind), dimension(num) :: help_vec1, help_vec2, help_vec3, help_vec4
    ncexps = unique_atoms(ua3)%l_ch(0)%n_exponents
    cexps => unique_atoms(ua3)%l_ch(0)%exponents(:)
    do k=1,ncexps

       gamma_arg2 = cexps(k)*fact0/(fact0 + cexps(k))* &
            (  (gamma_arg(:,1) - xc(1))**2 + &
            (gamma_arg(:,2) - xc(2))**2 + &
            (gamma_arg(:,3) - xc(3))**2)
       gamma_help(:,1:2) = gamma(2,gamma_arg2)

       help_vec1=four*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))*overlap1
       help_vec3=gamma_help(:,2)*cexps(k)/(fact0+cexps(k))
       if ( check_ab ) then ! position a = position b
          do i = 1,3
             help_vec2=help_vec1*help_vec3*(gamma_arg(:,i)-xc(i))
             if (moving_a) &
                  coul_int_a(i)%l(0)%m(:,k,1,m1,m2) = &
                  coul_int_a(i)%l(0)%m(:,k,1,m1,m2)- & 
                  aexp_arr*help_vec2

             if (moving_b) &
                  coul_int_b(i)%l(0)%m(:,k,1,m1,m2) = &
                  coul_int_b(i)%l(0)%m(:,k,1,m1,m2)- &
                  bexp_arr*help_vec2

             if (moving_c) &
                  help_arr(:,i) = &
                  fact0*help_vec2

          enddo
       else   ! position a != position b 
          do i = 1,3
             help_vec2=help_vec3*(gamma_arg(:,i)-xc(i))
             help_vec4=gamma_help(:,1)*fact2*dist(i)
             if (moving_a) &
                  coul_int_a(i)%l(0)%m(:,k,1,m1,m2) = &
                  coul_int_a(i)%l(0)%m(:,k,1,m1,m2)- &
                  help_vec1*(aexp_arr*help_vec2+help_vec4)

             if (moving_b) &
                  coul_int_b(i)%l(0)%m(:,k,1,m1,m2) = &
                  coul_int_b(i)%l(0)%m(:,k,1,m1,m2)- &
                  help_vec1*(bexp_arr*help_vec2-help_vec4)

             if (moving_c) &
                  help_arr(:,i) = & 
                  fact0*help_vec2*help_vec1
          enddo
       endif
       if (moving_c) then
          if(do_rotation) then
             ! make gradient totalsymmetric before adding
             do i_grad=1,grad_dim ! only if moving_c
                coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)=&
                     coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)+&
                     rotmat(i_grad,1)*help_arr(:,1)+&
                     rotmat(i_grad,2)*help_arr(:,2)+&
                     rotmat(i_grad,3)*help_arr(:,3)
             enddo
          else
             coul_int_c(1)%l(0)%m(:,k,1,1,1)=&
                  coul_int_c(1)%l(0)%m(:,k,1,1,1)+&
                  help_arr(:,1)
             coul_int_c(2)%l(0)%m(:,k,1,1,1)=&
!  integer(kind=i4_kind)::k_grad
                  coul_int_c(2)%l(0)%m(:,k,1,1,1)+&
                  help_arr(:,2)
             coul_int_c(3)%l(0)%m(:,k,1,1,1)=&
                  coul_int_c(3)%l(0)%m(:,k,1,1,1)+&
                  help_arr(:,3)
          end if
       endif
!  integer(kind=i4_kind)::k_grad
    enddo! s-exponents of third center 
  end subroutine s_coulomb
  !**************************************************************

  !**************************************************************
  subroutine calculate_nuc
    ! Purpose: calculate gradients of nuclear attraction
    !          and if necessary pvsp
!  integer(kind=i4_kind)::k_grad

    equal_3: do i_ea3 = 1,n_equal_c
       if (moving_c) then
          if(grad_dim==3) then
             if(sum((unique_atom_grad_info(imc)%m(:,:,i_ea3)-unity_matrix)**2)<&
                  1.0e-7_r8_kind) then
                do_rotation=.false.
             else
                do_rotation=.true.
                rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
             endif
          else
             do_rotation=.true.
             rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
          end if

       endif
       xc = unique_atoms(ua3)%position(:,i_ea3)

       check_ab =( (quadrupel%ua1.eq.quadrupel%ua2).and. &
            (i_ea1.eq.i_ea2))
       check_ac =( (quadrupel%ua1.eq.ua3).and. &
            (i_ea1.eq.i_ea3))
       check_bc =( (quadrupel%ua2.eq.ua3).and. &

            (i_ea2.eq.i_ea3))

       ! if all three centers are the same, the integrals
       !----------------------------------------------------------------
!       if(.true..and.cpks_energies) then ! no nuc grud but just nuc
!        print*,'nuclear ss grad',sum(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi))
!                 do k_grad=1,size(prim_int_3cob_grad)
!                 prim_int_3cob_grad(k_grad)%m(:,:,1,1)=&
!                      prim_int_3cob_grad(k_grad)%m(:,:,1,1) &
!                            +unpack(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi),cutoff,zero)
!                 enddo
!       print*,'nuclear ss grad',sum(prim_int_3cob_grad(2)%m(:,:,1,1))
!       endif
!      !--------------------------------------------------
!       !----------------------------------------------------------------
!       if(.true..and.cpks_energies) then ! no nuc grud but just nuc
!        print*,'nuclear ss grad',sum(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi))
!                 do k_grad=1,size(prim_int_3cob_grad)
!                 prim_int_3cob_grad(k_grad)%m(:,:,1,1)=&
!                      prim_int_3cob_grad(k_grad)%m(:,:,1,1) &
!                            +unpack(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi),cutoff,zero)
!                 enddo
!       print*,'nuclear ss grad',sum(prim_int_3cob_grad(2)%m(:,:,1,1))
!       endif
!      !--------------------------------------------------
       ! are zero
       if ( check_ab.and.check_bc ) then
          cycle equal_3
       endif
       gamma_arg2 = ( (gamma_arg(:,1)-xc(1))**2 + &
                      (gamma_arg(:,2)-xc(2))**2 + &
                      (gamma_arg(:,3)-xc(3))**2  )

       gamma_help(:,1:max_gamma) = gamma(max_gamma,fact0*gamma_arg2)

       ! grad_nuclear_a : derivative with respect to  center a:
       !                  < grad_a mu | V_nuc | nu >
       ! grad_nuclear_b :
       !                  < mu | V_nuc | grad_b nu >
       ! sign convention: '-' due to attractive interaction
       !                  between electrons and nuclei
       !                  is INCLUDED here.
       !----------------------------------------------------------------
!       if(.true..and.cpks_energies) then ! no nuc grud but just nuc
!        print*,'nuclear ss grad',sum(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi))
!                 do k_grad=1,size(prim_int_3cob_grad)
!                 prim_int_3cob_grad(k_grad)%m(:,:,1,1)=&
!                      prim_int_3cob_grad(k_grad)%m(:,:,1,1) &
!                            +unpack(two*zc*overlap1*gamma_help(:,1)*sqrt(fact0/pi),cutoff,zero)
!                 enddo
!       print*,'nuclear ss grad',sum(prim_int_3cob_grad(2)%m(:,:,1,1))
!       endif
      !--------------------------------------------------
       if ( .not. check_ab ) then
          do i = 1,3
            pprel1: if(iam_ppot.and.integralpar_relativistic) then
             if (moving_a) &
                  grad_nuc_pseudo_a(:,i) = grad_nuc_pseudo_a(:,i) - &
                  ( &
                  pack(spread(aexps,1,nbexps),cutoff)* &
                  (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) + &
                  fact2 * gamma_help(:,1) * dist(i) ) * &
                  fac*overlap1*four*zc
             if (moving_b) &
                  grad_nuc_pseudo_b(:,i) = grad_nuc_pseudo_b(:,i) - &
                  four*zc*fac*overlap1*( &
                  pack(spread(bexps,2,naexps),cutoff)* &
                  gamma_help(:,2)*(gamma_arg(:,i)-xc(i)) - &
                  fact2*gamma_help(:,1)*dist(i))
             
            else pprel1 !regular atom or any atom in non rel calc
             if (moving_a) &
                  grad_nuclear_a(:,i) = grad_nuclear_a(:,i) - &
                  ( &
                  pack(spread(aexps,1,nbexps),cutoff)* &
                  (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) + &
                  fact2 * gamma_help(:,1) * dist(i) ) * &
                  fac*overlap1*four*zc
             if (moving_b) &
                  grad_nuclear_b(:,i) = grad_nuclear_b(:,i) - &
                  four*zc*fac*overlap1*( &
                  pack(spread(bexps,2,naexps),cutoff)* &
                  gamma_help(:,2)*(gamma_arg(:,i)-xc(i)) - &
                  fact2*gamma_help(:,1)*dist(i))           
             endif pprel1

            if (moving_c) &
                  grad_nuclear_c(:,i) = four*zc*fac*overlap1*&
                  fact0*gamma_help(:,2)*(gamma_arg(:,i)-xc(i)) 
       end do
       
       
       else
          do i = 1,3
             if(iam_ppot.and.integralpar_relativistic) then
                if (moving_a .or. moving_b) then
                   grad_nuc_pseudo_a(:,i) = grad_nuc_pseudo_a(:,i) - &
                        fact0 * &
                        (gamma_arg(:,i)-xc(i)) * gamma_help(:,2)* &
                        fac*overlap1*2.0_r8_kind*zc
                   grad_nuc_pseudo_b(:,i) = grad_nuc_pseudo_a(:,i)
                   ! grad_c is needed only for ONE EQUAL
                endif
             else ! not iam_ppot
                if (moving_a .or. moving_b) then
                   grad_nuclear_a(:,i) = grad_nuclear_a(:,i) - &
                        fact0 * &
                        (gamma_arg(:,i)-xc(i)) * gamma_help(:,2)* &
                        fac*overlap1*2.0_r8_kind*zc
                   grad_nuclear_b(:,i) = grad_nuclear_a(:,i)
                   ! grad_c is needed only for ONE EQUAL
                endif
          
             ! Bei CO, Ni4 (TD) richtig
          end if
                  if (moving_c) &
                  grad_nuclear_c(:,i) = +four*zc*fac*overlap1*&
                  fact0*gamma_help(:,2)*(gamma_arg(:,i)-xc(i))            
          enddo
       endif
       ! make derivatives with respect to center c totalsymetric
       if (moving_c) then
          if(do_rotation) then
             counter=1
             do i_grad=1,grad_dim ! only if moving_c
                nuc_grad(:,i_grad)=nuc_grad(:,i_grad)+&
                     rotmat(i_grad,1)*grad_nuclear_c(:,1)+&
                     rotmat(i_grad,2)*grad_nuclear_c(:,2)+&
                     rotmat(i_grad,3)*grad_nuclear_c(:,3)
             enddo
             counter=counter+1
          else
             counter=1
             nuc_grad(:,1)=nuc_grad(:,1)+&
                  grad_nuclear_c(:,1)
             nuc_grad(:,2)=nuc_grad(:,2)+&
                  grad_nuclear_c(:,2)
             nuc_grad(:,3)=nuc_grad(:,3)+&
                  grad_nuclear_c(:,3)
             counter=counter+1
          end if
       endif
    enddo equal_3

  end subroutine calculate_nuc

  subroutine calculate_solv

    real(kind=r8_kind) :: z
    integer(kind=i4_kind) :: n_equal_solv,ism,N_pc
    integer(kind=i4_kind) :: i,j,k,l,m

    N_pc=0
    if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

    allocate( grad_solv_a(num,3),STAT=alloc_stat )  
    if (alloc_stat.ne.0) call error_handler &            
         ("SS_CALCULATE_GRAD: allocation grad_solv_a failed")  
    grad_solv_a = zero                                 
    allocate( grad_solv_b(num,3),STAT=alloc_stat )    
    if (alloc_stat.ne.0) call error_handler &             
         ("SS_CALCULATE_GRAD: allocation grad_solv_b failed")
    grad_solv_b = zero                                    
    allocate( grad_solv_c(num,3),grad_solv_c_help(num,3), &
                                        STAT=allocstat(7) )       
    if (allocstat(7).ne.0) call error_handler &        
         ("SS_CALCULATE_GRAD: allocation solv_c failed") 
        allocstat(7)=1
    grad_solv_c = zero                                  
    grad_solv_c_help = zero                               

    check_ab =( (quadrupel%ua1.eq.quadrupel%ua2).and. &
         (i_ea1.eq.i_ea2))

    n_size_s: do j=1,to_calc_grads%n_points 
       n_equal_solv=to_calc_grads%n_equal(j)
       z=to_calc_grads%Q(j)

       n_equal_s: do k=1,n_equal_solv
          xc =to_calc_grads%xyz(j,k,:)  
          gamma_arg2 = ( (gamma_arg(:,1)-xc(1))**2 + &
               (gamma_arg(:,2)-xc(2))**2 + &
               (gamma_arg(:,3)-xc(3))**2  )
          gamma_help(:,1:max_gamma) = gamma(max_gamma,fact0*gamma_arg2)	

          ! grad_solv_a : derivative with respect to  center a:
          !                  < grad_a mu | V_nuc | nu >
          ! grad_solv_b :
          !                  < mu | V_nuc | grad_b nu >
          ! sign convention: '-' due to attractive interaction
          !                  between electrons and nuclei
          !                  is INCLUDED here.
          if ( .not. check_ab ) then
             do i = 1,3
                if (moving_a) &
                     grad_solv_a(:,i) = grad_solv_a(:,i) - &
                     ( pack(spread(aexps,1,nbexps),cutoff)* &
                     (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) + &
                     fact2 * gamma_help(:,1) * dist(i) ) * &
                     fac*overlap1*four*z
                if (moving_b) &
                     grad_solv_b(:,i) = grad_solv_b(:,i) - &
                     ( pack(spread(bexps,2,naexps),cutoff)* &
                     (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) - &
                     fact2 * gamma_help(:,1) * dist(i) ) * &
                     fac*overlap1*four*z
             enddo
          else
             do i = 1,3
                if (moving_a .or. moving_b) then
                   grad_solv_a(:,i) = grad_solv_a(:,i) - &
                        fact0 * &
                        (gamma_arg(:,i)-xc(i)) * gamma_help(:,2)* &
                        fac*overlap1*2.0_r8_kind*z
                   grad_solv_b(:,i) = grad_solv_a(:,i)
                endif
             enddo
          endif

!*******************************************
          ism=to_calc_grads%i_symm_sort(j,k)
          unique_3: do ua3 = 1,N_moving_unique_atoms+N_pc
             if(ua3 <= N_moving_unique_atoms) then
                i=moving_unique_atom_index(ua3)
                n_equal_c = unique_atoms(i)%N_equal_atoms
                grad_dim = gradient_index(ua3+1) - gradient_index(ua3)
             else
                i=ua3-N_moving_unique_atoms
                n_equal_c=pointcharge_array(i)%N_equal_charges
                grad_dim = surf_points_grad_index(i+1)-surf_points_grad_index(i)
             end if

             grad_solv_c=zero
             do i = 1,3
                grad_solv_c_help(:,i) =four*fac*overlap1*&
                     fact0*gamma_help(:,2)* &
                     (gamma_arg(:,i)-xc(i))*z
             enddo
             do l=1,grad_dim
                do m=1,3
                   grad_solv_c(:,l)=grad_solv_c(:,l)+ &
                        grad_solv_c_help(:,m)*&
                        to_calc_grads%dxyz_totsyms(l,ua3)%m(m,ism)
                enddo
             enddo

             if(ua3 <= N_moving_unique_atoms) then
                index = gradient_index(ua3) - 1
                do i_grad=1,grad_dim 
                   prim_int_3cob_solv_grad(index+i_grad)%m(:,:,1,1)=&    
                        prim_int_3cob_solv_grad(index+i_grad)%m(:,:,1,1)&
			-unpack(grad_solv_c(:,i_grad),cutoff,zero)
                enddo
             else
                index=surf_points_grad_index(ua3-N_moving_unique_atoms) - 1
                do i_grad=1,grad_dim 
                   prim_int_3cob_solv_grad_pc(index+i_grad)%m(:,:,1,1)=&    
                        prim_int_3cob_solv_grad_pc(index+i_grad)%m(:,:,1,1)&
			+unpack(grad_solv_c(:,i_grad),cutoff,zero)
                enddo
             end if
    
      enddo unique_3

!******************************************
       enddo n_equal_s
    enddo n_size_s


    ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
    if (moving_a) then
       grad_mat_xa(:,:,1,1)=grad_mat_xa(:,:,1,1)-&
            unpack(grad_solv_a(:,1),cutoff,zero)
       grad_mat_ya(:,:,1,1)=grad_mat_ya(:,:,1,1)-&
            unpack(grad_solv_a(:,2),cutoff,zero)
       grad_mat_za(:,:,1,1)=grad_mat_za(:,:,1,1)-&
            unpack(grad_solv_a(:,3),cutoff,zero)
    endif
    if (moving_b) then
       grad_mat_xb(:,:,1,1)=grad_mat_xb(:,:,1,1)-&
            unpack(grad_solv_b(:,1),cutoff,zero)
       grad_mat_yb(:,:,1,1)=grad_mat_yb(:,:,1,1)-&
            unpack(grad_solv_b(:,2),cutoff,zero)
       grad_mat_zb(:,:,1,1)=grad_mat_zb(:,:,1,1)-&
            unpack(grad_solv_b(:,3),cutoff,zero)
    endif

    deallocate(grad_solv_a,grad_solv_b,stat=alloc_stat)
    if(alloc_stat/=0) call error_handler &
         ("SS_CALCULATE: deallocation grad_solv_a failed")
    
    pointer_prim_int=>prim_int_3cob_solv_grad
    call add_grads()
    if (model_density .and. spin_polarized) then
       spin_index = gradient_data_n_gradients + 1
       pointer_prim_int=>prim_int_3cob_solv_grad(spin_index:)
       call grad_spmat_to_grad_mat()
       call add_grads()
    endif
!!$do i=1,size(prim_int_3cob_solv_grad)
!!$print*,na,nb,'a',i,prim_int_3cob_solv_grad(i)%m(1,1,1,1)
!!$end do
!!$if(with_pc .and. .not.fixed_pc) then
!!$do i=1,size(prim_int_3cob_solv_grad_pc)
!!$print*,na,nb,'c',i,prim_int_3cob_solv_grad_pc(i)%m(1,1,1,1)
!!$end do
!!$end if

    deallocate(&
         grad_mat_xa,grad_mat_xb,&
         grad_mat_ya,grad_mat_yb,&
         grad_mat_za,grad_mat_zb,stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
    

    deallocate (aexp_arr,bexp_arr,fact0,fact1,fact2,fact3,fac,help_arr, &
                cutoff, &
                gamma_help,gamma_arg,gamma_arg2, &
                overlap1, grad_overlap_a, &
                grad_solv_c,grad_solv_c_help,STAT=allocstat(3))
    if (allocstat(1).ne.0) call error_handler &               ! aexp_arr bexp_arr fact0 fact1 
         ("SS_CALCULATE_GRAD: allocation of fact0 failed(1)") ! fact2 fact3 fac help_arr
        allocstat(2)=0 ! cutoff
        allocstat(4)=0 ! gamma_help gamma_arg gamma_arg2
        allocstat(5)=0 ! overlap1 grad_overlap_a
        allocstat(7)=0 ! grad_solv_c grad_solv_c_help

  end subroutine calculate_solv

  subroutine calculate_nuc_pc
    ! Purpose: calculate gradients of nuclear attraction on point charges
    !          and if necessary pvsp

    unique_3_PC: do ua3 = 1, pointcharge_N+n_timps
       if(ua3<=n_timps) then
          ua_pointer=> unique_timps(ua3)
          imc = ua_pointer%moving_atom
!!$        moving_c=.false. !!!! not false for EPE calculations
          moving_c = imc > 0
          if (moving_c) then
             grad_dim = gradient_index(N_moving_unique_atoms+imc+1) - &
                  gradient_index(N_moving_unique_atoms+imc)
             allocate(nuc_grad(num,grad_dim),stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("ss_calculate_grad: allocate nuc_grad failed")
             nuc_grad=0.0_r8_kind
          else
             grad_dim = 0
          endif
          zc = unique_timps(ua3)%Z - unique_timps(ua3)%zc
          n_equal_c = unique_timps(ua3)%N_equal_atoms          
       else
          cycle unique_3_PC ! pointcharges go to SHGI !!!!!!!!!!!AS
          moving_c=.false.
          zc = pointcharge_array(ua3-n_timps)%Z
          n_equal_c = pointcharge_array(ua3-n_timps)%N_equal_charges
       end if
       pointcharges: do i_ea3 = 1,n_equal_c
       if (moving_c) then
          if(grad_dim==3) then
             if(sum((unique_timp_grad_info(imc)%m(:,:,i_ea3)-unity_matrix)**2)<&
                  1.0e-7_r8_kind) then
                do_rotation=.false.
             else
                do_rotation=.true.
                rotmat=>unique_timp_grad_info(imc)%m(:,:,i_ea3)
             endif
          else
             do_rotation=.true.
             rotmat=>unique_timp_grad_info(imc)%m(:,:,i_ea3)
          end if
       endif
          if(ua3<=n_timps) then
             xc = unique_timps(ua3)%position(:,i_ea3)
          else
             xc = pointcharge_array(ua3-n_timps)%position(:,i_ea3)
          endif
          check_ab =( (quadrupel%ua1.eq.quadrupel%ua2).and. &
               (i_ea1.eq.i_ea2))

          gamma_arg2 = ( (gamma_arg(:,1)-xc(1))**2 + &
               (gamma_arg(:,2)-xc(2))**2 + &
               (gamma_arg(:,3)-xc(3))**2  )

          gamma_help(:,1:max_gamma) = gamma(max_gamma,fact0*gamma_arg2)

          ! grad_nuclear_a : derivative with respect to  center a:
          !                  < grad_a mu | V_nuc | nu >
          ! grad_nuclear_b :
          !                  < mu | V_nuc | grad_b nu >
          ! sign convention: '-' due to attractive interaction
          !                  between electrons and nuclei
          !                  is INCLUDED here.
          if ( .not. check_ab ) then
             do i = 1,3

                if(pseudopot_present.and.integralpar_relativistic) then 
                   if (moving_a) &
                        grad_nuc_pseudo_a(:,i) = grad_nuc_pseudo_a(:,i) - &
                        ( &
                        pack(spread(aexps,1,nbexps),cutoff)* &
                        (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) + &
                        fact2 * gamma_help(:,1) * dist(i) ) * &
                        fac*overlap1*four*zc
                   if (moving_b) &
                        grad_nuc_pseudo_b(:,i) = grad_nuc_pseudo_b(:,i) - &
                        four*zc*fac*overlap1*( &
                        pack(spread(bexps,2,naexps),cutoff)* &
                        gamma_help(:,2)*(gamma_arg(:,i)-xc(i)) - &
                        fact2*gamma_help(:,1)*dist(i)) 
                else
                   if (moving_a) &
                        grad_nuclear_a(:,i) = grad_nuclear_a(:,i) - &
                        ( &
                        pack(spread(aexps,1,nbexps),cutoff)* &
                        (gamma_arg(:,i)-xc(i)) * gamma_help(:,2) + &
                        fact2 * gamma_help(:,1) * dist(i) ) * &
                        fac*overlap1*four*zc
                   if (moving_b) &
                        grad_nuclear_b(:,i) = grad_nuclear_b(:,i) - &
                        four*zc*fac*overlap1*( &
                        pack(spread(bexps,2,naexps),cutoff)* &
                        gamma_help(:,2)*(gamma_arg(:,i)-xc(i)) - &
                        fact2*gamma_help(:,1)*dist(i))
                   if (moving_c) &
                        grad_nuclear_c(:,i) = four*zc*fac*overlap1* &
                        fact0*gamma_help(:,2)* (gamma_arg(:,i)-xc(i)) 
                end if

             enddo
          else
             if (moving_a .or. moving_b) then

                do i = 1,3
                   if(pseudopot_present.and.integralpar_relativistic) then
                      grad_nuc_pseudo_a(:,i) = grad_nuc_pseudo_a(:,i) - &
                           fact0 * &
                           (gamma_arg(:,i)-xc(i)) * gamma_help(:,2)* &
                           fac*overlap1*2.0_r8_kind*zc
                      grad_nuc_pseudo_b(:,i) = grad_nuc_pseudo_a(:,i)
                   else
                      grad_nuclear_a(:,i) = grad_nuclear_a(:,i) - &
                           fact0 * &
                           (gamma_arg(:,i)-xc(i)) * gamma_help(:,2)* &
                           fac*overlap1*2.0_r8_kind*zc
                      grad_nuclear_b(:,i) = grad_nuclear_a(:,i)
                   end if
                   if (moving_c) &
                        grad_nuclear_c(:,i) = +four*zc*fac*overlap1*&
                        fact0*gamma_help(:,2)*(gamma_arg(:,i)-xc(i))            

                enddo
             endif
          endif
       ! make derivatives with respect to center c totalsymetric
          if (moving_c) then
             if(do_rotation) then
                counter=1
                do i_grad=1,grad_dim ! only if moving_c
                   nuc_grad(:,i_grad)=nuc_grad(:,i_grad)+&
                        rotmat(i_grad,1)*grad_nuclear_c(:,1)+&
                        rotmat(i_grad,2)*grad_nuclear_c(:,2)+&
                        rotmat(i_grad,3)*grad_nuclear_c(:,3)
                enddo
                counter=counter+1
             else
                counter=1
                nuc_grad(:,1)=nuc_grad(:,1)+&
                     grad_nuclear_c(:,1)
                nuc_grad(:,2)=nuc_grad(:,2)+&
                     grad_nuclear_c(:,2)
                nuc_grad(:,3)=nuc_grad(:,3)+&
                     grad_nuclear_c(:,3)
                counter=counter+1
             end if
          endif
       enddo pointcharges

       if(moving_c) then
          index = gradient_index(imc+n_moving_unique_atoms) - 1
          do i_grad=1,grad_dim ! only if moving_c
           prim_int_3cob_grad(index+i_grad)%m(:,:,1,1)=&
               prim_int_3cob_grad(index+i_grad)%m(:,:,1,1)- &
                       unpack(nuc_grad(:,i_grad),cutoff,zero)
          end do

          deallocate(nuc_grad,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("SS_CALCULATE_GRAD: deallocation nuc_grad for nuc timps failed")
       end if

    enddo unique_3_PC

  end subroutine calculate_nuc_pc



  subroutine calculate_ol_and_kin
    ! overlap and gradient of overlap --------------------

    allocate(overlap1(num),grad_overlap_a(num,3), STAT=allocstat(5))
    ASSERT(allocstat(5).eq.0)
    allocstat(5)=1

    overlap1 = (two*sqrt(fact1)/fact0)* &
         sqrt((two*sqrt(fact1)/fact0))*exp(-fact2*arg)
    grad_overlap_a = zero

    allocate(grad_kinetic(num,3),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    grad_kinetic = zero

    gr_ind: do i=1,3
       dist(i) = xa(i)-xb(i)
       grad_overlap_a(:,i) = -two*fact2*dist(i)*overlap1
       grad_kinetic(:,i) = - (10.0_r8_kind*fact2**2 - four*fact2**2*fact2*arg) *dist(i) * overlap1 
       !                                                                  !^(a-b)**2
       ! grad_b overlap = -grad_overlap
       ! grad_b kinetic = -grad_kinetic
    enddo gr_ind


    deallocate(grad_kinetic,stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)
  end subroutine calculate_ol_and_kin

end subroutine ss_calculate_grads

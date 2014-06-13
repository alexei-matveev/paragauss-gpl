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
module ls_calculate_grads_module
!===============================================================
! Public interface of module
!===============================================================
  !
  !  Purpose: calculation of gradients of all primitive 2 center orbital
  !           and 3 center integrals for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,equal_atom2).
  !       For three center integrals, contraction and symmetry-
  !       adaption concerning fitfunctions is also performed.
  !
  !
  !  Author: MS
  !  Date:   8/96
  !
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Author: AH
  ! Date:   4/99
  ! Description: gradients for pseudopotentials has 
  !              been added
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
# include "def.h"
  use unique_atom_module, noname=>pseudopot_present
  use gamma_module
  use operations_module, only: operations_core_density
  use type_module
  use datatype
  use solid_harmonics_module, only : solid_harmonics_calc,solid_harmonics_scalar
  use int_data_2cob3c_module  
  use solhrules_module
  use fitcontract_module
  use integralpar_module
  use gradient_data_module
  use options_module, only: options_integral_expmax, options_split_gradients, &
                            options_xcmode, xcmode_model_density, &
                            options_spin_restricted
  use pointcharge_module
#ifdef WITH_EPE
  use ewaldpc_module
#endif
  use solv_cavity_module, only: to_calc_grads, &  !n_size,tessarea,cagr,i_symm_sort,Q_n,Q_e
       with_pc,fixed_pc
  use elec_static_field_module, only: totsym_field_length,surf_points_grad_index
  use calc3c_switches,only: old_3c_co, old_3c_fc, old_solv_grad,integralpar_dervs,new_3c_co_grad
  use cpksdervs_matrices,only: cpksalloc

  implicit none

  !== Interrupt end of public interface of module =================

  ! these globals are set from the sub:
  integer(kind=i4_kind) :: la
  integer(kind=i4_kind) :: equalb ! number of equal atom b

  !================================================================
  ! End of public interface of module
  !================================================================

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
! overlap_grad  (:,:,:,   1:3  ) :    d/dRa < xi_i | 1      | xi_j >
! kin_grad      (  :  ,   1:3  ) :    d/dRa < xi_i | T      | xi_j >
! nuc_grad_gr   (:,:,:,   1:3  ) :    d/dRa < xi_i | V_nuc  | xi_j >
! nuc_grad_gr   (:,:,:,   4:6  ) :    d/dRb < xi_i | V_nuc  | xi_j >
! nuc_grad      (:,:,:,grad_dim) : dsym/dRc < xi_i | V_nuc  | xi_j >
! pvsp_grad_gr  (:,:,:,   1:3  ) :    d/dRa < xi_i | V_pvsp | xi_j >
! pvsp_grad_gr  (:,:,:,   4:6  ) :    d/dRb < xi_i | V_pvsp | xi_j >
! pvsp_grad     (:,:,:,grad_dim) : dsym/dRc < xi_i | V_pvsp | xi_j >
! coul_int_grad       (   1:3  ) :    d/dRa [ xi_i | f_k    | xi_j ]
! coul_int_grad       (   4:6  ) :    d/dRb [ xi_i | f_k    | xi_j ]
! coul_int_grad_totsym(grad_dim) : dsym/dRc [ xi_i | f_k    | xi_j ]
!
! grad_mat_gr  (:,:,:,:,  1:3  ) :    d/dRa < xi_i | ...    | xi_j >
! grad_mat_spgr(:,:,:,:,  4:6  ) :    d/dRb < xi_i | ...    | xi_j >
!
! if ( la < lb ) then
! ~~~~~~~~~~~~~~~~~~~
! 1) overlap_grad and kin_grad contain d/dRb instead of d/dRa
! 2) nuc_grad_gr, pvsp_grad_gr, and coul_int_grad are stored in reverse order,
!    i.e. range 1:3 contains d/dRb and range 4:6 contains d/dRa
! 3) the storage order in grad_mat_gr and grad_mat_spgr is not altered
!..............................................................................

  integer(kind=i4_kind) :: naexps,nbexps,ncexps
  real(kind=r8_kind),pointer :: aexps(:),bexps(:)
  real(kind=r8_kind),pointer     :: cexps(:)
  real(kind=r8_kind) :: z  ! charge
  real(kind=r8_kind) :: zc ! core charge
  integer(kind=i4_kind)          :: max_order,grad_dim,index,i_grad,&
       l_index
  integer(kind=i4_kind)          :: ima, imb, imc, k_gr_0, k_gr_1, &
                                    spin_index, l_max
  logical                        :: moving_a, moving_b, moving_c, &
                                    split_gradients, model_density, &
                                    spin_polarized

logical                        :: iam_ppot
  ! constants
! real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
  real(kind=r8_kind),dimension(3,3) :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  real(kind=r8_kind),dimension(0:8) :: dfac=(/1.0_r8_kind,&
       1.0_r8_kind,3.0_r8_kind,15.0_r8_kind,105.0_r8_kind,945.0_r8_kind,&
       10395.0_r8_kind,135135.0_r8_kind,2027025.0_r8_kind/)
  integer(kind=i4_kind),dimension(3) :: xyz_map
  ! array help 
  logical,allocatable   :: cutoff(:,:)
  integer(kind=i4_kind) :: num,counter,m,ma,alloc_stat,ism
  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: fact0_arr, &
       fact1_arr,fact2_arr,fact10,help_mat
  real(kind=r8_kind),allocatable,dimension(:)  :: fact0,fact1, &
       fact2,fact4,fact5,fact6,fact7,fact8,rcsabc,tau,help_vec,rsaba

  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:)     :: gamma_arg(:,:),&
       gamma_arg2(:),gamma_help(:,:)


  ! help arrays for solid harmincs
  real(kind=r8_kind),allocatable                 :: yl_arr(:,:)
  real(kind=r8_kind),allocatable                 :: clmamb(:,:),clmamb_scalar(:)

  ! help arrays for product_rule and diff_rule
  real(kind=r8_kind),allocatable  :: prod_arr(:,:,:),prod_arr_ls_calc(:,:,:),&
       prod_arr_gr(:,:,:,:),prod_arr_rel_gr(:,:,:,:),&
       diff_arr0(:,:),diff_arr(:,:),&
       diff_arr_gr(:,:,:),&
       diff_arr_grad(:,:,:),&
       intermediate_gr(:,:,:,:,:),&
       cluster_epe(:,:)



  real(kind=r8_kind) :: arg

  ! cartesian coordinates
  real(kind=r8_kind),dimension(3)  :: xa,xb,xc,xd

  logical :: laltlb,lagtlb,&  ! flag to decide if la is lower then lb or not
       do_rotation ! flag to decide if gradients have to be rotated in order 
  ! to make it totalsymmetric

!!$  integer(kind=i4_kind)  :: i,i_la,i_lb,j,i_l,k,k_gr,i_ind,i_cnt,nm_la,nm_lb,&
!!$       la_org,lb_org,la,l,la2pm,k1,n,n1

  integer(kind=i4_kind)  :: nm_la,nm_lb,la_org,lb_org

  integer(kind=i4_kind)  :: lmax_ch,lmax_abs,ly_max
  integer(kind=i4_kind)  :: n_equals,n_independent_fcts, &
       n_contributing_fcts, &
       n_equal_solv !!!!!!!!!!!!!!!!!

  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coeff(:),rotmat(:,:)

  real(kind=r8_kind),allocatable :: &
       overlap(:,:), overlap_grad(:,:,:),&
       kinetic(:),kin_grad(:,:,:), &
       aexp_arr(:),bexp_arr(:)
  real(kind=r8_kind),allocatable :: nuc_grad(:,:,:),&
       nuc_grad_gr(:,:,:), &
       help_arr(:,:),help_arr2(:,:,:),&
       help_arr_gr(:,:,:,:), &
       grad_nuc_pseudo_ab(:,:,:),&
       grad_nuc_pseudo_c(:,:,:), &
       solv_grad(:,:,:),solv_grad_gr(:,:,:), &  
       help_arr_gr1(:,:,:,:)
  real(kind=r8_kind),allocatable :: pseudo_grad_gr(:,:,:), pseudo_grad(:,:,:)

  real(kind=r8_kind),pointer,dimension(:,:,:,:,:) :: & ! help pointers
       grad_mat_gr, grad_mat_spgr, pointer_s, pointer_r2, pointer_l

  real(kind=r8_kind),pointer,dimension(:,:,:,:,:,:) :: dervs_mat

  real(kind=r8_kind),pointer :: grad_mat_p (:,:,:,:) ! pointer on grad_mat_gr
  real(kind=r8_kind),pointer :: grad_mat_sp(:,:,:,:) ! pointer on grad_mat_spgr
!  type(three_center_l), save :: coul_int_grad(6)
  type(three_center_l) :: coul_int_grad(6)
  type(three_center_l),allocatable  :: coul_int_grad_totsym(:)
  ! totalsymmetric gradient
  type(arrmat4),pointer :: pointer_prim_int(:) ! help pointer to point on primitive integrals


  logical                           :: pseudopot_present ! same name as in UA module

contains

  subroutine calc_coul_gradients(index,k,i_ind)
    implicit none
    integer(kind=i4_kind) ::index ! angular momentum of fitting function
    integer(kind=i4_kind) :: k ! number of fitting function
    integer(kind=i4_kind) :: i_ind ! number of independent function
    real(kind=r8_kind),pointer,dimension(:,:,:,:,:) :: pointer_coul
    logical :: eff_moving_a, eff_moving_b
    integer(kind=i4_kind) :: k_gr

    eff_moving_a = (lagtlb .and. moving_a) .or. (laltlb .and. moving_b)
    eff_moving_b = (lagtlb .and. moving_b) .or. (laltlb .and. moving_a)
    ! derivative with respect to first center a
    ! LOAD [ d/dRa xi_i | f_k | xi_j ]
    if (eff_moving_a) then
    do k_gr=1,3
       pointer_coul=>coul_int_grad(k_gr)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)+help_arr_gr(:,k_gr,:,:)
    end do
    endif
 
    ! derivative with respect to third center c
    ! LOAD [ xi_i | dsym/dRc f_k | xi_j ]
    if (moving_c) then
    if(do_rotation) then
       ! make gradient totalsymmetric before adding
       do i_grad=1,grad_dim ! only if moving_c
          pointer_coul=>coul_int_grad_totsym(i_grad)%l(index)%m
          pointer_coul(:,k,i_ind,:,:)=&
               pointer_coul(:,k,i_ind,:,:)-&
               rotmat(i_grad,1)*help_arr_gr(:,4,:,:)-&
               rotmat(i_grad,2)*help_arr_gr(:,5,:,:)-&
               rotmat(i_grad,3)*help_arr_gr(:,6,:,:)
       enddo
    else
       pointer_coul=>coul_int_grad_totsym(1)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr(:,4,:,:)
       pointer_coul=>coul_int_grad_totsym(2)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr(:,5,:,:)
       pointer_coul=>coul_int_grad_totsym(3)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr(:,6,:,:)
    end if
    end if

    ! derivative with respect to second center b
    ! use that sum of derivatives is zero
    ! LOAD [ xi_i | f_k | d/dRb xi_j ]
    if (eff_moving_b) then
    do k_gr=1,3
       pointer_coul=>coul_int_grad(k_gr+3)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr(:,k_gr,:,:)+help_arr_gr(:,k_gr+3,:,:)
    end do
    end if
  end subroutine calc_coul_gradients

  subroutine add_grads(grad_mat_gr,compact_storage)
  real(kind=r8_kind) :: grad_mat_gr(:,:,:,:,:)
  logical, optional :: compact_storage
    ! LOAD dsym/dRa < xi_i | ... | xi_j >
    ! Input: grad_mat_gr(...,1:3) , Output: pointer_prim_int(ia:) or
    !                                       pointer_prim_int(1:N)
    if (moving_a) then
       if (present(compact_storage)) then
          index = 1
       else
          index = gradient_index(ima)
       endif
       grad_dim = gradient_index(ima+1) - gradient_index(ima)
    if(grad_dim==3) then
       if(sum((unique_atom_grad_info(ima)%m(:,:,1)-unity_matrix)**2)<&
            1.0e-7_r8_kind) then

          do_rotation=.false.
       else
          do_rotation=.true.
          rotmat=>unique_atom_grad_info(ima)%m(:,:,1)
       endif
    else
       do_rotation=.true.
       rotmat=>unique_atom_grad_info(ima)%m(:,:,1)
    end if
    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_a
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               rotmat(i_grad,1)*grad_mat_gr(:,:,:,:,1)+&
               rotmat(i_grad,2)*grad_mat_gr(:,:,:,:,2)+&
               rotmat(i_grad,3)*grad_mat_gr(:,:,:,:,3)
          index=index+1
       end do
    else
       pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
            grad_mat_gr(:,:,:,:,1)
       pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
            grad_mat_gr(:,:,:,:,2)
       pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
            grad_mat_gr(:,:,:,:,3)
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
       if(sum((unique_atom_grad_info(imb)%m(:,:,equalb)-unity_matrix)**2)<&
            1.0e-7_r8_kind) then
          do_rotation=.false.
       else
          do_rotation=.true.
          rotmat=>unique_atom_grad_info(imb)%m(:,:,equalb)
       endif
    else
       do_rotation=.true.
       rotmat=>unique_atom_grad_info(imb)%m(:,:,equalb)
    end if
    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_b
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               rotmat(i_grad,1)*grad_mat_gr(:,:,:,:,4)+&
               rotmat(i_grad,2)*grad_mat_gr(:,:,:,:,5)+&
               rotmat(i_grad,3)*grad_mat_gr(:,:,:,:,6)
          index=index+1
       end do
    else
       pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
            grad_mat_gr(:,:,:,:,4)
       pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
            grad_mat_gr(:,:,:,:,5)
       pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
            grad_mat_gr(:,:,:,:,6)
    endif
    endif
  end subroutine add_grads

  subroutine calculate_helpers
    implicit none
    ! purpose: calculate help quantity prod_arr_gr and 
    ! nuclear attraction gradients
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind) :: help_vec_arr(num,6),help_vec1(num),help_vec2(num)
    integer(kind=i4_kind) :: i,j,l_i
    integer(kind=i4_kind) :: i1,i2
    integer(kind=i4_kind), pointer :: index1_p(:),index2_p(:)
    real(kind=r8_kind) :: coeff
    real(kind=r8_kind), pointer :: coeff_p(:)
    integer(kind=i4_kind) :: k_gr,l,k,i_la,i_lb
    !------------ Executable code ------------------------------

    gamma_help(:,1:4+la)=gamma(4+la,gamma_arg2(:))
    prod_arr=0.0_r8_kind
    prod_arr_gr=0.0_r8_kind
    do i=1,(la+1)**2
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       coeff_p=>solhrules_product(i)%coef
       do j=1,solhrules_product(i)%n_summands
          i1=index1_p(j)
          l_i=solhrules_l_and_m_of_lm(1,i1)
          i2=index2_p(j)
          coeff=coeff_p(j)
          help_vec1=coeff*yl_arr(:,i1)
          prod_arr(:,i,l_i)=prod_arr(:,i,l_i)&
               +help_vec1*overlap(:,i2)
          help_vec2=coeff*overlap(:,i2)
          do k_gr=1,3
             prod_arr_gr(:,k_gr,i,l_i)=prod_arr_gr(:,k_gr,i,l_i)&
                  +help_vec1*overlap_grad(:,i2,k_gr)
             prod_arr_gr(:,k_gr+3,i,l_i)=prod_arr_gr(:,k_gr+3,i,l_i)&
                  +diff_arr_grad(:,i1,k_gr)*help_vec2
          end do
       end do
    enddo
    help_vec=(aexp_arr/fact0)
    do l=0,la
       do k=1,(la+1)**2
          do k_gr=1,3
             prod_arr_gr(:,k_gr,k,l)=prod_arr_gr(:,k_gr,k,l)+help_vec*&
                  prod_arr_gr(:,k_gr+3,k,l)
          end do
       end do
    end do
    do l=0,la
       do k=1,(la+1)**2
          do k_gr=1,3
             prod_arr_gr(:,k_gr,k,l+1) = prod_arr_gr(:,k_gr,k,l+1)+&
                  prod_arr(:,k,l)*(gamma_arg(:,k_gr)-xc(k_gr))
             prod_arr_gr(:,k_gr+3,k,l+1) = prod_arr_gr(:,k_gr+3,k,l+1)+&
                  prod_arr(:,k,l)*(gamma_arg(:,k_gr)-xc(k_gr))*rsaba
          end do
       end do
    end do
    ! now calculate nuclear attraction
    help_arr_gr=0.0_r8_kind
    do l=0,la+1
       help_vec=z*gamma_help(:,l+1)*&
            (-2.0_r8_kind*aexp_arr)**l*fact2
       do k_gr=1,6
          help_vec_arr(:,k_gr)=help_vec
       end do
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             l_index=la**2-1+i_la+i_lb 
             help_arr_gr(:,:,i_lb,i_la)=help_arr_gr(:,:,i_lb,i_la)&
                  +prod_arr_gr(:,:,l_index,l)*help_vec_arr
          enddo
       end do
    enddo
  end subroutine calculate_helpers


  subroutine calc_final_nuc_or_pvsp_grad(gradient_ab,gradient_c)
    implicit none
    real(kind=r8_kind) :: gradient_ab(num,2*la+1,6) !(:,:,:)
    ! derivatives with respect to centers a and b
    real(kind=r8_kind), optional :: gradient_c(num,2*la+1,grad_dim) !(:,:,:)
    logical :: eff_moving_a, eff_moving_b
    integer(kind=i4_kind) :: i_la,i_lb,k_gr

    eff_moving_a = (lagtlb .and. moving_a) .or. (laltlb .and. moving_b)
    eff_moving_b = (lagtlb .and. moving_b) .or. (laltlb .and. moving_a)
    ! derivative with respect to the centers a and b
    counter=1
    do i_la=1,nm_la
       do i_lb=1,nm_lb
          do k_gr=1,3
             ! derivative with respect to first center a
             if (eff_moving_a) then
             gradient_ab(:,counter,k_gr)=gradient_ab(:,counter,k_gr)+&
                  help_arr_gr(:,k_gr,i_lb,i_la)
             end if
             ! derivative with respect to second center b
             ! use that sum of derivatives is zero
             if (eff_moving_b) then
             gradient_ab(:,counter,k_gr+3)=gradient_ab(:,counter,k_gr+3)&
                  -help_arr_gr(:,k_gr,i_lb,i_la)+help_arr_gr(:,k_gr+3,i_lb,i_la)
             end if
          end do
          counter=counter+1
       end do
    end do

    if ( .not. present(gradient_c) .or. .not.moving_c) return

    ! derivative with respect to third center c
    if(do_rotation) then
       counter=1
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             ! make gradient totalsymmetric before adding
             do i_grad=1,grad_dim ! only if moving_c
                gradient_c(:,counter,i_grad)=gradient_c(:,counter,i_grad)-&
                     rotmat(i_grad,1)*help_arr_gr(:,4,i_lb,i_la)-&
                     rotmat(i_grad,2)*help_arr_gr(:,5,i_lb,i_la)-&
                     rotmat(i_grad,3)*help_arr_gr(:,6,i_lb,i_la)
             enddo
             counter=counter+1
          end do
       enddo
    else
       counter=1
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             gradient_c(:,counter,1)=gradient_c(:,counter,1)-&
                  help_arr_gr(:,4,i_lb,i_la)
             gradient_c(:,counter,2)=gradient_c(:,counter,2)-&
                  help_arr_gr(:,5,i_lb,i_la)
             gradient_c(:,counter,3)=gradient_c(:,counter,3)-&
                  help_arr_gr(:,6,i_lb,i_la)
             counter=counter+1
          end do
       enddo
    end if
  end subroutine calc_final_nuc_or_pvsp_grad



  subroutine s_coulomb(i)
    implicit none
    ! s-type coulomb fit integrals
    integer(kind=i4_kind), intent(in) :: i
    ! *** end of interface ***

    integer(kind=i4_kind) :: k,l,i_la,i_lb,k_gr
    real(kind=r8_kind)    :: help_vec_arr(num,6)

    ncexps = unique_atoms(i)%l_ch(0)%n_exponents
    cexps => unique_atoms(i)%l_ch(0)%exponents(:)
    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       ! precalculation of two factors
       rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
       fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))
       gamma_help(:,1:2+la)=gamma(2+la,gamma_arg2(:)*rcsabc(:))
       help_arr_gr=0.0_r8_kind
       do l=0,la+1
          help_vec=gamma_help(:,l+1)*&
               fact8*(-2.0_r8_kind*aexp_arr*rcsabc(:))**l
          do k_gr=1,6
             help_vec_arr(:,k_gr)=help_vec
          end do
          do i_la=1,nm_la
             do i_lb=1,nm_lb
                help_arr_gr(:,:,i_lb,i_la)=help_arr_gr(:,:,i_lb,i_la)+&
                     prod_arr_gr(:,:,la**2-1+i_la+i_lb,l)*help_vec_arr
             end do
          end do
       enddo! loop over l
       call calc_coul_gradients(0,k,1)
    end do! loop over fitfunctions
  end subroutine s_coulomb


  subroutine r2_coulomb(i)
    implicit none
    ! r2-type coloumb fit integral
    integer(kind=i4_kind), intent(in) :: i
    ! *** end of interface ***

    integer(kind=i4_kind) :: k,l,i_la,i_lb,k_gr
    real(kind=r8_kind) :: help_vec_arr(num,6)

    ncexps = unique_atoms(i)%r2_ch%n_exponents
    cexps => unique_atoms(i)%r2_ch%exponents(:)
    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       ! precalculation of two factors
       rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
       fact8=2.0_r8_kind*pi/(cexps(k)**2)*sqrt(fact0/(fact0+cexps(k)))&
            *(fact0/(fact0+cexps(k)))
       fact7=(fact0+1.5_r8_kind*cexps(k))/fact0
       gamma_help(:,1:3+la)=gamma(3+la,gamma_arg2(:)*rcsabc(:))
       help_arr_gr=0.0_r8_kind
       do l=0,la+1  
          help_vec=fact8*(-2.0_r8_kind*aexp_arr*rcsabc(:))**l*&
               ((fact7-real(l,r8_kind))*gamma_help(:,l+1)+&
               gamma_arg2(:)*rcsabc(:)*gamma_help(:,l+2))
          do k_gr=1,6
             help_vec_arr(:,k_gr)=help_vec
          end do
          do i_la=1,nm_la
             do i_lb=1,nm_lb
                help_arr_gr(:,:,i_lb,i_la)=help_arr_gr(:,:,i_lb,i_la)+&
                     prod_arr_gr(:,:,la**2-1+i_la+i_lb,l)*help_vec_arr
             end do
          end do
       enddo! loop over l
       call calc_coul_gradients(-1,k,1)

    end do! loop over fitfunctions
  end subroutine r2_coulomb


  subroutine l_coulomb(i,j)
    implicit none
    ! l-type coulomb fit integrals
    integer(kind=i4_kind), intent(in) :: i,j
    ! *** end of interface ***

!!$    integer(kind=i4_kind) :: k,l,i_la,i_lb,k_gr
    integer(kind=i4_kind) :: i_l,i_ind,i_cnt,k,l,k_gr,i_la,i_lb
    real(kind=r8_kind) :: help_vec_arr(num,6)

    do i_l=1,lmax_ch
       n_independent_fcts  = &
            unique_atoms(i)%symadapt_partner(1,i_l)%n_independent_fcts
       allocate(&
            intermediate_gr(num,6,2*la+1,n_independent_fcts,0:la+1),&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler('LS_CALCULATE_GRADS: allocation intermediate_gr failed')
       intermediate_gr=0.0_r8_kind
       !    
       ! do symmetry adaption
       independents: do i_ind=1,n_independent_fcts
          n_contributing_fcts = &
               unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%n_fcts
          coeff => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c
          magn => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m
          eq_atom => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%&
               I_equal_atom 
          contributing: do i_cnt=1,n_contributing_fcts
             if(eq_atom(i_cnt)==j) call calc_intermediate(i_l,i_ind,i_cnt)
          end do contributing
       end do independents
       ! end of symmetry adaption

       ncexps=  unique_atoms(i)%l_ch(i_l)%n_exponents
       cexps => unique_atoms(i)%l_ch(i_l)%exponents(:)

       do k=1,ncexps ! loop over fitexponents
          rcsabc=cexps(k)/(fact0+cexps(k))
          gamma_help(:,1:max_order)=&
               gamma(max_order,gamma_arg2(:)*rcsabc(:))
          fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))
          counter=1
          do i_ind=1,n_independent_fcts
             help_arr_gr=0.0_r8_kind
             do l=0,la+1
                help_vec=(2.0_r8_kind*fact0*rcsabc)**i_l*gamma_help(:,l+1+i_l)*&
                     (-2.0_r8_kind*aexp_arr*rcsabc(:))**l*fact8
                do k_gr=1,6
                   help_vec_arr(:,k_gr)=help_vec
                end do
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      help_arr_gr(:,:,i_lb,i_la)=help_arr_gr(:,:,i_lb,i_la)+&
                           intermediate_gr(:,:,i_la+i_lb-1,i_ind,l)*help_vec_arr
                   end do
                end do
             end do! loop over l
             call calc_coul_gradients(i_l,k,i_ind) 
          end do! loop over i_ind
       enddo! loop over k
       deallocate(intermediate_gr,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler &
            ("LS_CALCULATE_GRADS: deallocation intermediate_gr failed")
    end do! loop over lc
  end subroutine l_coulomb

  subroutine calc_intermediate(i_l,i_ind,i_cnt)
    implicit none
    ! Purpose: piece of the code
    !          calculate intermediate_gr, which is necesarry for
    !          three-center l-type coulomb integrals
    integer(kind=i4_kind), intent(in) :: i_l,i_ind,i_cnt
    ! *** end of interface ***

    integer(kind=i4_kind) :: i_p, i_p2, j_p
    integer(kind=i4_kind), pointer :: index1_p(:), index2_p(:)
    real(kind=r8_kind), pointer :: coef_p(:)
    real(kind=r8_kind) :: help_vec(num), help_vec2(num), help_vec2_arr(num,6), &
         help_vec3(num)
    integer(kind=i4_kind) :: k_gr,l

    help_vec=(aexp_arr/fact0)**i_l
    do k_gr=1,3
       diff_arr_gr(:,:,k_gr)=diff_rule(diff_arr_grad(:,:,k_gr)/fact10&
            ,1,(la+1)**2,i_l**2+magn(i_cnt))
       do l=1,(la+1)**2
          diff_arr_gr(:,l,k_gr)=diff_arr_gr(:,l,k_gr)*help_vec
       enddo
    end do
    diff_arr(:,:)=spread(help_vec,2,(la+1)**2)*&
         diff_rule(yl_arr(:,:)/&
         fact10(:,1:((ly_max+1)**2)),1,(la+1)**2,i_l**2+&
         magn(i_cnt))
    ! now apply product rule
    help_vec=aexp_arr/fact0
    do i_p=1,2*la+1
       i_p2=i_p+la*la
       index1_p=>solhrules_product(i_p2)%lm_sh1
       index2_p=>solhrules_product(i_p2)%lm_sh2
       coef_p=>solhrules_product(i_p2)%coef     
       do j_p=1,solhrules_product(i_p2)%n_summands
          help_vec2=diff_arr(:,index2_p(j_p))*coeff(i_cnt)*coef_p(j_p)
          do k_gr=1,6
             help_vec2_arr(:,k_gr)=help_vec2
          end do
          do l=0,la+1
             intermediate_gr(:,:,i_p,i_ind,l)=&
                  intermediate_gr(:,:,i_p,i_ind,l)+&
                  prod_arr_gr(:,:,index1_p(j_p),l)*help_vec2_arr
          end do
          do l=0,la
             help_vec2=coeff(i_cnt)*coef_p(j_p)*prod_arr(:,index1_p(j_p),l)
             help_vec3=help_vec2*help_vec
             do k_gr=1,3
                intermediate_gr(:,k_gr,i_p,i_ind,l)=&
                     intermediate_gr(:,k_gr,i_p,i_ind,l)+&
                     help_vec3*diff_arr_gr(:,index2_p(j_p),k_gr)
                intermediate_gr(:,k_gr+3,i_p,i_ind,l)=&
                     intermediate_gr(:,k_gr+3,i_p,i_ind,l)+&
                     help_vec2*diff_arr_gr(:,index2_p(j_p),k_gr)
             end do
          end do
       end do! 
    end do! i_p
  end subroutine calc_intermediate

end module ls_calculate_grads_module

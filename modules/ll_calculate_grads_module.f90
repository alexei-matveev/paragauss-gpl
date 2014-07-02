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
module ll_calculate_grads_module
  !
  !  Purpose: calculation of the gradints of all primitive 2 center orbital
  !           and 3 center integrals for a given set of indizes
  !       (unique_atom1,unique_atom2,la,lb).
  !       For three center integrals, contraction and symmetry-
  !       adaption concerning fitfunctions is also performed.
  !
  !
  !  Author: MS
  !  Date:  1/97
  !
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
  !-------------------------------------------------------------------

#include "def.h"
  use unique_atom_module, noname=>pseudopot_present
  use gamma_module
  use operations_module, only: operations_core_density
  use type_module
  use datatype
  use solid_harmonics_module, only: solid_harmonics_calc, &
       solid_harmonics_scalar
  use int_data_2cob3c_module, only: center1, center2,  &
       prim_int_2cob_ol_grad, prim_int_3cob_grad, prim_int_2cob_nuc_grad, &
       prim_int_2cob_ks_grad, prim_int_3cob_nuc_grad, prim_int_3cob_coul_grad, &
       prim_int_2cob_kin_grad, prim_int_2cob_pvsp_grad, prim_int_3cob_epe, &
       prim_int_3cob_solv_grad , prim_int_2cob_pseudo_grad, &
#ifndef no_cpks_coul_grads
       prim_int_cpks_coul_grad, &
#endif
       prim_int_3cob_solv_grad_pc, & !!!!!!!!!!!
       prim_int_coul_dervs,prim_int_2cob_ol_dervs
  use solhrules_module
  use fitcontract_module
  use integralpar_module
  use gradient_data_module, only: gradient_index, gradient_data_n_gradients, &
                                  gradient_data_n_spin_gradients,calc_cluster_epe_energy
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
  USE_MEMLOG

  implicit none

  !== Interrupt end of public interface of module ====================

  public :: l_coulomb
  public :: unity_matrix_p!(matrix) -> logical

!  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
!  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
!  integer(kind=i4_kind),intent(in) :: equalb ! number of equal atom b
!  integer(kind=i4_kind),intent(in) :: la ! angular momentum of unique atom a
!  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b

  ! This comment is not quite true, sadly everything in this module is
  ! public yet:

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
! overlap_xyz   (:,:,:,   1:3  ) :    d/dRa < xi_i | 1      | xi_j >
! kinetic_xyz   (  :  ,   1:3  ) :    d/dRa < xi_i | T      | xi_j >
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
! grad_mat     (:,:,:,:,  1:3  ) :    d/dRa < xi_i | ...    | xi_j >
! grad_mat     (:,:,:,:,  4:6  ) :    d/dRb < xi_i | ...    | xi_j >
! grad_mat_spin(:,:,:,:,  1:3  ) :    d/dRa < xi_i | ...    | xi_j >
! grad_mat_spin(:,:,:,:,  4:6  ) :    d/dRb < xi_i | ...    | xi_j >
! grad_mat_pvsp(:,:,:,:,  1:3  ) :    d/dRa < xi_i | ...    | xi_j >
! grad_mat_pvsp(:,:,:,:,  4:6  ) :    d/dRb < xi_i | ...    | xi_j >
!
!..............................................................................


  integer(kind=i4_kind) :: naexps,nbexps,ncexps
  real(kind=r8_kind),pointer :: aexps(:),bexps(:)
  real(kind=r8_kind),pointer     :: cexps(:)
  real(kind=r8_kind) :: z  ! charge
  real(kind=r8_kind) :: zc ! core charge
  integer(kind=i4_kind)          :: max_order,grad_dim,index,i_grad,&
       nm_la,nm_lb,la_index,lb_index,i_la,lasq,lbsq,laposq,lbposq,k_gr,k2dr
  integer(kind=i4_kind)          :: ima, imb, imc, k_gr_0, k_gr_1, &
                                    spin_index
  logical                        :: moving_a, moving_b, moving_c, &
                                    split_gradients, model_density, &
                                    spin_polarized
  ! constants
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
#if defined(OLD_LINUX)
  real(kind=r8_kind),dimension(0:8),parameter :: dfac=(/1.0_r8_kind,&
       1.0_r8_kind,3.0_r8_kind,15.0_r8_kind,105.0_r8_kind,945.0_r8_kind,&
       10395.0_r8_kind,135135.0_r8_kind,2027025.0_r8_kind/)
#else
  real(kind=r8_kind),dimension(0:8) :: dfac=(/1.0_r8_kind,&
       1.0_r8_kind,3.0_r8_kind,15.0_r8_kind,105.0_r8_kind,945.0_r8_kind,&
       10395.0_r8_kind,135135.0_r8_kind,2027025.0_r8_kind/)
#endif
  integer(kind=i4_kind),dimension(3) :: xyz_map
  ! array help
  logical,allocatable   :: cutoff(:,:)
  logical :: do_rotation ! flag to decide if gradients have to be rotated in order
  logical :: iam_ppot
  ! to make it totalsymmetric
  integer(kind=i4_kind) :: num,counter,m,ma,mb,lmax,l2pm
  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: fact0_arr, &
       fact1_arr,fact2_arr,fact10
  real(kind=r8_kind),allocatable,dimension(:)  :: fact0,fact1, &
       fact2,fact4,fact5,fact6,fact7,fact8,rcsabc,tau,help_vec,&
       help_vec0
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: help_mat
  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable     :: gamma_arg(:,:),gamma_arg2(:)
  real(kind=r8_kind),allocatable,dimension(:,:)     :: gamma_help


  ! help arrays for solid harmincs
  real(kind=r8_kind),allocatable     :: yl_arr(:,:),yl_arr_xyz(:,:,:),&
       yl_arr2(:,:),aqapb(:,:)
  real(kind=r8_kind),allocatable     :: &
       clmamb(:,:),clmamb_xyz(:,:,:),&
       clmamb2(:,:),clmamb2_xyz(:,:,:), &
       clmamb_scalar(:),clmamb_scalar_xyz(:,:)

  ! help arrays for product_rule and diff_rule
  real(kind=r8_kind),allocatable  :: &
       prod_arr(:,:,:,:),&
       prod_arr_gr(:,:,:,:,:),&
       prod_arr_gr_vec(:,:,:),&
       prod_arr_rel_gr(:,:,:,:,:),&
       diff_arr(:,:,:),&
       diff_arr_xyz(:,:,:,:),&
       diff_arr0(:,:),&
       diff_arr0_xyz(:,:,:), &
       intermediate(:,:,:,:,:),&
       intermediate_gr(:,:,:,:),&
       help_arr(:,:),&
       help_arr0(:,:),&
       help_arr2(:,:,:),&
       help_arr_2(:),&
       rk_bc_arg_help(:,:),&
       rk_ac_arg_help(:,:)

  real(kind=r8_kind) :: arg, help_real_0, help_real, help_real_2,help_real_3

  ! cartesian coordinates
  real(kind=r8_kind),dimension(3)  :: xa,xb,xc,xd
  ! exchange necessary or not

  integer(kind=i4_kind)  :: alloc_stat! (40)=0
  integer(kind=i4_kind)  :: i,j,i_l,j_l,i_lb,k,i_ind,i_cnt,l, &
       la2pm,lb2pm,i_i,i_ma,i_aexps,i_bexps
  integer(kind=i4_kind)  :: lmax_ch,lmax_xc,lmax_abs,ly_max
  integer(kind=i4_kind)  :: n_equals,n_independent_fcts, &
       n_contributing_fcts

  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coeff(:),rotmat(:,:)

  real(kind=r8_kind),pointer :: grad_mat_p(:,:,:,:)
  real(kind=r8_kind),allocatable :: grad_mat_sp(:,:,:,:)
  real(kind=r8_kind),allocatable :: grad_mat(:,:,:,:,:), &
        grad_mat_spin(:,:,:,:,:) !, grad_mat_pvsp(:,:,:,:,:)

  real(kind=r8_kind),pointer :: dervs_mat(:,:,:,:,:,:),ca_dervs_mat(:,:,:,:,:,:)

  real(kind=r8_kind),allocatable    :: &
       overlap(:,:,:), overlap_xyz(:,:,:,:), &
       aexp_arr(:),bexp_arr(:),twoa(:),twob(:)

  real(kind=r8_kind),allocatable    :: nested2_fac1(:,:),nested2_fac2(:,:)
  real(kind=r8_kind),allocatable    :: nested2_fac12(:,:,:)
  real(kind=r8_kind),allocatable    :: &
       nuc_grad(:,:,:,:),& ! nuclear gradients with respect to centers a and b
       nuc_grad_gr(:,:,:,:), & ! nuclear gradients with respect to center c
       help_arr_gr(:,:,:,:),&
       help_arr_gr_vec(:,:,:), &
       grad_nuc_pseudo_ab(:,:,:,:),&
       grad_nuc_pseudo_c(:,:,:,:)
  real(kind=r8_kind),allocatable    :: pseudo_grad_gr(:,:,:,:) ! arr for a &b components
  real(kind=r8_kind),allocatable  :: cluster_epe(:,:,:), &
       solv_grad(:,:,:,:), &   !!!!!!!!!!!!!
       solv_grad_gr(:,:,:,:), &   !!!!!!!!!!!!!
       help_arr_gr1(:,:,:,:)


#ifndef FPP_NOSAVE
  ! gfortran crashes on accessing components without "save" attribute:
  ! ( the vars with =>NULL() initialization of pointer components need it )
  type(three_center_l),dimension(6)  ,save, target   :: coul_int_grad
  type(three_center_l),dimension(6,6),save   :: coul_int_dervs
#else
  ! another fortran in Krsk crashes on allocation with "save" attribute:
  ! ( add FPPOPTIONS += -DFPP_NOSAVE to your machine.inc file )
  type(three_center_l),dimension(6)          :: coul_int_grad
  type(three_center_l),dimension(6,6)        :: coul_int_dervs
#endif

  type(three_center_l),allocatable, target  :: coul_int_grad_totsym(:)
  type(three_center_l),allocatable  :: coul_int_dervs_totsym(:,:)
  type(three_center_l),allocatable  :: coul_int_ca_dervs(:,:)

  logical                           :: pseudopot_present ! same name as in UA module

  ! totalsymmetric gradient

  type(arrmat4),pointer :: pointer_prim_int(:)     ! help pointer to point on primitive integrals

  type(arrmat5),pointer :: pointer_prim_3c(:) ! help pointer to point on primitive integrals
  real(kind=r8_kind),pointer,dimension(:,:,:,:,:) :: &
       pointer_s, pointer_r2, pointer_l

  type nested2_opt
   integer(kind=i4_kind):: n
   type(nested2_vars),pointer,dimension(:)::summands
  end type nested2_opt
  type(nested2_opt),pointer,dimension(:):: nested2_summands,nested2ne_summands
  integer(kind=i4_kind):: nested2_l1_max,nested2_l3_max,l1_max,l3_max,nne_summands
  logical:: nested2_optimized = .false.

  intrinsic max

contains

  subroutine calc_solv_grad(la,lb) !!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=i4_kind),intent(in):: la,lb
    !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    integer(kind=i4_kind) :: n_equal_solv,k1,m1
    !solv_grad_gr : orbital derivatives
    !solv_grad    : symmetrized operator derivaties

    fact8=2.0_r8_kind*sqrt(fact0/pi)

       ly_max=max(la,lb)
       max_order=4+la+lb
       allocate ( gamma_help(num,max_order), &
            gamma_arg2(num),&
            yl_arr(num,(ly_max+1)**2),&
            yl_arr_xyz(num,(ly_max+1)**2,3),&
            help_arr(num,(ly_max+1)**2),&
            help_arr0(num,laposq),&
            fact10(num,(ly_max+1)**2),&
            yl_arr2(num,(ly_max+1)**2),&
            aqapb(num,0:ly_max),&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LL_CALCULATE_GRADS : allocation 8 failed(solv)")
       aqapb(:,0)=1.0_r8_kind
       aqapb(:,1)=aexp_arr/fact0
       do i_l=2,ly_max
          aqapb(:,i_l)=aqapb(:,i_l-1)*aqapb(:,1)
       end do



       counter=1
       fact4=1.0_r8_kind
       do i_l=0,ly_max
          do ma=1,2*i_l+1
             fact10(:,counter)=fact4
             counter=counter+1
          enddo
          fact4=fact4*aqapb(:,1)
       enddo

       solv_grad_gr=0.0_r8_kind ! to provide summurizing of solv_grad_gr only one time
                                   ! over all surface point charges
       help_arr_gr1=0.0_r8_kind
       lab_k1: do k1=1,to_calc_grads%n_points
             call integral_interrupt_2cob3c()
             n_equal_solv=to_calc_grads%n_equal(k1)
             z=to_calc_grads%Q(k1)

             lab_m1: do m1=1,n_equal_solv
                yl_arr_xyz=0.0_r8_kind
                xc =to_calc_grads%xyz(k1,m1,:)
                yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                     spread(xc,1,num))
                gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                     (gamma_arg(:,3)-xc(3))**2)*fact0
                gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2(:))
                fact8=2.0_r8_kind*sqrt(fact0/pi)
                counter=1
                yl_arr2(:,:)=yl_arr(:,:)/fact10(:,1:((ly_max+1)**2))

                do l=0,ly_max
                   do m=1,2*l+1
                      l2pm=l*l+m
                      do k_gr=1,3
                         do k=1,solhrules_differential(xyz_map(k_gr),l2pm)%n_summands
                            yl_arr_xyz(:,l2pm,k_gr)=yl_arr_xyz(:,l2pm,k_gr)+&
                                 solhrules_differential(xyz_map(k_gr),l2pm)%coef(k)*&
                                 yl_arr(:,solhrules_differential(xyz_map(k_gr),l2pm)%lm_sh(k))
                         end do
                      end do
                   end do
                end do!loop over l

                do i_l=0,lb
                   help_arr0=spread(aqapb(:,i_l),2,laposq)
                   do mb=1,2*i_l+1
                      diff_arr(:,:,counter)=&
                           help_arr0*diff_rule(yl_arr2(:,:)&
                           ,1,laposq,counter)
                      do k_gr=1,3
                         diff_arr_xyz(:,:,counter,k_gr)=&
                              help_arr0*diff_rule(yl_arr_xyz(:,:,k_gr)/fact10&
                              ,1,laposq,counter)
                      end do
                      counter=counter+1
                   end do
                end do
                do mb=1,lbposq
                   do ma=1,laposq
                      help_mat(:,ma,mb)= overlap(:,ma,mb)*aqapb(:,1)
                   enddo
                enddo
                help_vec=-2.0_r8_kind*fact0
            if(nested2_optimized) then
             call calc_prod_arr_opt(la,lb,nested2_l1_max,nested2_l3_max, &
                       nested2_fac12)
            else
             call calc_prod_arr(la,lb)   !in  calc_solv_grad
            endif

                call calc_solv_gradients(la,lb,k1,m1)
             enddo lab_m1
       enddo lab_k1


       deallocate(yl_arr,yl_arr_xyz,yl_arr2,help_arr,&
            help_arr0,aqapb,gamma_help,gamma_arg2,fact10,&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LL_CALCULATE_GRADS : deallocation 17 failed(solv)")



    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  end subroutine calc_solv_grad    !!!!!!!!!!!!!!!!!!!!!!!!
  !**************************************************************


  subroutine calc_nuc_help_arr(la,lb)
  integer(kind=i4_kind),intent(in):: la,lb
    help_arr=0.0_r8_kind
    help_arr_gr=0.0_r8_kind
    do i_l =0,lb+la+1
       help_vec=fact8*gamma_help(:,i_l+1)*z
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             la_index=lasq+i_la
             lb_index=lbsq+i_lb
             do k_gr=1,6
                help_arr_gr(:,k_gr,i_lb,i_la)=help_arr_gr(:,k_gr,i_lb,i_la)&
                     +prod_arr_gr(:,k_gr,la_index,lb_index,i_l)*help_vec
             end do! loop over k_gr
          enddo! loop over i_lb
       end do! loop over i_la
    end do! loop over i_l
  end subroutine calc_nuc_help_arr

  subroutine calc_2center(la,lb,equalb)

  integer(kind=i4_kind),intent(in):: equalb
  integer(kind=i4_kind),intent(in):: la,lb


    ! Purpose:  calculate the gradients of overlap and kin

    ! first calculating 2-center integrals----------------
    ! fact5=fact2*(3.0_r8_kind-2.0_r8_kind*tau+2.0_r8_kind*la)
    ! a*b/(a+b)(3-2*tau+2*l)

    fact6=1.0_r8_kind/sqrt(aexp_arr**la*dfac(la))/sqrt(bexp_arr**lb*dfac(lb)) &
         *exp(-tau)*(4.0_r8_kind*fact2/fact0)**0.75_r8_kind

    fact5=fact2*fact6
    fact7=(fact2*2.0_r8_kind)**lb

    counter=1
    do i_l=0,lb
       do mb=1,2*i_l+1

          diff_arr0(:,counter)=reshape(diff_rule( &
                   spread(clmamb_scalar,1,1),1,laposq,counter),(/laposq/))

          ! now calculate derivatives of diff_arr0
          do k_gr=1,3
             diff_arr0_xyz(:,counter,k_gr)=reshape(diff_rule(&
                  spread(clmamb_scalar_xyz(:,k_gr),1,1),1,laposq,counter),(/laposq/))
                        !              ^
          enddo

          counter=counter+1
       end do
    end do


    diff_arr_xyz=spread(diff_arr0_xyz(:,:,:),1,num)

    diff_arr=spread(diff_arr0(:,:),1,num)

    counter=1
    do i_l=0,lb
       help_mat(:,:,1)=spread(fact6*(2.0_r8_kind*fact2)**i_l,2,laposq)
       magnetic_number_b: do mb=1,2*i_l+1
          ! overlap
          overlap(:,1:laposq,counter)=help_mat(:,:,1)*&
            prod_rule(spread(diff_arr0(:,counter),1,num),clmamb(:,:),1,laposq)
            ! d/da help_mat = help_mat*(-2.0_r8_kind*fact2)*xd(k_gr)

          ! now calculate derivatives of overlap

          do k_gr=1,3
             !                                    (d/da help_mat) * prod_rule(scal,scal)
             overlap_xyz(:,1:laposq,counter,k_gr)=overlap(:,1:laposq,counter)* &
                       (-2.0_r8_kind)*spread(fact2,2,laposq)*xd(k_gr)
             !      + help_mat * (d/da prod_rule)
             overlap_xyz(:,1:laposq,counter,k_gr)= overlap_xyz(:,1:laposq,counter,k_gr)&
                 +help_mat(:,:,1)*(prod_rule(diff_arr_xyz(:,:,counter,k_gr),clmamb(:,:),1,laposq) &
                                  +prod_rule(diff_arr(:,:,counter),clmamb_xyz(:,:,k_gr),1,laposq))
          end do

          counter=counter+1
       end do magnetic_number_b
    enddo

    call integral_interrupt_2cob3c()

    MEMLOG(-size(fact1)-size(tau))
    MEMLOG(-size(clmamb_scalar)*4-size(clmamb2)*4)
    deallocate(fact1,clmamb_scalar,clmamb_scalar_xyz,&
         clmamb2,clmamb2_xyz,tau,stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)

  end subroutine calc_2center


  subroutine calc_coul_gradients(index,k,i_ind)
    ! Purpose: Build final gradients of the coulomb integrals
    integer(kind=i4_kind) :: index ! angular momentum
    integer(kind=i4_kind) :: k ! number of fitfunction
    integer(kind=i4_kind) :: i_ind ! number of independent function
    !** End of interface ****************************************
    real(kind=r8_kind),pointer :: pointer_coul(:,:,:,:,:)
    ! derivative with respect to first center a
    ! LOAD [ d/dRa xi_i | f_k | xi_j ]
    if (moving_a) then
    do k_gr=1,3
       pointer_coul=>coul_int_grad(k_gr)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)+ &
            help_arr_gr_vec((k_gr-1)*num+1:k_gr*num,:,:)
    end do
    endif

    ! derivative with respect to third center c
    ! LOAD [ xi_i | dsym/dRc f_k | xi_j ]
    if (moving_c) then
    if(do_rotation) then
       ! make gradient totalsymmetric before adding
       do i_grad=1,grad_dim ! only if moving_c
          pointer_coul=>coul_int_grad_totsym(i_grad)%l(index)%m
          pointer_coul(:,k,i_ind,:,:) = &
               pointer_coul(:,k,i_ind,:,:) - &
               rotmat(i_grad,1)*help_arr_gr_vec(3*num+1:4*num,:,:) - &
               rotmat(i_grad,2)*help_arr_gr_vec(4*num+1:5*num,:,:) - &
               rotmat(i_grad,3)*help_arr_gr_vec(5*num+1:6*num,:,:)
       enddo
    else
       pointer_coul=>coul_int_grad_totsym(1)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr_vec(3*num+1:4*num,:,:)
       pointer_coul=>coul_int_grad_totsym(2)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr_vec(4*num+1:5*num,:,:)
       pointer_coul=>coul_int_grad_totsym(3)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr_vec(5*num+1:6*num,:,:)
    end if
    end if

    ! derivative with respect to second center b
    ! use that sum of derivatives is zero
    ! LOAD [ xi_i | f_k | d/dRb xi_j ]
    if (moving_b) then
    do k_gr=1,3
       pointer_coul=>coul_int_grad(k_gr+3)%l(index)%m
       pointer_coul(:,k,i_ind,:,:)=&
            pointer_coul(:,k,i_ind,:,:)-&
            help_arr_gr_vec((k_gr-1)*num+1:k_gr*num,:,:) + &
            help_arr_gr_vec((k_gr+2)*num+1:(k_gr+3)*num,:,:)
    end do
    end if
  end subroutine calc_coul_gradients
  !**************************************************************

  !**************************************************************
  subroutine calc_nuc_gradients(la,lb,third_center_required)
    ! Purpose: Calculate final gradient of nuclear attraction
    integer(kind=i4_kind),intent(in):: la,lb
    logical, intent(in) :: third_center_required
    !** End of interface ****************************************
    help_arr=0.0_r8_kind
    help_arr_gr=0.0_r8_kind
    do i_l =0,lb+la+1
       help_vec=fact8*gamma_help(:,i_l+1)*z
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             la_index=lasq+i_la
             lb_index=lbsq+i_lb
             do k_gr=1,6
                help_arr_gr(:,k_gr,i_lb,i_la)=help_arr_gr(:,k_gr,i_lb,i_la)&
                     +prod_arr_gr(:,k_gr,la_index,lb_index,i_l)*help_vec
             end do! loop over k_gr
          enddo! loop over i_lb
       end do! loop over i_la
    end do! loop over i_l
    if ( third_center_required ) then
       if(iam_ppot) then
          call calc_final_nuc_and_pvsp_grad(pseudo_grad_gr,nuc_grad)
       else
          call calc_final_nuc_and_pvsp_grad(nuc_grad_gr,nuc_grad)
       end if
    else
       if(iam_ppot) then
          call calc_final_nuc_and_pvsp_grad(pseudo_grad_gr)
       else
          call calc_final_nuc_and_pvsp_grad(nuc_grad_gr)
       end if
    endif
  end subroutine calc_nuc_gradients
  !**************************************************************

  !**************************************************************
  subroutine calc_solv_gradients(la,lb,kk,mm)
  integer(kind=i4_kind),intent(in):: la,lb

    integer(kind=i4_kind), intent(in) :: kk,mm
    integer(kind=i4_kind) :: nn,n1,ism,N_pc
    integer(kind=i4_kind) :: i,nequ,i_la,i_lb,km
!   integer(kind=i4_kind) :: ieq

    N_pc=0
    if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

    help_arr=0.0_r8_kind
    help_arr_gr=0.0_r8_kind
    do i_l =0,lb+la+1
       help_vec=fact8*gamma_help(:,i_l+1)*z
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             la_index=lasq+i_la
             lb_index=lbsq+i_lb
             do k_gr=1,6
                help_arr_gr(:,k_gr,i_lb,i_la)=help_arr_gr(:,k_gr,i_lb,i_la)&
                     +prod_arr_gr(:,k_gr,la_index,lb_index,i_l)*help_vec
             end do! loop over k_gr
          enddo! loop over i_lb
       end do! loop over i_la
    end do! loop over i_l
    call calc_final_nuc_and_pvsp_grad(solv_grad_gr) !only solv_grad_gr !!!!

    ism=to_calc_grads%i_symm_sort(kk,mm)


     do i=1,N_moving_unique_atoms+N_pc
        if( i <= N_moving_unique_atoms) then
           index = gradient_index(i) - 1
           grad_dim=gradient_index(i+1)-gradient_index(i)
           km = moving_unique_atom_index(i)
           nequ=unique_atoms(km)%n_equal_atoms
        else
           km=i-N_moving_unique_atoms
           nequ=pointcharge_array(km)%N_equal_charges
           index = surf_points_grad_index(km) - 1
           grad_dim=surf_points_grad_index(km+1)-surf_points_grad_index(km)
        end if

!      do ieq=1,nequ

         help_arr_gr1=0.0_r8_kind

         do nn=1,grad_dim
           do n1=4,6
            help_arr_gr1(:,nn,:,:)=help_arr_gr1(:,nn,:,:)+ &
                    help_arr_gr(:,n1,:,:)*&
                    to_calc_grads%dxyz_totsyms(nn,i)%m(n1-3,ism)
!                   to_calc_grads%dxyz_totsym(nn,i,ieq)%m(n1-3,ism)
           enddo
         enddo

!      enddo !ieq

       if(i <= N_moving_unique_atoms) then
          do i_grad=1,grad_dim
             grad_mat_p=>prim_int_3cob_solv_grad(index+i_grad)%m
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   grad_mat_p(:,:,i_lb,i_la)=grad_mat_p(:,:,i_lb,i_la)&   !!!!
                        +unpack(help_arr_gr1(:,i_grad,i_lb,i_la),cutoff,zero) !!!!
                end do
             end do
          end do
       else
          do i_grad=1,grad_dim
             grad_mat_p=>prim_int_3cob_solv_grad_pc(index+i_grad)%m
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   grad_mat_p(:,:,i_lb,i_la)=grad_mat_p(:,:,i_lb,i_la)&   !!!!
                        -unpack(help_arr_gr1(:,i_grad,i_lb,i_la),cutoff,zero) !!!!
                end do
             end do
          end do
       end if


  enddo
!       if(kk.eq.1) print*, sum(help_arr_gr1),'help_arr_grs',sum(help_arr_gr)

  end subroutine calc_solv_gradients
  !**************************************************************

  !**************************************************************

  subroutine calc_prod_arr(la,lb)
  integer(kind=i4_kind),intent(in):: la,lb
    ! Purpose: piece of the code
    !          calculate quantity prod_arr
    real(kind=r8_kind) :: help_arr(3*num,(ly_max+1)**2),help_arr2(3*num,0:la+lb,laposq),&
         coef
    integer(kind=i4_kind) :: l_sum,i3_p_j3,i3_p2_j3,i_p_j1,i2_p_j2,i0_p_j1,meta_min,meta_max
    real(kind=r8_kind),dimension(num) :: help_vec1,help_vec0,help_vec2,help_vec3
    integer(kind=i4_kind) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j3, l_j1
    integer(kind=i4_kind),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(kind=r8_kind),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    real(kind=r8_kind) :: &
         overlap_xyz_help(num*3), &
         overlap_help(num*6), &
         help_spread(num*6)

    prod_arr=0.0_r8_kind
    prod_arr_gr_vec=0.0_r8_kind
    do l=1,(ly_max+1)**2
       do k_gr=1,3
          help_arr((k_gr-1)*num+1:k_gr*num,l)=yl_arr(:,l)*(gamma_arg(:,k_gr)-xc(k_gr))
       end do
    end do
    help_vec=-2.0_r8_kind*fact0
    do i_l=0,lb
       do mb=1,2*i_l+1
          help_arr2=0.0_r8_kind
          lb_index=i_l**2+mb
          lm_in=lb_index
          lm_min=1
          lm_max=(la+1)**2
          index_p=>solhrules_product(lm_in)%lm_sh1
          index0_p=>solhrules_product(lm_in)%lm_sh2
          coef1_p=>solhrules_product(lm_in)%coef
          n_sum1=solhrules_product(lm_in)%n_summands
          do i=lm_min,lm_max
             index1_p=>solhrules_product(i)%lm_sh1
             index2_p=>solhrules_product(i)%lm_sh2
             n_sum2=solhrules_product(i)%n_summands
             coef2_p=>solhrules_product(i)%coef
             do j_1=1,n_sum1
                i_p_j1=index_p(j_1)
                i0_p_j1=index0_p(j_1)
                l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                help_vec2=twob**l1
                do j_2=1,n_sum2
                   i2_p_j2=index2_p(j_2)
                   index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
                   index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
                   n_sum3=solhrules_product(index1_p(j_2))%n_summands
                   coef3_p=>solhrules_product(index1_p(j_2))%coef
                   overlap_xyz_help(1:num) = overlap_xyz(:,i2_p_j2,i0_p_j1,1)
                   overlap_xyz_help(num+1:2*num) = overlap_xyz(:,i2_p_j2,i0_p_j1,2)
                   overlap_xyz_help(2*num+1:3*num) = overlap_xyz(:,i2_p_j2,i0_p_j1,3)
                   overlap_help(1:num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(num+1:2*num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(2*num+1:3*num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(3*num+1:4*num) = overlap(:,i2_p_j2,i0_p_j1)
                   overlap_help(4*num+1:5*num) = overlap(:,i2_p_j2,i0_p_j1)
                   overlap_help(5*num+1:6*num) = overlap(:,i2_p_j2,i0_p_j1)
                   do j_3=1,n_sum3
                      i3_p_j3=index3_p(j_3)
                      i3_p2_j3=index3_p2(j_3)
                      l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                      l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                      l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                      if(l_j3<=l_j1) then
                         coef=coef2_p(j_2)*coef3_p(j_3)*coef1_p(j_1)
                         help_vec1=coef*help_vec2*twoa**l3
                         l_sum=l1+l3
                         meta_min = l_sum * 6 * num + 1
                         meta_max = l_sum * 6 * num + 3 * num
                         help_vec0=help_vec1*&
                              yl_arr(:,i3_p_j3)*diff_arr(:,i3_p2_j3,i_p_j1)
                         prod_arr(:,i,lb_index,l_sum)=prod_arr(:,i,lb_index,l_sum)+&
                              help_vec0*overlap(:,i2_p_j2,i0_p_j1)
                         help_spread(1:num) = help_vec0
                         help_spread(num+1:2*num) = help_vec0
                         help_spread(2*num+1:3*num) = help_vec0
                         prod_arr_gr_vec(meta_min:meta_max,i,lb_index) = &
                              prod_arr_gr_vec(meta_min:meta_max,i,lb_index) + &
                              help_spread(1:3*num) * overlap_xyz_help
                         meta_max = (l_sum+1) * 6 * num
                         if(l_j3<l_j1) then
                            help_vec3=help_vec1*yl_arr(:,i3_p_j3)
                            help_spread(1:num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,1)
                            help_spread(num+1:2*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,2)
                            help_spread(2*num+1:3*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,3)
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_spread)
#endif
                            help_spread(3*num+1:6*num) = help_spread(1:3*num)
                            prod_arr_gr_vec(meta_min:meta_max,i,lb_index) = &
                                 prod_arr_gr_vec(meta_min:meta_max,i,lb_index) + &
                                 help_spread * overlap_help
                         endif
                         help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)
                         help_spread(1:num) = help_vec3*yl_arr_xyz(:,i3_p_j3,1)
                         help_spread(num+1:2*num) = help_vec3*yl_arr_xyz(:,i3_p_j3,2)
                         help_spread(2*num+1:3*num) = help_vec3*yl_arr_xyz(:,i3_p_j3,3)
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_spread)
#endif
                         help_spread(3*num+1:6*num) = help_spread(1:3*num)
                         prod_arr_gr_vec(meta_min:meta_max,i,lb_index) = &
                              prod_arr_gr_vec(meta_min:meta_max,i,lb_index) + &
                              help_spread * overlap_help
                         help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)*&
                              overlap(:,i2_p_j2,i0_p_j1)
                         help_spread(1:num) = help_vec3
                         help_spread(num+1:2*num) = help_vec3
                         help_spread(2*num+1:3*num) = help_vec3
                         help_arr2(:,l_sum,i) = help_arr2(:,l_sum,i) + &
                              help_spread(1:3*num) * help_arr(:,i3_p_j3)
                      endif
                   end do! j_3
                end do! j_2
             enddo! j_1
          end do! i=1
          help_spread(1:num) = twoa
          help_spread(num+1:2*num) = twoa
          help_spread(2*num+1:3*num) = twoa
          help_spread(3*num+1:4*num) = help_vec
          help_spread(4*num+1:5*num) = help_vec
          help_spread(5*num+1:6*num) = help_vec
          do i_ma=1,laposq
             do i_i=1,1+la+i_l
                prod_arr_gr_vec(i_i*6*num+1:(i_i*6+3)*num,i_ma,lb_index) = &
                     prod_arr_gr_vec(i_i*6*num+1:(i_i*6+3)*num,i_ma,lb_index) + &
                     help_spread(1:3*num) * help_arr2(:,i_i-1,i_ma)
                prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) = &
                     prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) + &
                     help_spread(3*num+1:6*num) * help_arr2(:,i_i-1,i_ma)
             enddo
          end do
       end do! mb
    end do! i_l
    call map_prod_arr_gr(la,lb)
  end subroutine calc_prod_arr
 !**************************************************************

!  subroutine calc_prod_arr_opt(la,lb,l1_max,l3_max,fact1,fact2)
  subroutine calc_prod_arr_opt(la,lb,l1_max,l3_max,fact12)
  integer(kind=i4_kind),intent(in):: la,lb
    ! Purpose: piece of the code
    !          calculate quantity prod_arr
    real(kind=r8_kind) :: help_arr(3*num,(ly_max+1)**2),help_arr2(3*num,0:la+lb,laposq)
    integer(kind=i4_kind) :: l_sum,i3_p_j3,i3_p2_j3,i_p_j1,i2_p_j2,i0_p_j1,meta_min,meta_max
    real(kind=r8_kind),dimension(num) :: help_vec1,help_vec0,help_vec3
    integer(kind=i4_kind) :: i
    real(kind=r8_kind) :: &
         overlap_xyz_help(num*3), &
         overlap_help(num*6), &
         help_spread(num*6)
    integer(kind=i4_kind),intent(in) :: l1_max,l3_max

!    real(kind=r8_kind),dimension(num,0:l3_max), intent(in) :: fact1
!    real(kind=r8_kind),dimension(num,0:l1_max), intent(in) :: fact2

    real(kind=r8_kind),dimension(num,0:l3_max,0:l1_max), intent(in) :: fact12

    type(nested2_vars),pointer, dimension(:) :: nested2

    prod_arr=0.0_r8_kind
    prod_arr_gr_vec=0.0_r8_kind

    do l=1,(ly_max+1)**2
       do k_gr=1,3
          help_arr((k_gr-1)*num+1:k_gr*num,l)=yl_arr(:,l)*(gamma_arg(:,k_gr)-xc(k_gr))
       end do
    end do

    help_vec=-2.0_r8_kind*fact0

    do i_l=0,lb
    do mb=1,2*i_l+1
    lb_index=i_l**2+mb
    help_arr2=0.0_r8_kind

    nested2=>nested2_summands(lb_index)%summands
    do i=1,nested2_summands(lb_index)%n
     l_sum=nested2(i)%l1+nested2(i)%l3

!     help_vec1=nested2(i)%coef_prod*fact1(:,nested2(i)%l3)*fact2(:,nested2(i)%l1)
     help_vec1=nested2(i)%coef_prod*fact12(:,nested2(i)%l3,nested2(i)%l1)

     i3_p_j3=nested2(i)%index3_p
     i3_p2_j3=nested2(i)%index3_p2
     i_p_j1=nested2(i)%index_p

     help_vec0=help_vec1*yl_arr(:,i3_p_j3)*diff_arr(:,i3_p2_j3,i_p_j1)

     i2_p_j2=nested2(i)%index2_p
     i0_p_j1=nested2(i)%index0_p

     prod_arr(:,nested2(i)%i,lb_index,l_sum)= &
          prod_arr(:,nested2(i)%i,lb_index,l_sum)+&
                              help_vec0* overlap(:,i2_p_j2,i0_p_j1)

     overlap_xyz_help(1:num) =       overlap_xyz(:,i2_p_j2,i0_p_j1,1) *help_vec0
     overlap_xyz_help(num+1:2*num)=  overlap_xyz(:,i2_p_j2,i0_p_j1,2) *help_vec0
     overlap_xyz_help(2*num+1:3*num)=overlap_xyz(:,i2_p_j2,i0_p_j1,3) *help_vec0

     meta_min = l_sum * 6 * num + 1
     meta_max = l_sum * 6 * num + 3 * num

!     help_spread(1:num) = help_vec0
!     help_spread(num+1:2*num) = help_vec0
!     help_spread(2*num+1:3*num) = help_vec0

  ! first contrib
    prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) = &
      prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index)+overlap_xyz_help
!                              overlap_xyz_help * help_spread(1:3*num)
!  help_spread is not required if only this contribution is in question

    meta_max = (l_sum+1) * 6 * num

    overlap_help(1:num) =         help_mat(:,i2_p_j2,i0_p_j1)
    overlap_help(num+1:2*num) =   help_mat(:,i2_p_j2,i0_p_j1)
    overlap_help(2*num+1:3*num) = help_mat(:,i2_p_j2,i0_p_j1)

    overlap_help(3*num+1:4*num) = overlap(:,i2_p_j2,i0_p_j1)
    overlap_help(4*num+1:5*num) = overlap(:,i2_p_j2,i0_p_j1)
    overlap_help(5*num+1:6*num) = overlap(:,i2_p_j2,i0_p_j1)

!   ! second contrib
!    help_vec3=help_vec1*yl_arr(:,i3_p_j3)
!    help_spread(1:num) =         help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,1)
!    help_spread(num+1:2*num) =   help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,2)
!    help_spread(2*num+1:3*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,3)
!#if defined(_VECTOR) && defined(_VPP)
!!OCL NOVREC(help_spread)
!#endif
!    help_spread(3*num+1:6*num) = help_spread(1:3*num)

!     prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) = &
!          prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) + &
!                                 help_spread * overlap_help
!

     help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)

     help_spread(1:num) =         help_vec3*yl_arr_xyz(:,i3_p_j3,1)
     help_spread(num+1:2*num) =   help_vec3*yl_arr_xyz(:,i3_p_j3,2)
     help_spread(2*num+1:3*num) = help_vec3*yl_arr_xyz(:,i3_p_j3,3)

#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_spread)
#endif
     help_spread(3*num+1:6*num) = help_spread(1:3*num)

     !third contrib
      prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) = &
      prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) + &
                              help_spread * overlap_help

      help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)*&
                              overlap(:,i2_p_j2,i0_p_j1)
      help_spread(1:num) = help_vec3
      help_spread(num+1:2*num) = help_vec3
      help_spread(2*num+1:3*num) = help_vec3

      help_arr2(:,l_sum,nested2(i)%i) = help_arr2(:,l_sum,nested2(i)%i) + &
                                 help_spread(1:3*num) * help_arr(:,i3_p_j3)

    end do

    nested2=>nested2ne_summands(lb_index)%summands
    do i=1,nested2ne_summands(lb_index)%n
     l_sum=nested2(i)%l1+nested2(i)%l3

!     help_vec1=nested2(i)%coef_prod*fact1(:,nested2(i)%l3)*fact2(:,nested2(i)%l1)
     help_vec1=nested2(i)%coef_prod*fact12(:,nested2(i)%l3,nested2(i)%l1)

     i3_p_j3=nested2(i)%index3_p
     i3_p2_j3=nested2(i)%index3_p2
     i_p_j1=nested2(i)%index_p
     help_vec0=help_vec1*yl_arr(:,i3_p_j3)*diff_arr(:,i3_p2_j3,i_p_j1)
     i2_p_j2=nested2(i)%index2_p
     i0_p_j1=nested2(i)%index0_p

     meta_min = l_sum * 6 * num + 1
     meta_max = (l_sum+1) * 6 * num
     help_vec3=help_vec1*yl_arr(:,i3_p_j3)

     help_spread(1:num) =         help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,1)
     help_spread(num+1:2*num) =   help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,2)
     help_spread(2*num+1:3*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,3)
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_spread)
#endif
     help_spread(3*num+1:6*num) = help_spread(1:3*num)

     overlap_help(1:num) =         help_mat(:,i2_p_j2,i0_p_j1)
     overlap_help(num+1:2*num) =   help_mat(:,i2_p_j2,i0_p_j1)
     overlap_help(2*num+1:3*num) = help_mat(:,i2_p_j2,i0_p_j1)

     overlap_help(3*num+1:4*num) = overlap(:,i2_p_j2,i0_p_j1)
     overlap_help(4*num+1:5*num) = overlap(:,i2_p_j2,i0_p_j1)
     overlap_help(5*num+1:6*num) = overlap(:,i2_p_j2,i0_p_j1)
   ! second contrib
     prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) = &
          prod_arr_gr_vec(meta_min:meta_max,nested2(i)%i,lb_index) + &
                                 help_spread * overlap_help

    end do


     help_spread(1:num) = twoa
     help_spread(num+1:2*num) = twoa
     help_spread(2*num+1:3*num) = twoa

     help_spread(3*num+1:4*num) = help_vec
     help_spread(4*num+1:5*num) = help_vec
     help_spread(5*num+1:6*num) = help_vec

      do i_ma=1,laposq
        do i_i=1,1+la+i_l
           !forth contrib
           prod_arr_gr_vec(i_i*6*num+1:(i_i*6+3)*num,i_ma,lb_index) = &
               prod_arr_gr_vec(i_i*6*num+1:(i_i*6+3)*num,i_ma,lb_index) + &
                     help_spread(1:3*num) * help_arr2(:,i_i-1,i_ma)
           prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) = &
               prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) + &
                     help_spread(3*num+1:6*num) * help_arr2(:,i_i-1,i_ma)
             enddo
          end do
        enddo
       enddo

    call map_prod_arr_gr(la,lb)

  end subroutine calc_prod_arr_opt


 !**************************************************************
 subroutine map_prod_arr_gr(la,lb)
  integer(kind=i4_kind),intent(in):: la,lb
   ! Purpose: piece of the code
   ! maps  prod_arr_gr_vec to  prod_arr_gr
   integer(kind=i4_kind) :: i_ma, i_mb, k_gr, l, i_min, i_max
   do i_mb=1,lbposq
      do i_ma=1,laposq
         i_min = 1
         i_max = num
         do l = 0, la+lb+1
            do k_gr=1,6
               prod_arr_gr(:,k_gr,i_ma,i_mb,l) = &
                    prod_arr_gr_vec(i_min:i_max,i_ma,i_mb)
               i_min = i_min + num
               i_max = i_max + num
            end do
         end do
      end do
   end do
 end subroutine map_prod_arr_gr
 !**************************************************************

  !**************************************************************
  subroutine calc_intermediate(la,lb)
  integer(kind=i4_kind),intent(in):: la,lb
    ! Purpose: piece of the code
    !          calculates intermediate, which is necesarry to calculate the
    !          3center integrals with l>0
    integer(kind=i4_kind) :: i,j_1,j_2,lm_min,lm_max,lm_in,n_sum1,n_sum2,l1,l2,&
         i2_2p_j2,i1_2p_j1,i2_1p_j2,i1_1p_j1,l_index,i_sp_min,i_sp_max, step
    integer(kind=i4_kind),pointer,dimension(:) :: index1_1p,&
         index1_2p,index2_1p,index2_2p
    real(kind=r8_kind),pointer :: coeff_p1(:),coeff_p2(:)
    real(kind=r8_kind) :: &
         coeff_loc, &
         help_arr0(num*3,0:(la+lb+1),2*la+1), &
         help_diff0(num*3), &
         help_prod0(num*3), &
         help_spread1(num*6), &
         help_spread2(num*6*(la+lb+2))

    mb_loop: do mb=1,2*lb+1
       help_arr0=0.0_r8_kind
       lm_in=(lb**2)+mb
       lm_min=la**2+1
       lm_max=(la+1)**2
       index1_1p=>solhrules_product(lm_in)%lm_sh1
       index1_2p=>solhrules_product(lm_in)%lm_sh2
       coeff_p1=>solhrules_product(lm_in)%coef
       n_sum1=solhrules_product(lm_in)%n_summands
       step = num * 6
       do i=lm_min,lm_max
          l_index=i+1-lm_min
          n_sum2=solhrules_product(i)%n_summands
          index2_1p=>solhrules_product(i)%lm_sh1
          index2_2p=>solhrules_product(i)%lm_sh2
          coeff_p2=>solhrules_product(i)%coef
          do j_1=1,n_sum1
             i1_2p_j1=index1_2p(j_1)
             i1_1p_j1=index1_1p(j_1)
             do j_2=1,n_sum2
                coeff_loc=coeff_p2(j_2)*coeff_p1(j_1)*coeff(i_cnt)
                i2_2p_j2=index2_2p(j_2)
                i2_1p_j2=index2_1p(j_2)
                l1=solhrules_l_and_m_of_lm(1,index2_1p(j_2))
                l2=solhrules_l_and_m_of_lm(1,index1_1p(j_1))
                if(l1+l2<=i_l) then
                ! sum up help_arr0
                if(l1+l2<i_l) then
                help_diff0(      1:  num) = diff_arr_xyz(:,i2_1p_j2,i1_1p_j1,1)
                help_diff0(  num+1:2*num) = diff_arr_xyz(:,i2_1p_j2,i1_1p_j1,2)
                help_diff0(2*num+1:3*num) = diff_arr_xyz(:,i2_1p_j2,i1_1p_j1,3)
                help_diff0 = coeff_loc * help_diff0
                do l=0,la+lb
                   help_prod0(      1:  num) = prod_arr(:,i2_2p_j2,i1_2p_j1,l)
                   help_prod0(  num+1:2*num) = prod_arr(:,i2_2p_j2,i1_2p_j1,l)
                   help_prod0(2*num+1:3*num) = prod_arr(:,i2_2p_j2,i1_2p_j1,l)
                   help_arr0(:,l,l_index) = help_arr0(:,l,l_index) + &
                        help_prod0 * help_diff0
                enddo
                endif
                ! sum up intermediate_gr
                help_vec=coeff_loc*diff_arr(:,i2_1p_j2,i1_1p_j1)
                i_sp_min = 1
                i_sp_max = num
                do k_gr=1,6
                   help_spread1(i_sp_min:i_sp_max) = help_vec
                   i_sp_min = i_sp_min + num
                   i_sp_max = i_sp_max + num
                enddo
                i_sp_min = 1
                i_sp_max = step
                do l=0,la+lb+1
                   help_spread2(i_sp_min:i_sp_max) = help_spread1
                   i_sp_min = i_sp_min + step
                   i_sp_max = i_sp_max + step
                enddo
                intermediate_gr(:,mb,l_index,i_ind) = intermediate_gr(:,mb,l_index,i_ind) + &
                     help_spread2 * prod_arr_gr_vec(:,i2_2p_j2,i1_2p_j1)
               endif
             end do! j_2
          enddo! j_1
       end do! i
       help_prod0(1:num) = aqapb(:,1)
       help_prod0(num+1:2*num) = aqapb(:,1)
       help_prod0(2*num+1:3*num) = aqapb(:,1)
       step = num * 3
       do i_ma=1,2*la+1
          i_sp_min = 1
          i_sp_max = step
          do l=0,la+lb
             intermediate_gr(i_sp_min:i_sp_max,mb,i_ma,i_ind) = &
                  intermediate_gr(i_sp_min:i_sp_max,mb,i_ma,i_ind) + &
                  help_arr0(:,l,i_ma) * help_prod0
             i_sp_min = i_sp_min + step
             i_sp_max = i_sp_max + step
             intermediate_gr(i_sp_min:i_sp_max,mb,i_ma,i_ind) = &
                  intermediate_gr(i_sp_min:i_sp_max,mb,i_ma,i_ind) + &
                  help_arr0(:,l,i_ma)
             i_sp_min = i_sp_min + step
             i_sp_max = i_sp_max + step
          end do
       end do
    end do mb_loop
  end subroutine calc_intermediate
  !****************************************************************

  !**************************************************************
  subroutine calc_final_nuc_and_pvsp_grad(gradient_ab,gradient_c)
    ! Purpose: piece of the code
    !          calculates the final gradient of nuc and pvsp matrixelements
    real(kind=r8_kind) :: gradient_ab(:,:,:,:)
    ! derivatives with respect to centers a and b
    real(kind=r8_kind), optional :: gradient_c(:,:,:,:)
    ! derivative with respect to the center c
    do i_la=1,nm_la
       do i_lb=1,nm_lb
          do k_gr=1,3
             ! derivative with respect to first center a
             if (moving_a) then
             gradient_ab(:,i_lb,i_la,k_gr)=gradient_ab(:,i_lb,i_la,k_gr)+&
                  help_arr_gr(:,k_gr,i_lb,i_la)
             end if
             ! derivative with respect to swecond center b
             ! use that sum of derivatives is zero
             if (moving_b) then
             gradient_ab(:,i_lb,i_la,3+k_gr)=gradient_ab(:,i_lb,i_la,3+k_gr)&
                  +help_arr_gr(:,3+k_gr,i_lb,i_la)-help_arr_gr(:,k_gr,i_lb,i_la)
             end if
          enddo
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
                gradient_c(:,i_lb,i_la,i_grad)=gradient_c(:,i_lb,i_la,i_grad)-&
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
             gradient_c(:,i_lb,i_la,1)=gradient_c(:,i_lb,i_la,1)-&
                  help_arr_gr(:,4,i_lb,i_la)
             gradient_c(:,i_lb,i_la,2)=gradient_c(:,i_lb,i_la,2)-&
                  help_arr_gr(:,5,i_lb,i_la)
             gradient_c(:,i_lb,i_la,3)=gradient_c(:,i_lb,i_la,3)-&
                  help_arr_gr(:,6,i_lb,i_la)
             counter=counter+1
          end do
       enddo
    end if
  end subroutine calc_final_nuc_and_pvsp_grad

  subroutine add_cpks_coul(equalb,grad_mat,compact_storage)
   use gradient_data_module,only:cpks_gradient_totalsym
#ifdef _COMPAC_FORTRAN
    use unique_atom_module
#endif
    ! Purpose: Add gradient of nuclear attraction or pvsp
    !          to the correct datatype pointer_prim_int
    integer(kind=i4_kind),intent(in):: equalb
#ifdef no_cpks_coul_grads
    real(kind=r8_kind) :: grad_mat(:,:)
#else
    real(kind=r8_kind) :: grad_mat(:,:,:,:,:,:)
#endif
    logical, optional :: compact_storage
    integer(i4_kind) :: index

    ! LOAD dsym/dRa < xi_i | ... | xi_j >
    ! Input: grad_mat(...,1:3) , Output: pointer_prim_int(ia:) or
    !                                    pointer_prim_int(1:N)
    if (moving_a) then
       if (present(compact_storage)) then
          index = 1
       else
          index = gradient_index(ima)
       endif

       grad_dim = gradient_index(ima+1) - gradient_index(ima)

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(ima) % m(:, :, 1)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_a
#ifdef no_cpks_coul_grads
       cpks_gradient_totalsym(:,index)=cpks_gradient_totalsym(:,index)+ &
               rotmat(i_grad,1)*grad_mat(:,1)+&
               rotmat(i_grad,2)*grad_mat(:,2)+&
               rotmat(i_grad,3)*grad_mat(:,3)
#else
          pointer_prim_3c(index)%m=pointer_prim_3c(index)%m+&   !(1)
               rotmat(i_grad,1)*grad_mat(:,:,:,:,:,1)+&
               rotmat(i_grad,2)*grad_mat(:,:,:,:,:,2)+&
               rotmat(i_grad,3)*grad_mat(:,:,:,:,:,3)
#endif

          index=index+1
       end do
    else
#ifdef no_cpks_coul_grads
       cpks_gradient_totalsym(:,index)=cpks_gradient_totalsym(:,index)+grad_mat(:,1)
       cpks_gradient_totalsym(:,index+1)=cpks_gradient_totalsym(:,index+1)+grad_mat(:,2)
       cpks_gradient_totalsym(:,index+2)=cpks_gradient_totalsym(:,index+2)+grad_mat(:,3)
#else
       pointer_prim_3c(index)%m=pointer_prim_3c(index)%m+&
            grad_mat(:,:,:,:,:,1)
       pointer_prim_3c(index+1)%m=pointer_prim_3c(index+1)%m+&
            grad_mat(:,:,:,:,:,2)
       pointer_prim_3c(index+2)%m=pointer_prim_3c(index+2)%m+&
            grad_mat(:,:,:,:,:,3)
#endif
    endif
    else
      grad_dim = 0 ! required for compact storage mode
    endif

    ! LOAD dsym/dRb < xi_i | ... | xi_j >
    ! Input: grad_mat(...,4:6) , Output: pointer_prim_int(ib:) or
    !                                    pointer_prim_int(N+1:M)
    if (moving_b) then
       if (present(compact_storage)) then
          index = grad_dim + 1
       else
          index = gradient_index(imb)
       endif
       grad_dim = gradient_index(imb+1) - gradient_index(imb)

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(imb) % m(:, :, equalb)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_b
#ifdef no_cpks_coul_grads
       cpks_gradient_totalsym(:,index)=cpks_gradient_totalsym(:,index)+ &
               rotmat(i_grad,1)*grad_mat(:,4)+&
               rotmat(i_grad,2)*grad_mat(:,5)+&
               rotmat(i_grad,3)*grad_mat(:,6)
#else
          pointer_prim_3c(index)%m=pointer_prim_3c(index)%m+&
               rotmat(i_grad,1)*grad_mat(:,:,:,:,:,4)+&
               rotmat(i_grad,2)*grad_mat(:,:,:,:,:,5)+&
               rotmat(i_grad,3)*grad_mat(:,:,:,:,:,6)
#endif
          index=index+1
       end do
    else
#ifdef no_cpks_coul_grads
     cpks_gradient_totalsym(:,index)=cpks_gradient_totalsym(:,index)+grad_mat(:,4)
     cpks_gradient_totalsym(:,index+1)=cpks_gradient_totalsym(:,index+1)+grad_mat(:,5)
     cpks_gradient_totalsym(:,index+2)=cpks_gradient_totalsym(:,index+2)+grad_mat(:,6)
#else
       pointer_prim_3c(index)%m=pointer_prim_3c(index)%m+&
            grad_mat(:,:,:,:,:,4)
       pointer_prim_3c(index+1)%m=pointer_prim_3c(index+1)%m+&
            grad_mat(:,:,:,:,:,5)
       pointer_prim_3c(index+2)%m=pointer_prim_3c(index+2)%m+&
            grad_mat(:,:,:,:,:,6)
#endif
    endif
    endif

  end subroutine add_cpks_coul


  subroutine add_nuc_or_pvsp(equalb,grad_mat,compact_storage)
#ifdef _COMPAC_FORTRAN
    use unique_atom_module
#endif
    ! Purpose: Add gradient of nuclear attraction or pvsp
    !          to the correct datatype pointer_prim_int
    integer(kind=i4_kind),intent(in):: equalb
    real(kind=r8_kind) :: grad_mat(:,:,:,:,:)
    logical, optional :: compact_storage

    ! LOAD dsym/dRa < xi_i | ... | xi_j >
    ! Input: grad_mat(...,1:3) , Output: pointer_prim_int(ia:) or
    !                                    pointer_prim_int(1:N)
    if (moving_a) then
       if (present(compact_storage)) then
          index = 1
       else
          index = gradient_index(ima)
       endif

       grad_dim = gradient_index(ima+1) - gradient_index(ima)

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(ima) % m(:, :, 1)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_a
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               rotmat(i_grad,1)*grad_mat(:,:,:,:,1)+&
               rotmat(i_grad,2)*grad_mat(:,:,:,:,2)+&
               rotmat(i_grad,3)*grad_mat(:,:,:,:,3)
          index=index+1
       end do
    else
       pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
            grad_mat(:,:,:,:,1)
       pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
            grad_mat(:,:,:,:,2)
       pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
            grad_mat(:,:,:,:,3)
    endif
    else
      grad_dim = 0 ! required for compact storage mode
    endif

    ! LOAD dsym/dRb < xi_i | ... | xi_j >
    ! Input: grad_mat(...,4:6) , Output: pointer_prim_int(ib:) or
    !                                    pointer_prim_int(N+1:M)
    if (moving_b) then
       if (present(compact_storage)) then
          index = grad_dim + 1
       else
          index = gradient_index(imb)
       endif
       grad_dim = gradient_index(imb+1) - gradient_index(imb)

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(imb) % m(:, :, equalb)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
       do i_grad=1,grad_dim ! only if moving_b
          pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
               rotmat(i_grad,1)*grad_mat(:,:,:,:,4)+&
               rotmat(i_grad,2)*grad_mat(:,:,:,:,5)+&
               rotmat(i_grad,3)*grad_mat(:,:,:,:,6)
          index=index+1
       end do
    else
       pointer_prim_int(index)%m=pointer_prim_int(index)%m+&
            grad_mat(:,:,:,:,4)
       pointer_prim_int(index+1)%m=pointer_prim_int(index+1)%m+&
            grad_mat(:,:,:,:,5)
       pointer_prim_int(index+2)%m=pointer_prim_int(index+2)%m+&
            grad_mat(:,:,:,:,6)
    endif
    endif

  end subroutine add_nuc_or_pvsp

  subroutine add_dervs(equalb,dervs_mat,ima,imb,pointer_prim_dervs,compact_storage)
  use cpksdervs_matrices
#ifdef _COMPAC_FORTRAN
    use unique_atom_module
#endif

    ! Purpose: Add 2nd dervs
    !          to the correct datatype pointer_prim_dervs
    implicit none
    integer(kind=i4_kind),intent(in):: equalb,ima,imb
    real(kind=r8_kind), intent(in)     :: dervs_mat(:,:,:,:,:,:)
    type(arrmat4)     , intent(inout)  :: pointer_prim_dervs(:,:) ! primitive integrals
    logical, optional :: compact_storage
    ! *** end of interface **

    ! local vars
    integer(kind=i4_kind):: sz(4),grad_dim_a,grad_dim_b,ind_a,ind_b,k2dr
    real(kind=r8_kind),allocatable:: dervs_mat_temp_a(:,:,:,:,:,:), &
                                     dervs_mat_temp_b(:,:,:,:,:,:)


    ! LOAD dsym/dRa < xi_i | ... | xi_j >
    ! Input: grad_mat(...,1:3) , Output: pointer_prim_int(ia:) or
    !                                    pointer_prim_int(1:N)

    if(moving_a) then
     grad_dim_a = gradient_index(ima+1) - gradient_index(ima)
          ind_a = gradient_index(ima)-1
    else
     grad_dim_a=0
    endif

    if(moving_b) then
     grad_dim_b = gradient_index(imb+1) - gradient_index(imb)
          ind_b= gradient_index(imb)-1
    else
     grad_dim_b=0
    endif

!    print*,' dervs_mat in add_dervs',ima,imb,equalb, &
!    print*, sum(dervs_mat(:,:,:,:,i_grad,1)), &
!            sum(dervs_mat(:,:,:,:,i_grad,2)), &
!            sum(dervs_mat(:,:,:,:,i_grad,3)), &
!            sum(dervs_mat(:,:,:,:,i_grad,4)), &
!            sum(dervs_mat(:,:,:,:,i_grad,5)), &
!            sum(dervs_mat(:,:,:,:,i_grad,6))
!   enddo
!   dervs_mat(:,:,:,:,4:6,4:6)=dervs_mat(:,:,:,:,1:3,1:3)
!   dervs_mat(:,:,:,:,1:3,4:6)=-dervs_mat(:,:,:,:,1:3,1:3)
!   dervs_mat(:,:,:,:,4:6,1:3)=-dervs_mat(:,:,:,:,1:3,1:3)


       sz(:)=shape(dervs_mat(:,:,:,:,1,1))
       allocate( dervs_mat_temp_a(sz(1),sz(2),sz(3),sz(4),6,grad_dim_a),&
                 dervs_mat_temp_b(sz(1),sz(2),sz(3),sz(4),6,grad_dim_b),&
                 stat=cpksalloc(121))
       ASSERT(cpksalloc(121).eq.0)
       dervs_mat_temp_a=0.0_r8_kind
       dervs_mat_temp_b=0.0_r8_kind

    if (moving_a) then

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(ima) % m(:, :, 1)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

     if(do_rotation) then
       do i_grad=1,grad_dim_a
!       print*,'dervs rotmat', rotmat(i_grad,:)
        dervs_mat_temp_a(:,:,:,:,1:6,i_grad)=dervs_mat_temp_a(:,:,:,:,1:6,i_grad) &
              +unique_atom_grad_info(ima)%m(i_grad,1,1)*dervs_mat(:,:,:,:,1:6,1) &
              +unique_atom_grad_info(ima)%m(i_grad,2,1)*dervs_mat(:,:,:,:,1:6,2) &
              +unique_atom_grad_info(ima)%m(i_grad,3,1)*dervs_mat(:,:,:,:,1:6,3)
       enddo

     else
        dervs_mat_temp_a(:,:,:,:,1:6,1:3)=dervs_mat_temp_a(:,:,:,:,1:6,1:3)+dervs_mat(:,:,:,:,1:6,1:3)
     endif
    else
      grad_dim_a = 0 ! required for compact storage mode
    endif

    ! This might be a 0 x 3 matrix. In general grad_dim x 3:
    rotmat => unique_atom_grad_info(imb) % m(:, :, equalb)

    ! Optimization only. No need to apply unity matrix:
    do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
       do i_grad=1,grad_dim_b
!        print*,' rotmat b 1',rotmat(i_grad,:)
        dervs_mat_temp_b(:,:,:,:,1:6,i_grad)= &
              dervs_mat_temp_b(:,:,:,:,1:6,i_grad) &
              +unique_atom_grad_info(imb)%m(i_grad,1,equalb)*dervs_mat(:,:,:,:,1:6,4) &
              +unique_atom_grad_info(imb)%m(i_grad,2,equalb)*dervs_mat(:,:,:,:,1:6,5) &
              +unique_atom_grad_info(imb)%m(i_grad,3,equalb)*dervs_mat(:,:,:,:,1:6,6)
       end do
    else
        dervs_mat_temp_b(:,:,:,:,1:6,4:6)= &
              dervs_mat_temp_b(:,:,:,:,1:6,4:6)+dervs_mat(:,:,:,:,1:6,4:6)

    endif




     ! ***  NOW TRANSFOR ON GRAD INDEX  AND REMAP

    if (moving_a) then

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(ima) % m(:, :, 1)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

     if(do_rotation) then
      do i_grad=1,grad_dim_a
!        print*,' rotmat a 2',rotmat(i_grad,:)
      do k2dr=1,grad_dim_a
          pointer_prim_dervs(ind_a+i_grad,ind_a+k2dr)%m= &
                  pointer_prim_dervs(ind_a+i_grad,ind_a+k2dr)%m&
         +unique_atom_grad_info(ima)%m(i_grad,1,1)*dervs_mat_temp_a(:,:,:,:,1,k2dr) &
         +unique_atom_grad_info(ima)%m(i_grad,2,1)*dervs_mat_temp_a(:,:,:,:,2,k2dr) &
         +unique_atom_grad_info(ima)%m(i_grad,3,1)*dervs_mat_temp_a(:,:,:,:,3,k2dr)
       enddo
      enddo

     do i_grad=1,grad_dim_a
!       print*,' rotmat b 2',rotmat(i_grad,:)
      do k2dr=1,grad_dim_b
          pointer_prim_dervs(ind_a+i_grad,ind_b+k2dr)%m= &
                             pointer_prim_dervs(ind_a+i_grad,ind_b+k2dr)%m&
              +rotmat(i_grad,1)*dervs_mat_temp_b(:,:,:,:,1,k2dr) &
              +rotmat(i_grad,2)*dervs_mat_temp_b(:,:,:,:,2,k2dr) &
              +rotmat(i_grad,3)*dervs_mat_temp_b(:,:,:,:,3,k2dr)
       end do
      enddo

     else
       do k2dr=1,grad_dim_a
       do i_grad=1,grad_dim_a
          pointer_prim_dervs(ind_a+i_grad,ind_a+k2dr)%m= &
             pointer_prim_dervs(ind_a+i_grad,ind_a+k2dr)%m&
            +dervs_mat_temp_a(:,:,:,:,i_grad,k2dr)
       enddo
       enddo

       do k2dr=1,grad_dim_b
       do i_grad=1,grad_dim_a
          pointer_prim_dervs(ind_a+i_grad,ind_b+k2dr)%m= &
              pointer_prim_dervs(ind_a+i_grad,ind_b+k2dr)%m&
             +dervs_mat_temp_b(:,:,:,:,i_grad,k2dr)
       enddo
       enddo

     endif
    else
      grad_dim_a = 0 ! required for compact storage mode
    endif


    mov_b: if (moving_b) then

       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(imb) % m(:, :, equalb)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
      do i_grad=1,grad_dim_b
      do k2dr=1,grad_dim_b
          pointer_prim_dervs(ind_b+i_grad,ind_b+k2dr)%m= &
                            pointer_prim_dervs(ind_b+i_grad,ind_b+k2dr)%m+&
          unique_atom_grad_info(imb)%m(i_grad,1,equalb)*dervs_mat_temp_b(:,:,:,:,4,k2dr)+&
          unique_atom_grad_info(imb)%m(i_grad,2,equalb)*dervs_mat_temp_b(:,:,:,:,5,k2dr)+&
          unique_atom_grad_info(imb)%m(i_grad,3,equalb)*dervs_mat_temp_b(:,:,:,:,6,k2dr)
       end do
    enddo

    do k2dr=1,grad_dim_a
!    print*,sum(dervs_mat_temp_a(:,:,:,:,4,k2dr)), &
!           sum(dervs_mat_temp_a(:,:,:,:,5,k2dr)), &
!           sum(dervs_mat_temp_a(:,:,:,:,6,k2dr))

    do i_grad=1,grad_dim_b
          pointer_prim_dervs(ind_b+i_grad,ind_a+k2dr)%m= &
                 pointer_prim_dervs(ind_b+i_grad,ind_a+k2dr)%m+&
               unique_atom_grad_info(imb)%m(i_grad,1,equalb)*dervs_mat_temp_a(:,:,:,:,4,k2dr)+&
               unique_atom_grad_info(imb)%m(i_grad,2,equalb)*dervs_mat_temp_a(:,:,:,:,5,k2dr)+&
               unique_atom_grad_info(imb)%m(i_grad,3,equalb)*dervs_mat_temp_a(:,:,:,:,6,k2dr)
       end do
    enddo

    else
      do k2dr=1,grad_dim_b
       do i_grad=1,grad_dim_b
       pointer_prim_dervs(ind_b+i_grad,ind_b+k2dr)%m= &
                         pointer_prim_dervs(ind_b+i_grad,ind_b+k2dr)%m+ &
               dervs_mat_temp_b(:,:,:,:,i_grad+3,k2dr)
       enddo
      enddo

      do k2dr=1,grad_dim_a
       do i_grad=1,grad_dim_b
       pointer_prim_dervs(ind_b+i_grad,ind_a+k2dr)%m= &
             pointer_prim_dervs(ind_b+i_grad,ind_a+k2dr)%m+ &
               dervs_mat_temp_a(:,:,:,:,i_grad+3,k2dr)
       enddo
      enddo

    endif
    endif mov_b
  deallocate(dervs_mat_temp_a,dervs_mat_temp_b,stat=cpksalloc(121))
  ASSERT(cpksalloc(121).eq.0)
  cpksalloc(121)=1

  end subroutine add_dervs

  subroutine add_ca_dervs(equalb,dervs_mat,ima,imb,imc,pointer_prim_dervs,compact_storage)
#ifdef _COMPAC_FORTRAN
    use unique_atom_module
#endif
    ! Purpose: Add 2nd dervs
    !          to the correct datatype pointer_prim_dervs
    implicit none
    integer(kind=i4_kind),intent(in):: equalb,ima,imb,imc
    real(kind=r8_kind), intent(in)     :: dervs_mat(:,:,:,:,:,:)
    type(arrmat4)     , intent(inout)  :: pointer_prim_dervs(:,:) ! primitive integrals
    logical, optional :: compact_storage
    ! *** end of interface **

    ! local vars
    integer(kind=i4_kind):: grad_dim_a,grad_dim_b,ind_a,ind_b,k2dr,ind_c,grad_dim_c

    if(moving_c) then
     ind_c=gradient_index(imc)-1
     grad_dim_c=size(ca_dervs_mat,5)
    else
     grad_dim_c=0
    endif

    if(moving_a) then
     grad_dim_a = gradient_index(ima+1) - gradient_index(ima)
       if (present(compact_storage)) then
          ind_a = 1
       else
          ind_a = gradient_index(ima)-1
       endif
    else
     grad_dim_a=0
    endif

    if(moving_b) then
     grad_dim_b = gradient_index(imb+1) - gradient_index(imb)

       if (present(compact_storage)) then
          ind_b = grad_dim_b + 1
       else
          ind_b= gradient_index(imb)-1
       endif
    else
     grad_dim_b=0
    endif


     ! ***  NOW TRANSFOR ON GRAD INDEXE a & b  AND REMAP

    if (moving_a) then
       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(ima) % m(:, :, 1)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

     if(do_rotation) then
      do i_grad=1,grad_dim_c
      do k2dr=1,grad_dim_a                                    ! transform a projections
          pointer_prim_dervs(ind_c+i_grad,ind_a+k2dr)%m= &
             pointer_prim_dervs(ind_c+i_grad,ind_a+k2dr)%m &
              +rotmat(k2dr,1)*dervs_mat(:,:,:,:,i_grad,1) &
              +rotmat(k2dr,2)*dervs_mat(:,:,:,:,i_grad,2) &
              +rotmat(k2dr,3)*dervs_mat(:,:,:,:,i_grad,3)

            ! this is a standard order,c goes after a
            ! should contain not transformed (-GGa-GaGb)
            pointer_prim_dervs(ind_a+k2dr,ind_c+i_grad)%m= &
            pointer_prim_dervs(ind_a+k2dr,ind_c+i_grad)%m &
              +rotmat(k2dr,1)*dervs_mat(:,:,:,:,i_grad,1) &
              +rotmat(k2dr,2)*dervs_mat(:,:,:,:,i_grad,2) &
              +rotmat(k2dr,3)*dervs_mat(:,:,:,:,i_grad,3)
       end do
      enddo

     else
       do i_grad=1,grad_dim_c
       do k2dr=1,3
          pointer_prim_dervs(ind_c+i_grad,ind_a+k2dr)%m= &
            pointer_prim_dervs(ind_c+i_grad,ind_a+k2dr)%m&
                           +dervs_mat(:,:,:,:,i_grad,k2dr)

          pointer_prim_dervs(ind_a+k2dr,ind_c+i_grad)%m= &
           pointer_prim_dervs(ind_a+k2dr,ind_c+i_grad)%m &
                           +dervs_mat(:,:,:,:,i_grad,k2dr)
       enddo
       enddo
     endif
    else
      grad_dim_a = 0 ! required for compact storage mode
    endif

    mov_b: if (moving_b) then
       ! This might be a 0 x 3 matrix. In general grad_dim x 3:
       rotmat => unique_atom_grad_info(imb) % m(:, :, equalb)

       ! Optimization only. No need to apply unity matrix:
       do_rotation = .not. unity_matrix_p (rotmat)

    if(do_rotation) then
      do k2dr=1,grad_dim_b      ! transform b projections
       do i_grad=1,grad_dim_c

           ! this corresponds to (7:9,4:6) block and thus we map
           ! untransposed -GaGb as stored in  dervs_mat(i_grad,4:6)

           pointer_prim_dervs(ind_c+i_grad,ind_b+k2dr)%m= &
           pointer_prim_dervs(ind_c+i_grad,ind_b+k2dr)%m  &
              +rotmat(k2dr,1)*dervs_mat(:,:,:,:,i_grad,4) &
              +rotmat(k2dr,2)*dervs_mat(:,:,:,:,i_grad,5) &
              +rotmat(k2dr,3)*dervs_mat(:,:,:,:,i_grad,6)

         ! this goes in stadard order b then c inexes
         ! (4:6,7:9) should containt transposed -GaGb
         ! and this trasposition is done here

          pointer_prim_dervs(ind_b+k2dr,ind_c+i_grad)%m= &
            pointer_prim_dervs(ind_b+k2dr,ind_c+i_grad)%m&
             +rotmat(k2dr,1)*dervs_mat(:,:,:,:,i_grad,4) &
             +rotmat(k2dr,2)*dervs_mat(:,:,:,:,i_grad,5) &
             +rotmat(k2dr,3)*dervs_mat(:,:,:,:,i_grad,6)

       end do
    enddo

    else
      do k2dr=1,3
       do i_grad=1,grad_dim_c
       pointer_prim_dervs(ind_c+i_grad,ind_b+k2dr)%m= &
         pointer_prim_dervs(ind_c+i_grad,ind_b+k2dr)%m+dervs_mat(:,:,:,:,i_grad,k2dr+3)

       pointer_prim_dervs(ind_b+k2dr,ind_c+i_grad)%m= &
         pointer_prim_dervs(ind_b+k2dr,ind_c+i_grad)%m+dervs_mat(:,:,:,:,i_grad,k2dr+3)

       enddo
      enddo

    endif
    endif mov_b

  end subroutine add_ca_dervs

  subroutine s_coulomb(la,lb)
  integer(kind=i4_kind),intent(in):: la,lb
    ! calculate s-type coulomb fit integrals
    integer(kind=i4_kind) :: i_min, i_max, step
    real(kind=r8_kind) :: help_loc(num*6)
    ncexps = unique_atoms(i)%l_ch(0)%n_exponents
    cexps => unique_atoms(i)%l_ch(0)%exponents(:)
    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       ! precalculation of two factors
       rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
       fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))
       gamma_help(:,1:2+la+lb)=gamma(2+la+lb,gamma_arg2(:)*rcsabc(:))
       help_arr_gr_vec=0.0_r8_kind
       step = num * 6
       i_min = 1
       i_max = step
       do l=0,la+lb+1
          help_vec=gamma_help(:,l+1)*fact8*rcsabc(:)**l
          help_loc(1:num) = help_vec
          help_loc(num+1:2*num) = help_vec
          help_loc(2*num+1:3*num) = help_vec
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_loc)
#endif
          help_loc(3*num+1:6*num) = help_loc(1:3*num)
          do i_la=1,nm_la
             la_index=lasq+i_la
             do i_lb=1,nm_lb
                lb_index=lbsq+i_lb
                help_arr_gr_vec(:,i_lb,i_la) = help_arr_gr_vec(:,i_lb,i_la) + &
                     prod_arr_gr_vec(i_min:i_max,la_index,lb_index) * help_loc
             end do! loop over i_lb
          end do! loop over i_lb
          i_min = i_min + step
          i_max = i_max + step
       enddo! loop over l
       call calc_coul_gradients(0,k,1)
    end do! loop over fitfunctions
   end subroutine s_coulomb
  !**************************************************************

  !**************************************************************
   subroutine r2_coulomb(la,lb)
     ! calculate r2-type coulomb fit integrals
   integer(kind=i4_kind),intent(in):: la,lb
     integer(kind=i4_kind) :: i_min, i_max, step
     real(kind=r8_kind) :: help_loc(num*6)
     ncexps = unique_atoms(i)%r2_ch%n_exponents
     cexps => unique_atoms(i)%r2_ch%exponents(:)
     call integral_interrupt_2cob3c()
     do k=1,ncexps ! loop over fitexponents
        ! precalculation of two factors
        rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
        fact8=2.0_r8_kind*pi/(cexps(k)**2)*sqrt(fact0/(fact0+cexps(k)))&
             *(fact0/(fact0+cexps(k)))
        fact7=(fact0+1.5_r8_kind*cexps(k))/fact0
        gamma_help(:,1:3+la+lb)=gamma(3+la+lb,gamma_arg2(:)*&
             rcsabc(:))
        help_arr_gr_vec=0.0_r8_kind
        step = num * 6
        i_min = 1
        i_max = step
        do i_l=0,la+lb+1
           help_vec=fact8*rcsabc(:)**i_l*((fact7-i_l)&
                *gamma_help(:,i_l+1)+gamma_arg2(:)*rcsabc(:)*&
                gamma_help(:,i_l+2))
           help_loc(1:num) = help_vec
           help_loc(num+1:2*num) = help_vec
           help_loc(2*num+1:3*num) = help_vec
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_loc)
#endif
           help_loc(3*num+1:6*num) = help_loc(1:3*num)
           do i_la=1,nm_la
              la_index=lasq+i_la
              do i_lb=1,nm_lb
                 lb_index=lbsq+i_lb
                 help_arr_gr_vec(:,i_lb,i_la)=help_arr_gr_vec(:,i_lb,i_la)+&
                      prod_arr_gr_vec(i_min:i_max,la_index,lb_index) * help_loc
              end do! loop over i_lb
           end do! loop over i_la
           i_min = i_min + step
           i_max = i_max + step
        enddo! loop over i_l
        call calc_coul_gradients(-1,k,1)
     end do! loop over fitfunctions
  end subroutine r2_coulomb
  !**************************************************************

  !**************************************************************
  subroutine l_coulomb(la,lb)
    ! calculate l-type coulomb fit integrals
    integer(kind=i4_kind),intent(in):: la,lb
    real(kind=r8_kind) :: help_loc(num*6)
    integer(kind=i4_kind) :: i_min,i_max,step

    do l=1,ly_max
       help_vec=1/aqapb(:,l)
       do m=1,2*l+1
          l2pm=l*l+m
          do k_gr=1,3
             yl_arr_xyz(:,l2pm,k_gr)=&
                  yl_arr_xyz(:,l2pm,k_gr)*help_vec
          end do
       end do
    end do

    do i_l=1,lmax_abs
       n_independent_fcts  = &
            unique_atoms(i)%symadapt_partner(1,i_l)%n_independent_fcts
       !
       allocate( &
            intermediate(num,nm_la,nm_lb,n_independent_fcts,0:la+lb),&
            intermediate_gr(num*6*(la+lb+2),nm_lb,nm_la,n_independent_fcts),&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler('LL_CALCULATE_GRADS: allocation 17 failed')
       intermediate=0.0_r8_kind
       intermediate_gr=0.0_r8_kind
       !
       ! do symmetry adaption
       independents: do i_ind=1,n_independent_fcts
          n_contributing_fcts = &
               unique_atoms(i)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%n_fcts
          coeff => unique_atoms(i)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          magn => unique_atoms(i)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          eq_atom => unique_atoms(i)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%&
               I_equal_atom
          contributing: do i_cnt=1,n_contributing_fcts
             if(eq_atom(i_cnt)==j) then
                counter=1
                help_vec0=(bexp_arr/aexp_arr)
                help_vec=1.0_r8_kind
                do i_lb=0,lb
                   help_arr0=spread(aqapb(:,i_l)*help_vec,2,laposq)
                   do mb=1,2*i_lb+1
                      diff_arr(:,:,counter)=help_arr0*&
                           diff_rule_nested(yl_arr2(:,:)&
                           ,1,laposq,counter,i_l**2+magn(i_cnt))
                      do k_gr=1,3
                         diff_arr_xyz(:,:,counter,k_gr)=help_arr0*&
                              diff_rule_nested(yl_arr_xyz(:,:,k_gr)&
                              ,1,laposq,counter,i_l**2+magn(i_cnt))
                      end do
                      counter=counter+1
                   end do
                   help_vec=help_vec*help_vec0
                end do
                call calc_intermediate(la,lb)
             end if
          end do contributing

       end do independents
       ! end of symmetry adaption

       ncexps=  unique_atoms(i)%l_ch(i_l)%n_exponents
       cexps => unique_atoms(i)%l_ch(i_l)%exponents(:)

       do k=1,ncexps ! loop over fitexponents

          rcsabc=cexps(k)/(fact0+cexps(k))
          gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2&
               (:)*rcsabc(:))
          fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/&
               (fact0+cexps(k)))
          do i_ind=1,n_independent_fcts
             help_arr_gr_vec = 0.0_r8_kind
             step = num * 6
             i_min = 1
             i_max = step
             do l=0,la+lb+1
                help_vec=(2.0_r8_kind*fact0*rcsabc)**i_l*gamma_help&
                     (:,l+1+i_l)*(rcsabc(:))**l*fact8
                help_loc(1:num) = help_vec
                help_loc(num+1:2*num) = help_vec
                help_loc(2*num+1:3*num) = help_vec
#if defined(_VECTOR) && defined(_VPP)
!OCL NOVREC(help_loc)
#endif
                help_loc(3*num+1:6*num) = help_loc(1:3*num)
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      help_arr_gr_vec(:,i_lb,i_la) = help_arr_gr_vec(:,i_lb,i_la) + &
                           intermediate_gr(i_min:i_max,i_lb,i_la,i_ind) * help_loc
                   enddo
                enddo
                i_min = i_min + step
                i_max = i_max + step
             enddo
             call calc_coul_gradients(i_l,k,i_ind)
          end do! loop over i_ind
       end do! loop over fitfunctions

       deallocate(intermediate,intermediate_gr,stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LL_CALCULATE_GRADS : deallocation 11 failed")
    end do
  end subroutine l_coulomb
  !**************************************************************

  function unity_matrix_p (m) result (yes)
    !
    ! Returns true iff  m is a 3x3 unity matrix.  In this case there
    ! is no need to transform cartesian gradients to the derivatives
    ! wrt the totally symmetric modes.
    !
    implicit none
    real(r8_kind), intent(in) :: m(:, :)
    logical :: yes
    ! *** end of interface ***

    integer, parameter :: I(3, 3) = reshape ([1, 0, 0, &
         0, 1, 0, &
         0, 0, 1], &
         [3, 3])
    yes = .false.
    if (size (m, 1) /= 3) return
    if (size (m, 2) /= 3) return

    ! FIXME: I assume it is safer to return false in case of doubts:
    yes = sum((m - I)**2) < 1.0e-7_r8_kind
  end function unity_matrix_p

end module ll_calculate_grads_module

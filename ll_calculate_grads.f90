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
subroutine ll_calculate_grads(na,nb,equalb,la,lb,imode)
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
 use ll_calculate_grads_module
 use calc3c_switches
 use shgi_cntrl, only: IPSEU
  implicit none


  !== Interrupt end of public interface of module ====================

  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
  integer(kind=i4_kind),intent(in) :: equalb ! number of equal atom b
  integer(kind=i4_kind),intent(in) :: la ! angular momentum of unique atom a
  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
  integer(kind=i8_kind),intent(in) :: imode ! for control

  integer(kind=i4_kind) alloc_status(10)
  
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

  ! was in ll_calculate_module before:
  type(unique_atom_type), pointer  :: ua_pointer

! MEMSET(0)
! write(*,'(A22,5I5)') 'll_calculate_grads:',na,nb,equalb,la,lb

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'll_grads: PP=',pseudopot_present,' imode=',imode

  nm_la=2*la+1
  nm_lb=2*lb+1
  lasq=la**2
  lbsq=lb**2
  laposq=(la+1)**2
  lbposq=(lb+1)**2

  naexps = unique_atoms(na)%l_ob(la)%n_exponents
  nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents
  ima = unique_atoms(na)%moving_atom
  imb = unique_atoms(nb)%moving_atom
  moving_a = ima > 0
  moving_b = imb > 0
  k_gr_0 = 4; if (moving_a) k_gr_0 = 1
  k_gr_1 = 3; if (moving_b) k_gr_1 = 6
  split_gradients = options_split_gradients()
  model_density = options_xcmode() == xcmode_model_density
  spin_polarized = .not. options_spin_restricted()
  lmax=max(la,lb)
  xyz_map(1)=3_i4_kind
  xyz_map(2)=4_i4_kind
  xyz_map(3)=2_i4_kind

  allocate(fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps),stat=alloc_stat)
       MEMLOG(size(fact0_arr)*4)
       ASSERT(alloc_stat.eq.0)


  xa = center1
  xb = center2
  xd =xa-xb
  aexps => unique_atoms(na)%l_ob(la)%exponents(:)
  bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)
  arg=sum(xd**2)
  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  where(fact2_arr*arg>options_integral_expmax()) ! cutoff: where almost no overlap
     cutoff=.false.              ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where
  num=count(cutoff)
  all_zero: if(num==0) then ! all integrals are equal zero

     if (integralpar_2cob_ol_grad) then
        do i_grad=1,size(prim_int_2cob_ol_grad)
           prim_int_2cob_ol_grad(i_grad)%m = 0.0_r8_kind
         ! FIXME: init prim_int_2cob_ol_dervs elsewhere:
         if(integralpar_dervs) then
          do k2dr=1,size(prim_int_2cob_ol_grad)
           prim_int_2cob_ol_dervs(i_grad,k2dr)%m = 0.0_r8_kind
          enddo
         endif
        end do
     end if

     if (integralpar_solv_grad) then
        do i_grad=1,gradient_data_n_spin_gradients
           prim_int_3cob_solv_grad(i_grad)%m = 0.0_r8_kind      !!!!!!!!!!!!!!
        enddo
        if (with_pc .and. .not.fixed_pc) then
           do i_grad=1,totsym_field_length
              prim_int_3cob_solv_grad_pc(i_grad)%m = 0.0_r8_kind      !!!!!!!!!!!!!!
           enddo
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
     end if

     MEMLOG(-size(fact0_arr)*4)
     deallocate(fact0_arr,fact1_arr,&
          fact2_arr,cutoff,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

     return
  end if all_zero

  allocate (fact0(num),fact1(num),fact2(num),fact4(num),fact5(num),&
       fact6(num),fact7(num),fact8(num),rcsabc(num),tau(num),&
       overlap(num,laposq,lbposq),&
       overlap_xyz(num,laposq,lbposq,3),&
       gamma_arg(num,3),aexp_arr(num),bexp_arr(num),&
       twoa(num),twob(num),&
       clmamb_scalar((lmax+1)**2),&
       clmamb_scalar_xyz((lmax+1)**2,3),&
       clmamb(num,laposq),&
       clmamb_xyz(num,laposq,3),&
       clmamb2(num,laposq),&
       clmamb2_xyz(num,laposq,3),&
       diff_arr(num,laposq,lbposq),&
       diff_arr_xyz(num,laposq,lbposq,3),&
       diff_arr0(laposq,lbposq),&
       diff_arr0_xyz(laposq,lbposq,3),&
       stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

       MEMLOG(size(fact0)*10)
       MEMLOG(size(fact0)*4)
       MEMLOG(size(overlap)*4+size(gamma_arg))
       MEMLOG(size(clmamb_scalar)*4+size(clmamb)*8)
       MEMLOG(size(diff_arr)*4+size(diff_arr0)*4)

  setup_grads: if (integralpar_3cob_grad) then

     allocate( nuc_grad_gr(num,nm_lb,nm_la,6),&
          help_arr_gr(num,6,nm_lb,nm_la),&
          help_arr_gr_vec(num*6,nm_lb,nm_la),&
          help_arr2(num,laposq,0:la+lb),&
          help_vec(num),help_vec0(num),&
          help_mat(num,laposq,lbposq),&
          grad_mat(nbexps,naexps,nm_lb,nm_la,6),&
          stat=alloc_stat)
          MEMLOG(size(nuc_grad_gr)+size(help_arr2))
          MEMLOG(size(help_arr_gr)+size(help_arr_gr_vec))
          MEMLOG(size(help_vec)*2+size(help_mat)+size(grad_mat))

          ASSERT(alloc_stat.eq.0)

     nuc_grad_gr=0.0_r8_kind
     grad_mat=0.0_r8_kind
     
     if (model_density .and. spin_polarized) then
        allocate( grad_mat_spin(nbexps,naexps,nm_lb,nm_la,6), stat=alloc_stat)
        if (alloc_stat /= 0) call error_handler &
             ("LL_CALCULATE_GRADS: allocation 3m failed")
        grad_mat_spin=0.0_r8_kind
     endif

  endif setup_grads

  if(integralpar_solv_grad) then  !!!!!!!!!!!!!!!!!!!!!!
     !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
     allocate( &
          solv_grad_gr(num,nm_lb,nm_la,6),&
          help_arr_gr(num,6,nm_lb,nm_la),&
          help_arr_gr1(num,3,nm_lb,nm_la),&
          help_arr_gr_vec(num*6,nm_lb,nm_la),&
          help_arr2(num,laposq,0:la+lb),&
          help_vec(num),help_vec0(num),&
          help_mat(num,laposq,lbposq),&
          grad_mat(nbexps,naexps,nm_lb,nm_la,6),&
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LL_CALCULATE_GRADS: allocation 3 failed")
     grad_mat=0.0_r8_kind
     if (model_density .and. spin_polarized) then
        allocate( grad_mat_spin(nbexps,naexps,nm_lb,nm_la,6), stat=alloc_stat)
        if (alloc_stat /= 0) call error_handler &
             ("LL_CALCULATE_GRADS: allocation 3m failed")
        grad_mat_spin=0.0_r8_kind
     endif
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  endif                             !!!!!!!!!!!!!!!!!!!!!!!

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  ! fact7= 1/sqrt(a**l*(2l-1)!!)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)
  fact0_arr=spread(aexps,1,nbexps)
  aexp_arr=pack(fact0_arr,cutoff)
  fact0_arr=spread(bexps,2,naexps)
  bexp_arr=pack(fact0_arr,cutoff)
  twoa=-2.0_r8_kind*aexp_arr
  twob=-2.0_r8_kind*bexp_arr

  MEMLOG(-size(fact0_arr)*3)
  deallocate(fact0_arr,fact1_arr,fact2_arr,stat=alloc_stat)
  ASSERT(alloc_stat.eq.0)

  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  gamma_arg(:,1)=(pack(spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps),cutoff))/fact0

  gamma_arg(:,2)=(pack(spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps),cutoff))/fact0

  gamma_arg(:,3)=(pack(spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps),cutoff))/fact0

  ! precalculation of solid harmonics
   clmamb_scalar=solid_harmonics_scalar(max(la,lb),xd)
   fact4=1.0_r8_kind
  counter=1

  tau=fact2*arg ! a*b/(a+b)*(A-B)**2

  ! calculate derivatives of solid harmonics
  clmamb_scalar_xyz=0.0_r8_kind

  do l=1,(lmax+1)**2
     do k_gr=1,3
        do k=1,solhrules_differential(xyz_map(k_gr),l)%n_summands
           clmamb_scalar_xyz(l,k_gr)=clmamb_scalar_xyz(l,k_gr)+&
                solhrules_differential(xyz_map(k_gr),l)%coef(k)*&
                clmamb_scalar(solhrules_differential(xyz_map(k_gr),l)%lm_sh(k))   
        end do
     end do
  end do!loop over l  

  fact4=1.0_r8_kind
  counter=1
  do l=0,la
     do m=1,2*l+1
        clmamb(:,counter)=clmamb_scalar(counter)*fact4
        clmamb2(:,counter)= &
            clmamb_scalar(counter)*fact4*(tau-real(l,kind=r8_kind))
        do k_gr=1,3

           clmamb_xyz(:,counter,k_gr)=clmamb_scalar_xyz(counter,k_gr)*fact4

           clmamb2_xyz(:,counter,k_gr)=clmamb_scalar_xyz(counter,k_gr)&
                                            *fact4*(tau-real(l,kind=r8_kind))&
                +fact4*clmamb_scalar(counter)*fact2*xd(k_gr)*2.0_r8_kind
                 !                            --------------------------  d(tau)/da
        end do
        counter=counter+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo

  ! PROCESS <xi_i|xi_j> and <xi_i|T|xi_j>
  call calc_2center(la,lb,equalb)

  ! PROCESS <xi_i|V_nuc|xi_j>, <xi_i|V_pvsp|xi_j>, and <xi_i|V_H|xi_j>

    if(integralpar_3cob_grad) then

       fact8=2.0_r8_kind*sqrt(fact0/pi)

       !  used for atoms and PC
       allocate(prod_arr(num,laposq,lbposq,0:la+lb),&
            prod_arr_gr(num,6,laposq,lbposq,0:la+lb+1),&
            prod_arr_gr_vec(num*6*(la+lb+2),laposq,lbposq),&
            stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       MEMLOG(size(prod_arr)+size(prod_arr_gr)+size(prod_arr_gr_vec))


       unique_atom_loop: do i=1,n_unique_atoms  ! loop over third center
          lmax_ch= unique_atoms(i)%lmax_ch      ! maximum l  for chargefit
          lmax_xc= unique_atoms(i)%lmax_xc      ! maximum l  for xcfit  
          ! determine the maximal angular momentum
          lmax_abs=lmax_ch
          ly_max=max(la,lb,lmax_ch)
          max_order=max(1+la+lb+lmax_abs,3+la+lb)+1
          z= unique_atoms(i)%z                  ! charge 
          zc=unique_atoms(i)%zc
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'll_grads: ua=',i,', zero its charge!'
        zc = zero
        z  = zero
          n_equals=unique_atoms(i)%n_equal_atoms
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
               ("LL_CALCULATE_GRADS : allocation 8 failed")
          aqapb(:,0)=1.0_r8_kind
          aqapb(:,1)=aexp_arr/fact0
          do i_l=2,ly_max
             aqapb(:,i_l)=aqapb(:,i_l-1)*aqapb(:,1)
          end do
          imc = unique_atoms(i)%moving_atom
          moving_c = imc > 0
          if (moving_c)then
             grad_dim = gradient_index(imc+1) - gradient_index(imc)
          else
             grad_dim = 0
          endif


          ! --- further allocation ----------------------------------
          ! num : number of pairs(a,b) which are inside the cutoff
          ! for s-and r2-type there is only 1 indep. fct

          allocate(coul_int_grad_totsym(grad_dim),&
               nuc_grad(num,nm_lb,nm_la,grad_dim),stat=alloc_stat)   ! Anguang
          ASSERT(alloc_stat.eq.0)
          nuc_grad=0.0_r8_kind

        if(.not.new_3c_co_grad) then
          do k_gr=1,6
             allocate (coul_int_grad(k_gr)%l(-1:lmax_ch),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate( coul_int_grad(k_gr)%l(-1)%m(num,ncexps,1,nm_lb,nm_la),&
                  stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)

             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_grad(k_gr)%l(0)%m(num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             coul_int_grad(k_gr)%l(0)%m=0.0_r8_kind
             coul_int_grad(k_gr)%l(-1)%m=0.0_r8_kind
          end do

          do i_grad=1,grad_dim ! only if moving_c
             allocate(coul_int_grad_totsym(i_grad)%l(-1:lmax_ch),&
                  stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_grad_totsym(i_grad)%l(-1)%m&
                  (num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_grad_totsym(i_grad)%l(0)%m&
                  (num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             coul_int_grad_totsym(i_grad)%l(0)%m=0.0_r8_kind
             coul_int_grad_totsym(i_grad)%l(-1)%m=0.0_r8_kind
          end do

          do i_l=1,lmax_ch
             ! allocation for coul_int
             n_independent_fcts  = &
                  unique_atoms(i)%symadapt_partner(1,i_l)%n_independent_fcts
             ncexps=  unique_atoms(i)%l_ch(i_l)%n_exponents
             cexps => unique_atoms(i)%l_ch(i_l)%exponents(:)
             do k_gr=1,6
                allocate(coul_int_grad(k_gr)%l(i_l)%m(num,ncexps,n_independent_fcts,nm_lb,nm_la),&
                         stat=alloc_stat)
                ASSERT(alloc_stat.eq.0)
                coul_int_grad(k_gr)%l(i_l)%m=0.0_r8_kind
             end do
             call integral_interrupt_2cob3c()
             do i_grad=1,grad_dim ! only if moving_c
                allocate(coul_int_grad_totsym(i_grad)%l(i_l)%m(num,ncexps,n_independent_fcts,nm_lb,nm_la),&
                     stat=alloc_stat)
                coul_int_grad_totsym(i_grad)%l(i_l)%m=0.0_r8_kind
                ASSERT(alloc_stat.eq.0)
             enddo
          end do
         endif

          ! do a precalculation of a factor needed for the
          ! diff rule
          counter=1
          fact4=1.0_r8_kind
          do i_l=0,ly_max
             do ma=1,2*i_l+1
                fact10(:,counter)=fact4
                counter=counter+1
             enddo
             fact4=fact4*aqapb(:,1)
          enddo

          equal_atoms: do j=1,n_equals 

             yl_arr_xyz=0.0_r8_kind
             call integral_interrupt_2cob3c()
             if(moving_c)then
                if(grad_dim==3) then
                   if(sum((unique_atom_grad_info(imc)%m(:,:,j)-unity_matrix)**2)<&
                        1.0e-7_r8_kind) then
                      do_rotation=.false.
                   else
                      do_rotation=.true.
                      rotmat=>unique_atom_grad_info(imc)%m(:,:,j)
                   endif
                else
                   do_rotation=.true.
                   rotmat=>unique_atom_grad_info(imc)%m(:,:,j)
                end if
             endif
             xc=unique_atoms(i)%position(:,j)
             yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                  spread(xc,1,num))
             gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                  (gamma_arg(:,3)-xc(3))**2)*fact0
             ! first calculation of the nuclear attraction
             gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2(:))
             fact8=2.0_r8_kind*sqrt(fact0/pi)
             counter=1
             yl_arr2(:,:)=yl_arr(:,:)/fact10(:,1:((ly_max+1)**2))
             ! now calculation of derivatives of yl_arr with respect to -c
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
                   ! now let us differentiate diff_arr
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
              call calc_prod_arr(la,lb)
            endif



             ! PROCESS <xi_i|Zc/|r-Rc||xi_j>
             iam_ppot=zc/=0.0_r8_kind &
                  .and.(.not.operations_core_density).and. &
                  pseudopot_present.and.integralpar_relativistic

!!$             call calc_nuc_gradients(third_center_required=.true.) 
             call calc_nuc_help_arr(la,lb)

             if(iam_ppot) then
                call calc_final_nuc_and_pvsp_grad(pseudo_grad_gr,nuc_grad) 
!!$                print*,'pseudo_grad_gr',sum(pseudo_grad_gr)
             else
                call calc_final_nuc_and_pvsp_grad(nuc_grad_gr,nuc_grad)
!!$                print*,'nuc_grad_gr',sum(nuc_grad_gr)
             end if
    
!!$        print*,'calc_nuc_gradients done pseudo_grad_gr',sum(pseudo_grad_gr)
             ! pseudo atom
             ! for PP and rel use other store pseudo_grad_gr and nuc_grad

!!$	if(ewpc_n.ne.0.and..false.) then
!!$        do mb=1,2*lb+1
!!$           do i_l =0,lb+la
!!$	     cluster_epe(:,:,mb)=cluster_epe(:,:,mb)+z*prod_arr(:,la**2+1:&
!!$                     (la+1)**2,lb**2+mb,i_l)*spread&
!!$                     (fact8*gamma_help(:,i_l+1),2,2*la+1)
!!$             end do
!!$          end do
!!$	endif

             ! PROCESS [xi_i|f_k|xi_j]
             ! first evaluate s-type coulomb integrals
             if(old_3c_co) call s_coulomb(la,lb)

             ! now treating r2-type coloumb integrals
             if(old_3c_co) call r2_coulomb(la,lb)

             ! finally calculate l-type coulomb integrals
             if(old_3c_co) call l_coulomb(la,lb)


          end do equal_atoms


          par_3cob_grad: if(integralpar_3cob_grad) then
             ! coulombic contributions of derivatives with respect
             ! to center a and b
             ! CONTRACT d/dRa [xi_i|f_k|xi_j] and d/dRb [xi_i|f_k|xi_j]

	 if(.false..and.old_3c_co) then
	     do i_l=-1,0 !!!! lmax_ch
		if(associated(coul_int_grad(3)%l(i_l)%m)) &
	        print*,sum(coul_int_grad(3)%l(i_l)%m),i_l,'ga ll'
		if(associated(coul_int_grad(6)%l(i_l)%m)) &
	        print*,sum(coul_int_grad(6)%l(i_l)%m),i_l,'gb ll'
	     enddo
         endif

             do k_gr=k_gr_0,k_gr_1
                grad_mat_p => grad_mat(:,:,:,:,k_gr)
                if (model_density .and. spin_polarized) then
                 if(old_3c_co) then
                   grad_mat_sp => grad_mat_spin(:,:,:,:,k_gr)
		   if(old_3c_fc) &
                   call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr),&
                        grad_mat_p,grad_mat_sp)
		 endif
                else ! standard SCF
	        if(old_3c_co) then
                grad_mat_p=>grad_mat(:,:,:,:,k_gr)
		if(old_3c_fc) &
                call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr),&
                     grad_mat_p)
		endif
                end if
             end do

           if(.not.new_3c_co_grad) then
             do k_gr=1,6
                do i_l = -1, lmax_ch
                   deallocate(coul_int_grad(k_gr)%l(i_l)%m,STAT=alloc_stat)
                   ASSERT(alloc_stat.eq.0)
                end do
                deallocate (coul_int_grad(k_gr)%l, STAT=alloc_stat)
                ASSERT(alloc_stat.eq.0)
             end do
           endif

             ! coulombic contributions of derivatives with respect
             ! to center c
             ! CONTRACT dsym/dRc [xi_i|f_k|xi_j]
         do i_grad=1,grad_dim ! only if moving_c
                index = gradient_index(imc) + i_grad - 1
                if (model_density) then
                   if (spin_polarized) then
                      spin_index = index + gradient_data_n_gradients
                      if (split_gradients) then
			if(old_3c_co.and.old_3c_fc) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad), &
                              prim_int_3cob_coul_grad(index)%m, & ! V_H
                              prim_int_3cob_grad(spin_index)%m, & ! V_X,spin
                              prim_int_3cob_grad(index)%m )       ! V_X,tot
                      else ! total_gradient only
			if(old_3c_co.and.old_3c_fc) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad), &
                              prim_int_3cob_grad(index)%m, &      ! V_H+V_X,tot
                              prim_int_3cob_grad(spin_index)%m )  ! V_X,spin
                      endif
                   else ! spin_restricted
                      if (split_gradients) then
			if(old_3c_co.and.old_3c_fc) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad), &
                              prim_int_3cob_coul_grad(index)%m, & ! V_H
                              mda_xcpot_gradients = &
                              prim_int_3cob_grad(index)%m )       ! V_X,tot
                      else ! total_gradient only
			if(old_3c_co.and.old_3c_fc) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad), &
                              prim_int_3cob_grad(index)%m )       ! V_H+V_X,tot
                      endif
                   endif
            else ! standard SCF

	     if(old_3c_co) then

!	     do i_l=-1,-1 !!!! lmax_ch
!		if(associated(coul_int_grad_totsym(i_grad)%l(i_l)%m)) &
!	        print*,sum(coul_int_grad_totsym(i_grad)%l(i_l)%m),i_l,'gc ll'
!	     enddo

		   if(old_3c_fc) &
                   call fitcontract('grad',num,i,cutoff, &
                        coul_int_grad_totsym(i_grad), &
                        prim_int_3cob_grad(index)%m )
	     endif
            endif

           if(.not.new_3c_co_grad) then
                do i_l = -1, lmax_ch
                   deallocate(coul_int_grad_totsym(i_grad)%l(i_l)%m,&
                        STAT=alloc_stat)
                   ASSERT(alloc_stat.eq.0)
                enddo
                deallocate (coul_int_grad_totsym(i_grad)%l,STAT=alloc_stat)
                ASSERT(alloc_stat.eq.0)
            endif
           end do

             deallocate(coul_int_grad_totsym,stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)

             ! nuclear charge contributions of derivatives with respect
             ! to center c
             ! ADD dsym/dRc <xi_i|V_nuc[c]|xi_j>, dsym/dRc <xi_i|V_pvsp[c]|xi_j>
             if (moving_c) index = gradient_index(imc) - 1
             do i_grad=1,grad_dim ! only if moving_c
                counter=1
                if(.not.integralpar_relativistic) then

                   if (split_gradients) then
                      grad_mat_p=>prim_int_3cob_nuc_grad(index+i_grad)%m
                   else
                      grad_mat_p=>prim_int_3cob_grad(index+i_grad)%m

                   endif

                else ! integralpar_relativistic
                   if(iam_ppot) then
                      grad_mat_p=>prim_int_2cob_pseudo_grad(index+i_grad)%m

                   else                           
                      grad_mat_p=>prim_int_2cob_nuc_grad(index+i_grad)%m                
                   end if

                end if
             
                   
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      grad_mat_p(:,:,i_lb,i_la) = grad_mat_p(:,:,i_lb,i_la)&
                           -unpack(nuc_grad(:,i_lb,i_la,i_grad),cutoff,zero)
                   end do
                end do
             end do
          endif par_3cob_grad

          deallocate(nuc_grad,yl_arr,yl_arr_xyz,yl_arr2,help_arr,&
               help_arr0,aqapb,gamma_help,gamma_arg2,fact10,&
               stat=alloc_stat)
          if (alloc_stat/=0) call error_handler &
               ("LL_CALCULATE_GRADS : deallocation 17 failed")
       end do unique_atom_loop

       ppot:if (pseudopot_present.and.(.not.operations_core_density)) then

          unqiue_pseudopot_loop: do i=1,n_unique_atoms+n_timps
             if(i<=n_unique_atoms) then
                imc = unique_atoms(i)%moving_atom
                moving_c = imc > 0
                ua_pointer=>unique_atoms(i)
                if (moving_c)then
                   grad_dim = gradient_index(imc+1) - gradient_index(imc)
                else
                   grad_dim = 0
                endif
             else
                ua_pointer=>unique_timps(i-n_unique_atoms)

!!$		moving_c=.false.
!!$		grad_dim = 0
		imc = unique_timps(i-n_moving_unique_atoms)%moving_atom
		moving_c = imc > 0
                if (moving_c)then
                   grad_dim = gradient_index(N_moving_unique_atoms+imc+1) - &
                        gradient_index(N_moving_unique_atoms+imc)
		else
                   grad_dim = 0
		endif
             end if
             
                 zc=ua_pointer%zc
          pp: if (zc/=0.0_r8_kind) then
            ABORT('not supported')
          endif pp! end of zc/=0.0_r8_kind
       enddo unqiue_pseudopot_loop
   
    endif ppot
       unique_charge_loop: do i=1,pointcharge_N+n_timps
          if(i<=n_timps) then
             z= unique_timps(i)%z -  unique_timps(i)%zc                ! charge 
             n_equals=unique_timps(i)%n_equal_atoms
             imc = unique_timps(i)%moving_atom
             moving_c = imc > 0
             if (moving_c)then
                grad_dim = gradient_index(N_moving_unique_atoms+imc+1) - &
                     gradient_index(N_moving_unique_atoms+imc)
                allocate( nuc_grad(num,nm_lb,nm_la,grad_dim),stat=alloc_stat)  
                if (alloc_stat/=0) call error_handler &
                     ("LL_CALCULATE_GRADS : allocation failed nuc_grad for timps")
                nuc_grad=0.0_r8_kind
             else
                grad_dim = 0
             endif
          else
	     cycle unique_charge_loop ! pointcharges go to SHGI !!!!!!!!!!!AS
             z= pointcharge_array(i-n_timps)%z                  ! charge 
             n_equals=pointcharge_array(i-n_timps)%n_equal_charges
             moving_c=.false.
          end if
         
          ly_max=max(la,lb)
          max_order=la+lb+4
          allocate ( &
               gamma_help(num,max_order), &
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
               ("LL_CALCULATE_GRADS : allocation pc failed")
          aqapb(:,0)=1.0_r8_kind
          aqapb(:,1)=aexp_arr/fact0
          do i_l=2,ly_max
             aqapb(:,i_l)=aqapb(:,i_l-1)*aqapb(:,1)
          end do

          ! do a precalculation of a factor needed for the diff rule
          counter=1
          fact4=1.0_r8_kind
          do i_l=0,ly_max
             do ma=1,2*i_l+1
                fact10(:,counter)=fact4
                counter=counter+1
             enddo
             fact4=fact4*aqapb(:,1)
          enddo

          equal_charges: do j=1,n_equals 


             if(i<=n_timps) then
                xc=unique_timps(i)%position(:,j)
                if (moving_c) then                
                   if(grad_dim==3) then
                      if(sum((unique_timp_grad_info(imc)%m(:,:,j)-unity_matrix)**2)<&
                           1.0e-7_r8_kind) then
                         do_rotation=.false.
                      else
                         do_rotation=.true.
                         rotmat=>unique_timp_grad_info(imc)%m(:,:,j)
                      endif
                   else
                      do_rotation=.true.
                      rotmat=>unique_timp_grad_info(imc)%m(:,:,j)
                   end if
                end if
                

             else
                xc=pointcharge_array(i-n_timps)%position(:,j)
             end if
             yl_arr_xyz=0.0_r8_kind

             call integral_interrupt_2cob3c()
             yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                  spread(xc,1,num))
             gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                  (gamma_arg(:,3)-xc(3))**2)*fact0
             ! first calculation of the nuclear attraction
             gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2(:))
             fact8=2.0_r8_kind*sqrt(fact0/pi)
             counter=1
             yl_arr2(:,:)=yl_arr(:,:)/fact10(:,1:((ly_max+1)**2))
             ! now calculation of derivatives of yl_arr with respect to -c
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
                   ! now let us differentiate diff_arr
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
              call calc_prod_arr(la,lb)
            endif


             ! PROCESS <xi_i|Zp/|r-Rp||xi_j>
             iam_ppot=pseudopot_present.and.integralpar_relativistic

!!$             call calc_nuc_gradients(third_center_required=moving_c)

             call calc_nuc_help_arr(la,lb)
             if ( moving_c ) then
                if(iam_ppot) then
                   call calc_final_nuc_and_pvsp_grad(pseudo_grad_gr,nuc_grad) 
!!$                   print*,'pseudo_grad_gr',sum(pseudo_grad_gr)
                else
                   call calc_final_nuc_and_pvsp_grad(nuc_grad_gr,nuc_grad)
!!$                   print*,'nuc_grad_gr',sum(nuc_grad_gr)
                end if

             else
                if(iam_ppot) then
                   call calc_final_nuc_and_pvsp_grad(pseudo_grad_gr)
                else
                   call calc_final_nuc_and_pvsp_grad(nuc_grad_gr)
                end if

             endif

          enddo equal_charges

          if (moving_c) then
             index = gradient_index(imc+n_moving_unique_atoms) - 1
             do i_grad=1,grad_dim ! only if moving_c
             grad_mat_p=>prim_int_3cob_grad(index+i_grad)%m
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      grad_mat_p(:,:,i_lb,i_la)=grad_mat_p(:,:,i_lb,i_la)&
                           -unpack(nuc_grad(:,i_lb,i_la,i_grad),cutoff,zero)             
                   end do
                end do
             end do
          end if
          
          deallocate ( &
               gamma_help, &
               gamma_arg2,&
               yl_arr,&
               yl_arr_xyz,&
               help_arr,&
               help_arr0,&
               fact10,&
               yl_arr2,&
               aqapb,&
               stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
 
          if(moving_c) then
             deallocate(nuc_grad,stat=alloc_stat)
             if (alloc_stat/=0) call error_handler &
               ("LL_CALCULATE_GRADS : deallocation of nuc_grad for timps failed")
          end if
          
       enddo unique_charge_loop
     endif

     if(old_solv_grad.and.integralpar_solv_grad) then
     allocate(&
         prod_arr(num,laposq,lbposq,0:la+lb),&
         prod_arr_gr(num,6,laposq,lbposq,0:la+lb+1),&
         prod_arr_gr_vec(num*6*(la+lb+2),laposq,lbposq),&
         stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)

     call calc_solv_grad(la,lb)
     endif

  if (integralpar_3cob_grad) then

     if((.not.integralpar_relativistic)) then
        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
        
        do k_gr=k_gr_0,k_gr_1
           do ma=1,2*la+1
              do mb=1,2*lb+1
                 grad_mat(:,:,mb,ma,k_gr)=grad_mat(:,:,mb,ma,k_gr)-&
                      unpack(nuc_grad_gr(:,mb,ma,k_gr),cutoff,zero)
              enddo
           enddo
        end do
        MEMLOG(-size(nuc_grad_gr)-size(help_arr2))
        deallocate(nuc_grad_gr,help_arr2,stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)

     end if! .not.integralpar_relativistic    

                  ! LOAD prim_int_3cob_grad(...,ia:) and prim_int_3cob_grad(...,ib: ) or
                  !      prim_int_2cob_ks_grad( 1:N) and prim_int_2cob_ks_grad(N+1:M)

     if (split_gradients) then
       pointer_prim_int=>prim_int_2cob_ks_grad
       call add_nuc_or_pvsp(equalb,grad_mat,compact_storage=.TRUE.) 
       if (model_density .and. spin_polarized) then
          spin_index = size(prim_int_2cob_ks_grad)/2 + 1
          pointer_prim_int => prim_int_2cob_ks_grad(spin_index:)
          call add_nuc_or_pvsp(equalb,grad_mat_spin,compact_storage=.TRUE.)
       endif

     else ! not split_gradients

       pointer_prim_int=>prim_int_3cob_grad
       call add_nuc_or_pvsp(equalb,grad_mat) !!!!!!!!!! this is not empty for relativistic case	

       if (model_density .and. spin_polarized) then
          spin_index = gradient_data_n_gradients + 1
          pointer_prim_int => prim_int_3cob_grad(spin_index:)
          call add_nuc_or_pvsp(equalb,grad_mat_spin)
       endif
     endif


     relmap: if(integralpar_relativistic) then

        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
        grad_mat=0.0_r8_kind
        do k_gr=k_gr_0,k_gr_1
           do ma=1,2*la+1
              do mb=1,2*lb+1
                 grad_mat(:,:,mb,ma,k_gr)=grad_mat(:,:,mb,ma,k_gr)-&
                      unpack(nuc_grad_gr(:,mb,ma,k_gr),cutoff,zero)
              enddo
           enddo
        end do

        MEMLOG(-size(nuc_grad_gr)-size(help_arr2))
        deallocate(nuc_grad_gr,help_arr2,stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)

        ! LOAD prim_int_2cob_nuc_grad(...,ia:), prim_int_2cob_nuc_grad(...,ib:)
        pointer_prim_int=>prim_int_2cob_nuc_grad
        call add_nuc_or_pvsp(equalb,grad_mat)

     end if relmap

  endif

  if(integralpar_solv_grad) then !!!!!!!!!!!!!!!!!!!!!!
     !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    if(old_solv_grad) then
     do k_gr=k_gr_0,k_gr_1
        do ma=1,2*la+1
           do mb=1,2*lb+1
              grad_mat(:,:,mb,ma,k_gr)=grad_mat(:,:,mb,ma,k_gr)-&
                   unpack(solv_grad_gr(:,mb,ma,k_gr),cutoff,zero)
           enddo
        enddo
     end do

     pointer_prim_int=>prim_int_3cob_solv_grad
     call add_nuc_or_pvsp(equalb,grad_mat)
     if (model_density .and. spin_polarized) then
        spin_index = gradient_data_n_gradients + 1
        pointer_prim_int => prim_int_3cob_solv_grad(spin_index:)
        call add_nuc_or_pvsp(equalb,grad_mat_spin)
     endif
    endif

     deallocate(solv_grad_gr,help_arr2,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ("LL_CALCULATE_GRADS: deallocation solv_grad_gr ... failed")   
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  endif                          !!!!!!!!!!!!!!!!!!!!!!

  MEMLOG(-size(grad_mat))
  deallocate(grad_mat,stat=alloc_stat)
  ASSERT(alloc_stat.eq.0)

  if (model_density .and. spin_polarized) then
     deallocate(grad_mat_spin,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ("LL_CALCULATE_GRADS: deallocation failed 5m grad_mat_spin")
  endif

  MEMLOG(-size(clmamb)*4)
  MEMLOG(-size(fact0)*8)
  deallocate(clmamb,clmamb_xyz,fact0,&
       fact2,fact4,fact5,fact6,fact7,fact8,rcsabc,&
       stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

  MEMLOG(-size(overlap)*4-size(gamma_arg))
  MEMLOG(-size(diff_arr0)*4-size(diff_arr)*4)
  MEMLOG(-size(help_vec)-size(help_vec0)-size(help_mat))
  MEMLOG(-size(help_arr_gr)-size(help_arr_gr_vec))
  deallocate(overlap,overlap_xyz,gamma_arg,diff_arr,diff_arr_xyz,diff_arr0,&
       diff_arr0_xyz,help_vec,help_vec0,&
       help_mat,help_arr_gr,help_arr_gr_vec,stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)


  if(integralpar_3cob_grad) then
  MEMLOG(-size(prod_arr)-size(prod_arr_gr)-size(prod_arr_gr_vec))
  deallocate(prod_arr,prod_arr_gr,prod_arr_gr_vec,stat=alloc_stat)
  ASSERT(alloc_stat.eq.0)
  endif

  if(integralpar_solv_grad) then
     deallocate(help_arr_gr1,stat=alloc_stat)  !!!!!!!!!!!!!!!!!
     if(alloc_stat/=0) call error_handler &    !!!!!!!!!!!!!!!!!
          ("LL_CALCULATE_GRADS: deallocation help_arr_gr1 failed") !!!!!!!!!
  endif
  if(old_solv_grad.and.integralpar_solv_grad) then
  deallocate( prod_arr,prod_arr_gr,prod_arr_gr_vec,stat=alloc_stat)
  if(alloc_stat/=0) call error_handler &
       ("LL_CALCULATE_GRADS: deallocation 7 failed")
  endif

  MEMLOG(-size(cutoff))
  MEMLOG(-size(aexp_arr)*4)
  deallocate(aexp_arr,bexp_arr,twoa,twob,cutoff,stat=alloc_stat)
  ASSERT(alloc_stat.eq.0)

!  MEMSET(0)
end subroutine ll_calculate_grads

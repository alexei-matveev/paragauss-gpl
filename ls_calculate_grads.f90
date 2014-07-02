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
subroutine ls_calculate_grads(na,nb,equalb_,la_in,lb,imode)
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

# include "def.h"
  use type_module
  use ls_calculate_grads_module
  use shgi_cntrl, only: IPSEU
  implicit none

  interface
     subroutine ls_grad_pseudo(i)
       use type_module, only: i4_kind
       implicit none
       integer(i4_kind), intent(in) :: i
     end subroutine ls_grad_pseudo
  end interface

  !== Interrupt end of public interface of module ====================

  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
  integer(kind=i4_kind),intent(in) :: la_in ! angular momentum of unique atom a
  integer(kind=i4_kind),intent(in) :: equalb_ ! number of equal atom b
  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
  integer(kind=i8_kind),intent(in) :: imode ! for control
  ! *** end of interface ***

  integer(i4_kind) :: i,j,l,k_gr,i_l,i_la,i_lb, k,k1,n,n1
  integer(i4_kind) :: la2pm, N_pc
  integer(i4_kind) :: k2dr

  ! was in ls_calculate_module before:
  type(unique_atom_type), pointer  :: ua_pointer

  intrinsic max

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'ls_grads: PP=',pseudopot_present,' imode=',imode

  ! set global in module:
  equalb = equalb_

  nm_la=2*la_in+1
  nm_lb=2*lb+1
  la_org=la_in
  lb_org=lb

  naexps = unique_atoms(na)%l_ob(la_in)%n_exponents
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
  xyz_map(1)=3_i4_kind
  xyz_map(2)=4_i4_kind
  xyz_map(3)=2_i4_kind

  ! set global LA:
  if(lb>la_in) then
     laltlb=.true.
     lagtlb=.false.
     la=lb
  else
     laltlb=.false.
     lagtlb=.true.
     la=la_in
  end if

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       diff_arr0((la+1)**2,3), &
       stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

  xa = center1  ! from int_data_module
  xb = center2  ! from int_data_module

  xd =xa-xb
  aexps => unique_atoms(na)%l_ob(la_in)%exponents(:)
  bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)

  arg=sum(xd**2)

  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  where(fact2_arr*arg>=options_integral_expmax()) ! cutoff: where almost no overlap
     cutoff=.false.              ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)

  all_zero: if(num==0) then    ! all integrals are equal zero

     if (integralpar_2cob_ol_grad) then
        do i_grad=1,size(prim_int_2cob_ol_grad)
           prim_int_2cob_ol_grad(i_grad)%m = 0.0_r8_kind
        if(integralpar_dervs) then
         ! FIXME: init prim_int_2cob_ol_dervs elsewhere
         do k2dr=1,size(prim_int_2cob_ol_grad)
           prim_int_2cob_ol_dervs(i_grad,k2dr)%m = 0.0_r8_kind
         enddo
        endif
        end do
     end if

     if (integralpar_solv_grad) then
        do i_grad=1,gradient_data_n_spin_gradients
           prim_int_3cob_solv_grad(i_grad)%m = 0.0_r8_kind      !!!!!!!!!!!!!!!
        enddo
        if(with_pc .and. .not.fixed_pc) then
        do i_grad=1,totsym_field_length
           prim_int_3cob_solv_grad_pc(i_grad)%m = 0.0_r8_kind      !!!!!!!!!!!!!!!
        enddo
        end if
     endif

     if (integralpar_3cob_grad) then

        do i_grad=1,gradient_data_n_spin_gradients

           prim_int_3cob_grad(i_grad)%m = 0.0_r8_kind

         if(integralpar_dervs) then
          ! FIXME: init prim_int_coul_dervs elsewhere
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

     deallocate(fact0_arr,fact1_arr,&
          fact2_arr,cutoff,diff_arr0,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

     return
  end if all_zero

  allocate ( &
       fact0(num), fact1(num), fact2(num), fact4(num), &
       fact5(num), fact6(num), fact7(num), fact8(num), &
       rcsabc(num), tau(num), rsaba(num), &
       aexp_arr(num), bexp_arr(num), &
       overlap(num,(la+1)**2), &
       overlap_grad(num,(la+1)**2,3), &
       kin_grad(num,2*la+1,3), &
       gamma_arg(num,3),&
       clmamb_scalar((la+1)**2), &
       clmamb(num,(la+1)**2), &
       diff_arr(num,(la+1)**2),&
       diff_arr_gr(num,(la+1)**2,3),&
       help_vec(num), &
       help_mat(num,2*la+1), &
       stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)


  if (integralpar_3cob_grad) then
     allocate( &
          kinetic(num), &
          nuc_grad_gr(num,2*la+1,6), &
          grad_mat_gr(nbexps,naexps,nm_lb,nm_la,6),&
          help_arr_gr(num,6,nm_lb,nm_la), &
          help_arr(num,(la+1)**2), &
          help_arr2(num,(la+1)**2,la+1),&
          stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
     nuc_grad_gr=0.0_r8_kind
     grad_mat_gr=0.0_r8_kind

     if (pseudopot_present) then
        allocate( grad_nuc_pseudo_ab(num,2*la+1,6), &
                  stat=alloc_stat)
        if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE_GRADS: allocation of pseudopotential failed")
        grad_nuc_pseudo_ab=0.0_r8_kind
        if(integralpar_relativistic) then
          allocate(pseudo_grad_gr(num,2*la+1,6),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          pseudo_grad_gr=0.0_r8_kind
        end if
     
     endif
     if (model_density .and. spin_polarized) then
        allocate( grad_mat_spgr(nbexps,naexps,nm_lb,nm_la,6), stat=alloc_stat)
        if (alloc_stat /= 0) call error_handler &
             ("LS_CALCULATE_GRADS: allocation 3m failed")
        grad_mat_spgr=0.0_r8_kind
     endif

  end if

  if (integralpar_solv_grad.and.old_solv_grad) then
     allocate( &
          solv_grad_gr(num,2*la+1,6), &
          grad_mat_gr(nbexps,naexps,nm_lb,nm_la,6),&
          help_arr_gr(num,6,nm_lb,nm_la), &
          help_arr_gr1(num,3,nm_lb,nm_la), &
          help_arr(num,(la+1)**2), &
          help_arr2(num,(la+1)**2,la+1),&
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE_GRADS: allocation 3 failed(s)")
     grad_mat_gr=0.0_r8_kind
     if (model_density .and. spin_polarized) then
        allocate( grad_mat_spgr(nbexps,naexps,nm_lb,nm_la,6), stat=alloc_stat)
        if (alloc_stat /= 0) call error_handler &
             ("LS_CALCULATE_GRADS: allocation 3m failed(s)")
        grad_mat_spgr=0.0_r8_kind
     endif
  endif

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  ! rsaba = (a+b)/a
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)
  
  if(lagtlb) then
     aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
     bexp_arr=pack(spread(bexps,2,naexps),cutoff)
  else
     aexp_arr=pack(spread(bexps,2,naexps),cutoff)
     bexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  end if
  rsaba=fact0/aexp_arr
  deallocate(fact0_arr,fact1_arr,fact2_arr,stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("LS_CALCULATE_GRADS: deallocation 1/2a failed")

  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  gamma_arg(:,1)=(pack(spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps),cutoff))/fact0

  gamma_arg(:,2)=(pack(spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps),cutoff))/fact0

  gamma_arg(:,3)=(pack(spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps),cutoff))/fact0

  if(laltlb) then
     xd=-xd
  end if


  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(la,xd)

  fact4=1.0_r8_kind
  counter=1
  do i=0,la
     do m=1,2*i+1
        clmamb(:,counter)=fact4
        counter=counter+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo


  ! first calculating 2-center integrals----------------
  tau=fact2*arg
  fact5=fact2*(3.0_r8_kind-2.0_r8_kind*tau+2*la)  ! a*b/(a+b)(3-2*tau+2*l)

  fact6=1.0_r8_kind/sqrt(aexp_arr**la*dfac(la))*exp(-tau)*&
       (4.0_r8_kind*fact2/fact0)**0.75_r8_kind ! 

  ! fact2=2.0_r8_kind*fact6*sqrt(fact0/pi)

  diff_arr0=0.0_r8_kind
  do l=1,(la+1)**2
     do k=1,solhrules_differential(3,l)%n_summands

        diff_arr0(l,1)=diff_arr0(l,1)+&
             solhrules_differential(3,l)%coef(k)*&
             clmamb_scalar(solhrules_differential(3,l)%lm_sh(k))   
     end do

     do k=1,solhrules_differential(4,l)%n_summands
        diff_arr0(l,2)=diff_arr0(l,2)+&
             solhrules_differential(4,l)%coef(k)*&
             clmamb_scalar(solhrules_differential(4,l)%lm_sh(k))   
     end do

     do k=1,solhrules_differential(2,l)%n_summands
        diff_arr0(l,3)=diff_arr0(l,3)+&
             solhrules_differential(2,l)%coef(k)*&
             clmamb_scalar(solhrules_differential(2,l)%lm_sh(k))   
     end do

  end do!loop over l

  ! PROCESS <xi_i|xi_j> and <xi_i|T|xi_j>
  standard_order: if(lagtlb) then

     standard_ol_grad: if (integralpar_2cob_ol_grad) then

        do l=0,la
           magnetic_number: do m=1,2*l+1
              ! overlap
              la2pm=l**2+m
              overlap(:,la2pm)=fact6*clmamb(:,la2pm)       *clmamb_scalar(la2pm)
              !                ^tau  ^not depend on a & b   ^depend on a & b
              ! gradient of overlap matrix element
              do k_gr=1,3
              ! d/da fact6 =-2.0_r8_kind*xd(k_gr)*fact2*fact6
              overlap_grad(:,la2pm,k_gr)=-2.0_r8_kind*xd(k_gr)*fact2*overlap(:,la2pm)+&
                                            fact6*clmamb(:,la2pm)*diff_arr0(la2pm,k_gr)
              ! d/da clmamb(:,la2pm) = clmamb(:,la2pm)*diff_arr0(la2pm,k_gr)
              enddo

           end do magnetic_number
        end do

     endif standard_ol_grad

     standard_kin_grad: if (integralpar_3cob_grad) then
        magnetic_number_kin: do m=1,2*la+1
           ! kinetic energy
           la2pm=la**2+m
           kinetic=fact5*overlap(:,la2pm) ! a*b/(a+b)(3-2*tau+2*l), tau=fact2*arg

           ! gradient of kinetic energy
          do k_gr=1,3   ! 1st of 2 contribs standard order
           kin_grad(:,m,k_gr)=-4.0_r8_kind*fact2*fact2*xd(k_gr)*overlap(:,la2pm)+&
                fact5*overlap_grad(:,la2pm,k_gr)
!           kin_grad(:,m,k_gr)=0.0_r8_kind !!! for check only
          enddo
        end do magnetic_number_kin

     endif standard_kin_grad

  else standard_order !i.e. reversed order

     reversed_order_olgrad: if (integralpar_2cob_ol_grad) then

        orbital_number: do l=0,la
           magnetic_number_2: do m=1,2*l+1
              ! overlap
              la2pm=l**2+m
              overlap(:,la2pm)=fact6*clmamb(:,la2pm)*clmamb_scalar(la2pm)
              ! gradient of overlap matrix element
              do k_gr=1,3
!              ! d/da fact6 =-2.0_r8_kind*xd(k_gr)*fact2*fact6
              overlap_grad(:,la2pm,k_gr)=-2.0_r8_kind*xd(k_gr)*fact2*overlap(:,la2pm)+&
                                            fact6*clmamb(:,la2pm)*diff_arr0(la2pm,k_gr)
!              ! d/da clmamb(:,la2pm) = clmamb(:,la2pm)*diff_arr0(la2pm,k_gr)
              enddo

           end do magnetic_number_2
        enddo orbital_number !(l)
     endif reversed_order_olgrad

     reversed_order_kin: if (integralpar_3cob_grad) then
        magnetic_number_kin_2: do m=1,2*la+1
           ! kinetic energy
           la2pm=la**2+m
           kinetic=fact5*overlap(:,m)
           ! re-map them to the int_data_2cob3c_stuff
           ! prim_int_2cob_kin(:,:,m,1)=unpack(kinetic,cutoff,zero)
           ! gradient of kinetic energy      2nd of 2 contribs
          do k_gr=1,3   ! 1st of 2 contribs standard order
           kin_grad(:,m,k_gr)=-4.0_r8_kind*fact2*fact2*xd(k_gr)*overlap(:,la2pm)+&
                fact5*overlap_grad(:,la2pm,k_gr)
          enddo
!!!           kin_grad(:,m,:)=0.0_r8_kind   !!!!
        end do magnetic_number_kin_2
     endif reversed_order_kin

  endif standard_order


  do l=1,(la+1)**2
     clmamb(:,l)=clmamb_scalar(l)* clmamb(:,l)
  end do
  fact2=2.0_r8_kind*sqrt(fact0/pi)

  allocate(&
       prod_arr(num,(la+1)**2,0:la),&
       prod_arr_gr(num,6,(la+1)**2,0:la+1),&
       stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("LS_CALCULATE_GRADS : allocation 5 failed")

  call integral_interrupt_2cob3c()


  ! calculate integrals that involve third center,
  ! i.e. gradient of coulomb fit integrals
  ! PROCESS <xi_i|V_nuc|xi_j>, <xi_i|V_pvsp|xi_j>, and <xi_i|V_H|xi_j>
  third_center_required: if (integralpar_3cob_grad) then

     unique_atom_loop: do i=1,n_unique_atoms +n_timps ! loop over third center
        ! allocate space for totalsymmetric gradient with respect
        ! to third center
        if(i<=n_unique_atoms) then
           ua_pointer => unique_atoms(i)
           do_rotation=.true.
           lmax_ch= unique_atoms(i)%lmax_ch      ! maximum l  for chargefit 
           ! determine the maximal angular momentum
           lmax_abs=lmax_ch
           ly_max=max(la,lmax_ch)
           max_order=max(2+la+lmax_abs,4+la)
           z= unique_atoms(i)%z                  ! charge 
           zc=unique_atoms(i)%zc
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'ls_grads: ua=',i,', zero its charge!'
        zc = zero
        z  = zero
           l_max = max((unique_atoms(i)%lmax_pseudo-1),la) ! ?may be deleted
           n_equals=unique_atoms(i)%n_equal_atoms
           imc = unique_atoms(i)%moving_atom
           moving_c = imc > 0
           if (moving_c)then
              grad_dim = gradient_index(imc+1) - gradient_index(imc)
           else
              grad_dim = 0
           endif
        else  ! I am timp
           ua_pointer => unique_timps(i-n_unique_atoms)
           n_equals=ua_pointer%n_equal_atoms
           z  = ua_pointer%z
           zc = ua_pointer%zc
           
           imc = unique_timps(i-n_unique_atoms)%moving_atom
           moving_c = imc > 0
           if (moving_c)then
              grad_dim = gradient_index(N_moving_unique_atoms+imc+1) - &
                   gradient_index(N_moving_unique_atoms+imc)

              allocate(grad_nuc_pseudo_c(num,2*la+1,grad_dim),&                  
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation gard_nuc_pseudo_c for timps failed ")
              grad_nuc_pseudo_c=0.0_r8_kind
           else
              grad_dim = 0
           endif
        end if
        l_max = max((ua_pointer%lmax_pseudo-1),la)
        iam_ppot=zc/=0.0_r8_kind.and.(.not.operations_core_density).and. &
                pseudopot_present
        atom: if(i<=n_unique_atoms) then
           ! --- further allocation ----------------------------------
           ! num : number of pairs(a,b) which are inside the cutoff
           ! for s-and r2-type there is only 1 indep. fct

           allocate ( &
                yl_arr(num,(ly_max+1)**2), &
                gamma_help(num,max_order), &
                gamma_arg2(num), &
                fact10(num,(ly_max+1)**2), &
                diff_arr_grad(num,(ly_max+1)**2,3),&
                nuc_grad(num,2*la+1,grad_dim), &
                stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CALCULATE_GRADS : allocation 7 failed")
           nuc_grad=0.0_r8_kind
           if (iam_ppot) then
              allocate(grad_nuc_pseudo_c(num,2*la+1,grad_dim),&                  
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation gard_nuc_pseudo_c failed ")
              grad_nuc_pseudo_c=0.0_r8_kind        
           endif ! iam_ppot

           allocate ( &
                coul_int_grad_totsym(grad_dim), &
                stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CALCULATE_GRADS : allocation coul_int_grad_totsym failed")

         if(.not.new_3c_co_grad) then
           do k_gr=1,6
              allocate (&
                   coul_int_grad(k_gr)%l(-1:lmax_ch),&
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation coul_int_grad(k_gr)%l failed")
              ncexps = unique_atoms(i)%r2_ch%n_exponents
              allocate(&
                   pointer_r2(num,ncexps,1,nm_lb,nm_la),&
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation coul_int_grad(k_gr)%l(-1)%m failed")
              coul_int_grad(k_gr)%l(-1)%m=>pointer_r2
              ncexps = unique_atoms(i)%l_ch(0)%n_exponents
              allocate(&
                   pointer_s(num,ncexps,1,nm_lb,nm_la),&
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation coul_int_grad(k_gr)%l(0)%m failed")
              coul_int_grad(k_gr)%l(0)%m=>pointer_s
              pointer_s=0.0_r8_kind
              pointer_r2=0.0_r8_kind
           end do

           do i_grad=1,grad_dim ! only if moving_c
              allocate(coul_int_grad_totsym(i_grad)%l(-1:lmax_ch),&
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation failed coul_int_grad_totsym(i_grads)%l")
              ncexps = unique_atoms(i)%r2_ch%n_exponents
              allocate(pointer_r2&
                   (num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
              coul_int_grad_totsym(i_grad)%l(-1)%m=>pointer_r2
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation failed coul_int_grad_totsym(i_grad)%l(-1)%m")
              ncexps = unique_atoms(i)%l_ch(0)%n_exponents
              allocate(pointer_s&
                   (num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
              coul_int_grad_totsym(i_grad)%l(0)%m=>pointer_s
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation failed coul_int_grad_totsym(i_grad)%l(0)%m")
              pointer_s=0.0_r8_kind
              pointer_r2=0.0_r8_kind
           end do

           do i_l=1,lmax_ch
              n_independent_fcts  = &
                   unique_atoms(i)%symadapt_partner(1,i_l)%n_independent_fcts
              ncexps=  unique_atoms(i)%l_ch(i_l)%n_exponents
              do k_gr=1,6
                 allocate(&
                      pointer_l(num,ncexps,n_independent_fcts,nm_lb,nm_la),&
                      stat=alloc_stat)
                 if(alloc_stat/=0) call error_handler &
                      ("LS_CALCULATE_GRADS: allocation coul_int_grad(k_gr)%l(i_l)%m failed")
                 coul_int_grad(k_gr)%l(i_l)%m=>pointer_l
                 pointer_l=0.0_r8_kind
              end do
              do i_grad=1,grad_dim ! only if moving_c
                 allocate(&
                      pointer_l(num,ncexps,n_independent_fcts,nm_lb,nm_la),&
                      stat=alloc_stat)
                 if(alloc_stat/=0) call error_handler &
                      ("LS_CALCULATE_GRADS: allocation coul_int_grad_totsym(i_grad)%l(i_l)%m failed")
                 coul_int_grad_totsym(i_grad)%l(i_l)%m=>pointer_l
                 pointer_l=0.0_r8_kind
              enddo
           end do
          endif


           counter=1
           fact4=1.0_r8_kind
           do i_l=0,ly_max
              do ma=1,2*i_l+1
                 fact10(:,counter)=fact4
                 counter=counter+1
              enddo
              fact4=fact4*aexp_arr/(fact0)
           enddo


           equal_atoms_c: do j=1,n_equals 

              diff_arr_grad=0.0_r8_kind
              if (moving_c) then
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
              xc=ua_pointer%position(:,j)
              yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                   spread(xc,1,num))
              gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2&
                   +(gamma_arg(:,3)-xc(3))**2)*fact0
              fact4=1.0_r8_kind
              do l=1,(ly_max+1)**2
                 do k=1,solhrules_differential(3,l)%n_summands
                    diff_arr_grad(:,l,1)=diff_arr_grad(:,l,1)+&
                         solhrules_differential(3,l)%coef(k)*&
                         yl_arr(:,solhrules_differential(3,l)%lm_sh(k))   
                 end do
                 do k=1,solhrules_differential(4,l)%n_summands
                    diff_arr_grad(:,l,2)=diff_arr_grad(:,l,2)+&
                         solhrules_differential(4,l)%coef(k)*&
                         yl_arr(:,solhrules_differential(4,l)%lm_sh(k))
                 end do
                 do k=1,solhrules_differential(2,l)%n_summands
                    diff_arr_grad(:,l,3)=diff_arr_grad(:,l,3)+&
                         solhrules_differential(2,l)%coef(k)*&
                         yl_arr(:,solhrules_differential(2,l)%lm_sh(k))
                 end do
              end do!loop over l

              ! derivative of prod_arr with respect to nucleus a and nucleus c
              call calculate_helpers()

              ! PROCESS <xi_i|Zc/|r-Rc||xi_j>
              if(integralpar_relativistic.and.iam_ppot) then
                call calc_final_nuc_or_pvsp_grad(pseudo_grad_gr,nuc_grad)
              else                 
                call calc_final_nuc_or_pvsp_grad(nuc_grad_gr,nuc_grad)
!              print*,'calc_final_nuc_or_pvsp_grad off'
              end if
           
              ! PROCESS [xi_i|f_k|xi_j]
              ! first evaluate s-type coulomb integrals
              if(old_3c_co) call s_coulomb(i)

              ! now treating r2-type coloumb integrals
              if(old_3c_co) call r2_coulomb(i)

              ! finally calculate l-type coulomb integrals
              if(old_3c_co) call l_coulomb(i,j)

           end do equal_atoms_c
        end if atom
        ! PROCESS  pseudopotentials
        pp: if (iam_ppot) then
          ABORT('not supported')
        endif pp
       
        grad_3cob_required: if(integralpar_3cob_grad.and.i<=n_unique_atoms) then

           ! coulombic contributions of derivatives with respect
           ! to center a and b
           ! CONTRACT d/dRa [xi_i|f_k|xi_j] and d/dRb [xi_i|f_k|xi_j]
           fit_stand_order: if(lagtlb) then
              do k_gr=k_gr_0,k_gr_1
                 grad_mat_p=>grad_mat_gr(:,:,:,:,k_gr)
                 if (model_density .and. spin_polarized) then
                    grad_mat_sp => grad_mat_spgr(:,:,:,:,k_gr)
                    if(old_3c_fc.and.old_3c_co) &
                    call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr),&
                         grad_mat_p,grad_mat_sp)

                 else ! standard SCF
                    if(old_3c_fc.and.old_3c_co) &
                    call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr),&
                      grad_mat_p)
                 endif
              end do
           else fit_stand_order
              do k_gr=1,3
                 if (moving_b) then ! load d/dRb from coul_int_grad(1:3)
                 grad_mat_p=>grad_mat_gr(:,:,:,:,k_gr+3)
                    if (model_density .and. spin_polarized) then
                       grad_mat_sp => grad_mat_spgr(:,:,:,:,k_gr+3)
                       if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff,coul_int_grad( &
                            k_gr),grad_mat_p,grad_mat_sp)
                    else ! standard SCF
                     if(old_3c_fc.and.old_3c_co) &
                     call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr),&
                      grad_mat_p)
                   endif
                 endif
                 if (moving_a) then ! load d/dRa from coul_int_grad(4:6)
                 grad_mat_p=>grad_mat_gr(:,:,:,:,k_gr)
                    if (model_density .and. spin_polarized) then
                       grad_mat_sp => grad_mat_spgr(:,:,:,:,k_gr)
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff,coul_int_grad( &
                            k_gr+3),grad_mat_p,grad_mat_sp)
                    else ! standard SCF
                     if(old_3c_fc.and.old_3c_co) &
                     call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr+3),&
                      grad_mat_p)
                    endif
                 endif
              end do
           endif fit_stand_order

         if(.not.new_3c_co_grad) then
           do k_gr=1,6
              do i_l = -1, lmax_ch
                 deallocate(coul_int_grad(k_gr)%l(i_l)%m,&
                   STAT=alloc_stat)
                 if(alloc_stat/=0) call error_handler &
                      ("LS_CALCULATE_GRADS : deallocation coul_int_grad_%l%m failed")
              enddo
              deallocate (coul_int_grad(k_gr)%l,&
                   STAT=alloc_stat)
              if(alloc_stat/=0) call error_handler &
                ("LS_CALCULATE_GRADS : deallocation coul_int_grad_%l failed")
           end do
          endif

           ! CONTRACT dsym/dRc [xi_i|f_k|xi_j]
           do i_grad=1,grad_dim ! only if moving_c
              index = gradient_index(imc) + i_grad - 1
              if (model_density) then
                 if (spin_polarized) then
                    spin_index = index + gradient_data_n_gradients
                    if (split_gradients) then
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff, &
                            coul_int_grad_totsym(i_grad), &
                            prim_int_3cob_coul_grad(index)%m, & ! V_H
                            prim_int_3cob_grad(spin_index)%m, & ! V_X,spin
                            prim_int_3cob_grad(index)%m )       ! V_X,tot
                    else ! total_gradient only
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff, &
                            coul_int_grad_totsym(i_grad), &
                            prim_int_3cob_grad(index)%m, &      ! V_H+V_X,tot
                            prim_int_3cob_grad(spin_index)%m )  ! V_X,spin
                    endif
                 else ! spin_restricted
                    if (split_gradients) then
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff, &
                            coul_int_grad_totsym(i_grad), &
                            prim_int_3cob_coul_grad(index)%m, & ! V_H
                            mda_xcpot_gradients = &                       
                            prim_int_3cob_grad(index)%m )       ! V_X,tot
                    else ! total_gradient only
                     if(old_3c_fc.and.old_3c_co) &
                       call fitcontract('grad',num,i,cutoff, &
                            coul_int_grad_totsym(i_grad), &
                            prim_int_3cob_grad(index)%m )       ! V_H+V_X,tot
                    endif
                 endif
              else ! standard SCF
                     if(old_3c_fc.and.old_3c_co) &
                      call fitcontract('grad',num,i,cutoff, &
                      coul_int_grad_totsym(i_grad), &
                      prim_int_3cob_grad(index)%m )
              endif

           if(.not.new_3c_co_grad) then
              do i_l = -1, lmax_ch
                 deallocate(coul_int_grad_totsym(i_grad)%l(i_l)%m,&
                      STAT=alloc_stat)
                 if(alloc_stat/=0) call error_handler &
                      ("LS_CALCULATE_GRADS : deallocation coul_int_grad_totsym%l%m failed")
              enddo
              deallocate (coul_int_grad_totsym(i_grad)%l,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("LS_CALCULATE_GRADS : deallocation coul_int_grad_totsym(i_grad)%l failed")
            endif
           end do
           deallocate(coul_int_grad_totsym,stat=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("LS_CALCULATE_GRADS : deallocation coul_int_grad_totsym failed")
 
           ! now add contribution of derivative of nuc with respect to third center
           ! ADD dsym/dRc <xi_i|V_nuc[c]|xi_j>, dsym/dRc <xi_i|V_pvsp[c]|xi_j>
           if (moving_c) index = gradient_index(imc) - 1
           do i_grad=1,grad_dim ! only if moving_c
              notrel: if(.not.integralpar_relativistic) then
                 if( split_gradients )then
                    counter=1
                    do i_la=1,nm_la
                       do i_lb=1,nm_lb
                          prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,i_lb,i_la) = &
                               prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,i_lb,i_la) - &
                               unpack(nuc_grad(:,counter,i_grad),cutoff,zero)
                          counter=counter+1
                       end do
                    end do
                    !:AH[
                    if (iam_ppot) then
                       counter=1
                       do i_la=1,nm_la
                          do i_lb=1,nm_lb
                             prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,i_lb,i_la) = &
                                  prim_int_3cob_nuc_grad(index+i_grad)%m(:,:,i_lb,i_la) + &
                                  unpack(grad_nuc_pseudo_c(:,counter,i_grad),cutoff,zero)
                             counter=counter+1
                          end do
                       end do
                    endif
                    !:AH]

                 else ! i.e. not split gradients
                    counter=1
                    do i_la=1,nm_la
                       do i_lb=1,nm_lb      ! nuc_grad 1st of 2 contribs
                          prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                               prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la) &
                                -unpack(nuc_grad(:,counter,i_grad),cutoff,zero)
                          counter=counter+1
                       end do
                    end do
                    !:AH[
                    if (iam_ppot) then
                       counter=1
                       do i_la=1,nm_la
                          do i_lb=1,nm_lb
                             prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                                  prim_int_3cob_grad(index+i_grad)%m&
                                  (:,:,i_lb,i_la)+unpack(grad_nuc_pseudo_c &
                                  (:,counter,i_grad),cutoff,zero)
                             counter=counter+1
                          end do
                       end do
                    endif
!:AH] 
              endif

           else notrel ! i.e. integralpar_relativistic!!!!
              counter=1
                 do i_la=1,nm_la
                    do i_lb=1,nm_lb
                    if (iam_ppot) then
                            prim_int_2cob_pseudo_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                            prim_int_2cob_pseudo_grad(index+i_grad)%m&
                            (:,:,i_lb,i_la) -unpack(nuc_grad &
                                      (:,counter,i_grad),cutoff,zero)

                    else
                       prim_int_2cob_nuc_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                            prim_int_2cob_nuc_grad(index+i_grad)%m&
                            (:,:,i_lb,i_la)-unpack(nuc_grad(:,counter,i_grad),cutoff,zero)
                    end if
                    
                       counter=counter+1
                    end do
                 end do
                if (iam_ppot) then ! this is true pp contrib
                    counter=1
                    do i_la=1,nm_la
                       do i_lb=1,nm_lb
                         prim_int_2cob_pseudo_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                            prim_int_2cob_pseudo_grad(index+i_grad)%m&
                            (:,:,i_lb,i_la) +unpack(grad_nuc_pseudo_c &
                         	(:,counter,i_grad),cutoff,zero)


                         counter=counter+1
                       end do
                    end do
                
                 endif
              endif notrel
           end do

        elseif(integralpar_3cob_grad.and.moving_c) then
           index = gradient_index(imc+n_moving_unique_atoms) - 1
           do i_grad=1,grad_dim ! only if moving_c
              counter=1
              do i_la=1,nm_la
                 do i_lb=1,nm_lb
                    prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la) = &
                         prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la) + &
                         unpack(grad_nuc_pseudo_c(:,counter,i_grad),cutoff,zero)
                    counter=counter+1
                 end do
              end do
           end do

        endif grad_3cob_required
        
        if(i<=n_unique_atoms) then
           deallocate( yl_arr, gamma_help, gamma_arg2, fact10,&
                diff_arr_grad, nuc_grad, stat=alloc_stat )
           if(alloc_stat/=0) call error_handler &
                ("LS_CALCULATE_GRADS: deallocation 7 failed ") 
           !:AH[
           if (iam_ppot) then
              deallocate(grad_nuc_pseudo_c,stat=alloc_stat )
              if(alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS: deallocation of grad_pseudo_c failed ")
           endif
           !:AH]          
        else 
           if(moving_c) then

              deallocate(grad_nuc_pseudo_c,stat=alloc_stat ) ! for timp
              if(alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS: deallocation of grad_pseudo_c for timps failed ")
           end if

        end if
     end do unique_atom_loop
     
     unique_charge_loop: do i=1,pointcharge_N+n_timps  ! loop over third center
        iam_ppot=i<=n_timps
        if(iam_ppot) then
           z= unique_timps(i)%z - unique_timps(i)%zc
           n_equals= unique_timps(i)%n_equal_atoms
           imc = unique_timps(i)%moving_atom
           moving_c = imc > 0
           if (moving_c)then
              grad_dim = gradient_index(N_moving_unique_atoms+imc+1) - &
                   gradient_index(N_moving_unique_atoms+imc)
             allocate ( nuc_grad(num,2*la+1,grad_dim), &
                   stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE_GRADS : allocation nuc_grad fot timps failed")
              nuc_grad=0.0_r8_kind


           else
              grad_dim = 0
              moving_c=.false.
           endif
        else
           cycle unique_charge_loop ! pointcharges go to SHGI !!!!!!!!!!AS
           moving_c=.false.
           z= pointcharge_array(i-n_timps)%z   ! charge 
           n_equals=pointcharge_array(i-n_timps)%n_equal_charges
        endif
        ly_max=la
        max_order=4+la

        allocate ( &
             yl_arr(num,(ly_max+1)**2), &
             gamma_help(num,max_order), &
             gamma_arg2(num), &
             diff_arr_grad(num,(ly_max+1)**2,3),&
             stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE_GRADS : allocation pc failed")
        equal_charges_c: do j=1,n_equals
           iam_ppot= i<=n_timps
           if(iam_ppot) then
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
              endif
              xc = unique_timps(i)%position(:,j)
           else
              xc = pointcharge_array(i-n_timps)%position(:,j)
           end if
           yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                spread(xc,1,num))
           gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2&
                +(gamma_arg(:,3)-xc(3))**2)*fact0
           fact4=1.0_r8_kind
           diff_arr_grad=0.0_r8_kind
           do l=1,(ly_max+1)**2
              do k=1,solhrules_differential(3,l)%n_summands
                 diff_arr_grad(:,l,1)=diff_arr_grad(:,l,1)+&
                      solhrules_differential(3,l)%coef(k)*&
                      yl_arr(:,solhrules_differential(3,l)%lm_sh(k))   
              end do
              do k=1,solhrules_differential(4,l)%n_summands
                 diff_arr_grad(:,l,2)=diff_arr_grad(:,l,2)+&
                      solhrules_differential(4,l)%coef(k)*&
                      yl_arr(:,solhrules_differential(4,l)%lm_sh(k))
              end do
              do k=1,solhrules_differential(2,l)%n_summands
                 diff_arr_grad(:,l,3)=diff_arr_grad(:,l,3)+&
                      solhrules_differential(2,l)%coef(k)*&
                      yl_arr(:,solhrules_differential(2,l)%lm_sh(k))
              end do
           end do!loop over l

           ! derivative of prod_arr with respect to nucleus a and nucleus c
           call calculate_helpers()

           ! ! PROCESS <xi_i|Zp/|r-Rp||xi_j>
           if(integralpar_relativistic.and.pseudopot_present) then
              call calc_final_nuc_or_pvsp_grad(pseudo_grad_gr)
           else
              if(moving_c) then
                 call calc_final_nuc_or_pvsp_grad(nuc_grad_gr,nuc_grad) !extended for timps
              else
              call calc_final_nuc_or_pvsp_grad(nuc_grad_gr,nuc_grad) !extended for timps
	      endif
           end if
           
        enddo equal_charges_c

        if(moving_c) then
           index = gradient_index(imc+n_moving_unique_atoms) - 1
           do i_grad=1,grad_dim ! only if moving_c
              counter=1
              do i_la=1,nm_la
                 do i_lb=1,nm_lb
                    prim_int_3cob_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                         prim_int_3cob_grad(index+i_grad)%m&
                         (:,:,i_lb,i_la)-unpack(nuc_grad(:,counter,i_grad),cutoff,zero)
                    counter=counter+1
                 end do
              end do
           end do
        endif
     
        deallocate ( &
             yl_arr, &
             gamma_help, &
             gamma_arg2, &
             diff_arr_grad,&
             stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE_GRADS : deallocation pc failed")
        if(moving_c) then
           deallocate (nuc_grad,stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CALCULATE_GRADS : deallocationi nuc_grad for timps failed")
        endif
        
     enddo unique_charge_loop

  end if third_center_required

  par_solv_grad: if(integralpar_solv_grad.and.old_solv_grad) then !!!!!!!!!!!!!!!!!!
     !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

        ly_max=la
        max_order=4+la
        
        allocate ( &
             yl_arr(num,(ly_max+1)**2), &
             gamma_help(num,max_order), &
             gamma_arg2(num), &
             fact10(num,(ly_max+1)**2), &
             diff_arr_grad(num,(ly_max+1)**2,3),&
             stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE_GRADS : allocation 7 failed")

        counter=1
        fact4=1.0_r8_kind
        do i_l=0,ly_max
           do ma=1,2*i_l+1
              fact10(:,counter)=fact4
              counter=counter+1
           enddo
           fact4=fact4*aexp_arr/(fact0)
        enddo

           solv_grad_gr=0.0_r8_kind 
           help_arr_gr1=0.0_r8_kind

           N_pc=0
           if(with_pc .and. .not.fixed_pc)  N_pc=pointcharge_N

           lab_k: do k1=1,to_calc_grads%n_points 
              call integral_interrupt_2cob3c()
              n_equal_solv=to_calc_grads%n_equal(k1)
              z=to_calc_grads%Q(k1) 
              
              lab_m: do m=1,n_equal_solv
                 xc =to_calc_grads%xyz(k1,m,:) 
                 yl_arr(:,:)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
                      spread(xc,1,num))
                 gamma_arg2(:)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2&
                      +(gamma_arg(:,3)-xc(3))**2)*fact0
                 fact4=1.0_r8_kind
                 diff_arr_grad=0.0_r8_kind
                 do l=1,(ly_max+1)**2
                    do k=1,solhrules_differential(3,l)%n_summands
                       diff_arr_grad(:,l,1)=diff_arr_grad(:,l,1)+&
                            solhrules_differential(3,l)%coef(k)*&
                            yl_arr(:,solhrules_differential(3,l)%lm_sh(k))   
                    end do
                    do k=1,solhrules_differential(4,l)%n_summands
                       diff_arr_grad(:,l,2)=diff_arr_grad(:,l,2)+&
                            solhrules_differential(4,l)%coef(k)*&
                            yl_arr(:,solhrules_differential(4,l)%lm_sh(k))
                    end do
                    do k=1,solhrules_differential(2,l)%n_summands
                       diff_arr_grad(:,l,3)=diff_arr_grad(:,l,3)+&
                            solhrules_differential(2,l)%coef(k)*&
                            yl_arr(:,solhrules_differential(2,l)%lm_sh(k))
                    end do
                 end do!loop over l

                 call calculate_helpers()
                 call calc_final_nuc_or_pvsp_grad(solv_grad_gr) ! only solv_grad_gr
                 ism=to_calc_grads%i_symm_sort(k1,m)

                 unique_moving_loop1: do i=1,N_moving_unique_atoms+N_pc
                    if(i <= N_moving_unique_atoms) then
                       grad_dim = gradient_index(i+1) - gradient_index(i)
                       l=moving_unique_atom_index(i)
                       n_equals=unique_atoms(l)%n_equal_atoms
                    else
                       l=i-N_moving_unique_atoms
                       n_equals=pointcharge_array(l)%N_equal_charges
                       grad_dim = surf_points_grad_index(l+1)-surf_points_grad_index(l)
                    end if

                   help_arr_gr1=0.0_r8_kind

!                  equal_atoms_cc: do j=1,n_equals 

                   do n=1,grad_dim
                    do n1=4,6
                       help_arr_gr1(:,n,:,:)=help_arr_gr1(:,n,:,:)+ &
!                           help_arr_gr(:,n1,:,:)*to_calc_grads%dxyz_totsym(n,i,j)%m(n1-3,ism)
                            help_arr_gr(:,n1,:,:)*to_calc_grads%dxyz_totsyms(n,i)%m(n1-3,ism)
                    enddo
                   enddo

!                  enddo equal_atoms_cc

                   if(i <= N_moving_unique_atoms) then
                      index = gradient_index(i) - 1
                      do i_grad=1,grad_dim 
                         do i_la=1,nm_la
                            do i_lb=1,nm_lb
                               prim_int_3cob_solv_grad(index+i_grad)%m(:,:,i_lb,i_la)=&
                                    prim_int_3cob_solv_grad(index+i_grad)%m(:,:,i_lb,i_la)+&
                                    unpack(help_arr_gr1(:,i_grad,i_lb,i_la),cutoff,zero)
                            end do
                         end do
                      enddo
                   else
                      index = surf_points_grad_index(l) - 1
                      do i_grad=1,grad_dim 
                         do i_la=1,nm_la
                            do i_lb=1,nm_lb
                               prim_int_3cob_solv_grad_pc(index+i_grad)%m(:,:,i_lb,i_la)=&
                                    prim_int_3cob_solv_grad_pc(index+i_grad)%m(:,:,i_lb,i_la)-&
                                    unpack(help_arr_gr1(:,i_grad,i_lb,i_la),cutoff,zero)
                            end do
                         end do
                      enddo
                   end if

                 enddo unique_moving_loop1
             
              enddo lab_m
           enddo lab_k

!!$if(with_pc .not. fixed_pc) then
!!$do i=1,size(prim_int_3cob_solv_grad_pc)
!!$print*,na,nb,'c',i,prim_int_3cob_solv_grad_pc(i)%m(1,1,1,1)
!!$end do
!!$end if
        deallocate( yl_arr, gamma_help, gamma_arg2, fact10,&
             diff_arr_grad, stat=alloc_stat )
        if(alloc_stat/=0) call error_handler &
             ("LS_CALCULATE_GRADS: deallocation 8 failed ")

     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  endif par_solv_grad     !!!!!!!!!!!!!!!!!!!!!

  ! remap nuc to the appropriate data structure
  grad_remap_required: if (integralpar_3cob_grad) then

     if(.not.integralpar_relativistic) then
        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>

       if(.true.) then   !!!!
        stand_nuc: if(lagtlb) then
           do k_gr=k_gr_0,k_gr_1
              do m=1,2*la+1 
                 grad_mat_gr(:,:,1,m,k_gr)=grad_mat_gr(:,:,1,m,k_gr)-&
                      unpack(nuc_grad_gr(:,m,k_gr),cutoff,zero)
              enddo
           end do
        else stand_nuc
           do k_gr=1,3
              do m=1,2*la+1 
                 if (moving_a) then ! load d/dRa from nuc_grad_gr(...,4:6)
                 grad_mat_gr(:,:,m,1,k_gr)=grad_mat_gr(:,:,m,1,k_gr)-&
                      unpack(nuc_grad_gr(:,m,k_gr+3),cutoff,zero)
                 endif
                 if (moving_b) then ! load d/dRb from nuc_grad_gr(...,1:3)
                 grad_mat_gr(:,:,m,1,k_gr+3)=grad_mat_gr(:,:,m,1,k_gr+3)-&
                      unpack(nuc_grad_gr(:,m,k_gr),cutoff,zero)
                 endif
              enddo
           end do
        endif stand_nuc
      endif

        if (pseudopot_present) then
           if(lagtlb) then
              do k_gr=k_gr_0,k_gr_1
                 do m=1,2*la+1
                    grad_mat_gr(:,:,1,m,k_gr)=grad_mat_gr(:,:,1,m,k_gr)+&
                    unpack(grad_nuc_pseudo_ab(:,m,k_gr),cutoff,zero)
                 enddo
              end do
           else
              do k_gr=1,3
                 do m=1,2*la+1
                    if (moving_a) then ! load d/dRa from nuc_grad_gr(...,4:6)
                       grad_mat_gr(:,:,m,1,k_gr)=grad_mat_gr(:,:,m,1,k_gr)+&
                       unpack(grad_nuc_pseudo_ab(:,m,k_gr),cutoff,zero)
                    endif
                    if (moving_b) then ! load d/dRb from nuc_grad_gr(...,1:3)
                    grad_mat_gr(:,:,m,1,k_gr+3)=grad_mat_gr(:,:,m,1,k_gr+3)+&
                        unpack(grad_nuc_pseudo_ab(:,m,k_gr+3),cutoff,zero)
                    endif
                 enddo
              end do
           end if
        end if! pseudopot_present
     endif! not integralpar_relativistic
     ! LOAD prim_int_3cob_grad(...,ia:) and prim_int_3cob_grad(...,ib: ) or
     !      prim_int_2cob_ks_grad( 1:N) and prim_int_2cob_ks_grad(N+1:M)

     if (split_gradients) then
        pointer_prim_int=>prim_int_2cob_ks_grad
        call add_grads(grad_mat_gr,compact_storage=.TRUE.)
        if (model_density .and. spin_polarized) then
           spin_index = size(prim_int_2cob_ks_grad)/2 + 1
           pointer_prim_int => prim_int_2cob_ks_grad(spin_index:)
           call add_grads(grad_mat_spgr,compact_storage=.TRUE.)
        endif

     else ! i.e. .not.split_gradients

        pointer_prim_int=>prim_int_3cob_grad
        call add_grads(grad_mat_gr)

        if (model_density .and. spin_polarized) then
           spin_index = gradient_data_n_gradients + 1
           pointer_prim_int => prim_int_3cob_grad(spin_index:)
           call add_grads(grad_mat_spgr)
        endif
     endif

     relmap: if(integralpar_relativistic) then

        ! ADD d/dRa <xi_i|V_nuc|xi_j> and d/dRb <xi_i|V_nuc|xi_j>
        grad_mat_gr=0.0_r8_kind
        if(lagtlb) then
           do k_gr=k_gr_0,k_gr_1
              do m=1,2*la+1 
                 grad_mat_gr(:,:,1,m,k_gr)=grad_mat_gr(:,:,1,m,k_gr)-&
                      unpack(nuc_grad_gr(:,m,k_gr)&
                             ,cutoff,zero)
               
              enddo
           end do
        else
           do k_gr=1,3
              do m=1,2*la+1 
                 if (moving_a) then ! load d/dRa from nuc_grad_gr(...,4:6)
                 grad_mat_gr(:,:,m,1,k_gr)=grad_mat_gr(:,:,m,1,k_gr)-&
                      unpack(nuc_grad_gr(:,m,k_gr+3)&
                               ,cutoff,zero)                                
                 endif
                 if (moving_b) then ! load d/dRb from nuc_grad_gr(...,1:3)
                 grad_mat_gr(:,:,m,1,k_gr+3)=grad_mat_gr(:,:,m,1,k_gr+3)-&
                      unpack(nuc_grad_gr(:,m,k_gr)&
                               ,cutoff,zero)
              endif
              enddo
           end do
        end if

        ! LOAD prim_int_2cob_nuc_grad(...,ia:), prim_int_2cob_nuc_grad(...,ib:)
        pointer_prim_int=>prim_int_2cob_nuc_grad
        call add_grads(grad_mat_gr)

        if (pseudopot_present) then
          grad_mat_gr=0.0_r8_kind
           if(lagtlb) then
              do k_gr=k_gr_0,k_gr_1
                 do m=1,2*la+1
                    grad_mat_gr(:,:,1,m,k_gr)=grad_mat_gr(:,:,1,m,k_gr)+&
                    unpack(grad_nuc_pseudo_ab(:,m,k_gr) &
                          - pseudo_grad_gr(:,m,k_gr) &
                            ,cutoff,zero)
                 enddo
              end do
           else
              do k_gr=1,3
                 do m=1,2*la+1
                    if (moving_a) then ! load d/dRa from nuc_grad_gr(...,4:6)
                       grad_mat_gr(:,:,m,1,k_gr)=grad_mat_gr(:,:,m,1,k_gr)+&
                       unpack(grad_nuc_pseudo_ab(:,m,k_gr)&
                       - pseudo_grad_gr(:,m,k_gr+3) &
                                       ,cutoff,zero)
                    endif
                    if (moving_b) then ! load d/dRb from nuc_grad_gr(...,1:3)
                    grad_mat_gr(:,:,m,1,k_gr+3)=grad_mat_gr(:,:,m,1,k_gr+3)+&
                        unpack(grad_nuc_pseudo_ab(:,m,k_gr+3)&
                        - pseudo_grad_gr(:,m,k_gr) &
                                        ,cutoff,zero)
                    endif
                 enddo
              end do
           end if
        pointer_prim_int=>prim_int_2cob_pseudo_grad
        call add_grads(grad_mat_gr)
        end if ! pseudopot_present

     end if relmap
  end if grad_remap_required

  grad_remap_solv: if(integralpar_solv_grad.and.old_solv_grad) then !!!!!!!!!!!!!!!!!!!!!!!
     !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
     if(lagtlb) then
        do k_gr=k_gr_0,k_gr_1
           do m=1,2*la+1 
              grad_mat_gr(:,:,1,m,k_gr)=grad_mat_gr(:,:,1,m,k_gr)-&
                   unpack(solv_grad_gr(:,m,k_gr),cutoff,zero)
           enddo
        end do
     else
        do k_gr=1,3
           do m=1,2*la +1 
              if (moving_a) then ! load d/dRa from nuc_grad_gr(...,4:6)
                 grad_mat_gr(:,:,m,1,k_gr)=grad_mat_gr(:,:,m,1,k_gr)-&
                      unpack(solv_grad_gr(:,m,k_gr+3),cutoff,zero)
              endif
              if (moving_b) then ! load d/dRb from nuc_grad_gr(...,1:3)
                 grad_mat_gr(:,:,m,1,k_gr+3)=grad_mat_gr(:,:,m,1,k_gr+3)-&
                      unpack(solv_grad_gr(:,m,k_gr),cutoff,zero)
              endif
           enddo
        end do
     end if

     pointer_prim_int=>prim_int_3cob_solv_grad
     call add_grads(grad_mat_gr)
     if (model_density .and. spin_polarized) then
        spin_index = gradient_data_n_gradients + 1
        pointer_prim_int => prim_int_3cob_solv_grad(spin_index:)
        call add_grads(grad_mat_spgr)
     endif
     !!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  endif grad_remap_solv  !!!!!!!!!!!!!!!!!

  if (integralpar_3cob_grad) then

     deallocate( kinetic, nuc_grad_gr, grad_mat_gr, &
          help_arr_gr, help_arr, help_arr2, stat=alloc_stat)
     ASSERT(alloc_stat.eq.0)
     
     if (pseudopot_present) then
        deallocate( grad_nuc_pseudo_ab, &
                  stat=alloc_stat)
        if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE_GRADS: deallocation of pseudopotential failed")
        if(integralpar_relativistic) then
          deallocate(pseudo_grad_gr,stat=alloc_stat)
         if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE_GRADS: deallocation of pseudo_grad_gr failed")
        end if
        
     endif
     if (model_density .and. spin_polarized) then
        deallocate( grad_mat_spgr, stat=alloc_stat)
        if (alloc_stat /= 0) call error_handler &
             ("LS_CALCULATE_GRADS: deallocation 3m failed")
     end if
  end if

  if(integralpar_solv_grad.and.old_solv_grad) then                               !!!!!!!!!!!!
     deallocate( solv_grad_gr, grad_mat_gr,help_arr_gr1, &                  !!!!!!!!!!!!
          help_arr_gr, help_arr, help_arr2, stat=alloc_stat)   !!!!!!!!!!!!
     if (alloc_stat.ne.0) call error_handler &                 !!!!!!!!!!!!
          ("LS_CALCULATE_GRADS: deallocation 3a failed")       !!!!!!!!!!!!
     if (model_density .and. spin_polarized) then              !!!!!!!!!!!!
        deallocate( grad_mat_spgr, stat=alloc_stat)            !!!!!!!!!!!!
        if (alloc_stat /= 0) call error_handler &              !!!!!!!!!!!!
             ("LS_CALCULATE_GRADS: deallocation 3ma failed")   !!!!!!!!!!!!
     end if                                                    !!!!!!!!!!!!
  endif                                                        !!!!!!!!!!!!

  deallocate(fact0,fact1,fact2,fact4,fact5,fact6,fact7,fact8,&
       rcsabc,tau,rsaba,aexp_arr,bexp_arr,overlap,overlap_grad,&
       kin_grad,gamma_arg,clmamb_scalar,clmamb,diff_arr,diff_arr_gr,&
       help_vec,help_mat,&
       cutoff,diff_arr0, &
       prod_arr,prod_arr_gr,&
       stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
end subroutine ls_calculate_grads

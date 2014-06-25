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
  subroutine calc_3center(na,la,nb,lb,equalb,equala,IMODE)
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
! End of public interface of module
!================================================================
#include "def.h"
    use type_module, only: i4_kind, r8_kind, i8_kind
    use calc_3center_module
    use cpksdervs_matrices
    use fitcontract_module
    use fit_coeff_module, only: fit_coeff_n_ch
#ifdef FPP_DEBUG
    use error_module, only: MyID
#endif
    use shgi_cntrl, only: setif, INUSD
    use gradient_data_module, only: cpks_gradient_totalsym
    use unique_atom_module, only: unique_atoms
    use int_data_2cob3c_module, only: quadrupel &
                                    , prim_int_3c_co &
                                    , center1 &
                                    , center2
    use integralpar_module, only: integralpar_gradients &
                                , integralpar_cpks_contribs &
                                , integralpar_solv_grad &
                                , integralpar_3c_co &
                                , integralpar_2cob_potential &
                                , integralpar_2cob_nuc &
                                , integralpar_dervs
    use calc3c_switches,    only: shift_xb &
                                , shift_xa

    ! FIXME: many of these are re-exported from other modules:
    use ll_calculate_grads_module, only: ima &
                                       , imb &
                                       , imc &
                                       , moving_a &
                                       , moving_b &
                                       , moving_c &
                                       , unity_matrix_p &
                                       , k_gr_0 &
                                       , k_gr_1 &
                                       , grad_mat_spin &
                                       , spin_index &
                                       , spin_polarized &
                                       , options_spin_restricted &
                                       , model_density &
                                       , options_xcmode &
                                       , xcmode_model_density &
                                       , split_gradients &
                                       , options_split_gradients &
                                       , add_nuc_or_pvsp &
                                       , add_cpks_coul &
                                       , add_dervs
    USE_MEMLOG
    implicit none

    integer(kind=i4_kind), optional,intent(in) :: equala
    integer(kind=i4_kind), optional,intent(in) :: equalb
    integer(kind=i4_kind), intent(in)     :: na ! number of unique atom a
    integer(kind=i4_kind), intent(in)     :: nb ! number of unique atom b
    integer(kind=i4_kind), intent(in)     :: la ! angular momentum of unique atom a
    integer(kind=i4_kind), intent(in)     :: lb ! angular momentum of unique atom b
    integer(kind=i8_kind), optional,intent(in) :: IMODE ! bitmask
    ! *** end of interface ***

    logical, save                         :: initialized = .false.

    real(kind=r8_kind), parameter         :: pi2 = pi*pi
    real(kind=r8_kind), parameter         :: zero = 0.0_r8_kind &
                                           , two  = 2.0_r8_kind
    character(len=50)                     :: quad

    ! declared in calc_3center_module:
!   type(arrmat5), allocatable :: coul_int(:) ! (-1:lmax_ch)

    integer (i4_kind) :: i_grad
    integer (i4_kind) :: k2dr
    integer (i4_kind) :: lp


      call setif( -1_i8_kind, .false. ) ! zero all bits
   if( present(IMODE) )then
      ! overwrite defaults:
      call setif( IMODE )
   else
      ! default values:
      call setif( INUSD )
   endif

  if(integralpar_gradients) then
   MEMSET(0)
  endif
  ! this logic is hopefully present in integralpar_setup():
! integralpar_dervs=integralpar_gradients &
!                   .and.integralpar_2dervs.and..not.integralpar_cpks_contribs
#ifdef no_cpks_coul_grads
    if(integralpar_cpks_contribs) cpks_gradient_totalsym=0.0_r8_kind ! no charge overlap contribs
#endif

    ! reset allocation flags on each entry:
    alloc_stat(:)=0


    ! tabulate factorials, once:
    if( .not.initialized )then
      call calc_3c_init()
      initialized = .true.
    endif


      if(l_print) then
         print*,integralpar_gradients,integralpar_solv_grad
         print*,integralpar_2cob_nuc,integralpar_3c_co,integralpar_2cob_potential
         write(quad,"(4i3)") quadrupel%ua1, quadrupel%l1, &
         quadrupel%ua2, quadrupel%l2
         call write_to_output_units("calc_3center : quadrupel "//trim(quad))
       endif

!          write(quad,"(4i3)") quadrupel%ua1, quadrupel%l1, &
!          quadrupel%ua2, quadrupel%l2
!         call write_to_output_units("calc_3center : quadrupel "//trim(quad))
!        if(integralpar_gradients) then
!       if(print_quad)  print*, 'quad 3c',quadrupel%ua1,quadrupel%l1, quadrupel%ua2, quadrupel%l2
!        endif

!      print*, 'quad 3c',quadrupel%ua1,quadrupel%l1, quadrupel%ua2, quadrupel%l2, &
!              'integralpar_cpksdervs', integralpar_cpksdervs,'integralpar_dervs',integralpar_dervs
!

    xa = center1
    if(integralpar_gradients.and.shift_xa) then
     xa(1)=xa(1)+0.001_r8_kind
     print*,'XA COORDINATE 0.001 SHIFTED'
    endif
    xb = center2
    if(integralpar_gradients.and.shift_xb) then
     xb(1)=xb(1)+0.001_r8_kind
     print*,'XB COORDINATE 0.001 SHIFTED'
    endif
    xd = xa-xb
    arg=sum(xd**2)
    aexps => unique_atoms(na)%l_ob(la)%exponents(:)
    naexps = unique_atoms(na)%l_ob(la)%n_exponents
    bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)
    nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents
    if(integralpar_gradients) then
     ima = unique_atoms(na)%moving_atom
     imb = unique_atoms(nb)%moving_atom
     moving_a = ima > 0
     moving_b = imb > 0
     k_gr_0 = 4; if (moving_a) k_gr_0 = 1
     k_gr_1 = 3; if (moving_b) k_gr_1 = 6
     split_gradients = options_split_gradients()
     model_density = options_xcmode() == xcmode_model_density
     spin_polarized = .not. options_spin_restricted()
    endif

    allocate( fact0_arr(nbexps,naexps), fact1_arr(nbexps,naexps), &
              fact2_arr(nbexps,naexps), cutoff(nbexps,naexps),    &
              stat=alloc_stat(2))
    MEMLOG(+size(fact0_arr)*4)
    ASSERT(alloc_stat(2).eq.0)
        alloc_stat(2)=1 ! fact0_arr fact1_arr fact2_arr
        alloc_stat(3)=1 ! cutoff

    fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
    fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

    where(fact0_arr>=very_small) ! prevent division by zero
       fact2_arr=fact1_arr/fact0_arr
    elsewhere
       fact2_arr=very_big
    end where

    where(fact2_arr*arg>options_integral_expmax()) ! cutoff: where almost no overlap
       cutoff=.false.
    elsewhere
       cutoff=.true.
    end where
    num=count(cutoff)

!-----------------------------------------------
!  if (integralpar_3cob_grad) then ! this is to zero all contyribs exept coul
!     do i_grad=1,gradient_data_n_spin_gradients
!          prim_int_3cob_grad(i_grad)%m = 0.0_r8_kind
!     enddo
!   endif
!----------------------------------------------------

    z_num: if(num==0) then ! all integrals are equal zero
      ! FIXME: Avoid special branches for num=0!
      !        Instead fix the main branch to handle that case
      !        properly!
      ! WARN('integrals screened')

      ! Assume these have been initialized to zero by allocate_primitives('initialize'):
      ! prim_int_2cob_poten, prim_int_3cob_solv_grad, prim_int_3cob_solv_grad_pc
     if(integralpar_3c_co) then
        ! FIXME: only importing for setting for zero
        prim_int_3c_co=0.0_r8_kind
     end if

     if (integralpar_3cob_grad) then

      do i_grad=1,gradient_data_n_spin_gradients
           prim_int_3cob_grad(i_grad)%m = 0.0_r8_kind
       if(integralpar_dervs) then
        do k2dr=1,gradient_data_n_gradients
         prim_int_coul_dervs(i_grad,k2dr)%m = 0.0_r8_kind
        enddo
       endif
      enddo

      spli: if(split_gradients) then
       mda_init:if (model_density) then
        do i_grad=1,gradient_data_n_gradients
           prim_int_3cob_coul_grad(i_grad)%m = 0.0_r8_kind
        end do
       endif mda_init

#ifndef no_cpks_coul_grads
      elseif(integralpar_cpksdervs) then
       do i_grad=1,gradient_data_n_gradients
           prim_int_cpks_coul_grad(i_grad)%m = 0.0_r8_kind
       end do
#endif

      endif spli
     endif

       MEMLOG(-size(fact0_arr)*4)
       deallocate(fact0_arr, fact1_arr, fact2_arr, cutoff,stat=alloc_stat(2))
       ASSERT(alloc_stat(2).eq.0)
       alloc_stat(3)=0 !cutoff
       !
       ! === RETURN POINT HERE ! ===
       !
       return
    end if z_num

    if(integralpar_gradients) then
       allocate(grad_mat(nbexps,naexps,2*lb+1,2*la+1,6),stat=alloc_stat(46))
          ASSERT(alloc_stat(46).eq.0)
                 alloc_stat(46)=1
       MEMLOG(size(grad_mat))
       grad_mat=0.0_r8_kind

     if(integralpar_dervs.or.ewpcdervs) then
       allocate(dervs_mat(nbexps,naexps,2*lb+1,2*la+1,6,6),stat=alloc_stat(72))
       MEMLOG(+size(dervs_mat))
          ASSERT(alloc_stat(72).eq.0)
                 alloc_stat(72)=1
       dervs_mat=0.0_r8_kind
     endif

     if(integralpar_cpksdervs.and.integralpar_3cob_grad) then
#ifdef no_cpks_coul_grads
       allocate( cpks_grad_mat(fit_coeff_n_ch(),6), &
#else
       allocate( cpks_grad_mat(nbexps,naexps,fit_coeff_n_ch(),2*lb+1,2*la+1,6), &
#endif
                 stat=cpksalloc(20))
          ASSERT(cpksalloc(20).eq.0)
          MEMLOG(+size(cpks_grad_mat))
          cpks_grad_mat=0.0_r8_kind
     endif

    endif

    if (integralpar_3cob_grad) then
     if (model_density .and. spin_polarized) then
        allocate( grad_mat_spin(nbexps,naexps,2*lb+1,2*la+1,6), stat=alloc_stat(64))
        MEMLOG(size(grad_mat_spin))
        ASSERT(alloc_stat(64).eq.0)
        alloc_stat(64)=1
        grad_mat_spin=0.0_r8_kind
     endif
    endif

    if (integralpar_2cob_nuc.or.integralpar_3cob_grad) then
       allocate(nuc(num,2*la+1,2*lb+1), stat=alloc_stat(6))
       MEMLOG(size(nuc))
       ASSERT(alloc_stat(6).eq.0)
              alloc_stat(6)=1
       nuc=0.0_r8_kind
    end if

    allocate (fact0(num),fact1(num),fact2(num), gamma_arg(num,3), &
              aexp_arr(num),bexp_arr(num), stat=alloc_stat(4))
    MEMLOG(+size(fact0)*8)
    ASSERT(alloc_stat(4).eq.0)
    alloc_stat(4)=1  ! fact0 fact1 fact2 gamma_arg
    alloc_stat(5)=1 ! aexp_arr bexp_arr


    fact0=pack(fact0_arr,cutoff)                 !  =>  a + b
    fact1=pack(fact1_arr,cutoff)                 !  =>  a * b
    fact2=pack(fact2_arr,cutoff)                 !  =>  a*b/(a+b)
    aexp_arr=pack(spread(aexps,1,nbexps),cutoff) !  => a
    bexp_arr=pack(spread(bexps,2,naexps),cutoff) !  => b
    MEMLOG(-size(fact0_arr)*3)
    deallocate(fact0_arr,fact1_arr,fact2_arr,stat=alloc_stat(2))
    ASSERT(alloc_stat(2).eq.0)

    if(integralpar_gradients) then
     allocate(two_aexp_arr(num),two_bexp_arr(num), &
                           two_fact2_vec(num,3),stat=alloc_stat(61))
     MEMLOG(+size(two_aexp_arr)*5)
     ASSERT(alloc_stat(61).eq.0)
     alloc_stat(61)=1
        two_aexp_arr=two*aexp_arr
        two_bexp_arr=two*bexp_arr
        two_fact2_vec(:,:)=spread(fact2(:),2,3)* &
                           spread(two*(xa(:)-xb(:)),1,num)
    endif
    ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
    ! d/da gamma_arg = a / (a + b) : aexp_arr/fact0
    do ix=1,3
       gamma_arg(:,ix)=(pack(spread(aexps*xa(ix),1,nbexps) + &
                             spread(bexps*xb(ix),2,naexps),cutoff))/fact0
    end do
    lmax_ch=0
    do i=1,n_unique_atoms
     l_max=max(unique_atoms(i)%lmax_ch,lmax_ch)
     lmax_ch=l_max
    enddo
    allocate(even_triangle(0:La,0:Lb,0:lmax_ch),stat=alloc_stat(7))
    MEMLOG(+size(even_triangle))
    ASSERT(alloc_stat(7).eq.0)

    do j1=0,la
     do j2=0,lb
      do j3=0,lmax_ch
       even_triangle(j1,j2,j3)=f_even_triangle(j1,j2,j3)
      enddo
     enddo
    enddo
    l_max=max(1,la,lb,lmax_ch)
    allocate(afac_sig(-la:la,-lb:lb,-lmax_ch:lmax_ch),stat=alloc_stat(7))
    MEMLOG(+size(afac_sig))
    ASSERT(alloc_stat(7).eq.0)
    alloc_stat(7)=1 ! even_triangle afac_sig

    do ma=-la,la
     do mb=-lb,lb
      do mc=-lmax_ch,lmax_ch
       afac_sig(ma,mb,mc)=1.0_r8_kind
        if(ma.gt.0) afac_sig(ma,mb,mc)=afac_sig(ma,mb,mc)*(-1)**ma
        if(mb.gt.0) afac_sig(ma,mb,mc)=afac_sig(ma,mb,mc)*(-1)**mb
        if(mc.gt.0) afac_sig(ma,mb,mc)=afac_sig(ma,mb,mc)*(-1)**mc
      enddo
     enddo
    enddo

    call csh_lm_map(l_max,"setup") ! that has to be made only once (latter) !

    allocate(csh_lm_of(0:l_max,-l_max:l_max),stat=alloc_stat(8))
    MEMLOG(+size(csh_lm_of))
    ASSERT(alloc_stat(8).eq.0)
    alloc_stat(8)=1
    call calc_csh_lm_of(l_max)

    call set_triplets()

    a=1; b=2; c=3

    ! CSH is used in PC contrib calcs as well
    allocate(CSH(3,3,(L_max+1)**2), &               ! CSH is used in both atomic and  PC contrib calcs
              rsh((L_max+1)**2), &                  ! temp in calculation of CSH products for one vector
              A_factor(3,(l_max+1)**2,0:l_max), &
    stat=alloc_stat(9)) ! actually more than needed !
    MEMLOG(size(CSH)+size(rsh)+size(A_factor))
    ASSERT(alloc_stat(9).eq.0)
    alloc_stat(9)=1

    grad_a: if(integralpar_gradients) then
            allocate(CSHG(3,3,(L_max+1)**2,3), &
                     CSHGb(3,3,(L_max+1)**2,3), &
                 rshg((L_max+1)**2,3), &
                 A_factorGa(3,(l_max+1)**2,0:l_max,3), &
                 A_factorGb(3,(l_max+1)**2,0:l_max,3), &
                 stat=alloc_stat(35))
        MEMLOG(+size(CSHG)*2+size(rshg)+size(A_factorGa)*2)
        ASSERT(alloc_stat(35).eq.0)
        alloc_stat(35)=1

     if(integralpar_dervs.or.ewpcdervs) then
      allocate(CSHGGa(3,3,(L_max+1)**2,3,3), &
               CSHGGb(3,3,(L_max+1)**2,3,3), &
               CSHGaGb(3,3,(L_max+1)**2,3,3), &
               A_factorGGa(3,(l_max+1)**2,0:l_max,3,3), &
               A_factorGGb(3,(l_max+1)**2,0:l_max,3,3), &
               A_factorGaGb(3,(l_max+1)**2,0:l_max,3,3), &
               rshgg((L_max+1)**2,3,3), stat=alloc_stat(65))
      MEMLOG(+size(CSHGGa)*3+size(A_factorGGa)*3+size(rshgg))
      ASSERT(alloc_stat(65).eq.0)
      alloc_stat(65)=1
    endif

    endif grad_a

    call calc_and_pack_3j()

           allocate(bin_fac(0:La+Lb,0:La+Lb),gamma_help_fac(num), r2_jjj_fac(num), &
                    gamma_help_fac_g(num),ghelp_shift_fac(num), stat=alloc_stat(11))
           MEMLOG(size(bin_fac)+size(gamma_help_fac)*4)
           ASSERT(alloc_stat(11).eq.0)
           alloc_stat(11)=1
           bin_fac=0.0_r8_kind

   csh_calcs: if(integralpar_2cob_nuc.or.integralpar_3c_co.or.integralpar_3cob_grad) then

    CSH(1:2,1:2,:)=czero
    CSH(a,b,:) = CSH_scalar(l_max,xa-xb)

    if(integralpar_3cob_grad) then

     CSHG (1:2,1:2,:,:)=czero ! this is (a) grad
     CSHG (3,3,:,:)=czero     ! this is (a) grad
     CSHGb(1:2,1:2,:,:)=czero
     CSHGb(3,3,:,:)=czero
     CSHG (a,b,:,:) =  CSHG_scalar(l_max,xa-xb)
     CSHGb(a,b,:,:) = -CSHG(a,b,:,:)
     A_factorGa=czero
     A_factorGb=czero

    if(integralpar_dervs.or.ewpcdervs) then

     CSHGGa(:2,:2,:,:,:)=czero
     CSHGGa(3,3,:,:,:)=czero
     CSHGGb(:2,:2,:,:,:)=czero
     CSHGGb(3,3,:,:,:)=czero
     CSHGaGb(:2,:2,:,:,:)=czero
     CSHGaGb(3,3,:,:,:)=czero

!!     CSHGGa(a,b,:,:,:) = CSHGG_scalar(l_max,xa-xb,i_p=1)
!!     variant for internal debug printing

!    according numerical test
!    both (GGa) and (GGb)  solid harmonic derivatives below have wrong sign
!    GaGb should be also of oposit sign

     CSHGGa(a,b,:,:,:) =  CSHGG_scalar(l_max,xa-xb)
!     CSHGGb(a,b,:,:,:) = -CSHGGa(a,b,:,:,:)   ! in two center integrals GGa and GGB should have tha same values
     CSHGGb(a,b,:,:,:) = CSHGGa(a,b,:,:,:)

!     CSHGaGb(a,b,:,:,:) = CSHGGa(a,b,:,:,:)
     CSHGaGb(a,b,:,:,:) = -CSHGGa(a,b,:,:,:)   ! in two center integrals GaGb=-GGa

!--------------------------------------------
!    print*, 'aa CSH 0 1 2',   sum(CSH(a,b,:)),sum(CSHG(a,b,:,1)),sum(CSHGGa(a,b,:,1,1))
!    print*, 'aa CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHG(a,b,:,2)),sum(CSHGGa(a,b,:,2,1))
!    print*, 'aa CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHG(a,b,:,3)),sum(CSHGGa(a,b,:,3,1))
!
!    print*, 'bb CSH 0 1 2',   sum(CSH(a,b,:)),sum(CSHGb(a,b,:,1)),sum(CSHGGb(a,b,:,1,1))
!    print*, 'bb CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHGb(a,b,:,2)),sum(CSHGGb(a,b,:,2,1))
!    print*, 'bb CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHGb(a,b,:,3)),sum(CSHGGb(a,b,:,3,1))
!
!    print*, 'a ab CSH 0 1 2',   sum(CSH(a,b,:)),sum(CSHG(a,b,:,1)),sum(CSHGaGb(a,b,:,1,1))
!    print*, 'a ab CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHG(a,b,:,2)),sum(CSHGaGb(a,b,:,2,1))
!    print*, 'a ab CSH21 0 1 2', sum(CSH(a,b,:)),sum(CSHG(a,b,:,3)),sum(CSHGaGb(a,b,:,3,1))
!----------------------------------------

    endif
    endif


          do lm=1,(l_max+1)**2
             CSH(b,a,lm) = (-1)**csh_map%l(lm)*CSH(a,b,lm)
          enddo

       if(integralpar_3cob_grad) then
          do lm=1,(l_max+1)**2
           CSHG(b,a,lm,:) = (-1)**csh_map%l(lm)*CSHG(a,b,lm,:)
           CSHGb(b,a,lm,:) = (-1)**csh_map%l(lm)*CSHGb(a,b,lm,:)
          enddo
          if(integralpar_dervs.or.ewpcdervs) then
           do lm=1,(l_max+1)**2
            CSHGGa(b,a,lm,:,:)=(-1)**csh_map%l(lm)*CSHGGa(a,b,lm,:,:)
            CSHGGb(b,a,lm,:,:)=(-1)**csh_map%l(lm)*CSHGGb(a,b,lm,:,:)
            CSHGaGb(b,a,lm,:,:)=(-1)**csh_map%l(lm)*CSHGaGb(a,b,lm,:,:)
           enddo
          endif
       endif

!    print*,'CSHGGb ba', sum(CSHGb(b,a,:,1)),sum(CSHGGb(b,a,:,1,1))
!    print*,'CSHGGb ab', sum(CSHGb(a,b,:,1)),sum(CSHGGb(a,b,:,1,1))
!    do lm=1,size(CSHGb,3)
!     print*,lm,CSHGb(b,a,lm,1),CSHGGb(b,a,lm,1,1)
!     print*,lm,CSHGb(a,b,lm,1),CSHGGb(a,b,lm,1,1)
!    enddo

    endif csh_calcs


      allocate(temp_overlap_nuc(num),gamma_argXC(num,3), stat=alloc_stat(12))
      MEMLOG(size(temp_overlap_nuc)*4)
       ASSERT(alloc_stat(12).eq.0)
       alloc_stat(12)=1

      temp_overlap_nuc = exp(-fact2*arg) / &
                    (fact0*sqrt(aexp_arr**La*df(2*La-1)*bexp_arr**Lb*df(2*Lb-1) )) &
                                    * 2*pi * (4*fact1(:)/pi2)**0.75_r8_kind
    if(integralpar_gradients) then
     allocate( radang_temp_gb(num,(2*La+1),(2*Lb+1),3), &
               radang_temp_ga(num,(2*La+1),(2*Lb+1),3), &
               Gb(num,3),Ga(num,3),GtFvec(num,2,3), &
               stat=alloc_stat(77))
      MEMLOG(size(radang_temp_gb)*2+size(Ga)*4)
      ASSERT(alloc_stat(77).eq.0)
      alloc_stat(77)=1
    if(integralpar_dervs) then
      allocate(radang_temp_gga(num,(2*La+1),(2*Lb+1),3,3),  &
               radang_temp_ggb(num,(2*La+1),(2*Lb+1),3,3),  &
               radang_temp_gagb(num,(2*La+1),(2*Lb+1),3,3), &
               GGt(num,(2*La+1),(2*Lb+1)), &
               GtGa(num,(2*La+1),(2*Lb+1),3), &
               GtGb(num,(2*La+1),(2*Lb+1),3), &
               GGa(num,3,3),GGb(num,3,3),GaGb(num,3,3), &
               stat=alloc_stat(76))
      MEMLOG(+size(radang_temp_gga)*3+size(GGt)*7+size(GGa)*3)
      ASSERT(alloc_stat(76).eq.0)
      alloc_stat(76)=1
    endif
    endif

!   print*, 'CALCULATE NUCLEAR ATRACTION AND COULOMB INTEGRALS    FOR ATOMS'

    nucfit_calc: if(integralpar_2cob_nuc.or.integralpar_3c_co.or.integralpar_3cob_grad) then


    ! NOW LOOP OVER THIRD INDEX for atoms
    nc = 0
    unique_atoms_fit_nuc : do i=1,n_unique_atoms
       uni_c=i

       ua_pointer => unique_atoms(i)
       lmax_ch= ua_pointer%lmax_ch
       n_equals = ua_pointer%n_equal_atoms

       z= ua_pointer%z                  ! charge
       if(ua_pointer%zc.ne.0.0_r8_kind) then
        z=0.0_r8_kind
       endif

       nc_fit_only=nc

 lcfit_calc: if(integralpar_3c_co.or.integralpar_3cob_grad) then

    if(integralpar_3c_co.or.integralpar_3cob_grad) then
        allocate(coul_int(-1:lmax_ch), stat=alloc_stat(13))
        ASSERT(alloc_stat(13).eq.0)
        alloc_stat(13)=1
    endif

       setup_gr: if(integralpar_3cob_grad) then

          imc = unique_atoms(i)%moving_atom
          moving_c = imc > 0

          if (moving_c)then
             grad_dim = gradient_index(imc+1) - gradient_index(imc)

          allocate (rotmat_eq(grad_dim, 3, n_equals), do_rotation_eq(n_equals), &
                   gamma_arg_xc(num,3,n_equals), stat=alloc_stat(62))
                     !used in  Radial_Ang_3cSA,  Radial_Ang_3cS
          MEMLOG(size(rotmat_eq)+size(do_rotation_eq)+size(gamma_arg_xc))
             ASSERT(alloc_stat(62).eq.0)
             alloc_stat(62)=1
             alloc_stat(63)=1

             do_rotation_eq = .false.
             if(grad_dim.gt.0) rotmat_eq=unique_atom_grad_info(imc)%m
          else
           allocate(gamma_arg_xc(num,3,n_equals),stat=alloc_stat(63))
           MEMLOG(size(gamma_arg_xc))
             ASSERT(alloc_stat(63).eq.0)
           alloc_stat(63)=1
             grad_dim = 0
          endif

          call calc_3c_codervs_setup(i,la,lb)
                                    !^ ca_dervs_mat is zero init for each unique (i)

         endif setup_gr

      if(integralpar_3cob_grad.and.integralpar_dervs) then
          allocate(gamma_arg2_vec(num,n_equals),temp_overlap(num), &
                   radang_temp(num,(2*La+1),(2*Lb+1)),c_exp_arg2(num), &
                   radangF(num), &
                   gamma_help_r2(num,La+Lb+4,n_equals), &
                   dervs_totsymM(num,2*Lb+1,2*La+1,grad_dim,grad_dim), &
                   dervs_totsymM_temp(num,2*Lb+1,2*La+1,grad_dim,3), &
                   dervsM(num,2*Lb+1,2*La+1,6,6), &
                   ca_dervsM(num,2*Lb+1,2*La+1,grad_dim,6), &
                                            stat=alloc_stat(15))
          MEMLOG(size(gamma_arg2_vec)+size(temp_overlap)*3)
          MEMLOG(size(radang_temp)+size(gamma_help_r2))
          MEMLOG(size(dervs_totsymM)+size(dervs_totsymM_temp))
          MEMLOG(size(dervsM)+size(ca_dervsM))
          ASSERT(alloc_stat(15).eq.0)
          alloc_stat(55)=1

      elseif(integralpar_3cob_grad) then
          allocate(gamma_arg2_vec(num,n_equals),temp_overlap(num), &
                   radang_temp(num,(2*La+1),(2*Lb+1)),c_exp_arg2(num), &
                   gamma_help_r2(num,La+Lb+3,n_equals), &
                   radangF(num), &
                   stat=alloc_stat(15)) ! for atoms (1)
          MEMLOG(size(gamma_arg2_vec)+size(temp_overlap)*3)
          MEMLOG(size(radang_temp)+size(gamma_help_r2))
          ASSERT(alloc_stat(15).eq.0)
          alloc_stat(55)=1 ! radangF
      else
          allocate(gamma_arg2_vec(num,n_equals),temp_overlap(num), &
                   radang_temp(num,(2*La+1),(2*Lb+1)),c_exp_arg2(num), &
                   gamma_help_r2(num,La+Lb+2,n_equals), &
                   stat=alloc_stat(15))
          MEMLOG(size(gamma_arg2_vec)+size(temp_overlap)*2)
          MEMLOG(size(radang_temp)+size(gamma_help_r2))
          ASSERT(alloc_stat(15).eq.0)
      endif
          alloc_stat(15)=1 ! gamma_arg2_vec
          alloc_stat(16)=1 ! temp_overlap
          alloc_stat(34)=1 ! radang_temp c_exp_arg2




      temp_overlap = exp(-fact2*arg) / &
                    (fact0*sqrt(aexp_arr**La*df(2*La-1)*bexp_arr**Lb*df(2*Lb-1) )) &
                                    * 2*pi2*sqrt(pi)  *  (4*fact1/pi2)**0.75_r8_kind

      lc_fit: do lc=0,lmax_ch

           n_jj=n_jj_terms(la,lb,lc)
           allocate(alph_pw(num,0:alph_pw_max),beta_pw(num,0:beta_pw_max) &
                   ,fact2_pw(num,0:fact2_pw_max), &
                    mup_llam(0:La+Lb+Lc),i3_fac(num,0:la+lb,0:Lc+La,0:Lb+Lc), &
                    i3_func_fac(num),                      stat=alloc_stat(17))
           MEMLOG(size(alph_pw)+size(beta_pw)+size(fact2_pw)+size(mup_llam)+size(i3_func_fac))
           MEMLOG(size(i3_fac))
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

        llam_m=llam_min(la,lb,lc)

                do m_up=0,la+lb
                 i3_func_fac(:)=temp_overlap(:)*fact2_pw(:,m_up)
                 if(overlap_var_fixed.and.integralpar_gradients)  &
                                                   i3_func_fac(:)=1.0_r8_kind
                   do i1=0,Lc+La
                   do i2=0,Lb+Lc
                   i3_fac(:,m_up,i1,i2)=i3_func_fac(:) *    alph_pw(:,i1) * beta_pw(:,i2)
                  enddo
                  enddo
                  enddo

           ncexps = ua_pointer%l_ch(lc)%n_exponents
           cexps => unique_atoms(i)%l_ch(lc)%exponents(:)

     if(integralpar_3cob_grad.and.integralpar_dervs) then
       !the only difference from gradients +2 but not +1 N_max increment

       N_max = la + lb + lc + 2  ! not used for r2_radial

       allocate ( gamma_help(num,1+N_max),  &                  ! not 3c procs
                  gamma_help_co(num,1+N_max,n_equals), &       ! 0:lc radial
                               !      ^
                  exp_arg(num,0:max(La+Lb+Lc,1)),gamma_arg2(num), &
                  jj_fac_SA(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                        !^   radial energy and up to second dervs temp, r2_radial
                  gamma_help_g(num,La+Lb+1,n_equals), &        ! used in radial_r2  energy
                  gamma_help_gG(num,La+Lb+1,n_equals), &       ! used in radial_r2  grads
                  gamma_help_gGG(num,La+Lb+1,n_equals), &      ! used in radial_r2  2nd dervs
                  jj_facG(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                       !^   -                                  ! first derivs r2_radial  temp
                  jj_facGG(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                       ! ^                                     ! 2nd dervs r2_radial temp
                  radial3cmat(num,n_jj,n_equals),  &           ! radial r2_radial
                  AngR3cMat(n_jj,2*la+1,2*lb+1,2*lc+1,n_equals), &
                  Gt(num,(2*La+1),(2*Lb+1)),&                  ! temp for rad_ang dervs routes
                                        stat=alloc_stat(18))
     MEMLOG(size(gamma_help)+size(gamma_help_co)+size(exp_arg)+size(gamma_arg2))
     MEMLOG(size(jj_fac_SA)+size(gamma_help_g)*3+size(jj_facG)*2)
     MEMLOG(size(radial3cmat)+size(AngR3cMat)+size(Gt))
     elseif(integralpar_3cob_grad) then
       N_max = la + lb + lc + 1
       allocate ( gamma_help(num,1+N_max), &             ! used to calc nuc
                  gamma_help_co(num,1+N_max,n_equals), &
                  exp_arg(num,0:max(La+Lb+Lc,1)),gamma_arg2(num), &
                  jj_fac_SA(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                  gamma_help_g(num,La+Lb+1,n_equals), & ! used in radial_r2
                  gamma_help_gG(num,La+Lb+1,n_equals), &
                  jj_facG(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                  radial3cmat(num,n_jj,n_equals),  &
                  AngR3cMat(n_jj,2*la+1,2*lb+1,2*lc+1,n_equals), &
                  Gt(num,(2*La+1),(2*Lb+1)),&                  ! temp for rad_ang dervs routes
                                                          stat=alloc_stat(18))
    MEMLOG(size(gamma_help)+size(gamma_help_co)+size(exp_arg)+size(gamma_arg2))
    MEMLOG(size(jj_fac_SA)+size(gamma_help_g)*2+size(jj_facG))
    MEMLOG(size(radial3cmat)+size(AngR3cMat)+size(Gt))
     else ! i.e. no derivs and grads
       N_max = la + lb + lc
       allocate ( gamma_help(num,1+N_max), &                 ! used to calc nuc
                  gamma_help_co(num,1+N_max,n_equals), &
                  exp_arg(num,0:max(La+Lb+Lc,1)),gamma_arg2(num), &
                  jj_fac_SA(num,0:la+lb,llam_m:la+lb+lc,n_equals), &
                        !^   - radial and r2_radial energy contris temp
                  gamma_help_g(num,La+Lb+1,n_equals), & ! used in radial_r2
                  radial3cmat(num,n_jj,n_equals),  &
                  AngR3cMat(n_jj,2*la+1,2*lb+1,2*lc+1,n_equals), &
                               !          ^
                                                          stat=alloc_stat(18))
    MEMLOG(size(gamma_help)+size(gamma_help_co)+size(exp_arg)+size(gamma_arg2))
    MEMLOG(size(jj_fac_SA)+size(gamma_help_g))
    MEMLOG(size(radial3cmat)+size(AngR3cMat))
     endif
                ASSERT(alloc_stat(18).eq.0)
                alloc_stat(18)=1 ! gamma_help gamma_help_co exp_arg gamma_arg2 jj_fac_SA
                alloc_stat(19)=1 ! radial3cmat
                alloc_stat(20)=1 ! AngR3cMat
                AngR3cMat = 0.0_r8_kind

      grads_AngR: if(integralpar_3cob_grad) then
       if(integralpar_dervs) then
        allocate(radial3cmatG(num,n_jj,n_equals), radial3cmatGG(num,n_jj,n_equals), &
                                                ! ^ radial and r2_radial second dervs
                AngR3cMatGa(n_jj,2*La+1,2*Lb+1,3,2*Lc+1,n_equals), &
                AngR3cMatGb(n_jj,2*La+1,2*Lb+1,3,2*Lc+1,n_equals), &
                AngR3cMatGGa(n_jj,2*La+1,2*Lb+1,3,3,2*Lc+1,n_equals), &
                AngR3cMatGGb(n_jj,2*La+1,2*Lb+1,3,3,2*Lc+1,n_equals), &
                AngR3cMatGaGb(n_jj,2*La+1,2*Lb+1,3,3,2*Lc+1,n_equals), &
                                                   stat=alloc_stat(50))

                ASSERT(alloc_stat(50).eq.0)
                MEMLOG(size(AngR3cMatGa)*2+size(AngR3cMatGGa)*3+size(radial3cmatG)*2)
                alloc_stat(50)=1

                AngR3cMatGGa=0.0_r8_kind
                AngR3cMatGGb=0.0_r8_kind
                AngR3cMatGaGb=0.0_r8_kind

       else
        allocate(radial3cmatG(num,n_jj,n_equals), &
                !           ^                          ! - radial and r2_radial  first dervs
                AngR3cMatGa(n_jj,2*La+1,2*Lb+1,3,2*Lc+1,n_equals), &
                AngR3cMatGb(n_jj,2*La+1,2*Lb+1,3,2*Lc+1,n_equals), &
                                     stat=alloc_stat(50))
                ASSERT(alloc_stat(50).eq.0)
                MEMLOG(size(AngR3cMatGa)*2+size(radial3cmatG))
       endif
                alloc_stat(50)=1 ! radial3cmatG radial3cmatGG
                alloc_stat(49)=1 ! AngR3cMatGa AngR3cMatGb
                AngR3cMatGa=0.0_r8_kind
                AngR3cMatGb=0.0_r8_kind
      endif grads_AngR


       lccase: if(lc.gt.0) then

       ! complex spherical harmonics and their derivatives
       ! are calculated separately for lc.eq.0 and lc.gt.0 cases

       lfit: if(new_3c_co) then

        allocate(m_c_required(-la:la,-lb:lb,-lc:lc,n_equals), &
                 mc_required(2*la+1,2*lb+1,2*lc+1,n_equals),stat=alloc_stat(21))
            ASSERT(alloc_stat(21).eq.0)
            MEMLOG(size(m_c_required)+size(mc_required))
                   alloc_stat(21)=1
        mc_required=.false.
        m_c_required=.false.
        call calc_mc_required(La,Lb,Lc)
    independents: if(n_independent_fcts.ne.0) then

       call calc_3c_colc_setup(la,lb)

       if(l_time) call start_timer(timer_ang3c)

       nc=nc_fit_only

       equal_atoms_angl  : do j=1,n_equals

         CSH(1:2,3,:)=czero
         CSH(3,1:2,:)=czero
         CSH(3,3,:)=czero

          A_factor = czero

          xc=unique_atoms(i)%position(:,j)
          nc=nc+1

          gamma_argXC(:,1)=gamma_arg(:,1)-xc(1)
          gamma_argXC(:,2)=gamma_arg(:,2)-xc(2)
          gamma_argXC(:,3)=gamma_arg(:,3)-xc(3)

          if(gamma_var_fixed.and.integralpar_gradients) gamma_argXC=1.0_r8_kind

          gamma_arg2_vec(:,j)=( gamma_argXC(:,1)**2 + &
                                gamma_argXC(:,2)**2 + &
                                gamma_argXC(:,3)**2   ) * fact0

          CSH(b,c,:) = CSH_scalar(l_max,xb-xc)
          CSH(c,a,:) = CSH_scalar(l_max,xc-xa)
          do lm=1,(l_max+1)**2
             lp = csh_map % l(lm)
             CSH(c,b,lm) = (-1)**lp * CSH(b,c,lm)
             CSH(a,c,lm) = (-1)**lp * CSH(c,a,lm)
          end do
          call calculate_factor(1,b,c,a,la)    ! A
          call calculate_factor(2,c,a,b,lb)    ! B
          call calculate_factor(3,a,b,c,lc)    ! C3c

           allocate(AngCMat(n_jj,-la:la, -lb:lb, -lc:lc), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
              alloc_stat(22)=1
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("calc_3c AngCMat allocation failed")
           endif

       gr1: if(integralpar_3cob_grad) then
        CSHG(b,c,:,:) = czero
        CSHGb(b,c,:,:) = CSHG_scalar(l_max,xb-xc)
        CSHGb(c,a,:,:) = czero
        CSHG(c,a,:,:) = -CSHG_scalar(l_max,xc-xa)
         do lm=1,(l_max+1)**2
            CSHG(c,b,lm,:) = (-1)**csh_map%l(lm) * CSHG(b,c,lm,:)
            CSHG(a,c,lm,:) = (-1)**csh_map%l(lm) * CSHG(c,a,lm,:)
            CSHGb(c,b,lm,:) = (-1)**csh_map%l(lm) * CSHGb(b,c,lm,:)
            CSHGb(a,c,lm,:) = (-1)**csh_map%l(lm) * CSHGb(c,a,lm,:)
         end do
          call calculate_factor_Ga(1,b,c,a,la)    ! A
          call calculate_factor_Ga(2,c,a,b,lb)    ! B
          call calculate_factor_Ga(3,a,b,c,lc)    ! C

          call calculate_factor_Gb(1,b,c,a,la)    ! A
          call calculate_factor_Gb(2,c,a,b,lb)    ! B
          call calculate_factor_Gb(3,a,b,c,lc)    ! C

          dr2els: if(integralpar_dervs) then
           allocate(AngCMatGa(n_jj,-la:la,-lb:lb,-lc:lc,3), &
                    AngCMatGb(n_jj,-la:la,-lb:lb,-lc:lc,3), &
                    AngCMatGGa(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                    AngCMatGGb(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                    AngCMatGaGb(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                                         stat=alloc_stat(48))
           MEMLOG(+size(AngCMatGa)*4+size(AngCMatGGa)*6)
           if(alloc_stat(48).eq.0) then
           alloc_stat(48)=1

           AngCMatGGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGGb = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGaGb = (0.0_r8_kind,0.0_r8_kind)

           AngCMatGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGb = (0.0_r8_kind,0.0_r8_kind)

           else
            call error_handler("48 1 AngCMat allocation failed")
           endif

          else dr2els
           allocate(AngCMatGa(n_jj,-la:la,-lb:lb,-lc:lc,3), &
                    AngCMatGb(n_jj,-la:la,-lb:lb,-lc:lc,3), stat=alloc_stat(48))
           MEMLOG(+size(AngCMatGa)*4)
           if(alloc_stat(48).eq.0) then
           alloc_stat(48)=1

           AngCMatGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGb = (0.0_r8_kind,0.0_r8_kind)
           else
            call error_handler("48 2 AngCMat allocation failed")
           endif
          endif dr2els

          lcdervs: if(integralpar_dervs) then
           CSHGGa(b,c,:,:,:) = czero ! (b,c) element do not depend on a
           CSHGGa(c,b,:,:,:) = czero


           CSHGGb(a,c,:,:,:) = czero ! (a,c) element do not depend on b
           CSHGGb(c,a,:,:,:) = czero

           CSHGaGb(a,c,:,:,:) = czero ! (a,c) element do not depend on b
           CSHGaGb(c,a,:,:,:) = czero

           CSHGaGb(c,b,:,:,:) = czero ! (c,b) element do not depend on a
           CSHGaGb(b,c,:,:,:) = czero

!          CSHGGa(c,a,:,:,:) = -CSHGG_scalar(l_max,xc-xa) !for second deriv (+) as in the next line
           CSHGGa(c,a,:,:,:) =  CSHGG_scalar(l_max,xc-xa)
           CSHGGb(c,b,:,:,:) =  CSHGG_scalar(l_max,xc-xb)


           do lm=1,(l_max+1)**2
            CSHGGa(a,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGa(c,a,lm,:,:)
            CSHGGb(b,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGb(c,b,lm,:,:)
!            print*,'j lm CSHGGa',i,lm,CSHG(a,c,lm,1),CSHGGa(a,c,lm,1,1)
!            print*,'j lm CSHGGb',i,lm,CSHGb(b,c,lm,1),CSHGGb(b,c,lm,1,1)
           enddo
!         print*,'i sum CSH11',i,sum(CSH(a,c,:)),sum(CSHG(a,c,:,1)),sum(CSHGGa(a,c,:,1,1))
!         print*,'i sum CSH11',i,sum(CSH(b,c,:)),sum(CSHGb(b,c,:,1)),sum(CSHGGb(b,c,:,1,1))

           call calculate_factor_GGa(1,b,c,a,la)    ! A
           call calculate_factor_GGa(2,c,a,b,lb)    ! B
           call calculate_factor_GGa(3,a,b,c,lc)    ! C

          ! similar statements for GGb counterpart follow

           call calculate_factor_GGb(1,b,c,a,la)    ! A
           call calculate_factor_GGb(2,c,a,b,lb)    ! B
           call calculate_factor_GGb(3,a,b,c,lc)    ! C

          ! now GaGb and GbGa cross therms should follow
           call calculate_factor_GaGb(1,b,c,a,la)    ! A
           call calculate_factor_GaGb(2,c,a,b,lb)    ! B
           call calculate_factor_GaGb(3,a,b,c,lc)    ! C

          endif lcdervs

       endif gr1

           call calc_angular_cmplx3cSA(la,lb,lc,AngCMat,n_jj) !lc>0

!-----------------------------------------------------
!        if(integralpar_2dervs.and.integralpar_gradients) then
!         do k=1,size(AngCMatGa(:,:,:,:,1),4)
!           print*,'AngCMatGGa',k,sum(AngCMatGa(:,:,:,k,1)),sum(AngCMatGGa(:,:,:,k,1,1))
!         enddo
!           print*,'AngCMatGGa',sum(AngCMatGa(:,:,:,:,1)),sum(AngCMatGGa(:,:,:,:,1,1))
!         endif
!---------------------------------------------------------

           AngRMat=>AngR3cMat(:,:,:,:,j)

           if(integralpar_3cob_grad) then
             AngRMatGa=>AngR3cMatGa(:,:,:,:,:,j)
             AngRMatGb=>AngR3cMatGb(:,:,:,:,:,j)

             if(integralpar_dervs) then
              AngRMatGGa=>AngR3cMatGGa(:,:,:,:,:,:,j)
              AngRMatGGb=>AngR3cMatGGb(:,:,:,:,:,:,j)
              AngRMatGaGb=>AngR3cMatGaGb(:,:,:,:,:,:,j)
             endif

            endif

            call Ang_Cmplx_to_RealSA(la,lb,lc)  ! AngRMatGGa added
                                                ! AngRMatGGb addeed
                                                ! AngRMatGaGb added

        if(harmonic_var_fixed.and.integralpar_gradients) then
           AngRMat=1.0_r8_kind     !!!!!
           AngRMatGa=0.0_r8_kind   !!!!!
           AngRMatGb=0.0_r8_kind   !!!!!
        endif

            MEMLOG(-size(AngCMat))
            deallocate(AngCMat, stat=alloc_stat(22))
            ASSERT(alloc_stat(22).eq.0)

        if(integralpar_3cob_grad) then

         if(integralpar_dervs) then
          MEMLOG(-size(AngCMatGa)*4-size(AngCMatGGa)*6)
          deallocate(AngCMatGGa,AngCMatGGb,AngCMatGaGb, &
                     AngCMatGa,AngCMatGb,stat=alloc_stat(48))
         else
          MEMLOG(-size(AngCMatGa)*4)
          deallocate(AngCMatGa,AngCMatGb,stat=alloc_stat(48))
         endif
         ASSERT(alloc_stat(48).eq.0)

          gamma_arg_xc(:,1,j)=gamma_arg(:,1)-unique_atoms(i)%position(1,j)
          gamma_arg_xc(:,2,j)=gamma_arg(:,2)-unique_atoms(i)%position(2,j)
          gamma_arg_xc(:,3,j)=gamma_arg(:,3)-unique_atoms(i)%position(3,j)

          if (moving_c) then  ! equal_atoms_angl
             ! If is  to guard against imc  < 0.  FIXME:  no reason to
             ! supply symmetry info only for moving atoms:
             do_rotation_eq(j) = .not. unity_matrix_p (unique_atom_grad_info(imc) % m(:, :, j))
          endif
       endif

         enddo equal_atoms_angl

         if(l_time) call stop_timer(timer_ang3c)
         if(l_time) call start_timer(timer_rad3c)

       radangT: if(integralpar_3cob_grad) then

        allocate(grad_totsymM(num,2*Lb+1,2*La+1,grad_dim), &
               gradM(num,2*Lb+1,2*La+1,6), &
               stat=alloc_stat(10))
        MEMLOG(size(grad_totsymM)+size(gradM))
               ASSERT(alloc_stat(10).eq.0)
                      alloc_stat(10)=1
       endif radangT

           pointer_coul => coul_int(lc)%m

        do k=1,ncexps

            c_exp_arg2(:)=(cexps(k)*sqrt(fact0(:)+cexps(k)))

                                      call calc_radial_3cSA(la,lb,lc,N_max) ! 1
            if(integralpar_3cob_grad) call calc_radial_3cGa(la,lb,lc)
            if(integralpar_dervs)    call calc_radial_3cGGa(la,lb,lc)

            call Radial_Ang_3cSA(la,lb,lc)   ! only call to this route

        enddo

       if(integralpar_3cob_grad) then
        MEMLOG(-size(grad_totsymM)-size(gradM))
        deallocate( grad_totsymM,gradM, stat=alloc_stat(10))
        ASSERT(alloc_stat(10).eq.0)
       endif

          if(l_time) call stop_timer(timer_rad3c)
       endif independents

            MEMLOG(-size(m_c_required)-size(mc_required))
         deallocate(mc_required,m_c_required,stat=alloc_stat(21))
         ASSERT(alloc_stat(21).eq.0)
       endif lfit

       else lccase ! i.e. now lc=0

           ! complex spherical harmonics and their derivatives
           ! are calculated separately for lc.eq.0 and lc.gt.0 cases

        ! nuc is also calc in lccase thus cal nucR_fac for radial routes

        allocate(nucR_fac(num,0:la+lb,0:Lc+La,0:Lb+Lc),stat=alloc_stat(25))
        ASSERT(alloc_stat(25).eq.0)
        alloc_stat(25)=1
        MEMLOG(size(nucR_fac))

                 do m_up=0,la+lb
                 i3_func_fac(:)=temp_overlap_nuc(:)*fact2_pw(:,m_up)
                 do i1=0,Lc+La
                 do i2=0,Lb+Lc
                  nucR_fac(:,m_up,i1,i2)=i3_func_fac(:) *    alph_pw(:,i1) * beta_pw(:,i2)
                 enddo
                 enddo
                 enddo



       nc=nc_fit_only
       if(l_time) call start_timer(timer_ang3c)

       nolc_equal_atoms  : do j=1,n_equals
!       equal_c=j


          CSH      = czero
          A_factor = czero
!          xc=unique_atoms(i)%position(:,equal_c)
          xc=unique_atoms(i)%position(:,j)
          nc=nc+1

          gamma_argXC(:,1)=gamma_arg(:,1)-xc(1)
          gamma_argXC(:,2)=gamma_arg(:,2)-xc(2)
          gamma_argXC(:,3)=gamma_arg(:,3)-xc(3)

          if(gamma_var_fixed.and.integralpar_gradients) gamma_argXC=1.0_r8_kind

          ! gamma_arg2_vec used in calc_radial
          gamma_arg2_vec(:,j)=( gamma_argXC(:,1)**2 + &
                                gamma_argXC(:,2)**2 + &
                                gamma_argXC(:,3)**2   ) * fact0

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
          call calculate_factor(3,a,b,c,lc )  ! C3c

       nolc_derivatives: if(integralpar_3cob_grad) then
          gamma_arg_xc(:,1,j)=gamma_arg(:,1)-unique_atoms(i)%position(1,j)
          gamma_arg_xc(:,2,j)=gamma_arg(:,2)-unique_atoms(i)%position(2,j)
          gamma_arg_xc(:,3,j)=gamma_arg(:,3)-unique_atoms(i)%position(3,j)

          if (moving_c) then ! equal_atoms_lfit
             ! If is  to guard against imc  < 0.  FIXME:  no reason to
             ! supply symmetry info only for moving atoms:
             do_rotation_eq(j) = .not. unity_matrix_p (unique_atom_grad_info(imc) % m(:, :, j))
          endif

        CSHG(b,c,:,:)  = czero ! (b,c) element do not depend on a
        CSHGb(c,a,:,:) = czero
        CSHGb(b,c,:,:) = CSHG_scalar(l_max,xb-xc)
        CSHG(c,a,:,:) = -CSHG_scalar(l_max,xc-xa)

         do lm=1,(l_max+1)**2
            lp = csh_map % l(lm)
            CSHG(c,b,lm,:) = (-1)**lp * CSHG(b,c,lm,:)
            CSHG(a,c,lm,:) = (-1)**lp * CSHG(c,a,lm,:)
            CSHGb(c,b,lm,:) = (-1)**lp * CSHGb(b,c,lm,:)
            CSHGb(a,c,lm,:) = (-1)**lp * CSHGb(c,a,lm,:)
         end do

          call calculate_factor_Ga(1,b,c,a,la)    ! A
          call calculate_factor_Ga(2,c,a,b,lb)    ! B
          call calculate_factor_Ga(3,a,b,c,lc)    ! C

          call calculate_factor_Gb(1,b,c,a,la)    ! A
          call calculate_factor_Gb(2,c,a,b,lb)    ! B
          call calculate_factor_Gb(3,a,b,c,lc)    ! C

          allocate(AngCMatGa(n_jj,-la:la,-lb:lb,-lc:lc,3), &
                   AngCMatGb(n_jj,-la:la,-lb:lb,-lc:lc,3), stat=alloc_stat(48))
          MEMLOG(+size(AngCMatGa)*4)
          if(alloc_stat(48).eq.0) then
          alloc_stat(48)=1
           AngCMatGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGb = (0.0_r8_kind,0.0_r8_kind)
          else
           call error_handler("48 AngCMat allocation failed")
          endif

         nolc_2dervs: if(integralpar_dervs) then
           allocate(AngCMatGGa(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                    AngCMatGGb(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                    AngCMatGaGb(n_jj,-la:la,-lb:lb,-lc:lc,3,3), &
                    stat=alloc_stat(48))
                    MEMLOG(+size(AngCMatGGa)*6)
                    ASSERT(alloc_stat(48).eq.0)
                           alloc_stat(48)=1

           AngCMatGGa = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGGb = (0.0_r8_kind,0.0_r8_kind)
           AngCMatGaGb = (0.0_r8_kind,0.0_r8_kind)

           CSHGGa(b,c,:,:,:) = czero ! (b,c) elements do not depend on a
           CSHGGa(c,b,:,:,:) = czero

           CSHGGb(a,c,:,:,:) = czero ! (a,c) element do not depend on b
           CSHGGb(c,a,:,:,:) = czero

           CSHGaGb(a,c,:,:,:) = czero ! (a,c) do not depend on b
           CSHGaGb(c,a,:,:,:) = czero

           CSHGaGb(b,c,:,:,:) = czero ! (b,c) do not depend on a
           CSHGaGb(c,b,:,:,:) = czero

           CSHGGa(c,a,:,:,:) =  CSHGG_scalar(l_max,xc-xa)
           CSHGGb(c,b,:,:,:) =  CSHGG_scalar(l_max,xc-xb)
!           CSHGGb(c,b,:,:,:) =  CSHGG_scalar(l_max,xc-xb,i_p=3)

           do lm=1,(l_max+1)**2
            CSHGGa(a,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGa(c,a,lm,:,:)
            CSHGGb(b,c,lm,:,:) = (-1)**csh_map%l(lm) * CSHGGb(c,b,lm,:,:)
           enddo

!          if(uni_c.eq.3.and.quadrupel%ua2.eq.1) then
!           print*,'CSHGGb bc 1 1',sum(CSHGb(c,b,:,1)),sum(CSHGGb(c,b,:,1,1))
!           print*,'CSHGGb bc 1 1',sum(CSHGb(b,c,:,1)),sum(CSHGGb(b,c,:,1,1))
!          do ma=1,size(CSHGb,3)
!           print*,ma,CSHGb(c,b,ma,1),CSHGGb(c,b,ma,1,1),'c,b'
!           print*,ma,CSHGb(b,c,ma,1),CSHGGb(b,c,ma,1,1),'b,c'
!           print*,ma,sum(CSHGb(:,:,ma,1)),sum(CSHGGb(:,:,ma,1,1)),':,:'
!          enddo
!          endif
!           if(uni_c.eq.3.and.quadrupel%ua2.eq.1) print*,sum(CSHGb(:,:,:,1)),sum(CSHGGb(:,:,:,1,1)),'CSHGb-CSHGGb'

           call calculate_factor_GGa(1,b,c,a,la)    ! A
           call calculate_factor_GGa(2,c,a,b,lb)    ! B
           call calculate_factor_GGa(3,a,b,c,lc)    ! C

          ! similar statements for GGb counterpart follow

           call calculate_factor_GGb(1,b,c,a,la)    ! A
           call calculate_factor_GGb(2,c,a,b,lb)    ! B
           call calculate_factor_GGb(3,a,b,c,lc)    ! C
           ! depens on CSHGGb CSH CSHGb only

!        if(uni_c.eq.3.and.quadrupel%ua2.eq.1) then
!         do ma=1,size(A_factorGb,2)
!
!             print*,'1 A_factorGGb 1 1', sum(A_factorGb(1,ma,:,1)),sum(A_factorGGb(1,ma,:,1,1)),ma
!             print*,'2 A_factorGGb 1 1', sum(A_factorGb(2,ma,:,1)),sum(A_factorGGb(2,ma,:,1,1))
!             print*,'3 A_factorGGb 1 1', sum(A_factorGb(3,ma,:,1)),sum(A_factorGGb(3,ma,:,1,1))
!          do lm=0,size(A_factorGb,3)-1
!            print*,'1 A_factorGGb 1 1', A_factorGb(1,ma,lm,1),A_factorGGb(1,ma,lm,1,1),ma,lm
!            print*,'2 A_factorGGb 1 1', A_factorGb(2,ma,lm,1),A_factorGGb(2,ma,lm,1,1)
!            print*,'3 A_factorGGb 1 1', A_factorGb(3,ma,lm,1),A_factorGGb(3,ma,lm,1,1)
!          enddo
!         enddo
!        endif

          ! now GaGb and GbGa cross therms should follow

           call calculate_factor_GaGb(1,b,c,a,la)    ! A
           call calculate_factor_GaGb(2,c,a,b,lb)    ! B
           call calculate_factor_GaGb(3,a,b,c,lc)    ! C

         endif nolc_2dervs

       endif nolc_derivatives

           allocate(AngCMat(n_jj,-la:la, -lb:lb, -lc:lc), stat=alloc_stat(22))
           MEMLOG(size(AngCMat))
           if(alloc_stat(22).eq.0) then
            AngCMat = (0.0_r8_kind,0.0_r8_kind)
           else
           call error_handler("calc_3c AngCMat allocation failed")
           endif
           alloc_stat(22)=1

         if(integralpar_gradients) then
           !    now extended with GaGb
           call calc_angular_cmplx3cGa(la,lb,lc,n_jj,3,integralpar_dervs, &
                                       AngCMatGb=AngCMatGb) !!!! 1 lc.eq.0 in nucfit
         else
           call calc_angular_cmplx3cGa(la,lb,lc,n_jj,3,integralpar_dervs)    !!!! 2 lc.eq.0 in nucfit
         endif

           AngRMat=>AngR3cMat(:,:,:,:,j)
           AngRMat=0.0_r8_kind

           if(integralpar_3cob_grad) then
              AngRMatGa=>AngR3cMatGa(:,:,:,:,:,j)
              AngRMatGb=>AngR3cMatGb(:,:,:,:,:,j)

             if(integralpar_dervs) then
              AngRMatGGa=>AngR3cMatGGa(:,:,:,:,:,:,j)
              AngRMatGGb=>AngR3cMatGGb(:,:,:,:,:,:,j)
              AngRMatGaGb=>AngR3cMatGaGb(:,:,:,:,:,:,j)
             endif

           endif

         if(integralpar_gradients) then
           call Ang_Cmplx_to_RealGa(la,lb,lc,integralpar_dervs,AngCMatGb=AngCMatGb)

!          if(uni_c.eq.3.and.quadrupel%ua2.eq.1) then
!           do ma=-la,la
!           do mb=-lb,lb
!           print*,'AngCMatGb 1 1',sum(AngCMatGb(:,ma,mb,0,1)),sum(AngCMatGGb(:,ma,mb,0,1,1)),ma,mb
!           enddo
!           enddo
!           print*,'AngCMatGb 1 1',sum(AngCMatGb(:,:,:,0,1)),sum(AngCMatGGb(:,:,:,0,1,1))
!           print*,'AngR3cMatGGb 1 1',sum(AngRMatGb(:,:,:,1,1)),sum(AngRMatGGb(:,:,:,1,1,1))
!          endif

          else
           call Ang_Cmplx_to_RealGa(la,lb,lc,integralpar_dervs)
          endif

           MEMLOG(-size(AngCMat))
           deallocate(AngCMat, stat=alloc_stat(22))
           ASSERT(alloc_stat(22).eq.0)

           if(integralpar_3cob_grad) then
             if(integralpar_dervs) then
             MEMLOG(-size(AngCMatGGa)*6-size(AngCMatGa)*4)
              deallocate(AngCMatGa,AngCMatGb, &
                AngCMatGGa,AngCMatGGb,AngCMatGaGb,stat=alloc_stat(48))
             else
              MEMLOG(-size(AngCMatGa)*4)
              deallocate(AngCMatGa,AngCMatGb,stat=alloc_stat(48))
             endif
              ASSERT(alloc_stat(48).eq.0)
           endif

      radial_nuc: if(integralpar_2cob_nuc.or.integralpar_3cob_grad) then

      rad_nuc_dervs_setup: if(integralpar_3cob_grad) then

        ! SETUP RADIAL_NUC CALCULATIONS

        if(integralpar_dervs) then
           allocate(radialNuc_matGGa(num,n_jj),radialNuc_matGa(num,n_jj), &
                    jj_facGa(num,0:la+lb,llam_m:la+lb),stat=alloc_stat(78))
           MEMLOG(+size(radialNuc_matGGa)*2+size(jj_facGa))
        else
           allocate(radialNuc_matGa(num,n_jj), &
                    jj_facGa(num,0:la+lb,llam_m:la+lb),stat=alloc_stat(78))
           MEMLOG(+size(radialNuc_matGa)+size(jj_facGa))
        endif
           ASSERT(alloc_stat(78).eq.0)
           alloc_stat(78)=1
       endif rad_nuc_dervs_setup

          allocate(radialNuc_mat(num,n_jj),jj_fac(num,0:la+lb,llam_m:la+lb),stat=alloc_stat(23))
            MEMLOG(size(radialNuc_mat)+size(jj_fac))
            ASSERT(alloc_stat(23).eq.0)
            alloc_stat(23)=1
            alloc_stat(43)=1

        ! DONE SETUP RADIAL_NUC CALCULATIONS

            call calc_radial_nucP(la,lb)

            call dgemm('n','n',num,(2*La+1)*(2*Lb+1),n_jj, &
                       1.0_r8_kind, radialNuc_mat(:,:),num,&
                       AngRMat(:,:,:,:),n_jj, &
                        1.0_r8_kind,nuc(:,:,:),num)
                       !^    - contribs of equal atoms summed up

            if(integralpar_3cob_grad) then
               call calc_radial_nucGa(La,Lb) ! first dervs as yet
               if(integralpar_dervs)  call calc_radial_nucGGa(la,lb)
               call radial_ang_nuc(La,Lb)
               ! result: nuclear attraction derivatives
            endif

       close_rad_nuc_dervs: if(integralpar_3cob_grad) then
        if(integralpar_dervs) then
           MEMLOG(-size(radialNuc_matGGa)*2-size(jj_facGa))
           deallocate(radialNuc_matGGa,radialNuc_matGa,jj_facGa, &
                                               stat=alloc_stat(78))
        else
           MEMLOG(-size(radialNuc_matGa)-size(jj_facGa))
           deallocate(radialNuc_matGa,jj_facGa,stat=alloc_stat(78))
        endif
           ASSERT(alloc_stat(78).eq.0)
       endif close_rad_nuc_dervs

            MEMLOG(-size(radialNuc_mat)-size(jj_fac))
            deallocate(radialNuc_mat,jj_fac,stat=alloc_stat(23))
            ASSERT(alloc_stat(23).eq.0)
            alloc_stat(43)=0 !       jj_fac

       endif radial_nuc

       enddo nolc_equal_atoms !(j)


       if(l_time) call stop_timer(timer_ang3c)
       if(l_time) call start_timer(timer_rad3c)

         nolc_calcrad_alloc: if(integralpar_gradients) then

!      print*,'ca_dervs_mat out prim_int_coul_dervs', &
!       sum(ca_dervs_mat(:,:,:,:,1,3)),imc,sum(prim_int_coul_dervs(1,1)%m)

          allocate(grad_totsymM(num,(2*Lb+1),(2*La+1),grad_dim), &
                   gradM(num,(2*Lb+1),(2*La+1),6), &
                   stat=alloc_stat(66))
           ASSERT(alloc_stat(66).eq.0)
           MEMLOG(size(grad_totsymM)+size(gradM))
                  alloc_stat(66)=1

         endif nolc_calcrad_alloc

          sfit: if(new_3c_co.and.(integralpar_3c_co.or.integralpar_3cob_grad)) then

           ncexps = unique_atoms(i)%l_ch(0)%n_exponents
           cexps => unique_atoms(i)%l_ch(0)%exponents(:)

           allocate(coul_int(0)%m(num, ncexps, 1, 2*lb+1, 2*la+1), stat=alloc_stat(14)) ! 2
           MEMLOG(+size(coul_int(0)%m))
            ASSERT(alloc_stat(14).eq.0)
                   alloc_stat(14)=1

            pointer_coul => coul_int(0)%m
            pointer_coul = 0.0_r8_kind

          do k=1,ncexps

            c_exp_arg2(:)=(cexps(k)*sqrt(fact0(:)+cexps(k)))

            call calc_radial_3cSA(la,lb,lc,N_max) ! (2) lc.eq.0 case

            if(integralpar_3cob_grad) then
              call calc_radial_3cGa(la,lb,lc)
              if(integralpar_dervs)   call calc_radial_3cGGa(la,lb,lc) ! res:radial3cmatGG
            endif


            call Radial_Ang_3cS(radial3cmat,la,lb,lc) ! Radial_Ang_3cSA for lccase, first call
                 ! depends on AngR3cMatGGb

           enddo ! exps
         endif sfit

        ncexps = ua_pointer%r2_ch%n_exponents
        r2fit: if(ncexps.ne.0.and.new_3c_co.and. &
                 (integralpar_3c_co.or.integralpar_3cob_grad)) then

           allocate(coul_int(-1)%m(num, ncexps, 1, 2*lb+1, 2*la+1), stat=alloc_stat(14)) ! 3
           MEMLOG(+size(coul_int(-1)%m))
              ASSERT(alloc_stat(14).eq.0)
                     alloc_stat(14)=1
           pointer_coul => coul_int(-1)%m
           pointer_coul = 0.0_r8_kind


        allocate(g_shift_fac(num),radial3cmat_g(num,n_jj,n_equals),stat=alloc_stat(24))
                 ASSERT(alloc_stat(24).eq.0)
         MEMLOG(size(g_shift_fac)+size(radial3cmat_g))
                        alloc_stat(24)=1

        cexps => unique_atoms(i)%r2_ch%exponents(:)
        do k=1,ncexps
         c_exp_arg2(:)=(cexps(k)*sqrt(fact0(:)+cexps(k)))
         g_shift_fac(:)=(2*fact0(:)+3*cexps(k))/(two*(cexps(k)+fact0(:))*cexps(k))


            call calc_radial_r2(la,lb,0) ! now before dervs extension this procedure
                                         ! treats integrals  gradients & dervs

             call Radial_Ang_3cS(radial3cmat_g,la,lb,-1)             ! second call
                 ! the same procedure called for s and r2 case but
                 ! with different first parameter
!!!      commented to leave only nuclear contrib
        enddo

        MEMLOG(-size(g_shift_fac)-size(radial3cmat_g))
        deallocate(g_shift_fac,radial3cmat_g,stat=alloc_stat(24))
        ASSERT(alloc_stat(24).eq.0)
       endif r2fit

         ! this block is here while nuc grads are not treated
         ! in this route

         if(integralpar_gradients) then
         MEMLOG(-size(gradM)-size(grad_totsymM))
            deallocate(gradM,grad_totsymM,stat=alloc_stat(66))
             ASSERT(alloc_stat(66).eq.0)
          endif

             if(l_time)       call stop_timer(timer_rad3c)

        MEMLOG(-size(nucR_fac))
        deallocate(nucR_fac,stat=alloc_stat(25))
        ASSERT(alloc_stat(25).eq.0)

     endif lccase !/else

     !--------------------------------------------------------------------------
     ! now all primitive calculated and one needs deallocate working hands

       if(integralpar_3cob_grad) then
        if(integralpar_dervs) then
                MEMLOG(-size(AngR3cMatGa)*2-size(AngR3cMatGGa)*3-size(radial3cmatG)*2)
            deallocate(radial3cmatG,radial3cmatGG, AngR3cMatGa,AngR3cMatGb, &
                   AngR3cMatGGa,AngR3cMatGGb,AngR3cMatGaGb, stat=alloc_stat(50))
        else
                MEMLOG(-size(AngR3cMatGa)*2-size(radial3cmatG))
            deallocate(radial3cmatG,AngR3cMatGa,AngR3cMatGb, stat=alloc_stat(50))
        endif
            ASSERT(alloc_stat(50).eq.0)
                   alloc_stat(49)=0          ! AngRMatGa
       endif

       if(integralpar_dervs) then
     MEMLOG(-size(gamma_help)-size(gamma_help_co)-size(exp_arg)-size(gamma_arg2))
     MEMLOG(-size(jj_fac_SA)-size(gamma_help_g)*3-size(jj_facG)*2)
     MEMLOG(-size(radial3cmat)-size(AngR3cMat)-size(Gt))
       MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw)-size(mup_llam)-size(i3_func_fac))
       MEMLOG(-size(i3_fac))
       deallocate(AngR3cMat,radial3cmat, &
                  gamma_help,gamma_help_co, gamma_arg2,exp_arg, &
                  jj_fac_SA, &
                  gamma_help_g, &
                  gamma_help_gG,gamma_help_gGG, &           ! r2_radial 1st and 2nd dervs
                  jj_facG, jj_facGG, &                      ! r2_radial 1st and 2nd dervs
                  i3_fac,i3_func_fac,mup_llam, &
                  alph_pw,beta_pw,fact2_pw, &
                  Gt,  stat=alloc_stat(20))

       elseif(integralpar_gradients) then
     MEMLOG(-size(gamma_help)-size(gamma_help_co)-size(exp_arg)-size(gamma_arg2))
     MEMLOG(-size(jj_fac_SA)-size(gamma_help_g)*2-size(jj_facG))
     MEMLOG(-size(radial3cmat)-size(AngR3cMat)-size(Gt))
       MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw)-size(mup_llam)-size(i3_func_fac))
       MEMLOG(-size(i3_fac))
       deallocate(AngR3cMat,radial3cmat, &
                  gamma_help,gamma_help_co, gamma_arg2,exp_arg, &
                  jj_fac_SA, &
                  gamma_help_g,gamma_help_gG,jj_facG, &
                  i3_fac,i3_func_fac,mup_llam, &
                  alph_pw,beta_pw,fact2_pw,Gt,  stat=alloc_stat(20))

       else ! i.e. no gradients and dervs
       MEMLOG(-size(alph_pw)-size(beta_pw)-size(fact2_pw)-size(mup_llam)-size(i3_func_fac))
       MEMLOG(-size(i3_fac))
       MEMLOG(-size(gamma_help)-size(gamma_help_co)-size(exp_arg)-size(gamma_arg2))
       MEMLOG(-size(jj_fac_SA)-size(gamma_help_g))
       MEMLOG(-size(radial3cmat)-size(AngR3cMat))
       deallocate(AngR3cMat,radial3cmat, &
                  gamma_help,gamma_help_co, gamma_arg2,exp_arg, &
                  jj_fac_SA,gamma_help_g, &
                  i3_fac, &
                  i3_func_fac,mup_llam, &
                  alph_pw,beta_pw,fact2_pw,  stat=alloc_stat(20))
       endif

       ASSERT(alloc_stat(20).eq.0)
       alloc_stat(19)=0 ! radial3cmat
       alloc_stat(18)=0 ! gamma_help,gamma_help_co, gamma_arg2,exp_arg ,Gt
       alloc_stat(17)=0 ! alph_pw,beta_pw,fact2_pw i3_fac,i3_func_fac,mup_llam
      enddo lc_fit

      MEMLOG(-size(radang_temp)-size(c_exp_arg2))
      deallocate(radang_temp,c_exp_arg2,stat=alloc_stat(34))
      ASSERT(alloc_stat(34).eq.0)

      if(integralpar_3cob_grad) then
       if(integralpar_dervs) then
          MEMLOG(-size(dervs_totsymM)-size(dervs_totsymM_temp))
          MEMLOG(-size(dervsM)-size(ca_dervsM))
          MEMLOG(-size(radangF))
        deallocate( radangF,dervsM,dervs_totsymM,dervs_totsymM_temp,ca_dervsM, &
                   stat=alloc_stat(55))
       else
        MEMLOG(-size(radangF))
        deallocate(radangF,stat=alloc_stat(55))
       endif
        ASSERT(alloc_stat(55).eq.0)
      endif
      !---------------------------------------------------

      call start_timer(timer_prod_nested)

     call calc_3c_fitcontract(i,equalb)


      gr_fitcont_clos: if(integralpar_3cob_grad) then
        if(moving_c) then
         MEMLOG(-size(rotmat_eq)-size(do_rotation_eq))
         deallocate (rotmat_eq, do_rotation_eq, stat=alloc_stat(62))
         ASSERT(alloc_stat(62).eq.0)
        endif
        MEMLOG(-size(gamma_arg_xc))
        deallocate(gamma_arg_xc,stat=alloc_stat(63))
          ASSERT(alloc_stat(63).eq.0)
      endif  gr_fitcont_clos

       call stop_timer(timer_prod_nested)

          MEMLOG(-size(gamma_arg2_vec)-size(temp_overlap))
          MEMLOG(-size(gamma_help_r2))
       deallocate(temp_overlap, gamma_arg2_vec, gamma_help_r2,stat=alloc_stat(16))
       ASSERT(alloc_stat(16).eq.0)
       alloc_stat(15)=0 ! gamma_arg2_vec,gamma_help_r2

     endif lcfit_calc

    end do unique_atoms_fit_nuc


!      deallocate(temp_overlap_nuc,gamma_argXC, stat=alloc_stat(12))
!      ASSERT(alloc_stat(12).eq.0)

       grad_3cob: if(integralpar_3cob_grad) then

                              pointer_prim_int=>prim_int_3cob_grad
         call add_nuc_or_pvsp(equalb,grad_mat)
                          ! ** here a & b coul grads are summed up
                          !    to pointer_prim_int

!         if(integralpar_dervs) then
!          call add_dervs(equalb,dervs_mat,ima,imb,prim_int_coul_dervs)  !(3) - third of 3 contribs
!         endif


         spmda: if (model_density .and. spin_polarized) then
          spin_index = gradient_data_n_gradients + 1
          pointer_prim_int => prim_int_3cob_grad(spin_index:)
          call add_nuc_or_pvsp(equalb,grad_mat_spin)
         endif spmda

         cpks1: if(integralpar_cpksdervs) then
#ifndef no_cpks_coul_grads
          pointer_prim_3c=>prim_int_cpks_coul_grad !(1)
#endif
          call add_cpks_coul(equalb,cpks_grad_mat)
          MEMLOG(-size(cpks_grad_mat))
          deallocate( cpks_grad_mat, stat=cpksalloc(20))
          ASSERT(cpksalloc(20).eq.0)
          cpksalloc(20)=1
         endif cpks1

        endif grad_3cob

      ! NUCL: NUCLEAR ATTRACTION HERE:
      if(new_nuc.and.integralpar_2cob_nuc) then
       do mb=1,2*lb+1
          do ma=1,2*la+1
           prim_int_2cob_nuc(:,:,mb,ma)=&!prim_int_2cob_nuc(:,:,mb,ma)+ &
                                         unpack(nuc(:,ma,mb),cutoff,zero)
          enddo
       end do
      endif

    endif nucfit_calc

    DPRINT 'DONE CALCULATE NUCLEAR ATRACTION AND COULOMB INTEGRALS'

#ifdef WITH_EPE
   if(new_ewpc.and. &
      (integralpar_2cob_nuc.or.integralpar_3cob_grad).and.ewpc_n.ne.0) then
      print*, 'calc_3c_ewpc call'
!    if((integralpar_3cob_grad.and.integralpar_dervs).and.ewpc_n.ne.0) then
      call calc_3c_ewpc(la,lb,equalb)
   end if
#endif

         if(integralpar_3cob_grad.and.integralpar_dervs) then
          call add_dervs(equalb,dervs_mat,ima,imb,prim_int_coul_dervs)  !(3) - third of 3 contribs
         endif

   if (integralpar_2cob_nuc.or.integralpar_3cob_grad) then ! nuc probable not used in add_ewpc sec hense this should
                                   ! stay before add_ewpc
       MEMLOG(-size(nuc))
       deallocate(nuc, stat=alloc_stat(6))
       ASSERT(alloc_stat(6).eq.0)
    endif
    if(integralpar_gradients) then
    MEMLOG(-size(radang_temp_ga)*2-size(Ga)*4)
    deallocate(radang_temp_ga,radang_temp_gb,Ga,Gb,GtFvec, &
               stat=alloc_stat(77))
      ASSERT(alloc_stat(77).eq.0)
    if(integralpar_dervs) then
      MEMLOG(-size(radang_temp_gga)*3-size(GGt)*7-size(GGa)*3)
      deallocate(radang_temp_gga,radang_temp_ggb,radang_temp_gagb, &
                 GGt,GtGa,GtGb,GGa,GGb,GaGb, &
                 stat=alloc_stat(76))
      ASSERT(alloc_stat(76).eq.0)
    endif
    endif

    !If you want to activate this subroutine again do not forget
    !to comment corresponding parts in calc_prim_shgi (integral_calc_quad_2cob3c)
!!$   if(integralpar_2cob_potential.or.integralpar_solv_grad) &
!!$    call calc_3c_solv(la,lb,equalb)

      MEMLOG(-size(temp_overlap_nuc)*4)
      deallocate(temp_overlap_nuc,gamma_argXC, stat=alloc_stat(12))
      ASSERT(alloc_stat(12).eq.0)

        MEMLOG(-size(bin_fac)-size(gamma_help_fac)*4)
        deallocate(bin_fac,gamma_help_fac, r2_jjj_fac,gamma_help_fac_g, &
                   ghelp_shift_fac,   stat=alloc_stat(11))
        ASSERT(alloc_stat(11).eq.0)

          do j1=0,l_max
             do j2=0,l_max
                do j3=0,l_max
                   if( jj_coeff(j1,j2,j3)%Num_M_terms > 0 ) then
                      MEMLOG(-size(jj_coeff(j1,j2,j3)%value)*4)
                      deallocate(jj_coeff(j1,j2,j3)%m1, jj_coeff(j1,j2,j3)%m2, &
                                 jj_coeff(j1,j2,j3)%m3, jj_coeff(j1,j2,j3)%value, stat=alloc_stat(31))
            ASSERT(alloc_stat(31).eq.0)
                   end if
                end do
             end do
          end do

          MEMLOG(-size(jj_coeff)-size(valjj)*4)
          deallocate(jj_coeff,m1jj,m2jj,m3jj,valjj, stat=alloc_stat(31))
          ASSERT(alloc_stat(31).eq.0)

       if(integralpar_gradients) then
        MEMLOG(-size(two_aexp_arr)*5)
        deallocate(two_fact2_vec,two_aexp_arr,two_bexp_arr,stat=alloc_stat(61))
         ASSERT(alloc_stat(61).eq.0)
        endif
    MEMLOG(-size(cutoff)-size(fact0)*6)
    deallocate(gamma_arg,fact0,fact1,fact2,cutoff, stat=alloc_stat(4))
    ASSERT(alloc_stat(4).eq.0)
    alloc_stat(3)=0 ! cutoff
    call csh_lm_map(l_max,"close")

    MEMLOG(-size(csh_lm_of))
    deallocate(csh_lm_of,stat=alloc_stat(8))
    ASSERT(alloc_stat(8).eq.0)
    MEMLOG(-size(afac_sig)-size(even_triangle))
    deallocate(afac_sig,even_triangle,stat=alloc_stat(7))
    call shutdown_triplets

    MEMLOG(-size(aexp_arr)*2)
    deallocate(aexp_arr,bexp_arr,stat=alloc_stat(5))
    ASSERT(alloc_stat(5).eq.0)

    MEMLOG(-size(CSH)-size(RSH)-size(A_factor))
    deallocate(CSH, RSH,A_factor, stat=alloc_stat(9))
    ASSERT(alloc_stat(9).eq.0)

    if(integralpar_gradients) then

        MEMLOG(-size(CSHG)*2-size(rshg)-size(A_factorGa)*2)
        deallocate(CSHG,CSHGb,rshg, A_factorGa,A_factorGb, stat=alloc_stat(35))
        ASSERT(alloc_stat(35).eq.0)

   if(integralpar_dervs.or.ewpcdervs) then

      MEMLOG(-size(CSHGGa)*3-size(rshgg)-size(A_factorGGa)*3)
      deallocate(CSHGGa,CSHGGb,CSHGaGb, rshgg, &
        A_factorGGa,A_factorGGb,A_factorGaGb,stat=alloc_stat(65))
       ASSERT(alloc_stat(65).eq.0)

       MEMLOG(-size(dervs_mat))
       deallocate(dervs_mat,stat=alloc_stat(72))
       ASSERT(alloc_stat(72).eq.0)
    endif

       MEMLOG(-size(grad_mat))
       deallocate(grad_mat,stat=alloc_stat(46))
       ASSERT(alloc_stat(46).eq.0)
    endif

    if (integralpar_3cob_grad) then
     if (model_density .and. spin_polarized) then
        MEMLOG(-size(grad_mat_spin))
        deallocate( grad_mat_spin, stat=alloc_stat(64))
        ASSERT(alloc_stat(64).eq.0)
     endif
    endif

#ifdef no_cpks_coul_grads
   if(.not.integralpar_pot_for_secderiv.and..not.integralpar_solv_grad) then
   if(allocated(cpks))  call calc_cpks_grad_totalsym()
   endif
#endif

     do i=1,size(alloc_stat)
      if(alloc_stat(i).ne.0) print*,i,' alloc_stat ne 0'
     enddo

  ! this logic is hopefully present in integralpar_setup():
! integralpar_dervs=integralpar_gradients &
!                   .and.integralpar_2dervs.and..not.integralpar_cpks_contribs

  if(integralpar_gradients) then
   MEMSET(0)
  endif
  contains

    subroutine calc_cpks_grad_totalsym()
    use gradient_data_module, only: cpks_gradient_totalsym, dervs_totalsym
    implicit none
    ! *** end of interface ***

    integer(i4_kind):: i_grad, k2dr

    if(integralpar_solv_grad.or.integralpar_pot_for_secderiv) return
    if(integralpar_cpks_contribs) then
!          cpks_gradient_totalsym=0.0_r8_kind ! no charge overlap contribs
          ! this cpks_gradient_totalsym is not used in cpks calcs
          ! but contribute to dervs_totalsym via product with cpks_fitcoeff_grads here
          do i_grad=1,size(cpks_gradient_totalsym,2)
           do k2dr=1,size(cpks_gradient_totalsym,2)
           dervs_totalsym(i_grad,k2dr)=dervs_totalsym(i_grad,k2dr) +  &
           dot_product(cpks_gradient_totalsym(:,i_grad),cpks_fitcoeff_grads(:,k2dr))
           enddo
          enddo
    endif

    end subroutine calc_cpks_grad_totalsym

  end subroutine Calc_3center



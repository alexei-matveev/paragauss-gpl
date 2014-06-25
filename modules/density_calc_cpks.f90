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
!================================================================
! Public interface of module
!================================================================
module density_calc_cpks
!---------------------------------------------------------------
!
!  Purpose: split off the cpksdervs_xc(..) sub from density_calc_module
!
!  Date: 1/07
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
! define FPP_TIMERS 2
# include "def.h"
!## define FPP_NOBLAS 0
  use type_module  ! type specification parameters
  use orbitalstore_module
! use orbital_module
! use iounitadmin_module
! use occupied_levels_module
  use datatype, only: arrmat4,arrmat5,arrmat6,arrmat7
! use options_module, only: options_spin_orbit
  use unique_atom_module
! use comm_module
  USE_MEMLOG
  implicit none
  private
  save

#ifdef WITH_EXPERIMENTAL
# undef WITH_SECDER
#endif

#ifdef WITH_SECDER
  interface transform_to_cart
   module procedure transform_to_cart_grarho
   module procedure transform_to_cart_gragamma
   module procedure transform_to_cart_dervs
   module procedure transform_to_cart_dervs_grarho
   module procedure transform_to_cart_dervs_grarho2
  end interface
  
!== Interrupt end of public interface of module =================

  integer(kind=i4_kind) :: alloc_stat(34)=1
! integer(kind=i4_kind) :: ispin,n_irrep,max_equal
! integer(kind=i4_kind),allocatable :: dims(:),partners(:),&
!      & dims_core(:),dims_fit(:)
! real(kind=r8_kind), allocatable :: coeff_core(:) !AM???

! integer(i4_kind), parameter, private ::&
!      & VALUE=0, X=1, Y=2, Z=3,&
!      & UP=1, DN=2,&
!      & UPUP=1, DNDN=2, UPDN=3,&
!      & RE=1, IM=2

  real(r8_kind), parameter ::&
       & ZERO = 0.0_r8_kind,&
       & ONE  = 1.0_r8_kind,&
       & TWO  = 2.0_r8_kind

  !------------ public functions and subroutines ------------------

  public :: cpksdervs_xc

!================================================================
! End of public interface of module
!================================================================

  !------------ public types of subroutines formal parameters ---

contains 

  subroutine cpksdervs_xc(full_vl, vec_length, ispin, dvdrho,orbs_ob,phi_ob,grdwts, & 
                               nuc_grads,dfdrho,nuc_grarho,graw,i_m,rho, &
                               nuc_grarho_imp, &
                               nuc_gragamma_imp, &
                               dervsrho_imp, &
                               nuc_sec_der, &
                               df_drhodgamma, &
                               nuc_sec_derrho, grarho,graphi, &
                               df_dgammadgamma,dfdgrarho, &
                               dervs_grarho_imp)

    use eigen_data_module
    use occupation_module
    use machineparameters_module
    use cpksdervs_matrices
    use calc3c_switches
    use density_calc_module, only: n_irrep, partners ! FIXME: historical global vars
    implicit none

    real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))

    integer(kind=i4_kind),intent(in) :: vec_length, full_vl ! vec_length <= full_vl
    integer(kind=i4_kind),intent(in) :: ispin
    real(kind=r8_kind), optional,intent(inout)   :: dfdgrarho(full_vl, 2 * ispin - 1) ! 1 or 3
    real(kind=r8_kind), optional,intent(inout)   :: df_drhodgamma(full_vl, 5 * ispin - 4) ! 1 or 6
    real(kind=r8_kind), optional,intent(inout)   :: df_dgammadgamma(full_vl, 5 * ispin - 4) ! 1 or 6
    real(kind=r8_kind), optional,intent(in)      :: grarho(full_vl, 3, ispin)
    type(arrmat5),optional,intent(inout) :: nuc_sec_derrho(:)
    type(arrmat4),optional,intent(in) :: graphi(:,:)

    real(kind=r8_kind),intent(in)   :: dvdrho(full_vl, 2 * ispin - 1) ! 1 or 3
    real(kind=r8_kind),intent(in)   :: dfdrho(full_vl, ispin) ! 1 or 2

    ! passed array may be longer, we only care about this much:
    real(kind=r8_kind),intent(in) :: grdwts(vec_length)
    type(orbital_type),pointer :: orbs_ob(:)
    type(orbital_spin_type),pointer :: phi_ob(:)
    type(orbital_nuclear_gradient_type), pointer :: nuc_grads(:,:)
    type(orbital_nuclear_sec_der_type),pointer, optional :: nuc_sec_der(:,:)
    type(arrmat4),allocatable      :: nuc_grarho(:) ! gradient of rho

    type(arrmat5),intent(inout),optional  :: nuc_gragamma_imp(:) 
    type(arrmat4),intent(inout),optional  :: nuc_grarho_imp(:) 
    type(arrmat5),intent(inout),optional  :: dervsrho_imp(:,:) 
    type(arrmat5),intent(inout),optional  :: dervs_grarho_imp(:,:) 
    real(kind=r8_kind),optional, intent(in)   :: rho(:,:) ! FIXME: UNUSED?
    type(arrmat3),allocatable :: graw(:)
    integer(kind=i4_kind),intent(in) ::i_m
    !** End of interface *****************************************

    real(kind=r8_kind),allocatable::nucgra_afac(:,:),nuc2der_afac(:,:,:), &
                    nucgra_efac(:,:,:),help_sd_eq(:,:,:)!!!,help_nuc_eq(:,:,:)
    real(kind=r8_kind),allocatable::ghelp_nuc_eq(:,:,:,:)

    integer(kind=i4_kind):: i_grad
    real(kind=r8_kind)  :: help_gp(3) 
    real(kind=r8_kind), allocatable :: PhiB(:,:),PhiSkl(:,:),graPhiB(:,:,:),graPhiSkl(:,:,:)
#if 0
    real(kind=r8_kind), pointer :: help_nuc_B(:,:,:,:), &
                                       help_2nd_Skl(:,:,:,:),help_2nd_B(:,:,:,:), &
                                       b_help(:,:,:)
#endif

    real(kind=r8_kind), allocatable :: nuc_imp_gragamma(:,:,:)

    real(kind=r8_kind):: nuc_imp_grarho(vec_length,size(cpks,1))

    real(kind=r8_kind),allocatable:: help_sd(:,:,:)
    real(kind=r8_kind),dimension(vec_length):: dvdrho_fac,wdvdrho,wdfdrho,graphi_fac
    real(kind=r8_kind),dimension(vec_length,3):: phi_p_fac

    real(kind=r8_kind),  allocatable :: drhodgamma_fac(:),phi_adfac(:), &
                         phi_p_adfac(:),graphi_adfac(:,:), &
                         dgammadgamma_fac(:,:),graphi_dgdg_fac(:),wdfdgrarho(:), &
                         dgrarho_fac(:,:,:),graphi_dgrarho_fac(:,:), &
                         phi_drdg_fac(:,:),grarho_drhodgamma(:,:),graw_adfac(:,:), &
                         phi_grarho_dgrarho(:,:),grarho_dgrarho(:,:),graphi_grarho_dgrarho(:), &
                         phi_wdgrarho(:,:),grarho_grarho(:,:,:)
    target:: grarho_drhodgamma,phi_drdg_fac
     real(kind=r8_kind), pointer:: grarho_dgammadrho(:,:),phi_dgdr_fac(:,:)

    type(arrmat4), allocatable ::imp_dervsrho(:)

    type(arrmat4), allocatable ::imp_dervs_grarho(:)
 
    real(kind=r8_kind),pointer :: ob(:,:,:),eigv(:,:,:),phi_p(:,:,:)
    integer(i4_kind) :: i,l,eig_dim
    integer(i4_kind) :: k,i_occ,j_occ,j_vir,i_ea,i_ua,counter,orb_dim,occ_dim, &
                        n_equal_max
    integer(i4_kind) :: kk
    integer(i4_kind) :: i_ua2
    real(kind=r8_kind),allocatable  :: help_nuc_arr(:,:,:,:),nuc_help(:,:,:,:), &
                                       sec_der_help(:,:,:,:)
#if 0
    real(kind=r8_kind),allocatable  :: help_nuc_arr_im(:,:,:,:)
#endif
    real(kind=r8_kind),allocatable  :: im_help_nuc_arr(:,:,:),vim_help_nuc_arr(:,:,:)
    real(kind=r8_kind),allocatable  :: im_help_sec_der_arr(:,:,:),vim_help_sec_der_arr(:,:,:)
    real(kind=r8_kind),allocatable,target  :: help_sec_der_arr(:,:,:,:)
    logical:: do_rotation
    logical:: nl_calc,Q_calc
    integer(i4_kind) :: index,grad_dim,i_irr,spin,n_ai,vir_dim
    real(kind=r8_kind)::grarho_dervs(3),rho_dervs(3)
    real(kind=r8_kind),pointer:: help_2nd(:,:,:,:)
    integer(i4_kind),parameter :: updn=3,uord=1,ud=2  
    integer(kind=i4_kind), parameter:: uup=1,ddn=2,udn=3,dup=4,udup=5,uddn=6
    integer(kind=i4_kind), parameter:: aa=1,bb=2,ab=3,ac=4,bc=5,cc=6

!    MEMSET(0) ! cpksdervs_xc
    ASSERT(vec_length<=full_vl)

    Q_calc=.not.cpks_Qxc_calculated
    nl_calc=present(df_drhodgamma)

    call allocate_intermediates()
    
    n_equal_max=maxval(unique_atoms(:)%n_equal_atoms) ! maximum number of equal_atoms

    FPP_TIMER_START(t_init_mocalc)
    do spin=1,ispin
     do i_irr=1,n_irrep

       !!!  in particular nuc_dervsrho_imp and dervs_grarho_imp  
       !!!  s1 contribs  are calculated

       i=i_irr
       ob=> orbs_ob(i_irr)%o
       eigv=>eigvec(i_irr)%m
       eig_dim=size(eigvec(i_irr)%m,1) ! sic!
       do l=1,partners(i_irr)
             call dgemm( 'n','n', vec_length,eig_dim,eig_dim,1.0_r8_kind,&
                         ob(1,1,l),full_vl, eigv(1,1,spin),eig_dim,&
                         0.0_r8_kind, phi_ob(i_irr)%o(:,:,l,spin),full_vl)
       enddo
     enddo
    enddo
   FPP_TIMER_STOP(t_init_mocalc)
   

  call setup_imp_dervs()

   global_spin_loop: do spin=1,ispin

   if(nl_calc)  nuc_imp_gragamma=0.0_r8_kind ! both for q and b calcs
                nuc_imp_grarho=0.0_r8_kind

   wdfdrho(:vec_length)=dfdrho(:vec_length,spin)*grdwts(:vec_length)
   wdvdrho(:vec_length)=dvdrho(:vec_length,spin)*grdwts(:vec_length)
   if(nl_calc)  then
    wdfdgrarho(:vec_length)=dfdgrarho(:vec_length,spin)*grdwts(:vec_length)
   endif


   ! s1 and B contributions to nuc_imp_grarho are calculated with
   ! separate calls and B contrib is real true for B(n_cpks+1)

    if(Q_calc) then

    init_irr: do i_irr=1,n_irrep

       !!!  in particular nuc_dervsrho_imp and dervs_grarho_imp  
       !!!  s1 contribs  are calculated

       i=i_irr
       ob=> orbs_ob(i_irr)%o
       eigv=>eigvec(i_irr)%m
       eig_dim=size(eigvec(i_irr)%m,1) ! sic!
        occ_dim=size(cpks(1,i_irr,spin)%s1,1)
        vir_dim=size(cpks(1,i_irr,spin)%hbh,2)
        n_ai=occ_dim*vir_dim

       occm: if (occ_dim > 0) then ! occ. orbitals

      
!             ** (1) s1 contrib to density  impicite nuclear gradients
!             ** (2) s1 contrib to density  impicite nuclear derivatives  

     phi_p=>phi_ob(i_irr)%o(:,:,:,spin)     ! to calc occ. orbitals for this spin

   partner:do l=1,partners(i_irr)


     FPP_TIMER_START(t_s1_imp_grarho)
      call s1_impgrads(nuc_imp_grarho) !!! result also nuc_imp_gragamma
     FPP_TIMER_STOP(t_s1_imp_grarho)

          ! ** s1 contrib to density  impicit nuclear derivatives
          ! ** requires calculation of nuc grads for occ MO
          !    and the same is required for calculations of Qai and h1
          !    thus dervs_s1 can be calculated together with Qai and h1
          !    using the same nuc grads for occ MO

#if 1
     FPP_TIMER_START(t_s1_imp_dervsrho)
     if(present(dervsrho_imp)) call s1_dervs()
     FPP_TIMER_STOP(t_s1_imp_dervsrho)
#endif 

 enddo partner

 endif occm

 enddo init_irr

      call  transform_to_cart(nuc_grarho_imp,nuc_imp_grarho,spin) !!! s1 contrib
                                   !-------------- cart 1st dervs for given spin
                                   !               used to calculate 2nd dervs_xc
                                   !            
                                   !               nuc_imp_grarho  - intermadiate 
                                   !               to calc spin dependent Qai and h1

      if(nl_calc) &
      call  transform_to_cart(nuc_gragamma_imp,nuc_imp_gragamma,spin) !!! s1 contrib

 endif

    QorB: if(q_calc) then
    ! *** calculate Qai, h1 and cpks_grad_xc 

    Q_irr: do i_irr=1,n_irrep

       i=i_irr
       ob=> orbs_ob(i_irr)%o
       eigv=>eigvec(i_irr)%m
       eig_dim=size(eigvec(i_irr)%m,1) ! sic!
        occ_dim=size(cpks(1,i_irr,spin)%Qai,1)
        vir_dim=size(cpks(1,i_irr,spin)%Qai,2)
        n_ai=occ_dim*vir_dim

           ! ** calculate  Q structure  s1 contrib 
           ! ** calculate h1 contibs

 if(occ_dim.gt.0) then

  partnerQ: do l=1,partners(i_irr)
    grads_s1: do i_grad=1,size(cpks,1)   

#if 1 /* here v_lambda_imp contrib due to s1 nuc_imp_grarho */ 
       call v_lambda_imp(cpks(i_grad,i_irr,spin)%Qai,dvdrho(:,spin),cpks_spin=spin) 
       if(ispin.eq.2) &               !!! s1 contrib !!!
       call v_lambda_imp(cpks(i_grad,i_irr,3-spin)%Qai,dvdrho(:,updn),cpks_spin=3-spin)
       !!! result also h1 complite hamiltonian derivative
#endif
    enddo grads_s1
  enddo partnerQ


     FPP_TIMER_START(t_xc_qh_expl)
       if(nl_calc) then  
        graphi_dgrarho_fac(:,uord)=two/partners(i_irr)*dfdgrarho(:vec_length,spin)*grdwts(:vec_length)          !!!(5)
        if(ispin.eq.2) &
        graphi_dgrarho_fac(:,ud)=grdwts(:vec_length)*dfdgrarho(:vec_length,updn)/partners(i_irr)
        
        !!! grarho_dgrarho factor to be used in all_expl_dervs
        grarho_dgrarho(:vec_length,:3)=two/partners(i_irr)*grarho(:vec_length,:3,spin)* &                       !!! (7)
                                       spread(dfdgrarho(:vec_length,spin),2,3)
        if(ispin.eq.2) grarho_dgrarho=grarho_dgrarho+grarho(:vec_length,:3,3-spin)*&
                       spread(dfdgrarho(:vec_length,updn),2,3)/partners(i_irr)
       endif

                  phi_p=>phi_ob(i_irr)%o(:,:,:,spin)
     q_partners2: do l=1,partners(i_irr)

     FPP_TIMER_START(t_h1q_dvdrho)

  FPP_TIMER_STOP(t_h1q_dvdrho)

    call all_expl_dervs()
  enddo q_partners2

        FPP_TIMER_STOP(t_xc_qh_expl)
endif
   enddo Q_irr


  else QorB !----------------------------------

    FPP_TIMER_START(t_xc_bimp_grarho)
    B_imp_grarho_irr: do i_irr=1,n_irrep
       i=i_irr
       ob=> orbs_ob(i_irr)%o
       eigv=>eigvec(i_irr)%m
       eig_dim=size(eigvec(i_irr)%m,1) ! sic!
       phi_p=>phi_ob(i_irr)%o(:,:,:,spin)     ! occ. orbitals for this spin
       occ_dim=size(cpks(1,i_irr,spin)%HBH,1)
       vir_dim=size(cpks(1,i_irr,spin)%HBH,2)
       phi_p=>phi_ob(i_irr)%o(:,:,:,spin)     ! to calc occ. orbitals for this spin

       n_ai=occ_dim*vir_dim

      if(n_ai.gt.0)  call impgrads(nuc_imp_grarho) !!! result also nuc_imp_gragamma

    enddo B_imp_grarho_irr

      if(i_cpks.gt.n_cpks) then                     ! 2nd of 2 nuc_grarho_imp contribs
       ! nuc_grarho_imp is requird only for calculation of xc_dervs when 
       ! cpks equations are solved
        call  transform_to_cart(nuc_grarho_imp,nuc_imp_grarho,spin)
        if(nl_calc) &
         call  transform_to_cart(nuc_gragamma_imp,nuc_imp_gragamma,spin)
       endif

    FPP_TIMER_STOP(t_xc_bimp_grarho)


    B_irr: do i_irr=1,n_irrep
    !!! B dependent nuc_imp_dervs  and 
    !!! B nuc_imp_grarho contribs to AB and H1
       i=i_irr
       ob=> orbs_ob(i_irr)%o
       eigv=>eigvec(i_irr)%m
       eig_dim=size(eigvec(i_irr)%m,1)
       phi_p=>phi_ob(i_irr)%o(:,:,:,spin)     ! occ. orbitals for this spin
       occ_dim=size(cpks(1,i_irr,spin)%HBH,1)
       vir_dim=size(cpks(1,i_irr,spin)%HBH,2)
       n_ai=occ_dim*vir_dim

   vir_GxBai: if(n_ai.gt.0) then

#if 0 /* if this PhiB used? */
       allocate( PhiB(vec_length,vir_dim),stat=cpksalloc(114))
           ASSERT(cpksalloc(114).eq.0)
           MEMLOG(size(PhiB))
#endif

    dervsrho_vir: if(i_cpks.gt.n_cpks) then 

    FPP_TIMER_START(t_dervs_imp)

         allocate(im_help_nuc_arr(vec_length,3,occ_dim),vim_help_nuc_arr(vec_length,3,vir_dim), &
                  stat=cpksalloc(155))
         ! vim_help_nuc_arr:
         ASSERT(cpksalloc(155).eq.0)
         alloc_stat(31)=0 ! im_help_nuc_arr
         MEMLOG(size(im_help_nuc_arr)+size(vim_help_nuc_arr))
         if(nl_calc) then
           allocate(im_help_sec_der_arr(vec_length,6,occ_dim),vim_help_sec_der_arr(vec_length,6,vir_dim), &
                    stat=cpksalloc(156))
         ASSERT(cpksalloc(156).eq.0)
         MEMLOG(size(im_help_sec_der_arr)+size(vim_help_sec_der_arr))
         endif

       dervsrho_vacmo: if (occ_dim < eig_dim) then 

       call all_imp_dervs() ! n_cpks+1

       MEMLOG(-size(im_help_nuc_arr)-size(vim_help_nuc_arr))
       deallocate(im_help_nuc_arr,vim_help_nuc_arr,stat=cpksalloc(155))
       ASSERT(cpksalloc(155).eq.0)
       alloc_stat(31)=1 ! im_help_nuc_arr
       cpksalloc(155)=1 ! vim_help_nuc_arr
       if(nl_calc) then
           MEMLOG(-size(im_help_sec_der_arr)-size(vim_help_sec_der_arr))
           deallocate(im_help_sec_der_arr,vim_help_sec_der_arr,stat=cpksalloc(156))
           ASSERT(cpksalloc(156).eq.0)
           cpksalloc(156)=1
       endif

    endif dervsrho_vacmo

    FPP_TIMER_STOP(t_dervs_imp)
    endif dervsrho_vir
      

    FPP_TIMER_START(t_xc_ab)
    ! input:  nuc_imp_grarho ! output: ABi
      do l=1,partners(i_irr)
       do i_grad=1,size(cpks,1)
#if 1 /* here v_lambda_imp contrib due to B nuc_imp_grarho */
       call v_lambda_imp(cpks(i_grad,i_irr,spin)%ABi,dvdrho(:,spin),cpks_spin=spin)    
       if(ispin.eq.2) &
       call v_lambda_imp(cpks(i_grad,i_irr,3-spin)%ABi,dvdrho(:,updn),cpks_spin=3-spin)
       !!! result also contrib to h1 complite hamiltonian derivatives for converged solutions
#endif
      enddo 
     enddo 
    FPP_TIMER_STOP(t_xc_ab)
   

#if 0 /* if this PhiB used? */
       MEMLOG(-size(PhiB))
       deallocate(PhiB,stat=cpksalloc(114))
       ASSERT(cpksalloc(114).eq.0)
       cpksalloc(114)=1
#endif

 endif vir_GxBai

 enddo B_irr

 endif QorB
 enddo global_spin_loop

 call stutdown_imp_dervs()

   call deallocate_intermediates()
!    print*,'cpksdervs_xc MEMSET'
!    MEMSET(0)

  contains

     subroutine s1_impgrads(nuc_imp_grarho)
     real(kind=r8_kind),intent(inout) :: nuc_imp_grarho(:,:)

           ! ** nuc_imp_grarho 

       allocate(PhiSkl(vec_length,occ_dim), stat=cpksalloc(96))
        ASSERT(cpksalloc(96).eq.0)
        MEMLOG(size(PhiSkl))
 
      if(nl_calc) then
       allocate(graPhiSkl(vec_length,3,occ_dim), stat=alloc_stat(30))
        ASSERT(alloc_stat(30).eq.0)
        MEMLOG(size(graPhiSkl))
      endif


        grads_s1_nuc_imp_grarho: do i_grad=1,size(cpks,1)   

        ! ** calculate phi_occ * s1

        call dgemm('n','n', vec_length,occ_dim,occ_dim, &
                   real(ispin-3,r8_kind),phi_p(:,:,l),full_vl,  &
                   cpks(i_grad,i_irr,spin)%s1,occ_dim, &       !%s1 used to calc nuc_imp_grarho
                   0.0_r8_kind,PhiSkl,vec_length)
       if(nl_calc) then
        call dgemm('n','n', 3*vec_length,occ_dim,occ_dim, &
                   real(ispin-3,r8_kind),graphi(i_irr,spin)%m(:vec_length,:3,:,l),3*vec_length,  &
                   cpks(i_grad,i_irr,spin)%s1,occ_dim, &       !%s1 used to calc nuc_imp_grarho
                   0.0_r8_kind,graPhiSkl,vec_length*3)
       endif

         
           !   here (igrad) in cpks(i_grad,i)%S1 is one of total symmetric projections
           !   and thus result summed up to related i_grad projectio of Qai as for 
           !   the coulomb  and other analitic contribs
           !   but everadging of XC contrib probable will be required as in 
           !   gradient posthoc calc (see post_hoc_average_gradient route)
           !   and may be note for this term as we calculate it as integral 
           !   and gradient   information is inhereted from %s1
        ! ** phi_occ * s1 calculated


       ! ** calculate phi_occ * s1 * phi_occ
       nuc_imp_grarho(:vec_length,i_grad)=  &            ! 1st of (2) nuc_imp_grarho contribs
       nuc_imp_grarho(:vec_length,i_grad) + &            ! contribs from both spin systems
           sum(PhiSkl(:vec_length,:occ_dim)*phi_p(:vec_length,:occ_dim,l),2)

        if(nl_calc) then
         nuc_imp_gragamma(:vec_length,:3,i_grad)= &
         nuc_imp_gragamma(:vec_length,:3,i_grad)  &
       + sum( spread(PhiSkl(:vec_length,:occ_dim),2,3) &
              *graphi(i_irr,spin)%m(:vec_length,:3,:occ_dim,l), 3) &
       + sum( graPhiSkl(:vec_length,:3,:occ_dim) * &
               spread(phi_p(:vec_length,:occ_dim,l),2,3),3)
        endif

       enddo grads_s1_nuc_imp_grarho

       MEMLOG(-size(PhiSkl))
       deallocate(PhiSkl,stat=cpksalloc(96))
       ASSERT(cpksalloc(96).eq.0)
       cpksalloc(96)=1
       if(nl_calc) then
        MEMLOG(-size(graPhiSkl))
        deallocate(graPhiSkl,stat=alloc_stat(30))
        ASSERT(alloc_stat(30).eq.0)
        alloc_stat(30)=1
       endif
     end subroutine s1_impgrads

      subroutine impgrads(nuc_imp_grarho)
      real(kind=r8_kind), intent(inout):: nuc_imp_grarho(:,:)

      integer(kind=i4_kind):: k,stat
#if 0
      real(kind=r8_kind):: BPhi(vec_length,occ_dim)
#else
      real(kind=r8_kind),allocatable:: BPhi(:,:)
      allocate(BPhi(vec_length,occ_dim),stat=stat)
      ASSERT(stat.eq.0)
#endif

       allocate( PhiB(vec_length,vir_dim),stat=cpksalloc(114))
           ASSERT(cpksalloc(114).eq.0)
           MEMLOG(size(PhiB))

       if(nl_calc) then
        allocate( graPhiB(vec_length,3,vir_dim),stat=alloc_stat(29))
        ASSERT(alloc_stat(29).eq.0)
        MEMLOG(size(graPhiB))
       endif


      b_imp_grarho_partners: do l=1,partners(i_irr)
      b_imp_grarho_grads: do i_grad=1,size(cpks,1)
       
       call dgemm( 'n','n', vec_length,vir_dim,occ_dim, &
                   two*(3-ispin),phi_p(:,:,l),full_vl, &
                   cpks(i_grad,i_irr,spin)%HBH,occ_dim, &
                   0.0_r8_kind,PhiB,vec_length)    ! PhiB(2)

        !!!       (Phi(:occ)*B(:occ,:vac)) * Phi(:vac)

       nuc_imp_grarho(:vec_length,i_grad)= &  ! (2) B contrib
       nuc_imp_grarho(:vec_length,i_grad) + sum( PhiB(:vec_length,:vir_dim)* &
                       phi_p(:vec_length,eig_dim-vir_dim+1:,l),2)


    if(nl_calc) then ! compare with s1 term

        !!!   (Phi(:occ)*B(:occ,:vac)) * graphi(:vac)
         nuc_imp_gragamma(:vec_length,:3,i_grad)=nuc_imp_gragamma(:vec_length,:3,i_grad) &
           +sum( spread(PhiB(:,:vir_dim),2,3) &
                *graphi(i_irr,spin)%m(:vec_length,:3,1+eig_dim-vir_dim:,l), 3)

        call dgemm( 'n','t', vec_length, occ_dim, vir_dim, real(2*(3-ispin),r8_kind), &
                    phi_p(1,1+eig_dim-vir_dim,l),full_vl, cpks(i_grad,i_irr,spin)%HBH,occ_dim, &
                    zero,BPhi,vec_length)
        do k=1,3
         nuc_imp_gragamma(:vec_length,k,i_grad)=nuc_imp_gragamma(:vec_length,k,i_grad) &
           +sum( graphi(i_irr,spin)%m(:vec_length,k,:occ_dim,l)*BPhi,2)
        enddo


      endif
        
    enddo b_imp_grarho_grads
  enddo b_imp_grarho_partners


       MEMLOG(-size(PhiB))
       deallocate(PhiB,stat=cpksalloc(114))
       ASSERT(cpksalloc(114).eq.0)
       cpksalloc(114)=1

       if(nl_calc) then
        MEMLOG(-size(graPhiB))
        deallocate( graPhiB,stat=alloc_stat(29))
        ASSERT(alloc_stat(29).eq.0)
        alloc_stat(29)=1
       endif
#if 1
      deallocate(BPhi,stat=stat)
      ASSERT(stat.eq.0)
#endif
  end subroutine impgrads

     subroutine s1_dervs()
#if 0
     real(kind=r8_kind):: im_help_nuc_occ(vec_length,3,unique_atoms(i_m)%n_equal_atoms,occ_dim)
     real(kind=r8_kind):: im_help_sec_der_occ(vec_length,6,unique_atoms(i_m)%n_equal_atoms,occ_dim)
#else
     integer(i4_kind):: stat
     real(kind=r8_kind),allocatable :: im_help_nuc_occ(:,:,:,:)
     real(kind=r8_kind),allocatable :: im_help_sec_der_occ(:,:,:,:)
     allocate(im_help_nuc_occ(vec_length,3,unique_atoms(i_m)%n_equal_atoms,occ_dim),&
              im_help_sec_der_occ(vec_length,6,unique_atoms(i_m)%n_equal_atoms,occ_dim), &
              stat=stat)
     ASSERT(stat.eq.0)
#endif

     im_help_nuc_occ=zero
     if(nl_calc) im_help_sec_der_occ=zero

         allocate(im_help_nuc_arr(vec_length,3,occ_dim),stat=alloc_stat(31))
         ASSERT(alloc_stat(31).eq.0)
         im_help_nuc_arr=0.0_r8_kind
         if(nl_calc) then
           allocate(im_help_sec_der_arr(vec_length,6,occ_dim),stat=alloc_stat(32))
           ASSERT(alloc_stat(32).eq.0)
           im_help_sec_der_arr=0.0_r8_kind
         endif

       counter=0

       nuc_gra_ua:  do i_ua=1,n_unique_atoms

       imp_dervsrho(i_ua)%m=0.0_r8_kind ! for s1 contrib
       if(nl_calc) imp_dervs_grarho(i_ua)%m=0.0_r8_kind

         orb_dim=size(nuc_grads(i_ua,i_irr)%o,4)

         orb_dim_blk: if (orb_dim > 0) then


        if(nl_calc) then ! compare with b part
         call alloc_dervs_helps(occ_dim)
         call calc_dervs_helps(help_sec_der_arr=help_sec_der_arr)
        else
         call alloc_dervs_helps(occ_dim)
         call calc_dervs_helps()
        endif

          if(i_ua.ne.i_m) then
           if(nl_calc) then
              call dervs_imp4(i_ua,help_nuc_arr,help_2nd=help_sec_der_arr)
              im_help_sec_der_occ(:vec_length,:6,1,:)= im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                     - sum(help_sec_der_arr(:vec_length,:6,:,:),3)
           else
              call dervs_imp4(i_ua,help_nuc_arr)
           endif
              im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                   -sum(help_nuc_arr(:vec_length,:3,:,:),3)
          elseif(size(im_help_nuc_occ,3).ne.1) then
            if(nl_calc) then
               im_help_sec_der_occ(:vec_length,:6,2:,:)=help_sec_der_arr(:vec_length,:6,2:,:) 
               im_help_sec_der_occ(:vec_length,:6, 1,:)=im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                       -sum(help_sec_der_arr(:vec_length,:6,2:,:),3)
            endif
                im_help_nuc_occ(:vec_length,:3,2:,:)=help_nuc_arr(:vec_length,:3,2:,:)
                im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                 -sum(help_nuc_arr(:vec_length,:3,2:,:),3)
          endif

         do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle
           im_help_nuc_arr(:vec_length,:3,:)=im_help_nuc_arr(:vec_length,:3,:) &
            - help_nuc_arr(:vec_length,:3,i_ea,:)

           if(nl_calc) then
           im_help_sec_der_arr(:vec_length,:6,:)=im_help_sec_der_arr(:vec_length,:6,:) &
           -help_sec_der_arr(:vec_length,:6,i_ea,:)
           endif

         enddo

         help_2nd=>help_sec_der_arr
         counter=counter+orb_dim

         call dealloc_dervs_helps()

         endif orb_dim_blk
 
      if(i_ua.ne.i_m) then
        call transform_to_cart(dervsrho_imp,imp_dervsrho,i_ua,spin) ! (1 - s1)
        if(nl_calc) call transform_to_cart(dervs_grarho_imp,imp_dervs_grarho,i_ua,spin) ! (1)
      endif

       enddo nuc_gra_ua 

     if(nl_calc) then
              call dervs_imp4(i_m,im_help_nuc_occ,help_2nd=im_help_sec_der_occ)
     else
              call dervs_imp4(i_m,im_help_nuc_occ)
     endif
       
         i_ua=i_m
         i_ea=1

        call transform_to_cart(dervsrho_imp,imp_dervsrho,i_m,spin) ! (1 - s1)

        if(nl_calc) then
         call transform_to_cart(dervs_grarho_imp,imp_dervs_grarho,i_m,spin) ! (1)
         deallocate(im_help_sec_der_arr,stat=alloc_stat(32))
         ASSERT(alloc_stat(32).eq.0)
         alloc_stat(32)=1
        endif
       deallocate(im_help_nuc_arr,stat=alloc_stat(31))
       ASSERT(alloc_stat(31).eq.0)
       alloc_stat(31)=1
#if 1
     deallocate(im_help_nuc_occ,im_help_sec_der_occ,stat=stat)
     ASSERT(stat.eq.0)
#endif
  end subroutine s1_dervs

  subroutine all_imp_dervs()
  
#if 0
    real(kind=r8_kind):: im_help_nuc_occ(vec_length,3,unique_atoms(i_m)%n_equal_atoms,occ_dim)
    real(kind=r8_kind):: im_help_nuc_vir(vec_length,3,unique_atoms(i_m)%n_equal_atoms,vir_dim)
    real(kind=r8_kind):: im_help_sec_der_occ(vec_length,6,unique_atoms(i_m)%n_equal_atoms,occ_dim)
    real(kind=r8_kind):: im_help_sec_der_vir(vec_length,6,unique_atoms(i_m)%n_equal_atoms,vir_dim)
#else
    real(kind=r8_kind),allocatable:: im_help_nuc_occ(:,:,:,:)
    real(kind=r8_kind),allocatable:: im_help_nuc_vir(:,:,:,:)
    real(kind=r8_kind),allocatable:: im_help_sec_der_occ(:,:,:,:)
    real(kind=r8_kind),allocatable:: im_help_sec_der_vir(:,:,:,:)
    integer(kind=i4_kind):: stat

    allocate(im_help_nuc_occ(vec_length,3,unique_atoms(i_m)%n_equal_atoms,occ_dim), &
             im_help_nuc_vir(vec_length,3,unique_atoms(i_m)%n_equal_atoms,vir_dim), &
             im_help_sec_der_occ(vec_length,6,unique_atoms(i_m)%n_equal_atoms,occ_dim), &
             im_help_sec_der_vir(vec_length,6,unique_atoms(i_m)%n_equal_atoms,vir_dim), &
             stat=stat)
    ASSERT(stat.eq.0)
#endif

 do l=1,partners(i_irr)
    im_help_nuc_occ=zero
    im_help_nuc_vir=zero
    if(nl_calc) then
     im_help_sec_der_occ=zero
     im_help_sec_der_vir=zero
    endif

         im_help_nuc_arr=0.0_r8_kind
         vim_help_nuc_arr=0.0_r8_kind
         if(nl_calc) then
          im_help_sec_der_arr=0.0_r8_kind
          vim_help_sec_der_arr=0.0_r8_kind
         endif

       counter=0
       b_gra_ua:  do i_ua=1,n_unique_atoms

         orb_dim=size(nuc_grads(i_ua,i_irr)%o,4)

         imp_dervsrho(i_ua)%m=0.0_r8_kind ! for b call

         if(nl_calc) imp_dervs_grarho(i_ua)%m=0.0_r8_kind

         b_orb_dim: if (orb_dim > 0) then

         occ_dim2: if(occ_dim.gt.0) then


         if(nl_calc) then ! compare with s1 part
          call alloc_dervs_helps(occ_dim)
          call calc_dervs_helps(help_sec_der_arr=help_sec_der_arr)      
         else
          call alloc_dervs_helps(occ_dim)
          call calc_dervs_helps()
         endif

          if(i_ua.ne.i_m) then
           if(nl_calc) then
              call dervs_imp2(i_ua,help_nuc_arr,help_2nd=help_sec_der_arr)
              im_help_sec_der_occ(:vec_length,:6,1,:)= im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                     - sum(help_sec_der_arr(:vec_length,:6,:,:),3)
           else
              call dervs_imp2(i_ua,help_nuc_arr)
           endif
              im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                   -sum(help_nuc_arr(:vec_length,:3,:,:),3)
          elseif(size(im_help_nuc_occ,3).ne.1) then
            if(nl_calc) then
               im_help_sec_der_occ(:vec_length,:6,2:,:)=help_sec_der_arr(:vec_length,:6,2:,:) 
               im_help_sec_der_occ(:vec_length,:6, 1,:)=im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                       -sum(help_sec_der_arr(:vec_length,:6,2:,:),3)
            endif
                im_help_nuc_occ(:vec_length,:3,2:,:)=help_nuc_arr(:vec_length,:3,2:,:)
                im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                 -sum(help_nuc_arr(:vec_length,:3,2:,:),3)
          endif

         if(nl_calc) help_2nd=>help_sec_der_arr

    do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle

           im_help_nuc_arr(:vec_length,:3,:)=im_help_nuc_arr(:vec_length,:3,:) &
            - help_nuc_arr(:vec_length,:3,i_ea,:)
           if(nl_calc) then
           im_help_sec_der_arr(:vec_length,:6,:)=im_help_sec_der_arr(:vec_length,:6,:) &
            - help_sec_der_arr(:vec_length,:6,i_ea,:) ! :6 is correct, I assume?
           endif
     enddo

     call dealloc_dervs_helps()

         if(nl_calc) then
         call alloc_dervs_helps(vir_dim) ! help_sec_der_arr(vir_dim)
         call calc_dervs_helps(eig_dim-vir_dim,help_sec_der_arr=help_sec_der_arr)
         else
         call alloc_dervs_helps(vir_dim) 
         call calc_dervs_helps(eig_dim-vir_dim)
         endif

          if(i_ua.ne.i_m) then
           if(nl_calc) then
              call dervs_imp3(i_ua,help_nuc_arr,help_2nd=help_sec_der_arr) !(1)
              im_help_sec_der_vir(:vec_length,:6,1,:)= im_help_sec_der_vir(:vec_length,:6,1,:) &
                                                     - sum(help_sec_der_arr(:vec_length,:6,:,:),3)
           else
              call dervs_imp3(i_ua,help_nuc_arr) !(2)
           endif
              im_help_nuc_vir(:vec_length,:3,1,:)=im_help_nuc_vir(:vec_length,:3,1,:)   &
                                                   -sum(help_nuc_arr(:vec_length,:3,:,:),3)
          elseif(size(im_help_nuc_vir,3).ne.1) then
            if(nl_calc) then
               im_help_sec_der_vir(:vec_length,:6,2:,:)=help_sec_der_arr(:vec_length,:6,2:,:) 
               im_help_sec_der_vir(:vec_length,:6, 1,:)=im_help_sec_der_vir(:vec_length,:6,1,:) &
                                                       -sum(help_sec_der_arr(:vec_length,:6,2:,:),3)
            endif
                im_help_nuc_vir(:vec_length,:3,2:,:)=help_nuc_arr(:vec_length,:3,2:,:)
                im_help_nuc_vir(:vec_length,:3,1,:)=im_help_nuc_vir(:vec_length,:3,1,:)   &
                                                 -sum(help_nuc_arr(:vec_length,:3,2:,:),3)
          endif

         do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle
           vim_help_nuc_arr(:vec_length,:3,:)=vim_help_nuc_arr(:vec_length,:3,:) &
            - help_nuc_arr(:vec_length,:3,i_ea,:)
           if(nl_calc) then
           vim_help_sec_der_arr(:vec_length,:6,:)=vim_help_sec_der_arr(:vec_length,:6,:) &
            - help_sec_der_arr(:vec_length,:6,i_ea,:) ! :6 is correct, I assume?
           endif
         enddo

         if(nl_calc) help_2nd=>help_sec_der_arr

      endif occ_dim2

                counter=counter+orb_dim

         call dealloc_dervs_helps()


          endif b_orb_dim


         if(i_m.ne.i_ua) then
          call transform_to_cart(dervsrho_imp,imp_dervsrho,i_ua,spin) !(2 b)
          if(nl_calc) call transform_to_cart(dervs_grarho_imp,imp_dervs_grarho,i_ua,spin) ! (1)
         endif

       enddo b_gra_ua

       call dervs_imp2(i_m,im_help_nuc_occ,help_2nd=im_help_sec_der_occ)
       call dervs_imp3(i_m,im_help_nuc_vir,help_2nd=im_help_sec_der_vir) !(3)

         i_ua=i_m
         i_ea=1

        call transform_to_cart(dervsrho_imp,imp_dervsrho,i_m,spin) 
        if(nl_calc) call transform_to_cart(dervs_grarho_imp,imp_dervs_grarho,i_m,spin) ! (1)

  enddo
    deallocate(im_help_nuc_occ,im_help_nuc_vir,im_help_sec_der_occ, &
               im_help_sec_der_vir, stat=stat)
    ASSERT(stat.eq.0)
 end subroutine all_imp_dervs

       subroutine dervs_imp3(i_ua,help_nuc_arr,help_2nd)
       integer(kind=i4_kind),intent(in)::i_ua
       real(kind=r8_kind), intent(in):: help_nuc_arr(:,:,:,:) 
       real(kind=r8_kind), intent(in),optional:: help_2nd(:,:,:,:) 

       integer(kind=i4_kind):: k,i_occ,j_vir,i_grad,off,i_ea
       integer(kind=i4_kind):: kk(3,3)=reshape( (/1,2,3, 2,4,5, 3,5,6/),(/3,3/) )
#if 0
       real(kind=r8_kind):: rho_dervs(occ_dim,vir_dim)
       real(kind=r8_kind):: phi_adfac(vir_dim,vec_length),occ(vec_length,occ_dim)
#else
      integer(i4_kind):: stat
      real(kind=r8_kind),allocatable:: phi_adfac(:,:),occ(:,:),rho_dervs(:,:)
      allocate(phi_adfac(vir_dim,vec_length),occ(vec_length,occ_dim), &
               rho_dervs(occ_dim,vir_dim),stat=stat)
      ASSERT(stat.eq.0)
#endif
       off=eig_dim-vir_dim

   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
!          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle


        do k=1,3
         do j_vir=1,vir_dim
         phi_adfac(j_vir,:)=real(3-ispin,r8_kind)*two* &
                            grdwts(:vec_length)*dfdrho(:vec_length,spin)* &
                            help_nuc_arr(:vec_length,k,i_ea,j_vir)
         enddo

          call dgemm('t', 't', occ_dim, vir_dim, vec_length, &
                   one,phi_p(1, 1, l), full_vl, phi_adfac, vir_dim, &
                   zero, rho_dervs, occ_dim)

         
         do i_grad=1,size(cpks,1)
          imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)+&
                                                   sum(cpks(i_grad,i_irr,spin)%HBH*rho_dervs) !(2 b)
         enddo

        enddo !k

        if(nl_calc) then
        do k=1,3
         do j_vir=1,vir_dim
           help_sd(:vec_length,k,1)=help_2nd(:vec_length,kk(1,k),i_ea,j_vir) 
           help_sd(:vec_length,k,2)=help_2nd(:vec_length,kk(2,k),i_ea,j_vir) 
           help_sd(:vec_length,k,3)=help_2nd(:vec_length,kk(3,k),i_ea,j_vir) 

         phi_adfac(j_vir,:)=real((3-ispin)*4,r8_kind)* &
                        dfdgrarho(:vec_length, spin) * grdwts(:vec_length) * &
                 sum(grarho(:vec_length,:,spin)*help_sd(:vec_length,k,:),2)
         if(ispin.eq.2) &
         phi_adfac(j_vir,:)=phi_adfac(j_vir,:)+real((3-ispin)*2,r8_kind)* &
                        dfdgrarho(:vec_length, updn) * grdwts(:vec_length) * &
                 sum(grarho(:vec_length,:,3-spin)*help_sd(:vec_length,k,:),2)

         enddo

          call dgemm('t','t',occ_dim,vir_dim,vec_length, &
                   one,phi_p(:,1:,l),full_vl, phi_adfac,vir_dim, &
                   zero,rho_dervs,occ_dim)
        do i_occ=1,occ_dim
         occ(:,i_occ)= real((3-ispin)*4,r8_kind)* &
          dfdgrarho(:vec_length,spin)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,i_occ,l)*grarho(:vec_length,:3,spin),2)
         if(ispin.eq.2) &
         occ(:,i_occ)=occ(:,i_occ)+real((3-ispin)*2,r8_kind)* &
          dfdgrarho(:vec_length,updn)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,i_occ,l)*grarho(:vec_length,:3,3-spin),2)
        enddo
          phi_adfac=transpose(help_nuc_arr(:vec_length,k,i_ea,:))
          call dgemm('t','t',occ_dim,vir_dim,vec_length, &
                   one,occ,vec_length, phi_adfac,vir_dim, &
                   one,rho_dervs,occ_dim)

         do i_grad=1,size(cpks,1)
          imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)+&
          sum(cpks(i_grad,i_irr,spin)%HBH*rho_dervs)
         enddo
         
        enddo !k
       endif

     enddo

#if 1
      deallocate(phi_adfac,occ,rho_dervs,stat=stat)
      ASSERT(stat.eq.0)
#endif

   end subroutine dervs_imp3

       subroutine dervs_imp2(i_ua,help_nuc_arr,help_2nd)
       integer(kind=i4_kind),intent(in)::i_ua
       real(kind=r8_kind), intent(in):: help_nuc_arr(:,:,:,:) 
       real(kind=r8_kind), intent(in),optional:: help_2nd(:,:,:,:) 

       integer(kind=i4_kind):: k,i_occ,j_vir,i_grad,off,i_ea
       integer(kind=i4_kind):: kk(3,3)=reshape( (/1,2,3, 2,4,5, 3,5,6/),(/3,3/) )
#if 0
       real(kind=r8_kind):: phi_adfac(occ_dim,vec_length),vir(vec_length,vir_dim)
       real(kind=r8_kind):: rho_dervs(occ_dim,vir_dim)
#else
       real(kind=r8_kind),allocatable:: phi_adfac(:,:),vir(:,:)
       real(kind=r8_kind),allocatable:: rho_dervs(:,:)
       integer(kind=i4_kind)::stat
       allocate(phi_adfac(occ_dim,vec_length),vir(vec_length,vir_dim), &
                rho_dervs(occ_dim,vir_dim),stat=stat)
       ASSERT(stat.eq.0)
#endif

       off=eig_dim-vir_dim

   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
!          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle


        do k=1,3
         do i_occ=1,occ_dim
         phi_adfac(i_occ,:)=real(3-ispin,r8_kind)*two* &
                            grdwts(:vec_length)*dfdrho(:vec_length,spin)* &
                            help_nuc_arr(:vec_length,k,i_ea,i_occ)
         enddo

          call dgemm('n', 'n', occ_dim, vir_dim, vec_length, &
                   one, phi_adfac, occ_dim, phi_p(1, off+1, l), full_vl, &
                   zero, rho_dervs, occ_dim)

         
         do i_grad=1,size(cpks,1)
          imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)+&
                                                   sum(cpks(i_grad,i_irr,spin)%HBH*rho_dervs) !(2 b)
         enddo

        enddo !k

        if(nl_calc) then
        do k=1,3
         do i_occ=1,occ_dim
           help_sd(:vec_length,k,1)=help_2nd(:vec_length,kk(1,k),i_ea,i_occ) 
           help_sd(:vec_length,k,2)=help_2nd(:vec_length,kk(2,k),i_ea,i_occ) 
           help_sd(:vec_length,k,3)=help_2nd(:vec_length,kk(3,k),i_ea,i_occ) 

         phi_adfac(i_occ,:)=real((3-ispin)*4,r8_kind)* &
                        dfdgrarho(:vec_length, spin) * grdwts(:vec_length) * &
                 sum(grarho(:vec_length,:,spin)*help_sd(:vec_length,k,:),2)
         if(ispin.eq.2) &
         phi_adfac(i_occ,:)=phi_adfac(i_occ,:)+real((3-ispin)*2,r8_kind)* &
                        dfdgrarho(:vec_length, updn) * grdwts(:vec_length) * &
                 sum(grarho(:vec_length,:,3-spin)*help_sd(:vec_length,k,:),2)

         enddo

          call dgemm('n','n',occ_dim,vir_dim,vec_length, &
                   one,phi_adfac,occ_dim, phi_p(:,off+1:,l),full_vl, &
                   zero,rho_dervs,occ_dim)
        do j_vir=1,vir_dim
         vir(:,j_vir)= real((3-ispin)*4,r8_kind)* &
          dfdgrarho(:vec_length,spin)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_vir+off,l)*grarho(:vec_length,:3,spin),2)
         if(ispin.eq.2) &
         vir(:,j_vir)=vir(:,j_vir)+real((3-ispin)*2,r8_kind)* &
          dfdgrarho(:vec_length,updn)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_vir+off,l)*grarho(:vec_length,:3,3-spin),2)
        enddo
          phi_adfac=transpose(help_nuc_arr(:vec_length,k,i_ea,:))
          call dgemm('n','n',occ_dim,vir_dim,vec_length, &
                   one,phi_adfac,occ_dim, vir,vec_length, &
                   one,rho_dervs,occ_dim)

         do i_grad=1,size(cpks,1)
          imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)+&
          sum(cpks(i_grad,i_irr,spin)%HBH*rho_dervs)
         enddo
         
        enddo !k
       endif

     enddo

#if 1
       deallocate(phi_adfac,vir,rho_dervs,stat=stat)
       ASSERT(stat.eq.0)
#endif
   end subroutine dervs_imp2

       subroutine dervs_imp4(i_ua,help_nuc_arr,help_2nd)
       integer(kind=i4_kind),intent(in)::i_ua
       real(kind=r8_kind), intent(in):: help_nuc_arr(:,:,:,:) 
       real(kind=r8_kind), intent(in),optional:: help_2nd(:,:,:,:) 

       integer(kind=i4_kind):: k,i_occ,i_grad,i_ea
       integer(kind=i4_kind):: kk(3,3)=reshape( (/1,2,3, 2,4,5, 3,5,6/),(/3,3/) )
#if 0
       real(kind=r8_kind):: rho_dervs(occ_dim,occ_dim)
       real(kind=r8_kind):: phi_adfac(occ_dim,vec_length),occ(vec_length,occ_dim)
#else
       integer(kind=i4_kind):: stat
       real(kind=r8_kind),allocatable:: rho_dervs(:,:)
       real(kind=r8_kind),allocatable:: phi_adfac(:,:),occ(:,:)
       allocate(rho_dervs(occ_dim,occ_dim),phi_adfac(occ_dim,vec_length), &
                occ(vec_length,occ_dim),stat=stat)
       ASSERT(stat.eq.0)
#endif

   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
!          if(i_ea.eq.1.and.i_m.eq.i_ua) cycle

        do k=1,3
         do i_occ=1,occ_dim
         phi_adfac(i_occ,:)=real((3-ispin)*2,r8_kind)* &
                            grdwts(:vec_length)*dfdrho(:vec_length,spin)* &
                            help_nuc_arr(:vec_length,k,i_ea,i_occ)
         enddo

          call dgemm('n','n',occ_dim,occ_dim,vec_length, &
                   one,phi_adfac,occ_dim, phi_p(:,:,l),full_vl, &
                   zero,rho_dervs,occ_dim)

         
         do i_grad=1,size(cpks,1)
          imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervsrho(i_ua)%m(k,i_ea,spin,i_grad)-&
                                                   sum(cpks(i_grad,i_irr,spin)%s1*rho_dervs) !(2 b)
         enddo

        enddo !k

        if(nl_calc) then
        do k=1,3
         do i_occ=1,occ_dim
           help_sd(:vec_length,k,1)=help_2nd(:vec_length,kk(1,k),i_ea,i_occ) 
           help_sd(:vec_length,k,2)=help_2nd(:vec_length,kk(2,k),i_ea,i_occ) 
           help_sd(:vec_length,k,3)=help_2nd(:vec_length,kk(3,k),i_ea,i_occ) 

         phi_adfac(i_occ, :) = (3 - ispin) * 4 * &
                        dfdgrarho(:vec_length, spin) * grdwts(:vec_length)* &
                 sum(grarho(:vec_length, :, spin) * help_sd(:vec_length, k, :), 2)
         if(ispin.eq.2) &
         phi_adfac(i_occ,:) = phi_adfac(i_occ,:) + (3 - ispin) * 2 * &
                        dfdgrarho(:vec_length, updn) * grdwts(:vec_length)* &
                 sum(grarho(:vec_length, :, 3 - spin) * help_sd(:vec_length, k, :), 2)

         enddo

          call dgemm('n','n',occ_dim,occ_dim,vec_length, &
                   one,phi_adfac,occ_dim, phi_p(:,1:,l),full_vl, &
                   zero,rho_dervs,occ_dim)
        do j_occ=1,occ_dim
         occ(:,j_occ)= real((3-ispin)*4,r8_kind)* &
          dfdgrarho(:vec_length,spin)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_occ,l)*grarho(:vec_length,:3,spin),2)
         if(ispin.eq.2) &
         occ(:,j_occ)=occ(:,j_occ)+real((3-ispin)*2,r8_kind)* &
          dfdgrarho(:vec_length,updn)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_occ,l)*grarho(:vec_length,:3,3-spin),2)
        enddo
          phi_adfac=transpose(help_nuc_arr(:vec_length,k,i_ea,:))
          call dgemm('n','n',occ_dim,occ_dim,vec_length, &
                   one,phi_adfac,occ_dim, occ,vec_length, &
                   one,rho_dervs,occ_dim)

         do i_grad=1,size(cpks,1)
          imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)=imp_dervs_grarho(i_ua)%m(k,i_ea,spin,i_grad)-&
          sum(cpks(i_grad,i_irr,spin)%s1*rho_dervs)
         enddo
         
        enddo !k
       endif

     enddo
#if 1
       deallocate(rho_dervs,phi_adfac,occ,stat=stat)
       ASSERT(stat.eq.0)
#endif
   end subroutine dervs_imp4

       subroutine dervs_imp()

       vir_mo: do j_vir=1,vir_dim
         phi_p_fac=real(3-ispin,r8_kind)*two* &
                   spread( grdwts(:vec_length)*dfdrho(:vec_length,spin)* &
                           phi_p(:vec_length,j_vir+eig_dim-vir_dim,l),2,3) 

        if(nl_calc) then
         nucgra_afac=real(3-ispin,r8_kind)*4.0_r8_kind* &
          spread(dfdgrarho(:vec_length,spin)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_vir+eig_dim-vir_dim,l)* &
              grarho(:vec_length,:3,spin),2),2,3)
         if(ispin.eq.2) &
         nucgra_afac=nucgra_afac+real(3-ispin,r8_kind)*2.0_r8_kind* &
          spread(dfdgrarho(:vec_length,updn)*grdwts(:vec_length)* &
          sum(graphi(i_irr,spin)%m(:vec_length,:3,j_vir+eig_dim-vir_dim,l)* &
              grarho(:vec_length,:3,3-spin),2),2,3)

         nuc2der_afac=real(3-ispin,r8_kind)*4.0_r8_kind* &
          spread(spread(phi_p(:vec_length,j_vir+eig_dim-vir_dim,l)* &
                        dfdgrarho(:vec_length, spin) * grdwts(:vec_length), 2, 3) * &
                 grarho(:vec_length,:,spin),2,3)
         if(ispin.eq.2) &
         nuc2der_afac=nuc2der_afac+real(3-ispin,r8_kind)*2.0_r8_kind* &
          spread(spread(phi_p(:vec_length,j_vir+eig_dim-vir_dim,l)* &
                        dfdgrarho(:vec_length, updn)*grdwts(:vec_length),2,3)* &
                 grarho(:vec_length,:,3-spin),2,3)
       endif

         occ_mo: do i_occ=1,occ_dim

         rho_dervs= sum(phi_p_fac*help_nuc_arr(:vec_length,:3,i_ea,i_occ),1)

         do i_grad=1,size(cpks,1)
          imp_dervsrho(i_ua)%m(:3,i_ea,spin,i_grad)=imp_dervsrho(i_ua)%m(:3,i_ea,spin,i_grad)+&
          cpks(i_grad,i_irr,spin)%HBH(i_occ,j_vir)*rho_dervs !(2 b)
         enddo

         if(nl_calc) then
         grarho_dervs= sum(help_nuc_arr(:,:3,i_ea,i_occ)*nucgra_afac,1)

                      help_sd(:,1,1)=help_2nd(:vec_length,1,i_ea,i_occ)
                      help_sd(:,1,2)=help_2nd(:vec_length,2,i_ea,i_occ)
                      help_sd(:,2,1)=help_sd(:,1,2)
                      help_sd(:,1,3)=help_2nd(:vec_length,3,i_ea,i_occ)
                      help_sd(:,3,1)=help_sd(:,1,3)
                      help_sd(:,2,2)=help_2nd(:vec_length,4,i_ea,i_occ)
                      help_sd(:,2,3)=help_2nd(:vec_length,5,i_ea,i_occ)
                      help_sd(:,3,2)=help_sd(:,2,3)
                      help_sd(:,3,3)=help_2nd(:vec_length,6,i_ea,i_occ)

                      grarho_dervs=grarho_dervs+sum(sum(help_sd*nuc2der_afac,3),1)
                      
         do i_grad=1,size(cpks,1)
          imp_dervs_grarho(i_ua)%m(:3,i_ea,spin,i_grad)=imp_dervs_grarho(i_ua)%m(:3,i_ea,spin,i_grad)+&
          cpks(i_grad,i_irr,spin)%HBH(i_occ,j_vir)*grarho_dervs
         enddo
         endif

        enddo occ_mo
        enddo vir_mo
       end subroutine dervs_imp

        subroutine grarho_rhodgamma(grarho_drhodgamma,grarho_dgdr,df_drhodgamma,vspin) 
        integer(kind=i4_kind),intent(in):: vspin
        real(kind=r8_kind),intent(in):: df_drhodgamma(:,:)
        real(kind=r8_kind),intent(inout)::grarho_drhodgamma(vec_length,3)
        real(kind=r8_kind),intent(inout)::grarho_dgdr(vec_length,3)

        integer(kind=i4_kind):: vl
#if 0
        real(kind=r8_kind)::grarho_dgammadrho(vec_length,3)
#else
        integer(kind=i4_kind):: stat
        real(kind=r8_kind),allocatable::grarho_dgammadrho(:,:)
        allocate(grarho_dgammadrho(vec_length,3),stat=stat)
        ASSERT(stat.eq.0)
#endif

        vl=vec_length

        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &        !!! (2) nl
                            spread(df_drhodgamma(:vl,spin)*grdwts(:vl),2,3)
       if(spin.ne.vspin) then
        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,3-spin)* &       !!! (2)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,5-spin),2,3)

        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)

        grarho_dgammadrho(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &        !!! (3)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,2+spin),2,3)

        grarho_dgammadrho=grarho_dgammadrho+grarho(:vl,:3,3-spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,7-spin),2,3)

        grarho_dgdr=grarho_dgammadrho                                            !!! (3)

       elseif(ispin.eq.2) then
        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,3-spin)/partners(i_irr)* & !!!(2)
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)
        grarho_dgdr=grarho_drhodgamma                                            !!! (3)
       else
        grarho_dgdr=grarho_drhodgamma                                            !!! (3)
       endif
#if 1
        deallocate(grarho_dgammadrho,stat=stat)
        ASSERT(stat.eq.0)
#endif
  end subroutine grarho_rhodgamma

       function h1q_grarho_grarho(df_dgammadgamma,vspin) result(grarho_grarho)
        real(kind=r8_kind),intent(in):: df_dgammadgamma(:,:)
        integer(kind=i4_kind),intent(in):: vspin

        integer(kind=i4_kind):: k,kk,vl
        real(kind=r8_kind):: grarho_grarho(vec_length,3,3)
        vl=vec_length

         do k=1,3
          do kk=1,3
          grarho_grarho(:,k,kk)=  4.0_r8_kind/partners(i_irr)* &        !!!(4)
          grarho(:vl,k,spin)*grarho(:vl,kk,spin)*df_dgammadgamma(:vl,spin)*grdwts(:vl)
          enddo
        enddo
          if(spin.ne.vspin) then
           do k=1,3
            do kk=1,3
            grarho_grarho(:vl,k,kk)=4.0_r8_kind/partners(i_irr)* &           !!!(4.2)
                  grarho(:vl,k,spin)*grarho(:vl,kk,vspin)*grdwts(:vl)*df_dgammadgamma(:vl,ab)
            grarho_grarho(:vl,k,kk)=grarho_grarho(:vl,k,kk) &    !!!(4)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,spin)* &
                                 grarho(:vl,kk,spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)
            grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,vspin)* &
                                 grarho(:vl,kk,vspin)*grdwts(:vl)*df_dgammadgamma(:vl,6-spin)
            grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk)+ &              !!!(4)
                                 grarho(:vl,k,vspin)*grarho(:vl,kk,spin)* &
                                 grdwts(:vl)*df_dgammadgamma(:vl,cc)/partners(i_irr)
            enddo
           enddo

          elseif(ispin.eq.2) then
           do k=1,3
            do kk=1,3
            grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,spin)* &
                                 grarho(:vl,kk,3-spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)
            grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,3-spin)* &
                                 grarho(:vl,kk,spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)

            grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &           !!!(4)
                                +grarho(:vl,k,3-spin)*grarho(:vl,kk,3-spin)* &
                                 grdwts(:vl)*df_dgammadgamma(:vl,cc)/partners(i_irr)
           enddo
          enddo
          endif

        end function h1q_grarho_grarho


   subroutine all_expl_dervs()
    integer(kind=i4_kind):: i
    real(kind=r8_kind):: im_help_nuc_occ(vec_length,3,unique_atoms(i_m)%n_equal_atoms,occ_dim)
    real(kind=r8_kind):: im_help_nuc_vir(vec_length,3,unique_atoms(i_m)%n_equal_atoms,vir_dim)
    real(kind=r8_kind):: im_help_sec_der_occ(vec_length,6,unique_atoms(i_m)%n_equal_atoms,occ_dim)
    real(kind=r8_kind):: im_help_sec_der_vir(vec_length,6,unique_atoms(i_m)%n_equal_atoms,vir_dim)
        real(kind=r8_kind):: grarho_drhodgamma(vec_length,3,ispin),grarho_dgdr(vec_length,3,ispin)
        real(kind=r8_kind):: grarho_grarho(vec_length,3,3,ispin)!!!,dgammadgamma_fac(vec_length,3) 

    im_help_nuc_occ=zero
    im_help_nuc_vir=zero

    if(nl_calc) then
     im_help_sec_der_occ=zero
     im_help_sec_der_vir=zero
      do i=1,ispin
       grarho_grarho(:,:,:,i)=h1q_grarho_grarho(df_dgammadgamma,spin*(3-2*i)+3*i-3)
      enddo

       do i=1,ispin
        call grarho_rhodgamma(grarho_drhodgamma(:,:,i),grarho_dgdr(:,:,i), &
                              df_drhodgamma,spin*(3-2*i)+3*i-3)
       enddo
    endif


    counter=0
    ua:  do i_ua=1,n_unique_atoms
      FPP_TIMER_START(t_dervs_helps)
     orb_dim=size(nuc_grads(i_ua,i_irr)%o,4)
     nonempty_irrua:if (orb_dim > 0) then

                ! ** calculate cpks_grad_xc contrib
                ! w * phi_occ * dvdrho*(d/dN rho) * phi_vac   

           call alloc_dervs_helps(vir_dim)

                ! (matrix of n_vac MO vecors ) x (matrix of basis vunctions on vector grid)

          if(eig_dim.gt.occ_dim) then
             if(nl_calc) then
              call calc_dervs_helps(eig_dim-vir_dim, & !off to take vir MO
                                    help_sec_der_arr=help_sec_der_arr)
             else
              call calc_dervs_helps(eig_dim-vir_dim) ! to calc Q
             endif

                ! *** now help_nuc_arr is gradients of virtual orbitals on the grid 
                !     calculated for given i_ua unique
           endif


          if(n_ai.gt.0) then   !!!! futher revision

          if(i_ua.ne.i_m) then
           if(nl_calc) then
              im_help_sec_der_vir(:vec_length,:6,1,:)= im_help_sec_der_vir(:vec_length,:6,1,:) &
                                                     -sum(help_sec_der_arr(:vec_length,:6,:,:),3)

           endif

                im_help_nuc_vir(:vec_length,:3,1,:)=im_help_nuc_vir(:vec_length,:3,1,:)   &
                                                   -sum(help_nuc_arr(:vec_length,:3,:,:),3)
          elseif(size(im_help_nuc_occ,3).ne.1) then
            if(nl_calc) then
               im_help_sec_der_vir(:vec_length,:6,2:,:)=help_sec_der_arr(:vec_length,:6,2:,:) 
               im_help_sec_der_vir(:vec_length,:6, 1,:)=im_help_sec_der_vir(:vec_length,:6,1,:) &
                                                   -sum(help_sec_der_arr(:vec_length,:6,2:,:),3)
            endif
                im_help_nuc_vir(:vec_length,:3,2:,:)=help_nuc_arr(:vec_length,:3,2:,:)
                im_help_nuc_vir(:vec_length,:3,1,:)=im_help_nuc_vir(:vec_length,:3,1,:)   &
                                                 -sum(help_nuc_arr(:vec_length,:3,2:,:),3)
          endif
    endif
  endif nonempty_irrua

 if(i_ua.ne.i_m) then
  if(nl_calc) then
    call q_dervs_phi_vir(orb_dim.gt.0,i_ua, help_nuc_arr,help_sec_der_arr=help_sec_der_arr, &
                                   nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho, &
                                   grarho_dgdr=grarho_dgdr,grarho_grarho=grarho_grarho)
  else ! i.e. lda
    call q_dervs_phi_vir(orb_dim.gt.0,i_ua, help_nuc_arr) !(1) lda
  endif
 endif

   if(orb_dim.gt.0)   call dealloc_dervs_helps()

   nonempty_irrua_occ: if(orb_dim.gt.0) then
              call alloc_dervs_helps(occ_dim)

           if(occ_dim.gt.0) then
             if(nl_calc) then
              call calc_dervs_helps(help_sec_der_arr=help_sec_der_arr)
             else
              call calc_dervs_helps() ! to calc h1 and Q
             endif
            endif

          if(i_ua.ne.i_m) then
           if(nl_calc) then
              im_help_sec_der_occ(:vec_length,:6,1,:)= im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                     - sum(help_sec_der_arr(:vec_length,:6,:,:),3)
           endif

                im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                   -sum(help_nuc_arr(:vec_length,:3,:,:),3)
          elseif(size(im_help_nuc_occ,3).ne.1) then
            if(nl_calc) then
               im_help_sec_der_occ(:vec_length,:6,2:,:)=help_sec_der_arr(:vec_length,:6,2:,:) 
               im_help_sec_der_occ(:vec_length,:6, 1,:)=im_help_sec_der_occ(:vec_length,:6,1,:) &
                                                   -sum(help_sec_der_arr(:vec_length,:6,2:,:),3)
            endif
                im_help_nuc_occ(:vec_length,:3,2:,:)=help_nuc_arr(:vec_length,:3,2:,:)
                im_help_nuc_occ(:vec_length,:3,1,:)=im_help_nuc_occ(:vec_length,:3,1,:)   &
                                                 -sum(help_nuc_arr(:vec_length,:3,2:,:),3)
          endif


                counter=counter+orb_dim
 endif nonempty_irrua_occ

    FPP_TIMER_STOP(t_dervs_helps)

   if(i_ua.ne.i_m) then
    FPP_TIMER_START(t_h1q_dvdrho)
    if(nl_calc) then
    call  nuc_dervs_blas_nl(orb_dim.gt.0,i_ua,grarho_drhodgamma,grarho_dgdr,grarho_grarho, &
                         graw=graw,nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho,   &
                         help_nuc_arr=help_nuc_arr,help_sec_der_arr=help_sec_der_arr)  !!! (6+7)
    call  h1_dervs_phi_occ(orb_dim.gt.0,i_ua, help_nuc_arr=help_nuc_arr,help_sec_der_arr=help_sec_der_arr, &
                        nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho, &
                        grarho_dgdr=grarho_dgdr,grarho_grarho=grarho_grarho)

    else ! i.e.lda
     call  nuc_dervs_blas_nl(orb_dim.gt.0,i_ua,  graw=graw,nuc_grarho=nuc_grarho, &
                            help_nuc_arr=help_nuc_arr)  !!! (6+7)
     call  h1_dervs_phi_occ(orb_dim.gt.0,i_ua, help_nuc_arr=help_nuc_arr)   !(1-lda)
    endif
    FPP_TIMER_STOP(t_h1q_dvdrho)
   endif
 
    if (orb_dim > 0) call dealloc_dervs_helps()

  enddo ua      

 FPP_TIMER_START(t_h1q_dvdrho)
   if(nl_calc) then
    call  q_dervs_phi_vir(.true.,i_m, im_help_nuc_vir,help_sec_der_arr=im_help_sec_der_vir, &
                        nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho, &
                        grarho_dgdr=grarho_dgdr,grarho_grarho=grarho_grarho)
    call  nuc_dervs_blas_nl(.true.,i_m, grarho_drhodgamma,grarho_dgdr,grarho_grarho, &
                         graw=graw,nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho,  &
                         help_nuc_arr=im_help_nuc_occ,help_sec_der_arr=im_help_sec_der_occ)  !!! (6+7)
    call  h1_dervs_phi_occ(.true.,i_m, im_help_nuc_occ,help_sec_der_arr=im_help_sec_der_occ, &
                        nuc_grarho=nuc_grarho,nuc_sec_derrho=nuc_sec_derrho, &
                        grarho_dgdr=grarho_dgdr,grarho_grarho=grarho_grarho)

   else ! i.e. lda

    call  q_dervs_phi_vir(.true.,i_m, im_help_nuc_vir)  !(2-lda)
    call  nuc_dervs_blas_nl(.true.,i_m, graw=graw,nuc_grarho=nuc_grarho,  &
                            help_nuc_arr=im_help_nuc_occ)  !!! (6+7)
    call  h1_dervs_phi_occ(.true.,i_m, help_nuc_arr=im_help_nuc_occ) !(2-lda)
   endif
  FPP_TIMER_STOP(t_h1q_dvdrho)
 end subroutine all_expl_dervs

 subroutine stutdown_imp_dervs()
       do i_ua=1,n_unique_atoms
        MEMLOG(-size(imp_dervsrho(i_ua)%m))
        deallocate(imp_dervsrho(i_ua)%m,stat=cpksalloc(98))
        ASSERT(cpksalloc(98).eq.0)
        cpksalloc(98)=1
        if(nl_calc) then
         MEMLOG(-size(imp_dervs_grarho(i_ua)%m))
         deallocate(imp_dervs_grarho(i_ua)%m,stat=cpksalloc(154))
         ASSERT(cpksalloc(154).eq.0)
         cpksalloc(154)=1
        endif
       enddo

       MEMLOG(-size(imp_dervsrho))
       deallocate(imp_dervsrho,stat=cpksalloc(99))
       ASSERT(cpksalloc(99).eq.0)
       cpksalloc(99)=1

  if(nl_calc) then
       MEMLOG(-size(imp_dervs_grarho))
       deallocate(imp_dervs_grarho,stat=cpksalloc(153))
       ASSERT(cpksalloc(153).eq.0)
       cpksalloc(153)=1

       MEMLOG(-size(nuc_imp_gragamma))
       deallocate(nuc_imp_gragamma,stat=cpksalloc(144))
       ASSERT(cpksalloc(144).eq.0)
       cpksalloc(144)=1
  endif
 end subroutine  stutdown_imp_dervs

  subroutine setup_imp_dervs()
   allocate(imp_dervsrho(n_unique_atoms),stat=cpksalloc(99))
   ASSERT(cpksalloc(99).eq.0)
   MEMLOG(size(imp_dervsrho))

   if(nl_calc) then
    allocate(imp_dervs_grarho(n_unique_atoms),stat=cpksalloc(153))
    ASSERT(cpksalloc(153).eq.0)
    MEMLOG(size(imp_dervs_grarho))
   endif

   do i_ua=1,n_unique_atoms

    nuc_grarho_imp(i_ua)%m=0.0_r8_kind
    if(nl_calc) nuc_gragamma_imp(i_ua)%m=0.0_r8_kind ! for both s1 and b calcs
    

    allocate( imp_dervsrho(i_ua)%m(3,unique_atoms(i_ua)%n_equal_atoms,ispin,size(cpks,1)), &
              stat=cpksalloc(98))
    MEMLOG(size(imp_dervsrho(i_ua)%m))

    imp_dervsrho(i_ua)%m=0.0_r8_kind ! for each s or b call

    if(nl_calc) then
     allocate(imp_dervs_grarho(i_ua)%m(3,unique_atoms(i_ua)%n_equal_atoms,ispin,size(cpks,1)), &
             stat=cpksalloc(154))
     ASSERT(cpksalloc(154).eq.0)
     MEMLOG(size(imp_dervs_grarho(i_ua)%m))
     imp_dervs_grarho(i_ua)%m=0.0_r8_kind
    endif

   if(present(dervsrho_imp)) then
    do i_ua2=1,n_unique_atoms
    dervsrho_imp(i_ua,i_ua2)%m=0.0_r8_kind
    if(nl_calc) dervs_grarho_imp(i_ua,i_ua2)%m=0.0_r8_kind
    enddo
   endif
   enddo

   if(nl_calc) then
     allocate(nuc_imp_gragamma(vec_length,3,size(cpks,1)),stat=cpksalloc(144))
       ASSERT(cpksalloc(144).eq.0)
       MEMLOG(size(nuc_imp_gragamma))
   endif

   !!! make nuc_grarho and nuc_sec_derrho correct
   do spin=1,ispin
             nuc_grarho(i_m)%m(:vec_length,:3,1,spin)=zero
             if(nl_calc) nuc_sec_derrho(i_m)%m(:vec_length,:3,:3,1,spin)=zero
             do i_ua=1,n_unique_atoms
              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
               if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
               nuc_grarho(i_m)%m(:vec_length,:,1,spin)= &
               nuc_grarho(i_m)%m(:vec_length,:,1,spin)-nuc_grarho(i_ua)%m(:vec_length,:,i_ea,spin)
               if(nl_calc)  nuc_sec_derrho(i_m)%m(:vec_length,:3,:3,1,spin)= &
                               nuc_sec_derrho(i_m)%m(:vec_length,:3,:3,1,spin)  &
                             - nuc_sec_derrho(i_ua)%m(:vec_length,:3,:3,i_ea,spin)
              enddo
             enddo
    enddo

  end subroutine setup_imp_dervs

   subroutine deallocate_intermediates()
    MEMLOG(-size(drhodgamma_fac)-size(phi_adfac))
    MEMLOG(-size(phi_p_adfac))
    MEMLOG(-size(graphi_adfac)-size(dgammadgamma_fac))
    MEMLOG(-size(graphi_dgdg_fac)-size(wdfdgrarho)-size(dgrarho_fac))
    MEMLOG(-size(graphi_dgrarho_fac)-size(grarho_drhodgamma))
    MEMLOG(-size(help_sd)-size(phi_drdg_fac)-size(nucgra_afac))
!    MEMLOG(-size(nuc2der_afac)-size(nuc_imp_grarho)-size(grarho_dgrarho))
    MEMLOG(-size(nuc2der_afac)-size(grarho_dgrarho))
    MEMLOG(-size(phi_grarho_dgrarho)-size(graphi_grarho_dgrarho))
    MEMLOG(-size(graw_adfac)-size(phi_wdgrarho))
    MEMLOG(-size(grarho_grarho))
    
    deallocate(drhodgamma_fac,phi_adfac, &
                         phi_p_adfac, &
                         graphi_adfac, &
                         dgammadgamma_fac, &
                         graphi_dgdg_fac,wdfdgrarho,dgrarho_fac,graphi_dgrarho_fac, &
                         grarho_drhodgamma, nucgra_afac,nuc2der_afac, &
                         phi_drdg_fac,help_sd,graw_adfac, &
                         grarho_dgrarho,phi_grarho_dgrarho,graphi_grarho_dgrarho,   &
                         phi_wdgrarho,grarho_grarho, &
                         stat=cpksalloc(145))
    ASSERT(cpksalloc(145).eq.0)
    cpksalloc(145)=1
    if(ispin.eq.2)  then
     MEMLOG(-size(grarho_dgammadrho)*2)
     deallocate(grarho_dgammadrho,phi_dgdr_fac,stat=cpksalloc(145))
     ASSERT(cpksalloc(145).eq.0)
     cpksalloc(145)=1
    endif
   end subroutine deallocate_intermediates

    subroutine allocate_intermediates()
    allocate( drhodgamma_fac(vec_length),phi_adfac(vec_length), &
                         phi_p_adfac(vec_length), &
                         graphi_adfac(vec_length,3), &
                         dgammadgamma_fac(vec_length,3), &
                         graphi_dgdg_fac(vec_length),wdfdgrarho(vec_length), &
                         dgrarho_fac(vec_length,3,ispin), &
                         graphi_dgrarho_fac(vec_length,ispin), &
                         phi_drdg_fac(vec_length,3), &
                         grarho_drhodgamma(vec_length,3), &
                         help_sd(vec_length,3,3), &
                         nucgra_afac(vec_length,3),nuc2der_afac(vec_length,3,3), &
!                         nuc_imp_grarho(vec_length,size(cpks,1)), &
                         graw_adfac(vec_length,3),&
                         grarho_dgrarho(vec_length,3),phi_grarho_dgrarho(vec_length,3), &
                         graphi_grarho_dgrarho(vec_length),phi_wdgrarho(vec_length,ispin), &
                         grarho_grarho(vec_length,3,3), &
                         stat=cpksalloc(145))
    ASSERT(cpksalloc(145).eq.0)
    MEMLOG(vec_length*62)
    if(ispin.eq.2)  then
     allocate(grarho_dgammadrho(vec_length,3),phi_dgdr_fac(vec_length,3),stat=cpksalloc(145))
     ASSERT(cpksalloc(145).eq.0)
     MEMLOG(2*size(grarho_dgammadrho))
    endif
   end subroutine allocate_intermediates

        subroutine nuc_dervs_blas_nl(orb_dim,i_ua,grarho_drhodgamma,grarho_dgdr,grarho_grarho, &
                                     nuc_grarho,nuc_sec_derrho,graw,help_nuc_arr,help_sec_der_arr)
      
        logical, intent(in):: orb_dim
        integer(kind=i4_kind),intent(in) ::  i_ua
        type(arrmat3),optional,intent(in) :: graw(:)
        type(arrmat4),optional,intent(in) :: nuc_grarho(:)
        type(arrmat5),optional,intent(in) :: nuc_sec_derrho(:)
        real(kind=r8_kind),optional,intent(in) :: help_nuc_arr(:,:,:,:)
        real(kind=r8_kind),optional,intent(in) :: help_sec_der_arr(:,:,:,:)
        real(kind=r8_kind),optional,intent(in) :: grarho_grarho(:,:,:,:)
        real(kind=r8_kind),optional,intent(inout),dimension(vec_length,3,2)::grarho_drhodgamma,grarho_dgdr

        integer(kind=i4_kind):: vl,off,k,i_ea,i_occ,j_vir,j_occ,kk,i,vspin
        real(kind=r8_kind):: phi_adfac(occ_dim,vec_length)
        real(kind=r8_kind):: graphi_vir_fac(vec_length,vir_dim),graphi_occ_fac(vec_length,occ_dim)
        real(kind=r8_kind):: ai(occ_dim,vir_dim,3,unique_atoms(i_ua)%n_equal_atoms)
        real(kind=r8_kind):: ii(occ_dim,occ_dim,3,unique_atoms(i_ua)%n_equal_atoms)
        integer(kind=i4_kind):: ind(2,2)=reshape( (/1,3,2,3/),(/2,2/) )

        
        off=eig_dim-vir_dim
        vl=vec_length

 do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      
         do k=1,3
          phi_adfac=zero !!!
  do i_occ=1,occ_dim

#if 1 /* dfdrho contrib */
          phi_adfac(i_occ,:)=phi_p(:vl,i_occ,l)* &                                          ! (6,qh)
          dfdrho(:vl,spin)/partners(i_irr)*graw(i_ua)%m(:vl,k,i_ea)
#endif
#if 1 /* dvdrho contrib */
         do i=1,ispin
          vspin=spin*(3-2*i)+3*i-3
          phi_adfac(i_occ,:)=phi_adfac(i_occ,:)-phi_p(:vl,i_occ,l)* &                       ! (1,lda)
            dvdrho(:vl,ind(i,spin))/partners(i_irr)*grdwts(:vl)*nuc_grarho(i_ua)%m(:vl,k,i_ea,vspin)
         enddo
#endif

         if(nl_calc) then  
#if 1 /* dgrarho contrib */
           phi_adfac(i_occ,:)=phi_adfac(i_occ,:)+graw(i_ua)%m(:vl,k,i_ea)* &                                !!! (7.1)
             sum(grarho_dgrarho*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)
#endif
          do i=1,ispin
            vspin=spin*(3-2*i)+3*i-3
            phi_adfac(i_occ,:)=phi_adfac(i_occ,:) - phi_p(:vl,i_occ,l)*&                                    !!!(2)
                               sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,spin*(3-2*i)+3*i-3)* &
                                   grarho_drhodgamma(:vl,:3,i),2)

            phi_adfac(i_occ,:)=phi_adfac(i_occ,:) &                                                         !!!(3)
                              -sum(grarho_dgdr(:,:,i)*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)* &
                               nuc_grarho(i_ua)%m(:vl,k,i_ea,vspin)
           do kk=1,3
           dgammadgamma_fac(:,kk)=sum(grarho_grarho(:,:,kk,i)*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)        !!!(4.1)
           enddo
           phi_adfac(i_occ,:)=phi_adfac(i_occ,:) &
              -sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)*dgammadgamma_fac(:,:),2)
#if 1 /* dgrarho contrib */
           phi_adfac(i_occ,:)=phi_adfac(i_occ,:)-sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)* &         !!! (5.1)
                 graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)*graphi_dgrarho_fac(:,uord+abs(spin-vspin))
#endif
          enddo
         endif

  if(orb_dim) then
           phi_adfac(i_occ,:)=phi_adfac(i_occ,:)-help_nuc_arr(:vl,k,i_ea,i_occ)* &       !!! (8,lda)
                              dfdrho(:vl,spin)*grdwts(:vl)/partners(i_irr)         

#if 1 /* dfdgrarho contrib */
           if(present(help_sec_der_arr)) &
           phi_adfac(i_occ,:)=phi_adfac(i_occ,:)-dfdgrarho_grads4(i_occ,k,i_ea,help_sec_der_arr)
#endif
  endif

  enddo

          call dgemm('n', 'n', occ_dim, vir_dim, vl, &
                   one, phi_adfac, occ_dim, phi_p(1, off+1, l), full_vl, &
                   zero, ai(:, :, k, i_ea), occ_dim)
          call dgemm('n','n',occ_dim,occ_dim,vl, &
                   one,phi_adfac,occ_dim, phi_p(:,:,l),full_vl, &
                   zero,ii(:,:,k,i_ea),occ_dim)

#if 1 /* dfdgrarho contrib */
       if(nl_calc.and.orb_dim) then
          graphi_occ_fac=zero
         do i_occ=1,occ_dim
          phi_adfac(i_occ,:)=help_nuc_arr(:vl,k,i_ea,i_occ)
         do i=1,ispin
          vspin=spin*(3-2*i)+3*i-3
          graphi_occ_fac(:,i_occ)=graphi_occ_fac(:,i_occ)-real(3-i,r8_kind)/partners(i_irr)* &
                 dfdgrarho(:vl,ind(i,spin))*grdwts(:vl)* &
                 sum(graphi(i_irr,spin)%m(:vl,:3,i_occ,l)*grarho(:vl,:3,vspin),2)
         enddo
         enddo
          graphi_vir_fac=zero
         do i=1,ispin
         vspin=spin*(3-2*i)+3*i-3
         do j_vir=1,vir_dim
          graphi_vir_fac(:,j_vir)=graphi_vir_fac(:,j_vir)-real(3-i,r8_kind)/partners(i_irr)*&
                 dfdgrarho(:vl,ind(i,spin))*grdwts(:vl)* &
                 sum(graphi(i_irr,spin)%m(:vl,:3,off+j_vir,l)*grarho(:vl,:3,vspin),2)
         enddo
         enddo
          call dgemm('n','n',occ_dim,vir_dim,vl, &
                     one,phi_adfac,occ_dim, graphi_vir_fac,vl, &
                     one,ai(:,:,k,i_ea),occ_dim)
          call dgemm('n','n',occ_dim,occ_dim,vl, &
                   one,phi_adfac,occ_dim, graphi_occ_fac,vl, &
                   one,ii(:,:,k,i_ea),occ_dim)
       endif
#endif

    enddo !k
  enddo 

#if 1
     do i_occ=1,occ_dim
      do j_vir=1,vir_dim
         cpks_grad_xc(i_ua)%m(:,:)=ai(i_occ,j_vir,:,:)
         call cpks_average_gradient(i_ua,cpks_grad_xc)           ! expl grad Vxc + graw  contribs
         call cpks_transform_to_totsymQai(i_ua,i_occ,j_vir,spin)         ! (1a)  Qai contrib
       enddo
     enddo
#endif

     do i_occ=1,occ_dim
      do j_occ=1,occ_dim
         cpks_grad_xc(i_ua)%m(:,:)=ii(i_occ,j_occ,:,:)
         call cpks_average_gradient(i_ua,cpks_grad_xc)           ! expl grad Vxc + graw  contribs
         call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)    
       enddo
     enddo

   end subroutine nuc_dervs_blas_nl

        subroutine h1_dervs_phi_occ(orb_dim,i_ua, help_nuc_arr,help_sec_der_arr, &
                                 nuc_grarho,nuc_sec_derrho,grarho_dgdr,grarho_grarho)
      
        logical, intent(in):: orb_dim
        integer(kind=i4_kind),intent(in) ::  i_ua
        real(kind=r8_kind),         intent(in) :: help_nuc_arr(:,:,:,:)
        real(kind=r8_kind),optional,intent(in) :: help_sec_der_arr(:,:,:,:)
        type(arrmat4),optional,intent(in) :: nuc_grarho(:)
        type(arrmat5),optional,intent(in) :: nuc_sec_derrho(:)
        real(kind=r8_kind),optional,intent(in) :: grarho_dgdr(:,:,:),grarho_grarho(:,:,:,:)

        integer(kind=i4_kind):: vl,k,i_ea,i_occ,j_occ,vspin,i
        real(kind=r8_kind):: ii(size(help_nuc_arr,4),size(help_nuc_arr,4),3,unique_atoms(i_ua)%n_equal_atoms)
        real(kind=r8_kind):: phi_occ_fac(vec_length,occ_dim),graphi_occ_fac(vec_length,occ_dim)
        real(kind=r8_kind):: dgammadgamma_fac(vec_length,3) 
        integer(kind=i4_kind):: ind(2,2)=reshape( (/1,3,2,3/),(/2,2/) )

        vl=vec_length

         graphi_h1eqa: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

         do k=1,3

         phi_occ_fac=zero

        if(nl_calc) then
         do j_occ=1,occ_dim
           phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ)+graw(i_ua)%m(:vl,k,i_ea)* &                            !!! (7.1)
             sum(grarho_dgrarho*graphi(i_irr,spin)%m(:vl,:3,j_occ,l),2)
         enddo
        endif
       if(present(grarho_dgdr)) then
        do i=1,ispin
         do j_occ=1,occ_dim
          phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ)   &
                              -sum(grarho_dgdr(:,:,i)*graphi(i_irr,spin)%m(:vl,:3,j_occ,l),2)* & ! (3.1)
                               nuc_grarho(i_ua)%m(:vl,k,i_ea,spin*(3-2*i)+3*i-3)
         enddo
        enddo
       endif

        if(present(nuc_sec_derrho)) then
         do j_occ=1,occ_dim
         do i=1,ispin
         vspin=spin*(3-2*i)+3*i-3
           do kk = 1, 3
              dgammadgamma_fac(:, kk) = &
                   sum(grarho_grarho(:, :, kk, i) * &
                   graphi(i_irr, spin)%m(:vl, :3, j_occ,l), 2) ! (4.1)
           enddo
           phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ) &                                                      !!! (4)
              -sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)*dgammadgamma_fac(:,:),2)
           phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ)-sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)* &     !!! (5.2)
                 graphi(i_irr,spin)%m(:vl,:3,j_occ,l),2)*graphi_dgrarho_fac(:,uord+abs(spin-vspin))
         enddo
         enddo
        endif

    if(orb_dim) then
         do j_occ=1,occ_dim
          phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ) &
            -wdfdrho(:)/partners(i_irr)*help_nuc_arr(:vl,k,i_ea,j_occ) !!! (8,h1)

         enddo
        if(present(help_sec_der_arr)) then
         do j_occ=1,occ_dim
          phi_occ_fac(:,j_occ)=phi_occ_fac(:,j_occ)-dfdgrarho_grads4(j_occ,k,i_ea,help_sec_der_arr)
         enddo
        endif

   endif !orb_dim

          call dgemm('t','n',occ_dim,occ_dim,vl, &
                   one,phi_p(:,:,l),full_vl, phi_occ_fac,vl, &
                   zero,ii(:,:,k,i_ea),occ_dim)

!        if(size(nuc_grads(i_ua,i_irr)%o,4).ne.0.and.nl_calc)  then
        if(nl_calc.and.orb_dim)  then
          graphi_occ_fac=zero
         do i_occ=1,occ_dim
          phi_occ_fac(:,i_occ)=help_nuc_arr(:vl,k,i_ea,i_occ)
         do i=1,ispin
          vspin=spin*(3-2*i)+3*i-3
          graphi_occ_fac(:,i_occ)=graphi_occ_fac(:,i_occ)-real(3-i,r8_kind)/partners(i_irr)* & !!!(9)
              dfdgrarho(:vl,ind(i,spin))*grdwts(:vl)* &
               sum(graphi(i_irr,spin)%m(:vl,:3,i_occ,l)*grarho(:vl,:3,vspin),2)
         enddo
         enddo
          call dgemm('t','n',occ_dim,occ_dim,vl, &
                   one, graphi_occ_fac,vl,phi_occ_fac,vl, &
                   one,ii(:,:,k,i_ea),occ_dim)
        endif

   enddo !k
  enddo graphi_h1eqa


     do i_occ=1,occ_dim
      do j_occ=1,occ_dim
         cpks_grad_xc(i_ua)%m(:,:)=ii(i_occ,j_occ,:,:)
         call cpks_average_gradient(i_ua,cpks_grad_xc)         
         call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)    
       enddo
     enddo

   end subroutine h1_dervs_phi_occ

        subroutine q_dervs_phi_vir(orb_dim,i_ua, help_nuc_arr,help_sec_der_arr, &
                                 nuc_grarho,nuc_sec_derrho, &
                                 grarho_dgdr,grarho_grarho)
      
        logical, intent(in):: orb_dim
        integer(kind=i4_kind),intent(in) ::  i_ua
        real(kind=r8_kind),optional,intent(in) :: help_nuc_arr(:,:,:,:)
        real(kind=r8_kind),optional,intent(in) :: help_sec_der_arr(:,:,:,:)
        real(kind=r8_kind),optional,intent(inout) :: grarho_grarho(:,:,:,:)
        real(kind=r8_kind),optional,intent(in) :: grarho_dgdr(:,:,:)
        type(arrmat4),optional,intent(in) :: nuc_grarho(:)
        type(arrmat5),optional,intent(in) :: nuc_sec_derrho(:)

        integer(kind=i4_kind):: vl,k,i_ea,i_occ,j_vir,off,vspin,i
        real(kind=r8_kind):: phi_vir_fac(vec_length,vir_dim),graphi_occ_fac(occ_dim,vec_length)
        real(kind=r8_kind):: ai(occ_dim,vir_dim,3,unique_atoms(i_ua)%n_equal_atoms)
        real(kind=r8_kind):: dgammadgamma_fac(vec_length,3) 
        integer(kind=i4_kind):: ind(2,2)=reshape( (/1,3,2,3/),(/2,2/) )
        
        off=eig_dim-vir_dim
        vl=vec_length

         graphi_h1eqa: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

         do k=1,3
         phi_vir_fac=zero
        

       if(nl_calc) then
         do j_vir=1,vir_dim
           phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir)+graw(i_ua)%m(:vl,k,i_ea)* &                            !!! (7.1)
             sum(grarho_dgrarho*graphi(i_irr,spin)%m(:vl,:3,off+j_vir,l),2)
         enddo
       endif
#if 1
       if(present(grarho_dgdr)) then
        do i=1,ispin
         do j_vir=1,vir_dim
            phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir) &                                                  !!! (3.2)
                              -sum(grarho_dgdr(:,:,i)*graphi(i_irr,spin)%m(:vl,:3,off+j_vir,l),2)* &
                               nuc_grarho(i_ua)%m(:vl,k,i_ea,spin*(3-2*i)+3*i-3)
          enddo
         enddo
        endif
#endif

        if(present(nuc_sec_derrho)) then
         do j_vir=1,vir_dim
           do i=1,ispin
           vspin=spin*(3-2*i)+3*i-3
           do kk=1,3
           dgammadgamma_fac(:,kk)=sum(grarho_grarho(:,:,kk,i)*graphi(i_irr,spin)%m(:vl,:3,off+j_vir,l),2)   !!!(4.2)
           enddo
           phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir) &
            -sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)*dgammadgamma_fac(:,:),2)
           phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir)-sum(nuc_sec_derrho(i_ua)%m(:vl,k,:3,i_ea,vspin)* &     !!! (5.2)
                 graphi(i_irr,spin)%m(:vl,:3,off+j_vir,l),2)*graphi_dgrarho_fac(:,uord+abs(spin-vspin))
           enddo
         enddo
        endif
    if(orb_dim) then
         do j_vir=1,vir_dim
                 phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir) &       !lda
                    -wdfdrho(:)/partners(i_irr)*help_nuc_arr(:vl,k,i_ea,j_vir)
         enddo

        if(present(help_sec_der_arr)) then
         do j_vir=1,vir_dim
          phi_vir_fac(:,j_vir)=phi_vir_fac(:,j_vir)-dfdgrarho_grads4(j_vir,k,i_ea,help_sec_der_arr) ! (9-help_sec_der_arr)
         enddo
        endif

  endif

          call dgemm('t','n',occ_dim,vir_dim,vl, &
                   one,phi_p(:,:,l),full_vl, phi_vir_fac,vl, &
                   zero,ai(:,:,k,i_ea),occ_dim)

        
#if 1
     
      if(nl_calc.and.orb_dim) then
          phi_vir_fac(:,:)=help_nuc_arr(:vl,k,i_ea,:)
          graphi_occ_fac=zero
         do i=1,ispin
         vspin=spin*(3-2*i)+3*i-3
         do i_occ=1,occ_dim
          graphi_occ_fac(i_occ,:)=graphi_occ_fac(i_occ,:)-real(3-i,r8_kind)/partners(i_irr)* &   ! (9-help_nuc_arr)
                 dfdgrarho(:vl,ind(i,spin))*grdwts(:vl)* &
                 sum(graphi(i_irr,spin)%m(:vl,:3,i_occ,l)*grarho(:vl,:3,vspin),2)
         enddo
         enddo
          call dgemm('n','n',occ_dim,vir_dim,vl, &
                     one,graphi_occ_fac,occ_dim, phi_vir_fac,vl, &
                     one,ai(:,:,k,i_ea),occ_dim)
      endif
#endif

    enddo !k
 enddo graphi_h1eqa

     do i_occ=1,occ_dim
      do j_vir=1,vir_dim
         cpks_grad_xc(i_ua)%m(:,:)=ai(i_occ,j_vir,:,:)
         call cpks_average_gradient(i_ua,cpks_grad_xc)           ! expl grad Vxc + graw  contribs
         call cpks_transform_to_totsymQai(i_ua,i_occ,j_vir,spin)         ! (1a)  Qai contrib
       enddo
     enddo

   end subroutine q_dervs_phi_vir

        subroutine nuc_dervs_afac(contrib,vspin,nuc_grarho,nuc_sec_derrho,graw)
        character(len=6),intent(in):: contrib
      
        integer(kind=i4_kind),intent(in) ::  vspin
        type(arrmat3),optional,intent(in) :: graw(:)
        type(arrmat4),optional,intent(in) :: nuc_grarho(:)
        type(arrmat5),optional,intent(in) :: nuc_sec_derrho(:)
        integer(kind=i4_kind):: vl

        vl=vec_length


    gw: if(present(graw)) then
     do i_ua=1,n_unique_atoms
                       cpks_grad_xc(i_ua)%m=0.0_r8_kind 
         graphi_h1eqa: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

                         cpks_grad_xc(i_ua)%m(:,i_ea)=cpks_grad_xc(i_ua)%m(:,i_ea)+ &
                                sum(graw(i_ua)%m(:vl,:3,i_ea)*graw_adfac(:vl,:3),1)
         enddo graphi_h1eqa

         call cpks_average_gradient(i_ua,cpks_grad_xc)           ! expl grad Vxc + graw  contribs
         select case(contrib)
         case('q_calc')
                  call cpks_transform_to_totsymQai(i_ua,i_occ,j_vir,spin)         ! (1a)  Qai contrib
         case('h1calc')
                call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)    
         case('h1gphi')
                call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)    
                call h1_transform_to_totsym(i_ua,j_occ,i_occ,spin)   
         end select
     enddo

    else  gw

         do i_ua=1,n_unique_atoms
         cpks_grad_xc(i_ua)%m=0.0_r8_kind 
         if(i_ua.ne.i_m) cpks_grad_xc(i_m)%m=0.0_r8_kind

          do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

            if(present(nuc_grarho)) then
             help_gp(:)=-sum(nuc_grarho(i_ua)%m(:vec_length,:,i_ea,vspin)*nucgra_afac,1)
            endif

             if(present(nuc_sec_derrho)) then
              help_gp(:)=-sum(sum(nuc_sec_derrho(i_ua)%m(:vec_length,:3,:3,i_ea,vspin) &
                                         *nuc2der_afac(:vec_length,:3,:3),3),1)
             endif

             cpks_grad_xc(i_ua)%m(:,i_ea)=cpks_grad_xc(i_ua)%m(:,i_ea)+help_gp
             cpks_grad_xc(i_m )%m(:,1)=cpks_grad_xc(i_m )%m(:,1)-help_gp

          enddo ! i_ea

                 call cpks_average_gradient(i_ua,cpks_grad_xc)           ! expl grad Vxc + graw  contribs
                 select case(contrib)
                 case('q_calc')
                  call cpks_transform_to_totsymQai(i_ua,i_occ,j_vir,spin)         ! (1a)  Qai contrib
                 case('h1calc')
                  call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                 case('h1gphi')
                  call h1_transform_to_totsym(i_ua,i_occ,j_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                  call h1_transform_to_totsym(i_ua,j_occ,i_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                 case default
                  print*,'no such case'
                  stop
                 end select
                 
                 

                 if(i_ua.ne.i_m) then
                  call cpks_average_gradient(i_m,cpks_grad_xc)           
                 select case(contrib)
                 case('q_calc')
                  call cpks_transform_to_totsymQai(i_m,i_occ,j_vir,spin)         ! (1b)  Qai contrib
                 case('h1calc')
                  call h1_transform_to_totsym(i_m,i_occ,j_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                 case('h1gphi')
                  call h1_transform_to_totsym(i_m,j_occ,i_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                  call h1_transform_to_totsym(i_m,i_occ,j_occ,spin)      ! 2nd of 4 contribs calculated for Q call
                 case default
                  print*,'no such case'
                  stop
                 end select
                 endif

         enddo ! i_ua
    endif gw
   end subroutine nuc_dervs_afac


#if 0
     subroutine h1q_func_dervs_reduced(graphi_dgrarho_fac,dvdrho,vspin,graw)
     !!! (2) d^2f/d_rho/d_rho
     real(kind=r8_kind),intent(in):: dvdrho(:),graphi_dgrarho_fac(:,:)
     integer(kind=i4_kind),intent(in):: vspin
     type(arrmat3),optional,intent(in):: graw(:)

     real(kind=r8_kind):: wdvdrho(vec_length),phi_grarho_grarho(vec_length,3,3)
     real(kind=r8_kind),pointer:: grarho_dgdr(:,:),phi_dgdr_f(:,:)

     integer(kind=i4_kind):: j_off,vl
     integer(kind=i4_kind):: ind,k,kk

     vl=vec_length
     j_off=eig_dim-vir_dim

!     wdvdrho(:)=dvdrho(:vl)/partners(i_irr)*grdwts(:vl)                           !!! (1)

 if(nl_calc) then  
        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &         !!! (2)
                                spread(df_drhodgamma(:vl,spin)*grdwts(:vl),2,3)
       if(spin.ne.vspin) then
        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,3-spin)* &       !!! (2)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,5-spin),2,3)

        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)

        grarho_dgammadrho(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &        !!! (3)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,2+spin),2,3)

        grarho_dgammadrho=grarho_dgammadrho+grarho(:vl,:3,3-spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,7-spin),2,3)

        grarho_dgdr=>grarho_dgammadrho                                            !!! (3)

       elseif(ispin.eq.2) then
        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,3-spin)/partners(i_irr)* & !!!(2)
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)
        grarho_dgdr=>grarho_drhodgamma                                            !!! (3)
       else
        grarho_dgdr=>grarho_drhodgamma                                            !!! (3)
       endif
        
       grarho_grarho=h1q_grarho_grarho(df_dgammadgamma,vspin)
     endif
       
      !__________________
       occ_i: do i_occ=1,occ_dim

!         phi_p_fac(:,1)=wdvdrho(:)*phi_p(:vl,i_occ,l)                                      !!!(1)

        if(nl_calc) then

          phi_drdg_fac(:,:3)=grarho_drhodgamma*spread(phi_p(:vl,i_occ,l),2,3)               !!!(2)

          if(spin.eq.vspin) then
           phi_dgdr_f=>phi_drdg_fac                                                         !!! (3)
          else
           phi_dgdr_fac(:,:3)=grarho_dgammadrho(:,:3)*spread(phi_p(:vl,i_occ,l),2,3)  
           phi_dgdr_f=>phi_dgdr_fac
          endif

          graphi_adfac=spread(sum(grarho_dgdr*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2),2,3)                 !!!(3)

!          phi_grarho_grarho=spread(spread(phi_p(:vl,i_occ,l),2,3),3,3)*grarho_grarho             !!!(4.2)
!          forall(k=1:3) &
!          dgammadgamma_fac(:,k)=sum(grarho_grarho(:,:,k)*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2) !!!(4.1)



          dgrarho_fac(:,:3,uord)=spread(graphi_dgrarho_fac(:,uord),2,3)* &                                 !!! (5.1)
                                  graphi(i_irr,spin)%m(:vl,:3,i_occ,l)
          phi_wdgrarho(:,uord)=graphi_dgrarho_fac(:,uord)*phi_p(:vl,i_occ,l)                               !!! (5.2)

          if(spin.ne.vspin) then
           dgrarho_fac(:,:3,ud)=spread(graphi_dgrarho_fac(:,ud),2,3)* &                                    !!! (5.1)
                                graphi(i_irr,spin)%m(:vl,:3,i_occ,l)
           phi_wdgrarho(:,ud)=graphi_dgrarho_fac(:,ud)*phi_p(:vl,i_occ,l)                                  !!! (5.2)
          endif

          phi_grarho_dgrarho(:vl,:3)=grarho_dgrarho(:vl,:3)* &                                             !!! (7.1)
                                     spread(phi_p(:vl,i_occ,l),2,3)
          graphi_grarho_dgrarho(:)=sum(grarho_dgrarho(:,:)* &                                              !!! (7.2)
                                       graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)
         endif

    !_________________________
     q_vir: do j_vir=1,vir_dim

!         nucgra_afac(:vl,:3)=spread(phi_p_fac(:vl,1)*phi_p(:vl,j_vir+j_off,l),2,3)                !!! (1)
!          nucgra_afac(:vl,:3)=zero

!     if(nl_calc) then
!         nucgra_afac=nucgra_afac+graphi_adfac(:vl,:3)*spread(phi_p(:vl,j_vir+j_off,l),2,3)        !!! (3.1)
!         nucgra_afac(:vl,:3)=nucgra_afac(:vl,:3)+ &                                               !!! (3.2) dgdr
!          spread(sum(phi_dgdr_f(:,:3)*graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l),2),2,3)
!      endif
                          ! ________________________
!      call  nuc_dervs_afac('q_calc',vspin,nuc_grarho)

!       if(nl_calc) then
!            nucgra_afac(:vl,:3)= spread(phi_p(:vl,j_vir+j_off,l),2,3)*( &
!              phi_drdg_fac(:vl,:3)+ &                                                          !!! (2)
!              dgammadgamma_fac(:vl,:3) &                                                       !!! (4) 
!               +dgrarho_fac(:vl,:3,uord+abs(spin-vspin)) )                                      !!! (5)

!            nucgra_afac=nucgra_afac+spread(phi_wdgrarho(:,uord+abs(spin-vspin)),2,3)* &         !!! (5)
!                                           graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l)

!            forall(k=1:3,kk=1:3) &
!            nuc2der_afac(:,k,kk)=nucgra_afac(:,kk)! &
!                +sum(graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l)*phi_grarho_grarho(:,:3,kk),2)   !!! (4)
                                !____________________________________________
!            call  nuc_dervs_afac('q_calc',vspin,nuc_sec_derrho=nuc_sec_derrho)

!       endif

!      if(present(graw)) then

!         graw_adfac=spread(dfdrho(:vl,spin)*phi_p(:vl,i_occ,l)* &                        !!! (6)
!                            phi_p(:vl,j_vir+j_off,l),2,3)/partners(i_irr)
!          graw_adfac=zero

!        if(nl_calc)  then
!         graw_adfac=graw_adfac+ spread(&                                                  !!! (7.1)
!          sum(phi_grarho_dgrarho*graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l),2),2,3)
!         graw_adfac=graw_adfac+spread(graphi_grarho_dgrarho*phi_p(:vl,j_vir+j_off,l),2,3) !!! (7.2)
!        endif
                            !________________________
!        call  nuc_dervs_afac('q_calc',vspin,graw=graw)  !!! (6+7)
!      endif

       enddo q_vir

       h1_occ_j: do j_occ=1,occ_dim
    
!         nucgra_afac(:vl,:3)=spread(phi_p_fac(:vl,1)*phi_p(:vl,j_occ,l),2,3)         !!! (1)
!                                !_________________________
!            call  nuc_dervs_afac('h1calc',vspin,nuc_grarho)

!       if(nl_calc) then
!         nucgra_afac=spread(phi_p(:vl,j_occ,l),2,3)*phi_drdg_fac                     !!! (2,h1)
!         nuc2der_afac(:vec_length,:3,:3)=spread(nucgra_afac,2,3)
!                             !____________________________________________
!         call  nuc_dervs_afac('h1calc',vspin,nuc_sec_derrho=nuc_sec_derrho)
!        endif 


!        if(present(graw)) then
!         graw_adfac(:vl,:3)=spread(dfdrho(:vl,spin) &                                 !!! (6,lda)
!           *phi_p(:vl,i_occ,l)*phi_p(:vl,j_occ,l),2,3)/partners(i_irr)
!          
!         call  nuc_dervs_afac('h1calc',vspin,graw=graw)
!        endif

   enddo h1_occ_j 

    if(nl_calc) then
       graphi_h1_occ: do j_occ=1,occ_dim
!        nucgra_afac= graphi_adfac*spread(phi_p(:vl,j_occ,l),2,3)                         !!! (3)
!        call  nuc_dervs_afac('h1gphi',vspin,nuc_grarho)

!        nuc2der_afac(:vl,:3,:3)=spread(spread(phi_p(:vl,j_occ,l),2,3)*&                   !!! (4.2)
!                                       dgammadgamma_fac(:,:3),2,3)                                            
        nuc2der_afac=zero

!        nuc2der_afac=nuc2der_afac+spread(spread(phi_p(:vl,j_occ,l),2,3) * &               !!! (5)
!                                         dgrarho_fac(:,:3,uord+abs(spin-vspin)),2,3)
!        call  nuc_dervs_afac('h1gphi',vspin,nuc_sec_derrho=nuc_sec_derrho)

!       if(present(graw)) then
!        graw_adfac(:vl,:3)= spread( sum(phi_grarho_dgrarho(:vl,:3)* &                     !!! (7)
!                                        graphi(i_irr,spin)%m(:vl,:3,j_occ,l),2),2,3)
!        call  nuc_dervs_afac('h1gphi',vspin,graw=graw)
!       endif

       enddo graphi_h1_occ
    endif

 enddo occ_i

 end subroutine h1q_func_dervs_reduced

     subroutine h1q_func_dervs(graphi_dgrarho_fac,dvdrho,vspin,graw)
     !!! (2) d^2f/d_rho/d_rho
     real(kind=r8_kind),intent(in):: dvdrho(:),graphi_dgrarho_fac(:,:)
     integer(kind=i4_kind),intent(in):: vspin
     type(arrmat3),optional,intent(in):: graw(:)

     real(kind=r8_kind):: wdvdrho(vec_length),phi_grarho_grarho(vec_length,3,3)
     real(kind=r8_kind),pointer:: grarho_dgdr(:,:),phi_dgdr_f(:,:)

     integer(kind=i4_kind):: j_off,vl
     integer(kind=i4_kind):: ind,k,kk

     vl=vec_length
     j_off=eig_dim-vir_dim

     wdvdrho(:)=dvdrho(:vl)/partners(i_irr)*grdwts(:vl)                           !!! (1)

 if(nl_calc) then  
        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &        !!! (2)
                                spread(df_drhodgamma(:vl,spin)*grdwts(:vl),2,3)
       if(spin.ne.vspin) then
        grarho_drhodgamma(:,:3)=two/partners(i_irr)*grarho(:vl,:3,3-spin)* &      !!! (2)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,5-spin),2,3)

        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)

        grarho_dgammadrho(:,:3)=two/partners(i_irr)*grarho(:vl,:3,spin)* &        !!! (3)
                                spread(grdwts(:vl)*df_drhodgamma(:vl,2+spin),2,3)

        grarho_dgammadrho=grarho_dgammadrho+grarho(:vl,:3,3-spin)/partners(i_irr)* &
                          spread(grdwts(:vl)*df_drhodgamma(:vl,7-spin),2,3)

        grarho_dgdr=>grarho_dgammadrho

       elseif(ispin.eq.2) then
        grarho_drhodgamma=grarho_drhodgamma+grarho(:vl,:3,3-spin)/partners(i_irr)* & !!!(2)
                          spread(grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)
        grarho_dgdr=>grarho_drhodgamma
       else
        grarho_dgdr=>grarho_drhodgamma
       endif

         forall(k=1:3,kk=1:3)  grarho_grarho(:,k,kk)=  4.0_r8_kind/partners(i_irr)* &        !!!(4.2)
                               grarho(:vl,k,spin)*grarho(:vl,kk,spin)*df_dgammadgamma(:vl,spin)*grdwts(:vl)
          if(spin.ne.vspin) then
            forall(k=1:3,kk=1:3) grarho_grarho(:vl,k,kk)=4.0_r8_kind/partners(i_irr)* &                                    !!!(4.2)
                  grarho(:vl,k,spin)*grarho(:vl,kk,vspin)*grdwts(:vl)*df_dgammadgamma(:vl,ab)
            forall(k=1:3,kk=1:3) grarho_grarho(:vl,k,kk)=grarho_grarho(:vl,k,kk) &    !!!(4.2)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,spin)* &
                                 grarho(:vl,kk,spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)
            forall(k=1:3,kk=1:3) grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4.2)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,vspin)* &
                                 grarho(:vl,kk,vspin)*grdwts(:vl)*df_dgammadgamma(:vl,6-spin)
            forall(k=1:3,kk=1:3) grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk)+ &              !!!(4.2)
                                 grarho(:vl,k,vspin)*grarho(:vl,kk,spin)* &
                                 grdwts(:vl)*df_dgammadgamma(:vl,cc)/partners(i_irr)

          elseif(ispin.eq.2) then
            forall(k=1:3,kk=1:3) grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4.2)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,spin)* &
                                 grarho(:vl,kk,3-spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)
            forall(k=1:3,kk=1:3) grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &        !!!(4.2)
                                +2.0_r8_kind/partners(i_irr)*grarho(:vl,k,3-spin)* &
                                 grarho(:vl,kk,spin)*grdwts(:vl)*df_dgammadgamma(:vl,3+spin)

            forall(k=1:3,kk=1:3) grarho_grarho(:,k,kk)=grarho_grarho(:,k,kk) &           !!!(4.2)
                                +grarho(:vl,k,3-spin)*grarho(:vl,kk,3-spin)* &
                                 grdwts(:vl)*df_dgammadgamma(:vl,cc)/partners(i_irr)
          endif
     endif
       
      !__________________
       occ_i: do i_occ=1,occ_dim

         phi_p_fac(:,1)=wdvdrho(:)*phi_p(:vl,i_occ,l)                                      !!!(1)

        if(nl_calc) then

          phi_drdg_fac(:,:3)=grarho_drhodgamma*spread(phi_p(:vl,i_occ,l),2,3)               !!!(2)

          if(spin.eq.vspin) then
           phi_dgdr_f=>phi_drdg_fac                                                         !!! (3)
          else
           phi_dgdr_fac(:,:3)=grarho_dgammadrho(:,:3)*spread(phi_p(:vl,i_occ,l),2,3)  
           phi_dgdr_f=>phi_dgdr_fac
          endif

          graphi_adfac=spread(sum(grarho_dgdr*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2),2,3)                 !!!(3)

          phi_grarho_grarho=spread(spread(phi_p(:vl,i_occ,l),2,3),3,3)*grarho_grarho             !!!(4.2)

          forall(k=1:3) &
          dgammadgamma_fac(:,k)=sum(grarho_grarho(:,:,k)*graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2) !!!(4.1)



          dgrarho_fac(:,:3,uord)=spread(graphi_dgrarho_fac(:,uord),2,3)* &                                 !!! (5.1)
                                  graphi(i_irr,spin)%m(:vl,:3,i_occ,l)
          phi_wdgrarho(:,uord)=graphi_dgrarho_fac(:,uord)*phi_p(:vl,i_occ,l)                               !!! (5.2)

          if(spin.ne.vspin) then
           dgrarho_fac(:,:3,ud)=spread(graphi_dgrarho_fac(:,ud),2,3)* &                                    !!! (5.1)
                                graphi(i_irr,spin)%m(:vl,:3,i_occ,l)
           phi_wdgrarho(:,ud)=graphi_dgrarho_fac(:,ud)*phi_p(:vl,i_occ,l)                                  !!! (5.2)
          endif

          phi_grarho_dgrarho(:vl,:3)=grarho_dgrarho(:vl,:3)* &                                             !!! (7.1)
                                     spread(phi_p(:vl,i_occ,l),2,3)
          graphi_grarho_dgrarho(:)=sum(grarho_dgrarho(:,:)* &                                              !!! (7.2)
                                       graphi(i_irr,spin)%m(:vl,:3,i_occ,l),2)
         endif

    !_________________________
     q_vir: do j_vir=1,vir_dim

         nucgra_afac(:vl,:3)=spread(phi_p_fac(:vl,1)*phi_p(:vl,j_vir+j_off,l),2,3)                !!! (1)

     if(nl_calc) then
         nucgra_afac=nucgra_afac+graphi_adfac(:vl,:3)*spread(phi_p(:vl,j_vir+j_off,l),2,3)        !!! (3.1)
         nucgra_afac(:vl,:3)=nucgra_afac(:vl,:3)+ &                                               !!! (3.2) dgdr
          spread(sum(phi_dgdr_f(:,:3)*graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l),2),2,3)
      endif
                          ! ________________________
      call  nuc_dervs_afac('q_calc',vspin,nuc_grarho)

       if(nl_calc) then
            nucgra_afac(:vl,:3)= spread(phi_p(:vl,j_vir+j_off,l),2,3)*( &
              phi_drdg_fac(:vl,:3)+ &                                                           !!! (2)
              dgammadgamma_fac(:vl,:3)+ &                                                       !!! (4) 
              dgrarho_fac(:vl,:3,uord+abs(spin-vspin)) )                                        !!! (5)

            nucgra_afac=nucgra_afac+spread(phi_wdgrarho(:,uord+abs(spin-vspin)),2,3)* &         !!! (5)
                                           graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l)

            forall(k=1:3,kk=1:3) &
            nuc2der_afac(:,k,kk)=nucgra_afac(:,kk)+sum(graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l)* &
                                              phi_grarho_grarho(:,:3,kk),2)
                                !____________________________________________
            call  nuc_dervs_afac('q_calc',vspin,nuc_sec_derrho=nuc_sec_derrho)

       endif

      if(present(graw)) then

         graw_adfac=spread(dfdrho(:vl,spin)*phi_p(:vl,i_occ,l)* &                        !!! (6)
                            phi_p(:vl,j_vir+j_off,l),2,3)/partners(i_irr)
        if(nl_calc)  then
         graw_adfac=graw_adfac+ spread(&                                                  !!! (7.1)
          sum(phi_grarho_dgrarho*graphi(i_irr,spin)%m(:vl,:3,j_vir+j_off,l),2),2,3)
         graw_adfac=graw_adfac+spread(graphi_grarho_dgrarho*phi_p(:vl,j_vir+j_off,l),2,3) !!! (7.2)
        endif
                            !________________________
        call  nuc_dervs_afac('q_calc',vspin,graw=graw)  !!! (6+7)
      endif

       enddo q_vir

       h1_occ_j: do j_occ=1,occ_dim
    
         nucgra_afac(:vl,:3)=spread(phi_p_fac(:vl,1)*phi_p(:vl,j_occ,l),2,3)         !!! (1)
                                !_________________________
            call  nuc_dervs_afac('h1calc',vspin,nuc_grarho)

       if(nl_calc) then
         nucgra_afac=spread(phi_p(:vl,j_occ,l),2,3)*phi_drdg_fac                     !!! (2,h1)
         nuc2der_afac(:vec_length,:3,:3)=spread(nucgra_afac,2,3)
                             !____________________________________________
         call  nuc_dervs_afac('h1calc',vspin,nuc_sec_derrho=nuc_sec_derrho)
        endif 


        if(present(graw)) then
         graw_adfac(:vl,:3)=spread(dfdrho(:vl,spin) &                                 !!! (6)
           *phi_p(:vl,i_occ,l)*phi_p(:vl,j_occ,l),2,3)/partners(i_irr)
          
         call  nuc_dervs_afac('h1calc',vspin,graw=graw)
        endif

   enddo h1_occ_j 

    if(nl_calc) then
       graphi_h1_occ: do j_occ=1,occ_dim
        nucgra_afac= graphi_adfac*spread(phi_p(:vl,j_occ,l),2,3)                         !!! (3)
        call  nuc_dervs_afac('h1gphi',vspin,nuc_grarho)

        nuc2der_afac(:vl,:3,:3)=spread(spread(phi_p(:vl,j_occ,l),2,3)*&                   !!! (4.2)
                                       dgammadgamma_fac(:,:3),2,3)                                            
        nuc2der_afac=nuc2der_afac+spread(spread(phi_p(:vl,j_occ,l),2,3) * &               !!! (5)
                                         dgrarho_fac(:,:3,uord+abs(spin-vspin)),2,3)
        call  nuc_dervs_afac('h1gphi',vspin,nuc_sec_derrho=nuc_sec_derrho)

       if(present(graw)) then
        graw_adfac(:vl,:3)= spread( sum(phi_grarho_dgrarho(:vl,:3)* &                     !!! (7)
                                        graphi(i_irr,spin)%m(:vl,:3,j_occ,l),2),2,3)
        call  nuc_dervs_afac('h1gphi',vspin,graw=graw)
       endif

       enddo graphi_h1_occ
    endif

 enddo occ_i

 end subroutine h1q_func_dervs
#endif

 subroutine v_lambda_imp(Qai,dvdrho,cpks_spin)

       !!! as now in each call  only up or down nuc_imp_grarho(spin) and nuc_imp_gragamma(spin) are known
       !!! thus cpks(cpks_spin=up) and cpks(cpks_spin=dn) dependent on these nuc_imp_grarho nuc_imp_gragamma 
       !!! are calculated 

       integer(kind=i4_kind), intent(in),optional:: cpks_spin
       real(kind=r8_kind),intent(in):: dvdrho(:)
       real(kind=r8_kind), intent(inout):: Qai(:,:)
       integer(kind=i4_kind):: occ_dim,vir_dim,vl
#if 0
       real(kind=r8_kind):: phi_p_adfac(size(Qai,1),vec_length),h1_phi_p_adfac(size(Qai,1),vec_length)
       real(kind=r8_kind):: wdvdrho(vec_length)
       occ_dim=size(Qai,1)
       vir_dim=size(Qai,2)
       vl=vec_length
#else
       integer(i4_kind):: stat
       real(kind=r8_kind),allocatable:: phi_p_adfac(:,:),h1_phi_p_adfac(:,:)
       real(kind=r8_kind),allocatable:: wdvdrho(:)
       occ_dim=size(Qai,1)
       vir_dim=size(Qai,2)
       vl=vec_length
       allocate(phi_p_adfac(occ_dim,vl),h1_phi_p_adfac(occ_dim,vl),wdvdrho(vl),stat=stat)
       ASSERT(stat.eq.0)
#endif


       phi_p=>phi_ob(i_irr)%o(:,:,:,cpks_spin)     ! occ. orbitals for this spin

       wdvdrho=grdwts(:vl)*dvdrho(:vl)
       dvdrho_fac(:vl)=nuc_imp_grarho(:vl,i_grad)*wdvdrho/partners(i_irr)        !!! (1)
      
       nl: if(nl_calc)  then
                                !!!contribs:__(4)___________(5)____________(2)_______
        call nuc_imp_gragamma_facs(cpks_spin,dgammadgamma_fac,dgrarho_fac,drhodgamma_fac) 


        if(spin.ne.cpks_spin) then   !! dgdr (3)
         graphi_adfac(:vl,:3)= two*grarho(:vl,:3,3-spin)/partners(i_irr) & !!!(3)
            *spread(nuc_imp_grarho(:vl,i_grad)*grdwts(:vl)*df_drhodgamma(:vl,5-spin),2,3)
         graphi_adfac(:vl,:3)=graphi_adfac(:vl,:3)+ &                      !!!(3)
         grarho(:vl,:3,spin)/partners(i_irr) &
            *spread(nuc_imp_grarho(:vl,i_grad)*grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3)
        else
         graphi_adfac(:vl,:3)=two*grarho(:vl,:3,spin)/partners(i_irr) &   !!!(3)dgdr
            *spread(nuc_imp_grarho(:vl,i_grad)*grdwts(:vl)*df_drhodgamma(:vl,spin),2,3) 
         if(ispin.eq.2) &
         graphi_adfac(:vl,:3)=graphi_adfac(:vl,:3)+ &                     !!!(3)
         grarho(:vl,:3,3-spin)/partners(i_irr) &
            *spread(nuc_imp_grarho(:vl,i_grad)*grdwts(:vl)*df_drhodgamma(:vl,4+spin),2,3) 
        endif

       endif nl


     do  i_occ=1,occ_dim

       phi_p_adfac(i_occ,:)=phi_p(:vl,i_occ,l)*dvdrho_fac     !!! (1)
 

        if(nl_calc) then
                                   !______(5)________ !______(2)________  ! _________(3)____
         call graphi_facs(cpks_spin,graphi_dgrarho_fac,phi_adfac=phi_adfac,graphi_fac=graphi_fac, &
                                    graphi_dgdg_fac=graphi_dgdg_fac) !!! (4)

         phi_p_adfac(i_occ,:)=phi_p_adfac(i_occ,:)+phi_adfac+graphi_fac+graphi_dgdg_fac &
                           +graphi_dgrarho_fac(:,1+abs(spin-cpks_spin))
         h1_phi_p_adfac(i_occ,:)=graphi_fac+graphi_dgdg_fac+graphi_dgrarho_fac(:,1+abs(spin-cpks_spin))
        endif

     enddo !occ

          call dgemm('n', 'n', occ_dim, vir_dim, vl, &
                   one, phi_p_adfac, occ_dim, phi_p(1, 1+eig_dim-vir_dim, l), full_vl, &
                   one,Qai(:,:), occ_dim)

    if(i_cpks.gt.n_cpks.or.Q_calc) then
          call dgemm('n','n',occ_dim,occ_dim,vl, &
                   one,phi_p_adfac,occ_dim, phi_p(:,:,l),full_vl, &
                   one,cpks(i_grad,i_irr,cpks_spin)%h1,occ_dim)
          if(nl_calc) &
          call dgemm('t','t',occ_dim,occ_dim,vl, &
                   one,phi_p(:,:,l),full_vl, h1_phi_p_adfac,occ_dim, &
                   one,cpks(i_grad,i_irr,cpks_spin)%h1,occ_dim)
    endif

    if(nl_calc) call nl_calcQai_grphi_vir(Qai,cpks_spin,graphi_adfac=graphi_adfac,  &    !!! (3)  n_cpks+1
                               dgrarho_fac=dgrarho_fac(:,:,1+abs(spin-cpks_spin)),  &    !!! (5)
                               dgammadgamma_fac=dgammadgamma_fac   )                     !!! (4)

#if 1
       deallocate(phi_p_adfac,h1_phi_p_adfac,wdvdrho,stat=stat)
       ASSERT(stat.eq.0)
#endif
  end subroutine v_lambda_imp

              subroutine calc_help_sd_eq(i_ua,j,help_sd_eq,ghelp_nuc_eq,help_sec_der_arr,help_nuc_arr)
              integer(kind=i4_kind),intent(in):: j,i_ua
              real(kind=r8_kind),intent(in):: help_sec_der_arr(:,:,:,:)
              real(kind=r8_kind),intent(in):: help_nuc_arr(:,:,:,:)
              real(kind=r8_kind), intent(out):: ghelp_nuc_eq(:,:,:,:),help_sd_eq(:,:,:)
              integer(kind=i4_kind)::k

                  do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                    help_sd(:vec_length,1,1)=help_sec_der_arr(:vec_length,1,i_ea,j) &
                                             *wdfdgrarho(:vec_length)

                    help_sd(:vec_length,1,2)=help_sec_der_arr(:vec_length,2,i_ea,j) &
                                             *wdfdgrarho(:vec_length)
                    help_sd(:vec_length,2,1)=help_sd(:vec_length,1,2)

                    help_sd(:vec_length,1,3)=help_sec_der_arr(:vec_length,3,i_ea,j) &
                                             *wdfdgrarho(:vec_length)
                    help_sd(:vec_length,3,1)=help_sd(:vec_length,1,3)

                    help_sd(:vec_length,2,2)=help_sec_der_arr(:vec_length,4,i_ea,j) &
                                             *wdfdgrarho(:vec_length)

                    help_sd(:vec_length,2,3)=help_sec_der_arr(:vec_length,5,i_ea,j) &
                                             *wdfdgrarho(:vec_length)
                    help_sd(:vec_length,3,2)=help_sd(:vec_length,2,3)

                    help_sd(:vec_length,3,3)=help_sec_der_arr(:vec_length,6,i_ea,j) &
                                             *wdfdgrarho(:vec_length)

                    help_sd_eq(:vec_length,:3,i_ea)=two/partners(i_irr)* &
                        sum(help_sd(:vec_length,:3,:3) &
                                *spread(grarho(:vec_length,:3,spin),2,3),3)
                   if(ispin.eq.2) then
                    help_sd(:vec_length,1,1)=help_sec_der_arr(:vec_length,1,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)

                    help_sd(:vec_length,1,2)=help_sec_der_arr(:vec_length,2,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)
                    help_sd(:vec_length,2,1)=help_sd(:vec_length,1,2)

                    help_sd(:vec_length,1,3)=help_sec_der_arr(:vec_length,3,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)
                    help_sd(:vec_length,3,1)=help_sd(:vec_length,1,3)

                    help_sd(:vec_length,2,2)=help_sec_der_arr(:vec_length,4,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)

                    help_sd(:vec_length,2,3)=help_sec_der_arr(:vec_length,5,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)
                    help_sd(:vec_length,3,2)=help_sd(:vec_length,2,3)

                    help_sd(:vec_length,3,3)=help_sec_der_arr(:vec_length,6,i_ea,j) &
                                      *grdwts(:vec_length)*dfdgrarho(:vec_length,updn)

                    help_sd_eq(:vec_length,:3,i_ea)=help_sd_eq(:vec_length,:3,i_ea)+ &
                        sum(help_sd(:vec_length,:3,:3) &
                         *spread(grarho(:vec_length,:3,3-spin),2,3),3)/partners(i_irr)
                   endif

                    forall(k=1:3)&
                    ghelp_nuc_eq(:vec_length,:3,k,i_ea)=two/partners(i_irr)* &
                    help_nuc_arr(:vec_length,:3,i_ea,j)*spread(wdfdgrarho(:vec_length)* &
                                                        grarho(:vec_length,k,spin),2,3)
                    if(ispin.eq.2) forall(k=1:3)&
                    ghelp_nuc_eq(:vec_length,:3,k,i_ea)=ghelp_nuc_eq(:vec_length,:3,k,i_ea) +  &
                    help_nuc_arr(:vec_length,:3,i_ea,j)*spread(grdwts(:vec_length)* &
                     dfdgrarho(:vec_length,updn)*grarho(:vec_length,k,3-spin),2,3)/partners(i_irr)
                    

                 enddo
              end subroutine calc_help_sd_eq


             function dfdgrarho_grads2(i,ghelp_nuc_eq) result(help_gp)
                !!! used to calc h1
                integer(kind=i4_kind), intent(in):: i
                real(kind=r8_kind),intent(in):: ghelp_nuc_eq(:,:,:,:)
                real(kind=r8_kind):: help_gp(3)

               forall(k=1:3) &
               help_gp(k)=-sum(sum(graphi(i_irr,spin)%m(:vec_length,:3,i,l)* &
                               ghelp_nuc_eq(:vec_length,k,:3,i_ea),2),1)
                    
              end function dfdgrarho_grads2


                function dfdgrarho_grads3(j,help_sec_der_arr) result(help_gp)
                real(kind=r8_kind),intent(in):: help_sec_der_arr(:,:,:,:)
                integer(kind=i4_kind), intent(in):: j
                real(kind=r8_kind):: help_gp(3)
                real(kind=r8_kind):: help_sd(vec_length,3,3)

                    help_sd(:vec_length,1,1)=help_sec_der_arr(:vec_length,1,i_ea,j)

                    help_sd(:vec_length,1,2)=help_sec_der_arr(:vec_length,2,i_ea,j)
                    help_sd(:vec_length,2,1)=help_sd(:vec_length,1,2)

                    help_sd(:vec_length,1,3)=help_sec_der_arr(:vec_length,3,i_ea,j)
                    help_sd(:vec_length,3,1)=help_sd(:vec_length,1,3)

                    help_sd(:vec_length,2,2)=help_sec_der_arr(:vec_length,4,i_ea,j)

                    help_sd(:vec_length,2,3)=help_sec_der_arr(:vec_length,5,i_ea,j)
                    help_sd(:vec_length,3,2)=help_sd(:vec_length,2,3)

                    help_sd(:vec_length,3,3)=help_sec_der_arr(:vec_length,6,i_ea,j) 

              help_gp(:)=-sum(sum(help_sd(:vec_length,:3,:3)*nuc2der_afac(:vec_length,:3,:3),3),1)
                    
        end function dfdgrarho_grads3

        function dfdgrarho_grads4(j,k,i_ea,help_sec_der_arr) result(help_gp)

                real(kind=r8_kind),intent(in):: help_sec_der_arr(:,:,:,:)
                integer(kind=i4_kind), intent(in):: j,k,i_ea
                real(kind=r8_kind):: help_gp(vec_length)
                real(kind=r8_kind):: help_sd(vec_length,3,3)
                integer(kind=i4_kind):: kk(3,3)=reshape( (/1,2,3, 2,4,5, 3,5,6/),(/3,3/) )

                    help_sd(:vec_length,k,1)=help_sec_der_arr(:vec_length,kk(1,k),i_ea,j) 
                    help_sd(:vec_length,k,2)=help_sec_der_arr(:vec_length,kk(2,k),i_ea,j)
                    help_sd(:vec_length,k,3)=help_sec_der_arr(:vec_length,kk(3,k),i_ea,j)

              help_gp(:)=sum(help_sd(:vec_length,k,:3)*grarho_dgrarho(:vec_length,:3),2)* &
                         grdwts(:vec_length)
                    
        end function dfdgrarho_grads4

         subroutine dealloc_dervs_helps()
                MEMLOG(-size(help_nuc_arr)-size(nucgra_efac)-size(help_sd_eq)-size(ghelp_nuc_eq))
                deallocate(help_nuc_arr,nucgra_efac,help_sd_eq,&
                           ghelp_nuc_eq,stat=cpksalloc(94))
                ASSERT(cpksalloc(94).eq.0)
                cpksalloc(94)=1


          if(nl_calc) then
                MEMLOG(-size(help_sec_der_arr))
                deallocate(help_sec_der_arr,stat=alloc_stat(21))
                ASSERT(alloc_stat(21).eq.0)
                alloc_stat(21)=1
          endif
         end subroutine dealloc_dervs_helps

        subroutine alloc_dervs_helps(mo_dim)
        integer(kind=i4_kind), intent(in):: mo_dim

         allocate(help_nuc_arr(vec_length,3,unique_atoms(i_ua)%n_equal_atoms,mo_dim),&
                  nucgra_efac(vec_length,3,unique_atoms(i_ua)%n_equal_atoms), &
                  help_sd_eq(vec_length,3,unique_atoms(i_ua)%n_equal_atoms), &
                  ghelp_nuc_eq(vec_length,3,3,unique_atoms(i_ua)%n_equal_atoms), &
                  stat=cpksalloc(94))
         ASSERT(cpksalloc(94).eq.0)
                MEMLOG(size(help_nuc_arr)+size(nucgra_efac)+size(ghelp_nuc_eq)+size(help_sd_eq))


        if(nl_calc) then ! compare with b part
          allocate(help_sec_der_arr(vec_length,6,unique_atoms(i_ua)%n_equal_atoms,mo_dim),&
                  stat=alloc_stat(21))
          ASSERT(alloc_stat(21).eq.0)
          MEMLOG(size(help_sec_der_arr))

       endif

       end subroutine alloc_dervs_helps

         subroutine calc_dervs_helps(eigv_offset,help_sec_der_arr)
         real(kind=r8_kind),intent(out),optional :: help_sec_der_arr(:,:,:,:)
         integer(kind=i4_kind), optional, intent(in):: eigv_offset
         integer(kind=i4_kind):: off

         if(present(eigv_offset)) then
          off=eigv_offset
         else
          off=0
         endif

       ! calculate ( d/dN Phi_occ)
                help_arr: if(vec_length==full_vl) then
                   call dgemm('n','n',vec_length*3*unique_atoms(i_ua)%n_equal_atoms, &
                        size(help_nuc_arr,4),orb_dim, &
                        1.0_r8_kind, &
                        nuc_grads(i_ua,i_irr)%o(:,:,:,:,l),vec_length*3*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1+off,spin),eig_dim, &
                        0.0_r8_kind,help_nuc_arr(1,1,1,1), &
                        vec_length*3*unique_atoms(i_ua)%n_equal_atoms)

!                 print*,'Q help_nuc_arr',sum(help_nuc_arr),cpks_test_counter,i_cpks

                 if(present(help_sec_der_arr)) then
                    call dgemm('n','n',vec_length*6*unique_atoms(i_ua)%n_equal_atoms,&
                        size(help_sec_der_arr,4),orb_dim, &
                        1.0_r8_kind,&
                        nuc_sec_der(i_ua,i_irr)%o(1,1,1,1,l),vec_length*6*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1+off,spin),eig_dim, &
                        0.0_r8_kind,help_sec_der_arr(1,1,1,1),&
                        vec_length*6*unique_atoms(i_ua)%n_equal_atoms)
                 endif

                else help_arr
                   ! copy data first
                   allocate(nuc_help(vec_length,3,unique_atoms(i_ua)%n_equal_atoms,orb_dim), &
                            stat=cpksalloc(95))
                   ASSERT(cpksalloc(95).eq.0)

                   nuc_help(:vec_length,:,:,:)=nuc_grads(i_ua,i_irr)%o(:vec_length,:,:,:,l)

                   call dgemm('n','n',vec_length*3*unique_atoms(i_ua)%n_equal_atoms,&
                                      size(help_nuc_arr,4),orb_dim, &
                        1.0_r8_kind,&
                        nuc_help(1,1,1,1),vec_length*3*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1+off,spin),eig_dim, &
                        0.0_r8_kind,help_nuc_arr(1,1,1,1),&
                                    vec_length*3*unique_atoms(i_ua)%n_equal_atoms)

                   deallocate(nuc_help, stat=cpksalloc(95))
                   ASSERT(cpksalloc(95).eq.0)
                   cpksalloc(95)=1

                 if(present(help_sec_der_arr)) then
                  allocate(sec_der_help(vec_length,6,unique_atoms(i_ua)%n_equal_atoms,orb_dim),&
                            stat=alloc_stat(22))
                    ASSERT(alloc_stat(22).eq.0)
                    MEMLOG(size(sec_der_help))
                 
                  sec_der_help(:vec_length,:,:,:)=nuc_sec_der(i_ua,i_irr)%o(:vec_length,:,:,:,l)

                   call dgemm('n','n',vec_length*6*unique_atoms(i_ua)%n_equal_atoms,&
                        size(help_sec_der_arr,4),orb_dim,1.0_r8_kind,&
                        sec_der_help(1,1,1,1),vec_length*6*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1+off,spin),eig_dim, &
                        0.0_r8_kind, help_sec_der_arr(1,1,1,1),&                !to calc nuc_dervsrho
                        vec_length*6*unique_atoms(i_ua)%n_equal_atoms)


                    MEMLOG(-size(sec_der_help))
                    deallocate(sec_der_help, stat=alloc_stat(22))
                    ASSERT(alloc_stat(22).eq.0)
                    alloc_stat(22)=1

                 endif
                endif  help_arr
       ! ** ( d/dN Phi_occ) calculated
       end subroutine calc_dervs_helps

       subroutine graphi_facs(cpks_spin,graphi_dgrarho_fac,phi_adfac,graphi_fac,graphi_dgdg_fac)
       integer(kind=i4_kind),intent(in):: cpks_spin
       real(kind=r8_kind), intent(out):: graphi_dgrarho_fac(:,:)
       real(kind=r8_kind), intent(out),optional:: phi_adfac(:),graphi_fac(:),graphi_dgdg_fac(:)
       integer(kind=i4_kind)::vl

       vl=vec_length

        phi_adfac=drhodgamma_fac(:vl) *phi_p(:vl,i_occ,l)                     !!! (2)

        graphi_fac(:vl)=sum(graphi(i_irr,cpks_spin)%m(:vl,:3,i_occ,l)* &      !!! (3)
                            graphi_adfac(:vl,:3),2) 

        graphi_dgdg_fac(:vl)=sum(graphi(i_irr,cpks_spin)%m(:vl,:3,i_occ,l)* & !!! (4)
                                 dgammadgamma_fac(:vl,:3),2) 

       if(spin.eq.cpks_spin)  then
         graphi_dgrarho_fac(:vl,uord)= &                                      !!! (5)
           sum(graphi(i_irr,cpks_spin)%m(:vl,:3,i_occ,l)*dgrarho_fac(:vl,:3,uord),2)
        else
         graphi_dgrarho_fac(:vl,ud)= &                                        !!! (5)
           sum(graphi(i_irr,cpks_spin)%m(:vl,:3,i_occ,l)*dgrarho_fac(:vl,:3,ud),2)
        endif

       end subroutine graphi_facs

       subroutine nuc_imp_gragamma_facs(cpks_spin,dgammadgamma_fac,dgrarho_fac,drhodgamma_fac)
       
       integer(kind=i4_kind),intent(in) :: cpks_spin
       real(kind=r8_kind), intent(out):: dgammadgamma_fac(:,:)
       real(kind=r8_kind), intent(out):: dgrarho_fac(:,:,:)
       real(kind=r8_kind), intent(out),optional:: drhodgamma_fac(:)
       
       integer(kind=i4_kind):: vl
       vl=vec_length

       if(spin.eq.cpks_spin) then
        drhodgamma_fac(:vec_length)= two*sum( &                           !!! (2)  drdg
               nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2) &
                   *grdwts(:vl)*df_drhodgamma(:vl,spin)/partners(i_irr)
        if(ispin.eq.2) drhodgamma_fac(:vl)=drhodgamma_fac(:vl)+sum( &
               nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,3-spin),2) &
                    *grdwts(:vl)*df_drhodgamma(:vl,4+spin)/partners(i_irr)

        dgammadgamma_fac=4.0_r8_kind/partners(i_irr)* &          !!! (4)
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2)* &
                grdwts(:vl)*df_dgammadgamma(:vl,spin),2,3)*grarho(:vl,:3,spin)

       if(ispin.eq.2) then
        dgammadgamma_fac=dgammadgamma_fac+2.0_r8_kind/partners(i_irr)* &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,3-spin),2)* &
            grdwts(:vl)*df_dgammadgamma(:vl,3+spin),2,3)*grarho(:vl,:3,spin)
        dgammadgamma_fac(:vl,:3)=dgammadgamma_fac(:vl,:3)+2.0_r8_kind/partners(i_irr)* &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2)* &
            grdwts(:vl)*df_dgammadgamma(:vl,3+spin),2,3)*grarho(:vl,:3,3-spin)

        dgammadgamma_fac=dgammadgamma_fac+ &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,3-spin),2)* &
         grdwts(:vl)*df_dgammadgamma(:vl,cc),2,3)*grarho(:vl,:3,3-spin)/partners(i_irr)
       endif

        dgrarho_fac(:vl,:3,uord)=two/partners(i_irr)* &                   !!! (5)
          nuc_imp_gragamma(:vl,:3,i_grad)*spread(grdwts(:vl)*dfdgrarho(:vl,spin),2,3)

       elseif(spin.ne.cpks_spin) then
        drhodgamma_fac=two*sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2) &    !!! (2)
                      *grdwts(:vl)*df_drhodgamma(:vl,2+spin)/partners(i_irr) !!! dup=4
        drhodgamma_fac=drhodgamma_fac+sum(nuc_imp_gragamma(:vl,:3,i_grad)* &
           grarho(:vl,:3,3-spin),2)*grdwts(:vl)*df_drhodgamma(:vl,7-spin)/partners(i_irr)

        dgrarho_fac(:vl,:3,ud)= &                     !!! (5)
          nuc_imp_gragamma(:vec_length,:3,i_grad)*spread(grdwts(:vl)* &
             dfdgrarho(:vl,updn),2,3)/partners(i_irr)

        dgammadgamma_fac=4.0_r8_kind/partners(i_irr)* &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2)* &
            grdwts(:vl)*df_dgammadgamma(:vl,ab),2,3)*grarho(:vl,:3,3-spin)

        dgammadgamma_fac=dgammadgamma_fac+2.0_r8_kind/partners(i_irr)* &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,3-spin),2)* &
            grdwts(:vl)*df_dgammadgamma(:vl,6-spin),2,3)*grarho(:vl,:3,3-spin)

        dgammadgamma_fac=dgammadgamma_fac+2.0_r8_kind/partners(i_irr)* &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,spin),2)* &
            grdwts(:vl)*df_dgammadgamma(:vl,3+spin),2,3)*grarho(:vl,:3,spin)

        dgammadgamma_fac=dgammadgamma_fac+ &
         spread(sum(nuc_imp_gragamma(:vl,:3,i_grad)*grarho(:vl,:3,3-spin),2)* &
         grdwts(:vl)*df_dgammadgamma(:vl,cc),2,3)*grarho(:vl,:3,spin)/partners(i_irr)
       endif

   end subroutine nuc_imp_gragamma_facs

   subroutine nl_calcQai_grphi_vir(Qai,cpks_spin,graphi_adfac,dgrarho_fac,dgammadgamma_fac)

       real(kind=r8_kind), intent(inout) :: Qai(:,:)
       integer(kind=i4_kind), intent(in) :: cpks_spin
       real(kind=r8_kind), intent(in),optional::graphi_adfac(:,:), &
                                                dgrarho_fac(:,:),&     !!! (5) uord / ud 
                                                dgammadgamma_fac(:,:) 
       
       real(kind=r8_kind):: graphi_fac(vec_length,size(Qai,2)),fac(vec_length,3)
       integer(kind=i4_kind):: off,occ_dim,vl,vir_dim
      
       occ_dim=size(Qai,1)
       vir_dim=size(Qai,2)
       off=eig_dim-vir_dim

       vl=vec_length
       fac=graphi_adfac+dgrarho_fac+dgammadgamma_fac !!! (3+5+4)
       do  j_vir=1,vir_dim
        graphi_fac(:,j_vir)=sum(graphi(i_irr,cpks_spin)%m(:vl,:3,j_vir+off,l)*fac,2)        !!! (3) df_dgdr contrib
!        Qai(:,j_vir)=Qai(:,j_vir)+sum(spread(graphi_fac,2,occ_dim)*phi_p(:vl,:occ_dim,l),1)
       enddo
          call dgemm('t','n',occ_dim,vir_dim,vl, &
                   one,phi_p(:,:,l),full_vl, graphi_fac,vl, &
                   one,Qai,occ_dim)

  end subroutine nl_calcQai_grphi_vir


  subroutine cpks_average_gradient(i_ua,grad_final)
    ! purpose : builds final gradient by averaging over all
    !           equal atoms. this is necesarry because the
    !           integration grid and the gradients are not totalsymmetric
    !
    integer,intent(in) :: i_ua
    type(arrmat2),intent(inout) :: grad_final(:)
    !** End of interface *****************************************

    real(kind=r8_kind) :: rotmat_tot(3,3),help_arr_local(3,n_equal_max)
    real(kind=r8_kind),pointer :: rotmat1(:,:),rotmat2(:,:)
    integer :: i_equal1,i_equal2,kk,n_ea

       n_ea = unique_atoms(moving_unique_atom_index(i_ua))%n_equal_atoms
       help_arr_local(:,1:n_ea)=0.0_r8_kind
       do i_equal1=1,n_ea
          do i_equal2=1,n_ea
             rotmat1=>unique_atom_grad_info(i_ua)%m(:,:,i_equal1)
             rotmat2=>unique_atom_grad_info(i_ua)%m(:,:,i_equal2)
             rotmat_tot=matmul(transpose(rotmat1),rotmat2)
             do kk=1,3
                help_arr_local(kk,i_equal1)=&
                     help_arr_local(kk,i_equal1)+&
                     sum(grad_final(i_ua)%m(:,i_equal2)*rotmat_tot(kk,:))
             end do
          end do
       end do
       grad_final(i_ua)%m=help_arr_local(:,:n_ea)
  end subroutine cpks_average_gradient

  subroutine cpks_transform_to_totsymQai(i_ua,mo_occ,mo_vir,spin)

    integer(kind=i4_kind),intent(in) :: i_ua,mo_occ,mo_vir,spin

             index = cpks_gradient_index(i_ua) - 1
             grad_dim = cpks_gradient_index(i_ua+1) - cpks_gradient_index(i_ua)
                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                if(grad_dim==3) then
                 if(sum((unique_atom_grad_info(i_ua)%m(:,:,i_ea)-unity_matrix)**2)<&
                    1.0e-7_r8_kind) then
                    do_rotation=.false.
                 else
                    do_rotation=.true.
                 endif
                else
                    do_rotation=.true.
                end if


                rot: if(do_rotation) then
                 gdim: do i_grad=1,grad_dim
                        cpks(index+i_grad,i,spin)%Qai(mo_occ,mo_vir)=cpks(index+i_grad,i,spin)%Qai(mo_occ,mo_vir) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,1,i_ea)*cpks_grad_xc(i_ua)%m(1,i_ea) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,2,i_ea)*cpks_grad_xc(i_ua)%m(2,i_ea) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,3,i_ea)*cpks_grad_xc(i_ua)%m(3,i_ea)
                 enddo gdim

                else rot
                Print*, ' NO ROTATION CASE IN CPKS DENSYTY MODULE PART'
                        cpks(index,  i,spin)%Qai(mo_occ,mo_vir)=cpks(index,  i,spin)%Qai(mo_occ,mo_vir) &
                                                          +cpks_grad_xc(i_ua)%m(1,i_ea) 
                        cpks(index+1,i,spin)%Qai(mo_occ,mo_vir)=cpks(index+1,i,spin)%Qai(mo_occ,mo_vir) &
                                                          +cpks_grad_xc(i_ua)%m(2,i_ea) 
                        cpks(index+2,i,spin)%Qai(mo_occ,mo_vir)=cpks(index+2,i,spin)%Qai(mo_occ,mo_vir) &
                                                          +cpks_grad_xc(i_ua)%m(3,i_ea) 
                 endif rot
               enddo
  end subroutine cpks_transform_to_totsymQai

  subroutine h1_transform_to_totsym(i_ua,occ1,occ2,spin)
    implicit none
    integer(kind=i4_kind),intent(in) :: i_ua,occ1,occ2,spin
    ! *** end of interface ***

             index = cpks_gradient_index(i_ua)-1
             grad_dim = cpks_gradient_index(i_ua+1) - cpks_gradient_index(i_ua)

                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

                if(grad_dim==3) then
                 if(sum((unique_atom_grad_info(i_ua)%m(:,:,i_ea)-unity_matrix)**2)<&
                    1.0e-7_r8_kind) then
                    do_rotation=.false.
                 else
                    do_rotation=.true.
                 endif
                else
                    do_rotation=.true.
                end if


                rot: if(do_rotation) then
                 gdim: do i_grad=1,grad_dim
                        cpks(index+i_grad,i_irr,spin)%h1(occ1,occ2)=cpks(index+i_grad,i_irr,spin)%h1(occ1,occ2) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,1,i_ea)*cpks_grad_xc(i_ua)%m(1,i_ea) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,2,i_ea)*cpks_grad_xc(i_ua)%m(2,i_ea) &
                        +unique_atom_grad_info(i_ua)%m(i_grad,3,i_ea)*cpks_grad_xc(i_ua)%m(3,i_ea)
                 enddo gdim

                else rot
                Print*, ' NO ROTATION CASE IN CPKS DENSYTY MODULE PART'
                STOP
                        cpks(index,  i,spin)%h1(occ1,occ2)=cpks(index,  i,spin)%h1(occ1,occ2) &
                                                          +cpks_grad_xc(abs(i_ua))%m(1,i_ea) 
                        cpks(index+1,i,spin)%h1(occ1,occ2)=cpks(index+1,i,spin)%h1(occ1,occ2) &
                                                          +cpks_grad_xc(abs(i_ua))%m(2,i_ea) 
                        cpks(index+2,i,spin)%h1(occ1,occ2)=cpks(index+2,i,spin)%h1(occ1,occ2) &
                                                          +cpks_grad_xc(abs(i_ua))%m(3,i_ea) 
                 endif rot
               enddo

    end subroutine h1_transform_to_totsym

  end subroutine cpksdervs_xc
#endif


#ifdef WITH_SECDER
 subroutine transform_to_cart_grarho(grad_cart,grad_totsym,spin)
   ! purpose: transform symmetry adapted gradient components to cartesian
   !          coordinates and add them to an array of cartesian gradients
 
  use cpksdervs_matrices, only:cpks_gradient_index
  integer(kind=i4_kind), intent(in):: spin
  type(arrmat4)     , intent(inout) :: grad_cart(:)   ! cartesian grad. array
  real(kind=r8_kind), intent(in   ) :: grad_totsym(:,:) ! symm. adapt. contr.
  !
  integer(kind=i4_kind) :: i_unique,i_center,i_equal,index,i,grad_dim,vl
  integer(kind=i4_kind) :: k
 
  real(kind=r8_kind),pointer :: rotmat(:,:)
  vl=size(grad_totsym,1)
 
  do i_unique=1,N_moving_unique_atoms
     i_center=moving_unique_atom_index(i_unique)
     index=cpks_gradient_index(i_unique)
     grad_dim=cpks_gradient_index(i_unique+1)-index
     index=index-1
     do i_equal=1,unique_atoms(i_center)%n_equal_atoms
        rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
        do i=1,grad_dim
        do k=1,3
           grad_cart(i_unique)%m(:vl,k,i_equal,spin) = &
           grad_cart(i_unique)%m(:vl,k,i_equal,spin) +  rotmat(i,k)*grad_totsym(:,index+i)
        enddo
        enddo
     end do
  enddo
 end subroutine transform_to_cart_grarho

 subroutine transform_to_cart_gragamma(grad_cart,grad_totsym,spin)
   ! purpose: transform symmetry adapted gradient components to cartesian
   !          coordinates and add them to an array of cartesian gradients
 
  use cpksdervs_matrices, only:cpks_gradient_index
  integer(kind=i4_kind), intent(in):: spin
  type(arrmat5)     , intent(inout) :: grad_cart(:)   ! cartesian grad. array
  real(kind=r8_kind), intent(in   ) :: grad_totsym(:,:,:) ! symm. adapt. contr.
  !
  integer(kind=i4_kind) :: i_unique,i_center,i_equal,index,i,grad_dim,vl,k
 
  real(kind=r8_kind),pointer :: rotmat(:,:)
  vl=size(grad_totsym,1)
 
  do i_unique=1,N_moving_unique_atoms
     i_center=moving_unique_atom_index(i_unique)
     index=cpks_gradient_index(i_unique)
     grad_dim=cpks_gradient_index(i_unique+1)-index
     index=index-1
     do i_equal=1,unique_atoms(i_center)%n_equal_atoms
        rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
        do i=1,grad_dim
          do k=1,3
           grad_cart(i_unique)%m(:vl,:,k,i_equal,spin) = &
             grad_cart(i_unique)%m(:vl,:,k,i_equal,spin) + &
                spread(rotmat(i,:),1,vl)*spread(grad_totsym(:vl,k,index+i),2,3)
          enddo
        end do
     end do
  enddo
 end subroutine transform_to_cart_gragamma

 subroutine transform_to_cart_dervs(dervs_cart,dervs_totsym,i_unique1,spin)
   ! purpose: transform symmetry adapted gradient components to cartesian
   !          coordinates and add them to an array of cartesian gradients
 
  use cpksdervs_matrices, only:cpks_gradient_index

  integer(kind=i4_kind), intent(in):: spin,i_unique1
  type(arrmat6), intent(inout) :: dervs_cart(:,:)   ! cartesian grad. array
  type(arrmat5), intent(in   ) :: dervs_totsym(:) ! symm. adapt. contr.
  !
  integer(kind=i4_kind) :: i_unique2,i_center,i_equal1,i_equal2,index,i,grad_dim,vl,k
 
  real(kind=r8_kind),pointer :: rotmat(:,:)
  vl=size(dervs_totsym(i_unique1)%m,1)

 
  ua2: do i_unique2=1,N_moving_unique_atoms
     i_center=moving_unique_atom_index(i_unique2)
     index=cpks_gradient_index(i_unique2)
     grad_dim=cpks_gradient_index(i_unique2+1)-index
     index=index-1
     eq2: do i_equal2=1,unique_atoms(i_center)%n_equal_atoms
        rotmat=>unique_atom_grad_info(i_unique2)%m(:,:,i_equal2)
        grad: do i=1,grad_dim
         eq1: do i_equal1=1,unique_atoms(i_unique1)%n_equal_atoms
          carts: do k=1,3
           dervs_cart(i_unique1,i_unique2)%m(:,k,i_equal1,:,i_equal2,spin) = &
             dervs_cart(i_unique1,i_unique2)%m(:,k,i_equal1,:,i_equal2,spin) + &
              spread(rotmat(i,:),1,vl)* &
                spread(dervs_totsym(i_unique1)%m(:,k,i_equal1,spin,index+i),2,3)
          enddo carts
          enddo eq1
        enddo grad
     enddo eq2
  enddo ua2
 end subroutine transform_to_cart_dervs


 subroutine transform_to_cart_dervs_grarho(dervs_cart,dervs_totsym,i_unique1,spin)
   ! purpose: transform symmetry adapted gradient components to cartesian
   !          coordinates and add them to an array of cartesian gradients
 
  use cpksdervs_matrices, only:cpks_gradient_index

  integer(kind=i4_kind), intent(in):: spin,i_unique1
  type(arrmat7), intent(inout) :: dervs_cart(:,:)   ! cartesian grad. array
  type(arrmat6), intent(in   ) :: dervs_totsym(:) ! symm. adapt. contr.
  !
  integer(kind=i4_kind) :: i_unique2,i_center,i_equal1,i_equal2,index,i,grad_dim,vl,k,kk
 
  real(kind=r8_kind),pointer :: rotmat(:,:)
  vl=size(dervs_totsym(i_unique1)%m,1)

 
  ua2: do i_unique2=1,N_moving_unique_atoms
     i_center=moving_unique_atom_index(i_unique2)
     index=cpks_gradient_index(i_unique2)
     grad_dim=cpks_gradient_index(i_unique2+1)-index
     index=index-1
     eq2: do i_equal2=1,unique_atoms(i_center)%n_equal_atoms
        rotmat=>unique_atom_grad_info(i_unique2)%m(:,:,i_equal2)
        grad: do i=1,grad_dim
         eq1: do i_equal1=1,unique_atoms(i_unique1)%n_equal_atoms
          carts: do k=1,3
                 do kk=1,3
           dervs_cart(i_unique1,i_unique2)%m(:,k,i_equal1,:,i_equal2,kk,spin) = &
           dervs_cart(i_unique1,i_unique2)%m(:,k,i_equal1,:,i_equal2,kk,spin) + &
              spread(rotmat(i,:),1,vl)* &
                spread(dervs_totsym(i_unique1)%m(:,k,i_equal1,kk,spin,index+i),2,3)
          enddo
          enddo carts
          enddo eq1
        enddo grad
     enddo eq2
  enddo ua2
 end subroutine transform_to_cart_dervs_grarho

 subroutine transform_to_cart_dervs_grarho2(dervs_cart,dervs_totsym,i_unique1,spin,i_eq)
   ! purpose: transform symmetry adapted gradient components to cartesian
   !          coordinates and add them to an array of cartesian gradients
 
  use cpksdervs_matrices, only:cpks_gradient_index

  integer(kind=i4_kind), intent(in):: spin,i_unique1
  integer(kind=i4_kind), optional, intent(in):: i_eq
  type(arrmat5), intent(inout) :: dervs_cart(:,:)   ! cartesian grad. array
  type(arrmat4), intent(in   ) :: dervs_totsym(:) ! symm. adapt. contr.
  !
  integer(kind=i4_kind) :: i_unique2,i_center,i_equal1,i_equal2,index,i,grad_dim,k
 
  real(kind=r8_kind),pointer :: rotmat(:,:)
 
  ua2: do i_unique2=1,N_moving_unique_atoms
     i_center=moving_unique_atom_index(i_unique2)
     index=cpks_gradient_index(i_unique2)
     grad_dim=cpks_gradient_index(i_unique2+1)-index
     index=index-1
     eq2: do i_equal2=1,unique_atoms(i_center)%n_equal_atoms
        rotmat=>unique_atom_grad_info(i_unique2)%m(:,:,i_equal2)
        grad: do i=1,grad_dim
         eq1: do i_equal1=1,unique_atoms(i_unique1)%n_equal_atoms
          carts: do k=1,3
           dervs_cart(i_unique1,i_unique2)%m(k,i_equal1,:,i_equal2,spin) = &
           dervs_cart(i_unique1,i_unique2)%m(k,i_equal1,:,i_equal2,spin) + &
              rotmat(i,:)* dervs_totsym(i_unique1)%m(k,i_equal1,spin,index+i)
          enddo carts
          if(present(i_eq)) exit 
          enddo eq1
        enddo grad
     enddo eq2
  enddo ua2
 end subroutine transform_to_cart_dervs_grarho2
#endif

end module density_calc_cpks

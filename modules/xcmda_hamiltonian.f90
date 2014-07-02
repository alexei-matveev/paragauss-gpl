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
  !===================================================================
! Public interface of module
  !===================================================================
module xcmda_hamiltonian
  !-------------------------------------------------------------------
  !
  !  Purpose: This module creates the XC-part of the hamiltionian
  !           within the model density approach (MDA)
  !
  !  Author: Uwe Birkenheuer
  !  Date: 6/97
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: Potential extended model density approach introduced
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  !------------ Modules used -----------------------------------------
# include "def.h"
  use type_module  ! type specification parameters
  use time_module, only: init_timer, start_timer, stop_timer
  use timer_module
  use orbitalstore_module
  use orbital_module
  use comm_module
  use msgtag_module
  use iounitadmin_module
  use density_calc_module
  use mat_charge_module
  use fit_coeff_module
  use linsys_module
  use unique_atom_module
  implicit none
  private
  save

  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------

  public :: xcmda_setup
  public :: xcmda_build
  public :: xcmda_close
  public :: mda_options_read
  public :: mda_options_write
  public :: mda_options_bcast
  public :: xcmda_get_exc
  public :: mda_constrain_rhospin
  public :: xcmda_coeff_store
  public :: xcmda_coeff_recover
  public :: mda_rho_shape_eps
  public :: mda_rho_cutoff
  public :: erf

  !===================================================================
  ! End of public interface of module
  !===================================================================

!..............................................................................
! model density approach (basic version)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! E_X[rho] = E_xc,num[rho_fit]
! with
! d/drho(r) rho_fit(s) = Sum(k,l) f_k(s) G_kl V_H[f_l](r)
! ==>
! V_X[rho](r) = Sum(k,l) < V_xc,num[rho_fit] | f_k > G_kl V_H[f_l](r)
!
! potential extended model density approach
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! E_X[rho] = E_xc,num[rho_fit] + < rho - rho_fit | V_xc,fit >
! with
! d/drho(r) rho_fit(s) = Sum(k,l) f_k(s) G_kl V_H[f_l](r)
! and
! d/dV_xc(r) V_xc,fit(s) = Sum(k,l) g_k(s) Y_kl g_l(r)
! where
! V_xc(r) = V_xc,num[rho_fit](r)
! ==>
! V_X[rho](r) = V_fit(r) + Sum(k,l) p_k G_kl V_H[f_l](r)
! with
! p_k = < V_xc,num[rho_fit] - V_xc,fit | f_k > +
!     + Sum(m,n) < rho - rho_fit | g_m > Y_mn { g_n | d/drho V_xc,num | f_k }
! where
! { g | d/drho V[rho] | f } := Int g(s) d/drho(r) V[rho](s) f(r) d^3s d^3r
!..............................................................................
  ! ----------- Default values for input parameters -------------
  real(kind=r8_kind) :: df_rho_shape_eps     = 1.0E-12_r8_kind
  real(kind=r8_kind) :: df_rho_cutoff        = 1.0E-50_r8_kind
  logical            :: df_constrain_rhospin = .true.
!:ANALYZE_MDA[
  logical            :: df_rho_exact = .false.
  logical            :: df_comp_exact = .false.
  logical            :: df_no_first = .false.
  real(kind=r8_kind) :: df_atomic_radius = 0.0_r8_kind
!:ANALYZE_MDA]

  ! ----------- MDA input parameters -----------------------------
  real(kind=r8_kind) :: rho_shape_eps, rho_cutoff
  logical            :: constrain_rhospin
!:ANALYZE_MDA[
  logical, public    :: rho_exact
  logical, public    :: comp_exact
  logical, public    :: no_first
  real(kind=r8_kind) :: atomic_radius, help_sum, help,help1
!:ANALYZE_MDA]

  !------------ private global arrays, switches and dimensions ---
  real(kind=r8_kind),allocatable :: rho(:,:)   , & ! rho(r,s)
                                    rho_ex(:,:), & ! rho_exact(r,s)
                                    fxc(:)     , & ! f_xc(rho(r),gamm(r))
                                    fxc_ex(:)  , & !
                                    dfdrho(:,:), & ! d/drho(r,s) f_xc(...)
                                    dfdrho_ex(:,:), & !
                                    dvdrho(:,:), & ! d/drho(r,s) V_xc(...,t)
                                    rhs_xc(:,:), & ! <g_l|V_xc(rho,s)>
!:TST[ testing
                                    rhs_lc(:,:), & !
!:TST]
                                    lhs_xc(:)  , & ! <g_k|g_l>
                                    rhs(:,:)   , & ! <f_k|V_xc(rho,s)>
                                    rhs_ex(:,:)    !
  real(kind=r8_kind),allocatable :: grad_rho(:,:,:), & ! ! Nabla rho(r,s)
                                    gamm(:,:)     , & ! gamm(r,t)
                                    dfdgrarho(:,:)     ! d/dgamm(r,t) f_xc(...)
  real(kind=r8_kind),allocatable :: grad_rho_ex(:,:,:), & ! storage for
                                    gamm_ex(:,:)     , & ! calc with
                                    dfdgrarho_ex(:,:)     ! exact dens
  ! private working arrays
  real(kind=r8_kind), allocatable :: tmp(:), hlp(:,:), theta(:,:), delta(:,:)
  real(kind=r8_kind), allocatable :: tmp_ex(:), hlp_ex(:,:)
  real(kind=r8_kind), allocatable :: tmp_st(:,:)

  logical :: &
       nl_calc, &
       becke88, &
       perdew, &
       xalpha, &
       vwn, &
       perdewwang91x, &
       perdewwang91c, &
       baerends94, &
       rxalpha, &
       rvwn, &
       ecmv92, &
       nrecmv92, &
       rbecke88, &
       rperdewwang91x, &
       rperdewwang91c, &
       ext_mda, &
       matinv_required, &
       matinv_exists, &
       eval_lincorr, &
       xalpha_only, &
       is_first_loop, &
       hcth_x, &
       hcth_c


  integer(kind=i4_kind) :: ispin,n_ch,n_xc,n_mat,vec_length
  real(kind=r8_kind)  :: charge_int(4), & ! Q(up) Q(dn) Qtrunc(up) Qtrunc(dn)
                         exc_int, exc_corr, exc_int_ex ! E_xc , E_xc(lin.corr.)
  type(orbital_type)          :: fcts_ch  ! charge fitting functions
!!! MF for testing
  type(orbital_type)          :: fcts_cl  ! charge fitting functions
  type(orbital_gradient_type) :: grads_ch ! charge fitting function gradients
  type(orbital_type)          :: fcts_xc  ! exchange fitting functions
!:ANALYZE_MDA[
  type(orbital_type)         ,pointer :: orbs_ob(:)  ! orbital basis functions
  type(orbital_gradient_type),pointer :: grads_ob(:) ! and their gradients
!:ANALYZE_MDA]

!:ANALYZE_MDA[
! namelist /mda_options/ rho_shape_eps, rho_cutoff, constrain_rhospin
  namelist /mda_options/ rho_shape_eps, rho_cutoff, constrain_rhospin, &
                         rho_exact, comp_exact, no_first, atomic_radius
!:ANALYZE_MDA]

contains

  subroutine xcmda_build (lh)
    !
    ! Wrapper for build_xcmda is called in every scf- cycle by
    ! main_scf().
    !
    implicit none
    logical, intent(in) :: lh
    !** End of interface *****************************************

    integer (i4_kind) :: info

    DPRINT 'in xcmda_build'
    print*,'xcmda_build'

    is_first_loop=lh

    if (ext_mda) call prepare_build_xcmda

    if (comm_parallel()) then
       if (comm_i_am_master()) then
          call comm_init_send (comm_all_other_hosts, msgtag_build_xcmda)
          call commpack(matinv_required,info)
          if(info/=0) call error_handler&
               ('Error packing matinv_required in xcmda_build')
          call commpack(matinv_exists,info)
          if(info/=0) call error_handler&
               ('Error packing matinv_exists in xcmda_build')
          call commpack(eval_lincorr,info)
          if(info/=0) call error_handler&
               ('Error packing eval_lincorr in xcmda_build')
          if (eval_lincorr) then
             call commpack(coeff_deltarho(1,1),n_xc*ispin,1,info)
             if(info/=0) call error_handler&
                  ('Error packing coeff_deltarho in xcmda_build')
          endif
          call commpack(is_first_loop,info)
          if(info/=0) call error_handler&
               ('Error packing matinv_required in xcmda_build')
          call comm_send()
       else
          call comm_save_recv (comm_master_host, msgtag_build_xcmda)
          call communpack(matinv_required,info)
          if(info/=0) call error_handler&
               ('Error unpacking matinv_required in build_xcmda')
          call communpack(matinv_exists,info)
          if(info/=0) call error_handler&
               ('Error unpacking matinv_exists in build_xcmda')
          call communpack(eval_lincorr,info)
          if(info/=0) call error_handler&
               ('Error unpacking eval_lincorr in build_xcmda')
          if (eval_lincorr) then
             allocate( coeff_deltarho(n_xc,ispin), stat=info )
             if (info /= 0) call error_handler&
                  ('Error allocating coeff_deltarho in build_xcmda')
             call communpack(coeff_deltarho(1,1),n_xc*ispin,1,info)
             if(info/=0) call error_handler&
                  ('Error unpacking coeff_deltarho in build_xcmda')
          endif
          call communpack(is_first_loop,info)
          if(info/=0) call error_handler&
               ('Error unpacking is_first_loop in build_xcmda')
       endif
    end if

    DPRINT 'call build_xcmda'
    call build_xcmda
  end subroutine xcmda_build

  !***********************************************************

  subroutine xcmda_allocate(soft_truncation)
    use machineparameters_module, only: machineparameters_veclen
    logical :: soft_truncation
    ! purpose: allocates non-permanent data
    integer :: status
    if (nl_calc) then ! GGA calculation
       call orbital_setup(machineparameters_veclen,do_gradients=.true.)
       call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch)
!      if(comp_exact) call fit_fct_allocate("ch",fcts=fcts_cl)
       if (ext_mda) call fit_fct_allocate("xc",fcts=fcts_xc)
       if (rho_exact) call orbital_allocate(orbs_ob=orbs_ob,grads=grads_ob)
!      if (comp_exact) call orbital_allocate(orbs_ob=orbs_ob,grads=grads_ob)
    else ! LDA calculation
       call orbital_setup(machineparameters_veclen)
       call fit_fct_allocate("ch",fcts=fcts_ch)
!      if(comp_exact) call fit_fct_allocate("ch",fcts=fcts_cl)
       if (ext_mda) call fit_fct_allocate("xc",fcts=fcts_xc)
!      if (comp_exact) call orbital_allocate(orbs_ob=orbs_ob)
       if (rho_exact) call orbital_allocate(orbs_ob=orbs_ob)
    endif
    allocate( tmp(vec_length), theta(vec_length,ispin), stat=status)
!   if(comp_exact) allocate(tmp_ex(vec_length), stat=status)
!   if(comp_exact) allocate(tmp_st(vec_length,2), stat=status)
    if (status /= 0) call error_handler&
         ('xcmda_allocate: allocation of tmp and theta failed')
    if (nl_calc) then
       allocate( hlp(vec_length,3), stat=status)
       if (status /= 0) call error_handler&
            ('xcmda_allocate: allocation of hlp failed')
!      if(comp_exact) allocate( hlp_ex(vec_length,3), stat=status)
    endif
    if (nl_calc .or. eval_lincorr) then
       if (soft_truncation) then
          allocate( delta(vec_length,ispin), stat=status)
          if (status /= 0) call error_handler&
               ('xcmda_allocate: allocation of delta failed')
       endif
    endif
    if (matinv_required) then
       allocate( lhs_xc(n_mat), stat=status)
       if (status /= 0) call error_handler&
            ('xcmda_allocate: allocation of lhs_xc failed')
    endif
  end subroutine xcmda_allocate

  !***********************************************************
  subroutine xcmda_deallocate(soft_truncation)
    logical :: soft_truncation
    ! purpose: deallocates non-permanent data
    integer :: status
    if (nl_calc) then ! GGA calculation
!:ANALYZE_MDA[
       if (rho_exact) call orbital_free(orbs_ob=orbs_ob,grads=grads_ob)
!      if (comp_exact) call orbital_free(orbs_ob=orbs_ob,grads=grads_ob)
!:ANALYZE_MDA]
       if (ext_mda) call fit_fct_free("xc",fcts=fcts_xc)
       call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch)
!      if(comp_exact) call fit_fct_free("ch",fcts=fcts_cl)
       call orbital_shutdown()
    else ! LDA calculation
!:ANALYZE_MDA[
!      if (comp_exact) call orbital_free(orbs_ob=orbs_ob)
       if (rho_exact) call orbital_free(orbs_ob=orbs_ob)
!:ANALYZE_MDA]
       if (ext_mda) call fit_fct_free("xc",fcts=fcts_xc)
       call fit_fct_free("ch",fcts=fcts_ch)
!      if(comp_exact) call fit_fct_free("ch",fcts=fcts_cl)
       call orbital_shutdown()
    endif
    deallocate( tmp, theta, stat=status)
    if (status /= 0) call error_handler&
         ('xcmda_deallocate: deallocation of tmp and theta failed')
!   if(comp_exact)  deallocate(tmp_ex,stat=status)
!   if(comp_exact)  deallocate(tmp_st,stat=status)
    if (nl_calc) then
       deallocate( hlp, stat=status)
       if (status /= 0) call error_handler&
            ('xcmda_deallocate: deallocation of hlp failed')
!      if(comp_exact) deallocate( hlp_ex, stat=status)
    endif
    if (nl_calc .or. eval_lincorr) then
       if (soft_truncation) then
          deallocate( delta, stat=status)
          if (status /= 0) call error_handler&
               ('xcmda_deallocate: deallocation of delta failed')
       endif
    endif
    if (matinv_required) then
       deallocate( lhs_xc, stat=status)
       if (status /= 0) call error_handler&
            ('xcmda_deallocate: deallocation of lhs_xc failed')
!:TST[
       write(1,*)'@@ deallocate lhs_xc in xcmda_deallocate'
!:TST]
    endif
  end subroutine xcmda_deallocate

  !***********************************************************

  subroutine xcmda_setup()
    ! purpose : routine performs the necessary allocations
    !           additionally some informations from
    !           symmetry_data_module are stored in module private
    !           variables
    !** End of interface **************************************
    use options_module          , only: options_n_spin, options_xcmode, &
                                        xcmode_extended_mda,            &
                                        options_orbitals_in_memory,     &
                                        options_orbitals_on_file
!!$    use xc_hamiltonian          , only: xc_give_control
    use xc_cntrl, only: xc_cntrl_give
    use machineparameters_module, only: machineparameters_veclen

    integer(kind=i4_kind) :: alloc_stat

    vec_length=machineparameters_veclen
    ! Initialize orbitalstore module
    if(options_orbitals_in_memory() .or. options_orbitals_on_file() ) then
       call orbitalstore_setup(vec_length)
!       print*,'done orbitalstore_setup'
    end if
    n_ch=fit_coeff_n_ch()
    n_xc=fit_coeff_n_xc()
    n_mat=(n_xc*(n_xc+1))/2
!!!MF test
    ext_mda=options_xcmode()==xcmode_extended_mda
    matinv_required=ext_mda ! dynamic control parameter reset by build_xcmda
    matinv_exists=.false.   ! dynamic control parameter reset by build_xcmda
    eval_lincorr=.false.    ! dynamic control parameter reset by build_xcmda
    ispin=options_n_spin()
!:UB[ added
    exc_corr = 0.0_r8_kind
!:UB]


!!$    call error_handler("xcmda/xcmda_setup: modification of xc_cntrl needed")
    call xc_cntrl_give(&
         lbecke88=becke88, &
         lperdew=perdew, &
         lxalpha=xalpha, &
         lrxalpha=rxalpha, &
         lvwn=vwn, &
         lperdewwang91x=perdewwang91x, &
         lperdewwang91c=perdewwang91c,&
         lrvwn=rvwn, &
         lrperdewwang91x=perdewwang91x, &
         lrperdewwang91c=perdewwang91c, &
         lbaerends94=baerends94, &
         lhcth_x=hcth_x, &
         lhcth_c=hcth_c, &
         lecmv92=ecmv92, &
         lnrecmv92=nrecmv92, &
         lrbecke88=rbecke88 &
         )

    nl_calc = perdew &
         .or. becke88 &
         .or. perdewwang91c &
         .or. perdewwang91x &
         .or. baerends94 &
         .or. ecmv92 &
         .or. nrecmv92 &
         .or. rbecke88 &
         .or. hcth_x &
         .or. hcth_c

!:UB[ removed
!:NO    ! allocation of coeff_xcmda moved to fit_coeff_allocate
!:NO    if( comm_i_am_master()) then
!:NO       allocate(coeff_xcmda(n_ch,ispin),stat=alloc_stat)
!:NO       if(alloc_stat/=0) call error_handler&
!:NO            ('allocation failed in su xcmda_setup')
!:NO    endif
!:UB]

    allocate(rho(vec_length,ispin),fxc(vec_length),dfdrho(vec_length,ispin), &
             rhs(n_ch,ispin),stat=alloc_stat)
!:MF[ for testing
!   if(comp_exact) then
!      allocate(rho_ex(vec_length,ispin),fxc_ex(vec_length),&
!            dfdrho_ex(vec_length,ispin), &
!             rhs_ex(n_ch,ispin),stat=alloc_stat)
!:MF]
!   endif
!:UB[ modified
    if (alloc_stat/=0) call error_handler(&
         'allocation (1) failed in xcmda_setup')
!:UB]
!:TST[ testing
    allocate(rhs_lc(n_ch,ispin),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler(&
         'allocation (1t) failed in xcmda_setup')
!:TST]
    if (nl_calc) then
       allocate(grad_rho(vec_length,3,ispin),gamm(vec_length,2*ispin-1), &
                dfdgrarho(vec_length,2*ispin-1),stat=alloc_stat)
!:MF[ for testing
!      if(comp_exact) then
!         allocate(grad_rho_ex(vec_length,3,ispin),&
!               gamm_ex(vec_length,2*ispin-1), &
!                dfdgrarho_ex(vec_length,2*ispin-1),stat=alloc_stat)
!      endif
!:MF]
       if (alloc_stat/=0) call error_handler(&
            'allocation (2) failed in xcmda_setup')
    endif
    if (ext_mda) then
       allocate(dvdrho(vec_length,2*ispin-1),rhs_xc(n_xc,ispin),stat=alloc_stat)
       if (alloc_stat/=0) call error_handler(&
            'allocation (3) failed in xcmda_setup')
    endif
    if (rho_exact) call density_calc_setup()
    if (.not.rho_exact) call fitted_density_calc_setup()
  end subroutine xcmda_setup

  !***************************************************************

  subroutine xcmda_close()
    ! purpose : perform the deallocation
!:UB[ removed
!:NO    !           keep coeff_xcmda for the analytical gradients
!:UB]
    !** End of interface *****************************************

    ! Close orbitalstore_module
    call orbitalstore_shutdown()
    deallocate(rho,fxc,dfdrho,rhs)

!:TST[ testing
!   if(comp_exact) deallocate(rho_ex,fxc_ex,dfdrho_ex,rhs_ex)

    deallocate(rhs_lc)
!:TST]
    if (nl_calc) deallocate(grad_rho,gamm,dfdgrarho)
!   if (nl_calc.and.comp_exact) deallocate(grad_rho_ex,gamm_ex,dfdgrarho_ex)
    if (ext_mda) deallocate(dvdrho,rhs_xc)
!:TST[
    if (ext_mda) &
    write(1,*)'@@ deallocate dvdrho in xcmda_close'
    if (ext_mda) &
    write(1,*)'@@ deallocate rhs_xc in xcmda_close'
!:TST]
!:ANALYZE_MDA[
!   if (rho_exact.or.comp_exact) call density_calc_close()
    if (rho_exact) call density_calc_close()
    if (.not.rho_exact) call fitted_density_calc_close()
!:ANALYZE_MDA]
  end subroutine xcmda_close

  !***************************************************************

  subroutine xcmda_clear()
    ! purpose : set some variables to zero
    !** End of interface *****************************************
    rho=0.0_r8_kind
    fxc=0.0_r8_kind
    dfdrho=0.0_r8_kind
!   if(comp_exact) then
!      rho_ex=0.0_r8_kind
!      fxc_ex=0.0_r8_kind
!      dfdrho_ex=0.0_r8_kind
!   endif
    if(nl_calc)then
      grad_rho=0.0_r8_kind
      gamm=0.0_r8_kind
      dfdgrarho=0.0_r8_kind

!     if(comp_exact) then
!       grad_rho_ex=0.0_r8_kind
!       gamm_ex=0.0_r8_kind
!       dfdgrarho_ex=0.0_r8_kind
!     endif

    endif
    if(eval_lincorr)then
      dvdrho=0.0_r8_kind
    endif
  end subroutine xcmda_clear

  !***************************************************************

   subroutine xcmda_send_tomaster()
     ! purpose : send results from numerical integration to the master
     !** End of interface ***************************************
     integer(kind=i4_kind) :: info

     call comm_init_send(comm_master_host,msgtag_xcmdaham_send)
     call commpack(charge_int,2*ispin,1,info)
     if(info/=0) call error_handler&
          ('Error packing the charge in su xcmda_send_tomaster')
     call commpack(exc_int,info)
     if(info/=0) call error_handler&
          ('Error packing the xc-energy in su xcmda_send_tomaster')
     call commpack(rhs(1,1),n_ch*ispin,1,info)
     if(info/=0) call error_handler&
          ('Error packing the rigth side of xc-fit')
!    if(comp_exact) then
!       call commpack(exc_int_ex,info)
!       if(info/=0) call error_handler&
!         ('Error packing the xc-energy_ex in su xcmda_send_tomaster')
!       call commpack(rhs_ex(1,1),n_ch*ispin,1,info)
!       if(info/=0) call error_handler&
!         ('Error packing the rigth side of xc-fit_ex')
!    endif
!:TST[ testing
     call commpack(rhs_lc(1,1),n_ch*ispin,1,info)
     if(info/=0) call error_handler&
          ('Error packing the rigth side of xc-fit (testing)')
!:TST]
      if (ext_mda) then
        call commpack(rhs_xc(1,1),n_xc*ispin,1,info)
        if(info/=0) call error_handler&
             ('Error packing the rigth side of extended xc-fit')
!:TST[
      write(1,*)'@@ send rhs_xc to the master in xcmda_send_tomaster'
      write(1,*)'@@      with n_xc,ispin = ',n_xc,ispin
!:TST]
        if (matinv_required) then
           call commpack(lhs_xc,n_mat,1,info)
           if(info/=0) call error_handler&
                ('Error packing the left side of extended xc-fit')
!:TST[
      write(1,*)'@@ send lhs_xc to the master in xcmda_send_tomaster'
      write(1,*)'@@      with n_mat = ',n_mat
!:TST]
        endif
     endif
     call comm_send()
   end subroutine xcmda_send_tomaster

   !***************************************************************

   subroutine xcmda_receive()
     ! purpose : receive results from numerical integration on the slaves
     !** End of interface *****************************************
     real(kind=r8_kind), allocatable :: help_arr(:,:), help_vec(:)
     real(kind=r8_kind)              :: help_ch(4), help_e
     integer(kind=i4_kind)           :: i,info,stat

     do i=1,comm_get_n_processors()-1
        call comm_save_recv(comm_all_other_hosts,msgtag_xcmdaham_send)
        if(comm_msgtag()/=msgtag_xcmdaham_send) call error_handler&
             ('Wrong msgtag in su xcmda_receive')
        call communpack(help_ch,2*ispin,1,info) ! unpacking integrated charge
        if(info/=0) call error_handler&
             ('Error unpacking the charge in su xcmda_receive')
        charge_int(:2*ispin) = charge_int(:2*ispin) + help_ch(:2*ispin)
        call communpack(help_e,info) ! unpacking integrated xc-energy
        if(info/=0) call error_handler&
             ('Error unpacking the xc energy in su xcmda_receive')
        exc_int = exc_int + help_e
        allocate( help_arr(n_ch,ispin), stat=stat )
        if (stat /= 0) call error_handler&
             ('Error allocating rhs_help for the xc(mda) fit')
        call communpack(help_arr(1,1),n_ch*ispin,1,info) ! unpacking rhs
        if(info/=0) call error_handler&
             ('Error unpacking the rhs for the xc(mda) fit')
        rhs = rhs + help_arr
!       if(comp_exact) then
!         call communpack(help_e,info) ! unpacking integrated xc-energy
!         if(info/=0) call error_handler&
!            ('Error unpacking the xc energy ex in su xcmda_receive')
!         exc_int_ex = exc_int_ex + help_e
!         help_arr=0.0_r8_kind
!         call communpack(help_arr(1,1),n_ch*ispin,1,info) ! unpacking rhs
!         if(info/=0) call error_handler&
!            ('Error unpacking the rhs for the xc(mda,ex) fit')
!         rhs_ex = rhs_ex + help_arr
!       endif
!:TST[ testing
        call communpack(help_arr(1,1),n_ch*ispin,1,info) ! unpacking rhs_lc
        if(info/=0) call error_handler&
             ('Error unpacking the rhs for the xc(mda) fit (testing)')
        rhs_lc = rhs_lc + help_arr
!:TST]
        deallocate( help_arr, stat=stat )
        if (stat /= 0) call error_handler&
             ('Error deallocating rhs_help for the xc(mda) fit')
        if (ext_mda) then
           allocate( help_arr(n_xc,ispin), stat=stat )
           if (stat /= 0) call error_handler&
                ('Error allocating rhs_xc_help for the xc(mda) fit')
           call communpack(help_arr(1,1),n_xc*ispin,1,info)
           if(info/=0) call error_handler&
                ('Error packing the rigth side of extended xc-fit')
           rhs_xc = rhs_xc + help_arr
           deallocate( help_arr, stat=stat )
           if (stat /= 0) call error_handler&
                ('Error deallocating rhs_xc_help for the xc(mda) fit')
           if (matinv_required) then
              allocate( help_vec(n_mat), stat=stat )
              if (stat /= 0) call error_handler&
                   ('Error allocating lhs_xc_help for the xc(mda) fit')
              call communpack(help_vec,n_mat,1,info)
              if(info/=0) call error_handler&
                   ('Error packing the left side of extended xc-fit')
              lhs_xc = lhs_xc + help_vec
              deallocate( help_vec, stat=stat )
              if (stat /= 0) call error_handler&
                   ('Error deallocating lhs_xc_help for the xc(mda) fit')
           endif
        endif
    end do
    1000 format("Num. integ. ",a14,":",F21.14:F22.14," (trunc.)")
    if (ispin == 1) then
       write(output_unit,1000)'charge density', &
             charge_int(1),charge_int(2)
    else
       write(output_unit,1000)'charge density', &
             charge_int(1)+charge_int(2),charge_int(3)+charge_int(4)
       write(output_unit,1000)'spin   density', &
             charge_int(1)-charge_int(2),charge_int(3)-charge_int(4)
    endif
    write(output_unit,1000)'xc-energy     ',exc_int
!    print *,'xc-energy     ',exc_int
   end subroutine xcmda_receive

   !***************************************************************

   subroutine prepare_build_xcmda
     ! purpose : performs preparatory processing of the extended
     !           charge density projection < rho - rho_fit | g_k >
     !
     ! only called on the master
     ! only called for the potential extended model density approach
     !
     !** End of interface *****************************************
     use output_module, only: output_chargefit
     integer(kind=i4_kind) :: i,s

     call init_timer(timer_grid_projections)
     call init_timer(timer_grid_fit_coeffs)

     eval_lincorr = matinv_exists .and. .not.pv_initialized

     if (.not.pv_initialized) then
        ! first load < rho_s - rho_fit,s | g_k > for s = up and down
        call start_timer(timer_grid_projections)
        coeff_deltarho = coeff_proj
        if (ispin == 1) then
           if (output_chargefit) then
              write(output_unit,*) 'COEFF_PROJ(tot)'
              write(output_unit,'(10F8.4)') coeff_deltarho(:,1)
           endif

           do i=1,n_ch
              coeff_deltarho(:,1) = coeff_deltarho(:,1) - &
                                    mat_xc_ch(:,i)*coeff_charge(i)
           end do
!+!:TST-TMP[
!+           if (n_xc == n_xc) then
!+              write(1,*)'@@ !! mat_xc_ch overwritten by mat_charge !!'
!+              coeff_deltarho = coeff_proj
!+              k = 0
!+              do i=1,n_ch
!+                 do j=1,i-1
!+                    k = k + 1
!+                    coeff_deltarho(j,1) = coeff_deltarho(j,1) - &
!+                                          mat_charge(k)*coeff_charge(i)
!+                    coeff_deltarho(i,1) = coeff_deltarho(i,1) - &
!+                                          mat_charge(k)*coeff_charge(j)
!+                 enddo
!+                 k = k + 1
!+                 coeff_deltarho(i,1) = coeff_deltarho(i,1) - &
!+                                       mat_charge(k)*coeff_charge(i)
!+              enddo
!+           endif
!+!:TST-TMP]
           if (output_chargefit) then
              write(output_unit,*) 'COEFF_DELTARHO(tot)'
              write(output_unit,'(10F8.4)') coeff_deltarho(:,1)
!+!:TST-TMP[
!+              coeff_xcmda(:,1) = 0.0_r8_kind
!+              where (charge_norm(:) /= 0.0_r8_kind)
!+                 coeff_xcmda(:,1) = coeff_deltarho(:,1) / charge_norm(:)
!+              end where
!+              write(output_unit,*) 'LAMBDA'
!+              write(output_unit,'(6F12.8)') coeff_xcmda(:,1)
!+!:TST-TMP]
           endif
        else
           ! load charge fit coefficients for spin up and down
           ! using coeff_xcmda as working array
           coeff_xcmda(:,1) = 0.5_r8_kind*( coeff_charge(:) + coeff_spin(:) )
           coeff_xcmda(:,2) = 0.5_r8_kind*( coeff_charge(:) - coeff_spin(:) )
           do s=1,ispin
              do i=1,n_ch
                 coeff_deltarho(:,s) = coeff_deltarho(:,s) - &
                                       mat_xc_ch(:,i)*coeff_xcmda(:,s)
              end do
           end do
           if (output_chargefit) then
              write(output_unit,*) 'COEFF_DELTARHO(up)'
              write(output_unit,'(10F8.4)') coeff_deltarho(:,1)
              write(output_unit,*) 'COEFF_DELTARHO(down)'
              write(output_unit,'(10F8.4)') coeff_deltarho(:,2)
           endif
        end if
        call stop_timer(timer_grid_projections)

        ! then evaluate the linear correction term to the xc energy
        ! exc_corr :=  Sum(s) < rho_s - rho_fit,s | V_xc,s,fit >
        do s=1,ispin
           exc_corr = exc_corr + sum( coeff_deltarho(:,s) * coeff_xc(:,s) )
        end do

        ! now determine the extended fit coefficients
        ! coeff_deltarho(l,s) := Sum(k) < rho_s - rho_fit,s | g_k > Y_kl
        call start_timer(timer_grid_fit_coeffs)
        if (matinv_exists) then
           do s=1,ispin
              call solve_linsys(n_xc,matinv_exchange,coeff_deltarho(:,s))
           end do
        endif
        call stop_timer(timer_grid_fit_coeffs)
     endif

   end subroutine prepare_build_xcmda

   !***************************************************************

   subroutine build_xcmda
     ! purpose : main routine for fiiting the xc-hamiltonian
     !           XC-Hamiltonian
     !
     use becke_perdew_module
     use perdew_wang_module
     use vwnc
     use output_module, only: output_chargefit
     use grid_module, only: more_grid, grid_loop_setup
     !** End of interface *****************************************
     logical :: need_zxc_off_diag
     real(kind=r8_kind),parameter :: zero=0.0_r8_kind, half=0.5_r8_kind, &
                                     two =2.0_r8_kind
     real(kind=r8_kind),pointer :: ob(:,:), gr(:,:,:)
     integer (i4_kind) :: i, j, k, l, s, t, vec_length_act
     real(kind=r8_kind) :: a2, zero_target, rho_min
     integer :: status
     real(kind=r8_kind), allocatable :: r1(:), r2(:)
     integer(kind=i4_kind) :: i_bas, i_exp, i_sa, n_sa, ih
     type(unique_atom_basis_type), pointer :: uab
     real(kind=r8_kind) :: alpha, alpha2, coeff, dr, pi, fact, fact2
     external error_handler
     real(kind=r8_kind),pointer :: grdpts(:,:), & ! grid points
                                   grdwts(:)      ! grid weights

     call init_timer(timer_grid_orbitals)
     call init_timer(timer_grid_density)
     call init_timer(timer_int_fit_density_calc) !AG
     call init_timer(timer_grid_trunc_dens)
     call init_timer(timer_grid_functionals)
     call init_timer(timer_grid_trunc_funcs)
     call init_timer(timer_grid_num_metric)
     if (.not.(comm_i_am_master() .and. ext_mda)) then
        call init_timer(timer_grid_projections)
        call init_timer(timer_grid_fit_coeffs)
     endif
     DPRINT 'timers initialized'

     a2 = rho_shape_eps*rho_shape_eps
     rho_min = 0.25_r8_kind * sqrt( epsilon(0.0_r8_kind) * a2 )
     need_zxc_off_diag = ispin > 1 .and. .not.(vwn.or.nl_calc)
     call xcmda_allocate(a2/=zero)

     charge_int=zero
     exc_int=zero
     exc_corr=zero
     exc_int_ex=zero
     rhs=zero
!    if(comp_exact) rhs_ex=zero
     rhs_lc=zero
     if (ext_mda) then
        rhs_xc=zero
        if (matinv_required) then
           lhs_xc=zero
        endif
     endif

     DPRINT 'loop over gridpoints'
     call grid_loop_setup()
     ! loop over gridpoints:
     do while(more_grid(vec_length,grdpts,grdwts) ) ! fetching part of the grid
        ! more_grid() will return
        !   grdpts => new batch of coordinates
        !   grdwts => corresponding weights
        ! not more than "vec_length" long

        ! this may be, in general, below vec_length, particularly in last
        ! iteration:
        vec_length_act = size(grdpts,1)

        ! zeroes global rho(:), fxc(:), etc, FIXME: make these vars local!
        call xcmda_clear()

        ! calculate the fitting functions on the grid
        call start_timer(timer_grid_orbitals)
        if (nl_calc) then ! GGA calculation
           call fit_fct_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
                "ch",fcts=fcts_ch,grads=grads_ch)
           DPRINT 'fit_fct_calculate done',shape(fcts_ch%o)
           if (ext_mda) & ! only LDA contributions so far (UB, 8/98)
           call fit_fct_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
                "xc",fcts=fcts_xc)
!:ANALYZE_MDA[
!          if (rho_exact.or.comp_exact) &
!          call orbital_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
!               orbs_ob=orbs_ob,grads=grads_ob)
!:ANALYZE_MDA]
        else ! LDA calculation
           call fit_fct_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
                "ch",fcts=fcts_ch)
!!! MF test
!          if(comp_exact) call fit_fct_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
!               "cl",fcts=fcts_cl)

           if (ext_mda) &
           call fit_fct_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
                "xc",fcts=fcts_xc)
!:ANALYZE_MDA[
!          if (rho_exact.or.comp_exact) &
           if (rho_exact) &
           call orbital_calculate(grdpts(1:vec_length_act,1:3),vec_length_act,&
                orbs_ob=orbs_ob)
!:ANALYZE_MDA]
        endif
        call stop_timer(timer_grid_orbitals)

        ! calculate the fitted density on the grid
        call start_timer(timer_grid_density)
        if (nl_calc) then ! GGA calculation
!:ANALYZE_MDA[
!AG[ text changes after merging (compare prev.vers !)
           if (rho_exact) then
             if (is_first_loop) then
                rho=zero
                gamm=zero
                grad_rho=zero
             else
                call density_calc_nl(vec_length_act,rho,gamm,&
                grad_rho,&
                orbs_ob,grads_ob)
             endif
           endif
!AG]
!          if (comp_exact) then
!            if(.not.rho_exact) then
!              call fitted_density_calc_close()
!              call density_calc_setup()
!            endif
!            if (is_first_loop) then
!               rho_ex=zero
!               gamm_ex=zero
!               grad_rho_ex=zero
!            else
!              call density_calc_nl(vec_length_act,rho_ex,gamm_ex,&
!               grad_rho_ex(:,1,:),grad_rho_ex(:,2,:),grad_rho_ex(:,3,:),&
!               orbs_ob,grads_ob)
!            endif
!            if(.not.rho_exact) then
!              call density_calc_close
!              call fitted_density_calc_setup()
!            endif
!          endif
!AG?(d)    endif
           DPRINT 'call fitted_density_calc'
           if (.not.rho_exact)&
                & call fitted_density_calc(vec_length_act,rho=rho,grad=grad_rho,&
                fcts=fcts_ch,grads=grads_ch)
           DPRINT 'call fitted_density_calc done'
!:ANALYZE_MDA]
        else ! LDA calculation
!:ANALYZE_MDA[
           if (rho_exact) then
             if (is_first_loop) then
                rho=zero
             else
               call density_calc(vec_length_act,rho,orbs_ob)
             endif
           endif
!          if (comp_exact) then
!             if(.not. rho_exact) then
!               call fitted_density_calc_close()
!               call density_calc_setup()
!             endif
!             if (is_first_loop) then
!               rho_ex=zero
!             else
!               call density_calc(vec_length_act,rho_ex,orbs_ob)
!             endif
!             if(.not. rho_exact) then
!               call density_calc_close
!               call fitted_density_calc_setup()
!             endif
!          endif
           if (.not.rho_exact) then
!:ANALYZE_MDA]
              call fitted_density_calc(vec_length_act,rho=rho,fcts=fcts_ch)
           endif
        endif
        call stop_timer(timer_grid_density)
!:ANALYZE_MDA]
        ! truncate the fitted density and
        ! integrate the density before and after truncation
        call start_timer(timer_grid_trunc_dens)
!..............................................................................
! density truncation
! ~~~~~~~~~~~~~~~~~~
! rho(trunc)_s(r) := S( rho_s(r) )  > 0
!
! soft truncation                             hard trunctation
! ~~~~~~~~~~~~~~~                             ~~~~~~~~~~~~~~~~
! S (x) = [ x + sqrt( a^2 + x^2 ) ] / 2       S (x) = max(x,b)
! S'(x) =    S(x) / sqrt( a^2 + x^2 )         S'(x) = Theta(x-b)
! S"(x) = 1/2 a^2 / sqrt( a^2 + x^2 )^3       S"(x) = delta(x-b) = 0 for x <> b
!
! Input : rho_s(r)
! Output: S(rho_(r)), S`(rho_s(r)), S``(rho_s(r))
!..............................................................................
        l = vec_length_act
        do s=1,ispin
           t = s+ispin
           charge_int(s) = charge_int(s) + sum( rho(:l,s) * grdwts(:l) )
!:ANALYZE_MDA]
           if (rho_exact) then
              theta = 1.0_r8_kind
              if ( (nl_calc .or. eval_lincorr) .and. a2 /= zero ) then
                 delta = 0.0_r8_kind
              endif
           else
!:ANALYZE_MDA]
           if (a2 == zero) then ! hard truncation
              rho(1:l,s) = rho(1:l,s) - rho_cutoff
              theta(:l,s) = half + sign(half,rho(:l,s))
              rho(:l,s) = rho_cutoff + theta(:l,s) * rho(:l,s)
           else ! soft truncation
!..............................................................................
! if rho < - a / sqrt(eps) with eps := min{ x | 1+x > 1 }_r8_kind
! the evaluation of S(rho) yields zero, numerically, though the correct
! value of S(rho) should be close to rho_min := sqrt(eps)*a/4
!..............................................................................
              tmp(:l  ) = sqrt( a2 + rho(:l,s)*rho(:l,s) )
              rho(:l,s) = half*( rho(:l,s) + tmp(:l) )
              if (nl_calc .or. eval_lincorr) then
                 delta(:l,s) = half * a2 / ( tmp(:l)*tmp(:l)*tmp(:l) )
              endif
              theta(:l,s) = rho(:l,s) / tmp(:l)
              rho(:l,s) = max( rho(:l,s) , rho_min ) ! to avoid numerical zero
           endif
!:ANALYZE_MDA[
           endif
!:ANALYZE_MDA]
           charge_int(t) = charge_int(t) + sum( rho(:l,s) * grdwts(:l) )
        enddo
        ! now build gamm of the truncated fitted density
!..............................................................................
! gamm_st(r) = < Nabla rho(trunc)_s(r) | Nabla rho(trunc)_t(r) >
! ==>
! gamm_st(r) = S'(rho_s(r)) < Nabla rho_s(r) | Nabla rho_t(r) > S'(rho_t(r))
!..............................................................................
        if (nl_calc) then
           !-> for exact rho already done
           gamm(:l,1) = grad_rho(:l,1,1)*grad_rho(:l,1,1) &
                       + grad_rho(:l,2,1)*grad_rho(:l,2,1) &
                       + grad_rho(:l,3,1)*grad_rho(:l,3,1)
           gamm(:l,1) = theta(:l,1) * gamm(:l,1) * theta(:l,1)
           if (ispin > 1) then
              gamm(:l,2) = grad_rho(:l,1,2)*grad_rho(:l,1,2) &
                          + grad_rho(:l,2,2)*grad_rho(:l,2,2) &
                          + grad_rho(:l,3,2)*grad_rho(:l,3,2)
              gamm(:l,2) = theta(:l,2) * gamm(:l,2) * theta(:l,2)
              gamm(:l,3) = grad_rho(:l,1,1)*grad_rho(:l,1,2) &
                          + grad_rho(:l,2,1)*grad_rho(:l,2,2) &
                          + grad_rho(:l,3,1)*grad_rho(:l,3,2)
              gamm(:l,3) = theta(:l,1) * gamm(:l,3) * theta(:l,2)
           endif
        endif
        call stop_timer(timer_grid_trunc_dens)

        ! evaluate the exchange correlation functionals
        call start_timer(timer_grid_functionals)
!       print*,'start grid_functionals'
!..............................................................................
! E_xc[rho_t] = Int(R^3) e_xc[rho_t](r) dV
! with
! e_xc[rho_t](r) = e_xc(rho_t(r),Nabla rho_t`(r))
! ==>
! < d/d rho_s E_xc[rho_t] | psi > =
! Int(R^3) V_xc,s[rho_t](r) psi(r) +  W_xc,s[rho_t](r) d/dr psi(r) dV
! with
! V_xc,s[rho_t](r) = [d/d       rho_s e_xc](rho_t(r),Nabla rho_t`(r))  and
! W_xc,s[rho_t](r) = [d/d Nabla rho_s e_xc](rho_t(r),Nabla rho_t`(r))
! ==> (for the LDA case)
! ( phi | d/d rho_s d/d rho_s` E_xc[rho_t] | psi ) =
! Int(R^3) phi(r) Z_xc,ss`[rho_t](r) psi(r) dV
! with
! Z_xc,ss`[rho_t](r) = [d/drho_s d/drho_s` e_xc](rho_t(r))
!
! e_xc(rho_t,Nabla rho_t') = f_xc(rho_t,gamm_tt')
! e_xc(rho_t,Nabla rho_t`) = f_xc(rho_t,gamma_tt`)
! with
! gamma_tt`(r) = < Nabla rho_t(r) | Nabla rho_t`(r) > ; tt` = 11,22,12
! ==>
! V_xc,s[rho_t](r) = [d/d rho_s f_xc](rho_t(r),gamma_tt`(r))  and
! gamm_tt`(r) = < Nabla rho_t(r) | Nabla rho_t'(r) > ; tt' = 11,22,12
! V_xc,s[rho_t](r) = [d/d rho_s f_xc](rho_t(r),gamm_tt`(r))  and
! W_xc,s[rho_t](r) = 2 G_xc,s s[rho_t](r) d/dr rho_ s(r)
!                  +   G_xc,s-s[rho_t](r) d/dr rho_-s(r)
! with
! G_xc,ss`[rho_t](r) := [d/d gamma_ss` f_xc](rho_t(r),gamma_tt`(r))
! and (for the LDA case)
! Z_xc,ss`[rho_t](r) = [d/drho_s d/drho_s` f_xc](rho(t))
!
! RKS : no spin dependencies _t, _tt`, _s, and _ss`
!       no mixed gamm derivative contribution G_xc,s-s[rho_t](r)
! LDA : no Nabla rho_t(r) and gamm_tt`(r) dependencies
!       no W_xc,s[rho_t] or W_xc[rho] potential contributions
!
! Input : S(rho_t(r)), gamm_tt`(r)
! Input : S(rho_t(r)), gamma_tt`(r)
! Output: e_xc[S(rho_t)](r), V_xc,s[S(rho_t)](r), G_xc,ss`[S(rho_t)](r), and
!         Z_xc,ss`[S(rho_t)](r)
!..............................................................................
! Up to now, only the Xalpha part of the xc energy is treated in
! potential extended fashion, during an extended MDA run. UB 8/98
        if ( vwn ) then
             call vwn_calc   (rho,dfdrho,ispin,fxc,l)
!            if(comp_exact) call vwn_calc(rho_ex,dfdrho_ex,ispin,fxc_ex,l)
             if ( perdewwang91c ) then
                ! keep fxc and dfdrho for use in perdewwang91c
                hlp(:l,1        ) = fxc   (:l  )
                hlp(:l,2:1+ispin) = dfdrho(:l,:)
                !!!test comp_exact ony for vwn and becke perdew
             end if
        endif
        if ( xalpha ) then
           if (eval_lincorr) then
              call xalpha_calc(rho,dfdrho,ispin,fxc,l,dvdrho=dvdrho, &
                               offset=.false.)
           else
              call xalpha_calc(rho,dfdrho,ispin,fxc,l,offset=.false.)
!             if(comp_exact) call xalpha_calc(rho_ex,dfdrho_ex,ispin,fxc_ex,l,offset=.false.)
           endif
        endif
        if ( perdew ) then
             call perdew_calcMDA(rho,gamm,dfdrho,ispin,fxc,dfdgrarho,l)
!            if(comp_exact) call perdew_calc(rho_ex,gamm_ex,dfdrho_ex,&
!                       ispin,fxc_ex,dfdgrarho_ex,l)
        endif
        if ( becke88 ) then
             call becke88_calcMDA(rho,gamm,dfdrho,ispin,fxc,dfdgrarho,l)
!             call becke88_calc(rho,gamm,dfdrho,ispin,fxc,dfdgrarho,l)
!            if(comp_exact) call becke88_calc(rho_ex,gamm_ex,dfdrho_ex,&
!                       ispin,fxc_ex,dfdgrarho_ex,l)
        endif
        if ( perdewwang91x ) then
             call pw91x_calc(rho,gamm,dfdrho,ispin,fxc,dfdgrarho,l)
        endif
        if ( perdewwang91c ) then
             call pw91c_calc(rho,gamm,dfdrho,ispin,fxc,dfdgrarho,l,&
                             fxc_lda=hlp(:,1),dfdrho_lda=hlp(:,2:1+ispin))
        endif
        exc_int = exc_int + sum( fxc(:l) * grdwts(:l) )
!       if(comp_exact) exc_int_ex = exc_int_ex + sum( fxc_ex(:l) * grdwts(:l) )
!print *,'ints by run',exc_int,exc_int_ex,comm_i_am_master()

        call stop_timer(timer_grid_functionals)

        ! start the exchange-correlation potential evaluation
        do s=1,ispin

           ! first load the weighted exchange-correlation potentials
           ! w_a V_xc,s[S(rho_t)](r_a) and w_a W_xc,s[S(rho_t)](r_a)
           call start_timer(timer_grid_functionals)
!..............................................................................
! W_xc,s[S(rho_t)](r) = 2 G_xc,s s[S(rho_t)](r) S`(rho_ s) d/dr rho_ s(r)
!                     +   G_xc,s-s[S(rho_t)](r) S`(rho_-s) d/dr rho_-s(r)
!..............................................................................
           if (nl_calc) then ! GGA co-potential
              ! 2 w_a G_xc,ss[rho_t](r_a) d/dr rho_s(r_a)
              tmp(1:l) = two * dfdgrarho(1:l,s) * grdwts(1:l) * theta(1:l,s)
              do k=1,3
                 hlp(1:l,k) = tmp(1:l) * grad_rho(1:l,k,s)
              end do
              if (ispin > 1) then
                 ! w_a G_xc,s-s[rho_t](r_a) d/dr rho_-s(r_a)
                 tmp(1:l) = dfdgrarho(1:l,3) * grdwts(1:l) * theta(1:l,3-s)
                 do k=1,3
                    hlp(1:l,k) = hlp(1:l,k) + &
                                 tmp(1:l) * grad_rho(1:l,k,3-s)
                 end do
              end if
!:MF[for testing
!             if(comp_exact) then
!               tmp_ex(1:l) = two * dfdgrarho_ex(1:l,s) * grdwts(1:l)
!               do k=1,3
!                hlp_ex(1:l,k) = tmp_ex(1:l) * grad_rho_ex(1:l,k,s)
!               end do
!               if (ispin > 1) then
!                tmp_ex(1:l) = dfdgrarho_ex(1:l,3) * grdwts(1:l)
!                do k=1,3
!                   hlp_ex(1:l,k) = hlp_ex(1:l,k) + &
!                                tmp_ex(1:l) * grad_rho_ex(1:l,k,3-s)
!                end do
!               end if
!             endif
!:MF]
           end if
           ! w_a V_xc,s[rho_t](r_a)

           tmp(1:l) = dfdrho(1:l,s) * grdwts(1:l)

!if(comp_exact.and.s==1) then
if(.false.) then
        !xc-density
        tmp_st(:,:)=0.0_r8_kind
        call fitted_density_calc(vec_length_act,rho=tmp_st,fcts=fcts_ch,&
                altern_coeff=coeff_xcmda_old)
        tmp_ex(:)=tmp_st(:,1)
        !mda xc-potential
        tmp_st(:,:)=0.0_r8_kind
        call fitted_density_calc(vec_length_act,rho=tmp_st,fcts=fcts_cl,&
                altern_coeff=coeff_xcmda_old)

        do ih=1,vec_length
!          if(abs(rhelp-sqrt(sum(grdpts(ih,1:3)**2)))>1.e-9) then
                write(output_unit,'(11es21.12E3)') &
                        grdpts(ih,1:3), &
                         grdwts(ih),&
                        !sqrt(sum(grdpts(ih,1:3)**2)),&
rho_ex(ih,s),rho(ih,s),tmp_ex(ih),dfdrho_ex(ih,s),dfdrho(ih,s),tmp_st(ih,1)
!                rhelp=sqrt(sum(grdpts(ih,1:3)**2))
!          endif
        enddo
endif
!          if(comp_exact) tmp_ex(1:l) = dfdrho_ex(1:l,s) * grdwts(1:l)


           call stop_timer(timer_grid_functionals)

           ! then perform the modifications due to the truncation
           call start_timer(timer_grid_trunc_funcs)
!..............................................................................
! E(trunc)_xc[rho_t] := E_xc[rho(trunc)_t] = E_xc[S(rho_t)]
! ==>
! e(trunc)_xc[rho_t](r) = e_xc(S(rho_t(r)),S`(rho_t`(r))*Nabla rho_t`(r))
! ==>
! V(trunc)_xc,s[rho_t](r) = V_xc,s[S(rho_t)](r) * S`(rho_s(r)) + ...
!                           W_xc,s[S(rho_t)](r) * S"(rho_s(r)) * Nabla rho_s(r)
! and
! W(trunc)_xc,s[rho_t](r) = W_xc,s[S(rho_t)](r) * S`(rho_s(r))
! and (for the LDA case)
! Z(trunc)_xc,ss`[rho_t](r) = Z_xc,ss`[S(rho_t)](r) * S`(rho_s(r)) S`(rho_s`(r))
!                           + V_xc,s  [S(rho_t)](r) * S"(rho_s(r)) delta_ss`
!
! RKS : no spin dependencies _t and _s
! LDA : no W_xc,s[rho_t] or W_xc[rho] contributions
!..............................................................................
           if (eval_lincorr) then
              ! at the moment only the Xalpha part is considered (UB,8/98) ==>
              need_zxc_off_diag = .false.
              ! w_a Z(trunc)_xc,ss`[rho_t](r_a)
              dvdrho(:l,s) = grdwts(:l) * dvdrho(:l,s) * theta(:l,s)*theta(:l,s)
              if (a2 /= zero) then ! soft truncation
                 dvdrho(:l,s) = dvdrho(:l,s) + tmp(:l) * delta(:l,s)
              endif
              if (s == 2 .and. need_zxc_off_diag) then
                 dvdrho(:l,3) = grdwts(:l)*dvdrho(:l,3)*theta(:l,1)*theta(:l,2)
              endif
           endif
           ! w_a V(trunc)_xc,s[rho_t](r_a)
           tmp(1:l) = tmp(1:l) * theta(:l,s)


           if (nl_calc) then ! GGA co-potential
              if (a2 /= zero) then ! soft truncation
                 tmp(1:l) = tmp(1:l) + &
                      hlp(1:l,1) * delta(:l,s) * grad_rho(1:l,1,s) + &
                      hlp(1:l,2) * delta(:l,s) * grad_rho(1:l,2,s) + &
                      hlp(1:l,3) * delta(:l,s) * grad_rho(1:l,3,s)
              end if
              ! w_a W(trunc)_xc,s[rho_t](r_a)
              do k=1,3
                 hlp(1:l,k) = hlp(1:l,k) * theta(:l,s)
              end do
!             if(comp_exact) -> no additional term, delta=zero
           end if
           call stop_timer(timer_grid_trunc_funcs)

           ! finally perform the potential projections < f_k | V_xc,s >num
           call start_timer(timer_grid_projections)
!          print*,'start grid_projections'
!..............................................................................
! < d/d rho_s E_xc[rho_t] | f_k > =
! Int(R^3) V_xc,s[rho_t](r) f_k(r) +  W_xc,s[rho_t](r) d/dr f_k(r) dV
!..............................................................................
           ob => fcts_ch%o(:,:,1)
           if (nl_calc) then ! GGA potential
              gr => grads_ch%o(:,:,:,1)
              do i=1,n_ch
                 rhs(i,s) = rhs(i,s) + sum( tmp(:l  ) * ob(:l,  i) + &
                                            hlp(:l,1) * gr(:l,1,i) + &
                                            hlp(:l,2) * gr(:l,2,i) + &
                                            hlp(:l,3) * gr(:l,3,i) )
!                if(comp_exact) then
!                   rhs_ex(i,s) = rhs_ex(i,s)+sum( tmp_ex(:l  ) * ob(:l,  i) + &
!                                           hlp_ex(:l,1) * gr(:l,1,i) + &
!                                           hlp_ex(:l,2) * gr(:l,2,i) + &
!                                           hlp_ex(:l,3) * gr(:l,3,i) )
!                endif
              end do
           else ! LDA potential
              do i=1,n_ch
                 rhs(i,s) = rhs(i,s) + sum( tmp(:l) * ob(:l,i) )
!                if(comp_exact) rhs_ex(i,s) = rhs_ex(i,s) + sum( tmp_ex(:l) * ob(:l,i) )
              end do
           endif
           ! and the same for < g_k | V_xc,s >num
           ! note, that up to now, only the LDA case is implemented (UB, 8/98)
           if (ext_mda) then
              ob => fcts_xc%o(:,:,1)
              do i=1,n_xc
                 rhs_xc(i,s) = rhs_xc(i,s) + sum( tmp(:l) * ob(:l,i) )
              end do
           endif
           call stop_timer(timer_grid_projections)
        enddo ! s

        ! add the extended projection contribution to rhs
        call start_timer(timer_grid_projections)
        ! ---------------------------------------------------------------------
        ! proj(r_a,s) = Sum(k) coeff_deltarho(k,s) * fcts_xc(k,r_a)
        ! and
        ! rhs(k,t) += Sum(a) w_a proj(r_a,s) dvdrho(r_a,s,t) fcts_ch(k,r_a)
        ! ---------------------------------------------------------------------
        if (eval_lincorr) then
           do s=1,ispin
              ! first load w_a * proj(r_a,s)
              ob => fcts_xc%o(:,:,1)
              theta(:l,s) = zero
              do i=1,n_xc
                 theta(:l,s) = theta(:l,s) + coeff_deltarho(i,s) * ob(:l,i)
              end do
              theta(:l,s) = theta(:l,s) * grdwts(:l)
              ! then add dvdrho(r_a,s,t) and perform the integration
              ob => fcts_ch%o(:,:,1)
              t = s
              tmp(:l) = theta(:l,s) * dvdrho(:l,s)
              do i=1,n_ch
                 rhs(i,t) = rhs(i,t) + sum( tmp(:l) * ob(:l,i) )
!:TST[ testing
                 rhs_lc(i,t) = rhs_lc(i,t) + sum( tmp(:l) * ob(:l,i) )
!:TST]
              end do
              if (need_zxc_off_diag) then
                 t = 3 - s
                 tmp(:l) = theta(:l,s) * dvdrho(:l,3)
                 do i=1,n_ch
                    rhs(i,t) = rhs(i,t) + sum( tmp(:l) * ob(:l,i) )
!:TST[ testing
                    rhs_lc(i,t) = rhs_lc(i,t) + sum( tmp(:l) * ob(:l,i) )
!:TST]
                 end do
              endif
           end do
        endif
        call stop_timer(timer_grid_projections)

        ! set up the numerical overlap matrix
        call start_timer(timer_grid_num_metric)
        if (matinv_required) then
           ob => fcts_xc%o(:,:,1)
           k = 0
           do i=1,n_xc
              tmp(:l) = ob(:l,i) * grdwts(:l)
              do j=1,i
                 k = k+1
                 lhs_xc(k) = lhs_xc(k) + sum( tmp(:l) * ob(:l,j) )
              end do
           end do
        endif
        call stop_timer(timer_grid_num_metric)
     enddo

     if(comm_i_am_master()) then
!:ANALYZE_MDA[
        if (rho_exact) &
        write(output_unit,*)'>>>TST: rho_exact is used for E_mda and V_mda !'
        pi = 4.0_r8_kind * atan( 1.0_r8_kind )
!:ANALYZE_MDA]
        ! receiving numerical integrals from the slaves
        call xcmda_receive()

        if (ext_mda) then
           ! solve the linear equation system to find V_xc,s,fit
           ! ------------------------------------------------------------------
           ! V_xc,s,fit(r) = Sum(k,l) = Sum(k) c_k,s g_k(r)
           ! with
           ! c_k,s = Sum(l) Y_kl < g_l | V_xc,s[rho_fit] >
           ! and
           ! Y_kl = {X^(-1)}_kl  where  X_kl = < g_k | g_l >_num
           ! ------------------------------------------------------------------
           if (matinv_required) then
              call start_timer(timer_grid_num_metric)
!:UB[ modified
!:TST[ <g_k|g_l> (please do not remove)
!             write(1,*)'compare <g_k|g_l>num and <g_k|g_l>anal'
!             write(1,1133)'i','j','lin','num.','anal.'
!             1122 format(2i3,i5,2f19.12)
!             1133 format(2a3,a5,2a19   )
!             k = 0
!             do i=1,n_xc
!                do j=1,i
!                   k = k+1
!                   zero_target = max( 1.0_r8_kind , abs(lhs_xc(k)) )*1.0E-2
!                   if ( abs(lhs_xc(k) - mat_exchange(k)) > zero_target ) then
!                      write(1,1122)i,j,k,lhs_xc(k),mat_exchange(k)
!                   endif
!                end do
!             end do
!:TST]
!:UB]
              call decompose_linsys(n_xc,lhs_xc)
              matinv_exchange = lhs_xc
!:TST[
              write(1,*)'@@ evaluate matinv_exchange[lhs_xc] in build_xcmda'
!:TST]
              call stop_timer(timer_grid_num_metric)

              ! don`t reset matinv_required here, it still controls deallocation
           endif
           call start_timer(timer_grid_fit_coeffs)
           do s=1,ispin
              call solve_linsys(n_xc,matinv_exchange,rhs_xc(:,s))
              coeff_xc(:,s) = rhs_xc(:,s)
           end do
!:TST[
           write(1,*)'@@ evaluate coeff_xc[rhs_xc] in build_xcmda'
!:TST]
           call stop_timer(timer_grid_fit_coeffs)
        endif

        ! solve the linear equation system to find V_X,s
        ! ---------------------------------------------------------------------
        ! -- Spin restricted MDA calculation
        !    V_mda,xc = Sum(k,l) <V_xc[rho]|f_k> G_kl V_H[f_l]
        ! -- Spin polarized calculation
        !    case 1: CONSTRAIN_RHOSPIN is turned on
        !       V_mda,xc(s) = Sum(k,l) <V_xc,s[rho]|f_k> G_kl V_H[f_l]
        !    case 2: CONSTRAIN_RHOSPIN is turned off
        !       V_mda,xc(up,down) = V_mda,xc[tot] +- V_mda,xc[spin]
        !       with
        !       V_mda,xc[t] = Sum(k,l) <V_xc:t[rho]|f_k> G_kl[t] V_H[f_l]
        !       and
        !       V_xc:tot,spin[rho] = ( V_xc,up[rho] +- V_xc,down[rho] ) / 2
        !
        ! In case of a potential extended MDA calculation the right-hand side
        ! <V_xc:t[rho]|f_k> of the equation system has to be replaced by
        ! <V_xc:t[rho]-V_xc,fit:t|f_k> + ...
        ! Sum(s) Sum(m,n) <rho_s-rho_fit,s|g_m> Y_mn (g_n|Z_xc,s:t[rho]|f_k)
        ! ---------------------------------------------------------------------
        if (ext_mda) then
           call start_timer(timer_grid_projections)
           ! add the analytical < V_xc,fit:s | f_k > contribution to rhs(:,s)
           do s=1,ispin
              if (output_chargefit) then
                 write(output_unit,*) '<f_k|V_xc,s> for spin',s
                 write(output_unit,'(10F8.4)') rhs(:,s)-rhs_lc(:,s)
              endif
              do i=1,n_xc
                 rhs(:,s) = rhs(:,s) - mat_xc_ch(i,:)*coeff_xc(i,s)
              end do
              if (output_chargefit) then
                 write(output_unit,*) '<f_k|V_xc,s-V_xc,fit,s> for spin',s
                 write(output_unit,'(10F8.4)') rhs(:,s)-rhs_lc(:,s)
                 write(output_unit,*) "Sum(n,m,s') (f_k|Z_xc,s,s'|g_n) Y_nm <g_m|rho_s'-rho_fit,s'>"
                 write(output_unit,'(10F8.4)') rhs_lc(:,s)
              endif
           end do
           call stop_timer(timer_grid_projections)
        endif
        call start_timer(timer_grid_fit_coeffs)
        ! transform rhs into V_tot and V_spin contributions (if required)

        if (ispin > 1 .and. .not.constrain_rhospin) then
           rhs(:,1) = half * ( rhs(:,1) + rhs(:,2) ) ! V_xc:tot
           rhs(:,2) =          rhs(:,1) - rhs(:,2)   ! V_xc:spin
        endif
        ! The sequence of linsys calls  m u s t  exactly correspond to
        ! the sequence of linsys calls used for the charge density fit, and
        ! the pre-decomposed unperturbted F matrix  m u s t  be used here.
        zero_target = zero
        call solve_linsys(n_ch,matinv_charge,rhs(:,1),dual_charge_norm, &
                          zero_target)
        if (ispin > 1) then
           if (constrain_rhospin) then
              zero_target = zero
              call solve_linsys(n_ch,matinv_charge,rhs(:,2),dual_charge_norm, &
                                zero_target)
           else
              call solve_linsys(n_ch,matinv_charge,rhs(:,2))
           endif
        endif
        ! back-transform coeff_xcmda into V_up and V_down contributions
        if (ispin > 1 .and. .not.constrain_rhospin) then
           coeff_xcmda(:,1) = rhs(:,1) + rhs(:,2) ! V_mda,xc(up)
           coeff_xcmda(:,2) = rhs(:,1) - rhs(:,2) ! V_mda,xc(down)
           coeff_xcmda_ts(:,1:2)=rhs(:,1:2)
        else
           coeff_xcmda = rhs
           if(ispin>1) then
                coeff_xcmda_ts(:,1)=half*(rhs(:,1) + rhs(:,2))
                coeff_xcmda_ts(:,2)=half*(rhs(:,1) - rhs(:,2))
           else
                coeff_xcmda_ts(:,:)=rhs(:,:)
           endif
        endif

!:MF[ for testing
!       if(comp_exact) then
!         if (ispin > 1 .and. .not.constrain_rhospin) then
!            rhs_ex(:,1) = half * ( rhs_ex(:,1) + rhs_ex(:,2) ) ! V_xc:tot
!            rhs_ex(:,2) =          rhs_ex(:,1) - rhs_ex(:,2)   ! V_xc:spin
!         endif
!         zero_target = zero
!         call solve_linsys(n_ch,matinv_charge,rhs_ex(:,1),dual_charge_norm, &
!                         zero_target)
!         if (ispin > 1) then
!             if (constrain_rhospin) then
!                zero_target = zero
!                call solve_linsys(n_ch,matinv_charge,rhs_ex(:,2),dual_charge_norm, &
!                               zero_target)
!             else
!                call solve_linsys(n_ch,matinv_charge,rhs_ex(:,2))
!             endif
!         endif
!         if (ispin > 1 .and. .not.constrain_rhospin) then
!             coeff_xcmda_old(:,1) = rhs_ex(:,1) + rhs_ex(:,2) ! V_mda,xc(up)
!             coeff_xcmda_old(:,2) = rhs_ex(:,1) - rhs_ex(:,2) ! V_mda,xc(down)
!         else
!             coeff_xcmda_old = rhs_ex
!         endif
!       endif
!:MF]

        call stop_timer(timer_grid_fit_coeffs)

!:ANALYZE_MDA[
        if(atomic_radius > 0.0_r8_kind)then
           ! set up the radial grid
           allocate(r1(vec_length),r2(vec_length))
           l = min(vec_length,1001)
           write(output_unit,*)'>>> l = ',l
           dr = atomic_radius / real(l-1,r8_kind)
           r1(1) = 1.0E-16_r8_kind
           do i=2,l
              r1(i) = real(i-1,r8_kind)*dr
           enddo
           r2(:l) = r1(:l)*r1(:l)
           write(output_unit,*)'>>>TST: rad(:)'
           write(output_unit,'(4es19.12)') r1(:l)
           ob => fcts_ch%o(:,:,1)

           ! evaluate the fit functions on the grid
!       print*,' start evaluate the fit functions on the grid'
           ob = 0.0_r8_kind
           i_sa = 0
           do i_bas=-1,min(unique_atoms(1)%lmax_ch,0) ! only s and r2 types
              if(i_bas==-1)uab => unique_atoms(1)%l_ch(0) ! s-type
              if(i_bas== 0)uab => unique_atoms(1)%r2_ch   ! r2-type
              do i_exp=1,uab%N_exponents
                 i_sa = i_sa + 1
                 alpha = uab%exponents(i_exp)
                 if(i_bas==-1)then ! s-type
                    ob(:l,i_sa) = exp( - alpha * r2(:l) )
                 endif
                 if(i_bas== 0)then ! r2-type
                    ob(:l,i_sa) = r2(:l) * exp( - alpha * r2(:l) )
                 endif
              enddo
           enddo
           n_sa = i_sa

           ! evaluate the Vxc(mda) charge density on the grid
           rho = 0.0_r8_kind
           do s=1,ispin
              do i=1,n_ch
                 rho(:l,s) = rho(:l,s) + coeff_xcmda(i,s) * ob(:l,i)
              enddo
              write(output_unit,'(1x,a,i1,a)')'>>>TST: rho_xc(:,',s,')'
              write(output_unit,'(4es19.12)') rho(:l,s)
           enddo

           ! evaluate the fitted charge density on the grid
           rho = 0.0_r8_kind
           if (ispin == 1) then
              do i=1,n_ch
                 coeff = coeff_charge(i)
                 rho(:l,1) = rho(:l,1) + coeff * ob(:l,i)
              enddo
           else
              do i=1,n_ch
                 coeff = half * ( coeff_charge(i) + coeff_spin(i) ) ! spin up
                 rho(:l,1) = rho(:l,1) + coeff * ob(:l,i)
                 coeff = half * ( coeff_charge(i) - coeff_spin(i) ) ! spin down
                 rho(:l,2) = rho(:l,2) + coeff * ob(:l,i)
              enddo
           endif

           ! truncate the fitted charge density on the grid
           do s=1,ispin
              if (a2 /= zero) then
                 rho(:l,s) = half*( rho(:l,s) + sqrt( a2+rho(:l,s)*rho(:l,s) ) )
              endif
              rho(:l,s) = max( rho(:l,s), rho_cutoff )
           enddo

           if (rho_exact) then
              allocate(grdpts(l,3))
              grdpts(:l,1) = unique_atoms(1)%position(1,1) + r1(:l)
              grdpts(:l,2) = unique_atoms(1)%position(2,1)
              grdpts(:l,3) = unique_atoms(1)%position(3,1)
              do i=1,SIZE(orbs_ob)
                 orbs_ob(i)%o = 0.0_r8_kind
              enddo
              call orbital_calculate(grdpts(1:l,1:3),l,orbs_ob=orbs_ob)
              rho = 0.0_r8_kind
              call density_calc(l,rho,orbs_ob)
              deallocate(grdpts)
           endif

           do s=1,ispin
              write(output_unit,'(1x,a,i1,a)')'>>>TST: rho(:,',s,')'
              write(output_unit,'(4es19.12)') rho(:l,s)
           enddo

           ! evaluate the exchange correlation potential  on the grid
           dfdrho = 0.0_r8_kind
           fxc    = 0.0_r8_kind
           if(vwn   ) call vwn_calc   (rho,dfdrho,ispin,fxc,l)
           if(xalpha) call xalpha_calc(rho,dfdrho,ispin,fxc,l)
           do s=1,ispin
              write(output_unit,'(1x,a,i1,a)')'>>>TST: Vxc(:,',s,')'
              write(output_unit,'(4es19.12)') dfdrho(:l,s)
           enddo

           ! evaluate the Coulomb potential of the fit functions on the grid
           ! s-type function
           ! rho(r) = exp(-a*r^2)
           ! phi(r) = sqrt(pi/a)^3 erf(sqrt(a)*r)/r
           ! r^2-type function
           ! rho(r) = r^2 exp(-a*r^2)
           ! phi(r) = 3/2a sqrt(pi/a)^3 erf(sqrt(a)*r)/r - pi/a^2 exp(-a*r^2)
           ob = 0.0_r8_kind
           i_sa = 0
           do i_bas=-1,min(unique_atoms(1)%lmax_ch,0) ! only s and r2 types
              if(i_bas==-1)uab => unique_atoms(1)%l_ch(0) ! s-type
              if(i_bas== 0)uab => unique_atoms(1)%r2_ch   ! r2-type
              do i_exp=1,uab%N_exponents
                 i_sa = i_sa + 1
                 alpha = uab%exponents(i_exp)
                 alpha2 = sqrt(alpha)
                 fact = sqrt(pi/alpha)*(pi/alpha)
                 ! s-type
                 do i=1,l
                   ob(i,i_sa) = fact * erf( alpha2 * r1(i) ) / r1(i)
                 enddo
                 if(i_bas== 0)then ! r2-type
                    fact  = 1.5_r8_kind/alpha
                    fact2 = pi/(alpha*alpha)
                    ob(:l,i_sa) = fact * ob(:l,i_sa) - &
                                  fact2 * exp( - alpha * r2(:l) )
                 endif
              enddo
           enddo

           ! evaluate the fitted Coulomb potential on the grid
           rho = 0.0_r8_kind
           if (ispin == 1) then
              do i=1,n_ch
                 coeff = coeff_charge(i)
                 rho(:l,1) = rho(:l,1) + coeff * ob(:l,i)
              enddo
           else
              do i=1,n_ch
                 coeff = half * ( coeff_charge(i) + coeff_spin(i) ) ! spin up
                 rho(:l,1) = rho(:l,1) + coeff * ob(:l,i)
                 coeff = half * ( coeff_charge(i) - coeff_spin(i) ) ! spin down
                 rho(:l,2) = rho(:l,2) + coeff * ob(:l,i)
              enddo
           endif
           do s=1,ispin
              write(output_unit,'(1x,a,i1,a)')'>>>TST: V_H[rho](:,',s,')'
              write(output_unit,'(4es19.12)') rho(:l,s)
           enddo

           ! evaluate the fitted exchange correlation potential on the grid
           rho = 0.0_r8_kind
           do s=1,ispin
              do i=1,n_ch
                 coeff = coeff_xcmda(i,s)
                 rho(:l,1) = rho(:l,1) + coeff * ob(:l,i)
              enddo
              write(output_unit,'(1x,a,i1,a)')'>>>TST: Vxc_fit(:,',s,')'
              write(output_unit,'(4es19.12)') rho(:l,s)
           enddo

           deallocate(r1,r2)
        endif !atomic radius
!:ANALYZE_MDA]
     else
        ! sending partial numerical integrals to the master
        call xcmda_send_tomaster()
        if (ispin > 1) then
           deallocate(coeff_spin,stat=status)
           if (status /= 0) call error_handler &
                ("BUILD_XCMDA: deallocation of coeff_spin failed")
!:TST[
           write(1,*)'@@ deallocate coeff_spin in build_xcmda'
!:TST]
        endif
        if (eval_lincorr) then
           deallocate(coeff_deltarho,stat=status)
           if (status /= 0) call error_handler &
                ("BUILD_XCMDA: deallocation of coeff_deltarho failed")
!:TST[
           write(1,*)'@@ deallocate coeff_deltarho in build_xcmda'
!:TST]
        endif
     end if
     call xcmda_deallocate(a2/=zero)
     if (matinv_required) then
        matinv_exists = .true.
        matinv_required = .false.
     endif

   end subroutine build_xcmda
!:ANALYZE_MDA[
  function erf(x)
    ! Purpose : calculate the complementary error-function
    !           using the chebyshev approximation
    !
    ! References: YL Luke 1975, pp 123-4
    !             Num. Rec. 1988 "CHEBEV"
    !
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in) :: x
    real(kind=r8_kind)             :: erf
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind),parameter :: na=25,nc=22
    real(kind=r8_kind)              :: a(0:na),c(0:nc)
    real(kind=r8_kind)              :: d,dd,alpha,z,sv
    integer(kind=i4_kind)           :: j
    real(kind=r8_kind),parameter    :: sqpi2_loc=1.128379167095513e0_r8_kind
    real(kind=r8_kind),parameter    :: zero  = 0.0_r8_kind
    real(kind=r8_kind),parameter    :: half  = 0.5_r8_kind
    real(kind=r8_kind),parameter    :: one   = 1.0_r8_kind
    real(kind=r8_kind),parameter    :: two   = 2.0_r8_kind
    real(kind=r8_kind),parameter    :: three = 3.0_r8_kind
    real(kind=r8_kind),parameter    :: four  = 4.0_r8_kind

    data a / &
         .109547129977762e+1_r8_kind, -.289175401126989e+0_r8_kind, &
         .110456398633795e+0_r8_kind, -.412531882278565e-1_r8_kind, &
         .140828380706516e-1_r8_kind, -.432929544743143e-2_r8_kind, &
         .119827190159228e-2_r8_kind, -.299972962353249e-3_r8_kind, &
         .683258603788747e-4_r8_kind, -.142469884548677e-4_r8_kind, &
         .273540877283989e-5_r8_kind, -.048619128719754e-5_r8_kind, &
         .008038727621172e-5_r8_kind, -.001241841831213e-5_r8_kind, &
         .000179953258879e-5_r8_kind, -.000024547948775e-5_r8_kind, &
         .000003162508603e-5_r8_kind, -.000000385902200e-5_r8_kind, &
         .000000044720291e-5_r8_kind, -.000000004933613e-5_r8_kind, &
         .000000000519303e-5_r8_kind, -.000000000052258e-5_r8_kind, &
         .000000000005037e-5_r8_kind, -.000000000000466e-5_r8_kind, &
         .000000000000041e-5_r8_kind, -.000000000000004e-5_r8_kind /
    DATA c / &
         .975083423708556e+0_r8_kind, -.240493938504146e-1_r8_kind, &
         .820452240880432e-3_r8_kind, -.434293081303427e-4_r8_kind, &
         .301844703403493e-5_r8_kind, -.025447331925082e-5_r8_kind, &
         .002485835302051e-5_r8_kind, -.000273172013238e-5_r8_kind, &
         .000033084722797e-5_r8_kind, -.000004350549080e-5_r8_kind, &
         .000000614121457e-5_r8_kind, -.000000092236928e-5_r8_kind, &
         .000000014635665e-5_r8_kind, -.000000002439278e-5_r8_kind, &
         .000000000424976e-5_r8_kind, -.000000000077084e-5_r8_kind, &
         .000000000014507e-5_r8_kind, -.000000000002824e-5_r8_kind, &
         .000000000000567e-5_r8_kind, -.000000000000117e-5_r8_kind, &
         .000000000000025e-5_r8_kind, -.000000000000005e-5_r8_kind, &
         .000000000000001e-5_r8_kind /

    d = zero
    dd = zero
    erf = zero
    if (abs(x).lt.three) then
       !         CALCULATE DIRECTLY
       z=x/three
       alpha=two-four*z*z
       do  J=NA,0,-1
          SV=D
          D=-ALPHA*D-DD+A(J)
          DD=SV
       enddo

       erf=sqpi2_loc*Z*(D-DD)
    else
       !      CALCULATE VIA ERFC
       Z=abs(THREE/X)
       ALPHA=TWO-FOUR*Z*Z
       do J=NC,0,-1
          SV=D
          D=-ALPHA*D-DD+C(J)
          DD=SV
       enddo
       if (X.GT.ZERO) then
          erf=ONE-HALF*exp(-X*X)/X*(D+HALF*ALPHA*DD)*sqpi2_loc
       else
          erf=(-ONE)+HALF*exp(-X*X)/(-X)*(D+HALF*ALPHA*DD)*sqpi2_loc
       endif
    endif

  end function erf
!:ANALYZE_MDA]

   !***************************************************************

!:UB[ modified
   function xcmda_get_exc(get_exc_corr)
!:UB]
     ! purpose : make the module private integrated xc_energy
     !           available to the calling unit
!:UB[ added
     real(kind=r8_kind), optional :: get_exc_corr
!:UB]
     real(kind=r8_kind) :: xcmda_get_exc
     !** End of interface *****************************************
!:UB[ modified
     if (present(get_exc_corr)) then
        xcmda_get_exc = exc_int
        get_exc_corr = exc_corr
     else
        xcmda_get_exc = exc_int + exc_corr
     endif
!:UB]
   end function xcmda_get_exc

   !***************************************************************

   subroutine mda_options_read()
     ! read in mda options from input file
     ! (called by read_input)
     !** End of interface *****************************************
     use input_module
     integer :: unit, status
     real(kind=r8_kind), parameter :: zero = 0.0_r8_kind
     !------------ Executable code --------------------------------
     rho_shape_eps     = df_rho_shape_eps
     rho_cutoff        = df_rho_cutoff
     constrain_rhospin = df_constrain_rhospin
!:ANALYZE_MDA[
     rho_exact         = df_rho_exact
     no_first          = df_no_first
     atomic_radius     = df_atomic_radius
!:ANALYZE_MDA]

     if ( input_line_is_namelist("mda_options") ) then
        call input_read_to_intermediate
        unit = input_intermediate_unit()
        read(unit, nml=mda_options, iostat=status)
        if (status .gt. 0) call input_error( &
             "mda_options_read: namelist occupation.")
        if (rho_shape_eps == zero .and. rho_cutoff <= zero) call input_error( &
             "mda_options_read: non-positive cutoff value encountered")
     endif
   end subroutine mda_options_read

   !***************************************************************

   subroutine mda_options_write(iounit)
     !
     ! Write namelist mda_options to iounit in input format
     ! (called by write_input)
     !
     use echo_input_module, only: start, real, flag, intg, strng, stop, &
          echo_level_full
     use operations_module, only: operations_echo_input_level
     implicit none
     integer, intent(in) :: iounit
     !** End of interface *****************************************


     call start("MDA_OPTIONS","MDA_OPTIONS_WRITE", &
          iounit,operations_echo_input_level)
     call flag("CONSTRAIN_RHOSPIN",constrain_rhospin,df_constrain_rhospin)
!:ANALYZE_MDA[
     call flag("RHO_EXACT        ",rho_exact        ,df_rho_exact        )
     call flag("COMP_EXACT       ",comp_exact       ,df_comp_exact       )
!:ANALYZE_MDA]
     call real("RHO_SHAPE_EPS    ",rho_shape_eps    ,df_rho_shape_eps    )
     call real("RHO_CUTOFF       ",rho_cutoff       ,df_rho_cutoff       )
!:ANALYZE_MDA[
     call real("ATOMIC_RADIUS    ",atomic_radius    ,df_atomic_radius    )
!:ANALYZE_MDA]
     call stop()

   end subroutine mda_options_write

   !****************************************************************************

   subroutine mda_options_bcast()
    !  Purpose: pack the structure mda_options into a comm buffer
    !** End of interface *****************************************
    use comm, only: comm_bcast
    !------------ Executable code ------------------------------------
    call comm_bcast(rho_shape_eps)
    call comm_bcast(rho_cutoff)
    call comm_bcast(constrain_rhospin)
    call comm_bcast(rho_exact)
    call comm_bcast(comp_exact)
    call comm_bcast(no_first)
   end subroutine mda_options_bcast

   !****************************************************************************

  real(kind=r8_kind) function mda_rho_shape_eps()
    !  Purpose: make private control variable rho_shape_eps accessible
    !           UB
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************

    !------------ Executable code ------------------------------------
    mda_rho_shape_eps = rho_shape_eps

  end function mda_rho_shape_eps

  !***************************************************************

  real(kind=r8_kind) function mda_rho_cutoff()
    !  Purpose: make private control variable rho_cutoff accessible
    !           UB
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************

    !------------ Executable code ------------------------------------
    mda_rho_cutoff = rho_cutoff

  end function mda_rho_cutoff

  !***************************************************************

  logical function mda_constrain_rhospin()
    !  Purpose: make private control variable constran_rhospin accessible
    !           UB
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************

    !------------ Executable code ------------------------------------
    mda_constrain_rhospin = constrain_rhospin

  end function mda_constrain_rhospin

  !*************************************************************

  subroutine xcmda_coeff_store(th)
    ! Purpose: stores the Coulomb-type XC coefficients on a
    !          readwriteblocked file in case the input switch
    !          "save_scfstate" is set.
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB, 8/97
    !------------ Modules ----------------------------------------
    use options_module
    use readwriteblocked_module
    use output_module, only: output_data_saved, output_main_scf
    !------------ Declaration of formal paramaters ---------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    !------------ Declaration of local variables   ---------------
    integer(kind=i4_kind) :: s, n_spin
    real(kind=r8_kind) :: yes(1) = (/1.0_r8_kind/), no(1) = (/0.0_r8_kind/)
    !------------ Executable code ------------------------------------
    n_spin = options_n_spin()

    if (output_main_scf) call write_to_output_units &
         ("XCMDA_COEFF_STORE: saving XC(MDA) coefficients")
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored XC(MDA) coefficients :'
    endif

    if ( options_xcmode() /= xcmode_model_density .and. &
         options_xcmode() /= xcmode_extended_mda ) then
       call readwriteblocked_write(no,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'model density used ?'
          write(output_unit,'(4es20.13)')no
       endif
    else
       call readwriteblocked_write(yes,th)
       call readwriteblocked_write((/real(n_spin,r8_kind)/),th)
       do s=1,n_spin
          call readwriteblocked_write(coeff_xcmda(:,s),th)
       enddo
       call readwriteblocked_write((/exc_int/),th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'model density used ?'
          write(output_unit,'(4es20.13)')yes
          write(output_unit,'(  a     )')'n_spin'
          write(output_unit,'(4es20.13)')(/real(n_spin,r8_kind)/)
          do s=1,n_spin
             write(output_unit,'( a,i1,a )')'coeff_xcmda(:,',s,')'
             write(output_unit,'(4es20.13)')coeff_xcmda(:,s)
          enddo
          write(output_unit,'(  a     )')'exc_int'
          write(output_unit,'(4es20.13)')(/exc_int/)
       endif
    endif

  end subroutine xcmda_coeff_store

!*************************************************************
! record  1: model_density
! record  2: n_spin               [if model_density]
! record  3: coeff_xcmda(spin=1)  [if model_density]
! record  4: coeff_xcmda(spin=2)  [if model_density & n_spin > 1]
! record  5: exc_int
!*************************************************************

  subroutine xcmda_coeff_recover(th)
    ! Purpose: recovers the Coulomb-type XC coefficients from a
    !          readwriteblocked file in case the input switch
    !          "read_scfstate" is set.
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB 8/97
    !------------ Modules ----------------------------------------
    use options_module
    use readwriteblocked_module
    use output_module, only: output_main_scf, output_data_read
    !------------ Declaration of formal_parameters ---------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    !------------ Declaration of local variables   ---------------
    allocatable           :: buffer
    integer(kind=i4_kind) :: n_spin, spin_stored
    real(kind=r8_kind)    :: dummy(1), spin(1), use_mda(1), buffer(:), &
                             zero = 0.0_r8_kind, half = 0.5_r8_kind
    logical               :: use_model_density
    integer               :: status
    external error_handler
    !------------ Executable code ------------------------------------
    n_spin = options_n_spin()
    use_model_density = options_xcmode() == xcmode_model_density .or. &
                        options_xcmode() == xcmode_extended_mda

    if (output_main_scf) call write_to_output_units &
         ("XCMDA_COEFF_RECOVER: reading XC(MDA) coefficients")

    call readwriteblocked_read(use_mda,th)
    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered XC(MDA) coefficients :'
       write(output_unit,'(  a     )')'model density used ?'
       write(output_unit,'(4es20.13)')use_mda(1)
    endif
    if (use_mda(1) /= zero) then
       call readwriteblocked_read(spin,th)
       spin_stored = int(spin(1),i4_kind)
       if (output_data_read) then
          write(output_unit,'(  a     )')'n_spin'
          write(output_unit,'(4es20.13)')spin(1)
       endif
       if (use_model_density) then
          if (output_main_scf) call write_to_output_units &
               ("XCMDA_COEFF_RECOVER: reading XC(MDA) coefficients")
       else
          if (output_main_scf) call write_to_output_units &
               ("XCMDA_COEFF_RECOVER: XC(MDA) coefficients ignored")
          call readwriteblocked_skipread(n_ch*n_spin+1,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xcmda skipped'
             write(output_unit,'(  a     )')'exc_int     skipped'
          endif
          return
       endif
    else
       if (use_model_density) then
          if (output_main_scf) call write_to_output_units &
               ("XCMDA_COEFF_RECOVER: XC(MDA) coefficients re-computed")
          call start_timer(timer_scf_xc)
          ABORT("SPMD?")
          call xcmda_build(.false.)
          call stop_timer(timer_scf_xc)
       endif
       return
    endif

    ! start recovering XC(MDA) coefficients now
    if (spin_stored > n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("XCMDA_COEFF_RECOVER: trying to convert spin-polarized data from")
       call write_to_output_units &
            ("                    tape into the spin-restricted data required")
       ! coeff_xcmda(:) := 1/2 * sum(s) coeff_xcmda(:,s)
    endif
    if (spin_stored < n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("XCMDA_COEFF_RECOVER: trying to convert spin-restricted data from")
       call write_to_output_units &
            ("                    tape into the spin-polarized data required")
       ! coeff_xcmda(:,s) := coeff_xcmda(:)
    endif

    call readwriteblocked_read(coeff_xcmda(:,1),th)
    if (output_data_read) then
       write(output_unit,'(  a     )')'coeff_xcmda(:,1)'
       write(output_unit,'(4es20.13)')coeff_xcmda(:,1)
    endif
    if (spin_stored > 1) then
       if (n_spin > 1) then
          call readwriteblocked_read(coeff_xcmda(:,2),th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'>>>STO: coeff_xcmda(:,2)'
             write(output_unit,'(4es20.13)')coeff_xcmda(:,2)
          endif
       else
          allocate(buffer(n_ch),stat=status)
          if (status.ne.0) call error_handler &
               ("XCMDA_COEFF_RECOVER: allocation of buffer failed")
          call readwriteblocked_read(buffer,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'>>>STO: coeff_xcmda(:,2)'
             write(output_unit,'(4es20.13)')buffer
          endif
          coeff_xcmda(:,1) = ( coeff_xcmda(:,1) + buffer ) * half
          deallocate(buffer,stat=status)
          if (status.ne.0) call error_handler &
               ("XCMDA_COEFF_RECOVER: deallocation of buffer failed")
       endif
    elseif (n_spin > 1) then
       coeff_xcmda(:,2) = coeff_xcmda(:,1)
    endif
    call readwriteblocked_read(dummy(1:1),th)
    exc_int = dummy(1)
    if (output_data_read) then
       write(output_unit,'(  a     )')'exc_int'
       write(output_unit,'(4es20.13)')dummy(1)
    endif

  end subroutine xcmda_coeff_recover

  !************************************************************
 end module xcmda_hamiltonian

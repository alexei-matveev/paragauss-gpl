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
module post_scf_module
  !
  !  Purpose: main module for calculation of the xc-energy
  !           on the post_scf_grid
  !           Every processor treats only one part of the grid.
  !           After the calculation the partial results are sent to the
  !           master.
  !           The structure of the module is similiar to the xc_module
  !           with two main differences:
  !             - the becke weights of the post_scf grid are multiplied
  !               in this module
  !             - the integration of matrix elements is not performed as it
  !               not nexessary
  !
  ! Module called by: main_master,main_scf
  !
  !
  !  Author: MS
  !  Date: 9/96
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   3/97
  ! Description: Adding routines for energy gradients
  !         post_scf_calc_xc_en_and_gr
  !         post_scf_calc_xc_en_and_gr_nl for GGA
  !         formulas for gradient terms can be found in
  !         Pople CPL 199 (1992) 557
  !
  ! Modification (Please copy before editing)
  ! Author: AM
  ! Date:   9/98
  ! Description: PBE functionals introduced
  !
  ! Modification
  ! Author: TS
  ! Date:   07/09
  ! Description: extensions for MGGA functionals
  !
  ! Modification (Please copy before editing)
  ! Author: AN
  ! Date:   10/10
  ! Description: use dynamic load balancing (dlb) for job distribution
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
! define FPP_TIMERS 3
# include "def.h"
  USE_MPI, only: MPI_WTIME
  use strings, only: itoa
  use symmetry_data_module ! description of irreps
  use type_module  ! type specification parameters
  use grid_module, only: atomicweight,grid_main,&
       grid_close,atomicweight_and_grad,weight_grads,atomicweight_and_dervs
  use grid_module, only: more_grid_atom, grid_loop_setup_atom
  use orbitalstore_module
  use orbital_module
  use density_data_module
  use occupied_levels_module, only: sndrcv_eigvec_occ
  use fit_coeff_module, only: fit_coeff_sndrcv, coeff_charge, coeff_spin, &
                              fit_coeff_shutdown, &
                              coeff_xcmda, fit_coeff_n_ch, coeff_xcmda_ts, ICHFIT
  use density_calc_module
  use unique_atom_module
  use machineparameters_module
  use spin_orbit_module
  use output_module
  use operations_module, only: operations_gradients, operations_response
  use options_module, only: options_split_gradients, options_xcmode, &
                            xcmode_model_density, xcmode_extended_mda, &
                            options_spin_orbit
  use xcmda_hamiltonian, only: mda_rho_shape_eps, mda_constrain_rhospin,comp_exact
  use datatype, only: arrmat2,arrmat4,arrmat5,arrmat6,arrmat3,arrmat7
  USE_MEMLOG
  use comm, only: comm_rank, comm_size
#ifdef FPP_TIMERS
  use comm, only: comm_reduce
#endif
#ifdef FPP_DEBUG
  use error_module, only: MyID
#endif
  implicit none
  private
  save

  !== Interrupt end of public interface of module =====================

  integer(i4_kind), parameter, private ::&
       & X=1, Y=2, Z=3,&
       & UP=1, DN=2,&
       & UPUP=1, DNDN=2, UPDN=3, &
       & XX = 1, &
       & XY = 2, &
       & XZ = 3, &
       & YY = 4, &
       & YZ = 5, &
       & ZZ = 6

  real(r8_kind), parameter :: ZERO = 0.0_r8_kind
  real(r8_kind), parameter :: rho_cutoff=1.0E-08_r8_kind

  !------------ public functions and subroutines ------------------
  public ::&
       & post_scf_main, &
       & post_scf_deallocate_grad_xc

#ifdef WITH_SECDER
  public :: cpks_xc_main
#endif


  !================================================================
  ! End of public interface of module
  !================================================================

!..............................................................................
! << OUTPUT ARRAYS >>
! ===================
!
! options_xcmode() /= xcmode_model_density
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! F1 = Sum(s) < V_xc,s[rho_up,rho_dn] | -d/dR rho_s >_num
!    = - 2 Sum(s,n) o_n < d/dR psi_n,s | V_xc,s[rho_up,rho_dn] | psi_n,s >
!    = - 2 Sum(i,j,s) P_ij,s < d/dR xi_i | V_xc,s[rho_up,rho_dn] | xi_j >
!
! F2 = - Sum(a) [ d/dR w_a ] e_xc[rho_up,rho_dn](r_a)
!      - Sum(a) w_a [ d/dR r_a ] Nabla e_xc[rho_up,rho_dn](r_a)
!
! options_xcmode() == xcmode_model_density
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! F1 = Sum(s) < V_xc,s[rhofit_up,rhofit_dn] | -d/dR rhofit_s >_num
!    = Sum(k,s) a_k,s < V_xc,s[rhofit_up,rhofit_dn] | -d/dR rho_k >
!    = Sum(k,t) a_k,t < V_xc,t[rhofit_up,rhofit_dn] | -d/dR rho_k >
!      with t = tot,spin
!
! F2 = - Sum(a) [ d/dR w_a ] e_xc[rhofit_up,rhofit_dn](r_a)
!      - Sum(a) w_a [ d/dR r_a ] Nabla e_xc[rhofit_up,rhofit_dn](r_a)
!
! Storage (negative d/dR E_tot contributions)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! options_split_gradients()       | FALSE  FALSE  TRUE   TRUE
! weight_grads                    | FALSE  TRUE   FALSE  TRUE
! --------------------------------+--------------------------
! grad_xc  (1:n_ua)%m(1:3,1:n_ea) | F1     F1+F2  F1     F1
! grad_grid(1:n_ua)%m(1:3,1:n_ea) | --     --     --     F2
!
! IF( model_density ) coeff_xcmda(1:n_ch,1:n_spin) := b_k,s
!
! << WORKING ARRAYS >>
! ====================
! grdpts         : r_a
! grdwts         : w_a
! graw           : d/dR w_a
!
! orbs_ob        :              xi_i(r_a)
! orbs_grads     :   d/dr       xi_i(r_a)  =  Sum(R) [-d/dR]      xi_i(r_a)
! nuc_grads      : [-d/dR]      xi_i(r_a)
! nuc_sec_der    : [-d/dR] d/dr xi_i(r_a)
!                    d/dr  d/dr xi_i(r_a)  =  Sum(R) [-d/dR] d/dr xi_i(r_a)
!
! fcts_ch        :              f_k(r_a)
! grads_ch       :   d/dr       f_k(r_a)  =  Sum(R) [-d/dR]      f_k(r_a)
! sec_ders_ch    :   d/dr  d/dr f_k(r_a)  =  Sum(R) [-d/dR] d/dr f_k(r_a)
! nuc_grads_ch   : [-d/dR]      f_k(r_a)
! sec_nuc_ders_ch: [-d/dR] d/dr f_k(r_a)
!
! rho            :              rho_s(r_a)
! grarho[xyz]    :   d/dr       rho_s(r_a)  =  Sum(R) [-d/dR]      rho_s(r_a)
! grad_rho       :   d/dr       rho_s(r_a)  =  Sum(R) [-d/dR]      rho_s(r_a)
! sec_der_rho    :   d/dr  d/dr rho_s(r_a)  =  Sum(R) [-d/dR] d/dr rho_s(r_a)
! nuc_grarho     : [-d/dR]      rho_s(r_a)
! nuc_grad_tau   : [-d/dR]      tau_s(r_a)
! nuc_sec_derrho : [-d/dR] d/dr rho_s(r_a)
! gamma          : < d/dr rho_s | d/dr rho_t > with (s,t) = (1,1),(2,2),(1,2)
!
! exc_xxx        : the various contributions to E_xc[...]
! fxc            : e_xc   [...](r_a)
! dfdrho         : V_xc,s [...](r_a) = d/drho_s    e_xc[...](r_a)
! dfdgrarho      : G_xc,st[...](r_a) = d/dgamma_st e_xc[...](r_a)
!
! rhs            : <V_xc,s[rho]|f_k> to evaluate the MDA coefficients b_k,s
!..............................................................................

  type(arrmat4),pointer :: graphi(:,:)
  real(kind=r8_kind),allocatable :: dvdrho(:,:) ! cpksdervs calculations
  real(kind=r8_kind),allocatable :: df_drhodgamma(:,:),df_dgammadrho(:,:)
  real(kind=r8_kind),allocatable :: df_dgammadgamma(:,:)

  real(kind=r8_kind),allocatable :: rho(:,:),&
       dfdrho(:,:),&  ! derivative of f with respect to rho
       fxc(:),&       ! funktion f according to JPG
       rhs(:,:)       ! XC-potential projections onto the fitting functions
#ifdef WITH_CORE_DENS
  real(kind=r8_kind),allocatable :: &
       rho_core(:,:), grad_rho_core(:,:,:), sec_grad_core(:,:,:,:)
#endif
  ! in case of spin polarized spin orbit we need
  real(kind=r8_kind),allocatable :: rho_ud(:,:),s_abs(:),s_dens(:,:),&
       & grarho_ud(:,:,:),& !(vl,X:Z,UP:DN) x(:,:),grarho_udy(:,:),grarho_udz(:,:),&
       & gras_dens(:,:) !(vl,X:Z) x(:),gras_densy(:),gras_densz(:)

  real(kind=r8_kind),allocatable :: rhs_intermed(:)

  ! FIXME: move phi_ob to type(arrmat4)
  type(orbital_spin_type),pointer :: phi_ob(:)    ! (n_irr) Molecular Orbitals (eigenfunctions) used in cpks
  type(arrmat5), pointer     :: phi_gr_nuc(:)     ! (n_irr) Molecular Orbitals Gradients wrt nuclear displacments
  real(r8_kind), allocatable :: rho_gr_nuc(:,:,:) ! (:vl,n_modes,ispin) Density Gradients wrt nuclear displacments
  real(r8_kind), allocatable :: wts_gr_nuc(:,:)   ! (:vl,n_modes)       Weight  Gradients wrt nuclear displacments
  real(r8_kind), allocatable :: rho_gr_imp(:,:,:) ! (:vl,n_modes,ispin) Density Impl-Grads wrt nuclear displacments
  real(r8_kind), allocatable :: eny_sd_imp(:,:,:) ! (n_modes,n_modes,ispin) Energy Impl-Dervs wrt nuclear displacments

  type(orbital_type),pointer :: orbs_ob(:)
  type(spinor_type),pointer  :: orbs_spinor_ob(:)
  type(orbital_gradient_type),pointer :: orbs_grads(:)
  type(spinor_gradient_type),pointer :: spinor_grads(:)
  type(orbital_nuclear_gradient_type),pointer :: nuc_grads(:,:)
  type(orbital_nuclear_sec_der_type),pointer :: nuc_sec_der(:,:)
  type(orbital_nuclear_sec_der_type),pointer :: nuc_3rd_der(:,:)

  type(orbital_type)                           :: fcts_ch
  type(orbital_gradient_type)                  :: grads_ch
  type(orbital_sec_der_type)                   :: sec_ders_ch
  type(orbital_nuclear_gradient_type), pointer :: nuc_grads_ch(:)
  type(orbital_nuclear_sec_der_type) , pointer :: sec_nuc_ders_ch(:)
#ifdef WITH_CORE_DENS
  type(core_orbital_type)                      :: orbs_ob_core
  type(core_orbital_gradient_type)             :: grads_core
  type(core_orbital_sec_der_type)              :: sec_der_core
  type(core_orbital_nuc_gradient_type),pointer :: nuc_grads_core(:)
  type(core_orbital_nuc_sec_der_type),pointer  :: nuc_sec_der_core(:)
  type(arrmat4), pointer                       :: nuc_grarho_core(:)
  type(arrmat5), pointer                       :: nuc_sec_derrho_core(:)
#endif
  ! nuc_sec_der(n_unique_atoms,n_irreps)
  ! this structure keeps the derivatives of the primitive orbitals
  ! with respect to nuclear displacements

  type(arrmat6),allocatable :: nuc_dervsrho(:,:)
  type(arrmat7),allocatable :: nuc_dervs_grarho(:,:)
  type(arrmat4),allocatable, target :: nuc_grarho(:)
  type(arrmat4),allocatable, target :: nuc_grad_tau(:)
  type(arrmat4),allocatable :: nuc_grarho_imp(:)
  type(arrmat5),allocatable :: nuc_gragamma_imp(:)
  type(arrmat5),allocatable :: dervsrho_imp(:,:)
  type(arrmat5),allocatable :: dervs_grarho_imp(:,:)
  ! nuc_grarho(n_unique_atoms)
  ! stores gradient of the density with respect to all nuclear
  ! displacements
  type(arrmat5),allocatable, target :: nuc_sec_derrho(:)

  type(arrmat5),allocatable :: dervsw(:,:) ! contains the derivatives of
  type(arrmat3),allocatable :: graw(:) ! contains the derivatives of
! type(arrmat2),pointer :: gradsw(:) ! contains the derivatives of
  ! of the integraton weights with respect to nuclear displeacements

  type(arrmat4), allocatable, public :: dervs_xc(:,:) ! modified from main_gradient()
  type(arrmat4), allocatable, private :: dervs_xc_tmp(:,:)
  type(arrmat2), allocatable, public, protected, target :: grad_xc(:)
  ! grad_xc(n_unique_atoms)
  ! contains xc-contribution to final energy gradient for all atoms
  type(arrmat2), allocatable, public, protected, target :: grad_grid(:)
  ! grad_grid(n_unique_atoms)
  ! contains the energy gradient contribution which arise from the finite
  ! accuracy of the nuclear position dependent integration nodes and weights
  real(kind=r8_kind),allocatable :: charge_int(:,:) ! integrated atomic charges
! MAKE THEM LOCAL: real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths

  logical :: do_grads,weight_grad,model_density,split_gradients
  ! steering parameters for control of exchange-correlation functionals
  integer(kind=i4_kind) :: ispin,n_equal_max,n_ch
  integer(kind=i4_kind) :: vec_length ,vec_length_act !uncommented vec_length_act (=> MDA grads.)
  integer(kind=i4_kind) :: charge_types
  real(kind=r8_kind)    :: rho_shape_eps
  logical               :: constrain_rhospin

  integer, parameter :: POST_SCF_CALL = 20120702 ! a date

FPP_TIMER_DECL(whole)
FPP_TIMER_DECL(loop)
FPP_TIMER_DECL(work)
FPP_TIMER_DECL(sched)

contains

  subroutine post_scf_main()
    ! Purpose : main routine for doing the post scf part
    use symmetry_data_module ! description of irreps
    use timer_module
    use time_module
    use unique_atom_methods, only: unique_atom_grad_information
    use iounitadmin_module
    use xc_cntrl, only: is_on, xc_mgga
    use ph_cntrl, only: nl_calc_ph
    use pointcharge_module
    use xc_func, only: xc_func_reset
    use integralpar_module, only: integralpar_2dervs,integralpar_cpksdervs
    use cpksdervs_matrices, only: cpksalloc
    use filename_module
    use eigen_data_module, only: eigen_data_bcast
#ifdef FPP_TIMERS
    use dlb, only: dlb_print_statistics
#endif
    !** End of interface **************************************

    integer(kind=i4_kind) :: alloc_stat,dummy,i,ua,n_ea
    integer(kind=i4_kind) :: i_ua,j_ua,ii
#ifdef FPP_TIMERS
    real(kind=r8_kind), allocatable :: times(:, :)
    integer(kind=i4_kind) :: ierr
    integer(kind=i4_kind) :: np, rank
#endif

    ! Write to trace unit to show progress of calculation:
    call dtrace ("Entering PostSCF part")

    dummy=1
    do_grads=operations_gradients
    split_gradients=options_split_gradients()
    model_density=options_xcmode()==xcmode_model_density .or. &
                  options_xcmode()==xcmode_extended_mda
    if (model_density) then
    endif
    n_equal_max=maxval(unique_atoms(:)%n_equal_atoms) ! maximum number of equal_atoms
    if(output_post_scf_main) &
         call write_to_output_units("POST_SCF_MAIN: start")
    if(do_grads) then
       if(output_post_scf_main) &
            call write_to_output_units("POST_SCF_MAIN: unique_atom_grad_information")
       ! set important informations for symmetrisation of gradients
       !MOVED to unique_atom_setup: call unique_atom_grad_information()
       if(n_moving_unique_timps.ne.0 .or. &
            (moving_pc.and.pointcharge_n> 0)) call unique_timp_grad_information()
    end if


    ! ===================================================================
    ! === FROM THIS POINT MASTER AND SLAVES RUN IN A PARALLEL CONTEXT ===
    ! ===================================================================

    if (comm_rank() == 0) call start_timer(timer_gridph)
    ! first send eigvecs and fitcoeffs to the slaves
    call sndrcv_eigvec_occ()
    if( integralpar_2dervs )then
      if(output_post_scf_main) &
         call write_to_output_units("POST_SCF_MAIN: send/recv eigenvectors")
      !
      ! send/recv eigenvectors ``eigvec'' (the whole bunch including virtual),
      ! some older code uses the structure ``eigvec_occ'' that does not
      ! include virtual orbitals
      !
      call eigen_data_bcast()
    endif

    ! sometimes there is also need of the fit_coeff
    if (model_density) then
      rho_shape_eps=mda_rho_shape_eps()
      constrain_rhospin=mda_constrain_rhospin()
      n_ch = fit_coeff_n_ch()
      if (comm_size() > 1) then
        if (output_post_scf_main) &
          call write_to_output_units &
          ("POST_SCF_MAIN: fit_coeff_sndrcv")
        ! send coeff_charge and coeff_spin to evaluate rho_fit on the grid
        ! coeff_xcmda has not yet been updated; they are sent later
        call fit_coeff_sndrcv(ICHFIT)
      endif
    endif !model_density


    ! write to trace unit to show progress of calculation
    call dtrace("Calculate PostSCF Exchange Energy")
    ! now prepare integration grid
    if (output_post_scf_main) &
         call write_to_output_units("POST_SCF_MAIN: grid_main")

    call grid_main(post_scf=.true.)

    if(output_post_scf_main) &
         call write_to_output_units("POST_SCF_MAIN: allocations")

    ! set up the density calculations
    ispin=ssym%n_spin
    vec_length=machineparameters_veclen

   if(integralpar_2dervs) then
    ASSERT(weight_grads)
   endif

    mda: if (model_density) then
       call fitted_density_calc_setup()
       if (pseudopot_present.and.core_density_setup) then
          ! BOMB HERE ->
          call error_handler("POST_SCF_MAIN: I am afraid PP+CoreDens needs some clean up")
#ifdef WITH_CORE_DENS
          call fitted_core_dens_calc_setup()
#endif
       end if
       ! allocation of orbitals and their derivatives for different cases
       nonlda: if (nl_calc_ph) then ! GGA calculation
          if (do_grads) then ! GGA gradients

             if(integralpar_2dervs) then
             call orbital_setup(vec_length,do_gradients=.true., &
                                           do_sec_der  =.true., &
                                           do_3rd_der  =.true.  )
             else
             call orbital_setup(vec_length,do_gradients=.true., &
                                           do_sec_der  =.true.)
             endif

             if (weight_grads) then ! d^2/dr^2 rho required
                call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch,&
                                           sec_ders=sec_ders_ch,&
                                           nuc_grads=nuc_grads_ch,&
                                           sec_nuc_ders=sec_nuc_ders_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup ) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           grads_core=grads_core, &
                                           sec_der_core=sec_der_core, &
                                           nuc_grads_core=nuc_grads_core,&
                                           nuc_sec_der_core=nuc_sec_der_core)
                end if
#endif
             else ! d^2/dr^2 rho not required
                call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch,&
                                           nuc_grads=nuc_grads_ch,&
                                           sec_nuc_ders=sec_nuc_ders_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           grads_core=grads_core, &
                                           nuc_grads_core=nuc_grads_core,&
                                           nuc_sec_der_core=nuc_sec_der_core)
                end if
#endif
             endif
          else if(comp_exact) then
             call orbital_setup(vec_length,do_gradients=.true.)
             call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch)
#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup ) &
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           grads_core=grads_core)
#endif
          else ! GGA energy only
             call orbital_setup(vec_length,do_gradients=.true.)
             call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           grads_core=grads_core)
                endif
#endif
          end if
       else nonlda ! i.e. lda
          if (do_grads) then ! LDA gradients
             call orbital_setup(vec_length,do_gradients=.true.)
             if (weight_grads) then ! d/dr rho required
                call fit_fct_allocate("ch",fcts=fcts_ch,grads=grads_ch,&
                                           nuc_grads=nuc_grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           grads_core=grads_core, &
                                           nuc_grads_core=nuc_grads_core)
                end if
#endif
             else ! d/dr rho not required
                call fit_fct_allocate("ch",fcts=fcts_ch,nuc_grads=nuc_grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                           nuc_grads_core=nuc_grads_core)
                end if
#endif
             endif
          else if(comp_exact) then
             call orbital_setup(vec_length)
             call fit_fct_allocate("ch",fcts=fcts_ch)
#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup ) &
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core)
#endif
          else ! LDA energy only
             call orbital_setup(vec_length)
             call fit_fct_allocate("ch",fcts=fcts_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_allocate(orbs_ob_core=orbs_ob_core)
                end if
#endif
          end if
       endif nonlda
       charge_types = 2*ispin

    else  mda !i.e. standard non mda procedure
       call density_calc_setup()

#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
          call fitted_core_dens_calc_setup()
       end if
#endif

       ! allocation of orbitals and their derivatives for different cases
       lda:if(.not.nl_calc_ph) then
          lda_nograds:if(.not.do_grads) then
             call orbital_setup(vec_length)

             if (options_spin_orbit) then
                call orbital_allocate(orbs_spinor_ob=orbs_spinor_ob) ! (1)
             else
                call orbital_allocate(orbs_ob) ! (2) LDA
             endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_allocate(orbs_ob_core=orbs_ob_core)
             endif
#endif

          else lda_nograds !i.e lda_grads

#ifdef WITH_SECDER
           if(integralpar_2dervs) then
             call secder_allocate(20)
           else
#endif
             call orbital_setup(vec_length,do_gradients=.true.)
             call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads) ! (3) LDA GRAD

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                        grads_core=grads_core,&
                                        nuc_grads_core=nuc_grads_core)
             endif
#endif
#ifdef WITH_SECDER
          endif ! if(integralpar_2dervs)
#endif

          endif lda_nograds

       else lda !i.e. not_lda
          no_grad: if(.not.do_grads) then
             call orbital_setup(vec_length,do_gradients=.true.)

             if (options_spin_orbit) then
                call orbital_allocate(orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
             else
                call orbital_allocate(orbs_ob,orbs_grads) ! (4)
             endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                        grads_core=grads_core)
             end if
#endif

          else no_grad ! i.e. calc grads

             if(integralpar_2dervs) then
             call orbital_setup(vec_length,do_gradients=.true.,do_sec_der=.true., &
                                                               do_3rd_der=.true.)
             call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads,& ! (5)
                  nuc_sec_der=nuc_sec_der,nuc_3rd_der=nuc_3rd_der)
             else
             call orbital_setup(vec_length,do_gradients=.true.,do_sec_der=.true.)
             call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads,& ! (5)
                  nuc_sec_der=nuc_sec_der)
             endif


#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                         grads_core=grads_core,&
                                         nuc_grads_core=nuc_grads_core, &
                                         nuc_sec_der_core=nuc_sec_der_core)
             end if
#endif

          endif no_grad
       endif lda
       charge_types = 1
    endif mda ! or  else

    so: if (options_spin_orbit) then

       sop: if (spin_orbit_polarized) then
          !
          ! OPEN SHELL
          !
          allocate(dfdrho(vec_length,2),s_dens(vec_length,3),s_abs(vec_length),rho_ud(vec_length,2), &
                   stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          if (nl_calc_ph) then
!!$             allocate(gras_densx(vec_length),gras_densy(vec_length),gras_densz(vec_length),stat=alloc_stat)
             allocate(gras_dens(vec_length,3),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
!!$             allocate(grarho_udx(vec_length,2),grarho_udy(vec_length,2),grarho_udz(vec_length,2),stat=alloc_stat)
             allocate(grarho_ud(vec_length,3,2),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
          endif
        else sop
          !
          ! CLOSED SHELL
          !
          allocate(dfdrho(vec_length,ispin),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        endif sop
     else so ! i.e. regular run
       allocate(dfdrho(vec_length,ispin),stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
     endif so

    allocate(rho(vec_length,ispin),&
         fxc(vec_length),charge_int(n_unique_atoms,charge_types),&
         stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

#ifdef WITH_CORE_DENS
    if (pseudopot_present.and.core_density_setup) then
      allocate(rho_core(vec_length,ispin), grad_rho_core(vec_length,3,ispin),&
              stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
    end if
#endif

    grads: if(do_grads) then
       allocate(nuc_grarho(n_unique_atoms),stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       if(is_on(xc_mgga)) then
          allocate(nuc_grad_tau(n_unique_atoms),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
       end if
    if(integralpar_2dervs) then
     allocate(nuc_dervsrho(n_unique_atoms,n_unique_atoms),stat=cpksalloc(86))
       ASSERT(cpksalloc(86).eq.0)
       MEMLOG(size(nuc_dervsrho))
       if(nl_calc_ph) then
       allocate(nuc_dervs_grarho(n_unique_atoms,n_unique_atoms),stat=cpksalloc(148))
       ASSERT(cpksalloc(148).eq.0)
       MEMLOG(size(nuc_dervs_grarho))
       endif
    endif

#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
         allocate(nuc_grarho_core(n_unique_atoms),stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
       end if
#endif

       do i=1,n_unique_atoms
          n_ea=unique_atoms(i)%n_equal_atoms
          allocate(nuc_grarho(i)%m(vec_length,3,n_ea,ispin),stat=cpksalloc(133))
          ASSERT(cpksalloc(133).eq.0)

          if(is_on(xc_mgga)) then
             allocate(nuc_grad_tau(i)%m(vec_length,3,n_ea,ispin),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
          end if

        if(integralpar_2dervs) then
         do ii=1,n_unique_atoms
          allocate(nuc_dervsrho(i,ii)%m(vec_length,3,n_ea, &
                                      3,unique_atoms(ii)%n_equal_atoms,ispin), &
                   stat=cpksalloc(87))
         ASSERT(cpksalloc(87).eq.0)
         MEMLOG(size(nuc_dervsrho(i,ii)%m))
         if(nl_calc_ph) then
          allocate(nuc_dervs_grarho(i,ii)%m(vec_length,3,n_ea, &
                   3,unique_atoms(ii)%n_equal_atoms,3,ispin), &
                   stat=cpksalloc(149))
          ASSERT(cpksalloc(149).eq.0)
          MEMLOG(size(nuc_dervs_grarho(i,ii)%m))
          nuc_dervs_grarho(i,ii)%m=0.0_r8_kind
         endif
         enddo
        endif

#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             allocate(nuc_grarho_core(i)%m(vec_length,3,n_ea,ispin),&
                      stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          end if
#endif
       end do

       allocate(grad_xc(N_moving_unique_atoms),stat=cpksalloc(21))
       ASSERT(cpksalloc(21).eq.0)
!       MEMLOG(size(grad_xc))

       init_dervs_xc_uniques: if(integralpar_2dervs)  then
        allocate( dervs_xc(N_moving_unique_atoms,N_moving_unique_atoms), &
                  dervs_xc_tmp(N_moving_unique_atoms,N_moving_unique_atoms), &
                  stat=cpksalloc(80))
        ASSERT(cpksalloc(80).eq.0)
        MEMLOG(size(dervs_xc)+size(dervs_xc_tmp))
       endif init_dervs_xc_uniques

       do i=1,N_moving_unique_atoms
          n_ea=unique_atoms(moving_unique_atom_index(i))%n_equal_atoms

          allocate(grad_xc(i)%m(3,n_ea),stat=cpksalloc(22)) !post_scf_main
          ASSERT(cpksalloc(22).eq.0)
          MEMLOG(size(grad_xc(i)%m))

        init_m_dervs_xc: if(integralpar_2dervs) then
          do ii=1,N_moving_unique_atoms
          allocate( dervs_xc(i,ii)%m(3,n_ea,3, &
                    unique_atoms(moving_unique_atom_index(ii))%n_equal_atoms), &
                    dervs_xc_tmp(i,ii)%m(3,n_ea,3, &
                    unique_atoms(moving_unique_atom_index(ii))%n_equal_atoms), &
                    stat=cpksalloc(81))
          ASSERT(cpksalloc(81).eq.0)
          MEMLOG(size(dervs_xc(i,ii)%m)+size(dervs_xc_tmp(i,ii)%m))
          dervs_xc(i,ii)%m=0.0_r8_kind
         enddo

         endif init_m_dervs_xc
       enddo


       if (model_density) then
          allocate(rhs(n_ch,ispin),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          allocate(rhs_intermed(n_ch),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
       endif

     if(weight_grads) call weight_dervs_alloc()

       if(nl_calc_ph) then
          ! allocate space for 2nd derivatives of density
          allocate(nuc_sec_derrho(n_unique_atoms),stat=cpksalloc(138))  ! post_scf_main
          MEMLOG(size(nuc_sec_derrho))
          ASSERT(cpksalloc(138).eq.0)

#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             allocate(nuc_sec_derrho_core(n_unique_atoms),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
          end if
#endif

          do i=1,n_unique_atoms
             n_ea=unique_atoms(i)%n_equal_atoms
             allocate(nuc_sec_derrho(i)%m(vec_length,3,3,n_ea,ispin)&    !post_scf_main
                     ,stat=cpksalloc(137))
             MEMLOG(size(nuc_sec_derrho(i)%m))
             ASSERT(cpksalloc(137).eq.0)

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
               allocate(nuc_sec_derrho_core(i)%m(vec_length,3,3,n_ea,ispin)&
                   ,stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             end if
#endif
          end do
       end if

    else if(model_density .and. comp_exact) then
       allocate(rhs(n_ch,ispin),stat=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('post_scf_main: allocation of rhs failed')
       allocate(rhs_intermed(n_ch),stat=alloc_stat)
    endif grads

    charge_int= 0.0_r8_kind

    call xc_func_reset() ! zero the energy contributions

#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif

    FPP_TIMER_START(whole)
    if (model_density) then
       ! calculate final MDA energy and gradients (if do_grads = .TRUE.)
       call post_scf_calc_xcmda_etot_grads()

    else ! i.e. standard non mda  SCF procedure
    nograd: if(.not.do_grads) then
       !
       ! Single-Point
       !
       ! (just calculate final total xc-energy)
       !
       if(output_post_scf_main) &
            call write_to_output_units("POST_SCF_MAIN: post_scf_calc_xc_energy")
       call post_scf_calc_xc_energy()
    else nograd !i.e grada
       !
       ! Gradients
       !
#ifndef WITH_EXPERIMENTAL
#ifdef WITH_SECDER
       !
       ! Branch 1. Regular case.
       !
       ! Main  source for  calculating integrals  and  gradients.  The
       ! source that  handles both Grads  and Ders, LDA and  GGA. Note
       ! that this very procedure (with different set of arguments) is
       ! (ab)used for other purposes in CPKS iterations.
       !
       ! In  more  recent  Fortran  standard passing  an  unassociated
       ! pointer  or  or an  unallocated  allocatable  as an  optional
       ! argument is the  same as not passing it  at all.  Presence of
       ! the optional  arguments used to define the  operation mode of
       ! post_scf_calc_xc_en_and_gr_nl().    Beware   that  in   plain
       ! gradient runs, the variable "dervsw" remains unassociated, so
       ! that it  would not matter  that it is literally  "present" in
       ! this call:
       !
       ASSERT(allocated(grad_xc))
       call post_scf_calc_xc_en_and_gr_nl(POST_SCF_CALL, dervsw=dervsw) ! (1) - POST_SCF_MAIN
#else /* of ifdef WITH_SECDER */
       !
       ! Branch 2. WITH_SECDER=0.
       !
       ! Historic  version of  post_scf_calc_xc_en_and_gr_nl segfaults
       ! in LDA case:
       !
       if(.not.nl_calc_ph) then
          call post_scf_calc_xc_en_and_gr()
       else
          call post_scf_calc_xc_en_and_gr_nl()           ! (1a) - POST_SCF_MAIN
       endif
#endif /* of ifdef WITH_SECDER */
#else /* of ifndef WITH_EXPERIMENTAL */
#ifdef WITH_SECDER
       !
       ! Branch 3. WITH_EXPERIMENTAL=1, not regularly tested.
       !
       ! Some experimental case uses the  old LDA code also for second
       ! derivatives.
       !
       if(.not.nl_calc_ph) then
          !
          ! LDA Gradients
          !
          if(output_post_scf_main)   &
              call write_to_output_units("post_scf_calc_xc_en_and_gr")
          if( .not. integralpar_2dervs )then
            ! compute energy and gradients:
            call post_scf_calc_xc_en_and_gr(10) ! gradients only
          else
            ! compute energy, gradients, and some more ...:
            call post_scf_calc_xc_en_and_gr(20,dervsw=dervsw)
          endif
       else
          !
          ! GGA Gradients
          !
          ! (Here there is no alternative)
          if(output_post_scf_main)  call write_to_output_units &
          call post_scf_calc_xc_en_and_gr_nl(dervsw=dervsw)
       endif
#else /* of ifdef WITH_SECDER */
       !
       ! Branch 1. WITH_EXPERIMENTAL=1, WITH_SECDER=0.
       !
       ! If  there are  no  second derivatives,  the main  calculation
       ! routines are still used after the choice of the gradients.
       !
       if(.not.nl_calc_ph) then
          call post_scf_calc_xc_en_and_gr()
       else
          call post_scf_calc_xc_en_and_gr_nl()
       endif
#endif /* of ifdef WITH_SECDER */
#endif /* of ifndef WITH_EXPERIMENTAL */


    endif nograd
    endif
    FPP_TIMER_STOP(whole)

#ifdef FPP_TIMERS
    np = comm_size()
    rank = comm_rank()
    allocate(times(0:np-1, 6), stat = ierr)
    ASSERT(ierr == 0)
    times = 0
    times(rank, 1) = FPP_TIMER_VALUE(whole)
    times(rank, 2) = FPP_TIMER_VALUE(loop)
    times(rank, 3) = FPP_TIMER_VALUE(work)
    times(rank, 4) = FPP_TIMER_VALUE(sched)
    times(rank, 5) = FPP_TIMER_SLICE(sched)

    call comm_reduce(times)
#endif

    master: if (comm_rank() == 0) then

       ! build the final xc-gradient
       ! this step is necessary because we only consider gridpoints for
       ! one equal atom and because the gradient is not totalsymmetric
       if ( do_grads ) then
          call post_scf_average_gradient(grad_xc)

          if (weight_grads .and. split_gradients) then
             call post_scf_average_gradient(grad_grid)
          end if
       endif

#ifdef FPP_TIMERS
        print *, ""
        print *, "TIMINGS STATISTICS FOR LOAD BALANCING GRID PART"
        print *, "All times are in seconds"
        print *, "Elapsed times              :", maxval(times(:, 1)) * np ! total
        print *, "|__Elapsed times loop      :", maxval(times(:, 2)) * np ! loop1
        print *, " |__Work                   :", sum(times(:, 3)) ! work2
        print *, " |__Scheduling             :", sum(times(:, 4) - times(:, 5)) ! sched2 w/o last interval
        print *, " |__Imbalance at end       :", sum(times(:, 5)) - np * minval(times(:, 5)) !sched2 last interval- term
        print *, "  \__Varianz at end        :", (maxval(times(:, 5)) - minval(times(:, 5))) * np ! sched2 last ineterval max -min
        print *, " |__Termination            :", minval(times(:, 5)) * np ! sched2 last interval
        print *, "|Efficiency (Work/Elapsed) :", sum(times(:, 3)) / (maxval(times(:, 2)) * np) ! work2 / loop2
        print *, ""
#endif
       call post_scf_write_results()

       call stop_timer(timer_gridph)

    else master

        if (model_density) then
           ! TODO: Check if this substitutes the deallocation below
           ! However, deallocations should be done collectively by
           ! finalization subroutine.
           call fit_coeff_shutdown()
    !      deallocate(coeff_charge,stat=alloc_stat)
    !      if (alloc_stat /= 0) call error_handler &
    !           ("POST_SCF_MAIN: deallocation of coeff_charge failed")
    !      if (ispin > 1) then
    !         deallocate(coeff_spin,stat=alloc_stat)
    !         if (alloc_stat /= 0) call error_handler &
    !              ("POST_SCF_MAIN: deallocation of coeff_spin failed")
    !      endif
        endif
    endif master ! else

#ifdef FPP_TIMERS
    deallocate(times, stat=ierr)
    ASSERT(ierr == 0)
    call dlb_print_statistics(0)
#endif

    if (output_post_scf_main) call write_to_output_units ("grid_close")
    if(.not. operations_response.and..not.integralpar_cpksdervs) then
       ! release ALL grid information
       call grid_close(.false.)
    else
       ! we have to keep some grid information which is needed to rebuild the
       ! post scf grid in the response module
        call grid_close(.true.)
    end if

    if(output_post_scf_main) &
         call write_to_output_units("POST_SCF_MAIN: deallocations")

    deallocate(dfdrho,rho,fxc,charge_int,stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

#ifdef WITH_CORE_DENS
    if (pseudopot_present.and.core_density_setup) then
       deallocate(rho_core,grad_rho_core, stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if
#endif

    if (options_spin_orbit.and.spin_orbit_polarized) then
       deallocate(s_dens,s_abs,rho_ud,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('Deallocation in post_scf_main failed')
       if (nl_calc_ph) then
!!$          deallocate(gras_densx,gras_densy,gras_densz,stat=alloc_stat)
          deallocate(gras_dens,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation of gras_dens failed in su xc_setup')
!!$          deallocate(grarho_udx,grarho_udy,grarho_udz,stat=alloc_stat)
          deallocate(grarho_ud,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation of grarho_ud failed in post_scf_main')
       endif
    endif

    mdadel: if (model_density) then
       if (nl_calc_ph) then ! GGA calculation
          if (do_grads) then ! GGA gradients
             if (weight_grads) then ! d^2/dr^2 rho required
                call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch,&
                                       sec_ders=sec_ders_ch,&
                                       nuc_grads=nuc_grads_ch,&
                                       sec_nuc_ders=sec_nuc_ders_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    grads_core=grads_core,&
                                    sec_der_core=sec_der_core,&
                                    nuc_grads_core=nuc_grads_core,&
                                    nuc_sec_der_core=nuc_sec_der_core)
                end if
#endif
             else ! d^2/dr^2 rho not required
                call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch,&
                                       nuc_grads=nuc_grads_ch,&
                                       sec_nuc_ders=sec_nuc_ders_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    grads_core=grads_core,&
                                    nuc_grads_core=nuc_grads_core,&
                                    nuc_sec_der_core=nuc_sec_der_core)
                end if
#endif
             endif
          elseif(comp_exact) then
                call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) &
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    grads_core=grads_core)
#endif
          else ! GGA energy only
             call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    grads_core=grads_core)
                end if
#endif
          endif
       else ! LDA calculation
          if (do_grads) then ! LDA gradients
             if (weight_grads) then ! d/dr rho required
                call fit_fct_free("ch",fcts=fcts_ch,grads=grads_ch,&
                                       nuc_grads=nuc_grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    grads_core=grads_core,&
                                    nuc_grads_core=nuc_grads_core)
                end if
#endif
             else ! d/dr rho not required
                call fit_fct_free("ch",fcts=fcts_ch,nuc_grads=nuc_grads_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                    nuc_grads_core=nuc_grads_core)
                end if
#endif
             endif
          elseif(comp_exact) then
                call fit_fct_free("ch",fcts=fcts_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) &
                   call core_orbital_free(orbs_ob_core=orbs_ob_core)
#endif
          else ! LDA energy only
             call fit_fct_free("ch",fcts=fcts_ch)
#ifdef WITH_CORE_DENS
                if (pseudopot_present.and.core_density_setup) then
                   call core_orbital_free(orbs_ob_core=orbs_ob_core)
                end if
#endif
          endif
       endif
       ! shutdown orbitals_module
       call orbital_shutdown()
       call fitted_density_calc_close()
#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
         call fitted_core_dens_calc_close()
       end if
#endif

    else mdadel !i.e. standard nonmda procedure

       nograds_del:if(.not.do_grads) then
          if(.not.nl_calc_ph) then
             if (options_spin_orbit) then
                call orbital_free(orbs_spinor_ob=orbs_spinor_ob)
             else
                call orbital_free(orbs_ob)
             endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_free(orbs_ob_core=orbs_ob_core)
             end if
#endif

          else  ! i.e. nl_calc_ph

             if (options_spin_orbit) then
                call orbital_free(orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
             else
                call orbital_free(orbs_ob,orbs_grads)
             endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                  grads_core=grads_core)
             end if
#endif

          endif

       else nograds_del ! i.e. gradients

          if(.not.nl_calc_ph) then  ! LDA case

#ifdef WITH_SECDER
          if(integralpar_2dervs) then
             call secder_free(20)
          else
#endif
             call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads)
#ifdef WITH_SECDER
          endif
#endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                  grads_core=grads_core,&
                                  nuc_grads_core=nuc_grads_core)
             endif
#endif

          else  ! i,e. non LDA

             if(integralpar_2dervs) then
             call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads,&
                  nuc_sec_der=nuc_sec_der,nuc_3rd_der=nuc_3rd_der)
             else
             call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads,&
                  nuc_sec_der=nuc_sec_der)
             endif

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                call core_orbital_free(orbs_ob_core=orbs_ob_core,&
                                  grads_core=grads_core,&
                                  nuc_grads_core=nuc_grads_core,&
                                  nuc_sec_der_core=nuc_sec_der_core)
             endif
#endif
          endif
       endif nograds_del !/else

       ! shutdown orbitals_module
       call orbital_shutdown()
       call density_calc_close()

#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
          call fitted_core_dens_calc_close
       endif
#endif

    endif mdadel ! /else

    if ((comm_rank() == 0) .and. (output_timing_post_scf &
         .or. output_timing_detailedpostscf)) then
       if(output_post_scf_main) &
            call write_to_output_units("POST_SCF_MAIN: timer_print_postscf()")
       call timer_print_postscf()
    endif

    if (do_grads) then

       if (nl_calc_ph) then
          do ua=1,N_unique_atoms
             MEMLOG(-size(nuc_sec_derrho(ua)%m))
             deallocate(nuc_sec_derrho(ua)%m,stat=cpksalloc(137))
             ASSERT(alloc_stat.eq.0)
             cpksalloc(137)=1

#ifdef WITH_CORE_DENS
             if (pseudopot_present.and.core_density_setup) then
                deallocate(nuc_sec_derrho_core(ua)%m,stat=alloc_stat)
                ASSERT(alloc_stat.eq.0)
             end if
#endif

          enddo
          MEMLOG(-size(nuc_sec_derrho))
          deallocate(nuc_sec_derrho,stat=cpksalloc(138))
          ASSERT(cpksalloc(138).eq.0)
          cpksalloc(138)=1

#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             deallocate(nuc_sec_derrho_core,stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
          end if
#endif
       endif

       if (weight_grads) then
          do i_ua=1,N_unique_atoms
             MEMLOG(-size(graw(i_ua)%m))
             deallocate(graw(i_ua)%m, stat=cpksalloc(120)) ! post_scf main
             ASSERT(cpksalloc(120).eq.0)
             cpksalloc(120)=1
            if(.true..and.integralpar_2dervs) then
             do j_ua=1,N_unique_atoms
              MEMLOG(-size(dervsw(i_ua,j_ua)%m))
              deallocate(dervsw(i_ua,j_ua)%m,stat=cpksalloc(83))
              ASSERT(cpksalloc(83).eq.0)
              cpksalloc(83)=1
             enddo
            endif
          enddo

          deallocate(graw, stat=cpksalloc(119))
          ASSERT(cpksalloc(119).eq.0)
          cpksalloc(119)=1

           if(integralpar_2dervs) then
            deallocate(dervsw,stat=cpksalloc(82))
           ASSERT(cpksalloc(82).eq.0)
           cpksalloc(82)=1
           endif

       endif

       if (model_density) then
          deallocate(rhs,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        deallocate(rhs_intermed,stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)
       endif

       do ua=1,N_unique_atoms
          deallocate(nuc_grarho(ua)%m,stat=cpksalloc(133))
          ASSERT(cpksalloc(133).eq.0)
          cpksalloc(133)=1

          if(integralpar_2dervs) then

           do ii=1,N_unique_atoms
            MEMLOG(-size(nuc_dervsrho(ua,ii)%m))
            deallocate(nuc_dervsrho(ua,ii)%m,stat=cpksalloc(87))
            ASSERT(cpksalloc(87).eq.0)
            cpksalloc(87)=1

            if(nl_calc_ph) then
             MEMLOG(-size(nuc_dervs_grarho(ua,ii)%m))
             deallocate(nuc_dervs_grarho(ua,ii)%m,stat=cpksalloc(149))
             ASSERT(cpksalloc(149).eq.0)
             cpksalloc(149)=1
            endif

           enddo
          endif

          !-------------------------------------------------------------+
          ! Deallocation of nuc_grad_tau - elements                     |
          !-------------------------------------------------------------+
          if(is_on(xc_mgga)) then
            deallocate(nuc_grad_tau(ua)%m,stat=alloc_stat)
            ASSERT(alloc_stat.eq.0)
          end if
#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             deallocate(nuc_grarho_core(ua)%m,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          end if
#endif
       enddo

       !----------------------------------------------------------------+
       ! Deallocation of nuc_grad_tau - array                           |
       !----------------------------------------------------------------+
       deallocate(nuc_grarho,stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

       if(is_on(xc_mgga)) then
         deallocate(nuc_grad_tau,stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
       end if

       if(integralpar_2dervs) then
        MEMLOG(-size(nuc_dervsrho))
        deallocate( nuc_dervsrho, stat=cpksalloc(86))
        ASSERT(cpksalloc(86).eq.0)
        cpksalloc(86)=1
       if(nl_calc_ph) then
        MEMLOG(-size(nuc_dervs_grarho))
        deallocate(nuc_dervs_grarho,stat=cpksalloc(148))
        ASSERT(cpksalloc(148).eq.0)
        cpksalloc(148)=1
       endif
       endif

#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             deallocate(nuc_grarho_core,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          end if
#endif

    elseif (model_density.and.comp_exact) then
          deallocate(rhs,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
     deallocate(rhs_intermed,stat=alloc_stat)
     ASSERT(alloc_stat.eq.0)
    endif

#if 0 /* probable this is not used in post_scf call similar to dvdrho */
   if(integralpar_2dervs) then
    if(nl_calc_ph) then
     MEMLOG(-size(df_drhodgamma)-size(df_dgammadgamma))
     deallocate(df_drhodgamma,df_dgammadgamma,stat=cpksalloc(142))
     ASSERT(cpksalloc(142).eq.0)
     cpksalloc(142)=1
    endif
   endif
#endif

    if(output_post_scf_main) call write_to_output_units("POST_SCF_MAIN: done")

    contains

        subroutine weight_dervs_alloc()
        integer(kind=i4_kind)::i_ua,ni_ea,j_ua,nj_ea,n_ea,i

          allocate(graw(n_unique_atoms), stat=cpksalloc(119))
          ASSERT(cpksalloc(119).eq.0)

       if(integralpar_2dervs) then
          allocate(dervsw(n_unique_atoms,n_unique_atoms),stat=cpksalloc(82))
          ASSERT(cpksalloc(82).eq.0)
       endif

          do i_ua=1,n_unique_atoms
             ni_ea=unique_atoms(i_ua)%n_equal_atoms
             allocate(graw(i_ua)%m(vec_length,3,ni_ea), stat=cpksalloc(120)) !post_scf_main
             ASSERT(cpksalloc(120).eq.0)
             MEMLOG(size(graw(i_ua)%m))

           if(integralpar_2dervs) then
            do j_ua=1,n_unique_atoms
              nj_ea=unique_atoms(j_ua)%n_equal_atoms
              allocate( dervsw(i_ua,j_ua)%m(vec_length,3,ni_ea,3,nj_ea), &
                        stat=cpksalloc(83))
              ASSERT(cpksalloc(83).eq.0)
              MEMLOG(size(dervsw(i_ua,j_ua)%m))
            enddo
           endif

          enddo

          if (split_gradients) then
             allocate(grad_grid(N_moving_unique_atoms),stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             do i=1,N_moving_unique_atoms
                n_ea=unique_atoms(moving_unique_atom_index(i))%n_equal_atoms
                allocate(grad_grid(i)%m(3,n_ea),stat=alloc_stat)
                ASSERT(alloc_stat.eq.0)
             end do
          endif
       end subroutine weight_dervs_alloc
  end subroutine post_scf_main

#ifdef WITH_SECDER
  subroutine cpks_xc_main(job)
    ! input:  cpks%s1  cpks%B
    ! output: cpks%Qai cpks%AB
    !         cpks%h1
    !         dervs_xc
    use symmetry_data_module ! description of irreps
    use timer_module
    use time_module
    use unique_atom_methods, only: unique_atom_grad_information
    use iounitadmin_module
    use ph_cntrl, only:  nl_calc_ph
    use pointcharge_module
    use integralpar_module, only: integralpar_2dervs,integralpar_cpksdervs
    use cpksdervs_matrices,only:cpksalloc,cpks,cpks_Qxc_calculated
#ifndef WITH_EXPERIMENTAL
    use cpksdervs_matrices,only: cpks_grad_xc
#endif
    use xc_func, only: xc_func_reset
    integer(i4_kind), intent(in) :: job ! 1,2,3,4
    !** End of interface **************************************
    integer(kind=i4_kind) :: i,ua,n_ea
    integer(kind=i4_kind) :: ua2

!!!    avi=0.0_r8_kind
     DPRINT MyID,'cpks_xc_main(',job,') start'
     DPRINT MyID,'cpks_xc_main(',job,') nl_calc_ph           =',nl_calc_ph
     DPRINT MyID,'cpks_xc_main(',job,') integralpar_2dervs   =',integralpar_2dervs
     DPRINT MyID,'cpks_xc_main(',job,') integralpar_cpksdervs=',integralpar_cpksdervs
     DPRINT MyID,'cpks_xc_main(',job,') cpks_Qxc_calculated  =',cpks_Qxc_calculated
!    MEMSET(0)


    do_grads=operations_gradients
    split_gradients=options_split_gradients()
    model_density=options_xcmode()==xcmode_model_density .or. &
                  options_xcmode()==xcmode_extended_mda
    n_equal_max=maxval(unique_atoms(:)%n_equal_atoms) ! maximum number of equal_atoms

    ! now prepare integration grid
    call grid_main(post_scf=.true.)

    master_only: if (comm_rank() == 0) then

       call start_timer(timer_gridph)
       ! first send eigvecs and fitcoeffs to the slaves
          ! send eigenvectors to evaluate rho on the grid
          ! do_grads: eigenvalues will be required for orbital Pulay correction

       ! write to trace unit to show progress of calculation
       call dtrace("Calculate cpks_xc_data ")
    endif master_only


    ! set up the density calculations
    ispin=ssym%n_spin
    vec_length=machineparameters_veclen

    allocate(dvdrho(vec_length,2*ispin-1),stat=cpksalloc(9)) !cpks_xc_main
    ASSERT(cpksalloc(9).eq.0)
    MEMLOG(size(dvdrho))
    dvdrho=0.0_r8_kind

    if(nl_calc_ph) then      !cpks_xc_main
     allocate(df_drhodgamma(vec_length,5*ispin-4),       &
              df_dgammadgamma(vec_length,5*ispin-4),stat=cpksalloc(142))
     ASSERT(cpksalloc(142).eq.0)
     MEMLOG(size(df_drhodgamma)+size(df_dgammadgamma))
     df_drhodgamma=0.0_r8_kind
     df_dgammadgamma=0.0_r8_kind
    endif

!       print*,'call density_calc_setup'
!       MEMSET(0)
       call density_calc_setup()
!       print*,'density_calc_setup done'
!       MEMSET(0)

       ! allocation of orbitals and their derivatives for different cases

       if(.not.nl_calc_ph) then !lda
          call secder_allocate(job)
       else ! i.e. gga
             call orbital_setup(vec_length,do_gradients=.true.,do_sec_der=.true.)
             call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads,& ! (6)
                                   nuc_sec_der=nuc_sec_der, &
                                   phi_ob=phi_ob) ! NONLDA grads cpks_xc_main
       endif
       charge_types = 1

    allocate(dfdrho(vec_length,ispin), rho(vec_length,ispin),&
             fxc(vec_length),charge_int(n_unique_atoms,charge_types),& !charge_int to be deleted
             stat=cpksalloc(10))
    ASSERT(cpksalloc(10).eq.0)
    MEMLOG(size(dfdrho)+size(rho)+size(fxc)+size(charge_int))

    if(do_grads) then
       allocate(nuc_grarho(n_unique_atoms),stat=cpksalloc(139))
       ASSERT(cpksalloc(139).eq.0)
       MEMLOG(size(nuc_grarho))
    if(integralpar_2dervs) then

       allocate(nuc_grarho_imp(n_unique_atoms),stat=cpksalloc(93))
       ASSERT(cpksalloc(93).eq.0)
       MEMLOG(size(nuc_grarho_imp))
       if(nl_calc_ph) then
        allocate(nuc_gragamma_imp(n_unique_atoms),stat=cpksalloc(146))
        ASSERT(cpksalloc(146).eq.0)
        MEMLOG(size(nuc_gragamma_imp))
       endif

     if(integralpar_cpksdervs) then
       allocate( dervsrho_imp(n_unique_atoms,n_unique_atoms),stat=cpksalloc(100))
       ASSERT(cpksalloc(100).eq.0)
       MEMLOG(size(dervsrho_imp))

       if(nl_calc_ph) then
       allocate( dervs_grarho_imp(n_unique_atoms,n_unique_atoms), stat=cpksalloc(151))
       ASSERT(cpksalloc(151).eq.0)
       MEMLOG(size(dervs_grarho_imp))
       endif

     endif
    endif

       do i=1,n_unique_atoms
          n_ea=unique_atoms(i)%n_equal_atoms
          allocate(nuc_grarho(i)%m(vec_length,3,n_ea,ispin),stat=cpksalloc(140))
          ASSERT(cpksalloc(140).eq.0)
          MEMLOG(size(nuc_grarho(i)%m))

        if(integralpar_2dervs) then

          allocate(nuc_grarho_imp(i)%m(vec_length,3,n_ea,ispin),stat=cpksalloc(92))
          ASSERT(cpksalloc(92).eq.0)
          MEMLOG(size(nuc_grarho_imp(i)%m))
          nuc_grarho_imp(i)%m=0.0_r8_kind

          if(nl_calc_ph) then
          allocate(nuc_gragamma_imp(i)%m(vec_length,3,3,n_ea,ispin),stat=cpksalloc(147))
          ASSERT(cpksalloc(147).eq.0)
          MEMLOG(size(nuc_gragamma_imp(i)%m))
          nuc_gragamma_imp(i)%m=0.0_r8_kind
          endif

         if(allocated(cpks)) then

          do ua=1,n_unique_atoms
          allocate(dervsrho_imp(i,ua)%m(3,n_ea,3,&
                                        unique_atoms(ua)%n_equal_atoms,ispin), &
                   stat=cpksalloc(101))
          ASSERT(cpksalloc(101).eq.0)
          MEMLOG(size(dervsrho_imp(i,ua)%m))
          dervsrho_imp(i,ua)%m=0.0_r8_kind

           if(nl_calc_ph) then
            allocate( dervs_grarho_imp(i,ua)%m(3,n_ea,3,&
                       unique_atoms(ua)%n_equal_atoms,ispin), &
            stat=cpksalloc(152))
            ASSERT(cpksalloc(152).eq.0)
            MEMLOG(size(dervs_grarho_imp(i,ua)%m))
            dervs_grarho_imp(i,ua)%m=0.0_r8_kind
           endif
          enddo
         endif

        endif

       end do

#ifndef WITH_EXPERIMENTAL
     if(.not.cpks_Qxc_calculated) then
       allocate(cpks_grad_xc(N_moving_unique_atoms),stat=cpksalloc(159))
       ASSERT(cpksalloc(159).eq.0)
       MEMLOG(size(cpks_grad_xc))

       do i=1,N_moving_unique_atoms
          n_ea=unique_atoms(moving_unique_atom_index(i))%n_equal_atoms
          allocate(cpks_grad_xc(i)%m(3,n_ea),stat=cpksalloc(14))
          ASSERT(cpksalloc(14).eq.0)
          MEMLOG(size(cpks_grad_xc(i)%m))
       end do
     endif
#endif

       if(weight_grads) then
          allocate(graw(n_unique_atoms),stat=cpksalloc(119))
          ASSERT(cpksalloc(119).eq.0)
          MEMLOG(size(graw))
          do i=1,n_unique_atoms
             n_ea=unique_atoms(i)%n_equal_atoms
             allocate(graw(i)%m(vec_length,3,n_ea),stat=cpksalloc(120))
          ASSERT(cpksalloc(120).eq.0)
          MEMLOG(size(graw(i)%m))
          enddo
       end if

       if(nl_calc_ph) then
          ! allocate space for 2nd derivatives of density
          allocate(nuc_sec_derrho(n_unique_atoms),stat=cpksalloc(138)) ! cpks_xc_main
          ASSERT(cpksalloc(138).eq.0)
          MEMLOG(size(nuc_sec_derrho))
          do i=1,n_unique_atoms
             n_ea=unique_atoms(i)%n_equal_atoms
             allocate(nuc_sec_derrho(i)%m(vec_length,3,3,n_ea,ispin)&  ! cpks_xc_main
                     ,stat=cpksalloc(137))
          ASSERT(cpksalloc(137).eq.0)
          MEMLOG(size(nuc_sec_derrho(i)%m))
          end do
       endif

    end if

    charge_int= 0.0_r8_kind

    call xc_func_reset()

       if(.not.nl_calc_ph) then !LDA
#ifdef WITH_EXPERIMENTAL
          WARN('calling LDA version')
          DPRINT MyID,'cpks_xc_main(',job,') call post_scf_calc_xc_en_and_gr(LDA)'
          call post_scf_calc_xc_en_and_gr(job,nuc_grarho_imp=nuc_grarho_imp) !(2) in cpks_main
#else
          call post_scf_calc_xc_en_and_gr_nl(job, nuc_grarho_imp=nuc_grarho_imp) ! in cpks_main
#endif
       else !i.e. nonlda
          DPRINT MyID,'cpks_xc_main(',job,') call post_scf_calc_xc_en_and_gr(GGA)'
          call post_scf_calc_xc_en_and_gr_nl(job, nuc_grarho_imp=nuc_grarho_imp,   &
                                             nuc_gragamma_imp=nuc_gragamma_imp) ! in cpks_main
       endif
       DPRINT MyID,'cpks_xc_main(',job,') done post_scf_calc_xc_en_and_gr(...)'
!       print*,'en_and_gr done'
!       MEMSET(0)

#ifndef WITH_EXPERIMENTAL
     if(.not.cpks_Qxc_calculated) then
       do i=1,N_moving_unique_atoms
          n_ea=unique_atoms(moving_unique_atom_index(i))%n_equal_atoms
          MEMLOG(-size(cpks_grad_xc(i)%m))
          deallocate(cpks_grad_xc(i)%m,stat=cpksalloc(14))
          ASSERT(cpksalloc(14).eq.0)
          cpksalloc(14)=1
       end do
       MEMLOG(-size(cpks_grad_xc))
       deallocate(cpks_grad_xc,stat=cpksalloc(159))
       ASSERT(cpksalloc(159).eq.0)
       cpksalloc(159)=1
     endif
#endif

    if (comm_rank() == 0) then
    ! no need to calc grad_xc with cpks_main
       ! build the final xc-gradient
       ! this step is necessary because we only consider gridpoints for
       ! one equal atom and because the gradient is not totalsymmetric

       call stop_timer(timer_gridph)

          ! we have to keep some grid information which is needed to rebuild the
          ! post scf grid in the response module
          ! (here grid_close only on the master wanted?)
          !call grid_close(.false.)
!          print*,'grid_close done'
!          MEMSET(0)
    endif
    call grid_close(.true.)

    MEMLOG(-size(dvdrho))
    deallocate(dvdrho,stat=cpksalloc(9))
    ASSERT(cpksalloc(9).eq.0)
    cpksalloc(9)=1
    if(nl_calc_ph) then
     MEMLOG(-size(df_drhodgamma)-size(df_dgammadgamma))
     deallocate(df_drhodgamma,df_dgammadgamma,stat=cpksalloc(142))
     ASSERT(cpksalloc(142).eq.0)
     cpksalloc(142)=1
    endif

    MEMLOG(-size(dfdrho)-size(rho)-size(fxc)-size(charge_int))
    deallocate(dfdrho,rho,fxc,charge_int,stat=cpksalloc(10))
    ASSERT(cpksalloc(10).eq.0)
    cpksalloc(10)=1



          if(.not.nl_calc_ph) then  ! LDA
            call secder_free(job)
          else  ! i.e. nonLDA
!             print*,'orbital_free stat'
!             MEMSET(0)
             call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads,&
                               nuc_sec_der=nuc_sec_der, &
                               phi_ob=phi_ob)  ! cpks_xc_main nonLDA
!             print*,'done'
!             MEMSET(0)
          endif

       call orbital_shutdown()
!       print*,'orbital_shutdown done'
!       MEMSET(0)
       call density_calc_close()
!       print*,'density_calc_close'
!       MEMSET(0)

    if ((comm_rank() == 0) .and. (output_timing_post_scf &
         .or. output_timing_detailedpostscf)) then
       call timer_print_postscf()
    endif

       if (nl_calc_ph) then
          do ua=1,N_unique_atoms
             MEMLOG(-size(nuc_sec_derrho(ua)%m))
             deallocate(nuc_sec_derrho(ua)%m,stat=cpksalloc(137))  ! cpks_main
             ASSERT(cpksalloc(137).eq.0)
             cpksalloc(137)=1
          enddo
          MEMLOG(-size(nuc_sec_derrho))
          deallocate(nuc_sec_derrho,stat=cpksalloc(138))           ! cpks_main
          ASSERT(cpksalloc(138).eq.0)
          cpksalloc(138)=1
       endif

       if (weight_grads) then
          do ua=1,N_unique_atoms
             MEMLOG(-size(graw(ua)%m))
             deallocate(graw(ua)%m,stat=cpksalloc(120)) ! cpks_xc_main
             ASSERT(cpksalloc(120).eq.0)
             cpksalloc(120)=1
          enddo

          MEMLOG(-size(graw))
          deallocate(graw,stat=cpksalloc(119)) ! cpks_xc_main
          ASSERT(cpksalloc(119).eq.0)
          cpksalloc(119)=1
       endif

       do ua=1,N_unique_atoms
          MEMLOG(-size(nuc_grarho(ua)%m))
          deallocate(nuc_grarho(ua)%m,stat=cpksalloc(140))
          ASSERT(cpksalloc(140).eq.0)
          cpksalloc(140)=1
       enddo
       MEMLOG(-size(nuc_grarho))
       deallocate(nuc_grarho,stat=cpksalloc(139))
       ASSERT(cpksalloc(139).eq.0)
       cpksalloc(139)=1

       if(integralpar_2dervs) then

   do ua=1,N_unique_atoms

          MEMLOG(-size(nuc_grarho_imp(ua)%m))
          deallocate(nuc_grarho_imp(ua)%m,stat=cpksalloc(92))
          ASSERT(cpksalloc(92).eq.0)
          cpksalloc(92)=1

         if(nl_calc_ph) then
          MEMLOG(-size(nuc_gragamma_imp(ua)%m))
          deallocate(nuc_gragamma_imp(ua)%m,stat=cpksalloc(147))
          ASSERT(cpksalloc(147).eq.0)
          cpksalloc(147)=1
         endif

       if(allocated(cpks)) then
        do ua2=1,N_unique_atoms
          MEMLOG(-size(dervsrho_imp(ua,ua2)%m))
          deallocate(dervsrho_imp(ua,ua2)%m,stat=cpksalloc(101))
          ASSERT(cpksalloc(101).eq.0)
          cpksalloc(101)=1
          if(nl_calc_ph) then
           MEMLOG(-size(dervs_grarho_imp(ua,ua2)%m))
           deallocate(dervs_grarho_imp(ua,ua2)%m,stat=cpksalloc(152))
           ASSERT(cpksalloc(152).eq.0)
           cpksalloc(152)=1
          endif
        enddo
       endif
   enddo

       MEMLOG(-size(nuc_grarho_imp))
       deallocate(nuc_grarho_imp,stat=cpksalloc(93))
       ASSERT(cpksalloc(93).eq.0)
       cpksalloc(93)=1
       if(nl_calc_ph) then
        MEMLOG(-size(nuc_gragamma_imp))
        deallocate(nuc_gragamma_imp,stat=cpksalloc(146))
        ASSERT(cpksalloc(146).eq.0)
        cpksalloc(146)=1
       endif
       if(integralpar_cpksdervs) then
        MEMLOG(-size(dervsrho_imp))
        deallocate(dervsrho_imp,stat=cpksalloc(100))
        ASSERT(cpksalloc(100).eq.0)
        cpksalloc(100)=1
       if(nl_calc_ph) then
        MEMLOG(-size(dervs_grarho_imp))
        deallocate(dervs_grarho_imp,stat=cpksalloc(151))
        ASSERT(cpksalloc(151).eq.0)
        cpksalloc(151)=1
       endif
       endif
       endif

       DPRINT MyID,'cpks_xc_main(',job,') END!!!!!!!!!!!!!!!!!!!'
!    MEMSET(0)

  end subroutine cpks_xc_main

  subroutine secder_allocate(job)
    use gradient_data_module, only: gradient_data_n_gradients
    implicit none
    integer(i4_kind), intent(in) :: job
    ! *** end of interface ***

    integer(i4_kind) :: istat

    DPRINT job,'N_GRADS=',gradient_data_n_gradients

    select case ( job )
    case ( 20 )
      call orbital_setup(vec_length,do_gradients=.true.,do_sec_der=.true.)
      call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                           , nuc_sec_der=nuc_sec_der               &
                           , phi_ob=phi_ob                         &
                           )
    case ( 21, 24 )
      call orbital_setup(vec_length,do_gradients=.true.)
      call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                           , nuc_sec_der=nuc_sec_der               &
                           , phi_ob=phi_ob                         &
                           , phi_gr_nuc=phi_gr_nuc                 &
                           , num_gr_nuc=gradient_data_n_gradients  &
                           )
    case ( 22, 23 )
      call orbital_setup(vec_length,do_gradients=.true.)
      call orbital_allocate(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                           , nuc_sec_der=nuc_sec_der               & ! ???
                           , phi_ob=phi_ob                         &
                           )
    case default
      print *,'ERROR: no such job',job,'in secder_allocate'
      ABORT('no such job')
    end select

    allocate( rho_gr_nuc(vec_length,gradient_data_n_gradients,ispin) &
            , rho_gr_imp(vec_length,gradient_data_n_gradients,ispin) &
            , wts_gr_nuc(vec_length,gradient_data_n_gradients)       &
            , eny_sd_imp(gradient_data_n_gradients,gradient_data_n_gradients,ispin) &
            , stat=istat )
    ASSERT(istat==0)
  end subroutine secder_allocate

  subroutine secder_free(job)
    implicit none
    integer(i4_kind), intent(in) :: job
    ! *** end of interface ***

    integer(i4_kind) :: istat

    DPRINT job,'call orbital_free()'

    select case ( job )
    case ( 20 )
      call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                       , nuc_sec_der=nuc_sec_der               &
                       , phi_ob=phi_ob                         &
                       )
    case ( 21, 24 )
      call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                       , nuc_sec_der=nuc_sec_der               &
                       , phi_ob=phi_ob                         &
                       , phi_gr_nuc=phi_gr_nuc                 &
                       )
    case ( 22, 23 )
      call orbital_free(orbs_ob,orbs_grads,nuc_grads=nuc_grads &
                       , nuc_sec_der=nuc_sec_der               &
                       , phi_ob=phi_ob                         &
                       )
    case default
      print *,'ERROR: no such job',job,'in secder_free'
      ABORT('no such job')
    end select

    deallocate( rho_gr_nuc &
              , rho_gr_imp &
              , wts_gr_nuc &
              , eny_sd_imp &
              , stat=istat )
    ASSERT(istat==0)
  end subroutine secder_free
#endif

  subroutine post_scf_deallocate_grad_xc()
    use cpksdervs_matrices,only: cpksalloc
    ! called by main_gradient
    !** End of interface *****************************************
    integer(kind=i4_kind) :: i_unique,alloc_stat,i_unique2

    if (allocated(grad_grid)) then
       do i_unique=1,N_moving_unique_atoms
          deallocate(grad_grid(i_unique)%m,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
       end do
       deallocate(grad_grid,stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if

    if (allocated(grad_xc)) then
        do i_unique=1,N_moving_unique_atoms
            MEMLOG(-size(grad_xc(i_unique)%m))
            deallocate(grad_xc(i_unique)%m,stat=cpksalloc(22))
            ASSERT(cpksalloc(22).eq.0)
            cpksalloc(22)=1
        enddo
        deallocate(grad_xc,stat=cpksalloc(21))
        ASSERT(cpksalloc(21).eq.0)
        cpksalloc(21)=1
    endif

       if (allocated(dervs_xc)) then
        do i_unique=1,N_moving_unique_atoms
           do i_unique2=1,N_moving_unique_atoms
               MEMLOG(-size(dervs_xc(i_unique,i_unique2)%m)-size(dervs_xc_tmp(i_unique,i_unique2)%m))
               deallocate( dervs_xc(i_unique,i_unique2)%m, &
                    dervs_xc_tmp(i_unique,i_unique2)%m, &
                    stat=cpksalloc(81))
           ASSERT(cpksalloc(81).eq.0)
           cpksalloc(81)=1
          enddo
        enddo
        MEMLOG(-size(dervs_xc)-size(dervs_xc_tmp))
        deallocate(dervs_xc,dervs_xc_tmp,stat=cpksalloc(80))
        ASSERT(cpksalloc(80).eq.0)
        cpksalloc(80)=1
    endif

!    MEMLOG(-size(grad_xc))

  end subroutine post_scf_deallocate_grad_xc

  subroutine post_scf_average_gradient(grad_final)
    ! purpose : builds final gradient by averaging over all
    !           equal atoms. this is necesarry because the
    !           integration grid and the gradients are not totalsymmetric
    !
    type(arrmat2),intent(inout) :: grad_final(:)
    !** End of interface *****************************************
    real(kind=r8_kind) :: rotmat_tot(3,3),help_arr_local(3,n_equal_max)
    real(kind=r8_kind),pointer :: rotmat1(:,:),rotmat2(:,:)
    integer :: i_unique,i_equal1,i_equal2,kk,n_ea

    do i_unique=1,N_moving_unique_atoms
       n_ea = unique_atoms(moving_unique_atom_index(i_unique))%n_equal_atoms
       help_arr_local(:,1:n_ea)=0.0_r8_kind
       do i_equal1=1,n_ea
       ! calc totsym contrib for i_equal1
          do i_equal2=1,n_ea
             rotmat1=>unique_atom_grad_info(i_unique)%m(:,:,i_equal1)
             rotmat2=>unique_atom_grad_info(i_unique)%m(:,:,i_equal2)
             rotmat_tot=matmul(transpose(rotmat1),rotmat2)
             do kk=1,3
                help_arr_local(kk,i_equal1)=&
                     help_arr_local(kk,i_equal1)+&
                     sum(grad_final(i_unique)%m(:,i_equal2)*rotmat_tot(kk,:))
             end do
          end do
       end do
       grad_final(i_unique)%m=help_arr_local(:,1:n_ea)
    end do
  end subroutine post_scf_average_gradient

  subroutine post_scf_average_dervs(dervs_final)
    ! purpose : builds final gradient by averaging over all
    !           equal atoms. this is necesarry because the
    !           integration grid and the gradients are not totalsymmetric
    !
    type(arrmat4),intent(inout) :: dervs_final(:,:)
    !** End of interface *****************************************

    real(kind=r8_kind) :: rotmat_tot(3,3),help_arr_local(3,n_equal_max,3,n_equal_max)
    real(kind=r8_kind),pointer :: rotmat1(:,:),rotmat2(:,:)
    integer :: i_unique,i_unique2,i_equal1,i_equal2,kk,n_ea,n_ea2,k

    uniques1: do i_unique=1,N_moving_unique_atoms
       n_ea = unique_atoms(moving_unique_atom_index(i_unique))%n_equal_atoms
     uniques2: do i_unique2=1,N_moving_unique_atoms
        n_ea2= unique_atoms(moving_unique_atom_index(i_unique2))%n_equal_atoms
        help_arr_local(:,1:n_ea,:,:n_ea2)=0.0_r8_kind
       do i_equal1=1,n_ea
          do i_equal2=1,n_ea
             rotmat1=>unique_atom_grad_info(i_unique)%m(:,:,i_equal1)
             rotmat2=>unique_atom_grad_info(i_unique)%m(:,:,i_equal2)
             rotmat_tot=matmul(transpose(rotmat1),rotmat2)
             do kk=1,3
              do k=1,3
                help_arr_local(kk,i_equal1,:,:n_ea2)=&
                   help_arr_local(kk,i_equal1,:,:n_ea2)+&
                      dervs_final(i_unique,i_unique2)%m(k,i_equal2,:,:n_ea2)*rotmat_tot(kk,k)
              enddo
             enddo
          end do
       end do
!       dervs_final(i_unique,i_unique2)%m=help_arr_local
       dervs_final(i_unique,i_unique2)%m=0.0_r8_kind
      do i_equal1=1,n_ea2
       do i_equal2=1,n_ea2
             rotmat1=>unique_atom_grad_info(i_unique2)%m(:,:,i_equal1)
             rotmat2=>unique_atom_grad_info(i_unique2)%m(:,:,i_equal2)
             rotmat_tot=matmul(transpose(rotmat1),rotmat2)
             do kk=1,3
              do k=1,3
             dervs_final(i_unique,i_unique2)%m(:,1:n_ea,kk,i_equal1)= &
               dervs_final(i_unique,i_unique2)%m(:,1:n_ea,kk,i_equal1) &
                  +help_arr_local(:,1:n_ea,kk,i_equal2)*rotmat_tot(kk,k)
              enddo
             enddo
       enddo
      enddo

     enddo uniques2
    enddo uniques1
  end subroutine post_scf_average_dervs


  subroutine ph_sndrcv(gradients_also, send_grad_grid, send_mda_coeff)
    ! Purpose : send results from slaves to the master
    !            - charge_int
    !            - xc_func_reduce
    !           if gradients_also:
    !            - grad_xc
    !            - grad_grid ( if send_grad_grid)
    !           if send_mda_coeff:
    !            - rhs
    !           if comm_dervs_xc:
    !           - dervs_xc
    !------------ Modules used ------------------- ---------------
    use xc_func, only: xc_func_reduce
    use cpksdervs_matrices,only: comm_dervs_xc
    use comm
    implicit none
    logical, intent(in)  :: gradients_also, send_grad_grid, send_mda_coeff
    !** End of interface *****************************************
    integer(kind=i4_kind) :: i_ua, ii_ua

    call comm_reduce(charge_int)
    call xc_func_reduce()

    if (gradients_also) then
      do i_ua=1,N_moving_unique_atoms
            call comm_reduce(grad_xc(i_ua)%m)
      enddo

      if (send_grad_grid) then
        do i_ua=1,N_moving_unique_atoms
          call comm_reduce(grad_grid(i_ua)%m)
        enddo
      endif

      if(comm_dervs_xc) then
        do i_ua=1,N_moving_unique_atoms
          do ii_ua=1,N_moving_unique_atoms
            call comm_reduce( dervs_xc(i_ua,ii_ua)%m)
          enddo
        enddo
      endif
    endif

    if (send_mda_coeff) then
      call comm_reduce(rhs)
    endif
  end subroutine ph_sndrcv

  subroutine post_scf_write_results
    use iounitadmin_module
    use ph_cntrl
    use xc_func, only: xc_func_write
    implicit none
    ! purpose : write the results of the post_scf part
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i,i_ua,i_ea,i_ca
    logical               :: sum_up_contributions
    real(r8_kind)         :: g(3) ! net xc-force


    write(output_unit,*) ''
    write(output_unit,*) ''
    write(output_unit,*) '************* POST SCF PART ***************'
    write(output_unit,*)
    write(output_unit,*)'**** charges ****'
    1000 format(A,I5,' : ',F25.12)
    1001 format(A,5X,' : ',F25.12)
    if (model_density) then
       if (ispin == 1) then
          write(output_unit,'( A)') &
               "The numerically integrated charges:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'charge for atom', i, charge_int(i,1)
          enddo
          write(output_unit,1001)'total charge   ',sum(charge_int(:,1))
          write(output_unit,'(/A)') &
               "The numerically integrated truncated charges:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'charge for atom', i, charge_int(i,2)
          enddo
          write(output_unit,1001)'total charge   ',sum(charge_int(:,2))
          write(output_unit,'()')
       else
          write(output_unit,'( A)') &
               "The numerically integrated charges:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'charge for atom', i, &
                                     charge_int(i,1) + charge_int(i,2)
          enddo
          write(output_unit,1001)'total charge   ', &
                                  sum(charge_int(:,1)) + sum(charge_int(:,2))
          write(output_unit,'(/A)') &
               "The numerically integrated truncated charges:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'charge for atom', i, &
                                     charge_int(i,3) + charge_int(i,4)
          enddo
          write(output_unit,1001)'total charge   ', &
                                  sum(charge_int(:,3)) + sum(charge_int(:,4))
          write(output_unit,'(/A)') &
               "The numerically integrated magn. moments:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'magn. moment for atom', i, &
                                     charge_int(i,1) - charge_int(i,2)
          enddo
          write(output_unit,1001)'total magn. moment   ', &
                                  sum(charge_int(:,1)) - sum(charge_int(:,2))
          write(output_unit,'(/A)') &
               "The numerically integrated truncated magn. moments:"
          do i=1,n_unique_atoms
             write(output_unit,1000)'magn. moment for atom', i, &
                                     charge_int(i,3) - charge_int(i,4)
          enddo
          write(output_unit,1001)'total magn. moment   ', &
                                  sum(charge_int(:,3)) - sum(charge_int(:,4))
          write(output_unit,'()')
       endif
    else
       do i=1,n_unique_atoms
          write(output_unit,'(A42,I5,A4,F25.12)') &
               'The numerically integrated charge for atom', i,'is:', &
               charge_int(i,1)
       enddo
       write(output_unit,*)
       write(output_unit,'(A50,F25.12)') &
            'The numerically integrated charge for all atoms is' , &
            sum(charge_int(:,1))
    endif

    call xc_func_write(output_unit)

    1004 FORMAT('PHXC: ',A15)
    if(do_grads) then
       sum_up_contributions = weight_grads .and. options_split_gradients()
       1003 FORMAT('PHXC: ',I15,3F18.12)
       1005 FORMAT('PHXC: ',A15,I5)
       write(output_unit,*)
       write(output_unit,1004) 'xc_gradients:'
       ! compute and print the net xc-force on the whole system as well
       g(1:3) = 0.0_r8_kind
       do i_ua=1,N_moving_unique_atoms
          i_ca = moving_unique_atom_index(i_ua)
          write(output_unit,1005) 'Unique_atom:',i_ca
          do i_ea=1,unique_atoms(i_ca)%n_equal_atoms
             if (sum_up_contributions) then
                write(output_unit,1003) i_ea, &
                     grad_xc(i_ua)%m(:,i_ea) + grad_grid(i_ua)%m(:,i_ea)
                g(1:3) = g(1:3) + grad_xc(i_ua)%m(:,i_ea) + grad_grid(i_ua)%m(:,i_ea)
             else
                write(output_unit,1003) i_ea, &
                     grad_xc(i_ua)%m(:,i_ea)
                g(1:3) = g(1:3) + grad_xc(i_ua)%m(:,i_ea)
             end if
          end do
       enddo
       if (sum_up_contributions) then
          write(output_unit,1005) 'SUM  (XC+G):',-1
       else
          write(output_unit,1005) 'SUM    (XC):',-1
       endif
          write(output_unit,1003)               -1, g(1:3)
    endif

  end subroutine post_scf_write_results

  !***************************************************************

  subroutine post_scf_calc_xc_energy()
    ! purpose : routine does the main calculations
    !           - calculation of orbitals
    !           - calculation of becke weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !** End of interface *****************************************
    use time_module
    use timer_module
    use ph_cntrl
    use iounitadmin_module
    use xc_func
    implicit none

    integer(i4_kind) :: vla, xyz

    real(kind=r8_kind),allocatable :: grarho(:,:,:) !(vl,X:Z,ISPIN)
    real(kind=r8_kind),allocatable :: gamma(:,:)
    real(kind=r8_kind),allocatable :: tau(:,:) !(vl,ISPIN)
    logical :: send_grad_grid, send_mda_coeff
    integer                        :: nspin
    integer :: i,alloc_stat
    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths

    ! maybe one should make *spin an input argument:
    nspin = ispin

    ! spin-orbit is usualy done in unpolarized fashion, unless ...
    if ( options_spin_orbit .and. spin_orbit_polarized ) then
      ASSERT(ispin==1)
      nspin = 2
    endif

    ! allocate gradients of density
    allocate(grarho(vec_length, X:Z, nspin), stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! for polarized calculations the last dimension is 3, otherwise 1:
    allocate(gamma(vec_length, 1 + 2*(nspin - 1)), stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! FIXME: not yet used in SO calculations:
    allocate(tau(vec_length, nspin), stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! prepare some timers!
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)
    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
    if(.not.nl_calc_ph) then        ! LDA calculation
       sp_orb: if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          if (spin_orbit_polarized) then
             !
             ! OPEN SHELL
             !
             ASSERT(nspin==2)
             FPP_TIMER_START(loop)
             FPP_TIMER_START(sched)
             do while (more_grid_atom(vec_length, i, grdpts, grdwts))    ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
                FPP_TIMER_STOP(sched)
                FPP_TIMER_START(work)
                rho=0.0_r8_kind
                dfdrho=0.0_r8_kind
                fxc=0.0_r8_kind
                vla=size(grdpts,1)
                call start_timer(timer_gridph_grdwts)
                grdwts(1:vla)=grdwts(1:vla)*atomicweight(i,grdpts(1:vla,:))*&
                              unique_atoms(i)%n_equal_atoms
                call stop_timer(timer_gridph_grdwts)
                call start_timer(timer_gridph_orbitals)
                call orbital_calculate(grdpts(1:vla,1:3),&
                     vla,orbs_spinor_ob=orbs_spinor_ob)
                call stop_timer(timer_gridph_orbitals)
                call start_timer(timer_gridph_density)
                call density_calc(vla,rho,orbs_spinor_ob=orbs_spinor_ob,s_dens=s_dens)
                call stop_timer(timer_gridph_density)
                call start_timer(timer_gridph_functionals)
                s_abs(1:vla) = sqrt(s_dens(1:vla,1)**2+s_dens(1:vla,2)**2+&
                     s_dens(1:vla,3)**2)
                rho_ud(1:vla,1) =  0.5_r8_kind*(rho(1:vla,1) + s_abs(1:vla)) ! spin up
                rho_ud(1:vla,2) =  0.5_r8_kind*(rho(1:vla,1) - s_abs(1:vla)) ! spin down

                call xc_functionals(vla, nspin, rho_ud, GRDWTS=grdwts)

                call stop_timer(timer_gridph_functionals)
                ! performing integration over the grid
                charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                     (grdwts(1:vla), 2, nspin))
                FPP_TIMER_STOP(work)
                FPP_TIMER_START(sched)
             end do! loop over gridpoints
          else
             !
             ! CLOSED SHELL
             !
             ASSERT(nspin==1)
             FPP_TIMER_START(loop)
             FPP_TIMER_START(sched)
             do while (more_grid_atom(vec_length, i, grdpts, grdwts))   ! loop over gridpoints
                FPP_TIMER_STOP(sched)
                FPP_TIMER_START(work)
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
                rho=0.0_r8_kind
                dfdrho=0.0_r8_kind
                fxc=0.0_r8_kind
                vla=size(grdpts,1)
                call start_timer(timer_gridph_grdwts)
                grdwts(1:vla)=grdwts(1:vla)*&
                     atomicweight(i,grdpts(1:vla,:))*&
                     unique_atoms(i)%n_equal_atoms
                call stop_timer(timer_gridph_grdwts)
                call start_timer(timer_gridph_orbitals)
                call orbital_calculate(grdpts(1:vla,1:3),&
                     vla,orbs_spinor_ob=orbs_spinor_ob)
                call stop_timer(timer_gridph_orbitals)
                call start_timer(timer_gridph_density)
                call density_calc(vla,rho,orbs_spinor_ob=orbs_spinor_ob)
                call stop_timer(timer_gridph_density)
                call start_timer(timer_gridph_functionals)

                call xc_functionals(vla, nspin, rho, GRDWTS=grdwts)

                call stop_timer(timer_gridph_functionals)
                ! performing integration over the grid
                charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                     (grdwts(1:vla), 2, nspin))
                FPP_TIMER_STOP(work)
                FPP_TIMER_START(sched)
             end do! loop over gridpoints
          endif! closed /open shell

       else sp_orb
          !
          ! STANDARD SCF, NONGGA
          !
          ASSERT(nspin==ispin)
          if (pseudopot_present.and.core_density_setup) then
             call error_handler("Post_SCF: clean LDA.and.pseudopot.and.core_density branch first, AM")
          endif
          FPP_TIMER_START(loop)
          FPP_TIMER_START(sched)
          do while (more_grid_atom(vec_length, i, grdpts, grdwts))  ! loop over gridpoints
             FPP_TIMER_STOP(sched)
             FPP_TIMER_START(work)
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore

             rho=0.0_r8_kind
             dfdrho=0.0_r8_kind
             fxc=0.0_r8_kind

             vla=size(grdpts,1)

             grdwts(1:vla)=grdwts(1:vla)*&
               atomicweight(i,grdpts(1:vla,:))*&
               unique_atoms(i)%n_equal_atoms

             call orbital_calculate(grdpts(1:vla,1:3),&
                  vla,orbs_ob)
             call dtrace("Post_SCF orbital_calculate <")

             call density_calc(vla,rho,orbs_ob)

             ! performing charge and exc integration over the grid
             charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                  (grdwts(1:vla), 2, nspin))

             call start_timer(timer_gridph_functionals)

             call xc_functionals(vla, nspin, rho, GRDWTS=grdwts)

             call stop_timer(timer_gridph_functionals)
             FPP_TIMER_STOP(work)
             FPP_TIMER_START(sched)
          end do! loop over gridpoints
       endif sp_orb

    else
       !
       !  NLDA calculation
       !
       FPP_TIMER_START(loop)
       FPP_TIMER_START(sched)
       do while (more_grid_atom(vec_length, i, grdpts, grdwts)) ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
          FPP_TIMER_STOP(sched)
          FPP_TIMER_START(work)
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             if (spin_orbit_polarized) then
                !
                ! OPEN SHELL
                !
                ASSERT(nspin==2)
                rho=0.0_r8_kind
                dfdrho=0.0_r8_kind
                fxc=0.0_r8_kind

                vla = size(grdpts,1)
                call start_timer(timer_gridph_grdwts)
                grdwts(1:vla)=grdwts(1:vla)*&
                     atomicweight(i,grdpts(1:vla,:))*&
                     unique_atoms(i)%n_equal_atoms
                call stop_timer(timer_gridph_grdwts)
                call start_timer(timer_gridph_orbitals)
                call orbital_calculate(grdpts(1:vla,1:3),&
                     vla,orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
                call stop_timer(timer_gridph_orbitals)
                call start_timer(timer_gridph_density)

                call density_calc_nl(vla,rho,gamma,grarho,&
                     & orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads,&
                     & s_abs=s_abs,s_dens=s_dens,&
                     & gras_dens=gras_dens)
                ! determine "spin up" and "spin down" densities
                ! following the transformations:
                ! rho_up   = (rho + |s|)/2
                ! rho_down = (rho - |s|)/2
                rho_ud(1:vla,1) =  0.5_r8_kind*(rho(1:vla,1) + s_abs(1:vla)) ! spin up
                rho_ud(1:vla,2) =  0.5_r8_kind*(rho(1:vla,1) - s_abs(1:vla)) ! spin down
                where (abs(rho_ud).lt.1.0e-40_r8_kind)
                   rho_ud = 0.0_r8_kind
                end where

                ! determine gradients of "spin up" and "spin down" densities
                do xyz=X,Z
                   grarho_ud(1:vla,xyz,UP) = 0.5_r8_kind*(grarho(1:vla,xyz,1) + gras_dens(1:vla,xyz))
                   grarho_ud(1:vla,xyz,DN) = 0.5_r8_kind*(grarho(1:vla,xyz,1) - gras_dens(1:vla,xyz))
                enddo

                gamma(1:vla,UPUP) = SUM(grarho_ud(1:vla,X:Z,UP)**2, DIM=2)
                gamma(1:vla,DNDN) = SUM(grarho_ud(1:vla,X:Z,DN)**2, DIM=2)
                gamma(1:vla,UPDN) = SUM(grarho_ud(1:vla,X:Z,UP)*grarho_ud(1:vla,X:Z,DN), DIM=2)

                where (abs(gamma).lt.1.0e-40_r8_kind)
                   gamma = 0.0_r8_kind
                end where

                call stop_timer(timer_gridph_density)
                call start_timer(timer_gridph_functionals)

                ! performing charge and exc integration over the grid
                charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                     (grdwts(1:vla), 2, nspin))

                ! now calculate functionals
                call xc_functionals(vla,2,rho_ud,GAM=gamma,GRDWTS=grdwts)

                call stop_timer(timer_gridph_functionals)
             else
                !
                ! CLOSED SHELL
                !
                ASSERT(nspin==1)
                rho=0.0_r8_kind
                dfdrho=0.0_r8_kind
                fxc=0.0_r8_kind
                vla=size(grdpts,1)
                call start_timer(timer_gridph_grdwts)
                grdwts(1:vla)=grdwts(1:vla)*&
                     atomicweight(i,grdpts(1:vla,:))*&
                     unique_atoms(i)%n_equal_atoms
                call stop_timer(timer_gridph_grdwts)
                call start_timer(timer_gridph_orbitals)
                call orbital_calculate(grdpts(1:vla,1:3),&
                     vla,orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
                call stop_timer(timer_gridph_orbitals)
                call start_timer(timer_gridph_density)

                call density_calc_nl(vla,rho,gamma,grarho,&
                     & orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads)

                call stop_timer(timer_gridph_density)
                call start_timer(timer_gridph_functionals)
                ! performing charge and exc integration over the grid
                charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                     (grdwts(1:vla), 2, nspin))

                ! now calculate functionals
                call xc_functionals(vla,1,rho,GAM=gamma,GRDWTS=grdwts)

                call stop_timer(timer_gridph_functionals)
             endif! open / closed shell

          else
             !
             ! STANDARD POST SCF, GGA
             !
             ASSERT(nspin==ispin)
             rho=ZERO
             dfdrho=ZERO
             fxc=ZERO

             gamma=ZERO

             grarho = ZERO

             vla=size(grdpts,1)
             call start_timer(timer_gridph_grdwts)

             ! atomic value modified with partition and symmetry factors
             grdwts(1:vla)=grdwts(1:vla)*&
                  atomicweight(i,grdpts(1:vla,:))*unique_atoms(i)%n_equal_atoms

             call stop_timer(timer_gridph_grdwts)
             call start_timer(timer_gridph_orbitals)
             call orbital_calculate(grdpts(1:vla,1:3),&
                  vla,orbs_ob=orbs_ob,grads=orbs_grads)

             call stop_timer(timer_gridph_orbitals)
             call start_timer(timer_gridph_density)

             call density_calc_nl(vla,rho,gamma,grarho,orbs_ob,orbs_grads,tau=tau)

             call stop_timer(timer_gridph_density)
             call start_timer(timer_gridph_functionals)
             ! performing charge and exc integration over the grid
             charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
                  (grdwts(1:vla), 2, nspin))

             ! now calculate functionals
             call xc_functionals(vla, nspin, rho, GAM=gamma, GRDWTS=grdwts, TAU=tau)

             call stop_timer(timer_gridph_functionals)

          endif! options_spin_orbit
          FPP_TIMER_STOP(work)
          FPP_TIMER_START(sched)

       end do! loop over gridpoints
    end if
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)
    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(.false., send_grad_grid, send_mda_coeff)

    ! deallocate gradients of density
    deallocate(grarho, stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    deallocate(gamma, stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    deallocate(tau, stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)
  end subroutine post_scf_calc_xc_energy


#ifdef WITH_SECDER
#ifdef WITH_EXPERIMENTAL
  subroutine post_scf_calc_xc_en_and_gr(job,dervsw,nuc_grarho_imp)

    ! XXX
    ! output: dervs_xc for given call

    ! purpose : routine does the main calculations
    !           - calculation of prim. orbitals and gradients of prim. orbitals
    !           - calculation of becke weights
    !           - calculation of gradients of weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !           - calculation of LDA-gradients
    !           - calculation of LDA-derivatives
    !** End of interface *****************************************
    use time_module
    use timer_module
    use ph_cntrl
    use xc_func ! exc_*, xc_functionals
    use cpks_grid_utils, only: symWtsNGra, cartNGra, cartNDer
    use cpks_xc_resp   , only: xc_resp_lda
    use cpks_dens_calc , only: density_calc_ph_v2, fix_rho_grads
    use cpksdervs_matrices, only: cpks,n_cpks,i_cpks,cpks_Qxc_calculated, &
                                  comm_dervs_xc
    use integralpar_module, only: integralpar_2dervs
#ifdef FPP_TIMERS
    use calc3c_switches !,only: t_ == timers
#endif
    implicit none
    integer(i4_kind), intent(in)         :: job ! 0,1,2,3,4
    type(arrmat5),intent(inout),optional :: dervsw(:,:) ! contains the derivatives of
    type(arrmat4),intent(inout),optional :: nuc_grarho_imp(:) ! contains the derivatives of
    ! *** end of interface ***

    integer(i4_kind) :: i,ii,i_ea,i_ua,i_spin,i_vec,i_ma,i_m
    integer(i4_kind) :: k,l
    integer(i4_kind) :: ii_ua,ii_ea,k_gr,ii_ma
    integer(i4_kind) :: vla
    logical :: split_gradients, moving_grid, moving_atom
    logical :: moving_atom2
    logical :: send_grad_grid, send_mda_coeff

    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths

    ! lot of help_arrays and pointers
    real(kind=r8_kind) :: g(3)   ! cartesian gradient
    real(kind=r8_kind) :: d(3,3) ! cartesian second derivative matrix
    real(kind=r8_kind) :: fxc_help(vec_length),&
         & grarho(vec_length,X:Z,ispin),&
         & help_vec(vec_length), help_vec2(vec_length,ispin)
    real(kind=r8_kind) :: gamma(vec_length,1+2*(ispin-1))
    real(kind=r8_kind),pointer :: gr(:),nuc(:,:,:)
    logical:: expl, impl, sder
    integer(i4_kind) :: gb

    FPP_TIMER_DECL(tot)
    FPP_TIMER_DECL(awts)
    FPP_TIMER_DECL(orb)
    FPP_TIMER_DECL(dens)
    FPP_TIMER_DECL(xcfun)
    FPP_TIMER_DECL(cpksxc)
    FPP_TIMER_DECL(zero)
    FPP_TIMER_DECL(grua)
    FPP_TIMER_DECL(dervs_im)

    FPP_TIMER_START(tot)
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)

    split_gradients = options_split_gradients()

    ! compute explicit derivatives (gradients):
    expl = job == 10 .or. job == 20

    ! this job is one of the sec. der. jobs:
    sder = job >= 20

    ! compute implicit derivatives:
    impl = present(nuc_grarho_imp) &
              .and.(i_cpks.gt.n_cpks.or..not.cpks_Qxc_calculated)

#ifdef  FPP_TIMERS
    print*,MyID,"xc_en_and_gr(",job,"): present(dervsw)        =",present(dervsw)
    print*,MyID,"xc_en_and_gr(",job,"): present(nuc_grarho_imp)=",present(nuc_grarho_imp)
    print*,MyID,"xc_en_and_gr(",job,"): integralpar_2dervs     =",integralpar_2dervs
    print*,MyID,"xc_en_and_gr(",job,"): allocated(cpks)        =",allocated(cpks)
    print*,MyID,"xc_en_and_gr(",job,"): i_cpks,n_cpks          =",i_cpks,n_cpks
    print*,MyID,"xc_en_and_gr(",job,"): impl                   =",impl
    print*,MyID,"xc_en_and_gr(",job,"): cpks_Qxc_calculated    =",cpks_Qxc_calculated
    print*,MyID,"xc_en_and_gr(",job,"): comm_dervs_xc          =",comm_dervs_xc
#endif

    if(present(dervsw))then
      ! only present when computing explicit dervs:
      ASSERT(job==20)
    endif

   if( expl )then
    FPP_TIMER_START(zero)
    do i=1,N_moving_unique_atoms

       grad_xc(i)%m=0.0_r8_kind

       if (weight_grads .and. split_gradients) then
          grad_grid(i)%m=0.0_r8_kind
       end if
    end do
    FPP_TIMER_STOP(zero)
   endif

   if( expl .and. sder ) then
    FPP_TIMER_START(zero)
    do i=1,N_moving_unique_atoms
      do ii=1,N_moving_unique_atoms
        dervs_xc(i,ii)%m=0.0_r8_kind !set to zero with post_scf call
      enddo
    end do
    FPP_TIMER_STOP(zero)
   endif

   gb = 0
    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
    FPP_TIMER_START(loop)
    FPP_TIMER_START(sched)
    gridpoints: do while (more_grid_atom(vec_length, i, grdpts, grdwts))
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
          FPP_TIMER_START(work)
          FPP_TIMER_STOP(sched)
          i_m = unique_atoms(i)%moving_atom
          moving_grid = i_m /= 0 .and. weight_grads
          ASSERT(i_m==i)

          FPP_TIMER_START(zero)
          do ii=1,n_unique_atoms
             nuc_grarho(ii)%m=0.0_r8_kind
          enddo

          rho=0.0_r8_kind
          dfdrho=0.0_r8_kind

          if( job >= 21 )then
            ASSERT(allocated(cpks))
            ! not needed for expl?
            dvdrho=0.0_r8_kind
          endif

          fxc=0.0_r8_kind
          fxc_help=0.0_r8_kind
          gamma=0.0_r8_kind
          grarho = ZERO
          FPP_TIMER_STOP(zero)

          vla=size(grdpts,1)
          gb = gb + 1

          call start_timer(timer_gridph_grdwts)

          if(weight_grads) then ! in post_scf_calc_xc_en_and_gr

             FPP_TIMER_START(zero)
             do ii=1,n_unique_atoms
                graw(ii)%m=0.0_r8_kind
             enddo
             FPP_TIMER_STOP(zero)

             FPP_TIMER_START(awts)
             if( job == 20 ) then
               ASSERT(present(dervsw))
               help_vec(1:vla)=atomicweight_and_dervs(i,grdpts(1:vla,:),graw,dervsw)
             else
               help_vec(1:vla)=atomicweight_and_grad(i,grdpts(1:vla,:),graw)
             endif
             FPP_TIMER_STOP(awts)

             ! grdwts modified with symmetry factor
             grdwts(1:vla)=grdwts(1:vla)*unique_atoms(i)%n_equal_atoms

             ! modyfy atomic weights with partition scheme contribs
             do i_ua=1,n_unique_atoms
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
               do k_gr=1,3
                graw(i_ua)%m(1:vla,k_gr,i_ea)=&
                     graw(i_ua)%m(1:vla,k_gr,i_ea)*grdwts(1:vla)
               enddo
             enddo
             enddo

           if( job == 20 ) then
             do i_ua=1,n_unique_atoms
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

               do ii_ua=1,n_unique_atoms
               do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms

                  do k_gr=1,3
                  dervsw(i_ua,ii_ua)%m(1:vla,k_gr,i_ea,:,ii_ea)= &
                   dervsw(i_ua,ii_ua)%m(1:vla,k_gr,i_ea,:,ii_ea) &
                    *spread(grdwts(1:vla),2,3)
                  enddo

               enddo
               enddo
             enddo
             enddo
           endif

             grdwts(1:vla)=grdwts(1:vla)*help_vec(1:vla)
             ! weight modified with partition factor

          else ! gradients of weights will not be used
             grdwts(1:vla)=grdwts(1:vla)*&
                  atomicweight(i,grdpts(1:vla,:))*unique_atoms(i)%n_equal_atoms
          end if

          call stop_timer(timer_gridph_grdwts)
          call start_timer(timer_gridph_orbitals)

         FPP_TIMER_START(orb)
        if( job == 20 )then
          ! one order more (nuc_sec_der) for second derivatives:
          call orbital_calculate(grdpts(1:vla,1:3),&
               vla,orbs_ob=orbs_ob,grads=orbs_grads,&
               nuc_grads=nuc_grads,nuc_sec_der=nuc_sec_der)
               !                   -----------------------
        else
          call orbital_calculate(grdpts(1:vla,1:3),&
               vla,orbs_ob=orbs_ob,grads=orbs_grads,&
               nuc_grads=nuc_grads)
        endif

!       ! MAGIC: ad-scf measure to ensure translational invariance:
!       if( gb == 1 )then
!         WARN('fix_NUC_grads on grid host!')
!       endif
!       do ii  =1,size(nuc_grads,2) ! irreps
!             nuc_grads(i_m ,ii)%o = 0.0D0
!       enddo

        FPP_TIMER_STOP(orb)

          call stop_timer(timer_gridph_orbitals)
          call start_timer(timer_gridph_density)

         FPP_TIMER_START(dens)
         select case ( job )
         case ( 10 ) ! gradients in a gradient run
           ! only density and density gradient wrt nuc are req:
           call density_calc_ph(vla                           &
                               , rho,gamma,grarho             &
                               , nuc_grarho                   &
                               , orbs_ob,orbs_grads,nuc_grads &
                               )
           ! MAGIC: ad-scf measure to ensure translational invariance:
           if( gb == 1 )then
             WARN('fix_RHO_grads on grid host!')
           endif
           call fix_rho_grads(i_m,vla,nuc_grarho)

         case ( 20 ) ! gradients (and more?) in a sec-der run
           call density_calc_ph_v2(vla, ispin                 & ! post_scf_calc_xc_en_and_gr
                              , rho,gamma,grarho              & ! output
                              , nuc_grarho                    & ! output rho grads wrt nuc
                              , orbs_ob,orbs_grads,nuc_grads  & ! input
                              , phi_ob=phi_ob                 & ! output MO-orbitals
                              , nuc_dervsrho=nuc_dervsrho     & ! output rho dervs wrt nuc
                              , nuc_sec_der=nuc_sec_der       & ! input orbital dervs
                              )
         case ( 21, 24 ) ! CPKS RHS and CPKS FINAL ITERATION
           ! MO-orbitals ``phi_ob'' and their gradients ``phi_gr_nuc''
           ! will be later used to compute CPKS RHS:
           call density_calc_ph_v2(vla, ispin                 & ! post_scf_calc_xc_en_and_gr
                              , rho,gamma,grarho              & ! output
                              , nuc_grarho                    & ! output    density grads wrt nuc (cartesian)
                              , orbs_ob,orbs_grads,nuc_grads  & ! input
                              , phi_ob=phi_ob                 & ! output MO-orbitals
                              , phi_gr_nuc=phi_gr_nuc         & ! output MO-orbital grads wrt nuc
                              , rho_gr_nuc=rho_gr_nuc         & ! output    density grads wrt nuc (symm modes)
                              )
         case ( 22, 23 ) ! CPKS iterations in a sec-der run
           ! MO-orbitals ``phi_ob'' will be later used to compute the density response due to
           ! (intemediate) CPKS solutions U(i,a):
           call density_calc_ph_v2(vla, ispin                 & ! post_scf_calc_xc_en_and_gr
                              , rho,gamma,grarho              & ! output rho, gamma, rho grads
                              , nuc_grarho                    & ! output             rho grads wrt nuc
                              , orbs_ob,orbs_grads,nuc_grads  & ! input
                              , phi_ob=phi_ob                 & ! output MO-orbitals
                              )
         case default
           print *,'ERROR: dont know such job=',job
           ABORT('dont know such job')
         end select

!        ! MAGIC: ad-scf measure to ensure translational invariance:
!        if( gb == 1 )then
!          WARN('fix_RHO_grads on grid host!')
!        endif
!        call fix_rho_grads(i_m,vla,nuc_grarho)

         FPP_TIMER_STOP(dens)

          ! performing charge and exc integration over the grid
          charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:) &
                         *spread(grdwts(1:vla),2,ispin))

          call stop_timer(timer_gridph_density)
          call start_timer(timer_gridph_functionals)


          FPP_TIMER_START(xcfun)
          if( job == 10 .or. job == 20 )then
            call xc_functionals(vla,ispin,rho,fxc,dfdrho,GRDWTS=grdwts)
          endif

          if( job >= 21 )then
            call xc_functionals( vla,ispin,rho,fxc,dfdrho,GRDWTS=grdwts, & ! lda cpks_xc_main
                                 FRR=dvdrho )
          endif
          FPP_TIMER_STOP(xcfun)

          call stop_timer(timer_gridph_functionals)

          FPP_TIMER_START(cpksxc)
          select case ( job )
          case ( 21, 22, 23, 24 ) ! RHS and FINAL (due to -s1(i,j)/2 and u(i,a))
            ! convert weight gradients from cartesian to symm modes:
            call symWtsNGra(vla,graw,wts_gr_nuc)

            call xc_resp_lda(job,gb,vla                      &
                            , ispin                          &
                            , dvdrho=dvdrho                  &
                            , dfdrho=dfdrho                  &
                            , wts=grdwts                     &
                            , wts_gr_nuc=wts_gr_nuc          &
                            , phi=phi_ob                     &
                            , phi_gr_nuc=phi_gr_nuc          &
                            , rho_gr_nuc=rho_gr_nuc          & ! input: density gradient wrt symm nuc
                            , rho_gr_imp=rho_gr_imp          & ! output: density impl-grad wrt symm nuc
!                           , dervsrho_imp=dervsrho_imp      &
                            , eny_sd_imp=eny_sd_imp          &
                            )
          end select

          select case ( job )
          case ( 21, 24 ) ! RHS and FINAL (due to -s1(i,j)/2 and u(i,a))
            ! convert implicit density gradients from symm modes to cartesian:
            call cartNGra(vla,rho_gr_imp,nuc_grarho_imp)
          end select

          select case ( job )
          case ( 24 ) ! FINAL
            ! convert implicit energy derivaties from symm modes to cartesian:
            call cartNDer(eny_sd_imp,dervsrho_imp)
          end select
          FPP_TIMER_STOP(cpksxc)

          ! now start building the gradient
          do i_spin=1,ispin
             help_vec2(1:vla,i_spin)= dfdrho(1:vla,i_spin)*grdwts(1:vla)
          end do

       if( sder )then
          FPP_TIMER_START(zero)
          do i_ua=1,n_unique_atoms
           do ii_ua=1,n_unique_atoms
            dervs_xc_tmp(i_ua,ii_ua)%m=0.0_r8_kind
           enddo
          enddo
          FPP_TIMER_STOP(zero)
       endif

       FPP_TIMER_START(grua)
       grua: do i_ua=1,n_unique_atoms

             i_ma = unique_atoms(i_ua)%moving_atom
             moving_atom = i_ma /= 0

          if( .not.moving_atom .and. .not.moving_grid) cycle

             equals: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

                g = 0.0_r8_kind
                nuc => nuc_grarho(i_ua)%m(:,:,i_ea,:)

             calc_grad_xc: if( expl ) then !calculated with main_postscf call
                do i_spin=1,ispin
                   do k=1,3
                      do i_vec=1,vla
                        !                      weigt*(d fxc/d rho)*(d rho/ dx)
                         g(k)=g(k) + help_vec2(i_vec,i_spin) &
                                   * nuc_grarho(i_ua)%m(i_vec,k,i_ea,i_spin)
                      end do
                   end do
                enddo

                ASSERT(moving_grid)
                if (moving_atom) then
                   ! here contrib from this part of grid to all muv atoms

                   ! In case of grids attached to atoms (as in weight_grads)
                   ! there is no force from the grid acting on its
                   ! own host atom due to displacements:
                   !
                   !   dE/dN += sum[P] dE/dRho * ( dRho/dP * dP/dN )
                   !
                   ! where some grid points P fall out as they dont
                   ! move wrt nucleus N:

!                  if( i_ma /= i_m .or. .not. moving_grid )then
                     grad_xc(i_ma)%m(:,i_ea) =  grad_xc(i_ma)%m(:,i_ea) + g
!                  endif

                   ! There might be still some forces from the own grid
                   ! due to changing integration weights -- see below.
                   ! They will be added to this pointer:
                   gr => grad_xc(i_ma)%m(:,i_ea)
                endif

                if (weight_grads) then ! add gradients of weights

                   if (split_gradients) then
                      ABORT('for the moment')
                      ! store two separate cancelling contributions:
                      if( moving_atom ) grad_xc(i_ma)%m(:,i_ea) = &
                                        grad_xc(i_ma)%m(:,i_ea) + g
                      if( moving_grid ) grad_grid(i_m)%m(:,1)   = &
                                        grad_grid(i_m)%m(:,1)   - g
                      ! add forces due to weight grads into:
                      if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                   else
                      ! here the same contrib is substructed on center posesing this grid points
                      if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                                        grad_xc(i_m)%m(:,1) - g
                   endif

                   g = 0.0_r8_kind
                   do k=1,3
                      do i_vec=1,vla
                         g(k)=g(k)-graw(i_ua)%m(i_vec,k,i_ea)*fxc(i_vec)
                      end do
                   end do
                   if( moving_atom ) gr = gr + g
                endif
              endif calc_grad_xc

              ![[=== second derivatives ===========================
              sdr9: if( sder ) then
                dervsua: do ii_ua=1,n_unique_atoms

                 ii_ma = unique_atoms(ii_ua)%moving_atom
                 moving_atom2 = ii_ma /= 0

                 if( .not.moving_atom2 .and. .not.moving_grid) cycle

                dervsea: do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms

                dervs_xc_expl: if( expl ) then

                ! nuc_dervsrho only contrib

                d=0.0_r8_kind
                do i_spin=1,ispin
                   do k=1,3
                    do l=1,3
                      do i_vec=1,vla
                        !                       weight X dfdrho X (d rho / d x)  (1 -sym)
                        !                   nuc_dervsrho = -(d/dx)(nuc_grarho)
                         d(k,l)=d(k,l) - help_vec2(i_vec,i_spin)* &
                          nuc_dervsrho(i_ua,ii_ua)%m(i_vec,k,i_ea,l,ii_ea,i_spin)
                         ! nuc_dervsrho is proportional to occupation numbers similar
                         ! to nuc_graro contribs
                      end do
                     end do
                   end do
                enddo

                if (moving_atom2) &
                   dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) = &
                       dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) + d

                if (weight_grads.and.moving_grid) then
                        dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) = &
                              dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) - d
                endif

               endif dervs_xc_expl


                impl=present(nuc_grarho_imp) &
                          .and.(i_cpks.gt.n_cpks.or..not.cpks_Qxc_calculated)

                dervs_xc_imp: if(impl) then

                d=0.0_r8_kind
                dervs_imp_spin: do i_spin=1,ispin


                        !                       weight * dfdrho * (d rho / d x)
                   if(.not.cpks_Qxc_calculated)  then

                      ! s1 contribs calculated

                   do k=1,3
                     do l=1,3
                      d(k,l) = d(k,l) + &
                        sum( dvdrho(:vla,i_spin) * grdwts(:vla)             &
                           * nuc_grarho(i_ua)%m(:vla,k,i_ea,i_spin)        &
                           * ( nuc_grarho_imp(ii_ua)%m(:vla,l,ii_ea,i_spin) &
                              -nuc_grarho(ii_ua)%m(:vla,l,ii_ea,i_spin)     &
                             )                                              &
                           )     ! (4)

                      ! FIXME: second derivatives are only computed AFTER CPKS:
                      ! d(k,l)=d(k,l) + dervsrho_imp(i_ua,ii_ua)%m(k,i_ea,l,ii_ea,i_spin)
                     enddo
                   enddo
                   else
                   ! B contrib calculated
                   do k=1,3
                    do l=1,3
                      d(k,l) = d(k,l) + &
                        sum( dvdrho(:vla,i_spin) * grdwts(:vla) * &
                         nuc_grarho(i_ua)%m(:vla,k,i_ea,i_spin) * &
                         nuc_grarho_imp(ii_ua)%m(:vla,l,ii_ea,i_spin))      ! (4 else)

                      ! nuc_grarho_imp is proportional to occupations via HBH

                      d(k,l)=d(k,l) + dervsrho_imp(i_ua,ii_ua)%m(k,i_ea,l,ii_ea,i_spin)
!                      d(k,:)=d(k,:) + &
!                         2.0_r8_kind*sum( spread( dfdrho(1:vla,i_spin)*grdwts(1:vla),2,3) * &
!                              nuc_dervsrho_imp(i_ua,ii_ua)%m(1:vla,k,i_ea,:,ii_ea,i_spin) ,1 )
                    enddo
                   enddo
                   endif
                enddo dervs_imp_spin

                if(moving_atom2) &
                  dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) = &
                       dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) + d

                if (weight_grads.and.moving_grid) then
                        dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) = &
                              dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) - d
                endif
                endif dervs_xc_imp

            enddo dervsea
         enddo dervsua

      wts9: if( weight_grads ) then

         dervs_func: if( expl ) then ! calculated with main_postscf call
           dervsw_ua: do ii_ua=1,n_unique_atoms
             ii_ma = unique_atoms(ii_ua)%moving_atom
             moving_atom2 = ii_ma /= 0
          dervs_ea: do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms
           do k=1,3

           ! dervsw only term   (2 symm)
           ! graw term goes with (-) and so do dervsw term

            dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea)= &
               dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea) - &
                 sum( dervsw(i_ua,ii_ua)%m(1:vla,k,i_ea,:,ii_ea)* &
                      spread(fxc,2,3), 1)

           spins: do i_spin=1,ispin

            ! the first of two cross terms d/dx (w (d fxc/d x) + (d w / d x) * fxc)
            ! the nuc_grarho * dfdrho * w contrib goes with (+) sign and so does this term

            dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea)= &                     !(3a+ pair)
               dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea) + sum(&
                 spread(nuc_grarho(i_ua)%m(:vla,k,i_ea,i_spin)* &
                   dfdrho(1:vla,i_spin),2,3)*graw(ii_ua)%m(:vla,:,ii_ea), 1)

            ! for atom i_m nuc_grarho is wrong and is substituted
            ! with sum of grads over all other atoms

            dervs_xc_tmp(i_m,ii_ma)%m(k,1,:,ii_ea)= &                         !(3a- pair)
               dervs_xc_tmp(i_m,ii_ma)%m(k,1,:,ii_ea) - sum(&
                 spread(nuc_grarho(i_ua)%m(:vla,k,i_ea,i_spin)* &
                   dfdrho(1:vla,i_spin),2,3)*graw(ii_ua)%m(:vla,:,ii_ea), 1)

            dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea)= &
               dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea) + sum(&
                 nuc_grarho(ii_ua)%m(:vla,:,ii_ea,i_spin)* &
                   spread(dfdrho(1:vla,i_spin)*graw(i_ua)%m(:vla,k,i_ea),2,3), 1)

            dervs_xc_tmp(i_ma,i_m)%m(k,i_ea,:,1)= &
               dervs_xc_tmp(i_ma,i_m)%m(k,i_ea,:,1) - sum(&
                 nuc_grarho(ii_ua)%m(:vla,:,ii_ea,i_spin)* &
                   spread(dfdrho(1:vla,i_spin)*graw(i_ua)%m(:vla,k,i_ea),2,3), 1)

           enddo spins

           enddo
          enddo dervs_ea
         enddo dervsw_ua

        endif dervs_func

         gra_imp: if(present(nuc_grarho_imp) &
             .and.(i_cpks.gt.n_cpks.or..not.cpks_Qxc_calculated)) then

         gra_imp_ua: do ii_ua=1,n_unique_atoms
             ii_ma = unique_atoms(ii_ua)%moving_atom
             moving_atom2 = ii_ma /= 0
          gra_imp_ea: do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms
           do k=1,3
           gra_imp_spins: do i_spin=1,ispin

            dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea)= &
               dervs_xc_tmp(i_ma,ii_ma)%m(k,i_ea,:,ii_ea) - sum(&
                 nuc_grarho_imp(ii_ua)%m(:vla,:,ii_ea,i_spin)* &                        !(2)
                   spread(dfdrho(1:vla,i_spin)*graw(i_ua)%m(:vla,k,i_ea),2,3), 1)
            dervs_xc_tmp(i_ma,i_m)%m(k,i_ea,:,1)= &
               dervs_xc_tmp(i_ma,i_m)%m(k,i_ea,:,1) + sum(&
                 nuc_grarho_imp(ii_ua)%m(:vla,:,ii_ea,i_spin)* &                        !(1)
                   spread(dfdrho(1:vla,i_spin)*graw(i_ua)%m(:vla,k,i_ea),2,3), 1)

           enddo gra_imp_spins

           enddo
          enddo gra_imp_ea
         enddo gra_imp_ua

        endif gra_imp

   else  wts9
         WARN('weights of dervs_xc off')
   endif wts9
   endif sdr9
   !]]=== EOF second derivatives =======================


     enddo equals ! (i_ea)
   enddo grua ! (i_ua)
   FPP_TIMER_STOP(grua)


   FPP_TIMER_START(dervs_im)
   dervs_im: if( sder ) then
     do ii_ua=1,n_unique_atoms
      dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,:)=0.0_r8_kind
      dervs_xc_tmp(ii_ua,i_m)%m(:,:,:,1)=0.0_r8_kind
     enddo

     do i_ua=1,n_unique_atoms
      i_ma = unique_atoms(i_ua)%moving_atom
         moving_atom = i_ma /= 0
          if( .not.moving_atom .and. .not.moving_grid) cycle

           do ii_ua=1,n_unique_atoms
            ii_ma = unique_atoms(ii_ua)%moving_atom
            if( ii_ma.eq.0 .and. .not.moving_grid) cycle

            do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
             if(i_m.eq.i_ma.and.i_ea.eq.1) cycle
             do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms

               if(i_m.eq.ii_ma.and.ii_ea.eq.1) cycle

               dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,ii_ea)= &
                           dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,ii_ea) &
                                           -dervs_xc_tmp(i_ua,ii_ua)%m(:,i_ea,:,ii_ea)
               dervs_xc_tmp(i_ua,i_m)%m(:,i_ea,:,1)= &
                           dervs_xc_tmp(i_ua,i_m)%m(:,i_ea,:,1) &
                                           -dervs_xc_tmp(i_ua,ii_ua)%m(:,i_ea,:,ii_ea)
             enddo
            enddo
           enddo
     enddo

     dervs_xc_tmp(i_m,i_m)%m(:,1,:,1)=0.0_r8_kind

       do i_ua=1,n_unique_atoms
        do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
          if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
               dervs_xc_tmp(i_m,i_m)%m(:,1,:,1) = &
               dervs_xc_tmp(i_m,i_m)%m(:,1,:,1)   &
                    -dervs_xc_tmp(i_m,i_ua)%m(:,1,:,i_ea)
        enddo
       enddo


      do i_ua=1,n_unique_atoms
      do ii_ua=1,n_unique_atoms
       dervs_xc(i_ua,ii_ua)%m = dervs_xc(i_ua,ii_ua)%m &
                              + dervs_xc_tmp(i_ua,ii_ua)%m
      enddo
      enddo

   endif dervs_im
   FPP_TIMER_STOP(dervs_im)
   FPP_TIMER_STOP(work)
   FPP_TIMER_START(sched)

    enddo  gridpoints
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)

  if( sder )then
    print*,MyID//"dervs_xc calculated impl comm_dervs_xc",sum(abs(dervs_xc(1,1)%m)), &
         impl,cpks_Qxc_calculated,comm_dervs_xc
  endif

    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(.true., send_grad_grid, send_mda_coeff)

    FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
    print *,MyID,'xc_en_and_gr(',job,'): TIMING total          =',FPP_TIMER_VALUE(tot)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-a-weights    =',FPP_TIMER_VALUE(awts)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-orbitals     =',FPP_TIMER_VALUE(orb)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-density      =',FPP_TIMER_VALUE(dens)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-xc func      =',FPP_TIMER_VALUE(xcfun)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-cpks_xc      =',FPP_TIMER_VALUE(cpksxc)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING   |-xc_qh_expl =',FPP_TIMER_VALUE(t_xc_qh_expl)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING   |-xc_ab      =',FPP_TIMER_VALUE(t_xc_ab)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING   |-dervs_imp  =',FPP_TIMER_VALUE(t_dervs_imp)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-grua         =',FPP_TIMER_VALUE(grua)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-dervs_im     =',FPP_TIMER_VALUE(dervs_im)
    print *,MyID,'xc_en_and_gr(',job,'): TIMING |-zero         =',FPP_TIMER_VALUE(zero)

!   print*,'t_grid_unique      =', FPP_TIMER_VALUE(t_grid_uni)
!   print*,'t_dervs_weight     =', FPP_TIMER_VALUE(t_dervs_weight)
!   print*,'t_orbital_calculate=', FPP_TIMER_VALUE(t_orbital_calc)
!   print*,'t_density_calculate=', FPP_TIMER_VALUE(t_density_calc)
!   print*,'t_nuc_dervsrho     =', FPP_TIMER_VALUE(t_nuc_dervsrho)
!   print*,'t_ph_sums          =', FPP_TIMER_VALUE(t_ph_sums)
!   if(.not.expl) then
!   print*,'t_xc_qh_expl       =', FPP_TIMER_VALUE(t_xc_qh_expl)
    print*,'t_s1_imp_grarho    =', FPP_TIMER_VALUE(t_s1_imp_grarho)
!   print*,'t_xc_ab            =', FPP_TIMER_VALUE(t_xc_ab)
    print*,'t_init_mocalc      =', FPP_TIMER_VALUE(t_init_mocalc)
!   print*,'t_dervs_imp        =', FPP_TIMER_VALUE(t_dervs_imp)
!   else
    print*,'t_nuc3rd           =', FPP_TIMER_VALUE(t_nuc3rd)
!   endif
#endif

  end subroutine post_scf_calc_xc_en_and_gr
#endif /* ifdef WITH_EXPERIMENTAL */

#else
  subroutine post_scf_calc_xc_en_and_gr()
    ! HISTORIC VERSION, FOR DEVELOPMENT USE THE ABOVE ONE!
    ! purpose : routine does the main calculations
    !           - calculation of prim. orbitals and gradients of prim. orbitals
    !           - calculation of becke weights
    !           - calculation of gradients of weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !           - calculation of LDA-gradients
    !** End of interface *****************************************
    use time_module
    use timer_module
    use ph_cntrl
    use xc_func ! exc_*, xc_functionals
    implicit none
    integer(i4_kind) :: i,ii,i_ea,i_ua,i_spin,k,i_vec,i_ma,i_m
    integer(i4_kind) :: vla
    logical :: split_gradients, moving_grid, moving_atom
    ! lot of help_arrays and pointers
    real(kind=r8_kind) :: helpvar(3),fxc_help(vec_length),&
         & grarho(vec_length,X:Z,ispin),&
         & help_vec(vec_length), help_vec2(vec_length,ispin)
    real(kind=r8_kind) :: gamma(vec_length,1+2*(ispin-1))
    real(kind=r8_kind),pointer :: gr(:),nuc(:,:,:)
    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths
    logical :: send_grad_grid, send_mda_coeff

    ! prepare some timers
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)

    split_gradients = options_split_gradients()
    do i=1,N_moving_unique_atoms
       grad_xc(i)%m=0.0_r8_kind
       if (weight_grads .and. split_gradients) then
          grad_grid(i)%m=0.0_r8_kind
       end if
    end do
    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
    FPP_TIMER_START(loop)
    FPP_TIMER_START(sched)
    do while (more_grid_atom(vec_length, i, grdpts, grdwts))
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
       FPP_TIMER_STOP(sched)
       FPP_TIMER_START(work)

       i_m = unique_atoms(i)%moving_atom
       moving_grid = i_m /= 0 .and. weight_grads
          do ii=1,n_unique_atoms
             nuc_grarho(ii)%m=0.0_r8_kind
          end do
          rho=0.0_r8_kind
          dfdrho=0.0_r8_kind

          fxc=0.0_r8_kind
          fxc_help=0.0_r8_kind
          gamma=0.0_r8_kind
          grarho = ZERO
          vla=size(grdpts,1)

          call start_timer(timer_gridph_grdwts)
          if(weight_grads) then ! gradients of weights will be used
             do ii=1,n_unique_atoms
                graw(ii)%m=0.0_r8_kind
             end do
             help_vec(1:vla)=&
                  atomicweight_and_grad(i,grdpts(1:vla,:),graw)
             grdwts(1:vla)=grdwts(1:vla)*&
                  unique_atoms(i)%n_equal_atoms
             do i_ua=1,n_unique_atoms
                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   graw(i_ua)%m(1:vla,1,i_ea)=&
                        graw(i_ua)%m(1:vla,1,i_ea)*&
                        grdwts(1:vla)
                   graw(i_ua)%m(1:vla,2,i_ea)=&
                        graw(i_ua)%m(1:vla,2,i_ea)*&
                        grdwts(1:vla)
                   graw(i_ua)%m(1:vla,3,i_ea)=&
                        graw(i_ua)%m(1:vla,3,i_ea)*&
                        grdwts(1:vla)
                end do
             end do
             grdwts(1:vla)=grdwts(1:vla)*&
                  help_vec(1:vla)
          else ! gradients of weights will not be used
             grdwts(1:vla)=grdwts(1:vla)*&
                  atomicweight(i,grdpts(1:vla,:))*&
                  unique_atoms(i)%n_equal_atoms
          end if
          call stop_timer(timer_gridph_grdwts)
          call start_timer(timer_gridph_orbitals)
          call orbital_calculate(grdpts(1:vla,1:3),&
               vla,orbs_ob=orbs_ob,grads=orbs_grads,&
               nuc_grads=nuc_grads)
          call stop_timer(timer_gridph_orbitals)
          call start_timer(timer_gridph_density)
          call density_calc_ph(vla,rho,gamma,grarho,&
               & nuc_grarho,orbs_ob,orbs_grads,nuc_grads)

          ! performing charge and exc integration over the grid
          charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
               (grdwts(1:vla),2,ispin))

          call stop_timer(timer_gridph_density)
          call start_timer(timer_gridph_functionals)

          call xc_functionals(vla,ispin,rho,fxc,dfdrho,GRDWTS=grdwts)

          call stop_timer(timer_gridph_functionals)
          ! now start building the gradient
          do i_spin=1,ispin
             help_vec2(1:vla,i_spin)=dfdrho(1:vla,i_spin)*grdwts(1:vla)
          end do

          grua: do i_ua=1,n_unique_atoms
             i_ma = unique_atoms(i_ua)%moving_atom
             moving_atom = i_ma /= 0
          if( .not.moving_atom .and. .not.moving_grid) cycle
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                helpvar=0.0_r8_kind
!!!MF is that correct?
                nuc=>nuc_grarho(i_ua)%m(:,:,i_ea,:)

                do i_spin=1,ispin
                   do k=1,3
                      do i_vec=1,vla
                        ! dfdrho X grdweight X grad nuc
                         helpvar(k)=helpvar(k)+help_vec2(i_vec,i_spin)*nuc(i_vec,k,i_spin)
                      end do
                   end do
                enddo

                if (moving_atom) then
                   gr => grad_xc(i_ma)%m(:,i_ea)
                   ! here contrib from this part of grid to all muv atoms
                   gr =  gr + helpvar
                endif
                if (weight_grads) then ! add gradients of weights
                   if (split_gradients) then
                      if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                      if( moving_grid ) grad_grid(i_m)%m(:,1) = &
                                        grad_grid(i_m)%m(:,1) - helpvar
                   else
                      ! here the same contrib is substructed on center posesing this grid points
                      if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                                        grad_xc(i_m)%m(:,1) - helpvar
                   endif
                   helpvar=0.0_r8_kind
                   do k=1,3
                      do i_vec=1,vla
                         helpvar(k)=helpvar(k)-fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                      end do
                   end do
                   if( moving_atom ) gr = gr + helpvar
                end if
             end do
          enddo grua
       FPP_TIMER_STOP(work)
       FPP_TIMER_START(sched)
    end do
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)

    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(.true., send_grad_grid, send_mda_coeff)
  end subroutine post_scf_calc_xc_en_and_gr
#endif

#ifdef WITH_SECDER
  subroutine post_scf_calc_xc_en_and_gr_nl(job, dervsw, nuc_grarho_imp, nuc_gragamma_imp)
    ! purpose : routine does the main calculations
    !           - calculation of prim. orbitals,gradients and sec. derivatives
    !           - calculation of becke weights
    !           - calculation of gradients of weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !           - calculation of GGA-gradients
    !           - calcularions of derivatives not implemented yet
    use time_module
    use timer_module
    use ph_cntrl
    use xc_cntrl, only: xc_hcth_version, is_on, xc_mgga
    use xc_func ! exc_*, xc_functionals

    use cpksdervs_matrices, only: cpks, n_cpks, i_cpks, cpks_Qxc_calculated, &
         cpksalloc
    use integralpar_module, only: integralpar_2dervs
#ifdef FPP_TIMERS
    use calc3c_switches !,only: t_ == timers
#endif
    use eigen_data_module
    use density_calc_cpks ! to know interface of cpksdervs_xc

    implicit none

    integer(i4_kind), intent(in) :: job
    type(arrmat5), intent(inout), optional :: dervsw(:,:) ! contains the derivatives of
    type(arrmat4), intent(inout), optional :: nuc_grarho_imp(:) ! contains the derivatives of
    type(arrmat5), intent(inout), optional :: nuc_gragamma_imp(:) ! contains the derivatives of
    !** End of interface *****************************************

    integer(kind=i4_kind), parameter:: uup=1,ddn=2,udn=3,dup=4,udup=5,uddn=6
    integer(kind=i4_kind), parameter:: aa=1,bb=2,ab=3,ac=4,bc=5,cc=6

    integer(i4_kind) :: ii_ua,ii_ea,k_gr,alloc_stat

    integer(i4_kind) :: i,ii,i_ea,i_ua,i_spin,k,i_vec,i_ma,i_m
    integer(i4_kind) :: ii_ma
    integer(i4_kind) :: spin,i_irr

    integer(i4_kind) :: vla,xyz
    logical :: split_gradients, moving_grid, moving_atom
    logical :: moving_atom2
    logical :: imp_dervs
    logical :: send_grad_grid, send_mda_coeff

    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths

    ! now some helparrays follow
    real(kind=r8_kind) :: helpvar2(3,3)
    real(kind=r8_kind) :: helpvar(3),fxc_help(vec_length),&
         & grarho(vec_length,X:Z,ispin),&
         & help_vec(vec_length),&
         & wdfdrho(vec_length,ispin), & ! instead of help_vec2 ien_and_gr
         & help_xyz_ud(vec_length,X:Z,UP:DN)

    ! arrays needed for mgga
    real(kind=r8_kind),allocatable :: tau(:,:)      , & ! (vl,ISPIN)
                                      dfdtau(:,:)   , & ! derivatives of f with respect to tau (vl,ISPIN)
                                      wdfdtau(:,:)      ! derivatives of f with respect to tau * grdwts

    real(kind=r8_kind) :: gamma(vec_length,1+2*(ispin-1)),  &
         dfdgrarho(vec_length,2*ispin-1)

    real(kind=r8_kind),allocatable :: wdfdgrarho(:),wdf_drhodgamma(:), &
         wdf_dgammadgamma(:),grarho_derrho(:,:)

    ! now pointer on gradients of weights, gradients of the energies and gradients and 2nd
    ! gradients of the density follow ===> that means:
    ! grw     = [-d/dR]w
    ! gr      = [-d/dR]f
    ! nuc     = [-d/dR](drho/dr)
    ! nuc_tau = [-d/dR]tau
    real(kind=r8_kind),pointer :: gr(:),nuc(:,:,:),sec_der(:,:,:,:),nuc_tau(:,:,:)

    logical ::  postscf_call

    ! These were variables at the  debugging stage which turned on and
    ! off specific second derivative  terms. In production all of them
    ! are on:
    logical, parameter :: t1_dervs = .true.
    logical, parameter :: t2_dervs = .true.
    logical, parameter :: t3_dervs = .true.

    real(kind=r8_kind), parameter:: two=2.0_r8_kind

!!$    real(kind=r8_kind),pointer :: nuc_core(:,:,:), sec_der_core(:,:,:,:)


    split_gradients = options_split_gradients()

    !
    ! There  is a  problem  with  making the  operation  mode of  this
    ! subroutine depend  on the value of  present(dervsw). In gradient
    ! runs it remains unassociated, but  is declared as a plain (not a
    ! pointer), though  optional, argument of  this procedure.  FIXME:
    ! This is illegal in Fortran  95, I think.  In later standards the
    ! argument  is assumed  to  be "not  present"  e.g. with  GFortran
    ! 4.6. Do not rely on that:
    !
    postscf_call = (job == POST_SCF_CALL)

    imp_dervs=present(nuc_grarho_imp) &
         .and.(i_cpks.gt.n_cpks.or..not.cpks_Qxc_calculated)

    if (is_on(xc_mgga)) then
       allocate(tau(vec_length,ispin), stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler&
            ('Allocation (1) in post_scf_calc_xc_en_and_gr_nl failed')
       allocate(dfdtau(vec_length,ispin), stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler&
            ('Allocation (2) in post_scf_calc_xc_en_and_gr_nl failed')
       allocate(wdfdtau(vec_length,ispin), stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler&
            ('Allocation (3) in post_scf_calc_xc_en_and_gr_nl failed')
    end if

    ! initialize timers
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)


    init_grad_xc: if (postscf_call) then
       do i=1,N_moving_unique_atoms

          grad_xc(i)%m = 0.0

          if(integralpar_2dervs) then
             do ii=1,N_moving_unique_atoms
                dervs_xc(i,ii)%m=0.0_r8_kind !set to zero with post_scf call
             enddo
          endif

          if (weight_grads .and. split_gradients) then
             grad_grid(i)%m=0.0_r8_kind
          end if

       end do
    endif init_grad_xc

    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
    FPP_TIMER_START(loop)
    FPP_TIMER_START(sched)
    do while (more_grid_atom(vec_length, i, grdpts, grdwts)) ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
       FPP_TIMER_STOP(sched)
       FPP_TIMER_START(work)
       FPP_TIMER_START(t_grid_uni)
       i_m = unique_atoms(i)%moving_atom
       moving_grid = i_m /= 0 .and. weight_grads
          do ii=1,n_unique_atoms
             nuc_grarho(ii)%m=0.0_r8_kind
             if(is_on(xc_mgga)) then
                nuc_grad_tau(ii)%m=0.0_r8_kind
             end if
             if(nl_calc_ph)  then
                nuc_sec_derrho(ii)%m=0.0_r8_kind
             endif

          enddo
          rho=0.0_r8_kind
          dfdrho=0.0_r8_kind

          if(allocated(cpks)) then
             dvdrho=0.0_r8_kind !nl call was allocated in cpks_xc main
             if(nl_calc_ph) then
                df_drhodgamma=0.0_r8_kind
                df_dgammadgamma=0.0_r8_kind
             endif
          endif

          fxc=0.0_r8_kind
          fxc_help=0.0_r8_kind
          gamma=0.0_r8_kind
          grarho = ZERO
          if(nl_calc_ph) dfdgrarho=0.0_r8_kind
          vla=size(grdpts,1)

          if(nl_calc_ph) then
             allocate(wdf_drhodgamma(vla),wdfdgrarho(vla), &
                  wdf_dgammadgamma(vla),grarho_derrho(vla,3), &
                  stat=cpksalloc(150))
             ASSERT(cpksalloc(150).eq.0)
             MEMLOG(vla*6)
          endif

          if(allocated(cpks).and.nl_calc_ph) then


             allocate(graphi(ssym%n_irrep,ssym%n_spin),stat=cpksalloc(160))
             ASSERT(cpksalloc(160).eq.0)
             MEMLOG(size(graphi))

             do spin=1,ssym%n_spin
                do i_irr=1,ssym%n_irrep
                   allocate(graphi(i_irr,spin)%m(vla,3,size(eigvec(i_irr)%m,1),ssym%partner(i_irr)), &
                        stat=cpksalloc(143))
                   ASSERT(cpksalloc(143).eq.0)
                   MEMLOG(size(graphi(i_irr,spin)%m))
                enddo
             enddo
          endif
#ifdef WITH_CORE_DENS
!!$          if (pseudopot_present.and.core_density_setup) then
!!$             do ii=1,n_unique_atoms
!!$                nuc_grarho_core(ii)%m=0.0_r8_kind
!!$                nuc_sec_derrho_core(ii)%m=0.0_r8_kind
!!$             end do
!!$          end if
#endif

          call start_timer(timer_gridph_grdwts)


          if(weight_grads) then
             FPP_TIMER_START(t_dervs_weight)
             do ii=1,n_unique_atoms
                graw(ii)%m=0.0_r8_kind
             end do

             if(present(dervsw).and.integralpar_2dervs) then
                help_vec(1:vla)=atomicweight_and_dervs(i,grdpts(1:vla,:),graw,dervsw) !main_postscf call
             else
                help_vec(1:vla)=atomicweight_and_grad(i,grdpts(1:vla,:),graw)
             endif
!!!           print*,'atomicweight_and_grad'
!!!           MEMSET(0)

             grdwts(1:vla)=grdwts(1:vla)*unique_atoms(i)%n_equal_atoms

             first_ua: do i_ua=1,n_unique_atoms
                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

                   do k_gr=1,3
                      graw(i_ua)%m(1:vla,k_gr,i_ea)=&
                           graw(i_ua)%m(1:vla,k_gr,i_ea)*grdwts(1:vla)
                   enddo

                   dervsw_fin: if(present(dervsw).and.integralpar_2dervs) then
                      do ii_ua=1,n_unique_atoms
                         do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms
                            do k_gr=1,3
                               dervsw(i_ua,ii_ua)%m(1:vla,k_gr,i_ea,:,ii_ea)= &
                                    dervsw(i_ua,ii_ua)%m(1:vla,k_gr,i_ea,:,ii_ea) &
                                    *spread(grdwts(1:vla),2,3)
                            enddo
                         enddo
                      enddo
                   endif dervsw_fin

                end do
             enddo first_ua

             grdwts(1:vla)=grdwts(1:vla)* help_vec(1:vla)
             FPP_TIMER_STOP(t_dervs_weight)
          else
             grdwts(1:vla)=grdwts(1:vla)*&
                  atomicweight(i,grdpts(1:vla,:))*unique_atoms(i)%n_equal_atoms
          endif

          call stop_timer(timer_gridph_grdwts)
          call start_timer(timer_gridph_orbitals)

          FPP_TIMER_START(t_orbital_calc)
          if(nl_calc_ph.or.(postscf_call.and.integralpar_2dervs)) then

             if(postscf_call.and.nl_calc_ph.and.integralpar_2dervs) then
                call orbital_calculate(grdpts(:vla,:3),&
                     vla,orbs_ob=orbs_ob,grads=orbs_grads,&
                     nuc_grads=nuc_grads,nuc_sec_der=nuc_sec_der, &
                     nuc_3rd_der=nuc_3rd_der)
             else
                call orbital_calculate(grdpts(1:vla,1:3),&
                     vla,orbs_ob=orbs_ob,grads=orbs_grads,&
                     nuc_grads=nuc_grads,nuc_sec_der=nuc_sec_der)
             endif

          else
             call orbital_calculate(grdpts(1:vla,1:3),&
                  vla,orbs_ob=orbs_ob,grads=orbs_grads,&
                  nuc_grads=nuc_grads)
          endif

          FPP_TIMER_STOP(t_orbital_calc)
          call stop_timer(timer_gridph_orbitals)
          call start_timer(timer_gridph_density)

          FPP_TIMER_START(t_density_calc)
          !
          ! FIXME: the comments on branch purpose may be wrong ...
          !
          if(nl_calc_ph .and. postscf_call .and. integralpar_2dervs) then
             !
             ! Case 1. GGA gradients and second energy derivatives.
             !
             call density_calc_ph_nl(vla, rho, gamma, grarho, nuc_grarho, &
                  orbs_ob, orbs_grads, nuc_grads, &
                  nuc_sec_derrho=nuc_sec_derrho, & ! (grad rho)^lambda
                  nuc_dervs_grarho=nuc_dervs_grarho, &
                  nuc_sec_der=nuc_sec_der, & ! input orbital dervs required for dervs and nonlda
                  nuc_3rd_der=nuc_3rd_der, & ! input orbital dervs required for dervs and nonlda
                  nuc_dervsrho=nuc_dervsrho) ! dervs only nuclear rho derivatives

          elseif (nl_calc_ph) then
             !
             ! GGA.   This branch  appears to  include  plain gradient
             ! calculations  and  CPKS  iterations  where  the  linear
             ! density  response to  nuclei displacemnts  are computed
             ! repeatedly.
             !
             if(allocated(cpks) ) then
                !
                ! Case 2. CPKS iterations.
                !
                select case(job)

                case(23)
                   call density_calc_ph_nl(vla, rho, gamma, grarho, nuc_grarho, &
                        orbs_ob, orbs_grads, nuc_grads, &
                        nuc_sec_der=nuc_sec_der,  & ! input orbital dervs for dervs and nonlda
                        graphi=graphi)

                case default
                   call density_calc_ph_nl(vla, rho, gamma, grarho, nuc_grarho, &
                        orbs_ob, orbs_grads, nuc_grads, &
                        nuc_sec_derrho=nuc_sec_derrho, & ! nonlda only (grad rho)^lambda
                        nuc_sec_der=nuc_sec_der, & ! input orbital dervs for dervs and nonlda
                        graphi=graphi)
                end select
             else ! that is not allocated(cpks) ...
                !
                ! Case 3.   This branch  is entered in  regular (M)GGA
                ! gradient runs.
                !
                if(is_on(xc_mgga)) then
                   !
                   ! 3.1  Meta GGA case,  need kinetic  energy density
                   ! TAU in addition to density and its gradients:
                   !
                   call density_calc_ph_nl(vla,rho,gamma,grarho,nuc_grarho, &
                        orbs_ob,orbs_grads,nuc_grads, &
                        nuc_sec_derrho=nuc_sec_derrho, & ! nonlda only
                        nuc_sec_der=nuc_sec_der, & ! input orbital dervs required for dervs and nonlda
                        tau=tau,nuc_grad_tau=nuc_grad_tau) ! mgga mode only
                else
                   !
                   ! 3.2  Plain   GGA  case,  need   density  and  its
                   ! gardients:
                   !
                   call density_calc_ph_nl(vla, rho, gamma, grarho, nuc_grarho, &
                        orbs_ob, orbs_grads, nuc_grads, &
                        nuc_sec_derrho=nuc_sec_derrho, & ! nonlda only
                        nuc_sec_der=nuc_sec_der) ! input orbital dervs required for dervs and nonlda
                end if
             endif
          elseif (postscf_call .and. integralpar_2dervs) then ! lda dervs
             !
             ! Case 4. LDA gradients and second energy derivatives.
             !
             call density_calc_ph_nl(vla, rho, gamma, grarho, &
                  nuc_grarho, &
                  orbs_ob, orbs_grads, nuc_grads, &
                  nuc_dervsrho=nuc_dervsrho, & ! output rho dervs
                  nuc_sec_der=nuc_sec_der)     ! input orbital dervs
          else ! no dervs lda
             !
             ! Case 5. The simplest case, LDA gradients.
             !
             call density_calc_ph(vla, rho, gamma, grarho, &
                  nuc_grarho, &
                  orbs_ob, orbs_grads, nuc_grads)
          endif

#ifdef WITH_CORE_DENS
!!$          if (pseudopot_present.and.core_density_setup) then
!!$            call fit_core_calculate(grdpts(1:vec_length_act,1:3),&
!!$                 vec_length_act,fcts=orbs_ob_core, &
!!$                 grads=grads_core,&
!!$                 nuc_grads=nuc_grads_core, &
!!$                 sec_nuc_ders=nuc_sec_der_core)
!!$          end if
#endif

          ! performing charge and exc integration over the grid
          charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:) &
               *spread(grdwts(1:vla),2,ispin))

#ifdef WITH_CORE_DENS
!!$          if (pseudopot_present.and.core_density_setup) then
!!$             call fitted_core_density_calc(vec_length_act,&
!!$                  rho=rho_core, grad=grad_rho_core,&
!!$                  nuc_grad=nuc_grarho_core,&
!!$                  sec_nuc_der=nuc_sec_derrho_core,&
!!$                  fcts=orbs_ob_core, grads=grads_core, &
!!$                  nuc_grads=nuc_grads_core,&
!!$                  sec_nuc_ders=nuc_sec_der_core)
!!$                if (ispin==1)then
!!$                  gamma(:,1) = gamma(:,1) + grad_rho_core(:,1,1)* &
!!$                                            grad_rho_core(:,1,1)&
!!$                                          + grad_rho_core(:,2,1)* &
!!$                                            grad_rho_core(:,2,1)&
!!$                                          + grad_rho_core(:,3,1)* &
!!$                                            grad_rho_core(:,3,1)
!!$                else
!!$                  gamma(:,1)= gamma(:,1) + grad_rho_core(:,1,1)* &
!!$                                           grad_rho_core(:,1,1)&
!!$                                         + grad_rho_core(:,2,1)* &
!!$                                           grad_rho_core(:,2,1)&
!!$                                         + grad_rho_core(:,3,1)* &
!!$                                           grad_rho_core(:,3,1)
!!$                  gamma(:,2)= gamma(:,2) + grad_rho_core(:,1,2)* &
!!$                                           grad_rho_core(:,1,2)&
!!$                                         + grad_rho_core(:,2,2)* &
!!$                                           grad_rho_core(:,2,2)&
!!$                                         + grad_rho_core(:,3,2)* &
!!$                                           grad_rho_core(:,3,2)
!!$                  gamma(:,3)= gamma(:,3) + grad_rho_core(:,1,1)* &
!!$                                           grad_rho_core(:,1,2)&
!!$                                         + grad_rho_core(:,2,1)* &
!!$                                           grad_rho_core(:,2,2)&
!!$                                         + grad_rho_core(:,3,1)* &
!!$                                           grad_rho_core(:,3,2)
!!$                end if
!!$                grarhox(:,:) = grarhox(:,:) + grad_rho_core(:,1,:)
!!$                grarhoy(:,:) = grarhoy(:,:) + grad_rho_core(:,2,:)
!!$                grarhoz(:,:) = grarhoz(:,:) + grad_rho_core(:,3,:)
!!$          end if
#endif
          FPP_TIMER_STOP(t_density_calc)
          call stop_timer(timer_gridph_density)

          if(integralpar_2dervs) then
             do i_ua=1,n_unique_atoms
                do ii_ua=1,n_unique_atoms
                   dervs_xc_tmp(i_ua,ii_ua)%m=0.0_r8_kind
                enddo
             enddo
          endif

          call start_timer(timer_gridph_functionals)



!!!        print*,'call cpksdervs_xc'
!!!        MEMSET(0)
          dervs_cpksgrads: if(allocated(cpks).and.nl_calc_ph) then ! non lda cpks call
#ifdef WITH_EXPERIMENTAL
             ABORT('recompile w/o -DWITH_EXPERIMENTAL')
#else

             call xc_functionals(vla,ispin,rho,fxc,dfdrho, &
                  gamma,dfdgrarho, &
                  GRDWTS=grdwts,FRR=dvdrho, &
                  FRG=df_drhodgamma, &
                  FGG=df_dgammadgamma) ! xc_functionals nonlda in cpks_xc_main

!                 print*,'outside cpksdervs_xc dfdgrarho',sum(dfdgrarho(:vla,1)),sum(dfdgrarho(:vla,2)),sum(dfdgrarho(:vla,3))
!                 print*,shape(grarho)

             call cpksdervs_xc(vec_length,vla,ispin,dvdrho,orbs_ob,phi_ob,grdwts, &      !non lda (1)
                  nuc_grads,dfdrho,nuc_grarho,graw,i_m,rho=rho, &
                  nuc_grarho_imp=nuc_grarho_imp, &
                  nuc_gragamma_imp=nuc_gragamma_imp, &
                  dervsrho_imp=dervsrho_imp, &
                  dervs_grarho_imp=dervs_grarho_imp, &
                  nuc_sec_der=nuc_sec_der, &
                  nuc_sec_derrho=nuc_sec_derrho, &
                  grarho=grarho, &
                  graphi=graphi, &
                  df_drhodgamma=df_drhodgamma, &
                  df_dgammadgamma=df_dgammadgamma, &
                  dfdgrarho=dfdgrarho)


             do spin=1,ispin
                do i_irr=1,ssym%n_irrep
                   MEMLOG(-size(graphi(i_irr,spin)%m))
                   deallocate(graphi(i_irr,spin)%m,stat=cpksalloc(143))
                   ASSERT(cpksalloc(143).eq.0)
                   cpksalloc(143)=1
                enddo
             enddo

             MEMLOG(-size(graphi))
             deallocate(graphi,stat=cpksalloc(160))
             ASSERT(cpksalloc(160).eq.0)
             cpksalloc(160)=1
#endif /* of ifdef WITH_EXPERIMENTAL */
          elseif(nl_calc_ph) then
             if (.not.is_on(xc_mgga)) then                   ! GGA case
                call xc_functionals(vla,ispin,rho,fxc,dfdrho,gamma,dfdgrarho, &
                     GRDWTS=grdwts)
             else                                            ! MGGA case
                call xc_functionals(vla,ispin,rho,fxc,dfdrho,gamma,dfdgrarho, &
                     GRDWTS=grdwts,tau=tau,ft=dfdtau)
             end if

          elseif(allocated(cpks)) then ! lda cpks call
#ifdef WITH_EXPERIMENTAL
             ABORT('call exp LDA vers')
#else
             call xc_functionals( vla,ispin,rho,fxc,dfdrho,GRDWTS=grdwts, & ! lda cpks_xc_main
                  FRR=dvdrho )

             call cpksdervs_xc(vec_length,vla,ispin,dvdrho,orbs_ob,phi_ob,grdwts, &    ! lda (2)
                  nuc_grads,dfdrho,nuc_grarho,graw,i_m,rho=rho, &
                  nuc_grarho_imp=nuc_grarho_imp, &
                  dervsrho_imp=dervsrho_imp)
#endif /* of ifdef WITH_EXPERIMENTAL */
          else ! lda no dervs
             call xc_functionals(vla,ispin,rho,fxc,dfdrho,GRDWTS=grdwts)

          endif dervs_cpksgrads
!!!print*,'cpksdervs_xc done'
!!!MEMSET(0)

          call stop_timer(timer_gridph_functionals)

          ! t1_dervs=.true.
          ! t2_dervs=.true.
          ! t3_dervs=.true.

!!!  make nuc_grarho and nuc_sec_derrho correct
          do i_spin=1,ispin
             nuc_grarho(i_m)%m(:vla,:3,1,i_spin)=zero
             if(nl_calc_ph) nuc_sec_derrho(i_m)%m(:vla,:3,:3,1,i_spin)=zero
             if(is_on(xc_mgga)) nuc_grad_tau(i_m)%m(:vla,:3,1,i_spin)=zero
             do i_ua=1,n_unique_atoms
                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
                   nuc_grarho(i_m)%m(:vla,:,1,i_spin)= &
                        nuc_grarho(i_m)%m(:vla,:,1,i_spin)-nuc_grarho(i_ua)%m(:vla,:,i_ea,i_spin)
                   if(nl_calc_ph)  nuc_sec_derrho(i_m)%m(:vla,:3,:3,1,i_spin)= &
                        nuc_sec_derrho(i_m)%m(:vla,:3,:3,1,i_spin)  &
                        - nuc_sec_derrho(i_ua)%m(:vla,:3,:3,i_ea,i_spin)
                   if(is_on(xc_mgga)) then
                     nuc_grad_tau(i_m)%m(:vla,:,1,i_spin) = &
                        nuc_grad_tau(i_m)%m(:vla,:,1,i_spin) - nuc_grad_tau(i_ua)%m(:vla,:,i_ea,i_spin)
                   end if
                enddo
             enddo
          enddo

!!! no dvdrho in post_scf call
#if 0 /* rho grarho */
          if(imp_dervs.and.ispin.eq.2.and.nl_calc_ph) then
             do i_spin=1,ispin
                if(cpks_Qxc_calculated) then
                   print*,'rho grarho',i_spin,sum(rho(:vla,i_spin)), -sum(nuc_grarho_imp(2)%m(:vla,1,1,i_spin))
                else
                   print*,'rho grarho',i_spin,sum(rho(:vla,i_spin)),sum(nuc_grarho(2)%m(:vla,1,1,i_spin)) &
                        -sum(nuc_grarho_imp(2)%m(:vla,1,1,i_spin))
                endif
             enddo
          endif
#endif
!!!            if(imp_dervs.and.ispin.eq.2.and.nl_calc_ph) then
!!!            print*,'2x dfdrho(up)^lambda',&
!!!             sum(dfdrho(:vla,up)),       &
!!!             sum( dvdrho(:vla,upup)*nuc_grarho(1)%m(:vla,1,1,up)     &
!!!                 +dvdrho(:vla,updn)*nuc_grarho(1)%m(:vla,1,1,dn) ) &
!!!                 +two*sum(df_drhodgamma(:vla,uup)* &
!!!                  sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2) )   &
!!!                 +sum(df_drhodgamma(:vla,udup)* &
!!!                  sum(grarho(:vla,:3,3-up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2) ) &
!!!                 +two*sum(df_drhodgamma(:vla,dup)* &
!!!                  sum(grarho(:vla,:3,dn)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2) )   &
!!!                 +sum(df_drhodgamma(:vla,udup)* &
!!!                  sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2) )
!!!            endif

!!!            if(imp_dervs.and.ispin.eq.2.and.nl_calc_ph) &
!!!            print*,'2x dfdrho(dn)^lambda',&
!!!             sum(dfdrho(:vla,dn)),       &
!!!             sum( dvdrho(:vla,dndn)*nuc_grarho(1)%m(:vla,1,1,dn)     &
!!!                 +dvdrho(:vla,updn)*nuc_grarho(1)%m(:vla,1,1,up) ) &
!!!                 +two*sum(df_drhodgamma(:vla,ddn)* &
!!!                  sum(grarho(:vla,:3,dn)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2) )   &
!!!                 +sum(df_drhodgamma(:vla,uddn)* &
!!!                  sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2) ) &
!!!                 +two*sum(df_drhodgamma(:vla,udn)* &
!!!                  sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2) )   &
!!!                 +sum(df_drhodgamma(:vla,uddn)* &
!!!                  sum(grarho(:vla,:3,dn)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2) )

!!!            if(imp_dervs.and.ispin.eq.1) &
!!!            print*,'2x dfdrho^lambda',&
!!!             sum(dfdrho(:vla,up)*nuc_grarho(1)%m(:vla,1,1,up)),     &
!!!             sum(dfdrho(:vla,up)),       &
!!!             sum( dvdrho(:vla,up)*nuc_grarho(1)%m(:vla,1,1,up)      )

!!!            if(imp_dervs.and.ispin.eq.1.and.nl_calc_ph) &
!!!            print*,'2x dfdgrarho^lamda', sum(dfdgrarho(:vla,up)),&
!!!             sum(df_drhodgamma(:vla,uup)*nuc_grarho(1)%m(:vla,1,1,up)) &
!!!            +sum(df_dgammadgamma(:vla,aa)*2* &
!!!              sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2))

!!!            if(imp_dervs.and.ispin.eq.2.and.nl_calc_ph) &
!!!            print*,'2x dfdgrarho^lamda', sum(dfdgrarho(:vla,up)),&
!!!             sum(df_drhodgamma(:vla,uup)*nuc_grarho(1)%m(:vla,1,1,up)) &
!!!            +sum(df_drhodgamma(:vla,udn)*nuc_grarho(1)%m(:vla,1,1,dn)) &
!!!            +sum(df_dgammadgamma(:vla,aa)*2* &
!!!              sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2)) &
!!!            +sum(df_dgammadgamma(:vla,ac)*1* &
!!!              sum(grarho(:vla,:3,dn)*nuc_sec_derrho(1)%m(:vla,1,:3,1,up),2)) &
!!!            +sum(df_dgammadgamma(:vla,ab)*2* &
!!!              sum(grarho(:vla,:3,dn)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2)) &
!!!            +sum(df_dgammadgamma(:vla,ac)*1* &
!!!              sum(grarho(:vla,:3,up)*nuc_sec_derrho(1)%m(:vla,1,:3,1,dn),2))

          if(integralpar_2dervs) call dervs_to_tmp() !1

          ! now start building the gradient
          restricted: if(ispin==1) then
             i_spin=1
             FPP_TIMER_START(t_ph_sums)
             wdfdrho(1:vla,1)=dfdrho(1:vla,i_spin)*grdwts(1:vla)

!-----------------------------------------------------------------------+
! Calculation of first part of gamma terms:                             |
!                                                                       |
!      2 df_dgamma * drho_dr * w                                        |
!-----------------------------------------------------------------------+
             if(nl_calc_ph) then
                wdfdgrarho(:vla)=dfdgrarho(:vla,i_spin)* grdwts(:vla)
                do xyz=X,Z
                   help_xyz_ud(:vla,xyz,UP) =  &
                        wdfdgrarho(:vla)*2.0_r8_kind*grarho(:vla,xyz,UP)
                enddo
             endif

!-----------------------------------------------------------------------+
! Calculation of first part of tau terms:                               |
!         df_dtau * w                                                   |
!-----------------------------------------------------------------------+
             if(is_on(xc_mgga)) then
               wdfdtau(1:vla,1)=dfdtau(1:vla,1)*grdwts(1:vla)
             end if

             grua: do i_ua=1,n_unique_atoms
                i_ma = unique_atoms(i_ua)%moving_atom
                moving_atom = i_ma /= 0
                if( .not.moving_atom .and. .not.moving_grid) cycle

                equals: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   helpvar=0.0_r8_kind
                   nuc=>nuc_grarho(i_ua)%m(:vla,:,i_ea,:)

                   if(nl_calc_ph) sec_der=>nuc_sec_derrho(i_ua)%m(:,:,:,i_ea,:)
                   if(is_on(xc_mgga)) nuc_tau=>nuc_grad_tau(i_ua)%m(:vla,:,i_ea,:)

!-----------------------------------------------------------------------+
! Complete calculation of rho terms:                                    |
!                                                                       |
! gr = df_drho * nuc * w                                                |
!-----------------------------------------------------------------------+
                   calc_grad_xc: if(postscf_call) then !calculated with main_postscf call
                      if(t1_dervs) then
                         do k=1,3
                            do i_vec=1,vla
                               helpvar(k)=helpvar(k)+ nuc(i_vec,k,1)*wdfdrho(i_vec,1) !!!(1)->(4,6,7)
                            end do
                         end do
                      endif

!-----------------------------------------------------------------------+
! Complete calculation of gamma terms to gradients                      |
!                                                                       |
!-----------------------------------------------------------------------+
                      if(nl_calc_ph) then
                         if(t2_dervs) then
                            do xyz=X,Z
                               do i_vec=1,vla
                                  helpvar(xyz)=helpvar(xyz)&
                                       +help_xyz_ud(i_vec,X,1)*sec_der(i_vec,xyz,X,1)&
                                       +help_xyz_ud(i_vec,Y,1)*sec_der(i_vec,xyz,Y,1)&
                                       +help_xyz_ud(i_vec,Z,1)*sec_der(i_vec,xyz,Z,1)
                               enddo
                            enddo
                         endif
                      endif

!-----------------------------------------------------------------------+
! Complete calculation of tau terms:                                    |
!                                                                       |
! gr = df_dtau * nuc_tau * w                                            |
!-----------------------------------------------------------------------+
                   if (is_on(xc_mgga)) then
                      do k=1,3
                         do i_vec=1,vla
                            helpvar(k)=helpvar(k)+wdfdtau(i_vec,1)*nuc_tau(i_vec,k,1)
                         end do
                      end do
                   end if

                      if( moving_atom )then
                         gr => grad_xc(i_ma)%m(:,i_ea)
                         gr =  gr + helpvar
                      endif

                      nucgr_weights: if (weight_grads) then
                         if (split_gradients) then
                            if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                            if( moving_grid ) grad_grid(i_m)%m(:,1) = &
                                 grad_grid(i_m)%m(:,1) - helpvar
                         else
                            if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                                 grad_xc(i_m)%m(:,1) - helpvar
                         endif

                         if(t3_dervs) then  !t3
                            helpvar=0.0_r8_kind
                            do k=1,3
                               do i_vec=1,vla  !!!!(3)->(1,2,3)
                                  helpvar(k)=helpvar(k)-fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                               enddo ! results in contribs (1,2,3)
                            end do
                            if( moving_atom ) gr = gr + helpvar
                         endif

                      endif nucgr_weights

                   endif calc_grad_xc


                   calc_dervs_xc:if(integralpar_2dervs) then


                      if(.false..and.i_ua.eq.2.and.i_ua.ne.i_m.and.imp_dervs) then
                         print*,'x dfdrho^lambda',sum(dfdrho(:vla,1)), &
                              sum(dvdrho(:vla,1)*nuc_grarho(i_ua)%m(:vla,1,1,1) &
                              +two*sum(nuc_sec_derrho(i_ua)%m(:vla,1,:3,1,1)*grarho(:vla,:3,1),2)*df_drhodgamma(:vla,1))
                         print*,'x dfdgrarho',sum(dfdgrarho(:vla,1)), sum(&
                              two*sum(nuc_sec_derrho(i_ua)%m(:vla,1,:3,1,1)*grarho(:vla,:3,1),2)*df_dgammadgamma(:vla,1)+ &
                              df_drhodgamma(:vla,1)*nuc_grarho(i_ua)%m(:vla,1,1,1) )
!!!             sum_rho=sum_rho+sum(grarho(:vla,:,1))
!!!             sum_grarho=sum_grarho+sum(nuc_sec_derrho(ii_ua)%m(:vla,:3,:3,1,1))
                      endif

                   endif calc_dervs_xc

                enddo equals
             enddo grua


             FPP_TIMER_STOP(t_ph_sums)


          else restricted          ! i.e. spinpolarized case

!-----------------------------------------------------------------------+
! Calculation of first part of rho terms:                               |
!     UP: df_drhoa * w                                                  |
!     DN: df_drhob * w                                                  |
!-----------------------------------------------------------------------+
             wdfdrho(1:vla,1)=dfdrho(1:vla,1)*grdwts(1:vla)
             wdfdrho(1:vla,2)=dfdrho(1:vla,2)*grdwts(1:vla)

!-----------------------------------------------------------------------+
! Calculation of first part of gamma terms:                             |
!     UP:                                                               |
!      (2 df_dgammaaa * drhoa_dr + df_dgammaab * drhob_dr) * w          |
!     DN:                                                               |
!      (2 df_dgammabb * drhob_dr + df_dgammaab * drhoa_dr) * w          |
!-----------------------------------------------------------------------+
             if(t2_dervs.and.nl_calc_ph) then
                do xyz=X,Z
                   help_xyz_ud(1:vla,xyz,UP) = (2.0_r8_kind*&
                        & dfdgrarho(1:vla,1)*grarho(1:vla,xyz,1)+&
                        & dfdgrarho(1:vla,3)*grarho(1:vla,xyz,2))*&
                        & grdwts(1:vla)
                   help_xyz_ud(1:vla,xyz,DN) = (2.0_r8_kind*&
                        & dfdgrarho(1:vla,2)*grarho(1:vla,xyz,2)+&
                        & dfdgrarho(1:vla,3)*grarho(1:vla,xyz,1))*&
                        & grdwts(1:vla)
                enddo
             endif

!-----------------------------------------------------------------------+
! Calculation of first part of tau terms - spin polarized:              |
!     UP: df_dtaua * w                                                  |
!     DN: df_dtaub * w                                                  |
!-----------------------------------------------------------------------+
             if(is_on(xc_mgga)) then
               wdfdtau(1:vla,1)=dfdtau(1:vla,1)*grdwts(1:vla)
               wdfdtau(1:vla,2)=dfdtau(1:vla,2)*grdwts(1:vla)
             end if

             do i_ua=1,n_unique_atoms
                i_ma = unique_atoms(i_ua)%moving_atom
                moving_atom = i_ma /= 0
                if( .not.moving_atom .and. .not.moving_grid) cycle
                do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

                   helpvar=0.0_r8_kind
                   nuc=>nuc_grarho(i_ua)%m(:,:,i_ea,:)
                   if(nl_calc_ph) sec_der=>nuc_sec_derrho(i_ua)%m(:,:,:,i_ea,:)
                   if(is_on(xc_mgga)) nuc_tau=>nuc_grad_tau(i_ua)%m(:,:,i_ea,:)

!-----------------------------------------------------------------------+
! Complete calculation of rho terms:                                    |
!                                                                       |
! gr = (df_drhoa * nuca + df_drhob * nucb) * w                          |
!-----------------------------------------------------------------------+
                   do i_spin=1,ispin
                      do k=1,3
                         do i_vec=1,vla
                            helpvar(k)=helpvar(k)+wdfdrho(i_vec,i_spin)*nuc(i_vec,k,i_spin)
                         end do
                      end do
                   enddo

!-----------------------------------------------------------------------+
! Complete calculation of gamma terms to gradients                      |
!                                                                       |
! gr = ((2 df_dgammaaa * drhoa_dr + df_dgammaab * drhob_dr) * sec_der_a |
!   +(2 df_dgammabb * drhob_dr + df_dgammaab * drhoa_dr) * sec_der_b)*w |
!                                                                       |
! with:         sec_der_a = [-d/dR]^T drhoa_dr                          |
!               sec_der_b = [-d/dR]^T drhob_dr                          |
!-----------------------------------------------------------------------+
                   if (t2_dervs .and. nl_calc_ph) then
                      do xyz = X, Z
                         helpvar(xyz) = helpvar(xyz) &
                              + sum(help_xyz_ud(:vla, X:Z, UP:DN) * sec_der(:vla, xyz, X:Z, UP:DN))
#if 0
                         ! below variant seems not to work with IFC version > 11.1 most likely due to
                         ! internal compiler thresholds. BUG seems not do occur for RKS path
                         ! FIXME: break down and simplify post_scf_calc_xc_en_and_gr_nl
                         do i_vec=1,vla
                            helpvar(xyz)=helpvar(xyz)&
                                 +help_xyz_ud(i_vec,X,UP)*sec_der(i_vec,xyz,X,UP)&
                                 +help_xyz_ud(i_vec,Y,UP)*sec_der(i_vec,xyz,Y,UP)&
                                 +help_xyz_ud(i_vec,Z,UP)*sec_der(i_vec,xyz,Z,UP)
                            helpvar(xyz)=helpvar(xyz)&
                                 +help_xyz_ud(i_vec,X,DN)*sec_der(i_vec,xyz,X,DN)&
                                 +help_xyz_ud(i_vec,Y,DN)*sec_der(i_vec,xyz,Y,DN)&
                                 +help_xyz_ud(i_vec,Z,DN)*sec_der(i_vec,xyz,Z,DN)
                         enddo
#endif
                      enddo
                   endif

!-----------------------------------------------------------------------+
! Complete calculation of tau terms:                                    |
!                                                                       |
! gr = (df_dtaua * nuc_grad_taua + df_dtaub * nuc_grad_taub) * w        |
!-----------------------------------------------------------------------+
                   if (is_on(xc_mgga)) then
                     do i_spin=1,ispin
                        do k=1,3
                           do i_vec=1,vla
                              helpvar(k)=helpvar(k)+wdfdtau(i_vec,i_spin)*nuc_tau(i_vec,k,i_spin)
                           end do
                        end do
                     enddo
                   end if

                   if( moving_atom )then
                      gr => grad_xc(i_ma)%m(:,i_ea)
                      gr =  gr + helpvar
                   endif
                   if (weight_grads) then
                      if (split_gradients) then
                         if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                         if( moving_grid ) grad_grid(i_m)%m(:,1) = &
                              grad_grid(i_m)%m(:,1) - helpvar
                      else
                         if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                              grad_xc(i_m)%m(:,1) - helpvar
                      endif

                      if(t3_dervs) then
                         helpvar=0.0_r8_kind
                         do k=1,3
                            do i_vec=1,vla
                               helpvar(k)=helpvar(k)-fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                            end do
                         end do
                         if( moving_atom ) gr = gr + helpvar
                      endif
                   end if

                enddo
             end do
          endif restricted

          if(integralpar_2dervs)  then
             call tmp_to_dervs()
             if(postscf_call.or.imp_dervs) then
                if (postscf_call) then
                   ASSERT(.not.imp_dervs)
                endif
                if (imp_dervs) then
                   ASSERT(.not.postscf_call)
                endif
                !
                ! Tell to compute explicit or implicit derivatives:
                !
                call to_dervs_xc(postscf_call)
             endif
          endif

          if(nl_calc_ph) then
             MEMLOG(-vla*6)
             deallocate(wdf_drhodgamma,wdfdgrarho,wdf_dgammadgamma,grarho_derrho, &
                  stat=cpksalloc(150))
             ASSERT(cpksalloc(150).eq.0)
             cpksalloc(150)=1
          endif
          FPP_TIMER_STOP(work)                                                    !
          FPP_TIMER_START(sched)
    end do
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)
    FPP_TIMER_STOP(t_grid_uni)
#ifdef FPP_TIMERS
    print *, MyID, 'ph_point ',i
#endif

    if (is_on(xc_mgga)) then
       deallocate(tau     , stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       deallocate(dfdtau  , stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       deallocate(wdfdtau , stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if


#ifdef FPP_TIMERS
    print*, MyID, 't_grid_unique           =', FPP_TIMER_VALUE(t_grid_uni)
    print*, MyID, 't_dervs_weight          =', FPP_TIMER_VALUE(t_dervs_weight)
    print*, MyID, 't_orbital_calculate     =', FPP_TIMER_VALUE(t_orbital_calc)
    print*, MyID, 't_density_calculate     =', FPP_TIMER_VALUE(t_density_calc)
    print*, MyID, 't_help_nuc_arr2         =', FPP_TIMER_VALUE(t_help_nuc_arr2)
    print*, MyID, 't_graphi_sec_der        =', FPP_TIMER_VALUE(t_graphi_sec_der)
    print*, MyID, 't_cross_nuc_dervsrho    =', FPP_TIMER_VALUE(t_cross_nuc_dervsrho)
    print*, MyID, 't_nuc3rd                =', FPP_TIMER_VALUE(t_nuc3rd)
    print*, MyID, 't_cross_nuc_dervs_grarho=', FPP_TIMER_VALUE(t_cross_dervs_grarho)
    print*, MyID, 't_nuc_dervsrho          =', FPP_TIMER_VALUE(t_nuc_dervsrho)
    print*, MyID, 't_ph_sums               =', FPP_TIMER_VALUE(t_ph_sums)
    if(.not.postscf_call) then
    print*, MyID, 't_xc_qh_expl            =', FPP_TIMER_VALUE(t_xc_qh_expl)
    print*, MyID, 't_s1_imp_grarho         =', FPP_TIMER_VALUE(t_s1_imp_grarho)
    print*, MyID, 't_xc_ab                 =', FPP_TIMER_VALUE(t_xc_ab)
    print*, MyID, 't_init_mocalc           =', FPP_TIMER_VALUE(t_init_mocalc)
    print*, MyID, 't_dervs_imp             =', FPP_TIMER_VALUE(t_dervs_imp)
    endif
#endif

    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(.true., send_grad_grid, send_mda_coeff)

!!!print*,'end cpks_xc_main'
!!!MEMSET(0)

  contains

  subroutine to_dervs_xc(explicit)
    implicit none
    logical, intent(in) :: explicit
    ! *** end of interface ***

    integer(i4_kind) :: i_ua, i_ma, i_ea, ii_ua, ii_ma, ii_ea
    integer(i4_kind) :: k, kk

    grua: do i_ua = 1, n_unique_atoms
       i_ma = unique_atoms(i_ua)%moving_atom
       moving_atom = i_ma /= 0

       if (.not. moving_atom .and. .not. moving_grid) cycle

       equals: do i_ea = 1, unique_atoms(i_ua)%n_equal_atoms
          dervsw_ua: do ii_ua = 1, n_unique_atoms
             ii_ma = unique_atoms(ii_ua)%moving_atom
             moving_atom2 = ii_ma /= 0
             dervs_ea: do ii_ea = 1, unique_atoms(ii_ua)%n_equal_atoms

                if (explicit) then
                   !
                   ! Explicit derivatives:
                   !

             if (t3_dervs) then
                ! (3.1f):
                forall (k = 1:3, kk = 1:3) &
                     dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &  ! (1)<-(3)
                     dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) - &
                     sum(dervsw(i_ua, ii_ua)%m(:vla, k, i_ea, kk, ii_ea) * fxc(:vla))
             endif

             if (t1_dervs) then
                do i_spin = 1, ispin
                   forall (k = 1:3, kk = 1:3) &
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = & ! (3a+ pair)
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                        sum(nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * dfdrho(:vla, i_spin) * graw(ii_ua)%m(:vla, kk, ii_ea))
                enddo
             endif

             if (t3_dervs) then
                do i_spin = 1, ispin
                   forall (k = 1:3, kk = 1:3) &  ! (3.2f)
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                        sum(nuc_grarho(ii_ua)%m(:vla, kk, ii_ea, i_spin) * dfdrho(:vla, i_spin) * graw(i_ua)%m(:vla, k, i_ea))
                enddo
             endif

            if(t2_dervs.and.nl_calc_ph)  then    ! d/dN w
               do i_spin = 1, ispin
                  ! forall ( k = 1:3, kk = 1:3) ! (5)<-(2)
                  do kk = 1, 3
                     do k = 1, 3
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &
                             dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                             two * sum(dfdgrarho(:vla, i_spin) * graw(ii_ua)%m(:vla, kk, ii_ea) * &
                             sum(grarho(:vla, :3, i_spin) * nuc_sec_derrho(i_ua)%m(:vla, k, :3, i_ea, i_spin), 2))
                     enddo
                  enddo
                  if (ispin .eq. 2) then
                     ! forall (k = 1:3, kk = 1:3) ! (5)<-(2)
                     do kk = 1, 3
                        do k = 1, 3
                           dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &
                                dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                                sum(dfdgrarho(:vla, updn) * graw(ii_ua)%m(:vla, kk, ii_ea) * &
                                sum(grarho(:vla, :3, 3 - i_spin) * nuc_sec_derrho(i_ua)%m(:vla, k, :3, i_ea, i_spin), 2) )
                        enddo
                     enddo
                  endif
               enddo
            endif

            if(nl_calc_ph.and.t3_dervs) then          !!! (3.3f)
               do i_spin = 1, ispin
                  ! forall( k = 1:3, kk = 1:3) ! (3)<-(3)
                  do kk = 1, 3
                     do k = 1, 3
                        dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &
                             dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                             two * sum(dfdgrarho(:vla, i_spin) * graw(i_ua)%m(:vla, k, i_ea) * &
                             sum(grarho(:vla, :3, i_spin) * nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea,i_spin), 2))
                     enddo
                  enddo
                  if (ispin .eq. 2) then
                     ! forall ( k = 1:3, kk = 1:3)
                     do kk = 1, 3
                        do k = 1, 3
                           dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) = &
                                dervs_xc(i_ma, ii_ma)%m(k, i_ea, kk, ii_ea) + &
                                sum(dfdgrarho(:vla, 3) * graw(i_ua)%m(:vla, k, i_ea) * &
                                sum(grarho(:vla, :3, 3 - i_spin) * nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin), 2))
                        enddo
                     enddo
                  endif
               enddo
             endif

    else ! that is .not. explicit
       !
       ! Implicit derivatives:
       !
              if(t1_dervs)then
                do i_spin=1,ispin
                forall(k=1:3) &
                dervs_xc(i_ua,ii_ua)%m(k,i_ea,:3,ii_ea)= &
                   dervs_xc(i_ua,ii_ua)%m(k,i_ea,:3,ii_ea) &
                 + dervsrho_imp(i_ua,ii_ua)%m(k,i_ea,:3,ii_ea,i_spin)
                enddo
              endif

            if(t2_dervs.and.nl_calc_ph) then
              do i_spin=1,ispin
              forall(k=1:3) & ! [d/dN nuc_sec_der]  (2) !restricted
              dervs_xc(i_ua,ii_ua)%m(k,i_ea,:3,ii_ea)= &
                 dervs_xc(i_ua,ii_ua)%m(k,i_ea,:3,ii_ea) &
               + dervs_grarho_imp(i_ua,ii_ua)%m(k,i_ea,:,ii_ea,i_spin)
              enddo
              endif



      if (t3_dervs) then
         do i_spin = 1, ispin
            forall (k = 1:3, kk = 1:3) & ! (3.2i)
                 dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                 dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) - &
                 sum(nuc_grarho_imp(ii_ua)%m(:vla, kk, ii_ea, i_spin) * dfdrho(:vla, i_spin) * graw(i_ua)%m(:vla, k, i_ea))
         enddo
      endif

      if(nl_calc_ph.and.t3_dervs) then              !!! (3.3i)
       do i_spin = 1, ispin
          ! forall(k = 1:3, kk = 1:3)
          do kk = 1, 3
             do k = 1, 3
                dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                     dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) - &
                     two * sum(sum(nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) * grarho(:vla, :3, i_spin), 2) * &
                     dfdgrarho(:vla, i_spin) * graw(i_ua)%m(:vla, k, i_ea))
             enddo
          enddo

          if (ispin .eq. 2) then
             ! forall (k = 1:3, kk = 1:3)
             do kk = 1, 3
                do k = 1, 3
                   dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                        dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) - &
                        sum(sum(nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) * grarho(:vla, :3, 3 - i_spin), 2) * &
                        dfdgrarho(:vla, 3) * graw(i_ua)%m(:vla, k, i_ea))
                enddo
             enddo
          endif
      enddo
      endif

         if(.not.cpks_Qxc_calculated) then

           if (t1_dervs) then
            do i_spin = 1, ispin
               ! forall (k = 1:3, kk = 1:3) & ! (6)<-(1)
               do kk = 1, 3
                  do k = 1, 3
                     dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                          dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                          sum(grdwts(:vla) * dvdrho(:vla, i_spin) * nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * &
                          (nuc_grarho_imp(ii_ua)%m(:vla, kk, ii_ea, i_spin) - &
                          nuc_grarho(ii_ua)%m(:vla, kk, ii_ea, i_spin)))
                  enddo
               enddo
              if (ispin .eq. 2) then
                 ! forall (k = 1:3, kk = 1:3)
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            sum(grdwts(:vla) * dvdrho(:vla, updn) * nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * &
                            (nuc_grarho_imp(ii_ua)%m(:vla, kk, ii_ea, 3 - i_spin) - &
                            nuc_grarho(ii_ua)%m(:vla, kk, ii_ea, 3 - i_spin)))
                    enddo
                 enddo
              endif
            enddo
            endif


             if(nl_calc_ph.and.t1_dervs) then   !!!(7)<-(1) d/dg dfdrho
              do i_spin=1,ispin
                 ! forall (k = 1:3, kk = 1:3) &
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * df_drhodgamma(:vla, i_spin) * nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * &
                            sum(grarho(:vla, :3, i_spin) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) - &
                            nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin)), 2))
                       enddo
                    enddo
              enddo
              if(ispin.eq.2) then
               ! forall (k = 1:3, kk = 1:3)
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, up) * &
                            (df_drhodgamma(:vla, udup) * &
                            sum(grarho(:vla, :3, dn) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up)  &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, up)), 2) &
                            + df_drhodgamma(:vla, udup) * &
                            sum(grarho(:vla, :3, up) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, dn)  &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, dn)), 2)))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, up) * &
                            df_drhodgamma(:vla, dup) * &
                            sum(grarho(:vla, :3, dn) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea,dn)  &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, dn)), 2))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, dn) * &
                            df_drhodgamma(:vla, udn) * &
                            sum(grarho(:vla, :3, up) * &
                            ( nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up)  &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, up)), 2))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, dn) * &
                            (df_drhodgamma(:vla, uddn) * &
                            sum(grarho(:vla, :3, up) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, dn)  &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, dn)), 2) &
                            + df_drhodgamma(:vla, uddn) * &
                            sum(grarho(:vla, :3, dn) * &
                            (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up)  &
                            -nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, up)), 2)))
                    enddo
                 enddo
                ! end forall
               endif
             endif


            if(nl_calc_ph.and.t2_dervs) then  !!!(9)<-(2) [d/dr df_dgamma]
             if(nl_calc_ph.and.t2_dervs) then   !!!(11)<-(2) [grarho^lambda]
              do i_spin = 1, ispin
                 ! forall (k = 1:3, kk = 1:3)
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * dfdgrarho(:vla, i_spin) * &
                            sum(nuc_sec_derrho(i_ua)%m(:vla, k, :3, i_ea, i_spin)* &
                            ( nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) &
                            - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) ), 2))
                    enddo
                 enddo
                 if (ispin .eq. 2) then
                    ! forall (k = 1:3, kk = 1:3)
                    do kk = 1, 3
                       do k = 1, 3
                          dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                               dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                               sum(grdwts(:vla) * dfdgrarho(:vla,updn) * &
                               sum(nuc_sec_derrho(i_ua)%m(:vla, k, :3, i_ea, i_spin) * &
                               (nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, 3 - i_spin) &
                               - nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, 3 - i_spin)), 2))
                       enddo
                    enddo
                 endif
              enddo
             endif

             do i_spin = 1, ispin
                ! forall(k=1:3,kk=1:3) &
                do kk = 1, 3
                   do k = 1, 3
                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)- &
                           two*sum( df_drhodgamma(:vla,i_spin)*(nuc_grarho(ii_ua)%m(:vla,kk,ii_ea,i_spin) &
                           -nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,i_spin) )*grdwts(:vla)* &
                           sum(grarho(:vla,:3,i_spin)*nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1)
                   enddo
                enddo
             enddo

             if (ispin .eq. 2) then
                do i_spin = 1, ispin
                   ! forall (k = 1:3, kk = 1:3)
                   do kk = 1, 3
                      do k = 1, 3
                         dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)- &
                              two*sum( df_drhodgamma(:vla,2+i_spin)*(nuc_grarho(ii_ua)%m(:vla,kk,ii_ea,3-i_spin) &
                              -nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,3-i_spin) )*grdwts(:vla)* &
                              sum(grarho(:vla,:3,i_spin)*nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1)

                         dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)- &
                              sum( df_drhodgamma(:vla,4+i_spin)*(nuc_grarho(ii_ua)%m(:vla,kk,ii_ea,i_spin) &
                              -nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,i_spin) )*grdwts(:vla)* &
                              sum(grarho(:vla,:3,3-i_spin)*nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1)

                         dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)- &
                              sum( df_drhodgamma(:vla,7-i_spin)*(nuc_grarho(ii_ua)%m(:vla,kk,ii_ea,3-i_spin) &
                              -nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,3-i_spin) )*grdwts(:vla)* &
                              sum(grarho(:vla,:3,3-i_spin)*nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1)
                      enddo
                   enddo
                   ! end forall
                enddo
             endif


             do i_spin=1,ispin
                ! forall (k = 1:3, kk = 1:3) ! t2_dervs [d/dg dfdgrarho]
                do kk = 1, 3
                   do k = 1, 3
                      dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                           dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                           4.0_r8_kind * sum(df_dgammadgamma(:vla, i_spin) * grdwts(:vla) * &
                           sum(nuc_sec_derrho(i_ua)%m(:vla, k, :3, i_ea, i_spin) &
                           *grarho(:vla, :3, i_spin), 2) * &
                           sum((nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin)) * &
                           grarho(:vla, :3, i_spin), 2), 1)
                   enddo
                enddo
             enddo

             if(ispin.eq.2) then
             do i_spin=1,ispin
             ! forall (k = 1:3, kk = 1:3)
                do kk = 1, 3
                   do k = 1, 3
                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           4.0_r8_kind*sum( df_dgammadgamma(:vla,ab)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin))* &
                           grarho(:vla,:3,3-i_spin),2), 1)

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           2.0_r8_kind*sum(df_dgammadgamma(:vla,3+i_spin)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,3-i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin))* &
                           grarho(:vla,:3,i_spin),2), 1)
                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           2.0_r8_kind*sum(df_dgammadgamma(:vla,3+i_spin)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin))* &
                           grarho(:vla,:3,3-i_spin),2), 1)

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           2.0_r8_kind*sum( df_dgammadgamma(:vla,6-i_spin)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,3-i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin))* &
                           grarho(:vla,:3,3-i_spin),2), 1)

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           2.0_r8_kind*sum( df_dgammadgamma(:vla,3+i_spin)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin))* &
                           grarho(:vla,:3,i_spin),2), 1)

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           sum( df_dgammadgamma(:vla,cc)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,3-i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin))* &
                           grarho(:vla,:3,3-i_spin),2), 1)

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                           sum( df_dgammadgamma(:vla,cc)*grdwts(:vla)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                           *grarho(:vla,:3,3-i_spin),2)* &
                           sum(( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin) &
                           -  nuc_sec_derrho(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin))* &
                           grarho(:vla,:3,i_spin),2), 1)
                   enddo
                enddo
             ! end forall
             enddo
             endif
             endif !t2_dervs

         else

           if(t1_dervs) then
             do i_spin=1,ispin
             forall(k=1:3,kk=1:3) &
              dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)= &
              dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+sum( &
                 grdwts(:vla)*dvdrho(:vla,i_spin)* &
                        nuc_grarho(i_ua)%m(:vla,k,i_ea,i_spin) &
                           * nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,i_spin),1)
             if (ispin .eq. 2) &
                  forall (k = 1:3, kk = 1:3) &
                  dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                  dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                  sum(grdwts(:vla) * dvdrho(:vla, updn) * &
                  nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * &
                  nuc_grarho_imp(ii_ua)%m(:vla, kk, ii_ea, 3 - i_spin))
            enddo
           endif


             if (nl_calc_ph .and. t1_dervs) then ! (7)<-(1) df_drho_dgamma * gamma^lambda
              do i_spin=1,ispin
                 ! forall (k = 1:3, kk = 1:3) &
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * df_drhodgamma(:vla, i_spin) * &
                            nuc_grarho(i_ua)%m(:vla, k, i_ea, i_spin) * &
                            sum(grarho(:vla, :3, i_spin) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, i_spin), 2))
                    enddo
                 enddo
              enddo

              if (ispin .eq. 2) then
                 ! forall (k = 1:3, kk = 1:3)
                 do kk = 1, 3
                    do k = 1, 3
                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua,ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, up) * &
                            (df_drhodgamma(:vla, udup) * &
                            sum(grarho(:vla, :3, dn) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up), 2) &
                            + df_drhodgamma(:vla, udup) * &
                            sum(grarho(:vla, :3, up) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, dn), 2)))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, up) * &
                            df_drhodgamma(:vla, dup) * &
                            sum(grarho(:vla, :3, dn) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, dn), 2))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            two * sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, dn) * &
                            df_drhodgamma(:vla, udn) * &
                            sum(grarho(:vla, :3, up) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up), 2))

                       dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) = &
                            dervs_xc(i_ua, ii_ua)%m(k, i_ea, kk, ii_ea) + &
                            sum(grdwts(:vla) * nuc_grarho(i_ua)%m(:vla, k, i_ea, dn) * &
                            (df_drhodgamma(:vla, uddn) * &
                            sum(grarho(:vla, :3, up) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, dn), 2) &
                            + df_drhodgamma(:vla, uddn) * &
                            sum(grarho(:vla, :3, dn) * &
                            nuc_gragamma_imp(ii_ua)%m(:vla, kk, :3, ii_ea, up), 2)))
                    enddo
                 enddo
                 ! end forall
               endif !ispin.eq.2
             endif !t1_dervs

            if(nl_calc_ph.and.t2_dervs) then ! [d/dr dfdgrarho]

            do i_spin=1,ispin
             ! forall( k = 1:3, kk = 1:3)   ! grarho^lambda
               do kk = 1, 3
                  do k = 1, 3
                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                          +two*sum( grdwts(:vla)*dfdgrarho(:vla,i_spin)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin)* &
                          nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin),2),1)
                  enddo
               enddo

             if (ispin .eq. 2) then
                ! forall(k=1:3,kk=1:3) &
                do kk = 1, 3
                   do k = 1, 3
                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                           +sum(grdwts(:vla)* dfdgrarho(:vla,updn)* &
                           sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin)* &
                           nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin),2),1)
                   enddo
                enddo
             endif

             ! forall(k=1:3,kk=1:3) &
             do kk = 1, 3
                do k = 1, 3
                   dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                        + two*sum(nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,i_spin)* &
                        grdwts(:vla)* df_drhodgamma(:vla,i_spin)*sum(grarho(:vla,:3,i_spin)* &
                        nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1 )
                enddo
             enddo

             if (ispin .eq. 2) then
              ! forall(k=1:3,kk=1:3)
                do kk = 1, 3
                   do k = 1, 3
                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                           + two*sum(df_drhodgamma(:vla,2+i_spin)*nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,3-i_spin)* &
                           grdwts(:vla)* sum(grarho(:vla,:3,i_spin)* &
                           nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1 )

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                           + sum(nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,i_spin)* &
                           grdwts(:vla)* df_drhodgamma(:vla,4+i_spin)*sum(grarho(:vla,:3,3-i_spin)* &
                           nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1 )

                      dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea) &
                           + sum(nuc_grarho_imp(ii_ua)%m(:vla,kk,ii_ea,3-i_spin)* &
                           grdwts(:vla)* df_drhodgamma(:vla,7-i_spin)*sum(grarho(:vla,:3,3-i_spin)* &
                           nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin),2),1 )
                   enddo
                enddo
              ! end forall
            endif

            ! forall(k=1:3,kk=1:3) &    !!! t2_dervs [d/dg dfdgrarho]
            do kk = 1, 3
               do k = 1, 3
                  dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                       4.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,i_spin)* &
                       sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                       *grarho(:vla,:3,i_spin),2)* &
                       sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin) * &
                       grarho(:vla,:3,i_spin),2), 1)
               enddo
            enddo
            if (ispin .eq. 2) then
               ! forall(k=1:3,kk=1:3)
               do kk = 1, 3
                  do k = 1, 3
                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          4.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,ab)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin)*grarho(:vla,:3,3-i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          2.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,3+i_spin)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,3-i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin)*grarho(:vla,:3,i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          2.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,3+i_spin)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin)*grarho(:vla,:3,3-i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          2.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,6-i_spin)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,3-i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin)*grarho(:vla,:3,3-i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          2.0_r8_kind*sum( grdwts(:vla)* df_dgammadgamma(:vla,3+i_spin)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin)*grarho(:vla,:3,i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          sum( grdwts(:vla)* df_dgammadgamma(:vla,cc)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,3-i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,i_spin)*grarho(:vla,:3,3-i_spin),2), 1)

                     dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)=dervs_xc(i_ua,ii_ua)%m(k,i_ea,kk,ii_ea)+ &
                          sum( grdwts(:vla)* df_dgammadgamma(:vla,cc)* &
                          sum(nuc_sec_derrho(i_ua)%m(:vla,k,:3,i_ea,i_spin) &
                          *grarho(:vla,:3,3-i_spin),2)* &
                          sum( nuc_gragamma_imp(ii_ua)%m(:vla,kk,:3,ii_ea,3-i_spin)*grarho(:vla,:3,i_spin),2), 1)
                  enddo
               enddo
               ! end forall
            endif
             enddo
             endif !t2_dervs
         endif !else
    endif ! of if (explicit) ... else ...

             enddo dervs_ea
          enddo dervsw_ua
       enddo equals
    enddo grua
  end subroutine to_dervs_xc

  subroutine dervs_to_tmp()
    implicit none
    ! *** end of interface ***

    integer :: k, kk

    grua: do i_ua=1,n_unique_atoms
                i_ma = unique_atoms(i_ua)%moving_atom
                moving_atom = i_ma /= 0

             if( .not.moving_atom .and. .not.moving_grid) cycle

         equals: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

           dervsw_ua: do ii_ua=1,n_unique_atoms
             ii_ma = unique_atoms(ii_ua)%moving_atom
             moving_atom2 = ii_ma /= 0
            dervs_ea: do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms

              dervs_xc_expl: if(postscf_call ) then
               helpvar2=0.0_r8_kind
              do i_spin=1,ispin
              if (t1_dervs) then
                 ! (8)<-(1) - wdfdrho * nuc^lambda
                 do kk = 1, 3
                    do k = 1, 3
                       helpvar2(k, kk) = helpvar2(k, kk) - &
                            sum(grdwts(:vla) * dfdrho(:vla, i_spin) * &
                            nuc_dervsrho(i_ua, ii_ua)%m(:vla, k, i_ea, kk, ii_ea, i_spin))
                    enddo
                 enddo
              endif

               if (t2_dervs .and. nl_calc_ph) forall(k=1:3,kk=1:3)  &
                 helpvar2(k,kk)=helpvar2(k,kk)-two*sum(grdwts(:vla)*dfdgrarho(:vla,i_spin)* &
                  sum( grarho(:vla,:3,i_spin)*nuc_dervs_grarho(i_ua,ii_ua)%m(:vla,k,i_ea,kk,ii_ea,:3,i_spin) ,2))

              if(ispin .eq. 2 .and. t2_dervs .and. nl_calc_ph) forall(k=1:3,kk=1:3) &
                    helpvar2(k,kk)=helpvar2(k,kk)-sum(grdwts(:vla)*dfdgrarho(:vla,updn)*sum(&
                       grarho(:vla,:3,3-i_spin)*nuc_dervs_grarho(i_ua,ii_ua)%m(:vla,k,i_ea,kk,ii_ea,:3,i_spin),2) )
              enddo

                if(moving_atom2) &
                   dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) = &
                       dervs_xc_tmp(i_ma,ii_ma)%m(:,i_ea,:,ii_ea) + helpvar2

                if (weight_grads.and.moving_grid) then
                        dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) = &
                              dervs_xc_tmp(i_ma,i_m)%m(:,i_ea,:,1) - helpvar2
                endif

              endif dervs_xc_expl

             enddo dervs_ea
             enddo dervsw_ua
                enddo equals
             enddo grua
          end subroutine dervs_to_tmp

   subroutine tmp_to_dervs()

     do ii_ua=1,n_unique_atoms
      dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,:)=0.0_r8_kind
      dervs_xc_tmp(ii_ua,i_m)%m(:,:,:,1)=0.0_r8_kind
     enddo

     dervs_im_ua: do i_ua=1,n_unique_atoms
      i_ma = unique_atoms(i_ua)%moving_atom
         moving_atom = i_ma /= 0
          if( .not.moving_atom .and. .not.moving_grid) cycle

              do ii_ua=1,n_unique_atoms
               ii_ma = unique_atoms(ii_ua)%moving_atom
               if( ii_ma.eq.0 .and. .not.moving_grid) cycle

              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
              if(i_m.eq.i_ma.and.i_ea.eq.1) cycle
              do ii_ea=1,unique_atoms(ii_ua)%n_equal_atoms

               if(i_m.eq.ii_ma.and.ii_ea.eq.1) cycle

               dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,ii_ea)= &
                           dervs_xc_tmp(i_m,ii_ua)%m(:,1,:,ii_ea) &
                                           -dervs_xc_tmp(i_ua,ii_ua)%m(:,i_ea,:,ii_ea)
               dervs_xc_tmp(i_ua,i_m)%m(:,i_ea,:,1)= &
                           dervs_xc_tmp(i_ua,i_m)%m(:,i_ea,:,1) &
                                           -dervs_xc_tmp(i_ua,ii_ua)%m(:,i_ea,:,ii_ea)

              enddo
              enddo
              enddo
     enddo dervs_im_ua

     dervs_xc_tmp(i_m,i_m)%m(:,1,:,1)=0.0_r8_kind

       do i_ua=1,n_unique_atoms
        do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
          if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
               dervs_xc_tmp(i_m,i_m)%m(:,1,:,1) = &
               dervs_xc_tmp(i_m,i_m)%m(:,1,:,1)   &
                    -dervs_xc_tmp(i_m,i_ua)%m(:,1,:,i_ea)
        enddo
       enddo

      do i_ua=1,n_unique_atoms
      do ii_ua=1,n_unique_atoms
       dervs_xc(i_ua,ii_ua)%m = dervs_xc(i_ua,ii_ua)%m &
                              + dervs_xc_tmp(i_ua,ii_ua)%m
      enddo
      enddo
    end subroutine tmp_to_dervs
  end subroutine post_scf_calc_xc_en_and_gr_nl
#else
  subroutine post_scf_calc_xc_en_and_gr_nl()
    ! HISTORIC VERSION, FOR DEVELOPMENT USE THE ABOVE ONE!
    ! purpose : routine does the main calculations
    !           - calculation of prim. orbitals,gradients and sec. derivatives
    !           - calculation of becke weights
    !           - calculation of gradients of weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !           - calculation of GGA-gradients
    !** End of interface *****************************************
    use time_module
    use timer_module
    use ph_cntrl
    use xc_cntrl, only: xc_hcth_version
    use xc_func ! exc_*, xc_functionals
    use integralpar_module, only: integralpar_2dervs
    implicit none
    integer(i4_kind) :: i,ii,i_ea,i_ua,i_spin,k,i_vec,i_ma,i_m
    integer(i4_kind) :: vla,xyz
    logical :: split_gradients, moving_grid, moving_atom
    logical :: send_grad_grid, send_mda_coeff

    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths

    ! now some helparrays follow
    real(kind=r8_kind) :: helpvar(3),fxc_help(vec_length),&
         & grarho(vec_length,X:Z,ispin),&
         & help_vec(vec_length),&
         & help_xyz_ud(vec_length,X:Z,UP:DN)
    real(kind=r8_kind) :: gamma(vec_length,1+2*(ispin-1)),dfdgrarho(vec_length,1+2*(ispin-1))
    ! now pointer on gradients of weights, gradients of the energies and gradients and 2nd
    ! gradients of the density follow
    real(kind=r8_kind),pointer :: gr(:),nuc(:,:,:),sec_der(:,:,:,:)
!!$    real(kind=r8_kind),pointer :: nuc_core(:,:,:), sec_der_core(:,:,:,:)

    split_gradients = options_split_gradients()

    ! initialize timers
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)

    ! set gradients to zero at the beginning
    do i=1,n_unique_atoms
       grad_xc(i)%m=0.0_r8_kind
       if (weight_grads .and. split_gradients) then
          grad_grid(i)%m=0.0_r8_kind
       end if
    end do
    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop

    FPP_TIMER_START(loop)
    FPP_TIMER_START(sched)
    do while (more_grid_atom(vec_length, i, grdpts, grdwts)) ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
       FPP_TIMER_STOP(sched)
       FPP_TIMER_START(work)
       i_m = unique_atoms(i)%moving_atom
       moving_grid = i_m /= 0 .and. weight_grads
       do ii=1,n_unique_atoms
          nuc_grarho(ii)%m=0.0_r8_kind
          nuc_sec_derrho(ii)%m=0.0_r8_kind
       end do
       rho=0.0_r8_kind
       dfdrho=0.0_r8_kind

       fxc=0.0_r8_kind
       fxc_help=0.0_r8_kind
       gamma=0.0_r8_kind
       grarho = ZERO
       dfdgrarho=0.0_r8_kind
       vla=size(grdpts,1)

       call start_timer(timer_gridph_grdwts)
       if(weight_grads) then
          do ii=1,n_unique_atoms
             graw(ii)%m=0.0_r8_kind
          end do
          help_vec(1:vla)=atomicweight_and_grad(i,grdpts(1:vla,:),graw)
          grdwts(1:vla)=grdwts(1:vla)*&
               unique_atoms(i)%n_equal_atoms
          do i_ua=1,n_unique_atoms
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                graw(i_ua)%m(1:vla,1,i_ea)=&
                     graw(i_ua)%m(1:vla,1,i_ea)*&
                     grdwts(1:vla)
                graw(i_ua)%m(1:vla,2,i_ea)=&
                     graw(i_ua)%m(1:vla,2,i_ea)*&
                     grdwts(1:vla)
                graw(i_ua)%m(1:vla,3,i_ea)=&
                     graw(i_ua)%m(1:vla,3,i_ea)*&
                     grdwts(1:vla)
             end do
          end do
          grdwts(1:vla)=grdwts(1:vla)*&
               help_vec(1:vla)
       else
          grdwts(1:vla)=grdwts(1:vla)*&
               atomicweight(i,grdpts(1:vla,:))*&
               unique_atoms(i)%n_equal_atoms
       end if
       call stop_timer(timer_gridph_grdwts)
       call start_timer(timer_gridph_orbitals)
       call orbital_calculate(grdpts(1:vla,1:3),&
            vla,orbs_ob=orbs_ob,grads=orbs_grads,&
            nuc_grads=nuc_grads,nuc_sec_der=nuc_sec_der)
       call stop_timer(timer_gridph_orbitals)
       call start_timer(timer_gridph_density)
       call density_calc_ph_nl(vla,rho,gamma,grarho &
            , nuc_grarho=nuc_grarho &
            , nuc_sec_derrho=nuc_sec_derrho &
            , orbs_ob=orbs_ob &
            , orbs_grads=orbs_grads &
            , nuc_grads=nuc_grads &
            , nuc_sec_der=nuc_sec_der )

       ! performing charge and exc integration over the grid
       charge_int(i,1)=charge_int(i,1)+sum(rho(1:vla,:)*spread&
            (grdwts(1:vla),2,ispin))

       call stop_timer(timer_gridph_density)
       call start_timer(timer_gridph_functionals)

       call xc_functionals(vla,ispin,rho,fxc,dfdrho,gamma,dfdgrarho,GRDWTS=grdwts)

       call stop_timer(timer_gridph_functionals)


       ! now start building the gradient
       if(ispin==1) then
          dfdgrarho(1:vla,1)=dfdgrarho(1:vla,1)*&
               grdwts(1:vla)
          dfdrho(1:vla,1)=dfdrho(1:vla,1)*&
               grdwts(1:vla)

          do xyz=X,Z
             help_xyz_ud(:vla,xyz,UP) = dfdgrarho(1:vla,UP)*2.0_r8_kind*&
                  & grarho(1:vla,xyz,UP)
          enddo

          grad: do i_ua=1,n_unique_atoms
             i_ma = unique_atoms(i_ua)%moving_atom
             moving_atom = i_ma /= 0
          if( .not.moving_atom .and. .not.moving_grid) cycle
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                helpvar=0.0_r8_kind
                nuc=>nuc_grarho(i_ua)%m(:,:,i_ea,:)
                sec_der=>nuc_sec_derrho(i_ua)%m(:,:,:,i_ea,:)
                do k=1,3
                   do i_vec=1,vla
                      helpvar(k)=helpvar(k)+dfdrho(i_vec,1)*&
                           nuc(i_vec,k,1)
                   end do
                end do

                ! dfdgrarho
                do xyz=X,Z
                   do i_vec=1,vla
                      helpvar(xyz)=helpvar(xyz)&
                           +help_xyz_ud(i_vec,X,1)*sec_der(i_vec,xyz,X,1)&
                           +help_xyz_ud(i_vec,Y,1)*sec_der(i_vec,xyz,Y,1)&
                           +help_xyz_ud(i_vec,Z,1)*sec_der(i_vec,xyz,Z,1)
                   enddo
                enddo

                if( moving_atom )then
                  gr => grad_xc(i_ma)%m(:,i_ea)
                  gr =  gr + helpvar
                endif
                if (weight_grads) then
                   if (split_gradients) then
                      if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                      if( moving_grid ) grad_grid(i_m)%m(:,1) = &
                                        grad_grid(i_m)%m(:,1) - helpvar
                   else
                      if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                                        grad_xc(i_m)%m(:,1) - helpvar
                   endif
                   helpvar=0.0_r8_kind
                   do k=1,3
                      do i_vec=1,vla
                         helpvar(k)=helpvar(k)-fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                      end do
                   end do
                   if( moving_atom ) gr = gr + helpvar
                end if
             enddo
          enddo grad
       else          ! spinpolarized case
!-----------------------------------------------------------------------+
! Calculation of first part of rho terms:                               |
!     UP: df_drhoa * w                                                  |
!     DN: df_drhob * w                                                  |
!-----------------------------------------------------------------------+
          dfdrho(1:vla,1)=dfdrho(1:vla,1)*&
               grdwts(1:vla)
          dfdrho(1:vla,2)=dfdrho(1:vla,2)*&
               grdwts(1:vla)
!-----------------------------------------------------------------------+
! Calculation of first part of gamma terms:                             |
!     UP:                                                               |
!      (2 df_dgammaaa * drhoa_dr + df_dgammaab * drhob_dr) * w          |
!     DN:                                                               |
!      (2 df_dgammabb * drhob_dr + df_dgammaab * drhoa_dr) * w          |
!-----------------------------------------------------------------------+
          do xyz=X,Z
             help_xyz_ud(1:vla,xyz,UP) = (2.0_r8_kind*&
                  & dfdgrarho(1:vla,1)*grarho(1:vla,xyz,1)+&
                  & dfdgrarho(1:vla,3)*grarho(1:vla,xyz,2))*&
                  & grdwts(1:vla)
             help_xyz_ud(1:vla,xyz,DN) = (2.0_r8_kind*&
                  & dfdgrarho(1:vla,2)*grarho(1:vla,xyz,2)+&
                  & dfdgrarho(1:vla,3)*grarho(1:vla,xyz,1))*&
                  & grdwts(1:vla)
          enddo
          do i_ua=1,n_unique_atoms
             i_ma = unique_atoms(i_ua)%moving_atom
             moving_atom = i_ma /= 0
          if( .not.moving_atom .and. .not.moving_grid) cycle
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                helpvar=0.0_r8_kind
                nuc=>nuc_grarho(i_ua)%m(:,:,i_ea,:)
                sec_der=>nuc_sec_derrho(i_ua)%m(:,:,:,i_ea,:)

!-----------------------------------------------------------------------+
! Complete calculation of rho terms:                                    |
!                                                                       |
! gr = (df_drhoa * nuca + df_drhob * nucb) * w                          |
!-----------------------------------------------------------------------+
                do i_spin=1,ispin
                   do k=1,3
                      do i_vec=1,vla
                         helpvar(k)=helpvar(k)+dfdrho(i_vec,i_spin)*&
                              nuc(i_vec,k,i_spin)
                      end do
                   end do
                end do

!-----------------------------------------------------------------------+
! Complete calculation of gamma terms to gradients                      |
!                                                                       |
! gr = ((2 df_dgammaaa * drhoa_dr + df_dgammaab * drhob_dr) * sec_der_a |
!   +(2 df_dgammabb * drhob_dr + df_dgammaab * drhoa_dr) * sec_der_b)*w |
!                                                                       |
! with:         sec_der_a = [-d/dR]^T drhoa_dr                          |
!               sec_der_b = [-d/dR]^T drhob_dr                          |
!                                                                       |
! note: sec_der=>nuc_sec_derrho                                         |
!-----------------------------------------------------------------------+
                do xyz=X,Z
                   do i_vec=1,vla
                      helpvar(xyz)=helpvar(xyz)&
                           +help_xyz_ud(i_vec,X,UP)*sec_der(i_vec,xyz,X,UP)&
                           +help_xyz_ud(i_vec,Y,UP)*sec_der(i_vec,xyz,Y,UP)&
                           +help_xyz_ud(i_vec,Z,UP)*sec_der(i_vec,xyz,Z,UP)
                      helpvar(xyz)=helpvar(xyz)&
                           +help_xyz_ud(i_vec,X,DN)*sec_der(i_vec,xyz,X,DN)&
                           +help_xyz_ud(i_vec,Y,DN)*sec_der(i_vec,xyz,Y,DN)&
                           +help_xyz_ud(i_vec,Z,DN)*sec_der(i_vec,xyz,Z,DN)
                   enddo
                enddo

                if( moving_atom )then
                   gr => grad_xc(i_ma)%m(:,i_ea)
                   gr =  gr + helpvar
                endif
                if (weight_grads) then
                   if (split_gradients) then
                      if( moving_atom ) gr => grad_grid(i_ma)%m(:,i_ea)
                      if( moving_grid ) grad_grid(i_m)%m(:,1) = &
                                        grad_grid(i_m)%m(:,1) - helpvar
                   else
                      if( moving_grid ) grad_xc(i_m)%m(:,1) = &
                                        grad_xc(i_m)%m(:,1) - helpvar
                   endif
                   helpvar=0.0_r8_kind
                   do k=1,3
                      do i_vec=1,vla
                         helpvar(k)=helpvar(k)-fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                      end do
                   end do
                   if( moving_atom ) gr = gr + helpvar
                end if
             enddo
          end do
       end if
       FPP_TIMER_STOP(work)
       FPP_TIMER_START(sched)
    end do
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)
    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(.true., send_grad_grid, send_mda_coeff)
  end subroutine post_scf_calc_xc_en_and_gr_nl
#endif


!AG [ Insert from pg_V2.1_MDA (subroutine post_scf_calc_xcmda_etot_grads)
!gamm => gamma IS NOT CHANGED

  subroutine post_scf_calc_xcmda_etot_grads
    ! purpose : routine does the main calculations
    !           - calculation of fit functions
    !           - calculation of becke weights
    !           - calculation of density
    !           - evaluation of functionals
    !           - integration of charge and energy
    !           - evaluation of gradients [if do_grads = .TRUE.]
    !** End of interface *****************************************
    use mat_charge_module, only: matinv_charge, dual_charge_norm
    use linsys_module, only: solve_linsys
    use time_module
    use timer_module
!temp comment    use rspace_plot
    use gamma_module
    use ph_cntrl
    use iounitadmin_module
    use xc_func ! exc_*, xc_functionals
    implicit none
    real(kind=r8_kind), allocatable :: grad_rho(:,:,:), sec_der_rho(:,:,:), &
                                       gamm(:,:), dfdgrarho(:,:), &
                                       tmp(:), theta(:,:), delta(:,:), hlp(:,:)
    target :: grad_rho, sec_der_rho
    real(kind=r8_kind), pointer :: nuc(:,:), sec_der(:,:,:), gr(:), &
                                   f(:,:), g(:,:,:)
    integer(kind=i4_kind) :: i,vl,s,t,i_m,i_ua,i_ea,i_ma,k,stat,i_vec,i_f
    real(kind=r8_kind)    :: a2,reg,zero_target
    real(kind=r8_kind), parameter :: zero = 0.0_r8_kind, half = 0.5_r8_kind, &
                                     two  = 2.0_r8_kind
    real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:) ! gridpoints and weigths
    logical               :: split_gradients, moving_grid
    logical :: send_grad_grid, send_mda_coeff
    ! prepare some timers!
    call init_timer(timer_gridph_orbitals)
    call init_timer(timer_gridph_density)
    call init_timer(timer_gridph_functionals)
    call init_timer(timer_gridph_grdwts)
    call init_timer(timer_gridph_gradients)

    ! allocation and initialization
    a2 = rho_shape_eps*rho_shape_eps
    if (nl_calc_ph) then ! GGA calculation
       allocate( theta    (vec_length,  ispin  ), &
                 grad_rho (vec_length,3,ispin  ), &
                 gamm    (vec_length,2*ispin-1), &
                 dfdgrarho(vec_length,2*ispin-1), stat = stat)
    ASSERT(stat.eq.0)

       if (do_grads.or.comp_exact .or. a2 /= zero ) then
          allocate( tmp(vec_length), stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: alloc of tmp failed')
       endif
       if (do_grads) then
          allocate( hlp(vec_length,3), stat = stat )
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: alloc of hlp failed')
          if (weight_grads) then
             allocate( sec_der_rho(vec_length,6,ispin), stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: alloc of sec_der_rho failed')
          endif
          if (a2 /= zero) then
             allocate( delta(vec_length,ispin), stat = stat)
          ASSERT(stat.eq.0)
          endif
       elseif(comp_exact) then
          allocate( hlp(vec_length,3), stat = stat )
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: alloc of hlp failed')
          if (a2 /= zero) then
             allocate( delta(vec_length,ispin), stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: alloc of delta failed')
          endif
       endif

    else ! LDA calculation
       if (do_grads) then
          allocate( theta(vec_length,ispin), &
                    tmp  (vec_length      ), stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: alloc of theta and tmp failed')
          if (weight_grads) then
             allocate( grad_rho(vec_length,3,ispin), stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: alloc of grad_rho failed')
          endif
       elseif (comp_exact) then
          allocate( theta(vec_length,ispin), &
                    tmp  (vec_length      ), stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: alloc of theta and tmp failed')
       endif
    endif
    if (do_grads) then
       split_gradients = options_split_gradients()
       do i=1,N_moving_unique_atoms
          grad_xc(i)%m=0.0_r8_kind
          if (weight_grads .and. split_gradients) then
             grad_grid(i)%m=0.0_r8_kind
          end if
       end do
       rhs = 0.0_r8_kind
    endif
    if(comp_exact) rhs=0.0_r8_kind

    ! main processing loop
    call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
    FPP_TIMER_START(loop)
    FPP_TIMER_START(sched)
    do while (more_grid_atom(vec_length, i, grdpts, grdwts)) ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
       FPP_TIMER_STOP(sched)
       FPP_TIMER_START(work)
       i_m = unique_atoms(i)%moving_atom
       moving_grid = do_grads .and. weight_grads .and. i_m > 0
       rho    = zero
       dfdrho = zero
       fxc    = zero
       if (do_grads .or. nl_calc_ph) then
          grad_rho  = zero
       endif
       if (nl_calc_ph) then
          gamm     = zero
          dfdgrarho = zero
          if (moving_grid .and. .not. comp_exact) then
             sec_der_rho = zero
          endif
       endif
       if (do_grads) then
          do i_ua=1,n_unique_atoms
             nuc_grarho(i_ua)%m = zero
          end do
          if (nl_calc_ph) then
             do i_ua=1,n_unique_atoms
                nuc_sec_derrho(i_ua)%m = zero
             end do
          endif
       endif
       ! fetching part of the grid
       vec_length_act=size(grdpts,1)
       vl = vec_length_act

       call start_timer(timer_gridph_grdwts)
       if (do_grads .and. weight_grads) then
          do i_ua=1,n_unique_atoms
             graw(i_ua)%m=0.0_r8_kind
          end do
          tmp   (1:vl) = atomicweight_and_grad(i,grdpts(1:vl,:),graw)
          grdwts(1:vl) = grdwts(1:vl) * unique_atoms(i)%n_equal_atoms
          do i_ua=1,n_unique_atoms
             do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                graw(i_ua)%m(1:vl,1,i_ea) = graw(i_ua)%m(1:vl,1,i_ea) * grdwts(1:vl)
                graw(i_ua)%m(1:vl,2,i_ea) = graw(i_ua)%m(1:vl,2,i_ea) * grdwts(1:vl)
                graw(i_ua)%m(1:vl,3,i_ea) = graw(i_ua)%m(1:vl,3,i_ea) * grdwts(1:vl)
             end do
          end do
          grdwts(1:vl) = grdwts(1:vl) * tmp(1:vl)
       else ! gradients of weights and nodes not required
          grdwts(1:vl) = grdwts(1:vl) * atomicweight(i,grdpts(1:vl,:)) &
                                      * unique_atoms(i)%n_equal_atoms
       end if
       call stop_timer(timer_gridph_grdwts)

       call start_timer(timer_gridph_orbitals)
       if (nl_calc_ph) then ! GGA calculation
          if (do_grads) then ! GGA gradients
             if (moving_grid) then ! d^2/dr^2 rho required
                call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch", &
                     fcts=fcts_ch,grads=grads_ch,sec_ders=sec_ders_ch,&
                     nuc_grads=nuc_grads_ch,sec_nuc_ders=sec_nuc_ders_ch)
             else ! d^2/dr^2 rho not required
                call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch", &
                     fcts=fcts_ch,grads=grads_ch,nuc_grads=nuc_grads_ch,&
                     sec_nuc_ders=sec_nuc_ders_ch)
             endif
          elseif(comp_exact) then
             call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch", &
                     fcts=fcts_ch,grads=grads_ch)
          else ! GGA energy only
             call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch", &
                  fcts=fcts_ch,grads=grads_ch)
          endif
       else ! LDA calculation
          if (do_grads) then ! LDA gradients
             call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch", &
                  fcts=fcts_ch,grads=grads_ch,nuc_grads=nuc_grads_ch)
          elseif(comp_exact) then
             call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch",fcts=fcts_ch)
          else ! LDA energy only
             call fit_fct_calculate(grdpts(1:vl,1:3),vl,"ch",fcts=fcts_ch)
          endif
       endif
       call stop_timer(timer_gridph_orbitals)

       call start_timer(timer_gridph_density)
       if (nl_calc_ph) then ! GGA calculation
          if (do_grads) then ! GGA gradients
             if (moving_grid) then ! d^2/dr^2 rho required
                call fitted_density_calc(vl,rho=rho,grad=grad_rho,&
                     sec_der=sec_der_rho,nuc_grad=nuc_grarho, &
                     sec_nuc_der=nuc_sec_derrho,&
                     fcts=fcts_ch,grads=grads_ch,sec_ders=sec_ders_ch,&
                     nuc_grads=nuc_grads_ch,sec_nuc_ders=sec_nuc_ders_ch)
             elseif (comp_exact) then
                call fitted_density_calc(vl,rho=rho,grad=grad_rho,&
                             fcts=fcts_ch,grads=grads_ch)
             else ! d^2/dr^2 rho not required
                call fitted_density_calc(vl,rho=rho,grad=grad_rho,&
                     nuc_grad=nuc_grarho,sec_nuc_der=nuc_sec_derrho,&
                     fcts=fcts_ch,grads=grads_ch,nuc_grads=nuc_grads_ch,&
                     sec_nuc_ders=sec_nuc_ders_ch)
             endif
          else ! GGA energy only
             call fitted_density_calc(vl,rho=rho,grad=grad_rho,&
                  fcts=fcts_ch,grads=grads_ch)
          endif
       else ! LDA calculation
          if (do_grads) then ! LDA gradients
             call fitted_density_calc(vl,rho=rho,grad=grad_rho,&
                  nuc_grad=nuc_grarho,fcts=fcts_ch,grads=grads_ch,&
                  nuc_grads=nuc_grads_ch)
          elseif (comp_exact) then
             call fitted_density_calc(vl,rho=rho,fcts=fcts_ch)
          else ! LDA energy only
             call fitted_density_calc(vl,rho=rho,fcts=fcts_ch)
          endif
       endif
       ! truncate the fitted density
!..............................................................................
! density truncation
! ~~~~~~~~~~~~~~~~~~
! rho(trunc)_s(r) := S( rho_s(r) )  >= 0
!
! soft truncation                             hard trunctation
! ~~~~~~~~~~~~~~~                             ~~~~~~~~~~~~~~~~
! S (x) = [ x + sqrt( a^2 + x^2 ) ] / 2       S (x) = max(0,x)
! S'(x) =    S(x) / sqrt( a^2 + x^2 )         S'(x) = Theta(x)
! S"(x) = 1/2 a^2 / sqrt( a^2 + x^2 )^3       S"(x) = delta(x) = 0 for x <> 0
!
! Input : rho_s(r)
! Output: S(rho_(r)), S'(rho_s(r)), S"(rho_s(r))
!..............................................................................
       do s=1,ispin
          t = s+ispin
          charge_int(i,s) = charge_int(i,s) + sum( rho(:vl,s)*grdwts(:vl) )
          if (nl_calc_ph .or. do_grads.or.comp_exact) then ! GGA calculations or gradients
             if (a2 == zero) then ! hard truncation
                theta(:vl,s) = half + sign(half,rho(:vl,s)) ! Theta(rho)
                rho  (:vl,s) = theta(:vl,s) * rho(:vl,s)    ! max(0,rho)
             else
                tmp(:vl  ) = sqrt( a2 + rho(:vl,s)*rho(:vl,s) )
                rho(:vl,s) = half*( rho(:vl,s) + tmp(:vl) )
                if ((do_grads.or.comp_exact) .and. nl_calc_ph) then
                   delta(:vl,s) = half * a2 / ( tmp(:vl)*tmp(:vl)*tmp(:vl) )
                endif
                theta(:vl,s) = rho(:vl,s) / tmp(:vl)
             endif
          else ! LDA energy only
             if (a2 == zero) then
                rho(:vl,s) = max( rho(:vl,s), zero )
             else
                rho(:vl,s) = half*( rho(:vl,s) + &
                                    sqrt( a2 + rho(:vl,s)*rho(:vl,s) ) )
             endif
          endif
          charge_int(i,t) = charge_int(i,t) + sum( rho(:vl,s)*grdwts(:vl) )
       enddo
       ! now build gamm of the truncated fitted density
!..............................................................................
! gamm_st(r) = < Nabla rho(trunc)_s(r) | Nabla rho(trunc)_t(r) >
! ==>
! gamm_st(r) = S'(rho_s(r)) < Nabla rho_s(r) | Nabla rho_t(r) > S'(rho_t(r))
!..............................................................................
       if (nl_calc_ph) then
          gamm(:vl,1) = grad_rho(:vl,1,1)*grad_rho(:vl,1,1) &
                       + grad_rho(:vl,2,1)*grad_rho(:vl,2,1) &
                       + grad_rho(:vl,3,1)*grad_rho(:vl,3,1)
          gamm(:vl,1) = theta(:vl,1) * gamm(:vl,1) * theta(:vl,1)
          if (ispin > 1) then
             gamm(:vl,2) = grad_rho(:vl,1,2)*grad_rho(:vl,1,2) &
                          + grad_rho(:vl,2,2)*grad_rho(:vl,2,2) &
                          + grad_rho(:vl,3,2)*grad_rho(:vl,3,2)
             gamm(:vl,2) = theta(:vl,2) * gamm(:vl,2) * theta(:vl,2)
             gamm(:vl,3) = grad_rho(:vl,1,1)*grad_rho(:vl,1,2) &
                          + grad_rho(:vl,2,1)*grad_rho(:vl,2,2) &
                          + grad_rho(:vl,3,1)*grad_rho(:vl,3,2)
             gamm(:vl,3) = theta(:vl,1) * gamm(:vl,3) * theta(:vl,2)
          endif
       endif
       call stop_timer(timer_gridph_density)

       ! evaluate the exchange correlation functionals
!..............................................................................
! E_xc[rho_t] = Int(R^3) e_xc[rho_t](r) dV
! with
! e_xc[rho_t](r) = e_xc(rho_t(r),Nabla rho_t'(r))
! ==>
! < d/d rho_s E_xc[rho_t] | psi > =
! Int(R^3) V_xc,s[rho_t](r) psi(r) +  W_xc,s[rho_t](r) d/dr psi(r) dV
! with
! V_xc,s[rho_t](r) = [d/d       rho_s e_xc](rho_t(r),Nabla rho_t'(r))  and
! W_xc,s[rho_t](r) = [d/d Nabla rho_s e_xc](rho_t(r),Nabla rho_t'(r))
!
! e_xc(rho_r,Nabla rho_t') = f_xc(rho_t,gamm_tt')
! with
! gamm_tt'(r) = < Nabla rho_t(r) | Nabla rho_t'(r) > ; tt' = 11,22,12
! ==>
! V_xc,s[rho_t](r) = [d/d rho_s f_xc](rho_t(r),gamm_tt'(r))  and
! W_xc,s[rho_t](r) = 2 G_xc,s s[rho_t](r) d/dr rho_ s(r)
!                  +   G_xc,s-s[rho_t](r) d/dr rho_-s(r)
! with
! G_xc,ss'[rho_t](r) := [d/d gamm_ss' f_xc](rho_t(r),gamm_tt'(r))
!
! RKS : no spin dependencies _t, _tt', _s, and _ss'
!       no mixed gamm derivative contribution G_xc,s-s[rho_t](r)
! LDA : no Nabla rho_t(r) and gamm_tt'(r) dependencies
!       no W_xc,s[rho_t] or W_xc[rho] potential contributions
!
! Input : S(rho_t(r)), gamm_tt'(r)
! Output: e_xc[S(rho_t)](r), V_xc,s[S(rho_t)](r), G_xc,ss'[S(rho_t)](r)
!..............................................................................
       call start_timer(timer_gridph_functionals)

       if(do_grads.or.comp_exact)then
          ! need Fxc,dFdn (and dFdg for gga)
          if(nl_calc_ph)then
             call xc_functionals(vl,ispin,rho,fxc,dfdrho,gamm,dfdgrarho,GRDWTS=grdwts)
          else
             call xc_functionals(vl,ispin,rho,fxc,dfdrho,GRDWTS=grdwts)
          endif
       else
          ! dont need need Fxc, only energies
          if(nl_calc_ph)then
             call xc_functionals(vl,ispin,rho,GAM=gamm,GRDWTS=grdwts)
          else
             call xc_functionals(vl,ispin,rho,GRDWTS=grdwts)
          endif
       endif

       call stop_timer(timer_gridph_functionals)

       ! start gradient evaluation
       call start_timer(timer_gridph_gradients)
       if (do_grads.or.comp_exact) then
          do s=1,ispin

             ! first load the weighted exchange-correlation potentials
             ! w_a V_xc,s[S(rho_t)](r_a) and w_a W_xc,s[S(rho_t)](r_a)
!..............................................................................
! W_xc,s[S(rho_t)](r) = 2 G_xc,s s[S(rho_t)](r) S'(rho_ s) d/dr rho_ s(r)
!                     +   G_xc,s-s[S(rho_t)](r) S'(rho_-s) d/dr rho_-s(r)
!..............................................................................
             if (nl_calc_ph) then ! GGA co-potential
                ! 2 w_a G_xc,ss[rho_t](r_a) d/dr rho_s(r_a)
                tmp(1:vl) = two*dfdgrarho(1:vl,s)*grdwts(1:vl)*theta(1:vl,s)
                do k=1,3
                   hlp(1:vl,k) = tmp(1:vl) * grad_rho(1:vl,k,s)
                end do
                if (ispin > 1) then
                   ! w_a G_xc,s-s[rho_t](r_a) d/dr rho_-s(r_a)
                   tmp(1:vl) = dfdgrarho(1:vl,3)*grdwts(1:vl)*theta(1:vl,3-s)
                   do k=1,3
                      hlp(1:vl,k) = hlp(1:vl,k) + &
                                    tmp(1:vl) * grad_rho(1:vl,k,3-s)
                   end do
                end if
             end if
             ! w_a V_xc,s[rho_t](r_a)
             tmp(1:vl) = dfdrho(1:vl,s) * grdwts(1:vl)

             ! then perform the modifications due to the truncation
!..............................................................................
! E(trunc)_xc[rho_t] := E_xc[rho(trunc)_t] = E_xc[S(rho_t)]
! ==>
! e(trunc)_xc[rho_t](r) = e_xc(S(rho_t(r)),S'(rho_t'(r))*Nabla rho_t'(r))
! ==>
! V(trunc)_xc,s[rho_t](r) = V_xc,s[S(rho_t)](r) * S'(rho_s(r)) + ...
!                           W_xc,s[S(rho_t)](r) * S"(rho_s(r)) * Nabla rho_s(r)
! and
! W(trunc)_xc,s[rho_t](r) = W_xc,s[S(rho_t)](r) * S'(rho_s(r))
!
! RKS : no spin dependencies _t and _s
! LDA : no W_xc,s[rho_t] or W_xc[rho] contributions
!..............................................................................
             tmp(1:vl) = tmp(1:vl) * theta(:vl,s)
             if (nl_calc_ph) then ! GGA co-potential
                if (a2 /= zero) then ! soft truncation
                   tmp(1:vl) = tmp(1:vl) + &
                        hlp(1:vl,1) * delta(:vl,s) * grad_rho(1:vl,1,s) + &
                        hlp(1:vl,2) * delta(:vl,s) * grad_rho(1:vl,2,s) + &
                        hlp(1:vl,3) * delta(:vl,s) * grad_rho(1:vl,3,s)
                end if
                do k=1,3
                   hlp(1:vl,k) = hlp(1:vl,k) * theta(:vl,s)
                end do
             end if

             if(do_grads) then
             ! finally sum up the individual gradient contributions
!..............................................................................
! -d/dR E_xc[rho_t] =
!       + Sum(a,s) w_a V_xc,s[rho_t](r_a) [-d/dR      rho_s](r_a)
!       + Sum(a,s) w_a W_xc,s[rho_t](r_a) [-d/dR d/dr rho_s](r_a)
!       - Sum(a,s) w_a V_xc,s[rho_t](r_a) [ d/dr      rho_s](r_a) (d/dR r_a)
!       - Sum(a,s) w_a W_xc,s[rho_t](r_a) [ d/dr d/dr rho_s](r_a) (d/dR r_a)
!       - Sum(a) [d/dR w_a] e_xc[rho_t](r_a)
! with
!       d/dR r_a = delta_R,R0  .
!
! RKS : no sum over s and no spin dependencies _t and _s
! LDA : no W_xc,s[rho_t] or W_xc[rho] contributions
!..............................................................................
             do i_ua=1,n_unique_atoms
                i_ma = unique_atoms(i_ua)%moving_atom
                if (i_ma > 0) then
                   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                      ! - Sum(a,s) w_a V_xc,s[rho_t](r_a) [d/dR rho_s](r_a)
                      gr => grad_xc(i_ma)%m(:,i_ea)
                      nuc => nuc_grarho(i_ua)%m(:,:,i_ea,s)
                      do k=1,3
                         reg = zero
                         do i_vec=1,vec_length_act
                            reg = reg + tmp(i_vec)*nuc(i_vec,k)
                         end do
                         gr(k) = gr(k) + reg
                      end do
                      if (nl_calc_ph) then
                         ! - Sum(a,s) w_a W_xc,s[rho_t](r_a) ...
                         !            ... [d/dR d/dr rho_s](r_a)
                         sec_der => nuc_sec_derrho(i_ua)%m(:,:,:,i_ea,s)
                         do k=1,3
                            reg = zero
                            do i_vec=1,vec_length_act
                               reg = reg + hlp(i_vec,1)*sec_der(i_vec,k,1) &
                                         + hlp(i_vec,2)*sec_der(i_vec,k,2) &
                                         + hlp(i_vec,3)*sec_der(i_vec,k,3)
                            end do
                            gr(k) = gr(k) + reg
                         end do
                      end if
                      if (weight_grads .and. s == 1) then
                         ! - Sum(a) [d/dR w_a] e_xc[rho_t](r_a) >>
                         if (split_gradients) then
                            gr => grad_grid(i_ma)%m(:,i_ea)
                         end if
                         do k=1,3
                            reg = zero
                            do i_vec=1,vec_length_act
                               reg = reg + fxc(i_vec)*graw(i_ua)%m(i_vec,k,i_ea)
                            end do
                            gr(k) = gr(k) - reg
                         end do
                      end if
                   end do ! i_ea
                end if ! i_ma > 0
             end do ! i_ua

             if (moving_grid) then
                ! - Sum(a,s) w_a V_xc,s[rho_t](r_a) [d/dr rho_s](r_a) >>
                if (split_gradients) then
                   gr => grad_grid(i_m)%m(:,1)
                else
                   gr => grad_xc(i_m)%m(:,1)
                end if
                do k=1,3
                   reg = zero
                   do i_vec=1,vec_length_act
                      reg = reg + tmp(i_vec)*grad_rho(i_vec,k,s)
                   end do
                   gr(k) = gr(k) - reg
                end do
                if (nl_calc_ph) then
                   ! - Sum(a,s) w_a W_xc,s[rho_t](r_a) ...
                   !            ... [d/dr d/dr rho_s](r_a) >>

                   ! FIXME: find a better way to multiply
                   ! a vector by the matrix stored in a
                   ! triangular mode:
                   ! k=1; sum(J) XJ >>>
                   reg = zero
                   do i_vec=1,vec_length_act
                      reg = reg &
                           + hlp(i_vec,1)*sec_der_rho(i_vec,XX,s) &
                           + hlp(i_vec,2)*sec_der_rho(i_vec,XY,s) &
                           + hlp(i_vec,3)*sec_der_rho(i_vec,XZ,s)
                   end do
                   gr(1) = gr(1) - reg
                   ! k=2; sum(J) YJ >>>
                   reg = zero
                   do i_vec=1,vec_length_act
                      reg = reg &
                           + hlp(i_vec,1)*sec_der_rho(i_vec,XY,s) & ! YX
                           + hlp(i_vec,2)*sec_der_rho(i_vec,YY,s) &
                           + hlp(i_vec,3)*sec_der_rho(i_vec,YZ,s)
                   end do
                   gr(2) = gr(2) - reg
                   ! k=3; sum(J) ZJ >>>
                   reg = zero
                   do i_vec=1,vec_length_act
                      reg = reg &
                           + hlp(i_vec,1)*sec_der_rho(i_vec,XZ,s) & ! ZX
                           + hlp(i_vec,2)*sec_der_rho(i_vec,YZ,s) & ! ZY
                           + hlp(i_vec,3)*sec_der_rho(i_vec,ZZ,s)
                   end do
                   gr(3) = gr(3) - reg
                end if
             end if
             endif !do_grads in do_grads.or.comp_exact

                ! and perform the potential projections < f_k | V_xc,s >num
!..............................................................................
! < d/d rho_s E_xc[rho_t] | f_k > =
! Int(R^3) V_xc,s[rho_t](r) f_k(r) +  W_xc,s[rho_t](r) d/dr f_k(r) dV
!..............................................................................
             f => fcts_ch%o(:,:,1)
             if (nl_calc_ph) then ! GGA potential
                g => grads_ch%o(:,:,:,1)
                do i_f=1,n_ch
!!! MF Bug?! vgl. xcmda_hamiltonian
!                     rhs(i_f,s) = sum( tmp(:vl  ) * f(:vl,  i_f) + &
                   rhs(i_f,s) = rhs(i_f,s) + sum( tmp(:vl  ) * f(:vl,  i_f) + &
                                     hlp(:vl,1) * g(:vl,1,i_f) + &
                                     hlp(:vl,2) * g(:vl,2,i_f) + &
                                     hlp(:vl,3) * g(:vl,3,i_f) )
                end do
             else ! LDA potential
                do i_f=1,n_ch
!!! MF Bug?! vgl. xcmda_hamiltonian
!                 rhs(i_f,s) = sum( tmp(:vl) * f(:vl,i_f) )
                  rhs(i_f,s) = rhs(i_f,s) + sum( tmp(:vl) * f(:vl,i_f) )
               end do
            endif

         end do ! s
      end if ! do_grads
      call stop_timer(timer_gridph_gradients)

      FPP_TIMER_STOP(work)
      FPP_TIMER_START(sched)
    end do! loop over gridpoints
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
    FPP_TIMER_STOP(sched)
    FPP_TIMER_STOP(loop)

    if (nl_calc_ph) then ! GGA calculation
       deallocate( theta, grad_rho, gamm, dfdgrarho, stat = stat)
       if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
            &_grads: dealloc of theta, grad_rho, gamm and dfdgrarho failed')

       if (do_grads .or. a2 /= zero .or. comp_exact) then
          deallocate( tmp, stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: dealloc of tmp failed')
       endif
       if (do_grads) then
          deallocate( hlp, stat = stat )
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: dealloc of hlp failed')
          if (weight_grads) then
             deallocate( sec_der_rho, stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: dealloc of sec_der_rho failed')
          endif
          if (a2 /= zero) then
             deallocate( delta, stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: dealloc of delta failed')
          endif
       elseif (comp_exact) then
          deallocate( hlp, stat = stat )
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: dealloc of hlp failed')
          if (a2 /= zero) then
             deallocate( delta, stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: dealloc of delta failed')
          endif
       endif
    else ! LDA calculation
       if (do_grads) then
          deallocate( theta, tmp, stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: dealloc of tmp failed')
          if (weight_grads) then
             deallocate( grad_rho, stat = stat)
             if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
                  &_grads: dealloc of grad_rho failed')
          endif
       else if(comp_exact) then
          deallocate( theta, tmp, stat = stat)
          if (stat /= 0) call error_handler('post_scf_calc_xcmda_energy_and&
               &_grads: dealloc of tmp failed')
       endif
    endif

    send_grad_grid = weight_grads .and. options_split_gradients()
    send_mda_coeff = options_xcmode() == xcmode_model_density .or. &
                     options_xcmode() == xcmode_extended_mda

    call ph_sndrcv(do_grads, send_grad_grid, send_mda_coeff)

    if (comm_rank() == 0) then
      if(do_grads.or.comp_exact) then
       ! solve linear equation systems to get b_k,s
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
       ! ---------------------------------------------------------------------
       if (ispin > 1 .and. .not.constrain_rhospin) then
          rhs(:,1) = half * ( rhs(:,1) + rhs(:,2) ) ! V_xc:tot
          rhs(:,2) =          rhs(:,1) - rhs(:,2)   ! V_xc:spin
       endif
       ! The sequence of linsys calls  m u s t  exactly correspond to
       ! the sequence of linsys calls used for the charge density fit, and
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
       if (ispin > 1 .and. .not.constrain_rhospin) then
          coeff_xcmda(:,1) = rhs(:,1) + rhs(:,2) ! V_mda,xc(up)
          coeff_xcmda(:,2) = rhs(:,1) - rhs(:,2) ! V_mda,xc(down)
          coeff_xcmda_ts(:,:) = rhs(:,:)
       else
          coeff_xcmda = rhs
          if(ispin>1) then
                coeff_xcmda_ts(:,1)=half*(rhs(:,1) + rhs(:,2))
                coeff_xcmda_ts(:,2)=half*(rhs(:,1) - rhs(:,2))
          else
                coeff_xcmda_ts(:,:)=rhs(:,:)
          endif
       endif

       if(comp_exact) then

        allocate(theta(n_ch,2))
        if(ispin>1) then
                theta(:,1)=(coeff_charge(:)+coeff_spin(:))/2.0_r8_kind
                theta(:,2)=(coeff_charge(:)-coeff_spin(:))/2.0_r8_kind
        endif
        do i=1,n_ch
        write(output_unit,*) coeff_xcmda(i,:ispin),theta(i,:ispin)
        enddo

            ! FIXME: why do we need incomplete gamma function in
            !        numerical XC code?
            if(gamma_is_closed()) call gamma_setup(16)
            ! cp. 16 with default numj=17 in old vers. of gamma_module
            call grid_loop_setup_atom() !initalizing jobs distribution
                                ! needed before every loop
            do while (more_grid_atom(vec_length, i, grdpts, grdwts)) ! loop over gridpoints
                ! fetching part of the grid on i, only gets false if dlb cannot steal anymore
               i_m = unique_atoms(i)%moving_atom
               moving_grid = do_grads .and. weight_grads .and. i_m > 0
                    vec_length_act=size(grdpts,1)
                    vl = vec_length_act
                    ! FIXME: there seems to be no usefull code in this loop?
!AG temp comment                        call rspace_compare(grdpts,vl,coeff_xcmda,n_ch,"xc ",grdwts)
!AG temp comment                        call rspace_compare(grdpts,vl,theta,n_ch,"rho",grdwts)
            enddo
            deallocate(theta)
       endif
      end if !do_grads
    end if

  end subroutine post_scf_calc_xcmda_etot_grads

!AG ] End Insert from pg_V2.1_MDA (subroutine post_scf_calc_xcmda_etot_grads)

  subroutine dtrace (msg)
    use comm, only: comm_rank
    use iounitadmin_module, only: write_to_trace_unit
    implicit none
    character(len=*), intent(in) :: msg
    ! *** end of interface ***

    if (comm_rank() == 0) then
       call write_to_trace_unit (msg)
    endif
  end subroutine dtrace

end module post_scf_module

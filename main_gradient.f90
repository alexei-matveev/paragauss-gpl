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
!===============================================================
! Public interface of module
!===============================================================
subroutine main_gradient(loop)
  !----------------------------------------------------------------
  !
  !  Purpose: This is the main routine of the integral part
  !           for calculation of energy gradients
  !           executed by the master.
  !
  !  Subroutine called by: main_master
  !
  !  References: Publisher Document: Concepts of Integral Part
  !
  !
  !  Author: MS
  !  Date: 12/96
  !
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: Uwe Birkenheuer
  ! Date:   June 1998
  ! Description: Moving_Unique_Atom concept introduced
  !              Split_Gradients concept introduced
  !              Gradients for Model_Density_Approach introduced
  !
  !================================================================
  ! End of public interface of module
  !================================================================
!..............................................................................
!
! Individual force contributions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! F_grid     = Sum(a) [ d/dR w_a ] e_xc[rho_up,rho_dn](r_a)
!              Sum(a) w_a [ d/dR r_a ] Nabla e_xc[rho_up,rho_dn](r_a)
!
! F_grid_mda = Sum(a) [ d/dR w_a ] e_xc[rhofit_up,rhofit_dn](r_a)
!              Sum(a) w_a [ d/dR r_a ] Nabla e_xc[rhofit_up,rhofit_dn](r_a)
!
! F_ob_num   = Sum(s) < d/dR rho_s | V_xc,s[rho_up,rho_dn] >_num
!
! F_rho_num  = Sum(s) < d/dR rho_s | V_xc,s[rhofit_up,rhofit_dn] >_num
!
! F_ch_fit   = - [ rhofit | d/dR rhofit ]
!
! F_rho_fit  = - Sum(s) < d/dR rhofit_s | V_X,s >
!
! F_vxc_fit  = - Sum(s) < rhofit_s | d/dR V_X,s >
!
! F_ob_3c    = Sum(n,s) < d/dR psi_n,s | T+V_nuc+V_H - eps_n,s | psi_n,s >
!
! F_ob_mda   = Sum(n,s) < d/dR psi_n,s | H_KS,s - eps_n,s | psi_n,s >
!
! F_ch_3c    = [ rho | d/dR rho_fit ]
!
! F_vxc_3c   = Sum(s) < rho_s | d/dR V_X,s >
!
! F_hf_3c    = < rho_tot | d/dR V_nuc >
!
! F_nn       = d/dR E_nn
!
! Physical force contributions
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 1. Hellmann-Feyman force              : F_hf  = F_hf_3c + F_nn
! 2. Orbital Pulay correction           : F_ob  = F_ob_3c + F_ob_num [non-MDA]
!                                               = F_ob_mda           [    MDA]
! 3. Charge fit Pulay correction        : F_ch  = F_ch_3c + F_ch_fit
! 4. XCMDA correction due to d/dR V_X,s : F_vcx = F_vxc_3c + F_vxc_fit
! 5. XCMDA correction due to d/dR rho_s : F_rho = F_rho_num + F_rho_fit
! 6. Numerical grid correction          : F_num = F_grid             [non-MDA]
!                                                 F_grid_mda         [    MDA]
!
! Storage
! ~~~~~~~
!
!              total_gradients          split_gradients              type
! ----------------------------------------------------------------------------
! F_grid       -grad_xc                 -grad_grid                   non-MDA
! F_grid_mda   -grad_xc                 -grad_grid                   MDA
! F_ob_num     -grad_xc                 -grad_xc                     non-MDA
! F_rho_num    -grad_xc                 -grad_xc                     MDA
! F_ch_fit     -grad_fit_ch_cartesian   -grad_fit_ch_cartesian       both
! F_rho_fit    -grad_fit_ch_cartesian   -grad_mda_rhofit_cartesian   MDA
! F_vxc_fit    -grad_fit_ch_cartesian   -grad_mda_xcpot_cartesian    MDA
! F_ob_3c       gradient_totalsym        gradient_ob_pulay           non-MDA
! F_ob_mda      gradient_totalsym        gradient_ob_pulay           MDA
! F_ch_3c       gradient_totalsym        gradient_ch_pulay           both
! F_vxc_3c      gradient_totalsym        gradient_mda_vxc            MDA
! F_hf_3c       gradient_totalsym        gradient_totalsym           both
! F_nn          gradient_totalsym        gradient_totalsym           both
!
! F_hf          gradient_totalsym        gradient_totalsym           both
! F_ob          gradient_totalsym        gradient_ob_pulay           both
! F_ch          gradient_totalsym        gradient_ch_pulay           both
! F_vxc         gradient_totalsym        gradient_mda_vxc            MDA
! F_rho         gradient_totalsym        gradient_mda_rho            MDA
! F_num        -grad_xc                 -grad_grid                   both
!..............................................................................
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: A.Hu
  ! Date:   4/99
  ! Description: added pseudopotentials contributions to gradients
  !
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   6/97
  ! Description: pvm -> comm
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
  use type_module          ! type specification parameters
  use efield_module, only: efield_gradient, efield_applied
  use efield_module, only: efield_intensity, efield_gradient_store
  use efield_module, only: efield_send_recv
#ifdef _COMPAC_FORTRAN
  use datatype
#endif
  use output_module        ! defines amount of output
  use iounitadmin_module   ! to open output units
  ! comm related information and routines
  use comm, only: comm_rank, comm_size, comm_reduce
  use integralpar_module   ! steering information for integral part
  use operations_module, only: operations_geo_opt,operations_post_scf, &
                               operations_core_density,operations_solvation_effect, &
#ifdef WITH_MOLMECH
                               operations_qm_mm_new, &
#endif
                               operations_qm_mm, operations_integral
  use gradient_data_module, only : gradient_data_n_gradients                   &
                                 , gradient_totalsym                           &
                                 , gradient_cartesian                          &
                                 , gradient_ob_pulay                           &
                                 , gradient_ch_pulay                           &
                                 , gradient_mda_vxc                            &
                                 , gradient_mda_rho                            &
                                 , gradient_index                              &
                                 , grad_fit_ch_cartesian                       &
                                 , grad_mda_rhofit_cartesian                   &
                                 , grad_mda_xcpot_cartesian                    &
                                 , dervs_totalsym                              &
                                 , dervs_cartesian                             &
                                 , cpks_gradient_totalsym                      &
                                 , cpks_grad_fit_totasym                       &
                                 , cpks_grad_fit_ch_cartesian                  &
                                 , partial_cartesian                           &
                                 , gradient_data_setup                         &
                                 , gradient_data_write_cartesians              &
                                 , gradient_data_write_cart_hess               &
                                 , gradient_sndrcv_fit_ch                      &
                                 , gradient_sndrcv_3c                          &
                                 , send_receive_q_grads                        &
                                 , gradient_sndrcv_3c                          &
                                 , gradient_data_write_gxfile                  &
                                 , gradient_data_shutdown                      &
                                 , cpks_add_fitmat_ch_grads
#ifdef WITH_EPE
  use gradient_data_module, only: calc_cluster_epe_energy, cluster_epe_energy
#endif
#ifdef WITH_MOLMECH
  use gradient_data_module, only: qm_grads_to_qmmm
#endif

  use occupied_levels_module
  use xc_cntrl, only: xc_is_on=>is_on, xc_ANY
  use post_scf_module, only: post_scf_deallocate_grad_xc,grad_xc,grad_grid, &
                             dervs_xc
  use grid_module, only: weight_grads
  use fit_coeff_module, only: fit_coeff_sndrcv, ICHFIT, IXCFIT, coeff_charge
  use options_module, only : options_relativistic, options_split_gradients, &
                             options_xcmode, xcmode_model_density           &
                           , options_debug_key
  use unique_atom_module, only: unique_atoms, N_unique_atoms, &
      N_moving_unique_atoms, moving_unique_atom_index, &
      unique_atom_grad_info
  use pointcharge_module, only: moving_pc, print_pc_grad, &
       pc_grad_cart_write, transform_pc_grad_to_cart
#ifdef WITH_EPE
  use ewaldpc_module, only: cluster_nuc_epe_en, epe_relaxation
#endif
  use point_dqo_module, only: dealloc_dqo,calc_nuc_dqo_grads,dealloc_center_inform,moving_X_centers, &
       calc_X_grads,print_X_grad,transform_x_grad_to_cart,x_grad_cart_write
  use point_dqo_module, only: moving_R_centers,print_R_grad
  use induced_dipoles_module, only: moving_Pol_centers, print_id_grad
  use calc_id_module, only: calc_nuc_id_grads,calc_id_grads
  use calc_id_module, only: transform_id_grad_to_cart,id_grad_cart_write
#ifdef WITH_EFP
  use efp_rep_module, only: dealloc_repf,dealloc_repf_inform,transform_repf_grad_to_cart, &
       repf_grad_cart_write
  use efp_efp_module, only: efp_efp_gradients
  use efp_module, only: efp, n_efp, efp_sum_up_grads, efp_grad_cart_write, efp_write_gxfile
  use efp_module, only: qm_fixed
  use efp_only_opt_module, only: efp_opt_main
  use energy_calc_module, only: get_energy
  use efp_solv_grad_module, only: transform_X_solv_grad_to_cart, X_solv_grad_cart_write
#endif
  use solv_cavity_module
  use solv_electrostat_module
  use solv_2nd_deriv_module, only: nuc_solv_2nd_deriv,charge_solv_2nd_deriv,matrix_2nd_deriv !!!!!!!!!AS
#ifdef WITH_MOLMECH
  use qmmm_interface_module  !!!!!!!!!!!!!AS
  use qmmm1_interface_module !!!!!!!!!!!!AS
#endif
  use datatype
  use calc3c_switches
  use cpksdervs_matrices, only: cpksalloc
#ifdef WITH_SECDER
  use cpks_utils, only: cpks_bcast_hr1 &
                      , cpks_free_hr1
#endif
#ifdef WITH_DFTPU
  use dft_plus_u_module,   only: dft_plus_u_grad_init                          &
                               , dft_plus_u_grad_finalize                      &
                               , dft_plus_u_mo_grad_init
#endif
  use density_data_module, only: gendensmat_occ, density_data_free
#ifdef WITH_ERI4C
  use density_data_module, only: densmat
#endif
  use overlap_module, only: overlap
  use error_module, only: MyID
  use interfaces, only: main_integral, integral_trafo, RELGRAD, RELSDER
  use interfaces, only: potential_calculate
  use interfaces, only: grad_solv_calculate
  use occupation_module, only: dealloc_occ_num
#ifdef WITH_EXPERIMENTAL
  use symmetry_data_module, only: irrep_dimensions
  use cpks_common, only: cpks_alloc_p1w1, cpks_free_p1w1, cpks_p1, cpks_w1
  use cpks_utils, only: cpks_calc_p1w1
#endif
#if WITH_ERI4C == 1
  use eri4c_options, only: J_exact, K_exact
#endif
  implicit none


  integer(i4_kind), intent(in) :: loop ! used only on master so far
  ! *** end of interface ***

  integer :: rank
  logical :: model_density, split_gradients
  integer (kind=i4_kind) :: i, n_equal_max
! real(kind=r8_kind),allocatable :: gradient_totalsym_fit(:)
  integer(i4_kind) :: IFIT
  integer(i4_kind) :: IREL
  logical :: tty                ! used in say() subprogram
  logical :: xc ! indicates if XC is required at all
  integer, parameter :: grads=1, deriv=2
#ifdef WITH_EFP
  real (r8_kind) :: energy
#endif

  FPP_TIMER_DECL(tot)
  FPP_TIMER_DECL(int1)
  FPP_TIMER_DECL(rel)
  FPP_TIMER_DECL(cpks)
  FPP_TIMER_DECL(int2)
  !----------------------------------------------------------------

  !------------ Executable code -----------------------------------

  ! ================ CONTEXT: ALL PROCESSORS ===============
  FPP_TIMER_START(tot)
  DPRINT MyID,"main_gradient:  start"

  rank = comm_rank ()

  ! Tty is  used in subrogram  say() to determine whether  the current
  ! processor should print the given output:
  tty = output_main_gradient .and. rank == 0

  ! indicates if XC is required at all:
  xc  = xc_is_on(xc_ANY) .and. operations_post_scf

#ifdef WITH_EFP
  if(efp .and. n_efp > 1 .and. qm_fixed) then
     DPRINT MyID,"main_gradient: starting EFP optimization"
     call efp_opt_main()
     DPRINT MyID,"main_gradient: finishing EFP optimization"
     return
  end if
#endif

  ! ================ CONTEXT: ALL PROCESSORS ===============
  split_gradients = options_split_gradients()
  model_density = options_xcmode() == xcmode_model_density

  call say ("start")

  DPRINT MyID,"main_gradient: call gradient_data_setup() "
  call say ("gradient_data_setup")
  !
  ! See matching gradient_data_shutdown() below:
  !
  call gradient_data_setup()
  DPRINT MyID,"done"

  ! ================ CONTEXT: MASTER ONLY ==================
  master1: if (rank == 0) then

     solv_eff: if(operations_solvation_effect) then
        ! Calculations  of the additional  contributions to  the total
        ! molecular  gradients  due  to  the  solvent  (nonelectronic)
        ! effect
        do_gradients=.true.
        if(disp_rep_energy) then
           do_cavitation=.false.
           do_disp_rep=.true.
           call disp_rep_wrap()
           if(output_solv_grads) then
              call transform_to_cart_grads(grad_solv_cart,grad_solv_totsym)
              call gradient_data_write_cartesians( &
                   ' Cartesian Solvation Gradients - Disp - Rep term:',grad_solv_cart)
              call add_solv_grads()
              grad_solv_totsym=0.0_r8_kind
              do i=1,N_moving_unique_atoms
                 grad_solv_cart(i)%m=0.0_r8_kind
              enddo
           endif
        endif

        if(cavitation_energy) then
           do_cavitation=.true.
           do_disp_rep=.false.
           call points_on_cavity_surface()

           ! This sets  energy_cav in solv_cavity_module,  among other
           ! things.  If  any  of the  operations_solvation_effect  or
           ! cavitation_energy   is  not   set  that   energy  remains
           ! uninitialized:
           call energy_and_grad_of_cavity()

           if(output_solv_grads) then
              call transform_to_cart_grads(grad_solv_cart,grad_solv_totsym)
              call gradient_data_write_cartesians( &
                   ' Cartesian Solvation Gradients - Cavitation term:',grad_solv_cart)
              call add_solv_grads()
              grad_solv_totsym=0.0_r8_kind
              do i=1,N_moving_unique_atoms
                 grad_solv_cart(i)%m=0.0_r8_kind
              enddo
           endif
        endif
        call say ("solvation setup")
        do_cavitation=.false.
        do_disp_rep=.false.
        call points_on_cavity_surface
        call say ("solvation setup done")
     endif solv_eff

#ifdef WITH_EFP
     if(efp .and. n_efp > 0) call efp_efp_gradients()
#endif
  endif master1
  ! ================ CONTEXT: ALL PROCESSORS ===============
  if(.not. operations_integral) goto 1000 ! not to calculate QM gradients

     ! FIXME:  yet better, make  sure that  fit_coeff_sndrcv() handles
     ! the serial case gracefully:
     if (comm_size () > 1) then
        call say ("fit_coeff_sndrcv")
        ! MDA:  update coeff_xcmda must  be sent  to each  slave else:
        ! coeff_charge has not yet been send to the slaves
        IFIT = ICHFIT
        if (model_density) IFIT = IFIT + IXCFIT
        call fit_coeff_sndrcv (IFIT)
     endif

     call say ("main_integral(1)")

     if(options_relativistic) then
        call integralpar_set('RelGrads')
     else
        call integralpar_set('Gradients')
     end if

     !********Initial steps before computing the gradients in presence of electric field*****!
     ! Passing the value of eletric field strength to master and slave.
     !
     call efield_send_recv()


     ! Computation of 2e integral derivatives with ERI4C library
     ! TODO: AVOID GENERATING DENSITY MATRIX AGAIN....
#if WITH_ERI4C == 1
     IF ( K_exact .or. J_exact ) THEN
       call gendensmat_occ()
       !
       call exact_2e_grads( densmat, gradient_totalsym )
       ! TODO: in the case of enabled DFT+U densmat is generated again by init
       call density_data_free()
     ENDIF
#endif

#ifdef WITH_DFTPU
     !
     ! Inital steps before computing DFT+U gradients:
     !
     call dft_plus_u_grad_init()
     call dft_plus_u_mo_grad_init()
#endif

     ! DONE by integralpar_set(...):
     ! integralpar_cpks_contribs=.false. ! explicit functional dervs only

     FPP_TIMER_START(int1)
     call main_integral () ! (1) First call: gradients, maybe second derivatives
     FPP_TIMER_STOP(int1)
     DPRINT MyID,'main_gradient:  calc_3c=',FPP_TIMER_VALUE(t_calc_3center)
     DPRINT MyID,'main_gradient:  int1=',FPP_TIMER_VALUE(int1)

#ifdef WITH_DFTPU
     !
     ! Final steps after computing DFT+U gradients,
     ! Deallocate internal structures and 'densmat'
     !
     call dft_plus_u_grad_finalize()
#endif

     if (integralpar_2dervs .and. rank == 0) then
        ! === context: master only ===
          ! orbital gradient_totalsym contribs calculated
          ! to be summed up with fit_charge contribs and
          ! to be transformed to cartesians
        call secder_add_ch_grads() !ito calculate cpks_gradient_totalsym,cpks_grad_fit_totasym, see contains
     endif

     call say ("gradient_sndrcv_fit_ch")
     ! GET grad_fit_ch_cartesian     : < rhofit | d/dR V_H[rhofit] >
     ! GET grad_mda_rhofit_cartesian : < d/dR rhofit | V_X[rhofit] >
     ! GET grad_mda_xcpot_cartesian  : < rhofit | d/dR V_X[rhofit] >
     call gradient_sndrcv_fit_ch()   !  msgtag_grad_ch

   ! ================ CONTEXT: MASTER ONLY ==================
   master2: if (rank == 0) then

     call say ("add_fit_ch_grads")

     sum_fit_ch_grads: if (split_gradients) then
        call add_fit_ch_grads(gradient_ch_pulay,grad_fit_ch_cartesian)
        if (model_density) then
           call add_fit_ch_grads(gradient_mda_rho,grad_mda_rhofit_cartesian)
           call add_fit_ch_grads(gradient_mda_vxc,grad_mda_xcpot_cartesian)
        endif

     else ! i.e. .not.split_gradients
#define IMPL_CONTRIBS_TO_DERVS
#ifdef IMPL_CONTRIBS_TO_DERVS
      call add_fit_ch_grads(gradient_totalsym,grad_fit_ch_cartesian)  !fit(1)
#endif
               ! totalsym func contrib, for coeff grad part see below
               ! rho contrib only if commented
     endif sum_fit_ch_grads

     call say ("core_gradient_calc")
     call core_gradient_calc() !GET d/dR E_nn
     DCALL print_totsym_grads('GRAD: gradient_totalsym fin orbital and fit contribs')
   endif master2
   ! ================ CONTEXT: ALL PROCESSORS ===============

   ![[=============== Relativistic Transformations: =========
   if(options_relativistic) then

     IREL = RELGRAD

     if(integralpar_2dervs)  then
       IREL = IREL + RELSDER
     endif

     DCALL print_totsym_grads('GRAD: gradient_totalsym before integral_trafo')
     call say ("start rel trafos")

     FPP_TIMER_START(rel)
     call integral_trafo(IREL)
     FPP_TIMER_STOP(rel)
     DPRINT MyID,'main_gradient:   rel=',FPP_TIMER_VALUE(rel)
     DCALL print_totsym_grads('GRAD: gradient_totalsym after integral_trafo')
   endif
   !]]=======================================================

1000 continue ! .not. operations_integral
     if(operations_solvation_effect) then
        call say ("grad_solv_calculate ")
        call grad_solv_calculate('SolvGrads')
        if(integralpar_2dervs) then
           call say ("grad_solv_calculate - Q deriv ")
           call grad_solv_calculate('Q_SolvGrads')
        end if

        if (rank == 0) then
           if(VTN) then
              call say ("solvation nuc_grad_vtn")
              call nuc_grad_vtn(gradient_index)
           else
              call say ("solvation nuc_grad")
              call nuc_grad(gradient_index,integralpar_2dervs)
           end if
        end if
        call say ("solvation matrix_grad")
        if(VTN) then
           call matrix_grad_vtn(gradient_index)
        else
           call matrix_grad(gradient_index)
        end if

        if(integralpar_2dervs) then
           if (rank == 0) then
              call say ("nuc_solv_2nd_deriv ")
              call nuc_solv_2nd_deriv()
           end if
        end if
     endif
  if(.not. operations_integral) goto 1001 ! not to calculate QM gradients

     call say ("gradient_sndrcv_3c")
     ! GET gradient_ob_pulay : Sum(n) < d/dR psi_n | H_KS - eps_n | psi_n >
     !                         without the numerical V_xc contribution
     ! GET gradient_ch_pulay : < rho | d/dR V_H[rhofit] >
     ! GET gradient_mda_vxc  : < rho | d/dR V_X[rhofit] >
     ! GET gradient_totalsym : < rho | d/dR V_nuc >
     call gradient_sndrcv_3c(grads)

     if (operations_solvation_effect) then
        if (integralpar_2dervs) then
#if 1 /* def IMPL_CONTRIBS_TO_DERVS */
           if (rank == 0) then
              call charge_solv_2nd_deriv()
           end if
#endif
           ! FIXME: better  yet make sure  that send_receive_Q_grads()
           ! handles the serial case:
           if (comm_size () > 1) then
              call send_receive_Q_grads ()
           endif
           call matrix_2nd_deriv()
           call potential_calculate ('SolvDervs')
        end if
     endif

     call gradient_sndrcv_3c(deriv)


   ! ================ CONTEXT: MASTER ONLY ==================
   master3: if (rank == 0) then

#ifdef IMPL_CONTRIBS_TO_DERVS
     if( xc )then
       call add_xc_grads(gradient_cartesian,grad_xc)
       DCALL print_cart_grads('GRAD: gradient_cartesian only xc contribs')
     endif
#endif

     call say ("transform_to_cart_grads")
     call transform_to_cart_grads(gradient_cartesian,gradient_totalsym)

     ! below xc gradients were previously added  directly to gradient_cartesians
     ! now gradient calc finished but gradient_totalsym will be modified futher
     ! in cpks calcs

     DCALL print_cart_grads('GRAD: gradient_cartesian with xc contribs')
  endif master3
  ! ================ CONTEXT: ALL PROCESSORS ===============

#ifdef WITH_SECDER
    if(integralpar_cpksdervs)  then
      DPRINT 't_cpksdervs_quad_sums= ',FPP_TIMER_VALUE(t_dervs_quad_sums)

      DPRINT MyID//'cpks_g4constructs_______________________________________________'
      FPP_TIMER_START(cpks)
      call cpks_g4constructs() ! -> parametrization matrix
      FPP_TIMER_STOP(cpks)
      DPRINT MyID,'main_gradient: cpks=',FPP_TIMER_VALUE(cpks)

      ! distribute saved rel-ham:
      if(options_relativistic) then
         ! do it before second call to main_integral()
         ! in a sec-der run.
         call cpks_bcast_hr1()
      endif
#if 0
      ! Add fit contribution aX(k) * gY(k,l) * a(l) (moved from main_integral)
      call impl_fit_dervs(dervs_totalsym)
#endif

#ifdef WITH_EXPERIMENTAL
      ! allocate cpks_p1/cpks_w1 matrices from cpks_common:
      call cpks_alloc_p1w1(irrep_dimensions,n_spn=1,n_gra=size(cpks,1))
      ! cf. cpks_free_p1w1

      ! Compute (yet another) density matrix and energy weited density matrix:
      call cpks_calc_p1w1(1,p1=cpks_p1,w1=cpks_w1)
#endif

      ! FIXME: slave have it allocated?
      cpks_gradient_totalsym=0.0_r8_kind
                            ! this term will be recalculated
                            ! to be multiplied on fitcoeff gradients or
                            ! it should not be modified after initial
                            ! calculation

      ! no need for RelGrads as they were saved before CPKS:
      call integralpar_set('Gradients2') ! sets integralpar_cpks_contribs=.true.

      call say ("main_integral(2)")

      FPP_TIMER_START(int2)
      call main_integral ()
      FPP_TIMER_STOP(int2)
      DPRINT MyID,'main_gradient:  calc_3c=',FPP_TIMER_VALUE(t_calc_3center)
      DPRINT MyID,'main_gradient:  int2=',FPP_TIMER_VALUE(int2)
      DPRINT MyID,'main_gradient:  w1_sumup=',FPP_TIMER_VALUE(t_w1_sumup)
      DPRINT MyID,'main_gradient:  p1_sumup=',FPP_TIMER_VALUE(t_p1_sumup)
      ! no regular gradient_totalsym contibs are required for this call

#if 1
      ! Add fit contribution aX(k) * gY(k,l) * a(l) (moved from main_integral)
      call impl_fit_dervs(dervs_totalsym)
      DPRINT MyID//'impl_fit_dervs done'
#endif

      if (operations_solvation_effect) then
         call say ("grad_solv_calculate: cpks contribution")
         call grad_solv_calculate('SolvGrads2')
         call potential_calculate ('SolvDervs2')
         ! NOTE: potential_calculate also calls main_integral()
         !       that does solvation integrals
      endif

#ifdef WITH_EXPERIMENTAL
      ! deallocate cpks_p1/cpks_w1 matrices from cpks_common:
      call cpks_free_p1w1() ! cf. cpks_alloc_p1w1
#endif

      ! deallocate saved rel-ham:
      if(options_relativistic) then
         ! do it after second call to main_integral(Gradients2)
         ! in a sec-der run.
         call cpks_free_hr1()
      endif

      ! Sum  partial contributions  to second  derivatives  on master.
      ! FIXME:  should  we  do   allreduce  instead  to  minimize  the
      ! differences between workers?
      call comm_reduce (dervs_totalsym)
   endif ! integralpar_cpksdervs
#endif


1001 continue ! .not. operations_integral !!!!!!!!!!!!!
   ! ================ CONTEXT: MASTER ONLY ==================
   master4: if (rank == 0) then

    QM2: if(operations_integral) then ! if QM gradients exist
       !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#ifdef WITH_EPE
     if(epe_relaxation.and.calc_cluster_epe_energy) then
        DPRINT 'main_gradient: cluster_epe_energy',cluster_epe_energy
        write(output_unit,*)'energy of interaction of the exact orbital electronic'
        write(output_unit,*)'density with epe centers only,cluster_epe_energy'
        write(output_unit,*)'main_gradient: cluster_epe_energy',cluster_epe_energy
        write(output_unit,*)
        write(output_unit,*)'energy of interaction of cluster and epe centers only'
        write(output_unit,*)cluster_epe_energy- cluster_nuc_epe_en
        write(output_unit,*)
     end if
#endif


     if(integralpar_2dervs)  then
      call transform_to_cart_dervs(dervs_cartesian,dervs_totalsym)
      DPRINT MyID//'transform_to_cart_dervs done'
     endif

!    if(integralpar_cpksdervs) &
!     call cpks_transform_to_cart_grads(cpks_gradient_cartesian,cpks_gradient_totalsym)

     if(operations_solvation_effect) then
        call transform_to_cart_grads(grad_solv_cart,grad_solv_totsym)
        if(output_solv_grads) then
           call gradient_data_write_cartesians( &
                ' Cartesian Solvation Gradients - Electrostatic term:',grad_solv_cart)
        endif
        call add_solv_grads()
        DCALL print_cart_grads('GRAD: gradient_cartesian with solv contribs')
        if(with_pc .and. .not.fixed_pc) call solv_forces_on_pc()
     endif

     splig_cart: if (split_gradients) then
        ! 1. Hellmann-Feyman contribution
        call gradient_data_write_cartesians( &
             ' 1. The Hellmann-Feyman force contributions', &
             gradient_cartesian)
        ! 2. orbital Pulay correction
        call transform_to_cart_grads(partial_cartesian,gradient_ob_pulay)
        if (.not.model_density) then
           ! GET grad_xc : - < d/dR rho | Vxc[rho] >_num
           if(xc) call add_xc_grads(partial_cartesian,grad_xc)
        endif
        call gradient_data_write_cartesians( &
             ' 2. The orbital Pulay correction to the forces', &
             partial_cartesian)
        call sum_up_and_reset_part_grads()
        ! 3. charge fit Pulay correction
        call transform_to_cart_grads(partial_cartesian,gradient_ch_pulay)
        call gradient_data_write_cartesians( &
             ' 3. The charge fit Pulay correction to the forces', &
             partial_cartesian)
        call sum_up_and_reset_part_grads()
        if (weight_grads) then
           ! 4. numerical grid correction
           ! GET grad_grid : - Sum(a) (d/dR w_a) exc[rho](r_a)
           !                 - Sum(a) w_a (d/dR r_a) Nabla exc[rho](r_a)
           if(xc) call add_xc_grads(partial_cartesian,grad_grid)
           call gradient_data_write_cartesians( &
                ' 4. The numerical grid correction to the forces', &
                partial_cartesian)
           call sum_up_and_reset_part_grads()
        endif

        if (model_density) then
           ! 5. XCMDA correction : d/dR V_X[rhofit]
           call transform_to_cart_grads(partial_cartesian,gradient_mda_vxc)
           call gradient_data_write_cartesians( &
                ' 5. The XC force corrections due to d/dR V_xc[rhofit]', &
                partial_cartesian)
           call sum_up_and_reset_part_grads()
           ! 6. XCMDA correction : d/dR rhofit
           call transform_to_cart_grads(partial_cartesian,gradient_mda_rho)
           ! GET grad_xc : - < d/dR rhofit | Vxc[rhofit] >_num
           if(xc) call add_xc_grads(partial_cartesian,grad_xc)
           call gradient_data_write_cartesians( &
                ' 6. The XC force corrections due to d/dR rhofit', &
                partial_cartesian)
           call sum_up_and_reset_part_grads()
        endif
     else splig_cart
      call say ("add_xc_grads")
     ! GET grad_xc   : - < d/dR rho | Vxc[rho] >_num
     ! GET grad_grid : - Sum(a) (d/dR w_a) exc[rho](r_a)
     !                 - Sum(a) w_a (d/dR r_a) Nabla exc[rho](r_a)

!     call add_xc_grads(gradient_cartesian,grad_xc)
     !   grad_xc contribs are directly added to  other cartesian contribs
     !   while othe contribs transforme from the total symmetric contribs
#ifdef WITH_SECDER
     if(integralpar_2dervs) then
      n_equal_max=maxval(unique_atoms(:)%n_equal_atoms)
      if( xc )then
        call average_xc_dervs(dervs_xc)
#if 0
       print *,'XC hessian cartesian'
       call gradient_data_write_cart_hess(stdout_unit,dervs_xc)
#endif
        call add_xc_dervs(dervs_cartesian,dervs_xc)
      DPRINT MyID//'add_xc_dervs done'
      endif
      call store_dervs_cart(dervs_cartesian)
      call gradient_data_write_cart_hess(stdout_unit,dervs_cartesian)
      call gradient_data_write_cart_hess(output_unit,dervs_cartesian)
      DPRINT 't_cpksdervs_quad_sums= ',FPP_TIMER_VALUE(t_dervs_quad_sums)
     endif
#endif
     endif splig_cart

     ! total force due to presence of electric field (nuclear + electronic) contribution
     if (efield_applied()) then
       ! I guess this adds the nuclear contribution to the gradient_cartesian:
       call efield_gradient(gradient_cartesian)
     endif

     call say ("gradient_data_write_cartesians")
     call gradient_data_write_cartesians(' Final Cartesian Gradients:', &
                                            gradient_cartesian)
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     end if QM2

     if(moving_pc) then
        call say ("PC_grad_to_cart")
        call transform_PC_grad_to_cart()
        if(print_pc_grad) call pc_grad_cart_write()
     end if
     if(moving_X_centers .or. moving_R_centers) then
        call say ("X_grad_to_cart")
        call transform_X_grad_to_cart()
        if(print_X_grad .or. print_R_grad) call X_grad_cart_write()
     end if
#ifdef WITH_EFP
     if(moving_R_centers) then
        call say ("repf_grad_to_cart")
        call transform_repf_grad_to_cart()
        if(print_R_grad) call repf_grad_cart_write()
     end if
#endif
     if(moving_Pol_centers) then
        call say ("id_grad_to_cart")
        call transform_id_grad_to_cart()
        if(print_id_grad) call id_grad_cart_write()
     end if
#ifdef WITH_EFP
     if(operations_solvation_effect) then
        if(moving_pc) then
           call say ("X_solv_grad_to_cart")
           call transform_X_solv_grad_to_cart()
           if(output_solv_grads .and. print_pc_grad) call X_solv_grad_cart_write()
        end if
     end if
     if(efp) then
        call say ("efp_sum_up_grads")
        call efp_sum_up_grads()
        call efp_grad_cart_write()
     end if
#endif
  endif master4
  ! ================ CONTEXT: ALL PROCESSORS ===============

     !
     ! FIXME: does it need to be here? Was the allocation intiated
     !        from this scope? Maybe finalize_geometry() is a better
     !        place?
     !
     ! make it unconditional and deallocate only what is there
     call post_scf_deallocate_grad_xc()

  ! ================ CONTEXT: MASTER ONLY ==================
  master5: if (rank == 0) then

#ifdef WITH_MOLMECH
     if(operations_qm_mm_new) then
#ifdef WITH_EFP
        if(efp .and. operations_geo_opt) then
           call get_energy (tot=energy)
           call efp_write_gxfile (energy, gradient_cartesian)
        endif
#endif
        if(imomm) call qm_grads_to_qmmm()
        if(qm_mm .or. qm_mm_1) call qm_grads_to_qmmm1()
     else
#endif
        if (operations_geo_opt .or. operations_qm_mm) then
           call say ("call gradient_data_write_gxfile()")
           call gradient_data_write_gxfile (loop)
           call say ("done gradient_data_write_gxfile()")
        end if
#ifdef WITH_MOLMECH
     end if
#endif
  endif master5
  ! ================ CONTEXT: ALL PROCESSORS ===============

  DPRINT MyID,"main_gradient: call gradient_data_shutdown()"
  !
  ! See matching gradient_data_setup() above:
  !
  call gradient_data_shutdown()

  FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
  DPRINT MyID,'main_gradient: TIMING: total  =',FPP_TIMER_VALUE(tot)
  DPRINT MyID,'main_gradient: TIMING: |-int1 =',FPP_TIMER_VALUE(int1)
  DPRINT MyID,'main_gradient: TIMING: |-rel  =',FPP_TIMER_VALUE(rel)
  DPRINT MyID,'main_gradient: TIMING: |-cpks =',FPP_TIMER_VALUE(cpks)
  DPRINT MyID,'main_gradient: TIMING: |-int2 =',FPP_TIMER_VALUE(int2)
#endif

  call say ("end")

  DPRINT MyID,"main_gradient: end"

contains

#ifdef _DPRINT
     subroutine print_cart_grads(header)
       implicit none
       character(len=*), intent(in), optional :: header
       ! *** end of interface ***

       integer(i4_kind) :: u,e

       if(present(header)) &
         print*,header
         print '(2A4,3A20)','ua','ea','x','y','z'
       do u=1,size(gradient_cartesian)
       do e=1,size(gradient_cartesian(u)%m,2)
         print '(2I4,3F20.10)',u,e,gradient_cartesian(u)%m(:,e)
       enddo
       enddo
     end subroutine print_cart_grads

     subroutine print_totsym_grads(header)
       implicit none
       character(len=*), intent(in), optional :: header
       ! *** end of interface ***

       integer(i4_kind) :: n,a

       if(present(header)) &
         print*,header
         print '(2A4,3A20)','from','to','+2','+3','+1'
       n = size(gradient_totalsym)/3
       do a=0,n-1
         print '(2I4,3F20.10)',3*a+1,3*a+3,gradient_totalsym(3*a+2) &
                                          ,gradient_totalsym(3*a+3) &
                                          ,gradient_totalsym(3*a+1)
       enddo
       if( 3*n < size(gradient_totalsym) )&
         print '(2I4,3F20.10)',3*n+1,3*n+3,gradient_totalsym(3*n+1:)
     end subroutine print_totsym_grads
#endif

  subroutine say (msg)
    implicit none
    character (len=*), intent (in) :: msg
    ! *** end of interface ***

    ! tty is variable of main_gradient (thus global in this module)
    ! it is used to generate only on master output and only if
    ! output is demanded
    if (tty) then
        call write_to_output_units ("main_gradient: " // msg)
        call write_to_trace_unit ("main_gradient: " // msg)
    endif
  end subroutine say

  subroutine secder_add_ch_grads()
    implicit none
    ! *** end of interface ***

     DPRINT 'main_integral(1) cpks_totalsym_fitcharge and dervs'

       if(integralpar_cpksdervs) &
        call cpks_add_fitmat_ch_grads( cpks_gradient_totalsym,cpks_grad_fit_totasym, &  !result matrixes
                                       cpks_grad_fit_ch_cartesian)
  end subroutine secder_add_ch_grads

  subroutine store_dervs_cart(dervs_cartesian)
    use filename_module, only: inpfile
    implicit none
    type(arrmat4), intent(in) :: dervs_cartesian(:,:) ! (NUA,NUA)
    ! *** end of interface ***

    integer(kind=i4_kind) :: io_hesse_cart,n_atoms,i_unique,j_unique,i_equal,j_equal, &
    i_off,j_off
    real(kind=r8_kind), allocatable:: hesse_cartesian(:,:)

!    print*, 'openget_iounit hesse_cartesian.dat in main_gra'
    io_hesse_cart=openget_iounit(status='unknown',form='formatted',file=&
                                 trim(inpfile('hesse_cartesian.dat')))
!    print*, 'done openget_iounit hesse_cart'
    n_atoms=0
    uniques: do i_unique=1,N_moving_unique_atoms
       n_atoms = n_atoms +  unique_atoms(moving_unique_atom_index(i_unique))%n_equal_atoms
    enddo uniques
    allocate(hesse_cartesian(3*n_atoms,3*n_atoms),stat=cpksalloc(132))
    ASSERT(cpksalloc(132).eq.0)
    hesse_cartesian=0.0_r8_kind

    i_off=0
    uniques_1: do i_unique=1,N_moving_unique_atoms
    do i_equal=1,unique_atoms(moving_unique_atom_index(i_unique))%n_equal_atoms
    j_off=0
    uniques_2: do j_unique=1,N_moving_unique_atoms
    do j_equal=1,unique_atoms(moving_unique_atom_index(j_unique))%n_equal_atoms
    hesse_cartesian(j_off+1:j_off+3,i_off+1:i_off+3)= &
          dervs_cartesian(j_unique,i_unique)%m(:,j_equal,:,i_equal)
    j_off=j_off+3
    enddo
    enddo uniques_2
    i_off=i_off+3
    enddo
    enddo uniques_1

    write(io_hesse_cart,*) hesse_cartesian
    call returnclose_iounit(io_hesse_cart)

   deallocate(hesse_cartesian,stat=cpksalloc(132))
   ASSERT(cpksalloc(132).eq.0)
   cpksalloc(132)=1

  end subroutine store_dervs_cart

  subroutine average_xc_dervs(dervs_final)
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
                  +help_arr_local(:,1:n_ea,k,i_equal2)*rotmat_tot(kk,k)
!                  +help_arr_local(:,1:n_ea,kk,i_equal2)*rotmat_tot(kk,k)
              enddo
             enddo
       enddo
      enddo

     enddo uniques2
    enddo uniques1
  end subroutine average_xc_dervs



subroutine transform_to_cart_grads(grad_cart,grad_totsym)
  use pointcharge_module
  ! purpose: transform symmetry adapted gradient components to cartesian
  !          coordinates and add them to an array of cartesian gradients
 type(arrmat2)     , intent(inout) :: grad_cart(:)   ! cartesian grad. array
 real(kind=r8_kind), intent(in   ) :: grad_totsym(:) ! symm. adapt. contr.
 !
 integer(kind=i4_kind) :: i_unique,i_center,i_equal,index,i,grad_dim
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    index=index-1
    do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
       do i=1,grad_dim
          grad_cart(i_unique)%m(:,i_equal) = &
               grad_cart(i_unique)%m(:,i_equal) + &
               rotmat(i,:)*grad_totsym(index+i)
!       print*,grad_cart(1)%m(:,1),'in'
       end do
    end do
 enddo

 if(N_moving_unique_timps.ne.0) then
    do i_unique=1,N_moving_unique_timps
       i_center=moving_unique_timp_index(i_unique)
       index=gradient_index(N_moving_unique_atoms+i_unique)
       grad_dim=gradient_index(N_moving_unique_atoms+i_unique+1)-index
       index=index-1
       do i_equal=1,unique_timps(i_center)%n_equal_atoms
          rotmat=>unique_timp_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             grad_cart(i_unique+N_moving_unique_atoms)%m(:,i_equal) = &
                  grad_cart(i_unique+N_moving_unique_atoms)%m(:,i_equal) + &
                  rotmat(i,:)*grad_totsym(index+i)
          end do
       end do
    enddo
 end if

end subroutine transform_to_cart_grads

subroutine transform_to_cart_dervs(dervs_cart,dervs_totsym)

  use pointcharge_module

  ! purpose: transform symmetry adapted gradient components to cartesian
  !          coordinates and add them to an array of cartesian gradients
 type(arrmat4)     , intent(inout) :: dervs_cart(:,:)   ! cartesian dervs. array
 real(kind=r8_kind), intent(in   ) :: dervs_totsym(:,:) ! symm. adapt. contr.
 !
 integer(kind=i4_kind) :: i_unique,i_center,i_equal,index,i,grad_dim
 integer(kind=i4_kind) :: grad_dim_tot,k,i_unique2
 real(kind=r8_kind),pointer :: rotmat(:,:)

 type(arrmat3) :: help_totsym_cart(N_moving_unique_atoms)


!allocate(help_totsym_cart(N_moving_unique_atoms),stat=cpksalloc(88))
!ASSERT(cpksalloc(88).eq.0)
 grad_dim_tot=size(dervs_totsym,1)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    index=index-1

    allocate( help_totsym_cart(i_unique)%m(3,unique_atoms(i_center)%n_equal_atoms, &
                                          grad_dim_tot),  stat=cpksalloc(89))
    ASSERT(cpksalloc(89).eq.0)
    help_totsym_cart(i_unique)%m=0.0_r8_kind

   ! transform by first dervs_totsym index to uniques and equals
    do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
       do i=1,grad_dim
          help_totsym_cart(i_unique)%m(:,i_equal,:) = &
          help_totsym_cart(i_unique)%m(:,i_equal,:) + &
               spread(rotmat(i,:),2,grad_dim_tot)* &
                    spread(dervs_totsym(index+i,:),1,3)
       end do
    end do
!    print*,sum(help_totsym_cart(i_unique)%m),i_unique,'help_totsym_cart,i_unique'
 enddo

 uas: do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    index=index-1


    equals: do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
       totsyms: do i=1,grad_dim
       ua2: do i_unique2=1,N_moving_unique_atoms
            cart: do k=1,3
            dervs_cart(i_unique2,i_unique)%m(:,:,k,i_equal) = &
               dervs_cart(i_unique2,i_unique)%m(:,:,k,i_equal) + &
            rotmat(i,k)*help_totsym_cart(i_unique2)%m(:,:,index+i)
            enddo cart
       enddo ua2
       enddo totsyms
    enddo equals
 enddo uas

 do i_unique=1,N_moving_unique_atoms
  deallocate(help_totsym_cart(i_unique)%m,stat=cpksalloc(89))
  ASSERT(cpksalloc(89).eq.0)
  cpksalloc(89)=1
 enddo

!deallocate(help_totsym_cart,stat=cpksalloc(88))
!ASSERT(cpksalloc(88).eq.0)
!cpksalloc(88)=1
end subroutine transform_to_cart_dervs

subroutine cpks_transform_to_cart_grads(grad_cart,grad_totsym)
  use pointcharge_module
  ! purpose: transform symmetry adapted gradient components to cartesian
  !          coordinates and add them to an array of cartesian gradients
 type(arrmat3)     , intent(inout) :: grad_cart(:)   ! cartesian grad. array
 real(kind=r8_kind), intent(in   ) :: grad_totsym(:,:) ! symm. adapt. contr.
 !
 integer(kind=i4_kind) :: i_unique,i_center,i_equal,index,i,grad_dim,k,n_ch
 real(kind=r8_kind),pointer :: rotmat(:,:)
 n_ch=size(coeff_charge)
 ua1: do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    index=index-1
    do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       rotmat=>unique_atom_grad_info(i_unique)%m(:,:,i_equal)
       do i=1,grad_dim
        do k=1,n_ch
          grad_cart(i_unique)%m(k,:,i_equal) = &
            grad_cart(i_unique)%m(k,:,i_equal) + &
                rotmat(i,:)*grad_totsym(k,index+i)
        enddo

       end do
    end do
 enddo ua1

 if(N_moving_unique_timps.ne.0) then
    do i_unique=1,N_moving_unique_timps
       i_center=moving_unique_timp_index(i_unique)
       index=gradient_index(N_moving_unique_atoms+i_unique)
       grad_dim=gradient_index(N_moving_unique_atoms+i_unique+1)-index
       index=index-1
       do i_equal=1,unique_timps(i_center)%n_equal_atoms
          rotmat=>unique_timp_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
           do k=1,n_ch
             grad_cart(i_unique+N_moving_unique_atoms)%m(k,:,i_equal) = &
                  grad_cart(i_unique+N_moving_unique_atoms)%m(k,:,i_equal) + &
                  rotmat(i,:)*grad_totsym(k,index+i)
           enddo
          end do
       end do
    enddo
 endif

end subroutine cpks_transform_to_cart_grads


subroutine  cpks_add_fit_ch_grads(grad_totsym,grad_fit)
  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components
  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom
 real(kind=r8_kind), intent(inout) :: grad_totsym(:,:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:,:)  ! fitfct. contributions
 integer(kind=i4_kind) :: i_unique,i_center,index,i,grad_dim,k
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
     do k=1,size(grad_totsym,1)
       grad_totsym(k,index) = grad_totsym(k,index) - &     ! sign changed
            weight * sum( rotmat(i,:) * grad_fit(k,:,i_unique) )
     enddo
       index=index+1
    end do
 end do
end subroutine cpks_add_fit_ch_grads

subroutine  add_fit_ch_grads(grad_totsym,grad_fit)

  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components

  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom

 real(kind=r8_kind), intent(inout) :: grad_totsym(:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:)  ! fitfct. contributions
 integer(kind=i4_kind) :: i_unique,i_center,index,i,grad_dim
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
!   print*,'rotmat index i_unique i',index,i_unique,i,rotmat(i,:)
       grad_totsym(index) = grad_totsym(index) - &
            weight * sum( rotmat(i,:) * grad_fit(:,i_unique) )
       index=index+1
    enddo
 enddo
end subroutine add_fit_ch_grads


subroutine  add_xc_grads(grad_cart,grad_xc)
  ! purpose: add xc gradient contributions to a cartesian gradient array
  !    note: xc gradient contributions contain negative gradients -d/dR !
  type(arrmat2), intent(inout) :: grad_cart(:) ! cartesian gradient array
  type(arrmat2), intent(in   ) :: grad_xc(:)   ! xc gradient contribution
 integer(kind=i4_kind) :: i_unique,i_center,i_equal

 do i_unique=1,N_moving_unique_atoms
    i_center = moving_unique_atom_index(i_unique)
    do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       grad_cart(i_unique)%m(:,i_equal) = &
            grad_cart(i_unique)%m(:,i_equal) - &
            grad_xc(i_unique)%m(:,i_equal)
    end do
 end do
end subroutine add_xc_grads

subroutine  add_xc_dervs(dervs_cart,dervs_xc)
  ! purpose: add xc gradient contributions to a cartesian gradient array
  !    note: xc gradient contributions contain negative gradients -d/dR !
  type(arrmat4), intent(inout) :: dervs_cart(:,:) ! cartesian gradient array
  type(arrmat4), intent(in   ) :: dervs_xc(:,:)   ! xc dervs contribution
 integer(kind=i4_kind) :: i_unique,i_center,i_equal
 integer(kind=i4_kind) :: i_unique2,i_center2,i_equal2

 do i_unique=1,N_moving_unique_atoms
    i_center = moving_unique_atom_index(i_unique)
    do i_equal=1,unique_atoms(i_center)%n_equal_atoms
       do i_unique2=1,N_moving_unique_atoms
       i_center2=moving_unique_atom_index(i_unique2)
        do i_equal2=1,unique_atoms(i_center2)%n_equal_atoms
        dervs_cart(i_unique,i_unique2)%m(:,i_equal,:,i_equal2) = &
          dervs_cart(i_unique,i_unique2)%m(:,i_equal,:,i_equal2) - &
            dervs_xc(i_unique,i_unique2)%m(:,i_equal,:,i_equal2)
        enddo
       enddo
    end do
 enddo
end subroutine add_xc_dervs

subroutine  sum_up_and_reset_part_grads()
  ! purpose: sum up the partial cartesian gradient contributions
  !          and re-initialize the contribution array to zero
  integer(kind=i4_kind) :: i_unique,i_center,i_equal

  do i_unique=1,N_moving_unique_atoms
     i_center = moving_unique_atom_index(i_unique)
     do i_equal=1,unique_atoms(i_center)%n_equal_atoms
        gradient_cartesian(i_unique)%m(:,i_equal) = &
        gradient_cartesian(i_unique)%m(:,i_equal) + &
             partial_cartesian(i_unique)%m(:,i_equal)
     end do
     partial_cartesian(i_unique)%m = 0.0_r8_kind
  end do
end subroutine sum_up_and_reset_part_grads

subroutine core_gradient_calc()
  !  Purpose: This routine actually calculates the
  !           gradient of the electrostatic core-core interaction.
  !           The gradient is added to gradient_totalsym
  !** End of interface *****************************************
  !------------ Modules used ----------------------------------
  use pointcharge_module
#ifdef WITH_EPE
  use ewaldpc_module
#endif
  use population_module, only: m_charge
  use energy_calc_module, only: get_energy
#ifdef WITH_EPE
   use energy_calc_module, only: addional_core_ewpc_interaction
#endif
  implicit none
  !------------ Declaration of local variables -----------------
  real(kind=r8_kind)          :: za,zb,dist,wa,dist2
  real(kind=r8_kind)          :: zca,zcb,C,A
  integer(kind=i4_kind)       :: ma,na,nb,eq_b,i, n_equal_charges
! integer(kind=i4_kind)       :: na1,na2
  integer(kind=i4_kind)       :: index_a,index_b,grad_dim_a,grad_dim_b
  real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:),rotmat(:,:)
  real(kind=r8_kind) :: gradient_a(3),derivative_aa(3,3)
! real(kind=r8_kind) :: dervs_core(9,9)=0.0_r8_kind
  real(kind=r8_kind), allocatable,dimension(:,:) :: totalsym_cartesian_aa, totalsym_cartesian_ba
  real(kind=r8_kind), allocatable,dimension(:,:,:) :: totalsym_cartesian_bb,totalsym_cartesian_ab
  integer (i4_kind) :: k
#ifdef WITH_EPE
  real (r8_kind) :: tot_en
  character(len=300) :: inp_dir !!!!!!!!!AS
  logical :: lcomp_ch !!!!!!!!!!!AS
  integer (i4_kind) :: funit
#endif

  !------------ Executable code --------------------------------

!  print*
!  print*,'gradient_totalsym core contrib ONLY !!!!!!!!!'
!  gradient_totalsym=0.0_r8_kind          !!!!!!!
!  dervs_totalsym=0.0_r8_kind               !!!!!!!

!!   this is for check only, to be commented
!   dervs_totalsym=0.0_r8_kind               !!!!!!!
!   do na1=1,N_unique_atoms
!    do na2=1,N_unique_atoms
!     dervs_cartesian(na1,na2)%m=0.0_r8_kind
!    enddo
!   enddo
!      call transform_to_cart_dervs(dervs_cartesian,dervs_totalsym)        !!!!!
!    print*
!!------------------------------------------

#ifdef WITH_EPE
  if(ewpc_n > 0 .and. nzepe > 0) then
     call getenv("TTFSINPUTDIR",inp_dir)
     inquire(file=trim(inp_dir)//'/compl_charges',exist=lcomp_ch)
     funit=openget_iounit(trim(inp_dir)//'/compl_charges',  &
          form='formatted', status='unknown')
     if(lcomp_ch) then
        do i=1,N_unique_atoms
           read(funit,*) m_charge(i)
        enddo
     else
        do i=1,N_unique_atoms
           m_charge(i)=zepe(i)-m_charge(i)
           if(zepe(i) == 0.0_r8_kind) m_charge(i)=0.0_r8_kind
           write(funit,'(f16.12)') m_charge(i)
        enddo
     endif
     call returnclose_iounit(funit)
  endif
#endif

  unique1: do ma=1,N_moving_unique_atoms+n_moving_unique_timps
     if(ma.gt.N_moving_unique_atoms) then
        na = moving_unique_timp_index(ma-N_moving_unique_atoms)
        wa = unique_timps(na)%n_equal_atoms
        za = unique_timps(na)%Z
        xa => unique_timps(na)%position
        zca = unique_timps(na)%ZC
     else

        na = moving_unique_atom_index(ma)
        wa = unique_atoms(na)%n_equal_atoms
        za = unique_atoms(na)%Z
        xa => unique_atoms(na)%position
        zca = unique_atoms(na)%ZC
     end if

     grad_dim_a=gradient_index(ma+1)-gradient_index(ma)
     index_a=gradient_index(ma)-1

     if(integralpar_2dervs) then
      allocate(totalsym_cartesian_aa(grad_dim_a,3),stat=cpksalloc(110))
      ASSERT(cpksalloc(110).eq.0)
     endif

     if(ma.gt.N_moving_unique_atoms) then
        rotmat=>unique_timp_grad_info(ma-N_moving_unique_atoms)%m(:,:,1)
     else
        rotmat=>unique_atom_grad_info(ma)%m(:,:,1)
     end if


     if (zca/=0.0_r8_kind.and.(.not.operations_core_density))za=za-zca
     gradient_a=0.0_r8_kind

!    if(integralpar_2dervs) derivative_aa=0.0_r8_kind
                            ! each derivative_aa contrib is remaped
                            ! to  dervs totalsym at the and of  unique2 loop

!   << gradient = - d/dRa E_nuc = Sum(Rb/=Ra) Za Zb (Ra - Rb) / |Ra - Rb|^3 >>
     unique2: do nb=1,N_unique_atoms

#if 1
     if(integralpar_2dervs) then
!!     dervs_totalsym=0.0_r8_kind !!!   this is only to see individual components
!!                                !!!   to be commented
        derivative_aa=0.0_r8_kind
                            ! each derivative_aa contrib is remaped
                            ! to  dervs totalsym at the and of  unique2 loop
        grad_dim_b=gradient_index(nb+1)-gradient_index(nb)
        index_b=gradient_index(nb)-1
      allocate(totalsym_cartesian_bb(grad_dim_b,3,unique_atoms(nb)%n_equal_atoms), &
               totalsym_cartesian_ba(grad_dim_b,3), &
               totalsym_cartesian_ab(grad_dim_a,3,unique_atoms(nb)%n_equal_atoms), &
               stat=cpksalloc(111))
      ASSERT(cpksalloc(111).eq.0)

      totalsym_cartesian_bb=0.0_r8_kind
      totalsym_cartesian_ba=0.0_r8_kind
      totalsym_cartesian_ab=0.0_r8_kind

      totalsym_cartesian_aa=0.0_r8_kind

     endif
#endif

        zb = unique_atoms(nb)%Z
        xb => unique_atoms(nb)%position
        zcb = unique_atoms(nb)%ZC
        if (zcb/=0.0_r8_kind.and.(.not.operations_core_density))zb=zb-zcb

        equal2: do eq_b=1,unique_atoms(nb)%N_equal_atoms
        derivative_aa=0.0_r8_kind
                            ! each derivative_aa contrib is remaped
                            ! to  dervs totalsym at the and of  unique2 loop

           if(ma.le.N_moving_unique_atoms .and. &
                (na.eq.nb).and.(1.eq.eq_b)) then
              cycle equal2
           endif

           dist = sqrt(sum((xa(:,1)-xb(:,eq_b))**2))
           gradient_a = gradient_a + za*zb/dist**3*(xa(:,1)-xb(:,eq_b))

  if(integralpar_2dervs) then
           do k=1,3
            derivative_aa(k,k)=derivative_aa(k,k)+za*zb/dist**3
           enddo
            derivative_aa = derivative_aa - &
           spread(3.0_r8_kind*za*zb/dist**5*(xa(:,1)-xb(:,eq_b)),1,3)* &
           spread((xa(:,1)-xb(:,eq_b)),2,3)

!     print*,'derivative_aa',derivative_aa(1,:)
!     print*,'derivative_aa',derivative_aa(2,:)
!     print*,'derivative_aa',derivative_aa(3,:)

     do i=1,grad_dim_a
        totalsym_cartesian_aa(i,:) =  &
            totalsym_cartesian_aa(i,:) + 0.5_r8_kind* &
             wa * ( rotmat(i,1) * derivative_aa(1,:)   &
                   +rotmat(i,2) * derivative_aa(2,:)   &
                   +rotmat(i,3) * derivative_aa(3,:)   )
        totalsym_cartesian_ab(i,:,eq_b)= &
            totalsym_cartesian_ab(i,:,eq_b) + 0.5_r8_kind* &
             wa * ( rotmat(i,1) * derivative_aa(1,:)   &
                   +rotmat(i,2) * derivative_aa(2,:)   &
                   +rotmat(i,3) * derivative_aa(3,:)   )
     enddo

     do i=1,grad_dim_b
        totalsym_cartesian_bb(i,:,eq_b) =  &
            totalsym_cartesian_bb(i,:,eq_b) + 0.5_r8_kind* &
             wa * ( unique_atom_grad_info(nb)%m(i,1,eq_b) * derivative_aa(1,:)   &
                   +unique_atom_grad_info(nb)%m(i,2,eq_b) * derivative_aa(2,:)   &
                   +unique_atom_grad_info(nb)%m(i,3,eq_b) * derivative_aa(3,:)   )
        totalsym_cartesian_ba(i,:) =  &
            totalsym_cartesian_ba(i,:) +  0.5_r8_kind* &
             wa * ( unique_atom_grad_info(nb)%m(i,1,eq_b) * derivative_aa(1,:)   &
                   +unique_atom_grad_info(nb)%m(i,2,eq_b) * derivative_aa(2,:)   &
                   +unique_atom_grad_info(nb)%m(i,3,eq_b) * derivative_aa(3,:)   )
     enddo
!  print*,ma,nb,eq_b,'core contrib'
!  print*,derivative_aa(1,:),-derivative_aa(1,:)
!  print*,derivative_aa(2,:),-derivative_aa(2,:)
!  print*,derivative_aa(3,:),-derivative_aa(3,:)

!  print*,-derivative_aa(1,:),derivative_aa(1,:)
!  print*,-derivative_aa(2,:),derivative_aa(2,:)
!  print*,-derivative_aa(3,:),derivative_aa(3,:)


   endif

  enddo equal2

!   print*,'uniques',ma,nb

#if 1
  dervs: if(integralpar_2dervs) then
   if(.true.) then

!    print*,index_a+1,index_a+grad_dim_a,index_a+1,index_a+grad_dim_a

     do i=1,grad_dim_a
                !     print*,'totalsym_cartesian_a',totalsym_cartesian_aa(i,:),i
     do k=1,3
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i) = &
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i) - &   !(1)
     unique_atom_grad_info(ma)%m(i,k,1)* &
                        totalsym_cartesian_aa(:,k)
          !   this is aa dia block thus - sign here
     enddo
              !   print*,dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i)
     enddo

              !    print*,index_a+1,index_a+grad_dim_a,index_b+1,index_b+grad_dim_b
     do i=1,grad_dim_b
     do eq_b=1,unique_atoms(nb)%N_equal_atoms
     do k=1,3
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_b+i) = &
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_b+i) + &   !(2)
     unique_atom_grad_info(nb)%m(i,k,eq_b)* &
                        totalsym_cartesian_ab(:,k,eq_b)
            ! thid is off dia block, + sign here

     enddo
     enddo
              !     print*,dervs_totalsym(index_a+1:index_a+grad_dim_a,index_b+i)
     enddo

     do i=1,grad_dim_a
     do k=1,3
     dervs_totalsym(index_b+1:index_b+grad_dim_b,index_a+i) = &
     dervs_totalsym(index_b+1:index_b+grad_dim_b,index_a+i) + &   !(3)
     unique_atom_grad_info(ma)%m(i,k,1)* &
                        totalsym_cartesian_ba(:,k)
     enddo
     enddo

     do i=1,grad_dim_b
     do eq_b=1,unique_atoms(nb)%N_equal_atoms
     do k=1,3
     dervs_totalsym(index_b+1:index_b+grad_dim_b,index_b+i) = &
     dervs_totalsym(index_b+1:index_b+grad_dim_b,index_b+i) - &   !(4)
     unique_atom_grad_info(nb)%m(i,k,eq_b)* &
                        totalsym_cartesian_bb(:,k,eq_b)
     enddo
     enddo
     enddo
   endif

      deallocate(totalsym_cartesian_bb,totalsym_cartesian_ab, &
                 totalsym_cartesian_ba,stat=cpksalloc(111))
      ASSERT(cpksalloc(111).eq.0)
      cpksalloc(111)=1

  endif dervs
#endif

!   print*,' core dervs ma, nb', ma,nb
!   do i=1,size(dervs_totalsym,1)
!    write(*,'(10f8.4)') dervs_totalsym(i,:),sum(dervs_totalsym(i,:))
!   enddo
!    write(*,'(10f8.4)') sum(dervs_totalsym,1)

!!   this is for check only, to be commented
!   do na1=1,N_unique_atoms
!    do na2=1,N_unique_atoms
!     dervs_cartesian(na1,na2)%m=0.0_r8_kind
!    enddo
!   enddo
!      call transform_to_cart_dervs(dervs_cartesian,dervs_totalsym)        !!!!!
!    print*
!!------------------------------------------

  enddo unique2

!     print*,'gradient core',gradient,ma
!     print*, derivative(:,1)
!     print*, derivative(:,2)
!     print*, derivative(:,3)

!   << gradient += - d/dRa E_PC = Sum(Rc) Za Zc (Ra - Rc) / |Ra-Rc|^3 >>
     if(integralpar_2dervs) derivative_aa=0.0_r8_kind
     pointcharge: do nb=1,pointcharge_N + n_timps
        if(nb<=n_timps) then
           zb = unique_timps(nb)%Z - unique_timps(nb)%ZC
           xb => unique_timps(nb)%position
           n_equal_charges = unique_timps(nb)%n_equal_atoms
        else
           zb = pointcharge_array(nb-n_timps)%Z
           xb => pointcharge_array(nb-n_timps)%position
           n_equal_charges = pointcharge_array(nb-n_timps)%n_equal_charges
           C=pointcharge_array(nb-n_timps)%C
           A=pointcharge_array(nb-n_timps)%A
        end if
        equal_pc: do eq_b=1,n_equal_charges
           if(ma.gt.N_moving_unique_atoms &
                .and.na.eq.nb.and.eq_b.eq.1) cycle equal_pc
           dist2 = sum((xa(:,1)-xb(:,eq_b))**2)
           dist = sqrt(dist2)
           if(dist <= 1.0e-10_r8_kind) cycle equal_pc
           gradient_a = gradient_a + za*zb/dist**3*(xa(:,1)-xb(:,eq_b))
           if(C /= 0.0_r8_kind) then
              gradient_a = gradient_a - &
                   C*exp(-A*dist2)*za*zb*(1.0+2.0_r8_kind*A*dist2)*(xa(:,1)-xb(:,eq_b))/dist**3
           end if
           if(integralpar_2dervs) then
             do k=1,3
              derivative_aa(k,k)=derivative_aa(k,k)+za*zb/dist**3
             enddo
              derivative_aa = derivative_aa - &
             spread(3.0_r8_kind*za*zb/dist**5*(xa(:,1)-xb(:,eq_b)),1,3)* &
             spread((xa(:,1)-xb(:,eq_b)),2,3)
           endif
        enddo equal_pc
     enddo pointcharge

#ifdef WITH_EPE
     if(EWPC_N.ne.0) then
!!!        gradient_a=0.0_r8_kind   !!!!!!!!!!!!!!!! degug only stat
        if(nzepe > 0) za=za+m_charge(ma) !!!!!!!!!!!!!AS
        ewpc: do nb=1,EWPC_N
           zb = ewpc_array(nb)%Z
           xb => ewpc_array(nb)%position
           equal_ewpc: do eq_b=1,ewpc_array(nb)%n_equal_charges
              dist = sqrt(sum((xa(:,1)-xb(:,eq_b))**2))
              gradient_a = gradient_a + za*zb/dist**3*(xa(:,1)-xb(:,eq_b))
           if(integralpar_2dervs) then
             do k=1,3
              derivative_aa(k,k)=derivative_aa(k,k)+za*zb/dist**3
             enddo
              derivative_aa = derivative_aa - &
             spread(3.0_r8_kind*za*zb/dist**5*(xa(:,1)-xb(:,eq_b)),1,3)* &
             spread((xa(:,1)-xb(:,eq_b)),2,3)
           endif
           enddo equal_ewpc
        enddo ewpc
     endif ! EWPC_N.ne.0
#endif


!   << gradient_totalsym += d/dRa E_nuc + E_PC >>

     do i=1,grad_dim_a
        gradient_totalsym(index_a+i) = gradient_totalsym(index_a+i) - &
                                  wa * sum( rotmat(i,:) * gradient_a(:) )

!  print*, gradient_totalsym(index_a+i),' nuc gradient_totalsym'
!  print*, dervs_totalsym(index_a+i,:)
     enddo

    if(integralpar_2dervs) then
#ifdef WITH_EPE
    if(EWPC_N.ne.0 .or. pointcharge_N+n_timps.ne.0) then

     totalsym_cartesian_aa=0.0_r8_kind
     do i=1,grad_dim_a
        totalsym_cartesian_aa(i,:) =  &
            totalsym_cartesian_aa(i,:) + &! 0.5_r8_kind* &
             wa * ( rotmat(i,1) * derivative_aa(1,:)   &
                   +rotmat(i,2) * derivative_aa(2,:)   &
                   +rotmat(i,3) * derivative_aa(3,:)   )
     enddo
     do i=1,grad_dim_a
                !     print*,'totalsym_cartesian_a',totalsym_cartesian_aa(i,:),i
     do k=1,3
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i) = &
     dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i) - &   !(1)
     unique_atom_grad_info(ma)%m(i,k,1)* &
                        totalsym_cartesian_aa(:,k)
     enddo
!                 print*,dervs_totalsym(index_a+1:index_a+grad_dim_a,index_a+i)*2.0_r8_kind
     enddo
    endif
#endif
      deallocate(totalsym_cartesian_aa,stat=cpksalloc(110))
      ASSERT(cpksalloc(110).eq.0)
      cpksalloc(110)=1
    endif

  enddo unique1

#ifdef WITH_EPE
  if (ewpc_n > 0 .and. nzepe > 0) then
     call addional_core_ewpc_interaction(m_charge) !!!!!!!!!!!!!AS
     call get_energy(tot=tot_en)
     call write_to_trace_unit("MAIN_GRADIENT: corrected final total energy (core-ewpc): ",real=tot_en)
     print *, 'MAIN_GRADIENT: corrected final total energy (core-ewpc): ',tot_en
  endif

  if(allocated(zepe)) deallocate(zepe) !!!!!!!!!!!!!!!AS
#endif

  ! FIXME: WITH_EPE?
  if(allocated(m_charge)) deallocate(m_charge) !!!!!!!!!!!!!!!AS

  nullify(xa)
  nullify(xb)

  call calc_nuc_dqo_grads(gradient_totalsym,gradient_index)

  call calc_nuc_id_grads(gradient_totalsym,gradient_index)

  if(moving_pc) call calc_PC_grads()

  if(moving_X_centers) call calc_X_grads()

  if(moving_Pol_centers) call calc_id_grads()

!     if(integralpar_2dervs) then                                             !!!!!
!   do na=1,N_unique_atoms
!    do nb=1,N_unique_atoms
!     dervs_cartesian(na,nb)%m=0.0_r8_kind
!    enddo
!   enddo
!      call transform_to_cart_dervs(dervs_cartesian,dervs_totalsym)        !!!!!
!    print*
! endif

end subroutine core_gradient_calc

subroutine  add_solv_grads
  integer(kind=i4_kind) :: i_unique,i_center,i_equal

  do i_unique=1,N_moving_unique_atoms
     i_center = moving_unique_atom_index(i_unique)
     do i_equal=1,unique_atoms(i_center)%n_equal_atoms
        gradient_cartesian(i_unique)%m(:,i_equal) = &
             gradient_cartesian(i_unique)%m(:,i_equal) + &
             grad_solv_cart(i_unique)%m(:,i_equal)
     end do
  end do

end subroutine add_solv_grads

  subroutine impl_fit_dervs(dervs)
    !
    ! seems to compute fit contribution to second derivatives:
    !
    ! dervs(X,Y) += SUM(kl) aY(k) * GX(k,l) * a(l)
    !
    ! input: cpks_fitcoeff_grads calculated in cpks_g4_constructs
    use cpksdervs_matrices
    implicit none
    real(r8_kind), intent(inout) :: dervs(:,:)
    ! *** end of interface ***

    integer(i4_kind) :: x,y
    real(r8_kind)    :: d(size(dervs,1),size(dervs,2))

    do x=1,size(dervs,1)
      do y=1,size(dervs,2)
        d(x,y)= dot_product( cpks_fitcoeff_grads(:,y),cpks_grad_fit_totasym(:,x))
       enddo
    enddo

    ! OUTPUT: add my contribution:
    dervs = dervs + d

    deallocate(cpks_fitcoeff_grads,stat=cpksalloc(108))
    ASSERT(cpksalloc(108).eq.0)
    cpksalloc(108)=1
  end subroutine impl_fit_dervs

#if WITH_ERI4C == 1
  subroutine exact_2e_grads( densmat,                        gradient_totalsym )
    !
    ! Contains calling sequence for eri4c gradient wrapper
    !
    use datatype,             only : arrmat3
    use eri4c_options,        only : J_exact                                   &
                                   , K_exact                                   &
                                   , QQP_tol                                   &
                                   , EXX_factor
    use eri_main,             only : erg_blk_main
    use symm_positions,       only : symm_mapping_shells
    use group_module,         only : group_num_el
    use symmetry_data_module, only : symmetry_data_n_spin                      &
                                   , symmetry_data_n_irreps
    use comm,                 only : comm_bcast                                &
                                   , comm_rank
    !
    type(arrmat3),                 intent(in)    :: densmat(:) ! dmat(irr)%m   !
    real(r8_kind),                 intent(inout) :: gradient_totalsym(:)
    !
    real(r8_kind), allocatable                   :: jexact_grad(:)
    real(r8_kind), allocatable                   :: kexact_grad(:)
    integer(i4_kind)                             :: n_spin
    integer(i4_kind)                             :: n_atm, n_sym
    integer(i4_kind)                             :: n_shl, n_irr
    integer(i4_kind)                             :: n_ua1
    !
    !--------------------------------------------------------------------------+
    !
    n_atm  = sum( unique_atoms(:)%n_equal_atoms )
    n_ua1  = N_unique_atoms + 1
    n_shl  = sum( unique_atoms(:)%n_equal_atoms                                &
                * ( unique_atoms(:)%lmax_ob + 1 ) )
    n_irr  = symmetry_data_n_irreps()
    n_spin = symmetry_data_n_spin()
if(gradient_data_n_gradients/=size(gradient_totalsym))stop 'gradient_data_n_gradients/=size(gradient_totalsym'
    !
    ! GET SYMMETRY MAPPING INDICES FOR INDIVIDUAL SHELLS
    n_sym = group_num_el
    call comm_bcast( n_sym )
    !
    allocate( kexact_grad(gradient_data_n_gradients) )
    kexact_grad = 0.0_r8_kind
    IF ( K_exact .and. J_exact ) THEN
      !
      stop 'j-exact gradients not supported'
      !
    ELSEIF ( K_exact .and. .not. J_exact ) THEN
      !
      call erg_blk_main( QQP_tol, n_sym, n_atm, n_ua1                          &
                       , gradient_data_n_gradients                             &
                       , n_spin, n_shl, n_irr, -1,                 kexact_grad &
                       , unique_atom_grad_info, densmat                        &
                       , gradient_index, unique_atoms                          )
      !
    ELSEIF ( J_exact .and. .not. K_exact ) THEN
      !
      stop 'j-exact gradients not supported'
      !
    ENDIF
    !
    IF( comm_rank() == 0 ) gradient_totalsym = gradient_totalsym - kexact_grad * EXX_factor
    !
  end subroutine exact_2e_grads
#endif

end subroutine main_gradient

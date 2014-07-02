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
subroutine  integral_calc_quad_2cob3c()
!----------------------------------------------------------------
!
!  Purpose: This routine is the main routine for the
!           integral calculation of one quadrupel
!           ( unique atom 1, l 1, unique atom 2, l 2 )
!           for 2 center orbital and three center integral
!           calculation. The corresponding data (including
!           the quadrupel to be calculated) are stored in
!           int_data_2cob3c_module.
!           Contraction and symmetry-adaption are directly
!           done within this subroutine.
!
!  Subroutine called by: main_slave integral_main_2cob3c
!
!  References: Publisher Document: Concepts of Integral Part
!
!  Author: TB
!  Date: 6/96
!
  !===================================================================
  ! End of public interface of module
  !===================================================================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: MS
! Date:   3/97
! Description: Necessary changes and supplements for calculation of
! energy gradients
!
! Modification (Please copy before editing)
! Author: MS
! Date:  7/97
! Description: Changes for  relativistic and relativistic gradients have
!        been added
!
! Modification (Please copy before editing)
! Author: MS
! Date:   8/97
! Description: The symmetryadaption was extended and is now also
!              able to treat pseudo 2D irreps
!
! Modification (Please copy before editing)
! Author: Uwe Birkenheuer
! Date:   June 1998
! Description: Moving_Unique_Atom concept introduced
!              Split_Gradients concept introduced
!              Gradients for Model_Density_Approach introduced
!
! Modification (Please copy before editing)
! Author: MM
! Date:   8/97
! Description: Symmetryadaption of spin orbit matrix elements
!
! Modification (Please copy before editing)
! Author: AS
! Date:   11-12/99
! Description: integrals of electrostatic potential are added
!
! Modification (Please copy before editing)
! Modification
! Author: AM
! Date:   first: 02.03.1999
! Description: ...
!
! Modification use DLB instead of master/slave for 2cop3c integrals
! Author: AN
! Date:   4/11
! Description: use DLB for scheduling batches of 2cop3c integrals
!              remove reporting back of slaves, they get their new
!              tasks via DLB
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!----------------------------------------------------------------
!..............................................................................
! F1 = Sum(n,s) < d/dR psi_n,s | T+V_nuc+V_H       - eps_n,s | psi_n,s >
!      or
!      Sum(n,s) < d/dR psi_n,s | T+V_nuc+V_H+V_X,s - eps_n,s | psi_n,s >
!
! F2 = < rho_tot | d/dR V_nuc >
!
! F3 = < rho_tot | d/dR V_H >
!
! F4 = Sum(s) < rho_s | d/dR V_X,s >
!
!
! Storage
! ~~~~~~~
! options_split_gradients      | FALSE     FALSE        TRUE   TRUE
! model_density                ! FALSE     TRUE         FALSE  TRUE
! -----------------------------+------------------------------------
! gradient_totalsym(1:n_grads) | F1+F2+F3  F1+F2+F3+F4  F2     F2
! gradient_ob_pulay(1:n_grads) | --        --           F1     F1
! gradient_ch_pulay(1:n_grads) | --        --           F3     F3
! gradient_mda_vxc (1:n_grads) | --        --           --     F4
!..............................................................................

!.......................... STORAGE FOR GRADIENTS .............................
!  this a copy of comments from LL_CALCULATE_GRADS:
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

  !------------ Modules used ------------------------------------
! use CPU_TIME:
! define FPP_TIMERS 2
# include "def.h"
! define _NUSD_BY_CALC_3CENTER
  use type_module ! type specification parameters
  use datatype
  use strings, only: itoa
  use timer_module
  use time_module
  use integralpar_module,        only: integralpar_gradients                   &
                                     , integralpar_dervs                       &
                                     , integralpar_pot_for_secderiv            &
                                     , integralpar_totalsymmetric              &
                                     , integralpar_3cob_grad                   &
                                     , integralpar_3c_co                       &
                                     , integralpar_3c_co_resp                  &
                                     , integralpar_3c_r2_pvsp                  &
                                     , integralpar_3c_r2_pvxp                  &
                                     , integralpar_3c_xc                       &
                                     , integralpar_3c_rcoul_pvsp               &
                                     , integralpar_3c_rcoul_pvxp               &
                                     , integralpar_send_3c                     &
                                     , integralpar_i_int_part                  &
                                     , integralpar_relativistic                &
                                     , integralpar_2cob_pvec                   &
                                     , integralpar_2cob_potential              &
                                     , integralpar_2cob_field                  &
                                     , integralpar_solv_grad                   &
                                     , integralpar_Q_solv_grad                 &
                                     , integralpar_2cob_pc_grad                &
                                     , integralpar_2cob_X_grad                 &
                                     , integralpar_2cob_ipd_grad               &
                                     , integralpar_2cob_kin                    &
                                     , integralpar_2cob_kin_grad               &
                                     , integralpar_2cob_nuc                    &
                                     , integralpar_2cob_nuc_grad               &
                                     , integralpar_2cob_pvsp                   &
                                     , integralpar_2cob_pvsp_grad              &
                                     , integralpar_2cob_pvxp                   &
                                     , integralpar_2cob_ol                     &
                                     , integralpar_2cob_ol_grad                &
                                     , integralpar_cpks_contribs               &
                                     , integralpar_cpksdervs                   &
                                     , integralpar_pseudo
#ifdef WITH_EFP
  use integralpar_module,        only: integralpar_2cob_pc                     &
                                     , integralpar_2cob_X                      &
                                     , integralpar_2cob_field_efp              &
                                     , integralpar_efp_gradients
#endif
  use int_data_2cob3c_module
  use integral_calc_quad_module, only: &
      integral_calc_quad_densmat, quad_pmat, quad_wmat  &
    , integral_calc_quad_close &
    , IDONT_CONTRACT, IDONT_SUM_OVER_PARTNERS &
    , contrsym2c, contrsym3c, contr3c &
    , quad_density_mat
  use output_module
  use int_send_2cob3c_module, only : int_send_2cob3c_send
  use int_send_2cob3c_spor_module, only : int_send_2cob3c_spor_send
  use iounitadmin_module
  use gradient_data_module,      only: gradient_totalsym                       &
                                     , gradient_data_n_gradients               &
                                     , gradient_data_n_spin_gradients          &
                                     , dervs_totalsym                          &
                                     , gradient_ob_pulay                       &
                                     , gradient_ch_pulay                       &
                                     , gradient_mda_vxc                        &
                                     , gradient_index
#ifdef WITH_EPE
  use gradient_data_module, only: calc_cluster_epe_energy, cluster_epe_energy
#endif
  use unique_atom_module,        only: unique_atom_symequivatoms_type          &
                                     , unique_atom_partner_type                &
                                     , unique_atom_sa_int_type                 &
                                     , unique_atom_symequiv                    &
                                     , unique_atoms                            &
                                     , N_unique_atoms                          &
                                     , pseudopot_present
  use symmetry_data_module
  use comm_module
  use msgtag_module
  use options_module, only: options_split_gradients, options_n_spin, &
                            xcmode_model_density, xcmode_extended_mda, &
                            options_xcmode, options_spin_orbit
#ifdef WITH_EPE
  use ewaldpc_module,only: ewpc_n
#endif
  use pointcharge_module, only: totsym_grad_pc_length
  use point_dqo_module, only: totsym_grad_dip_length, &
       totsym_grad_quad_length,totsym_grad_oct_length, totsym_grad_rep_length
  use induced_dipoles_module, only: totsym_grad_idip_length
  use elec_static_field_module, only: totsym_field_length
  use solv_cavity_module, only: with_pc, fixed_pc
  use spin_orbit_module, only: is_on,op_FitTrafo
  use error_module,      only: MyID !<<< debug
  use calc3c_switches!!!,only: integralpar_cpksdervs,integralpar_2dervs,integralpar_dervs
  USE_MEMLOG
  use efield_module, only: efield_applied

#ifdef WITH_RESPONSE
  use int_send_3c_resp, only: int_send_3c_resp_save
#endif

  use interfaces ! FIXME: need it here?

#ifdef WITH_DFTPU
  use dft_plus_u_module, only: dft_plus_u_grad, dft_plus_u_in_use &! for DFT+U gradients
                             , dft_plus_u_mo_in_use
#endif
  implicit none

  !
  ! interfaces moved to
  !     modules/interfaces.f90
  !

  !------------ Declaration of local variables ------------------
  type(unique_atom_symequivatoms_type), pointer :: symequivatoms
    ! different sym.equiv. atom pairs to be processed
  integer(kind=i4_kind)  :: i_symequiv
    ! index of not-symmetry-equivalent connection vector
  real(r8_kind)          :: weight
    ! weight of not-symmetry-equivalent connection vector
  integer(kind=i4_kind)  :: i_ua1, i_ua2, i_l1, i_l2
    ! local vars for for quadrupel%ua1, quadrupel%ua2, ...
  integer(kind=i4_kind)  :: i_ea1, i_ea2
    ! loop indices for equal atoms
  logical :: split_gradients, model_density

  real(r8_kind), parameter :: zero = 0.0_r8_kind

  logical                  :: &
       diagonal_ua, & ! if ( ua1 == ua2 )
       diagonal_ea, & ! if ( ua1 == ua2 && ea1 = ea2 )
       diagonal_sh    ! if ( ua1 == ua2 && ea1 = ea2 && l1 == l2 )

  logical :: do_sum_up
#ifdef WITH_EFP
  logical :: do_sum_up_efp
#endif

#ifdef FPP_TIMERS
  real(r8_kind), parameter :: interval=300.0
  real(r8_kind)            :: last_time=-2*interval
  FPP_TIMER_DECL(i2c3c)
  FPP_TIMER_DECL(oli)
  FPP_TIMER_DECL(olg)
  FPP_TIMER_DECL(shgi)
  FPP_TIMER_DECL(cntr3)
  FPP_TIMER_DECL(sym3c)
  FPP_TIMER_DECL(sumup)
  FPP_TIMER_DECL(tomo)
#endif

 !--------------------------------------------------------------
 !------------ Executable code ---------------------------------

  FPP_TIMER_START(i2c3c)

  i_ua1 = quadrupel%ua1
  i_l1  = quadrupel%l1
  i_ua2 = quadrupel%ua2
  i_l2  = quadrupel%l2

  call dtrace("start with quadrupel "//itoa(i_ua1)//" " &
                                     //itoa(i_l1) //" " &
                                     //itoa(i_ua2)//" " &
                                     //itoa(i_l2)       &
                                     )
  call dtrace("setup")
  call setup()
  split_gradients = options_split_gradients()
  model_density = options_xcmode() == xcmode_model_density .or. &
                  options_xcmode() == xcmode_extended_mda

  ! decide whether or not this is a "property" run that will
  ! require taking traces with density matrix (see sum_up_gradient)
  do_sum_up = integralpar_gradients .or. integralpar_pot_for_secderiv
#ifdef WITH_EFP
  do_sum_up = do_sum_up .or. integralpar_2cob_field_efp .or. &
       integralpar_2cob_pc .or. integralpar_2cob_X .or. integralpar_efp_gradients
  do_sum_up_efp = integralpar_2cob_field_efp.or. &
       integralpar_2cob_pc .or. integralpar_2cob_X .or. integralpar_efp_gradients
#endif

  if( do_sum_up )then
    ! if this is a property run, then prepare the section of the
    ! densmat (and energy-weighted densmat) corresponding to the current
    ! quadrupel:
    call integral_calc_quad_densmat(i_ua1, i_l1, i_ua2, i_l2, quad_pmat, quad_wmat)
  endif

  ! allocate primitives now here, since same for all ea's
  ! call allocate_primitives()
  ! MOVED TO int_data_2cob3c_setup()
  !
  ! FIXME: do the same with contracted
  ! currently cont_* are allocated in
  ! contract_2center() and contract_3center()
  ! for each ea


  totalsymmetric: if (integralpar_totalsymmetric) then
     ! In case only total symmertic integrals are calculated,
     ! integrals are only necessary between selected equal atoms

     symequivatoms => unique_atom_symequiv(quadrupel%ua1,quadrupel%ua2)

     ! loop over not symmetry equivalent pairs of equal atoms
     symequiv: do i_symequiv = 1, symequivatoms%n_equiv_dist

        call dtrace("processing pair of non-symmetryequivalent atoms: "//itoa(i_symequiv))

        i_ea1 = 1
        i_ea2 = symequivatoms%index_eq_atom2(i_symequiv)
        center1 = ua1%position(:,i_ea1)
        center2 = ua2%position(:,i_ea2)

        ! weight of a pair (number of symm-equiv distances):
        weight = symequivatoms%weight(i_symequiv)

        diagonal_ua = ( i_ua1 == i_ua2 )
        diagonal_ea = ( diagonal_ua .and. (i_ea1 == i_ea2) )
        diagonal_sh = ( diagonal_ea .and. (i_l1  == i_l2 ) )

        call dtrace("index of secound equal atom: "//itoa(i_ea2) )

        ! allocation used to have been done inside of
        ! call calc_primitives_and_contract()
        ! it was moved outside the loop over ea-ea distances
        ! some ?buggy? ??_calculate_*() progs reguire re-initialization
        ! of primitives to zero
        ! FIXME: better fix the actual subs not to depend on input state
        call allocate_primitives('initialize')

        ! calculate primitive integrals and contract them in non_relativistic
        ! case or else map contracted integrals to primitives
        call dtrace("calc_primitives")

#ifdef no_cpks_coul_grads
    if(allocated(prim_int_dens_mat))then
      ! it is allocated in the (first?) gradient run
      ! for computing (explicit?) gradients due to fit-functions
      ! See int_data_2cob3c_module.

      ! only in first gradient run of the sec-der calculations:
      !   Compute the section of the density matrix in (uncontracted?) basis
      !   that corresponds to the curent quadrupel.
      call quad_density_mat(i_ua1,i_ea1,i_l1,i_ua2,i_ea2,i_l2,weight,prim_int_dens_mat)
    endif
#endif
        call calc_primitives(i_ua1,i_ea1,i_l1,i_ua2,i_ea2,i_l2)

        ! add contribution of one not symmetry equivalent pair
        ! of equal atom to symmetry adapted integrals
        call dtrace("adding results to symmetry adapted arrays(a)")

        call contract_and_symadapt(i_ua1,i_ea1,i_l1,i_ua2,i_ea2,i_l2,0) ! sum over partners totalsym

        call dtrace("1 pair of non-symmetryequivalent atoms done ")

     enddo symequiv

  else ! now not total symmetric

     eqal_atom_1: do i_ea1 = 1, ua1%N_equal_atoms
        eqal_atom_2: do i_ea2 = 1, ua2%N_equal_atoms

           call dtrace("processing pair of equal atoms "//itoa(i_ea1)//" and "//itoa(i_ea2))

           center1 = ua1%position(:,i_ea1)
           center2 = ua2%position(:,i_ea2)

           ! weight of a pair (all ea1-ea2 distances):
           weight = 1.0_r8_kind


           call allocate_primitives('initialize')

           ! calculate primitive integrals and contract them in non_relativistic
           ! case or else map contracted integrals to primitives
           call dtrace("calc_primitives_and_contract")

           call calc_primitives(i_ua1,i_ea1,i_l1,i_ua2,i_ea2,i_l2)


           ! add contribution of one pair
           ! of equal atom to symmetry adapted integrals
           call dtrace("adding results to symmetry adapted arrays(b)")

           ! not a totally symmetric case:
           call contract_and_symadapt(i_ua1, i_ea1, i_l1, i_ua2, i_ea2, i_l2, IDONT_SUM_OVER_PARTNERS)

#ifdef WITH_RESPONSE
           if (integralpar_3c_co_resp) then
              call symadapt_add_nottotalsym_co(i_l1, i_l2, i_ea1, i_ea2, contr3c(prim_int_3c_co))
           end if
#endif

           call dtrace("1 pair of equal atoms done")

        end do eqal_atom_2
     end do eqal_atom_1


  endif totalsymmetric !/else

  ! FIXME: also in not-totalsymmetric case?:
  if(.not.options_spin_orbit)then
     ! in the case of groups with 2D pseudo irreps we have to make
     ! an adittional symmetry adaption

     call symadp_pseudo2D() !!! the only call
  endif

  call dtrace("calculation done")

  if(do_sum_up) then
#ifdef WITH_EFP
     if(.not. do_sum_up_efp) then
#endif
        call dtrace("sum_up_gradient")
        call sum_up_gradient()

#ifdef WITH_DFTPU
        !*************************************************************************************
        !                  DFT+U Gradient subroutine is called here
        !*************************************************************************************
        ! DFT+U Gradient is calculated only when both "integralpar_3cob_grad" and
        ! "dft_plus_u_in_use" are "t". When only "dft_plus_u_in_use" was used, gradient
        ! calculations with both DFT+U and SOLVATION terminated with error.
        if( integralpar_3cob_grad .and. (dft_plus_u_in_use .or. dft_plus_u_mo_in_use) )then
          call dft_plus_u_grad(i_ua1, i_l1, i_ua2, i_l2, symadapt_int_2cob_ol_grad, gradient_totalsym)
        end if
#endif

#ifdef WITH_EFP
     else
        call dtrace("sum_up_grads_on_efp")
        call sum_up_grads_on_efp()
     end if
#endif
     ! free the storage in integral_calc_quad_module
     ! that was allocated by integral_calc_quad_densmat(..)
     ! before proceeding to the next quadrupel:
     call integral_calc_quad_close()
  endif

  if ( integralpar_send_3c .and. .not. integralpar_pot_for_secderiv) then
     call start_timer(timer_int_send_2cob3c(integralpar_i_int_part))
     call dtrace("int_send_2cob3c_send")
     if( .not. options_spin_orbit )then
       call int_send_2cob3c_send()
     else
       call int_send_2cob3c_spor_send()
     endif
     call stop_timer(timer_int_send_2cob3c(integralpar_i_int_part))
  endif

  ! store relativistic gradients/derivatives for future transformations:
  if( integralpar_relativistic .and..not. options_spin_orbit )then
    call relstore(i_ua1,i_l1,i_ua2,i_l2)
  endif
#ifdef WITH_RESPONSE
     if ( integralpar_3c_co_resp ) then
        call start_timer(timer_int_send_2cob3c(integralpar_i_int_part))
        call dtrace("int_send_2cob3c_send")
        call int_send_3c_resp_save(quadrupel%ua1, quadrupel%ua2, &
             quadrupel%l1, quadrupel%l2, symadapt_int_3c_co_resp)
        call stop_timer(timer_int_send_2cob3c(integralpar_i_int_part))
     endif
#endif

  call dtrace("shutdown")
  call shutdown()

  call dtrace("end")

  FPP_TIMER_STOP(i2c3c)
#ifdef FPP_TIMERS
  if( FPP_TIMER_VALUE(i2c3c) - last_time > interval )then
    last_time = FPP_TIMER_VALUE(i2c3c)
    print*,MyID,'TIMER: i2c3c tot =',FPP_TIMER_VALUE(i2c3c)
    print*,MyID,'TIMER:  old ints =',FPP_TIMER_VALUE(oli)
    if( integralpar_gradients )then
    print*,MyID,'TIMER:  old grads=',FPP_TIMER_VALUE(olg)
    endif
    print*,MyID,'TIMER: shgi ints =',FPP_TIMER_VALUE(shgi)
    print*,MyID,'TIMER: contr  3c =',FPP_TIMER_VALUE(cntr3)
    print*,MyID,'TIMER: symadp 3c =',FPP_TIMER_VALUE(sym3c)
    if( integralpar_gradients )then
    print*,MyID,'TIMER: sum up    =',FPP_TIMER_VALUE(sumup)
    print*,MyID,'TIMER:  |- to MO =',FPP_TIMER_VALUE(tomo)
    endif
  endif
#endif


 !--------------------------------------------------------------
 !------------ Private Subroutines -----------------------------
contains ! integral_calc_quad_2cob3c

  subroutine dtrace(msg)
    implicit none
    character(len=*), intent(in) :: msg
    ! *** end of interface ***

    if ( output_int_loops ) &
      call write_to_output_units(MyID//"integral_calc_quad_2cob3c: "//msg)
  end subroutine dtrace


  !*************************************************************
  subroutine setup
    !  Purpose: setup routines of various modules
    !------------ Executable code ------------------------------

    call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
    call init_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_quadrupel_2cob3c(integralpar_i_int_part))

    call int_data_2cob3c_setup()

  end subroutine setup


  subroutine shutdown
    !  Purpose: shutdown routines of various modules
    !------------ Executable code ------------------------------

    call int_data_2cob3c_shutdown()

    call stop_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call stop_timer(timer_int_quadrupel_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
    call timer_small_to_large( &
         timer_int_calc_2cob3c(integralpar_i_int_part), &
         timer_int_calcsum_2cob3c(integralpar_i_int_part) )

  end subroutine shutdown
  !**************************************************************


  !**************************************************************
  !
  ! allocate_primitives()
  ! deallocate_primitives()
  ! deallocate_contracted()
  ! all moved to int_data_2cob3c_module
  !
  !**************************************************************

  !**************************************************************
  ! These were moved (with slight modifications) to
  ! integral_calc_quad_module:
  !
  !   subroutine contrsym2c *
  !   subroutine contrsym3c *
  !   subroutine contract_2center
  !   subroutine uncontract_2center
  !   subroutine contract_3center
  !   function contr2c *
  !   function contr3c
  !
  ! Only those marked by star are exported back to this sub.
  !**************************************************************

  !**************************************************************
  subroutine calc_primitives(U1,E1,L1,U2,E2,L2)
    ! calculate primitive integrals and do contraction
    use dip_prim_module, only: ll_momnt
    use spin_orbit_module !<<< options only
    use cpksdervs_matrices
    use calc3c_switches
    use interfaces
    use shgi_cntrl, only: IPSEU
    use shgi_cntrl, only: INUSD
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind), parameter  :: SS=0,LS=1,LL=2
    integer(kind=i4_kind)             :: ILL ! one of the above
    integer(kind=i8_kind)             :: IMODE
    !------------ Executable code ------------------------------

    ILL = LS
    if( L1==0 .and. L2==0 ) ILL = SS
    if( L1/=0 .and. L2/=0 ) ILL = LL

    IMODE = 0
! no more supported:
!ifdef WITH_OLD_PSEUDO
!   WARN('WITH_OLD_PSEUDO')
!   if( pseudopot_present ) IMODE = IOR(IMODE,IPSEU)
!   ! otherwise it will be computed by SHGI
!endif

#ifdef _NUSD_BY_CALC_3CENTER
    ! tell calc_3center() explicitly to add sec der of nuc:
    if( integralpar_dervs ) IMODE = IOR(IMODE,INUSD)
#endif

    ! calculate primitive integrals
    call start_timer(timer_int_prim_2cob3c(integralpar_i_int_part))
#ifdef WITH_EFP
    if (.not.integralpar_gradients .and. .not.integralpar_efp_gradients) then
#else
    if (.not.integralpar_gradients) then
#endif
       FPP_TIMER_START(oli)
#ifdef NEW_INTEGRALS
       call start_timer(timer_new_ll_integral)
       call calc_3center(U1, L1, U2, L2, IMODE=IMODE)
       call stop_timer(timer_new_ll_integral)
#endif

       ! FIXME: is this debugging?
       if(cpks_noco) then
         prim_int_3c_co = 0.0_r8_kind
       endif
       FPP_TIMER_STOP(oli)

       if (integralpar_2cob_pvec) then
          call ll_momnt(&
               & ua2_basis%exponents,   ua1_basis%exponents,&    !from int_data_2cob_
               & L2,L1,& !...
               & center2,     center1, &     !...
               & prim_2cv_ints(int_prim_pvec)%prim &
               & )
       endif

    else ! ========= GRADIENT PART ==========
       FPP_TIMER_START(olg)
       select case(ILL)
       case (SS)
          call dtrace("calc_primitives: calling ss_calculate_grads")
#ifndef NEW_INTEGRALS
          call start_timer(timer_ss_integral)
          call ss_calculate_grads(i_ea1,i_ea2,IMODE)
          call stop_timer(timer_ss_integral)
#else
          call start_timer(timer_new_ss_integral)
FPP_TIMER_START(t_calc_3center)

#if 1 /* commented to see memory use */
          call calc_3center(U1,L1,U2,L2, &
                            equalb=i_ea2,equala=i_ea1,IMODE=IMODE)
#endif

FPP_TIMER_STOP(t_calc_3center)
          call stop_timer(timer_new_ss_integral)
#endif
       case(LS)
          call dtrace("calc_primitives: calling ls_calculate_grads")
#ifndef NEW_INTEGRALS
          call start_timer(timer_ls_integral)
          call ls_calculate_grads(U1,U2,i_ea2,L1,L2,IMODE)
          call stop_timer(timer_ls_integral)
#else
          call start_timer(timer_new_ls_integral)
FPP_TIMER_START(t_calc_3center)

#if 1 /* commented to see memory use */
          call calc_3center(U1,L1,U2,L2,equalb=i_ea2,IMODE=IMODE)
#endif

FPP_TIMER_STOP(t_calc_3center)
          call stop_timer(timer_new_ls_integral)
#endif
       case (LL)
          call dtrace("calc_primitives: calling ll_calculate_grads")
#ifndef NEW_INTEGRALS
          call start_timer(timer_ll_integral)
          call ll_calculate_grads(U1,U2,i_ea2,L1,L2,IMODE)
          call stop_timer(timer_ll_integral)
#else
          call start_timer(timer_new_ll_integral)
FPP_TIMER_START(t_calc_3center)
#if 1 /* commented to see memory use */
          call calc_3center(U1,L1,U2,L2,equalb=i_ea2,IMODE=IMODE)
#endif
FPP_TIMER_STOP(t_calc_3center)
          call stop_timer(timer_new_ll_integral)
#endif
       case default
          call error_handler( &
               "calc_quad_2cob3c: calc_primitives:&
               & a wonder")
       end select
       FPP_TIMER_STOP(olg)
    endif

    DPRINT   'call calc_prim_shgi(',U1,E1,L1,U2,E2,L2,')'
    FPP_TIMER_START(shgi)
    if(integralpar_2cob_potential.or.integralpar_2cob_field) then
       call calc_prim_pot_and_field(U1,E1,L1,U2,E2,L2)
    else if (integralpar_solv_grad .or.integralpar_Q_solv_grad)then
       call calc_prim_solv(U1,E1,L1,U2,E2,L2)
    else

#ifdef WITH_EFP
      if(.not. integralpar_2cob_pc .and. &
         .not. integralpar_2cob_X  .and. &
         .not.integralpar_efp_gradients)then
#endif
         call calc_prim_shgi(U1,E1,L1,U2,E2,L2)
#ifdef WITH_EFP
      endif
#endif

       if( &
#ifdef WITH_EFP
          integralpar_2cob_pc .or. &
          integralpar_2cob_X .or. &
#endif
          integralpar_2cob_pc_grad.or. &
          integralpar_2cob_X_grad .or. &
          integralpar_2cob_ipd_grad) then
          call dtrace("calc_primitives: calling calc_prim_shgi_X")
          call calc_prim_shgi_X(U1,E1,L1,U2,E2,L2)
       end if

       ! calculation of forces in presence of electric field
       if (efield_applied() .and. integralpar_gradients) then
          call calc_grad_efield(U1,E1,L1,U2,E2,L2)
       endif
    endif
    FPP_TIMER_STOP(shgi)

    call dtrace("calc_primitives: primitive integrals done")

    call stop_timer(timer_int_prim_2cob3c(integralpar_i_int_part))
  end subroutine calc_primitives

!**************************************************************
  subroutine contract_and_symadapt(U1, E1, L1, U2, E2, L2, imode)
    ! calculate primitive integrals and do contraction
    use fit_coeff_module, only: fit_coeff_n_ch
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    integer(i4_kind), intent(in) :: imode
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)  :: GRAD_DIM,SPIN_GRAD_DIM
    integer(kind=i4_kind)  :: i_grad,k2dr
    integer(kind=i4_kind)  :: irel2c
    integer(kind=i4_kind)  :: irel3c
    integer(kind=i4_kind)  :: n_ff_ch ! fit_coeff_n_ch()
    !------------ Executable code ------------------------------

    call start_timer(timer_int_cont_2cob3c(integralpar_i_int_part))

    if (options_spin_orbit) then
       ! contraction of 3c_co is done inside:
       call symadapt_add_totalsym_spor(weight)
       GOTO 999 ! finalize and exit
    endif

    ! === NORMAL INTEGRALS ===
    ! default is imode, in relativistic calculations
    ! some 2c/3c integrals require special treatment:
    irel2c = imode
    if ( integralpar_relativistic ) irel2c = IOR(irel2c,IDONT_CONTRACT)

    irel3c = imode
    if ( is_on(op_FitTrafo)       ) irel3c = IOR(irel3c,IDONT_CONTRACT)

    if ( integralpar_2cob_kin ) then
       call contrsym2c(                      &
                    U1,E1,L1,U2,E2,L2,weight &
                   ,prim_int_2cob_kin        &
                   ,symadapt_int_2cob_kin    &
                   ,irel2c                   &
                   )
    endif
    if ( integralpar_2cob_nuc ) then
       call contrsym2c(                      &
                    U1,E1,L1,U2,E2,L2,weight &
                   ,prim_int_2cob_nuc        &
                   ,symadapt_int_2cob_nuc    &
                   ,irel2c                   &
                   )
       if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
       call contrsym2c(                           &
                    U1,E1,L1,U2,E2,L2,weight      &
                   ,prim_int_2cob_nuc_pseudo      &
                   ,symadapt_int_2cob_nuc_pseudo  &
                   ,irel2c                        &
                   )
       endif
    endif
    if ( integralpar_2cob_pvsp ) then
       call contrsym2c(                      &
                    U1,E1,L1,U2,E2,L2,weight &
                   ,prim_int_2cob_pvsp       &
                   ,symadapt_int_2cob_pvsp   &
                   ,irel2c                   &
                   )
    endif
    if ( integralpar_2cob_ol ) then
       call  contrsym2c(                      &
                     U1,E1,L1,U2,E2,L2,weight &
                    ,prim_int_2cob_ol         &
                    ,symadapt_int_2cob_ol     &
                    ,irel2c                   &
                    )
    endif
    if ( integralpar_3c_xc ) then
       call  contrsym3c(                      &
                     U1,E1,L1,U2,E2,L2,weight &
                    ,prim_int_3c_xc           &
                    ,symadapt_int_3c_xc       &
                    ,irel3c                   &
                    )
    endif
    if ( integralpar_3c_co ) then
       !
       ! FIXME: int TDDFT (response) calculations, the shape of
       ! the primitive fit integrals, prim_int_3c_co, deviates
       ! from normal. The middle axis has the length of a total
       ! number of fit functions including the non-totally symmetric
       ! once. As the totally symmetric irrep is usually the first
       ! the leading once are those used for density expansion:
       !
       n_ff_ch  = fit_coeff_n_ch()
       call  contrsym3c(                      &
                     U1,E1,L1,U2,E2,L2,weight &
                    ,prim_int_3c_co(:, :, :n_ff_ch, :, :) &
                    ,symadapt_int_3c_co       &
                    ,irel3c                   &
                    )
    endif
    if ( integralpar_2cob_potential ) then
       call  contrsym3c(                      &
                     U1,E1,L1,U2,E2,L2,weight &
                    ,prim_int_2cob_poten      &
                    ,symadapt_int_2cob_poten  &
                    ,imode                    &
                    )
    endif
    if ( integralpar_2cob_field ) then
       call  contrsym3c(                      &
                     U1,E1,L1,U2,E2,L2,weight &
                    ,prim_int_2cob_field      &
                    ,symadapt_int_2cob_field  &
                    ,imode                    &
                    )
    endif
    ! ========================

    ! ===    GRADIENTS     ===
    grad_dim=size(symadapt_int_2cob_ol_grad,1)

!   split_gradients = options_split_gradients()
!   model_density = options_xcmode() == xcmode_model_density .or. &
!                   options_xcmode() == xcmode_extended_mda
    if (model_density) then
       spin_grad_dim = grad_dim * options_n_spin()
    else
       spin_grad_dim = grad_dim
    endif

#ifdef WITH_EPE
    if ( integralpar_3cob_grad .and. ewpc_n.ne.0.and.calc_cluster_epe_energy) then
       call contrsym2c(                      &
                  U1,E1,L1,U2,E2,L2,weight   &
                 ,prim_int_3cob_epe%m        &
                 ,symadapt_int_3cob_epe(1,:) &
                 ,imode                      &
                 )
    endif
#endif

    do i_grad=1,GRAD_DIM
       if ( integralpar_2cob_ol_grad ) then
          call contrsym2c(                                  &
                     U1,E1,L1,U2,E2,L2,weight               &
                    ,prim_int_2cob_ol_grad(i_grad)%m        &
                    ,symadapt_int_2cob_ol_grad(i_grad,:)    &
                    ,imode                                  &
                    )

          if( integralpar_relativistic ) then
             call contrsym2c(                                  &
                       U1,E1,L1,U2,E2,L2,weight                &
                      ,prim_int_2cob_ol_grad(i_grad)%m         &
                      ,symadapt_int_2cob_ol_rel_grad(i_grad,:) &
                      ,irel2c                                  &
                      )
          endif
       endif
       if ( integralpar_2cob_kin_grad ) then
          call contrsym2c(                                &
                     U1,E1,L1,U2,E2,L2,weight             &
                    ,prim_int_2cob_kin_grad(i_grad)%m     &
                    ,symadapt_int_2cob_kin_grad(i_grad,:) &
                    ,irel2c                               &
                    )
       endif
    enddo

    do i_grad=1,gradient_data_n_gradients
       if ( integralpar_2cob_nuc_grad ) then
          call contrsym2c(                                &
                     U1,E1,L1,U2,E2,L2,weight             &
                    ,prim_int_2cob_nuc_grad(i_grad)%m     &
                    ,symadapt_int_2cob_nuc_grad(i_grad,:) &
                    ,irel2c                               &
                    )
          if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
             call contrsym2c(                                   &
                        U1,E1,L1,U2,E2,L2,weight                &
                       ,prim_int_2cob_pseudo_grad(i_grad)%m     &
                       ,symadapt_int_2cob_pseudo_grad(i_grad,:) &
                       ,irel2c                                  &
                       )
          endif
       endif
       if ( integralpar_2cob_pvsp_grad ) then
          call contrsym2c(                                 &
                     U1,E1,L1,U2,E2,L2,weight              &
                    ,prim_int_2cob_pvsp_grad(i_grad)%m     &
                    ,symadapt_int_2cob_pvsp_grad(i_grad,:) &
                    ,irel2c                                &
                    )
       endif
       if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
          call contrsym2c(                                 &
                     U1,E1,L1,U2,E2,L2,weight              &
                    ,prim_int_3cob_solv_grad(i_grad)%m     &
                    ,symadapt_int_3cob_solv_grad(i_grad,:) &
                    ,imode                                 &
                    )
       endif
       if ( integralpar_3cob_grad ) then
          call contrsym2c(                            &
                     U1,E1,L1,U2,E2,L2,weight         &
                    ,prim_int_3cob_grad(i_grad)%m     &
                    ,symadapt_int_3cob_grad(i_grad,:) &
                    ,imode                            &
                    )
       endif
    enddo

    if ( integralpar_2cob_pc_grad ) then
       do i_grad=1,totsym_grad_pc_length
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_pc_grad(i_grad)%m     &
                    ,symadapt_int_2cob_pc_grad(i_grad,:) &
                    ,imode                               &
                    )
       end do
    endif

    if ( integralpar_2cob_X_grad ) then
       do i_grad=1,totsym_grad_dip_length
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_pd_grad(i_grad)%m     &
                    ,symadapt_int_2cob_pd_grad(i_grad,:) &
                    ,imode                               &
                    )
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_pd_torq(i_grad)%m     &
                    ,symadapt_int_2cob_pd_torq(i_grad,:) &
                    ,imode                               &
                    )
       end do

       do i_grad=1,totsym_grad_quad_length
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_pq_grad(i_grad)%m     &
                    ,symadapt_int_2cob_pq_grad(i_grad,:) &
                    ,imode                               &
                    )
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_pq_torq(i_grad)%m     &
                    ,symadapt_int_2cob_pq_torq(i_grad,:) &
                    ,imode                               &
                    )
       end do

       do i_grad=1,totsym_grad_oct_length
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_po_grad(i_grad)%m     &
                    ,symadapt_int_2cob_po_grad(i_grad,:) &
                    ,imode                               &
                    )
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_po_torq(i_grad)%m     &
                    ,symadapt_int_2cob_po_torq(i_grad,:) &
                    ,imode                               &
                    )
       end do

       do i_grad=1,totsym_grad_rep_length
          call contrsym2c(                               &
                     U1,E1,L1,U2,E2,L2,weight            &
                    ,prim_int_2cob_rc_grad(i_grad)%m     &
                    ,symadapt_int_2cob_rc_grad(i_grad,:) &
                    ,imode                               &
                    )
       end do
    endif

    if ( integralpar_2cob_ipd_grad ) then
       do i_grad=1,totsym_grad_idip_length
          call contrsym2c(                                &
                     U1,E1,L1,U2,E2,L2,weight             &
                    ,prim_int_2cob_ipd_grad(i_grad)%m     &
                    ,symadapt_int_2cob_ipd_grad(i_grad,:) &
                    ,imode                                &
                    )
          call contrsym2c(                                &
                     U1,E1,L1,U2,E2,L2,weight             &
                    ,prim_int_2cob_ipd_torq(i_grad)%m     &
                    ,symadapt_int_2cob_ipd_torq(i_grad,:) &
                    ,imode                                &
                    )
       end do
    endif

    if ( integralpar_solv_grad .and. integralpar_cpksdervs .and. with_pc .and. .not.fixed_pc ) then
      do i_grad=1,TOTSYM_FIELD_LENGTH
          call contrsym2c(                                    &
                     U1,E1,L1,U2,E2,L2,weight                 &
                    ,prim_int_3cob_solv_grad_pc(i_grad)%m     &
                    ,symadapt_int_3cob_solv_grad_pc(i_grad,:) &
                    ,imode                                    &
                    )
      enddo
    endif

    if (integralpar_3cob_grad) then
       if (split_gradients) then
          do i_grad=1,SPIN_GRAD_DIM
             call contrsym2c(                               &
                        U1,E1,L1,U2,E2,L2,weight            &
                       ,prim_int_2cob_ks_grad(i_grad)%m     &
                       ,symadapt_int_2cob_ks_grad(i_grad,:) &
                       ,imode                               &
                       )
          enddo
          do i_grad=1,gradient_data_n_gradients
             call contrsym2c(                                &
                        U1,E1,L1,U2,E2,L2,weight             &
                       ,prim_int_3cob_nuc_grad(i_grad)%m     &
                       ,symadapt_int_3cob_nuc_grad(i_grad,:) &
                       ,imode                                &
                       )

             if ( model_density ) then
                call contrsym2c(                                 &
                           U1,E1,L1,U2,E2,L2,weight              &
                          ,prim_int_3cob_coul_grad(i_grad)%m     &
                          ,symadapt_int_3cob_coul_grad(i_grad,:) &
                          ,imode                                 &
                          )
             endif
          enddo
#ifndef no_cpks_coul_grads
       elseif(integralpar_cpksdervs) then
          ! used to calculate both cpks_matrixes(int1) and impl contribs(int2)
          do i_grad=1,gradient_data_n_gradients
             call contrsym3c(                                 &
                        U1,E1,L1,U2,E2,L2,weight              &
                       ,prim_int_cpks_coul_grad(i_grad)%m     &
                       ,symadapt_int_cpks_coul_grad(i_grad,:) &
                       ,imode                                 &
                       )
          enddo
#endif

       endif
    endif

    if (integralpar_solv_grad.and.integralpar_dervs) then     !!!!!!!!!!!!!!AS
      !
      ! FIXME: is this empty branch of any value?
      !
    else if(integralpar_Q_solv_grad) then
      ! FIXME: this empty branch is needed, otherwise
      !        in Q_SolvGrads the execution enters into
      !        the next "elseif" branch.
    elseif(integralpar_dervs) then

       do i_grad=1,gradient_data_n_spin_gradients
          do k2dr=1,gradient_data_n_spin_gradients
             call contrsym2c(                                  &
                        U1,E1,L1,U2,E2,L2,weight               &
                       ,prim_int_coul_dervs(i_grad,k2dr)%m     &
                       ,symadapt_int_coul_dervs(i_grad,k2dr,:) &
                       ,imode                                  &
                       )
          enddo
       enddo

       if (integralpar_2cob_ol_grad) then
          do i_grad=1,GRAD_DIM
             do k2dr=1,GRAD_DIM
                call contrsym2c(                                     &
                           U1,E1,L1,U2,E2,L2,weight                  &
                          ,prim_int_2cob_ol_dervs(i_grad,k2dr)%m     &
                          ,symadapt_int_2cob_ol_dervs(i_grad,k2dr,:) &
                          ,imode                                     &
                          )
             enddo
          enddo
       endif

       if ( integralpar_relativistic ) then
!      print *,'2cob3c: contrsym nucl_dervs, pvsp_dervs!'
       do i_grad=1,size(prim_int_nucl_dervs,1) ! gradient_data_n_spin_gradients
          do k2dr=1,size(prim_int_nucl_dervs,2) ! gradient_data_n_spin_gradients
             call contrsym2c(                                  &
                        U1,E1,L1,U2,E2,L2,weight               &
                       ,prim_int_nucl_dervs(i_grad,k2dr)%m     &
                       ,symadapt_int_nucl_dervs(i_grad,k2dr,:) &
                       ,irel2c                                 &
                       )
             call contrsym2c(                                  &
                        U1,E1,L1,U2,E2,L2,weight               &
                       ,prim_int_pvsp_dervs(i_grad,k2dr)%m     &
                       ,symadapt_int_pvsp_dervs(i_grad,k2dr,:) &
                       ,irel2c                                 &
                       )
          enddo
       enddo
       do i_grad=1,size(prim_int_2cob_kin_dervs,1) ! GRAD_DIM
          do k2dr=1,size(prim_int_2cob_kin_dervs,2) ! GRAD_DIM
             ! (uncontracted) overlap...
             call contrsym2c(                                      &
                        U1,E1,L1,U2,E2,L2,weight                   &
                       ,prim_int_2cob_ol_dervs(i_grad,k2dr)%m      &
                       ,symadapt_int_2cob_olu_dervs(i_grad,k2dr,:) &
                       ,irel2c                                     &
                       )
             ! ...and kinetic energy:
             call contrsym2c(                                      &
                        U1,E1,L1,U2,E2,L2,weight                   &
                       ,prim_int_2cob_kin_dervs(i_grad,k2dr)%m     &
                       ,symadapt_int_2cob_kin_dervs(i_grad,k2dr,:) &
                       ,irel2c                                     &
                       )
          enddo
       enddo
       endif
    endif
    ! ========================

999 CONTINUE ! FINILIZE AND EXIT!
    call stop_timer(timer_int_cont_2cob3c(integralpar_i_int_part))

    call dtrace("calc_primitives: contraction done")
  end subroutine contract_and_symadapt
  !**************************************************************

  !**************************************************************
  subroutine relstore(U1,L1,U2,L2)
    use gradient_data_module, only: gradient_number
    use int_data_2cob3c_module, only: g=>grad_dim_index
    use relgrads_store
    implicit none
    integer(IK), intent(in) :: U1,L1,U2,L2
    ! *** end of interface ***

    integer(IK) :: irr,i,j
    integer(IK) :: ng1,ng2,ng12

    if( .not. integralpar_2cob_nuc ) goto 111 ! jump to grads

    ! save only one triangle:
    ASSERT(U1>=U2)
    if( U1==U2 )then
    ASSERT(L1>=L2)
    endif

    do irr=1,size(symadapt_int_2cob_nuc)  ! N_IRREPS
      call rg_put( RGNUCL,irr,U2,L2,U1,L1, symadapt_int_2cob_nuc(irr)%int )
      call rg_put( RGPVSP,irr,U2,L2,U1,L1, symadapt_int_2cob_pvsp(irr)%int )
      call rg_put( RGOVRL,irr,U2,L2,U1,L1, symadapt_int_2cob_ol(irr)%int )
      call rg_put( RGKNTC,irr,U2,L2,U1,L1, symadapt_int_2cob_kin(irr)%int )
      if( pseudopot_present .or. integralpar_pseudo) then
      call rg_put( RGPSEU,irr,U2,L2,U1,L1, symadapt_int_2cob_nuc_pseudo(irr)%int )
      endif
    enddo

111 if( .not. integralpar_gradients .and. .not. integralpar_dervs ) RETURN

    ! FIXME: this is a tribute to the special way
    !        to store 2c kin and ovrl grads:
    ng1 = gradient_number(U1)
    ng2 = gradient_number(U2)

    ng12= ng1 + ng2
!   if( U1==U2 ) ng12 = ng12/2

    if( .not. integralpar_gradients ) GOTO 222 ! jump to sec-ders.
!   print *,'2cob3c: SAVE grads ua=',U1,' la=',L1,' ub=',U2,' lb=',L2

    i = size(symadapt_int_2cob_kin_grad,1)
    ASSERT(i==ng1+ng2)
    i = size(symadapt_int_2cob_ol_rel_grad,1)
    ASSERT(i==ng1+ng2)

    do irr=1,size(symadapt_int_2cob_nuc_grad,2)  ! N_IRREPS
       do i=1,size(symadapt_int_2cob_nuc_grad,1) ! gradient_data_n_spin_gradients
         call rg_gr_put( RGNUCL,irr,U2,L2,U1,L1,i,symadapt_int_2cob_nuc_grad(i,irr)%int    )
         call rg_gr_put( RGPVSP,irr,U2,L2,U1,L1,i,symadapt_int_2cob_pvsp_grad(i,irr)%int   )
         if( pseudopot_present .or. integralpar_pseudo) then
         call rg_gr_put( RGPSEU,irr,U2,L2,U1,L1,i,symadapt_int_2cob_pseudo_grad(i,irr)%int )
         endif
       enddo

       ! FIXME: i = i + offset(U1)
       do i=1,size(symadapt_int_2cob_kin_grad,1) !ng12 ! GRAD_DIM or GRAD_DIM/2
         call rg_gr_put( RGKNTC,irr,U2,L2,U1,L1,g(i),symadapt_int_2cob_kin_grad(i,irr)%int    )
         call rg_gr_put( RGOVRL,irr,U2,L2,U1,L1,g(i),symadapt_int_2cob_ol_rel_grad(i,irr)%int )
       enddo
    enddo

222 if( .not. integralpar_dervs ) RETURN ! === RETURN POINT ===
!   print *,'2cob3c: dervs ua=',U1,' la=',L1,' ub=',U2,' lb=',L2

    do irr=1,size(symadapt_int_nucl_dervs,3)     ! N_IRREPS

       do j=1,size(symadapt_int_nucl_dervs,2) ! gradient_data_n_spin_gradients
       do i=1,j ! UPPER TRIANGLE (or size(symadapt_int_nucl_dervs,1) if you dont care)
         ASSERT(i<=j)
         call rg_sd_put( RGNUCL,irr,U2,L2,U1,L1,i,j,symadapt_int_nucl_dervs(i,j,irr)%int )
         call rg_sd_put( RGPVSP,irr,U2,L2,U1,L1,i,j,symadapt_int_pvsp_dervs(i,j,irr)%int )
       enddo
       enddo

       do j=1,ng12 ! GRAD_DIM or GRAD_DIM/2
       do i=1,ng12 ! GRAD_DIM or GRAD_DIM/2
         if( g(i) > g(j) ) CYCLE ! UPPER TRIANGLE
         ASSERT(g(i)<=g(j))
         call rg_sd_put( RGKNTC,irr,U2,L2,U1,L1,g(i),g(j),symadapt_int_2cob_kin_dervs(i,j,irr)%int )
         call rg_sd_put( RGOVRL,irr,U2,L2,U1,L1,g(i),g(j),symadapt_int_2cob_olu_dervs(i,j,irr)%int )
       enddo
       enddo
    enddo
  end subroutine relstore
  !**************************************************************

  subroutine calc_prim_relfit(U1,E1,L1,U2,E2,L2)
    use spin_orbit_module, only: &
         whatis, op_RelFit, SO_RELFIT_1C
    use shgi, only: shgi_drv
    use shgi_cntrl, only: IOVRL,IKNTC,INUCL,INUSR,INUSO &
                                           ,ICHSR,ICHSO
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------

    logical                           :: diag
    integer(i8_kind) :: IMODE
    !------------ Executable code ------------------------------

    DPRINT '2cob3c::calc_prim_relfit: ',U1,E1,L1,U2,E2,L2

    if (.not.is_on(op_RelFit)) then
       WARN('calc_prim_relfit: why calling?')
       return
    endif

    ! by default:
    IMODE = IOVRL+IKNTC+INUCL+INUSR+INUSO + ICHSR+ICHSO

    if(whatis(op_RelFit).eq.SO_RELFIT_1C)then
       diag = (U1==U2) .and. (E1==E2)
       if(.not.diag)then
          ! S :
          prim_3c_ints( int_prim_rcoul_pvsp)%prim = 0.0
          print *,'EE: zero <',U1,E1,L1,'|ps S  p|',U2,E2,L2,'>'

          prim_3cv_ints(int_prim_rcoul_pvxp)%prim = 0.0
          print *,'EE: zero <',U1,E1,L1,'|px S  p|',U2,E2,L2,'>'

          ! R2:
          prim_3c_ints( int_prim_r2_pvsp   )%prim = 0.0
          print *,'EE: zero <',U1,E1,L1,'|ps R2 p|',U2,E2,L2,'>'

          prim_3cv_ints(int_prim_r2_pvxp   )%prim = 0.0
          print *,'EE: zero <',U1,E1,L1,'|px R2 p|',U2,E2,L2,'>'

          IMODE = 0
       else
          print *,'EE: KEEP <',U1,E1,L1,'|XXXXXXX|',U2,E2,L2,'> (same atom)'

          ! will be computed only if same center (A==B):
          IMODE = ICHSR+ICHSO
       endif
       ! always:
       IMODE = IMODE + IOVRL+IKNTC+INUCL+INUSR+INUSO
    endif

    ! FinNuc is controlled by UA%NUCLEAR_RADIUS > 0:
    call shgi_drv(U1,E1,L1,U2,E2,L2 &
          , unique_atoms &
          , OVRL=prim_2c_ints( int_prim_ol        )%prim &
          , KNTC=prim_2c_ints( int_prim_kin       )%prim &
          , NUCL=prim_2c_ints( int_prim_nuc       )%prim &
          , NUSR=prim_2c_ints( int_prim_pvsp      )%prim &
          , NUSO=prim_2cv_ints(int_prim_pvxp      )%prim &
          , CHSR=prim_3c_ints( int_prim_rcoul_pvsp)%prim &
          , R2SR=prim_3c_ints( int_prim_r2_pvsp   )%prim &
          , CHSO=prim_3cv_ints(int_prim_rcoul_pvxp)%prim &
          , R2SO=prim_3cv_ints(int_prim_r2_pvxp   )%prim &
          , IMOD=IMODE                                   &
         )
  end subroutine calc_prim_relfit

  subroutine calc_prim_shgi(U1,E1,L1,U2,E2,L2)
    use spin_orbit_module, only: is_on, op_RelFit
    use shgi, only: shgi_drv, shgi_gr_drv, shgi_sd_drv
    use shgi_cntrl, only: INUCL,INUGR,INUSD, &
                          INUSR,ISRGR,ISRSD, &
                          IPSEU,IPSGR,IPSSD, &
                          IOVRL,IOVGR,IOVSD, &
                          IKNTC,IKNGR,IKNSD, &
                          IPCSD,INUSO,IADKH
#ifdef WITH_EPE
    use ewaldpc_module, only: ewpc_n, ewpc_array
#endif
    use datatype, only: pct=>pointcharge_type !!!!!!!!!!!!!AS
    use pointcharge_module, only: pointcharge_array, pointcharge_N !!!!!!!!!AS
    use options_module, only: option!(adkh)
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------
    type(pct), pointer :: pc(:) => NULL()
    integer(i8_kind) :: IMODE
    !------------ Executable code ------------------------------

    IMODE = 0

#ifdef WITH_EPE
    ! FIXME: workaround for ewpc_n==0 special case:
    if (ewpc_n /= 0 .and. pointcharge_N /= 0) &
         call error_handler("Ewald and Regular PC still cannot be used simultaneously")
    if(ewpc_n /= 0) pc=>ewpc_array
#endif

    ! This array is zero-sized when there are no point charges:
    pc => pointcharge_array
    ASSERT(associated(pc))

    ! jump to gradients:
    if( integralpar_gradients ) GOTO 111

    !========= INTEGRAL ENTRY POINT =========================

    IMODE = 0

    ! overlap:
    if(integralpar_2cob_ol)then
       IMODE = IOR(IMODE,IOVRL)
    endif

    ! kinetic:
    if(integralpar_2cob_kin)then
       IMODE = IOR(IMODE,IKNTC)
    endif

    ! nuclear attraction now from here:
    if(integralpar_2cob_nuc)then
       IMODE = IOR(IMODE,INUCL)

       if(integralpar_2cob_pvsp)then ! FIXME: SR w/o NUC?
          IMODE = IOR(IMODE,INUSR)
       endif

       ! spin orbit integrals
       if (integralpar_2cob_pvxp) then
          IMODE = IOR(IMODE,INUSO)
       endif

       ! new PP code called from here
       if(pseudopot_present .or. integralpar_pseudo)then ! FIXME: PP w/o NUC?
          IMODE = IOR(IMODE,IPSEU)
       endif

       if( option('adkh') )then
          IMODE = IOR(IMODE,IADKH)
       endif
    endif

    if ( is_on(op_RelFit) ) then
       call calc_prim_relfit(U1,E1,L1,U2,E2,L2)
       ! RETURN POINT FOR RELFIT, FIXME: merge case list instead!
       GOTO 999 ! clean up and exit
    endif

    SELECT CASE( IAND(IMODE,IOVRL+IKNTC+INUCL+INUSR+INUSO+IPSEU) )
    CASE (0)
       ! do nothing here, but see below after END SELECT !
    CASE (IOVRL)
       ! I guess it is unnecessary
       !WARN('calc_prim_shgi: do you really need OVRL?')
    CASE (IOVRL+IKNTC+INUCL, IOVRL+IKNTC+INUCL+IPSEU)
       ! in non-rel. NUCL and PSEU go together:
       call shgi_drv(U1,E1,L1,U2,E2,L2 &
             , unique_atoms &
             , pc &
             , OVRL=prim_int_2cob_ol       &
             , KNTC=prim_int_2cob_kin      &
             , NUCL=prim_int_2cob_nuc      &
             , IMOD=IMODE                  &
             )
    CASE (IOVRL+IKNTC+INUCL+INUSR)
       call shgi_drv(U1,E1,L1,U2,E2,L2 &
             , unique_atoms &
             , pc  &
             , OVRL=prim_int_2cob_ol          &
             , KNTC=prim_int_2cob_kin         &
             , NUCL=prim_int_2cob_nuc         &
             , NUSR=prim_int_2cob_pvsp        &
             , IMOD=IMODE                     &
             )
    CASE (IOVRL+IKNTC+INUCL+INUSR+INUSO)
       call shgi_drv(U1,E1,L1,U2,E2,L2 &
             , unique_atoms &
             , pc  &
             , OVRL=prim_2c_ints( int_prim_ol  )%prim &
             , KNTC=prim_2c_ints( int_prim_kin )%prim &
             , NUCL=prim_2c_ints( int_prim_nuc )%prim &
             , NUSR=prim_2c_ints( int_prim_pvsp)%prim &
             , NUSO=prim_2cv_ints(int_prim_pvxp)%prim &
             , IMOD=IMODE                             &
             )
    CASE (IOVRL+IKNTC+INUCL+INUSR+IPSEU)
       ! PP and Nuclear attraction separated,
       ! PP goes into separate prim_int_2cob_nuc_pseudo
       call shgi_drv(U1,E1,L1,U2,E2,L2 &
             , unique_atoms &
             , pc &
             , OVRL=prim_int_2cob_ol          &
             , KNTC=prim_int_2cob_kin         &
             , NUCL=prim_int_2cob_nuc         &
             , NUSR=prim_int_2cob_pvsp        &
             , PSEU=prim_int_2cob_nuc_pseudo  &
             , IMOD=IMODE                     &
             )
    CASE DEFAULT
       WARN('no such in ints')
       goto 666 ! print bits and abort
!      print '(" ERROR: in ints IMODE=",B32)',IMODE
!      ABORT('no such in ints')
    END SELECT

    ! RETURN POINT FOR ENERGY CALC (NO GRADS):
    GOTO 999 ! clean up and exit

    !========= GRADIENTS ENTRY POINT ===============
111 CONTINUE

    IMODE = 0

    ! overlap:
    if(integralpar_2cob_ol_grad)then
       IMODE = IOR(IMODE,IOVGR)
    endif

    if(integralpar_3cob_grad)then
       ! nuclear attraction:
       IMODE = IOR(IMODE,INUGR)
       ! kinetic:
       IMODE = IOR(IMODE,IKNGR)

       if(pseudopot_present .or. integralpar_pseudo)then
          ! PP goes only with NUC:
          IMODE = IOR(IMODE,IPSGR)
       endif

       if(integralpar_relativistic)then
          ! PVSP goes only with NUC:
          IMODE = IOR(IMODE,ISRGR)
       endif
    endif

    if( option('adkh') )then
       IMODE = IOR(IMODE,IADKH)
    endif

     ! Of course it is faster to share code between
     ! gradients and second derivatives ...
     ! As soon as the gradient driver shgi_gr_drv()
     ! will cope with ALL second derivatives we may think to
     ! merge the case lists. So far the the sec der driver
     ! is called anyway ...

    SELECT CASE( IAND(IMODE,IOVGR+IKNGR+INUGR+IPSGR+ISRGR) )
    CASE (0)
       ! do nothing here, but see below after END SELECT !
    CASE (IOVGR)
       ! I guess it is unnecessary
       !WARN('calc_prim_shgi: do you really need OVGR?')
    CASE (IOVGR+IKNGR+INUGR, IOVGR+IKNGR+INUGR+IPSGR)
       ! in non-rel. NUCL and PSEU go together:
       call shgi_gr_drv(U1,E1,L1,U2,E2,L2 &    !!!(2)
             , unique_atoms &
             , pc &
             , OVGR=prim_int_2cob_ol_grad      &
             , NUGR=prim_int_3cob_grad         & ! WARNING: non-split grads
             , IMOD=IMODE                      &
             )
    CASE (IOVGR+IKNGR+INUGR+ISRGR)
       ! in rel. KIN, NUCL (and PSEU) must be separated:
       call shgi_gr_drv(U1,E1,L1,U2,E2,L2 &    !!!(3)
             , unique_atoms &
             , pc &
             , OVGR=prim_int_2cob_ol_grad      &
             , KNGR=prim_int_2cob_kin_grad     &
             , NUGR=prim_int_2cob_nuc_grad     & ! WARNING: different from the above!
             , SRGR=prim_int_2cob_pvsp_grad    &
             , IMOD=IMODE                      &
             )
    CASE (IOVGR+IKNGR+INUGR+ISRGR+IPSGR)
       ! PP and Nuclear attraction separated,
       ! PP goes into separate prim_int_2cob_pseudo_grad
       call shgi_gr_drv(U1,E1,L1,U2,E2,L2 &    !!!(4)
             , unique_atoms &
             , pc &
             , OVGR=prim_int_2cob_ol_grad      &
             , KNGR=prim_int_2cob_kin_grad     &
             , NUGR=prim_int_2cob_nuc_grad     &
             , SRGR=prim_int_2cob_pvsp_grad    &
             , PSGR=prim_int_2cob_pseudo_grad  &
             , IMOD=IMODE                      &
             )
    CASE DEFAULT
       WARN('no such in grads')
       goto 666 ! print bits and abort
    END SELECT

    ! RETURN POINT FOR GRADIENT CALC:
    if( .not. integralpar_dervs ) GOTO 999 ! clean up and exit

    !========= SECOND DERIVATIVES ENTRY POINT ===============

    IMODE = 0

       ! overlap:
       IMODE = IOR(IMODE,IOVSD)

       ! kinetic:
       IMODE = IOR(IMODE,IKNSD)

#ifndef _NUSD_BY_CALC_3CENTER
       ! nuclear attraction, tell shgi explicitly to add sec der of nuc:
       IMODE = IOR(IMODE,INUSD)
#endif

       if(pseudopot_present)then
          ! PP goes only with NUC:
          IMODE = IOR(IMODE,IPSSD)
       endif

       if(integralpar_relativistic)then
          ! PVSP goes only with NUC:
          IMODE = IOR(IMODE,ISRSD)
       endif

       if( option('adkh') )then
          IMODE = IOR(IMODE,IADKH)
       endif

#ifdef WITH_EPE
       if( ewpc_n.ne.0 ) then
          IMODE = IOR(IMODE,IPCSD)
!         if(.not.pseudopot_present) IMODE = IOR(IMODE,IPSSD) !AS:correct me if it is bad
!                                                             !but now we can calculate PC secnd deriv
!                                                             !if no pseudopotentials in QM system
       end if
#endif
       if( pointcharge_N.ne.0 ) then
          IMODE = IOR(IMODE,IPCSD)
!         if(.not.pseudopot_present) IMODE = IOR(IMODE,IPSSD) !The same as above
       end if

   SELECT CASE( IAND( IMODE, IPSSD+ISRSD ) )
    CASE ( 0, IPSSD )
       ASSERT(IMODE/=0)
       ! in non-rel. KIN, NUCL, and PSEU go together:
       call shgi_sd_drv(U1,E1,L1,U2,E2,L2      &
             , unique_atoms                    &
             , pcs=pc &
             , OVSD=prim_int_2cob_ol_dervs     &
             , NUSD=prim_int_coul_dervs        &
             , IMOD=IMODE                      &
             )
    CASE ( ISRSD )
       ! in REL case NUCL and KIN go separate:
       call shgi_sd_drv(U1,E1,L1,U2,E2,L2      &
             , unique_atoms                    &
             , OVSD=prim_int_2cob_ol_dervs     &
             , KNSD=prim_int_2cob_kin_dervs    &
             , NUSD=prim_int_nucl_dervs        &
             , SRSD=prim_int_pvsp_dervs        &
             , IMOD=IMODE                      &
             )
    CASE ( ISRSD + IPSSD )
       ! PP and Nuclear attraction separated,
       ! PP goes into separate prim_int_2cob_??
       call shgi_sd_drv(U1,E1,L1,U2,E2,L2      &
             , unique_atoms                    &
             , pcs=pc &
             , OVSD=prim_int_2cob_ol_dervs     &
             , KNSD=prim_int_2cob_kin_dervs    &
             , NUSD=prim_int_nucl_dervs        &
             , SRSD=prim_int_pvsp_dervs        &
             , PSSD=prim_int_coul_dervs        &
             , IMOD=IMODE                      &
             )
    CASE DEFAULT
       WARN('no such in dervs')
       goto 666 ! print bits and abort
    END SELECT

999 CONTINUE ! clean up and exit:
    nullify(pc)
    RETURN

666 CONTINUE ! error occured
       print '(" ERROR: wrong IMODE=",B32)',IMODE

       print '("              IKNTC=",B32)',IKNTC
       print '("              INUCL=",B32)',INUCL
       print '("              IPSEU=",B32)',IPSEU
       print '("              IOVRL=",B32)',IOVRL
       print '("              INUSR=",B32)',INUSR

       print '("              IKNGR=",B32)',IKNGR
       print '("              INUGR=",B32)',INUGR
       print '("              IPSGR=",B32)',IPSGR
       print '("              IOVGR=",B32)',IOVGR
       print '("              ISRGR=",B32)',ISRGR

       print '("              IKNSD=",B32)',IKNSD
       print '("              INUSD=",B32)',INUSD
       print '("              IPSSD=",B32)',IPSSD
       print '("              IOVSD=",B32)',IOVSD
       print '("              ISRSD=",B32)',ISRSD
       ABORT('no such mode')
  end subroutine calc_prim_shgi
  !**************************************************************
  subroutine calc_prim_pot_and_field(U1,E1,L1,U2,E2,L2)
    use shgi_ep_ef, only: shgi_pot_drv,shgi_field_drv
    use potential_module, only: point_in_space
    use elec_static_field_module, only: surface_points
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------
    !------------ Executable code ------------------------------

#ifdef NEW_INTEGRALS
    if(integralpar_2cob_potential)then
       ! Electrostatic potential
       call shgi_pot_drv(U1,E1,L1,U2,E2,L2, &
            unique_atoms,point_in_space,    &
            prim_int_2cob_poten)
    endif
    if(integralpar_2cob_field)then
       ! Electrostatic field
       call shgi_field_drv(U1,E1,L1,U2,E2,L2, &
            unique_atoms,surface_points,    &
            prim_int_2cob_field)
    endif
#else
    ! FIXME: rm e-field from ??_calculate (WHY NOT)
#endif

  end subroutine calc_prim_pot_and_field
  !**************************************************************

  !**************************************************************
  subroutine calc_prim_solv(U1,E1,L1,U2,E2,L2)
    use shgi, only: shgi_gr_solv_drv, shgi_gr_Q_solv_drv, shgi_gr_solv_drv_vtn, shgi_sd_solv_drv
    use solv_cavity_module, only: grad_solv_totsym, grad_solv_totsym_tes,VTN
#ifdef WITH_EFP
    use efp_solv_grad_module, only: gradient_mpole_solv_totalsym
    use efp_solv_grad_module, only: torque_mpole_solv_totalsym
#endif
    use elec_static_field_module, only: totalsym_field
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------
    !------------ Executable code ------------------------------

#ifdef NEW_INTEGRALS
    if(integralpar_solv_grad)then
       ! Gradients of solvation
       if(integralpar_cpksdervs)then
         !
         ! a) In a SecDer calculation, one needs matrices for
         ! RHS of CPKS equation. b) In a post-CPKS run
         ! the traces tr(P^X,I^Y) are handled the old way,
         ! through matrices. FIXME the case (b)!
         !
         DPRINT  '2cob3c: shgi_gr_solv_drv(matrix)'
         call shgi_gr_solv_drv(U1,E1,L1,U2,E2,L2, &
              unique_atoms,                       &
              MATUAGR=prim_int_3cob_solv_grad,    &
              MATPCGR=prim_int_3cob_solv_grad_pc  )
       else
         !
         ! In calculation of forces pass the uncontracted
         ! AO density matrix and form the traces tr(P,I^X) early.
         ! Return just the scalar gradients:
         !
         if(VTN) then
            DPRINT  '2cob3c: shgi_gr_solv_drv_vtn(vector)'
#ifdef WITH_EFP
            call shgi_gr_solv_drv_vtn(U1,E1,L1,U2,E2,L2, &
                 unique_atoms,                           &
                 prim_int_dens_mat,                      &
                 grad_solv_totsym,                       &
                 gradient_mpole_solv_totalsym,           &
                 torque_mpole_solv_totalsym              )
#else
            call shgi_gr_solv_drv_vtn(U1,E1,L1,U2,E2,L2, &
                 unique_atoms,                           &
                 prim_int_dens_mat,                      &
                 grad_solv_totsym                        )
#endif
         else
            DPRINT  '2cob3c: shgi_gr_solv_drv(vector)'
            call shgi_gr_solv_drv(U1,E1,L1,U2,E2,L2, &
                 unique_atoms,                       &
                 DENSMAT=prim_int_dens_mat,          &
                 VECUAGR=grad_solv_totsym,           &
                 VECPCGR=totalsym_field              )
         end if
       endif
    endif
#else
    ! FIXME: rm solvation gradients from ??_calculate_grads (WHY NOT)
    !        The storage prim_int_3cob_solv_grad/symadapt_* is only
    !        Allocated if( integralpar_cpksdervs .and. integralpar_cpksdervs )
    !        That is, only in secder calculations, but old ints dont consider this.
    ABORT('needs fixes')
#endif

    ! RETURN POINT FOR GRADIENT CALC:
    if( .not. integralpar_dervs ) RETURN

    if(integralpar_solv_grad)then
       !
       ! Explicit contribution of solvation to second derivatives
       ! Pass AO density matrix and let the driver form/return
       ! the traces tr(P,I^XY) immediately:
       !
       call shgi_sd_solv_drv(U1,E1,L1,U2,E2,L2      &
             , unique_atoms                         &
             , prim_int_dens_mat                    &
             , dervs_totalsym                       &
             )
    else if(integralpar_Q_solv_grad) then
       !
       ! Gradients wrt PCs of solvation tessarea (== electric field there)
       ! Pass AO density matrix and let the driver form/return
       ! the traces tr(P,I^X) immediately:
       !
       call shgi_gr_Q_solv_drv(U1,E1,L1,U2,E2,L2             &
                              , unique_atoms                 &
                              , prim_int_dens_mat            &
                              , grad_solv_totsym_tes         &
                              )
    end if
  end subroutine calc_prim_solv
  !**************************************************************

  !**************************************************************
  subroutine calc_grad_efield(U1,E1,L1,U2,E2,L2)
    use shgi_dip, only: shgi_gr_efield_drv
    use efield_module, only:  efield_field
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! uniqe atom, equiv. atom, and ang momenta for center 1 and 2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------
    !------------ Executable code ------------------------------

    print*, "calc_grad_efield enter"
    !
    ! to calculate the forces in presence of electric field
    !
ASSERT(allocated(prim_int_dens_mat))
ASSERT(allocated(gradient_totalsym))
    call shgi_gr_efield_drv(U1,E1,L1,U2,E2,L2  &
                           , unique_atoms      &
                           , prim_int_dens_mat &
                           , efield_field()    &
                           , gradient_totalsym )
    print*, "calc_grad_efield exit"
  end subroutine calc_grad_efield
  !**************************************************************

  !**************************************************************
  subroutine calc_prim_shgi_X(U1,E1,L1,U2,E2,L2)
    use pointcharge_module, only: pointcharge_array
    use point_dqo_module, only: pd_array,pq_array,po_array,rc_array
    use induced_dipoles_module, only: ipd_array
    use shgi_pcm, only: shgi_pc_grad
    use shgi_ext_c, only: shgi_X_grad,shgi_X_torq
#ifdef WITH_EFP
    use shgi_pcm, only: shgi_PCs_wrap
    use shgi_ext_c, only: shgi_X_wrap
#endif
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! *** end of interface ***
    !------------ Declaration of local variables ---------------
    !------------ Executable code ------------------------------

#ifdef NEW_INTEGRALS

#ifdef WITH_EFP
    if(associated(pointcharge_array)) then
    if(integralpar_2cob_pc) &
         call shgi_PCs_wrap(U1,E1,L1,U2,E2,L2  &
                            ,unique_atoms      &
                            ,prim_int_2cob_nuc)
    endif

    if(associated(pd_array) .or. associated(pq_array) .or. &
       associated(po_array) .or. associated(rc_array)) then
    if(integralpar_2cob_X)  &
         call shgi_X_wrap(U1,E1,L1,U2,E2,L2  &
                          ,unique_atoms      &
                          ,prim_int_2cob_nuc)
    endif
#endif

    if(integralpar_2cob_pc_grad) &
         call shgi_pc_grad(U1,E1,L1,U2,E2,L2      &
                          ,unique_atoms           &
                          ,pointcharge_array      &
                          ,prim_int_2cob_pc_grad)

    if(integralpar_2cob_X_grad.or.integralpar_2cob_ipd_grad) then
       call shgi_X_grad(U1,E1,L1,U2,E2,L2       &
                        ,unique_atoms           &
                        ,pd_array               &
                        ,pq_array               &
                        ,po_array               &
                        ,rc_array               &
                        ,ipd_array              &
                        ,prim_int_2cob_pd_grad  &
                        ,prim_int_2cob_pq_grad  &
                        ,prim_int_2cob_po_grad  &
                        ,prim_int_2cob_rc_grad  &
                        ,prim_int_2cob_ipd_grad)

       call shgi_X_torq(U1,E1,L1,U2,E2,L2       &
                        ,unique_atoms           &
                        ,pd_array               &
                        ,pq_array               &
                        ,po_array               &
                        ,ipd_array              &
                        ,prim_int_2cob_pd_torq  &
                        ,prim_int_2cob_pq_torq  &
                        ,prim_int_2cob_po_torq  &
                        ,prim_int_2cob_ipd_torq)
    end if
#endif

  end subroutine calc_prim_shgi_X
  !**************************************************************

  !**************************************************************
  ! These were moved (with slight modifications) to
  ! integral_calc_quad_module:
  !
  !   subroutine symadp2c
  !   subroutine unsymadp2c
  !   subroutine symadp3c
  !
  ! None of them is directly called from this sub.
  !**************************************************************

  !**************************************************************
  subroutine symadapt_add_totalsym_spor(weight)
    use type_module,only: IK=>i4_kind,RK=>r8_kind
    use symm_adapt_struct, only: LSymAdp, LSymAdpL, LSymAdpS
    use symm_adapt_module, only: symm_adapt_2c, symm_adapt_2cv, &
                                 symm_adapt_3c, symm_adapt_3cv
    use spin_orbit_module, only: is_on, op_FitTrafo
    implicit none
    real(RK), intent(in) :: weight
    ! *** end of interface ***

    integer(IK)                            :: L1,L2,ea1,ea2
    integer(IK)                            :: ua1,ua2
!   real(RK)                               :: weight

    L1  = quadrupel%l1; ! public ?
    L2  = quadrupel%l2
    ua1 = quadrupel%ua1
    ua2 = quadrupel%ua2
    ea1 = i_ea1
    ea2 = i_ea2

!   weight = symequivatoms%weight(i_symequiv); ! public ?

    ! TWO CENTER, SCALARS:
    if(integralpar_2cob_kin)then
      call symm_adapt_2c( prim_2c_ints(int_prim_kin)%prim, sa_2c_ints(:, int_sa_kin), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_2cob_nuc)then
      call symm_adapt_2c( prim_2c_ints(int_prim_nuc)%prim, sa_2c_ints(:, int_sa_nuc), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_2cob_pvsp)then
      call symm_adapt_2c( prim_2c_ints(int_prim_pvsp)%prim, sa_2c_ints(:, int_sa_pvsp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_2cob_ol)then
      call symm_adapt_2c( prim_2c_ints(int_prim_ol)%prim, sa_2c_ints(:, int_sa_ol), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif

    ! TWO CENTER, VECTORS
    if(integralpar_2cob_pvxp)then
      call symm_adapt_2cv( prim_2cv_ints(int_prim_pvxp)%prim, sa_2c_ints(:, int_sa_pvxp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_2cob_pvec)then
      call symm_adapt_2cv( prim_2cv_ints(int_prim_pvec)%prim, sa_2c_ints(:, int_sa_sigp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif

    ! THREE CENTER, SCALARS:
    if(integralpar_3c_xc)then
      ABORT('FIXME!')
      call symm_adapt_3c( prim_3c_ints(int_prim_xc)%prim, sa_3c_ints(:, int_sa_xc), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_3c_co)then
      ! FIXME: contraction if no FitTrafo!?
      if( is_on(op_FitTrafo) )then
        ! DONT contract 3c_co:
        call symm_adapt_3c(         prim_3c_ints(int_prim_co)%prim , sa_3c_ints(:, int_sa_co), &
                         ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
      else
        ! contract 3c_co:
        call symm_adapt_3c( contr3c(prim_3c_ints(int_prim_co)%prim), sa_3c_ints(:, int_sa_co), &
                         ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
      endif
    endif

    ! THREE CENTER, SCALARS, continued
    if(integralpar_3c_rcoul_pvsp)then
      call symm_adapt_3c( prim_3c_ints(int_prim_rcoul_pvsp)%prim, sa_3c_ints(:, int_sa_rcoul_pvsp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_3c_r2_pvsp)then
      call symm_adapt_3c( prim_3c_ints(int_prim_r2_pvsp)%prim, sa_3c_ints(:, int_sa_r2_pvsp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif

    ! Three CENTER, VECTORS :
    if(integralpar_3c_rcoul_pvxp)then
      call symm_adapt_3cv( prim_3cv_ints(int_prim_rcoul_pvxp)%prim, sa_3c_ints(:, int_sa_rcoul_pvxp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif
    if(integralpar_3c_r2_pvxp)then
      call symm_adapt_3cv( prim_3cv_ints(int_prim_r2_pvxp)%prim, sa_3c_ints(:, int_sa_r2_pvxp), &
                       ea2,L2,lsymadp(ua2),ea1,L1,lsymadp(ua1),weight)
    endif

    if(is_on(op_FitTrafo))then
       !
       ! Alpha * P needs to be calculated
       !
       call symm_adapt_2cv( prim_2cv_ints(int_prim_pvec)%prim, sa_alphap(:, 1),&
            & ea2, L2, lsymadpL(ua2),& !<<< L
            & ea1, L1, lsymadpS(ua1),& !<<< S
            & weight)

       call symm_adapt_2cv( prim_2cv_ints(int_prim_pvec)%prim, sa_alphap(:, 2),&
            & ea2, L2, lsymadpS(ua2),& !<<< S
            & ea1, L1, lsymadpL(ua1),& !<<< L
            & weight)
    endif
  end subroutine symadapt_add_totalsym_spor

#ifdef WITH_RESPONSE
  subroutine symadapt_add_nottotalsym_co(la, lb, i_ea1, i_ea2, contr)
    use clebsch_gordan,     only: cg =>cg_eliminated, prod_bas
    use ch_response_module, only: dimension_of_fit_ch
    use debug
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind), intent(in) :: la, lb, i_ea1, i_ea2
    real(kind=r8_kind),    intent(in) :: contr(:,:,:,:,:)


    real(kind=r8_kind), pointer, dimension(:,:,:,:,:) :: &
         saint_3c_co_resp
    type(unique_atom_partner_type), pointer :: sap_a, sap_b
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_ir_c, i_ir_b, i_ir_a
    integer(kind=i4_kind) :: i_if1, i_if2, &
         i_cf1, i_cf2, m1, m2, mult_irrep

    real(kind=r8_kind)    :: coef1, coef, coef2, coeffcg

    integer(kind=i4_kind) :: i_pa_a, i_pa_b, i_pa_c

    integer(kind=i4_kind) :: start

    integer(kind=i4_kind) :: n_ir, n_ua, dim
    integer(kind=i4_kind) :: n_pa_a, n_pa_b, n_pa_c, imult

    type(prod_bas),pointer  :: pcg
    integer(i4_kind) :: off(20)

    !------------ Executable code ------------------------------
    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))

    n_ir = symmetry_data_n_irreps()
    ASSERT(n_ir<=20)
    n_ua = N_unique_atoms

    !       _________________________________
    ! irr   |   irr1  |    irr2   |   irr3  |
    ! pa    |   pa1   | pa1 | pa2 |   pa1   |
    ! ff    | 1 2 3 4 | 1 2 | 1 2 | 1 2 3 4 |
    ! indx :  1 2 3 4   5 6   7 8   9 0 1 2
    !                 ^^^ off(irr2) returns the first index of irr2
    start = 1
    do i_ir_c = 1,n_ir
       off(i_ir_c) = start
       start = start + dimension_of_fit_ch(i_ir_c) * symmetry_data_n_partners(i_ir_c)
    enddo

    irrep_c_: do i_ir_c = 1, n_ir
       n_pa_c = symmetry_data_n_partners(i_ir_c)
       if (dimension_of_fit_ch(i_ir_c) == 0) cycle
       irrep_a_: do i_ir_a = 1, n_ir
          n_pa_a = symmetry_data_n_partners(i_ir_a)
          sap_a => ua1%symadapt_partner(i_ir_a, la)
          if (sap_a%N_independent_fcts == 0) cycle
          irrep_b_: do i_ir_b = 1, n_ir
             n_pa_b = symmetry_data_n_partners(i_ir_b)
             sap_b => ua2%symadapt_partner(i_ir_b, lb)
             if (sap_b%N_independent_fcts == 0) cycle

             mult_irrep = cg(i_ir_c,i_ir_a,i_ir_b)%mult

!!$             if (mult_irrep .eq. 0) cycle

             i_mlt_: do imult = 1, mult_irrep

                pcg              => cg(i_ir_c,i_ir_a,i_ir_b)%sub(imult)
                saint_3c_co_resp => symadapt_int_3c_co_resp(i_ir_c,i_ir_a,i_ir_b)%mult(imult)%int

                pa_c_: do i_pa_c = 1, n_pa_c
                   pa_a_: do i_pa_a = 1, n_pa_a
                      pa_b_: do i_pa_b = 1, n_pa_b

                         ind_fct_1: do i_if1 = 1, sap_a%N_independent_fcts
                            sat1 => sap_a%sa_int(i_ea1,i_if1,i_pa_a)
                            ind_fct_2: do i_if2 = 1, sap_b%N_independent_fcts
                               sat2 => sap_b%sa_int(i_ea2,i_if2,i_pa_b)
                               m_sum_1: do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1 = sat1%c(i_cf1)
                                  m_sum_2: do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2 = sat2%c(i_cf2)
                                     coef  = coef1 * coef2

                                     dim      = dimension_of_fit_ch(i_ir_c)
                                     start    = off(i_ir_c)
                                     start    = start + dim * (i_pa_c - 1)

                                     coeffcg = pcg%c(i_pa_c,i_pa_a,i_pa_b)
                                     saint_3c_co_resp(:,:,:,i_if2,i_if1) =  saint_3c_co_resp(:,:,:,i_if2,i_if1) &
                                          + coef &
                                          * pcg%c(i_pa_c,i_pa_a,i_pa_b) &
                                          * contr(:,:,start:(start+dim)-1,m2,m1) &
                                          / symmetry_data_n_partners(i_ir_c)

                                  end do m_sum_2
                               end do m_sum_1
                            end do ind_fct_2
                         end do ind_fct_1

                      end do pa_b_
                   end do pa_a_
                end do pa_c_

             end do i_mlt_

          enddo irrep_b_
       end do irrep_a_
    enddo irrep_c_

    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
  end subroutine symadapt_add_nottotalsym_co
#endif

  !**************************************************************

  ! RENORMALIZATION STAFF DELETED (AM)

  !**************************************************************
  subroutine sum_up_gradient()
    ! Purpose: Sums up the contributions of one quadrupel to the final
    !          gradient. For this purpose the necessary parts of density and
    !          energy weighted density matrix have to built first
    use orbitalprojection_module
    use occupied_levels_module
    use eigen_data_module
    use symmetry_data_module ! description of irreps
    use occupied_levels_module
    use solv_cavity_module, only: grad_solv_totsym,to_calc_grads
    use elec_static_field_module, only: totalsym_field !!!!!!!!!!!
    use solv_electrostat_module, only: Q_grad
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: gradient_dip_totalsym,gradient_quad_totalsym,gradient_oct_totalsym
    use point_dqo_module, only: torque_dip_totalsym,torque_quad_totalsym,torque_oct_totalsym
    use point_dqo_module, only: gradient_rep_totalsym
    use induced_dipoles_module, only: grad_idip_totalsym
    use induced_dipoles_module, only: torque_idip_totalsym
    use calc3c_switches
    use cpksdervs_matrices
    use virtual_levels_module
    use error_module
#ifdef WITH_EXPERIMENTAL
    use cpks_common, only: cpks_p1, cpks_w1
#endif
   implicit none
   integer(kind=i4_kind) :: i,j,i_ir,i_spin,l_bound_1,l_bound_2,&
         u_bound_1,u_bound_2,i_occ,n_spin,nu_dim,mu_dim,alloc_stat,&
         n_if1,n_if2,n_exp1,n_exp2,s,n_spin_dens,&
         ima,imb,ind_a,ind_b,dim_a,dim_b,off,dim,n_grads
!  logical :: model_density, split_gradients, use_spin_dens, &
   logical :: use_spin_dens, &
              moving_a, moving_b
    real(kind=r8_kind) :: fact
    real(kind=r8_kind), parameter::two=2.0_r8_kind
    real(kind=r8_kind), pointer :: pmunu(:,:) ! density matrix
    real(kind=r8_kind), pointer :: smunu(:,:) ! spin-density matrix
    real(kind=r8_kind), pointer :: wmunu(:,:) ! energy-weighted density matrix
#ifndef NO_SECOND_INTEGRAL_RUN
    real(kind=r8_kind),allocatable :: p1munu(:,:,:),w1munu(:,:)
#ifndef WITH_EXPERIMENTAL
    real(kind=r8_kind),allocatable :: p1munu_temp_occ(:,:),p1munu_temp_vir(:,:), &
                                      w1munu_temp_occ(:,:)
#endif
   integer(kind=i4_kind) :: j_occ!,i_vir,j_vir
#endif /* NO_SECOND_INTEGRAL_RUN */
    real(kind=r8_kind),pointer::eigv_nu(:,:),eigv_mu(:,:)!,eigenvalues(:),occ(:)

    integer(kind=i4_kind):: spin_n_occ,spin_n_vir,eig_dim
    type(symadapt_totsym_2c_int_type), allocatable :: pot_buff(:) !!!!!!!!!!!!!!!AS
!    real(kind=r8_kind):: fixed_occ(2)=(/0.01923418346262114,-0.644173682947655/)
!    real(kind=r8_kind):: fixed_vir(2)=(/-0.654386635734533,  0.116756537457598/)
!    real(kind=r8_kind):: fixed_occ(2)=(/0.00022583688826112,-0.640512045259315/)
!    real(kind=r8_kind):: fixed_vir(2)=(/-0.654669207950957, 0.1354110140245738/)
!   real(kind=r8_kind):: fixed_occ(2)=(/2.258368913317099E-004, -0.640512049198151/)
!   real(kind=r8_kind):: fixed_vir(2)=(/-0.654669225215912, 0.135411009192467/)

    integer(i4_kind) :: i_grad,j_grad
    real(kind=r8_kind), allocatable:: occao(:,:),s1_weighted(:,:)

    !------------ Executable code ------------------------------
    FPP_TIMER_START(sumup)

!   if(cpks_fixed_coeff) then
!    eigvec_occ(1)%m(1:2,1,1)=fixed_occ
!    eigvec_vir(1)%m(1:2,1,1)=fixed_vir
!   endif

    ima = unique_atoms(quadrupel%ua1)%moving_atom
    imb = unique_atoms(quadrupel%ua2)%moving_atom
    moving_a = ima > 0
    moving_b = imb > 0
    if (moving_a) then
      ind_a = gradient_index(ima  ) - 1
      dim_a = gradient_index(ima+1) - ind_a - 1
    else
      dim_a = 0
    endif
    if (moving_b) then
      ind_b = gradient_index(imb  ) - 1
      dim_b = gradient_index(imb+1) - ind_b - 1
    else
      dim_b = 0
    endif

    n_spin=ssym%n_spin
!   split_gradients = options_split_gradients()
!   model_density = options_xcmode() == xcmode_model_density .or. &
!   options_xcmode() == xcmode_extended_mda

    n_spin_dens = 1
    if (model_density) n_spin_dens = n_spin
    use_spin_dens = n_spin_dens > 1

    if(.not.integralpar_pot_for_secderiv) then
       if(      .not. integralpar_solv_grad      &
          .and. .not. integralpar_Q_solv_grad    &
#ifdef WITH_EFP
          .and. .not. integralpar_2cob_field_efp &
#endif
         ) then
          if(size(symadapt_int_3cob_grad,1) ==0) goto 999 ! clean up and exit
       else if( integralpar_solv_grad .and. integralpar_cpksdervs ) then
          if(size(symadapt_int_3cob_solv_grad,1) ==0) goto 999 ! clean up and exit
       endif
    end if

        irrs: do i_ir=1,ssym%n_irrep
              eig_dim=size(eigvec(i_ir)%m,1)

       !
       ! Build parts of density matrix pmunu and energy weigthed
       ! density matrices wmunu which are necesarry for this quadrupel
       !
       ! quad_pmat and quad_wmat were allocated and computed by
       ! integral_calc_quad_densmat() in integral_calc_quad_module:
       !
ASSERT(allocated(quad_pmat(i_ir,1)%m))

       !
       ! Section of the density, energy-density, and spin-density
       ! matrices corresponding to shells pair < 2 | ... | 1 >:
       !
       pmunu => quad_pmat(i_ir, 1)%m(:, :)
       wmunu => quad_wmat(i_ir)%m(:, :)
       if( use_spin_dens )then
         smunu => quad_pmat(i_ir, 2)%m(:, :)
       endif

       !
       ! The latter are not defined, in general, but only the product
       ! counts:
       !
       ! mu_dim = n_exp2 * n_if2
       ! nu_dim = n_exp1 * n_if1
       !
       mu_dim = size(pmunu, 1)
       nu_dim = size(pmunu, 2)

       !
       ! Since not all atomic shells contribute to all irreps,
       ! the number of independent functions n_if1 or n_if2
       ! may turn up to be zero. Then the corresponding section
       ! of the matrices is zero-sized:
       !
       !       (0 x 0) or (nu_dim x 0) or (0 x mu_dim)
       !
       ! This may cause out-of bouns errors in calls to "dgemm".
       ! These out-of-bounds access maybe harmless, since the extents
       ! are zero, but they may cause crashes if bound checking is
       ! enabled.
       !
       ! There are no apparent counters running inside of the
       ! loop over irreps, so if a batch is empty, just cycle:
       if( nu_dim * mu_dim == 0 ) CYCLE irrs !!! *** CYCLE POINT *** !!!

       !
       ! moved from within the branchy "if" at the beginning of the loop,
       ! this is intermediate storage, will be deallocated before going into
       ! next irrep:
       !
       if ( integralpar_pot_for_secderiv ) then
          ! FIXME: lay out this copying into a subprogram, these
          !        varaiables are only used in this branch:
          n_exp2 = size(symadapt_int_2cob_poten(i_ir)%int, 1)
          n_exp1 = size(symadapt_int_2cob_poten(i_ir)%int, 2)
          n_if2  = size(symadapt_int_2cob_poten(i_ir)%int, 4)
          n_if1  = size(symadapt_int_2cob_poten(i_ir)%int, 5)

          n_grads=size(Q_grad,1)
          allocate(pot_buff(n_grads),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

          do i=1,n_grads
             allocate(pot_buff(i)%int(n_exp2, n_exp1, n_if2, n_if1), stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             pot_buff(i)%int=0.0_r8_kind
          end do
          do i=1,n_grads
             do j=1,to_calc_grads%n_points
                pot_buff(i)%int=pot_buff(i)%int-Q_grad(i,j)* &
                     symadapt_int_2cob_poten(i_ir)%int(:,:,j,:,:)
             end do
          end do
       endif

       !
       ! These are offsets into the subblocks of the full matrices
       ! (e.g. Hamiltonian or Densmat) corresponding to the current
       ! irrep i_ir:
       !
       l_bound_1 = orbitalprojection_ob(i_ir, quadrupel%l1, quadrupel%ua1)
       u_bound_1 = l_bound_1 + nu_dim - 1
       l_bound_2 = orbitalprojection_ob(i_ir, quadrupel%l2, quadrupel%ua2)
       u_bound_2 = l_bound_2 + mu_dim - 1

!     print*, 'quad:', quadrupel%ua1,quadrupel%l1,quadrupel%ua2,quadrupel%l2,diagonal
!     print*, l_bound_1,u_bound_1,l_bound_2,u_bound_2,' l1 u1 l2 u2'

#ifdef WITH_SECDER
  cpks_mats: if(integralpar_cpksdervs) then
      !!! calculates cpks matrixes and their contributions to dervs
FPP_TIMER_START(t_dervs_quad_sums)

     DPRINT MyID//"integralpar_cpks_contribs sum eigval", integralpar_cpks_contribs

          impl_dervs: if(integralpar_cpks_contribs) then

          impl_dervs_occ: if(size(eigvec_occ(i_ir)%m,2).gt.0) then

#ifndef NO_SECOND_INTEGRAL_RUN
          ![[=== Contributions due to PX * FY - WX * SY ===========
          !  Note that
          !  a) to compute energy-weighted gradient WX
          !     we take %h1,, %h1ai (explicit gradient of the Fock matrix FX)
          !     from CPKS storage computed in a previous grad+secder run, but
          !  b) we multiply density-matrix gradient PX with the
          !     Fock-matrix gradient FY computed ``on the fly''
          !     in a second gradient run!
          !  c) The latter must be the ``full'' derivative,
          !     otherwise I dont see a point. AM

          allocate( p1munu(mu_dim,nu_dim,n_spin_dens), &
                    w1munu(mu_dim,nu_dim), stat=cpksalloc(107))
          ASSERT(cpksalloc(107).eq.0)

          p1_grads: do i_grad=1,size(cpks,1)
#ifdef WITH_EXPERIMENTAL
       ! Copy the density matrix gradient:
       p1munu(:,:,1) =   cpks_p1(i_ir,i_grad)%m(l_bound_2:u_bound_2,l_bound_1:u_bound_1,1)
       w1munu(:,:)   = - cpks_w1(i_ir,i_grad)%m(l_bound_2:u_bound_2,l_bound_1:u_bound_1)
       ! FIXME: note minus for W1!
       if(.not.diagonal) then
         p1munu=p1munu*2
         w1munu=w1munu*2
       endif
#else
          w1munu=0.0_r8_kind
          p1munu=0.0_r8_kind

          cpks_mats_spins: do i_spin=1,n_spin
                  s = min(i_spin,n_spin_dens) ! for non MDA case n_spin_dens and s = 1
            spin_n_occ=size(cpks(i_grad,i_ir,i_spin)%s1,1)
            spin_n_vir=size(cpks(i_grad,i_ir,i_spin)%HBH,2)

          allocate( p1munu_temp_occ(spin_n_occ,nu_dim), &
                    p1munu_temp_vir(spin_n_vir,nu_dim), &
                    w1munu_temp_occ(spin_n_occ,nu_dim), &
                    stat=cpksalloc(136))
          ASSERT(cpksalloc(136).eq.0)

FPP_TIMER_START(t_p1_sumup)
          ! hole corrected spks%s1 dependent contrib
#if 0
          p1munu_temp_occ= matmul(cpks(i_grad,i_ir,i_spin)%s1, &
                                  transpose(eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)))
#else
          call dgemm('n','t',spin_n_occ,nu_dim,spin_n_occ, &
                      1.0_r8_kind,cpks(i_grad,i_ir,i_spin)%s1,spin_n_occ, &
                                  eigvec(i_ir)%m(l_bound_1,1,i_spin),eig_dim, &
                      0.0_r8_kind,p1munu_temp_occ,spin_n_occ)
#endif

#if 0
          p1munu(:,:,s)=p1munu(:,:,s) - &
                 matmul( eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin), &
                         p1munu_temp_occ)
#else
         if(spin_n_occ*size(p1munu(:,:,s)).gt.0) &
          call dgemm('n','n',mu_dim,nu_dim,spin_n_occ,  -1.0_r8_kind, &
                      eigvec(i_ir)%m(l_bound_2,1,i_spin),eig_dim, &
                      p1munu_temp_occ,spin_n_occ, &
                      1.0_r8_kind,p1munu(1,1,s),mu_dim)
#endif

!          p1vir: if(size(cpks3c(i_ir,i_spin)%eigvec_vir_and_holes,2).gt.0) then
          p1vir: if(spin_n_vir.gt.0) then

#if 0
!          p1munu_temp_occ= matmul( cpks(i_grad,i_ir,i_spin)%HBH(:,:), &
!               transpose(cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_1:u_bound_1,:)))
          p1munu_temp_occ= matmul( cpks(i_grad,i_ir,i_spin)%HBH(:,:), &
               transpose(eigvec(i_ir)%m(l_bound_1:u_bound_1,1+eig_dim-spin_n_vir:,i_spin)))
#else
          if(spin_n_vir*size(p1munu_temp_occ).gt.0) &
          call dgemm('n','t',spin_n_occ,nu_dim,spin_n_vir, 1.0_r8_kind, &
                     cpks(i_grad,i_ir,i_spin)%HBH(:,:),spin_n_occ, &
                     eigvec(i_ir)%m(l_bound_1,1+eig_dim-spin_n_vir,i_spin), eig_dim,  &
                     0.0_r8_kind,p1munu_temp_occ,spin_n_occ )
#endif
#if 0
          p1munu(:,:,s)=p1munu(:,:,s) + &
                 matmul( eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin), &
                         p1munu_temp_occ)
#else
          if(spin_n_occ*size(p1munu(:,:,s)).gt.0) &
          call dgemm('n','n',mu_dim,nu_dim,spin_n_occ, 1.0_r8_kind, &
                     eigvec(i_ir)%m(l_bound_2,1,i_spin), eig_dim, &
                     p1munu_temp_occ, spin_n_occ, &
                     1.0_r8_kind,p1munu(:,:,s),mu_dim)
#endif

!          p1munu_temp_vir= matmul( transpose(cpks(i_grad,i_ir,i_spin)%HBH), &
!                                   transpose(eigvec_occ(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)))
!          p1munu(:,:,s)=p1munu(:,:,s) + &
!                     matmul( cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_2:u_bound_2,:), &
!                             p1munu_temp_vir)
#if 0
          p1munu_temp_vir= matmul( transpose(cpks(i_grad,i_ir,i_spin)%HBH), &
                                   transpose(eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)))
#else
          if(spin_n_occ*size(p1munu_temp_vir).gt.0) &
          call  dgemm('t','t',spin_n_vir,nu_dim,spin_n_occ, 1.0_r8_kind, &
                      cpks(i_grad,i_ir,i_spin)%HBH,spin_n_occ, &
                      eigvec(i_ir)%m(l_bound_1,1,i_spin),eig_dim, &
                      0.0_r8_kind,p1munu_temp_vir,spin_n_vir)
#endif
#if 0
          p1munu(:,:,s)=p1munu(:,:,s) + &
           matmul( eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin), &
                   p1munu_temp_vir)
#else
        if(spin_n_vir*size(p1munu(:,:,s)).gt.0) &
        call dgemm('n','n',mu_dim,nu_dim,spin_n_vir, 1.0_r8_kind, &
                   eigvec(i_ir)%m(l_bound_2,1+eig_dim-spin_n_vir,i_spin),eig_dim, &
                   p1munu_temp_vir,spin_n_vir, &
                   1.0_r8_kind,p1munu(:,:,s),mu_dim)
#endif
          endif p1vir
FPP_TIMER_STOP(t_p1_sumup)

!          if(diagonal) then
!           p1munu(:,:,s)=p1munu(:,:,s)*2.0_r8_kind
!          else
!           p1munu(:,:,s)=p1munu(:,:,s)*4.0_r8_kind
!          endif

!          print*,sum(pmunu(:,:,s)),sum(p1munu(:,:,s)),i_grad,'pmunu p1munu',diagonal

!          if(.true.) then   ! 1st of 3 implicit contribs
!             call trace( ssym%partner(i_ir)*p1munu(:,:,s), &
!                         symadapt_int_3cob_grad(:,i_ir), &
!                         dervs_totalsym(:,i_grad) )             ! 1st of 5 dervs_totalsym contribs
!          else
!           print*,'p1 cpks contrib to totalsym_dervs is off'
!          endif

FPP_TIMER_START(t_w1_sumup)
#if 0
           w1munu_temp_occ(:,:)=matmul( cpks(i_grad,i_ir,i_spin)%h1, &
                                transpose(eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)))
#else
           if(spin_n_occ*size(w1munu_temp_occ).gt.0) &
           call dgemm('n','t',spin_n_occ,nu_dim,spin_n_occ, 1.0_r8_kind, &
                      cpks(i_grad,i_ir,i_spin)%h1,spin_n_occ,&
                      eigvec(i_ir)%m(l_bound_1,1,i_spin),eig_dim,&
                      0.0_r8_kind,w1munu_temp_occ,spin_n_occ)
#endif
#if 0
           w1munu(:,:)=w1munu(:,:)-matmul(eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin), &
                                          w1munu_temp_occ)
#else
          if(spin_n_occ*size(w1munu).gt.0) &
          call dgemm('n','n',mu_dim,nu_dim,spin_n_occ, -1.0_r8_kind, &
                     eigvec(i_ir)%m(l_bound_2,1,i_spin),eig_dim, &
                     w1munu_temp_occ,spin_n_occ, &
                     1.0_r8_kind,w1munu,mu_dim)
#endif

#if 0
          w1_s1: do i_occ=1,spin_n_occ
          do j_occ=1,spin_n_occ
          do  nu=1,nu_dim
          do  mu=1,mu_dim
           w1munu(mu,nu)=w1munu(mu,nu)+ &
                 (eigval_occ(i_ir)%m(i_occ,i_spin)+eigval_occ(i_ir)%m(j_occ,i_spin))* &
                   eigvec(i_ir)%m(mu+l_bound_2-1,i_occ,i_spin)* &
                     eigvec(i_ir)%m(nu+l_bound_1-1,j_occ,i_spin) &
                       *cpks(i_grad,i_ir,i_spin)%s1(i_occ,j_occ)
           enddo
           enddo
           enddo
           enddo w1_s1
#else
     eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
     eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
     allocate(s1_weighted(spin_n_occ,spin_n_occ),occao(spin_n_occ,size(eigv_mu,1)),stat=cpksalloc(157))
     ASSERT(cpksalloc(157).eq.0)
     do i_occ=1,spin_n_occ
     do j_occ=1,spin_n_occ
     s1_weighted(i_occ,j_occ)=cpks(i_grad,i_ir,i_spin)%s1(i_occ,j_occ)* &
              (eigval_occ(i_ir)%m(i_occ,i_spin)+eigval_occ(i_ir)%m(j_occ,i_spin))
     enddo
     enddo
#if 0
     occao=matmul(s1_weighted,transpose(eigv_mu))
#else
     if(spin_n_occ*size(occao).gt.0) &
     call dgemm('n','t',spin_n_occ,mu_dim,spin_n_occ, 1.0_r8_kind, &
                s1_weighted,spin_n_occ, &
                eigv_mu,mu_dim, &
                0.0_r8_kind, occao, spin_n_occ) ! FIXME: temp copy
#endif
#if 0
     w1munu=w1munu+matmul(transpose(occao),transpose(eigv_nu))
#else
     if(spin_n_occ*size(w1munu).gt.0) &
     call dgemm('t','t',mu_dim,nu_dim,spin_n_occ, 1.0_r8_kind, &
                occao, spin_n_occ, eigv_nu, nu_dim, &
                1.0_r8_kind, w1munu, mu_dim) ! FIXME: temp copy
#endif

     deallocate(s1_weighted,occao,stat=cpksalloc(157))
     ASSERT(cpksalloc(157).eq.0)
     cpksalloc(157)=1
#endif

   vir_w1: if(spin_n_vir.gt.0) then
!                eigv_mu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_2:u_bound_2,:)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin)
#if 0
   occvir_w1:  do i_occ=1,size(eigv_nu,2)
               do j_vir=1,size(eigv_mu,2)
          do  nu=1,nu_dim
          do  mu=1,mu_dim
           w1munu(mu,nu)=w1munu(mu,nu)-eigval_occ(i_ir)%m(i_occ,i_spin)* &
             eigv_mu(mu,j_vir)*eigv_nu(nu,i_occ) &
                    *cpks(i_grad,i_ir,i_spin)%HBH(i_occ,j_vir)
!                    *cpks(i_grad,i_ir,1)%b(i_occ,j_vir,n_cpks+1)*ssym%partner(i_ir)
           enddo
           enddo


     enddo! j mo loop
    enddo occvir_w1
#else
  allocate(occao(size(eigv_nu,2),size(eigv_mu,1)),stat=cpksalloc(157))
  ASSERT(cpksalloc(157).eq.0)
#if 0
      occao=matmul(cpks(i_grad,i_ir,i_spin)%HBH,transpose(eigv_mu))
#else
    if(spin_n_vir*size(occao).gt.0) &
    call dgemm('n','t',spin_n_occ,mu_dim,spin_n_vir, 1.0_r8_kind, &
                cpks(i_grad,i_ir,i_spin)%HBH,spin_n_occ, &
                eigv_mu, mu_dim, &
                0.0_r8_kind, occao, spin_n_occ) ! FIXME: temp copy
#endif
#if 0
   do nu=1,size(eigv_nu,1)
    do i_occ=1,size(eigv_nu,2)
     do mu=1,size(eigv_mu,1)
      w1munu(mu,nu)=w1munu(mu,nu)- &
       occao(i_occ,mu)*eigval_occ(i_ir)%m(i_occ,i_spin)*eigv_nu(nu,i_occ)
     enddo
    enddo
   enddo
#else
   do i_occ=1,spin_n_occ
    occao(i_occ,:)=occao(i_occ,:)*eigval_occ(i_ir)%m(i_occ,i_spin)
   enddo
   if(spin_n_occ*size(w1munu).gt.0) &
   call dgemm('t','t',mu_dim,nu_dim,spin_n_occ, -1.0_r8_kind, &
               occao,spin_n_occ,  eigv_nu,nu_dim, &
               1.0_r8_kind, w1munu, mu_dim ) ! FIXME: temp copy
#endif
  deallocate(occao,stat=cpksalloc(157))
  ASSERT(cpksalloc(157).eq.0)
  cpksalloc(157)=1
#endif

!                eigv_nu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_1:u_bound_1,:)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,1+eig_dim-spin_n_vir:,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
#if 0
    virocc_w1:  do i_vir=1,size(eigv_nu,2)
                        do j_occ=1,size(eigv_mu,2)

          do  nu=1,nu_dim
          do  mu=1,mu_dim
           w1munu(mu,nu)=w1munu(mu,nu)- &
                   eigval_occ(i_ir)%m(j_occ,i_spin)* &
                   eigv_mu(mu,j_occ)*eigv_nu(nu,i_vir) &
                    *cpks(i_grad,i_ir,i_spin)%HBH(j_occ,i_vir)
          enddo
          enddo


     enddo! j mo loop
    enddo virocc_w1
#else
  allocate(occao(size(eigv_mu,2),size(eigv_nu,1)),stat=cpksalloc(157))
  ASSERT(cpksalloc(157).eq.0)
#if 0
    occao=matmul(cpks(i_grad,i_ir,i_spin)%HBH,transpose(eigv_nu))
#else
    if(spin_n_vir*size(occao).gt.0) &
    call dgemm('n','t',spin_n_occ,nu_dim,spin_n_vir, 1.0_r8_kind, &
                cpks(i_grad,i_ir,i_spin)%HBH, spin_n_occ, &
                eigv_nu, nu_dim, &
                0.0_r8_kind, occao,spin_n_occ) ! FIXME: temp copy
#endif
#if 0
   do nu=1,size(eigv_nu,1)
    do j_occ=1,size(eigv_mu,2)
     do mu=1,size(eigv_mu,1)
      w1munu(mu,nu)=w1munu(mu,nu)- &
       occao(j_occ,nu)*eigval_occ(i_ir)%m(j_occ,i_spin)*eigv_mu(mu,j_occ)
     enddo
    enddo
   enddo
#else
   do j_occ=1,spin_n_occ
    occao(j_occ,:)=occao(j_occ,:)*eigval_occ(i_ir)%m(j_occ,i_spin)
   enddo
   if(spin_n_occ*size(w1munu).gt.0) &
   call dgemm('n','n', mu_dim,nu_dim, spin_n_occ, -1.0_r8_kind, &
               eigv_mu, mu_dim,  occao,spin_n_occ, &
               1.0_r8_kind, w1munu, mu_dim) ! FIXME: temp copy
#endif
  deallocate(occao,stat=cpksalloc(157))
  ASSERT(cpksalloc(157).eq.0)
  cpksalloc(157)=1
#endif

   endif vir_w1
          deallocate( p1munu_temp_occ, p1munu_temp_vir, w1munu_temp_occ,  &
                      stat=cpksalloc(136))
                      ASSERT(cpksalloc(136).eq.0)
                      cpksalloc(136)=1

  enddo cpks_mats_spins

          if(diagonal) then
           p1munu=p1munu*(3-n_spin) * ssym%partner(i_ir)
           w1munu=w1munu*(3-n_spin) * ssym%partner(i_ir)
          else
           p1munu=p1munu*2.0_r8_kind*(3-n_spin) * ssym%partner(i_ir)
           w1munu=w1munu*2.0_r8_kind*(3-n_spin) * ssym%partner(i_ir)
          endif
FPP_TIMER_STOP(t_w1_sumup)
#endif /* of ifdef WITH_EXPERIMENTAL */

    no_pot: if(.not.integralpar_pot_for_secderiv) then !!!!!!!!!!!!!!!!!
      no_solv: if(.not.integralpar_solv_grad .and. .not.integralpar_Q_solv_grad) then     !!!!!!!!!!!!!!!!!!
#define IMPL_CONTRIBS_TO_DERVS
          ![[=== Add terms due to P^i * F^j ======================
#ifdef IMPL_CONTRIBS_TO_DERVS
             call trace( p1munu(:,:,1), &
                         symadapt_int_3cob_grad(:,i_ir), &
                         dervs_totalsym(:,i_grad) )             ! 1st of 5 dervs_totalsym contribs
#endif

             ![[=== Add terms due to 1st derivatives of rel. ham.:
             if( integralpar_relativistic )then
               ! Hrel = Trel + Vrel, was saved by cpks_h1_store before CPKS
               ! in %hr1
               do j_grad=1,size(cpks,1)
                 ! dervs(j,i) += trace( P^i * Hrel^j ):
                 dervs_totalsym(j_grad,i_grad) = dervs_totalsym(j_grad,i_grad) &
                   + SUM(   p1munu(:,:,1)                                        &
                          * cpks(j_grad,i_ir,1)%hr1(l_bound_2:u_bound_2,l_bound_1:u_bound_1) &
                        )
               enddo
             endif
             !]]==================================================
          !]]=====================================================


         w1_ol_grad: if(totsym_w_ol_grad) then !(3) 2nd of 3 implicit contibs
             if (moving_a) then
                 call trace( w1munu(:,:), &
                      symadapt_int_2cob_ol_grad(1:dim_a,i_ir), &
                      dervs_totalsym(ind_a+1:ind_a+dim_a,i_grad) )
             endif

             if( moving_b) then                                  ! 1st of 2  b contrib
                 call trace( w1munu(:,:), &
                      symadapt_int_2cob_ol_grad(dim_a+1:dim_a+dim_b,i_ir), &
                      dervs_totalsym(ind_b+1:ind_b+dim_b,i_grad) )
             endif
         endif w1_ol_grad

      else if(integralpar_solv_grad) then !no_solv
       ASSERT(integralpar_cpksdervs )
       !!! here two  implicite contribs to dervs follow
         call trace( p1munu(:,:,1), &
              symadapt_int_3cob_solv_grad(:,i_ir), &  !1st impl solv contrib to dervs
              dervs_totalsym(:,i_grad) )                  !correct becouse P1 is OK
      end if no_solv

    else no_pot !i.e. potential for 2nd dervs
       call trace( p1munu(:,:,1), &
           pot_buff(:),dervs_totalsym(:,i_grad) )     !2nd impl solv contrib
    end if no_pot !!!!!!!!!!!!!AS

   enddo p1_grads


          deallocate(p1munu, w1munu,stat=cpksalloc(107))
          ASSERT(cpksalloc(107).eq.0)
          cpksalloc(107)=1
          !]]=== EOF Contributions due to PX * FY - WX * SY ========
#endif /* NO_SECOND_INTEGRAL_RUN */


#ifndef no_cpks_coul_grads
        orbgr_fitgr: if(.not.integralpar_solv_grad .and. &
             .not.integralpar_Q_solv_grad .and. &
             .not.integralpar_pot_for_secderiv) then ! 2nd of 5 dervs_totalsym contribs
          cpks_gradient_totalsym=0.0_r8_kind ! no charge overlap contribs
          call cpks_trace( pmunu(:,:), &
                           symadapt_int_cpks_coul_grad(:,i_ir), &   !(2nd use) impl contribs
                           cpks_gradient_totalsym )            ! 1st contrib
          ! this cpks_gradient_totalsym is not used in cpks calcs
          ! but contribute to dervs_totalsym via product with cpks_fitcoeff_grads here

          do i_grad=1,size(cpks_gradient_totalsym,2)
           do j_grad=1,size(cpks_gradient_totalsym,2)
           dervs_totalsym(i_grad,j_grad)=dervs_totalsym(i_grad,j_grad) +  &
           dot_product(cpks_gradient_totalsym(:,i_grad),cpks_fitcoeff_grads(:,j_grad))
           enddo
          enddo
        endif orbgr_fitgr
#endif


       endif impl_dervs_occ

       else impl_dervs ! actually matixes to solve CPKS
                       ! no direct contribs to dervs

             if(diagonal) then
                fact=0.5_r8_kind
             else
                fact=1.0_r8_kind
             endif

   no_cpks_solv: if(.not.integralpar_solv_grad .and. &
          .not.integralpar_Q_solv_grad .and. &
          .not.integralpar_pot_for_secderiv) then     !!!!!!!!!!!!!!!!!!
#ifndef no_cpks_coul_grads
             call cpks_trace( pmunu(:,:), &
                              symadapt_int_cpks_coul_grad(:,i_ir), &   !(1st use) cpks matrixes
                              cpks_gradient_totalsym ) ! other contrib is from charge overlap matrix
#endif

      s1w_spins: do i_spin=1,n_spin
      spin_n_occ=size(cpks(1,i_ir,i_spin)%s1,1)
      spin_n_vir=size(cpks(1,i_ir,i_spin)%h1ai,2)

!      do i_grad=1,size(cpks,1)
!      DPRINT MyID//'sum %s1 start',sum(cpks(i_grad,i_ir,i_spin)%s1),i_grad,i_ir,i_spin
!      enddo

!                eigv_nu=>eigvec_occ(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
!                eigv_mu=>eigvec_occ(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
!            print*,'eigvec_occ'
!            print*,eigvec(i_ir)%m(:,:spin_n_occ,i_spin)
             DPRINT MyID//'s1 eigvec sum',sum(eigvec(i_ir)%m)

            mova: if (moving_a) then
              cpks_ga: do i_grad=1,dim_a
                call MO( eigv_mu, eigv_nu                                 &
                       , symadapt_int_2cob_ol_grad(i_grad+0,i_ir)%int &
                       , cpks(i_grad+ind_a,i_ir,i_spin)%s1                &
                       , fact, "s" )
      enddo cpks_ga
      endif mova

            movb: if (moving_b) then
             cpks_gb: do i_grad=1,dim_b
                call MO( eigv_mu, eigv_nu                                 &
                       , symadapt_int_2cob_ol_grad(i_grad+dim_a,i_ir)%int &
                       , cpks(i_grad+ind_b,i_ir,i_spin)%s1                &
                       , fact, "s" )
      enddo cpks_gb
      endif movb

      h1_cpks_w: do i_grad=1,size(cpks,1)
                !  calculated for explicit calc call
                !  explicit contrib to be augmented with XC part
                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_3cob_grad(i_grad,i_ir)%int &
                       , cpks(i_grad,i_ir,i_spin)%h1             &
                       , fact, "s" )
      enddo h1_cpks_w

#if 1
   h1ai_cpks_g: do i_grad=1,size(cpks,1)

!          cpks(i_grad,i_ir,1)%h1ai=0.0_r8_kind   !!!!!!

!                eigv_nu=>eigvec_occ(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
!                eigv_mu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_2:u_bound_2,:)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_3cob_grad(i_grad,i_ir)%int &
                       , cpks(i_grad,i_ir,i_spin)%h1ai           &    ! (1st-contrib)
                       , fact, "t" )

!                eigv_nu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_1:u_bound_1,:)
!                eigv_mu=>eigvec_occ(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,1+eig_dim-spin_n_vir:,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_3cob_grad(i_grad,i_ir)%int &
                       , cpks(i_grad,i_ir,i_spin)%h1ai           &    !(2nd-contib)
                       , fact, "n" )

    enddo h1ai_cpks_g
#endif

!      do i_grad=1,size(cpks,1)
!       print*, MyID//'sum %h1ai fin',sum(cpks(i_grad,i_ir,1)%h1ai),i_grad,i_ir,&
!               sum(symadapt_int_3cob_grad(i_grad,i_ir)%int), &
!               sum(eigval_vir(i_ir)%m),sum(eigvec_vir(i_ir)%m), &
!               l_bound_1,u_bound_1,l_bound_2,u_bound_2, &
!               sum(eigval(i_ir)%m),sum(eigvec_occ(i_ir)%m)
!      enddo


           if (moving_a) then
           cpks_s1aiga: do i_grad=1,dim_a

!                eigv_nu=>eigvec_occ(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
!                eigv_mu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_2:u_bound_2,:)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_2cob_ol_grad(i_grad+0,i_ir)%int &
                       , cpks(i_grad+ind_a,i_ir,i_spin)%s1ai           &
                       , fact, "t" )

!                eigv_nu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_1:u_bound_1,:)
!                eigv_mu=>eigvec_occ(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,1+eig_dim-spin_n_vir:,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_2cob_ol_grad(i_grad+0,i_ir)%int &
                       , cpks(i_grad+ind_a,i_ir,i_spin)%s1ai           &
                       , fact, "n" )
    enddo cpks_s1aiga
    endif

    if (moving_b) then
     cpks_s1aigb: do i_grad=1,dim_b

!                eigv_nu=>eigvec_occ(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
!                eigv_mu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_2:u_bound_2,:)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_2cob_ol_grad(i_grad+dim_a,i_ir)%int &
                       , cpks(i_grad+ind_b,i_ir,i_spin)%s1ai           &
                       , fact, "t" )

!                eigv_nu=>cpks3c(i_ir,i_spin)%eigvec_vir_and_holes(l_bound_1:u_bound_1,:)
!                eigv_mu=>eigvec_occ(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
                eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,eig_dim-spin_n_vir+1:,i_spin)
                eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)

                call MO( eigv_mu, eigv_nu                        &
                       , symadapt_int_2cob_ol_grad(i_grad+dim_a,i_ir)%int &
                       , cpks(i_grad+ind_b,i_ir,i_spin)%s1ai           &
                       , fact, "n" )
    enddo cpks_s1aigb
    endif

   enddo s1w_spins

  else if(integralpar_solv_grad .or. integralpar_pot_for_secderiv) then !no_cpks_solv
      !!! matrixes for CPKS no direct contributions to dervs
      !!! solvation contribs to h1 and h1ai which are correct
      !!! becouse CPKS solutions are exact and h1 us exact
      s1w_spins_solv: do i_spin=1,n_spin
         spin_n_occ=size(cpks(1,i_ir,i_spin)%s1,1)
         spin_n_vir=size(cpks(1,i_ir,i_spin)%h1ai,2)

         eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
         eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
         h1_cpks_w_solv: do i_grad=1,size(cpks,1)
            if(.not.integralpar_pot_for_secderiv) then
               ASSERT(integralpar_cpksdervs )
               call MO( eigv_mu, eigv_nu                        &
                    , symadapt_int_3cob_solv_grad(i_grad,i_ir)%int &
                    , cpks(i_grad,i_ir,i_spin)%h1             &         ! tested OK
                    , fact, "s" )
            else !.i.e. potential contrib
               call MO( eigv_mu, eigv_nu                        &
                    , pot_buff(i_grad)%int &
                    , cpks(i_grad,i_ir,i_spin)%h1             &         ! tested OK
                    , fact, "s" )
            end if
         enddo h1_cpks_w_solv

         h1ai_cpks_g_solv: do i_grad=1,size(cpks,1)
            eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,:spin_n_occ,i_spin)
            eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,1+eig_dim-spin_n_vir:,i_spin)
            ! for cpks equation one need here one partical part
            if(.not.integralpar_pot_for_secderiv) then
               ASSERT(integralpar_cpksdervs )
               call MO( eigv_mu, eigv_nu                        &
                    , symadapt_int_3cob_solv_grad(i_grad,i_ir)%int &
                    , cpks(i_grad,i_ir,i_spin)%h1ai           &       ! tested OK
                    , fact, "t" )
            else
               call MO( eigv_mu, eigv_nu                        &
                    , pot_buff(i_grad)%int &                          !tested OK
                    , cpks(i_grad,i_ir,i_spin)%h1ai           &
                    , fact, "t" )
            end if
            eigv_nu=>eigvec(i_ir)%m(l_bound_1:u_bound_1,1+eig_dim-spin_n_vir:,i_spin)
            eigv_mu=>eigvec(i_ir)%m(l_bound_2:u_bound_2,:spin_n_occ,i_spin)
            if(.not.integralpar_pot_for_secderiv) then
               ASSERT(integralpar_cpksdervs )
               call MO( eigv_mu, eigv_nu                        &
                    , symadapt_int_3cob_solv_grad(i_grad,i_ir)%int &   !tested OK
                    , cpks(i_grad,i_ir,i_spin)%h1ai           &
                    , fact, "n" )
            else
               call MO( eigv_mu, eigv_nu                        &       !tested OK
                    , pot_buff(i_grad)%int &
                    , cpks(i_grad,i_ir,i_spin)%h1ai           &
                    , fact, "n" )
            end if
         end do h1ai_cpks_g_solv
      end do s1w_spins_solv
   endif no_cpks_solv !!!!!!!!!!!!!!!!!!!!!AS

  endif impl_dervs !else
FPP_TIMER_STOP(t_dervs_quad_sums)
  endif cpks_mats
#endif

       ! building of density matrices finished

#ifdef WITH_EPE
       if(ewpc_n.ne.0.and.calc_cluster_epe_energy)  then !exact density and epe PCs
        call trace( pmunu(:,:),symadapt_int_3cob_epe(:,i_ir),cluster_epe_energy)
       endif
#endif

       ! now sum up contribution for every gradient
       nopot: if(.not.integralpar_pot_for_secderiv ) then
       notsolv: if (.not.integralpar_solv_grad .and. .not.integralpar_Q_solv_grad )then

          spli: if (split_gradients) then
             ! LOAD orbital Pulay corrections
             if (moving_a) then

                ! ADD tr( W d/dRab S )
                call trace( wmunu(:,:),symadapt_int_2cob_ol_grad(1:dim_a,i_ir), &
                     gradient_ob_pulay(ind_a+1:ind_a+dim_a) )

                ! ADD tr( P_tot d/dRab H_RKS )
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_ks_grad(1:dim_a,i_ir), &
                     gradient_ob_pulay(ind_a+1:ind_a+dim_a) )
                if (use_spin_dens) then
                   ! ADD tr( P_spin d/dRab V_xc,spin )
                   off = dim_a + dim_b
                   call trace( smunu(:,:), & ! < xi_i | V_xc,spin | xi_j >
                        symadapt_int_2cob_ks_grad(off+1:off+dim_a,i_ir), &
                        gradient_ob_pulay(ind_a+1:ind_a+dim_a) )
                endif
             endif
             if( moving_b) then
                off = dim_a
                ! ADD tr( W d/dRab S )
                call trace( wmunu(:,:),symadapt_int_2cob_ol_grad(off+1:off+dim_b,i_ir), &
                     gradient_ob_pulay(ind_b+1:ind_b+dim_b) )
                ! ADD tr( P_tot d/dRab H_RKS )
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_ks_grad(off+1:off+dim_b,i_ir), &
                     gradient_ob_pulay(ind_b+1:ind_b+dim_b) )
                if (use_spin_dens) then
                   ! ADD tr( P_spin d/dRab V_xc,spin )
                   off = off + dim_a + dim_b
                   call trace( smunu(:,:), & ! < xi_i | V_xc,spin | xi_j >
                        symadapt_int_2cob_ks_grad(off+1:off+dim_b,i_ir), &
                        gradient_ob_pulay(ind_b+1:ind_b+dim_b) )
                endif
             endif

             ! LOAD charge fit Pulay correction contribution
             if (model_density) then
                ! ADD tr( P_tot d/dRc V_H )
                call trace( pmunu(:,:), &
                     symadapt_int_3cob_coul_grad(:,i_ir), &
                     gradient_ch_pulay(:) )

             else ! standard SCF
                ! ADD tr( P_tot d/dRc V_H )
                call trace( pmunu(:,:), symadapt_int_3cob_grad(:,i_ir), &
                     gradient_ch_pulay(:) )
             endif

             ! LOAD MDA XC-potential Pulay correction contribution
             if (model_density) then
                ! ADD tr( P_tot d/dRc V_X,tot )
                dim = gradient_data_n_gradients
                call trace( pmunu(:,:), &
                     symadapt_int_3cob_grad(1:dim,i_ir), &
                     gradient_mda_vxc(:) )
                if (use_spin_dens) then
                   ! ADD tr( P_spin d/dRc V_X,spin )
                   off = dim
                   call trace( smunu(:,:), &
                        symadapt_int_3cob_grad(off+1:off+dim,i_ir), &
                        gradient_mda_vxc(:) )
                endif
             endif

             ! LOAD Hellman-Feyman force contributions
             ! ADD tr( P_tot d/dRc V_nuc )
             call trace( pmunu(:,:), &
                  symadapt_int_3cob_nuc_grad(:,i_ir), &
                  gradient_totalsym(:) )

          else spli ! i.e. total gradient only

!            pmunu=1.0_r8_kind   ! to check bare functional part

             ! LOAD weighted density matrix contributions

            w_ol_grad: if(totsym_w_ol_grad) then ! (1)
             if (moving_a) then                                  ! 1st of 2  a contrib
                 call trace( wmunu(:,:), &
                      symadapt_int_2cob_ol_grad(1:dim_a,i_ir), &
                      gradient_totalsym(ind_a+1:ind_a+dim_a) )
             endif
             if( moving_b) then                                  ! 1st of 2  b contrib
                 call trace( wmunu(:,:), &
                      symadapt_int_2cob_ol_grad(dim_a+1:dim_a+dim_b,i_ir), &
                      gradient_totalsym(ind_b+1:ind_b+dim_b) )
             endif
             else w_ol_grad
             print*, ' wmunu x ol_grad contribs off'
             endif w_ol_grad

             ! LOAD standard density matrix contributions
             dim = gradient_data_n_gradients

  not_impl_dervs: if(.not.integralpar_cpks_contribs) then
!            do i_grad=1,size(symadapt_int_3cob_grad,1)
!             print*,sum(symadapt_int_3cob_grad(i_grad,i_ir)%int), &
!                      'sum symadapt_int_3cob_grad(i_grad,i_ir)',i_grad,i_ir
!            enddo

#if 0 /* uncomment to cut off fit contribs */
            gradient_totalsym(:)=0.0_r8_kind
             print*,'!!!!!!!!!! fit contribs off'
#endif

             call trace( pmunu(:,:),symadapt_int_3cob_grad(1:dim,i_ir), &
                         gradient_totalsym(:) )     ! this is one of two contribs
#if 0
  print*,'sumup 3cob_grad', gradient_totalsym(1)
#endif



       sc_dervs: if(integralpar_dervs) then
             ! LOAD standard density matrix contributions
             dim = gradient_data_n_gradients

            orb_dervs_contribs: if(totsym_w_ol_grad) then     !(2)  3rd of 5 dervs_totalsym contribs

            ! this is 2nd of two explicit dervs contribs

             if (moving_a) then

                 call trace_dervs( wmunu(:,:), &
                      symadapt_int_2cob_ol_dervs(1:dim_a,1:dim_a,i_ir), &        !!! sum up  (int1)
                      dervs_totalsym(ind_a+1:ind_a+dim_a,ind_a+1:ind_a+dim_a) )  !(1)

              if(moving_b) then
                 call trace_dervs( wmunu(:,:), &
                      symadapt_int_2cob_ol_dervs(1:dim_a,dim_a+1:dim_a+dim_b,i_ir), &
                      dervs_totalsym(ind_a+1:ind_a+dim_a,ind_b+1:ind_b+dim_b) )  !(2)
                 call trace_dervs( wmunu(:,:), &
                      symadapt_int_2cob_ol_dervs(dim_a+1:dim_a+dim_b,1:dim_a,i_ir), &
                      dervs_totalsym(ind_b+1:ind_b+dim_b,ind_a+1:ind_a+dim_a) )  !(3)
              endif

             endif

             if (moving_b) then
              call trace_dervs( wmunu(:,:), &
               symadapt_int_2cob_ol_dervs(dim_a+1:dim_a+dim_b,dim_a+1:dim_a+dim_b,i_ir), &
                      dervs_totalsym(ind_b+1:ind_b+dim_b,ind_b+1:ind_b+dim_b) )  !(4)
             endif

             else orb_dervs_contribs
             print*,'ol_dervs * wmunu contribs off'
            endif orb_dervs_contribs


            ! 1st of two explicite dervs contribs

!           dervs_totalsym(:,:)=0.0_r8_kind  !!!!! uncomment to see individual contribs
             call trace_dervs( pmunu(:,:), &                              !(int1)
                               symadapt_int_coul_dervs(1:dim,1:dim,i_ir), & !4th of 5 dervs_totalsym contribs
                               dervs_totalsym(:,:) )  !(5)                  !actually contains con3c dervs

!             print*,'COUL + H1 SUM  COUL_DERVS (2,2) '

!!             do i_grad=1,size(dervs_totalsym,1)
!!             do k2dr=1,size(dervs_totalsym,1)
!!              print*,sum(symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int), &
!!              dervs_totalsym(i_grad,k2dr),i_grad,k2dr,i_ir
!!             enddo
!!             enddo
!    print*,'dervs_totalsym(1,1) individual',dervs_totalsym(1,1)
!    do i_grad=1,size(dervs_totalsym,1)
!      print*,dervs_totalsym(i_grad,:)
!    write(*,'(10f8.4)') dervs_totalsym(i_grad,:),sum(dervs_totalsym(i_grad,:))
!    enddo
!    write(*,'(10f8.4)') sum(dervs_totalsym(:,:),1)


             endif sc_dervs
          endif not_impl_dervs

             if (use_spin_dens) then
                off = dim
                call trace( smunu(:,:), &
                     symadapt_int_3cob_grad(off+1:off+dim,i_ir), &
                     gradient_totalsym(:) )
             endif

          endif spli

          if(.not.integralpar_cpks_contribs) then
             if(integralpar_2cob_pc_grad) then
                dim=totsym_grad_pc_length

                call trace( pmunu(:,:), &
                     symadapt_int_2cob_pc_grad(1:dim,i_ir), &
                     gradient_pc_totalsym(:))
             end if

             if(integralpar_2cob_X_grad) then
                dim=totsym_grad_dip_length
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_pd_grad(1:dim,i_ir), &
                     gradient_dip_totalsym(:))
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_pd_torq(1:dim,i_ir), &
                     torque_dip_totalsym(:))

                dim=totsym_grad_quad_length
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_pq_grad(1:dim,i_ir), &
                     gradient_quad_totalsym(:))
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_pq_torq(1:dim,i_ir), &
                     torque_quad_totalsym(:))

                dim=totsym_grad_oct_length
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_po_grad(1:dim,i_ir), &
                     gradient_oct_totalsym(:))
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_po_torq(1:dim,i_ir), &
                     torque_oct_totalsym(:))

                dim=totsym_grad_rep_length
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_rc_grad(1:dim,i_ir), &
                     gradient_rep_totalsym(:))
             end if

             if(integralpar_2cob_ipd_grad) then
                dim=totsym_grad_idip_length
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_ipd_grad(1:dim,i_ir), &
                     grad_idip_totalsym(:))
                call trace( pmunu(:,:), &
                     symadapt_int_2cob_ipd_torq(1:dim,i_ir), &
                     torque_idip_totalsym(:))
             end if
            ! here a strage conflict is marked by darcs 2, hope this resolves it
          end if

       else notsolv ! i.e. solvation grads and 2nd dervs
          dim = gradient_data_n_gradients
          no_impl_2d:if(.not.integralpar_cpks_contribs) then
             if( integralpar_solv_grad .and. integralpar_cpksdervs ) then
#if 1 /* p0 solv_grad contrib */
                call trace( pmunu(:,:), &
                     symadapt_int_3cob_solv_grad(1:dim,i_ir), &
                     grad_solv_totsym(:) )
#endif
                if (use_spin_dens) then
                   off = dim
                   call trace( smunu(:,:), &
                        symadapt_int_3cob_solv_grad(off+1:off+dim,i_ir), &
                        grad_solv_totsym(:) )
                endif
                if (with_pc .and. .not.fixed_pc) then
                   dim = totsym_field_length
                   call trace( pmunu(:,:), &
                        symadapt_int_3cob_solv_grad_pc(1:dim,i_ir), &
                        totalsym_field(:) )
                end if
             endif
          end if no_impl_2d
       endif notsolv
       end if nopot

       if(integralpar_pot_for_secderiv) then
          do i=1,n_grads
             deallocate(pot_buff(i)%int,stat=alloc_stat)
             ASSERT(alloc_stat.eq.0)
          end do
          deallocate(pot_buff,stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
       end if
    enddo irrs

999 CONTINUE ! clean up and exit

    FPP_TIMER_STOP(sumup)

  end subroutine sum_up_gradient

  !**************************************************************
  ! This sub was moved (with slight modifications) to
  ! integral_calc_quad_module:
  !
  !   subroutine quad_density_mat
  !
  ! It is exported and called from here.
  !**************************************************************

#ifdef WITH_EFP
  subroutine sum_up_grads_on_efp()
    ! Purpose: Sums up the contributions of one quadrupel to the final
    !          gradient. For this purpose the necessary parts of density and
    !          energy weighted density matrix have to built first
    use symmetry_data_module ! description of irreps
    use elec_static_field_module, only: totalsym_field !!!!!!!!!!!
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: gradient_dip_totalsym,gradient_quad_totalsym,gradient_oct_totalsym
    use point_dqo_module, only: torque_dip_totalsym,torque_quad_totalsym,torque_oct_totalsym
    use point_dqo_module, only: gradient_rep_totalsym
    use induced_dipoles_module, only: grad_idip_totalsym
    use induced_dipoles_module, only: torque_idip_totalsym
    use efp_efp_module, only: qm_efp_energy
    use error_module
    implicit none
    integer(kind=i4_kind) :: i,i_ir,n_spin,n_spin_dens,dim
    logical :: use_spin_dens
    real(kind=r8_kind), pointer :: pmunu(:,:) ! density matrix section

    !------------ Executable code ------------------------------
    FPP_TIMER_START(sumup)


    n_spin=ssym%n_spin

    n_spin_dens = 1
    if (model_density) n_spin_dens = n_spin
    use_spin_dens = n_spin_dens > 1

    if( integralpar_2cob_field_efp) then
       if(size(symadapt_int_2cob_field(1)%int,3) ==0) goto 999 ! clean up and exit
    end if

    irrs: do i_ir=1,ssym%n_irrep
       ! quad_pmat and quad_wmat were allocated and computed by
       ! integral_calc_quad_densmat() in integral_calc_quad_module:
       pmunu => quad_pmat(i_ir,1)%m(:,:)

       ! building of density matrices finished

       if(integralpar_2cob_pc .or. integralpar_2cob_X) then
          call trace_energy( pmunu(:,:), &
               symadapt_int_2cob_nuc(i_ir), &
               qm_efp_energy)
       end if

       if(integralpar_2cob_pc_grad) then
          dim=totsym_grad_pc_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_pc_grad(1:dim,i_ir), &
               gradient_pc_totalsym(:))
       end if

       if(integralpar_2cob_X_grad) then
          dim=totsym_grad_dip_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_pd_grad(1:dim,i_ir), &
               gradient_dip_totalsym(:))
          call trace( pmunu(:,:), &
               symadapt_int_2cob_pd_torq(1:dim,i_ir), &
               torque_dip_totalsym(:))

          dim=totsym_grad_quad_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_pq_grad(1:dim,i_ir), &
               gradient_quad_totalsym(:))
          call trace( pmunu(:,:), &
               symadapt_int_2cob_pq_torq(1:dim,i_ir), &
               torque_quad_totalsym(:))

          dim=totsym_grad_oct_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_po_grad(1:dim,i_ir), &
               gradient_oct_totalsym(:))
          call trace( pmunu(:,:), &
               symadapt_int_2cob_po_torq(1:dim,i_ir), &
               torque_oct_totalsym(:))

          dim=totsym_grad_rep_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_rc_grad(1:dim,i_ir), &
               gradient_rep_totalsym(:))
       end if

       if(integralpar_2cob_ipd_grad) then
          dim=totsym_grad_idip_length
          call trace( pmunu(:,:), &
               symadapt_int_2cob_ipd_grad(1:dim,i_ir), &
               grad_idip_totalsym(:))
          call trace( pmunu(:,:), &
               symadapt_int_2cob_ipd_torq(1:dim,i_ir), &
               torque_idip_totalsym(:))
       end if

       if(integralpar_2cob_field_efp) then
          do i=1,totsym_field_length
            totalsym_field(i) = totalsym_field(i)                                &
                              + tr( pmunu(:,:)                                   &
                                  , symadapt_int_2cob_field(i_ir)%int(:,:,i,:,:) &
                                  )
          enddo
       end if
    enddo irrs

999 CONTINUE ! clean up and exit

    FPP_TIMER_STOP(sumup)

  end subroutine sum_up_grads_on_efp
#endif

  subroutine MO(M,N,p,h,f,t)
    !
    ! Purpose: Transforms prim. ints ``p'' in legacy storage
    !          to MO-basis and increments ``h''
    !
    ! case ``n''
    !   h(a,i) = h(a,i) + f * SUM(mu,nu) M(mu,a) * p(mu,nu) * N(nu,i)
    !
    ! case ``t''
    !   h(i,a) = h(i,a) +  ... the same ...
    !
    ! case ``s''
    !   h(a,i) = h(a,i) +  ... the same ... + (a<->i)
    !
    real(r8_kind), intent(in)    :: M(:,:)     ! (NMU,NA)
    real(r8_kind), intent(in)    :: N(:,:)     ! (NNU,NI)
    real(r8_kind), intent(in)    :: p(:,:,:,:) ! (n_exp2,n_exp1,n_if2,n_if1)
    real(r8_kind), intent(inout) :: h(:,:)     ! (NA,NI)
    real(r8_kind), intent(in)    :: f
    character(len=1), intent(in) :: t          ! "n", "t", or "s"
    ! NOTE: n2==n_exp2*n_if2, n1==n_exp1*n_if1
    !------------ Local Variables ------------------------------
    integer(i4_kind) :: n_exp1, n_if1, n_exp2, n_if2
    integer(i4_kind) :: i_exp1, i_if1, i_exp2, i_if2
    integer(i4_kind) :: NNU, NMU
    integer(i4_kind) :: nu, mu
    integer(i4_kind) :: NI, NA
    integer(i4_kind) :: i, a
#if _SLOW_REFERENCE_CODE
    real(r8_kind)    :: hai
#else
    real(r8_kind)    :: mna(size(N,1),size(M,2)) ! ( NNU x NA )
    real(r8_kind)    :: mia(size(N,2),size(M,2)) ! ( NI  x NA )
#endif
    !------------ Executable code ------------------------------

    FPP_TIMER_START(tomo)

    n_exp2 = size(p,1)
    n_exp1 = size(p,2)
    n_if2  = size(p,3)
    n_if1  = size(p,4)

    NA  = size(M,2)
    NMU = size(M,1)
    NI  = size(N,2)
    NNU = size(N,1)

    ASSERT(NMU==n_exp2*n_if2)
    ASSERT(NNU==n_exp1*n_if1)

    select case ( t )
    case ( "n", "N" )
      ASSERT(size(h,1)==NA)
      ASSERT(size(h,2)==NI)
    case ( "t", "T" )
      ASSERT(size(h,1)==NI)
      ASSERT(size(h,2)==NA)
    case ( "s", "S" )
      ASSERT(NA==NI)
      ASSERT(size(h,1)==size(h,2))
      ASSERT(size(h,1)==NA)
      ASSERT(size(h,2)==NI)
    case default
      ABORT('no such case')
    end select

#if _SLOW_REFERENCE_CODE
    ! FIXME: slow reference implementation:
      do i=1,NI
      do a=1,NA
        hai = 0.0_r8_kind
        nu = 0
        do i_if1=1,n_if1
          do i_exp1=1,n_exp1
            nu = nu + 1
            mu = 0
            do i_if2=1,n_if2
              do i_exp2=1,n_exp2
                mu = mu + 1
                hai = hai + M(mu,a)                      &
                          * p(i_exp2,i_exp1,i_if2,i_if1) &
                          * N(nu,i)
              enddo
            enddo
          enddo
        enddo
        select case ( t )
        case ( "n", "N" )
          h(a,i) = h(a,i) + hai * f
        case ( "t", "T" )
          h(i,a) = h(i,a) + hai * f
        case ( "s", "S" )
          h(a,i) = h(a,i) + hai * f
          h(i,a) = h(i,a) + hai * f
        end select
      enddo ! a MO loop
      enddo ! i MO loop
#else
    ! FIXME: depending on NA <> NI sum first over the biggest
    !        of them.

    ! MM1: compute mna(nu,a) = SUM(mu) M(mu,a) * p(mu,nu):
    ! loop over a:
    do a=1,NA

      ! loop over nu:
      nu = 0
      do i_if1=1,n_if1
        do i_exp1=1,n_exp1
          nu = nu + 1

          ! compute SUM(mu) M(mu,a) * p(mu,nu):
          mna(nu,a) = 0.0_r8_kind

          ! loop over mu:
          mu = 0
          do i_if2=1,n_if2
            do i_exp2=1,n_exp2
              mu = mu + 1

              mna(nu,a) = mna(nu,a) + M(mu,a) * p(i_exp2,i_exp1,i_if2,i_if1)
            enddo
          enddo
        enddo
      enddo
    enddo

    ! MM2: compute mia(i,a) = SUM(mu,nu) M(mu,a) * p(mu,nu) * N(nu,i):
    !                                  ^^^ mna(nu,a) ^^^
    do a=1,NA
      do i=1,NI
        mia(i,a) = 0.0_r8_kind
        do nu=1,NNU
          mia(i,a) = mia(i,a) + N(nu,i) * mna(nu,a)
        enddo
        ! NOTE the overall factor f:
        mia(i,a) = mia(i,a) * f
      enddo
    enddo

    select case ( t )
    case ( "n", "N" )
      h = h + transpose(mia)
    case ( "t", "T" )
      h = h + mia
    case ( "s", "S" )
      h = h + mia + transpose(mia)
    end select
#endif
    FPP_TIMER_STOP(tomo)
  end subroutine MO

  function tr(dmat,int) result(res)
    !
    ! Purpose: Retruns trace(P*I) = Sum(i,j) P_ij I_ij
    !
    real(r8_kind), intent(in) :: dmat(:,:)    ! (n_exp2 * n_if2, n_exp1 * n_if1)
    real(r8_kind), intent(in) :: int(:,:,:,:) ! (n_exp2, n_exp1, n_if2, n_if1)
    real(r8_kind)             :: res
    !------------ Local Variables ------------------------------
    integer(i4_kind) :: n_exp1, n_if1, n_exp2, n_if2
    integer(i4_kind) :: i_exp1, i_if1, i_exp2, i_if2
    integer(i4_kind) :: nu, mu
    !------------ Executable code ------------------------------
    n_exp2 = size(int,1)
    n_exp1 = size(int,2)
    n_if2  = size(int,3)
    n_if1  = size(int,4)

    res = 0.0_r8_kind

    nu = 0
    do i_if1=1,n_if1
       do i_exp1=1,n_exp1
          nu = nu + 1
          mu = 0
          do i_if2=1,n_if2
             do i_exp2=1,n_exp2
                mu = mu + 1
                res  = res  + dmat(mu,nu) * int(i_exp2,i_exp1,i_if2,i_if1)
             end do
          end do
       end do
    end do
  end function tr

  function tr3(dmat,int) result(res)
    !
    ! Purpose: Retruns trace(P*I) = Sum(i,j) P_ij I_ij
    !
    real(r8_kind), intent(in) :: dmat(:,:)      ! (n_exp2*n_if2,n_exp1*n_if1)
    real(r8_kind), intent(in) :: int(:,:,:,:,:) ! (n_exp2,n_exp1,NFF,n_if2,n_if1)
    real(r8_kind)             :: res(size(int,3)) !             (NFF)
    !------------ Local Variables ------------------------------
    integer(i4_kind) :: n_exp1, n_if1, n_exp2, n_if2
    integer(i4_kind) :: i_exp1, i_if1, i_exp2, i_if2
    integer(i4_kind) :: nu, mu
    !------------ Executable code ------------------------------
    n_exp2 = size(int,1)
    n_exp1 = size(int,2)
    !NFF   = size(int,3)
    n_if2  = size(int,4)
    n_if1  = size(int,5)

    res = 0.0_r8_kind

    nu = 0
    do i_if1=1,n_if1
       do i_exp1=1,n_exp1
          nu = nu + 1
          mu = 0
          do i_if2=1,n_if2
             do i_exp2=1,n_exp2
                mu = mu + 1
                res  = res  + dmat(mu,nu) * int(i_exp2,i_exp1,:,i_if2,i_if1)
             end do
          end do
       end do
    end do
  end function tr3

  subroutine trace(densmat,integral,gradient)
    !
    ! Purpose: Performs F(k) = F(k) + Sum(i,j) P_ij I_ij(k)
    !
    ! Input
    ! I_ij(k) : integral(n_grad)%int(n_exp2,n_exp1,n_if2,n_if1)
    ! P_ij    : densmat (n_exp2*n_if2,n_exp1*n_if1)
    !
    ! Output
    ! F(k)    : gradient(n_grad)
    !
    real(kind=r8_kind)               , intent(in   ) :: densmat(:,:)
    type(symadapt_totsym_2c_int_type), intent(in   ) :: integral(:)
    real(kind=r8_kind)               , intent(inout) :: gradient(:)
    !------------ Local Variables ------------------------------
    integer(kind=i4_kind) :: i_grad
    integer(kind=i4_kind) :: n_grad
    !------------ Executable code ------------------------------

    n_grad = size(gradient)
    ASSERT(n_grad==size(integral))
    do i_grad=1,n_grad
       gradient(i_grad) = gradient(i_grad) + tr(densmat,integral(i_grad)%int)
    end do
  end subroutine trace

#ifdef WITH_EFP
  subroutine trace_energy(densmat,integral,energy)
    !
    ! Purpose: Performs E = E + Sum(i,j) P_ij I_ij
    !
    ! Input
    ! I_ij : integral%int(n_exp2,n_exp1,n_if2,n_if1)
    ! P_ij : densmat (n_exp2*n_if2,n_exp1*n_if1)
    !
    ! Output
    ! E    : energy
    !
    real(kind=r8_kind)               , intent(in   ) :: densmat(:,:)
    type(symadapt_totsym_2c_int_type), intent(in   ) :: integral
    real(kind=r8_kind)               , intent(inout) :: energy
    !------------ Local Variables ------------------------------
    !------------ Executable code ------------------------------

    energy = energy + tr(densmat,integral%int)

  end subroutine trace_energy
#endif

  subroutine trace_dervs(densmat,integral,dervs)
    !
    ! Purpose: Performs F(k) = F(k) + Sum(i,j) P_ij I_ij(k)
    !
    ! Input
    ! I_ij(k) : integral(n_grad)%int(n_exp2,n_exp1,n_if2,n_if1)
    ! P_ij    : densmat (n_exp2*n_if2,n_exp1*n_if1)
    !
    ! Output
    ! F(k)    : dervs(n_grad,n_grad)
    !
    real(kind=r8_kind)               , intent(in   ) :: densmat(:,:)
    type(symadapt_totsym_2c_int_type), intent(in   ) :: integral(:,:)
    real(kind=r8_kind)               , intent(inout) :: dervs(:,:)
    !------------ Local Variables ------------------------------
    integer(kind=i4_kind) :: i_grad, k2dr
    !------------ Executable code ------------------------------

   do k2dr=1,size(dervs,2)
    do i_grad=1,size(dervs,1)
       dervs(i_grad,k2dr) = dervs(i_grad,k2dr) + tr(densmat,integral(i_grad,k2dr)%int)
    end do
   enddo
  end subroutine trace_dervs

  subroutine cpks_trace(densmat,integral,gradient)
    !
    ! Purpose: Performs F(k) = F(k) + Sum(i,j) P_ij I_ij(k)
    !
    ! Input
    ! I_ij(k) : integral(n_grad)%int(n_exp2,n_exp1,n_if2,n_if1)
    ! P_ij    : densmat (n_exp2*n_if2,n_exp1*n_if1)
    !
    ! Output
    ! F(k)    : gradient(n_grad)
    !
    real(kind=r8_kind)               , intent(in   ) :: densmat(:,:)
    type(symadapt_totsym_3c_int_type), intent(in   ) :: integral(:)
    real(kind=r8_kind)               , intent(inout) :: gradient(:,:)
    !------------ Local Variables ------------------------------
    integer(kind=i4_kind) :: i_grad
    !------------ Executable code ------------------------------

    do i_grad=1,size(gradient,2)
       gradient(:,i_grad) = gradient(:,i_grad) + tr3(densmat,integral(i_grad)%int)
    end do
  end subroutine cpks_trace

  subroutine symadp_pseudo2D()
    ! Purpose: In case of pseudo 2D irreps the symmetry adaption for
    !          for total symmetric operator ( Goerling, Thesis, 2.73)
    !          doesnt apply. An additional symmitrization is necessary:
    !
    !          <ax|O|bx>:=0.5*(<ax|O|bx>+<ay|O|by>)
    !          <ay|O|by>=<ax|O|bx>
    !          <ax|O|by>:=0.5*(<ax|O|by>-<ay|O|bx>)
    !          <ay|O|bx>=-<ax|O|by>
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: i_ir
    integer(kind=i4_kind) :: n_independent1_1dim, n_independent2_1dim
    ! multiplicities of the original irreps, only used for pseudo 2D
    integer(i4_kind) :: i_grad,k2dr

    ! in the case of groups with 2D pseudo irreps we have to make
    ! an adittional symmetry adaption
    do i_ir=1,symmetry_data_n_irreps()
       if( .not. symmetry_data_pseudo(i_ir) ) CYCLE

        n_independent1_1dim=ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts/2
        n_independent2_1dim=ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts/2
        if ( integralpar_2cob_kin ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,&
                saint_p2=symadapt_int_2cob_kin(i_ir)%int)
        endif
        if ( integralpar_2cob_nuc ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,&
                saint_p2=symadapt_int_2cob_nuc(i_ir)%int)
           if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) &
                call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,&
                saint_p2=symadapt_int_2cob_nuc_pseudo(i_ir)%int)
        endif
        if ( integralpar_2cob_pvsp ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                symadapt_int_2cob_pvsp(i_ir)%int)
        endif
        if ( integralpar_2cob_ol ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,&
                saint_p2=symadapt_int_2cob_ol(i_ir)%int)
        endif
        if ( integralpar_3c_xc ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,saint_p3=&
                symadapt_int_3c_xc(i_ir)%int)
        endif
        if ( integralpar_3c_co ) then
           call pseudo_2d_symadapt&
                (n_independent1_1dim,n_independent2_1dim,saint_p3=&
                symadapt_int_3c_co(i_ir)%int)
        endif
        if ( integralpar_2cob_potential ) then
           call pseudo_2d_symadapt&                           !!!!!!!!!!!!!
                (n_independent1_1dim,n_independent2_1dim,&    !!!!!!!!!!!!!
                saint_p3=symadapt_int_2cob_poten(i_ir)%int)   !!!!!!!!!!!!!
        endif
        if ( integralpar_2cob_field ) then
           call pseudo_2d_symadapt&                           !!!!!!!!!!!!!
                (n_independent1_1dim,n_independent2_1dim,&    !!!!!!!!!!!!!
                saint_p3=symadapt_int_2cob_field(i_ir)%int)   !!!!!!!!!!!!!
        endif

        if ( integralpar_solv_grad .and. integralpar_cpksdervs  ) then
           do i_grad=1,gradient_data_n_spin_gradients                    !!!!!!!!!!!!!!
              call pseudo_2d_symadapt&                                   !!!!!!!!!!!!!!
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&   !!!!!!!!!!!!!!
                   symadapt_int_3cob_solv_grad(i_grad,i_ir)%int)         !!!!!!!!!!!!!!
           end do!!!!!!!!!!!!!!
           if(with_pc .and. .not.fixed_pc) then
              do i_grad=1,totsym_field_length                    !!!!!!!!!!!!!!
                 call pseudo_2d_symadapt&                                   !!!!!!!!!!!!!!
                      (n_independent1_1dim,n_independent2_1dim,saint_p2=&   !!!!!!!!!!!!!!
                      symadapt_int_3cob_solv_grad_pc(i_grad,i_ir)%int)         !!!!!!!!!!!!!!
              end do!!!!!!!!!!!!!!
           end if
        endif

        int3cob_grad: if ( integralpar_3cob_grad ) then
#ifdef WITH_EPE
          if( ewpc_n.ne.0.and.calc_cluster_epe_energy) call pseudo_2d_symadapt&
               (n_independent1_1dim,n_independent2_1dim,saint_p2=&
               symadapt_int_3cob_epe(1,i_ir)%int)
#endif
           do i_grad=1,gradient_data_n_spin_gradients
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_3cob_grad(i_grad,i_ir)%int)
           enddo
#ifndef   no_cpks_coul_grads
        if(integralpar_cpksdervs) then
           do i_grad=1,size(symadapt_int_cpks_coul_grad,1)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p3=&
                   symadapt_int_cpks_coul_grad(i_grad,i_ir)%int)
           enddo
        endif
#endif

         if(.not.integralpar_cpks_contribs.and.integralpar_dervs) then

           do i_grad=1,size(symadapt_int_coul_dervs,1)
            do k2dr=1,size(symadapt_int_coul_dervs,2)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int)
            enddo
           enddo
           do i_grad=1,size(symadapt_int_2cob_ol_dervs,1)
            do k2dr=1,size(symadapt_int_2cob_ol_dervs,2)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int)
            enddo
           enddo
         endif

        if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
           do i_grad=1,gradient_data_n_spin_gradients
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_pseudo_grad(i_grad,i_ir)%int)
           end do
         endif

        split_gr:   if (split_gradients) then
              do i_grad=1,size(symadapt_int_2cob_ks_grad,1)
                 call pseudo_2d_symadapt( &
                      n_independent1_1dim,n_independent2_1dim,saint_p2= &
                      symadapt_int_2cob_ks_grad(i_grad,i_ir)%int)
              end do
              do i_grad=1,gradient_data_n_gradients
                 call pseudo_2d_symadapt( &
                      n_independent1_1dim,n_independent2_1dim,saint_p2= &
                      symadapt_int_3cob_nuc_grad(i_grad,i_ir)%int)
              end do
              if (model_density) then
                 do i_grad=1,gradient_data_n_gradients
                    call pseudo_2d_symadapt( &
                         n_independent1_1dim,n_independent2_1dim,saint_p2= &
                         symadapt_int_3cob_coul_grad(i_grad,i_ir)%int)
                 end do
              endif
           endif split_gr

        endif int3cob_grad

        if ( integralpar_2cob_nuc_grad ) then
           do i_grad=1,gradient_data_n_gradients
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_nuc_grad(i_grad,i_ir)%int)
           end do
        endif
        if ( integralpar_2cob_pvsp_grad ) then
           do i_grad=1,gradient_data_n_gradients
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_pvsp_grad(i_grad,i_ir)%int)
           end do
        endif
        if ( integralpar_2cob_kin_grad ) then
           do i_grad=1,size(symadapt_int_2cob_kin_grad,1)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_kin_grad(i_grad,i_ir)%int)
           end do
        endif

        ol_grad: if ( integralpar_2cob_ol_grad ) then
           do i_grad=1,size(symadapt_int_2cob_ol_grad,1)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_ol_grad(i_grad,i_ir)%int)
           end do

           if(integralpar_relativistic) then
              do i_grad=1,size(symadapt_int_2cob_ol_grad,1)
                 call pseudo_2d_symadapt&
                      (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                      symadapt_int_2cob_ol_rel_grad(i_grad,i_ir)%int)
#if 0 /* to be treated not only for relativistic if calculated */
            if(integralpar_dervs) then
            do k2dr=1,size(symadapt_int_2cob_ol_grad,1)
              call pseudo_2d_symadapt&
                   (n_independent1_1dim,n_independent2_1dim,saint_p2=&
                   symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int)     !relativistic ???
            enddo
            WARN('and pseudo_2d of coul_dervs?')
            endif
#endif
              end do
           end if

        endif ol_grad
    end do! loop over irreps
  end subroutine symadp_pseudo2D

  subroutine pseudo_2d_symadapt(n_independent1_1dim,n_independent2_1dim,saint_p2,&
       saint_p3)
    ! Purpose: In case of pseudo 2D irreps the symmetry adaption for
    !          for total symmetric operator ( Goerling, Thesis, 2.73)
    !          doesnt apply. An additional symmitrization is necessary:
    !
    !          <ax|O|bx>:=0.5*(<ax|O|bx>+<ay|O|by>)
    !          <ay|O|by>=<ax|O|bx>
    !          <ax|O|by>:=0.5*(<ax|O|by>-<ay|O|bx>)
    !          <ay|O|bx>=-<ax|O|by>
    integer(kind=i4_kind) :: n_independent1_1dim,n_independent2_1dim
    ! number of independent funcions in the original 1-dim irreps =
    ! half of the number of independent functions of the 2d pseudo irrep
    real(kind=r8_kind), optional, pointer :: saint_p2(:,:,:,:),&
         saint_p3(:,:,:,:,:)
    ! symmetry adapted integrals
    !------------ Declaration of local variables ------------------
    integer(kind=i4_kind) :: i_if1,i_if2
    if(present(saint_p2)) then
       do i_if1=1,n_independent1_1dim
          do i_if2=1,n_independent2_1dim
             saint_p2(:,:,i_if2,i_if1) = 0.5_r8_kind* &
                  ( saint_p2(:,:,i_if2,i_if1)+ &
                  saint_p2(:,:,i_if2+n_independent2_1dim, i_if1+n_independent1_1dim))
             saint_p2(:,:,i_if2+n_independent2_1dim, i_if1+n_independent1_1dim) = &
                  saint_p2(:,:,i_if2,i_if1)
             saint_p2(:,:,i_if2,i_if1+n_independent1_1dim) = &
                  0.5_r8_kind*(saint_p2(:,:,i_if2,i_if1+n_independent1_1dim)-&
                  saint_p2(:,:,i_if2+n_independent2_1dim,i_if1))
             saint_p2(:,:,i_if2+n_independent2_1dim,i_if1) = - &
                  saint_p2(:,:,i_if2,i_if1+n_independent1_1dim)
          end do
       end do
    end if

    if(present(saint_p3)) then
       do i_if1=1,n_independent1_1dim
          do i_if2=1,n_independent2_1dim
             saint_p3(:,:,:,i_if2,i_if1) = 0.5_r8_kind* &
                  ( saint_p3(:,:,:,i_if2,i_if1)+ &
                  saint_p3(:,:,:,i_if2+n_independent2_1dim, i_if1+n_independent1_1dim))
             saint_p3(:,:,:,i_if2+n_independent2_1dim, i_if1+n_independent1_1dim) = &
                  saint_p3(:,:,:,i_if2,i_if1)
             saint_p3(:,:,:,i_if2,i_if1+n_independent1_1dim) =  &
                  0.5_r8_kind*(saint_p3(:,:,:,i_if2,i_if1+n_independent1_1dim)-&
                  saint_p3(:,:,:,i_if2+n_independent2_1dim,i_if1))
             saint_p3(:,:,:,i_if2+n_independent2_1dim,i_if1) = - &
                  saint_p3(:,:,:,i_if2,i_if1+n_independent1_1dim)
          end do
       end do
    end if
  end subroutine pseudo_2d_symadapt
  !**************************************************************

end subroutine integral_calc_quad_2cob3c

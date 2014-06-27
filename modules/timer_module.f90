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
module  timer_module
!---------------------------------------------------------------
!
!  Purpose: Timing. Contains  all timers used in program.   Thus it is
!  USEd at quite a few places. Any import of another module here has a
!  huge   potential   for    cyclic   dpendencies.    Dont   do   that
!  (unnecessarily).
!
!  Module called by: main_scf, main_master, ...
!
!  see time_module for decription of timing concept!
!
!  Author: TB
!  Date: 2/96
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: TB
! Date:   11/96
! Description: Complete restructuring:
!              + print routines for seperate program parts
!              + output parameters to control amount of output
!              + lots of new integral timings
!              + routines to collect slave integral timings
!              + integral timers are now arrays to allow
!                multiple executions of integral part
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification (Please copy before editing)
! Author: ...
! Date:   9/99
! Description: +EPE-timing
!
!----------------------------------------------------------------

# include "def.h"
  use type_module, only: r8_kind, i4_kind ! type specification parameters
  use time_module, only: timer_type
  ! FIXME: by commenting this a cyclic dependency is broken:
!!$  use integralpar_module, only: n_int_parts => integralpar_n_int_parts
!:UB[ modified
use options_module, only: options_relativistic, options_xcmode, &
                          xcmode_model_density, xcmode_extended_mda
!:UB]
implicit none
save ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================


!------------ Declaration of constants and variables ------------
public :: timer_type

! Total number  of integral parts (for timing).   FIXME: when changing
! this adpat  a copy of  this number in integralpar_module.   There is
! some  cyclic  dependency  if  one  imports  integralpar_module  into
! timer_module:
integer, parameter :: n_int_parts = 4

type(timer_type), public, dimension(n_int_parts) :: &
     timer_int, &
     timer_int_2cff, &
     timer_int_interrupt_2cff, &
     timer_int_calcsum_2cff, &
     timer_int_primcont_2cff, &
     timer_int_prim_2cff, &
     timer_int_cont_2cff, &
     timer_int_idle_2cff, &
     timer_int_send_2cff, &
     timer_int_sendnorm_2cff, &
     timer_int_waitnorm_2cff, &
     timer_int_communication_2cff, &
     timer_int_quadrupel_2cff, &
     timer_int_receive_2cff, &
     timer_int_norm_2cff, &
     timer_int_rewrite_2cff, &
     timer_int_newtask_2cff, &
     timer_int_2cob3c, &
     timer_int_interrupt_2cob3c, &
     timer_int_calc_2cob3c, &
     timer_int_calcsum_2cob3c, &
     timer_int_prim_2cob3c, &
     timer_int_cont_2cob3c, &
     timer_int_symadapt_2cob3c, &
     timer_int_symadapt_pvxp_2cob3c, &
     timer_int_symadapt_npxp_2cob3c, &
     timer_int_idle_2cob3c, &
     timer_int_write_2cob3c, &
     timer_int_send_2cob3c, &
     timer_int_commsend_2cob3c, &
     timer_int_pack_2cob3c, &
     timer_int_quadrupel_2cob3c, &
     timer_int_receive_2cob3c, &
     timer_int_norm_2cob3c, &
     timer_int_rewrite_2cob3c, &
     timer_int_newtask_2cob3c, &
     timer_int_rel
type(timer_type), public :: &
     timer_efield, &
     timer_scf, &
     timer_scf_ham, &
     timer_scf_eigen, &
     timer_scf_reoc, &
     timer_scf_chfit, &
     timer_scf_cycle, &
     timer_scf_preparations, &
     timer_scf_xc, &
     timer_initialisation, &
     timer_first_delay, &
     timer_1, &
     timer_2, &
     timer_3, &
     timer_4, &
     timer_5, &
     timer_symm, &
     timer_grid_setup, &
     timer_grid_setup_ph, &
     timer_grid_orbitals_sum, &
     timer_grid_density_sum, &
     timer_grid_functionals_sum, &
     timer_grid_xcbuild_sum, &
!:UB[ added
     timer_grid_trunc_dens_sum, &
     timer_grid_trunc_funcs_sum, &
     timer_grid_num_metric_sum, &
     timer_grid_projections_sum, &
     timer_grid_fit_coeffs_sum
type(timer_type), public :: &
!:UB]
     timer_grid_orbitals, &
     timer_grid_density, &
     timer_grid_functionals, &
     timer_grid_xcbuild, &
     timer_int_fit_density_calc, &     !AG
     timer_int_fit_density_calc_sum, & !AG
!:UB[ added
     timer_grid_trunc_dens, &
     timer_grid_trunc_funcs, &
     timer_grid_num_metric, &
     timer_grid_projections, &
     timer_grid_fit_coeffs, &
!:UB]
     timer_gridph_orbitals, &
     timer_gridph_density, &
     timer_gridph_functionals, &
     timer_gridph_grdwts, &
     timer_gridph_gradients, &
!:UB[ added
     timer_gridph_trunc_dens, &
     timer_gridph_trunc_funcs, &
     timer_gridph_projections, &
     timer_gridph_fit_coeffs, &
!:UB]
     timer_gridph, &
     timer_dipole, &
     timer_dipole_integral, &
     timer_dipole_calculate, &
     timer_dipole_print, &
     timer_rel_total, &
     timer_rel_diag, &
     timer_rel_trafo, &
     timer_rel_gradients
type(timer_type), public :: &
     timer_resp,&
     timer_resp_preparations,&
     timer_resp_setup,&
     timer_resp_header,&
     timer_resp_dipole,&
     timer_resp_2index,&
     timer_resp_2index_integration,&
     timer_resp_3index,&
     timer_resp_3index_integration,&
     timer_resp_Coulomb,&
     timer_resp_Coulomb_reading,&
     timer_resp_Coulomb_writing,&
     timer_resp_Coulomb_trafo,&
     timer_resp_allClb,&
     timer_resp_calcFG,&
     timer_resp_addCoul,&
     timer_resp_all3ind,&
     timer_resp_rhsh

type(timer_type), public :: &
     timer_epe,                       &
     timer_epe_cycle,                 &
     timer_epe_forces,                &
     timer_get_epe_energy,            &
     timer_collect_grads,             &
     timer_send_epe_ref,              &
     timer_main_epe,                  &
     timer_latt_epe_forces,           &
     timer_defect_contributions,      &
     timer_collect_def_contributions, &
     timer_defect_ext

type(timer_type), public :: timer_ll_integral, &
                            timer_ls_integral, &
                            timer_ss_integral, &
                            timer_new_ll_integral, &
                            timer_new_ls_integral, &
                            timer_new_ss_integral, &
                            timer_ang3c, &
                            timer_rad3c, &
                            timer_prod_nested

!------------ public functions and subroutines ------------------
public timer_setup, timer_print_summary, timer_print_scfcycle, &
     timer_grid_small_to_large, timer_print_slavetiming, &
     timer_print_integral, timer_print_scf, timer_print_postscf, &
     timer_print_response

public :: timer_gather_slave_int_timing!()


!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of private constants and variables ----
type(timer_type), private, dimension(n_int_parts) :: &
     t_sl_int, &
     t_sl_int_2cff, &
     t_sl_int_interrupt_2cff, &
     t_sl_int_calcsum_2cff, &
     t_sl_int_idle_2cff, &
     t_sl_int_norm_2cff, &
     t_sl_int_2cob3c, &
     t_sl_int_interrupt_2cob3c, &
     t_sl_int_calcsum_2cob3c, &
     t_sl_int_idle_2cob3c, &
     t_sl_int_receive_2cob3c, &
     t_sl_int_norm_2cob3c, &
     t_sl_int_rewrite_2cob3c, &
     t_sl_int_rel


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

#define init_timer(t) init_timer_v2(t, __STRING(t))


   !*************************************************************
   subroutine timer_setup()
   !  Purpose: initialize all timers in module
   use time_module, only: timer_type, INIT_TIMER, init_timer_v2
   implicit none
   !** End of interface *****************************************
   integer :: i_int

   do i_int = 1, size(timer_int) ! n_int_parts
      call init_timer(timer_int(i_int))
      call init_timer(timer_int_2cff(i_int))
      call init_timer(timer_int_interrupt_2cff(i_int))
      call init_timer(timer_int_calcsum_2cff(i_int))
      call init_timer(timer_int_prim_2cff(i_int))
      call init_timer(timer_int_primcont_2cff(i_int))
      call init_timer(timer_int_cont_2cff(i_int))
      call init_timer(timer_int_idle_2cff(i_int))
      call init_timer(timer_int_send_2cff(i_int))
      call init_timer(timer_int_sendnorm_2cff(i_int))
      call init_timer(timer_int_waitnorm_2cff(i_int))
      call init_timer(timer_int_communication_2cff(i_int))
      call init_timer(timer_int_quadrupel_2cff(i_int))
      call init_timer(timer_int_receive_2cff(i_int))
      call init_timer(timer_int_norm_2cff(i_int))
      call init_timer(timer_int_rewrite_2cff(i_int))
      call init_timer(timer_int_newtask_2cff(i_int))
      call init_timer(timer_int_2cob3c(i_int))
      call init_timer(timer_int_interrupt_2cob3c(i_int))
      call init_timer(timer_int_calc_2cob3c(i_int))
      call init_timer(timer_int_calcsum_2cob3c(i_int))
      call init_timer(timer_int_prim_2cob3c(i_int))
      call init_timer(timer_int_cont_2cob3c(i_int))
      call init_timer(timer_int_symadapt_2cob3c(i_int))
      call init_timer(timer_int_symadapt_pvxp_2cob3c(i_int))
      call init_timer(timer_int_symadapt_npxp_2cob3c(i_int))
      call init_timer(timer_int_idle_2cob3c(i_int))
      call init_timer(timer_int_write_2cob3c(i_int))
      call init_timer(timer_int_send_2cob3c(i_int))
      call init_timer(timer_int_pack_2cob3c(i_int))
      call init_timer(timer_int_commsend_2cob3c(i_int))
      call init_timer(timer_int_quadrupel_2cob3c(i_int))
      call init_timer(timer_int_receive_2cob3c(i_int))
      call init_timer(timer_int_norm_2cob3c(i_int))
      call init_timer(timer_int_rewrite_2cob3c(i_int))
      call init_timer(timer_int_newtask_2cob3c(i_int))
      call init_timer(timer_int_rel(i_int))
   enddo
   call init_timer(timer_efield)
   call init_timer(timer_scf)
   call init_timer(timer_scf_ham)
   call init_timer(timer_scf_eigen)
   call init_timer(timer_scf_reoc)
   call init_timer(timer_scf_chfit)
   call init_timer(timer_scf_cycle)
   call init_timer(timer_scf_preparations)
   call init_timer(timer_grid_setup)
   call init_timer(timer_scf_xc)
   call init_timer(timer_initialisation)
   call init_timer(timer_first_delay)
   call init_timer(timer_1)
   call init_timer(timer_2)
   call init_timer(timer_3)
   call init_timer(timer_4)
   call init_timer(timer_5)
   call init_timer(timer_symm)
   call init_timer(timer_grid_orbitals_sum)
   call init_timer(timer_grid_density_sum)
   call init_timer(timer_grid_functionals_sum)
!:UB[ added
   call init_timer(timer_grid_xcbuild_sum)
   call init_timer(timer_grid_trunc_dens_sum)
   call init_timer(timer_grid_trunc_funcs_sum)
   call init_timer(timer_grid_num_metric_sum)
   call init_timer(timer_grid_projections_sum)
   call init_timer(timer_grid_fit_coeffs_sum)
!:UB]
!AG[---------------------------------------------
   call init_timer(timer_int_fit_density_calc)
   call init_timer(timer_int_fit_density_calc_sum)
!AG]---------------------------------------------
   call init_timer(timer_gridph)
   call init_timer(timer_rel_total)
   call init_timer(timer_rel_diag)
   call init_timer(timer_rel_trafo)
   call init_timer(timer_rel_gradients)
   call init_timer(timer_dipole)
   call init_timer(timer_dipole_integral)
   call init_timer(timer_dipole_calculate)
   call init_timer(timer_dipole_print)
   call init_timer(timer_resp)
   call init_timer(timer_resp_preparations)
   call init_timer(timer_resp_setup)
   call init_timer(timer_resp_header)
   call init_timer(timer_resp_dipole)
   call init_timer(timer_resp_2index)
   call init_timer(timer_resp_2index_integration)
   call init_timer(timer_resp_3index)
   call init_timer(timer_resp_3index_integration)
   call init_timer(timer_resp_Coulomb)
   call init_timer(timer_resp_Coulomb_reading)
   call init_timer(timer_resp_Coulomb_writing)
   call init_timer(timer_resp_Coulomb_trafo)
   call init_timer(timer_resp_allClb)
   call init_timer(timer_resp_calcFG)
   call init_timer(timer_resp_addCoul)
   call init_timer(timer_resp_all3ind)
   call init_timer(timer_resp_rhsh)
!AG[
   call init_timer(timer_epe)
   call init_timer(timer_epe_cycle)
   call init_timer(timer_epe_forces)
   call init_timer(timer_get_epe_energy)
   call init_timer(timer_collect_grads)
   call init_timer(timer_send_epe_ref)
   call init_timer(timer_main_epe)
   call init_timer(timer_latt_epe_forces)
   call init_timer(timer_defect_contributions)
   call init_timer(timer_collect_def_contributions)
   call init_timer(timer_defect_ext)
!AG]
   call init_timer(timer_prod_nested)
   call init_timer(timer_ll_integral)
   call init_timer(timer_ls_integral)
   call init_timer(timer_ss_integral)
   call init_timer(timer_new_ll_integral)
   call init_timer(timer_new_ls_integral)
   call init_timer(timer_new_ss_integral)
   call init_timer(timer_ang3c)
   call init_timer(timer_rad3c)

   ! initialise slve timers as well on slaves
   ! in order to make timer_print_summary() to work on them
   do i_int = 1, size(t_sl_int) ! n_int_parts
      call init_timer( t_sl_int(i_int) )
      call init_timer( t_sl_int_2cff(i_int) )
      call init_timer( t_sl_int_interrupt_2cff(i_int) )
      call init_timer( t_sl_int_calcsum_2cff(i_int) )
      call init_timer( t_sl_int_idle_2cff(i_int) )
      call init_timer( t_sl_int_norm_2cff(i_int) )
      call init_timer( t_sl_int_2cob3c(i_int) )
      call init_timer( t_sl_int_interrupt_2cob3c(i_int) )
      call init_timer( t_sl_int_calcsum_2cob3c(i_int) )
      call init_timer( t_sl_int_idle_2cob3c(i_int) )
      call init_timer( t_sl_int_receive_2cob3c(i_int) )
      call init_timer( t_sl_int_norm_2cob3c(i_int) )
      call init_timer( t_sl_int_rewrite_2cob3c(i_int) )
      call init_timer( t_sl_int_rel(i_int) )
   enddo
   end subroutine timer_setup
   !*************************************************************


   !*************************************************************
   subroutine timer_print_summary (integralpar_int_part_name)
   !  Purpose: prints summary of timing of program
     use time_module, only: timer_type,print_timer,print_time
     use iounitadmin_module, only: output_unit
     use operations_module, only: operations_gradients
     use options_module, only: options_relativistic
     use output_module
     implicit none
     character(len=*), intent(in) :: integralpar_int_part_name(:)
     !** End of interface *****************************************

   integer :: i_int
!:UB[ added
   logical :: extended_mda, xcmda_called
!:UB]

   if (output_unit <= 0) return ! FIXME: maybe dont call at all?

!:UB[ added
   extended_mda = options_xcmode() == xcmode_extended_mda
   xcmda_called = options_xcmode() == xcmode_model_density .or. extended_mda
!:UB]
   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "#######################"
   write(output_unit,*) "## Summary of Timing ##"
   write(output_unit,*) "#######################"
   write(output_unit,*) ""
   write(output_unit,*) ""
   call print_timer(timer_initialisation,output_unit, &
        "Initialisation", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_first_delay,output_unit, &
        "first Delay", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_1,output_unit,  "1", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_2,output_unit,  "2", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_3,output_unit,  "3", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_4,output_unit,  "4", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_5,output_unit,  "5", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_symm,output_unit, &
        "Entire Symmetry part", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   do i_int = 1, size(timer_int) ! n_int_parts
      call print_timer(timer_int(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: Total", &
           nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
      call print_timer(timer_int_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Total", &
           nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
      call print_timer(timer_int_quadrupel_2cff(i_int),output_unit,integralpar_int_part_name(i_int)//  &
           "Integral Part: 2 center fitfunction integrals: Quadrupel", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_prim_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Primitive Integrals and Symmetry Adaption", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_cont_2cff(i_int),output_unit,integralpar_int_part_name(i_int)//  &
           "Integral Part: 2 center fitfunction integrals: Contractionn and Renormalisation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_send_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Send", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_calcsum_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_interrupt_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_idle_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_receive_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: receive", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_norm_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_rewrite_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: rewrite", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_newtask_2cff(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: newtask", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Total", &
           nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
      call print_timer(timer_int_quadrupel_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Quadrupel", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_prim_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Primitive Integrals", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_cont_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Contraction", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_symadapt_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_symadapt_pvxp_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption of PVXP", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_symadapt_npxp_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption except PVXP", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_send_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Send", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_commsend_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: comm Send", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_pack_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Pack", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_write_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Write", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_calcsum_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_interrupt_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_idle_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_receive_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: receive", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_norm_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_rewrite_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: rewrite", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(timer_int_newtask_2cob3c(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: newtask", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_rel(i_int),output_unit,integralpar_int_part_name(i_int)// &
           "Integral Part: relativistic transformations: Total", &
           nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   enddo
   call print_timer(timer_efield,output_unit, &
        "Calculation of integrals of external electrical field", &
        diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_scf,output_unit, &
        "Entire SCF part", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_scf_preparations,output_unit, &
        "Preparations in SCF part")
   if ( output_timing_detailedsummary ) &
   call print_timer(timer_grid_setup,output_unit, &
        "Grid setup")
   call print_timer(timer_scf_cycle,output_unit, &
        "Entire SCF cycle",diff=.false. &
        ,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_scf_xc,output_unit, &
        "Numeric Build up of Exchange Part of Hamiltonian", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_ham,output_unit, &
        "Build up of Hamiltonian", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_eigen,output_unit, &
        "Eigensolver", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_scf_reoc,output_unit, &
        "Reoccupation of Orbitals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_chfit,output_unit, &
        "Chargefit", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_grid_orbitals_sum,output_unit, &
!:UB[ modified
        "Calculation of orbitals", &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_grid_density_sum,output_unit, &
!:UB[ modified
        "Calculation of density",  &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_grid_trunc_dens_sum,output_unit, &
        "Truncation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_grid_functionals_sum,output_unit, &
!:UB[ modified
        "Calculation of xc-functionals",  &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
!AG[
   if( output_timing_detailedsummary .and. xcmda_called ) &
       call print_timer(timer_int_fit_density_calc_sum, output_unit, &
       "Fitted density calc. (int.)", &
       diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!AG]
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_grid_trunc_funcs_sum,output_unit, &
        "Truncation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_grid_projections_sum,output_unit, &
        "Calculation of numerical projections",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary .and. extended_mda ) &
        call print_timer(timer_grid_num_metric_sum,output_unit, &
        "Calculation of numerical metric",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_grid_fit_coeffs_sum,output_unit, &
        "Calculation of fitting coefficients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
!:UB[ modified
   if ( output_timing_detailedsummary .and. .not.xcmda_called ) &
!:UB]
        call print_timer(timer_grid_xcbuild_sum,output_unit, &
        "Construction of the numerical XC-Hamiltonian",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_gridph,output_unit, &
        "Post SCF: Complete grid_part", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_gridph_orbitals,output_unit, &
        "Post scf: calculation of orbitals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_gridph_density,output_unit, &
        "Post scf: calculation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_gridph_trunc_dens,output_unit, &
        "Post scf: truncation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_gridph_functionals,output_unit, &
        "Post scf: calculation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_gridph_trunc_funcs,output_unit, &
        "Post scf: truncation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_gridph_projections,output_unit, &
        "Post scf: calculation of projections",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary .and. xcmda_called ) &
        call print_timer(timer_gridph_fit_coeffs,output_unit, &
        "Post scf: calculation of fit coefficients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_gridph_grdwts,output_unit, &
        "Post scf: calculation of grdwts",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_gridph_gradients,output_unit, &
        "Post scf: calculation of gradients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_dipole,output_unit, &
        "Dipole: Complete part", &
        diff=.false.,absolut=.false.,sum=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_dipole_integral,output_unit, &
        "Dipole: Integrals", &
        diff=.false.,absolut=.false.,sum=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_dipole_calculate,output_unit, &
        "Dipole: Calculation", &
        diff=.false.,absolut=.false.,sum=.true.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_dipole_print,output_unit, &
        "Dipole: Output", &
        diff=.false.,absolut=.false.,sum=.true.)
   if(operations_gradients.and.options_relativistic) &
        call print_timer(timer_rel_gradients,output_unit, &
        "Calculation of rel. gradients", &
        nbr=.true.,diff=.false.,sum=.true.,average=.true.,absolut=.false.)
   call print_timer(timer_resp,output_unit, &
        "Entire RESPONSE part", diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_resp_preparations,output_unit, &
        "Response: Preparations", diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_setup,output_unit, &
        "Response: Setup", diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_header,output_unit, &
        "Response: Writing header files", diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_resp_dipole,output_unit, &
        "Response: Rewriting transition dipole moments",&
        & diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_resp_2index,output_unit, &
        "Response: Numerical 2-index-integrals", diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_2index_integration,output_unit, &
        &"Response: Numerical 2-index-integrals: integration",&
        & diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_resp_3index,output_unit, &
        "Response: Numerical 3-index-integrals", diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_3index_integration,output_unit, &
        &"Response: Numerical 3-index-integrals: integration",&
        & diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_resp_Coulomb,output_unit, &
        "Response: Transformation of Coulomb integrals",&
        & diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_Coulomb_reading,output_unit, &
        "Response: Transformation of Coulomb integrals: reading tmp file",&
        & diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_Coulomb_writing,output_unit, &
        "Response: Transformation of Coulomb integrals: writing interfacefile",&
        & diff=.false.,sum=.true.,absolut=.false.)
   if ( output_timing_detailedsummary ) &
        call print_timer(timer_resp_Coulomb_trafo,output_unit, &
        "Response: Transformation of Coulomb integrals: transformation",&
        & diff=.false.,sum=.true.,absolut=.false.)

   call print_timer(timer_epe,output_unit,"EPE: Relaxation",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_epe_cycle,output_unit,"EPE: Relaxation cycle",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_epe_forces,output_unit,"EPE: forces",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_get_epe_energy,output_unit,"EPE: get_epe_energy",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_collect_grads,output_unit,"EPE: collect_grads",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_send_epe_ref,output_unit,"EPE: send_epe_ref",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_main_epe,output_unit,"EPE: main_epe",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_latt_epe_forces,output_unit,"EPE: latt_epe_forces",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_defect_contributions,output_unit,"EPE: defect_contrib.",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_collect_def_contributions,output_unit,"EPE: collect.def.",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_defect_ext,output_unit,"EPE: def.ext",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)

   call print_timer(timer_prod_nested,output_unit,"ll_fit_cont : ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_ll_integral,output_unit,"LL_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_ls_integral,output_unit,"LS_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_ss_integral,output_unit,"SS_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_new_ll_integral,output_unit,"LL_new_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_new_ls_integral,output_unit,"LS_new_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_new_ss_integral,output_unit,"SS_new_integral: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_ang3c,output_unit,"LL_ang3c: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)
   call print_timer(timer_rad3c,output_unit,"LL_rad3c: ",&
        & nbr=.true.,diff=.false.,sum=.true.,absolut=.false.)

   call print_time(output_unit,"Time at End of Program")
   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) ""
   end subroutine timer_print_summary
   !*************************************************************


   !*************************************************************
   subroutine timer_print_integral(i_int, int_part_name)
     !
     ! Purpose:  prints summary  of  timing for  actual integral  part
     ! timing for quadrupels and parts  are given for all hosts timing
     ! for  different  tasks  (calculate, receive  integrals,  receive
     ! norms, distribute  new quadrupel,  idle) are kept  separare for
     ! master and slaves
     !
     use time_module, only: timer_type,print_timer,print_time
     use iounitadmin_module, only: output_unit
     use output_module
     use comm, only: comm_size
     implicit none
     integer, intent(in) :: i_int
     character(len=*), intent(in) :: int_part_name
     !** End of interface *****************************************

   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "######################################"
   write(output_unit,*) "## Timing of "//int_part_name//" Integral Part ##"
   write(output_unit,*) "######################################"
   write(output_unit,*) ""
   write(output_unit,*) ""
   call print_timer(timer_int(i_int), output_unit, int_part_name// &
        "Integral Part: Total")

   if (timer_int_2cff(i_int)%nbr_of_calls .gt. 0) then
      write(output_unit,*) ""
      write(output_unit,*) ""
      write(output_unit,*) "### Two center fitfunction part ###"
      write(output_unit,*) ""
      call print_timer(timer_int_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Total")
      write(output_unit,*) ""
      write(output_unit,*) "## Master and Slaves together ##"
      write(output_unit,*) ""
      call print_timer(timer_int_quadrupel_2cff(i_int), output_unit, int_part_name//  &
           "Integral Part: 2 center fitfunction integrals: Quadrupel", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_primcont_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Primitive Integrals and Contraction", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_prim_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Primitive Integrals and Symmetry Adaption", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_cont_2cff(i_int), output_unit, int_part_name//  &
           "Integral Part: 2 center fitfunction integrals: Contractionn and Renormalisation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_communication_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Communication", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_send_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Send", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_sendnorm_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Send Norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_waitnorm_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Wait Norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      write(output_unit,*) ""
      write(output_unit,*) "## Master ##"
      write(output_unit,*) ""
      call print_timer(timer_int_calcsum_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_interrupt_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_idle_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(timer_int_receive_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Receive", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_norm_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_rewrite_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Rewrite", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(timer_int_newtask_2cff(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center fitfunction integrals: Newtask", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)

      if (comm_size() > 1) then
         !
         ! FIXME:  Otherwise it is  not meaningfull,  and, which  is more
         ! important there is a division by zero somewhere in this case.
         !
         write(output_unit,*) ""
         write(output_unit,*) "## Average per slave ##"
         write(output_unit,*) ""
         call print_timer(t_sl_int_calcsum_2cff(i_int), output_unit, int_part_name// &
              "Integral Part: 2 center fitfunction integrals: Calculation", &
              diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
         if ( output_timing_detailedintegrals ) &
              call print_timer(t_sl_int_interrupt_2cff(i_int), output_unit, int_part_name// &
              "Integral Part: 2 center fitfunction integrals: Interrupts", &
              diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
         call print_timer(t_sl_int_idle_2cff(i_int), output_unit, int_part_name// &
              "Integral Part: 2 center fitfunction integrals: Idle", &
              diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
         if ( output_timing_detailedintegrals ) &
              call print_timer(t_sl_int_norm_2cff(i_int), output_unit, int_part_name// &
              "Integral Part: 2 center fitfunction integrals: Norm", &
              diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
         write(output_unit,*) ""
         write(output_unit,*) ""
      endif
   endif

   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "### Two center orbital and three center part ###"
   write(output_unit,*) ""
   call print_timer(timer_int_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Total")
   write(output_unit,*) ""
   write(output_unit,*) "## Master and Slaves together ##"
   write(output_unit,*) ""
   call print_timer(timer_int_quadrupel_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Quadrupel", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_prim_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Primitive Integrals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_cont_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Contraction", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_symadapt_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_symadapt_pvxp_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption of PVXP", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_symadapt_npxp_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Symmetry Adaption except PVXP", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
!!$   call print_timer(timer_int_renorm_2cob3c(i_int), output_unit, int_part_name// &
!!$        "Integral Part: 2 center orbital and 3 center integrals: Renormalisation", &
!!$        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_send_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Send", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_commsend_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: comm Send", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_pack_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Pack", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_write_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Write", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   write(output_unit,*) ""
   write(output_unit,*) "## Master ##"
   write(output_unit,*) ""
   call print_timer(timer_int_calcsum_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Calculation", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_interrupt_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Interrupts", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_idle_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: Idle", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_receive_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: receive", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_norm_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: norm", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_int_rewrite_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: rewrite", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   if ( output_timing_detailedintegrals ) &
        call print_timer(timer_int_newtask_2cob3c(i_int), output_unit, int_part_name// &
        "Integral Part: 2 center orbital and 3 center integrals: newtask", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)

   if (comm_size() > 1) then
      !
      ! FIXME:  Otherwise it is  not meaningfull,  and, which  is more
      ! important there is a division by zero somewhere in this case.
      !
      write(output_unit,*) ""
      write(output_unit,*) "## Average per slave ##"
      write(output_unit,*) ""
      call print_timer(t_sl_int_calcsum_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(t_sl_int_interrupt_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_idle_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_receive_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: receive", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedintegrals ) &
           call print_timer(t_sl_int_norm_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_rewrite_2cob3c(i_int), output_unit, int_part_name// &
           "Integral Part: 2 center orbital and 3 center integrals: rewrite", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   endif

   if (timer_int_rel(i_int)%nbr_of_calls .gt. 0) then
      write(output_unit,*) ""
      write(output_unit,*) ""
      write(output_unit,*) "###  relativistic transformations ###"
      write(output_unit,*) ""
      call print_timer(timer_int_rel(i_int), output_unit, int_part_name// &
           "Integral Part: relativistic transformations: Total")
   endif
   write(output_unit,*) ""
   write(output_unit,*) ""
   call print_time(output_unit,"Time at End of "//int_part_name//" Integral Part")
   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) ""
   end subroutine timer_print_integral
   !*************************************************************


   !*************************************************************
   subroutine timer_print_scfcycle()
   !  Purpose: prints summary of timing of one SCF cycle
     use time_module, only: timer_type,print_timer
     use iounitadmin_module, only: output_unit
     implicit none
   !** End of interface *****************************************

   logical :: extended_mda, xcmda_called

   ! FIXME: what should  slaves do when they dont  have an output file
   ! open? For now bail out:
   if (output_unit <= 0) RETURN

   extended_mda = options_xcmode() == xcmode_extended_mda
   xcmda_called = options_xcmode() == xcmode_model_density .or. extended_mda

   write(output_unit,*) ""
   write(output_unit,*) "Summary of Timing for SCF Cycle"
   write(output_unit,*) ""

   call print_timer(timer_scf_cycle,output_unit, &
        "Entire SCF cycle",start=.true.,sum=.true.,average=.true.,nbr=.true.)

   call print_timer(timer_scf_xc,output_unit, &
        "Numeric Build up of Exchange Part of Hamiltonian")

   call print_timer(timer_scf_ham,output_unit, &
        "Build up of Hamiltonian")

   call print_timer(timer_scf_eigen,output_unit, &
        "Eigensolver")

   call print_timer(timer_scf_reoc,output_unit, &
        "Reoccupation of Orbitals")

   call print_timer(timer_scf_chfit,output_unit, &
        "Chargefit")

   call print_timer(timer_grid_orbitals,output_unit, &
        "Calculation of orbitals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)

   call print_timer(timer_grid_density,output_unit, &
        "Calculation of density", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( xcmda_called ) &
        call print_timer(timer_grid_trunc_dens,output_unit, &
        "Truncation of density", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   call print_timer(timer_grid_functionals,output_unit, &
        "Calculation of xc-functionals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( xcmda_called ) &
        call print_timer(timer_grid_trunc_funcs,output_unit, &
        "Truncation of xc-functionals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( xcmda_called ) &
        call print_timer(timer_grid_projections,output_unit, &
        "Calculation of numerical projections", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( extended_mda ) &
        call print_timer(timer_grid_num_metric,output_unit, &
        "Calculation of numericl metric", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( xcmda_called ) &
        call print_timer(timer_grid_fit_coeffs,output_unit, &
        "Calculation of fitting coefficients", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   if ( .not.xcmda_called ) &
        call print_timer(timer_grid_xcbuild,output_unit, &
        "Construction of the numerical XC-Hamiltonian", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)

   write(output_unit,*) ""
   end subroutine timer_print_scfcycle
   !*************************************************************


   !*************************************************************
   subroutine timer_print_scf ()
   !  Purpose: prints summary of timing of the SCF part
     use time_module, only: timer_type,print_timer
     use iounitadmin_module, only: output_unit
     use output_module, only: output_timing_detailedscf
     implicit none
     !** End of interface *****************************************

!:UB[ added
!  implicit none
   logical :: extended_mda, xcmda_called
!:UB]

   if (output_unit <= 0) return ! FIXME: maybe dont call at all?

!:UB[ added
   extended_mda = options_xcmode() == xcmode_extended_mda
   xcmda_called = options_xcmode() == xcmode_model_density .or. extended_mda
!:UB]

   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "########################"
   write(output_unit,*) "## Timing of SCF Part ##"
   write(output_unit,*) "########################"
   write(output_unit,*) ""
   write(output_unit,*) ""
   call print_timer(timer_scf,output_unit, &
        "Entire SCF part")
   call print_timer(timer_scf_preparations,output_unit, &
        "Preparations in SCF part")
   if ( output_timing_detailedscf ) &
        call print_timer(timer_grid_setup,output_unit, &
        "Grid setup")
   call print_timer(timer_scf_cycle,output_unit, &
        "Entire SCF cycle",diff=.false. &
        ,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
   call print_timer(timer_scf_xc,output_unit, &
        "Numeric Build up of Exchange Part of Hamiltonian", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_ham,output_unit, &
        "Build up of Hamiltonian", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_eigen,output_unit, &
        "Eigensolver", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf ) &
        call print_timer(timer_scf_reoc,output_unit, &
        "Reoccupation of Orbitals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   call print_timer(timer_scf_chfit,output_unit, &
        "Chargefit", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf ) &
        call print_timer(timer_grid_orbitals_sum,output_unit, &
!:UB[ modified
        "Calculation of orbitals", &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf ) &
        call print_timer(timer_grid_density_sum,output_unit, &
!:UB[ modified
        "Calculation of density",  &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedscf .and. xcmda_called ) &
        call print_timer(timer_grid_trunc_dens_sum,output_unit, &
        "Truncation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedscf ) &
        call print_timer(timer_grid_functionals_sum,output_unit, &
!:UB[ modified
        "Calculation of xc-functionals",  &
!:UB]
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedscf .and. xcmda_called ) &
        call print_timer(timer_grid_trunc_funcs_sum,output_unit, &
        "Truncation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf .and. xcmda_called ) &
        call print_timer(timer_grid_projections_sum,output_unit, &
        "Calculation of numerical projections",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf .and. extended_mda ) &
        call print_timer(timer_grid_num_metric_sum,output_unit, &
        "Calculation of numerical metric",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedscf .and. xcmda_called ) &
        call print_timer(timer_grid_fit_coeffs_sum,output_unit, &
        "Calculation of fitting coefficients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
!:UB[ modified
   if ( output_timing_detailedscf .and. .not.xcmda_called) &
!:UB]
        call print_timer(timer_grid_xcbuild_sum,output_unit, &
        "Construction of the numerical XC-Hamiltonian",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   write(output_unit,*) ""
   write(output_unit,*) ""
   end subroutine timer_print_scf
   !*************************************************************


   !*************************************************************
   subroutine timer_print_postscf()
   !  Purpose: prints summary of timing of the PostSCF part
     use time_module, only: timer_type,print_timer
     use iounitadmin_module, only: output_unit
     use output_module, only: output_timing_detailedpostscf
     implicit none
!:UB[ added
   logical :: extended_mda, xcmda_called
!:UB]
   !** End of interface *****************************************
!:UB[ added
   extended_mda = options_xcmode() == xcmode_extended_mda
   xcmda_called = options_xcmode() == xcmode_model_density .or. extended_mda
!:UB]
     !** End of interface *****************************************
   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "############################"
   write(output_unit,*) "## Timing of PostSCF Part ##"
   write(output_unit,*) "############################"
   write(output_unit,*) ""
   write(output_unit,*) ""
   call print_timer(timer_gridph,output_unit, &
        "Post SCF: Complete grid_part")
   if ( output_timing_detailedpostscf ) &
        call print_timer(timer_gridph_orbitals,output_unit, &
        "Post scf: calculation of orbitals", &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedpostscf ) &
        call print_timer(timer_gridph_density,output_unit, &
        "Post scf: calculation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedpostscf .and. xcmda_called ) &
        call print_timer(timer_gridph_trunc_dens,output_unit, &
        "Post scf: truncation of density",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedpostscf ) &
        call print_timer(timer_gridph_functionals,output_unit, &
        "Post scf: calculation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB[ added
   if ( output_timing_detailedpostscf .and. xcmda_called ) &
        call print_timer(timer_gridph_trunc_funcs,output_unit, &
        "Post scf: truncation of xc-functionals",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedpostscf .and. xcmda_called ) &
        call print_timer(timer_gridph_projections,output_unit, &
        "Post scf:  calculation of projections",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedpostscf .and. xcmda_called ) &
        call print_timer(timer_gridph_fit_coeffs,output_unit, &
        "Post scf: calculation of fit coefficients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
!:UB]
   if ( output_timing_detailedpostscf ) &
        call print_timer(timer_gridph_grdwts,output_unit, &
        "Post scf: calculation of grdwts",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   if ( output_timing_detailedpostscf ) &
        call print_timer(timer_gridph_gradients,output_unit, &
        "Post scf: calculation of gradients",  &
        diff=.false.,absolut=.false.,sum=.true.,average=.true.)
   write(output_unit,*) ""
   write(output_unit,*) ""
   end subroutine timer_print_postscf
   !*************************************************************


   !*************************************************************
   subroutine timer_print_response(hyper,method_RI1)
     !  Purpose: prints summary of timing of the Response part
     use time_module, only: timer_type,print_timer
     use iounitadmin_module, only: output_unit
     implicit none
     logical,intent(in) :: hyper, method_RI1
     !** End of interface *****************************************
     write(output_unit,*) ""
     write(output_unit,*) ""
     write(output_unit,*) "#############################"
     write(output_unit,*) "## Timing of Response Part ##"
     write(output_unit,*) "#############################"
     write(output_unit,*) ""
     write(output_unit,*) ""
     call print_timer(timer_resp,output_unit, "Entire RESPONSE part")
     call print_timer(timer_resp_preparations,output_unit, &
          "Response: Preparations", diff=.false.,sum=.true.,absolut=.false.)
     call print_timer(timer_resp_setup,output_unit, &
          "Response: Setup", diff=.false.,sum=.true.,absolut=.false.)
     call print_timer(timer_resp_header,output_unit, &
          "Response: Writing header files", diff=.false.,sum=.true.,absolut=.false.)
     if(.NOT.hyper) then
        call print_timer(timer_resp_dipole,output_unit, &
             "Response: Rewriting transition dipole moments",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_2index,output_unit, &
             "Response: Numerical 2-index-integrals", diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_2index_integration,output_unit, &
             &"Response: Numerical 2-index-integrals: integration",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_3index,output_unit, &
             "Response: Numerical 3-index-integrals", diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_3index_integration,output_unit, &
             &"Response: Numerical 3-index-integrals: integration",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_Coulomb,output_unit, &
             "Response: Transformation of Coulomb integrals",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_Coulomb_reading,output_unit, &
             "Response: Transformation of Coulomb integrals: reading tmp file",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_Coulomb_writing,output_unit, &
             "Response: Transformation of Coulomb integrals: writing interfacefile",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_Coulomb_trafo,output_unit, &
             "Response: Transformation of Coulomb integrals: transformation",&
             & diff=.false.,sum=.true.,absolut=.false.)
     else
        call print_timer(timer_resp_allClb,output_unit, &
             "Response: Transformation of all Coulomb integrals",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_calcFG,output_unit, &
             "Response: Calculation of 2-index matrices F and G",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_addCoul,output_unit, &
             "Response: Evaluation of Coulomb contribution to matrix F",&
             & diff=.false.,sum=.true.,absolut=.false.)
        if(method_RI1) call print_timer(timer_resp_all3ind,output_unit, &
             "Response: Calculation of all numerical 3-index integrals",&
             & diff=.false.,sum=.true.,absolut=.false.)
        call print_timer(timer_resp_rhsh,output_unit, &
             "Response: Evaluation of the final righthandside",&
             & diff=.false.,sum=.true.,absolut=.false.)
     end if
     write(output_unit,*) ""
     write(output_unit,*) ""
   end subroutine timer_print_response
   !*************************************************************


   !*************************************************************
   subroutine timer_grid_small_to_large()
   !  Purpose: Stores timing on grid for one SCF cycle
     use time_module, only: timer_type,timer_small_to_large
     implicit none
     !** End of interface *****************************************
!:UB[ added
   logical :: extended_mda, xcmda_called
!:UB]
   !------------ Executable code --------------------------------
!:UB[ added
   extended_mda = options_xcmode() == xcmode_extended_mda
   xcmda_called = options_xcmode() == xcmode_model_density .or. extended_mda
!:UB]
   call timer_small_to_large(timer_grid_orbitals,timer_grid_orbitals_sum)
   call timer_small_to_large(timer_grid_density,timer_grid_density_sum)
!:UB[ added
   if ( xcmda_called ) &
   call timer_small_to_large(timer_grid_trunc_dens,timer_grid_trunc_dens_sum)
!:UB]
   call timer_small_to_large(timer_grid_functionals,timer_grid_functionals_sum)
!:UB[ added
   if ( xcmda_called ) &
   call timer_small_to_large(timer_grid_trunc_funcs,timer_grid_trunc_funcs_sum)
   if ( xcmda_called ) &
   call timer_small_to_large(timer_grid_projections,timer_grid_projections_sum)
   if ( extended_mda ) &
   call timer_small_to_large(timer_grid_num_metric,timer_grid_num_metric_sum)
   if ( xcmda_called ) &
   call timer_small_to_large(timer_grid_fit_coeffs,timer_grid_fit_coeffs_sum)
!:UB]
!AG
   if ( xcmda_called ) &
   call timer_small_to_large(timer_int_fit_density_calc,timer_int_fit_density_calc_sum)
!AG
!:UB[ modified
   if ( .not.xcmda_called ) &
   call timer_small_to_large(timer_grid_xcbuild,timer_grid_xcbuild_sum)
!:UB]
   end subroutine timer_grid_small_to_large
   !*************************************************************

   !*************************************************************
   subroutine timer_send_slave_int_timing(i_int)
   !  Purpose: sends timing of integral part of one slave to master
     use time_module, only: timer_type,pack_timer
     use comm_module
     use msgtag_module, only: msgtag_slavetiming
     implicit none
     integer, intent(in) :: i_int
     !** End of interface *****************************************

   call comm_init_send(comm_master_host, msgtag_slavetiming)
   call pack_timer(timer_int(i_int))
   call pack_timer(timer_int_2cff(i_int))
   call pack_timer(timer_int_interrupt_2cff(i_int))
   call pack_timer(timer_int_calcsum_2cff(i_int))
   call pack_timer(timer_int_primcont_2cff(i_int))
   call pack_timer(timer_int_prim_2cff(i_int))
   call pack_timer(timer_int_cont_2cff(i_int))
   call pack_timer(timer_int_idle_2cff(i_int))
   call pack_timer(timer_int_send_2cff(i_int))
   call pack_timer(timer_int_sendnorm_2cff(i_int))
   call pack_timer(timer_int_waitnorm_2cff(i_int))
   call pack_timer(timer_int_communication_2cff(i_int))
   call pack_timer(timer_int_quadrupel_2cff(i_int))
   call pack_timer(timer_int_norm_2cff(i_int))
   call pack_timer(timer_int_2cob3c(i_int))
   call pack_timer(timer_int_interrupt_2cob3c(i_int))
   call pack_timer(timer_int_calcsum_2cob3c(i_int))
   call pack_timer(timer_int_prim_2cob3c(i_int))
   call pack_timer(timer_int_cont_2cob3c(i_int))
   call pack_timer(timer_int_symadapt_2cob3c(i_int))
   call pack_timer(timer_int_symadapt_pvxp_2cob3c(i_int))
   call pack_timer(timer_int_symadapt_npxp_2cob3c(i_int))
   call pack_timer(timer_int_idle_2cob3c(i_int))
   call pack_timer(timer_int_write_2cob3c(i_int))
   call pack_timer(timer_int_send_2cob3c(i_int))
   call pack_timer(timer_int_commsend_2cob3c(i_int))
   call pack_timer(timer_int_pack_2cob3c(i_int))
   call pack_timer(timer_int_quadrupel_2cob3c(i_int))
   call pack_timer(timer_int_receive_2cob3c(i_int))
   call pack_timer(timer_int_norm_2cob3c(i_int))
   call pack_timer(timer_int_rewrite_2cob3c(i_int))
   call pack_timer(timer_int_rel(i_int))
   call comm_send()
   end subroutine timer_send_slave_int_timing
   !*************************************************************


   !*************************************************************
   subroutine timer_receive_slave_int_timing(i_int)
   ! Purpose: receive timings of slves at end of integral part.
   ! timing for quadrupels and parts are collected for all hosts
   ! timing for different tasks (calculate, receive integrals,
   ! receive norms, distribute new quadrupel, idle) are kept
   ! separare for master and slaves
     use time_module, only: timer_type,unpack_add_timer,divide_timer
     use comm_module, only: comm_get_n_processors, comm_save_recv, &
          comm_all_other_hosts
     use msgtag_module, only: msgtag_slavetiming
     implicit none
     integer, intent(in) :: i_int
     !** End of interface *****************************************

     integer (i4_kind) :: i, np
     real (r8_kind) :: n_slaves
   !------------ Executable code --------------------------------

     np = comm_get_n_processors()
   do i = 1, np - 1
      call comm_save_recv (comm_all_other_hosts, msgtag_slavetiming)
      call unpack_add_timer(t_sl_int(i_int))
      call unpack_add_timer(t_sl_int_2cff(i_int))
      call unpack_add_timer(t_sl_int_interrupt_2cff(i_int))
      call unpack_add_timer(t_sl_int_calcsum_2cff(i_int))
      call unpack_add_timer(timer_int_primcont_2cff(i_int))
      call unpack_add_timer(timer_int_prim_2cff(i_int))
      call unpack_add_timer(timer_int_cont_2cff(i_int))
      call unpack_add_timer(t_sl_int_idle_2cff(i_int))
      call unpack_add_timer(timer_int_send_2cff(i_int))
      call unpack_add_timer(timer_int_sendnorm_2cff(i_int))
      call unpack_add_timer(timer_int_waitnorm_2cff(i_int))
      call unpack_add_timer(timer_int_communication_2cff(i_int))
      call unpack_add_timer(timer_int_quadrupel_2cff(i_int))
      call unpack_add_timer(t_sl_int_norm_2cff(i_int))
      call unpack_add_timer(t_sl_int_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_interrupt_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_calcsum_2cob3c(i_int))
      call unpack_add_timer(timer_int_prim_2cob3c(i_int))
      call unpack_add_timer(timer_int_cont_2cob3c(i_int))
      call unpack_add_timer(timer_int_symadapt_2cob3c(i_int))
      call unpack_add_timer(timer_int_symadapt_pvxp_2cob3c(i_int))
      call unpack_add_timer(timer_int_symadapt_npxp_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_idle_2cob3c(i_int))
      call unpack_add_timer(timer_int_write_2cob3c(i_int))
      call unpack_add_timer(timer_int_send_2cob3c(i_int))
      call unpack_add_timer(timer_int_commsend_2cob3c(i_int))
      call unpack_add_timer(timer_int_pack_2cob3c(i_int))
      call unpack_add_timer(timer_int_quadrupel_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_receive_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_norm_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_rewrite_2cob3c(i_int))
      call unpack_add_timer(t_sl_int_rel(i_int))
   enddo
   ! calculate average for one slave
   n_slaves = real (np - 1, r8_kind)
   if (n_slaves > 0.0) then
      call divide_timer( t_sl_int(i_int), n_slaves )
      call divide_timer( t_sl_int_2cff(i_int), n_slaves )
      call divide_timer( t_sl_int_interrupt_2cff(i_int), n_slaves )
      call divide_timer( t_sl_int_calcsum_2cff(i_int), n_slaves )
      call divide_timer( t_sl_int_idle_2cff(i_int), n_slaves )
      call divide_timer( t_sl_int_norm_2cff(i_int), n_slaves )
      call divide_timer( t_sl_int_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_interrupt_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_calcsum_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_idle_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_receive_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_norm_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_rewrite_2cob3c(i_int), n_slaves )
      call divide_timer( t_sl_int_rel(i_int), n_slaves )
   endif
   end subroutine timer_receive_slave_int_timing
   !*************************************************************

   subroutine timer_gather_slave_int_timing(i_int)
     use comm, only: comm_rank
     implicit none
     integer(i4_kind), intent(in) :: i_int
     ! *** end of interface ***

     if ( comm_rank() == 0 ) then
        call timer_receive_slave_int_timing(i_int)
     else
        call timer_send_slave_int_timing(i_int)
     endif
   end subroutine timer_gather_slave_int_timing

   !*************************************************************
   subroutine timer_print_slavetiming (integralpar_int_part_name)
   !  Purpose: prints summary of timing of slaves
     use time_module, only: timer_type,print_timer
     use iounitadmin_module, only: output_unit
     use output_module
     implicit none
     character(len=*), intent(in) :: integralpar_int_part_name(:)
     !** End of interface *****************************************

   integer :: i_int
   !------------ Executable code --------------------------------

   if (output_unit <= 0) return ! FIXME: maybe dont call at all?

   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) "#################################"
   write(output_unit,*) "## Summary of Timing of Slaves ##"
   write(output_unit,*) "#################################"
   write(output_unit,*) ""
   write(output_unit,*) ""
   do i_int = 1, size(t_sl_int) ! n_int_parts
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int(i_int),output_unit, &
           "Integral Part: Total", &
           diff=.false.,absolut=.false.,sum=.true.,average=.false.,nbr=.false.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_2cff(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Total", &
           diff=.false.,absolut=.false.,sum=.true.,average=.false.,nbr=.false.)
      call print_timer(t_sl_int_calcsum_2cff(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_interrupt_2cff(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_idle_2cff(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_norm_2cff(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center fitfunction integrals: Norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Total", &
           diff=.false.,absolut=.false.,sum=.true.,average=.false.,nbr=.false.)
      call print_timer(t_sl_int_calcsum_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Calculation", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_interrupt_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Interrupts", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_idle_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: Idle", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_receive_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: receive", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_norm_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: norm", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      call print_timer(t_sl_int_rewrite_2cob3c(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: 2 center orbital and 3 center integrals: rewrite", &
           diff=.false.,absolut=.false.,sum=.true.,average=.true.,nbr=.true.)
      if ( output_timing_detailedsummary ) &
           call print_timer(t_sl_int_rel(i_int),output_unit, integralpar_int_part_name(i_int)// &
           "Integral Part: relativistic transformations: Total", &
           diff=.false.,absolut=.false.,sum=.true.,average=.false.,nbr=.false.)
   enddo
   write(output_unit,*) ""
   write(output_unit,*) ""
   write(output_unit,*) ""
   end subroutine timer_print_slavetiming
   !*************************************************************


end module timer_module

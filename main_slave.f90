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
subroutine main_slave()
!
! This subroutine tries to receive messages from master in an infinite
! loop. Depending  on msgtag  of received message,  it goes  to switch
! statement of  all known msgtags to  decide what the  slave should do
! and  invokes  appropriate methods  of  modules  or subroutines  that
! unpack the  received data and perform  calculations.  The subroutine
! returns when receiving messagetag msgtag_finito.
!
! Subroutine called by: main
!
! Author: Folke
! Date:   7/95
!
!================================================================
! End of public interface of module
!================================================================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: pvm -> comm
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

#include "def.h"
use type_module
use comm_module, only: comm_save_recv, comm_msgtag, &
     comm_master_host, comm_any_message
use msgtag_module
use density_data_module, only: density_data_free
use pointcharge_module
use induced_dipoles_module, only: send_receive_id
use calc_id_module, only: calc_Pol_ham
use fit_coeff_module, only: fit_coeff_receive
#ifdef WITH_RESPONSE
  use response_module, only: response_setup, &
       & response_calc_2index_int_v2
  use int_resp_module, only: int_resp_Clb_3c
  use init_tddft_module, only: init_tddft_start
  use tddft_diag, only: diag_init
  use global_module, only: global_dealloc_M
  use int_send_2c_resp, only: int_send_2c_resp_rewrite
  use noRI_module, only: noRI_2c
#endif

#ifdef WITH_EPE
use epe_module
use main_epe_module, only : main_epe,                      &
                            epe_receive_shells_and_cores,  &
                            cons_latt_gradients,           &
                            init_lat_gradients,            &
                            epe_init_slave,                &
                            epe_finish_slave,              &
                            cons_defect_contrubutions,     &
                            defect_contributions,          &
                            defect_contributions_fin
#endif
use potential_module, only: bounds_free_poten, start_read_poten_e
use elec_static_field_module, only: receive_surf_point, &
     surf_points_gradinfo_dealloc,surf_points_grad_information,read_field_e,send_receive_field, &
     bounds_free_field
use interfaces, only: ISLAV
use calc3c_switches
#ifdef WITH_MOLMECH
use molmech_msgtag_module, only: msgtag_start_molmech
use qmmm_interface_module, only: slave_run
#endif

implicit none
integer (i4_kind) :: msgtag

do ! while comm_msgtag() /= msgtag_finito, then RETURN

   call comm_save_recv(comm_master_host, comm_any_message)

   msgtag = comm_msgtag()
   DPRINT 'main_slave: received msgtag=', msgtag
   select case (msgtag)
   case( msgtag_free_bnds_ptn )
      call say("free_bounds_poten")
      call bounds_free_poten()
   case (msgtag_start_poten)
      call say("start_read_poten_e")
      call start_read_poten_e()
   case (msgtag_density_data_free)
      ! FIXME: density_data_free1() called from
      ! elec_pot_and_field/potential_calc_module.f90 uses this tag.
      call say("density_data_free")
      call density_data_free()
   case(msgtag_fit_coeff_send)
      call say("fit_coeff_receive")
      call fit_coeff_receive()
#ifdef WITH_RESPONSE
     case(msgtag_response_setup)
        call say("response_setup")
        call response_setup()
     case(msgtag_response_2index)
        call say("response_calc_2index_int_v2")
        call response_calc_2index_int_v2()
     case(msgtag_response_3Clb_start)
        call say("int_resp_Clb_3c")
        call int_resp_Clb_3c()
     case(msgtag_tddft_eps_eta)
        call say("init_tddft_eps_eta")
        call init_tddft_start()
     case(msgtag_tddft_clshell)
        call say("closed_shell_init")
        call diag_init()
     case(msgtag_tddft_dealloM)
        call say("global_dealloc_M")
        call global_dealloc_M()
     case(msgtag_send_2c_start)
        call say("int_send_2c_resp_rewrite")
        call int_send_2c_resp_rewrite()

     case(msgtag_nori_2c_send)
        call say("noRI_2c start")
        call noRI_2c()
#endif

#ifdef WITH_EPE
!AG ============================ EPE-distribution =========================
  case(msgtag_epe_data_sent)
     call say("receiving epe_data")
     call epe_receive_data()
  case(msgtag_epe_do_gradients)
     call say("receiving reference")
     call epe_receive_reference()
     call say("epe_field_and_forces")
     call epe_field_and_forces_par()
  case(msgtag_epe_init_slave)
     call say("epe_init")
     call epe_init_slave()
  case(msgtag_epe_grads_init)
     call say("receiving shell and cores (init)")
     call epe_receive_shells_and_cores()
     call init_lat_gradients()
  case(msgtag_epe_grads_cons)
     call say("receiving shell and cores (cons)")
     call epe_receive_shells_and_cores()
     call cons_latt_gradients()
  case(msgtag_epe_send_only)
     call say("receive shels and cores (ONLY)")
     call epe_receive_shells_and_cores()
  case(msgtag_epe_consdef)
     call say("begin cons.defect contributions")
     call cons_defect_contrubutions()
  case(msgtag_epe_defects)
     call say("begin defect contributions")
     call defect_contributions()
  case(msgtag_epe_def_fin)
     call say("begin defect contrib.fin")
     call defect_contributions_fin()
  case(msgtag_epe_finish_slave)
     call say("epe_finish")
     call epe_finish_slave()
!AG ========================= end of EPE-distribution ======================
#endif
  case (msgtag_surf_point)
     call say("receive_surf_point")
     call receive_surf_point()
  case (msgtag_surf_point_sa)
     call say("surf_points_symadapt_unpack")
     call surf_points_grad_information()
  case (msgtag_sp_grinfo_dealloc)
     call say("surf_points_gradinfo_dealloc")
     call surf_points_gradinfo_dealloc()
   case( msgtag_start_field )
      call say("read_field")
      call read_field_e
      call send_receive_field
   case( msgtag_free_bnds_fld )
      call say("free_bounds_field")
      call bounds_free_field()
   case(msgtag_ind_dipmom)
      call say("send_receive_id")
      call send_receive_id()
   case(msgtag_Pol_ham)
      call say("calc_Pol_ham")
      call calc_Pol_ham()

#ifdef WITH_MOLMECH
  case(msgtag_start_molmech)
     ! molecular mechanics part
     call main_molmech(slave_run)
#endif
  case (msgtag_finito)
     call say("exiting")
     return
  case default
     print *,'main_slave: received msgtag=',msgtag,': ',msgtag_name(msgtag)
     call error_handler('main_slave: wrong message tag: '//msgtag_name(msgtag))
  end select

enddo
   !
   ! FIXME:   this   section  is   not   reachable,   see  RETURN   on
   ! msgtag_finito!   Gfortran was  smart  enough to  notice that  and
   ! optimize away the call to print_alloc_grid() from here.

contains

  subroutine say(phrase)
    use output_module, only: output_slaveoperations
    use iounitadmin_module, only: write_to_output_units
    implicit none
    character(len=*), intent(in) :: phrase
    ! *** end of interface ***

    if( output_slaveoperations ) then
        call write_to_output_units("main_slave: "//phrase)
    endif
  end subroutine say

end subroutine main_slave


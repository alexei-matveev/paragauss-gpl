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
  !===================================================================
  ! End of public interface of module
  !===================================================================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
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
!---------------------------------------------------------------------

#include "def.h"
use type_module
use comm_module, only: comm_save_recv, comm_msgtag, &
     comm_master_host, comm_any_message
use msgtag_module
use pointcharge_module
use induced_dipoles_module, only: send_receive_id
use calc_id_module, only: calc_Pol_ham
use fit_coeff_module, only: fit_coeff_receive
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
use elec_static_field_module, only: receive_surf_point, &
     surf_points_gradinfo_dealloc,surf_points_grad_information,read_field_e,send_receive_field, &
     bounds_free_field
use calc3c_switches

implicit none
integer (i4_kind) :: msgtag

do ! while comm_msgtag() /= msgtag_finito, then RETURN

   call comm_save_recv(comm_master_host, comm_any_message)

   msgtag = comm_msgtag()
   DPRINT 'main_slave: received msgtag=', msgtag
   select case (msgtag)
   case (msgtag_fit_coeff_send)
      !
      ! FIXME: this is here  because of the unmatched fit_coeff_send()
      ! main_epe_block() only. Make it run in parallel context and get
      ! rid of this:
      !
      call say("fit_coeff_receive")
      call fit_coeff_receive()

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


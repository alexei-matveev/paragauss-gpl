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
module molmech_slave_module
#include "def.h"
  !------------ Modules used --------------------------------------
  use type_module
  use comm_module, only: &
       comm_save_recv,   &
       comm_master_host, &
       comm_any_message, &
       comm_msgtag
  use molmech_msgtag_module
  use tasks_main_options_module, only : &
       send_receive_tasks_and_options,  &
       calc_gradients,calc_hessian
  use slab_module, only : &
       slab_calc
  use species_module, only :       &
       send_receive_species,       &
       shutdown_species_on_slaves, &
       lattice_calc 
  use ewald_module, only : &
       send_receive_ewald, &
       recipr_ew_E_and_F,  &
       calc_gauss,         &
       shutdown_ewald
  use ewald2d_module, only : &
       send_receive_ewald2d, &
       calc_2d_slab_vec,     &
       recipr_ew2d_E_self,   &
       recipr_ew2d_E_and_F,  &
       shutdown_ewald2d
  use energy_and_forces_module, only : &
       init_energy_and_forces,         &
       shutdown_gradients,             &
       send_receive_e_g_h,             &
       send_receive_e_g_h_cascad
  use hess_and_opt_module, only : &
       alloc_H_on_slaves,         &
       shutdown_H_on_slaves
  use n_body_lists_module, only : &
       send_receive_n_total,      &
       shutdown_interacting_lists

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ------------------
  public molmech_slave
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine molmech_slave()

    integer(i4_kind) :: msg_tag

    do ! while comm_msgtag() /= msgtag_mm_shutdown, then RETURN

       call comm_save_recv(comm_master_host, comm_any_message)

       msg_tag = comm_msgtag()
       DPRINT 'main_molmech_slave: received msgtag=',msg_tag

       select case(msg_tag)
       case(msgtag_mm_initslaves)
          call send_receive_tasks_and_options()
       case(msgtag_mm_send_ntotal)
          call send_receive_n_total()
       case(msgtag_mm_send_species)
          call send_receive_species()
          call init_energy_and_forces()
          if(calc_hessian) call alloc_H_on_slaves()
       case(msgtag_mm_send_ewald)
          call send_receive_ewald()
          call calc_gauss()
       case(msgtag_mm_send_ewald2d)
          call send_receive_ewald2d()
          call calc_2d_slab_vec()
          call recipr_ew2d_E_self()
       case(msgtag_mm_calc_start)
          if(lattice_calc) then
             call recipr_ew_E_and_F()
             call shutdown_ewald()
          else if(slab_calc) then
             call recipr_ew2d_E_and_F()
             call shutdown_ewald2d()
          end if
!!$       call send_receive_e_g_h
          call send_receive_e_g_h_cascad
          call shutdown_species_on_slaves()
          if(calc_gradients) call shutdown_gradients()
          if(calc_hessian) call shutdown_H_on_slaves()
       case(msgtag_mm_shutdown)
          call shutdown_interacting_lists()
          return
       end select
    end do

  end subroutine molmech_slave
  !****************************************************************
end module molmech_slave_module

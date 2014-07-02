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
subroutine main_molmech (job, iwork, geo_steps)
  !
  ! Executed by main_master(). Runs on all workers.
  !
  use type_module
  use comm, only: comm_rank
!!$  use common_data_module, only: infinity,ten
  use inp_out_module
  use mm_timer_module
  use tasks_main_options_module
  use slab_module, only: slab_calc
  use species_module
  use external_field_module
  use n_body_lists_module
  use energy_and_forces_module
  use potentials_module
  use covalent_module
  use van_der_waals_module
  use coulomb_module
  use ewald_module
  use ewald2d_module
  use calc_energy_module
  use hess_and_opt_module
  use comm_module, only: comm_all_other_hosts, comm_init_send,comm_send
  use molmech_msgtag_module, only: msgtag_mm_shutdown
  use pc_array_module
  use qmmm1_interface_module, only: mm_grads_to_qmmm1
  use molmech_slave_module, only: molmech_slave !!!!!!!!!!

  !=========================================================
  implicit none
  save
  !=========================================================
  integer (i4_kind) :: job
  integer (i4_kind), optional :: iwork, geo_steps

  !*********************************************************

  ! Slaves  will  spin  in  molmech_slave() waiting  for  orders  from
  ! master:
  if (comm_rank() /= 0) then
     ! Exit upon receival of  msgtag_mm_shutdown sent by master at the
     ! end of this sub:
     call molmech_slave()
     return
  endif

  call init_timer()
  call start_mm_timer(full_time)

  qmmm=job
  
  if(qmmm == 4) goto 4
  
  call read_molmech_input()
!!$print*,'READ INPUT DONE'
  if(qmmm == -1) return

  iwork=0
  geo_steps=1
  if(with_optimizer) then 
     iwork=int(-gx_work,i4_kind)
     geo_steps=n_iterations
  end if

4 if(calc_energy) then
!!$  call species_resorting()
!!$  call resorting_Nb_lists()

     if(coulomb) then 
        if(trim(coulomb_type) == "DIPOLE-DIPOLE") then
           call dip_atom_pairs()
        end if
     end if

     call start_mm_timer(nbl_time)
     call build_2b_nonbonded_lists(.true.,.true.)
     call stop_mm_timer(nbl_time)
!!$print*,'N-BODY LISTS DONE'

     call init_energy_and_forces() 
!!$print*,'INIT ENERGY AND GRAD DONE'

     if(calc_optimization .or. calc_hessian) then
        if(calc_optimization) then 
           call init_hess()
           call minimize()
           if(to_pg_hess) then
              calc_optimization=.false.
              call calc_hess()
           end if
        else
           call calc_hess()
        end if
     else
        call total_energy_and_grad(do_corr=.true.)
     end if
!!$print*,'DONE MAIN TASK'

     if(calc_optimization) call write_final_geometry()
     if(calc_strain .and. calc_gradients) call print_lat_param_and_grad()
     call write_energy_and_gradients()
!!$print*,'DONE FINAL_GEOM'

     if(to_xyz_file .or. to_xmkmol_file) call save_xyz()
     if(to_viewmol) call save_vm()
     if(write_gx .and. calc_gradients .and. read_gx) call write_gxfile()
     if(qmmm==2 .or. qmmm==3) call grads_to_qmmm()
     if(qmmm==4) call mm_grads_to_qmmm1()
     if(to_epe .and. (lattice_calc .or. slab_calc)) call epe_interface()
!!$print*,'ALL DONE'
  end if

  if(calc_ppc_array) then
     call gen_pc_array()
  end if

  if(calc_energy) then
     call shutdown_gradients()
     call shutdown_species_data()
     if(ext_field .and. ext_pc) call shutdown_ext_pc()
     call shutdown_poten_data()
     call shutdown_interacting_lists()
     call shutdown_electrostatic()
     call shutdown_ewald()
     call shutdown_ewald2d()
     call shutdown_hess()
  end if

  call start_mm_timer(comm_time)

  ! Tells slaves to exit molmech_slave():
  if (comm_rank() == 0) then
    call comm_init_send (comm_all_other_hosts, msgtag_mm_shutdown)
    call comm_send()
  endif
  call stop_mm_timer(comm_time)

  call stop_mm_timer(full_time)
  call print_timing()

  call close_output_device()
!!$print*,'FINISH'

end subroutine main_molmech





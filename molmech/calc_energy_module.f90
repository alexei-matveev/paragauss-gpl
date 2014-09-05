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
module calc_energy_module
  !------------ Modules used -----------------------------------------
  use type_module
  use common_data_module
  use tasks_main_options_module
  use slab_module
  use species_module
  use n_body_lists_module
  use covalent_module
  use van_der_waals_module
  use coulomb_module
  use ewald_module
  use ewald2d_module
  use energy_and_forces_module
  use molmech_msgtag_module
  use comm_module
  use mm_timer_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of public constants and variables -----
  
  !------------ public functions and subroutines ---------------------
  public total_energy_and_grad
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ------------------------------------------
contains
  !****************************************************************
  subroutine total_energy_and_grad(do_corr)

    logical :: do_corr
    logical :: do_ewald
    real(r8_kind) :: cutr
    integer(i4_kind) :: i

    do_ewald=(coulomb .and. trim(coulomb_type) =="EWALD")

    call start_mm_timer(single_point_time)

    if(comm_parallel() .and. do_ewald) then
       call start_mm_timer(comm_time)
       call send_receive_tasks_and_options()
       call stop_mm_timer(comm_time)
    end if

    E_total=zero
    E=zero
    E_coulomb=zero
    E_solv_tot=zero
    E_solv_el=zero
    E_solv_dr=zero
    E_solv_cav=zero
    if(calc_gradients) then
       Grad=zero
       if(calc_strain)  Grad_s=zero
    end if
    E_ew_d=zero
    E_ew_r=zero

    if(do_ewald) then
       if(lattice_calc) then
          call shutdown_ewald()
          call init_ewald()
       else if(slab_calc) then
          call shutdown_ewald2d()
          call init_ewald2d()
       end if
    end if
!!$    print*,'INIT  DONE'

    if(do_corr) then
       if(lattice_calc) then
          call frac2cart(n_species,vect,atoms_frac,atoms_cart)
          cutr=max(rcut_vdw,rcut_ew)
          if(.not.minimal_image) then
             call start_mm_timer(cis_time)
             call calc_images_species(cutr)
             call stop_mm_timer(cis_time)
          end if
       else if(slab_calc) then
          do i=1,n_species
!!$             call cart2frac_slab(atoms_cart(i)%r,atoms_frac(i)%r)
             call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
          end do
          cutr=max(rcut_vdw,rcut_ew)
          if(.not.minimal_image) then
             call start_mm_timer(cis_time)
             call calc_images_slab_species(cutr)
             call stop_mm_timer(cis_time)
          end if
       end if
!!$       if(cutr /= infinity) then
          call start_mm_timer(nbl_time)
          call build_2b_nonbonded_lists(.false.,.false.)
          call stop_mm_timer(nbl_time)
!!$       end if
    end if
!!$    print*,'CORRECTION  DONE'

    if(comm_parallel() .and. do_ewald) then
       call start_mm_timer(comm_time)
       call send_receive_species()
       if(coulomb.and.slab_calc) call send_receive_n_total()
       call stop_mm_timer(comm_time)
    end if

    if(coulomb) then
       if(trim(coulomb_type) =="DIRECT") then
          call start_mm_timer(coulomb_time)
          call direct_coulomb_E_and_F()
          call stop_mm_timer(coulomb_time)
       else if(trim(coulomb_type) =="DIPOLE-DIPOLE") then
          call start_mm_timer(dipdip_time)
          call dipol_dipol_E_and_F()
          call stop_mm_timer(dipdip_time)
       else if(trim(coulomb_type) =="EWALD") then
          call start_mm_timer(ewald_r_time)
          if(lattice_calc) then
             call start_mm_timer(comm_time)
             call send_receive_ewald()
             call stop_mm_timer(comm_time)
             if(comm_parallel()) then
                call start_mm_timer(comm_time)
                call comm_init_send(comm_all_other_hosts,msgtag_mm_calc_start)
                call comm_send()
                call stop_mm_timer(comm_time)
             end if
             call calc_gauss()
             call recipr_ew_E_and_F()
!!$             call recipr_ew_E_and_F_t()
          else if(slab_calc) then
             call start_mm_timer(comm_time)
             call send_receive_ewald2d()
             call stop_mm_timer(comm_time)
             if(comm_parallel()) then
                call start_mm_timer(comm_time)
                call comm_init_send(comm_all_other_hosts,msgtag_mm_calc_start)
                call comm_send()
                call stop_mm_timer(comm_time)
             end if
             call calc_2d_slab_vec()
             call recipr_ew2d_E_self()
             call recipr_ew2d_E_and_F()
          end if
          call stop_mm_timer(ewald_r_time)

          call start_mm_timer(ewald_d_time)
          if(lattice_calc) then
             if(do_corr) then
                call self_ewald_energy()
             end if
             call direct_ew_E_and_F()
          else if(slab_calc) then
             call self_ewald2d_energy()
             call direct_ew2d_E_and_F()
          end if
          call stop_mm_timer(ewald_d_time)
       end if
!!$       print*,'COULOMB DONE'
    end if

    if(comm_parallel() .and. do_ewald) then
       call start_mm_timer(comm_time)
!!$       call send_receive_e_g_h()
       call send_receive_e_g_h_cascad()
       call stop_mm_timer(comm_time)
    end if

    call start_mm_timer(ebonding_time)
    call two_body_E_and_F()
!!$    print*,'2-BODY  DONE'
    call three_body_E_and_F()
!!$    print*,'3-BODY  DONE'
    call four_body_E_and_F()
!!$    print*,'4-BODY  DONE'
    call many_body_E_and_F()
!!$    print*,'MANY-BODY  DONE'
    call core_shell_E_and_F()
!!$    print*,'CORE-SHELL DONE'
    call stop_mm_timer(ebonding_time)
    if(van_der_waals) then
       call start_mm_timer(vdw_time)
       call vdw_E_and_F()
       call stop_mm_timer(vdw_time)
!!$    print*,'VDW DONE'
    end if

    if(solvent) then
!!$grad=zero
!!$E_total=zero
       call calc_solv_eff()
!!$       print*,'SOLV DONE'
    end if

    call stop_mm_timer(single_point_time)

  end subroutine total_energy_and_grad

end module calc_energy_module

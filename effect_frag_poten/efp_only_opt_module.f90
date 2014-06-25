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
module efp_only_opt_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
#include "def.h"
  use type_module ! type specification parameters
  use elec_static_field_module, only: totalsym_field, E_ele, E_nuc,totsym_field_length
  use elec_static_field_module, only: transform_to_cart_field, get_field_nuc
  use efp_module, only: print_id, do_pol, energy_conv, n_density_updates
  use efp_efp_module, only: efp_efp_energy,efp_efp_en,efp_efp_gradients,qm_efp_energy
  use efp_polar_module, only: calc_polar
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public efp_opt_main, geom_converged
  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  real(r8_kind) :: qm_qm_energy,total_energy,total_energy_save

  integer(i4_kind), parameter :: snd_epe=499_i4_kind
  integer(i4_kind), parameter :: get_en=498_i4_kind
  integer(i4_kind), parameter :: get_fld=497_i4_kind
  integer(i4_kind), parameter :: get_grd=496_i4_kind

  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine efp_opt_main()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only:  comm_i_am_master,comm_parallel,comm_myindex
    use operations_module, only: operations_gradients,operations_geo_opt
    use pointcharge_module, only: print_energy
    use iounitadmin_module, only: write_to_trace_unit,no_trace_output, &
         no_output_unit_output,output_unit
    use convergence_module, only: convergence_max_geo_iteration
    use efp_module, only: efp_max_grad,efp_write_gxfile,calc_efield_points1
    use energy_calc_module, only: get_energy
    use hesse_module, only: step_counter
    use elec_static_field_module, only: receive_surf_point1,surf_points_grad_information1
    use induced_dipoles_module, only: send_receive_id1
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i, n_iter
    character(len=80) :: trace_buff
    real(r8_kind) :: max_grad
    logical :: geo_conv
    real(r8_kind) :: previous_t,next_t,sum_t
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       call write_to_trace_unit("======================================================")
       call write_to_trace_unit("OPTIMIZATION OF EFP REGION - FIXED QM DENSITY APPROACH")
       call write_to_trace_unit("------------------------------------------------------")
       print*,"======================================================"
       print*,"OPTIMIZATION OF EFP REGION - FIXED QM DENSITY APPROACH"
       print*,"------------------------------------------------------"
       if( .not. operations_geo_opt) then
          call write_to_trace_unit("Just Energy and Gradients Calculation")
          call write_to_trace_unit("------------------------------------------------------")
          print*,"Just Energy and Gradients Calculation"
          print*,"------------------------------------------------------"
       else
          call write_to_trace_unit("Iter  Energy            Max grad    Time")
          call write_to_trace_unit("------------------------------------------")
       end if
    end if

    if(comm_i_am_master()) then
       call get_energy(tot=total_energy)
       total_energy_save=total_energy
    end if

    if(comm_i_am_master()) then
       n_iter=convergence_max_geo_iteration()
       step_counter=0
       previous_t=0.0_r8_kind
       call time_period(previous_t,next_t,sum_t)
    end if

    if(comm_parallel()) call start_cycle(n_iter)

    i=0
    geo_run:do
       !I   - ALLOCATE GRADIENTS
       call alloc_EFP_grads()

       i=i+1
       if(comm_i_am_master()) then
          if(i > 1) step_counter=step_counter+1
          print*,''
          print*,'======================================'
          print*,'EFP GEOMETRY ITERATION: ',i
          print*,'======================================'
          print*,''

          if(i > 1) call read_efp_coor()
       end if
       if(i > 1) then
          if(comm_parallel()) call send_receive_efp()
          if(do_pol) then
             if(comm_i_am_master()) then
                call calc_efield_points1()
             else
                call receive_surf_point1
             end if
             call surf_points_grad_information1()
          end if
       end if
       if(comm_i_am_master()) no_output_unit_output=.false.

       !II  - ENERGY CALCULATION
       no_trace_output=.true.
       !1. Calculation of QM-EFP energy
       call calc_QM_EFP_energy()

       !2. Calculation of electrostatic field at the current positions
       !   of EFP centers
       if(do_pol) call calc_QM_field()

       !3. Calculation EFP-EFP and total energy
       if(comm_i_am_master()) then
          calc_polar=.true.
          call efp_efp_energy(print_id)
          calc_polar=.false.
          qm_efp_energy=qm_efp_energy+efp_efp_en
          if(i==1) then
             qm_qm_energy=total_energy-qm_efp_energy
          end if
          total_energy=qm_qm_energy+qm_efp_energy
          if(print_energy) print*,'TOTAL EFP energy',qm_efp_energy
          if(print_energy) print*,'TOTAL energy',total_energy
       end if

       if(comm_parallel() .and. do_pol) call send_receive_id1()

      !III - GRADIENTS  CALCULATION
       !4. QM-EFP gradients
       call calc_QM_EFP_grads()

       !5. EFP-EFP gradients and total gradients
       if(comm_i_am_master()) then
          call efp_efp_gradients()

          write(output_unit,'(a28)') 'Independent EFP optimization'
          write(output_unit,'(a18,f19.10)') 'TOTAL EFP energy: ' ,qm_efp_energy
          write(output_unit,'(a18,f19.10)') 'TOTAL energy:     ',total_energy
          write(output_unit,*) ''
          call transform_EFP_grad_to_cart()
       end if

       no_trace_output=.false.
       if(.not. operations_geo_opt) then
          if(comm_i_am_master()) then
             max_grad=efp_max_grad()
             write(trace_buff,'(a15,e16.8,a12,f10.7)') &
                  "QM-EFP Energy: ",qm_efp_energy, "  Max Grad: ",max_grad
             call write_to_trace_unit(trim(trace_buff))
          end if
          call dealloc_EFP_grads()
          if(comm_parallel()) call exit_cycle_grad()
          exit geo_run
       else
          if(comm_i_am_master()) then
             call time_period(previous_t,next_t,sum_t)
             max_grad=efp_max_grad()
             write(trace_buff,'(i4,2x,e16.8,2x,f10.7,1x,f7.2)') &
                  i,qm_efp_energy,max_grad,next_t
             call write_to_trace_unit(trim(trace_buff))

             call efp_write_gxfile(qm_efp_energy)

             call optimization_step(geo_conv)
          end if
          call dealloc_EFP_grads()
          if(comm_parallel()) call exit_cycle_opt(geo_conv,i)
          if(geo_conv) exit geo_run
          if(i == n_iter) exit geo_run
       end if
       if(comm_parallel()) call check_cycle(i)
    end do geo_run

    if(comm_i_am_master()) then
       call time_period(previous_t,next_t,sum_t)
       call write_to_trace_unit("------------------------------------------------------")
       call write_to_trace_unit("OPTIMIZATION OF EFP REGION - FINISHING")
       call write_to_trace_unit("======================================================")
       print*,"------------------------------------------------------"
       print*,"OPTIMIZATION OF EFP REGION - FINISHING"
       print*,"======================================================"
    end if

    call main_shutdown()

    if(comm_parallel()) call check_exit()

  end subroutine efp_opt_main
  !*************************************************************

  !*************************************************************
  subroutine start_cycle(n_iterations)
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master,comm_get_n_processors, &
         comm_save_recv,comm_init_send,comm_send,comm_all_other_hosts, &
         comm_master_host, commpack, communpack
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(inout) :: n_iterations
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind),parameter :: i_am_ready=500
    integer(i4_kind),parameter :: begin_cycle=501
    integer(i4_kind) :: n_procs,i,info
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       n_procs=comm_get_n_processors()

       do i=2,n_procs
          call comm_save_recv(i,i_am_ready)
       end do

       call comm_init_send(comm_all_other_hosts,begin_cycle)

       call commpack(n_iterations,info)
       ASSERT(info==0)

       call commpack(do_pol,info)
       ASSERT(info==0)

       call comm_send()
    else
       call comm_init_send(comm_master_host,i_am_ready)
       call comm_send()

       call comm_save_recv(comm_master_host,begin_cycle)

       call communpack(n_iterations,info)
       ASSERT(info==0)

       call communpack(do_pol,info)
       ASSERT(info==0)
    end if

  end subroutine  start_cycle
  !*************************************************************

  !*************************************************************
  subroutine exit_cycle_grad()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master,comm_get_n_processors, &
         comm_save_recv,comm_init_send,comm_send,comm_all_other_hosts, &
         comm_master_host
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind),parameter :: i_am_ready=510
    integer(i4_kind),parameter :: leave_cycle=511
    integer(i4_kind) :: n_procs,i
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       n_procs=comm_get_n_processors()

       do i=2,n_procs
          call comm_save_recv(i,i_am_ready)
       end do

       call comm_init_send(comm_all_other_hosts,leave_cycle)
       call comm_send()
    else
       call comm_init_send(comm_master_host,i_am_ready)
       call comm_send()

       call comm_save_recv(comm_master_host,leave_cycle)
    end if

  end subroutine exit_cycle_grad
  !*************************************************************

  !*************************************************************
  subroutine exit_cycle_opt(conv,iter)
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master,comm_get_n_processors, &
         comm_save_recv,comm_init_send,comm_send,comm_all_other_hosts, &
         comm_master_host, commpack, communpack
    !------------ Declaration of formal parameters ---------------
    logical,intent(inout) :: conv
    integer(i4_kind),intent(inout) :: iter
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind),parameter :: i_am_ready=520
    integer(i4_kind),parameter :: leave_cycle=521
    integer(i4_kind) :: n_procs,i,info
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       n_procs=comm_get_n_processors()

       do i=2,n_procs
          call comm_save_recv(i,i_am_ready)
       end do

       call comm_init_send(comm_all_other_hosts,leave_cycle)

       call commpack(conv,info)
       ASSERT(info==0)

       call commpack(iter,info)
       ASSERT(info==0)

       call comm_send()
    else
       call comm_init_send(comm_master_host,i_am_ready)
       call comm_send()

       call comm_save_recv(comm_master_host,leave_cycle)

       call communpack(conv,info)
       ASSERT(info==0)

       call communpack(iter,info)
       ASSERT(info==0)
    end if

  end subroutine exit_cycle_opt
  !*************************************************************

  !*************************************************************
  subroutine check_cycle(iter)
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master,comm_get_n_processors, &
         comm_save_recv,comm_init_send,comm_send,comm_all_other_hosts, &
         comm_master_host, commpack, communpack
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent(in) :: iter
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind),parameter :: i_am_ready=530
    integer(i4_kind),parameter :: next_cycle=531
    integer(i4_kind) :: n_procs,i,info,iter_slave
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       n_procs=comm_get_n_processors()

       do i=2,n_procs
          call comm_save_recv(i,i_am_ready)

          call communpack(iter_slave,info)
          ASSERT(info==0)

          ASSERT(iter==iter_slave)
       end do

       call comm_init_send(comm_all_other_hosts,next_cycle)
       call comm_send()
    else
       call comm_init_send(comm_master_host,i_am_ready)

       call commpack(iter,info)
       ASSERT(info==0)

       call comm_send()

       call comm_save_recv(comm_master_host,next_cycle)
    end if

  end subroutine check_cycle
  !*************************************************************

  !*************************************************************
  subroutine check_exit()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master,comm_get_n_processors, &
         comm_save_recv,comm_init_send,comm_send,comm_all_other_hosts, &
         comm_master_host, commpack, communpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind),parameter :: i_am_ready=540
    integer(i4_kind),parameter :: stop_opt=541
    integer(i4_kind) :: n_procs, i
    !------------ Executable code --------------------------------

    if(comm_i_am_master()) then
       n_procs=comm_get_n_processors()

       do i=2,n_procs
          call comm_save_recv(i,i_am_ready)
       end do

       call comm_init_send(comm_all_other_hosts,stop_opt)
       call comm_send()
    else
       call comm_init_send(comm_master_host,i_am_ready)
       call comm_send()

       call comm_save_recv(comm_master_host,stop_opt)
    end if

  end subroutine  check_exit
  !*************************************************************

  !*************************************************************
  subroutine read_efp_coor()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use efp_module, only: read_gx_qm,def_efp_arrays,calc_X_points, &
         calc_efield_points
    use induced_dipoles_module, only: calc_Pol_centers
    use iounitadmin_module, only: no_output_unit_output
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    call read_gx_qm()
    call def_efp_arrays()
    no_output_unit_output=.true.
    call calc_X_points()
!!$    if(calc_Pol_centers()) call calc_efield_points()

  end subroutine read_efp_coor
  !*************************************************************

  !*************************************************************
  subroutine send_receive_efp()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use pointcharge_module,     only: pointcharge_bcast
    use point_dqo_module,       only: external_centers_bcast
    use induced_dipoles_module, only: pol_centers_bcast
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------
    !
    call pointcharge_bcast()
    !
    call external_centers_bcast()
    !
    call pol_centers_bcast()
    !
  end subroutine send_receive_efp
  !*************************************************************

  !*************************************************************
  subroutine calc_QM_EFP_energy()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use integralpar_module, only: integralpar_set
    use interfaces, only: main_integral
    use comm_module, only:  comm_i_am_master, comm_parallel
    use pointcharge_module, only : calc_nuc_pc_energy, print_energy
    use point_dqo_module, only : calc_nuc_dqo_energy
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    qm_efp_energy=0.0_r8_kind

    call integralpar_set("QM_EFP_energy")

    call main_integral ()

    if(comm_parallel()) then
       if(comm_i_am_master()) then
          call receive_energy()
       else
          call send_energy()
       end if
    end if

    if(comm_i_am_master()) then
       if(print_energy) print*,'ELECTRON-EFP energy',qm_efp_energy
       qm_efp_energy=qm_efp_energy+calc_nuc_pc_energy()
       qm_efp_energy=qm_efp_energy+calc_nuc_dqo_energy()
       if(print_energy) print*,'QM-EFP energy',qm_efp_energy
    end if

  end subroutine calc_QM_EFP_energy
  !*************************************************************

  !*************************************************************
  subroutine receive_energy()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_get_n_processors, comm_save_recv, communpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: info,n_procs,i
    real(r8_kind) :: e_buf
    !------------ Executable code --------------------------------

    n_procs=comm_get_n_processors()

    do i=2,n_procs
       call comm_save_recv(i,get_en)

       call communpack(e_buf,info)
       ASSERT(info==0)

       qm_efp_energy=qm_efp_energy+e_buf
    end do

  end subroutine receive_energy
  !*************************************************************

  !*************************************************************
  subroutine send_energy()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_init_send,comm_master_host,comm_send, commpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: info
    !------------ Executable code --------------------------------

    call comm_init_send(comm_master_host,get_en)

    call commpack(qm_efp_energy,info)
    ASSERT(info==0)

    call comm_send()

  end subroutine send_energy
  !*************************************************************

  !*************************************************************
  subroutine calc_QM_field()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use integralpar_module, only: integralpar_set
    use interfaces, only: main_integral
    use comm_module, only:  comm_i_am_master, comm_parallel
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    call integralpar_set("Field_at_EFP")

    call main_integral ()

    if(comm_parallel()) then
       if(comm_i_am_master()) then
          call receive_field()
       else
          call send_field()
       end if
    end if

    if(comm_i_am_master()) then
       call transform_to_cart_field(E_ele,totalsym_field)
       totalsym_field=0.0_r8_kind

       call get_field_nuc()
    end if


  end subroutine calc_QM_field
  !*************************************************************

  !*************************************************************
  subroutine receive_field()
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_get_n_processors,comm_save_recv, communpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: info,n_procs,i
    real(r8_kind) :: help_arr(totsym_field_length)
    !------------ Executable code --------------------------------

    n_procs=comm_get_n_processors()

    do i=2,n_procs
       call comm_save_recv(i,get_fld)

       call communpack(help_arr,totsym_field_length,1,info)
       ASSERT(info==0)

       totalsym_field=totalsym_field+help_arr
    end do

  end subroutine receive_field
  !*************************************************************

  !*************************************************************
  subroutine send_field()
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_init_send,comm_master_host,comm_send, commpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: info
    !------------ Executable code --------------------------------

    call comm_init_send(comm_master_host,get_fld)

    call commpack(totalsym_field,totsym_field_length,1,info)
    ASSERT(info==0)

    call comm_send()

  end subroutine send_field
  !*************************************************************

  !*************************************************************
  subroutine alloc_EFP_grads()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use pointcharge_module, only: moving_pc, init_pointcharges_grads
    use point_dqo_module, only: moving_X_centers, moving_R_centers, init_X_centers_grads
    use induced_dipoles_module, only: moving_Pol_centers, init_Pol_centers_grads
    use efp_rep_module, only: init_repf_grads
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    if (moving_pc) call init_pointcharges_grads()
    if (moving_X_centers .or. moving_R_centers) call init_X_centers_grads()
    if (moving_R_centers) call init_repf_grads()
    if (moving_Pol_centers) call init_Pol_centers_grads()

  end subroutine alloc_EFP_grads
  !*************************************************************

  !*************************************************************
  subroutine dealloc_EFP_grads()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use pointcharge_module, only: moving_pc, pc_grads_shutdown
    use point_dqo_module, only: moving_X_centers, moving_R_centers, X_centers_grads_shutdown
    use induced_dipoles_module, only: moving_Pol_centers, Pol_centers_grads_shutdown
    use efp_rep_module, only: repf_grads_shutdown
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    if (moving_pc) call pc_grads_shutdown()
    if (moving_X_centers .or. moving_R_centers) call X_centers_grads_shutdown()
    if (moving_R_centers) call repf_grads_shutdown()
    if (moving_Pol_centers) call Pol_centers_grads_shutdown()

  end subroutine dealloc_EFP_grads
  !*************************************************************

  !*************************************************************
  subroutine calc_QM_EFP_grads()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use integralpar_module, only: integralpar_set
    use interfaces, only: main_integral
    use comm_module, only:  comm_i_am_master, comm_parallel
    use pointcharge_module, only: moving_pc, calc_PC_grads
    use point_dqo_module, only: moving_X_centers, calc_X_grads
    use induced_dipoles_module, only: moving_Pol_centers
    use calc_id_module, only: calc_id_grads
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    call integralpar_set("EFP_gradients")

    call main_integral ()

    if(comm_parallel()) then
       if(comm_i_am_master()) then
          call receive_efp_grads()
       else
          call send_efp_grads()
       end if
    end if

    if(comm_i_am_master()) then
       if(moving_pc) call calc_PC_grads()
       if(moving_X_centers) call calc_X_grads()
       if(moving_Pol_centers) call calc_id_grads()
    end if

  end subroutine calc_QM_EFP_grads
  !*************************************************************

  !*************************************************************
  subroutine receive_efp_grads()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_get_n_processors,comm_save_recv
    use pointcharge_module, only: moving_pc,totsym_PC_grad_unpack
    use point_dqo_module, only: moving_X_centers,moving_R_centers, &
         totsym_X_grad_unpack
    use induced_dipoles_module, only: moving_Pol_centers,totsym_id_grad_unpack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: n_procs,i
    !------------ Executable code --------------------------------

    n_procs=comm_get_n_processors()

    do i=2,n_procs
       call comm_save_recv(i,get_grd)

       if(moving_pc) call totsym_PC_grad_unpack()
       if(moving_X_centers .or. moving_R_centers) call totsym_X_grad_unpack()
       if(moving_Pol_centers) call totsym_id_grad_unpack()

    end do

  end subroutine receive_efp_grads
  !*************************************************************

  !*************************************************************
  subroutine send_efp_grads()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_init_send,comm_master_host,comm_send
    use pointcharge_module, only: moving_pc,totsym_PC_grad_pack
    use point_dqo_module, only: moving_X_centers,moving_R_centers, &
         totsym_X_grad_pack
    use induced_dipoles_module, only: moving_Pol_centers,totsym_id_grad_pack
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    call comm_init_send(comm_master_host,get_grd)

    if(moving_pc) call totsym_PC_grad_pack()
    if(moving_X_centers .or. moving_R_centers) call totsym_X_grad_pack()
    if(moving_Pol_centers) call totsym_id_grad_pack()

    call comm_send()

  end subroutine send_efp_grads
  !*************************************************************

  !*************************************************************
  subroutine transform_EFP_grad_to_cart()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use pointcharge_module, only: transform_PC_grad_to_cart, &
         pc_grad_cart_write,print_pc_grad, moving_pc
    use point_dqo_module, only: transform_X_grad_to_cart, &
         X_grad_cart_write, print_X_grad, print_R_grad, &
         moving_X_centers, moving_R_centers
    use induced_dipoles_module, only: moving_Pol_centers, &
         print_id_grad
    use calc_id_module,only: transform_id_grad_to_cart, &
         id_grad_cart_write
    use efp_rep_module, only: transform_repf_grad_to_cart, &
         repf_grad_cart_write
    use efp_module, only: efp_sum_up_grads, efp_grad_cart_write
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    if(moving_pc) then
       call transform_PC_grad_to_cart()
       if(print_pc_grad) call pc_grad_cart_write()
    end if
    if(moving_X_centers .or. moving_R_centers) then
       call transform_X_grad_to_cart()
       if(print_X_grad .or. print_R_grad) call X_grad_cart_write()
    end if
    if(moving_R_centers) then
       call transform_repf_grad_to_cart()
       if(print_R_grad) call repf_grad_cart_write()
    end if
    if(moving_Pol_centers) then
       call transform_id_grad_to_cart()
       if(print_id_grad) call id_grad_cart_write()
    end if

    call efp_sum_up_grads()
    call efp_grad_cart_write()

  end subroutine transform_EFP_grad_to_cart
  !*************************************************************

  !*************************************************************
  subroutine optimization_step(geo_conv)
    !  Purpose:
    !------------ Modules used -----------------------------------
    use optimizer, only: main_opt
    !------------ Declaration of formal parameters ---------------
    logical         , intent(out) :: geo_conv
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    logical :: conv
    character(len=32) :: optimizer_task='GeoOpt'
    !------------ Executable code --------------------------------

    call main_opt( task=optimizer_task,converged=conv )

    geo_conv = conv

  end subroutine optimization_step
  !*************************************************************

  !*************************************************************
  function geom_converged(qm_cycle) result(conv)
    !  Purpose:
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    logical :: conv
    integer(i4_kind) :: qm_cycle
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    conv=.false.

print*,total_energy_save,total_energy,abs(total_energy-total_energy_save),'@@@@@@@@@@@@@@'
    if(abs(total_energy-total_energy_save) <= energy_conv) conv=.true.
    if(qm_cycle >= n_density_updates) conv=.true.

  end function geom_converged
  !*************************************************************

  !*************************************************************
  subroutine main_shutdown()
    !  Purpose:
    !------------ Modules used -----------------------------------
    use interfaces, only: IPARA
    use unique_atom_methods, only: unique_atom_close
    use pointcharge_module, only: pointcharge_close
    use point_dqo_module, only: dealloc_dqo,moving_X_centers,moving_R_centers, &
         dealloc_center_inform
    use induced_dipoles_module, only: dealloc_id,dealloc_Pol_center_inform
    use efp_rep_module, only: dealloc_repf,dealloc_repf_inform
    use occupied_levels_module, only: eigvec_occ_dealloc
    use efp_module, only: shutdown_efp_arrays
    use comm_module, only:  comm_i_am_master
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    call unique_atom_close()

    call pointcharge_close()

    call dealloc_dqo()

    if(moving_X_centers .or. moving_R_centers) call dealloc_center_inform()

    call dealloc_id()

    call dealloc_Pol_center_inform()

    if(comm_i_am_master())then
       call dealloc_repf()
       if(moving_R_centers) call dealloc_repf_inform()

       call shutdown_efp_arrays()
    end if

    call eigvec_occ_dealloc()
  end subroutine main_shutdown
  !*************************************************************

  !*************************************************************
  subroutine time_period(previous,next,sum)
    !  Purpose:
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(inout) ::previous,next,sum
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: previous1
    character(len=10) :: date,time,zone
    integer(i4_kind) :: values(8)
    !------------ Executable code --------------------------------

    previous1=previous
    call date_and_time(date,time,zone,values)
    next=values(8)*0.001_r8_kind+values(7)*1.0_r8_kind+ &
         values(6)*60.0_r8_kind+values(5)*3600.0_r8_kind+ &
         values(3)*86400.0_r8_kind
    if(previous == 0.0_r8_kind) then
       previous=next
       next=0.0_r8_kind
       sum=0.0_r8_kind
    else
       previous=next
       next=next-previous1
    endif
    sum=sum+next

  end subroutine time_period
  !*************************************************************

  !*************************************************************
end module efp_only_opt_module

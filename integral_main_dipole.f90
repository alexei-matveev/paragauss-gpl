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
subroutine  integral_main_dipole
  !-------------------------------------------------------------------
  !
  !  Purpose: This is the main routine of the dipole integral part
  !           executed by the master.
  !
  !  Subroutine called by: main_integral
  !
  !  References: Publisher Document: Concepts of Integral Part
  !
  !
  !  Author: TB
  !  Date: 7/97
  !
  !
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   6/97
  ! Description: ...
  !
  ! Modification master/slave -> DLB
  ! Author: AN
  ! Date:   4/11
  ! Description: use DLB for scheduling of dipole integrals
  !              master/slave use the same code to run
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
# include "def.h"
  !------------ Modules used -----------------------------------------
  use type_module          ! type specification parameters
  use output_module        ! defines amount of output
  use iounitadmin_module   ! to open output units
  use comm_module           ! comm related information and routines
  use msgtag_module
  use integralpar_module   ! steering information for integral part
  use int_distribute_module
  use int_data_dipole_module, only: quadrupel
  use int_send_dipole_module
  use timer_module
  use time_module

  implicit none

  integer(i4_kind) :: n_quad

  !-------------------------------------------------------------------
  !------------ Executable code -----------------------------------

! ================ CONTEXT: EVERY PROC ================

  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_main_dipole: start")

  if( comm_i_am_master() )then
    ! ================ CONTEXT: MASTER ONLY ================
    call start_timer(timer_int_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))


    ! prepare distribution of integral quadrupels
    if (output_int_detailedprogress) call write_to_output_units( &
         "integral_main_dipole: int_distribute_setup")
  endif ! i am master
  ! ================ CONTEXT: EVERY PROC ================
    n_quad = int_distribute_setup(int_distribute_dipole_run)

  ! do setup work, starting the idle timer
  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_main_dipole: setup")
  call integral_setup_dipole()

  if( comm_i_am_master() )then
    ! ================ CONTEXT: MASTER ONLY ================
    ! prepare receiving of results
    if (output_int_detailedprogress) call write_to_output_units( &
         "integral_main_dipole: int_send_dipole_setup")
    call int_send_dipole_setup(n_quad)
    ! start calculation
    if (output_int_detailedprogress) call write_to_output_units( &
         "integral_main_dipole: int_distribute_start")
  endif
  ! ================ CONTEXT: EVERY PROC ================
    ! last preparations for actual start, for example DLB setup
    call int_distribute_start()

    ! The main loop for calculation, int_distribute_next_job tells
    ! if there is still someting to do for the current processor, and if
    ! yes on which quadrupel to work on, integral_calc_quad_dipole works
    ! on quadrupel stored in int_data_dipole_module
    do while ( int_distribute_next_job(quadrupel))
         call integral_calc_quad_dipole()
        if (comm_i_am_master()) then
            ! dipole results will be sended to master
            ! which will non-blocking check for arrival
            call try_get_results()
        endif
    enddo

  ! do shutdown work
  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_main_dipole: shutdown")

    if (comm_i_am_master() .and. output_int_detailedprogress) &
          call write_to_output_units( &
         "integral_main_dipole: int_distribute_shutdown")
    call int_distribute_shutdown()

  ! ================ CONTEXT: MASTER ONLY ================
  if( comm_i_am_master() )then
    ! shutdown int_send_module
    if (output_int_detailedprogress) call write_to_output_units( &
         "integral_main_dipole: int_send_dipole_shutdown")
    ! Master will receive all missing dipoles
    call int_send_dipole_shutdown()
  endif ! i am master

  ! stop timer
  call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
  call stop_timer(timer_int_2cob3c(integralpar_i_int_part))

  ! ================ CONTEXT: EVERY PROC ================
  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_main_dipole: done")

contains

  subroutine try_get_results()
    !Purpose: master tries to get the results from the slaves
    !         uses nonblocking mpi
    use comm_module, only: comm_save_recv_nonblocking
    implicit none
    do while (comm_save_recv_nonblocking &
                (comm_all_other_hosts, msgtag_int_dipole_result))
             call int_send_dipole_receive()
    enddo
  end subroutine try_get_results

end subroutine integral_main_dipole

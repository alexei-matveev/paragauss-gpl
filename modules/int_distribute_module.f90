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
module  int_distribute_module
!---------------------------------------------------------------
!
!  Purpose: This module does the distribution of quadrupels
!           ( unique atom 1, l 1, unique atom 2, l 2 )
!           on hosts for dynamic load balancing.
!
!           It can be used both for 2 center fitfunction
!           integral calculation and for 2 center orbital
!           and three center integral calculation.
!
!           In case norm integrals are calculated, it
!           takes into account that the integrals for
!           diagonal quadrupels (ua1,l1,ua1,l1) have to be
!           calculated before any other integrals of (ua1,l1)
!           can be calculated.
!
!           The module tries to start calculations with higher
!           priority (that take longer) first to minimize waiting
!           time at end. However, in case norms are calculated,
!           diagonal quadrupels with lowest priority are dispatched
!           first to obtain non-diagonal quadrupels that can be 
!           calculated as soon as possible.
!
!           The module dynamically decides if master should
!           join in calculating integrals. It takes its initial
!           strategy from the setting in comm-module but then
!           watches if the sum of idle_times of the slaves is higher
!           than the time the master spends upon calculating
!           quadrupels.
!
!
!           Call int_distribute_setup(is_2cff) before using module
!           and int_distribute_shutdown() afterwards.
!           Call int_distribute_start() to start calculation.
!           Call int_distribute_newtask(host,idle_time) when host
!              finished work on quadruple.
!           Call int_distribute_workleft() to inquire if there
!              are any tasks left to be done.
!           Call int_distribute_mastertask() to inquire if there
!              is a new task for master
!
!
!  Module called by: integral_main_2zff, integral_interrupt_2zff,
!                    integral_main_2zob3z, integral_interrupt_2zob3z
!  (executed only by master)
!
!  References: Publisher Document: Concepts of Integral Part
! 
!
!  Author: TB
!  Date: 5/96
!
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing) 
! Author:      Uwe Birkenheuer
! Date:        8/98
! Description: distribution of quadrupel tasks adapted to the 
!              new control parameter INTEGRALPAR_2C_MIXED
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification master/slave -> DLB
! Author: AN
! Date:   4/11
! Description: the complex algorithm with master/slave concept
!              with variants like master_work and so on is replaced by
!              simple call to DLB library which will distribute and
!              redistribute the tasks accordingly. Therefore many of the
!              old routines are not needed any more and omitted.
!              int_distribute_setup, int_distribute_shutdown and
!              int_distribute_start are still for the same tasks, only with
!              different content, note that int_distribute_start does not
!              start the actual calculation, this would be only done by
!              a following call of int_distribute_next_job
!              the function int_distribute_next_job(quadrupels) is the new function
!              it might be used to get a valid quadrupel to work on, consider
!              that int_distribute_next_job has two outputs: the return value is
!              an logical telling if there are any jobs left, and qudrupels is
!              (but only if return value was true) a valid set of not yet
!              calculated quadrupel, which is supposed to be calculated the next
!              time the processor re-asks for a new job
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  use quadrupel_module
  use quadrupelstore_module
  use integerstack_module
  use integralpar_module
  use unique_atom_module
  use comm_module
  use comm, only: comm_size
  use msgtag_module
  use time_module
  use timer_module
  use output_module
  use iounitadmin_module

  implicit none

  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module =================

  !------------ Declaration of public constants ----------
  integer, parameter, public :: &
       int_distribute_2cob3c_run = 1, &
       int_distribute_2cff_run = 2, &
       int_distribute_dipole_run = 3

  !------------ public functions and subroutines ----------------
  public int_distribute_setup, int_distribute_shutdown, &
         int_distribute_start

public int_distribute_next_job
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ public entitities from used modules -------------
  public quadrupel_type


  !------------ Declaration of private types --------------------

  !------------ Declaration of constants and variables ----------
  type(quadrupel_type), allocatable :: quadrupels(:)
  ! quadrupels:               quadrupel to calculate (in case
  !                          norms are calculated except for diagonal)
  integer(i4_kind) :: n_procs ! for not having to calculate it over and
                              ! over again
  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains


  !**************************************************************
  integer function int_distribute_setup(run_type)
    !  Purpose: setup routine that has to be called before module
    !           can be used
    implicit none
    !------------ Declaration of formal parameters --------------
    integer, intent(in) :: run_type
    ! legal values are: int_distribute_2cob3c_run, 
    !    int_distribute_2cff_run, int_distribute_dipole_run
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: i_ua
    integer(kind=i4_kind) :: min_l, max_l(N_unique_atoms)
    integer(kind=i4_kind) :: which_2cff
    logical :: need_ch, need_xc
    logical :: fitfunction, quadrupels_triangular
    ! fitfunction: true if fitfunction integrals are to be calculated
    ! quadrupels_triangular: only quadrupels of upper ua l triangle required
    !------------ Executable code -------------------------------

    n_procs = comm_size()
!..............................................................................
! If only charge fit integrals have to be evaluated the 2-center integral
! quadrupels are distributed for the lower triangle of the corresponding
! n_ch x n_ch matrix (which_2cff = 1). 
!
! If, on the other hand, only echange fit integrals have to be evaluated 
! the quadrupels are distributed for a n_xc x n_xc matrix (which_2cff = 2)
!
! If, however, both types 2-center integrals and/or mixed overlap integrals
! are to be evaluated the 2-center integral quadrupels are distributed such
! as to evaluate a quadratic matrix of order max(n_ch,n_xc).
!..............................................................................

    select case (run_type)
    case(int_distribute_2cob3c_run)
       fitfunction = .false.
       quadrupels_triangular = .true.

    case(int_distribute_2cff_run)
       fitfunction = .true.
       quadrupels_triangular = .true.
    case(int_distribute_dipole_run)
       fitfunction = .false.
       ! FIXME: why so?
       quadrupels_triangular = .not. integralpar_offdiag_dipoles
    case default
       call error_handler("int_distribute_setup: illegal value for run_type")
    end select

    if (integralpar_gradients) then
       need_ch = .true.
       need_xc = .false.
    else
       need_ch = integralpar_2cch_no .or. integralpar_2c_mixed
       need_xc = integralpar_2cxc_no .or. integralpar_2c_mixed
    endif
    if (need_ch) then
       if (need_xc) then
          which_2cff = 3
       else
          which_2cff = 1
       endif
    else
       if (need_xc) then
          which_2cff = 2
       else
          which_2cff = 0
       endif
    endif


    if ( fitfunction ) then
       min_l = -1
    else
       min_l = 0
    endif

    do i_ua = 1, N_unique_atoms
        if ( fitfunction ) then
           select case (which_2cff)
           case (1)
              max_l(i_ua) = unique_atoms(i_ua)%lmax_ch
           case (2)
              max_l(i_ua) = unique_atoms(i_ua)%lmax_xc
           case (3)
              max_l(i_ua) = max( unique_atoms(i_ua)%lmax_ch, &
                           unique_atoms(i_ua)%lmax_xc )
           case default
              call error_handler( &
                   "int_distribute_setup: invalid fitfunction case")
           end select
        else
           max_l(i_ua) = unique_atoms(i_ua)%lmax_ob
        endif
    enddo


    call quadrupelstore_construct(quadrupels, min_l, max_l, quadrupels_triangular, fitfunction)
    int_distribute_setup = size(quadrupels)

    ! master only print output and have timer:
    if (comm_i_am_master()) then
        call write_to_trace_unit("int_disribute_setup: &
            &total number of quadrupels: ",inte= size(quadrupels))
    endif

  end function int_distribute_setup
  !**************************************************************


  !**************************************************************
  subroutine int_distribute_shutdown()
    !  Purpose: shutdown routine that has to be called after
    !           module has been used
    !** End of interface ****************************************
    implicit none
    integer(kind=i4_kind) :: ierr
    !------------ Executable code -------------------------------

    deallocate(quadrupels, stat = ierr)
    ASSERT(ierr == 0)
  end subroutine int_distribute_shutdown
  !**************************************************************


  !**************************************************************
  subroutine int_distribute_start()
    !  Purpose: This subroutine distributes tasks statically and
    !           thus makes last preparations for run
    !** End of interface ****************************************
    use dlb, only: dlb_setup_rr, idlb_kind
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=idlb_kind) :: N
    !------------ Executable code -------------------------------
        N = size(quadrupels)

         call dlb_setup_rr(N)
  end subroutine int_distribute_start
  !**************************************************************



  !**************************************************************
  logical function int_distribute_next_job(quad)
    !  Purpose: returns if any task s available.
    !           In case there is, quadrupel should be calculated.
    !           Else quadrupel is undefined.
    use dlb, only: dlb_give_more_rr, idlb_kind
    implicit none
    !------------ Declaration of formal parameters --------------
    type(quadrupel_type), intent(out)  :: quad
    !** End of interface ****************************************
    integer(idlb_kind) :: jobs(3), length
    integer(i4_kind) :: q_num
    !------------ Executable code -------------------------------
    length = 1
    int_distribute_next_job = dlb_give_more_rr(length, jobs)

    if (int_distribute_next_job) then
        ! transform them in guaranteed i4_kind
        q_num = jobs(1) + 1
        quad = quadrupelstore_give_num(quadrupels, q_num)
    endif

  end function int_distribute_next_job

  !**************************************************************

end module int_distribute_module

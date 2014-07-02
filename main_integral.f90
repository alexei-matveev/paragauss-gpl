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
subroutine main_integral ()
!----------------------------------------------------------------
!
!  Purpose: This is the main routine of the integral part
!           executed by the master.
!
!  Subroutine called by: main_mater
!
!  References: Publisher Document: Concepts of Integral Part
! 
!
!  Author: TB
!  Date: 5/96
!
!
  !===================================================================
  ! End of public interface of module
  !===================================================================
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
!------------ Modules used --------------------------------------
use type_module          ! type specification parameters
use output_module        ! defines amount of output
use iounitadmin_module   ! to open output units
use comm, only: comm_barrier
use comm_module           ! comm related information and routines
use msgtag_module
use integralpar_module   ! steering information for integral part
use fit_coeff_module, only: fit_coeff_calc_chargefitnorm !, fit_coeff_read_chargefitnorm
use time_module, only: clktime
use timer_module
USE_MEMLOG
use xpack, only: pck
use operations_module, only: operations_gradients
! possible values of argument CONTEXT are combinations of:
use interfaces, only: integral_trafo, RELENRG, RELGRAD
use shgi_utils, only: shgi_timing
implicit none
! *** end of interface ***


!----------------------------------------------------------------
logical :: io
!integer(kind=i4_kind):: i_grad
integer(kind=i4_kind):: iop
real(r8_kind)        :: start_time, stop_time, comp_time
!------------ Executable code -----------------------------------

! ================ CONTEXT: EVERY PROC ================

start_time = clktime()

! only master does IO:
io = comm_i_am_master()

if (integralpar_cpks_contribs) then
   call say("integralpar_cpks_contribs to be calculated")
else
   call say("no integralpar_cpks_contribs")
endif

! ================ CONTEXT: EVERY PROC ================
if (output_int_progress) call write_to_output_units( &
     "main_integral:  start")
call say("start")

if (output_int_progress) call write_to_output_units( &
     "main_integral: setup")
call integral_setup()

! calculate 2 center fitfunction integrals
if ( integralpar_2cff ) then
   if (output_int_progress) call write_to_output_units( &
        "main_integral: 2 center fitfunction integrals")
   call say("2 center fitfunction integrals")
   call integral_main_2cff()
endif

if ( comm_i_am_master() ) then
! ================ CONTEXT: MASTER ONLY ================

! calculate 1 center charge fitfunction norm integral
if ( integralpar_1cch_no ) then
   if (output_int_progress) call write_to_output_units( &
        "main_integral: 1 center charge fitfunction norm integral")
   call fit_coeff_calc_chargefitnorm()
endif
endif ! i am master

! ================ CONTEXT: EVERY PROC ================
! calculate 2 center orbital and 3 center integrals
if (integralpar_2cob3c) then
   if (output_int_progress) call write_to_output_units( &
        "main_integral: 2 center orbital and 3 center integrals")
   call say("2 center orbital and 3 center integrals")
   call integral_main_2cob3c()
endif


! calculate dipole integrals
if (integralpar_2cob_dipole) then
   if (output_int_progress) call write_to_output_units( &
        "main_integral: dipole integrals")
   call say("dipole integrals")
   DPRINT 'main_integral: call integral_main_dipole()'
   call integral_main_dipole()
   DPRINT 'main_integral: done integral_main_dipole()'
endif

! do shutdown work
call say("shutdown")

if (output_int_progress) call write_to_output_units("main_integral: shutdown")
call integral_shutdown()

! perform relativistic transformations
!
! make all procs go to integral_trafo() entree
! integral_main_rel() (see above) is called from there
!
if ( integralpar_relativistic.and.(.not.integralpar_rel_gradients) ) then
   ! note: the calculation of relativistic gradients is driven by main_gradient
   call say("relativistic transformations")

   MEMSET(0)
   iop = RELENRG
   if( operations_gradients ) iop = iop + RELGRAD
   ! do PREPARATIONS for gradients that will follow SCF!
   call integral_trafo(iop)
endif

if ( comm_i_am_master() ) then
! ================ CONTEXT: MASTER ONLY ================

! print summary of timings
if (output_timing_integrals .or. output_timing_detailedintegrals) then
   if (output_int_progress) call write_to_output_units( &
        "main_integral: timer_print_integral()")
   call timer_print_integral(integralpar_i_int_part, &
        integralpar_int_part_name(integralpar_i_int_part))
endif
endif ! i am master

! print summary of timings in SHGI integrals:
DCALL shgi_timing(0.0D0) ! from shgi_utils

! ================ CONTEXT: EVERY PROC ================

stop_time = clktime()
comp_time = stop_time - start_time

call say("done")
if (output_int_progress) call write_to_output_units("main_integral: done")

contains

  subroutine say(phrase)
    use iounitadmin_module, only: write_to_trace_unit
    implicit none
    character(len=*), intent(in) :: phrase
    ! *** end of interface ***

    !
    ! Logical io is true on master:
    !
    if ( io ) then
        call write_to_trace_unit("main_integral: "//phrase)
    endif
  end subroutine say

end subroutine main_integral

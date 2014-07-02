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
  !===================================================================
! Public interface of module
  !===================================================================
subroutine  integral_calc_quad_2cff()
  !-------------------------------------------------------------------
  !
  !  Purpose: This routine is the main routine for the
  !           integral calculation of one quadrupel
  !           ( unique atom 1, l 1, unique atom 2, l 2 )
  !           for 2 center fitfunction integral
  !           calculation. The corresponding data (including
  !           the quadrupel to be calculated) are stored in
  !           int_data_2cff_module.
  !
  !  Subroutine called by: main_slave integral_main_2cff
  !
  !  References: Publisher Document: Concepts of Integral Part
  !
  !  Author: TB
  !  Date: 6/96
  !
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   3/97
  ! Description: Necesary chances for gradient run implemented
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: evaluation of mixed charge-fit/exchange-fit integrals
  !              introduced. To this end, new global control parameters
  !              NEED_CH_CH, NEED_CH_XC, NEED_XC_CH, and NEED_XC_XC as
  !              well as a new input flag INT_TYPE to subroutine XC-OVERLAP
  !              is introduced. (Gradients not yet considered)
  !
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------

  !------------ Modules used ------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use timer_module
  use time_module
  use integralpar_module
  use int_data_2cff_module
  use integral_2c_fit_ch_module
  use integral_2c_fit_xc_module
  use int_send_2cff_module
  use output_module, only: output_int_loops, output_int_taskdistribution
  use iounitadmin_module, only: stdout_unit
  use gradient_2c_fit_ch_module
  implicit none
  ! *** end of interface ***

  if ( output_int_loops .or. output_int_taskdistribution) then
     write(stdout_unit,*) "integral_calc_quad_2cff: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
  endif

  call start_timer(timer_int_primcont_2cff(integralpar_i_int_part))

  call say("setup")
  call setup()

  ! calculate primitive integrals, do symmetry adaption and contract
  if(.not.integralpar_gradients) then
     if ( (integralpar_2cch_no .or. diagonal) .and. need_ch_ch ) then
        call say("calling charge_overlap")
        call charge_overlap()
     endif
  else
     call say("calling charge_overlap_grad")
     call charge_overlap_grad()
  end if

  if ( (integralpar_2cxc_no .or. diagonal) .and. need_xc_xc ) then
     call say("calling xc_overlap('xc_xc')")
     call xc_overlap('xc_xc')
  endif

  if ( integralpar_2c_mixed .and. need_ch_xc ) then
     call say("calling xc_overlap('ch_xc')")
     call xc_overlap('ch_xc')
  endif

  if ( integralpar_2c_mixed .and. need_xc_ch ) then
     call say("calling xc_overlap('xc_ch')")
     call xc_overlap('xc_ch')
  endif

  call switch_timer(timer_int_primcont_2cff(integralpar_i_int_part), &
       timer_int_communication_2cff(integralpar_i_int_part))

  call say("report_back")

  if ( integralpar_send_2c) then
     call say("int_send_2cff_send")
     call int_send_2cff_send()
  endif

  call stop_timer(timer_int_communication_2cff(integralpar_i_int_part))

  call say("shutdown")
  call shutdown()

  call say("end")

  !------------ Private Subroutines -----------------------------
contains


  !*************************************************************
  subroutine setup
    !  Purpose: contains setup routines of various modules
    !------------ Executable code ------------------------------

    call stop_timer(timer_int_idle_2cff(integralpar_i_int_part))
    call start_timer(timer_int_calcsum_2cff(integralpar_i_int_part))
    call start_timer(timer_int_quadrupel_2cff(integralpar_i_int_part))

    call int_data_2cff_setup()

  end subroutine setup
  !**************************************************************


  !*************************************************************
  subroutine shutdown
    !  Purpose: contains shutdown routines of various modules
    !------------ Executable code ------------------------------

    call int_data_2cff_shutdown()

    call stop_timer(timer_int_calcsum_2cff(integralpar_i_int_part))
    call stop_timer(timer_int_quadrupel_2cff(integralpar_i_int_part))
    call start_timer(timer_int_idle_2cff(integralpar_i_int_part))

  end subroutine shutdown
  !**************************************************************

  subroutine say(phrase)
    use output_module, only: output_int_loops
    use iounitadmin_module, only: write_to_output_units
    implicit none
    character(len=*), intent(in) :: phrase
    ! *** end of interface ***

    if( output_int_loops ) then
        call write_to_output_units("integral_calc_quad_2cff: "//phrase)
    endif
  end subroutine say

end subroutine integral_calc_quad_2cff

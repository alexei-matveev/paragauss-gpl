!
! ParaGauss, a program package for high-performance computations
! of molecular systems
! Copyright (C) 2014
! T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
! M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
! A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
! T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
! M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
! M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License version 2 as published
! by the Free Software Foundation [1].
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
module mm_timer_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: timer_data
     real(r8_kind) :: tot_rt,tot_ut,tot_st
     real(r8_kind) :: rt,ut,st
     real(r8_kind) :: p_rt,p_ut,p_st
     integer(i4_kind) :: N_calls, N_calls_p
  end type timer_data


  type(timer_data), public :: full_time !Full time of job
  type(timer_data), public :: nbl_time  !generation of non-bonded lists
  type(timer_data), public :: cis_time  !calculation of image species
  type(timer_data), public :: optim_time  !full optimization time
  type(timer_data), public :: single_point_time !full single_point time
  type(timer_data), public :: ebonding_time !calculation of bonded interactions
  type(timer_data), public :: vdw_time !calculation of van der Waals interactions
  type(timer_data), public :: coulomb_time !calculation of interaction of non-periodic
                                           !point charges
  type(timer_data), public :: dipdip_time  !calculation time of dipole-dipole interactions
  type(timer_data), public :: ewald_r_time !Ewald method: real-space contribution
  type(timer_data), public :: ewald_d_time !Ewald method: reciprocal-space contribution
  type(timer_data), public :: comm_time !communication time
  type(timer_data), public :: hinv_time !time of hessian inversion

  !------------ public functions and subroutines ------------------
  public init_timer, start_mm_timer, stop_mm_timer, print_timing, print_time
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  type :: curr_time
     real(r8_kind) :: real_time
     real(r8_kind) :: user_time
     real(r8_kind) :: syst_time
  end type curr_time

  type(curr_time) :: initial_values

  type :: date_time
     integer(r4_kind) :: value(8)
     character(len=10) :: real_clock(3)
  end type date_time

  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine init_timer()

    initial_values=get_time()

    full_time%N_calls=0; full_time%N_calls_p=0
    full_time%tot_rt=zero
    full_time%tot_ut=zero
    full_time%tot_st=zero

    optim_time%N_calls=0; optim_time%N_calls_p=0
    optim_time%tot_rt=zero
    optim_time%tot_ut=zero
    optim_time%tot_st=zero

    single_point_time%N_calls=0; single_point_time%N_calls_p=0
    single_point_time%tot_rt=zero
    single_point_time%tot_ut=zero
    single_point_time%tot_st=zero

    ewald_r_time%N_calls=0; ewald_r_time%N_calls_p=0
    ewald_r_time%tot_rt=zero
    ewald_r_time%tot_ut=zero
    ewald_r_time%tot_st=zero

    ewald_d_time%N_calls=0; ewald_d_time%N_calls_p=0
    ewald_d_time%tot_rt=zero
    ewald_d_time%tot_ut=zero
    ewald_d_time%tot_st=zero

    vdw_time%N_calls=0; vdw_time%N_calls_p=0
    vdw_time%tot_rt=zero
    vdw_time%tot_ut=zero
    vdw_time%tot_st=zero

    coulomb_time%N_calls=0; coulomb_time%N_calls_p=0
    coulomb_time%tot_rt=zero
    coulomb_time%tot_ut=zero
    coulomb_time%tot_st=zero

    dipdip_time%N_calls=0; dipdip_time%N_calls_p=0
    dipdip_time%tot_rt=zero
    dipdip_time%tot_ut=zero
    dipdip_time%tot_st=zero

    ebonding_time%N_calls=0; ebonding_time%N_calls_p=0
    ebonding_time%tot_rt=zero
    ebonding_time%tot_ut=zero
    ebonding_time%tot_st=zero

    nbl_time%N_calls=0; nbl_time%N_calls_p=0
    nbl_time%tot_rt=zero
    nbl_time%tot_ut=zero
    nbl_time%tot_st=zero

    cis_time%N_calls=0; cis_time%N_calls_p=0
    cis_time%tot_rt=zero
    cis_time%tot_ut=zero
    cis_time%tot_st=zero

    comm_time%N_calls=0; comm_time%N_calls_p=0
    comm_time%tot_rt=zero
    comm_time%tot_ut=zero
    comm_time%tot_st=zero

    hinv_time%N_calls=0; hinv_time%N_calls_p=0
    hinv_time%tot_rt=zero
    hinv_time%tot_ut=zero
    hinv_time%tot_st=zero

  end subroutine init_timer
  !****************************************************************

  !****************************************************************
  function get_time()

#ifdef FPP_AIX_XLF
   use xlfutility
#endif

    type(curr_time) :: get_time

    type(date_time) :: dt
#ifdef FPP_AIX_XLF
    type(tb_type) :: buffer
    real(r8_kind) :: rt
#else
    real(r4_kind) :: tarray(2),tsum
#endif

#ifndef INTRINSIC_ETIME
  real(r4_kind), external :: etime
#endif

    call date_and_time(dt%real_clock(1),dt%real_clock(2),dt%real_clock(3), &
         dt%value)
    get_time%real_time=dt%value(8)*0.001_r8_kind+ &
         dt%value(7)*1.0_r8_kind+dt%value(6)*60.0_r8_kind+ &
         dt%value(5)*3600.0_r8_kind+dt%value(3)*86400.0_r8_kind
#ifdef FPP_AIX_XLF
    rt=real(etime_(buffer) ,kind=r8_kind)
    get_time%user_time=real(buffer%usrtime ,kind=r8_kind)
    get_time%syst_time=real(buffer%systime ,kind=r8_kind)
#else
    tsum=etime(tarray)
    get_time%user_time=tarray(1)
    get_time%syst_time=tarray(2)
#endif

  end function get_time
  !****************************************************************

  !****************************************************************
  subroutine start_mm_timer(mm_timer)

    type(timer_data) :: mm_timer
    type(curr_time) :: ct

    if(mm_timer%N_calls /= mm_timer%N_calls_p) call error_handler &
         ("MolMech: start_mm_timer could not be started")
    ct=get_time()

    mm_timer%p_rt=ct%real_time
    mm_timer%p_ut=ct%user_time
    mm_timer%p_st=ct%syst_time
    mm_timer%N_calls=mm_timer%N_calls+1

  end subroutine start_mm_timer
  !****************************************************************

  !****************************************************************
  subroutine stop_mm_timer(mm_timer)

    type(timer_data) :: mm_timer
    type(curr_time) :: ct

    ct=get_time()

    mm_timer%rt=ct%real_time-mm_timer%p_rt
    mm_timer%ut=ct%user_time-mm_timer%p_ut
    mm_timer%st=ct%syst_time-mm_timer%p_st

    mm_timer%tot_rt=mm_timer%tot_rt+mm_timer%rt
    mm_timer%tot_ut=mm_timer%tot_ut+mm_timer%ut
    mm_timer%tot_st=mm_timer%tot_st+mm_timer%st

    mm_timer%N_calls_p=mm_timer%N_calls

  end subroutine stop_mm_timer
  !****************************************************************

  !****************************************************************
  subroutine print_timing()

    write(output_device,'(/,a13)') "TIMING OF JOB"
    write(output_device,'(80("-"))')
    write(output_device,'(a4,40x,a6,3x,a4,6x,a4,5x,a6)') "Task","Number","Real","User","System"
    write(output_device,'(80("-"))')

    if(nbl_time%N_calls > 0) call print_time(nbl_time,"Generation of nonbonded lists")
    if(cis_time%N_calls > 0) call print_time(cis_time,"Calculation of image species")
    if(ebonding_time%N_calls > 0) call print_time(vdw_time,"Bonding interactions")
    if(vdw_time%N_calls > 0) call print_time(vdw_time,"Van der Waals interaction")
    if(coulomb_time%N_calls > 0) call print_time(coulomb_time,"Coulomb interaction")
    if(dipdip_time%N_calls > 0) call print_time(dipdip_time,"Dipole-dipole interaction")
    if(ewald_d_time%N_calls > 0) call print_time(ewald_d_time,"Ewald: direct part")
    if(ewald_r_time%N_calls > 0) call print_time(ewald_r_time,"Ewald: reciprocal part")
    if(hinv_time%N_calls > 0) call print_time(hinv_time,"Calculation of inverse Hessian")
    if(single_point_time%N_calls > 0) call print_time(single_point_time, &
         "Single point: energy, gradients and hessian")
    if(calc_optimization) call print_time(optim_time,"Optimization")
    if(comm_time%N_calls > 0) call print_time(comm_time,"Communication time")
    call print_time(full_time,"WHOLE JOB")

  end subroutine print_timing
  !****************************************************************

  !****************************************************************
  subroutine print_time(mm_timer,timer_name)

    type(timer_data) :: mm_timer
    character(len=*) :: timer_name
    
    character(len=44) :: message

    message=trim(timer_name)

    write(output_device,'(a44,a1,i5,$)') message,":",mm_timer%N_calls
    write(output_device,'(3f10.2)') mm_timer%tot_rt, mm_timer%tot_ut, mm_timer%tot_st
    if(mm_timer%N_calls > 1) then
       write(output_device,'(a50,$)') "per task:     "
       write(output_device,'(3f10.2)') mm_timer%tot_rt/mm_timer%N_calls, &
            mm_timer%tot_ut/mm_timer%N_calls, mm_timer%tot_st/mm_timer%N_calls
    end if

  end subroutine print_time
  !****************************************************************

  !****************************************************************
end module mm_timer_module

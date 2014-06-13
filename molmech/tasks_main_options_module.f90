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
module tasks_main_options_module

  !------------ Modules used --------------------------------------
  use type_module
  use inp_out_module
  use common_data_module
  use molmech_msgtag_module
  use comm_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  logical, public :: calc_energy, calc_gradients, calc_optimization, &
       calc_hessian, calc_ppc_array, solvent
  logical, public :: to_xyz_file,to_viewmol,to_xmkmol_file
  logical, public :: read_gx
  logical, public :: write_gx
  logical, public :: to_epe
  logical, public :: to_pg_hess
  logical, public :: with_optimizer

  character(len=100), public :: job_topic

  real(kind=r8_kind), public :: rcut_ew
  !list of main options
  logical, public :: calc_strain
  real(kind=r8_kind), public :: rcut_vdw,eps,rcut_cs
  logical, public :: coulomb,van_der_waals
  character(len=10), public :: energy_unit
  character(len=10), public :: length_unit
  character(len=10), public :: angle_unit
  logical, public :: minimal_image
  real(kind=r8_kind), public :: rc_inc

  integer(i4_kind), public :: qmmm
  !------------ public functions and subroutines ------------------
  public read_topic, write_topic_to_output, read_tasks, &
       write_tasks_to_output,read_main_options,write_options_to_output, &
       send_receive_tasks_and_options
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  logical, parameter :: df_calc_energy = .false.
  logical, parameter :: df_calc_gradients = .false.
  logical, parameter :: df_calc_optimization = .false.
  logical, parameter :: df_calc_hessian = .false.
  logical, parameter :: df_calc_ppc_array = .false.
  logical, parameter :: df_solvent = .false.
  logical, parameter :: df_to_xyz_file = .false.
  logical, parameter :: df_to_xmkmol_file = .false.
  logical, parameter :: df_to_viewmol = .false.
  logical, parameter :: df_read_gx = .false.
  logical, parameter :: df_write_gx = .false.
  logical, parameter :: df_to_epe = .false.
  logical, parameter :: df_to_pg_hess = .false.
  logical, parameter :: df_with_optimizer = .false.

  character(len=100) :: job_tasks
  namelist /tasks/ job_tasks 

  namelist /topic/ job_topic

  logical, parameter :: df_calc_strain=.false.
  real(kind=r8_kind), parameter :: df_rcut_vdw=infinity
  real(kind=r8_kind), parameter :: df_rcut_cs=one
  real(kind=r8_kind), parameter :: df_eps=one
  logical, parameter :: df_coulomb=.true. 
  logical, parameter :: df_van_der_waals=.true. 
  character(len=10), parameter :: df_energy_unit="KJ/MOL" !KCAL/MOL, EV
  character(len=10), parameter :: df_length_unit="ANGSTROM" !BOHR
  character(len=10), parameter :: df_angle_unit="DEGREE" ! RADIAN
  logical, parameter :: df_minimal_image=.true.
  real(kind=r8_kind), parameter :: df_rc_inc=one    !incriments sum of covalent radii of atom
                                                    !to define chemical bond

  namelist /main_options/ coulomb,van_der_waals,rcut_vdw,eps, &
       energy_unit,length_unit,angle_unit,rc_inc, &
       minimal_image,rcut_cs,calc_strain
  !------------ Subroutines -----------------------------------------
contains
  !******************************************************************
  function read_topic()

    logical :: read_topic

    integer(i4_kind) :: i

    job_topic="..."
    call  go_to_first_input_line
    read_topic=find_namelist("&TOPIC",i)
    if(.not.read_topic) return
    read(input_device,nml=topic, end=100, err=200)
    if(trim(job_topic)=="...") then
       read_topic = .false.
    else
       read_topic = .true.
    end if
    return
100 read_topic = .false.
    return
200 call input_nm_error(0,"TOPIC")

  end function  read_topic
  !******************************************************************
  
  !******************************************************************
  subroutine write_topic_to_output()

    write(output_device,'(80("*"))') 
    write(output_device,'(a2,a76,a2)') '* ',job_topic,' *'
    write(output_device,'(80("*"))') 

  end subroutine write_topic_to_output
  !******************************************************************
  
  !******************************************************************
  function read_tasks()

    logical :: read_tasks

    integer(i4_kind) :: i

    calc_energy=df_calc_energy
    calc_gradients=df_calc_gradients
    calc_optimization=df_calc_optimization
    calc_hessian=df_calc_hessian
    calc_ppc_array=df_calc_ppc_array
    solvent=df_solvent
    to_xyz_file=df_to_xyz_file
    to_xmkmol_file=df_to_xmkmol_file
    to_viewmol=df_to_viewmol
    read_gx=df_read_gx
    write_gx=df_write_gx
    to_epe=df_to_epe
    to_pg_hess=df_to_pg_hess
    with_optimizer=df_with_optimizer

    if(qmmm == 2 .or. qmmm == 3) then
       calc_energy=.true.
       calc_gradients=.true.
       read_tasks = .true.
       return
    end if
    if(qmmm == -1) then
       calc_energy=.true.
       calc_gradients=.true.
       solvent=.false.
       with_optimizer=.false.
       read_tasks = .true.
       return
    end if
    
    call go_to_first_input_line
    read_tasks=find_namelist("&TASKS",i)
    if(.not.read_tasks) return
    read(input_device,nml=tasks, end=100, err=200)

    call upcase(job_tasks)
    if(check_string(job_tasks,'ENERGY')) then
       calc_energy=.true.
    elseif(check_string(job_tasks,'GRADIENTS')) then
       calc_energy=.true.
       calc_gradients=.true.
    elseif(check_string(job_tasks,'OPTIMIZATION')) then
       calc_energy=.true.
       calc_gradients=.true.
       calc_optimization=.true.
    elseif(check_string(job_tasks,'HESSIAN')) then
       calc_energy=.true.
       calc_gradients=.true.
       calc_hessian=.true.
    end if

    if(calc_gradients) then
       if(check_string(job_tasks,'WITH_OPTIM')) with_optimizer=.true.
    end if

    if(check_string(job_tasks,'PERIODIC_PC_ARRAY')) then
       calc_ppc_array=.true.
    end if

    if(.not.calc_energy .and. .not.calc_ppc_array) then
       call error_handler("MolMech: There are no main TASKS. What shuold I do? Think twice")
    end if

    if(check_string(job_tasks,'TO_XYZ')) then
       to_xyz_file=.true.
    end if
    if(check_string(job_tasks,'TO_XMAKEMOL')) then
       to_xmkmol_file=.true.
    end if
    if(check_string(job_tasks,'TO_VIEWMOL')) then
       to_viewmol=.true.
    end if
    if(check_string(job_tasks,'READ_GX')) then
       read_gx=.true.
    end if
    if(check_string(job_tasks,'WRITE_GX')) then
       write_gx=.true.
    end if
    if(check_string(job_tasks,'TO_EPE')) then
       to_epe=.true.
    end if
    if(check_string(job_tasks,'TO_PG_HESS')) then
       to_pg_hess=.true.
    end if
    if(check_string(job_tasks,'SOLVENT')) then
       solvent=.true.
    end if

    if(solvent) then
       calc_hessian=.false.
       calc_ppc_array=.false.
       to_pg_hess=.false.
    end if

    if(with_optimizer) then
       calc_optimization=.false.
       calc_hessian=.false.
       read_gx=.true.
       write_gx=.true.
    end if

    if((to_pg_hess .and. .not.calc_optimization) .and. &
         (to_pg_hess .and. .not.calc_hessian)) then
       write(output_device,*) "TO_PG_HESS option can be used only in combination with" &
            //achar(10)//"either OPTIMIZATION or HESSIAN options. Now ignored!"
       to_pg_hess=.false.
    end if

    read_tasks = .true.
    return

100 read_tasks = .false.
    return
200 call input_nm_error(0,"TASKS") 

  end function  read_tasks
  !******************************************************************
  
  !******************************************************************
  subroutine write_tasks_to_output()

    character(len=76) :: message

    write(output_device,'(80("*"))') 
    message='Current actions:'
    write(output_device,'(a2,a76,a2)') '* ',message,' *'
    if (calc_optimization) then
       message='OPTIMIZATION:    geometry optimization'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    else if (calc_hessian) then
       message='HESSIAN:       Single point calculation - Energy + forces + hessian'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    else if (calc_gradients) then
       message='GRADIENTS:       Single point calculation - Energy + forces'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    else if (calc_energy) then
       message='ENERGY:          Single point calculation - only Energy'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(solvent) then
       message='Calulation in solvent'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if (calc_ppc_array) then
       message='POINT CHARGE ARRAY: 3-D potential reproducing'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(read_gx) then
       write(output_device,'(a45,34x,a1)') "* Initial geometry will be taken from GX file","*"
    end if
    if(write_gx .and. calc_gradients) then
       write(output_device,'(a63,16x,a1)') &
            "* Final geometry, energy and gradients will be saved in GX file","*"
    end if
    if(to_xyz_file) then
       write(output_device,'(a42,37x,a1)') "* Final geometry will be saved in XYZ file","*"
    end if
    if(to_viewmol) then
       write(output_device,'(a48,31x,a1)') "* Final geometry will be saved in Viewmol format","*"
    end if
    if(to_epe) then
       write(output_device,'(a51,28x,a1)') "* Final geometry will be saved to be read in by EPE","*"
    end if
    write(output_device,'(a47,32x,a1)') "* Internal units: kJ/mol and angstrom - FOREVER","*"
    write(output_device,'(80("*"),/)') 
    
  end subroutine write_tasks_to_output
  !******************************************************************
  
  !******************************************************************
  function read_main_options()

    logical :: read_main_options

    integer(i4_kind) :: i

    calc_strain=df_calc_strain
    rcut_vdw=df_rcut_vdw
    rcut_cs=df_rcut_cs
    rc_inc=df_rc_inc
    eps=df_eps
    coulomb=df_coulomb
    van_der_waals=df_van_der_waals
    energy_unit=df_energy_unit
    length_unit=df_length_unit
    angle_unit=df_angle_unit
    minimal_image=df_minimal_image

    call  go_to_first_input_line
    read_main_options=find_namelist("&MAIN_OPTIONS",i)
    if(.not.read_main_options) return
    read(input_device,nml=main_options, end=100, err=200)
    call upcase(energy_unit)
    if(trim(energy_unit) /= "KJ/MOL" .and.  &
       trim(energy_unit) /= "KCAL/MOL" .and. &
       trim(energy_unit) /= "EV") goto 200
    call upcase(length_unit)
    if(trim(length_unit) /= "ANGSTROM" .and.  trim(length_unit) /= "BOHR") goto 200
    call upcase(angle_unit)
    if(trim(angle_unit) /= "DEGREE" .and.  trim(angle_unit) /= "RADIAN") goto 200

    if(solvent .and. calc_strain) calc_strain=.false.

    read_main_options=.true.
    return

100 read_main_options=.false.
    return
200 call input_nm_error(0,"MAIN_OPTIONS")

  end function read_main_options
  !******************************************************************
  
  !******************************************************************
  subroutine write_options_to_output()

    character(len=76) :: message

    write(output_device,'(80("*"))') 
    message='Main options:'
    write(output_device,'(a2,a76,a2)') '* ',message,' *'
    write(output_device,'(a23,f12.6,44x,a1)') '* Dielectric constant= ',eps,'*'
    write(output_device,'(a38,f10.4,a9,22x,a1)') '* Global van der Waals cutoff radius= ',rcut_vdw, &
         ' angstrom','*'
    write(output_device,'(a45,f10.2,a9,15x,a1)') '* Cutoff radius to define core-shell pairs is ', &
         rcut_cs, ' angstrom','*'
    write(output_device,'(a42,f5.2,32x,a1)') '* Sum of atomic covalent radii increment= ', &
         rc_inc,'*'
    if(.not.coulomb) then
       message='No electrostatic interaction'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(.not.van_der_waals) then
       message='No van der Waals interaction'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    write(output_device,'(a15,a10,a2,a10,a5,a10,27x,a1)') "* Input units: ",energy_unit,", ", &
         length_unit," and ",angle_unit,"*"
    write(output_device,'(80("*"))') 

  end subroutine write_options_to_output
  !******************************************************************
  
  !******************************************************************
  subroutine send_receive_tasks_and_options()

    integer(i4_kind) :: info

    if(comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_mm_initslaves)
       call commpack(calc_gradients,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_gradients pack failed")
       call commpack(calc_hessian,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_hessian pack failed")
       call commpack(calc_strain,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_strain pack failed")
       call comm_send()
    else
       call communpack(calc_gradients,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_gradients unpack failed")
       call communpack(calc_hessian,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_hessian unpack failed")
       call communpack(calc_strain,info)
       if( info /= 0) call error_handler &
            ("send_receive_tasks_and_options: calc_strain unpack failed")
    end if

  end subroutine send_receive_tasks_and_options
  !******************************************************************
  
  !******************************************************************
end module tasks_main_options_module



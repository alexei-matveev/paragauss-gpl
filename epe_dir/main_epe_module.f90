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
module main_epe_module
#include <def.h>
  use epecom_module
  use culon_module
  use mol_module
  use str_module 
  use epe_pg_module
  use epe_pot_module
  use iounitadmin_module
  use time_module
  use timer_module
  use comm_module
  use msgtag_module, only: &
                     msgtag_epe_consdef &
                   , msgtag_epe_defects &
                   , msgtag_epe_def_fin &
                   , msgtag_epe_finish_slave &
                   , msgtag_epe_forces_done &
                   , msgtag_epe_grads_cons &
                   , msgtag_epe_grads_init &
                   , msgtag_epe_init_slave &
                   , msgtag_epe_send_def &
                   , msgtag_epe_send_only &
                   , msgtag_imp_intermediate &
                   , msgtag_vac_intermediate
  
  implicit none
  private
  save
  public:: main_epe, start_epe, n_epe_vacancies, &
       finish_epe,                           &
       epe_send_init_to_slave,               &
       epe_send_finish_to_slave,             &
       epe_send_shells_and_cores,            &
       epe_receive_shells_and_cores,         & 
       init_lat_gradients,                   &
       cons_latt_gradients,                  &
       cons_defect_contrubutions,            &
       defect_contributions,                 &
       epe_init_slave, epe_finish_slave,     &
       start_slave_defect,                   &
       epe_collect_def_contributions,        &
       print_energy,                         &
       defect_contributions_fin

        type imp_grads
        real(kind=r8_kind), pointer, dimension(:,:) :: imp_core_grad,imp_shell_grad
        end type imp_grads
        type(imp_grads) prev,curr
  logical:: epein,not_first_step   
  real(kind=r8_kind), allocatable,dimension(:,:)        :: short_range_grad
  real(kind=r8_kind), allocatable,dimension(:,:)        :: reg1_shells_coordinates, &
                                                           reg1_cores_coordinates
  real(kind=r8_kind),dimension(:,:),allocatable:: d_r_sh,d_sh_nuc
  real(kind=r8_kind) :: etotal_gopt,ETM,EPL,E3B,dels,DISPM,DISTM,pol
  real(kind=r8_kind) :: CKALP,CKALT,UGOLT,UGOLP,RLIND,SPX,SPY,SPZ,CPX,CPY,CPZ
  integer(kind=i4_kind) :: io,ia,k,in,i,j,ity,ind,l
  integer(kind=i4_kind) :: L1,l3
  integer(kind=i4_kind) :: hdscom_unit,in_unit,out_unit,ucellxyz
  integer(kind=i4_kind) :: n_epe_vacancies
  integer(kind=i4_kind) :: alloc_status ,    count1=0
  integer(kind=i4_kind), parameter ::  msg_init=0,     & 
                                       msg_cons=1,     &
                                       msg_defects=2,  &
                                       msg_consdef=3,  &
                                       msg_def_fin=4
!******************************************************
contains

  subroutine calc_lattice
    CALL trans_vectors
 !  print*,'CALL epe_generater(N_GEN_IONS,.true.)'
    CALL epe_generater(N_GEN_IONS,.true.)
 !  print*,'epe_generater dene'
    if(.not. core_shell) then
       R_NUC_ION(1:reg_2a_n_ions,1)=epe(1:reg_2a_n_ions)%r(1)
       R_NUC_ION(1:reg_2a_n_ions,2)=epe(1:reg_2a_n_ions)%r(2)
       R_NUC_ION(1:reg_2a_n_ions,3)=epe(1:reg_2a_n_ions)%r(3)
       R_SH_ION(1:reg_2a_n_ions,:)=R_NUC_ION(1:reg_2a_n_ions,:)
    else
       R_NUC_ION(1:reg_2a_n_ions,1)=epe(1:reg_2a_n_ions)%c(1)
       R_NUC_ION(1:reg_2a_n_ions,2)=epe(1:reg_2a_n_ions)%c(2)
       R_NUC_ION(1:reg_2a_n_ions,3)=epe(1:reg_2a_n_ions)%c(3)
       R_SH_ION(1:reg_2a_n_ions,1)=epe(1:reg_2a_n_ions)%s(1)
       R_SH_ION(1:reg_2a_n_ions,2)=epe(1:reg_2a_n_ions)%s(2)
       R_SH_ION(1:reg_2a_n_ions,3)=epe(1:reg_2a_n_ions)%s(3)
    endif
    
       CALL madelung
    if(.not. periodic_optimization) then
       CALL mott_litleton
       if(region_I_output) call regionI2xyz()
    endif

  end subroutine calc_lattice
!**********************************************************

!**********************************************************
  subroutine ucell2xyz

    integer(kind=i4_kind) ::i
    logical :: lucell

    inquire (file=trim(epe_input_dir)//'/xyz_ucell',exist=lucell)
    if(.not.lucell) then
       ucellxyz=openget_iounit(trim(epe_input_dir)//'/xyz_ucell',  &
          form='formatted', status='unknown')
    end if

    write(ucellxyz,'(i4)') n_ions_cell
    write(ucellxyz,'(a9)') 'UNIT_CELL'
    do i=1,n_ions_cell
       write(ucellxyz,'(a2,3(1x,f12.8))') NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
    end do

  end subroutine ucell2xyz
!**********************************************************

!**********************************************************
  subroutine regionI2xyz

    integer(kind=i4_kind) :: rIxyz,i,k

    rIxyz=openget_iounit(trim(epe_input_dir)//'/xyz_reg_I',  &
          form='formatted', status='unknown')

    write(rIxyz,'(i4)') reg_I_n_ions+n_fixed
    write(rIxyz,'(a8)') 'region_I'
    do i=1,reg_I_n_ions+n_fixed
       k=epe(i)%k
       if(periodic_optimization) then
          write(rIxyz,'(a2,3(1x,f12.8))') name_of_type(k), epe(i)%r
       else
          write(rIxyz,'(a2,3(1x,f12.8))') name_of_type(k), epe(i)%c
       endif
    end do

    call returnclose_iounit(rIxyz)

  end subroutine regionI2xyz
!**********************************************************
!**********************************************************
  subroutine regionIs2xyz

    integer(kind=i4_kind) :: rIxyz,i,k

    rIxyz=openget_iounit(trim(epe_input_dir)//'/xyz_reg_Is',  &
          form='formatted', status='unknown')

    write(rIxyz,'(i4)') reg_I_n_ions+n_fixed
    write(rIxyz,'(a8)') 'region_I'
    do i=1,reg_I_n_ions+n_fixed
       k=epe(i)%k
       write(rIxyz,'(a2,3(1x,f12.8))') name_of_type(k), epe(i)%s
    end do

    call returnclose_iounit(rIxyz)

  end subroutine regionIs2xyz
!**********************************************************

!**********************************************************
   subroutine start_epe
     integer(kind=i4_kind)::reg2a_unit,counter
     ! part of code: generates and allocates space for EPE
     ! if exist and asked then previous configuration of reg. I
     ! is read

        ! on the very first step calc of fefect contributions bypassed thus
        defect_energy_short=0.0_r8_kind
        defect_energy_coul=0.0_r8_kind

     if(periodic_optimization) read_configuration=.false.
     if(periodic_optimization) write_configuration=.false.
     make_epe_reference=make_epe_reference 
     pg_interfaced_mode=0

     if(qm_interfaced_mode.or.operations_make_reg_reference) call pggener()

     not_first_step=.not.qm_interfaced_mode

     call calc_lattice !(1)
 !   print*,'calc_lattice 1 done'

     if(reg_2a_treated) then
       allocate (reg_I_pg(reg_2a_n_ions),stat=epealloc_stat(21))
     else
        allocate (reg_I_pg(reg_I_n_ions),stat=epealloc_stat(21))
     end if
     ASSERT(epealloc_stat(21).eq.0)
     epealloc_stat(21)=1
  
     allocate (point_0_core_grad(reg_I_n_ions,3) &
          ,point_0_shell_grad(reg_I_n_ions,3) &
          ,point_1_core_grad(reg_I_n_ions,3) &
          ,point_1_shell_grad(reg_I_n_ions,3) &
          ,short_range_grad(reg_I_n_ions,3), stat=epealloc_stat(23) )
     ASSERT(epealloc_stat(23).eq.0)
     epealloc_stat(23)=1

     short_range_grad=0.0_r8_kind
     point_0_core_grad=0.0_r8_kind
     point_0_shell_grad=0.0_r8_kind
     point_1_core_grad=0.0_r8_kind
     point_1_shell_grad=0.0_r8_kind
     allocate(reg1_shells_coordinates(reg_I_n_ions,3), &
              reg1_cores_coordinates(reg_I_n_ions,3) ,stat=epealloc_stat(24))
     ASSERT(epealloc_stat(24).eq.0)
     epealloc_stat(24)=1

     if(.not.allocated(regI_previous)) then 
        allocate(regI_previous(1:reg_I_n_ions),stat=alloc_status)
        if(alloc_status.ne.0) call error_handler(" allocate regI_previous failed")
     endif

     if(n_vacancies.ne.0) then
        regI_previous(:)%r_core(1)=epe(1:reg_I_n_ions)%r(1)
        regI_previous(:)%r_core(2)=epe(1:reg_I_n_ions)%r(2)
        regI_previous(:)%r_core(3)=epe(1:reg_I_n_ions)%r(3)
        regI_previous(:)%r_shell(1)=epe(1:reg_I_n_ions)%r(1)
        regI_previous(:)%r_shell(2)=epe(1:reg_I_n_ions)%r(2)
        regI_previous(:)%r_shell(3)=epe(1:reg_I_n_ions)%r(3)
     end if

 !  print*,' read /epeinout'
     inquire(file=trim(epe_input_dir)//'/epeinout',exist=epein)
      read_configuration=read_configuration.and.epein
      if(read_configuration)  call read_config
 !  print*,'done'

     select case(basic_action)
     case(0:1)
        if(qm_interfaced_mode) call get_pgcent
        ! in particular reads regular and reference configurations
        ! r_imp is defined

        n_epe_vacancies=n_vacancies

        if(n_impurities.ne.0) then
              allocate(impurities_old_conf(n_impurities),&
             & impurities_new_conf(n_impurities),stat=epealloc_stat(22))
        ASSERT(epealloc_stat(22).eq.0)

        allocate(curr%imp_core_grad(n_impurities,3), curr%imp_shell_grad(n_impurities,3), &
                 prev%imp_core_grad(n_impurities,3), prev%imp_shell_grad(n_impurities,3), &
                 stat=epealloc_stat(15))
                 ASSERT(epealloc_stat(15).eq.0)
                        epealloc_stat(15)=1

     curr%imp_shell_grad=0.0_r8_kind
     curr%imp_core_grad=0.0_r8_kind
     prev%imp_shell_grad=0.0_r8_kind
     prev%imp_core_grad=0.0_r8_kind

        else 
        allocate(curr%imp_core_grad(0,0), curr%imp_shell_grad(0,0), &
                 prev%imp_core_grad(0,0), prev%imp_shell_grad(0,0), &
                 stat=epealloc_stat(15))
                 ASSERT(epealloc_stat(15).eq.0)
                        epealloc_stat(15)=1
           
     endif

        IF(read_configuration) then
         inquire (file=trim(epe_input_dir)//'/reg2a_in',exist=ml_displacements_fixed)  
        if(ml_displacements) then
           if(ml_cluster_simulated) then
!!$     print*,'start_epe: calculate reg2a_ml_displacements ...'
              CALL reg2a_ml_displacements
              elseif(ml_displacements_fixed) then
!!$     print*,'start_epe: read already fixed ml displacements...'              

                   reg2a_unit=openget_iounit( &
                      trim(epe_input_dir)//'/reg2a_in',status='old') 
                 do counter=reg_I_n_ions+1,reg_2a_n_ions
                    read(reg2a_unit,*) epe(counter)%s
                    read(reg2a_unit,*) epe(counter)%c
                 end do
                 call returnclose_iounit(reg2a_unit)  
                 reg_2a_treated=.false.
              end if
              endif
           ! after reading configuration of reg I make ml_displacements 
        endif

     endselect

!--------------------------------------------
 ! print*,'start EPE done'
   contains

     subroutine read_config
       integer(kind=i4_kind) :: n_saved_ions,i
       logical :: c3_symm,yes
! **make binary input of previosly saved configuration
       inquire(file=trim(epe_input_dir)//'/epeinout',exist=epein)
 !     print*,'read_config'
       if(.not.epein) then
          print*,'file ',   trim(epe_input_dir)//'/epeinout','not found'
       end if
  
       in_unit=openget_iounit(trim(epe_input_dir)//'/epeinout', &
                                 form='unformatted',status='old') 
       DPRINT 'epeinout unit assigned to ',in_unit
       READ(in_unit,err=100)n_saved_ions
       if(n_saved_ions.eq.reg_I_n_ions) then
          READ(in_unit,err=100)R_SH_ION(:reg_I_n_ions,:),R_NUC_ION(:reg_I_n_ions,:)
#ifdef NEW_EPE
          READ(in_unit,err=200,end=200) c3_symm
          yes=(c3_symm.and.option_c3_symm).or.(.not.c3_symm.and..not.option_c3_symm)
          ASSERT(yes)
200       continue
          READ(in_unit,err=300,end=300)(epe(i)%o,i=1,reg_I_n_ions)
#endif

          call sorting_epeinout(reg_I_n_ions)
#ifdef NEW_EPE
300       continue
#endif
          epe(:reg_I_n_ions)%s(1)=R_SH_ION(:reg_I_n_ions,1)
          epe(:reg_I_n_ions)%s(2)=R_SH_ION(:reg_I_n_ions,2)
          epe(:reg_I_n_ions)%s(3)=R_SH_ION(:reg_I_n_ions,3)
          epe(:reg_I_n_ions)%c(1)=R_NUC_ION(:reg_I_n_ions,1)
          epe(:reg_I_n_ions)%c(2)=R_NUC_ION(:reg_I_n_ions,2)
          epe(:reg_I_n_ions)%c(3)=R_NUC_ION(:reg_I_n_ions,3)
      if(.not.allocated(epe_at_start)) then
           allocate(epe_at_start(reg_I_n_ions),stat=epealloc_stat(6))
           ASSERT(epealloc_stat(6).eq.0)
           epealloc_stat(6)=1
      endif
          epe_at_start(:reg_I_n_ions)%s(1)=R_SH_ION(:reg_I_n_ions,1)
          epe_at_start(:reg_I_n_ions)%s(2)=R_SH_ION(:reg_I_n_ions,2)
          epe_at_start(:reg_I_n_ions)%s(3)=R_SH_ION(:reg_I_n_ions,3)
          epe_at_start(:reg_I_n_ions)%c(1)=R_NUC_ION(:reg_I_n_ions,1)
          epe_at_start(:reg_I_n_ions)%c(2)=R_NUC_ION(:reg_I_n_ions,2)
          epe_at_start(:reg_I_n_ions)%c(3)=R_NUC_ION(:reg_I_n_ions,3)

          DPRINT 'read_config: configuration of ions in region I is read from file'
          call returnclose_iounit(in_unit)
       else
          call returnclose_iounit(in_unit)  
          print*,n_saved_ions,reg_I_n_ions
          ASSERT(n_saved_ions==reg_I_n_ions)
       endif
       if(basic_action.eq.0.and.moving_epecenter.ne.0) then
       DPRINT  'initial cordinates of moving_epecenter' 
       DPRINT epe(moving_epecenter)%s
       DPRINT epe(moving_epecenter)%c
        epe(moving_epecenter)%s=coordinates_moving_epecenter
        R_SH_ION(moving_epecenter,:)=coordinates_moving_epecenter
        epe(moving_epecenter)%c=coordinates_moving_epecenter
        R_NUC_ION(moving_epecenter,:)=coordinates_moving_epecenter
        DPRINT 'cordinates of moving_epecenter',moving_epecenter,' are set to',epe(moving_epecenter)%s
       endif
       return
100    call error_handler("epein(epeinout) has wrong format for this system")
     end subroutine read_config
   end subroutine start_epe
   !**************************************************************
!AG 
   !*************************************************************
   subroutine epe_send_init_to_slave()
     !  Purpose: Tells to slaves to begin epe-initialization.  
     !           Called by Master (epe-driver).
     !------------ Modules used ----------------------------------
        use epepar_module, only: ec_max_type_ions=>max_type_ions
     implicit none
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind)            :: status, i_ion, size_reg_reference, &
          size_epe_reference
     integer(kind=i4_kind) :: n_hosts,n_reg,i,nt
     logical             :: allocated_epe_reference
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
!!$     print*, 'Master: => EPE-init message'
     n_hosts=comm_get_n_processors()
     if(n_types_central_atoms_3body > 0) then
        n_reg=ceiling(n_tetrahedrons*1.0_r8_kind/n_hosts)
        first_ind=1
        last_ind=first_ind+n_reg-1
     endif

     do i=2,n_hosts
        call comm_init_send(i, msgtag_epe_init_slave)
        call commpack(n_gen_ions, status)
DPRINT 'n_gen_ions MASTER',n_gen_ions
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [1]")
        call commpack(reg_I_n_ions, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [2]")
        call commpack(n_vacancies, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [2]")
        call commpack(n_impurities, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [2]")   
        call commpack(reg_2a_n_ions, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [2]")
        call commpack(PIS, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [2]")
        call commpack(erfo, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [3]")
        call commpack(et2, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [4]")
        call commpack(error_function_parameter, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [5]")

        do i_ion=1,n_gen_ions
           call commpack(epe(i_ion)%r(1),3,1,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [6]")
           call commpack(epe(i_ion)%k, status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [7]")
           call commpack(epe(i_ion)%ec, status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [7]")
           call commpack(epe(i_ion)%m, status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [7]")
        end do

        call commpack(radius_long_interact, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [8]")
        call commpack(max_type_ions, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [11]")
        call commpack(fixed_reg2a_relcontribs, status)
        if(status.ne.0) call error_handler &
        ("epe_send_init_to_slave : error fixed_reg2a_relcontribs")        
        call commpack(explicit_coupling, status)
        if(status.ne.0) call error_handler &
        ("epe_send_init_to_slave : error explicit_coupling")
        if(explicit_coupling) then
        call commpack(ec_max_type_ions, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error ec_max_type_ions")
        call commpack(explicit_coupling_type(1),n_impurities,1,status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error explicit_coupling_type")
        endif
        call commpack(q_shell(1), ndt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [12]")
        call commpack(q_nuclear(1), ndt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [13]")
        call commpack(q_ion(1), ndt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [14]")
        call commpack(pk(1), ndt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [15]")
        if(explicit_coupling) then
        nt=max(ec_max_type_ions,max_type_ions)
        else
        nt=max_type_ions
        endif
        call commpack(host%ro(1,1,1),nt*nt*(nt+1) , 1, status) !!!!!!!!!AS
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [16]")
        call commpack(host%b(1,1,1),nt*nt*(nt+1) , 1, status)  !!!!!!!!!!AS
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [17]") 
        call commpack(host%c(1,1,1),nt*nt*(nt+1) , 1, status)  !!!!!!!!!!!AS
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [18]")
        call commpack(host%d(1,1,1),nt*nt*(nt+1) , 1, status)  !!!!!!!!!!!AS
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [19]") 
        call commpack(host%sr1(1,1,1),nt*nt*(nt+1),1,status)  !!!!!!!!!!!AS
        call commpack(host%sr2(1,1,1),nt*nt*(nt+1),1,status)  !!!!!!!!!!!AS

        call commpack(host%k(1,1,1),nt*nt*(nt+1),1,status)  
        call commpack(host%r0(1,1,1),nt*nt*(nt+1),1,status)
        call commpack(host%k1(1,1,1),nt*nt*(nt+1),1,status)
        call commpack(host%r1(1,1,1),nt*nt*(nt+1),1,status)

        call commpack(n_bs_points, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [20]")
        call commpack(GSTR(1,1), 3*ndngv, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [21]")
        call commpack(RSIN(1), ndngv, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [22]")
        call commpack(RCOS(1), ndngv, 1, status) 
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [23]")
        call commpack(madc%emd(1), nmadt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [24]") 
        call commpack(use_epe_reference, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [26]") 
        if(use_epe_reference) then
           size_reg_reference=size(reg_reference)
           call commpack(size_reg_reference, status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [27]") 
           do i_ion=1, size_reg_reference
              call commpack(reg_reference(i_ion)%rc(1), 3,1, status)
              if(status.ne.0) call error_handler &
                   ("epe_send_init_to_slave : error [28]") 
              call commpack(reg_reference(i_ion)%rs(1), 3,1, status)
              if(status.ne.0) call error_handler &
                   ("epe_send_init_to_slave : error [29]") 
           end do
        end if

        allocated_epe_reference=allocated( epe_reference )
        call commpack(allocated_epe_reference, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [30]")
        if( allocated_epe_reference ) then
           size_epe_reference=size( epe_reference )
           call commpack(size_epe_reference, status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [31]")
           do i_ion=1,size_epe_reference
              call commpack(epe_reference(i_ion)%rc(1), 3,1, status)
              if(status.ne.0) call error_handler &
                   ("epe_send_init_to_slave : error [32]") 
              call commpack(epe_reference(i_ion)%rs(1), 3,1, status)
              if(status.ne.0) call error_handler &
                   ("epe_send_init_to_slave : error [33]") 
           end do
        end if ! if( allocated_epe_reference )
        
        call commpack(type_impurity(1), ndpt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [34]")
        call commpack(q_sh_impurity(1), ndpt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [35]")
        call commpack(q_nuc_impurity(1), ndpt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [36]")
        call commpack(imp2center(1), ndpt, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [37]")
        call commpack(r_imp(1,1),ndpt*3,1,status) 
        ! I did not find other place where r_imp communicated to slaves
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [37a]")
        DPRINT 'r_imp commpack done'

        call commpack(n_types_central_atoms_3body,status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [41]")

        if(n_types_central_atoms_3body > 0) then
           first_ind=last_ind+1
           last_ind=first_ind+n_reg-1
           if(last_ind > n_tetrahedrons) last_ind=n_tetrahedrons

           call commpack(n_tetrahedrons,status)       !!!!!!!!!!!!!!!!!AS
           call commpack(first_ind,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [37]")
           call commpack(last_ind,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [38]")
!!$           allocate(tet(last_ind-first_ind+1,5))
!!$           tet=tetra_atoms(first_ind:last_ind,1:5)
!!$           call commpack(tet(1,1),5*(last_ind-first_ind+1),1,status)
           call commpack(tetra_atoms(1,1),5*n_tetrahedrons,1,status) !!!!!!!!!!!!!!!!!!AS
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [38]")
!!$           deallocate(tet)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [39]")
           call commpack(ki(1,1,1),max_type_ions*max_type_ions*max_type_ions,1,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [40]")
           call commpack(theta_0(1,1,1),max_type_ions*max_type_ions*max_type_ions,1,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [41]")

           if(explicit_coupling) then
           call commpack(ec%ki(1,1,1),ec_max_type_ions*ec_max_type_ions*ec_max_type_ions,1,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [40]")
           call commpack(ec%theta_0(1,1,1),ec_max_type_ions*ec_max_type_ions*ec_max_type_ions,1,status)
           if(status.ne.0) call error_handler &
                ("epe_send_init_to_slave : error [42]")
           endif

        endif

        call comm_send()
     enddo

     if(n_types_central_atoms_3body > 0) then
        first_ind=1
        last_ind=first_ind+n_reg-1
     endif
   end subroutine epe_send_init_to_slave
   !*************************************************************

   !*************************************************************
   subroutine epe_send_finish_to_slave()
     !  Purpose: Tells to slaves to finish epe-calculations  
     !           Called by Master (epe-driver).
     !------------ Modules used ----------------------------------
     implicit none
     !------------ Executable code --------------------------------
!!$     print*, 'Master: => EPE-finish message'
     call comm_init_send(comm_all_other_hosts, msgtag_epe_finish_slave)
     call comm_send()
!!$     print*, 'Master: EPE-finish message, Ok' 
   end subroutine epe_send_finish_to_slave
   !*************************************************************
   
   !*************************************************************
   subroutine epe_init_slave()
     !  Purpose: Initialization of slaves for parallel 
     !           epe_calculations. Called by slave.
     !------------ Modules used ----------------------------------
        use epepar_module, only: ec_max_type_ions=>max_type_ions
     implicit none
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind) :: status, i_ion, size_reg_reference, &
          size_epe_reference,nt
     logical               :: allocated_epe_reference
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
     call communpack(n_gen_ions, status)
!!$print*,n_gen_ions,'n_gen_ions SLAVE',comm_myindex()
     if(status.ne.0) call error_handler("epe_init_slave : error [1]")
     call communpack(reg_I_n_ions, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [2]")
     call communpack(n_vacancies, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [3]")
     call communpack(n_impurities, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [2]")
     call communpack(reg_2a_n_ions, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [4]")

     allocate( r_nuc_ion(n_gen_ions,3), stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_init_slave : error [5]")
        epealloc_stat(9)=1
     allocate( r_sh_ion(n_gen_ions,3), stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_init_slave : error [6]")
     allocate( epe(n_gen_ions), stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_init_slave : error [7]")

     allocate( point_0_core_grad(reg_I_n_ions,3), stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [8]")
     allocate( point_0_shell_grad(reg_I_n_ions,3), stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [9]")
     allocate( point_1_core_grad(reg_I_n_ions,3), stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [10]")
     allocate( point_1_shell_grad(reg_I_n_ions,3), stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [11]")
     allocate( short_range_grad(reg_I_n_ions,3),  stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [12]")

     if(n_impurities.ne.0) then
     allocate( curr%imp_core_grad(n_impurities,3),curr%imp_shell_grad(n_impurities,3) &
              ,prev%imp_core_grad(n_impurities,3),prev%imp_shell_grad(n_impurities,3) &
              ,stat=epealloc_stat(15))
     ASSERT(epealloc_stat(15).eq.0)
            epealloc_stat(15)=1
endif

     call communpack(PIS, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [14]")
     call communpack(erfo, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [15]")
     call communpack(et2, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [16]")
     call communpack(error_function_parameter, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [17]")
     do i_ion=1, n_gen_ions
        call communpack(epe(i_ion)%r(1), 3, 1, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [18]")
        call communpack(epe(i_ion)%k, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [19]")
        call communpack(epe(i_ion)%ec, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [19]")
        call communpack(epe(i_ion)%m, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [19]")
     end do
     call communpack(radius_long_interact, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [20]")
     call communpack(max_type_ions, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [23]")
     call communpack(fixed_reg2a_relcontribs, status)
     if(status.ne.0) call error_handler &
        ("epe_init_slave : error fixed_reg2a_relcontribs")
     call communpack(explicit_coupling, status)
     if(status.ne.0) call error_handler &
        ("epe_init_slave : error explicit_coupling")
     if(explicit_coupling) then
        call communpack(ec_max_type_ions, status)
     if(status.ne.0) call error_handler &
        ("epe_init_slave : error ec_max_type_ions")
        call communpack(explicit_coupling_type(1),n_impurities,1,status)
     if(status.ne.0) call error_handler &
        ("epe_init_slave : error explicit_coupling_type")        
     endif
     call communpack(q_shell(1), ndt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [24]")
     call communpack(q_nuclear(1), ndt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [25]")
     call communpack(q_ion(1), ndt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [26]")
     call communpack(pk(1), ndt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [27]")
        if(explicit_coupling) then
        nt=max(max_type_ions,ec_max_type_ions)
        else
        nt=max_type_ions
        endif
     allocate(host%sr1(nt,nt,nt+1), host%sr2(nt,nt,nt+1), &                !!!!!!!!!!!!!!!!!!AS
          host%r0(nt,nt,nt+1), host%k(nt,nt,nt+1), &
          host%r1(nt,nt,nt+1), host%k1(nt,nt,nt+1), &
          host%ro(nt,nt,nt+1), host%b(nt,nt,nt+1), host%c(nt,nt,nt+1), &
          host%d(nt,nt,nt+1), &
         stat=epealloc_stat(4))
     if(epealloc_stat(4).ne.0) call error_handler("allocate host  failed")     
        epealloc_stat(4)=1
     call communpack(host%ro(1,1,1),nt*nt*(nt+1) , 1, status)         !!!!!!!!!!!!!!!AS
     if(status.ne.0) call error_handler("epe_init_slave : error [28]")
     call communpack(host%b(1,1,1),nt*nt*(nt+1) , 1, status)         !!!!!!!!!!!!!!!AS
     if(status.ne.0) call error_handler("epe_init_slave : error [29]")
     call communpack(host%c(1,1,1),nt*nt*(nt+1) , 1, status)         !!!!!!!!!!!!!!!AS
     if(status.ne.0) call error_handler("epe_init_slave : error [30]")
     call communpack(host%d(1,1,1),nt*nt*(nt+1) , 1, status)         !!!!!!!!!!!!!!!AS

     call communpack(host%sr1(1,1,1),nt*nt*(nt+1),1,status)         !!!!!!!!!!!!!!!AS
     call communpack(host%sr2(1,1,1),nt*nt*(nt+1),1,status)

     call communpack(host%k(1,1,1),nt*nt*(nt+1),1,status)
     call communpack(host%r0(1,1,1),nt*nt*(nt+1),1,status)
     call communpack(host%k1(1,1,1),nt*nt*(nt+1),1,status)
     call communpack(host%r1(1,1,1),nt*nt*(nt+1),1,status)
     if(status.ne.0) call error_handler("epe_init_slave : error [31]")

     call communpack(n_bs_points, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [32]")
     call communpack(GSTR(1,1), 3*ndngv, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [33]")
     call communpack(RSIN(1), ndngv, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [34]")
     call communpack(RCOS(1), ndngv, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [35]")
     call communpack(madc%emd(1), nmadt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [36]")
     call communpack(use_epe_reference, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [38]")   
     
     if(use_epe_reference) then
        call communpack(size_reg_reference, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [39]")
        if( .not. allocated(reg_reference) ) & 
             allocate( reg_reference(size_reg_reference),stat=status )
        if(status.ne.0) call error_handler(" allocate reg_reference failed")
        do i_ion=1, size_reg_reference
           call communpack(reg_reference(i_ion)%rc(1), 3, 1, status)
           if(status.ne.0) call error_handler("epe_init_slave : error [40]")
           
           call communpack(reg_reference(i_ion)%rs(1), 3, 1, status)
           if(status.ne.0) call error_handler("epe_init_slave : error [41]")
        end do
     end if

     call communpack(allocated_epe_reference, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [42]")

     if( allocated_epe_reference ) then
        call communpack(size_epe_reference, status)
        if(status.ne.0) call error_handler("epe_init_slave : error [43]")
        if(.not.allocated(epe_reference)) then
           allocate( epe_reference(size_epe_reference), stat= status)
           if(status.ne.0) then
              call error_handler("allocate epe_reference failed")
!!$           else
!!$              call   write_to_output_units(" epe_reference on slave allocated")
           end if
        end if

        do i_ion=1,size_epe_reference
           call communpack(epe_reference(i_ion)%rc(1), 3,1, status)
           if(status.ne.0) call error_handler("epe_init_slave : error [44]")
           call communpack(epe_reference(i_ion)%rs(1), 3,1, status)
           if(status.ne.0) call error_handler("epe_init_slave : error [45]")
        end do
     end if! if( allocated_epe_reference )

     call communpack(type_impurity(1), ndpt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [46]")   
     call communpack(q_sh_impurity(1), ndpt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [47]")
     call communpack(q_nuc_impurity(1), ndpt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [48]") 
     call communpack(imp2center(1), ndpt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [49]")
     call communpack(r_imp(1,1),3*ndpt, 1, status)
     if(status.ne.0) call error_handler("epe_init_slave : error [49a]")
     DPRINT 'r_imp communpack done'

     call communpack(n_types_central_atoms_3body,status)
     if(status.ne.0) call error_handler("epe_init_slave : error [54]")

     if(n_types_central_atoms_3body > 0) then

        call communpack(n_tetrahedrons,status)       !!!!!!!!!!!!!!!!!AS
        call communpack(first_ind,status)
        if(status.ne.0) call error_handler &
             ("epe_init_slave : error [49]")
        call communpack(last_ind,status)
        if(status.ne.0) call error_handler &
             ("epe_init_slave : error [50]")
!!$        allocate(tetra_atoms(last_ind-first_ind+1,5))
        allocate(tetra_atoms(n_tetrahedrons,5))          !!!!!!!!!!!!!!!AS
!!$        call communpack(tetra_atoms(1,1),5*(last_ind-first_ind+1),1,status)
        call communpack(tetra_atoms(1,1),5*n_tetrahedrons,1,status)    !!!!!!!!!!!!!AS
        if(status.ne.0) call error_handler &
             ("epe_init_slave : error [51]")
        allocate(ki(max_type_ions,max_type_ions,max_type_ions))
        call communpack(ki(1,1,1),max_type_ions*max_type_ions*max_type_ions,1,status)
        if(status.ne.0) call error_handler &
             ("epe_init_slave : error [53]")

        allocate(theta_0(max_type_ions,max_type_ions,max_type_ions),stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
               epealloc_stat(13)=1
        call communpack(theta_0(1,1,1),max_type_ions*max_type_ions*max_type_ions,1,status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [54]")
       
           if(explicit_coupling) then
        allocate(ec%ki(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions))
        call communpack(ec%ki(1,1,1),ec_max_type_ions*ec_max_type_ions*ec_max_type_ions,1,status)
        if(status.ne.0) call error_handler &
             ("epe_init_slave : error [53]")
        allocate(ec%theta_0(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions),stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
               epealloc_stat(13)=1
        call communpack(ec%theta_0(1,1,1),ec_max_type_ions*ec_max_type_ions*ec_max_type_ions,1,status)
        if(status.ne.0) call error_handler &
             ("epe_send_init_to_slave : error [54]")
           endif
     endif
     
   end subroutine epe_init_slave
   !*************************************************************

   !*************************************************************
   subroutine epe_finish_slave()
     !  Purpose: Finish epe on slaves.
     !           Called by slave.
     !------------ Modules used ----------------------------------
     implicit none
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind) :: status,i_alloc
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
     deallocate( r_nuc_ion, stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_finish_slave : error [1]")
     deallocate( r_sh_ion, stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_finish_slave : error [2]")
     deallocate( epe, stat=epealloc_stat(9))
     if(epealloc_stat(9).ne.0) call error_handler("epe_finish_slave : error [3]")

     deallocate( point_0_core_grad, stat=status)
     if(status.ne.0) call error_handler("epe_finish_slave : error [4]")
     deallocate( point_0_shell_grad, stat=status)
     if(status.ne.0) call error_handler("epe_finish_slave : error [5]")
     deallocate( point_1_core_grad, stat=status)
     if(status.ne.0) call error_handler("epe_finish_slave : error [6]")
     deallocate( point_1_shell_grad, stat=status)
     if(status.ne.0) call error_handler("epe_finish_slave : error [7]")
     deallocate( short_range_grad,  stat=status)
     if(status.ne.0) call error_handler("epe_init_slave : error [12]")

     if(n_impurities.ne.0) then
     deallocate( curr%imp_core_grad,curr%imp_shell_grad, &
                 prev%imp_core_grad,prev%imp_shell_grad,stat=epealloc_stat(15))
     ASSERT(epealloc_stat(15).eq.0)
     endif

     if(associated(host%sr1)) then
        deallocate(host%sr1,host%sr2,stat=epealloc_stat(1))
        if(epealloc_stat(1).ne.0)  call error_handler(" deallocate host%sr1 failed")
     end if
     if(associated(host%b)) then
        deallocate(host%b,host%ro,host%c,host%d,host%k,host%r0,host%k1,host%r1, &
                   stat=epealloc_stat(2))
      if(epealloc_stat(2).ne.0) call error_handler(" deallocate host%b failed")
     endif

     if(allocated(reg_reference)) then
        deallocate(reg_reference,stat=status)
        if(status.ne.0) call error_handler(" deallocate reg_reference failed")
     end if
     if(allocated(epe_reference)) then
        deallocate(epe_reference,stat=status)
        if(status.ne.0) call error_handler(" deallocate epe_reference failed")
     end if
     if(allocated(tetra_atoms)) then
        deallocate(tetra_atoms,stat=status)
        if(status.ne.0) call error_handler(" deallocate tetra_atoms failed")
     end if
     if(allocated(ki)) then
        deallocate(ki,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
     end if
     if(allocated(theta_0)) then
        deallocate(theta_0,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
     end if

     if(associated(ec%ki)) then
        deallocate(ec%ki,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
     end if
     if(associated(ec%theta_0)) then
        deallocate(ec%theta_0,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
     end if

     
   
    do i_alloc=1,size(epealloc_stat)
     if(epealloc_stat(i_alloc).ne.0) print*,'epealloc_stat not dealloc',i_alloc
    enddo 
   end subroutine epe_finish_slave

   subroutine main_epe
     real(kind=r8_kind)               :: var_contrib,e_lattice
     !integer(kind=i4_kind), parameter :: msg_init=0,     & 
     !                                    msg_cons=1,     &
     !                                    msg_defects=2,  &
     !                                    msg_consdef=3,  &
     !                                    msg_def_fin=4
     real(kind=r8_kind) :: previous_t,next_t,sum_t
     logical :: stop_opt=.false.
 
     integer(kind=i4_kind):: reg2a_unit,counter

     if(lpotcalc) then
        print*,'THE CALCULATION OF AN ELECTROSTATIC POTENTIAL 2-D MAP'
        write(output_epe,*) 'THE CALCULATION OF AN ELECTROSTATIC POTENTIAL 2-D MAP'
        call get_plane_grid()
        call calc_pot()
        print*,'THE CALCULATION OF AN ELECTROSTATIC POTENTIAL 2-D MAP HAS BEEN DONE SUCCESSFULLY'
        write(output_epe,*) 'THE CALCULATION OF AN ELECTROSTATIC POTENTIAL 2-D MAP HAS BEEN DONE SUCCESSFULLY'
        return
     endif

     if(periodic_optimization) then

        call periodic_optim_actions() !=>gxcell

        if(unit_cell_output) call ucell2xyz()   
        write(output_epe,*) '------------PERIODIC OPTIMIZATION---------------'
        print*, '------------PERIODIC OPTIMIZATION---------------'
        stop_opt=.false.
        previous_t=zero
        call time_period(previous_t,next_t,sum_t)
        do
           if(stop_opt) exit
           IF(epeit_count.ge.N_ITERATIONS) exit ! else continue optimization
           epeit_count=epeit_count+1

           call gopt_main_action
           if(unit_cell_output .and. mod(epeit_count,output_step)==0) call ucell2xyz()  

           call time_period(previous_t,next_t,sum_t)
           write(*,'(a37,i4,a9,f15.8)') '<<<<< next_step_to_minimum_gopt >>>>>',epeit_count,', time = ',next_t
           write(output_epe,*) '==========================================='
           write(output_epe,'(a37,i4,a9,f15.8)') '<<<<< NEXT STEP TO MINIMUM GOPT >>>>>',epeit_count,', time = ',next_t
           CALL next_step_BFGS_gopt(stop_opt)
           write(*,*) '==========================================================='
           if(n_vacancies.ne.0) then
              regI_previous(:)%r_core(1)=epe(1:reg_I_n_ions)%r(1)
              regI_previous(:)%r_core(2)=epe(1:reg_I_n_ions)%r(2)
              regI_previous(:)%r_core(3)=epe(1:reg_I_n_ions)%r(3)
              regI_previous(:)%r_shell(1)=epe(1:reg_I_n_ions)%r(1)
              regI_previous(:)%r_shell(2)=epe(1:reg_I_n_ions)%r(2)
              regI_previous(:)%r_shell(3)=epe(1:reg_I_n_ions)%r(3)
           end if
        enddo

        if(unit_cell_output .and. mod(epeit_count,output_step) == 0) call ucell2xyz()
        if(region_I_output) call regionI2xyz()

        print*,'-----------periodic optimization finished------------------'
        write(*,'(a25,f15.8)') 'Full optimization time = ',sum_t
        write(output_epe,*) '---------------periodic optimization finished-------------------'
        write(output_epe,'(a25,f15.8)') 'Full optimization time = ',sum_t
        return
     endif ! periodic_optimization

     DPRINT 'use_epe_reference in main_epe ', use_epe_reference
     use_epe_pgdata=use_pgdata.and.ex_pgdata

 !   print*,'basic_action', basic_action
     bas_ac: select case(basic_action)
     case(0:1)

 !     print*, 'pg_interfaced_mode', pg_interfaced_mode
       pg_int1: selectcase(pg_interfaced_mode) 
        case(0:1)
           if( comm_parallel() )  call epe_send_shells_and_cores(msg_init) ! + allocation
 !         print*, 'init_lat_gradients'
           call init_lat_gradients() ! in particular in initial epe call
 !         print*, 'start_slave_defect'
           call start_slave_defect ( comm_parallel(), msg_defects ) 
           call start_timer(timer_defect_ext)
           DPRINT 'before defect_contributions s c',defect_energy_short,defect_energy_coul
 !         print*, 'defect_contributions'
           call defect_contributions()
 !         print*,'defect_contributions done 1',energy_shortrange(.false.)

           DPRINT 'after  defect_contributions s c',defect_energy_short,defect_energy_coul
           call stop_timer(timer_defect_ext)
           call init_e_calc !+ save initial var contrib
           call get_epe_energies(etot_epe_corrected=e_lattice)
 !         print*,'get_epe_energies 1 done'
           if(pg_interfaced_mode.eq.1) &
                call get_epe_energies(var_epe_contrib=var_contrib)
 !         print*,'get_epe_energies 2 done'

           if(pg_interfaced_mode.ne.0)    n_iterations=1

 !         print*,'n_iterations.eq.0',n_iterations
           if(n_iterations.eq.0) then
              if(qm_interfaced_mode)  then
                 imp_conf=>impurities_old_conf
                 call epe2pg(prev%imp_core_grad,prev%imp_shell_grad)
              endif
              if(write_configuration) call write_config
 !            print*,'write_config done 2'
              pg_interfaced_mode=1
              return
           endif ! n_iterations.eq.0


 !print*,
           !** continue initialization
           if(qm_interfaced_mode) then ! no cluster 
              if(read_configuration.and.not_first_step) then ! second and futher steps
                 call swap_regI_data 
                 DELTA=DEL
                 DPRINT 'continue initialization, next_step_to_minimum'
!!$                 if (basic_action.ne.0) CALL next_step_to_minimum()  !!!!!!!!!!!!!!?????????????

              elseif(basic_action.ne.0) then
!!!             CALL reg2a_ml_displacements
                 if(.not.read_configuration) then
!!$                    print*, '** make first step mott-litleton displacements',&
!!$                         ' init_regI_displacements'
!!$                    call init_regI_displacements  !!!!!!!!!!!!!!????????????
                 endif
              endif
           endif   !interfaced_mode/else

           pg_int2: selectcase(pg_interfaced_mode)
           case(1)
              if(reg_2a_treated.and..not.ml_displacements_fixed) then   
                 call reg2a_ml_displacements
                 reg2a_unit=openget_iounit( &
                      trim(epe_input_dir)//'/reg2a_out',status='unknown')

                 do counter=reg_I_n_ions+1,reg_2a_n_ions
                    write(reg2a_unit,*) epe(counter)%s
                    write(reg2a_unit,*) epe(counter)%c
                 end do
                 call returnclose_iounit(reg2a_unit)
                 ml_displacements_fixed=.true.
              end if
              
               if(basic_action.ne.0) then
!!$              print*,'swap and next_step_to_minimum'
                 call swap_regI_data 
                 CALL next_step_to_minimum()  !!!!!!!!!!!!!!!!!?????????????
               endif
!!$              print *, 'switching to pg_interfaced_mode=2 stop_opt ',stop_opt
              imp_conf=>impurities_old_conf
              call epe2pg(prev%imp_core_grad,prev%imp_shell_grad)
              pg_interfaced_mode=2
              stop_opt=.false.
              return
           endselect pg_int2
        endselect pg_int1

        if(pg_interfaced_mode.eq.2.and..not.stop_opt) then
           if( comm_parallel() )  call epe_send_shells_and_cores(msg_cons) 
           call cons_latt_gradients()
           call start_slave_defect( comm_parallel(), msg_consdef )
           call start_timer(timer_defect_ext)
           call cons_defect_contrubutions()
           call stop_timer(timer_defect_ext) 
           if(basic_action.ne.0) then
              if(epeit_count.eq.0) then
                 write(output_epe,*) ' -------------PG epe cycles------------'
                 print*,' -------------PG epe cycles------------ '
                 call write_to_trace_unit("-----------------------------------")
                 call write_to_trace_unit("epe_loop   epe_energy      gnorm")
                 call write_to_trace_unit("-----------------------------------")
                 CALL next_step_to_minimum() !!!!!!!!!!!!???????????????
                 epeit_count=epeit_count+1
              else
                 CALL next_step_DFPPG(stop_opt)
                 epeit_count=epeit_count+1
                 write(output_epe,*) ' --------------------------------------'
                 print*,' -------------------------------------- '
              endif
           else
              dg_convergence_reached=.true.
              stop_opt=.true.
                 write(output_epe,*) ' --------------- END ------------------'
                 print*,' --------------- END ------------------ '
           endif

           if(read_configuration)  then
!!$     print*,'check_for_changes for pg_interfaced mode 2'
        call check_for_changes()
        endif
           call epe2pg(prev%imp_core_grad,prev%imp_shell_grad) !->return


        elseif(stop_opt) then
!!$           print* ,'no actions for flag stop_opt on'
           dg_convergence_reached=.true.
           if(read_configuration) call check_for_changes()
           call epe2pg(prev%imp_core_grad,prev%imp_shell_grad) !->return

        else
           call write_to_trace_unit("======== Start EPE relaxation ==========")
           call write_to_trace_unit("-------------------------------------------")
           call write_to_trace_unit("epe_loop   epe_energy      gnorm      time ")
           call write_to_trace_unit("-------------------------------------------")
           print*,'==================Start EPE relaxation======================='
           write(output_epe,*) '===================Start EPE relaxation======================='
                 call swap_regI_data 
                if (basic_action.ne.0) CALL next_step_to_minimum() !!!!!!!!!!!!!!!!!??????????????
           stop_opt=.false.
           previous_t=zero
           call time_period(previous_t,next_t,sum_t)

           do
              if(stop_opt) exit 
              if( comm_parallel() )  call epe_send_shells_and_cores(msg_cons) ! + allocation 
              IF(epeit_count.ge.N_ITERATIONS) exit ! else continue optimization
              epeit_count=epeit_count+1
              call cons_latt_gradients()
              call start_slave_defect( comm_parallel(), msg_consdef )
              call start_timer(timer_defect_ext)
              call cons_defect_contrubutions()
              call stop_timer(timer_defect_ext)
              call time_period(previous_t,next_t,sum_t)
              write(*,'(1x,a7,f15.8)') 'Time = ',next_t
              write(output_epe,'(1x,a7,f15.8)') 'Time = ',next_t
              if(basic_action.ne.0) CALL next_step_DFP(stop_opt,next_t)

              print*,'============================================================='
              write(output_epe,*) '============================================================='

              if(.not.(qm_interfaced_mode.or.operations_make_reg_reference)) cycle 
  !             print*,'call epe2pg 1'
              if(N_IMPURITIES.gt.0) then
                   call epe2pg(prev%imp_core_grad,prev%imp_shell_grad)
               else
                   call epe2pg_epe()
               endif
   !             print*,'done'

           enddo
           call write_to_trace_unit("========= End of EPE relaxation ===========")
        endif
        write(output_epe,'(a25,f15.8)') 'Full optimization time = ',sum_t
        write(*,'(a25,f15.8)') 'Full optimization time = ',sum_t
        if(print_gradients) call print_gradients_r()
        if(region_I_output) call regionIs2xyz
        not_first_step=.true. 

     case(2) ! basic_action
        if( comm_parallel() )  call epe_send_shells_and_cores(msg_defects)
        call start_slave_defect ( comm_parallel(), msg_defects ) 
        call start_timer(timer_defect_ext)
        call defect_contributions()

        call stop_timer(timer_defect_ext) 
        call init_e_calc
        call init_imp_displacements
        if(ml_displacements) CALL reg2a_ml_displacements
        do
           IF(epeit_count.ge.N_ITERATIONS) exit
           epeit_count=epeit_count+1
           if( comm_parallel() )  call epe_send_shells_and_cores(msg_defects)
           call start_slave_defect( comm_parallel(), msg_consdef )
           call start_timer(timer_defect_ext)
           call cons_defect_contrubutions()
           call stop_timer(timer_defect_ext)
           CALL next_step_to_minimum() !!!!!!!!!!!!!!!!!???????????
!!!    CALL reg2a_ml_displacements 
        enddo
        if(print_gradients) call print_gradients_imp()
        call output_r_imp
     endselect bas_ac
  !    print*,'main_epe_done'

!-------------------------------------------------
   contains

     subroutine time_period(previous,next,sum)

       real(kind=r8_kind), intent(inout) ::previous,next,sum
       real(kind=r8_kind) :: previous1
!!$       integer(kind=i4_kind) :: count,count_rate,count_max
       character(len=10) :: date,time,zone
       integer(kind=i4_kind) :: values(8)
       
       previous1=previous
!!$       call system_clock(count,count_rate,count_max)
!!$       next=(count*1.0_r8_kind/count_rate)
       call date_and_time(date,time,zone,values)
       next=values(8)*0.001_r8_kind+values(7)*1.0_r8_kind+ &
            values(6)*60.0_r8_kind+values(5)*3600.0_r8_kind+ &
            values(3)*86400.0_r8_kind
       if(previous == zero) then
          previous=next
          next=zero
          sum=zero
       else
          previous=next
          next=next-previous1
       endif
       sum=sum+next

     end subroutine time_period
!-----------------------------------------------------------

!-----------------------------------------------------------
     subroutine print_gradients_r

       integer(kind=i4_kind) :: i
       
       write(output_epe,*) '============================================================================'
       write(output_epe,*) '                                     GRADIENTS'
       write(output_epe,*) '                  SHELL                                 CORE'
       do i=n_vacancies+1,reg_I_n_ions
          write(output_epe,'(i3,6(1x,f12.8))')i,point_1_shell_grad(i,:),point_1_core_grad(i,:)
       enddo
       write(output_epe,*) '============================================================================'

     end subroutine print_gradients_r
!-----------------------------------------------------------

!-----------------------------------------------------------
     subroutine print_gradients_imp

       integer(kind=i4_kind) :: i
       
       write(output_epe,*) '============================================================================'
       write(output_epe,*) '                             IMPURITY  GRADIENTS'
       write(output_epe,*) '                  SHELL                                 CORE'
       do i=1,N_IMPURITIES
          write(output_epe,'(i3,6(1x,f12.8))')i,curr%imp_shell_grad(i,:),curr%imp_core_grad(i,:)
       enddo
       write(output_epe,*) '============================================================================'

     end subroutine print_gradients_imp
!-----------------------------------------------------------

!-----------------------------------------------------------
     subroutine swap_regI_data
       IF(N_IMPURITIES.ne.0) then
          impurities_new_conf=impurities_old_conf !!!!
          curr%imp_shell_grad(1:N_IMPURITIES,1:3)=prev%imp_shell_grad(1:N_IMPURITIES,1:3)
          curr%imp_core_grad(1:N_IMPURITIES,1:3)=prev%imp_core_grad(1:N_IMPURITIES,1:3)
       endif ! N_IMPURITIES.ne.0
       regI_previous(:)%r_shell(1)=R_SH_ION(1:reg_I_n_ions,1)
       regI_previous(:)%r_shell(2)=R_SH_ION(1:reg_I_n_ions,2)
       regI_previous(:)%r_shell(3)=R_SH_ION(1:reg_I_n_ions,3)
       regI_previous(:)%r_core(1)=R_NUC_ION(1:reg_I_n_ions,1)
       regI_previous(:)%r_core(2)=R_NUC_ION(1:reg_I_n_ions,2)
       regI_previous(:)%r_core(3)=R_NUC_ION(1:reg_I_n_ions,3)
       point_1_core_grad(1:reg_I_n_ions,1:3)=point_0_core_grad(1:reg_I_n_ions,1:3)
       point_1_shell_grad(1:reg_I_n_ions,1:3)=point_0_shell_grad(1:reg_I_n_ions,1:3)
     end subroutine swap_regI_data
!----------------------------------------------------

!----------------------------------------------------
     subroutine init_imp_displacements
       R_SH_IMPO(:N_IMPURITIES,:)=R_SH_IMP(:N_IMPURITIES,:)
       R_NUC_IMPO(:N_IMPURITIES,:)=R_NUC_IMP(:N_IMPURITIES,:)
       do i=1,n_impurities
          R_SH_IMP(i,:)=R_SH_IMPO(i,:)-DEL*prev%imp_shell_grad(i,:)/(ABS(prev%imp_shell_grad(i,:))+0.01)
          R_NUC_IMP(i,:)=R_SH_IMP(i,:)-prev%imp_core_grad(i,:)/PK_IMPURITY(i)
       enddo ! I=1,N_IMPURITIES
     end subroutine init_imp_displacements
!----------------------------------------------------

!----------------------------------------------------
     subroutine next_step_to_minimum()
       ! **this program makes displacement of ions

       use epecom_module
       use culon_module
       use mol_module, only: dielectric_const
       real(kind=r8_kind) :: regI_displ_convergence,imp_displ_convergence
       real(kind=r8_kind) :: sx,sy,AMG,DELGRA,DRS,RCPROM(3),RSPROM(3),DELX
       real(kind=r8_kind) :: DELTAR,DPROM,epe_ml_dprime_factors(reg_I_n_ions,3),displ_vector(reg_I_n_ions,3)
       integer(kind=i4_kind) :: i,ind,IGR
       real(kind=r8_kind),dimension(reg_I_n_ions,3)::delta_g

       regI_displ_convergence=zero
       imp_displ_convergence=zero

       IF(basic_action.EQ.1)  then
          do I=n_vacancies+1,reg_I_n_ions
             epe_ml_dprime_factors(i,:)=ml_dprime_factors(epe(i)%k)
          enddo

if(relax_shells_only) then
          delta_g(n_vacancies+1:,:)=point_1_shell_grad(n_vacancies+1:,:) &
               -point_0_shell_grad(n_vacancies+1:,:) 
          displ_vector(n_vacancies+1:,:)=-0.05*(point_1_shell_grad(n_vacancies+1:,:)) &
               *(dielectric_const+2.0_r8_kind)/(3.0_r8_kind*dielectric_const)/epe_ml_dprime_factors( &
               n_vacancies+1:,:)
else
          delta_g(n_vacancies+1:,:)=point_1_shell_grad(n_vacancies+1:,:) &
               -point_0_shell_grad(n_vacancies+1:,:) &
               +point_1_core_grad(n_vacancies+1:,:)-point_0_core_grad(n_vacancies+1:,:)
          displ_vector(n_vacancies+1:,:)=-0.05*(point_1_shell_grad(n_vacancies+1:,:) & 
               +point_1_core_grad(n_vacancies+1:,:)) &
               *(dielectric_const+2.0_r8_kind)/(3.0_r8_kind*dielectric_const)/epe_ml_dprime_factors( &
               n_vacancies+1:,:)
endif
!        print*, 'next step: point_1_shell_grad ',sum(abs(point_1_shell_grad(n_vacancies+1:,:))), &
!                                                 sum(abs(point_0_shell_grad(n_vacancies+1:,:)))


          regI_previous(n_vacancies+1:)%r_core(1)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,1)
          regI_previous(n_vacancies+1:)%r_core(2)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,2)
          regI_previous(n_vacancies+1:)%r_core(3)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,3)
          regI_previous(n_vacancies+1:)%r_shell(1)=R_SH_ION(n_vacancies+1:reg_I_n_ions,1)
          regI_previous(n_vacancies+1:)%r_shell(2)=R_SH_ION(n_vacancies+1:reg_I_n_ions,2)
          regI_previous(n_vacancies+1:)%r_shell(3)=R_SH_ION(n_vacancies+1:reg_I_n_ions,3)

          point_0_core_grad(n_vacancies+1:reg_I_n_ions,:)=point_1_core_grad(n_vacancies+1:reg_I_n_ions,:)
          point_0_shell_grad(n_vacancies+1:reg_I_n_ions,:)=point_1_shell_grad(n_vacancies+1:reg_I_n_ions,:)

         if(relax_shells_only) then
#if 1
             R_SH_ION(I,:)=R_SH_ION(I,:)+displ_vector(I,:)
#endif
         else
          DO I=n_vacancies+1,reg_I_n_ions
             RCPROM=R_NUC_ION(I,:)-R_SH_ION(I,:) 
             RSPROM=regI_previous(i)%r_shell(:)+displ_vector(I,:)
             R_SH_ION(I,:)=RSPROM+0.5_r8_kind*point_1_core_grad(I,:)/PK(epe(I)%k)
             R_NUC_ION(I,:)=RSPROM-0.5_r8_kind*point_1_core_grad(I,:)/PK(epe(I)%k)+RCPROM
          enddo
         endif

          if(option_c3_symm) then
             call symmetrize_forces(R_sh_ION(1:reg_I_n_ions,:))  !!!
             call symmetrize_forces(R_nuc_ION(1:reg_I_n_ions,:))  !!!
          endif
       endif ! basic_action.eq.1

! **minimization by lattice coordinates
       IF(basic_action.eq.2.and.N_IMPURITIES.ne.0) then
          ! **get new coordinates of impurities (R_SH_IMP,R_NUC_IMP)
          DO IND=1,3
             DO I=1,N_IMPURITIES
                AMG=0.
                DELGRA =curr%imp_shell_grad(I,IND)-prev%imp_shell_grad (I,IND)+ &
                        curr%imp_core_grad(I,IND)-prev%imp_core_grad (I,IND)
                IF(ABS(DELGRA ).GE.1.0E-5)AMG=(R_SH_IMP(I,IND)-R_SH_IMPO(I,IND))/DELGRA
                DRS=-AMG*(curr%imp_shell_grad  (I,IND)+curr%imp_core_grad  (I,IND) )
                IF(AMG.LE.0.0)  DRS=-DRS/(ABS(DRS)+0.01)*0.01*scaling_factor
                IF(ABS(DRS).GT.imp_displ_convergence)  imp_displ_convergence=ABS(DRS)
                IF(ABS(DRS).LT.DELTA ) then
                   DELTAR=DELTA/ABS(DRS)
                   IF(DELTAR.LE.1.0)  DRS=DRS*DELTAR
                endif
                R_NUC_IMPO(i,ind)=R_NUC_IMP(i,ind)
                R_SH_IMPO(i,ind)=R_SH_IMP(i,ind)
                prev%imp_core_grad(i,ind)=curr%imp_core_grad(i,ind)
                prev%imp_shell_grad(i,ind)=curr%imp_shell_grad(i,ind)
                ! **RCPROM - displaysment of nucleus with respect to shell
                RCPROM=R_NUC_IMP(I,IND)-R_SH_IMP(I,IND)
                ! **RSPROM - intermediate value R_SH_IMP
                RSPROM(ind)=   R_SH_IMPO(I,IND)+DRS
                IGR=-1
                IF(ABS(curr%imp_core_grad(I,IND)).GT.ABS(curr%imp_shell_grad(I,IND))) IGR=+1
                DPROM=curr%imp_core_grad(I,IND)
                IF(IGR.EQ.+1) DPROM=curr%imp_shell_grad(I,IND)
                ! **DELX - additional displacement of nucleus with respect to shell
                DELX=+DPROM/PK_IMPURITY(I)
                R_SH_IMP(I,IND)=RSPROM(ind)-IGR*DELX/2.
                R_NUC_IMP(I,IND)=RSPROM(ind)+IGR*DELX/2.+RCPROM(ind)
             enddo !i=1,n_impurities
          enddo !IND=1,3
          ! **minimization by impurity coordinates
          ! **done get new coordinates of impurities (R_SH_IMP,R_NUC_IMP)
       endif !basic_action.ne.1
       SX=regI_displ_convergence/scaling_factor
       SY=imp_displ_convergence/scaling_factor
     END subroutine next_step_to_minimum
!------------------------------------------------------

!--------------------------------------------------------------
     subroutine next_step_DFP(stop_opt,time_cycle)

       use epecom_module
       use culon_module

       real(kind=r8_kind) :: time_cycle
       logical :: stop_opt
       real(kind=r8_kind), allocatable :: A_m(:,:)
       real(kind=r8_kind), allocatable, dimension(:) :: g_dim,g0_dim, &
            rs_dim,rs0_dim,rn_dim,rn0_dim,dsn, &
            step,g_step,store_v0,store_v1,r_step
       real(kind=r8_kind) :: store,delta_g,r1,r2,r3,rm,rm0,en1,en2,en3,enm
       real(kind=r8_kind) :: alpha,aa,bb,rnorm,cosm,gnorm

       integer(kind=i4_kind) :: i,j,k,icount1,status
       logical :: resett
       character(len=44) :: trace_message


       allocate(g_dim(3*(reg_I_n_ions-n_vacancies)), &
            g0_dim(3*(reg_I_n_ions-n_vacancies)), &
            rs_dim(3*(reg_I_n_ions-n_vacancies)), & 
            rs0_dim(3*(reg_I_n_ions-n_vacancies)), &
            rn_dim(3*(reg_I_n_ions-n_vacancies)), &
            rn0_dim(3*(reg_I_n_ions-n_vacancies)), &
            dsn(3*(reg_I_n_ions-n_vacancies)), &
            step(3*(reg_I_n_ions-n_vacancies)),stat=status)
       if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [1]") 

       resett=.false.
       icount1=0

1      continue

!        print*, 'next step DFP: point_1_shell_grad ',sum(abs(point_1_shell_grad(n_vacancies+1:,:))), &
!                                                 sum(abs(point_0_shell_grad(n_vacancies+1:,:)))
       ! swap data to one dimentional store
if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             g_dim(k)=point_1_shell_grad(i,j)
             g0_dim(k)=point_0_shell_grad(i,j)
             rs_dim(k)=r_sh_ion(i,j)
             rs0_dim(k)=regI_previous(i)%r_shell(j)
             rn_dim(k)=r_nuc_ion(i,j)
             rn0_dim(k)=regI_previous(i)%r_core(j)
          enddo
       enddo
else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             g_dim(k)=point_1_shell_grad(i,j)+point_1_core_grad(i,j)
             g0_dim(k)=point_0_shell_grad(i,j)+point_0_core_grad(i,j)
             rs_dim(k)=r_sh_ion(i,j)
             rs0_dim(k)=regI_previous(i)%r_shell(j)
             rn_dim(k)=r_nuc_ion(i,j)
             rn0_dim(k)=regI_previous(i)%r_core(j)
          enddo
       enddo
endif

       ! calculate delta_g
       delta_g=zero
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             if(abs(point_1_shell_grad(i,j)) > delta_g) &
                  delta_g=abs(point_1_shell_grad(i,j))
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             if(abs(point_1_shell_grad(i,j)+point_1_core_grad(i,j)) > delta_g) &
                  delta_g=abs(point_1_shell_grad(i,j)+point_1_core_grad(i,j))
          enddo
       enddo
      endif

       ! control direction of search
       rnorm=sqrt(dot_product(rs_dim-rs0_dim,rs_dim-rs0_dim))
       gnorm=sqrt(dot_product(g_dim,g_dim))
       cosm=abs(dot_product(rs_dim-rs0_dim,g_dim))/(rnorm*gnorm)
       if(.not.resett) then
          write (output_epe,*) 'Gnorm= ',gnorm/(3*(reg_I_n_ions-n_vacancies))
          write (*,*) 'DFP: Gnorm= ',gnorm/(3*(reg_I_n_ions-n_vacancies))
       endif
       write(trace_message,'(i5,5x,f11.7,3x,1p,e11.5e2,0p,3x,f5.2)') epeit_count, &
            etot_epe/eau_ev,gnorm/(3*(reg_I_n_ions-n_vacancies)),time_cycle
       call write_to_trace_unit(trace_message)

       ! check convirgence by delta_g criterium
       if (epeit_count > 1 .and. .not. resett) then
          if (delta_g <= abs_g) then
             write (output_epe,*) '================================='
             write (output_epe,*) '     Minimun has been found'
             write (output_epe,*) '================================='
             stop_opt=.true.
             dg_convergence_reached=.true.
!!$             print*, '     Minimun has been found'
             if(allocated(H_m)) deallocate(H_m, &
                  g_dim, g0_dim, &
                  rs_dim, rs0_dim, &
                  rn_dim, rn0_dim, &
                  step,dsn,stat=status)
             if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [2]")
             return
          endif
       endif
       
       if(cosm < 0.01_r8_kind .and. epeit_count > 1) then
          resett=.true.
       endif

    !initialize or reset Hessian
       alpha=0.5_r8_kind
       if ( epeit_count == 1 .or. &
            (mod(epeit_count-1,n_hess_update_cycles) == 0 &
            .and. epeit_count > 1) .or.  resett ) then
          if(.not. allocated(H_m)) then
             allocate(H_m(3*(reg_I_n_ions-n_vacancies), &
                  3*(reg_I_n_ions-n_vacancies)),stat=status)
             if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [3]")
          endif
          if(epeit_count /=1) then
             write (*,*) '     Hessian reset'
             write (output_epe,*) '================================='
             write (output_epe,*) '        Hessian reset'
             write (output_epe,*) '================================='
             else
             write (*,*) ' Hessian initialized next_step_DFP'  
          endif
          H_m=zero
          do i=1,3*(reg_I_n_ions-n_vacancies)
             H_m(i,i)=0.5_r8_kind*weight_hess
          enddo
          etot_epe_0=etot_epe
          goto 2
       endif

       !define new search direction
       if (epeit_count > 1 .and. mod(epeit_count-1,n_hess_update_cycles) /= 0 &
            .or. .not.resett) then
          allocate(A_m(3*(reg_I_n_ions-n_vacancies),3*(reg_I_n_ions-n_vacancies)), &
               store_v0(3*(reg_I_n_ions-n_vacancies)), &
               store_v1(3*(reg_I_n_ions-n_vacancies)), &
               r_step(3*(reg_I_n_ions-n_vacancies)), &
               g_step(3*(reg_I_n_ions-n_vacancies)),stat=status)
          if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [4]") 


          r_step=rs_dim-rs0_dim
          g_step=g_dim-g0_dim
!vvvvvvvvvvvvvvvvvvvvvvv
          store=dot_product(matmul(g_step,H_m),g_step)
          store_v0=matmul(H_m,g_step)
          store_v1=matmul(g_step,H_m)
          do i=1,3*(reg_I_n_ions-n_vacancies)
             do j=1,3*(reg_I_n_ions-n_vacancies)
                A_m(i,j)=store_v0(i)*store_v1(j)/store
             enddo
          enddo
          H_m=H_m-A_m
!^^^^^^^^^^^^^^^^^^^^^^^^^^^

          store=dot_product(r_step,g_step)


          do i=1,3*(reg_I_n_ions-n_vacancies)
             do j=1,3*(reg_I_n_ions-n_vacancies)
                A_m(i,j)=r_step(i)*r_step(j)/store
             enddo
          enddo
          H_m=H_m+A_m

          deallocate(A_m, &
               store_v0, store_v1, &
               r_step, g_step,stat=status)
          if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [5]") 

       endif
       ! done H_m is calculated

2      rs0_dim=rs_dim
       rn0_dim=rn_dim
       dsn=rn0_dim-rs0_dim

!..............................................
       resett=.false.

       r1=zero
       en1=etot_epe

       do
          step=-alpha*matmul(H_m,g_dim)
          r3=sqrt(dot_product(step,step))
          if (r3 <= 1.0_r8_kind) exit
          alpha=alpha*1.0_r8_kind/r3
       enddo
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k) 
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k) &
                  +0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
             rn_dim(k)=rs0_dim(k)+step(k)+dsn(k) &
                  -0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
          enddo
       enddo
      endif

       !new positions
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             r_nuc_ion(i,j)=rn_dim(k)
          enddo
       enddo
      endif
       
       ! calc energy for new positions of lattice
       if(n_types_central_atoms_3body > 0) then
          en3=energy_shortrange(.false.)+energy_coulomb(.false.)+energy_3_body()
         DPRINT 'energy_shortrange energy_coulomb energy_3_body',en3
       else
          en3=energy_shortrange(.false.)+energy_coulomb(.false.)
       endif

       ! add contributions due to vacancies and impurities
       if(n_vacancies > 0 .or. n_impurities >0) then
          call energy_vac_imp(use_epe_reference,n_vacancies,n_impurities) !1
          en3=en3+defect_energy_short+defect_energy_coul
       DPRINT 'en3=',en3
       endif

       step=-alpha*matmul(H_m,g_dim)/5.0_r8_kind
       r2=sqrt(dot_product(step,step))
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k)
          enddo
       enddo
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k)+ &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
             rn_dim(k)=rs0_dim(k)+step(k)+dsn(k)- &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
          enddo
       enddo
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             r_nuc_ion(i,j)=rn_dim(k)
          enddo
       enddo
      endif
       if(n_types_central_atoms_3body > 0) then
          en2=energy_shortrange(.false.)+energy_coulomb(.false.)+energy_3_body()
          DPRINT 'energy_shortrange energy_coulomb energy_3_body',en2
       else
          en2=energy_shortrange(.false.)+energy_coulomb(.false.)
       endif
       if(n_vacancies > 0 .or. n_impurities >0) then
          call energy_vac_imp(use_epe_reference,n_vacancies,n_impurities) !2
          en2=en2+defect_energy_short+defect_energy_coul
          DPRINT 'en2=',en2,'defect_energies',defect_energy_short,defect_energy_coul
       endif

       rm0=r2

       aa=((en1-en2)*(r1-r3)-(en1-en3)*(r1-r2))/ &
            ((r1-r3)*(r1*r1-r2*r2)-(r1-r2)*(r1*r1-r3*r3))
       if(aa > zero) then
          bb=((en1-en2)-aa*(r1*r1-r2*r2))/(r1-r2)
          rm=-bb/(2.0_r8_kind*aa)
       else
          rm=min(abs(r2),abs(r3))/5.0_r8_kind
       endif
       if(abs(rm) > 1.0_r8_kind) rm=rm0/2.0_r8_kind
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+(rm/rm0)*step(k)
          enddo
       enddo
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+(rm/rm0)*step(k)+ &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
             rn_dim(k)=rs0_dim(k)+(rm/rm0)*step(k)+dsn(k)- &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
          enddo
       enddo
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             r_nuc_ion(i,j)=rn_dim(k)
          enddo
       enddo
      endif
       if(n_types_central_atoms_3body > 0) then
          enm=energy_shortrange(.false.)+energy_coulomb(.false.)+energy_3_body()
       else
          enm=energy_shortrange(.false.)+energy_coulomb(.false.)
       endif
       if(n_vacancies > 0 .or. n_impurities >0) then
          call energy_vac_imp(use_epe_reference,n_vacancies,n_impurities) !3
          enm=enm+defect_energy_short+defect_energy_coul
       DPRINT 'enm etot_epe', enm,etot_epe
       endif

       if(enm < etot_epe) then
       else
          resett=.true.
          icount1=icount1+1
          if(icount1 ==1) then
         if(relax_shells_only) then
             do i=n_vacancies+1,reg_I_n_ions
                do j=1,3
                   k=3*(i-n_vacancies-1)+j
                   r_sh_ion(i,j)=rs0_dim(k)
                enddo
             enddo
         else
             do i=n_vacancies+1,reg_I_n_ions
                do j=1,3
                   k=3*(i-n_vacancies-1)+j
                   r_sh_ion(i,j)=rs0_dim(k)
                   r_nuc_ion(i,j)=rn0_dim(k)
                enddo
             enddo
          endif
             goto 1
          else
             write (output_epe,*) '================================='
             write (output_epe,*) '   Next step cannot be done.'
             write (output_epe,*) '================================='
             write (*,*) '     Next step cannot be done 1'
             stop_opt=.true.
             dg_convergence_reached=.true.
             do i=n_vacancies+1,reg_I_n_ions
                do j=1,3
                   k=3*(i-n_vacancies-1)+j
                   r_sh_ion(i,j)=rs0_dim(k)
                   r_nuc_ion(i,j)=rn0_dim(k)
                enddo
             enddo
             deallocate(H_m,g_dim, g0_dim, &
                  rs_dim, rs0_dim, &
                  rn_dim, rn0_dim, &
                  step,dsn,stat=status)
             if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [6]")
             return           
          endif
       endif
!...........................................

      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             regI_previous(i)%r_shell(j)=rs0_dim(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             regI_previous(i)%r_shell(j)=rs0_dim(k)
             r_nuc_ion(i,j)=rn_dim(k)
             regI_previous(i)%r_core(j)=rn0_dim(k)
          enddo
       enddo
       endif

       etot_epe_0=etot_epe
       etot_epe=enm

       point_0_shell_grad=point_1_shell_grad
       point_0_core_grad=point_1_core_grad

       if(option_c3_symm) then
          call symmetrize_forces(R_sh_ION(1:reg_I_n_ions,:))  
          call symmetrize_forces(R_nuc_ION(1:reg_I_n_ions,:))  
       endif

       deallocate(g_dim, g0_dim, &
            rs_dim, rs0_dim, &
            rn_dim, rn0_dim, &
            step, dsn,stat=status)
       if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [7]") 

     end subroutine next_step_DFP
!------------------------------------------------------
!--------------------------------------------------------------
     subroutine next_step_DFPPG(stop_opt)

       use epecom_module
       use culon_module

       logical :: stop_opt

       real(kind=r8_kind), allocatable :: A_m(:,:)
       real(kind=r8_kind), allocatable, dimension(:) :: g_dim,g0_dim, &
            rs_dim,rs0_dim,rn_dim,rn0_dim,dsn, &
            step,g_step,store_v0,store_v1,r_step
       real(kind=r8_kind) :: store,delta_g
       real(kind=r8_kind) :: alpha,rnorm,cosm,gnorm,step_length,max_step=0.01

       integer(kind=i4_kind) :: i,j,k,icount1,status
       character(len=36) :: trace_message

       if(use_lin_search) then
          if(etot_epe_0 <= etot_epe .and. n_ls < 2) return
       end if

       allocate(g_dim(3*(reg_I_n_ions-n_vacancies)), &
            g0_dim(3*(reg_I_n_ions-n_vacancies)), &
            rs_dim(3*(reg_I_n_ions-n_vacancies)), & 
            rs0_dim(3*(reg_I_n_ions-n_vacancies)), &
            rn_dim(3*(reg_I_n_ions-n_vacancies)), &
            rn0_dim(3*(reg_I_n_ions-n_vacancies)), &
            dsn(3*(reg_I_n_ions-n_vacancies)), &
            step(3*(reg_I_n_ions-n_vacancies)),stat=status)
       if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [1]") 

       icount1=0

1      continue

      ! swap data to one dimentional store
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             g_dim(k)=point_1_shell_grad(i,j)
             g0_dim(k)=point_0_shell_grad(i,j)
             rs_dim(k)=r_sh_ion(i,j)
             rs0_dim(k)=regI_previous(i)%r_shell(j)
             rn_dim(k)=r_nuc_ion(i,j)
             rn0_dim(k)=regI_previous(i)%r_core(j)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             g_dim(k)=point_1_shell_grad(i,j)+point_1_core_grad(i,j)
             g0_dim(k)=point_0_shell_grad(i,j)+point_0_core_grad(i,j)
             rs_dim(k)=r_sh_ion(i,j)
             rs0_dim(k)=regI_previous(i)%r_shell(j)
             rn_dim(k)=r_nuc_ion(i,j)
             rn0_dim(k)=regI_previous(i)%r_core(j)
          enddo
       enddo
      endif

       ! calculate delta_g
      delta_g=zero
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             if(abs(point_1_shell_grad(i,j)) > delta_g) &
                  delta_g=abs(point_1_shell_grad(i,j))
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             if(abs(point_1_shell_grad(i,j)+point_1_core_grad(i,j)) > delta_g) &
                  delta_g=abs(point_1_shell_grad(i,j)+point_1_core_grad(i,j))
          enddo
       enddo
      endif

       rnorm=sqrt(dot_product(rs_dim-rs0_dim,rs_dim-rs0_dim))
       gnorm=sqrt(dot_product(g_dim,g_dim))
       print*,cosm,'cosm',gnorm/(3*(reg_I_n_ions-n_vacancies)),'Gnorm'
       if(gnorm/(3*(reg_I_n_ions-n_vacancies)) <= 2.0e-5_r8_kind .and. use_lin_search) &
            use_lin_search=.false.
       cosm=abs(dot_product(rs_dim-rs0_dim,g_dim))/(rnorm*gnorm+1.0E-10)
       print*,cosm,'cosm',gnorm/(3*(reg_I_n_ions-n_vacancies)),'Gnorm'
       if(.not.reset) then
          write (output_epe,*) 'Gnorm= ',gnorm/(3*(reg_I_n_ions-n_vacancies))
          write (*,*) 'Gnorm= ',gnorm/(3*(reg_I_n_ions-n_vacancies))
       endif
       write(trace_message,'(i5,5x,f11.8,3x,1p,e11.5e2,0p)') epeit_count, &
            etot_epe/eau_ev,gnorm/(3*(reg_I_n_ions-n_vacancies))
       call write_to_trace_unit(trace_message)


       ! check convirgence by delta_g criterium
       if (epeit_count > 1 .and. .not. reset) then
          if (delta_g <= abs_g) then
             write (output_epe,*) '================================='
             write (output_epe,*) '     Minimun has been found'
             write (output_epe,*) '================================='
             stop_opt=.true.
       dg_convergence_reached=.true.
!!$             print*, '     Minimun has been found, dg_convergence_reached '
             if(allocated(H_m)) deallocate(H_m, &
                  g_dim, g0_dim, &
                  rs_dim, rs0_dim, &
                  rn_dim, rn0_dim, &
                  step,dsn,stat=status)
             if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [2]")
             return
          endif
       endif
       
       ! control direction of search
       if(cosm < 0.01_r8_kind .and. epeit_count > 1) then
          reset=.true.
       endif
    

    !initialize or reset Hessian
       alpha=0.5_r8_kind
       if ( epeit_count == 1 .or.(mod(epeit_count-1,n_hess_update_cycles) == 0 &
            .and. epeit_count > 1) .or.  reset ) then
          if(.not. allocated(H_m)) then
             allocate(H_m(3*(reg_I_n_ions-n_vacancies), &
                  3*(reg_I_n_ions-n_vacancies)),stat=status)
             if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [3]")
          endif
          if(epeit_count /=1) then
             write (*,*) '     Hessian reset'
             write (output_epe,*) '================================='
             write (output_epe,*) '        Hessian reset next_step_DFPPG'
             write (output_epe,*) '================================='
             else
             write (*,*) ' Hessian initialized'  
          endif
          H_m=zero
          do i=1,3*(reg_I_n_ions-n_vacancies)
             H_m(i,i)=0.1_r8_kind*weight_hess
          enddo
          goto 2
       endif

       !define new search direction
       if (epeit_count > 1 .and. mod(epeit_count-1,n_hess_update_cycles) /= 0 &
            .or. .not.reset) then
          allocate(A_m(3*(reg_I_n_ions-n_vacancies),3*(reg_I_n_ions-n_vacancies)), &
               store_v0(3*(reg_I_n_ions-n_vacancies)), &
               store_v1(3*(reg_I_n_ions-n_vacancies)), &
               r_step(3*(reg_I_n_ions-n_vacancies)), &
               g_step(3*(reg_I_n_ions-n_vacancies)),stat=status)
          if(status.ne.0) call error_handler("next_step_DFP: wrong allocation [4]") 


          r_step=rs_dim-rs0_dim
          g_step=g_dim-g0_dim

          store=dot_product(matmul(g_step,H_m),g_step)

          
          store_v0=matmul(H_m,g_step)
          store_v1=matmul(g_step,H_m)
          do i=1,3*(reg_I_n_ions-n_vacancies)
             do j=1,3*(reg_I_n_ions-n_vacancies)
                A_m(i,j)=store_v0(i)*store_v1(j)/(store+1.0e-10)
             enddo
          enddo
          H_m=H_m-A_m


          store=dot_product(r_step,g_step)


          do i=1,3*(reg_I_n_ions-n_vacancies)
             do j=1,3*(reg_I_n_ions-n_vacancies)
                A_m(i,j)=r_step(i)*r_step(j)/(store+1.0e-10)
             enddo
          enddo
          H_m=H_m+A_m

          deallocate(A_m, store_v0, store_v1, r_step, g_step,stat=status)
          if(status.ne.0) call error_handler("next_step_DFP: wrong deallocation [5]") 

       endif
       ! done H_m is calculated

2      rs0_dim=rs_dim
       rn0_dim=rn_dim
       dsn=rn0_dim-rs0_dim

!..............................................
       reset=.false.
       step=-alpha*matmul(H_m,g_dim)
       step_length=sqrt(dot_product(step,step))/(3*(reg_I_n_ions-n_vacancies))
       write(*,*) step_length,max_step , ' next_step_DFPPG:   step_length'
       if(step_length.gt.max_step) step=step*max_step/step_length

      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             rs_dim(k)=rs0_dim(k)+step(k)+ &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
             rn_dim(k)=rs0_dim(k)+step(k)+dsn(k)- &
                  0.5_r8_kind*point_1_core_grad(i,j)/pk(epe(i)%k)
          enddo
       enddo
      endif

       !new positions
      if(relax_shells_only) then
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             regI_previous(i)%r_shell(j)=rs0_dim(k)
          enddo
       enddo
      else
       do i=n_vacancies+1,reg_I_n_ions
          do j=1,3
             k=3*(i-n_vacancies-1)+j
             r_sh_ion(i,j)=rs_dim(k)
             regI_previous(i)%r_shell(j)=rs0_dim(k)
             r_nuc_ion(i,j)=rn_dim(k)
             regI_previous(i)%r_core(j)=rn0_dim(k)
          enddo
       enddo
     endif

       point_0_shell_grad=point_1_shell_grad
       point_0_core_grad=point_1_core_grad

       if(option_c3_symm) then
          call symmetrize_forces(R_sh_ION(1:reg_I_n_ions,:))  
          call symmetrize_forces(R_nuc_ION(1:reg_I_n_ions,:))  
       endif

       deallocate(g_dim, g0_dim, rs_dim, rs0_dim,  rn_dim, rn0_dim, &
            step, dsn,stat=status)
       if(status.ne.0) call error_handler("next_step_DFPG: wrong deallocation [7]") 

     end subroutine next_step_DFPPG
!------------------------------------------------------

!------------------------------------------------------
     subroutine init_regI_displacements
       real(kind=r8_kind):: DELS(3),delmax,gradmax,two=2.0_r8_kind,three=3.0_r8_kind, &
            maxfac
       maxfac=0.0_r8_kind
       delmax=0.0_r8_kind
       gradmax=0.0_r8_kind
       regI_previous(n_vacancies+1:)%r_core(1)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,1)
       regI_previous(n_vacancies+1:)%r_core(2)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,2)
       regI_previous(n_vacancies+1:)%r_core(3)=R_NUC_ION(n_vacancies+1:reg_I_n_ions,3)
       regI_previous(n_vacancies+1:)%r_shell(1)=R_SH_ION(n_vacancies+1:reg_I_n_ions,1)
       regI_previous(n_vacancies+1:)%r_shell(2)=R_SH_ION(n_vacancies+1:reg_I_n_ions,2)
       regI_previous(n_vacancies+1:)%r_shell(3)=R_SH_ION(n_vacancies+1:reg_I_n_ions,3)
       DO i=n_vacancies+1,reg_I_n_ions
          DELS=-0.50*(point_0_shell_grad(I,:)+ point_0_core_grad(I,:))* &
               (dielectric_const+two)/(three*dielectric_const)/(ml_dprime_factors(epe(I)%k))
          R_SH_ION(I,:)=epe(I)%r(:)+DELS
         if(.not.relax_shells_only) R_NUC_ION(I,:)=epe(I)%r(:)+ &
               (DELS-(dielectric_const+two)/(three*dielectric_const)*point_0_core_grad(I,:)/PK(epe(I)%k))
          if(maxval(dels).gt.delmax) delmax=maxval(dels)
          if(maxval(point_0_shell_grad(I,:)+point_0_core_grad(I,:)).gt.gradmax) &
               gradmax = maxval(point_0_shell_grad(I,:)+point_0_core_grad(I,:))
          if( (dielectric_const+two)/(three*dielectric_const)/(ml_dprime_factors(epe(I)%k)).gt. &
               maxfac)  &
               maxfac= (dielectric_const+two)/(three*dielectric_const)/(ml_dprime_factors(epe(I)%k))  
       enddo ! 
       write(output_epe,*) 'init_regI_displacements delmax, gradmax', delmax,gradmax

       if(option_c3_symm) then
          call symmetrize_forces(r_sh_ion(:reg_I_n_ions,:))
          call symmetrize_forces(r_nuc_ion(:reg_I_n_ions,:))
       endif
     end subroutine init_regI_displacements
!----------------------------------------------

!----------------------------------------------
     subroutine ml_reg_I_displacements
!!$       print*,'ml_reg_I_displacements'
       DO i=n_vacancies+1,reg_I_n_ions
          regI_previous(i)%r_shell=R_SH_ION(I,1:3)
          regI_previous(i)%r_core=R_NUC_ION(I,1:3)

          R_SH_ION(I,:)=regI_previous(i)%r_shell-0.1_r8_kind* &
               point_0_shell_grad(I,:)*(dielectric_const+ &
               2.)/(3.*dielectric_const)/(2.*ml_dprime_factors(epe(I)%k))
          if(.not.relax_shells_only) &
          R_NUC_ION(I,:)=regI_previous(i)%r_core- 0.1_r8_kind* &
               (dielectric_const+2.)/(3.*dielectric_const)*point_0_core_grad(I,:)/PK(epe(I)%k)
       enddo
     end subroutine ml_reg_I_displacements
!-----------------------------------------------

!-----------------------------------------------
     subroutine init_e_calc
       write(output_epe, *) 'init_e_calc : start'
       select case(basic_action)
       case(0:1)
          if(n_types_central_atoms_3body > 0) then
             E3B=energy_3_body()
          else
             E3B=zero
          endif
          ETM=energy_shortrange(.true.)
          EPL=energy_coulomb(.true.)
          write(output_epe,*) 'shortrange, 3-body  and coulomb energies'
          write(output_epe,*) ETM,E3B,EPL
       case(2)
          E3B=zero
          ETM=zero
          EPL=zero
       endselect
       IF(epeit_count.EQ.0)  then
          etot_epe=ETM+E3B+EPL+defect_energy_short+defect_energy_coul
          DPRINT 'etot_epe comp:ETM+E3B+EPL+defect_energy_short+defect_energy_cou'
          DPRINT ETM,E3B,EPL,defect_energy_short,defect_energy_coul

          call get_epe_energies(var_epe_contrib=var_contrib)
  !          print*,'get_epe_energies 3 done'
          write(output_epe, *) 'epeit_count= ', epeit_count
          write(output_epe, *) 'etot_epe (a.u.)= ', etot_epe/eau_ev
          var_epe_contrib_atstart=var_contrib 
       endif
       write(output_epe, *) 'init_e_calc : finish'
     end subroutine init_e_calc
!----------------------------------------

!----------------------------------------
     subroutine output_r_imp
       DO  I=1,N_IMPURITIES
          write(output_epe, 13)I,R_SH_IMP(I,1)*RECA,R_SH_IMP(I,2)*RECA,R_SH_IMP(I,3)*RECA, &  
               R_NUC_IMP(I,1)*RECA,R_NUC_IMP(I,2)*RECA,R_NUC_IMP(I,3)*RECA     
13        format(1X,'impurity',i3,' coordinates of shell',3f9.4,  &
               /10X,                 ' coordinates of nuclear',3f9.4)
       enddo ! I=1,N_IMPURITIES
     end subroutine output_r_imp
!----------------------------------------

!---------------------------------------
     subroutine periodic_optim_actions

       CALL lattice_grad_gopt(r_sh_ion(1,1),r_sh_ion(1,2),r_sh_ion(1,3)  &
            ,r_nuc_ion(1,1),r_nuc_ion(1,2),r_nuc_ion(1,3)  &
            ,point_0_shell_grad(1,1),point_0_shell_grad(1,2),point_0_shell_grad(1,3)  &
            ,point_0_core_grad(1,1),point_0_core_grad(1,2),point_0_core_grad(1,3)  &
            ,short_range_grad(1,1),short_range_grad(1,2),short_range_grad(1,3))
        
       ETM=energy_shortrange_po(.true.)
       EPL=energy_coulomb_po(.true.)
       write(output_epe,*) 'coulomb energy of lattice ',EPL
       write(output_epe,*) 'short-range energy of lattice',ETM
       if(n_types_central_atoms_3body>0) then 
          E3B=energy_3_body_po()
          write(output_epe,*) 'three-body energy of lattice',E3B
       else
          E3B=zero
       endif
       etotal_gopt=ETM+EPL+E3B

       write(output_epe,*) 'energy of lattice ', ETM+EPL+E3B
       print*,'energy of lattice ',etotal_gopt ! ETM+EPL+E3B
       return
       write(gxcell_unit,'(2f20.7,3i5)')ETM+EPL+E3B,EPL,81,1,0

!!$ia=which_epe_ion(35)%new
!!$print*,short_range_grad(ia,2)+point_0_shell_grad(ia,2)+point_0_core_grad(ia,2)
!!$stop

       do io=1,n_ions_cell
          ia=which_epe_ion(io)%new
          write(gxcell_unit,'((i5,5x,3f15.6))') &
               io,scaling_factor*(short_range_grad(ia,1)+point_0_shell_grad(ia,1)+point_0_core_grad(ia,1))  &
               ,scaling_factor*(short_range_grad(ia,2)+point_0_shell_grad(ia,2)+point_0_core_grad(ia,2))  &
               ,scaling_factor*(short_range_grad(ia,3)+point_0_shell_grad(ia,3)+point_0_core_grad(ia,3))
       enddo ! i=1,n_ions_cell
!!$       write(output_epe,*) 'point_0_shell_grad+point_0_core_grad+short_range_grad'
!!$       write(output_epe,*) 1, point_0_shell_grad(1,:)+point_0_core_grad(1,:)+short_range_grad(1,:)
       point_0_shell_grad(1:reg_I_n_ions,:)=point_0_shell_grad(1:reg_I_n_ions,:)+short_range_grad(1:reg_I_n_ions,:)
     end subroutine periodic_optim_actions
!-----------------------------------

!-----------------------------------
     subroutine write_config
 !   print*, 'write_config'
       !write configuration for further use
       
       DPRINT 'write_config epeinout unit assigned to',out_unit
       out_unit=openget_iounit (trim(epe_input_dir)//'/epeinout' &
            ,form='unformatted',status='unknown')
       write(output_epe,*) ' save configuration of ions in region A  '
       WRITE(out_unit) reg_I_n_ions
       WRITE(out_unit) R_SH_ION(1:reg_I_n_ions,:),R_NUC_ION(1:reg_I_n_ions,:)
#ifdef NEW_EPE
       WRITE(out_unit) option_c3_symm
       WRITE(out_unit) (epe(i)%r,i=1,reg_I_n_ions)
       print*,'EPE CONFIGURATION WRITTEN'
#else
       WRITE(out_unit) point_0_shell_grad(1:reg_I_n_ions,1:3), &
            point_0_core_grad(1:reg_I_n_ions,1:3) 
#endif

       ! **step for interface with lcgto

       write(output_epe,*) 'gradients of initial point are writen on tape'      
       call returnclose_iounit(out_unit)
     end subroutine write_config
!--------------------------------------------

!--------------------------------------------------------------
     subroutine next_step_BFGS_gopt(stop_opt)

       use epecom_module
       use culon_module

       logical :: stop_opt
       real(kind=r8_kind), allocatable :: A_m(:,:),B_m(:,:)
       real(kind=r8_kind), allocatable, dimension(:) :: g_dim,g0_dim, &
            r_dim,r0_dim,step,g_step,r_step,buffer
       real(kind=r8_kind) :: store,delta_g,r1,r2,r3,rm,rm0,en1,en2,en3,enm
       real(kind=r8_kind) :: alpha,aa,bb,rnorm,gnorm,cosm

       integer(kind=i4_kind) :: i,j,k,l,l1,icount1,status
       logical :: resett
       real(kind=r8_kind), parameter :: small=1.0e-10_r8_kind

       allocate(buffer(n_ions_cell),stat=status)
       if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong allocation [1]") 
       buffer=which_epe_ion(:)%new
       if(r_max > 0.5_r8_kind) r_max=0.5_r8_kind

       allocate(g_dim(3*n_ions_cell), &
            g0_dim(3*n_ions_cell), &
            r_dim(3*n_ions_cell), & 
            r0_dim(3*n_ions_cell), &
            step(3*n_ions_cell),stat=status)
       if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong allocation [2]")

       resett=.false.
       icount1=0

1      continue

       do i=1,n_ions_cell
          l=buffer(i)
          if(epeit_count == 1 ) then
             l1=l
          else
             l1=which_epe_ion(i)%old
          endif
          do j=1,3
             k=3*(i-1)+j
             g_dim(k)=point_1_shell_grad(l,j)+point_1_core_grad(l,j)
             g0_dim(k)=point_0_shell_grad(l1,j)+point_0_core_grad(l1,j)
             r_dim(k)=r_ion_in_cell(i,j)
             r0_dim(k)=r_ion_in_cell_old(i,j)
          enddo
       enddo   
  
       delta_g=zero
       do i=1,n_ions_cell
          l=which_epe_ion(i)%new
          do j=1,3
             if(abs(point_1_shell_grad(l,j)+point_1_core_grad(l,j)) > delta_g) &
                  delta_g=abs(point_1_shell_grad(l,j)+point_1_core_grad(l,j))
          enddo
       enddo

       if (epeit_count > 1 .and. .not. resett) then
          if (delta_g <= abs_g) then
             write (output_epe,*) '================================='
             write (output_epe,*) '     Minimun has been found'
             write (output_epe,*) '================================='
             stop_opt=.true.
             deallocate(H_m1, &
                  g_dim, g0_dim, &
                  r_dim, r0_dim, &
                  step,buffer,stat=status)
             if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong deallocation [3]")
             return
          endif
       endif

       rnorm=sqrt(dot_product(r_dim-r0_dim,r_dim-r0_dim))
       gnorm=sqrt(dot_product(g_dim,g_dim))
       cosm=abs(dot_product(r_dim-r0_dim,g_dim))/(rnorm*gnorm)
!!$print*,cosm,'cosm',gnorm/(3*n_ions_cell),'Gnorm'
       if(.not.resett) then
          write (output_epe,*) 'Gnorm= ',gnorm/(3*n_ions_cell)
          write (*,*) 'Gnorm= ',gnorm/(3*n_ions_cell)
       endif
       if(cosm < 0.01_r8_kind .and. epeit_count > 1) then
          resett=.true.
       endif

       alpha=1.0_r8_kind
       if ( epeit_count == 1 .or. &
            (mod(epeit_count-1,n_hess_update_cycles) == 0 .and. epeit_count > 1).or. &
            resett ) then
          if(.not. allocated(H_m1)) then
             allocate(H_m1(3*n_ions_cell,3*n_ions_cell),stat=status)
             if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong allocation [4]")
          endif
          if(epeit_count /=1) then
             write (*,*) '     Hessian reset'
             write (output_epe,*) '================================='
             write (output_epe,*) '          Hessian reset'
             write (output_epe,*) '================================='
          endif
          H_m1=zero
          do i=1,3*n_ions_cell
             H_m1(i,i)=0.5_r8_kind*weight_hess
          enddo
          etotal_gopt_0=etotal_gopt
          goto 2
       endif
       
       if (epeit_count > 1 .and. mod(epeit_count-1,n_hess_update_cycles) /= 0 .and. .not.resett) then
          allocate(A_m(3*n_ions_cell,3*n_ions_cell),B_m(3*n_ions_cell,3*n_ions_cell), &
!!$               C_m(3*n_ions_cell,3*n_ions_cell), &
!!$               lda(3*n_ions_cell),ldb(3*n_ions_cell),ldc(3*n_ions_cell), &
!!$               store_v0(3*n_ions_cell), &
!!$               store_v1(3*n_ions_cell), &
               r_step(3*n_ions_cell), &
               g_step(3*n_ions_cell),stat=status)
          if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong allocation [5]")

          r_step=r_dim-r0_dim
          g_step=g_dim-g0_dim
!vvvvvvvvvvvvvvvvvvvvvvv       
!!$          store=dot_product(matmul(g_step,H_m1),g_step)
!!$          store_v0=matmul(H_m1,g_step)
!!$          store_v1=matmul(g_step,H_m1)
!!$          do i=1,3*n_ions_cell
!!$             do j=1,3*n_ions_cell
!!$                A_m(i,j)=store_v0(i)*store_v1(j)/store
!!$             enddo
!!$          enddo
!!$          H_m1=H_m1-A_m
!^^^^^^^^^^^^^^^^^^^^^^^^^^^
          
          store=dot_product(r_step,g_step)

!vvvvvvvvvvvvvvvvvvvvvvv
          do i=1,3*n_ions_cell
             do j=1,3*n_ions_cell
                A_m(i,j)=r_step(i)*g_step(j)/store
                if(i==j) then 
                   A_m(i,j)=1.0_r8_kind-A_m(i,j)
                else
                   A_m(i,j)=-A_m(i,j)
                endif
             enddo
          enddo
!!$          B_m=matmul(A_m,H_m1)
!!$          call dgemm('n','n',3*n_ions_cell,3*n_ions_cell,3*n_ions_cell, &
!!$               alph,A_m,lda,H_m1,ldb,bet,B_m,ldc)
          do i=1,3*n_ions_cell
             do j=i,3*n_ions_cell
                B_m(i,j)=0.0_r8_kind
                do k=1,3*n_ions_cell
                   B_m(i,j)=B_m(i,j)+A_m(i,k)*H_m1(k,j)
                enddo
                B_m(j,i)=B_m(i,j)
             enddo
          enddo


!!$          C_m=transpose(A_m)

!!$          H_m1=matmul(B_m,C_m)
          do i=1,3*n_ions_cell
             do j=i,3*n_ions_cell
                H_m1(i,j)=0.0_r8_kind
                do k=1,3*n_ions_cell
                   H_m1(i,j)=H_m1(i,j)+B_m(i,k)*A_m(j,k)
                enddo
                H_m1(j,i)=H_m1(i,j)
             enddo
          enddo
!^^^^^^^^^^^^^^^^^^^^^^^^^^^

          do i=1,3*n_ions_cell
             do j=1,3*n_ions_cell
                A_m(i,j)=r_step(i)*r_step(j)/store
             enddo
          enddo
          H_m1=H_m1+A_m

          deallocate(A_m, B_m, & !C_m, lda, ldb, ldc, &
!!$               store_v0, store_v1, &
               r_step, g_step,stat=status)
          if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong deallocation [6]")
       endif

2      r0_dim=r_dim

!..............................................
       resett=.false.
       
       r1=zero
       en1=etotal_gopt

       do
          step=-alpha*matmul(H_m1,g_dim)
          r3=sqrt(dot_product(step,step))
          if (r3 <= 2.0_r8_kind*r_max+small) exit
          alpha=alpha*2.0_r8_kind*r_max/r3
       enddo

       r_dim=r0_dim+step
       do i=1,n_ions_cell
          do j=1,3
             k=3*(i-1)+j
             r_ion_in_cell(i,j)=r_dim(k)
          enddo
       enddo
       call calc_lattice !(2)
 !   print*,'calc_lattice 2 done'
       if(n_types_central_atoms_3body > 0) then
          en3=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)+energy_3_body_po()
       else
          en3=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)
       endif

       step=-alpha*matmul(H_m1,g_dim)/5.0_r8_kind
       r2=sqrt(dot_product(step,step))
       r_dim=r0_dim+step
       do i=1,n_ions_cell
          do j=1,3
             k=3*(i-1)+j
             r_ion_in_cell(i,j)=r_dim(k)
          enddo
       enddo
       call calc_lattice !(3)
 !   print*,'calc_lattice 3 done'
       if(n_types_central_atoms_3body > 0) then
          en2=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)+energy_3_body_po()
       else
          en2=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)
       endif

       rm0=r2

       aa=((en1-en2)*(r1-r3)-(en1-en3)*(r1-r2))/ &
            ((r1-r3)*(r1*r1-r2*r2)-(r1-r2)*(r1*r1-r3*r3))
       if(aa > zero) then
          bb=((en1-en2)-aa*(r1*r1-r2*r2))/(r1-r2)
          rm=-bb/(2.0_r8_kind*aa)
       else
          rm=min(abs(r2),abs(r3))/5.0_r8_kind
       endif
       if(abs(rm) > 2.0_r8_kind) rm=rm0/2.0_r8_kind
!!$       r_dim=r0_dim+(rm/rm0)*(r_dim-r0_dim)
       r_dim=r0_dim+(rm/rm0)*step
       do i=1,n_ions_cell
          do j=1,3
             k=3*(i-1)+j
             r_ion_in_cell(i,j)=r_dim(k)
          enddo
       enddo
       call calc_lattice !(4)
 !   print*,'calc_lattice 4 done'
       if(n_types_central_atoms_3body > 0) then
          enm=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)+energy_3_body_po()
       else
          enm=energy_shortrange_po(.false.)+energy_coulomb_po(.false.)
       endif

!!$       print*,r1,r2,r3,rm
!!$       print*,en1,en2,en3,enm

!!$       if(enm < etotal_gopt .and. rm > zero) then
       if(enm <= etotal_gopt) then
       else       
          resett=.true.
          icount1=icount1+1
          if(icount1 == 1) then
             do i=1,n_ions_cell
                do j=1,3
                   k=3*(i-1)+j
                   r_ion_in_cell(i,j)=r0_dim(k)
                enddo
             enddo
             goto 1
          else
             write (output_epe,*) '================================='
             write (output_epe,*) '   Next step cannot be done.'
             write (output_epe,*) '================================='
             write (*,*) '     Next step cannot be done.'
             deallocate(H_m1,g_dim, g0_dim, &
                  r_dim, r0_dim, &
                  step,buffer,stat=status)
             if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong deallocation [6]")
             stop_opt=.true.
             return  
          endif
       endif
!...........................................

       do i=1,n_ions_cell
          do j=1,3
             k=3*(i-1)+j
             r_ion_in_cell(i,j)=r_dim(k)
             r_ion_in_cell_old(i,j)=r0_dim(k)
          enddo
       enddo

       point_0_shell_grad=point_1_shell_grad
       point_0_core_grad=point_1_core_grad

       which_epe_ion(:)%old=buffer


       etotal_gopt_0=etotal_gopt
       etotal_gopt=enm
       r_max=abs(rm)

       deallocate(g_dim, g0_dim, &
            r_dim, r0_dim, &
            step,buffer,stat=status)
       if(status.ne.0) call error_handler("next_step_BFGS_gopt: wrong deallocation [7]")

     end subroutine next_step_BFGS_gopt

   end subroutine main_epe

! !*************************************************************
! subroutine start_slave_lattice( parallel_calculations, message_num)
!   !  Purpose: Start lattice grads. calculations on slaves.
!   !           Called by master.
!   !------------ Modules used ----------------------------------
!   implicit none
!  !------------ Declaration of local variables -----------------
!  integer, intent(in)  :: message_num
!  logical, intent(in)  :: parallel_calculations
!  !------------ Declaration of subroutines used ----------------
!  external error_handler
!  !------------ Executable code --------------------------------
!   if(parallel_calculations) then
!   call epe_send_shells_and_cores(message_num) ! + allocation 
!   end if
! end subroutine start_slave_lattice
! !*************************************************************
        
!AG (moved form main_epe)
   subroutine cons_latt_gradients()
       
     call start_timer(timer_latt_epe_forces)
     CALL epe_forces(r_sh_ion(:,1),r_sh_ion(:,2),r_sh_ion(:,3)  &
          ,r_nuc_ion(:,1),r_nuc_ion(:,2),r_nuc_ion(:,3)  &
          ,point_1_shell_grad(:,1),point_1_shell_grad(:,2),point_1_shell_grad(:,3)      &
          ,point_1_core_grad(:,1),point_1_core_grad(:,2),point_1_core_grad(:,3) )
     call stop_timer(timer_latt_epe_forces) 
     if( comm_i_am_master() .and. comm_parallel() ) &
          call epe_collect_forces( &
          point_1_shell_grad(:,1),point_1_shell_grad(:,2),point_1_shell_grad(:,3), &
          point_1_core_grad(:,1),point_1_core_grad(:,2),point_1_core_grad(:,3))
   end subroutine cons_latt_gradients
!*****************************************************

!****************************************************
   subroutine init_lat_gradients()
     call start_timer(timer_latt_epe_forces)
     CALL epe_forces(r_sh_ion(:,1),r_sh_ion(:,2),r_sh_ion(:,3)  &
          ,r_nuc_ion(:,1),r_nuc_ion(:,2),r_nuc_ion(:,3)  &
          ,point_0_shell_grad(:,1),point_0_shell_grad(:,2),point_0_shell_grad(:,3)  &
          ,point_0_core_grad(:,1), point_0_core_grad(:,2), point_0_core_grad(:,3)  )
     call stop_timer(timer_latt_epe_forces)
     if( comm_i_am_master() .and. comm_parallel() ) then
        call epe_collect_forces( &
             point_0_shell_grad(:,1),point_0_shell_grad(:,2),point_0_shell_grad(:,3), &
             point_0_core_grad(:,1), point_0_core_grad(:,2), point_0_core_grad(:,3)  )
     end if
        
   end subroutine init_lat_gradients
!**************************************************************

!*************************************************************
   subroutine start_slave_defect( parallel_calculations, message_num)
     !  Purpose: Start calculations of defect contributions
     !           on slaves. Called by master.
     !------------ Modules used ----------------------------------
     implicit none
     !------------ Declaration of local variables -----------------
     integer, intent(in)  :: message_num
     logical, intent(in)  :: parallel_calculations
     !------------ Executable code --------------------------------
     if(parallel_calculations) then
        select case(message_num)
        case(2)
           call comm_init_send(comm_all_other_hosts, msgtag_epe_defects)
        case(3)
           call comm_init_send(comm_all_other_hosts, msgtag_epe_consdef)
        case(4)
           call comm_init_send(comm_all_other_hosts, msgtag_epe_def_fin)
        end select
        call comm_send() 
     end if
   end subroutine start_slave_defect
!*************************************************************

!*************************************************************
!    (moved from subroutine main_epe) 
   subroutine cons_defect_contrubutions()
   if(n_impurities.ne.0.or.n_vacancies.ne.0) then
     if( comm_i_am_master() ) call start_timer(timer_defect_contributions)
     imp_conf=>impurities_new_conf
     select case(pg_interfaced_mode)
     case(0)
        CALL calc_def_contributions(point_1_core_grad,point_1_shell_grad,short_range_grad, &
             curr%imp_core_grad,curr%imp_shell_grad,use_epe_reference)
     case(1:2)
        CALL calc_def_contributions(point_1_core_grad,point_1_shell_grad,short_range_grad, &
             curr%imp_core_grad,curr%imp_shell_grad,use_epe_reference)
     endselect
     if( comm_i_am_master() ) call stop_timer(timer_defect_contributions)
     !AG[
     if( comm_i_am_master() .and. comm_parallel() ) then
        call start_timer(timer_collect_def_contributions)
        call epe_collect_def_contributions(point_1_core_grad,point_1_shell_grad,short_range_grad, &
             curr%imp_core_grad,curr%imp_shell_grad)
        call stop_timer(timer_collect_def_contributions)
     end if
     !AG]
  endif
!!$     call write_to_output_units('point_1_shell_grad completed ')              
!!$     print*,'point_1_shell_grad completed '
     point_1_shell_grad(1:reg_I_n_ions,:)=point_1_shell_grad(1:reg_I_n_ions,:)+ &
          short_range_grad(1:reg_I_n_ions,:)
     if(option_c3_symm) then
        call symmetrize_forces(point_1_core_grad)
        call symmetrize_forces(point_1_shell_grad)
     endif
     if(comm_i_am_master() ) call print_energy()

   end subroutine cons_defect_contrubutions
!***********************************************************

!***********************************************************
   subroutine defect_contributions()
     if( comm_i_am_master() ) call start_timer(timer_defect_contributions)
     imp_conf=>impurities_old_conf

     DPRINT 'calc_def_cont for pg_interface and use_epe_reference', pg_interfaced_mode,use_epe_reference
     select case(pg_interfaced_mode)
     case(0)
 !      print*,'CALL calc_def_contributions'
        CALL calc_def_contributions(point_0_core_grad,point_0_shell_grad, &
             short_range_grad,prev%imp_core_grad,prev%imp_shell_grad,use_epe_reference) 
 !      print*,'calc_def_contributions done'
     case(1:2)
        CALL calc_def_contributions(point_0_core_grad,point_0_shell_grad, &
             short_range_grad,prev%imp_core_grad,prev%imp_shell_grad,use_epe_reference) 
     endselect

     if( comm_i_am_master() ) call stop_timer(timer_defect_contributions)
     if( comm_i_am_master() .and. comm_parallel() ) then
        call start_timer(timer_collect_def_contributions)
        call epe_collect_def_contributions(point_0_core_grad,point_0_shell_grad,short_range_grad, &
        prev%imp_core_grad,prev%imp_shell_grad)
        call stop_timer(timer_collect_def_contributions)
     end if

     point_0_shell_grad(:reg_I_n_ions,:)=point_0_shell_grad(:reg_I_n_ions,:)+ &
          short_range_grad(:reg_I_n_ions,:)
     if(option_c3_symm) then
        call symmetrize_forces(point_0_core_grad)
        call symmetrize_forces(point_0_shell_grad)
     endif
   end subroutine defect_contributions
!*************************************************

!*************************************************************
   subroutine defect_contributions_fin()
     !  Purpose: Final call of "calc_def_contributions"
     !           on slaves. Called both by master and slave.
     !------------ Modules used ----------------------------------
     implicit none
     !------------ Declaration of local variables -----------------
     !------------ Executable code --------------------------------
     if( comm_i_am_master() ) call start_timer(timer_defect_contributions)

     CALL calc_def_contributions(point_1_core_grad,point_1_shell_grad, &
          short_range_grad,curr%imp_core_grad, curr%imp_shell_grad,use_epe_reference)
     if( comm_i_am_master() ) call stop_timer(timer_defect_contributions)
     if( comm_i_am_master() .and. comm_parallel() ) then
        call start_timer(timer_collect_def_contributions)
        call epe_collect_def_contributions(point_1_core_grad,point_1_shell_grad,short_range_grad, &
             curr%imp_core_grad,curr%imp_shell_grad)
        call stop_timer(timer_collect_def_contributions)
     end if
   end subroutine defect_contributions_fin

   subroutine print_energy
     real(kind=r8_kind) :: var_contrib
     select case(basic_action)
     case(0:1)
        ETM=energy_shortrange(.true.)
        EPL=energy_coulomb(.true.)
        if(n_types_central_atoms_3body>0) then
           E3B=energy_3_body()
           write(output_epe, *) 'REGULAR LATTICE THREE BODY ENERGY   =',E3B
        else
           E3B=zero
        endif
     end select

#if 1
     etot_epe=ETM+EPL+E3B+defect_energy_short+defect_energy_coul
                          ! here contrib from imp_3_body_1
#else
     etot_epe=defect_energy_short+defect_energy_coul+E3B+ETM
#endif
     
     call get_epe_energies(var_epe_contrib=var_contrib)
 !         print*,'get_epe_energies 4 done'
     write(output_epe, *) 'epeit_count= ', epeit_count
     write(output_epe, *) 'etot_epe (a.u.)= ', etot_epe/eau_ev
     print*, 'epeit_count= ', epeit_count
     print*, 'etot_epe (a.u.)= ', etot_epe/eau_ev
     if(epeit_count.eq.4) e_epe_atstart=etot_epe
   end subroutine print_energy
!**********************************************************

!*************************************************************
   subroutine epe_collect_forces(DSX,DSY,DSZ,DCX,DCY,DCZ)
     !  Purpose: Receives epe_forces calculated by slaves.
     !           Called by master.
     !------------ Modules used ----------------------------------
     implicit none
     real(kind=r8_kind), dimension(reg_I_n_ions)   :: DCX, DCY, DCZ, &
                                                    DSX, DSY, DSZ
     !------------ Declaration of local variables -----------------
     real(kind=r8_kind), allocatable, dimension(:)   :: DCX_send, DCY_send, DCZ_send, &
          DSX_send, DSY_send, DSZ_send  
     !                                                     temp. receiving arrays
     real(kind=r8_kind), allocatable, dimension(:,:) :: short_range_grad_send
     integer             :: status, n_received
     logical             :: all_epe_forces_received
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
!!$     call write_to_output_units("Master: collect epe_forces")
     n_received=0
    
     all_epe_forces_received=.false.
     do while(.not.all_epe_forces_received )
        
        do while ( .not. comm_save_recv_nonblocking(comm_all_other_hosts, &
                                         msgtag_epe_forces_done) )
        end do
        n_received=n_received+1

        if(.not.allocated(DCX_send)) then
           allocate ( DCX_send(reg_I_n_ions), & 
                DCY_send(reg_I_n_ions), &
                DCZ_send(reg_I_n_ions), &
                DSX_send(reg_I_n_ions), &
                DSY_send(reg_I_n_ions), &
                DSZ_send(reg_I_n_ions),stat=status  ) 
           if(status.ne.0) then
              call error_handler &
                   ("epe_collect_forces :DC DS alloc failed  ")
           end if
        end if
    
        if(.not.allocated(short_range_grad_send)) then
           allocate( short_range_grad_send(reg_I_n_ions,3),stat=status)
           if(status.ne.0) call error_handler &
                   ("epe_collect_forces :short_range_grad_send alloc failed ")
        end if
    
        call communpack(DCX_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [3]")
        call communpack(DCY_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [4]")
        call communpack(DCZ_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [5]")
        call communpack(DSX_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [6]")
        call communpack(DSY_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [7]")
        call communpack(DSZ_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [8]")
        call communpack(short_range_grad_send(1,1), reg_I_n_ions*3, 1, status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces :  error [9]")

        DCX(:) = DCX(:) + DCX_send(:)
        DCY(:) = DCY(:) + DCY_send(:)
        DCZ(:) = DCZ(:) + DCZ_send(:)
        DSX(:) = DSX(:) + DSX_send(:)
        DSY(:) = DSY(:) + DSY_send(:)
        DSZ(:) = DSZ(:) + DSZ_send(:)
        short_range_grad(:,:)=short_range_grad(:,:) + short_range_grad_send(:,:)

        deallocate ( DCX_send, DCY_send, DCZ_send, &
             DSX_send, DSY_send, DSZ_send,stat=status  )
        if(status.ne.0) call error_handler &
             ("epe_collect_forces : deallocate DCX_send failed")   
        
        deallocate ( short_range_grad_send ,stat=status)
        if(status.ne.0) call error_handler &
             ("epe_collect_forces : deallocate short_range_grad_send failed " )
        
        all_epe_forces_received=n_received==comm_get_n_processors()-1
     end do
!!$     call write_to_output_units('Master: epe_forces received')                    
   end subroutine epe_collect_forces
!*************************************************************

!*************************************************************
   subroutine epe_collect_def_contributions(point_1_core_grad,point_1_shell_grad,  &
        short_range_grad1, &
        impurity_core_grad,impurity_shell_grad )
     !  Purpose: Receives epe_defect contributions calculated 
     !           by slaves. Called by master.
     !------------ Modules used ----------------------------------
     implicit none 
     real(kind=r8_kind), intent(inout), &
          dimension(reg_I_n_ions,3)   :: point_1_core_grad, &
          point_1_shell_grad
     real(kind=r8_kind), intent(inout), &
          dimension(reg_I_n_ions,3)   :: short_range_grad1
     real(kind=r8_kind), intent(inout), &
          dimension(n_impurities,3)           :: impurity_core_grad, &
          impurity_shell_grad
     !------------ Declaration of local variables -----------------
     real(kind=r8_kind), allocatable, dimension(:,:) :: impurity_core_grad_send,        &
                                                      impurity_shell_grad_send,       &
                                                      point_1_core_grad_send,         &
                                                      point_1_shell_grad_send,        &
                                                      short_range_grad_send
     real(kind=r8_kind), allocatable, dimension(:)   :: temp_array
     real(kind=r8_kind), parameter                   :: zero=0.0_r8_kind
     integer                                         :: status, n_received
     logical                                         :: all_def_contr_received
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
!!$     call write_to_output_units("Master: collect epe_def_contrib.")
     n_received=0
     all_def_contr_received=.false.
     allocate( impurity_core_grad_send(n_impurities,3) , impurity_shell_grad_send(n_impurities,3) &
          ,temp_array(19) , point_1_core_grad_send(reg_I_n_ions,3) &
          ,point_1_shell_grad_send(reg_I_n_ions,3), short_range_grad_send(reg_I_n_ions,3) , &
          stat=status )
     if(status.ne.0) call error_handler("allocate impurity_core_grad_send failed")

     do while( .not.all_def_contr_received )

        do while ( .not. comm_save_recv_nonblocking(comm_all_other_hosts, &
                                         msgtag_epe_send_def) )
        end do
        n_received=n_received+1

        impurity_core_grad_send(:,:)  = zero
        impurity_shell_grad_send(:,:) = zero 
        point_1_core_grad_send(:,:)   = zero
        point_1_shell_grad_send(:,:)  = zero
        short_range_grad_send(:,:)    = zero
        temp_array(:)                 = zero   
        call communpack(point_1_core_grad_send(1,1),  3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_rec_defect_contrib:  error [1]")
        call communpack(point_1_shell_grad_send(1,1), 3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_rec_defect_contrib:  error [2]")
        call communpack(short_range_grad_send(1,1),   3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_rec_defect_contrib:  error [3]")    
        call communpack(impurity_core_grad_send(1,1), 3*n_impurities, 1, status)
        if(status.ne.0) call error_handler("epe_rec_defect_contrib:  error [4]")
        call communpack(impurity_shell_grad_send(1,1),3*n_impurities, 1, status)
        if(status.ne.0) call error_handler("epe_rec_defect_contrib:  error [5]")
 
        point_1_core_grad(:,:)   = point_1_core_grad(:,:)   + point_1_core_grad_send(:,:)
        point_1_shell_grad(:,:)  = point_1_shell_grad(:,:)  + point_1_shell_grad_send(:,:)
        short_range_grad1(:,:)    = short_range_grad1(:,:)    + short_range_grad_send(:,:)
        impurity_core_grad(:,:)  = impurity_core_grad(:,:)  + impurity_core_grad_send(:,:)
        impurity_shell_grad(:,:) = impurity_shell_grad(:,:) + impurity_shell_grad_send(:,:)
        all_def_contr_received=n_received==comm_get_n_processors()-1
     end do
     deallocate( impurity_core_grad_send, impurity_shell_grad_send, temp_array, &
          point_1_core_grad_send , point_1_shell_grad_send , short_range_grad_send ,stat=status )
     if(status.ne.0) call error_handler("epe_rec_defect_contrib: deallocate failed")
!!$     call write_to_output_units('Master: epe_def_contrib. received')              
   end subroutine epe_collect_def_contributions

   subroutine finish_epe(do_deallocate)
     use epe_pg_module 

     logical, intent(in) :: do_deallocate
     real(kind=r8_kind)::e_lattice,var_contrib
     integer(kind=i4_kind):: allocation_status,i_alloc
     
     print*, 'finish_epe:  optimization is to be  finished'
     write(output_epe,*) 'finish_epe:  optimization is finished'
        if(periodic_optimization) return

! **DISPM,DISTM - modulus of DISP and d_sh_nuc
! **CKALP,CKALT - scalar products of the vectors R and DISP, R and d_sh_nuc
! **UGOLT,UGOLP - cosine of angles between vectors R and d_sh_nuc, R and DISP
! **POL - energy of electron polarization of L - ion
     deallocate(regI_previous,stat=allocation_status)
     if(allocation_status.ne.0) call error_handler("deallocate regI_previous failed")


     IF(print_ions_from_first_sphere) then
        allocate(d_r_sh(reg_I_n_ions,3),d_sh_nuc(reg_I_n_ions,3),stat=allocation_status)
        if(allocation_status.ne.0) call error_handler('d_r_sh allocation failed')

        d_r_sh(1:reg_I_n_ions,1)=R_SH_ION(1:reg_I_n_ions,1)-epe(1:reg_I_n_ions)%r(1)
        d_r_sh(1:reg_I_n_ions,2)=R_SH_ION(1:reg_I_n_ions,2)-epe(1:reg_I_n_ions)%r(2)
        d_r_sh(1:reg_I_n_ions,3)=R_SH_ION(1:reg_I_n_ions,3)-epe(1:reg_I_n_ions)%r(3)
        d_sh_nuc(1:reg_I_n_ions,:)=R_SH_ION(1:reg_I_n_ions,:)-R_NUC_ION(1:reg_I_n_ions,:)

!!$    do i=1,n_ions_cell
!!$1116   FORMAT (I3,I2,1X,A4,3F12.6)
!!$       write(output_epe,1116) i,type_of_ion(i),NAME_OF_IONS(i),r_ion_in_cell(i,:)+d_r_sh(which_epe_ion(i)%new,:)
!!$       write(output_epe,1116) i,type_of_ion(i),NAME_OF_IONS(i),r_ion_in_cell(i,:)+ &
!!$             d_r_sh(which_epe_ion(i)%new,:)-d_sh_nuc(which_epe_ion(i)%new,:)
!!$     enddo

        L1=n_vacancies+1
        write(output_epe, 775)
775     FORMAT(1X,78('-'))
        write(output_epe, 773)
773     FORMAT(1X,' L Q_ZL epe(L)%r                         ! R_SH_ION-epe%r  ',  &
             '     d_r_sh   UGOLP')
        write(output_epe,1773)
1773    FORMAT(41X,'  R_SH_ION-R_NUC_ION   DISTM   UGOLT   POL')
        write(output_epe, 775)
        DO l=n_vacancies+1,reg_I_n_ions
           DISPM=0.0
           DISTM=0.0
           CKALP=0.0
           CKALT=0.0
           DISPM=dot_product(d_r_sh(l,:),d_r_sh(l,:))
           DISTM=dot_product(d_sh_nuc(L,:),d_sh_nuc(L,:))
           CKALP=dot_product(epe(L)%r-R_CENT_GENER(epe(l)%gc,:),d_r_sh(L,:))
           CKALT=dot_product(epe(L)%r-R_CENT_GENER(epe(l)%gc,:),d_sh_nuc(l,:))
           POL=DISTM*PK(epe(L)%k)*0.5
           DISTM=SQRT(DISTM)
           DISPM=SQRT(DISPM)
           UGOLT=1.
           UGOLP=1.
           IF(DISTM.GT.1.E-6.and.epe(l)%d.gt.1.E-6) UGOLT=CKALT/epe(L)%d/DISTM
           IF(DISPM.GT.1.E-6.and.epe(l)%d.gt.1.E-6) UGOLP=CKALP/epe(L)%d/DISPM
           DISPM=DISPM*RECA
           DISTM=DISTM*RECA
           write(output_epe, 776)L,Q_ZL(L),epe(L)%r*RECA, d_r_sh(L,:)*RECA, &
                DISPM,UGOLP
776        FORMAT(1X,I3,4F8.3,' !',3F8.4,F7.4,F8.4,6x,' !')
           write(output_epe,1776)d_sh_nuc(L,1)*RECA,d_sh_nuc(L,2)*RECA,d_sh_nuc(L,3)*RECA, & !STX,STY,STZ
                DISTM,UGOLT,POL
1776       FORMAT(36X,' !',3F8.4,F7.4,F8.4,F7.4,'!')
           IF(L.NE.epe(L1)%u) cycle
           write(output_epe, 775)
           L1=L1+1
        enddo ! l=n_vacancies+1,reg_I_n_ions
        write(output_epe, 775)
        deallocate(d_r_sh,d_sh_nuc,stat=allocation_status)
        if(allocation_status.ne.0) call error_handler("deallocate d_r_sh failed")
     endif !print_ions_from_first_sphere

     IF(N_IMPURITIES.ne.0) then
        write(output_epe, 116)
116     format(1x,' c o o r d i n a t e s  of  i m p u r i t i e s')
        DO I=1,N_IMPURITIES
           write(output_epe, 13)I,R_SH_IMP(I,1)*RECA,R_SH_IMP(I,2)*RECA,R_SH_IMP(I,3)*RECA, & !SPX,SPY,SPZ
                R_NUC_IMP(I,1)*RECA,R_NUC_IMP(I,2)*RECA,R_NUC_IMP(I,3)*RECA    !CPX,CPY,CPZ
13         format(1X,'impurity',i3,' coordinates of shell',3f9.4,  &
                /10X,                 ' coordinates of core   ',3f9.4)
        enddo
     endif

     imp_conf=>impurities_new_conf
     if( comm_parallel() )  call epe_send_shells_and_cores(msg_defects)
     call start_slave_defect( comm_parallel(), msg_def_fin )
     call start_timer(timer_defect_ext)
     call defect_contributions_fin()
     call stop_timer(timer_defect_ext)
 !print*,'finish epe energy_shortrange 1'
     ETM=energy_shortrange(.true.)
     EPL=energy_coulomb(.true.)
     if(n_types_central_atoms_3body>0) then
        E3B=energy_3_body()
        write(output_epe, *) 'THREE BODY ENERGY   =',E3B
     else
        E3B=zero
     endif
     etot_epe=ETM+EPL+E3B+defect_energy_short+defect_energy_coul
     call get_epe_energies(etot_epe_corrected=e_lattice)
 !          print*,'get_epe_energies 5 done'
     call get_epe_energies(var_epe_contrib=var_contrib)
 !          print*,'get_epe_energies 6 done'
     write(output_epe, *)'etot_epe,var_contrib (a.u.)  epeit_count ',epeit_count
     write(output_epe, *) etot_epe/eau_ev,var_contrib
     write(output_epe, *) 'end of epe minimization '

     out_unit=openget_iounit(trim(epe_input_dir)//'/epeinout' &
          ,form='unformatted',status='unknown')
     WRITE(out_unit)reg_I_n_ions
     WRITE(out_unit) R_SH_ION(1:reg_I_n_ions,1:3),R_NUC_ION(1:reg_I_n_ions,1:3)
#ifdef NEW_EPE
     WRITE(out_unit) option_c3_symm
     WRITE(out_unit) (epe(i)%r,i=1,reg_I_n_ions)
     print*, 'finish: EPE CONFIGURATION IS WRITTEN'
#else
     WRITE(out_unit) point_0_shell_grad,point_0_core_grad
#endif 
     ! write gradients of the current
     ! step for interface with lhdscom
     call returnclose_iounit(out_unit)
     
     if(qm_interfaced_mode) then
        imp_conf=>impurities_new_conf
        call epe2pg(curr%imp_core_grad,curr%imp_shell_grad)
     endif
     if(use_epe_reference.and.qm_interfaced_mode)  then 
        if(allocated(reg_reference) ) then
           deallocate(reg_reference,stat=allocation_status)
           if(allocation_status.ne.0) call error_handler(" reg_reference deallocate failed")
        end if

        if(allocated(epe_reference)) then
           deallocate(epe_reference,stat=allocation_status)
           if(allocation_status.ne.0) call error_handler(" epe_reference deallocate failed")
        end if
        
     end if
     
     deallocate(reg_I_pg,point_0_core_grad,point_0_shell_grad, &
          point_1_core_grad,point_1_shell_grad,reg1_shells_coordinates, &
          reg1_cores_coordinates,&
          short_range_grad,stat=epealloc_stat(21))
          ! reg_I_pg:
          ASSERT(epealloc_stat(21).eq.0)
          epealloc_stat(23)=0 ! point_0_core_grad point_0_shell_grad point_1_core_grad point_1_shell_grad short_range_grad
          epealloc_stat(24)=0 ! reg1_shells_coordinates reg1_cores_coordinates 

     if(do_deallocate) then
        deallocate(r_nuc_ion,r_sh_ion,epe,q_zl,stat=epealloc_stat(9))
        ASSERT(epealloc_stat(9).eq.0)
        epealloc_stat(10)=0 ! q_zl
     end if

     if(qm_interfaced_mode) then
        deallocate(lato,nhds,gx_imp_r,stat=epealloc_stat(5)) ! 5 - gx_imp_r
           ASSERT(epealloc_stat(5).eq.0)
           epealloc_stat(11)=0 ! lato
           epealloc_stat(12)=0
     else 
       if(allocated(gx_imp_r)) then
        deallocate(gx_imp_r,stat=epealloc_stat(5))
        ASSERT(epealloc_stat(5).eq.0)
       endif
     DPRINT 'lato epealloc_stat(11)',allocated(lato),allocated(nhds)
        if(allocated(lato)) then
         deallocate(lato,nhds,stat=epealloc_stat(11))
         ASSERT(epealloc_stat(11).eq.0)
         epealloc_stat(12)=0  ! nhds
        endif
      
     end if
        
     if(do_deallocate) call returnclose_iounit(output_epe)

     if(unit_cell_output)call returnclose_iounit(ucellxyz)

     if(n_impurities.ne.0) then
        deallocate(impurities_old_conf,impurities_new_conf,stat=epealloc_stat(22))
        ASSERT(epealloc_stat(22).eq.0)

        deallocate(curr%imp_core_grad, curr%imp_shell_grad, & !finish
                   prev%imp_core_grad, prev%imp_shell_grad,stat=epealloc_stat(15))
        ASSERT(epealloc_stat(15).eq.0)
     end if
        
     if(n_types_central_atoms_3body >0  .and. do_deallocate) then
        deallocate(ki,theta_0,types,r3b,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
        epealloc_stat(14)=0 ! types
     endif

     if(n_ec_types_central_atoms_3body >0  .and. do_deallocate) then
        deallocate(ec%ki,ec%theta_0,ec%types,ec%r3b,stat=epealloc_stat(13))
        ASSERT(epealloc_stat(13).eq.0)
        epealloc_stat(14)=0 ! types
     endif

      if(allocated(ml_dprime_factors) .and. do_deallocate)  then
#ifdef NEW_EPE
       deallocate(ml_dprime_factors,ml_fac,stat=epealloc_stat(7))
#else
       deallocate(ml_dprime_factors,ml_fac,which_epe_ion,stat=epealloc_stat(7))
#endif
       ! which_epe_ion:
       ASSERT(epealloc_stat(7).eq.0)
       epealloc_stat(8)=0 ! ml_dprime_factors ml_fac
       endif

#ifdef NEW_EPE
      deallocate(which_epe_ion,stat=allocation_status)
      ASSERT(allocation_status==0)
#endif

        if(associated(host%b) .and. do_deallocate) then
           deallocate(&
                host%b, host%c, host%ro, host%d, &
                host%k, host%r0, host%k1, host%r1, &
                host%sr1,host%sr2, stat=epealloc_stat(2))
           if(epealloc_stat(2).ne.0)  call error_handler("deallocate host failed")
           epealloc_stat(1)=0 !sr1 sr2
        end if

     if(option_c3_symm) deallocate(c3,permut_ind)
     if(ml_cluster_simulated) then
        deallocate(ml_cluster_simulators,stat=alloc_status)
        if(alloc_status.ne.0) &
             call error_handler("deallocate ml_cluster_simulators failed") 
     end if

     if(allocated(epe_at_start)) then
      deallocate(epe_at_start,stat=epealloc_stat(6))
      ASSERT(epealloc_stat(6).eq.0)
     endif

     do i_alloc=1,size(epealloc_stat)
     if(epealloc_stat(i_alloc).ne.0) print*, 'epealloc_stat ne 0',i_alloc
     enddo
   end subroutine finish_epe

   subroutine gopt_main_action
     
     short_range_grad(1:reg_I_n_ions,:)=zero
     point_1_core_grad(1:reg_I_n_ions,:)=zero
     point_1_shell_grad(1:reg_I_n_ions,:)=zero
     CALL lattice_grad_gopt(r_sh_ion(1,1),r_sh_ion(1,2),r_sh_ion(1,3)  &
          ,r_nuc_ion(1,1),r_nuc_ion(1,2),r_nuc_ion(1,3)  &
          ,point_1_shell_grad(1,1),point_1_shell_grad(1,2),point_1_shell_grad(1,3)      &
          ,point_1_core_grad(1,1), point_1_core_grad(1,2), point_1_core_grad(1,3)  &
          ,short_range_grad(1,1),short_range_grad(1,2),short_range_grad(1,3))
     ETM=energy_shortrange_po(.true.)
     EPL=energy_coulomb_po(.true.)
     write(output_epe,*) 'coulomb energy of lattice ',EPL
     write(output_epe,*) 'short-range energy of lattice ',ETM
     if(n_types_central_atoms_3body>0) then
        E3B=energy_3_body_po()
        write(output_epe,*) 'three-body energy of lattice ',E3B
     else
        E3B=zero
     endif
     etotal_gopt=ETM+EPL+E3B
     write(output_epe,*) 'energy of lattice ', ETM+EPL+E3B
     print*,'energy of lattice ', ETM+EPL+E3B

     do i=1,n_ions_cell
        write(gxcell_unit,203) an(i),r_ion_in_cell(i,:)  &
             ,(indgx(k,i),k=1,8),impu(i)
203     FORMAT((f5.2,3(2X,f13.7),2i4,2x,3I3,2X,3I3,2x,2i3))
     enddo      ! i=1,n_ions_cell
     write(gxcell_unit,203) zero,zero,zero,zero,0,0,0,0,0,0,0,0,0
     write(gxcell_unit,'(2f20.7,3i5)')ETM+EPL+E3B,EPL,81,1,0

     do io=1,n_ions_cell
        ia=which_epe_ion(io)%new
        write(gxcell_unit,'((i5,5x,3f15.6))') io,&
             scaling_factor* &
             (short_range_grad(ia,1)+point_1_shell_grad(ia,1)+point_1_core_grad(ia,1))  &
             ,scaling_factor* &
             (short_range_grad(ia,2)+point_1_shell_grad(ia,2)+point_1_core_grad(ia,2))  &
             ,scaling_factor* &
             (short_range_grad(ia,3)+point_1_shell_grad(ia,3)+point_1_core_grad(ia,3))
     enddo ! i=1,n_ions_cell

     point_1_shell_grad= &
          point_1_shell_grad+short_range_grad

   end subroutine gopt_main_action
!****************************************

!****************************************
   subroutine gopt_displacements
     real(r8_kind)::     displ_vector(3)
     write(output_epe,*) 'displacements'
     do io=1,n_ions_cell
        ia=which_epe_ion(io)%new
        which_epe_ion(io)%old=which_epe_ion(io)%new
        R_ION_IN_CELL_OLD(io,:)=r_ion_in_cell(io,:)
        displ_vector=-0.05*(point_0_shell_grad(Ia,:)+point_0_core_grad(Ia,:)) &
             *(dielectric_const+2.)/(3.*dielectric_const)/(2.*ml_dprime_factors(epe(ia)%k))
        DELS=sqrt(dot_product(displ_vector,displ_vector))
        if(DELS.GT.DMAX) displ_vector=displ_vector*DMAX/DELS
        r_ion_in_cell(io,:)=r_ion_in_cell(io,:)+displ_vector
     enddo ! i=1,n_ions_cell
   end subroutine gopt_displacements
!****************************************

!****************************************
   SUBROUTINE epe_forces(XS,YS,ZS,XC,YC,ZC,DSX,DSY,DSZ, &
        DCX,DCY,DCZ)
     save
! **calculation of relaxation gradients of lattice ions 
! **(DS - shell, DC - core)

     real(kind=r8_kind),intent(out),&
          dimension(reg_I_n_ions) :: DSX,DSY,DSZ,DCX,DCY,DCZ
! real(kind=r8_kind),intent(out),dimension(ndr1) :: DX,DY,DZ !short range contrib.
     real(kind=r8_kind),intent(in),dimension(n_gen_ions) :: XS,YS,ZS,XC,YC,ZC
     real(kind=r8_kind) :: DIST,fF1,ERF1,RRR2,RSSX,RSSY,RSSZ,CSIX,CSIY,CSIZ
     real(kind=r8_kind) :: RSS1,RSS2,RSSR2,RSSR8,RSSD,ESS1,ESS2,RSSRO
     real(kind=r8_kind) :: BRO,OSO,OSS,RSCX,RSCY,RSCZ,RCSX,RCSY,RCSZ,RCCX,RCCY
     real(kind=r8_kind) :: RCCZ,RSRX,RSRY,RSRZ,RCRX,RCRY,RCRZ,CCIX,CCIY,CCIZ
     real(kind=r8_kind) :: RSS3
     real(kind=r8_kind) :: RCS3,RSC3,RCC3,RSR1,RCR1,OSC,OSR,OCS,OCC,OCR,OSSCSI
     real(kind=r8_kind) :: OCSCSI,OSCCCI,OCCCCI
     real(kind=r8_kind) :: RRSX,RRSY,RRSZ,RRCX,RRCY,RRCZ,RRS2,RRC2,RRS1,RRC1
     real(kind=r8_kind) :: OSX,OSY,OSZ,ORS,ORC,OSSX,OSSY,OSSZ,OSCX,OSCY,OSCZ
     real(kind=r8_kind) :: OCSX,OCSY,OCSZ,OCCX,OCCY,OCCZ,CMCSX,CMCSY,CMCSZ
     real(kind=r8_kind) :: CSIS2,CSIC2,ET2S,ET2C,CSIS,OEDS,CSIC,OEDC,GS(3)
     real(kind=r8_kind) :: GCX,GCY,GCZ,GIS2,GIC2,PSIS,PSIC,PI2S,PI2C,OEPKX,OEPKY,OEPKZ
     real(kind=r8_kind),   allocatable, dimension(:) :: DCX_send, DCY_send, DCZ_send, &
          DSX_send, DSY_send, DSZ_send
     integer(kind=i4_kind) :: n_proc, n_forces_total, n_forces_master, n_forces_slave 
     integer(kind=i4_kind) :: f_1, f_2, status,ic

     type epe_vectors
        real(kind=r8_kind),dimension(3)::ss,cs,sc,cc,sr,cr,ds,dc
     end type epe_vectors
     type(epe_vectors) vec
     type epe_dot_products
        real(kind=r8_kind)::ss,cs,sc,cc,sr,cr,ss_ds,cc_dc,sc_dc,cs_ds
     end type epe_dot_products
     type(epe_dot_products)::dp
     type force_contrib
        real(kind=r8_kind),dimension(3)::s,c
     end type force_contrib
     type(force_contrib) f1,f2,f3

     integer(kind=i4_kind) :: i,ig,j,n,k,ia_1,ia_2,ia_3,i1,j1,k1,l
     real(kind=r8_kind) :: e1(3),e2(3),g3b,theta,cos_th,sin_th,scal
     real(kind=r8_kind) :: deg2rad

     epe(:)%s(1)=r_sh_ion(:,1)
     epe(:)%s(2)=r_sh_ion(:,2)
     epe(:)%s(3)=r_sh_ion(:,3)
     epe(:)%c(1)=r_nuc_ion(:,1)
     epe(:)%c(2)=r_nuc_ion(:,2)
     epe(:)%c(3)=r_nuc_ion(:,3)

     DIST=7.0_r8_kind
     ff1=0.026_r8_kind
     ERF1=4.0_r8_kind/3.0_r8_kind*ERROR_FUNCTION_PARAMETER**3/PIS
     short_range_grad(:,:)=zero
     DSX=zero
     DSY=zero
     DSZ=zero
     DCX=zero
     DCY=zero
     DCZ=zero
!AG ===================================================================
!AG                        Distribute force calculations
     n_proc=comm_get_n_processors()
     n_forces_total  = reg_I_n_ions
     n_forces_master = n_forces_total - &
          (n_forces_total/n_proc)*(n_proc-1)

     if(n_proc > 1) then
        n_forces_slave = (n_forces_total-n_forces_master)/(n_proc-1)
     else
        n_forces_slave=0
     end if
!                      forces_s=>( . . . f_1........f_2 . . . )     
     if( comm_i_am_master() ) then
        f_1=( comm_get_n_processors()-1 )*n_forces_slave + 1
        f_2=f_1 + n_forces_master - 1 
     else 
        f_1=( comm_myindex()-2 )*n_forces_slave + 1
        f_2=f_1 + n_forces_slave - 1
        if(.not.allocated(DCX_send))then
           allocate(DCX_send(reg_I_n_ions), &
                DCY_send(reg_I_n_ions), &
                DCZ_send(reg_I_n_ions), &
                DSX_send(reg_I_n_ions), & 
                DSY_send(reg_I_n_ions), &
                DSZ_send(reg_I_n_ions) ,stat=status )
           if(status.ne.0) call error_handler("allocate DCX_send failed")
        end if
    
        DSX_send=zero
        DSY_send=zero
        DSZ_send=zero
        DCX_send=zero
        DCY_send=zero
        DCZ_send=zero
     end if

! ** contribution to gradients due to 3-body interaction into I region
     if(n_types_central_atoms_3body > 0 ) then
        deg2rad=pi/180.0_r8_kind
!!$        fst: do i=1,last_ind-first_ind+1
        fst: do i=first_ind,last_ind
           do l=1,5
              if (tetra_atoms(i,l) <= reg_I_n_ions) goto 1
           enddo
           cycle fst
1          ia_1=tetra_atoms(i,1)
           i1=epe(ia_1)%k
           scnd :do j=2,4
              ia_2=tetra_atoms(i,j)
              if(ia_2 == 0) cycle scnd
              j1=epe(ia_2)%k
              thrd: do k=j+1,5
                 ia_3=tetra_atoms(i,k)
                 if(ia_3 == 0) cycle thrd
                 if(ia_1 > reg_I_n_ions.and.ia_2 > reg_I_n_ions.and.ia_3 > reg_I_n_ions) cycle thrd
                 k1=epe(ia_3)%k
                 rss1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))
                 rss2=dot_product(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
                 scal=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
                 rss1=sqrt(rss1)
                 rss2=sqrt(rss2)
                 cos_th=scal/(rss1*rss2)
                 theta=acos(cos_th)
                 sin_th=sin(theta)
                 g3b=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
                 e1=(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))/rss1
                 e2=(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))/rss2
                 if(ia_3 <= reg_I_n_ions) &
                      short_range_grad(ia_3,:)=short_range_grad(ia_3,:)+ & 
                      g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                 if(ia_2 <= reg_I_n_ions) &
                 short_range_grad(ia_2,:)=short_range_grad(ia_2,:)+ & 
                      g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
                 if(ia_1 <= reg_I_n_ions) &
                 short_range_grad(ia_1,:)=short_range_grad(ia_1,:)+ & 
                      g3b*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/ &
                      (rss2*rss1*sin_th)
              enddo thrd!k=j+1,4
           enddo scnd!j=1,3
        enddo fst!i=1,n_tetra_atoms
     endif

! ** area 1 - start
!AG  DO I=1,reg_I_n_ions (replaced)
     do I=f_1, f_2
        K=epe(I)%k
! ** area 2a - start
        DO J=reg_I_n_ions+1,reg_2a_n_ions       ! contrib due to reg. II discplacements

           RRR2=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)
           IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle
           N=epe(J)%k
           vec%ss=epe(I)%s-epe(J)%s
           vec%ds=epe(J)%s-epe(J)%r
           dp%ss=dot_product(vec%ss,vec%ss)
           dp%ss_ds=dot_product(vec%ss,vec%ds)
           RSS1=SQRT(dp%ss)
           if(RRR2.GT.16.0_r8_kind) then
              ic=1
           else
              ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
              if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
           endif

           IF(RRR2.le.host%sr1(K,N,ic)**2) then
              RSSR2=1.0_r8_kind/dp%ss
              RSSR8=(RSSR2**2)**2
              RSSD=RSSR2*host%D(K,N,ic)
              ESS1=(6.0_r8_kind*host%C(K,N,ic) + 8.0_r8_kind*RSSD)*RSSR8
              ESS2=(48.0_r8_kind*host%C(K,N,ic) +80.0_r8_kind*RSSD)*(RSSR8*RSSR2)

              IF(RRR2.le.host%sr2(K,N,ic)**2) then
                 RSSRO=1.0_r8_kind/(host%RO(K,N,ic)*RSS1)
                 BRO=host%B(K,N,ic)*RSSRO*EXP(-RSS1/host%RO(K,N,ic))
                 ESS1=-BRO+ESS1
                 ESS2=BRO*(RSSR2+RSSRO)-ESS2
              endif! 

              ESS2=0.5_r8_kind*ESS2*dp%ss_ds
              OSO=0.5_r8_kind/dp%ss
              RSSR8=(RSSR2**2)**2
              RSSD=RSSR2*host%D(K,N,ic)
              ESS1=(6.0_r8_kind*host%C(K,N,ic) + 8.0_r8_kind*RSSD)*RSSR8
              ESS2=(48.0_r8_kind*host%C(K,N,ic) +80.0_r8_kind*RSSD)*(RSSR8*RSSR2)
              IF(RRR2.le.host%sr2(K,N,ic)**2) then
                 RSSRO=1.0_r8_kind/(host%RO(K,N,ic)*RSS1)
                 BRO=host%B(K,N,ic)*RSSRO*EXP(-RSS1/host%RO(K,N,ic))
                 ESS1=-BRO+ESS1
                 ESS2=BRO*(RSSR2+RSSRO)-ESS2
              endif! 
              ESS2=0.5_r8_kind*ESS2*dp%ss_ds
              OSO=0.5_r8_kind*ESS1
              OSS=ESS1+ESS2
              short_range_grad(i,:)=short_range_grad(i,:)+OSS*vec%ss+OSO*vec%ds
           endif! (RRR2.le.T
           ! r in a.u , energy in eV
           if(host%k(K,N,ic).ne.0.0_r8_kind.and. &
              abs(rss1-host%r0(K,N,ic)/auangs).lt.0.1) then
!           ess1=2
!               eau_ev=27.21.....
           endif

           vec%sc=epe(i)%s-epe(j)%c
           vec%cs=epe(i)%c-epe(j)%s
           vec%cc=epe(i)%c-epe(j)%c
           vec%sr=epe(i)%s-epe(j)%r
           vec%cr=epe(i)%c-epe(j)%r
           vec%dc=epe(j)%c-epe(j)%r

           dp%ss=dot_product(vec%ss,vec%ss)
           dp%sc=dot_product(vec%sc,vec%sc)
           dp%cs=dot_product(vec%cs,vec%cs)
           dp%cc=dot_product(vec%cc,vec%cc)
           dp%sr=dot_product(vec%sr,vec%sr)
           dp%cr=dot_product(vec%cr,vec%cr)
           dp%cs_ds=dot_product(vec%cs,vec%ds)
           dp%sc_dc=dot_product(vec%sc,vec%dc)
           dp%cc_dc=dot_product(vec%cc,vec%dc)

           RSS3=dp%ss*RSS1
           RCS3=dp%cs*SQRT(dp%cs)
           RSC3=dp%sc*SQRT(dp%sc)
           RCC3=dp%cc*SQRT(dp%cc)
           RSR1=SQRT(dp%sr)
           RCR1=SQRT(dp%cr)

           OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
           OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
           OSR=Q_SHELL(K)*Q_ion(epe(J)%k)*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*dp%sr))/dp%sr
           OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
           OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
           OCR=Q_NUCLEAR(K)*Q_ion(epe(J)%k)*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*dp%cr))/dp%cr
           OSSCSI=-3.0_r8_kind*OSS*dp%ss_ds/dp%ss
           OCSCSI=-3.0_r8_kind*OCS*dp%cs_ds/dp%cs
           OSCCCI=-3.0_r8_kind*OSC*dp%sc_dc/dp%sc
           OCCCCI=-3.0_r8_kind*OCC*dp%cc_dc/dp%cc

           f1%s=OSS*vec%ss+OSC*vec%sc+OSR*vec%sr
           f1%c=OCC*vec%cc+OCS*vec%cs+OCR*vec%cr
           f2%s=OSS*vec%ds+OSC*vec%dc+OSCCCI*vec%sc+OSSCSI*vec%ss
           f2%c=OCC*vec%dc+OCS*vec%ds+OCSCSI*vec%cs+OCCCCI*vec%cc
           DSX(I)=DSX(I)+f1%s(1)+0.5_r8_kind*f2%s(1)
           DSy(I)=DSy(I)+f1%s(2)+0.5_r8_kind*f2%s(2)
           DSz(I)=DSz(I)+f1%s(3)+0.5_r8_kind*f2%s(3)
           DCx(I)=DCx(I)+f1%c(1)+0.5_r8_kind*f2%c(1)
           DCy(I)=DCy(I)+f1%c(2)+0.5_r8_kind*f2%c(2)
           DCz(I)=DCz(I)+f1%c(3)+0.5_r8_kind*f2%c(3)
        enddo
! **done, contrib. due to region 2a 

! **treat contrib. from polarized region I
        DO J=I+1,reg_I_n_ions
           N=epe(J)%k
           
           RSSX=XS(I)-XS(J)
           RSSY=YS(I)-YS(J)
           RSSZ=ZS(I)-ZS(J)
           vec%ss=epe(i)%s-epe(j)%s

           RSCX=XS(I)-XC(J)
           RSCY=YS(I)-YC(J)
           RSCZ=ZS(I)-ZC(J)

           RCSX=XC(I)-XS(J)
           RCSY=YC(I)-YS(J)
           RCSZ=ZC(I)-ZS(J)

           RCCX=XC(I)-XC(J)
           RCCY=YC(I)-YC(J)
           RCCZ=ZC(I)-ZC(J)
 
           RSRX=XS(I)-epe(j)%r(1)
           RSRY=YS(I)-epe(j)%r(2)
           RSRZ=ZS(I)-epe(j)%r(3)

           RCRX=XC(I)-epe(j)%r(1)
           RCRY=YC(I)-epe(j)%r(2)
           RCRZ=ZC(I)-epe(j)%r(3)

           RRSX=epe(I)%r(1)-XS(J)
           RRSY=epe(I)%r(2)-YS(J)
           RRSZ=epe(I)%r(3)-ZS(J)
           
           RRCX=epe(I)%r(1)-XC(J)
           RRCY=epe(I)%r(2)-YC(J)
           RRCZ=epe(I)%r(3)-ZC(J)

           RRR2=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)
           dp%ss=dot_product(vec%ss,vec%ss)
           dp%sc=RSCX**2+RSCY**2+RSCZ**2
           dp%cs=RCSX**2+RCSY**2+RCSZ**2
           dp%cc=RCCX**2+RCCY**2+RCCZ**2
           dp%sr=RSRX**2+RSRY**2+RSRZ**2
           dp%cr=RCRX**2+RCRY**2+RCRZ**2
           RRS2=RRSX**2+RRSY**2+RRSZ**2
           RRC2=RRCX**2+RRCY**2+RRCZ**2

           RSS1=SQRT(dp%ss)
           RSS3=dp%ss*RSS1
           RCC3=dp%cc*SQRT(dp%cc)
           RSC3=dp%sc*SQRT(dp%sc)
           RCS3=dp%cs*SQRT(dp%cs)
           RSR1=SQRT(dp%sr)
           RCR1=SQRT(dp%cr)
           RRS1=SQRT(RRS2)
           RRC1=SQRT(RRC2)

        ! **short region interacting part of area 1 - start
           if(RRR2.GT.16.0_r8_kind) then
              ic=1
           else
              ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
              if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
           endif
           IF(RRR2.le.host%sr1(K,N,ic)**2) then
              RSSR2=1.0_r8_kind/dp%ss
              ESS1=(6.0_r8_kind*host%C(K,N,ic)+8.0_r8_kind*host%D(K,N,ic)*RSSR2)*((RSSR2**2)**2)
              IF(RRR2.LE.host%sr1(K,N,ic)**2) &
                   ESS1=ESS1-host%B(K,N,ic)/(host%RO(K,N,ic)*RSS1)*EXP(-RSS1/host%RO(K,N,ic))
              OSX=ESS1*RSSX
              OSY=ESS1*RSSY
              OSZ=ESS1*RSSZ
              f3%s=ESS1*vec%ss
              short_range_grad(i,:)=short_range_grad(i,:)+f3%s
              short_range_grad(j,:)=short_range_grad(j,:)-f3%s
           endif
        ! **short region interacting part of area 1 - end

           OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
           OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
           OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
           OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
           IF(ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST) then  
              OSR=Q_SHELL(K)*Q_ion(epe(J)%k)*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*dp%sr))/dp%sr
              OCR=Q_NUCLEAR(K)*Q_ion(epe(J)%k)*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*dp%cr))/dp%cr
              ORS=(Q_ion(epe(i)%k)*Q_SHELL(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRS1)/RRS1+ERFO*EXP(-ET2*RRS2))/RRS2
              ORC=(Q_ion(epe(i)%k)*Q_NUCLEAR(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRC1)/RRC1+ERFO*EXP(-ET2*RRC2))/RRC2
           else
              OSR=Q_SHELL(K)*Q_ion(epe(J)%k)/(dp%sr*RSR1)
              OCR=Q_NUCLEAR(K)*Q_ion(epe(J)%k)/(dp%cr*RCR1)
              ORS=Q_ion(epe(i)%k)*Q_SHELL(N)/(RRS2*RRS1)
              ORC=Q_ion(epe(i)%k)*Q_NUCLEAR(N)/(RRC2*RRC1)
           endif!  ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST/else

           OSSX=OSS*RSSX
           OSSY=OSS*RSSY
           OSSZ=OSS*RSSZ
           OSCX=OSC*RSCX
           OSCY=OSC*RSCY
           OSCZ=OSC*RSCZ
           OCSX=OCS*RCSX
           OCSY=OCS*RCSY
           OCSZ=OCS*RCSZ
           OCCX=OCC*RCCX
           OCCY=OCC*RCCY
           OCCZ=OCC*RCCZ
           DSX(I)=DSX(I)+OSSX+OSCX+OSR*RSRX
           DSY(I)=DSY(I)+OSSY+OSCY+OSR*RSRY
           DSZ(I)=DSZ(I)+OSSZ+OSCZ+OSR*RSRZ
           DSX(J)=DSX(J)-OSSX-OCSX-ORS*RRSX
           DSY(J)=DSY(J)-OSSY-OCSY-ORS*RRSY
           DSZ(J)=DSZ(J)-OSSZ-OCSZ-ORS*RRSZ
           DCX(I)=DCX(I)+OCSX+OCCX+OCR*RCRX
           DCY(I)=DCY(I)+OCSY+OCCY+OCR*RCRY
           DCZ(I)=DCZ(I)+OCSZ+OCCZ+OCR*RCRZ
           DCX(J)=DCX(J)-OSCX-OCCX-ORC*RRCX
           DCY(J)=DCY(J)-OSCY-OCCY-ORC*RRCY
           DCZ(J)=DCZ(J)-OSCZ-OCCZ-ORC*RRCZ
        !AG[
           if(.not.comm_i_am_master() ) then
              DSX_send(J)=DSX(J)
              DSY_send(J)=DSY(J)
              DSZ_send(J)=DSZ(J)
              DCX_send(J)=DCX(J)
              DCY_send(J)=DCY(J)
              DCZ_send(J)=DCZ(J)
           end if
        !AG]
        enddo

        ! **polarized area 1 - end
        
        CSIX=XS(I)-epe(I)%r(1)
        CSIY=YS(I)-epe(I)%r(2)
        CSIZ=ZS(I)-epe(I)%r(3)
        CCIX=XC(I)-epe(I)%r(1)
        CCIY=YC(I)-epe(I)%r(2)
        CCIZ=ZC(I)-epe(I)%r(3)
        CMCSX=XC(I)-XS(I)
        CMCSY=YC(I)-YS(I)
        CMCSZ=ZC(I)-ZS(I)
        CSIS2=CSIX**2+CSIY**2+CSIZ**2
        CSIC2=CCIX**2+CCIY**2+CCIZ**2
        ET2S=ET2*CSIS2
        ET2C=ET2*CSIC2

        IF(ET2S.GT.ff1) then
           CSIS=SQRT(CSIS2)
           OEDS=(ERF(ERROR_FUNCTION_PARAMETER*CSIS)/CSIS+ERFO*EXP(-ET2S))/CSIS2
        else
           OEDS=ERF1*(1.-ET2S*(0.6-0.2142856*ET2S))
        endif! ET2S.GT.ff1)/else
    
        IF(ET2C.GT.ff1) then
           CSIC=SQRT(CSIC2)
           OEDC=(ERF(ERROR_FUNCTION_PARAMETER*CSIC)/CSIC+ERFO*EXP(-ET2C))/CSIC2
        else
           OEDC=ERF1*(1.-ET2C*(0.6-0.2142856*ET2C))
        endif! ET2C.GT.ff1

        OEDS=Q_SHELL(K)*Q_ion(epe(i)%k)*OEDS
        OEDC=Q_NUCLEAR(K)*Q_ion(epe(i)%k)*OEDC
        GS=zero
        GCX=zero
        GCY=zero
        GCZ=zero
        DO IG=1,n_bs_points
           GIS2=PI2*(GSTR(IG,1)*XS(I)+GSTR(IG,2)*YS(I)+GSTR(IG,3)*ZS(I))
           GIC2=PI2*(GSTR(IG,1)*XC(I)+GSTR(IG,2)*YC(I)+GSTR(IG,3)*ZC(I))
           PSIS=RSIN(IG)*COS(GIS2)-RCOS(IG)*SIN(GIS2)
           PSIC=RSIN(IG)*COS(GIC2)-RCOS(IG)*SIN(GIC2)
           GS=GS+PSIS*GSTR(IG,:)
           GCX=GCX+PSIC*GSTR(IG,1)
           GCY=GCY+PSIC*GSTR(IG,2)
           GCZ=GCZ+PSIC*GSTR(IG,3)
        enddo!IG=1,n_bs_points

        PI2S=PI2*Q_SHELL(K)
        PI2C=PI2*Q_NUCLEAR(K)
        GS=GS*PI2S
        GCX=GCX*PI2C
        GCY=GCY*PI2C
        GCZ=GCZ*PI2C
        OEPKX=PK(K)*CMCSX
        OEPKY=PK(K)*CMCSY
        OEPKZ=PK(K)*CMCSZ
        DSX(I)=DSX(I)+CSIX*OEDS+GS(1)-OEPKX
        DSY(I)=DSY(I)+CSIY*OEDS+GS(2)-OEPKY
        DSZ(I)=DSZ(I)+CSIZ*OEDS+GS(3)-OEPKZ
        DCX(I)=DCX(I)+CCIX*OEDC+GCX+OEPKX
        DCY(I)=DCY(I)+CCIY*OEDC+GCY+OEPKY
        DCZ(I)=DCZ(I)+CCIZ*OEDC+GCZ+OEPKZ

!!$    write(output_epe,*) i,max(DSX(I),DSy(I),Dsz(I))
     enddo!reg_I_n_ions
!!$  write(output_epe,*) short_range_grad(1,:),' reg_I contrib',radius_first_short_interact,host%sr1
!AG=========================================================================
!AG                  Prepare array to send forces to master (slave only)
     if( .not.comm_i_am_master() ) then
        call comm_init_send(comm_master_host, msgtag_epe_forces_done)

        DCX_send(f_1:f_2)=DCX(f_1:f_2)
        DCY_send(f_1:f_2)=DCY(f_1:f_2)
        DCZ_send(f_1:f_2)=DCZ(f_1:f_2)
        DSX_send(f_1:f_2)=DSX(f_1:f_2)
        DSY_send(f_1:f_2)=DSY(f_1:f_2)
        DSZ_send(f_1:f_2)=DSZ(f_1:f_2)
!AG                         Now send forces for master to gather
        call commpack(DCX_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [1]")
        call commpack(DCY_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [2]")
        call commpack(DCZ_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [3]")
        call commpack(DSX_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [4]")
        call commpack(DSY_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [5]")
        call commpack(DSZ_send(1), reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [6]")
        call commpack(short_range_grad(1,1), 3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_forces:  error [7]")
        call comm_send()
        deallocate(DCX_send,DCY_send,DCZ_send,DSX_send,DSY_send,DSZ_send ,stat=status)  
        if(status.ne.0) call error_handler("deallocate DCX_send failed")
     end if! (slave work)
!AG =========================================================================
!!$call cpu_time(t2)
!!$print*,'>>>>>>>>>> epe_forces ',t2-t1
   END SUBROUTINE epe_forces
!***************************************************

!***************************************************
   subroutine symmetrize_forces(init_forces)
     real(kind=r8_kind), dimension(:,:),intent(inout):: init_forces
     real(kind=r8_kind),dimension(reg_I_n_ions,3):: symmetrized_forces, &
          first_vec(3),second_vec(3)
     symmetrized_forces=zero
     do i=1,reg_I_n_ions
        first_vec=init_forces(i,:)
!if(i.le.10)    print*,'f',i,dot_product(first_vec,first_vec)
        do j=1,3
           second_vec=matmul(c3(j)%rotmat,first_vec)
           symmetrized_forces(permut_ind(i,j),:)=symmetrized_forces(permut_ind(i,j),:) &
                + second_vec/3.0_r8_kind
        enddo
     enddo
     init_forces=symmetrized_forces
   end subroutine symmetrize_forces

   subroutine calc_def_contributions(point_x_core_grad,point_x_shell_grad, &
        short_range_grad1,impurity_core_grad,impurity_shell_grad,use_ref_data)

! **  energies and gradients due to 
! **  defect - defect and defect - lattice interactions
! **a) all impurities are atoms of the cluster

     use epecom_module
     use culon_module

     logical, intent(in):: use_ref_data 
     real(kind=r8_kind),intent(inout),dimension(:,:) ::   & !(reg_I_n_ions,3)
          point_x_core_grad,point_x_shell_grad                                  
     real(kind=r8_kind),intent(inout),dimension(:,:) :: short_range_grad1
     real(kind=r8_kind),intent(inout),dimension(:,:) :: impurity_core_grad,impurity_shell_grad !(n_impurities,3)
     real(kind=r8_kind), allocatable, dimension(:,:) :: send_intermediate_array
     real(kind=r8_kind), allocatable, dimension(:,:) :: impurity_core_grad_send,  &
                                                     impurity_shell_grad_send
     real(kind=r8_kind), allocatable, dimension(:)   :: reg_2a_ec, &
                                                     reg_2a_ec_send
     real(kind=r8_kind) :: vs_ewaback,vc_ewaback
     real(kind=r8_kind) :: eshort_lat_relaxation=0.d0, &
          ecoul_vaccluster_epe,ec_cluster_epe
     real(kind=r8_kind) :: ec_reg2a,ec_current,ec_reg2a_ref,eshort_epecluster_epe
     real(kind=r8_kind) :: angle_fac,gc_ewaback(3),gs_ewaback(3)
     real(kind=r8_kind) :: ec_ewadir_cluster, eshort_vacepe, &
        e_vac_3b=0.0_r8_kind, &
        e_imp_3b=0.0_r8_kind ! contrib to eshort_epecluster
     integer(kind=i4_kind) :: i,ind,j,k,ii1,jj1
     integer(kind=i4_kind) :: n_proc, i_1, i_2, status
     integer(kind=i4_kind) :: n_vacancies_master,  n_vacancies_slave
     integer(kind=i4_kind) :: n_impurities_master, n_impurities_slave  
     type epe_dot_products
        real(kind=r8_kind)::ss,sq_ss,sq_sc,sq_cs,sq_cc,ss3,ss6,ss8,ss10 &
             ,cs,cc,sc,sr,cr,rr,rc,rs_sr,rc_cr,ss_sr,sc_cr,cc_cr,cs_sr &
             ,rs
     end type epe_dot_products
     type factors
        real(kind=r8_kind)::ss,sc,cs,cc,sewa,cewa,gsewa,gcewa, &
             rs,rc
     end type factors
     real(kind=r8_kind),parameter::third=0.33333333_r8_kind,tenth=0.1_r8_kind,one=1.0_r8_kind
     real(kind=r8_kind)::serfunc_arg(3),cerfunc_arg(3),gewa_const
     type(epe_dot_products),target::reg,ref,var
     type(epe_dot_products)::sim
     type(epe_dot_products),pointer::l
     type each_center_contributions
        real(kind=r8_kind) coul,short,tot,regI,reg2a,eshort,gshort
        type(factors):: fac
     end type each_center_contributions
     type(each_center_contributions),allocatable,dimension(:)::clus,vaca
     type(each_center_contributions),allocatable,dimension(:)::ml_simulator
     real(kind=r8_kind)::ecoul_lat_relcontrib,eshort_lat_relcontrib
     real(kind=r8_kind), parameter::tail_exp_dist=7.0_r8_kind, &
          erfarg_th=0.035_r8_kind,gerfarg_th=0.026_r8_kind
     type epe_imp_tet
        real(kind=r8_kind) :: coord(5,3)
        character(len=3) :: center(5)
        integer(kind=i4_kind) :: number(5)
     end type epe_imp_tet
     real(r8_kind) :: short_forces(3)
     real(r8_kind) :: e_imp_3bT
     short_forces=0.0_r8_kind

!!$call cpu_time(t1)
!!$    print*, '=> calc_def_contributions', comm_myindex(),t1
 !print*, 'epe(:)%s='
     epe(:)%s(1)=R_sh_ION(:,1)
     epe(:)%s(2)=R_sh_ION(:,2)
     epe(:)%s(3)=R_sh_ION(:,3)
     epe(:)%c(1)=R_nuc_ION(:,1)
     epe(:)%c(2)=R_nuc_ION(:,2)
     epe(:)%c(3)=R_nuc_ION(:,3)

! **E - energy, V - vacancy, P - impurity, R - lattice, O - shortaction
! **C - coulomb

  !    print*,'init'
     eshort_vaccluster=zero
     eshort_vacepe=zero
     eshort_epecluster=zero
     eshort_epecluster_epe=zero
     eshort_lat_relaxation=zero
     eshort_lat_relcorr=zero
     eshort_reg_I=zero
     eshort_reg2a=zero
     e_imp_3b=zero

     ecoul_vaccluster=zero
     ecoul_vaccluster_epe=zero
     ecoul_epecluster=zero
     ec_cluster_epe=zero

     ecoul_lat_relaxation=zero
     ecoul_lat_relcorr=zero

     ec_cluster_reg2a=zero

     gewa_const=4.0_r8_kind/3.0_r8_kind*ERROR_FUNCTION_PARAMETER**3/PIS

   !   print*,size(impurity_core_grad), size(impurity_shell_grad),N_IMPURITIES
     if(N_IMPURITIES.gt.0) then
      impurity_core_grad(1:N_IMPURITIES,1:3)=zero
      impurity_shell_grad(1:N_IMPURITIES,1:3)=zero
     endif

     IF(N_VACANCIES.ne.0.AND.N_IMPURITIES.ne.0 .and.comm_i_am_master() ) then
        WRITE(output_epe, 108)
108     format(1x,'      energy  /ev/             ', &
             '      coulomn      short range      total')
        WRITE(output_epe, 106)
106     FORMAT(1X,78('-'))
     endif
!   print*,'  **vacancies start monitoring'
     allocate(vaca(N_VACANCIES), stat=status)
     if(status.ne.0 ) call error_handler("allocate vaca failed")
     if(n_vacancies.gt.0) then
      vaca(1:n_vacancies)%coul  = zero ! ?
      vaca(1:n_vacancies)%short = zero
      vaca(1:n_vacancies)%tot   = zero
     endif
![AG                       go in parallel : for each vacancy
     n_proc=comm_get_n_processors()
     if(n_proc > n_vacancies) then
        n_vacancies_master = 1
        n_vacancies_slave = 1
     else
        n_vacancies_master = n_vacancies - &
             (n_vacancies/n_proc)*(n_proc-1)
        if(n_proc > 1) then
           n_vacancies_slave = (n_vacancies-n_vacancies_master)/(n_proc-1)
        else
           n_vacancies_slave=0
        end if
     endif

     if(comm_i_am_master() ) then
        i_1 = ( comm_get_n_processors() - 1 )*n_vacancies_slave +1
        i_2 = i_1 + n_vacancies_master - 1
     else
        i_1 = 1 + ( comm_myindex() - 2 )*n_vacancies_slave
        i_2 = i_1 + n_vacancies_slave - 1
        point_x_core_grad   = zero
        point_x_shell_grad  = zero
        short_range_grad1    = zero
        ec_cluster_epe      = zero 
        ec_cluster_reg2a    = zero
        diffpg_ec_ecref     = zero
        defect_energy_short = zero
        defect_energy_coul  = zero
        ecrr                = zero
     end if
!AG]

     if(i_1 <= n_vacancies) then
        vac: do I=i_1, i_2  !  changed <= vac:    DO I=1,N_VACANCIES
           k=epe(I)%k
           vaca(i)%reg2a=zero
    
!!$   vaca(i)%regi=-madc%emd(epe(i)%m)
   ! for vacancy in the shifted position one need recalculate
   ! Madelung energy contribution


           vaca(i)%regi=zero

!           print*,'direct space ewald contributions'
           do j=1,reg_2a_n_ions

              var%sr=dot_product(epe(i)%s-epe(j)%r,epe(i)%s-epe(j)%r)
              var%cr=dot_product(epe(i)%c-epe(j)%r,epe(i)%c-epe(j)%r)
              
              serfunc_arg(1)=error_function_parameter*sqrt(var%sr)
              serfunc_arg(2)=serfunc_arg(1)**2
              serfunc_arg(3)=serfunc_arg(2)**2
              cerfunc_arg(1)=error_function_parameter*sqrt(var%cr)
              cerfunc_arg(2)=cerfunc_arg(1)**2
              cerfunc_arg(3)=cerfunc_arg(2)**2

              vaca(i)%fac%sewa=erfo*(one-third*serfunc_arg(2)+tenth*serfunc_arg(3))
              vaca(i)%fac%cewa=erfo*(one-third*cerfunc_arg(2)+tenth*cerfunc_arg(3))
              if(serfunc_arg(2).gt.erfarg_th) vaca(i)%fac%sewa=-erf(serfunc_arg(1))/sqrt(var%sr)
              if(cerfunc_arg(2).gt.erfarg_th) vaca(i)%fac%cewa=-erf(cerfunc_arg(1))/sqrt(var%cr)

              vaca(i)%regI=vaca(i)%regI - q_ion(epe(j)%k)* &
                   (q_shell(epe(i)%k)*vaca(i)%fac%sewa+ &
                   q_nuclear(epe(i)%k)*vaca(i)%fac%cewa) 
           end do
!           print*,'back space ewald contributions'
           vs_ewaback=zero
           do ind=1,n_bs_points
              angle_fac=dot_product(gstr(ind,:),epe(i)%s)*pi2
              vs_ewaback=vs_ewaback+rsin(ind)*sin(angle_fac)+rcos(ind)*cos(angle_fac)
           enddo!i=1,n_bs_points
           vaca(i)%regi=vaca(i)%regi-q_shell(epe(i)%k)*vs_ewaback

           vc_ewaback=zero
           do ind=1,n_bs_points
              angle_fac=dot_product(gstr(ind,:),epe(i)%c)*pi2
              vc_ewaback=vc_ewaback+rsin(ind)*sin(angle_fac)+rcos(ind)*cos(angle_fac)
           enddo!i=1,n_bs_points
           vaca(i)%regi=vaca(i)%regi-q_nuclear(epe(i)%k)*vc_ewaback   

           ! shell model intra ionic energy contrib (4)
           vaca(i)%regi=vaca(i)%regi - pk(epe(i)%k)* & 
                dot_product(epe(i)%c-epe(i)%s,epe(i)%c-epe(i)%s)/2.0_r8_kind

           call vac_to_vac_contibs !in calc_def_cont

           ! **vacancy - lattice 
           ! **treat contibutions from region I 
 
           call vac_to_latt_contribs  !in calc_def_cont
           vaca(i)%coul=vaca(i)%reg2a+vaca(i)%regI
           vaca(i)%tot=vaca(i)%short+vaca(i)%coul
           eshort_vacepe=eshort_vacepe+vaca(i)%short
           ecoul_vaccluster_epe=ecoul_vaccluster_epe+vaca(i)%coul
        enddo vac
     endif

     if(n_types_central_atoms_3body > 0 .and. comm_i_am_master()) then
        call vac_3_body(e_vac_3b)
        eshort_vacepe=eshort_vacepe+e_vac_3b  !!  ??
     end if

!AG[
     if(comm_parallel() .and. n_vacancies > 0) then
        if( comm_i_am_master() ) then
           call receive_intermediate_vac() 
        else
           call send_intermediate_vac()
        end if
     end if
     
!AG]
     if(comm_i_am_master() ) then
        do I=1,n_vacancies
           WRITE(output_epe,102) I, vaca(i)%coul, vaca(i)%short, vaca(i)%tot
102        FORMAT(1X,'  vacancy  ',I2,' - lattice        !',F16.8, &
                ' !',F16.8,' !',F16.8)
        end do
        if(n_vacancies > 0) then
           WRITE(output_epe,*) 'sum_over_vacancies_short ',sum(vaca(:)%short)
           WRITE(output_epe,*) 'sum_over_vacancies_3body ',e_vac_3b
           WRITE(output_epe,*) 'sum_over_vacancies_coulomb ', sum( vaca(:)%coul)
           WRITE(output_epe,*) 'sum_over_vacancies_total ',sum( vaca(:)%tot)+e_vac_3b
        endif
     end if! master print only
     deallocate(vaca,stat=status)
     if(status.ne.0)  call error_handler(" deallocate vaca failed")

! **start treatment of impurities
     diffpg_ec_ecref=zero
     diff_reg2a_ec_ecref=zero

    imp_present: if(n_impurities.ne.0) then
     if(.not.allocated(clus)) then
     allocate(clus(N_IMPURITIES),stat=status)
     if(status.ne.0)   call error_handler("main_epe allocate clus failed")
     endif
     if(.not.allocated(reg_2a_ec)) then
     allocate(reg_2a_ec(n_impurities),stat=status)
     if(status.ne.0)   call error_handler("main_epe allocate reg_2a_ec failed")
     endif
     clus(1:n_impurities)%coul  = zero
     clus(1:n_impurities)%short = zero
     clus(1:n_impurities)%tot   = zero 
     reg_2a_ec(:)               = zero
![AG                   go in parallel : for each impurity
     if(n_proc > n_impurities) then
        n_impurities_master = 1
        n_impurities_slave = 1
     else
        n_impurities_master = n_impurities - &
             (n_impurities/n_proc)*(n_proc-1)
        if(n_proc > 1) then
           n_impurities_slave = (n_impurities-n_impurities_master)/(n_proc-1)
        else
           n_impurities_slave=0
        end if
        if(comm_i_am_master() ) then
           i_1 = ( comm_get_n_processors() - 1 )*n_impurities_slave + 1
           i_2 = i_1 + n_impurities_master - 1
        else
           i_1 = 1 + ( comm_myindex() - 2 )*n_impurities_slave
           i_2 = i_1 + n_impurities_slave - 1
        end if
     endif
!AG]

     i_1isimp: if(i_1 <= n_impurities) then
     imp2rest: do I=i_1, i_2
        K=TYPE_IMPURITY(I)
        DO J=I+1,N_IMPURITIES
           ! **coulomb interaction of imurities : energies and gradients

              ref_data1: if(use_ref_data) then !only regular contribution (a constant)

                 ii1=imp2center(i) 
                 jj1=imp2center(j)
                 reg%ss=dot_product(reg_reference(i)%rs-reg_reference(j)%rs, &
                      reg_reference(i)%rs-reg_reference(j)%rs)
                 reg%cs=dot_product(reg_reference(i)%rc-reg_reference(j)%rs, &
                      reg_reference(i)%rc-reg_reference(j)%rs)
                 reg%sc=dot_product(reg_reference(i)%rs-reg_reference(j)%rc, &
                      reg_reference(i)%rs-reg_reference(j)%rc)
                 reg%cc=dot_product(reg_reference(i)%rc-reg_reference(j)%rc, &
                      reg_reference(i)%rc-reg_reference(j)%rc)
                 reg%rr=dot_product(R_IMP(I,:)-R_IMP(J,:),R_IMP(I,:)-R_IMP(J,:))

                 clus(i)%coul=Q_SH_IMPURITY(I)*Q_SH_IMPURITY(J)/sqrt(reg%ss) &
                      +Q_SH_IMPURITY(I)*Q_NUC_IMPURITY(J)/sqrt(reg%sc) &
                      +Q_NUC_IMPURITY(I)*Q_SH_IMPURITY(J)/sqrt(reg%cs) &
                      +Q_NUC_IMPURITY(I)*Q_NUC_IMPURITY(J)/sqrt(reg%cc)

                 l=>reg ! rr calculated
                 CALL calc_shortrange(TYPE_IMPURITY(I),TYPE_IMPURITY(J),ii1,jj1, & ! def contr 1
                      clus(i)%short,clus(i)%gshort)

              else ref_data1 ! calculate energy of incluster interactions
                 ii1=imp2center(i); jj1=imp2center(j)
                 var%ss=dot_product(R_SH_IMP(I,:) &
                      -R_SH_IMP(J,:),R_SH_IMP(I,:)-R_SH_IMP(J,:))
                 var%cs=dot_product(R_NUC_IMP(I,:) &
                      -R_SH_IMP(J,:),R_NUC_IMP(I,:)-R_SH_IMP(J,:))
                 var%sc=dot_product(R_SH_IMP(I,:) &
                      -R_NUC_IMP(J,:),R_SH_IMP(I,:)-R_NUC_IMP(J,:))
                 var%cc=dot_product(R_NUC_IMP(I,:) &
                      -R_NUC_IMP(J,:),R_NUC_IMP(I,:)-R_NUC_IMP(J,:))
                 var%rr=dot_product(R_IMP(I,:)-R_IMP(J,:), &
                      R_IMP(I,:)-R_IMP(J,:))

                 ! **short range interaction : energies and gradients

                 l=>var !rr calculated
                 CALL calc_shortrange(TYPE_IMPURITY(I),TYPE_IMPURITY(J), ii1,jj1,& ! def contr 2
                      clus(i)%short,clus(i)%gshort)
                 impurity_shell_grad(I,:)=impurity_shell_grad(I,:)+ &
                      clus(i)%gshort*(R_SH_IMP(I,:)-R_SH_IMP(J,:))
                 impurity_shell_grad(J,:)=impurity_shell_grad(J,:)- &
                      clus(i)%gshort*(R_SH_IMP(I,:)-R_SH_IMP(J,:))
                 clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SH_IMPURITY(J)/SQRT(var%ss)
                 clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SH_IMPURITY(J)/SQRT(var%cs)
                 clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUC_IMPURITY(J)/SQRT(var%sc)
                 clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUC_IMPURITY(J)/SQRT(var%cc)

                 impurity_shell_grad(I,:)=impurity_shell_grad(I,:)- &
                      clus(i)%fac%ss/var%ss*(R_SH_IMP(I,:)-R_SH_IMP(J,:))- &
                      clus(i)%fac%sc/var%sc*(R_SH_IMP(I,:)-R_NUC_IMP(J,:))
                 impurity_shell_grad(J,:)=impurity_shell_grad(J,:)+ &
                      clus(i)%fac%ss/var%ss*(R_SH_IMP(I,:)-R_SH_IMP(J,:))+ &
                      clus(i)%fac%cs/var%cs*(R_NUC_IMP(I,:)-R_SH_IMP(J,:))
                 impurity_core_grad(I,:)=impurity_core_grad(I,:)- & 
                      clus(i)%fac%cs/var%cs*(R_NUC_IMP(I,:)-R_SH_IMP(J,:))- &
                      clus(i)%fac%cc/var%cc*(R_NUC_IMP(I,:)-R_NUC_IMP(J,:))
                 impurity_core_grad(J,:)=impurity_core_grad(J,:)+ &
                      clus(i)%fac%sc/var%sc*(R_SH_IMP(I,:)-R_NUC_IMP(J,:))+ &
                      clus(i)%fac%cc/var%cc*(R_NUC_IMP(I,:)-R_NUC_IMP(J,:))
                 ! calc  energy contrib (1)
                 clus(i)%coul=clus(i)%fac%ss+clus(i)%fac%sc+clus(i)%fac%cs+clus(i)%fac%cc

              endif ref_data1

           clus(i)%tot=clus(i)%short+clus(i)%coul
           eshort_epecluster=eshort_epecluster+clus(i)%short
           ecoul_epecluster=ecoul_epecluster+clus(i)%coul

        end DO

! **done impurity - impurity iteraction 

! **impurity - lattice interaction 
        clus(i)%coul=zero
        clus(i)%short=zero

! **coulomb energy & gradients
! **a) calculate contribution to the coulomb energy 
! **   due to the evald sumation in the coordinate space
! **b) calculate contributions to the gradients of impurities (imp_shell_grad,DPC)
! **   due to the Evald sumation in the coordinate space

        ec_ewadir_cluster=zero
        imp2vac: DO J=1,N_VACANCIES
           if(use_ref_data) then ! only regular data 
              reg%sr=dot_product(reg_reference(I)%rs-epe(J)%r, &
                   reg_reference(I)%rs-epe(J)%r)
              reg%cr=dot_product(reg_reference(I)%rc-epe(J)%r, &
                   reg_reference(I)%rc-epe(J)%r)
           else
              reg%sr=dot_product(R_SH_IMP(I,:)-epe(J)%r, &
                   R_SH_IMP(I,:)-epe(J)%r )
              reg%cr=dot_product(R_NUC_IMP(I,:)-epe(J)%r, &
                   R_NUC_IMP(I,:)-epe(J)%r)
           endif

           serfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(reg%sr)
           serfunc_arg(2)=serfunc_arg(1)**2
           serfunc_arg(3)=serfunc_arg(2)**2
           cerfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(reg%cr)
           cerfunc_arg(2)=cerfunc_arg(1)**2
           cerfunc_arg(3)=cerfunc_arg(2)**2
           clus(i)%fac%sewa=ERFO*(one-third*serfunc_arg(2)+tenth*serfunc_arg(3))
           clus(i)%fac%cewa=ERFO*(one-third*cerfunc_arg(2)+tenth*cerfunc_arg(3))
           IF(serfunc_arg(2).GT.erfarg_th) clus(i)%fac%sewa=-ERF(serfunc_arg(1))/SQRT(reg%sr)
           if(cerfunc_arg(2).gt.erfarg_th) clus(i)%fac%cewa=-erf(cerfunc_arg(1))/sqrt(reg%cr)

           !calc ewald dir summation energy contrib (2)          
           ec_ewadir_cluster=ec_ewadir_cluster+ Q_ion(epe(J)%k)*( &
                Q_SH_IMPURITY(I)*clus(i)%fac%sewa+ &
                Q_NUC_IMPURITY(I)*clus(i)%fac%cewa)

           ! treat gradients
           clus(i)%fac%gsewa=gewa_const* &
                (1.0_r8_kind-0.6_r8_kind*serfunc_arg(2)+0.214286_r8_kind*serfunc_arg(3))
           clus(i)%fac%gcewa=gewa_const* &
                (1.0_r8_kind-0.6_r8_kind*cerfunc_arg(2)+0.214286_r8_kind*cerfunc_arg(3))
           IF(serfunc_arg(2).GT.gerfarg_th)clus(i)%fac%gsewa= &
                (ERF(serfunc_arg(1))/SQRT(reg%sr)+ERFO*EXP(-serfunc_arg(2)))/reg%sr
           IF(cerfunc_arg(2).GT.gerfarg_th)clus(i)%fac%gcewa= &
                (ERF(cerfunc_arg(1))/SQRT(reg%cr)+ERFO*EXP(-cerfunc_arg(2)))/reg%cr

           impurity_shell_grad(I,:)=impurity_shell_grad(I,:)+(R_SH_IMP(I,:)- &
                epe(j)%r)*Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%gsewa
           impurity_core_grad(I,:)=impurity_core_grad(I,:)+(R_NUC_IMP(I,:)- &
                epe(j)%r)*Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%gcewa
        enddo imp2vac

        clus(i)%coul=clus(i)%coul+ec_ewadir_cluster

! **Evald method reciprocial space contribution to energy (clus(i)%coul)
! **and gradients of impurities (imp_core_grad,imp_shell_grad) 
! **acctually  clus(i)%coul contribution for use_ref_data.eq.true have to be calculated
! **quantumchemically
        gc_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),R_NUC_IMP(i,:))*PI2
           gc_ewaback=gc_ewaback+(RSIN(ind)*COS(angle_fac)-RCOS(ind)*SIN(angle_fac))*GSTR(ind,:)
        enddo
        gc_ewaback=gc_ewaback*PI2
        gs_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),R_SH_IMP(i,:))*PI2
           gs_ewaback=gs_ewaback+(RSIN(ind)*COS(angle_fac)-RCOS(ind)*SIN(angle_fac))*GSTR(ind,:)
        enddo
        gs_ewaback=gs_ewaback*PI2
        impurity_core_grad(I,:)=impurity_core_grad(I,:)+PK_IMPURITY(I)*(R_NUC_IMP(I,:)- &
             R_SH_IMP(I,:))+Q_NUC_IMPURITY(I)*gc_ewaback
        impurity_shell_grad(I,:)=impurity_shell_grad(I,:)-PK_IMPURITY(I)*(R_NUC_IMP(I,:)- &
             R_SH_IMP(I,:))+Q_SH_IMPURITY(I)*gs_ewaback

   !** treat energy contributions
        if(use_ref_data) then
           vs_ewaback=zero
           DO ind=1,n_bs_points
              angle_fac=dot_product(GSTR(ind,:),reg_reference(i)%rs)*PI2
              vs_ewaback=vs_ewaback+ &
                   RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
           enddo
           clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*vs_ewaback
           vc_ewaback=zero
           DO ind=1,n_bs_points
              angle_fac=dot_product(GSTR(ind,:),reg_reference(i)%rc)*PI2
              vc_ewaback=vc_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
           enddo!I=1,n_bs_points
           clus(i)%coul=clus(i)%coul+Q_NUC_IMPURITY(I)*vc_ewaback
           clus(i)%coul=clus(i)%coul + PK_IMPURITY(I)* &
                dot_product(reg_reference(i)%rc-reg_reference(i)%rs, &
                reg_reference(i)%rc-reg_reference(i)%rs)/2.0_r8_kind

        else ! regular run

           ! ewald backspace energy contrib (3)
           vs_ewaback=zero
           DO ind=1,n_bs_points
              angle_fac=dot_product(GSTR(ind,:),R_SH_IMP(i,:))*PI2
              vs_ewaback=vs_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
           enddo!I=1,n_bs_points
           clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*vs_ewaback

           vc_ewaback=zero
           DO ind=1,n_bs_points
              angle_fac=dot_product(GSTR(ind,:),R_nuc_IMP(i,:))*PI2
              vc_ewaback=vc_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
           enddo!I=1,n_bs_points
           clus(i)%coul=clus(i)%coul+Q_NUC_IMPURITY(I)*vc_ewaback


           ! shell model intra ionic energy contrib (4)
           clus(i)%coul=clus(i)%coul + pk_impurity(i)* & !!!!!
                dot_product(r_nuc_imp(i,:)-r_sh_imp(i,:), &
                r_nuc_imp(i,:)-r_sh_imp(i,:))/2.0_r8_kind
        endif
   ! **done Ewald method reciprocial space contributions

! **a)calculate contributions to energy (clus(i)%coul) due to interaction
! **  of impurities with shell model ions of region A and 
! **  corresponding contributions of Ewald summation in the 
! **  coordinate space 
! **b)calculate contributions to gradients on atoms of region A 
! **  (D1CUL,D1SUL) due to interactions with impurities (lcgto-cluster)
! **c)calculate contributions to gradients on impurities (imp_core_grad,imp_shell_grad)
! **  due to interaction with ions of region 'A' and due to 
! **  Ewald summation in the coordinate space 


        ! for use_ref_data mode add regular and substract reference contribs
        imp2regI: DO J=n_vacancies+1,reg_I_n_ions

           !calc dot products
           var%sr=dot_product(R_SH_IMP(I,:)-epe(J)%r,R_SH_IMP(I,:)-epe(J)%r)
           var%cr=dot_product(R_NUC_IMP(I,:)-epe(J)%r,R_NUC_IMP(I,:)-epe(J)%r)
           var%ss=dot_product(R_SH_IMP(I,:)-epe(J)%s,R_SH_IMP(I,:)-epe(J)%s)
           var%cc=dot_product(R_NUC_IMP(I,:)-epe(J)%c,R_NUC_IMP(I,:)-epe(J)%c)
           var%sc=dot_product(R_SH_IMP(I,:)-epe(J)%c,R_SH_IMP(I,:)-epe(J)%c)
           var%cs=dot_product(R_NUC_IMP(I,:)-epe(J)%s,R_NUC_IMP(I,:)-epe(J)%s)
           var%rr=dot_product(R_IMP(I,:)-epe(J)%r,R_IMP(I,:)-epe(J)%r)

           ! error_function contribs
           serfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(var%sr)
           serfunc_arg(2)=serfunc_arg(1)**2
           cerfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(var%cr)
           cerfunc_arg(2)=cerfunc_arg(1)**2
           clus(i)%fac%sewa=ERF(serfunc_arg(1))/SQRT(var%sr)
           clus(i)%fac%cewa=ERF(cerfunc_arg(1))/SQRT(var%cr)

           ! energy prefactors  
           clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%cc)
           clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%cs)
           clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%sc)
           clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%ss)

           if(use_ref_data) then
              reg%ss=dot_product(reg_reference(i)%rs-epe(j)%s, &
                   reg_reference(i)%rs- epe(j)%s)
              reg%sc=sqrt(dot_product(reg_reference(i)%rs-epe(j)%c, &
                   reg_reference(i)%rs-epe(j)%c ))
              reg%sr=sqrt(dot_product(reg_reference(i)%rs-epe(J)%r,  &
                   reg_reference(i)%rs-epe(J)%r  ))
              reg%cr=sqrt(dot_product(reg_reference(i)%rc-epe(J)%r,  &
                   reg_reference(i)%rc-epe(J)%r  ))
              reg%cs=sqrt(dot_product(reg_reference(i)%rc-epe(j)%s, &
                   reg_reference(i)%rc-epe(j)%s ))
              reg%cc=sqrt(dot_product(reg_reference(i)%rc-epe(j)%c, &
                   reg_reference(i)%rc-epe(j)%c ))
              reg%rr=var%rr
              clus(i)%fac%sewa= ERF(ERROR_FUNCTION_PARAMETER*reg%sr)/reg%sr
              clus(i)%fac%cewa= ERF(ERROR_FUNCTION_PARAMETER*reg%cr)/reg%cr
              !*** coulon contrib regular data
              clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*( &
                   Q_SHELL(epe(j)%k)/sqrt(reg%ss)+ Q_NUCLEAR(epe(j)%k)/reg%sc &
                   -Q_ion(epe(J)%k)*clus(i)%fac%sewa)&
                   +Q_NUC_IMPURITY(I)*( &
                   Q_SHELL(epe(j)%k)/reg%cs+ Q_NUCLEAR(epe(j)%k)/reg%cc- &
                   Q_ion(epe(J)%k)*clus(i)%fac%cewa)

              if(i.eq.n_impurities.and.use_epe_pgdata) diffpg_ec_ecref=diffpg_ec_ecref &
                   -(reg_I_pg(j)%vs*Q_SHELL(epe(j)%k)+ &
                   reg_I_pg(j)%vc*Q_NUCLEAR(epe(j)%k))*7.17115929581458_r8_kind*scale_factor
              ! add (PG_var-PG_ref) contrib

           else ! regular run not use_ref_data mode

              ! qq/r + ewald err function contrib
              clus(i)%coul=clus(i)%coul+ &
                   clus(i)%fac%ss+clus(i)%fac%sc &
                   -Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &   
                   + clus(i)%fac%cc+clus(i)%fac%cs &
                   -Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa    
           endif        !       use_ref_data/else

           ! prefactors for gradients
           clus(i)%fac%cc=clus(i)%fac%cc/(var%cc)
           clus(i)%fac%cs=clus(i)%fac%cs/(var%cs)
           clus(i)%fac%sc=clus(i)%fac%sc/(var%sc)
           clus(i)%fac%ss=clus(i)%fac%ss/(var%ss)
           clus(i)%fac%gsewa=Q_SH_IMPURITY(I)*Q_ion(epe(J)%k) &
                *(clus(i)%fac%sewa+ERFO*EXP(-serfunc_arg(2)))/var%sr
           clus(i)%fac%gcewa=Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k) &
                *(clus(i)%fac%cewa+ERFO*EXP(-cerfunc_arg(2)))/var%cr

           if(use_ref_data) then
              point_x_core_grad(j,:)=point_x_core_grad(j,:)+& 
                   (reg_reference(i)%rc-R_NUC_ION(J,:))*Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/reg%cc**3 &
                   +(reg_reference(i)%rs-R_NUC_ION(J,:))*Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/reg%sc**3
              point_x_shell_grad(j,:)=point_x_shell_grad(j,:)+& !***  coulon contribution for epe shells
                   (reg_reference(i)%rc-R_SH_ION(J,:))*Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/reg%cs**3 + &
                   (reg_reference(i)%rs-R_SH_ION(J,:))*Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/ &
                   (reg%ss*sqrt(reg%ss))
              if(i.eq.n_impurities.and.use_epe_pgdata) then !*** delta contribution from PG
                 point_x_core_grad(J,:)=point_x_core_grad(J,:)+ &
                      reg_I_pg(j)%gc(:)*13.5515309044_r8_kind*Q_NUCLEAR(epe(j)%k)*scale_factor
                 point_x_shell_grad(J,:)=point_x_shell_grad(J,:)+ &
                      reg_I_pg(j)%gs(:)*13.5515309044_r8_kind*Q_SHELL(epe(j)%k)*scale_factor
              endif! i.eq.n_impurities
           else
              point_x_core_grad(J,:)=point_x_core_grad(J,:) &
                   +clus(i)%fac%cc*(R_NUC_IMP(I,:)-R_NUC_ION(J,:)) &
                   +clus(i)%fac%sc*(R_SH_IMP(I,:)-R_NUC_ION(J,:))

              point_x_shell_grad(J,:)=point_x_shell_grad(J,:) &
                   +clus(i)%fac%cs*(R_NUC_IMP(I,:)-R_SH_ION(J,:)) &
                   +clus(i)%fac%ss*(R_SH_IMP(I,:)-R_SH_ION(J,:))
           endif

           impurity_core_grad(I,:)=impurity_core_grad(I,:)+clus(i)%fac%gcewa* &
                (R_NUC_IMP(I,:)-epe(J)%r) &
                -clus(i)%fac%cc*(R_NUC_IMP(I,:)-R_NUC_ION(J,:)) &
                -clus(i)%fac%cs*(R_NUC_IMP(I,:)-R_SH_ION(J,:))
           impurity_shell_grad(I,:)=impurity_shell_grad(I,:)+clus(i)%fac%gsewa* &
                (R_SH_IMP(I,:)-epe(J)%r) &
                -clus(i)%fac%ss*(R_SH_IMP(I,:)-R_SH_ION(J,:)) &
                -clus(i)%fac%sc*(R_SH_IMP(I,:)-R_NUC_ION(J,:))

           ! **calculate short-range interaction contributions to 
           ! **energy (clus(i)%short) and gradient on ions of region A (short_range_grad)
           ! **and on impurities impurity_shell_grad()
           ! **short-range interaction impurity-lattice
           
           l=>var !rr calculated

           if(use_ref_data.and.explicit_coupling) then
              ii1=imp2center(i)
              CALL calc_shortrange(explicit_coupling_type(i),epe(j)%k,ii1,j, &   !def contr 3
                   clus(i)%eshort,clus(i)%gshort,harmonic_bond_fix=.true.)
                   ! bonds can be fixed only for explicit_coupling 
                   ! and it is done for var configuration only


           else
              ii1=imp2center(i)
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j, &
                   clus(i)%eshort,clus(i)%gshort)
           end if                       

           clus(i)%short=clus(i)%short+clus(i)%eshort ! No 1 (reg1 var)
           eshort_reg_I=eshort_reg_I+clus(i)%eshort 
 
           short_range_grad1(J,:)=short_range_grad1(J,:)- &
                clus(i)%gshort*(R_SH_IMP(I,:)-epe(j)%s)
           impurity_shell_grad(I,:)=impurity_shell_grad(I,:)+ &
                clus(i)%gshort*(R_SH_IMP(I,:)-epe(j)%s)

!            if(i.eq.70.and.sqrt(sum( (R_SH_IMP(I,:)-epe(j)%s)**2)).lt.2.00) then
!             print*, R_SH_IMP(I,:), i,j,'R_SH_IMP', &
!                     sqrt(sum( (R_SH_IMP(I,:)-epe(j)%s)**2))
!             print*, clus(i)%gshort*(R_SH_IMP(I,:)-epe(j)%s)
!             print*
!            endif

           refdat2: if(use_ref_data) then
           ! two reference configurations
              ii1=imp2center(i)
              l=>reg !rr calculated
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j, &   !def contr 4
                   clus(i)%eshort,clus(i)%gshort)
                   ! for regular configuration parameters are not redefined
                   ! with expicit coupling


!              short_forces=short_forces-clus(i)%gshort*(reg_reference(i)%rs-epe(j)%s)

!            if(i.eq.70.and.sqrt(sum( (reg_reference(i)%rs-epe(j)%s)**2)).lt.2.00) then
!             print*, reg_reference(i)%rs, i,j,'reg_reference', &
!                     sqrt(sum( (reg_reference(i)%rs-epe(j)%s)**2))
!             print*, clus(i)%gshort*(reg_reference(i)%rs-epe(j)%s)
!             print*
!            endif
                   
              short_range_grad1(J,:)=short_range_grad1(J,:)- &
                   clus(i)%gshort*(reg_reference(i)%rs-epe(j)%s)
              clus(i)%short=clus(i)%short+clus(i)%eshort ! (reg1 reg)

              ref%ss=dot_product(epe_reference(i)%rs-epe(j)%s, &
                   epe_reference(i)%rs-epe(j)%s )
              ref%rr=var%rr

              l=>ref !rr calculated
              if(explicit_coupling) then
                 ii1=imp2center(i)
                 CALL calc_shortrange(explicit_coupling_type(i),epe(j)%k,ii1,j, & !def contr 5
                      clus(i)%eshort,clus(i)%gshort)
              else
                 CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,i,j, &    ! def contr 6
                      clus(i)%eshort,clus(i)%gshort)
              end if
              clus(i)%short=clus(i)%short-clus(i)%eshort ! (reg1 ref)
              short_range_grad1(J,:)=short_range_grad1(J,:)+ &
                   clus(i)%gshort*(epe_reference(i)%rs-epe(j)%s)
!            if(i.eq.70.and.sqrt(sum( (epe_reference(i)%rs-epe(j)%s)**2)).lt.2.00) then
!             print*, epe_reference(i)%rs, i,j,'epe_reference', &
!                     sqrt(sum( (epe_reference(i)%rs-epe(j)%s)**2))
!             print*, clus(i)%gshort*(epe_reference(i)%rs-epe(j)%s)
!             print*
!            endif
         endif refdat2
        enddo imp2regI

        ! ** 2a-region start
        ec_reg2a=zero
        imp2reg2a: DO J=reg_I_n_ions+1,reg_2a_n_ions

           ! calc dot products
           var%rr=dot_product(r_imp(i,:)-epe(j)%r,r_imp(i,:)-epe(j)%r)
           var%sr=sqrt(dot_product(R_SH_IMP(I,:)-epe(j)%r,R_SH_IMP(I,:)-epe(j)%r))
           var%cr=sqrt(dot_product(R_NUC_IMP(I,:)-epe(j)%r,R_NUC_IMP(I,:)-epe(j)%r))
           var%ss=dot_product(R_SH_IMP(I,:)-epe(j)%s,R_SH_IMP(I,:)-epe(j)%s)
           var%cc=dot_product(R_NUC_IMP(I,:)-epe(j)%c,R_NUC_IMP(I,:)-epe(j)%c)
           var%sc=dot_product(R_SH_IMP(I,:)-epe(j)%c,R_SH_IMP(I,:)-epe(j)%c)
           var%cs=dot_product(R_NUC_IMP(I,:)-epe(j)%s,R_NUC_IMP(I,:)-epe(j)%s)
           var%ss_sr=dot_product(R_SH_IMP(I,:)-epe(j)%s,epe(j)%s-epe(j)%r)
           var%sc_cr=dot_product(R_SH_IMP(I,:)-epe(j)%c,epe(j)%c-epe(j)%r)
           var%cc_cr=dot_product(R_NUC_IMP(I,:)-epe(j)%c,epe(j)%c-epe(j)%r)
           var%cs_sr=dot_product(R_NUC_IMP(I,:)-epe(j)%s,epe(j)%s-epe(j)%r)

           ! var shortrange contribs
           l=>var !rr calculated
           if(use_ref_data.and.explicit_coupling) then
              ii1=imp2center(i)
              CALL calc_shortrange(explicit_coupling_type(I),epe(j)%k,ii1,j, & !def contr 7
                   clus(i)%eshort,clus(i)%gshort)
           else
              ii1=imp2center(i)
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j, & ! def contr 8
                   clus(i)%eshort,clus(i)%gshort)
           end if

           clus(i)%short=clus(i)%short+clus(i)%eshort  ! No 2 (reg2 var)

           eshort_reg2a=eshort_reg2a+clus(i)%eshort   

           eshort_lat_relcontrib=0.5_r8_kind*clus(i)%gshort*var%ss_sr !(var)

           if(.not.fixed_reg2a_relcontribs) then
              !switched off  as  no dipole contributions is calculated in PG
              eshort_lat_relaxation=eshort_lat_relaxation+eshort_lat_relcontrib
              eshort_lat_relcorr=eshort_lat_relcorr+eshort_lat_relcontrib
              ! contrib No 2 cluster, var configuration
              !else contribs only due cluster simulator
           endif

         ! var longrange prefactors
           clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%cc) 
           clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%cs)  
           clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%sc)
           clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%ss)
           clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*var%sr)/var%sr
           clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*var%cr)/var%cr 

        !var  qq/r + ewald error function contrib
        ! the same contribution is calculated with PG
        ! when cluster is embedded in epe PC
           ec_current= &
                clus(i)%fac%ss+clus(i)%fac%sc + &
                clus(i)%fac%cc+clus(i)%fac%cs  &
                -Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &
                -Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa

           ec_reg2a=ec_reg2a+ec_current
           ec_cluster_reg2a=ec_cluster_reg2a+ec_current 
           diff_reg2a_ec_ecref=diff_reg2a_ec_ecref+ec_current

           clus(i)%fac%cc=clus(i)%fac%cc/var%cc 
           clus(i)%fac%cs=clus(i)%fac%cs/var%cs
           clus(i)%fac%sc=clus(i)%fac%sc/var%sc
           clus(i)%fac%ss=clus(i)%fac%ss/var%ss

        !(1) var rel energy contrib --------------------------------------------------
           if(.not.fixed_reg2a_relcontribs) then
           ecoul_lat_relcontrib= &
                0.5*( clus(i)%fac%ss*var%ss_sr+ clus(i)%fac%sc*var%sc_cr + &
                clus(i)%fac%cc*var%cc_cr + clus(i)%fac%cs*var%cs_sr )

           ecoul_lat_relcorr=ecoul_lat_relcorr-ecoul_lat_relcontrib

           ecoul_lat_relaxation=ecoul_lat_relaxation- ecoul_lat_relcontrib
        !var contribution due to ML dicpacements 
        !these contributions are equal to zero for frozen reg2A-------------
           !else contribution only due cluster simulator, see below
        end if
        
        !** now treat var configuration for forces
           IF(ERROR_FUNCTION_PARAMETER*var%sr.LT.tail_exp_dist) then
              clus(i)%fac%gsewa=Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*(clus(i)%fac%sewa+ &
                   ERFO*EXP(-ET2*var%sr**2))/var%sr**2
              clus(i)%fac%gcewa=Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*(clus(i)%fac%cewa+ &
                   ERFO*EXP(-ET2*var%cr**2))/var%cr**2
           else
              clus(i)%fac%gsewa=Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)/(var%sr**3)
              clus(i)%fac%gcewa=Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)/(var%cr**3)
           endif

           impurity_core_grad(I,:)=impurity_core_grad(I,:)+ &
                clus(i)%fac%gcewa*(R_NUC_IMP(I,:)-epe(J)%r) &
                -clus(i)%fac%cc*(R_NUC_IMP(I,:)-R_NUC_ION(J,:)) &
                -clus(i)%fac%cs*(R_SH_IMP(I,:)-R_NUC_ION(J,:)) 
           if(.not.fixed_reg2a_relcontribs) then
              impurity_core_grad(I,:)=impurity_core_grad(I,:) &
                   -0.5_r8_kind*(clus(i)%fac%cs*((epe(J)%r-R_SH_ION(J,:)) &
                   -3.0_r8_kind*(R_NUC_IMP(I,:)-R_SH_ION(J,:))*var%cs_sr/var%cs) &
                        +clus(i)%fac%cc*((epe(J)%r-R_NUC_ION(J,:)) &
                   -3.0_r8_kind*(R_NUC_IMP(I,:)-R_NUC_ION(J,:))*var%cc_cr/var%cc))
           end if
        
           impurity_shell_grad(I,:)=impurity_shell_grad(I,:)+ &
                clus(i)%fac%gsewa*(R_SH_IMP(I,:)-epe(J)%r) &
                +clus(i)%gshort*(R_SH_IMP(I,:)-R_SH_ION(J,:)) &
                -clus(i)%fac%ss*(R_SH_IMP(I,:)-R_SH_ION(J,:)) &
                -clus(i)%fac%sc*(R_SH_IMP(I,:)-R_NUC_ION(J,:)) 
           if(.not.fixed_reg2a_relcontribs) then
              impurity_shell_grad(I,:)=impurity_shell_grad(I,:) &  
                   -0.5_r8_kind*(clus(i)%fac%ss*((epe(J)%r-R_SH_ION(J,:)) &
                   -3.0_r8_kind*(R_SH_IMP(I,:)-R_SH_ION(J,:))*var%ss_sr/var%ss) &
                        +clus(i)%fac%sc*((epe(J)%r-R_NUC_ION(J,:)) &
                   -3.0_r8_kind*(R_SH_IMP(I,:)-R_NUC_ION(J,:))*var%sc_cr/var%sc))
           end if
             

           if(use_ref_data) then ! modify code for energy only
              reg%ss=dot_product(reg_reference(i)%rs-epe(J)%s, &
                   reg_reference(i)%rs-epe(J)%s)
              reg%sr=sqrt(dot_product(reg_reference(i)%rs-epe(j)%r, &
                   reg_reference(i)%rs-epe(j)%r))
              reg%cr=sqrt(dot_product(reg_reference(i)%rc-epe(j)%r, &
                   reg_reference(i)%rc-epe(j)%r))
              reg%cc=dot_product(reg_reference(i)%rc-epe(J)%c, &
                   reg_reference(i)%rc-epe(J)%c)
              reg%sc=dot_product(reg_reference(i)%rs-epe(J)%c, &
                   reg_reference(i)%rs-epe(J)%c)
              reg%cs=dot_product(reg_reference(i)%rc-epe(J)%s, &
                   reg_reference(i)%rc-epe(J)%s)
              reg%rr=dot_product(r_imp(i,:)-epe(j)%r,r_imp(i,:)-epe(j)%r)

              clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(reg%cc)
              clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(reg%cs)
              clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(reg%sc)
              clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(reg%ss)
              clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*reg%sr)/reg%sr
              clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*reg%cr)/reg%cr

           ! var longrange reg contribs for use_ref_data mode 
              ec_reg2a=ec_reg2a+clus(i)%fac%ss+clus(i)%fac%sc- &
                   Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &
                   +clus(i)%fac%cc+clus(i)%fac%cs- &
                   Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa
           
           ! reg shortrange contribs for use_ref_data mode
              reg%ss_sr=dot_product(reg_reference(i)%rs-epe(j)%s,epe(j)%s-epe(j)%r)
              l=>reg !rr calculated
              ii1=imp2center(i)
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(J)%k,ii1,j,clus(i)%eshort,clus(i)%gshort) !def contr 8
              clus(i)%short=clus(i)%short+clus(i)%eshort ! (reg2 reg)

              if(.not.fixed_reg2a_relcontribs) then
              eshort_lat_relaxation=eshort_lat_relaxation+ &
                   0.5_r8_kind*clus(i)%gshort*reg%ss_sr
              ! contrib No 3 cluster, regular configuration

              
              !*** (2) reg rel contrib for use_ref_data mode
                 ecoul_lat_relaxation=ecoul_lat_relaxation- &
                      0.5_r8_kind*( clus(i)%fac%ss*reg%ss_sr/reg%ss+ &
                      clus(i)%fac%sc*dot_product(reg_reference(i)%rs-epe(j)%c,epe(j)%c- &
                      epe(j)%r)/reg%sc+ &
                      clus(i)%fac%cc*dot_product(reg_reference(i)%rc-epe(j)%c,epe(j)%c- &
                      epe(j)%r)/reg%cc+ &
                      clus(i)%fac%cs*dot_product(reg_reference(i)%rc-epe(j)%s,epe(j)%s- &
                      epe(j)%r)/reg%cs ) 
              end if
           
        !** now treat reference
              ref%ss=dot_product(epe_reference(i)%rs-epe(j)%s,epe_reference(i)%rs-epe(j)%s)
              ref%sr=sqrt(dot_product(epe_reference(i)%rs-epe(j)%r,epe_reference(i)%rs-epe(j)%r))
              ref%cr=sqrt(dot_product(epe_reference(i)%rc-epe(j)%r,epe_reference(i)%rc-epe(j)%r))
              ref%cc=dot_product(epe_reference(i)%rc-epe(j)%c,epe_reference(i)%rc-epe(j)%c)
              ref%sc=dot_product(epe_reference(i)%rs-epe(j)%c,epe_reference(i)%rs-epe(j)%c)
              ref%cs=dot_product(epe_reference(i)%rc-epe(j)%s,epe_reference(i)%rc-epe(j)%s)
              ref%rr=var%rr

              clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(ref%cc) 
              clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(ref%cs)   
              clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(ref%sc)
              clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(ref%ss)
              clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*ref%sr)/ref%sr
              clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*ref%cr)/ref%cr 

           ! ref lonrange contrib for use_ref_data mode
              ec_reg2a_ref=     clus(i)%fac%ss+clus(i)%fac%sc- &
                   Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa  &
                   +clus(i)%fac%cs+clus(i)%fac%cc- &
                   Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa
              ec_reg2a=ec_reg2a-ec_reg2a_ref

              diff_reg2a_ec_ecref=diff_reg2a_ec_ecref-ec_reg2a_ref 

              ! ref shortrange contribs for use_ref_data mode ----------------------------
!new              ref%ss_sr=dot_product(epe_reference(i)%rs-epe(j)%s,epe(j)%s-epe(j)%r)
              ref%ss_sr=dot_product(epe_reference(i)%rs-R_SH_ION(J,:),R_SH_ION(J,:)-epe(j)%r)
              l=>ref !rr calculated

              if(explicit_coupling) then
                 ii1=imp2center(i)
                 CALL calc_shortrange(explicit_coupling_type(I),epe(J)%k,ii1,j, & !def cont 9
                      clus(i)%eshort,clus(i)%gshort)
              else
                 CALL calc_shortrange(TYPE_IMPURITY(I),epe(J)%k,i,j, &  ! def cont 10
                      clus(i)%eshort,clus(i)%gshort)
              endif
              clus(i)%short=clus(i)%short-clus(i)%eshort  ! (reg2 ref)

              if(.not.fixed_reg2a_relcontribs) then
                 eshort_lat_relaxation=eshort_lat_relaxation- &
                      0.5_r8_kind*clus(i)%gshort*ref%ss_sr 
                 ! contrib No 4, reference configuration

                 ! (3) rel ref long-range contrib
                 ecoul_lat_relaxation=ecoul_lat_relaxation+0.5_r8_kind*( &
                      clus(i)%fac%ss*ref%ss_sr/ref%ss  & 
                      +clus(i)%fac%sc*dot_product(epe_reference(i)%rs-epe(j)%c, &
                      epe(j)%c-epe(j)%r)/ref%sc  &
                      +clus(i)%fac%cc*dot_product(epe_reference(i)%rc-epe(j)%c, &
                      epe(j)%c-epe(j)%r)/ref%cc &
                      +clus(i)%fac%cs*dot_product(epe_reference(i)%rc-epe(j)%s, &
                      epe(j)%s-epe(j)%r)/ref%cs ) 
                 ! contributions due to ML dicplacenents for reference configuration
              end if
           
           endif
        enddo imp2reg2a

        clus(i)%tot=clus(i)%coul+ec_reg2a+clus(i)%short
        eshort_epecluster_epe=eshort_epecluster_epe+clus(i)%short
        ec_cluster_epe=ec_cluster_epe+ec_reg2a+clus(i)%coul
        reg_2a_ec(i)=ec_reg2a

     enddo imp2rest
     endif i_1isimp
     endif imp_present

     if(qm_interfaced_mode.and. comm_i_am_master()) then ! treat impurity contrib
        if(n_types_central_atoms_3body > 0) then
         if(use_ref_data) then
           e_imp_3b=0.0_r8_kind
           call imp_3_body_1(e_imp_3bT,reg=.true.) ! in calc_def_contributions
           DPRINT 'e_imp_3bT', e_imp_3bT
           e_imp_3b=e_imp_3b+e_imp_3bT
           call imp_3_body_1(e_imp_3bT,ref=.true.) ! in calc_def_contributions
           DPRINT 'e_imp_3bT', e_imp_3bT
           e_imp_3b=e_imp_3b+e_imp_3bT     !calculated with -1 prefactor inside imp_3_body_1
           call imp_3_body_1(e_imp_3bT)            ! in calc_def_contributions
           DPRINT 'e_imp_3bT', e_imp_3bT
           e_imp_3b=e_imp_3b+e_imp_3bT
           eshort_reg_I=eshort_reg_I+e_imp_3bT
           cross_boundary_3b=e_imp_3bT
           ! this term eshort_reg_I will be substructed from energy in gxfile
           ! to be calculated  again in optimizer run
           print*, '3b term to be recalculated to the same value in optimizer', e_imp_3bT
         else
           call imp_3_body_1(e_imp_3b,reg=.true.) ! in calc_def_contributions
         endif
        else
           e_imp_3b=zero
        endif
     else
        if(comm_i_am_master()) then
           if(n_types_central_atoms_3body_im > 0) then
              call imp_3_body(e_imp_3b)
           else
              e_imp_3b=0.0_r8_kind
           endif
        endif
     endif
     DPRINT 'e_imp_3b in calc_def_contributions', e_imp_3b
!AG[
     if(comm_parallel() .and. (n_impurities > 0 .or. n_vacancies >0)) then
        if( comm_i_am_master() ) then
           call receive_intermediate_imp()
        else
           call send_intermediate_imp()
        end if
     end if

     if(comm_i_am_master() ) then  
        do I=1,n_impurities
           ec_reg2a = reg_2a_ec(i)
           WRITE(output_epe, 104) I, clus(i)%coul+ec_reg2a, clus(i)%short, clus(i)%tot
104        FORMAT(1X,'  impurity ',I2,' - lattice        !',F16.8, ' !',F16.8,' !',F16.8)
        end do
        if(n_impurities > 0) then
           WRITE(output_epe,*) 'sum_over_impurities_short ',sum(clus(:)%short)
           WRITE(output_epe,*) 'sum_over_impurities_coul ',sum(clus(:)%coul+reg_2a_ec(:))
        IF(n_types_central_atoms_3body_im.gt.0.or.qm_interfaced_mode) &
           WRITE(output_epe,*) 'sum_over_impurities_3_body ',e_imp_3b
           WRITE(output_epe,*) 'sum_over_impurities_total ',sum(clus(:)%tot)+e_imp_3b
        endif
        if(n_impurities /= 0) deallocate( reg_2a_ec ,stat=status)
        if(status.ne.0) call error_handler(" deallocate reg_2a_ec failed")
     end if! (only master prints)

!AG]
     if(n_impurities /= 0) then 
        deallocate(clus,stat=status)
        if(status.ne.0) call error_handler("deallocate clus failed")
     end if
     ec_cluster_epe=ec_cluster_epe+diffpg_ec_ecref
     ! **region 2a finished
     ! **impurities done

!AG[
!AG        Only master prints
!AG]
     if( comm_i_am_master() ) then
        defect_energy_short=eshort_vaccluster+eshort_vacepe+ &
             eshort_epecluster+eshort_epecluster_epe+e_imp_3b
        defect_energy_coul=ec_cluster_epe &
             +ecoul_vaccluster_epe+ecoul_vaccluster+ecoul_epecluster+E2BDEF

        if(ml_cluster_simulated.and.n_ml_cluster_simulators.ne.0) then
           allocate(ml_simulator(n_ml_cluster_simulators),stat=status)
           if(status.ne.0) call error_handler("allocate ml_simulator failed")
           do i=1,n_ml_cluster_simulators
              do j=reg_I_n_ions+1,reg_2a_n_ions
                 sim%rs=dot_product(ml_cluster_simulators(i)%r-epe(j)%s, &
                                    ml_cluster_simulators(i)%r-epe(j)%s)
                 sim%rc=dot_product(ml_cluster_simulators(i)%r-epe(j)%c, &
                                    ml_cluster_simulators(i)%r-epe(j)%c)
                 
                 ml_simulator(i)%fac%rs= &
                      ml_cluster_simulators(i)%q*Q_SHELL(epe(j)%k)/sqrt(sim%rs) 
                 ml_simulator(i)%fac%rc= &
                      ml_cluster_simulators(i)%q*Q_NUCLEAR(epe(j)%k)/sqrt(sim%rc) 

                    ! (4) rel  long-range contrib
                 ecoul_lat_relaxation=ecoul_lat_relaxation+ &
                      0.5_r8_kind*( &
                      ml_simulator(i)%fac%rs* &
                      dot_product(ml_cluster_simulators(i)%r-epe(j)%s, &
                                                    epe(j)%r-epe(j)%s)/sim%rs &
                     +ml_simulator(i)%fac%rc* &
                      dot_product(ml_cluster_simulators(i)%r-epe(j)%c, &
                                                    epe(j)%r-epe(j)%c)/sim%rc)
              end do

           end do
           deallocate(ml_simulator,stat=status) 
           if(status.ne.0) call error_handler("deallocate ml_simulator failed")
           
        end if

        IF(N_VACANCIES.ne.0 .or. N_IMPURITIES.ne.0) then
           WRITE(output_epe, 109)E2BDEF,E2BDEF
109        format(1x,' constitutiens of  2b          !',F16.8,' !',16X, &
                ' !',F16.8)
           defect_contrib=defect_energy_short+defect_energy_coul
           WRITE(output_epe &
                &, "(1X,'  defect-defect+defect-lattice !',F16.8, &
                &' !',F16.8,' !',F16.8)") &
                defect_energy_coul,defect_energy_short, &
                defect_contrib
           WRITE(output_epe, 106)
        endif
!!$        WRITE(output_epe,*) 'eshort_lat_relaxation,ecoul_lat_relaxation '
!!$        WRITE(output_epe,*) eshort_lat_relaxation,ecoul_lat_relaxation
        !  e_lat_relaxation=eshort_lat_relaxation+ecoul_lat_relaxation
        ecrr=eshort_lat_relaxation+ecoul_lat_relaxation
        !AG
     end if! master print

!                         Now return results to master 
     if( .not.comm_i_am_master() ) then
        call comm_init_send(comm_master_host, msgtag_epe_send_def)
        call commpack(point_x_core_grad(1,1), 3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_send_defect_contrib:  error [1]")
        call commpack(point_x_shell_grad(1,1), 3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_send_defect_contrib:  error [2]")
        call commpack(short_range_grad1(1,1),   3*reg_I_n_ions, 1, status)
        if(status.ne.0) call error_handler("epe_send_defect_contrib:  error [3]")   
        allocate( impurity_core_grad_send(n_impurities,3) , &
             impurity_shell_grad_send(n_impurities,3),stat=status )
        if(status.ne.0) call error_handler("allocate impurity_core_grad_send failed")
        impurity_core_grad_send(1:n_impurities,:) = zero
        impurity_shell_grad_send(1:n_impurities,:)= zero
        impurity_core_grad_send(1:n_impurities,:) = impurity_core_grad(1:n_impurities,:)
        impurity_shell_grad_send(1:n_impurities,:)= impurity_shell_grad(1:n_impurities,:)
        call commpack(impurity_core_grad_send(1,1), 3*n_impurities, 1, status)
        if(status.ne.0) call error_handler("epe_send_defect_contrib:  error [3]")
        call commpack(impurity_shell_grad_send(1,1), 3*n_impurities, 1, status)
        if(status.ne.0) call error_handler("epe_send_defect_contrib:  error [4]")
        call comm_send()
        deallocate(impurity_core_grad_send, &
             impurity_shell_grad_send,stat=status)   
        if(status.ne.0) call error_handler(" deallocate impurity_core_grad_send failed")
     end if! (slave work)
!AG]
   contains ! of calc_def_contributions

     SUBROUTINE calc_shortrange(K,N,ik,jn,EP,GP,harmonic_bond_fix)
       ! **potential of interaction of ions with a lattice and
       ! **each other.

       use epecom_module, only: host
    
       real(kind=r8_kind) :: EP,GP
       integer(kind=i4_kind) :: M,ik,jn
       integer(kind=i4_kind),intent(in) :: K,N
       logical,optional, intent(in):: harmonic_bond_fix
  ASSERT(K.ne.0)
  ASSERT(K.le.size(host%ro,2))
  ASSERT(N.ne.0)
  ASSERT(N.le.size(host%ro,2))

       l%sq_ss=sqrt(l%ss)
       l%ss3=l%sq_ss*l%ss
       l%ss6=l%ss3**2
       l%ss8=l%ss6*l%ss
       l%ss10=l%ss8*l%ss

       ! **calculation of energy and gradients for remote neighbours

!       if(l%ss .gt. 16.0_r8_kind) then
       if(l%rr .gt. 16.0_r8_kind) then
          M=1
       else
          M=common_atom(ik,jn)+1                 !!!!!!!!!!!!!!AS
          if(host%ro(K,N,M)==0.0_r8_kind) M=1    !!!!!!!!!!!!!!AS
       end if
     ASSERT(M.le.size(host%ro,3))

!       if(l%ss.le.host%sr1(K,N,M)**2) then
       if(l%rr.le.host%sr1(K,N,M)**2) then
          EP =-host%C(K,N,M)/(l%ss6)-host%D(K,N,M)/(l%ss8)
          GP = 6.0_r8_kind*host%C(K,N,M)/(l%ss8)+8.0_r8_kind*host%D(K,N,M)/(l%ss10)

          IF(l%rr.le.host%sr2(K,N,M)**2) then
             EP=EP+host%B(K,N,M)*EXP(-l%sq_ss/host%RO(K,N,M))
             GP=GP-host%B(K,N,M)/(host%RO(K,N,M)*l%sq_ss)*EXP(-l%sq_ss/host%RO(K,N,M))
          endif
       else
          EP = 0.0_r8_kind
          GP= 0.0_r8_kind
       end if
           if(present(harmonic_bond_fix)) then
           ! r in a.u , energy in eV
           if(host%k(K,N,m).ne.0.0.and.abs(sqrt(l%rr)-host%r0(K,N,m)).lt.0.1) then
           DPRINT 'main_epe harmonic_bond_fix',l%sq_ss
           EP=EP+host%k(K,N,m)*(eau_ev/auangs**2)*(l%sq_ss-host%r0(K,N,m))**2
           GP=GP+2.0_r8_kind*(eau_ev/auangs**2)*host%k(K,N,m)*(l%sq_ss-host%r0(K,N,m))/l%sq_ss
!               eau_ev=27.21..... energy converded to eV as alwais in basic code
           endif
           endif


     END SUBROUTINE calc_shortrange

     subroutine vac_to_vac_contibs

       vaca(i)%coul=zero

       do j=1,n_vacancies
          if(i.eq.j) cycle
          var%rr=dot_product(epe(I)%r-epe(J)%r,epe(I)%r-epe(J)%r)
          var%ss=dot_product(epe(I)%s-epe(J)%s,epe(I)%s-epe(J)%s)
          var%cs=dot_product(epe(I)%c-epe(J)%s,epe(I)%c-epe(J)%s)
          var%sc=dot_product(epe(I)%s-epe(J)%c,epe(I)%s-epe(J)%c)
          var%cc=dot_product(epe(I)%c-epe(J)%c,epe(I)%c-epe(J)%c)

          vaca(i)%fac%cc=q_nuclear(epe(i)%k)*q_nuclear(epe(j)%k)/sqrt(var%cc)
          vaca(i)%fac%cs=q_nuclear(epe(i)%k)*q_shell(epe(j)%k)/sqrt(var%cs)
          vaca(i)%fac%sc=q_shell(epe(i)%k)*q_nuclear(epe(j)%k)/sqrt(var%sc)
          vaca(i)%fac%ss=q_shell(epe(i)%k)*q_shell(epe(j)%k)/sqrt(var%ss)

!!$        vaca(i)%regI=vaca(i)%regI+vaca(i)%fac%ss+ &
!!$             vaca(i)%fac%sc+vaca(i)%fac%cs+vaca(i)%fac%cc 
          if(i.ge.j) cycle
          vaca(i)%coul=-(vaca(i)%fac%ss+vaca(i)%fac%sc+vaca(i)%fac%cs+vaca(i)%fac%cc) 
          l=>var !rr calculated
          call calc_shortrange(epe(i)%k,epe(j)%k,i,j,vaca(i)%short,vaca(i)%gshort)
          vaca(i)%tot= vaca(i)%coul-vaca(i)%short
          eshort_vaccluster=eshort_vaccluster-vaca(i)%short
          ecoul_vaccluster=ecoul_vaccluster+vaca(i)%coul
!!$        write(output_epe, "(1x,'! vacancy  ',i2,' - vacancy  ',i2,' !',f12.4,  &
!!$             ' !',f12.4,' ! ',f12.4)")i,j,vaca(i)%coul,-vaca(i)%short,vaca(i)%tot
       enddo
!!$
!!$   do j=1,n_vacancies
!!$      IF(I.EQ.J) cycle
!!$      var%rr=dot_product(epe(I)%r-epe(J)%r,epe(I)%r-epe(J)%r)
!!$      vaca(i)%regI=vaca(i)%regI+q_ion(epe(I)%k)*q_ion(epe(J)%k)/SQRT(var%rr)
!!$      IF(I.GE.J) cycle
!!$      vaca(i)%coul=-q_ion(epe(I)%k)*q_ion(epe(J)%k)/SQRT(var%rr)
!!$      var%ss=var%rr
!!$      l=>var
!!$      CALL calc_shortrange(epe(I)%k,epe(J)%k,i,j,vaca(i)%short,vaca(i)%gshort)
!!$      vaca(i)%tot= vaca(i)%coul-vaca(i)%short
!!$      eshort_vaccluster=eshort_vaccluster-vaca(i)%short
!!$      ecoul_vaccluster=ecoul_vaccluster+vaca(i)%coul
!!$      WRITE(output_epe, "(1X,'! vacancy  ',I2,' - vacancy  ',I2,' !',F12.4,  &
!!$           ' !',F12.4,' ! ',F12.4)")I,J,vaca(i)%coul,-vaca(i)%short,vaca(i)%tot
!!$   enddo
     end subroutine vac_to_vac_contibs

     subroutine vac_to_latt_contribs

       vaca(i)%short=zero

       do j=n_vacancies+1,reg_i_n_ions

      ! dot products
          var%rr=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)
          var%ss=dot_product(epe(i)%s-epe(j)%s,epe(i)%s-epe(j)%s)
          var%cs=dot_product(epe(i)%c-epe(j)%s,epe(i)%c-epe(j)%s)
          var%sc=dot_product(epe(i)%s-epe(j)%c,epe(i)%s-epe(j)%c)
          var%cc=dot_product(epe(i)%c-epe(j)%c,epe(i)%c-epe(j)%c)

          l=>var !rr calculated
          call calc_shortrange(epe(i)%k,epe(j)%k,i,j,vaca(i)%eshort,vaca(i)%gshort)
          vaca(i)%short=vaca(i)%short-vaca(i)%eshort
          short_range_grad1(j,:)=short_range_grad1(j,:)+vaca(i)%gshort*(epe(i)%s-epe(j)%s)
          var%rc=dot_product(epe(i)%r-epe(j)%c, epe(i)%r-epe(j)%c)

!!$!      vaca(i)%regi=vaca(i)%regi+q_ion(epe(i)%k)*(q_ion(epe(j)%k)/sqrt(var%rr) &
!!$!           -q_shell(epe(j)%k)/var%sq_ss-q_nuclear(epe(j)%k)/sqrt(var%rc))
          var%sq_sc=sqrt(var%sc)
          var%sq_cs=sqrt(var%cs)
          var%sq_cc=sqrt(var%cc)

          vaca(i)%regI=vaca(i)%regI &
               - q_shell(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_ss+q_nuclear(epe(j)%k)/var%sq_sc) &
               - q_nuclear(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_cs+q_nuclear(epe(j)%k)/var%sq_cc)
          
          point_x_core_grad(j,:)=point_x_core_grad(j,:)-q_nuclear(epe(j)%k)*( &
               q_nuclear(epe(i)%k)/(var%sq_cc*var%cc) *(epe(i)%c-epe(j)%c)+ &
               q_shell(epe(i)%k)/(var%sq_sc*var%sc) *(epe(i)%s-epe(j)%c) ) 

          point_x_shell_grad(j,:)=point_x_shell_grad(j,:)-q_shell(epe(j)%k)*( &
               q_shell(epe(i)%k)/(var%ss3) *(epe(i)%s-epe(j)%s)+ &
               q_nuclear(epe(i)%k)/(var%sq_cs*var%cs)*(epe(i)%c-epe(j)%s) )
       enddo
   ! ** area    i  finished 
 
! **treat region 2a contributions
      DO J=reg_I_n_ions+1,reg_2a_n_ions  

          var%rr=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)

          var%ss=dot_product(epe(i)%s-epe(j)%s,epe(i)%s-epe(j)%s)
          var%cs=dot_product(epe(i)%c-epe(j)%s,epe(i)%c-epe(j)%s)
          var%sc=dot_product(epe(i)%s-epe(j)%c,epe(i)%s-epe(j)%c)
          var%cc=dot_product(epe(i)%c-epe(j)%c,epe(i)%c-epe(j)%c)
          var%ss_sr=dot_product(epe(i)%s-epe(j)%s,epe(j)%s-epe(j)%r)

          l=>var !var rr calculated
          call calc_shortrange(epe(I)%k,epe(J)%k,i,j,vaca(i)%eshort,vaca(i)%gshort)
          vaca(i)%short=vaca(i)%short-vaca(i)%eshort

          if(.not.fixed_reg2a_relcontribs) then
             eshort_lat_relaxation=eshort_lat_relaxation &
                  -0.5_r8_kind*vaca(i)%gshort*var%ss_sr
             ! contrib No 1 vacancy
          end if
       

          var%rc=dot_product(epe(i)%r-epe(J)%c,epe(i)%r-epe(J)%c)
          var%rs_sr=dot_product(epe(i)%r-epe(J)%s,epe(j)%s-epe(J)%r)
          var%rc_cr=dot_product(epe(i)%r-epe(j)%c,epe(j)%c-epe(j)%r)

! **monopol contribution        core and shell j with the vacancy

          var%sq_cc=sqrt(var%cc)
          var%sq_sc=sqrt(var%sc)
          var%sq_cs=sqrt(var%cs)
          vaca(i)%reg2a=vaca(i)%reg2a &
             !!$            + q_ion(epe(i)%k)*q_ion(epe(j)%k)/sqrt(var%rr) &
               - q_shell(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_ss+q_nuclear(epe(j)%k)/var%sq_sc) &
               - q_nuclear(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_cs+q_nuclear(epe(j)%k)/var%sq_cc)

          if(.not.fixed_reg2a_relcontribs) then
             ! **(5) dipol contribution, j-dipol - i-vacanvy
             ecoul_lat_relaxation=ecoul_lat_relaxation &
                  + 0.5_r8_kind*Q_ion(epe(I)%k) * &
                  (Q_SHELL(epe(J)%k)*var%rs_sr/(var%ss3) &
                  +Q_NUCLEAR(epe(J)%k)*var%rc_cr/(var%rc*SQRT(var%rc)))
          end if
       
       enddo
     ! **area 2a finished
 
     end subroutine vac_to_latt_contribs
!-------------------------------------------------

!-------------------------------------------------
     subroutine imp_3_body_1(e3bi,reg,ref)
       ! use calc3c_switches, only: print_epe

!    called in qm_interfaced_mode from  defect contib
!    calculates energy and gradients 
!    depends on r_sh_imp and  reg_reference(if(use_ref_data))

       logical, optional,intent(in):: reg
       logical, optional,intent(in):: ref
       real(kind=r8_kind), intent(out) :: e3bi
       real(kind=r8_kind) :: e1(3),e2(3),g3bi,theta,cos_th,sin_th,rss1,rss2,scal1
       real(kind=r8_kind) :: r1(3),r2(3),r3(3)
       real(kind=r8_kind) :: deg2rad
       integer(kind=i4_kind) :: i,j,k,ia_1,ia_2,ia_3,i1,j1,k1,l,ii
       integer(kind=i4_kind) :: at1,at2,at3
       integer(kind=i4_kind) :: num_3b_terms
       real(kind=r8_kind) :: g3bi_sum

       deg2rad=pi/180.0_r8_kind
       e3bi=zero
       g3bi_sum=zero

       if( n_vacancies == 0 ) return

         if(use_ref_data.and.present(reg)) then
          print*,'reg reference 3b  treated'
         elseif(use_ref_data.and.present(ref)) then
          print*,'qm reference 3b  treated'
         else
          print*,'var 3b  treated'
         endif

       num_3b_terms=0
       fst: do i=1,n_tetrahedrons
          do l=1,5
             if (tetra_atoms(i,l) <= n_vacancies) goto 1
             ! at least one atom belong to cluster
          enddo
          cycle fst
1         do l=1,5                                      !!!!!AS
             !if (tetra_atoms(i,l) > n_vacancies) goto 2 !!!!!AS
             ! at least one atom out of cluster         !!!!!AS
          enddo                                         !!!!!AS
2         ia_1=tetra_atoms(i,1)
           i1=epe(ia_1)%k !type of ia_1
          if(explicit_coupling) then
          if(use_ref_data.and.present(reg)) then
           i1=epe(ia_1)%k !type of ia_1
          elseif(use_ref_data.and.present(ref)) then
           i1=epe(ia_1)%k
          else
           i1=epe(ia_1)%ec
          endif
         endif
          
          scnd: do j=2,4
             ia_2=tetra_atoms(i,j)
             if(ia_2 == 0) cycle scnd
             j1=epe(ia_2)%k
          if(explicit_coupling) then
          if(use_ref_data.and.present(reg)) then
           j1=epe(ia_2)%k !type of ia_1
          elseif(use_ref_data.and.present(ref)) then
           j1=epe(ia_2)%k
          else
           j1=epe(ia_2)%ec
          endif
         endif
             thrd: do k=j+1,5
                ia_3=tetra_atoms(i,k)
                if(ia_3 == 0) cycle thrd
                k1=epe(ia_3)%k
          if(explicit_coupling) then
          if(use_ref_data.and.present(reg)) then
           k1=epe(ia_3)%k !type of ia_1
          elseif(use_ref_data.and.present(ref)) then
           k1=epe(ia_3)%k
          else
           k1=epe(ia_3)%ec
          endif
         endif
                if(ia_1 > n_vacancies.and.ia_2 >n_vacancies .and.ia_3 > n_vacancies) cycle thrd

                if(ia_1 > n_vacancies) then
                   r1=epe(ia_1)%s(:)
                else
                   ii=epe(ia_1)%gc
                    if(use_ref_data.and.present(reg)) then
                      r1=reg_reference(ii)%rs
                    elseif(use_ref_data.and.present(ref)) then
                      r1=epe_reference(ii)%rs
                   else
                      r1=r_sh_imp(ii,:)
                   endif
                endif
                if(ia_2 > n_vacancies) then
                   r2=epe(ia_2)%s(:)
                else
                   ii=epe(ia_2)%gc
                   if(use_ref_data.and.present(reg)) then
                      r2=reg_reference(ii)%rs
                   elseif(use_ref_data.and.present(ref)) then
                      r2=epe_reference(ii)%rs
                   else
                      r2=r_sh_imp(ii,:)
                   endif
                endif
                if(ia_3 > n_vacancies) then
                   r3=epe(ia_3)%s(:)
                else
                   ii=epe(ia_3)%gc
                   if(use_ref_data.and.present(reg)) then
                      r3=reg_reference(ii)%rs
                   elseif(use_ref_data.and.present(ref)) then
                      r3=epe_reference(ii)%rs
                   else
                      r3=r_sh_imp(ii,:)
                   endif
                endif

                 if((use_ref_data.and..not.present(reg)).and. &
                    (ia_1.le.n_vacancies.and.ia_2.le.n_vacancies.and.ia_3.le.n_vacancies)) cycle thrd
               

                rss1=dot_product(r2-r1,r2-r1)
                rss2=dot_product(r3-r1,r3-r1)
                scal1=dot_product(r2-r1,r3-r1)

                rss1=sqrt(rss1)
                rss2=sqrt(rss2)
                cos_th=scal1/(rss1*rss2)
                theta=acos(cos_th)
                num_3b_terms=num_3b_terms+1
!                if(print_epe) print*, 'theta',num_3b_terms,i1, theta/deg2rad,theta_0(j1,i1,k1),ia_1,ia_2,ia_3
                if(ia_1> n_vacancies) at1=ia_1
                if(ia_1.le. n_vacancies) at1= epe(ia_1)%gc
                if(ia_2> n_vacancies) at2= ia_2
                if(ia_2.le. n_vacancies) at2= epe(ia_2)%gc
                if(ia_3> n_vacancies) at3= ia_3
                if(ia_3.le. n_vacancies) at3= epe(ia_3)%gc
                sin_th=sin(theta)

               if(use_ref_data.and.present(ref)) then
! print*, 'theta',num_3b_terms,i1, theta/deg2rad,theta_0(j1,i1,k1),ki(j1,i1,k1),at1,at2,at3
                e3bi=e3bi-0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2
                g3bi=-ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
               elseif(use_ref_data.and.present(reg)) then
! print*, 'theta',num_3b_terms,i1, theta/deg2rad,theta_0(j1,i1,k1),ki(j1,i1,k1),at1,at2,at3
                e3bi=e3bi+0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2
                g3bi=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
               else
! print*, 'theta',num_3b_terms,i1, theta/deg2rad,ec%theta_0(j1,i1,k1),ec%ki(j1,i1,k1),at1,at2,at3
                e3bi=e3bi+0.5*ec%ki(j1,i1,k1)*(theta-ec%theta_0(j1,i1,k1)*deg2rad)**2
                g3bi=ec%ki(j1,i1,k1)*(theta-ec%theta_0(j1,i1,k1)*deg2rad)
               endif
                g3bi_sum=g3bi_sum+g3bi
                
                e1=(r2-r1)/rss1
                e2=(r3-r1)/rss2

                if (ia_3 > n_vacancies) then
                 short_range_grad1(ia_3,:)=short_range_grad1(ia_3,:)+g3bi*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                else
                 ii=epe(ia_3)%gc
                 impurity_shell_grad(ii,:)=impurity_shell_grad(ii,:)+g3bi*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                endif

                if (ia_2 > n_vacancies) then
                 short_range_grad1(ia_2,:)=short_range_grad1(ia_2,:)+g3bi*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
                else
                 ii=epe(ia_2)%gc
                 impurity_shell_grad(ii,:)=impurity_shell_grad(ii,:)+g3bi*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
                endif

                if (ia_1 > n_vacancies) then
                   short_range_grad1(ia_1,:)=short_range_grad1(ia_1,:)+ &
                        g3bi*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
                else
                   ii=epe(ia_1)%gc
                   impurity_shell_grad(ii,:)=impurity_shell_grad(ii,:)+ &
                        g3bi*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
                endif

             enddo thrd!k=j+1,4
          enddo scnd!j=1,3
       enddo fst!i=1,n_tetra_atoms
       DPRINT 'g3bi_sum ',g3bi_sum

     end subroutine imp_3_body_1


     subroutine imp_3_body(e3bi)

!    called for .not.qm_interfaced_mode

       real(kind=r8_kind), intent(out) :: e3bi
       type(epe_imp_tet) :: tetra_imp

       integer(kind=i4_kind) :: i,j,k,i1,jj

       e3bi=zero
       if( n_impurities == 0 ) return
       fst: do i=1,n_types_central_atoms_3body_im
          i1=types_im(i,1)

          ! impurity-impurity + impurity-lattice
          jj=0
          scnd: do j=1,n_impurities

             if(i1 == type_impurity(j)) then
                call select_centers(tetra_imp,j,'imp',i)
                call energy_and_grad(e3bi,tetra_imp)
             else
                jj=jj+1
             endif
          enddo scnd

          ! lattice-impurity
          if(jj == n_impurities) then
             thrd: do j=n_vacancies+1,reg_I_n_ions
                if(i1 ==  epe(j)%k) then
                   call select_centers(tetra_imp,j,'epe',i)
                   do k=2,5
                      if(tetra_imp%center(k) == 'imp') goto 1
                   enddo
                   cycle thrd
1                  call energy_and_grad(e3bi,tetra_imp)
                endif
             enddo thrd
          endif
       enddo fst

     end subroutine imp_3_body

     subroutine energy_and_grad(e3bi,tetra_imp)

       real(kind=r8_kind), intent(inout) :: e3bi
       type(epe_imp_tet),intent(in) :: tetra_imp

       real(kind=r8_kind) :: gg(3),e1(3),e2(3),g3bi,theta,cos_th,sin_th,rss1,rss2,scal1
       real(kind=r8_kind) :: deg2rad
       integer(kind=i4_kind) :: ia_1,ia_2,ia_3,i1,j1,k1,j,k

       deg2rad=pi/180.0_r8_kind
       ia_1=tetra_imp%number(1)
       if(tetra_imp%center(1) == 'imp') then 
          i1=type_impurity(ia_1)
       else
          i1=epe(ia_1)%k
       endif

       do j=2,4
          ia_2=tetra_imp%number(j)
          if(ia_2 == 0) cycle
          if(tetra_imp%center(j) == 'imp') then 
             j1=type_impurity(ia_2)
          else
             j1=epe(ia_2)%k
          endif

          do k=j+1,5
             ia_3=tetra_imp%number(k)
             if(ia_3 == 0) cycle
             if(tetra_imp%center(k) == 'imp') then 
                k1=type_impurity(ia_3)
             else
                k1=epe(ia_3)%k
             endif
             if(tetra_imp%center(1) /='imp'.and. tetra_imp%center(j) /='imp'.and. &
                  tetra_imp%center(k) /= 'imp') cycle

             rss1=dot_product(tetra_imp%coord(j,:)-tetra_imp%coord(1,:), &
                  tetra_imp%coord(j,:)-tetra_imp%coord(1,:))
             rss2=dot_product(tetra_imp%coord(k,:)-tetra_imp%coord(1,:), &
                  tetra_imp%coord(k,:)-tetra_imp%coord(1,:))
             scal1=dot_product(tetra_imp%coord(j,:)-tetra_imp%coord(1,:), &
                  tetra_imp%coord(k,:)-tetra_imp%coord(1,:))

             rss1=sqrt(rss1)
             rss2=sqrt(rss2)
             cos_th=scal1/(rss1*rss2)
             theta=acos(cos_th)
             sin_th=sin(theta)

             e3bi=e3bi+0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2

             g3bi=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
             e1=(tetra_imp%coord(j,:)-tetra_imp%coord(1,:))/rss1
             e2=(tetra_imp%coord(k,:)-tetra_imp%coord(1,:))/rss2

             gg=g3bi*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
             if(tetra_imp%center(1)=='imp') then
                impurity_shell_grad(ia_1,:)=impurity_shell_grad(ia_1,:)+gg
             else
                short_range_grad1(ia_1,:)=short_range_grad1(ia_1,:)+gg
             endif
             gg=g3bi*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
             if(tetra_imp%center(j)=='imp') then
                impurity_shell_grad(ia_2,:)=impurity_shell_grad(ia_2,:)+gg
             else
                short_range_grad1(ia_2,:)=short_range_grad1(ia_2,:)+gg
             endif
             gg=g3bi*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
             if(tetra_imp%center(k)=='imp') then
                impurity_shell_grad(ia_3,:)=impurity_shell_grad(ia_3,:)+gg
             else
                short_range_grad1(ia_3,:)=short_range_grad1(ia_3,:)+gg
             endif
          enddo
       enddo

     end subroutine energy_and_grad

     subroutine select_centers(tetra_imp,num,who_am_i,n_type_3b)

       type(epe_imp_tet),intent(out) :: tetra_imp
       integer(kind=i4_kind),intent(in) :: num,n_type_3b
       character(len=3),intent(in) :: who_am_i

       real(kind=r8_kind) :: dist,rrr(3)
       real(kind=r8_kind), dimension(4) :: buf_dist
       integer(kind=i4_kind), dimension(4) :: buf_index
       character(len=3), dimension(4) :: buf_type
       integer(kind=i4_kind) :: max_ind(1)
       real(kind=r8_kind) :: max_val
       integer(kind=i4_kind) :: i,j,k

       tetra_imp%coord=0.0_r8_kind
       tetra_imp%number(1)=num
       if(who_am_i == 'imp') then
          tetra_imp%center(1)='imp'
          if(use_ref_data) then
             tetra_imp%coord(1,:)=reg_reference(num)%rs
             rrr=reg_reference(num)%rs
          else
             tetra_imp%coord(1,:)=r_sh_imp(num,:)
             rrr=r_sh_imp(num,:)
          endif
       else
          tetra_imp%center(1)='epe'
          tetra_imp%coord(1,:)=epe(num)%r
          rrr=epe(num)%r
       endif

       buf_dist=0.0_r8_kind
       buf_index=0
       fst: do i=1,n_impurities
          do j=2,5
             if(type_impurity(i)==types_im(n_type_3b,j)) goto 1
          enddo
          cycle fst

1         if(use_ref_data) then
             dist=dot_product(reg_reference(i)%rs-rrr, &
                  reg_reference(i)%rs-rrr)
!                  reg_reference(j)%rs-rrr) ! bug fixed
          else
             dist=dot_product(r_sh_imp(i,:)-rrr, &
                  r_sh_imp(i,:)-rrr)
          endif
          dist=sqrt(dist)
          if(dist > r3b_im(n_type_3b)) cycle fst

          do k=1,4
             if(buf_dist(k) == 0.0_r8_kind) then
                buf_dist(k)=dist
                buf_index(k)=i
                buf_type(k)='imp'
                cycle fst
             endif
          enddo
          max_ind=maxloc(buf_dist)
          max_val=maxval(buf_dist)
          if(dist < max_val) then
             buf_dist(max_ind)=dist
             buf_index(max_ind)=i
             buf_type(max_ind)='imp'
          endif
       enddo fst

       scnd: do i=n_vacancies+1,reg_I_n_ions
          if(i == num) cycle scnd
          do j=2,5
             if(epe(i)%k == types_im(n_type_3b,j)) goto 2
          enddo
          cycle scnd

2         dist=dot_product(epe(i)%r-rrr, &
               epe(i)%r-rrr)
          dist=sqrt(dist)
          if(dist > r3b_im(n_type_3b)) cycle scnd

          do k=1,4
             if(buf_dist(k) == 0.0_r8_kind) then
                buf_dist(k)=dist
                buf_index(k)=i
                buf_type(k)='epe'
                cycle scnd
             endif
          enddo
          max_ind=maxloc(buf_dist)
          max_val=maxval(buf_dist)
          if(dist < max_val) then
             buf_dist(max_ind)=dist
             buf_index(max_ind)=i
             buf_type(max_ind)='epe'
          endif
       enddo scnd
       
       tetra_imp%number(2:5)=buf_index
       tetra_imp%center(2:5)=buf_type
       do i=2,5
          j=tetra_imp%number(i)
          if(tetra_imp%center(i) == 'imp') then
             if(use_ref_data) then
                tetra_imp%coord(i,:)=reg_reference(j)%rs
             else
                tetra_imp%coord(i,:)=r_sh_imp(j,:)
             endif
          else
             tetra_imp%coord(i,:)=epe(j)%s
          endif
       enddo

     end subroutine select_centers

     subroutine vac_3_body(e3bv)

       real(kind=r8_kind), intent(out) :: e3bv
       real(kind=r8_kind) :: e1(3),e2(3),g3bv,theta,cos_th,sin_th,rss1,rss2,scal1
       real(kind=r8_kind) :: deg2rad

       integer(kind=i4_kind) :: i,j,k,ia_1,ia_2,ia_3,i1,j1,k1,l

       deg2rad=pi/180.0_r8_kind
       e3bv=zero
       if( n_vacancies == 0 .or. n_types_central_atoms_3body == 0 ) return
!!$       fst:  do i=1,last_ind-first_ind+1
!!$       fst: do i=first_ind,last_ind
       fst: do i=1,n_tetrahedrons
          do l=1,5
             if (tetra_atoms(i,l) <= n_vacancies) goto 1
          enddo
          cycle fst
1         ia_1=tetra_atoms(i,1)
          i1=epe(ia_1)%k
          scnd: do j=2,4
             ia_2=tetra_atoms(i,j)
             if(ia_2 == 0) cycle scnd
             j1=epe(ia_2)%k
             thrd: do k=j+1,5
                ia_3=tetra_atoms(i,k)
                if(ia_3 == 0) cycle thrd
                if(ia_1 > n_vacancies.and.ia_2 >n_vacancies .and.ia_3 > n_vacancies) cycle thrd
                k1=epe(ia_3)%k
   
                rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
                     epe(ia_2)%s(:)-epe(ia_1)%s(:))
!!$                if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies) &
!!$                     rss1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_2)%r(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies) &
!!$                     rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_2)%s(:)-epe(ia_1)%r(:))
!!$                   
!!$                if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies) &
!!$                     rss1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_2)%r(:)-epe(ia_1)%s(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_2 > n_vacancies) &
!!$                     rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_2)%s(:)-epe(ia_1)%s(:))

                rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%s(:), &
                     epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$                if(ia_1 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     rss2=dot_product(epe(ia_3)%r(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                     rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$                 
!!$                if(ia_1 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     rss2=dot_product(epe(ia_3)%r(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                     rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%s(:))

                scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
                     epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$                if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_2 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$                if(ia_1 > n_vacancies.and.ia_2 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                     scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                     epe(ia_3)%s(:)-epe(ia_1)%s(:))
      
                rss1=sqrt(rss1)
                rss2=sqrt(rss2)
                cos_th=scal1/(rss1*rss2)
                theta=acos(cos_th)
                sin_th=sin(theta)

                e3bv=e3bv-0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2

                g3bv=-ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
                
                e1=(epe(ia_2)%s-epe(ia_1)%s)/rss1
!!$                if (ia_1 <= n_vacancies.and.ia_2 <= n_vacancies ) &
!!$                     e1=(epe(ia_2)%r-epe(ia_1)%r)/rss1
!!$
!!$                if (ia_1 <= n_vacancies.and.ia_2 > n_vacancies ) &
!!$                     e1=(epe(ia_2)%s-epe(ia_1)%r)/rss1
!!$
!!$                if (ia_1 > n_vacancies.and.ia_2 <= n_vacancies ) &
!!$                     e1=(epe(ia_2)%r-epe(ia_1)%s)/rss1
!!$
!!$                if (ia_1 > n_vacancies.and.ia_2 > n_vacancies ) &
!!$                     e1=(epe(ia_2)%s-epe(ia_1)%s)/rss1

                e2=(epe(ia_3)%s-epe(ia_1)%s)/rss2
!!$                if (ia_1 <= n_vacancies.and.ia_3 <= n_vacancies ) &
!!$                     e2=(epe(ia_3)%r-epe(ia_1)%r)/rss2
!!$
!!$                if (ia_1 <= n_vacancies.and.ia_3 > n_vacancies ) &
!!$                     e2=(epe(ia_3)%s-epe(ia_1)%r)/rss2
!!$
!!$                if (ia_1 > n_vacancies.and.ia_3 <= n_vacancies ) &
!!$                     e2=(epe(ia_3)%r-epe(ia_1)%s)/rss2
!!$
!!$                if (ia_1 > n_vacancies.and.ia_3 > n_vacancies ) &
!!$                     e2=(epe(ia_3)%s-epe(ia_1)%s)/rss2
                
                if (ia_3 > n_vacancies) &
                     short_range_grad1(ia_3,:)=short_range_grad1(ia_3,:)+g3bv*(cos_th*e2(:)-e1(:))/(rss2*sin_th)

                if (ia_2 > n_vacancies) &
                     short_range_grad1(ia_2,:)=short_range_grad1(ia_2,:)+g3bv*(cos_th*e1(:)-e2(:))/(rss1*sin_th)

                if (ia_1 > n_vacancies) &
                     short_range_grad1(ia_1,:)=short_range_grad1(ia_1,:)+ &
                     g3bv*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
             enddo thrd!k=j+1,4
          enddo scnd!j=1,3
       enddo fst!i=1,n_tetra_atoms

     end subroutine vac_3_body

     subroutine send_intermediate_vac
       !  Purpose: Sends a part of results to master. 
       !  Called by slave.
       !------------ Modules used ----------------------------------
       !------------ Declaration of local variables -----------------
       integer  :: status
       !------------ Executable code --------------------------------

       allocate( send_intermediate_array(n_vacancies,3) ,stat=status)
       if(status.ne.0) call error_handler("allocate send_intermediate_array failed(1)")

       send_intermediate_array=zero
       send_intermediate_array(i_1:i_2,1) = vaca(i_1:i_2)%coul
       send_intermediate_array(i_1:i_2,2) = vaca(i_1:i_2)%short
       send_intermediate_array(i_1:i_2,3) = vaca(i_1:i_2)%tot

       call comm_init_send(comm_master_host, msgtag_vac_intermediate)
       call commpack(send_intermediate_array(1,1),3*n_vacancies, 1, status)
       if(status.ne.0) & 
            call error_handler("epe_send_intermediate: error [1]")
       call commpack(e_vac_3b,status)
       if(status.ne.0) & 
            call error_handler("e_vac_3b - pack: error [1]")

       call comm_send()

       deallocate(send_intermediate_array,stat=status)
       if(status.ne.0) call error_handler("deallocate  send_intermediate_array failed")

     end subroutine send_intermediate_vac

     subroutine receive_intermediate_vac
       !  Purpose: Receives a part of results from slaves. 
       !  Called by master.
       !------------ Modules used ----------------------------------
       !------------ Declaration of local variables -----------------
       integer             :: n_received, status
       logical             :: all_results_received
       real(kind=r8_kind)  :: e_var_sl
       !------------ Executable code --------------------------------

       allocate( send_intermediate_array(n_vacancies,3) ,stat=status)
       if(status.ne.0) call error_handler("allocate send_intermediate_array failed")
       send_intermediate_array=zero

       all_results_received=.false.
       n_received=0
       do while(.not.all_results_received )
          do while ( .not. comm_save_recv_nonblocking(comm_all_other_hosts, &
               msgtag_vac_intermediate) )
          end do
          n_received=n_received+1 
          call communpack(send_intermediate_array(1,1),3*n_vacancies, 1, status)
          if(status.ne.0) & 
               call error_handler("epe_receive_intermediate: error [1]")
          call communpack(e_var_sl, status)
          if(status.ne.0) & 
               call error_handler("e_var_sl - unpack: error [1]")

          vaca(1:n_vacancies)%coul  = vaca(1:n_vacancies)%coul  +  send_intermediate_array(:,1)
          vaca(1:n_vacancies)%short = vaca(1:n_vacancies)%short +  send_intermediate_array(:,2)
          vaca(1:n_vacancies)%tot   = vaca(1:n_vacancies)%tot   +  send_intermediate_array(:,3)

          e_vac_3b=e_vac_3b+e_var_sl

          all_results_received=n_received==comm_get_n_processors()-1
       end do

       deallocate( send_intermediate_array, stat=status )
       if(status.ne.0) call error_handler("deallocate  send_intermediate_array failed")

     end subroutine receive_intermediate_vac

     subroutine send_intermediate_imp
       !  Purpose: Sends a part of results to master. 
       !  Called by slave.
       !------------ Modules used ----------------------------------
       !------------ Declaration of local variables -----------------
       integer  :: status
       !------------ Executable code --------------------------------

       allocate( send_intermediate_array(n_impurities+19,3) ,& 
            reg_2a_ec_send(n_impurities) ,stat=status)
       if(status.ne.0) call error_handler(" allocate send_intermediate_array failed(2)")
       send_intermediate_array(1:n_impurities+19,:) = zero
       if(n_impurities > 0) then
          reg_2a_ec_send(1:n_impurities)               = zero
          send_intermediate_array(i_1:i_2,1) = clus(i_1:i_2)%coul
          send_intermediate_array(i_1:i_2,2) = clus(i_1:i_2)%short
          send_intermediate_array(i_1:i_2,3) = clus(i_1:i_2)%tot
          reg_2a_ec_send(1:n_impurities)     = reg_2a_ec(1:n_impurities)
       endif
!
!       ... and partial energy contributions to print on master !
!
       send_intermediate_array(n_impurities+1, 1) = eshort_vaccluster
       send_intermediate_array(n_impurities+2, 1) = eshort_vacepe
       send_intermediate_array(n_impurities+3, 1) = eshort_epecluster
       send_intermediate_array(n_impurities+4, 1) = eshort_epecluster_epe
       send_intermediate_array(n_impurities+5, 1) = eshort_lat_relaxation
       send_intermediate_array(n_impurities+6, 1) = eshort_lat_relcorr
       send_intermediate_array(n_impurities+7, 1) = eshort_reg_I
       send_intermediate_array(n_impurities+8, 1) = eshort_reg2a  
       send_intermediate_array(n_impurities+9, 1) = ecoul_vaccluster
       send_intermediate_array(n_impurities+10,1) = ecoul_vaccluster_epe
       send_intermediate_array(n_impurities+11,1) = ecoul_epecluster
       send_intermediate_array(n_impurities+12,1) = ec_cluster_epe 
       send_intermediate_array(n_impurities+13,1) = ecoul_lat_relaxation
       send_intermediate_array(n_impurities+14,1) = ecoul_lat_relcorr
       send_intermediate_array(n_impurities+15,1) = ec_cluster_reg2a
       send_intermediate_array(n_impurities+16,1) = diffpg_ec_ecref
       send_intermediate_array(n_impurities+17,1) = defect_energy_short
       send_intermediate_array(n_impurities+18,1) = defect_energy_coul
       send_intermediate_array(n_impurities+19,1) = ecrr

       call comm_init_send(comm_master_host, msgtag_imp_intermediate)
       call commpack(send_intermediate_array(1,1),3*(n_impurities+19), 1, status)
       if(status.ne.0) & 
            call error_handler("epe_send_intermediate: error [1a]")
       if(n_impurities > 0) then
          call commpack(reg_2a_ec_send, n_impurities, 1, status)
          if(status.ne.0) & 
               call error_handler("epe_send_intermediate: error [2a]")
       endif
       call comm_send()

       deallocate(send_intermediate_array,&
            reg_2a_ec_send ,stat=status)
       if(status.ne.0) call error_handler("deallocate send_intermediate_array failed")

     end subroutine send_intermediate_imp

     subroutine receive_intermediate_imp
       !  Purpose: Receives a part of results from slaves. 
       !  Called by master.
       !------------ Modules used ----------------------------------
       !------------ Declaration of local variables -----------------
       integer             :: n_received, status
       logical             :: all_results_received
       !------------ Executable code --------------------------------

       allocate( send_intermediate_array(n_impurities+19,3),reg_2a_ec_send(n_impurities),&
            stat=status )
       if(status.ne.0) call error_handler("allocate send_intermediate_array failed(3)")
       send_intermediate_array=zero

       all_results_received=.false.
       n_received=0
       do while(.not.all_results_received )
          do while ( .not. comm_save_recv_nonblocking(comm_all_other_hosts, &
               msgtag_imp_intermediate) )
          end do
          n_received=n_received+1 
          call communpack(send_intermediate_array(1,1),3*(n_impurities+19) , 1, status)
          if(status.ne.0) & 
               call error_handler("epe_receive_intermediate: error [1a]")
          if(n_impurities > 0) then
             call communpack(reg_2a_ec_send, n_impurities, 1, status)
             if(status.ne.0) & 
                  call error_handler("epe_receive_intermediate: error [2a]")
          endif

          if(n_impurities > 0) then
             clus(1:n_impurities)%coul  = clus(1:n_impurities)%coul+ send_intermediate_array(1:n_impurities,1)
             clus(1:n_impurities)%short = clus(1:n_impurities)%short+send_intermediate_array(1:n_impurities,2)
             clus(1:n_impurities)%tot   = clus(1:n_impurities)%tot + send_intermediate_array(1:n_impurities,3)
          endif
          eshort_vaccluster     = eshort_vaccluster      +  send_intermediate_array(n_impurities+1, 1) 
          eshort_vacepe         = eshort_vacepe          +  send_intermediate_array(n_impurities+2, 1)
          eshort_epecluster     = eshort_epecluster      +  send_intermediate_array(n_impurities+3, 1)
          eshort_epecluster_epe = eshort_epecluster_epe  +  send_intermediate_array(n_impurities+4, 1)
          eshort_lat_relaxation = eshort_lat_relaxation  +  send_intermediate_array(n_impurities+5, 1) 
          eshort_lat_relcorr    = eshort_lat_relcorr     +  send_intermediate_array(n_impurities+6, 1)
          eshort_reg_I          = eshort_reg_I           +  send_intermediate_array(n_impurities+7, 1) 
          eshort_reg2a          = eshort_reg2a           +  send_intermediate_array(n_impurities+8, 1)   
          ecoul_vaccluster      = ecoul_vaccluster       +  send_intermediate_array(n_impurities+9, 1)
          ecoul_vaccluster_epe  = ecoul_vaccluster_epe   +  send_intermediate_array(n_impurities+10,1)
          ecoul_epecluster      = ecoul_epecluster       +  send_intermediate_array(n_impurities+11,1) 
          ec_cluster_epe        = ec_cluster_epe         +  send_intermediate_array(n_impurities+12,1) 
          ecoul_lat_relaxation  = ecoul_lat_relaxation   +  send_intermediate_array(n_impurities+13,1) 
          ecoul_lat_relcorr     = ecoul_lat_relcorr      +  send_intermediate_array(n_impurities+14,1) 
          ec_cluster_reg2a      = ec_cluster_reg2a       +  send_intermediate_array(n_impurities+15,1) 
          diffpg_ec_ecref       = diffpg_ec_ecref        +  send_intermediate_array(n_impurities+16,1)
          defect_energy_short   = defect_energy_short    +  send_intermediate_array(n_impurities+17,1) 
          defect_energy_coul    = defect_energy_coul     +  send_intermediate_array(n_impurities+18,1)
          ecrr                  = ecrr                   +  send_intermediate_array(n_impurities+19,1)
          if(n_impurities > 0) then
             reg_2a_ec(1:n_impurities) = reg_2a_ec(1:n_impurities) + reg_2a_ec_send(1:n_impurities)
          endif

          all_results_received=n_received==comm_get_n_processors()-1
       end do

       deallocate( send_intermediate_array , reg_2a_ec_send, stat=status )
       if(status.ne.0) call error_handler("deallocate send_intermediate_array failed")

     end subroutine receive_intermediate_imp

   END SUBROUTINE calc_def_contributions

   subroutine epe_send_shells_and_cores(msg_type)
   !  Purpose: sends coordinates of shells and cores to be
   !           used in latt.grad. calculations to all slaves. 
   !           Called by master (subroutine main_epe)
   !           
   !------------ Modules used ----------------------------------    
     implicit none
     integer, intent(in) :: msg_type
     !------------ Declaration of local variables ----------------
     integer             :: status
     !------------ Declaration of subroutines used ---------------
     external error_handler
     !------------ Executable code -------------------------------
 
     select case(msg_type)
     case(0) ! init_latt_gradients
        call comm_init_send(comm_all_other_hosts, msgtag_epe_grads_init)
     case(1) ! cons_lattice_gradients
        call comm_init_send(comm_all_other_hosts, msgtag_epe_grads_cons)
     case(2)
        call comm_init_send(comm_all_other_hosts, msgtag_epe_send_only)
     case(3)
        call comm_init_send(comm_all_other_hosts, msgtag_epe_send_only)
     end select
  
! pack shells and cores to be sent

     call commpack(r_sh_ion(1,1),3*n_gen_ions,1,status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [1]")  
     call commpack(r_nuc_ion(1,1),3*n_gen_ions,1,status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [2]")   
     
     call commpack(r_sh_imp(1,1), 3*ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [3]")
     call commpack(r_nuc_imp(1,1), 3*ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [4]")
     call commpack(use_epe_pgdata,    status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [5]")
     call commpack(pk_impurity(1), ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_send_shells_and_cores : error [6]")     
     call comm_send()
   end subroutine epe_send_shells_and_cores
!*************************************************************

!*************************************************************
   subroutine epe_receive_shells_and_cores()
    !  Purpose: Complementary subroutine to
    !           epe_send_shells_and_cores 
    !           Called by slaves (subroutine main_epe)
    !           
    !------------ Modules used ----------------------------------    
     implicit none

     !------------ Declaration of local variables -----------------
     integer             :: status
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------

! unpack shells and cores received

     call communpack(r_sh_ion(1,1),3*n_gen_ions,1,status)
!!$     print*,'epe_receive_shells_and_cores epe%s :', sum(r_sh_ion(:reg_2a_n_ions,1)**2 + &
!!$          r_sh_ion(:reg_2a_n_ions,2)**2 + r_sh_ion(:reg_2a_n_ions,3)**2 )
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [1]")  
     call communpack(r_nuc_ion(1,1),3*n_gen_ions,1,status)
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [2]")  
     call communpack(r_sh_imp(1,1), 3*ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [3]") 
     call communpack(r_nuc_imp(1,1), 3*ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [4]") 
     call communpack(use_epe_pgdata,    status)
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [5]") 
     call communpack(pk_impurity(1), ndpt, 1, status)
     if(status.ne.0) call error_handler &
          ("epe_receive_shells_and_cores : error [6]")
   end subroutine epe_receive_shells_and_cores
!*************************************************************

!*************************************************************
   SUBROUTINE lattice_grad_gopt(XS,YS,ZS,XC,YC,ZC,DSX,DSY,DSZ, &
        DCX,DCY,DCZ,DX,DY,DZ)
! **calculation of relaxation gradients of lattice ions 
! **(DS - shell, DC - core)

     use epecom_module
     use culon_module

     real(kind=r8_kind),intent(out), &
          dimension(reg_I_n_ions) :: DSX,DSY,DSZ,DCX,DCY,DCZ 
     real(kind=r8_kind),intent(out),dimension(reg_I_n_ions) :: DX,DY,DZ 
     real(kind=r8_kind),intent(in),dimension(n_gen_ions) :: XS,YS,ZS,XC,YC,ZC
     real(kind=r8_kind),dimension(3) :: gs,gc,gpsis,gpsic
     real(kind=r8_kind) :: DIST,F1,ERF1,RRR2,RSSX,RSSY,RSSZ,CSIX,CSIY,CSIZ
     real(kind=r8_kind) :: RSS2,RSSCSI,RSS1,RSSR2,RSSR8,RSSD,ESS1,ESS2,RSSRO
     real(kind=r8_kind) :: BRO,OSO,OSS,RSCX,RSCY,RSCZ,RCSX,RCSY,RCSZ,RCCX,RCCY
     real(kind=r8_kind) :: RCCZ,RSRX,RSRY,RSRZ,RCRX,RCRY,RCRZ,CCIX,CCIY,CCIZ
     real(kind=r8_kind) :: RSC2,RCS2,RCC2,RSR2,RCR2,RCSCSI,RSCCCI,RCCCCI,RSS3
     real(kind=r8_kind) :: RCS3,RSC3,RCC3,RSR1,RCR1,OSC,OSR,OCS,OCC,OCR,OSSCSI
     real(kind=r8_kind) :: OCSCSI,OSCCCI,OCCCCI,DISX,DISY,DISZ,DICX,DICY,DICZ
     real(kind=r8_kind) :: DIS3X,DIS3Y,DIS3Z,DIC3X,DIC3Y,DIC3Z,RRR(3)
     real(kind=r8_kind) :: RRSX,RRSY,RRSZ,RRCX,RRCY,RRCZ,RRS2,RRC2,RRS1,RRC1
     real(kind=r8_kind) :: OSX,OSY,OSZ,ORS,ORC,OSSX,OSSY,OSSZ,OSCX,OSCY,OSCZ
     real(kind=r8_kind) :: OCSX,OCSY,OCSZ,OCCX,OCCY,OCCZ,CMCSX,CMCSY,CMCSZ
     real(kind=r8_kind) :: CSIS2,CSIC2,ET2S,ET2C,CSIS,OEDS,CSIC,OEDC
     real(kind=r8_kind) :: GIS2,GIC2,PSIS,PSIC,PI2S,PI2C,OEPKX,OEPKY,OEPKZ
     real(kind=r8_kind) :: est,psi,gig,res,gcel,gcelpi,qns,qnc,prsin,prcos
     real(kind=r8_kind) :: theta,scal,g3b,cos_th,sin_th
     real(kind=r8_kind),dimension(3) :: e1,e2
     real(kind=r8_kind) :: deg2rad
     integer(kind=i4_kind) :: i,io,ig,jo,k,n,j,i1,j1,k1,ia_1,ia_2,ia_3,icl,m,ic

     DIST=7.d0
     F1=0.026d0
     ERF1=4.d0/3.d0*ERROR_FUNCTION_PARAMETER**3/PIS
     DX=zero
     DY=zero
     DZ=zero
     DSX=zero
     DSY=zero
     DSZ=zero
     DCX=zero
     DCY=zero
     DCZ=zero

     if(n_types_central_atoms_3body > 0 ) then
        deg2rad=pi/180.0_r8_kind
        fst: do io=1,n_ions_cell
           icl=which_epe_ion(io)%new
           do i=1,n_tetrahedrons
              if(icl==tetra_atoms(i,1)) goto 1
           enddo
           cycle fst
1          ia_1=tetra_atoms(i,1)
           i1=epe(ia_1)%k
           scnd: do j=2,4
              ia_2=tetra_atoms(i,j)
              if(ia_2==0_i4_kind) cycle scnd
              j1=epe(ia_2)%k
              thrd: do k=j+1,5
                 ia_3=tetra_atoms(i,k)
                 if(ia_3==0_i4_kind) cycle thrd
                 k1=epe(ia_3)%k
                 rss1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))
                 rss2=dot_product(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
                 scal=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                      r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
                 rss1=sqrt(rss1)
                 rss2=sqrt(rss2)
                 cos_th=scal/(rss1*rss2)
                 theta=acos(cos_th)
                 sin_th=sin(theta)
                 g3b=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
                 e1=(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))/rss1
                 e2=(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))/rss2

                 m=epe(ia_2)%m
                 ia_2=which_epe_ion(m)%new
                 m=epe(ia_3)%m
                 ia_3=which_epe_ion(m)%new

                 DX(ia_3)=DX(ia_3)+g3b*(cos_th*e2(1)-e1(1))/(rss2*sin_th)
                 DY(ia_3)=DY(ia_3)+g3b*(cos_th*e2(2)-e1(2))/(rss2*sin_th)
                 DZ(ia_3)=DZ(ia_3)+g3b*(cos_th*e2(3)-e1(3))/(rss2*sin_th)
                 DX(ia_2)=DX(ia_2)+g3b*(cos_th*e1(1)-e2(1))/(rss1*sin_th)
                 DY(ia_2)=DY(ia_2)+g3b*(cos_th*e1(2)-e2(2))/(rss1*sin_th)
                 DZ(ia_2)=DZ(ia_2)+g3b*(cos_th*e1(3)-e2(3))/(rss1*sin_th)
                 DX(ia_1)=DX(ia_1)+g3b*((rss1-rss2*cos_th)*e1(1)+(rss2-rss1*cos_th)*e2(1))/ &
                      (rss2*rss1*sin_th)
                 DY(ia_1)=DY(ia_1)+g3b*((rss1-rss2*cos_th)*e1(2)+(rss2-rss1*cos_th)*e2(2))/ &
                      (rss2*rss1*sin_th)
                 DZ(ia_1)=DZ(ia_1)+g3b*((rss1-rss2*cos_th)*e1(3)+(rss2-rss1*cos_th)*e2(3))/ &
                      (rss2*rss1*sin_th)
              enddo thrd!k=j+1,5
           enddo scnd!j=2,4
        enddo fst!io=1,n_ions_cell
     endif

     est=zero
     do io=1,n_ions_cell 
        i=which_epe_ion(io)%new
        K=epe(I)%k
        DO J=reg_I_n_ions+1,reg_2a_n_ions       ! contrib due to reg. II discplacements
           if(epe(j)%m.eq.io) cycle
           RRR2=dot_product(epe(I)%r-epe(j)%r,epe(I)%r-epe(j)%r)
           IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle
           N=epe(J)%k
           RSSX=XS(I)-XS(J)
           RSSY=YS(I)-YS(J)
           RSSZ=ZS(I)-ZS(J)
           CSIX=XS(J)-epe(J)%r(1)
           CSIY=YS(J)-epe(J)%r(2)
           CSIZ=ZS(J)-epe(J)%r(3)
           RSS2=RSSX**2+RSSY**2+RSSZ**2
           RSSCSI=RSSX*CSIX+RSSY*CSIY+RSSZ*CSIZ
           RSS1=SQRT(RSS2)
           if(RRR2.GT.16.0_r8_kind) then
              ic=1
           else
              ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
              if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
           endif
           IF(RRR2.le.host%sr1(K,N,ic)**2) then
              RSSR2=1.0_r8_kind/RSS2
              RSSR8=(RSSR2**2)**2
              RSSD=RSSR2*host%D(K,N,ic)
              ESS1=(6.0_r8_kind *host%C(K,N,ic) + 8.0_r8_kind*RSSD)*RSSR8
              ESS2=(48.0_r8_kind*host%C(K,N,ic) +80.0_r8_kind*RSSD)*(RSSR8*RSSR2)
              IF(RRR2.le.host%sr2(K,N,ic)**2) then
                 RSSRO=1.0_r8_kind/(host%RO(K,N,ic)*RSS1)
                 BRO=host%B(K,N,ic)*RSSRO*EXP(-RSS1/host%RO(K,N,ic))
                 ESS1=-BRO+ESS1
                 ESS2= BRO*(RSSR2+RSSRO)-ESS2
              endif! RRR2.le.RADIUS_SECOND_SHORT_INTERACT
              ESS2=0.5_r8_kind*ESS2*RSSCSI
              OSO=0.5_r8_kind*ESS1
              OSS=ESS1+ESS2

              DX(I)=DX(I)+OSS*RSSX+OSO*CSIX
              DY(I)=DY(I)+OSS*RSSY+OSO*CSIY
              DZ(I)=DZ(I)+OSS*RSSZ+OSO*CSIZ
           endif! (RRR2.le.T

           IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle
           RSCX=XS(I)-XC(J)
           RSCY=YS(I)-YC(J)
           RSCZ=ZS(I)-ZC(J)

           RCSX=XC(I)-XS(J)
           RCSY=YC(I)-YS(J)
           RCSZ=ZC(I)-ZS(J)

           RCCX=XC(I)-XC(J)
           RCCY=YC(I)-YC(J)
           RCCZ=ZC(I)-ZC(J)

           RSRX=XS(I)-epe(J)%r(1)
           RSRY=YS(I)-epe(J)%r(2)
           RSRZ=ZS(I)-epe(J)%r(3)

           RCRX=XC(I)-epe(J)%r(1)
           RCRY=YC(I)-epe(J)%r(2)
           RCRZ=ZC(I)-epe(J)%r(3)

           CCIX=XC(J)-epe(J)%r(1)
           CCIY=YC(J)-epe(J)%r(2)
           CCIZ=ZC(J)-epe(J)%r(3)

           RSS2=RSSX**2+RSSY**2+RSSZ**2
           RSC2=RSCX**2+RSCY**2+RSCZ**2
           RCS2=RCSX**2+RCSY**2+RCSZ**2
           RCC2=RCCX**2+RCCY**2+RCCZ**2
           RSR2=RSRX**2+RSRY**2+RSRZ**2
           RCR2=RCRX**2+RCRY**2+RCRZ**2

           RCSCSI=RCSX*CSIX+RCSY*CSIY+RCSZ*CSIZ
           RSCCCI=RSCX*CCIX+RSCY*CCIY+RSCZ*CCIZ
           RCCCCI=RCCX*CCIX+RCCY*CCIY+RCCZ*CCIZ

           RSS3=RSS2*RSS1
           RCS3=RCS2*SQRT(RCS2)
           RSC3=RSC2*SQRT(RSC2)
           RCC3=RCC2*SQRT(RCC2)
           RSR1=SQRT(RSR2)
           RCR1=SQRT(RCR2)

           OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
           OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
           OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*RSR2))/RSR2
           OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
           OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
           OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*RCR2))/RCR2
           OSSCSI=-3.*OSS*RSSCSI/RSS2
           OCSCSI=-3.*OCS*RCSCSI/RCS2
           OSCCCI=-3.*OSC*RSCCCI/RSC2
           OCCCCI=-3.*OCC*RCCCCI/RCC2

           DISX=OSS*RSSX+OSC*RSCX+OSR*RSRX
           DISY=OSS*RSSY+OSC*RSCY+OSR*RSRY
           DISZ=OSS*RSSZ+OSC*RSCZ+OSR*RSRZ
           DICX=OCS*RCSX+OCC*RCCX+OCR*RCRX
           DICY=OCS*RCSY+OCC*RCCY+OCR*RCRY
           DICZ=OCS*RCSZ+OCC*RCCZ+OCR*RCRZ
           DIS3X=OSS*CSIX+OSC*CCIX+OSCCCI*RSCX+OSSCSI*RSSX
           DIS3Y=OSS*CSIY+OSC*CCIY+OSCCCI*RSCY+OSSCSI*RSSY
           DIS3Z=OSS*CSIZ+OSC*CCIZ+OSCCCI*RSCZ+OSSCSI*RSSZ
           DIC3X=OCS*CSIX+OCC*CCIX+OCSCSI*RCSX+OCCCCI*RCCX
           DIC3Y=OCS*CSIY+OCC*CCIY+OCSCSI*RCSY+OCCCCI*RCCY
           DIC3Z=OCS*CSIZ+OCC*CCIZ+OCSCSI*RCSZ+OCCCCI*RCCZ
           DSX(I)=DSX(I)+DISX+0.5*DIS3X
           DSY(I)=DSY(I)+DISY+0.5*DIS3Y
           DSZ(I)=DSZ(I)+DISZ+0.5*DIS3Z
           DCX(I)=DCX(I)+DICX+0.5*DIC3X
           DCY(I)=DCY(I)+DICY+0.5*DIC3Y
           DCZ(I)=DCZ(I)+DICZ+0.5*DIC3Z
        enddo!J=NN2A,NK2A       
! **done, contrib. due to region IIa 

! **treat contrib. from polarized region I
        DO J=1,reg_I_n_ions     
           if(epe(j)%m.eq.io.or.j.eq.i)  cycle
           N=epe(J)%k
           rrr=epe(i)%r-epe(j)%r

           RSSX=XS(I)-XS(J)
           RSSY=YS(I)-YS(J)
           RSSZ=ZS(I)-ZS(J)

           RSCX=XS(I)-XC(J)
           RSCY=YS(I)-YC(J)
           RSCZ=ZS(I)-ZC(J)

           RCSX=XC(I)-XS(J)
           RCSY=YC(I)-YS(J)
           RCSZ=ZC(I)-ZS(J)

           RCCX=XC(I)-XC(J)
           RCCY=YC(I)-YC(J)
           RCCZ=ZC(I)-ZC(J)

           RSRX=XS(I)-epe(J)%r(1)
           RSRY=YS(I)-epe(J)%r(2)
           RSRZ=ZS(I)-epe(J)%r(3)

           RCRX=XC(I)-epe(J)%r(1)
           RCRY=YC(I)-epe(J)%r(2)
           RCRZ=ZC(I)-epe(J)%r(3)

           RRSX=epe(I)%r(1)-XS(J)
           RRSY=epe(I)%r(2)-YS(J)
           RRSZ=epe(I)%r(3)-ZS(J)

           RRCX=epe(I)%r(1)-XC(J)
           RRCY=epe(I)%r(2)-YC(J)
           RRCZ=epe(I)%r(3)-ZC(J)
           RRR2=dot_product(rrr,rrr)
           RSS2=RSSX**2+RSSY**2+RSSZ**2
           RSC2=RSCX**2+RSCY**2+RSCZ**2
           RCS2=RCSX**2+RCSY**2+RCSZ**2
           RCC2=RCCX**2+RCCY**2+RCCZ**2
           RSR2=RSRX**2+RSRY**2+RSRZ**2
           RCR2=RCRX**2+RCRY**2+RCRZ**2
           RRS2=RRSX**2+RRSY**2+RRSZ**2
           RRC2=RRCX**2+RRCY**2+RRCZ**2

           RSS1=SQRT(RSS2)
           RSS3=RSS2*RSS1
           RCC3=RCC2*SQRT(RCC2)
           RSC3=RSC2*SQRT(RSC2)
           RCS3=RCS2*SQRT(RCS2)
           RSR1=SQRT(RSR2)
           RCR1=SQRT(RCR2)
           RRS1=SQRT(RRS2)
           RRC1=SQRT(RRC2)

! **short region interacting part of area 1 - start
           if(RRR2.GT.16.0_r8_kind) then
              ic=1
           else
              ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
              if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
           endif
           IF(RRR2.le.host%sr1(K,N,ic)**2) then
              RSSR2=1.0_r8_kind/RSS2
              ESS1=(6.0_r8_kind*host%C(K,N,ic)+8.0_r8_kind*host%D(K,N,ic)*RSSR2)*((RSSR2**2)**2)
              IF(RRR2.LE.host%sr2(K,N,ic)**2)  &
                   ESS1=ESS1-host%B(K,N,ic)/(host%RO(K,N,ic)*RSS1)*EXP(-RSS1/host%RO(K,N,ic))
              OSX=ESS1*RSSX
              OSY=ESS1*RSSY
              OSZ=ESS1*RSSZ
              DX(I)=DX(I)+OSX
              DY(I)=DY(I)+OSY
              DZ(I)=DZ(I)+OSZ
           endif! RRR2.le.T
! **short region interacting part of area 1 - end

           OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
           OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
           OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
           OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
           IF(ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST) then  
              OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*RSR2))/RSR2
              OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*RCR2))/RCR2
              ORS=(Q_ion(epe(i)%k)*Q_SHELL(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRS1)/RRS1+ERFO*EXP(-ET2*RRS2))/RRS2
              ORC=(Q_ion(epe(i)%k)*Q_NUCLEAR(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRC1)/RRC1+ERFO*EXP(-ET2*RRC2))/RRC2
           else
              OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))/(RSR2*RSR1)
              OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))/(RCR2*RCR1)
              ORS=(Q_ion(epe(i)%k)*Q_SHELL(N))/(RRS2*RRS1)
              ORC=(Q_ion(epe(i)%k)*Q_NUCLEAR(N))/(RRC2*RRC1)
           endif!ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST/else

           OSSX=OSS*RSSX
           OSSY=OSS*RSSY
           OSSZ=OSS*RSSZ
           OSCX=OSC*RSCX
           OSCY=OSC*RSCY
           OSCZ=OSC*RSCZ
           OCSX=OCS*RCSX
           OCSY=OCS*RCSY
           OCSZ=OCS*RCSZ
           OCCX=OCC*RCCX
           OCCY=OCC*RCCY
           OCCZ=OCC*RCCZ
           DSX(I)=DSX(I)+OSSX+OSCX+OSR*RSRX
           DSY(I)=DSY(I)+OSSY+OSCY+OSR*RSRY
           DSZ(I)=DSZ(I)+OSSZ+OSCZ+OSR*RSRZ
           DCX(I)=DCX(I)+OCSX+OCCX+OCR*RCRX
           DCY(I)=DCY(I)+OCSY+OCCY+OCR*RCRY
           DCZ(I)=DCZ(I)+OCSZ+OCCZ+OCR*RCRZ
        enddo!J=1,NNA
! **done polirized region I

        CSIX=XS(I)-epe(I)%r(1)
        CSIY=YS(I)-epe(I)%r(2)
        CSIZ=ZS(I)-epe(I)%r(3)
        CCIX=XC(I)-epe(I)%r(1)
        CCIY=YC(I)-epe(I)%r(2)
        CCIZ=ZC(I)-epe(I)%r(3)
        CMCSX=XC(I)-XS(I)
        CMCSY=YC(I)-YS(I)
        CMCSZ=ZC(I)-ZS(I)
        CSIS2=CSIX**2+CSIY**2+CSIZ**2
        CSIC2=CCIX**2+CCIY**2+CCIZ**2
        ET2S=ET2*CSIS2
        ET2C=ET2*CSIC2

        IF(ET2S.GT.F1) then
           CSIS=SQRT(CSIS2)
           OEDS=(ERF(ERROR_FUNCTION_PARAMETER*CSIS)/CSIS+ERFO*EXP(-ET2S))/CSIS2
        else
           OEDS=ERF1*(1.0_r8_kind-ET2S*(0.6_r8_kind-0.2142856_r8_kind*ET2S))
        endif! ET2S.GT.F1)/else
    
        IF(ET2C.GT.F1) then
           CSIC=SQRT(CSIC2)
           OEDC=(ERF(ERROR_FUNCTION_PARAMETER*CSIC)/CSIC+ERFO*EXP(-ET2C))/CSIC2
        else
           OEDC=ERF1*(1.0_r8_kind-ET2C*(0.6_r8_kind-0.2142856_r8_kind*ET2C))
        endif! ET2C.GT.F1

        OEDS=Q_SHELL(K)*Q_ion(epe(i)%k)*OEDS
        OEDC=Q_NUCLEAR(K)*Q_ion(epe(i)%k)*OEDC
        
        gs(:)=zero
        gc(:)=zero
        psi=zero
        DO IG=1,n_bs_points
           gig=dot_product(gstr(ig,:),gstr(ig,:))
           res=exp(-pta2*gig)/gig*qpivc
           GIS2=PI2*(GSTR(IG,1)*XS(I)+GSTR(IG,2)*YS(I)+GSTR(IG,3)*ZS(I))
           GIC2=PI2*(GSTR(IG,1)*XC(I)+GSTR(IG,2)*YC(I)+GSTR(IG,3)*ZC(I))
           PSIS=RSIN(IG)*COS(GIS2)-RCOS(IG)*SIN(GIS2)
           PSIC=RSIN(IG)*COS(GIC2)-RCOS(IG)*SIN(GIC2)
           gs(:)=gs(:)+psis*gstr(ig,:)
           gc(:)=gc(:)+psic*gstr(ig,:)
           psi=psi+RSIN(IG)*SIN(GIC2)+RCOS(IG)*COS(GIC2)
           do jo=1,n_ions_cell  ! treat region one
              j=which_epe_ion(jo)%new
              gcel=dot_product(gstr(ig,:),R_ION_IN_CELL(jo,:))
              gcelpi=pi2*gcel
              qns=q_z(jo)*sin(gcelpi)
              qnc=q_z(jo)*cos(gcelpi)
              prsin=res*qns
              prcos=res*qnc
              gpsis(:)=gstr(ig,:)*(prcos*sin(gis2)-prsin*cos(gis2))
              gpsic(:)=gstr(ig,:)*(prcos*sin(gic2)-prsin*cos(gic2))
              DSX(j)=DSX(j)+gpsis(1)*PI*q_shell(k)
              DSY(j)=DSY(j)+gpsis(2)*PI*q_shell(k)
              DSZ(j)=DSZ(j)+gpsis(3)*PI*q_shell(k)
              DCX(j)=DCX(j)+gpsic(1)*PI*q_nuclear(k)
              DCY(j)=DCY(j)+gpsic(2)*PI*q_nuclear(k)
              DCZ(j)=DCZ(j)+gpsic(3)*PI*q_nuclear(k)
           enddo! jo=1,n_ions_cell
        enddo!IG=1,n_bs_points

        PI2S=PI*Q_SHELL(K)      ! divided by two for gradients as for energy
        PI2C=PI*Q_NUCLEAR(K)
        est=est+psi*Q_NUCLEAR(K)/two
        OEPKX=PK(K)*CMCSX
        OEPKY=PK(K)*CMCSY
        OEPKZ=PK(K)*CMCSZ

        DSX(I)=DSX(I)+CSIX*OEDS+GS(1)*PI2S-OEPKX
        DSY(I)=DSY(I)+CSIY*OEDS+GS(2)*PI2S-OEPKY
        DSZ(I)=DSZ(I)+CSIZ*OEDS+GS(3)*PI2S-OEPKZ
        DCX(I)=DCX(I)+OEPKX+CCIX*OEDC+GC(1)*PI2C
        DCY(I)=DCY(I)+OEPKY+CCIY*OEDC+GC(2)*PI2C
        DCZ(I)=DCZ(I)+OEPKZ+CCIZ*OEDC+GC(3)*PI2C
     enddo!io=1,n_ions_cell

   END SUBROUTINE lattice_grad_gopt

end module main_epe_module

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
subroutine main_epe_block()
  !-------------------------------------------------------------------
  !
  !  Purpose: subroutine starts EPE calculations:
  !           - call epe_draver
  !           - send epe_data to slaves
  !           - make epe- optimization loop (main_epe)
  !
  !  Subroutine called by: main_master
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  !------------ Modules used -----------------------------------------
#include "def.h"
  use type_module ! type specification parameters
  use filename_module, only: data_dir
  use time_module
  use timer_module
  use iounitadmin_module
  use unique_atom_module
  use orbitalprojection_module
  use options_module,   only: options_xcmode, xcmode_model_density
  use fit_coeff_module, only: fit_coeff_send, ICHFIT, IXCFIT
  use symmetry_data_module, only: get_totalsymmetric_irrep
  use comm_module
  use epe_module
  use main_epe_module, only : main_epe,                &
                              finish_epe,              &
                              n_epe_vacancies,         &
                              epe_send_init_to_slave,  &
                              epe_send_finish_to_slave

  use epecom_module, only: reg_I_pg, reg_I_n_ions, etot_epe, &
       epe_rel_converged, n_pgepe_iterations, epg_cluster_reg_I, &
       get_epe_energies, &
                            qau_qepe,rel_converged_unit, &
                            make_epe_reference,                           &
                            dealloc_epe_ref,                              &
                            eperef_unit,                                  &
                            epe_iter,                                     &
                            n_grads_master,                               &
                            n_grads_slave,                                &
                            n_grads_total,                                &
                            var_epe_contrib_atstart,                      &
                            embed_convergence_check,                      &
                            embed_convergence_limit,                      &
                            dg_convergence_reached, &
                            reg_2a_treated, reg_2a_n_ions,                &
                            end_treated_region,                           &
                            n_ls,                                         &
                            reset,                                        &
                            use_lin_search,                               &
                            basic_action,                                 &
                            etot_epe_0

  use energy_calc_module, only: get_energy
  implicit none

  !== Interrupt end of public interface of module ====================

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of subroutines ------------------------
  !------------ Declaration of local constants --------------------
  !------------ Declaration of local variables --------------------
  real(kind=r8_kind)    :: var_contrib
  real(kind=r8_kind)    :: energy,e_sh_regIau
  integer(kind=i4_kind) :: n_proc
  logical               :: model_density, conv_achived
  integer(kind=i4_kind) :: IFIT
  !-------------------------------------------------------------------

  !------------ Executable code -----------------------------------
  call write_to_trace_unit("main_epe_block: start of epe relaxation")
    call start_timer(timer_epe)
    call write_to_output_units(" main_master => epe_driver")
    call epe_driver
    call write_to_output_units(" main_master <= epe_driver")

! number of gradients to be calculated by each processor
! master is also involved

    n_proc=comm_get_n_processors()

    if(reg_2a_treated) then
       n_grads_total  = (reg_2a_n_ions-n_epe_vacancies)
       end_treated_region = reg_2a_n_ions
    else
       n_grads_total  = (reg_I_n_ions-n_epe_vacancies)
       end_treated_region = reg_I_n_ions
    end if

    n_grads_master = n_grads_total - &
                     (n_grads_total/n_proc)*(n_proc-1)

    if(n_proc > 1) then
    n_grads_slave = (n_grads_total-n_grads_master)/(n_proc-1)
    else
    n_grads_slave=0
    end if
!!$    call write_to_output_units('n_proc', n_proc)
!!$    call write_to_output_units('n_grads_total',  n_grads_total )
!!$    call write_to_output_units('n_grads_slave',  n_grads_slave )
!!$    call write_to_output_units('n_grads_master', n_grads_master)

    if(comm_parallel()) then
       model_density=options_xcmode()==xcmode_model_density
       IFIT = ICHFIT
       if(model_density) IFIT = IFIT + IXCFIT
       call fit_coeff_send(IFIT)
       call epe_send_data()
    end if

    dg_convergence_reached = .false.
    n_ls=0; reset=.false.;use_lin_search=.true.
    epeit: do  epe_iter=1,n_pgepe_iterations
       call start_timer(timer_epe_cycle)

       dealloc_epe_ref=epe_iter.eq.n_pgepe_iterations.or.dg_convergence_reached
!!$           call write_to_output_units("main_epe_block: epe_read_write_reference")
       call epe_read_write_reference()
!!$           call write_to_output_units("epe_read_write_reference done")

       if(comm_parallel() )  then
          call start_timer(timer_send_epe_ref)
          call epe_send_reference()
          call stop_timer(timer_send_epe_ref)
       end if

       call start_timer(timer_epe_forces)
!!$        call write_to_output_units("main_epe_block: epe_field_and_forces_par")
       call epe_field_and_forces_par()

       call stop_timer(timer_epe_forces)

       if(comm_parallel() ) then
          call start_timer(timer_collect_grads)
          call epe_collect_gradients()
          call stop_timer(timer_collect_grads)
       end if

       if(dealloc_epe_ref) then
          call stop_timer(timer_epe_cycle)
          exit
       endif

       if(epe_iter.eq.1) epg_cluster_reg_I_at_start=epg_cluster_reg_I

       call start_timer(timer_main_epe)
       call main_epe()
       if(use_lin_search.and.basic_action/=0) then
          if(etot_epe_0 <= etot_epe .and. n_ls < 2) then
             call lin_search_epe()
          else
             n_ls=0
             reset=.false.
             etot_epe_0=etot_epe
          end if
       else
          reset=.false.
          etot_epe_0=etot_epe
       end if
       call stop_timer(timer_main_epe)

       call stop_timer(timer_epe_cycle)
    enddo epeit

        ! check  for convergence
        call get_epe_energies(var_epe_contrib=var_contrib)
        epe_rel_converged= &
             abs(var_contrib-var_epe_contrib_atstart).lt.embed_convergence_limit &
             .or. .not.embed_convergence_check .or. make_epe_reference

        ! etot_epe - energy minimized with EPE optimization
        print*,etot_epe,' etot_epe eV, converged',epe_rel_converged
        write(output_unit,*) ' etot_epe eV, converged',epe_rel_converged

!---------------------------------------------------------------------------------------
        call get_epe_energies(diff_ec_ecref=ec_ecref, &
             etot_epe_corrected=etot_epe_corr)
        etot_epe_corr=etot_epe_corr+epg_cluster_reg_I
                                   ! fitted density + nucleaR
        ! EPE energy with some contributions deleted to evoid double
        ! counting as they are in PG total energy
!---------------------------------------------------------------------------------------

        print*,'pg_energies (ref,var)', &
             epg_reference_reg_I,epg_cluster_reg_I,etot_epe_corr

        write(output_unit,*) ' energies of interaction of the cluster and epe  region I '
        write(output_unit,*) 'calculated with use of fitted charge density for reference and actual configurations'
        write(output_unit,*) epg_reference_reg_I, epg_cluster_reg_I
        write(output_unit,*)

        print*,'epeg_nuc (var,rs,rc)',epg_nuc
        write(output_unit,*) ' energy of interaction of the claster nuclei with region I epe centers '
        write(output_unit,*) epg_nuc
        write(output_unit,*)
        print*,dot_product(reg_I_pg(1+n_epe_vacancies)%rs,reg_I_pg(1+n_epe_vacancies)%rs)
        print*,dot_product(reg_I_pg(1+n_epe_vacancies)%rc,reg_I_pg(1+n_epe_vacancies)%rc)
        print*,'finish_epe done'

    call finish_epe(.true.)
    if(comm_parallel() ) call epe_send_finish_to_slave()
    call returnclose_iounit(eperef_unit)



    if(.not.epe_rel_converged)  then
        call get_epe_energies(eshort_reg_Iau=e_sh_regIau)
        print*,'Relaxation is not converged',etot_epe_corr,e_sh_regIau
        rel_converged_unit=&
        openget_iounit(file=trim(data_dir)//"/epe_rel_unconverged",&
                                form='unformatted',status='unknown')
        write(rel_converged_unit) etot_epe_corr,epg_cluster_reg_I
        inquire(file=trim(data_dir)//"/epe_rel_converged",exist=conv_achived)
      if(conv_achived)  then
        print*,'file epe_rel_converged exist'
        call system(" rm "//trim(data_dir)//"/epe_rel_converged")
        print*," rm "//trim(data_dir)//"/epe_rel_converged"
      else
        print*,'file epe_rel_converged do not exist'
      endif
        call returnclose_iounit(rel_converged_unit)

    else ! epe_rel_converged

        inquire(file=trim(data_dir)//"/epe_rel_unconverged",exist=conv_achived)
        if(conv_achived) &
             call system(" rm "//trim(data_dir)//"/epe_rel_unconverged")
        rel_converged_unit=&
             openget_iounit(file=trim(data_dir)//"/epe_rel_converged",&
             form='unformatted',status='unknown')
        write(rel_converged_unit) etot_epe_corr,epg_cluster_reg_I
        call returnclose_iounit(rel_converged_unit)
        call get_epe_energies(eshort_reg_Iau=e_sh_regIau)
        print*,'epe relaxation is claimed to be converded', etot_epe_corr,e_sh_regIau
    endif


        call get_energy(tot=energy)
        print*,'PG QM totel energy', energy

    print*,'main_epe_block: lattice_energy', etot_epe_corr

    call write_to_trace_unit("main_epe_block: end of epe relaxation")
    call stop_timer(timer_epe)
end subroutine main_epe_block

!             Added from new epe_field_and_forces.f90

subroutine gxepe_allocate

   use type_module
   use pointcharge_module
   use unique_atom_module
   use ewaldpc_module
   use iounitadmin_module

   integer(kind=i4_kind)::i_ua,n_equal_atoms
   write(output_unit,*) &
        'unique_atom_unique_read: N_unique_atoms,n_timps', &
        N_unique_atoms,n_timps

       if(associated(gxepe_array)) then
      do i_ua=1,N_unique_atoms
         deallocate(gxepe_array(i_ua)%position,stat=ewa_allocstat(23))
         ASSERT(ewa_allocstat(23).eq.0)
      end do

      deallocate(gxepe_array,gxepe_impu,stat=ewa_allocstat(22))
                 ASSERT(ewa_allocstat(22).eq.0)
   end if


   allocate(gxepe_array(N_unique_atoms+n_timps), &
        gxepe_impu(N_unique_atoms+n_timps),stat=ewa_allocstat(22))
       if(ewa_allocstat(22).ne.0) call error_handler('allocation of gxepe_array failed')
          ewa_allocstat(22)=1
   do i_ua=1,N_unique_atoms+n_timps
      if(i_ua.gt.N_unique_atoms) then
         n_equal_atoms=unique_timps(i_ua-N_unique_atoms)%n_equal_atoms
      else
         n_equal_atoms=unique_atoms(i_ua)%n_equal_atoms
      end if

      allocate(gxepe_array(i_ua)%position(3,n_equal_atoms),stat=ewa_allocstat(23))
      if(ewa_allocstat(23).ne.0) &
           call error_handler('allocation of gxepe_array(i_ua)%position failed')
         ewa_allocstat(23)=1
           gxepe_array(i_ua)%position(:,:)=0.0_r8_kind
   enddo ! i_ua=1,N_unique_atoms+

 end subroutine gxepe_allocate

 subroutine read_gxtimps(io_u,io_gxepe)
  use type_module
  use pointcharge_module
  use unique_atom_module
  use ewaldpc_module
  use operations_module
  use filename_module
  integer(kind=i4_kind), intent(in)::io_u,io_gxepe
  integer             ::  status,i_ua,counter_equal
  integer             ::  ieq_dummy,indexes(7),impu
  real(kind=r8_kind) :: x_coord,y_coord,z_coord,z_dummy
  real(kind=r8_kind), dimension(3) :: r_gxepe

  do i_ua=1,N_unique_atoms
     counter_equal=1

     do
        if(counter_equal>unique_atoms(i_ua)%n_equal_atoms) exit
        read(io_u,*, iostat=status) z_dummy,x_coord,y_coord,z_coord,ieq_dummy, &
             indexes,impu
        if(ex_gxepe.and.impu.ne.0) read(io_gxepe,*) r_gxepe
        if (status .gt. 0) call error_handler &
             ("unique_atom_unique_read: "//trim(data_dir)//'/gxfile')
        if(ieq_dummy.ne.0) then
           if(ex_gxepe) then
              gxepe_impu(i_ua)=impu
              if(impu.ne.0) gxepe_array(i_ua)%position(:,counter_equal)=r_gxepe
           endif
           counter_equal=counter_equal+1
        endif
     enddo
  enddo ! i_ua=1,N_unique_atoms

  do i_ua=1,n_timps
     counter_equal=1
     do
        if(counter_equal>unique_timps(i_ua)%n_equal_atoms) exit
        if(operations_gx_epeformat) then
          read(io_u,*) z_dummy,x_coord,y_coord,z_coord,ieq_dummy, &
               indexes,impu
          if(ex_gxepe.and.impu.ne.0) read(io_gxepe,*) r_gxepe
          else
        if(operations_gx_highprec) then
           read(io_u,'(F5.2,3(2x,F21.12),i3)', iostat=status) &
                z_dummy,x_coord,y_coord,z_coord,ieq_dummy
        else
           !               read(io_u,'(F5.2,3(2x,F13.7),i3)', iostat=status) &
           !                        z_dummy,x_coord,y_coord,z_coord,ieq_dummy
           read(io_u,*)  z_dummy,x_coord,y_coord,z_coord,ieq_dummy,indexes,impu
           if(ex_gxepe.and.impu.ne.0) read(io_gxepe,*) r_gxepe
        endif
     end if


        if (status .gt. 0) call error_handler &
             ("unique_atom_unique_read: "//trim(data_dir)//'/gxfile')
        if(ieq_dummy.ne.0) then
           if(counter_equal==1) then
              unique_timps(i_ua)%position_first_ea(1)=x_coord
              unique_timps(i_ua)%position_first_ea(2)=y_coord
              unique_timps(i_ua)%position_first_ea(3)=z_coord
           endif
           if(ex_gxepe) then
              gxepe_impu(i_ua+N_unique_atoms)=impu
              if(impu.ne.0) gxepe_array &
                   (N_unique_atoms+i_ua)%position(:,counter_equal)=r_gxepe
           endif
           counter_equal=counter_equal+1
        endif
     end do

  end do
    end subroutine read_gxtimps

   subroutine unique_atom_make_epegx()
     ! Purpose: write a template for EPE related gx file
     !          the gx file is used as a input for the geometry optimizer
     ! Subroutine called by: main_gradient ?
     !** End of interface *****************************************
     use type_module
     use pointcharge_module
     use unique_atom_module
     use iounitadmin_module
     use filename_module
     use operations_module, only: operations_gx_highprec
     ! ----------- declaration of local variables -------------
     integer(kind=i4_kind) :: i_ua,i_ea,io_u,counter, &
        deck(7)=(/0,0,0,0,0,0,-1/)
     external error_handler
     !
     ! ----------- executable code -----------------------------
     !
!     deck=0
     io_u=get_iounit()
     open (io_u,status='unknown',form='formatted',file=&
          trim(outfile('gxfile')))
     counter=1
     do i_ua=1,n_unique_atoms
        do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
           if (operations_gx_highprec) then
              write(io_u,'(F5.2,3(2x,F21.12),2I3,2X,3I3,2X,3I3,i5)') &
                   unique_atoms(i_ua)%z,unique_atoms(i_ua)%position(:,i_ea),&
                   i_ua,counter,deck
           else
              write(io_u,'(F5.2,3(2x,F13.7),2I4,2X,3I3,2X,3I3,I5)') &
                   unique_atoms(i_ua)%z,unique_atoms(i_ua)%position(:,i_ea),&
                   i_ua,counter,deck
           endif
           counter=counter+1
        end do
     end do
     do i_ua=1,n_timps
        do i_ea=1,unique_timps(i_ua)%n_equal_atoms
           if (operations_gx_highprec) then
              write(io_u,'(F5.2,3(2x,F21.12),2I3,2X,3I3,2X,3I3,i5)') &
                   unique_timps(i_ua)%z,unique_timps(i_ua)%position(:,i_ea),&
                   i_ua,counter,deck
           else
              write(io_u,'(F5.2,3(2x,F13.7),2I4,2X,3I3,2X,3I3,I5)') &
                   unique_timps(i_ua)%z,unique_timps(i_ua)%position(:,i_ea),&
                   i_ua,counter,deck
           end if
        end do
     end do

        if (operations_gx_highprec) then
!       write(io_u,'(F5.1)') -7.0_r8_kind
        write(io_u,'(F5.2,3(2x,F13.7),2I4,2X,3I3,2X,3I3,I5)')  &
        -7.0_r8_kind,0.0,0.0,0.0, 0,0, deck
        else
     write(io_u,'(F5.1,3(2x,F13.7),2I4,2X,3I3,2X,3I3,I5)') &
     -7.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0,0,deck(1:6),0
        endif
     close(io_u)
     call return_iounit(io_u)
   end subroutine unique_atom_make_epegx








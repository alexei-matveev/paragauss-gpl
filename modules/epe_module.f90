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
!===============================================================
! Public interface of module
!===============================================================
module epe_module
  !---------------------------------------------------------------
  !       
  !  Purpose: ......
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
  !------------ Modules used --------------------------------------
#include "def.h"
  use type_module ! type specification parameters
  use unique_atom_module
  use gamma_module
  use datatype
  use msgtag_module
  use solid_harmonics_module, only : solid_harmonics_calc
  use solhrules_module, only: solhrules_differential
  use fitcontract_module, only: fitcontract,fitcontract_eperef
  use integralpar_module
  use iounitadmin_module
  use filename_module
  use gradient_data_module
  use integralpar_module
  use options_module, only: options_integral_expmax, options_split_gradients, &
                            options_xcmode, xcmode_model_density, &
                            options_spin_restricted
  use int_data_2cob3c_module, only:                              &
                            ua1,ua2,&    ! indices of unique atoms of quadrupel
                            ua1_basis,ua2_basis, & ! => uai%l_ob(0)
                            quadrupel,&
                            prim_int_3cob_grad,&
                            prim_int_3cob_coul_grad,&
                            prim_int_2cob_ks_grad
  use epecom_module, only:  reg_I_pg,reg_I_n_ions,output_epe,i_ir_eperef, &
                            e_epe_atstart,e_epe_final,etot_epe, &
                            epe_rel_converged,use_epe_reference,&
                            n_pgepe_iterations, ex_pgdata, &
                            ecoul_epecluster,eshort_epecluster, &
                            ecoul_vaccluster, &
                            epe, q_shell,q_nuclear,epg_cluster_reg_I, &
                            ndt, n_gen_ions,&
                            get_epe_energies, &
                            qau_qepe,rel_converged_unit,diffpg_ec_ecref, &
                            make_epe_reference, &
                            dealloc_epe_ref,        &
                            eperef_unit,            &
                            epe_iter,               &
                            n_grads_master,         &
                            n_grads_slave,          &
                            n_grads_total,          &
                            reg_2a_treated,         &
                            reg_2a_n_ions
  use comm_module
  use main_epe_module, only: main_epe,n_epe_vacancies,finish_epe
  use orbitalprojection_module
  use fit_coeff_module, only: coeff_charge, coeff_charge_eperef                &
                            , fit_coeff_shutdown
  use symmetry_data_module, only: get_totalsymmetric_irrep
  use energy_calc_module, only: get_energy
  use pointcharge_module
  use msgtag_module
!

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  intrinsic max,maxval
   type(unique_atom_type), public, pointer  :: unique_timps_eperef(:)
   type(unique_atom_type), pointer  :: ua
   type(unique_atom_type), public, pointer  :: ua_eperef 
   real(kind=r8_kind), allocatable, dimension(:,:) :: epe_send_array
   real(kind=r8_kind), allocatable, dimension(:)   :: epe_convergence_energies
   real(kind=r8_kind), public    :: pg_energy,ec_ecref,etot_epe_corr, &
                                    epg_reference_reg_I,epg_nuc,      &
                                    epg_cluster_reg_I_at_start
   real(kind=r8_kind)    :: zcc ! core charge of PP atom
   integer(kind=i4_kind) :: status, length
   integer(kind=i4_kind) :: N_partner         
   integer(kind=i4_kind) :: n_irreps_eperef,lmax_sym_eperef
   integer(kind=i4_kind) :: n_unique_atoms_eperef!!!,n_pgepe_iterations
   integer(kind=i4_kind) :: n_timps_eperef
   integer(kind=i4_kind) :: independent_fct, n_coeff_charge_eperef 
   integer(kind=i4_kind), public :: epemalloc_stat(51)=0
   integer(kind=i4_kind) :: i,ua3,max_order,                          &
                            max_gamma,lmax_ch,lmax_ch_eperef,k,                  & 
                            n_equal_c,lm,n_indep_max,n_indep,i_l,i_sum,          &
                            n_indep_fcts,n_contributing_fcts,i_ind,i_cont,i_ea3
   integer(kind=i4_kind), allocatable,public, dimension(:,:) :: N_partner_ua
!!$   integer(kind=i4_kind), allocatable, public,dimension(:,:) :: n_c_expua
   integer(kind=i4_kind), allocatable, public, dimension(:)   :: n_irreps_eperef_ua
   integer(kind=i4_kind), allocatable, public, dimension(:)   :: lmax_sym_eperef_ua
   integer(kind=i4_kind), allocatable, dimension(:),public   :: i_ir_eperef_ua
   logical,public  :: get_epe_reference 
   logical :: do_rotation
   integer(kind=i4_kind) :: ind

  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public ::  epe_send_data,                 & 
             epe_receive_data,              &
             epe_read_write_reference,      &
             epe_send_reference,            &
             epe_receive_reference,         &
             epe_field_and_forces_par,      &
             epe_collect_gradients
  public :: print_epemod_alloc

  !------------ Declaration of constants and variables ----
  !================================================================
  ! End of public interface of module
  !================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine epe_send_data()
    !  Purpose: sends epe_data(coordinates of epe-centers) 
    !           to all slaves. Called by master
    !           
    !------------ Modules used ------------------- ---------------    
    implicit none

   !------------ Declaration of local variables -----------------
   integer             :: status, index
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
 
  call comm_init_send(comm_all_other_hosts, msgtag_epe_data_sent)
  
! pack epe_data to be sent
  call commpack(reg_I_n_ions,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at reg_I_n_ion")
  call commpack(reg_2a_n_ions,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at reg_2a_n_ion")  
  call commpack(output_epe,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at output_epe")
  call commpack(i_ir_eperef,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at i_ir_eperef")  
  call commpack(e_epe_atstart,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at e_epe_atstart")  
  call commpack(e_epe_final,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at ,e_epe_final")  
  call commpack(etot_epe,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at etot_epe")  
  call commpack(epe_rel_converged,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at epe_rel_converged")  
  call commpack(use_epe_reference,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at use_epe_reference")  
  call commpack(n_pgepe_iterations,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_pgepe_iterations")  
  call commpack(n_epe_vacancies,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_epe_vacancies")
  call commpack(make_epe_reference,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at make_epe_reference")
! call commpack(n_gen_ions,status)
! if(status.ne.0) call error_handler &
!                 ("epe_send_data : error at n_gener_ions")
! call commpack(q_nuclear(1),ndt,1,status)
! if(status.ne.0) call error_handler &
!                 ("epe_send_data : error at q_nuclear")
! call commpack(q_shell(1),ndt,1,status)
! if(status.ne.0) call error_handler &
!                 ("epe_send_data : error at q_shell")
  do index=1,n_gen_ions
  call commpack(epe(index)%k,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at epe")
  end do
  call commpack(reg_2a_treated,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at reg_2a_treated")
  call commpack(n_grads_master,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_grads_master")  
  call commpack(n_grads_slave,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_grads_slave")  
  call comm_send()
!!$  print*, 'Master: epe_data sent'
!!$  call write_to_output_units( 'Master: epe_data sent' )
  end subroutine epe_send_data
  !*************************************************************

  !*************************************************************
  subroutine epe_receive_data()
    !  Purpose: Receives epe_data(coordinates of epe_centers)
    !           from master. Called by slaves.
    !           gamma-function is also initialized
    !------------ Modules used ------------------- ---------------
    implicit none
   !------------ Declaration of local variables -----------------
   integer             :: status, index
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
  call write_to_output_units("Slave : Receiving epe_data")
!!$print*, 'Slave: receiving epe_data'
! unpack epe_data received
  call communpack(reg_I_n_ions,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at reg_I_n_ion")
  call communpack(reg_2a_n_ions,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at reg_2a_n_ion")
  call communpack(output_epe,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at output_epe")
  call communpack(i_ir_eperef,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at i_ir_eperef")  
  call communpack(e_epe_atstart,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at e_epe_atstart")  
  call communpack(e_epe_final,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at ,e_epe_final")  
  call communpack(etot_epe,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at etot_epe")  
  call communpack(epe_rel_converged,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at epe_rel_converged")  
  call communpack(use_epe_reference,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at use_epe_reference")  
  call communpack(n_pgepe_iterations,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at n_pgepe_iterations")  
  call communpack(n_epe_vacancies,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at n_epe_vacancies")
  call communpack(make_epe_reference,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at make_epe_reference")
!-----------------------------------------------------------------
! call communpack(n_gen_ions,status)
! if(status.ne.0) call error_handler &
!                 ("epe_receive_data : error at n_gener_ions")
! call communpack(q_nuclear(1),ndt,1,status)
! if(status.ne.0) call error_handler &
!                 ("epe_receive_data : error at q_nuclear")
! call communpack(q_shell(1),ndt,1,status)
! if(status.ne.0) call error_handler &
!                 ("epe_receive_data : error at q_shell")
! if( allocated(epe) ) deallocate( epe )
! allocate(epe(n_gen_ions),stat=status)
! if(status.ne.0) call error_handler &
!                 ("epe_receive_data : error at epe allocation")
!---------------------------------------------------------------------
  do index=1,n_gen_ions  
  call communpack(epe(index)%k,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at epe")
  end do 
  
 call communpack(reg_2a_treated,status)
  if(status.ne.0) call error_handler &
                  ("epe_receive_data : error at ml_displacements")
 
  if(reg_2a_treated ) then
     allocate(reg_I_pg(reg_2a_n_ions),stat=epemalloc_stat(1) )
  else
     allocate( reg_I_pg(reg_I_n_ions), stat=epemalloc_stat(1) )
  end if
  ASSERT(epemalloc_stat(1).eq.0)
  epemalloc_stat(1)=1

  call communpack(n_grads_master,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_grads_master")  
  call communpack(n_grads_slave,status)
  if(status.ne.0) call error_handler &
                  ("epe_send_data : error at n_grads_slave")  

!!$print*, 'Slave: epe_data_received, n_grads_slave=', n_grads_slave
!!$  call write_to_output_units("Slave : Receiving epe_data - OK")
  end subroutine epe_receive_data

  subroutine epe_read_write_reference()
    !  Purpose: Writes and reads epe_reference configuration.
    !           Called by master.
    !------------ Modules used ------------------- ---------------
    use epecom_module, only: epedata_dir
    implicit none
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind) :: status
   integer(kind=i4_kind) :: lmax_ch
   !------------ Declaration of subroutines used ----------------
   external error_handler
   intrinsic max, maxval
   !------------ Executable code --------------------------------

    inquire(file=trim(inpfile("pgepe_reference")), exist=get_epe_reference)
    make_epe_reference=make_epe_reference.and..not.get_epe_reference 
!!$print*, 'epe_read_write_reference: MAKE_EPE_REFERENCE :', make_epe_reference
!!$    if(get_epe_reference) &
!!$             print*,'________epe_field_and_forces: ** using epe reference **' 
!!$    if( make_epe_reference) print*, & 
!!$             'epe_field_and_forces: ** making epe reference **'
 if(make_epe_reference) then    !   Store epe_reference configuration
   eperef_unit=openget_iounit(file=trim(inpfile("pgepe_reference")), &
                                  form='unformatted',status='unknown')
 
         write(eperef_unit)  n_unique_atoms, &
                             n_timps, &
                             size(coeff_charge)      !(1)
         write(eperef_unit)  unique_atoms(:)%lmax_ch,            &
                             unique_atoms(:)%z,                  &
                             unique_atoms(:)%zc, &
                             unique_atoms(:)%N_equal_atoms,      &  !(2)
                             unique_atoms(:)%N_glob_cons_ch,     &
                             coeff_charge(:)
         if(n_timps.ne.0) write(eperef_unit) unique_timps(:)%z,  &
                                             unique_timps(:)%zc, &
                                             unique_timps(:)%N_equal_atoms
         lmax_ch=maxval(unique_atoms(:)%lmax_ch)
         write(eperef_unit) orbitalprojection_ch
   unique_3: do ua3 = 1,N_unique_atoms+n_timps
     if(ua3.gt.N_unique_atoms) then
        ua=>unique_timps(ua3-N_unique_atoms)
     else
        ua=>unique_atoms(ua3)
     endif
         write(eperef_unit)  ua%position(:,:) !(3)
    if(ua3.le.N_unique_atoms) then

         lmax_ch = int(unique_atoms(ua3)%lmax_ch,kind=i4_kind)
         max_order = max(3,lmax_ch+2)

         write(eperef_unit)  unique_atoms(ua3)%r2_ch%n_exponents, &
                             unique_atoms(ua3)%r2_ch%N_uncontracted_fcts, & !(4)
                             unique_atoms(ua3)%r2_ch%N_contracted_fcts 
         write(eperef_unit)  unique_atoms(ua3)%r2_ch%exponents 
         if(unique_atoms(ua3)%r2_ch%N_contracted_fcts.ne.0) &
         write(eperef_unit)  unique_atoms(ua3)%r2_ch%contractions


         write(eperef_unit)  unique_atoms(ua3)%l_ch(0)%n_exponents,&
                           unique_atoms(ua3)%l_ch(0)%N_uncontracted_fcts,&
                           unique_atoms(ua3)%l_ch(0)%N_contracted_fcts

         write(eperef_unit) unique_atoms(ua3)%l_ch(0)%exponents 
         if(unique_atoms(ua3)%l_ch(0)%N_contracted_fcts.ne.0) &
         write(eperef_unit) unique_atoms(ua3)%l_ch(0)%contractions

   if (lmax_ch.gt.0) then
        i_ir_eperef=get_totalsymmetric_irrep()
        n_irreps_eperef=size(unique_atoms(ua3)%symadapt_partner(:,:),dim=1)
        lmax_sym_eperef=size(unique_atoms(ua3)%symadapt_partner(:,:),dim=2)
          write(eperef_unit) n_irreps_eperef,lmax_sym_eperef,i_ir_eperef   !(8)
          write(eperef_unit) &                                  !(9)
          unique_atoms(ua3)%symadapt_partner(:,:)%n_independent_fcts

       do i_l = 1,lmax_ch

         write(eperef_unit) unique_atoms(ua3)%l_ch(i_l)%n_exponents, &
                         unique_atoms(ua3)%l_ch(i_l)%n_uncontracted_fcts , &
                            unique_atoms(ua3)%l_ch(i_l)%N_contracted_fcts  
         write(eperef_unit) unique_atoms(ua3)%l_ch(i_l)%exponents !(11)
         if(unique_atoms(ua3)%l_ch(i_l)%N_contracted_fcts.ne.0) &
         write(eperef_unit) unique_atoms(ua3)%l_ch(i_l)%contractions

      enddo  !  (i_l)?x AG
   endif ! lmax_ch.gt.0
!======================= Store data to be used in  subroutine "calc_sym_coef_eperef"
 ang_momentum_symadapt: do i_l=1,lmax_ch
    n_indep_fcts =  &
                    unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts

         N_partner=size(unique_atoms(ua3)%symadapt_partner(1,i_l)%symadapt,dim=2)
         write(eperef_unit) N_partner
         write(eperef_unit) unique_atoms(ua3)%symadapt_partner(1,i_l)%&
                            symadapt(:,:)%N_fcts
         do i_ind=1,n_indep_fcts
        write(eperef_unit) unique_atoms(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m, &
                           unique_atoms(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c, &
                           unique_atoms(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)% &
                           I_equal_atom
         enddo
 enddo ang_momentum_symadapt
end if ! (ua3.le.N_unique_atoms)
  enddo unique_3
  call returnclose_iounit(eperef_unit)
!!$print*, 'epe_reference ----> written'
  get_epe_reference=.not.get_epe_reference
  end if   ! (make_epe_reference), now read,pack and send epe_reference configuration
         
  get_epe_reference=get_epe_reference .and. epe_iter.eq.1

! FIRST READ EPE_REFERENCE
  if(get_epe_reference) then
         eperef_unit=openget_iounit(file=trim(inpfile("pgepe_reference")), &
                                    form='unformatted',status='old')

         read(eperef_unit)  n_unique_atoms_eperef, &
                            n_timps_eperef, &
                            n_coeff_charge_eperef  !(1)
         allocate(unique_atoms_eperef(n_unique_atoms_eperef),  &     !(2)
                  coeff_charge_eperef(n_coeff_charge_eperef),STAT=epemalloc_stat(2)) !done 
         ASSERT(epemalloc_stat(2).eq.0)
         epemalloc_stat(2)=1
            
         
         if(n_timps_eperef.ne.0) then
            allocate(unique_timps_eperef(n_timps_eperef),STAT=epemalloc_stat(3))
            ASSERT(epemalloc_stat(3).eq.0)
            epemalloc_stat(3)=1
         endif

         read(eperef_unit)  unique_atoms_eperef(:)%lmax_ch,           &
                            unique_atoms_eperef(:)%z,                 &
                            unique_atoms_eperef(:)%zc,                &
                            unique_atoms_eperef(:)%N_equal_atoms,     &  !(2)
                            unique_atoms_eperef(:)%N_glob_cons_ch,    &
                            coeff_charge_eperef(:)

         if(n_timps_eperef.ne.0) &
          read(eperef_unit) unique_timps_eperef(:)%z, &
                            unique_timps_eperef(:)%zc, &
                            unique_timps_eperef(:)%N_equal_atoms

         lmax_ch = maxval(unique_atoms_eperef(:)%lmax_ch)
         allocate(orbitalprojection_ch_eperef( & !(4)
                  -1:lmax_ch, n_unique_atoms_eperef), stat=epemalloc_stat(4)) ! done
         ASSERT(epemalloc_stat(4).eq.0)
         epemalloc_stat(4)=1

         read(eperef_unit ) orbitalprojection_ch_eperef(-1:lmax_ch,:)


! auxiliary arrays to distinguish ua-dependent data
!!$   allocate( n_c_expua(N_unique_atoms_eperef,lmax_ch+2),stat=epemalloc_stat(5) )
!!$   ASSERT(epemalloc_stat(5).eq.0)
!!$   epemalloc_stat(5)=1
   allocate( N_partner_ua(N_unique_atoms_eperef,lmax_ch),stat=epemalloc_stat(6) )
   ASSERT(epemalloc_stat(6).eq.0)
   epemalloc_stat(6)=1
   allocate( n_irreps_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(7) )
   ASSERT(epemalloc_stat(7).eq.0)
   epemalloc_stat(7)=1
   allocate( lmax_sym_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(8) )
   ASSERT(epemalloc_stat(8).eq.0)
   epemalloc_stat(8)=1
   allocate( i_ir_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(9) )
   ASSERT(epemalloc_stat(9).eq.0)
   epemalloc_stat(9)=1

   unique_3_ref: do ua3 = 1,N_unique_atoms_eperef+n_timps_eperef

         if(ua3.gt.N_unique_atoms_eperef) then
            ua_eperef=>unique_timps_eperef(ua3-N_unique_atoms_eperef)
         else
            ua_eperef=>unique_atoms_eperef(ua3)
         end if

         n_equal_c = ua_eperef%N_equal_atoms
         allocate(ua_eperef%position(3,n_equal_c),stat=epemalloc_stat(10))
   ASSERT(epemalloc_stat(10).eq.0)
   epemalloc_stat(10)=1
      
         read(eperef_unit) ua_eperef%position(:,:)  !(3)

   if(ua3.le.N_unique_atoms_eperef) then
         lmax_ch = int(unique_atoms_eperef(ua3)%lmax_ch,kind=i4_kind)
         max_order = max(3,lmax_ch+2)

         if(unique_atoms_eperef(ua3)%N_glob_cons_ch.gt.0) then 
          print*,'alloc(28)', unique_atoms_eperef(ua3)%N_glob_cons_ch
          allocate(unique_atoms_eperef(ua3)%glob_con_ch( &     !(28)
           unique_atoms_eperef(ua3)%N_glob_cons_ch),stat=epemalloc_stat(28)) ! %glob_con_ch
           ASSERT(epemalloc_stat(28).eq.0)
                epemalloc_stat(28)=1
         endif

         read(eperef_unit)  unique_atoms_eperef(ua3)%r2_ch%n_exponents, &
                            unique_atoms_eperef(ua3)%r2_ch%N_uncontracted_fcts, & !(4)
                            unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts

         allocate(unique_atoms_eperef(ua3)%r2_ch%exponents(& !(25)
                         unique_atoms_eperef(ua3)%r2_ch%n_exponents),&
                                                  stat=epemalloc_stat(25)) !done          
         ASSERT(epemalloc_stat(25).eq.0)
         epemalloc_stat(25)=1

         read(eperef_unit)  unique_atoms_eperef(ua3)%r2_ch%exponents !!!&

       if(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts.ne.0) then
         allocate(unique_atoms_eperef(ua3)%r2_ch%contractions(&
                  unique_atoms_eperef(ua3)%r2_ch%n_exponents, &
                  unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts), &
                                                  stat=epemalloc_stat(24))         
         ASSERT(epemalloc_stat(24).eq.0)
                epemalloc_stat(24)=1

         read(eperef_unit) unique_atoms(ua3)%r2_ch%contractions
       end if


         allocate(unique_atoms_eperef(ua3)%l_ch(0:lmax_ch),stat=epemalloc_stat(32))
         ASSERT(epemalloc_stat(32).eq.0)
                epemalloc_stat(32)=1
         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(0)%n_exponents,&
                           unique_atoms_eperef(ua3)%l_ch(0)%N_uncontracted_fcts, &
                           unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts
! the same for r2

         allocate(unique_atoms_eperef(ua3)%l_ch(0)%exponents(&
                  unique_atoms_eperef(ua3)%l_ch(0)%n_exponents) &
                                            ,stat=epemalloc_stat(27)) !l_ch(0)%exponents
         ASSERT(epemalloc_stat(27).eq.0)
                epemalloc_stat(27)=1

         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(0)%exponents

         if(unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts.ne.0) then
          allocate(unique_atoms_eperef(ua3)%l_ch(0)%contractions(&
                             unique_atoms_eperef(ua3)%l_ch(0)%n_exponents, &
                             unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts), &
                             stat=epemalloc_stat(41)) ! %contractions l_ch(0)
          ASSERT(epemalloc_stat(41).eq.0)
                epemalloc_stat(41)=1

         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(0)%contractions
         end if

   if (lmax_ch.gt.0) then
         read(eperef_unit) n_irreps_eperef_ua(ua3),lmax_sym_eperef_ua(ua3), &
                                                   i_ir_eperef_ua(ua3)   !(8)
         n_irreps_eperef=n_irreps_eperef_ua(ua3)
         lmax_sym_eperef=lmax_sym_eperef_ua(ua3)
         i_ir_eperef=i_ir_eperef_ua(ua3)
         allocate(unique_atoms_eperef(ua3)%symadapt_partner(1:n_irreps_eperef,&
                  0:lmax_sym_eperef-1),stat=epemalloc_stat(42)) ! %symadapt_partner
         ASSERT(epemalloc_stat(42).eq.0)
                epemalloc_stat(42)=1

         read(eperef_unit) & !(9)
         unique_atoms_eperef(ua3)%symadapt_partner(:,:)%n_independent_fcts

       do i_l = 1,lmax_ch
         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents, &
                           unique_atoms_eperef(ua3)%l_ch(i_l)%n_uncontracted_fcts, &
                           unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts

         allocate(unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(&
                  unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents), &
                                                stat=epemalloc_stat(45)) ! l_ch(i_l)%exponents
         ASSERT(epemalloc_stat(45).eq.0)
                epemalloc_stat(45)=1


         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(i_l)%exponents !(11)

       if(unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts.ne.0) then
         allocate(unique_atoms_eperef(ua3)%l_ch(i_l)%contractions(&
                           unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents, &
                           unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts), &
                           stat=epemalloc_stat(47)) ! %contractions l_ch(i_l)
         ASSERT(epemalloc_stat(47).eq.0)
                epemalloc_stat(47)=1

         read(eperef_unit) unique_atoms_eperef(ua3)%l_ch(i_l)%contractions
       end if

      enddo  !  (i_l) 
   endif ! lmax_ch.gt.0

 ang_momentum_store: do i_l=1,lmax_ch
   n_indep_fcts=unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
     read(eperef_unit) N_partner_ua(ua3,i_l)
      allocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l) &
          %symadapt(n_indep_fcts,N_partner_ua(ua3,i_l)), stat=epemalloc_stat(35) )  ! symadapt_partner
         ASSERT(epemalloc_stat(35).eq.0)
                epemalloc_stat(35)=1
         read(eperef_unit, iostat=status) unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                           symadapt(:,:)%N_fcts
       do i_ind=1,n_indep_fcts
         n_contributing_fcts = &
                          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%N_fcts
         allocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%m(n_contributing_fcts) &
                         ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%c(n_contributing_fcts) &
                         ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%I_equal_atom(n_contributing_fcts), &
                          stat=epemalloc_stat(36) ) ! %m %c %I_equal_atom
        ASSERT(epemalloc_stat(36).eq.0)
               epemalloc_stat(36)=1

        read(eperef_unit) unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                                                   symadapt(i_ind,1)%m, &
                          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                                                   symadapt(i_ind,1)%c, &
                          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                                                   symadapt(i_ind,1)%I_equal_atom
       enddo
 enddo ang_momentum_store
 end if !(ua3.le.N_unique_atoms)
 enddo unique_3_ref
 if(.not.comm_parallel()) then
!!$   deallocate( n_c_expua, stat=epemalloc_stat(5) ) !
!!$   ASSERT(epemalloc_stat(5).eq.0)

   deallocate( N_partner_ua,n_irreps_eperef_ua, &
        lmax_sym_eperef_ua,i_ir_eperef_ua,stat=epemalloc_stat(6) )  
   ASSERT(epemalloc_stat(6).eq.0)
   epemalloc_stat(7)=0 ! n_irreps_eperef_ua
   epemalloc_stat(8)=0 ! lmax_sym_eperef_ua
   epemalloc_stat(9)=0 ! i_ir_eperef_ua
 end if 
end if ! get_epe_reference
!  print*,'associated(orbitalprojection_ch_eperef)', associated(orbitalprojection_ch_eperef)
  end subroutine epe_read_write_reference

  subroutine epe_send_reference()
    !  Purpose: Sends epe_reference configuration to slaves.
    !           Called by master.
    !------------ Modules used ------------------- ---------------
    implicit none
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind) :: status, length, ion
   integer(kind=i4_kind) :: i_pack,j_pack
   integer(kind=i4_kind) :: end_treated_region
   !------------ Declaration of subroutines used ----------------
   external error_handler
   intrinsic max,maxval
   !------------ Executable code --------------------------------

  call comm_init_send(comm_all_other_hosts,msgtag_epe_do_gradients)

!          Send current ccordinates and gradients to slave
  if(reg_2a_treated) then
     end_treated_region= reg_2a_n_ions  
  else
     end_treated_region=reg_I_n_ions
  end if

  do ion=1,end_treated_region 
         call commpack(reg_I_pg(ion)%vs,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [1]")
         call commpack(reg_I_pg(ion)%gs(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [2]")
         call commpack(reg_I_pg(ion)%vc,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [3]")
         call commpack(reg_I_pg(ion)%gc(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [4]")
         call commpack(reg_I_pg(ion)%rs(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [5]")
         call commpack(reg_I_pg(ion)%rc(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [6]")
  end do
         call commpack(dealloc_epe_ref,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [7]")         
         call commpack(get_epe_reference,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [8]")         
  if(get_epe_reference) then
         call commpack(n_unique_atoms_eperef,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [9]")
! + TIMPs
         call commpack(n_timps_eperef, status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [91]")
         call commpack(n_coeff_charge_eperef,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference : error [10]")
         do i_pack=1,n_unique_atoms_eperef
         call commpack(unique_atoms_eperef(i_pack)%lmax_ch,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [11]")
         call commpack(unique_atoms_eperef(i_pack)%z,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [12]")
         call commpack(unique_atoms_eperef(i_pack)%zc,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [121]")
         call commpack(unique_atoms_eperef(i_pack)%N_equal_atoms,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [13]")
         call commpack(unique_atoms_eperef(i_pack)%N_glob_cons_ch,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [14]")
         end do
!                                  pack TIMPs 
         do i_pack=1, n_timps_eperef
         call commpack(unique_timps_eperef(i_pack)%z, status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [141]")
         call commpack(unique_timps_eperef(i_pack)%zc, status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [142]")
         call commpack(unique_timps_eperef(i_pack)%N_equal_atoms, status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [143]")
         end do
! end TIMPS packing
         call commpack(coeff_charge_eperef(1),n_coeff_charge_eperef,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [15]")
         length= ( maxval(unique_atoms_eperef(:)%lmax_ch)+2 ) * n_unique_atoms_eperef
         call commpack(orbitalprojection_ch_eperef(-1,1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [16]")

   unique_3_ref: do ua3 = 1,N_unique_atoms_eperef+n_timps_eperef

         if(ua3.gt.N_unique_atoms) then
            ua_eperef=>unique_timps_eperef(ua3-N_unique_atoms)
         else
            ua_eperef=>unique_atoms_eperef(ua3)
         endif    

         n_equal_c = ua_eperef%N_equal_atoms
         call commpack(ua_eperef%position(1,1),3*n_equal_c,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [17]") 
    if(ua3.le.N_unique_atoms_eperef) then
         lmax_ch = int(unique_atoms_eperef(ua3)%lmax_ch,kind=i4_kind)
         call commpack(unique_atoms_eperef(ua3)%r2_ch%n_exponents,status)
         if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [18]")
         call commpack(unique_atoms_eperef(ua3)%r2_ch%N_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [19]")
         call commpack(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [20]")
         length=unique_atoms_eperef(ua3)%r2_ch%n_exponents         
         call commpack(unique_atoms_eperef(ua3)%r2_ch%exponents(1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [21]")
  
       if(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts.ne.0) then
         length=unique_atoms_eperef(ua3)%r2_ch%n_exponents * &
                unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts
         call commpack(unique_atoms_eperef(ua3)%r2_ch%contractions(1,1),&
                                                         length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [22]")
       end if
!!$         call commpack(n_c_expua(ua3,1),status)
!!$         if(status.ne.0) call error_handler &
!!$                         ("store_and_send_epe_reference :  error [23]")

         call commpack(unique_atoms_eperef(ua3)%l_ch(0)%n_exponents,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [25]")
         call commpack(unique_atoms_eperef(ua3)%l_ch(0)%N_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [26]")
         call commpack(unique_atoms_eperef(ua3)%l_ch(0)%n_contracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [27]")
! the same for r2
         length=unique_atoms_eperef(ua3)%l_ch(0)%n_exponents
         call commpack(unique_atoms_eperef(ua3)%l_ch(0)%exponents(1),&
                                                      length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [28]")
     if(unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts.ne.0) then
         length=unique_atoms_eperef(ua3)%l_ch(0)%n_exponents * &
                unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts
         call commpack(unique_atoms_eperef(ua3)%l_ch(0)%contractions(1,1),&
                        length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [29]")
     end if
!!$         call commpack(n_c_expua(ua3,2),status)
!!$         if(status.ne.0) call error_handler &
!!$                         ("store_and_send_epe_reference :  error [30]")

   if (lmax_ch.gt.0) then
         call commpack(n_irreps_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [32]")
         n_irreps_eperef=n_irreps_eperef_ua(ua3)
         call commpack(lmax_sym_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [33]")
         lmax_sym_eperef=lmax_sym_eperef_ua(ua3)
         call commpack(i_ir_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [34]")
         i_ir_eperef=i_ir_eperef_ua(ua3)
         do i_pack=1,n_irreps_eperef
         do j_pack=0,lmax_sym_eperef-1
         call commpack(unique_atoms_eperef(ua3) &
                     %symadapt_partner(i_pack,j_pack)%n_independent_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [35]")
         end do  ! i_pack
         end do  ! j_pack

       do i_l = 1,lmax_ch
         call commpack(unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [36]")
         call commpack(unique_atoms_eperef(ua3)%l_ch(i_l)%n_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [37]")
         call commpack(unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [38]")
         length=unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents
         call commpack(unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [39]")
      if(unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts.ne.0) then
         length=unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents * &
                unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts
         call commpack(unique_atoms_eperef(ua3)%l_ch(i_l)%contractions(1,1), &
                       length,1,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [40]")
      end if
!!$       do independent_fct=1,unique_atoms_eperef(ua3)&
!!$                             %symadapt_partner(1,i_l)%N_independent_fcts
!!$         call commpack(n_c_expua(ua3,i_l+2),status)
!!$         if(status.ne.0) call error_handler &
!!$                         ("store_and_send_epe_reference :  error [41]")
!!$
!!$!         call commpack(unique_atoms_eperef(ua3)% &
!!$!                       renormaliation_partner_ch(i_l)% &
!!$!                       renorm(independent_fct)%c_exp(1),n_c_expua(ua3,i_l+2),1,status)
!!$!         if(status.ne.0) call error_handler &
!!$!                         ("store_and_send_epe_reference :  error [42]")         
!!$
!!$       enddo
      enddo  !  (i_l) 
    endif ! lmax_ch.gt.0
!======================= Store data to be used in  subroutine "calc_sym_coef_eperef"
 ang_momentum_store: do i_l=1,lmax_ch
    n_indep_fcts =  &
              unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
         N_partner=N_partner_ua(ua3,i_l) 
         call commpack(N_partner,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [43]")
         do i_pack=1,n_indep_fcts
         do j_pack=1,N_partner
         call commpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                       symadapt(i_pack,j_pack)%N_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [44]")
         end do ! i_pack
         end do ! j_pack
       do i_ind=1,n_indep_fcts
         n_contributing_fcts = &
                          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%N_fcts
        length = n_contributing_fcts
        call commpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                                     symadapt(i_ind,1)%m(1), length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [45]")                                     
        call commpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                                     symadapt(i_ind,1)%c(1), length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [46]")                                     
        call commpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                           symadapt(i_ind,1)%I_equal_atom(1),length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [47]")
       enddo
 enddo ang_momentum_store
 end if ! (ua3.le.N_unique_atoms_eperef) then
 enddo unique_3_ref
!!$   deallocate( n_c_expua ,stat=epemalloc_stat(5) )
!!$   ASSERT(epemalloc_stat(5).eq.0)

   deallocate( N_partner_ua,stat=epemalloc_stat(6) )
   ASSERT(epemalloc_stat(6).eq.0)
   deallocate( n_irreps_eperef_ua,stat=epemalloc_stat(7) )
   ASSERT(epemalloc_stat(7).eq.0)
   deallocate( lmax_sym_eperef_ua,stat=epemalloc_stat(8) )
   ASSERT(epemalloc_stat(8).eq.0)
   deallocate( i_ir_eperef_ua,stat=epemalloc_stat(9) )  
   ASSERT(epemalloc_stat(9).eq.0)
 end if ! get_epe_reference
  call comm_send()

  end subroutine epe_send_reference
  !*************************************************************

  !*************************************************************
  subroutine epe_receive_reference()
    !  Purpose: Receives stored reference configuration 
    !           from master. Called by slaves.
    !------------ Modules used ------------------- ---------------
    implicit none
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind)              :: status, length, ion
   integer(kind=i4_kind)              :: end_treated_region
   integer(kind=i4_kind)              :: N_partner,i_pack,j_pack
   !------------ Declaration of subroutines used ----------------
   external error_handler
   intrinsic max, maxval
   !------------ Executable code --------------------------------
!!$  call write_to_output_units("==================================")
!!$  call write_to_output_units(" Slave is receiving epe_reference ")
!!$  call write_to_output_units("==================================")
  if(reg_2a_treated) then
     end_treated_region=reg_2a_n_ions
     else
     end_treated_region= reg_I_n_ions
  end if
  
  do ion=1,end_treated_region
         call communpack(reg_I_pg(ion)%vs,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [1]")
         call communpack(reg_I_pg(ion)%gs(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [2]")
         call communpack(reg_I_pg(ion)%vc,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [3]")
         call communpack(reg_I_pg(ion)%gc(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [4]")
         call communpack(reg_I_pg(ion)%rs(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [5]")
         call communpack(reg_I_pg(ion)%rc(1),3,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [6]")
  end do
         call communpack(dealloc_epe_ref,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [7]")
         call communpack(get_epe_reference,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [8]")   
    if(get_epe_reference) then 
         call communpack(n_unique_atoms_eperef,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [9]")
!+ TIMPs
         call communpack(n_timps_eperef, status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference : error [91]")
         call communpack(n_coeff_charge_eperef,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [10]")
         allocate(unique_atoms_eperef(n_unique_atoms_eperef),  &
                  coeff_charge_eperef(n_coeff_charge_eperef),stat=epemalloc_stat(2)) !done
         ASSERT(epemalloc_stat(2).eq.0)
                epemalloc_stat(2)=1
     do i_pack=1,n_unique_atoms_eperef        
         call communpack(unique_atoms_eperef(i_pack)%lmax_ch,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [11]")
         call communpack(unique_atoms_eperef(i_pack)%z,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [12]")
         call communpack(unique_atoms_eperef(i_pack)%zc,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [121]")
         call communpack(unique_atoms_eperef(i_pack)%N_equal_atoms,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [13]")
         call communpack(unique_atoms_eperef(i_pack)%N_glob_cons_ch,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [14]")
     end do ! i_pack
!                                unpack TIMPs 
        if(n_timps_eperef.ne.0) then
         allocate( unique_timps_eperef(n_timps_eperef),stat=epemalloc_stat(3) )
         ASSERT(epemalloc_stat(3).eq.0)
                epemalloc_stat(3)=1
        endif

         do i_pack=1, n_timps_eperef
         call communpack(unique_timps_eperef(i_pack)%z, status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [141]")
         call communpack(unique_timps_eperef(i_pack)%zc, status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [142]")
         call communpack(unique_timps_eperef(i_pack)%N_equal_atoms, status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [143]")
         end do
! end TIMPS unpacking
         call communpack(coeff_charge_eperef(1),n_coeff_charge_eperef,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [15]")
         allocate(orbitalprojection_ch_eperef( &
                  -1:maxval(unique_atoms_eperef(:)%lmax_ch),&
                   n_unique_atoms_eperef),stat=epemalloc_stat(4)) ! done        
         ASSERT(epemalloc_stat(4).eq.0)
                epemalloc_stat(4)=1

         length= ( maxval(unique_atoms_eperef(:)%lmax_ch)+2 ) * n_unique_atoms_eperef
         call communpack(orbitalprojection_ch_eperef(-1,1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [16]")
   lmax_ch = maxval(unique_atoms_eperef(:)%lmax_ch)
! auxiliary arrays to distinguish ua-dependent data
!!$   allocate( n_c_expua(N_unique_atoms_eperef,lmax_ch+2),stat=epemalloc_stat(5) )
!!$   ASSERT(epemalloc_stat(5).eq.0)
!!$          epemalloc_stat(5)=1
   allocate( N_partner_ua(N_unique_atoms_eperef,lmax_ch),stat=epemalloc_stat(6))
   ASSERT(epemalloc_stat(6).eq.0)
   epemalloc_stat(6)=1
   allocate( n_irreps_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(7) )
   ASSERT(epemalloc_stat(7).eq.0)
   epemalloc_stat(7)=1
   allocate( lmax_sym_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(8) )
   ASSERT(epemalloc_stat(8).eq.0)
   epemalloc_stat(8)=1

   allocate( i_ir_eperef_ua(N_unique_atoms_eperef),stat=epemalloc_stat(9) )
   ASSERT(epemalloc_stat(9).eq.0)
   epemalloc_stat(9)=1
        
   unique_3_ref: do ua3 = 1,N_unique_atoms_eperef


         if(ua3.gt.N_unique_atoms) then
            ua_eperef=>unique_timps_eperef(ua3-N_unique_atoms)
         else
            ua_eperef=>unique_atoms_eperef(ua3)
         endif   

         n_equal_c = ua_eperef%N_equal_atoms

         allocate(unique_atoms_eperef(ua3)%position(3,n_equal_c),stat=epemalloc_stat(10))
         ASSERT(epemalloc_stat(10).eq.0)
         epemalloc_stat(10)=1

         call communpack(unique_atoms_eperef(ua3)%position(1,1), &
                                          3*n_equal_c,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [17]")
  if(ua3.le.N_unique_atoms_eperef) then
   lmax_ch = int(unique_atoms_eperef(ua3)%lmax_ch,kind=i4_kind)
   if(unique_atoms_eperef(ua3)%N_glob_cons_ch.gt.0) then
         allocate(unique_atoms_eperef(ua3)%glob_con_ch( &
                  unique_atoms_eperef(ua3)%N_glob_cons_ch), &
                                         stat=epemalloc_stat(28))   ! %glob_con_ch
         ASSERT(epemalloc_stat(28).eq.0)
                epemalloc_stat(28)=1
    endif

         call communpack(unique_atoms_eperef(ua3)%r2_ch%n_exponents,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [18]")
         allocate(unique_atoms_eperef(ua3)%r2_ch%exponents(&
                            unique_atoms_eperef(ua3)%r2_ch%n_exponents), &
                                                      stat=epemalloc_stat(25)) !done
         ASSERT(epemalloc_stat(25).eq.0)
         epemalloc_stat(25)=1

         call communpack(unique_atoms_eperef(ua3)%r2_ch%N_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [19]")
         call communpack(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts,status)         
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [20]")
         length=unique_atoms_eperef(ua3)%r2_ch%n_exponents         
         call communpack(unique_atoms_eperef(ua3)%r2_ch%exponents(1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [21]")         
         if(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts.ne.0) then
         allocate(unique_atoms_eperef(ua3)%r2_ch%contractions(&
                  unique_atoms_eperef(ua3)%r2_ch%n_exponents, &
                  unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts), &
                  stat=epemalloc_stat(24))         
         ASSERT(epemalloc_stat(24).eq.0)
         epemalloc_stat(24)=1

         length=unique_atoms_eperef(ua3)%r2_ch%n_exponents * &
                unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts
         call communpack(unique_atoms_eperef(ua3)%r2_ch%contractions(1,1), &
                                                         length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [22]")
         end if

         allocate(unique_atoms_eperef(ua3)%l_ch(0:lmax_ch),stat=epemalloc_stat(32))
         ASSERT(epemalloc_stat(32).eq.0)
                epemalloc_stat(32)=1
         call communpack(unique_atoms_eperef(ua3)%l_ch(0)%n_exponents,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [25]")
         call communpack(unique_atoms_eperef(ua3)%l_ch(0)%N_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [26]")
         call communpack(unique_atoms_eperef(ua3)%l_ch(0)%n_contracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [27]")
! the same for r2
         allocate(unique_atoms_eperef(ua3)%l_ch(0)%exponents(&
                  unique_atoms_eperef(ua3)%l_ch(0)%n_exponents), &
                   stat=epemalloc_stat(27)) !%exponents l_ch(0)
         ASSERT(epemalloc_stat(27).eq.0)
                epemalloc_stat(27)=1

         length=unique_atoms_eperef(ua3)%l_ch(0)%n_exponents
         call communpack(unique_atoms_eperef(ua3)%l_ch(0)%exponents(1),&
                                                        length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [28]")
         if(unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts.ne.0) then
         allocate(unique_atoms_eperef(ua3)%l_ch(0)%contractions(&
                             unique_atoms_eperef(ua3)%l_ch(0)%n_exponents, &
                             unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts) , &
                             stat=epemalloc_stat(41)) ! %contractions
         ASSERT(epemalloc_stat(41).eq.0)
                epemalloc_stat(41)=1

         length=unique_atoms_eperef(ua3)%l_ch(0)%n_exponents * &
                unique_atoms_eperef(ua3)%l_ch(0)%N_contracted_fcts
         call communpack(unique_atoms_eperef(ua3)%l_ch(0)%contractions(1,1),&
                        length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [29]")
         end if

   if (lmax_ch.gt.0) then
         call communpack(n_irreps_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [32]")
         n_irreps_eperef=n_irreps_eperef_ua(ua3)
         call communpack(lmax_sym_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [33]")
         lmax_sym_eperef=lmax_sym_eperef_ua(ua3)
         call communpack(i_ir_eperef_ua(ua3),status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [34]")
         i_ir_eperef=i_ir_eperef_ua(ua3)
         allocate(unique_atoms_eperef(ua3)&
          %symadapt_partner(1:n_irreps_eperef,0:lmax_sym_eperef-1),stat=epemalloc_stat(42))
         ASSERT(epemalloc_stat(42).eq.0)
                epemalloc_stat(42)=1

         length=n_irreps_eperef * lmax_sym_eperef
         do i_pack=1,n_irreps_eperef
         do j_pack=0,lmax_sym_eperef-1         
         call communpack(unique_atoms_eperef(ua3) &
                     %symadapt_partner(i_pack,j_pack)%n_independent_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [35]")
         end do ! i_pack
         end do ! j_pack
       do i_l = 1,lmax_ch
         call communpack(unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [36]")
         call communpack(unique_atoms_eperef(ua3)%l_ch(i_l)%n_uncontracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [37]")
         call communpack(unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [38]")
         allocate(unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(&
                  unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents), &
                  stat=epemalloc_stat(45)) ! l_ch(i_l)%exponents
         ASSERT(epemalloc_stat(45).eq.0)
                epemalloc_stat(45)=1

         length=unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents
         call communpack(unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(1),length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [39]")
       if(unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts.ne.0) then
         allocate(unique_atoms_eperef(ua3)%l_ch(i_l)%contractions &
                                (unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents, &
                           unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts), &
                           stat=epemalloc_stat(47))  ! %contractions  l_ch(i_l)
         ASSERT(epemalloc_stat(47).eq.0)
                epemalloc_stat(47)=1
         length=unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents * &
                unique_atoms_eperef(ua3)%l_ch(i_l)%N_contracted_fcts
         call communpack(unique_atoms_eperef(ua3)%l_ch(i_l)%contractions(1,1), &
                       length,1,status)
         if(status.ne.0) call error_handler &
                         ("receive_epe_reference :  error [40]")
       end if
      enddo  !  (i_l) 
    endif ! lmax_ch.gt.0
!======================= Receive data to be used in  subroutine "calc_sym_coef_eperef"
 ang_momentum_symadapt: do i_l=1,lmax_ch
    n_indep_fcts =  &
                  unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
         call communpack(N_partner_ua(ua3,i_l),status)
         N_partner=N_partner_ua(ua3,i_l)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [43]")
          
         allocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                  symadapt(n_indep_fcts,N_partner),stat=epemalloc_stat(35) ) !  %symadapt
         ASSERT(epemalloc_stat(35).eq.0)
                epemalloc_stat(35)=1

         do i_pack=1,n_indep_fcts
         do j_pack=1, N_partner
         call communpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                       symadapt(i_pack,j_pack)%N_fcts,status)
         if(status.ne.0) call error_handler &
                         ("store_and_send_epe_reference :  error [44]")
         end do ! i_pack
         end do ! j_pack                        
       do i_ind=1,n_indep_fcts
         n_contributing_fcts = &
                          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%N_fcts
         allocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%m(n_contributing_fcts) &
                         ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%c(n_contributing_fcts) &
                         ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
                          symadapt(i_ind,1)%I_equal_atom(n_contributing_fcts), &
               stat=epemalloc_stat(36) ) ! symadapt(i_ind,1)%m, symadapt(i_ind,1)%c, symadapt(i_ind,1)%I_equal_atom
               ASSERT(epemalloc_stat(36).eq.0)
                      epemalloc_stat(36)=1
        length = n_contributing_fcts 
        call communpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                                     symadapt(i_ind,1)%m(1), length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [45]")                                     
        call communpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                                     symadapt(i_ind,1)%c(1), length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [46]")                                     
        call communpack(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)% &
                           symadapt(i_ind,1)%I_equal_atom(1),length,1,status)
        if(status.ne.0) call error_handler &
                        ("store_and_send_epe_reference :  error [47]")
       enddo
 enddo ang_momentum_symadapt
 end if ! (ua3.le.N_unique_atoms_eperef) then
 enddo unique_3_ref
!!$   deallocate( n_c_expua , N_partner_ua , n_irreps_eperef_ua ,stat=epemalloc_stat(5))
   deallocate( N_partner_ua , n_irreps_eperef_ua ,stat=epemalloc_stat(5))
   ASSERT(epemalloc_stat(5).eq.0)
   epemalloc_stat(6)=0 ! N_partner_ua
   epemalloc_stat(7)=0 ! n_irreps_eperef_ua
   deallocate( lmax_sym_eperef_ua,stat=epemalloc_stat(8))
   ASSERT(epemalloc_stat(8).eq.0)
   deallocate( i_ir_eperef_ua,stat=epemalloc_stat(9) )  
   ASSERT(epemalloc_stat(9).eq.0)
!!$          print*, 'Slave : EPE_reference  is  received, OK' 
!!$          call write_to_output_units("Slave : EPE_reference  is  received, OK ")
 end if ! << get_epe_reference >>

  end subroutine epe_receive_reference
  !*************************************************************

  !*************************************************************
  subroutine epe_collect_gradients()
    !  Purpose: Receives epe_gradients and epe_energies calculated
    !           by slaves  to use in "main_epe" subroutine.
    !           Called by master.
    !------------ Modules used ----------------------------------
    implicit none
   !------------ Declaration of local variables -----------------
   integer :: n_treated_ions 
   integer             :: status, n_received
   logical             :: all_gradients_received
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
!!$    call write_to_output_units("Master: collect gradients")
    n_received=0
    all_gradients_received=.false.
    do while(.not.all_gradients_received )
       if(reg_2a_treated) then
        n_treated_ions= reg_2a_n_ions
        else
        n_treated_ions = reg_I_n_ions
     end if
     
    allocate ( epe_send_array(n_treated_ions,8),stat=epemalloc_stat(39) )
    ASSERT(epemalloc_stat(39).eq.0)
           epemalloc_stat(39)=1
    allocate ( epe_convergence_energies(3),stat=epemalloc_stat(38) )
    ASSERT(epemalloc_stat(38).eq.0)
           epemalloc_stat(38)=1

    do while ( .not. comm_save_recv_nonblocking(comm_all_other_hosts, &
                                         msgtag_epe_grad_done) )
    end do
    n_received=n_received+1
    call communpack(epe_send_array(1,1),8*n_treated_ions,1,status)
          if(status.ne.0) call error_handler &
                         ("epe_collect_gradients :  error [1]")
    call communpack(epe_convergence_energies,3,1,status)
          if(status.ne.0) call error_handler &
                         ("epe_collect_gradients :  error [2]")
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vs = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vs +    &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,1)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(1) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(1) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,2)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(2) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(2) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,3)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(3) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(3) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,4)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vc = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vc +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,5)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(1) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(1) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,6)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(2) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(2) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,7)
    reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(3) = &
                  reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(3) +  &
                  epe_send_array(1+n_epe_vacancies:n_treated_ions,8)
    epg_reference_reg_I = epg_reference_reg_I +                     &
                          epe_convergence_energies(1)
    epg_cluster_reg_I  =  epg_cluster_reg_I   +                     &
                          epe_convergence_energies(2)
    epg_nuc            =  epg_nuc +                                 &
                          epe_convergence_energies(3)
    deallocate( epe_send_array,stat=epemalloc_stat(39) )
    ASSERT(epemalloc_stat(39).eq.0)
    deallocate( epe_convergence_energies,stat=epemalloc_stat(38) )
    ASSERT(epemalloc_stat(38).eq.0)
    all_gradients_received=n_received==comm_get_n_processors()-1
    end do
!!$       call write_to_output_units('Master: gradients received')                
  end subroutine epe_collect_gradients 
  !*************************************************************


  !*************************************************************
  subroutine epe_field_and_forces_par()
  !----------------------------------------------------------------
  ! Purpose: calculate coulomb  field and gradients of the field 
  !          in positions of EPE centers. On the first iteration
  !          of EPE relaxation loop the data storadge for pgepe_reference
  !          is allocated  and the pgepe_reference data are read. If
  !          pgepe_reference file is missing then it is created and
  !        the current charge density distribution  is considered 
  !          as EPE reference charge distribution.
  !          Symmetry adaption for fitfunctions
  !          is included.
  !
  !  Subroutine called by: main_epe, main_slave
  !
  !  References: ...
  !
  !  Author: VN
  !  Date: ...
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !== Interrupt end of public interface of module =================

  !================================================================
  ! End of public interface of module
  !================================================================
!..............................................................................
! << OUTPUT ARRAYS >>
! ===================
! prim_int_2cob_ol_grad  (  1:N) : dsym/dRa < xi_i | 1      | xi_j >
! prim_int_2cob_ol_grad  (N+1:M) : dsym/dRb < xi_i | 1      | xi_j >
! << relativistic calculation >>
! prim_int_2cob_kin_grad (  1:N) : dsym/dRa < xi_i | T      | xi_j >
! prim_int_2cob_kin_grad (N+1:M) : dsym/dRb < xi_i | T      | xi_j >
! prim_int_2cob_nuc_grad ( ia: ) : dsym/dRa < xi_i | V_nuc  | xi_j >
! prim_int_2cob_nuc_grad ( ib: ) : dsym/dRb < xi_i | V_nuc  | xi_j >
! prim_int_2cob_nuc_grad ( ic: ) : dsym/dRc < xi_i | V_nuc  | xi_j >
! prim_int_2cob_pvsp_grad( ia: ) : dsym/dRa < xi_i | V_pvsp | xi_j >
! prim_int_2cob_pvsp_grad( ib: ) : dsym/dRb < xi_i | V_pvsp | xi_j >
! prim_int_2cob_pvsp_grad( ic: ) : dsym/dRc < xi_i | V_pvsp | xi_j >
! prim_int_3cob_grad     ( ia: ) : dsym/dRa < xi_i | V_H    | xi_j >
! prim_int_3cob_grad     ( ib: ) : dsym/dRb < xi_i | V_H    | xi_j >
! prim_int_3cob_grad     ( ic: ) : dsym/dRc < xi_i | V_H    | xi_j >
! << non-relativistic calculation with total gradients >>
! prim_int_3cob_grad     ( ia: ) : dsym/dRa < xi_i | T + V_nuc + V_H | xi_j >
! prim_int_3cob_grad     ( ib: ) : dsym/dRb < xi_i | T + V_nuc + V_H | xi_j >
! prim_int_3cob_grad     ( ic: ) : dsym/dRc < xi_i |     V_nuc + V_H | xi_j >
!
! Model_Density_Approach
! ~~~~~~~~~~~~~~~~~~~~~~
! << relativistic calculation >>
! prim_int_3cob_grad     ( toff+ia: ) : dsym/dRa <xi_i|        V_H+V_X,t|xi_j>
! prim_in   t_3cob_grad     ( toff+ib: ) : dsym/dRb <xi_i|        V_H+V_X,t|xi_j>
! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|        V_H+V_X,t|xi_j>
! << non-relativistic calculation with total gradients >>
! prim_int_3cob_grad     ( toff+ia: ) : dsym/dRa <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
! prim_int_3cob_grad     ( toff+ib: ) : dsym/dRb <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|  V_nuc+V_H+V_X,t|xi_j>
! << non-relativistic calculation with split gradients >>
! prim_int_2cob_ks_grad  (toff+  1:N) : dsym/dRa <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
! prim_int_2cob_ks_grad  (toff+N+1:M) : dsym/dRb <xi_i|T+V_nuc+V_H+V_X,t|xi_j>
! prim_int_3cob_nuc_grad (      ic: ) : dsym/dRc <xi_i|  V_nuc          |xi_j>
! prim_int_3cob_coul_grad(      ic: ) : dsym/dRc <xi_i|        V_H      |xi_j>
! prim_int_3cob_grad     ( toff+ic: ) : dsym/dRc <xi_i|            V_X,t|xi_j>
!
! << WORKING ARRAYS >>
! ====================
! nuc_grad      (:,1:grad_dim) : dsym/dRc < xi_i | V_nuc  | xi_j >
! help_arr_c    (:,1:3       ) :    d/dRc < xi_i | V_pvsp | xi_j >
! help_arr      (:,1:3       ) :    d/dRb [ xi_i | f_k    | xi_j ]
! coul_int_c    (  1:grad_dim) : dsym/dRb [ xi_i | f_k    | xi_j ]
!
!..............................................................................

  !------------ Declaration of subroutines ------------------------
  external error_handler
  intrinsic max,maxval
  !------------ Declaration of local constants --------------------
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: one=1.0_r8_kind,&
       two=2.0_r8_kind,&
       three=3.0_r8_kind, &
       four=4.0_r8_kind
  real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  !------------ Declaration of local variables --------------------
  type(unique_atom_type), pointer :: unique_timps_eperef(:)
  type(unique_atom_type), pointer  :: ua,ua_eperef 
  integer(kind=i4_kind)          :: ncexps,grad_dim
  real(kind=r8_kind),pointer     :: cexps(:)
  ! mapping of exponents to one dimension and cutoff of small integrals
  integer(kind=i4_kind) :: num 
  integer(kind=i4_kind) :: n_treated_ions
  ! help factors
  integer(kind=i4_kind) :: i_grad
  logical                        :: model_density, &
                                    spin_polarized, dealloc_coeff_charge
  logical,allocatable, dimension(:,:):: cutoff
  ! help factors for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:):: gamma_help
  real(kind=r8_kind),allocatable :: gamma_arg2(:),gamma_arg(:,:),help_arr(:,:)
  ! help arrays for solid harmonics
  real(kind=r8_kind),allocatable,dimension(:,:)   :: yl_arg
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: yl_arr
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: yl_arr_grad
  ! help arrays for symmetry adaption
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: sym_coef1
  real(kind=r8_kind),allocatable,dimension(:,:,:,:,:) :: sym_coef2
  type(arrmat4),pointer ::prim_epe_grad(:)
  type(arrmat4),pointer ::prim_epe_gref(:)

  !  ! 3-center integrals
  !  type(three_center_l_grad)  :: grad_coul_1, grad_coul_2
  type(three_center_l),allocatable  :: coul_int_c(:)

  ! help pointers to unique atoms module data
  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coef(:),rotmat(:,:)
  real(kind=r8_kind),allocatable :: nuc_grad(:,:),nuc_v(:)
  real(kind=r8_kind),allocatable :: nuc_gref(:,:),nuc_vref(:)

  integer(kind=i4_kind) :: i,ua3,max_order,&
       max_gamma,lmax_ch,k, & 
       n_equal_c,lm,n_indep_max,n_indep,i_l,i_sum, &
       n_indep_fcts,n_contributing_fcts,i_ind,i_cont,i_ea3
  real(kind=r8_kind),dimension(3)   :: xc 
  real(kind=r8_kind)    :: zc
  integer(kind=i4_kind) :: g_1,g_2  ! bounds for gradients to store
  integer(kind=4)       :: ix, i_epe  !, i_g
  !----------------------------------------------------------------

  !------------ Executable code -----------------------------------

  call gamma_setup(16) ! cp. default numj=17 in old vers. of gamma_module

  grad_dim=3   
  model_density = options_xcmode() == xcmode_model_density
  spin_polarized = .not. options_spin_restricted()

!  grad_s=>( . . . g_1........g_2 . . . )     
     if( comm_i_am_master() ) then
       num=2*n_grads_master
       g_1=n_epe_vacancies+1+( comm_get_n_processors()-1 )*n_grads_slave
       g_2=g_1 + n_grads_master - 1 
     else 
       num=2*n_grads_slave
       g_1=n_epe_vacancies+1+( comm_myindex()-2 )*n_grads_slave
       g_2=g_1 + n_grads_slave - 1
     end if

     allocate(cutoff(num,1),STAT=epemalloc_stat(11))                 !deal done
     ASSERT(epemalloc_stat(11).eq.0)
     epemalloc_stat(11)=1

     cutoff=.true.
     allocate(prim_epe_grad(0:3),prim_epe_gref(0:3),STAT=epemalloc_stat(12)) ! 0-v;1,2,3-grads
     ASSERT(epemalloc_stat(12).eq.0)
     epemalloc_stat(12)=1
     allocate(prim_epe_grad(0)%m(num,1,1,1), &
              prim_epe_grad(1)%m(num,1,1,1), &
              prim_epe_grad(2)%m(num,1,1,1), &
              prim_epe_grad(3)%m(num,1,1,1), &
                          STAT=epemalloc_stat(13)) !done
     ASSERT(epemalloc_stat(13).eq.0)
     epemalloc_stat(13)=1
     allocate(prim_epe_gref(0)%m(num,1,1,1), &
              prim_epe_gref(1)%m(num,1,1,1), &
              prim_epe_gref(2)%m(num,1,1,1), &
              prim_epe_gref(3)%m(num,1,1,1), &
                          STAT=epemalloc_stat(14)) !done
     ASSERT(epemalloc_stat(14).eq.0)
     epemalloc_stat(14)=1
     allocate (gamma_arg(num,3),gamma_arg2(num),help_arr(num,3),&
                                             STAT=epemalloc_stat(15))     !done
     ASSERT(epemalloc_stat(15).eq.0)
     epemalloc_stat(15)=1
     allocate(nuc_grad(num,3),nuc_v(num),nuc_gref(num,3),nuc_vref(num), &
                                                     stat=epemalloc_stat(16)) !done
     ASSERT(epemalloc_stat(16).eq.0)
     epemalloc_stat(16)=1
     max_gamma=2
     allocate(coul_int_c(0:3),STAT=epemalloc_stat(17)) !done
     ASSERT(epemalloc_stat(17).eq.0)
     epemalloc_stat(17)=1

  lmax_ch = maxval(unique_atoms_eperef(:)%lmax_ch)
  max_order = max(lmax_ch+2,3)
  allocate (gamma_help(num,max_order),STAT=epemalloc_stat(18)) !done
     ASSERT(epemalloc_stat(18).eq.0)
     epemalloc_stat(18)=1
     prim_epe_grad(0)%m=0.0_r8_kind
     prim_epe_gref(0)%m=0.0_r8_kind
     do i_grad=1,3
        prim_epe_grad(i_grad)%m=0.0_r8_kind
        prim_epe_gref(i_grad)%m=0.0_r8_kind
        gamma_arg(1:num/2,i_grad)     = reg_I_pg(g_1:g_2)%rs(i_grad)
        gamma_arg(num/2+1:num,i_grad) = reg_I_pg(g_1:g_2)%rc(i_grad)
     end do

        help_arr=zero
        nuc_grad=zero   
        nuc_v=zero      
        nuc_vref=zero   
        nuc_gref=zero

!!$     print*,'loop unique_3_ref'
  unique_3_ref: do ua3 = 1,N_unique_atoms_eperef+n_timps_eperef

!!$if( comm_i_am_master() ) then
!!$  if(associated(unique_timps_eperef) ) &
!!$       print*, 'Master: unique_timps size ', size(unique_timps_eperef), n_timps_eperef
!!$else
!!$  if(associated(unique_timps_eperef) ) & 
!!$       print*, 'Slave:  unique_timps size ', size(unique_timps_eperef), n_timps_eperef
!!$end if

     if(ua3.gt.N_unique_atoms_eperef) then
        ua_eperef=>unique_timps_eperef(ua3-N_unique_atoms_eperef)
     else
        ua_eperef=>unique_atoms_eperef(ua3)
     end if 

     zc = unique_atoms_eperef(ua3)%Z
     zcc = ua_eperef%ZC
     n_equal_c = ua_eperef%N_equal_atoms
     
  equal_3_ref: do i_ea3 = 1,n_equal_c
         xc = ua_eperef%position(:,i_ea3)

         gamma_arg2 = ( (gamma_arg(:,1)-xc(1))**2 + &
                        (gamma_arg(:,2)-xc(2))**2 + &
                        (gamma_arg(:,3)-xc(3))**2   )
          nuc_vref(:)=nuc_vref(:)+ (zc-zcc)/sqrt(gamma_arg2)
          do i_grad=1,3
          nuc_gref(:,i_grad)=nuc_gref(:,i_grad) + &
                             (zc-zcc)*(gamma_arg(:,i_grad)-xc(i_grad)) / &
                             sqrt(gamma_arg2)/gamma_arg2
          end do        
  enddo equal_3_ref

     ! now fitfunctions  ***************
     uaeperef: if(ua3.le.N_unique_atoms_eperef) then

     ! Orbital gradients:
     lmax_ch = int(unique_atoms_eperef(ua3)%lmax_ch,kind=i4_kind)
     max_order = max(3,lmax_ch+2)

  do i_grad=0,3     
        allocate(coul_int_c(i_grad)%l(-1:lmax_ch),stat=epemalloc_stat(19))
        ASSERT(epemalloc_stat(19).eq.0)
        epemalloc_stat(19)=1
        ncexps = unique_atoms_eperef(ua3)%r2_ch%n_exponents
        allocate(coul_int_c(i_grad)%l(-1)%m&
                           (num,ncexps,1,1,1),stat=epemalloc_stat(20))
        ASSERT(epemalloc_stat(20).eq.0)
        epemalloc_stat(20)=1

        ncexps = unique_atoms_eperef(ua3)%l_ch(0)%n_exponents
        allocate(coul_int_c(i_grad)%l(0)%m&
                           (num,ncexps,1,1,1),stat=epemalloc_stat(20))
        ASSERT(epemalloc_stat(20).eq.0)
        epemalloc_stat(20)=1
        coul_int_c(i_grad)%l(0)%m=0.0_r8_kind
        coul_int_c(i_grad)%l(-1)%m=0.0_r8_kind
  enddo ! i_grad=0,3

     if (lmax_ch.gt.0) then

        allocate(yl_arr(num,(lmax_ch+1)**2,n_equal_c),  &
                 yl_arr_grad(num,n_equal_c,(lmax_ch+1)**2,3),  &
                                            STAT=epemalloc_stat(21))
        ASSERT(epemalloc_stat(21).eq.0)
        epemalloc_stat(21)=1
        yl_arr = zero
        yl_arr_grad = zero

        n_indep_max = maxval(unique_atoms_eperef(ua3)%symadapt_partner&
             (1,:)%n_independent_fcts)
 
        allocate(sym_coef1(num,n_equal_c,n_indep_max,lmax_ch),&
                 sym_coef2(num,n_equal_c,n_indep_max,lmax_ch,3),&
                                             STAT=epemalloc_stat(22))
        ASSERT(epemalloc_stat(22).eq.0)
        epemalloc_stat(22)=1

        sym_coef1 = zero
        sym_coef2 = zero

        do i_l = 1,lmax_ch
           ncexps = unique_atoms_eperef(ua3)%l_ch(i_l)%n_exponents
           n_indep =  &
           unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
           do i_grad=0,3   
              allocate (coul_int_c(i_grad)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                                                         STAT=epemalloc_stat(20))
        ASSERT(epemalloc_stat(20).eq.0)
        epemalloc_stat(20)=1
              coul_int_c(i_grad)%l(i_l)%m = zero
           enddo ! i_grad=0,3
        enddo
     endif ! lmax_ch.gt.0

     equal_3_coul_ref: do i_ea3 = 1,n_equal_c
        do_rotation=.false. !!!!
        xc = unique_atoms_eperef(ua3)%position(:,i_ea3)

        ! now do a precalculation of solid harmonics
        if (lmax_ch.gt.0) then

           allocate(yl_arg(num,3),STAT=epemalloc_stat(23))
           ASSERT(epemalloc_stat(23).eq.0)
           epemalloc_stat(23)=1

           yl_arg(:,1) = gamma_arg(:,1) - xc(1)
           yl_arg(:,2) = gamma_arg(:,2) - xc(2)
           yl_arg(:,3) = gamma_arg(:,3) - xc(3)

           yl_arr(:,:,i_ea3) = solid_harmonics_calc(lmax_ch,yl_arg)
           deallocate(yl_arg,STAT=epemalloc_stat(23)) 
           ASSERT(epemalloc_stat(23).eq.0)

           do lm=1,(lmax_ch+1)**2
              do i_sum=1,solhrules_differential(3,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,1) = &
                      yl_arr_grad(:,i_ea3,lm,1) + &
                      solhrules_differential(3,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(3,lm)%lm_sh(i_sum),i_ea3)
              enddo
              do i_sum=1,solhrules_differential(4,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,2) = &
                      yl_arr_grad(:,i_ea3,lm,2) + &
                      solhrules_differential(4,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(4,lm)%lm_sh(i_sum),i_ea3)
              enddo
              do i_sum=1,solhrules_differential(2,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,3) = &
                      yl_arr_grad(:,i_ea3,lm,3) + &
                      solhrules_differential(2,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(2,lm)%lm_sh(i_sum),i_ea3)
              enddo
           enddo ! lm=1,(lmax_ch+1)**2

        endif ! lmax_ch.gt.0

! PROCESS [xi_i|f_k|xi_j] **********************
    ncexps = unique_atoms_eperef(ua3)%l_ch(0)%n_exponents
    cexps => unique_atoms_eperef(ua3)%l_ch(0)%exponents(:)
         call s_coulomb()
    ncexps = unique_atoms_eperef(ua3)%r2_ch%n_exponents
    cexps => unique_atoms_eperef(ua3)%r2_ch%exponents(:)
         call r2_coulomb()

     enddo equal_3_coul_ref

     if (lmax_ch.gt.0) then
        call calc_sym_coef_eperef()
        ! finally calculate l-type coulomb integrals
        call l_coulomb(.true.)
     endif

        do i_grad=0,3
           call fitcontract_eperef&
              ('use_eperef',num,ua3,cutoff,coul_int_c(i_grad)%l,&
                                               prim_epe_gref(i_grad)%m)
           do i_l = -1, lmax_ch
              deallocate(coul_int_c(i_grad)%l(i_l)%m,STAT=epemalloc_stat(20))
              ASSERT(epemalloc_stat(20).eq.0)
           enddo
        do i_l =  1, lmax_ch
           if(dealloc_epe_ref) then
                do independent_fct=1,unique_atoms_eperef(ua3)&
                                %symadapt_partner(1,i_l)%N_independent_fcts
!          deallocate( unique_atoms_eperef(ua3)%&
!        renormaliation_partner_ch(i_l)%renorm(independent_fct)%c_exp,stat=)
                enddo
           endif
        enddo
!!!!!!!!!!!!!!!
           deallocate (coul_int_c(i_grad)%l,STAT=epemalloc_stat(19))
           ASSERT(epemalloc_stat(19).eq.0)
        enddo

        if (lmax_ch.gt.0) then
          deallocate(sym_coef2,sym_coef1,STAT=epemalloc_stat(22))
          ASSERT(epemalloc_stat(22).eq.0) 
          deallocate(yl_arr_grad,yl_arr,STAT=epemalloc_stat(21))
          ASSERT(epemalloc_stat(21).eq.0)
       endif
    deperef: if(dealloc_epe_ref) then 
          if(unique_atoms_eperef(ua3)%r2_ch%N_contracted_fcts.ne.0) then
             deallocate(unique_atoms_eperef(ua3)%r2_ch%contractions,stat=epemalloc_stat(24))
             ASSERT(epemalloc_stat(24).eq.0)
          endif

          deallocate(unique_atoms_eperef(ua3)%r2_ch%exponents,stat=epemalloc_stat(25))
          ASSERT(epemalloc_stat(25).eq.0)

!!$          deallocate(unique_atoms_eperef(ua3)%r2_ch,stat=epemalloc_stat(30))
!!$          ASSERT(epemalloc_stat(30).eq.0)

          deallocate(unique_atoms_eperef(ua3)%l_ch(0)%exponents,stat=epemalloc_stat(27)) 
          ASSERT(epemalloc_stat(27).eq.0)
 
         do i_l = 1,lmax_ch
          deallocate(unique_atoms_eperef(ua3)%l_ch(i_l)%exponents,stat=epemalloc_stat(45)) 
          ASSERT(epemalloc_stat(45).eq.0)
         enddo

          deallocate(unique_atoms_eperef(ua3)%l_ch,stat=epemalloc_stat(32))  !!!!!
          ASSERT(epemalloc_stat(32).eq.0)

          if (lmax_ch.gt.0) then 
             deallocate(unique_atoms_eperef(ua3)%symadapt_partner,stat=epemalloc_stat(42)) 
             ASSERT(epemalloc_stat(42).eq.0)
          endif

          if(unique_atoms_eperef(ua3)%N_glob_cons_ch.gt.0) then
           deallocate(unique_atoms_eperef(ua3)%glob_con_ch,stat=epemalloc_stat(28))
           ASSERT(epemalloc_stat(28).eq.0)
          endif
    endif deperef
 endif uaeperef
           if(dealloc_epe_ref) then
               deallocate(ua_eperef%position,stat=epemalloc_stat(10))
              ASSERT(epemalloc_stat(10).eq.0)
           endif
  enddo unique_3_ref

   if(dealloc_epe_ref) then
!     print*,'assiciated orbitalprojection_ch_eperef',associated(orbitalprojection_ch_eperef)
          if ( associated(orbitalprojection_ch_eperef) ) &
          deallocate(orbitalprojection_ch_eperef,stat=epemalloc_stat(4))
          ASSERT(epemalloc_stat(4).eq.0)
          deallocate(unique_atoms_eperef,coeff_charge_eperef,stat=epemalloc_stat(2))
          ASSERT(epemalloc_stat(2).eq.0)
   endif 
!MMMMMMMMMMMMMMMMMMMMMMM end REFERENCE MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

!WWWWWWWWWWWWWWWWWWWWWWW   begin VAR   WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  deallocate(gamma_help,STAT=epemalloc_stat(18))
          ASSERT(epemalloc_stat(18).eq.0)
  dealloc_coeff_charge=dealloc_epe_ref
  lmax_ch = maxval(unique_atoms(:)%lmax_ch)
  max_order = max(lmax_ch+2,3)
  allocate (gamma_help(num,max_order),STAT=epemalloc_stat(18)) !done
          ASSERT(epemalloc_stat(18).eq.0)
          epemalloc_stat(18)=1

  unique_3: do ua3 = 1,N_unique_atoms+n_timps
     if(ua3.gt.N_unique_atoms) then
        ua=>unique_timps(ua3-N_unique_atoms)
     else
        ua=>unique_atoms(ua3)
     endif

   zc = ua%Z
   zcc= ua%ZC
   n_equal_c = ua%N_equal_atoms

   equal_3: do i_ea3 = 1,n_equal_c
        xc = ua%position(:,i_ea3)
        gamma_arg2 = ( (gamma_arg(:,1)-xc(1))**2 + &
                       (gamma_arg(:,2)-xc(2))**2 + &
                       (gamma_arg(:,3)-xc(3))**2  )
     do i_grad=1,3
        nuc_grad(:,i_grad)=nuc_grad(:,i_grad)+&
        (zc-zcc)*(gamma_arg(:,i_grad)-xc(i_grad))/sqrt(gamma_arg2)/gamma_arg2
     end do
        nuc_v(:)=nuc_v(:)+ (zc-zcc)/sqrt(gamma_arg2)
   enddo equal_3


     ! now fitfunctions  ***************

   if(ua3.le.N_unique_atoms) then
     ! Orbital gradients:
     lmax_ch = int(unique_atoms(ua3)%lmax_ch,kind=i4_kind)
     max_order = max(3,lmax_ch+2)

     do i_grad=0,3     
        allocate(coul_int_c(i_grad)%l(-1:lmax_ch), stat=epemalloc_stat(19)) !done
        ASSERT(epemalloc_stat(19).eq.0)
        epemalloc_stat(19)=1

        ncexps = unique_atoms(ua3)%r2_ch%n_exponents
        allocate(coul_int_c(i_grad)%l(-1)%m&
                           (num,ncexps,1,1,1),stat=epemalloc_stat(20)) !done
        ASSERT(epemalloc_stat(20).eq.0)
        epemalloc_stat(20)=1
        ncexps = unique_atoms(ua3)%l_ch(0)%n_exponents
        allocate(coul_int_c(i_grad)%l(0)%m&
                           (num,ncexps,1,1,1),stat=epemalloc_stat(20)) !done
        ASSERT(epemalloc_stat(20).eq.0)
        epemalloc_stat(20)=1

        coul_int_c(i_grad)%l(0)%m=0.0_r8_kind
        coul_int_c(i_grad)%l(-1)%m=0.0_r8_kind
     enddo ! i_grad=0,3

     if (lmax_ch.gt.0) then

        allocate(yl_arr(num,(lmax_ch+1)**2,n_equal_c),  &
                 yl_arr_grad(num,n_equal_c,(lmax_ch+1)**2,3), &
                                           STAT=epemalloc_stat(21)) !done
        ASSERT(epemalloc_stat(21).eq.0)
        epemalloc_stat(21)=1
        yl_arr = zero
        yl_arr_grad = zero

        n_indep_max = maxval(unique_atoms(ua3)%symadapt_partner&
             (1,:)%n_independent_fcts)
        allocate(sym_coef1(num,n_equal_c,n_indep_max,lmax_ch),&
                 sym_coef2(num,n_equal_c,n_indep_max,lmax_ch,3),&
                                           STAT=epemalloc_stat(22)) !done
        ASSERT(epemalloc_stat(22).eq.0)
        epemalloc_stat(22)=1

        sym_coef1 = zero
        sym_coef2 = zero

        do i_l = 1,lmax_ch
           ncexps = unique_atoms(ua3)%l_ch(i_l)%n_exponents
           n_indep =  &
                unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
           do i_grad=0,3   
              allocate (coul_int_c(i_grad)%l(i_l)%m(num,ncexps,n_indep,1,1),&
                                             STAT=epemalloc_stat(20)) !done
              ASSERT(epemalloc_stat(20).eq.0)
              epemalloc_stat(20)=1
              coul_int_c(i_grad)%l(i_l)%m = zero
           enddo ! i_grad=0,3
        enddo  !  (i_l)?x AG
     endif ! lmax_ch.gt.0

     equal_3_coul: do i_ea3 = 1,n_equal_c
        do_rotation=.false. !!!!
        xc = unique_atoms(ua3)%position(:,i_ea3)

        ! now do a precalculation of solid harmonics
        if (lmax_ch.gt.0) then

           allocate(yl_arg(num,3),STAT=epemalloc_stat(23)) !done
              ASSERT(epemalloc_stat(23).eq.0)
              epemalloc_stat(23)=1

           yl_arg(:,1) = gamma_arg(:,1) - xc(1)
           yl_arg(:,2) = gamma_arg(:,2) - xc(2)
           yl_arg(:,3) = gamma_arg(:,3) - xc(3)

           yl_arr(:,:,i_ea3) = solid_harmonics_calc(lmax_ch,yl_arg)

           deallocate(yl_arg,STAT=epemalloc_stat(23)) 
              ASSERT(epemalloc_stat(23).eq.0)

           do lm=1,(lmax_ch+1)**2
              do i_sum=1,solhrules_differential(3,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,1) = &
                      yl_arr_grad(:,i_ea3,lm,1) + &
                      solhrules_differential(3,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(3,lm)%lm_sh(i_sum),i_ea3)
              enddo
              do i_sum=1,solhrules_differential(4,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,2) = &
                      yl_arr_grad(:,i_ea3,lm,2) + &
                      solhrules_differential(4,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(4,lm)%lm_sh(i_sum),i_ea3)
              enddo
              do i_sum=1,solhrules_differential(2,lm)%n_summands
                 yl_arr_grad(:,i_ea3,lm,3) = &
                      yl_arr_grad(:,i_ea3,lm,3) + &
                      solhrules_differential(2,lm)%coef(i_sum)* &
                      yl_arr(:,solhrules_differential(2,lm)%lm_sh(i_sum),i_ea3)
              enddo
           enddo ! lm=1,(lmax_ch+1)**2

        endif ! lmax_ch.gt.0

    ncexps = unique_atoms(ua3)%l_ch(0)%n_exponents
    cexps => unique_atoms(ua3)%l_ch(0)%exponents(:)
         call s_coulomb()
    ncexps = unique_atoms(ua3)%r2_ch%n_exponents
    cexps => unique_atoms(ua3)%r2_ch%exponents(:)
         call r2_coulomb()
     enddo equal_3_coul

     if (lmax_ch.gt.0) then
        call calc_sym_coef()
        ! finally calculate l-type coulomb integrals
        call l_coulomb(.false.)
     endif

        do i_grad=0,3
           call fitcontract('grad',num,ua3,cutoff,coul_int_c(i_grad)%l,&
                                              prim_epe_grad(i_grad)%m)
!           call fitcontract_eperef &
!            ('default case',num,ua3,cutoff,coul_int_c(i_grad)%l,&
!                                               prim_epe_grad(i_grad)%m)
           do i_l = -1, lmax_ch
              deallocate(coul_int_c(i_grad)%l(i_l)%m,STAT=epemalloc_stat(20))
              ASSERT(epemalloc_stat(20).eq.0)
           enddo

           deallocate (coul_int_c(i_grad)%l,STAT=epemalloc_stat(19))
           ASSERT(epemalloc_stat(19).eq.0)
        enddo

        if (lmax_ch.gt.0) then
           deallocate(sym_coef2,sym_coef1,STAT=epemalloc_stat(22))
           ASSERT(epemalloc_stat(22).eq.0)
           deallocate(yl_arr_grad,yl_arr,STAT=epemalloc_stat(21))
           ASSERT(epemalloc_stat(21).eq.0)
        endif
     end if
  enddo unique_3
!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO End VAR OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

!!$  if( comm_i_am_master() ) then
!!$     print*,'Master: epe_field_and reference:,Del(Vel-Vnuc),Vel,Vnuc'
!!$     print*,'M',prim_epe_grad(0)%m(num,1,1,1)-nuc_v(num) &
!!$                     ,prim_epe_grad(0)%m(num,1,1,1),nuc_v(num), num
! old print out
!!$     print*,'M',prim_epe_grad(1)%m(1,1,1,1)-nuc_grad(1,1)&
!!$              ,'M',prim_epe_grad(2)%m(1,1,1,1)-nuc_grad(1,2)&
!!$              ,'M',prim_epe_grad(3)%m(1,1,1,1)-nuc_grad(1,3)
!!$
!!$     print*,'Master:reference Delta,El,nuc_vref'
!!$     print*,'M',prim_epe_gref(0)%m(num,1,1,1)-nuc_vref(num) , &
!!$            'M',prim_epe_gref(0)%m(num,1,1,1),nuc_vref(num)
!!$
!!$     print*,'M',prim_epe_gref(1)%m(1,1,1,1)-nuc_gref(1,1)&
!!$              ,'M',prim_epe_gref(2)%m(1,1,1,1)-nuc_gref(1,2)&
!!$              ,'M',prim_epe_gref(3)%m(1,1,1,1)-nuc_gref(1,3)
!!$   else
!       print*,'Slave: epe_field_and reference:,Del(Vel-Vnuc),Vel,Vnuc'
!       print*,'S',prim_epe_grad(0)%m(num,1,1,1)-nuc_v(num) &
!                       ,prim_epe_grad(0)%m(num,1,1,1),nuc_v(num), num

!       print*,'S',prim_epe_grad(1)%m(1,1,1,1)-nuc_grad(1,1)&
!             ,'S',prim_epe_grad(2)%m(1,1,1,1)-nuc_grad(1,2)&
!             ,'S',prim_epe_grad(3)%m(1,1,1,1)-nuc_grad(1,3)

!       print*,'Slave:reference Delta,El,nuc_vref'
!       print*,'S',prim_epe_gref(0)%m(num,1,1,1)-nuc_vref(num) , &
!              'S',prim_epe_gref(0)%m(num,1,1,1),nuc_vref(num)

!       print*,'S',prim_epe_gref(1)%m(1,1,1,1)-nuc_gref(1,1)&
!             ,'S',prim_epe_gref(2)%m(1,1,1,1)-nuc_gref(1,2)&
!             ,'S',prim_epe_gref(3)%m(1,1,1,1)-nuc_gref(1,3)
!!$   end if
       if(reg_2a_treated) then
        n_treated_ions= reg_2a_n_ions
        else
        n_treated_ions = reg_I_n_ions
     end if
  if(comm_i_am_master() ) then
       reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vs=zero
       reg_I_pg(1+n_epe_vacancies:n_treated_ions)%vc=zero
     do ix=1,3
       reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gc(ix)=zero
       reg_I_pg(1+n_epe_vacancies:n_treated_ions)%gs(ix)=zero  
     end do
  else ! slave
       reg_I_pg(1:n_treated_ions)%vs=zero
       reg_I_pg(1:n_treated_ions)%vc=zero
     do ix=1,3
       reg_I_pg(1:n_treated_ions)%gc(ix)=zero
       reg_I_pg(1:n_treated_ions)%gs(ix)=zero  
     end do
  end if

   reg_I_pg(g_1:g_2)%vs= &
             prim_epe_grad(0)%m(1:num/2,1,1,1)-nuc_v(1:num/2)-&
            ( prim_epe_gref(0)%m(1:num/2,1,1,1)-nuc_vref(1:num/2))
   reg_I_pg(g_1:g_2)%vc= &
            prim_epe_grad(0)%m(1+num/2:num,1,1,1)-nuc_v(1+num/2:num)-&
           (prim_epe_gref(0)%m(1+num/2:num,1,1,1)-nuc_vref(1+num/2:num))

   do ix=1,3
     reg_I_pg(g_1:g_2)%gs(ix)= &
          prim_epe_grad(ix)%m(1:num/2,1,1,1)-nuc_grad(1:num/2,ix)-&
        ( prim_epe_gref(ix)%m(1:num/2,1,1,1)-nuc_gref(1:num/2,ix) )
     reg_I_pg(g_1:g_2)%gc(ix)= &
         prim_epe_grad(ix)%m(1+num/2:num,1,1,1)-nuc_grad(1+num/2:num,ix)-&
       ( prim_epe_gref(ix)%m(1+num/2:num,1,1,1)-nuc_gref(1+num/2:num,ix) )
   end do

        epg_cluster_reg_I=zero
        epg_nuc=zero
        epg_reference_reg_I=zero

        do i_epe=1,num/2 
        epg_reference_reg_I=epg_reference_reg_I &
               +(prim_epe_gref(0)%m(i_epe,1,1,1)-nuc_vref(i_epe))* &    
               q_shell(epe(g_1-1+i_epe)%k)/qau_qepe &
               +(prim_epe_gref(0)%m(num/2+i_epe,1,1,1)  &
               -nuc_vref(num/2+i_epe))*&
               q_nuclear(epe(g_1-1+i_epe)%k)/qau_qepe

        epg_cluster_reg_I=epg_cluster_reg_I+ &
               (prim_epe_grad(0)%m(i_epe,1,1,1)-nuc_v(i_epe))* &        
               q_shell(epe(g_1-1+i_epe)%k)/qau_qepe &
               +(prim_epe_grad(0)%m(num/2+i_epe,1,1,1)  &
               -nuc_v(num/2+i_epe)) * &
               q_nuclear(epe(g_1-1+i_epe)%k)/qau_qepe

        epg_nuc=epg_nuc+nuc_v(i_epe)*q_shell(epe(g_1-1+i_epe)%k)/qau_qepe+ &
               nuc_v(num/2+i_epe)*q_nuclear(epe(g_1-1+i_epe)%k)/qau_qepe
        enddo
!!$     if(comm_i_am_master()) then
!!$     print*,'Master: epg_cluster_reg_I,epg_nuc',epg_cluster_reg_I,epg_nuc
!!$     else
!        print*,'Slave: epg_cluster_reg_I,epg_nuc',epg_cluster_reg_I,epg_nuc
!!$     end if

!   Sending the results to collect them in subroutine 
!  "epe_collect_gradients()" (slaves only)

 slav: if(.not. comm_i_am_master() ) then
    call comm_init_send(comm_master_host, msgtag_epe_grad_done)
    allocate ( epe_send_array(n_treated_ions,8) ,stat=epemalloc_stat(39))
    ASSERT(epemalloc_stat(39).eq.0)
           epemalloc_stat(39)=1
    allocate ( epe_convergence_energies(3),stat=epemalloc_stat(38) )
    ASSERT(epemalloc_stat(38).eq.0)
           epemalloc_stat(38)=1
    epe_send_array(:,:)=zero
    epe_send_array(1:n_treated_ions,1) = reg_I_pg(1:n_treated_ions)%vs
    epe_send_array(1:n_treated_ions,2) = reg_I_pg(1:n_treated_ions)%gs(1)
    epe_send_array(1:n_treated_ions,3) = reg_I_pg(1:n_treated_ions)%gs(2)
    epe_send_array(1:n_treated_ions,4) = reg_I_pg(1:n_treated_ions)%gs(3)
    epe_send_array(1:n_treated_ions,5) = reg_I_pg(1:n_treated_ions)%vc
    epe_send_array(1:n_treated_ions,6) = reg_I_pg(1:n_treated_ions)%gc(1)
    epe_send_array(1:n_treated_ions,7) = reg_I_pg(1:n_treated_ions)%gc(2)
    epe_send_array(1:n_treated_ions,8) = reg_I_pg(1:n_treated_ions)%gc(3)

    epe_convergence_energies(1) = epg_reference_reg_I
    epe_convergence_energies(2) = epg_cluster_reg_I
    epe_convergence_energies(3) = epg_nuc
    call commpack(epe_send_array(1,1),8*n_treated_ions,1,status)
          if(status.ne.0) call error_handler &
                         ("epe_field_and_forces :  error [1]")
    call commpack(epe_convergence_energies(1),3,1,status)
          if(status.ne.0) call error_handler &
                         ("epe_field_and_forces :  error [2]")
    call comm_send()
    deallocate( epe_send_array,stat=epemalloc_stat(39) )
    ASSERT(epemalloc_stat(39).eq.0)

    deallocate( epe_convergence_energies,stat=epemalloc_stat(38) )
    ASSERT(epemalloc_stat(38).eq.0)

!!$     print*,'epe_field_and_forces: dealloc_coeff_charge', dealloc_coeff_charge
     if(dealloc_coeff_charge) then
        ! FIXME: not sure if shutdown call is ok for deallocating coeff_charge
        call fit_coeff_shutdown()
!       deallocate(coeff_charge,stat=epemalloc_stat(37))
!       ASSERT(epemalloc_stat(37).eq.0)

        deallocate(reg_I_pg,stat=epemalloc_stat(1))
           if (epemalloc_stat(1) /= 0) call error_handler &
                ("epe_field_and_forces_par: deallocation of reg_I_pg failed")
     end if
 endif slav
        call gamma_close() 
        deallocate (gamma_help,STAT=epemalloc_stat(18))  
        ASSERT(epemalloc_stat(18).eq.0)
        deallocate(coul_int_c,STAT=epemalloc_stat(17))
        ASSERT(epemalloc_stat(17).eq.0)
        deallocate(prim_epe_grad(0)%m, &
                   prim_epe_grad(1)%m, &
                   prim_epe_grad(2)%m, &
                   prim_epe_grad(3)%m,STAT=epemalloc_stat(13))
        ASSERT(epemalloc_stat(13).eq.0)
        deallocate(prim_epe_gref(0)%m, &
                   prim_epe_gref(1)%m, &
                   prim_epe_gref(2)%m, &
                   prim_epe_gref(3)%m,STAT=epemalloc_stat(14))
        ASSERT(epemalloc_stat(14).eq.0)
        deallocate(prim_epe_grad,prim_epe_gref,stat=epemalloc_stat(12))
        ASSERT(epemalloc_stat(12).eq.0)
        deallocate(gamma_arg,gamma_arg2,help_arr,stat=epemalloc_stat(15))
        ASSERT(epemalloc_stat(15).eq.0)
        deallocate(nuc_grad,nuc_v,nuc_gref,nuc_vref,stat=epemalloc_stat(16))
        ASSERT(epemalloc_stat(16).eq.0)
        deallocate(cutoff,stat=epemalloc_stat(11))
        ASSERT(epemalloc_stat(11).eq.0)

contains
!       subroutine define_do_rotation
!        if(grad_dim==3) then
!           if(sum((unique_atom_grad_info(imc)%m(:,:,i_ea3)-unity_matrix)**2)<&
!                                                          1.0e-7_r8_kind) then
!              do_rotation=.false.
!           else
!              do_rotation=.true.
!              rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
!           endif
!        else
!           do_rotation=.true.
!           rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
!        endif
!       end subroutine define_do_rotation

  !**************************************************************

  subroutine calc_sym_coef_eperef()
    ! Purpose: piece of the code
    !          calculate symmetry coefficients for the l-type 3 center
    !          fitintegrals

    ang_momentum_symadapt: do i_l=1,lmax_ch
       n_indep_fcts =  &
          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       do i_ind=1,n_indep_fcts
          n_contributing_fcts = &
               unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%N_fcts
          magn => unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          coef =>  unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          eq_atom => unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%I_equal_atom

          do i_cont=1,n_contributing_fcts
             lm = i_l**2 + magn(i_cont)
             sym_coef1(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coef1(:,eq_atom(i_cont),i_ind,i_l) + &
                  yl_arr(:,lm,eq_atom(i_cont))*&
                  coef(i_cont)
             do i = 1,3
                sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) = &
                     sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) + &
                     yl_arr_grad(:,eq_atom(i_cont),lm,i)*coef(i_cont)
             enddo
          enddo   ! i_cont=1,n_contributing_fcts

     if(dealloc_epe_ref) then
      deallocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m&
       ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c&
       ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%I_equal_atom, &
        stat=epemalloc_stat(36))
       ASSERT(epemalloc_stat(36).eq.0)
      endif   ! dealloc_epe_ref
       enddo
       if(dealloc_epe_ref) then
        deallocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt, &
                                                         stat=epemalloc_stat(35)) ! %symadapt
        ASSERT(epemalloc_stat(35).eq.0)
       endif
    enddo ang_momentum_symadapt
  end subroutine calc_sym_coef_eperef
  !**************************************************************
  subroutine calc_sym_coef()
    ! Purpose: piece of the code
    !          calculate symmetry coefficients for the l-type 3 center
    !          fitintegrals
!  integer(kind=i4_kind)          :: N_partner

    ang_momentum_symadapt: do i_l=1,lmax_ch
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
 
       do i_ind=1,n_indep_fcts
          n_contributing_fcts = &
               unique_atoms(ua3)%symadapt_partner(1,i_l)%&
                                  symadapt(i_ind,1)%N_fcts
          magn => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          coef =>  unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          eq_atom => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%I_equal_atom
          do i_cont=1,n_contributing_fcts
             lm = i_l**2 + magn(i_cont)
             sym_coef1(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coef1(:,eq_atom(i_cont),i_ind,i_l) + &
                  yl_arr(:,lm,eq_atom(i_cont))*&
                  coef(i_cont)
             do i = 1,3
                sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) = &
                     sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) + &
                     yl_arr_grad(:,eq_atom(i_cont),lm,i)*coef(i_cont)
             enddo
          enddo
       enddo
    enddo ang_momentum_symadapt
  end subroutine calc_sym_coef


  subroutine l_coulomb(epe_reference)
    ! Purpose: claculate gradients of l-type fit integrals
    real(kind=r8_kind), dimension(num) :: help_vec1, help_vec2, &
         help_vec3, help_vec4,  help_vec7
         logical::epe_reference

    ang_momentum: do i_l=1,lmax_ch
!       n_indep_fcts =  &
!            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
        if(epe_reference) then
       n_indep_fcts =  &
            unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms_eperef(ua3)%l_ch(i_l)%N_exponents
       cexps => unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(:)
        else
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms(ua3)%l_ch(i_l)%N_exponents
       cexps => unique_atoms(ua3)%l_ch(i_l)%exponents(:)
        endif ! get_epe_reference
       do i_ind=1,n_indep_fcts
          equal_3_l: do i_ea3=1,n_equal_c
!             if(grad_dim==3) then
!                if(sum((unique_atom_grad_info(imc)%m(:,:,i_ea3)&
!                              -unity_matrix)**2)< 1.0e-7_r8_kind) then
!                   do_rotation=.false.
!                else
!                   do_rotation=.true.
!                   rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
!                endif
!             else
!                do_rotation=.true.
!                rotmat=>unique_atom_grad_info(imc)%m(:,:,i_ea3)
!             end if
        do_rotation=.false.
        if(epe_reference) then
             xc =  unique_atoms_eperef(ua3)%position(:,i_ea3)
        else
             xc =  unique_atoms(ua3)%position(:,i_ea3)
        endif
             do k = 1,ncexps
                gamma_arg2 = cexps(k)* &
                     (  (gamma_arg(:,1) - xc(1))**2 + &
                        (gamma_arg(:,2) - xc(2))**2 + &
                        (gamma_arg(:,3) - xc(3))**2) 
                gamma_help(:,1:lmax_ch+2) = gamma(lmax_ch+2,gamma_arg2)
                help_vec1=two*pi/cexps(k)*(two*cexps(k))**i_l
                help_vec2=gamma_help(:,i_l+1)
                help_vec3=two*cexps(k)* &
                     gamma_help(:,i_l+2)*sym_coef1(:,i_ea3,i_ind,i_l)
                do i = 1,3
                   help_vec4=help_vec3*(gamma_arg(:,i)-xc(i))
                   help_vec7=help_vec2*sym_coef2(:,i_ea3,i_ind,i_l,i)
                   help_arr(:,i) =  help_vec1* ( -help_vec7 +  help_vec4)
                enddo
        coul_int_c(0)%l(i_l)%m(:,k,i_ind,1,1)=&
                coul_int_c(0)%l(i_l)%m(:,k,i_ind,1,1) &
                  +help_vec1*sym_coef1(:,i_ea3,i_ind,i_l)*gamma_help(:,i_l+1)
                if(do_rotation) then
                   ! make gradient totalsymmetric before adding
                   do i_grad=1,grad_dim ! only if moving_c
                      coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)=&
                           coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)+&
                           rotmat(i_grad,1)*help_arr(:,1)+&
                           rotmat(i_grad,2)*help_arr(:,2)+&
                           rotmat(i_grad,3)*help_arr(:,3)
                   enddo
                else
                   coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,1)
                   coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,2)
                   coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,3)
                end if
             enddo! loop over k
          enddo equal_3_l
       enddo
    enddo ang_momentum
  end subroutine l_coulomb
  !**************************************************************

  !**************************************************************
  subroutine r2_coulomb()
    ! Purpose: calculate gradients of r2 coulomb integrals
    !---** r2-type Fitfct. **---
    real(kind=r8_kind), dimension(num) :: help_vec1,  help_vec3, &
         help_vec4

    do k=1,ncexps

       gamma_arg2 = cexps(k)* &
            (  (gamma_arg(:,1) - xc(1))**2 + &
               (gamma_arg(:,2) - xc(2))**2 + &
               (gamma_arg(:,3) - xc(3))**2)
       gamma_help(:,1:3) = gamma(3,gamma_arg2)

       help_vec1=two*pi/(cexps(k)**2)
        coul_int_c(0)%l(-1)%m(:,k,1,1,1)=coul_int_c(0)%l(-1)%m(:,k,1,1,1)&
                +help_vec1*( gamma_help(:,1)+gamma_arg2*gamma_help(:,2) )
       help_vec3= two*((-gamma_arg2*gamma_help(:,3))*cexps(k))
       help_vec4=(gamma_help(:,1)+gamma_arg2*gamma_help(:,2))
          help_vec1=help_vec1*help_vec3
          do i = 1,3
             help_arr(:,i) = -help_vec1*(gamma_arg(:,i)-xc(i))
          enddo

       if(do_rotation) then
          ! make gradient totalsymmetric before adding
          do i_grad=1,grad_dim ! only if moving_c
             coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)=&
                  coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)+&
                  rotmat(i_grad,1)*help_arr(:,1)+&
                  rotmat(i_grad,2)*help_arr(:,2)+&
                  rotmat(i_grad,3)*help_arr(:,3)
          enddo
       else
          coul_int_c(1)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(1)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,1)
          coul_int_c(2)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(2)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,2)
          coul_int_c(3)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(3)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,3)
       end if
    enddo! r2-exponents, third center
  end subroutine r2_coulomb
  !**************************************************************

  !**************************************************************
  subroutine s_coulomb
    ! Purpose: calculate gradients of s type coulomb fitintegrals
    ! loop over exponents of third center ---** s-type Fitfct. **---
    real(kind=r8_kind), dimension(num) :: help_vec1

    do k=1,ncexps
       gamma_arg2 = cexps(k)* &
           ((gamma_arg(:,1) - xc(1))**2 + &
            (gamma_arg(:,2) - xc(2))**2 + &
            (gamma_arg(:,3) - xc(3))**2)

       gamma_help(:,1:2) = gamma(2,gamma_arg2)
       help_vec1=four*pi*gamma_help(:,2)
        coul_int_c(0)%l(0)%m(:,k,1,1,1)=coul_int_c(0)%l(0)%m(:,k,1,1,1)&
                        +two*pi/cexps(k)*gamma_help(:,1)
          do i = 1,3
             help_arr(:,i) =  help_vec1*(gamma_arg(:,i)-xc(i))
          enddo
       if(do_rotation) then
          ! make gradient totalsymmetric before adding
          do i_grad=1,grad_dim ! only if moving_c
             coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)=&
                  coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)+&
                  rotmat(i_grad,1)*help_arr(:,1)+&
                  rotmat(i_grad,2)*help_arr(:,2)+&
                  rotmat(i_grad,3)*help_arr(:,3)
          enddo
       else
          coul_int_c(1)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(1)%l(0)%m(:,k,1,1,1)+help_arr(:,1)
          coul_int_c(2)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(2)%l(0)%m(:,k,1,1,1)+help_arr(:,2)
          coul_int_c(3)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(3)%l(0)%m(:,k,1,1,1)+help_arr(:,3)
       endif
    enddo! s-exponents of third center 
  end subroutine s_coulomb
  !**************************************************************
  end subroutine epe_field_and_forces_par

  subroutine print_epemod_alloc()
  integer(kind=i4_kind):: i
  do i=1,size(epemalloc_stat)
   if(epemalloc_stat(i).ne.0) print*,i, 'epemod_alloc ne 0'
  enddo
  end subroutine print_epemod_alloc
  
  
end module epe_module

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
module  gradient_data_module
  !---------------------------------------------------------------
  !
  !  Purpose: This module contains the data for the
  !           nuclear displacement gradients, and
  !           the corresponding subroutines for allocation, sending, ...
  !
  !  Module called by: main_gradient
  !
  !
  !
  !
  !
  !  Author: MS
  !  Date: 12/96
  !
  ! Modification (Please copy before editing)
  ! Author: Uwe Birkenheuer
  ! Date:   June 1998
  ! Description: Moving_Unique_Atom concept introduced
  !              Split_Gradients concept introduced
  !              Gradients for Model_Density_Approach introduced
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
  !== Interrupt end of public interface of module =================
#ifdef _LINUX
# ifdef _UNDERSCORE
#  define GETENV getenv_
# else
#  define GETENV getenv
# endif
#endif

  !------------ all modules used are public !! ------------------
# include "def.h"
  use type_module
  use unique_atom_module
  use datatype
  use iounitadmin_module
  use comm_module
  use msgtag_module
  use symmetry_data_module
  use operations_module
  use options_module, only: options_relativistic, options_split_gradients, &
                            options_xcmode, xcmode_model_density, &
                            options_n_spin
  use output_module
  use dipole_module
  use solv_cavity_module, only: grad_solv_totsym,grad_solv_totsym_tes, &
       grad_solv_cart,dealloc_geom_grad
  use solv_electrostat_module, only: init_forces_on_pc,Q_grad
  use solv_cavity_module, only: with_pc,fixed_pc,to_calc_grads
  use potential_module, only: N_points
  use elec_static_field_module, only: totsym_field_length,totalsym_field
  use point_dqo_module
  use induced_dipoles_module
#ifdef WITH_EFP
  use efp_rep_module, only: init_repf_grads, repf_grads_shutdown
#endif
#ifdef WITH_EPE
  use ewaldpc_module,only: ewpc_n
#endif
#ifdef FPP_DEBUG
    use error_module, only: MyID
#endif
    USE_MEMLOG

  !== Interrupt of public interface of module =====================
  implicit none
  save ! save all variables defined in this module
  private

  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ----------------------------



  !------------ Declaration of constants and variables ----------

  integer(kind=i4_kind),public :: gradient_data_n_gradients
  ! total number of symmetry adapted Etot gradient components
  integer(kind=i4_kind),public :: gradient_data_n_spin_gradients
  ! total number of symmetry adapted Etot gradient components times the
  ! spin-degeneracy required for evaluating the 3-center contributions
  integer(kind=i4_kind),public,allocatable :: gradient_rel_index(:,:)
  ! (gradient_data_n_gradients,n_irreps)
  ! contains a processor index p = gradient_rel_index(x,irr)
  ! that will transform x-gradient of irrep irr
  integer(kind=i4_kind),public,allocatable :: gradient_der_index(:,:,:)
  ! (gradient_data_n_gradients,gradient_data_n_gradients,n_irreps)
  ! contains a processor index p = gradient_der_index(x,y,irr)
  ! that will transform (x,y)-derivative of irrep irr

  integer(kind=i4_kind),public :: gradient_data_rel_blocklength
  ! blocklength used for quadrupel files

  real(kind=r8_kind),public,allocatable,target :: dervs_totalsym(:,:)
  real(kind=r8_kind),public,allocatable :: cpks_grad_fit_totasym(:,:)
  real(kind=r8_kind),public,allocatable,target :: cpks_gradient_totalsym(:,:)
  real(kind=r8_kind),public,allocatable,target :: cpks_grad_totalsym(:,:)

  real(kind=r8_kind),public,allocatable,target :: gradient_totalsym(:)
  ! gradient_totalsym(gradient_data_n_gradients) : stores all
  ! totalsymmetric first energy gradients
  ! (not including the contributions from the numerical grid, UB)

  real(kind=r8_kind),public,allocatable :: gradient_ob_pulay(:)
  ! gradient_ob_pulay(gradient_data_n_gradients)
  ! Pulay correction to e_tot derivatives due to finite orbital basis sets

  logical,public:: calc_cluster_epe_energy=.true.

  real(kind=r8_kind),public,allocatable :: cluster_epe_energy(:)
  ! cluster_epe_energy(1)
  ! Energy of interaction of cluster and epe centers

  real(kind=r8_kind),public,allocatable :: gradient_ch_pulay(:)
  ! gradient_ch_pulay(gradient_data_n_gradients)
  ! Pulay correction to e_tot derivatives due to finite charge fit basis sets
  real(kind=r8_kind),public,allocatable :: gradient_mda_rho(:)
  ! gradient_mda_rho(gradient_data_n_gradients)
  ! MDA Pulay correction to e_tot derivatives due to Sum(k) a_k d/dR f_k
  real(kind=r8_kind),public,allocatable :: gradient_mda_vxc(:)
  ! gradient_mda_vxc(gradient_data_n_gradients)
  ! MDA Pulay correction to e_tot derivatives due to Sum(k) b_k d/dR V_H[f_k]
  type(arrmat4),public,allocatable :: dervs_cartesian(:,:)
  type(arrmat3),public,allocatable :: cpks_gradient_cartesian(:)
  type(arrmat2),public,allocatable :: gradient_cartesian(:)
  ! gradient_cartesian(N_moving_unique_atoms
  !                           +N_moving_unique_timps)%m(3,ua%N_equal_atom)
  ! the cartesian coordinates of the gradients of each non-fixed atom
  type(arrmat2),public,allocatable :: partial_cartesian(:)
  ! partial_cartesian(N_moving_unique_atoms)%m(3,ua%N_equal_atom)
  ! cartesian coordinate contributions to the gradients on each non-fixed atom

  integer(kind=i4_kind),public,allocatable,target :: gradient_index(:)
  ! gradient_index(N_moving_unique_atoms +
  !                               N_moving_unique_timps +  1 )
  ! 1:N_moving_unique_atoms :
  ! pointer on the first symmetry adapted gradient of each non-fixed
  ! unique atom within the array of linearly packed symmetry adapted
  ! gradients
  ! N_moving_unique_atoms + 1 :
  ! pointer  a f t e r  the last symmetry adapted gradient
  real(kind=r8_kind),public,allocatable   :: dervs_fit_ch_cartesian(:,:,:,:)
  real(kind=r8_kind),public,allocatable   :: grad_fit_ch_cartesian(:,:)
  real(kind=r8_kind),public,allocatable   :: cpks_grad_fit_ch_cartesian(:,:,:)
  ! grad_fit_ch_cartesian(N_moving_unique_atoms,3)
  ! The cartesian coordinates of the gradients of the first equal atom of
  ! each non-fixed unique atom as arising from the charge fit functions.
  real(kind=r8_kind),public,allocatable   :: grad_mda_rhofit_cartesian(:,:)
  ! grad_mda_rhofit_cartesian(N_moving_unique_atoms,3)
  ! The cartesian coordinates of the gradients of the first equal atom of
  ! each non-fixed unique atom as arising from the derivative d/dR rho_fit
  ! = Sum(k) a_k d/dR f_k(r) entering the MDA XC force contributions.
  real(kind=r8_kind),public,allocatable   :: grad_mda_xcpot_cartesian(:,:)
  ! grad_mda_xcpot_cartesian(N_moving_unique_atoms,3)
  ! The cartesian coordinates of the gradients of the first equal atom of
  ! each non-fixed unique atom as arising from the derivative d/dR Vxc,mda
  ! = Sum(k) b_k V_H[d/dR f_k] enerting the MDA XC force contributions.

  !------------ public functions and subroutines ----------------
  public gradient_data_setup,gradient_data_shutdown,&
       gradient_data_write_cartesians,&
       gradient_data_write_gxfile

  public :: gradient_sndrcv_fit_ch
  public :: gradient_sndrcv_3c, send_receive_Q_grads

#ifdef WITH_MOLMECH
  public :: qm_grads_to_qmmm
#endif
  public cpks_add_fitmat_ch_grads

  public :: gradient_number

  public :: optimizer_write_cart_hess

#ifdef WITH_SECDER
  public :: gradient_data_write_cart_hess
#endif



  !================================================================
  ! End of public interface of module
  !================================================================



  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains

  function gradient_number(UA) result(num)
    implicit none
    integer(i4_kind), intent(in) :: UA
    integer(i4_kind)             :: num ! result
    ! *** end of interface ***

    ASSERT(allocated(gradient_index))
    ASSERT(UA>0)
    ASSERT(UA+1<=size(gradient_index))

    num = gradient_index(UA+1) - gradient_index(UA)
  end function gradient_number

  subroutine gradient_data_setup()
    !  Purpose: allocation of necesarry arrays
    !
    !  Called by: main_gradient
    !** End of interface ****************************************
    use pointcharge_module
!   use calc3c_switches
    use integralpar_module, only: integralpar_2dervs,integralpar_cpksdervs
    use cpksdervs_matrices
    use occupied_levels_module
    use virtual_levels_module
    use occupation_module, only: occ_num,alloc_occ_num
    use eigen_data_module,only:eigvec,eigval
    use fit_coeff_module, only: get_fit,fit
#ifdef WITH_EPE
    use epecom_module, only: epealloc_stat
#endif
#ifdef WITH_EFP
    use efp_solv_grad_module, only: init_X_solv_grads
#endif
!    use epe_module, only:epemalloc_stat
    implicit none
    !------------ Declaration of local variables ----------------
    type(unique_atom_type), pointer :: ua
    integer(kind=i4_kind)           :: counter,alloc_stat,i_unique,ts
    integer(kind=i4_kind)           :: i_unique2
    logical                         :: split_gradients, model_density
    integer(kind=i4_kind)           :: i_ir,i_grad
!    integer(kind=i4_kind)           :: n_occmo
    integer(kind=i4_kind)           :: n_occs,n_virs
    type(fit)                       :: n_fit
    integer(kind=i4_kind)           :: n_holes,i_occ,i_vir
    integer(kind=i4_kind)           :: i_spn!,i
    real(kind=r8_kind)              :: dgen
    integer(kind=i4_kind)           :: eig_dim
    integer(kind=i4_kind)           :: info

    !------------ Executable code -------------------------------
    split_gradients = options_split_gradients()
    model_density = options_xcmode() == xcmode_model_density

    allocate(gradient_index(N_moving_unique_atoms+N_moving_unique_timps+1), &
         gradient_cartesian(N_moving_unique_atoms+N_moving_unique_timps), &
         stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)

    if(integralpar_2dervs) then
     allocate(dervs_cartesian(N_moving_unique_atoms,N_moving_unique_atoms), &
              stat=cpksalloc(90))
     ASSERT(cpksalloc(90).eq.0)
     MEMLOG(size(dervs_cartesian))
    endif

    if(integralpar_cpksdervs) then
     call get_fit(n_fit)
     ! to be calculated when cpks equation are solved
     allocate(cpks_gradient_cartesian(size(gradient_cartesian)),&
                                             stat=cpksalloc(28))
     ASSERT(cpksalloc(28).eq.0)
    endif

    if (split_gradients) then
       allocate(partial_cartesian(N_moving_unique_atoms+&
       N_moving_unique_timps),stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler( &
            "gradient_data_setup: allocate of gradient_index (s) failed" )
    endif

    if(operations_solvation_effect) then
       allocate(grad_solv_cart(N_moving_unique_atoms),stat=alloc_stat) !!!!!!!!!
       if ( alloc_stat/=0 ) call error_handler( &                      !!!!!!!!!
            "gradient_data_setup: allocate of grad_solv_cart  failed" )!!!!!!!!!
    endif

    ts = get_totalsymmetric_irrep()
    counter = 1
    uniques: do i_unique=1,N_moving_unique_atoms
       ua => unique_atoms(moving_unique_atom_index(i_unique))

       allocate(gradient_cartesian(i_unique)%m(3,ua%N_equal_atoms),&
            stat=alloc_stat)
            ASSERT(alloc_stat.eq.0)
            gradient_cartesian(i_unique)%m=0.0_r8_kind

       if(integralpar_2dervs) then
        do  i_unique2=1,N_moving_unique_atoms
         allocate(dervs_cartesian(i_unique,i_unique2)%m(3,ua%N_equal_atoms,3, &
            unique_atoms(moving_unique_atom_index(i_unique2))%N_equal_atoms), &
            stat=cpksalloc(91) )
          ASSERT(cpksalloc(91).eq.0)
          dervs_cartesian(i_unique,i_unique2)%m=0.0_r8_kind
        enddo
       endif

       if(integralpar_cpksdervs) then
        allocate(cpks_gradient_cartesian(i_unique)%m(n_fit%n_ch,3,ua%N_equal_atoms),&
            stat=cpksalloc(29))
            ASSERT(cpksalloc(29).eq.0)
            cpks_gradient_cartesian(i_unique)%m=0.0_r8_kind
        endif

       if(operations_solvation_effect) then
          allocate(grad_solv_cart(i_unique)%m(3,ua%N_equal_atoms),&         !!!!!!!!!!!
               stat=alloc_stat)                                             !!!!!!!!!!!
          if ( alloc_stat/=0 ) call error_handler( &                        !!!!!!!!!!!
               "gradient_data_setup: allocate of grad_solv_cart%m failed" ) !!!!!!!!!!!
          grad_solv_cart(i_unique)%m=0.0_r8_kind                            !!!!!!!!!!!
       endif

       if (split_gradients) then
          allocate(partial_cartesian(i_unique)%m(3,ua%N_equal_atoms),&
               stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "gradient_data_setup: allocate of partial_cartesian failed" )
          partial_cartesian(i_unique)%m=0.0_r8_kind
       endif
       gradient_index(i_unique)=counter
       counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
    enddo uniques

    gradient_index(N_moving_unique_atoms+1)=counter

    timps: if(n_moving_unique_timps.ne.0_i4_kind) then
    do i_unique=1,N_moving_unique_timps
       ua => unique_timps(moving_unique_timp_index(i_unique))

       allocate(gradient_cartesian(&
            i_unique+N_moving_unique_atoms)%m(3,ua%N_equal_atoms),&
            stat=alloc_stat)
       if ( alloc_stat/=0 ) call error_handler( &
            "gradient_data_setup: allocate of gradient_cartesian failed" )
       gradient_cartesian(i_unique+N_moving_unique_atoms)%m=0.0_r8_kind
       if (split_gradients) then
          allocate(partial_cartesian(&
               i_unique+N_moving_unique_atoms)%m(3,ua%N_equal_atoms),&
               stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "gradient_data_setup: allocate of partial_cartesian failed" )
          partial_cartesian(i_unique+&
               N_moving_unique_atoms)%m=0.0_r8_kind
       endif
       gradient_index(i_unique+N_moving_unique_atoms)=counter
       counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
    enddo
    gradient_index(N_moving_unique_atoms+&
         N_moving_unique_timps+1)=counter
    endif timps

    gradient_data_n_gradients=counter-1
    if (model_density) then
       gradient_data_n_spin_gradients = gradient_data_n_gradients * &
                                        options_n_spin()
    else
       gradient_data_n_spin_gradients = gradient_data_n_gradients
    end if

#ifdef WITH_EPE
    if(ewpc_n.ne.0) then
        allocate(cluster_epe_energy(1),stat=epealloc_stat(26))
        ASSERT(epealloc_stat(26).eq.0)
        epealloc_stat(26)=1
        cluster_epe_energy=0.0_r8_kind
    endif
#endif

    allocate(gradient_totalsym(gradient_data_n_gradients),stat=alloc_stat)
     ASSERT(alloc_stat.eq.0)
    gradient_totalsym=0.0_r8_kind

    if(integralpar_2dervs) then
     allocate(dervs_totalsym(gradient_data_n_gradients,gradient_data_n_gradients),&
                             stat=cpksalloc(39))
     ASSERT(cpksalloc(39).eq.0)
     MEMLOG(size(dervs_totalsym))
     dervs_totalsym=0.0_r8_kind
    endif

    cpks1: if(integralpar_cpksdervs) then

     allocate(cpks(gradient_data_n_gradients,symmetry_data_n_irreps(), &
                   options_n_spin()), stat=cpksalloc(1))
     ASSERT(cpksalloc(1).eq.0)
     allocate(cpks3c(symmetry_data_n_irreps(),options_n_spin()), &
              stat=cpksalloc(7))
      ASSERT(cpksalloc(7).eq.0)
      MEMLOG(size(cpks)+size(cpks3c))


     allocate(cpks_gradient_totalsym(n_fit%n_ch,gradient_data_n_gradients),&
              cpks_grad_totalsym(n_fit%n_ch,gradient_data_n_gradients),&
              stat=cpksalloc(27))
     ASSERT(cpksalloc(27).eq.0)
     MEMLOG(size(cpks_gradient_totalsym))
     MEMLOG(size(cpks_grad_totalsym))
     cpks_gradient_totalsym=0.0_r8_kind
     cpks_grad_totalsym=0.0_r8_kind

     allocate(cpks_grad_fit_totasym(n_fit%n_ch,gradient_data_n_gradients), &
              stat=cpksalloc(164))
     ASSERT(cpksalloc(164).eq.0)
     cpksalloc(164)=0  ! cpks_grad_fit_totasym

     cpks_grad_fit_totasym=0
     MEMLOG(size(cpks_grad_fit_totasym))
     cpks_grad_fit_totasym=0.0_r8_kind

    !!! what are dimentions of up and down cpks equation systems

    !!! first dimension  = occ + holes
    !!! second dimension = vir + holes

   if(comm_parallel()) then
    if(.not.comm_i_am_master())  then
    call alloc_occ_num()
    call comm_save_recv(comm_master_host,msgtag_packed_message)
    do i_ir=1,symmetry_data_n_irreps()
       call communpack(occ_num(i_ir)%m(1,1),size(occ_num(i_ir)%m),1,info)
     ASSERT(info.eq.0)
    enddo

    else
    call comm_init_send(comm_all_other_hosts,msgtag_packed_message)
    do i_ir=1,symmetry_data_n_irreps()
      call commpack(occ_num(i_ir)%m(1,1),size(occ_num(i_ir)%m),1,info)
     ASSERT(info.eq.0)
    enddo
    call comm_send
    endif
   endif

cpks_size=0
irr_cpks: do i_ir=1,symmetry_data_n_irreps()

     eig_dim=size(eigvec(i_ir)%m,1)
     dgen=symmetry_data_n_partners(i_ir)*(2/options_n_spin())

  spn: do i_spn=1,options_n_spin()
       n_holes=0
       n_occs=0
       n_virs=0
       DPRINT  MyID//'occ_num print',associated(occ_num(i_ir)%m)

      do i_occ=1,size(occ_num(i_ir)%m,1) !!!n_occmo
       DPRINT occ_num(i_ir)%m(i_occ,i_spn),dgen
       if(occ_num(i_ir)%m(i_occ,i_spn).lt.0.00001_r8_kind) then
        n_virs=n_virs+1
       elseif(dgen-occ_num(i_ir)%m(i_occ,i_spn).gt.0.00001_r8_kind) then
        n_holes=n_holes+1
       else
        n_occs=n_occs+1
        endif
      enddo

      DPRINT  MyID//' spin i_ir',i_spn,i_ir
      DPRINT  MyID//' n_occs n_holes n_virs ',n_occs,n_holes,n_virs

     cpks3c(i_ir,i_spn)%n_holes=n_holes
     allocate(cpks3c(i_ir,i_spn)%eigvec_vir_and_holes(eig_dim,n_virs+n_holes), &
              cpks3c(i_ir,i_spn)%eigval_vir_and_holes(n_virs+n_holes), &
              cpks3c(i_ir,i_spn)%parocc_vir_and_holes(n_virs+n_holes), &
              cpks3c(i_ir,i_spn)%parocc_occ(n_occs+n_holes), &
              stat=cpksalloc(135))
     MEMLOG(size(cpks3c(i_ir,i_spn)%eigvec_vir_and_holes))
     MEMLOG(size(cpks3c(i_ir,i_spn)%eigval_vir_and_holes))
     MEMLOG(size(cpks3c(i_ir,i_spn)%parocc_vir_and_holes+size(cpks3c(i_ir,i_spn)%parocc_occ)))
     ASSERT(cpksalloc(135).eq.0)

 if(n_virs.ge.1) then
  DPRINT  'dgen for iir',i_ir,dgen
     cpks3c(i_ir,i_spn)%eigvec_vir_and_holes(:,n_holes+1:n_virs+n_holes)= &
                        eigvec(i_ir)%m(:,1+eig_dim-n_virs:,i_spn)
     cpks3c(i_ir,i_spn)%eigval_vir_and_holes(n_holes+1:n_virs+n_holes)= &
                        eigval(i_ir)%m(1+eig_dim-n_virs:,i_spn)
     cpks3c(i_ir,i_spn)%parocc_vir_and_holes(n_holes+1:n_virs+n_holes)= &
                         (dgen-occ_num(i_ir)%m(1+eig_dim-n_virs:,i_spn))/dgen
  endif

 if(n_holes.ne.0)  then

     cpks3c(i_ir,i_spn)%eigvec_vir_and_holes(:,:n_holes)=eigvec(i_ir)%m(:,n_occs+1:n_occs+n_holes,i_spn)
     cpks3c(i_ir,i_spn)%eigval_vir_and_holes(:n_holes)=eigval(i_ir)%m(n_occs+1:n_occs+n_holes,i_spn)

     do  i_occ=1,n_holes
      dgen=symmetry_data_n_partners(i_ir)*(2/options_n_spin())
      cpks3c(i_ir,i_spn)%parocc_vir_and_holes(i_occ)= &
               (dgen-occ_num(i_ir)%m(n_occs+i_occ,i_spn))/dgen
      cpks3c(i_ir,i_spn)%parocc_occ(n_occs+i_occ)= &
               (occ_num(i_ir)%m(n_occs+i_occ,i_spn))/dgen

!    DPRINT  MyID//'i_ir i_hole parocc_occ parocc_vir_and_holes', i_ir,i_occ, &
!             cpks3c(i_ir,i_spn)%parocc_occ(n_occs+i_occ), &
!             cpks3c(i_ir,i_spn)%parocc_vir_and_holes(i_occ)
     enddo
  endif

     cpks3c(i_ir,i_spn)%parocc_occ(:n_occs)=1.0_r8_kind
     cpks3c(i_ir,i_spn)%n_ai=(n_occs+n_holes)*(n_virs+n_holes)

   do i_grad=1,gradient_data_n_gradients
      cpks(i_grad,i_ir,i_spn)%occ_dim=n_occs+n_holes
      allocate(cpks(i_grad,i_ir,i_spn)%s1(n_occs+n_holes,n_occs+n_holes), &
               cpks(i_grad,i_ir,i_spn)%h1(n_occs+n_holes,n_occs+n_holes), &
               cpks(i_grad,i_ir,i_spn)%s1ai(n_occs+n_holes,n_virs+n_holes), &
               cpks(i_grad,i_ir,i_spn)%h1ai(n_occs+n_holes,n_virs+n_holes), &
               stat=cpksalloc(2))
      ASSERT(cpksalloc(2).eq.0)
      cpksalloc(162)=0 !!! s1ai to be deallocated sepatately
      cpksalloc(167)=0 !h1ai
      MEMLOG(size(cpks(i_grad,i_ir,i_spn)%s1)*2)
      MEMLOG(size(cpks(i_grad,i_ir,i_spn)%s1ai)*2)
      cpks_size=cpks_size+size(cpks(i_grad,i_ir,i_spn)%s1)*2 &
                         +size(cpks(i_grad,i_ir,i_spn)%s1ai)*2

             cpks(i_grad,i_ir,i_spn)%s1=0.0_r8_kind
             cpks(i_grad,i_ir,i_spn)%h1=0.0_r8_kind
             cpks(i_grad,i_ir,i_spn)%s1ai=0.0_r8_kind
             cpks(i_grad,i_ir,i_spn)%h1ai=0.0_r8_kind
   enddo

       DPRINT  'occp occ  and vir',i_ir,i_spn
       do i_occ=1,n_occs+n_holes
        DPRINT i_occ,cpks3c(i_ir,i_spn)%parocc_occ(i_occ)
       enddo
       do i_vir=1,n_holes+n_virs
        DPRINT i_vir,cpks3c(i_ir,i_spn)%parocc_vir_and_holes(i_vir)
       enddo
  enddo spn

     enddo irr_cpks
    endif cpks1

    if (operations_solvation_effect) then                                     !!!!!!!!!!!!!
       allocate(grad_solv_totsym(gradient_data_n_gradients),stat=alloc_stat)  !!!!!!!!!!!!!
       if ( alloc_stat/=0 ) call error_handler( &                             !!!!!!!!!!!!!
            "gradient_data_setup: allocate of grad_solv_totsym failed" )      !!!!!!!!!!!!!
       grad_solv_totsym=0.0_r8_kind                                           !!!!!!!!!!!!!
       if(with_pc .and. .not.fixed_pc) call init_forces_on_pc()               !!!!!!!!!!!!!
       if(integralpar_2dervs) then                                            !!!!!!!!!!!!!
          allocate(grad_solv_totsym_tes(gradient_data_n_gradients,N_points), &
               Q_grad(gradient_data_n_gradients,N_points), stat=alloc_stat)
          if ( alloc_stat/=0 ) call error_handler( &
               "gradient_data_setup: allocate of grad_solv_totsym_tes or Q_grad failed" )
          grad_solv_totsym_tes=0.0_r8_kind
          Q_grad=0.0_r8_kind
       end if
    endif                                                                     !!!!!!!!!!!!!

    if (split_gradients) then
       allocate(gradient_ob_pulay(gradient_data_n_gradients), &
                gradient_ch_pulay(gradient_data_n_gradients), stat=alloc_stat)
       if ( alloc_stat/=0 ) call error_handler( &
            "gradient_data_setup: allocate of gradient_xx_pulay failed" )
       gradient_ob_pulay=0.0_r8_kind
       gradient_ch_pulay=0.0_r8_kind
       if (model_density) then
          allocate(gradient_mda_rho(gradient_data_n_gradients),&
                   gradient_mda_vxc(gradient_data_n_gradients),stat=alloc_stat)
          if ( alloc_stat/=0 ) call error_handler( &
               "gradient_data_setup: allocate of gradient_mda_xxx failed" )
          gradient_mda_rho=0.0_r8_kind
          gradient_mda_vxc=0.0_r8_kind
       endif
    endif

    ! allocate  data structure for Fitfunction-gradients
    allocate(grad_fit_ch_cartesian(3,N_moving_unique_atoms),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    grad_fit_ch_cartesian = 0.0_r8_kind
    if(integralpar_cpksdervs) then
     allocate(cpks_grad_fit_ch_cartesian(n_fit%n_ch,3,N_moving_unique_atoms),&
                                         stat=cpksalloc(26))
     ASSERT(cpksalloc(26).eq.0)
     cpks_grad_fit_ch_cartesian=0.0_r8_kind

    endif
    if (model_density .and. split_gradients) then
       allocate(grad_mda_rhofit_cartesian(3,N_moving_unique_atoms), &
                grad_mda_xcpot_cartesian(3,N_moving_unique_atoms), &
                STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("gradient_data_setup: allocation (3m) failed")
       grad_mda_rhofit_cartesian = 0.0_r8_kind
       grad_mda_xcpot_cartesian = 0.0_r8_kind
    endif
    if(options_relativistic) then
       call gradient_data_rel_setup()
    end if

    if (moving_pc) then
       call init_pointcharges_grads()
#ifdef WITH_EFP
       if(operations_solvation_effect) call init_X_solv_grads()
#endif
    end if
    if (moving_X_centers .or. moving_R_centers) call init_X_centers_grads()
#ifdef WITH_EFP
    if (moving_R_centers) call init_repf_grads()
#endif
    if (moving_Pol_centers) then
       call init_Pol_centers_grads()
    end if

  end subroutine gradient_data_setup

   subroutine cpks_data_shutdown()
    use cpksdervs_matrices
   integer(kind=i4_kind):: i_ir,i_spn,i_grad
     do i_ir=1,size(cpks,2)
      do i_spn=1,size(cpks,3)

     ! see equivalent block in cpks_gvec calculations
!      MEMLOG(-size(cpks3c(i_ir,i_spn)%co_ai))
!      deallocate(cpks3c(i_ir,i_spn)%co_ai,stat=cpksalloc(8))
!      ASSERT(cpksalloc(8).eq.0)

     MEMLOG(-size(cpks3c(i_ir,i_spn)%eigvec_vir_and_holes))
     MEMLOG(-size(cpks3c(i_ir,i_spn)%eigval_vir_and_holes))
     MEMLOG(-size(cpks3c(i_ir,i_spn)%parocc_vir_and_holes-size(cpks3c(i_ir,i_spn)%parocc_occ)))
     deallocate(cpks3c(i_ir,i_spn)%eigvec_vir_and_holes, &
                cpks3c(i_ir,i_spn)%eigval_vir_and_holes, &
                cpks3c(i_ir,i_spn)%parocc_vir_and_holes, &
                cpks3c(i_ir,i_spn)%parocc_occ, &
                stat=cpksalloc(135) )
     ASSERT(cpksalloc(135).eq.0)
     cpksalloc(135)=1

      do i_grad=1,size(cpks,1)

      if(associated(cpks(i_grad,i_ir,i_spn)%s1)) then
       MEMLOG(-size(cpks(i_grad,i_ir,i_spn)%s1)-size(cpks(i_grad,i_ir,i_spn)%h1))
       deallocate(cpks(i_grad,i_ir,i_spn)%s1, &
                  cpks(i_grad,i_ir,i_spn)%h1,stat=cpksalloc(2))
       ASSERT(cpksalloc(2).eq.0)
       cpksalloc(2)=1 ! s1 h1
      endif

       MEMLOG(-size(cpks(i_grad,i_ir,i_spn)%HBH))
       deallocate( cpks(i_grad,i_ir,i_spn)%HBH,   &
                   stat=cpksalloc(166))
      ASSERT(cpksalloc(166).eq.0)
      cpksalloc(166)=1  !HBH

      ! in relativistic case only:
      if( i_spn == 0 )then
        if( associated(cpks(i_grad,i_ir,i_spn)%hr1) )then
          MEMLOG(-size(cpks(i_grad,i_ir,i_spn)%hr1))
          deallocate( cpks(i_grad,i_ir,i_spn)%hr1, stat=cpksalloc(2) )
          ASSERT(cpksalloc(2).eq.0)
        endif
      endif

      enddo
      enddo

     enddo

     MEMLOG(-size(cpks)-size(cpks3c))
     deallocate(cpks,cpks3c,stat=cpksalloc(1))
      ASSERT(cpksalloc(1).eq.0)
             cpksalloc(1)=1
             cpksalloc(7)=1

      call cpksalloc_stat()
   end subroutine cpks_data_shutdown

  subroutine gradient_data_shutdown()
    !  Purpose: deallocating
    ! ---------- modules ---------------------------------------
    use pointcharge_module, only: moving_pc,pc_grads_shutdown
!   use calc3c_switches
    use integralpar_module, only: integralpar_2dervs,integralpar_cpksdervs
    use cpksdervs_matrices
#ifdef WITH_EFP
    use efp_solv_grad_module, only: X_solv_grads_shutdown
#endif
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: alloc_stat,i_unique!,i_ir,i_grad
!   integer(kind=i4_kind) :: i_spn,i
    integer(kind=i4_kind) :: i_unique2
    logical               :: split_gradients, model_density
    !------------ Executable code -------------------------------
    split_gradients = options_split_gradients()
    model_density = options_xcmode() == xcmode_model_density

    do i_unique=1,N_moving_unique_atoms
       deallocate(gradient_cartesian(i_unique)%m,stat=alloc_stat)
       ASSERT(alloc_stat.eq.0)

     if(integralpar_2dervs) then
      do i_unique2=1,N_moving_unique_atoms
       deallocate(dervs_cartesian(i_unique,i_unique2)%m,stat=cpksalloc(91))
       ASSERT(cpksalloc(91).eq.0)
       cpksalloc(91)=1
      enddo
     endif

       if(integralpar_cpksdervs) then
       deallocate(cpks_gradient_cartesian(i_unique)%m,stat=cpksalloc(29))
       ASSERT(cpksalloc(29).eq.0)
       cpksalloc(29)=1
       endif

       if (operations_solvation_effect) then                                       !!!!!!!!!!!!!
          deallocate(grad_solv_cart(i_unique)%m,stat=alloc_stat)
          if ( alloc_stat/=0 ) call error_handler( &
               "gradient_data_shutdown: deallocate of grad_solv_cart%m failed" )
       endif
       if (split_gradients) then
          deallocate(partial_cartesian(i_unique)%m,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler("gradient_data_shutdown: &
               &deallocate of partial_cartesian%m failed" )
       endif
    enddo
    if (split_gradients) then
       deallocate(gradient_ob_pulay,gradient_ch_pulay,stat=alloc_stat)
       if ( alloc_stat/=0 ) call error_handler &
            ("gradient_data_shutdown: deallocate of gradient_xx_pulay failed")
       if (model_density) then
          deallocate(gradient_mda_rho,gradient_mda_vxc,stat=alloc_stat)
          if ( alloc_stat/=0 ) call error_handler &
               ("gradient_data_shutdown: deallocate of gradient_mda_xxx failed")
       endif
       deallocate(partial_cartesian,stat=alloc_stat)
       if ( alloc_stat/=0 ) call error_handler &
            ("gradient_data_shutdown: deallocate of partial_cartesian failed")
    endif

#ifdef WITH_EPE
    if(ewpc_n.ne.0) then
        deallocate(cluster_epe_energy)
    endif
#endif

    deallocate(gradient_totalsym,gradient_index,gradient_cartesian,&
               stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

    if(integralpar_2dervs) then
    MEMLOG(-size(dervs_totalsym)-size(dervs_cartesian))
    deallocate(dervs_totalsym, dervs_cartesian,&
               stat=cpksalloc(39))
             ASSERT(cpksalloc(39).eq.0)
             cpksalloc(39)=1
             cpksalloc(90)=1
    endif

     if(integralpar_cpksdervs) then
      deallocate(cpks_gradient_cartesian,stat=cpksalloc(28))
      ASSERT(cpksalloc(28).eq.0)
      cpksalloc(28)=1
     endif


    if (operations_solvation_effect) then                                       !!!!!!!!!!!!!
       deallocate(grad_solv_totsym,grad_solv_cart,stat=alloc_stat)              !!!!!!!!!!!!!
       if ( alloc_stat/=0 ) call error_handler( &                               !!!!!!!!!!!!!
            "gradient_data_shutdown: deallocate of grad_solv_totsym failed" )   !!!!!!!!!!!!!
       if(integralpar_2dervs) then                                              !!!!!!!!!!!!!
          deallocate(grad_solv_totsym_tes,Q_grad,stat=alloc_stat)
          if ( alloc_stat/=0 ) call error_handler( &
               "gradient_data_shutdown: deallocate of grad_solv_totsym_tes or Q_grad failed" )
       end if
       call dealloc_geom_grad()                                                 !!!!!!!!!!!!!
    endif                                                                       !!!!!!!!!!!!!
    deallocate(grad_fit_ch_cartesian,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    if(integralpar_cpksdervs) then
    MEMLOG(-size(cpks_gradient_totalsym)-size(cpks_grad_fit_ch_cartesian))
    MEMLOG(-size(cpks_grad_totalsym))
    deallocate(cpks_grad_fit_ch_cartesian,cpks_gradient_totalsym, &
               cpks_grad_totalsym, &
                                                STAT=cpksalloc(26))
    ASSERT(cpksalloc(26).eq.0)
           cpksalloc(26)=1
           cpksalloc(27)=1 ! cpks_gradient_totalsym
    MEMLOG(-size(cpks_grad_fit_totasym))
    deallocate(cpks_grad_fit_totasym, &
               STAT=cpksalloc(164))
    ASSERT(cpksalloc(164).eq.0)
    cpksalloc(164)=1 ! cpks_grad_fit_totasym
    endif

    if (model_density .and. split_gradients) then
       deallocate(grad_mda_rhofit_cartesian, &
                  grad_mda_xcpot_cartesian, STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    endif
    if(options_relativistic) then
       call gradient_data_rel_shutdown()
    end if

    if(allocated(cpks).and.integralpar_cpksdervs)  call cpks_data_shutdown()

    if(moving_pc) then
       call pc_grads_shutdown()
#ifdef WITH_EFP
       if(operations_solvation_effect) call X_solv_grads_shutdown()
#endif
    end if
    if(moving_X_centers .or. moving_R_centers) call X_centers_grads_shutdown()
#ifdef WITH_EFP
    if(moving_R_centers) call repf_grads_shutdown()
#endif
    if(moving_Pol_centers) then
       call Pol_centers_grads_shutdown()
    end if

  end subroutine gradient_data_shutdown

  subroutine gradient_data_rel_setup()
    !  Purpose: allocation of necesarry arrays
    !** End of interface ****************************************
!   use calc3c_switches, only: integralpar_2dervs
    use integralpar_module, only: integralpar_2dervs
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: n_irrep,n_host,i_host,i_irrep,alloc_stat
    integer(kind=i4_kind) :: x,y
    !------------ Executable code -------------------------------
    n_host=comm_get_n_processors()
    n_irrep=symmetry_data_n_irreps()
    allocate( gradient_rel_index(gradient_data_n_gradients,n_irrep),&
          stat=alloc_stat)
    ASSERT(alloc_stat==0)

    ! assign work by round-robin over procs:
    i_host = 0
    do i_irrep=1,n_irrep
      do x=1,gradient_data_n_gradients
        gradient_rel_index(x,i_irrep) = 1 + mod(i_host,n_host)
        i_host = i_host + 1
#ifdef FPP_DEBUG
        if( comm_i_am_master() )then
          WRITE(*,'(X,"gradient_data_rel_setup: grad(",I6,") of irrep",I3," on proc",I3)') &
               x, i_irrep, gradient_rel_index(x,i_irrep)
        endif
#endif
      enddo
    enddo

    if(.not.integralpar_2dervs) RETURN

    ! continue with assignment of second derivatives
    allocate(gradient_der_index( gradient_data_n_gradients &
                               , gradient_data_n_gradients &
                               , n_irrep                   &
                               )                           &
            ,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    ! assign work:
    i_host = 0
    do i_irrep=1,n_irrep
      do y=1,gradient_data_n_gradients
        do x=1,y
          gradient_der_index(x,y,i_irrep) = 1 + mod(i_host,n_host)
          gradient_der_index(y,x,i_irrep) = gradient_der_index(x,y,i_irrep)
          i_host = i_host + 1
!         if( comm_i_am_master() )then
!           WRITE(*,'(X,"gradient_data_rel_setup: derv(",2I3,") of irrep",I3," on proc",I3)') &
!                x, y, i_irrep, gradient_der_index(x,y,i_irrep)
!         endif
        enddo
      enddo
    enddo
  end subroutine gradient_data_rel_setup
  !**************************************************************

  subroutine gradient_data_rel_shutdown()
    !  Purpose: deallocating
    ! ---------- modules ---------------------------------------
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: alloc_stat !,i_irrep,n_irrep
    !------------ Executable code -------------------------------

    deallocate( gradient_rel_index,stat=alloc_stat)
    if(alloc_stat/=0) call error_handler(&
         'gradient_data_rel_shutdown: deallocating failed')

    if(allocated(gradient_der_index))then
      ! only with second derivatives:
      deallocate(gradient_der_index,stat=alloc_stat)
      ASSERT(alloc_stat==0)
    endif
  end subroutine gradient_data_rel_shutdown

  subroutine gradient_data_write_cartesians(header,grad_final)
    use pointcharge_module
    !  Purpose: writing the cartesian gradients
    implicit none
    character(len=*) :: header
    type(arrmat2)    :: grad_final(:)
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: i_unique,i_center,i_equal
    !------------ Executable code -------------------------------
    write(output_unit,'(/A)') header
    do i_unique=1,N_moving_unique_atoms
       i_center = moving_unique_atom_index(i_unique)
       write(output_unit,*) 'Unique Center:',i_center
       do i_equal=1,unique_atoms(i_center)%n_equal_atoms
          write(output_unit,'(A20,3F20.12)') 'Equal Center:',&
               grad_final(i_unique)%m(:,i_equal)
       enddo
    end do
    if(N_moving_unique_timps.ne.0) then
       do i_unique=1,N_moving_unique_timps
          i_center = moving_unique_timp_index(i_unique)
          write(output_unit,*) 'Unique Center:',i_center
          do i_equal=1,unique_timps(i_center)%n_equal_atoms
             write(output_unit,'(A20,3F20.12)') 'Equal Center:',&
                  grad_final(i_unique+N_moving_unique_atoms)%m(:,i_equal)
          enddo
       end do
    end if

  end subroutine gradient_data_write_cartesians

#ifdef WITH_SECDER
  subroutine gradient_data_write_cart_hess(iou,h)
    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    integer(i4_kind)  , intent(in) :: iou
    type(arrmat4), intent(in) :: h(:,:) ! (NUA,NUA)
    ! *** end of interface ***

    integer(i4_kind)           :: na
    real(r8_kind), allocatable :: k(:,:), m(:)
    real(r8_kind), allocatable :: r(:)
    real(r8_kind), allocatable :: w(:), mef(:), ang(:,:), rms(:,:)
    real(r8_kind), allocatable :: u(:,:), v(:,:), intensity(:)

    na = sum( uas(:)%n_equal_atoms )

    allocate( k(3*na,3*na), m(3*na), r(3*na) )
    allocate( w(3*na), mef(3*na), ang(3,3*na), rms(3,3*na) )
    allocate( u(3*na,3*na), v(3*na,3*na), intensity(3*na) )

    ! ``unwind'' legacy hessian into plain matrix:
    call cart_hess_mat(iou,h,k)

    ! build (diagonal) mass matrix:
    call cart_mass_mat(  0,uas,m) ! dont print masses

    ! build cartesian coordinate vector (optional):
    call cart_coor_vec (uas, r) ! does not print

    ! compute the frequencies:
    call freq(k,m,r,w,mef,ang,rms,u,v,intensity)
    call freq_print(iou,m,r,w,mef,ang,rms,u,v,intensity)

    deallocate( k, m, r )
    deallocate( w, mef, ang, rms, u, v, intensity )
  end subroutine gradient_data_write_cart_hess
#endif

  subroutine optimizer_write_cart_hess(iou,k,charge,xyz,na)
!    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    integer(i4_kind)  , intent(in) :: iou,na
!    type(arrmat4), intent(in) :: h(:,:) ! (NUA,NUA)
    real(r8_kind),intent(in):: k(:,:),charge(:),xyz(:,:)
    ! *** end of interface ***

    real(r8_kind), allocatable :: m(:)
    real(r8_kind), allocatable :: r(:)
    real(r8_kind), allocatable :: w(:), mef(:), ang(:,:), rms(:,:)
    real(r8_kind), allocatable :: u(:,:), v(:,:), intensity(:)

    integer(i4_kind)::i

    allocate(  m(3*na), r(3*na) )
    allocate( w(3*na), mef(3*na), ang(3,3*na), rms(3,3*na) )
    allocate( u(3*na,3*na), v(3*na,3*na), intensity(3*na) )

!    ! ``unwind'' legacy hessian into plain matrix:
!    call cart_hess_mat(iou,h,k)

    ! build (diagonal) mass matrix:
    call optimizer_cart_mass_mat(  0,charge,na,m) ! dont print masses

    ! build cartesian coordinate vector (optional):
     do i = 1, na
        r(3*i-2)=xyz(1,i)
        r(3*i-1)=xyz(2,i)
        r(3*i  )=xyz(3,i)
     enddo


    ! compute the frequencies:
    call freq(k,m,r,w,mef,ang,rms,u,v,intensity)
    call freq_print(iou,m,r,w,mef,ang,rms,u,v,intensity)

    deallocate( m, r )
    deallocate( w, mef, ang, rms, u, v, intensity )
  end subroutine optimizer_write_cart_hess
  !**************************************************************

  !**************************************************************
  subroutine gradient_data_write_gxfile(loop, coordinates_only)
    !  Purpose: writing the cartesian gradients to gxfile
    !** End of interface ****************************************
    use iounitadmin_module
    use filename_module, strlen=>filename_namelengthmax
    use energy_calc_module
    use fermi_module
    use operations_module, only: OPERATIONS_gx_epeformat
    implicit none
    integer(kind=i4_kind), intent(in) :: loop
    logical, optional, intent(in) :: coordinates_only ! default = .false.
    ! *** end of interface ***

    logical :: write_gradients

    DPRINT 'gradient_data_write_gxfile: entered'
    if (present(coordinates_only)) then
       write_gradients = .not. coordinates_only
    else
       write_gradients = .true.
    endif

#ifdef WITH_SIMOL /* never, actually ...*/
    DPRINT 'gradient_data_write_gxfile: call simol_write()'
    call simol_write()
#else
    DPRINT 'gradient_data_write_gxfile: call optimizer_write()'
    call optimizer_write()
#endif

    DPRINT 'done gradient_data_write_gxfile'

  contains

#ifdef WITH_SIMOL /* never, actually ...*/
    subroutine simol_write()
        use pointcharge_module
#ifdef WITH_EPE
        use epecom_module
        use ewaldpc_module, only: cluster_nuc_epe_en,epe_relaxation
#endif
      integer(kind=i4_kind) :: i_unique,i_moving,i_equal,io_u,counter
      real(kind=r8_kind)    :: energy,entropy,gradient(3)
      ! the following variables are used to read  in the whole gxfile
      integer(kind=i4_kind),parameter :: max_atoms=1000
      integer(kind=i4_kind)    :: zmat(3,max_atoms),numx(3,max_atoms) &
                                  ,iepe(max_atoms) !EPE index
      integer(kind=i4_kind)    :: index_unique(max_atoms),index_eq(max_atoms)
      real(kind=r8_kind)       :: x(max_atoms),y(max_atoms),z(max_atoms),&
           charge(max_atoms)
      integer(kind=i4_kind)    :: n_atoms,i,j,alloc_stat
      ! help variables in case that the stepsize is specified for
      ! a simol(frequency) calculation
      integer(kind=i4_kind),allocatable :: sym_equiv(:)
      logical,allocatable               :: change(:)
      real(kind=r8_kind),allocatable    :: step_size(:)
      integer(kind=i4_kind)             :: num_change,min_num,max_num
      ! 'sym_equiv'-array implemented for later use.

      real(kind=r8_kind):: energy2
#ifdef WITH_EPE
      real(kind=r8_kind):: epe_latt_energy,cluster_regI,eshort
#endif

      ! initialize help variables
      x=0.0_r8_kind
      y=0.0_r8_kind
      z=0.0_r8_kind
      charge=0.0_r8_kind
      zmat=0_i4_kind
      numx=0_i4_kind
      iepe=0_i4_kind
      index_unique=0_i4_kind
      index_eq=0_i4_kind
      num_change=0_i4_kind

      io_u=openget_iounit(status='old',form='formatted',file=&
           trim(inpfile('gxfile')))

      ! now read the gxfile
      n_atoms=0
        call write_to_trace_unit('simol_write gxfile recreation')
        if(OPERATIONS_gx_epeformat) &
         call write_to_trace_unit('free epe format  read')
      input_loop: do i=1,max_atoms
         if (operations_gx_highprec) then
            if(OPERATIONS_gx_epeformat) then
              read(io_u,*) charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)

               else
            read(io_u,1020)charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i)
            endif
         else
          if(OPERATIONS_gx_epeformat) then
           read(io_u,*) charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)
          else
            read(io_u,1000)charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)
          endif
         endif
         if (charge(i)<=0.5_r8_kind) then
            exit input_loop
         endif
         n_atoms=n_atoms+1
      enddo input_loop
         call write_to_trace_unit('done')
      backspace io_u

      max_num=maxval(numx)
      min_num=minval(numx)
      ! find out the highest absolute value of numx
      if (min_num<0) then
         min_num=-min_num
         if (min_num>max_num) then
            max_num=min_num
         endif
         !only if there is a negative number in the numx-array we
         ! have to handle step-sizes:
         allocate(sym_equiv(max_num),change(max_num),stat=alloc_stat)
         if (alloc_stat/=0) call error_handler &
              ("gradient_data_writegxfile: allocation (1) failed")
         sym_equiv=0_i4_kind
         change=.false.
         do j=1,3
            do i=1,n_atoms
               sym_equiv(abs(numx(j,i))) = sym_equiv(abs(numx(j,i)))+1
               if (numx(j,i)<0) change(abs(numx(j,i)))=.true.
            enddo
         enddo
         num_change=count(change)
         allocate(step_size(num_change),stat=alloc_stat)
         if (alloc_stat/=0) call error_handler &
              ("gradient_data_writegxfile: allocation (2) failed")
         do i=1,num_change
            read(io_u,*)step_size(i)
         enddo

      endif

      close(io_u,status='delete')

      ! now re-create the gxfile
      open(io_u,status='new',form='formatted',file=&
           trim(inpfile('gxfile')))
      do i=1,n_atoms
         if (operations_gx_highprec) then
            if(operations_gx_epeformat) then
               print*,'simol_write '
              write(io_u,1021 )charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)
               else
            write(io_u,1020)charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i)
            endif
         else
            if(operations_gx_epeformat) then
            write(io_u,1001)charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)
            else
            write(io_u,1000)charge(i),x(i),y(i),z(i),&
                 index_unique(i),index_eq(i),&
                 zmat(1,i),zmat(2,i),zmat(3,i),&
                 numx(1,i),numx(2,i),numx(3,i),iepe(i)
             endif ! operations_ewpc
         endif
      enddo

      if (operations_gx_highprec) then
         if(operations_gx_epeformat) then
           write(io_u,1031) charge(n_atoms+1),x(n_atoms+1), &
                y(n_atoms+1),z(n_atoms+1),&
                index_unique(n_atoms+1),index_eq(n_atoms+1),&
                zmat(:,n_atoms+1),numx(:,n_atoms+1),iepe(n_atoms+1)
           print *,' gradient_data'
           print 1031,charge(n_atoms+1),x(n_atoms+1), &
                y(n_atoms+1),z(n_atoms+1),&
                index_unique(n_atoms+1),index_eq(n_atoms+1),&
                zmat(:,n_atoms+1),numx(:,n_atoms+1),iepe(n_atoms+1)

            else
         write(io_u,1030)charge(n_atoms+1),x(n_atoms+1),y(n_atoms+1),z(n_atoms+1),&
                 index_unique(n_atoms+1),index_eq(n_atoms+1),&
                 zmat(1,n_atoms+1)
         endif
      else
         write(io_u,1010)charge(n_atoms+1),x(n_atoms+1),y(n_atoms+1),z(n_atoms+1),&
                 index_unique(n_atoms+1),index_eq(n_atoms+1),&
                 zmat(:,n_atoms+1),numx(:,n_atoms+1),iepe(n_atoms+1)
      endif
         call write_to_trace_unit('gx-recreated')

      ! if changed step-sizes are present write them now
      if (num_change/=0) then
         do i=1,num_change
            write(io_u,*)step_size(i)
         enddo
         deallocate(step_size,change,sym_equiv,stat=alloc_stat)
         if (alloc_stat/=0) call error_handler &
              ("gradient_data_writegxfile: deallocation (1) failed")
      endif

      if (write_gradients) then
      counter=1
      call get_energy(tot=energy)
      energy2=energy

#ifdef WITH_EPE
      if(epe_relaxation) then
         call write_to_trace_unit('get_epe_energies')
         call get_epe_energies( &
              lattice_energy=epe_latt_energy,epg_cluster_reg_I=cluster_regI,eshort_coupling_au=eshort)
         energy=energy+epe_latt_energy
                      !in epe_latt_energy the term  cluster_regI is substrected from etot_epe
         energy2=energy-cluster_regI+(cluster_epe_energy(1)-cluster_nuc_epe_en)
                       !coulomb orbital interaction
                                                           !coulomb nuclear interact
                                    !coulomb fit and nuclear
         print*,'energy,energy2,eshort',energy,energy2,eshort
         print*,'epe_side_optimized_energy', energy+eshort
         epe_side_optimized_energy=energy+eshort
         print*,'SCF cluster_regI contribution', cluster_regI
         print*,'epe_latt_energy',epe_latt_energy
         print*,'cluster_epe_energy-cluster_nuc_epe_en', &
              (cluster_epe_energy(1)-cluster_nuc_epe_en)
      endif
#endif

      entropy=fermi_get_entropy()
      if (.not.operations_gx_highprec) then
         write(io_u,'(2F20.8,2x,3I5)')energy-entropy,&
              energy-entropy,counter,counter,counter
      else
         write(io_u,'(2F24.12,2x,3I5)') energy-entropy,&
              energy-entropy,counter,counter,counter
      endif
      do i_unique=1,n_unique_atoms
         i_moving = unique_atoms(i_unique)%moving_atom
         do i_equal=1,unique_atoms(i_unique)%n_equal_atoms
            if (i_moving == 0) then
               gradient = 0.0_r8_kind
            else
               gradient = gradient_cartesian(i_moving)%m(:,i_equal)
            endif
            if (.not.operations_gx_highprec) then
               write(io_u,'(I5,5X,3F15.6)') counter,gradient
            else
               write(io_u,'(I5,5x,3F17.12)') counter,gradient
            endif
            counter=counter+1
         enddo
      end do

      if(n_moving_unique_timps.ne.0) then
         do i_unique=1,n_timps
            i_moving = unique_timps(i_unique)%moving_atom
            do i_equal=1,unique_timps(i_unique)%n_equal_atoms
               if (i_moving == 0) then
                  gradient = 0.0_r8_kind
               else
                  gradient = gradient_cartesian( &
                       n_unique_atoms+i_moving)%m(:,i_equal)
               endif
               if (.not.operations_gx_highprec) then
                  write(io_u,'(I5,5X,3F15.6)') counter,gradient
               else
                  write(io_u,'(I5,5x,3F17.12)') counter,gradient
               endif
               counter=counter+1
            enddo
         end do
      end if

      endif

      close(io_u)
      call return_iounit(io_u)

1000  format((f5.2,3(2X,f13.7),2i3,2x,3I3,2X,3I3,i5)) ! simol format
1001  format((f5.2,3(f15.7),2i4,i5,2I4,i5,2I4,i5)) ! gx_epe_f extended
1020  format((f5.2,3(2X,f21.12),2i3,2x,3I3,2X,3I3)) ! simol format
1021  format((f5.2,3(2X,f21.12),2i4,2x,3I4,2X,3I4,i5)) ! gx_epe_f extended
1010  format(f5.1,3(2X,f13.7),2i3,2x,3I3,2X,3I3,i5)
1030  format(f5.1,3(2X,f21.12),2i3,2x,I3)
1031  format(f5.1,3(2X,f21.12),2i4,2x,3I4,2X,3I4,i5) ! gx_epe_format
    end subroutine simol_write

#else /* i.e. WITH_SIMOL is not def, always actually ... */

    subroutine optimizer_write()
#ifdef WITH_EPE
      use epecom_module
      use ewaldpc_module, only: cluster_nuc_epe_en,epe_relaxation
#endif
      use gxfile, only: gxfile_read
      use filename_module, only: max_path=>filename_namelengthmax
      use unique_atom_methods, only: unique_atom_make_gx
      implicit none
      ! *** end of interface ***

      integer (i4_kind) :: i, i_unique, i_moving, counter, i_equal, n_atoms, &
           counter_equal
      real(kind=r8_kind)      :: energy,entropy,gradient(3)
      real(kind=r8_kind),parameter   :: dummy_charge=99.0_r8_kind
      ! the following variables are used to read  in the whole gxfile
      integer(kind=i4_kind),parameter :: max_atoms=1000
      integer(kind=i4_kind)    :: zmat(3,max_atoms),numx(3,max_atoms)
      integer(kind=i4_kind)    :: index_unique(max_atoms),index_eq(max_atoms)
      integer(kind=i4_kind)    :: iepe(max_atoms)
      real(kind=r8_kind)       :: xyz(3,max_atoms),charge(max_atoms)
      ! help variables in case that the stepsize is specified for
      ! a simol(frequency) calculation
      integer(kind=i4_kind)             :: num_change

      real(kind=r8_kind):: energy2
#ifdef WITH_EPE
      real(kind=r8_kind):: epe_latt_energy,cluster_regI,eshort
#else
      !FIXME: unexport it from epecom_module, and use constants:
      real(r8_kind), parameter :: zero=0.0_r8_kind
#endif

      integer(kind=i4_kind) :: n_centers
      real(kind=r8_kind)    :: geo_loop
      logical               :: ex_gx
      character(len=1) :: cha_loop

      ! Path and the descriptor for gx file:
      character (len=max_path) :: file
      integer (i4_kind) :: io_u

      DPRINT 'optimizer_write: entered, epefmt=',operations_gx_epeformat
      ! initialize help variables
      charge=0.0_r8_kind
      zmat=0_i4_kind
      numx=0_i4_kind
      index_unique=0_i4_kind
      index_eq=0_i4_kind
      num_change=0_i4_kind

      if (operations_transit) then
         write (cha_loop, '(i1)') loop
         file = trim (inpfile ('gx.' // cha_loop))
      else
         file = trim (inpfile ('gxfile'))
      endif

      ! See if the file exists in the input directory.  If it did not,
      ! create one. Older  version of this sub created  gx file in the
      ! output directory instead.
      inquire (file= file, exist= ex_gx)
      if (.not. ex_gx) then
         call unique_atom_make_gx (iloop= 1)

         ! FIXME: I set ex_gx to true  now (even though it is not used
         ! anywhere  below).  The  file  does exists  now.  The  older
         ! logic only executed the bulk  of the code below only if the
         ! file existed before, otherwise  it just created a fresh one
         ! and did  very little (apart from trying  to write gradients
         ! to a non-initialized unit).
         ex_gx = .true.         ! not used below
      endif

      ! Make sure it is initalized before use:
      io_u = openget_iounit (file, status= 'old', form= 'formatted')

      ! now read the gxfile
      if(operations_gx_epeformat)then
         call gxfile_read( io_u, n_centers, geo_loop, charge, xyz, &
              & index_unique, index_eq, zmat, numx, &
              & iepe )
      else
         call gxfile_read( io_u, n_centers, geo_loop, charge, xyz, &
              & index_unique, index_eq, zmat, numx)
      endif
      n_atoms = n_centers ! counting dummies?

      ! Now re-adjust the positions of the equal atoms in case 'optimizer' moved
      ! them out of their symmetry-equivalent positions due to numerical
      ! precision limits, the positions of dummy atoms remain fixed
      counter=0
      do i=1,N_unique_atoms
         counter_equal=0
         equals: do
            counter=counter+1
            if (charge(counter)==dummy_charge.or.index_unique(counter).eq.0) cycle equals
            counter_equal=counter_equal+1
            if (counter_equal==1) then  !this is the 'unique'
               xyz(:,counter) = unique_atoms(i)%position_first_ea(:)
            else
               xyz(:,counter) = unique_atoms(i)%position(:,counter_equal)
            endif
            if (counter_equal==unique_atoms(i)%N_equal_atoms) exit equals
         enddo equals
      enddo

      if (operations_transit) then
         call returnclose_iounit (io_u, status= 'keep')
      else
         call returnclose_iounit (io_u, status= 'delete')
      endif
      io_u = -1

      DPRINT 're-write the gxfile with its re-adjusted atomic positions'
      ! Write  unconditionally  the staff  into  gxfile  in the  input
      ! directory. FIXME: should it be "file" again?
      io_u = openget_iounit (file= trim (inpfile('gxfile')), &
           status= 'unknown', form= 'formatted')

      DPRINT 'operations_gx_epeformat=',operations_gx_epeformat
      do i=1,n_atoms
         if(operations_gx_epeformat) then
            write(io_u,1021)charge(i),&
                 & xyz(1:3,i), &
              index_unique(i),index_eq(i),&
              zmat(1,i),zmat(2,i),zmat(3,i),&
              numx(1,i),numx(2,i),numx(3,i), iepe(i)
            else
         write(io_u,1020)charge(i),&
              & xyz(1:3,i), &
              index_unique(i),index_eq(i),&
              zmat(1,i),zmat(2,i),zmat(3,i),&
              numx(1,i),numx(2,i),numx(3,i)
         endif
      enddo

      if(operations_gx_epeformat) then
        print*,'optimizer_write: instruction line', charge(n_atoms+1)

         write(io_u,1031)charge(n_atoms+1), &
              zero,zero,zero,&
              index_unique(n_atoms+1),index_eq(n_atoms+1),&
              0,0,0, 0,0,0, 0

         else
      write(io_u,1030)charge(n_atoms+1), &
            zero,zero,zero,&
              index_unique(n_atoms+1),index_eq(n_atoms+1),&
              0,0,0, 0,0,0
      end if

      if (.not.write_gradients) goto 900 ! skip writing gradients ...

      counter=1
      call get_energy(tot=energy)
      energy2=energy

#ifdef WITH_EPE
      if(epe_relaxation) then
         call write_to_trace_unit('optimizer_write get_epe_energies')
         call get_epe_energies( &
              lattice_energy=epe_latt_energy,epg_cluster_reg_I=cluster_regI, &
              eshort_coupling_au=eshort)
         energy=energy+epe_latt_energy
                      !in epe_latt_energy the term  cluster_regI is substrected from etot_epe
         print*,epe_latt_energy, 'epe_latt_energy to be summed with SCF energy'
         print*,'cluster_regI', cluster_regI
         print*,'cluster_nuc_epe_en',cluster_nuc_epe_en

        if(calc_cluster_epe_energy) then
         print*,'cluster_epe_energy(1)',cluster_epe_energy(1)
         energy2=energy-cluster_regI+(cluster_epe_energy(1)-cluster_nuc_epe_en)
         print*,'energy,energy2,eshort',energy,energy2,eshort
         print*,'cluster_epe_energy-cluster_nuc_epe_en', &
              cluster_epe_energy(1) - cluster_nuc_epe_en
       endif

         print*,'epe_side_optimized_energy', energy+eshort
         epe_side_optimized_energy=energy+eshort
      endif
#endif

      entropy=fermi_get_entropy()

      write(io_u,'(2F24.12,2x,3I5)')energy-entropy,&
           energy-entropy,counter,counter,counter
      do i_unique=1,n_unique_atoms
         i_moving = unique_atoms(i_unique)%moving_atom
         do i_equal=1,unique_atoms(i_unique)%n_equal_atoms
            if (i_moving == 0) then
               gradient = 0.0_r8_kind
            else
               gradient = gradient_cartesian(i_moving)%m(:,i_equal)
            endif
            write(io_u,'(I5,5x,3F17.12)') counter,gradient
            counter=counter+1
         enddo
      end do

900   CONTINUE ! by writing dipole moments ...
      if ( output_dipole_optimizer ) then
         ! Append dipole moments to gxfile. They are used in
         ! numerical frequency calculations to estimate the intensities:
         call dipole_write_optimizer(io_u)
      endif

      ! close gxfile and exit
      call returnclose_iounit (io_u)
      io_u = -1
      RETURN

1020  format((f7.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4)) ! simol format
1021  format((f7.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)) ! simol format
1030  format(f7.1,3(2X,f21.12),2i4,2x,3I4,2X,3I4)
1031  format(f7.1,3(2X,f21.12),2i4,2x,3I4,2X,3I4,i5)
    end subroutine optimizer_write

#endif /* ifdef WITH_SIMOL */

  end subroutine gradient_data_write_gxfile
  !**************************************************************

  !**************************************************************
#ifdef WITH_MOLMECH
  subroutine qm_grads_to_qmmm()

    use energy_calc_module
    use qmmm_interface_module

    integer(i4_kind) :: i,i_unique,i_equal,i_moving
    real(r8_kind) :: energy

    call get_energy(tot=energy)
    energy_qm=energy

    i=0
    do i_unique=1,n_unique_atoms
       i_moving = unique_atoms(i_unique)%moving_atom
       do i_equal=1,unique_atoms(i_unique)%n_equal_atoms
          i=i+1
          if (i_moving == 0) then
             grad_qm(i)%x = 0.0_r8_kind
             grad_qm(i)%y = 0.0_r8_kind
             grad_qm(i)%z = 0.0_r8_kind
          else
             grad_qm(i)%x = gradient_cartesian(i_moving)%m(1,i_equal)
             grad_qm(i)%y = gradient_cartesian(i_moving)%m(2,i_equal)
             grad_qm(i)%z = gradient_cartesian(i_moving)%m(3,i_equal)
          endif
       enddo
    end do

  end subroutine qm_grads_to_qmmm
#endif

  subroutine gradient_sndrcv_3c(what_to_do)
    ! to be called from parallel context
    use comm_module
    use msgtag_module, only: msgtag_grad_3c
    implicit none
    ! *** end of interface ***
    integer :: what_to_do

    if( comm_i_am_master() )then
       call gradient_receive_3c(what_to_do)
    else
       call gradient_send_3c(what_to_do)
    endif

  end subroutine gradient_sndrcv_3c

  subroutine gradient_receive_3c(what_to_do)
    !  Purpose: receive 3 center contribution to the gradient from the
    !           slave
    use cpksdervs_matrices,only:cpks,cpksalloc
    use integralpar_module, only: integralpar_2dervs
    use pointcharge_module, only: moving_pc,totsym_PC_grad_unpack
#ifdef WITH_EFP
    use efp_solv_grad_module, only: totsym_X_solv_grad_unpack
#endif
    implicit none

    integer :: what_to_do
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: n_procs,i,info
    integer(kind=i4_kind) :: i_spn,i1
    integer(kind=i4_kind) :: n_occ,i_g,i_grad,n_vac
    real(kind=r8_kind) :: help_arr(gradient_data_n_gradients)
    real(r8_kind) :: help_arr1(totsym_field_length)
#ifdef WITH_EPE
    real(r8_kind) :: buffer
#endif
    real(kind=r8_kind),allocatable :: sh1_help_arr(:,:),ai_help_arr(:,:)
    real(kind=r8_kind),allocatable:: help_cpks_grad_totalsym(:,:)
    integer, parameter :: grads=1, deriv=2
    !------------ Executable code -------------------------------

    n_procs=comm_get_n_processors()
    if(what_to_do==deriv .and. allocated(cpks)) then
        allocate(help_cpks_grad_totalsym(size(cpks_gradient_totalsym,1),size(cpks,1)), &
                 stat=cpksalloc(123))
        ASSERT(cpksalloc(123).eq.0)
    endif

    do i=2,n_procs
       call comm_save_recv(i,msgtag_grad_3c)

     if(what_to_do==grads) then
       if(comm_msgtag()/=msgtag_grad_3c) &
          call error_handler('Wrong msgtag in subroutine gradient_receive_3c')
       call communpack(help_arr,gradient_data_n_gradients,1,info)
       ASSERT(info.eq.0)
       gradient_totalsym=gradient_totalsym+help_arr
     end if

     if(what_to_do==deriv .and. allocated(cpks)) then

       call communpack(help_cpks_grad_totalsym(1,1),size(cpks_gradient_totalsym),1,info)
       ASSERT(info.eq.0)
       cpks_gradient_totalsym=cpks_gradient_totalsym+help_cpks_grad_totalsym

     do i_g=1,symmetry_data_n_irreps()
        do i_spn=1,size(cpks,3)
        n_occ=size(cpks(1,i_g,i_spn)%s1,1)
        n_vac=size(cpks(1,i_g,i_spn)%s1ai,2)

       if(n_occ.gt.0) then
        allocate(sh1_help_arr(n_occ,n_occ), &
                 ai_help_arr(n_occ,n_vac),stat=cpksalloc(122))
         ASSERT(cpksalloc(122).eq.0)

        gr: do i_grad=1,size(cpks,1)

         call communpack(sh1_help_arr(1,1),n_occ**2,1,info)
         ASSERT(info.eq.0)
         cpks(i_grad,i_g,i_spn)%s1=cpks(i_grad,i_g,i_spn)%s1+sh1_help_arr

         call communpack(ai_help_arr(1,1),n_occ*n_vac,1,info)
         ASSERT(info.eq.0)
         cpks(i_grad,i_g,i_spn)%s1ai=cpks(i_grad,i_g,i_spn)%s1ai+ai_help_arr

         call communpack(ai_help_arr(1,1),n_occ*n_vac,1,info)
         ASSERT(info.eq.0)
         cpks(i_grad,i_g,i_spn)%h1ai=cpks(i_grad,i_g,i_spn)%h1ai+ai_help_arr

        enddo gr

        deallocate(sh1_help_arr,ai_help_arr,stat=cpksalloc(122))
         ASSERT(cpksalloc(122).eq.0)
         cpksalloc(122)=1
       endif
       enddo
      enddo
     endif

     if(what_to_do==grads) then
       if (operations_solvation_effect) then                                !!!!!!!!!!!!
          call communpack(help_arr,gradient_data_n_gradients,1,info)        !!!!!!!!!!!!
          if(info/=0) call error_handler&                                   !!!!!!!!!!!!
               ('gradient_receive_3c: unpacking grad_solv_totsym failed')   !!!!!!!!!!!!
          grad_solv_totsym=grad_solv_totsym+help_arr                        !!!!!!!!!!!!
          if(with_pc .and. .not.fixed_pc) then
             call communpack(help_arr1,totsym_field_length,1,info)          !!!!!!!!!!!!
             if(info/=0) call error_handler&                                !!!!!!!!!!!!
                  ('gradient_receive_3c: unpacking totalsym_field failed')  !!!!!!!!!!!!
             totalsym_field=totalsym_field+help_arr1                        !!!!!!!!!!!!
          end if                                                            !!!!!!!!!!!!
          if(integralpar_2dervs) then                                       !!!!!!!!!!!!
             do i1=1,to_calc_grads%n_points
                call communpack(help_arr,gradient_data_n_gradients,1,info)
                if(info/=0) call error_handler&
                     ('gradient_receive_3c: unpacking grad_solv_totsym_tes failed')
                grad_solv_totsym_tes(:,i1)=grad_solv_totsym_tes(:,i1)+help_arr
             end do
          end if
       endif                                                                !!!!!!!!!!!!

       if (options_split_gradients()) then
           call communpack(help_arr,gradient_data_n_gradients,1,info)
           if(info/=0) call error_handler&
                ('gradient_receive_3c: unpacking gradient_ob_pulay failed')
           gradient_ob_pulay=gradient_ob_pulay+help_arr
           call communpack(help_arr,gradient_data_n_gradients,1,info)
           if(info/=0) call error_handler&
                ('gradient_receive_3c: unpacking gradient_ch_pulay failed')
           gradient_ch_pulay=gradient_ch_pulay+help_arr
           if (options_xcmode() == xcmode_model_density) then
              call communpack(help_arr,gradient_data_n_gradients,1,info)
              if(info/=0) call error_handler&
                   ('gradient_receive_3c: unpacking gradient_mda_vxc failed')
              gradient_mda_vxc=gradient_mda_vxc+help_arr
              ! gradient_mda_rho is not sent/received !
           endif
       endif

#ifdef WITH_EPE
       if(ewpc_n.ne.0) then
          call communpack (buffer, info)
          if(info/=0) call error_handler&
               ('gradient_receive_3c: unpacking cluster_epe_energy failed')
          cluster_epe_energy = cluster_epe_energy + buffer
       end if
#endif

       if(moving_pc) call totsym_PC_grad_unpack()
       if(moving_X_centers .or. moving_R_centers) call totsym_X_grad_unpack()
       if(moving_Pol_centers) call totsym_id_grad_unpack()
#ifdef WITH_EFP
       if(operations_solvation_effect) then
          if(moving_pc) call totsym_X_solv_grad_unpack()
       end if
#endif

      end if
     end do
    if(what_to_do==deriv .and. allocated(cpks)) then
     deallocate(help_cpks_grad_totalsym,stat=cpksalloc(123))
     ASSERT(cpksalloc(123).eq.0)
     cpksalloc(123)=1
    endif
    DPRINT 'gradient_receive_3c: exit'
  end subroutine gradient_receive_3c

  subroutine gradient_send_3c(what_to_do)
    !  Purpose: send 3 center contribution to the gradient from the
    !           slave to the master
    !** End of interface ****************************************
  use cpksdervs_matrices, only:cpks
  use integralpar_module, only: integralpar_2dervs
  use pointcharge_module, only: moving_pc,totsym_PC_grad_pack
#ifdef WITH_EFP
  use efp_solv_grad_module, only: totsym_X_solv_grad_pack
#endif
    implicit none

    integer :: what_to_do
    !------------ Declaration of local variables ----------------
#ifdef WITH_EPE
    real(r8_kind) :: buffer
#endif
    integer(kind=i4_kind) :: info
    integer(kind=i4_kind) :: i_spn,i
    integer(kind=i4_kind) :: i_g,i_grad,n_occ,n_vac
    integer, parameter :: grads=1, deriv=2

    !------------ Executable code -------------------------------
    DPRINT 'gradient_send_3c: entered'
    call comm_init_send(comm_master_host,msgtag_grad_3c)
    if(what_to_do==grads) then
       call commpack(gradient_totalsym,gradient_data_n_gradients,1,info)
       ASSERT(info.eq.0)
    end if

    if(what_to_do==deriv .and. allocated(cpks)) then

        call cpks_add_fitmat_ch_grads( cpks_gradient_totalsym,cpks_grad_fit_totasym, &
                                       cpks_grad_fit_ch_cartesian) ! used in gradient_send_3c

         call commpack(cpks_gradient_totalsym(1,1),size(cpks_gradient_totalsym),1,info)
        do i_g=1,symmetry_data_n_irreps()
        do i_spn=1,size(cpks,3)
        n_occ=size(cpks(1,i_g,i_spn)%s1,1)
        n_vac=size(cpks(1,i_g,i_spn)%s1ai,2)
        if(n_occ.gt.0) then
        do i_grad=1,size(cpks,1)
        DPRINT  MyId//'cpks(i_grad,i_g,i_spn)%h1ai in send',sum(cpks(i_grad,i_g,i_spn)%h1ai),i_g,i_grad
         call commpack(cpks(i_grad,i_g,i_spn)%s1(1,1),n_occ**2,1,info)
         call commpack(cpks(i_grad,i_g,i_spn)%s1ai(1,1),n_occ*n_vac,1,info)
         call commpack(cpks(i_grad,i_g,i_spn)%h1ai(1,1),n_occ*n_vac,1,info)
        enddo
        endif
        enddo
        enddo
       endif ! cpks

    if(what_to_do==grads) then
    if (operations_solvation_effect) then                               !!!!!!!!!!!!!
       call commpack(grad_solv_totsym,gradient_data_n_gradients,1,info) !!!!!!!!!!!!!
       if(info/=0) call error_handler&                                  !!!!!!!!!!!!!
         ('gradient_send_3c: packing grad_solv_totsym failed')          !!!!!!!!!!!!!
       if(with_pc .and. .not.fixed_pc) then
          call commpack(totalsym_field,totsym_field_length,1,info)
          if(info/=0) call error_handler&
               ('gradient_send_3c: packing totalsym_field failed')
       end if
       if(integralpar_2dervs) then                                       !!!!!!!!!!!!
          do i=1,to_calc_grads%n_points
             call commpack(grad_solv_totsym_tes(1,i),gradient_data_n_gradients,1,info)
             if(info/=0) call error_handler&
                  ('gradient_send_3c: packing grad_solv_totsym_tes failed')
          end do
       end if
    endif                                                               !!!!!!!!!!!!!
    if (options_split_gradients()) then
       call commpack(gradient_ob_pulay,gradient_data_n_gradients,1,info)
       if(info/=0) call error_handler&
            ('gradient_send_3c: packing gradient_ob_pulay failed')
       call commpack(gradient_ch_pulay,gradient_data_n_gradients,1,info)
       if(info/=0) call error_handler&
            ('gradient_send_3c: packing gradient_ch_pulay failed')
       if (options_xcmode() == xcmode_model_density) then
          call commpack(gradient_mda_vxc,gradient_data_n_gradients,1,info)
          if(info/=0) call error_handler&
               ('gradient_send_3c: packing gradient_mda_vxc failed')
          ! gradient_mda_rho is not sent/received !
       endif
    endif

#ifdef WITH_EPE
    if(ewpc_n.ne.0) then
       buffer = cluster_epe_energy(1)
       call commpack (buffer, info)
       ASSERT(info.eq.0)
    endif
#endif

    if(moving_pc) call totsym_PC_grad_pack()
    if(moving_X_centers .or. moving_R_centers) call totsym_X_grad_pack()
    if(moving_Pol_centers) call totsym_id_grad_pack()
#ifdef WITH_EFP
    if(operations_solvation_effect) then
       if(moving_pc) call totsym_X_solv_grad_pack()
    end if
#endif

    end if

    call comm_send()
    DPRINT 'gradient_send_3c: exit'
  end subroutine gradient_send_3c

  subroutine send_receive_Q_grads() !!!!!!!!!!!!!!AS
    !purpose : destribute solvation gradients of charges among
    !          all nodes to calculate solvation second
    !          derivatives
    use comm_module
    use msgtag_module, only: msgtag_grad_3c

    integer :: info

    if( comm_i_am_master() )then
       call comm_init_send(comm_all_other_hosts,msgtag_grad_3c)
       call commpack(Q_grad(1,1),gradient_data_n_gradients*to_calc_grads%n_points,1,info)
       if(info/=0) call error_handler&
            ('send_receive_tes_grads: packing Q_grad failed')
       call comm_send()
    else
       call comm_save_recv(comm_master_host,msgtag_grad_3c)
       call communpack(Q_grad(1,1),gradient_data_n_gradients*to_calc_grads%n_points,1,info)
       if(info/=0) call error_handler&
            ('send_receive_tes_grads: packing Q_grad failed')
    endif

  end subroutine send_receive_Q_grads !!!!!!!!!!!!!AS

#if 0 /* old cpks_add_fitmat_ch_grads */
subroutine  cpks_add_fitmat_ch_grads(grad_totsym,grad_fit)

  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components
  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom
 use unique_atom_module, only : unique_atom_grad_info !It needs for Compac Linux Fortran

 real(kind=r8_kind), intent(inout) :: grad_totsym(:,:,:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:,:,:)  ! fitfct. contributions
 integer(kind=i4_kind) :: i_unique,i_center,index,i,grad_dim,k,l
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
     do k=1,size(grad_totsym,1)
     do l=1,size(grad_totsym,1)
       grad_totsym(k,l,index) = grad_totsym(k,l,index) - &     ! sign changed
            weight * sum( rotmat(i,:) * grad_fit(k,l,:,i_unique) )
     enddo
     enddo
       index=index+1
    enddo
 end do
end subroutine cpks_add_fitmat_ch_grads
#else
#if 0
subroutine  cpks_add_fitmat_ch_grads (add_grad_totsym, grad_totsym, grad_fit)

  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components
  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom
 use unique_atom_module, only : unique_atom_grad_info !It needs for Compac Linux Fortran
 use fit_coeff_module, only: coeff_charge
 implicit none
 real(kind=r8_kind), intent(inout) :: grad_totsym(:,:),add_grad_totsym(:,:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:,:,:)  ! fitfct. contributions
 ! *** end of interface ***

 integer(kind=i4_kind) :: i_unique,i_center,index,i,grad_dim,k,l
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind)    :: totsym_mat(size(grad_totsym,1),size(grad_totsym,1))
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
     totsym_mat=0.0_r8_kind
     do k=1,size(grad_totsym,1)
     do l=1,size(grad_totsym,1)
       totsym_mat(k,l) = totsym_mat(k,l) - &     ! sign changed
            weight * sum( rotmat(i,:) * grad_fit(k,l,:,i_unique) )
     enddo
     enddo
      add_grad_totsym(:,index)=add_grad_totsym(:,index)+ &
                  matmul(coeff_charge,totsym_mat(:,:))
      grad_totsym(:,index)=matmul(totsym_mat(:,:),coeff_charge)
       index=index+1
    enddo
 end do
end subroutine cpks_add_fitmat_ch_grads
#else
subroutine  cpks_add_fitmat_ch_grads (add_grad_totsym, grad_totsym, grad_fit)

  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components
  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom
 use unique_atom_module, only : unique_atom_grad_info !It needs for Compac Linux Fortran
 implicit none
 real(kind=r8_kind), intent(inout) :: grad_totsym(:,:),add_grad_totsym(:,:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:,:)  ! fitfct. contributions
 ! *** end of interface ***

 integer(kind=i4_kind) :: i_unique,i_center,index,i,grad_dim,k!,l
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind)    :: totsym_mat(size(grad_totsym,1))
 real(kind=r8_kind),pointer :: rotmat(:,:)

 do i_unique=1,N_moving_unique_atoms
    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
     totsym_mat=0.0_r8_kind
     do k=1,size(grad_totsym,1)
       totsym_mat(k) = totsym_mat(k) - &     ! sign changed
            weight * sum( rotmat(i,:) * grad_fit(k,:,i_unique) )
     enddo
      add_grad_totsym(:,index)=add_grad_totsym(:,index)+totsym_mat
      grad_totsym(:,index)=totsym_mat
      index=index+1
    enddo
 end do
end subroutine cpks_add_fitmat_ch_grads
#endif
#if 0 /* new code working with one row of grad_fit matrix */
subroutine  cpks_add_fitmat_ch_grads (add_grad_totsym, grad_totsym, grad_fit, i_unique)

  ! purpose: transform cartesian gradients which arise from ch fitfunctions
  !          into symmetry adapted gradient components and add them to
  !          an array of symmetry adapted gradient components
  !    note: fit contributions contain negative gradients -d/dR !
  !          and are only loaded for the first equal atom of each unique atom
 use unique_atom_module, only : unique_atom_grad_info !It needs for Compac Linux Fortran
 use fit_coeff_module, only: coeff_charge
 implicit none
 real(kind=r8_kind), intent(inout) :: grad_totsym(:,:),add_grad_totsym(:,:) ! symm. adapt. grad. array
 real(kind=r8_kind), intent(in)    :: grad_fit(:,:,:)  ! fitfct. contributions
 ! *** end of interface ***

 integer(kind=i4_kind),intent(in) :: i_unique
 integer(kind=i4_kind) :: i_center,index,i,grad_dim,k,l
 real(kind=r8_kind)    :: weight
 real(kind=r8_kind)    :: totsym_mat(size(grad_totsym,1),size(grad_totsym,1))
 real(kind=r8_kind),pointer :: rotmat(:,:)

    i_center=moving_unique_atom_index(i_unique)
    weight=unique_atoms(i_center)%n_equal_atoms
    index=gradient_index(i_unique)
    grad_dim=gradient_index(i_unique+1)-index
    rotmat=>unique_atom_grad_info(i_unique)%m(:,:,1)
    do i=1,grad_dim
     totsym_mat=0.0_r8_kind
     do k=1,size(grad_totsym,1)
     do l=1,size(grad_totsym,1)
       totsym_mat(k,l) = totsym_mat(k,l) - &     ! sign changed
            weight * sum( rotmat(i,:) * grad_fit(k,l,:) )
     enddo
     enddo
      add_grad_totsym(:,index)=add_grad_totsym(:,index)+ &
                  matmul(coeff_charge,totsym_mat(:,:))
      grad_totsym(:,index)=add_grad_totsym(:,index)+matmul(totsym_mat(:,:),coeff_charge)
       index=index+1
    enddo
end subroutine cpks_add_fitmat_ch_grads
#endif
#endif

  subroutine gradient_sndrcv_fit_ch()
    ! to be called from parallel context
    use comm_module
    use msgtag_module, only: msgtag_grad_ch
    implicit none
    ! *** end of interface ***

    if( comm_i_am_master() )then
       call gradient_receive_fit_ch()
    else
       call gradient_send_fit_ch()
    endif
  end subroutine gradient_sndrcv_fit_ch

  subroutine gradient_receive_fit_ch()
    !  Purpose: receive ch fit contribution to the gradient from the
    !           slave
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: n_procs,i,info
    real(kind=r8_kind) :: help_arr(3*N_moving_unique_atoms)
    !------------ Executable code -------------------------------

    n_procs=comm_get_n_processors()
    do i=2,n_procs
       call comm_save_recv(i,msgtag_grad_ch)
       if(comm_msgtag()/=msgtag_grad_ch) call error_handler&
            ('Wrong msgtag in subroutine gradient_receive_fit_ch')

       call communpack(help_arr,3*N_moving_unique_atoms,1,info)
       ASSERT(info.eq.0)
       grad_fit_ch_cartesian=grad_fit_ch_cartesian+&
            reshape(help_arr,(/3,N_moving_unique_atoms/))

       if ( options_xcmode() == xcmode_model_density .and. &
            options_split_gradients()) then
          call communpack(help_arr,3*N_moving_unique_atoms,1,info)
          if(info/=0) call error_handler&
               ('gradient_receive_fit_ch: unpacking grad_mda_rhofit failed')
          grad_mda_rhofit_cartesian=grad_mda_rhofit_cartesian+&
                reshape(help_arr,(/3,N_moving_unique_atoms/))
          call communpack(help_arr,3*N_moving_unique_atoms,1,info)
          if(info/=0) call error_handler&
               ('gradient_receive_fit_ch: unpacking grad_mda_xcpot failed')
          grad_mda_xcpot_cartesian=grad_mda_xcpot_cartesian+&
                reshape(help_arr,(/3,N_moving_unique_atoms/))
       endif
    end do
  end subroutine gradient_receive_fit_ch

  subroutine gradient_send_fit_ch()
    !  Purpose: send 3 center contribution to the gradient from the
    !           slave to the master
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: info

    !------------ Executable code -------------------------------

    call comm_init_send(comm_master_host,msgtag_grad_ch)
    call commpack(reshape(grad_fit_ch_cartesian,&
         (/3*N_moving_unique_atoms/)),3*N_moving_unique_atoms,1,info)
    ASSERT(info.eq.0)

    if ( options_xcmode() == xcmode_model_density .and. &
         options_split_gradients()) then
       call commpack(reshape(grad_mda_rhofit_cartesian,&
            (/3*N_moving_unique_atoms/)),3*N_moving_unique_atoms,1,info)
       if(info/=0) call error_handler&
            ('gradient_send_fit_ch: packing grad_mda_rhofit_cartesian failed')
       call commpack(reshape(grad_mda_xcpot_cartesian,&
            (/3*N_moving_unique_atoms/)),3*N_moving_unique_atoms,1,info)
       if(info/=0) call error_handler&
            ('gradient_send_fit_ch: packing grad_mda_xcpot_cartesian failed')
    endif
    call comm_send()

  end subroutine gradient_send_fit_ch
  !**************************************************************

  subroutine freq(k,m,r,w,mef,ang,rms,u,v,intensity)
    !
    ! Solves the eigenvalue equation:
    !
    !     ( K - w2 * M ) X = 0
    !
    ! with force matrix K and mass matrix M
    !
    ! Outputs result on "iou"
    !
    use matrix_module, only: eigs
    use efield_module, only: efield_der, efield_intensity
    implicit none
    real(r8_kind)   , intent(in)    :: m(:)     ! (3*NA)
    real(r8_kind)   , intent(in)    :: k(:,:)   ! (3*NA,3*NA)
    real(r8_kind)   , intent(in)    :: r(:)     ! (3*NA)
    real(r8_kind)   , intent(out)   :: w(:)     ! frequency
    real(r8_kind)   , intent(out)   :: mef(:)   ! ``effective'' mass
    real(r8_kind)   , intent(out)   :: ang(:,:) ! rotation    speed
    real(r8_kind)   , intent(out)   :: rms(:,:) ! mass center speed
    real(r8_kind)   , intent(out)   :: u(:,:)   ! modes  in m-frame
    real(r8_kind)   , intent(out)   :: v(:,:)
    real(r8_kind)   , intent(out)   :: intensity(:)
    optional                        :: r,intensity
    ! *** end of interface ***

    integer(i4_kind) :: i,a,j,i_mode,N
    real(r8_kind)    :: h(size(m),size(m))
    real(r8_kind)    :: nrm(size(m))             ! lenght-norm sqare of mode vectors
    real(r8_kind)    :: rm(3), mtot              ! mass center, mass
    real(r8_kind)    :: itot(3)                  ! inertia moments
    real(r8_kind)    :: x(size(m))               ! coords in m-frame
    real(r8_kind)    :: der(3,3,size(m)/3)       ! double derivative of energy in cartesian coordinates
    real(r8_kind)    :: der_q(3,size(m))         ! double derivative of energy in normal mode
    logical          :: have_coordinates
    logical          :: error

    !------------ Executable code -------------------------------

    have_coordinates = present(r)


    ! ``orthonormalize'' generalized eigenvalue problem:
    do i=1,size(m)
         h(:,i) = k(:,i) / sqrt(m(i))
    enddo
    do i=1,size(m)
         h(i,:) = h(i,:) / sqrt(m(i))
    enddo

    ! solve eigenvalue problem:
    call eigs(h,w,v)

    ! square roots of w^2, negative==imaginary:
    where( w > 0.0_r8_kind )
       w = sqrt(w)
    elsewhere
       w = - sqrt(-w)
    endwhere

    ! back-transform eigenvectors:
    do i=1,size(m)
      v(i,:) = ( 1 / sqrt(m(i)) ) * v(i,:)
    enddo


    ! for IR intensity calculation
    if ( efield_intensity() ) then
      call efield_der(der,error)

      do i_mode = 1, size(m)
         der_q = 0.0
         do i = 1, 3              ! cartesian component of electric field
           do N = 1, size(m)/3    ! atomic indices
             do j = 1, 3          ! cartesian component of atomic position
                der_q(i,i_mode) = der_q(i,i_mode) + der(i,j,N) * v((N-1)*3+j,i_mode)
             enddo
           enddo
         enddo
         intensity(i_mode) = sum(der_q(:,i_mode)**2)
       enddo

    endif

    ! evaluate their cartesian length-norm:
    do i=1,size(m)
      nrm(i) = dot_product( v(:,i), v(:,i) )
    enddo

    ! and re-norm mode vectors to unit length:
    do i=1,size(m)
      v(:,i) = v(:,i) / sqrt(nrm(i))
    enddo

    ! compute the effective masses (only diagonal):
    do i=1,size(m)
      mef(i) = dot_product( v(:,i), m(:)*v(:,i) )
    enddo

    ! compute the mass center ``velocity'' for the mode:

    mtot = sum(m)/3 ! total mass

    do i=1,size(m)
      rms(:,i) = 0.0_r8_kind
      do a=0,size(m)/3 - 1
        j = 3*a+1 ! take e.g. m(x) of atom a
        rms(1,i) = rms(1,i) + m(j) * v(3*a+1,i)
        rms(2,i) = rms(2,i) + m(j) * v(3*a+2,i)
        rms(3,i) = rms(3,i) + m(j) * v(3*a+3,i)
      enddo
      rms(:,i) = rms(:,i) / mtot
    enddo

    ! subtract translation part of the mode:
    do i=1,size(m)
      do a=0,size(m)/3 - 1
        u(3*a+1,i) = v(3*a+1,i) - rms(1,i)
        u(3*a+2,i) = v(3*a+2,i) - rms(2,i)
        u(3*a+3,i) = v(3*a+3,i) - rms(3,i)
      enddo
    enddo

    ![[=== otional, if coordinates are available ===
    if( have_coordinates )then

      mtot = sum(m)/3 ! total mass

      ! coordinates specified => compute mass center
      rm(:) = 0.0_r8_kind
      do a=0,size(m)/3 - 1
        j = 3*a+1 ! take e.g. m(x) of atom a
        rm(1) = rm(1) + m(j) * r(3*a+1)
        rm(2) = rm(2) + m(j) * r(3*a+2)
        rm(3) = rm(3) + m(j) * r(3*a+3)
      enddo
      rm(:) = rm(:) / mtot

      ! go to center mass frame:
      do a=0,size(m)/3 - 1
        x(3*a+1) = r(3*a+1) - rm(1)
        x(3*a+2) = r(3*a+2) - rm(2)
        x(3*a+3) = r(3*a+3) - rm(3)
      enddo

      ! compute inertia moments:
      itot(:) = 0.0
      do a=0,size(m)/3 - 1
        j = 3*a+1 ! take e.g. m(x) of atom a
        itot(3) = itot(3) + m(j) * x(3*a+1)**2 &
                          + m(j) * x(3*a+2)**2
        itot(1) = itot(1) + m(j) * x(3*a+2)**2 &
                          + m(j) * x(3*a+3)**2
        itot(2) = itot(2) + m(j) * x(3*a+3)**2 &
                          + m(j) * x(3*a+1)**2
      enddo

      ! compute angular velocity of mode:
      do i=1,size(m)
        ang(:,i) = 0.0_r8_kind
        do a=0,size(m)/3 - 1
          j = 3*a+1 ! take e.g. m(x) of atom a
          ang(3,i) = ang(3,i) + m(j) * x(3*a+1) * u(3*a+2,i) &
                              - m(j) * x(3*a+2) * u(3*a+1,i)
          ang(1,i) = ang(1,i) + m(j) * x(3*a+2) * u(3*a+3,i) &
                              - m(j) * x(3*a+3) * u(3*a+2,i)
          ang(2,i) = ang(2,i) + m(j) * x(3*a+3) * u(3*a+1,i) &
                              - m(j) * x(3*a+1) * u(3*a+3,i)
        enddo
      enddo

      do i=1,3
        if( abs(itot(i)) > 1.0D-7 )then
          ang(i,:) = ang(i,:) / itot(i)
        else
          ang(i,:) = 0.0_r8_kind
        endif
      enddo

      ! subtract rotation part of the mode:
      do i=1,size(m)
        do a=0,size(m)/3 - 1
          u(3*a+3,i) = u(3*a+3,i) - ang(1,i) * x(3*a+2) &
                                  + ang(2,i) * x(3*a+1)
          u(3*a+1,i) = u(3*a+1,i) - ang(2,i) * x(3*a+3) &
                                  + ang(3,i) * x(3*a+2)
          u(3*a+2,i) = u(3*a+2,i) - ang(3,i) * x(3*a+1) &
                                  + ang(1,i) * x(3*a+3)
        enddo
      enddo
    endif
    !]]=============================================

  end subroutine freq
!***************************************************************

!**************************************************************
  subroutine freq_print(iou,m,r,w,mef,ang,rms,u,v,intensity)
    use atom_data_module, only: nuc_mass!(1&6)
    use efield_module, only: efield_intensity
    implicit none
    integer(i4_kind), intent(in)   :: iou
    real(r8_kind)   , intent(in)   :: m(:)           ! (3*NA)
    real(r8_kind)   , intent(in)   :: r(:)           ! (3*NA)
    real(r8_kind)   , intent(in)   :: w(:)           ! frequency
    real(r8_kind)   , intent(in)   :: mef(:)         ! ``effective'' mass
    real(r8_kind)   , intent(in)   :: ang(:,:)       ! rotation    speed
    real(r8_kind)   , intent(in)   :: rms(:,:)       ! mass center speed
    real(r8_kind)   , intent(in)   :: u(:,:)         ! modes  in m-frame
    real(r8_kind)   , intent(in)   :: v(:,:)         ! modes  in m-frame
    real(r8_kind)   , intent(in)   :: intensity(:)   ! modes  in m-frame
    optional                       :: r, intensity

    ! *** end of interface ***
    real(r8_kind), parameter   :: au2ev = 27.211658_r8_kind    ! eV/au
    real(r8_kind), parameter   :: ev2cm = 8065.5410_r8_kind    ! cm-1/eV
    real(r8_kind), parameter   :: au2cm = au2ev * ev2cm        ! cm-1/au
    real(r8_kind), parameter   :: mp = 1836.15267261_r8_kind   ! Mp/Me
    real(r8_kind), parameter   :: au2D  = 2.541584_r8_kind     ! debye/au
    real(r8_kind), parameter   :: au2A  = 0.529177_r8_kind     ! angstram/aueV
    real(r8_kind), parameter   :: D_A   = 23.071130228_r8_kind ! (au2D/au2A)^2
    ! for conversion between different unit of Intensity, please follow JCP 84, 2262 (1986)
    real(r8_kind)              :: m1
    integer(i4_kind)           :: i,a
    logical                    :: have_coordinates

    !------------ Executable code -------------------------------

    have_coordinates = present(r)

    ! atomic mass unit in (electronic) a.u. (~1836):
    m1 = (nuc_mass(6)/12) * (mp/nuc_mass(1))


  if( iou > 0 )then
      write(iou,'("SD:")')
      write(iou,'("SD: MASS MATRIX:")')
      write(iou,'("SD:",A4,2A16)') "I","MASS [C12]","MASS [Me]"
      do a=0,size(m)/3 - 1
      write(iou,'("SD:",I4,2F16.6)') a+1,m(3*a+1)/m1,m(3*a+1)
      enddo

      write(iou,'("SD:")')
      write(iou,'("SD: FREQUENCIES:")')
      write(iou,'("SD:",A4,4A16)') "I","OMEGA [au]"  ,"MEFF [au]"  &
                                      ,"OMEGA [cm-1]","MEFF [C12]"
      do i=size(w),1,-1
      write(iou,'("SD:",I4,F16.6,F16.1,2F16.3)') i, w(i)      , mef(i)    &
                                                  , w(i)*au2cm, mef(i)/m1
      enddo

      write(iou,'("SD:")')
      write(iou,'("SD: MODE VECTORS (IN M-FRAME):")')
      write(iou,'("SD:",A4,A16,A38)') "I","OMEGA [cm-1]","X(:,I) [length units]"
      do i=size(w),1,-1
        write(iou,'("SD:",I4,F16.6,3F16.6,/:("SD:",20X,3F16.6))') &
             i, w(i)*au2cm, v(:,i)
        write(iou,'("SD:",I4,F16.6,3F16.6,/:("SD:",20X,3F16.6))') &
             i, w(i)*au2cm, u(:,i)
        write(iou,'("SD:",4X,A16  ,3F16.6)') "m-frame speed", rms(:,i)
        if( have_coordinates )then
        write(iou,'("SD:",4X,A16  ,3F16.6)') "angular speed", ang(:,i)
        endif
        write(iou,'("SD:")')
      enddo
      if (efield_intensity()) then
        write(iou,'("SD:",A4,5A20)') "I","OMEGA [au]","OMEGA [cm-1]", &
                                  "[(D/A)^2amu^-1]", "[Km/mol]", "[cm^-2atm^-1]"
        do i=1,size(w)
        write(iou,'("SD:",I4,5F20.12)') i, w(i), w(i)*au2cm, &
           intensity(i)*D_A*1000, intensity(i)*974.8706077*1000, &
           intensity(i)*3960.159*1000
       ! 1000 is for conversion of kilometer ==> meter
      enddo
      endif
    endif
  end subroutine freq_print
!**********************************************************************

  subroutine cart_coor_vec (uas, r)
    use unique_atom_module, only: uat=>unique_atom_type
    implicit none
    type(uat), intent(in) :: uas(:)    ! (NUA)
    real(r8_kind), intent(out) :: r(:) ! (3*NA)
    ! *** end of interface ***

    integer(i4_kind) :: u,e,a,i

    ! construct linear vector of coordinates:
    a = 0
    do u=1,size(uas)
    do e=1,uas(u)%n_equal_atoms
       a = a + 1
       do i=1,3
         r(i+3*(a-1)) = uas(u)%position(i,e)
       enddo
    enddo
    enddo
  end subroutine cart_coor_vec


  subroutine cart_mass_mat(iou,uas,m)
    use unique_atom_module, only: uat=>unique_atom_type
    use atom_data_module, only: nuc_mass
    implicit none
    integer(i4_kind), intent(in)  :: iou
    type(uat)       , intent(in)  :: uas(:) ! (NUA)
    real(r8_kind)   , intent(out) :: m(:)   ! (3*NA)
    ! *** end of interface ***

    integer(i4_kind)           :: u,e,a,i
    real(r8_kind)              :: mc,mm,z
    real(r8_kind), parameter   :: mp = 1836.15267261_r8_kind ! Mp/Me

    ! construct (diagonal) mass matrix:
    if( iou > 0 )then
       write(iou,'("SD:")')
       write(iou,'("SD: MASS MATRIX:")')
       write(iou,'("SD:",A4,A6,A20,A20)') "UA","Z(UA)","MASS(UA) [C12]","MASS(UA) [Me]"
    endif
    a = 0
    do u=1,size(uas)
       z  = uas(u)%Z
       mc = nuc_mass( NINT(z) )
       mm = mc * ( mp / nuc_mass(1) )
       if( iou > 0 )then
          write(iou,'("SD:",I4,F6.2,2F20.12)') u, z, mc, mm
       endif
       do e=1,uas(u)%n_equal_atoms
          a = a + 1
          do i=1,3
            m(i+3*(a-1)) = mm
          enddo
       enddo
    enddo
  end subroutine cart_mass_mat

  subroutine optimizer_cart_mass_mat(iou,charge,na,m)
    use unique_atom_module, only: uat=>unique_atom_type
    use atom_data_module, only: nuc_mass
    implicit none
    integer(i4_kind), intent(in)  :: iou,na
    real(r8_kind)       , intent(in)  :: charge(:) ! (NUA)
    real(r8_kind)   , intent(out) :: m(:)   ! (3*NA)
    ! *** end of interface ***

    integer(i4_kind)           :: u,a,i!,e
    real(r8_kind)              :: mc,mm
    real(r8_kind), parameter   :: mp = 1836.15267261_r8_kind ! Mp/Me

    ! construct (diagonal) mass matrix:
    if( iou > 0 )then
       write(iou,'("SD:")')
       write(iou,'("SD: MASS MATRIX:")')
       write(iou,'("SD:",A4,A6,A20,A20)') "UA","Z(UA)","MASS(UA) [C12]","MASS(UA) [Me]"
    endif
    a = 0
    do u=1,na
       mc = nuc_mass( NINT(charge(u)) )
       mm = mc * ( mp / nuc_mass(1) )
       if( iou > 0 )then
          write(iou,'("SD:",I4,F6.2,2F20.12)') u, charge(u), mc, mm
       endif
       a = a + 1
       do i=1,3
         m(i+3*(a-1)) = mm
       enddo
    enddo
  end subroutine optimizer_cart_mass_mat

  subroutine cart_hess_mat(iou,h,hh)
    implicit none
    integer(i4_kind), intent(in)  :: iou
    type(arrmat4)   , intent(in)  :: h(:,:)  ! (NUA,NUA)
    real(r8_kind)   , intent(out) :: hh(:,:) ! (3*NA,3*NA)
    ! *** end of interface ***

    integer(i4_kind) :: u1,u2,e1,e2,a1,a2
    integer(i4_kind) :: i,j
    integer(i4_kind) :: n2,n1
    character(len=1), parameter  :: x(3) = (/ 'X', 'Y', 'Z' /)

    n2 = size(hh,1) / 3
    n1 = size(hh,2) / 3
    ASSERT(n1==n2)

    if( iou > 0 )then
    write(iou,'("SD: "," ========================= CARTESIAN HESSIAN ===================================")')
    write(iou,'("SD: ",4A4,4X,3A20)') 'U1', 'E1', 'U2', 'E2', (x(i), i=1,3)
    endif

    a2 = 0
    do u2=1,size(h,2)
      do e2=1,size(h(1,u2)%m,4)
         a2 = a2 + 1
         ASSERT(a2<=n2)

    a1 = 0
    do u1=1,size(h,1)
      do e1=1,size(h(u1,1)%m,2)
        a1 = a1 + 1
        ASSERT(a1<=n1)

        do j=1,3
        do i=1,3
           hh(i+3*(a1-1),j+3*(a2-1)) = h(u1,u2)%m(i,e1,j,e2)
        enddo
        enddo

        if( iou > 0 )then
        write(iou,'("SD: ",4I4,A4,3F20.12)') &
                   u1, e1, u2, e2, &
                   x(1), ( h(u1,u2)%m(1,e1,j,e2), j=1,3 )
        write(iou,'("SD: ",16X,A4,3F20.12)') &
                   x(2), ( h(u1,u2)%m(2,e1,j,e2), j=1,3 )
        write(iou,'("SD: ",16X,A4,3F20.12)') &
                   x(3), ( h(u1,u2)%m(3,e1,j,e2), j=1,3 )
        ! FIXME: replace the three statements with one (check Intel):
!       write(iou,'("SD: ",4I4,A4,3F20.12,/:("SD: ",16X,A4,3F20.12))') &
!                  u1, e1, u2, e2, &
!                  ( ( x(i), ( h(u1,u2)%m(i,e1,j,e2), j=1,3 ) ), i=1,3 )
        endif
    enddo
    enddo
    enddo
    enddo

    if( iou > 0 )then
    write(iou,'("SD: "," ===============================================================================")')
    endif
  end subroutine cart_hess_mat

end module gradient_data_module

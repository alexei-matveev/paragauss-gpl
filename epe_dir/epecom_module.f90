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
module epecom_module
!
!== Interrupt of public interface of module ====================
!
!---------------------------------------------------------------
! Modifications
!---------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!--------------------------------------------------------------
!== Interrupt end of public interface of module ===============

# include "def.h"
  use type_module
  use filename_module, only: max_path => filename_namelengthmax
  implicit none
  save
  public
  real(kind=r8_kind):: cross_boundary_3b=0.0_r8_kind
  integer(kind=i4_kind):: epealloc_stat(26)=0
  logical::pc_aswritemode=.false.
  real(kind=r8_kind) :: coordinates_moving_epecenter(3)
  integer(kind=i4_kind):: moving_epecenter=0

  logical:: fixed_reg2a_relcontribs=.true.
  logical:: ml_displacements_fixed=.false.
  logical:: explicit_coupling=.false.
  logical:: explicit_coupling_3body=.true.
  logical :: reg_2a_treated=.false.
  integer(kind=i4_kind), public :: N_VACANCIES=0,N_IMPURITIES=0
        logical:: dg_convergence_reached=.false.
  logical :: ml_displacements=.false.
  logical :: GX_HIGHPREC = .true.
  integer(kind=i4_kind):: gxcell_unit,output_epe,input_epe,eperef_unit
  logical :: relax_shells_only=.false.

  integer(kind=i4_kind), parameter :: ndr1=2500 &   ! limit for   size of region I
                                   ,ndpt=200 &   ! posible number of impurtities
                                   ,ndnc=210 &   ! possible No of generation centers
                                   ,ndt=20 &    ! possible No of ion types in the cell
                                   ,ndcell=400 & ! possible No of ions in cell
                                   ,ndngv=48000 & ! possible No of reciprocial translations
                                   ,ngxat=300 &  ! No of atoms in gx file
                                   ,n_hds_cycles=6 &
                                   ,nmadt=300 &  ! number of types with respect to field
                                   ,klgen=250    ! number of reg. 1 ions for init  gener
  integer(kind=i4_kind)  :: n_grads_master
  integer(kind=i4_kind)  :: n_grads_slave
  integer(kind=i4_kind)  :: n_grads_total
  integer(kind=i4_kind)  :: end_treated_region

  real(kind=r8_kind), public :: DELTA
  real(kind=r8_kind)::epe_side_optimized_energy
  real(kind=r8_kind)::epe_side_optimized_energy_prev

logical :: ml_tensors=.false.
logical :: independent_epe_charges=.false.
logical :: pseudobond_approach=.false.
  logical :: ml_cluster_simulated=.false.
  integer(kind=i4_kind)  ::n_ml_cluster_simulators=0
  type cluster_simulator
     real(kind=r8_kind)::r(3)
     real(kind=r8_kind)::q
  end type cluster_simulator
  type(cluster_simulator), allocatable:: ml_cluster_simulators(:)

  real(kind=r8_kind)::scaling_factor
  type pg_field_at
     real(kind=r8_kind)::vs
     real(kind=r8_kind), dimension(3)::gs
     real(kind=r8_kind)::vc
     real(kind=r8_kind), dimension(3)::gc
     real(kind=r8_kind), dimension(3)::rs
     real(kind=r8_kind), dimension(3)::rc
  end type pg_field_at
  type (pg_field_at), allocatable, dimension(:) :: reg_I_pg

  real(kind=r8_kind), dimension(3,3):: rot_gto_to_epe,gmat
  real(kind=r8_kind) :: E2BIND,E2BDIS,E2BDEF,ECRR,shft_gto_to_epe(3)=(/0.,0.,0./)

  real(kind=r8_kind) :: shft_gto_to_hds(3) ! duumy input for consistency
  integer(kind=i4_kind), dimension(ndpt)::      TYPE_IMPURITY
  integer(kind=i4_kind), dimension(ndpt)::      imp2center
  integer(kind=i4_kind), dimension(ndpt)::      explicit_coupling_type
        ! to be used in calculations of var and ref configuration
  integer(kind=i4_kind), dimension(:),allocatable:: ml_dprime_factors
  integer(kind=i4_kind) :: ihds,nucen_pc,nmvmax,max_type_ions
  logical :: periodic_optimization=.false.      ! make global oprimization
  logical :: dealloc_epe_ref
  logical :: use_pgdata=.true.,use_epe_pgdata
  logical :: read_configuration=.true.,write_configuration=.true., &
           write_config_to_lcgto_format,print_deform, &
           print_gradients=.false., &
           print_ions_from_first_sphere=.false.
  logical :: option_c3_symm=.false.
  real(kind=r8_kind):: epeside_energy_limit = 0.00005_r8_kind
  logical :: embed_convergence_check = .true.
  logical,parameter:: df_embed_convergence_check= .true.
  real(kind=r8_kind):: embed_convergence_limit = 0.00001_r8_kind
  real(kind=r8_kind),parameter:: df_embed_convergence_limit= 0.00001_r8_kind

  type symm_operations
     real(kind=r8_kind), dimension(3,3)::rotmat
  end type symm_operations
  type(symm_operations),dimension(:),allocatable::c3
  integer(kind=i4_kind),dimension(:,:),allocatable:: permut_ind

  type ep_enviroment
     real(kind=r8_kind),dimension(3)::r,s,c,o
     integer(kind=i4_kind) ::l,u,k,m,gc,sym_st_ind,ec
     real(kind=r8_kind) ::d,q,qs,qc
     logical :: fix
  end type ep_enviroment
  type(ep_enviroment), allocatable,target, dimension(:)::epe
  type(ep_enviroment):: varepe
  type(ep_enviroment),allocatable,dimension(:) :: epe_at_start

  real(kind=r8_kind),dimension(:),allocatable:: Q_zl

  integer(kind=i4_kind) :: n_centres_of_generation
  real(kind=r8_kind), dimension(ndnc,3)::r_cent_gener

  real(kind=r8_kind), allocatable, dimension(:,:):: R_NUC_ION,r_sh_ion

  type short_range_par
     real(kind=r8_kind), pointer, dimension(:,:,:)::b,ro,c,d,sr1,sr2,k,r0,k1,r1 !!!!!!!!!!!!!AS
  endtype short_range_par
  type (short_range_par)::host,host_tmp

  real(kind=r8_kind), dimension(3,3):: VECTORS_TRANS
  real(kind=r8_kind), dimension(ndcell,3):: R_ION_IN_CELL,R_ION_IN_CELL_OLD
  type epe_ions
     real(kind=r8_kind) :: r_core(3),r_shell(3)
  end type epe_ions
  type(epe_ions), allocatable,dimension(:) :: regI_previous
  logical, dimension(ndcell)::fixed
  integer(kind=i4_kind) :: n_fixed
  real(kind=r8_kind), dimension(ndcell,3):: R_SHELL_IN_CELL,R_CORE_IN_CELL
  real(kind=r8_kind), allocatable,target,dimension(:,:) :: point_1_core_grad, &
                                                           point_1_shell_grad,&
                                                           point_0_core_grad, &
                                                           point_0_shell_grad
  real(kind=r8_kind), dimension(3,ndcell):: gpsi ! grad of reciprocial contributions
  real(kind=r8_kind), dimension(ndcell):: Q_z
  integer(kind=i4_kind), dimension(ndcell)::type_of_ion
  integer(i4_kind) :: n_ls
  logical :: reset, use_lin_search

  type ref_to_epe
     integer(kind=i4_kind)::new,old
  end type ref_to_epe
  type(ref_to_epe),allocatable,dimension(:):: which_epe_ion

  integer(kind=i4_kind) :: n_trans_primitive_cell(3),n_ions_cell
  real(kind=r8_kind), dimension(ndngv,3):: GSTR
  real(kind=r8_kind), dimension(ndngv):: rsin,rcos
  integer(kind=i4_kind) :: n_bs_points

  type madi
     real(kind=r8_kind),dimension(nmadt):: EMD,FMD,AMD,RMD
  end type madi
  type(madi) :: mad,madc
  type ml_tensors_type
   real(kind=r8_kind),dimension(3,3)::c,s
  end type ml_tensors_type
  type (ml_tensors_type),dimension(:), allocatable:: ml_ten
  type ml_factors
     real(kind=r8_kind)::s,c
  end type ml_factors
  type(ml_factors), dimension(:), allocatable:: ml_fac
  real(kind=r8_kind), dimension(ndt):: Q_NUCLEAR,q_shell,pk,q_ion

  real(kind=r8_kind), dimension(ndt)::q_epecl_coupl
  logical :: qm_interfaced_mode=.false.
  integer(kind=i4_kind) :: &
                        n_gen_ions & ! limit in number of ions in the dir. latt. sum
                        ,n_print_additional_ions, &
                         in2
  integer(kind=i4_kind) :: basic_action
! 0-no epe displacements, 1- lattice displacements,2-displacements of impurities
  real(kind=r8_kind) :: del,dmax,reca

  real(kind=r8_kind), parameter :: qau_qepe=3.794813_r8_kind, auangs=.529177_r8_kind,zero=0.0_r8_kind,two=2.0_r8_kind
  real(kind=r8_kind), parameter :: pi=3.14159265_r8_kind, PI2=two*pi, eau_ev=27.211652_r8_kind

  real(kind=r8_kind), parameter :: df_radius_first_sphere=5.0_r8_kind, df_radius_2A_sphere=15.0_r8_kind
  real(kind=r8_kind) :: radius_first_sphere, radius_2A_sphere

  real(kind=r8_kind) :: qpivc,pta2,vc
  real(kind=r8_kind) :: ET2,ERFO,PIS

  real(kind=r8_kind),parameter :: df_ERROR_FUNCTION_PARAMETER=1.9177190_r8_kind
  real(kind=r8_kind) ::  ERROR_FUNCTION_PARAMETER

  real(kind=r8_kind) :: sescp,cescp,etot_epe,an(ngxat),vshf(ngxat)
  real(kind=r8_kind), allocatable,dimension(:,:) :: gx_imp_r ! reg position of gx impurities
  real(kind=r8_kind), dimension(3,ngxat):: gx_centers  ! positons of gx file
  integer(kind=i4_kind) :: jps(ngxat),indgx(8,ngxat),impu(ngxat)
  integer(kind=i4_kind) :: natoms,isym

  real(kind=r8_kind), dimension(ndpt,3) ::R_NUC_IMP,R_SH_IMP,R_IMP,R_NUC_IMPO,R_SH_IMPO,RPO,dir

  real(kind=r8_kind),parameter :: df_radius_long_interact=100.0_r8_kind, &
  df_radius_first_short_interact=10.0, &
  df_radius_second_short_interact=10.0
  real(kind=r8_kind) :: radius_long_interact


  real(kind=r8_kind), dimension(ndpt)::Q_NUC_IMPURITY,Q_SH_IMPURITY,Q_IMPURITY,pk_impurity
  real(kind=r8_kind), dimension(ndpt,ndpt):: B_IM,RO_IM,C_IM,D_IM
  character (len=max_path), public :: epedata_dir, epe_input_dir
  integer(kind=i4_kind),public:: reg_I_n_ions,n_pgepe_iterations=95,reg_2a_n_ions
  integer(kind=i4_kind), public :: epe_iter, epeit_count
  logical::pg_cfield_calculated

  integer(kind=i4_kind) ::pg_interfaced_mode=0,i_ir_eperef
  integer(kind=i4_kind),parameter:: df_pg_interfaced_mode=0
  real(kind=r8_kind) e_epe_atstart,e_epe_final
  logical, parameter:: df_epe_rel_converged=.false.
  logical:: epe_rel_converged=.false.
  logical:: make_epe_reference=.false.
  logical, parameter:: df_make_epe_reference=.false.
  logical:: use_epe_reference=.false.
  logical,parameter:: df_use_epe_reference=.false.

  logical:: extended_epe=.false.
  logical,parameter:: df_extended_epe=.false.

  logical:: ex_pgdata=.false.
  logical,parameter:: df_ex_pgdata=.false.

  logical:: fixed_dielectric_const=.false.
  logical,parameter:: df_fixed_dielectric_const=.false.

  real(kind=r8_kind) :: fx_dielectric_const= 1.d0
  real(kind=r8_kind),parameter:: df_fx_dielectric_const= 1.d0

  type,public:: epe_centers
     real(kind=r8_kind),dimension(3):: rs,rc
  end type epe_centers
  type(epe_centers),allocatable,dimension(:),public :: reg_reference,epe_reference
  type,public::impurity
     real(kind=r8_kind):: gc(3),gs(3),rc(3),rs(3)
  end type impurity
  type(impurity), allocatable,dimension(:),target,public:: impurities_new_conf
  type(impurity), allocatable,dimension(:),target,public:: impurities_old_conf
  type(impurity), pointer,dimension(:),public:: imp_conf
  real(kind=r8_kind):: defect_energy_coul=0.0_r8_kind,defect_energy_short=0.0_r8_kind &
                      ,ecoul_epecluster,eshort_epecluster &
                      ,ecoul_vaccluster,eshort_vaccluster &
                      ,ec_cluster_reg2a,epg_cluster_reg_I,eshort_reg_I &
                      ,eshort_reg2a,diffpg_ec_ecref,diff_reg2a_ec_ecref &
                      ,ecoul_lat_relaxation,ecoul_lat_relcorr,eshort_lat_relcorr &
                      ,defect_contrib=0.0_r8_kind
  real(kind=r8_kind),parameter :: df_defect_contrib=0.0_r8_kind
  integer(kind=i4_kind)::rel_converged_unit
  real(kind=r8_kind):: var_epe_contrib_atstart

! public variables to calculate 3 body interaction
  integer(kind=i4_kind) :: n_types_central_atoms_3body = 0
  integer(kind=i4_kind) :: n_ec_types_central_atoms_3body = 0
  real(kind=r8_kind), allocatable :: ki(:,:,:),theta_0(:,:,:),r3b(:),r3b_im(:)

  type, public:: expicite_coupling_3body
   real(kind=r8_kind), pointer :: ki(:,:,:),theta_0(:,:,:),r3b(:)
   integer(kind=i4_kind),pointer :: types(:,:),n_3b(:),index(:,:,:)
  end type expicite_coupling_3body
  type(expicite_coupling_3body)::ec


  integer(kind=i4_kind), allocatable :: types(:,:)
  integer(kind=i4_kind) :: n_types_central_atoms_3body_im = 0
  integer(kind=i4_kind), allocatable :: types_im(:,:)
  integer(kind=i4_kind) :: n_tetrahedrons=0, first_ind, last_ind
  integer(kind=i4_kind), allocatable :: tetra_atoms(:,:)

! public variables to perform periodic optimization and epe relaxation
  real(kind=r8_kind), allocatable :: H_m(:,:),H_m1(:,:)
  real(kind=r8_kind) :: etot_epe_0,etotal_gopt_0, r_max=0.5_r8_kind
  real(kind=r8_kind) :: abs_g=0.00002_r8_kind
  integer(kind=i4_kind) :: n_iterations=200, n_hess_update_cycles=20
  real(kind=r8_kind) :: weight_hess=1.0_r8_kind

! public variables what cantrol output of XYZ files
  logical :: unit_cell_output=.false. ,region_I_output=.false.
  integer :: output_step=1
  character(len=3), dimension(ndt):: name_of_type
  character(len=3), dimension(ndcell):: name_of_ions


! this variable informs program that unit cell was taken from GULP calculation
  logical :: core_shell


! variables controlling a calculation of electrostatic potential
  logical :: lpotcalc
  real(kind=r8_kind) :: point00(3),point01(3),point02(3),dxaxis(2),dyaxis(2),v_abs_limit
  integer(kind=i4_kind) :: nxgrid,nygrid
  character(len=4) :: output_units
  character(len=9) :: output_format

! this variables redefines epe charges to use them in qm-epe calculation
  real(kind=r8_kind) :: out_charges(ndt)

  real(kind=r8_kind) :: scale_factor=1.0_r8_kind

  logical :: do_print=.true.
contains

 subroutine get_epe_defaults()
   embed_convergence_limit=df_embed_convergence_limit
   pg_interfaced_mode=df_pg_interfaced_mode
   defect_contrib=df_defect_contrib
   fx_dielectric_const  =df_fx_dielectric_const
   fixed_dielectric_const=df_fixed_dielectric_const
   ex_pgdata= df_ex_pgdata
   extended_epe= df_extended_epe
   make_epe_reference=df_make_epe_reference
   use_epe_reference=df_use_epe_reference
   epe_rel_converged= df_epe_rel_converged
   radius_first_sphere=df_radius_first_sphere
   radius_2A_sphere=df_radius_2A_sphere
   ERROR_FUNCTION_PARAMETER=df_ERROR_FUNCTION_PARAMETER
   radius_long_interact=df_radius_long_interact
   embed_convergence_check= df_embed_convergence_check
 end subroutine get_epe_defaults

 subroutine get_epe_energies &
      (diff_ec_ecref,etot_epe_corrected,lattice_energy,epg_cluster_reg_I,&
      var_epe_contrib,eshort_reg_Iau,eshort_coupling_au )
   !** purpes: gives epe energies in AU

   use filename_module
   use iounitadmin_module
   use time_module
   use timer_module

   real(kind=r8_kind),intent(out),optional:: diff_ec_ecref,etot_epe_corrected &
                  ,lattice_energy,epg_cluster_reg_I,var_epe_contrib &
                  ,eshort_reg_Iau, eshort_coupling_au

   call start_timer(timer_get_epe_energy)

   if(present(eshort_reg_Iau)) eshort_reg_Iau=eshort_reg_I/eau_ev

   if(present(diff_ec_ecref)) &
        diff_ec_ecref=(diffpg_ec_ecref+diff_reg2a_ec_ecref)/eau_ev
   if(present(var_epe_contrib)) then
      var_epe_contrib=(etot_epe-ec_cluster_reg2a)/eau_ev
   endif

   if(present(etot_epe_corrected)) then
        !(1)  Intra-cluster interaction never goes in etot_epe except
        !     for regular reference contribution. Thus  it is constant on EPE side
        !     and on QM side it is present in SCF energy.
        !     No correction for this interaction is required.
        !(2)  Cluster  is linked to ions of region I. Such contribution is not
        !     present in PG SCF energy but it is treated with optimizer.
        !     Thus  related contrib eshort_reg_I is substracted.
        !(3)  The coulomb interaction of var cluster  ions with  unmovable ions i
        !     of Reg 2 is in etot_epe. The same interaction is in PG SCF energy.
        !     Thus  term ec_cluster_reg2a responsible for it on EPE side should
        !     be substracted.
        !(4)  cluster - reg1 coulomb is substracted from etot_epe when epe_lat_energy
        !     is calculated in main_epe_block

      etot_epe_corrected=(etot_epe-eshort_reg_I-ec_cluster_reg2a)/eau_ev

      ! to exclude double counting one also needs to substruct cluster_regI
      ! coulomb contrib see main_epe_block


        print*,'get_epe_energies etot_epe', etot_epe
        print*,'ec_cluster_reg2a/eau_ev', ec_cluster_reg2a/eau_ev
        print*,'eshort_reg_I/eau_ev',eshort_reg_I/eau_ev

      if(.true.) then
         !this contribution is calculated in pure epe treatment
         !therefore it should not be substructed to prevent double
         !counting
         etot_epe_corrected=etot_epe_corrected+ecoul_lat_relcorr/eau_ev ! a negative correction
         etot_epe_corrected=etot_epe_corrected-eshort_lat_relcorr/eau_ev
      endif

      if(extended_epe) &
           etot_epe_corrected=etot_epe_corrected-eshort_reg2a/eau_ev
   endif

        if(present(eshort_coupling_au)) then
        ! var contributions only as they will be recalculate in simol/optimizer
        eshort_coupling_au=eshort_reg_I/eau_ev
        if(extended_epe) eshort_coupling_au=eshort_coupling_au+eshort_reg2a/eau_ev
        endif

   if(present(lattice_energy)) then
      lattice_energy=0.0_r8_kind
      inquire(file=trim(data_dir)//"/epe_rel_converged",exist=epe_rel_converged)
      if(epe_rel_converged) then
         rel_converged_unit=&
              openget_iounit(file=trim(data_dir)//"/epe_rel_converged",&
              form='unformatted',status='unknown')
      else
         rel_converged_unit=&
              openget_iounit(file=trim(data_dir)//"/epe_rel_unconverged",&
              form='unformatted',status='unknown')
      endif
      read(rel_converged_unit) lattice_energy,epg_cluster_reg_I
      !                        etot_epe_corr,epg_cluster_reg_I
      ! saved in epe_field and forces
      print*,'with epe relaxation converged the lattice lattice interaction energy is' &
           ,lattice_energy,epg_cluster_reg_I

      call returnclose_iounit(rel_converged_unit)
      DPRINT 'returnclose_iounit for rel_converged_unit done'
   endif

   call stop_timer(timer_get_epe_energy)
 end subroutine get_epe_energies

 subroutine check_for_changes()
 real(kind=r8_kind)::lsd
   lsd=sum((epe(:reg_I_n_ions)%s(1)-epe_at_start(:reg_I_n_ions)%s(1))**2)
   lsd=lsd+sum((epe(:reg_I_n_ions)%s(2)-epe_at_start(:reg_I_n_ions)%s(2))**2)
   lsd=lsd+sum((epe(:reg_I_n_ions)%s(3)-epe_at_start(:reg_I_n_ions)%s(3))**2)
   lsd=lsd+sum((epe(:reg_I_n_ions)%c(1)-epe_at_start(:reg_I_n_ions)%c(1))**2)
   lsd=lsd+sum((epe(:reg_I_n_ions)%c(2)-epe_at_start(:reg_I_n_ions)%c(2))**2)
   lsd=lsd+sum((epe(:reg_I_n_ions)%c(3)-epe_at_start(:reg_I_n_ions)%c(3))**2)
!!$ if(lsd.ne.0.0_r8_kind) print*, '****LSD*******',lsd
 end subroutine check_for_changes

end module epecom_module




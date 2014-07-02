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
module efp_module
  !-------------------------------------------------------------------
  !
  !  Purpose: ...
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
#include "def.h"
  use type_module ! type specification parameters
  use efp_data_module
  use common_data_module
  use qmmm_interface_module, only: efp
  use datatype, only: arrmat2
  use pointcharge_module, only: rcm, no_fixed, qm_fixed, efp_fixed, present_X_centers
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  public efp, no_fixed, qm_fixed, efp_fixed
  integer(i4_kind), public :: n_efp
  logical, public :: print_id, do_pol
  real(r8_kind), public :: energy_conv
  integer(i4_kind), public :: n_density_updates

  type, public :: EFP_atoms
     character(len=3) :: name(3)
     real(r8_kind) :: coor(3,3)
  end type EFP_atoms
  type(EFP_atoms), public, allocatable :: efp_H2O_atoms(:)
  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: efp_read_input
  public :: efp_write_input
  public :: def_efp_arrays
  public :: read_gx_qm
  public :: calc_X_points
  public :: calc_efield_points
  public :: shutdown_efp_arrays
  public :: efp_sum_up_grads
  public :: efp_grad_cart_write
  public :: efp_write_gxfile
  public :: efp_max_grad
  public :: update_efps
  public :: efp_bcast
  public :: calc_efield_points1
#if 0
  public :: def_efp_arrays1
  public :: read_gx_qm1
#endif

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  integer(i4_kind) :: n_efp_uniq
  real(r8_kind), allocatable :: r_uniq(:,:)

  real(r8_kind), allocatable :: mp(:,:) !multipoles
  real(r8_kind), allocatable :: lmo(:,:) ! polarizable points
  real(r8_kind), allocatable :: rp_q(:,:) ! repulsive points
  real(r8_kind), allocatable :: rp_f(:,:) ! repulsive points

  real(r8_kind), allocatable :: Z(:) !charges
  real(r8_kind), allocatable :: C_z(:), A_z(:) !charge screening parameters
  real(r8_kind), allocatable :: C_zf(:), A_zf(:) !charge screening parameters
  real(r8_kind), allocatable :: D(:,:) !dipoles
  real(r8_kind), allocatable :: Q(:,:,:) !quadrupoles
  real(r8_kind), allocatable :: O(:,:,:,:) !octopoles
  real(r8_kind), allocatable :: pol_ten(:,:,:) !tensor of polarization
  real(r8_kind), allocatable :: C_q(:), A_q(:) !qm-efp repulsive parameters
  real(r8_kind), allocatable :: C_f(:,:), A_f(:,:) !efp-efp repulsive parameters

  type(arrmat2),allocatable :: grad_efp_cartesian(:)
  type(arrmat2),allocatable :: torque_efp_cartesian(:)

  real(r8_kind),allocatable :: grad_efp_cart_total(:,:),torque_efp_cart_total(:,:)

  character(len=12), allocatable :: mname(:), lname(:), rname_q(:), rname_f(:)

  integer(i4_kind), allocatable :: num(:,:)

  integer(i4_kind) :: gx_id, gx_line
  real(r8_kind) :: iwork
  logical :: ex_gx

  logical :: dftwater_efp 
  logical :: angs
  real(r8_kind) :: C1,C2,C3,C4,C5

  character(len=10) :: efp_solvent
  character(len=8) :: units !angstrom, au
  character(len=3) :: fixed_region
  logical :: do_micro_iter
  integer(i4_kind) :: update_ID_every
  logical :: calc_pc,calc_pd,calc_pq,calc_po,calc_pol,calc_rep !for test aims only
  logical :: calc_pol_pcm         !not for test aims only
  logical :: calc_scr             !for test aims only
  logical :: print_energy_contrib !for test aims only
  logical :: print_grad_contrib   !for test aims only
  logical :: print_grad_centers   !for test aims only
  logical :: print_efp_grad       !for test aims only
  logical :: print_efp_centers    !for test aims only
  logical :: print_id_cycles      !for test aims only
  logical :: print_as_in_input    !for test aims only

  character(len=5) :: name1, name2, name3
  real(r8_kind) :: coor1(3), coor2(3), coor3(3)

  character(len=10), parameter :: df_efp_solvent="DFT_WATER"
  integer(i4_kind), parameter :: df_n_efp=0
  character(len=8), parameter :: df_units="AU"
  character(len=3), parameter :: df_fixed_region ="NO" !"QM","EFP" 
  logical, parameter :: df_do_micro_iter=.false.
  integer(i4_kind) :: df_update_ID_every=1
  logical, parameter :: df_calc_pc=.true.
  logical, parameter :: df_calc_pd=.true.
  logical, parameter :: df_calc_pq=.true.
  logical, parameter :: df_calc_po=.true.
  logical, parameter :: df_calc_pol=.true.
  logical, parameter :: df_calc_pol_pcm=.true.
  logical, parameter :: df_calc_rep=.true.
  logical, parameter :: df_calc_scr=.true.
  logical, parameter :: df_print_energy_contrib = .false.
  logical, parameter :: df_print_grad_contrib = .false.
  logical, parameter :: df_print_grad_centers = .false.
  logical, parameter :: df_print_efp_grad = .false.
  logical, parameter :: df_print_efp_centers = .false.
  logical, parameter :: df_print_id_cycles = .false.  
  logical, parameter :: df_print_as_in_input = .true.
  real(r8_kind), parameter :: df_energy_conv = 1.0e-7_r8_kind
  integer(i4_kind), parameter :: df_n_density_updates = 15

  character(len=5), parameter :: df_name="XXXXX"
  real(r8_kind), parameter :: df_coor(3)=0.0_r8_kind


  namelist /efp_data/ n_efp, efp_solvent, units, fixed_region, energy_conv, n_density_updates, &
       do_micro_iter, update_ID_every, &
       calc_pc, calc_scr, calc_pd, calc_pq, &
       calc_po, calc_pol, calc_pol_pcm, calc_rep, &
       print_grad_contrib,print_efp_centers,print_id_cycles, &
       print_grad_centers,print_efp_grad,print_as_in_input,print_energy_contrib
  namelist/efp_water/ name1,coor1,name2,coor2,name3,coor3
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  !*************************************************************
  subroutine efp_read_input()
    !  Purpose: read namelist EFP_DATA
    !------------ Modules used ------------------- ---------------
    use input_module
    use iounitadmin_module, only : write_to_trace_unit,write_to_output_units
    use inp_out_module, only : upcase, check_string
    use pointcharge_module, only : print_energy
    use induced_dipoles_module, only: do_iterations, do_pol_pcm, n_update
    use operations_module, only: operations_solvation_effect
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i, j, unit, status
    !------------ Executable code ------------------------------------

    dftwater_efp=.false.
    angs=.false.

    n_efp=df_n_efp
    efp_solvent=trim(df_efp_solvent)
    fixed_region=trim(df_fixed_region)
    units=trim(df_units)
    do_micro_iter=df_do_micro_iter
    update_ID_every=df_update_ID_every
    calc_pc=df_calc_pc
    calc_pd=df_calc_pd
    calc_pq=df_calc_pq
    calc_po=df_calc_po
    calc_pol=df_calc_pol
    calc_pol_pcm=df_calc_pol_pcm
    calc_rep=df_calc_rep
    calc_scr=df_calc_scr
    print_energy_contrib=df_print_energy_contrib
    print_grad_contrib=df_print_grad_contrib
    print_grad_centers=df_print_grad_centers
    print_efp_grad=df_print_efp_grad
    print_efp_centers=df_print_efp_centers
    print_id_cycles=df_print_id_cycles
    print_as_in_input=df_print_as_in_input
    energy_conv=df_energy_conv
    n_density_updates=df_n_density_updates

    if(input_line_is_namelist("efp_data")) then
       unit = input_intermediate_unit()
       call input_read_to_intermediate
       read(unit, nml=efp_data, iostat=status)
       if (status .gt. 0) call input_error("efp_read_input: namelist efp_data")
    endif

    do_pol=calc_pol
    do_pol_pcm=calc_pol.and.calc_pol_pcm.and.(n_efp>0)
    do_iterations=do_micro_iter
    n_update=update_ID_every

    no_fixed=.false.; qm_fixed=.false.; efp_fixed=.false.
    call upcase(fixed_region)
    select case(trim(fixed_region))
    case("NO")
       no_fixed  =.true.
    case("QM")
       qm_fixed  =.true.
    case("EFP")
       efp_fixed =.true.
    case default
       no_fixed  =.true.
    end select
    if(operations_solvation_effect) no_fixed  =.true.

    if(efp_fixed) then
       print_grad_contrib=.false.
       print_grad_centers=.false.
       print_efp_grad    =.false.
    end if

    if(.not.calc_pc) calc_scr=.false.
    print_id=print_id_cycles
    print_energy=print_energy_contrib

    call upcase(efp_solvent)
    dftwater_efp=check_string(efp_solvent,"DFT_WATER")
    
    call upcase(units)
    angs=check_string(efp_solvent,"ANGSTROM")

    if(n_efp <= 0 .and. efp) then 
       n_efp=0
       call write_to_trace_unit("NO EFP DEFINED. EFP METHOD IS DISABLE !!!!!!!!!!")
       call write_to_output_units("NO EFP DEFINED. EFP METHOD IS DISABLE !!!!!!!!!!")
       print*,"NO EFP DEFINED. EFP METHOD IS DISABLE !!!!!!!!!!"
       efp=.false.
    end if

    if(.not.efp) return

    call init_efp_arrays()

    do i=1,n_efp
       if(dftwater_efp) then
          if ( .not. input_line_is_namelist("efp_water")) &
               call input_error("efp_read_input: namelist EFP_WATER expected")
          name1=df_name; name2=df_name; name3=df_name
          coor1=df_coor; coor1=df_coor; coor1=df_coor

          call input_read_to_intermediate()
          read(unit, nml=efp_water, iostat=status)
          if (status .gt. 0) call input_error("efp_read_input: namelist efp_water")

          call upcase(name1); call upcase(name2); call upcase(name3);

          if(.not.check_string(trim(name1),"O1") .and. .not.check_string(trim(name1),"H2")) &
               call input_error("efp_read_input: efp_water: wrong name for the first atom")
          if(.not.check_string(trim(name2),"O1") .and. .not.check_string(trim(name2),"H2")) &
               call input_error("efp_read_input: efp_water: wrong name for the second atom")
          if(.not.check_string(trim(name3),"H3")) &
               call input_error("efp_read_input: efp_water: wrong name for the third atom")
          if(trim(name1) == trim(name2)) &
               call input_error("efp_read_input: efp_water: the first and second atoms has to be different")

          num(1,i)=1; num(2,i)=2; num(3,i)=3

          j=1+(n_nuc+n_emp)*(i-1); mp(:,j)=coor1; mname(j)=trim(name1)
          if(check_string(mname(j),"H2")) num(1,i)=2
          efp_H2O_atoms(i)%name(1)=mname(j)
          efp_H2O_atoms(i)%coor(:,1)=mp(:,j)

          j=2+(n_nuc+n_emp)*(i-1); mp(:,j)=coor2; mname(j)=trim(name2)
          if(check_string(mname(j),"O1")) num(2,i)=1
          efp_H2O_atoms(i)%name(2)=mname(j)
          efp_H2O_atoms(i)%coor(:,2)=mp(:,j)

          j=3+(n_nuc+n_emp)*(i-1); mp(:,j)=coor3; mname(j)=trim(name3)
          efp_H2O_atoms(i)%name(3)=mname(j)
          efp_H2O_atoms(i)%coor(:,3)=mp(:,j)
       end if
    end do

  end subroutine efp_read_input
  !*************************************************************

  !*************************************************************
  subroutine efp_write_input(iounit)
    !  Purpose: write namelist EFP_DATA in input.out
    !------------ Modules used ------------------- ---------------
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i,j,k
    !------------ Executable code ------------------------------------

    call start("EFP_DATA","EFP_WRITE_INPUT",iounit,operations_echo_input_level)
    call intg("N_EFP            ",n_efp            ,df_n_efp      )
    call word("UNITS            ",units            ,df_units      )
    call word("EFP_SOLVENT      ",efp_solvent      ,df_efp_solvent)
    call word("FIXED_REGION     ",fixed_region     ,df_fixed_region)
    call real("ENERGY_CONV      ",energy_conv      ,df_energy_conv)
    call intg("N_DENSITY_UPDATES",n_density_updates,df_n_density_updates)
    call flag("DO_MICRO_ITER    ",do_micro_iter    ,df_do_micro_iter)
    call intg("UPDATE_ID_EVERY  ",update_ID_every  ,df_update_ID_every)
    call stop()

    if(.not.efp) return

    do i=1,n_efp
       if(dftwater_efp) then
          call start("EFP_WATER","EFP_WRITE_INPUT",iounit,operations_echo_input_level)
          k=1
          if(num(1,i)==2) k=2
          j=k+(n_nuc+n_emp)*(i-1)
          call   word("NAME1",mname(j),df_name  )
          call real_d("COOR1",mp(:,j) ,df_coor,3)

          k=2
          if(num(2,i)==1) k=1
          j=k+(n_nuc+n_emp)*(i-1)
          call   word("NAME2",mname(j),df_name  )
          call real_d("COOR2",mp(:,j) ,df_coor,3)

          k=3
          j=k+(n_nuc+n_emp)*(i-1)
          call   word("NAME3",mname(j),df_name  )
          call real_d("COOR3",mp(:,j) ,df_coor,3)

          call stop()
       end if
    end do

  end subroutine efp_write_input
  !*************************************************************

  !*****************************************************************************
  subroutine efp_bcast()
    !------------ Modules used ------------------- ---------------
    use comm,                only: comm_bcast                                  &
                                 , comm_rank
    use operations_module,   only: operations_solvation_effect
    !------------ Declaration of local variables ---------------------
    integer(i4_kind)   :: status
    !------------ Executable code ----------------------------------------------
    !
    call comm_bcast( efp       )
    call comm_bcast( n_efp     )
    call comm_bcast( no_fixed  )
    call comm_bcast( qm_fixed  )
    call comm_bcast( efp_fixed )
    call comm_bcast( do_pol    )
    !
    if(operations_solvation_effect) then
      if(allocated(rcm)) then
        if ( comm_rank() /= 0 ) then
          deallocate(rcm,stat=status)
          allocate( rcm(3,n_efp), stat=status )
          ASSERT(status==0)
        endif
        call comm_bcast( rcm )
      end if
    end if
    !
  end subroutine efp_bcast
  !*****************************************************************************

  !*************************************************************
  subroutine init_efp_arrays()
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n_mp,n_lmo,n_rp_q,n_rp_f
    integer(i4_kind) :: status
    !------------ Executable code ------------------------------------

    if(dftwater_efp) then
       n_mp=(n_nuc+n_emp)*n_efp
       n_lmo=n_pol*n_efp
       n_rp_q=(n_rep_q+n_rep_q)*n_efp
       n_rp_f=n_rep_f*n_efp
    end if

    allocate(efp_H2O_atoms(n_efp),stat=status)
    ASSERT(status==0)

    allocate(num(n_nuc,n_efp),stat=status)
    ASSERT(status==0)

    allocate(rcm(3,n_efp),stat=status)
    ASSERT(status==0)

    allocate(mp(3,n_mp),stat=status)
    ASSERT(status==0)
    mp=zero
    allocate(Z(n_mp),C_z(n_mp),A_z(n_mp),stat=status)
    ASSERT(status==0)
    Z=zero; C_z=zero; A_z=zero
    allocate(C_zf(n_mp),A_zf(n_mp),stat=status)
    ASSERT(status==0)
    C_zf=zero; A_zf=zero
    allocate(D(3,n_mp),stat=status)
    ASSERT(status==0)
    D=zero
    allocate(Q(3,3,n_mp),stat=status)
    ASSERT(status==0)
    Q=zero
    allocate(O(3,3,3,n_mp),stat=status)
    ASSERT(status==0)
    O=zero
    allocate(mname(n_mp),stat=status)
    ASSERT(status==0)

    allocate(lmo(3,n_lmo),stat=status)
    ASSERT(status==0)
    allocate(pol_ten(3,3,n_lmo),stat=status)
    ASSERT(status==0)
    pol_ten=zero
    allocate(lname(n_lmo),stat=status)
    ASSERT(status==0)

    allocate(rp_q(3,n_rp_q),C_q(n_rp_q),A_q(n_rp_q),stat=status)
    ASSERT(status==0)
    C_q=zero; A_q=zero
    allocate(rp_f(3,n_rp_f),C_f(n_rep_f,n_rp_f),A_f(n_rep_f,n_rp_f),stat=status)
    ASSERT(status==0)
    C_f=zero; A_f=zero
    allocate(rname_q(n_rp_q),rname_f(n_rp_f),stat=status)
    ASSERT(status==0)

  end subroutine init_efp_arrays
  !*************************************************************

  !*************************************************************
  subroutine shutdown_efp_arrays()
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: status
    !------------ Executable code ------------------------------------

    deallocate(num,stat=status)
    ASSERT(status==0)

    if(qm_fixed) then
       deallocate(rcm,stat=status)
       ASSERT(status==0)
    end if

    deallocate(efp_H2O_atoms,stat=status)
    ASSERT(status==0)
    deallocate(mp,stat=status)
    ASSERT(status==0)
    deallocate(Z,C_z,A_z,stat=status)
    ASSERT(status==0)
    deallocate(C_zf,A_zf,stat=status)
    ASSERT(status==0)
    deallocate(D,stat=status)
    ASSERT(status==0)
    deallocate(Q,stat=status)
    ASSERT(status==0)
    deallocate(O,stat=status)
    ASSERT(status==0)
    deallocate(mname,stat=status)
    ASSERT(status==0)

    deallocate(lmo,stat=status)
    ASSERT(status==0)
    deallocate(pol_ten,stat=status)
    ASSERT(status==0)
    deallocate(lname,stat=status)
    ASSERT(status==0)

    deallocate(rp_q,C_q,A_q,stat=status)
    ASSERT(status==0)
    deallocate(rp_f,C_f,A_f,stat=status)
    ASSERT(status==0)
    deallocate(rname_q,rname_f,stat=status)
    ASSERT(status==0)

  end subroutine shutdown_efp_arrays
  !*************************************************************

  !*************************************************************
  subroutine correct_efp_water
    ! Purpose: correct (if necessary) structure of read in efp waters
    !          to correspond tabulated parameters
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: OH,HOH
    real(r8_kind) :: r(3,3),r12(3),r13(3),d12,d13,cos_a,angle
    real(r8_kind) :: rm(3,3),zmat(3),rq(3),r1q(3),d1q
    integer(i4_kind) :: i,j,k
    !------------ Executable code ------------------------------------

    OH=get_dft_OH(); HOH=get_dft_HOH()

!print*,OH*b2a,HOH*rad2deg

    do i=1,n_efp
       j=1+(n_nuc+n_emp)*(i-1)
       k=num(1,i)
       r(:,k)=mp(:,j)

       j=2+(n_nuc+n_emp)*(i-1)
       k=num(2,i)
       r(:,k)=mp(:,j)

       j=3+(n_nuc+n_emp)*(i-1); r(:,3)=mp(:,j)

       r12=r(:,2)-r(:,1); d12=sqrt(dot_product(r12,r12))
       r13=r(:,3)-r(:,1); d13=sqrt(dot_product(r13,r13))

       !bond correction
       if(d12 /= OH) then
          r12=r12*OH/d12; r(:,2)=r(:,1)+r12
       end if
       if(d13 /= OH) then
          r13=r13*OH/d13; r(:,3)=r(:,1)+r13
       end if

       !angle correction
       cos_a=dot_product(r12,r13)/(d12*d13); angle=acos(cos_a)
       if(angle /= HOH) then
          rq=(r(:,2)+r(:,3))/two
          r1q=rq-r(:,1); d1q=sqrt(dot_product(r1q,r1q))

          rm(:,1)=r(:,3); rm(:,2)=rq; rm(:,3)=r(:,1)
          zmat(1)=OH; zmat(2)=HOH/two; zmat(3)=pi
          call add_1_4_atom(r(:,2),rm,zmat)

          rm(:,1)=r(:,2); rm(:,2)=rq; rm(:,3)=r(:,1)
          zmat(1)=OH; zmat(2)=HOH/two; zmat(3)=pi
          call add_1_4_atom(r(:,3),rm,zmat)
       end if

       j=1+(n_nuc+n_emp)*(i-1); mp(:,j)=r(:,1)
       j=2+(n_nuc+n_emp)*(i-1); mp(:,j)=r(:,2)
       j=3+(n_nuc+n_emp)*(i-1); mp(:,j)=r(:,3)

       r12=r(:,2)-r(:,1); d12=sqrt(dot_product(r12,r12))
       r13=r(:,3)-r(:,1); d13=sqrt(dot_product(r13,r13))
       cos_a=dot_product(r12,r13)/(d12*d13); angle=acos(cos_a)
!print*,d12*b2a,d13*b2a,angle*rad2deg
    end do

  end subroutine correct_efp_water
  !*************************************************************

  !*************************************************************
  function get_dft_OH() result(length)
    ! Purpose : calculate length of OH bond
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: length
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: r1(3),r2(3),r3(3)
    !------------ Executable code ------------------------------------

    r1=dft_water_nuc(1)%coor
    r2=dft_water_nuc(2)%coor
    r3=dft_water_nuc(3)%coor

    length=sqrt(dot_product(r2-r1,r2-r1))
    
  end function get_dft_OH
  !*************************************************************

  !*************************************************************
  function get_dft_HOH() result(angle)
    ! Purpose : calculate HOH angle
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: angle
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: r1(3),r2(3),r3(3),d12,d13,cos_a
    !------------ Executable code ------------------------------------

    r1=dft_water_nuc(1)%coor
    r2=dft_water_nuc(2)%coor
    r3=dft_water_nuc(3)%coor

    d12=sqrt(dot_product(r2-r1,r2-r1))
    d13=sqrt(dot_product(r3-r1,r3-r1))

    cos_a=dot_product(r2-r1,r3-r1)/(d12*d13)
    angle=acos(cos_a)
    
  end function get_dft_HOH
  !*************************************************************

  !*************************************************************
  subroutine struct_H2O_param()
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: d16,d171,d172,d191,d192,d1n,dnm,d1l
    real(r8_kind) :: r1(3),r2(3),r3(3)
    real(r8_kind) :: p1(3),p2(3),p3(3),p4(3),p5(3)
    real(r8_kind) :: r1n(3),rnm(3),r1l(3)
    !------------ Executable code ------------------------------------

    r1=dft_water_nuc(1)%coor
    r2=dft_water_nuc(2)%coor
    r3=dft_water_nuc(3)%coor

    r1n=(r2+r3)/two-r1; d1n=sqrt(dot_product(r1n,r1n))
    rnm=(r3-r2)/two;    dnm=sqrt(dot_product(rnm,rnm))
    r1l=vector_product(r2-r1,r3-r1); d1l=sqrt(dot_product(r1l,r1l))

    p1=dft_water_pol(1)%coor
    p2=dft_water_pol(2)%coor
    p3=dft_water_pol(3)%coor
    p4=dft_water_pol(4)%coor
    p5=dft_water_pol(5)%coor

    d16 =abs(r1(3)-p1(3))
    d171=abs(r1(3)-p2(3))
    d172=abs(r1(1)-p2(1))
    d191=abs(r1(3)-p4(3))
    d192=abs(r1(2)-p4(2))
    
    C1=d16/d1n
    C2=d171/d1n
    C3=d172/dnm
    C4=d191/d1n
    C5=d192/d1l

  end subroutine struct_H2O_param
  !*************************************************************

  !*************************************************************
  subroutine add_1_4_atom(ra,rb,internals)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: ra(3),rb(3,3)
    real(r8_kind) :: internals(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: bond,angle,dihedral
    real(r8_kind) :: cos_a,sin_a,cos_d,sin_d
    real(r8_kind) :: r21(3),d21,r32(3),d32,r1(3),d1,r2(3),d2,r3(3),r4(3),d4
    !------------ Executable code ------------------------------------

    bond=internals(1)
    angle=internals(2)
    dihedral=internals(3)

    cos_a=cos(angle); sin_a=sin(angle)
    cos_d=cos(dihedral); sin_d=sin(dihedral)

    r21=rb(:,2)-rb(:,1); r32=rb(:,3)-rb(:,2)
    d21=sqrt(dot_product(r21,r21))
    r21=r21/d21
    d32=sqrt(dot_product(r32,r32))
    r32=r32/d32

    r1=vector_product(r21,r32)
    d1=sqrt(dot_product(r1,r1))
    r1=r1/d1

    r2=vector_product(r1,r32)
    d2=sqrt(dot_product(r2,r2))
    r2=r2/d2

    r3=r2*cos_d+r1*sin_d

    r4=-r32*cos_a+r3*sin_a
    d4=sqrt(dot_product(r4,r4))
    r4=r4*bond/d4

    ra=rb(:,3)+r4

  end subroutine add_1_4_atom
  !*************************************************************

  !*************************************************************
  function vector_product(v1,v2)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: vector_product(3)
    real(r8_kind) :: v1(3),v2(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !*************************************************************

  !*************************************************************
  subroutine def_efp_arrays()
    !------------ Modules used -----------------------------------
    use operations_module, only: operations_geo_opt, operations_read_gx
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,l,m,status,n_points
    real(r8_kind) :: r(3,3),rmat1(3,3),rmat2(3,3),rmat3(3,3),rmat_3(3,3),rbuf(3)
    real(r8_kind), allocatable :: rpf(:,:), rcf(:,:)
    real(r8_kind) :: dip(3,n_emp),quad(3,3,n_emp),oct(3,3,3,n_emp),pol(3,3,n_pol)
    real(r8_kind) :: qa(3,3),oc(3,3,3),summ(3)
    !------------ Executable code ------------------------------------

    water: if(dftwater_efp) then
       call struct_H2O_param()

       if(operations_geo_opt .or. operations_read_gx) call read_gx_efp_water()
       call correct_efp_water()

       n_points=n_nuc+n_emp+n_pol+n_rep_q+n_rep_q+n_rep_f
       ASSERT(23==n_points)
       allocate(rpf(3,n_nuc),rcf(3,n_points),stat=status)
       ASSERT(status==0)

       !define coordinate system of predefined water fragment
       do i=1,3
          r(:,i)=dft_water_nuc(i)%coor
       end do

       call rotmat(r,rmat1)
       rmat1=transpose(rmat1)

       j=0
       do i=1,n_nuc
          j=j+1
          rpf(:,j)=dft_water_nuc(i)%coor
       end do

       rbuf=rpf(:,1)
       do i=1,n_nuc
          rpf(:,i)=rpf(:,i)-rbuf(:)
       end do

       n_eff_frag: do i=1,n_efp

          !define coordinate system of the i-th fragment
          do j=1,3
             k=j+(n_nuc+n_emp)*(i-1)
             l=j
             r(:,l)=mp(:,k)
          end do
          call rotmat(r,rmat2)

          !find rotation matrix connecting two coordinate system
          rmat3=matmul(rmat2,rmat1)
          rmat_3=transpose(rmat3)

          !rotate and translate predefined fragment to coinside with i-th fragment
          rcf=matmul(rmat3,rpf)
          do j=1,n_nuc
             rcf(:,j)=rcf(:,j)+r(:,1)
          end do
!!$if(i==1 ) rcf(2,1)=rcf(2,1)-0.0001_r8_kind
          ! multipoles
          rcf(:,4) =rcf(:,1) !1
          rcf(:,5) =rcf(:,2) !2
          rcf(:,6) =rcf(:,3) !3
          rcf(:,7) =(rcf(:,1)+rcf(:,2))/two !4
          rcf(:,8) =(rcf(:,1)+rcf(:,3))/two !5
          ! polarizable points
          rcf(:,9) =(one-C1)*rcf(:,1)+C1*(rcf(:,2)+rcf(:,3))/two !6
          rcf(:,10)=(one-C2)*rcf(:,1)+C2*(rcf(:,2)+rcf(:,3))/two+C3*(rcf(:,3)-rcf(:,2))/two !7
          rcf(:,11)=(one-C2)*rcf(:,1)+C2*(rcf(:,2)+rcf(:,3))/two-C3*(rcf(:,3)-rcf(:,2))/two !8
          rcf(:,12)=(one+C4)*rcf(:,1)-C4*(rcf(:,2)+rcf(:,3))/two- &
               C5*vector_product(rcf(:,2)-rcf(:,1),rcf(:,3)-rcf(:,1)) !9
          rcf(:,13)=(one+C4)*rcf(:,1)-C4*(rcf(:,2)+rcf(:,3))/two+ &
               C5*vector_product(rcf(:,2)-rcf(:,1),rcf(:,3)-rcf(:,1)) !10
!!$print*,'6',rcf(:,9)
!!$print*,'7',rcf(:,10)
!!$print*,'8',rcf(:,11)
!!$print*,'9',rcf(:,12)
!!$print*,'10',rcf(:,13)
          ! qm - efp repulsive points
          rcf(:,14)=rcf(:,1) !1
          rcf(:,15)=rcf(:,2) !2
          rcf(:,16)=rcf(:,3) !3
          rcf(:,17)=rcf(:,1) !4
          rcf(:,18)=rcf(:,2) !5
          rcf(:,19)=rcf(:,3) !6
          ! efp - efp repulsive points
          rcf(:,20)=rcf(:,1) !1
          rcf(:,21)=rcf(:,2) !2
          rcf(:,22)=rcf(:,3) !3
          rcf(:,23)=(rcf(:,1)*dft_water_nuc(1)%mass+ &
                     rcf(:,2)*dft_water_nuc(2)%mass+ &
                     rcf(:,3)*dft_water_nuc(3)%mass)/ &
                     (dft_water_nuc(1)%mass+dft_water_nuc(2)%mass+dft_water_nuc(3)%mass) !11

          rcm(:,i)=rcf(:,23)
!          rcm(:,i)=(rcf(:,1)+rcf(:,2)+rcf(:,3))/three

          !rotating dipoles
          do j=1,n_emp
             dip(:,j)=dft_water_mp(j)%dip
          end do
          dip=matmul(rmat3,dip)

          !rotating quadrupoles and transformation to traceless form
          do j=1,n_emp
             qa(:,:)=dft_water_mp(j)%quad
             qa(:,:)=matmul(matmul(rmat3,qa(:,:)),rmat_3)

             summ(1)=zero
             do k=1,3
                summ(1)=summ(1)+qa(k,k)
             end do

             do k=1,3
                do l=1,3
                   quad(k,l,j)=three*qa(k,l)
                   if(k==l) quad(k,l,j)=quad(k,l,j)-summ(1)
                end do
             end do
             quad(:,:,j)=quad(:,:,j)*half
          end do

          !rotating octopoles and transformation to traceless form
          do j=1,n_emp
             oc(:,:,:)=dft_water_mp(j)%oct
             do k=1,3
                oc(:,:,k)=matmul(matmul(rmat3,oc(:,:,k)),rmat_3)
             end do
             do k=1,3
                oc(:,k,:)=matmul(oc(:,k,:),rmat_3)
             end do

             summ=zero
             do k=1,3
                summ(1)=summ(1)+oc(1,k,k)
             end do
             do k=1,3
                summ(2)=summ(2)+oc(2,k,k)
             end do
             do k=1,3
                summ(3)=summ(3)+oc(3,k,k)
             end do

             do k=1,3
                do l=1,3
                   do m=1,3
                      oct(k,l,m,j)=five*oc(k,l,m)
                      if(l==m) oct(k,l,m,j)=oct(k,l,m,j)-summ(k)
                      if(k==m) oct(k,l,m,j)=oct(k,l,m,j)-summ(l)
                      if(k==l) oct(k,l,m,j)=oct(k,l,m,j)-summ(m)
                   end do
                end do
             end do
             oct(:,:,:,j)=oct(:,:,:,j)*half
          end do

          !rotating tensors of polarization
          do j=1,n_pol
             pol(:,:,j)=dft_water_pol(j)%alpha
             pol(:,:,j)=matmul(matmul(rmat3,pol(:,:,j)),rmat_3)
          end do

          call fill_efp_arrays(rcf,dip,quad,oct,pol,i)

       end do n_eff_frag

       deallocate(rpf,rcf,stat=status)
       ASSERT(status==0)
    end if water

    call print_efp_coor()

  end subroutine def_efp_arrays
  !*************************************************************

  !*************************************************************
#if 0
  subroutine def_efp_arrays1()
    !------------ Modules used -----------------------------------
    use operations_module, only: operations_geo_opt, operations_read_gx
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,l,m,status,n_points
    real(r8_kind) :: r(3,3),rmat1(3,3),rmat2(3,3),rmat3(3,3),rmat_3(3,3),rbuf(3)
    real(r8_kind), allocatable :: rpf(:,:), rcf(:,:)
    real(r8_kind) :: dip(3,n_emp),quad(3,3,n_emp),oct(3,3,3,n_emp),pol(3,3,n_pol)
    real(r8_kind) :: qa(3,3),oc(3,3,3),summ(3)
    !------------ Executable code ------------------------------------

    water: if(dftwater_efp) then
       call struct_H2O_param()

       if(operations_geo_opt .or. operations_read_gx) call read_gx_efp_water1()
       call correct_efp_water()

       n_points=n_nuc+n_emp+n_pol+n_rep_q+n_rep_q+n_rep_f
       ASSERT(23==n_points)
       allocate(rpf(3,n_nuc),rcf(3,n_points),stat=status)
       ASSERT(status==0)

       !define coordinate system of predefined water fragment
       do i=1,3
          r(:,i)=dft_water_nuc(i)%coor
       end do

       call rotmat(r,rmat1)
       rmat1=transpose(rmat1)

       j=0
       do i=1,n_nuc
          j=j+1
          rpf(:,j)=dft_water_nuc(i)%coor
       end do

       rbuf=rpf(:,1)
       do i=1,n_nuc
          rpf(:,i)=rpf(:,i)-rbuf(:)
       end do

       n_eff_frag: do i=1,n_efp

          !define coordinate system of the i-th fragment
          do j=1,3
             k=j+(n_nuc+n_emp)*(i-1)
             l=j
             r(:,l)=mp(:,k)
          end do
          call rotmat(r,rmat2)

          !find rotation matrix connecting two coordinate system
          rmat3=matmul(rmat2,rmat1)
          rmat_3=transpose(rmat3)

          !rotate and translate predefined fragment to coinside with i-th fragment
          rcf=matmul(rmat3,rpf)
          do j=1,n_nuc
             rcf(:,j)=rcf(:,j)+r(:,1)
          end do
!!$if(i==1 ) rcf(2,1)=rcf(2,1)-0.0001_r8_kind
          ! multipoles
          rcf(:,4) =rcf(:,1) !1
          rcf(:,5) =rcf(:,2) !2
          rcf(:,6) =rcf(:,3) !3
          rcf(:,7) =(rcf(:,1)+rcf(:,2))/two !4
          rcf(:,8) =(rcf(:,1)+rcf(:,3))/two !5
          ! polarizable points
          rcf(:,9) =(one-C1)*rcf(:,1)+C1*(rcf(:,2)+rcf(:,3))/two !6
          rcf(:,10)=(one-C2)*rcf(:,1)+C2*(rcf(:,2)+rcf(:,3))/two+C3*(rcf(:,3)-rcf(:,2))/two !7
          rcf(:,11)=(one-C2)*rcf(:,1)+C2*(rcf(:,2)+rcf(:,3))/two-C3*(rcf(:,3)-rcf(:,2))/two !8
          rcf(:,12)=(one+C4)*rcf(:,1)-C4*(rcf(:,2)+rcf(:,3))/two- &
               C5*vector_product(rcf(:,2)-rcf(:,1),rcf(:,3)-rcf(:,1)) !9
          rcf(:,13)=(one+C4)*rcf(:,1)-C4*(rcf(:,2)+rcf(:,3))/two+ &
               C5*vector_product(rcf(:,2)-rcf(:,1),rcf(:,3)-rcf(:,1)) !10
!!$print*,'6',rcf(:,9)
!!$print*,'7',rcf(:,10)
!!$print*,'8',rcf(:,11)
!!$print*,'9',rcf(:,12)
!!$print*,'10',rcf(:,13)
          ! qm - efp repulsive points
          rcf(:,14)=rcf(:,1) !1
          rcf(:,15)=rcf(:,2) !2
          rcf(:,16)=rcf(:,3) !3
          rcf(:,17)=rcf(:,1) !4
          rcf(:,18)=rcf(:,2) !5
          rcf(:,19)=rcf(:,3) !6
          ! efp - efp repulsive points
          rcf(:,20)=rcf(:,1) !1
          rcf(:,21)=rcf(:,2) !2
          rcf(:,22)=rcf(:,3) !3
          rcf(:,23)=(rcf(:,1)*dft_water_nuc(1)%mass+ &
                     rcf(:,2)*dft_water_nuc(2)%mass+ &
                     rcf(:,3)*dft_water_nuc(3)%mass)/ &
                     (dft_water_nuc(1)%mass+dft_water_nuc(2)%mass+dft_water_nuc(3)%mass) !11

          rcm(:,i)=rcf(:,23)
!          rcm(:,i)=(rcf(:,1)+rcf(:,2)+rcf(:,3))/three

          !rotating dipoles
          do j=1,n_emp
             dip(:,j)=dft_water_mp(j)%dip
          end do
          dip=matmul(rmat3,dip)

          !rotating quadrupoles and transformation to traceless form
          do j=1,n_emp
             qa(:,:)=dft_water_mp(j)%quad
             qa(:,:)=matmul(matmul(rmat3,qa(:,:)),rmat_3)

             summ(1)=zero
             do k=1,3
                summ(1)=summ(1)+qa(k,k)
             end do

             do k=1,3
                do l=1,3
                   quad(k,l,j)=three*qa(k,l)
                   if(k==l) quad(k,l,j)=quad(k,l,j)-summ(1)
                end do
             end do
             quad(:,:,j)=quad(:,:,j)*half
          end do

          !rotating octopoles and transformation to traceless form
          do j=1,n_emp
             oc(:,:,:)=dft_water_mp(j)%oct
             do k=1,3
                oc(:,:,k)=matmul(matmul(rmat3,oc(:,:,k)),rmat_3)
             end do
             do k=1,3
                oc(:,k,:)=matmul(oc(:,k,:),rmat_3)
             end do

             summ=zero
             do k=1,3
                summ(1)=summ(1)+oc(1,k,k)
             end do
             do k=1,3
                summ(2)=summ(2)+oc(2,k,k)
             end do
             do k=1,3
                summ(3)=summ(3)+oc(3,k,k)
             end do

             do k=1,3
                do l=1,3
                   do m=1,3
                      oct(k,l,m,j)=five*oc(k,l,m)
                      if(l==m) oct(k,l,m,j)=oct(k,l,m,j)-summ(k)
                      if(k==m) oct(k,l,m,j)=oct(k,l,m,j)-summ(l)
                      if(k==l) oct(k,l,m,j)=oct(k,l,m,j)-summ(m)
                   end do
                end do
             end do
             oct(:,:,:,j)=oct(:,:,:,j)*half
          end do

          !rotating tensors of polarization
          do j=1,n_pol
             pol(:,:,j)=dft_water_pol(j)%alpha
             pol(:,:,j)=matmul(matmul(rmat3,pol(:,:,j)),rmat_3)
          end do

          call fill_efp_arrays(rcf,dip,quad,oct,pol,i)

       end do n_eff_frag

       deallocate(rpf,rcf,stat=status)
       ASSERT(status==0)
    end if water

    call print_efp_coor()

  end subroutine def_efp_arrays1
#endif
  !*************************************************************

  !*************************************************************
  subroutine print_efp_coor()
    ! Purpose : print coordinates of the main EFP centers
    !------------ Modules used -----------------------------------
    use iounitadmin_module, only: output_unit
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,n,k
    !------------ Executable code ------------------------------------

    write(output_unit,'(a50)') '=========== EFFECTIVE FRAGMENTS (AU) ============'
    write(output_unit,'(a1)') ''
    do i=1,n_efp
       write(output_unit,'(a9,3x,i4)') 'Fragment ',i
       write(output_unit,'(a50)') '--------------------------------------------------'
       n=(i-1)*(n_nuc+n_emp)
       do j=1,n_gx_points
          k=j
          if(print_as_in_input) k=num(j,i)
          write(output_unit,'(a4,3x,3f15.10)') mname(k+n), mp(:,k+n)
       end do
       write(output_unit,'(a1)') ''
    end do
    write(output_unit,'(a1)') ''
    write(output_unit,'(a50)') '========== EFFECTIVE FRAGMENTS (ANG) ============'
    write(output_unit,'(a1)') ''
    do i=1,n_efp
       write(output_unit,'(a9,3x,i4)') 'Fragment ',i
       write(output_unit,'(a50)') '--------------------------------------------------'
       n=(i-1)*(n_nuc+n_emp)
       do j=1,n_gx_points
          k=j
          if(print_as_in_input) k=num(j,i)
          write(output_unit,'(a4,3x,3f15.10)') mname(k+n), mp(:,k+n)*b2a
       end do
       write(output_unit,'(a1)') ''
    end do
    write(output_unit,'(a1)') ''

  end subroutine print_efp_coor
  !*************************************************************

  !*************************************************************
  subroutine read_gx_qm()
    ! Purpose : read coordinates of QM atoms from GX file
    !------------ Modules used -----------------------------------
    use iounitadmin_module, only: openget_iounit
    use filename_module, only: inpfile
    use unique_atom_module, only: N_unique_atoms, unique_atoms
    use operations_module, only: operations_geo_opt, operations_read_gx
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iuniq
    real(r8_kind) :: at_num,coor(3)
    character(len=5) ::  gx_buff
    integer(i4_kind) :: i,counter_equal
    !------------ Executable code ------------------------------------

    inquire(file=inpfile('gxfile'), exist=ex_gx)
    if((operations_geo_opt .or. operations_read_gx) .and. ex_gx) then
       gx_id = openget_iounit(file=inpfile('gxfile'), status='old', &
            form='formatted')

       gx_line=0
       do i=1,N_unique_atoms
          counter_equal=1
          do 
             if(counter_equal>unique_atoms(i)%n_equal_atoms) exit
             gx_line=gx_line+1
             read(gx_id,*,end=100,err=101) at_num,coor,iuniq
             if(iuniq /= 0) then
                if(counter_equal==1) then
                   unique_atoms(i)%position_first_ea(1)=coor(1)
                   unique_atoms(i)%position_first_ea(2)=coor(2)
                   unique_atoms(i)%position_first_ea(3)=coor(3)
                endif
                counter_equal=counter_equal+1
             end if
          end do
       end do
    else
       gx_line=0
       do i=1,N_unique_atoms
          gx_line=gx_line+unique_atoms(i)%n_equal_atoms
       end do
    end if

    return

100 call error_handler("efp_module: read_gx_qm:  attempt to read after end of GX file")
101 write(gx_buff,'(i5)') gx_line
    call error_handler("efp_module: read_gx_qm: wrong GX file. Line - "//trim(gx_buff))

  end subroutine read_gx_qm
  !*************************************************************

  !*************************************************************
#if 0
  subroutine read_gx_qm1()
    ! Purpose : read coordinates of QM atoms from GX file
    !------------ Modules used -----------------------------------
    use iounitadmin_module, only: openget_iounit
    use filename_module, only: inpfile
    use unique_atom_module, only: N_unique_atoms, unique_atoms
    use operations_module, only: operations_geo_opt, operations_read_gx
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iuniq
    real(r8_kind) :: at_num,coor(3)
    character(len=5) ::  gx_buff
    integer(i4_kind) :: i,counter_equal
    !------------ Executable code ------------------------------------

    inquire(file=inpfile('gxfile1'), exist=ex_gx)
    if((operations_geo_opt .or. operations_read_gx) .and. ex_gx) then
       gx_id = openget_iounit(file=inpfile('gxfile1'), status='old', &
            form='formatted')

       gx_line=0
       do i=1,N_unique_atoms
          counter_equal=1
          do 
             if(counter_equal>unique_atoms(i)%n_equal_atoms) exit
             gx_line=gx_line+1
             read(gx_id,*,end=100,err=101) at_num,coor,iuniq
             if(iuniq /= 0) then
                if(counter_equal==1) then
                   unique_atoms(i)%position_first_ea(1)=coor(1)
                   unique_atoms(i)%position_first_ea(2)=coor(2)
                   unique_atoms(i)%position_first_ea(3)=coor(3)
                endif
                counter_equal=counter_equal+1
             end if
          end do
       end do
    else
       gx_line=0
       do i=1,N_unique_atoms
          gx_line=gx_line+unique_atoms(i)%n_equal_atoms
       end do
    end if

    return

100 call error_handler("efp_module: read_gx_qm:  attempt to read after end of GX file")
101 write(gx_buff,'(i5)') gx_line
    call error_handler("efp_module: read_gx_qm: wrong GX file. Line - "//trim(gx_buff))

  end subroutine read_gx_qm1
#endif
  !*************************************************************

  !*************************************************************
  subroutine read_gx_efp_water()
    ! Purpose : read coordinates of fragments from GX file
    !------------ Modules used ------------------- ---------------
    use iounitadmin_module, only: returnclose_iounit
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iuniq
    real(r8_kind) :: at_num,coor(3)
    character(len=5) ::  gx_buff
    integer(i4_kind) :: i,j,k,l,m
    !------------ Executable code ------------------------------------

    if(ex_gx) then
    do i=1,n_efp
       k=0
       do j=1,n_gx_points
1         gx_line=gx_line+1
          read(gx_id,*,end=100,err=101) at_num,coor,iuniq
          if(int(at_num) == 99 .or. iuniq == 0) goto 1
          k=k+1
          m=num(k,i)
          if(abs(at_num-gx_num(m)) .ge. 1.0e-2_r8_kind) goto 102
          l=k+(n_nuc+n_emp)*(i-1); mp(:,l)=coor
       end do
    end do

    gx_line=gx_line+1
    read(gx_id,*,end=100,err=101) iwork

!!$    call returnclose_iounit(gx_id)

    else
       gx_line=gx_line+n_efp*3+1
       iwork=-1.0_r8_kind
    endif

    return

100 call error_handler("efp_module: read_gx_efp_water:  attempt to read after end of GX file")
101 write(gx_buff,'(i5)') gx_line
    call error_handler("efp_module: read_gx_efp_water: wrong GX file. Line - "//trim(gx_buff))
102 write(gx_buff,'(i5)') gx_line
    call error_handler( &
         "efp_module: read_gx_efp_water: wrong point order within water fragment. Line - "//trim(gx_buff))

  end subroutine read_gx_efp_water
  !*************************************************************

  !*************************************************************
#if 0
  subroutine read_gx_efp_water1()
    ! Purpose : read coordinates of fragments from GX file
    !------------ Modules used ------------------- ---------------
    use iounitadmin_module, only: returnclose_iounit
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iuniq
    real(r8_kind) :: at_num,coor(3)
    character(len=5) ::  gx_buff
    integer(i4_kind) :: i,j,k,l,m
    !------------ Executable code ------------------------------------

    if(ex_gx) then
    do i=1,n_efp
       k=0
       do j=1,n_gx_points
1         gx_line=gx_line+1
          read(gx_id,*,end=100,err=101) at_num,coor,iuniq
!if(i==1) coor(1)=coor(1)-0.0001_r8_kind
          if(int(at_num) == 99 .or. iuniq == 0) goto 1
          k=k+1
          m=num(k,i)
          if(abs(at_num-gx_num(m)) .ge. 1.0e-2_r8_kind) goto 102
          l=k+(n_nuc+n_emp)*(i-1); mp(:,l)=coor
       end do
    end do

    gx_line=gx_line+1
    read(gx_id,*,end=100,err=101) iwork

!!$    call returnclose_iounit(gx_id)

    else
       gx_line=gx_line+n_efp*3+1
       iwork=-1.0_r8_kind
    endif

    return

100 call error_handler("efp_module: read_gx_efp_water:  attempt to read after end of GX file")
101 write(gx_buff,'(i5)') gx_line
    call error_handler("efp_module: read_gx_efp_water: wrong GX file. Line - "//trim(gx_buff))
102 write(gx_buff,'(i5)') gx_line
    call error_handler( &
         "efp_module: read_gx_efp_water: wrong point order within water fragment. Line - "//trim(gx_buff))

  end subroutine read_gx_efp_water1
#endif
  !*************************************************************

  !*************************************************************
  subroutine rotmat(r,rotm)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(in)  :: r(3,3)
    real(r8_kind), intent(out) :: rotm(3,3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: r12(3),r13(3),d12,d13
    real(r8_kind) :: e1(3),e2(3),e3(3)
    !------------ Executable code ------------------------------------

    r12=r(:,2)-r(:,1); r13=r(:,3)-r(:,1)
    d12=sqrt(dot_product(r12,r12)); d13=sqrt(dot_product(r13,r13))

    e1=r12/d12; e2=r13/d13; e2=e2-dot_product(e1,e2)*e1
    e3=vector_product(e1,e2)

    e1=e1/sqrt(dot_product(e1,e1))
    e2=e2/sqrt(dot_product(e2,e2))
    e3=e3/sqrt(dot_product(e3,e3))

    rotm(:,1)=e1; rotm(:,2)=e2; rotm(:,3)=e3

  end subroutine rotmat
  !*************************************************************

  !*************************************************************
  subroutine fill_efp_arrays(coor,dp,qp,op,pt,i_efp)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(in) :: coor(:,:),dp(:,:),qp(:,:,:),op(:,:,:,:),pt(:,:,:)
    integer(i4_kind), intent(in) :: i_efp
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,n
    real(r8_kind) :: rmass(3),frag_mass
    !------------ Executable code ------------------------------------

!!$print*,'FRAGMENT :',i_efp
    j=0
    n=(i_efp-1)*(n_nuc+n_emp)
    rmass=zero; frag_mass=zero
    do i=1,n_nuc+n_emp
       k=i
       j=j+1
       mp(:,n+k)=coor(:,j)
       if(i <= n_nuc) then
          mname(n+i) = dft_water_nuc(k)%name
          Z(n+i) = dft_water_nuc(k)%z%q
          D(:,n+i) = dft_water_nuc(k)%dip
          Q(:,:,n+i) = dft_water_nuc(k)%quad
          O(:,:,:,n+i) = dft_water_nuc(k)%oct
          rmass=rmass+mp(:,n+i)*dft_water_nuc(k)%mass
          frag_mass=frag_mass+dft_water_nuc(k)%mass
       else
          mname(n+i) = dft_water_mp(k-n_nuc)%name
          Z(n+i) = dft_water_mp(k-n_nuc)%z%q
          if(calc_scr) then
             C_z(n+i) = dft_water_mp(k-n_nuc)%z%C
             A_z(n+i) = dft_water_mp(k-n_nuc)%z%A
             C_zf(n+i) = dft_water_mp_f(k-n_nuc)%z%C
             A_zf(n+i) = dft_water_mp_f(k-n_nuc)%z%A
          end if
          D(:,n+i) = dp(:,k-n_nuc)
          Q(:,:,n+i) = qp(:,:,k-n_nuc)
          O(:,:,:,n+i) = op(:,:,:,k-n_nuc)
       end if
    end do
    rmass=rmass/frag_mass

!!$print*,'MULTIPOLES'
!do i=1,n_nuc+n_emp
!!$   print*,i+n,mname(i+n),mp(:,i+n)
!!$   print*,'----------------------------------------------------'
!!$   print*,'CHARGE',Z(i),C_z(i),A_z(i),C_zf(i),A_zf(i)
!!$   print*,'----------------------------------------------------'
!!$   print*,'DIPOLE',D(:,i)
!!$   print*,'----------------------------------------------------'
!!$   print*,'QUAD',Q(:,:,i)
!!$   print*,'----------------------------------------------------'
!!$   print*,'OCTO',O(:,:,:,i)
!!$   print*,'====================================================='
!if((i==2 .or. i==5) .and. i_efp==1) mp(1,i+n)=mp(1,i+n)-0.0001_r8_kind
!end do

    n=(i_efp-1)*n_pol
    do i=1,n_pol
       k=i
       j=j+1
       lmo(:,n+k) = coor(:,j)
       lname(n+i) = dft_water_pol(k)%name
       pol_ten(:,:,n+i) = pt(:,:,k)
    end do
!!$print*,'POLARIZATION'
!do i=1,n_pol
!!$   print*,i+n,lname(i+n),lmo(:,i+n)
!!$   print*,'----------------------------------------------------'
!!$   print*,'POL',pol_ten(:,:,i+n)
!!$   print*,'====================================================='
!if(i==5 .and. i_efp==1 ) lmo(2,i+n)=lmo(2,i+n)-0.0001_r8_kind
!end do

    n=(i_efp-1)*(n_rep_q+n_rep_q)
    do i=1,n_rep_q+n_rep_q
       k=i
       j=j+1
       rp_q(:,n+k) = coor(:,j)
       if(i <= n_rep_q) then
          rname_q(n+i) = dft_water_qrep(k)%name
          C_q(n+i) = dft_water_qrep(k)%C(1)
          A_q(n+i) = dft_water_qrep(k)%A(1)
       else if(i > n_rep_q) then
          rname_q(n+i) = dft_water_qrep(k-n_rep_q)%name
          C_q(n+i) = dft_water_qrep(k-n_rep_q)%C(2)
          A_q(n+i) = dft_water_qrep(k-n_rep_q)%A(2)
       end if
    end do
!print*,'REPULSION QM'
!do i=1,n_rep_q+n_rep_q
!!$   print*,i,rname_q(i+n),rp_q(:,i+n)
!!$   print*,'----------------------------------------------------'
!!$   print*,'REP_QM',C_q(i),A_q(i)
!!$   print*,'====================================================='
!if(i==1 .or. i==4 .and. i_efp==1) rp_q(1,i+n)=rp_q(1,i+n)-0.0001_r8_kind
!end do
    
    n=(i_efp-1)*n_rep_f
    do i=1,n_rep_f
       k=i
       j=j+1
       rp_f(:,n+k) = coor(:,j)
       rname_f(n+i) = dft_water_frep(k)%name
       C_f(:,n+i) = dft_water_frep(k)%C(:)
       A_f(:,n+i) = dft_water_frep(k)%A(:)
    end do
!!$print*,'REPULSION EFP'
!do i=1,n_rep_f
!!$   print*,i+n,rname_f(i+n),rp_f(:,i+n)
!!$   print*,'----------------------------------------------------'
!!$   print*,'REP_EFP',C_f(:,i+n),A_f(:,i+n)
!!$   print*,'====================================================='
!if(i==1 .and. i_efp==1) rp_f(1,i+n)=rp_f(1,i+n)-0.0001_r8_kind
!end do
    
  end subroutine fill_efp_arrays
  !*************************************************************

  !*************************************************************
  subroutine calc_X_points()
    !------------ Modules used -----------------------------------
    use point_dqo_module, only: symm_external_points, moving_X_centers, moving_R_centers
    use point_dqo_module, only: print_X_grad, print_R_grad,output_geometry_ec
    use pointcharge_module, only: moving_pc, print_pc_grad
    use induced_dipoles_module, only: symm_Pol_points, moving_Pol_centers, print_id_grad
    use induced_dipoles_module, only: output_geometry_id
    use efp_rep_module, only: symm_repf_points
    use unique_atom_module, only: N_unique_atoms
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n_points,i,i1,i2,status,k1,k2
    integer(i4_kind), allocatable :: frag_group(:),rc_type(:)
    integer(i4_kind), allocatable :: i_group(:)
    !------------ Executable code ------------------------------------

    allocate(i_group(n_efp*10),stat=status)
    ASSERT(status==0)

    present_X_centers=0

    n_points=(n_nuc+n_emp)*n_efp
    k1=1; k2=n_nuc+n_emp
    do i=1,n_efp
       k1=1+(n_nuc+n_emp)*(i-1)
       k2=1+(n_nuc+n_emp)*i
       i_group(k1:k2)=(/0,0,0,1,1,1,1,1/)
    end do

    allocate(frag_group(n_points),stat=status)
    ASSERT(status==0)

    do i=1,n_efp
       i1=(i-1)*(n_nuc+n_emp)+1
       i2=i1+(n_nuc+n_emp)-1
       frag_group(i1:i2)=i
    end do

    ! Point charges
    if(calc_pc) then
       present_X_centers=present_X_centers+1
       if(.not. efp_fixed) moving_pc=.true.
       print_pc_grad=print_grad_contrib
       call symm_external_points(N_total=n_points, &
                                 P=mp,ig=i_group(1:n_points),ef_grp=frag_group, &
                                 Name=mname,Z=Z,Cz=C_z,Az=A_z, &
                                 Cz_f=C_zf,Az_f=A_zf)
       if(print_efp_centers) call output_geometry_ec("PC")
    end if

    if(calc_pd .or. calc_pq .or. calc_po) then
       if(.not. efp_fixed) moving_X_centers=.true.
       print_X_grad=print_grad_contrib
    end if
    ! Point dipoles
    if(calc_pd) then
       present_X_centers=present_X_centers+1
       call symm_external_points(N_total=n_points, &
                                 P=mp,ig=i_group(1:n_points),ef_grp=frag_group, &
                                 Name=mname,D=D)
       if(print_efp_centers) call output_geometry_ec("PD")
    end if

    ! Point quadrupoles
    if(calc_pq) then
       present_X_centers=present_X_centers+1
       call symm_external_points(N_total=n_points, &
                                 P=mp,ig=i_group(1:n_points),ef_grp=frag_group, &
                                 Name=mname,Q=Q)
       if(print_efp_centers) call output_geometry_ec("PQ")
    end if

    ! Point octopoles
    if(calc_po) then
       present_X_centers=present_X_centers+1
       call symm_external_points(N_total=n_points, &
                                 P=mp,ig=i_group(1:n_points),ef_grp=frag_group, &
                                 Name=mname,O=O)
       if(print_efp_centers) call output_geometry_ec("PO")
    end if

    deallocate(frag_group,stat=status)
    ASSERT(status==0)

    n_points=n_pol*n_efp

    allocate(frag_group(n_points),stat=status)
    ASSERT(status==0)

    do i=1,n_efp
       i1=(i-1)*n_pol+1
       i2=i1+n_pol-1
       frag_group(i1:i2)=i
    end do
    ! Point induced dipoles
    if(calc_pol) then
       if(.not. efp_fixed) moving_Pol_centers=.true.
       print_id_grad=print_grad_contrib
       call symm_Pol_points(n_points,lmo,frag_group,lname,pol_ten)
       if(print_efp_centers) call output_geometry_id()
    end if

    deallocate(frag_group,stat=status)
    ASSERT(status==0)

    ! Repulsive centers
    if(calc_rep) then
       present_X_centers=present_X_centers+1
       if(.not. efp_fixed) moving_R_centers=.true.
       print_R_grad=print_grad_contrib

       if(N_unique_atoms > 0) then
          n_points=2*n_rep_q*n_efp
          k1=1; k2=2*n_rep_q
          do i=1,n_efp
             k1=1+2*n_rep_q*(i-1)
             k2=1+2*n_rep_q*i
             i_group(k1:k2)=(/0,0,0,1,1,1/)
          end do

          allocate(frag_group(n_points),stat=status)
          ASSERT(status==0)

          do i=1,n_efp
             i1=(i-1)*n_rep_q*2+1
             i2=i1+2*n_rep_q-1
             frag_group(i1:i2)=i
          end do

          call symm_external_points(N_total=n_points, &
                                    P=rp_q(:,1:n_points),ig=i_group(1:n_points),ef_grp=frag_group, &
                                    Name=rname_q(1:n_points),C=C_q(1:n_points),A=A_q(1:n_points))

          deallocate(frag_group,stat=status)
          ASSERT(status==0)
       end if

       if(n_efp > 1) then
          n_points=n_rep_f*n_efp

          allocate(frag_group(n_points),rc_type(n_points),stat=status)
          ASSERT(status==0)

          do i=1,n_efp
             i1=(i-1)*n_rep_f+1
             i2=i1+n_rep_f-1
             frag_group(i1:i2)=i
             rc_type(i1:i2)=(/1,2,3,4/)
          end do

          call symm_repf_points(n_points,rp_f,frag_group,rc_type,rname_f,C_f,A_f)

          deallocate(frag_group,rc_type,stat=status)
          ASSERT(status==0)
       end if
    end if

    deallocate(i_group,stat=status)
    ASSERT(status==0)

  end subroutine calc_X_points
  !*************************************************************

  !*************************************************************
  subroutine calc_efield_points
    !------------ Modules used -----------------------------------
    use elec_static_field_module, only: fill_surf_points, surf_points_grad_information
    use induced_dipoles_module, only: N_ipd
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    call fill_surf_points(lmo,N_ipd)
    call surf_points_grad_information()

  end subroutine calc_efield_points
  !*************************************************************

  !*************************************************************
  subroutine calc_efield_points1
    !------------ Modules used -----------------------------------
    use elec_static_field_module, only: fill_surf_points
    use induced_dipoles_module, only: N_ipd
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    call fill_surf_points(lmo,N_ipd)

  end subroutine calc_efield_points1
  !*************************************************************

  !*************************************************************
  subroutine efp_sum_up_grads()
    ! Purpose: summing up contributions of different EFP centers
    !------------ Modules used -----------------------------------
    use unique_atom_module, only: N_unique_atoms
    use pointcharge_module, only: gradient_pc_cartesian,torque_pc_cartesian
    use point_dqo_module, only: gradient_dip_cartesian,torque_dip_cartesian
    use point_dqo_module, only: gradient_quad_cartesian,torque_quad_cartesian
    use point_dqo_module, only: gradient_oct_cartesian, torque_oct_cartesian
    use point_dqo_module, only: gradient_rep_cartesian, torque_rep_cartesian
    use induced_dipoles_module, only: grad_idip_cartesian, torque_idip_cartesian
    use efp_rep_module, only: gradient_repf_cartesian,torque_repf_cartesian
    use iounitadmin_module, only: output_unit
    use operations_module, only: operations_solvation_effect
    use efp_solv_grad_module, only: gradient_mpole_solv_cartesian,torque_mpole_solv_cartesian
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,n,n1,n2,n3,n4,n5,status,n_points_total,n_points
    type(arrmat2),allocatable :: grad_efp_center_cart(:)
    type(arrmat2),allocatable :: torque_efp_center(:)
    real(r8_kind) :: vect(3)
    !------------ Executable code ------------------------------------

    if(n_efp > 0) then
       allocate(grad_efp_cart_total(n_efp,3),torque_efp_cart_total(n_efp,3),stat=status)
       ASSERT(status==0)
       grad_efp_cart_total=zero; torque_efp_cart_total=zero

       if(.not.efp_fixed) then
          n_points=n_efp*n_gx_points
          allocate(grad_efp_cartesian(n_points),stat=status)
          ASSERT(status==0)
          allocate(torque_efp_cartesian(n_points),stat=status)
          ASSERT(status==0)
          do i=1,n_points
             allocate(grad_efp_cartesian(i)%m(3,1),stat=status)
             ASSERT(status==0)
             grad_efp_cartesian(i)%m=zero
             allocate(torque_efp_cartesian(i)%m(3,1),stat=status)
             ASSERT(status==0)
             torque_efp_cartesian(i)%m=zero
          end do
       end if
    end if

    water: if(dftwater_efp) then
       efp_fix: if(.not.efp_fixed) then
       if(n_efp > 0) then
          n_points_total=n_efp*n_water_centers
          allocate(grad_efp_center_cart(n_points_total),stat=status)
          ASSERT(status==0)
          allocate(torque_efp_center(n_points_total),stat=status)
          ASSERT(status==0)
          do i=1,n_points_total
             allocate(grad_efp_center_cart(i)%m(3,1),stat=status)
             ASSERT(status==0)
             grad_efp_center_cart(i)%m=zero
             allocate(torque_efp_center(i)%m(3,1),stat=status)
             ASSERT(status==0)
             torque_efp_center(i)%m=zero
          end do
       end if
       !  pc_n_q+pc_e_q+pd_e_q+pq_e_q+po_e_q+rep_q  O1
       !  pc_n_q+pc_e_q+pd_e_q+pq_e_q+po_e_q+rep_q  H2
       !  pc_n_q+pc_e_q+pd_e_q+pq_e_q+po_e_q+rep_q  H3
       !  pc_e_q+pd_e_q+pq_e_q+po_e_q               B12
       !  pc_e_q+pd_e_q+pq_e_q+po_e_q               B13
       !  pol_q                                     LMO1
       !  pol_q                                     LMO2
       !  pol_q                                     LMO3
       !  pol_q                                     LMO4
       !  pol_q                                     LMO5
       !                                            CMS
       k=0
       do i=1,n_efp
          n =(i-1)*n_water_centers
          n1=(i-1)*(n_nuc+n_emp)
          n2=(i-1)*(n_rep_q+n_rep_q)
          n4=(i-1)*n_rep_f
          n3=(i-1)*n_pol

          !collect Cartesian gradients on EFP centers
          if(calc_pc) then
             grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                  gradient_pc_cartesian(1+n1)%m(:,1) + gradient_pc_cartesian(1+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                  gradient_pc_cartesian(2+n1)%m(:,1) + gradient_pc_cartesian(2+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                  gradient_pc_cartesian(3+n1)%m(:,1) + gradient_pc_cartesian(3+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(4+n)%m(:,1) = grad_efp_center_cart(4+n)%m(:,1) + &
                  gradient_pc_cartesian(4+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(5+n)%m(:,1) = grad_efp_center_cart(5+n)%m(:,1) + &
                  gradient_pc_cartesian(5+n_nuc+n1)%m(:,1)

             torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                  torque_pc_cartesian(1+n1)%m(:,1) + torque_pc_cartesian(1+n_nuc+n1)%m(:,1)
             torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                  torque_pc_cartesian(2+n1)%m(:,1) + torque_pc_cartesian(2+n_nuc+n1)%m(:,1)
             torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                  torque_pc_cartesian(3+n1)%m(:,1) + torque_pc_cartesian(3+n_nuc+n1)%m(:,1)
             torque_efp_center(4+n)%m(:,1) = torque_efp_center(4+n)%m(:,1) + &
                  torque_pc_cartesian(4+n_nuc+n1)%m(:,1)
             torque_efp_center(5+n)%m(:,1) = torque_efp_center(5+n)%m(:,1) + &
                  torque_pc_cartesian(5+n_nuc+n1)%m(:,1)
          end if
          if(calc_pd) then
             grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                  gradient_dip_cartesian(1+n_nuc+n1)%m(:,1)+gradient_dip_cartesian(1+n1)%m(:,1)
             grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                  gradient_dip_cartesian(2+n_nuc+n1)%m(:,1)+gradient_dip_cartesian(2+n1)%m(:,1)
             grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                  gradient_dip_cartesian(3+n_nuc+n1)%m(:,1)+gradient_dip_cartesian(3+n1)%m(:,1)
             grad_efp_center_cart(4+n)%m(:,1) = grad_efp_center_cart(4+n)%m(:,1) + &
                  gradient_dip_cartesian(4+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(5+n)%m(:,1) = grad_efp_center_cart(5+n)%m(:,1) + &
                  gradient_dip_cartesian(5+n_nuc+n1)%m(:,1)

             torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                  torque_dip_cartesian(1+n_nuc+n1)%m(:,1)+torque_dip_cartesian(1+n1)%m(:,1)
             torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                  torque_dip_cartesian(2+n_nuc+n1)%m(:,1)+torque_dip_cartesian(2+n1)%m(:,1)
             torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                  torque_dip_cartesian(3+n_nuc+n1)%m(:,1)+torque_dip_cartesian(3+n1)%m(:,1)
             torque_efp_center(4+n)%m(:,1) = torque_efp_center(4+n)%m(:,1) + &
                  torque_dip_cartesian(4+n_nuc+n1)%m(:,1)
             torque_efp_center(5+n)%m(:,1) = torque_efp_center(5+n)%m(:,1) + &
                  torque_dip_cartesian(5+n_nuc+n1)%m(:,1)
          end if
          if(calc_pq) then
             grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                  gradient_quad_cartesian(1+n_nuc+n1)%m(:,1)+gradient_quad_cartesian(1+n1)%m(:,1)
             grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                  gradient_quad_cartesian(2+n_nuc+n1)%m(:,1)+gradient_quad_cartesian(2+n1)%m(:,1)
             grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                  gradient_quad_cartesian(3+n_nuc+n1)%m(:,1)+gradient_quad_cartesian(3+n1)%m(:,1)
             grad_efp_center_cart(4+n)%m(:,1) = grad_efp_center_cart(4+n)%m(:,1) + &
                  gradient_quad_cartesian(4+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(5+n)%m(:,1) = grad_efp_center_cart(5+n)%m(:,1) + &
                  gradient_quad_cartesian(5+n_nuc+n1)%m(:,1)

             torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                  torque_quad_cartesian(1+n_nuc+n1)%m(:,1)+torque_quad_cartesian(1+n1)%m(:,1)
             torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                  torque_quad_cartesian(2+n_nuc+n1)%m(:,1)+torque_quad_cartesian(2+n1)%m(:,1)
             torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                  torque_quad_cartesian(3+n_nuc+n1)%m(:,1)+torque_quad_cartesian(3+n1)%m(:,1)
             torque_efp_center(4+n)%m(:,1) = torque_efp_center(4+n)%m(:,1) + &
                  torque_quad_cartesian(4+n_nuc+n1)%m(:,1)
             torque_efp_center(5+n)%m(:,1) = torque_efp_center(5+n)%m(:,1) + &
                  torque_quad_cartesian(5+n_nuc+n1)%m(:,1)
          end if
          if(calc_po) then
             grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                  gradient_oct_cartesian(1+n_nuc+n1)%m(:,1)+gradient_oct_cartesian(1+n1)%m(:,1)
             grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                  gradient_oct_cartesian(2+n_nuc+n1)%m(:,1)+gradient_oct_cartesian(2+n1)%m(:,1)
             grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                  gradient_oct_cartesian(3+n_nuc+n1)%m(:,1)+gradient_oct_cartesian(3+n1)%m(:,1)
             grad_efp_center_cart(4+n)%m(:,1) = grad_efp_center_cart(4+n)%m(:,1) + &
                  gradient_oct_cartesian(4+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(5+n)%m(:,1) = grad_efp_center_cart(5+n)%m(:,1) + &
                  gradient_oct_cartesian(5+n_nuc+n1)%m(:,1)

             torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                  torque_oct_cartesian(1+n_nuc+n1)%m(:,1)+torque_oct_cartesian(1+n1)%m(:,1)
             torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                  torque_oct_cartesian(2+n_nuc+n1)%m(:,1)+torque_oct_cartesian(2+n1)%m(:,1)
             torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                  torque_oct_cartesian(3+n_nuc+n1)%m(:,1)+torque_oct_cartesian(3+n1)%m(:,1)
             torque_efp_center(4+n)%m(:,1) = torque_efp_center(4+n)%m(:,1) + &
                  torque_oct_cartesian(4+n_nuc+n1)%m(:,1)
             torque_efp_center(5+n)%m(:,1) = torque_efp_center(5+n)%m(:,1) + &
                  torque_oct_cartesian(5+n_nuc+n1)%m(:,1)
          end if

          if(operations_solvation_effect) then
             grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                  gradient_mpole_solv_cartesian(1+n_nuc+n1)%m(:,1)+gradient_mpole_solv_cartesian(1+n1)%m(:,1)
             grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                  gradient_mpole_solv_cartesian(2+n_nuc+n1)%m(:,1)+gradient_mpole_solv_cartesian(2+n1)%m(:,1)
             grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                  gradient_mpole_solv_cartesian(3+n_nuc+n1)%m(:,1)+gradient_mpole_solv_cartesian(3+n1)%m(:,1)
             grad_efp_center_cart(4+n)%m(:,1) = grad_efp_center_cart(4+n)%m(:,1) + &
                  gradient_mpole_solv_cartesian(4+n_nuc+n1)%m(:,1)
             grad_efp_center_cart(5+n)%m(:,1) = grad_efp_center_cart(5+n)%m(:,1) + &
                  gradient_mpole_solv_cartesian(5+n_nuc+n1)%m(:,1)

             torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                  torque_mpole_solv_cartesian(1+n_nuc+n1)%m(:,1)+torque_mpole_solv_cartesian(1+n1)%m(:,1)
             torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                  torque_mpole_solv_cartesian(2+n_nuc+n1)%m(:,1)+torque_mpole_solv_cartesian(2+n1)%m(:,1)
             torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                  torque_mpole_solv_cartesian(3+n_nuc+n1)%m(:,1)+torque_mpole_solv_cartesian(3+n1)%m(:,1)
             torque_efp_center(4+n)%m(:,1) = torque_efp_center(4+n)%m(:,1) + &
                  torque_mpole_solv_cartesian(4+n_nuc+n1)%m(:,1)
             torque_efp_center(5+n)%m(:,1) = torque_efp_center(5+n)%m(:,1) + &
                  torque_mpole_solv_cartesian(5+n_nuc+n1)%m(:,1)
          end if

          if(calc_rep) then
             if(N_unique_atoms > 0) then
                grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                     gradient_rep_cartesian(1+n2)%m(:,1) + gradient_rep_cartesian(1+n_rep_q+n2)%m(:,1)
                grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                     gradient_rep_cartesian(2+n2)%m(:,1) + gradient_rep_cartesian(2+n_rep_q+n2)%m(:,1)
                grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                     gradient_rep_cartesian(3+n2)%m(:,1) + gradient_rep_cartesian(3+n_rep_q+n2)%m(:,1)

                torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                     torque_rep_cartesian(1+n2)%m(:,1) + torque_rep_cartesian(1+n_rep_q+n2)%m(:,1)
                torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                     torque_rep_cartesian(2+n2)%m(:,1) + torque_rep_cartesian(2+n_rep_q+n2)%m(:,1)
                torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                     torque_rep_cartesian(3+n2)%m(:,1) + torque_rep_cartesian(3+n_rep_q+n2)%m(:,1)
             end if
             if(n_efp > 1) then
                grad_efp_center_cart(1+n)%m(:,1) = grad_efp_center_cart(1+n)%m(:,1) + &
                     gradient_repf_cartesian(1+n4)%m(:,1)
                grad_efp_center_cart(2+n)%m(:,1) = grad_efp_center_cart(2+n)%m(:,1) + &
                     gradient_repf_cartesian(2+n4)%m(:,1)
                grad_efp_center_cart(3+n)%m(:,1) = grad_efp_center_cart(3+n)%m(:,1) + &
                     gradient_repf_cartesian(3+n4)%m(:,1)
                grad_efp_center_cart(11+n)%m(:,1) = grad_efp_center_cart(11+n)%m(:,1) + &
                     gradient_repf_cartesian(4+n4)%m(:,1)

                torque_efp_center(1+n)%m(:,1) = torque_efp_center(1+n)%m(:,1) + &
                     torque_repf_cartesian(1+n4)%m(:,1)
                torque_efp_center(2+n)%m(:,1) = torque_efp_center(2+n)%m(:,1) + &
                     torque_repf_cartesian(2+n4)%m(:,1)
                torque_efp_center(3+n)%m(:,1) = torque_efp_center(3+n)%m(:,1) + &
                     torque_repf_cartesian(3+n4)%m(:,1)
                torque_efp_center(11+n)%m(:,1) = torque_efp_center(11+n)%m(:,1) + &
                     torque_repf_cartesian(4+n4)%m(:,1)
             end if
          end if
          if(calc_pol) then
             grad_efp_center_cart(6+n)%m(:,1) =grad_efp_center_cart(6+n)%m(:,1) +grad_idip_cartesian(1+n3)%m(:,1)
             grad_efp_center_cart(7+n)%m(:,1) =grad_efp_center_cart(7+n)%m(:,1) +grad_idip_cartesian(2+n3)%m(:,1)
             grad_efp_center_cart(8+n)%m(:,1) =grad_efp_center_cart(8+n)%m(:,1) +grad_idip_cartesian(3+n3)%m(:,1)
             grad_efp_center_cart(9+n)%m(:,1) =grad_efp_center_cart(9+n)%m(:,1) +grad_idip_cartesian(4+n3)%m(:,1)
             grad_efp_center_cart(10+n)%m(:,1)=grad_efp_center_cart(10+n)%m(:,1)+grad_idip_cartesian(5+n3)%m(:,1)

             torque_efp_center(6+n)%m(:,1) =torque_efp_center(6+n)%m(:,1) +torque_idip_cartesian(1+n3)%m(:,1)
             torque_efp_center(7+n)%m(:,1) =torque_efp_center(7+n)%m(:,1) +torque_idip_cartesian(2+n3)%m(:,1)
             torque_efp_center(8+n)%m(:,1) =torque_efp_center(8+n)%m(:,1) +torque_idip_cartesian(3+n3)%m(:,1)
             torque_efp_center(9+n)%m(:,1) =torque_efp_center(9+n)%m(:,1) +torque_idip_cartesian(4+n3)%m(:,1)
             torque_efp_center(10+n)%m(:,1)=torque_efp_center(10+n)%m(:,1)+torque_idip_cartesian(5+n3)%m(:,1)
          end if

          if(print_grad_centers) then
             write(output_unit,'(a22,i4,a37)') 'Water fragment number ',i, &
                  ',gradients and torques on EFP centers'
             write(output_unit,'(a40)') '........................................'
             do j=1,n_water_centers
                write(output_unit,'(i3,3F15.10,a3,3F15.10)') j, &
                     grad_efp_center_cart(j+n)%m(:,1),' / ', &
                     torque_efp_center(j+n)%m(:,1)
             end do
             write(output_unit,'(a1)') ' '
          end if

          !Destribution of efp gradients between the efp atoms
          n5=(i-1)*n_gx_points
          do j=1,n_water_centers
             grad_efp_cart_total(i,:)=grad_efp_cart_total(i,:)+grad_efp_center_cart(j+n)%m(:,1)
             torque_efp_cart_total(i,:)=torque_efp_cart_total(i,:)+torque_efp_center(j+n)%m(:,1)

             if(print_efp_grad) then
                vect=grad_efp_center_cart(j+n)%m(:,1)
                grad_efp_cartesian(1+n5)%m(:,1)=grad_efp_cartesian(1+n5)%m(:,1)+matmul(vect,dRi_dRj(j,1,i))
!!$if(j>=6 .and. j<= 10 .and. i==1) print*,i,j,'1',dRi_dRj(j,1,i)
!!$if(j>=6 .and. j<= 10 .and. i==1) print*,i,j,'2',dRi_dRj(j,2,i)
!!$if(j>=6 .and. j<= 10 .and. i==1) print*,i,j,'3',dRi_dRj(j,3,i)
                grad_efp_cartesian(2+n5)%m(:,1)=grad_efp_cartesian(2+n5)%m(:,1)+matmul(vect,dRi_dRj(j,2,i))
                grad_efp_cartesian(3+n5)%m(:,1)=grad_efp_cartesian(3+n5)%m(:,1)+matmul(vect,dRi_dRj(j,3,i))
             end if
          end do
       end do

       if(n_efp > 0) then
          do i=1,n_points_total
             deallocate(grad_efp_center_cart(i)%m,stat=status)
             ASSERT(status==0)
             deallocate(torque_efp_center(i)%m,stat=status)
             ASSERT(status==0)
          end do
          deallocate(grad_efp_center_cart,stat=status)
          ASSERT(status==0)
          deallocate(torque_efp_center,stat=status)
          ASSERT(status==0)
       end if
       endif efp_fix
    end if water

  end subroutine efp_sum_up_grads
  !*************************************************************

  !*************************************************************
  function dRi_dRj(i,j,i_efp) result(deriv)
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: i,j,i_efp
    real(r8_kind) :: deriv(3,3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: k,n
    real(r8_kind), parameter :: e1(3)=(/one,zero,zero/) 
    real(r8_kind), parameter :: e2(3)=(/zero,one,zero/) 
    real(r8_kind), parameter :: e3(3)=(/zero,zero,one/)
    real(r8_kind) :: r1(3),r2(3)
    !------------ Executable code ------------------------------------

    ASSERT(j<=3)
    ASSERT(i<=n_water_centers)

    n=(i_efp-1)*(n_nuc+n_emp)

    
    r1=mp(:,2+n)-mp(:,1+n); r2=mp(:,3+n)-mp(:,1+n)

    deriv=zero

    !dr1/dr1, dr2/dr2, dr3/dr3
    if(i==j) then 
       do k=1,3
          deriv(k,k)=one
       end do
    end if
    !dr4/dr1, dr4/dr2
    if(i==4) then
       if(j==1 .or. j==2) then
          do k=1,3
             deriv(k,k)=half
          end do
       end if
    end if
    !dr5/dr1, dr5/dr3
    if(i==5) then
       if(j==1 .or. j==3) then
          do k=1,3
             deriv(k,k)=half
          end do
       end if
    end if
    !dr6/dr1, dr6/dr2, dr6/dr3
    if(i==6) then
       if(j==1) then
          do k=1,3
             deriv(k,k)=(one-C1)
          end do
       end if
       if(j==2 .or. j==3) then
          do k=1,3
             deriv(k,k)=C1/two
          end do
       end if
    end if
    !dr7/dr1, dr7/dr2, dr7/dr3
    if(i==7) then
        if(j==1) then
          do k=1,3
             deriv(k,k)=(one-C2)
          end do
       end if
       if(j==2) then
          do k=1,3
             deriv(k,k)=(C2-C3)/two
          end do
       end if
       if(j==3) then
          do k=1,3
             deriv(k,k)=(C2+C3)/two
          end do
       end if
    end if
    !dr8/dr1, dr8/dr2, dr8/dr3
    if(i==8) then
        if(j==1) then
          do k=1,3
             deriv(k,k)=(one-C2)
          end do
       end if
       if(j==2) then
          do k=1,3
             deriv(k,k)=(C2+C3)/two
          end do
       end if
       if(j==3) then
          do k=1,3
             deriv(k,k)=(C2-C3)/two
          end do
       end if
    end if
    !dr9/dr1, dr9/dr2, dr9/dr3
    if(i==9) then
       if(j==1) then
          do k=1,3
             deriv(k,k)=(one+C4)
          end do
          deriv(:,1)=deriv(:,1)+C5*(vector_product(e1,r2)+vector_product(r1,e1))
          deriv(:,2)=deriv(:,2)+C5*(vector_product(e2,r2)+vector_product(r1,e2))
          deriv(:,3)=deriv(:,3)+C5*(vector_product(e3,r2)+vector_product(r1,e3))
       end if
       if(j==2) then
          do k=1,3
             deriv(k,k)=-C4/two
          end do
          deriv(:,1)=deriv(:,1)-C5*vector_product(e1,r2)
          deriv(:,2)=deriv(:,2)-C5*vector_product(e2,r2)
          deriv(:,3)=deriv(:,3)-C5*vector_product(e3,r2)
       end if
       if(j==3) then
          do k=1,3
             deriv(k,k)=-C4/two
          end do
          deriv(:,1)=deriv(:,1)-C5*vector_product(r1,e1)
          deriv(:,2)=deriv(:,2)-C5*vector_product(r1,e2)
          deriv(:,3)=deriv(:,3)-C5*vector_product(r1,e3)
       end if
    end if
    !dr10/dr1, dr10/dr2, dr10/dr3
    if(i==10) then
       if(j==1) then
          do k=1,3
             deriv(k,k)=(one+C4)
          end do
          deriv(:,1)=deriv(:,1)-C5*(vector_product(e1,r2)+vector_product(r1,e1))
          deriv(:,2)=deriv(:,2)-C5*(vector_product(e2,r2)+vector_product(r1,e2))
          deriv(:,3)=deriv(:,3)-C5*(vector_product(e3,r2)+vector_product(r1,e3))
       end if
       if(j==2) then
          do k=1,3
             deriv(k,k)=-C4/two
          end do
          deriv(:,1)=deriv(:,1)+C5*vector_product(e1,r2)
          deriv(:,2)=deriv(:,2)+C5*vector_product(e2,r2)
          deriv(:,3)=deriv(:,3)+C5*vector_product(e3,r2)
       end if
       if(j==3) then
          do k=1,3
             deriv(k,k)=-C4/two
          end do
          deriv(:,1)=deriv(:,1)+C5*vector_product(r1,e1)
          deriv(:,2)=deriv(:,2)+C5*vector_product(r1,e2)
          deriv(:,3)=deriv(:,3)+C5*vector_product(r1,e3)
       end if
    end if
    !dr11/dr1, dr11/dr2, dr11/dr3
    if(i==11) then
       if(j==1 .or. j==2 .or. j==3) then
          do k=1,3
             deriv(k,k)=dft_water_nuc(j)%mass/ &
                  (dft_water_nuc(1)%mass+dft_water_nuc(2)%mass+dft_water_nuc(3)%mass)
          end do
       end if
    end if

    return

  end function dRi_dRj
  !*************************************************************

  !*************************************************************
  function efp_max_grad() result(max_grad)
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: max_grad
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j
    !------------ Executable code ------------------------------------

   max_grad =zero

    do i=1,n_efp
       do j=1,3
          if(abs(grad_efp_cart_total(i,j)) > max_grad) max_grad=abs(grad_efp_cart_total(i,j))
          if(abs(torque_efp_cart_total(i,j)) > max_grad) max_grad=abs(torque_efp_cart_total(i,j))
       end do
    end do

  end function efp_max_grad
  !*************************************************************

  !*************************************************************
  subroutine efp_grad_cart_write()
    !------------ Modules used -----------------------------------
    use iounitadmin_module, only: output_unit
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,n,n5
    !------------ Executable code ------------------------------------

    if(print_efp_grad) then
       write(output_unit,'(/A)') 'EFP Cartesian gradients'
       write(output_unit,'(A)')  '----------------------------------'
       do i=1,n_efp
          n=(i-1)*(n_nuc+n_emp)
          n5=(i-1)*n_gx_points
          write(output_unit,'(a22,i4)') 'Water fragment number ',i
          write(output_unit,'(a26)') '..........................'
          do j=1,n_gx_points
             k=j
             if(print_as_in_input) k=num(j,i)
             write(output_unit,'(a4,3F20.12)') mname(n+k),grad_efp_cartesian(k+n5)%m(:,1)
          end do
          write(output_unit,'(a1)') ' '
       end do
    end if

    write(output_unit,'(/A)') 'Total EFP Cartesian gradients and torques'
    write(output_unit,'(A)')  '-----------------------------------------'
    do i=1,n_efp
       write(output_unit,'(a22,i4)') 'Water fragment number ',i
       write(output_unit,'(a40)') '........................................'
       write(output_unit,'(3F15.10,a3,3F15.10)') &
            grad_efp_cart_total(i,:),' / ',torque_efp_cart_total(i,:)
       write(output_unit,'(a1)') ' '
    end do
    write(output_unit,'(A)')  '-----------------------------------------'

  end subroutine efp_grad_cart_write
  !*************************************************************

  !*************************************************************
  subroutine efp_write_gxfile(energy,grad_cart)
    !------------ Modules used -----------------------------------
    use iounitadmin_module, only: openget_iounit,returnclose_iounit
    use filename_module, only: inpfile
    use unique_atom_module, only: N_unique_atoms, unique_atoms
    use point_dqo_module, only: pd_array,N_pd
    use induced_dipoles_module, only: ipd_array, N_ipd
!!$    use energy_calc_module, only: get_energy
!!$    use gradient_data_module, only: gradient_cartesian
!!$    use efp_efp_module, only: qm_efp_energy
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(in) :: energy
    type(arrmat2), optional, intent(in) :: grad_cart(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,k,l,status,n1,natoms,n_particls,l1
    real(r8_kind), allocatable :: at_num(:),coor(:,:)
    integer(i4_kind), allocatable :: iuniq(:),inum(:),zmat(:,:),numx(:,:)
    character(len=5) ::  gx_buff
!!$    real(r8_kind) :: energy
    real(r8_kind) :: grad(3)
    !------------ Executable code ------------------------------------

    allocate(at_num(gx_line-1),coor(3,gx_line-1),iuniq(gx_line-1),inum(gx_line-1), &
         zmat(3,gx_line-1),numx(3,gx_line-1),stat=status)
    ASSERT(status==0)

    if(ex_gx) then
    rewind(gx_id)

    do i=1,gx_line-1
       read(gx_id,*,err=100,end=101) at_num(i),coor(:,i),iuniq(i),inum(i),zmat(:,i),numx(:,i)
    end do
    zmat=0; numx=0

    j=0; l=0
    do i=1,N_unique_atoms
       k=0
       equals: do 
          j=j+1
          if (at_num(j)==99.0_r8_kind .or. iuniq(j)==0) cycle equals
          l=l+1
          k=k+1
          coor(:,l) = unique_atoms(i)%position(:,k)
          at_num(l)=at_num(j)
          iuniq(l)=i
          inum(l)=l
          if (k==unique_atoms(i)%N_equal_atoms) exit equals
       enddo equals
    enddo
    natoms=l

    do i=1,n_efp
       if(dftwater_efp) then
          n1=(i-1)*(n_nuc+n_emp)

1         j=j+1
          if(at_num(j) == 99.0_r8_kind .or. iuniq(j) == 0) goto 1
          k=num(1,i)
          l=l+1
          at_num(l)=at_num(j)
          iuniq(l)=natoms+i
          inum(l)=natoms+i
          coor(:,l)=mp(:,k+n1)
2         j=j+1
          if(at_num(j) == 99.0_r8_kind .or. iuniq(j) == 0) goto 2
          k=num(2,i)
          l=l+1
          at_num(l)=at_num(j)
          iuniq(l)=natoms+i
          inum(l)=natoms+i
          coor(:,l)=mp(:,k+n1)
3         j=j+1
          if(at_num(j) == 99.0_r8_kind .or. iuniq(j) == 0) goto 3
          k=num(3,i)
          l=l+1
          at_num(l)=at_num(j)
          iuniq(l)=natoms+i
          inum(l)=natoms+i
          coor(:,l)=mp(:,k+n1)
       end if
    end do
    n_particls=l
    ASSERT(j==gx_line-1)

    rewind(gx_id)
    else
       zmat=0; numx=0

       j=0; l=0
       do i=1,N_unique_atoms
          do j=1,unique_atoms(i)%N_equal_atoms 
             l=l+1
             coor(:,l) = unique_atoms(i)%position(:,j)
             at_num(l)= unique_atoms(i)%Z
             iuniq(l)=i
             inum(l)=l
          enddo
       enddo
       natoms=l

       do i=1,n_efp
          if(dftwater_efp) then
             n1=(i-1)*(n_nuc+n_emp)
             do j=1,3
                k=num(j,i)
                l=l+1
                if(k==1) at_num(l)=8.01
                if(k==2) at_num(l)=1.02
                if(k==3) at_num(l)=1.03
                iuniq(l)=natoms+i
                inum(l)=natoms+i
                coor(:,l)=mp(:,k+n1)
             end do
          end if
       end do
       n_particls=l

       gx_id = openget_iounit(file=inpfile('gxfile'), status='new', &
            form='formatted')
    endif

    do i=1,n_particls
      write(gx_id,1020) at_num(i),coor(:,i),iuniq(i),inum(i),zmat(:,i),numx(:,i)
    end do
    write(gx_id,1031) iwork, zero,zero,zero,0,0,0,0,0,0,0,0

!!$    if(qm_fixed) then
!!$       energy=qm_efp_energy
!!$    else
!!$       call get_energy(tot=energy)
!!$    end if
    write(gx_id,'(2F24.12,2x,3I5)') energy,energy,gx_line-1,gx_line-1,gx_line-1

    l=0
    do i=1,n_unique_atoms   
       j = unique_atoms(i)%moving_atom
       do k=1,unique_atoms(i)%n_equal_atoms
          if (j == 0 .or. qm_fixed) then
             grad = zero
          else
!!$             grad = gradient_cartesian(j)%m(:,k)
             grad = grad_cart(j)%m(:,k)
          endif
          l=l+1
          write(gx_id,'(I5,5x,3F17.12)') l,grad
       enddo
    end do
    ASSERT(l==natoms)

    l1=l
    do i=1,n_efp
       l=l+1
       grad=grad_efp_cart_total(i,:)
       write(gx_id,'(I5,5x,3F17.12)') l,grad
       grad=torque_efp_cart_total(i,:)
       write(gx_id,'(I5,5x,3F17.12)') l,grad
    end do

    call returnclose_iounit(gx_id)

    deallocate(at_num,coor,iuniq,inum,zmat,numx,stat=status)
    ASSERT(status==0)

    if(.not. efp_fixed) then
       do i=1,n_efp*n_gx_points
          deallocate(grad_efp_cartesian(i)%m,stat=status)
          ASSERT(status==0)
          deallocate(torque_efp_cartesian(i)%m,stat=status)
          ASSERT(status==0)
       end do
       deallocate(grad_efp_cartesian,stat=status)
       ASSERT(status==0)
       deallocate(torque_efp_cartesian,stat=status)
       ASSERT(status==0)
    end if
    deallocate(grad_efp_cart_total,torque_efp_cart_total,stat=status)
    ASSERT(status==0)

    if(.not. qm_fixed) call shutdown_efp_arrays()

    return
    
100 call error_handler("efp_module: efp_write_gxfile: attempt to read after end of GX file")
101 write(gx_buff,'(i5)') i
    call error_handler("efp_module: efp_write_gxfile: wrong GX file. Line - "//trim(gx_buff))

1020  format((f7.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4))
1031  format(f7.1,3(2X,f21.12),2i4,2x,3I4,2X,3I4,i5)

  end subroutine efp_write_gxfile
  !*************************************************************
  
  !*************************************************************
  subroutine update_efps(x,y,z,ds)
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(inout) :: x(:),y(:),z(:) ! Coordinates of EFP centers
                                                   ! Size == n_efp
    real(r8_kind), intent(in)    :: ds(:)          ! Displacements of EFP centers
                                                   ! Include translations and rotations
                                                   ! Size == 6*n_efp
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,l,k,m
    real(r8_kind) :: rotmat(3,3),r_tmp(3)
    !------------ Executable code ------------------------------------

    ASSERT(size(x)==3*n_efp)
    ASSERT(size(ds)==6*n_efp)

    m=0
    l=1
    do i=1,n_efp
       k=l+1
       call calc_rotmat(ds(3*k-2),ds(3*k-1),ds(3*k),rotmat)

       do j=1,3
          m=m+1
          r_tmp(1)=x(m)-rcm(1,i); r_tmp(2)=y(m)-rcm(2,i); r_tmp(3)=z(m)-rcm(3,i)
          r_tmp=matmul(rotmat,r_tmp)
          x(m)=r_tmp(1)+rcm(1,i); y(m)=r_tmp(2)+rcm(2,i); z(m)=r_tmp(3)+rcm(3,i)
          x(m)=x(m)+ds(3*l-2)   ; y(m)=y(m)+ds(3*l-1)   ; z(m)=z(m)+ds(3*l)
       end do

       l=l+2
    end do

  contains
    subroutine calc_rotmat(dx,dy,dz,rmat)
      !------------ Modules used -----------------------------------
      !------------ Declaration of formal parameters ---------------
      real(r8_kind), intent(in)  :: dx,dy,dz
      real(r8_kind), intent(out) :: rmat(3,3)
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      real(r8_kind) :: dx2,dy2,dz2,dxdy,dxdz,dydz,dr,dr2,r1,r2
      !------------ Executable code --------------------------------

      dx2 =dx*dx; dy2 =dy*dy; dz2 =dz*dz
      dxdy=dx*dy; dxdz=dx*dy; dydz=dy*dz
      
      dr2=dx2+dy2+dz2; dr=sqrt(dr2)

      if(dr < small) then
         r1=one; r2=half
      else
         r1=sin(dr)/dr; r2=(one-cos(dr))/dr2
      end if

      rmat(1,1)= one-(dy2+dz2)*r2
      rmat(1,2)=-dz*r1+dxdy*r2
      rmat(1,3)= dy*r1+dxdz*r2
      rmat(2,1)= dz*r1+dxdy*r2
      rmat(2,2)= one-(dx2+dz2)*r2
      rmat(2,3)=-dx*r1+dydz*r2
      rmat(3,1)=-dy*r1+dxdz*r2
      rmat(3,2)= dx*r1+dydz*r2
      rmat(3,3)= one-(dx2+dy2)*r2

    end subroutine calc_rotmat
  end subroutine update_efps
  !*************************************************************
!--------------- End of module ----------------------------------
end module efp_module

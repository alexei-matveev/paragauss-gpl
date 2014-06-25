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
!===============================================================
! Public interface of module
!===============================================================
module induced_dipoles_module
  !---------------------------------------------------------------
  !
  !  Purpose: induced point dipoles
  !           (accounting for polarization of eXternal centers)
  !
  !
  !  Module called by:
  !
  !
  !  References: ...
  !
  !
  !  Author: AS
  !  Date: 24.08.2007
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
#include <def.h>

  use type_module ! type specification parameters
  use datatype
  use unique_atom_module, only: unique_atom_type
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use pointcharge_module, only: rcm
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------
  type, public :: ind_dipole_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_dipoles
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind), pointer :: pol_tensor(:,:,:)
     real(r8_kind), pointer :: idipole(:,:)
     real(r8_kind), pointer :: idipole1(:,:)
     integer(i4_kind), pointer :: group(:)
  end type ind_dipole_type

  !------------ Declaration of constants and variables ------------
  integer(i4_kind), public :: N_ipd=0_i4_kind
  type(ind_dipole_type), pointer, public :: ipd_array(:)
  real(r8_kind),public,allocatable :: totsym_ind_dip(:)
  type(arrmat2),public,allocatable :: ham_Pol(:)
  real(r8_kind),public :: en_ind_dipole

  type(arrmat3), pointer, public :: ind_dip_grad_info(:)
  integer(i4_kind), allocatable, public :: ind_dip_grad_index(:)

  logical, public :: moving_Pol_centers=.false.
  integer(i4_kind), public :: totsym_grad_idip_length=0
  real(r8_kind),public,allocatable,target :: grad_idip_totalsym(:)
  type(arrmat2),public,allocatable :: grad_idip_cartesian(:)
  real(r8_kind),public,allocatable,target :: torque_idip_totalsym(:)
  type(arrmat2),public,allocatable :: torque_idip_cartesian(:)
  logical, public :: print_id_grad=.false.
  logical, public :: do_iterations, do_pol_pcm
  integer(i4_kind), public :: n_update

  !------------ Interface statements -----------------------------

  !------------ public functions and subroutines ------------------
  public :: dealloc_id
  public :: dealloc_Pol_center_inform
  public :: symm_Pol_points
  public :: calc_Pol_centers
  public :: output_geometry_id
  public :: pol_centers_bcast
  public :: init_Pol_centers_grads
  public :: Pol_centers_grads_shutdown
  public :: totsym_id_grad_unpack
  public :: totsym_id_grad_pack
  public :: send_receive_id
  public :: free_field_arrays
  public :: send_receive_id1

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------
  !------------ Declaration of constants and variables ------------
  type(unique_atom_type), pointer :: unique_Pol_points(:)
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  function calc_Pol_centers() result (do_pol)
    !Calculate or not QM-EFP polarizabilities
    !** End of interface *****************************************

    logical :: do_pol

    do_pol=.false.
!!$    if(N_ipd > 0 .and. N_unique_atoms > 0) do_pol=.true.
    if(N_ipd > 0) do_pol=.true.

  end function calc_Pol_centers
  !*************************************************************

  !*************************************************************
  subroutine symm_Pol_points (N_total, P, ef_grp, name, Alpha)
    !
    ! Symmetrize polarizable point centers (dipoles)
    !
    use group_module, only: symm_transformation_int, sub_group, &
         group_coset, ylm_trafos, group_num_el, group_coset_decomp
    use symm_module, only: symm_adapt_centers
    implicit none
    integer(i4_kind) :: N_total
    real(r8_kind) :: P(3,N_total)
    integer(i4_kind) :: ef_grp(N_total)
    character(len=12)  :: name(N_total)
    real(r8_kind) :: Alpha(3,3,N_total)
    !** End of interface *****************************************

    real(r8_kind),parameter :: small = 1.e-10_r8_kind
    ! very small value
    !------------ Declaration of local variables -----------------
    integer (i4_kind) :: i, j, k, n, l, alloc_stat, n_size
    ! counters
    integer(i4_kind) :: n_equal,n_equal_chck,N_eq
    ! number of equal atoms
    real(r8_kind) :: position(3),position2(3),position3(3)
    ! coordinates of unique atom in the order (x,z,y)
    type(symm_transformation_int),allocatable :: point_trafos_buf(:)
    type(symm_transformation_int), pointer :: point_trafos_id(:)
    type(sub_group) :: local_groups_id
    ! local symmetry group of the unique external centers
    type(group_coset) :: cosets_id
    ! coset of the local symmetry group of the unique external centers
    logical, allocatable :: help_dim(:)
    type(ind_dipole_type), allocatable :: buffer_id(:)
    logical :: do_alloc
    !------------ Executable code --------------------------------

    allocate(help_dim(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    help_dim=.false.

    allocate(point_trafos_buf(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    allocate(buffer_id(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    n_size=0
    l=0
    i_N_tot: do i=1,N_total
       if(help_dim(i)) cycle i_N_tot
       n_size=n_size+1

       ! reorder coordinates of unique atom as
       ! (x,y,z) --> (x,z,y) in order to comply with the
       ! convention for angular momentum l=1
       position(1) = P(1,i)
       position(2) = P(3,i)
       position(3) = P(2,i)

       !
       ! determine local symmetry groups
       !

       ! now apply all symmetry operations to the position of the
       ! unique atom
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
          end if
       end do

       ! allocate group elements
       local_groups_id%num_el = n
       allocate(local_groups_id%elements(n))

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups_id%elements(n) = j
          end if
       end do

       !
       ! now determine symmetry equivalent atoms
       !
       call group_coset_decomp(n_equal,local_groups_id,&
            cosets_id,point_trafos_buf(n_size)%matrix)

       allocate(buffer_id(n_size)%position(3,n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_id(n_size)%idipole(3,n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_id(n_size)%idipole1(3,n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_id(n_size)%pol_tensor(3,3,n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_id(n_size)%group(n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)

       n_equal_chck=0
       j_n_equal: do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_id%elements(1,j)),position)
          position3(1) = position2(1)
          position3(2) = position2(3)
          position3(3) = position2(2)
          k_N_tot: do k=1,N_total
             if(.not.help_dim(k)) then
                if (sqrt(dot_product(position3-P(:,k),position3-P(:,k))) <= small) then
                   help_dim(k)=.true.
                   n_equal_chck=n_equal_chck+1
                   buffer_id(n_size)%position(:,j)=P(:,k)
                   buffer_id(n_size)%group(j)=ef_grp(k)
                   buffer_id(n_size)%idipole(:,j)=0.0_r8_kind
                   buffer_id(n_size)%idipole1(:,j)=0.0_r8_kind
                   buffer_id(n_size)%pol_tensor(:,:,j)=Alpha(:,:,k)
                   buffer_id(n_size)%n_equal_dipoles=n_equal
                   buffer_id(n_size)%name=trim(name(k))
                endif
             endif
          end do k_N_tot
       end do j_n_equal
       ASSERT(n_equal_chck==n_equal)

       deallocate(local_groups_id%elements, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do i_N_tot

    do_alloc=.false.
    N_ipd=n_size
    if(.not. associated(ipd_array)) then
       allocate(ipd_array(N_ipd),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do_alloc=.true.
    end if
    do i=1,N_ipd
       ipd_array(i)%name=buffer_id(i)%name
       ipd_array(i)%N_equal_dipoles=buffer_id(i)%n_equal_dipoles
       N_eq=ipd_array(i)%N_equal_dipoles
       if(do_alloc) then
          allocate(ipd_array(i)%position(3,N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(ipd_array(i)%group(N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(ipd_array(i)%idipole(3,N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(ipd_array(i)%idipole1(3,N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(ipd_array(i)%pol_tensor(3,3,N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end if
       ipd_array(i)%position=buffer_id(i)%position
       ipd_array(i)%group=buffer_id(i)%group
       ipd_array(i)%idipole=buffer_id(i)%idipole
       ipd_array(i)%idipole1=buffer_id(i)%idipole1
       ipd_array(i)%pol_tensor=buffer_id(i)%pol_tensor

       deallocate(buffer_id(i)%position,stat=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(buffer_id(i)%idipole,stat=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(buffer_id(i)%pol_tensor,stat=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(buffer_id(i)%group,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do
    deallocate(buffer_id,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    deallocate(help_dim,stat=alloc_stat)
    ASSERT(alloc_stat==0)

!    if(operations_gradients .and. moving_Pol_centers) then
       allocate(point_trafos_id(n_size),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          n_equal=ipd_array(i)%n_equal_dipoles
          allocate(point_trafos_id(i)%matrix(n_equal,n_equal,group_num_el),stat=alloc_stat)
          ASSERT(alloc_stat==0)

          point_trafos_id(i)%matrix=point_trafos_buf(i)%matrix

          deallocate(point_trafos_buf(i)%matrix,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
!    end if
    deallocate(point_trafos_buf, stat=alloc_stat)
    ASSERT(alloc_stat==0)

!    if(operations_gradients .and. moving_Pol_centers) then
       allocate(unique_Pol_points(n_size), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          unique_Pol_points(i)%N_equal_atoms=ipd_array(i)%N_equal_dipoles
          allocate(unique_Pol_points(i)%position(3,ipd_array(i)%N_equal_dipoles), &
               stat=alloc_stat)
          ASSERT(alloc_stat==0)

          unique_Pol_points(i)%position=ipd_array(i)%position
          unique_Pol_points(i)%lmax_all=1
       enddo

       call symm_adapt_centers(unique_Pol_points, point_trafos_id, L_MAX=1)

       do i=1,n_size
          deallocate(point_trafos_id(i)%matrix, stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
       deallocate(point_trafos_id, stat=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(point_trafos_id)

       call pol_points_grad_inform()
!    end if

  end subroutine symm_Pol_points
  !*************************************************************

  !*************************************************************
  subroutine pol_points_grad_inform()
    use unique_atom_module, only: unique_atom_type, unique_atom_partner_type, &
         unique_atom_symadapt_type
    use symmetry_data_module, only: get_totalsymmetric_irrep,find_totalsymmetric_irrep
    use uatom_symmadapt, only: sa_free
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    type(unique_atom_type)         , pointer :: ua
    type(unique_atom_partner_type) , pointer :: sap
    type(unique_atom_symadapt_type), pointer :: sa
    type(arrmat3)                  , pointer :: gi
    integer(i4_kind) :: i,j,k
    integer(i4_kind) :: i_ua,i_if,i_cd,i_ea,ts, counter, &
         n_indep,n_equals,alloc_stat
    logical :: do_alloc
    !----------- executable code -----------------------------

    if(N_ipd > 0) then
       do_alloc=.false.
       if(.not.associated(ind_dip_grad_info)) then
          allocate(ind_dip_grad_info(N_ipd),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if

       call find_totalsymmetric_irrep()
       ts = get_totalsymmetric_irrep()

       do i_ua=1,N_ipd
          ua => unique_Pol_points(i_ua)
          sap => ua%symadapt_partner(ts,1) ! l = 1
          gi => ind_dip_grad_info(i_ua)
          n_indep = sap%N_independent_fcts
          n_equals = ua%N_equal_atoms
          if(do_alloc) then
             allocate ( gi%m(n_indep,3,n_equals),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          gi%m = 0.0_r8_kind

          do i_if=1,n_indep ! loop over all sym. adapted gradients of ua
             sa => sap%symadapt(i_if,1) ! first partner
             do i_cd=1,sa%N_fcts ! loop over all contributing derivatives
                i_ea = sa%I_equal_atom(i_cd)

                select case (sa%m(i_cd))
                case (2); gi%m(i_if,1,i_ea) = sa%c(i_cd)
                case (3); gi%m(i_if,2,i_ea) = sa%c(i_cd)
                case (1); gi%m(i_if,3,i_ea) = sa%c(i_cd)
                case default
                   call error_handler("external_point_grad_information: sth. wrong")
                end select

             enddo
          enddo
       end do

       if(do_alloc) then
          allocate(ind_dip_grad_index(N_ipd+1),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end if

       counter = 1
       do i_ua=1,N_ipd
          ua => unique_Pol_points(i_ua)
          ind_dip_grad_index(i_ua)=counter
          counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
       enddo
       ind_dip_grad_index(N_ipd+1)=counter
       totsym_grad_idip_length=counter-1

       do i=1,size(unique_Pol_points)
          ua=>unique_Pol_points(i)
          deallocate(ua%position,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          do j=1,size(ua%symadapt_partner,1)
             do k=1,ua%lmax_all
                call sa_free(ua%symadapt_partner(j,k))
             end do
          end do
          deallocate(ua%symadapt_partner,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(unique_Pol_points,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       if(.not. allocated(totsym_ind_dip)) then
          allocate(totsym_ind_dip(totsym_grad_idip_length),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end if
       totsym_ind_dip=0.0_r8_kind
    end if

  end subroutine pol_points_grad_inform
  !*************************************************************

  !*************************************************************
  subroutine dealloc_id
    !  Purpose: deallocate allocated external centers
    !------------ Modules used ------------------- ---------------
    use comm_module, only: comm_i_am_master
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,status
    !------------ Executable code --------------------------------
    if(associated(ipd_array)) then

       do i=1,N_ipd
          deallocate(ipd_array(i)%position,stat=status)
          ASSERT(status==0)
          deallocate(ipd_array(i)%idipole,stat=status)
          ASSERT(status==0)
          deallocate(ipd_array(i)%idipole1,stat=status)
          ASSERT(status==0)
          deallocate(ipd_array(i)%pol_tensor,stat=status)
          ASSERT(status==0)
          if(comm_i_am_master()) then
             deallocate(ipd_array(i)%group,stat=status)
             ASSERT(status==0)
          end if
       end do

       deallocate(ipd_array,stat=status)
       ASSERT(status==0)

       deallocate(totsym_ind_dip,STAT=status)
       ASSERT(status==0)
    end if

  end subroutine dealloc_id
  !*************************************************************

  !*************************************************************
  subroutine dealloc_Pol_center_inform()
    ! Purpose: deallocate the variables 'ind_dip_grad_info'
    !** End of interface *****************************************
    ! ----------- declaration of local variables -------------
    integer(kind=i4_kind) :: i,alloc_stat
    ! ----------- executable code -----------------------------

    if(associated(ind_dip_grad_info)) then
       do i=1,size(ind_dip_grad_info)
          deallocate(ind_dip_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(ind_dip_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(ind_dip_grad_info)

       deallocate(ind_dip_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

  end subroutine dealloc_Pol_center_inform
  !*************************************************************

  !*************************************************************
  subroutine output_geometry_id()
    ! Purpose: write the geomtry of external polarizable centers to output
    use iounitadmin_module, only : output_unit
    implicit none
    ! *** end of interface ***

    integer(i4_kind)            :: i_eq,i_ua
    real(r8_kind),parameter     :: d2au=2.541766_r8_kind

    if(N_ipd > 0) then
       write(output_unit,*)"-- Geometry of Polarizable Centers in au --------------------------------"
       write(output_unit,*)"          pol. tensor                                      coordinates"

       do i_ua=1,N_ipd
          do i_eq=1,ipd_array(i_ua)%N_equal_dipoles
             write(output_unit,1000) ipd_array(i_ua)%pol_tensor(1,1:3,i_eq),&
                  trim(ipd_array(i_ua)%name),&
                  ipd_array(i_ua)%position(1,i_eq),&
                  ipd_array(i_ua)%position(2,i_eq),&
                  ipd_array(i_ua)%position(3,i_eq)
             write(output_unit,1200) ipd_array(i_ua)%pol_tensor(2,1:3,i_eq)
             write(output_unit,1200) ipd_array(i_ua)%pol_tensor(3,1:3,i_eq)
          enddo
       enddo
    end if

1000 format('',3(F9.4,2X),2X,A4,2X,3(F11.6,7X))
1200 format('',3(F9.4,2X))
  end subroutine output_geometry_id
  !*************************************************************

  !*****************************************************************************
  subroutine pol_centers_bcast()
    !  Purpose:  broadcasts all information about external polarizable points
    !  called by send_recv_init_options
    use comm,               only: comm_bcast                                   &
                                , comm_rank
    !** End of interface *******************************************************
    !----------- declaration of local variables --------------------------------
    type(ind_dipole_type), pointer :: ipd
    integer(i4_kind)               :: i, status
    logical                        :: do_alloc
    !--- executable code--------------------------------------------------------
    !
    call comm_bcast( moving_Pol_centers )
#ifdef WITH_EFP
    call comm_bcast( efp                )
#endif
    call comm_bcast( N_ipd              )
    call comm_bcast( do_pol_pcm         )
    !
    if( N_ipd > 0 ) then
      !
      do_alloc=.false.
      if( comm_rank() /= 0 .and. .not. associated(ipd_array) ) then
        allocate(ipd_array(N_ipd), stat=status)
        ASSERT(status==0)
        do_alloc=.true.
      end if
      !
      do i = 1, N_ipd
        !
        ipd => ipd_array(i)
        !
        call comm_bcast( ipd%name            )
        call comm_bcast( ipd%N_equal_dipoles )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate(ipd%position(3,ipd%N_equal_dipoles),     stat=status)
          ASSERT(status==0)
          allocate(ipd%idipole(3,ipd%N_equal_dipoles),      stat=status)
          ASSERT(status==0)
          allocate(ipd%idipole1(3,ipd%N_equal_dipoles),     stat=status)
          ASSERT(status==0)
          allocate(ipd%pol_tensor(3,3,ipd%N_equal_dipoles), stat=status)
          ASSERT(status==0)
        end if
        call comm_bcast( ipd%position   )
        call comm_bcast( ipd%idipole    )
        call comm_bcast( ipd%idipole1   )
        call comm_bcast( ipd%pol_tensor )
        !
      end do
      !
      call pol_points_info_bcast()
      !
      if( comm_rank() /= 0 .and. .not.allocated(totsym_ind_dip) ) then
        allocate( totsym_ind_dip(totsym_grad_idip_length), stat=status)
        ASSERT(status==0)
      end if
      call comm_bcast( totsym_ind_dip )
      !
    end if
    !
    contains
    !
    subroutine pol_points_info_bcast()
      !------------ Declaration of local variables -----------------------------
      integer(i4_kind) :: j, status
      integer(i4_kind) :: n1, n2, n3
      logical          :: do_alloc
      !------------ Executable code --------------------------------------------
      !
      do_alloc=.false.
      !
      if( comm_rank() /= 0 .and. .not. associated(ind_dip_grad_info) ) then
        !
        allocate(ind_dip_grad_info(N_ipd),stat=status)
        ASSERT(status==0)
        do_alloc=.true.
        !
      end if
      !
      do j = 1, N_ipd
        if ( comm_rank() == 0 ) then
          ! calculate dimensions
          n1=size(ind_dip_grad_info(j)%m,1)
          n2=size(ind_dip_grad_info(j)%m,2)
          n3=size(ind_dip_grad_info(j)%m,3)
        endif
        !
        call comm_bcast( n1 )
        call comm_bcast( n2 )
        call comm_bcast( n3 )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( ind_dip_grad_info(j)%m(n1,n2,n3), stat=status )
          ASSERT(status==0)
        end if
        !
        call comm_bcast( ind_dip_grad_info(j)%m )
        !
      end do
      !
      if( comm_rank() /= 0 .and. do_alloc ) then
        allocate( ind_dip_grad_index(N_ipd+1), stat=status )
        ASSERT(status==0)
      end if
      !
      call comm_bcast( ind_dip_grad_index      )
      call comm_bcast( totsym_grad_idip_length )
      !
    end subroutine pol_points_info_bcast
    !
  end subroutine pol_centers_bcast
  !*****************************************************************************

  !*************************************************************
  subroutine init_Pol_centers_grads()
    ! Purpose: allocate arrays to calculate gradients on pol. centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code --------------------------------

    if(N_ipd > 0) then
       allocate(grad_idip_cartesian(N_ipd),stat=status)
       ASSERT(status==0)
#ifdef WITH_EFP
       if(efp) then
          allocate(torque_idip_cartesian(N_ipd),stat=status)
          ASSERT(status==0)
       end if
#endif

       do i=1,N_ipd
          n_eq=ipd_array(i)%N_equal_dipoles
          allocate(grad_idip_cartesian(i)%m(3,n_eq),stat=status)
          ASSERT(status==0)
          grad_idip_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
          if(efp) then
             allocate(torque_idip_cartesian(i)%m(3,n_eq),stat=status)
             ASSERT(status==0)
             torque_idip_cartesian(i)%m=0.0_r8_kind
          end if
#endif
       end do

       allocate(grad_idip_totalsym(totsym_grad_idip_length),stat=status)
       ASSERT(status==0)
       grad_idip_totalsym=0.0_r8_kind
#ifdef WITH_EFP
       if(efp) then
          allocate(torque_idip_totalsym(totsym_grad_idip_length),stat=status)
          ASSERT(status==0)
          torque_idip_totalsym=0.0_r8_kind
       end if
#endif
    end if

  end subroutine init_Pol_centers_grads
  !*************************************************************

  !*************************************************************
  subroutine Pol_centers_grads_shutdown()
    !  Purpose: free arrays for calculating gradients on pol. centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i
    !------------ Executable code --------------------------------

    if(N_ipd > 0) then
       do i=1,N_ipd
          deallocate(grad_idip_cartesian(i)%m,stat=status)
          ASSERT(status==0)
#ifdef WITH_EFP
          if(efp) then
             deallocate(torque_idip_cartesian(i)%m,stat=status)
             ASSERT(status==0)
          end if
#endif
       end do
       deallocate(grad_idip_cartesian,stat=status)
       ASSERT(status==0)
#ifdef WITH_EFP
       if(efp) then
          deallocate(torque_idip_cartesian,stat=status)
          ASSERT(status==0)
       end if
#endif

       deallocate(grad_idip_totalsym,stat=status)
       ASSERT(status==0)
#ifdef WITH_EFP
       if(efp) then
          deallocate(torque_idip_totalsym,stat=status)
          ASSERT(status==0)
       end if
#endif
    end if

  end subroutine Pol_centers_grads_shutdown
  !*************************************************************

  !*************************************************************
  subroutine totsym_id_grad_unpack()
    !  Purpose: unpack 2 center contribution to the Pol. center gradient from the
    !           slave
    use comm_module, only: communpack
    !------------ Declaration of local variables ----------------
    real(r8_kind) :: help_arr_ipd(totsym_grad_idip_length)
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    if(N_ipd > 0) then
       call communpack(help_arr_ipd,totsym_grad_idip_length,1,info)
       ASSERT(info==0)
       grad_idip_totalsym=grad_idip_totalsym+help_arr_ipd
#ifdef WITH_EFP
       if(efp) then
          call communpack(help_arr_ipd,totsym_grad_idip_length,1,info)
          ASSERT(info==0)
          torque_idip_totalsym=torque_idip_totalsym+help_arr_ipd
       end if
#endif
    end if

  end subroutine totsym_id_grad_unpack
  !*************************************************************

  !*************************************************************
  subroutine totsym_id_grad_pack()
    !  Purpose: pack 2 center contribution to the Pol center gradient on the
    !           slave
    use comm_module, only: commpack
    !------------ Declaration of local variables ----------------
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    if(N_ipd > 0) then
       call commpack(grad_idip_totalsym,totsym_grad_idip_length,1,info)
       ASSERT(info==0)
#ifdef WITH_EFP
       if(efp) then
          call commpack(torque_idip_totalsym,totsym_grad_idip_length,1,info)
          ASSERT(info==0)
       end if
#endif
    end if

  end subroutine totsym_id_grad_pack
  !*************************************************************

  !*************************************************************
  subroutine send_receive_id()
    ! Purpose: Master distributes calculated induced dipoles
    ! between slaves
    ! called by calc_induced_dipmom and main_slave
    use comm_module, only: comm_i_am_master, comm_init_send, comm_send, &
         comm_all_other_hosts, commpack, communpack
    use msgtag_module, only : msgtag_ind_dipmom
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i,n_eq,info
    real(r8_kind), pointer :: id(:,:)
    !--- executable code-------------------------------------

    if(comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_ind_dipmom)
       do i=1,N_ipd
          n_eq=ipd_array(i)%N_equal_dipoles
          id=>ipd_array(i)%idipole
          call commpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
          id=>ipd_array(i)%idipole1
          call commpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
       end do
       call commpack(totsym_ind_dip(1),totsym_grad_idip_length,1,info)
       ASSERT(info==0)
       call comm_send()
    else
       do i=1,N_ipd
          n_eq=ipd_array(i)%N_equal_dipoles
          id=>ipd_array(i)%idipole
          call communpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
          id=>ipd_array(i)%idipole1
          call communpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
       end do
       call communpack(totsym_ind_dip(1),totsym_grad_idip_length,1,info)
       ASSERT(info==0)
    end if

  end subroutine send_receive_id
  !*************************************************************

  !*************************************************************
  subroutine send_receive_id1()
    ! Purpose: Master distributes calculated induced dipoles
    ! between slaves
    ! called by calc_induced_dipmom and main_slave
    use comm_module, only : comm_i_am_master, comm_init_send, comm_send, &
         comm_all_other_hosts,comm_save_recv,comm_master_host, commpack, communpack
    use msgtag_module, only : msgtag_ind_dipmom
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i,n_eq,info
    real(r8_kind), pointer :: id(:,:)
    !--- executable code-------------------------------------

    if(comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_ind_dipmom)
       do i=1,N_ipd
          n_eq=ipd_array(i)%N_equal_dipoles
          id=>ipd_array(i)%idipole
          call commpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
          id=>ipd_array(i)%idipole1
          call commpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
       end do
!!$       call commpack(totsym_ind_dip(1),totsym_grad_idip_length,1,info)
!!$       ASSERT(info==0)
       call comm_send()
    else
       call comm_save_recv(comm_master_host,msgtag_ind_dipmom)
       do i=1,N_ipd
          n_eq=ipd_array(i)%N_equal_dipoles
          id=>ipd_array(i)%idipole
          call communpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
          id=>ipd_array(i)%idipole1
          call communpack(id(1,1),3*n_eq,1,info)
          ASSERT(info==0)
       end do
!!$       call communpack(totsym_ind_dip(1),totsym_grad_idip_length,1,info)
!!$       ASSERT(info==0)
    end if

  end subroutine send_receive_id1
  !*************************************************************

  !*************************************************************
  subroutine free_field_arrays()
    !
    ! This looks  like a clean-up  function. It should be  executed by
    ! all  workers. It  is called  from main_scf()  that will  be soon
    ! running on all workers. Check this next time it is used.
    !
    use operations_module, only: operations_integral
    use elec_static_field_module, only: deallocate_field
    use elec_static_field_module, only: bounds_free_field,destroy_field_file
    use elec_static_field_module, only: surf_points_gradinfo_dealloc,dealloc_surf_points
    use options_module, only: options_integrals_on_file
    use integralstore_module, only: integralstore_deallocate_pcm
    implicit none
    !** End of interface *****************************************

    ABORT("check!")
    if (.not. options_integrals_on_file() .and. operations_integral) then
       call integralstore_deallocate_pcm()
    endif

    call surf_points_gradinfo_dealloc()

    if (operations_integral) then
       call bounds_free_field()
       call destroy_field_file()
    end if
    call dealloc_surf_points()
    call deallocate_field()
  end subroutine free_field_arrays
  !*************************************************************
  !--------------- End of module ----------------------------------
end module induced_dipoles_module

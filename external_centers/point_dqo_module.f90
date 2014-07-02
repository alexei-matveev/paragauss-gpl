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
module point_dqo_module
  !-------------------------------------------------------------------
  !
  !  Purpose: keep information on eXternal point dipoles,
  !           point quadrupoles, point octopoles and repulsive
  !           centers
  !
  !  Module called by: integral and gradient parts
  !
  !
  !  References: ...
  !
  !
  !  Author: AS
  !  Date: 31.05.2007
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
#include <def.h>

  use type_module ! type specification parameters
  use datatype
  use unique_atom_module, only: unique_atom_type
  use pointcharge_module
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of types ---------------------------------
  type, public :: pointdipole_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_dipoles
     real(r8_kind) :: position_1(3)
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind) :: dipole_1(3)
     real(r8_kind), pointer :: dipole(:,:)
     integer(i4_kind), pointer :: group(:)
  end type pointdipole_type

  type, public :: pointquadrupole_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_qpoles
     real(r8_kind) :: position_1(3)
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind) :: quadrupole_1(3,3)
     real(r8_kind), pointer :: quadrupole(:,:,:)
     integer(i4_kind), pointer :: group(:)
  end type pointquadrupole_type

  type, public :: pointoctopole_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_opoles
     real(r8_kind) :: position_1(3)
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind) :: octopole_1(3,3,3)
     real(r8_kind), pointer :: octopole(:,:,:,:)
     integer(i4_kind), pointer :: group(:)
  end type pointoctopole_type

  type, public :: repulsion_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_centers
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind) :: C
     real(r8_kind) :: a
     integer(i4_kind), pointer :: group(:)
  end type repulsion_type

  !------------ Declaration of constants and variables ---------------
  integer(i4_kind), public :: N_pd=0_i4_kind
  integer(i4_kind), public :: N_pq=0_i4_kind
  integer(i4_kind), public :: N_po=0_i4_kind
  integer(i4_kind), public :: N_po1=0_i4_kind
  integer(i4_kind), public :: N_rc=0_i4_kind
  type(pointdipole_type), pointer, public :: pd_array(:)
  type(pointquadrupole_type), pointer, public :: pq_array(:)
  type(pointoctopole_type), pointer, public :: po_array(:)
  type(repulsion_type), pointer, public :: rc_array(:)

  type(arrmat3), pointer, public :: dipoles_grad_info(:)
  integer(i4_kind), allocatable, public :: dipoles_grad_index(:)
  type(arrmat3), pointer, public :: qpoles_grad_info(:)
  integer(i4_kind), allocatable, public :: qpoles_grad_index(:)
  type(arrmat3), pointer, public :: opoles_grad_info(:)
  integer(i4_kind), allocatable, public :: opoles_grad_index(:)
  type(arrmat3), pointer, public :: rep_grad_info(:)
  integer(i4_kind), allocatable, public :: rep_grad_index(:)

  logical, public :: moving_X_centers=.false.
  logical, public :: moving_R_centers=.false.
  integer(i4_kind), public :: totsym_grad_dip_length=0
  integer(i4_kind), public :: totsym_grad_quad_length=0
  integer(i4_kind), public :: totsym_grad_oct_length=0
  integer(i4_kind), public :: totsym_grad_rep_length=0

  real(r8_kind),public,allocatable,target :: gradient_dip_totalsym(:)
  real(r8_kind),public,allocatable,target :: torque_dip_totalsym(:)
  real(r8_kind),public,allocatable,target :: gradient_quad_totalsym(:)
  real(r8_kind),public,allocatable,target :: torque_quad_totalsym(:)
  real(r8_kind),public,allocatable,target :: gradient_oct_totalsym(:)
  real(r8_kind),public,allocatable,target :: torque_oct_totalsym(:)
  real(r8_kind),public,allocatable,target :: gradient_rep_totalsym(:)

  type(arrmat2),public,allocatable :: gradient_dip_cartesian(:)
  type(arrmat2),public,allocatable :: torque_dip_cartesian(:)
  type(arrmat2),public,allocatable :: gradient_quad_cartesian(:)
  type(arrmat2),public,allocatable :: torque_quad_cartesian(:)
  type(arrmat2),public,allocatable :: gradient_oct_cartesian(:)
  type(arrmat2),public,allocatable :: torque_oct_cartesian(:)
  type(arrmat2),public,allocatable :: gradient_rep_cartesian(:)
  type(arrmat2),public,allocatable :: torque_rep_cartesian(:)

  logical, public :: print_X_grad=.false.
  logical, public :: print_R_grad=.false.

  !------------ Interface statements -----------------------------

  !------------ public functions and subroutines ---------------------
  public :: dealloc_dqo
  public :: dealloc_center_inform
  public :: symm_external_points
  public :: output_geometry_ec
  public :: external_centers_bcast
  public :: calc_nuc_dqo_energy
  public :: calc_nuc_dqo_grads
  public :: init_X_centers_grads
  public :: X_centers_grads_shutdown
  public :: totsym_X_grad_unpack
  public :: totsym_X_grad_pack
  public :: transform_X_grad_to_cart
  public :: X_grad_cart_write
  public :: calc_X_grads

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of types ---------------------------------
  !------------ Declaration of constants and variables ---------------
  type(unique_atom_type), pointer :: unique_external_points(:)
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  subroutine symm_external_points (N_total, P, ig, ef_grp, name, &
       Z, Cz, Az, Cz_f, Az_f, D, Q, O, C, A)
    !
    ! Symmetrize    external   point   centers    (charges,   dipoles,
    ! quadrupoles, octopoles, repulsive centers)
    !
    use operations_module, only: operations_gradients
    use group_module, only: group_num_el, symm_transformation_int, &
         sub_group, group_coset, ylm_trafos, group_coset_decomp
    use symm_module, only: symm_adapt_centers
    implicit none
    integer(i4_kind) :: N_total
    real(r8_kind) :: P(3,N_total)
    integer(i4_kind) :: ig(N_total)
    integer(i4_kind) :: ef_grp(N_total)
    character(len=12)  :: name(N_total)
    real(r8_kind), optional :: Z(N_total)
    real(r8_kind), optional :: Cz(N_total)
    real(r8_kind), optional :: Az(N_total)
    real(r8_kind), optional :: Cz_f(N_total)
    real(r8_kind), optional :: Az_f(N_total)
    real(r8_kind), optional :: D(3,N_total)
    real(r8_kind), optional :: Q(3,3,N_total)
    real(r8_kind), optional :: O(3,3,3,N_total)
    real(r8_kind), optional :: C(N_total)
    real(r8_kind), optional :: A(N_total)
    !** End of interface *****************************************

    real(r8_kind),parameter :: small = 1.e-10_r8_kind
    ! very small value

    logical :: do_Z,do_D, do_Q, do_O, do_R
    integer (i4_kind) :: i, j, k, n, l, alloc_stat, n_size
    ! counters
    integer(i4_kind) :: n_equal,n_equal_check,N_eq,i_group
    ! number of equal atoms
    real(r8_kind) :: position(3),position2(3),position3(3)
    ! coordinates of unique atom in the order (x,z,y)
    type(symm_transformation_int),allocatable :: point_trafos_buf(:)
    type(symm_transformation_int), pointer :: point_trafos_ec(:)
    type(sub_group) :: local_groups_ec
    ! local symmetry group of the unique external centers
    type(group_coset) :: cosets_ec
    ! coset of the local symmetry group of the unique external centers
    logical, allocatable :: help_dim(:)
    type(pointcharge_type), allocatable :: buffer_z(:)
    type(pointdipole_type), allocatable :: buffer_d(:)
    type(pointquadrupole_type), allocatable :: buffer_q(:)
    type(pointoctopole_type), allocatable :: buffer_o(:)
    type(repulsion_type), allocatable :: buffer_r(:)
    logical :: do_alloc
    !------------ Executable code ------------------------------------

    n=0
    do_Z=.false.; do_D=.false.; do_Q=.false.; do_O=.false.; do_R=.false.

    if(present(Z) .and. present(Cz) .and. present(Az) &
         .and. present(Cz_f) .and. present(Az_f)) then
       do_Z=.true.; n=n+1
    end if
    if(present(D)) then
       do_D=.true.; n=n+1
    end if
    if(present(Q)) then
       do_Q=.true.; n=n+1
    end if
    if(present(O)) then
       do_O=.true.; n=n+1
    end if
    if(present(C) .and. present(A)) then
       do_R=.true.; n=n+1
    end if
    if(n /= 1) &
         call error_handler("Point_dqo_module:symm_external_points: just one type of external points can be treated")

    allocate(help_dim(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    help_dim=.false.

    allocate(point_trafos_buf(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    if(do_D) then
       allocate(buffer_d(N_total),stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_Q) then
       allocate(buffer_q(N_total),stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_O) then
       allocate(buffer_o(N_total),stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_R) then
       allocate(buffer_r(N_total),stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_Z) then
       allocate(buffer_z(N_total),stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

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

       i_group=ig(i)

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
       local_groups_ec%num_el = n
       allocate(local_groups_ec%elements(n))

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups_ec%elements(n) = j
          end if
       end do

       !
       ! now determine symmetry equivalent atoms
       !
       call group_coset_decomp(n_equal,local_groups_ec,&
            cosets_ec,point_trafos_buf(n_size)%matrix)

       if(do_D) then
          allocate(buffer_d(n_size)%position(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_d(n_size)%dipole(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_d(n_size)%group(n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       else if(do_Q) then
          allocate(buffer_q(n_size)%position(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_q(n_size)%quadrupole(3,3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_q(n_size)%group(n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       else if(do_O) then
          allocate(buffer_o(n_size)%position(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_o(n_size)%octopole(3,3,3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_o(n_size)%group(n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       else if(do_R) then
          allocate(buffer_r(n_size)%position(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_r(n_size)%group(n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       else if(do_Z) then
          allocate(buffer_z(n_size)%position(3,n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(buffer_z(n_size)%group(n_equal),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end if

       n_equal_check=0
       j_n_equal: do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_ec%elements(1,j)),position)
          position3(1) = position2(1)
          position3(2) = position2(3)
          position3(3) = position2(2)
          k_N_tot: do k=1,N_total
             if(.not.help_dim(k)) then
                if(i_group /= ig(k)) cycle k_N_tot
                if (sqrt(dot_product(position3-P(:,k),position3-P(:,k))) <= small) then
                   help_dim(k)=.true.
                   n_equal_check=n_equal_check+1
                   if(do_D) then
                      buffer_d(n_size)%position(:,j)=P(:,k)
                      buffer_d(n_size)%dipole(:,j)=D(:,k)
                      buffer_d(n_size)%group(j)=ef_grp(k)
                      buffer_d(n_size)%n_equal_dipoles=n_equal
                      buffer_d(n_size)%name=trim(name(k))
                   else if(do_Q) then
                      buffer_q(n_size)%position(:,j)=P(:,k)
                      buffer_q(n_size)%quadrupole(:,:,j)=Q(:,:,k)
                      buffer_q(n_size)%group(j)=ef_grp(k)
                      buffer_q(n_size)%n_equal_qpoles=n_equal
                      buffer_q(n_size)%name=trim(name(k))
                   else if(do_O) then
                      buffer_o(n_size)%position(:,j)=P(:,k)
                      buffer_o(n_size)%octopole(:,:,:,j)=O(:,:,:,k)
                      buffer_o(n_size)%group(j)=ef_grp(k)
                      buffer_o(n_size)%n_equal_opoles=n_equal
                      buffer_o(n_size)%name=trim(name(k))
                   else if(do_R) then
                      buffer_r(n_size)%position(:,j)=P(:,k)
                      buffer_r(n_size)%group(j)=ef_grp(k)
                      buffer_r(n_size)%n_equal_centers=n_equal
                      buffer_r(n_size)%name=trim(name(k))
                      buffer_r(n_size)%C=C(k)
                      buffer_r(n_size)%A=A(k)
                   else if(do_Z) then
                      buffer_z(n_size)%position(:,j)=P(:,k)
                      buffer_z(n_size)%Z=Z(k)
                      buffer_z(n_size)%group(j)=ef_grp(k)
                      buffer_z(n_size)%n_equal_charges=n_equal
                      buffer_z(n_size)%name=trim(name(k))
                      buffer_z(n_size)%C=Cz(k)
                      buffer_z(n_size)%A=Az(k)
                      buffer_z(n_size)%Cf=Cz_f(k)
                      buffer_z(n_size)%Af=Az_f(k)
                   end if
                end if
             endif
          end do k_N_tot
       end do j_n_equal
       ASSERT(n_equal_check==n_equal)

       deallocate(local_groups_ec%elements, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do i_N_tot

    if(do_D) then
       N_pd=n_size
       do_alloc=.false.
       if(.not. associated(pd_array)) then
          allocate(pd_array(N_pd),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if
       do i=1,N_pd
          pd_array(i)%name=buffer_d(i)%name
          pd_array(i)%N_equal_dipoles=buffer_d(i)%n_equal_dipoles
          N_eq=pd_array(i)%N_equal_dipoles
          if(do_alloc) then
             allocate(pd_array(i)%position(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(pd_array(i)%dipole(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(pd_array(i)%group(N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          pd_array(i)%position=buffer_d(i)%position
          pd_array(i)%dipole=buffer_d(i)%dipole
          pd_array(i)%group=buffer_d(i)%group

          deallocate(buffer_d(i)%position,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_d(i)%dipole,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_d(i)%group,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(buffer_d,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_Q) then
       N_pq=n_size
       do_alloc=.false.
       if(.not. associated(pq_array)) then
          allocate(pq_array(N_pq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if
       do i=1,N_pq
          pq_array(i)%name=buffer_q(i)%name
          pq_array(i)%N_equal_qpoles=buffer_q(i)%n_equal_qpoles
          N_eq=pq_array(i)%N_equal_qpoles
          if(do_alloc) then
             allocate(pq_array(i)%position(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(pq_array(i)%quadrupole(3,3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(pq_array(i)%group(N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          pq_array(i)%position=buffer_q(i)%position
          pq_array(i)%quadrupole=buffer_q(i)%quadrupole
          pq_array(i)%group=buffer_q(i)%group

          deallocate(buffer_q(i)%position,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_q(i)%quadrupole,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_q(i)%group,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(buffer_q,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_O) then
       N_po=n_size
       N_po1=0
       do_alloc=.false.
       if(.not. associated(po_array)) then
          allocate(po_array(N_po),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if
       do i=1,N_po
          po_array(i)%name=buffer_o(i)%name
          po_array(i)%N_equal_opoles=buffer_o(i)%n_equal_opoles
          N_eq=po_array(i)%N_equal_opoles
          if(do_alloc) then
             allocate(po_array(i)%position(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(po_array(i)%octopole(3,3,3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(po_array(i)%group(N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          po_array(i)%position=buffer_o(i)%position
          po_array(i)%octopole=buffer_o(i)%octopole
          po_array(i)%group=buffer_o(i)%group

          deallocate(buffer_o(i)%position,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_o(i)%octopole,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_o(i)%group,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(buffer_o,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_R) then
       N_rc=n_size
       do_alloc=.false.
       if(.not. associated(rc_array)) then
          allocate(rc_array(N_rc),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if
       do i=1,N_rc
          rc_array(i)%name=buffer_r(i)%name
          rc_array(i)%N_equal_centers=buffer_r(i)%n_equal_centers
          N_eq=rc_array(i)%N_equal_centers
          if(do_alloc) then
             allocate(rc_array(i)%position(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(rc_array(i)%group(N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          rc_array(i)%position=buffer_r(i)%position
          rc_array(i)%C=buffer_r(i)%C
          rc_array(i)%A=buffer_r(i)%A
          rc_array(i)%group=buffer_r(i)%group

          deallocate(buffer_r(i)%position,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_r(i)%group,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(buffer_r,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    else if(do_Z) then
       pointcharge_N=n_size
       do_alloc=.false.
       if(.not. associated(pointcharge_array)) then
          allocate(pointcharge_array(pointcharge_N),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if
       do i=1,pointcharge_N
          pointcharge_array(i)%name=buffer_z(i)%name
          pointcharge_array(i)%N_equal_charges=buffer_z(i)%n_equal_charges
          N_eq=pointcharge_array(i)%N_equal_charges
          if(do_alloc) then
             allocate(pointcharge_array(i)%position(3,N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             allocate(pointcharge_array(i)%group(N_eq),stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          pointcharge_array(i)%position=buffer_z(i)%position
          pointcharge_array(i)%group=buffer_z(i)%group
          pointcharge_array(i)%Z=buffer_z(i)%Z
          pointcharge_array(i)%C=buffer_z(i)%C
          pointcharge_array(i)%A=buffer_z(i)%A
          pointcharge_array(i)%Cf=buffer_z(i)%Cf
          pointcharge_array(i)%Af=buffer_z(i)%Af

          deallocate(buffer_z(i)%position,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(buffer_z(i)%group,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(buffer_z,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    deallocate(help_dim,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    if(operations_gradients .and. (moving_X_centers .or. moving_R_centers .or. moving_pc)) then
       if((do_D .or. do_Q .or. do_O) .and. moving_X_centers ) allocate(point_trafos_ec(n_size),stat=alloc_stat)
       if(do_R .and. moving_R_centers ) allocate(point_trafos_ec(n_size),stat=alloc_stat)
       if(do_Z .and. moving_pc ) allocate(point_trafos_ec(n_size),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          if (do_D .and. moving_X_centers) n_equal=pd_array(i)%n_equal_dipoles
          if (do_Q .and. moving_X_centers) n_equal=pq_array(i)%n_equal_qpoles
          if (do_O .and. moving_X_centers) n_equal=po_array(i)%n_equal_opoles
          if (do_R .and. moving_R_centers) n_equal=rc_array(i)%n_equal_centers
          if (do_Z .and. moving_pc)        n_equal=pointcharge_array(i)%n_equal_charges

          allocate(point_trafos_ec(i)%matrix(n_equal,n_equal,group_num_el),stat=alloc_stat)
          ASSERT(alloc_stat==0)

          point_trafos_ec(i)%matrix=point_trafos_buf(i)%matrix

          deallocate(point_trafos_buf(i)%matrix,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
    end if
    deallocate(point_trafos_buf, stat=alloc_stat)
    ASSERT(alloc_stat==0)

    if(operations_gradients .and. (moving_X_centers .or. moving_R_centers .or. moving_pc)) then
       if((do_D .or. do_Q .or. do_O) .and. moving_X_centers ) allocate(unique_external_points(n_size), stat=alloc_stat)
       if(do_R .and. moving_R_centers ) allocate(unique_external_points(n_size), stat=alloc_stat)
       if(do_Z .and. moving_pc ) allocate(unique_external_points(n_size), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          if(do_D .and. moving_X_centers) then
             unique_external_points(i)%N_equal_atoms=pd_array(i)%N_equal_dipoles
             allocate(unique_external_points(i)%position(3,pd_array(i)%N_equal_dipoles), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             unique_external_points(i)%position=pd_array(i)%position
          else if(do_Q .and. moving_X_centers) then
             unique_external_points(i)%N_equal_atoms=pq_array(i)%N_equal_qpoles
             allocate(unique_external_points(i)%position(3,pq_array(i)%N_equal_qpoles), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             unique_external_points(i)%position=pq_array(i)%position
          else if(do_O .and. moving_X_centers) then
             unique_external_points(i)%N_equal_atoms=po_array(i)%N_equal_opoles
             allocate(unique_external_points(i)%position(3,po_array(i)%N_equal_opoles), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             unique_external_points(i)%position=po_array(i)%position
          else if(do_R .and. moving_R_centers) then
             unique_external_points(i)%N_equal_atoms=rc_array(i)%N_equal_centers
             allocate(unique_external_points(i)%position(3,rc_array(i)%N_equal_centers), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             unique_external_points(i)%position=rc_array(i)%position
          else if(do_Z .and. moving_pc) then
             unique_external_points(i)%N_equal_atoms=pointcharge_array(i)%N_equal_charges
             allocate(unique_external_points(i)%position(3,pointcharge_array(i)%N_equal_charges), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             unique_external_points(i)%position=pointcharge_array(i)%position
          end if
          unique_external_points(i)%lmax_all=1
       enddo

       call symm_adapt_centers(unique_external_points, point_trafos_ec, L_MAX=1)

       do i=1,n_size
          deallocate(point_trafos_ec(i)%matrix, stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
       deallocate(point_trafos_ec, stat=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(point_trafos_ec)

       if(do_D .and. moving_X_centers)  call external_points_grad_inform("PD")
       if(do_Q .and. moving_X_centers)  call external_points_grad_inform("PQ")
       if(do_O .and. moving_X_centers)  call external_points_grad_inform("PO")
       if(do_R .and. moving_R_centers)  call external_points_grad_inform("RC")
       if(do_Z .and. moving_pc)         call external_points_grad_inform("PC")
    end if

  end subroutine symm_external_points
  !*************************************************************

  !*************************************************************
  subroutine external_points_grad_inform(type_of_center)
    use unique_atom_module, only: unique_atom_type, unique_atom_partner_type, &
         unique_atom_symadapt_type
    use symmetry_data_module, only: get_totalsymmetric_irrep,find_totalsymmetric_irrep
    use uatom_symmadapt, only: sa_free
    !------------ Declaration of formal parameters ---------------
    character(len=2) :: type_of_center
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    type(arrmat3), pointer :: points_grad_info(:)
    integer(i4_kind), allocatable :: points_grad_index(:)
    type(unique_atom_type)         , pointer :: ua
    type(unique_atom_partner_type) , pointer :: sap
    type(unique_atom_symadapt_type), pointer :: sa
    type(arrmat3)                  , pointer :: gi
    integer(i4_kind) :: N_centers,i,j,k,m1,m2,m3
    integer(i4_kind) :: i_ua,i_if,i_cd,i_ea,ts, counter, &
         n_indep,n_equals,alloc_stat
    logical :: do_D, do_Q, do_O, do_R, do_Z
    logical :: do_alloc
    !----------- executable code -----------------------------

    do_D=.false.; do_Q=.false.; do_O=.false.; do_R=.false.; do_Z=.false.

    N_centers=0
    if(type_of_center=='PD') then
       N_centers=N_pd; do_D=.true.
    else if(type_of_center=='PQ') then
       N_centers=N_pq; do_Q=.true.
    else if(type_of_center=='PO') then
       N_centers=N_po; do_O=.true.
    else if(type_of_center=='RC') then
       N_centers=N_rc; do_R=.true.
    else if(type_of_center=='PC') then
       N_centers=pointcharge_N; do_Z=.true.
    end if
    if(N_centers==0) return

    if(N_centers > 0) then
       allocate(points_grad_info(N_centers),stat=alloc_stat)
       ASSERT(alloc_stat==0)

       call find_totalsymmetric_irrep()
       ts = get_totalsymmetric_irrep()

       do i_ua=1,N_centers
          ua => unique_external_points(i_ua)
          sap => ua%symadapt_partner(ts,1) ! l = 1
          gi => points_grad_info(i_ua)
          n_indep = sap%N_independent_fcts
          n_equals = ua%N_equal_atoms
          allocate ( gi%m(n_indep,3,n_equals),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
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

       allocate(points_grad_index(N_centers+1),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       counter = 1
       do i_ua=1,N_centers
          ua => unique_external_points(i_ua)
          points_grad_index(i_ua)=counter
          counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
       enddo
       points_grad_index(N_centers+1)=counter

       do_alloc=.false.
       if (do_D) then
          if(.not. associated(dipoles_grad_info)) then
             allocate(dipoles_grad_info(N_centers),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             do_alloc=.true.
          end if
          do i=1,N_pd
             m1=size(points_grad_info(i)%m,1)
             m2=size(points_grad_info(i)%m,2)
             m3=size(points_grad_info(i)%m,3)
             if(do_alloc) then
                allocate(dipoles_grad_info(i)%m(m1,m2,m3),stat=alloc_stat)
                ASSERT(alloc_stat==0)
             end if
             dipoles_grad_info(i)%m=points_grad_info(i)%m
          end do

          if(do_alloc) then
             allocate(dipoles_grad_index(N_centers+1),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          dipoles_grad_index=points_grad_index
          totsym_grad_dip_length=counter-1
       else if (do_Q) then
          if(.not. associated(qpoles_grad_info)) then
             allocate(qpoles_grad_info(N_centers),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             do_alloc=.true.
          end if
          do i=1,N_pq
             m1=size(points_grad_info(i)%m,1)
             m2=size(points_grad_info(i)%m,2)
             m3=size(points_grad_info(i)%m,3)
             if(do_alloc) then
                allocate(qpoles_grad_info(i)%m(m1,m2,m3),stat=alloc_stat)
                ASSERT(alloc_stat==0)
             end if
             qpoles_grad_info(i)%m=points_grad_info(i)%m
          end do

          if(do_alloc) then
             allocate(qpoles_grad_index(N_centers+1),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          qpoles_grad_index=points_grad_index
          totsym_grad_quad_length=counter-1
       else if (do_O) then
          if(.not. associated(opoles_grad_info)) then
             allocate(opoles_grad_info(N_centers),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             do_alloc=.true.
          end if
          do i=1,N_po
             m1=size(points_grad_info(i)%m,1)
             m2=size(points_grad_info(i)%m,2)
             m3=size(points_grad_info(i)%m,3)
             if(do_alloc) then
                allocate(opoles_grad_info(i)%m(m1,m2,m3),stat=alloc_stat)
                ASSERT(alloc_stat==0)
             end if
             opoles_grad_info(i)%m=points_grad_info(i)%m
          end do

          if(do_alloc) then
             allocate(opoles_grad_index(N_centers+1),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          opoles_grad_index=points_grad_index
          totsym_grad_oct_length=counter-1
       else if (do_R) then
          if(.not. associated(rep_grad_info)) then
             allocate(rep_grad_info(N_centers),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             do_alloc=.true.
          end if
          do i=1,N_rc
             m1=size(points_grad_info(i)%m,1)
             m2=size(points_grad_info(i)%m,2)
             m3=size(points_grad_info(i)%m,3)
             if(do_alloc) then
                allocate(rep_grad_info(i)%m(m1,m2,m3),stat=alloc_stat)
                ASSERT(alloc_stat==0)
             end if
             rep_grad_info(i)%m=points_grad_info(i)%m
          end do

          if(do_alloc) then
             allocate(rep_grad_index(N_centers+1),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          rep_grad_index=points_grad_index
          totsym_grad_rep_length=counter-1
       else if (do_Z) then
          if(.not. associated(unique_pc_grad_info)) then
             allocate(unique_pc_grad_info(N_centers),stat=alloc_stat)
             ASSERT(alloc_stat==0)
             do_alloc=.true.
          end if
          do i=1,pointcharge_N
             m1=size(points_grad_info(i)%m,1)
             m2=size(points_grad_info(i)%m,2)
             m3=size(points_grad_info(i)%m,3)
             if(do_alloc) then
                allocate(unique_pc_grad_info(i)%m(m1,m2,m3),stat=alloc_stat)
                ASSERT(alloc_stat==0)
             end if
             unique_pc_grad_info(i)%m=points_grad_info(i)%m
          end do

          if(do_alloc) then
             allocate(unique_pc_index(N_centers+1),STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          unique_pc_index=points_grad_index
          totsym_grad_pc_length=counter-1
       end if

       do i=1,N_centers
          deallocate(points_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(points_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       deallocate(points_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    do i=1,size(unique_external_points)
       ua=>unique_external_points(i)
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
    deallocate(unique_external_points,STAT=alloc_stat)
    ASSERT(alloc_stat==0)

  end subroutine external_points_grad_inform
  !*************************************************************

  !*************************************************************
  subroutine dealloc_dqo
    !  Purpose: deallocate allocated external centers
    !------------ Modules used -----------------------------------
    use comm_module, only: comm_i_am_master
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,status
    !------------ Executable code ------------------------------------
    if(associated(pd_array)) then

       do i=1,N_pd
          deallocate(pd_array(i)%position,stat=status)
          ASSERT(status==0)
          deallocate(pd_array(i)%dipole,stat=status)
          ASSERT(status==0)
          if(comm_i_am_master()) then
             deallocate(pd_array(i)%group,stat=status)
             ASSERT(status==0)
          end if
       end do

       deallocate(pd_array,stat=status)
       ASSERT(status==0)
    end if

    if(associated(pq_array)) then

       do i=1,N_pq
          deallocate(pq_array(i)%position,stat=status)
          ASSERT(status==0)
          deallocate(pq_array(i)%quadrupole,stat=status)
          ASSERT(status==0)
          if(comm_i_am_master()) then
             deallocate(pq_array(i)%group,stat=status)
             ASSERT(status==0)
          end if
       end do

       deallocate(pq_array,stat=status)
       ASSERT(status==0)
    end if

    if(associated(po_array)) then

       do i=1,N_po
          deallocate(po_array(i)%position,stat=status)
          ASSERT(status==0)
          deallocate(po_array(i)%octopole,stat=status)
          ASSERT(status==0)
          if(comm_i_am_master()) then
             deallocate(po_array(i)%group,stat=status)
             ASSERT(status==0)
          end if
       end do

       deallocate(po_array,stat=status)
       ASSERT(status==0)
    end if

    if(associated(rc_array)) then

       do i=1,N_rc
          deallocate(rc_array(i)%position,stat=status)
          ASSERT(status==0)
          if(comm_i_am_master()) then
             deallocate(rc_array(i)%group,stat=status)
             ASSERT(status==0)
          end if
       end do

       deallocate(rc_array,stat=status)
       ASSERT(status==0)
    end if

  end subroutine dealloc_dqo
  !*************************************************************

  !*************************************************************
  subroutine dealloc_center_inform()
    ! Purpose: deallocate the variables 'point_grad_info'
    !** End of interface *****************************************
    ! ----------- declaration of local variables -------------
    integer(kind=i4_kind) :: i,alloc_stat
    ! ----------- executable code -----------------------------

    if(associated(dipoles_grad_info)) then
       do i=1,size(dipoles_grad_info)
          deallocate(dipoles_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(dipoles_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(dipoles_grad_info)

       deallocate(dipoles_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if(associated(qpoles_grad_info)) then
       do i=1,size(qpoles_grad_info)
          deallocate(qpoles_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(qpoles_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(qpoles_grad_info)

       deallocate(qpoles_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if(associated(opoles_grad_info)) then
       do i=1,size(opoles_grad_info)
          deallocate(opoles_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(opoles_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(opoles_grad_info)

       deallocate(opoles_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if(associated(rep_grad_info)) then
       do i=1,size(rep_grad_info)
          deallocate(rep_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(rep_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(rep_grad_info)

       deallocate(rep_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

  end subroutine dealloc_center_inform
  !*************************************************************

  !*************************************************************
  subroutine output_geometry_ec(type_of_center)
    ! Purpose: write the geomtry of external centers to output
    use iounitadmin_module, only : output_unit
    use constants, only: angstrom
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=2) :: type_of_center
    !--------------------------------------------------
    integer(i4_kind)            :: i_eq,i_ua
    real(r8_kind),parameter     :: d2au=2.541766_r8_kind
    !--- executable code--------------------------------

    if(type_of_center=='PC' .and. pointcharge_N > 0) then
       write(output_unit,*)"-- Geometry of Point Charges in au ------------------------------------"

       do i_ua=1,pointcharge_N
          do i_eq=1,pointcharge_array(i_ua)%N_equal_charges
             write(output_unit,900)pointcharge_array(i_ua)%Z,&
                  trim(pointcharge_array(i_ua)%name),&
                  pointcharge_array(i_ua)%position(1,i_eq),&
                  pointcharge_array(i_ua)%position(2,i_eq),&
                  pointcharge_array(i_ua)%position(3,i_eq)
          enddo
       enddo
    end if

    if(type_of_center=='PD' .and. N_pd > 0) then
       write(output_unit,*)"-- Geometry of Point Dipoles in au ------------------------------------"
       write(output_unit,*)"            dipoles                                        coordinates"

       do i_ua=1,N_pd
          do i_eq=1,pd_array(i_ua)%N_equal_dipoles
             write(output_unit,1000) pd_array(i_ua)%dipole(1:3,i_eq),&
                  trim(pd_array(i_ua)%name),&
                  pd_array(i_ua)%position(1,i_eq),&
                  pd_array(i_ua)%position(2,i_eq),&
                  pd_array(i_ua)%position(3,i_eq)
          enddo
       enddo

       write(output_unit,*)" "
       write(output_unit,*)"-- Geometry of Point Dipoles in Debye, Angstroms ----------------------"
       write(output_unit,*)"            dipoles                                        coordinates"
       do i_ua=1,N_pd
          do i_eq=1,pd_array(i_ua)%N_equal_dipoles
             write(output_unit,1100)pd_array(i_ua)%dipole(1:3,i_eq)/d2au,&
                  trim(pd_array(i_ua)%name),&
                  pd_array(i_ua)%position(1, i_eq) / angstrom,&
                  pd_array(i_ua)%position(2, i_eq) / angstrom,&
                  pd_array(i_ua)%position(3, i_eq) / angstrom
          enddo
       enddo
       write(output_unit,*)" "
       write(output_unit,*)"-----------------------------------------------------------------------"
    end if

    if(type_of_center=='PQ' .and. N_pq > 0) then
       write(output_unit,*)"-- Geometry of Point Quadrupoles in au --------------------------------"
       write(output_unit,*)"          quadrupoles                                      coordinates"

       do i_ua=1,N_pq
          do i_eq=1,pq_array(i_ua)%N_equal_qpoles
             write(output_unit,1000) pq_array(i_ua)%quadrupole(1,1:3,i_eq),&
                  trim(pq_array(i_ua)%name),&
                  pq_array(i_ua)%position(1,i_eq),&
                  pq_array(i_ua)%position(2,i_eq),&
                  pq_array(i_ua)%position(3,i_eq)
             write(output_unit,1200) pq_array(i_ua)%quadrupole(2,1:3,i_eq)
             write(output_unit,1200) pq_array(i_ua)%quadrupole(3,1:3,i_eq)
          enddo
       enddo
    end if

    if(type_of_center=='PO' .and. N_po > 0) then
       write(output_unit,*)"-- Geometry of Point Octopoles in au ----------------------------------"
       write(output_unit,*)"           octopoles                                       coordinates"

       do i_ua=1,N_po
          do i_eq=1,po_array(i_ua)%N_equal_opoles
             write(output_unit,1000) po_array(i_ua)%octopole(1,1:3,1,i_eq),&
                  trim(po_array(i_ua)%name),&
                  po_array(i_ua)%position(1,i_eq),&
                  po_array(i_ua)%position(2,i_eq),&
                  po_array(i_ua)%position(3,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(2,1:3,1,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(3,1:3,1,i_eq)

             write(output_unit,1200) po_array(i_ua)%octopole(1,1:3,2,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(2,1:3,2,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(3,1:3,2,i_eq)

             write(output_unit,1200) po_array(i_ua)%octopole(1,1:3,3,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(2,1:3,3,i_eq)
             write(output_unit,1200) po_array(i_ua)%octopole(3,1:3,3,i_eq)
          enddo
       enddo
    end if

900 format('au         ',F5.1,2X,A4,2X,3(F11.6,7X))
1000 format('',3(F9.4,2X),2X,A4,2X,3(F11.6,7X))
1100 format('',3(F9.4,2X),2X,A4,2X,3(F11.6,7X))
1200 format('',3(F9.4,2X))
  end subroutine output_geometry_ec
  !*************************************************************

  !*****************************************************************************
  subroutine external_centers_bcast()
    !--------------------------------------------------------------------------+
    !  Purpose: broadcasts all information about external points
    !  called by send_recv_init_options
    use comm,               only: comm_bcast                                   &
                                , comm_rank
    !** End of interface *******************************************************
    !----------- declaration of local variables -------------------------------+
    type(pointdipole_type),     pointer :: pd
    type(pointquadrupole_type), pointer :: pq
    type(pointoctopole_type),   pointer :: po
    type(repulsion_type),       pointer :: rc
    integer(i4_kind)                    :: i, status
    logical                             :: do_alloc
    !--- executable code ------------------------------------------------------+
    !
    call comm_bcast( moving_X_centers )
    call comm_bcast( moving_R_centers )
#ifdef WITH_EFP
    call comm_bcast( efp              )
#endif
    !
    call comm_bcast( N_pd             )
    if( N_pd > 0 ) then
      !
      do_alloc=.false.
      !
      if( comm_rank() /= 0 .and. .not. associated(pd_array) ) then
        allocate( pd_array(N_pd), stat=status )
        ASSERT(status==0)
        do_alloc=.true.
      end if
      !
      do i = 1, N_pd
        pd => pd_array(i)
        !
        call comm_bcast( pd%name            )
        call comm_bcast( pd%N_equal_dipoles )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( pd%position(3,pd%N_equal_dipoles), stat=status )
          ASSERT(status==0)
          allocate( pd%dipole(3,pd%N_equal_dipoles)  , stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( pd%position        )
        call comm_bcast( pd%dipole          )
       end do
       !
       if(moving_X_centers) call external_points_info_bcast( "PD" )
       !
    end if
    !
    call comm_bcast( N_pq )
    if( N_pq > 0 ) then
      !
      do_alloc=.false.
      !
      if( comm_rank() /= 0 .and. .not. associated(pq_array) ) then
        allocate( pq_array(N_pq), stat=status )
        ASSERT(status==0)
        do_alloc=.true.
      end if
      !
      do i = 1, N_pq
        pq => pq_array(i)
        !
        call comm_bcast( pq%name           )
        call comm_bcast( pq%N_equal_qpoles )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( pq%position(3,pq%N_equal_qpoles)    , stat=status )
          ASSERT(status==0)
          allocate( pq%quadrupole(3,3,pq%N_equal_qpoles), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( pq%position       )
        call comm_bcast( pq%quadrupole     )
      end do
      !
      if(moving_X_centers) call external_points_info_bcast( "PQ" )
      !
    end if
    !
    call comm_bcast( N_po )
    if( N_po > 0 ) then
      !
      do_alloc=.false.
      !
      if( comm_rank() /= 0 .and. .not. associated(po_array) ) then
        allocate( po_array(N_po), stat=status )
        ASSERT(status==0)
        do_alloc=.true.
      end if
      !
      do i = 1, N_po
        po => po_array(i)
        !
        call comm_bcast( po%name           )
        call comm_bcast( po%N_equal_opoles )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( po%position(3,po%N_equal_opoles)    , stat=status )
          ASSERT(status==0)
          allocate( po%octopole(3,3,3,po%N_equal_opoles), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( po%position       )
        call comm_bcast( po%octopole       )
      end do
      !
      if(moving_X_centers) call external_points_info_bcast( "PO" )
      !
    end if
    !
    call comm_bcast( N_rc )
    !
    if( N_rc > 0 ) then
      !
      do_alloc=.false.
      !
      if( comm_rank() /= 0 .and. .not. associated(rc_array) ) then
        allocate( rc_array(N_rc), stat=status )
        ASSERT(status==0)
        do_alloc=.true.
      end if
      !
      do i = 1, N_rc
        rc => rc_array(i)
        !
        call comm_bcast( rc%name            )
        call comm_bcast( rc%N_equal_centers )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( rc%position(3,rc%N_equal_centers), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( rc%position        )
        !
        call comm_bcast( rc%C,status        )
        call comm_bcast( rc%A,status        )
       end do
       !
       if(moving_R_centers) call external_points_info_bcast( "RC" )
       !
    end if
    !
    contains
    !
    subroutine external_points_info_bcast( type_of_centers )
      !------------ Declaration of formal parameters ---------------------------
      character(len=2)      :: type_of_centers
      !** End of interface *****************************************************
      !------------ Declaration of local variables -----------------------------
      integer(i4_kind) :: i, status
      integer(i4_kind) :: n1, n2, n3
      logical          :: do_alloc
      !------------ Executable code --------------------------------------------
      !
      if( type_of_centers == "PD" ) then
        !
        do_alloc=.false.
        !
        if( comm_rank() /= 0 .and. .not. associated(dipoles_grad_info) ) then
          allocate( dipoles_grad_info(N_pd), stat=status )
          ASSERT(status==0)
          do_alloc=.true.
        end if
        !
        do i = 1, N_pd
          !
          if( comm_rank() == 0 ) then
            n1 = size( dipoles_grad_info(i)%m, 1 )
            n2 = size( dipoles_grad_info(i)%m, 2 )
            n3 = size( dipoles_grad_info(i)%m, 3 )
          endif
          !
          call comm_bcast( n1 )
          call comm_bcast( n2 )
          call comm_bcast( n3 )
          !
          if( comm_rank() /= 0 .and. do_alloc ) then
            allocate( dipoles_grad_info(i)%m(n1,n2,n3), stat=status)
            ASSERT(status==0)
          end if
          call comm_bcast( dipoles_grad_info(i)%m )
          !
        end do
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( dipoles_grad_index(N_pd+1), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( dipoles_grad_index     )
        call comm_bcast( totsym_grad_dip_length )
        !
      elseif(type_of_centers == "PQ") then
        !
        do_alloc=.false.
        !
        if( comm_rank() /= 0 .and. .not. associated(qpoles_grad_info) ) then
          allocate( qpoles_grad_info(N_pq), stat=status )
          ASSERT(status==0)
          do_alloc=.true.
        end if
        !
        do i = 1, N_pq
          !
          if( comm_rank() == 0 ) then
            n1 = size( qpoles_grad_info(i)%m, 1 )
            n2 = size( qpoles_grad_info(i)%m, 2 )
            n3 = size( qpoles_grad_info(i)%m, 3 )
          endif
          !
          call comm_bcast( n1 )
          call comm_bcast( n2 )
          call comm_bcast( n3 )
          !
          if( comm_rank() /= 0 .and. do_alloc ) then
            allocate( qpoles_grad_info(i)%m(n1,n2,n3), stat=status )
            ASSERT(status==0)
          end if
          !
          call comm_bcast( qpoles_grad_info(i)%m )
          !
        end do
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( qpoles_grad_index(N_pq+1), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( qpoles_grad_index       )
        call comm_bcast( totsym_grad_quad_length )
        !
      elseif(type_of_centers == "PO") then
        !
        do_alloc=.false.
        !
        if( comm_rank() /= 0 .and. .not. associated(opoles_grad_info) ) then
          allocate( opoles_grad_info(N_po), stat=status )
          ASSERT(status==0)
          do_alloc=.true.
        end if
        !
        do i = 1, N_po
          !
          if( comm_rank() == 0 ) then
            n1 = size( opoles_grad_info(i)%m, 1 )
            n2 = size( opoles_grad_info(i)%m, 2 )
            n3 = size( opoles_grad_info(i)%m, 3 )
          endif
          !
          call comm_bcast( n1 )
          call comm_bcast( n2 )
          call comm_bcast( n3 )
          !
          if( comm_rank() /= 0 .and. do_alloc ) then
            allocate( opoles_grad_info(i)%m(n1,n2,n3), stat=status )
            ASSERT(status==0)
          end if
          !
          call comm_bcast( opoles_grad_info(i)%m )
          !
        end do
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( opoles_grad_index(N_po+1), stat=status )
          ASSERT(status==0)
        end if
        call comm_bcast( opoles_grad_index      )
        call comm_bcast( totsym_grad_oct_length )
        !
      elseif(type_of_centers == "RC") then
        !
        do_alloc=.false.
        !
        if( comm_rank() /= 0 .and. .not. associated(rep_grad_info) ) then
          allocate( rep_grad_info(N_rc), stat=status )
          ASSERT(status==0)
          do_alloc=.true.
        end if
        !
        do i = 1, N_rc
          !
          if( comm_rank() == 0 ) then
            n1 = size( rep_grad_info(i)%m, 1 )
            n2 = size( rep_grad_info(i)%m, 2 )
            n3 = size( rep_grad_info(i)%m, 3 )
          endif
          !
          call comm_bcast( n1 )
          call comm_bcast( n2 )
          call comm_bcast( n3 )
          !
          if( comm_rank() /= 0 .and. do_alloc ) then
            allocate( rep_grad_info(i)%m(n1,n2,n3), stat=status )
            ASSERT(status==0)
          end if
          !
          call comm_bcast( rep_grad_info(i)%m )
          !
        end do
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( rep_grad_index(N_rc+1), stat=status )
          ASSERT(status==0)
        end if
        !
         call comm_bcast( rep_grad_index         )
         call comm_bcast( totsym_grad_rep_length )
         !
      end if
      !
    end subroutine external_points_info_bcast
    !
  end subroutine external_centers_bcast
  !*****************************************************************************

  !*************************************************************
  subroutine init_X_centers_grads()
    ! Purpose: allocate arrays to calculate garadients on X centers
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: status,i,n_eq
     !------------ Executable code --------------------------------

     if(moving_X_centers) then
     if(N_pd > 0) then
        allocate(gradient_dip_cartesian(N_pd),stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_dip_cartesian(N_pd),stat=status)
           ASSERT(status==0)
        end if
#endif

        do i=1,N_pd
           n_eq=pd_array(i)%N_equal_dipoles
           allocate(gradient_dip_cartesian(i)%m(3,n_eq),stat=status)
           ASSERT(status==0)
           gradient_dip_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
           if(efp) then
              allocate(torque_dip_cartesian(i)%m(3,n_eq),stat=status)
              ASSERT(status==0)
              torque_dip_cartesian(i)%m=0.0_r8_kind
           end if
#endif
        end do

        allocate(gradient_dip_totalsym(totsym_grad_dip_length),stat=status)
        ASSERT(status==0)
        gradient_dip_totalsym=0.0_r8_kind
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_dip_totalsym(totsym_grad_dip_length),stat=status)
           ASSERT(status==0)
           torque_dip_totalsym=0.0_r8_kind
        end if
#endif
     end if

     if(N_pq > 0) then
        allocate(gradient_quad_cartesian(N_pq),stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_quad_cartesian(N_pq),stat=status)
           ASSERT(status==0)
        end if
#endif

        do i=1,N_pq
           n_eq=pq_array(i)%N_equal_qpoles
           allocate(gradient_quad_cartesian(i)%m(3,n_eq),stat=status)
           ASSERT(status==0)
           gradient_quad_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
           if(efp) then
              allocate(torque_quad_cartesian(i)%m(3,n_eq),stat=status)
              ASSERT(status==0)
              torque_quad_cartesian(i)%m=0.0_r8_kind
           end if
#endif
        end do

        allocate(gradient_quad_totalsym(totsym_grad_quad_length),stat=status)
        ASSERT(status==0)
        gradient_quad_totalsym=0.0_r8_kind
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_quad_totalsym(totsym_grad_quad_length),stat=status)
           ASSERT(status==0)
           torque_quad_totalsym=0.0_r8_kind
        end if
#endif
     end if

     if(N_po > 0) then
        allocate(gradient_oct_cartesian(N_po),stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_oct_cartesian(N_po),stat=status)
        end if
#endif

        do i=1,N_po
           n_eq=po_array(i)%N_equal_opoles
           allocate(gradient_oct_cartesian(i)%m(3,n_eq),stat=status)
           ASSERT(status==0)
           gradient_oct_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
           if(efp) then
              allocate(torque_oct_cartesian(i)%m(3,n_eq),stat=status)
              ASSERT(status==0)
              torque_oct_cartesian(i)%m=0.0_r8_kind
           end if
#endif
        end do

        allocate(gradient_oct_totalsym(totsym_grad_oct_length),stat=status)
        ASSERT(status==0)
        gradient_oct_totalsym=0.0_r8_kind
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_oct_totalsym(totsym_grad_oct_length),stat=status)
           ASSERT(status==0)
           torque_oct_totalsym=0.0_r8_kind
        end if
#endif
     end if
     end if

     if(moving_R_centers) then
     if(N_rc > 0) then
        allocate(gradient_rep_cartesian(N_rc),stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_rep_cartesian(N_rc),stat=status)
           ASSERT(status==0)
        end if
#endif

        do i=1,N_rc
           n_eq=rc_array(i)%N_equal_centers
           allocate(gradient_rep_cartesian(i)%m(3,n_eq),stat=status)
           ASSERT(status==0)
           gradient_rep_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
           if(efp) then
              allocate(torque_rep_cartesian(i)%m(3,n_eq),stat=status)
              ASSERT(status==0)
              torque_rep_cartesian(i)%m=0.0_r8_kind
           end if
#endif
        end do

        allocate(gradient_rep_totalsym(totsym_grad_rep_length),stat=status)
        ASSERT(status==0)
        gradient_rep_totalsym=0.0_r8_kind
     end if
     end if

  end subroutine init_X_centers_grads
  !*************************************************************

  !*************************************************************
  subroutine X_centers_grads_shutdown()
     !  Purpose: free arrays for calculating gradients on X centers
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: status,i
     !------------ Executable code --------------------------------

     if(moving_X_centers) then
     if(N_pd > 0) then
        do i=1,N_pd
           deallocate(gradient_dip_cartesian(i)%m,stat=status)
           ASSERT(status==0)
#ifdef WITH_EFP
           if(efp) then
              deallocate(torque_dip_cartesian(i)%m,stat=status)
              ASSERT(status==0)
           end if
#endif
        end do
        deallocate(gradient_dip_cartesian,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_dip_cartesian,stat=status)
           ASSERT(status==0)
        end if
#endif

        deallocate(gradient_dip_totalsym,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_dip_totalsym,stat=status)
           ASSERT(status==0)
        end if
#endif
     end if

     if(N_pq > 0) then
        do i=1,N_pq
           deallocate(gradient_quad_cartesian(i)%m,stat=status)
           ASSERT(status==0)
#ifdef WITH_EFP
           if(efp) then
              deallocate(torque_quad_cartesian(i)%m,stat=status)
              ASSERT(status==0)
           end if
#endif
        end do
        deallocate(gradient_quad_cartesian,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_quad_cartesian,stat=status)
           ASSERT(status==0)
        end if
#endif

        deallocate(gradient_quad_totalsym,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_quad_totalsym,stat=status)
           ASSERT(status==0)
        end if
#endif
     end if

     if(N_po > 0) then
        do i=1,N_po
           deallocate(gradient_oct_cartesian(i)%m,stat=status)
           ASSERT(status==0)
#ifdef WITH_EFP
           if(efp) then
              deallocate(torque_oct_cartesian(i)%m,stat=status)
              ASSERT(status==0)
           end if
#endif
        end do
        deallocate(gradient_oct_cartesian,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_oct_cartesian,stat=status)
           ASSERT(status==0)
        end if
#endif

        deallocate(gradient_oct_totalsym,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_oct_totalsym,stat=status)
           ASSERT(status==0)
        end if
#endif
     end if
     end if

     if(moving_R_centers) then
     if(N_rc > 0) then
        do i=1,N_rc
           deallocate(gradient_rep_cartesian(i)%m,stat=status)
           ASSERT(status==0)
#ifdef WITH_EFP
           if(efp) then
              deallocate(torque_rep_cartesian(i)%m,stat=status)
              ASSERT(status==0)
           end if
#endif
        end do
        deallocate(gradient_rep_cartesian,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_rep_cartesian,stat=status)
           ASSERT(status==0)
        end if
#endif

        deallocate(gradient_rep_totalsym,stat=status)
        ASSERT(status==0)
     end if
     end if

  end subroutine X_centers_grads_shutdown
  !*************************************************************

  !*************************************************************
  subroutine totsym_X_grad_unpack()
    !  Purpose: unpack 2 center contribution to the X center gradient from the
    !           slave
    use comm_module, only: communpack
    !------------ Declaration of local variables ----------------
    real(r8_kind) :: help_arr_pd(totsym_grad_dip_length)
    real(r8_kind) :: help_arr_pq(totsym_grad_quad_length)
    real(r8_kind) :: help_arr_po(totsym_grad_oct_length)
    real(r8_kind) :: help_arr_rc(totsym_grad_rep_length)
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    if(moving_X_centers) then
    if(N_pd > 0) then
       call communpack(help_arr_pd,totsym_grad_dip_length,1,info)
       ASSERT(info==0)
       gradient_dip_totalsym=gradient_dip_totalsym+help_arr_pd
#ifdef WITH_EFP
       if(efp) then
          call communpack(help_arr_pd,totsym_grad_dip_length,1,info)
          ASSERT(info==0)
          torque_dip_totalsym=torque_dip_totalsym+help_arr_pd
       end if
#endif
    end if
    if(N_pq > 0) then
       call communpack(help_arr_pq,totsym_grad_quad_length,1,info)
       ASSERT(info==0)
       gradient_quad_totalsym=gradient_quad_totalsym+help_arr_pq
#ifdef WITH_EFP
       if(efp) then
          call communpack(help_arr_pq,totsym_grad_quad_length,1,info)
          ASSERT(info==0)
          torque_quad_totalsym=torque_quad_totalsym+help_arr_pq
       end if
#endif
    end if
    if(N_po > 0) then
       call communpack(help_arr_po,totsym_grad_oct_length,1,info)
       ASSERT(info==0)
       gradient_oct_totalsym=gradient_oct_totalsym+help_arr_po
#ifdef WITH_EFP
       if(efp) then
          call communpack(help_arr_po,totsym_grad_oct_length,1,info)
          ASSERT(info==0)
          torque_oct_totalsym=torque_oct_totalsym+help_arr_po
       end if
#endif
    end if
    end if
    if(moving_R_centers) then
    if(N_rc > 0) then
       call communpack(help_arr_rc,totsym_grad_rep_length,1,info)
       ASSERT(info==0)
       gradient_rep_totalsym=gradient_rep_totalsym+help_arr_rc
    end if
    end if

  end subroutine totsym_X_grad_unpack
  !*************************************************************

  !*************************************************************
  subroutine totsym_X_grad_pack()
    !  Purpose: pack 2 center contribution to the X center gradient on the
    !           slave
    use comm_module, only: commpack
    !------------ Declaration of local variables ----------------
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    if(moving_X_centers) then
    if(N_pd > 0) then
       call commpack(gradient_dip_totalsym,totsym_grad_dip_length,1,info)
       ASSERT(info==0)
#ifdef WITH_EFP
       if(efp) then
          call commpack(torque_dip_totalsym,totsym_grad_dip_length,1,info)
          ASSERT(info==0)
       end if
#endif
    end if
    if(N_pq > 0) then
       call commpack(gradient_quad_totalsym,totsym_grad_quad_length,1,info)
       ASSERT(info==0)
#ifdef WITH_EFP
       if(efp) then
          call commpack(torque_quad_totalsym,totsym_grad_quad_length,1,info)
          ASSERT(info==0)
       end if
#endif
    end if
    if(N_po > 0) then
       call commpack(gradient_oct_totalsym,totsym_grad_oct_length,1,info)
       ASSERT(info==0)
#ifdef WITH_EFP
       if(efp) then
          call commpack(torque_oct_totalsym,totsym_grad_oct_length,1,info)
          ASSERT(info==0)
       end if
#endif
    end if
    end if
    if(moving_R_centers) then
    if(N_rc > 0) then
       call commpack(gradient_rep_totalsym,totsym_grad_rep_length,1,info)
       ASSERT(info==0)
    end if
    end if

  end subroutine totsym_X_grad_pack
  !*************************************************************

  !*************************************************************
  function calc_nuc_dqo_energy() result (e_nuc_ext)
    ! Purpose: Calculate interaction energy between atomic nuclears
    ! end external centers (dipoles, quadrupoles, octopoles)
    use unique_atom_module, only : unique_atoms, N_unique_atoms, &
         pseudopot_present
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: e_nuc_ext
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ua,i_ea,i_uc,i_ec,i,j,k
    integer(i4_kind) :: n_ec,n_ea
    real(r8_kind),pointer  :: xa(:,:),xc(:,:)
    real(r8_kind),pointer  :: dipole(:,:)
    real(r8_kind),pointer  :: quadrupole(:,:,:)
    real(r8_kind),pointer  :: octopole(:,:,:,:)
    real(r8_kind) :: Z, Zc,r_ac(3),d_ac,Eij,Eijk
    real(r8_kind) :: e_nuc_dip,e_nuc_quad,e_nuc_oct
    !--- executable code-------------------------------------

    e_nuc_ext=0.0_r8_kind

    e_nuc_dip=0.0_r8_kind
    if(N_pd > 0) then
       do i_uc=1,N_pd
          n_ec=pd_array(i_uc)%N_equal_dipoles
          xc=>pd_array(i_uc)%position
          dipole=>pd_array(i_uc)%dipole

          do i_ec=1,n_ec

             do i_ua=1,N_unique_atoms
                Z=unique_atoms(i_ua)%Z
                Zc=unique_atoms(i_ua)%Zc
                if (.not.pseudopot_present) Zc = 0.0_r8_kind
                Z=Z-Zc
                xa=>unique_atoms(i_ua)%position
                n_ea=unique_atoms(i_ua)%N_equal_atoms

                do i_ea=1,n_ea
                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind
                   e_nuc_dip=e_nuc_dip+ &
                        Z*dot_product(r_ac,dipole(:,i_ec))/d_ac**3

                end do
             end do
          end do
       end do
    end if

    e_nuc_quad=0.0_r8_kind
    if(N_pq > 0) then
       do i_uc=1,N_pq
          n_ec=pq_array(i_uc)%N_equal_qpoles
          xc=>pq_array(i_uc)%position
          quadrupole=>pq_array(i_uc)%quadrupole

          do i_ec=1,n_ec

             do i_ua=1,N_unique_atoms
                Z=unique_atoms(i_ua)%Z
                Zc=unique_atoms(i_ua)%Zc
                if (.not.pseudopot_present) Zc = 0.0_r8_kind
                Z=Z-Zc
                xa=>unique_atoms(i_ua)%position
                n_ea=unique_atoms(i_ua)%N_equal_atoms

                do i_ea=1,n_ea
                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                   do i=1,3
                      do j=1,3
                         Eij=3.0_r8_kind*r_ac(i)*r_ac(j)/d_ac**5
                         if(i==j) Eij=Eij-1.0_r8_kind/d_ac**3

                         e_nuc_quad=e_nuc_quad+ &
                              Z*quadrupole(i,j,i_ec)*Eij/3.0_r8_kind
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if

    e_nuc_oct=0.0_r8_kind
    if(N_po > 0) then
       do i_uc=1,N_po
          n_ec=po_array(i_uc)%N_equal_opoles
          xc=>po_array(i_uc)%position
          octopole=>po_array(i_uc)%octopole

          do i_ec=1,n_ec

             do i_ua=1,N_unique_atoms
                Z=unique_atoms(i_ua)%Z
                Zc=unique_atoms(i_ua)%Zc
                if (.not.pseudopot_present) Zc = 0.0_r8_kind
                Z=Z-Zc
                xa=>unique_atoms(i_ua)%position
                n_ea=unique_atoms(i_ua)%N_equal_atoms

                do i_ea=1,n_ea
                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                   do i=1,3
                      do j=1,3
                         do k=1,3
                            Eijk=15.0_r8_kind*r_ac(i)*r_ac(j)*r_ac(k)/d_ac**7
                            if(i==j) Eijk=Eijk-3.0_r8_kind*r_ac(k)/d_ac**5
                            if(i==k) Eijk=Eijk-3.0_r8_kind*r_ac(j)/d_ac**5
                            if(j==k) Eijk=Eijk-3.0_r8_kind*r_ac(i)/d_ac**5

                            e_nuc_oct=e_nuc_oct+ &
                                 Z*octopole(i,j,k,i_ec)*Eijk/15.0_r8_kind
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if
#ifdef WITH_EFP
    if((N_pd > 0) .and. efp .and. print_energy) print*,'NUC-DIP  EFP ',e_nuc_dip
    if((N_pq > 0) .and. efp .and. print_energy) print*,'NUC-QUAD EFP ',e_nuc_quad
    if((N_po > 0) .and. efp .and. print_energy) print*,'NUC-OCT  EFP ',e_nuc_oct
#endif
    e_nuc_ext=e_nuc_dip+e_nuc_quad+e_nuc_oct
!print*,e_nuc_ext,'###############'

  end function calc_nuc_dqo_energy
  !*************************************************************

  !*************************************************************
  subroutine calc_nuc_dqo_grads(totsym_grad,gradient_index)
    ! Purpose: Calculate gradients of energy between atomic nuclears
    ! end external centers (dipoles, quadrupoles, octopoles)
    use unique_atom_module, only : unique_atoms, &
         pseudopot_present,N_moving_unique_atoms,moving_unique_atom_index, &
         unique_atom_grad_info
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(inout) :: totsym_grad(:)
    integer(i4_kind), intent(in) :: gradient_index(N_moving_unique_atoms + 1)
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ua,i_ea,i_uc,i_ec,ma,i,j,k,l
    integer(i4_kind) :: n_ec,n_ea
    real(r8_kind),pointer  :: xa(:,:),xc(:,:),rotmat(:,:)
    real(r8_kind),pointer  :: dipole(:,:)
    real(r8_kind),pointer  :: quadrupole(:,:,:)
    real(r8_kind),pointer  :: octopole(:,:,:,:)
    real(r8_kind) :: gradient(3),Gik,Gkj,Gii
    real(r8_kind) :: Gijk,Gijkl,G1ijk,G1ijkl
    real(r8_kind) :: Z, Zc,r_ac(3),d_ac
    integer(i4_kind) :: grad_dim,index
    !--- executable code-------------------------------------

    uniq_at: do ma=1,N_moving_unique_atoms
       i_ua=moving_unique_atom_index(ma)
       n_ea=unique_atoms(i_ua)%n_equal_atoms
       Z=unique_atoms(i_ua)%Z
       Zc=unique_atoms(i_ua)%Zc
       if (.not.pseudopot_present) Zc = 0.0_r8_kind
       Z=Z-Zc
       xa=>unique_atoms(i_ua)%position

       grad_dim=gradient_index(ma+1)-gradient_index(ma)
       index=gradient_index(ma)-1

       eq_at: do i_ea=1,n_ea
          rotmat=>unique_atom_grad_info(ma)%m(:,:,i_ea)
          gradient=0.0_r8_kind

          Point_Dipoles: if(N_pd > 0) then
             uniq_pd: do i_uc=1,N_pd
                n_ec=pd_array(i_uc)%N_equal_dipoles
                xc=>pd_array(i_uc)%position
                dipole=>pd_array(i_uc)%dipole

                eq_pd: do i_ec=1,n_ec
                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                   gradient=gradient+(Z*dipole(:,i_ec)/d_ac**3- &
                        3.0_r8_kind*Z*dot_product(r_ac,dipole(:,i_ec))*r_ac/d_ac**5)
                end do eq_pd
             end do uniq_pd
          end if Point_Dipoles

          Point_Quadrupoles: if(N_pq > 0) then
             uniq_pq: do i_uc=1,N_pq
                n_ec=pq_array(i_uc)%N_equal_qpoles
                xc=>pq_array(i_uc)%position
                quadrupole=>pq_array(i_uc)%quadrupole

                eq_pq: do i_ec=1,n_ec

                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                   do k=1,3
                      do i=1,3
                         do j=1,3
                            Gkj=0.0_r8_kind
                            if(k==i) Gkj=r_ac(j)*quadrupole(k,j,i_ec)
                            Gik=0.0_r8_kind
                            if(k==j) Gik=r_ac(i)*quadrupole(i,k,i_ec)
                            Gii=0.0_r8_kind
                            if(i==j) Gii=r_ac(k)*quadrupole(i,j,i_ec)
                            gradient(k)=gradient(k)+Z*( &
                                 -5.0_r8_kind*r_ac(k)*r_ac(i)*r_ac(j)*quadrupole(i,j,i_ec)/d_ac**7+ &
                                 (Gkj+Gik+Gii)/d_ac**5)
                         end do
                      end do
                   end do
                end do eq_pq
             end do uniq_pq
          end if Point_Quadrupoles

          Point_Octopoles: if(N_po > 0) then
             uniq_po: do i_uc=1,N_po
                n_ec=po_array(i_uc)%N_equal_opoles
                xc=>po_array(i_uc)%position
                octopole=>po_array(i_uc)%octopole

                eq_po: do i_ec=1,n_ec

                   r_ac=xa(:,i_ea)-xc(:,i_ec)
                   d_ac=sqrt(dot_product(r_ac,r_ac))
                   if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                   do l=1,3
                      do i=1,3
                         do j=1,3
                            do k=1,3
                               Gijk=r_ac(i)*r_ac(j)*r_ac(k)*octopole(i,j,k,i_ec)

                               Gijkl=0.0_r8_kind
                               if(l==i) Gijkl=Gijkl+r_ac(j)*r_ac(k)
                               if(l==j) Gijkl=Gijkl+r_ac(i)*r_ac(k)
                               if(l==k) Gijkl=Gijkl+r_ac(i)*r_ac(j)
                               Gijkl=Gijkl*octopole(i,j,k,i_ec)

                               G1ijk=0.0_r8_kind
                               if(i==j) G1ijk=G1ijk+r_ac(k)
                               if(i==k) G1ijk=G1ijk+r_ac(j)
                               if(j==k) G1ijk=G1ijk+r_ac(i)
                               G1ijk=G1ijk*octopole(i,j,k,i_ec)

                               G1ijkl=0.0_r8_kind
                               if(i==j .and. k==l) G1ijkl=G1ijkl+1.0_r8_kind
                               if(i==k .and. j==l) G1ijkl=G1ijkl+1.0_r8_kind
                               if(j==k .and. i==l) G1ijkl=G1ijkl+1.0_r8_kind
                               G1ijkl=G1ijkl*octopole(i,j,k,i_ec)

                               gradient(l)=gradient(l)+Z*( &
                                   -35_r8_kind*r_ac(l)*Gijk/d_ac**9+ &
                                   5_r8_kind*Gijkl/d_ac**7+ &
                                   5_r8_kind*r_ac(l)*G1ijk/d_ac**7- &
                                   G1ijkl/d_ac**5)/5_r8_kind
                            end do
                         end do
                      end do
                   end do
                end do eq_po
             end do uniq_po
          end if Point_Octopoles

!print*,ma,i_ea,gradient
          do i=1,grad_dim
             totsym_grad(index+i) = totsym_grad(index+i)+ &
                  sum(rotmat(i,:)*gradient(:))

          enddo
       end do eq_at
    end do uniq_at

  end subroutine calc_nuc_dqo_grads
  !*************************************************************


  !*************************************************************
  subroutine calc_X_grads()
    ! Purpose: Calculate gradients on X centers (interaction with atomic nuclears)
    use unique_atom_module, only : unique_atoms, N_unique_atoms, &
         pseudopot_present
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ux,i_ex,n_ex,i_ua,i_ea,n_ea,i,j,k,l
    real(r8_kind), pointer :: xx(:,:),xa(:,:),rotmat(:,:)
    integer(i4_kind) :: grad_dim,index
    real(r8_kind) :: gradient(3),torque(3),r_ax(3),d_ax
    real(r8_kind),pointer  :: dipole(:,:)
    real(r8_kind),pointer  :: quadrupole(:,:,:)
    real(r8_kind),pointer  :: octopole(:,:,:,:)
    real (r8_kind) :: Z_a, Orr(3)
#ifdef WITH_EFP
    real (r8_kind) :: Qr(3), EE(3)
#endif
    real(r8_kind) :: Gik,Gkj,Gii
    real(r8_kind) :: Gijk,Gijkl,G1ijk,G1ijkl
    !---executable code--------------------------------------

    ! Point dipoles
    do i_ux=1,N_pd
       n_ex=pd_array(i_ux)%N_equal_dipoles
       xx=>pd_array(i_ux)%position
       dipole=>pd_array(i_ux)%dipole

       grad_dim=dipoles_grad_index(i_ux+1)-dipoles_grad_index(i_ux)
       index=dipoles_grad_index(i_ux)-1

       do i_ex=1,n_ex
          rotmat=>dipoles_grad_info(i_ux)%m(:,:,i_ex)
          gradient=0.0_r8_kind
          torque=0.0_r8_kind

          do i_ua=1,N_unique_atoms
             n_ea=unique_atoms(i_ua)%N_equal_atoms
             Z_a=unique_atoms(i_ua)%Z
             if(pseudopot_present) Z_a=Z_a-unique_atoms(i_ua)%Zc
             xa=>unique_atoms(i_ua)%position

             do i_ea=1,n_ea
                r_ax=xa(:,i_ea)-xx(:,i_ex)
                d_ax=sqrt(dot_product(r_ax,r_ax))
                if(d_ax < 0.001_r8_kind) d_ax=0.001_r8_kind

                gradient=gradient+(Z_a*dipole(:,i_ex)/d_ax**3- &
                        3.0_r8_kind*Z_a*dot_product(r_ax,dipole(:,i_ex))*r_ax/d_ax**5)

#ifdef WITH_EFP
                torque=torque+vector_product(dipole(:,i_ex),r_ax)*Z_a/d_ax**3
#endif
             end do
          end do

          do i=1,grad_dim
             gradient_dip_totalsym(index+i)=gradient_dip_totalsym(index+i)- &
                  sum(rotmat(i,:)*gradient(:))

#ifdef WITH_EFP
             torque_dip_totalsym(index+i)=torque_dip_totalsym(index+i)+ &
                  sum(rotmat(i,:)*torque(:))
#endif
          enddo

       end do
    end do

    ! Point quadrupoles
    do i_ux=1,N_pq
       n_ex=pq_array(i_ux)%N_equal_qpoles
       xx=>pq_array(i_ux)%position
       quadrupole=>pq_array(i_ux)%quadrupole

       grad_dim=qpoles_grad_index(i_ux+1)-qpoles_grad_index(i_ux)
       index=qpoles_grad_index(i_ux)-1

       do i_ex=1,n_ex

          rotmat=>qpoles_grad_info(i_ux)%m(:,:,i_ex)
          gradient=0.0_r8_kind
          torque=0.0_r8_kind

          do i_ua=1,N_unique_atoms
             n_ea=unique_atoms(i_ua)%N_equal_atoms
             Z_a=unique_atoms(i_ua)%Z
             if(pseudopot_present) Z_a=Z_a-unique_atoms(i_ua)%Zc
             xa=>unique_atoms(i_ua)%position

             do i_ea=1,n_ea
                r_ax=xa(:,i_ea)-xx(:,i_ex)
                d_ax=sqrt(dot_product(r_ax,r_ax))
                if(d_ax < 0.001_r8_kind) d_ax=0.001_r8_kind

                do k=1,3
                   do i=1,3
                      do j=1,3
                         Gkj=0.0_r8_kind
                         if(k==i) Gkj=r_ax(j)*quadrupole(k,j,i_ex)
                         Gik=0.0_r8_kind
                         if(k==j) Gik=r_ax(i)*quadrupole(i,k,i_ex)
                         Gii=0.0_r8_kind
                         if(i==j) Gii=r_ax(k)*quadrupole(i,j,i_ex)
                         gradient(k)=gradient(k)+Z_a*( &
                              -5.0_r8_kind*r_ax(k)*r_ax(i)*r_ax(j)*quadrupole(i,j,i_ex)/d_ax**7+ &
                              (Gkj+Gik+Gii)/d_ax**5)
                      end do
                   end do
                end do

#ifdef WITH_EFP
                Qr=matmul(quadrupole(:,:,i_ex),r_ax)
                EE=2.0_r8_kind*Z_a*r_ax/d_ax**5
                torque=torque+vector_product(Qr,EE)
#endif
             end do
          end do

          do i=1,grad_dim
             gradient_quad_totalsym(index+i)=gradient_quad_totalsym(index+i)- &
                  sum(rotmat(i,:)*gradient(:))

#ifdef WITH_EFP
             torque_quad_totalsym(index+i)=torque_quad_totalsym(index+i)+ &
                  sum(rotmat(i,:)*torque(:))
#endif
          enddo

       end do
    end do

    ! Point octopoles
    do i_ux=1,N_po
       n_ex=po_array(i_ux)%N_equal_opoles
       xx=>po_array(i_ux)%position
       octopole=>po_array(i_ux)%octopole

       grad_dim=opoles_grad_index(i_ux+1)-opoles_grad_index(i_ux)
       index=opoles_grad_index(i_ux)-1

       do i_ex=1,n_ex

          rotmat=>opoles_grad_info(i_ux)%m(:,:,i_ex)
          gradient=0.0_r8_kind
          torque=0.0_r8_kind

          do i_ua=1,N_unique_atoms
             n_ea=unique_atoms(i_ua)%N_equal_atoms
             Z_a=unique_atoms(i_ua)%Z
             if(pseudopot_present) Z_a=Z_a-unique_atoms(i_ua)%Zc
             xa=>unique_atoms(i_ua)%position

             do i_ea=1,n_ea
                r_ax=xa(:,i_ea)-xx(:,i_ex)
                d_ax=sqrt(dot_product(r_ax,r_ax))
                if(d_ax < 0.001_r8_kind) d_ax=0.001_r8_kind

                Orr=0.0_r8_kind
                do l=1,3
                   do i=1,3
                      do j=1,3
                         do k=1,3
                            Gijk=r_ax(i)*r_ax(j)*r_ax(k)*octopole(i,j,k,i_ex)

                            Gijkl=0.0_r8_kind
                            if(l==i) Gijkl=Gijkl+r_ax(j)*r_ax(k)
                            if(l==j) Gijkl=Gijkl+r_ax(i)*r_ax(k)
                            if(l==k) Gijkl=Gijkl+r_ax(i)*r_ax(j)
                            Gijkl=Gijkl*octopole(i,j,k,i_ex)

                            G1ijk=0.0_r8_kind
                            if(i==j) G1ijk=G1ijk+r_ax(k)
                            if(i==k) G1ijk=G1ijk+r_ax(j)
                            if(j==k) G1ijk=G1ijk+r_ax(i)
                            G1ijk=G1ijk*octopole(i,j,k,i_ex)

                            G1ijkl=0.0_r8_kind
                            if(i==j .and. k==l) G1ijkl=G1ijkl+1.0_r8_kind
                            if(i==k .and. j==l) G1ijkl=G1ijkl+1.0_r8_kind
                            if(j==k .and. i==l) G1ijkl=G1ijkl+1.0_r8_kind
                            G1ijkl=G1ijkl*octopole(i,j,k,i_ex)

                            gradient(l)=gradient(l)+Z_a*( &
                                 -35_r8_kind*r_ax(l)*Gijk/d_ax**9+ &
                                 5_r8_kind*Gijkl/d_ax**7+ &
                                 5_r8_kind*r_ax(l)*G1ijk/d_ax**7- &
                                 G1ijkl/d_ax**5)/5_r8_kind

#ifdef WITH_EFP
                            if(l==i) Orr(l)=Orr(l)+octopole(i,j,k,i_ex)*r_ax(j)*r_ax(k)
                            if(l==j) Orr(l)=Orr(l)+octopole(i,j,k,i_ex)*r_ax(i)*r_ax(k)
                            if(l==k) Orr(l)=Orr(l)+octopole(i,j,k,i_ex)*r_ax(i)*r_ax(j)
#endif
                         end do
                      end do
                   end do
                end do

#ifdef WITH_EFP
                EE=Z_a*r_ax/d_ax**7
                torque=torque+vector_product(Orr,EE)
#endif
             end do
          end do

          do i=1,grad_dim
             gradient_oct_totalsym(index+i)=gradient_oct_totalsym(index+i)- &
                  sum(rotmat(i,:)*gradient(:))

#ifdef WITH_EFP
             torque_oct_totalsym(index+i)=torque_oct_totalsym(index+i)+ &
                  sum(rotmat(i,:)*torque(:))
#endif
          enddo

       end do
    end do

  end subroutine calc_X_grads
  !*************************************************************

  !*************************************************************
  subroutine transform_X_grad_to_cart
    ! purpose: transform symmetry adapted gradient components to cartesian
    !          coordinates and add them to an array of cartesian gradients
    !
    implicit none
    ! *** end of interface ***

    integer (i4_kind) :: i_unique, i_equal, index, i, grad_dim
    real(kind=r8_kind),pointer :: rotmat(:,:)
#ifdef WITH_EFP
    integer (i4_kind) :: i_group
    real (r8_kind) :: r12(3), r1(3), r2(3), f(3)
#endif

    if(moving_X_centers) then
    do i_unique=1,N_pd
       index=dipoles_grad_index(i_unique)
       grad_dim=dipoles_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,pd_array(i_unique)%n_equal_dipoles
          rotmat=>dipoles_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             gradient_dip_cartesian(i_unique)%m(:,i_equal) = &
                  gradient_dip_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*gradient_dip_totalsym(index+i)
#ifdef WITH_EFP
             if(efp) then
                torque_dip_cartesian(i_unique)%m(:,i_equal) = &
                     torque_dip_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*torque_dip_totalsym(index+i)
             end if
#endif
          end do
       end do
    enddo
#ifdef WITH_EFP
    if(efp) then
       do i_unique=1,N_pd
          do i_equal=1,pd_array(i_unique)%n_equal_dipoles
             r1=pd_array(i_unique)%position(:,i_equal)
             i_group=pd_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-gradient_dip_cartesian(i_unique)%m(:,i_equal)

             torque_dip_cartesian(i_unique)%m(:,i_equal) = &
                  torque_dip_cartesian(i_unique)%m(:,i_equal) + vector_product(r12,f)
          end do
       end do
    end if
#endif

    do i_unique=1,N_pq
       index=qpoles_grad_index(i_unique)
       grad_dim=qpoles_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,pq_array(i_unique)%n_equal_qpoles
          rotmat=>qpoles_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             gradient_quad_cartesian(i_unique)%m(:,i_equal) = &
                  gradient_quad_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*gradient_quad_totalsym(index+i)
#ifdef WITH_EFP
             if(efp) then
                torque_quad_cartesian(i_unique)%m(:,i_equal) = &
                     torque_quad_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*torque_quad_totalsym(index+i)
             end if
#endif
          end do
       end do
    enddo
#ifdef WITH_EFP
    if(efp) then
       do i_unique=1,N_pq
          do i_equal=1,pq_array(i_unique)%n_equal_qpoles
             r1=pq_array(i_unique)%position(:,i_equal)
             i_group=pq_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-gradient_quad_cartesian(i_unique)%m(:,i_equal)

             torque_quad_cartesian(i_unique)%m(:,i_equal) = &
                  torque_quad_cartesian(i_unique)%m(:,i_equal) + vector_product(r12,f)
          end do
       end do
    end if
#endif

    do i_unique=1,N_po
       index=opoles_grad_index(i_unique)
       grad_dim=opoles_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,po_array(i_unique)%n_equal_opoles
          rotmat=>opoles_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             gradient_oct_cartesian(i_unique)%m(:,i_equal) = &
                  gradient_oct_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*gradient_oct_totalsym(index+i)
#ifdef WITH_EFP
             if(efp) then
                torque_oct_cartesian(i_unique)%m(:,i_equal) = &
                     torque_oct_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*torque_oct_totalsym(index+i)
             end if
#endif
          end do
       end do
    enddo
#ifdef WITH_EFP
    if(efp) then
       do i_unique=1,N_po
          do i_equal=1,po_array(i_unique)%n_equal_opoles
             r1=po_array(i_unique)%position(:,i_equal)
             i_group=po_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-gradient_oct_cartesian(i_unique)%m(:,i_equal)

             torque_oct_cartesian(i_unique)%m(:,i_equal) = &
                  torque_oct_cartesian(i_unique)%m(:,i_equal) + vector_product(r12,f)
          end do
       end do
    end if
#endif
    end if

    if(moving_R_centers) then
    do i_unique=1,N_rc
       index=rep_grad_index(i_unique)
       grad_dim=rep_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,rc_array(i_unique)%n_equal_centers
          rotmat=>rep_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             gradient_rep_cartesian(i_unique)%m(:,i_equal) = &
                  gradient_rep_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*gradient_rep_totalsym(index+i)
          end do
       end do
    enddo
#ifdef WITH_EFP
    if(efp) then
       do i_unique=1,N_rc
          do i_equal=1,rc_array(i_unique)%n_equal_centers
             r1=rc_array(i_unique)%position(:,i_equal)
             i_group=rc_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-gradient_rep_cartesian(i_unique)%m(:,i_equal)

             torque_rep_cartesian(i_unique)%m(:,i_equal) = &
                  vector_product(r12,f)
          end do
       end do
    end if
#endif
    end if

  end subroutine transform_X_grad_to_cart
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
  subroutine X_grad_cart_write()
    !  Purpose: writing the cartesian gradients
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
     use iounitadmin_module, only: output_unit
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j
    !------------ Executable code ------------------------------------

    if(print_X_grad) then
    if(N_pd > 0) then
#ifdef WITH_EFP
       if(efp) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Point Dipoles'
          do i=1,N_pd
             write(output_unit,*) 'Unique Point Dipole:',i
             do j=1,pd_array(i)%n_equal_dipoles
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     gradient_dip_cartesian(i)%m(:,j),' / ',torque_dip_cartesian(i)%m(:,j)
             enddo
          end do
       else
#endif
          write(output_unit,'(/A)') 'Cartesian gradients on Point Dipoles'
          do i=1,N_pd
             write(output_unit,*) 'Unique Point Dipole:',i
             do j=1,pd_array(i)%n_equal_dipoles
                write(output_unit,'(A14,3F15.10)') 'Equal Center: ',&
                     gradient_dip_cartesian(i)%m(:,j)
             enddo
          end do
#ifdef WITH_EFP
       end if
#endif
    end if

    if(N_pq > 0) then
#ifdef WITH_EFP
       if(efp) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Point Quadrupoles'
          do i=1,N_pq
             write(output_unit,*) 'Unique Point Quadrupole:',i
             do j=1,pq_array(i)%n_equal_qpoles
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     gradient_quad_cartesian(i)%m(:,j),' / ',torque_quad_cartesian(i)%m(:,j)
             enddo
          end do
#endif
       else
          write(output_unit,'(/A)') 'Cartesian gradients on Point Quadrupoles'
          do i=1,N_pq
             write(output_unit,*) 'Unique Point Quadrupole:',i
             do j=1,pq_array(i)%n_equal_qpoles
                write(output_unit,'(A14,3F15.10)') 'Equal Center:',&
                     gradient_quad_cartesian(i)%m(:,j)
             enddo
          end do
#ifdef WITH_EFP
       end if
#endif
    end if

    if(N_po > 0) then
#ifdef WITH_EFP
       if(efp) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Point Octopoles'
          do i=1,N_po
             write(output_unit,*) 'Unique Point Octopole:',i
             do j=1,po_array(i)%n_equal_opoles
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     gradient_oct_cartesian(i)%m(:,j),' / ',torque_oct_cartesian(i)%m(:,j)
             enddo
          end do
       else
#endif
          write(output_unit,'(/A)') 'Cartesian gradients on Point Octopoles'
          do i=1,N_po
             write(output_unit,*) 'Unique Point Octopole:',i
             do j=1,po_array(i)%n_equal_opoles
                write(output_unit,'(A14,3F15.10)') 'Equal Center:',&
                     gradient_oct_cartesian(i)%m(:,j)
             enddo
          end do
#ifdef WITH_EFP
       end if
#endif
    end if
    end if

    if(print_R_grad) then
    if(N_rc > 0) then
#ifdef WITH_EFP
       if(efp) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Repulsion Centers (QM-EFP)'
          do i=1,N_rc
             write(output_unit,*) 'Unique Repulsion Center:',i
             do j=1,rc_array(i)%n_equal_centers
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     gradient_rep_cartesian(i)%m(:,j),' / ',torque_rep_cartesian(i)%m(:,j)
             enddo
          end do
       else
#endif
          write(output_unit,'(/A)') 'Cartesian gradients on Repulsion Centers (QM-EFP)'
          do i=1,N_rc
             write(output_unit,*) 'Unique Repulsion Center:',i
             do j=1,rc_array(i)%n_equal_centers
                write(output_unit,'(A14,3F15.10)') 'Equal Center: ',&
                     gradient_rep_cartesian(i)%m(:,j)
             enddo
          end do
#ifdef WITH_EFP
       end if
#endif
    end if
    end if

  end subroutine X_grad_cart_write
  !*************************************************************
  !--------------- End of module -------------------------------------
end module point_dqo_module

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
module efp_rep_module
  !---------------------------------------------------------------
  !
  !  Purpose: Calculate interfragment interactions
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: AS
  !  Date: 11.07
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
  use datatype, only: arrmat2, arrmat3
  use unique_atom_module, only: unique_atom_type
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------
  type, public :: repulsive_f_type
     character(len=12) :: name
     integer(i4_kind) :: N_equal_centers
     real(r8_kind), pointer :: position(:,:)
     real(r8_kind) :: C(4)
     real(r8_kind) :: a(4)
     integer(i4_kind), pointer :: group(:)
     integer(i4_kind), pointer :: type(:)
  end type repulsive_f_type

  !------------ Declaration of constants and variables ------------
  integer(i4_kind), public :: N_rcf=0_i4_kind
  type(repulsive_f_type), pointer, public :: rcf_array(:)
  type(arrmat3), pointer, public :: repf_grad_info(:)
  integer(i4_kind), allocatable, public :: repf_grad_index(:)
  integer(i4_kind), public :: totsym_grad_repf_length=0
  real(r8_kind), public,allocatable,target :: gradient_repf_totalsym(:)
  type(arrmat2), public,allocatable :: gradient_repf_cartesian(:)
  type(arrmat2), public,allocatable :: torque_repf_cartesian(:)

  !------------ Interface statements ------------------------------
  !------------ public functions and subroutines ------------------
  public symm_repf_points, dealloc_repf, init_repf_grads, repf_grads_shutdown
  public dealloc_repf_inform, transform_repf_grad_to_cart, repf_grad_cart_write

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------
  !------------ Declaration of constants and variables ------------
  type(unique_atom_type), pointer :: unique_external_points(:)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains
  subroutine symm_repf_points (N_total, P, ef_grp, type, name, C, A)
    ! Symmetrize repulsive centers
    !** End of interface *****************************************
    !------------ Modules used -----------------------------------
    use operations_module, only: operations_gradients
    use point_dqo_module, only: moving_R_centers
    use group_module, only: symm_transformation_int, group_num_el, &
         sub_group, group_coset, ylm_trafos, group_coset_decomp
    use symm_module, only: symm_adapt_centers
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind) :: N_total
    real(r8_kind) :: P(3,N_total)
    integer(i4_kind) :: ef_grp(N_total)
    integer(i4_kind) :: type(N_total)
    character(len=12)  :: name(N_total)
    real(r8_kind) :: C(4,N_total)
    real(r8_kind) :: A(4,N_total)
    !------------ Declaration of local constants  ----------------
    real(r8_kind),parameter :: small = 1.e-10_r8_kind
    ! very small value
    !------------ Declaration of local variables -----------------
    integer (i4_kind) :: i, j, k, n, l, alloc_stat, n_size
    ! counters
    integer (i4_kind) :: n_equal, n_equal_check, N_eq
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
    type(repulsive_f_type), allocatable :: buffer_r(:)
    logical :: do_alloc
    !------------ Executable code --------------------------------

    allocate(help_dim(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    help_dim=.false.

    allocate(point_trafos_buf(N_total),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    allocate(buffer_r(N_total),stat=alloc_stat)
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

       allocate(buffer_r(n_size)%position(3,n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_r(n_size)%group(n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       allocate(buffer_r(n_size)%type(n_equal),stat=alloc_stat)
       ASSERT(alloc_stat==0)

       n_equal_check=0
       j_n_equal: do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_ec%elements(1,j)),position)
          position3(1) = position2(1)
          position3(2) = position2(3)
          position3(3) = position2(2)
          k_N_tot: do k=1,N_total
             if(.not.help_dim(k)) then
                if (sqrt(dot_product(position3-P(:,k),position3-P(:,k))) <= small) then
                   help_dim(k)=.true.
                   n_equal_check=n_equal_check+1
                   buffer_r(n_size)%position(:,j)=P(:,k)
                   buffer_r(n_size)%group(j)=ef_grp(k)
                   buffer_r(n_size)%type(j)=type(k)
                   buffer_r(n_size)%n_equal_centers=n_equal
                   buffer_r(n_size)%name=trim(name(k))
                   buffer_r(n_size)%C=C(:,k)
                   buffer_r(n_size)%A=A(:,k)
                end if
             endif
          end do k_N_tot
       end do j_n_equal
       ASSERT(n_equal_check==n_equal)

       deallocate(local_groups_ec%elements, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do i_N_tot

    N_rcf=n_size
    do_alloc=.false.
    if(.not. associated(rcf_array)) then
       allocate(rcf_array(N_rcf),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do_alloc=.true.
    end if
    do i=1,N_rcf
       rcf_array(i)%name=buffer_r(i)%name
       rcf_array(i)%N_equal_centers=buffer_r(i)%n_equal_centers
       N_eq=rcf_array(i)%N_equal_centers
       if(do_alloc) then
          allocate(rcf_array(i)%position(3,N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(rcf_array(i)%group(N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(rcf_array(i)%type(N_eq),stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end if
       rcf_array(i)%position=buffer_r(i)%position
       rcf_array(i)%group=buffer_r(i)%group
       rcf_array(i)%type=buffer_r(i)%type
       rcf_array(i)%C=buffer_r(i)%C
       rcf_array(i)%A=buffer_r(i)%A

       deallocate(buffer_r(i)%position,stat=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(buffer_r(i)%group,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do
    deallocate(buffer_r,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    deallocate(help_dim,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    if(operations_gradients .and. moving_R_centers) then
       allocate(point_trafos_ec(n_size),stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          n_equal=rcf_array(i)%n_equal_centers
          allocate(point_trafos_ec(i)%matrix(n_equal,n_equal,group_num_el),stat=alloc_stat)
          ASSERT(alloc_stat==0)

          point_trafos_ec(i)%matrix=point_trafos_buf(i)%matrix

          deallocate(point_trafos_buf(i)%matrix,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
    end if
    deallocate(point_trafos_buf, stat=alloc_stat)
    ASSERT(alloc_stat==0)

    if(operations_gradients .and. moving_R_centers) then
       allocate(unique_external_points(n_size), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_size
          unique_external_points(i)%N_equal_atoms=rcf_array(i)%N_equal_centers
          allocate(unique_external_points(i)%position(3,rcf_array(i)%N_equal_centers), &
               stat=alloc_stat)
          ASSERT(alloc_stat==0)

          unique_external_points(i)%position=rcf_array(i)%position
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

       call repf_points_grad_inform()
    end if

  end subroutine symm_repf_points
  !*************************************************************

  !*************************************************************
  subroutine repf_points_grad_inform()
    use unique_atom_module, only: unique_atom_type, unique_atom_partner_type, &
         unique_atom_symadapt_type
    use symmetry_data_module, only: get_totalsymmetric_irrep,find_totalsymmetric_irrep
    use uatom_symmadapt, only: sa_free
    !** End of interface *****************************************

    !----------- declaration of local variables -------------
    type(unique_atom_type)         , pointer :: ua
    type(unique_atom_partner_type) , pointer :: sap
    type(unique_atom_symadapt_type), pointer :: sa
    type(arrmat3)                  , pointer :: gi
    integer (i4_kind) :: i, j ,k
    integer(i4_kind) :: i_ua,i_if,i_cd,i_ea,ts, counter, &
         n_indep,n_equals,alloc_stat
    logical :: do_alloc
    !----------- executable code -----------------------------

    if(N_rcf > 0) then
       do_alloc=.false.
       if(.not. associated(repf_grad_info)) then
          allocate(repf_grad_info(N_rcf),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          do_alloc=.true.
       end if

       call find_totalsymmetric_irrep()
       ts = get_totalsymmetric_irrep()

       do i_ua=1,N_rcf
          ua => unique_external_points(i_ua)
          sap => ua%symadapt_partner(ts,1) ! l = 1
          gi => repf_grad_info(i_ua)
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
          allocate(repf_grad_index(N_rcf+1),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end if

       counter = 1
       do i_ua=1,N_rcf
          ua => unique_external_points(i_ua)
          repf_grad_index(i_ua)=counter
          counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
       enddo
       repf_grad_index(N_rcf+1)=counter

       totsym_grad_repf_length=counter-1

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
    end if

  end subroutine repf_points_grad_inform
  !*************************************************************

  !*************************************************************
  subroutine dealloc_repf
    !  Purpose: deallocate allocated repalsive centers
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,status
    !------------ Executable code --------------------------------

    if(associated(rcf_array)) then

       do i=1,N_rcf
          deallocate(rcf_array(i)%position,stat=status)
          ASSERT(status==0)
          deallocate(rcf_array(i)%group,stat=status)
          ASSERT(status==0)
          deallocate(rcf_array(i)%type,stat=status)
          ASSERT(status==0)
       end do

       deallocate(rcf_array,stat=status)
       ASSERT(status==0)
    end if

  end subroutine dealloc_repf
  !*************************************************************

  !*************************************************************
  subroutine dealloc_repf_inform()
    ! Purpose: deallocate the variables 'repf_grad_info'
    !** End of interface *****************************************
    ! ----------- declaration of local variables -------------
    integer(kind=i4_kind) :: i,alloc_stat
    ! ----------- executable code -----------------------------

    if(associated(repf_grad_info)) then
       do i=1,size(repf_grad_info)
          deallocate(repf_grad_info(i)%m,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
       deallocate(repf_grad_info,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(repf_grad_info)

       deallocate(repf_grad_index,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

  end subroutine dealloc_repf_inform
  !*************************************************************

  !*************************************************************
  subroutine init_repf_grads()
    ! Purpose: allocate arrays to calculate garadients at repalsive centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use point_dqo_module, only: moving_R_centers
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code --------------------------------

    if(N_rcf > 0) then
       allocate(gradient_repf_cartesian(N_rcf),stat=status)
       ASSERT(status==0)
       allocate(torque_repf_cartesian(N_rcf),stat=status)
       ASSERT(status==0)

       do i=1,N_rcf
          n_eq=rcf_array(i)%N_equal_centers
          allocate(gradient_repf_cartesian(i)%m(3,n_eq),stat=status)
          ASSERT(status==0)
          gradient_repf_cartesian(i)%m=0.0_r8_kind
          allocate(torque_repf_cartesian(i)%m(3,n_eq),stat=status)
          ASSERT(status==0)
          torque_repf_cartesian(i)%m=0.0_r8_kind
       end do

       allocate(gradient_repf_totalsym(totsym_grad_repf_length),stat=status)
       ASSERT(status==0)
       gradient_repf_totalsym=0.0_r8_kind
    end if

  end subroutine init_repf_grads
  !*************************************************************

  !*************************************************************
  subroutine repf_grads_shutdown()
    !  Purpose: free arrays for calculating gradients at repalsive centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use point_dqo_module, only: moving_R_centers
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i
    !------------ Executable code --------------------------------

    if(moving_R_centers) then
       if(N_rcf > 0) then
          do i=1,N_rcf
             deallocate(gradient_repf_cartesian(i)%m,stat=status)
             ASSERT(status==0)
             deallocate(torque_repf_cartesian(i)%m,stat=status)
             ASSERT(status==0)
          end do
          deallocate(gradient_repf_cartesian,stat=status)
          ASSERT(status==0)
          deallocate(torque_repf_cartesian,stat=status)
          ASSERT(status==0)

          deallocate(gradient_repf_totalsym,stat=status)
          ASSERT(status==0)
       end if
    end if

  end subroutine repf_grads_shutdown
  !*************************************************************

  !*************************************************************
  subroutine transform_repf_grad_to_cart
    ! purpose: transform symmetry adapted gradient components to cartesian
    !          coordinates and add them to an array of cartesian gradients
    !
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use point_dqo_module, only: moving_R_centers
    use pointcharge_module, only: rcm
    !------------ Declaration of local variables ----------------
    integer (i4_kind) :: i_unique, i_equal, index, i, grad_dim, i_group
    real(kind=r8_kind),pointer :: rotmat(:,:)
    real(kind=r8_kind) :: r12(3),r1(3),r2(3),f(3)
    !------------ Executable code -------------------------------

    if(moving_R_centers) then
       do i_unique=1,N_rcf
          index=repf_grad_index(i_unique)
          grad_dim=repf_grad_index(i_unique+1)-index
          index=index-1
          do i_equal=1,rcf_array(i_unique)%n_equal_centers
             rotmat=>repf_grad_info(i_unique)%m(:,:,i_equal)
             do i=1,grad_dim
                gradient_repf_cartesian(i_unique)%m(:,i_equal) = &
                     gradient_repf_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*gradient_repf_totalsym(index+i)
             end do
          end do
       enddo

       do i_unique=1,N_rcf
          do i_equal=1,rcf_array(i_unique)%n_equal_centers
             r1=rcf_array(i_unique)%position(:,i_equal)
             i_group=rcf_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-gradient_repf_cartesian(i_unique)%m(:,i_equal)

             torque_repf_cartesian(i_unique)%m(:,i_equal) = &
                  vector_product(r12,f)
          end do
       end do
    end if

  end subroutine transform_repf_grad_to_cart
  !*************************************************************

  !*************************************************************
  function vector_product(v1,v2)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: vector_product(3)
    real(r8_kind) :: v1(3),v2(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !*************************************************************

  !*************************************************************
  subroutine repf_grad_cart_write()
    !  Purpose: writing the cartesian gradients
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use iounitadmin_module, only: output_unit
    use point_dqo_module, only: print_R_grad
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,j
    !------------ Executable code --------------------------------

    if(print_R_grad) then
       if(N_rcf > 0) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Repulsive Centers (EFP-EFP)'
          do i=1,N_rcf
             write(output_unit,*) 'Unique Repulsive Center:',i
             do j=1,rcf_array(i)%n_equal_centers
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     gradient_repf_cartesian(i)%m(:,j),' / ',torque_repf_cartesian(i)%m(:,j)
             enddo
          end do
       end if
    end if

  end subroutine repf_grad_cart_write
  !*************************************************************

  !--------------- End of module ----------------------------------
end module efp_rep_module

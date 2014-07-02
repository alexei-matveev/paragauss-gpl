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
  !===================================================================
! Public interface of module
  !===================================================================
module elec_static_field_module
!
!  Calculate electrostatic field
!
!  Author: AS
!  Date: 03/00
!
!== Interrupt of public interface of module ====
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: MF
! Date:   7/2000
! Description: adaption to VPP
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

# include "def.h"
  use type_module          ! contains standard data types
  use datatype    ! user defined types
  use options_module
  use msgtag_module
  use comm_module
  use unique_atom_module, only: unique_atom_type
  use symmetry_data_module  ! provide ssym

  implicit none

  private
  save

!------------ public functions and subroutines ------------------
  public send_surf_point,receive_surf_point, dealloc_surf_points, surf_points_grad_information, &
       surf_points_gradinfo_dealloc,start_read_field_e,get_field_nuc,get_field_pc, &
       read_field_e,send_receive_field,destroy_field_file,deallocate_field,bounds_calc_field, &
       bounds_free_field, fill_surf_points, &
       transform_to_cart_field, field_integral_open, field_integral_close
  public receive_surf_point1,surf_points_grad_information1

  type, public :: outward_normal
     integer(kind=i4_kind)                         :: N_equal_points
     ! Number of partners of unique point charge
     real(kind=r8_kind), pointer                :: position(:,:)
     ! positions of first partner ("equal point")
     real(kind=r8_kind)                :: out_normal(3)
     real(kind=r8_kind)                :: area
  end type outward_normal

  ! FIXME: %lower1 and %upper1 are obscure and may be removed,
  !        see bounds_module.f90. After that this struct is just
  !        a plain array.
  type, public ::  field_bounds
     integer(kind=i4_kind)                :: lower1,upper1
     ! start bounds
     integer(kind=i4_kind),pointer        :: item_arr(:)
  end type field_bounds

  type(outward_normal), allocatable, target, public :: surface_points(:) ! to calculate E
  integer(kind=i4_kind), public :: N_surface_points  ! to calculate E
  type(unique_atom_type), pointer, public :: unique_surf_points(:)
  type(arrmat3),pointer,public :: surf_points_grad_info(:)
  integer(kind=i4_kind),public,allocatable :: surf_points_grad_index(:)
  logical,public :: calc_normal = .false.

  type(arrmat2), allocatable, public :: E_ele(:)
  type(arrmat2), allocatable, public :: E_nuc(:)
  type(arrmat2), allocatable, public :: E_pc(:)
#ifdef WITH_EFP
  type(arrmat2), allocatable, public :: E_mp(:) !field coming from different multipoles
  type(arrmat2), allocatable, public :: E_id(:) !field coming from induced dipoles
  type(arrmat2), allocatable, public :: E_id1(:) !field coming from induced dipoles
  type(arrmat2), allocatable, public :: E_cav(:) !field coming from cavity surface charges
                                                 !produced by efp multipoles and induced dipoles
  type(arrmat2), allocatable, public :: E_cav1(:) !field coming from cavity surface charges
                                                 !produced by efp multipoles and induced dipoles
#endif
  real(r8_kind),allocatable, public :: totalsym_field(:)
  integer(i4_kind), public :: totsym_field_length

  real(kind=r8_kind),allocatable, public:: E_n(:) ! projection of electrostatic field on
                                                  ! normal of solvation cavity surface
                                                  ! only electron part
  real(kind=r8_kind),allocatable :: E_n_tmp(:)

  type(field_bounds), public      :: bounds
  !===================================================================
  ! End of public interface of module
  !===================================================================
!To calculate electrostatic field:
!1. Define points where electrostatic field has to be calculated P(3,N)
!2. calc_normal = .false.
!3. call fill_surf_points(P,N)
!4. call surf_points_grad_information()
!5. call field_calculate() - integrals
!6. call main_scf() or use saved density matrix (main_scf can be called
!   in any place before calculation of electrostatic field)
!7. call start_read_field_e()
!8. call get_field_nuc()
!9. call get_field_pc()
!10. call surf_points_gradinfo_dealloc()
!11.Here you can use E_ele, E_nuc, E_pc
!12.call deallocate_field()
!13.call dealloc_surf_points()
!14.call bounds_free_field()
!15.call post_dens_master()
!16.If you want to keep integrals in memory, do not forget to deallocate them in the end
!           if ( .not. options_integrals_on_file() ) then !!!!!!!!!!!!!!!!!
!              call integralstore_deallocate_pcm()
!              if ( comm_parallel() ) then
!                 call comm_init_send(comm_all_other_hosts,msgtag_intstore_dealloc)
!                 call comm_send()
!              endif
!            end if

contains

  !******************************************************
  subroutine fill_surf_points(point_array,N_total)
    !------------ Modules used --------------------------------------
    use group_module, only : ylm_trafos,sub_group,group_coset, &
         symm_transformation_int,group_num_el,group_coset_decomp
    use symm_module, only : symm_adapt_centers
    !== End of interface ==========================================
    real(r8_kind)    :: point_array(3,N_total)
    integer(i4_kind) :: N_total
    !------------ Declaration of local variables ------------------
    integer(kind=i4_kind) :: n_equal,n_equal_check,n_size
    real(kind=r8_kind)    :: position(3),position2(3),position3(3)
    type(sub_group) :: local_groups
    type(group_coset) :: cosets
    type(symm_transformation_int),allocatable :: point_trafos_buf(:)
    type(symm_transformation_int), pointer :: point_trafos(:)
    logical, allocatable :: help_dim(:)
    type(outward_normal),allocatable :: buffer(:)

    real(kind=r8_kind),parameter :: small = 1.0e-10_r8_kind
    integer(kind=i4_kind) :: i,n,j,k,l,status
    logical :: do_alloc
    !------------ Executable code ---------------------------------

    allocate(help_dim(N_total),stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: allocation of help_dim is failed")
    help_dim=.false.

    allocate(buffer(N_total),point_trafos_buf(N_total),stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: allocation of buffer is failed")

    n_size=0
    l=0
    i_N_tot: do i=1,N_total
       if(help_dim(i)) cycle i_N_tot
       n_size=n_size+1
       ! reorder coordinates of points as
       ! (x,y,z) --> (x,z,y) in order to comply with the
       ! convention for angular momentum l=1
       position(1) = point_array(1,i)
       position(2) = point_array(3,i)
       position(3) = point_array(2,i)
       !
       ! determine local symmetry groups
       !
       ! now apply all symmetry operations to the position of the
       ! surface point
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
          endif
       enddo

       ! allocate group elements
       local_groups%num_el = n
       allocate(local_groups%elements(n))
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: allocation LOCAL_GROUPS  is failed")

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups%elements(n) = j
          end if
       enddo
       !
       ! now determine symmetry equivalent atoms
       !
       call group_coset_decomp(n_equal,local_groups,&
            cosets,point_trafos_buf(n_size)%matrix)
       !
       ! search of positions of equal atoms
       !
       allocate(buffer(n_size)%position(3,n_equal), stat=status)
       if (status .ne. 0) call error_handler( &
            "elec_static_field_module: allocation of BUFFER%POSITION  failed")

       n_equal_check=0
       j_n_equal: do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets%elements(1,j)),position)
          position3(1) = position2(1)
          position3(2) = position2(3)
          position3(3) = position2(2)
          k_N_tot: do k=1,N_total
             if(.not.help_dim(k)) then
                if (sqrt(dot_product(position3-point_array(:,k), &
                     position3-point_array(:,k))) <= small) then

                   help_dim(k)=.true.
                   n_equal_check=n_equal_check+1
                   buffer(n_size)%position(:,j)=point_array(:,k)
                   buffer(n_size)%n_equal_points=n_equal
                endif
             endif
          end do k_N_tot
       end do j_n_equal
       if (n_equal_check /= n_equal) call error_handler( &
            "elec_static_field_module: The point array to calculate electrostatic &
            & field is not corresponded to the chosen symmetry")

       deallocate(local_groups%elements,stat=status)
       if (status .ne. 0 ) call error_handler( &
            "elec_static_field_module: deallocation of points HELPERS is failed")
    end do i_N_tot

    do_alloc=.false.
    N_surface_points=n_size
    if(.not. allocated(surface_points)) then
       allocate(surface_points(N_surface_points),stat=status)
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: allocation SURFACE_POINTS is  failed")
       do_alloc=.true.
    end if
    do i=1,n_size
       surface_points(i)%n_equal_points=buffer(i)%n_equal_points

       if(do_alloc) then
          allocate(surface_points(i)%position(3,buffer(i)%n_equal_points), stat=status)
          if (status .ne. 0 ) call error_handler( &
               "elec_static_field_module: allocation of SURFACE_POINTS(i)%position is failed")
       end if

       surface_points(i)%position=buffer(i)%position

       deallocate(buffer(i)%position, stat=status)
       if (status .ne. 0 ) call error_handler( &
            "elec_static_field_module: deallocation of BUFFER(i)%position is failed")
    enddo
    deallocate(buffer,stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: deallocation of BUFFER is failed")
    deallocate(help_dim,stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: deallocation of HELP_DIM  is failed")

    allocate(point_trafos(n_size),  stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: allocation of point_trafos is failed")
    do i=1,n_size
       n_equal=surface_points(i)%n_equal_points
       allocate(point_trafos(i)%matrix(n_equal,n_equal,group_num_el), stat=status)
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: allocation of point_trafos%matix is failed")

       point_trafos(i)%matrix=point_trafos_buf(i)%matrix

       deallocate(point_trafos_buf(i)%matrix, stat=status)
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: deallocation of point_trafos_buf%matrix is failed")
    enddo

    deallocate(point_trafos_buf, stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: deallocation of point_trafos_buf is failed")

    allocate(unique_surf_points(n_size),  stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: allocation of unique_surf_charges failed")

    do i=1,n_size
       unique_surf_points(i)%N_equal_atoms=surface_points(i)%n_equal_points
       allocate(unique_surf_points(i)%position(3,surface_points(i)%n_equal_points), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: allocation of unique_surf_points%position failed")
       do j=1,surface_points(i)%n_equal_points
          unique_surf_points(i)%position(:,j)=surface_points(i)%position(:,j)
          unique_surf_points(i)%lmax_all=1
       enddo
    enddo

    call symm_adapt_centers(unique_surf_points, point_trafos, L_MAX=1)

    do i=1,n_size
       deallocate(point_trafos(i)%matrix, stat=status)
       if ( status /= 0) call error_handler( &
            "elec_static_field_module: deallocation point_trafos%matix is failed")
    enddo
    deallocate(point_trafos, stat=status)
    if ( status /= 0) call error_handler( &
         "elec_static_field_module: deallocation point_trafos is failed")

    if (comm_parallel()) call send_surf_point()

  end subroutine fill_surf_points
  !********************** ********************************

  !******************************************************
  subroutine send_surf_point
   !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, length, status
    type(outward_normal), pointer :: ps

    call comm_init_send(comm_all_other_hosts,msgtag_surf_point)

    call commpack(calc_normal,status)
    call commpack(N_surface_points,status)
    do i=1,N_surface_points
       ps => surface_points(i)
       call commpack(ps%N_equal_points,status)
       length = ps%N_equal_points*3
       call commpack(ps%position(1,1),length , 1, status)
       call commpack(ps%out_normal(1), 3, 1, status)
       call commpack(ps%area,status)
    enddo

    call comm_send()

  end subroutine send_surf_point
  !******************************************************

  !******************************************************
  subroutine receive_surf_point
   !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, length, status
    logical :: do_alloc

    call communpack(calc_normal,status)
    call communpack(N_surface_points,status)
    do_alloc=.false.
    if(.not. allocated(surface_points)) then
       allocate(surface_points(N_surface_points),stat=status)
       if ( status .ne. 0) call error_handler( &
            "elec_ststic_field_module:receive_surf_point: allocated surface_points failed" )
       do_alloc=.true.
    end if

    do i=1,N_surface_points
       call communpack(surface_points(i)%N_equal_points,status)

       if(do_alloc) then
          allocate(surface_points(i)%position(3,surface_points(i)%N_equal_points),stat=status)
          if ( status .ne. 0) call error_handler( &
               "elec_ststic_field_module:receive_surf_point: allocated surface_points%position failed" )
       end if

       length = surface_points(i)%N_equal_points*3
       call communpack(surface_points(i)%position(1,1),length , 1, status)

       call communpack(surface_points(i)%out_normal(1), 3, 1, status)
       call communpack(surface_points(i)%area,status)
    enddo

  end subroutine receive_surf_point
  !******************************************************

  !******************************************************
  subroutine receive_surf_point1
   !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, length, status
    logical :: do_alloc

    call comm_save_recv(comm_master_host,msgtag_surf_point)

    call communpack(calc_normal,status)
    call communpack(N_surface_points,status)
    do_alloc=.false.
    if(.not. allocated(surface_points)) then
       allocate(surface_points(N_surface_points),stat=status)
       ASSERT(status==0)
       do_alloc=.true.
    end if

    do i=1,N_surface_points
       call communpack(surface_points(i)%N_equal_points,status)

       if(do_alloc) then
          allocate(surface_points(i)%position(3,surface_points(i)%N_equal_points),stat=status)
          ASSERT(status==0)
       end if

       length = surface_points(i)%N_equal_points*3
       call communpack(surface_points(i)%position(1,1),length , 1, status)

       call communpack(surface_points(i)%out_normal(1), 3, 1, status)
       call communpack(surface_points(i)%area,status)
    enddo

  end subroutine receive_surf_point1
  !******************************************************

  !******************************************************
  subroutine dealloc_surf_points
   !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, status

    do i=1,N_surface_points
       deallocate(surface_points(i)%position, stat=status)
       if ( status .ne. 0) call error_handler( &
         "elec_ststic_field_module:dealloc_surf_points: deallocation of surface_points%position are failed" )
    enddo

    deallocate(surface_points,stat=status)
    if ( status /= 0) call error_handler( &
            "elec_ststic_field_module: deallocation of unique_surf_points are failed")

  end subroutine dealloc_surf_points
  !******************************************************

  !******************************************************
  subroutine surf_points_grad_information
    !** End of interface *****************************************
    use unique_atom_module, only: unique_atom_type, unique_atom_partner_type, &
         unique_atom_symadapt_type
    use symmetry_data_module, only: get_totalsymmetric_irrep,find_totalsymmetric_irrep
    use uatom_symmadapt, only: sa_free
    ! ----------- declaration of local variables -------------
    type(unique_atom_type)         , pointer :: ua
    type(unique_atom_partner_type) , pointer :: sap
    type(unique_atom_symadapt_type), pointer :: sa
    type(arrmat3)                  , pointer :: gi
    integer(kind=i4_kind) :: i_ua,i_if,i_cd,i_ea,ts,counter,i,j,k, &
         n_indep,n_equals,alloc_stat
    logical :: do_alloc

    external error_handler
    ! ----------- executable code -----------------------------

    if(comm_i_am_master()) then
       do_alloc=.false.
       if(.not.associated(surf_points_grad_info)) then
          allocate( surf_points_grad_info(N_surface_points),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation (1) failed")
          do_alloc=.true.
       end if

       call find_totalsymmetric_irrep()
       ts = get_totalsymmetric_irrep()

       do i_ua=1,N_surface_points
          ua => unique_surf_points(i_ua)
          sap => ua%symadapt_partner(ts,1) ! l = 1
          gi => surf_points_grad_info(i_ua)
          n_indep = sap%N_independent_fcts
          n_equals = ua%N_equal_atoms
          if(do_alloc) then
             allocate ( gi%m(n_indep,3,n_equals),STAT=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("surf_point_grad_information: allocation (2) failed")
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
                   call error_handler("surf_point_grad_information: sth. wrong")
                end select

             enddo
          enddo
       end do

       if(do_alloc) then
          allocate(surf_points_grad_index(N_surface_points+1),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation (2) failed")
       end if

       counter = 1
       do i_ua=1,N_surface_points
          ua => unique_surf_points(i_ua)
          surf_points_grad_index(i_ua)=counter
          counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
       enddo
       surf_points_grad_index(N_surface_points+1)=counter
       totsym_field_length=counter-1

       do i=1,size(unique_surf_points)
          ua=>unique_surf_points(i)
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
       deallocate(unique_surf_points,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       if(.not.allocated(totalsym_field)) then
          allocate(totalsym_field(totsym_field_length),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation Totalsym_field failed")
       end if
       totalsym_field=0.0_r8_kind

       call comm_init_send(comm_all_other_hosts,msgtag_surf_point_sa)
       call surf_points_info_pack()
       call comm_send()
    else
       call surf_points_info_unpack()
    end if

  end subroutine surf_points_grad_information
  !******************************************************

  !******************************************************
  subroutine surf_points_grad_information1()
    !** End of interface *****************************************
    use unique_atom_module, only: unique_atom_type, unique_atom_partner_type, &
         unique_atom_symadapt_type
    use symmetry_data_module, only: get_totalsymmetric_irrep,find_totalsymmetric_irrep
    use uatom_symmadapt, only: sa_free
    ! ----------- declaration of local variables -------------
    type(unique_atom_type)         , pointer :: ua
    type(unique_atom_partner_type) , pointer :: sap
    type(unique_atom_symadapt_type), pointer :: sa
    type(arrmat3)                  , pointer :: gi
    integer(kind=i4_kind) :: i_ua,i_if,i_cd,i_ea,ts,counter,i,j,k, &
         n_indep,n_equals,alloc_stat
    logical :: do_alloc

    external error_handler
    ! ----------- executable code -----------------------------

    if(comm_i_am_master()) then
       do_alloc=.false.
       if(.not.associated(surf_points_grad_info)) then
          allocate( surf_points_grad_info(N_surface_points),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation (1) failed")
          do_alloc=.true.
       end if

       call find_totalsymmetric_irrep()
       ts = get_totalsymmetric_irrep()

       do i_ua=1,N_surface_points
          ua => unique_surf_points(i_ua)
          sap => ua%symadapt_partner(ts,1) ! l = 1
          gi => surf_points_grad_info(i_ua)
          n_indep = sap%N_independent_fcts
          n_equals = ua%N_equal_atoms
          if(do_alloc) then
             allocate ( gi%m(n_indep,3,n_equals),STAT=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("surf_point_grad_information: allocation (2) failed")
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
                   call error_handler("surf_point_grad_information: sth. wrong")
                end select

             enddo
          enddo
       end do

       if(do_alloc) then
          allocate(surf_points_grad_index(N_surface_points+1),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation (2) failed")
       end if

       counter = 1
       do i_ua=1,N_surface_points
          ua => unique_surf_points(i_ua)
          surf_points_grad_index(i_ua)=counter
          counter=counter+ua%symadapt_partner(ts,1)%N_independent_fcts
       enddo
       surf_points_grad_index(N_surface_points+1)=counter
       totsym_field_length=counter-1

       do i=1,size(unique_surf_points)
          ua=>unique_surf_points(i)
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
       deallocate(unique_surf_points,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       if(.not.allocated(totalsym_field)) then
          allocate(totalsym_field(totsym_field_length),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("surf_points_grad_information: allocation Totalsym_field failed")
       end if
       totalsym_field=0.0_r8_kind

       call comm_init_send(comm_all_other_hosts,msgtag_surf_point_sa)
       call surf_points_info_pack()
       call comm_send()
    else
       call comm_save_recv(comm_master_host,msgtag_surf_point_sa)
       call surf_points_info_unpack()
    end if

  end subroutine surf_points_grad_information1
  !******************************************************

  !*************************************************************
  subroutine surf_points_info_pack()
    !** End of interface *****************************************
    !------------ modules used ------------------- ---------------
    use comm_module, only: commpack
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,status
    integer(i4_kind) :: n1,n2,n3,nn
    !------------ Declaration of subroutines used ----------------
    !------------ Executable code ------------------------------------

    do i=1,N_surface_points
       n1=size(surf_points_grad_info(i)%m,1)
       n2=size(surf_points_grad_info(i)%m,2)
       n3=size(surf_points_grad_info(i)%m,3)
       nn=n1*n2*n3
       call commpack(n1,status)
       ASSERT(status==0)
       call commpack(n2,status)
       ASSERT(status==0)
       call commpack(n3,status)
       ASSERT(status==0)
       call commpack(surf_points_grad_info(i)%m(1,1,1),nn,1,status)
       ASSERT(status==0)
    end do

    call commpack(surf_points_grad_index(1),N_surface_points+1,1,status)
    ASSERT(status==0)
    call commpack(totsym_field_length,status)
    ASSERT(status==0)

    call commpack(totalsym_field(1),totsym_field_length,1,status)
    ASSERT(status==0)

  end subroutine surf_points_info_pack
  !*************************************************************

  !*************************************************************
  subroutine surf_points_info_unpack()
    !** End of interface *****************************************
    !------------ modules used ------------------- ---------------
    use comm_module, only: communpack
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,status
    integer(i4_kind) :: n1,n2,n3,nn
    logical :: do_alloc
    !------------ Executable code ------------------------------------

    do_alloc=.false.
    if(.not. associated(surf_points_grad_info)) then
       allocate(surf_points_grad_info(N_surface_points),stat=status)
       ASSERT(status==0)
       do_alloc=.true.
    end if

    do i=1,N_surface_points
       call communpack(n1,status)
       ASSERT(status==0)
       call communpack(n2,status)
       ASSERT(status==0)
       call communpack(n3,status)
       ASSERT(status==0)
       nn=n1*n2*n3

       if(do_alloc) then
          allocate(surf_points_grad_info(i)%m(n1,n2,n3),stat=status)
          ASSERT(status==0)
       end if

       call communpack(surf_points_grad_info(i)%m(1,1,1),nn,1,status)
       ASSERT(status==0)
    end do

    if(do_alloc) then
       allocate(surf_points_grad_index(N_surface_points+1),stat=status)
       ASSERT(status==0)
    end if

    call communpack(surf_points_grad_index(1),N_surface_points+1,1,status)
    ASSERT(status==0)
    call communpack(totsym_field_length,status)
    ASSERT(status==0)

    if(.not.allocated(totalsym_field)) then
       allocate(totalsym_field(totsym_field_length),stat=status)
       ASSERT(status==0)
    end if
    call communpack(totalsym_field(1),totsym_field_length,1,status)
    ASSERT(status==0)

  end subroutine surf_points_info_unpack
  !*************************************************************

  !*************************************************************
  subroutine surf_points_gradinfo_dealloc()
    ! Purpose: deallocate the variable 'surf_point_grad_info'
    ! subroutine called by: main_gradient
    !** End of interface *****************************************
    ! ----------- declaration of local variables -------------
    integer(kind=i4_kind) :: i, alloc_stat
    external error_handler
    ! ----------- executable code -----------------------------

    if (comm_i_am_master() .and. comm_parallel()) then
       call comm_init_send(comm_all_other_hosts,msgtag_sp_grinfo_dealloc)
       call comm_send()
    endif

    do i=1,N_surface_points
       deallocate(surf_points_grad_info(i)%m,STAT=alloc_stat)
       if (alloc_stat /= 0 ) call error_handler &
            ("surf_points_gradinfo_dealloc: deallocation (1) failed")
    enddo
    deallocate(surf_points_grad_info,STAT=alloc_stat)
    if (alloc_stat /= 0 ) call error_handler &
         ("surf_point_gradinfo_dealloc: deallocation (2) failed")

    deallocate(surf_points_grad_index,STAT=alloc_stat)
    if (alloc_stat /= 0 ) call error_handler &
         ("surf_point_gradinfo_dealloc: deallocation (3) failed")

    deallocate(totalsym_field,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("surf_points_grad_information: deallocation Totalsym_field failed")

  end subroutine surf_points_gradinfo_dealloc
  !*************************************************************

  !*********************************************************
  subroutine start_read_field_e
   !** End of interface *****************************************

    if(comm_parallel()) then
       call comm_init_send(comm_all_other_hosts,msgtag_start_field)
       call comm_send()
    endif

    call read_field_e
    call send_receive_field

  end subroutine start_read_field_e
  !*********************************************************

  !*********************************************************
  subroutine read_field_e
    !** End of interface *****************************************

    use density_data_module,only: densmat
    use readwriteblocked_module
    use integralstore_module, only : integralstore_3c_field
    !** End of interface *****************************************

    type(readwriteblocked_tapehandle) :: th_field
    real(kind=r8_kind),pointer :: field_int(:)
    integer (i4_kind) :: item_arr_field, my_ind, n_irrep, i_gamma, mul
    logical ::  spin_polarized,integrals_on_file
    integer(kind=i4_kind) :: n, status, m, k, i_meta, i_last
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
#ifdef _VECTOR
    integer(kind=i4_kind) :: i_mn, n_mn
    integer(kind=i4_kind), allocatable :: n_i(:),m_i(:)
#endif

    integrals_on_file=options_integrals_on_file()
    spin_polarized = symmetry_data_n_spin() > 1

    n_irrep = symmetry_data_n_irreps()

    allocate(dim_irrep(n_irrep),stat=status)
    if ( status .ne. 0) call error_handler( &
         "elec_static_field_module:read_field_e: allocated dim_irrep failed" )

    do n=1,n_irrep
       dim_irrep(n)  = symmetry_data_dimension(n)
    enddo

    my_ind=comm_myindex()
    item_arr_field=bounds%item_arr(my_ind)

    if (item_arr_field == 0) then
       deallocate(dim_irrep,stat=status)
       if ( status .ne. 0) call error_handler( &
            "elec_static_field_module:read_field_e: deallocated dim_irrep failed" )
       return
    endif

    allocate(E_n_tmp(item_arr_field), stat=status)
    if ( status .ne. 0) call error_handler( &
         "elec_static_field_module:read_field_e: allocated of E_n_tmp failed" )
    E_n_tmp=0.0_r8_kind

    if ( integrals_on_file ) then
       ! open the integral files ----------------------
       call field_integral_open(th_field)
    else
       i_meta=1
    endif

#ifndef _VECTOR
    if (integrals_on_file) then
       allocate( field_int(item_arr_field),STAT=status)
       if(status.ne.0) call error_handler&
            ("elec_static_field_module: allocation of field_int failed")
    endif
#endif

    do i_gamma = 1,n_irrep

#ifndef _VECTOR
       do m=1,dim_irrep(i_gamma)
          do n=1,m
             mul = 2
             if (n==m) mul = 1
             ! read in integrals file
             if (integrals_on_file) then
                call readwriteblocked_read(field_int,th_field)
             else
                i_last=i_meta+item_arr_field-1
                field_int=> integralstore_3c_field(i_meta:i_last)
                i_meta=i_last+1
             endif
             do k=1,item_arr_field
                if(spin_polarized) then
                   E_n_tmp(k)=E_n_tmp(k)+field_int(k)* &
                        mul * (densmat(i_gamma) % m(n, m, 1) + densmat(i_gamma) % m(n, m, 2))
                else
                   E_n_tmp(k)=E_n_tmp(k)+field_int(k)* &
                        mul *densmat(i_gamma) % m(n, m, 1)
                endif
             enddo
          enddo
       enddo
#else
       n_mn = ( dim_irrep(i_gamma) * (dim_irrep(i_gamma) + 1) ) / 2
       allocate(m_i(n_mn),n_i(n_mn))
       i_mn = 1
       do m=1,dim_irrep(i_gamma)
          do n=1,m
             m_i(i_mn)=m
             n_i(i_mn)=n
             i_mn=i_mn+1
          enddo
       enddo
       if ( integrals_on_file ) then
          allocate( field_int(item_arr_field*n_mn),STAT=status )
          if(status.ne.0) call error_handler &
               ("elec_static_module: allocation of field_int failed")
          call readwriteblocked_read(field_int,th_field)
          do k=1,item_arr_field
             do i_mn = 1, n_mn
                mul = 2
                if ( n_i(i_mn).eq.m_i(i_mn) ) mul = 1
                if( spin_polarized ) then
                   E_n_tmp(k)=E_n_tmp(k) +&
                        field_int(k+(i_mn-1)*item_arr_field)* &
                        mul * (densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1) + &
                        densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 2))
                else
                   E_n_tmp(k)=E_n_tmp(k) +&
                        field_int(k+(i_mn-1)*item_arr_field)* &
                        mul * densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)
                endif !spin_polarized
             enddo
          enddo
          deallocate(field_int,STAT=status)
          if(status.ne.0) call error_handler &
               ("elec_static_modul: deallocation of field_int failed")
       else ! integrals not on file
          do k=1,item_arr_fiel
             do i_mn = 1, n_mn
                mul = 2
                if (n_i(i_mn) .eq. m_i(i_mn)) mul = 1
                if(spin_polarized) then
                   E_n_tmp(k)=E_n_tmp(k)-&
                        integralstore_3c_field(i_meta-1+k+(i_mn-1)*item_arr_field)* &
                        mul * (densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)&
                        +densmat(i_gamma)%m(n_i(i_mn),m_i(i_mn),2))
                else
                   E_n_tmp(k)=E_n_tmp(k)-&
                        integralstore_3c_field(i_meta-1+k+(i_mn-1)*item_arr_field)* &
                        mul * densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)
                endif
             enddo !i_mn
          enddo !k
          i_meta=i_meta+item_arr_field*n_mn
       endif ! intgrals on file
       deallocate(m_i,n_i)
#endif
    enddo ! i_gamma

#ifndef _VECTOR
    if (integrals_on_file) then
       deallocate(field_int,STAT=status)
       if(status.ne.0) call error_handler&
            ("elec_static_field_module: deallocation of field_int failed")
    endif
#endif

    if ( integrals_on_file ) then
       ! close the integral files ----------------------
            call field_integral_close(th_field)
    endif

    deallocate(dim_irrep,stat=status)
    if ( status .ne. 0) call error_handler( &
         "elec_static_field_module:read_field_e: deallocated dim_irrep failed(1)" )

  end subroutine read_field_e
  !*********************************************************

  !*********************************************************
  subroutine send_receive_field
   !** End of interface *****************************************

    integer(kind=i4_kind) :: my_ind,field_item_arr,bound
    integer(kind=i4_kind) :: status,info,i_pr,wa
    integer(kind=i4_kind) :: i

    if(comm_parallel()) then
       if(comm_i_am_master()) then
          if(calc_normal.and. .not.allocated(E_n)) then
             allocate(E_n(N_surface_points),stat=status)
             if ( status .ne. 0) call error_handler( &
                  "field: allocated of E_n failed" )
          else
                 totalsym_field=0.0_r8_kind
          endif

          i_pr=comm_get_n_processors()
          bound=0
          do i=1,i_pr
             if (i == 1) then
                field_item_arr=bounds%item_arr(1)
                if(field_item_arr == 0) cycle
                if(calc_normal) then
                   E_n(bound+1:bound+field_item_arr)= E_n_tmp
                else
                   totalsym_field(bound+1:bound+field_item_arr)=E_n_tmp
                end if
                deallocate(E_n_tmp,stat=status)
                if ( status .ne. 0) call error_handler( &
                     "field: deallocated of E_n_tmp failed(1)" )
             else
                call comm_save_recv(i,msgtag_finish_field)

                call communpack(field_item_arr,info)
                if( info.ne.0) call error_handler &
                     ("field:unpack of field_item_arr failed")
                if (field_item_arr == 0) cycle

                if(calc_normal) then
                   call communpack(E_n(bound+1:bound+field_item_arr),field_item_arr,1,info)
                else
                   call communpack(totalsym_field(bound+1:bound+field_item_arr),field_item_arr,1,info)
                end if
                if( info.ne.0) call error_handler &
                     ("field:unpack of E_n failed")
             endif
             bound=bound+field_item_arr
          enddo
       else
          call comm_init_send(comm_master_host,msgtag_finish_field)

          my_ind=comm_myindex()

          call commpack(bounds%item_arr(my_ind),info)
          if( info.ne.0) call error_handler &
               ("field:pack of item_arr failed")

          if(bounds%item_arr(my_ind) /= 0) then
             call commpack(E_n_tmp,bounds%item_arr(my_ind),1,info)
             if( info.ne.0) call error_handler &
                  ("field:pack of E_n_tmp failed")

             deallocate(E_n_tmp,stat=status)
             if ( status .ne. 0) call error_handler( &
                  "field: deallocated of E_n_tmp failed(2)" )
          endif
          call comm_send()
       endif
    else
       if (calc_normal.and. .not.allocated(E_n)) then
          allocate(E_n(N_surface_points),stat=status)
          if ( status .ne. 0) call error_handler( &
               "field: allocated of E_n failed(1)" )
       endif
       if (calc_normal) then
          E_n=E_n_tmp
       else
          totalsym_field=E_n_tmp
       end if

       deallocate(E_n_tmp,stat=status)
       if ( status .ne. 0) call error_handler( &
            "field: deallocated of E_n_tmp failed(3)" )

    endif

    if(comm_i_am_master() .and. .not. calc_normal)then
       if(.not. allocated(E_ele)) then
          allocate(E_ele(N_surface_points),stat=status)
          if ( status .ne. 0) call error_handler( &
               "field:get_field_nuc: allocated E_ele failed" )
       end if

       do i=1,N_surface_points
          wa=surface_points(i)%N_equal_points
          allocate(E_ele(i)%m(3,wa),stat=status)
          if (status .ne. 0) call error_handler( &
               "field:send_receive_field: allocated E_ele%m failed" )
       end do
       call transform_to_cart_field(E_ele,totalsym_field)
       totalsym_field=0.0_r8_kind
    end if

  end subroutine send_receive_field
  !*********************************************************

  !*********************************************************
  subroutine get_field_nuc()
    use unique_atom_module, only : unique_atoms, N_unique_atoms, &
         pseudopot_present
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind)       :: f_dim,index,wa
    real(kind=r8_kind)          :: z1,dist,zc1
    integer(kind=i4_kind)       :: na, eq_a, i, j, status
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:),rotmat(:,:)
    real(r8_kind) :: EE(3),rr(3)

    if(.not. allocated(E_nuc)) then
       allocate(E_nuc(N_surface_points),stat=status)
       if ( status .ne. 0) call error_handler( &
            "field:get_field_nuc: allocated E_nuc failed" )
    end if

    totalsym_field=0.0_r8_kind

    do i=1,N_surface_points
       EE=0.0_r8_kind
       xb => surface_points(i)%position
       wa=surface_points(i)%N_equal_points
       if(.not. allocated(E_nuc(i)%m)) then
          allocate(E_nuc(i)%m(3,wa),stat=status)
          if (status .ne. 0) call error_handler( &
               "field:get_field_nuc: allocated E_nuc%m failed" )
       end if

       unique_a: do na=1,N_unique_atoms
          z1 = unique_atoms(na)%Z
          zc1 = unique_atoms(na)%zc
          if (.not.pseudopot_present) zc1 = 0.0_r8_kind
          xa => unique_atoms(na)%position

          equal_a: do eq_a=1,unique_atoms(na)%N_equal_atoms
             rr=xa(:,eq_a)-xb(:,1)
             dist=sqrt(dot_product(rr,rr))
             if(dist < 0.001_r8_kind) dist=0.001_r8_kind
             EE = EE + (z1-zc1)*rr/dist**3
          enddo equal_a
       enddo unique_a

       f_dim=surf_points_grad_index(i+1)-surf_points_grad_index(i)
       rotmat=>surf_points_grad_info(i)%m(:,:,1)
       index=surf_points_grad_index(i)
       do j=1,f_dim
          totalsym_field(index) = totalsym_field(index) - &
               wa*sum(rotmat(j,:)*EE(:))
          index=index+1
       end do
    end do

    call transform_to_cart_field(E_nuc,totalsym_field)

    totalsym_field=0.0_r8_kind

  end subroutine get_field_nuc
  !*********************************************************

  !*********************************************************
  subroutine get_field_pc()
   !** End of interface *****************************************

    use pointcharge_module, only : pointcharge_array, pointcharge_N
#ifdef WITH_EPE
    use ewaldpc_module, only : ewpc_array, EWPC_N
#endif

    integer(kind=i4_kind)       :: f_dim,index,wa
    real(kind=r8_kind)          :: z1,dist
    integer(kind=i4_kind)       :: na, eq_a, i, j, status
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:),rotmat(:,:)
    real(r8_kind) :: EE(3),rr(3)

    allocate(E_pc(N_surface_points),stat=status)
    if ( status .ne. 0) call error_handler( &
         "field:get_field_pc: allocated E_pc failed" )

    totalsym_field=0.0_r8_kind

    do i=1,N_surface_points
       EE=0.0_r8_kind
       xb => surface_points(i)%position
       wa=surface_points(i)%N_equal_points
       allocate(E_pc(i)%m(3,wa),stat=status)
       if (status .ne. 0) call error_handler( &
            "field:get_field_pc: allocated E_pc%m failed" )

       unique_a: do na=1,pointcharge_N
          z1 = pointcharge_array(na)%Z
          xa => pointcharge_array(na)%position

          equal_a: do eq_a=1,pointcharge_array(na)%N_equal_charges
             rr=xa(:,eq_a)-xb(:,1)
             dist=sqrt(dot_product(rr,rr))
             if(dist < 0.001_r8_kind) dist=0.001_r8_kind
             EE = EE + z1*rr/dist**3
          enddo equal_a
       enddo unique_a

#ifdef WITH_EPE
       unique_ae: do na=1,EWPC_N
          z1 = ewpc_array(na)%Z
          xa => ewpc_array(na)%position

          equal_ae: do eq_a=1,ewpc_array(na)%N_equal_charges
             rr=xa(:,eq_a)-xb(:,1)
             dist=sqrt(dot_product(rr,rr))
             if(dist < 0.001_r8_kind) dist=0.001_r8_kind
             EE = EE + z1*rr/dist**3
          enddo equal_ae
       enddo unique_ae
#endif

       f_dim=surf_points_grad_index(i+1)-surf_points_grad_index(i)
       rotmat=>surf_points_grad_info(i)%m(:,:,1)
       index=surf_points_grad_index(i)
       do j=1,f_dim
          totalsym_field(index) = totalsym_field(index) - &
               wa*sum( rotmat(j,:)*EE(:) )
          index=index+1
       end do
    end do

    call transform_to_cart_field(E_pc,totalsym_field)

    totalsym_field=0.0_r8_kind

  end subroutine get_field_pc
  !*********************************************************

  !*********************************************************
  subroutine transform_to_cart_field(cart_field,totsym_field)
    !** End of interface *****************************************

    type(arrmat2), intent(inout) :: cart_field(:)   ! cartesian field array
    real(r8_kind), intent(in)    :: totsym_field(:) ! symm. adapt. contr.
    integer(kind=i4_kind) :: i_unique,i_equal,index,i,f_dim
    real(kind=r8_kind),pointer :: rotmat(:,:)

    do i_unique=1,N_surface_points
       index=surf_points_grad_index(i_unique)
       f_dim=surf_points_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,surface_points(i_unique)%N_equal_points
          rotmat=>surf_points_grad_info(i_unique)%m(:,:,i_equal)
          cart_field(i_unique)%m(:,i_equal) = 0.0_r8_kind
          do i=1,f_dim
             cart_field(i_unique)%m(:,i_equal) = &
                  cart_field(i_unique)%m(:,i_equal) - &  !+
                  rotmat(i,:)*totsym_field(index+i)
          end do
       end do
    enddo

  end subroutine transform_to_cart_field
  !*********************************************************

  !*********************************************************
  subroutine deallocate_field
    !** End of interface *****************************************

    integer(kind=i4_kind) :: status,i

    if(allocated(E_n)) then
       deallocate(E_n,stat=status)
       if ( status .ne. 0) call error_handler( &
            "field:deallocation of E_n is failed" )
    end if
    if(allocated(E_ele)) then
       do i=1,N_surface_points
          deallocate(E_ele(i)%m,stat=status)
          if ( status .ne. 0) call error_handler( &
               "field:deallocation of E_ele%m is failed" )
       end do
       deallocate(E_ele,stat=status)
       if ( status .ne. 0) call error_handler( &
            "field:deallocation of E_ele is failed" )
    end if
    if(allocated(E_nuc)) then
       do i=1,N_surface_points
          deallocate(E_nuc(i)%m,stat=status)
          if ( status .ne. 0) call error_handler( &
               "field:deallocation of E_nuc%m is failed" )
       end do
       deallocate(E_nuc,stat=status)
       if ( status .ne. 0) call error_handler( &
            "field:deallocation of E_nuc is failed" )
    end if
    if(allocated(E_pc)) then
       do i=1,N_surface_points
          deallocate(E_pc(i)%m,stat=status)
          if ( status .ne. 0) call error_handler( &
               "field:deallocation of E_pc%m is failed" )
       end do
       deallocate(E_pc,stat=status)
       if ( status .ne. 0) call error_handler( &
            "field:deallocation of E_pc is failed" )
    end if

  end subroutine deallocate_field
  !*********************************************************

  !**********************************************************
  subroutine field_integral_open(th_field2)
   !** End of interface *****************************************
    !  Purpose: decide which files to open, perform 'startread'
    ! Modules:
    use filename_module, only: tmpfile
    use readwriteblocked_module
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle),intent(out) :: th_field2
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    call readwriteblocked_startread(trim(tmpfile('field.dat')), th_field2)

  end subroutine field_integral_open
  !*************************************************************

  !*************************************************************
  subroutine field_integral_close(th_field1)
    !  Purpose: performs a 'stopread' on the files contained
    !           in the specified tapehandle,
    !           returns the corresponding io_unit and
    !           closes the files
    use readwriteblocked_module
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle),intent(inout) :: th_field1
    !** End of interface *****************************************
    call readwriteblocked_stopread(th_field1)
  end subroutine field_integral_close
  !*************************************************************

  !*************************************************************
  subroutine destroy_field_file
    ! Purpose: to delete files of electric field matrix elements
    !          from disk
   !** End of interface *****************************************
    use filename_module, only: tmpfile
    use readwriteblocked_module

    logical :: yes
    type(readwriteblocked_tapehandle) :: th_field3

    inquire(file=trim(tmpfile('field.dat')), exist=yes)
    if (yes) then
       call readwriteblocked_startread(trim(tmpfile('field.dat')), th_field3)
       call readwriteblocked_stopread(th_field3,status='delete')
    endif

  end subroutine destroy_field_file
  !*************************************************************

  !*********************************************************
  subroutine bounds_calc_field ()
    use comm, only: comm_size, comm_same
    !** End of interface *****************************************

    real(kind=r8_kind)      :: rhelp
    integer(kind=i4_kind)   :: i,alloc_stat,i_pr,ind,i_size

    ! Make sure the input is the same on all workers:
    if (.not. comm_same(N_surface_points)) stop "FIX!"
    if (.not. comm_same(totsym_field_length)) stop "FIX!"

    i_pr = comm_size ()

    allocate (bounds%item_arr(i_pr),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler&
         ("field: allocation of BOUNDS failed ")
    bounds%item_arr=0

    if (comm_parallel()) then
       rhelp = real(N_surface_points,kind=r8_kind)
       if(calc_normal) then
             do i=1,i_pr
                bounds%item_arr(i) = ceiling(rhelp/i_pr)
                rhelp = rhelp - real(bounds%item_arr(i),kind=r8_kind)
                i_pr = i_pr - 1
             enddo

             bounds%lower1 = bounds%item_arr(1)+1
             bounds%upper1 = bounds%item_arr(1) + bounds%item_arr(2)
       else
             ind=1
             do i=1,i_pr
                i_size=ceiling(rhelp/i_pr)
                bounds%item_arr(i)=surf_points_grad_index(i_size+ind)- &
                     surf_points_grad_index(ind)
                rhelp = rhelp - real(i_size,kind=r8_kind)
                i_pr = i_pr - 1
                ind=ind+i_size
             end do
             bounds%lower1 = bounds%item_arr(1)+1
             bounds%upper1 = bounds%item_arr(1) + bounds%item_arr(2)
       end if
    else  ! no parallel processing
       if(calc_normal) then
          bounds%item_arr(1) = N_surface_points
       else
          bounds%item_arr(1) = totsym_field_length
       end if
       bounds%lower1 = 1
       bounds%upper1 = 0
    endif! end of if( comm_parallel )

  end subroutine bounds_calc_field
  !******************************************************************

  !*************************************************************
  subroutine bounds_free_field()
    !  Purpose: Release the private bounds pointer
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: status
    !------------ Executable code ------------------------------------

    if (comm_parallel() .and. comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_free_bnds_fld)
       call comm_send()
    endif

    deallocate (bounds%item_arr,STAT=status)
    if (status /= 0) call error_handler&
         ("field: deallocation of BOUNDS failed ")
  end subroutine bounds_free_field
  !*************************************************************

end module elec_static_field_module

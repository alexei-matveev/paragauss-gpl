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
! Public interface of module
!=====================================================================
module solv_electrostat_module
!
!  This module are used to calculate electrostatic
!  contribution to the solvent effect of molecule
!
!  The module was prepared by extracting corresponding routines from
!  old solvation_module

!== Interrupt of public interface of module =========
!  Author: AS
!  Date: 07/69
!
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!== Interrupt end of public interface of module =====

!----modules used ------------------
! use FPP_USR_TIMERs:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only: i4_kind, r8_kind ! type specification parameters
  use datatype, only: arrmat2
  use filename_module, only: inpfile
  use iounitadmin_module, only: output_unit, write_to_output_units, &
       openget_iounit, returnclose_iounit
  use solv_cavity_module, only: tessarea, to_calc_grads, Q_e, Q_n, &
       Q_e_old, n_size, center2sphere, with_pc, fixed_pc, &
       cor_el, cor_nuc, &
       grad_solv_totsym, grad_solv_totsym_tes, n_Q_update
#ifdef WITH_EFP
  use efp_module, only: efp, n_efp
#endif

  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------

  public :: surface_charge_moments       ! (charge, dipole), both out
  public :: const_part_of_ham_and_energy ! () is called by all workers
  public :: build_solv_ham               ! () is called by all workers

  public solv_poten_transfer_data, &
       calc_Q_e, &
       sphere2center, &
       alloc_ham_solv, &
       dealloc_ham_solv, solv_energy_el,&
       matrix_generation, matrix_grad,nuc_grad,matrix_grad_vtn,nuc_grad_vtn, &
       shutdown_solvation, &
       init_forces_on_pc,solv_forces_on_pc,dealloc_solv_pc

  !===================================================================
!== Interrupt end of public interface of module =====
!--public variables ---

  ! (A  half)  of  the  interaction energy  between  molecular  charge
  ! density including electrons and nuclei and the polarized solvent:
  real (r8_kind), public, protected :: energy_solv_el

  ! Additional   part   of  hamiltonian   due   to  interaction   with
  ! polarization charges of the solvent:
  type (arrmat2), allocatable, public, protected :: ham_solv_el(:)
#if 0
  type(arrmat2),allocatable,public :: ham_solv_el_keep(:)
#endif
  ! End of public interface of module

  real(kind=r8_kind), allocatable :: A_matrix(:,:)
  real (r8_kind), public, allocatable :: A_matrix_inv(:,:) ! not protected

  real (r8_kind), public, allocatable :: Q_grad(:,:) ! not protected
  real (r8_kind) :: energy_core_n

  type(arrmat2), allocatable, public, protected :: F_solv_pc(:)

  real(kind=r8_kind) :: correct_factor_n
  external error_handler

contains
#if 0
  !********************************************************************
   subroutine fix_A_mat_inv()
    integer(kind=i4_kind):: fix_mat_unit
    logical:: fix_mat
     inquire(file=trim(inpfile("fixmat.dat")), exist=fix_mat)
     if(fix_mat) then
     print*, 'solv matrix fixed'
               fix_mat_unit=openget_iounit(file=trim(inpfile('fixmat.dat')), &
                               form='FORMATTED',status='old')

           read(fix_mat_unit,*) A_matrix_inv

           call returnclose_iounit(fix_mat_unit,status='keep')
    endif
    fix_mat_unit=openget_iounit(file=trim(inpfile('fixmat.dat')), &
                               form='FORMATTED',status='unknown')
           write(fix_mat_unit,*) A_matrix_inv
           call returnclose_iounit(fix_mat_unit,status='keep')
    end subroutine fix_A_mat_inv

   subroutine fix_cavity()
   integer(kind=i4_kind):: fix_cav_unit,i,l
   logical:: fix_cav

     inquire(file=trim(inpfile("fixcav.dat")), exist=fix_cav)
     if(fix_cav) then
     print*, 'cavity fixed'
               fix_cav_unit=openget_iounit(file=trim(inpfile('fixcav.dat')), &
                               form='FORMATTED',status='old')

      do i=1,size(tessarea)
      do l=1,tessarea(i)%n_equal
           read(fix_cav_unit,*) tessarea(i)%area,tessarea(i)%xyz(l,:)
      enddo
      enddo
           call returnclose_iounit(fix_cav_unit,status='keep')
    endif
    fix_cav_unit=openget_iounit(file=trim(inpfile('fixcav.dat')), &
                               form='FORMATTED',status='unknown')
      do i=1,size(tessarea)
      do l=1,tessarea(i)%n_equal
           write(fix_cav_unit,*) tessarea(i)%area,tessarea(i)%xyz(l,:)
      enddo
      enddo
           call returnclose_iounit(fix_cav_unit,status='keep')
  end subroutine fix_cavity
#endif

  !********************************************************************

  subroutine matrix_generation ()
    !
    ! Generate the "interaction" matrix of charged surface areas.  The
    ! inverse = matrix for calculating the induced surface charge from
    ! the molecule potential on the surface. Does no communucation.
    !
    use math_module, only: invert_matrix
    implicit none
    !** End of interface *****************************************

    real (r8_kind), parameter :: pi = 3.14159265355897932368_r8_kind
    integer (i4_kind) :: i, k, l, status
    integer (i4_kind) :: first_index, second_index, first_equal, second_equal
    real (r8_kind) :: distance, A_matrix_fs
    real (r8_kind) :: vect(3)

    ! Generating direct matrix A for the main COSMO equation
    allocate (A_matrix_inv(n_size, n_size), stat=status)
    if (status /= 0) call error_handler &
         ("matrix_generation: allocation of A_MATRIX is failed")

!!$call cpu_time(tt)
!!$print*,tt
#if 0
   call fix_cavity()
#endif

    first_index = 1
    do i = 1, n_size
       first_equal = tessarea(i) % n_equal
       second_index = 1
       do k = 1, n_size
          second_equal = tessarea(k) % n_equal
          A_matrix_inv(i, k) = 0.0
          do l = 1, second_equal
             f_s: if (first_index == second_index) then
                A_matrix_fs = &
                     1.07_r8_kind * sqrt (4 * pi / tessarea(i) % area)
             else
                vect = tessarea(i) % xyz(1, :) - tessarea(k) % xyz(l, :)
                distance = sqrt (dot_product (vect, vect))
                A_matrix_fs = 1 / distance
             endif f_s
             A_matrix_inv(i, k) = A_matrix_inv(i, k) + A_matrix_fs
             second_index = second_index + 1
          enddo
       enddo
       first_index = first_index + first_equal
    enddo
!!$print*,sum(A_matrix_inv(:,:)),'AAAAAAAAAAAA'


!!$call cpu_time(tt)
!!$print*,tt

    call invert_matrix (n_size, A_matrix_inv)
#if 0
   call fix_A_mat_inv()
#endif

!!$call cpu_time(tt)
!!$print*,tt

!     deallocate(A_matrix,stat=status)
!     if ( status /= 0) call error_handler( &
!          "matrix_generation: deallocation of A_MATRIX is failed")
  end subroutine matrix_generation
  !********************************************************************

  !***************************************************
  subroutine dealloc_A_inv
   ! deallocates the inverse of the cavity matrix A
   ! (A_inv is needed to calculate the induced surface charges)
   !** End of interface *****************************************
    integer(kind=i4_kind) :: status

    deallocate(A_matrix_inv,stat=status)
    if ( status /= 0) call error_handler( &
         "dealloc_A_inv: deallocation of A_MATRIX_INV is failed")

  end subroutine dealloc_A_inv
  !***************************************************

  !***************************************************
  subroutine matrix_grad(grad_index)
    !gradient of the cavity matrix A
    use unique_atom_module
    use pointcharge_module
    use elec_static_field_module
    use comm, only: comm_size, comm_rank
    implicit none
    integer(kind=i4_kind), intent(in) :: grad_index(N_moving_unique_atoms + 1)
    !** End of interface *****************************************

    real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind
    integer(kind=i4_kind) :: i, i1, j, k, l, ng, ng1, na, ea !, m
    integer(kind=i4_kind) :: first_index,second_index,first_equal,second_equal
    real(kind=r8_kind) :: distance,help,eps_help,help_d,ds(3),dxyz1(3,3),dxyz2(3,3)
    real(kind=r8_kind) :: vect(3), Qi, Q1i, Qk !,tt
    integer(kind=i4_kind) :: index,grad_dim,N_pc
    integer(i4_kind) :: n_proc,my_index,N_species,n_start,n_length

    eps_help=to_calc_grads%dielconst/(2.0_r8_kind*(to_calc_grads%dielconst-1.0_r8_kind))

    N_pc=0
    if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

    DPRINT 'matrix_grad: entered'
    DPRINT 'matrix_grad: grad_index=',grad_index

    N_species=N_moving_unique_atoms+N_pc
    n_proc = comm_size()
    my_index = 1 + comm_rank()
    n_start=1
    n_length=0
    do i=1,my_index
       n_start=n_start+n_length
       N_species=N_species-n_length
       n_length=int(N_species/n_proc)
       n_proc=n_proc-1
    end do

    if(n_length == 0) return
!!$call cpu_time(tt)
!!$print*,tt
    ng_N: do ng=n_start,n_start+n_length-1
       if(ng <= N_moving_unique_atoms) then
          na = moving_unique_atom_index(ng)
          ea = unique_atoms(na)%n_equal_atoms
       else
          na=ng-N_moving_unique_atoms
          ea=pointcharge_array(na)%N_equal_charges
       end if
!!$       m_ea: do m=1,ea
          first_index=0
          i_n_size: do i=1,to_calc_grads%n_points
             first_equal=to_calc_grads%n_equal(i)

             Qi=to_calc_grads%Q(i)
             Q1i=Qi
             help_d=-1.07_r8_kind*Qi*Q1i*sqrt(pi/to_calc_grads%s(i))/&
                  to_calc_grads%s(i)

             j_first_equal: do j=1,first_equal
                first_index=first_index+1

                second_index=0
                k_n_size: do k=1,to_calc_grads%n_points
                   second_equal=to_calc_grads%n_equal(k)

                   Qk=to_calc_grads%Q(k)

                   l_second_equal: do l=1,second_equal
                      second_index=second_index+1
                      if(second_index < first_index) cycle l_second_equal
                      f_s: if (first_index==second_index) then
                         ds(1)=to_calc_grads%ds_totsyms(1,ng)%m(to_calc_grads%i_symm_sort(i,j))
                         ds(2)=to_calc_grads%ds_totsyms(2,ng)%m(to_calc_grads%i_symm_sort(i,j))
                         ds(3)=to_calc_grads%ds_totsyms(3,ng)%m(to_calc_grads%i_symm_sort(i,j))
                         if(dot_product(ds,ds) < 1.0e-12_r8_kind) cycle

                         if(ng <= N_moving_unique_atoms) then
                            grad_dim=grad_index(ng+1)-grad_index(ng)
                            index=grad_index(ng)
                            do i1=1,grad_dim
                               grad_solv_totsym(index)=grad_solv_totsym(index) + &
                                    help_d*eps_help*ds(i1)
                               index=index+1
                            end do
                         else
                            ng1=ng-N_moving_unique_atoms
                            grad_dim=surf_points_grad_index(ng1+1)-surf_points_grad_index(ng1)
                            index=surf_points_grad_index(ng1)
                            do i1=1,grad_dim
                               totalsym_field(index)=totalsym_field(index) - &
                                    help_d*eps_help*ds(i1)
                               index=index+1
                            end do
                         end if
                      else
                         dxyz1(1,:)=to_calc_grads%dxyz_totsyms(1,ng)%m(:,to_calc_grads%i_symm_sort(i,j))
                         dxyz1(2,:)=to_calc_grads%dxyz_totsyms(2,ng)%m(:,to_calc_grads%i_symm_sort(i,j))
                         dxyz1(3,:)=to_calc_grads%dxyz_totsyms(3,ng)%m(:,to_calc_grads%i_symm_sort(i,j))
                         dxyz2(1,:)=to_calc_grads%dxyz_totsyms(1,ng)%m(:,to_calc_grads%i_symm_sort(k,l))
                         dxyz2(2,:)=to_calc_grads%dxyz_totsyms(2,ng)%m(:,to_calc_grads%i_symm_sort(k,l))
                         dxyz2(3,:)=to_calc_grads%dxyz_totsyms(3,ng)%m(:,to_calc_grads%i_symm_sort(k,l))
                         if(sum((dxyz1-dxyz2)*(dxyz1-dxyz2)) < 1.0e-12_r8_kind) cycle

                         vect=to_calc_grads%xyz(i,j,:)-to_calc_grads%xyz(k,l,:)
                         distance=sqrt(dot_product(vect,vect))
                         help=-1.0_r8_kind/distance**3
                         if(ng <= N_moving_unique_atoms) then
                            grad_dim=grad_index(ng+1)-grad_index(ng)
                            index=grad_index(ng)
                            do i1=1,grad_dim
                               grad_solv_totsym(index)=grad_solv_totsym(index) + &
                                    help*2.0_r8_kind*eps_help*Q1i*Qk* &
                                    dot_product(vect, &
                                    to_calc_grads%dxyz_totsyms(i1,ng)%m(:,to_calc_grads%i_symm_sort(i,j))- &
                                    to_calc_grads%dxyz_totsyms(i1,ng)%m(:,to_calc_grads%i_symm_sort(k,l)))
                               index=index+1
                            end do
                         else
                            ng1=ng-N_moving_unique_atoms
                            grad_dim=surf_points_grad_index(ng1+1)-surf_points_grad_index(ng1)
                            index=surf_points_grad_index(ng1)
                            do i1=1,grad_dim
                               totalsym_field(index)=totalsym_field(index) - &
                                    help*2.0_r8_kind*eps_help*Q1i*Qk* &
                                    dot_product(vect, &
                                    to_calc_grads%dxyz_totsyms(i1,ng)%m(:,to_calc_grads%i_symm_sort(i,j))- &
                                    to_calc_grads%dxyz_totsyms(i1,ng)%m(:,to_calc_grads%i_symm_sort(k,l)))
                               index=index+1
                            end do
                         endif
                      end if f_s
                   enddo l_second_equal
                enddo k_n_size
             enddo j_first_equal
          enddo i_n_size
!!$       enddo m_ea
    enddo ng_N
!!$call cpu_time(tt)
!!$print*,tt

  end subroutine matrix_grad
  !***************************************************

  !***************************************************
  subroutine matrix_grad_vtn(grad_index)
    !gradient of the cavity matrix A (variable tesserae number approximation)
    !------------ Modules used -----------------------------------
    use unique_atom_module
#ifdef WITH_EFP
    use pointcharge_module, only: pointcharge_N,pointcharge_array, &
         unique_pc_index,unique_pc_grad_info,rcm
    use efp_solv_grad_module, only: gradient_mpole_solv_totalsym, &
         torque_mpole_solv_totalsym
#endif
    use comm, only: comm_size, comm_rank
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: grad_index(N_moving_unique_atoms + 1)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind
    integer(kind=i4_kind) :: i, i1, j, k, l, m, ng, na, na1, ea
    integer(kind=i4_kind) :: first_index,second_index,first_equal,second_equal
    real(kind=r8_kind) :: distance,help,eps_help
    real(kind=r8_kind) :: vect(3), gradient(3), Qi, Q1i, Qk, grad(3)
    real(kind=r8_kind),pointer  :: rotmat(:,:)
    integer(kind=i4_kind) :: index,grad_dim,N_xc
    integer(i4_kind) :: n_proc,my_index,N_species,n_start,n_length
    integer(i4_kind) :: sphere_c, sphere_t1, sphere_t2
#ifdef WITH_EFP
    integer(kind=i4_kind) :: ng1
    real(kind=r8_kind) :: torque(3),rc(3),vect_c(3)
    integer(i4_kind) :: center1(2), center2(2),grp,grp1,grp2
#endif
    !------------ Executable code ------------------------------------

    eps_help=to_calc_grads%dielconst/(2.0_r8_kind*(to_calc_grads%dielconst-1.0_r8_kind))

    N_xc=0
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) N_xc=pointcharge_N
#endif

    N_species=N_moving_unique_atoms+N_xc
    n_proc = comm_size()
    my_index = 1 + comm_rank()
    n_start=1
    n_length=0
    do i=1,my_index
       n_start=n_start+n_length
       N_species=N_species-n_length
       n_length=int(N_species/n_proc)
       n_proc=n_proc-1
    end do

    if(n_length == 0) return
    ng_N: do ng=n_start,n_start+n_length-1
       if (ng <= N_moving_unique_atoms) then
          na = moving_unique_atom_index(ng)
          ea = unique_atoms(na) % n_equal_atoms
          na1 = na
#ifdef WITH_EFP
       elseif (ng > N_moving_unique_atoms) then
          na = ng - N_moving_unique_atoms
          ea = pointcharge_array(na) % N_equal_charges
          na1 = na + N_unique_atoms
#endif
       else
          na = -1
          ea = -1
          na1 = -1
          ABORT ("should not happen!")
       end if
       m_ea: do m=1,ea
          gradient=0.0_r8_kind
#ifdef WITH_EFP
          if(ng > N_moving_unique_atoms) then
             torque=0.0_r8_kind
             grp=pointcharge_array(na)%group(m)
             rc=rcm(:,grp)
          end if
#endif
          sphere_c=center2sphere(na1,m)

          first_index=0
          i_n_size: do i=1,to_calc_grads%n_points
             first_equal=to_calc_grads%n_equal(i)

             Qi=to_calc_grads%Q(i)
#ifdef WITH_EFP
             Q1i=to_calc_grads%Q1(i)
#else
             Q1i=Qi
#endif
             j_first_equal: do j=1,first_equal
                first_index=first_index+1

                second_index=0
                k_n_size: do k=1,to_calc_grads%n_points
                   second_equal=to_calc_grads%n_equal(k)

                   Qk=to_calc_grads%Q(k)

                   l_second_equal: do l=1,second_equal
                      second_index=second_index+1
                      if(second_index < first_index) cycle l_second_equal
                      f_s: if (first_index==second_index) then
                         cycle l_second_equal
                      else
                         sphere_t1=to_calc_grads%sphere(i,j)
                         sphere_t2=to_calc_grads%sphere(k,l)

                         if(sphere_c==sphere_t1 .and. sphere_c==sphere_t2) cycle l_second_equal
                         if(sphere_c/=sphere_t1 .and. sphere_c/=sphere_t2) cycle l_second_equal

#ifdef WITH_EFP
                         center1=sphere2center(sphere_t1)
                         center2=sphere2center(sphere_t2)
                         if(ng > N_moving_unique_atoms .and. &
                              center1(1) > N_unique_atoms .and. &
                              center2(1) > N_unique_atoms) then
                            grp1=pointcharge_array(center1(1)-N_unique_atoms)%group(center1(2))
                            grp2=pointcharge_array(center2(1)-N_unique_atoms)%group(center2(2))
                            if(grp1 == grp2) cycle l_second_equal
                         end if
#endif
                         vect=to_calc_grads%xyz(i,j,:)-to_calc_grads%xyz(k,l,:)
                         distance=sqrt(dot_product(vect,vect))
                         if(sphere_c==sphere_t2) distance=-distance
                         help=-1.0_r8_kind/distance**3
                         grad=help*2.0_r8_kind*eps_help*Q1i*Qk*vect
                         gradient=gradient+grad

#ifdef WITH_EFP
                         if(ng > N_moving_unique_atoms) then
                           if(sphere_c==sphere_t1) then
                               vect_c=to_calc_grads%xyz(i,j,:)-rc
                            elseif(sphere_c==sphere_t2) then
                               vect_c=to_calc_grads%xyz(k,l,:)-rc
                            end if
                            torque=torque-vector_product(grad,vect_c)
                         end if
#endif
                      end if f_s
                   enddo l_second_equal
                enddo k_n_size
             enddo j_first_equal
          enddo i_n_size
          if(ng <= N_moving_unique_atoms) then
             grad_dim=grad_index(ng+1)-grad_index(ng)
             index=grad_index(ng)
             rotmat=>unique_atom_grad_info(ng)%m(:,:,m)
             do i1=1,grad_dim
                grad_solv_totsym(index)=grad_solv_totsym(index) + &
                     sum(rotmat(i1,:)*gradient(:))
                index=index+1
             end do
#ifdef WITH_EFP
          elseif(ng > N_moving_unique_atoms .and. ng <= N_moving_unique_atoms+pointcharge_N) then
             ng1=ng-N_moving_unique_atoms
             grad_dim=unique_pc_index(ng1+1)-unique_pc_index(ng1)
             index=unique_pc_index(ng1)
             rotmat=>unique_pc_grad_info(ng1)%m(:,:,m)
             do i1=1,grad_dim
                gradient_mpole_solv_totalsym(index)=gradient_mpole_solv_totalsym(index) + &
                     sum(rotmat(i1,:)*gradient(:))
                torque_mpole_solv_totalsym(index)=torque_mpole_solv_totalsym(index) + &
                     sum(rotmat(i1,:)*torque(:))
                index=index+1
             end do
#endif
          endif
       enddo m_ea
    enddo ng_N

  end subroutine matrix_grad_vtn
  !***************************************************

  !*************************************************************
  function sphere2center(i_sphere) result (center)

    use unique_atom_module, only: N_unique_atoms, unique_atoms
    use pointcharge_module, only: pointcharge_N, pointcharge_array

    integer(i4_kind) :: center(2)
    integer(i4_kind) :: i_sphere

    integer(i4_kind) :: is,js,ks

    do is=1,N_unique_atoms+pointcharge_N
       if(is <= N_unique_atoms) then
          js=unique_atoms(is)%n_equal_atoms
       else
          js=pointcharge_array(is-N_unique_atoms)%n_equal_charges
       end if
       ks=center2sphere(is,js)
       if(ks==i_sphere) then
          center(1)=is; center(2)=js
          exit
       end if
    end do

  end function sphere2center
  !***************************************************

  !*************************************************************
  function vector_product(v1, v2)
    implicit none
    real(r8_kind), intent(in) :: v1(:), v2(:) ! (3)
    real(r8_kind) :: vector_product(3)
    !** End of interface *****************************************

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product
  !*************************************************************

  !***************************************************
  subroutine solv_poten_transfer_data()
    !
    ! Copies data from solv_cavity_module to potential_module. Does no
    ! communication.
    !
    use potential_module, only: N_points, point_in_space
    use solv_cavity_module, only: tessarea, dealloc_cavity
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: status
    integer (i4_kind) :: i, j
    N_points = n_size

    allocate (point_in_space(N_points), stat=status)
    if (status /= 0) call error_handler &
         ("solv_poten_transfer_data: allocation of point_in_space is failed")

    do i = 1, N_points
       point_in_space(i) % N_equal_points = tessarea(i) % n_equal
       allocate (point_in_space(i) % position(3, point_in_space(i) % N_equal_points), &
            stat=status)
       if (status /= 0) call error_handler &
            ("solv_poten_transfer_data: allocation of point_in_space%N_equal_points is failed")
       do j = 1, point_in_space(i) % N_equal_points
          point_in_space(i) % position(:,j) = tessarea(i) % xyz(j, :)
       enddo
    enddo

    call dealloc_cavity()       ! no comm
  end subroutine solv_poten_transfer_data
  !******************************************************

  subroutine const_part_of_ham_and_energy()
    !
    ! See do_const_part_of_ham_and_energy() below.
    ! Executed in parallel context.
    !
    use potential_module, only: N_points
    use comm, only: comm_rank, comm_bcast
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: status

    !
    ! Assuming N_points is meaningfull already here:
    !
    allocate(Q_n(N_points), stat=status)
    ASSERT(status==0)

    if ( comm_rank() == 0 ) then
        !
        ! Master does the real work alone:
        !
        call do_const_part_of_ham_and_energy()
    endif

    !
    ! Slaves get the results, this does nothing usefull in serial runs:
    !
    call comm_bcast(Q_n)
  end subroutine const_part_of_ham_and_energy

  subroutine do_const_part_of_ham_and_energy()
    !
    ! This subroutine calculates constant contributions to Hamiltonian
    ! and  energy discribing interaction  between solute  molecule and
    ! point charges destributed on cavity surface
    !
    ! 1) The contribution to energy: interaction between solute nuclei
    !    and surface charges due to solute nuclei
    !
    ! 2) The contributions to Hamiltonian:
    !
    !   a) J-matrix (interaction between electronic part of the solute
    !     potential and point  charges due to the nuclear  part of the
    !     charge distribution of the solute molecule)
    !
    !   b) Y-matrix  (interaction between  nuclear part of  the solute
    !      potential and  point charges due to the  electronic part of
    !      the charge distribution of the solute molecule)
    !
    ! because of both J and Y matrices have to be the same the program
    ! use only J matrix.
    !
    use solv_cavity_module, only: e => dielectric_constant ! and more
    use potential_module, only: N_points, point_in_space, V_pot_n, &
         V_pot_pc, get_poten_n, get_poten_pc
#ifdef WITH_EFP
    use potential_module, only: V_pot_mp
#endif
    use unique_atom_module, only: unique_atoms
    use unique_atom_methods, only: core_charge => unique_atom_core_charge
    use pointcharge_module, only: pointcharge_N,pointcharge_array
    implicit none
    ! *** end of interface ***

    real (r8_kind) :: V_buf(N_points) ! assuming N_points is defined on entry
    real (r8_kind) :: Q_n_sum, Z_sum
    integer (i4_kind) :: n_eql
    integer (i4_kind) :: i

    !
    ! 1) The contribution to  energy: interaction between solute cores
    !    and surface charges due to solute cores.
    !
    ! This fills V_pot_n(:) with  potentials of the atomic cores.  For
    ! ECPs the charge of the nuclei is shielded by the core electrons.
    !
    call get_poten_n()               ! get V_pot_n
    if (with_pc) call get_poten_pc() ! get V_pot_pc

    !
    ! V_pot_n and V_pot_e are point_in_space(i) % N_equal_points times
    ! the  potential, because:  The integrals  contain  these factors,
    ! because the totally symmetric sum of all equivalent space points
    ! is used as operand.
    !

    do i = 1, N_points
       n_eql = point_in_space(i) % N_equal_points

       V_buf(i) = V_pot_n(i) / n_eql

       if (with_pc) then
          V_buf(i) = V_buf(i) + V_pot_pc(i) / n_eql
       endif

#ifdef WITH_EFP
       if (efp .and. n_efp > 0) then
          V_buf(i) = V_buf(i) + V_pot_mp(i) / n_eql
       endif
#endif
    enddo

    ! Calculation  of point  charges due  to the  nuclear part  of the
    ! solute molecule:
    Q_n = -((e - 1) / e) * MATMUL (A_matrix_inv, V_buf)

    ! Correction of point charges respect to Gauss equation. Count the
    ! total charge of the atomic  cores shielded by core electrons (in
    ! case of ECPs):
    Z_sum = dot_product (core_charge (unique_atoms(:)), unique_atoms(:) % N_equal_atoms)

    if (with_pc) then
       do i = 1, pointcharge_N
           Z_sum = Z_sum + pointcharge_array(i) % N_equal_charges * pointcharge_array(i) % Z
       end do
    end if

    Z_sum = -Z_sum * ((e - 1) / e)

    Q_n_sum = 0.0
    do i = 1, N_points
       n_eql = point_in_space(i) % N_equal_points
       Q_n_sum = Q_n_sum + n_eql * Q_n(i)
    enddo

    write (output_unit,*) 'SOLVENT EFFECT: Discretization error of nuclear part'
    write (output_unit,*) '(Q_n-Z_n)(e-1)/e', Q_n_sum - Z_sum
    write (output_unit,*) 'Q_n*(e-1)/e', Q_n_sum, 'Z_n*(e-1)/e', Z_sum

    ! Issue a warning in any case ...
    if ((Z_sum / Q_n_sum - 1)**2 > 1.e-4_r8_kind) then
       call write_to_output_units ("solvation module: warning: large &
            & discretization error. It is recommended to choose smaller &
            & tessera size")
    endif

#ifdef WITH_EFP
    ! FIXME: this changes a variable in another module:
    if (efp .and. n_efp > 0) cor_nuc = .false.
#endif

    ! Apply correction factor only if requested. For gradients (MF):
    if (cor_nuc) then
       correct_factor_n = Z_sum / Q_n_sum
    else
       correct_factor_n = 1.0
    endif

    Q_n = correct_factor_n * Q_n

    !      1
    ! E = --- <Q * V >
    !  n   2    n   n
    energy_core_n = 0.0
    do i = 1, N_points
       energy_core_n = energy_core_n + Q_n(i) * V_pot_n(i) / 2

       if (with_pc) then
          energy_core_n = energy_core_n + Q_n(i) * V_pot_pc(i) / 2
       endif

#ifdef WITH_EFP
       if (efp .and. n_efp > 0) then
          energy_core_n = energy_core_n + Q_n(i) * V_pot_mp(i) / 2
       endif
#endif
    enddo

    write (output_unit,*) 'CONSTANT PART OF SOLVATION ENERGY'
    write (output_unit,*) '(included in energy, not in spectra)'
    write (output_unit,*) energy_core_n, 'Hartree'
    write (*,*) 'CONSTANT PART OF SOLVATION ENERGY'
    write (*,*) '(included in energy, not in spectra)'
    write (*,*) energy_core_n, 'Hartree'
  end subroutine do_const_part_of_ham_and_energy
  !**********************************************************

  !******************************************************
  subroutine calc_Q_e()
    !
    ! Computes  the charges  on the  surface  of the  cavity from  the
    ! values of electrostatic  potential there. Called from main_scf()
    ! and starts on on all workers, but does not do the same on all of
    ! them.
    !
    use solv_cavity_module, only: e => dielectric_constant ! and more
    use potential_module, only: point_in_space, N_points, V_pot_e, &
         start_read_poten_e
    use occupation_module, only: get_n_elec
    use solv_charge_mixing_module, only: mix_charges
    use comm, only: comm_rank
    implicit none
    !** End of interface *****************************************

    real (r8_kind) :: V_buf(N_points)
    real (r8_kind) :: Q_e_sum, ne_sum , correct_factor_e
    real (r8_kind) :: N_electrons
    integer (i4_kind) :: status, n_eql
    integer (i4_kind) :: i

    DPRINT 'slv:build_solv_ham: eneterd'

    ! To be executed on all workers:
    call start_read_poten_e()

    ! FIXME: why not let all of the workers do the rest?
    if (comm_rank() /= 0) return

    call get_n_elec(N_electrons)

    do i = 1, N_points
       n_eql = point_in_space(i) % N_equal_points
       V_buf(i) = V_pot_e(i) / n_eql
    enddo

    if (.not. allocated (Q_e)) then
       allocate (Q_e(N_points), stat=status)
       ASSERT(status==0)

       if (mix_charges()) then
          allocate (Q_e_old (N_points), stat=status)
          ASSERT(status==0)
          Q_e_old = 0.0
       end if
    endif

    Q_e = -((e - 1) / e) * MATMUL (A_matrix_inv, V_buf)

    !
    ! V_buf(:) is no more used ...
    !

    ne_sum = N_electrons * ((e - 1) / e)

    Q_e_sum = 0.0
    do i = 1, N_points
       n_eql = point_in_space(i) % N_equal_points
       Q_e_sum = Q_e_sum + n_eql * Q_e(i)
    enddo
    print *, Q_e_sum - ne_sum, 'Q_e_sum-ne_sum...', Q_e_sum, ne_sum
    correct_factor_e = ne_sum / Q_e_sum

#ifdef WITH_EFP
    if (efp .and. n_efp > 0) cor_el=.false.
#endif
    if(.not.cor_el) then
        if ((correct_factor_e - correct_factor_n)**2 > 1.e-4_r8_kind) then
                call write_to_output_units ("solvation module: warning: large &
                & discretization error. It is recommended to choose smaller &
                & tessera size")
        endif
        correct_factor_e = correct_factor_n
    endif

    Q_e = correct_factor_e * Q_e

#if 0
   call fix_solvq()

   contains

     subroutine fix_solvq()
       integer(kind=i4_kind):: fix_q_unit
       logical:: fix_q

       inquire(file=trim(inpfile("fix_solvq.dat")), exist=fix_q)
       if(fix_q) then
          print*, 'solv q fixed'
          fix_q_unit=openget_iounit(file=trim(inpfile('fix_solvq.dat')), &
               form='FORMATTED',status='unknown')

          read(fix_q_unit,*) Q_n,Q_e
          ! Q_e = -Q_n
          call returnclose_iounit(fix_q_unit,status='keep')
       endif
     end subroutine fix_solvq
#endif
  end subroutine calc_Q_e
  !*********************************************************

  subroutine do_surface_charge_moments (q, pts, charge, dipole)
    !
    ! Returns the total  charge and the dipole moment  of the apparent
    ! surface charge distribution.
    !
    use potential_module, only: spacepoint_type
    implicit none
    real (r8_kind), intent (in) :: q(:)           ! (n_points)
    type (spacepoint_type), intent (in) :: pts(:) ! (n_points)
    real (r8_kind), intent (out) :: charge, dipole(3)
    ! *** end of interface ***

    integer :: i

    ASSERT(size(q)==size(pts))

    charge = 0.0
    dipole = 0.0
    do i = 1, size (q)
       ! Second  axis of  %position(:, :)  is  for symmetry-equivalent
       ! points, all of them having the same charge:
       charge = charge + q(i) * size (pts(i) % position(:, :), 2)
       dipole = dipole + q(i) * sum (pts(i) % position(:, :), 2)
    enddo
  end subroutine do_surface_charge_moments


  subroutine surface_charge_moments (charge, dipole)
    !
    ! Lazy version that fetches the arguments from the modules.
    !
    use potential_module, only: point_in_space
    implicit none
    real (r8_kind), intent (out) :: charge, dipole(3)
    ! *** end of interface ***

    if (allocated (Q_e) .and. allocated (point_in_space)) then
       call do_surface_charge_moments (Q_e + Q_n, point_in_space, charge, dipole)
    else
       ! This  calculation  is  probably  without PCM.  This  is  true
       ! nevertheless:
       charge = 0.0
       dipole = 0.0
    endif
  end subroutine surface_charge_moments

  !******************************************************
  subroutine charge_mix_wrapper (scf_iter, first_iter)
    !
    ! Called on  master from build_solv_ham().  Slaves do  not seem to
    ! execute this sub.
    !
    use solv_charge_mixing_module, only: solv_charge_mix, mix_charges
    use potential_module, only: N_points
    implicit none
    integer (i4_kind), intent(in) :: scf_iter, first_iter
    !** End of interface *****************************************

    if (mix_charges() .and. n_Q_update == 1) then
       call solv_charge_mix (Q_e, Q_e_old, N_points, scf_iter, first_iter)
    endif
  end subroutine charge_mix_wrapper
  !*********************************************************

  !******************************************************
  subroutine build_solv_ham (loop, first_loop)
    !
    ! Setup  the part  of  hamiltonian due  to the  electron-dependent
    ! electrostatic interactions between solute and solvent.
    !
    ! Runs in parallel context. Called from main_scf().
    !
    use comm, only: comm_rank
    implicit none
    integer (i4_kind), intent (in) :: loop, first_loop
    !** End of interface *****************************************

    ! Was previousely called from main_scf() on master only:
    if (comm_rank() == 0) then
       call charge_mix_wrapper (loop, first_loop)
    endif

    ! Broadcast Q_e, Q_id, etc. to all workers:
    call send_receive_Q_e()

    call do_build_solv_ham()
  end subroutine build_solv_ham
  !******************************************************

  !******************************************************
  subroutine alloc_ham_solv
   !** End of interface *****************************************

    use symmetry_data_module  ! provide ssym

    integer(kind=i4_kind), allocatable :: dim_irrep(:)
    integer(kind=i4_kind)       :: n_irrep,status
    integer(kind=i4_kind)       :: n,i

    n_irrep = symmetry_data_n_irreps()

    allocate(dim_irrep(n_irrep),stat=status)
    if ( status .ne. 0) call error_handler( &
         "alloc_ham_solv: allocated dim_irrep failed" )

    do n=1,n_irrep
       dim_irrep(n) = symmetry_data_dimension(n)
    enddo

    allocate (ham_solv_el(n_irrep),STAT=status)
    if (status.ne.0) call error_handler &
         ("alloc_ham_solv: allocation of ham_solv_el failed ")

    do i=1,n_irrep
       allocate( ham_solv_el(i)%m(dim_irrep(i),dim_irrep(i)),STAT=status)
       if (status.ne.0) call error_handler &
            ("alloc_ham_solv: allocation of ham_solv_el(1) failed")

       ham_solv_el(i)%m=0.0_r8_kind
    enddo

    deallocate(dim_irrep,stat=status)
    if ( status .ne. 0) call error_handler( &
         "alloc_ham_solv: deallocated dim_irrep failed" )

  end subroutine alloc_ham_solv
  !******************************************************

  !******************************************************
  subroutine dealloc_ham_solv()
    !
    ! Does no communication. Idempotent.
    !
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: i, status

    ! To make it idempotent:
    if (.not. allocated (ham_solv_el)) return

#if 0
       if(.not.allocated(ham_solv_el_keep)) &
       allocate( ham_solv_el_keep(size(ham_solv_el)))
#endif
    do i = 1, size (ham_solv_el) ! n_irrep
#if 0
       if(.not.associated(ham_solv_el_keep(i)%m)) &
       allocate(ham_solv_el_keep(i)%m(size(ham_solv_el(i)%m,1),size(ham_solv_el(i)%m,2)))
       ham_solv_el_keep(i)%m=ham_solv_el(i)%m
    print*,'scf ham_solv',sum(ham_solv_el_keep(i)%m)
#endif
       deallocate (ham_solv_el(i) %m, STAT=status)
       if (status.ne.0) call error_handler &
            ("dealloc_ham_solv: deallocation of ham_solv_el(1) failed")
    enddo

    deallocate (ham_solv_el, STAT=status)
    if (status.ne.0) call error_handler &
         ("dealloc_ham_solv: deallocation of ham_solv_el failed ")
  end subroutine dealloc_ham_solv
  !******************************************************

  !******************************************************
  subroutine send_receive_Q_e ()
    !
    ! Broadcast  Q_e, Q_id,  and  Q_id1.  Slaves  will allocate  them,
    ! unless they already have been allocated.
    !
    use comm, only: comm_bcast
    use potential_module, only: N_points
#ifdef WITH_EFP
    use solv_cavity_module, only: Q_id, Q_id1
    use induced_dipoles_module, only: do_pol_pcm
#endif
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: status

    ! FIXME: The storage on slaves is not allocated?
    if (.not. allocated (Q_e)) then
       allocate (Q_e(N_points), stat=status)
       ASSERT(status==0)
    endif

#ifdef WITH_EFP
       if (do_pol_pcm) then
          if (.not. allocated (Q_id)) then
             allocate (Q_id(N_points), Q_id1(N_points), stat=status)
             ASSERT(status==0)
          endif
       endif
#endif

       ! Master sends, slaves receive ...
       ASSERT(size(Q_e)==N_points)
       call comm_bcast (Q_e)

#ifdef WITH_EFP
       if (do_pol_pcm) then
          ASSERT(size(Q_id)==N_points)
          call comm_bcast (Q_id)

          ASSERT(size(Q_id1)==N_points)
          call comm_bcast (Q_id1)
       endif
#endif
  end subroutine send_receive_Q_e
  !******************************************************

#if 0
  !******************************************************
  subroutine dealloc_Q
   !** End of interface *****************************************

    integer(kind=i4_kind) :: status

    if(allocated(Q_n)) then
       deallocate(Q_n,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_e)) then
       deallocate(Q_e,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_mp)) then
       deallocate(Q_n,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_id)) then
       deallocate(Q_id,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_id1)) then
       deallocate(Q_id1,stat=status)
       ASSERT(status == 0)
    endif

  end subroutine dealloc_Q
  !******************************************************
#endif

  !******************************************************
  subroutine do_build_solv_ham ()
    !
    ! The  contribution  to  the  Hamiltonian due  to  interaction  of
    ! electrons  with the  surface  charges induced  by the  electrons
    ! themselves.
    !
    ! Runs in parallel context.
    !
    use options_module, only: options_integrals_on_file
    use potential_module, only: poten_bounds, get_bounds_poten, &
         poten_integral_open, poten_integral_close
    use readwriteblocked_module
    use symmetry_data_module, only: symmetry_data_n_irreps, &
         symmetry_data_dimension
    use comm, only: comm_rank, comm_reduce
    use integralstore_module, only: integralstore_3c_poten
#ifdef WITH_EFP
    use solv_cavity_module, only: Q_id
    use induced_dipoles_module, only: do_pol_pcm
#endif
    implicit none
    !** End of interface *****************************************

    type(readwriteblocked_tapehandle) :: th_poten
    type(poten_bounds)      :: bounds
    real(kind=r8_kind),pointer :: poten_int(:)
    real(kind=r8_kind), allocatable :: sum_x(:)
    integer(kind=i4_kind), allocatable :: dim_irrep(:)
    integer (i4_kind) :: item_arr_poten, my_ind, n_irrep, i_gamma, &
         first_i
    logical ::  integrals_on_file
    integer(kind=i4_kind) :: n, status, i_mn, m, n_mn, k, i, i_meta, i_last

    integrals_on_file=options_integrals_on_file()

    n_irrep = symmetry_data_n_irreps()

    allocate(dim_irrep(n_irrep),stat=status)
    if ( status .ne. 0) call error_handler( &
         "do_build_solv_ham: allocated dim_irrep failed" )

    do n=1,n_irrep
       dim_irrep(n) = symmetry_data_dimension(n)
    enddo

    call get_bounds_poten(bounds)

    my_ind = 1 + comm_rank()
    item_arr_poten=bounds%item_arr(my_ind)

    first_i=1
    do i=1,my_ind-1
       first_i=first_i+bounds%item_arr(i)
    enddo

    call alloc_ham_solv()

    if ( integrals_on_file ) then
       ! open the integral files ----------------------
       if (item_arr_poten /= 0) call poten_integral_open(th_poten)
    else
       i_meta=1
    endif

#ifndef _VECTOR
    if (integrals_on_file) then
       allocate( poten_int(item_arr_poten),STAT=status)
       if(status.ne.0) call error_handler&
            ("do_build_solv_ham: allocation of poten_int, poten_q failed")
    endif
#endif
    i_gamma_lab: do i_gamma = 1,n_irrep

       n_mn = ( dim_irrep(i_gamma) * (dim_irrep(i_gamma) + 1) ) / 2

       if (item_arr_poten /= 0) then
          allocate( sum_x(n_mn), STAT=status )
          if (status.ne.0) call error_handler &
               ("do_build_solv_ham: allocation of sum_x failed")
          sum_x = 0.0_r8_kind

#ifndef _VECTOR
          i_mn = 1
          do m=1,dim_irrep(i_gamma)
             do n=1,m
                ! read in integrals file
                if (integrals_on_file) then
                   call readwriteblocked_read(poten_int,th_poten)
                else
                   i_last=i_meta+item_arr_poten-1
                   poten_int => integralstore_3c_poten(i_meta:i_last)
                   i_meta=i_last+1
                endif
                do k=1,item_arr_poten
                   sum_x(i_mn) =sum_x(i_mn)-poten_int(k)*Q_e(k+first_i-1)
                   sum_x(i_mn) =sum_x(i_mn)-poten_int(k)*Q_n(k+first_i-1)
#ifdef WITH_EFP
                   if(do_pol_pcm) &
                        sum_x(i_mn) =sum_x(i_mn)+poten_int(k)*Q_id(k+first_i-1) !???
#endif
                enddo
                i_mn = i_mn + 1
             enddo
          enddo
#else
          if ( integrals_on_file ) then
             allocate( poten_int(item_arr_poten*n_mn),STAT=status )
             if(status.ne.0) call error_handler &
                  ("potential_module: allocation of poten_int failed")
             call readwriteblocked_read(poten_int,th_poten)
             do k=1,item_arr_poten
                do i_mn = 1, n_mn
                   sum_x(i_mn) =sum_x(i_mn)-poten_int(k+(i_mn-1)*item_arr_poten)*Q_e(k+first_i-1)
                   sum_x(i_mn) =sum_x(i_mn)-poten_int(k+(i_mn-1)*item_arr_poten)*Q_n(k+first_i-1)
                enddo
             enddo
             deallocate(poten_int,STAT=status)
             if(status.ne.0) call error_handler &
                  ("potential_modul: deallocation of poten_int failed")
          else
             do k=1,item_arr_poten
                do i_mn = 1, n_mn
                   sum_x(i_mn) =sum_x(i_mn)-integralstore_3c_poten(i_meta-1+k+(i_mn-1)*item_arr_poten)*Q_e(k+first_i-1)
                   sum_x(i_mn) =sum_x(i_mn)-integralstore_3c_poten(i_meta-1+k+(i_mn-1)*item_arr_poten)*Q_n(k+first_i-1)
                enddo
             enddo
             i_meta=i_meta+item_arr_poten*n_mn
          endif
#endif

          i_mn = 1
          do m=1,dim_irrep(i_gamma)
             do n=1,m-1
                ham_solv_el(i_gamma)%m(m,n) = sum_x(i_mn)
                ham_solv_el(i_gamma)%m(n,m) = sum_x(i_mn)
                i_mn = i_mn + 1
             enddo
             ham_solv_el(i_gamma)%m(m,m) = sum_x(i_mn)
             i_mn = i_mn + 1
          enddo

          deallocate( sum_x, STAT=status )
          if (status.ne.0) call error_handler &
               ("J_Y_part_of_ham: deallocation of sum_j or sum_y failed")
       endif

       ! Reducing Hamiltonian on master:
       call comm_reduce (ham_solv_el(i_gamma) % m)
    enddo i_gamma_lab

#ifndef _VECTOR
    if (integrals_on_file) then
       deallocate(poten_int,STAT=status)
       if(status.ne.0) call error_handler&
            ( "do_build_solv_ham: deallocation of poten_int failed")
    endif
#endif

    if ( integrals_on_file ) then
       ! close the integral files ----------------------
       if (item_arr_poten /= 0) call poten_integral_close(th_poten)
    endif

    ! FIXME:  master does  it later  after  ham_calc_main() presumably
    ! added all cotributions to  Fock matrix, see main_scf(). It looks
    ! like the slaves try to free O(N^2) storage as soon as possible:
    if (my_ind /= 1) then
       call dealloc_ham_solv()  ! no comm, idempotent
    endif

    deallocate(dim_irrep,stat=status)
    if ( status .ne. 0) call error_handler( &
         "do_build_solv_ham: deallocated dim_irrep failed" )
  end subroutine do_build_solv_ham
  !**************************************************

  !**************************************************
  subroutine solv_energy_el()
    !
    ! Addition  of  different electrostatic  parts  of  energy due  to
    ! solvent polarization.   Sets module global  var "energy_solv_el"
    ! to this expression (in the most common case):
    !
    !        1
    !   E = --- <(Q + Q ) * (V + V)>
    !        2     e   n      e   n
    !
    ! This  is a  HALF of  the interaction  energy of  the  solute and
    ! solvent charge distributions. Note that the nn-term
    !
    !        1
    !       --- <Q * V >
    !        2    n   n
    !
    ! stored  in  the  module  private  variable  "energy_core_n"  was
    ! pre-computed earlier by do_const_part_of_ham_and_energy().
    !
    ! Does no communication.  Is executed  on all workers, but see the
    ! body. Historically was executed by master only.
    !
    use comm, only: comm_rank
    use potential_module, only: N_points, V_pot_e, V_pot_n, V_pot_pc
#ifdef WITH_EFP
    use solv_cavity_module, only: Q_id
    use potential_module, only: V_pot_mp
    use induced_dipoles_module, only: do_pol_pcm
#endif
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: i

    ! FIXME: do slaves have all the data to perform computation?
    if (comm_rank() /= 0) return

    energy_solv_el = 0.0
    do i = 1, N_points
       energy_solv_el = energy_solv_el + (Q_e(i) * V_pot_n(i) + Q_e(i) * V_pot_e(i) + &
            Q_n(i) * V_pot_e(i)) * 0.5_r8_kind

       if (with_pc) then
          energy_solv_el = energy_solv_el + Q_e(i) * V_pot_pc(i) * 0.5_r8_kind
       endif

#ifdef WITH_EFP
       if (efp .and. n_efp > 0) then
          energy_solv_el = energy_solv_el + Q_e(i) * V_pot_mp(i) * 0.5_r8_kind

          if (do_pol_pcm) then
             energy_solv_el = energy_solv_el + &
                  Q_id(i) * (V_pot_n(i) + V_pot_e(i) + V_pot_mp(i)) * 0.5_r8_kind
          endif
       end if
#endif
    enddo

    energy_solv_el = energy_solv_el + energy_core_n

#ifdef WITH_EFP
    if (.not. do_pol_pcm) then
       DPRINT 'ASC ', sum (Q_e) + sum (Q_n)
    else
       DPRINT 'ASC ', sum (Q_e) + sum (Q_n) + sum (Q_id)
    endif
#endif
  end subroutine solv_energy_el
  !**************************************************

  !**************************************************
  subroutine nuc_grad (grad_index, keep_tes_grads)
    !
    ! Calculates  and adds the  part of  gradient (induced  charges) *
    ! (gradient of electrostatic potential due to nuclei)
    !
    use unique_atom_module
    use unique_atom_methods, only: core_charge => unique_atom_core_charge
    use pointcharge_module, only : pointcharge_array, pointcharge_N
    use elec_static_field_module, only: surf_points_grad_info,surf_points_grad_index, &
         totalsym_field
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind), intent(in) :: grad_index(N_moving_unique_atoms + 1)
    logical :: keep_tes_grads

    real (r8_kind) :: za, dist, Q, QQ
    integer(kind=i4_kind)       :: ma, ma1, na, na1, index, grad_dim, j, k, l, m, i1, ea, nn, kk
    integer(kind=i4_kind)       :: n_equal_p,n_equal_a,ism
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:),rotmat(:,:)
    real(kind=r8_kind) :: dr_n(3,3),dr_ts(3,3),dxyz(3,3)
    real (r8_kind), parameter :: unit(3, 3) = reshape ([1, 0, 0, &
                                                        0, 1, 0, &
                                                        0, 0, 1], &
                                                       [3, 3])
    integer :: N_pc,N_pc1

    ! FIXME: sanity check, remove when confident:
    if (.not. pseudopot_present) then
ASSERT(all(unique_atoms(:)%zc==0))
    endif

    N_pc=0; N_pc1=0
    if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
    if(with_pc) N_pc1=pointcharge_N

    unique1: do ma=1,N_moving_unique_atoms+N_pc
       if(ma <= N_moving_unique_atoms) then
          na = moving_unique_atom_index(ma)
          ea = unique_atoms(na)%n_equal_atoms
          na1=na
       else
          na=ma-N_moving_unique_atoms
          ea=pointcharge_array(na)%N_equal_charges
          na1=na+N_unique_atoms
       end if

       unique3: do j=1,to_calc_grads%n_points
          xb => to_calc_grads%xyz(j,:,:)
          n_equal_p=to_calc_grads%n_equal(j)
          Q=to_calc_grads%Q(j)
          QQ=Q

          unique4: do k=1,n_equal_p
             ism=to_calc_grads%i_symm_sort(j,k)
             dxyz(1,:)=to_calc_grads%dxyz_totsyms(1,ma)%m(:,ism)
             dxyz(2,:)=to_calc_grads%dxyz_totsyms(2,ma)%m(:,ism)
             dxyz(3,:)=to_calc_grads%dxyz_totsyms(3,ma)%m(:,ism)

             unique5: do l=1,N_unique_atoms+N_pc1
                if(l <= N_unique_atoms) then
                   n_equal_a=unique_atoms(l)%n_equal_atoms
                   za = core_charge (unique_atoms(l))
                   xa => unique_atoms(l)%position
                else
                   n_equal_a=pointcharge_array(l-N_unique_atoms)%N_equal_charges
                   za = pointcharge_array(l-N_unique_atoms)%Z
                   xa => pointcharge_array(l-N_unique_atoms)%position
                end if

                unique6: do m=1,n_equal_a
                   if(ma <= N_moving_unique_atoms) then
                      grad_dim=grad_index(ma+1)-grad_index(ma)
                      rotmat=>unique_atom_grad_info(ma)%m(:,:,1)
                      index=grad_index(ma)
                   else
                      ma1=ma-N_moving_unique_atoms
                      grad_dim=surf_points_grad_index(ma1+1)-surf_points_grad_index(ma1)
                      rotmat=>surf_points_grad_info(ma1)%m(:,:,1)
                      index=surf_points_grad_index(ma1)
                   end if

                   dr_n = 0.0_r8_kind
                   if(l == na1 .and. m == 1) dr_n=unit
                   dr_ts=0.0_r8_kind
                   do nn=1,grad_dim
                      do kk=1,3
                         dr_ts(nn,kk)=dr_ts(nn,kk)+dot_product(rotmat(nn,:),dr_n(:,kk))*ea
                      end do
                   end do

                   if(sum((dr_ts-dxyz)*(dr_ts-dxyz)) < 1.0e-12_r8_kind) cycle unique6

                   dist = sqrt(sum((xa(:,m)-xb(k,:))**2))

                   if(ma <= N_moving_unique_atoms) then
                      do i1=1,grad_dim
                         grad_solv_totsym(index) = grad_solv_totsym(index) - &
                              (za*QQ/dist**3)*dot_product((xa(:,m)-xb(k,:)), &
                              (dr_ts(i1,:) - to_calc_grads%dxyz_totsyms(i1,ma)%m(:,ism)))
                         if(keep_tes_grads) then
                            grad_solv_totsym_tes(index,j) = grad_solv_totsym_tes(index,j) - &
                                 (za/dist**3)*dot_product((xa(:,m)-xb(k,:)), &
                                 (dr_ts(i1,:) - to_calc_grads%dxyz_totsyms(i1,ma)%m(:,ism)))
                         end if
                         index=index+1
                      end do
                   else
                      do i1=1,grad_dim
                         totalsym_field(index) = totalsym_field(index) + &
                              (za*QQ/dist**3)*dot_product((xa(:,m)-xb(k,:)), &
                              (dr_ts(i1,:) - to_calc_grads%dxyz_totsyms(i1,ma)%m(:,ism)))
                         index=index+1
                      end do
                   endif
                enddo unique6
             enddo unique5
          enddo unique4
       enddo unique3
    enddo unique1

  end subroutine nuc_grad
  !****************************************************************

  !****************************************************************
  subroutine nuc_grad_vtn (grad_index)
    use unique_atom_module
    use unique_atom_methods, only: core_charge => unique_atom_core_charge
#ifdef WITH_EFP
    use pointcharge_module, only : pointcharge_array, pointcharge_N, unique_pc_grad_info, &
         unique_pc_index, rcm
    use point_dqo_module, only: pd_array, N_pd, pq_array, N_pq
    use induced_dipoles_module, only: ipd_array, N_ipd
    use efp_solv_grad_module, only: dV_dr, dV_dtet, gradient_mpole_solv_totalsym, torque_mpole_solv_totalsym
#endif
    integer(kind=i4_kind), intent(in) :: grad_index(N_moving_unique_atoms + 1)
    !** End of interface *****************************************


    real (r8_kind) :: dist, Q, QQ
    real (r8_kind) :: Z
    real(kind=r8_kind)          :: grad1(3)
    integer(i4_kind)            :: N_efp_centers
#ifdef WITH_EFP
    real(kind=r8_kind)          :: QU(3,3)
    real(kind=r8_kind)          :: D1(3), D2(3)
    real(kind=r8_kind)          :: grad2(3)
    real(kind=r8_kind)          :: torque(3)
    real(r8_kind)               :: Q1,r_rcm(3)
    integer(i4_kind)            :: grp_c1,center_t(2),grp_t,grp_c2,grp_efp
#endif
    integer(kind=i4_kind)       :: ma, na, index, grad_dim, i, j, k, l, ll, m, i1, ea
    integer(kind=i4_kind)       :: n_equal_p,n_equal_a
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:),rotmat(:,:)
    real(kind=r8_kind) :: gradient(3)
    integer :: sphere_c,sphere_t
    !------------ Executable code ------------------------------------

    ! FIXME: sanity check, remove when confident:
    if (.not. pseudopot_present) then
ASSERT(all(unique_atoms(:)%zc==0))
    endif

    N_efp_centers=0
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) N_efp_centers=pointcharge_N+N_ipd
#endif

    unique1: do ma=1,N_moving_unique_atoms
       na = moving_unique_atom_index(ma)
       ea = unique_atoms(na)%n_equal_atoms

       unique2: do i=1,ea
          gradient=0.0_r8_kind
          sphere_c=center2sphere(na,i)

          unique3: do j=1,to_calc_grads%n_points
             xb => to_calc_grads%xyz(j,:,:)
             n_equal_p=to_calc_grads%n_equal(j)
             Q=to_calc_grads%Q(j)
#ifdef WITH_EFP
             Q1=to_calc_grads%Q1(j)
             QQ=(Q+Q1)*0.5_r8_kind
#else
             QQ=Q
#endif

             unique4: do k=1,n_equal_p
                sphere_t=to_calc_grads%sphere(j,k)

                unique5: do l=1,N_unique_atoms+N_efp_centers
                   if (l <= N_unique_atoms) then
                      ll = l
                      n_equal_a = unique_atoms(l) % n_equal_atoms
                      xa => unique_atoms(l)%position
                      Z = core_charge (unique_atoms(l))
#ifdef WITH_EFP
                   elseif (l > N_unique_atoms .and. l <= N_unique_atoms+pointcharge_N) then
                      ll = l - N_unique_atoms
                      n_equal_a = pointcharge_array(ll) % n_equal_charges
                      xa => pointcharge_array(ll) % position
                      Z = pointcharge_array(ll) % Z
                   elseif (l > N_unique_atoms + pointcharge_N) then
                      ll = l - (N_unique_atoms + pointcharge_N)
                      n_equal_a = ipd_array(ll) % n_equal_dipoles
                      xa => ipd_array(ll) % position
                      ! FIXME: and what is Z in this case?
#endif
                   else
                      ll = -1
                      n_equal_a = -1
                      xa => NULL()
                      Z = -1.0
                      ABORT ("should not happen!")
                   end if

                   unique6: do m=1,n_equal_a
                      if(l <= N_unique_atoms) then
#ifdef WITH_EFP
                      elseif(l > N_unique_atoms .and. l <= N_unique_atoms+pointcharge_N) then
                         if(N_pd > 0) D1=pd_array(ll)%dipole(:,m)
                         if(N_pq > 0) QU=pq_array(ll)%quadrupole(:,:,m)
                      elseif(l > N_unique_atoms+pointcharge_N) then
                         D1=ipd_array(ll)%idipole(:,m)
                         D2=ipd_array(ll)%idipole1(:,m)
#endif
                      end if

                      if((na==l .and. i==m) .and. sphere_c==sphere_t) cycle unique6
                      if((na/=l  .or. i/=m) .and. sphere_c/=sphere_t) cycle unique6

                      dist = sqrt(sum((xa(:,m)-xb(k,:))**2))
                      if(sphere_c==sphere_t) dist=-dist
#ifdef WITH_EFP
                      grad1=dV_dr(l,xa(:,m)-xb(k,:),dist,QQ,Z,D1,QU,D2,Q,Q1)
#else
                      grad1=-(QQ*Z/dist**3)*(xa(:,m)-xb(k,:))
#endif
                      gradient=gradient+grad1

                   enddo unique6
                enddo unique5
             enddo unique4
          enddo unique3

          grad_dim=grad_index(ma+1)-grad_index(ma)
          rotmat=>unique_atom_grad_info(ma)%m(:,:,1)
          index=grad_index(ma)
          do i1=1,grad_dim
             grad_solv_totsym(index)=grad_solv_totsym(index) + &
                  sum(rotmat(i1,:)*gradient(:))
             index=index+1
          end do
       end do unique2
    enddo unique1

#ifdef WITH_EFP
    if(efp .and. n_efp > 0) then
       grp_efp=0
       unique_efp1: do ma=1,pointcharge_N
          na = ma
          ea = pointcharge_array(na)%n_equal_charges

          unique_efp2: do i=1,ea
             grp_c1=pointcharge_array(na)%group(i)

             if(grp_c1 /= grp_efp) then
                grp_efp=grp_c1
             else
                cycle unique_efp2
             end if

             gradient=0.0_r8_kind
             torque=0.0_r8_kind

             unique_efp3: do j=1,to_calc_grads%n_points
                xb => to_calc_grads%xyz(j,:,:)
                n_equal_p=to_calc_grads%n_equal(j)
                Q=to_calc_grads%Q(j)
                Q1=to_calc_grads%Q1(j)
                QQ=(Q+Q1)*0.5_r8_kind

                unique_efp4: do k=1,n_equal_p
                   sphere_t=to_calc_grads%sphere(j,k)
                   center_t=sphere2center(sphere_t)
                   if(center_t(1) <= N_unique_atoms) then
                      grp_t=-1
                   else
                      grp_t=pointcharge_array(center_t(1)-N_unique_atoms)%group(center_t(2))
                   end if

                   unique_efp5: do l=1,N_unique_atoms+pointcharge_N+N_ipd
                      if(l <= N_unique_atoms) then
                         ll=l
                         n_equal_a=unique_atoms(ll)%n_equal_atoms
                         xa => unique_atoms(ll)%position
                         Z = core_charge (unique_atoms(ll))
                      elseif(l > N_unique_atoms .and. l <= N_unique_atoms+pointcharge_N) then
                         ll=l-N_unique_atoms
                         n_equal_a=pointcharge_array(ll)%n_equal_charges
                         xa => pointcharge_array(ll)%position
                         Z = pointcharge_array(ll)%Z
                      elseif(l > N_unique_atoms+pointcharge_N) then
                         ll=l-(N_unique_atoms+pointcharge_N)
                         n_equal_a=ipd_array(ll)%n_equal_dipoles
                         xa => ipd_array(ll)%position
                      end if

                      unique_efp6: do m=1,n_equal_a
                         if(l <= N_unique_atoms) then
                            grp_c2=-1
                         elseif(l > N_unique_atoms .and. l <= N_unique_atoms+pointcharge_N) then
                            grp_c2=pointcharge_array(ll)%group(m)
                            if(N_pd > 0) D1=pd_array(ll)%dipole(:,m)
                            if(N_pq > 0) QU=pq_array(ll)%quadrupole(:,:,m)
                         elseif(l > N_unique_atoms+pointcharge_N) then
                            grp_c2=ipd_array(ll)%group(m)
                            D1=ipd_array(ll)%idipole(:,m)
                            D2=ipd_array(ll)%idipole1(:,m)
                         end if

                         dist = sqrt(sum((xa(:,m)-xb(k,:))**2))

                         if(grp_c1==grp_t .and. grp_c1==grp_c2) cycle unique_efp6
                         if(grp_c1/=grp_t .and. grp_c1/=grp_c2) cycle unique_efp6
                         if(grp_c1==grp_t .and. grp_c1/=grp_c2) then !B
                            grad1=dV_dr(l,xa(:,m)-xb(k,:),-dist,QQ,Z,D1,QU,D2,Q,Q1)
                            gradient=gradient+grad1
                            r_rcm=xb(k,:)-rcm(:,grp_c1)
                            torque=torque-vector_product(grad1,r_rcm)
                         elseif(grp_c1/=grp_t .and. grp_c1==grp_c2) then !A
                            grad2=dV_dr(l,xa(:,m)-xb(k,:),dist,QQ,Z,D1,QU,D2,Q,Q1)
                            gradient=gradient+grad2
                            r_rcm=xa(:,m)-rcm(:,grp_c1)
                            torque=torque-vector_product(grad2,r_rcm)
                            torque=torque+dV_dtet(l,xa(:,m)-xb(k,:),dist,QQ,Z,D1,QU,D2,Q,Q1)
                         end if

                      enddo unique_efp6
                   enddo unique_efp5
                enddo unique_efp4
             enddo unique_efp3

             grad_dim=unique_pc_index(na+1)-unique_pc_index(na)
             index=unique_pc_index(na)
             rotmat=>unique_pc_grad_info(na)%m(:,:,i)
             do i1=1,grad_dim
                gradient_mpole_solv_totalsym(index)=gradient_mpole_solv_totalsym(index) + &
                     sum(rotmat(i1,:)*gradient(:))
                torque_mpole_solv_totalsym(index)=torque_mpole_solv_totalsym(index) + &
                     sum(rotmat(i1,:)*torque(:))
                index=index+1
             end do
          end do unique_efp2
       enddo unique_efp1
    end if
#endif

  end subroutine nuc_grad_vtn
  !****************************************************************

  !****************************************************************
  subroutine init_forces_on_pc()
    use elec_static_field_module, only: totalsym_field

    totalsym_field=0.0_r8_kind

  end subroutine init_forces_on_pc
  !****************************************************************

  !****************************************************************
  subroutine solv_forces_on_pc()

    use elec_static_field_module, only: N_surface_points,surface_points,transform_to_cart_field
    use elec_static_field_module, only: totalsym_field

    integer :: status,i !,j

    allocate(F_solv_pc(N_surface_points),stat=status)
    if ( status .ne. 0) call error_handler( &
         "solv_forces_on_pc: allocated  F_solv_pc failed" )

    do i=1,N_surface_points
       allocate(F_solv_pc(i)%m(3,surface_points(i)%N_equal_points),stat=status)
       if (status .ne. 0) call error_handler( &
            "solv_forces_on_pc: allocated F_solv_pc%m failed" )
    end do

    call transform_to_cart_field(F_solv_pc,totalsym_field)
!!$do i=1,N_surface_points
!!$do j=1,surface_points(i)%N_equal_points
!!$write(output_unit,'(a20,3F20.12)') ' ',F_solv_pc(i)%m(:,j)
!!$write(*,'(a20,3F20.12)') ' ',F_solv_pc(i)%m(:,j)
!!$end do
!!$end do

    totalsym_field=0.0_r8_kind

  end subroutine solv_forces_on_pc
  !****************************************************************

  !****************************************************************
  subroutine dealloc_solv_pc()

    use elec_static_field_module, only: N_surface_points

    integer :: status,i

    if(allocated(F_solv_pc)) then
       do i=1,N_surface_points
          deallocate(F_solv_pc(i)%m,stat=status)
          if ( status .ne. 0) call error_handler( &
               "solv_forces_on_pc:deallocation of F_solv_pc%m is failed" )
       end do
       deallocate(F_solv_pc,stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_forces_on_pc: deallocated  F_solv_pc failed" )
    end if

  end subroutine dealloc_solv_pc
  !****************************************************************

  subroutine shutdown_solvation()
    !
    ! Executed  by all  workers  in a  parallel  context. Called  from
    ! finalize_geometry() every  geometry iteration.  Deallocates some
    ! used variables. Let us make an entry point for a full shutdown.
    !
    use comm, only: comm_rank
    use potential_module, only: deallocate_pot, bounds_free_poten, &
        destroy_poten_file, dealloc_space_points
    use unique_atom_module, only: N_unique_atoms
    implicit none
    !** End of interface *****************************************

    !
    ! Delete everything related to solvation effect
    !
    if (comm_rank() == 0) then
        ! FIXME: was historically different for master and slaves
        call deallocate_pot()

        if (N_unique_atoms > 0) call bounds_free_poten()

        call dealloc_A_inv()
    endif

    if (N_unique_atoms > 0) call destroy_poten_file()

    call dealloc_space_points()
  end subroutine shutdown_solvation

end module solv_electrostat_module

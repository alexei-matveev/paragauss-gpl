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
!**************************************************************
subroutine disp_rep_wrap()
  !
  ! Calls the routines for calculating the dispersion/repulsion energy
  ! or its gradients  and 2nd derivatives. Does not  look like it does
  ! any communication.
  !
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use solv_cavity_module
  use output_module ,only: output_cavity_data
  use iounitadmin_module, only: output_unit
  use disp_rep_module
  use unique_atom_module, only: N_unique_atoms,unique_atoms,N_moving_unique_atoms, &
       unique_atom_grad_info,moving_unique_atom_index
  use pointcharge_module, only: pointcharge_N,pointcharge_array
  use elec_static_field_module
  use integralpar_module, only: integralpar_2dervs
  use gradient_data_module, only: grad_index=>gradient_index
  use gradient_data_module, only : dervs_totsym=>dervs_totalsym
  implicit none
  !** End of interface *****************************************

  real(kind=r8_kind) :: Vsol,K_slv,R_slv
  character(len=2) :: Name_slv
  real(kind=r8_kind),pointer  :: rotmat(:,:),rotmat1(:,:)
  real(kind=r8_kind) :: gradient(3),hessian(3,3),buffer(3,3)
  integer(kind=i4_kind) :: grad_dim,index,na,ea
  integer(kind=i4_kind) :: grad_dim1,index1,na1,ea1
  integer :: N_pc,N_pc2
  integer(kind=i4_kind) :: i,status,ma,ma1,mb,mb1,mc,i1,j,k,l
  real(kind=r8_kind) , parameter :: ang_au = 0.529177249_r8_kind
  real(kind=r8_kind) , parameter :: avogadro = 6.0221367e23_r8_kind  !mol^(-1)

  Vsol=solvent_volume*1.0e24_r8_kind/(ang_au*ang_au*ang_au)
  num_dens_dr=avogadro/Vsol

  N_pc=0; N_pc2=0
  if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
  if(with_pc) N_pc2=pointcharge_N

  allocate(nm_solute_at(N_unique_atoms+N_pc2),stat=status)
  if ( status /= 0) call error_handler( &
       "disp_rep_wrap: allocation of nm_solute_at is failed")

  if(.not. allocated(skip_short)) then
     allocate(skip_short(N_unique_atoms+N_pc2),stat=status)
     if ( status /= 0) call error_handler( &
          "disp_rep_wrap: allocation of skip_short is failed")
  endif

  if(.not. use_rappe_data) then
     call initialize_data(skip_short,with_pc)
  else
     call get_rappe_data(with_pc)
     skip_short(:)=.false.
  endif

  if (.not.do_gradients) then
     ! Public vars of disp_rep_module:
     E_disp = 0.0
     E_rep = 0.0

     do i=1,ndsa
        if(output_cavity_data) write (output_unit,*) 'THE CAVITY -',i
        call get_k_and_r(i,k_slv,R_slv,Name_slv)
        if(K_slv==0.0_r8_kind .or. R_slv==0.0_r8_kind) then
           cycle
        endif

        R_access=r_sphere_at(i)
        call points_on_cavity_surface()

        call solv_disprep_transer_data()

        ! Increments E_disp, E_rep among other things. DO_GRADIENTS is
        ! false here:
        call disp_rep_energies_and_grad(slv_cycle=i, &
                                        do_grad  =do_gradients, &
                                        with_pc  =with_pc)

        call deallocate_D_R_points(do_gradients)
     enddo
  else

     do i=1,ndsa
        if(output_cavity_data) write (output_unit,*) 'THE CAVITY -',i
        call get_k_and_r(i,k_slv,R_slv,Name_slv)
        if(K_slv==0.0_r8_kind .or. R_slv==0.0_r8_kind) then
           cycle
        endif

        R_access=r_sphere_at(i)
        call points_on_cavity_surface()
        call solv_disprep_transer_data()
        unique1: do ma=1,N_moving_unique_atoms+N_pc
           if(ma <= N_moving_unique_atoms) then
              na = moving_unique_atom_index(ma)
              ea = unique_atoms(na)%n_equal_atoms
           else
              na=ma-N_moving_unique_atoms
              ea=pointcharge_array(na)%N_equal_charges
           end if
           unique2: do mb=1,ea
              gradient=0.0_r8_kind
              call disp_rep_energies_and_grad(slv_cycle=i, &
                                              do_grad  =do_gradients, &
                                              with_pc  =with_pc, &
                                              gradient =gradient, &
                                              ma       =ma, &
                                              ea       =mb)
              if(ma <= N_moving_unique_atoms) then
                 grad_dim=grad_index(ma+1)-grad_index(ma)
                 rotmat=>unique_atom_grad_info(ma)%m(:,:,mb)
                 index=grad_index(ma)
                 do i1=1,grad_dim
                    grad_solv_totsym(index) = grad_solv_totsym(index) + &
                         sum( rotmat(i1,:) * gradient(:) )
                    index=index+1
                 enddo
              else
                 mc=ma-N_moving_unique_atoms
                 grad_dim=surf_points_grad_index(mc+1)-surf_points_grad_index(mc)
                 rotmat=>surf_points_grad_info(mc)%m(:,:,mb)
                 index=surf_points_grad_index(mc)
                 do i1=1,grad_dim
                    totalsym_field(index) = totalsym_field(index) - &
                         sum( rotmat(i1,:) * gradient(:) )
                    index=index+1
                 end do
              end if

              if(integralpar_2dervs) then
                 unique1a: do ma1=1,N_moving_unique_atoms
                    na1 = moving_unique_atom_index(ma1)
                    ea1 = unique_atoms(na1)%n_equal_atoms
                    grad_dim1=grad_index(ma1+1)-grad_index(ma1)
                    unique2a: do mb1=1,ea1
                       rotmat1=>unique_atom_grad_info(ma1)%m(:,:,mb1)

                       hessian=0.0_r8_kind
                       call disp_rep_energies_and_grad(slv_cycle=i, &
                                                       do_grad  =do_gradients, &
                                                       with_pc  =with_pc, &
                                                       ma       =ma, &
                                                       ea       =mb, &
                                                       hessian  =hessian, &
                                                       ma1      =ma1, &
                                                       ea1      =mb1)

                       buffer=0.0_r8_kind
                       do j=1,grad_dim
                          do k=1,3
                             do l=1,3
                                buffer(j,l)=buffer(j,l)+rotmat(j,k)*hessian(k,l)
                             end do
                          end do
                       end do

                       index=grad_index(ma)
                       do j=1,grad_dim
                          index1=grad_index(ma1)
                          do k=1,grad_dim1
                             do l=1,3
                                dervs_totsym(index,index1)=dervs_totsym(index,index1)+ &
                                     rotmat1(k,l)*buffer(j,l)
                             end do
                             index1=index1+1
                          end do
                          index=index+1
                       end do
                    end do unique2a
                 end do unique1a
              end if
           enddo unique2
        enddo unique1

        call deallocate_D_R_points(do_gradients)
     enddo
  endif

  call deallocate_parameters(do_gradients)
  deallocate(nm_solute_at,stat=status)
  if ( status /= 0) call error_handler( &
       "disp_rep_wrap: deallocation of nm_solute_at is failed")

contains ! of disp_rep_wrap

  subroutine solv_disprep_transer_data()
    !stores data needed for calculating the dispersion/repulsion
    !energy or its gradients in variables
    !D_R_points
    !D_R_grad
    !** End of interface *****************************************

    integer(kind=i4_kind) :: status,nat,natm,eat
    integer (i4_kind) :: j, j1, k, k1, l, m, n, is
    integer(i4_kind)      :: j11,k11,n1,m1,nat1,eat1
    integer(kind=i4_kind) :: i_sym,N_pc1

    DPRINT 'solv_disprep_t_d: do_gradients=',do_gradients

    N_pc1=0
    if(with_pc) N_pc1=pointcharge_N
    N_points_dr=n_size

    allocate(D_R_points(N_points_dr),stat=status)
    ASSERT(status==0)
    if (do_gradients) then
       allocate(D_R_grad(N_points_dr),stat=status)
       ASSERT(status==0)
       if(integralpar_2dervs) then
          allocate(D_R_hess(N_points_dr),stat=status)
          ASSERT(status==0)
       end if
    endif

    ASSERT(size(tessarea)==n_size)
    ASSERT(ua_dim_max.ne.-1)

    if(do_gradients)then
       DPRINT 'solv_disprep_t_d: shape(i_symm_sort)=',shape(i_symm_sort)
    endif

    do l=1,N_points_dr
       D_R_points(l)%N_equal_points = tessarea(l)%n_equal

       allocate(D_R_points(l)%position(3,D_R_points(l)%N_equal_points),stat=status)
       ASSERT(status==0)
       allocate(D_R_points(l)%out_normal(3,D_R_points(l)%N_equal_points),stat=status)
       ASSERT(status==0)

       if(do_gradients) then
          allocate(D_R_grad(l)%position &
               (3,N_moving_unique_atoms+N_pc1,ua_dim_max,3,tessarea(l)%n_equal),&
               D_R_grad(l)%out_normal &
               (3,N_moving_unique_atoms+N_pc1,ua_dim_max,3,tessarea(l)%n_equal), &
               D_R_grad(l)%area &
               (3,N_moving_unique_atoms+N_pc1,ua_dim_max,tessarea(l)%n_equal), &
               D_R_grad(l)%distance &
               (3,N_moving_unique_atoms+N_pc1,ua_dim_max,N_unique_atoms+N_pc1,ua_dim_max, &
               3,tessarea(l)%n_equal), &
               stat=status)
          ASSERT(status==0)
          if(integralpar_2dervs) then
             allocate(D_R_hess(l)%area(tessarea(l)%n_equal), &
                  D_R_hess(l)%out_normal(3,tessarea(l)%n_equal), &
                  D_R_hess(l)%distance(3,tessarea(l)%n_equal), &
                  stat=status)
             ASSERT(status==0)
          end if
       endif

       do j=1,D_R_points(l)%N_equal_points

          if(do_gradients .and. integralpar_2dervs) then
             allocate(D_R_hess(l)%area(j)%m(3,N_moving_unique_atoms,ua_dim_max, &
                  3,N_moving_unique_atoms,ua_dim_max),stat=status)
             ASSERT(status==0)
             do is=1,3
                allocate(D_R_hess(l)%out_normal(is,j)%m(3,N_moving_unique_atoms,ua_dim_max, &
                     3,N_moving_unique_atoms,ua_dim_max),stat=status)
                ASSERT(status==0)
                allocate(D_R_hess(l)%distance(is,j)%m(3,N_moving_unique_atoms,ua_dim_max, &
                     3,N_moving_unique_atoms,ua_dim_max),stat=status)
                ASSERT(status==0)
             end do
          end if

          D_R_points(l)%position(:,j) = tessarea(l)%xyz(j,:)
          k=tessarea(l)%sphere(j)
          D_R_points(l)%out_normal(:,j)=(tessarea(l)%xyz(j,:)-xyz_sphere(k,:))/ &
               sqrt(dot_product(tessarea(l)%xyz(j,:)-xyz_sphere(k,:),tessarea(l)%xyz(j,:)-xyz_sphere(k,:)))

          if(do_gradients) then
             i_sym = i_symm_sort(l,j)
             do j1=1,N_moving_unique_atoms+N_pc1 ! gradient coodinates
                if(j1 <= N_moving_unique_atoms) then
                   nat = moving_unique_atom_index(j1)
                   eat = unique_atoms(nat)%n_equal_atoms
                   natm=nat
                else
                   nat=j1-N_moving_unique_atoms
                   eat=pointcharge_array(nat)%N_equal_charges
                   natm=nat+N_unique_atoms
                end if
                do k1=1,eat
                   do n=1,3
                      ! SIGSEGV here:
                      D_R_grad(l)%position(n,j1,k1,:,j)=&
                           cagr%dcenter(n,j1,k1)%m(:,i_sym)
                      D_R_grad(l)%area(n,j1,k1,j)=&
                           cagr%darea(n,j1,k1)%m(i_sym)

                      do m=1,3 ! coordinate(x,y,z) of point j
                         D_R_grad(l)%out_normal(n,j1,k1,m,j)=&
                              (cagr%dcenter(n,j1,k1)%m(m,i_sym)-&
                              cagr%dc(m,k)%xyz_grad(n,j1,k1))/r_sphere(k)

                         D_R_grad(l)%distance(n,j1,k1,:,:,m,j)=&
                              D_R_grad(l)%position(n,j1,k1,m,j)
                      enddo

                      D_R_grad(l)%distance(n,j1,k1,natm,k1,n,j)=&
                           D_R_grad(l)%distance(n,j1,k1,natm,k1,n,j)-1.0_r8_kind

                      if(integralpar_2dervs) then
                         do j11=1,N_moving_unique_atoms
                            nat1 = moving_unique_atom_index(j11)
                            eat1 = unique_atoms(nat1)%n_equal_atoms
                            do k11=1,eat1
                               do n1=1,3
                                  D_R_hess(l)%area(j)%m(n,j1,k1,n1,j11,k11)= &
                                       cagr%d2area(n,j1,k1,n1,j11,k11)%m(i_sym)

                                  do m1=1,3
                                     D_R_hess(l)%out_normal(m1,j)%m(n,j1,k1,n1,j11,k11)= &
                                          (cagr%d2center(n,j1,k1,n1,j11,k11)%m(m1,i_sym)-&
                                          cagr%dc(m1,k)%xyz_hess(n,j1,k1,n1,j11,k11))/r_sphere(k)
                                     D_R_hess(l)%distance(m1,j)%m(n,j1,k1,n1,j11,k11)= &
                                          cagr%d2center(n,j1,k1,n1,j11,k11)%m(m1,i_sym)
                                  end do
                               end do
                            end do
                         end do
                      end if

                   enddo
                enddo
             enddo
          endif
       enddo
       D_R_points(l)%area=tessarea(l)%area
    enddo

    call dealloc_cavity()
    if(do_gradients) then
       call dealloc_geom_deriv_part1(cagr)
       call dealloc_geom_deriv_part2()
    endif

  end subroutine solv_disprep_transer_data

end subroutine disp_rep_wrap
!**************************************************************

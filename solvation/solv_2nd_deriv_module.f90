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
!===============================================================
module solv_2nd_deriv_module
!== Interrupt of public interface of module =========
!  Author: AS
!  Date: 03/06
!
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
!== Interrupt end of public interface of module =====              

!----modules used ------------------
# include <def.h>
  use type_module ! type specification parameters
  use datatype
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public nuc_solv_2nd_deriv, charge_solv_2nd_deriv, matrix_2nd_deriv
  public cpks_solv_main, cpks_solv_Qai, cpks_solv_H1, cpks_solv_H1_1
!  public cpks_solv_ham,calc_solv_ham_mo 

!== Interrupt of public interface of module =========

  external error_handler

contains
#if 0  /* to calc solv hamiltonian */
   subroutine fix_solv_q(Q_e)
use solv_cavity_module, only: to_calc_grads
use filename_module, only: input_dir
use iounitadmin_module
   integer(kind=i4_kind):: fix_q_unit
   logical:: fix_q
   real(r8_kind), dimension(to_calc_grads%n_points):: Q_n,Q_e
    

     inquire(file=trim(input_dir)//"/fix_solvq.dat",exist=fix_q)
     if(fix_q) then
     print*, 'solv q Q_n added'
               fix_q_unit=openget_iounit(file=trim(input_dir)//'/fix_solvq.dat',&
                               form='FORMATTED',status='unknown')

           read(fix_q_unit,*) Q_n
           print*,'Q_n',sum(Q_n)
           Q_e=Q_e+Q_n
!           Q_e=-Q_n
           call returnclose_iounit(fix_q_unit,status='keep')
     endif
   end subroutine fix_solv_q
#endif

  !********************************************************************
  subroutine nuc_solv_2nd_deriv
    !Nuclear contribution to explicit solvation second 
    !derivatives (Q*V_nuc_xy) 
    use solv_cavity_module, only:to_calc_grads
    use unique_atom_module, only:N_moving_unique_atoms,moving_unique_atom_index,pseudopot_present
    use unique_atom_module, only:unique_atoms,N_unique_atoms,unique_atom_grad_info
    use gradient_data_module, only:gradient_index,dervs_totalsym 

    real(i8_kind), pointer :: xa(:,:),xt(:,:)
    real(i8_kind) :: za,zca,dist,dist3,dist5,Q,QZ
    integer(i4_kind) :: n_equal_t,n_equal_a,ism
    integer(i4_kind) :: ma,na,ea,ma1,na1,ea1,i,j,k,l,m,n,i1,n1,kk,kk1
    real(i8_kind) :: dR1(3,3),dR2(3,3),r_at(3)!,dR1_ts(3),dR2_ts(3)
    real(i8_kind) :: dR1_ts(3,3),dR2_ts(3,3)
    real(i8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    real(i8_kind) :: unit(3,3)=reshape((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind, &
                                         0.0_r8_kind,1.0_r8_kind,0.0_r8_kind, &
                                         0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/), &
                                       (/3,3/))
    integer(i4_kind) :: grad_dim,grad_dim1,index,index1

    unique1: do ma=1,N_moving_unique_atoms
       na = moving_unique_atom_index(ma)
       ea = unique_atoms(na)%n_equal_atoms
       grad_dim=gradient_index(ma+1)-gradient_index(ma)
       unique2: do i=1,1!ea
          rotmat1=>unique_atom_grad_info(ma)%m(:,:,i)

          unique1a: do ma1=1,N_moving_unique_atoms
             na1 = moving_unique_atom_index(ma1)
             ea1 = unique_atoms(na1)%n_equal_atoms
             grad_dim1=gradient_index(ma1+1)-gradient_index(ma1)

             unique2a: do i1=1,1!ea1
                rotmat2=>unique_atom_grad_info(ma1)%m(:,:,i1)

                unique3: do j=1,to_calc_grads%n_points
                   n_equal_t=to_calc_grads%n_equal(j)
                   xt => to_calc_grads%xyz(j,:,:)
                   Q=to_calc_grads%Q(j)

                   unique4: do k=1,n_equal_t
                      ism=to_calc_grads%i_symm_sort(j,k)

                      unique5: do l=1,N_unique_atoms
                         n_equal_a=unique_atoms(l)%n_equal_atoms
                         za = unique_atoms(l)%Z
                         xa => unique_atoms(l)%position
                         zca = unique_atoms(l)%ZC
                         if (.not.pseudopot_present) zca = 0.0_r8_kind
                         za=za-zca
                         QZ=Q*za

                         unique6: do m=1,1!n_equal_a
                            r_at=xa(1:3,m)-xt(k,1:3)
                            dist = sqrt(dot_product(r_at,r_at))
                            dist3 = dist*dist*dist
                            dist5 = dist3*dist*dist
                            dR1=0.0_r8_kind; dR2=0.0_r8_kind
                            if(l==ma .and. m==i) dR1=unit
                            if(l==ma1 .and. m==i1) dR2=unit
                            index=gradient_index(ma)

                            dR1_ts=0.0_r8_kind
                            do n=1,grad_dim
                               do kk=1,3
                                  dR1_ts(n,kk)=dR1_ts(n,kk)+dot_product(rotmat1(n,:),dR1(:,kk))
                               end do
                            end do
                            do n=1,grad_dim
                               index1=gradient_index(ma1)

                               dR2_ts=0.0_r8_kind
                               do n1=1,grad_dim1
                                  do kk1=1,3
                                     dR2_ts(n1,kk1)=dR2_ts(n1,kk1)+dot_product(rotmat2(n1,:),dR2(:,kk1))
                                  end do
                               end do
                               do n1=1,grad_dim1
                                  dervs_totalsym(index,index1)=dervs_totalsym(index,index1)+QZ*( &
                                       3.0_r8_kind* &
                                       dot_product(r_at,dR1_ts(n,:)- &
                                                   to_calc_grads%dxyz_totsyms(n,ma)%m(:,ism))* &
                                       dot_product(r_at,dR2_ts(n1,:)- &
                                                   to_calc_grads%dxyz_totsyms(n1,ma1)%m(:,ism))/dist5- &
                                       dot_product(dR1_ts(n,:)- &
                                                   to_calc_grads%dxyz_totsyms(n,ma)%m(:,ism), &
                                                   dR2_ts(n1,:)- &
                                                   to_calc_grads%dxyz_totsyms(n1,ma1)%m(:,ism))/dist3- &
                                       dot_product(r_at, &
                                                  -to_calc_grads%d2xyz_totsyms(n,ma,n1,ma1)%m(:,ism))/dist3)
                                  index1=index1+1
                               end do
                               index=index+1
                            end do
                         end do unique6
                      end do unique5
                   end do unique4
                end do unique3
             end do unique2a
          end do unique1a
       end do unique2
    end do unique1

  end subroutine nuc_solv_2nd_deriv
  !********************************************************************

  !********************************************************************
  subroutine matrix_2nd_deriv
    !Calculation of A_xy*q1*q2+2.0*A_x*q1_y*q2 contribution to the second derivatives
    use solv_cavity_module, only:to_calc_grads
    use solv_electrostat_module, only: Q_grad
    use unique_atom_module, only: N_moving_unique_atoms
    use gradient_data_module, only:gradient_index,dervs_totalsym 
    use comm_module, only: comm_get_n_processors,comm_myindex

    integer(i4_kind) :: ma,ma1,i,j,k,l,n,n1
    integer(i4_kind) :: grad_dim,grad_dim1,index,index1,ism,ism1
    integer(i4_kind) :: first_index,first_equal,second_index,second_equal
    real(i8_kind) :: eps_help,Qi,QiQi,Qj,QiQj,Aii,s,ss
    integer(i4_kind) :: n_proc,my_index,N_species,n_start,n_length
    real(i8_kind) :: Rij(3),Dij,Dij3,Dij5,dRij1(3),dRij2(3),d2Rij(3)
    real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind

    N_species=N_moving_unique_atoms
    n_proc=comm_get_n_processors()
    my_index=comm_myindex()
    n_start=1
    n_length=0
    do i=1,my_index
       n_start=n_start+n_length
       N_species=N_species-n_length
       n_length=int(N_species/n_proc)
       n_proc=n_proc-1
    end do

    if(n_length == 0) return

    eps_help=to_calc_grads%dielconst/(2.0_r8_kind*(to_calc_grads%dielconst-1.0_r8_kind))

    unique1: do ma=n_start,n_start+n_length-1
       grad_dim=gradient_index(ma+1)-gradient_index(ma)

       unique1a: do ma1=1,N_moving_unique_atoms
          grad_dim1=gradient_index(ma1+1)-gradient_index(ma1)

          call calc_deriv()

       end do unique1a
    end do unique1

  contains

    subroutine calc_deriv()

          first_index=0
          i_n_size: do i=1,to_calc_grads%n_points
             first_equal=to_calc_grads%n_equal(i)

             Qi=to_calc_grads%Q(i)
             QiQi=Qi*Qi
             s=to_calc_grads%s(i)
             ss=s**2
             Aii=-1.07_r8_kind*sqrt(4.0_r8_kind*pi/s)

             j_first_equal: do j=1,first_equal
                first_index=first_index+1
                ism=to_calc_grads%i_symm_sort(i,j)

                second_index=0
                k_n_size: do k=1,to_calc_grads%n_points
                   second_equal=to_calc_grads%n_equal(k)

                   Qj=to_calc_grads%Q(k)
                   QiQj=Qi*Qj

                   l_second_equal: do l=1,second_equal
                      second_index=second_index+1
                      ism1=to_calc_grads%i_symm_sort(k,l)

                      f_s: if (first_index==second_index) then
                         index=gradient_index(ma)
                         do n=1,grad_dim

                            index1=gradient_index(ma1)
                            do n1=1,grad_dim1
                               dervs_totalsym(index,index1)=dervs_totalsym(index,index1)- &
                                    QiQi*eps_help*(0.75_r8_kind*(Aii/ss)* &
                                    to_calc_grads%ds_totsyms(n,ma)%m(ism)* &
                                    to_calc_grads%ds_totsyms(n1,ma1)%m(ism)- &
                                    (Aii*0.5_r8_kind/s)* &
                                    to_calc_grads%d2s_totsyms(n,ma,n1,ma1)%m(ism))

                               dervs_totalsym(index,index1)=dervs_totalsym(index,index1)+ &
                                    2.0_r8_kind*Qi*eps_help*(Aii/(2.0_r8_kind*s))* &
                                    to_calc_grads%ds_totsyms(n,ma)%m(ism)* &
                                    Q_grad(index1,i)

                               index1=index1+1
                            end do
                            index=index+1
                         end do
                      else
                         Rij=to_calc_grads%xyz(i,j,:)-to_calc_grads%xyz(k,l,:)
                         Dij=sqrt(dot_product(Rij,Rij))
                         Dij3=Dij*Dij*Dij
                         Dij5=Dij3*Dij*Dij

                         index=gradient_index(ma)
                         do n=1,grad_dim

                            dRij1=to_calc_grads%dxyz_totsyms(n,ma)%m(:,ism)- &
                                 to_calc_grads%dxyz_totsyms(n,ma)%m(:,ism1)
                            index1=gradient_index(ma1)
                            do n1=1,grad_dim1
                               dRij2=to_calc_grads%dxyz_totsyms(n1,ma1)%m(:,ism)- &
                                    to_calc_grads%dxyz_totsyms(n1,ma1)%m(:,ism1)

                               d2Rij=to_calc_grads%d2xyz_totsyms(n,ma,n1,ma1)%m(:,ism)- &
                                    to_calc_grads%d2xyz_totsyms(n,ma,n1,ma1)%m(:,ism1)

                               dervs_totalsym(index,index1)=dervs_totalsym(index,index1)- &
                                    QiQj*eps_help*( &
                                    dot_product(dRij2,dRij1)/Dij3+ &
                                    dot_product(Rij,d2Rij)/Dij3- &
                                    3.0_r8_kind*dot_product(Rij,dRij1)*dot_product(Rij,dRij2)/Dij5)

                               dervs_totalsym(index,index1)=dervs_totalsym(index,index1)- &
                                    2.0_r8_kind*Qj*eps_help*Q_grad(index1,i)*dot_product(Rij,dRij1)/Dij3

                               index1=index1+1
                            end do
                            index=index+1
                         end do
                      end if f_s
                   end do l_second_equal
                end do k_n_size
             end do j_first_equal
          end do i_n_size
    end subroutine calc_deriv

  end subroutine matrix_2nd_deriv
  !********************************************************************

  !********************************************************************
  subroutine charge_solv_2nd_deriv
    !Charge contribution to explicit solvation second 
    !derivatives (Q_y*V_x=-A^-1*(V_y+A_y*Q)*V_x) 
    use solv_cavity_module, only: grad_solv_totsym_tes, to_calc_grads
    use solv_electrostat_module, only: A_matrix_inv, Q_grad
    use unique_atom_module, only:N_moving_unique_atoms
    use gradient_data_module, only:gradient_index,dervs_totalsym 

    type(arrmat1), allocatable :: A_matrix_inv_grad(:,:)
    real(r8_kind) :: coeff
    real(r8_kind), pointer :: Q(:)
    integer(i4_kind) :: ma,ma1,i,n,n1,j,j1,status
    integer(i4_kind) :: index,index1,grad_dim,grad_dim1,N_p

    allocate(A_matrix_inv_grad(3,N_moving_unique_atoms),stat=status)
    ASSERT(status.eq.0)

    N_p=to_calc_grads%n_points
    do i=1,N_moving_unique_atoms
       do j=1,3
          allocate(A_matrix_inv_grad(j,i)%m(N_p),stat=status)
          ASSERT(status.eq.0)
          A_matrix_inv_grad(j,i)%m=0.0_r8_kind
       end do
    end do

    Q=>to_calc_grads%Q
    call inverse_matrix_grad()

    coeff=(to_calc_grads%dielconst-1.0_r8_kind)/to_calc_grads%dielconst

    !Calculation of Q_grad
    unique_q: do ma=1,N_moving_unique_atoms
       grad_dim=gradient_index(ma+1)-gradient_index(ma)

       unique_q1: do j1=1,to_calc_grads%n_points
          index=gradient_index(ma)

          do n=1,grad_dim
             Q_grad(index,:)=Q_grad(index,:)-A_matrix_inv(:,j1)* &
                  (coeff*grad_solv_totsym_tes(index,j1)+A_matrix_inv_grad(n,ma)%m(j1))
             index=index+1
          end do
       end do unique_q1
    end do unique_q

    !Charge contribution to second derivatives
    unique1: do ma=1,N_moving_unique_atoms
       grad_dim=gradient_index(ma+1)-gradient_index(ma)

       unique1a: do ma1=1,N_moving_unique_atoms
          grad_dim1=gradient_index(ma1+1)-gradient_index(ma1)

          unique3: do j1=1,to_calc_grads%n_points
             index=gradient_index(ma)

             do n=1,grad_dim
                index1=gradient_index(ma1)

                do n1=1,grad_dim1
                   dervs_totalsym(index,index1)=dervs_totalsym(index,index1)+ &
                        Q_grad(index1,j1)*grad_solv_totsym_tes(index,j1)
                   index1=index1+1
                end do
                index=index+1
             end do
          end do unique3
       end do unique1a
    end do unique1

    do i=1,N_moving_unique_atoms
       do j=1,3
          deallocate(A_matrix_inv_grad(j,i)%m,stat=status)
          ASSERT(status.eq.0)
       end do
    end do
    deallocate(A_matrix_inv_grad,stat=status)
    ASSERT(status.eq.0)

  contains

    subroutine inverse_matrix_grad()

      real(r8_kind), allocatable :: A_matrix_grad(:,:)
      real(r8_kind) :: A_matrix_fs,vect(3),dvect(3),distance,distance3
      integer(i4_kind) :: N_pi,ii,ji,ki,li,mi,statusi
      integer(i4_kind) :: first_index,first_equal,second_index,second_equal
      integer(i4_kind) :: ism,ism1,grad_dimi
      real(r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind

      N_pi=to_calc_grads%n_points
 
      allocate(A_matrix_grad(N_pi,N_pi),stat=statusi)
      ASSERT(statusi.eq.0)

      do ii=1,N_moving_unique_atoms
         grad_dimi=gradient_index(ii+1)-gradient_index(ii)
         do ji=1,grad_dimi

            A_matrix_grad=0.0_r8_kind
            first_index=1
            do li=1,N_pi
               first_equal=to_calc_grads%n_equal(li)
               ism=to_calc_grads%i_symm_sort(li,1)

               second_index=1
               do ki=1,N_pi
                  second_equal=to_calc_grads%n_equal(ki)

                  do mi=1,second_equal
                     ism1=to_calc_grads%i_symm_sort(ki,mi)

                     f_s: if (first_index==second_index) then
                        A_matrix_fs = &
                             (1.07_r8_kind*sqrt(pi/to_calc_grads%s(li))/to_calc_grads%s(li))* &
                             to_calc_grads%ds_totsyms(ji,ii)%m(ism)
                     else
                        vect=to_calc_grads%xyz(li,1,:)-to_calc_grads%xyz(ki,mi,:)
                        distance=sqrt(dot_product(vect,vect))
                        distance3=distance*distance*distance
                        dvect=to_calc_grads%dxyz_totsyms(ji,ii)%m(:,ism)- &
                             to_calc_grads%dxyz_totsyms(ji,ii)%m(:,ism1)

                        A_matrix_fs = dot_product(vect,dvect)/distance3
                     endif f_s

                     A_matrix_grad(li,ki)=A_matrix_grad(li,ki)-A_matrix_fs
                     second_index=second_index+1
                  enddo
               enddo
               first_index=first_index+first_equal
            enddo
            A_matrix_inv_grad(ji,ii)%m=matmul(A_matrix_grad,Q)

         end do
      end do

      deallocate(A_matrix_grad,stat=statusi)
      ASSERT(statusi.eq.0)

    end subroutine inverse_matrix_grad

  end subroutine charge_solv_2nd_deriv
  !********************************************************************

  !********************************************************************
  subroutine calc_cpks_V_mo(oo_ou)
    !Transform the electrostatic potential matrix elements from
    !atomic orbital basis to molecular one. Parts of the final vector cpks_V_vec
    !are kept on corresponding hosts
    use readwriteblocked_module
    use integralstore_module, only: integralstore_3c_poten
    use options_module, only: options_integrals_on_file
    use symmetry_data_module, only: symmetry_data_n_spin,symmetry_data_n_irreps
    use symmetry_data_module, only: symmetry_data_n_partners,symmetry_data_dimension
    use comm_module, only: comm_myindex
    use potential_module, only: get_bounds_poten,poten_bounds
    use potential_module, only: poten_integral_open,poten_integral_close
    use cpksdervs_matrices, only: cpks,cpks3c
    use eigen_data_module, only : eigvec
    !** End of interface *****************************************
    character(len=2) :: oo_ou

    type(readwriteblocked_tapehandle) :: th_poten
    type(poten_bounds) :: bounds
    real(r8_kind),pointer :: poten_int(:),eigvec1(:),eigvec2(:,:,:)
    integer(i4_kind) :: item_arr_poten,my_ind,n_irrep,i_gamma,n_spin
    logical ::  spin_polarized,integrals_on_file
    integer(i4_kind) :: n,status,m,k,i_meta,i_last,i_s,n1,n2
    integer(i4_kind),allocatable :: dim_irrep(:)
    integer(i4_kind) :: dim,dim1,dim2,i
    real(r8_kind) :: mul
    logical :: ai,mn
    real(r8_kind),allocatable :: V_moao(:,:,:)
    real(r8_kind),pointer :: buffer(:,:,:)

    ai=.false.; mn=.false.
    ai=(index(trim(oo_ou),'ai')==1)
    mn=(index(trim(oo_ou),'mn')==1)

    integrals_on_file=options_integrals_on_file()

    n_spin =symmetry_data_n_spin()
    spin_polarized = n_spin > 1

    my_ind=comm_myindex()

    call get_bounds_poten(bounds)
    item_arr_poten=bounds%item_arr(my_ind)
    if (item_arr_poten == 0) then 
       return
    else
       n_irrep = symmetry_data_n_irreps()
       allocate(dim_irrep(n_irrep),stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_2nd_deriv_module:calc_cpks_V_vec: allocation of dim_irrep failed" )
       do n=1,n_irrep
          dim_irrep(n)  = symmetry_data_dimension(n)
       enddo
    endif

    do i_gamma=1,n_irrep
       do i_s=1,n_spin
          dim1=size(cpks(1,i_gamma,i_s)%Qai,1)
          if(ai) then
             dim2=size(cpks(1,i_gamma,i_s)%Qai,2)
             allocate(cpks3c(i_gamma,i_s)%V_ai(item_arr_poten,dim1,dim2),stat=status)
          else if(mn) then
             dim2=dim1
             allocate(cpks3c(i_gamma,i_s)%V_mn(item_arr_poten,dim1,dim2),stat=status)
          end if
          if ( status .ne. 0) call error_handler( &
               "solv_2nd_deriv_module:calc_cpks_V_vec: allocation of cpks3c%V_ai  failed" )
          if(ai) then
             cpks3c(i_gamma,i_s)%V_ai=0.0_r8_kind
          else if(mn) then
             cpks3c(i_gamma,i_s)%V_mn=0.0_r8_kind
          end if
       enddo
    enddo

    if (integrals_on_file) then 
       allocate( poten_int(item_arr_poten),STAT=status)
       if(status.ne.0) call error_handler&
            ("solv_2nd_deriv_module:calc_cpks_V_vec: allocation of poten_int failed")
    endif

    do i_s=1,n_spin
       if ( integrals_on_file ) then
          call poten_integral_open(th_poten)
       else
          i_meta=1
       endif
       do i_gamma = 1,n_irrep
          dim=size(eigvec(i_gamma)%m,1)
          dim1=size(cpks(1,i_gamma,i_s)%Qai,1)
          if(mn) then 
             dim2=dim1
          else if(ai) then
             dim2=size(cpks(1,i_gamma,i_s)%Qai,2)
          end if

          allocate(V_moao(item_arr_poten,dim1,dim))
          V_moao=0.0_r8_kind 

          do m=1,dim_irrep(i_gamma)
             do n=1,m
                mul=2.0_r8_kind 
                if (n==m) mul=1.0_r8_kind

                if (integrals_on_file) then
                   call readwriteblocked_read(poten_int,th_poten)
                else
                   i_last=i_meta+item_arr_poten-1
                   poten_int=> integralstore_3c_poten(i_meta:i_last)
                   i_meta=i_last+1
                end if

                do n1=1,dim1
                   eigvec1 => eigvec(i_gamma)%m(:,n1,i_s)
                   do i=1,item_arr_poten
                      V_moao(i,n1,m)=V_moao(i,n1,m)+mul*poten_int(i)*eigvec1(n)
                      V_moao(i,n1,n)=V_moao(i,n1,n)+mul*poten_int(i)*eigvec1(m)
                   end do
                enddo
             enddo
          end do

          if(mn) then
             buffer=> cpks3c(i_gamma,i_s)%V_mn
          else if(ai) then 
             buffer=> cpks3c(i_gamma,i_s)%V_ai
          end if

          if(mn) then
             eigvec2 => eigvec(i_gamma)%m
             k=0
          else if(ai) then
             eigvec2 => eigvec(i_gamma)%m
             k=dim-dim2
          end if

          do i=1,item_arr_poten
             do n1=1,dim
                do n2=1,dim2
                   buffer(i,:,n2) = buffer(i,:,n2)- &
                        V_moao(i,:,n1)*eigvec2(n1,k+n2,i_s)
                end do
             end do
          end do
          deallocate(V_moao)
       end do
       if (integrals_on_file) then 
          call poten_integral_close(th_poten)
       endif
    end do

    if (integrals_on_file) then 
       deallocate(poten_int,STAT=status)
       if(status.ne.0) call error_handler&
            ("solv_2nd_deriv_module:calc_cpks_V_vec: deallocation of poten_int failed")
    endif
       
    deallocate(dim_irrep,stat=status)
    if ( status .ne. 0) call error_handler( &
         "solv_2nd_deriv_module:calc_cpks_V_vec: deallocation of dim_irrep failed" )

  end subroutine calc_cpks_V_mo
  !********************************************************************
#if 0 /*  subroutine calc_solv_ham_mo(oo_ou) */
    !Transform the electrostatic potential matrix elements from
    !atomic orbital basis to molecular one. Parts of the final vector cpks_V_vec
    !are kept on corresponding hosts
    use symmetry_data_module, only: symmetry_data_n_spin,symmetry_data_n_irreps
    use symmetry_data_module, only: symmetry_data_n_partners,symmetry_data_dimension
    use cpksdervs_matrices, only: cpks,cpks3c
    use eigen_data_module, ONLY : eigvec
    use solv_electrostat_module, only:ham_solv_el_keep 
    !** End of interface *****************************************
    character(len=2) :: oo_ou

    real(r8_kind),pointer :: eigvec1(:),eigvec2(:)
    integer(i4_kind) :: n_irrep,i_gamma,n_spin
    logical ::  spin_polarized
    integer(i4_kind) :: n,status,i_mn,m,n_mn,k,i_meta,i_last,i_s,n1,n2,i_grad,n_partners
    integer(i4_kind),allocatable :: dim_irrep(:)
    integer(i4_kind) :: dim1,dim2,eig_dim
    real(r8_kind) :: mul
    logical :: ai,mn
    real(r8_kind),pointer :: buffer(:,:)


    ai=.false.; mn=.false.
    ai=(index(trim(oo_ou),'ai')==1)
    mn=(index(trim(oo_ou),'mn')==1)

    n_spin =symmetry_data_n_spin()
    spin_polarized = n_spin > 1


       n_irrep = symmetry_data_n_irreps()
       allocate(dim_irrep(n_irrep),stat=status)
       do n=1,n_irrep
          dim_irrep(n)  = symmetry_data_dimension(n)
       enddo

    do i_gamma=1,n_irrep
       do i_s=1,n_spin
       dim1=size(cpks(1,i_gamma,i_s)%h1ai,1)
          if(ai) then
             dim2=size(cpks(1,i_gamma,i_s)%h1ai,2)
             allocate(cpks3c(i_gamma,i_s)%solv_ham(dim1,dim2),stat=status)
          end if
          if(mn) then
             dim2=dim1
             allocate(cpks3c(i_gamma,i_s)%solv_ham(dim1,dim2),stat=status)
          end if
          cpks3c(i_gamma,i_s)%solv_ham=0.0_r8_kind
       enddo
    enddo

    do i_gamma = 1,n_irrep
       eig_dim=size(eigvec(i_gamma)%m,1)
       do m=1,dim_irrep(i_gamma)
          do n=1,m
             mul=2.0_r8_kind 
             if (n==m) mul=1.0_r8_kind

             do i_s=1,n_spin
                buffer=> cpks3c(i_gamma,i_s)%solv_ham

                dim1=size(buffer,1)
                dim2=size(buffer,2)
                do n1=1,dim1
                   eigvec1 => eigvec(i_gamma)%m(:,n1,i_s)
                   do n2=1,dim2
                      if(mn) eigvec2 => eigvec(i_gamma)%m(:,n2,i_s)
                      if(ai) eigvec2 => eigvec(i_gamma)%m(:,n2-dim2+eig_dim,i_s)
                      buffer(n1,n2) = buffer(n1,n2)- &
                           mul*ham_solv_el_keep(i_gamma)%m(m,n)*(eigvec1(m)*eigvec2(n)+eigvec1(n)*eigvec2(m))
                   enddo
                enddo
             enddo
          end do
       end do
     print*,'ham_solv_el_keep MO',sum(cpks3c(i_gamma,1)%solv_ham)
    end do
       
    deallocate(dim_irrep,stat=status)

#endif /*  end subroutine calc_solv_ham_mo */

  !********************************************************************
  subroutine shutdown_V_vec(oo_ou)
    use symmetry_data_module, only: symmetry_data_n_spin,symmetry_data_n_irreps
    use cpksdervs_matrices, only: cpks3c

    character(len=2) :: oo_ou
    integer(i4_kind) :: n_irrep,i_gamma,n_spin,i_s,status
    logical :: ai,mn
    real(r8_kind),pointer :: buffer(:,:,:)

    ai=.false.; mn=.false.
    ai=(index(trim(oo_ou),'ai')==1)
    mn=(index(trim(oo_ou),'mn')==1)

    n_irrep = symmetry_data_n_irreps()
    n_spin =symmetry_data_n_spin()
    do i_gamma=1,n_irrep
       do i_s=1,n_spin
          if(mn) then
             buffer=> cpks3c(i_gamma,i_s)%V_mn
          else if(ai) then
             buffer=> cpks3c(i_gamma,i_s)%V_ai
          end if
          deallocate(buffer,stat=status)
          if ( status .ne. 0) call error_handler( &
               "solv_2nd_deriv_module:shutdown_V_vec: deallocation of cpks3c%V  failed" )
       enddo
    enddo
  end subroutine shutdown_V_vec
  !********************************************************************

  !********************************************************************
  subroutine calc_P1(oo_ou)
    use symmetry_data_module, only: symmetry_data_n_spin,symmetry_data_n_irreps
    use symmetry_data_module, only: symmetry_data_dimension
    use cpksdervs_matrices, only: cpks
    use potential_module, only: get_bounds_poten,poten_bounds
    use comm_module, only: comm_myindex,comm_i_am_master,comm_init_send
    use comm_module, only: comm_send,comm_get_n_processors,comm_save_recv,comm_master_host
    use comm_module, only: commpack, communpack
    use msgtag_module, only: msgtag_cpks_solv1

    character(len=2) :: oo_ou

    type(poten_bounds) :: bounds
    integer(i4_kind) :: item_arr_poten,my_ind,n_irrep,i_gamma,n_spin
    integer(i4_kind) :: i_spin,i_grad,status,n_pr,info
    logical ::  spin_polarized
    integer(i4_kind) :: dim1,dim2,i

    logical :: ai,mn
    logical :: p0

    ai=.false.; mn=.false.
    p0=.false.
    ai=(index(trim(oo_ou),'ai')==1)
    mn=(index(trim(oo_ou),'mn')==1)
    p0=(index(trim(oo_ou),'p0')==1)

    my_ind=comm_myindex()

    call get_bounds_poten(bounds)
    item_arr_poten=bounds%item_arr(my_ind)
    if (item_arr_poten == 0 .and. .not.comm_i_am_master()) then 
       return
    end if
    n_spin =symmetry_data_n_spin()
    spin_polarized = n_spin > 1
    n_irrep = symmetry_data_n_irreps()

    do i_gamma=1,n_irrep
       do i_spin=1,n_spin
          dim1=size(cpks(1,i_gamma,i_spin)%Qai,1)
          if(mn.or.p0) then
             dim2=dim1
          else if(ai) then
             dim2=size(cpks(1,i_gamma,i_spin)%Qai,2)
          end if
          do i_grad=1,size(cpks,1)
             allocate(cpks(i_grad,i_gamma,i_spin)%P1(dim1,dim2),stat=status)
             if ( status .ne. 0) call error_handler( &
                  "solv_2nd_deriv_module:calc_P1: allocation of P1 failed" )
             if(comm_i_am_master()) then
                if(ai) then 
                   cpks(i_grad,i_gamma,i_spin)%P1=cpks(i_grad,i_gamma,i_spin)%HBH
                else if(mn) then
                   cpks(i_grad,i_gamma,i_spin)%P1=-cpks(i_grad,i_gamma,i_spin)%S1/2.0_r8_kind
                elseif(p0) then
                 cpks(i_grad,i_gamma,i_spin)%P1=0.0_r8_kind
                 do i=1,size(cpks(i_grad,i_gamma,i_spin)%P1,1)
                  cpks(i_grad,i_gamma,i_spin)%P1(i,i)=1.0_r8_kind
                 enddo
                end if
             endif
          end do
       end do
    end do

    if(comm_i_am_master()) then
       n_pr=comm_get_n_processors()
       do i=2,n_pr
          item_arr_poten=bounds%item_arr(i)
          if(item_arr_poten==0) cycle
          call comm_init_send(i,msgtag_cpks_solv1)
          do i_gamma=1,n_irrep
             do i_spin=1,n_spin
                do i_grad=1,size(cpks,1)
                   dim1=size(cpks(i_grad,i_gamma,i_spin)%P1,1)
                   dim2=size(cpks(i_grad,i_gamma,i_spin)%P1,2)
                   call commpack(cpks(i_grad,i_gamma,i_spin)%P1(1,1), &
                       dim1*dim2,1,info)
                end do
             end do
          end do
          call comm_send()
       end do
    else
       call comm_save_recv(comm_master_host,msgtag_cpks_solv1)
       do i_gamma=1,n_irrep
          do i_spin=1,n_spin
             do i_grad=1,size(cpks,1)
                dim1=size(cpks(i_grad,i_gamma,i_spin)%P1,1)
                dim2=size(cpks(i_grad,i_gamma,i_spin)%P1,2)
                call communpack(cpks(i_grad,i_gamma,i_spin)%P1(1,1), &
                     dim1*dim2,1,info)
             end do
          end do
       end do
    end if

  end subroutine calc_P1
  !********************************************************************

  !********************************************************************
  subroutine calc_ABai_Qai_solv(ABai_Qai_H1,h0)
    ! Calculation of X_mn=-Sum(A_inv_ij*Sum(P1_kl*V_klj)*V_mni)
    use symmetry_data_module, only: symmetry_data_n_spin,symmetry_data_n_irreps
    use symmetry_data_module, only: symmetry_data_dimension
    use cpksdervs_matrices, only: cpks,cpks3c
    use potential_module, only: get_bounds_poten,poten_bounds
    use comm_module, only: comm_myindex,comm_i_am_master,comm_init_send,comm_send
    use comm_module, only: comm_get_n_processors,comm_save_recv,comm_master_host
    use comm_module, only: commpack, communpack
    use msgtag_module, only: msgtag_cpks_solv1,msgtag_cpks_solv2,msgtag_cpks_solv3,msgtag_cpks_solv4
    use solv_cavity_module, only: to_calc_grads
    use solv_electrostat_module, only: A_matrix_inv
    implicit none
    character(len=*), intent(in)    :: ABai_Qai_H1
    character(len=2),optional,intent(in) :: h0
    ! *** end of interface ***
    
    logical :: ABai,Qai,H1,H1_1
    real(r8_kind), pointer :: buffer(:,:),buffer1(:,:,:),buffer2(:,:,:),buffer3(:)
    type(poten_bounds) :: bounds
    integer(i4_kind) :: bnd,n_points
    real(r8_kind), allocatable :: P1xV(:,:),P1xV_buf(:,:),P1xV_tot(:,:),AxP1xVxV(:,:)
    integer(i4_kind) :: item_arr_poten,my_ind,n_irrep,i_gamma,n_spin,item_arr_poten1
    integer(i4_kind) :: i_spin,i_grad,status,n_pr,info
    logical ::  spin_polarized
    integer(i4_kind) :: n_occ,n_vir,i,j,k,dim1,dim2,n_grads
    real(r8_kind) :: coeff
#if 0
    real(r8_kind):: q_solv(to_calc_grads%n_points)
#endif

    ABai=.false.; Qai=.false.; H1=.false.; H1_1=.false.
    ABai=(index(trim(ABai_Qai_H1),'ABai')==1)
    Qai=(index(trim(ABai_Qai_H1),'Qai')==1)
    H1=(index(trim(ABai_Qai_H1),'H1')==1)
    H1_1=(index(trim(ABai_Qai_H1),'H1_1')==1)

    coeff=(to_calc_grads%dielconst-1.0_r8_kind)/to_calc_grads%dielconst

    my_ind=comm_myindex()

    call get_bounds_poten(bounds)
    item_arr_poten=bounds%item_arr(my_ind)
    if (item_arr_poten == 0 .and. .not.comm_i_am_master()) then 
       return
    end if
    n_points=to_calc_grads%n_points
    n_spin =symmetry_data_n_spin()
    spin_polarized = n_spin > 1
    n_irrep = symmetry_data_n_irreps()
    n_grads=size(cpks,1)

    !Store for Global Sum(P1_kl*V_klj) 
    allocate(P1xV_tot(n_points,n_grads),stat=status)
    if ( status .ne. 0) call error_handler( &
         "solv_2nd_deriv_module:calc_ABai_Qai_solv: allocation of P1xV_tot failed" )
    P1xV_tot=0.0_r8_kind

    if (item_arr_poten /= 0) then
       !Allocation of local array  P1xV  on each node
       allocate(P1xV(item_arr_poten,n_grads),stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_2nd_deriv_module:calc_ABai_Qai_solv: allocation of P1xV failed" )
       P1xV=0.0_r8_kind
    end if

    !calculation of P1xV on each node
    do i_gamma=1,n_irrep
       do i_spin=1,n_spin
          if(ABai.or.H1_1) then 
             buffer1=>cpks3c(i_gamma,i_spin)%V_ai
          else if(Qai.or.H1) then 
             buffer1=>cpks3c(i_gamma,i_spin)%V_mn
          end if

          do i_grad=1,n_grads
             dim1=size(cpks(i_grad,i_gamma,i_spin)%P1,1)
             dim2=size(cpks(i_grad,i_gamma,i_spin)%P1,2)
             do k=1,item_arr_poten
                do i=1,dim1
                   do j=1,dim2
                      P1xV(k,i_grad)=P1xV(k,i_grad)+cpks(i_grad,i_gamma,i_spin)%P1(i,j)*buffer1(k,i,j)
                   end do
                end do
             end do

             deallocate(cpks(i_grad,i_gamma,i_spin)%P1,stat=status)
             if ( status .ne. 0) call error_handler( &
                  "solv_2nd_deriv_module:calc_ABai_Qai_solv: deallocation of P1 failed" )
          end do
       end do
    end do

    !Sending local P1xV arrays to  master
    if(comm_i_am_master()) then
       !Gathering local P1xV arrays on master mode
       if (item_arr_poten /= 0) P1xV_tot(1:item_arr_poten,1:n_grads)=P1xV
       n_pr=comm_get_n_processors()
       bnd=item_arr_poten
       do i=2,n_pr
          item_arr_poten1=bounds%item_arr(i)
          if(item_arr_poten==0) cycle
          allocate(P1xV_buf(item_arr_poten1,n_grads))
          call comm_save_recv(i,msgtag_cpks_solv2)
          call communpack(P1xV_buf(1,1),item_arr_poten1*n_grads,1,info)
          P1xV_tot(bnd+1:bnd+item_arr_poten1,1:n_grads)=P1xV_buf
          deallocate(P1xV_buf)
          bnd=bnd+item_arr_poten1
       end do

       !Sending global P1xV_tot array to each node participated
       do i=2,n_pr
          item_arr_poten1=bounds%item_arr(i)
          if(item_arr_poten1==0) cycle
          call comm_init_send(i,msgtag_cpks_solv3)
          call commpack(P1xV_tot(1,1),N_points*n_grads,1,info)
          call comm_send()
       end do
    else
       call comm_init_send(comm_master_host,msgtag_cpks_solv2)
       call commpack(P1xV(1,1),item_arr_poten*n_grads,1,info)
       call comm_send()
       
       call comm_save_recv(comm_master_host,msgtag_cpks_solv3)
       call communpack(P1xV_tot(1,1),N_points*n_grads,1,info)
    end if
    if (item_arr_poten /= 0) then
       deallocate(P1xV,stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_2nd_deriv_module:calc_ABai_Qai_solv: deallocation of P1xV failed" )
    end if

    !Distribution of the matrix A_inv between slaves
    if(comm_i_am_master()) then
       n_pr=comm_get_n_processors()
       bnd=item_arr_poten
       do i=2,n_pr
          item_arr_poten1=bounds%item_arr(i)
          if(item_arr_poten==0) cycle
          call comm_init_send(i,msgtag_cpks_solv1)
          call commpack(A_matrix_inv(1,bnd+1),N_points*item_arr_poten1,1,info)
          call comm_send()
          bnd=bnd+item_arr_poten1
       end do
    else
       call comm_save_recv(comm_master_host,msgtag_cpks_solv1)
       allocate(A_matrix_inv(N_points,item_arr_poten),stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_2nd_deriv_module:calc_ABai_Qai_solv: allocation of A_matrix_inv failed" )
       call communpack(A_matrix_inv(1,1),N_points*item_arr_poten,1,info)
    end if

    !calculation of ABi solv
    !Cycles i_gamma, i_spin and i_grad run simulteniously on all nodes participated
#if 1
    g1: do i_gamma=1,n_irrep
       do i_spin=1,n_spin
          n_occ=size(cpks(1,i_gamma,i_spin)%Qai,1)
          n_vir=size(cpks(1,i_gamma,i_spin)%Qai,2)
          if(ABai.or.Qai) then 
             buffer2=>cpks3c(i_gamma,i_spin)%V_ai
          else if(H1.or.H1_1) then 
             buffer2=>cpks3c(i_gamma,i_spin)%V_mn
             n_vir=n_occ
          end if
          do i_grad=1,n_grads

             if(comm_i_am_master()) then
                if(ABai) then 
                   buffer => cpks(i_grad,i_gamma,i_spin)%ABi
                else if(Qai) then
                   buffer => cpks(i_grad,i_gamma,i_spin)%Qai
                else if(H1.or.H1_1) then
                   buffer => cpks(i_grad,i_gamma,i_spin)%H1
                end if
             end if



             !Calculation of local solvent contribution of each node to ABi matrix
             allocate(AxP1xVxV(n_occ,n_vir),stat=status)
             if ( status .ne. 0) call error_handler( &
                  "solv_2nd_deriv_module:calc_ABai_Qai_solv: allocation of AxP1xVxV failed" )
             AxP1xVxV=0.0_r8_kind
             if (item_arr_poten /= 0) then
                allocate(buffer3(item_arr_poten))
                buffer3=0.0_r8_kind
                do i=1,n_points
                   buffer3(:)=buffer3(:)+coeff*A_matrix_inv(i,1:item_arr_poten)*P1xV_tot(i,i_grad)
                enddo
                do j=1,item_arr_poten
                   AxP1xVxV=AxP1xVxV-buffer3(j)*buffer2(j,:,:)
                enddo
                deallocate(buffer3)
             end if

#if 0
            if(i_grad.eq.1.and.present(h0)) then
             q_solv=-matmul(P1xV_tot,A_matrix_inv)
             q_solv=q_solv*coeff
             print*,'q_solv A_matrix_inv', &
                    q_solv(1),to_calc_grads%Q(1),sum(q_solv),sum(to_calc_grads%Q), &
                    sum(A_matrix_inv)
            q_solv=0.0_r8_kind
            call fix_solv_q(q_solv)
            do i=1,n_points
              AxP1xVxV=AxP1xVxV+q_solv(i)*buffer2(i,:,:)
            enddo
            endif
#endif

             if(comm_i_am_master()) then
                !Summing up all local solvent contribution of ABi on master node
                buffer=buffer+AxP1xVxV
                n_pr=comm_get_n_processors()
                do i=2,n_pr
                   item_arr_poten1=bounds%item_arr(i)
                   if(item_arr_poten==0) cycle
                   call comm_save_recv(i,msgtag_cpks_solv4)
                   call communpack(AxP1xVxV(1,1),n_occ*n_vir,1,info)
                   buffer=buffer+AxP1xVxV
                end do
             else
                !Sending local solvent contribution to ABi from slaves to master
                call comm_init_send(comm_master_host,msgtag_cpks_solv4)
                call commpack(AxP1xVxV(1,1),n_occ*n_vir,1,info)
                call comm_send()
             end if

             deallocate(AxP1xVxV,stat=status)
             if ( status .ne. 0) call error_handler( &
                  "solv_2nd_deriv_module:calc_ABai_Qai_solv: deallocation of AxP1xVxV failed" )
          end do
       end do
    enddo g1
#endif

    deallocate(P1xV_tot,stat=status)
    if ( status .ne. 0) call error_handler( &
         "solv_2nd_deriv_module:calc_ABai_solv: deallocation of P1xV_tot failed" )

    if(.not.comm_i_am_master()) then
       deallocate(A_matrix_inv,stat=status)
       if ( status .ne. 0) call error_handler( &
            "solv_2nd_deriv_module:calc_ABai_Qai_solv: deallocation of A_matrix_inv failed" )
    end if
  end subroutine calc_ABai_Qai_solv
  !********************************************************************

  !********************************************************************
  subroutine cpks_solv_main()

    call calc_P1('ai')
    call calc_ABai_Qai_solv('ABai')

  end subroutine cpks_solv_main
  !********************************************************************

  !********************************************************************
  subroutine cpks_solv_Qai()

    call calc_cpks_V_mo('ai')
    call calc_cpks_V_mo('mn')
    call calc_P1('mn')
    call calc_ABai_Qai_solv('Qai')
    call shutdown_V_vec('mn')


  end subroutine cpks_solv_Qai
  !********************************************************************
#if 0 
  subroutine cpks_solv_ham()

    use comm_module, only : comm_myindex

    call calc_cpks_V_mo('ai')
    call calc_cpks_V_mo('mn')
    call calc_P1('p0')
    call calc_ABai_Qai_solv('Qai','h0')
    call shutdown_V_vec('mn')
    call shutdown_V_vec('ai')
  end subroutine cpks_solv_ham
#endif 

  !********************************************************************
  subroutine cpks_solv_H1()

    call calc_cpks_V_mo('mn')
    call calc_P1('mn')
    call calc_ABai_Qai_solv('H1')

  end subroutine cpks_solv_H1
  !********************************************************************

  !********************************************************************
  subroutine cpks_solv_H1_1()

    use integralstore_module, only: integralstore_deallocate_pcm
    use options_module, only: options_integrals_on_file

    call calc_P1('ai')
    call calc_ABai_Qai_solv('H1_1')
    call shutdown_V_vec('mn')
    call shutdown_V_vec('ai')
    if (.not.options_integrals_on_file()) then 
       call integralstore_deallocate_pcm()
    endif

  end subroutine cpks_solv_H1_1
  !********************************************************************

  !********************************************************************
end module solv_2nd_deriv_module


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
  subroutine calc_sym_coef_eperef(ua3,lmax_ch)
#include <def.h>
  use epe_module
  use type_module
  use unique_atom_module
  use epecom_module, only:  dealloc_epe_ref
  implicit none
   integer(kind=i4_kind), intent(in)::ua3,lmax_ch

   integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
   real(kind=r8_kind),pointer      :: coef(:)
   integer(kind=i4_kind):: n_indep_fcts,i_l,i_ind,i_cont,lm,i,n_contributing_fcts
    ! Purpose: piece of the code
    !          calculate symmetry coefficients for the l-type 3 center
    !          fitintegrals
    ang_momentum_symadapt: do i_l=1,lmax_ch
       n_indep_fcts =  &
          unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       do i_ind=1,n_indep_fcts
          n_contributing_fcts = &
               unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%N_fcts
          magn => unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          coef =>  unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          eq_atom => unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%I_equal_atom

          do i_cont=1,n_contributing_fcts
             lm = i_l**2 + magn(i_cont)
             sym_coef1(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coef1(:,eq_atom(i_cont),i_ind,i_l) + &
                  yl_arr(:,lm,eq_atom(i_cont))*&
                  coef(i_cont)
             do i = 1,3
                sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) = &
                     sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) + &
                     yl_arr_grad(:,eq_atom(i_cont),lm,i)*coef(i_cont)
             enddo
          enddo   ! i_cont=1,n_contributing_fcts

     if(dealloc_epe_ref) then
      deallocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m&
       ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c&
       ,unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%I_equal_atom, &
        stat=epemalloc_stat(36))
       ASSERT(epemalloc_stat(36).eq.0)
      endif   ! dealloc_epe_ref
       enddo
       if(dealloc_epe_ref) then
        deallocate(unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%symadapt, &
                                                         stat=epemalloc_stat(35)) ! %symadapt
        ASSERT(epemalloc_stat(35).eq.0)
       endif
    enddo ang_momentum_symadapt
  end subroutine calc_sym_coef_eperef

  subroutine calc_sym_coef(ua3,lmax_ch)
  use epe_module
  use type_module
  use unique_atom_module
  implicit none
   integer(kind=i4_kind), intent(in)::ua3,lmax_ch

   integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
   real(kind=r8_kind),pointer      :: coef(:)
   integer(kind=i4_kind):: n_indep_fcts,i_l,i_ind,i_cont,lm,i,n_contributing_fcts
    !          calculates symmetry coefficients for the l-type 3 center
    !          fitintegrals

    ang_momentum_symadapt: do i_l=1,lmax_ch
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
 
       do i_ind=1,n_indep_fcts
          n_contributing_fcts = &
               unique_atoms(ua3)%symadapt_partner(1,i_l)%&
                                  symadapt(i_ind,1)%N_fcts
          magn => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%m
          coef =>  unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%c
          eq_atom => unique_atoms(ua3)%symadapt_partner(1,i_l)%&
               symadapt(i_ind,1)%I_equal_atom
          do i_cont=1,n_contributing_fcts
             lm = i_l**2 + magn(i_cont)
             sym_coef1(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coef1(:,eq_atom(i_cont),i_ind,i_l) + &
                  yl_arr(:,lm,eq_atom(i_cont))*&
                  coef(i_cont)
             do i = 1,3
                sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) = &
                     sym_coef2(:,eq_atom(i_cont),i_ind,i_l,i) + &
                     yl_arr_grad(:,eq_atom(i_cont),lm,i)*coef(i_cont)
             enddo
          enddo
       enddo
    enddo ang_momentum_symadapt
  end subroutine calc_sym_coef
  subroutine epe_l_coulomb(ua3,epe_reference,grad_dim,lmax_ch)
    ! Purpose: claculate gradients of l-type fit integrals
  use epe_module
  use type_module
  use gamma_module
  use unique_atom_module
  implicit none
    integer(kind=i4_kind),intent(in):: grad_dim,ua3,lmax_ch
    integer(kind=i4_kind):: n_indep_fcts,ncexps,n_equal_c,i_ind,i_ea3, &
                            i_grad,i_l,k,i
	 logical::epe_reference
     real(kind=r8_kind):: xc(3)

   ang_momentum: do i_l=1,lmax_ch
   if(epe_reference) then
       n_indep_fcts =  &
            unique_atoms_eperef(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms_eperef(ua3)%l_ch(i_l)%N_exponents
       cexps => unique_atoms_eperef(ua3)%l_ch(i_l)%exponents(:)
       n_equal_c=unique_atoms_eperef(ua3)%n_equal_atoms
    else
       n_indep_fcts =  &
            unique_atoms(ua3)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms(ua3)%l_ch(i_l)%N_exponents
       cexps => unique_atoms(ua3)%l_ch(i_l)%exponents(:)
       n_equal_c=unique_atoms(ua3)%n_equal_atoms
    endif ! get_epe_reference
       do i_ind=1,n_indep_fcts
          equal_3_l: do i_ea3=1,n_equal_c
        do_rotation=.false.
        if(epe_reference) then
             xc =  unique_atoms_eperef(ua3)%position(:,i_ea3)
        else
             xc =  unique_atoms(ua3)%position(:,i_ea3)
        endif
             do k = 1,ncexps
                gamma_arg2 = cexps(k)* &
                     (  (gamma_arg(:,1) - xc(1))**2 + &
                        (gamma_arg(:,2) - xc(2))**2 + &
                        (gamma_arg(:,3) - xc(3))**2) 
                gamma_help(:,1:lmax_ch+2) = gamma(lmax_ch+2,gamma_arg2)
                help_vec1=two*pi/cexps(k)*(two*cexps(k))**i_l
                help_vec2=gamma_help(:,i_l+1)
                help_vec3=two*cexps(k)* &
                     gamma_help(:,i_l+2)*sym_coef1(:,i_ea3,i_ind,i_l)
                do i = 1,3
                   help_vec4=help_vec3*(gamma_arg(:,i)-xc(i))
                   help_vec7=help_vec2*sym_coef2(:,i_ea3,i_ind,i_l,i)
                   help_arr(:,i) =  help_vec1* ( -help_vec7 +  help_vec4)
                enddo
         coul_int_c(0)%l(i_l)%m(:,k,i_ind,1,1)=&
         coul_int_c(0)%l(i_l)%m(:,k,i_ind,1,1) &
                  +help_vec1*sym_coef1(:,i_ea3,i_ind,i_l)*gamma_help(:,i_l+1)
                if(do_rotation) then
                   ! make gradient totalsymmetric before adding
                   do i_grad=1,grad_dim ! only if moving_c
                      coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)=&
                           coul_int_c(i_grad)%l(i_l)%m(:,k,i_ind,1,1)+&
                           rotmat(i_grad,1)*help_arr(:,1)+&
                           rotmat(i_grad,2)*help_arr(:,2)+&
                           rotmat(i_grad,3)*help_arr(:,3)
                   enddo
                else
                   coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(1)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,1)
                   coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(2)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,2)
                   coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)=&
                        coul_int_c(3)%l(i_l)%m(:,k,i_ind,1,1)+&
                        help_arr(:,3)
                end if
             enddo! loop over k
          enddo equal_3_l
       enddo
    enddo ang_momentum
  end subroutine epe_l_coulomb
  subroutine epe_s_coulomb(ncexps,grad_dim)
   use epe_module
   use type_module
   use gamma_module
   implicit none
    ! Purpose: calculate gradients of s type coulomb fitintegrals
    ! loop over exponents of third center ---** s-type Fitfct. **---
    real(kind=r8_kind):: xc(3)
    integer(kind=i4_kind),intent(in):: ncexps,grad_dim
    integer(kind=i4_kind):: k,i,i_grad
!    real(kind=r8_kind), dimension(num) :: help_vec1

    do k=1,ncexps
       gamma_arg2 = cexps(k)* &
           ((gamma_arg(:,1) - xc(1))**2 + &
            (gamma_arg(:,2) - xc(2))**2 + &
            (gamma_arg(:,3) - xc(3))**2)

       gamma_help(:,1:2) = gamma(2,gamma_arg2)
       help_vec1=four*pi*gamma_help(:,2)
       coul_int_c(0)%l(0)%m(:,k,1,1,1)=coul_int_c(0)%l(0)%m(:,k,1,1,1)&
                                       +two*pi/cexps(k)*gamma_help(:,1)
          do i = 1,3
             help_arr(:,i) =  help_vec1*(gamma_arg(:,i)-xc(i))
          enddo
       if(do_rotation) then
          ! make gradient totalsymmetric before adding
          do i_grad=1,grad_dim ! only if moving_c
             coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)=&
                  coul_int_c(i_grad)%l(0)%m(:,k,1,1,1)+&
                  rotmat(i_grad,1)*help_arr(:,1)+&
                  rotmat(i_grad,2)*help_arr(:,2)+&
                  rotmat(i_grad,3)*help_arr(:,3)
          enddo
       else
          coul_int_c(1)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(1)%l(0)%m(:,k,1,1,1)+help_arr(:,1)
          coul_int_c(2)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(2)%l(0)%m(:,k,1,1,1)+help_arr(:,2)
          coul_int_c(3)%l(0)%m(:,k,1,1,1)=&
               coul_int_c(3)%l(0)%m(:,k,1,1,1)+help_arr(:,3)
       endif
    enddo! s-exponents of third center 
  end subroutine epe_s_coulomb

  subroutine epe_r2_coulomb(ncexps,grad_dim) ! 1
   use epe_module
   use type_module
   use gamma_module
!   use unique_atom_module
   implicit none

    real(kind=r8_kind):: xc(3)
    integer(kind=i4_kind), intent(in):: ncexps,grad_dim
    integer(kind=i4_kind):: k,i,i_grad
    ! Purpose: calculate gradients of r2 coulomb integrals
    !---** r2-type Fitfct. **---

    do k=1,ncexps

       gamma_arg2 = cexps(k)* &
            (  (gamma_arg(:,1) - xc(1))**2 + &
               (gamma_arg(:,2) - xc(2))**2 + &
               (gamma_arg(:,3) - xc(3))**2)
       gamma_help(:,1:3) = gamma(3,gamma_arg2)

       help_vec1=two*pi/(cexps(k)**2)
       coul_int_c(0)%l(-1)%m(:,k,1,1,1)=coul_int_c(0)%l(-1)%m(:,k,1,1,1)&
                +help_vec1*( gamma_help(:,1)+gamma_arg2*gamma_help(:,2) )
       help_vec3= two*((-gamma_arg2*gamma_help(:,3))*cexps(k))
       help_vec4=(gamma_help(:,1)+gamma_arg2*gamma_help(:,2))
          help_vec1=help_vec1*help_vec3
          do i = 1,3
             help_arr(:,i) = -help_vec1*(gamma_arg(:,i)-xc(i))
          enddo

       if(do_rotation) then
          ! make gradient totalsymmetric before adding
          do i_grad=1,grad_dim ! only if moving_c
             coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)=&
                  coul_int_c(i_grad)%l(-1)%m(:,k,1,1,1)+&
                  rotmat(i_grad,1)*help_arr(:,1)+&
                  rotmat(i_grad,2)*help_arr(:,2)+&
                  rotmat(i_grad,3)*help_arr(:,3)
          enddo
       else
          coul_int_c(1)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(1)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,1)
          coul_int_c(2)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(2)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,2)
          coul_int_c(3)%l(-1)%m(:,k,1,1,1)=&
               coul_int_c(3)%l(-1)%m(:,k,1,1,1)+&
               help_arr(:,3)
       end if
    enddo! r2-exponents, third center
  end subroutine epe_r2_coulomb

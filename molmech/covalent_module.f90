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
module covalent_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use slab_module
  use species_module
  use external_field_module
  use tasks_main_options_module
  use energy_and_forces_module
  use potentials_module
  use n_body_lists_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ------------------
  public two_body_E_and_F, three_body_E_and_F, four_body_E_and_F,  &
       many_body_E_and_F, core_shell_E_and_F
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine core_shell_E_and_F()

    integer(kind=i4_kind) :: i,j,k,length,ia1,ia2,id,k1,k2,l,m,m1,ns
    real(kind=r8_kind) :: r(3),dr,dr2,dr3,E_buf
    real(kind=r8_kind) :: f1(3),f2(3),dE_dr,V,VV,HH
    real(kind=r8_kind) :: d2E_dr2,ff(6,6),fstore
    real(kind=r8_kind) :: p(n_parameter)

    ns=3*n_species

    i_n_2b: do i=1,n_cs
       k=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+n_2b_n+i
       length=size(core_shell(i)%list,2)

       j_length: do j=1,length
          ia1=core_shell(i)%list(1,j)
          ia2=core_shell(i)%list(2,j)
          r=atoms_cart(ia1)%r-atoms_cart(ia2)%r
          if(lattice_calc) then 
             call image(r)
          else if(slab_calc) then
             call image_slab(r)
          end if
          dr=sqrt(dot_product(r,r))
          if(dr < small2) dr=small2
!!$          if(dr < small0*ten) cycle
          dr2=dr*dr
          dr3=dr2*dr

          id=poten(k)%id
          p=poten(k)%param
          select case (id)
          case (15)  !core_shell
             E_buf=(p(1)/two)*dr2
!!$             if(abs(E_buf) < small1) E_buf=zero
             E(15)=E(15)+E_buf
             if(calc_gradients) dE_dr=p(1)*dr
             if(calc_hessian) d2E_dr2=p(1)
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             f1=r/dr
             f2=-f1
             Grad(:,ia1)=Grad(:,ia1)+dE_dr*f1
             Grad(:,ia2)=Grad(:,ia2)+dE_dr*f2
             if(calc_strain) then
                V=dE_dr/dr
                if(lattice_calc) then
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                   Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                   Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                   Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
                else if(slab_calc) then
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                end if
             end if
          end if

          if(calc_hessian) then
             ff(1,1)=one/dr-r(1)*r(1)/dr3
             ff(1,2)=-r(1)*r(2)/dr3
             ff(1,3)=-r(1)*r(3)/dr3
             ff(1,4)=-ff(1,1)
             ff(1,5)=-ff(1,2)
             ff(1,6)=-ff(1,3)
             ff(2,1)=-r(2)*r(1)/dr3
             ff(2,2)=one/dr-r(2)*r(2)/dr3
             ff(2,3)=-r(2)*r(3)/dr3
             ff(2,4)=-ff(2,1)
             ff(2,5)=-ff(2,2)
             ff(2,6)=-ff(2,3)
             ff(3,1)=-r(3)*r(1)/dr3
             ff(3,2)=-r(3)*r(2)/dr3
             ff(3,3)=one/dr-r(3)*r(3)/dr3
             ff(3,4)=-ff(3,1)
             ff(3,5)=-ff(3,2)
             ff(3,6)=-ff(3,3)
             ff(4,1:3)=ff(1:3,4)
             ff(5,1:3)=ff(1:3,5)
             ff(6,1:3)=ff(1:3,6)
             ff(4,4:6)=ff(1,1:3)
             ff(5,4:6)=ff(2,1:3)
             ff(6,4:6)=ff(3,1:3)
             
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             do l=1,3
                do m=1,6
                   if(m <= 3) then
                      fstore=f1(m)
                      m1=k1+m
                   elseif(m > 3) then 
                      m1=k2+(m-3)
                      fstore=f2(m-3)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
                end do
             end do
             
             if(calc_strain) then
                VV=(d2E_dr2-V)/(dr*dr)

                !d2E_dek_del
                if(lattice_calc) then
                   H(ns+1,ns+1)=H(ns+1,ns+1)+V*two*r(1)*r(1)+VV*r(1)*r(1)*r(1)*r(1)
                   HH=VV*r(1)*r(1)*r(2)*r(2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=VV*r(1)*r(1)*r(3)*r(3)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH
                   HH=VV*r(1)*r(1)*r(2)*r(3)
                   H(ns+1,ns+4)=H(ns+1,ns+4)+HH; H(ns+4,ns+1)=H(ns+4,ns+1)+HH
                   HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(1)*r(1)*r(1)*r(3)
                   H(ns+1,ns+5)=H(ns+1,ns+5)+HH; H(ns+5,ns+1)=H(ns+5,ns+1)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                   H(ns+1,ns+6)=H(ns+1,ns+6)+HH; H(ns+6,ns+1)=H(ns+6,ns+1)+HH
                
                   H(ns+2,ns+2)=H(ns+2,ns+2)+V*two*r(2)*r(2)+VV*r(2)*r(2)*r(2)*r(2)
                   HH=VV*r(2)*r(2)*r(3)*r(3)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH
                   HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(2)*r(2)*r(2)*r(3)
                   H(ns+2,ns+4)=H(ns+2,ns+4)+HH; H(ns+4,ns+2)=H(ns+4,ns+2)+HH
                   HH=VV*r(2)*r(2)*r(1)*r(3)
                   H(ns+2,ns+5)=H(ns+2,ns+5)+HH; H(ns+5,ns+2)=H(ns+5,ns+2)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                   H(ns+2,ns+6)=H(ns+2,ns+6)+HH; H(ns+6,ns+2)=H(ns+6,ns+2)+HH

                   H(ns+3,ns+3)=H(ns+3,ns+3)+V*two*r(3)*r(3)+VV*r(3)*r(3)*r(3)*r(3)
                   HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(3)*r(3)*r(3)*r(2)
                   H(ns+3,ns+4)=H(ns+3,ns+4)+HH; H(ns+4,ns+3)=H(ns+4,ns+3)+HH
                   HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(3)*r(3)*r(3)*r(1)
                   H(ns+3,ns+5)=H(ns+3,ns+5)+HH; H(ns+5,ns+3)=H(ns+5,ns+3)+HH
                   HH=VV*r(3)*r(3)*r(1)*r(2)
                   H(ns+3,ns+6)=H(ns+3,ns+6)+HH; H(ns+6,ns+3)=H(ns+6,ns+3)+HH
                
                   H(ns+4,ns+4)=H(ns+4,ns+4)+ &
                        V*quarter*(r(2)*r(2)+r(3)*r(3))+V*quarter*(r(2)*r(2)+r(3)*r(3))+ &
                        VV*r(2)*r(3)*r(2)*r(3)
                   HH=V*quarter*r(1)*r(2)+V*quarter*r(1)*r(2)+VV*r(2)*r(3)*r(3)*r(1)
                   H(ns+4,ns+5)=H(ns+4,ns+5)+HH; H(ns+5,ns+4)=H(ns+5,ns+4)+HH
                   HH=V*quarter*r(1)*r(3)+V*quarter*r(1)*r(3)+VV*r(2)*r(3)*r(2)*r(1)
                   H(ns+4,ns+6)=H(ns+4,ns+6)+HH; H(ns+6,ns+4)=H(ns+6,ns+4)+HH
                
                   H(ns+5,ns+5)=H(ns+5,ns+5)+&
                        V*quarter*(r(1)*r(1)+r(3)*r(3))+V*quarter*(r(1)*r(1)+r(3)*r(3))+ &
                        VV*r(1)*r(3)*r(1)*r(3)
                   HH=V*quarter*r(2)*r(3)+V*quarter*r(2)*r(3)+VV*r(1)*r(2)*r(1)*r(3)
                   H(ns+5,ns+6)=H(ns+5,ns+6)+HH; H(ns+6,ns+5)=H(ns+6,ns+5)+HH

                   H(ns+6,ns+6)=H(ns+6,ns+6)+&
                        V*quarter*(r(1)*r(1)+r(2)*r(2))+V*quarter*(r(1)*r(1)+r(2)*r(2))+ &
                        VV*r(1)*r(2)*r(1)*r(2)
                else if(slab_calc) then
                   H(ns+1,ns+1)=H(ns+1,ns+1)+V*two*r(1)*r(1)+VV*r(1)*r(1)*r(1)*r(1)
                   HH=VV*r(1)*r(1)*r(2)*r(2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH

                   H(ns+2,ns+2)=H(ns+2,ns+2)+V*two*r(2)*r(2)+VV*r(2)*r(2)*r(2)*r(2)
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH

                   H(ns+3,ns+3)=H(ns+3,ns+3)+&
                        V*quarter*(r(1)*r(1)+r(2)*r(2))+V*quarter*(r(1)*r(1)+r(2)*r(2))+ &
                        VV*r(1)*r(2)*r(1)*r(2)
                end if

                !d2E_dri_dek
                k1=3*(ia1-1); k2=3*(ia2-1)

                if(lattice_calc) then
                   HH=V*two*r(1)+VV*r(1)*r(1)*r(1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(2)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                   HH=VV*r(1)*r(3)*r(3)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+1,ns+4)=H(k1+1,ns+4)+HH; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                   H(k2+1,ns+4)=H(k2+1,ns+4)-HH; H(ns+4,k2+1)=H(ns+4,k2+1)-HH
                   HH=V*r(3)+VV*r(1)*r(1)*r(3)
                   H(k1+1,ns+5)=H(k1+1,ns+5)+HH; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                   H(k2+1,ns+5)=H(k2+1,ns+5)-HH; H(ns+5,k2+1)=H(ns+5,k2+1)-HH
                   HH=V*r(2)+VV*r(1)*r(1)*r(2)
                   H(k1+1,ns+6)=H(k1+1,ns+6)+HH; H(ns+6,k1+1)=H(ns+6,k1+1)+HH
                   H(k2+1,ns+6)=H(k2+1,ns+6)-HH; H(ns+6,k2+1)=H(ns+6,k2+1)-HH

                   HH=VV*r(2)*r(1)*r(1)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                   HH=V*two*r(2)+VV*r(2)*r(2)*r(2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                   HH=VV*r(2)*r(3)*r(3)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH
                   HH=V*r(3)+VV*r(2)*r(2)*r(3)
                   H(k1+2,ns+4)=H(k1+2,ns+4)+HH; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                   H(k2+2,ns+4)=H(k2+2,ns+4)-HH; H(ns+4,k2+2)=H(ns+4,k2+2)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+2,ns+5)=H(k1+2,ns+5)+HH; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                   H(k2+2,ns+5)=H(k2+2,ns+5)-HH; H(ns+5,k2+2)=H(ns+5,k2+2)-HH
                   HH=V*r(1)+VV*r(2)*r(1)*r(2)
                   H(k1+2,ns+6)=H(k1+2,ns+6)+HH; H(ns+6,k1+2)=H(ns+6,k1+2)+HH
                   H(k2+2,ns+6)=H(k2+2,ns+6)-HH; H(ns+6,k2+2)=H(ns+6,k2+2)-HH
                
                   HH=VV*r(3)*r(1)*r(1)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                   HH=VV*r(3)*r(2)*r(2)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                   HH=V*two*r(3)+VV*r(3)*r(3)*r(3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                   HH=V*r(2)+VV*r(3)*r(2)*r(3)
                   H(k1+3,ns+4)=H(k1+3,ns+4)+HH; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                   H(k2+3,ns+4)=H(k2+3,ns+4)-HH; H(ns+4,k2+3)=H(ns+4,k2+3)-HH
                   HH=V*r(1)+VV*r(3)*r(1)*r(3)
                   H(k1+3,ns+5)=H(k1+3,ns+5)+HH; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                   H(k2+3,ns+5)=H(k2+3,ns+5)-HH; H(ns+5,k2+3)=H(ns+5,k2+3)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+3,ns+6)=H(k1+3,ns+6)+HH; H(ns+6,k1+3)=H(ns+6,k1+3)+HH
                   H(k2+3,ns+6)=H(k2+3,ns+6)-HH; H(ns+6,k2+3)=H(ns+6,k2+3)-HH
                else if(slab_calc) then
                   HH=V*two*r(1)+VV*r(1)*r(1)*r(1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(2)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                   HH=V*r(2)+VV*r(1)*r(1)*r(2)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH

                   HH=VV*r(2)*r(1)*r(1)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                   HH=V*two*r(2)+VV*r(2)*r(2)*r(2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                   HH=V*r(1)+VV*r(2)*r(1)*r(2)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH

                   HH=VV*r(3)*r(1)*r(1)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                   HH=VV*r(3)*r(2)*r(2)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                end if
             end if
          end if
       end do j_length
    end do i_n_2b
    
!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine core_shell_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine two_body_E_and_F()

    integer(kind=i4_kind) :: i,j,length,ia1,ia2,id,k1,k2,l,m,m1,ns,l1,l2
    real(kind=r8_kind) :: r(3),dr,dr1,E_buf,dr3,dr6,dr7,dr8,r1(3),r2(3)
    real(kind=r8_kind) :: f1(3),f2(3),dE_dr,V,VV
    real(kind=r8_kind) :: d2E_dr2,ff(6,6),fstore,HH
    real(kind=r8_kind) :: p(n_parameter)

    ns=3*n_species

    i_n_2b: do i=1,n_2b
       length=size(two_body(i)%list,2)

       j_length: do j=1,length
          ia1=two_body(i)%list(1,j)
          ia2=two_body(i)%list(2,j)
          if(ia1 <= n_species) then
             r1=atoms_cart(ia1)%r
          else
             r1=pc_cart(ia1-n_species)%r
          end if
          if(ia2 <= n_species) then
             r2=atoms_cart(ia2)%r
          else
             r2=pc_cart(ia2-n_species)%r
          end if
          r=r1-r2
          if(lattice_calc) then
             call image(r)
          else if(slab_calc) then
             call image_slab(r)
          end if
          dr=sqrt(dot_product(r,r))

          id=poten(i)%id
          p=poten(i)%param
          select case (id)
          case (1)  !harm_str
             dr1=dr-p(2)
             E_buf=(p(1)/two)*dr1**2  
             E(1)=E(1)+E_buf
             if(calc_gradients) dE_dr=p(1)*dr1
             if(calc_hessian) d2E_dr2=p(1)
          case (2)  !morze
             dr1=dr-p(3)
             E_buf=p(1)*((one-exp(-p(2)*dr1))**2-one)
             E(2)=E(2)+E_buf
             if(calc_gradients) dE_dr=two*p(1)*p(2)*(one-exp(-p(2)*dr1))*exp(-p(2)*dr1)
             if(calc_hessian) d2E_dr2=two*p(1)*p(2)**2*(two*exp(-p(2)*dr1)-one)*exp(-p(2)*dr1)
          case (3)  !quart_str
             dr1=dr-p(4)
             E_buf=(p(1)/two)*dr1**2+(p(2)/three)*dr1**3+(p(3)/four)*dr1**4
             E(3)=E(3)+E_buf
             if(calc_gradients) dE_dr=p(1)*dr1+p(2)*dr1**2+p(3)*dr1**3
             if(calc_hessian) d2E_dr2=p(1)+two*p(2)*dr1+three*p(3)*dr1**2
          case (16) !bck_short
             dr3=dr*dr*dr
             dr6=dr3*dr3
             dr7=dr6*dr
             dr8=dr7*dr
             E_buf=p(1)*exp(-dr/p(2))-p(3)/dr6
             E(16)=E(16)+E_buf
             if(calc_gradients) dE_dr=-(p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7
             if(calc_hessian) d2E_dr2=(p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8
          case (n_poten+1 :) !user_def
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             f1=r/dr
             f2=-f1
             if(ia1 <= n_species) Grad(:,ia1)=Grad(:,ia1)+dE_dr*f1
             if(ia2 <= n_species) Grad(:,ia2)=Grad(:,ia2)+dE_dr*f2
             if(calc_strain) then
                V=dE_dr/dr
                if(lattice_calc) then
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                   Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                   Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                   Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
                else if(slab_calc) then
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                end if
             end if
          end if

          if(calc_hessian) then
             ff(1,1)=one/dr-r(1)*r(1)/dr**3
             ff(1,2)=-r(1)*r(2)/dr**3
             ff(1,3)=-r(1)*r(3)/dr**3
             ff(1,4)=-ff(1,1)
             ff(1,5)=-ff(1,2)
             ff(1,6)=-ff(1,3)
             ff(2,1)=-r(2)*r(1)/dr**3
             ff(2,2)=one/dr-r(2)*r(2)/dr**3
             ff(2,3)=-r(2)*r(3)/dr**3
             ff(2,4)=-ff(2,1)
             ff(2,5)=-ff(2,2)
             ff(2,6)=-ff(2,3)
             ff(3,1)=-r(3)*r(1)/dr**3
             ff(3,2)=-r(3)*r(2)/dr**3
             ff(3,3)=one/dr-r(3)*r(3)/dr**3
             ff(3,4)=-ff(3,1)
             ff(3,5)=-ff(3,2)
             ff(3,6)=-ff(3,3)
             ff(4,1:3)=ff(1:3,4)
             ff(5,1:3)=ff(1:3,5)
             ff(6,1:3)=ff(1:3,6)
             ff(4,4:6)=ff(1,1:3)
             ff(5,4:6)=ff(2,1:3)
             ff(6,4:6)=ff(3,1:3)
             
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             do l=1,3
                do m=1,6
                   if(m <= 3) then
                      fstore=f1(m)
                      m1=k1+m
                   elseif(m > 3) then 
                      m1=k2+(m-3)
                      fstore=f2(m-3)
                   end if
                   l1=k1+l; l2=k2+l
                   if(l1<=ns .and. m1<=ns) H(l1,m1)=H(l1,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                   if(l2<=ns .and. m1<=ns) H(l2,m1)=H(l2,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
                end do
             end do

             if(calc_strain) then
                VV=(d2E_dr2-V)/(dr*dr)

                !d2E_dek_del
                if(lattice_calc) then
                   H(ns+1,ns+1)=H(ns+1,ns+1)+V*two*r(1)*r(1)+VV*r(1)*r(1)*r(1)*r(1)
                   HH=VV*r(1)*r(1)*r(2)*r(2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=VV*r(1)*r(1)*r(3)*r(3)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH
                   HH=VV*r(1)*r(1)*r(2)*r(3)
                   H(ns+1,ns+4)=H(ns+1,ns+4)+HH; H(ns+4,ns+1)=H(ns+4,ns+1)+HH
                   HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(1)*r(1)*r(1)*r(3)
                   H(ns+1,ns+5)=H(ns+1,ns+5)+HH; H(ns+5,ns+1)=H(ns+5,ns+1)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                   H(ns+1,ns+6)=H(ns+1,ns+6)+HH; H(ns+6,ns+1)=H(ns+6,ns+1)+HH

                   H(ns+2,ns+2)=H(ns+2,ns+2)+V*two*r(2)*r(2)+VV*r(2)*r(2)*r(2)*r(2)
                   HH=VV*r(2)*r(2)*r(3)*r(3)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH
                   HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(2)*r(2)*r(2)*r(3)
                   H(ns+2,ns+4)=H(ns+2,ns+4)+HH; H(ns+4,ns+2)=H(ns+4,ns+2)+HH
                   HH=VV*r(2)*r(2)*r(1)*r(3)
                   H(ns+2,ns+5)=H(ns+2,ns+5)+HH; H(ns+5,ns+2)=H(ns+5,ns+2)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                   H(ns+2,ns+6)=H(ns+2,ns+6)+HH; H(ns+6,ns+2)=H(ns+6,ns+2)+HH

                   H(ns+3,ns+3)=H(ns+3,ns+3)+V*two*r(3)*r(3)+VV*r(3)*r(3)*r(3)*r(3)
                   HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(3)*r(3)*r(3)*r(2)
                   H(ns+3,ns+4)=H(ns+3,ns+4)+HH; H(ns+4,ns+3)=H(ns+4,ns+3)+HH
                   HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(3)*r(3)*r(3)*r(1)
                   H(ns+3,ns+5)=H(ns+3,ns+5)+HH; H(ns+5,ns+3)=H(ns+5,ns+3)+HH
                   HH=VV*r(3)*r(3)*r(1)*r(2)
                   H(ns+3,ns+6)=H(ns+3,ns+6)+HH; H(ns+6,ns+3)=H(ns+6,ns+3)+HH

                   H(ns+4,ns+4)=H(ns+4,ns+4)+ &
                        V*quarter*(r(2)*r(2)+r(3)*r(3))+V*quarter*(r(2)*r(2)+r(3)*r(3))+ &
                        VV*r(2)*r(3)*r(2)*r(3)
                   HH=V*quarter*r(1)*r(2)+V*quarter*r(1)*r(2)+VV*r(2)*r(3)*r(3)*r(1)
                   H(ns+4,ns+5)=H(ns+4,ns+5)+HH; H(ns+5,ns+4)=H(ns+5,ns+4)+HH
                   HH=V*quarter*r(1)*r(3)+V*quarter*r(1)*r(3)+VV*r(2)*r(3)*r(2)*r(1)
                   H(ns+4,ns+6)=H(ns+4,ns+6)+HH; H(ns+6,ns+4)=H(ns+6,ns+4)+HH

                   H(ns+5,ns+5)=H(ns+5,ns+5)+&
                        V*quarter*(r(1)*r(1)+r(3)*r(3))+V*quarter*(r(1)*r(1)+r(3)*r(3))+ &
                        VV*r(1)*r(3)*r(1)*r(3)
                   HH=V*quarter*r(2)*r(3)+V*quarter*r(2)*r(3)+VV*r(1)*r(2)*r(1)*r(3)
                   H(ns+5,ns+6)=H(ns+5,ns+6)+HH; H(ns+6,ns+5)=H(ns+6,ns+5)+HH

                   H(ns+6,ns+6)=H(ns+6,ns+6)+&
                        V*quarter*(r(1)*r(1)+r(2)*r(2))+V*quarter*(r(1)*r(1)+r(2)*r(2))+ &
                        VV*r(1)*r(2)*r(1)*r(2)
                else if(slab_calc) then
                   H(ns+1,ns+1)=H(ns+1,ns+1)+V*two*r(1)*r(1)+VV*r(1)*r(1)*r(1)*r(1)
                   HH=VV*r(1)*r(1)*r(2)*r(2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH

                   H(ns+2,ns+2)=H(ns+2,ns+2)+V*two*r(2)*r(2)+VV*r(2)*r(2)*r(2)*r(2)
                   HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH

                   H(ns+3,ns+3)=H(ns+3,ns+3)+&
                        V*quarter*(r(1)*r(1)+r(2)*r(2))+V*quarter*(r(1)*r(1)+r(2)*r(2))+ &
                        VV*r(1)*r(2)*r(1)*r(2)
                end if

                !d2E_dri_dek
                k1=3*(ia1-1); k2=3*(ia2-1)

                if(lattice_calc) then
                   HH=V*two*r(1)+VV*r(1)*r(1)*r(1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(2)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                   HH=VV*r(1)*r(3)*r(3)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+1,ns+4)=H(k1+1,ns+4)+HH; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                   H(k2+1,ns+4)=H(k2+1,ns+4)-HH; H(ns+4,k2+1)=H(ns+4,k2+1)-HH
                   HH=V*r(3)+VV*r(1)*r(1)*r(3)
                   H(k1+1,ns+5)=H(k1+1,ns+5)+HH; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                   H(k2+1,ns+5)=H(k2+1,ns+5)-HH; H(ns+5,k2+1)=H(ns+5,k2+1)-HH
                   HH=V*r(2)+VV*r(1)*r(1)*r(2)
                   H(k1+1,ns+6)=H(k1+1,ns+6)+HH; H(ns+6,k1+1)=H(ns+6,k1+1)+HH
                   H(k2+1,ns+6)=H(k2+1,ns+6)-HH; H(ns+6,k2+1)=H(ns+6,k2+1)-HH

                   HH=VV*r(2)*r(1)*r(1)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                   HH=V*two*r(2)+VV*r(2)*r(2)*r(2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                   HH=VV*r(2)*r(3)*r(3)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH
                   HH=V*r(3)+VV*r(2)*r(2)*r(3)
                   H(k1+2,ns+4)=H(k1+2,ns+4)+HH; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                   H(k2+2,ns+4)=H(k2+2,ns+4)-HH; H(ns+4,k2+2)=H(ns+4,k2+2)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+2,ns+5)=H(k1+2,ns+5)+HH; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                   H(k2+2,ns+5)=H(k2+2,ns+5)-HH; H(ns+5,k2+2)=H(ns+5,k2+2)-HH
                   HH=V*r(1)+VV*r(2)*r(1)*r(2)
                   H(k1+2,ns+6)=H(k1+2,ns+6)+HH; H(ns+6,k1+2)=H(ns+6,k1+2)+HH
                   H(k2+2,ns+6)=H(k2+2,ns+6)-HH; H(ns+6,k2+2)=H(ns+6,k2+2)-HH

                   HH=VV*r(3)*r(1)*r(1)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                   HH=VV*r(3)*r(2)*r(2)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                   HH=V*two*r(3)+VV*r(3)*r(3)*r(3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                   HH=V*r(2)+VV*r(3)*r(2)*r(3)
                   H(k1+3,ns+4)=H(k1+3,ns+4)+HH; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                   H(k2+3,ns+4)=H(k2+3,ns+4)-HH; H(ns+4,k2+3)=H(ns+4,k2+3)-HH
                   HH=V*r(1)+VV*r(3)*r(1)*r(3)
                   H(k1+3,ns+5)=H(k1+3,ns+5)+HH; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                   H(k2+3,ns+5)=H(k2+3,ns+5)-HH; H(ns+5,k2+3)=H(ns+5,k2+3)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+3,ns+6)=H(k1+3,ns+6)+HH; H(ns+6,k1+3)=H(ns+6,k1+3)+HH
                   H(k2+3,ns+6)=H(k2+3,ns+6)-HH; H(ns+6,k2+3)=H(ns+6,k2+3)-HH
                else if(slab_calc) then
                   HH=V*two*r(1)+VV*r(1)*r(1)*r(1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                   HH=VV*r(1)*r(2)*r(2)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                   HH=V*r(2)+VV*r(1)*r(1)*r(2)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH

                   HH=VV*r(2)*r(1)*r(1)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                   HH=V*two*r(2)+VV*r(2)*r(2)*r(2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                   HH=V*r(1)+VV*r(2)*r(1)*r(2)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH

                   HH=VV*r(3)*r(1)*r(1)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                   HH=VV*r(3)*r(2)*r(2)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                   HH=VV*r(1)*r(2)*r(3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                end if
             end if
          end if
       end do j_length
    end do i_n_2b
    
!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do


  end subroutine two_body_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine three_body_E_and_F()

    integer(kind=i4_kind) :: i,j,k,length,ia1,ia2,ia3,id,l,m,m1,k1,k2,k3,l1,l2,l3,ns
    real(kind=r8_kind) :: r1(3),r2(3),r3(3),dr1,dr2,dr3,dot_prod,theta,dtheta
    real(kind=r8_kind) :: sin_th,cos_th,r12,r22,r32,r13,r23,r1_buf(3),r2_buf(3),r3_buf(3)
    real(kind=r8_kind) :: f1(3),f2(3),f3(3),dE_dt,e1(3),e2(3)
    real(kind=r8_kind) :: p(n_parameter),E_buf
    real(r8_kind) :: d2E_dt2,ff(9,9),fstore,HH
    real(r8_kind) :: vc,V1,V2,V3,A(6),rr1(6),rr2(6),rr3(6)
    real(r8_kind) :: C
    real(r8_kind) :: B1(3),B2(3),B3(2),AA(3,6),BB(3,9)

    ns=3*n_species

    i_n_3b: do i=1,n_3b
       k=n_2b+i
       length=size(three_body(i)%list,2)
       j_length: do j=1,length
          ia1=three_body(i)%list(1,j)
          ia2=three_body(i)%list(2,j)
          ia3=three_body(i)%list(3,j)
          if(ia1 <= n_species) then
             r1_buf=atoms_cart(ia1)%r
          else
             r1_buf=pc_cart(ia1-n_species)%r
          end if
          if(ia2 <= n_species) then
             r2_buf=atoms_cart(ia2)%r
          else
             r2_buf=pc_cart(ia2-n_species)%r
          end if
          if(ia3 <= n_species) then
             r3_buf=atoms_cart(ia3)%r
          else
             r3_buf=pc_cart(ia3-n_species)%r
          end if
          r1=r1_buf-r2_buf
          r2=r3_buf-r2_buf
!!$          r3=atoms_cart(ia3)%r-atoms_cart(ia1)%r
          if(lattice_calc) then 
             call image(r1)
             call image(r2)
!!$             call image(r3)
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
!!$             call image_slab(r3)
          end if
          r3=r2-r1
          r12=dot_product(r1,r1)
          r22=dot_product(r2,r2)
          r32=dot_product(r3,r3)
          dr1=sqrt(r12)
          dr2=sqrt(r22)
          dr3=sqrt(r32)
          r13=dr1*r12
          r23=dr2*r22
          dot_prod=dot_product(r1,r2)
          cos_th=dot_prod/(dr1*dr2)
          theta=acos(cos_th)

          id=poten(k)%id
          p=poten(k)%param
          select case(id)
          case (6) ! harm_bnd
             dtheta=theta-p(2)*deg2rad
             E_buf=(p(1)/two)*dtheta**2
             E(6)=E(6)+E_buf
             if(calc_gradients) dE_dt=p(1)*dtheta
             if(calc_hessian) d2E_dt2=p(1)
          case (7) ! quart_bnd
             dtheta=theta-p(4)*deg2rad
             E_buf=(p(1)/two)*dtheta**2+(p(2)/three)*dtheta**3+(p(3)/four)*dtheta**4
             E(7)=E(7)+E_buf
             if(calc_gradients) dE_dt=p(1)*dtheta+p(2)*dtheta**2+p(3)*dtheta**3
             if(calc_hessian) d2E_dt2=p(1)+two*p(2)*dtheta+three*p(3)*dtheta**2
          case (8) ! six_bnd
             dtheta=theta-p(6)*deg2rad
             E_buf=(p(1)/two)*dtheta**2+(p(2)/three)*dtheta**3+(p(3)/four)*dtheta**4+ &
                  (p(4)/five)*dtheta**5+(p(5)/six)*dtheta**6
             E(8)=E(8)+E_buf
             if(calc_gradients) dE_dt=p(1)*dtheta+p(2)*dtheta**2+p(3)*dtheta**3+ &
                  p(4)*dtheta**4+p(5)*dtheta**5
             if(calc_hessian) d2E_dt2=p(1)+two*p(2)*dtheta+three*p(3)*dtheta**2+ &
                  four*p(4)*dtheta**3+five*p(5)*dtheta**4
          case (n_poten+1 :) !user_def
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             sin_th=sin(theta)
             sin_th = sign(max(small,abs(sin_th)),sin_th)

             e1=r1/dr1
             e2=r2/dr2
             f1=(cos_th*e1-e2)/(dr1*sin_th)
             f3=(cos_th*e2-e1)/(dr2*sin_th)
             f2=((dr1-dr2*cos_th)*e1+(dr2-dr1*cos_th)*e2)/(dr1*dr2*sin_th)
             if(ia1 <= n_species) Grad(:,ia1)=Grad(:,ia1)+dE_dt*f1
             if(ia2 <= n_species) Grad(:,ia2)=Grad(:,ia2)+dE_dt*f2
             if(ia3 <= n_species) Grad(:,ia3)=Grad(:,ia3)+dE_dt*f3
             if(calc_strain) then
                vc=dE_dt/sin_th
                V1=-vc*(r12-r22+r32)/(two*r13*dr2)
                V2=-vc*(r22-r12+r32)/(two*r23*dr1)
                V3= vc/(dr1*dr2)

                if(lattice_calc) then
                   rr1(1)=r1(1)*r1(1); rr2(1)=r2(1)*r2(1); rr3(1)=r3(1)*r3(1)
                   rr1(2)=r1(2)*r1(2); rr2(2)=r2(2)*r2(2); rr3(2)=r3(2)*r3(2)
                   rr1(3)=r1(3)*r1(3); rr2(3)=r2(3)*r2(3); rr3(3)=r3(3)*r3(3)
                   rr1(4)=r1(2)*r1(3); rr2(4)=r2(2)*r2(3); rr3(4)=r3(2)*r3(3)
                   rr1(5)=r1(1)*r1(3); rr2(5)=r2(1)*r2(3); rr3(5)=r3(1)*r3(3)
                   rr1(6)=r1(1)*r1(2); rr2(6)=r2(1)*r2(2); rr3(6)=r3(1)*r3(2)

                   A(1)=V1*rr1(1)+V2*rr2(1)+V3*rr3(1)
                   A(2)=V1*rr1(2)+V2*rr2(2)+V3*rr3(2)
                   A(3)=V1*rr1(3)+V2*rr2(3)+V3*rr3(3)
                   A(4)=V1*rr1(4)+V2*rr2(4)+V3*rr3(4)
                   A(5)=V1*rr1(5)+V2*rr2(5)+V3*rr3(5)
                   A(6)=V1*rr1(6)+V2*rr2(6)+V3*rr3(6)

                   Grad_s(1)=Grad_s(1)+A(1)
                   Grad_s(2)=Grad_s(2)+A(2)
                   Grad_s(3)=Grad_s(3)+A(3)
                   Grad_s(4)=Grad_s(4)+A(4)
                   Grad_s(5)=Grad_s(5)+A(5)
                   Grad_s(6)=Grad_s(6)+A(6)
                else if(slab_calc) then
                   rr1(1)=r1(1)*r1(1); rr2(1)=r2(1)*r2(1); rr3(1)=r3(1)*r3(1)
                   rr1(2)=r1(2)*r1(2); rr2(2)=r2(2)*r2(2); rr3(2)=r3(2)*r3(2)
                   rr1(3)=r1(1)*r1(2); rr2(3)=r2(1)*r2(2); rr3(3)=r3(1)*r3(2)

                   A(1)=V1*rr1(1)+V2*rr2(1)+V3*rr3(1)
                   A(2)=V1*rr1(2)+V2*rr2(2)+V3*rr3(2)
                   A(3)=V1*rr1(3)+V2*rr2(3)+V3*rr3(3)

                   Grad_s(1)=Grad_s(1)+A(1)
                   Grad_s(2)=Grad_s(2)+A(2)
                   Grad_s(3)=Grad_s(3)+A(3)
                end if
             end if
          end if

          if(calc_hessian) then
             ff(1,1)=(-sin_th*f1(1)*e1(1)+cos_th*(one/dr1-r1(1)*r1(1)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-r1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

             ff(1,2)=(-sin_th*f1(2)*e1(1)+cos_th*(-r1(1)*r1(2)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-r1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

             ff(1,3)=(-sin_th*f1(3)*e1(1)+cos_th*(-r1(1)*r1(3)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-r1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

             ff(1,4)=(-sin_th*f2(1)*e1(1)+cos_th*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                  one/dr2-r2(1)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(r1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

             ff(1,5)=(-sin_th*f2(2)*e1(1)+cos_th*(r1(1)*r1(2)/dr1**3)- &
                  r2(1)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(r1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

             ff(1,6)=(-sin_th*f2(3)*e1(1)+cos_th*(r1(1)*r1(3)/dr1**3)- &
                  r2(1)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(r1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

             ff(1,7)=(-sin_th*f3(1)*e1(1)-one/dr2+r2(1)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-cos_th*f3(1)/(dr1*sin_th**2))

             ff(1,8)=(-sin_th*f3(2)*e1(1)+r2(1)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-cos_th*f3(2)/(dr1*sin_th**2))

             ff(1,9)=(-sin_th*f3(3)*e1(1)+r2(1)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(1)-e2(1))*(-cos_th*f3(3)/(dr1*sin_th**2))
             !........................................................
             ff(2,1)=(-sin_th*f1(1)*e1(2)+cos_th*(-r1(2)*r1(1)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-r1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

             ff(2,2)=(-sin_th*f1(2)*e1(2)+cos_th*(one/dr1-r1(2)*r1(2)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-r1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

             ff(2,3)=(-sin_th*f1(3)*e1(2)+cos_th*(-r1(2)*r1(3)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-r1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

             ff(2,4)=(-sin_th*f2(1)*e1(2)+cos_th*(r1(2)*r1(1)/dr1**3)- &
                  r2(2)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(r1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

             ff(2,5)=(-sin_th*f2(2)*e1(2)+cos_th*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                  one/dr2-r2(2)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(r1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

             ff(2,6)=(-sin_th*f2(3)*e1(2)+cos_th*(r1(2)*r1(3)/dr1**3)- &
                  r2(2)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(r1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

             ff(2,7)=(-sin_th*f3(1)*e1(2)+r2(2)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-cos_th*f3(1)/(dr1*sin_th**2))

             ff(2,8)=(-sin_th*f3(2)*e1(2)-one/dr2+r2(2)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-cos_th*f3(2)/(dr1*sin_th**2))

             ff(2,9)=(-sin_th*f3(3)*e1(2)+r2(2)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(2)-e2(2))*(-cos_th*f3(3)/(dr1*sin_th**2))
             !........................................................
             ff(3,1)=(-sin_th*f1(1)*e1(3)+cos_th*(-r1(3)*r1(1)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-r1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

             ff(3,2)=(-sin_th*f1(2)*e1(3)+cos_th*(-r1(3)*r1(2)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-r1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

             ff(3,3)=(-sin_th*f1(3)*e1(3)+cos_th*(one/dr1-r1(3)*r1(3)/dr1**3))/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-r1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

             ff(3,4)=(-sin_th*f2(1)*e1(3)+cos_th*(r1(3)*r1(1)/dr1**3)- &
                  r2(3)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(r1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

             ff(3,5)=(-sin_th*f2(2)*e1(3)+cos_th*(r1(3)*r1(2)/dr1**3)- &
                  r2(3)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(r1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

             ff(3,6)=(-sin_th*f2(3)*e1(3)+cos_th*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                  one/dr2-r2(3)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(r1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

             ff(3,7)=(-sin_th*f3(1)*e1(3)+r2(3)*r2(1)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-cos_th*f3(1)/(dr1*sin_th**2))

             ff(3,8)=(-sin_th*f3(2)*e1(3)+r2(3)*r2(2)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-cos_th*f3(2)/(dr1*sin_th**2))

             ff(3,9)=(-sin_th*f3(3)*e1(3)-one/dr2+r2(3)*r2(3)/dr2**3)/(dr1*sin_th)+ &
                  (cos_th*e1(3)-e2(3))*(-cos_th*f3(3)/(dr1*sin_th**2))
             !........................................................
             ff(7,1:3)=ff(1:3,7)

             ff(7,4)=(-sin_th*f2(1)*e2(1)+cos_th*(-one/dr2+r2(1)*r2(1)/dr2**3)+ &
                  one/dr1-r1(1)*r1(1)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(r2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

             ff(7,5)=(-sin_th*f2(2)*e2(1)+cos_th*(r2(1)*r2(2)/dr2**3)- &
                  r1(1)*r1(2)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(r2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

             ff(7,6)=(-sin_th*f2(3)*e2(1)+cos_th*(r2(1)*r2(3)/dr2**3)- &
                  r1(1)*r1(3)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(r2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

             ff(7,7)=(-sin_th*f3(1)*e2(1)+cos_th*(one/dr2-r2(1)*r2(1)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(-r2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

             ff(7,8)=(-sin_th*f3(2)*e2(1)+cos_th*(-r2(1)*r2(2)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(-r2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

             ff(7,9)=(-sin_th*f3(3)*e2(1)+cos_th*(-r2(1)*r2(3)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(1)-e1(1))*(-r2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
             !........................................................
             ff(8,1:3)=ff(1:3,8)

             ff(8,4)=(-sin_th*f2(1)*e2(2)+cos_th*(r2(2)*r2(1)/dr2**3)- &
                  r1(2)*r1(1)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(r2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

             ff(8,5)=(-sin_th*f2(2)*e2(2)+cos_th*(-one/dr2+r2(2)*r2(2)/dr2**3)+ &
                  one/dr1-r1(2)*r1(2)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(r2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

             ff(8,6)=(-sin_th*f2(3)*e2(2)+cos_th*(r2(2)*r2(3)/dr2**3)- &
                  r1(2)*r1(3)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(r2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

             ff(8,7)=(-sin_th*f3(1)*e2(2)+cos_th*(-r2(2)*r2(1)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(-r2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

             ff(8,8)=(-sin_th*f3(2)*e2(2)+cos_th*(one/dr2-r2(2)*r2(2)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(-r2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

             ff(8,9)=(-sin_th*f3(3)*e2(2)+cos_th*(-r2(2)*r2(3)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(2)-e1(2))*(-r2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
             !........................................................
             ff(9,1:3)=ff(1:3,9)

             ff(9,4)=(-sin_th*f2(1)*e2(3)+cos_th*(r2(3)*r2(1)/dr2**3)- &
                  r1(3)*r1(1)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(r2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

             ff(9,5)=(-sin_th*f2(2)*e2(3)+cos_th*(r2(3)*r2(2)/dr2**3)- &
                  r1(3)*r1(2)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(r2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

             ff(9,6)=(-sin_th*f2(3)*e2(3)+cos_th*(-one/dr2+r2(3)*r2(3)/dr2**3)+ &
                  one/dr1-r1(3)*r1(3)/dr1**3)/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(r2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

             ff(9,7)=(-sin_th*f3(1)*e2(3)+cos_th*(-r2(3)*r2(1)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(-r2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

             ff(9,8)=(-sin_th*f3(2)*e2(3)+cos_th*(-r2(3)*r2(2)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(-r2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

             ff(9,9)=(-sin_th*f3(3)*e2(3)+cos_th*(one/dr2-r2(3)*r2(3)/dr2**3))/(dr2*sin_th)+ &
                  (cos_th*e2(3)-e1(3))*(-r2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
             !........................................................
             ff(4,1:3)=ff(1:3,4)

             ff(4,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(1)+ &
                      (dr1-dr2*cos_th)*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                      (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(1)+ &
                      (dr2-dr1*cos_th)*(-one/dr2+r2(1)*r2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                      (r1(1)/(dr1**3*dr2*sin_th)+r2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

             ff(4,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(1)+ &
                      (dr1-dr2*cos_th)*(r1(1)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(1)+ &
                      (dr2-dr1*cos_th)*(r2(1)*r2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                      (r1(2)/(dr1**3*dr2*sin_th)+r2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

             ff(4,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(1)+ &
                      (dr1-dr2*cos_th)*(r1(1)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(1)+ &
                      (dr2-dr1*cos_th)*(r2(1)*r2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                      (r1(3)/(dr1**3*dr2*sin_th)+r2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))

             ff(4,7:9)=ff(7:9,4)
             !........................................................
             ff(5,1:3)=ff(1:3,5)

             ff(5,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(2)+ &
                      (dr1-dr2*cos_th)*(r1(2)*r1(1)/dr1**3)+ &
                      (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(2)+ &
                      (dr2-dr1*cos_th)*(r2(2)*r2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                      (r1(1)/(dr1**3*dr2*sin_th)+r2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

             ff(5,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(2)+ &
                      (dr1-dr2*cos_th)*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(2)+ &
                      (dr2-dr1*cos_th)*(-one/dr2+r2(2)*r2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                      (r1(2)/(dr1**3*dr2*sin_th)+r2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

             ff(5,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(2)+ &
                      (dr1-dr2*cos_th)*(r1(2)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(2)+ &
                      (dr2-dr1*cos_th)*(r2(2)*r2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                      (r1(3)/(dr1**3*dr2*sin_th)+r2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))

             ff(5,7:9)=ff(7:9,5)
             !........................................................
             ff(6,1:3)=ff(1:3,6)

             ff(6,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(3)+ &
                      (dr1-dr2*cos_th)*(r1(3)*r1(1)/dr1**3)+ &
                      (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(3)+ &
                      (dr2-dr1*cos_th)*(r2(3)*r2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                      (r1(1)/(dr1**3*dr2*sin_th)+r2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

             ff(6,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(3)+ &
                      (dr1-dr2*cos_th)*(r1(3)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(3)+ &
                      (dr2-dr1*cos_th)*(r2(3)*r2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                      (r1(2)/(dr1**3*dr2*sin_th)+r2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

             ff(6,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(3)+ &
                      (dr1-dr2*cos_th)*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(3)+ &
                      (dr2-dr1*cos_th)*(-one/dr2+r2(3)*r2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                      ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                      (r1(3)/(dr1**3*dr2*sin_th)+r2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))

             ff(6,7:9)=ff(7:9,6)
             !........................................................

             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             do l=1,3
                do m=1,9
                   if(m <= 3) then
                      fstore=f1(m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fstore=f2(m-3)
                   else if(m > 6 .and. m <= 9) then
                      m1=k3+(m-6)
                      fstore=f3(m-6)
                   end if
                   l1=k1+l; l2=k2+l; l3=k3+l
                   if(l1<=ns .and. m1<=ns) H(l1,m1)=H(l1,m1)+d2E_dt2*fstore*f1(l)+dE_dt*ff(l,m)
                   if(l2<=ns .and. m1<=ns) H(l2,m1)=H(l2,m1)+d2E_dt2*fstore*f2(l)+dE_dt*ff(l+3,m)
                   if(l3<=ns .and. m1<=ns) H(l3,m1)=H(l3,m1)+d2E_dt2*fstore*f3(l)+dE_dt*ff(l+6,m)
                end do
             end do

             if(calc_strain) then
                !d2E_dek_del
                C=(d2E_dt2*sin_th-dE_dt*cos_th)/(sin_th*dE_dt*dE_dt)

                B1(1)=(r12-three*r22+three*r32)/(two*r13*r12*dr2)
                B1(2)=(r12+r22+r32)/(two*r13*r23)
                B1(3)=-one/(r13*dr2)

                B2(1)=(r12+r22+r32)/(two*r13*r23)
                B2(2)=(r22-three*r12+three*r32)/(two*r23*r22*dr1)
                B2(3)=-one/(r23*dr1)

                B3(1)=B1(3)
                B3(2)=B2(3)

                if(lattice_calc) then
                   AA(1,1)=C*V1*A(1)+vc*(B1(1)*rr1(1)+B1(2)*rr2(1)+B1(3)*rr3(1))
                   AA(2,1)=C*V2*A(1)+vc*(B2(1)*rr1(1)+B2(2)*rr2(1)+B2(3)*rr3(1))
                   AA(3,1)=C*V3*A(1)+vc*(B3(1)*rr1(1)+B3(2)*rr2(1))
                   AA(1,2)=C*V1*A(2)+vc*(B1(1)*rr1(2)+B1(2)*rr2(2)+B1(3)*rr3(2))
                   AA(2,2)=C*V2*A(2)+vc*(B2(1)*rr1(2)+B2(2)*rr2(2)+B2(3)*rr3(2))
                   AA(3,2)=C*V3*A(2)+vc*(B3(1)*rr1(2)+B3(2)*rr2(2))
                   AA(1,3)=C*V1*A(3)+vc*(B1(1)*rr1(3)+B1(2)*rr2(3)+B1(3)*rr3(3))
                   AA(2,3)=C*V2*A(3)+vc*(B2(1)*rr1(3)+B2(2)*rr2(3)+B2(3)*rr3(3))
                   AA(3,3)=C*V3*A(3)+vc*(B3(1)*rr1(3)+B3(2)*rr2(3))
                   AA(1,4)=C*V1*A(4)+vc*(B1(1)*rr1(4)+B1(2)*rr2(4)+B1(3)*rr3(4))
                   AA(2,4)=C*V2*A(4)+vc*(B2(1)*rr1(4)+B2(2)*rr2(4)+B2(3)*rr3(4))
                   AA(3,4)=C*V3*A(4)+vc*(B3(1)*rr1(4)+B3(2)*rr2(4))
                   AA(1,5)=C*V1*A(5)+vc*(B1(1)*rr1(5)+B1(2)*rr2(5)+B1(3)*rr3(5))
                   AA(2,5)=C*V2*A(5)+vc*(B2(1)*rr1(5)+B2(2)*rr2(5)+B2(3)*rr3(5))
                   AA(3,5)=C*V3*A(5)+vc*(B3(1)*rr1(5)+B3(2)*rr2(5))
                   AA(1,6)=C*V1*A(6)+vc*(B1(1)*rr1(6)+B1(2)*rr2(6)+B1(3)*rr3(6))
                   AA(2,6)=C*V2*A(6)+vc*(B2(1)*rr1(6)+B2(2)*rr2(6)+B2(3)*rr3(6))
                   AA(3,6)=C*V3*A(6)+vc*(B3(1)*rr1(6)+B3(2)*rr2(6))

                   HH=rr1(1)*AA(1,1)+rr2(1)*AA(2,1)+rr3(1)*AA(3,1)+two*A(1)
                   H(ns+1,ns+1)=H(ns+1,ns+1)+HH
                   HH=rr1(1)*AA(1,2)+rr2(1)*AA(2,2)+rr3(1)*AA(3,2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=rr1(1)*AA(1,3)+rr2(1)*AA(2,3)+rr3(1)*AA(3,3)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH
                   HH=rr1(1)*AA(1,4)+rr2(1)*AA(2,4)+rr3(1)*AA(3,4)
                   H(ns+1,ns+4)=H(ns+1,ns+4)+HH; H(ns+4,ns+1)=H(ns+4,ns+1)+HH
                   HH=rr1(1)*AA(1,5)+rr2(1)*AA(2,5)+rr3(1)*AA(3,5)+A(5)
                   H(ns+1,ns+5)=H(ns+1,ns+5)+HH; H(ns+5,ns+1)=H(ns+5,ns+1)+HH
                   HH=rr1(1)*AA(1,6)+rr2(1)*AA(2,6)+rr3(1)*AA(3,6)+A(6)
                   H(ns+1,ns+6)=H(ns+1,ns+6)+HH; H(ns+6,ns+1)=H(ns+6,ns+1)+HH

                   HH=rr1(2)*AA(1,2)+rr2(2)*AA(2,2)+rr3(2)*AA(3,2)+two*A(2)
                   H(ns+2,ns+2)=H(ns+2,ns+2)+HH
                   HH=rr1(2)*AA(1,3)+rr2(2)*AA(2,3)+rr3(2)*AA(3,3)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH
                   HH=rr1(2)*AA(1,4)+rr2(2)*AA(2,4)+rr3(2)*AA(3,4)+A(4)
                   H(ns+2,ns+4)=H(ns+2,ns+4)+HH; H(ns+4,ns+2)=H(ns+4,ns+2)+HH
                   HH=rr1(2)*AA(1,5)+rr2(2)*AA(2,5)+rr3(2)*AA(3,5)
                   H(ns+2,ns+5)=H(ns+2,ns+5)+HH; H(ns+5,ns+2)=H(ns+5,ns+2)+HH
                   HH=rr1(2)*AA(1,6)+rr2(2)*AA(2,6)+rr3(2)*AA(3,6)+A(6)
                   H(ns+2,ns+6)=H(ns+2,ns+6)+HH; H(ns+6,ns+2)=H(ns+6,ns+2)+HH

                   HH=rr1(3)*AA(1,3)+rr2(3)*AA(2,3)+rr3(3)*AA(3,3)+two*A(3)
                   H(ns+3,ns+3)=H(ns+3,ns+3)+HH
                   HH=rr1(3)*AA(1,4)+rr2(3)*AA(2,4)+rr3(3)*AA(3,4)+A(4)
                   H(ns+3,ns+4)=H(ns+3,ns+4)+HH; H(ns+4,ns+3)=H(ns+4,ns+3)+HH
                   HH=rr1(3)*AA(1,5)+rr2(3)*AA(2,5)+rr3(3)*AA(3,5)+A(5)
                   H(ns+3,ns+5)=H(ns+3,ns+5)+HH; H(ns+5,ns+3)=H(ns+5,ns+3)+HH
                   HH=rr1(3)*AA(1,6)+rr2(3)*AA(2,6)+rr3(3)*AA(3,6)
                   H(ns+3,ns+6)=H(ns+3,ns+6)+HH; H(ns+6,ns+3)=H(ns+6,ns+3)+HH

                   HH=rr1(4)*AA(1,4)+rr2(4)*AA(2,4)+rr3(4)*AA(3,4)+half*(A(2)+A(3))
                   H(ns+4,ns+4)=H(ns+4,ns+4)+HH
                   HH=rr1(4)*AA(1,5)+rr2(4)*AA(2,5)+rr3(4)*AA(3,5)+half*A(6)
                   H(ns+4,ns+5)=H(ns+4,ns+5)+HH; H(ns+5,ns+4)=H(ns+5,ns+4)+HH
                   HH=rr1(4)*AA(1,6)+rr2(4)*AA(2,6)+rr3(4)*AA(3,6)+half*A(5)
                   H(ns+4,ns+6)=H(ns+4,ns+6)+HH; H(ns+6,ns+4)=H(ns+6,ns+4)+HH

                   HH=rr1(5)*AA(1,5)+rr2(5)*AA(2,5)+rr3(5)*AA(3,5)+half*(A(1)+A(3))
                   H(ns+5,ns+5)=H(ns+5,ns+5)+HH
                   HH=rr1(5)*AA(1,6)+rr2(5)*AA(2,6)+rr3(5)*AA(3,6)+half*A(4)
                   H(ns+5,ns+6)=H(ns+5,ns+6)+HH; H(ns+6,ns+5)=H(ns+6,ns+5)+HH

                   HH=rr1(6)*AA(1,6)+rr2(6)*AA(2,6)+rr3(6)*AA(3,6)+half*(A(1)+A(2))
                   H(ns+6,ns+6)=H(ns+6,ns+6)+HH
                else if(slab_calc) then
                   AA(1,1)=C*V1*A(1)+vc*(B1(1)*rr1(1)+B1(2)*rr2(1)+B1(3)*rr3(1))
                   AA(2,1)=C*V2*A(1)+vc*(B2(1)*rr1(1)+B2(2)*rr2(1)+B2(3)*rr3(1))
                   AA(3,1)=C*V3*A(1)+vc*(B3(1)*rr1(1)+B3(2)*rr2(1))
                   AA(1,2)=C*V1*A(2)+vc*(B1(1)*rr1(2)+B1(2)*rr2(2)+B1(3)*rr3(2))
                   AA(2,2)=C*V2*A(2)+vc*(B2(1)*rr1(2)+B2(2)*rr2(2)+B2(3)*rr3(2))
                   AA(3,2)=C*V3*A(2)+vc*(B3(1)*rr1(2)+B3(2)*rr2(2))
                   AA(1,3)=C*V1*A(3)+vc*(B1(1)*rr1(3)+B1(2)*rr2(3)+B1(3)*rr3(3))
                   AA(2,3)=C*V2*A(3)+vc*(B2(1)*rr1(3)+B2(2)*rr2(3)+B2(3)*rr3(3))
                   AA(3,3)=C*V3*A(3)+vc*(B3(1)*rr1(3)+B3(2)*rr2(3))

                   HH=rr1(1)*AA(1,1)+rr2(1)*AA(2,1)+rr3(1)*AA(3,1)+two*A(1)
                   H(ns+1,ns+1)=H(ns+1,ns+1)+HH
                   HH=rr1(1)*AA(1,2)+rr2(1)*AA(2,2)+rr3(1)*AA(3,2)
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=rr1(1)*AA(1,3)+rr2(1)*AA(2,3)+rr3(1)*AA(3,3)+A(3)
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH

                   HH=rr1(2)*AA(1,2)+rr2(2)*AA(2,2)+rr3(2)*AA(3,2)+two*A(2)
                   H(ns+2,ns+2)=H(ns+2,ns+2)+HH
                   HH=rr1(2)*AA(1,3)+rr2(2)*AA(2,3)+rr3(2)*AA(3,3)+A(3)
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH

                   HH=rr1(3)*AA(1,3)+rr2(3)*AA(2,3)+rr3(3)*AA(3,3)+half*(A(1)+A(2))
                   H(ns+3,ns+3)=H(ns+3,ns+3)+HH
                end if

                !d2E_dri_dek
                C=C*dE_dt

                BB(1,1)=C*V1*f1(1)+vc*(B1(1)*r1(1)-B1(3)*r3(1))
                BB(2,1)=C*V2*f1(1)+vc*(B2(1)*r1(1)-B2(3)*r3(1))
                BB(3,1)=C*V3*f1(1)+vc*(B3(1)*r1(1))
                BB(1,2)=C*V1*f1(2)+vc*(B1(1)*r1(2)-B1(3)*r3(2))
                BB(2,2)=C*V2*f1(2)+vc*(B2(1)*r1(2)-B2(3)*r3(2))
                BB(3,2)=C*V3*f1(2)+vc*(B3(1)*r1(2))
                BB(1,3)=C*V1*f1(3)+vc*(B1(1)*r1(3)-B1(3)*r3(3))
                BB(2,3)=C*V2*f1(3)+vc*(B2(1)*r1(3)-B2(3)*r3(3))
                BB(3,3)=C*V3*f1(3)+vc*(B3(1)*r1(3))
                BB(1,4)=C*V1*f2(1)+vc*(-B1(1)*r1(1)-B1(2)*r2(1))
                BB(2,4)=C*V2*f2(1)+vc*(-B2(1)*r1(1)-B2(2)*r2(1))
                BB(3,4)=C*V3*f2(1)+vc*(-B3(1)*r1(1)-B3(2)*r2(1))
                BB(1,5)=C*V1*f2(2)+vc*(-B1(1)*r1(2)-B1(2)*r2(2))
                BB(2,5)=C*V2*f2(2)+vc*(-B2(1)*r1(2)-B2(2)*r2(2))
                BB(3,5)=C*V3*f2(2)+vc*(-B3(1)*r1(2)-B3(2)*r2(2))
                BB(1,6)=C*V1*f2(3)+vc*(-B1(1)*r1(3)-B1(2)*r2(3))
                BB(2,6)=C*V2*f2(3)+vc*(-B2(1)*r1(3)-B2(2)*r2(3))
                BB(3,6)=C*V3*f2(3)+vc*(-B3(1)*r1(3)-B3(2)*r2(3))
                BB(1,7)=C*V1*f3(1)+vc*(B1(2)*r2(1)+B1(3)*r3(1))
                BB(2,7)=C*V2*f3(1)+vc*(B2(2)*r2(1)+B2(3)*r3(1))
                BB(3,7)=C*V3*f3(1)+vc*(B3(2)*r2(1))
                BB(1,8)=C*V1*f3(2)+vc*(B1(2)*r2(2)+B1(3)*r3(2))
                BB(2,8)=C*V2*f3(2)+vc*(B2(2)*r2(2)+B2(3)*r3(2))
                BB(3,8)=C*V3*f3(2)+vc*(B3(2)*r2(2))
                BB(1,9)=C*V1*f3(3)+vc*(B1(2)*r2(3)+B1(3)*r3(3))
                BB(2,9)=C*V2*f3(3)+vc*(B2(2)*r2(3)+B2(3)*r3(3))
                BB(3,9)=C*V3*f3(3)+vc*(B3(2)*r2(3))

                if(lattice_calc) then
                   HH=V1*two*r1(1)-V3*two*r3(1)+rr1(1)*BB(1,1)+rr2(1)*BB(2,1)+rr3(1)*BB(3,1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   HH=rr1(2)*BB(1,1)+rr2(2)*BB(2,1)+rr3(2)*BB(3,1)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   HH=rr1(3)*BB(1,1)+rr2(3)*BB(2,1)+rr3(3)*BB(3,1)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   HH=rr1(4)*BB(1,1)+rr2(4)*BB(2,1)+rr3(4)*BB(3,1)
                   H(k1+1,ns+4)=H(k1+1,ns+4)+HH; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                   HH=V1*r1(3)-V3*r3(3)+rr1(5)*BB(1,1)+rr2(5)*BB(2,1)+rr3(5)*BB(3,1)
                   H(k1+1,ns+5)=H(k1+1,ns+5)+HH; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                   HH=V1*r1(2)-V3*r3(2)+rr1(6)*BB(1,1)+rr2(6)*BB(2,1)+rr3(6)*BB(3,1)
                   H(k1+1,ns+6)=H(k1+1,ns+6)+HH; H(ns+6,k1+1)=H(ns+6,k1+1)+HH

                   HH=rr1(1)*BB(1,2)+rr2(1)*BB(2,2)+rr3(1)*BB(3,2)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   HH=V1*two*r1(2)-V3*two*r3(2)+rr1(2)*BB(1,2)+rr2(2)*BB(2,2)+rr3(2)*BB(3,2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   HH=rr1(3)*BB(1,2)+rr2(3)*BB(2,2)+rr3(3)*BB(3,2)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   HH=V1*r1(3)-V3*r3(3)+rr1(4)*BB(1,2)+rr2(4)*BB(2,2)+rr3(4)*BB(3,2)
                   H(k1+2,ns+4)=H(k1+2,ns+4)+HH; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                   HH=rr1(5)*BB(1,2)+rr2(5)*BB(2,2)+rr3(5)*BB(3,2)
                   H(k1+2,ns+5)=H(k1+2,ns+5)+HH; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                   HH=V1*r1(1)-V3*r3(1)+rr1(6)*BB(1,2)+rr2(6)*BB(2,2)+rr3(6)*BB(3,2)
                   H(k1+2,ns+6)=H(k1+2,ns+6)+HH; H(ns+6,k1+2)=H(ns+6,k1+2)+HH

                   HH=rr1(1)*BB(1,3)+rr2(1)*BB(2,3)+rr3(1)*BB(3,3)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   HH=rr1(2)*BB(1,3)+rr2(2)*BB(2,3)+rr3(2)*BB(3,3)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   HH=V1*two*r1(3)-V3*two*r3(3)+rr1(3)*BB(1,3)+rr2(3)*BB(2,3)+rr3(3)*BB(3,3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   HH=V1*r1(2)-V3*r3(2)+rr1(4)*BB(1,3)+rr2(4)*BB(2,3)+rr3(4)*BB(3,3)
                   H(k1+3,ns+4)=H(k1+3,ns+4)+HH; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                   HH=V1*r1(1)-V3*r3(1)+rr1(5)*BB(1,3)+rr2(5)*BB(2,3)+rr3(5)*BB(3,3)
                   H(k1+3,ns+5)=H(k1+3,ns+5)+HH; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                   HH=rr1(6)*BB(1,3)+rr2(6)*BB(2,3)+rr3(6)*BB(3,3)
                   H(k1+3,ns+6)=H(k1+3,ns+6)+HH; H(ns+6,k1+3)=H(ns+6,k1+3)+HH

                   HH=-V1*two*r1(1)-V2*two*r2(1)+rr1(1)*BB(1,4)+rr2(1)*BB(2,4)+rr3(1)*BB(3,4)
                   H(k2+1,ns+1)=H(k2+1,ns+1)+HH; H(ns+1,k2+1)=H(ns+1,k2+1)+HH
                   HH=rr1(2)*BB(1,4)+rr2(2)*BB(2,4)+rr3(2)*BB(3,4)
                   H(k2+1,ns+2)=H(k2+1,ns+2)+HH; H(ns+2,k2+1)=H(ns+2,k2+1)+HH
                   HH=rr1(3)*BB(1,4)+rr2(3)*BB(2,4)+rr3(3)*BB(3,4)
                   H(k2+1,ns+3)=H(k2+1,ns+3)+HH; H(ns+3,k2+1)=H(ns+3,k2+1)+HH
                   HH=rr1(4)*BB(1,4)+rr2(4)*BB(2,4)+rr3(4)*BB(3,4)
                   H(k2+1,ns+4)=H(k2+1,ns+4)+HH; H(ns+4,k2+1)=H(ns+4,k2+1)+HH
                   HH=-V1*r1(3)-V2*r2(3)+rr1(5)*BB(1,4)+rr2(5)*BB(2,4)+rr3(5)*BB(3,4)
                   H(k2+1,ns+5)=H(k2+1,ns+5)+HH; H(ns+5,k2+1)=H(ns+5,k2+1)+HH
                   HH=-V1*r1(2)-V2*r2(2)+rr1(6)*BB(1,4)+rr2(6)*BB(2,4)+rr3(6)*BB(3,4)
                   H(k2+1,ns+6)=H(k2+1,ns+6)+HH; H(ns+6,k2+1)=H(ns+6,k2+1)+HH

                   HH=rr1(1)*BB(1,5)+rr2(1)*BB(2,5)+rr3(1)*BB(3,5)
                   H(k2+2,ns+1)=H(k2+2,ns+1)+HH; H(ns+1,k2+2)=H(ns+1,k2+2)+HH
                   HH=-V1*two*r1(2)-V2*two*r2(2)+rr1(2)*BB(1,5)+rr2(2)*BB(2,5)+rr3(2)*BB(3,5)
                   H(k2+2,ns+2)=H(k2+2,ns+2)+HH; H(ns+2,k2+2)=H(ns+2,k2+2)+HH
                   HH=rr1(3)*BB(1,5)+rr2(3)*BB(2,5)+rr3(3)*BB(3,5)
                   H(k2+2,ns+3)=H(k2+2,ns+3)+HH; H(ns+3,k2+2)=H(ns+3,k2+2)+HH
                   HH=-V1*r1(3)-V2*r2(3)+rr1(4)*BB(1,5)+rr2(4)*BB(2,5)+rr3(4)*BB(3,5)
                   H(k2+2,ns+4)=H(k2+2,ns+4)+HH; H(ns+4,k2+2)=H(ns+4,k2+2)+HH
                   HH=rr1(5)*BB(1,5)+rr2(5)*BB(2,5)+rr3(5)*BB(3,5)
                   H(k2+2,ns+5)=H(k2+2,ns+5)+HH; H(ns+5,k2+2)=H(ns+5,k2+2)+HH
                   HH=-V1*r1(1)-V2*r2(1)+rr1(6)*BB(1,5)+rr2(6)*BB(2,5)+rr3(6)*BB(3,5)
                   H(k2+2,ns+6)=H(k2+2,ns+6)+HH; H(ns+6,k2+2)=H(ns+6,k2+2)+HH

                   HH=rr1(1)*BB(1,6)+rr2(1)*BB(2,6)+rr3(1)*BB(3,6)
                   H(k2+3,ns+1)=H(k2+3,ns+1)+HH; H(ns+1,k2+3)=H(ns+1,k2+3)+HH
                   HH=rr1(2)*BB(1,6)+rr2(2)*BB(2,6)+rr3(2)*BB(3,6)
                   H(k2+3,ns+2)=H(k2+3,ns+2)+HH; H(ns+2,k2+3)=H(ns+2,k2+3)+HH
                   HH=-V1*two*r1(3)-V2*two*r2(3)+rr1(3)*BB(1,6)+rr2(3)*BB(2,6)+rr3(3)*BB(3,6)
                   H(k2+3,ns+3)=H(k2+3,ns+3)+HH; H(ns+3,k2+3)=H(ns+3,k2+3)+HH
                   HH=-V1*r1(2)-V2*r2(2)+rr1(4)*BB(1,6)+rr2(4)*BB(2,6)+rr3(4)*BB(3,6)
                   H(k2+3,ns+4)=H(k2+3,ns+4)+HH; H(ns+4,k2+3)=H(ns+4,k2+3)+HH
                   HH=-V1*r1(1)-V2*r2(1)+rr1(5)*BB(1,6)+rr2(5)*BB(2,6)+rr3(5)*BB(3,6)
                   H(k2+3,ns+5)=H(k2+3,ns+5)+HH; H(ns+5,k2+3)=H(ns+5,k2+3)+HH
                   HH=rr1(6)*BB(1,6)+rr2(6)*BB(2,6)+rr3(6)*BB(3,6)
                   H(k2+3,ns+6)=H(k2+3,ns+6)+HH; H(ns+6,k2+3)=H(ns+6,k2+3)+HH

                   HH=V2*two*r2(1)+V3*two*r3(1)+rr1(1)*BB(1,7)+rr2(1)*BB(2,7)+rr3(1)*BB(3,7)
                   H(k3+1,ns+1)=H(k3+1,ns+1)+HH; H(ns+1,k3+1)=H(ns+1,k3+1)+HH
                   HH=rr1(2)*BB(1,7)+rr2(2)*BB(2,7)+rr3(2)*BB(3,7)
                   H(k3+1,ns+2)=H(k3+1,ns+2)+HH; H(ns+2,k3+1)=H(ns+2,k3+1)+HH
                   HH=rr1(3)*BB(1,7)+rr2(3)*BB(2,7)+rr3(3)*BB(3,7)
                   H(k3+1,ns+3)=H(k3+1,ns+3)+HH; H(ns+3,k3+1)=H(ns+3,k3+1)+HH
                   HH=rr1(4)*BB(1,7)+rr2(4)*BB(2,7)+rr3(4)*BB(3,7)
                   H(k3+1,ns+4)=H(k3+1,ns+4)+HH; H(ns+4,k3+1)=H(ns+4,k3+1)+HH
                   HH=V2*r2(3)+V3*r3(3)+rr1(5)*BB(1,7)+rr2(5)*BB(2,7)+rr3(5)*BB(3,7)
                   H(k3+1,ns+5)=H(k3+1,ns+5)+HH; H(ns+5,k3+1)=H(ns+5,k3+1)+HH
                   HH=V2*r2(2)+V3*r3(2)+rr1(6)*BB(1,7)+rr2(6)*BB(2,7)+rr3(6)*BB(3,7)
                   H(k3+1,ns+6)=H(k3+1,ns+6)+HH; H(ns+6,k3+1)=H(ns+6,k3+1)+HH

                   HH=rr1(1)*BB(1,8)+rr2(1)*BB(2,8)+rr3(1)*BB(3,8)
                   H(k3+2,ns+1)=H(k3+2,ns+1)+HH; H(ns+1,k3+2)=H(ns+1,k3+2)+HH
                   HH=V2*two*r2(2)+V3*two*r3(2)+rr1(2)*BB(1,8)+rr2(2)*BB(2,8)+rr3(2)*BB(3,8)
                   H(k3+2,ns+2)=H(k3+2,ns+2)+HH; H(ns+2,k3+2)=H(ns+2,k3+2)+HH
                   HH=rr1(3)*BB(1,8)+rr2(3)*BB(2,8)+rr3(3)*BB(3,8)
                   H(k3+2,ns+3)=H(k3+2,ns+3)+HH; H(ns+3,k3+2)=H(ns+3,k3+2)+HH
                   HH=V2*r2(3)+V3*r3(3)+rr1(4)*BB(1,8)+rr2(4)*BB(2,8)+rr3(4)*BB(3,8)
                   H(k3+2,ns+4)=H(k3+2,ns+4)+HH; H(ns+4,k3+2)=H(ns+4,k3+2)+HH
                   HH=rr1(5)*BB(1,8)+rr2(5)*BB(2,8)+rr3(5)*BB(3,8)
                   H(k3+2,ns+5)=H(k3+2,ns+5)+HH; H(ns+5,k3+2)=H(ns+5,k3+2)+HH
                   HH=V2*r2(1)+V3*r3(1)+rr1(6)*BB(1,8)+rr2(6)*BB(2,8)+rr3(6)*BB(3,8)
                   H(k3+2,ns+6)=H(k3+2,ns+6)+HH; H(ns+6,k3+2)=H(ns+6,k3+2)+HH

                   HH=rr1(1)*BB(1,9)+rr2(1)*BB(2,9)+rr3(1)*BB(3,9)
                   H(k3+3,ns+1)=H(k3+3,ns+1)+HH; H(ns+1,k3+3)=H(ns+1,k3+3)+HH
                   HH=rr1(2)*BB(1,9)+rr2(2)*BB(2,9)+rr3(2)*BB(3,9)
                   H(k3+3,ns+2)=H(k3+3,ns+2)+HH; H(ns+2,k3+3)=H(ns+2,k3+3)+HH
                   HH=V2*two*r2(3)+V3*two*r3(3)+rr1(3)*BB(1,9)+rr2(3)*BB(2,9)+rr3(3)*BB(3,9)
                   H(k3+3,ns+3)=H(k3+3,ns+3)+HH; H(ns+3,k3+3)=H(ns+3,k3+3)+HH
                   HH=V2*r2(2)+V3*r3(2)+rr1(4)*BB(1,9)+rr2(4)*BB(2,9)+rr3(4)*BB(3,9)
                   H(k3+3,ns+4)=H(k3+3,ns+4)+HH; H(ns+4,k3+3)=H(ns+4,k3+3)+HH
                   HH=V2*r2(1)+V3*r3(1)+rr1(5)*BB(1,9)+rr2(5)*BB(2,9)+rr3(5)*BB(3,9)
                   H(k3+3,ns+5)=H(k3+3,ns+5)+HH; H(ns+5,k3+3)=H(ns+5,k3+3)+HH
                   HH=rr1(6)*BB(1,9)+rr2(6)*BB(2,9)+rr3(6)*BB(3,9)
                   H(k3+3,ns+6)=H(k3+3,ns+6)+HH; H(ns+6,k3+3)=H(ns+6,k3+3)+HH
                else if(slab_calc) then
                   HH=V1*two*r1(1)-V3*two*r3(1)+rr1(1)*BB(1,1)+rr2(1)*BB(2,1)+rr3(1)*BB(3,1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   HH=rr1(2)*BB(1,1)+rr2(2)*BB(2,1)+rr3(2)*BB(3,1)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   HH=V1*r1(2)-V3*r3(2)+rr1(3)*BB(1,1)+rr2(3)*BB(2,1)+rr3(3)*BB(3,1)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH

                   HH=rr1(1)*BB(1,2)+rr2(1)*BB(2,2)+rr3(1)*BB(3,2)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   HH=V1*two*r1(2)-V3*two*r3(2)+rr1(2)*BB(1,2)+rr2(2)*BB(2,2)+rr3(2)*BB(3,2)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   HH=V1*r1(1)-V3*r3(1)+rr1(3)*BB(1,2)+rr2(3)*BB(2,2)+rr3(3)*BB(3,2)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH

                   HH=rr1(1)*BB(1,3)+rr2(1)*BB(2,3)+rr3(1)*BB(3,3)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   HH=rr1(2)*BB(1,3)+rr2(2)*BB(2,3)+rr3(2)*BB(3,3)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   HH=rr1(3)*BB(1,3)+rr2(3)*BB(2,3)+rr3(3)*BB(3,3)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH

                   HH=-V1*two*r1(1)-V2*two*r2(1)+rr1(1)*BB(1,4)+rr2(1)*BB(2,4)+rr3(1)*BB(3,4)
                   H(k2+1,ns+1)=H(k2+1,ns+1)+HH; H(ns+1,k2+1)=H(ns+1,k2+1)+HH
                   HH=rr1(2)*BB(1,4)+rr2(2)*BB(2,4)+rr3(2)*BB(3,4)
                   H(k2+1,ns+2)=H(k2+1,ns+2)+HH; H(ns+2,k2+1)=H(ns+2,k2+1)+HH
                   HH=-V1*r1(2)-V2*r2(2)+rr1(3)*BB(1,4)+rr2(3)*BB(2,4)+rr3(3)*BB(3,4)
                   H(k2+1,ns+3)=H(k2+1,ns+3)+HH; H(ns+3,k2+1)=H(ns+3,k2+1)+HH

                   HH=rr1(1)*BB(1,5)+rr2(1)*BB(2,5)+rr3(1)*BB(3,5)
                   H(k2+2,ns+1)=H(k2+2,ns+1)+HH; H(ns+1,k2+2)=H(ns+1,k2+2)+HH
                   HH=-V1*two*r1(2)-V2*two*r2(2)+rr1(2)*BB(1,5)+rr2(2)*BB(2,5)+rr3(2)*BB(3,5)
                   H(k2+2,ns+2)=H(k2+2,ns+2)+HH; H(ns+2,k2+2)=H(ns+2,k2+2)+HH
                   HH=-V1*r1(1)-V2*r2(1)+rr1(3)*BB(1,5)+rr2(3)*BB(2,5)+rr3(3)*BB(3,5)
                   H(k2+2,ns+3)=H(k2+2,ns+3)+HH; H(ns+3,k2+2)=H(ns+3,k2+2)+HH

                   HH=rr1(1)*BB(1,6)+rr2(1)*BB(2,6)+rr3(1)*BB(3,6)
                   H(k2+3,ns+1)=H(k2+3,ns+1)+HH; H(ns+1,k2+3)=H(ns+1,k2+3)+HH
                   HH=rr1(2)*BB(1,6)+rr2(2)*BB(2,6)+rr3(2)*BB(3,6)
                   H(k2+3,ns+2)=H(k2+3,ns+2)+HH; H(ns+2,k2+3)=H(ns+2,k2+3)+HH
                   HH=rr1(3)*BB(1,6)+rr2(3)*BB(2,6)+rr3(3)*BB(3,6)
                   H(k2+3,ns+3)=H(k2+3,ns+3)+HH; H(ns+3,k2+3)=H(ns+3,k2+3)+HH

                   HH=V2*two*r2(1)+V3*two*r3(1)+rr1(1)*BB(1,7)+rr2(1)*BB(2,7)+rr3(1)*BB(3,7)
                   H(k3+1,ns+1)=H(k3+1,ns+1)+HH; H(ns+1,k3+1)=H(ns+1,k3+1)+HH
                   HH=rr1(2)*BB(1,7)+rr2(2)*BB(2,7)+rr3(2)*BB(3,7)
                   H(k3+1,ns+2)=H(k3+1,ns+2)+HH; H(ns+2,k3+1)=H(ns+2,k3+1)+HH
                   HH=V2*r2(2)+V3*r3(2)+rr1(3)*BB(1,7)+rr2(3)*BB(2,7)+rr3(3)*BB(3,7)
                   H(k3+1,ns+3)=H(k3+1,ns+3)+HH; H(ns+3,k3+1)=H(ns+3,k3+1)+HH

                   HH=rr1(1)*BB(1,8)+rr2(1)*BB(2,8)+rr3(1)*BB(3,8)
                   H(k3+2,ns+1)=H(k3+2,ns+1)+HH; H(ns+1,k3+2)=H(ns+1,k3+2)+HH
                   HH=V2*two*r2(2)+V3*two*r3(2)+rr1(2)*BB(1,8)+rr2(2)*BB(2,8)+rr3(2)*BB(3,8)
                   H(k3+2,ns+2)=H(k3+2,ns+2)+HH; H(ns+2,k3+2)=H(ns+2,k3+2)+HH
                   HH=V2*r2(1)+V3*r3(1)+rr1(3)*BB(1,8)+rr2(3)*BB(2,8)+rr3(3)*BB(3,8)
                   H(k3+2,ns+3)=H(k3+2,ns+3)+HH; H(ns+3,k3+2)=H(ns+3,k3+2)+HH

                   HH=rr1(1)*BB(1,9)+rr2(1)*BB(2,9)+rr3(1)*BB(3,9)
                   H(k3+3,ns+1)=H(k3+3,ns+1)+HH; H(ns+1,k3+3)=H(ns+1,k3+3)+HH
                   HH=rr1(2)*BB(1,9)+rr2(2)*BB(2,9)+rr3(2)*BB(3,9)
                   H(k3+3,ns+2)=H(k3+3,ns+2)+HH; H(ns+2,k3+3)=H(ns+2,k3+3)+HH
                   HH=rr1(3)*BB(1,9)+rr2(3)*BB(2,9)+rr3(3)*BB(3,9)
                   H(k3+3,ns+3)=H(k3+3,ns+3)+HH; H(ns+3,k3+3)=H(ns+3,k3+3)+HH
                end if
             end if
          end if

       end do j_length
    end do i_n_3b
!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine three_body_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine four_body_E_and_F()

    integer(kind=i4_kind) :: i,j,i1,j1,k,m,length,ia1,ia2,ia3,ia4,id,kk,k1,k2,k3,k4,l,m1
    real(kind=r8_kind) :: r1(3),r2(3),r3(3),r4(3),r5(3),r6(3),dr1,dr2,dot_prod,phi,dphi
    real(kind=r8_kind) :: r21(3),r32(3),dr21,dr32,dr3,dr4,dr5,dr6,dr2s,da1,da2,dh1,dh2
    real(kind=r8_kind) :: sin_phi,cos_phi
    real(kind=r8_kind) :: f1(3),f2(3),f3(3),f4(3),dE_dp
    real(kind=r8_kind) :: A(4,3),B(4,3),C(4,3)
    real(kind=r8_kind) :: p(n_parameter),E_buf
    real(r8_kind) :: d2E_dp2,ff(12,12),fstore
    real(r8_kind) :: dcos(12),d_sin(12),m2121,m2132,m3232,g2121(12),g2132(12),g3232(12)
    real(r8_kind) :: aa(12,12),bb(12,12),cc(12,12),fs(3)

    i_n_4b: do i=1,n_4b
       k=n_2b+n_3b+i
       length=size(four_body(i)%list,2)
       if(length == 0) cycle i_n_4b

       j_length: do j=1,length

          ia1=four_body(i)%list(1,j)
          ia2=four_body(i)%list(2,j)
          ia3=four_body(i)%list(3,j)
          ia4=four_body(i)%list(4,j)
          r1=atoms_cart(ia2)%r-atoms_cart(ia1)%r
          r2=atoms_cart(ia3)%r-atoms_cart(ia2)%r
          r3=atoms_cart(ia4)%r-atoms_cart(ia3)%r
          r4=atoms_cart(ia3)%r-atoms_cart(ia1)%r
          r5=atoms_cart(ia4)%r-atoms_cart(ia2)%r
          r6=atoms_cart(ia4)%r-atoms_cart(ia1)%r
          if(lattice_calc) then
             call image(r1)
             call image(r2)
             call image(r3)
             call image(r4)  !!!
             call image(r5)  !!! 
             call image(r6)  !!!
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
             call image_slab(r3)
             call image_slab(r4) !!!
             call image_slab(r5) !!!
             call image_slab(r6) !!!
          end if

          dr1=dot_product(r1,r1)
          dr2=dot_product(r2,r2)
          dr3=dot_product(r3,r3)
          dr4=dot_product(r4,r4)
          dr5=dot_product(r5,r5)
          dr6=dot_product(r6,r6)
          dr2s=sqrt(dr2)
          da1=(dr1-dr4+dr2)/(two*dr2s)
          da2=(dr5-dr3+dr2)/(two*dr2s)
          dh1=dr1-da1*da1
          dh2=dr5-da2*da2
!!$          cos_phi=(dh1+dh2-dr6+(da2-da1)*(da2-da1))/(two*sqrt(dh1*dh2))

          r21=vector_product(r1,r2)
          r32=vector_product(r2,r3)
          dr21=sqrt(dot_product(r21,r21))
          dr32=sqrt(dot_product(r32,r32))
          dot_prod=dot_product(r21,r32)
          cos_phi=dot_prod/(dr21*dr32)
          if(cos_phi < mone) cos_phi=mone
          if(cos_phi > one) cos_phi=one
          phi=acos(cos_phi)

          id=poten(k)%id
          p=poten(k)%param
          select case (id)
          case (9)  ! harm_trs 
             dphi=phi-p(2)*deg2rad !may be convert to degree???
             E_buf=(p(1)/two)*dphi**2
             E(9)=E(9)+E_buf
             if(calc_gradients) dE_dp=p(1)*dphi
             if(calc_hessian) d2E_dp2=p(1)
          case (10) ! tripl_cos
             E_buf=(p(1)*(one+cos(phi))+p(2)*(one-cos(two*phi))+ &
                  p(3)*(one+cos(three*phi)))/two
             E(10)=E(10)+E_buf
             if(calc_gradients) dE_dp=(-p(1)*sin(phi)+two*p(2)*sin(two*phi)- &
                  three*p(3)*sin(three*phi))/two
             if(calc_hessian) d2E_dp2=(-p(1)*cos(phi)+four*p(2)*cos(two*phi)- &
                  nine*p(3)*cos(three*phi))/two
          case (n_poten+1 :) !user_def
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             sin_phi=sin(phi)
             sin_phi = sign(max(small,abs(sin_phi)),sin_phi)

             A(1,1)=-r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))+r3(1)*(r2(2)*r2(2)+r2(3)*r2(3))
             A(1,2)=-r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))+r3(2)*(r2(1)*r2(1)+r2(3)*r2(3))
             A(1,3)=-r2(3)*(r2(1)*r3(1)+r2(2)*r3(2))+r3(3)*(r2(1)*r2(1)+r2(2)*r2(2))

             A(2,1)=-r1(1)*(r2(2)*r3(2)+r2(3)*r3(3))+r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))- &
                  r3(1)*(r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))+ &
                  two*r2(1)*(r1(2)*r3(2)+r1(3)*r3(3))
             A(2,2)=-r1(2)*(r2(1)*r3(1)+r2(3)*r3(3))+r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))- &
                  r3(2)*(r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))+ &
                  two*r2(2)*(r1(1)*r3(1)+r1(3)*r3(3))
             A(2,3)=-r1(3)*(r2(1)*r3(1)+r2(2)*r3(2))+r2(3)*(r2(1)*r3(1)+r2(2)*r3(2))- &
                  r3(3)*(r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))+ &
                  two*r2(3)*(r1(1)*r3(1)+r1(2)*r3(2))

             A(3,1)=r1(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))- &
                  r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))+r3(1)*(r1(2)*r2(2)+r1(3)*r2(3))- &
                  two*r2(1)*(r1(2)*r3(2)+r1(3)*r3(3))
             A(3,2)=r1(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))- &
                  r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))+r3(2)*(r1(1)*r2(1)+r1(3)*r2(3))- &
                  two*r2(2)*(r1(1)*r3(1)+r1(3)*r3(3))
             A(3,3)=r1(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r2(2)*r3(2)+r2(1)*r3(1))- &
                  r2(3)*(r1(2)*r2(2)+r1(1)*r2(1))+r3(3)*(r1(2)*r2(2)+r1(1)*r2(1))- &
                  two*r2(3)*(r1(2)*r3(2)+r1(1)*r3(1))

             A(4,1)=-r1(1)*(r2(2)*r2(2)+r2(3)*r2(3))+r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))
             A(4,2)=-r1(2)*(r2(1)*r2(1)+r2(3)*r2(3))+r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))
             A(4,3)=-r1(3)*(r2(1)*r2(1)+r2(2)*r2(2))+r2(3)*(r1(1)*r2(1)+r1(2)*r2(2))

             B(1,1)=-two*r1(1)*(r2(2)*r2(2)+r2(3)*r2(3))+two*r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))
             B(1,2)=-two*r1(2)*(r2(1)*r2(1)+r2(3)*r2(3))+two*r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))
             B(1,3)=-two*r1(3)*(r2(2)*r2(2)+r2(1)*r2(1))+two*r2(3)*(r1(2)*r2(2)+r1(1)*r2(1))

             B(2,1)=two*r1(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r1(2)*r2(2)+r1(3)*r2(3))- &
                  two*r2(1)*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             B(2,2)=two*r1(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r1(1)*r2(1)+r1(3)*r2(3))- &
                  two*r2(2)*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             B(2,3)=two*r1(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r1(2)*r2(2)+r1(1)*r2(1))- &
                  two*r2(3)*(r1(2)*r1(2)+r1(1)*r1(1)+r1(2)*r2(2)+r1(1)*r2(1))

             B(3,1)=-two*r1(1)*(r1(2)*r2(2)+r1(3)*r2(3))+two*r2(1)*(r1(2)*r1(2)+r1(3)*r1(3))
             B(3,2)=-two*r1(2)*(r1(1)*r2(1)+r1(3)*r2(3))+two*r2(2)*(r1(1)*r1(1)+r1(3)*r1(3))
             B(3,3)=-two*r1(3)*(r1(2)*r2(2)+r1(1)*r2(1))+two*r2(3)*(r1(2)*r1(2)+r1(1)*r1(1))

             B(4,:)=zero

             C(1,:)=zero

             C(2,1)=two*r3(1)*(r2(2)*r3(2)+r2(3)*r3(3))-two*r2(1)*(r3(2)*r3(2)+r3(3)*r3(3))
             C(2,2)=two*r3(2)*(r2(1)*r3(1)+r2(3)*r3(3))-two*r2(2)*(r3(1)*r3(1)+r3(3)*r3(3))
             C(2,3)=two*r3(3)*(r2(2)*r3(2)+r2(1)*r3(1))-two*r2(3)*(r3(2)*r3(2)+r3(1)*r3(1))

             C(3,1)=-two*r3(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))+ &
                  two*r2(1)*(r3(2)*r3(2)+r3(3)*r3(3)+r2(2)*r3(2)+r2(3)*r3(3))
             C(3,2)=-two*r3(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))+ &
                  two*r2(2)*(r3(1)*r3(1)+r3(3)*r3(3)+r2(1)*r3(1)+r2(3)*r3(3))
             C(3,3)=-two*r3(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r2(2)*r3(2)+r2(1)*r3(1))+ &
                  two*r2(3)*(r3(2)*r3(2)+r3(1)*r3(1)+r2(2)*r3(2)+r2(1)*r3(1))

             C(4,1)=two*r3(1)*(r2(2)*r2(2)+r2(3)*r2(3))-two*r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))
             C(4,2)=two*r3(2)*(r2(1)*r2(1)+r2(3)*r2(3))-two*r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))
             C(4,3)=two*r3(3)*(r2(2)*r2(2)+r2(1)*r2(1))-two*r2(3)*(r2(2)*r3(2)+r2(1)*r3(1))

             f1=-(A(1,:)/(dr21*dr32)-(cos_phi/two)*(B(1,:)/(dr21*dr21)+C(1,:)/(dr32*dr32)))/sin_phi
             f2=-(A(2,:)/(dr21*dr32)-(cos_phi/two)*(B(2,:)/(dr21*dr21)+C(2,:)/(dr32*dr32)))/sin_phi
             f3=-(A(3,:)/(dr21*dr32)-(cos_phi/two)*(B(3,:)/(dr21*dr21)+C(3,:)/(dr32*dr32)))/sin_phi
             f4=-(A(4,:)/(dr21*dr32)-(cos_phi/two)*(B(4,:)/(dr21*dr21)+C(4,:)/(dr32*dr32)))/sin_phi

             Grad(:,ia1)=Grad(:,ia1)+dE_dp*f1
             Grad(:,ia2)=Grad(:,ia2)+dE_dp*f2
             Grad(:,ia3)=Grad(:,ia3)+dE_dp*f3
             Grad(:,ia4)=Grad(:,ia4)+dE_dp*f4
          end if
          
          if(calc_hessian) then
             dcos(1:3)=-sin_phi*f1
             dcos(4:6)=-sin_phi*f2
             dcos(7:9)=-sin_phi*f3
             dcos(10:12)=-sin_phi*f4

             d_sin(1:3)=-f1*cos_phi/sin_phi**2
             d_sin(4:6)=-f2*cos_phi/sin_phi**2
             d_sin(7:9)=-f3*cos_phi/sin_phi**2
             d_sin(10:12)=-f4*cos_phi/sin_phi**2

             m2121=one/(dr21*dr21)
             m2132=one/(dr21*dr32)
             m3232=one/(dr32*dr32)

             g2121(1:3)=-B(1,1:3)/(dr21**4)
             g2121(4:6)=-B(2,1:3)/(dr21**4)
             g2121(7:9)=-B(3,1:3)/(dr21**4)
             g2121(10:12)=-B(4,1:3)/(dr21**4)

             g3232(1:3)=-C(1,1:3)/(dr32**4)
             g3232(4:6)=-C(2,1:3)/(dr32**4)
             g3232(7:9)=-C(3,1:3)/(dr32**4)
             g3232(10:12)=-C(4,1:3)/(dr32**4)

             g2132(1:3)=-B(1,1:3)/(two*dr21**3*dr32)-C(1,1:3)/(two*dr21*dr32**3)
             g2132(4:6)=-B(2,1:3)/(two*dr21**3*dr32)-C(2,1:3)/(two*dr21*dr32**3)
             g2132(7:9)=-B(3,1:3)/(two*dr21**3*dr32)-C(3,1:3)/(two*dr21*dr32**3)
             g2132(10:12)=-B(4,1:3)/(two*dr21**3*dr32)-C(4,1:3)/(two*dr21*dr32**3)
             !......................................
             aa(1,1:3)=zero
             aa(1,4)=r2(2)*r3(2)+r2(3)*r3(3)
             aa(1,5)=-r2(1)*(-r3(2))+r3(1)*(-two*r2(2))
             aa(1,6)=-r2(1)*(-r3(3))+r3(1)*(-two*r2(3))
             aa(1,7)=-(r2(2)*r3(2)+r2(3)*r3(3))-(r2(2)*r2(2)+r2(3)*r2(3))
             aa(1,8)=-r2(1)*(r3(2)-r2(2))+r3(1)*(two*r2(2))
             aa(1,9)=-r2(1)*(r3(3)-r2(3))+r3(1)*(two*r2(3))
             aa(1,10)=r2(2)*r2(2)+r2(3)*r2(3)
             aa(1,11)=-r2(1)*(r2(2))
             aa(1,12)=-r2(1)*(r2(3))
             !......................................
             aa(2,1:3)=zero
             aa(2,4)=-r2(2)*(-r3(1))+r3(2)*(-two*r2(1))
             aa(2,5)=r2(1)*r3(1)+r2(3)*r3(3)
             aa(2,6)=-r2(2)*(-r3(3))+r3(2)*(-two*r2(3))
             aa(2,7)=-r2(2)*(r3(1)-r2(1))+r3(2)*(two*r2(1))
             aa(2,8)=-(r2(1)*r3(1)+r2(3)*r3(3))-(r2(1)*r2(1)+r2(3)*r2(3))
             aa(2,9)=-r2(2)*(r3(3)-r2(3))+r3(2)*(two*r2(3))
             aa(2,10)=-r2(2)*(r2(1))
             aa(2,11)=r2(1)*r2(1)+r2(3)*r2(3)
             aa(2,12)=-r2(2)*(r2(3))
             !......................................
             aa(3,1:3)=zero
             aa(3,4)=-r2(3)*(-r3(1))+r3(3)*(-two*r2(1))
             aa(3,5)=-r2(3)*(-r3(2))+r3(3)*(-two*r2(2))
             aa(3,6)=r2(1)*r3(1)+r2(2)*r3(2)
             aa(3,7)=-r2(3)*(r3(1)-r2(1))+r3(3)*(two*r2(1))
             aa(3,8)=-r2(3)*(r3(2)-r2(2))+r3(3)*(two*r2(2))
             aa(3,9)=-(r2(1)*r3(1)+r2(2)*r3(2))-(r2(1)*r2(1)+r2(2)*r2(2))
             aa(3,10)=-r2(3)*(r2(1))
             aa(3,11)=-r2(3)*(r2(2))
             aa(3,12)=r2(1)*r2(1)+r2(2)*r2(2)
             !......................................
             aa(4,1:3)=aa(1:3,4)
             aa(4,4)=-two*(r2(2)*r3(2)+r2(3)*r3(3))-two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(4,5)=-r1(1)*(-r3(2))+r2(1)*(-r3(2))-r3(1)*(r2(2)-r1(2)-two*r2(2))+ &
                     two*r2(1)*(r3(2))
             aa(4,6)=-r1(1)*(-r3(3))+r2(1)*(-r3(3))-r3(1)*(r2(3)-r1(3)-two*r2(3))+ &
                     two*r2(1)*(r3(3))
             aa(4,7)=(r2(2)*r3(2)+r2(3)*r3(3))+ &
                     (r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))+two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(4,8)=-r1(1)*(r3(2)-r2(2))+r2(1)*(r3(2)-r2(2))-r3(1)*(r1(2)+two*r2(2))+ &
                     two*r2(1)*(-r1(2))
             aa(4,9)=-r1(1)*(r3(3)-r2(3))+r2(1)*(r3(3)-r2(3))-r3(1)*(r1(3)+two*r2(3))+ &
                     two*r2(1)*(-r1(3))
             aa(4,10)=-(r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))
             aa(4,11)=-r1(1)*(r2(2))+r2(1)*(r2(2))+two*r2(1)*(r1(2))
             aa(4,12)=-r1(1)*(r2(3))+r2(1)*(r2(3))+two*r2(1)*(r1(3))
             !......................................
             aa(5,1:3)=aa(1:3,5)
             aa(5,4)=-r1(2)*(-r3(1))+r2(2)*(-r3(1))-r3(2)*(r2(1)-r1(1)-two*r2(1))+ &
                     two*r2(2)*(r3(1))             
             aa(5,5)=-two*(r2(1)*r3(1)+r2(3)*r3(3))-two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(5,6)=-r1(2)*(-r3(3))+r2(2)*(-r3(3))-r3(2)*(r2(3)-r1(3)-two*r2(3))+ &
                     two*r2(2)*(r3(3))             
             aa(5,7)=-r1(2)*(r3(1)-r2(1))+r2(2)*(r3(1)-r2(1))-r3(2)*(r1(1)+two*r2(1))+ &
                     two*r2(2)*(-r1(1))
             aa(5,8)=(r2(1)*r3(1)+r2(3)*r3(3))+ &
                     (r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))+two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(5,9)=-r1(2)*(r3(3)-r2(3))+r2(2)*(r3(3)-r2(3))-r3(2)*(r1(3)+two*r2(3))+ &
                     two*r2(2)*(-r1(3))
             aa(5,10)=-r1(2)*(r2(1))+r2(2)*(r2(1))+two*r2(2)*(r1(1))
             aa(5,11)=-(r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))
             aa(5,12)=-r1(2)*(r2(3))+r2(2)*(r2(3))+two*r2(2)*(r1(3))
             !......................................
             aa(6,1:3)=aa(1:3,6)
             aa(6,4)=-r1(3)*(-r3(1))+r2(3)*(-r3(1))-r3(3)*(r2(1)-r1(1)-two*r2(1))+ &
                     two*r2(3)*(r3(1))             
             aa(6,5)=-r1(3)*(-r3(2))+r2(3)*(-r3(2))-r3(3)*(r2(2)-r1(2)-two*r2(2))+ &
                     two*r2(3)*(r3(2))             
             aa(6,6)=-two*(r2(1)*r3(1)+r2(2)*r3(2))-two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(6,7)=-r1(3)*(r3(1)-r2(1))+r2(3)*(r3(1)-r2(1))-r3(3)*(r1(1)+two*r2(1))+ &
                     two*r2(3)*(-r1(1))
             aa(6,8)=-r1(3)*(r3(2)-r2(2))+r2(3)*(r3(2)-r2(2))-r3(3)*(r1(2)+two*r2(2))+ &
                     two*r2(3)*(-r1(2))
             aa(6,9)=(r2(1)*r3(1)+r2(2)*r3(2))+ &
                     (r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))+two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(6,10)=-r1(3)*(r2(1))+r2(3)*(r2(1))+two*r2(3)*(r1(1))
             aa(6,11)=-r1(3)*(r2(2))+r2(3)*(r2(2))+two*r2(3)*(r1(2))
             aa(6,12)=-(r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))
             !......................................
             aa(7,1:3)=aa(1:3,7)
             aa(7,4:6)=aa(4:6,7)
             aa(7,7)=-two*(r1(2)*r2(2)+r1(3)*r2(3))-two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(7,8)=r1(1)*(two*r2(2)+r3(2)-r2(2))-r2(1)*(r1(2))+r3(1)*(r1(2))- &
                     two*r2(1)*(-r1(2))
             aa(7,9)=r1(1)*(two*r2(3)+r3(3)-r2(3))-r2(1)*(r1(3))+r3(1)*(r1(3))- &
                     two*r2(1)*(-r1(3))
             aa(7,10)=r1(2)*r2(2)+r1(3)*r2(3)
             aa(7,11)=r1(1)*(r2(2))-two*r2(1)*(r1(2))
             aa(7,12)=r1(1)*(r2(3))-two*r2(1)*(r1(3))
             !......................................
             aa(8,1:3)=aa(1:3,8)
             aa(8,4:6)=aa(4:6,8)
             aa(8,7)=r1(2)*(two*r2(1)+r3(1)-r2(1))-r2(2)*(r1(1))+r3(2)*(r1(1))- &
                     two*r2(2)*(-r1(1))
             aa(8,8)=-two*(r1(1)*r2(1)+r1(3)*r2(3))-two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(8,9)=r1(2)*(two*r2(3)+r3(3)-r2(3))-r2(2)*(r1(3))+r3(2)*(r1(3))- &
                     two*r2(2)*(-r1(3))
             aa(8,10)=r1(2)*(r2(1))-two*r2(2)*(r1(1))
             aa(8,11)=r1(1)*r2(1)+r1(3)*r2(3)
             aa(8,12)=r1(2)*(r2(3))-two*r2(2)*(r1(3))
             !......................................
             aa(9,1:3)=aa(1:3,9)
             aa(9,4:6)=aa(4:6,9)
             aa(9,7)=r1(3)*(two*r2(1)+r3(1)-r2(1))-r2(3)*(r1(1))+r3(3)*(r1(1))- &
                     two*r2(3)*(-r1(1))
             aa(9,8)=r1(3)*(two*r2(2)+r3(2)-r2(2))-r2(3)*(r1(2))+r3(3)*(r1(2))- &
                     two*r2(3)*(-r1(2))
             aa(9,9)=-two*(r1(1)*r2(1)+r1(2)*r2(2))-two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(9,10)=r1(3)*(r2(1))-two*r2(3)*(r1(1))
             aa(9,11)=r1(3)*(r2(2))-two*r2(3)*(r1(2))
             aa(9,12)=r1(1)*r2(1)+r1(2)*r2(2)
             !......................................
             aa(10,1:3)=aa(1:3,10)
             aa(10,4:6)=aa(4:6,10)
             aa(10,7:9)=aa(7:9,10)
             aa(10,10:12)=zero
             !......................................
             aa(11,1:3)=aa(1:3,11)
             aa(11,4:6)=aa(4:6,11)
             aa(11,7:9)=aa(7:9,11)
             aa(11,10:12)=zero
             !......................................
             aa(12,1:3)=aa(1:3,12)
             aa(12,4:6)=aa(4:6,12)
             aa(12,7:9)=aa(7:9,12)
             aa(12,10:12)=zero
             !......................................
             !......................................
             bb(1,1)=two*(r2(2)*r2(2)+r2(3)*r2(3))
             bb(1,2)=two*r2(1)*(-r2(2))
             bb(1,3)=two*r2(1)*(-r2(3))
             bb(1,4)=-two*(r2(2)*r2(2)+r2(3)*r2(3))-two*(r1(2)*r2(2)+r1(3)*r2(3))
             bb(1,5)=-two*r1(1)*(-two*r2(2))+two*r2(1)*(r2(2)-r1(2))
             bb(1,6)=-two*r1(1)*(-two*r2(3))+two*r2(1)*(r2(3)-r1(3))
             bb(1,7)=two*(r1(2)*r2(2)+r1(3)*r2(3))
             bb(1,8)=-two*r1(1)*(two*r2(2))+two*r2(1)*(r1(2))
             bb(1,9)=-two*r1(1)*(two*r2(3))+two*r2(1)*(r1(3))
             bb(1,10:12)=zero
             !......................................
             bb(2,1)=two*r2(2)*(-r2(1))
             bb(2,2)=two*(r2(1)*r2(1)+r2(3)*r2(3))
             bb(2,3)=two*r2(2)*(-r2(3))
             bb(2,4)=-two*r1(2)*(-two*r2(1))+two*r2(2)*(r2(1)-r1(1))
             bb(2,5)=-two*(r2(1)*r2(1)+r2(3)*r2(3))-two*(r1(1)*r2(1)+r1(3)*r2(3))
             bb(2,6)=-two*r1(2)*(-two*r2(3))+two*r2(2)*(r2(3)-r1(3))
             bb(2,7)=-two*r1(2)*(two*r2(1))+two*r2(2)*(r1(1))
             bb(2,8)=two*(r1(1)*r2(1)+r1(3)*r2(3))
             bb(2,9)=-two*r1(2)*(two*r2(3))+two*r2(2)*(r1(3))
             bb(2,10:12)=zero
             !......................................
             bb(3,1)=two*r2(3)*(-r2(1))
             bb(3,2)=two*r2(3)*(-r2(2))
             bb(3,3)=two*(r2(1)*r2(1)+r2(2)*r2(2))
             bb(3,4)=-two*r1(3)*(-two*r2(1))+two*r2(3)*(r2(1)-r1(1))
             bb(3,5)=-two*r1(3)*(-two*r2(2))+two*r2(3)*(r2(2)-r1(2))
             bb(3,6)=-two*(r2(1)*r2(1)+r2(2)*r2(2))-two*(r1(1)*r2(1)+r1(2)*r2(2))
             bb(3,7)=-two*r1(3)*(two*r2(1))+two*r2(3)*(r1(1))
             bb(3,8)=-two*r1(3)*(two*r2(2))+two*r2(3)*(r1(2))
             bb(3,9)=two*(r1(2)*r2(2)+r1(1)*r2(1))
             bb(3,10:12)=zero
             !......................................
             bb(4,1:3)=bb(1:3,4)
             bb(4,4)=two*(r2(2)*r2(2)+r2(3)*r2(3)+r1(2)*r2(2)+r1(3)*r2(3))+ &
                     two*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             bb(4,5)=two*r1(1)*(-two*r2(2)+r2(2)-r1(2))-two*r2(1)*(two*r1(2)+r2(2)-r1(2))
             bb(4,6)=two*r1(1)*(-two*r2(3)+r2(3)-r1(3))-two*r2(1)*(two*r1(3)+r2(3)-r1(3))
             bb(4,7)=-two*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             bb(4,8)=two*r1(1)*(two*r2(2)+r1(2))-two*r2(1)*(r1(2))
             bb(4,9)=two*r1(1)*(two*r2(3)+r1(3))-two*r2(1)*(r1(3))
             bb(4,10:12)=zero
             !......................................
             bb(5,1:3)=bb(1:3,5)
             bb(5,4)=two*r1(2)*(-two*r2(1)+r2(1)-r1(1))-two*r2(2)*(two*r1(1)+r2(1)-r1(1))
             bb(5,5)=two*(r2(1)*r2(1)+r2(3)*r2(3)+r1(1)*r2(1)+r1(3)*r2(3))+ &
                     two*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             bb(5,6)=two*r1(2)*(-two*r2(3)+r2(3)-r1(3))-two*r2(2)*(two*r1(3)+r2(3)-r1(3))
             bb(5,7)=two*r1(2)*(two*r2(1)+r1(1))-two*r2(2)*(r1(1))
             bb(5,8)=-two*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             bb(5,9)=two*r1(2)*(two*r2(3)+r1(3))-two*r2(2)*(r1(3))
             bb(5,10:12)=zero
             !......................................
             bb(6,1:3)=bb(1:3,6)
             bb(6,4)=two*r1(3)*(-two*r2(1)+r2(1)-r1(1))-two*r2(3)*(two*r1(1)+r2(1)-r1(1))
             bb(6,5)=two*r1(3)*(-two*r2(2)+r2(2)-r1(2))-two*r2(3)*(two*r1(2)+r2(2)-r1(2))
             bb(6,6)=two*(r2(1)*r2(1)+r2(2)*r2(2)+r1(1)*r2(1)+r1(2)*r2(2))+ &
                     two*(r1(1)*r1(1)+r1(2)*r1(2)+r1(1)*r2(1)+r1(2)*r2(2))
             bb(6,7)=two*r1(3)*(two*r2(1)+r1(1))-two*r2(3)*(r1(1))
             bb(6,8)=two*r1(3)*(two*r2(2)+r1(2))-two*r2(3)*(r1(2))
             bb(6,9)=-two*(r1(1)*r1(1)+r1(2)*r1(2)+r1(1)*r2(1)+r1(2)*r2(2))
             bb(6,10:12)=zero
             !......................................
             bb(7,1:3)=bb(1:3,7)
             bb(7,4:6)=bb(4:6,7)
             bb(7,7)=two*(r1(2)*r1(2)+r1(3)*r1(3))
             bb(7,8)=-two*r1(1)*(r1(2))
             bb(7,9)=-two*r1(1)*(r1(3))
             bb(7,10:12)=zero
             !......................................
             bb(8,1:3)=bb(1:3,8)
             bb(8,4:6)=bb(4:6,8)
             bb(8,7)=-two*r1(2)*(r1(1))
             bb(8,8)=two*(r1(1)*r1(1)+r1(3)*r1(3))
             bb(8,9)=-two*r1(2)*(r1(3))
             bb(8,10:12)=zero
             !......................................
             bb(9,1:3)=bb(1:3,9)
             bb(9,4:6)=bb(4:6,9)
             bb(9,7)=-two*r1(3)*(r1(1))
             bb(9,8)=-two*r1(3)*(r1(2))
             bb(9,9)=two*(r1(1)*r1(1)+r1(2)*r1(2))
             bb(9,10:12)=zero
             !......................................
             bb(10:12,1:12)=zero
             !......................................
             !......................................
             cc(1:3,1:12)=zero
             !......................................
             cc(4,1:3)=zero
             cc(4,4)=two*(r3(2)*r3(2)+r3(3)*r3(3))
             cc(4,5)=two*r3(1)*(-r3(2))
             cc(4,6)=two*r3(1)*(-r3(3))
             cc(4,7)=-two*(r2(2)*r3(2)+r2(3)*r3(3))-two*(r3(2)*r3(2)+r3(3)*r3(3))
             cc(4,8)=two*r3(1)*(r3(2)-r2(2))-two*r2(1)*(-two*r3(2))
             cc(4,9)=two*r3(1)*(r3(3)-r2(3))-two*r2(1)*(-two*r3(3))
             cc(4,10)=two*(r2(2)*r3(2)+r2(3)*r3(3))
             cc(4,11)=two*r3(1)*(r2(2))-two*r2(1)*(two*r3(2))
             cc(4,12)=two*r3(1)*(r2(3))-two*r2(1)*(two*r3(3))
             !......................................
             cc(5,1:3)=zero
             cc(5,4)=two*r3(2)*(-r3(1))
             cc(5,5)=two*(r3(1)*r3(1)+r3(3)*r3(3))
             cc(5,6)=two*r3(2)*(-r3(3))
             cc(5,7)=two*r3(2)*(r3(1)-r2(1))-two*r2(2)*(-two*r3(1))
             cc(5,8)=-two*(r2(1)*r3(1)+r2(3)*r3(3))-two*(r3(1)*r3(1)+r3(3)*r3(3))
             cc(5,9)=two*r3(2)*(r3(3)-r2(3))-two*r2(2)*(-two*r3(3))
             cc(5,10)=two*r3(2)*(r2(1))-two*r2(2)*(two*r3(1))
             cc(5,11)=two*(r2(1)*r3(1)+r2(3)*r3(3))
             cc(5,12)=two*r3(2)*(r2(3))-two*r2(2)*(two*r3(3))
             !......................................
             cc(6,1:3)=zero
             cc(6,4)=two*r3(3)*(-r3(1))
             cc(6,5)=two*r3(3)*(-r3(2))
             cc(6,6)=two*(r3(1)*r3(1)+r3(2)*r3(2))
             cc(6,7)=two*r3(3)*(r3(1)-r2(1))-two*r2(3)*(-two*r3(1))
             cc(6,8)=two*r3(3)*(r3(2)-r2(2))-two*r2(3)*(-two*r3(2))
             cc(6,9)=-two*(r2(1)*r3(1)+r2(2)*r3(2))-two*(r3(1)*r3(1)+r3(2)*r3(2))
             cc(6,10)=two*r3(3)*(r2(1))-two*r2(3)*(two*r3(1))
             cc(6,11)=two*r3(3)*(r2(2))-two*r2(3)*(two*r3(2))
             cc(6,12)=two*(r2(1)*r3(1)+r2(2)*r3(2))
             !......................................
             cc(7,1:3)=zero
             cc(7,4:6)=cc(4:6,7)
             cc(7,7)=two*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))+ &
                     two*(r3(2)*r3(2)+r3(3)*r3(3)+r2(2)*r3(2)+r2(3)*r3(3))
             cc(7,8)=-two*r3(1)*(two*r2(2)+r3(2)-r2(2))+two*r2(1)*(-two*r3(2)+r3(2)-r2(2))
             cc(7,9)=-two*r3(1)*(two*r2(3)+r3(3)-r2(3))+two*r2(1)*(-two*r3(3)+r3(3)-r2(3))
             cc(7,10)=-two*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))
             cc(7,11)=-two*r3(1)*(r2(2))+two*r2(1)*(two*r3(2)+r2(2))
             cc(7,12)=-two*r3(1)*(r2(3))+two*r2(1)*(two*r3(3)+r2(3))
             !......................................
             cc(8,1:3)=zero
             cc(8,4:6)=cc(4:6,8)
             cc(8,7)=-two*r3(2)*(two*r2(1)+r3(1)-r2(1))+two*r2(2)*(-two*r3(1)+r3(1)-r2(1))
             cc(8,8)=two*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))+ &
                     two*(r3(1)*r3(1)+r3(3)*r3(3)+r2(1)*r3(1)+r2(3)*r3(3))
             cc(8,9)=-two*r3(2)*(two*r2(3)+r3(3)-r2(3))+two*r2(2)*(-two*r3(3)+r3(3)-r2(3))
             cc(8,10)=-two*r3(2)*(r2(1))+two*r2(2)*(two*r3(1)+r2(1))
             cc(8,11)=-two*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))
             cc(8,12)=-two*r3(2)*(r2(3))+two*r2(2)*(two*r3(3)+r2(3))
             !......................................
             cc(9,1:3)=zero
             cc(9,4:6)=cc(4:6,9)
             cc(9,7)=-two*r3(3)*(two*r2(1)+r3(1)-r2(1))+two*r2(3)*(-two*r3(1)+r3(1)-r2(1))
             cc(9,8)=-two*r3(3)*(two*r2(2)+r3(2)-r2(2))+two*r2(3)*(-two*r3(2)+r3(2)-r2(2))
             cc(9,9)=two*(r2(1)*r2(1)+r2(2)*r2(2)+r2(1)*r3(1)+r2(2)*r3(2))+ &
                     two*(r3(1)*r3(1)+r3(2)*r3(2)+r2(1)*r3(1)+r2(2)*r3(2))
             cc(9,10)=-two*r3(3)*(r2(1))+two*r2(3)*(two*r3(1)+r2(1))
             cc(9,11)=-two*r3(3)*(r2(2))+two*r2(3)*(two*r3(2)+r2(2))
             cc(9,12)=-two*(r2(1)*r2(1)+r2(2)*r2(2)+r2(1)*r3(1)+r2(2)*r3(2))
             !......................................
             cc(10,1:3)=zero
             cc(10,4:6)=cc(4:6,10)
             cc(10,7:9)=cc(7:9,10)
             cc(10,10)=two*(r2(2)*r2(2)+r2(3)*r2(3))
             cc(10,11)=-two*r2(1)*(r2(2))
             cc(10,12)=-two*r2(1)*(r2(3))
             !......................................
             cc(11,1:3)=zero
             cc(11,4:6)=cc(4:6,11)
             cc(11,7:9)=cc(7:9,11)
             cc(11,10)=-two*r2(2)*(r2(1))
             cc(11,11)=two*(r2(1)*r2(1)+r2(3)*r2(3))
             cc(11,12)=-two*r2(2)*(r2(3))
             !......................................
             cc(12,1:3)=zero
             cc(12,4:6)=cc(4:6,12)
             cc(12,7:9)=cc(7:9,12)
             cc(12,10)=-two*r2(3)*(r2(1))
             cc(12,11)=-two*r2(3)*(r2(2))
             cc(12,12)=two*(r2(1)*r2(1)+r2(2)*r2(2))
             !......................................

             do i1=1,12
                kk=int((i1-1)/3)+1
                m=i1-3*(kk-1)
                if(kk==1) then
                   fs=f1
                else if(kk==2) then
                   fs=f2
                else if(kk==3) then
                   fs=f3
                else if(kk==4) then
                   fs=f4
                end if
                do j1=1,12
                   ff(i1,j1)= fs(m)*sin_phi*d_sin(j1)- &
                        (aa(i1,j1)*m2132+A(kk,m)*g2132(j1)- &
                        (dcos(j1)/two)*(B(kk,m)*m2121+C(kk,m)*m3232)- &
                        (cos_phi/two)*(bb(i1,j1)*m2121+B(kk,m)*g2121(j1)+ &
                        cc(i1,j1)*m3232+C(kk,m)*g3232(j1)))/sin_phi
                end do
             end do

             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             k4=3*(ia4-1)
             do l=1,3
                do m=1,12
                   if(m <= 3) then
                      fstore=f1(m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fstore=f2(m-3)
                   else if(m > 6 .and. m <= 9) then
                      m1=k3+(m-6)
                      fstore=f3(m-6)
                   else if(m > 9 .and. m <= 12) then
                      m1=k4+(m-9)
                      fstore=f4(m-9)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+d2E_dp2*fstore*f1(l)+dE_dp*ff(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+d2E_dp2*fstore*f2(l)+dE_dp*ff(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+d2E_dp2*fstore*f3(l)+dE_dp*ff(l+6,m)
                   H(k4+l,m1)=H(k4+l,m1)+d2E_dp2*fstore*f4(l)+dE_dp*ff(l+9,m)
                end do
             end do

          end if
             
       end do j_length
    end do i_n_4b
!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine four_body_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine many_body_E_and_F()

    integer(kind=i4_kind) :: i,j,k,length,ia1,ia2,ia3,ia4,ia5,id,i1,j1
    real(kind=r8_kind) :: r1(3),r2(3),r3(3),r4(3),dr1,dr2,dr3,dr4,dr1a,dr2a
    real(kind=r8_kind) :: sin_th1,sin_th2,e1(3),e2(3),e3(3),e4(3)
    real(kind=r8_kind) :: dot_prod1,dot_prod2,cos_th1,cos_th2,theta1,theta2,dtheta1,dtheta2
    real(r8_kind) :: f_1(3,3),f_2(3,3),f_3(3,3),fb1(5,3),fb2(5,3),ft1(4,3),ft2(4,3)
    real(kind=r8_kind) :: dE_dr1,dE_dr2,E_buf,dE_dt1,dE_dt2,dE_dp
    real(kind=r8_kind) :: r21(3),r32(3),dr21,dr32,cos_phi,sin_phi,phi
    real(kind=r8_kind) :: A(4,3),B(4,3),C(4,3)
    real(kind=r8_kind) :: p(n_parameter)
    integer(i4_kind) :: kk,k1,k2,k3,k4,k5,l,m,m1
    real(r8_kind) :: d2E_dr12,d2E_dr1dr2,d2E_dr2dr1,d2E_dr22
    real(r8_kind) :: d2E_dt12,d2E_dt1dt2,d2E_dt2dt1,d2E_dt22
    real(r8_kind) :: d2E_dr1dt1,d2E_dr2dt1,d2E_dt1dr1,d2E_dt1dr2
    real(r8_kind) :: d2E_dr2dp,d2E_dpdr2,d2E_dp2
    real(r8_kind) :: ff1(9,9),ff2(9,9),ff3(9,9),ffb1(15,15),ffb2(15,15),fft1(12,12),fft2(12,12)
    real(r8_kind) :: fst_1,fst_2,fst_3
    real(r8_kind) :: dcos(12),d_sin(12),m2121,m2132,m3232,g2121(12),g2132(12),g3232(12)
    real(r8_kind) :: aa(12,12),bb(12,12),cc(12,12),fs(3)
    real(r8_kind) :: rad2deg1

    if(trim(angle_unit) == "DEGREE") rad2deg1=rad2deg
    if(trim(angle_unit) == "RADIAN") rad2deg1=one

    i_n_ss: do i=1,n_ss
       k=n_2b+n_3b+n_4b+i
       length=size(str_str(i)%list,2)

       j_lenght: do j=1,length
          ia1=str_str(i)%list(1,j)
          ia2=str_str(i)%list(2,j)
          ia3=str_str(i)%list(3,j)
          r1=atoms_cart(ia1)%r-atoms_cart(ia2)%r
          r2=atoms_cart(ia3)%r-atoms_cart(ia2)%r
          if(lattice_calc) then 
             call image(r1)
             call image(r2)
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
          end if
          dr1=sqrt(dot_product(r1,r1))
          dr2=sqrt(dot_product(r2,r2))

          id=poten(k)%id
          p=poten(k)%param
          select case(id)
          case(11) ! str_str
             dr1a=dr1-p(2)
             dr2a=dr2-p(3)
             E_buf=(p(1)/two)*dr1a*dr2a
             E(11)=E(11)+E_buf
             if(calc_gradients) then
                dE_dr1=(p(1)/two)*dr2a
                dE_dr2=(p(1)/two)*dr1a
             end if
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             f_1(1,1:3)=r1/dr1
             f_1(2,1:3)=-r1/dr1
             f_1(3,1:3)=zero
             f_2(1,1:3)=zero
             f_2(2,1:3)=-r2/dr2
             f_2(3,1:3)=r2/dr2
             Grad(:,ia1)=Grad(:,ia1)+dE_dr1*f_1(1,:)+dE_dr2*f_2(1,:)
             Grad(:,ia2)=Grad(:,ia2)+dE_dr1*f_1(2,:)+dE_dr2*f_2(2,:)
             Grad(:,ia3)=Grad(:,ia3)+dE_dr1*f_1(3,:)+dE_dr2*f_2(3,:)
          end if

          if(calc_hessian) then
             d2E_dr12=zero
             d2E_dr1dr2=p(1)/two
             d2E_dr2dr1=p(1)/two
             d2E_dr22=zero

             ff1(1,1)=one/dr1-r1(1)*r1(1)/dr1**3
             ff1(1,2)=-r1(1)*r1(2)/dr1**3
             ff1(1,3)=-r1(1)*r1(3)/dr1**3
             ff1(1,4:6)=-ff1(1,1:3)
             ff1(1,7:9)=zero
             !...................................
             ff1(2,1)=ff1(1,2)
             ff1(2,2)=one/dr1-r1(2)*r1(2)/dr1**3
             ff1(2,3)=-r1(2)*r1(3)/dr1**3
             ff1(2,4:6)=-ff1(2,1:3)
             ff1(2,7:9)=zero
             !...................................
             ff1(3,1:2)=ff1(1:2,3)
             ff1(3,3)=one/dr1-r1(3)*r1(3)/dr1**3
             ff1(3,4:6)=-ff1(3,1:3)
             ff1(3,7:9)=zero
             !...................................
             ff1(4,1:3)=ff1(1:3,4)
             ff1(4,4:6)=ff1(1,1:3)
             ff1(4,7:9)=zero
             !...................................
             ff1(5,1:3)=ff1(1:3,5)
             ff1(5,4:6)=ff1(2,1:3)
             ff1(5,7:9)=zero
             !...................................
             ff1(6,1:3)=ff1(1:3,6)
             ff1(6,4:6)=ff1(3,1:3)
             ff1(6,7:9)=zero
             !...................................
             ff1(7:9,1:9)=zero
             !...................................
             !...................................
             ff2(1:3,1:9)=zero
             !...................................
             ff2(4,1:3)=ff2(1:3,4)
             ff2(4,4)=one/dr2-r2(1)*r2(1)/dr2**3
             ff2(4,5)=-r2(1)*r2(2)/dr2**3
             ff2(4,6)=-r2(1)*r2(3)/dr2**3
             ff2(4,7:9)=-ff2(4,4:6)
             !...................................
             ff2(5,1:3)=ff2(1:3,5)
             ff2(5,4)=ff2(4,5)
             ff2(5,5)=one/dr2-r2(2)*r2(2)/dr2**3
             ff2(5,6)=-r2(2)*r2(3)/dr2**3
             ff2(5,7:9)=-ff2(5,4:6)
             !...................................
             ff2(6,1:3)=ff2(1:3,6)
             ff2(6,4:5)=ff2(4:5,6)
             ff2(6,6)=one/dr2-r2(3)*r2(3)/dr2**3
             ff2(6,7:9)=-ff2(6,4:6)
             !...................................
             ff2(7,1:6)=ff2(1:6,7)
             ff2(7,7:9)=ff2(4,4:6)
             !...................................
             ff2(8,1:6)=ff2(1:6,8)
             ff2(8,7:9)=ff2(5,4:6)
             !...................................
             ff2(9,1:6)=ff2(1:6,9)
             ff2(9,7:9)=ff2(6,4:6)
             !...................................
             !...................................
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             do l=1,3
                do m=1,9
                   if(m <= 3) then
                      fst_1=f_1(1,m)
                      fst_2=f_2(1,m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fst_1=f_1(2,m-3)
                      fst_2=f_2(2,m-3)
                   else if(m > 6) then
                      m1=k3+(m-6)
                      fst_1=f_1(3,m-6)
                      fst_2=f_2(3,m-6)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2)*f_1(1,l)+ &
                        (d2E_dr2dr1*fst_1+d2E_dr22*fst_2)*f_2(1,l)+dE_dr1*ff1(l,m)+dE_dr2*ff2(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2)*f_1(2,l)+ &
                        (d2E_dr2dr1*fst_1+d2E_dr22*fst_2)*f_2(2,l)+dE_dr1*ff1(l+3,m)+dE_dr2*ff2(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2)*f_1(3,l)+ &
                        (d2E_dr2dr1*fst_1+d2E_dr22*fst_2)*f_2(3,l)+dE_dr1*ff1(l+6,m)+dE_dr2*ff2(l+6,m)
                end do
             end do
          end if

       end do j_lenght
    end do i_n_ss
!!$print*,'111111111111'

    i_n_bb: do i=1,n_bb
       k=n_2b+n_3b+n_4b+n_ss+i
       length=size(bnd_bnd(i)%list,2)

       j_lenght1: do j=1,length
          ia1=bnd_bnd(i)%list(1,j)
          ia2=bnd_bnd(i)%list(2,j)
          ia3=bnd_bnd(i)%list(3,j)
          ia4=bnd_bnd(i)%list(4,j)
          ia5=bnd_bnd(i)%list(5,j)
          r1=atoms_cart(ia1)%r-atoms_cart(ia3)%r
          r2=atoms_cart(ia2)%r-atoms_cart(ia3)%r
          r3=atoms_cart(ia4)%r-atoms_cart(ia3)%r
          r4=atoms_cart(ia5)%r-atoms_cart(ia3)%r
          if(lattice_calc) then 
             call image(r1)
             call image(r2)
             call image(r3)
             call image(r4)
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
             call image_slab(r3)
             call image_slab(r4)
          end if
          dr1=sqrt(dot_product(r1,r1))
          dr2=sqrt(dot_product(r2,r2))
          dr3=sqrt(dot_product(r3,r3))
          dr4=sqrt(dot_product(r4,r4))
          dot_prod1=dot_product(r1,r2)
          dot_prod2=dot_product(r3,r4)
          cos_th1=dot_prod1/(dr1*dr2)
          cos_th2=dot_prod2/(dr3*dr4)
          theta1=acos(cos_th1)
          theta2=acos(cos_th2)

          id=poten(k)%id
          p=poten(k)%param
          select case(id)
          case (12) ! bnd_bnd
             dtheta1=(theta1-p(2)*deg2rad)
             dtheta2=(theta2-p(3)*deg2rad)
             E_buf=(p(1)/two)*dtheta1*dtheta2*rad2deg1*rad2deg1
             E(12)=E(12)+E_buf
!!$print*,ia1,ia2,ia3,ia4,ia5,E_buf*j2c
             if(calc_gradients) then
                dE_dt1=(p(1)/two)*dtheta2*rad2deg1*rad2deg1
                dE_dt2=(p(1)/two)*dtheta1*rad2deg1*rad2deg1
             end if
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             sin_th1=sin(theta1)
             sin_th1 = sign(max(small,abs(sin_th1)),sin_th1)
             sin_th2=sin(theta2)
             sin_th2 = sign(max(small,abs(sin_th2)),sin_th2)

             e1=r1/dr1
             e2=r2/dr2
             e3=r3/dr3
             e4=r4/dr4
             fb1(1,1:3)=(cos_th1*e1-e2)/(dr1*sin_th1)
             fb1(2,1:3)=(cos_th1*e2-e1)/(dr2*sin_th1)
             fb1(3,1:3)=((dr1-dr2*cos_th1)*e1+(dr2-dr1*cos_th1)*e2)/(dr1*dr2*sin_th1)
             fb1(4:5,1:3)=zero
             fb2(1:2,1:3)=zero
             fb2(3,1:3)=((dr3-dr4*cos_th2)*e3+(dr4-dr3*cos_th2)*e4)/(dr4*dr3*sin_th2)
             fb2(4,1:3)=(cos_th2*e3-e4)/(dr3*sin_th2)
             fb2(5,1:3)=(cos_th2*e4-e3)/(dr4*sin_th2)
             
             Grad(:,ia1)=Grad(:,ia1)+dE_dt1*fb1(1,:)+dE_dt2*fb2(1,:)
             Grad(:,ia2)=Grad(:,ia2)+dE_dt1*fb1(2,:)+dE_dt2*fb2(2,:)
             Grad(:,ia3)=Grad(:,ia3)+dE_dt1*fb1(3,:)+dE_dt2*fb2(3,:)
             Grad(:,ia4)=Grad(:,ia4)+dE_dt1*fb1(4,:)+dE_dt2*fb2(4,:)
             Grad(:,ia5)=Grad(:,ia5)+dE_dt1*fb1(5,:)+dE_dt2*fb2(5,:)
          end if

          if(calc_hessian) then
             d2E_dt12=zero
             d2E_dt1dt2=(p(1)/two)*rad2deg1*rad2deg1
             d2E_dt2dt1=(p(1)/two)*rad2deg1*rad2deg1
             d2E_dt22=zero
             !.......................................
             ffb1(1,1)=(-sin_th1*fb1(1,1)*e1(1)+cos_th1*(one/dr1-r1(1)*r1(1)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(1)/(dr1**3*sin_th1)-cos_th1*fb1(1,1)/(dr1*sin_th1**2))
             ffb1(1,2)=(-sin_th1*fb1(1,2)*e1(1)+cos_th1*(-r1(1)*r1(2)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(2)/(dr1**3*sin_th1)-cos_th1*fb1(1,2)/(dr1*sin_th1**2))
             ffb1(1,3)=(-sin_th1*fb1(1,3)*e1(1)+cos_th1*(-r1(1)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(1,3)/(dr1*sin_th1**2))
             ffb1(1,4)=(-sin_th1*fb1(2,1)*e1(1)-one/dr2+r2(1)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*fb1(2,1)/(dr1*sin_th1**2))
             ffb1(1,5)=(-sin_th1*fb1(2,2)*e1(1)+r2(1)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*fb1(2,2)/(dr1*sin_th1**2))
             ffb1(1,6)=(-sin_th1*fb1(2,3)*e1(1)+r2(1)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*fb1(2,3)/(dr1*sin_th1**2))
             ffb1(1,7)=(-sin_th1*fb1(3,1)*e1(1)+cos_th1*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                  one/dr2-r2(1)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(1)/(dr1**3*sin_th1)-cos_th1*fb1(3,1)/(dr1*sin_th1**2))
             ffb1(1,8)=(-sin_th1*fb1(3,2)*e1(1)+cos_th1*(r1(1)*r1(2)/dr1**3)- &
                  r2(1)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(2)/(dr1**3*sin_th1)-cos_th1*fb1(3,2)/(dr1*sin_th1**2))
             ffb1(1,9)=(-sin_th1*fb1(3,3)*e1(1)+cos_th1*(r1(1)*r1(3)/dr1**3)- &
                  r2(1)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*sin_th1**2))
             ffb1(1,10:15)=zero
             !.......................................
             ffb1(2,1)=ffb1(1,2)
             ffb1(2,2)=(-sin_th1*fb1(1,2)*e1(2)+cos_th1*(one/dr1-r1(2)*r1(2)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-r1(2)/(dr1**3*sin_th1)-cos_th1*fb1(1,2)/(dr1*sin_th1**2))
             ffb1(2,3)=(-sin_th1*fb1(1,3)*e1(2)+cos_th1*(-r1(2)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(1,3)/(dr1*sin_th1**2))
             ffb1(2,4)=(-sin_th1*fb1(2,1)*e1(2)+r2(2)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*fb1(2,1)/(dr1*sin_th1**2))
             ffb1(2,5)=(-sin_th1*fb1(2,2)*e1(2)-one/dr2+r2(2)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*fb1(2,2)/(dr1*sin_th1**2))
             ffb1(2,6)=(-sin_th1*fb1(2,3)*e1(2)+r2(2)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*fb1(2,3)/(dr1*sin_th1**2))
             ffb1(2,7)=(-sin_th1*fb1(3,1)*e1(2)+cos_th1*(r1(2)*r1(1)/dr1**3)- &
                  r2(2)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(1)/(dr1**3*sin_th1)-cos_th1*fb1(3,1)/(dr1*sin_th1**2))
             ffb1(2,8)=(-sin_th1*fb1(3,2)*e1(2)+cos_th1*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                  one/dr2-r2(2)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(2)/(dr1**3*sin_th1)-cos_th1*fb1(3,2)/(dr1*sin_th1**2))
             ffb1(2,9)=(-sin_th1*fb1(3,3)*e1(2)+cos_th1*(r1(2)*r1(3)/dr1**3)- &
                  r2(2)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*sin_th1**2))
             ffb1(2,10:15)=zero
             !.......................................
             ffb1(3,1:2)=ffb1(1:2,3)
             ffb1(3,3)=(-sin_th1*fb1(1,3)*e1(3)+cos_th1*(one/dr1-r1(3)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(1,3)/(dr1*sin_th1**2))
             ffb1(3,4)=(-sin_th1*fb1(2,1)*e1(3)+r2(3)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*fb1(2,1)/(dr1*sin_th1**2))
             ffb1(3,5)=(-sin_th1*fb1(2,2)*e1(3)+r2(3)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*fb1(2,2)/(dr1*sin_th1**2))
             ffb1(3,6)=(-sin_th1*fb1(2,3)*e1(3)-one/dr2+r2(3)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*fb1(2,3)/(dr1*sin_th1**2))
             ffb1(3,7)=(-sin_th1*fb1(3,1)*e1(3)+cos_th1*(r1(3)*r1(1)/dr1**3)- &
                  r2(3)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(1)/(dr1**3*sin_th1)-cos_th1*fb1(3,1)/(dr1*sin_th1**2))
             ffb1(3,8)=(-sin_th1*fb1(3,2)*e1(3)+cos_th1*(r1(3)*r1(2)/dr1**3)- &
                  r2(3)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(2)/(dr1**3*sin_th1)-cos_th1*fb1(3,2)/(dr1*sin_th1**2))
             ffb1(3,9)=(-sin_th1*fb1(3,3)*e1(3)+cos_th1*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                  one/dr2-r2(3)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(3)/(dr1**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*sin_th1**2))
             ffb1(3,10:15)=zero
             !.......................................
             ffb1(4,1:3)=ffb1(1:3,4)
             ffb1(4,4)=(-sin_th1*fb1(2,1)*e2(1)+cos_th1*(one/dr2-r2(1)*r2(1)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(1)/(dr2**3*sin_th1)-cos_th1*fb1(2,1)/(dr2*sin_th1**2))
             ffb1(4,5)=(-sin_th1*fb1(2,2)*e2(1)+cos_th1*(-r2(1)*r2(2)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(2)/(dr2**3*sin_th1)-cos_th1*fb1(2,2)/(dr2*sin_th1**2))
             ffb1(4,6)=(-sin_th1*fb1(2,3)*e2(1)+cos_th1*(-r2(1)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(2,3)/(dr2*sin_th1**2))
             ffb1(4,7)=(-sin_th1*fb1(3,1)*e2(1)+cos_th1*(-one/dr2+r2(1)*r2(1)/dr2**3)+ &
                  one/dr1-r1(1)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(1)/(dr2**3*sin_th1)-cos_th1*fb1(3,1)/(dr2*sin_th1**2))
             ffb1(4,8)=(-sin_th1*fb1(3,2)*e2(1)+cos_th1*(r2(1)*r2(2)/dr2**3)- &
                  r1(1)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(2)/(dr2**3*sin_th1)-cos_th1*fb1(3,2)/(dr2*sin_th1**2))
             ffb1(4,9)=(-sin_th1*fb1(3,3)*e2(1)+cos_th1*(r2(1)*r2(3)/dr2**3)- &
                  r1(1)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr2*sin_th1**2))
             ffb1(4,10:15)=zero
             !.......................................
             ffb1(5,1:4)=ffb1(1:4,5)
             ffb1(5,5)=(-sin_th1*fb1(2,2)*e2(2)+cos_th1*(one/dr2-r2(2)*r2(2)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(-r2(2)/(dr2**3*sin_th1)-cos_th1*fb1(2,2)/(dr2*sin_th1**2))
             ffb1(5,6)=(-sin_th1*fb1(2,3)*e2(2)+cos_th1*(-r2(2)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(2,3)/(dr2*sin_th1**2))
             ffb1(5,7)=(-sin_th1*fb1(3,1)*e2(2)+cos_th1*(r2(2)*r2(1)/dr2**3)- &
                  r1(2)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(1)/(dr2**3*sin_th1)-cos_th1*fb1(3,1)/(dr2*sin_th1**2))
             ffb1(5,8)=(-sin_th1*fb1(3,2)*e2(2)+cos_th1*(-one/dr2+r2(2)*r2(2)/dr2**3)+ &
                  one/dr1-r1(2)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(2)/(dr2**3*sin_th1)-cos_th1*fb1(3,2)/(dr2*sin_th1**2))
             ffb1(5,9)=(-sin_th1*fb1(3,3)*e2(2)+cos_th1*(r2(2)*r2(3)/dr2**3)- &
                  r1(2)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr2*sin_th1**2))
             ffb1(5,10:15)=zero
             !.......................................
             ffb1(6,1:5)=ffb1(1:5,6)
             ffb1(6,6)=(-sin_th1*fb1(2,3)*e2(3)+cos_th1*(one/dr2-r2(3)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(2,3)/(dr2*sin_th1**2))
             ffb1(6,7)=(-sin_th1*fb1(3,1)*e2(3)+cos_th1*(r2(3)*r2(1)/dr2**3)- &
                  r1(3)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(1)/(dr2**3*sin_th1)-cos_th1*fb1(3,1)/(dr2*sin_th1**2))
             ffb1(6,8)=(-sin_th1*fb1(3,2)*e2(3)+cos_th1*(r2(3)*r2(2)/dr2**3)- &
                  r1(3)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(2)/(dr2**3*sin_th1)-cos_th1*fb1(3,2)/(dr2*sin_th1**2))
             ffb1(6,9)=(-sin_th1*fb1(3,3)*e2(3)+cos_th1*(-one/dr2+r2(3)*r2(3)/dr2**3)+ &
                  one/dr1-r1(3)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(3)/(dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr2*sin_th1**2))
             ffb1(6,10:15)=zero
             !.......................................
             ffb1(7,1:6)=ffb1(1:6,7)
             ffb1(7,7)=((-e1(1)+e2(1)*cos_th1+dr2*sin_th1*fb1(3,1))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                      (-e2(1)+e1(1)*cos_th1+dr1*sin_th1*fb1(3,1))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(1)*r2(1)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(1)/(dr1**3*dr2*sin_th1)+r2(1)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,1)/(dr1*dr2*sin_th1**2))
             ffb1(7,8)=((-e1(2)+e2(2)*cos_th1+dr2*sin_th1*fb1(3,2))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(r1(1)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th1+dr1*sin_th1*fb1(3,2))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(r2(1)*r2(2)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(2)/(dr1**3*dr2*sin_th1)+r2(2)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,2)/(dr1*dr2*sin_th1**2))
             ffb1(7,9)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*fb1(3,3))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(r1(1)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*fb1(3,3))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(r2(1)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*dr2*sin_th1**2))
             ffb1(7,10:15)=zero
             !.......................................
             ffb1(8,1:7)=ffb1(1:7,8)
             ffb1(8,8)=((-e1(2)+e2(2)*cos_th1+dr2*sin_th1*fb1(3,2))*e1(2)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th1+dr1*sin_th1*fb1(3,2))*e2(2)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(2)*r2(2)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(2)+(dr2-dr1*cos_th1)*e2(2))* &
                      (r1(2)/(dr1**3*dr2*sin_th1)+r2(2)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,2)/(dr1*dr2*sin_th1**2))
             ffb1(8,9)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*fb1(3,3))*e1(2)+ &
                      (dr1-dr2*cos_th1)*(r1(2)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*fb1(3,3))*e2(2)+ &
                      (dr2-dr1*cos_th1)*(r2(2)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(2)+(dr2-dr1*cos_th1)*e2(2))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*dr2*sin_th1**2))
             ffb1(8,10:15)=zero
             !.......................................
             ffb1(9,1:8)=ffb1(1:8,9)
             ffb1(9,9)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*fb1(3,3))*e1(3)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*fb1(3,3))*e2(3)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(3)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(3)+(dr2-dr1*cos_th1)*e2(3))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*fb1(3,3)/(dr1*dr2*sin_th1**2))
             ffb1(9,10:15)=zero
             !.......................................
             ffb1(10:15,1:15)=zero
             !.......................................
             !.......................................
             ffb2(1:6,1:15)=zero
             !.......................................
             ffb2(10,1:6)=zero
             ffb2(10,10)=(-sin_th2*fb2(4,1)*e3(1)+cos_th2*(one/dr3-r3(1)*r3(1)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-r3(1)/(dr3**3*sin_th2)-cos_th2*fb2(4,1)/(dr3*sin_th2**2))
             ffb2(10,11)=(-sin_th2*fb2(4,2)*e3(1)+cos_th2*(-r3(1)*r3(2)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-r3(2)/(dr3**3*sin_th2)-cos_th2*fb2(4,2)/(dr3*sin_th2**2))
             ffb2(10,12)=(-sin_th2*fb2(4,3)*e3(1)+cos_th2*(-r3(1)*r3(3)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(4,3)/(dr3*sin_th2**2))
             ffb2(10,13)=(-sin_th2*fb2(5,1)*e3(1)-one/dr4+r4(1)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-cos_th2*fb2(5,1)/(dr3*sin_th2**2))
             ffb2(10,14)=(-sin_th2*fb2(5,2)*e3(1)+r4(1)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-cos_th2*fb2(5,2)/(dr3*sin_th2**2))
             ffb2(10,15)=(-sin_th2*fb2(5,3)*e3(1)+r4(1)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(-cos_th2*fb2(5,3)/(dr3*sin_th2**2))
             ffb2(10,7)=(-sin_th2*fb2(3,1)*e3(1)+cos_th2*(-one/dr3+r3(1)*r3(1)/dr3**3)+ &
                  one/dr4-r4(1)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(r3(1)/(dr3**3*sin_th2)-cos_th2*fb2(3,1)/(dr3*sin_th2**2))
             ffb2(10,8)=(-sin_th2*fb2(3,2)*e3(1)+cos_th2*(r3(1)*r3(2)/dr3**3)- &
                  r4(1)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(r3(2)/(dr3**3*sin_th2)-cos_th2*fb2(3,2)/(dr3*sin_th2**2))
             ffb2(10,9)=(-sin_th2*fb2(3,3)*e3(1)+cos_th2*(r3(1)*r3(3)/dr3**3)- &
                  r4(1)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(1)-e4(1))*(r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*sin_th2**2))
             !.......................................
             ffb2(11,1:6)=zero
             ffb2(11,10)=ffb2(10,11)
             ffb2(11,11)=(-sin_th2*fb2(4,2)*e3(2)+cos_th2*(one/dr3-r3(2)*r3(2)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(-r3(2)/(dr3**3*sin_th2)-cos_th2*fb2(4,2)/(dr3*sin_th2**2))
             ffb2(11,12)=(-sin_th2*fb2(4,3)*e3(2)+cos_th2*(-r3(2)*r3(3)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(-r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(4,3)/(dr3*sin_th2**2))
             ffb2(11,13)=(-sin_th2*fb2(5,1)*e3(2)+r4(2)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(-cos_th2*fb2(5,1)/(dr3*sin_th2**2))
             ffb2(11,14)=(-sin_th2*fb2(5,2)*e3(2)-one/dr4+r4(2)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(-cos_th2*fb2(5,2)/(dr3*sin_th2**2))
             ffb2(11,15)=(-sin_th2*fb2(5,3)*e3(2)+r4(2)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(-cos_th2*fb2(5,3)/(dr3*sin_th2**2))
             ffb2(11,7)=(-sin_th2*fb2(3,1)*e3(2)+cos_th2*(r3(2)*r3(1)/dr3**3)- &
                  r4(2)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(r3(1)/(dr3**3*sin_th2)-cos_th2*fb2(3,1)/(dr3*sin_th2**2))
             ffb2(11,8)=(-sin_th2*fb2(3,2)*e3(2)+cos_th2*(-one/dr3+r3(2)*r3(2)/dr3**3)+ &
                  one/dr4-r4(2)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(r3(2)/(dr3**3*sin_th2)-cos_th2*fb2(3,2)/(dr3*sin_th2**2))
             ffb2(11,9)=(-sin_th2*fb2(3,3)*e3(2)+cos_th2*(r3(2)*r3(3)/dr3**3)- &
                  r4(2)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(2)-e4(2))*(r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*sin_th2**2))             
             !.......................................
             ffb2(12,1:6)=zero
             ffb2(12,10:11)=ffb2(10:11,12)
             ffb2(12,12)=(-sin_th2*fb2(4,3)*e3(3)+cos_th2*(one/dr3-r3(3)*r3(3)/dr3**3))/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(-r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(4,3)/(dr3*sin_th2**2))
             ffb2(12,13)=(-sin_th2*fb2(5,1)*e3(3)+r4(3)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(-cos_th2*fb2(5,1)/(dr3*sin_th2**2))
             ffb2(12,14)=(-sin_th2*fb2(5,2)*e3(3)+r4(3)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(-cos_th2*fb2(5,2)/(dr3*sin_th2**2))
             ffb2(12,15)=(-sin_th2*fb2(5,3)*e3(3)-one/dr4+r4(3)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(-cos_th2*fb2(5,3)/(dr3*sin_th2**2))
             ffb2(12,7)=(-sin_th2*fb2(3,1)*e3(3)+cos_th2*(r3(3)*r3(1)/dr3**3)- &
                  r4(3)*r4(1)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(r3(1)/(dr3**3*sin_th2)-cos_th2*fb2(3,1)/(dr3*sin_th2**2))
             ffb2(12,8)=(-sin_th2*fb2(3,2)*e3(3)+cos_th2*(r3(3)*r3(2)/dr3**3)- &
                  r4(3)*r4(2)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(r3(2)/(dr3**3*sin_th2)-cos_th2*fb2(3,2)/(dr3*sin_th2**2))
             ffb2(12,9)=(-sin_th2*fb2(3,3)*e3(3)+cos_th2*(-one/dr3+r3(3)*r3(3)/dr3**3)+ &
                  one/dr4-r4(3)*r4(3)/dr4**3)/(dr3*sin_th2)+ &
                  (cos_th2*e3(3)-e4(3))*(r3(3)/(dr3**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*sin_th2**2))
             !.......................................
             ffb2(13,1:6)=zero
             ffb2(13,10:12)=ffb2(10:12,13)
             ffb2(13,13)=(-sin_th2*fb2(5,1)*e4(1)+cos_th2*(one/dr4-r4(1)*r4(1)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(-r4(1)/(dr4**3*sin_th2)-cos_th2*fb2(5,1)/(dr4*sin_th2**2))
             ffb2(13,14)=(-sin_th2*fb2(5,2)*e4(1)+cos_th2*(-r4(1)*r4(2)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(-r4(2)/(dr4**3*sin_th2)-cos_th2*fb2(5,2)/(dr4*sin_th2**2))
             ffb2(13,15)=(-sin_th2*fb2(5,3)*e4(1)+cos_th2*(-r4(1)*r4(3)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(-r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(5,3)/(dr4*sin_th2**2))
             ffb2(13,7)=(-sin_th2*fb2(3,1)*e4(1)+cos_th2*(-one/dr4+r4(1)*r4(1)/dr4**3)+ &
                  one/dr3-r3(1)*r3(1)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(r4(1)/(dr4**3*sin_th2)-cos_th2*fb2(3,1)/(dr4*sin_th2**2))
             ffb2(13,8)=(-sin_th2*fb2(3,2)*e4(1)+cos_th2*(r4(1)*r4(2)/dr4**3)- &
                  r3(1)*r3(2)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(r4(2)/(dr4**3*sin_th2)-cos_th2*fb2(3,2)/(dr4*sin_th2**2))
             ffb2(13,9)=(-sin_th2*fb2(3,3)*e4(1)+cos_th2*(r4(1)*r4(3)/dr4**3)- &
                  r3(1)*r3(3)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(1)-e3(1))*(r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr4*sin_th2**2))
             !.......................................
             ffb2(14,1:6)=zero
             ffb2(14,10:13)=ffb2(10:13,14)
             ffb2(14,14)=(-sin_th2*fb2(5,2)*e4(2)+cos_th2*(one/dr4-r4(2)*r4(2)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(2)-e3(2))*(-r4(2)/(dr4**3*sin_th2)-cos_th2*fb2(5,2)/(dr4*sin_th2**2))
             ffb2(14,15)=(-sin_th2*fb2(5,3)*e4(2)+cos_th2*(-r4(2)*r4(3)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(2)-e3(2))*(-r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(5,3)/(dr4*sin_th2**2))
             ffb2(14,7)=(-sin_th2*fb2(3,1)*e4(2)+cos_th2*(r4(2)*r4(1)/dr4**3)- &
                  r3(2)*r3(1)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(2)-e3(2))*(r4(1)/(dr4**3*sin_th2)-cos_th2*fb2(3,1)/(dr4*sin_th2**2))
             ffb2(14,8)=(-sin_th2*fb2(3,2)*e4(2)+cos_th2*(-one/dr4+r4(2)*r4(2)/dr4**3)+ &
                  one/dr3-r3(2)*r3(2)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(2)-e3(2))*(r4(2)/(dr4**3*sin_th2)-cos_th2*fb2(3,2)/(dr4*sin_th2**2))
             ffb2(14,9)=(-sin_th2*fb2(3,3)*e4(2)+cos_th2*(r4(2)*r4(3)/dr4**3)- &
                  r3(2)*r3(3)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(2)-e3(2))*(r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr4*sin_th2**2))
             !.......................................
             ffb2(15,1:6)=zero
             ffb2(15,10:14)=ffb2(10:14,15)
             ffb2(15,15)=(-sin_th2*fb2(5,3)*e4(3)+cos_th2*(one/dr4-r4(3)*r4(3)/dr4**3))/(dr4*sin_th2)+ &
                  (cos_th2*e4(3)-e3(3))*(-r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(5,3)/(dr4*sin_th2**2))
             ffb2(15,7)=(-sin_th2*fb2(3,1)*e4(3)+cos_th2*(r4(3)*r4(1)/dr4**3)- &
                  r3(3)*r3(1)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(3)-e3(3))*(r4(1)/(dr4**3*sin_th2)-cos_th2*fb2(3,1)/(dr4*sin_th2**2))
             ffb2(15,8)=(-sin_th2*fb2(3,2)*e4(3)+cos_th2*(r4(3)*r4(2)/dr4**3)- &
                  r3(3)*r3(2)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(3)-e3(3))*(r4(2)/(dr4**3*sin_th2)-cos_th2*fb2(3,2)/(dr4*sin_th2**2))
             ffb2(15,9)=(-sin_th2*fb2(3,3)*e4(3)+cos_th2*(-one/dr4+r4(3)*r4(3)/dr4**3)+ &
                  one/dr3-r3(3)*r3(3)/dr3**3)/(dr4*sin_th2)+ &
                  (cos_th2*e4(3)-e3(3))*(r4(3)/(dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr4*sin_th2**2))
             !.......................................
             ffb2(7,1:6)=zero
             ffb2(7,7)=((-e3(1)+e4(1)*cos_th2+dr4*sin_th2*fb2(3,1))*e3(1)+ &
                      (dr3-dr4*cos_th2)*(-one/dr3+r3(1)*r3(1)/dr3**3)+ &
                      (-e4(1)+e3(1)*cos_th2+dr3*sin_th2*fb2(3,1))*e4(1)+ &
                      (dr4-dr3*cos_th2)*(-one/dr4+r4(1)*r4(1)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(1)+(dr4-dr3*cos_th2)*e4(1))* &
                      (r3(1)/(dr3**3*dr4*sin_th2)+r4(1)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,1)/(dr3*dr4*sin_th2**2))
             ffb2(7,8)=((-e3(2)+e4(2)*cos_th2+dr4*sin_th2*fb2(3,2))*e3(1)+ &
                      (dr3-dr4*cos_th2)*(r3(1)*r3(2)/dr3**3)+ &
                      (-e4(2)+e3(2)*cos_th2+dr3*sin_th2*fb2(3,2))*e4(1)+ &
                      (dr4-dr3*cos_th2)*(r4(1)*r4(2)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(1)+(dr4-dr3*cos_th2)*e4(1))* &
                      (r3(2)/(dr3**3*dr4*sin_th2)+r4(2)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,2)/(dr3*dr4*sin_th2**2))
             ffb2(7,9)=((-e3(3)+e4(3)*cos_th2+dr4*sin_th2*fb2(3,3))*e3(1)+ &
                      (dr3-dr4*cos_th2)*(r3(1)*r3(3)/dr3**3)+ &
                      (-e4(3)+e3(3)*cos_th2+dr3*sin_th2*fb2(3,3))*e4(1)+ &
                      (dr4-dr3*cos_th2)*(r4(1)*r4(3)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(1)+(dr4-dr3*cos_th2)*e4(1))* &
                      (r3(3)/(dr3**3*dr4*sin_th2)+r4(3)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*dr4*sin_th2**2))
             ffb2(7,10:15)=ffb2(10:15,7)
             !.......................................
             ffb2(8,1:6)=zero
             ffb2(8,7)=ffb2(7,8)
             ffb2(8,8)=((-e3(2)+e4(2)*cos_th2+dr4*sin_th2*fb2(3,2))*e3(2)+ &
                      (dr3-dr4*cos_th2)*(-one/dr3+r3(2)*r3(2)/dr3**3)+ &
                      (-e4(2)+e3(2)*cos_th2+dr3*sin_th2*fb2(3,2))*e4(2)+ &
                      (dr4-dr3*cos_th2)*(-one/dr4+r4(2)*r4(2)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(2)+(dr4-dr3*cos_th2)*e4(2))* &
                      (r3(2)/(dr3**3*dr4*sin_th2)+r4(2)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,2)/(dr3*dr4*sin_th2**2))
             ffb2(8,9)=((-e3(3)+e4(3)*cos_th2+dr4*sin_th2*fb2(3,3))*e3(2)+ &
                      (dr3-dr4*cos_th2)*(r3(2)*r3(3)/dr3**3)+ &
                      (-e4(3)+e3(3)*cos_th2+dr3*sin_th2*fb2(3,3))*e4(2)+ &
                      (dr4-dr3*cos_th2)*(r4(2)*r4(3)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(2)+(dr4-dr3*cos_th2)*e4(2))* &
                      (r3(3)/(dr3**3*dr4*sin_th2)+r4(3)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*dr4*sin_th2**2))
             ffb2(8,10:15)=ffb2(10:15,8)
             !.......................................
             ffb2(9,1:6)=zero
             ffb2(9,7:8)=ffb2(7:8,9)
             ffb2(9,9)=((-e3(3)+e4(3)*cos_th2+dr4*sin_th2*fb2(3,3))*e3(3)+ &
                      (dr3-dr4*cos_th2)*(-one/dr3+r3(3)*r3(3)/dr3**3)+ &
                      (-e4(3)+e3(3)*cos_th2+dr3*sin_th2*fb2(3,3))*e4(3)+ &
                      (dr4-dr3*cos_th2)*(-one/dr4+r4(3)*r4(3)/dr4**3))/(dr3*dr4*sin_th2)+ &
                      ((dr3-dr4*cos_th2)*e3(3)+(dr4-dr3*cos_th2)*e4(3))* &
                      (r3(3)/(dr3**3*dr4*sin_th2)+r4(3)/(dr3*dr4**3*sin_th2)-cos_th2*fb2(3,3)/(dr3*dr4*sin_th2**2))
             ffb2(9,10:15)=ffb2(10:15,9)
             !.......................................
             !.......................................
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             k4=3*(ia4-1)
             k5=3*(ia5-1)
             do l=1,3
                do m=1,15
                   if(m <= 3) then
                      fst_1=fb1(1,m)
                      fst_2=fb2(1,m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fst_1=fb1(2,m-3)
                      fst_2=fb2(2,m-3)
                   else if(m > 6 .and. m <= 9) then
                      m1=k3+(m-6)
                      fst_1=fb1(3,m-6)
                      fst_2=fb2(3,m-6)
                   else if(m > 9 .and. m <= 12) then
                      m1=k4+(m-9)
                      fst_1=fb1(4,m-9)
                      fst_2=fb2(4,m-9)
                   else if(m > 12) then
                      m1=k5+(m-12)
                      fst_1=fb1(5,m-12)
                      fst_2=fb2(5,m-12)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+(d2E_dt12*fst_1+d2E_dt1dt2*fst_2)*fb1(1,l)+ &
                        (d2E_dt2dt1*fst_1+d2E_dt22*fst_2)*fb2(1,l)+dE_dt1*ffb1(l,m)+dE_dt2*ffb2(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+(d2E_dt12*fst_1+d2E_dt1dt2*fst_2)*fb1(2,l)+ &
                        (d2E_dt2dt1*fst_1+d2E_dt22*fst_2)*fb2(2,l)+dE_dt1*ffb1(l+3,m)+dE_dt2*ffb2(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+(d2E_dt12*fst_1+d2E_dt1dt2*fst_2)*fb1(3,l)+ &
                        (d2E_dt2dt1*fst_1+d2E_dt22*fst_2)*fb2(3,l)+dE_dt1*ffb1(l+6,m)+dE_dt2*ffb2(l+6,m)
                   H(k4+l,m1)=H(k4+l,m1)+(d2E_dt12*fst_1+d2E_dt1dt2*fst_2)*fb1(4,l)+ &
                        (d2E_dt2dt1*fst_1+d2E_dt22*fst_2)*fb2(4,l)+dE_dt1*ffb1(l+9,m)+dE_dt2*ffb2(l+9,m)
                   H(k5+l,m1)=H(k5+l,m1)+(d2E_dt12*fst_1+d2E_dt1dt2*fst_2)*fb1(5,l)+ &
                        (d2E_dt2dt1*fst_1+d2E_dt22*fst_2)*fb2(5,l)+dE_dt1*ffb1(l+12,m)+dE_dt2*ffb2(l+12,m)
                end do
             end do
          end if

       end do j_lenght1
    end do i_n_bb
!!$print*,'22222222222'

    i_n_sb: do i=1,n_sb
       k=n_2b+n_3b+n_4b+n_ss+n_bb+i
       length=size(str_bnd(i)%list,2)

       j_lenght2: do j=1,length
          ia1=str_bnd(i)%list(1,j)
          ia2=str_bnd(i)%list(2,j)
          ia3=str_bnd(i)%list(3,j)
          r1=atoms_cart(ia1)%r-atoms_cart(ia2)%r
          r2=atoms_cart(ia3)%r-atoms_cart(ia2)%r
          if(lattice_calc) then 
             call image(r1)
             call image(r2)
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
          end if
          dr1=sqrt(dot_product(r1,r1))
          dr2=sqrt(dot_product(r2,r2))
          dot_prod1=dot_product(r1,r2)
          cos_th1=dot_prod1/(dr1*dr2)
          theta1=acos(cos_th1)
          
          id=poten(k)%id
          p=poten(k)%param
          select case(id)
          case (13) !str_bnd
             dr1a=dr1-p(2)
             dr2a=dr2-p(3)
             dtheta1=(theta1-p(4)*deg2rad)
             E_buf=p(1)*dr1a*dtheta1*rad2deg1+p(1)*dr2a*dtheta1*rad2deg1
             E(13)=E(13)+E_buf
             if(calc_gradients) then
                dE_dr1=p(1)*dtheta1*rad2deg1
                dE_dr2=p(1)*dtheta1*rad2deg1
                dE_dt1=p(1)*dr1a*rad2deg1+p(1)*dr2a*rad2deg1
             end if
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             sin_th1=sin(theta1)
             sin_th1 = sign(max(small,abs(sin_th1)),sin_th1)

             e1=r1/dr1
             e2=r2/dr2

             f_1(1,1:3)=r1/dr1
             f_1(2,1:3)=-r1/dr1
             f_1(3,1:3)=zero
             f_2(1,1:3)=zero
             f_2(2,1:3)=-r2/dr2
             f_2(3,1:3)=r2/dr2
             f_3(1,1:3)=(cos_th1*e1-e2)/(dr1*sin_th1)
             f_3(2,1:3)=((dr1-dr2*cos_th1)*e1+(dr2-dr1*cos_th1)*e2)/(dr1*dr2*sin_th1)
             f_3(3,1:3)=(cos_th1*e2-e1)/(dr2*sin_th1)

             Grad(:,ia1)=Grad(:,ia1)+dE_dr1*f_1(1,:)+dE_dr2*f_2(1,:)+dE_dt1*f_3(1,:)
             Grad(:,ia2)=Grad(:,ia2)+dE_dr1*f_1(2,:)+dE_dr2*f_2(2,:)+dE_dt1*f_3(2,:)
             Grad(:,ia3)=Grad(:,ia3)+dE_dr1*f_1(3,:)+dE_dr2*f_2(3,:)+dE_dt1*f_3(3,:)
          end if

          if(calc_hessian) then
             d2E_dr12=zero
             d2E_dr1dr2=zero
             d2E_dr1dt1=p(1)*rad2deg1
             d2E_dr2dr1=zero
             d2E_dr22=zero
             d2E_dr2dt1=p(1)*rad2deg1
             d2E_dt1dr1=p(1)*rad2deg1
             d2E_dt1dr2=p(1)*rad2deg1
             d2E_dt12=zero
             !...................................
             ff1(1,1)=one/dr1-r1(1)*r1(1)/dr1**3
             ff1(1,2)=-r1(1)*r1(2)/dr1**3
             ff1(1,3)=-r1(1)*r1(3)/dr1**3
             ff1(1,4:6)=-ff1(1,1:3)
             ff1(1,7:9)=zero
             !...................................
             ff1(2,1)=ff1(1,2)
             ff1(2,2)=one/dr1-r1(2)*r1(2)/dr1**3
             ff1(2,3)=-r1(2)*r1(3)/dr1**3
             ff1(2,4:6)=-ff1(2,1:3)
             ff1(2,7:9)=zero
             !...................................
             ff1(3,1:2)=ff1(1:2,3)
             ff1(3,3)=one/dr1-r1(3)*r1(3)/dr1**3
             ff1(3,4:6)=-ff1(3,1:3)
             ff1(3,7:9)=zero
             !...................................
             ff1(4,1:3)=ff1(1:3,4)
             ff1(4,4:6)=ff1(1,1:3)
             ff1(4,7:9)=zero
             !...................................
             ff1(5,1:3)=ff1(1:3,5)
             ff1(5,4:6)=ff1(2,1:3)
             ff1(5,7:9)=zero
             !...................................
             ff1(6,1:3)=ff1(1:3,6)
             ff1(6,4:6)=ff1(3,1:3)
             ff1(6,7:9)=zero
             !...................................
             ff1(7:9,1:9)=zero
             !...................................
             !...................................
             ff2(1:3,1:9)=zero
             !...................................
             ff2(4,1:3)=ff2(1:3,4)
             ff2(4,4)=one/dr2-r2(1)*r2(1)/dr2**3
             ff2(4,5)=-r2(1)*r2(2)/dr2**3
             ff2(4,6)=-r2(1)*r2(3)/dr2**3
             ff2(4,7:9)=-ff2(4,4:6)
             !...................................
             ff2(5,1:3)=ff2(1:3,5)
             ff2(5,4)=ff2(4,5)
             ff2(5,5)=one/dr2-r2(2)*r2(2)/dr2**3
             ff2(5,6)=-r2(2)*r2(3)/dr2**3
             ff2(5,7:9)=-ff2(5,4:6)
             !...................................
             ff2(6,1:3)=ff2(1:3,6)
             ff2(6,4:5)=ff2(4:5,6)
             ff2(6,6)=one/dr2-r2(3)*r2(3)/dr2**3
             ff2(6,7:9)=-ff2(6,4:6)
             !...................................
             ff2(7,1:6)=ff2(1:6,7)
             ff2(7,7:9)=ff2(4,4:6)
             !...................................
             ff2(8,1:6)=ff2(1:6,8)
             ff2(8,7:9)=ff2(5,4:6)
             !...................................
             ff2(9,1:6)=ff2(1:6,9)
             ff2(9,7:9)=ff2(6,4:6)
             !...................................
             !...................................
             ff3(1,1)=(-sin_th1*f_3(1,1)*e1(1)+cos_th1*(one/dr1-r1(1)*r1(1)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(1)/(dr1**3*sin_th1)-cos_th1*f_3(1,1)/(dr1*sin_th1**2))
             ff3(1,2)=(-sin_th1*f_3(1,2)*e1(1)+cos_th1*(-r1(1)*r1(2)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(2)/(dr1**3*sin_th1)-cos_th1*f_3(1,2)/(dr1*sin_th1**2))
             ff3(1,3)=(-sin_th1*f_3(1,3)*e1(1)+cos_th1*(-r1(1)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(1,3)/(dr1*sin_th1**2))
             ff3(1,4)=(-sin_th1*f_3(2,1)*e1(1)+cos_th1*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                  one/dr2-r2(1)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(1)/(dr1**3*sin_th1)-cos_th1*f_3(2,1)/(dr1*sin_th1**2))
             ff3(1,5)=(-sin_th1*f_3(2,2)*e1(1)+cos_th1*(r1(1)*r1(2)/dr1**3)- &
                  r2(1)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(2)/(dr1**3*sin_th1)-cos_th1*f_3(2,2)/(dr1*sin_th1**2))
             ff3(1,6)=(-sin_th1*f_3(2,3)*e1(1)+cos_th1*(r1(1)*r1(3)/dr1**3)- &
                  r2(1)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*sin_th1**2))
             ff3(1,7)=(-sin_th1*f_3(3,1)*e1(1)-one/dr2+r2(1)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*f_3(3,1)/(dr1*sin_th1**2))
             ff3(1,8)=(-sin_th1*f_3(3,2)*e1(1)+r2(1)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*f_3(3,2)/(dr1*sin_th1**2))
             ff3(1,9)=(-sin_th1*f_3(3,3)*e1(1)+r2(1)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(1)-e2(1))*(-cos_th1*f_3(3,3)/(dr1*sin_th1**2))
             !........................................................
             ff3(2,1)=ff3(1,2)
             ff3(2,2)=(-sin_th1*f_3(1,2)*e1(2)+cos_th1*(one/dr1-r1(2)*r1(2)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-r1(2)/(dr1**3*sin_th1)-cos_th1*f_3(1,2)/(dr1*sin_th1**2))
             ff3(2,3)=(-sin_th1*f_3(1,3)*e1(2)+cos_th1*(-r1(2)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(1,3)/(dr1*sin_th1**2))
             ff3(2,4)=(-sin_th1*f_3(2,1)*e1(2)+cos_th1*(r1(2)*r1(1)/dr1**3)- &
                  r2(2)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(1)/(dr1**3*sin_th1)-cos_th1*f_3(2,1)/(dr1*sin_th1**2))
             ff3(2,5)=(-sin_th1*f_3(2,2)*e1(2)+cos_th1*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                  one/dr2-r2(2)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(2)/(dr1**3*sin_th1)-cos_th1*f_3(2,2)/(dr1*sin_th1**2))
             ff3(2,6)=(-sin_th1*f_3(2,3)*e1(2)+cos_th1*(r1(2)*r1(3)/dr1**3)- &
                  r2(2)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*sin_th1**2))
             ff3(2,7)=(-sin_th1*f_3(3,1)*e1(2)+r2(2)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*f_3(3,1)/(dr1*sin_th1**2))
             ff3(2,8)=(-sin_th1*f_3(3,2)*e1(2)-one/dr2+r2(2)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*f_3(3,2)/(dr1*sin_th1**2))
             ff3(2,9)=(-sin_th1*f_3(3,3)*e1(2)+r2(2)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(2)-e2(2))*(-cos_th1*f_3(3,3)/(dr1*sin_th1**2))
             !........................................................
             ff3(3,1:2)=ff3(1:2,3)
             ff3(3,3)=(-sin_th1*f_3(1,3)*e1(3)+cos_th1*(one/dr1-r1(3)*r1(3)/dr1**3))/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(1,3)/(dr1*sin_th1**2))
             ff3(3,4)=(-sin_th1*f_3(2,1)*e1(3)+cos_th1*(r1(3)*r1(1)/dr1**3)- &
                  r2(3)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(1)/(dr1**3*sin_th1)-cos_th1*f_3(2,1)/(dr1*sin_th1**2))
             ff3(3,5)=(-sin_th1*f_3(2,2)*e1(3)+cos_th1*(r1(3)*r1(2)/dr1**3)- &
                  r2(3)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(2)/(dr1**3*sin_th1)-cos_th1*f_3(2,2)/(dr1*sin_th1**2))
             ff3(3,6)=(-sin_th1*f_3(2,3)*e1(3)+cos_th1*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                  one/dr2-r2(3)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(r1(3)/(dr1**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*sin_th1**2))
             ff3(3,7)=(-sin_th1*f_3(3,1)*e1(3)+r2(3)*r2(1)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*f_3(3,1)/(dr1*sin_th1**2))
             ff3(3,8)=(-sin_th1*f_3(3,2)*e1(3)+r2(3)*r2(2)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*f_3(3,2)/(dr1*sin_th1**2))
             ff3(3,9)=(-sin_th1*f_3(3,3)*e1(3)-one/dr2+r2(3)*r2(3)/dr2**3)/(dr1*sin_th1)+ &
                  (cos_th1*e1(3)-e2(3))*(-cos_th1*f_3(3,3)/(dr1*sin_th1**2))
             !........................................................
             ff3(7,1:3)=ff3(1:3,7)
             ff3(7,4)=(-sin_th1*f_3(2,1)*e2(1)+cos_th1*(-one/dr2+r2(1)*r2(1)/dr2**3)+ &
                  one/dr1-r1(1)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(1)/(dr2**3*sin_th1)-cos_th1*f_3(2,1)/(dr2*sin_th1**2))
             ff3(7,5)=(-sin_th1*f_3(2,2)*e2(1)+cos_th1*(r2(1)*r2(2)/dr2**3)- &
                  r1(1)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(2)/(dr2**3*sin_th1)-cos_th1*f_3(2,2)/(dr2*sin_th1**2))
             ff3(7,6)=(-sin_th1*f_3(2,3)*e2(1)+cos_th1*(r2(1)*r2(3)/dr2**3)- &
                  r1(1)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr2*sin_th1**2))
             ff3(7,7)=(-sin_th1*f_3(3,1)*e2(1)+cos_th1*(one/dr2-r2(1)*r2(1)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(1)/(dr2**3*sin_th1)-cos_th1*f_3(3,1)/(dr2*sin_th1**2))
             ff3(7,8)=(-sin_th1*f_3(3,2)*e2(1)+cos_th1*(-r2(1)*r2(2)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(2)/(dr2**3*sin_th1)-cos_th1*f_3(3,2)/(dr2*sin_th1**2))
             ff3(7,9)=(-sin_th1*f_3(3,3)*e2(1)+cos_th1*(-r2(1)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(1)-e1(1))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(3,3)/(dr2*sin_th1**2))
             !........................................................
             ff3(8,1:3)=ff3(1:3,8)
             ff3(8,4)=(-sin_th1*f_3(2,1)*e2(2)+cos_th1*(r2(2)*r2(1)/dr2**3)- &
                  r1(2)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(1)/(dr2**3*sin_th1)-cos_th1*f_3(2,1)/(dr2*sin_th1**2))
             ff3(8,5)=(-sin_th1*f_3(2,2)*e2(2)+cos_th1*(-one/dr2+r2(2)*r2(2)/dr2**3)+ &
                  one/dr1-r1(2)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(2)/(dr2**3*sin_th1)-cos_th1*f_3(2,2)/(dr2*sin_th1**2))
             ff3(8,6)=(-sin_th1*f_3(2,3)*e2(2)+cos_th1*(r2(2)*r2(3)/dr2**3)- &
                  r1(2)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr2*sin_th1**2))
             ff3(8,7)=ff3(7,8)
             ff3(8,8)=(-sin_th1*f_3(3,2)*e2(2)+cos_th1*(one/dr2-r2(2)*r2(2)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(-r2(2)/(dr2**3*sin_th1)-cos_th1*f_3(3,2)/(dr2*sin_th1**2))
             ff3(8,9)=(-sin_th1*f_3(3,3)*e2(2)+cos_th1*(-r2(2)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(2)-e1(2))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(3,3)/(dr2*sin_th1**2))
             !........................................................
             ff3(9,1:3)=ff3(1:3,9)
             ff3(9,4)=(-sin_th1*f_3(2,1)*e2(3)+cos_th1*(r2(3)*r2(1)/dr2**3)- &
                  r1(3)*r1(1)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(1)/(dr2**3*sin_th1)-cos_th1*f_3(2,1)/(dr2*sin_th1**2))
             ff3(9,5)=(-sin_th1*f_3(2,2)*e2(3)+cos_th1*(r2(3)*r2(2)/dr2**3)- &
                  r1(3)*r1(2)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(2)/(dr2**3*sin_th1)-cos_th1*f_3(2,2)/(dr2*sin_th1**2))
             ff3(9,6)=(-sin_th1*f_3(2,3)*e2(3)+cos_th1*(-one/dr2+r2(3)*r2(3)/dr2**3)+ &
                  one/dr1-r1(3)*r1(3)/dr1**3)/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr2*sin_th1**2))
             ff3(9,7:8)=ff3(7:8,9)
             ff3(9,9)=(-sin_th1*f_3(3,3)*e2(3)+cos_th1*(one/dr2-r2(3)*r2(3)/dr2**3))/(dr2*sin_th1)+ &
                  (cos_th1*e2(3)-e1(3))*(-r2(3)/(dr2**3*sin_th1)-cos_th1*f_3(3,3)/(dr2*sin_th1**2))
             !........................................................
             ff3(4,1:3)=ff3(1:3,4)
             ff3(4,4)=((-e1(1)+e2(1)*cos_th1+dr2*sin_th1*f_3(2,1))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(1)*r1(1)/dr1**3)+ &
                      (-e2(1)+e1(1)*cos_th1+dr1*sin_th1*f_3(2,1))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(1)*r2(1)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(1)/(dr1**3*dr2*sin_th1)+r2(1)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,1)/(dr1*dr2*sin_th1**2))
             ff3(4,5)=((-e1(2)+e2(2)*cos_th1+dr2*sin_th1*f_3(2,2))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(r1(1)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th1+dr1*sin_th1*f_3(2,2))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(r2(1)*r2(2)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(2)/(dr1**3*dr2*sin_th1)+r2(2)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,2)/(dr1*dr2*sin_th1**2))
             ff3(4,6)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*f_3(2,3))*e1(1)+ &
                      (dr1-dr2*cos_th1)*(r1(1)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*f_3(2,3))*e2(1)+ &
                      (dr2-dr1*cos_th1)*(r2(1)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(1)+(dr2-dr1*cos_th1)*e2(1))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*dr2*sin_th1**2))
             ff3(4,7:9)=ff3(7:9,4)
             !........................................................
             ff3(5,1:4)=ff3(1:4,5)
             ff3(5,5)=((-e1(2)+e2(2)*cos_th1+dr2*sin_th1*f_3(2,2))*e1(2)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(2)*r1(2)/dr1**3)+ &
                      (-e2(2)+e1(2)*cos_th1+dr1*sin_th1*f_3(2,2))*e2(2)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(2)*r2(2)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(2)+(dr2-dr1*cos_th1)*e2(2))* &
                      (r1(2)/(dr1**3*dr2*sin_th1)+r2(2)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,2)/(dr1*dr2*sin_th1**2))
             ff3(5,6)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*f_3(2,3))*e1(2)+ &
                      (dr1-dr2*cos_th1)*(r1(2)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*f_3(2,3))*e2(2)+ &
                      (dr2-dr1*cos_th1)*(r2(2)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(2)+(dr2-dr1*cos_th1)*e2(2))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*dr2*sin_th1**2))
             ff3(5,7:9)=ff3(7:9,5)
             !........................................................
             ff3(6,1:5)=ff3(1:5,6)
             ff3(6,6)=((-e1(3)+e2(3)*cos_th1+dr2*sin_th1*f_3(2,3))*e1(3)+ &
                      (dr1-dr2*cos_th1)*(-one/dr1+r1(3)*r1(3)/dr1**3)+ &
                      (-e2(3)+e1(3)*cos_th1+dr1*sin_th1*f_3(2,3))*e2(3)+ &
                      (dr2-dr1*cos_th1)*(-one/dr2+r2(3)*r2(3)/dr2**3))/(dr1*dr2*sin_th1)+ &
                      ((dr1-dr2*cos_th1)*e1(3)+(dr2-dr1*cos_th1)*e2(3))* &
                      (r1(3)/(dr1**3*dr2*sin_th1)+r2(3)/(dr1*dr2**3*sin_th1)-cos_th1*f_3(2,3)/(dr1*dr2*sin_th1**2))
             ff3(6,7:9)=ff3(7:9,6)
             !........................................................
             !........................................................
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             do l=1,3
                do m=1,9
                   if(m <= 3) then
                      fst_1=f_1(1,m)
                      fst_2=f_2(1,m)
                      fst_3=f_3(1,m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fst_1=f_1(2,m-3)
                      fst_2=f_2(2,m-3)
                      fst_3=f_3(2,m-3)
                   else if(m > 6) then
                      m1=k3+(m-6)
                      fst_1=f_1(3,m-6)
                      fst_2=f_2(3,m-6)
                      fst_3=f_3(3,m-6)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2+d2E_dr1dt1*fst_3)*f_1(1,l)+ &
                                         (d2E_dr2dr1*fst_1+d2E_dr22*fst_2+d2E_dr2dt1*fst_3)*f_2(1,l)+ &
                                         (d2E_dt1dr1*fst_1+d2E_dt1dr2*fst_2+d2E_dt12*fst_3)*f_3(1,l)+ &
                                         dE_dr1*ff1(l,m)+dE_dr2*ff2(l,m)+dE_dt1*ff3(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2+d2E_dr1dt1*fst_3)*f_1(2,l)+ &
                                         (d2E_dr2dr1*fst_1+d2E_dr22*fst_2+d2E_dr2dt1*fst_3)*f_2(2,l)+ &
                                         (d2E_dt1dr1*fst_1+d2E_dt1dr2*fst_2+d2E_dt12*fst_3)*f_3(2,l)+ &
                                         dE_dr1*ff1(l+3,m)+dE_dr2*ff2(l+3,m)+dE_dt1*ff3(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+(d2E_dr12*fst_1+d2E_dr1dr2*fst_2+d2E_dr1dt1*fst_3)*f_1(3,l)+ &
                                         (d2E_dr2dr1*fst_1+d2E_dr22*fst_2+d2E_dr2dt1*fst_3)*f_2(3,l)+ &
                                         (d2E_dt1dr1*fst_1+d2E_dt1dr2*fst_2+d2E_dt12*fst_3)*f_3(3,l)+ &
                                         dE_dr1*ff1(l+6,m)+dE_dr2*ff2(l+6,m)+dE_dt1*ff3(l+6,m)
                end do
             end do
             
          end if

       end do j_lenght2
    end do i_n_sb
!!$print*,'3333333333333'

    i_n_st: do i=1,n_st
       k=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+i
       length=size(str_trs(i)%list,2)

       j_lenght3: do j=1,length
          ia1=str_trs(i)%list(1,j)
          ia2=str_trs(i)%list(2,j)
          ia3=str_trs(i)%list(3,j)
          ia4=str_trs(i)%list(4,j)
          r1=atoms_cart(ia2)%r-atoms_cart(ia1)%r
          r2=atoms_cart(ia3)%r-atoms_cart(ia2)%r
          r3=atoms_cart(ia4)%r-atoms_cart(ia3)%r
          if(lattice_calc) then 
             call image(r1)
             call image(r2)
             call image(r3)
          else if(slab_calc) then
             call image_slab(r1)
             call image_slab(r2)
             call image_slab(r3)
          end if
          r21=vector_product(r1,r2)
          r32=vector_product(r2,r3)
          dr2=sqrt(dot_product(r2,r2))
          dr21=sqrt(dot_product(r21,r21))
          dr32=sqrt(dot_product(r32,r32))
          dot_prod1=dot_product(r21,r32)
          cos_phi=dot_prod1/(dr21*dr32)
          if(cos_phi < mone) cos_phi=mone
          if(cos_phi > one) cos_phi=one
          phi=acos(cos_phi)

          id=poten(k)%id
          p=poten(k)%param
          select case (id)
          case (14) ! str_trs
             dr2a=dr2-p(2)
             E_buf=(p(1)/two)*dr2a*(one+cos(three*phi))
             E(14)=E(14)+E_buf
             if(calc_gradients) then
                dE_dr2=(p(1)/two)*(one+cos(three*phi))
                dE_dp=-(three*p(1)/two)*dr2a*sin(three*phi)
             end if
          end select
          E_total=E_total+E_buf

          if(calc_gradients) then
             sin_phi=sin(phi)
             sin_phi = sign(max(small,abs(sin_phi)),sin_phi)

             A(1,1)=-r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))+r3(1)*(r2(2)*r2(2)+r2(3)*r2(3))
             A(1,2)=-r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))+r3(2)*(r2(1)*r2(1)+r2(3)*r2(3))
             A(1,3)=-r2(3)*(r2(1)*r3(1)+r2(2)*r3(2))+r3(3)*(r2(1)*r2(1)+r2(2)*r2(2))

             A(2,1)=-r1(1)*(r2(2)*r3(2)+r2(3)*r3(3))+r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))- &
                  r3(1)*(r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))+ &
                  two*r2(1)*(r1(2)*r3(2)+r1(3)*r3(3))
             A(2,2)=-r1(2)*(r2(1)*r3(1)+r2(3)*r3(3))+r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))- &
                  r3(2)*(r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))+ &
                  two*r2(2)*(r1(1)*r3(1)+r1(3)*r3(3))
             A(2,3)=-r1(3)*(r2(1)*r3(1)+r2(2)*r3(2))+r2(3)*(r2(1)*r3(1)+r2(2)*r3(2))- &
                  r3(3)*(r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))+ &
                  two*r2(3)*(r1(1)*r3(1)+r1(2)*r3(2))

             A(3,1)=r1(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))- &
                  r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))+r3(1)*(r1(2)*r2(2)+r1(3)*r2(3))- &
                  two*r2(1)*(r1(2)*r3(2)+r1(3)*r3(3))
             A(3,2)=r1(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))- &
                  r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))+r3(2)*(r1(1)*r2(1)+r1(3)*r2(3))- &
                  two*r2(2)*(r1(1)*r3(1)+r1(3)*r3(3))
             A(3,3)=r1(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r2(2)*r3(2)+r2(1)*r3(1))- &
                  r2(3)*(r1(2)*r2(2)+r1(1)*r2(1))+r3(3)*(r1(2)*r2(2)+r1(1)*r2(1))- &
                  two*r2(3)*(r1(2)*r3(2)+r1(1)*r3(1))

             A(4,1)=-r1(1)*(r2(2)*r2(2)+r2(3)*r2(3))+r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))
             A(4,2)=-r1(2)*(r2(1)*r2(1)+r2(3)*r2(3))+r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))
             A(4,3)=-r1(3)*(r2(1)*r2(1)+r2(2)*r2(2))+r2(3)*(r1(1)*r2(1)+r1(2)*r2(2))

             B(1,1)=-two*r1(1)*(r2(2)*r2(2)+r2(3)*r2(3))+two*r2(1)*(r1(2)*r2(2)+r1(3)*r2(3))
             B(1,2)=-two*r1(2)*(r2(1)*r2(1)+r2(3)*r2(3))+two*r2(2)*(r1(1)*r2(1)+r1(3)*r2(3))
             B(1,3)=-two*r1(3)*(r2(2)*r2(2)+r2(1)*r2(1))+two*r2(3)*(r1(2)*r2(2)+r1(1)*r2(1))

             B(2,1)=two*r1(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r1(2)*r2(2)+r1(3)*r2(3))- &
                  two*r2(1)*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             B(2,2)=two*r1(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r1(1)*r2(1)+r1(3)*r2(3))- &
                  two*r2(2)*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             B(2,3)=two*r1(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r1(2)*r2(2)+r1(1)*r2(1))- &
                  two*r2(3)*(r1(2)*r1(2)+r1(1)*r1(1)+r1(2)*r2(2)+r1(1)*r2(1))

             B(3,1)=-two*r1(1)*(r1(2)*r2(2)+r1(3)*r2(3))+two*r2(1)*(r1(2)*r1(2)+r1(3)*r1(3))
             B(3,2)=-two*r1(2)*(r1(1)*r2(1)+r1(3)*r2(3))+two*r2(2)*(r1(1)*r1(1)+r1(3)*r1(3))
             B(3,3)=-two*r1(3)*(r1(2)*r2(2)+r1(1)*r2(1))+two*r2(3)*(r1(2)*r1(2)+r1(1)*r1(1))

             B(4,:)=zero

             C(1,:)=zero

             C(2,1)=two*r3(1)*(r2(2)*r3(2)+r2(3)*r3(3))-two*r2(1)*(r3(2)*r3(2)+r3(3)*r3(3))
             C(2,2)=two*r3(2)*(r2(1)*r3(1)+r2(3)*r3(3))-two*r2(2)*(r3(1)*r3(1)+r3(3)*r3(3))
             C(2,3)=two*r3(3)*(r2(2)*r3(2)+r2(1)*r3(1))-two*r2(3)*(r3(2)*r3(2)+r3(1)*r3(1))

             C(3,1)=-two*r3(1)*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))+ &
                  two*r2(1)*(r3(2)*r3(2)+r3(3)*r3(3)+r2(2)*r3(2)+r2(3)*r3(3))
             C(3,2)=-two*r3(2)*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))+ &
                  two*r2(2)*(r3(1)*r3(1)+r3(3)*r3(3)+r2(1)*r3(1)+r2(3)*r3(3))
             C(3,3)=-two*r3(3)*(r2(2)*r2(2)+r2(1)*r2(1)+r2(2)*r3(2)+r2(1)*r3(1))+ &
                  two*r2(3)*(r3(2)*r3(2)+r3(1)*r3(1)+r2(2)*r3(2)+r2(1)*r3(1))

             C(4,1)=two*r3(1)*(r2(2)*r2(2)+r2(3)*r2(3))-two*r2(1)*(r2(2)*r3(2)+r2(3)*r3(3))
             C(4,2)=two*r3(2)*(r2(1)*r2(1)+r2(3)*r2(3))-two*r2(2)*(r2(1)*r3(1)+r2(3)*r3(3))
             C(4,3)=two*r3(3)*(r2(2)*r2(2)+r2(1)*r2(1))-two*r2(3)*(r2(2)*r3(2)+r2(1)*r3(1))

             ft1(1,1:3)=zero
             ft1(2,1:3)=-r2/dr2
             ft1(3,1:3)=r2/dr2
             ft1(4,1:3)=zero
             ft2(1,1:3)=-(A(1,:)/(dr21*dr32)-(cos_phi/two)*(B(1,:)/(dr21*dr21)+C(1,:)/(dr32*dr32)))/sin_phi
             ft2(2,1:3)=-(A(2,:)/(dr21*dr32)-(cos_phi/two)*(B(2,:)/(dr21*dr21)+C(2,:)/(dr32*dr32)))/sin_phi
             ft2(3,1:3)=-(A(3,:)/(dr21*dr32)-(cos_phi/two)*(B(3,:)/(dr21*dr21)+C(3,:)/(dr32*dr32)))/sin_phi
             ft2(4,1:3)=-(A(4,:)/(dr21*dr32)-(cos_phi/two)*(B(4,:)/(dr21*dr21)+C(4,:)/(dr32*dr32)))/sin_phi

             Grad(:,ia1)=Grad(:,ia1)+dE_dr2*ft1(1,:)+dE_dp*ft2(1,:)
             Grad(:,ia2)=Grad(:,ia2)+dE_dr2*ft1(2,:)+dE_dp*ft2(2,:)
             Grad(:,ia3)=Grad(:,ia3)+dE_dr2*ft1(3,:)+dE_dp*ft2(3,:)
             Grad(:,ia4)=Grad(:,ia4)+dE_dr2*ft1(4,:)+dE_dp*ft2(4,:)
          end if
          
          if(calc_hessian) then
             d2E_dr22=zero
             d2E_dr2dp=-(three*p(1)/two)*sin(three*phi)
             d2E_dpdr2=-(three*p(1)/two)*sin(three*phi)
             d2E_dp2=-(nine*p(1)/two)*dr2a*cos(three*phi)
             !......................................
             !......................................
             fft1(1:3,1:12)=zero
             !......................................
             fft1(4,1:3)=zero
             fft1(4,4)=one/dr2-r2(1)*r2(1)/dr2**3
             fft1(4,5)=-r2(1)*r2(2)/dr2**3
             fft1(4,6)=-r2(1)*r2(3)/dr2**3
             fft1(4,7:9)=-fft1(4,4:6)
             fft1(4,10:12)=zero
             !......................................
             fft1(5,1:3)=zero
             fft1(5,4)=fft1(4,5)
             fft1(5,5)=one/dr2-r2(2)*r2(2)/dr2**3
             fft1(5,6)=-r2(2)*r2(3)/dr2**3
             fft1(5,7:9)=-fft1(5,4:6)
             fft1(4,10:12)=zero
             !......................................
             fft1(6,1:3)=zero
             fft1(6,4:5)=fft1(4:5,6)
             fft1(6,6)=one/dr2-r2(3)*r2(3)/dr2**3
             fft1(6,7:9)=-fft1(6,4:6)
             fft1(6,10:12)=zero
             !......................................
             fft1(7,1:3)=zero
             fft1(7,4:6)=fft1(4:6,7)
             fft1(7,7:9)=fft1(4,4:6)
             fft1(7,10:12)=zero
             !......................................
             fft1(8,1:3)=zero
             fft1(8,4:6)=fft1(4:6,8)
             fft1(8,7:9)=fft1(5,4:6)
             fft1(8,10:12)=zero
             !......................................
             fft1(9,1:3)=zero
             fft1(9,4:6)=fft1(4:6,9)
             fft1(9,7:9)=fft1(6,4:6)
             fft1(9,10:12)=zero
             !......................................
             fft1(10:12,1:12)=zero
             !......................................
             !......................................
             dcos(1:3)=-sin_phi*ft2(1,:)
             dcos(4:6)=-sin_phi*ft2(2,:)
             dcos(7:9)=-sin_phi*ft2(3,:)
             dcos(10:12)=-sin_phi*ft2(4,:)

             d_sin(1:3)=-ft2(1,:)*cos_phi/sin_phi**2
             d_sin(4:6)=-ft2(2,:)*cos_phi/sin_phi**2
             d_sin(7:9)=-ft2(3,:)*cos_phi/sin_phi**2
             d_sin(10:12)=-ft2(4,:)*cos_phi/sin_phi**2

             m2121=one/(dr21*dr21)
             m2132=one/(dr21*dr32)
             m3232=one/(dr32*dr32)

             g2121(1:3)=-B(1,1:3)/(dr21**4)
             g2121(4:6)=-B(2,1:3)/(dr21**4)
             g2121(7:9)=-B(3,1:3)/(dr21**4)
             g2121(10:12)=-B(4,1:3)/(dr21**4)

             g3232(1:3)=-C(1,1:3)/(dr32**4)
             g3232(4:6)=-C(2,1:3)/(dr32**4)
             g3232(7:9)=-C(3,1:3)/(dr32**4)
             g3232(10:12)=-C(4,1:3)/(dr32**4)

             g2132(1:3)=-B(1,1:3)/(two*dr21**3*dr32)-C(1,1:3)/(two*dr21*dr32**3)
             g2132(4:6)=-B(2,1:3)/(two*dr21**3*dr32)-C(2,1:3)/(two*dr21*dr32**3)
             g2132(7:9)=-B(3,1:3)/(two*dr21**3*dr32)-C(3,1:3)/(two*dr21*dr32**3)
             g2132(10:12)=-B(4,1:3)/(two*dr21**3*dr32)-C(4,1:3)/(two*dr21*dr32**3)
             !......................................
             aa(1,1:3)=zero
             aa(1,4)=r2(2)*r3(2)+r2(3)*r3(3)
             aa(1,5)=-r2(1)*(-r3(2))+r3(1)*(-two*r2(2))
             aa(1,6)=-r2(1)*(-r3(3))+r3(1)*(-two*r2(3))
             aa(1,7)=-(r2(2)*r3(2)+r2(3)*r3(3))-(r2(2)*r2(2)+r2(3)*r2(3))
             aa(1,8)=-r2(1)*(r3(2)-r2(2))+r3(1)*(two*r2(2))
             aa(1,9)=-r2(1)*(r3(3)-r2(3))+r3(1)*(two*r2(3))
             aa(1,10)=r2(2)*r2(2)+r2(3)*r2(3)
             aa(1,11)=-r2(1)*(r2(2))
             aa(1,12)=-r2(1)*(r2(3))
             !......................................
             aa(2,1:3)=zero
             aa(2,4)=-r2(2)*(-r3(1))+r3(2)*(-two*r2(1))
             aa(2,5)=r2(1)*r3(1)+r2(3)*r3(3)
             aa(2,6)=-r2(2)*(-r3(3))+r3(2)*(-two*r2(3))
             aa(2,7)=-r2(2)*(r3(1)-r2(1))+r3(2)*(two*r2(1))
             aa(2,8)=-(r2(1)*r3(1)+r2(3)*r3(3))-(r2(1)*r2(1)+r2(3)*r2(3))
             aa(2,9)=-r2(2)*(r3(3)-r2(3))+r3(2)*(two*r2(3))
             aa(2,10)=-r2(2)*(r2(1))
             aa(2,11)=r2(1)*r2(1)+r2(3)*r2(3)
             aa(2,12)=-r2(2)*(r2(3))
             !......................................
             aa(3,1:3)=zero
             aa(3,4)=-r2(3)*(-r3(1))+r3(3)*(-two*r2(1))
             aa(3,5)=-r2(3)*(-r3(2))+r3(3)*(-two*r2(2))
             aa(3,6)=r2(1)*r3(1)+r2(2)*r3(2)
             aa(3,7)=-r2(3)*(r3(1)-r2(1))+r3(3)*(two*r2(1))
             aa(3,8)=-r2(3)*(r3(2)-r2(2))+r3(3)*(two*r2(2))
             aa(3,9)=-(r2(1)*r3(1)+r2(2)*r3(2))-(r2(1)*r2(1)+r2(2)*r2(2))
             aa(3,10)=-r2(3)*(r2(1))
             aa(3,11)=-r2(3)*(r2(2))
             aa(3,12)=r2(1)*r2(1)+r2(2)*r2(2)
             !......................................
             aa(4,1:3)=aa(1:3,4)
             aa(4,4)=-two*(r2(2)*r3(2)+r2(3)*r3(3))-two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(4,5)=-r1(1)*(-r3(2))+r2(1)*(-r3(2))-r3(1)*(r2(2)-r1(2)-two*r2(2))+ &
                     two*r2(1)*(r3(2))
             aa(4,6)=-r1(1)*(-r3(3))+r2(1)*(-r3(3))-r3(1)*(r2(3)-r1(3)-two*r2(3))+ &
                     two*r2(1)*(r3(3))
             aa(4,7)=(r2(2)*r3(2)+r2(3)*r3(3))+ &
                     (r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))+two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(4,8)=-r1(1)*(r3(2)-r2(2))+r2(1)*(r3(2)-r2(2))-r3(1)*(r1(2)+two*r2(2))+ &
                     two*r2(1)*(-r1(2))
             aa(4,9)=-r1(1)*(r3(3)-r2(3))+r2(1)*(r3(3)-r2(3))-r3(1)*(r1(3)+two*r2(3))+ &
                     two*r2(1)*(-r1(3))
             aa(4,10)=-(r1(2)*r2(2)+r1(3)*r2(3)+r2(2)*r2(2)+r2(3)*r2(3))
             aa(4,11)=-r1(1)*(r2(2))+r2(1)*(r2(2))+two*r2(1)*(r1(2))
             aa(4,12)=-r1(1)*(r2(3))+r2(1)*(r2(3))+two*r2(1)*(r1(3))
             !......................................
             aa(5,1:3)=aa(1:3,5)
             aa(5,4)=-r1(2)*(-r3(1))+r2(2)*(-r3(1))-r3(2)*(r2(1)-r1(1)-two*r2(1))+ &
                     two*r2(2)*(r3(1))             
             aa(5,5)=-two*(r2(1)*r3(1)+r2(3)*r3(3))-two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(5,6)=-r1(2)*(-r3(3))+r2(2)*(-r3(3))-r3(2)*(r2(3)-r1(3)-two*r2(3))+ &
                     two*r2(2)*(r3(3))             
             aa(5,7)=-r1(2)*(r3(1)-r2(1))+r2(2)*(r3(1)-r2(1))-r3(2)*(r1(1)+two*r2(1))+ &
                     two*r2(2)*(-r1(1))
             aa(5,8)=(r2(1)*r3(1)+r2(3)*r3(3))+ &
                     (r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))+two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(5,9)=-r1(2)*(r3(3)-r2(3))+r2(2)*(r3(3)-r2(3))-r3(2)*(r1(3)+two*r2(3))+ &
                     two*r2(2)*(-r1(3))
             aa(5,10)=-r1(2)*(r2(1))+r2(2)*(r2(1))+two*r2(2)*(r1(1))
             aa(5,11)=-(r1(1)*r2(1)+r1(3)*r2(3)+r2(1)*r2(1)+r2(3)*r2(3))
             aa(5,12)=-r1(2)*(r2(3))+r2(2)*(r2(3))+two*r2(2)*(r1(3))
             !......................................
             aa(6,1:3)=aa(1:3,6)
             aa(6,4)=-r1(3)*(-r3(1))+r2(3)*(-r3(1))-r3(3)*(r2(1)-r1(1)-two*r2(1))+ &
                     two*r2(3)*(r3(1))             
             aa(6,5)=-r1(3)*(-r3(2))+r2(3)*(-r3(2))-r3(3)*(r2(2)-r1(2)-two*r2(2))+ &
                     two*r2(3)*(r3(2))             
             aa(6,6)=-two*(r2(1)*r3(1)+r2(2)*r3(2))-two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(6,7)=-r1(3)*(r3(1)-r2(1))+r2(3)*(r3(1)-r2(1))-r3(3)*(r1(1)+two*r2(1))+ &
                     two*r2(3)*(-r1(1))
             aa(6,8)=-r1(3)*(r3(2)-r2(2))+r2(3)*(r3(2)-r2(2))-r3(3)*(r1(2)+two*r2(2))+ &
                     two*r2(3)*(-r1(2))
             aa(6,9)=(r2(1)*r3(1)+r2(2)*r3(2))+ &
                     (r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))+two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(6,10)=-r1(3)*(r2(1))+r2(3)*(r2(1))+two*r2(3)*(r1(1))
             aa(6,11)=-r1(3)*(r2(2))+r2(3)*(r2(2))+two*r2(3)*(r1(2))
             aa(6,12)=-(r1(1)*r2(1)+r1(2)*r2(2)+r2(1)*r2(1)+r2(2)*r2(2))
             !......................................
             aa(7,1:3)=aa(1:3,7)
             aa(7,4:6)=aa(4:6,7)
             aa(7,7)=-two*(r1(2)*r2(2)+r1(3)*r2(3))-two*(r1(2)*r3(2)+r1(3)*r3(3))
             aa(7,8)=r1(1)*(two*r2(2)+r3(2)-r2(2))-r2(1)*(r1(2))+r3(1)*(r1(2))- &
                     two*r2(1)*(-r1(2))
             aa(7,9)=r1(1)*(two*r2(3)+r3(3)-r2(3))-r2(1)*(r1(3))+r3(1)*(r1(3))- &
                     two*r2(1)*(-r1(3))
             aa(7,10)=r1(2)*r2(2)+r1(3)*r2(3)
             aa(7,11)=r1(1)*(r2(2))-two*r2(1)*(r1(2))
             aa(7,12)=r1(1)*(r2(3))-two*r2(1)*(r1(3))
             !......................................
             aa(8,1:3)=aa(1:3,8)
             aa(8,4:6)=aa(4:6,8)
             aa(8,7)=r1(2)*(two*r2(1)+r3(1)-r2(1))-r2(2)*(r1(1))+r3(2)*(r1(1))- &
                     two*r2(2)*(-r1(1))
             aa(8,8)=-two*(r1(1)*r2(1)+r1(3)*r2(3))-two*(r1(1)*r3(1)+r1(3)*r3(3))
             aa(8,9)=r1(2)*(two*r2(3)+r3(3)-r2(3))-r2(2)*(r1(3))+r3(2)*(r1(3))- &
                     two*r2(2)*(-r1(3))
             aa(8,10)=r1(2)*(r2(1))-two*r2(2)*(r1(1))
             aa(8,11)=r1(1)*r2(1)+r1(3)*r2(3)
             aa(8,12)=r1(2)*(r2(3))-two*r2(2)*(r1(3))
             !......................................
             aa(9,1:3)=aa(1:3,9)
             aa(9,4:6)=aa(4:6,9)
             aa(9,7)=r1(3)*(two*r2(1)+r3(1)-r2(1))-r2(3)*(r1(1))+r3(3)*(r1(1))- &
                     two*r2(3)*(-r1(1))
             aa(9,8)=r1(3)*(two*r2(2)+r3(2)-r2(2))-r2(3)*(r1(2))+r3(3)*(r1(2))- &
                     two*r2(3)*(-r1(2))
             aa(9,9)=-two*(r1(1)*r2(1)+r1(2)*r2(2))-two*(r1(1)*r3(1)+r1(2)*r3(2))
             aa(9,10)=r1(3)*(r2(1))-two*r2(3)*(r1(1))
             aa(9,11)=r1(3)*(r2(2))-two*r2(3)*(r1(2))
             aa(9,12)=r1(1)*r2(1)+r1(2)*r2(2)
             !......................................
             aa(10,1:3)=aa(1:3,10)
             aa(10,4:6)=aa(4:6,10)
             aa(10,7:9)=aa(7:9,10)
             aa(10,10:12)=zero
             !......................................
             aa(11,1:3)=aa(1:3,11)
             aa(11,4:6)=aa(4:6,11)
             aa(11,7:9)=aa(7:9,11)
             aa(11,10:12)=zero
             !......................................
             aa(12,1:3)=aa(1:3,12)
             aa(12,4:6)=aa(4:6,12)
             aa(12,7:9)=aa(7:9,12)
             aa(12,10:12)=zero
             !......................................
             !......................................
             bb(1,1)=two*(r2(2)*r2(2)+r2(3)*r2(3))
             bb(1,2)=two*r2(1)*(-r2(2))
             bb(1,3)=two*r2(1)*(-r2(3))
             bb(1,4)=-two*(r2(2)*r2(2)+r2(3)*r2(3))-two*(r1(2)*r2(2)+r1(3)*r2(3))
             bb(1,5)=-two*r1(1)*(-two*r2(2))+two*r2(1)*(r2(2)-r1(2))
             bb(1,6)=-two*r1(1)*(-two*r2(3))+two*r2(1)*(r2(3)-r1(3))
             bb(1,7)=two*(r1(2)*r2(2)+r1(3)*r2(3))
             bb(1,8)=-two*r1(1)*(two*r2(2))+two*r2(1)*(r1(2))
             bb(1,9)=-two*r1(1)*(two*r2(3))+two*r2(1)*(r1(3))
             bb(1,10:12)=zero
             !......................................
             bb(2,1)=two*r2(2)*(-r2(1))
             bb(2,2)=two*(r2(1)*r2(1)+r2(3)*r2(3))
             bb(2,3)=two*r2(2)*(-r2(3))
             bb(2,4)=-two*r1(2)*(-two*r2(1))+two*r2(2)*(r2(1)-r1(1))
             bb(2,5)=-two*(r2(1)*r2(1)+r2(3)*r2(3))-two*(r1(1)*r2(1)+r1(3)*r2(3))
             bb(2,6)=-two*r1(2)*(-two*r2(3))+two*r2(2)*(r2(3)-r1(3))
             bb(2,7)=-two*r1(2)*(two*r2(1))+two*r2(2)*(r1(1))
             bb(2,8)=two*(r1(1)*r2(1)+r1(3)*r2(3))
             bb(2,9)=-two*r1(2)*(two*r2(3))+two*r2(2)*(r1(3))
             bb(2,10:12)=zero
             !......................................
             bb(3,1)=two*r2(3)*(-r2(1))
             bb(3,2)=two*r2(3)*(-r2(2))
             bb(3,3)=two*(r2(1)*r2(1)+r2(2)*r2(2))
             bb(3,4)=-two*r1(3)*(-two*r2(1))+two*r2(3)*(r2(1)-r1(1))
             bb(3,5)=-two*r1(3)*(-two*r2(2))+two*r2(3)*(r2(2)-r1(2))
             bb(3,6)=-two*(r2(1)*r2(1)+r2(2)*r2(2))-two*(r1(1)*r2(1)+r1(2)*r2(2))
             bb(3,7)=-two*r1(3)*(two*r2(1))+two*r2(3)*(r1(1))
             bb(3,8)=-two*r1(3)*(two*r2(2))+two*r2(3)*(r1(2))
             bb(3,9)=two*(r1(2)*r2(2)+r1(1)*r2(1))
             bb(3,10:12)=zero
             !......................................
             bb(4,1:3)=bb(1:3,4)
             bb(4,4)=two*(r2(2)*r2(2)+r2(3)*r2(3)+r1(2)*r2(2)+r1(3)*r2(3))+ &
                     two*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             bb(4,5)=two*r1(1)*(-two*r2(2)+r2(2)-r1(2))-two*r2(1)*(two*r1(2)+r2(2)-r1(2))
             bb(4,6)=two*r1(1)*(-two*r2(3)+r2(3)-r1(3))-two*r2(1)*(two*r1(3)+r2(3)-r1(3))
             bb(4,7)=-two*(r1(2)*r1(2)+r1(3)*r1(3)+r1(2)*r2(2)+r1(3)*r2(3))
             bb(4,8)=two*r1(1)*(two*r2(2)+r1(2))-two*r2(1)*(r1(2))
             bb(4,9)=two*r1(1)*(two*r2(3)+r1(3))-two*r2(1)*(r1(3))
             bb(4,10:12)=zero
             !......................................
             bb(5,1:3)=bb(1:3,5)
             bb(5,4)=two*r1(2)*(-two*r2(1)+r2(1)-r1(1))-two*r2(2)*(two*r1(1)+r2(1)-r1(1))
             bb(5,5)=two*(r2(1)*r2(1)+r2(3)*r2(3)+r1(1)*r2(1)+r1(3)*r2(3))+ &
                     two*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             bb(5,6)=two*r1(2)*(-two*r2(3)+r2(3)-r1(3))-two*r2(2)*(two*r1(3)+r2(3)-r1(3))
             bb(5,7)=two*r1(2)*(two*r2(1)+r1(1))-two*r2(2)*(r1(1))
             bb(5,8)=-two*(r1(1)*r1(1)+r1(3)*r1(3)+r1(1)*r2(1)+r1(3)*r2(3))
             bb(5,9)=two*r1(2)*(two*r2(3)+r1(3))-two*r2(2)*(r1(3))
             bb(5,10:12)=zero
             !......................................
             bb(6,1:3)=bb(1:3,6)
             bb(6,4)=two*r1(3)*(-two*r2(1)+r2(1)-r1(1))-two*r2(3)*(two*r1(1)+r2(1)-r1(1))
             bb(6,5)=two*r1(3)*(-two*r2(2)+r2(2)-r1(2))-two*r2(3)*(two*r1(2)+r2(2)-r1(2))
             bb(6,6)=two*(r2(1)*r2(1)+r2(2)*r2(2)+r1(1)*r2(1)+r1(2)*r2(2))+ &
                     two*(r1(1)*r1(1)+r1(2)*r1(2)+r1(1)*r2(1)+r1(2)*r2(2))
             bb(6,7)=two*r1(3)*(two*r2(1)+r1(1))-two*r2(3)*(r1(1))
             bb(6,8)=two*r1(3)*(two*r2(2)+r1(2))-two*r2(3)*(r1(2))
             bb(6,9)=-two*(r1(1)*r1(1)+r1(2)*r1(2)+r1(1)*r2(1)+r1(2)*r2(2))
             bb(6,10:12)=zero
             !......................................
             bb(7,1:3)=bb(1:3,7)
             bb(7,4:6)=bb(4:6,7)
             bb(7,7)=two*(r1(2)*r1(2)+r1(3)*r1(3))
             bb(7,8)=-two*r1(1)*(r1(2))
             bb(7,9)=-two*r1(1)*(r1(3))
             bb(7,10:12)=zero
             !......................................
             bb(8,1:3)=bb(1:3,8)
             bb(8,4:6)=bb(4:6,8)
             bb(8,7)=-two*r1(2)*(r1(1))
             bb(8,8)=two*(r1(1)*r1(1)+r1(3)*r1(3))
             bb(8,9)=-two*r1(2)*(r1(3))
             bb(8,10:12)=zero
             !......................................
             bb(9,1:3)=bb(1:3,9)
             bb(9,4:6)=bb(4:6,9)
             bb(9,7)=-two*r1(3)*(r1(1))
             bb(9,8)=-two*r1(3)*(r1(2))
             bb(9,9)=two*(r1(1)*r1(1)+r1(2)*r1(2))
             bb(9,10:12)=zero
             !......................................
             bb(10:12,1:12)=zero
             !......................................
             !......................................
             cc(1:3,1:12)=zero
             !......................................
             cc(4,1:3)=zero
             cc(4,4)=two*(r3(2)*r3(2)+r3(3)*r3(3))
             cc(4,5)=two*r3(1)*(-r3(2))
             cc(4,6)=two*r3(1)*(-r3(3))
             cc(4,7)=-two*(r2(2)*r3(2)+r2(3)*r3(3))-two*(r3(2)*r3(2)+r3(3)*r3(3))
             cc(4,8)=two*r3(1)*(r3(2)-r2(2))-two*r2(1)*(-two*r3(2))
             cc(4,9)=two*r3(1)*(r3(3)-r2(3))-two*r2(1)*(-two*r3(3))
             cc(4,10)=two*(r2(2)*r3(2)+r2(3)*r3(3))
             cc(4,11)=two*r3(1)*(r2(2))-two*r2(1)*(two*r3(2))
             cc(4,12)=two*r3(1)*(r2(3))-two*r2(1)*(two*r3(3))
             !......................................
             cc(5,1:3)=zero
             cc(5,4)=two*r3(2)*(-r3(1))
             cc(5,5)=two*(r3(1)*r3(1)+r3(3)*r3(3))
             cc(5,6)=two*r3(2)*(-r3(3))
             cc(5,7)=two*r3(2)*(r3(1)-r2(1))-two*r2(2)*(-two*r3(1))
             cc(5,8)=-two*(r2(1)*r3(1)+r2(3)*r3(3))-two*(r3(1)*r3(1)+r3(3)*r3(3))
             cc(5,9)=two*r3(2)*(r3(3)-r2(3))-two*r2(2)*(-two*r3(3))
             cc(5,10)=two*r3(2)*(r2(1))-two*r2(2)*(two*r3(1))
             cc(5,11)=two*(r2(1)*r3(1)+r2(3)*r3(3))
             cc(5,12)=two*r3(2)*(r2(3))-two*r2(2)*(two*r3(3))
             !......................................
             cc(6,1:3)=zero
             cc(6,4)=two*r3(3)*(-r3(1))
             cc(6,5)=two*r3(3)*(-r3(2))
             cc(6,6)=two*(r3(1)*r3(1)+r3(2)*r3(2))
             cc(6,7)=two*r3(3)*(r3(1)-r2(1))-two*r2(3)*(-two*r3(1))
             cc(6,8)=two*r3(3)*(r3(2)-r2(2))-two*r2(3)*(-two*r3(2))
             cc(6,9)=-two*(r2(1)*r3(1)+r2(2)*r3(2))-two*(r3(1)*r3(1)+r3(2)*r3(2))
             cc(6,10)=two*r3(3)*(r2(1))-two*r2(3)*(two*r3(1))
             cc(6,11)=two*r3(3)*(r2(2))-two*r2(3)*(two*r3(2))
             cc(6,12)=two*(r2(1)*r3(1)+r2(2)*r3(2))
             !......................................
             cc(7,1:3)=zero
             cc(7,4:6)=cc(4:6,7)
             cc(7,7)=two*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))+ &
                     two*(r3(2)*r3(2)+r3(3)*r3(3)+r2(2)*r3(2)+r2(3)*r3(3))
             cc(7,8)=-two*r3(1)*(two*r2(2)+r3(2)-r2(2))+two*r2(1)*(-two*r3(2)+r3(2)-r2(2))
             cc(7,9)=-two*r3(1)*(two*r2(3)+r3(3)-r2(3))+two*r2(1)*(-two*r3(3)+r3(3)-r2(3))
             cc(7,10)=-two*(r2(2)*r2(2)+r2(3)*r2(3)+r2(2)*r3(2)+r2(3)*r3(3))
             cc(7,11)=-two*r3(1)*(r2(2))+two*r2(1)*(two*r3(2)+r2(2))
             cc(7,12)=-two*r3(1)*(r2(3))+two*r2(1)*(two*r3(3)+r2(3))
             !......................................
             cc(8,1:3)=zero
             cc(8,4:6)=cc(4:6,8)
             cc(8,7)=-two*r3(2)*(two*r2(1)+r3(1)-r2(1))+two*r2(2)*(-two*r3(1)+r3(1)-r2(1))
             cc(8,8)=two*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))+ &
                     two*(r3(1)*r3(1)+r3(3)*r3(3)+r2(1)*r3(1)+r2(3)*r3(3))
             cc(8,9)=-two*r3(2)*(two*r2(3)+r3(3)-r2(3))+two*r2(2)*(-two*r3(3)+r3(3)-r2(3))
             cc(8,10)=-two*r3(2)*(r2(1))+two*r2(2)*(two*r3(1)+r2(1))
             cc(8,11)=-two*(r2(1)*r2(1)+r2(3)*r2(3)+r2(1)*r3(1)+r2(3)*r3(3))
             cc(8,12)=-two*r3(2)*(r2(3))+two*r2(2)*(two*r3(3)+r2(3))
             !......................................
             cc(9,1:3)=zero
             cc(9,4:6)=cc(4:6,9)
             cc(9,7)=-two*r3(3)*(two*r2(1)+r3(1)-r2(1))+two*r2(3)*(-two*r3(1)+r3(1)-r2(1))
             cc(9,8)=-two*r3(3)*(two*r2(2)+r3(2)-r2(2))+two*r2(3)*(-two*r3(2)+r3(2)-r2(2))
             cc(9,9)=two*(r2(1)*r2(1)+r2(2)*r2(2)+r2(1)*r3(1)+r2(2)*r3(2))+ &
                     two*(r3(1)*r3(1)+r3(2)*r3(2)+r2(1)*r3(1)+r2(2)*r3(2))
             cc(9,10)=-two*r3(3)*(r2(1))+two*r2(3)*(two*r3(1)+r2(1))
             cc(9,11)=-two*r3(3)*(r2(2))+two*r2(3)*(two*r3(2)+r2(2))
             cc(9,12)=-two*(r2(1)*r2(1)+r2(2)*r2(2)+r2(1)*r3(1)+r2(2)*r3(2))
             !......................................
             cc(10,1:3)=zero
             cc(10,4:6)=cc(4:6,10)
             cc(10,7:9)=cc(7:9,10)
             cc(10,10)=two*(r2(2)*r2(2)+r2(3)*r2(3))
             cc(10,11)=-two*r2(1)*(r2(2))
             cc(10,12)=-two*r2(1)*(r2(3))
             !......................................
             cc(11,1:3)=zero
             cc(11,4:6)=cc(4:6,11)
             cc(11,7:9)=cc(7:9,11)
             cc(11,10)=-two*r2(2)*(r2(1))
             cc(11,11)=two*(r2(1)*r2(1)+r2(3)*r2(3))
             cc(11,12)=-two*r2(2)*(r2(3))
             !......................................
             cc(12,1:3)=zero
             cc(12,4:6)=cc(4:6,12)
             cc(12,7:9)=cc(7:9,12)
             cc(12,10)=-two*r2(3)*(r2(1))
             cc(12,11)=-two*r2(3)*(r2(2))
             cc(12,12)=two*(r2(1)*r2(1)+r2(2)*r2(2))
             !......................................

             do i1=1,12
                kk=int((i1-1)/3)+1
                m=i1-3*(kk-1)
                if(kk==1) then
                   fs=ft2(1,:)
                else if(kk==2) then
                   fs=ft2(2,:)
                else if(kk==3) then
                   fs=ft2(3,:)
                else if(kk==4) then
                   fs=ft2(4,:)
                end if
                do j1=1,12
                   fft2(i1,j1)= fs(m)*sin_phi*d_sin(j1)- &
                        (aa(i1,j1)*m2132+A(kk,m)*g2132(j1)-(dcos(j1)/two)*(B(kk,m)*m2121+C(kk,m)*m3232)- &
                        (cos_phi/two)*(bb(i1,j1)*m2121+B(kk,m)*g2121(j1)+cc(i1,j1)*m3232+C(kk,m)*g3232(j1)))/sin_phi
                end do
             end do
             !......................................
             !......................................
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             k4=3*(ia4-1)
             do l=1,3
                do m=1,12
                   if(m <= 3) then
                      fst_1=ft1(1,m)
                      fst_2=ft2(1,m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      fst_1=ft1(2,m-3)
                      fst_2=ft2(2,m-3)
                      m1=k2+(m-3)
                   else if(m > 6 .and. m <= 9) then
                      fst_1=ft1(3,m-6)
                      fst_2=ft2(3,m-6)
                      m1=k3+(m-6)
                   else if(m > 9 .and. m <= 12) then
                      fst_1=ft1(4,m-9)
                      fst_2=ft2(4,m-9)
                      m1=k4+(m-9)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+(d2E_dr22*fst_1+d2E_dr2dp*fst_2)*ft1(1,l)+ &
                        (d2E_dpdr2*fst_1+d2E_dp2*fst_2)*ft2(1,l)+ &
                        dE_dr2*fft1(l,m)+dE_dp*fft2(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+(d2E_dr22*fst_1+d2E_dr2dp*fst_2)*ft1(2,l)+ &
                        (d2E_dpdr2*fst_1+d2E_dp2*fst_2)*ft2(2,l)+ &
                        dE_dr2*fft1(l+3,m)+dE_dp*fft2(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+(d2E_dr22*fst_1+d2E_dr2dp*fst_2)*ft1(3,l)+ &
                        (d2E_dpdr2*fst_1+d2E_dp2*fst_2)*ft2(3,l)+ &
                        dE_dr2*fft1(l+6,m)+dE_dp*fft2(l+6,m)
                   H(k4+l,m1)=H(k4+l,m1)+(d2E_dr22*fst_1+d2E_dr2dp*fst_2)*ft1(4,l)+ &
                        (d2E_dpdr2*fst_1+d2E_dp2*fst_2)*ft2(4,l)+ &
                        dE_dr2*fft1(l+9,m)+dE_dp*fft2(l+9,m)
                end do
             end do

          end if

       end do j_lenght3
    end do i_n_st
!!$print*,'4444444444'

!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine many_body_E_and_F
  !****************************************************************

  !****************************************************************
  function vector_product(v1,v2)

    real(kind=r8_kind) :: vector_product(3)
    real(kind=r8_kind) :: v1(3),v2(3)

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !****************************************************************

  !****************************************************************
end module covalent_module

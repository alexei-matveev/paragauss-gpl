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
module van_der_waals_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use tasks_main_options_module
  use slab_module
  use species_module
  use external_field_module
  use potentials_module
  use n_body_lists_module
  use energy_and_forces_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ------------------
  public vdw_E_and_F
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine vdw_E_and_F()

    integer(kind=i4_kind) :: i,j,jj,length,length1,ind,id,i1,k,kk,ns
    integer(kind=i4_kind) :: img,jj1,k1,k2,k3,k4,l,m,m1,l1,l2,l3,l4
    real(kind=r8_kind) :: r(3),dr,r1(3),r2(3),red_fac(2),rj(3),rj1(3),E_buf
    real(kind=r8_kind) :: p(n_parameter),f1(3),f2(3),E_r,E_d
    real(kind=r8_kind) :: dE_dr,dE_r_dr,dE_d_dr,V,VV,HH
    real(kind=r8_kind) :: d2E_dr2,ff(12,12),fstore,a1,a2,d2E_r_dr2,d2E_d_dr2
    integer(kind=i4_kind) :: iup,in,ty,count
    real(kind=r8_kind) :: eps,eps2,eps3,y1,y2,y1_2,y2_2,st,sfac
    real(r8_kind) :: dr3,dr6,dr7,dr8,dr12,dr13,dr14,st2

    ns=3*n_species

    if(n_2b_n==0) return
    do i=1,n_species
       if(.not. associated(non_bonded(i)%list)) cycle
       length=size(non_bonded(i)%list,2)
       if(length==0) cycle
       
       ty=atoms_cart(i)%type
       red_fac(1)=atoms(ty)%redfac
       i1=atoms_cart(i)%neigh_sp
       r1=red_fac(1)*(atoms_cart(i)%r-atoms_cart(i1)%r)+atoms_cart(i1)%r
       do j=1,length
          jj=non_bonded(i)%list(1,j)
          if(jj <= n_species) then
             ty=atoms_cart(jj)%type
             red_fac(2)=atoms(ty)%redfac
             jj1=atoms_cart(jj)%neigh_sp
             rj=atoms_cart(jj)%r; rj1=atoms_cart(jj1)%r
          else
             ty=pc_cart(jj-n_species)%type
             red_fac(2)=one
             jj1=jj
             rj=pc_cart(jj-n_species)%r; rj1=pc_cart(jj1-n_species)%r
          end if

          if(non_bonded(i)%first_image(j)) then
             !calculation vdw intractions between atoms in real unit cell
             !if minimal image = true, there are not interactions between
             !exlcuded atoms
             r2=red_fac(2)*(rj-rj1)+rj1

!!$          r=atoms_cart(jj)%r-atoms_cart(i)%r
             r=r2-r1
             if(lattice_calc .and. minimal_image) then 
                call image(r)
             else if(slab_calc .and. minimal_image) then
                call image_slab(r)
             end if
             dr=sqrt(dot_product(r,r))

             dr3=dr*dr*dr
             dr6=dr3*dr3
             dr7=dr6*dr
             dr8=dr7*dr
             dr12=dr6*dr6
             dr13=dr12*dr
             dr14=dr13*dr

             ind=non_bonded(i)%list(2,j)
             id=poten(ind)%id
             p=poten(ind)%param
             select case(id)
             case (4) !buck
                E_buf=p(1)*exp(-dr/p(2))-p(3)/dr6
                E(4)=E(4)+E_buf
!!$print*,E_buf*j2c,E(4)*j2c,'=============='
                if(calc_gradients) dE_dr=-(p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7
                if(calc_hessian) d2E_dr2=(p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8
             case (5) !l_j
                E_buf=p(1)/dr12-p(2)/dr6
                E(4)=E(4)+E_buf
                if(calc_gradients) dE_dr=-twelve*p(1)/dr13+six*p(2)/dr7
                if(calc_hessian) d2E_dr2=156.0_r8_kind*p(1)/dr14-42.0_r8_kind*p(3)/dr8
             case (n_poten+1 :) !user_def
                iup=id-n_poten
                st=user_poten(iup)%step
                st2=st*st
                in=int((dr-user_poten(iup)%r(1))/st)+1
                eps=(dr-user_poten(iup)%r(in))/st
                eps2=eps*eps
                eps3=eps2*eps

                y1=user_poten(iup)%f(in)
                y2=user_poten(iup)%f(in+1)
                y1_2=user_poten(iup)%d2f(in)
                y2_2=user_poten(iup)%d2f(in+1)
                E_r=p(1)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                     eps3*st2*(y2_2-y1_2)/six)
                if(calc_gradients) dE_r_dr=p(1)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                     two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                if(calc_hessian) d2E_r_dr2= &
                     p(1)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                y1=user_poten(iup)%f1(in)
                y2=user_poten(iup)%f1(in+1)
                y1_2=user_poten(iup)%d2f1(in)
                y2_2=user_poten(iup)%d2f1(in+1)
                E_d=p(2)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                     eps3*st2*(y2_2-y1_2)/six)
                if(calc_gradients) dE_d_dr=p(2)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                     two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                if(calc_hessian) d2E_d_dr2= &
                     p(2)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                E_buf=E_r+E_d
                E(id)=E(id)+E_buf
                if(calc_gradients) dE_dr=dE_r_dr+dE_d_dr
                if(calc_hessian) d2E_dr2=d2E_r_dr2+d2E_d_dr2
             end select
             E_total=E_total+E_buf

             if(calc_gradients) then
                f1=r/dr
                f2=-f1
                Grad(:,i)=Grad(:,i)+dE_dr*f2*(red_fac(1))
                if(jj <= n_species) Grad(:,jj)=Grad(:,jj)+dE_dr*f1*(red_fac(2))
                Grad(:,i1)=Grad(:,i1)+dE_dr*f2*(-red_fac(1)+one)
                if(jj1 <= n_species) Grad(:,jj1)=Grad(:,jj1)+dE_dr*f1*(-red_fac(2)+one)
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
                a1=red_fac(1); a2=red_fac(2)

                ff(1,1)=a1*a1/dr-a1*a1*r(1)*r(1)/dr3
                ff(1,2)=-a1*a1*r(1)*r(2)/dr3
                ff(1,3)=-a1*a1*r(1)*r(3)/dr3
                ff(1,4)=-a1*a2/dr+a1*a2*r(1)*r(1)/dr3
                ff(1,5)=a1*a2*r(1)*r(2)/dr3
                ff(1,6)=a1*a2*r(1)*r(3)/dr3
                ff(1,7)=-a1*(a1-one)/dr+a1*(a1-one)*r(1)*r(1)/dr3
                ff(1,8)=a1*(a1-one)*r(1)*r(2)/dr3
                ff(1,9)=a1*(a1-one)*r(1)*r(3)/dr3
                ff(1,10)=-a1*(one-a2)/dr+a1*(one-a2)*r(1)*r(1)/dr3
                ff(1,11)=a1*(one-a2)*r(1)*r(2)/dr3
                ff(1,12)=a1*(one-a2)*r(1)*r(3)/dr3

                ff(2,1)=-a1*a1*r(2)*r(1)/dr3
                ff(2,2)=a1*a1/dr-a1*a1*r(2)*r(2)/dr3
                ff(2,3)=-a1*a1*r(2)*r(3)/dr3
                ff(2,4)=a1*a2*r(2)*r(1)/dr3
                ff(2,5)=-a1*a2/dr+a1*a2*r(2)*r(2)/dr3
                ff(2,6)=a1*a2*r(2)*r(3)/dr3
                ff(2,7)=a1*(a1-one)*r(2)*r(1)/dr3
                ff(2,8)=-a1*(a1-one)/dr+a1*(a1-one)*r(2)*r(2)/dr3
                ff(2,9)=a1*(a1-one)*r(2)*r(3)/dr3
                ff(2,10)=a1*(one-a2)*r(2)*r(1)/dr3
                ff(2,11)=-a1*(one-a2)/dr+a1*(one-a2)*r(2)*r(2)/dr3
                ff(2,12)=a1*(one-a2)*r(2)*r(3)/dr3

                ff(3,1)=-a1*a1*r(3)*r(1)/dr3
                ff(3,2)=-a1*a1*r(3)*r(2)/dr3
                ff(3,3)=a1*a1/dr-a1*a1*r(3)*r(3)/dr3
                ff(3,4)=a1*a2*r(3)*r(1)/dr3
                ff(3,5)=a1*a2*r(3)*r(2)/dr3
                ff(3,6)=-a1*a2/dr+a1*a2*r(3)*r(3)/dr3
                ff(3,7)=a1*(a1-one)*r(3)*r(1)/dr3
                ff(3,8)=a1*(a1-one)*r(3)*r(2)/dr3
                ff(3,9)=-a1*(a1-one)/dr+a1*(a1-one)*r(3)*r(3)/dr3
                ff(3,10)=a1*(one-a2)*r(3)*r(1)/dr3
                ff(3,11)=a1*(one-a2)*r(3)*r(2)/dr3
                ff(3,12)=-a1*(one-a2)/dr+a1*(one-a2)*r(3)*r(3)/dr3

                ff(4,1:3)=ff(1:3,4)
                ff(4,4)=a2*a2/dr-a2*a2*r(1)*r(1)/dr3
                ff(4,5)=-a2*a2*r(1)*r(2)/dr3
                ff(4,6)=-a2*a2*r(1)*r(3)/dr3
                ff(4,7)=a2*(a1-one)/dr-a2*(a1-one)*r(1)*r(1)/dr3
                ff(4,8)=-a2*(a1-one)*r(1)*r(2)/dr3
                ff(4,9)=-a2*(a1-one)*r(1)*r(3)/dr3
                ff(4,10)=a2*(one-a2)/dr-a2*(one-a2)*r(1)*r(1)/dr3
                ff(4,11)=-a2*(one-a2)*r(1)*r(2)/dr3
                ff(4,12)=-a2*(one-a2)*r(1)*r(3)/dr3

                ff(5,1:3)=ff(1:3,5)
                ff(5,4)=-a2*a2*r(2)*r(1)/dr3
                ff(5,5)=a2*a2/dr-a2*a2*r(2)*r(2)/dr3
                ff(5,6)=-a2*a2*r(2)*r(3)/dr3
                ff(5,7)=-a2*(a1-one)*r(2)*r(1)/dr3
                ff(5,8)=a2*(a1-one)/dr-a2*(a1-one)*r(2)*r(2)/dr3
                ff(5,9)=-a2*(a1-one)*r(2)*r(3)/dr3
                ff(5,10)=-a2*(one-a2)*r(2)*r(1)/dr3
                ff(5,11)=a2*(one-a2)/dr-a2*(one-a2)*r(2)*r(2)/dr3
                ff(5,12)=-a2*(one-a2)*r(2)*r(3)/dr3

                ff(6,1:3)=ff(1:3,6)
                ff(6,4)=-a2*a2*r(3)*r(1)/dr3
                ff(6,5)=-a2*a2*r(3)*r(2)/dr3
                ff(6,6)=a2*a2/dr-a2*a2*r(3)*r(3)/dr3
                ff(6,7)=-a2*(a1-one)*r(3)*r(1)/dr3
                ff(6,8)=-a2*(a1-one)*r(3)*r(2)/dr3
                ff(6,9)=a2*(a1-one)/dr-a2*(a1-one)*r(3)*r(3)/dr3
                ff(6,10)=-a2*(one-a2)*r(3)*r(1)/dr3
                ff(6,11)=-a2*(one-a2)*r(3)*r(2)/dr3
                ff(6,12)=a2*(one-a2)/dr-a2*(one-a2)*r(3)*r(3)/dr3

                ff(7,1:3)=ff(1:3,7)
                ff(7,4:6)=ff(4:6,7)
                ff(7,7)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(1)*r(1)/dr3
                ff(7,8)=-(a1-one)*(a1-one)*r(1)*r(2)/dr3
                ff(7,9)=-(a1-one)*(a1-one)*r(1)*r(3)/dr3
                ff(7,10)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(1)*r(1)/dr3
                ff(7,11)=-(a1-one)*(one-a2)*r(1)*r(2)/dr3
                ff(7,12)=-(a1-one)*(one-a2)*r(1)*r(3)/dr3

                ff(8,1:3)=ff(1:3,8)
                ff(8,4:6)=ff(4:6,8)
                ff(8,7)=-(a1-one)*(a1-one)*r(2)*r(1)/dr3
                ff(8,8)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(2)*r(2)/dr3
                ff(8,9)=-(a1-one)*(a1-one)*r(2)*r(3)/dr3
                ff(8,10)=-(a1-one)*(one-a2)*r(2)*r(1)/dr3
                ff(8,11)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(2)*r(2)/dr3
                ff(8,12)=-(a1-one)*(one-a2)*r(2)*r(3)/dr3

                ff(9,1:3)=ff(1:3,9)
                ff(9,4:6)=ff(4:6,9)
                ff(9,7)=-(a1-one)*(a1-one)*r(3)*r(1)/dr3
                ff(9,8)=-(a1-one)*(a1-one)*r(3)*r(2)/dr3
                ff(9,9)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(3)*r(3)/dr3
                ff(9,10)=-(a1-one)*(one-a2)*r(3)*r(1)/dr3
                ff(9,11)=-(a1-one)*(one-a2)*r(3)*r(2)/dr3
                ff(9,12)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(3)*r(3)/dr3
             
                ff(10,1:3)=ff(1:3,10)
                ff(10,4:6)=ff(4:6,10)
                ff(10,7:9)=ff(7:9,10)
                ff(10,10)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(1)*r(1)/dr3
                ff(10,11)=-(one-a2)*(one-a2)*r(1)*r(2)/dr3
                ff(10,12)=-(one-a2)*(one-a2)*r(1)*r(3)/dr3

                ff(11,1:3)=ff(1:3,11)
                ff(11,4:6)=ff(4:6,11)
                ff(11,7:9)=ff(7:9,11)
                ff(11,10)=-(one-a2)*(one-a2)*r(2)*r(1)/dr3
                ff(11,11)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(2)*r(2)/dr3
                ff(11,12)=-(one-a2)*(one-a2)*r(2)*r(3)/dr3

                ff(12,1:3)=ff(1:3,12)
                ff(12,4:6)=ff(4:6,12)
                ff(12,7:9)=ff(7:9,12)
                ff(12,10)=-(one-a2)*(one-a2)*r(3)*r(1)/dr3
                ff(12,11)=-(one-a2)*(one-a2)*r(3)*r(2)/dr3
                ff(12,12)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(3)*r(3)/dr3

                k1=3*(i-1)
                k2=3*(jj-1)
                k3=3*(i1-1)
                k4=3*(jj1-1)
                do l=1,3
                   do m=1,12
                      if(m <= 3) then
                         fstore=f2(m)*a1
                         m1=k1+m
                      else if(m > 3 .and. m <= 6) then 
                         m1=k2+(m-3)
                         fstore=f1(m-3)*a2
                      else if(m > 6 .and. m <= 9) then
                         m1=k3+(m-6)
                         fstore=f2(m-6)*(one-a1)
                      else if(m > 9 .and. m <= 12) then
                         m1=k4+(m-9)
                         fstore=f1(m-9)*(one-a2)
                      end if
                      l1=k1+l; l2=k2+l; l3=k3+l; l4=k4+l
                      if(l1<=ns .and. m1<=ns) H(l1,m1)=H(l1,m1)+d2E_dr2*fstore*f2(l)*a1+dE_dr*ff(l,m)
                      if(l2<=ns .and. m1<=ns) H(l2,m1)=H(l2,m1)+d2E_dr2*fstore*f1(l)*a2+dE_dr*ff(l+3,m)
                      if(l3<=ns .and. m1<=ns) H(l3,m1)=H(l3,m1)+d2E_dr2*fstore*f2(l)*(one-a1)+dE_dr*ff(l+6,m)
                      if(l4<=ns .and. m1<=ns) H(l4,m1)=H(l4,m1)+d2E_dr2*fstore*f1(l)*(one-a2)+dE_dr*ff(l+9,m)
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
                   k1=3*(i-1); k2=3*(jj-1)

                   if(lattice_calc) then
                      HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                      H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                      H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(2)
                      H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                      H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                      HH=-VV*r(1)*r(3)*r(3)
                      H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                      H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+1,ns+4)=H(k1+1,ns+4)+HH; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                      H(k2+1,ns+4)=H(k2+1,ns+4)-HH; H(ns+4,k2+1)=H(ns+4,k2+1)-HH
                      HH=-(V*r(3)+VV*r(1)*r(1)*r(3))
                      H(k1+1,ns+5)=H(k1+1,ns+5)+HH; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                      H(k2+1,ns+5)=H(k2+1,ns+5)-HH; H(ns+5,k2+1)=H(ns+5,k2+1)-HH
                      HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                      H(k1+1,ns+6)=H(k1+1,ns+6)+HH; H(ns+6,k1+1)=H(ns+6,k1+1)+HH
                      H(k2+1,ns+6)=H(k2+1,ns+6)-HH; H(ns+6,k2+1)=H(ns+6,k2+1)-HH
                
                      HH=-VV*r(2)*r(1)*r(1)
                      H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                      H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                      HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                      H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                      H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                      HH=-VV*r(2)*r(3)*r(3)
                      H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                      H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH
                      HH=-(V*r(3)+VV*r(2)*r(2)*r(3))
                      H(k1+2,ns+4)=H(k1+2,ns+4)+HH; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                      H(k2+2,ns+4)=H(k2+2,ns+4)-HH; H(ns+4,k2+2)=H(ns+4,k2+2)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+2,ns+5)=H(k1+2,ns+5)+HH; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                      H(k2+2,ns+5)=H(k2+2,ns+5)-HH; H(ns+5,k2+2)=H(ns+5,k2+2)-HH
                      HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                      H(k1+2,ns+6)=H(k1+2,ns+6)+HH; H(ns+6,k1+2)=H(ns+6,k1+2)+HH
                      H(k2+2,ns+6)=H(k2+2,ns+6)-HH; H(ns+6,k2+2)=H(ns+6,k2+2)-HH
                
                      HH=-VV*r(3)*r(1)*r(1)
                      H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                      H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                      HH=-VV*r(3)*r(2)*r(2)
                      H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                      H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                      HH=-(V*two*r(3)+VV*r(3)*r(3)*r(3))
                      H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                      H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                      HH=-(V*r(2)+VV*r(3)*r(2)*r(3))
                      H(k1+3,ns+4)=H(k1+3,ns+4)+HH; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                      H(k2+3,ns+4)=H(k2+3,ns+4)-HH; H(ns+4,k2+3)=H(ns+4,k2+3)-HH
                      HH=-(V*r(1)+VV*r(3)*r(1)*r(3))
                      H(k1+3,ns+5)=H(k1+3,ns+5)+HH; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                      H(k2+3,ns+5)=H(k2+3,ns+5)-HH; H(ns+5,k2+3)=H(ns+5,k2+3)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+3,ns+6)=H(k1+3,ns+6)+HH; H(ns+6,k1+3)=H(ns+6,k1+3)+HH
                      H(k2+3,ns+6)=H(k2+3,ns+6)-HH; H(ns+6,k2+3)=H(ns+6,k2+3)-HH
                   else if(slab_calc) then
                      HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                      H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                      H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(2)
                      H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                      H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                      HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                      H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                      H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH

                      HH=-VV*r(2)*r(1)*r(1)
                      H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                      H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                      HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                      H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                      H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                      HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                      H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                      H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH

                      HH=-VV*r(3)*r(1)*r(1)
                      H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                      H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                      HH=-VV*r(3)*r(2)*r(2)
                      H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                      H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                      H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                   end if
                end if
             end if
          end if

          if(.not.minimal_image) then
             !calculation vdw intractions between atoms in real unit cell and unit cell imeges
             sfac=one; if(i==jj) sfac=half
             count=0
             do 
                count=count+1
                img=add_list(i)%list(count,j)
                if(img==0) exit
                r2=red_fac(2)*(im_coor(jj)%r(:,img)-im_coor(jj1)%r(:,img))+im_coor(jj1)%r(:,img)
!!$                r=im_coor(jj)%r(:,img)-atoms_cart(i)%r
                r=r2-r1
                dr=sqrt(dot_product(r,r))

                dr3=dr*dr*dr
                dr6=dr3*dr3
                dr7=dr6*dr
                dr8=dr7*dr
                dr12=dr6*dr6
                dr13=dr12*dr
                dr14=dr13*dr

                ind=non_bonded(i)%list(2,j)
                id=poten(ind)%id
                p=poten(ind)%param
                select case(id)
                case (4) !buck
                   E_buf=sfac*(p(1)*exp(-dr/p(2))-p(3)/dr6)
                   E(4)=E(4)+E_buf
!!$print*,E_buf*j2c,E(4)*j2c,'=============='
                   if(calc_gradients) dE_dr=sfac*(-(p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7)
                   if(calc_hessian) d2E_dr2=sfac*((p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8)
                case (5) !l_j
                   E_buf=sfac*(p(1)/dr12-p(2)/dr6)
                   E(4)=E(4)+E_buf
                   if(calc_gradients) dE_dr=sfac*(-twelve*p(1)/dr13+six*p(2)/dr7)
                   if(calc_hessian) d2E_dr2=sfac*(156.0_r8_kind*p(1)/dr14-42.0_r8_kind*p(3)/dr8)
                case (n_poten+1 :) !user_def
                   iup=id-n_poten
                   st=user_poten(iup)%step
                   st2=st*st
                   in=int((dr-user_poten(iup)%r(1))/st)+1
                   eps=(dr-user_poten(iup)%r(in))/st
                   eps2=eps*eps
                   eps3=eps2*eps

                   y1=user_poten(iup)%f(in)
                   y2=user_poten(iup)%f(in+1)
                   y1_2=user_poten(iup)%d2f(in)
                   y2_2=user_poten(iup)%d2f(in+1)
                   E_r=p(1)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                        eps3*st2*(y2_2-y1_2)/six)
                   if(calc_gradients) dE_r_dr=p(1)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                        two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                   if(calc_hessian) d2E_r_dr2=p(1)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                   y1=user_poten(iup)%f1(in)
                   y2=user_poten(iup)%f1(in+1)
                   y1_2=user_poten(iup)%d2f1(in)
                   y2_2=user_poten(iup)%d2f1(in+1)
                   E_d=p(2)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                        eps3*st2*(y2_2-y1_2)/six)
                   if(calc_gradients) dE_d_dr=p(2)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                        two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                   if(calc_hessian) d2E_d_dr2=p(2)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                   E_buf=sfac*(E_r+E_d)
                   E(id)=E(id)+E_buf
                   if(calc_gradients) dE_dr=sfac*(dE_r_dr+dE_d_dr)
                   if(calc_hessian) d2E_dr2=sfac*(d2E_r_dr2+d2E_d_dr2)
                end select
                E_total=E_total+E_buf

                if(i /= jj) then
                   if(calc_gradients) then
                      f1=r/dr
                      f2=-f1
                      Grad(:,i)=Grad(:,i)+dE_dr*f2*(red_fac(1))
                      Grad(:,jj)=Grad(:,jj)+dE_dr*f1*(red_fac(2))
                      Grad(:,i1)=Grad(:,i1)+dE_dr*f2*(-red_fac(1)+one)
                      Grad(:,jj1)=Grad(:,jj1)+dE_dr*f1*(-red_fac(2)+one)
                   end if

                   if(calc_hessian) then
                      a1=red_fac(1); a2=red_fac(2)

                      ff(1,1)=a1*a1/dr-a1*a1*r(1)*r(1)/dr3
                      ff(1,2)=-a1*a1*r(1)*r(2)/dr3
                      ff(1,3)=-a1*a1*r(1)*r(3)/dr3
                      ff(1,4)=-a1*a2/dr+a1*a2*r(1)*r(1)/dr3
                      ff(1,5)=a1*a2*r(1)*r(2)/dr3
                      ff(1,6)=a1*a2*r(1)*r(3)/dr3
                      ff(1,7)=-a1*(a1-one)/dr+a1*(a1-one)*r(1)*r(1)/dr3
                      ff(1,8)=a1*(a1-one)*r(1)*r(2)/dr3
                      ff(1,9)=a1*(a1-one)*r(1)*r(3)/dr3
                      ff(1,10)=-a1*(one-a2)/dr+a1*(one-a2)*r(1)*r(1)/dr3
                      ff(1,11)=a1*(one-a2)*r(1)*r(2)/dr3
                      ff(1,12)=a1*(one-a2)*r(1)*r(3)/dr3

                      ff(2,1)=-a1*a1*r(2)*r(1)/dr3
                      ff(2,2)=a1*a1/dr-a1*a1*r(2)*r(2)/dr3
                      ff(2,3)=-a1*a1*r(2)*r(3)/dr3
                      ff(2,4)=a1*a2*r(2)*r(1)/dr3
                      ff(2,5)=-a1*a2/dr+a1*a2*r(2)*r(2)/dr3
                      ff(2,6)=a1*a2*r(2)*r(3)/dr3
                      ff(2,7)=a1*(a1-one)*r(2)*r(1)/dr3
                      ff(2,8)=-a1*(a1-one)/dr+a1*(a1-one)*r(2)*r(2)/dr3
                      ff(2,9)=a1*(a1-one)*r(2)*r(3)/dr3
                      ff(2,10)=a1*(one-a2)*r(2)*r(1)/dr3
                      ff(2,11)=-a1*(one-a2)/dr+a1*(one-a2)*r(2)*r(2)/dr3
                      ff(2,12)=a1*(one-a2)*r(2)*r(3)/dr3

                      ff(3,1)=-a1*a1*r(3)*r(1)/dr3
                      ff(3,2)=-a1*a1*r(3)*r(2)/dr3
                      ff(3,3)=a1*a1/dr-a1*a1*r(3)*r(3)/dr3
                      ff(3,4)=a1*a2*r(3)*r(1)/dr3
                      ff(3,5)=a1*a2*r(3)*r(2)/dr3
                      ff(3,6)=-a1*a2/dr+a1*a2*r(3)*r(3)/dr3
                      ff(3,7)=a1*(a1-one)*r(3)*r(1)/dr3
                      ff(3,8)=a1*(a1-one)*r(3)*r(2)/dr3
                      ff(3,9)=-a1*(a1-one)/dr+a1*(a1-one)*r(3)*r(3)/dr3
                      ff(3,10)=a1*(one-a2)*r(3)*r(1)/dr3
                      ff(3,11)=a1*(one-a2)*r(3)*r(2)/dr3
                      ff(3,12)=-a1*(one-a2)/dr+a1*(one-a2)*r(3)*r(3)/dr3

                      ff(4,1:3)=ff(1:3,4)
                      ff(4,4)=a2*a2/dr-a2*a2*r(1)*r(1)/dr3
                      ff(4,5)=-a2*a2*r(1)*r(2)/dr3
                      ff(4,6)=-a2*a2*r(1)*r(3)/dr3
                      ff(4,7)=a2*(a1-one)/dr-a2*(a1-one)*r(1)*r(1)/dr3
                      ff(4,8)=-a2*(a1-one)*r(1)*r(2)/dr3
                      ff(4,9)=-a2*(a1-one)*r(1)*r(3)/dr3
                      ff(4,10)=a2*(one-a2)/dr-a2*(one-a2)*r(1)*r(1)/dr3
                      ff(4,11)=-a2*(one-a2)*r(1)*r(2)/dr3
                      ff(4,12)=-a2*(one-a2)*r(1)*r(3)/dr3

                      ff(5,1:3)=ff(1:3,5)
                      ff(5,4)=-a2*a2*r(2)*r(1)/dr3
                      ff(5,5)=a2*a2/dr-a2*a2*r(2)*r(2)/dr3
                      ff(5,6)=-a2*a2*r(2)*r(3)/dr3
                      ff(5,7)=-a2*(a1-one)*r(2)*r(1)/dr3
                      ff(5,8)=a2*(a1-one)/dr-a2*(a1-one)*r(2)*r(2)/dr3
                      ff(5,9)=-a2*(a1-one)*r(2)*r(3)/dr3
                      ff(5,10)=-a2*(one-a2)*r(2)*r(1)/dr3
                      ff(5,11)=a2*(one-a2)/dr-a2*(one-a2)*r(2)*r(2)/dr3
                      ff(5,12)=-a2*(one-a2)*r(2)*r(3)/dr3

                      ff(6,1:3)=ff(1:3,6)
                      ff(6,4)=-a2*a2*r(3)*r(1)/dr3
                      ff(6,5)=-a2*a2*r(3)*r(2)/dr3
                      ff(6,6)=a2*a2/dr-a2*a2*r(3)*r(3)/dr3
                      ff(6,7)=-a2*(a1-one)*r(3)*r(1)/dr3
                      ff(6,8)=-a2*(a1-one)*r(3)*r(2)/dr3
                      ff(6,9)=a2*(a1-one)/dr-a2*(a1-one)*r(3)*r(3)/dr3
                      ff(6,10)=-a2*(one-a2)*r(3)*r(1)/dr3
                      ff(6,11)=-a2*(one-a2)*r(3)*r(2)/dr3
                      ff(6,12)=a2*(one-a2)/dr-a2*(one-a2)*r(3)*r(3)/dr3

                      ff(7,1:3)=ff(1:3,7)
                      ff(7,4:6)=ff(4:6,7)
                      ff(7,7)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(1)*r(1)/dr3
                      ff(7,8)=-(a1-one)*(a1-one)*r(1)*r(2)/dr3
                      ff(7,9)=-(a1-one)*(a1-one)*r(1)*r(3)/dr3
                      ff(7,10)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(1)*r(1)/dr3
                      ff(7,11)=-(a1-one)*(one-a2)*r(1)*r(2)/dr3
                      ff(7,12)=-(a1-one)*(one-a2)*r(1)*r(3)/dr3

                      ff(8,1:3)=ff(1:3,8)
                      ff(8,4:6)=ff(4:6,8)
                      ff(8,7)=-(a1-one)*(a1-one)*r(2)*r(1)/dr3
                      ff(8,8)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(2)*r(2)/dr3
                      ff(8,9)=-(a1-one)*(a1-one)*r(2)*r(3)/dr3
                      ff(8,10)=-(a1-one)*(one-a2)*r(2)*r(1)/dr3
                      ff(8,11)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(2)*r(2)/dr3
                      ff(8,12)=-(a1-one)*(one-a2)*r(2)*r(3)/dr3

                      ff(9,1:3)=ff(1:3,9)
                      ff(9,4:6)=ff(4:6,9)
                      ff(9,7)=-(a1-one)*(a1-one)*r(3)*r(1)/dr3
                      ff(9,8)=-(a1-one)*(a1-one)*r(3)*r(2)/dr3
                      ff(9,9)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(3)*r(3)/dr3
                      ff(9,10)=-(a1-one)*(one-a2)*r(3)*r(1)/dr3
                      ff(9,11)=-(a1-one)*(one-a2)*r(3)*r(2)/dr3
                      ff(9,12)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(3)*r(3)/dr3
             
                      ff(10,1:3)=ff(1:3,10)
                      ff(10,4:6)=ff(4:6,10)
                      ff(10,7:9)=ff(7:9,10)
                      ff(10,10)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(1)*r(1)/dr3
                      ff(10,11)=-(one-a2)*(one-a2)*r(1)*r(2)/dr3
                      ff(10,12)=-(one-a2)*(one-a2)*r(1)*r(3)/dr3

                      ff(11,1:3)=ff(1:3,11)
                      ff(11,4:6)=ff(4:6,11)
                      ff(11,7:9)=ff(7:9,11)
                      ff(11,10)=-(one-a2)*(one-a2)*r(2)*r(1)/dr3
                      ff(11,11)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(2)*r(2)/dr3
                      ff(11,12)=-(one-a2)*(one-a2)*r(2)*r(3)/dr3

                      ff(12,1:3)=ff(1:3,12)
                      ff(12,4:6)=ff(4:6,12)
                      ff(12,7:9)=ff(7:9,12)
                      ff(12,10)=-(one-a2)*(one-a2)*r(3)*r(1)/dr3
                      ff(12,11)=-(one-a2)*(one-a2)*r(3)*r(2)/dr3
                      ff(12,12)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(3)*r(3)/dr3

                      k1=3*(i-1)
                      k2=3*(jj-1)
                      k3=3*(i1-1)
                      k4=3*(jj1-1)
                      do l=1,3
                         do m=1,12
                            if(m <= 3) then
                               fstore=f2(m)*a1
                               m1=k1+m
                            else if(m > 3 .and. m <= 6) then 
                               m1=k2+(m-3)
                               fstore=f1(m-3)*a2
                            else if(m > 6 .and. m <= 9) then
                               m1=k3+(m-6)
                               fstore=f2(m-6)*(one-a1)
                            else if(m > 9 .and. m <= 12) then
                               m1=k4+(m-9)
                               fstore=f1(m-9)*(one-a2)
                            end if
                            H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)*a1+dE_dr*ff(l,m)
                            H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)*a2+dE_dr*ff(l+3,m)
                            H(k3+l,m1)=H(k3+l,m1)+d2E_dr2*fstore*f2(l)*(one-a1)+dE_dr*ff(l+6,m)
                            H(k4+l,m1)=H(k4+l,m1)+d2E_dr2*fstore*f1(l)*(one-a2)+dE_dr*ff(l+9,m)
                         end do
                      end do
                   end if
                end if

                if(calc_strain .and. calc_gradients) then
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
                
                if(calc_strain .and. calc_hessian) then
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
                   k1=3*(i-1); k2=3*(jj-1)

                   if(lattice_calc) then
                      HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                      H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                      H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(2)
                      H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                      H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                      HH=-VV*r(1)*r(3)*r(3)
                      H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                      H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+1,ns+4)=H(k1+1,ns+4)+HH; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                      H(k2+1,ns+4)=H(k2+1,ns+4)-HH; H(ns+4,k2+1)=H(ns+4,k2+1)-HH
                      HH=-(V*r(3)+VV*r(1)*r(1)*r(3))
                      H(k1+1,ns+5)=H(k1+1,ns+5)+HH; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                      H(k2+1,ns+5)=H(k2+1,ns+5)-HH; H(ns+5,k2+1)=H(ns+5,k2+1)-HH
                      HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                      H(k1+1,ns+6)=H(k1+1,ns+6)+HH; H(ns+6,k1+1)=H(ns+6,k1+1)+HH
                      H(k2+1,ns+6)=H(k2+1,ns+6)-HH; H(ns+6,k2+1)=H(ns+6,k2+1)-HH
                
                      HH=-VV*r(2)*r(1)*r(1)
                      H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                      H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                      HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                      H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                      H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                      HH=-VV*r(2)*r(3)*r(3)
                      H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                      H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH
                      HH=-(V*r(3)+VV*r(2)*r(2)*r(3))
                      H(k1+2,ns+4)=H(k1+2,ns+4)+HH; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                      H(k2+2,ns+4)=H(k2+2,ns+4)-HH; H(ns+4,k2+2)=H(ns+4,k2+2)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+2,ns+5)=H(k1+2,ns+5)+HH; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                      H(k2+2,ns+5)=H(k2+2,ns+5)-HH; H(ns+5,k2+2)=H(ns+5,k2+2)-HH
                      HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                      H(k1+2,ns+6)=H(k1+2,ns+6)+HH; H(ns+6,k1+2)=H(ns+6,k1+2)+HH
                      H(k2+2,ns+6)=H(k2+2,ns+6)-HH; H(ns+6,k2+2)=H(ns+6,k2+2)-HH
                
                      HH=-VV*r(3)*r(1)*r(1)
                      H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                      H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                      HH=-VV*r(3)*r(2)*r(2)
                      H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                      H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                      HH=-(V*two*r(3)+VV*r(3)*r(3)*r(3))
                      H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                      H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                      HH=-(V*r(2)+VV*r(3)*r(2)*r(3))
                      H(k1+3,ns+4)=H(k1+3,ns+4)+HH; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                      H(k2+3,ns+4)=H(k2+3,ns+4)-HH; H(ns+4,k2+3)=H(ns+4,k2+3)-HH
                      HH=-(V*r(1)+VV*r(3)*r(1)*r(3))
                      H(k1+3,ns+5)=H(k1+3,ns+5)+HH; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                      H(k2+3,ns+5)=H(k2+3,ns+5)-HH; H(ns+5,k2+3)=H(ns+5,k2+3)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+3,ns+6)=H(k1+3,ns+6)+HH; H(ns+6,k1+3)=H(ns+6,k1+3)+HH
                      H(k2+3,ns+6)=H(k2+3,ns+6)-HH; H(ns+6,k2+3)=H(ns+6,k2+3)-HH
                   else if(slab_calc) then
                      HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                      H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                      H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                      HH=-VV*r(1)*r(2)*r(2)
                      H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                      H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                      HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                      H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                      H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH

                      HH=-VV*r(2)*r(1)*r(1)
                      H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                      H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                      HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                      H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                      H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                      HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                      H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                      H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH

                      HH=-VV*r(3)*r(1)*r(1)
                      H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                      H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                      HH=-VV*r(3)*r(2)*r(2)
                      H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                      H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                      HH=-VV*r(1)*r(2)*r(3)
                      H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                      H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                   end if
                end if

             end do

             !subsraction of interactions between non interacted atoms 
             length1=size(bonded_spcs(i)%list,1)
             if(length1==0) cycle
             do k=1,length1
                kk=bonded_spcs(i)%list(k,1)
                if(kk == i) cycle
                if(kk /= jj) cycle

!!$                r2=red_fac(2)*(atoms_cart(jj)%r-atoms_cart(jj1)%r)+atoms_cart(jj1)%r
                r2=red_fac(2)*(rj-rj1)+rj1
                r=r2-r1
                if(lattice_calc) then 
                   call image(r)
                else if(slab_calc) then
                   call image_slab(r)
                end if
                dr=sqrt(dot_product(r,r))

                dr3=dr*dr*dr
                dr6=dr3*dr3
                dr7=dr6*dr
                dr8=dr7*dr
                dr12=dr6*dr6
                dr13=dr12*dr
                dr14=dr13*dr

                ind=non_bonded(i)%list(2,j)
                id=poten(ind)%id
                p=poten(ind)%param
                select case(id)
                case (4) !buck
                   E_buf=p(1)*exp(-dr/p(2))-p(3)/dr6
                   E(4)=E(4)-E_buf
!!$print*,E_buf*j2c,E(4)*j2c,'=============='
                   if(calc_gradients) dE_dr=-(p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7
                   if(calc_hessian) d2E_dr2=(p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8
                case (5) !l_j
                   E_buf=p(1)/dr12-p(2)/dr6
                   E(4)=E(4)-E_buf
                   if(calc_gradients) dE_dr=-twelve*p(1)/dr13+six*p(2)/dr7
                   if(calc_hessian) d2E_dr2=156.0_r8_kind*p(1)/dr14-42.0_r8_kind*p(3)/dr8
                case (n_poten+1 :) !user_def
                   iup=id-n_poten
                   st=user_poten(iup)%step
                   st2=st*st
                   in=int((dr-user_poten(iup)%r(1))/st)+1
                   eps=(dr-user_poten(iup)%r(in))/st
                   eps2=eps*eps
                   eps3=eps2*eps

                   y1=user_poten(iup)%f(in)
                   y2=user_poten(iup)%f(in+1)
                   y1_2=user_poten(iup)%d2f(in)
                   y2_2=user_poten(iup)%d2f(in+1)
                   E_r=p(1)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                        eps3*st2*(y2_2-y1_2)/six)
                   if(calc_gradients) dE_r_dr=p(1)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                        two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                   if(calc_hessian) d2E_r_dr2= &
                        p(1)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                   y1=user_poten(iup)%f1(in)
                   y2=user_poten(iup)%f1(in+1)
                   y1_2=user_poten(iup)%d2f1(in)
                   y2_2=user_poten(iup)%d2f1(in+1)
                   E_d=p(2)*(y1+eps*(y2-y1-st2*(two*y1_2+y2_2)/six)+eps2*st2*y1_2/two+ &
                        eps3*st2*(y2_2-y1_2)/six)
                   if(calc_gradients) dE_d_dr=p(2)*(y2-y1-st2*(two*y1_2+y2_2)/six+ &
                        two*eps*st2*y1_2/two+three*eps2*st2*(y2_2-y1_2)/six)/st
                   if(calc_hessian) d2E_d_dr2= &
                        p(2)*(two*st2*y1_2/two+six*eps*st2*(y2_2-y1_2)/six)/st2

                   E_buf=E_r+E_d
                   E(id)=E(id)-E_buf
                   if(calc_gradients) dE_dr=dE_r_dr+dE_d_dr
                   if(calc_hessian) d2E_dr2=d2E_r_dr2+d2E_d_dr2
                end select
                E_total=E_total-E_buf

                if(calc_gradients) then
                   f1=r/dr
                   f2=-f1
                   Grad(:,i)=Grad(:,i)-dE_dr*f2*(red_fac(1))
                   if(jj <= n_species) Grad(:,jj)=Grad(:,jj)-dE_dr*f1*(red_fac(2))
                   Grad(:,i1)=Grad(:,i1)-dE_dr*f2*(-red_fac(1)+one)
                   if(jj1 <= n_species) Grad(:,jj1)=Grad(:,jj1)-dE_dr*f1*(-red_fac(2)+one)
                   if(calc_strain) then
                      V=dE_dr/dr
                      if(lattice_calc) then
                         Grad_s(1)=Grad_s(1)-V*r(1)*r(1)
                         Grad_s(2)=Grad_s(2)-V*r(2)*r(2)
                         Grad_s(3)=Grad_s(3)-V*r(3)*r(3)
                         Grad_s(4)=Grad_s(4)-V*r(2)*r(3)
                         Grad_s(5)=Grad_s(5)-V*r(1)*r(3)
                         Grad_s(6)=Grad_s(6)-V*r(1)*r(2)
                      else if(slab_calc) then
                         Grad_s(1)=Grad_s(1)-V*r(1)*r(1)
                         Grad_s(2)=Grad_s(2)-V*r(2)*r(2)
                         Grad_s(3)=Grad_s(3)-V*r(1)*r(2)
                      end if
                   end if
                end if
                
                if(calc_hessian) then
                   a1=red_fac(1); a2=red_fac(2)

                   ff(1,1)=a1*a1/dr-a1*a1*r(1)*r(1)/dr3
                   ff(1,2)=-a1*a1*r(1)*r(2)/dr3
                   ff(1,3)=-a1*a1*r(1)*r(3)/dr3
                   ff(1,4)=-a1*a2/dr+a1*a2*r(1)*r(1)/dr3
                   ff(1,5)=a1*a2*r(1)*r(2)/dr3
                   ff(1,6)=a1*a2*r(1)*r(3)/dr3
                   ff(1,7)=-a1*(a1-one)/dr+a1*(a1-one)*r(1)*r(1)/dr3
                   ff(1,8)=a1*(a1-one)*r(1)*r(2)/dr3
                   ff(1,9)=a1*(a1-one)*r(1)*r(3)/dr3
                   ff(1,10)=-a1*(one-a2)/dr+a1*(one-a2)*r(1)*r(1)/dr3
                   ff(1,11)=a1*(one-a2)*r(1)*r(2)/dr3
                   ff(1,12)=a1*(one-a2)*r(1)*r(3)/dr3

                   ff(2,1)=-a1*a1*r(2)*r(1)/dr3
                   ff(2,2)=a1*a1/dr-a1*a1*r(2)*r(2)/dr3
                   ff(2,3)=-a1*a1*r(2)*r(3)/dr3
                   ff(2,4)=a1*a2*r(2)*r(1)/dr3
                   ff(2,5)=-a1*a2/dr+a1*a2*r(2)*r(2)/dr3
                   ff(2,6)=a1*a2*r(2)*r(3)/dr3
                   ff(2,7)=a1*(a1-one)*r(2)*r(1)/dr3
                   ff(2,8)=-a1*(a1-one)/dr+a1*(a1-one)*r(2)*r(2)/dr3
                   ff(2,9)=a1*(a1-one)*r(2)*r(3)/dr3
                   ff(2,10)=a1*(one-a2)*r(2)*r(1)/dr3
                   ff(2,11)=-a1*(one-a2)/dr+a1*(one-a2)*r(2)*r(2)/dr3
                   ff(2,12)=a1*(one-a2)*r(2)*r(3)/dr3

                   ff(3,1)=-a1*a1*r(3)*r(1)/dr3
                   ff(3,2)=-a1*a1*r(3)*r(2)/dr3
                   ff(3,3)=a1*a1/dr-a1*a1*r(3)*r(3)/dr3
                   ff(3,4)=a1*a2*r(3)*r(1)/dr3
                   ff(3,5)=a1*a2*r(3)*r(2)/dr3
                   ff(3,6)=-a1*a2/dr+a1*a2*r(3)*r(3)/dr3
                   ff(3,7)=a1*(a1-one)*r(3)*r(1)/dr3
                   ff(3,8)=a1*(a1-one)*r(3)*r(2)/dr3
                   ff(3,9)=-a1*(a1-one)/dr+a1*(a1-one)*r(3)*r(3)/dr3
                   ff(3,10)=a1*(one-a2)*r(3)*r(1)/dr3
                   ff(3,11)=a1*(one-a2)*r(3)*r(2)/dr3
                   ff(3,12)=-a1*(one-a2)/dr+a1*(one-a2)*r(3)*r(3)/dr3

                   ff(4,1:3)=ff(1:3,4)
                   ff(4,4)=a2*a2/dr-a2*a2*r(1)*r(1)/dr3
                   ff(4,5)=-a2*a2*r(1)*r(2)/dr3
                   ff(4,6)=-a2*a2*r(1)*r(3)/dr3
                   ff(4,7)=a2*(a1-one)/dr-a2*(a1-one)*r(1)*r(1)/dr3
                   ff(4,8)=-a2*(a1-one)*r(1)*r(2)/dr3
                   ff(4,9)=-a2*(a1-one)*r(1)*r(3)/dr3
                   ff(4,10)=a2*(one-a2)/dr-a2*(one-a2)*r(1)*r(1)/dr3
                   ff(4,11)=-a2*(one-a2)*r(1)*r(2)/dr3
                   ff(4,12)=-a2*(one-a2)*r(1)*r(3)/dr3

                   ff(5,1:3)=ff(1:3,5)
                   ff(5,4)=-a2*a2*r(2)*r(1)/dr3
                   ff(5,5)=a2*a2/dr-a2*a2*r(2)*r(2)/dr3
                   ff(5,6)=-a2*a2*r(2)*r(3)/dr3
                   ff(5,7)=-a2*(a1-one)*r(2)*r(1)/dr3
                   ff(5,8)=a2*(a1-one)/dr-a2*(a1-one)*r(2)*r(2)/dr3
                   ff(5,9)=-a2*(a1-one)*r(2)*r(3)/dr3
                   ff(5,10)=-a2*(one-a2)*r(2)*r(1)/dr3
                   ff(5,11)=a2*(one-a2)/dr-a2*(one-a2)*r(2)*r(2)/dr3
                   ff(5,12)=-a2*(one-a2)*r(2)*r(3)/dr3

                   ff(6,1:3)=ff(1:3,6)
                   ff(6,4)=-a2*a2*r(3)*r(1)/dr3
                   ff(6,5)=-a2*a2*r(3)*r(2)/dr3
                   ff(6,6)=a2*a2/dr-a2*a2*r(3)*r(3)/dr3
                   ff(6,7)=-a2*(a1-one)*r(3)*r(1)/dr3
                   ff(6,8)=-a2*(a1-one)*r(3)*r(2)/dr3
                   ff(6,9)=a2*(a1-one)/dr-a2*(a1-one)*r(3)*r(3)/dr3
                   ff(6,10)=-a2*(one-a2)*r(3)*r(1)/dr3
                   ff(6,11)=-a2*(one-a2)*r(3)*r(2)/dr3
                   ff(6,12)=a2*(one-a2)/dr-a2*(one-a2)*r(3)*r(3)/dr3

                   ff(7,1:3)=ff(1:3,7)
                   ff(7,4:6)=ff(4:6,7)
                   ff(7,7)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(1)*r(1)/dr3
                   ff(7,8)=-(a1-one)*(a1-one)*r(1)*r(2)/dr3
                   ff(7,9)=-(a1-one)*(a1-one)*r(1)*r(3)/dr3
                   ff(7,10)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(1)*r(1)/dr3
                   ff(7,11)=-(a1-one)*(one-a2)*r(1)*r(2)/dr3
                   ff(7,12)=-(a1-one)*(one-a2)*r(1)*r(3)/dr3
                
                   ff(8,1:3)=ff(1:3,8)
                   ff(8,4:6)=ff(4:6,8)
                   ff(8,7)=-(a1-one)*(a1-one)*r(2)*r(1)/dr3
                   ff(8,8)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(2)*r(2)/dr3
                   ff(8,9)=-(a1-one)*(a1-one)*r(2)*r(3)/dr3
                   ff(8,10)=-(a1-one)*(one-a2)*r(2)*r(1)/dr3
                   ff(8,11)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(2)*r(2)/dr3
                   ff(8,12)=-(a1-one)*(one-a2)*r(2)*r(3)/dr3
                   
                   ff(9,1:3)=ff(1:3,9)
                   ff(9,4:6)=ff(4:6,9)
                   ff(9,7)=-(a1-one)*(a1-one)*r(3)*r(1)/dr3
                   ff(9,8)=-(a1-one)*(a1-one)*r(3)*r(2)/dr3
                   ff(9,9)=(a1-one)*(a1-one)/dr-(a1-one)*(a1-one)*r(3)*r(3)/dr3
                   ff(9,10)=-(a1-one)*(one-a2)*r(3)*r(1)/dr3
                   ff(9,11)=-(a1-one)*(one-a2)*r(3)*r(2)/dr3
                   ff(9,12)=(a1-one)*(one-a2)/dr-(a1-one)*(one-a2)*r(3)*r(3)/dr3

                   ff(10,1:3)=ff(1:3,10)
                   ff(10,4:6)=ff(4:6,10)
                   ff(10,7:9)=ff(7:9,10)
                   ff(10,10)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(1)*r(1)/dr3
                   ff(10,11)=-(one-a2)*(one-a2)*r(1)*r(2)/dr3
                   ff(10,12)=-(one-a2)*(one-a2)*r(1)*r(3)/dr3

                   ff(11,1:3)=ff(1:3,11)
                   ff(11,4:6)=ff(4:6,11)
                   ff(11,7:9)=ff(7:9,11)
                   ff(11,10)=-(one-a2)*(one-a2)*r(2)*r(1)/dr3
                   ff(11,11)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(2)*r(2)/dr3
                   ff(11,12)=-(one-a2)*(one-a2)*r(2)*r(3)/dr3

                   ff(12,1:3)=ff(1:3,12)
                   ff(12,4:6)=ff(4:6,12)
                   ff(12,7:9)=ff(7:9,12)
                   ff(12,10)=-(one-a2)*(one-a2)*r(3)*r(1)/dr3
                   ff(12,11)=-(one-a2)*(one-a2)*r(3)*r(2)/dr3
                   ff(12,12)=(one-a2)*(one-a2)/dr-(one-a2)*(one-a2)*r(3)*r(3)/dr3

                   k1=3*(i-1)
                   k2=3*(jj-1)
                   k3=3*(i1-1)
                   k4=3*(jj1-1)
                   do l=1,3
                      do m=1,12
                         if(m <= 3) then
                            fstore=f2(m)*a1
                            m1=k1+m
                         else if(m > 3 .and. m <= 6) then
                            m1=k2+(m-3)
                            fstore=f1(m-3)*a2
                         else if(m > 6 .and. m <= 9) then
                            m1=k3+(m-6)
                            fstore=f2(m-6)*(one-a1)
                         else if(m > 9 .and. m <= 12) then
                            m1=k4+(m-9)
                            fstore=f1(m-9)*(one-a2)
                         end if
                         l1=k1+l; l2=k2+l; l3=k3+l; l4=k4+l
                         if(l1<=ns .and. m1<=ns) H(l1,m1)=H(l1,m1)-d2E_dr2*fstore*f2(l)*a1-dE_dr*ff(l,m)
                         if(l2<=ns .and. m1<=ns) H(l2,m1)=H(l2,m1)-d2E_dr2*fstore*f1(l)*a2-dE_dr*ff(l+3,m)
                         if(l3<=ns .and. m1<=ns) H(l3,m1)=H(l3,m1)-d2E_dr2*fstore*f2(l)*(one-a1)-dE_dr*ff(l+6,m)
                         if(l4<=ns .and. m1<=ns) H(l4,m1)=H(l4,m1)-d2E_dr2*fstore*f1(l)*(one-a2)-dE_dr*ff(l+9,m)
                      end do
                   end do

                   if(calc_strain) then
                      VV=(d2E_dr2-V)/(dr*dr)

                      !d2E_dek_del
                      if(lattice_calc) then
                         H(ns+1,ns+1)=H(ns+1,ns+1)-V*two*r(1)*r(1)-VV*r(1)*r(1)*r(1)*r(1)
                         HH=VV*r(1)*r(1)*r(2)*r(2) 
                         H(ns+1,ns+2)=H(ns+1,ns+2)-HH; H(ns+2,ns+1)=H(ns+2,ns+1)-HH
                         HH=VV*r(1)*r(1)*r(3)*r(3)
                         H(ns+1,ns+3)=H(ns+1,ns+3)-HH; H(ns+3,ns+1)=H(ns+3,ns+1)-HH
                         HH=VV*r(1)*r(1)*r(2)*r(3)
                         H(ns+1,ns+4)=H(ns+1,ns+4)-HH; H(ns+4,ns+1)=H(ns+4,ns+1)-HH
                         HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(1)*r(1)*r(1)*r(3)
                         H(ns+1,ns+5)=H(ns+1,ns+5)-HH; H(ns+5,ns+1)=H(ns+5,ns+1)-HH
                         HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                         H(ns+1,ns+6)=H(ns+1,ns+6)-HH; H(ns+6,ns+1)=H(ns+6,ns+1)-HH

                         H(ns+2,ns+2)=H(ns+2,ns+2)-V*two*r(2)*r(2)-VV*r(2)*r(2)*r(2)*r(2)
                         HH=VV*r(2)*r(2)*r(3)*r(3)
                         H(ns+2,ns+3)=H(ns+2,ns+3)-HH; H(ns+3,ns+2)=H(ns+3,ns+2)-HH
                         HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(2)*r(2)*r(2)*r(3)
                         H(ns+2,ns+4)=H(ns+2,ns+4)-HH; H(ns+4,ns+2)=H(ns+4,ns+2)-HH
                         HH=VV*r(2)*r(2)*r(1)*r(3)
                         H(ns+2,ns+5)=H(ns+2,ns+5)-HH; H(ns+5,ns+2)=H(ns+5,ns+2)-HH
                         HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                         H(ns+2,ns+6)=H(ns+2,ns+6)-HH; H(ns+6,ns+2)=H(ns+6,ns+2)-HH

                         H(ns+3,ns+3)=H(ns+3,ns+3)+V*two*r(3)*r(3)+VV*r(3)*r(3)*r(3)*r(3)
                         HH=V*half*r(2)*r(3)+V*half*r(2)*r(3)+VV*r(3)*r(3)*r(3)*r(2)
                         H(ns+3,ns+4)=H(ns+3,ns+4)-HH; H(ns+4,ns+3)=H(ns+4,ns+3)-HH
                         HH=V*half*r(1)*r(3)+V*half*r(1)*r(3)+VV*r(3)*r(3)*r(3)*r(1)
                         H(ns+3,ns+5)=H(ns+3,ns+5)-HH; H(ns+5,ns+3)=H(ns+5,ns+3)-HH
                         HH=VV*r(3)*r(3)*r(1)*r(2)
                         H(ns+3,ns+6)=H(ns+3,ns+6)-HH; H(ns+6,ns+3)=H(ns+6,ns+3)-HH

                         H(ns+4,ns+4)=H(ns+4,ns+4)- &
                              V*quarter*(r(2)*r(2)+r(3)*r(3))-V*quarter*(r(2)*r(2)+r(3)*r(3))- &
                              VV*r(2)*r(3)*r(2)*r(3)
                         HH=V*quarter*r(1)*r(2)+V*quarter*r(1)*r(2)+VV*r(2)*r(3)*r(3)*r(1)
                         H(ns+4,ns+5)=H(ns+4,ns+5)-HH; H(ns+5,ns+4)=H(ns+5,ns+4)-HH
                         HH=V*quarter*r(1)*r(3)+V*quarter*r(1)*r(3)+VV*r(2)*r(3)*r(2)*r(1)
                         H(ns+4,ns+6)=H(ns+4,ns+6)-HH; H(ns+6,ns+4)=H(ns+6,ns+4)-HH
                   
                         H(ns+5,ns+5)=H(ns+5,ns+5)- &
                              V*quarter*(r(1)*r(1)+r(3)*r(3))-V*quarter*(r(1)*r(1)+r(3)*r(3))- &
                              VV*r(1)*r(3)*r(1)*r(3)
                         HH=V*quarter*r(2)*r(3)+V*quarter*r(2)*r(3)+VV*r(1)*r(2)*r(1)*r(3)
                         H(ns+5,ns+6)=H(ns+5,ns+6)-HH; H(ns+6,ns+5)=H(ns+6,ns+5)-HH
                         
                         H(ns+6,ns+6)=H(ns+6,ns+6)- &
                              V*quarter*(r(1)*r(1)+r(2)*r(2))-V*quarter*(r(1)*r(1)+r(2)*r(2))- &
                              VV*r(1)*r(2)*r(1)*r(2)
                      else if(slab_calc) then
                         H(ns+1,ns+1)=H(ns+1,ns+1)-V*two*r(1)*r(1)-VV*r(1)*r(1)*r(1)*r(1)
                         HH=VV*r(1)*r(1)*r(2)*r(2) 
                         H(ns+1,ns+2)=H(ns+1,ns+2)-HH; H(ns+2,ns+1)=H(ns+2,ns+1)-HH
                         HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(1)*r(1)*r(1)*r(2)
                         H(ns+1,ns+3)=H(ns+1,ns+3)-HH; H(ns+3,ns+1)=H(ns+3,ns+1)-HH

                         H(ns+2,ns+2)=H(ns+2,ns+2)-V*two*r(2)*r(2)-VV*r(2)*r(2)*r(2)*r(2)
                         HH=V*half*r(1)*r(2)+V*half*r(1)*r(2)+VV*r(2)*r(2)*r(2)*r(1)
                         H(ns+2,ns+3)=H(ns+2,ns+3)-HH; H(ns+3,ns+2)=H(ns+3,ns+2)-HH

                         H(ns+3,ns+3)=H(ns+3,ns+3)- &
                              V*quarter*(r(1)*r(1)+r(2)*r(2))-V*quarter*(r(1)*r(1)+r(2)*r(2))- &
                              VV*r(1)*r(2)*r(1)*r(2)
                      end if

                      !d2E_dri_dek
                      k1=3*(i-1); k2=3*(jj-1)

                      if(lattice_calc) then
                         HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                         H(k1+1,ns+1)=H(k1+1,ns+1)-HH; H(ns+1,k1+1)=H(ns+1,k1+1)-HH
                         H(k2+1,ns+1)=H(k2+1,ns+1)+HH; H(ns+1,k2+1)=H(ns+1,k2+1)+HH
                         HH=-VV*r(1)*r(2)*r(2)
                         H(k1+1,ns+2)=H(k1+1,ns+2)-HH; H(ns+2,k1+1)=H(ns+2,k1+1)-HH
                         H(k2+1,ns+2)=H(k2+1,ns+2)+HH; H(ns+2,k2+1)=H(ns+2,k2+1)+HH
                         HH=-VV*r(1)*r(3)*r(3)
                         H(k1+1,ns+3)=H(k1+1,ns+3)-HH; H(ns+3,k1+1)=H(ns+3,k1+1)-HH
                         H(k2+1,ns+3)=H(k2+1,ns+3)+HH; H(ns+3,k2+1)=H(ns+3,k2+1)+HH
                         HH=-VV*r(1)*r(2)*r(3)
                         H(k1+1,ns+4)=H(k1+1,ns+4)-HH; H(ns+4,k1+1)=H(ns+4,k1+1)-HH
                         H(k2+1,ns+4)=H(k2+1,ns+4)+HH; H(ns+4,k2+1)=H(ns+4,k2+1)+HH
                         HH=-(V*r(3)+VV*r(1)*r(1)*r(3))
                         H(k1+1,ns+5)=H(k1+1,ns+5)-HH; H(ns+5,k1+1)=H(ns+5,k1+1)-HH
                         H(k2+1,ns+5)=H(k2+1,ns+5)+HH; H(ns+5,k2+1)=H(ns+5,k2+1)+HH
                         HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                         H(k1+1,ns+6)=H(k1+1,ns+6)-HH; H(ns+6,k1+1)=H(ns+6,k1+1)-HH
                         H(k2+1,ns+6)=H(k2+1,ns+6)+HH; H(ns+6,k2+1)=H(ns+6,k2+1)+HH
                
                         HH=-VV*r(2)*r(1)*r(1)
                         H(k1+2,ns+1)=H(k1+2,ns+1)-HH; H(ns+1,k1+2)=H(ns+1,k1+2)-HH
                         H(k2+2,ns+1)=H(k2+2,ns+1)+HH; H(ns+1,k2+2)=H(ns+1,k2+2)+HH
                         HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                         H(k1+2,ns+2)=H(k1+2,ns+2)-HH; H(ns+2,k1+2)=H(ns+2,k1+2)-HH
                         H(k2+2,ns+2)=H(k2+2,ns+2)+HH; H(ns+2,k2+2)=H(ns+2,k2+2)+HH
                         HH=-VV*r(2)*r(3)*r(3)
                         H(k1+2,ns+3)=H(k1+2,ns+3)-HH; H(ns+3,k1+2)=H(ns+3,k1+2)-HH
                         H(k2+2,ns+3)=H(k2+2,ns+3)+HH; H(ns+3,k2+2)=H(ns+3,k2+2)+HH
                         HH=-(V*r(3)+VV*r(2)*r(2)*r(3))
                         H(k1+2,ns+4)=H(k1+2,ns+4)-HH; H(ns+4,k1+2)=H(ns+4,k1+2)-HH
                         H(k2+2,ns+4)=H(k2+2,ns+4)+HH; H(ns+4,k2+2)=H(ns+4,k2+2)+HH
                         HH=-VV*r(1)*r(2)*r(3)
                         H(k1+2,ns+5)=H(k1+2,ns+5)-HH; H(ns+5,k1+2)=H(ns+5,k1+2)-HH
                         H(k2+2,ns+5)=H(k2+2,ns+5)+HH; H(ns+5,k2+2)=H(ns+5,k2+2)-HH
                         HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                         H(k1+2,ns+6)=H(k1+2,ns+6)-HH; H(ns+6,k1+2)=H(ns+6,k1+2)-HH
                         H(k2+2,ns+6)=H(k2+2,ns+6)+HH; H(ns+6,k2+2)=H(ns+6,k2+2)+HH
                
                         HH=-VV*r(3)*r(1)*r(1)
                         H(k1+3,ns+1)=H(k1+3,ns+1)-HH; H(ns+1,k1+3)=H(ns+1,k1+3)-HH
                         H(k2+3,ns+1)=H(k2+3,ns+1)+HH; H(ns+1,k2+3)=H(ns+1,k2+3)+HH
                         HH=-VV*r(3)*r(2)*r(2)
                         H(k1+3,ns+2)=H(k1+3,ns+2)-HH; H(ns+2,k1+3)=H(ns+2,k1+3)-HH
                         H(k2+3,ns+2)=H(k2+3,ns+2)+HH; H(ns+2,k2+3)=H(ns+2,k2+3)+HH
                         HH=-(V*two*r(3)+VV*r(3)*r(3)*r(3))
                         H(k1+3,ns+3)=H(k1+3,ns+3)-HH; H(ns+3,k1+3)=H(ns+3,k1+3)-HH
                         H(k2+3,ns+3)=H(k2+3,ns+3)+HH; H(ns+3,k2+3)=H(ns+3,k2+3)+HH
                         HH=-(V*r(2)+VV*r(3)*r(2)*r(3))
                         H(k1+3,ns+4)=H(k1+3,ns+4)-HH; H(ns+4,k1+3)=H(ns+4,k1+3)-HH
                         H(k2+3,ns+4)=H(k2+3,ns+4)+HH; H(ns+4,k2+3)=H(ns+4,k2+3)+HH
                         HH=-(V*r(1)+VV*r(3)*r(1)*r(3))
                         H(k1+3,ns+5)=H(k1+3,ns+5)-HH; H(ns+5,k1+3)=H(ns+5,k1+3)-HH
                         H(k2+3,ns+5)=H(k2+3,ns+5)+HH; H(ns+5,k2+3)=H(ns+5,k2+3)+HH
                         HH=-VV*r(1)*r(2)*r(3)
                         H(k1+3,ns+6)=H(k1+3,ns+6)-HH; H(ns+6,k1+3)=H(ns+6,k1+3)-HH
                         H(k2+3,ns+6)=H(k2+3,ns+6)+HH; H(ns+6,k2+3)=H(ns+6,k2+3)+HH
                      else if(slab_calc) then
                         HH=-(V*two*r(1)+VV*r(1)*r(1)*r(1))
                         H(k1+1,ns+1)=H(k1+1,ns+1)-HH; H(ns+1,k1+1)=H(ns+1,k1+1)-HH
                         H(k2+1,ns+1)=H(k2+1,ns+1)+HH; H(ns+1,k2+1)=H(ns+1,k2+1)+HH
                         HH=-VV*r(1)*r(2)*r(2)
                         H(k1+1,ns+2)=H(k1+1,ns+2)-HH; H(ns+2,k1+1)=H(ns+2,k1+1)-HH
                         H(k2+1,ns+2)=H(k2+1,ns+2)+HH; H(ns+2,k2+1)=H(ns+2,k2+1)+HH
                         HH=-(V*r(2)+VV*r(1)*r(1)*r(2))
                         H(k1+1,ns+3)=H(k1+1,ns+3)-HH; H(ns+3,k1+1)=H(ns+3,k1+1)-HH
                         H(k2+1,ns+3)=H(k2+1,ns+3)+HH; H(ns+3,k2+1)=H(ns+3,k2+1)+HH

                         HH=-VV*r(2)*r(1)*r(1)
                         H(k1+2,ns+1)=H(k1+2,ns+1)-HH; H(ns+1,k1+2)=H(ns+1,k1+2)-HH
                         H(k2+2,ns+1)=H(k2+2,ns+1)+HH; H(ns+1,k2+2)=H(ns+1,k2+2)+HH
                         HH=-(V*two*r(2)+VV*r(2)*r(2)*r(2))
                         H(k1+2,ns+2)=H(k1+2,ns+2)-HH; H(ns+2,k1+2)=H(ns+2,k1+2)-HH
                         H(k2+2,ns+2)=H(k2+2,ns+2)+HH; H(ns+2,k2+2)=H(ns+2,k2+2)+HH
                         HH=-(V*r(1)+VV*r(2)*r(1)*r(2))
                         H(k1+2,ns+3)=H(k1+2,ns+3)-HH; H(ns+3,k1+2)=H(ns+3,k1+2)-HH
                         H(k2+2,ns+3)=H(k2+2,ns+3)+HH; H(ns+3,k2+2)=H(ns+3,k2+2)+HH

                         HH=-VV*r(3)*r(1)*r(1)
                         H(k1+3,ns+1)=H(k1+3,ns+1)-HH; H(ns+1,k1+3)=H(ns+1,k1+3)-HH
                         H(k2+3,ns+1)=H(k2+3,ns+1)+HH; H(ns+1,k2+3)=H(ns+1,k2+3)+HH
                         HH=-VV*r(3)*r(2)*r(2)
                         H(k1+3,ns+2)=H(k1+3,ns+2)-HH; H(ns+2,k1+3)=H(ns+2,k1+3)-HH
                         H(k2+3,ns+2)=H(k2+3,ns+2)+HH; H(ns+2,k2+3)=H(ns+2,k2+3)+HH
                         HH=-VV*r(1)*r(2)*r(3)
                         H(k1+3,ns+3)=H(k1+3,ns+3)-HH; H(ns+3,k1+3)=H(ns+3,k1+3)-HH
                         H(k2+3,ns+3)=H(k2+3,ns+3)+HH; H(ns+3,k2+3)=H(ns+3,k2+3)+HH
                      end if
                   end if
                end if
             end do
          end if
       end do
    end do
!!$print*,"E_total=",E_total
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine vdw_E_and_F
  !****************************************************************

  !****************************************************************
end module van_der_waals_module

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
module ewald2d_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module
  use species_module
  use slab_module
  use energy_and_forces_module
  use n_body_lists_module
  use potentials_module, only : n_cs
  use coulomb_module
  use ewald_module, only : erf_erfc
  use molmech_msgtag_module
  use comm_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ------------------
  public init_ewald2d, self_ewald2d_energy,direct_ew2d_E_and_F, &
       recipr_ew2d_E_and_F, recipr_ew2d_E_self, &
       calc_2d_slab_vec, shutdown_ewald2d, send_receive_ewald2d, &
       calc_2d_poten
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  real(r8_kind) :: accuracy=small
  real(r8_kind) :: ewald_param ! alpha
  real(r8_kind) :: ewp3
  real(r8_kind) :: kcut_ew
  real(r8_kind) :: sqrt_pi
  integer(i4_kind), parameter :: max_kv = 200 

  integer(i4_kind), allocatable :: k_indexes(:)
  integer(i4_kind) :: n_kvec
  real(r8_kind), allocatable :: kvector(:,:),factor(:),dk(:)

  integer(i4_kind) :: first_index,last_index

  real(r8_kind) :: E_self,E_self0,E_self1
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine init_ewald2d()

    real(r8_kind) :: f

    ewald_param=1.2_r8_kind/sqrt(area)
    ewp3=ewald_param*ewald_param*ewald_param
    sqrt_pi=sqrt(pi)
    f=sqrt(-log(accuracy))
    rcut_ew=f/ewald_param
    kcut_ew=two*f*ewald_param

    call calc_kindexes()

  end subroutine init_ewald2d
  !****************************************************************

  !****************************************************************
  subroutine calc_kindexes()

    integer(i4_kind) :: max1,max2
    real(r8_kind) :: drk,rk1(2),rk2(2),rk(2)
    integer(i4_kind) :: i,j,status

    type store
       integer(i4_kind) :: index_store
       type(store), pointer :: next_data
    end type store
    type(store), target :: first_data
    type(store), pointer ::   current_data, tmp_data

    rk1=kvect_s%v1
    drk=sqrt(dot_product(rk1,rk1))
    max1=int(kcut_ew/drk)+1

    rk2=kvect_s%v2
    drk=sqrt(dot_product(rk2,rk2))
    max2=int(kcut_ew/drk)+1

    if(max1 > max_kv .or. max2 > max_kv) &
         call error_handler("MolMech: Number of K vectors is too large(2d)")

    current_data=>first_data
    nullify(current_data%next_data)

    n_kvec=0
    do i=0,max1
       rk1=i*kvect_s%v1
       do j=-max2,max2
          if(i==0.and.j==0) cycle
          rk2=j*kvect_s%v2
          rk=rk1+rk2
          drk=sqrt(dot_product(rk,rk))
          if(drk > kcut_ew) cycle
          n_kvec=n_kvec+1
          allocate(tmp_data,stat=status)
          if(status /= 0) call error_handler("MolMech: failed tmp_data allocation (ew2)")
          tmp_data%index_store=1000*(i+max_kv)+(j+max_kv)
          nullify(tmp_data%next_data)
          current_data%next_data =>tmp_data
          current_data => tmp_data
       end do
    end do

    allocate(k_indexes(n_kvec),stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Failed allocation of K_INDEXES")
    current_data=>first_data
    do i=1,n_kvec
       current_data => current_data%next_data
       k_indexes(i)=current_data%index_store
    end do
    nullify(current_data)

  end subroutine calc_kindexes
  !****************************************************************

  !****************************************************************
  subroutine send_receive_ewald2d()

    integer(i4_kind) :: info,status,i,myindex,n_proc,nn


    n_proc=comm_get_n_processors()
    myindex=comm_myindex()
    nn = nint (real (n_species) / real (n_proc))
    first_index=nn*(myindex-1)+1
    last_index=first_index+nn-1
    if(myindex==n_proc) last_index=n_species

    if(comm_i_am_master()) then
       do i=2,n_proc
          call comm_init_send(i,msgtag_mm_send_ewald2d)

          call commpack(ewald_param,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: ewald_param pack failed(2d)")
          call commpack(ewp3,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: ewp3 pack failed(2d)")
          call commpack(sqrt_pi,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: sqrt_pi pack failed(2d)")
          call commpack(kvect_s%v1(1),2,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v1 pack failed(2d)")
          call commpack(kvect_s%v2(1),2,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v2 pack failed(2d)")
          call commpack(area,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: area pack failed")

          call commpack(n_kvec,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: nkv pack failed(2d)")
          call commpack(k_indexes(1),n_kvec,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: k_indexes pack failed(2d)")

          call comm_send()
       end do
    else
       call communpack(ewald_param,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: ewald_param unpack failed(2d)")
       call communpack(ewp3,info)
       if( info /= 0) call error_handler &
               ("send_receive_ewald: ewp3 unpack failed(2d)")
       call communpack(sqrt_pi,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: sqrt_pi unpack failed(2d)")
       call communpack(kvect_s%v1(1),2,1,info)
       if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v1 unpack failed(2d)")
       call communpack(kvect_s%v2(1),2,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: kvect%v2 unpack failed(2d)")
       call communpack(area,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: area unpack failed")

       call communpack(n_kvec,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: nkv unpack failed(2d)")
       allocate(k_indexes(n_kvec),stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: Failed allocation of K_INDEXES on slave(2d)")
       
       call communpack(k_indexes(1),n_kvec,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: k_indexes unpack failed(2d)")

    end if

  end subroutine send_receive_ewald2d
  !****************************************************************

  !****************************************************************
  subroutine self_ewald2d_energy()

    real(r8_kind) :: self_energy
    integer(i4_kind) :: i
    integer(kind=i4_kind) :: a_type1
    real(r8_kind) :: q1

    self_energy=zero
    do i=1,n_species
       a_type1=atoms_cart(i)%type
       q1=atoms(a_type1)%charge  
       self_energy=self_energy-q1*q1
    end do

    E_self0=(coulomb_factor*ewald_param/(sqrt_pi))*self_energy

  end subroutine  self_ewald2d_energy
  !****************************************************************

  !****************************************************************
  subroutine direct_ew2d_E_and_F()

    integer(i4_kind) :: i,j,jj,k,kk,l,m,k1,k2,m1,length,length1,count
    integer(i4_kind) :: img,i1,j1
    integer(i4_kind) :: a_type1,a_type2,ns
    real(r8_kind) :: q1,q2,r(3),dr,sfac,erfun,erfunc,er
    real(r8_kind) :: dr2,dr3,cqq,er2
    real(r8_kind) :: f1(3),f2(3),ff(6,6),fstore,V,VV
    real(r8_kind) :: E_buf,dE_dr,d2E_dr2,HH
    real(r8_kind) :: self_energy,rsm

    rsm=small0*ten

    !taking into account that core and shell are the same atom
    self_energy=zero
    do i=1,n_cs
       length=size(core_shell(i)%list,2)
       do j=1,length
          i1=core_shell(i)%list(1,j)
          j1=core_shell(i)%list(2,j)
          r=atoms_cart(i1)%r-atoms_cart(j1)%r
          dr=sqrt(dot_product(r,r))
          if(dr < rsm) then
             a_type1=atoms_cart(i1)%type
             a_type2=atoms_cart(j1)%type
             q1=atoms(a_type1)%charge  
             q2=atoms(a_type2)%charge
             self_energy=self_energy-two*q1*q2
          end if
       end do
    end do

    E_self=E_self0+(coulomb_factor*ewald_param/(sqrt_pi))*self_energy

    E_coulomb=E_coulomb+E_self
    E_total=E_total+E_self
    E_ew_r=E_ew_r+E_self

    ns=3*n_species

    do i=1,n_species
       length=size(coulombic(i)%list,1)
       if(length==0) cycle
       a_type1=atoms_cart(i)%type
       q1=atoms(a_type1)%charge

       do j=1,length
          k=coulombic(i)%list(j,1)
          a_type2=atoms_cart(k)%type
          q2=atoms(a_type2)%charge
          cqq=coulomb_factor*q1*q2

          if(coulombic(i)%first_image(j)) then
             !calculation direct contributions to Ewald  between atoms in real unit cell
             r=atoms_cart(k)%r-atoms_cart(i)%r
             if(slab_calc .and. minimal_image) call image_slab(r)
             dr=sqrt(dot_product(r,r))
             if(dr < rsm) goto 1
             dr2=dr*dr; dr3=dr2*dr
             er=ewald_param*dr; er2=er*er
             erfunc=erf_erfc(er,ierfc)

             E_buf=cqq*erfunc/dr
             E_coulomb=E_coulomb+E_buf
             E_total=E_total+E_buf
             E_ew_d=E_ew_d+E_buf

             if(calc_gradients) dE_dr=cqq*(-erfunc/dr2-two*ewald_param*exp(-er2)/(dr*sqrt_pi))

             if(calc_hessian) d2E_dr2=cqq*two*(erfunc/dr3+ &
                  two*ewald_param*exp(-er2)/(dr2*sqrt_pi)+ &
                  two*ewp3*exp(-er2)/sqrt_pi)

             if(calc_gradients) then
                f1=r/dr
                f2=-f1
                Grad(:,i)=Grad(:,i)+dE_dr*f2
                Grad(:,k)=Grad(:,k)+dE_dr*f1
                if(calc_strain) then
                   V=dE_dr/dr
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                end if
             endif

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
             
                k1=3*(i-1)
                k2=3*(k-1)
                do l=1,3
                   do m=1,6
                      fstore=f1(m)
                      m1=k1+m
                      if(m > 3) then 
                         m1=k2+(m-3)
                         fstore=f2(m-3)
                      end if
!                      H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!                      H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                      H(k1+l,m1)=H(k1+l,m1)-d2E_dr2*fstore*f2(l)-dE_dr*ff(l+3,m)  !??????
                      H(k2+l,m1)=H(k2+l,m1)-d2E_dr2*fstore*f1(l)-dE_dr*ff(l,m)    !?????
                   end do
                end do

                if(calc_strain) then
                   VV=(d2E_dr2-V)/(dr*dr)

                   !d2E_dek_del
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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

1         if(.not.minimal_image) then
             !calculation direct contributions to Ewald  between atoms in real unit cell 
             !and unit cell images
             sfac=one; if(i==k) sfac=half
             count=0
             do
                count=count+1
                img=add_list_c(i)%list(count,j)
                if(img==0) exit
                r=im_coor(k)%r(:,img)-atoms_cart(i)%r
                dr=sqrt(dot_product(r,r))
                dr2=dr*dr; dr3=dr2*dr
                er=ewald_param*dr; er2=er*er
                erfunc=erf_erfc(er,ierfc)

                E_buf=sfac*cqq*erfunc/dr
                E_coulomb=E_coulomb+E_buf
                E_total=E_total+E_buf
                E_ew_d=E_ew_d+E_buf

                if(calc_gradients) dE_dr=sfac*cqq*(-erfunc/dr2- &
                     two*ewald_param*exp(-er2)/(dr*sqrt_pi))
                if(calc_hessian) d2E_dr2=sfac*cqq*two*(erfunc/dr3+ &
                     two*ewald_param*exp(-er2)/(dr2*sqrt_pi)+ &
                     two*ewp3*exp(-er2)/sqrt_pi)

                if(calc_gradients) then
                   f1=r/dr
                   f2=-f1
                   Grad(:,i)=Grad(:,i)+dE_dr*f2
                   Grad(:,k)=Grad(:,k)+dE_dr*f1
                   if(calc_strain) then
                      V=dE_dr/dr
                      Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                      Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                      Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                   end if
                endif

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
             
                   k1=3*(i-1)
                   k2=3*(k-1)
                   do l=1,3
                      do m=1,6
                         fstore=f1(m)
                         m1=k1+m
                         if(m > 3) then 
                            m1=k2+(m-3)
                            fstore=f2(m-3)
                         end if
!                         H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!                         H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                         H(k1+l,m1)=H(k1+l,m1)-d2E_dr2*fstore*f2(l)-dE_dr*ff(l+3,m)  !?????
                         H(k2+l,m1)=H(k2+l,m1)-d2E_dr2*fstore*f1(l)-dE_dr*ff(l,m)    !?????
                      end do
                   end do

                   if(calc_strain) then
                      VV=(d2E_dr2-V)/(dr*dr)

                      !d2E_dek_del
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

                      !d2E_dri_dek
                      k1=3*(i-1); k2=3*(k-1)

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
          end if

          length1=size(bonded_spcs_c(i)%list,1)
          if(length1==0) cycle
          do jj=1,length1
             kk=bonded_spcs_c(i)%list(jj,1)
             if(kk == i) cycle
             if(kk /= k) cycle

             !subsraction interactions between non interacted atoms in direct space 
             r=atoms_cart(k)%r-atoms_cart(i)%r
             if(slab_calc) call image_slab(r)
             dr=sqrt(dot_product(r,r))
             if(dr < rsm) cycle
             dr2=dr*dr; dr3=dr2*dr
             er=ewald_param*dr; er2=er*er
             erfunc=erf_erfc(er,ierfc)

             E_buf=-cqq*erfunc/dr
             E_coulomb=E_coulomb+E_buf
             E_total=E_total+E_buf
             E_ew_d=E_ew_d+E_buf

             if(calc_gradients) dE_dr=-cqq*(-erfunc/dr2- &
                  two*ewald_param*exp(-er2)/(dr*sqrt_pi))
             if(calc_hessian) d2E_dr2=-cqq*two*(erfunc/dr3+ &
                  two*ewald_param*exp(-er2)/(dr2*sqrt_pi)+ &
                  two*ewp3*exp(-er2)/sqrt_pi)

             if(calc_gradients) then
                f1=r/dr
                f2=-f1
                Grad(:,i)=Grad(:,i)+dE_dr*f2
                Grad(:,k)=Grad(:,k)+dE_dr*f1
                if(calc_strain) then
                   V=dE_dr/dr
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                end if
             endif

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
             
                k1=3*(i-1)
                k2=3*(k-1)
                do l=1,3
                   do m=1,6
                      fstore=f1(m)
                      m1=k1+m
                      if(m > 3) then 
                         m1=k2+(m-3)
                         fstore=f2(m-3)
                      end if
!                      H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)    
!                      H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)      
                      H(k1+l,m1)=H(k1+l,m1)-d2E_dr2*fstore*f2(l)-dE_dr*ff(l+3,m)  !?????
                      H(k2+l,m1)=H(k2+l,m1)-d2E_dr2*fstore*f1(l)-dE_dr*ff(l,m)   !?????
                   end do
                end do

                if(calc_strain) then
                   VV=(d2E_dr2-V)/(dr*dr)

                   !d2E_dek_del
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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

             !subsraction interactions between non interacted atoms in reciprocal space 
             erfun=erf_erfc(er,ierf)

             E_buf=-cqq*erfun/dr
             E_ew_d=E_ew_d+E_buf
             E_coulomb=E_coulomb+E_buf
             E_total=E_total+E_buf

             if(calc_gradients) dE_dr=cqq*(erfun/dr2- &
                  two*ewald_param*exp(-er2)/(dr*sqrt_pi))

             if(calc_hessian) d2E_dr2=cqq*two*(-erfun/dr3+ &
                  two*ewald_param*exp(-er2)/(dr2*sqrt_pi)+ &
                  two*ewp3*exp(-er2)/sqrt_pi)

             if(calc_gradients) then
                Grad(:,i)=Grad(:,i)+dE_dr*f2
                Grad(:,k)=Grad(:,k)+dE_dr*f1
                if(calc_strain) then
                   V=dE_dr/dr
                   Grad_s(1)=Grad_s(1)+V*r(1)*r(1)
                   Grad_s(2)=Grad_s(2)+V*r(2)*r(2)
                   Grad_s(3)=Grad_s(3)+V*r(1)*r(2)
                end if
             endif

             if(calc_hessian) then
                k1=3*(i-1)
                k2=3*(k-1)
                do l=1,3
                   do m=1,6
                      fstore=f1(m)
                      m1=k1+m
                      if(m > 3) then 
                         m1=k2+(m-3)
                         fstore=f2(m-3)
                      end if
!                      H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!                      H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                      H(k1+l,m1)=H(k1+l,m1)-d2E_dr2*fstore*f2(l)-dE_dr*ff(l+3,m)  !?????
                      H(k2+l,m1)=H(k2+l,m1)-d2E_dr2*fstore*f1(l)-dE_dr*ff(l,m)    !????
                   end do
                end do

                if(calc_strain) then
                   VV=(d2E_dr2-V)/(dr*dr)

                   !d2E_dek_del
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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
       end do
    end do

  end subroutine direct_ew2d_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine calc_2d_slab_vec()

    integer(i4_kind) :: ii,i,j,k_ind,status
    real(r8_kind) :: kv(2)

    allocate(kvector(2,n_kvec),dk(n_kvec),factor(n_kvec),stat=status)
    if(status /= 0) call error_handler("MolMech: failed allocation of KVECTOR for slab")

    do ii=1,n_kvec
       factor(ii)=two
       k_ind=k_indexes(ii)
       i=k_ind/1000-max_kv
       j=(k_ind-(i+max_kv)*1000)-max_kv
       if(i==0) factor(ii)=one
       kv=i*kvect_s%v1+j*kvect_s%v2
       kvector(:,ii)=kv
       dk(ii)=sqrt(dot_product(kv,kv))
    end do

  end subroutine calc_2d_slab_vec
  !****************************************************************

  !****************************************************************
  subroutine recipr_ew2d_E_self()

    integer(i4_kind) :: i,ik,a_type,ns
    real(r8_kind) :: ep2,erfunc,er,cqq,cf,cf1,dk2,HH
    real(r8_kind) :: E_buf_z,E_buf,qi,ffff,kv11,kv12,kv22,iiii

    ns=3*n_species

    ep2=ewald_param*ewald_param
    cf=coulomb_factor*pi/area
    cf1=one/(ewald_param*sqrt_pi)

    E_self1=zero
    do i=first_index,last_index
       a_type=atoms_cart(i)%type
       qi=atoms(a_type)%charge
       cqq=cf*qi*qi

       E_buf_z=-cqq*cf1
       E_self1=E_self1+E_buf_z
       
       if(calc_gradients .and. calc_strain) then
          Grad_s(1)=Grad_s(1)-E_buf_z
          Grad_s(2)=Grad_s(2)-E_buf_z

          if(calc_hessian) then
             !d2E_dek_del
             H(ns+1,ns+1)=H(ns+1,ns+1)+E_buf_z
             H(ns+1,ns+2)=H(ns+1,ns+2)+E_buf_z; H(ns+2,ns+1)=H(ns+2,ns+1)+E_buf_z

             H(ns+2,ns+2)=H(ns+2,ns+2)+E_buf_z
          end if
       end if

       do ik=1,n_kvec
          dk2=dk(ik)*dk(ik)
          er=dk(ik)/(two*ewald_param)
          erfunc=erf_erfc(er,ierfc)

          E_buf=factor(ik)*cqq*erfunc/dk(ik)
          E_self1=E_self1+E_buf

          if(calc_gradients .and. calc_strain) then
             kv11=kvector(1,ik)*kvector(1,ik)
             kv22=kvector(2,ik)*kvector(2,ik)
             kv12=kvector(1,ik)*kvector(2,ik)

             ffff=cf1*exp(-er*er)/(erfunc*dk(ik))+one/dk2

             Grad_s(1)=Grad_s(1)-E_buf*(one-kv11*ffff)
             Grad_s(2)=Grad_s(2)-E_buf*(one-kv22*ffff)
             Grad_s(3)=Grad_s(3)+E_buf*kv12*ffff

             if(calc_hessian) then
                !d2E_dek_del
                iiii=(cf1*exp(-er*er)/((erfunc*dk(ik))))*(one/(four*ep2)- &
                     cf1*exp(-er*er)/(erfunc*dk(ik))+one/dk2)+two/(dk2*dk2)

                H(ns+1,ns+1)=H(ns+1,ns+1)+E_buf*(one-kv11*ffff)*(one-kv11*ffff)- &
                     E_buf*ffff*two*kv11+E_buf*kv11*kv11*iiii
                HH=E_buf*(one-kv11*ffff)*(one-kv22*ffff)+E_buf*kv11*kv22*iiii
                H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                HH=-E_buf*(one-kv11*ffff)*kv12*ffff-E_buf*ffff*kv12+ &
                     E_buf*kv11*kv12*iiii
                H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH

                H(ns+2,ns+2)=H(ns+2,ns+2)+E_buf*(one-kv22*ffff)*(one-kv22*ffff)- &
                     E_buf*ffff*two*kv22+E_buf*kv22*kv22*iiii
                HH=-E_buf*(one-kv22*ffff)*kv12*ffff-E_buf*ffff*kv12+ &
                     E_buf*kv22*kv12*iiii
                H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH

                H(ns+3,ns+3)=H(ns+3,ns+3)+E_buf*kv12*ffff*kv12*ffff- &
                     E_buf*ffff*half*(kv11+kv22)+E_buf*kv12*kv12*iiii
             end if
          end if
       end do
    end do

  end subroutine recipr_ew2d_E_self
  !****************************************************************

  !****************************************************************
  subroutine recipr_ew2d_E_and_F()

    integer(i4_kind) :: i,j,jj,ik,a_type,k1,k2,l,m,ns
    real(r8_kind) :: erfunc_z,ez,ez2,er1,er2,erfunc1,erfunc2
    real(r8_kind) :: ep2,twoep,cf,cf1,cf2,cqq,qi,qj,ri(3),rj(3)
    real(r8_kind) :: z,dz,dz2,r(2),krd,kzd,dk2,kv11,kv12,kv22
    real(r8_kind) :: E_buf_z,E_buf,dE_dz,dE_dr
    real(r8_kind) :: d2E_dz2,d2E_dr2,d2E_dzdr,d2E_drdz,exp12,exp22
    real(r8_kind) :: expez2,expkzd,exp_kzd,coskrd,sinkrd
    real(r8_kind) :: eeee,hhhh,gggg,ffff,ffff1,iiii,jjjj
    real(r8_kind) :: ff(6,6),HH

    ns=3*n_species

    ep2=ewald_param*ewald_param
    twoep=two*ewald_param
    cf=coulomb_factor*pi/area
    cf1=one/(ewald_param*sqrt_pi)

    E_total=E_total+E_self1
    E_coulomb=E_coulomb+E_self1
    E_ew_r=E_ew_r+E_self1

    do i=first_index,last_index
       a_type=atoms_cart(i)%type
       qi=atoms(a_type)%charge

       ri=atoms_cart(i)%r

       do j=1,N_total(i)
          jj=i+j
          if(jj > n_species) jj=jj-n_species
          a_type=atoms_cart(jj)%type
          qj=atoms(a_type)%charge
          cqq=cf*qi*qj

          rj=atoms_cart(jj)%r

          z=rj(3)-ri(3); dz=abs(z); dz2=dz*dz
          ez=ewald_param*z; ez2=ez*ez
          erfunc_z=erf_erfc(ez,ierf)
          expez2=exp(-ez2)

          E_buf_z=-two*cqq*(z*erfunc_z+cf1*expez2)
          E_total=E_total+E_buf_z
          E_coulomb=E_coulomb+E_buf_z
          E_ew_r=E_ew_r+E_buf_z

          if(calc_gradients) dE_dz=-two*cqq*(erfunc_z+ &
               z*twoep*expez2/sqrt_pi- &
               two*cf1*ep2*z*expez2)

          if(calc_hessian) d2E_dz2=-two*cqq*expez2* &
               (twoep/sqrt_pi+twoep/sqrt_pi-four*ep2*dz2/sqrt_pi- &
               two*cf1*ep2+four*cf1*ewp3*dz2)

          if(calc_gradients) then
!!$             f1_z=one
!!$             f2_z=-f1_z
!!$             Grad(3,i)=Grad(3,i)+dE_dz*f2_z
!!$             Grad(3,jj)=Grad(3,jj)+dE_dz*f1_z
             Grad(3,i)=Grad(3,i)-dE_dz
             Grad(3,jj)=Grad(3,jj)+dE_dz
             if(calc_strain) then
                Grad_s(1)=Grad_s(1)-E_buf_z
                Grad_s(2)=Grad_s(2)-E_buf_z
             end if
          end if

          if(calc_hessian) then
             k1=3*(i-1)
             k2=3*(jj-1)
!!$             H(k1+3,k1+3)=H(k1+3,k1+3)+d2E_dz2*f2_z*f2_z
!!$             H(k1+3,k2+3)=H(k1+3,k2+3)+d2E_dz2*f1_z*f2_z
!!$             H(k2+3,k1+3)=H(k2+3,k1+3)+d2E_dz2*f2_z*f1_z
!!$             H(k2+3,k2+3)=H(k2+3,k2+3)+d2E_dz2*f1_z*f1_z
             H(k1+3,k1+3)=H(k1+3,k1+3)+d2E_dz2
             H(k1+3,k2+3)=H(k1+3,k2+3)-d2E_dz2
             H(k2+3,k1+3)=H(k2+3,k1+3)-d2E_dz2
             H(k2+3,k2+3)=H(k2+3,k2+3)+d2E_dz2

             if(calc_strain) then
                !d2E_dek_del
                H(ns+1,ns+1)=H(ns+1,ns+1)+E_buf_z
                H(ns+1,ns+2)=H(ns+1,ns+2)+E_buf_z; H(ns+2,ns+1)=H(ns+2,ns+1)+E_buf_z
                   
                H(ns+2,ns+2)=H(ns+2,ns+2)+E_buf_z

                !d2E_drk_del
                k1=3*(i-1); k2=3*(jj-1)

!!$                H(k1+1,ns+1)=H(k1+1,ns+1)+dE_dz; H(ns+1,k1+1)=H(ns+1,k1+1)+dE_dz
!!$                H(k2+1,ns+1)=H(k2+1,ns+1)-dE_dz; H(ns+1,k2+1)=H(ns+1,k2+1)-dE_dz
!!$                H(k1+1,ns+2)=H(k1+1,ns+2)+dE_dz; H(ns+2,k1+1)=H(ns+2,k1+1)+dE_dz
!!$                H(k2+1,ns+2)=H(k2+1,ns+2)-dE_dz; H(ns+2,k2+1)=H(ns+2,k2+1)-dE_dz
!!$
!!$                H(k1+2,ns+1)=H(k1+2,ns+1)+dE_dz; H(ns+1,k1+2)=H(ns+1,k1+2)+dE_dz
!!$                H(k2+2,ns+1)=H(k2+2,ns+1)-dE_dz; H(ns+1,k2+2)=H(ns+1,k2+2)-dE_dz
!!$                H(k1+2,ns+2)=H(k1+2,ns+2)+dE_dz; H(ns+2,k1+2)=H(ns+2,k1+2)+dE_dz
!!$                H(k2+2,ns+2)=H(k2+2,ns+2)-dE_dz; H(ns+2,k2+2)=H(ns+2,k2+2)-dE_dz

                H(k1+3,ns+1)=H(k1+3,ns+1)+dE_dz; H(ns+1,k1+3)=H(ns+1,k1+3)+dE_dz
                H(k2+3,ns+1)=H(k2+3,ns+1)-dE_dz; H(ns+1,k2+3)=H(ns+1,k2+3)-dE_dz
                H(k1+3,ns+2)=H(k1+3,ns+2)+dE_dz; H(ns+2,k1+3)=H(ns+2,k1+3)+dE_dz
                H(k2+3,ns+2)=H(k2+3,ns+2)-dE_dz; H(ns+2,k2+3)=H(ns+2,k2+3)-dE_dz
             end if
          end if

          r=rj(1:2)-ri(1:2)

          do ik=1,n_kvec
             dk2=dk(ik)*dk(ik)
             krd=dot_product(kvector(:,ik),r)
             kzd=dk(ik)*z
             expkzd=exp(kzd); exp_kzd=exp(-kzd)
             coskrd=cos(krd)/dk(ik); sinkrd=sin(krd)/dk(ik)

             er1=dk(ik)/(twoep)+ez
             erfunc1=erf_erfc(er1,ierfc)
             er2=dk(ik)/(twoep)-ez
             erfunc2=erf_erfc(er2,ierfc)

             eeee=expkzd*erfunc1+exp_kzd*erfunc2
             cf2=cqq*factor(ik)

             E_buf=cf2*coskrd*eeee
             E_total=E_total+E_buf
             E_coulomb=E_coulomb+E_buf
             E_ew_r=E_ew_r+E_buf

             if(calc_gradients) then
                exp12=exp(-er1*er1); exp22=exp(-er2*er2)
                dE_dr=-cf2*sinkrd*eeee
                gggg=dk(ik)*expkzd*erfunc1- &
                     expkzd*twoep*exp12/sqrt_pi- &
                     dk(ik)*exp_kzd*erfunc2+ &
                     exp_kzd*twoep*exp22/sqrt_pi
                dE_dz=cf2*coskrd*gggg
             end if
             
             if(calc_hessian) then
                hhhh=dk2*expkzd*erfunc1-dk(ik)*expkzd*(twoep/sqrt_pi)*exp12- &
                     dk(ik)*expkzd*(twoep/sqrt_pi)*exp12+ &
                     expkzd*(twoep/sqrt_pi)*two*er1*ewald_param*exp12+ &
                     dk2*exp_kzd*erfunc2-dk(ik)*exp_kzd*(twoep/sqrt_pi)*exp22- &
                     dk(ik)*exp_kzd*(twoep/sqrt_pi)*exp22+ &
                     exp_kzd*(twoep/sqrt_pi)*two*er2*ewald_param*exp22

                d2E_dz2=cf2*coskrd*hhhh
                d2E_drdz=-cf2*sinkrd*gggg
                d2E_dzdr=d2E_drdz
                d2E_dr2=-cf2*coskrd*eeee
             end if

             if(calc_gradients) then
!!$                f1_z=one
!!$                f2_z=-f1_z
                Grad(1:2,i)=Grad(1:2,i)-dE_dr*kvector(:,ik)
                Grad(1:2,jj)=Grad(1:2,jj)+dE_dr*kvector(:,ik)
!!$                Grad(3,i)=Grad(3,i)+dE_dz*f2_z
!!$                Grad(3,jj)=Grad(3,jj)+dE_dz*f1_z
                Grad(3,i)=Grad(3,i)-dE_dz
                Grad(3,jj)=Grad(3,jj)+dE_dz

                kv11=kvector(1,ik)*kvector(1,ik)
                kv12=kvector(1,ik)*kvector(2,ik)
                kv22=kvector(2,ik)*kvector(2,ik)

                if(calc_strain) then
                   ffff=-z*expkzd*erfunc1+expkzd*cf1*exp12+ &
                        z*exp_kzd*erfunc2+exp_kzd*cf1*exp22
                   ffff1=ffff/(eeee*dk(ik))+one/dk2

                   Grad_s(1)=Grad_s(1)-E_buf*(one-kv11*ffff1)
                   Grad_s(2)=Grad_s(2)-E_buf*(one-kv22*ffff1)
                   Grad_s(3)=Grad_s(3)+E_buf*kv12*ffff1
                end if
             end if

             if(calc_hessian) then
                ff(1,1)= kv11
                ff(1,2)= kv12
!!$                ff(1,3)= f2_z*(-kvector(1,ik))
                ff(1,3)= kvector(1,ik)
                ff(1,4)=-kv11
                ff(1,5)=-kv12
!!$                ff(1,6)= f1_z*(-kvector(1,ik))
                ff(1,6)=-kvector(1,ik)

                ff(2,1)= kv12
                ff(2,2)= kv22
!!$                ff(2,3)= f2_z*(-kvector(2,ik))
                ff(2,3)= kvector(2,ik)
                ff(2,4)=-kv12
                ff(2,5)=-kv22
!!$                ff(2,6)= f1_z*(-kvector(2,ik))
                ff(2,6)=-kvector(2,ik)

!!$                ff(3,1)=-kvector(1,ik)*f2_z
!!$                ff(3,2)=-kvector(2,ik)*f2_z
!!$                ff(3,3)= f2_z*f2_z
!!$                ff(3,4)= kvector(1,ik)*f2_z
!!$                ff(3,5)= kvector(2,ik)*f2_z
!!$                ff(3,6)= f1_z*f2_z
                ff(3,1)= kvector(1,ik)
                ff(3,2)= kvector(2,ik)
                ff(3,3)= one
                ff(3,4)=-kvector(1,ik)
                ff(3,5)=-kvector(2,ik)
                ff(3,6)=-one

                ff(4,1)=-kv11
                ff(4,2)=-kv12
!!$                ff(4,3)= f2_z*kvector(1,ik)
                ff(4,3)=-kvector(1,ik)
                ff(4,4)= kv11
                ff(4,5)= kv12
!!$                ff(4,6)= f1_z*kvector(1,ik)
                ff(4,6)= kvector(1,ik)

                ff(5,1)=-kv12
                ff(5,2)=-kv22
!!$                ff(5,3)= f2_z*kvector(2,ik)
                ff(5,3)=-kvector(2,ik)
                ff(5,4)= kv12
                ff(5,5)= kv22
!!$                ff(5,6)= f1_z*kvector(2,ik)
                ff(5,6)= kvector(2,ik)

!!$                ff(6,1)=-kvector(1,ik)*f1_z
!!$                ff(6,2)=-kvector(2,ik)*f1_z
!!$                ff(6,3)= f2_z*f1_z
!!$                ff(6,4)= kvector(1,ik)*f1_z
!!$                ff(6,5)= kvector(2,ik)*f1_z
!!$                ff(6,6)= f1_z*f1_z
                ff(6,1)=-kvector(1,ik)
                ff(6,2)=-kvector(2,ik)
                ff(6,3)=-one
                ff(6,4)= kvector(1,ik)
                ff(6,5)= kvector(2,ik)
                ff(6,6)= one

                k1=3*(i-1)
                k2=3*(jj-1)
                do l=1,2
                   do m=1,2
                      H(k1+l,k1+m)=H(k1+l,k1+m)+d2E_dr2*ff(l,m)
                      H(k1+l,k2+m)=H(k1+l,k2+m)+d2E_dr2*ff(l,m+3)
                      H(k2+l,k1+m)=H(k2+l,k1+m)+d2E_dr2*ff(l+3,m)
                      H(k2+l,k2+m)=H(k2+l,k2+m)+d2E_dr2*ff(l+3,m+3)
                   end do
                   H(k1+l,k1+3)=H(k1+l,k1+3)+d2E_drdz*ff(l,3)
                   H(k1+l,k2+3)=H(k1+l,k2+3)+d2E_drdz*ff(l,6)
                   H(k2+l,k1+3)=H(k2+l,k1+3)+d2E_drdz*ff(l+3,3)
                   H(k2+l,k2+3)=H(k2+l,k2+3)+d2E_drdz*ff(l+3,6)
                end do
                do m=1,2
                   H(k1+3,k1+m)=H(k1+3,k1+m)+d2E_dzdr*ff(3,m)
                   H(k1+3,k2+m)=H(k1+3,k2+m)+d2E_dzdr*ff(3,m+3)
                   H(k2+3,k1+m)=H(k2+3,k1+m)+d2E_dzdr*ff(6,m)
                   H(k2+3,k2+m)=H(k2+3,k2+m)+d2E_dzdr*ff(6,m+3)
                end do
                H(k1+3,k1+3)=H(k1+3,k1+3)+d2E_dz2*ff(3,3)
                H(k1+3,k2+3)=H(k1+3,k2+3)+d2E_dz2*ff(3,6)
                H(k2+3,k1+3)=H(k2+3,k1+3)+d2E_dz2*ff(6,3)
                H(k2+3,k2+3)=H(k2+3,k2+3)+d2E_dz2*ff(6,6)

                if(calc_strain) then
                   !d2E_dek_del
                   iiii=(dz2*expkzd*erfunc1-two*cf1*expkzd*exp12+ &
                        two*cf1*expkzd*exp12*er1+ &
                        dz2*exp_kzd*erfunc2+two*cf1*exp_kzd*exp22+ &
                        two*cf1*exp_kzd*exp22*er2)/(eeee*dk2)- &
                        ffff*ffff/(eeee*eeee*dk2)+ffff/(eeee*dk(ik)*dk2)+ &
                        two/(dk2*dk2)

                   H(ns+1,ns+1)=H(ns+1,ns+1)+E_buf*(one-kv11*ffff1)*(one-kv11*ffff1)- &
                        E_buf*ffff1*two*kv11+E_buf*kv11*kv11*iiii
                   HH=E_buf*(one-kv11*ffff1)*(one-kv22*ffff1)+E_buf*kv11*kv22*iiii
                   H(ns+1,ns+2)=H(ns+1,ns+2)+HH; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
                   HH=-E_buf*(one-kv11*ffff1)*kv12*ffff-E_buf*ffff1*kv12+ &
                        E_buf*kv11*kv12*iiii
                   H(ns+1,ns+3)=H(ns+1,ns+3)+HH; H(ns+3,ns+1)=H(ns+3,ns+1)+HH

                   H(ns+2,ns+2)=H(ns+2,ns+2)+E_buf*(one-kv22*ffff1)*(one-kv22*ffff1)- &
                        E_buf*ffff1*two*kv22+E_buf*kv22*kv22*iiii
                   HH=-E_buf*(one-kv22*ffff1)*kv12*ffff1-E_buf*ffff1*kv12+ &
                        E_buf*kv22*kv12*iiii
                   H(ns+2,ns+3)=H(ns+2,ns+3)+HH; H(ns+3,ns+2)=H(ns+3,ns+2)+HH

                   H(ns+3,ns+3)=H(ns+3,ns+3)+E_buf*kv12*ffff1*kv12*ffff1- &
                        E_buf*ffff1*half*(kv11+kv22)+E_buf*kv12*kv12*iiii

                   !d2E_drk_del
                   HH=dE_dr*(one-kv11/dk2+kv11*ffff/(eeee*dk(ik)))*kvector(1,ik)+dE_dr*kvector(1,ik)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   H(k2+1,ns+1)=H(k2+1,ns+1)-HH; H(ns+1,k2+1)=H(ns+1,k2+1)-HH
                   HH=dE_dr*(one-kv22/dk2+kv22*ffff/(eeee*dk(ik)))*kvector(1,ik)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   H(k2+1,ns+2)=H(k2+1,ns+2)-HH; H(ns+2,k2+1)=H(ns+2,k2+1)-HH
                   HH=dE_dr*(-kv12/dk2+kv12*ffff/(eeee*dk(ik)))*kvector(1,ik)+dE_dr*kvector(2,ik)*half
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   H(k2+1,ns+3)=H(k2+1,ns+3)-HH; H(ns+3,k2+1)=H(ns+3,k2+1)-HH

                   HH=dE_dr*(one-kv11/dk2+kv11*ffff/(eeee*dk(ik)))*kvector(2,ik)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   H(k2+2,ns+1)=H(k2+2,ns+1)-HH; H(ns+1,k2+2)=H(ns+1,k2+2)-HH
                   HH=dE_dr*(one-kv22/dk2+kv22*ffff/(eeee*dk(ik)))*kvector(2,ik)+dE_dr*kvector(2,ik)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   H(k2+2,ns+2)=H(k2+2,ns+2)-HH; H(ns+2,k2+2)=H(ns+2,k2+2)-HH
                   HH=dE_dr*(-kv12/dk2+kv12*ffff/(eeee*dk(ik)))*kvector(2,ik)+dE_dr*kvector(1,ik)*half
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   H(k2+2,ns+3)=H(k2+2,ns+3)-HH; H(ns+3,k2+2)=H(ns+3,k2+2)-HH

                   if(gggg /= zero) then
                      jjjj=kzd*expkzd*erfunc1-cf1*z*expkzd*exp12- &
                           (twoep/sqrt_pi)*z*expkzd*exp12+(twoep/sqrt_pi)*two*er1*expkzd*exp12+ &
                           kzd*exp_kzd*erfunc2+cf1*z*exp_kzd*exp22- &
                           (twoep/sqrt_pi)*z*exp_kzd*exp22-(twoep/sqrt_pi)*two*er2*exp_kzd*exp22

                      HH=-dE_dz*(-one+kv11/dk2-kv11*jjjj/(gggg*dk(ik)))
                      H(k1+3,ns+1)=H(k1+3,ns+1)+HH; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                      H(k2+3,ns+1)=H(k2+3,ns+1)-HH; H(ns+1,k2+3)=H(ns+1,k2+3)-HH
                      HH=-dE_dz*(-one+kv22/dk2-kv22*jjjj/(gggg*dk(ik)))
                      H(k1+3,ns+2)=H(k1+3,ns+2)+HH; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                      H(k2+3,ns+2)=H(k2+3,ns+2)-HH; H(ns+2,k2+3)=H(ns+2,k2+3)-HH
                      HH=-dE_dz*(kv12/dk2-kv12*jjjj/(gggg*dk(ik)))
                      H(k1+3,ns+3)=H(k1+3,ns+3)+HH; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                      H(k2+3,ns+3)=H(k2+3,ns+3)-HH; H(ns+3,k2+3)=H(ns+3,k2+3)-HH
                   end if
                end if
             end if
          end do
       end do
    end do

  end subroutine recipr_ew2d_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine calc_2d_poten(n_points,rj,Vj,dVj_dri,do_grads)

    integer(i4_kind), intent(in) :: n_points       !number of points
    real(r8_kind), intent(in)    :: rj(:,:)        !coordinates of points potential to be calclulated
    real(r8_kind), intent(inout) :: Vj(:)          !calculated potential
    real(r8_kind), intent(inout) :: dVj_dri(:,:,:) !calculated potential first derivatives with respect
                                                   !to atomic positions
    logical, intent(in)          :: do_grads

    integer(i4_kind) :: j,i,im,a_type_i,ik
    real(r8_kind) :: qi,cqi,rji(3),drji,drji2,er,er2,ri(3)
    real(r8_kind) :: dV_drji,fi(3),erfunc,dk2,zji,ez
    real(r8_kind) :: ez2,erfunc_z,expez2,dV_dzji,xyji(2)
    real(r8_kind) :: coskrd,sinkrd,dV_dxyji,eeee,ep2,er1
    real(r8_kind) :: erfunc1,erfunc2,exp12,exp22,expkzd,exp_kzd
    real(r8_kind) :: gggg,krd,kzd

    !direct space
    if(.not.do_grads) then
       Vj=zero
    else
       dVj_dri=zero
    end if

    do j=1,n_points
       lab_i:do i=1,n_species
          a_type_i=atoms_cart(i)%type
          qi=atoms(a_type_i)%charge
          cqi=qi

          !calculating potential from atoms located in the real unit cell
          rji=rj(j,:)-atoms_cart(i)%r
          if(minimal_image) call image_slab(rji)
          drji=sqrt(dot_product(rji,rji))
          if(drji > rcut_ew) cycle lab_i
          drji2=drji*drji
          er=ewald_param*drji; er2=er*er
          erfunc=erf_erfc(er,ierfc)

          if(.not.do_grads) Vj(j)=Vj(j)+cqi*erfunc/drji
          
          if(do_grads) dV_drji=cqi*(-erfunc/drji2-two*ewald_param*exp(-er2)/(drji*sqrt_pi))
          if(do_grads) then
             fi=-rji/drji
             dVj_dri(j,i,:)=dVj_dri(j,i,:)+dV_drji*fi
          end if

          if(.not.minimal_image) then
             !calculating potential from atoms located in unit cell images
             lab_im:do im=1,n_images
                rji=rj(j,:)-im_coor(i)%r(:,im)
                drji=sqrt(dot_product(rji,rji))
                if(drji > rcut_ew) cycle lab_im
                drji2=drji*drji
                er=ewald_param*drji; er2=er*er
                erfunc=erf_erfc(er,ierfc)

                if(.not.do_grads) Vj(j)=Vj(j)+cqi*erfunc/drji
                
                if(do_grads) dV_drji=cqi*(-erfunc/drji2-two*ewald_param*exp(-er2)/(drji*sqrt_pi))
                if(do_grads) then
                   fi=-rji/drji
                   dVj_dri(j,i,:)=dVj_dri(j,i,:)+dV_drji*fi
                end if
             end do lab_im
          end if
       end do lab_i
    end do

    !reciprocal space
    ep2=ewald_param*ewald_param
    do j=1,n_points
       lab1_i:do i=1,n_species
          a_type_i=atoms_cart(i)%type
          qi=atoms(a_type_i)%charge
          cqi=(pi/area)*qi

          ri=atoms_cart(i)%r

          zji=rj(j,3)-ri(3)
          ez=ewald_param*zji;ez2=ez*ez
          erfunc_z=erf_erfc(ez,ierf)
          expez2=exp(-ez2)

          if(.not.do_grads) Vj(j)=Vj(j)-two*cqi*(zji*erfunc_z+(one/(ewald_param*sqrt_pi))*expez2)

          if(do_grads) dV_dzji=-two*cqi*(erfunc_z+ &
               zji*(two*ewald_param)*expez2/sqrt_pi- &
               two*(one/(ewald_param*sqrt_pi))*ep2*zji*expez2)

          if(do_grads) then
             dVj_dri(j,i,3)=dVj_dri(j,i,3)-dV_dzji
          end if

          xyji=rj(j,1:2)-ri(1:2)

          do ik=1,n_kvec
             dk2=dk(ik)*dk(ik)
             krd=dot_product(kvector(:,ik),xyji)
             kzd=dk(ik)*zji
             expkzd=exp(kzd); exp_kzd=exp(-kzd)
             coskrd=cos(krd)/dk(ik); sinkrd=sin(krd)/dk(ik)

             er1=dk(ik)/(two*ewald_param)+ez
             erfunc1=erf_erfc(er1,ierfc)
             er2=dk(ik)/(two*ewald_param)-ez
             erfunc2=erf_erfc(er2,ierfc)

             eeee=expkzd*erfunc1+exp_kzd*erfunc2

             if(.not.do_grads) Vj(j)=Vj(j)+cqi*factor(ik)*coskrd*eeee

             if(do_grads) then
                exp12=exp(-er1*er1); exp22=exp(-er2*er2)
                dV_dxyji=-cqi*factor(ik)*sinkrd*eeee
                gggg=dk(ik)*expkzd*erfunc1- &
                     expkzd*(two*ewald_param)*exp12/sqrt_pi- &
                     dk(ik)*exp_kzd*erfunc2+ &
                     exp_kzd*(two*ewald_param)*exp22/sqrt_pi
                dV_dzji=cqi*factor(ik)*coskrd*gggg
             end if

             if(do_grads) then
                dVj_dri(j,i,1:2)=dVj_dri(j,i,1:2)-dV_dxyji*kvector(:,ik)
                dVj_dri(j,i,3)=dVj_dri(j,i,3)-dV_dzji
             end if
          end do
       end do lab1_i
    end do
    
  end subroutine calc_2d_poten
  !****************************************************************

  !****************************************************************
  subroutine shutdown_ewald2d()

    integer(i4_kind) :: status

    if(allocated(k_indexes)) then 
       deallocate(k_indexes,stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: Failed deallocation of K_INDEXES (2D)")
    end if

    if(allocated(kvector)) then
       deallocate(kvector,dk,factor,stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: failed deallocation of KVECTOR (2D)")
    end if

  end subroutine shutdown_ewald2d
  !****************************************************************

  !****************************************************************
end module ewald2d_module

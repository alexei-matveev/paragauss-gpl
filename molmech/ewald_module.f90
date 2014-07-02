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
module ewald_module
  !------------ Modules used -----------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module
  use species_module
  use energy_and_forces_module
  use n_body_lists_module
  use potentials_module, only : n_cs
  use coulomb_module
  use molmech_msgtag_module
  use comm_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ---------------------
  public init_ewald, self_ewald_energy,direct_ew_E_and_F,recipr_ew_E_and_F, &
       send_receive_ewald, calc_gauss, shutdown_ewald,recipr_ew_E_and_F_t, &
       erf_erfc,calc_3d_poten
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  real(r8_kind) :: w_factor=one
  real(r8_kind) :: accuracy=small
  real(r8_kind) :: ewald_param ! beta
  real(r8_kind) :: ewp3
  real(r8_kind) :: kcut_ew
  real(r8_kind) :: sqrt_pi
  integer(i4_kind), parameter :: max_kv = 200 

  integer(i4_kind), allocatable :: k_indexes(:)
  integer(i4_kind) :: n_kvec
  integer(i4_kind) :: n_kvec_on_node
  real(r8_kind), allocatable :: gauss_multiplier(:)
  real(r8_kind), allocatable :: kvector(:,:)

  integer(i4_kind) :: n_proc

  real(r8_kind) :: E_self,E_self0
  !------------ Subroutines ------------------------------------------
contains
  !****************************************************************
  subroutine init_ewald()

    real(r8_kind) :: f

    ewald_param=sqrt(((n_species*w_factor*pi**3)/volume**2)**(one/three))
    ewp3=ewald_param*ewald_param*ewald_param
    sqrt_pi=sqrt(pi)
    f=sqrt(-log(accuracy))
    rcut_ew=f/ewald_param
    kcut_ew=two*f*ewald_param

    call calc_kindexes()

    n_proc=comm_get_n_processors()

  end subroutine init_ewald
  !****************************************************************

  !****************************************************************
  subroutine calc_kindexes()

    integer(i4_kind) :: max1,max2,max3
    real(r8_kind) :: drk,rk1(3),rk2(3),rk3(3),rk(3),cut2
    integer(i4_kind) :: i,j,k,status

    type store
       integer(i4_kind) :: index_store
       type(store), pointer :: next_data
    end type store
    type(store), target :: first_data
    type(store), pointer ::   current_data, tmp_data

    rk1=kvect%v1
    drk=sqrt(dot_product(rk1,rk1))
    max1=int(kcut_ew/drk)+1

    rk2=kvect%v2
    drk=sqrt(dot_product(rk2,rk2))
    max2=int(kcut_ew/drk)+1

    rk3=kvect%v3
    drk=sqrt(dot_product(rk3,rk3))
    max3=int(kcut_ew/drk)+1

    if(max1 > max_kv .or. max2 > max_kv .or. max3 > max_kv) &
         call error_handler("MolMech: Number of K vectors is too large")

    current_data=>first_data
    nullify(current_data%next_data)

    cut2=kcut_ew*kcut_ew

    n_kvec=0
    do i=0,max1
       rk1=i*kvect%v1
       do j=-max2,max2
          rk2=j*kvect%v2
          do k=-max3,max3
             if(i==0.and.j==0.and.k==0) cycle
             rk3=k*kvect%v3
             rk=rk1+rk2+rk3
             drk=dot_product(rk,rk)
             if(drk > cut2) cycle
             n_kvec=n_kvec+1
             allocate(tmp_data,stat=status)
             if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(ew3)")
             tmp_data%index_store=1000000*(i+max_kv)+1000*(j+max_kv)+(k+max_kv)
             nullify(tmp_data%next_data)
             current_data%next_data =>tmp_data
             current_data => tmp_data
          end do
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
  subroutine send_receive_ewald()

    integer(i4_kind) :: info,status,i,nkv,ilow,iup

    if (comm_i_am_master()) then
       n_kvec_on_node = nint (real (n_kvec) / real (n_proc))
       ilow=1
       do i=2,n_proc
          call comm_init_send(i,msgtag_mm_send_ewald)

          call commpack(ewald_param,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: ewald_param pack failed")
          call commpack(kvect%v1(1),3,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v1 pack failed")
          call commpack(kvect%v2(1),3,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v2 pack failed")
          call commpack(kvect%v3(1),3,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v3 pack failed")
          call commpack(volume,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: volume pack failed")

          ilow=ilow+n_kvec_on_node
          iup=ilow+n_kvec_on_node-1
          if(i==n_proc) iup=n_kvec
          nkv=iup-ilow+1
          call commpack(nkv,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: nkv pack failed")
          call commpack(k_indexes(ilow),nkv,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_ewald: k_indexes pack failed")

          call comm_send()
       end do
    else
       call communpack(ewald_param,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: ewald_param unpack failed")
       call communpack(kvect%v1(1),3,1,info)
       if( info /= 0) call error_handler &
               ("send_receive_ewald: kvect%v1 unpack failed")
       call communpack(kvect%v2(1),3,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: kvect%v2 unpack failed")
       call communpack(kvect%v3(1),3,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: kvect%v3 unpack failed")
       call communpack(volume,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: volume unpack failed")

       call communpack(n_kvec_on_node,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: nkv unpack failed")
       allocate(k_indexes(n_kvec_on_node),stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: Failed allocation of K_INDEXES on slave")
       
       call communpack(k_indexes(1),n_kvec_on_node,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_ewald: k_indexes unpack failed")
    end if

  end subroutine send_receive_ewald
  !****************************************************************

  !****************************************************************
  subroutine self_ewald_energy()

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

  end subroutine  self_ewald_energy
  !****************************************************************

  !****************************************************************
  subroutine direct_ew_E_and_F()

    integer(i4_kind) :: i,j,jj,k,kk,l,m,k1,k2,m1,length,length1,count
    integer(i4_kind) :: img,i1,j1
    integer(i4_kind) :: a_type1,a_type2,ns
    real(r8_kind) :: q1,q2,r(3),dr,sfac,erfun,erfunc,er
    real(r8_kind) :: dr2,dr3,cqq,er2
    real(r8_kind) :: f1(3),f2(3),ff(6,6),fstore
    real(r8_kind) :: E_buf,dE_dr,d2E_dr2,V,VV,HH
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
             if(lattice_calc .and. minimal_image) call image(r)
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
                   Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                   Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                   Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                   Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
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
!!$                      H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!!$                      H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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
                      Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                      Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                      Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                      Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
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
!!$                         H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!!$                         H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
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

                      !d2E_dri_dek
                      k1=3*(i-1); k2=3*(k-1)

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
             if(lattice_calc) call image(r)
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
                   Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                   Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                   Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                   Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
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
!!$                      H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)    
!!$                      H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)      
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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
                   Grad_s(3)=Grad_s(3)+V*r(3)*r(3)
                   Grad_s(4)=Grad_s(4)+V*r(2)*r(3)
                   Grad_s(5)=Grad_s(5)+V*r(1)*r(3)
                   Grad_s(6)=Grad_s(6)+V*r(1)*r(2)
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
!!$                   H(k1+l,m1)=H(k1+l,m1)+d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
!!$                   H(k2+l,m1)=H(k2+l,m1)+d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
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

                   !d2E_dri_dek
                   k1=3*(i-1); k2=3*(k-1)

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
                end if
             end if
          end do
       end do
    end do

  end subroutine direct_ew_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine calc_gauss()

    integer(i4_kind) :: ii,i,j,k,k_ind,status
    real(r8_kind) :: kv(3),dkv,ep2,factor

    allocate(gauss_multiplier(n_kvec_on_node),kvector(3,n_kvec_on_node),stat=status)
    if(status /= 0) call error_handler("MolMech: failed allocation of GAUSS_MULTIPLIER")

    ep2=ewald_param*ewald_param

    do ii=1,n_kvec_on_node
       factor=two
       k_ind=k_indexes(ii)
       i=k_ind/1000000-max_kv
       j=(k_ind-(i+max_kv)*1000000)/1000-max_kv
       k=k_ind-(i+max_kv)*1000000-(j+max_kv)*1000-max_kv
       if(i==0) factor=one

       kv=i*kvect%v1+j*kvect%v2+k*kvect%v3
       kvector(:,ii)=kv
       dkv=dot_product(kv,kv)

       gauss_multiplier(ii)=factor*coulomb_factor*four*pi*exp(-dkv/(four*ep2))/(dkv*volume)
    end do

  end subroutine calc_gauss
  !****************************************************************

  !****************************************************************
  subroutine recipr_ew_E_and_F()

    integer(i4_kind) :: ii,jj,kk,i,a_type,k1
    integer(i4_kind) :: i1,i2,i3,j1,j2,j3,ns,status
    real(r8_kind) :: qcos,qsin,qcos2,qsin2,q,q1
    real(r8_kind) :: r(3),kr,dkv,kv11,kv12,kv13,kv22,kv23,kv33
    real(r8_kind) :: E_buf,HH,aaa,aaa1,bbb,g_buf,qcs,cs
    real(r8_kind) :: qgm,sin2,cos2,ep2,coss,sinn,coss1,sinn1
    real(r8_kind),allocatable :: krd(:),cosd(:),sind(:)

    ns=3*n_species

    ep2=ewald_param*ewald_param

    allocate(krd(N_species),cosd(N_species),sind(N_species),stat=status)
    if(status /= 0) call error_handler("MolMech: failed KRD allocation(ew3)")

    do ii=1,n_kvec_on_node
       qcos=zero; qsin=zero
       do jj=1,N_species
          a_type=atoms_cart(jj)%type
          q=atoms(a_type)%charge
          r=atoms_cart(jj)%r
          krd(jj)=dot_product(kvector(:,ii),r)
          kr=krd(jj) !dot_product(kvector(:,ii),r)
          cosd(jj)=cos(kr); sind(jj)=sin(kr)
          qcos=qcos+q*cosd(jj) !cos(kr)
          qsin=qsin+q*sind(jj) !sin(kr)
       end do
       qcos2=qcos*qcos
       qsin2=qsin*qsin

       E_buf=gauss_multiplier(ii)*(qcos2+qsin2)*half
       E_ew_r=E_ew_r+E_buf
       E_coulomb=E_coulomb+E_buf
       E_total=E_total+E_buf

       if(calc_gradients) then
          kv11=kvector(1,ii)*kvector(1,ii)
          kv12=kvector(1,ii)*kvector(2,ii)
          kv13=kvector(1,ii)*kvector(3,ii)
          kv22=kvector(2,ii)*kvector(2,ii)
          kv23=kvector(2,ii)*kvector(3,ii)
          kv33=kvector(3,ii)*kvector(3,ii)
       end if

       if(calc_gradients .and. calc_strain) then
          dkv=dot_product(kvector(:,ii),kvector(:,ii))
          aaa=two*(one/dkv+one/(four*ep2)); aaa1=aaa

          Grad_s(1)=Grad_s(1)-E_buf*(one-kv11*aaa)
          Grad_s(2)=Grad_s(2)-E_buf*(one-kv22*aaa)
          Grad_s(3)=Grad_s(3)-E_buf*(one-kv33*aaa)
          Grad_s(4)=Grad_s(4)+E_buf*kv23*aaa
          Grad_s(5)=Grad_s(5)+E_buf*kv13*aaa
          Grad_s(6)=Grad_s(6)+E_buf*kv12*aaa
          
          if(calc_hessian) then
             !d2E_dek_del
             bbb=aaa*aaa+four/(dkv*dkv)
             aaa=E_buf*aaa
             bbb=E_buf*bbb

             H(ns+1,ns+1)=H(ns+1,ns+1)+E_buf*two-five*kv11*aaa+kv11*kv11*bbb- &
                  (E_buf-kv11*aaa)
             HH=E_buf-(kv11+kv22)*aaa+kv11*kv22*bbb
             H(ns+1,ns+2)=H(ns+1,ns+2)+HH!; H(ns+2,ns+1)=H(ns+2,ns+1)+HH
             HH=E_buf-(kv11+kv33)*aaa+kv11*kv33*bbb
             H(ns+1,ns+3)=H(ns+1,ns+3)+HH!; H(ns+3,ns+1)=H(ns+3,ns+1)+HH
             HH=-kv23*aaa+kv11*kv23*bbb
             H(ns+1,ns+4)=H(ns+1,ns+4)+HH!; H(ns+4,ns+1)=H(ns+4,ns+1)+HH
             HH=-two*kv13*aaa+kv11*kv13*bbb
             H(ns+1,ns+5)=H(ns+1,ns+5)+HH!; H(ns+5,ns+1)=H(ns+5,ns+1)+HH
             HH=-two*kv12*aaa+kv11*kv12*bbb
             H(ns+1,ns+6)=H(ns+1,ns+6)+HH!; H(ns+6,ns+1)=H(ns+6,ns+1)+HH

             H(ns+2,ns+2)=H(ns+2,ns+2)+E_buf*two-five*kv22*aaa+kv22*kv22*bbb- &
                  (E_buf-kv22*aaa)
             HH=E_buf-(kv22+kv33)*aaa+kv22*kv33*bbb
             H(ns+2,ns+3)=H(ns+2,ns+3)+HH!; H(ns+3,ns+2)=H(ns+3,ns+2)+HH
             HH=-two*kv23*aaa+kv22*kv23*bbb
             H(ns+2,ns+4)=H(ns+2,ns+4)+HH!; H(ns+4,ns+2)=H(ns+4,ns+2)+HH
             HH=-kv13*aaa+kv22*kv13*bbb     
             H(ns+2,ns+5)=H(ns+2,ns+5)+HH!; H(ns+5,ns+2)=H(ns+5,ns+2)+HH
             HH=-two*kv12*aaa+kv22*kv12*bbb
             H(ns+2,ns+6)=H(ns+2,ns+6)+HH!; H(ns+6,ns+2)=H(ns+6,ns+2)+HH

             H(ns+3,ns+3)=H(ns+3,ns+3)+E_buf*two-five*kv33*aaa+kv33*kv33*bbb- &
                  (E_buf-kv33*aaa)
             HH=-two*kv23*aaa+kv33*kv23*bbb
             H(ns+3,ns+4)=H(ns+3,ns+4)+HH!; H(ns+4,ns+3)=H(ns+4,ns+3)+HH
             HH=-two*kv13*aaa+kv33*kv13*bbb
             H(ns+3,ns+5)=H(ns+3,ns+5)+HH!; H(ns+5,ns+3)=H(ns+5,ns+3)+HH
             HH=-kv12*aaa+kv33*kv12*bbb
             H(ns+3,ns+6)=H(ns+3,ns+6)+HH!; H(ns+6,ns+3)=H(ns+6,ns+3)+HH

!!$             H(ns+4,ns+4)=H(ns+4,ns+4)+E_buf*half-quarter3*(kv33+kv22)*aaa+kv23*kv23*bbb- &
!!$                  quarter*(E_buf-kv33*aaa+E_buf-kv22*aaa)
             H(ns+4,ns+4)=H(ns+4,ns+4)-half*(kv33+kv22)*aaa+kv23*kv23*bbb
             HH=-half*kv12*aaa+kv23*kv13*bbb
             H(ns+4,ns+5)=H(ns+4,ns+5)+HH!; H(ns+5,ns+4)=H(ns+5,ns+4)+HH
             HH=-half*kv13*aaa+kv23*kv12*bbb
             H(ns+4,ns+6)=H(ns+4,ns+6)+HH!; H(ns+6,ns+4)=H(ns+6,ns+4)+HH
             
!!$             H(ns+5,ns+5)=H(ns+5,ns+5)+E_buf*half-quarter3*(kv33+kv11)*aaa+kv13*kv13*bbb- &
!!$                  quarter*(E_buf-kv33*aaa+E_buf-kv11*aaa)
             H(ns+5,ns+5)=H(ns+5,ns+5)-half*(kv33+kv11)*aaa+kv13*kv13*bbb
             HH=-half*kv23*aaa+kv13*kv12*bbb
             H(ns+5,ns+6)=H(ns+5,ns+6)+HH!; H(ns+6,ns+5)=H(ns+6,ns+5)+HH
             
!!$             H(ns+6,ns+6)=H(ns+6,ns+6)+E_buf*half-quarter3*(kv22+kv11)*aaa+kv12*kv12*bbb- &
!!$                  quarter*(E_buf-kv22*aaa+E_buf-kv11*aaa)
             H(ns+6,ns+6)=H(ns+6,ns+6)-half*(kv22+kv11)*aaa+kv12*kv12*bbb
          end if
       end if

       if(calc_gradients .or. calc_hessian) then
          do jj=1,N_species
             a_type=atoms_cart(jj)%type
             q=atoms(a_type)%charge
!!$             r=atoms_cart(jj)%r
!!$             kr=krd(jj) !dot_product(kvector(:,ii),r)

             qgm=q*gauss_multiplier(ii)
             coss=cosd(jj) !cos(kr)
             sinn=sind(jj) !sin(kr)
             g_buf=qgm*(sinn*qcos-coss*qsin)

             Grad(:,jj)=Grad(:,jj)-g_buf*kvector(:,ii)

             if(calc_hessian) then
                if(calc_strain) then
                   !d2E_dri_dek
                   k1=3*(jj-1)

                   HH=g_buf*(kvector(1,ii)-kv11*kvector(1,ii)*aaa1)
                   H(k1+1,ns+1)=H(k1+1,ns+1)+HH!; H(ns+1,k1+1)=H(ns+1,k1+1)+HH
                   HH=g_buf*(kvector(1,ii)-kv22*kvector(1,ii)*aaa1)
                   H(k1+1,ns+2)=H(k1+1,ns+2)+HH!; H(ns+2,k1+1)=H(ns+2,k1+1)+HH
                   HH=g_buf*(kvector(1,ii)-kv33*kvector(1,ii)*aaa1)
                   H(k1+1,ns+3)=H(k1+1,ns+3)+HH!; H(ns+3,k1+1)=H(ns+3,k1+1)+HH
                   HH=g_buf*(-kv23*kvector(1,ii)*aaa1)
                   H(k1+1,ns+4)=H(k1+1,ns+4)+HH!; H(ns+4,k1+1)=H(ns+4,k1+1)+HH
                   HH=g_buf*(-kv13*kvector(1,ii)*aaa1)
                   H(k1+1,ns+5)=H(k1+1,ns+5)+HH!; H(ns+5,k1+1)=H(ns+5,k1+1)+HH
                   HH=g_buf*(-kv12*kvector(1,ii)*aaa1)
                   H(k1+1,ns+6)=H(k1+1,ns+6)+HH!; H(ns+6,k1+1)=H(ns+6,k1+1)+HH
                
                   HH=g_buf*(kvector(2,ii)-kv11*kvector(2,ii)*aaa1)
                   H(k1+2,ns+1)=H(k1+2,ns+1)+HH!; H(ns+1,k1+2)=H(ns+1,k1+2)+HH
                   HH=g_buf*(kvector(2,ii)-kv22*kvector(2,ii)*aaa1)
                   H(k1+2,ns+2)=H(k1+2,ns+2)+HH!; H(ns+2,k1+2)=H(ns+2,k1+2)+HH
                   HH=g_buf*(kvector(2,ii)-kv33*kvector(2,ii)*aaa1)
                   H(k1+2,ns+3)=H(k1+2,ns+3)+HH!; H(ns+3,k1+2)=H(ns+3,k1+2)+HH
                   HH=g_buf*(-kv23*kvector(2,ii)*aaa1)
                   H(k1+2,ns+4)=H(k1+2,ns+4)+HH!; H(ns+4,k1+2)=H(ns+4,k1+2)+HH
                   HH=g_buf*(-kv13*kvector(2,ii)*aaa1)
                   H(k1+2,ns+5)=H(k1+2,ns+5)+HH!; H(ns+5,k1+2)=H(ns+5,k1+2)+HH
                   HH=g_buf*(-kv12*kvector(2,ii)*aaa1)
                   H(k1+2,ns+6)=H(k1+2,ns+6)+HH!; H(ns+6,k1+2)=H(ns+6,k1+2)+HH
                
                   HH=g_buf*(kvector(3,ii)-kv11*kvector(3,ii)*aaa1)
                   H(k1+3,ns+1)=H(k1+3,ns+1)+HH!; H(ns+1,k1+3)=H(ns+1,k1+3)+HH
                   HH=g_buf*(kvector(3,ii)-kv22*kvector(3,ii)*aaa1)
                   H(k1+3,ns+2)=H(k1+3,ns+2)+HH!; H(ns+2,k1+3)=H(ns+2,k1+3)+HH
                   HH=g_buf*(kvector(3,ii)-kv33*kvector(3,ii)*aaa1)
                   H(k1+3,ns+3)=H(k1+3,ns+3)+HH!; H(ns+3,k1+3)=H(ns+3,k1+3)+HH
                   HH=g_buf*(-kv23*kvector(3,ii)*aaa1)
                   H(k1+3,ns+4)=H(k1+3,ns+4)+HH!; H(ns+4,k1+3)=H(ns+4,k1+3)+HH
                   HH=g_buf*(-kv13*kvector(3,ii)*aaa1)
                   H(k1+3,ns+5)=H(k1+3,ns+5)+HH!; H(ns+5,k1+3)=H(ns+5,k1+3)+HH
                   HH=g_buf*(-kv12*kvector(3,ii)*aaa1)
                   H(k1+3,ns+6)=H(k1+3,ns+6)+HH!; H(ns+6,k1+3)=H(ns+6,k1+3)+HH
                end if

                sin2=sinn*sinn
                cos2=coss*coss
                qcs=qgm*(coss*qcos-q*sin2+sinn*qsin-q*cos2)

                i1=3*(jj-1)+1; i2=i1+1; i3=i1+2

                H(i1,i1)=H(i1,i1)-kv11*qcs
                HH=-kv12*qcs
                H(i1,i2)=H(i1,i2)+HH
                !H(i2,i1)=H(i2,i1)+HH
                HH=-kv13*qcs
                H(i1,i3)=H(i1,i3)+HH
                !H(i3,i1)=H(i3,i1)+HH
                
                H(i2,i2)=H(i2,i2)-kv22*qcs
                HH=-kv23*qcs
                H(i2,i3)=H(i2,i3)+HH
                !H(i3,i2)=H(i3,i2)+HH

                H(i3,i3)=H(i3,i3)-kv33*qcs

                do kk=jj+1,N_species
                   j1=3*(kk-1)+1; j2=j1+1; j3=j1+2
                   a_type=atoms_cart(kk)%type
                   q1=atoms(a_type)%charge
!!$                   r=atoms_cart(kk)%r
!!$                   kr1=krd(kk) !dot_product(kvector(:,ii),r)
                   coss1=cosd(kk) !cos(kr1)
                   sinn1=sind(kk) !sin(kr1)
                   cs=q1*qgm*(-sinn*sinn1-coss*coss1)
            
                   HH=-kv11*cs
                   H(i1,j1)=H(i1,j1)+HH
                   !H(j1,i1)=H(j1,i1)+HH
                   HH=-kv12*cs
                   H(i1,j2)=H(i1,j2)+HH
                   !H(j2,i1)=H(j2,i1)+HH
                   H(i2,j1)=H(i2,j1)+HH
                   !H(j1,i2)=H(j1,i2)+HH
                   HH=-kv13*cs
                   H(i1,j3)=H(i1,j3)+HH
                   !H(j3,i1)=H(j3,i1)+HH
                   H(i3,j1)=H(i3,j1)+HH
                   !H(j1,i3)=H(j1,i3)+HH
            
                   HH=-kv22*cs
                   H(i2,j2)=H(i2,j2)+HH
                   !H(j2,i2)=H(j2,i2)+HH
                   HH=-kv23*cs
                   H(i2,j3)=H(i2,j3)+HH
                   !H(j3,i2)=H(j3,i2)+HH
                   H(i3,j2)=H(i3,j2)+HH
                   !H(j2,i3)=H(j2,i3)+HH
            
                   HH=-kv33*cs
                   H(i3,j3)=H(i3,j3)+HH
                   !H(j3,i3)=H(j3,i3)+HH
                end do
                
!!$                do l=1,3
!!$                   do n=l,3
!!$                      HH=-kvector(l,ii)*kvector(n,ii)*qcs
!!$                      H(i1+l,i1+n)=H(i1+l,i1+n)+HH
!!$                      if(l /= n) H(i1+n,i1+l)=H(i1+n,i1+l)+HH
!!$                   end do
!!$                   do kk=jj+1,N_species
!!$                      j1=3*(kk-1)
!!$                      a_type=atoms_cart(kk)%type
!!$                      q1=atoms(a_type)%charge
!!$                      r=atoms_cart(kk)%r
!!$                      kr1=krd(kk) !dot_product(kvector(:,ii),r)
!!$                      coss1=cosd(kk) !cos(kr1)
!!$                      sinn1=sind(kk) !sin(kr1)
!!$                      cs=qgm*(-sinn*sinn1-coss*coss1)
!!$                      do m=1,3
!!$                         HH=-q1*kvector(l,ii)*kvector(m,ii)*cs
!!$                         H(i1+l,j1+m)=H(i1+l,j1+m)+HH
!!$                         H(j1+m,i1+l)=H(j1+m,i1+l)+HH
!!$                      end do
!!$                   end do
!!$                end do
             end if
          end do
       end if
    end do

    if(calc_hessian) then
       if(calc_strain) ns=ns+6
       do i=1,ns-1
          H(i+1:ns,i)=H(i,i+1:ns)
       end do
    end if

    deallocate(krd,cosd,sind,stat=status)
    if(status /= 0) call error_handler("MolMech: failed KRD deallocation(ew3)")

  end subroutine recipr_ew_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine recipr_ew_E_and_F_t
    
    integer(i4_kind) :: i,j,jj,k,k1,k2,l,m,m1,a_type
    real(r8_kind) :: qi,qj,ri(3),rj(3),E_buf,kri,krj
    real(r8_kind) :: ff(6,6)
    real(r8_kind) :: cossin,qqgm,qqgmc

    E_buf=zero
    do i=1,n_species
       a_type=atoms_cart(i)%type
       qi=atoms(a_type)%charge
       do k=1,n_kvec
          E_buf=E_buf+qi*qi*gauss_multiplier(k)*half
       end do
    end do
    E_total=E_total+E_buf
    E_coulomb=E_coulomb+E_buf
    E_ew_r=E_ew_r+E_buf
      
    do i=1,n_species
       a_type=atoms_cart(i)%type
       qi=atoms(a_type)%charge
       ri=atoms_cart(i)%r
       do j=1,N_total(i)
          jj=i+j
          if(jj > n_species) jj=jj-n_species
          a_type=atoms_cart(jj)%type
          qj=atoms(a_type)%charge
          rj=atoms_cart(jj)%r

          E_buf=zero
          do k=1,n_kvec
             kri=dot_product(kvector(:,k),ri)
             krj=dot_product(kvector(:,k),rj)
             cossin=cos(kri)*cos(krj)+sin(kri)*sin(krj)
             qqgm=qi*qj*gauss_multiplier(k)
             qqgmc=qqgm*cossin

             E_buf=E_buf+qi*qj*gauss_multiplier(k)*cossin

             if(calc_gradients) then
                Grad(:,i) =Grad(:,i) + &
                     qqgm*kvector(:,k)*(-sin(kri)*cos(krj)+cos(kri)*sin(krj))
                Grad(:,jj)=Grad(:,jj)+ &
                     qqgm*kvector(:,k)*(-cos(kri)*sin(krj)+sin(kri)*cos(krj))
             end if

             if(calc_hessian) then
                ff(1,1:3)=-qqgmc*kvector(1,k)*kvector(1:3,k)
                ff(1,4:6)= qqgmc*kvector(1,k)*kvector(1:3,k)
                
                ff(2,1)  = ff(1,2)
                ff(2,2:3)=-qqgmc*kvector(2,k)*kvector(2:3,k)
                ff(2,4:6)= qqgmc*kvector(2,k)*kvector(1:3,k)

                ff(3,1:2)= ff(1:2,3)
                ff(3,3)  =-qqgmc*kvector(3,k)*kvector(3,k)
                ff(3,4:6)= qqgmc*kvector(3,k)*kvector(1:3,k)

                ff(4,1:3)= ff(1:3,4)
                ff(4,4:6)= -qqgmc*kvector(1,k)*kvector(1:3,k)

                ff(5,1:4)= ff(1:4,5)
                ff(5,5:6)= -qqgmc*kvector(2,k)*kvector(2:3,k)

                ff(6,1:5)= ff(1:5,6)
                ff(6,6)  = -qqgmc*kvector(3,k)*kvector(3,k)

                k1=3*(i-1)
                k2=3*(jj-1)
                do l=1,3
                   do m=1,6
                      m1=k1+m
                      if(m > 3) then
                         m1=k2+(m-3)
                      end if
                      H(k1+l,m1)=H(k1+l,m1)+ff(l,m)
                      H(k2+l,m1)=H(k2+l,m1)+ff(l+3,m)
                   end do
                end do
             end if
          end do
          E_total=E_total+E_buf
          E_coulomb=E_coulomb+E_buf
          E_ew_r=E_ew_r+E_buf

       end do
    end do

  end subroutine recipr_ew_E_and_F_t
  !****************************************************************

  !****************************************************************
  subroutine calc_3d_poten(n_points,rj,Vj,dVj_dri,do_grads)

    integer(i4_kind), intent(in) :: n_points       !number of points
    real(r8_kind), intent(in)    :: rj(:,:)        !coordinates of points potential to be calclulated
    real(r8_kind), intent(inout) :: Vj(:)          !calculated potential
    real(r8_kind), intent(inout) :: dVj_dri(:,:,:) !calculated potential first derivatives with respect
                                                   !to atomic positions
    logical, intent(in)          :: do_grads

    integer :: i,j,a_type_i,im,k
    real(r8_kind) :: qi,cqi,rji(3),drji,drji2,er,erfunc,dV_drji
    real(r8_kind) :: fi(3),er2,krji

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

          !Direct space
          !calculating potential from atoms located in the real unit cell
          rji=rj(j,:)-atoms_cart(i)%r
          if(minimal_image) call image(rji)
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

          !Reciprocal space
          lab_K: do k=1,n_kvec
             krji=dot_product(kvector(:,k),rji)

             if(.not.do_grads) Vj(j)=Vj(j)+(qi/coulomb_factor)*cos(krji)*gauss_multiplier(k)

             if(do_grads) then
                dV_drji=(qi/coulomb_factor)*sin(krji)*gauss_multiplier(k)
                dVj_dri(j,i,:)=dVj_dri(j,i,:)+dV_drji*kvector(:,k)
             end if
          end do lab_K
       end do lab_i
    end do

  end subroutine calc_3d_poten
  !****************************************************************

  !****************************************************************
  function erf_erfc(arg_inp,func)

    real(r8_kind) :: erf_erfc

    real(r8_kind) :: arg_inp
    integer(i4_kind) :: func

    integer(i4_kind), parameter :: eerf=0
    integer(i4_kind), parameter :: eerfc=1
    !Mathematical constants
    real(r8_kind), parameter :: ezero=0.0_r8_kind
    real(r8_kind), parameter :: ehalf=0.5_r8_kind
    real(r8_kind), parameter :: eone =1.0_r8_kind
    real(r8_kind), parameter :: etwo =2.0_r8_kind
    real(r8_kind), parameter :: efour=4.0_r8_kind
    real(r8_kind), parameter :: esixteen=16.0_r8_kind
    real(r8_kind), parameter :: thresh=0.46875_r8_kind
    real(r8_kind), parameter :: esqrt_pi=5.6418958354775628695e-1_r8_kind
    !Machine-dependent constants
    real(r8_kind), parameter :: xinf=1.79e308_r8_kind
    real(r8_kind), parameter :: xneg=-26.628_r8_kind
    real(r8_kind), parameter :: xsmall=1.11e-16_r8_kind
    real(r8_kind), parameter :: xbig=26.543_r8_kind
    real(r8_kind), parameter :: xhuge=6.71e7_r8_kind
    real(r8_kind), parameter :: xmax=2.53e307_r8_kind
    !Coefficients for approximation to  erf  in first interval
    real(r8_kind) :: a(5)=(/3.16112374387056560_r8_kind, &
                            1.13864154151050156e2_r8_kind, &
                            3.77485237685302021e2_r8_kind, &
                            3.20937758913846947e3_r8_kind, &
                            1.85777706184603153e-1_r8_kind/)
    real(r8_kind) :: b(4)=(/2.36012909523441209e1_r8_kind, &
                            2.44024637934444173e2_r8_kind, &
                            1.28261652607737228e3_r8_kind, &
                            2.84423683343917062e3_r8_kind/)
    !Coefficients for approximation to  erfc  in second interval
    real(r8_kind) :: c(9)=(/5.64188496988670089e-1_r8_kind, &
                            8.88314979438837594_r8_kind, &
                            6.61191906371416295e1_r8_kind, &
                            2.98635138197400131e2_r8_kind, &
                            8.81952221241769090e2_r8_kind, &
                            1.71204761263407058e3_r8_kind, &
                            2.05107837782607147e3_r8_kind, &
                            1.23033935479799725e3_r8_kind, &
                            2.15311535474403846e-8_r8_kind/)
    real(r8_kind) :: d(8)=(/1.57449261107098347e1_r8_kind, &
                            1.17693950891312499e2_r8_kind, &
                            5.37181101862009858e2_r8_kind, &
                            1.62138957456669019e3_r8_kind, &
                            3.29079923573345963e3_r8_kind, &
                            4.36261909014324716e3_r8_kind, &
                            3.43936767414372164e3_r8_kind, &
                            1.23033935480374942e3_r8_kind/)
    !Coefficients for approximation to  erfc  in third interval
    real(r8_kind) :: p(6)=(/3.05326634961232344e-1_r8_kind, &
                            3.60344899949804439e-1_r8_kind, &
                            1.25781726111229246e-1_r8_kind, &
                            1.60837851487422766e-2_r8_kind, &
                            6.58749161529837803e-4_r8_kind, &
                            1.63153871373020978e-2_r8_kind/)
    real(r8_kind) :: q(5)=(/2.56852019228982242_r8_kind, &
                            1.87295284992346047_r8_kind, &
                            5.27905102951428412e-1_r8_kind, &
                            6.05183413124413191e-2_r8_kind, &
                            2.33520497626869185e-3_r8_kind/)

    integer(i4_kind) :: i
    real(r8_kind) :: del,x,y,ysq,xden,xnum

    x=arg_inp
    y=abs(x)
    if(y <= thresh) then
       !Evaluate  erf  for  |X| <= 0.46875
       ysq=ezero
       if(y > xsmall) ysq=y*y
       xnum=a(5)*ysq
       xden=ysq
       do i=1,3
          xnum=(xnum+a(i))*ysq
          xden=(xden+b(i))*ysq
       end do
       erf_erfc=x*(xnum+a(4))/(xden+b(4))
       if(func /= eerf) erf_erfc=eone-erf_erfc
       return
    elseif(y > thresh .and. y <= efour) then
       !Evaluate  erfc  for 0.46875 <= |X| <= 4.0
       xnum=c(9)*y
       xden=y
       do i=1,7
          xnum=(xnum+c(i))*y
          xden=(xden+d(i))*y
       end do
       erf_erfc=(xnum+c(8))/(xden+d(8))
       ysq=aint(y*esixteen)/esixteen
       del=(y-ysq)*(y+ysq)
       erf_erfc=exp(-ysq*ysq)*exp(-del)*erf_erfc
    else
       !Evaluate  erfc  for |X| > 4.0
       erf_erfc=ezero
       if(y >= xbig) then
          if(y >= xmax) goto 1
          if(y >= xhuge) then 
             erf_erfc=esqrt_pi/y
             goto 1
          end if
       end if
       ysq=eone/(y*y)
       xnum=p(6)*ysq
       xden=ysq
       do i=1,4
          xnum=(xnum+p(i))*y
          xden=(xden+q(i))*y
       end do
       erf_erfc=ysq*(xnum+p(5))/(xden+q(5))
       erf_erfc=(esqrt_pi-erf_erfc)/y
       ysq=aint(y*esixteen)/esixteen
       del=(y-ysq)*(y+ysq)
       erf_erfc=exp(-ysq*ysq)*exp(-del)*erf_erfc
    end if

1   if(func == eerf) then
       erf_erfc=(ehalf-erf_erfc)+ehalf
       if(x < ezero) erf_erfc=-erf_erfc
    elseif(func == eerfc) then
       if(x < ezero) erf_erfc=etwo-erf_erfc
    end if

  end function erf_erfc
  !****************************************************************

  !****************************************************************
  subroutine shutdown_ewald()

    integer(i4_kind) :: status

    if(allocated(k_indexes)) then 
       deallocate(k_indexes,stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: Failed deallocation of K_INDEXES")
    end if

    if(allocated(gauss_multiplier)) then
       deallocate(gauss_multiplier,kvector,stat=status)
       if(status /= 0) call error_handler &
            ("MolMech: failed deallocation of GAUSS_MULTIPLIER")
    end if

  end subroutine shutdown_ewald
  !****************************************************************

  !****************************************************************
end module ewald_module

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
!**************************************************************
subroutine energy_and_grad_of_cavity()
  !
  ! Calculates the cavitation energy or its gradients (scaled particle
  ! theory) formulas from
  !
# include "def.h"
  use type_module, only: i4_kind, r8_kind ! type specification parameters
  use unique_atom_module, only: N_moving_unique_atoms, unique_atoms, &
       moving_unique_atom_index, unique_atom_grad_info
  use solv_cavity_module, only: energy_cav, tessarea, grad_solv_totsym, &
       skip_short, abs_temperature, cavitation_all, do_gradients, fixed_pc, &
       n_size, n_spheres, solvent_radius, solvent_volume, with_pc, iuniq, &
       r_sphere, i_symm_sort, cagr, dealloc_geom_deriv_part2, dealloc_cavity
  use pointcharge_module, only: pointcharge_N, pointcharge_array
  use elec_static_field_module, only: surf_points_grad_info, surf_points_grad_index, &
       totalsym_field
  use integralpar_module, only: integralpar_2dervs
  use gradient_data_module, only: grad_index=>gradient_index
  use gradient_data_module, only: dervs_totsym=>dervs_totalsym
  implicit none
  !** End of interface *****************************************

  real(kind=r8_kind) :: Kb, P, Vsol, num_dens, rsolv
  real(kind=r8_kind) :: K0, K1, K2, K3, dd, Y, YY, Ecav, area_total, area_s
  integer(kind=i4_kind) :: i,i1,j,k,l,m,ma,ma1,na,ea,grad_dim,index
  integer(i4_kind) :: mb,nb,eb,n,n1,l1,grad_dim1,index1
  real(kind=r8_kind),pointer  :: rotmat(:,:),rotmat1(:,:)
  real (kind=r8_kind) :: gradient(3),hessian(3,3),buffer(3,3)
  integer :: N_pc,ism

  real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind
  real(kind=r8_kind) , parameter :: ang_au = 0.529177249_r8_kind
  real(kind=r8_kind) , parameter :: au_to_erg = 4.35926565522e-11_r8_kind
  real(kind=r8_kind) , parameter :: au_to_dyn = 8.237e-3_r8_kind
  real(kind=r8_kind) , parameter :: avogadro = 6.0221367e23_r8_kind  !mol^(-1)
  real(kind=r8_kind) , parameter :: pressure = 1.01325e6_r8_kind      !dyn/cm^2
  real(kind=r8_kind) , parameter :: boltzmann = 1.380658e-16_r8_kind   !erg/K

  DPRINT 'energy_and_grad_of_cavity: entered'

  !convert some parameters to atomic units
  Kb=boltzmann/au_to_erg
  P=(pressure*ang_au*ang_au/au_to_dyn)*1.0e-16_r8_kind
  Vsol=solvent_volume*1.0e24_r8_kind/(ang_au*ang_au*ang_au)
  rsolv=solvent_radius/ang_au

  num_dens=avogadro/Vsol

  if(.not.do_gradients) then
     energy_cav=0.0_r8_kind
     
     Y=4.0_r8_kind*pi*num_dens*rsolv**3/3.0_r8_kind
     YY=Y/(1.0_r8_kind-Y)
     K0= -log(1.0_r8_kind-Y)
     K1= 3.0_r8_kind*YY
     K2= 3.0_r8_kind*YY+4.5_r8_kind*YY**2
!!$          K3= Y*P/(num_dens*Kb*abs_temperature)
     K3=0.0_r8_kind

     do_nsp_i: do i=1,N_spheres
        if(.not. cavitation_all) then
           !use skip_short only when allocated!
           if(skip_short(iuniq(i))) cycle do_nsp_i
        endif
        dd=r_sphere(i)/rsolv
        
        Ecav=(K0+K1*dd+K2*dd**2+K3*dd**3)*Kb*abs_temperature
        
        area_total=4.0_r8_kind*pi*r_sphere(i)**2
        
        area_s=0.0_r8_kind
        do j=1,n_size
           do k=1,tessarea(j)%n_equal
              if(tessarea(j)%sphere(k) == i) area_s=area_s+tessarea(j)%area
           enddo
        enddo
        energy_cav=energy_cav+Ecav*area_s/area_total
     enddo do_nsp_i
  else
     N_pc=0
     if(with_pc .and. .not.fixed_pc)  N_pc=pointcharge_N
     
     Y=4.0_r8_kind*pi*num_dens*rsolv**3/3.0_r8_kind
     YY=Y/(1.0_r8_kind-Y)
     K0= -log(1.0_r8_kind-Y)
     K1= 3.0_r8_kind*YY
     K2= 3.0_r8_kind*YY+4.5_r8_kind*YY**2
!!$             K3= Y*P/(num_dens*Kb*abs_temperature)
     K3=0.0_r8_kind             
     
     unique1: do ma=1,N_moving_unique_atoms+N_pc
        if(ma <= N_moving_unique_atoms) then
           na=moving_unique_atom_index(ma)
           ea=unique_atoms(na)%n_equal_atoms
        else
           na=ma-N_moving_unique_atoms
           ea=pointcharge_array(na)%N_equal_charges
        end if
        
        unique2: do m=1,ea
           gradient=0.0_r8_kind
           
           i_N: do i=1,N_spheres
              if(.not. cavitation_all) then
                 !use skip_short only when allocated!
                 if(skip_short(iuniq(i))) cycle i_N
              endif
              dd=r_sphere(i)/rsolv
              
              Ecav=(K0+K1*dd+K2*dd**2+K3*dd**3)*Kb*abs_temperature
              
              area_total=4.0_r8_kind*pi*r_sphere(i)**2
              
              DPRINT 'i_symm_sort used here, shape=',shape(i_symm_sort)
              do j=1,n_size
                 do k=1,tessarea(j)%n_equal
                    ism=i_symm_sort(j,k)
                    if(tessarea(j)%sphere(k) == i) then
                       do l=1,3
                          gradient(l)=gradient(l)+&
                               Ecav*cagr%darea(l,ma,m)%m(ism)/area_total
!!$print*,l,ma,m,j,k,cagr%darea(l,ma,m)%m(i_symm_sort(j,k)),ism
                       enddo
                    endif
                 enddo
              enddo
           enddo i_N

           if(ma <= N_moving_unique_atoms) then
              grad_dim=grad_index(ma+1)-grad_index(ma)
              rotmat=>unique_atom_grad_info(ma)%m(:,:,m)
              index=grad_index(ma)
              do i1=1,grad_dim
!!! ??? Sign ???
!               grad_solv_totsym(index) = grad_solv_totsym(index) - &
!                    sum( rotmat(i1,:) * gradient(:) )
                 grad_solv_totsym(index) = grad_solv_totsym(index) + &
                      sum( rotmat(i1,:) * gradient(:) )
                 index=index+1
              end do
           else
              ma1=ma-N_moving_unique_atoms
              grad_dim=surf_points_grad_index(ma1+1)-surf_points_grad_index(ma1)
              rotmat=>surf_points_grad_info(ma1)%m(:,:,m)
              index=surf_points_grad_index(ma1)
              do i1=1,grad_dim
                 totalsym_field(index) = totalsym_field(index) - &
                      sum( rotmat(i1,:) * gradient(:) )
                 index=index+1
              end do
           end if
        enddo unique2
     enddo unique1
     
     if(integralpar_2dervs) then
        unique_a: do ma=1,N_moving_unique_atoms
           na=moving_unique_atom_index(ma)
           ea=unique_atoms(na)%n_equal_atoms
           grad_dim=grad_index(ma+1)-grad_index(ma)
           
           unique_b: do mb=1,N_moving_unique_atoms
              nb=moving_unique_atom_index(mb)
              eb=unique_atoms(nb)%n_equal_atoms
              grad_dim1=grad_index(mb+1)-grad_index(mb)
              
              do n=1,ea
                 rotmat=>unique_atom_grad_info(ma)%m(:,:,n)
                 do n1=1,eb
                    rotmat1=>unique_atom_grad_info(mb)%m(:,:,1)
                    
                    hessian=0.0_r8_kind
                    
                    i_Ns:do i=1,N_spheres
                       if(.not. cavitation_all) then
                          !use skip_short only when allocated!
                          if(skip_short(iuniq(i))) cycle i_Ns
                       endif
                       dd=r_sphere(i)/rsolv
                       
                       Ecav=(K0+K1*dd+K2*dd**2+K3*dd**3)*Kb*abs_temperature
                       
                       area_total=4.0_r8_kind*pi*r_sphere(i)**2
                       
                       do j=1,n_size
                          do k=1,tessarea(j)%n_equal
                             ism=i_symm_sort(j,k)
                             if(tessarea(j)%sphere(k) == i) then
                                do l=1,3
                                   do l1=1,3
                                      hessian(l,l1)=hessian(l,l1)+&
                                           Ecav*cagr%d2area(l,ma,n,l1,mb,n1)%m(ism)/area_total
                                   end do
                                enddo
                             endif
                          enddo
                       enddo
                    end do i_Ns

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
                       index1=grad_index(mb)
                       do k=1,grad_dim1
                          do l=1,3
                             dervs_totsym(index,index1)=dervs_totsym(index,index1)+ &
                                  rotmat1(k,l)*buffer(j,l)
                          end do
                          index1=index1+1
                       end do
                       index=index+1
                    end do
                 end do
              end do
           end do unique_b
        end do unique_a
     end if
  endif

  if(do_gradients) then
     call dealloc_geom_deriv_part2
  endif

  call dealloc_cavity()

  if(allocated(skip_short)) then
     !maybe not allocated if dispersion repulsion not calculated
     deallocate(skip_short,stat=i)
     if(i/=0) call error_handler("energy_and_grad_of_cavity: dealloc of skip short failed")
  endif

end subroutine energy_and_grad_of_cavity
!**************************************************************

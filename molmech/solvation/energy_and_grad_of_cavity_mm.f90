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
subroutine energy_and_grad_of_cavity_mm()
  !calculates the cavitation energy or its gradients 
  !(scaled particle theory) formulas from
  !
  !** End of interface *****************************************
# include <def.h>
  use type_module ! type specification parameters
  use datatype
  use cavity_module
  use energy_and_forces_module
  use common_data_module
  use tasks_main_options_module
  use species_module

  implicit none

  real(kind=r8_kind) :: Kb, Vsol, num_dens, rsolv
  real(kind=r8_kind) :: K0, K1, K2, K3, dd, Y, YY, Ecav, area_total, area_s
  integer(kind=i4_kind) :: i,j,k,l,ma
  real (kind=r8_kind) :: gradient(3)

  Kb=(boltzmann/h2erg)*h2kJm ! (in kJm/K)
  Vsol=solvent_volume*1.0e24_r8_kind
  rsolv=solvent_radius

  num_dens=avogadro/Vsol

  E_solv_cav=zero
  do_nsp_i: do i=1,N_spheres
     Y=four*pi*num_dens*rsolv**3/three
     YY=Y/(one-Y)
     
     dd=r_sphere(i)/rsolv

     K0= -log(one-Y)
     K1= three*YY
     K2= three*YY+4.5_r8_kind*YY**2
     K3= zero

     Ecav=(K0+K1*dd+K2*dd**2+K3*dd**3)*Kb*abs_temperature
       
     area_total=four*pi*r_sphere(i)**2

     area_s=zero
     do j=1,n_size
        do k=1,tessarea(j)%n_equal
           if(tessarea(j)%sphere(k) == i) area_s=area_s+tessarea(j)%area
        enddo
     enddo
     E_solv_cav=E_solv_cav+Ecav*area_s/area_total
  enddo do_nsp_i
  E_solv_tot=E_solv_tot+E_solv_cav
  E_total=E_total+E_solv_cav
  if(calc_gradients) then
     unique1: do ma=1,n_species       
        gradient=zero
        
        i_N: do i=1,N_spheres
           Y=four*pi*num_dens*rsolv**3/three
           YY=Y/(one-Y)
           
           dd=r_sphere(i)/rsolv
           
           K0= -log(one-Y)
           K1= three*YY
           K2= three*YY+4.5_r8_kind*YY**2
           K3= zero
           
           Ecav=(K0+K1*dd+K2*dd**2+K3*dd**3)*Kb*abs_temperature

           area_total=four*pi*r_sphere(i)**2

           do j=1,n_size
              do k=1,tessarea(j)%n_equal
                 if(tessarea(j)%sphere(k) == i) then
                    do l=1,3
                       gradient(l)=gradient(l)+&
                            Ecav*cagr%darea(l,ma,1)%m(j)/area_total
                    enddo
                 endif
              enddo
           enddo
        enddo i_N
!!$print*,ma,gradient/(h2kJm*a2b)
        Grad(:,ma)=Grad(:,ma)+gradient
     enddo unique1
  endif
  if(calc_gradients) then
     call dealloc_geom_deriv_part2
  endif

  call dealloc_cavity_mm()

end subroutine energy_and_grad_of_cavity_mm
!**************************************************************

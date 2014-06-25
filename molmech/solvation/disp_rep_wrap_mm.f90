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
# include <def.h>
  subroutine disp_rep_wrap_mm()
   !calls the routines for calculating the dispersion/repulsion
   !energy or its gradients
   !** End of interface *****************************************
    
  use type_module ! type specification parameters
  use datatype
  use common_data_module
  use tasks_main_options_module
  use species_module
  use energy_and_forces_module
  use cavity_module
  use vdwcm_module
  use inp_out_module, only: output_device

  implicit none

  real(kind=r8_kind) :: Vsol
  real(kind=r8_kind) :: gradient(3) !,gr(9,3)
  integer(kind=i4_kind) :: i,status,ma

  call disp_rep_read() !!!!!!!!!!!!1

  Vsol=solvent_volume*1.0e24_r8_kind
  num_dens_dr=avogadro/Vsol

  allocate(nm_solute_at(n_species),stat=status)
  if ( status /= 0) call error_handler( &
       "MolMech:disp_rep_wrap: allocation of nm_solute_at is failed")

  if(.not. allocated(skip_short)) then
     allocate(skip_short(n_species),stat=status)
     if ( status /= 0) call error_handler( &
          "MolMech:disp_rep_wrap: allocation of skip_short is failed")
  endif

  call initialize_data(skip_short)

  E_disp=zero
  E_rep=zero

!!$gr=zero

  do i=1,ndsa
     if(output_level >=1) write (output_device,*) 'THE CAVITY -',i
     R_access=r_sphere_at(i)
     call points_on_cavity_surface()
     call solv_disprep_transer_data_mm()
     call disp_rep_energies_and_grad(i,do_energy=.true.,do_grad=.false.)
     if(calc_gradients) then
        unique2: do ma=1,n_species
           gradient=0.0_r8_kind
           call disp_rep_energies_and_grad(i,do_energy=.false.,do_grad=.true.,gradient=gradient,ma=ma)
           Grad(:,ma)=Grad(:,ma)+gradient
!!$             gr(ma,:)=gr(ma,:)+gradient
        enddo unique2
     end if
     call deallocate_D_R_points(calc_gradients)
  enddo
  E_solv_dr=E_disp+E_rep
  E_solv_tot=E_solv_tot+E_solv_dr
  E_total=E_total+E_solv_dr

!!$do ma=1,n_species
!!$print*,ma,gr(ma,:)/(h2kJm*a2b)
!!$end do

  call deallocate_parameters(calc_gradients)
  deallocate(nm_solute_at,stat=status)
  if ( status /= 0) call error_handler( &
       "disp_rep_wrap: deallocation of nm_solute_at is failed")

contains

  subroutine solv_disprep_transer_data_mm()
    !stores data needed for calculating the dispersion/repulsion
    !energy or its gradients in variables
    !D_R_points
    !D_R_grad
    !** End of interface *****************************************

    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: j1,k,l,m,n
    integer(kind=i4_kind) :: i_sym

    N_points_dr=n_size

    allocate(D_R_points(N_points_dr),stat=status)
    if (calc_gradients) then
       allocate(D_R_grad(N_points_dr),stat=status)
    endif

    do l=1,N_points_dr
       if(calc_gradients) then
          allocate(D_R_grad(l)&
               & %position  (3,n_species,3),&
               D_R_grad(l)&
               & %out_normal(3,n_species,3), &
               D_R_grad(l)&
               & %area      (3,n_species), &
               D_R_grad(l)&
               & %distance  (3,n_species,n_species,3), &
               stat=status)
       endif

       D_R_points(l)%position(:)=tessarea(l)%xyz(1,:)
       k=tessarea(l)%sphere(1)
       D_R_points(l)%out_normal(:)=(tessarea(l)%xyz(1,:)-xyz_sphere(k,:))/ &
            sqrt(dot_product(tessarea(l)%xyz(1,:)-xyz_sphere(k,:), &
            tessarea(l)%xyz(1,:)-xyz_sphere(k,:)))

       if(calc_gradients) then
          i_sym = l
          do j1=1,n_species
             do n=1,3
                D_R_grad(l)%position(n,j1,:)=&
                     cagr%dcenter(n,j1,1)%m(:,i_sym)
                D_R_grad(l)%area(n,j1)=&
                     cagr%darea(n,j1,1)%m(i_sym)
                    
                do m=1,3 ! coordinate(x,y,z) of point j
                   D_R_grad(l)%out_normal(n,j1,m)=&
                        (cagr%dcenter(n,j1,1)%m(m,i_sym)-&
                        cagr%dc(m,k)%xyz_grad(n,j1,1))/r_sphere(k)
                       
                   D_R_grad(l)%distance(n,j1,:,m)=&
                        D_R_grad(l)%position(n,j1,m)
                enddo
                       
                D_R_grad(l)%distance(n,j1,j1,n)=&
                     D_R_grad(l)%distance(n,j1,j1,n)-1.0_r8_kind
             enddo
          enddo
       endif
       D_R_points(l)%area=tessarea(l)%area
    enddo

    call dealloc_cavity_mm()
    if(calc_gradients) then
       call dealloc_geom_deriv_part1(cagr)
       call dealloc_geom_deriv_part2()
    endif

  end subroutine solv_disprep_transer_data_mm

end subroutine disp_rep_wrap_mm

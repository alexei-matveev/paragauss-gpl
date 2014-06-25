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
!********************************************************************
subroutine calc_solv_eff()

# include <def.h>
  use type_module ! type specification parameters
  use datatype
  use cavity_module
  use solv_elec_stat_module
  use tasks_main_options_module
  use species_module, only: lattice_calc
  use slab_module, only: slab_calc
  use ewald_solv_module

  implicit none

  save            ! save all variables defined in this module

  if(lattice_calc .or. slab_calc) then
     disp_rep_energy=.false.
     cavitation_energy=.false.
  end if

  if(cavitation_energy) then
     do_cavitation=.true.
     do_disp_rep=.false.
     do_electr=.false.
!!$print*,'cav cav'
     call points_on_cavity_surface()
!!$print*,'cav energy'
     call energy_and_grad_of_cavity_mm()
  endif
  if(disp_rep_energy) then
     do_cavitation=.false.
     do_disp_rep=.true.
     do_electr=.false.
!!$print*,'disp-rep cav and energy'
     call disp_rep_wrap_mm()
  endif
  if(electr_energy) then
     do_cavitation=.false.
     do_disp_rep=.false.
     do_electr=.true.
!!$print*,'cav'
     call points_on_cavity_surface()
     if(lattice_calc .and. use_ewald) then
        call init_ewald_3D_solv()
        call calc_gauss_solv()
        call calc_number_of_3D_images()
     end if
     call matrix_generation()
!!$print*,'energy'
     call el_solv_energy()
     if(calc_gradients) then
!!$print*,'grad mat'
        call matrix_grad()
!!$print*,'grad elec'
        call solv_el_grad()
     end if
     call dealloc_cavity_mm()
     call shutdown_ewald_3D_solv()
     if(calc_gradients) then
        call dealloc_geom_deriv_part1(cagr)
        call dealloc_geom_deriv_part2()
     endif
  end if

end subroutine calc_solv_eff
!********************************************************************

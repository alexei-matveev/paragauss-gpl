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
!===============================================================
! Public interface of module
!===============================================================
subroutine lin_search_epe
  !----------------------------------------------------------------
  !
  !  Purpose: Perform primitive line search during epe relaxation
  !           around QM cluster 
  !
  !  Subroutine called by: main_epe_block
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AS
  !  Date: 23.11.09
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
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
#include <def.h>
  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use epecom_module, only: n_ls, output_epe, reg_I_n_ions
  use epecom_module, only: n_vacancies,regI_previous,relax_shells_only
  use epecom_module, only: r_sh_ion,r_nuc_ion,pk,epe,point_0_core_grad,reset

  implicit none

  !== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ------------------
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of subroutines ------------------------
  external error_handler

  !------------ Declaration of local constants --------------------

  !------------ Declaration of local variables --------------------
  real(r8_kind), allocatable, dimension(:) :: rs_dim,rs0_dim,rn_dim,rn0_dim,dsn,step
  integer(i4_kind) :: status,i,j,k

  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------

  n_ls=n_ls+1
  print*,'Line search',n_ls
  write (output_epe,*) 'Line search',n_ls

  allocate(rs_dim(3*(reg_I_n_ions-n_vacancies)), & 
       rs0_dim(3*(reg_I_n_ions-n_vacancies)), &
       rn_dim(3*(reg_I_n_ions-n_vacancies)), &
       rn0_dim(3*(reg_I_n_ions-n_vacancies)), &
       step(3*(reg_I_n_ions-n_vacancies)), &
       dsn(3*(reg_I_n_ions-n_vacancies)),stat=status)
  ASSERT(status==0) 

  ! swap data to one dimentional store
  if(relax_shells_only) then
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           rs_dim(k)=r_sh_ion(i,j)
           rs0_dim(k)=regI_previous(i)%r_shell(j)
           rn_dim(k)=r_nuc_ion(i,j)
           rn0_dim(k)=regI_previous(i)%r_core(j)
        enddo
     enddo
  else
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           rs_dim(k)=r_sh_ion(i,j)
           rs0_dim(k)=regI_previous(i)%r_shell(j)
           rn_dim(k)=r_nuc_ion(i,j)
           rn0_dim(k)=regI_previous(i)%r_core(j)
        enddo
     enddo
  endif

  step=(rs_dim-rs0_dim)*0.25_r8_kind
  dsn=rn0_dim-rs0_dim

  if(relax_shells_only) then
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           rs_dim(k)=rs0_dim(k)+step(k)
        enddo
     enddo
  else
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           rs_dim(k)=rs0_dim(k)+step(k)+ &
                0.5_r8_kind*point_0_core_grad(i,j)/pk(epe(i)%k)
           rn_dim(k)=rs0_dim(k)+step(k)+dsn(k)- &
                0.5_r8_kind*point_0_core_grad(i,j)/pk(epe(i)%k)
        enddo
     enddo
  endif

  if(relax_shells_only) then
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           r_sh_ion(i,j)=rs_dim(k)
        enddo
     enddo
  else
     do i=n_vacancies+1,reg_I_n_ions
        do j=1,3
           k=3*(i-n_vacancies-1)+j
           r_sh_ion(i,j)=rs_dim(k)
           r_nuc_ion(i,j)=rn_dim(k)
        enddo
     enddo
  endif

  deallocate(rs_dim,rs0_dim,rn_dim,rn0_dim,dsn,stat=status)
  ASSERT(status==0) 

  if(n_ls == 2) then
     reset=.true.
  end if

end subroutine lin_search_epe

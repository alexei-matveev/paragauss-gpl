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
!===============================================================
! Public interface of module
!===============================================================
module vdwcm_module

# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use iounitadmin_module, only: output_unit,write_to_output_units
  use common_data_module
  use inp_out_module, only: output_device
  use element_data_module
  use species_module
  use energy_and_forces_module

  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public disp_rep_read, initialize_data, disp_rep_energies_and_grad, deallocate_D_R_points, &
       deallocate_parameters


  ! data type to collect unique surface points, normals and areas
  ! of a cavity tessellation.
  ! idices: (coordinates,number_of_equal_partner)
  type, public :: disp_rep_points
     ! center of tessera
     real(kind=r8_kind)               :: position(3)
     ! outward normal in tessera center
     real(kind=r8_kind)               :: out_normal(3)
     ! area of tessera
     real(kind=r8_kind)               :: area
  end type disp_rep_points

  ! type to collect the gradients of the surface point positions,
  ! the tessera normals, the areas and the distance to all sphere
  ! centers.
  type, public :: disp_rep_points_grad
     real(kind=r8_kind), pointer                :: position(:,:,:)
     real(kind=r8_kind), pointer                :: out_normal(:,:,:)
     real(kind=r8_kind), pointer                :: distance(:,:,:,:)
     real(kind=r8_kind), pointer                :: area(:,:)
  end type disp_rep_points_grad

  ! energy contributions form dispersion and repulsion interaction
  real(kind=r8_kind), public :: E_disp, E_rep
  ! number of unique surface points of the actual cavity
  integer(kind=i4_kind), public :: N_points_dr
  ! position/normals/areas of surface tesserae of the actual cavity
  ! for dispersion/repulsion contribution
  type(disp_rep_points), allocatable, target, public :: D_R_points(:)
  ! gradients of position/normals/areas of surface tesserae of the 
  ! actual cavity
  type(disp_rep_points_grad), allocatable, target, public :: D_R_grad(:)

  ! Names of the solute atom (atoms of the molecule)
  character*2, public, allocatable :: nm_solute_at(:)
  real(kind=r8_kind), public, allocatable :: r_sphere_at(:)
  ! number density of special solvent atom
  real(kind=r8_kind), public :: num_dens_dr
  ! number of different elements in a solvent molecule
  integer(kind=i4_kind), public :: ndsa
  ! if I shall use data from Rappe
  logical, public :: use_rappe_data
!===============================================================
! End of public interface of module
!===============================================================

  real(kind=r8_kind), allocatable :: K_solv(:),R_solv(:),K_solu(:),R_solu(:)

!!$  integer(kind=i4_kind) :: n_differ_solv_atoms
!!$  integer(kind=i4_kind) :: n_external_pot_par
!!$
!!$  namelist /disper_repal/ &
!!$       n_differ_solv_atoms, &
!!$       n_external_pot_par, &
!!$       use_rappe_data
!!$
!!$  integer(kind=i4_kind) :: df_n_differ_solv_atoms = 0
!!$  integer(kind=i4_kind) :: df_n_external_pot_par = 0
!!$  logical :: df_use_rappe_data = .false.
!!$
!!$  character*2 :: name_solv_atom
!!$  integer(kind=i4_kind) :: n_solv_atom
  character*2, allocatable :: nm_slv_at(:)
  integer(kind=i4_kind), allocatable :: n_slv_at(:)
!!$
!!$  namelist /solvent_atom/ &
!!$       name_solv_atom, &
!!$       n_solv_atom


!!$  character*2 :: atom
!!$  real(kind=r8_kind) :: K_atom
!!$  real(kind=r8_kind) :: R0_atom
  character*2, allocatable :: atm(:)
  real(kind=r8_kind), allocatable :: K_at(:)
  real(kind=r8_kind), allocatable :: R0_at(:)

!!$  namelist /potential_param/ &
!!$       atom, &
!!$       K_atom, &
!!$       R0_atom
!!$
!!$  logical :: yes_no = .false.

  real(kind=r8_kind) , allocatable :: AA(:,:),CC(:,:),zeta_av(:,:)
  real(kind=r8_kind) , parameter :: A_p=0.214_r8_kind ! kcal/mol
  real(kind=r8_kind) , parameter :: C_p=47000.0_r8_kind ! kcal/mol
  real(kind=r8_kind) , parameter :: alpha_p=12.35_r8_kind


  external error_handler

contains

  !********************************************************************
  subroutine disp_rep_read
  ! read in namelist disper_repal and if not default solvent (water)
  ! read in namelists solvent_atom to specify the solvent
  !** End of interface *****************************************

    integer :: status

    ndsa=2
    if(.not. allocated(nm_slv_at)) then
       allocate(nm_slv_at(2),n_slv_at(2),stat=status)
       if ( status /= 0) call error_handler( &
            "disp_rep_read: allocation of nm_slv_at is failed")
    endif
    nm_slv_at(1)="H_"
    n_slv_at(1)=2
    nm_slv_at(2)="O_"
    n_slv_at(2)=1       

  end subroutine disp_rep_read
  !*********************************************************

  !*********************************************************
  subroutine initialize_data(skip_short)

  ! set the default () values for the dispersion repulsion data
  ! or read in user defined parameters if given for an (unique) atom
  !** End of interface *****************************************

    logical, intent(inout) :: skip_short(:)
    integer(kind=i4_kind) :: i,j,status,typ

    !initialization of solvent atoms
    allocate(K_solv(ndsa) , &
         R_solv(ndsa)     , &
         r_sphere_at(ndsa), &
         stat=status)
       if ( status /= 0) call error_handler( &
            "MolMech:initialize_data: allocation of K_solv,R_solv,r_sphere_at are failed")

    m1: do i=1,ndsa
       do j=1,98
          if(nm_slv_at(i)==atom_data(j)%name) then
             K_solv(i)=atom_data(j)%k_dr
             R_solv(i)=atom_data(j)%r0_dr
             r_sphere_at(i)=atom_data(j)%r0_dr
             cycle m1
          endif
       enddo
    enddo m1

    !initialization of solute atoms
    allocate(K_solu(n_species),R_solu(n_species),stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:initialize_data: allocation of K_solu,R_solu are failed")

    do i=1,n_species
       typ=atoms_cart(i)%type
       nm_solute_at(i)=atoms(typ)%name
       K_solu(i)=atoms(typ)%k_dr
       R_solu(i)=atoms(typ)%r0_dr
    enddo

    skip_short(1:n_species) = .false.
    do i=1,n_species
        if(R_solu(i)<1.e-10) skip_short(i)=.true.
    enddo

    allocate(AA(ndsa,n_species),stat=status)
    if(status/=0) call error_handler("MolMech:dis_rep initialize_data, alloc (1) failed")
    allocate(CC(ndsa,n_species),stat=status)
    if(status/=0) call error_handler("MolMech:dis_rep initialize_data, alloc (2) failed")
    allocate(zeta_av(ndsa,n_species),stat=status)
    if(status/=0) call error_handler("MolMech:dis_rep initialize_data, alloc (3) failed")

    do i=1,ndsa
       AA(i,:)=A_p*c2J
       CC(i,:)=C_p*c2J
       zeta_av(i,:)=alpha_p
    enddo

  end subroutine initialize_data
  !*********************************************************

  !*********************************************************
!!$  subroutine get_rappe_data(with_pc)
!!$  ! do the same as initialize_data but use the data from
!!$  ! Rappe as default.
!!$  ! The different exponential factors in the Rappe data set
!!$  ! are averaged geometrical for each pair of atoms.
!!$  !** End of interface *****************************************
!!$
!!$    use unique_atom_module, only: N_unique_atoms,unique_atoms
!!$    use pointcharge_module, only: pointcharge_N,pointcharge_array
!!$
!!$    logical :: with_pc
!!$    integer(kind=i4_kind) :: i,j,k,ind,status
!!$    integer :: N_pc
!!$    real(kind=r8_kind), allocatable :: z_solv(:),z_solu(:)
!!$
!!$    !initialization of solvent atoms
!!$    allocate(K_solv(ndsa) , &
!!$         R_solv(ndsa)     , &
!!$         z_solv(ndsa)     , &
!!$         r_sphere_at(ndsa), &
!!$         stat=status)
!!$       if ( status /= 0) call error_handler( &
!!$            "initialize_data: allocation of K_solv,R_solv,z_solv, r_sphere_at are failed")
!!$
!!$    m1: do i=1,ndsa
!!$       do k=1,n_external_pot_par
!!$          if(nm_slv_at(i)==atm(k)) then
!!$             K_solv(i)=K_at(k)
!!$             R_solv(i)=R0_at(k)/ang_au       !convert to atomic units
!!$             r_sphere_at(i)=R0_at(k)/ang_au  !convert to atomic units
!!$             z_solv(i)=alpha_p
!!$             cycle m1
!!$          endif
!!$       enddo
!!$       do j=1,98
!!$          if(nm_slv_at(i)==atom_name(j)) then
!!$             K_solv(i)=sqrt(D_def_rap(j)/au_kcalmol)
!!$             R_solv(i)=R_def_rap(j)/(2.0_r8_kind*ang_au)     !convert to atomic units
!!$             z_solv(i)=zeta_def_rap(j)
!!$             r_sphere_at(i)=R_solv(i)
!!$             cycle m1
!!$          endif
!!$       enddo
!!$    enddo m1
!!$
!!$    !initialization of solute atoms
!!$    N_pc=0
!!$    if(with_pc) N_pc=pointcharge_N
!!$
!!$    do i=1,N_unique_atoms
!!$       ind=int(unique_atoms(i)%Z)
!!$       nm_solute_at(i)=atom_name(ind)
!!$    enddo
!!$    if(with_pc) then
!!$       do i=1,pointcharge_N
!!$          nm_solute_at(i+N_unique_atoms)=trim(pointcharge_array(i)%name)
!!$       end do
!!$    end if
!!$
!!$    allocate(K_solu(N_unique_atoms+N_pc), &
!!$         R_solu(N_unique_atoms+N_pc)    , &
!!$      z_solu(N_unique_atoms+N_pc)    , &
!!$         stat=status)
!!$    if ( status /= 0) call error_handler( &
!!$         "initialize_data: allocation of K_solu,R_solu, z_solu are failed")
!!$
!!$    m2: do i=1,N_unique_atoms+N_pc
!!$       do k=1,n_external_pot_par
!!$          if(nm_solute_at(i)==atm(k)) then
!!$             K_solu(i)=K_at(k)
!!$             R_solu(i)=R0_at(k)/ang_au    !convert to atomic units
!!$          z_solu(i)=alpha_p
!!$             cycle m2
!!$          endif
!!$       enddo
!!$       do j=1,98
!!$          if(nm_solute_at(i)==atom_name(j)) then
!!$             K_solu(i)=sqrt(D_def_rap(j)/au_kcalmol)
!!$             R_solu(i)=R_def_rap(j)/(2.0_r8_kind*ang_au)     !convert to atomic units
!!$             z_solu(i)=zeta_def_rap(j)
!!$             cycle m2
!!$          endif
!!$       enddo
!!$    enddo m2      
!!$
!!$    allocate(AA(ndsa,N_unique_atoms+N_pc),stat=status)
!!$    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (1) failed")
!!$    allocate(CC(ndsa,N_unique_atoms+N_pc),stat=status)
!!$    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (2) failed")
!!$    allocate(zeta_av(ndsa,N_unique_atoms+N_pc),stat=status)
!!$    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (3) failed")
!!$
!!$    do i=1,ndsa
!!$       do j=1,N_unique_atoms+N_pc
!!$          zeta_av(i,j)=sqrt(z_solu(j)*z_solv(i))
!!$          AA(i,j)=zeta_av(i,j)/(zeta_av(i,j)-6.0_r8_kind)
!!$          CC(i,j)=6.0_r8_kind/(zeta_av(i,j)-6.0_r8_kind)*exp(zeta_av(i,j))
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine get_rappe_data
  !*********************************************************

  !*********************************************************
  subroutine disp_rep_energies_and_grad(slv_cycle,do_energy,do_grad,gradient,ma)
    ! calculate dispersion and repulsion energies or their
    ! gradients, respectively.
    ! Called for each type of atom in the solvent molecule
    ! (with an other cavity each time). Based on boundary element
    ! method. The parametrized pair potential of both, dispersion
    ! and repulsion,has the form
    ! V_ij=-A /r_ij**6 + B exp(-zeta*r_ij)
    ! the solvent atoms are assumed to be distribouted homogenously
    ! outside the cavity.
    ! If parameters are not available the contribution is skipped.

    integer(kind=i4_kind) :: slv_cycle
    logical :: do_grad,do_energy
    real(kind=r8_kind),optional :: gradient(3)
    integer(kind=i4_kind),optional :: ma
    !** End of interface *****************************************

    real(kind=r8_kind) :: D_solv_solu,B_solv_solu,gamma_solv_solu,g_r
    real(kind=r8_kind) :: disp_p, rep_p, R_abs, R_xyz(3),dot_prod
    real(kind=r8_kind):: grad_disp,grad_rep,grad_abs,grad_p_no
    integer :: nslv
    integer(kind=i4_kind) :: i,k,l2
    real(r8_kind) :: xyz(3)


    nslv=n_slv_at(slv_cycle)
    if(K_solv(slv_cycle)==0.0_r8_kind .or. R_solv(slv_cycle)==0.0_r8_kind) then
       write(output_device,*) "------------------------------------------------------"
       write(output_device,*) "Solvent effect: Dispersion-repulsion term: WARNING!!!"
       write(output_device,*) "The contribution of the solvent atom "//nm_slv_at(slv_cycle)//" has been skiped,"
       write(output_device,*) "because of either K or R0 parameter of this atom is "
       write(output_device,*) "equal ZERO"
       write(output_device,*) "------------------------------------------------------"
       return
    endif

    m1: do i=1,n_species
       if(K_solu(i)==0.0_r8_kind .or. R_solu(i)==0.0_r8_kind) then
          write(output_device,*) "------------------------------------------------------"
          write(output_device,*) "Solvent effect: Dispersion-repulsion term: WARNING!!!"
          write(output_device,*) "The contribution of the solute atom "//nm_solute_at(i)//" has been skiped,"
          write(output_device,*) "because of either K or R0 parameter of this atom is "
          write(output_device,*) "equal ZERO"
          write(output_device,*) "------------------------------------------------------"
          cycle m1
       endif
       D_solv_solu=K_solv(slv_cycle)*K_solu(i)*AA(slv_cycle,i)* &
            64.0_r8_kind*(R_solv(slv_cycle)*R_solu(i))**3
       B_solv_solu=K_solv(slv_cycle)*K_solu(i)*CC(slv_cycle,i)

       gamma_solv_solu=zeta_av(slv_cycle,i)/&
                (two*sqrt(R_solv(slv_cycle)*R_solu(i)))

       xyz(:)=atoms_cart(i)%r
       m3: do k=1,N_points_dr

          R_xyz(:)=D_R_points(k)%position(:)-xyz(:)
          R_abs=sqrt(dot_product(R_xyz,R_xyz))
          disp_p=-D_solv_solu/(three*R_abs**6)
          g_r=gamma_solv_solu*R_abs
          rep_p=B_solv_solu*exp(-g_r)*(one/g_r+two/g_r**2+two/g_r**3)
          dot_prod=dot_product(R_xyz,D_R_points(k)%out_normal(:))

          if(do_energy) E_disp=E_disp+nslv*num_dens_dr*disp_p*D_R_points(k)%area*dot_prod

          if(do_energy) E_rep=E_rep+nslv*num_dens_dr*rep_p*D_R_points(k)%area*dot_prod

          if(do_grad) then
             call calc_grads()
          endif
       enddo m3
    enddo m1

  contains

    subroutine calc_grads()

       do l2=1,3
         grad_abs=dot_product(R_xyz,D_R_grad(k)%distance(l2,ma,i,:))/R_abs

         grad_disp=two*D_solv_solu/(R_abs**7)*grad_abs
         
         grad_rep=gamma_solv_solu*grad_abs*(&
              -rep_p+ B_solv_solu*exp(-g_r)*(&
              -one/g_r**2-four/g_r**3-six/g_r**4))

         grad_p_no=&
              dot_product(R_xyz,&
              D_R_grad(k)%out_normal(l2,ma,:))&
              +dot_product(D_R_grad(k)%distance(l2,ma,i,:),&
              D_R_points(k)%out_normal(:))

         gradient(l2)=gradient(l2)+&
              nslv*num_dens_dr*((grad_disp+grad_rep)*&
              D_R_points(k)%area*dot_prod+(disp_p+rep_p)*&
              D_R_grad(k)%area(l2,ma)*dot_prod+(disp_p+rep_p)*&
              D_R_points(k)%area*grad_p_no)
      enddo
    end subroutine calc_grads
  end subroutine disp_rep_energies_and_grad
  !*********************************************************
 
  !*********************************************************
  subroutine deallocate_D_R_points(do_grad)
    !** End of interface *****************************************

    logical, intent(in) :: do_grad
    integer(kind=i4_kind) :: i,status

    if(do_grad) then
       do i=1,N_points_dr
          deallocate(D_R_grad(i)%position,D_R_grad(i)%out_normal,&
               D_R_grad(i)%area,D_R_grad(i)%distance,stat=status)
          if ( status /= 0) call error_handler( &
               "MolMech:deallocate_D_R_points: deallocation of D_R_grad%position failed")
       enddo
    endif

    deallocate(D_R_points,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:deallocate_D_R_points: deallocation of D_R_points failed")

    if(do_grad) then
       deallocate(D_R_grad,stat=status)
       if ( status /= 0) call error_handler( &
            "MolMech:deallocate_D_R_points: deallocation of D_R_grad failed")
    endif

  end subroutine deallocate_D_R_points
  !*********************************************************

  !*********************************************************
  subroutine deallocate_parameters(do_grad)
    !** End of interface *****************************************

    logical, intent(in) :: do_grad
    integer(kind=i4_kind) :: status

    deallocate(K_solv,R_solv,r_sphere_at,K_solu,R_solu,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:deallocate_parameters: deallocation of K_solv is failed")

    deallocate(AA,CC,zeta_av,stat=status)
    if(status/=0) call error_handler("MolMech:deallocate_parameters, dealloc (1) failed")

   if(do_grad) then
      deallocate(nm_slv_at,n_slv_at,stat=status)
      if ( status /= 0) call error_handler( &
         "MolMech:disp_rep_write: deallocation of nm_slv_at is failed")
   endif

  end subroutine deallocate_parameters

end module vdwcm_module

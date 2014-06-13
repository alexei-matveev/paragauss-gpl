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
module disp_rep_module
!
!  Author: AS,MF
!  Date: 4/5 2000
!

# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use iounitadmin_module, only: output_unit,write_to_output_units
  use integralpar_module, only: integralpar_2dervs
  use atoms_data_module

  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public disp_rep_read, disp_rep_write, check_read_dr, initialize_data, &
       disp_rep_energies_and_grad, deallocate_D_R_points,&
       deallocate_parameters, &
       get_rappe_data, get_k_and_r


  ! data type to collect unique surface points, normals and areas
  ! of a cavity tessellation.
  ! idices: (coordinates,number_of_equal_partner)
  type, public :: disp_rep_points
     ! Number of partners of unique point charge
     integer(kind=i4_kind)                         :: N_equal_points
     ! center of tessera
     real(kind=r8_kind), pointer                :: position(:,:)
     ! outward normal in tessera center
     real(kind=r8_kind), pointer                :: out_normal(:,:)
     ! area of tessera
     real(kind=r8_kind)                :: area
  end type disp_rep_points

  ! type to collect the gradients of the surface point positions,
  ! the tessera normals, the areas and the distance to all sphere
  ! centers.
  type, public :: disp_rep_points_grad
     real(kind=r8_kind), pointer                :: position(:,:,:,:,:)
     real(kind=r8_kind), pointer                :: out_normal(:,:,:,:,:)
     real(kind=r8_kind), pointer                :: distance(:,:,:,:,:,:,:)
     real(kind=r8_kind), pointer                :: area(:,:,:,:)
  end type disp_rep_points_grad

  ! type to collect the second derivatives of
  ! the tessera normals, the areas and the distance to all sphere
  ! centers.
  type, public :: disp_rep_second_deriv
     type(arrmat6), pointer :: out_normal(:,:)
     type(arrmat6), pointer :: distance(:,:)
     type(arrmat6), pointer :: area(:)
  end type disp_rep_second_deriv

  !
  ! Energy  contributions form  dispersion and  repulsion interaction.
  ! FIXME: these  energies are computed (incremented)  in this module.
  ! However they are reset  to zero from outside, see disp_rep_wrap(),
  ! so that they cannot be made "protected":
  !
  real (r8_kind), public :: E_disp, E_rep

  ! number of unique surface points of the actual cavity
  integer(kind=i4_kind), public :: N_points_dr
  ! position/normals/areas of surface tesserae of the actual cavity
  ! for dispersion/repulsion contribution
  type(disp_rep_points), allocatable, target, public :: D_R_points(:)
  ! gradients of position/normals/areas of surface tesserae of the
  ! actual cavity
  type(disp_rep_points_grad), allocatable, target, public :: D_R_grad(:)
  ! 2nd derivatives of surface tesserae of the actual cavity
  type(disp_rep_second_deriv), allocatable, target, public :: D_R_hess(:)

  ! Names of the solute atom (atoms of the molecule)
  character(len=2), public, allocatable :: nm_solute_at(:)
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

  integer(kind=i4_kind) :: n_differ_solv_atoms
  integer(kind=i4_kind) :: n_external_pot_par

  namelist /disper_repal/ &
       n_differ_solv_atoms, &
       n_external_pot_par, &
       use_rappe_data

  integer(kind=i4_kind) :: df_n_differ_solv_atoms = 0
  integer(kind=i4_kind) :: df_n_external_pot_par = 0
  logical :: df_use_rappe_data = .false.

  character(len=2) :: name_solv_atom
  integer(kind=i4_kind) :: n_solv_atom
  character(len=2), allocatable :: nm_slv_at(:)
  integer(kind=i4_kind), allocatable :: n_slv_at(:)

  namelist /solvent_atom/ &
       name_solv_atom, &
       n_solv_atom


  character(len=2) :: atom
  real(kind=r8_kind) :: K_atom
  real(kind=r8_kind) :: R0_atom
  character(len=2), allocatable :: atm(:)
  real(kind=r8_kind), allocatable :: K_at(:)
  real(kind=r8_kind), allocatable :: R0_at(:)

  namelist /potential_param/ &
       atom, &
       K_atom, &
       R0_atom

  logical :: yes_no = .false.

  real(kind=r8_kind) , allocatable :: AA(:,:),CC(:,:),zeta_av(:,:)

  ! FIXME: at other places in PG we use 0.52917706 instead, see module
  ! constants:
  real(kind=r8_kind) , parameter :: ang_au = 0.529177249_r8_kind
  real(kind=r8_kind) , parameter :: A_p = 0.214_r8_kind   ! kcal/mol
  real(kind=r8_kind) , parameter :: C_p = 47000.0_r8_kind ! kcal/mol
  real(kind=r8_kind) , parameter :: alpha_p=12.35_r8_kind


  external error_handler

contains

#define DOT3(u, v)      dot_product(u, v)

#ifndef DOT3
  pure function DOT3(u, v) result(uv)
    !
    ! scalar product for 3-vectors, hope for inlining
    !
    implicit none
    real(r8_kind), intent(in) :: u(:)
    real(r8_kind), intent(in) :: v(:)
    real(r8_kind)             :: uv ! result
    ! *** end of interface ***

    uv = u(1) * v(1) + u(2) * v(2) + u(3) * v(3)
  end function DOT3
#endif

#define LEN3(v)         sqrt(DOT3(v, v))

#ifndef LEN3
  pure function LEN3(v) result(r)
    !
    ! scalar product for 3-vectors, hope for inlining
    !
    implicit none
    real(r8_kind), intent(in) :: v(:)
    real(r8_kind)             :: r ! result
    ! *** end of interface ***

    r = SQRT(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
  end function LEN3
#endif

  !********************************************************************
  subroutine disp_rep_read
  ! read in namelist disper_repal and if not default solvent (water)
  ! read in namelists solvent_atom to specify the solvent
  !** End of interface *****************************************

    use input_module

    integer(kind=i4_kind) :: unit,status,i,j

    n_differ_solv_atoms=df_n_differ_solv_atoms
    n_external_pot_par=df_n_external_pot_par
    use_rappe_data=df_use_rappe_data

    unit = input_intermediate_unit()
    call input_read_to_intermediate
    read(unit, nml=disper_repal, iostat=status)
    if (status .ne. 0) call input_error( &
         "disp_rep_read: namelist disper_repal.")

    yes_no = .true.

    if(n_differ_solv_atoms<=0) then
       ndsa=2
       if(.not. allocated(nm_slv_at)) then
        allocate(nm_slv_at(2),n_slv_at(2),stat=status)
        if ( status /= 0) call error_handler( &
            "disp_rep_read: allocation of nm_slv_at is failed")
       endif
       nm_slv_at(1)="H "
       n_slv_at(1)=2
       nm_slv_at(2)="O "
       n_slv_at(2)=1
    else
       ndsa=n_differ_solv_atoms
       if(.not. allocated(nm_slv_at)) then
         allocate(nm_slv_at(ndsa),n_slv_at(ndsa),stat=status)
         if ( status /= 0) call error_handler( &
            "disp_rep_read: allocation of nm_slv_at is failed(1)")
       endif
       do i=1,ndsa
          call input_read_to_intermediate()
          read(unit, nml=solvent_atom, iostat=status)
          if (status .ne. 0) call input_error( &
               "disp_rep_read: namelist solvent_atom.")
          do j=1,98
             if(name_solv_atom==atom_name(j)) exit
             if (j==98) call input_error( &
                  "You have specified the wrong name of atom in SOLVENT_ATOM namelist")
          enddo
          nm_slv_at(i)=name_solv_atom
          n_slv_at(i)=n_solv_atom
       enddo
    endif

    if(n_external_pot_par > 0) then
       if(.not. allocated(atm)) then
         allocate(atm(n_external_pot_par),K_at(n_external_pot_par), &
            R0_at(n_external_pot_par),stat=status)
         if ( status /= 0) call error_handler( &
            "disp_rep_read: allocation of atm, K_at, R0_at are failed")
       endif
       do i=1,n_external_pot_par
          call input_read_to_intermediate()
          read(unit, nml=potential_param, iostat=status)
          if (status .ne. 0) call input_error( &
               "disp_rep_read: namelist potential_param.")
          do j=1,98
             if(atom==atom_name(j)) exit
             if (j==98) call input_error( &
                  "You have specified the wrong name of atom in POTENTIAL_PARAM namelist")
          enddo
          atm(i)=atom
          K_at(i)=K_atom
          R0_at(i)=R0_atom
       enddo
    endif

  end subroutine disp_rep_read
  !*********************************************************

  !*********************************************************
  subroutine check_read_dr
  ! set defaults if no namlist disper_repal is present
  !** End of interface *****************************************


    integer(kind=i4_kind) :: status

    if (.not.yes_no) then
       n_differ_solv_atoms=df_n_differ_solv_atoms
       n_external_pot_par=df_n_external_pot_par
       use_rappe_data=df_use_rappe_data
       ndsa=2
       if(.not. allocated(nm_slv_at)) then
         allocate(nm_slv_at(2),n_slv_at(2),stat=status)
         if ( status /= 0) call error_handler( &
            "check_read_dr: allocation of nm_slv_at is failed")
       endif
       nm_slv_at(1)="H "
       n_slv_at(1)=2
       nm_slv_at(2)="O "
       n_slv_at(2)=1
    endif

  end subroutine check_read_dr
  !********************************************************

  !*********************************************************
  subroutine disp_rep_write(unit)
    !
    ! Write input/default values of namelist disper_repal and if given
    ! by the user also external potential parameters.
    !
    use echo_input_module, only: start, real, flag, intg, word, stop, &
         echo_level_full, word_format
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(kind=i4_kind), intent(in) :: unit
    ! *** end of interface ***

    integer(kind=i4_kind) :: i,status

    word_format = '("    ",a," = ",a7  :" # ",a)'

    call start("DISPER_REPAL","DISP_REP_WRITE",unit,operations_echo_input_level)
    call intg("N_DIFFER_SOLV_ATOMS  ",n_differ_solv_atoms  ,df_n_differ_solv_atoms  )
    call intg("N_EXTERNAL_POT_PAR   ",n_external_pot_par   ,df_n_external_pot_par   )
    call flag("USE_RAPPE_DATA       ",use_rappe_data       ,df_use_rappe_data       )
    call stop()

    do i=1,n_differ_solv_atoms
       call start("SOLVENT_ATOM","DISP_REP_WRITE",unit,operations_echo_input_level)
       call word("NAME_SOLV_ATOM       ",nm_slv_at(i)         ,"  "                    )
       call intg("N_SOLV_ATOM          ",n_slv_at(i)          ,0                       )
       call stop()
    enddo

    do i=1,n_external_pot_par
       call start("POTENTIAL_PARAM","DISP_REP_WRITE",unit,operations_echo_input_level)
       call word("ATOM                 ",atm(i)               ,"  "                    )
       call real("K_ATOM               ",K_at(i)              ,0.0_r8_kind             )
       call real("R0_ATOM              ",R0_at(i)             ,0.0_r8_kind             )
       call stop()
    enddo

!   deallocate(nm_slv_at,n_slv_at,stat=status)
!   if ( status /= 0) call error_handler( &
!        "disp_rep_write: deallocation of nm_slv_at is failed")

    if(n_external_pot_par > 0) then
       deallocate(atm,K_at,R0_at,stat=status)
       if ( status /= 0) call error_handler( &
            "disp_rep_write: deallocation of atm is failed")
    endif

  end subroutine disp_rep_write
  !*********************************************************

  !*********************************************************
  subroutine initialize_data(skip_short,with_pc)
    !
    ! Set the default  () values for the dispersion  repulsion data or
    ! read in user defined parameters if given for an (unique) atom
    !
    use constants, only: kcal ! to convert A_p and C_p
    use unique_atom_module, only: N_unique_atoms,unique_atoms
    use pointcharge_module, only: pointcharge_N,pointcharge_array
    implicit none
    ! *** end of interface ***

    logical, intent(inout) :: skip_short(:)
    logical :: with_pc
    integer :: N_pc
    integer(kind=i4_kind) :: i,j,k,ind,status

    !initialization of solvent atoms
    allocate(K_solv(ndsa) , &
         R_solv(ndsa)     , &
         r_sphere_at(ndsa), &
         stat=status)
    if ( status /= 0) call error_handler( &
         "initialize_data: allocation of K_solv,R_solv,r_sphere_at are failed")
    K_solv=0.0_r8_kind; R_solv=0.0_r8_kind; r_sphere_at=0.0_r8_kind

    m1: do i=1,ndsa
       do k=1,n_external_pot_par
          if(nm_slv_at(i)==atm(k)) then
             K_solv(i)=K_at(k)
             R_solv(i)=R0_at(k)/ang_au       !convert to atomic units
             r_sphere_at(i)=R0_at(k)/ang_au  !convert to atomic units
             cycle m1
          endif
       enddo
       do j=1,98
          if(nm_slv_at(i)==atom_name(j)) then
             K_solv(i)=K_def(j)
             R_solv(i)=R0_def(j)/ang_au      !convert to atomic units
             r_sphere_at(i)=R0_def(j)/ang_au !convert to atomic units
             cycle m1
          endif
       enddo
       if(K_solv(i)==0.0_r8_kind .or.R_solv(i)==0.0_r8_kind) then
          call write_to_output_units( &
               "------------------------------------------------------")
          call write_to_output_units( &
               "Solvent effect: Dispersion-repulsion term: WARNING!!!")
          call write_to_output_units( &
               "The contribution of the solvent atom (pc) "//nm_slv_at(i)//" has been skiped,")
          call write_to_output_units( &
               "because of either K or R0 parameter of this atom is ")
          call write_to_output_units("equal ZERO")
       endif
    enddo m1

    !initialization of solute atoms
    N_pc=0
    if(with_pc) N_pc=pointcharge_N

    do i=1,N_unique_atoms
       ind=int(unique_atoms(i)%Z)
       nm_solute_at(i)=atom_name(ind)
    enddo
    if(with_pc) then
       do i=1,pointcharge_N
          nm_solute_at(i+N_unique_atoms)=trim(pointcharge_array(i)%name)
       end do
    end if

    allocate(K_solu(N_unique_atoms+N_pc), &
         R_solu(N_unique_atoms+N_pc)    , &
         stat=status)
    if ( status /= 0) call error_handler( &
         "initialize_data: allocation of K_solu,R_solu are failed")
    K_solu=0.0_r8_kind; R_solu=0.0_r8_kind

    m2: do i=1,N_unique_atoms+N_pc
       do k=1,n_external_pot_par
          if(nm_solute_at(i)==atm(k)) then
             K_solu(i)=K_at(k)
             R_solu(i)=R0_at(k)/ang_au    !convert to atomic units
             cycle m2
          endif
       enddo
       do j=1,98
          if(nm_solute_at(i)==atom_name(j)) then
             K_solu(i)=K_def(j)
             R_solu(i)=R0_def(j)/ang_au   !convert to atomic units
             cycle m2
          endif
       enddo
       if(K_solu(i)==0.0_r8_kind .or. R_solu(i)==0.0_r8_kind) then
          call write_to_output_units( &
               "------------------------------------------------------")
          call write_to_output_units( &
               "Solvent effect: Dispersion-repulsion term: WARNING!!!")
          call write_to_output_units( &
               "The contribution of the solute atom "//nm_solute_at(i)//" has been skiped,")
          call write_to_output_units( &
               "because of either K or R0 parameter of this atom is ")
          call write_to_output_units("equal ZERO")
       endif
    enddo m2

    skip_short(1:N_unique_atoms+N_pc) = .false.
    do i=1,N_unique_atoms+N_pc
        if(R_solu(i)<1.e-10) skip_short(i)=.true.
    enddo

    allocate(AA(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep initialize_data, alloc (1) failed")
    allocate(CC(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep initialize_data, alloc (2) failed")
    allocate(zeta_av(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep initialize_data, alloc (3) failed")

    ! A_p and C_p global module parameters are in kcal/mol, convert to
    ! atomic units:
    do i = 1, ndsa
       AA(i, :) = A_p * kcal
       CC(i, :) = C_p * kcal
       zeta_av(i,:)=alpha_p
    enddo

  end subroutine initialize_data
  !*********************************************************

  !*********************************************************
  subroutine get_rappe_data(with_pc)
    !
    ! Do the same as initialize_data()  but use the data from Rappe as
    ! default.  The  different exponential  factors in the  Rappe data
    ! set are averaged geometrical for each pair of atoms.
    !
    use constants, only: kcal   ! D_def_rap(:) is in kcals
    use unique_atom_module, only: N_unique_atoms,unique_atoms
    use pointcharge_module, only: pointcharge_N,pointcharge_array
    implicit none
    ! *** end of interface ***

    logical :: with_pc
    integer(kind=i4_kind) :: i,j,k,ind,status
    integer :: N_pc
    real(kind=r8_kind), allocatable :: z_solv(:),z_solu(:)

    !initialization of solvent atoms
    allocate(K_solv(ndsa) , &
         R_solv(ndsa)     , &
         z_solv(ndsa)     , &
         r_sphere_at(ndsa), &
         stat=status)
    if ( status /= 0) call error_handler( &
         "initialize_data: allocation of K_solv,R_solv,z_solv, r_sphere_at are failed")
    K_solv=0.0_r8_kind; R_solv=0.0_r8_kind; r_sphere_at=0.0_r8_kind

    m1: do i=1,ndsa
       do k=1,n_external_pot_par
          if(nm_slv_at(i)==atm(k)) then
             K_solv(i)=K_at(k)
             R_solv(i)=R0_at(k)/ang_au       !convert to atomic units
             r_sphere_at(i)=R0_at(k)/ang_au  !convert to atomic units
             z_solv(i)=alpha_p
             cycle m1
          endif
       enddo

       ! The data  in D_def_rap(:) array  (see atoms_data_module) must
       ! be in kcal/mol, I assume:
       do j=1,98
          if(nm_slv_at(i)==atom_name(j)) then
             K_solv(i) = sqrt (D_def_rap(j) * kcal)
             R_solv(i)=R_def_rap(j)/(2.0_r8_kind*ang_au)     !convert to atomic units
             z_solv(i)=zeta_def_rap(j)
             r_sphere_at(i)=R_solv(i)
             cycle m1
          endif
       enddo
       if(K_solv(i)==0.0_r8_kind .or.R_solv(i)==0.0_r8_kind) then
          call write_to_output_units( &
               "------------------------------------------------------")
          call write_to_output_units( &
               "Solvent effect: Dispersion-repulsion term: WARNING!!!")
          call write_to_output_units( &
               "The contribution of the solvent atom (pc) "//nm_slv_at(i)//" has been skiped,")
          call write_to_output_units( &
               "because of either K or R0 parameter of this atom is ")
          call write_to_output_units("equal ZERO")
       endif
    enddo m1

    !initialization of solute atoms
    N_pc=0
    if(with_pc) N_pc=pointcharge_N

    do i=1,N_unique_atoms
       ind=int(unique_atoms(i)%Z)
       nm_solute_at(i)=atom_name(ind)
    enddo
    if(with_pc) then
       do i=1,pointcharge_N
          nm_solute_at(i+N_unique_atoms)=trim(pointcharge_array(i)%name)
       end do
    end if

    allocate(K_solu(N_unique_atoms+N_pc), &
         R_solu(N_unique_atoms+N_pc)    , &
         z_solu(N_unique_atoms+N_pc)    , &
         stat=status)
    if ( status /= 0) call error_handler( &
         "initialize_data: allocation of K_solu,R_solu, z_solu are failed")
    K_solu=0.0_r8_kind; R_solu=0.0_r8_kind

    m2: do i=1,N_unique_atoms+N_pc
       do k=1,n_external_pot_par
          if(nm_solute_at(i)==atm(k)) then
             K_solu(i)=K_at(k)
             R_solu(i)=R0_at(k)/ang_au    !convert to atomic units
             z_solu(i)=alpha_p
             cycle m2
          endif
       enddo

       ! The data  in D_def_rap(:) array  (see atoms_data_module) must
       ! be in kcal/mol, I assume:
       do j=1,98
          if(nm_solute_at(i)==atom_name(j)) then
             K_solu(i) = sqrt (D_def_rap(j) * kcal)
             R_solu(i)=R_def_rap(j)/(2.0_r8_kind*ang_au)     !convert to atomic units
             z_solu(i)=zeta_def_rap(j)
             cycle m2
          endif
       enddo
       if(K_solu(i)==0.0_r8_kind .or. R_solu(i)==0.0_r8_kind) then
          call write_to_output_units( &
               "------------------------------------------------------")
          call write_to_output_units( &
               "Solvent effect: Dispersion-repulsion term: WARNING!!!")
          call write_to_output_units( &
               "The contribution of the solute atom "//nm_solute_at(i)//" has been skiped,")
          call write_to_output_units( &
               "because of either K or R0 parameter of this atom is ")
          call write_to_output_units("equal ZERO")
       endif
    enddo m2

    allocate(AA(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (1) failed")
    allocate(CC(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (2) failed")
    allocate(zeta_av(ndsa,N_unique_atoms+N_pc),stat=status)
    if(status/=0) call error_handler("dis_rep get_rappe_data, alloc (3) failed")

    do i=1,ndsa
       do j=1,N_unique_atoms+N_pc
          zeta_av(i,j)=sqrt(z_solu(j)*z_solv(i))
          AA(i,j)=zeta_av(i,j)/(zeta_av(i,j)-6.0_r8_kind)
          CC(i,j)=6.0_r8_kind/(zeta_av(i,j)-6.0_r8_kind)*exp(zeta_av(i,j))
       enddo
    enddo

  end subroutine get_rappe_data
  !*********************************************************

  !*********************************************************
  subroutine get_k_and_r(slv_cycle,K_slv,R_slv,Name_slv)
    integer(i4_kind), intent(in) :: slv_cycle
    real(r8_kind), intent(out) :: K_slv,R_slv
    character(len=2), intent(out) :: Name_slv

    K_slv=K_solv(slv_cycle)
    R_slv=R_solv(slv_cycle)
    Name_slv=nm_slv_at(slv_cycle)

  end subroutine get_k_and_r
  !*********************************************************

  !*********************************************************
  subroutine disp_rep_energies_and_grad (slv_cycle, do_grad, with_pc, gradient, ma, ea, hessian, ma1, ea1)
    ! calculate dispersion and repulsion energies or their
    ! gradients, respectively.
    ! Called for each type of atom in the solvent molecule
    ! (with an other cavity each time). Based on boundary element
    ! method. The parametrized pair potential of both, dispersion
    ! and repulsion,has the form
    ! V_ij=-A /r_ij**6 + B exp(-zeta*r_ij)
    ! the solvent atoms are assumed to be distributed homogenously
    ! outside the cavity.
    ! If parameters are not available the contribution is skipped.
    use unique_atom_module, only: N_unique_atoms,unique_atoms
    use pointcharge_module, only: pointcharge_N,pointcharge_array

    integer(kind=i4_kind), intent(in) :: slv_cycle
    logical              , intent(in) :: do_grad
    logical :: with_pc
    real(kind=r8_kind),intent(inout),optional :: gradient(3)
    integer(kind=i4_kind), intent(in),optional :: ma,ea
    real(kind=r8_kind),intent(inout),optional :: hessian(3,3)
    integer(kind=i4_kind), intent(in),optional :: ma1,ea1
    !** End of interface *****************************************

    real(kind=r8_kind) :: D_solv_solu,B_solv_solu,gamma_solv_solu,g_r
    real(kind=r8_kind) :: disp_p, rep_p, R_abs, R_xyz(3), g_rep1, g_rep2
    real(kind=r8_kind) :: grad_disp,grad_rep,grad_abs,grad_p_no
    real(r8_kind)      :: grad_disp1,grad_rep1,grad_abs1,grad_p_no1
    real(r8_kind)      :: hess_disp,hess_rep,hess_abs,hess_p_no

    integer(kind=i4_kind) :: i,j,k,l,l2,l3,m
    integer :: N_pc,n_equal
    real(r8_kind) :: xyz(3)

    N_pc=0
    if(with_pc) N_pc=pointcharge_N

    DPRINT 'disp_rep_energies_and_grad: entered'
    DPRINT '_and_grad: slv_cycle=',slv_cycle
    if(do_grad)then
       ASSERT(present(ma))
       ASSERT(present(ea))
       DPRINT '_and_grad: ma=',ma,' ea=',ea
    endif
    DPRINT '_and_grad: shape(K_solv)=', shape(K_solv)
    DPRINT '_and_grad: shape(R_solv)=', shape(R_solv)
    DPRINT '_and_grad: shape(nm_slv_at)=',shape(nm_slv_at)
    DPRINT '_and_grad: shape(K_solu)=', shape(K_solu)
    DPRINT '_and_grad: shape(R_solu)=', shape(R_solu)
    DPRINT '_and_grad: shape(nm_solute_at)=',shape(nm_solute_at)
    DPRINT '_and_grad: N_unique_atoms=',N_unique_atoms
    m1: do i=1,N_unique_atoms+N_pc
       DPRINT '_and_grad: i=',i
       if(K_solu(i)==0.0_r8_kind .or. R_solu(i)==0.0_r8_kind) then
          cycle m1
       endif
       DPRINT '_and_grad: shape(AA)=',shape(AA)
       D_solv_solu=K_solv(slv_cycle)*K_solu(i)*AA(slv_cycle,i)* &
            64.0_r8_kind*(R_solv(slv_cycle)*R_solu(i))**3
       DPRINT '_and_grad: shape(CC)=',shape(CC)
       B_solv_solu=K_solv(slv_cycle)*K_solu(i)*CC(slv_cycle,i)

       DPRINT '_and_grad: shape(zeta_av)=',shape(zeta_av)
       gamma_solv_solu=zeta_av(slv_cycle,i)/&
                (2.0_r8_kind*sqrt(R_solv(slv_cycle)*R_solu(i)))

       !SSERT(N_points_dr==size(D_R_points))
       if(i <= N_unique_atoms) then
          n_equal=unique_atoms(i)%N_equal_atoms
       else
          n_equal=pointcharge_array(i-N_unique_atoms)%N_equal_charges
       end if
       m2: do j=1,n_equal
          if(i <= N_unique_atoms) then
             xyz(:)=unique_atoms(i)%position(:,j)
          else
             xyz(:)=pointcharge_array(i-N_unique_atoms)%position(:,j)
          end if
          m3: do k=1,N_points_dr
             do l=1,D_R_points(k)%N_equal_points
                R_xyz(:)=D_R_points(k)%position(:,l)-xyz(:)

                R_abs=LEN3(R_xyz)

                disp_p=-D_solv_solu/(3.0_r8_kind*R_abs**6)

                g_r=gamma_solv_solu*R_abs
                rep_p=B_solv_solu*exp(-g_r)*(1.0_r8_kind/g_r+ &
                     2.0_r8_kind/g_r**2+2.0_r8_kind/g_r**3)

                if(.not.do_grad) then
                   E_disp = E_disp + n_slv_at(slv_cycle) * num_dens_dr * disp_p * D_R_points(k) % area * &
                        DOT3 (R_xyz, D_R_points(k) % out_normal(:, l))

                   E_rep = E_rep + n_slv_at(slv_cycle) * num_dens_dr * rep_p * D_R_points(k) % area * &
                        DOT3 (R_xyz, D_R_points(k) % out_normal(:, l))
                else
                   call calc_grads()
                endif
             enddo
          enddo m3
       enddo m2
    enddo m1

  contains
    subroutine calc_grads()

      g_rep1=B_solv_solu*exp(-g_r)*(1.0_r8_kind/g_r**2+4.0_r8_kind/g_r**3 &
           +6.0_r8_kind/g_r**4)
      g_rep2=rep_p+g_rep1

      do l2=1,3
         grad_abs=DOT3(R_xyz, D_R_grad(k)%distance(l2,ma,ea,i,j,:,l))/R_abs

         grad_disp=2.0_r8_kind*D_solv_solu/&
              (R_abs**7)* grad_abs

         grad_rep=-gamma_solv_solu*g_rep2*grad_abs

         grad_p_no=&
              DOT3(R_xyz, D_R_grad(k)%out_normal(l2,ma,ea,:,l))+ &
              DOT3(D_R_grad(k)%distance(l2,ma,ea,i,j,:,l), D_R_points(k)%out_normal(:,l))

         if(present(gradient)) then
            gradient(l2)=gradient(l2)+&
                 n_slv_at(slv_cycle)*num_dens_dr*(&
                 (grad_disp+grad_rep)*&
                 D_R_points(k)%area* &
                 DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+&
                 (disp_p+rep_p)*&
                 D_R_grad(k)%area(l2,ma,ea,l)* &
                 DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+&
                 (disp_p+rep_p)*&
                 D_R_points(k)%area * grad_p_no)
!!$if(slv_cycle==1.and.ma==3.and.ea==1)print*,grad_abs,grad_disp,grad_rep,grad_p_no,gradient(l2),l2,l,k,j,i
         end if
         if(present(hessian)) then
            do l3=1,3
               grad_abs1=DOT3(R_xyz, D_R_grad(k)%distance(l3,ma1,ea1,i,j,:,l))/R_abs

               grad_disp1=2.0_r8_kind*D_solv_solu/&
                    (R_abs**7)* grad_abs1

               grad_rep1=-gamma_solv_solu*g_rep2*grad_abs1

               grad_p_no1=&
                    DOT3(R_xyz, D_R_grad(k)%out_normal(l3,ma1,ea1,:,l))+ &
                    DOT3(D_R_grad(k)%distance(l3,ma1,ea1,i,j,:,l), D_R_points(k)%out_normal(:,l))

               hess_abs=DOT3(D_R_grad(k)%distance(l3,ma1,ea1,i,j,:,l), D_R_grad(k)%distance(l2,ma,ea,i,j,:,l)) &
                    / R_abs - &
                    DOT3(R_xyz,D_R_grad(k)%distance(l2,ma,ea,i,j,:,l))* &
                    DOT3(R_xyz,D_R_grad(k)%distance(l3,ma1,ea1,i,j,:,l))/R_abs**3
               do m=1,3
                  hess_abs=hess_abs+R_xyz(m)*D_R_hess(k)%distance(m,l)%m(l2,ma,ea,l3,ma1,ea1)/R_abs
               end do

               hess_p_no=DOT3(D_R_grad(k)%distance(l2,ma,ea,i,j,:,l), D_R_grad(k)%out_normal(l3,ma1,ea1,:,l))+ &
                    DOT3(D_R_grad(k)%distance(l3,ma1,ea1,i,j,:,l), D_R_grad(k)%out_normal(l2,ma,ea,:,l))
               do m=1,3
                  hess_p_no=hess_p_no+R_xyz(m)*D_R_hess(k)%out_normal(m,l)%m(l2,ma,ea,l3,ma1,ea1)
               end do
               do m=1,3
                  hess_p_no=hess_p_no+ &
                       D_R_hess(k)%distance(m,l)%m(l2,ma,ea,l3,ma1,ea1)*D_R_points(k)%out_normal(m,l)
               end do

               hess_disp=-14.0_r8_kind*(D_solv_solu/R_abs**8)*grad_abs1*grad_abs+ &
                    2.0_r8_kind*(D_solv_solu/R_abs**7)*hess_abs

               hess_rep=-gamma_solv_solu*( &
                    grad_rep1)*grad_abs+ &
                    gamma_solv_solu*( &
                    gamma_solv_solu*g_rep1+ &
                    gamma_solv_solu*B_solv_solu*exp(-g_r)* &
                    (2.0_r8_kind/g_r**3+12.0_r8_kind/g_r**4+24.0_r8_kind/g_r**5) &
                    )*grad_abs1*grad_abs- &
                    gamma_solv_solu*g_rep2*hess_abs

               hessian(l2,l3)=hessian(l2,l3)+n_slv_at(slv_cycle)*num_dens_dr*( &
                    (hess_disp+hess_rep)* &
                    D_R_points(k)%area*DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+ &
                    (grad_disp+grad_rep)* &
                    D_R_grad(k)%area(l3,ma1,ea1,l)*DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+ &
                    (grad_disp+grad_rep)* &
                    D_R_points(k)%area*grad_p_no1+ &
                    (grad_disp1+grad_rep1)* &
                    D_R_grad(k)%area(l2,ma,ea,l)*DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+ &
                    (disp_p+rep_p)* &
                    D_R_hess(k)%area(l)%m(l2,ma,ea,l3,ma1,ea1)*DOT3(R_xyz,D_R_points(k)%out_normal(:,l))+ &
                    (disp_p+rep_p)* &
                    D_R_grad(k)%area(l2,ma,ea,l)*grad_p_no1+ &
                    (grad_disp1+grad_rep1)* &
                    D_R_points(k)%area*grad_p_no+ &
                    (disp_p+rep_p)* &
                    D_R_grad(k)%area(l3,ma1,ea1,l)*grad_p_no+ &
                    (disp_p+rep_p)* &
                    D_R_points(k)%area*hess_p_no)
            end do
         end if
      enddo
    end subroutine calc_grads

  end subroutine disp_rep_energies_and_grad
  !*********************************************************

  !*********************************************************
  subroutine deallocate_D_R_points(do_grad)
    !** End of interface *****************************************

    logical, intent(in) :: do_grad
    integer(kind=i4_kind) :: i,j,k,status

    do i=1,N_points_dr
       deallocate(D_R_points(i)%position,D_R_points(i)%out_normal,stat=status)
       ASSERT(status==0)
       if(do_grad) then
          deallocate(D_R_grad(i)%position,D_R_grad(i)%out_normal,&
               D_R_grad(i)%area,D_R_grad(i)%distance,stat=status)
          ASSERT(status==0)
          if(integralpar_2dervs) then
             do j=1,D_R_points(i)%N_equal_points
                deallocate(D_R_hess(i)%area(j)%m,stat=status)
                ASSERT(status==0)
                do k=1,3
                   deallocate(D_R_hess(i)%out_normal(k,j)%m,stat=status)
                   ASSERT(status==0)
                   deallocate(D_R_hess(i)%distance(k,j)%m,stat=status)
                   ASSERT(status==0)
                end do
             end do
             deallocate(D_R_hess(i)%area,D_R_hess(i)%out_normal, &
                  D_R_hess(i)%distance,stat=status)
             ASSERT(status==0)
          end if
       endif
    enddo

    deallocate(D_R_points,stat=status)
    ASSERT(status==0)

    if(do_grad) then
       deallocate(D_R_grad,stat=status)
       ASSERT(status==0)
       if(integralpar_2dervs) then
          deallocate(D_R_hess,stat=status)
          ASSERT(status==0)
       end if
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
         "deallocate_parameters: deallocation of K_solv is failed")

    deallocate(AA,CC,zeta_av,stat=status)
    if(status/=0) call error_handler("deallocate_parameters, dealloc (1) failed")

   if(do_grad) then
      deallocate(nm_slv_at,n_slv_at,stat=status)
      if ( status /= 0) call error_handler( &
         "disp_rep_write: deallocation of nm_slv_at is failed")
   endif

  end subroutine deallocate_parameters

end module disp_rep_module

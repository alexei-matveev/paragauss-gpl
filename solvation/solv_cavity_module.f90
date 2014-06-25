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
module solv_cavity_module
!
!  This module are used to
!  1) read in input data
!  2) build the molecular cavity and generate
!  the set of cavity surface points using GEPOL_GB
!  algorithm
!  3) calculate of the cavity building energy
!  4) call routines to calculate the dispersion/repulsion energy
!  5) calculation of gradients of that parts of energy
!     (resp the call of appropriate routines in case of dispersion/repulsion)
!
!  The module was prepared by extracting corresponding routines from
!  old solvation_module
!
!== Interrupt of public interface of module =========
!  Author: AS
!  Date: 07/06
!
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
!== Interrupt end of public interface of module =====

!----modules used ------------------
! use FPP_USR_TIMERs:
! define FPP_TIMERS 2
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use filename_module, only: inpfile
  use iounitadmin_module, only: output_unit,write_to_output_units,openget_iounit, &
       returnclose_iounit
  use group_module ,only: symm_transformation_int
  use output_module ,only: output_cavity_data, output_cavity_long
  use atoms_data_module
  use integralpar_module, only: integralpar_2dervs
  use polyhedron_module
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
#endif
#ifdef WITH_EFP
  use efp_module, only: efp, n_efp
#endif
  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public solvation_read, solvation_write, check_read, &
       points_on_cavity_surface, &
       send_receive_geom_grad, dealloc_geom_deriv_part2, &
       dealloc_geom_deriv_part1, dealloc_geom_grad, dealloc_cavity, &
       correction_param

!================================================================
!================================================================
!== Interrupt of public interface of module =========
!----private types
# define MAX_POL_VER 30
  type, private :: poligon
     integer(i4_kind) :: sphere
     integer(kind=i4_kind) :: n_vertises
     real(kind=r8_kind) :: xyz_vertex(MAX_POL_VER,3)
     integer(kind=i4_kind) :: bounds(MAX_POL_VER,2)
     real(kind=r8_kind) :: r_bound(MAX_POL_VER,2)
     real(kind=r8_kind) :: xyz_bound(MAX_POL_VER,6)
     integer(kind=i4_kind) :: n_sphere(MAX_POL_VER,2)
  end type poligon

  type, public :: grad_atomic_center
     real(kind=r8_kind), pointer :: xyz_grad(:,:,:)
     real(kind=r8_kind), pointer :: xyz_hess(:,:,:,:,:,:)
  end type grad_atomic_center

  type, public :: geom_deriv
     type(grad_atomic_center), pointer :: dc(:,:)
     type(grad_atomic_center), pointer :: dR(:)
     type(arrmat1), pointer :: darea(:,:,:)
     type(arrmat1), pointer :: d2area(:,:,:,:,:,:)
     type(arrmat2), pointer :: dcenter(:,:,:)
     type(arrmat2), pointer :: d2center(:,:,:,:,:,:)
  end type geom_deriv

!== Interrupt end of public interface of module =====
!--public types ---
  type, public :: cavity_data
     integer(kind=i4_kind) :: n_equal
     real(kind=r8_kind),pointer :: xyz(:,:)
     real(kind=r8_kind) :: area
     real(kind=r8_kind) :: r_tes
     integer(kind=i4_kind),pointer :: sphere(:)
     logical :: cut_off
  end type cavity_data

!--type collecting information needed for calculating gradients
  type, public:: grad_data
     real(kind=r8_kind) :: dielconst
     integer(kind=i4_kind) :: n_points
     integer(kind=i4_kind),pointer :: n_equal(:)
     real(kind=r8_kind),pointer :: xyz(:,:,:)
     real(kind=r8_kind),pointer :: s(:)
     integer(kind=i4_kind),pointer :: i_symm_sort(:,:)
     type(arrmat1), pointer :: ds_totsyms(:,:)
     type(arrmat1), pointer :: d2s_totsyms(:,:,:,:)
     type(arrmat2), pointer :: dxyz_totsyms(:,:)
     type(arrmat2), pointer :: d2xyz_totsyms(:,:,:,:)
     real(kind=r8_kind),pointer :: Q(:)
     integer(kind=i4_kind),pointer :: sphere(:,:)
#ifdef WITH_EFP
     real(kind=r8_kind),pointer :: Q1(:)
#endif
  end type grad_data

!--public variables ---
  logical, public :: stop_solv
  logical, public :: with_pc=.false.
  logical, public :: fixed_pc=.false.

  integer(kind=i4_kind), public :: n_size
  real(kind=r8_kind), public, allocatable :: Q_n(:), Q_e(:), Q_e_old(:)
#ifdef WITH_EFP
  real(r8_kind), public, allocatable :: Q_mp(:), Q_id(:), Q_id1(:)
#endif
  !if I will ever calculate this:
  logical, public :: cavitation_energy
  logical, public :: disp_rep_energy

  !if I will calculate this now:
  logical , public :: do_cavitation,do_disp_rep
  logical, public :: do_gradients

  ! Caviation energy, not accessed in this module:
  real(r8_kind), public :: energy_cav

  !collected data for calculating the solvent specific
  !gradients:
  type(grad_data),public :: to_calc_grads
  !totsym and cartesian gradient of energies speficic to
  !solvent effects:
  real(kind=r8_kind),public,allocatable,target :: grad_solv_totsym(:)
  real(kind=r8_kind),public,allocatable,target :: grad_solv_totsym_tes(:,:) !!!!!!!!!!AS
  type(arrmat2),public,allocatable :: grad_solv_cart(:)
  type(cavity_data), public, allocatable :: tessarea(:)
  integer(i4_kind),public,allocatable :: center2sphere(:,:)

! End of public interface of module

  logical :: save_cav
  type(symm_transformation_int), pointer :: point_trafos(:)
  integer(kind=i4_kind) :: N_total

  integer(kind=i4_kind),public :: N_spheres
  integer(kind=i4_kind) :: N_atom_spheres,N_bulk_spheres
#ifdef WITH_EFP
  integer(i4_kind) :: N_qm_atom_spheres
#endif
  real(kind=r8_kind), public, allocatable :: r_sphere(:)
  real(kind=r8_kind), public, allocatable :: xyz_sphere(:,:)
  integer(kind=i4_kind), allocatable :: parents(:,:)
  logical, allocatable :: zero_area(:)
  real(kind=r8_kind),public :: R_access
  real(kind=r8_kind) :: m_t_a

  type(geom_deriv),public :: cagr
  integer(kind=i4_kind),public, allocatable :: i_symm_sort(:,:)
  integer(kind=i4_kind), allocatable :: i_unique_at_index(:)
  integer(kind=i4_kind),public :: ua_dim_max = -1

  character(len=5) :: solvation_model
  character(len=4) :: spec_point_group
  character(len=5) :: atomic_radii
  real(kind=r8_kind),public :: solvent_radius
  real(kind=r8_kind) :: scaled_factor
  real(kind=r8_kind), public :: dielectric_constant
  real(kind=r8_kind),public :: abs_temperature
  real(kind=r8_kind),public :: solvent_volume
  integer(kind=i4_kind) :: point_factor
  logical :: get_vdwr
  logical :: no_hydrogen_sphere
  integer(kind=i4_kind) :: external_vdwr
  real(kind=r8_kind) :: max_tes_area
  logical, public :: hydrogen_no_scale
  integer(kind=i4_kind), public :: sol_start_cycle
  integer(kind=i4_kind), public :: n_Q_update
  character(len=13) :: correction_type
  character(len=13) :: cent_type
  logical,public :: cor_el, cor_nuc
  logical :: orig_cent,weight_cent_sin,weight_cent_mass,weight_cent
  real(kind=r8_kind) :: fradio,overlap_angle !old overlap factor for Gepol87
  real (r8_kind) :: rmin_gepol, min_area
  logical :: fix_number_add
  logical,public :: cavitation_all
  !new
  logical :: only_cavities  !cancel calculations just after cavity calculations
  integer :: gepol          !choose gepol algorithm (Gepol87 or Gepol93)
  real(r8_kind) :: overlap_factor !new overlap factor for Gepol93
  logical :: view_cavity
  logical, public :: VTN   !Variable tesserae number algorithm - recomended only for EFP
  logical :: AreaScaling   !Scaling tessera area to avoid strong iteractions between
                           !close located surface point charges - recomended only with EFP
                           !Both VTN and AreaScaling provide just approximated gradients
                           !due to skiped contribution of intersected tesserae
  !Non tessera smoothing partition of cavity surface
  character(len=9) :: smoothing
  integer(i4_kind) :: smooth

  integer(i4_kind) :: atom_rad

  logical :: yes_no = .false.

  namelist /solvation/ &
       solvation_model, &
       smoothing, &
       atomic_radii, &
       VTN, &
       AreaScaling, &
       cavitation_energy, &
       disp_rep_energy, &
       dielectric_constant, &
       abs_temperature, &
       solvent_volume, &
       solvent_radius, &
       scaled_factor, &
       point_factor, &
       max_tes_area, &
       get_vdwr, &
       external_vdwr, &
       no_hydrogen_sphere, &
       hydrogen_no_scale, &
       sol_start_cycle, &
       n_Q_update, &
       correction_type, &
       spec_point_group, &
       orig_cent, &
       weight_cent, &
       cent_type, &
       cavitation_all, &
       fradio, &
       overlap_angle, &
       rmin_gepol, &
       fix_number_add, &
       min_area, &
       only_cavities, &
       gepol, &
       overlap_factor, &
       view_cavity

  integer(kind=i4_kind) :: number_unique_atom
  real(kind=r8_kind) :: vdw_rad
  integer(kind=i4_kind), allocatable :: n_uq_at(:)
  real(kind=r8_kind), allocatable :: vdwr(:)

  namelist /van_der_Waals_radius/ &
       number_unique_atom, &
       vdw_rad


  character(len=5) :: df_solvation_model = "COSMO"
  character(len=4) :: df_spec_point_group = "C1  "
  character(len=5) :: df_atomic_radii = "BONDI" !"UFF"
  logical :: df_cavitation_energy = .true.
  logical :: df_disp_rep_energy = .true.
  real(kind=r8_kind) :: df_solvent_radius = 1.400_r8_kind   !angstrom
  real(kind=r8_kind) :: df_scaled_factor = 1.2_r8_kind
  real(kind=r8_kind) :: df_dielectric_constant = 78.4_r8_kind
  real(kind=r8_kind) :: df_abs_temperature = 298.0_r8_kind    !K
  real(kind=r8_kind) :: df_solvent_volume = 18.07_r8_kind   !cm^3/mol
  integer(kind=i4_kind) :: df_point_factor = 1
  logical :: df_get_vdwr = .false.
  logical :: df_no_hydrogen_sphere = .false.
  integer(kind=i4_kind) :: df_external_vdwr = 0
  real(kind=r8_kind) :: df_max_tes_area = 0.5_r8_kind !angstrom^2
  integer(kind=i4_kind) :: df_sol_start_cycle = 1
  integer(kind=i4_kind) :: df_n_Q_update = 1
  logical :: df_hydrogen_no_scale = .true.
  character(len=13) :: df_correction_type = "Scale_Nuc_El "
  character(len=13) :: df_cent_type =       "default      "
  logical :: df_orig_cent = .false.
  logical :: df_weight_cent = .false.
  real(kind=r8_kind) :: df_fradio = 0.7_r8_kind           !Gepol 87
  real(kind=r8_kind) :: df_overlap_angle = 40.0_r8_kind   !Gepol 87
  real(kind=r8_kind) :: df_rmin_gepol = 0.2_r8_kind
  real (r8_kind) :: df_min_area = 1.0e-7_r8_kind
  logical :: df_fix_number_add = .false.
  logical :: df_cavitation_all = .false.
  !new
  logical :: df_only_cavities = .false.
  integer :: df_gepol=93 !87
  real(r8_kind) :: df_overlap_factor=0.77_r8_kind         !Gepol 93
  logical :: df_view_cavity=.false.
  logical :: df_VTN=.false.
  logical :: df_AreaScaling=.false.
  !smoothing
  character(len=9) :: df_smoothing="NO" ! "TesLess","TesLess_a","FIXPVA"
  real(r8_kind),parameter :: d=0.15_r8_kind
  integer(i4_kind) :: Nd
  real(r8_kind) :: a(7)
  !parameters defining behavior FIXPVA
  real(r8_kind), parameter :: mm11=0.02_r8_kind
  real(r8_kind), parameter :: mm12=0.55_r8_kind
  real(r8_kind), parameter :: nn11=1.0_r8_kind
  real(r8_kind), parameter :: nn12=2.0_r8_kind
!!$  real(r8_kind), parameter :: mm11=0.038_r8_kind  !(0.02A)
!!$  real(r8_kind), parameter :: mm12=0.567_r8_kind  !(0.30A)
!!$  real(r8_kind), parameter :: nn11=1.890_r8_kind  !(1.00A)
!!$  real(r8_kind), parameter :: nn12=2.835_r8_kind  !(1.50A)

  logical, allocatable,public :: skip_short(:)
  integer, allocatable,public :: iuniq(:)

  !to provide graphical representation of tessera
  real(r8_kind), allocatable :: tess_export(:,:,:)
  integer(i4_kind), allocatable :: tess_exp_n(:)
  integer(i4_kind), allocatable :: vert_bounds(:,:,:)

  real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind
  real(kind=r8_kind) , parameter :: ang_au = 0.529177249_r8_kind
  integer(i4_kind), parameter :: N_max_spheres=1000 !maximal number of spheres

  integer(kind=i4_kind) :: counter_n
  external error_handler

contains ! of module

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

    r = SQRT(v(1) * v(1) + v(2) * v(2) + v(3) * v(3))
  end function LEN3
#endif

  !********************************************************************
  subroutine solvation_read
   !reads in the namelist "solvation", sets defaults of elements
   !and checks input values.
   !reads also van der Waals radia if given in input
   !called by read_input
   !** End of interface *****************************************

    use input_module
    use unique_atom_module, only: N_unique_atoms
    use operations_module, only: operations_gradients,operations_geo_opt

    integer(kind=i4_kind) :: unit,status,i

    solvation_model=df_solvation_model
    cavitation_energy=df_cavitation_energy
    disp_rep_energy=df_disp_rep_energy
    dielectric_constant=df_dielectric_constant
    abs_temperature=df_abs_temperature
    solvent_radius=df_solvent_radius
    solvent_volume=df_solvent_volume
    scaled_factor=df_scaled_factor
    point_factor=df_point_factor
    get_vdwr=df_get_vdwr
    no_hydrogen_sphere=df_no_hydrogen_sphere
    hydrogen_no_scale=df_hydrogen_no_scale
    external_vdwr=df_external_vdwr
    max_tes_area=df_max_tes_area
    sol_start_cycle=df_sol_start_cycle
    n_Q_update=df_n_Q_update
    correction_type=df_correction_type
    cent_type=df_cent_type
    spec_point_group=df_spec_point_group
    orig_cent=df_orig_cent
    weight_cent=df_weight_cent
    overlap_angle=df_overlap_angle
    fradio=df_fradio
    rmin_gepol=df_rmin_gepol
    min_area = df_min_area
    fix_number_add = df_fix_number_add
    if(operations_geo_opt) fix_number_add=.true.
    cavitation_all=df_cavitation_all
    only_cavities=df_only_cavities
    gepol=df_gepol
    overlap_factor=df_overlap_factor
    view_cavity=df_view_cavity
    VTN=df_VTN
    AreaScaling=df_AreaScaling
    smoothing=df_smoothing
    smooth=0
    atomic_radii=df_atomic_radii
    atom_rad=0

    unit = input_intermediate_unit()
    call input_read_to_intermediate
    read(unit, nml=solvation, iostat=status)
    if (status .ne. 0) call input_error( &
         "solvation_read: namelist solvation.")
    yes_no = .true.

    stop_solv=only_cavities
    save_cav=view_cavity
    if ( solvation_model/='COSMO' .and. &
         solvation_model/='cosmo')  call input_error( &
         "solvation_read: "// &
         "I do not know anything about choosen solvation model. I deal with COSMO(cosmo) only.")
    if (solvent_radius <= 0.0) call input_error( &
         "solvation_read: You have put wrong solvent_radius. Sorry")
    if (scaled_factor < 1.0) call input_error( &
         "solvation_read: You have put wrong scaled_factor. Sorry")
!!$    if (point_factor > 2 .or. point_factor <= 0) call input_error( &
!!$         "solvation_read: point_factor can not be more then 2 and less 1. Sorry")
    if (dielectric_constant < 1.0 ) call input_error( &
         "solvation_read: You have put wrong dielectric_constant. Sorry")

    if(trim(atomic_radii)=='BONDI' .or. trim(atomic_radii)=='bondi') atom_rad=0
    if(trim(atomic_radii)=='UFF' .or. trim(atomic_radii)=='uff') atom_rad=1

    if(trim(smoothing)=='NO' .or. trim(smoothing)=='no') smooth=0
    if(trim(smoothing)=='TESLESS' .or. trim(smoothing)=='tesless') smooth=1
    if(trim(smoothing)=='TESLESS_A' .or. trim(smoothing)=='tesless_a') smooth=2
    if(trim(smoothing)=='FIXPVA' .or. trim(smoothing)=='fixpva') smooth=3
    if(smooth /= 0) then
       gepol = 93
       VTN=.false.
       AreaScaling=.false.
    end if
    if(smooth == 3) then
       if(point_factor < 2) point_factor=2
       overlap_factor = 0.0_r8_kind
       scaled_factor = 1.125_r8_kind
       no_hydrogen_sphere = .false.
       hydrogen_no_scale = .true.
       atom_rad=0
    endif

#ifdef WITH_EFP
    if(efp) then
       smooth = 0
       VTN=.true.
       AreaScaling=.true.
       no_hydrogen_sphere = .false.
    end if
#endif
    if(VTN) then
       cavitation_energy = .false.
       disp_rep_energy = .false.
    end if
    if(cavitation_energy .and. .not.  disp_rep_energy) cavitation_all=.true.
    if(gepol /= 87 .and. gepol /= 93) gepol = 93
    if(overlap_factor < 0.0_r8_kind .or. overlap_factor >= 1.0_r8_kind) &
         overlap_factor=df_overlap_factor
    if(n_Q_update < 1) n_Q_update=1
    if(VTN) then
       gepol = 93
       overlap_factor = 0.0_r8_kind
    end if
    if(.not.VTN) AreaScaling=.false.
    if(overlap_factor == 0.0_r8_kind) rmin_gepol=0.8_r8_kind

    !type of surface tessera center (representative point) definition
    select case(cent_type)
     case("corners      ")
        !
        ! FIXME: Same as cent_type="triangles"?
        !
        weight_cent=.false.
        weight_cent_mass=.false.
        weight_cent_sin=.false.
        orig_cent=.false.
     case("sin          ")
        weight_cent=.true.
        weight_cent_mass=.false.
        weight_cent_sin=.true.
        orig_cent=.false.
     case("triangles    ")
        !
        ! FIXME: Same as cent_type="corners"?
        !
        weight_cent=.false.
        weight_cent_mass=.false.
        weight_cent_sin=.false.
        orig_cent=.false.
     case("mass         ")
        weight_cent=.true.
        weight_cent_mass=.true.
        weight_cent_sin=.false.
        orig_cent=.false.
     case("default      ")
        !for old kind input options
        weight_cent_mass=weight_cent
        weight_cent_sin=.false.
     case default
       call input_error("solvation_read: cent_type: allowed values are&
                        & corners,mass,triangles,sin,default")
    end select

    select case(correction_type)
     case("None         ")
        cor_el           =.false.
        cor_nuc          =.false.
     case("Scale_Nuc_El ")
        cor_el           =.true.
        cor_nuc          =.true.
     case("Scale_Nuclear")
        cor_el           =.false.
        cor_nuc          =.true.
     case default
       call input_error("solvation_read: correction_type: allowed values are&
                        & None,Scale_Nuc_El,Scale_Nuclear,Charge_Inside ")
    end select
    if(operations_gradients.and.(cor_el.or.cor_nuc)) then
      call write_to_output_units("Warning: In connection with the calculation &
        & of gradients only correction_type=None is allowed. I reset this &
        & value ")
      cor_el =.false.
      cor_nuc=.false.
      correction_type="None         "
    endif

    if(external_vdwr > 0) then
       if(.not. allocated(n_uq_at)) then
          allocate(n_uq_at(external_vdwr),vdwr(external_vdwr),stat=status)
          if ( status /= 0) call error_handler( &
               "solvation_read: allocation of n_uq_at is failed")
       endif
    endif

    do i=1,external_vdwr
       call input_read_to_intermediate()
       read(unit, nml=van_der_Waals_radius, iostat=status)
       n_uq_at(i)=number_unique_atom
       if (number_unique_atom > N_unique_atoms .or. &
            number_unique_atom <= 0 ) call input_error( &
            "solvation_read: You have put wrong number_unique_atom . Sorry")

       vdwr(i)=vdw_rad
       if( n_uq_at(i) == 0 .or. vdwr(i) == 0.0_r8_kind) call input_error( &
            "solvation_read: Please specify something reasonable")
    enddo

  end subroutine solvation_read
  !********************************************************************

  !********************************************************************
  subroutine check_read
   ! sets defaults if namelist "solvation" is not present
   ! called by main_master
   !** End of interface *****************************************
   use operations_module, only: operations_geo_opt

    if (.not.yes_no) then
       solvation_model = df_solvation_model
       cavitation_energy = df_cavitation_energy
       disp_rep_energy = df_disp_rep_energy
       dielectric_constant = df_dielectric_constant
       abs_temperature = df_abs_temperature
       solvent_radius = df_solvent_radius
       solvent_volume=df_solvent_volume
       scaled_factor = df_scaled_factor
       n_Q_update = df_n_Q_update
       point_factor = df_point_factor
       get_vdwr=df_get_vdwr
       no_hydrogen_sphere=df_no_hydrogen_sphere
       hydrogen_no_scale=df_hydrogen_no_scale
       external_vdwr=df_external_vdwr
!!! MF >>>>
       max_tes_area=df_max_tes_area
       sol_start_cycle=df_sol_start_cycle
       correction_type=df_correction_type
       cent_type=df_cent_type
       spec_point_group=df_spec_point_group
       orig_cent=df_orig_cent
       weight_cent=df_weight_cent
       weight_cent_mass=df_weight_cent
       weight_cent_sin=.false.
       fradio=df_fradio
       rmin_gepol=df_rmin_gepol
       min_area = df_min_area
       overlap_angle=df_overlap_angle
       fix_number_add = df_fix_number_add
       if(operations_geo_opt) fix_number_add=.true.
       cavitation_all=df_cavitation_all
       !!! MF <<<<
       only_cavities=df_only_cavities
       stop_solv=only_cavities
       gepol=df_gepol
       overlap_factor=df_overlap_factor
       view_cavity=df_view_cavity
       save_cav=view_cavity
       VTN=df_VTN
       AreaScaling=df_AreaScaling
       smoothing=df_smoothing
       smooth=0
       atomic_radii=df_atomic_radii
       atom_rad=0
#ifdef WITH_EFP
      if(efp) then
         smooth = 0
         cavitation_energy = .false.
         disp_rep_energy = .false.
         VTN=.true.
         AreaScaling=.true.
         no_hydrogen_sphere = .false.
      end if
      if(VTN) then
         gepol = 93
         overlap_factor = 0.0_r8_kind
      end if
#endif
    endif

  end subroutine check_read
  !********************************************************************

  !********************************************************************
  subroutine solvation_write(unit)
    !
    ! Writes  the used  parameters of  namelist "solvation"  called by
    ! write_input.
    !
    use echo_input_module, only: start, real, flag, intg, word, stop, &
         echo_level_full, real_format1, real_format2, word_format
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(kind=i4_kind), intent(in) :: unit
    ! *** end of interface ***

    integer(kind=i4_kind) :: i

    real_format1 = '("    ",a," = ",f7.3:" # ",a)'
    real_format2 = '("    ",a," = ",es12.3:" # ",a)'
    word_format = '("    ",a," = ",a15 :" # ",a)'

    call start("SOLVATION","SOLVATION_WRITE",unit,operations_echo_input_level)
    call word("SOLVATION_MODEL    ",solvation_model    ,df_solvation_model    )
!    VTN and AreaScaling options are removed from input.out not to confuse users
!    as they are used for EFP only (default values) and are not recomended for standard
!    calculations
!    call flag("VTN                ",VTN                ,df_VTN                )
!    call flag("AREASCALING        ",AreaScaling        ,df_AreaScaling        )
    call flag("CAVITATION_ENERGY  ",cavitation_energy  ,df_cavitation_energy  )
    call flag("DISP_REP_ENERGY    ",disp_rep_energy    ,df_disp_rep_energy    )
    call flag("ONLY_CAVITIES      ",only_cavities      ,df_only_cavities      )
    call flag("VIEW_CAVITY        ",view_cavity        ,df_view_cavity        )
    call real("DIELECTRIC_CONSTANT",dielectric_constant,df_dielectric_constant,1)
    call real("ABS_TEMPERATURE    ",abs_temperature    ,df_abs_temperature    ,1)
    call word("ATOMIC_RADII       ",atomic_radii       ,df_atomic_radii       )
    call word("SMOOTHING          ",smoothing          ,df_smoothing          )
    call intg("GEPOL              ",gepol              ,df_gepol              )
    call real("SOLVENT_VOLUME     ",solvent_volume     ,df_solvent_volume     ,1)
    call real("SOLVENT_RADIUS     ",solvent_radius     ,df_solvent_radius     ,1)
    call real("SCALED_FACTOR      ",scaled_factor      ,df_scaled_factor      ,1)
    call intg("POINT_FACTOR       ",point_factor       ,df_point_factor       )
    call real("MAX_TES_AREA       ",max_tes_area       ,df_max_tes_area       ,2)
!!$    call flag("GET_VDWR           ",get_vdwr           ,df_get_vdwr           )
    call intg("EXTERNAL_VDWR      ",external_vdwr      ,df_external_vdwr      )
    call flag("NO_HYDROGEN_SPHERE ",no_hydrogen_sphere ,df_no_hydrogen_sphere )
    call flag("HYDROGEN_NO_SCALE  ",hydrogen_no_scale  ,df_hydrogen_no_scale  )
    call intg("SOL_START_CYCLE    ",sol_start_cycle    ,df_sol_start_cycle    )
    call intg("N_Q_UPDATE         ",n_Q_update         ,df_n_Q_update         )
    call word("CORRECTION_TYPE    ",correction_type    ,df_correction_type    )
       write (unit,'(A100)') &
            "   # Allowed values are: None Scale_Nuc_El Scale_Nuclear  "
    call word("CENT_TYPE          ",cent_type          ,df_cent_type          )
       write (unit,'(A100)') &
            "   # Allowed values are: corners, mass, sin, triangles, default "
    call flag("WEIGHT_CENT        ",weight_cent        ,df_weight_cent        )
       write (unit,'(A100)') &
            "   # (only for old fashion input (without cent_type) "
    call real("FRADIO             ",fradio             ,df_fradio,1)
    call real("OVERLAP_ANGLE      ",overlap_angle      ,df_overlap_angle,1)
    call real("OVERLAP_FACTOR     ",overlap_factor     ,df_overlap_factor,1)
    call real("RMIN_GEPOL         ",rmin_gepol         ,df_rmin_gepol   ,2)
    call word("SPEC_POINT_GROUP   ",spec_point_group   ,df_spec_point_group   )
    call flag("FIX_NUMBER_ADD     ",fix_number_add     ,df_fix_number_add)
    call flag("CAVITATION_ALL     ",cavitation_all     ,df_cavitation_all)
    call real("MIN_AREA           ",min_area           ,df_min_area   ,2)
    call stop()

    do i=1,external_vdwr
       call start("VAN_DER_WAALS_RADIUS","VAN_DER_WAALS_RADIUS_WRITE",unit,operations_echo_input_level)
       call intg("NUMBER_UNIQUE_ATOM",n_uq_at(i)        ,0                     )
       call real("VDW_RAD           ",vdwr(i)           ,0.0_r8_kind           )
       call stop()
    enddo

  end subroutine solvation_write
  !**************************************************************

  !**************************************************************
  subroutine correction_param()
    use solv_charge_mixing_module, only: mix_charges

    if(mix_charges()) then
       cor_el           =.false.
       cor_nuc          =.false.
    end if

  end subroutine correction_param
  !**************************************************************

  !**************************************************************
  subroutine points_on_cavity_surface
    !
    ! Calls  routines  for producing  the  appropriate cavity  surface
    ! grid:
    !
    ! 1) calc_cavity() or calc_caity_1(): find set of spheres building
    !    the cavity
    !
    ! 2) denerate_cube(), generate_dodecahedron(), or
    !    generate_doublepyramide(): produce points on a spherical
    !    surface with local symmetry, starting from a appropriate
    !    polihedron.
    !
    ! 3) tesselation(): calculates areas and centers of cutted (by
    !    other spheres) and non cutted surface pieces (tesserae)
    !
    ! 4) symm_sorted_centers(): find and store symmetry equivalent
    !    surface tesserae
    !
    ! Called by main_master(), main_gradient() and (intern)
    ! disp_rep_wrap().
    !
    use unique_atom_module
    use pointcharge_module
#ifdef WITH_EFP
    use point_dqo_module
    use induced_dipoles_module
#endif
    use symmetry_data_module, only : symmetry_data_point_group
    use group_module, only : ylm_trafos,sub_group,group_coset, &
         symm_transformation_int,group_num_el,group_coset_decomp
    use symm_module, only : symm_adapt_centers
    use help_cavity_module
    use cavity_image_module
    implicit none
    !** End of interface *****************************************

    character(len=4) :: name_point_group
    real(kind=r8_kind), allocatable :: xyz_tes_c(:,:),area_tes(:),r_tes(:)           !! for tesselation
    integer(kind=i4_kind), allocatable :: sphere(:)                                  !! for tesselation
    logical, allocatable :: cuttt(:)                                                 !! for tesselation
    type(poligon),allocatable :: data_tes(:)                                         !! for tesselation
!!! MF >>>>
    integer (kind=i4_kind), allocatable :: cut_rad_sort(:,:)
!!! MF <<<<

    integer(kind=i4_kind) :: i,j !,k
    integer(kind=i4_kind) :: i_poly
    integer(kind=i4_kind) :: ret_stat,status

    external error_handler
    if(.not.do_cavitation .and. .not.do_disp_rep) then
       if(output_cavity_data .and. .not.do_gradients) then
          write (output_unit,*) '*******************************************************************'
          write (output_unit,*) '           SOLVATION EFFECT - ELECTROSTATIC CONTRIBUTION           '
          write (output_unit,*) '                       C O S M O - model                           '
          write (output_unit,*) '==================================================================='
          write (output_unit,*) 'Definition of cavity and coordinats of point charges on its surface'
          if(smooth == 0 .and. gepol == 93) &
          write (output_unit,*) '                        GEPOL-93 scheme                            '
          if(smooth == 1) &
          write (output_unit,*) '                    TesLess smoothing scheme                       '
          if(smooth == 3) &
          write (output_unit,*) '                     FIXPVA smoothing scheme                       '
          if(atom_rad == 0) &
          write (output_unit,*) '                    Bondi + UFF atomic radii                       '
          if(atom_rad == 1) &
          write (output_unit,*) '                        UFF atomic radii                           '
          write (output_unit,*) '==================================================================='
       endif
!!$       if(get_vdwr.and. .not.do_gradients ) then
!!$          write (output_unit,*) '       Default values of van der Waals radii     '
!!$          write (output_unit,*) '-------------------------------------------------'
!!$          write (output_unit,*) 'Element      van der Waals radius(A)'
!!$          do i=1,98
!!$          if(vdW_radius(i)/=0.0_r8_kind) then
!!$                write (output_unit,'(3x,a2,17x,f5.2,1x,a5)')atom_name(i),vdW_radius(i),'Bondi'
!!$             else
!!$                write (output_unit,'(3x,a2,17x,f5.2,1x,a5)')atom_name(i),R_def_rap(i)/2.0_r8_kind,'Rappe'
!!$             endif
!!$          enddo
!!$          write (output_unit,*) '-------------------------------------------------'
!!$       endif
    endif

    if(do_cavitation .and. output_cavity_data .and. .not.do_gradients) then
       write (output_unit,*) '*******************************************************************'
       write (output_unit,*) '          SOLVATION EFFECT - CAVITATION CONTRIBUTION'
       write (output_unit,*) '                    Definition of cavity'
       write (output_unit,*) '==================================================================='
    endif

    if(do_disp_rep .and. output_cavity_data .and. .not.do_gradients) then
       write (output_unit,*) '*******************************************************************'
       write (output_unit,*) '       SOLVATION EFFECT - DISPERSION-REPULSION CONTRIBUTION'
       write (output_unit,*) '                    Definition of cavity'
       write (output_unit,*) '==================================================================='
    endif

    ua_dim_max = maxval(unique_atoms(:)%N_equal_atoms)
#ifndef WITH_EFP
    if(with_pc) ua_dim_max = max(maxval(pointcharge_array(:)%N_equal_charges),ua_dim_max)
#else
    if(with_pc .or. efp) ua_dim_max = max(maxval(pointcharge_array(:)%N_equal_charges),ua_dim_max)
#endif
    DPRINT 'points_on_cavity_surface: ua_dim_max=',ua_dim_max

    if(do_disp_rep .or. do_cavitation) then
       if(do_cavitation) R_access=0.0_r8_kind
       call calc_cavity_1()
!!$       print*,'calc_cavity_1 N_spheres', N_spheres
       if(do_gradients) call calc_cavity_gradient(.false.)
    else
       if(gepol == 87) then
          call calc_cavity()
       elseif(gepol == 93) then
          if(symmetry_data_point_group() =='C1  ') NDIV=4
          name_point_group='C1  '
          call generate_dodecahedron(1,point_factor,name_point_group, &
               do_cavitation,do_gradients)
          call calc_cavity_93()
          deallocate(surf_elem%xyz_centers,surf_elem%xyz,surf_elem%index)
!!$          print*,'calc_cavity N_spheres', N_spheres
       end if
       if(save_cav .and. .not.do_gradients) call save_cavity_image( &
            gepol,N_spheres,xyz_sphere,r_sphere,zero_area)
       if(do_gradients) call calc_cavity_gradient(.true.)
    endif

!   name_point_group=symmetry_data_point_group()
    if(spec_point_group==df_spec_point_group) then
       name_point_group =symmetry_data_point_group()
    else
        name_point_group=spec_point_group
    endif

!!! MF >>>>
    m_t_a=max_tes_area/(ang_au**2)
    call find_max_dim_surf_elem()
!!! MF <<<<

    i_poly=0
    if ( N_atom_spheres >= MAX_ATOMS_QMMM .and. &
         name_point_group=='C1  ') then
       if(output_cavity_data .and. .not.do_gradients) &
            write (output_unit,*) 'The OCTAHEDRON is inscribed into each sphere'
       call generate_octahedron(N_atom_spheres,do_cavitation,do_gradients)
       i_poly=1
    elseif ( (name_point_group=='C1  ')  .or. &
         (name_point_group=='C2  ')  .or. &
         (name_point_group=='C3  ')  .or. &
         (name_point_group=='C5  ')  .or. &
         (name_point_group=='CS  ')  .or. &
         (name_point_group=='Ci  ')  .or. &
         (name_point_group=='C2V ')  .or. &
         (name_point_group=='C2H ')  .or. &
         (name_point_group=='C3V ')  .or. &
         (name_point_group=='C5V ')  .or. &
         (name_point_group=='S6  ')  .or. &
         (name_point_group=='S10 ')  .or. &
         (name_point_group=='D2  ')  .or. &
         (name_point_group=='D2H ')  .or. &
         (name_point_group=='D3  ')  .or. &
         (name_point_group=='D3D ')  .or. &
         (name_point_group=='D5  ')  .or. &
         (name_point_group=='D5D ')  .or. &
         (name_point_group=='I   ')  .or. &
         (name_point_group=='IH  ')  ) then
       if(output_cavity_data .and. .not.do_gradients) &
            write (output_unit,*) 'The DODECAHEDRON is inscribed into each sphere'
       call generate_dodecahedron(0,point_factor,name_point_group, &
               do_cavitation,do_gradients)
    elseif ( (name_point_group=='C4  ') .or. &
         (name_point_group=='C4V ')  .or. &
         (name_point_group=='C4H ')  .or. &
         (name_point_group=='D4  ')  .or. &
         (name_point_group=='D4H ')  .or. &
         (name_point_group=='S4  ')  .or. &
         (name_point_group=='D2D ')  .or. &
         (name_point_group=='O   ')  .or. &
         (name_point_group=='OH  ')  .or. &
         (name_point_group=='T   ')  .or. &
         (name_point_group=='TH  ')  .or. &
         (name_point_group=='TD  ') ) then
       if(output_cavity_data .and. .not.do_gradients) &
            write (output_unit,*) 'The CUBE is inscribed into each sphere'
       call generate_cube(point_factor,do_cavitation,do_gradients)
    else
!!! MF >>>>
       n_rotations=0
       if( (name_point_group=='C3H ')  .or. &
            (name_point_group=='D3H ') ) then
          n_rotations=3
       elseif( (name_point_group=='C5H ')  .or. &
            (name_point_group=='D5H ') ) then
          n_rotations=5
       elseif( (name_point_group=='C6  ')  .or. &
            (name_point_group=='C6V ')  .or. &
            (name_point_group=='C6H ')  .or. &
            (name_point_group=='D6  ')  .or. &
            (name_point_group=='D6H ') ) then
          n_rotations=6
       elseif( (name_point_group=='C7  ')  .or. &
            (name_point_group=='C7V ')  .or. &
            (name_point_group=='C7H ')  .or. &
            (name_point_group=='D7  ')  .or. &
            (name_point_group=='D7H ') ) then
          n_rotations=7
       elseif( (name_point_group=='C8  ')  .or. &
            (name_point_group=='C8V ')  .or. &
            (name_point_group=='C8H ')  .or. &
            (name_point_group=='D4D ')  .or. &
            (name_point_group=='S8  ')  .or. &
            (name_point_group=='D8  ')  .or. &
            (name_point_group=='D8H ') ) then
          n_rotations=8
       elseif( (name_point_group=='C9  ')  .or. &
            (name_point_group=='C9V ')  .or. &
            (name_point_group=='C9H ')  .or. &
            (name_point_group=='D9  ')  .or. &
            (name_point_group=='D9H ') ) then
          n_rotations=9
       elseif( (name_point_group=='C10 ')  .or. &
            (name_point_group=='C10V')  .or. &
            (name_point_group=='C10H')  .or. &
            (name_point_group=='D10 ')  .or. &
            (name_point_group=='D10H') ) then
          n_rotations=10
       elseif( (name_point_group=='S12 ')  .or. &
            (name_point_group=='D6D ')) then
          n_rotations=12
       elseif( (name_point_group=='S14 ')  .or. &
            (name_point_group=='D7D ')) then
          n_rotations=14
       elseif( (name_point_group=='S16 ')  .or. &
            (name_point_group=='D8D ')) then
          n_rotations=16
       elseif( (name_point_group=='S18 ')  .or. &
            (name_point_group=='D9D ')) then
          n_rotations=18
       elseif( (name_point_group=='S20 ')  .or. &
            (name_point_group=='D10D')) then
          n_rotations=20
       endif
       if(n_rotations==0) then
          call error_handler( &
               "points_on_cavity_surface: S12 ... S20 and  &
               & D6D ... D10D are yet missing &
               & in the solvent model.Sorry")
       endif

       if(output_cavity_data .and. .not.do_gradients) &
            write (output_unit,*) 'The DOUBLEPIRAMIDE is inscribed into each sphere'
       call generate_doublepyramide(point_factor,do_cavitation,do_gradients)
!!! MF <<<<
    endif

    call init_cut_rad_sort
    counter_n=0
!!$    print*,'tesselation'
    call tesselation

    if (.not.do_cavitation .and. .not.do_disp_rep) then
       if(save_cav .and. .not.do_gradients .and. smooth == 0) then
          allocate(tess_export(N_total,MAX_POL_VER,3), tess_exp_n(N_total), &
               vert_bounds(N_total,MAX_POL_VER,2),stat=status)
          ASSERT(status==0)

          do i=1,N_total
             tess_exp_n(i)=data_tes(i)%n_vertises
             do j=1,data_tes(i)%n_vertises
                tess_export(i,j,:)=data_tes(i)%xyz_vertex(j,:)
                vert_bounds(i,j,:)=data_tes(i)%bounds(j,:)
             end do
          end do

          call save_tess_image(N_total,tess_exp_n,tess_export,vert_bounds)

          deallocate(tess_export,tess_exp_n,vert_bounds,stat=status)
          ASSERT(status==0)
       end if
    end if

!!$    print *,'          call symm_sorted_centers(ret_stat,verbose=.false.)'
    call symm_sorted_centers(ret_stat,verbose=.false.)
!!$    if(ret_stat/=0) then
!!$       print *,'retrying  call symm_sorted_centers(ret_stat,verbose=.true.)'
!!$       call symm_sorted_centers(ret_stat,verbose=.true.)
!!$       if(ret_stat/=0)then
!!$          stop 'symm_sorted_centers failed'
!!$       endif
!!$    endif
    DPRINT 'end of points of cavity'
    if(.not.do_cavitation .and. .not.do_disp_rep) then
       DPRINT 'call matrix'
!       if(.not. do_gradients) call matrix_generation
       DPRINT 'end matrix'
       if(output_cavity_data .and. .not.do_gradients ) then
          write (output_unit,*) '        THE END OF THE SOLUTE CAVITY POINTS GENERATION'
          write (output_unit,*) '*****************************************************************'
       endif
    endif
    if(do_cavitation .or. do_disp_rep .and. output_cavity_data .and. .not.do_gradients ) then
       write (output_unit,*) '*****************************************************************'
    endif

    DPRINT 'end of points of cavity'
  contains ! of points_on_cavity_surface


    !------------------------------------------------------------
    ! Public interface of module
    subroutine calc_cavity_1
      !private subroutine
      !solvent accessible surface (only for dispersion-repulsion and cavitation terms)
      ! End of public interface of module

      real(kind=r8_kind) :: v_d_w_r
      integer(kind=i4_kind) :: vdW_index
      logical :: ext
      integer(kind=i4_kind) :: i,j,k,status

      N_spheres = 0
      do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1_i4_kind) cycle
         N_spheres = N_spheres + unique_atoms(i)%N_equal_atoms
      enddo
      if(with_pc) then
         do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle
            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 100
            end do
            cycle
100         N_spheres = N_spheres + pointcharge_array(i)%N_equal_charges
         end do
      end if
      N_atom_spheres=N_spheres

      allocate(r_sphere(N_spheres),iuniq(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface(1): allocation of r_sphere is failed")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface(1): allocation of xyz_sphere is failed")

      if(do_gradients) then
         allocate(i_unique_at_index(N_spheres),stat=status)
         if(status/=0) call error_handler("calc cavity_1: allocation &
              & of i_unique_at_index is failed")
         i_unique_at_index(:)=0
      endif

      k=0
      n_uni: do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1) cycle n_uni
         ext=.false.
         if(external_vdwr /= 0) then
            do j=1,external_vdwr
               if(i ==  n_uq_at(j)) then
                  vdW_index=j
                  ext=.true.
               endif
            enddo
         endif
         if(external_vdwr /= 0 .and. ext) then
            v_d_w_r=vdwr(vdW_index)
         else
            vdW_index=int(unique_atoms(i)%Z)
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            if(v_d_w_r == 0.0_r8_kind ) then
               write (output_unit,*) i,' -unique_atom has no default value of van der Waals radius'
               call error_handler( &
                    "points_on_cavity_surface: one of atoms in the solute has no &
                    & default value of van der Waals radius (See output file). If &
                    & You want to know no default van der Waals radii please specify &
                    & option GET_VDWR")
            endif
         endif
         do j=1,unique_atoms(i)%N_equal_atoms
            k=k+1
            r_sphere(k)=v_d_w_r/ang_au+R_access
            iuniq(k)=i
            xyz_sphere(k,:)=unique_atoms(i)%position(:,j)
            ASSERT(ua_dim_max.ne.-1)
            if(do_gradients)then
               i_unique_at_index(k)= (i - 1) * ua_dim_max + (j - 1)
               ASSERT(j.le.ua_dim_max)
            endif
         enddo
      enddo n_uni

      if(with_pc) then
         n_uni_pc: do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle n_uni_pc

            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 200
            end do
            cycle n_uni_pc
200         vdW_index=j
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            do j=1,pointcharge_array(i)%N_equal_charges
               k=k+1
               r_sphere(k)=v_d_w_r/ang_au+R_access
               iuniq(k)=i
               xyz_sphere(k,:)=pointcharge_array(i)%position(:,j)
               ASSERT(ua_dim_max.ne.-1)
               if(do_gradients)then
                  i_unique_at_index(k)= (i+N_unique_atoms - 1) * ua_dim_max + (j - 1)
                  ASSERT(j.le.ua_dim_max)
               endif
            enddo
         end do n_uni_pc
      end if

      N_bulk_spheres=N_atom_spheres

      if(output_cavity_data .and. .not.do_gradients ) then
         write (output_unit,'(a23,i4,a19)') 'The cavity consists of ',N_spheres,' overlaping spheres'
         write (output_unit,*) 'number       radius(a.u.)                 coordinates(a.u)'
         do i=1,N_spheres
            write (output_unit,'(1x,i4,6x,f11.8,6x,3(f13.9,1x))') i,r_sphere(i),xyz_sphere(i,:)
         enddo
         write (output_unit,*) '---------------------------------------------------'
      endif

    end subroutine calc_cavity_1
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine calc_cavity
      ! the solvent excluding surface
      !** End of interface *****************************************

      real(kind=r8_kind) :: r_mini, r_new
      real(kind=r8_kind), allocatable :: r_sphere_buf(:)
      real(kind=r8_kind), allocatable :: xyz_sphere_buf(:,:)
      integer(kind=i4_kind), allocatable :: parents_buf(:,:)
      integer(kind=i4_kind), pointer :: fix_par(:,:)
      integer(kind=i4_kind) :: n_pair, N_spheres_old, n_par_fix,n_max_par

      real(kind=r8_kind) :: distance12, distance13, distance23, d_13, d_23, &
           dist_overlap, dist_h, dist_hiding !, sp
      real(kind=r8_kind) :: delta_xyz12(3), delta_xyz13(3), delta_xyz23(3)

      real(kind=r8_kind) :: r_solv,v_d_w_r,s_f,help,rmin,rmax
      real(kind=r8_kind) :: overlap_angle_rad,coso,sino

      integer(kind=i4_kind) :: vdW_index, one_more_sphere
      logical :: ext,use_stored,parents_file_exists
      integer(kind=i4_kind) :: i,j,k,l,status !,ii

      !converting parameters to radian and atomic units
      overlap_angle_rad= overlap_angle*pi/180.0_r8_kind
      sino=sin(overlap_angle_rad)
      coso=cos(overlap_angle_rad)
      r_solv=solvent_radius/ang_au !!!!!!!
      r_mini=rmin_gepol/ang_au     !!!!!!
      if(do_cavitation) then
         s_f=1.0_r8_kind
      else
         s_f=scaled_factor
      endif

      parents_file_exists=.false.
      if(fix_number_add) call read_fixed_parents(parents_file_exists,n_max_par,fix_par)
      use_stored=(fix_number_add .and. parents_file_exists)
      if(use_stored) write(output_unit,*) "THE CAVITY CONFIGURATIONS IS TAKEN FROM FIX_PAR FILE"

      N_spheres = 0
      do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1_i4_kind) cycle
         N_spheres = N_spheres + unique_atoms(i)%N_equal_atoms
      enddo
      if(with_pc) then
         do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle
            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 100
            end do
            cycle
100         N_spheres = N_spheres + pointcharge_array(i)%N_equal_charges
         end do
      end if
      N_atom_spheres=N_spheres

      if(do_gradients) then
         allocate(i_unique_at_index(N_spheres),stat=status)
         if(status/=0) call error_handler("calc cavity: allocation &
              & of i_unique_at_index is failed")
         i_unique_at_index(:)=0
      endif

      allocate(r_sphere(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of r_sphere is failed")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of xyz_sphere is failed")
      allocate(parents(N_spheres,2), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of parents is failed")

      !definition of initial spheres
      k=0
      n_uni: do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1) cycle n_uni
         ext=.false.
         if(external_vdwr /= 0) then
            do j=1,external_vdwr
               if(i ==  n_uq_at(j)) then
                  vdW_index=j
                  ext=.true.
               endif
            enddo
         endif
         if(external_vdwr /= 0 .and. ext) then
            v_d_w_r=vdwr(vdW_index)
         else
            vdW_index=int(unique_atoms(i)%Z)
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            if(v_d_w_r == 0.0_r8_kind ) then
               write (output_unit,*) i,' -unique_atom has no default value of van der Waals radius'
               call error_handler( &
                    "points_on_cavity_surface: one of atoms in the solute has no &
                    & default value of van der Waals radius (See output file). If &
                    & You want to know no default van der Waals radii please specify &
                    & option GET_VDWR")
            endif
         endif
         do j=1,unique_atoms(i)%N_equal_atoms
            k=k+1
            if((.not.hydrogen_no_scale).or.(unique_atoms(i)%Z/=1.).or.&
                 (do_cavitation)) then
               r_sphere(k)=v_d_w_r*s_f/ang_au  !!!!!
            else
               r_sphere(k)=v_d_w_r/ang_au      !!!!!!!!
            endif
            xyz_sphere(k,:)=unique_atoms(i)%position(:,j) !*ang_au !!!!!!!
            parents(k,:)=0
            ASSERT(ua_dim_max.ne.-1)
            if(do_gradients)then
               i_unique_at_index(k)= (i - 1) * ua_dim_max + (j - 1)
               ASSERT(j.le.ua_dim_max)
            endif
         enddo
      enddo n_uni

      if(with_pc) then
         n_uni_pc: do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle n_uni_pc

            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 200
            end do
            cycle n_uni_pc
200         vdW_index=j
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            do j=1,pointcharge_array(i)%N_equal_charges
               k=k+1
               if((.not.hydrogen_no_scale).or.(trim(pointcharge_array(i)%name) /= "H").or.&
                    (do_cavitation)) then
                  r_sphere(k)=v_d_w_r*s_f/ang_au  !!!!!
               else
                  r_sphere(k)=v_d_w_r/ang_au
               end if
               xyz_sphere(k,:)=pointcharge_array(i)%position(:,j)
               ASSERT(ua_dim_max.ne.-1)
               parents(k,:)=0
               if(do_gradients)then
                  i_unique_at_index(k)= (i+N_unique_atoms - 1) * ua_dim_max + (j - 1)
                  ASSERT(j.le.ua_dim_max)
               endif
            enddo
         end do n_uni_pc
      end if

      N_bulk_spheres=N_atom_spheres

      !start main cycle
!       print*,'start main cycle use stored',use_stored
      N_spheres_old=0
      n_par_fix=1

      main_cycle_0: do
         n_pair = (N_spheres*(N_spheres-1_i4_kind))/2_i4_kind
         allocate(r_sphere_buf(N_spheres+n_pair),stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of r_sphere_buf is failed")
         allocate(xyz_sphere_buf(N_spheres+n_pair,3), stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of xyz_sphere is failed")
         allocate(parents_buf(N_spheres+n_pair,2), stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of parants_buf is failed")

         r_sphere_buf(1:N_spheres)=r_sphere(:)
         xyz_sphere_buf(1:N_spheres,:)=xyz_sphere(:,:)
         parents_buf(1:N_spheres,:)=parents(:,:)

         one_more_sphere=0
         main_cycle_1: do i=1,N_spheres-1
            l=i+1
            if(l < N_spheres_old+1) l= N_spheres_old+1
            main_cycle_2: do j=l,N_spheres

               if(use_stored ) then
!                   print*,'fix_par', fix_par(n_par_fix,1),fix_par(n_par_fix,2)
                    if(fix_par(n_par_fix,1)==i .and. fix_par(n_par_fix,2)==j) then
                        n_par_fix=n_par_fix+1
                    else
                        cycle main_cycle_2
                    endif
               endif

               delta_xyz12(:) = xyz_sphere(j,:) - xyz_sphere(i,:)
               distance12 = LEN3(delta_xyz12)

               if(.not. use_stored) then
                 !checking distance between two spheres
                 if (distance12 >= (r_sphere(i)+r_sphere(j)+2.0_r8_kind*r_solv)) cycle  main_cycle_2

                 !checking overlaping of two spheres
                 if (r_sphere(i) <= r_sphere(j)) then
                    rmin=r_sphere(i); rmax=r_sphere(j)
                 else
                    rmin=r_sphere(j); rmax=r_sphere(i)
                 endif
                 help=rmax**2-(rmin*sino)**2
                 if(help<0.0_r8_kind) help=0.0_r8_kind
                 dist_overlap=rmin*coso+sqrt(help)
!                 help=rmax**2-(rmin*coso)**2
!                 if(help<0.0_r8_kind) help=0.0_r8_kind
!                 dist_overlap=rmin*sino+sqrt(help)

                 if (distance12 <= dist_overlap) cycle  main_cycle_2
!                 help=((rmax+r_solv)**2+rmax**2-(r_solv+r_mini)**2)/rmax      !??????????
!                 if (distance12 >= help) cycle  main_cycle_2                  !??????????

                 !are there any other spheres in between two spheres
                 check_hiding: do k=1,N_spheres
                   if (k == i .or. k == j) cycle

                   delta_xyz13(:) = xyz_sphere(k,:) - xyz_sphere(i,:)
                   distance13 = LEN3(delta_xyz13)
                   delta_xyz23(:) = xyz_sphere(k,:) - xyz_sphere(j,:)
                   distance23 = LEN3(delta_xyz23)

!                   d_13=abs(DOT3(delta_xyz13,delta_xyz12)/distance12)
                   d_13=DOT3(delta_xyz13,delta_xyz12)/distance12
                   if(abs(d_13) <= 1.0e-8_r8_kind) d_13=0.0_r8_kind
!                   d_23=abs(DOT3(delta_xyz23,delta_xyz12)/distance12)
                   d_23=DOT3(delta_xyz23,-delta_xyz12)/distance12
                   if(abs(d_23) <= 1.0e-8_r8_kind) d_23=0.0_r8_kind
                   if(d_13 <= 0.0_r8_kind .or. d_23 <= 0.0_r8_kind) cycle

                   help=1.0_r8_kind-(DOT3(delta_xyz13,delta_xyz12)/(distance13*distance12))**2
                   if(help<0.0_r8_kind) help=0.0_r8_kind
                   dist_h=distance13*sqrt(help)

!                   sp = (distance12+distance13+distance23)/2.0_r8_kind      !old fashion
!                   help=4.0_r8_kind*(sp*(sp-distance12)*(sp-distance13)*(sp-distance23))/distance12**2
!                   dist_h=sqrt(help)

                   dist_hiding = r_sphere(k)*fradio
!                   if( dist_h <= dist_hiding .and. abs(distance12-d_13-d_23)<=1.0e-8_r8_kind ) then
                   if( dist_h <= dist_hiding) then
                      cycle main_cycle_2
                   end if
                 enddo check_hiding

               endif !use_stored

               ! creation of new sphere

               if ( distance12 <= (r_sphere(i)+r_sphere(j))) then
                  !case A of GEPOL algorithm
                  help= (r_sphere(j)+r_solv)**2+0.5_r8_kind*(distance12+r_sphere(j)- &
                       r_sphere(i))*(0.5_r8_kind*(distance12+r_sphere(j)-r_sphere(i))- &
                       ((r_sphere(j)+r_solv)**2+distance12**2-(r_sphere(i)+r_solv)**2)/distance12)
                  if(help<0.0_r8_kind) help=0.0_r8_kind
                  r_new= sqrt(help)-r_solv
                  if ( r_new < r_mini .and. .not. use_stored ) cycle main_cycle_2
!                 print*,'case A of GEPOL algorithm',distance12,r_sphere(i),r_sphere(j),use_stored
                  one_more_sphere=one_more_sphere+1
                  r_sphere_buf(N_spheres+one_more_sphere)=r_new
                  xyz_sphere_buf(N_spheres+one_more_sphere,:)=&
                       0.5_r8_kind*(xyz_sphere(j,:)+xyz_sphere(i,:))- &
                       (r_sphere(j)-r_sphere(i))*(xyz_sphere(j,:)- &
                       xyz_sphere(i,:))/(2.0_r8_kind*distance12)
               else
                  help= (r_sphere(j)+r_solv)**2+0.5_r8_kind*(distance12+r_sphere(j)- &
                       r_sphere(i))*(0.5_r8_kind*(distance12+r_sphere(j)-r_sphere(i))- &
                       ((r_sphere(j)+r_solv)**2+distance12**2-(r_sphere(i)+r_solv)**2)/distance12)
                  if(help<0.0_r8_kind) help=0.0_r8_kind
                  r_new= sqrt(help)-r_solv

                  if ( r_new >= r_mini .or. use_stored) then
                     !case B of GEPOL algorithm
!                    print*,'case B of GEPOL algorithm'
                     one_more_sphere=one_more_sphere+1
                     r_sphere_buf(N_spheres+one_more_sphere)=r_new
                     xyz_sphere_buf(N_spheres+one_more_sphere,:)= &
                          0.5_r8_kind*(xyz_sphere(j,:)+xyz_sphere(i,:))- &
                          (r_sphere(j)-r_sphere(i))*(xyz_sphere(j,:)- &
                          xyz_sphere(i,:))/(2.0_r8_kind*distance12)
                  else
                     !case C of GEPOL algorithm
                     if( r_sphere(j) >= r_sphere(i) ) then
                        help=(r_sphere(j)+r_solv)**2+r_sphere(j)**2- &
                             r_sphere(j)*((r_sphere(j)+r_solv)**2+distance12**2- &
                             (r_sphere(i)+r_solv)**2)/distance12
                        if(help<0.0_r8_kind) help=0.0_r8_kind
                        r_new= sqrt(help)-r_solv
                        if ( r_new < r_mini .and. .not. use_stored ) cycle main_cycle_2
!                       print*,'case C of GEPOL algorithm'
                        one_more_sphere=one_more_sphere+1
                        r_sphere_buf(N_spheres+one_more_sphere)=r_new
                        xyz_sphere_buf(N_spheres+one_more_sphere,:)= &
                             xyz_sphere(j,:)-r_sphere(j)*(xyz_sphere(j,:)-xyz_sphere(i,:))/distance12
                     else

                        help=(r_sphere(i)+r_solv)**2+r_sphere(i)**2- &
                             r_sphere(i)*((r_sphere(i)+r_solv)**2+distance12**2- &
                             (r_sphere(j)+r_solv)**2)/distance12
                        r_new= sqrt(help)-r_solv
                        if ( r_new < r_mini .and. .not. use_stored) cycle main_cycle_2
!                      print*,'last case of GEPOL algorithm'
                        one_more_sphere=one_more_sphere+1
                        r_sphere_buf(N_spheres+one_more_sphere)=r_new
                        xyz_sphere_buf(N_spheres+one_more_sphere,:)= &
                             xyz_sphere(i,:)-r_sphere(i)*(xyz_sphere(i,:)-xyz_sphere(j,:))/distance12
                     endif
                  endif
               endif
               parents_buf(N_spheres+one_more_sphere,1)=i
               parents_buf(N_spheres+one_more_sphere,2)=j

            enddo main_cycle_2
         enddo main_cycle_1

         deallocate(r_sphere, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of r_sphere is failed")
         deallocate(xyz_sphere, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of xyz_sphere is failed")
         deallocate(parents, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of parents is failed")

         N_spheres_old=N_spheres
         N_spheres= N_spheres+one_more_sphere

         allocate(r_sphere(N_spheres), stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of r_sphere is failed(1)")
         allocate(xyz_sphere(N_spheres,3), stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of xyz_sphere is failed(1)")
         allocate(parents(N_spheres,2), stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: allocation of parents is failed(1)")

         r_sphere(:)=r_sphere_buf(1:N_spheres)
         xyz_sphere(:,:)=xyz_sphere_buf(1:N_spheres,:)
         parents(:,:)=parents_buf(1:N_spheres,:)
         do j=1,N_spheres
            do i=1,3
               if(abs(xyz_sphere(j,i)) <= 1.0E-10_r8_kind) xyz_sphere(j,i)=0.0_r8_kind
            enddo
         enddo

         deallocate(r_sphere_buf, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation r_sphere_buf is failed")
         deallocate(xyz_sphere_buf, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation xyz_sphere_buf is failed")
         deallocate(parents_buf, stat=status)
         if ( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of parents_buf is failed")

         if (one_more_sphere == 0 ) exit main_cycle_0
      enddo main_cycle_0

      if(fix_number_add .and. .not. parents_file_exists) then
                call write_fixed_parents(N_spheres,N_atom_spheres,parents)
      else if(fix_number_add) then
                deallocate(fix_par,stat=status)
                if( status /= 0) call error_handler( &
                        "points_on_cavity_surface: deallocation of fix_par is failed")
      endif
      if(output_cavity_data .and. .not.do_gradients) then
         write (output_unit,'(a23,i4,a19)') 'The cavity consists of ',N_spheres,' overlaping spheres'
         write (output_unit,*) 'number       radius(a.u.)                 coordinates(a.u)'
         do i=1,N_spheres
            write (output_unit,'(1x,i4,6x,f11.8,6x,3(f13.9,1x))') i,r_sphere(i),xyz_sphere(i,:)
         enddo
         write (output_unit,*) '---------------------------------------------------'
      endif

    end subroutine calc_cavity
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine read_fixed_parents(f_exists, n_max_par, fix_par)
      use iounitadmin_module
      use filename_module
      implicit none
      logical,intent(out) :: f_exists
      integer, intent(out) :: n_max_par
      integer, pointer :: fix_par(:,:)
      ! *** end of interface ***

      integer(kind=i4_kind) :: st, iounit, ii

      if(gepol==87) then
         inquire(file=trim(inpfile('fix_par')), EXIST=f_exists)
         if(f_exists) then
            iounit = openget_iounit(file=trim(inpfile('fix_par')), form='FORMATTED', &
                 action='READ',status='OLD')
         else
            inquire(file=trim(data_dir)//'/fix_par',EXIST=f_exists)
            if(f_exists) then
               iounit = openget_iounit(file=trim(data_dir)//'/fix_par',form='FORMATTED', &
                    action='READ',status='OLD')
            else
               return
            endif
         endif
      elseif(gepol==93) then
         inquire(file=trim(inpfile('fix_par.93a')), EXIST=f_exists)
         if(f_exists) then
            iounit = openget_iounit(file=trim(inpfile('fix_par.93a')), form='FORMATTED', &
                 action='READ',status='OLD')
         else
            inquire(file=trim(data_dir)//'/fix_par.93a',EXIST=f_exists)
            if(f_exists) then
               iounit = openget_iounit(file=trim(data_dir)//'/fix_par.93a',form='FORMATTED', &
                    action='READ',status='OLD')
            else
               return
            endif
         endif
      end if

      read(iounit, fmt=*, iostat=st) n_max_par
      if ( st .ne. 0 ) call error_handler("read_fixed_parents: read n_max_par failed")

      if(gepol==87) then
         allocate(fix_par(n_max_par+1,2),stat=st)
      elseif(gepol==93) then
         allocate(fix_par(n_max_par+1,3),stat=st)
      end if
      if ( st .ne. 0 ) call error_handler("read_fixed_parents: allocation failed")
      fix_par=0

      do ii=1,n_max_par
         if(gepol==87) then
            read(iounit, fmt=*, iostat=st) fix_par(ii,1:2)
         elseif(gepol==93) then
            read(iounit, fmt=*, iostat=st) fix_par(ii,1:3)
         end if
         if ( st .ne. 0 ) call error_handler("read_fixed_parents: read fix_par failed")
      enddo

      call returnclose_iounit(iounit,status='KEEP')
    end subroutine read_fixed_parents
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine write_fixed_parents(n_sph,n_at_sph,fixed_par)
      use iounitadmin_module
      use filename_module
      implicit none
      integer(kind=i4_kind), intent(in) :: n_sph,n_at_sph,fixed_par(:,:)
      ! *** end of interface ***

      integer(kind=i4_kind) :: iostat, iounit,ii, n_par,iounit1

      n_par=n_sph-n_at_sph

      !
      ! I assume  one can safely overwrite old  "fix_par" files?  When
      ! using  status="NEW"  if  would  fail  when  the  file  already
      ! exists. The  "runpg" script used to  delete output directories
      ! prior to  starting the  program. This is  not always  the case
      ! anymore.
      !
      if(gepol == 87) then
         iounit = openget_iounit(file=trim(data_dir)//'/fix_par',form='FORMATTED', &
              action='WRITE') !,status='NEW')
         iounit1 = openget_iounit(file=trim(outfile('fix_par')), form='FORMATTED', &
              action='WRITE') !,status='NEW')
      elseif(gepol==93) then
         iounit = openget_iounit(file=trim(data_dir)//'/fix_par.93a',form='FORMATTED', &
              action='WRITE') !,status='NEW')
         iounit1 = openget_iounit(file=trim(outfile('fix_par.93a')), form='FORMATTED', &
              action='WRITE') !,status='NEW')
      end if

      write(iounit,iostat=iostat,fmt='(i5)') n_par
      if ( iostat .ne. 0 ) call error_handler("write_fixed_parents: write failed (data_dir)" )
      write(iounit1,iostat=iostat,fmt='(i5)') n_par
      if ( iostat .ne. 0 ) call error_handler("write_fixed_parents: write failed (output dir)" )

      do ii=1,n_par
         if(gepol == 87) then
            write(iounit,iostat=iostat,fmt='(2i5)') fixed_par(ii+n_at_sph,:)
         elseif(gepol==93) then
            write(iounit,iostat=iostat,fmt='(3i5)') fixed_par(ii+n_at_sph,:)
         end if
         if ( iostat .ne. 0 ) call error_handler("write_fixed_parents: write failed (data)" )
         if(gepol == 87) then
            write(iounit1,iostat=iostat,fmt='(2i5)') fixed_par(ii+n_at_sph,:)
         elseif(gepol==93) then
            write(iounit1,iostat=iostat,fmt='(3i5)') fixed_par(ii+n_at_sph,:)
         end if
         if ( iostat .ne. 0 ) call error_handler("write_fixed_parents: write failed (output)" )
      enddo

      call returnclose_iounit(iounit,status='KEEP')
      call returnclose_iounit(iounit1,status='KEEP')
    end subroutine write_fixed_parents
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine calc_cavity_93
      !Solvent excluding surface (Gepol93 algorithm)
      !** End of interface *****************************************
      real(r8_kind), allocatable :: r_sphere_buf(:)
      real(r8_kind), allocatable :: xyz_sphere_buf(:,:)
      integer(i4_kind), allocatable :: parents_buf(:,:)
      integer(i4_kind), allocatable :: pairs(:,:)
      integer(i4_kind) :: n_pairs
      integer(i4_kind) :: sym_pairs(500,2)
      integer(i4_kind), pointer :: fix_par(:,:)
      integer(i4_kind) :: n_par_fix,n_max_par

      real(r8_kind) :: r_mini,r_max,r_min,r_solv,v_d_w_r,s_f,ofac
      integer(i4_kind) :: vdW_index,N_spheres_old,one_more_sphere
      real(r8_kind) :: distance12,dist13,dist23
      real(r8_kind) :: delta_xyz12(3),delta_xyz13(3)
      real(r8_kind) :: xyz_new(3),r_new,dist_1new,dist_2new,dist_3new
      real(r8_kind) :: xyz_surf(3),dist_sk
      real(r8_kind) :: xyz_symm(500,4)
      integer(i4_kind) :: max_i,min_i,n_sym

      logical :: ext,parents_file_exists,use_stored
      integer(i4_kind) :: i,j,k,l,m,kk,status !,n
      real(r8_kind), parameter :: smalll=1.0e-8_r8_kind
#ifdef WITH_EFP
      character(len=1) :: efp_name
#endif

      !converting parameters to radian and atomic units
      r_solv=solvent_radius/ang_au !!!!!!!
      r_mini=rmin_gepol/ang_au  !!!!!!!!!
      if(do_cavitation) then
         s_f=1.0_r8_kind
      else
         s_f=scaled_factor
      endif
      ofac=1.0_r8_kind-2.0_r8_kind*overlap_factor

      allocate(r_sphere_buf(N_max_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of r_sphere_buf is failed (93)")
      allocate(xyz_sphere_buf(N_max_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of xyz_sphere_buf is failed (93)")
      allocate(parents_buf(N_max_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of parents_buf is failed (93)")
      allocate(pairs(N_max_spheres*N_max_spheres,2), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of pairs is failed (93)")

      parents_file_exists=.false.
      if(fix_number_add) call read_fixed_parents(parents_file_exists,n_max_par,fix_par)
      use_stored=(fix_number_add .and. parents_file_exists)
      if(use_stored) write(output_unit,*) "THE CAVITY CONFIGURATIONS IS TAKEN FROM FIX_PAR.93 FILE"

      N_spheres = 0
      do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1_i4_kind) cycle
         N_spheres = N_spheres + unique_atoms(i)%N_equal_atoms
      enddo
#ifdef WITH_EFP
      N_qm_atom_spheres=N_spheres
#endif
      if(with_pc) then
         do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle
            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 100
            end do
            cycle
100         N_spheres = N_spheres + pointcharge_array(i)%N_equal_charges
         end do
      end if

#ifdef WITH_EFP
      if(efp .and. n_efp > 0) then
         if(pointcharge_N > 0) then
            do i=1,pointcharge_N
               if((trim(pointcharge_array(i)%name) /= "O1") .and. &
                  (trim(pointcharge_array(i)%name) /= "H2") .and. &
                  (trim(pointcharge_array(i)%name) /= "H3")) cycle
               N_spheres = N_spheres + pointcharge_array(i)%N_equal_charges
            end do
         else
            ASSERT(pointcharge_N > 0)
         end if
      end if
#endif
      N_atom_spheres=N_spheres

      if(do_gradients .and. .not. VTN) then
         allocate(i_unique_at_index(N_spheres),stat=status)
         if(status/=0) call error_handler("calc cavity: allocation &
              & of i_unique_at_index is failed (93)")
         i_unique_at_index(:)=0
      endif
      if(do_gradients .and. VTN) then
         allocate(center2sphere(N_unique_atoms+pointcharge_N,ua_dim_max),stat=status)
         ASSERT(status==0)
         center2sphere=0
      end if

      !definition of initial spheres
      k=0
      n_uni: do i=1,N_unique_atoms
         if(no_hydrogen_sphere .and. int(unique_atoms(i)%Z)==1) cycle n_uni
         ext=.false.
         if(external_vdwr /= 0) then
            do j=1,external_vdwr
               if(i ==  n_uq_at(j)) then
                  vdW_index=j
                  ext=.true.
               endif
            enddo
         endif
         if(external_vdwr /= 0 .and. ext) then
            v_d_w_r=vdwr(vdW_index)
         else
            vdW_index=int(unique_atoms(i)%Z)
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
               if(smooth==3 .and. vdW_index==1) v_d_w_r=1.125_r8_kind
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            if(v_d_w_r == 0.0_r8_kind ) then
               write (output_unit,*) i,' -unique_atom has no default value of van der Waals radius'
               call error_handler( &
                    "points_on_cavity_surface: one of atoms in the solute has no &
                    & default value of van der Waals radius (See output file). If &
                    & You want to know no default van der Waals radii please specify &
                    & option GET_VDWR")
            endif
         endif

         do j=1,unique_atoms(i)%N_equal_atoms
            k=k+1
            if((.not.hydrogen_no_scale).or.(unique_atoms(i)%Z/=1.0_r8_kind).or.&
                 (do_cavitation)) then
               r_sphere_buf(k)=v_d_w_r*s_f/ang_au  !!!!!!!!!!!
            else
               r_sphere_buf(k)=v_d_w_r/ang_au  !!!!!!!!
            endif
            xyz_sphere_buf(k,:)=unique_atoms(i)%position(:,j) !*ang_au !!!!!!!
            parents_buf(k,:)=0
            if(do_gradients .and. .not. VTN) i_unique_at_index(k)= (i - 1) * ua_dim_max + (j - 1) !i*10+j
            if(do_gradients .and. VTN) center2sphere(i,j)=k
         enddo
      enddo n_uni

      if(with_pc) then
         n_uni_pc: do i=1,pointcharge_N
            if(no_hydrogen_sphere .and. (trim(pointcharge_array(i)%name) == "H")) cycle n_uni_pc

            do j=1,98
               if(trim(pointcharge_array(i)%name) == trim(atom_name(j))) goto 200
            end do
            cycle n_uni_pc
200         vdW_index=j
            if(atom_rad == 0) then
               if(vdW_radius(vdW_index) /= 0.0_r8_kind) then
                  v_d_w_r=vdW_radius(vdW_index)
               else
                  v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
               endif
               if(smooth==3 .and. vdW_index==1) v_d_w_r=1.125_r8_kind
            else if(atom_rad == 1) then
               v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
            end if
            do j=1,pointcharge_array(i)%N_equal_charges
               k=k+1
               if((.not.hydrogen_no_scale).or.(trim(pointcharge_array(i)%name) /= "H").or.&
                    (do_cavitation)) then
                  r_sphere_buf(k)=v_d_w_r*s_f/ang_au  !!!!!
               else
                  r_sphere_buf(k)=v_d_w_r/ang_au
               end if
               xyz_sphere_buf(k,:)=pointcharge_array(i)%position(:,j)
               parents_buf(k,:)=0
               if(do_gradients) i_unique_at_index(k)= (i+N_unique_atoms - 1) * ua_dim_max + (j - 1)
            enddo
         end do n_uni_pc
      end if

#ifdef WITH_EFP
      if(efp .and. n_efp > 0) then
         if(pointcharge_N > 0) then
            do i=1,pointcharge_N
               if((trim(pointcharge_array(i)%name) /= "O1") .and. &
                  (trim(pointcharge_array(i)%name) /= "H2") .and. &
                  (trim(pointcharge_array(i)%name) /= "H3")) cycle
               if(trim(pointcharge_array(i)%name) == "H2" .or. &
                    trim(pointcharge_array(i)%name) == "H3") then
                  vdW_index=1
                  efp_name="H"
               end if
               if(trim(pointcharge_array(i)%name) == "O1") then
                  vdW_index=8
                  efp_name="O"
               end if
               v_d_w_r=vdW_radius(vdW_index)
               do j=1,pointcharge_array(i)%N_equal_charges
                  k=k+1
                  if((.not.hydrogen_no_scale) .or. (efp_name /= "H") .or. (do_cavitation)) then
                     r_sphere_buf(k)=v_d_w_r*s_f/ang_au
                  else
                     r_sphere_buf(k)=v_d_w_r/ang_au
                  end if

                  xyz_sphere_buf(k,:)=pointcharge_array(i)%position(:,j)
                  parents_buf(k,:)=0
                  if(do_gradients) center2sphere(i+N_unique_atoms,j)=k
               end do
            end do
         else
            ASSERT(pointcharge_N > 0)
         end if
      end if
#endif

      do j=1,N_spheres
         do i=1,3
            if(abs(xyz_sphere_buf(j,i)) <= 1.0E-10_r8_kind) xyz_sphere_buf(j,i)=0.0_r8_kind
         enddo
      enddo

      pairs=0; n_pairs=0

      if(ofac == 1.0_r8_kind) goto 1 !$$$$$$$$$$$$$$$$!!!!!!!

      n_par_fix=1
      !start BULK mode: fill internal space by spheres which have no contact to
      !solvent (solvent accessible surface are zero)
      N_spheres_old=0
      bulk: do
         one_more_sphere=0

         !choose the first sphere of the pair
         first: do i=1,N_spheres-1
            l=i+1
            if(l < N_spheres_old+1) l= N_spheres_old+1

            !choose the second sphere of the pair
            second: do j=l,N_spheres

               do k=1,n_pairs
                  if(i==pairs(k,1) .and. j==pairs(k,2)) cycle second
               end do

               if(use_stored ) then
                    if(fix_par(n_par_fix,1)==i .and. &
                       fix_par(n_par_fix,2)==j .and. &
                       fix_par(n_par_fix,3)==1) then
                     else
                        cycle second
                    endif
               endif

               !calculate the distance between the spheres
               delta_xyz12(:) = xyz_sphere_buf(j,:) - xyz_sphere_buf(i,:)
               distance12 = LEN3(delta_xyz12)

               if(.not.use_stored) then
                  !can the solvent pass between spheres? if Yes choose new second sphere
                  if (distance12 >= (r_sphere_buf(i)+r_sphere_buf(j)+2.0_r8_kind*r_solv)) then
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second
                  end if

                  !choose the largest sphere
                  if(r_sphere_buf(i) >= r_sphere_buf(j)) then
                     r_max=r_sphere_buf(i)
                     r_min=r_sphere_buf(j)
                  else
                     r_max=r_sphere_buf(j)
                     r_min=r_sphere_buf(i)
                  end if

                  !check overlapping two spheres. If spheres too close to each other
                  !choose new second sphere
                  if(distance12 <= r_max+r_min*ofac) then
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second
                  end if

                  !compute coordinates, radius of new sphere and its distance to
                  !the first sphere
                  xyz_new=0.5_r8_kind*(xyz_sphere_buf(j,:)+xyz_sphere_buf(i,:))- &
                       (r_sphere_buf(j)-r_sphere_buf(i))*(xyz_sphere_buf(j,:)-xyz_sphere_buf(i,:))/ &
                       (2.0_r8_kind*distance12)
                  r_new=r_max
                  dist_1new=0.5_r8_kind*(distance12-r_sphere_buf(j)+r_sphere_buf(i))


                  !check overlapping the new sphere
                  third: do k=1,N_spheres+one_more_sphere
                     delta_xyz13(:) = xyz_new(:) - xyz_sphere_buf(k,:)
                     dist_3new = LEN3(delta_xyz13)

                     !if no overlapping choose new the third sphere
                     if(dist_3new >= r_sphere_buf(k)+r_new) cycle third
                     !if overlapping
                     if(r_new > r_sphere_buf(k)) cycle third
                     !if overlapping too strong cancel the process and choose new second sphere
                     if(dist_3new <= r_sphere_buf(k)+r_new*ofac) then
                        sym_pairs=0
                        call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                             sym_pairs,n_sym)
                        do kk=1,n_sym
                           n_pairs=n_pairs+1
                           pairs(n_pairs,:)=sym_pairs(kk,:)
                        end do
                        cycle second
                     end if
                  end do third

                  !Determine if new sphere contact with solvent
                  !(has non zero solvent accessible surface)
                  !if Yes choose new second sphere (define new pair)
                  surf_points:do k=1,N_centers_on_sphere
                     xyz_surf=surf_elem%xyz_centers(k,:)*(r_new+r_solv)/radius+xyz_new

                     atoms:do m=1,N_atom_spheres
                        delta_xyz13(:) = xyz_surf(:) - xyz_sphere_buf(m,:)
                        dist_sk = LEN3(delta_xyz13)

                        if(dist_sk < r_sphere_buf(m)+r_solv) cycle surf_points
                     end do atoms
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second
                  end do surf_points
               else
                  !choose the largest sphere
                  if(r_sphere_buf(i) >= r_sphere_buf(j)) then
                     r_max=r_sphere_buf(i)
                     r_min=r_sphere_buf(j)
                  else
                     r_max=r_sphere_buf(j)
                     r_min=r_sphere_buf(i)
                  end if
                  !compute coordinates, radius of new sphere and its distance to
                  !the first sphere
                  xyz_new=0.5_r8_kind*(xyz_sphere_buf(j,:)+xyz_sphere_buf(i,:))- &
                       (r_sphere_buf(j)-r_sphere_buf(i))*(xyz_sphere_buf(j,:)-xyz_sphere_buf(i,:))/ &
                       (2.0_r8_kind*distance12)
                  r_new=r_max
               end if

               sym_pairs=0
               call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                    sym_pairs,n_sym,xyz_new,xyz_symm)
               l_kk1:do kk=1,n_sym
                  n_pairs=n_pairs+1
                  pairs(n_pairs,:)=sym_pairs(kk,:)
                  if(int(xyz_symm(kk,4))==0) cycle
                  one_more_sphere=one_more_sphere+1
                  xyz_sphere_buf(N_spheres+one_more_sphere,:)=xyz_symm(kk,1:3)
                  r_sphere_buf(N_spheres+one_more_sphere)=r_new
                  parents_buf(N_spheres+one_more_sphere,1)=sym_pairs(kk,1)
                  parents_buf(N_spheres+one_more_sphere,2)=sym_pairs(kk,2)
                  parents_buf(N_spheres+one_more_sphere,3)=1
                  if(use_stored ) n_par_fix=n_par_fix+1
               end do l_kk1
            end do second
         end do first

         N_spheres_old=N_spheres
         N_spheres=N_spheres+one_more_sphere
         if(N_spheres > N_max_spheres) call error_handler("calc_cavity 93: number of &
              & spheres have exceeded maximal value 1000")

         if(one_more_sphere == 0) exit bulk

      end do bulk

      N_bulk_spheres=N_spheres

      do j=1,N_spheres
         do i=1,3
            if(abs(xyz_sphere_buf(j,i)) <= 1.0E-10_r8_kind) xyz_sphere_buf(j,i)=0.0_r8_kind
         enddo
      enddo

      pairs=0; n_pairs=0

      !start CREA mode: refining the surface built after BULK mode run
      N_spheres_old=0
      crea: do
         one_more_sphere=0

         !choose the first sphere of the pair
         first1: do i=1,N_spheres-1
            l=i+1
            if(l < N_spheres_old+1) l= N_spheres_old+1

            !choose the second sphere of the pair
            second1: do j=l,N_spheres
               do k=1,n_pairs
                  if(i==pairs(k,1) .and. j==pairs(k,2)) cycle second1
               end do

               if(use_stored ) then
                    if(fix_par(n_par_fix,1)==i .and. &
                       fix_par(n_par_fix,2)==j .and. &
                       fix_par(n_par_fix,3)==2) then
                    else
                        cycle second1
                    endif
               endif

               !calculate the distance between the spheres
               delta_xyz12(:) = xyz_sphere_buf(j,:) - xyz_sphere_buf(i,:)
               distance12 = LEN3(delta_xyz12)

               if(.not.use_stored) then
                  !can the solvent pass between spheres? if Yes choose new second sphere
                  if (distance12 >= (r_sphere_buf(i)+r_sphere_buf(j)+2.0_r8_kind*r_solv)) then
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second1
                  end if

                  !choose the largest sphere
                  if(r_sphere_buf(i) >= r_sphere_buf(j)) then
                     r_max=r_sphere_buf(i)
                     r_min=r_sphere_buf(j)
                     max_i=i; min_i=j
                  else
                     r_max=r_sphere_buf(j)
                     r_min=r_sphere_buf(i)
                     max_i=j; min_i=i
                  end if

                  !check overlapping two spheres. If spheres too close to each other
                  !choose new second sphere
                  if(distance12 < r_max+r_min*ofac) then
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second1
                  end if
               else
                  !choose the largest sphere
                  if(r_sphere_buf(i) >= r_sphere_buf(j)) then
                     r_max=r_sphere_buf(i)
                     r_min=r_sphere_buf(j)
                     max_i=i; min_i=j
                  else
                     r_max=r_sphere_buf(j)
                     r_min=r_sphere_buf(i)
                     max_i=j; min_i=i
                  end if
               end if

               dist13=r_min+r_solv
               dist23=r_max+r_solv
               dist_2new=0.5_r8_kind*(distance12-r_min+r_max)

               if(distance12 <= r_sphere_buf(i)+r_sphere_buf(j)) then
                  !variant A of sphere building
                  r_new=sqrt(dist23**2+dist_2new* &
                       (dist_2new-(dist23**2+distance12**2-dist13**2)/distance12))- &
                       r_solv
                  if(r_new <= r_mini) then
                     sym_pairs=0
                     call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                          sym_pairs,n_sym)
                     do kk=1,n_sym
                        n_pairs=n_pairs+1
                        pairs(n_pairs,:)=sym_pairs(kk,:)
                     end do
                     cycle second1
                  end if
                  xyz_new=0.5_r8_kind*(xyz_sphere_buf(max_i,:)+xyz_sphere_buf(min_i,:))- &
                       (r_max-r_min)*(xyz_sphere_buf(max_i,:)-xyz_sphere_buf(min_i,:))/ &
                       (2.0_r8_kind*distance12)
                  dist_1new=0.5_r8_kind*(distance12-r_sphere_buf(j)+r_sphere_buf(i))
               else
                  !variants B and C
                  r_new=sqrt(dist23**2+dist_2new* &
                       (dist_2new-(dist23**2+distance12**2-dist13**2)/distance12))- &
                       r_solv
                  if(r_new > r_mini) then
                     !variant B
                     xyz_new=0.5_r8_kind*(xyz_sphere_buf(max_i,:)+xyz_sphere_buf(min_i,:))- &
                          (r_max-r_min)*(xyz_sphere_buf(max_i,:)-xyz_sphere_buf(min_i,:))/ &
                          (2.0_r8_kind*distance12)
                     dist_1new=0.5_r8_kind*(distance12-r_sphere_buf(j)+r_sphere_buf(i))
                  else
                     !variant C
                     r_new=sqrt(dist23**2+r_max**2- &
                          (r_max/distance12)*(dist23**2+distance12**2-dist13**2))- &
                          r_solv
                     if(r_new <= r_mini) then
                        sym_pairs=0
                        call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                             sym_pairs,n_sym)
                        do kk=1,n_sym
                           n_pairs=n_pairs+1
                           pairs(n_pairs,:)=sym_pairs(kk,:)
                        end do
                        cycle second1
                     end if
!!$                     xyz_new=0.5_r8_kind*(xyz_sphere_buf(max_i,:)+xyz_sphere_buf(min_i,:))- &
!!$                          (r_max-r_min)*(xyz_sphere_buf(max_i,:)-xyz_sphere_buf(min_i,:))/ &
!!$                          (2.0_r8_kind*distance12)
                     xyz_new=xyz_sphere_buf(max_i,:)- &
                          r_max*(xyz_sphere_buf(max_i,:)-xyz_sphere_buf(min_i,:))/distance12
                     if(max_i == i) then
                        dist_1new=dist_2new
                     else
                        dist_1new=distance12-dist_2new
                     end if
                  end if
               end if

               if(.not.use_stored) then
                  !check overlapping the new sphere
                  third1: do k=1,N_spheres+one_more_sphere
                     delta_xyz13(:) = xyz_new(:) - xyz_sphere_buf(k,:)
                     dist_3new = LEN3(delta_xyz13)

                     !if no overlapping choose new the third sphere
                     if(dist_3new >= r_sphere_buf(k)+r_new) cycle third1
                     !if overlapping
                     if(r_new > r_sphere_buf(k) .and.  &
                          abs(r_new-r_sphere_buf(k)) > smalll) cycle third1
                     !if overlapping too strong cancel the process and choose new second sphere
                     if(dist_3new <= r_sphere_buf(k)+r_new*ofac) then
                        sym_pairs=0
                        call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                             sym_pairs,n_sym)
                        do kk=1,n_sym
                           n_pairs=n_pairs+1
                           pairs(n_pairs,:)=sym_pairs(kk,:)
                        end do
                        cycle second1
                     end if
                  end do third1
               end if

               sym_pairs=0
               call calc_symm_pairs(i,j,distance12,xyz_sphere_buf(1:N_spheres,1:3),N_spheres, &
                    sym_pairs,n_sym,xyz_new,xyz_symm)
               l_kk:do kk=1,n_sym
                  n_pairs=n_pairs+1
                  pairs(n_pairs,:)=sym_pairs(kk,:)
                  if(int(xyz_symm(kk,4))==0) cycle
                  one_more_sphere=one_more_sphere+1
                  xyz_sphere_buf(N_spheres+one_more_sphere,:)=xyz_symm(kk,1:3)
                  r_sphere_buf(N_spheres+one_more_sphere)=r_new
                  parents_buf(N_spheres+one_more_sphere,1)=sym_pairs(kk,1)
                  parents_buf(N_spheres+one_more_sphere,2)=sym_pairs(kk,2)
                  parents_buf(N_spheres+one_more_sphere,3)=2
                  if(use_stored ) n_par_fix=n_par_fix+1
               end do l_kk
            end do second1
         end do first1

         N_spheres_old=N_spheres
         N_spheres=N_spheres+one_more_sphere
         if(N_spheres > N_max_spheres) call error_handler("calc_cavity 93: number of &
              & spheres have exceeded maximal value 1000 (1)")

         if(one_more_sphere == 0) exit crea

      end do crea

1     allocate(r_sphere(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of r_sphere is failed (93)")
      allocate(zero_area(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of zero_area is failed (93)")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of xyz_sphere is failed (93)")
      allocate(parents(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of parents is failed (93)")

      r_sphere(1:N_spheres)=r_sphere_buf(1:N_spheres)
      xyz_sphere(1:N_spheres,:)=xyz_sphere_buf(1:N_spheres,:)
      parents(1:N_spheres,:)=parents_buf(1:N_spheres,:)

      do j=1,N_spheres
         do i=1,3
            if(abs(xyz_sphere(j,i)) <= 1.0E-10_r8_kind) xyz_sphere(j,i)=0.0_r8_kind
         enddo
      enddo

      !marking spheres with zero contribution to solvent excluding surface
      zero_area=.true.
      one:do i=1,N_spheres
         three:do k=1,N_centers_on_sphere
            xyz_surf=surf_elem%xyz_centers(k,:)*r_sphere(i)/radius+xyz_sphere(i,:)
            two:do j=1,N_spheres
               if(j==i) cycle two

               delta_xyz13(:) = xyz_surf(:) - xyz_sphere(j,:)
               dist_sk = LEN3(delta_xyz13)

               if(dist_sk <= r_sphere_buf(j)) cycle three
            end do two
            zero_area(i)=.false.
            cycle one
         enddo three
      enddo one

      deallocate(r_sphere_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation of r_sphere_buf is failed (93)")
      deallocate(xyz_sphere_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation of xyz_sphere_buf is failed (93)")
      deallocate(parents_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation of parents_buf is failed (93)")
      deallocate(pairs, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation of pairs is failed (93)")

      if(fix_number_add .and. .not. parents_file_exists) then
         call write_fixed_parents(N_spheres,N_atom_spheres,parents)
      else if(fix_number_add) then
         deallocate(fix_par,stat=status)
         if( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of fix_par is failed (93)")
      end if

      if(output_cavity_data .and. .not.do_gradients) then
         write (output_unit,'(a23,i4,a19)') 'The cavity consists of ',N_spheres,' overlaping spheres'
         write (output_unit,*) 'number       radius(a.u.)                 coordinates(a.u)'
         do i=1,N_spheres
            write (output_unit,'(1x,i4,6x,f11.8,6x,3(f13.9,1x),l1)') i,r_sphere(i),xyz_sphere(i,:),zero_area(i)
         enddo
         write (output_unit,*) '---------------------------------------------------'
      endif

    end subroutine calc_cavity_93
    !------------------------------------------------------------
!!!MF - gradients >>>>
    !------------------------------------------------------------
    subroutine calc_cavity_gradient(cavity_flag)
      ! formulas see Cossi 1996
      use constants
#ifdef WITH_EFP
      use efp_data_module, only: n_nuc
#endif
      !** End of interface *****************************************
      logical, intent(in) :: cavity_flag ! spheres with parents?

      integer(kind=i4_kind) :: i,j,k,l,l0,l1,l2,l3,status,na,ea
      integer(i4_kind) :: m,m0,m1,m2,m3,na1,ea1
      real(kind=r8_kind) :: r1s,r2s,ras,d12,vbuf(3),help,da1,da2,da3
      real(kind=r8_kind) :: drdr1,drdr2,drda(3),dadr(3),dbdadiag1,dbdadiag2,dada1(3,3),dada2(3,3), &
           dbdanondiag,r_solv

      real(r8_kind) :: d2r_dr1dr1,d2r_dr1dr2,d2r_dr2dr2,d2r_dr2dr1
      real(r8_kind) :: d2r_dr1da1(3),d2r_dr1da2(3),d2r_dr2da1(3),d2r_dr2da2(3)
      real(r8_kind) :: d2r_da1dr1(3),d2r_da1dr2(3),d2r_da2dr1(3),d2r_da2dr2(3)
      real(r8_kind) :: d2r_da1da1(3,3),d2r_da1da2(3,3),d2r_da2da1(3,3),d2r_da2da2(3,3)
      real(r8_kind) :: d2a_dr1dr1(3),d2a_dr1dr2(3),d2a_dr2dr1(3),d2a_dr2dr2(3)
      real(r8_kind) :: d2a_dr1da1(3,3),d2a_dr1da2(3,3),d2a_dr2da1(3,3),d2a_dr2da2(3,3)
      real(r8_kind) :: d2a_da1dr1(3,3),d2a_da1dr2(3,3),d2a_da2dr1(3,3),d2a_da2dr2(3,3)
      real(r8_kind) :: d2a_da1da1(3,3,3),d2a_da1da2(3,3,3),d2a_da2da1(3,3,3),d2a_da2da2(3,3,3)
      real(r8_kind) :: dr1_dai(3),dr1_daj(3),da1_dai(3,3),da1_daj(3,3)
      real(r8_kind) :: dr2_dai(3),dr2_daj(3),da2_dai(3,3),da2_daj(3,3)
      real(r8_kind) :: d2r1_daidaj(3,3),d2a1_daidaj(3,3,3)
      real(r8_kind) :: d2r2_daidaj(3,3),d2a2_daidaj(3,3,3)
      real(r8_kind) :: d2r_dr1daj(3),d2r_dr2daj(3)
      real(r8_kind) :: d2r_da1daj(3,3),d2r_da2daj(3,3)
      real(r8_kind) :: d2a_dr1daj(3,3),d2a_dr2daj(3,3)
      real(r8_kind) :: d2a_da1daj(3,3,3),d2a_da2daj(3,3,3)

      logical :: dist_check_j,dist_check_k
      integer :: N_pc
      real(kind=r8_kind), parameter :: small=1.e-8_r8_kind
      real(i8_kind) :: unit(3,3)=reshape((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind, &
                                           0.0_r8_kind,1.0_r8_kind,0.0_r8_kind, &
                                           0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/), &
                                         (/3,3/))

      if(VTN) return

      N_pc=0
!!$      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
      if(with_pc) N_pc=pointcharge_N

      DPRINT 'calc_cavity_gradient(',cavity_flag,'): enetered'
      r_solv=solvent_radius/ang_au
      call alloc_geom_deriv_part1(cagr)

      ! gradients of atomic centers and radia
      ASSERT(ua_dim_max.ne.-1)
      do i=1,N_atom_spheres
         DPRINT '_gradient: i_unique_at_index(',i,')=',i_unique_at_index(i)
         na= 1 +     i_unique_at_index(i) / ua_dim_max
         ea= 1 + mod(i_unique_at_index(i),  ua_dim_max)
         DPRINT '_gradient: na=',na,' ea=',ea
         do j=1,3
            cagr%dc(j,i)%xyz_grad(j,na,ea)=1.0_r8_kind
         enddo
      enddo

      deallocate(i_unique_at_index,stat=status)
      if(status/=0) call error_handler("calc_cavity_gradient: &
                & deallocation of i_unique_at_index is failed")

      ! gradients of additional centers and radii
      ! spheres are ordered such that parents always before childs
      cav_f: if(cavity_flag) then
         do i=N_atom_spheres+1,N_spheres
            j=parents(i,1)
            k=parents(i,2)
            if(j>=i .or. k>=i) call error_handler("calc_cavity_gradient: &
                 & wrong order of spheres")

            vbuf=xyz_sphere(j,:)-xyz_sphere(k,:)

            d12=LEN3(vbuf)
            r1s=r_sphere(j)+r_solv
            r2s=r_sphere(k)+r_solv
            ras=r_sphere(i)+r_solv

            vbuf=xyz_sphere(j,:)-r_sphere(j)/d12*(xyz_sphere(j,:)-xyz_sphere(k,:))&
                 -xyz_sphere(i,:)
            dist_check_j = DOT3(vbuf,vbuf)<small
            if(dist_check_j) then
               dist_check_k = .false.
            else
               vbuf=xyz_sphere(k,:)-r_sphere(k)/d12* &
                    (xyz_sphere(k,:)-xyz_sphere(j,:))&
                    -xyz_sphere(i,:)
               dist_check_k = DOT3(vbuf,vbuf)<small
            endif

            if((.not.dist_check_j .and. .not. dist_check_k) .or. i <= N_bulk_spheres) then
               ! case A/B formulas (18)-(21)
               if(i > N_bulk_spheres) then
                  drdr1=(3.0_r8_kind*r1s*(d12-r1s)-r2s*(d12-r2s)&
                       +2.0_r8_kind*r1s*r2s)/(4.0_r8_kind*d12*ras)
                  drdr2=(3.0_r8_kind*r2s*(d12-r2s)-r1s*(d12-r1s)&
                       +2.0_r8_kind*r2s*r1s)/(4.0_r8_kind*d12*ras)
                  drda(:)= (-d12**3+(r1s+r2s)*(r1s-r2s)**2)/&
                       (4.0_r8_kind*d12**3*ras) &
                       *(xyz_sphere(j,:)-xyz_sphere(k,:))
               else
                  drdr1=zero; drdr2=zero; drda=zero
               end if
               dadr(:)=-(xyz_sphere(j,:)-xyz_sphere(k,:))/(2.0_r8_kind*d12)
               dbdadiag1= 0.5_r8_kind - (r_sphere(j)-r_sphere(k))/(2.0_r8_kind*d12)
               dbdadiag2= 0.5_r8_kind + (r_sphere(j)-r_sphere(k))/(2.0_r8_kind*d12)
               dbdanondiag = (r_sphere(j)-r_sphere(k))/(2.0_r8_kind*d12**3)

               if(integralpar_2dervs)then
                  if(i <= N_bulk_spheres) then
                     d2r_dr1dr1=zero; d2r_dr1dr2=zero; d2r_dr2dr1=zero; d2r_dr2dr2=zero
                     d2r_dr1da1=zero; d2r_dr1da2=zero; d2r_dr2da1=zero; d2r_dr2da2=zero
                     d2r_da1dr1=zero; d2r_da1dr2=zero; d2r_da2dr1=zero; d2r_da2dr2=zero
                     d2r_da1da1=zero; d2r_da1da2=zero; d2r_da2da1=zero; d2r_da2da2=zero
                  else
                     d2r_dr1dr1=(three*d12-six*r1s+two*r2s)/(four*d12*ras)-drdr1**2/ras
                     d2r_dr1dr2=(two*r2s-d12+two*r1s)/(four*d12*ras)-drdr1*drdr2/ras
                     d2r_dr2dr1=(two*r1s-d12+two*r2s)/(four*d12*ras)-drdr2*drdr1/ras
                     d2r_dr2dr2=(three*d12-six*r2s+two*r1s)/(four*d12*ras)-drdr2**2/ras

                     d2r_dr1da1=(three*r1s-r2s)*(xyz_sphere(j,:)-xyz_sphere(k,:))/(four*d12**2*ras)- &
                          drdr1*(xyz_sphere(j,:)-xyz_sphere(k,:))/d12**2-drdr1*drda/ras
                     d2r_dr1da2=-d2r_dr1da1
                     d2r_dr2da1=(three*r2s-r1s)*(xyz_sphere(j,:)-xyz_sphere(k,:))/(four*d12**2*ras)- &
                          drdr2*(xyz_sphere(j,:)-xyz_sphere(k,:))/d12**2+drdr2*drda/ras
                     d2r_dr2da2=-d2r_dr2da1

                     d2r_da1dr1=(xyz_sphere(j,:)-xyz_sphere(k,:))*((r1s-r2s)**2+two*r1s**2-two*r2s**2)/ &
                          (four*d12**3*ras)-drda*drdr1/ras
                     d2r_da1dr2=(xyz_sphere(j,:)-xyz_sphere(k,:))*((r2s-r1s)**2+two*r2s**2-two*r1s**2)/ &
                          (four*d12**3*ras)-drda*drdr2/ras
                     d2r_da2dr1=-d2r_da1dr1
                     d2r_da2dr2=-d2r_da1dr2

                     do l1=1,3
                        da1=xyz_sphere(j,l1)-xyz_sphere(k,l1)
                        do l2=1,3
                           da2=xyz_sphere(j,l2)-xyz_sphere(k,l2)

                           d2r_da1da1(l1,l2)=-(three/four)*(da1*da2)/(d12**2*ras)-three*drda(l1)*da2/d12**2- &
                                drda(l1)*drda(l2)/ras
                           if(da1>= 1.0e-10_r8_kind) d2r_da1da1(l1,l2)=d2r_da1da1(l1,l2)+drda(l1)*unit(l1,l2)/da1
                           d2r_da1da2(l1,l2)=-d2r_da1da1(l1,l2)
                           d2r_da2da1(l1,l2)=-d2r_da1da1(l1,l2)
                           d2r_da2da2(l1,l2)=d2r_da1da1(l1,l2)
                        end do
                     end do
                  end if

                  dada1(:,:)=(half-(r_sphere(j)-r_sphere(k))/(two*d12))*unit(:,:)
                  dada2(:,:)=(half-(r_sphere(j)+r_sphere(k))/(two*d12))*unit(:,:)

                  d2a_dr1dr1=zero; d2a_dr1dr2=zero; d2a_dr2dr1=zero; d2a_dr2dr2=zero

                  do l1=1,3
                     da1=xyz_sphere(j,l1)-xyz_sphere(k,l1)
                     do l2=1,3
                        da2=xyz_sphere(j,l2)-xyz_sphere(k,l2)


                        dada1(l1,l2)=dada1(l1,l2)+(r_sphere(j)-r_sphere(k))*da1*da2/(two*d12**3)
                        dada2(l1,l2)=dada2(l1,l2)-(r_sphere(j)-r_sphere(k))*da1*da2/(two*d12**3)

                        d2a_dr1da1(l1,l2)=-unit(l1,l2)/(two*d12)+da1*da2/(two*d12**3)
                        d2a_dr1da2(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_dr2da1(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_dr2da2(l1,l2)=d2a_dr1da1(l1,l2)

                        d2a_da1dr1(l1,l2)=d2a_dr1da1(l1,l2)
                        d2a_da1dr2(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_da2dr1(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_da2dr2(l1,l2)=d2a_dr1da1(l1,l2)

                        do l3=1,3
                           da3=xyz_sphere(j,l3)-xyz_sphere(k,l3)

                           d2a_da1da1(l1,l2,l3)=(r_sphere(j)-r_sphere(k))*unit(l1,l2)*da3/(two*d12**3)+ &
                                (r_sphere(j)-r_sphere(k))*da2*unit(l1,l3)/(two*d12**3)+ &
                                (r_sphere(j)-r_sphere(k))*da1*unit(l2,l3)/(two*d12**3)- &
                                three*(r_sphere(j)-r_sphere(k))*da1*da2*da3/(two*d12**5)
                           d2a_da1da2(l1,l2,l3)=-d2a_da1da1(l1,l2,l3)
                           d2a_da2da1(l1,l2,l3)=-d2a_da1da1(l1,l2,l3)
                           d2a_da2da2(l1,l2,l3)=d2a_da1da1(l1,l2,l3)
                        end do
                     end do
                  end do
               end if

               do l=1,N_moving_unique_atoms+N_pc
                  if(l <= N_moving_unique_atoms) then
                     na=moving_unique_atom_index(l)
                     ea=unique_atoms(na)%n_equal_atoms
                  else
                     na=l-N_moving_unique_atoms
                     if(with_pc) ea=pointcharge_array(na)%N_equal_charges
                  end if
                  do l0=1,ea
                     do l1=1,3
                        help=xyz_sphere(j,l1)-xyz_sphere(k,l1) !beta_1-beta_2
                        do l2=1,3
                           cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                                cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                                + dbdanondiag*help* &
                                (xyz_sphere(j,l2)-xyz_sphere(k,l2)) * &
                                cagr%dc(l2,j)%xyz_grad(:,l,l0) &
                                + dbdanondiag*help* &
                                (xyz_sphere(k,l2)-xyz_sphere(j,l2)) * &
                                cagr%dc(l2,k)%xyz_grad(:,l,l0)
                        enddo

                        cagr%dR(i)%xyz_grad(:,l,l0) = cagr%dR(i)%xyz_grad(:,l,l0) &
                             + drda(l1)*cagr%dc(l1,j)%xyz_grad(:,l,l0) &
                             - drda(l1)*cagr%dc(l1,k)%xyz_grad(:,l,l0)

                        cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                             + dadr(l1)*cagr%dR(j)%xyz_grad(:,l,l0) &
                             - dadr(l1)*cagr%dR(k)%xyz_grad(:,l,l0)

                        cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                             + dbdadiag1*cagr%dc(l1,j)%xyz_grad(:,l,l0) &
                             + dbdadiag2*cagr%dc(l1,k)%xyz_grad(:,l,l0)

                     enddo
                     cagr%dR(i)%xyz_grad(:,l,l0) = cagr%dR(i)%xyz_grad(:,l,l0) + &
                          drdr1*cagr%dR(j)%xyz_grad(:,l,l0) + &
                          drdr2*cagr%dR(k)%xyz_grad(:,l,l0)
!!$print*,'A-B dr/dx  ',i,l,l0,cagr%dR(i)%xyz_grad(:,l,l0)
!!$print*,'A-B dx/dx x',i,l,l0,cagr%dc(1,i)%xyz_grad(:,l,l0)
!!$print*,'A-B dx/dx y',i,l,l0,cagr%dc(2,i)%xyz_grad(:,l,l0)
!!$print*,'A-B dx/dx z',i,l,l0,cagr%dc(3,i)%xyz_grad(:,l,l0)
                  enddo
               enddo
            else
               ! case C
               ! test if j/k resp r1 and r2 have to been swaped !!!
               if(dist_check_k) then
                  l=k
                  k=j
                  j=l
                  help=r1s
                  r1s=r2s
                  r2s=help
               endif

               drdr1=(2.0_r8_kind*(d12*r1s+d12*r_sphere(j)- &             !??????
                    r1s*r_sphere(j)) - r1s**2 -r2s**2 +d12**2)/ &
                    (2.0_r8_kind*d12*ras)
               drdr2= (r2s*r_sphere(j))/(d12*ras)
               drda(:) = r_sphere(j)*(-d12**2+r1s**2-r2s**2)/ &
                    (2.0_r8_kind*d12**3*ras) * &
                    (xyz_sphere(j,:)-xyz_sphere(k,:))
               dadr(:)=-(xyz_sphere(j,:)-xyz_sphere(k,:))/d12
               dbdadiag2=  r_sphere(j)/d12
               dbdadiag1= 1.0_r8_kind - dbdadiag2
               dbdanondiag = r_sphere(j)/d12**3

               if(integralpar_2dervs)then
                  d2r_dr1dr1=(two*d12-two*r1s-r_sphere(j))/(d12*ras)- &
                       drdr1**2/ras
                  d2r_dr1dr2=-r2s/(d12*ras)-drdr1*drdr2/ras
                  d2r_dr2dr1= r2s/(d12*ras)-drdr2*drdr1/ras
                  d2r_dr2dr2= r_sphere(j)/(d12*ras)-drdr2**2/ras

                  d2r_dr1da1=(r1s+r_sphere(j)+two*d12)*(xyz_sphere(j,:)-xyz_sphere(k,:))/(d12**2*ras)- &
                       drdr1*(xyz_sphere(j,:)-xyz_sphere(k,:))/d12**2-drdr1*drda/ras
                  d2r_dr1da2=-d2r_dr1da1
                  d2r_dr2da1=-r2s*r_sphere(j)/(d12**3*ras)-drdr2*drda/ras
                  d2r_dr2da2=-d2r_dr2da1

                  d2r_da1dr1=(xyz_sphere(j,:)-xyz_sphere(k,:))*(-d12**2+r1s**2-r2s**2)/(two*d12**3*ras)+ &
                       (xyz_sphere(j,:)-xyz_sphere(k,:))*r_sphere(j)*r1s/(d12**3*ras)-drda*drdr1/ras
                  d2r_da1dr2=-(xyz_sphere(j,:)-xyz_sphere(k,:))*r_sphere(j)*r2s/(d12**3*ras)-drda*drdr1/ras
                  d2r_da2dr1=-d2r_da1dr1
                  d2r_da2dr2=-d2r_da1dr2

                  do l1=1,3
                     da1=xyz_sphere(j,l1)-xyz_sphere(k,l1)
                     do l2=1,3
                        da2=xyz_sphere(j,l2)-xyz_sphere(k,l2)

                        d2r_da1da1(l1,l2)=-r_sphere(j)*(da1*da2)/(d12**3*ras)-three*drda(l1)*da2/d12**2- &
                             drda(l1)*drda(l2)/ras
                        if(da1>= 1.0e-10_r8_kind) d2r_da1da1(l1,l2)=d2r_da1da1(l1,l2)+drda(l1)*unit(l1,l2)/da1
                        d2r_da1da2(l1,l2)=-d2r_da1da1(l1,l2)
                        d2r_da2da1(l1,l2)=-d2r_da1da1(l1,l2)
                        d2r_da2da2(l1,l2)=d2r_da1da1(l1,l2)
                     end do
                  end do

                  d2a_dr1dr1=zero; d2a_dr1dr2=zero; d2a_dr2dr1=zero; d2a_dr2dr2=zero

                  dada1(:,:)=(one-r_sphere(j)/d12)*unit(:,:)
                  dada2(:,:)=(r_sphere(j)/d12)*unit(:,:)

                  do l1=1,3
                     da1=xyz_sphere(j,l1)-xyz_sphere(k,l1)
                     do l2=1,3
                        da2=xyz_sphere(j,l2)-xyz_sphere(k,l2)

                        dada1(l1,l2)=dada1(l1,l2)+r_sphere(j)*da1*da2/d12**3
                        dada2(l1,l2)=dada2(l1,l2)-r_sphere(j)*da1*da2/d12**3

                        d2a_dr1da1(l1,l2)=-unit(l1,l2)/d12+da1*da2/d12**3
                        d2a_dr1da2(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_dr2da1(l1,l2)=zero
                        d2a_dr2da2(l1,l2)=zero

                        d2a_da1dr1(l1,l2)=d2a_dr1da1(l1,l2)
                        d2a_da1dr2(l1,l2)=zero
                        d2a_da2dr1(l1,l2)=-d2a_dr1da1(l1,l2)
                        d2a_da2dr2(l1,l2)=zero

                        do l3=1,3
                           da3=xyz_sphere(j,l3)-xyz_sphere(k,l3)

                           d2a_da1da1(l1,l2,l3)=r_sphere(j)*unit(l1,l2)*da3/d12**3+ &
                                r_sphere(j)*da2*unit(l1,l3)/d12**3+ &
                                r_sphere(j)*da1*unit(l2,l3)/d12**3- &
                                three*r_sphere(j)*da1*da2*da3/d12**5
                           d2a_da1da2(l1,l2,l3)=-d2a_da1da1(l1,l2,l3)
                           d2a_da2da1(l1,l2,l3)=-d2a_da1da1(l1,l2,l3)
                           d2a_da2da2(l1,l2,l3)=d2a_da1da1(l1,l2,l3)
                        end do
                     end do
                  end do
               end if

               do l=1,N_moving_unique_atoms+N_pc
                  if(l <= N_moving_unique_atoms) then
                     na=moving_unique_atom_index(l)
                     ea=unique_atoms(na)%n_equal_atoms
                  else
                     na=l-N_moving_unique_atoms
                     if(with_pc) ea=pointcharge_array(na)%N_equal_charges
                  end if
                  do l0=1,ea
                     do l1=1,3
                        help=xyz_sphere(j,l1)-xyz_sphere(k,l1) !beta_1-beta_2
                        do l2=1,3
                           cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                                cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                                + dbdanondiag*help* &
                                (xyz_sphere(j,l2)-xyz_sphere(k,l2)) * &
                                cagr%dc(l2,j)%xyz_grad(:,l,l0) &
                                + dbdanondiag*help* &
                                (xyz_sphere(k,l2)-xyz_sphere(j,l2)) * &
                                cagr%dc(l2,k)%xyz_grad(:,l,l0)
                        enddo

                        cagr%dR(i)%xyz_grad(:,l,l0) = cagr%dR(i)%xyz_grad(:,l,l0) &
                             + drda(l1)*cagr%dc(l1,j)%xyz_grad(:,l,l0) &
                             - drda(l1)*cagr%dc(l1,k)%xyz_grad(:,l,l0)

                        cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                             + dadr(l1)*cagr%dR(j)%xyz_grad(:,l,l0)

                        cagr%dc(l1,i)%xyz_grad(:,l,l0) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,l0) &
                             + dbdadiag1*cagr%dc(l1,j)%xyz_grad(:,l,l0) &
                             + dbdadiag2*cagr%dc(l1,k)%xyz_grad(:,l,l0)

                     enddo
                     cagr%dR(i)%xyz_grad(:,l,l0) = cagr%dR(i)%xyz_grad(:,l,l0) + &
                          drdr1*cagr%dR(j)%xyz_grad(:,l,l0) + &
                          drdr2*cagr%dR(k)%xyz_grad(:,l,l0)
!!$print*,'C   dr/dx  ',i,l,l0,cagr%dR(i)%xyz_grad(:,l,l0)
!!$print*,'C   dx/dx x',i,l,l0,cagr%dc(1,i)%xyz_grad(:,l,l0)
!!$print*,'C   dx/dx y',i,l,l0,cagr%dc(2,i)%xyz_grad(:,l,l0)
!!$print*,'C   dx/dx z',i,l,l0,cagr%dc(3,i)%xyz_grad(:,l,l0)
                  enddo
               enddo
            endif

            if(integralpar_2dervs)then
               mov_atom1: do l=1,N_moving_unique_atoms+N_pc
                  if(l <= N_moving_unique_atoms) then
                     na=moving_unique_atom_index(l)
                     ea=unique_atoms(na)%n_equal_atoms
                  else
                     na=l-N_moving_unique_atoms
                     ea=pointcharge_array(na)%N_equal_charges
                  end if
                  uniq_atom1: do l0=1,ea
                     dr1_dai(:)=cagr%dR(j)%xyz_grad(:,l,l0); dr2_dai(:)=cagr%dR(k)%xyz_grad(:,l,l0)
                     do l1=1,3
                        da1_dai(l1,:)=cagr%dc(l1,j)%xyz_grad(:,l,l0)
                        da2_dai(l1,:)=cagr%dc(l1,k)%xyz_grad(:,l,l0)

                     end do

                     mov_atom2: do m=1,N_moving_unique_atoms+N_pc
                        if(m <= N_moving_unique_atoms) then
                           na1=moving_unique_atom_index(m)
                           ea1=unique_atoms(na1)%n_equal_atoms
                        else
                           na1=m-N_moving_unique_atoms
                           ea1=pointcharge_array(na1)%N_equal_charges
                        end if
                        uniq_atom2: do m0=1,ea1
                           dr1_daj(:)=cagr%dR(j)%xyz_grad(:,m,m0); dr2_daj(:)=cagr%dR(k)%xyz_grad(:,m,m0)
                           d2r1_daidaj(:,:)=cagr%dR(j)%xyz_hess(:,l,l0,:,m,m0)
                           d2r2_daidaj(:,:)=cagr%dR(k)%xyz_hess(:,l,l0,:,m,m0)
                           do m1=1,3
                              da1_daj(m1,:)=cagr%dc(m1,j)%xyz_grad(:,m,m0)
                              da2_daj(m1,:)=cagr%dc(m1,k)%xyz_grad(:,m,m0)
                              d2a1_daidaj(m1,:,:)=cagr%dc(m1,j)%xyz_hess(:,l,l0,:,m,m0)
                              d2a2_daidaj(m1,:,:)=cagr%dc(m1,k)%xyz_hess(:,l,l0,:,m,m0)
                           end do

                           d2r_dr1daj(:)=d2r_dr1dr1*dr1_daj(:)+d2r_dr1dr2*dr2_daj(:)+ &
                                matmul(d2r_dr1da1(:),da1_daj(:,:))+matmul(d2r_dr1da2(:),da2_daj(:,:))
                           d2r_dr2daj(:)=d2r_dr2dr1*dr1_daj(:)+d2r_dr2dr2*dr2_daj(:)+ &
                                matmul(d2r_dr2da1(:),da1_daj(:,:))+matmul(d2r_dr2da2(:),da2_daj(:,:))
                           do m2=1,3
                              d2r_da1daj(m2,:)=d2r_da1dr1(m2)*dr1_daj(:)+d2r_da1dr2(m2)*dr2_daj(:)+ &
                                   matmul(d2r_da1da1(m2,:),da1_daj(:,:))+matmul(d2r_da1da2(m2,:),da2_daj(:,:))
                              d2r_da2daj(m2,:)=d2r_da2dr1(m2)*dr1_daj(:)+d2r_da2dr2(m2)*dr2_daj(:)+ &
                                   matmul(d2r_da2da1(m2,:),da1_daj(:,:))+matmul(d2r_da2da2(m2,:),da2_daj(:,:))
                              d2a_dr1daj(m2,:)=d2a_dr1dr1(m2)*dr1_daj(:)+d2a_dr1dr2(m2)*dr2_daj(:)+ &
                                   matmul(d2a_dr1da1(m2,:),da1_daj(:,:))+matmul(d2a_dr1da2(m2,:),da2_daj(:,:))
                              d2a_dr2daj(m2,:)=d2a_dr2dr1(m2)*dr1_daj(:)+d2a_dr2dr2(m2)*dr2_daj(:)+ &
                                   matmul(d2a_dr2da1(m2,:),da1_daj(:,:))+matmul(d2a_dr2da2(m2,:),da2_daj(:,:))
                              do m3=1,3
                                 d2a_da1daj(m3,m2,:)=d2a_da1dr1(m3,m2)*dr1_daj(:)+d2a_da1dr2(m3,m2)*dr2_daj(:)+ &
                                      matmul(d2a_da1da1(m3,m2,:),da1_daj(:,:))+matmul(d2a_da1da2(m3,m2,:),da2_daj(:,:))
                                 d2a_da2daj(m2,m3,:)=d2a_da2dr1(m3,m2)*dr1_daj(:)+d2a_da2dr2(m3,m2)*dr2_daj(:)+ &
                                      matmul(d2a_da2da1(m3,m2,:),da1_daj(:,:))+matmul(d2a_da2da2(m3,m2,:),da2_daj(:,:))
                              end do
                           end do

                           do m1=1,3
                              cagr%dR(i)%xyz_hess(:,l,l0,m1,m,m0)= &
                                   d2r_dr1daj(m1)*dr1_dai(:)+drdr1*d2r1_daidaj(:,m1)+ &
                                   d2r_dr2daj(m1)*dr2_dai(:)+drdr2*d2r2_daidaj(:,m1)+ &
                                   matmul(d2r_da1daj(:,m1),da1_dai(:,:))+matmul(drda(:),d2a1_daidaj(:,:,m1))+ &
                                   matmul(d2r_da2daj(:,m1),da2_dai(:,:))-matmul(drda(:),d2a2_daidaj(:,:,m1))
!!$if(i==3) print*,'ABC d2r/dx2  ',i,l,l0,m,m0,cagr%dR(i)%xyz_hess(:,l,l0,m1,m,m0)
                              do m2=1,3
                                 cagr%dc(m2,i)%xyz_hess(:,l,l0,m1,m,m0)= &
                                      d2a_dr1daj(m2,m1)*dr1_dai(:)+dadr(m2)*d2r1_daidaj(:,m1)+ &
                                      d2a_dr2daj(m2,m1)*dr2_dai(:)-dadr(m2)*d2r2_daidaj(:,m1)+ &
                                      matmul(d2a_da1daj(m2,:,m1),da1_dai(:,:))+matmul(dada1(m2,:),d2a1_daidaj(:,:,m1))+ &
                                      matmul(d2a_da2daj(m2,:,m1),da2_dai(:,:))-matmul(dada2(m2,:),d2a2_daidaj(:,:,m1))
                              end do

!!$if(i==3) print*,'Atm d2x/dx2 x',i,l,l0,m,m0,cagr%dc(1,i)%xyz_hess(:,l,l0,m1,m,m0)
!!$if(i==3) print*,'Atm d2x/dx2 y',i,l,l0,m,m0,cagr%dc(2,i)%xyz_hess(:,l,l0,m1,m,m0)
!!$if(i==3) print*,'Atm d2x/dx2 z',i,l,l0,m,m0,cagr%dc(3,i)%xyz_hess(:,l,l0,m1,m,m0)
                           end do
                        end do uniq_atom2
                     end do mov_atom2
                  end do uniq_atom1
               end do mov_atom1
            end if
         enddo
      endif cav_f
!!$stop

    end subroutine calc_cavity_gradient
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine alloc_geom_deriv_part1(ca)
      !allocation of dR,dc gradients of geom_deriv ca
      type(geom_deriv) :: ca
      !** End of interface *****************************************

      integer(kind=i4_kind) :: i,alloc_stat,j
      integer :: N_pc

      if(VTN) return

      N_pc=0
!!$      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
      if(with_pc) N_pc=pointcharge_N

      allocate(ca%dR(N_spheres),stat=alloc_stat)
      if ( alloc_stat /= 0) call error_handler( &
           "alloc_geom_deriv_part: allocation of grad_atomic_center type  array is failed (dR)")
      allocate(ca%dc(3,N_spheres),stat=alloc_stat)
      if ( alloc_stat /= 0) call error_handler( &
           "alloc_geom_deriv_part: allocation of grad_atomic_center type  array is failed (dc)")

      do i=1,N_spheres
         allocate(ca%dR(i)%xyz_grad(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
         if ( alloc_stat /= 0) call error_handler( &
              "alloc_geom_deriv_part: allocation of xyz_grad is failed (dR)")
         ca%dR(i)%xyz_grad(:,:,:)=0.0_r8_kind
         if(integralpar_2dervs) then
            allocate(ca%dR(i)%xyz_hess(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
                                       3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
            if ( alloc_stat /= 0) call error_handler( &
                 "alloc_geom_deriv_part: allocation of xyz_hess is failed (dR)")
            ca%dR(i)%xyz_hess=0.0_r8_kind
         end if
         do j=1,3
            allocate(ca%dc(j,i)%xyz_grad(3,N_moving_unique_atoms+N_pc,ua_dim_max),&
                 stat=alloc_stat)
            if ( alloc_stat /= 0) call error_handler( "alloc_geom_&
                 &deriv_part: allocation of xyz_grad is failed (dc)")
            ca%dc(j,i)%xyz_grad(:,:,:)=0.0_r8_kind
            if(integralpar_2dervs) then
               allocate(ca%dc(j,i)%xyz_hess(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
                                       3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
               if ( alloc_stat /= 0) call error_handler( &
                    "alloc_geom_deriv_part: allocation of xyz_hess is failed (dc)")
               ca%dc(j,i)%xyz_hess=0.0_r8_kind
         end if
         enddo
      enddo

    end subroutine alloc_geom_deriv_part1
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine alloc_geom_deriv_part2(ca,n_max)
      !allocation vo dare,dcenter of geom_deriv ca
      !** End of interface *****************************************
      type(geom_deriv) :: ca
      integer(kind=i4_kind) ::n_max

      integer(kind=i4_kind) :: i,l,alloc_stat,j !,na,ea
      integer(i4_kind) :: i1,l1,j1
      integer :: N_pc

      if(VTN) return

      N_pc=0
!!$      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
      if(with_pc) N_pc=pointcharge_N

      allocate(ca%darea(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      allocate(ca%dcenter(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      if(integralpar_2dervs) then
         allocate(ca%d2area(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
              3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
         allocate(ca%d2center(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
              3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
      end if

      do i=1,N_moving_unique_atoms+N_pc
         do l=1,ua_dim_max
            do j=1,3
               allocate(ca%darea(j,i,l)%m(n_max),stat=alloc_stat)
               ASSERT(alloc_stat.eq.0)
               allocate(ca%dcenter(j,i,l)%m(3,n_max),stat=alloc_stat)
               ASSERT(alloc_stat.eq.0)
               ca%darea(j,i,l)%m=0.0_r8_kind
               ca%dcenter(j,i,l)%m=0.0_r8_kind
               if(integralpar_2dervs)then
                  do i1=1,N_moving_unique_atoms+N_pc
                     do l1=1,ua_dim_max
                        do j1=1,3
                           allocate(ca%d2area(j,i,l,j1,i1,l1)%m(n_max),stat=alloc_stat)
                           ASSERT(alloc_stat.eq.0)
                           allocate(ca%d2center(j,i,l,j1,i1,l1)%m(3,n_max),stat=alloc_stat)
                           ASSERT(alloc_stat.eq.0)
                           ca%d2area(j,i,l,j1,i1,l1)%m=0.0_r8_kind
                           ca%d2center(j,i,l,j1,i1,l1)%m=0.0_r8_kind
                        end do
                     end do
                  end do
               end if
            enddo
         enddo
      enddo

    end subroutine alloc_geom_deriv_part2
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine find_max_dim_surf_elem
      !find maximal numbers of centers needed for a subdivision
      !with average areas less than max_tes_area
      !** End of interface *****************************************
      integer (kind=i4_kind) :: i,l
      real (kind=r8_kind) :: area_comp

      max_cent=0
      do i=1,N_spheres
         area_comp=4.0_r8_kind*pi*r_sphere(i)**2
!        l=int(area_comp/m_t_a,kind=i4_kind)
         l=ceiling(area_comp/m_t_a)
         if(l > max_cent) max_cent=l
      enddo
DPRINT 'max cent',max_cent
    end subroutine find_max_dim_surf_elem
    !------------------------------------------------------------
!!!MF <<<<

    !------------------------------------------------------------
    subroutine tesselation
      !creates points on each sphere
      !calls def_poligon_tesselation ("cutting routine")
      !stores resulting surface tesselation in data_tes
      !** End of interface *****************************************

      type(triangles), allocatable :: triangles_on_sphere(:)
      type(poligon) :: tes_temp
      integer(i4_kind), allocatable :: int_sec_spheres(:)
      real(kind=r8_kind) :: xyz_pol_c(3),area_pol,area_tot,normal(3),volume_tot
      integer(kind=i4_kind) :: status,i,j,k,l,i1,j1,k1,l1
      integer(kind=i4_kind), allocatable :: i_sorted_spheres(:),N_cent_loc(:)
      integer(kind=i4_kind) :: ihelp,N_new_cent,N_new_xyz,n_int_sec_spheres
      real(kind=r8_kind) :: area_comp,xyz_buf(3)
      logical :: not_cav_dis_rep

      not_cav_dis_rep=(.not. do_cavitation .and. .not. do_disp_rep)

      allocate(i_sorted_spheres(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of i_sorted_spheres is failed")
      allocate(N_cent_loc(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation of N_cent_loc is failed")

      i_sorted_spheres(1)=1
      do i=2,N_spheres
         i_sorted_spheres(i)=i
         labj: do j=1,i-1
            if(r_sphere(i)<r_sphere(i_sorted_spheres(j))) then
               do k=i-1,j,-1
                  i_sorted_spheres(k+1)=i_sorted_spheres(k)
               enddo
               i_sorted_spheres(j)=i
               exit  labj
            endif
         enddo labj
      enddo

      allocate(triangles_on_sphere(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "tesselation: allocation of triangles_on_sphere is failed")

      !projecting set of points of arbitrary sphere on each sphere forming the cavity
      labisorted: do ihelp=1,N_spheres
         i=i_sorted_spheres(ihelp)

         if(.not. do_cavitation) then
!!$            if(do_disp_rep .or. .not.with_pc) then
#ifdef WITH_EFP
            if((.not.with_pc .or. (i_poly /= 1)) .and. .not. efp) then
#else
            if(.not.with_pc .or. (i_poly /= 1)) then
#endif
            area_comp=4.0_r8_kind*pi*r_sphere(i)**2
            local_point_factor=ceiling(log(&
                 area_comp/(N_centers_on_sphere*m_t_a))/log(4.0_r8_kind))
            do k=1,local_point_factor
               N_new_xyz=N_points_of_triangles+N_centers_on_sphere*3/2
               N_new_cent=N_centers_on_sphere*4
#ifdef _LINUX1
               call more_triangles(N_new_cent,N_points_of_triangles,&
                        N_centers_on_sphere)
#else
               call more_triangles(surf_elem%xyz,surf_elem%index,&
                    N_new_cent,N_points_of_triangles,N_centers_on_sphere)
#endif
               N_centers_on_sphere=N_new_cent
               N_points_of_triangles=N_new_xyz
            end do
            end if
            if(output_cavity_data .and. .not.do_gradients) &
                 write (output_unit,'(a27,i3,a11,i5)') &
                 'the number of triangles on ',i,' sphere is ', N_centers_on_sphere
         endif

         N_cent_loc(i)=N_centers_on_sphere
         do k=1,N_centers_on_sphere
            i1=surf_elem%index(k,1)
            l1=surf_elem%index(k,2)
            k1=surf_elem%index(k,3)
            xyz_buf=(surf_elem%xyz(i1,:)+surf_elem%xyz(l1,:)+surf_elem%xyz(k1,:))/3.0_r8_kind
            surf_elem%xyz_centers(k,:)=(radius/LEN3(xyz_buf))*xyz_buf
         enddo
         allocate(triangles_on_sphere(i)%xyz(N_points_of_triangles,3), &
              triangles_on_sphere(i)%index(N_centers_on_sphere,3), &
              triangles_on_sphere(i)%xyz_centers(N_centers_on_sphere,3), &
              stat=status)
         if ( status /= 0) then
            call error_handler( &
                 "tesselation: allocation of triangles_on_sphere(1) is failed")
         endif
         do j=1,N_points_of_triangles
            triangles_on_sphere(i)%xyz(j,:)= surf_elem%xyz(j,:)*r_sphere(i)/radius+ &
                 xyz_sphere(i,:)
         enddo
         triangles_on_sphere(i)%index(1:N_centers_on_sphere,:)= &
              surf_elem%index(1:N_centers_on_sphere,:)
         do k=1,N_centers_on_sphere
            triangles_on_sphere(i)%xyz_centers(k,:)= surf_elem%xyz_centers(k,:)*r_sphere(i)/radius+ &
                 xyz_sphere(i,:)
         enddo
         triangles_on_sphere(i)%radius=r_sphere(i)
         if(smooth /= 0 .and. not_cav_dis_rep) then
            triangles_on_sphere(i)%area=(4.0_r8_kind*pi*r_sphere(i)**2)/N_centers_on_sphere
         else
            triangles_on_sphere(i)%area=0.0_r8_kind
         end if
!       if(i.eq.1) print*,'lab ed',r_sphere(i),triangles_on_sphere(1)%xyz_centers(3,:)
!       if(i.eq.1) print*,triangles_on_sphere(1)%index(3,:)
      enddo labisorted

!!! MF: why not?!
!     deallocate(i_sorted_spheres, stat=status)
      deallocate(i_sorted_spheres,surf_elem%xyz_centers,surf_elem%xyz,surf_elem%index, stat=status)
      if ( status /= 0) call error_handler( &
           "tesselation: deallocation of i_sorted_spheres is failed")

      allocate(xyz_tes_c(N_spheres*N_centers_on_sphere,3), &
           area_tes(N_spheres*N_centers_on_sphere), &
           r_tes(N_spheres*N_centers_on_sphere), &
           sphere(N_spheres*N_centers_on_sphere), &
           cuttt(N_spheres*N_centers_on_sphere), &
           data_tes(N_spheres*N_centers_on_sphere),stat=status)
      if ( status /= 0) call error_handler( &
           "tesselation: allocation of xyz_tes_c, area_tes, r_tes are failed")

!!! MF 5/2000 gradients
      if(do_gradients) then
        call alloc_geom_deriv_part2(cagr,N_spheres*N_centers_on_sphere)
      endif

      allocate(int_sec_spheres(N_spheres),stat=status)
      if ( status /= 0) call error_handler( &
           "tesselation: failed int_sec_spheres allocation")

      if(not_cav_dis_rep) then
         if(smooth == 1 .or. smooth == 2) call calc_smooth_function()
      end if

      !definition of the cavity elements
      area_tot=0.0_r8_kind
      volume_tot=0.0_r8_kind
      N_total=0
      lab1: do i=1,N_spheres
         if(.not. do_cavitation .and. .not. do_disp_rep) then
            if(gepol == 93) then
               if(zero_area(i)) cycle lab1
            end if
         end if

         call calc_intersect_spheres(i,xyz_sphere,r_sphere,N_spheres, &
              n_int_sec_spheres,int_sec_spheres)

         lab2: do j=1,N_cent_loc(i)
           j1=triangles_on_sphere(i)%index(j,1)
           k1=triangles_on_sphere(i)%index(j,2)
           l1=triangles_on_sphere(i)%index(j,3)

           if(smooth /= 0 .and. not_cav_dis_rep) then
!!$              if(center_inside_cav(triangles_on_sphere(i)%xyz_centers(j,:), &
!!$                   xyz_sphere,r_sphere,N_spheres, &
!!$                   n_int_sec_spheres,int_sec_spheres)) cycle lab2
           else
              if(if_tess_inside_cav(triangles_on_sphere(i)%xyz(j1,:), &
                   triangles_on_sphere(i)%xyz(k1,:), &
                   triangles_on_sphere(i)%xyz(l1,:), &
                   i,N_cent_loc(i),xyz_sphere,r_sphere,N_spheres, &
                   n_int_sec_spheres,int_sec_spheres)) cycle lab2
           end if

!       if(i.eq.1.and.j.eq.3) print*
!       if(i.eq.1.and.j.eq.3) print*,triangles_on_sphere(i)%xyz(j1,:)
!       if(i.eq.1.and.j.eq.3) print*,triangles_on_sphere(i)%xyz(k1,:)
!       if(i.eq.1.and.j.eq.3) print*,triangles_on_sphere(i)%xyz(l1,:)
!
!print*,'sphere- ',i,' tess- ',N_total
           if((smooth == 1 .or. smooth == 2) .and. not_cav_dis_rep) then
              xyz_pol_c=triangles_on_sphere(i)%xyz_centers(j,:)
              area_pol=triangles_on_sphere(i)%area
              call smooth_tess_area(i,j,xyz_pol_c,area_pol,tes_temp)
           else if(smooth == 3 .and. not_cav_dis_rep) then
              xyz_pol_c=triangles_on_sphere(i)%xyz_centers(j,:)
              area_pol=triangles_on_sphere(i)%area
              call smooth_tess_area1(i,j,xyz_pol_c,area_pol,tes_temp)
           else if(smooth == 0 .or. .not. not_cav_dis_rep) then
              call def_poligon_tes(i, &
                   triangles_on_sphere(i)%xyz(j1,:), &
                   triangles_on_sphere(i)%xyz(k1,:), &
                   triangles_on_sphere(i)%xyz(l1,:), &
                   xyz_pol_c,area_pol,tes_temp)
           end if

           if(area_pol /= 0.0_r8_kind) then
!          if(N_total.eq.2) print*,'def pol tes',j1,k1,l1,i,j
                  N_total=N_total+1
                  xyz_tes_c(N_total,:)=xyz_pol_c
                  area_tes(N_total)=area_pol
                  area_tot=area_tot+area_tes(N_total)
                  r_tes(N_total)=triangles_on_sphere(i)%radius
                  sphere(N_total)=i
                  normal(:)=(xyz_pol_c(:)-xyz_sphere(i,:))/r_sphere(i)
                  volume_tot=volume_tot+dot_product(normal,xyz_pol_c)*area_pol/3.0_r8_kind

                  cuttt(N_total)=.false.
                  do l=1,tes_temp%n_vertises
                        if(tes_temp%n_sphere(l,1)/=0) cuttt(N_total) = .true.
                  enddo

                  data_tes(N_total)=tes_temp
           endif
         enddo lab2
      enddo lab1

      if(AreaScaling) call scale_tess_area()

      if(output_cavity_data .and. .not.do_gradients) then
         write(output_unit,*) '-------------------------------------------------------------'
         write(output_unit,'(i5,a27)') N_total,' points have been generated'
         write(output_unit,'(f18.9,a27)') area_tot,' a_u^2 total surface area  '
         write(output_unit,'(f18.9,a27)') volume_tot,' a_u^3 total cavity volume '
         if(output_cavity_long) then
           write(output_unit,*)'number                   coordinates(a.u.)               area   n_vert sphere'
           do i=1,N_total
!!$              write(output_unit,*) data_tes(i)%sphere
            write(output_unit,'(1x,i4,3x,3(1x,f13.9),1x,f14.9,1x,i3,1x,i3)') i,xyz_tes_c(i,:), &
                 area_tes(i),data_tes(i)%n_vertises,data_tes(i)%sphere
!!$            do j=1,data_tes(i)%n_vertises
!!$               write(output_unit,*) data_tes(i)%xyz_vertex(j,:)
!!$            end do
           enddo
         endif
         write(output_unit,*) '-------------------------------------------------------------'
      endif
!       print*, 'ter ter' ,2, xyz_tes_c(2,:)

      do i=1,N_spheres
         deallocate(triangles_on_sphere(i)%xyz, &
              triangles_on_sphere(i)%index, &
              triangles_on_sphere(i)%xyz_centers, &
              stat=status)
         if ( status /= 0) then
            call error_handler( &
                 "tesselation: deallocation of triangles_on_sphere(1) is failed")
         endif

      enddo
      deallocate(triangles_on_sphere,stat=status)
      if ( status /= 0) call error_handler( &
         "tesselation: deallocation of triangles_on_sphere(2) is failed")

      deallocate(N_cent_loc,cut_rad_sort, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation of N_cent_loc is failed")

      deallocate(int_sec_spheres,stat=status)
      if ( status /= 0) call error_handler( &
           "tesselation: failed int_sec_spheres deallocation")

      if(do_gradients .and. .not. do_disp_rep) then
         if(.not. VTN) call dealloc_geom_deriv_part1(cagr)
      endif

    end subroutine tesselation
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine smooth_tess_area(current_sphere,current_tess,xyz_center,area,poli)
      !  Purpose: TsLess smooth partitioning cavity surface
      !  C. S.Pomelli J Comput Chem 25, 1532, 2004
      !------------ Modules used ------------------- ---------------
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind), intent(in) :: current_sphere,current_tess
      real(r8_kind), intent(in)    :: xyz_center(3)
      real(r8_kind), intent(inout) :: area
      type(poligon), intent(out)   :: poli
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      integer(i4_kind)         :: i,j
      real(r8_kind)            :: rho,area0,x,sigma,sigma_tot
      logical                  :: reduced
      !------------ Executable code --------------------------------

      reduced=.false.
      area0=area
      Rho=sqrt(area0/pi)
      if(smooth == 2) Rho=Rho/4.0_r8_kind

      poli%n_vertises=0

      sigma_tot=1.0_r8_kind
      a1:do i=1,N_spheres
!         if(zero_area(i)) cycle a1
         if(i==current_sphere) cycle a1

         x=DOT3(xyz_center(:)-xyz_sphere(i,:), xyz_center(:)-xyz_sphere(i,:))

         x=(sqrt(x)-r_sphere(i))/rho

         if(x <= d) then
            area=0.0_r8_kind
            return
         elseif (x > d .and. x < 1.0_r8_kind+d) then
            sigma=0.0_r8_kind
            do j=1,Nd
               sigma=sigma+a(j)*x**(j-1)
            end do
         else
            sigma=1.0_r8_kind
         end if
         sigma_tot=sigma_tot*sigma
      end do a1

      area=area0*sigma_tot
      if(sigma_tot /= 1.0_r8_kind) then
         reduced=.true.
         poli%n_vertises=1
         poli%n_sphere(1,1)=1
      end if
      poli%sphere=current_sphere

      if(do_gradients) then
         if(reduced) then
            call geom_grad_reduced_triang(current_sphere,area0,&
                                xyz_center,Rho)
         else
            call geom_grad_orig_triang(current_sphere,area0,&
                                xyz_center)
         endif
      endif

    end subroutine smooth_tess_area
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine calc_smooth_function()
      !  Purpose: Calulation of switching function for TsLess
      !------------ Modules used ------------------- ---------------
      implicit none
      !------------ Declaration of formal parameters ---------------
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      real(r8_kind)    :: mat(7,7),vec(7,1),ipiv(7)
      real(r8_kind)    :: mat1(6,6),vec1(6,1),ipiv1(6)
      integer(i4_kind) :: i,info !,status,j,k
      !------------ Executable code --------------------------------

      if(smooth == 1) then
         Nd=7
         do i=1,7
            mat(1,i)=d**(i-1)
            mat(2,i)=1.0_r8_kind
            mat(3,i)=(i-1)*d**(i-2)
            mat(4,i)=(i-1)*1.0_r8_kind
            mat(5,i)=(i-1)*(i-2)*d**(i-3)
            mat(6,i)=(i-1)*(i-2)*1.0_r8_kind
            mat(7,i)=(1-d**i)/i
         end do

         vec(1,1)=0.0_r8_kind
         vec(2,1)=1.0_r8_kind
         vec(3:6,1)=0.0_r8_kind
         vec(7,1)=1.0_r8_kind-d

         call dgesv(7,1,mat,7,ipiv,vec,7,info)
         ASSERT(info==0)

         a(1:7)=vec(1:7,1)
      else if(smooth == 2) then
         Nd=6
         do i=1,6
            mat1(1,i)=d**(i-1)
            mat1(2,i)=1.0_r8_kind
            mat1(3,i)=(i-1)*d**(i-2)
            mat1(4,i)=(i-1)*1.0_r8_kind
            mat1(5,i)=(i-1)*(i-2)*d**(i-3)
            mat1(6,i)=(i-1)*(i-2)*1.0_r8_kind
         end do

         vec1(1,1)=0.0_r8_kind
         vec1(2,1)=1.0_r8_kind
         vec1(3:6,1)=0.0_r8_kind

         call dgesv(6,1,mat1,6,ipiv1,vec1,6,info)
         ASSERT(info==0)

         a(1:6)=vec1(1:6,1)
      end if

    end subroutine calc_smooth_function
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine smooth_tess_area1(current_sphere,current_tess,xyz_center,area,poli)
      !  Purpose: FIXPVA smooth partitioning cavity surface
      !  P. Su, H. Li  J Chem Phys 130 074109 2009
      !------------ Modules used ------------------- ---------------
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind), intent(in) :: current_sphere,current_tess
      real(r8_kind), intent(in)    :: xyz_center(3)
      real(r8_kind), intent(inout) :: area
      type(poligon), intent(out)   :: poli
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      real(r8_kind)            :: f1,f2,nn,mm !n,m,
      integer(i4_kind)         :: i
      real(r8_kind)            :: r1(3),r2(3),r3(3),r5(3) !,r4(3)
      real(r8_kind)            :: r2r3(3),dr2r3,r1r5(3),r1r3(3),dr1r3
      real(r8_kind)            :: rc,cosa,sina,cosb,sinb,acb,ac,acd
      real(r8_kind)            :: n1n1,n2n2,m1m1,m2m2
      real(r8_kind)            :: sigma,sigma_tot,area0
      logical                  :: reduced
      !------------ Executable code --------------------------------

      reduced=.false.
      area0=area
      poli%n_vertises=0

      n1n1=nn11*nn11; n2n2=nn12*nn12; m1m1=mm11*mm11; m2m2=mm12*mm12

      r2=xyz_sphere(current_sphere,:)
      r1=xyz_center
      rc=r_sphere(current_sphere)

      sigma_tot=1.0_r8_kind
      a1:do i=1,N_spheres
!         if(zero_area(i)) cycle a1
         if(i==current_sphere) cycle a1

         r3=xyz_sphere(i,:)
         r2r3=r2-r3
         dr2r3=LEN3(r2r3)
         r1r3=r1-r3
         dr1r3=LEN3(r1r3)

         acb=rc*rc+dr2r3*dr2r3-r_sphere(i)*r_sphere(i)
         ac=2.0_r8_kind*rc*dr2r3
         cosa=acb/ac
         sina=sqrt(ac*ac-acb*acb)/ac

         acd=rc*rc+dr2r3*dr2r3-dr1r3*dr1r3
         cosb=acd/ac
         sinb=sqrt(ac*ac-acd*acd)/ac

         mm=2.0_r8_kind*rc*rc*(1.0_r8_kind-cosb*cosa-sinb*sina)

         if(dr2r3 >= rc+r_sphere(i)) then
            f1=1.0_r8_kind
         else if(dr1r3 <= r_sphere(i)) then
            f1=0.0_r8_kind
         else if(mm > m2m2) then
            f1=1.0_r8_kind
         else if(mm <= m1m1) then
            f1=0.0_r8_kind
         else if(mm > m1m1 .and. mm<= m2m2) then
            f1=10.0_r8_kind*((mm-m1m1)/(m2m2-m1m1))**3- &
               15.0_r8_kind*((mm-m1m1)/(m2m2-m1m1))**4+ &
                6.0_r8_kind*((mm-m1m1)/(m2m2-m1m1))**5
         end if

         r5=r2r3*r_sphere(i)/dr2r3+r3
         r1r5=r1-r5
         nn=DOT3(r1r5,r1r5)

         if(nn > n2n2) then
            f2=1.0_r8_kind
         else if(nn <= n1n1) then
            f2=0.0_r8_kind
         else if(nn > n1n1 .and. nn<= n2n2) then
            f2=10.0_r8_kind*((nn-n1n1)/(n2n2-n1n1))**3- &
               15.0_r8_kind*((nn-n1n1)/(n2n2-n1n1))**4+ &
                6.0_r8_kind*((nn-n1n1)/(n2n2-n1n1))**5
         end if

         if(f1==0.0_r8_kind .or. f2==0.0_r8_kind) then
            area=0.0_r8_kind
            return
         end if

         sigma=f1*f2

         sigma_tot=sigma_tot*sigma
      end do a1

      area=area0*sigma_tot
      if(sigma_tot /= 1.0_r8_kind) then
         reduced=.true.
         poli%n_vertises=1
         poli%n_sphere(1,1)=1
      end if
      poli%sphere=current_sphere

      if(do_gradients) then
         if(reduced) then
            call geom_grad_reduced_triang1(current_sphere,area0,&
                                xyz_center)
         else
            call geom_grad_orig_triang(current_sphere,area0,&
                                xyz_center)
         endif
      endif

    end subroutine smooth_tess_area1
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine scale_tess_area()
      !HUI LI, JAN H. JENSEN, J Comput Chem 25: 1449-1462, 2004
      !used in combination with VTN approach and for EFP calculations

      integer(i4_kind) :: i,j,k,l
      real(r8_kind) :: P1,P2,r1(3),r2(3),norm1(3),norm2(3)
      real(r8_kind) :: dr12,area1,area2,scale_factor

      P1=1.0_r8_kind/sqrt(48.0_r8_kind)

      do i=1,N_total
         k=sphere(i)
         r1=xyz_tes_c(i,:)-xyz_sphere(k,:)
         norm1=r1/LEN3(r1)
         area1=area_tes(i); area1=sqrt(area1)

         do j=i,N_total
            l=sphere(j)
            if(k==l) cycle

            r2=xyz_tes_c(j,:)-xyz_sphere(l,:)
            norm2=r2/LEN3(r2)
            area2=area_tes(j); area2=sqrt(area2)

            P2=1.0_r8_kind+(1.0_r8_kind-DOT3(norm1,norm2))**4

            dr12=LEN3(xyz_tes_c(i,:)-xyz_tes_c(j,:))

            if(dr12 < P1*P2*(area1+area2)) then
               scale_factor=(dr12/(P1*P2*(area1+area2)))**2
               area_tes(i)=area_tes(i)*scale_factor
               area_tes(j)=area_tes(j)*scale_factor
            end if
         end do
      end do

    end subroutine scale_tess_area
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine init_cut_rad_sort()
      integer(kind=i4_kind) :: i,j,status,k,l
      real (kind=r8_kind),allocatable :: rcut(:)
      real (kind=r8_kind) :: help,vd(3)

      allocate(cut_rad_sort(N_spheres,N_spheres),rcut(N_spheres),stat=status)
      if(status/=0) call error_handler("alloc of cut_rad_sort failed")
      cut_rad_sort=0

      do i=1,N_spheres
!       print*,xyz_sphere(i,:),i,'xyz_sphere'
         rcut=-1.0_r8_kind
         labj : do j=1,N_spheres
            if(i==j) cycle labj
            vd(:)=xyz_sphere(j,:)-xyz_sphere(i,:)
            help=LEN3(vd)
            if (r_sphere(i)+r_sphere(j) < help) cycle
            help=(help**2+r_sphere(i)**2-r_sphere(j)**2)/&
                                (2.0_r8_kind*help)
            help=r_sphere(i)**2-help**2
            if(help<0.0_r8_kind) help=0.0_r8_kind
            rcut(j)=sqrt(help)
         enddo labj

         do j=1,N_spheres
            cut_rad_sort(i,j)=j
            labk: do k=1,j-1
               if(rcut(cut_rad_sort(i,j))>rcut(cut_rad_sort(i,k))) then
                  do l=j-1,k,-1
                     cut_rad_sort(i,l+1)=cut_rad_sort(i,l)
                  enddo
                  cut_rad_sort(i,k)=j
                  exit  labk
               endif
            enddo labk
         enddo
      enddo
!       print*,sum(abs(rcut)),'rcut sum abs',sum(abs(xyz_sphere(1:N_spheres,:))),sum(r_sphere(1:N_spheres))

      deallocate(rcut)
    end subroutine init_cut_rad_sort
    !------------------------------------------------------------
    !------------------------------------------------------------
    subroutine def_poligon_tes(sphere_number, &
         xyz_triangl_1, xyz_triangl_2, xyz_triangl_3, &
         xyz_poligon_centre, poligon_area, poli)
      ! This is very important procedure calculating area and coordinates
      ! of geometrical center of tessarea of polygonal shape.
      ! The polygonal type of tessarea is formed as result of intersection
      ! between two or some spheres.
      implicit none
      integer(kind=i4_kind), intent(in) :: sphere_number
      real(kind=r8_kind), intent(in) :: xyz_triangl_1(:), xyz_triangl_2(:), xyz_triangl_3(:)
      real(kind=r8_kind), intent(out) :: xyz_poligon_centre(:), poligon_area
      type(poligon), intent(out) :: poli
      !** End of interface *****************************************

      real(kind=r8_kind) :: new_vertex(6)
      real(kind=r8_kind) :: alpha, beta, cos_gamma,v1(3),v2(3),v_buf(3)
      logical :: combination(MAX_POL_VER,MAX_POL_VER)
      logical :: inside_flag(MAX_POL_VER),cutted
      real(kind=r8_kind) :: d_rad,d_rad1,d_vv
      real(kind=r8_kind), parameter :: small_d=1.0e-8_r8_kind

      integer(kind=i4_kind) :: ncut,k_new,l1,l2,n_extra
      integer(kind=i4_kind) :: i,i1,ii,j,k,jh,k1
      real(kind=r8_kind) :: tn(MAX_POL_VER,3),tns(MAX_POL_VER,3),vn(MAX_POL_VER,3),vns(MAX_POL_VER,3)
      logical :: extra_point(MAX_POL_VER),h_case(MAX_POL_VER,2),h_case_sphere
      real(kind=r8_kind), dimension(MAX_POL_VER) :: costhn,phin
      real(kind=r8_kind) :: cut_center(3),cut_radius,help,tanphi2,center_factor
      integer(kind=i4_kind) :: is(MAX_POL_VER),num_v(MAX_POL_VER),n_tes_part,max_tes_part

      character(len=54) :: message

      tn=0.0_r8_kind
      tns=0.0_r8_kind
      vn=0.0_r8_kind
      vns=0.0_r8_kind

      poli%xyz_vertex=0.0_r8_kind
      poli%bounds=0
      poli%r_bound=0.0_r8_kind
      poli%xyz_bound=0.0_r8_kind
      poli%n_sphere=0

      counter_n=counter_n+1

      poli%sphere=sphere_number
      poli%n_vertises=3

      poli%xyz_vertex(1,:)=xyz_triangl_1
      poli%xyz_vertex(2,:)=xyz_triangl_2
      poli%xyz_vertex(3,:)=xyz_triangl_3

      ! ensure anticlockwise sort
      v_buf=vector_product(xyz_triangl_1-xyz_sphere(sphere_number,:),&
        xyz_triangl_2-xyz_sphere(sphere_number,:))

      if(DOT3(v_buf,xyz_triangl_3-xyz_sphere(sphere_number,:))&
                                                <0.0_r8_kind) then
        poli%xyz_vertex(2,:)=xyz_triangl_3
        poli%xyz_vertex(3,:)=xyz_triangl_2
      endif

      poli%bounds(1,1)=2
      poli%bounds(1,2)=3
      poli%bounds(2,1)=3
      poli%bounds(2,2)=1
      poli%bounds(3,1)=1
      poli%bounds(3,2)=2

      poli%r_bound(1:3,:)=r_sphere(sphere_number)

      poli%xyz_bound(1,1:3)=xyz_sphere(sphere_number,:)
      poli%xyz_bound(2,1:3)=xyz_sphere(sphere_number,:)
      poli%xyz_bound(3,1:3)=xyz_sphere(sphere_number,:)
      poli%xyz_bound(1,4:6)=xyz_sphere(sphere_number,:)
      poli%xyz_bound(2,4:6)=xyz_sphere(sphere_number,:)
      poli%xyz_bound(3,4:6)=xyz_sphere(sphere_number,:)

      cutted=.false.
      n_extra=0_i4_kind
      extra_point=.false.
      h_case=.false.

      j_spheres: do jh=1,N_spheres
         j=cut_rad_sort(sphere_number,jh)
         inside_flag(:) = .false.
         h_case_sphere=.false.
         if (j==sphere_number) cycle j_spheres
         if(.not. do_cavitation .and. .not. do_disp_rep) then
            if(gepol == 93) then
               if(zero_area(j)) cycle j_spheres
            end if
         end if
         v1=xyz_sphere(j,:)-xyz_sphere(sphere_number,:)
         d_rad1=LEN3(v1)
         if (r_sphere(sphere_number)+r_sphere(j) < d_rad1) cycle j_spheres

         d_vv=(d_rad1**2+r_sphere(sphere_number)**2-r_sphere(j)**2)/&
                                (2.0_r8_kind*d_rad1)

         cut_center=xyz_sphere(sphere_number,:)+v1/d_rad1*d_vv
         help=r_sphere(sphere_number)**2-d_vv**2
         if(help<0.0_r8_kind) help=0.0_r8_kind
         cut_radius=sqrt(help)

         !h_case: sphere more than half in an other, typically
         !for H-Atoms
         if(DOT3(v1,cut_center(:)-xyz_sphere(sphere_number,:))<0.0_r8_kind) then
            h_case_sphere=.true.
         endif

         ! definition of the vertises which are inside of J-th sphere
         k1=0
         do i=1,poli%n_vertises
            d_rad=LEN3(poli%xyz_vertex(i,:)-xyz_sphere(j,:))
            ! is vertex inside of sphere?
            if (d_rad < r_sphere(j) - small_d) then
               k1=k1+1
               inside_flag(i)= .true.
               cutted=.true.
            endif
         enddo
         ! stop tesselation because the tessarea is complitely within the J-th sphere
         if(k1==poli%n_vertises .and. d_rad1>=r_sphere(j)) then
            poligon_area=0.0_r8_kind
            return
         end if

         k=poli%n_vertises+1
         i_vertex: do i=1,poli%n_vertises
            i1=poli%bounds(i,1)

!!$print*,sphere_number,j,inside_flag(i),inside_flag(i1),'!!!!!!!!!!!!!'
            call line_cut(poli%xyz_vertex(i,:),poli%xyz_vertex(i1,:),&
                poli%xyz_bound(i,1:3),poli%r_bound(i,1),&
                xyz_sphere(j,:),r_sphere(j),&
                inside_flag(i),inside_flag(i1), &
                xyz_sphere(sphere_number,:),h_case(i,1),&
                ncut,new_vertex)
!!$print*,sphere_number,j,inside_flag(i),inside_flag(i1),ncut,'!!!!!!!!!!!!!'

            if(k+ncut>MAX_POL_VER) call error_handler("to many vert")

            select case(ncut)
              case (1)
                  poli%xyz_vertex(k,:)=new_vertex(1:3)
!       if(sphere_number.eq.2) then
!          print*,'case 1',poli%xyz_vertex(k,:),jh,i,i1,j
!          print*,sum(abs(poli%xyz_bound(i,1:3))), &
!                   poli%r_bound(i,1),r_sphere(j),h_case(i,1),inside_flag(i),inside_flag(i1)
!       endif

                  poli%bounds(k,1)=i1
                  poli%bounds(k,2)=i
                  poli%bounds(i1,2)=k
                  poli%bounds(i,1)=k
                  if(inside_flag(i) .and. .not. inside_flag(i1) ) then
                     l1=1
                     l2=2
                  else if(.not. inside_flag(i) .and. inside_flag(i1) ) then
                     l1=2
                     l2=1
                  else
                     write(message,'(a19,i3,a17,i3,a9,i3)' ) &
                          'tessarea on sphere ',sphere_number,', cutting sphere ',j,', vertex ',k
                     call error_handler("def_poligon_tes: this case (1) &
                                & should not occure"//achar(10)//trim(message))
                  endif

                  poli%r_bound(k,l1) = poli%r_bound(i,1)
                  poli%xyz_bound(k,3*l1-2:3*l1) = poli%xyz_bound(i,1:3)
                  poli%n_sphere(k,l1)=poli%n_sphere(i,1)
                  h_case(k,l1)=h_case(i,1)
                  poli%r_bound(k,l2) = cut_radius
                  poli%xyz_bound(k,3*l2-2:3*l2) = cut_center(:)
                  poli%n_sphere(k,l2)=j
                  if(h_case_sphere) h_case(k,l2)=.true.
                  k=k+1
                  cutted=.true.

              case (11)
                  ! add a point, which is not really a vertex but a point on an
                  ! edge
                  ! only, if no h_case edge
                  if(.not. h_case(i,1)) then
                    poli%xyz_vertex(k,:)=new_vertex(1:3)
!            if(sphere_number.eq.2) print*,'case 11',poli%xyz_vertex(k,:)

                    poli%bounds(k,1)=i1
                    poli%bounds(k,2)=i
                    poli%bounds(i1,2)=k
                    poli%bounds(i,1)=k

                    poli%r_bound(k,1) = poli%r_bound(i,1)
                    poli%r_bound(k,2) = poli%r_bound(i,1)
                    poli%xyz_bound(k,1:3) = poli%xyz_bound(i,1:3)
                    poli%xyz_bound(k,4:6) = poli%xyz_bound(i,1:3)
                    poli%n_sphere(k,1)=poli%n_sphere(i,1)
                    poli%n_sphere(k,2)=poli%n_sphere(i,1)
                    h_case(k,1)=h_case(i,1)
                    h_case(k,2)=h_case(i,1)
                    n_extra=n_extra+1
                    extra_point(k)=.true.
                    k=k+1
                    cutted=.true.
                  endif

              case (2)
!AS 28.12.2005
!Here we try to avoid situation when the j-th sphere does not cross any edges of the tessarea
!but only touchs one of them. Othewise we can have situation of creating pseudo vertex
!with angle exactly equal 180 degree
                  v_buf=new_vertex(1:3)-new_vertex(4:6)
                  if(LEN3(v_buf)<small_d .and. &
                       .not.inside_flag(i) .and. .not. inside_flag(i1)) cycle i_vertex

                  poli%xyz_vertex(k,:)=new_vertex(1:3)
                  poli%xyz_vertex(k+1,:)=new_vertex(4:6)
!       if(sphere_number.eq.2) print*,'case (2)',poli%xyz_vertex(k,:)

                  poli%bounds(i,1)=k
                  poli%bounds(k,1)=k+1
                  poli%bounds(k+1,1)=i1
                  poli%bounds(k,2)=i
                  poli%bounds(k+1,2)=k
                  poli%bounds(i1,2)=k+1
                  if(inside_flag(i) .and. inside_flag(i1)) then
                     l1=1
                     l2=2
                  else if(.not.inside_flag(i) .and. .not. inside_flag(i1)) then
                     l1=2
                     l2=1
                  else
                     write(message,'(a19,i3,a17,i3,a9,i3)' ) &
                          'tessarea on sphere ',sphere_number,', cutting sphere ',j,', vertex ',k
                     call error_handler("def_poligon_tes: this case (2) &
                                & should not occure"//achar(10)//trim(message))
                  endif

                  poli%r_bound(k,l1)=poli%r_bound(i,1)
                  poli%r_bound(k+1,l2)=poli%r_bound(i,1)
                  poli%r_bound(k,l2)=cut_radius
                  poli%r_bound(k+1,l1)=cut_radius
                  poli%n_sphere(k,l1)=poli%n_sphere(i,1)
                  h_case(k,l1)=h_case(i,1)
                  poli%n_sphere(k+1,l2)=poli%n_sphere(i,1)
                  h_case(k+1,l2)=h_case(i,1)
                  poli%n_sphere(k,l2)=j
                  poli%n_sphere(k+1,l1)=j
                  poli%xyz_bound(k,3*l1-2:3*l1)=poli%xyz_bound(i,1:3)
                  poli%xyz_bound(k+1,3*l2-2:3*l2)=poli%xyz_bound(i,1:3)
                  poli%xyz_bound(k,3*l2-2:3*l2)=cut_center(:)
                  poli%xyz_bound(k+1,3*l1-2:3*l1)=cut_center(:)
                  if(h_case_sphere) then
                        h_case(k,l2)=.true.
                        h_case(k+1,l1)=.true.
                  endif
                  k=k+2
                  cutted=.true.
                case default
                   if(inside_flag(i)  .neqv. inside_flag(i1)) then
                      write(message,'(a19,i3,a17,i3,a9,i3)' ) &
                           'tessarea on sphere ',sphere_number,', cutting sphere ',j,', vertex ',k
                      call error_handler("this case (0) &
                           & should not occure"//achar(10)//trim(message))
                   end if
                   ! do nothing in cases n_cut=0 n_cut=3
                   ! n_cut=0 appears also, if both vertises are inside.
                   ! in this case, cutted=.true. is already set above
            end select

         enddo i_vertex
         k_new=k-1-poli%n_vertises
         ! sort in and out new / old corners (only old corners, new are on
         ! cutting circle by construction
         ! delete very close vertises
         do i=1, poli%n_vertises
                if(.not. inside_flag(i)) then
                        lab_l1 : do l1=1,2
                           i1=poli%bounds(i,l1)
                           ii=1
                           if(l1==1) ii=2
                           v_buf=poli%xyz_vertex(i,:)-&
                                poli%xyz_vertex(i1,:)
                           if(LEN3(v_buf)<small_d) then
                                poli%bounds(i1,ii)=poli%bounds(i,ii)
                                poli%r_bound(i1,ii)=poli%r_bound(i,ii)
                                poli%xyz_bound(i1,ii*3-2:ii*3)=&
                                        poli%xyz_bound(i,3*ii-2:3*ii)
                                poli%n_sphere(i1,ii)=poli%n_sphere(i,ii)
                                h_case(i1,ii)=h_case(i,ii)
                                inside_flag(i)=.true.
                                exit lab_l1
                           endif
                        enddo lab_l1
                endif
         enddo

         do i=1, poli%n_vertises
                if(inside_flag(i)) then
                        poli%bounds(poli%bounds(i,1),2)=poli%bounds(i,2)
                        poli%bounds(poli%bounds(i,2),1)=poli%bounds(i,1)
                endif
          enddo
         num_v=0
         if(k_new-n_extra>3 .and. .not. h_case_sphere) then
           !re-connect vertises, in case, tessera is cutted in several pieces
           n_tes_part=1
           is=-1
           ii=0
           new_points: do k=poli%n_vertises+1,poli%n_vertises+k_new
            if(extra_point(k).or.inside_flag(k)) cycle new_points
            do i=1,ii
              if(is(i)==k) cycle new_points
            enddo
            !k belongs to an separate piece of tessera
            ii=ii+1
            is(ii)=k
            !search for new vertexes belonging to the same piece of tessera
            i1=k
            do
               i1=poli%bounds(i1,1)
               if (i1==k) exit
               if (.not. extra_point(i1) .and. i1>poli%n_vertises) then
                  ii=ii+1
                  is(ii)=i1
               endif
            enddo
            num_v(n_tes_part)=ii
            n_tes_part=n_tes_part+1
           enddo new_points
           max_tes_part=n_tes_part-1
           do k=max_tes_part,2,-1
            num_v(k)=num_v(k)-num_v(k-1)
           enddo
           !connect vertices from j_sphere cuttings which
           !are the next ones in the opposite (cyclic) direction
           !over the vertices from where
           !the j_sphere cutting edge joins two of them
           !the is-field is constructed in 1 direction
           !therefore in 1 direction is(2) has to be connected with
           !is(3) ... is(1) with is(num), whereas in the 2 direction
           !is(1) has to be connected with is(2) and so on
           ii=1
           tessera_parts: do n_tes_part=1,max_tes_part
             if(num_v(n_tes_part)<4) then
               ! part is not devided in several new parts
               ii=ii+num_v(n_tes_part)
               cycle tessera_parts
             endif

             !find direction
             l1=2
             l2=1
             if(poli%n_sphere(is(ii),1)==j) then
                l1=1
                l2=2
             endif
             do i=ii-1+l1,ii-1+num_v(n_tes_part),2
                i1=i-1
                if(i1<ii) i1=ii+num_v(n_tes_part)-1
                poli%bounds(is(i1),2)=is(i)
                poli%bounds(is(i),1)=is(i1)
                poli%n_sphere(is(i1),2)=j
                poli%n_sphere(is(i),1)=j
                if(h_case_sphere) then
                        h_case(is(i1),2)=.true.
                        h_case(is(i),1)=.true.
                endif
                poli%xyz_bound(is(i),1:3)=cut_center(:)
                poli%xyz_bound(is(i1),4:6)=cut_center(:)
                poli%r_bound(is(i),1)=cut_radius
                poli%r_bound(is(i1),2)=cut_radius
             enddo
             ii=ii+num_v(n_tes_part)
           enddo tessera_parts
         endif

         ! contract poligon (remove empty space)
         ! also only over old corners
!       if(sphere_number.eq.2.and.poli%n_vertises.eq.4) print*,'n_vertises', poli%n_vertises
         i=1
         do
             if(i>poli%n_vertises) exit
!               if(i.eq.0) print*,'i eq 0'
             if(inside_flag(i)) then
                do k=i+1,poli%n_vertises+k_new
                  if(.not. inside_flag(k)) then
                   poli%xyz_vertex(k-1,:)=poli%xyz_vertex(k,:)
                   poli%r_bound(k-1,:)=poli%r_bound(k,:)
                   poli%xyz_bound(k-1,:)=poli%xyz_bound(k,:)
                   poli%n_sphere(k-1,:)=poli%n_sphere(k,:)
                   h_case(k-1,:)=h_case(k,:)
                   poli%bounds(poli%bounds(k,1),2)=k-1
                   poli%bounds(poli%bounds(k,2),1)=k-1
                  endif
                  poli%bounds(k-1,:)=poli%bounds(k,:)
                  inside_flag(k-1)=inside_flag(k)
                  extra_point(k-1)=extra_point(k)
                enddo
                poli%n_vertises=poli%n_vertises-1
                i=i-1
             endif
             i=i+1
         enddo
         poli%n_vertises=poli%n_vertises+k_new
      enddo j_spheres

!       if(sphere_number.eq.2.and.poli%n_vertises.eq.4)  &
!            print*, sum(abs(poli%xyz_vertex(1:poli%n_vertises,:))), &
!             'just defined',k_new,sum(abs(poli%xyz_vertex(1:3,:)))

      !count again number of tessera pieces, because some
      !may be totally inside an other sphere
      ii=0
      is=-1
      n_tes_part=1

      count_pieces: do k=1,poli%n_vertises
            do i=1,ii
              if(is(i)==k) cycle count_pieces
            enddo
            ii=ii+1
            is(ii)=k
            i1=k
            do l1=1,poli%n_vertises
               i1=poli%bounds(i1,1)
               if (i1==k) exit
               ii=ii+1
               is(ii)=i1
            enddo
            n_tes_part=n_tes_part+1
      enddo count_pieces
      max_tes_part=n_tes_part-1


      !test for half-circle edges or edges with angle > pi
      k=poli%n_vertises+1
      n_extra=0
      do i=1,poli%n_vertises
         if(poli%n_sphere(i,1)==0) cycle
         if(h_case(i,1)) cycle
         i1=poli%bounds(i,1)
         v_buf=vector_product(poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:),&
                        poli%xyz_vertex(i1,:)-xyz_sphere(sphere_number,:))
         d_vv= DOT3(v_buf,poli%xyz_vertex(i,:)-poli%xyz_bound(i,1:3))

         if(d_vv <small_d .and. -d_vv<small_d) then
           v_buf=v_buf/LEN3(v_buf)*poli%r_bound(i,1)+&
                        poli%xyz_bound(i,1:3)
         else if( &
           DOT3(v_buf, 0.5_r8_kind*(poli%xyz_vertex(i,:)+poli%xyz_vertex(i1,:))-poli%xyz_bound(i,1:3)) &
           < small_d ) then
           v_buf= -0.5_r8_kind*poli%xyz_vertex(i,:) &
                   -0.5_r8_kind*poli%xyz_vertex(i1,:) &
                   +poli%xyz_bound(i,1:3)
           v_buf= v_buf/LEN3(v_buf)*poli%r_bound(i,1)&
                   +poli%xyz_bound(i,1:3)
         else
           cycle
         endif
         poli%xyz_vertex(k,:)=v_buf


         poli%bounds(k,1)=i1
         poli%bounds(k,2)=i
         poli%bounds(i1,2)=k
         poli%bounds(i,1)=k

         poli%r_bound(k,1) = poli%r_bound(i,1)
         poli%r_bound(k,2) = poli%r_bound(i,1)
         poli%xyz_bound(k,1:3) = poli%xyz_bound(i,1:3)
         poli%xyz_bound(k,4:6) = poli%xyz_bound(i,1:3)
         poli%n_sphere(k,1)=poli%n_sphere(i,1)
         poli%n_sphere(k,2)=poli%n_sphere(i,1)
         h_case(k,1)=h_case(i,1)
         h_case(k,2)=h_case(i,1)
         extra_point(k)=.true.

         k=k+1
      enddo
      poli%n_vertises=k-1

      ! sorting already enshured
      ! test sorting
      do i=1,poli%n_vertises
        i1=poli%bounds(i,1)
        if(poli%bounds(i1,2)/=i) print *,"wrong sorting"
        i1=poli%bounds(i,2)
        if(poli%bounds(i1,1)/=i) print *,"wrong sorting"
      enddo

      if(poli%n_vertises==0) then
        poligon_area=0.0_r8_kind
        return
      endif

      ! computing poligone center
      xyz_poligon_centre=0.0_r8_kind
      if(.not. weight_cent) then
       if(.not. orig_cent .and. .not. weight_cent) then
!       if(sphere_number.eq.2) print*,'case 1 xyz_poligon_centre'
         do i=1,poli%n_vertises
          if(.not. extra_point(i)) then
            xyz_poligon_centre=xyz_poligon_centre+ & !(case 1)
              (poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:))
          endif
         enddo
        else if(orig_cent) then
!       if(sphere_number.eq.2) print*,'case 2 xyz_poligon_centre'
        ! use the original triangle vertices -> no problem with gradient
            xyz_poligon_centre=xyz_triangl_1+xyz_triangl_2+xyz_triangl_3 & !(case 2)
                        -3.0_r8_kind*xyz_sphere(sphere_number,:)
        endif
        center_factor=LEN3(xyz_poligon_centre)
        xyz_poligon_centre=xyz_poligon_centre*(r_sphere(sphere_number))/ &
           center_factor+xyz_sphere(sphere_number,:)
!!!! MF bug fix 9/2001 (for further using)
        center_factor=center_factor/r_sphere(sphere_number)
      endif
!     else  if weigt_cent -> calculate after area computation

      ! computing poligon area (Gauss-Bonet theorem)
      ! attention: there are lots of errors in the paper of Cossi & al: J Comp Chem 17,57-73 (1996)
      ! read (p 60):
      ! phi_n=arccos(hat_v(n) hat_v*(n+1))
      ! cos(theta_n)= hat_P(n) hat_T(n), where
      ! T(n) is the vector pointing form the center of actual sphere to
      ! the center of the cutting sphere (the  difference to the definition in the paper
      ! is a sign in the case, one sphere is more then half in the other)

      ! more or less pathologic cases lead to negative areas
      ! in case, when one geodetic bonding edge of th poligon
      ! is cutted twice by the same sphere
      ! (e.g. if this cutting sphere is very small, or the poligon is
      ! very flat
      ! 4/2000 neglegt this areas (better: contract them with
      ! an neigbouring poligon
      combination=.false.

      poligon_area=0.0_r8_kind
      phin=0.0_r8_kind
      costhn=0.0_r8_kind

!!!MF >>>>
      first_vertex: do i=1,poli%n_vertises
         j=1
            i1=poli%bounds(i,j)

            if(poli%n_sphere(i,j)/=0) then
               v_buf=xyz_sphere(poli%n_sphere(i,j),:)-xyz_sphere(sphere_number,:)
               d_rad=LEN3(v_buf)
               cos_gamma=DOT3(poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:),v_buf)/&
                       (r_sphere(sphere_number)*d_rad)

               alpha=acosnum( &
                 DOT3(poli%xyz_vertex(i,:)-poli%xyz_bound(i,3*j-2:3*j), poli%xyz_vertex(i1,:)-poli%xyz_bound(i,3*j-2:3*j)) &
                 / poli%r_bound(i,j)**2)

               ! because of the new sorting of vertices, always with neigbour
               !        number 1
               !if(j==1) then
                  costhn(i)=cos_gamma
                  phin(i)=alpha
               !else
               !   costhn(i1)=cos_gamma
               !   phin(i1)=alpha
               !endif
               poligon_area=poligon_area+alpha*cos_gamma
            else if(weight_cent_mass) then !also phin needed for uncut bounds
                                      !i.e. bound= sphere centre
               phin(i)=acosnum( &
                 DOT3(poli%xyz_vertex(i,:)-poli%xyz_bound(i,3*j-2:3*j), poli%xyz_vertex(i1,:)-poli%xyz_bound(i,3*j-2:3*j)) &
                 / poli%r_bound(i,j)**2)
            endif
      enddo first_vertex

!       if(sphere_number.eq.2.and.poli%n_vertises.eq.4) print*,sum(phin(1:poli%n_vertises)), &
!               sum(abs(poli%xyz_vertex(1:poli%n_vertises,:)))
!!!MF <<<<

      do i=1,poli%n_vertises
         i1=poli%bounds(i,1)
         vn(i,:)=poli%xyz_vertex(i,:)-poli%xyz_bound(i,1:3)
         v_buf=vector_product(vn(i,:), &
              poli%xyz_vertex(i1,:)-poli%xyz_bound(i,1:3))
         v1=vector_product(vn(i,:),v_buf)

         v1=v1/LEN3(v1)
         i1=poli%bounds(i,2)
         vns(i,:)=poli%xyz_vertex(i,:)-poli%xyz_bound(i,4:6)

         v_buf=vector_product(vns(i,:), &
              poli%xyz_vertex(i1,:)-poli%xyz_bound(i,4:6))

         v2=vector_product(vns(i,:),v_buf)

         v2=v2/LEN3(v2)
         if(.not. extra_point(i)) then
                beta=pi-acosnum(DOT3(v1,v2))

                poligon_area=poligon_area-beta
         endif

         !v1,v2 are the tn,tn+1 in Cossi96
         if (do_gradients .or. weight_cent_sin) then
                tn(i,:)=v1
                tns(i,:)=v2
         endif
      enddo

      poligon_area=(poligon_area+2.0_r8_kind*pi*max_tes_part)*r_sphere(sphere_number)**2

      if(weight_cent_sin) then
!         if(sphere_number.eq.2) print*,'case 3 xyz_poligon_centre'
         do i=1,poli%n_vertises
          if(.not. extra_point(i)) then
            xyz_poligon_centre=xyz_poligon_centre+ &    ! (case 3)
              (poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:)) * &
              (1.0_r8_kind-(DOT3(tn(i,:),tns(i,:)))**2)
          endif
         enddo
         center_factor=LEN3(xyz_poligon_centre)
         xyz_poligon_centre=xyz_poligon_centre*(r_sphere(sphere_number))/ &
           center_factor+ xyz_sphere(sphere_number,:)
!!!!! MF 9/2001 bug fix (for further using)
         center_factor=center_factor/r_sphere(sphere_number)

      else if(weight_cent_mass) then
         !mass=0.0_r8_kind
         !real total mass not needed, for another scaling is anyhow
         !done to put the representative point (center) on the sphere
         !again
         do i=1,poli%n_vertises
            i1=poli%bounds(i,1)
            !run also over extra points which are constructed
            !to avoid sign problems with angles larger than pi
            !mass=mass+poli%r_bound(i,1)*phin(i)

            tanphi2=tan(phin(i)/2.0_r8_kind)

!           xyz_poligon_centre=xyz_poligon_centre+ &
!               poli%r_bound(i,1)*(&
!                       (phin(i)-2.0_r8_kind*tanphi2)*&
!                          (poli%xyz_bound(i,1:3)-xyz_sphere(sphere_number,:))+&
!                       tanphi2*(poli%xyz_vertex(i,:)+poli%xyz_vertex(i1,:)-&
!                               2.0_r8_kind*xyz_sphere(sphere_number,:))&
!               )
!           easier:

            xyz_poligon_centre=xyz_poligon_centre+ &   !(case 4)
                poli%r_bound(i,1)*(&
                        phin(i)*(poli%xyz_bound(i,1:3)-xyz_sphere(sphere_number,:))+&
                        tanphi2*(vn(i,:)+vns(i1,:))&
                )
         enddo
         center_factor=LEN3(xyz_poligon_centre)
!       if(sphere_number.eq.2) then
!         print*,'case 4 xyz_poligon_centre',poli%n_vertises,center_factor
!        print*,sum(abs(poli%xyz_bound(1:poli%n_vertises,:))),sum(abs(phin(1:poli%n_vertises)))
!       endif

         xyz_poligon_centre=xyz_poligon_centre*(r_sphere(sphere_number))/ &
                center_factor+xyz_sphere(sphere_number,:)
!!!!! MF 9/2001 bug fix (for further using)
         center_factor=center_factor/r_sphere(sphere_number)
      endif

!!! MF 4/2000 cancle negative areas
      if (poligon_area < min_area) then
        DPRINT 'area small/negative (',poligon_area,') , canceled'
        DPRINT 'at center',xyz_poligon_centre
        poligon_area=0.0_r8_kind
      endif

!!! MF 5/2000
      if(do_gradients .and. poligon_area>0.0_r8_kind ) then
         if(cutted) then
            call geom_grad_cutted_triang(poli,sphere_number,tn,tns,&
                vn,vns,costhn,phin,poligon_area,extra_point,&
                xyz_poligon_centre,center_factor)
         else
            call geom_grad_orig_triang(sphere_number,poligon_area,&
                                xyz_poligon_centre)
         endif
      endif
!if(poligon_area>0.0_r8_kind .and. .not.do_cavitation .and. .not. do_disp_rep ) then
!n_extra=0
!i1=0
!do i=1,poli%n_vertises
!if(extra_point(i)) n_extra=n_extra+1
!if(poli%n_sphere(i,1)/=0) i1=i1+1
! enddo
! DPRINT 'poligon ' ,N_total+1,poli%n_vertises,i1,n_extra,max_tes_part
!endif

    end subroutine def_poligon_tes
!---------------------------------------------------------------------------
    !------------------------------------------------------------
!!!MF numerical save acos function >>>>
    function acosnum(cos_t)
      real(kind=r8_kind), intent(in) :: cos_t
      real(kind=r8_kind) :: acosnum

!!$print*,cos_t
      if(cos_t*cos_t-1.0_r8_kind > 1.0e-15_r8_kind) then
        DPRINT cos_t
        DPRINT counter_n
        call error_handler("acos of value >1 or <-1 demanded?!")
      else if(cos_t>1.0_r8_kind) then
        acosnum=acos(1.0_r8_kind)
      else if(cos_t<-1.0_r8_kind) then
        acosnum=acos(-1.0_r8_kind)
      else
        acosnum=acos(cos_t)
      end if

    end function acosnum
    !------------------------------------------------------------

!---------------------------------------------------------------------------
    subroutine line_cut(p1, p2, cent1, r1, cent2, r2, ins1i, ins2i, cent0, h_case, ncut, pcut)
      ! idea: any point of the edge considered (center at origin)
      ! can be written as:
      ! r=(v1+lam (v2-v1))/|(v1+lam (v2-v1))| *r_bound
      ! with lam in [0,1]
      ! look for solutions of |r-d|=R, where d is the vector to the
      ! cutting sphere center and R is its radius
      implicit none
      real(kind=r8_kind), intent(in) :: p1(:), p2(:), cent1(:), cent2(:), cent0(:) ! (3)
      real(kind=r8_kind), intent(in) :: r1, r2
      logical, intent(in) :: ins1i, ins2i, h_case
      integer (kind=i4_kind), intent(out) :: ncut
      real(kind=r8_kind), intent(out), optional :: pcut(6)
      !** End of interface *****************************************

      integer (kind=i4_kind) :: i,ncut_o,ih
      real(kind=r8_kind), dimension(3) :: v1,v12,vd
      real(kind=r8_kind) :: d,del,v1v12,v1vd,v12vd,ha,hb,hc,lam(2),v12v12,v1v1,aha,ahb,ahc
      real(kind=r8_kind), parameter :: small = 1.0e-8_r8_kind
      real (kind=r8_kind) :: vb(3),help,lower_bound,pcut_help(6)
      logical :: near(2,2),l_help,twocase,colin_case,smaller_pi,second_run,ins1,ins2,interchange

      smaller_pi=.true.
      second_run=.false.
      colin_case=.false.
      twocase=.false.
      interchange=.false.
      ins1=ins1i
      ins2=ins2i
      lower_bound=0.0_r8_kind
      lam=-100.0_r8_kind
      v1=p1-cent1
      v12=p2-p1

      !test if vectors are linear dependent
      vd=vector_product(v1,v12)
!!$print*,v1,v12
      if(DOT3(vd,vd)<small .and. .not. DOT3(v12,v12) < small)&
                colin_case=.true. !vectors parallel!

      !test if edge has angle larger then pi -> only for cutting edges
      if(.not.colin_case) then
        vd=cent1-cent0
        if(DOT3(vd,vd)>small .and. .not. h_case ) then
           vd=vector_product(p1-cent0,p2-cent0)
!          larger pi case for h_case
           if(h_case) vd=-vd
           if(DOT3(vd,0.5_r8_kind*(p1-cent1)+0.5_r8_kind*(p2-cent1))&
                        <0.0_r8_kind ) then
              smaller_pi=.false.
              vb= -0.5_r8_kind*p1 -0.5_r8_kind*p2 +cent1
              vb= vb/LEN3(vb) *r1
              v12= vb -v1
              vd=vb-(cent2-cent1)
              if(DOT3(vd,vd)<r2**2) then
                        ins2=.true.
              else
                        ins2=.false.
              endif
           endif
        endif
      endif

      if(colin_case) then
        v12=vector_product(p1-cent0,p2-cent0)
        v12=v12/LEN3(v12) * r1
        ! opposite for h_case (convex bond)
!       if (h_case) v12=-v12
        lower_bound=-1.0_r8_kind
      endif

      vd=cent2-cent1
      ncut=0

      d=LEN3(vd)
      del=r1**2+d**2-r2**2

1000  v1v12=DOT3(v1,v12)
      v1vd=DOT3(v1,vd)
      v12vd=DOT3(v12,vd)
      v12v12=DOT3(v12,v12)
      v1v1=DOT3(v1,v1)
      near=.false.

      if(.not. colin_case) then
         ha=4.0_r8_kind*r1**2*v12vd**2-del**2*v12v12
         hb=8.0_r8_kind*r1**2*v12vd*v1vd-2.0_r8_kind*del**2*v1v12
         hc=4.0_r8_kind*r1**2*v1vd**2-del**2*v1v1
      else
         ha=4.0_r8_kind*(v1vd**2+v12vd**2)
         hb=4.0_r8_kind*(-del*v1vd)
         hc=del**2-4.0_r8_kind*v12vd**2
      endif

      aha=abs(ha)
      ahb=abs(hb)
      ahc=abs(hc)
!!$print*,aha,ahb,ahc,colin_case
      ! this a littel bit complcated formulation I do because of
      ! not clear accuracy of intrinsic functions
      ! solutions of: ha*lam**2+hb*lam+hc==0
      if(aha<small .and. ahb<small .and. ahc<small) then
         !all lambdas fullfill the equation (esotheric case)
         ncut=3
         return
      else if (aha<small .and. ahb<small ) then
         !no solution
         ncut=0
         if (.not. colin_case .and. smaller_pi) return
      else if (aha<small ) then
         !exactly one solution
         ncut=1
         lam(1)=-hc/hb
         if(lam(1)<lower_bound-small .or. lam(1)>1.0_r8_kind+small) ncut=0
         if(abs(lam(1)-lower_bound)<10.0_r8_kind*small) then
                if(ins1.neqv.ins2) then
                    lam(1)=lower_bound
                    ncut=1
                else
                    ncut=0
                endif
         else if(abs(lam(1)-1.0_r8_kind)<10.0_r8_kind*small) then
                if(ins1.neqv.ins2) then
                    lam(1)=1.0_r8_kind
                    ncut=1
                else
                    ncut=0
                endif
         endif
      else if ((hb/ha)**2-4.0_r8_kind*hc/ha> -small) then
!!$print*,'2'
         !two possible solutions
         twocase = .true.
         ncut=2
         if((hb/(2.0_r8_kind*ha))**2-hc/ha > small*10.0_r8_kind) then
            lam(1)=sqrt((hb/(2.0_r8_kind*ha))**2 - hc/ha)
         else
            lam(1)=0.0_r8_kind
         endif
         lam(2)=-hb/(2.0_r8_kind*ha)+lam(1)
         lam(1)=-hb/(2.0_r8_kind*ha)-lam(1)
         near(1,1)=(abs(lam(1)-lower_bound)<small*10.0_r8_kind)
         near(2,1)=(abs(lam(2)-lower_bound)<small*10.0_r8_kind)
         near(1,2)=(abs(lam(1)-1.0_r8_kind)<small*10.0_r8_kind)
         near(2,2)=(abs(lam(2)-1.0_r8_kind)<small*10.0_r8_kind)
         ! always lam(1)<lam(2)
!!$print*,lam(1),lam(2),'lam(1),lam(2)',lower_bound
         if(lam(1)>1.0_r8_kind+small .or. lam(2)<lower_bound-small) then
            ! no solution in the range [0,1]
            ncut=0
         else if(lam(2)>1.0_r8_kind+small .and. lam(1)<lower_bound-small)  then
            ! no solution in the range [0,1]
            ncut=0
         else if(lam(2)>1.0_r8_kind+small) then! .or.
            ! only lam(1) in the range [0,1]
            ncut=1
         else if(lam(1)<lower_bound-small) then!.or. ((ins1i.neqv.ins2i).and.lam(2)>0.0_r8_kind)) then
            ! only lam(2) solution in range [0,1]
            ncut=1
            help=lam(2)
            lam(2)=lam(1)
            lam(1)=help
            l_help=near(1,1)
            near(1,1)=near(2,1)
            near(2,1)=l_help
            l_help=near(1,2)
            near(1,2)=near(2,2)
            near(2,2)=l_help
         endif
!!$print*,ncut,'ncut'
      else
         !no real solution
         ncut=0
         if(.not. colin_case .and. smaller_pi) return
      endif
      i=1
      do
         if(i>ncut) exit

         if(.not.colin_case) then
           vb=v1+lam(i)*v12
           vb=vb/(LEN3(vb))
         else
           if(lam(i)**2<1.0_r8_kind) then
              vb=lam(i)*v1+sqrt(1.0_r8_kind-lam(i)**2)*v12
           else
              vb=lam(i)*v1
           endif
           vb=vb/r1
         endif
         ! test if not solution of - vd*pcut = d**2+r1**2-r2**2
!!$         if(abs(-2.0_r8_kind*DOT3(vd,vb)*r1-del)<sqrt(small) .and. &
         if(abs(2.0_r8_kind*DOT3(vd,vb)*r1-del)>=sqrt(small) .and. &
                abs(DOT3(vd,vb))>sqrt(small)) then
            ! non valid solution
            if(i==1 .and. ncut==1) then
               ncut=0
               if(.not. colin_case .and. smaller_pi) return
            else if(i==1 .and. ncut==2) then
               i=0
               ncut=1
               help=lam(1)
               lam(1)=lam(2)
               lam(2)=help
               l_help=near(1,1)
               near(1,1)=near(2,1)
               near(2,1)=l_help
               l_help=near(1,2)
               near(1,2)=near(2,2)
               near(2,2)=l_help
               interchange=.true.
            else
               ncut=1
            endif
         else
            ! valid solution
            pcut(3*i-2:3*i) = cent1 + vb*r1
         endif
         i=i+1
      enddo

!!$print*,ncut,'ncut1'
      if(twocase) then
         if(ncut==0 .and. (ins1 .neqv. ins2)) then
            if(near(1,1) .or. near(2,1)) then
               lam(1)=lower_bound
            else if(near(1,2) .or. near(2,2)) then
               lam(1)=1.0_r8_kind
            else
               ncut=-1
            endif
            ncut=ncut+1
         else if(ncut==2 .and. (ins1 .neqv. ins2)) then
            if(ins1 .and. near(1,1)) then
               lam(1)=lam(2)
            else if (ins1 .and. near(2,1)) then
            else if (ins2 .and. near(2,2)) then
            else if (ins2 .and. near(1,2)) then
               lam(1)=lam(2)
            else if (near(1,1) .and. near(2,1)) then
            else if (near(1,2) .and. near(2,2)) then
               lam(1)=lam(2)
            else if(abs(lam(1)-lam(2))<sqrt(small)) then
            else
               ncut=3
            endif
            ncut=ncut-1
         else if(ncut ==1 .and. (ins1 .eqv. ins2) .and. ins1) then
            if (near(1,1) .or. near(1,2) .or. near (2,1) .or. near (2,2))ncut=0
         else if(ncut ==1 .and. (ins1 .eqv. ins2) .and. .not. ins1) then
            if (near(1,1) .or. near(1,2)) then
               ncut=0
            else if (near(2,1) .or. near(2,2)) then
               ncut=2
               if(near(2,1)) then
                  lam(2)=lower_bound
                  interchange=.true.
               else if(near(2,2)) then
                  lam(2)=1.0_r8_kind
               endif
            endif
         endif
!!$print*,ncut,'ncut2'

!        recalculate in any case
         do i=1,ncut
           if(.not.colin_case) then
             vb=v1+lam(i)*v12
             vb=vb/(LEN3(vb))
           else
             if(lam(i)**2<1.0_r8_kind) then
              vb=lam(i)*v1+sqrt(1.0_r8_kind-lam(i)**2)*v12
             else
              vb=lam(i)*v1
             endif
           vb=vb/r1
           endif
           pcut(3*i-2:3*i) = cent1 + vb*r1
         enddo
      endif

      if(.not.smaller_pi .and. ncut<2 .and. .not. second_run) then
           pcut_help=pcut
           ncut_o=ncut
           v1=-0.5_r8_kind*p1 -0.5_r8_kind*p2 + cent1
           v1=v1/LEN3(v1)*r1
           v12=p2-cent1-v1
           ins1=ins2
           ins2=ins2i
           lam=-100.0_r8_kind
           second_run=.true.
           goto 1000
      else if (.not. smaller_pi .and. second_run) then
           ncut=ncut+ncut_o
           do i=ncut_o+1,ncut
             ih=i-ncut_o
             pcut_help(3*i-2:3*i)=pcut(3*ih-2:3*ih)
           enddo
           pcut=pcut_help
      endif
!!$print*,ncut,'ncut3'

!     if(ncut==0 .and. colin_case) then
!        pcut(1:3)=cent1+v12
!        ncut=11
!     endif
!     if(ncut==0 .and. .not. smaller_pi) then
!        vb= -0.5_r8_kind*p1 -0.5_r8_kind*p2 +cent1
!        vb= vb/DOT3(vb,vb) *r1
!        pcut(1:3)=cent1+vb
!        ncut=11
!     endif

    end subroutine line_cut
    !------------------------------------------------------------

    !------------------------------------------------------------
    function vector_product(vector1, vector2)
      implicit none
      real(kind=r8_kind), intent(in) :: vector1(:), vector2(:)
      real(kind=r8_kind) :: vector_product(3)
      ! *** end of interface ***

      vector_product(1)=vector1(2)*vector2(3)-vector2(2)*vector1(3)
      vector_product(2)=vector1(3)*vector2(1)-vector2(3)*vector1(1)
      vector_product(3)=vector1(1)*vector2(2)-vector2(1)*vector1(2)
    end function vector_product
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine geom_grad_cutted_triang(poli,nsp,tn,tns,vn,vns,costhn,phin,area,&
         extra_point,xyz_cent,cent_f)
      !gradients of areas and centers of cutted triangles
      !** End of interface *****************************************

      type(poligon), intent(in) :: poli
      integer(kind=i4_kind), intent(in) :: nsp
      real(kind=r8_kind), intent(in),dimension(MAX_POL_VER,3) :: tn,tns,vn,vns
      real(kind=r8_kind), intent(in),dimension(MAX_POL_VER) :: phin,costhn
      real(kind=r8_kind), intent(in) :: area,xyz_cent(3),cent_f
      logical, intent(in) :: extra_point(MAX_POL_VER)

      real(kind=r8_kind), dimension(3) :: vc,vc1,vd,vd1,vcc,vdd,tnsg,tng!,aa
      real(kind=r8_kind), dimension(3) :: ep,ep1,no,nos,vbuf,vbufs,vbuf1,dno,dnos
      real(kind=r8_kind), dimension(3) :: dik,diks,epk,dep,d_no,d_nos,dvbuf
      real(kind=r8_kind) :: dp1,dpnoep,dpnosep,dpnotns,dpnostn,mu,mus,dmu,dmus,f,df,df1,d2f
      real(kind=r8_kind) :: dist1,dist2,help,help1,tanphi2,rr2,drr2,drr21
      real(r8_kind) :: d2phi,d2costh,d2omeg
      real(kind=r8_kind), dimension(3) :: dP,dTc,dTcs,d2P,e_cent,dtn,dtns
      real(kind=r8_kind), dimension(3) :: vvv,dvvv,dvvv1,d2vvv,vv,dvv,d2Tc,d2Tcs,d2no,d2nos
      real(kind=r8_kind), dimension(3) :: v1,v2,dv11,dv12,dv21,dv22,d2v1,d2v2,d2vb1,d2vb2
      integer(kind=i4_kind) :: i,j,l,i1,i2,k,status,na,ea,l1 !,k1
      integer(i4_kind) :: m,m1,n,na1,ea1 !,ii
      logical :: b_cut,bs_cut,change_sign,change_signs
      real (kind=r8_kind), parameter :: small=1.0e-11_r8_kind
      integer :: N_pc

      real(r8_kind), allocatable :: mu_d(:,:,:,:),mus_d(:,:,:,:),dTc_d(:,:,:,:,:),dTcs_d(:,:,:,:,:)
      type(arrmat2), allocatable :: cgdP(:,:,:),cgdv(:,:,:),cgdvs(:,:,:)
      type(arrmat2), allocatable :: cgd2T(:,:,:,:,:,:),cgd2Ts(:,:,:,:,:,:),cgd2P(:,:,:,:,:,:)
      type(arrmat1), allocatable :: cgdomeg(:,:,:),cgdcosth(:,:,:),cgdphi(:,:,:)

      !bound is the edge and vertex the corner of the poligon
      !all derivatives relative to sphere center
      !dP: vertex relative to sphere center
      !dv: vertex relative to bound center
      !dTc: bound center
      !tn.. bound tangents at vertices
      !phin: bound arc angle

      if(VTN) return

      if(integralpar_2dervs) then
          if ( weight_cent .or. weight_cent_mass .or. weight_cent_sin .or. orig_cent) then
              print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
              print *, "XX                                                            XX"
              print *, "XX       F E A T U R E  N O T  I M P L E M E N T E D !        XX"
              print *, "XX                                                            XX"
              print *, "XX Use of weighted center of tessarea is not yet implemented  XX"
              print *, "XX     in combination with second energy derivatives.         XX"
              print *, "XX                                                            XX"
              print *, "XX Use cent_type='corners' or cent_type='triangles' instead!  XX"
              print *, "XX                                                            XX"
              print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
              call error_handler('cavity 2nd dervs: switch off all weight_cent_... and orig_cent parameters ')
          endif
      endif

      N_pc=0
!!$      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
      if(with_pc) N_pc=pointcharge_N

      allocate(cgdP(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang0")
      allocate(cgdv(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang1")
      allocate(cgdvs(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang2")
      allocate(cgdomeg(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang3")
      allocate(cgdcosth(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang4")
      allocate(cgdphi(3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
      if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang5")
      if(integralpar_2dervs)then
         allocate(mu_d(poli%n_vertises,3,N_moving_unique_atoms+N_pc,ua_dim_max), &
              mus_d(poli%n_vertises,3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
         ASSERT(status.eq.0)
         allocate(dTc_d(poli%n_vertises,3,3,N_moving_unique_atoms+N_pc,ua_dim_max), &
              dTcs_d(poli%n_vertises,3,3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
         ASSERT(status.eq.0)
         allocate(cgd2T(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
              3,N_moving_unique_atoms+N_pc,ua_dim_max), &
              cgd2Ts(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
              3,N_moving_unique_atoms+N_pc,ua_dim_max), &
              cgd2P(3,N_moving_unique_atoms+N_pc,ua_dim_max, &
              3,N_moving_unique_atoms+N_pc,ua_dim_max),stat=status)
         ASSERT(status.eq.0)
      end if

      do i=1,N_moving_unique_atoms+N_pc
         if(i <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(i)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=i-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l=1,ea
            do j=1,3
               allocate(cgdP(j,i,l)%m(3,poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang6")
               allocate(cgdv(j,i,l)%m(3,poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang7")
               allocate(cgdvs(j,i,l)%m(3,poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang8")
               allocate(cgdomeg(j,i,l)%m(poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang9")
               allocate(cgdcosth(j,i,l)%m(poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang10")
               allocate(cgdphi(j,i,l)%m(poli%n_vertises),stat=status)
               if(status/=0) call error_handler("alloc fail geom_grad_cutted_triang11")
               cgdP(j,i,l)%m=0.0_r8_kind
               cgdv(j,i,l)%m=0.0_r8_kind
               cgdvs(j,i,l)%m=0.0_r8_kind
               cgdcosth(j,i,l)%m=0.0_r8_kind
               cgdphi(j,i,l)%m=0.0_r8_kind
               cgdomeg(j,i,l)%m=0.0_r8_kind
               if(integralpar_2dervs)then
                  do k=1,N_moving_unique_atoms+N_pc
                     if(i <= N_moving_unique_atoms) then
                        na1=moving_unique_atom_index(k)
                        ea1=unique_atoms(na1)%n_equal_atoms
                     else
                        na1=k-N_moving_unique_atoms
                        ea1=pointcharge_array(na1)%N_equal_charges
                     end if
                     do m=1,ea1
                        do n=1,3
                           allocate(cgd2T(j,i,l,n,k,m)%m(3,poli%n_vertises),stat=status)
                           ASSERT(status.eq.0)
                           allocate(cgd2Ts(j,i,l,n,k,m)%m(3,poli%n_vertises),stat=status)
                           ASSERT(status.eq.0)
                           allocate(cgd2P(j,i,l,n,k,m)%m(3,poli%n_vertises),stat=status)
                           ASSERT(status.eq.0)
                           cgd2T(j,i,l,n,k,m)%m=0.0_r8_kind
                           cgd2Ts(j,i,l,n,k,m)%m=0.0_r8_kind
                           cgd2P(j,i,l,n,k,m)%m=0.0_r8_kind
                        end do
                     end do
                  end do
               end if
            enddo
         enddo
      enddo

      do i=1,poli%n_vertises
         !P here relative vertex vector
         ep1=poli%xyz_vertex(i,:)-xyz_sphere(nsp,:)
         ep=ep1/LEN3(ep1)

         i1=poli%bounds(i,1)
         i2=poli%bounds(i,2)

         vbuf=vector_product(vn(i,:),vns(i1,:))
         no=vbuf/LEN3(vbuf)


         vbuf=vector_product(vns(i,:),vn(i2,:))
         nos=vbuf/LEN3(vbuf)

         b_cut=.false.
         if(poli%n_sphere(i,1)/=0) b_cut=.true.
         bs_cut=.false.
         if(poli%n_sphere(i,2)/=0) bs_cut=.true.

         if(b_cut) then
            dik=xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:)
            dist1=LEN3(dik)
            if(DOT3(no,dik)<0.0_r8_kind) no=-no
         endif
         if(bs_cut) then
            diks=xyz_sphere(poli%n_sphere(i,2),:)-xyz_sphere(nsp,:)
            dist2=LEN3(diks)
            if(DOT3(nos,diks)<0.0_r8_kind) nos=-nos
         endif

         ! compute products after prospective sign change!!!
         dpnotns=DOT3(no,tns(i,:))
         dpnostn=DOT3(nos,tn(i,:))
         dpnoep =DOT3(no ,ep)
         dpnosep=DOT3(nos,ep)

         do l=1,N_moving_unique_atoms+N_pc
            if(l <= N_moving_unique_atoms) then
               na=moving_unique_atom_index(l)
               ea=unique_atoms(na)%n_equal_atoms
            else
               na=l-N_moving_unique_atoms
               if(with_pc) ea=pointcharge_array(na)%N_equal_charges
            end if
            do l1=1,ea
               do j=1,3
                  dTc(:)=0.0_r8_kind
                  dTcs(:)=0.0_r8_kind
                  dno=0.0_r8_kind
                  dnos=0.0_r8_kind

                  !calculate dno and dnos
                  !calculate dTc and dTcs
                  if(b_cut) then
                     do k=1,3
                        vbuf(k)= cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,l1) -&
                             cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                     enddo

                     help1=DOT3(dik,vbuf)
                     dno=vbuf/dist1-help1* dik/dist1**3

                     help=(r_sphere(nsp)**2- &
                          r_sphere(poli%n_sphere(i,1))**2)/dist1**2
                     dTc(:)=0.5_r8_kind*(1.0_r8_kind+help)*vbuf + &
                          dik*(r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,l1)-  &
                          r_sphere(poli%n_sphere(i,1)) *                &
                          cagr%dR(poli%n_sphere(i,1))%xyz_grad(j,l,l1)  &
                          -help1*help)/dist1**2
                  endif
                  if(bs_cut) then
                     do k=1,3
                        vbuf(k)= cagr%dc(k,poli%n_sphere(i,2))%xyz_grad(j,l,l1)-&
                             cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                     enddo

                     help1=DOT3(diks,vbuf)
                     dnos=vbuf/dist2-help1* diks/dist2**3

                     help=(r_sphere(nsp)**2- &
                          r_sphere(poli%n_sphere(i,2))**2)/dist2**2
                     dTcs(:)=0.5_r8_kind*(1.0_r8_kind+help)*vbuf + &
                          diks*(r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,l1)-  &
                          r_sphere(poli%n_sphere(i,2))*               &
                          cagr%dR(poli%n_sphere(i,2))%xyz_grad(j,l,l1)  -   &
                          help1*help)/dist2**2
                  endif

                  if(integralpar_2dervs)then
                     dTc_d(i,:,j,l,l1)=dTc
                     dTcs_d(i,:,j,l,l1)=dTcs
                  end if

                  if(.not.extra_point(i)) then
                    mus=(DOT3(no,dTc(:))-                        &
                       DOT3(dno,vn(i,:)) -                       &
                       dpnoep*cagr%dR(nsp)%xyz_grad(j,l,l1))/dpnotns
                    mu=(DOT3(nos,dTcs(:))-                       &
                       DOT3(dnos,vns(i,:)) -                     &
                       dpnosep*cagr%dR(nsp)%xyz_grad(j,l,l1))/dpnostn

                    dP(:)= ep * cagr%dR(nsp)%xyz_grad(j,l,l1) + &
                       tn(i,:) * mu + tns(i,:) *mus
                    if(integralpar_2dervs)then
                       mu_d(i,j,l,l1)=mu
                       mus_d(i,j,l,l1)=mus
                    end if
                  else
                    vbuf=vector_product(ep,tn(i,:))
                    mu=(DOT3(no,dTc(:))-&
                       DOT3(dno,vn(i,:)) - &
                       dpnoep*cagr%dR(nsp)%xyz_grad(j,l,l1))/&
                          DOT3(no,vbuf)
                    dP(:)= ep * cagr%dR(nsp)%xyz_grad(j,l,l1) + &
                       vbuf(:) * mu
                  endif

                  do k=1,3
                     cgdP(j,l,l1)%m(k,i)=dP(k)
                     if(.not.extra_point(i)) then
                       if(.not. weight_cent) then
                         cagr%dcenter(j,l,l1)%m(k,N_total+1)=dP(k) +&
                          cagr%dcenter(j,l,l1)%m(k,N_total+1)
                       else if(weight_cent_sin) then
                         cagr%dcenter(j,l,l1)%m(k,N_total+1)= &
                             cagr%dcenter(j,l,l1)%m(k,N_total+1)+&
                             dP(k)*(1.0_r8_kind-(DOT3(tn(i,:),tns(i,:)))**2)
                       endif
                     endif

                     cgdv(j,l,l1)%m(k,i)=dP(k)-dTc(k) !! use for recalculating dTc
                     cgdvs(j,l,l1)%m(k,i)=dP(k)-dTcs(k)

                  enddo
                  if(costhn(i)*costhn(i)>small) then
                     do k=1,3
                        vd(k)= cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,l1)-&
                             cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                     enddo
                     cgdcosth(j,l,l1)%m(i)= &
                          deriv_cos(poli%xyz_vertex(i,:)-xyz_sphere(nsp,:),&
                          xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:),&
                          dP,vd,.false.,.false.)
                  else
                     cgdcosth(j,l,l1)%m(i)= 0.0_r8_kind
                  endif
               enddo
            enddo
         enddo
      end do

      secnd_d: if(integralpar_2dervs)then
         do i=1,poli%n_vertises
            ep1=poli%xyz_vertex(i,:)-xyz_sphere(nsp,:)
            ep=ep1/LEN3(ep1)

            i1=poli%bounds(i,1)
            i2=poli%bounds(i,2)

            vbuf=vector_product(vn(i,:),vns(i1,:))
            no=vbuf/LEN3(vbuf)
            change_sign=.false.

            vbuf=vector_product(vns(i,:),vn(i2,:))
            nos=vbuf/LEN3(vbuf)
            change_signs=.false.

            b_cut=.false.
            if(poli%n_sphere(i,1)/=0) b_cut=.true.
            bs_cut=.false.
            if(poli%n_sphere(i,2)/=0) bs_cut=.true.

            if(b_cut) then
               dik=xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:)
               dist1=LEN3(dik)
               if(DOT3(no,dik)<0.0_r8_kind) then
                  no=-no; change_sign=.true.
               end if
            endif
            if(bs_cut) then
               diks=xyz_sphere(poli%n_sphere(i,2),:)-xyz_sphere(nsp,:)
               dist2=LEN3(diks)
               if(DOT3(nos,diks)<0.0_r8_kind) then
                  nos=-nos; change_signs=.true.
               end if
            endif

            ! compute products after prospective sign change!!!
            dpnotns=DOT3(no,tns(i,:))
            dpnostn=DOT3(nos,tn(i,:))
            dpnoep =DOT3(no ,ep)
            dpnosep=DOT3(nos,ep)

            do l=1,N_moving_unique_atoms+N_pc
               if(l <= N_moving_unique_atoms) then
                  na=moving_unique_atom_index(l)
                  ea=unique_atoms(na)%n_equal_atoms
               else
                  na=l-N_moving_unique_atoms
                  ea=pointcharge_array(na)%N_equal_charges
               end if
               do l1=1,ea
                  do j=1,3
                     dno=0.0_r8_kind
                     if(b_cut) then
                        do k=1,3
                           vbuf(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,l1)-cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                        end do
                        dno=vbuf/dist1-dik*DOT3(dik,vbuf)/dist1**3
                     end if
                     dnos=0.0_r8_kind
                     if(bs_cut) then
                        do k=1,3
                           vbufs(k)=cagr%dc(k,poli%n_sphere(i,2))%xyz_grad(j,l,l1)-cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                        end do
                        dnos=vbufs/dist2-diks*DOT3(diks,vbufs)/dist2**3
                     end if

!!$if(N_total+1 == 274)then
!!$   print*,i,l,l1,j,mu_d(i,j,l,l1),cgdP(j,l,l1)%m(:,i)
!!$endif
                     do m=1,N_moving_unique_atoms+N_pc
                        if(m <= N_moving_unique_atoms) then
                           na1=moving_unique_atom_index(m)
                           ea1=unique_atoms(na1)%n_equal_atoms
                        else
                           na1=m-N_moving_unique_atoms
                           ea1=pointcharge_array(na1)%N_equal_charges
                        end if
                        do m1=1,ea1
                           do n=1,3
                              dP(:)=cgdP(n,m,m1)%m(:,i)
                              dep=unit_vector_deriv(ep1,dP)
                              !--------------------------------------
                              vvv=vector_product(vns(i,:),vector_product(vns(i,:),vn(i2,:)))
                              dvvv=cgdvs(n,m,m1)%m(:,i)*DOT3(vns(i,:),vn(i2,:))+ &
                                   vns(i,:)*DOT3(cgdvs(n,m,m1)%m(:,i),vn(i2,:))+ &
                                   vns(i,:)*DOT3(vns(i,:),cgdv(n,m,m1)%m(:,i2))- &
                                   cgdv(n,m,m1)%m(:,i2)*DOT3(vns(i,:),vns(i,:))- &
                                   2.0_r8_kind*vn(i2,:)*DOT3(vns(i,:),cgdvs(n,m,m1)%m(:,i))
                              dtns=unit_vector_deriv(vvv,dvvv)
                              !--------------------------------------
                              vv=vector_product(vn(i,:),vns(i1,:))
                              dvv=vector_product(cgdv(n,m,m1)%m(:,i),vns(i1,:))+vector_product(vn(i,:),cgdvs(n,m,m1)%m(:,i1))
                              d_no=unit_vector_deriv(vv,dvv); if(change_sign) d_no=-d_no
                              !--------------------------------------
                              d2Tc(:)=0.0_r8_kind; d2no=0.0_r8_kind
                              if(b_cut) then
                                 rr2=r_sphere(nsp)**2-r_sphere(poli%n_sphere(i,1))**2

                                 f=(dist1**2+rr2)/(2.0_r8_kind*dist1**2)

                                 do k=1,3
                                    dvbuf(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_hess(j,l,l1,n,m,m1)- &
                                         cagr%dc(k,nsp)%xyz_hess(j,l,l1,n,m,m1)
                                    vbuf1(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(n,m,m1)-cagr%dc(k,nsp)%xyz_grad(n,m,m1)
                                 end do

                                 d2Tc(:)=d2Tc(:)+f*dvbuf(:)

                                 drr2=r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,l1)-  &
                                      r_sphere(poli%n_sphere(i,1))*cagr%dR(poli%n_sphere(i,1))%xyz_grad(j,l,l1)
                                 drr21=r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(n,m,m1)-  &
                                      r_sphere(poli%n_sphere(i,1))*cagr%dR(poli%n_sphere(i,1))%xyz_grad(n,m,m1)

                                 df=(drr2-rr2*DOT3(dik,vbuf)/dist1**2)/dist1**2
                                 df1=(drr21-rr2*DOT3(dik,vbuf1)/dist1**2)/dist1**2

                                 d2Tc(:)=d2Tc(:)+df1*vbuf(:)+df*vbuf1(:)

                                 d2f=(cagr%dR(nsp)%xyz_grad(n,m,m1)*cagr%dR(nsp)%xyz_grad(j,l,l1)+ &
                                      r_sphere(nsp)*cagr%dR(nsp)%xyz_hess(j,l,l1,n,m,m1)- &
                                      cagr%dR(poli%n_sphere(i,1))%xyz_grad(n,m,m1)*cagr%dR(poli%n_sphere(i,1))%xyz_grad(j,l,l1)- &
                                      r_sphere(poli%n_sphere(i,1))*cagr%dR(poli%n_sphere(i,1))%xyz_hess(j,l,l1,n,m,m1)- &
                                      2.0_r8_kind*drr21*DOT3(dik,vbuf)/dist1**2- &
                                      rr2*DOT3(vbuf1,vbuf)/dist1**2- &
                                      rr2*DOT3(dik,dvbuf)/dist1**2+ &
                                      2.0_r8_kind*rr2*DOT3(dik,vbuf)*DOT3(dik,vbuf1)/dist1**4)/dist1**2- &
                                      2.0_r8_kind*df*DOT3(dik,vbuf1)/dist1**2

                                 d2Tc(:)=d2Tc(:)+dik(:)*d2f

                                 cgd2T(j,l,l1,n,m,m1)%m(:,i)=d2Tc(:)
                                 !----------------------------------------
                                 d2no=dvbuf(:)/dist1-vbuf*DOT3(dik,vbuf1)/dist1**3- &
                                      vbuf1*DOT3(dik,vbuf)/dist1**3- &
                                      dik*DOT3(vbuf1,vbuf)/dist1**3-dik*DOT3(dik,dvbuf)/dist1**3+ &
                                      3.0_r8_kind*dik*DOT3(dik,vbuf)*DOT3(dik,vbuf1)/dist1**5
                              end if

                              dmus=(DOT3(d_no,dTc_d(i,:,j,l,l1))+DOT3(no,d2Tc)-DOT3(d2no,vn(i,:))- &
                                   DOT3(dno,cgdv(n,m,m1)%m(:,i))-DOT3(d_no,ep)*cagr%dR(nsp)%xyz_grad(j,l,l1)- &
                                   DOT3(no,dep)*cagr%dR(nsp)%xyz_grad(j,l,l1)- &
                                   DOT3(no,ep)*cagr%dR(nsp)%xyz_hess(j,l,l1,n,m,m1))/dpnotns

                              !===========================================
                              vvv=vector_product(vn(i,:),vector_product(vn(i,:),vns(i1,:)))
                              dvvv=cgdv(n,m,m1)%m(:,i)*DOT3(vn(i,:),vns(i1,:))+ &
                                   vn(i,:)*DOT3(cgdv(n,m,m1)%m(:,i),vns(i1,:))+ &
                                   vn(i,:)*DOT3(vn(i,:),cgdvs(n,m,m1)%m(:,i1))- &
                                   cgdvs(n,m,m1)%m(:,i1)*DOT3(vn(i,:),vn(i,:))- &
                                   2.0_r8_kind*vns(i1,:)*DOT3(vn(i,:),cgdv(n,m,m1)%m(:,i))
                              dtn=unit_vector_deriv(vvv,dvvv)
                              !--------------------------------------
                              vv=vector_product(vns(i,:),vn(i2,:))
                              dvv=vector_product(cgdvs(n,m,m1)%m(:,i),vn(i2,:))+vector_product(vns(i,:),cgdv(n,m,m1)%m(:,i2))
                              d_nos=unit_vector_deriv(vv,dvv); if(change_signs) d_no=-d_no
                              !--------------------------------------
                              d2Tcs(:)=0.0_r8_kind; d2nos=0.0_r8_kind
                              if(bs_cut) then
                                 rr2=r_sphere(nsp)**2-r_sphere(poli%n_sphere(i,2))**2

                                 f=(dist1**2+rr2)/(2.0_r8_kind*dist2**2)

                                 do k=1,3
                                    dvbuf(k)=cagr%dc(k,poli%n_sphere(i,2))%xyz_hess(j,l,l1,n,m,m1)- &
                                         cagr%dc(k,nsp)%xyz_hess(j,l,l1,n,m,m1)
                                    vbuf1(k)=cagr%dc(k,poli%n_sphere(i,2))%xyz_grad(n,m,m1)-cagr%dc(k,nsp)%xyz_grad(n,m,m1)
                                 end do

                                 d2Tcs(:)=d2Tcs(:)+f*dvbuf(:)

                                 drr2=r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,l1)-  &
                                      r_sphere(poli%n_sphere(i,2))*cagr%dR(poli%n_sphere(i,2))%xyz_grad(j,l,l1)
                                 drr21=r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(n,m,m1)-  &
                                      r_sphere(poli%n_sphere(i,2))*cagr%dR(poli%n_sphere(i,2))%xyz_grad(n,m,m1)

                                 df=(drr2-rr2*DOT3(diks,vbufs)/dist2**2)/dist2**2
                                 df1=(drr21-rr2*DOT3(diks,vbuf1)/dist2**2)/dist2**2

                                 d2Tcs(:)=d2Tcs(:)+df1*vbufs(:)+df*vbuf1(:)

                                 d2f=(cagr%dR(nsp)%xyz_grad(n,m,m1)*cagr%dR(nsp)%xyz_grad(j,l,l1)+ &
                                      r_sphere(nsp)*cagr%dR(nsp)%xyz_hess(j,l,l1,n,m,m1)- &
                                      cagr%dR(poli%n_sphere(i,2))%xyz_grad(n,m,m1)*cagr%dR(poli%n_sphere(i,2))%xyz_grad(j,l,l1)- &
                                      r_sphere(poli%n_sphere(i,2))*cagr%dR(poli%n_sphere(i,2))%xyz_hess(j,l,l1,n,m,m1)- &
                                      2.0_r8_kind*drr21*DOT3(diks,vbufs)/dist2**2- &
                                      rr2*DOT3(vbuf1,vbufs)/dist2**2- &
                                      rr2*DOT3(diks,dvbuf)/dist2**2+ &
                                      2.0_r8_kind*rr2*DOT3(diks,vbufs)*DOT3(diks,vbuf1)/dist2**4)/dist2**2- &
                                      2.0_r8_kind*df*DOT3(diks,vbuf1)/dist2**2


                                 d2Tcs(:)=d2Tcs(:)+diks(:)*d2f

                                 cgd2Ts(j,l,l1,n,m,m1)%m(:,i)=d2Tcs(:)
                                 !----------------------------------------
                                 d2nos=dvbuf(:)/dist2-vbufs*DOT3(diks,vbuf1)/dist2**3- &
                                      vbuf1*DOT3(diks,vbufs)/dist2**3- &
                                      diks*DOT3(vbuf1,vbufs)/dist2**3-diks*DOT3(diks,dvbuf)/dist2**3+ &
                                      3.0_r8_kind*diks*DOT3(diks,vbufs)*DOT3(diks,vbuf1)/dist2**5
                              end if

                              dmu=(DOT3(d_nos,dTcs_d(i,:,j,l,l1))+DOT3(nos,d2Tcs)-DOT3(d2nos,vns(i,:))- &
                                   DOT3(dnos,cgdvs(n,m,m1)%m(:,i))-DOT3(d_nos,ep)*cagr%dR(nsp)%xyz_grad(j,l,l1)- &
                                   DOT3(nos,dep)*cagr%dR(nsp)%xyz_grad(j,l,l1)- &
                                   DOT3(nos,ep)*cagr%dR(nsp)%xyz_hess(j,l,l1,n,m,m1))/dpnostn


                              d2P(:)=dep(:)*cagr%dR(nsp)%xyz_grad(j,l,l1)+ep(:)*cagr%dR(nsp)%xyz_hess(j,l,l1,n,m,m1) + &
                                   mus_d(i,j,l,l1)*dtns(:)+tns(i,:)*dmus- &
                                   tns(i,:)*mus_d(i,j,l,l1)*(DOT3(d_no,tns(i,:))+DOT3(no,dtns))/dpnotns+ &
                                   mu_d(i,j,l,l1)*dtn(:)+tn(i,:)*dmu- &
                                   tn(i,:)*mu_d(i,j,l,l1)*(DOT3(d_nos,tn(i,:))+DOT3(nos,dtn))/dpnostn

                              cgd2P(j,l,l1,n,m,m1)%m(:,i)=d2P(:)

!!$if(N_total+1 == 274)then
!!$   print*,i,l,l1,j,m,m1,d2P,b_cut,bs_cut
!!$endif
                              if(.not.extra_point(i)) then
                                 do k=1,3
                                    cagr%d2center(j,l,l1,n,m,m1)%m(k,N_total+1)=d2P(k)+ &
                                         cagr%d2center(j,l,l1,n,m,m1)%m(k,N_total+1)
                                 end do
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         enddo
      end if secnd_d

!print*,'n_vertises',poli%n_vertises
      do i=1,poli%n_vertises
!print*,'vertex- ',i,poli%xyz_vertex(i,:)
         i1=poli%bounds(i,1)
         i2=poli%bounds(i,2)
!print*,'1- ',poli%xyz_vertex(i1,:)
!print*,'2- ',poli%xyz_vertex(i2,:)

         do l=1,N_moving_unique_atoms+N_pc
!print*,'Atom- ',l
            if(l <= N_moving_unique_atoms) then
               na=moving_unique_atom_index(l)
               ea=unique_atoms(na)%n_equal_atoms
            else
               na=l-N_moving_unique_atoms
               if(with_pc) ea=pointcharge_array(na)%N_equal_charges
            end if
            do l1=1,ea
               do j=1,3

                  vc=cgdv(j,l,l1)%m(:,i)
                  vd=cgdvs(j,l,l1)%m(:,i1)

                  cgdphi(j,l,l1)%m(i)=  deriv_cos(vn(i,:),vns(i1,:),vc,vd,.true.,.false.)

                  if(weight_cent_mass) then
!!!! MF bug fix 12/2001
                         ! recalculate dTc
                         dTc(:)=cgdP(j,l,l1)%m(:,i)-cgdv(j,l,l1)%m(:,i)

                         dp1=phin(i)*0.5_r8_kind
                         tanphi2=tan(dp1)
                         help=DOT3(vn(i,:),vc)/poli%r_bound(i,1)
                                !derivative of r_bound(i,1)
                         ep=poli%xyz_bound(i,1:3)-xyz_sphere(nsp,:)
                         epk=vn(i,:)+vns(i1,:)
                         cagr%dcenter(j,l,l1)%m(:,N_total+1)= &
                             cagr%dcenter(j,l,l1)%m(:,N_total+1)+&
                             help*(phin(i)*ep+tanphi2*epk) + & !deriv of r_bound
                             cgdphi(j,l,l1)%m(i)*poli%r_bound(i,1)* &
                                (ep+0.5_r8_kind/(cos(dp1)**2)*epk) + &
!!!! MF bug fix 9/2001
!wrong:                      poli%r_bound(i,1)*(phin(i)*dTc(i)+help*(vc+vd))
                             poli%r_bound(i,1)*(phin(i)*dTc(:)+tanphi2*(vc+vd))
                  endif

                  ! expansion of the double vector product terms ax(bxc)=b(ac)-c(ab)
                  ! and setting (vv) = r_bound**2
                  ! derivative of not normalized tangent
                  dp1=DOT3(vn(i,:),vns(i1,:))
                  tng =vn(i,:)*(DOT3(vn(i,:),vd)+DOT3(vc,vns(i1,:))) &
                       + vc *  dp1                                   &
                       - vd * poli%r_bound(i,1)**2                  &
                       - 2.0_r8_kind * vns(i1,:) * DOT3(vc,vn(i,:))

                  ! norm of double vector product vx(vxv*)
                  help=poli%r_bound(i,1)**2 * &
                       (poli%r_bound(i,1)**4 - dp1**2)
                  if(help<0.0_r8_kind) help=0.0_r8_kind
                  help=sqrt(help)

                  tng=tng/help
                  tng=tng-tn(i,:)*DOT3(tn(i,:),tng)

                  vc=cgdvs(j,l,l1)%m(:,i)
                  vd=cgdv(j,l,l1)%m(:,i2)

                  dp1=DOT3(vns(i,:),vn(i2,:))

                  tnsg=vns(i,:)*(DOT3(vns(i,:),vd)+DOT3(vc,vn(i2,:))) &
                       + vc *  dp1                                  &
                       - vd * poli%r_bound(i,2)**2                  &
                       - 2.0_r8_kind *vn(i2,:)* DOT3(vc,vns(i,:))

                  help=poli%r_bound(i,2)**2 * &
                       (poli%r_bound(i,2)**4 - dp1**2)
                  if(help<0.0_r8_kind) help=0.0_r8_kind
                  help=sqrt(help)
                  tnsg=tnsg/help
                  tnsg=tnsg-tns(i,:)*DOT3(tns(i,:),tnsg)

                  cgdomeg(j,l,l1)%m(i)=-deriv_cos(tn(i,:),tns(i,:),tng,tnsg,.true.,extra_point(i))

                  if(i==1) cagr%darea(j,l,l1)%m(N_total+1) = 2.0_r8_kind * &
                       cagr%dR(nsp)%xyz_grad(j,l,l1) *area / r_sphere(nsp)
                  cagr%darea(j,l,l1)%m(N_total+1)= &
                       cagr%darea(j,l,l1)%m(N_total+1) +&
                       (phin(i)* cgdcosth(j,l,l1)%m(i)    + &
                       costhn(i)* cgdphi(j,l,l1)%m(i)     - &
                       cgdomeg(j,l,l1)%m(i)) *r_sphere(nsp)**2
                  if(weight_cent_sin) then
                      cagr%dcenter(j,l,l1)%m(:,N_total+1)= &
                             cagr%dcenter(j,l,l1)%m(:,N_total+1)-&
                             (poli%xyz_vertex(i,:)-xyz_sphere(nsp,:))*2.0_r8_kind*( &
                                DOT3(tn(i,:),tnsg(:))+ &
                                DOT3(tng(:),tns(i,:)))
                  endif
               enddo
            enddo
         enddo
      enddo

      if(integralpar_2dervs) then
         do i=1,poli%n_vertises
            i1=poli%bounds(i,1)
            i2=poli%bounds(i,2)

            b_cut=.false.
            if(poli%n_sphere(i,1)/=0) b_cut=.true.
            bs_cut=.false.
            if(poli%n_sphere(i,2)/=0) bs_cut=.true.

            do l=1,N_moving_unique_atoms+N_pc
               if(l <= N_moving_unique_atoms) then
                  na=moving_unique_atom_index(l)
                  ea=unique_atoms(na)%n_equal_atoms
               else
                  na=l-N_moving_unique_atoms
                  ea=pointcharge_array(na)%N_equal_charges
               end if
               do l1=1,ea
                  do j=1,3

                     do m=1,N_moving_unique_atoms+N_pc
                        if(m <= N_moving_unique_atoms) then
                           na1=moving_unique_atom_index(m)
                           ea1=unique_atoms(na1)%n_equal_atoms
                        else
                           na1=m-N_moving_unique_atoms
                           ea1=pointcharge_array(na1)%N_equal_charges
                        end if
                        do m1=1,ea1
                           do n=1,3
                              !----------------------------------------
                              if(b_cut) then
                                 v1(:)=poli%xyz_vertex(i,:)-xyz_sphere(nsp,:)
                                 dv11(:)=cgdP(j,l,l1)%m(:,i)
                                 dv12(:)=cgdP(n,m,m1)%m(:,i)
                                 d2v1(:)=cgd2P(j,l,l1,n,m,m1)%m(:,i)

                                 v2(:)=xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:)
                                 do k=1,3
                                    dv21(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,l1)-&
                                         cagr%dc(k,nsp)%xyz_grad(j,l,l1)
                                    dv22(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(n,m,m1)-&
                                         cagr%dc(k,nsp)%xyz_grad(n,m,m1)
                                    d2v2(k)=cagr%dc(k,poli%n_sphere(i,1))%xyz_hess(j,l,l1,n,m,m1)-&
                                         cagr%dc(k,nsp)%xyz_hess(j,l,l1,n,m,m1)
                                 end do
                                 d2costh=deriv2_cos_and_angle(v1,v2,dv11,dv12,dv21,dv22,d2v1,d2v2,.false.,.false.)
                              else
                                 d2costh=0.0_r8_kind
                              end if
                              !----------------------------------------
                              v1(:)=vn(i,:)
                              v2(:)=vns(i1,:)
                              dv11(:)=cgdv(j,l,l1)%m(:,i)
                              dv12(:)=cgdv(n,m,m1)%m(:,i)
                              dv21(:)=cgdvs(j,l,l1)%m(:,i1)
                              dv22(:)=cgdvs(n,m,m1)%m(:,i1)
                              d2v1(:)=cgd2P(j,l,l1,n,m,m1)%m(:,i)-cgd2T(j,l,l1,n,m,m1)%m(:,i)
                              d2v2(:)=cgd2P(j,l,l1,n,m,m1)%m(:,i1)-cgd2Ts(j,l,l1,n,m,m1)%m(:,i1)

                              d2vb1=d2v1; d2vb2=d2v2

                              d2phi=deriv2_cos_and_angle(v1,v2,dv11,dv12,dv21,dv22,d2v1,d2v2,.true.,extra_point(i))
                              !----------------------------------------
                              v1(:)=tn(i,:)

                              vvv=vector_product(vn(i,:),vector_product(vn(i,:),vns(i1,:)))
                              dvvv=cgdv(j,l,l1)%m(:,i)*DOT3(vn(i,:),vns(i1,:))+ &
                                   vn(i,:)*DOT3(cgdv(j,l,l1)%m(:,i),vns(i1,:))+ &
                                   vn(i,:)*DOT3(vn(i,:),cgdvs(j,l,l1)%m(:,i1))- &
                                   cgdvs(j,l,l1)%m(:,i1)*DOT3(vn(i,:),vn(i,:))- &
                                   2.0_r8_kind*vns(i1,:)*DOT3(vn(i,:),cgdv(j,l,l1)%m(:,i))
                              dv11(:)=unit_vector_deriv(vvv,dvvv)

                              dvvv1=cgdv(n,m,m1)%m(:,i)*DOT3(vn(i,:),vns(i1,:))+ &
                                   vn(i,:)*DOT3(cgdv(n,m,m1)%m(:,i),vns(i1,:))+ &
                                   vn(i,:)*DOT3(vn(i,:),cgdvs(n,m,m1)%m(:,i1))- &
                                   cgdvs(n,m,m1)%m(:,i1)*DOT3(vn(i,:),vn(i,:))- &
                                   2.0_r8_kind*vns(i1,:)*DOT3(vn(i,:),cgdv(n,m,m1)%m(:,i))
                              dv12(:)=unit_vector_deriv(vvv,dvvv1)

                              d2vvv=d2vb1(:)*DOT3(vn(i,:),vns(i1,:))+ &
                                   cgdv(j,l,l1)%m(:,i)*DOT3(cgdv(n,m,m1)%m(:,i),vns(i1,:))+ &
                                   cgdv(j,l,l1)%m(:,i)*DOT3(vn(i,:),cgdvs(n,m,m1)%m(:,i1))+ &
                                   cgdv(n,m,m1)%m(:,i)*DOT3(cgdv(j,l,l1)%m(:,i),vns(i1,:))+ &
                                   vn(i,:)*DOT3(d2vb1,vns(i1,:))+ &
                                   vn(i,:)*DOT3(cgdv(j,l,l1)%m(:,i),cgdvs(n,m,m1)%m(:,i1))+ &
                                   cgdv(n,m,m1)%m(:,i)*DOT3(vn(i,:),cgdvs(j,l,l1)%m(:,i1))+ &
                                   vn(i,:)*DOT3(cgdv(n,m,m1)%m(:,i),cgdvs(j,l,l1)%m(:,i1))+ &
                                   vn(i,:)*DOT3(vn(i,:),d2vb2)- &
                                   d2vb2*DOT3(vn(i,:),vn(i,:))- &
                                   2.0_r8_kind*cgdvs(j,l,l1)%m(:,i1)*DOT3(vn(i,:),cgdv(n,m,m1)%m(:,i))- &
                                   2.0_r8_kind*cgdvs(n,m,m1)%m(:,i1)*DOT3(vn(i,:),cgdv(j,l,l1)%m(:,i))- &
                                   2.0_r8_kind*vns(i1,:)*DOT3(cgdv(n,m,m1)%m(:,i),cgdv(j,l,l1)%m(:,i))- &
                                   2.0_r8_kind*vns(i1,:)*DOT3(vn(i,:),d2vb1)
                              d2v1(:)=unit_vector_2deriv(vvv,dvvv,dvvv1,d2vvv)
                              !.............................................
                              v2(:)=tns(i,:)

                              vvv=vector_product(vns(i,:),vector_product(vns(i,:),vn(i2,:)))
                              dvvv=cgdvs(j,l,l1)%m(:,i)*DOT3(vns(i,:),vn(i2,:))+ &
                                   vns(i,:)*DOT3(cgdvs(j,l,l1)%m(:,i),vn(i2,:))+ &
                                   vns(i,:)*DOT3(vns(i,:),cgdv(j,l,l1)%m(:,i2))- &
                                   cgdv(j,l,l1)%m(:,i2)*DOT3(vns(i,:),vns(i,:))- &
                                   2.0_r8_kind*vn(i2,:)*DOT3(vns(i,:),cgdvs(j,l,l1)%m(:,i))
                              dv21(:)=unit_vector_deriv(vvv,dvvv)

                              dvvv1=cgdvs(n,m,m1)%m(:,i)*DOT3(vns(i,:),vn(i2,:))+ &
                                   vns(i,:)*DOT3(cgdvs(n,m,m1)%m(:,i),vn(i2,:))+ &
                                   vns(i,:)*DOT3(vns(i,:),cgdv(n,m,m1)%m(:,i2))- &
                                   cgdv(n,m,m1)%m(:,i2)*DOT3(vns(i,:),vns(i,:))- &
                                   2.0_r8_kind*vn(i2,:)*DOT3(vns(i,:),cgdvs(n,m,m1)%m(:,i))
                              dv22(:)=unit_vector_deriv(vvv,dvvv1)

                              d2vb1(:)=cgd2P(j,l,l1,n,m,m1)%m(:,i)-cgd2Ts(j,l,l1,n,m,m1)%m(:,i)
                              d2vb2(:)=cgd2P(j,l,l1,n,m,m1)%m(:,i2)-cgd2T(j,l,l1,n,m,m1)%m(:,i2)

                              d2vvv=d2vb1(:)*DOT3(vns(i,:),vn(i2,:))+ &
                                   cgdvs(j,l,l1)%m(:,i)*DOT3(cgdvs(n,m,m1)%m(:,i),vn(i2,:))+ &
                                   cgdvs(j,l,l1)%m(:,i)*DOT3(vns(i,:),cgdv(n,m,m1)%m(:,i2))+ &
                                   cgdvs(n,m,m1)%m(:,i)*DOT3(cgdvs(j,l,l1)%m(:,i),vn(i2,:))+ &
                                   vns(i,:)*DOT3(d2vb1,vn(i2,:))+ &
                                   vns(i,:)*DOT3(cgdvs(j,l,l1)%m(:,i),cgdv(n,m,m1)%m(:,i2))+ &
                                   cgdvs(n,m,m1)%m(:,i)*DOT3(vns(i,:),cgdv(j,l,l1)%m(:,i2))+ &
                                   vns(i,:)*DOT3(cgdvs(n,m,m1)%m(:,i),cgdv(j,l,l1)%m(:,i2))+ &
                                   vns(i,:)*DOT3(vns(i,:),d2vb2)- &
                                   d2vb2*DOT3(vns(i,:),vns(i,:))- &
                                   2.0_r8_kind*cgdv(j,l,l1)%m(:,i2)*DOT3(vns(i,:),cgdvs(n,m,m1)%m(:,i))- &
                                   2.0_r8_kind*cgdv(n,m,m1)%m(:,i2)*DOT3(vns(i,:),cgdvs(j,l,l1)%m(:,i))- &
                                   2.0_r8_kind*vn(i2,:)*DOT3(cgdvs(n,m,m1)%m(:,i),cgdvs(j,l,l1)%m(:,i))- &
                                   2.0_r8_kind*vn(i2,:)*DOT3(vns(i,:),d2vb1)
                              d2v2(:)=unit_vector_2deriv(vvv,dvvv,dvvv1,d2vvv)

                              d2omeg=-deriv2_cos_and_angle(v1,v2,dv11,dv12,dv21,dv22,d2v1,d2v2,.true.,extra_point(i))
                              !----------------------------------------

                              cagr%d2area(j,l,l1,n,m,m1)%m(N_total+1)= &
                                   cagr%d2area(j,l,l1,n,m,m1)%m(N_total+1)+ &
                                   (phin(i)* cgdcosth(j,l,l1)%m(i)+ &
                                   costhn(i)* cgdphi(j,l,l1)%m(i)- &
                                   cgdomeg(j,l,l1)%m(i))*r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(n,m,m1)* &
                                   2.0_r8_kind

                              cagr%d2area(j,l,l1,n,m,m1)%m(N_total+1)= &
                                   cagr%d2area(j,l,l1,n,m,m1)%m(N_total+1)+ &
                                   (cgdphi(n,m,m1)%m(i)*cgdcosth(j,l,l1)%m(i)+ &
                                   phin(i)*d2costh+ &
                                   cgdcosth(n,m,m1)%m(i)*cgdphi(j,l,l1)%m(i)+ &
                                   costhn(i)*d2phi- &
                                   d2omeg)*r_sphere(nsp)**2
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end if !integralpar_2derivs

!     now terms because of the scaling of the center to the sphere surface
!
!     ep=0.0_r8_kind
!     help=1.0_r8_kind
!     do i=1,poli%n_vertises
!        if(.not. extra_point(i)) then
!           if(weight_cent_sin) help=1.0_r8_kind-(DOT3(tn(i,:),tns(i,:)))**2
!           ep=ep+(poli%xyz_vertex(i,:)-xyz_sphere(nsp,:))*help
!        endif
!     enddo
      !representative point before scaling to sphere surface
      ep=(xyz_cent-xyz_sphere(nsp,:))*cent_f

      help=LEN3(ep)
      e_cent=ep/help

      do l=1,N_moving_unique_atoms+N_pc
         if(l <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(l)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=l-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l1=1,ea
            do j=1,3
             if(.not.orig_cent) then
               do k=1,3
                  vd(k)=cagr%dcenter(j,l,l1)%m(k,N_total+1)
               enddo
               vc=vd/help-ep*DOT3(ep,vd)/help**3
               do k=1,3
                  cagr%dcenter(j,l,l1)%m(k,N_total+1)= vc(k)*r_sphere(nsp) + &
                       cagr%dc(k,nsp)%xyz_grad(j,l,l1) + &
                       ep(k)/help*cagr%dR(nsp)%xyz_grad(j,l,l1)
               enddo
             else if(orig_cent) then !orig_cent
               do k=1,3
                  cagr%dcenter(j,l,l1)%m(k,N_total+1) = &
                       cagr%dc(k,nsp)%xyz_grad(j,l,l1) + &
                       (xyz_cent(k) - xyz_sphere(nsp,k))* &
                       cagr%dR(nsp)%xyz_grad(j,l,l1)/r_sphere(nsp)
               enddo
             endif
!!$if(N_total+1 == 274) then
!!$   print*, l,l1,j,cagr%darea(j,l,l1)%m(N_total+1)
!!$endif
             if(integralpar_2dervs)then
                do m=1,N_moving_unique_atoms+N_pc
                   if(m <= N_moving_unique_atoms) then
                      na1=moving_unique_atom_index(m)
                      ea1=unique_atoms(na1)%n_equal_atoms
                   else
                      na1=m-N_moving_unique_atoms
                      ea1=pointcharge_array(na1)%N_equal_charges
                   end if
                   do m1=1,ea1
                      do i=1,3
!!$                         vd1=sum(cgdP(i,m,m1)%m(:,:),2)
                         vd1=0.0_r8_kind
                         do k=1,poli%n_vertises
                            if(extra_point(i)) cycle
                            vd1=vd1+cgdP(i,m,m1)%m(:,k)
                         end do
                         vc1=vd1/help-ep*DOT3(ep,vd1)/help**3
                         vdd=cagr%d2center(j,l,l1,i,m,m1)%m(:,N_total+1)
                         vcc=vdd/help-vd*DOT3(ep,vd1)/help**3-vd1*DOT3(ep,vd)/help**3- &
                              ep*DOT3(vd1,vd)/help**3-ep*DOT3(ep,vdd)/help**3+ &
                              3.0_r8_kind*ep*DOT3(ep,vd)*DOT3(ep,vd1)/help**5
                         do k=1,3
                            cagr%d2center(j,l,l1,i,m,m1)%m(k,N_total+1)= &
                                 cagr%dc(k,nsp)%xyz_hess(j,l,l1,i,m,m1)+ &
                                 e_cent(k)*cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)+ &
                                 vc(k)*cagr%dR(nsp)%xyz_grad(i,m,m1)+vc1(k)*cagr%dR(nsp)%xyz_grad(j,l,l1)+ &
                                 vcc(k)*r_sphere(nsp)
                         end do
                         cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)=cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)+ &
                              2.0_r8_kind*(cagr%darea(i,m,m1)%m(N_total+1)*cagr%dR(nsp)%xyz_grad(j,l,l1)/r_sphere(nsp)- &
                              area*cagr%dR(nsp)%xyz_grad(j,l,l1)*cagr%dR(nsp)%xyz_grad(i,m,m1)/r_sphere(nsp)**2+ &
                              area*cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)/r_sphere(nsp))
!!$if(N_total+1 == 274) then
!!$   print*, l,l1,j,m,m1,i,cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)
!!$endif
                      end do
                   end do
                end do
             end if
            end do
         enddo
      enddo
!!$if(N_total+1 == 274) stop

      do i=1,N_moving_unique_atoms+N_pc
         if(i <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(i)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=i-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l=1,ea
            do j=1,3
               deallocate(cgdP(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang1")
               deallocate(cgdv(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang2")
               deallocate(cgdvs(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang3")
               deallocate(cgdomeg(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang4")
               deallocate(cgdcosth(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang5")
               deallocate(cgdphi(j,i,l)%m,stat=status)
               if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang6")
               if(integralpar_2dervs)then
                  do k=1,N_moving_unique_atoms+N_pc
                     if(i <= N_moving_unique_atoms) then
                        na1=moving_unique_atom_index(k)
                        ea1=unique_atoms(na1)%n_equal_atoms
                     else
                        na1=k-N_moving_unique_atoms
                        ea1=pointcharge_array(na1)%N_equal_charges
                     end if
                     do m=1,ea1
                        do n=1,3
                           deallocate(cgd2T(j,i,l,n,k,m)%m,cgd2Ts(j,i,l,n,k,m)%m,stat=status)
                           ASSERT(status.eq.0)
                           deallocate(cgd2P(j,i,l,n,k,m)%m,stat=status)
                           ASSERT(status.eq.0)
                        end do
                     end do
                  end do
               end if
            enddo
         enddo
      enddo

      deallocate(cgdP,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang7")
      deallocate(cgdv,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang8")
      deallocate(cgdvs,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang9")
      deallocate(cgdomeg,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang10")
      deallocate(cgdcosth,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang11")
      deallocate(cgdphi,stat=status)
      if(status/=0) call error_handler("dealloc fail geom_grad_cutted_triang12")
      if(integralpar_2dervs)then
         deallocate(mu_d,mus_d,stat=status)
         ASSERT(status.eq.0)
         deallocate(dTc_d,dTcs_d,stat=status)
         ASSERT(status.eq.0)
         deallocate(cgd2T,cgd2Ts,cgd2P,stat=status)
         ASSERT(status.eq.0)
      end if


    end subroutine geom_grad_cutted_triang
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine geom_grad_orig_triang(nsp,area,xyz_cent)
      !gradients of non cutted surface triangles
      !** End of interface *****************************************
      integer(kind=i4_kind) :: l,j,k,l1,na,ea
      integer(kind=i4_kind) :: m,i,m1,na1,ea1 !,n
      integer(kind=i4_kind) , intent(in):: nsp
      real(kind=r8_kind) , intent(in):: area,xyz_cent(3)
      real(r8_kind) :: e_cent(3)
      integer :: N_pc

      if(VTN) return

      N_pc=0
!!$      if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N
      if(with_pc) N_pc=pointcharge_N

      e_cent=(xyz_cent(:)-xyz_sphere(nsp,:))/r_sphere(nsp)

      do l=1,N_moving_unique_atoms+N_pc
         if(l <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(l)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=l-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l1=1,ea
            do j=1,3
               cagr%darea(j,l,l1)%m(N_total+1) = 2.0_r8_kind * area * &
                    cagr%dR(nsp)%xyz_grad(j,l,l1) / r_sphere(nsp)
               do k=1,3
                    !also valid for weight_cent_sin and weight_cent_mass
                    !movement only by sphere movement
                    cagr%dcenter(j,l,l1)%m(k,N_total+1) = &
                       cagr%dc(k,nsp)%xyz_grad(j,l,l1) + &
                       e_cent(k)*cagr%dR(nsp)%xyz_grad(j,l,l1)
               enddo
               if(integralpar_2dervs)then
                  do m=1,N_moving_unique_atoms+N_pc
                     if(m <= N_moving_unique_atoms) then
                        na1=moving_unique_atom_index(m)
                        ea1=unique_atoms(na1)%n_equal_atoms
                     else
                        na1=m-N_moving_unique_atoms
                        ea1=pointcharge_array(na1)%N_equal_charges
                     end if
                     do m1=1,ea1
                        do i=1,3
                           cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)= &
                                (2.0_r8_kind*area/r_sphere(nsp))*(cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)+ &
                                cagr%dR(nsp)%xyz_grad(j,l,l1)*cagr%dR(nsp)%xyz_grad(i,m,m1)/r_sphere(nsp))
                           do k=1,3
                              cagr%d2center(j,l,l1,i,m,m1)%m(k,N_total+1)= &
                                   cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)*e_cent(k)+ &
                                   cagr%dc(k,nsp)%xyz_hess(j,l,l1,i,m,m1)
                           end do
                        end do
                     end do
                  end do
               end if
            end do
         enddo
      enddo

    end subroutine geom_grad_orig_triang
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine geom_grad_reduced_triang(nsp,area,xyz_cent,rho)
      ! Purpose : Calculation of gradients and hessians of
      !           cavity surface elements - TsLess partition
      !------------ Modules used ------------------- ---------------
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind) , intent(in):: nsp
      real(r8_kind) , intent(in):: area,xyz_cent(3),rho
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      integer(i4_kind) :: l,j,k,l1,na,ea
      integer(i4_kind) :: m,i,m1,na1,ea1
      integer(i4_kind) :: n,n1,n2,n3,n4
      real(r8_kind) :: e_cent(3),x,x1,ms_r1,s_r1(3),ms_r2,s_r2(3),sigma,sigma_tot
      real(r8_kind) :: dPdx,dxdr1(3),dSigma_dr1(3),ds_r1(3,3),d2s_r(3,3,3),d2xdr2(3,3)
      real(r8_kind) :: d2Pdx2,dxdr2(3),dSigma_dr2(3),ds_r2(3,3),d2Sigma_dr2(3,3)
      real(r8_kind) :: v1(3),v2(3),v3(3)
      integer :: N_pc
      !------------ Executable code --------------------------------

      N_pc=0
      if(with_pc) N_pc=pointcharge_N

      e_cent=(xyz_cent(:)-xyz_sphere(nsp,:))/r_sphere(nsp)

      do l=1,N_moving_unique_atoms+N_pc
         if(l <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(l)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=l-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l1=1,ea
            do j=1,3
               do k=1,3
                    cagr%dcenter(j,l,l1)%m(k,N_total+1) = &
                       cagr%dc(k,nsp)%xyz_grad(j,l,l1) + &
                       e_cent(k)*cagr%dR(nsp)%xyz_grad(j,l,l1)
               enddo

               d_area:do n=1,N_spheres
!                  if(zero_area(n)) cycle d_area
                  if(n == nsp) cycle d_area
                  s_r1=xyz_cent(:)-xyz_sphere(n,:)
                  ms_r1=LEN3(s_r1)

                  x=(ms_r1-r_sphere(n))/rho

                  if(x <= d) then
                     dSigma_dr1(j)=0.0_r8_kind
                  else if(x >= 1.0_r8_kind) then
                     dSigma_dr1(j)=0.0_r8_kind
                  else if(x > d .and. x < 1.0_r8_kind) then
                     dPdx=0.0_r8_kind
                     do n1=1,Nd
                        dPdx=dPdx+(n1-1)*a(n1)*x**(n1-2)
                     end do

                     do n1=1,3
                        ds_r1(j,n1)=cagr%dcenter(j,l,l1)%m(n1,N_total+1)- &
                                   cagr%dc(n1,n)%xyz_grad(j,l,l1)
                     end do

                     dxdr1(j)=(DOT3(s_r1(:),ds_r1(j,:))/ms_r1- &
                              cagr%dR(n)%xyz_grad(j,l,l1))/rho

                     dSigma_dr1(j)=dPdx*dxdr1(j)
                  end if

                  sigma_tot=1.0_r8_kind
                  a1:do n2=1,N_spheres
!                     if(zero_area(n2)) cycle a1
                     if(n2==nsp) cycle a1
                     if(n2==n) cycle a1
                     x1=LEN3(xyz_cent(:)-xyz_sphere(n2,:))
                     x1=(x1-r_sphere(n2))/rho
                     if(x1 <= d) then
                        sigma=0.0_r8_kind
                     elseif(x1 >= 1.0_r8_kind) then
                        sigma=1.0_r8_kind
                     elseif(x1 > d .and. x1 < 1.0_r8_kind) then
                        sigma=0.0_r8_kind
                        do n1=1,Nd
                           sigma=sigma+a(n1)*x1**(n1-1)
                        end do
                     end if
                     sigma_tot=sigma_tot*sigma
                  end do a1
                  cagr%darea(j,l,l1)%m(N_total+1)= &
                       cagr%darea(j,l,l1)%m(N_total+1)+dSigma_dr1(j)*sigma_tot
               end do d_area
               cagr%darea(j,l,l1)%m(N_total+1)=cagr%darea(j,l,l1)%m(N_total+1)*area

            end do
         enddo
      enddo

      if(integralpar_2dervs)then
         do l=1,N_moving_unique_atoms+N_pc
            if(l <= N_moving_unique_atoms) then
               na=moving_unique_atom_index(l)
               ea=unique_atoms(na)%n_equal_atoms
            else
               na=l-N_moving_unique_atoms
               if(with_pc) ea=pointcharge_array(na)%N_equal_charges
            end if
            do l1=1,ea
               do j=1,3
                  do m=1,N_moving_unique_atoms+N_pc
                     if(m <= N_moving_unique_atoms) then
                        na1=moving_unique_atom_index(m)
                        ea1=unique_atoms(na1)%n_equal_atoms
                     else
                        na1=m-N_moving_unique_atoms
                        ea1=pointcharge_array(na1)%N_equal_charges
                     end if
                     do m1=1,ea1
                        do i=1,3
                           do k=1,3
                              cagr%d2center(j,l,l1,i,m,m1)%m(k,N_total+1)= &
                                   cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)*e_cent(k)+ &
                                   cagr%dc(k,nsp)%xyz_hess(j,l,l1,i,m,m1)
                           end do

                           area1:do n=1,N_spheres
!                              if(zero_area(n)) cycle area1
                              if(n==nsp) cycle area1

                              s_r1=xyz_cent(:)-xyz_sphere(n,:)
                              ms_r1=LEN3(s_r1)

                              x=(ms_r1-r_sphere(n))/rho

                              if(x <= d) then
                                 dSigma_dr1(j)=0.0_r8_kind
                                 d2Sigma_dr2(j,i)=0.0_r8_kind
                              else if(x >= 1.0_r8_kind) then
                                 dSigma_dr1(j)=0.0_r8_kind
                                 d2Sigma_dr2(j,i)=0.0_r8_kind
                              else if(x > d .and. x < 1.0_r8_kind) then
                                 dPdx=0.0_r8_kind
                                 d2Pdx2=0.0_r8_kind
                                 do n1=1,Nd
                                    dPdx=dPdx+(n1-1)*a(n1)*x**(n1-2)
                                    d2Pdx2=d2Pdx2+(n1-1)*(n1-2)*a(n1)*x**(n1-3)
                                 end do

                                 do n1=1,3
                                    ds_r1(j,n1)=cagr%dcenter(j,l,l1)%m(n1,N_total+1)- &
                                         cagr%dc(n1,n)%xyz_grad(j,l,l1)
                                    ds_r2(i,n1)=cagr%dcenter(i,l,l1)%m(n1,N_total+1)- &
                                         cagr%dc(n1,n)%xyz_grad(i,m,m1)
                                    d2s_r(j,i,n1)=cagr%d2center(j,l,l1,i,m,m1)%m(n1,N_total+1)- &
                                         cagr%dc(n1,n)%xyz_hess(j,l,l1,i,m,m1)
                                 end do

                                 d2xdr2(j,i)=DOT3(ds_r1(j,:),ds_r2(i,:))/ms_r1
                                 d2xdr2(j,i)=d2xdr2(j,i)+DOT3(s_r1(:),d2s_r(:,j,i))/ms_r1

                                 v1=DOT3(s_r1(:),ds_r1(j,:))
                                 v2=DOT3(s_r1(:),ds_r2(i,:))
                                 d2xdr2(j,i)=d2xdr2(j,i)-v1(j)*v2(i)/ms_r1**3

                                 d2xdr2(j,i)=d2xdr2(j,i)-cagr%dR(n)%xyz_hess(j,l,l1,i,m,m1)
                                 d2xdr2(j,i)=d2xdr2(j,i)/rho

                                 dxdr1(j)=(DOT3(s_r1(:),ds_r1(j,:))/ms_r1- &
                                      cagr%dR(n)%xyz_grad(j,l,l1))/rho
                                 dxdr2(i)=(DOT3(s_r1(:),ds_r2(i,:))/ms_r1- &
                                      cagr%dR(n)%xyz_grad(i,m,m1))/rho

                                 dSigma_dr1(j)=dPdx*dxdr1(j)

                                 d2Sigma_dr2(j,i)=d2Pdx2*dxdr1(j)*dxdr2(i)+dPdx*d2xdr2(j,i)
                              end if

                              sigma_tot=1.0_r8_kind
                              area2:do n2=1,N_spheres
!                                 if(zero_area(n2)) cycle area2
                                 if(n2==nsp) cycle area2
                                 if(n2==n) cycle area2
                                 x1=LEN3(xyz_cent(:)-xyz_sphere(n2,:))
                                 x1=(x1-r_sphere(n2))/rho
                                 if(x1 <= d) then
                                    sigma=0.0_r8_kind
                                 elseif(x1 >= 1.0_r8_kind) then
                                    sigma=1.0_r8_kind
                                 elseif(x1 > d .and. x1 < 1.0_r8_kind) then
                                    sigma=0.0_r8_kind
                                    do n1=1,Nd
                                       sigma=sigma+a(n1)*x1**(n1-1)
                                    end do
                                 end if
                                 sigma_tot=sigma_tot*sigma
                              end do area2

                              cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)= &
                                   cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)+sigma_tot*d2Sigma_dr2(j,i)

                              v3=0.0_r8_kind
                              area3:do n3=1,N_spheres
!                                 if(zero_area(n3)) cycle area3
                                 if(n3==nsp) cycle area3
                                 if(n3==n) cycle area3

                                 s_r2=xyz_cent(:)-xyz_sphere(n3,:)
                                 ms_r2=LEN3(s_r2)

                                 x=(ms_r2-r_sphere(n3))/rho

                                 if(x <= d) then
                                    dSigma_dr2(i)=0.0_r8_kind
                                 else if(x >= 1.0_r8_kind) then
                                    dSigma_dr2(i)=0.0_r8_kind
                                 else if(x > d .and. x < 1.0_r8_kind) then
                                    dPdx=0.0_r8_kind
                                    do n1=1,Nd
                                       dPdx=dPdx+(n1-1)*a(n1)*x**(n1-2)
                                    end do

                                    do n1=1,3
                                       ds_r2(i,n1)=cagr%dcenter(i,l,l1)%m(n1,N_total+1)- &
                                            cagr%dc(n1,n3)%xyz_grad(i,l,l1)
                                    end do

                                    dxdr2(i)=(DOT3(s_r2(:),ds_r2(i,:))/ms_r2- &
                                         cagr%dR(n3)%xyz_grad(i,l,l1))/rho

                                    dSigma_dr2(i)=dPdx*dxdr2(i)
                                 end if

                                 sigma_tot=1.0_r8_kind
                                 area4:do n4=1,N_spheres
!                                    if(zero_area(n4)) cycle area4
                                    if(n4==nsp) cycle area4
                                    if(n4==n) cycle area4
                                    if(n4==n3) cycle area4
                                    x1=LEN3(xyz_cent(:)-xyz_sphere(n4,:))
                                    x1=(x1-r_sphere(n4))/rho
                                    if(x1 <= d) then
                                       sigma=0.0_r8_kind
                                    elseif(x1 >= 1.0_r8_kind) then
                                       sigma=1.0_r8_kind
                                    elseif(x1 > d .and. x1 < 1.0_r8_kind) then
                                       sigma=0.0_r8_kind
                                       do n1=1,Nd
                                          sigma=sigma+a(n1)*x1**(n1-1)
                                       end do
                                    end if
                                    sigma_tot=sigma_tot*sigma
                                 end do area4

                                 v3(i)=v3(i)+sigma_tot*dSigma_dr2(i)
                              end do area3

                              cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)= &
                                   cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)+dSigma_dr1(j)*v3(i)
                           end do area1
                        end do
                     end do
                  end do
               end do
            enddo
         enddo
      end if

    end subroutine geom_grad_reduced_triang
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine geom_grad_reduced_triang1(nsp,area,xyz_cent)
      ! Purpose : Calculation of gradients of
      !           cavity surface elements - FIXVPN partition
      !------------ Modules used ------------------- ---------------
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind) , intent(in):: nsp
      real(r8_kind) , intent(in):: area,xyz_cent(3)
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      integer(i4_kind) :: l,j,k,l1,na,ea
      integer(i4_kind) :: m,i,m1,na1,ea1
      integer(i4_kind) :: n,n1,n2,n3,n4
      real(r8_kind) :: e_cent(3)
      real(r8_kind) :: f1,f2,df1dx(3),df2dx(3),df1dy(3),df2dy(3),d2f1dxdy(3,3),d2f2dxdy(3,3)
      real(r8_kind) :: s1(3),rj(3),ri(3),r5(3),r15(3),nn(3),nn2,mm2
      real(r8_kind) :: rji(3),drji,r1i(3),dr1i,drjidx(3,3),drjidy(3,3),dr1idx(3,3),dr1idy(3,3)
      real(r8_kind) :: Rsj,Rsi,cosa,sina,cosb,sinb,aji,bji,cji,dji,eji,fji
      real(r8_kind) :: dnn2dx(3),dmm2dx(3),ds1dx(3,3),dridx(3,3),drjdx(3,3)
      real(r8_kind) :: dnn2dy(3),dmm2dy(3),ds1dy(3,3),dridy(3,3),drjdy(3,3)
      real(r8_kind) :: d2nn2dxdy(3,3),d2mm2dxdy(3,3),d2s1dxdy(3,3,3)
      real(r8_kind) :: d2ridxdy(3,3,3),d2rjdxdy(3,3,3),d2rjidxdy(3,3,3),d2r1idxdy(3,3,3)
      real(r8_kind) :: dr5dx(3,3),dRsidx(3),dRsjdx(3),rbuf,rbuf1,rbuf2
      real(r8_kind) :: dr5dy(3,3),dRsidy(3),dRsjdy(3),dr15dx(3,3),dr15dy(3,3),d2r15dxdy(3,3,3)
      real(r8_kind) :: d2r5dxdy(3,3,3),d2Rsidxdy(3,3),d2Rsjdxdy(3,3)
      real(r8_kind) :: dajidx(3),dbjidx(3),dcjidx(3),ddjidx(3),dejidx(3),dfjidx(3)
      real(r8_kind) :: dajidy(3),dbjidy(3),dcjidy(3),ddjidy(3),dejidy(3),dfjidy(3)
      real(r8_kind) :: d2ajidxdy(3,3),d2bjidxdy(3,3),d2cjidxdy(3,3),d2djidxdy(3,3),d2ejidxdy(3,3),d2fjidxdy(3,3)
      real(r8_kind) :: nn21,nn22,mm21,mm22
      real(r8_kind) :: sigma,sigma_tot,sigma_tot1,dsigmadx(3),dsigmady(3),d2sigmadxdy(3,3),sigma_buf
      integer :: N_pc
      !------------ Executable code --------------------------------

      N_pc=0
      if(with_pc) N_pc=pointcharge_N

      e_cent=(xyz_cent(:)-xyz_sphere(nsp,:))/r_sphere(nsp)

      mm21=mm11*mm11; mm22=mm12*mm12; nn21=nn11*nn11; nn22=nn12*nn12

      s1=xyz_cent
      rj=xyz_sphere(nsp,:)
      Rsj=r_sphere(nsp)

      do l=1,N_moving_unique_atoms+N_pc
         if(l <= N_moving_unique_atoms) then
            na=moving_unique_atom_index(l)
            ea=unique_atoms(na)%n_equal_atoms
         else
            na=l-N_moving_unique_atoms
            if(with_pc) ea=pointcharge_array(na)%N_equal_charges
         end if
         do l1=1,ea
            do j=1,3
               do k=1,3
                    cagr%dcenter(j,l,l1)%m(k,N_total+1) = &
                       cagr%dc(k,nsp)%xyz_grad(j,l,l1) + &
                       e_cent(k)*cagr%dR(nsp)%xyz_grad(j,l,l1)
               enddo

               cagr%darea(j,l,l1)%m(N_total+1)=0.0_r8_kind
               d_area:do n=1,N_spheres
!                  if(zero_area(n)) cycle d_area
                  if(n == nsp) cycle d_area

                  dRsjdx(j)=cagr%dR(nsp)%xyz_grad(j,l,l1)

                  Rsi=r_sphere(n)
                  dRsidx(j)=cagr%dR(n)%xyz_grad(j,l,l1)

                  ri=xyz_sphere(n,:)
                  rji=rj-ri
                  drji=LEN3(rji)

                  do n1=1,3
                     ds1dx(j,n1)=cagr%dcenter(j,l,l1)%m(n1,N_total+1)
                     drjdx(j,n1)=cagr%dc(n1,nsp)%xyz_grad(j,l,l1)
                     dridx(j,n1)=cagr%dc(n1,n)%xyz_grad(j,l,l1)
                  end do
                  drjidx(j,:)=drjdx(j,:)-dridx(j,:)

                  r5=ri+rji*Rsi/drji

                  rbuf=DOT3(rji,drjidx(j,:))
                  do n1=1,3
                     dr5dx(j,n1)=dridx(j,n1)+drjidx(j,n1)*Rsi/drji+ &
                               rji(j)*dRsidx(n1)/drji-rji(n1)*rbuf*Rsi/drji**3
                  end do
                  dr15dx(j,:)=ds1dx(j,:)-dr5dx(j,:)

                  nn=s1-r5; r15=s1-r5
                  nn2=DOT3(nn,nn)

                  r1i=s1-ri
                  dr1i=LEN3(r1i)
                  dr1idx(j,:)=ds1dx(j,:)-dridx(j,:)

                  aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                  bji=2.0_r8_kind*Rsj*drji
                  cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                  dji=sqrt(bji*bji-aji*aji)
                  eji=sqrt(bji*bji-cji*cji)
                  fji=1.0_r8_kind/(bji*bji)

                  dajidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+DOT3(rji,drjidx(j,:))- &
                                       DOT3(r1i,dr1idx(j,:)))
                  dbjidx(j)=2.0_r8_kind*(dRsjdx(j)*drji+ &
                                       Rsj*DOT3(rji,drjidx(j,:))/drji)
                  dcjidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+DOT3(rji,drjidx(j,:))- &
                                       Rsi*dRsidx(j))
                  dfjidx(j)=-2.0_r8_kind*dbjidx(j)/bji**3

                  cosb=aji/bji; sinb=dji/bji
                  cosa=cji/bji; sina=eji/bji

                  ddjidx(j)=(bji*dbjidx(j)-aji*dajidx(j))/dji
                  dejidx(j)=(bji*dbjidx(j)-cji*dcjidx(j))/eji

                  mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                  if(drji >= Rsj+Rsi) then
                     df1dx=0.0_r8_kind
                     f1=1.0_r8_kind
                  else if(dr1i <= Rsi) then
                     df1dx=0.0_r8_kind
                     f1=0.0_r8_kind
                  else if(mm2 > mm22) then
                     df1dx=0.0_r8_kind
                     f1=1.0_r8_kind
                  else if(mm2 <= mm21) then
                     df1dx=0.0_r8_kind
                     f1=0.0_r8_kind
                  else if(mm2 > mm21 .and. mm2 <= mm22) then
                     f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                        15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                         6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5

                     dmm2dx(j)=4.0_r8_kind*Rsj*dRsjdx(j)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)- &
                               2.0_r8_kind*Rsj*Rsj*(dajidx(j)*cji*fji+      &
                                                    aji*dcjidx(j)*fji+      &
                                                    aji*cji*dfjidx(j)+      &
                                                    ddjidx(j)*eji*fji+      &
                                                    dji*dejidx(j)*fji+      &
                                                    dji*eji*dfjidx(j))
                     df1dx(j)=(30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2-   &
                               60.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3+   &
                               30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4)*  &
                               dmm2dx(j)/(mm22-mm21)
                  end if

                  if(nn2 > nn22) then
                     df2dx=0.0_r8_kind
                     f2=1.0_r8_kind
                  else if(nn2 <= nn21) then
                     df2dx=0.0_r8_kind
                     f2=0.0_r8_kind
                  else if(nn2 > nn21 .and. nn2 <= nn22) then
                     f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                        15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                         6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5

                     dnn2dx(j)=2.0_r8_kind*(DOT3(r15,dr15dx(j,:)))
                     df2dx(j)=(30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2-   &
                               60.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3+   &
                               30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4)*  &
                               dnn2dx(j)/(nn22-nn21)
                  end if

                  dsigmadx(j)=df1dx(j)*f2+df2dx(j)*f1

                  sigma_tot=1.0_r8_kind
                  a1:do n2=1,N_spheres
!                     if(zero_area(n2)) cycle a1
                     if(n2==nsp) cycle a1
                     if(n2==n) cycle a1

                     Rsi=r_sphere(n2)

                     ri=xyz_sphere(n2,:)
                     rji=rj-ri
                     drji=LEN3(rji)

                     r5=ri+rji*Rsi/drji

                     nn=s1-r5
                     nn2=DOT3(nn,nn)

                     r1i=s1-ri
                     dr1i=LEN3(r1i)

                     aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                     bji=2.0_r8_kind*Rsj*drji
                     cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                     dji=sqrt(bji*bji-aji*aji)
                     eji=sqrt(bji*bji-cji*cji)
                     fji=1.0_r8_kind/(bji*bji)

                     cosb=aji/bji; sinb=dji/bji
                     cosa=cji/bji; sina=eji/bji

                     mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                     if(drji >= Rsj+Rsi) then
                        f1=1.0_r8_kind
                     else if(dr1i <= Rsi) then
                        f1=0.0_r8_kind
                     else if(mm2 > mm22) then
                        f1=1.0_r8_kind
                     else if(mm2 <= mm21) then
                        f1=0.0_r8_kind
                     else if(mm2 > mm21 .and. mm2 <= mm22) then
                        f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                           15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                            6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5
                     end if

                     if(nn2 > nn22) then
                        f2=1.0_r8_kind
                     else if(nn2 <= nn21) then
                        f2=0.0_r8_kind
                     else if(nn2 > nn21 .and. nn2 <= nn22) then
                        f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                           15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                            6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5
                     end if
                     sigma=f1*f2

                     sigma_tot=sigma_tot*sigma
                  end do a1
                  cagr%darea(j,l,l1)%m(N_total+1)= &
                       cagr%darea(j,l,l1)%m(N_total+1)+dsigmadx(j)*sigma_tot
               end do d_area
               cagr%darea(j,l,l1)%m(N_total+1)=cagr%darea(j,l,l1)%m(N_total+1)*area

            end do
         enddo
      enddo

      if(integralpar_2dervs)then
         do l=1,N_moving_unique_atoms+N_pc
            if(l <= N_moving_unique_atoms) then
               na=moving_unique_atom_index(l)
               ea=unique_atoms(na)%n_equal_atoms
            else
               na=l-N_moving_unique_atoms
               if(with_pc) ea=pointcharge_array(na)%N_equal_charges
            end if
            do l1=1,ea
               do j=1,3
                  do m=1,N_moving_unique_atoms+N_pc
                     if(m <= N_moving_unique_atoms) then
                        na1=moving_unique_atom_index(m)
                        ea1=unique_atoms(na1)%n_equal_atoms
                     else
                        na1=m-N_moving_unique_atoms
                        ea1=pointcharge_array(na1)%N_equal_charges
                     end if
                     do m1=1,ea1
                        do i=1,3
                           do k=1,3
                              cagr%d2center(j,l,l1,i,m,m1)%m(k,N_total+1)= &
                                   cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)*e_cent(k)+ &
                                   cagr%dc(k,nsp)%xyz_hess(j,l,l1,i,m,m1)
                           end do

                           cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)=0.0_r8_kind
                           area1:do n=1,N_spheres
!                              if(zero_area(n)) cycle area1
                              if(n==nsp) cycle area1

                              dRsjdx(j)=cagr%dR(nsp)%xyz_grad(j,l,l1)
                              dRsjdy(i)=cagr%dR(nsp)%xyz_grad(i,m,m1)

                              d2Rsjdxdy(j,i)=cagr%dR(nsp)%xyz_hess(j,l,l1,i,m,m1)

                              Rsi=r_sphere(n)
                              dRsidx(j)=cagr%dR(n)%xyz_grad(j,l,l1)
                              dRsidy(i)=cagr%dR(n)%xyz_grad(i,m,m1)

                              d2Rsidxdy(j,i)=cagr%dR(n)%xyz_hess(j,l,l1,i,m,m1)

                              ri=xyz_sphere(n,:)
                              rji=rj-ri
                              drji=LEN3(rji)

                              do n1=1,3
                                 ds1dx(j,n1)=cagr%dcenter(j,l,l1)%m(n1,N_total+1)
                                 drjdx(j,n1)=cagr%dc(n1,nsp)%xyz_grad(j,l,l1)
                                 dridx(j,n1)=cagr%dc(n1,n)%xyz_grad(j,l,l1)
                                 ds1dy(i,n1)=cagr%dcenter(i,m,m1)%m(n1,N_total+1)
                                 drjdy(i,n1)=cagr%dc(n1,nsp)%xyz_grad(i,m,m1)
                                 dridy(i,n1)=cagr%dc(n1,n)%xyz_grad(i,m,m1)
                              end do
                              drjidx(j,:)=drjdx(j,:)-dridx(j,:); drjidy(i,:)=drjdy(i,:)-dridy(i,:)

                              do n1=1,3
                                 d2s1dxdy(j,i,n1)=cagr%d2center(j,l,l1,i,m,m1)%m(n1,N_total+1)
                                 d2rjdxdy(j,i,n1)=cagr%dc(n1,nsp)%xyz_hess(j,l,l1,i,m,m1)
                                 d2ridxdy(j,i,n1)=cagr%dc(n1,n)%xyz_hess(j,l,l1,i,m,m1)
                              end do
                              d2rjidxdy(j,i,:)=d2rjdxdy(j,i,:)-d2ridxdy(j,i,:)

                              r5=ri+rji*Rsi/drji

                              rbuf1=DOT3(rji,drjidx(j,:))
                              do n1=1,3
                                 dr5dx(j,n1)=dridx(j,n1)+drjidx(j,n1)*Rsi/drji+ &
                                             rji(j)*dRsidx(n1)/drji-rji(n1)*rbuf1*Rsi/drji**3
                              end do
                              dr15dx(j,:)=ds1dx(j,:)-dr5dx(j,:)

                              rbuf2=DOT3(rji,drjidy(i,:))
                              do n1=1,3
                                 dr5dy(i,n1)=dridy(i,n1)+drjidy(i,n1)*Rsi/drji+ &
                                             rji(i)*dRsidy(n1)/drji-rji(n1)*rbuf2*Rsi/drji**3
                              end do
                              dr15dy(i,:)=ds1dy(i,:)-dr5dy(i,:)

                              do n1=1,3
                                 d2r5dxdy(j,i,n1)=d2ridxdy(j,i,n1)+                                          &
                                                  d2rjidxdy(j,i,n1)*Rsi/drji+                                &
                                                  drjidx(j,n1)*dRsidy(i)/drji-                               &
                                                  drjidx(j,n1)*Rsi*rbuf2/drji**3+                            &
                                                  drjidy(i,n1)*dRsidx(j)/drji+                               &
                                                  rji(j)*d2Rsidxdy(i,n1)/drji-                               &
                                                  rji(j)*dRsidx(n1)*rbuf2/drji**3-                           &
                                                  drjidy(i,n1)*Rsi*rbuf1/drji**3-                            &
                                                  rji(n1)*dRsidy(i)*rbuf1/drji**3-                           &
                                                  rji(n1)*Rsi*DOT3(drjidy(i,:),drjidx(j,:))/drji**3-  &
                                                  rji(n1)*Rsi*DOT3(rji,d2rjidxdy(j,i,:))/drji**3+     &
                                                  3.0_r8_kind*rji(n1)*Rsi*rbuf1*rbuf2/drji**5
                              end do
                              d2r15dxdy(j,i,:)=d2s1dxdy(j,i,:)-d2r5dxdy(j,i,:)

                              nn=s1-r5; r15=s1-r5
                              nn2=DOT3(nn,nn)

                              r1i=s1-ri
                              dr1i=LEN3(r1i)
                              dr1idx(j,:)=ds1dx(j,:)-dridx(j,:); dr1idy(i,:)=ds1dy(i,:)-dridy(i,:)
                              d2r1idxdy(j,i,:)=d2s1dxdy(j,i,:)-d2ridxdy(j,i,:)

                              aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                              bji=2.0_r8_kind*Rsj*drji
                              cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                              dji=sqrt(bji*bji-aji*aji)
                              eji=sqrt(bji*bji-cji*cji)
                              fji=1.0_r8_kind/(bji*bji)

                              dajidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+DOT3(rji,drjidx(j,:))- &
                                   DOT3(r1i,dr1idx(j,:)))
                              dbjidx(j)=2.0_r8_kind*(dRsjdx(j)*drji+ &
                                       Rsj*DOT3(rji,drjidx(j,:))/drji)
                              dcjidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+DOT3(rji,drjidx(j,:))- &
                                   Rsi*dRsidx(j))
                              dfjidx(j)=-2.0_r8_kind*dbjidx(j)/bji**3

                              dajidy(i)=2.0_r8_kind*(Rsj*dRsjdy(i)+DOT3(rji,drjidy(i,:))- &
                                   DOT3(r1i,dr1idy(i,:)))
                              dbjidy(i)=2.0_r8_kind*(dRsjdy(i)*drji+ &
                                       Rsj*DOT3(rji,drjidy(i,:))/drji)
                              dcjidy(i)=2.0_r8_kind*(Rsj*dRsjdy(i)+DOT3(rji,drjidy(i,:))- &
                                   Rsi*dRsidy(i))
                              dfjidy(i)=-2.0_r8_kind*dbjidy(i)/bji**3

                              d2ajidxdy(j,i)=2.0_r8_kind*(dRsjdy(i)*dRsjdx(j)+Rsj*d2Rsjdxdy(j,i)+ &
                                                          DOT3(drjidy(i,:),drjidx(j,:))+   &
                                                          DOT3(rji(:),d2rjidxdy(j,i,:))-   &
                                                          DOT3(dr1idy(i,:),dr1idx(j,:))-   &
                                                          DOT3(r1i(:),d2r1idxdy(j,i,:)))

                              d2bjidxdy(j,i)=2.0_r8_kind*(d2Rsjdxdy(j,i)*drji+                           &
                                                          dRsjdx(j)*DOT3(rji,drjidy(i,:))/drji+   &
                                                          dRsjdy(i)*DOT3(rji,drjidx(j,:))/drji+   &
                                                          Rsj*DOT3(drjidy(i,:),drjidx(j,:))/drji+ &
                                                          Rsj*DOT3(rji,d2rjidxdy(j,i,:))/drji-    &
                                                          Rsj*DOT3(rji,drjidx(j,:))*              &
                                                              DOT3(rji,drjidy(i,:))/drji**3)

                              d2cjidxdy(j,i)=2.0_r8_kind*(dRsjdy(i)*dRsjdx(j)+Rsj*d2Rsjdxdy(j,i)+ &
                                                          DOT3(drjidy(i,:),drjidx(j,:))+   &
                                                          DOT3(rji(:),d2rjidxdy(j,i,:))-   &
                                                          dRsidy(i)*dRsidx(j)-Rsi*d2Rsidxdy(j,i))

                              d2fjidxdy(j,i)=-2.0_r8_kind*d2bjidxdy(j,i)/bji**3+ &
                                              6.0_r8_kind*dbjidx(j)*dbjidy(i)/bji**4

                              cosb=aji/bji; sinb=dji/bji
                              cosa=cji/bji; sina=eji/bji

                              ddjidx(j)=(bji*dbjidx(j)-aji*dajidx(j))/dji
                              dejidx(j)=(bji*dbjidx(j)-cji*dcjidx(j))/eji
                              ddjidy(i)=(bji*dbjidy(i)-aji*dajidy(i))/dji
                              dejidy(i)=(bji*dbjidy(i)-cji*dcjidy(i))/eji

                              d2djidxdy(j,i)=(dbjidy(i)*dbjidx(j)+bji*d2bjidxdy(j,i)-      &
                                              dajidy(i)*dajidx(j)-aji*d2ajidxdy(j,i))/dji- &
                                              ddjidx(j)*ddjidy(i)/dji
                              d2ejidxdy(j,i)=(dbjidy(i)*dbjidx(j)+bji*d2bjidxdy(j,i)-      &
                                              dcjidy(i)*dcjidx(j)-aji*d2cjidxdy(j,i))/eji- &
                                              dejidx(j)*dejidy(i)/eji

                              mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                              if(drji >= Rsj+Rsi) then
                                 d2f1dxdy=0.0_r8_kind
                                 df1dx=0.0_r8_kind; df1dy=0.0_r8_kind
                                 f1=1.0_r8_kind
                              else if(dr1i <= Rsi) then
                                 d2f1dxdy=0.0_r8_kind
                                 df1dx=0.0_r8_kind; df1dy=0.0_r8_kind
                                 f1=0.0_r8_kind
                              else if(mm2 > mm22) then
                                 d2f1dxdy=0.0_r8_kind
                                 df1dx=0.0_r8_kind; df1dy=0.0_r8_kind
                                 f1=1.0_r8_kind
                              else if(mm2 <= mm21) then
                                 d2f1dxdy=0.0_r8_kind
                                 df1dx=0.0_r8_kind; df1dy=0.0_r8_kind
                                 f1=0.0_r8_kind
                              else if(mm2 > mm21 .and. mm2 <= mm22) then
                                 f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                                    15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                                     6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5

                                 dmm2dx(j)=4.0_r8_kind*Rsj*dRsjdx(j)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)- &
                                           2.0_r8_kind*Rsj*Rsj*(dajidx(j)*cji*fji+      &
                                                                aji*dcjidx(j)*fji+      &
                                                                aji*cji*dfjidx(j)+      &
                                                                ddjidx(j)*eji*fji+      &
                                                                dji*dejidx(j)*fji+      &
                                                                dji*eji*dfjidx(j))

                                 dmm2dy(i)=4.0_r8_kind*Rsj*dRsjdy(i)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)- &
                                           2.0_r8_kind*Rsj*Rsj*(dajidy(i)*cji*fji+      &
                                                                aji*dcjidy(i)*fji+      &
                                                                aji*cji*dfjidy(i)+      &
                                                                ddjidy(i)*eji*fji+      &
                                                                dji*dejidy(i)*fji+      &
                                                                dji*eji*dfjidy(i))

                                 d2mm2dxdy(j,i)=4.0_r8_kind*dRsjdy(i)*dRsjdx(j)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)+ &
                                                4.0_r8_kind*Rsj*d2Rsjdxdy(j,i)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)- &
                                                4.0_r8_kind*Rsj*dRsjdx(j)*(dajidy(i)*cji*fji+   &
                                                                           aji*dcjidy(i)*fji+   &
                                                                           aji*cji*dfjidy(i)+   &
                                                                           ddjidy(i)*eji*fji+   &
                                                                           dji*dejidy(i)*fji+   &
                                                                           dji*eji*dfjidy(i))-  &
                                                4.0_r8_kind*Rsj*dRsjdy(i)*(dajidx(j)*cji*fji+   &
                                                                           aji*dcjidx(j)*fji+   &
                                                                           aji*cji*dfjidx(j)+   &
                                                                           ddjidx(j)*eji*fji+   &
                                                                           dji*dejidx(j)*fji+   &
                                                                           dji*eji*dfjidx(j))-  &
                                                2.0_r8_kind*Rsj*Rsj*(d2ajidxdy(j,i)*cji*fji+    &
                                                                     dajidx(j)*dcjidy(i)*fji+   &
                                                                     dajidx(j)*cji*dfjidy(i)+   &
                                                                     dajidy(i)*dcjidx(j)*fji+   &
                                                                     aji*d2cjidxdy(j,i)*fji+    &
                                                                     aji*dcjidx(j)*dfjidy(i)+   &
                                                                     dajidy(i)*cji*dfjidx(j)+   &
                                                                     aji*dcjidy(i)*dfjidx(j)+   &
                                                                     aji*cji*d2fjidxdy(j,i)+    &
                                                                     d2djidxdy(j,i)*eji*fji+    &
                                                                     ddjidx(j)*dejidy(i)*fji+   &
                                                                     ddjidx(j)*eji*dfjidy(i)+   &
                                                                     ddjidy(i)*dejidx(j)*fji+   &
                                                                     dji*d2ejidxdy(j,i)*fji+    &
                                                                     dji*dejidx(j)*dfjidy(i)+   &
                                                                     ddjidy(i)*eji*dfjidx(j)+   &
                                                                     dji*dejidy(i)*dfjidx(j)+   &
                                                                     dji*eji*d2fjidxdy(j,i))

                                 df1dx(j)=(30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2-   &
                                           60.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3+   &
                                           30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4)*  &
                                          dmm2dx(j)/(mm22-mm21)
                                 df1dy(i)=(30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2-   &
                                           60.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3+   &
                                           30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4)*  &
                                          dmm2dy(i)/(mm22-mm21)
                                 d2f1dxdy(j,i)=(60.0_r8_kind*((mm2-mm21)/(mm22-mm21))-      &
                                               180.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2+   &
                                               120.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3)*  &
                                               dmm2dx(j)*dmm2dy(i)/(mm22-mm21)**2+          &
                                               (30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2-   &
                                                60.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3+   &
                                                30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4)*  &
                                               d2mm2dxdy(j,i)/(mm22-mm21)
                              end if

                              if(nn2 > nn22) then
                                 d2f2dxdy=0.0_r8_kind
                                 df2dx=0.0_r8_kind; df2dy=0.0_r8_kind
                                 f2=1.0_r8_kind
                              else if(nn2 <= nn21) then
                                 d2f2dxdy=0.0_r8_kind
                                 df2dx=0.0_r8_kind; df2dy=0.0_r8_kind
                                 f2=0.0_r8_kind
                              else if(nn2 > nn21 .and. nn2 <= nn22) then
                                 f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                                    15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                                     6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5

                                 dnn2dx(j)=2.0_r8_kind*(DOT3(r15,dr15dx(j,:)))
                                 dnn2dy(i)=2.0_r8_kind*(DOT3(r15,dr15dy(i,:)))

                                 df2dx(j)=(30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2-   &
                                           60.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3+   &
                                           30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4)*  &
                                           dnn2dx(j)/(nn22-nn21)
                                 df2dy(i)=(30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2-   &
                                           60.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3+   &
                                           30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4)*  &
                                           dnn2dy(i)/(nn22-nn21)

                                 d2nn2dxdy(j,i)=2.0_r8_kind*(DOT3(dr15dy(i,:),dr15dx(j,:))+ &
                                                             DOT3(r15(:),d2r15dxdy(j,i,:)))
                                 d2f2dxdy(j,i)=(60.0_r8_kind*((nn2-nn21)/(nn22-nn21))-      &
                                               180.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2+   &
                                               120.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3)*  &
                                                dnn2dx(j)*dnn2dy(i)/(nn22-nn21)**2+         &
                                               (30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2-   &
                                                60.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3+   &
                                                30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4)*  &
                                                d2nn2dxdy(j,i)/(nn22-nn21)
                              end if

                              dsigmadx(j)=df1dx(j)*f2+df2dx(j)*f1
                              d2sigmadxdy(j,i)=d2f1dxdy(j,i)*f2+df1dx(j)*df2dy(i)+df1dy(i)*df2dx(j)+d2f2dxdy(j,i)*f1

                              sigma_tot=1.0_r8_kind
                              area2:do n2=1,N_spheres
!                                 if(zero_area(n2)) cycle area2
                                 if(n2==nsp) cycle area2
                                 if(n2==n) cycle area2

                                 Rsi=r_sphere(n2)

                                 ri=xyz_sphere(n2,:)
                                 rji=rj-ri
                                 drji=LEN3(rji)

                                 r5=ri+rji*Rsi/drji

                                 nn=s1-r5
                                 nn2=DOT3(nn,nn)

                                 r1i=s1-ri
                                 dr1i=LEN3(r1i)

                                 aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                                 bji=2.0_r8_kind*Rsj*drji
                                 cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                                 dji=sqrt(bji*bji-aji*aji)
                                 eji=sqrt(bji*bji-cji*cji)
                                 fji=1.0_r8_kind/(bji*bji)

                                 cosb=aji/bji; sinb=dji/bji
                                 cosa=cji/bji; sina=eji/bji

                                 mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                                 if(drji >= Rsj+Rsi) then
                                    f1=1.0_r8_kind
                                 else if(dr1i <= Rsi) then
                                    f1=0.0_r8_kind
                                 else if(mm2 > mm22) then
                                    f1=1.0_r8_kind
                                 else if(mm2 <= mm21) then
                                    f1=0.0_r8_kind
                                 else if(mm2 > mm21 .and. mm2 <= mm22) then
                                    f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                                       15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                                        6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5
                                 end if

                                 if(nn2 > nn22) then
                                    f2=1.0_r8_kind
                                 else if(nn2 <= nn21) then
                                    f2=0.0_r8_kind
                                 else if(nn2 > nn21 .and. nn2 <= nn22) then
                                    f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                                       15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                                        6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5
                                 end if
                                 sigma=f1*f2

                                 sigma_tot=sigma_tot*sigma
                              end do area2

                              cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)=cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)+ &
                                   d2sigmadxdy(j,i)*sigma_tot

                              sigma_buf=0.0_r8_kind
                              area3:do n3=1,N_spheres
!                                 if(zero_area(n3)) cycle area3
                                 if(n3==nsp) cycle area3
                                 if(n3==n) cycle area3

                                 dRsjdy(i)=cagr%dR(nsp)%xyz_grad(i,m,m1)

                                 Rsi=r_sphere(n3)
                                 dRsidy(i)=cagr%dR(n3)%xyz_grad(i,m,m1)

                                 ri=xyz_sphere(n3,:)
                                 rji=rj-ri
                                 drji=LEN3(rji)

                                 do n1=1,3
                                    ds1dy(i,n1)=cagr%dcenter(i,m,m1)%m(n1,N_total+1)
                                    drjdy(i,n1)=cagr%dc(n1,nsp)%xyz_grad(i,m,m1)
                                    dridy(i,n1)=cagr%dc(n1,n3)%xyz_grad(i,m,m1)
                                 end do
                                 drjidy(i,:)=drjdy(i,:)-dridy(i,:)

                                 r5=ri+rji*Rsi/drji

                                 rbuf=DOT3(rji,drjidy(i,:))
                                 do n1=1,3
                                    dr5dy(i,n1)=dridy(i,n1)+drjidy(i,n1)*Rsi/drji+ &
                                              rji(i)*dRsidy(n1)/drji-rji(n1)*rbuf*Rsi/drji**3
                                 end do
                                 dr15dy(i,:)=ds1dy(i,:)-dr5dy(i,:)

                                 nn=s1-r5; r15=s1-r5
                                 nn2=DOT3(nn,nn)

                                 r1i=s1-ri
                                 dr1i=LEN3(r1i)
                                 dr1idy(i,:)=ds1dy(i,:)-dridy(i,:)

                                 aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                                 bji=2.0_r8_kind*Rsj*drji
                                 cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                                 dji=sqrt(bji*bji-aji*aji)
                                 eji=sqrt(bji*bji-cji*cji)
                                 fji=1.0_r8_kind/(bji*bji)

                                 dajidy(i)=2.0_r8_kind*(Rsj*dRsjdy(i)+DOT3(rji,drjidy(i,:))- &
                                                      DOT3(r1i,dr1idy(i,:)))
                                 dbjidy(i)=2.0_r8_kind*(dRsjdy(i)*drji+ &
                                                      Rsj*DOT3(rji,drjidy(i,:))/drji)
                                 dcjidy(i)=2.0_r8_kind*(Rsj*dRsjdy(i)+DOT3(rji,drjidy(i,:))- &
                                                      Rsi*dRsidy(i))
                                 dfjidy(i)=-2.0_r8_kind*dbjidy(i)/bji**3

                                 cosb=aji/bji; sinb=dji/bji
                                 cosa=cji/bji; sina=eji/bji

                                 ddjidy(i)=(bji*dbjidy(i)-aji*dajidy(i))/dji
                                 dejidy(i)=(bji*dbjidy(i)-cji*dcjidy(i))/eji

                                 mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                                 if(drji >= Rsj+Rsi) then
                                    df1dy=0.0_r8_kind
                                    f1=1.0_r8_kind
                                 else if(dr1i <= Rsi) then
                                    df1dy=0.0_r8_kind
                                    f1=0.0_r8_kind
                                 else if(mm2 > mm22) then
                                    df1dy=0.0_r8_kind
                                    f1=1.0_r8_kind
                                 else if(mm2 <= mm21) then
                                    df1dy=0.0_r8_kind
                                    f1=0.0_r8_kind
                                 else if(mm2 > mm21 .and. mm2 <= mm22) then
                                    f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                                       15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                                        6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5

                                    dmm2dy(i)=4.0_r8_kind*Rsj*dRsjdy(i)*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)- &
                                              2.0_r8_kind*Rsj*Rsj*(dajidy(i)*cji*fji+      &
                                                                   aji*dcjidy(i)*fji+      &
                                                                   aji*cji*dfjidy(i)+      &
                                                                   ddjidy(i)*eji*fji+      &
                                                                   dji*dejidy(i)*fji+      &
                                                                   dji*eji*dfjidy(i))
                                    df1dy(i)=(30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**2-   &
                                              60.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3+   &
                                              30.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4)*  &
                                              dmm2dy(i)/(mm22-mm21)
                                 end if

                                 if(nn2 > nn22) then
                                    df2dx=0.0_r8_kind
                                    f2=1.0_r8_kind
                                 else if(nn2 <= nn21) then
                                    df2dx=0.0_r8_kind
                                    f2=0.0_r8_kind
                                 else if(nn2 > nn21 .and. nn2 <= nn22) then
                                    f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                                       15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                                        6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5

                                    dnn2dy(i)=2.0_r8_kind*(DOT3(r15,dr15dy(i,:)))

                                    df2dy(i)=(30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**2-   &
                                              60.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3+   &
                                              30.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4)*  &
                                              dnn2dy(i)/(nn22-nn21)
                                 end if

                                 dsigmady(i)=df1dy(i)*f2+df2dy(i)*f1

                                 sigma_tot1=1.0_r8_kind
                                 area4:do n4=1,N_spheres
!                                    if(zero_area(n4)) cycle area4
                                    if(n4==nsp) cycle area4
                                    if(n4==n) cycle area4
                                    if(n4==n3) cycle area4

                                    Rsi=r_sphere(n4)

                                    ri=xyz_sphere(n4,:)
                                    rji=rj-ri
                                    drji=LEN3(rji)

                                    r5=ri+rji*Rsi/drji

                                    nn=s1-r5
                                    nn2=DOT3(nn,nn)

                                    r1i=s1-ri
                                    dr1i=LEN3(r1i)

                                    aji=Rsj*Rsj+drji*drji-dr1i*dr1i
                                    bji=2.0_r8_kind*Rsj*drji
                                    cji=Rsj*Rsj+drji*drji-Rsi*Rsi
                                    dji=sqrt(bji*bji-aji*aji)
                                    eji=sqrt(bji*bji-cji*cji)
                                    fji=1.0_r8_kind/(bji*bji)

                                    cosb=aji/bji; sinb=dji/bji
                                    cosa=cji/bji; sina=eji/bji

                                    mm2=2.0_r8_kind*Rsj*Rsj*(1.0_r8_kind-aji*cji*fji-dji*eji*fji)

                                    if(drji >= Rsj+Rsi) then
                                       f1=1.0_r8_kind
                                    else if(dr1i <= Rsi) then
                                       f1=0.0_r8_kind
                                    else if(mm2 > mm22) then
                                       f1=1.0_r8_kind
                                    else if(mm2 <= mm21) then
                                       f1=0.0_r8_kind
                                    else if(mm2 > mm21 .and. mm2 <= mm22) then
                                       f1=10.0_r8_kind*((mm2-mm21)/(mm22-mm21))**3- &
                                          15.0_r8_kind*((mm2-mm21)/(mm22-mm21))**4+ &
                                           6.0_r8_kind*((mm2-mm21)/(mm22-mm21))**5
                                    end if

                                    if(nn2 > nn22) then
                                       f2=1.0_r8_kind
                                    else if(nn2 <= nn21) then
                                       f2=0.0_r8_kind
                                    else if(nn2 > nn21 .and. nn2 <= nn22) then
                                       f2=10.0_r8_kind*((nn2-nn21)/(nn22-nn21))**3- &
                                          15.0_r8_kind*((nn2-nn21)/(nn22-nn21))**4+ &
                                           6.0_r8_kind*((nn2-nn21)/(nn22-nn21))**5
                                    end if
                                    sigma=f1*f2

                                    sigma_tot1=sigma_tot1*sigma
                                 end do area4

                                 sigma_buf=sigma_buf+dsigmady(i)*sigma_tot1

                              end do area3

                              cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)=cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)+ &
                                   dsigmadx(j)*sigma_buf
                           end do area1

                           cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)=cagr%d2area(j,l,l1,i,m,m1)%m(N_total+1)*area

                        end do
                     end do
                  end do
               end do
            enddo
         enddo
      end if

    end subroutine geom_grad_reduced_triang1
    !------------------------------------------------------------

    !------------------------------------------------------------
    function unit_vector_deriv(vect,vect_deriv)

      real(r8_kind) :: unit_vector_deriv(3)
      real(r8_kind), intent(in) :: vect(3),vect_deriv(3)
      real(r8_kind) :: dvect,vvd

      dvect=LEN3(vect)

      vvd=DOT3(vect,vect_deriv)

      unit_vector_deriv=vect_deriv/dvect-vvd*vect/dvect**3

    end function unit_vector_deriv
    !------------------------------------------------------------

    !------------------------------------------------------------
    function unit_vector_2deriv(vect,vect_deriv1,vect_deriv2,vect_2deriv)

      real(r8_kind) :: unit_vector_2deriv(3)
      real(r8_kind), intent(in) :: vect(3),vect_deriv1(3),vect_deriv2(3),vect_2deriv(3)
      real(r8_kind) :: dvect,vvd1,vvd2,vd1vd2,vv2d

      dvect=LEN3(vect)

      vvd1=DOT3(vect,vect_deriv1)
      vvd2=DOT3(vect,vect_deriv2)
      vd1vd2=DOT3(vect_deriv1,vect_deriv2)
      vv2d=DOT3(vect,vect_2deriv)

      unit_vector_2deriv=vect_2deriv/dvect-vect_deriv1*vvd2/dvect**3- &
           vect_deriv2*vvd1/dvect**3-vect*vd1vd2/dvect**3-vect*vv2d/dvect**3+ &
           3.0_r8_kind*vect*vvd1*vvd2/dvect**5

    end function unit_vector_2deriv
    !------------------------------------------------------------

    !------------------------------------------------------------
    function deriv_cos(a, b, da, db, do_acos, extra_p)
      implicit none
      real(r8_kind), intent(in) :: a(:), b(:), da(:), db(:) ! (3)
      real(r8_kind) :: deriv_cos
      ! *** end of interface ***

     real(kind=r8_kind) :: costh,absa,absb,vbuf(3),sinth
     logical, intent(in) :: do_acos, extra_p
     real(kind=r8_kind), parameter :: small=1.e-9
     ! if do_acos == .true. the derivative of the angle is computed
     ! else the derivative of the cosine of the angle is computed

     ! see formula (31) in Cossi 96

     absa= LEN3(a)
     absb= LEN3(b)
     costh= DOT3(a,b)
     if(absa>small .and. absb>small) then
        costh=costh/(absa*absb)
        if(do_acos) then
           if(costh**2<1.0_r8_kind) then
              sinth = sqrt(1.0_r8_kind-costh**2)
           else
              sinth = 0.0_r8_kind !sign(max(small,abs(sinth)),sinth) !0.0_r8_kind
           endif
           if(sinth<small .and. .not. extra_p) then
              ABORT("deriv_cos: angle to small")
           else if (sinth<small) then
              deriv_cos=0.0_r8_kind
              return
           endif
        endif
        vbuf(:)=a(:)-costh*absa/absb* b(:)
        deriv_cos=DOT3(vbuf,db)
        vbuf(:)=b(:)-costh*absb/absa* a(:)
        deriv_cos=deriv_cos+DOT3(vbuf,da)
        deriv_cos=deriv_cos/(absa*absb)
        if(do_acos) deriv_cos = - deriv_cos/sinth
     else
        deriv_cos=0.0_r8_kind
     endif

    end function deriv_cos
 !!!MF <<<<
    !------------------------------------------------------------

    !------------------------------------------------------------
    function deriv2_cos_and_angle(v1,v2,dv11,dv12,dv21,dv22,d2v1,d2v2,do_angle,extra_p) result(deriv2)

      real(r8_kind) :: deriv2
      real(r8_kind) :: v1(3),v2(3),dv11(3),dv12(3),dv21(3),dv22(3),d2v1(3),d2v2(3)
      logical :: do_angle,extra_p

      real(r8_kind) :: cos_a,sin_a,dv1,dv2,aaa,daaa2,vect1(3),vect2(3),dvect(3)
      real(r8_kind) :: d2cos_a,dcos1,dcos2,dsin2
      real(kind=r8_kind), parameter :: small=1.e-9

      dv1=LEN3(v1); dv2=LEN3(v2)
      if(dv1 > small .and. dv2 > small) then
         cos_a=DOT3(v1,v2)/(dv1*dv2)
         if(cos_a**2 < 1.0_r8_kind) then
            sin_a=sqrt(1.0_r8_kind-cos_a**2)
         else
            sin_a=0.0_r8_kind
         end if

         dcos1=deriv_cos(v1,v2,dv11,dv21,.false.,extra_p)
         dcos2=deriv_cos(v1,v2,dv12,dv22,.false.,extra_p)

         vect1(:)=v1(:)-cos_a*v2(:)*dv1/dv2
         aaa=DOT3(vect1,dv21)
         vect2(:)=v2(:)-cos_a*v1(:)*dv2/dv1
         aaa=aaa+DOT3(vect2,dv11)

         dvect(:)=dv12(:)-dcos2*v2(:)*dv1/dv2-cos_a*DOT3(v1,dv12)*v2/(dv1*dv2)- &
              cos_a*dv1*unit_vector_deriv(v2,dv22)
         daaa2=DOT3(dvect,dv21)+DOT3(vect1,d2v2)
         dvect(:)=dv22(:)-dcos2*v1(:)*dv2/dv1-cos_a*DOT3(v2,dv22)*v1/(dv1*dv2)- &
              cos_a*dv2*unit_vector_deriv(v1,dv12)
         daaa2=daaa2+DOT3(dvect,dv11)+DOT3(vect2,d2v1)

         d2cos_a=-DOT3(v1,dv12)*aaa/(dv1**3*dv2)-DOT3(v2,dv22)*aaa/(dv1*dv2**3)+ &
              daaa2/(dv1*dv2)

         if(do_angle) then
            if(sin_a < small .and. .not. extra_p) then
               call error_handler("deriv2_cos: angle to small")
            else if(sin_a < small) then
               deriv2=0.0_r8_kind
            else
               dsin2=-cos_a*dcos2/sin_a
               deriv2=-d2cos_a/sin_a+dsin2*dcos1/sin_a**2
            end if
         else
            deriv2=d2cos_a
         end if
      else
         deriv2=0.0_r8_kind
      end if

    end function deriv2_cos_and_angle
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine symm_sorted_centers(ret_stat,verbose)
      ! find and sort symmetry equivalent surface tessera,
      ! store in tessarea
      implicit none
      integer(kind=i4_kind), intent(out) :: ret_stat
      logical, intent(in)                :: verbose
      !** End of interface *****************************************

      integer(kind=i4_kind) :: n_equal,n_equal_check
      real(kind=r8_kind)    :: this(3),that(3),cartes(3),check(3)
      type(sub_group) :: local_groups
      type(group_coset) :: cosets
      type(symm_transformation_int), allocatable :: point_trafos_buf(:)
      logical, allocatable :: checked(:)
      type(cavity_data),allocatable :: tessarea_buf(:)

      real(kind=r8_kind),parameter :: small = 1.0e-7_r8_kind
      integer(kind=i4_kind) :: i,n,j,k,status,l

      integer(i4_kind), parameter :: max_num_el = 256
      integer(i4_kind)            :: local_group_el(max_num_el)
      integer(i4_kind)            :: uniq
      real(r8_kind)               :: diff,norma

# define VPRINT if(verbose) print *,

      FPP_TIMER_DECL(tot)
      FPP_TIMER_DECL(prep)
      FPP_TIMER_DECL(cmp)

      VPRINT 'symm_sorted_centers: entered'
      FPP_TIMER_START(tot)

      ret_stat = 0

      allocate(checked(N_total),point_trafos_buf(N_total),stat=status)
      ASSERT(status==0)
      checked=.false.

      allocate(tessarea_buf(N_total),stat=status)
      ASSERT(status==0)

      if(do_gradients) then
         DPRINT 'symm_sorted_centers: allocate(i_symm_sort)'
         allocate(i_symm_sort(N_total,group_num_el),stat=status)
         ASSERT(status==0)
      endif

      VPRINT 'total=',n_total
      uniq=0
      l=0
      i_N_total: do i=1,N_total
         if(checked(i)) cycle
            VPRINT 'i=',i
            FPP_TIMER_START(prep)
            ! hopefully one more Stuek surface:
            uniq = uniq + 1
            ! reorder coordinates of points as
            ! (x,y,z) --> (x,z,y) in order to comply with the
            ! convention for angular momentum l=1
            this(1) = xyz_tes_c(i,1)
            this(2) = xyz_tes_c(i,3)
            this(3) = xyz_tes_c(i,2)
!!$            call normalize(this)
            !
            ! determine local symmetry groups
            !
            ! now apply all symmetry operations to the position of the
            ! surface point
            norma = max(norm(this),1.0_r8_kind)
            n = 0
            do j=1,group_num_el
               that = MATMUL(ylm_trafos(1)%matrix(:,:,j),this)
               VPRINT '1) this=',this
               VPRINT '1) that=',that
!!$               diff = DOT3(that-this,that-this)
               diff  = norm(that-this)
               if (diff < small*norma ) then
                  VPRINT '1) that==this: diff=',diff,'norma=',norma
                  n = n+1
                  ASSERT(n<=max_num_el)
                  local_group_el(n) = j
               else
                  VPRINT '1) that/=this: diff=',diff,'norma=',norma
               endif
            enddo

            ! allocate group elements
            local_groups%num_el = n
            allocate(local_groups%elements(n))
            ASSERT(status==0)
            ! and fill up group elements
            local_groups%elements(:) = local_group_el(:n)

            !
            ! now determine symmetry equivalent atoms
            !

            call group_coset_decomp(n_equal,local_groups,&
                 cosets,point_trafos_buf(uniq)%matrix)

             VPRINT 'n_equal=',n_equal


            ! search of positions of equal atoms

            allocate(tessarea_buf(uniq)%xyz(n_equal,3), &
                 tessarea_buf(uniq)%sphere(n_equal), stat=status)
            ASSERT(status==0)

            FPP_TIMER_STOP(prep)
            FPP_TIMER_START(cmp)

            n_equal_check=0
            j_n_equal: do j=1,n_equal
               that = MATMUL(ylm_trafos(1)%matrix(:,:,cosets%elements(1,j)),this)
               cartes(1) = that(1)
               cartes(2) = that(3)
               cartes(3) = that(2)
               k_N_total: do k=1,N_total
                  if(.not.checked(k)) then
!!$               if (LEN3(cartes-xyz_tes_c(k,:)) <= small*10.0_r8_kind) then
!!$               if (DOT3(cartes-xyz_tes_c(k,:), cartes-xyz_tes_c(k,:)) <= small*10.0_r8_kind) then
                     check = xyz_tes_c(k,:)
                     VPRINT '2) that  =',cartes,'(',i,j,')'
                     VPRINT '2) check =',check,'(',k,')'
!!$                     diff = DOT3(cartes-check,cartes-check)
                     diff = norm(cartes-check)
                     if (diff < small*norma) then
                        VPRINT '2) that==check: diff=',diff,'norma=',norma
                        checked(k) = .true.
                        n_equal_check = n_equal_check + 1
                        l=l+1
                        tessarea_buf(uniq)%xyz(j,:)  = xyz_tes_c(k,:)
                        tessarea_buf(uniq)%area      = area_tes(k)
                        tessarea_buf(uniq)%r_tes     = r_tes(k)
                        tessarea_buf(uniq)%sphere(j) = sphere(k)
                        tessarea_buf(uniq)%n_equal   = n_equal
!!$                        tessarea_buf(uniq)%tessar_data(j)=data_tes(k)
                        tessarea_buf(uniq)%cut_off   = cuttt(k)

                        if(do_gradients) then
                           i_symm_sort(uniq,j)=k
                        endif
                     else
                        VPRINT '2) that/=check: diff=',diff,'norma=',norma
                     endif
                  endif
               enddo k_N_total
            enddo j_n_equal
            if (n_equal_check /= n_equal) then
                print *,'failure at tes/area/cut',i,area_tes(i),cuttt(i)
                print *,'n_equal=',n_equal
                print *,'n_equal_check=',n_equal_check
                print *,'count(.not.checked)=',count(.not.checked)
                ret_stat = 1
                call error_handler( &
!!$                print *,&
                 "The program has calculated more equivalent points than &
                 & it had been able to find on the cavity surface. Who fib from us? &
                 & I do not know")
                print *, 'retrying...'
            endif
            deallocate(local_groups%elements, stat=status)
            ASSERT(status==0)
            if(ret_stat/=0) goto 1001 ! cleanup and return error
      enddo i_N_total
      FPP_TIMER_STOP(cmp)

      ! set the global var:
      n_size = uniq
      print *,'solv::symm_sorted_centers: uniq=',n_size,' total=',n_total

      allocate(tessarea(n_size),stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: allocation TESSAREA is failed")

      do i=1,n_size
         tessarea(i)%n_equal=tessarea_buf(i)%n_equal
         tessarea(i)%area=tessarea_buf(i)%area
         tessarea(i)%r_tes=tessarea_buf(i)%r_tes
         tessarea(i)%cut_off=tessarea_buf(i)%cut_off

         allocate(tessarea(i)%xyz(tessarea(i)%n_equal,3), &
              tessarea(i)%sphere(tessarea(i)%n_equal), stat=status)
         ASSERT(status==0)

         tessarea(i)%xyz    = tessarea_buf(i)%xyz
         tessarea(i)%sphere = tessarea_buf(i)%sphere

!!$         deallocate(tessarea_buf(i)%xyz, &
!!$              tessarea_buf(i)%sphere, stat=status)
!!$         ASSERT(status==0)
      enddo

      if(output_cavity_data .and. .not.do_gradients) then
         write (output_unit,*) '--------------------------------------------------------------'
         write (output_unit,*) &
              'The symmetrization of surface points has been done succesfully, number of unique points:', &
              n_size
         if(output_cavity_long) then
           write(output_unit,*)'number n_equal                  coordinates(a.u.)               area'
           do i=1,n_size
            write(output_unit,'(1x,i4,2x,i4,3x,3(f13.9,1x),f14.9)') i,tessarea(i)%n_equal, &
                 tessarea(i)%xyz(1,:),tessarea(i)%area
           enddo
          endif
         write(output_unit,*) '---------------------------------------------------------------'
      endif

      do i=1,n_size
         if(tessarea(i)%area <= 0.0_r8_kind)  then
            call write_to_output_units("---------------------------------------------------------")
            call write_to_output_units("because  of  some  tessarea  have  very  small  size  the")
            call write_to_output_units("calculating   areas   of  such  surface   elements   give")
            call write_to_output_units("unstable  results. Sometimes  the  value of area  can  be")
            call write_to_output_units("equal ZERO or  even be NEGATIVE. To  continue  your  work")
            call write_to_output_units("You should  either change  SCALED-FACTOR  to use  another")
            call write_to_output_units("size  of cavity  or choose  another  POINT-FACTOR to  use")
            call write_to_output_units("different number  of surface point  or if is possible  of")
            call write_to_output_units("course change symmetry of system to take another type  of")
            call write_to_output_units("poligone inscribed  to each sphere.  Sorry  unconvinient.")
            call write_to_output_units("---------------------------------------------------------")
            call error_handler("symm_sorted_centers: bad area. See description and output file!")
         endif
      enddo

      deallocate(xyz_tes_c, &
           area_tes, &
           r_tes, &
           sphere, &
           cuttt, &
           data_tes, stat=status)
      if ( status /= 0) call error_handler( &
           "points_on_cavity_surface: deallocation xyz_tes_c, area_tes, r_tes are failed")

1001  continue ! exeption enters here:

      do i=1,uniq ! could be less than n_size
         deallocate(&
              & tessarea_buf(i)%xyz, &
              & tessarea_buf(i)%sphere, &
              & stat=status)
         ASSERT(status==0)
      enddo

      deallocate(tessarea_buf,stat=status)
      ASSERT(status==0)
      deallocate(checked,stat=status)
      ASSERT(status==0)

      if(do_gradients.and.ret_stat/=0) then
         deallocate(i_symm_sort, stat=status)
         ASSERT(status==0)
      endif

      FPP_TIMER_STOP(tot)
      FPP_TIMER_PRINT(tot)
      FPP_TIMER_PRINT(prep)
      FPP_TIMER_PRINT(cmp)

    end subroutine symm_sorted_centers
    !--------------------------------------------------------

    real(r8_kind) function norm(x)
      implicit none
      real(r8_kind), intent(in) :: x(:) ! (3)
      ! *** end of interface ***

      norm = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
    end function norm

  end subroutine points_on_cavity_surface
  !****************************************************

  !************************************************************
  subroutine dealloc_geom_deriv_part1(cagr)
    ! dealloc dR and dc of cagr
    type(geom_deriv), intent(inout) :: cagr
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i,alloc_stat,j

    DPRINT '_part1: entered, N_spheres=',N_spheres
    DPRINT '_part1: shape(cagr%dR)=',shape(cagr%dR)
    DPRINT '_part1: shape(cagr%dc)=',shape(cagr%dc)
    ASSERT(size(cagr%dR)==n_spheres)
    ASSERT(size(cagr%dc,2)==n_spheres)
    ASSERT(size(cagr%dc,1)==3)

    if(VTN) return

    do i=1,N_spheres !size(cagr%dR)
       deallocate(cagr%dR(i)%xyz_grad,&
            stat=alloc_stat)
       if ( alloc_stat /= 0) call error_handler( &
            "dealloc_geom_deriv_part: deallocation of xyz_grad is failed (dR)")
       if(integralpar_2dervs) then
          deallocate(cagr%dR(i)%xyz_hess,&
               stat=alloc_stat)
          if ( alloc_stat /= 0) call error_handler( &
               "dealloc_geom_deriv_part: deallocation of xyz_hess is failed (dR)")
       end if
       do j=1,3
          DPRINT '_part1:',i,j,associated(cagr%dc(j,i)%xyz_grad)
          deallocate(cagr%dc(j,i)%xyz_grad,&
               stat=alloc_stat)
          if ( alloc_stat /= 0) call error_handler( &
               "dealloc_geom_deriv_part: deallocation of xyz_grad is failed (dc)")
          if(integralpar_2dervs) then
             deallocate(cagr%dc(j,i)%xyz_hess,&
                  stat=alloc_stat)
             if ( alloc_stat /= 0) call error_handler( &
                  "dealloc_geom_deriv_part: deallocation of xyz_hess is failed (dc)")
          end if
       enddo
    enddo
    deallocate(cagr%dR,stat=alloc_stat)
    if ( alloc_stat /= 0) call error_handler( &
         "dealloc_geom_deriv_part: deallocation of grad_atomic_center failed (dR)")
    deallocate(cagr%dc,stat=alloc_stat)
    if ( alloc_stat /= 0) call error_handler( &
         "dealloc_geom_deriv_part: deallocation of grad_atomic_center failed (dc)")

  end subroutine dealloc_geom_deriv_part1
  !************************************************************

  !******************************************************
  subroutine dealloc_geom_deriv_part2
   !dealloc darea,dcenter of cagr
   !** End of interface *****************************************

    use unique_atom_module, only : N_moving_unique_atoms,unique_atoms,moving_unique_atom_index
    use pointcharge_module, only : pointcharge_array, pointcharge_N

    integer(kind=i4_kind) :: status,N_pc
    integer(kind=i4_kind) :: i,j,l,na,ea
    integer(kind=i4_kind) :: i1,j1,l1,na1,ea1

    if(VTN) goto 1

    N_pc=0
    if(with_pc) N_pc=pointcharge_N

    do i=1,N_moving_unique_atoms+N_pc
       if(i <= N_moving_unique_atoms) then
          na = moving_unique_atom_index(i)
          ea = unique_atoms(na)%n_equal_atoms
       else
          na=i-N_moving_unique_atoms
          if(with_pc) ea=pointcharge_array(na)%N_equal_charges
       end if
       do l=1,ea
          do j=1,3
             deallocate(cagr%dcenter(j,i,l)%m, &
                  cagr%darea(j,i,l)%m,stat=status)
             ASSERT(status.eq.0)
             if(integralpar_2dervs) then
                do i1=1,N_moving_unique_atoms+N_pc
                   if(i1 <= N_moving_unique_atoms) then
                      na1 = moving_unique_atom_index(i1)
                      ea1 = unique_atoms(na1)%n_equal_atoms
                   else
                      na1=i1-N_moving_unique_atoms
                      ea1=pointcharge_array(na1)%N_equal_charges
                   end if
                   do l1=1,ea1
                      do j1=1,3
                         deallocate(cagr%d2center(j,i,l,j1,i1,l1)%m, &
                              cagr%d2area(j,i,l,j1,i1,l1)%m,stat=status)
                         ASSERT(status.eq.0)
                      end do
                   end do
                end do
             end if
          enddo
       enddo
    enddo

    deallocate(cagr%dcenter,cagr%darea,stat=status)
    ASSERT(status.eq.0)
    if(integralpar_2dervs) then
       deallocate(cagr%d2center,cagr%d2area,stat=status)
       ASSERT(status.eq.0)
    end if

1   deallocate(i_symm_sort,stat=status)
    ASSERT(status.eq.0)

  end subroutine dealloc_geom_deriv_part2
  !*****************************************************

  !***************************************************
  subroutine send_receive_geom_grad(grad_index)
   !** End of interface *****************************************

    use group_module, only : group_num_el
    use unique_atom_module, only: N_moving_unique_atoms,unique_atoms,moving_unique_atom_index, &
         unique_atom_grad_info,N_unique_atoms
    use elec_static_field_module, only: surf_points_grad_index,surf_points_grad_info
    use pointcharge_module, only: pointcharge_N,pointcharge_array
#ifdef WITH_EFP
    use induced_dipoles_module, only: N_ipd, do_pol_pcm
#endif
    use comm_module
    use msgtag_module, only: msgtag_geom_grad
    implicit none
    integer(kind=i4_kind), intent(in),optional:: grad_index(N_moving_unique_atoms + 1)
    integer(kind=i4_kind) :: dim,dim1,status,info,N_pc
    integer(kind=i4_kind) :: i,j,l,na,ea,ma,grad_dim,k
    integer(kind=i4_kind) :: i1,j1,l1,na1,ea1,grad_dim1 !, ma1
    real(kind=r8_kind),pointer :: rotmat(:,:), rotmat1(:,:)
    type(arrmat2) :: buffer(3,3)
    type(arrmat1) :: buffer1(3,3)

    if(comm_i_am_master()) then
       if(.not. present(grad_index)) call error_handler(&
          "send_receive_geom_grad: grad_index missing when master calls")

       to_calc_grads%dielconst=dielectric_constant

       N_pc=0
       if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

       to_calc_grads%n_points=n_size

       allocate(to_calc_grads%n_equal(to_calc_grads%n_points), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%n_equal on master failed")
       do i=1,to_calc_grads%n_points
          to_calc_grads%n_equal(i)=tessarea(i)%n_equal
       enddo

       dim=maxval(to_calc_grads%n_equal)
       allocate(to_calc_grads%xyz(to_calc_grads%n_points,dim,3), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%xyz on master failed")

       allocate(to_calc_grads%s(to_calc_grads%n_points), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%s on master failed")

       do i=1,to_calc_grads%n_points
          if(.not. VTN) to_calc_grads%s(i)=tessarea(i)%area
          do j=1,to_calc_grads%n_equal(i)
             to_calc_grads%xyz(i,j,:)=tessarea(i)%xyz(j,:)
          enddo
       enddo

       if(VTN) then
          allocate(to_calc_grads%sphere(to_calc_grads%n_points,dim), stat=status)
          ASSERT(status==0)
          do i=1,to_calc_grads%n_points
             do j=1,to_calc_grads%n_equal(i)
                to_calc_grads%sphere(i,j)=tessarea(i)%sphere(j)
             end do
          end do
       end if

#if 0
   call fix_solvq()
#endif
       allocate(to_calc_grads%Q(to_calc_grads%n_points), &
            stat=status)
       ASSERT(status==0)
#ifdef WITH_EFP
       allocate(to_calc_grads%Q1(to_calc_grads%n_points), &
            stat=status)
       ASSERT(status==0)
#endif
       do i=1,to_calc_grads%n_points
          to_calc_grads%Q(i)=0.0_r8_kind
          if(N_unique_atoms > 0) to_calc_grads%Q(i)=to_calc_grads%Q(i)+Q_n(i)+Q_e(i)
#ifdef WITH_EFP
          to_calc_grads%Q1(i)=to_calc_grads%Q(i)
          if(efp .and. n_efp > 0) then
             if((pointcharge_N > 0) .and. &
                  N_unique_atoms == 0) then
                to_calc_grads%Q(i)=to_calc_grads%Q(i)+Q_mp(i)
                to_calc_grads%Q1(i)=to_calc_grads%Q1(i)+Q_mp(i)
             end if
             if(N_ipd > 0 .and. do_pol_pcm) then
                to_calc_grads%Q1(i)=to_calc_grads%Q1(i)+Q_id1(i)
                to_calc_grads%Q(i)=to_calc_grads%Q(i)+Q_id(i)
             end if
          end if
#endif
       enddo
       call dealloc_Q()

       allocate(to_calc_grads%i_symm_sort(N_total,group_num_el), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%i_symm_sort on master failed")
       to_calc_grads%i_symm_sort=i_symm_sort

       if(.not. VTN) then
       allocate(to_calc_grads%dxyz_totsyms(3,N_moving_unique_atoms+N_pc), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%dxyz_totsyms on master failed")

       allocate(to_calc_grads%ds_totsyms(3,N_moving_unique_atoms+N_pc), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%ds_totsyms on master failed")
       endif

       if(integralpar_2dervs)then
          allocate(to_calc_grads%d2xyz_totsyms(3,N_moving_unique_atoms+N_pc, &
               3,N_moving_unique_atoms+N_pc), stat=status)
          if ( status /= 0) call error_handler( &
               "send_receive_geom_grad:allocation to_calc_grads%d2xyz_totsyms on master failed")
          allocate(to_calc_grads%d2s_totsyms(3,N_moving_unique_atoms+N_pc, &
               3,N_moving_unique_atoms+N_pc), stat=status)
          if ( status /= 0) call error_handler( &
               "send_receive_geom_grad:allocation to_calc_grads%d2s_totsyms on master failed")
       end if

       if(.not. VTN) then
       dim1=N_spheres*N_centers_on_sphere
       do i=1,N_moving_unique_atoms+N_pc
          if(i <= N_moving_unique_atoms) then
             na = moving_unique_atom_index(i)
             ea = unique_atoms(na)%n_equal_atoms
             grad_dim = grad_index(i+1) - grad_index(i)
          else
             if(with_pc) then
                na=i-N_moving_unique_atoms
                ea=pointcharge_array(na)%N_equal_charges
                grad_dim = surf_points_grad_index(na+1)-surf_points_grad_index(na)
             end if
          end if
          do j=1,3
             allocate(to_calc_grads%dxyz_totsyms(j,i)%m(3,dim1), &
                  stat=status)
             if ( status /= 0) call error_handler( &
                  "send_receive_geom_grad:allocation to_calc_grads%dxyz_totsyms%m on master failed")
             to_calc_grads%dxyz_totsyms(j,i)%m=0.0_r8_kind

             allocate(to_calc_grads%ds_totsyms(j,i)%m(dim1), &
                  stat=status)
             if ( status /= 0) call error_handler( &
                  "send_receive_geom_grad:allocation to_calc_grads%ds_totsyms%m on master failed")
             to_calc_grads%ds_totsyms(j,i)%m=0.0_r8_kind

          enddo
          do l=1,ea
             if(i <= N_moving_unique_atoms) then
                rotmat=>unique_atom_grad_info(i)%m(:,:,l)
             else
                if(with_pc) rotmat=>surf_points_grad_info(na)%m(:,:,l)
             end if
             do j=1,grad_dim
                do ma=1,3
                   to_calc_grads%dxyz_totsyms(j,i)%m(:,:)=&
                        to_calc_grads%dxyz_totsyms(j,i)%m(:,:)+&
                        rotmat(j,ma)*cagr%dcenter(ma,i,l)%m(:,:)

                   to_calc_grads%ds_totsyms(j,i)%m(:)=&
                        to_calc_grads%ds_totsyms(j,i)%m(:)+&
                        rotmat(j,ma)*cagr%darea(ma,i,l)%m(:)
                enddo
             enddo
          enddo
       enddo
       endif

       if(integralpar_2dervs)then
          do i=1,N_moving_unique_atoms+N_pc
             if(i <= N_moving_unique_atoms) then
                na = moving_unique_atom_index(i)
                ea = unique_atoms(na)%n_equal_atoms
                grad_dim = grad_index(i+1) - grad_index(i)
             else
                na=i-N_moving_unique_atoms
                ea=pointcharge_array(na)%N_equal_charges
                grad_dim = surf_points_grad_index(na+1)-surf_points_grad_index(na)
             end if
             do i1=1,N_moving_unique_atoms+N_pc
                if(i1 <= N_moving_unique_atoms) then
                   na1 = moving_unique_atom_index(i1)
                   ea1 = unique_atoms(na1)%n_equal_atoms
                   grad_dim1 = grad_index(i1+1) - grad_index(i1)
                else
                   na1=i1-N_moving_unique_atoms
                   ea1=pointcharge_array(na1)%N_equal_charges
                   grad_dim1 = surf_points_grad_index(na1+1)-surf_points_grad_index(na1)
                end if
                do j=1,3
                   do j1=1,3
                      allocate(to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(3,dim1), &
                           buffer(j,j1)%m(3,dim1),stat=status)
                      if ( status /= 0) call error_handler( &
                           "send_receive_geom_grad:allocation to_calc_grads%d2xyz_totsyms%m on master failed")
                      to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m=0.0_r8_kind
                      buffer(j,j1)%m=0.0_r8_kind

                      allocate(to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(dim1), &
                           buffer1(j,j1)%m(dim1),stat=status)
                      if ( status /= 0) call error_handler( &
                           "send_receive_geom_grad:allocation to_calc_grads%d2s_totsyms%m on master failed")
                      to_calc_grads%d2s_totsyms(j,i,j1,i1)%m=0.0_r8_kind
                      buffer1(j,j1)%m=0.0_r8_kind
                   end do
                end do

                do l=1,ea
                   if(i <= N_moving_unique_atoms) then
                      rotmat=>unique_atom_grad_info(i)%m(:,:,l)
                   else
                      rotmat=>surf_points_grad_info(na)%m(:,:,l)
                   end if
                   do l1=1,ea1
                      if(i1 <= N_moving_unique_atoms) then
                         rotmat1=>unique_atom_grad_info(i1)%m(:,:,l1)
                      else
                         rotmat1=>surf_points_grad_info(na1)%m(:,:,l1)
                      end if

                      do j=1,3
                         do j1=1,3
                            buffer(j,j1)%m=0.0_r8_kind
                            buffer1(j,j1)%m=0.0_r8_kind
                         end do
                      end do

                      do j=1,grad_dim
                         do ma=1,3
                            do k=1,3
                               buffer(j,k)%m(:,:)=buffer(j,k)%m(:,:)+&
                                    rotmat(j,ma)*cagr%d2center(ma,i,l,k,i1,l1)%m(:,:)
                               buffer1(j,k)%m(:)=buffer1(j,k)%m(:)+&
                                    rotmat(j,ma)*cagr%d2area(ma,i,l,k,i1,l1)%m(:)
                            end do
                         end do
                      end do

                      do j=1,grad_dim
                         do j1=1,grad_dim1
                            do k=1,3
                               to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(:,:)= &
                                    to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(:,:)+ &
                                    rotmat1(j1,k)*buffer(j,k)%m(:,:)
                               to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(:)= &
                                    to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(:)+ &
                                    rotmat1(j1,k)*buffer1(j,k)%m(:)
                            end do
                         end do
                      end do
                   end do
                end do

                do j=1,3
                   do j1=1,3
                      deallocate(buffer(j,j1)%m,buffer1(j,j1)%m,stat=status)
                      if ( status /= 0) call error_handler( &
                           "send_receive_geom_grad:deallocation buffer%m on master failed")
                   end do
                end do
             end do
          end do
       end if

       call dealloc_geom_deriv_part2()
       call dealloc_cavity()

    endif

    if(.not.comm_parallel()) RETURN
    ! RETURN POINT IF NOT PARALLEL

       if(comm_i_am_master()) then
          ! === context: master
          call comm_init_send(comm_all_other_hosts,msgtag_geom_grad)

          call commpack(VTN,info)
          ASSERT(info==0)

          call commpack(to_calc_grads%dielconst,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:dielconst pack failed")

          call commpack(with_pc,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:with_pc pack failed")
          call commpack(fixed_pc,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:fixed_pc pack failed")
          call commpack(to_calc_grads%n_points,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:n_points pack failed")
          call commpack(to_calc_grads%n_equal(1),to_calc_grads%n_points,1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:n_equal pack failed")

          call commpack(dim,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:dim pack failed")
          call commpack(to_calc_grads%xyz(1,1,1),to_calc_grads%n_points*dim*3,1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:n_equal xyz failed")

          if(.not. VTN) then
             call commpack(to_calc_grads%s(1),to_calc_grads%n_points,1,info)
             if(info/=0) call error_handler &
                  ("send_receive_geom_grad:n_equal s failed")
          end if

          if(VTN) then
             call commpack(to_calc_grads%sphere(1,1),to_calc_grads%n_points*dim,1,info)
             ASSERT(info==0)
          end if

          call commpack(to_calc_grads%Q(1),to_calc_grads%n_points,1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:n_equal Q failed")
#ifdef WITH_EFP
          call commpack(to_calc_grads%Q1(1),to_calc_grads%n_points,1,info)
          ASSERT(info==0)
#endif
          call commpack(N_total,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:N_total pack failed")
          call commpack(group_num_el,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:group_num_el pack failed")
          call commpack(to_calc_grads%i_symm_sort(1,1),N_total*group_num_el,1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:i_symm_sort pack failed")

          call commpack(ua_dim_max,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:ua_dim_max pack failed")

          if(.not. VTN) then
          call commpack(dim1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:dim1 pack failed")
          do i=1,N_moving_unique_atoms+N_pc
!!$             na = moving_unique_atom_index(i)
             do j=1,3
!!!?dim1
                   call commpack(to_calc_grads%dxyz_totsyms(j,i)%m(1,1),3*dim1,1,info)
                   if(info/=0) call error_handler &
                        ("send_receive_geom_grad:dxyz pack failed")
             enddo
          enddo

          do i=1,N_moving_unique_atoms+N_pc
             do j=1,3
                call commpack(to_calc_grads%ds_totsyms(j,i)%m(1),dim1,1,info)
                if(info/=0) call error_handler &
                     ("send_receive_geom_grad:ds pack failed")
             enddo
          enddo
          endif

          if(integralpar_2dervs)then
             do i=1,N_moving_unique_atoms+N_pc
                do j=1,3
                   do i1=1,N_moving_unique_atoms+N_pc
                      do j1=1,3
                         call commpack(to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(1,1),3*dim1,1,info)
                         if(info/=0) call error_handler &
                              ("send_receive_geom_grad:d2xyz pack failed")
                      end do
                   end do
                end do
             end do
             do i=1,N_moving_unique_atoms+N_pc
                do j=1,3
                   do i1=1,N_moving_unique_atoms+N_pc
                      do j1=1,3
                         call commpack(to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(1),dim1,1,info)
                         if(info/=0) call error_handler &
                              ("send_receive_geom_grad:d2s pack failed")
                      end do
                   end do
                end do
             end do
          end if

          if(VTN) then
             call commpack(center2sphere(1,1),(N_unique_atoms+pointcharge_N)*ua_dim_max,1,info)
             ASSERT(info==0)
          end if

          call comm_send()

       else
       ! === context: slaves
       call comm_save_recv(comm_master_host,msgtag_geom_grad)

       call communpack(VTN,info)
       ASSERT(info==0)

       call communpack(to_calc_grads%dielconst,info)
       if(info/=0) call error_handler &
               ("send_receive_geom_grad:dielconst unpack failed")

       call communpack(with_pc,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:with_pc unpack failed")

       call communpack(fixed_pc,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:fixed_pc unpack failed")

       N_pc=0
       if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

       call communpack(to_calc_grads%n_points,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:n_points unpack failed")

       allocate(to_calc_grads%n_equal(to_calc_grads%n_points), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%n_equal on slave failed")
       call communpack(to_calc_grads%n_equal(1),to_calc_grads%n_points,1,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:n_equal unpack failed")

       call communpack(dim,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:dim unpack failed")

       allocate(to_calc_grads%xyz(to_calc_grads%n_points,dim,3), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%xyz on slave failed")
       call communpack(to_calc_grads%xyz(1,1,1),to_calc_grads%n_points*dim*3,1,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:xyz unpack failed")

       if(.not. VTN) then
          allocate(to_calc_grads%s(to_calc_grads%n_points), &
               stat=status)
          if ( status /= 0) call error_handler( &
               "send_receive_geom_grad:allocation to_calc_grads%s on slave failed")
          call communpack(to_calc_grads%s(1),to_calc_grads%n_points,1,info)
          if(info/=0) call error_handler &
               ("send_receive_geom_grad:s unpack failed")
       end if

       if(VTN) then
          allocate(to_calc_grads%sphere(to_calc_grads%n_points,dim), stat=status)
          ASSERT(status==0)
          call communpack(to_calc_grads%sphere(1,1),to_calc_grads%n_points*dim,1,info)
          ASSERT(info==0)
       end if

       allocate(to_calc_grads%Q(to_calc_grads%n_points), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%Q on slave failed")
       call communpack(to_calc_grads%Q(1),to_calc_grads%n_points,1,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:Q unpack failed")
#ifdef WITH_EFP
       allocate(to_calc_grads%Q1(to_calc_grads%n_points), &
            stat=status)
       ASSERT(status==0)
       call communpack(to_calc_grads%Q1(1),to_calc_grads%n_points,1,info)
       ASSERT(info==0)
#endif
       call dealloc_Q()

       call communpack(N_total,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:N_total unpack failed")

       call communpack(group_num_el,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:group_num_el unpack failed")

       allocate(to_calc_grads%i_symm_sort(N_total,group_num_el), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%i_symm_sort on slave failed")
       call communpack(to_calc_grads%i_symm_sort(1,1),N_total*group_num_el,1,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:i_symm_sort unpack failed")

       call communpack(ua_dim_max,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:ua_dim_max unpack failed")

       if(.not. VTN) then
       call communpack(dim1,info)
       if(info/=0) call error_handler &
            ("send_receive_geom_grad:dim1 unpack failed")

       allocate(to_calc_grads%dxyz_totsyms(3,N_moving_unique_atoms+N_pc), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%dxyz_totsyms on slave failed")

       allocate(to_calc_grads%ds_totsyms(3,N_moving_unique_atoms+N_pc), &
            stat=status)
       if ( status /= 0) call error_handler( &
            "send_receive_geom_grad:allocation to_calc_grads%ds_totsyms on slave failed")
       endif

       if(integralpar_2dervs)then
          allocate(to_calc_grads%d2xyz_totsyms(3,N_moving_unique_atoms+N_pc, &
               3,N_moving_unique_atoms+N_pc),stat=status)
          ASSERT(status.eq.0)
          allocate(to_calc_grads%d2s_totsyms(3,N_moving_unique_atoms+N_pc, &
               3,N_moving_unique_atoms+N_pc),stat=status)
          ASSERT(status.eq.0)
       end if

       if(.not. VTN) then
       do i=1,N_moving_unique_atoms+N_pc
!!$          na = moving_unique_atom_index(i)
          do j=1,3
                allocate(to_calc_grads%dxyz_totsyms(j,i)%m(3,dim1), &
                     stat=status)
                if ( status /= 0) call error_handler( &
                     "send_receive_geom_grad:allocation to_calc_grads%dxyz_totsyms%m on slave failed")
                call communpack(to_calc_grads%dxyz_totsyms(j,i)%m(1,1),3*dim1,1,info)
                if(info/=0) call error_handler &
                     ("send_receive_geom_grad:dxyz_totsyms unpack failed")

          enddo
       enddo

       do i=1,N_moving_unique_atoms+N_pc
          do j=1,3
             allocate(to_calc_grads%ds_totsyms(j,i)%m(dim1), &
                  stat=status)
             if ( status /= 0) call error_handler( &
                  "send_receive_geom_grad:allocation to_calc_grads%ds_totsyms%m on slave failed")
             call communpack(to_calc_grads%ds_totsyms(j,i)%m(1),dim1,1,info)
             if(info/=0) call error_handler &
                  ("send_receive_geom_grad:ds_totsyms unpack failed")
          enddo
       enddo
       endif

       if(integralpar_2dervs)then
          do i=1,N_moving_unique_atoms+N_pc
             do j=1,3
                do i1=1,N_moving_unique_atoms+N_pc
                   do j1=1,3
                      allocate(to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(3,dim1), &
                           stat=status)
                      if ( status /= 0) call error_handler( &
                           "send_receive_geom_grad:allocation to_calc_grads%d2xyz_totsyms%m on slave failed")
                      call communpack(to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m(1,1),3*dim1,1,info)
                      if(info/=0) call error_handler &
                           ("send_receive_geom_grad:d2xyz unpack failed")
                   end do
                end do
             end do
          end do
          do i=1,N_moving_unique_atoms+N_pc
             do j=1,3
                do i1=1,N_moving_unique_atoms+N_pc
                   do j1=1,3
                      allocate(to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(dim1), &
                           stat=status)
                      if ( status /= 0) call error_handler( &
                           "send_receive_geom_grad:allocation to_calc_grads%d2s_totsyms%m on slave failed")
                      call communpack(to_calc_grads%d2s_totsyms(j,i,j1,i1)%m(1),dim1,1,info)
                      if(info/=0) call error_handler &
                           ("send_receive_geom_grad:d2s unpack failed")
                   end do
                end do
             end do
          end do
       end if

       if(VTN) then
          allocate(center2sphere(N_unique_atoms+pointcharge_N,ua_dim_max), stat=status)
          ASSERT(status==0)
          call communpack(center2sphere(1,1),(N_unique_atoms+pointcharge_N)*ua_dim_max,1,info)
          ASSERT(info==0)
       end if
    endif
#if 0
  contains ! of send_receive_geom_grad

   subroutine fix_solvq()
   integer(kind=i4_kind):: fix_q_unit,i,l
   logical:: fix_q

     print*, 'solv q writen sum(Q_n)',sum(Q_n)
               fix_q_unit=openget_iounit(file=trim(inpfile('fix_solvq.dat')), &
                               form='FORMATTED',status='unknown')

           write(fix_q_unit,*) Q_n,Q_e
           call returnclose_iounit(fix_q_unit,status='keep')
   end subroutine fix_solvq
#endif

  end subroutine send_receive_geom_grad
  !*****************************************************

  !******************************************************
  subroutine dealloc_Q
   !** End of interface *****************************************
    use solv_charge_mixing_module, only: mix_charges

    integer(kind=i4_kind) :: status

    if(allocated(Q_n)) then
       deallocate(Q_n,stat=status)
       if(status /= 0) call error_handler("dealloc_Q:deallocation Q_n failed")
    endif
    if(allocated(Q_e)) then
       deallocate(Q_e,stat=status)
       if(status /= 0) call error_handler("dealloc_Q:deallocation Q_e failed")
       if(mix_charges()) then
          deallocate(Q_e_old,stat=status)
          ASSERT(status==0)
       end if
    endif
#ifdef WITH_EFP
    if(allocated(Q_mp)) then
       deallocate(Q_mp,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_id)) then
       deallocate(Q_id,stat=status)
       ASSERT(status == 0)
    endif
    if(allocated(Q_id1)) then
       deallocate(Q_id1,stat=status)
       ASSERT(status == 0)
    endif
#endif

  end subroutine dealloc_Q
  !******************************************************

  !****************************************************
  subroutine dealloc_geom_grad
   !deallocates to_calc_grads
   !** End of interface *****************************************

    use unique_atom_module, only: N_moving_unique_atoms
    use pointcharge_module, only: pointcharge_N
    use comm_module
    use commpack_module

    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: i,i1,j,j1,N_pc !l,ea,

    N_pc=0
    if(with_pc .and. .not.fixed_pc) N_pc=pointcharge_N

    if(.not. VTN) then
       do i=1,N_moving_unique_atoms+N_pc
          do j=1,3
             deallocate(to_calc_grads%dxyz_totsyms(j,i)%m, &
                  stat=status)
             if ( status /= 0) call error_handler( &
                  "dealloc_geom_grad:deallocation to_calc_grads%dxyz_totsyms%m failed")

             deallocate(to_calc_grads%ds_totsyms(j,i)%m, &
                  stat=status)
             if ( status /= 0) call error_handler( &
                  "dealloc_geom_grad:deallocation to_calc_grads%ds_totsyms%m failed")

             if(integralpar_2dervs)then
                do i1=1,N_moving_unique_atoms+N_pc
                   do j1=1,3
                      deallocate(to_calc_grads%d2xyz_totsyms(j,i,j1,i1)%m, &
                           stat=status)
                      if ( status /= 0) call error_handler( &
                           "dealloc_geom_grad:deallocation to_calc_grads%d2xyz_totsyms%m failed")
                   end do
                end do
             end if
          enddo
       enddo
    end if

    deallocate(to_calc_grads%n_equal, &
         to_calc_grads%xyz, &
         to_calc_grads%Q, &
         to_calc_grads%i_symm_sort, &
         stat=status)
    if ( status /= 0) call error_handler( &
            "dealloc_geom_grad:deallocation to_calc_grad failed")
    if(.not. VTN) then
       deallocate(to_calc_grads%dxyz_totsyms, &
            to_calc_grads%ds_totsyms, &
            stat=status)
       ASSERT(status==0)
    end if
    if(VTN) then
       deallocate(to_calc_grads%sphere,stat=status)
       ASSERT(status==0)
       deallocate(center2sphere,stat=status)
       ASSERT(status==0)
    end if

    if(integralpar_2dervs)then
       deallocate(to_calc_grads%d2xyz_totsyms,stat=status)
       if ( status /= 0) call error_handler( &
            "dealloc_geom_grad:deallocation to_calc_grad%d2xyz_totsyms failed")
       deallocate(to_calc_grads%d2s_totsyms,stat=status)
       if ( status /= 0) call error_handler( &
            "dealloc_geom_grad:deallocation to_calc_grad%d2s_totsyms failed")
    end if

  end subroutine dealloc_geom_grad
  !****************************************************

  !******************************************************
  subroutine dealloc_cavity
   !deallocate tesserea
   !** End of interface *****************************************

    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: i

    DPRINT 'dealloc_cavity: enetered'
    DPRINT 'dealloc_cavity: size(tessarea)=',size(tessarea)
    ASSERT(n_size==size(tessarea))
    do i=1,size(tessarea)
       deallocate(tessarea(i)%xyz, &
            tessarea(i)%sphere, stat=status)
       if (status .ne. 0 ) call error_handler( &
            "dealloc_cavity: deallocation of TESSAREA%XYZ  failed")
    enddo
    deallocate(tessarea,stat=status)
    if ( status /= 0) call error_handler( &
         "dealloc_cavity: deallocation TESSAREA is failed")

    DPRINT 'dealloc_cavity: shape(r_sphere)=', shape(r_sphere)
    deallocate(r_sphere, stat=status)
    if ( status /= 0) call error_handler( &
         "dealloc_cavity: deallocation r_sphere is failed")
    if(allocated(iuniq))  then
        deallocate(iuniq,stat=status)
       if ( status /= 0) call error_handler( &
            "dealloc_cavity: deallocation of iuniq is failed")
    endif

    DPRINT 'dealloc_cavity: shape(xyz_sphere)=',shape(xyz_sphere)
    deallocate(xyz_sphere, stat=status)
    if ( status /= 0) call error_handler( &
         "dealloc_cavity: deallocation xyz_sphere is failed")

    if(allocated(parents)) then
       deallocate(parents, stat=status)
       if ( status /= 0) call error_handler( &
            "dealloc_cavity: deallocation of parents is failed")
    endif

    if(allocated(zero_area)) then
       deallocate(zero_area, stat=status)
       if ( status /= 0) call error_handler( &
            "dealloc_cavity: deallocation of zero_area is failed")
    endif

  end subroutine dealloc_cavity
  !******************************************************
end module solv_cavity_module


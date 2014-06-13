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
! Public interface of module
!===============================================================
module cavity_module

# include <def.h>
  use type_module ! type specification parameters
  use datatype
  use common_data_module
  use species_module
  use energy_and_forces_module
  use filename_module, only: inpfile
  use iounitadmin_module, only: output_unit,write_to_output_units,openget_iounit, &
       returnclose_iounit
  use inp_out_module
  use tasks_main_options_module
  use slab_module, only: slab_calc

  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private
!------------ public functions and subroutines ------------------
  public read_solv,write_solv,points_on_cavity_surface
  public dealloc_cavity_mm
  public dealloc_geom_deriv_part1,dealloc_geom_deriv_part2

!================================================================
  !if I will ever calculate this:
  logical, public :: electr_energy
  logical, public :: cavitation_energy
  logical, public :: disp_rep_energy

  !if I will calculate this now:
  logical , public :: do_electr, do_cavitation,do_disp_rep

!== Interrupt of public interface of module =========
!----private types
  type, private :: triangles
     real(kind=r8_kind),pointer :: xyz(:,:)
     real(kind=r8_kind),pointer :: xyz_centers(:,:)
     integer(kind=i4_kind),pointer :: index(:,:)
     real(kind=r8_kind) :: radius
     real(kind=r8_kind) :: area
  end type triangles

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

  type, public :: cavity_data
     integer (i4_kind) :: n_equal
     real(kind=r8_kind),pointer :: xyz(:,:)
     real(kind=r8_kind) :: area
     real(kind=r8_kind) :: r_tes
     integer(kind=i4_kind),pointer :: sphere(:)
     logical :: cut_off
  end type cavity_data

  type, public :: grad_atomic_center
     real(kind=r8_kind), pointer :: xyz_grad(:,:,:)
  end type grad_atomic_center

  type, public :: geom_deriv
     type(grad_atomic_center), pointer :: dc(:,:)
     type(grad_atomic_center), pointer :: dR(:)
     type(arrmat1), pointer :: darea(:,:,:)
     type(arrmat2), pointer :: dcenter(:,:,:)
  end type geom_deriv

  integer(kind=i4_kind) :: N_total
  integer(kind=i4_kind),public :: n_size

  real(kind=r8_kind), public, allocatable :: Q_at(:)

  type(cavity_data),public,allocatable :: tessarea(:)
  integer(kind=i4_kind),public :: N_spheres
  integer(kind=i4_kind) :: N_atom_spheres
  real(kind=r8_kind), public, allocatable :: r_sphere(:)
  real(kind=r8_kind), public, allocatable :: xyz_sphere(:,:)
  integer(kind=i4_kind), allocatable :: parents(:,:)
  logical, allocatable :: zero_area(:)
  integer(kind=i4_kind) :: N_points_of_triangles
  integer(kind=i4_kind) :: N_centers_on_sphere
  real(kind=r8_kind),public :: R_access
  real(kind=r8_kind) :: m_t_a

  integer(kind=i4_kind), allocatable :: at_index(:)

  type(geom_deriv),public :: cagr

  logical :: fix_number_add = .false.
  real(kind=r8_kind) :: scaled_factor=1.2_r8_kind

  real(kind=r8_kind),public :: solvent_radius
  real(kind=r8_kind),public :: dielectric_constant
  real(kind=r8_kind),public :: abs_temperature
  real(kind=r8_kind),public :: solvent_volume
  logical :: no_hydrogen_sphere
  logical, public :: hydrogen_no_scale
  real(kind=r8_kind) :: fradio,overlap_angle,rmin_gepol,min_area
  integer :: gepol          !choose gepol algorithm (Gepol87 or Gepol93)
  logical :: save_cav !output cavity and surface points into VRML format
  real(r8_kind) :: overlap_factor !new overlap factor for Gepol93
  integer,public :: output_level
  logical :: smoothing
  logical,public :: use_ewald

  namelist /solv_options/ &
       smoothing, &
       electr_energy, &
       cavitation_energy, &
       disp_rep_energy, &
       dielectric_constant, &
       abs_temperature, &
       solvent_volume, &
       solvent_radius, &
       no_hydrogen_sphere, &
       hydrogen_no_scale, &
       fradio, &
       overlap_angle, &
       rmin_gepol, &
       min_area, &
       gepol, &
       save_cav, &
       overlap_factor, &
       output_level, &
       use_ewald

  logical :: df_smoothing = .false.
  logical :: df_electr_energy = .true.
  logical :: df_cavitation_energy = .true.
  logical :: df_disp_rep_energy = .true.
  real(kind=r8_kind) :: df_dielectric_constant = eps_h2o
  real(kind=r8_kind) :: df_solvent_radius = 1.400_r8_kind   !angstrom
  real(kind=r8_kind) :: df_abs_temperature = 298.0_r8_kind    !K
  real(kind=r8_kind) :: df_solvent_volume = 18.07_r8_kind   !cm^3/mol
  logical :: df_no_hydrogen_sphere = .false.
  logical :: df_hydrogen_no_scale = .true.
  real(kind=r8_kind) :: df_fradio = 0.7_r8_kind
  real(kind=r8_kind) :: df_overlap_angle = 40.0_r8_kind
  real(kind=r8_kind) :: df_rmin_gepol = 0.2_r8_kind
  real(kind=r8_kind) :: df_min_area = 1.0e-7_r8_kind
  integer :: df_gepol=93
  logical :: df_save_cav = .false.
  real(r8_kind) :: df_overlap_factor=0.77_r8_kind
  integer :: df_output_level=0
  logical :: df_use_ewald = .false.

  !parameters defining behavior FIXPVA
  real(r8_kind), parameter :: mm11=0.01058354498_r8_kind
  real(r8_kind), parameter :: mm12=0.29104748695_r8_kind
  real(r8_kind), parameter :: nn11=0.529177249_r8_kind
  real(r8_kind), parameter :: nn12=1.0583544980_r8_kind

  logical, public, allocatable :: skip_short(:)
!!$  integer, allocatable :: iuniq(:)

  !to provide graphical representation of tessera
  real(r8_kind), allocatable :: tess_export(:,:,:)
  integer(i4_kind), allocatable :: tess_exp_n(:)
  integer(i4_kind), allocatable :: vert_bounds(:,:,:)
  !to provide graphical representation of tessera centers (point charge positions)
  real(r8_kind), allocatable :: tess_centers(:,:)

  integer(i4_kind), parameter :: N_max_spheres=5000 !maximal number of spheres
  integer(i4_kind) :: NDIV=4

  integer(i4_kind),parameter :: MAX_ATOMS=40
  integer(i4_kind),parameter :: MAX_ATOMS1=100

  integer(kind=i4_kind) :: counter_n
  external error_handler

contains

  !********************************************************************
  function read_solv()

    logical :: read_solv

    integer(kind=i4_kind) :: i

    smoothing=df_smoothing
    electr_energy=df_electr_energy
    cavitation_energy=df_cavitation_energy
    disp_rep_energy=df_disp_rep_energy
    dielectric_constant=df_dielectric_constant
    abs_temperature=df_abs_temperature
    solvent_radius=df_solvent_radius
    solvent_volume=df_solvent_volume
    no_hydrogen_sphere=df_no_hydrogen_sphere
    hydrogen_no_scale=df_hydrogen_no_scale
    overlap_angle=df_overlap_angle
    fradio=df_fradio
    rmin_gepol=df_rmin_gepol
    min_area=df_min_area
    gepol=df_gepol
    save_cav=df_save_cav
    overlap_factor=df_overlap_factor
    output_level=df_output_level
    use_ewald=df_use_ewald

    if(calc_optimization .or. with_optimizer) fix_number_add=.true.

    call  go_to_first_input_line
    read_solv=find_namelist("&SOLV_OPTIONS",i)
    if(.not.read_solv) return
    read(input_device,nml=solv_options, end=100, err=200)

    if (solvent_radius <= 0.0)  goto 301
    if (dielectric_constant < 1.0 ) goto 302
    if(gepol /= 87 .and. gepol /= 93) gepol = 93
    if(overlap_factor < 0.0_r8_kind .or. overlap_factor >= 1.0_r8_kind) &
         overlap_factor=df_overlap_factor
    if(smoothing) then
       gepol = 93
       overlap_factor = 0.0_r8_kind
       scaled_factor = 1.125_r8_kind
       no_hydrogen_sphere = .false.
       hydrogen_no_scale = .true.
       fix_number_add = .false.
    endif
    if(overlap_factor == 0.0_r8_kind) rmin_gepol=0.8_r8_kind

    if(.not.electr_energy .and. .not.cavitation_energy .and. .not.disp_rep_energy) solvent=.false.
    read_solv=.true.
    return

100 read_solv=.false.
    return
200 call input_nm_error(0,"SOLV_OPTIONS")
301 call error_handler("MolMech: You have put wrong solvent_radius in namelist SOLV_OPTIONS")
302 call error_handler("MolMech: Wrong dielectric_constant in namelist SOLV_OPTION")

  end function  read_solv
  !********************************************************************

  !**************************************************************
  subroutine write_solv()

    character(len=76) :: message

    write(output_device,'(80("*"))')
    message='Solvent effect options:'
    write(output_device,'(a2,a76,a2)') '* ',message,' *'
    message='Calculated contributions:'
    write(output_device,'(a2,a76,a2)') '* ',message,' *'
    message=''
    if(electr_energy) message='Electrostatic'
    if(cavitation_energy) message=trim(message)//' Cavitation'
    if(disp_rep_energy) message=trim(message)//' Dispersion-Repulsion'
    write(output_device,'(a2,a76,a2)') '* ',message,' *'
    write(output_device,'(80("*"))')

  end subroutine write_solv
  !********************************************************************

  !**************************************************************
  subroutine points_on_cavity_surface
   !** End of interface *****************************************

    use cavity_image_module

    real(kind=r8_kind) :: radius
    type(triangles), target :: surf_elem                 !! for generate_dodecahedron
    real(kind=r8_kind), allocatable :: xyz_tes_c(:,:),area_tes(:),r_tes(:)           !! for tesselation
    integer(kind=i4_kind), allocatable :: sphere(:)                                  !! for tesselation
    logical, allocatable :: cuttt(:)                                                 !! for tesselation
    type(poligon),allocatable :: data_tes(:)                                         !! for tesselation
    integer (kind=i4_kind), allocatable :: cut_rad_sort(:,:)
    integer (i4_kind) :: i,j,status

    external error_handler

    if((slab_calc.or.lattice_calc) .and. .not.smoothing) then
       smoothing = .true.
       gepol = 93
       overlap_factor = 0.0_r8_kind
       scaled_factor = 1.125_r8_kind
       no_hydrogen_sphere = .false.
       hydrogen_no_scale = .true.
       fix_number_add = .false.
    endif

    if(do_electr .and. output_level >= 1) then
          write (output_device,*) '*******************************************************************'
          write (output_device,*) '          SOLVATION EFFECT - ELECTROSTATIC CONTRIBUTION          '
          write (output_device,*) '                     C O S M O - model                           '
          write (output_device,*) '==================================================================='
          write (output_device,*) 'Definition of cavity and coordinats of point charges on its surface'
          write (output_device,*) '==================================================================='
    endif

    if(do_cavitation .and. output_level >= 1) then
       write (output_device,*) '*******************************************************************'
       write (output_device,*) '          SOLVATION EFFECT - CAVITATION CONTRIBUTION'
       write (output_device,*) '                    Definition of cavity'
       write (output_device,*) '==================================================================='
    endif

    if(do_disp_rep .and. output_level >= 1) then
       write (output_device,*) '*******************************************************************'
       write (output_device,*) '       SOLVATION EFFECT - DISPERSION-REPULSION CONTRIBUTION'
       write (output_device,*) '                    Definition of cavity'
       write (output_device,*) '==================================================================='
    endif

    if(do_disp_rep .or. do_cavitation) then
       if(do_cavitation) R_access=0.0_r8_kind
       call calc_cavity_1()
       if(calc_gradients) call calc_cavity_gradient(.false.)
    else if(do_electr) then
       if(gepol == 87) then
          call calc_cavity()
       else if(gepol == 93) then
          call generate_dodecahedron(1)
          call calc_cavity_93()
          deallocate(surf_elem%xyz_centers,surf_elem%xyz,surf_elem%index)
       end if
       if(calc_gradients) then
          if(smoothing) then
             call calc_cavity_gradient(.false.)
          else
             call calc_cavity_gradient(.true.)
          end if
       end if
    end if

    if ( N_spheres >= MAX_ATOMS) then
       if( output_level >= 1)write (output_device,*) 'The OCTAHEDRON is inscribed into each sphere'
       call generate_octahedron()
    else
       if( output_level >= 1)write (output_device,*) 'The DODECAHEDRON is inscribed into each sphere'
       call generate_dodecahedron(0)
    end if

    if(.not. smoothing) then
       call init_cut_rad_sort()
    elseif(smoothing .and. .not.do_electr) then
       call init_cut_rad_sort()
    end if

    counter_n=0
    call tesselation

    if(do_electr.and.save_cav) then
       if(.not.smoothing) then
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
       call save_cavity_image(gepol,N_spheres,xyz_sphere,r_sphere,zero_area)

       allocate(tess_centers(N_total,3),stat=status)
       ASSERT(status==0)

       tess_centers=xyz_tes_c(1:N_total,:)
       call save_points_image(N_total,tess_centers)

       if(slab_calc.or.lattice_calc) then
          call save_cavity_and_points_image &
               (gepol,N_spheres,xyz_sphere,r_sphere,zero_area,N_total,tess_centers,N_atom_spheres)
       else
          call save_cavity_and_points_image &
               (gepol,N_spheres,xyz_sphere,r_sphere,zero_area,N_total,tess_centers)
       end if

       deallocate(tess_centers,stat=status)
       ASSERT(status==0)
    end if

    call symm_sorted_centers()
    if((do_electr .or. do_cavitation .or. do_disp_rep) .and. output_level >= 1) then
       write (output_device,*) '*****************************************************************'
    endif

  contains

    !------------------------------------------------------------
    ! Public interface of module
    subroutine calc_cavity_1
      !private subroutine
      !solvent accessible surface (only for dispersion-repulsion and cavitation terms)
      ! End of public interface of module

      integer(kind=i4_kind) :: typ
      integer(kind=i4_kind) :: i,k,status

      N_spheres = 0
      do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         N_spheres = N_spheres + 1
      enddo
      N_atom_spheres=N_spheres

      allocate(r_sphere(N_spheres),at_index(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface(1): allocation of r_sphere is failed")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface(1): allocation of xyz_sphere is failed")

      at_index(:)=0

      k=0
      n_uni: do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         k=k+1
         r_sphere(k)=atoms(typ)%r_vdw+R_access
         xyz_sphere(k,:)=atoms_cart(i)%r
         at_index(k)=i
      enddo n_uni

      if(output_level >= 1) then
         write (output_device,'(a23,i4,a19)') 'The cavity consists of ',N_spheres,' overlaping spheres'
         write (output_device,*) 'number       radius(angs)                 coordinates(angs)'
         do i=1,N_spheres
            write (output_device,'(1x,i4,6x,f11.8,6x,3(f13.9,1x))') i,r_sphere(i),xyz_sphere(i,:)
         enddo
         write (output_device,*) '---------------------------------------------------'
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
           dist_overlap, dist_h, dist_hiding
      real(kind=r8_kind) :: delta_xyz12(3), delta_xyz13(3), delta_xyz23(3)

      real(kind=r8_kind) :: r_solv,s_f,help,rmin,rmax
      real(kind=r8_kind) :: overlap_angle_rad,coso,sino

      integer(kind=i4_kind) :: one_more_sphere
      logical :: use_stored,parents_file_exists
      integer(kind=i4_kind) :: i,j,k,l,status,typ

      !converting parameters to radian and atomic units
      overlap_angle_rad= overlap_angle*deg2rad
      sino=sin(overlap_angle_rad)
      coso=cos(overlap_angle_rad)
      r_solv=solvent_radius !!!!!!!
      r_mini=rmin_gepol     !!!!!!
      s_f=scaled_factor

      parents_file_exists=.false.
      if(fix_number_add) call read_fixed_parents(parents_file_exists,n_max_par,fix_par)
      use_stored=(fix_number_add .and. parents_file_exists)
!!$      if(use_stored) write(output_unit,*) "THE CAVITY CONFIGURATIONS IS TAKEN FROM FIX_PAR FILE"

      N_spheres = 0
      do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         N_spheres = N_spheres + 1
      enddo
      N_atom_spheres=N_spheres

      allocate(r_sphere(N_spheres),at_index(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of r_sphere is failed")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of xyz_sphere is failed")
      allocate(parents(N_spheres,2), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of parents is failed")

      at_index(:)=0

      !definition of initial spheres
      k=0
      n_uni: do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         k=k+1
         r_sphere(k)=atoms(typ)%r_vdw*s_f
         xyz_sphere(k,:)=atoms_cart(i)%r
         parents(k,:)=0
         at_index(k)=i
      enddo n_uni

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
                  if(fix_par(n_par_fix,1)==i .and. fix_par(n_par_fix,2)==j) then
                     n_par_fix=n_par_fix+1
                  else
                     cycle main_cycle_2
                  endif
               endif

               delta_xyz12(:) = xyz_sphere(j,:) - xyz_sphere(i,:)
               distance12 = sqrt(dot_product(delta_xyz12,delta_xyz12))

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

                 if (distance12 <= dist_overlap) cycle  main_cycle_2

                 !are there any other spheres in between two spheres
                 check_hiding: do k=1,N_spheres
                   if (k == i .or. k == j) cycle

                   delta_xyz13(:) = xyz_sphere(k,:) - xyz_sphere(i,:)
                   distance13 = sqrt(dot_product(delta_xyz13,delta_xyz13))
                   delta_xyz23(:) = xyz_sphere(k,:) - xyz_sphere(j,:)
                   distance23 = sqrt(dot_product(delta_xyz23,delta_xyz23))

                   d_13=dot_product(delta_xyz13,delta_xyz12)/distance12
                   if(abs(d_13) <= 1.0e-8_r8_kind) d_13=0.0_r8_kind
                   d_23=dot_product(delta_xyz23,-delta_xyz12)/distance12
                   if(abs(d_23) <= 1.0e-8_r8_kind) d_23=0.0_r8_kind
                   if(d_13 <= 0.0_r8_kind .or. d_23 <= 0.0_r8_kind) cycle

                   help=1.0_r8_kind-(dot_product(delta_xyz13,delta_xyz12)/(distance13*distance12))**2
                   if(help<0.0_r8_kind) help=0.0_r8_kind
                   dist_h=distance13*sqrt(help)

                   dist_hiding = r_sphere(k)*fradio
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
              "MolMech:points_on_cavity_surface: deallocation of r_sphere is failed")
         deallocate(xyz_sphere, stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: deallocation of xyz_sphere is failed")
         deallocate(parents, stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: deallocation of parents is failed")

         N_spheres_old=N_spheres
         N_spheres= N_spheres+one_more_sphere

         allocate(r_sphere(N_spheres), stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: allocation of r_sphere is failed(1)")
         allocate(xyz_sphere(N_spheres,3), stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: allocation of xyz_sphere is failed(1)")
         allocate(parents(N_spheres,2), stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: allocation of parents is failed(1)")

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
              "MolMech:points_on_cavity_surface: deallocation r_sphere_buf is failed")
         deallocate(xyz_sphere_buf, stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: deallocation xyz_sphere_buf is failed")
         deallocate(parents_buf, stat=status)
         if ( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: deallocation of parents_buf is failed")

         if (one_more_sphere == 0 ) exit main_cycle_0
      enddo main_cycle_0

      if(fix_number_add .and. .not. parents_file_exists) then
         call write_fixed_parents(N_spheres,N_atom_spheres,parents)
      else if(fix_number_add) then
         deallocate(fix_par,stat=status)
         if( status /= 0) call error_handler( &
              "points_on_cavity_surface: deallocation of fix_par is failed")
      endif
!!$      if(output_cavity_data) then
!!$         write (output_unit,'(a23,i4,a19)') 'The cavity consists of ',N_spheres,' overlaping spheres'
!!$         write (output_unit,*) 'number       radius(a.u.)                 coordinates(a.u)'
!!$         do i=1,N_spheres
!!$            write (output_unit,'(1x,i4,6x,f11.8,6x,3(f13.9,1x))') i,r_sphere(i),xyz_sphere(i,:)
!!$         enddo
!!$         write (output_unit,*) '---------------------------------------------------'
!!$      endif

    end subroutine calc_cavity
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine read_fixed_parents(f_exists,n_max_par,fix_par)
      !** End of interface *****************************************
      use iounitadmin_module
      use filename_module
      integer(kind=i4_kind) :: st,ii
      integer(kind=i4_kind) :: iounit
      logical,intent(out) :: f_exists
      integer, intent(out) :: n_max_par
      integer, pointer :: fix_par(:,:)

      if(gepol==87) then
         inquire(file=trim(inpfile('fix_par.mm')), EXIST=f_exists)
         if(f_exists) then
            iounit = openget_iounit(file=trim(inpfile('fix_par.mm')), form='FORMATTED', &
                 action='READ',status='OLD')
         else
            return
         endif
      elseif(gepol==93) then
         inquire(file=trim(inpfile('fix_par.93.mm')), EXIST=f_exists)
         if(f_exists) then
            iounit = openget_iounit(file=trim(inpfile('fix_par.93.mm')), form='FORMATTED', &
                 action='READ',status='OLD')
         else
            return
         endif
      end if

      read(iounit, fmt=*, iostat=st) n_max_par
      if ( st .ne. 0 ) call error_handler("MolMech:read_fixed_parents: read n_max_par failed")

      if(gepol==87) then
         allocate(fix_par(n_max_par+1,2),stat=st)
      elseif(gepol==93) then
         allocate(fix_par(n_max_par+1,3),stat=st)
      end if
      if ( st .ne. 0 ) call error_handler("MolMech:read_fixed_parents: allocation failed")
      fix_par=0

      do ii=1,n_max_par
         if(gepol==87) then
            read(iounit, fmt=*, iostat=st) fix_par(ii,1:2)
         elseif(gepol==93) then
            read(iounit, fmt=*, iostat=st) fix_par(ii,1:3)
         end if
         if ( st .ne. 0 ) call error_handler("MolMech:read_fixed_parents: read fix_par failed")
      enddo

      call returnclose_iounit(iounit,status='KEEP')
    end subroutine read_fixed_parents
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine write_fixed_parents(n_sph,n_at_sph,fixed_par)
    !** End of interface *****************************************
      use iounitadmin_module
      use filename_module
      integer(kind=i4_kind) :: iostat, iounit,ii, n_par
      integer(kind=i4_kind), intent(in) :: n_sph,n_at_sph,fixed_par(:,:)

      n_par=n_sph-n_at_sph

      if(gepol == 87) then
         iounit = openget_iounit(file=trim(inpfile('fix_par.mm')), form='FORMATTED', &
              action='WRITE',status='NEW')
      elseif(gepol==93) then
         iounit = openget_iounit(file=trim(inpfile('fix_par.93.mm')), form='FORMATTED', &
              action='WRITE',status='NEW')
      end if

      write(iounit,iostat=iostat,fmt='(i5)') n_par
      if ( iostat .ne. 0 ) call error_handler("MolMech:write_fixed_parents: write failed (data_dir)" )

      do ii=1,n_par
         if(gepol == 87) then
            write(iounit,iostat=iostat,fmt='(2i5)') fixed_par(ii+n_at_sph,:)
         elseif(gepol==93) then
            write(iounit,iostat=iostat,fmt='(3i5)') fixed_par(ii+n_at_sph,:)
         end if
         if ( iostat .ne. 0 ) call error_handler("MolMech:write_fixed_parents: write failed (data)" )
      enddo

      call returnclose_iounit(iounit,status='KEEP')
    end subroutine write_fixed_parents
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine calc_cavity_93
      !Solvent excluding surface (Gepol93 algorithm)
      !** End of interface *****************************************
      real(r8_kind), allocatable :: r_sphere_buf(:)
      real(r8_kind), allocatable :: xyz_sphere_buf(:,:)
      integer(i4_kind), allocatable :: parents_buf(:,:)
      integer(i4_kind), pointer :: fix_par(:,:)
      integer(i4_kind) :: n_par_fix,n_max_par

      real(r8_kind) :: r_mini,r_max,r_min,r_solv,s_f,ofac
      integer(i4_kind) :: N_spheres_old,one_more_sphere
      real(r8_kind) :: distance12,dist13,dist23
      real(r8_kind) :: delta_xyz12(3),delta_xyz13(3)
      real(r8_kind) :: xyz_new(3),r_new,dist_1new,dist_2new,dist_3new
      real(r8_kind) :: xyz_surf(3),dist_sk
      integer(i4_kind) :: max_i,min_i

      logical :: parents_file_exists,use_stored
      integer(i4_kind) :: i,j,k,l,m,status,typ
      integer(i4_kind) :: N_S,K_S
      real(r8_kind) :: d_min,d_current2,r_current(3)

      !converting parameters to radian and atomic units
      r_solv=solvent_radius
      r_mini=rmin_gepol
      ofac=1.0_r8_kind-2.0_r8_kind*overlap_factor

      allocate(r_sphere_buf(N_max_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of r_sphere_buf is failed (93)")
      allocate(xyz_sphere_buf(N_max_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of xyz_sphere_buf is failed (93)")
      allocate(parents_buf(N_max_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of parents_buf is failed (93)")

      parents_file_exists=.false.
      if(fix_number_add) call read_fixed_parents(parents_file_exists,n_max_par,fix_par)
      use_stored=(fix_number_add .and. parents_file_exists)
      if(use_stored .and. output_level >=1 ) &
           write(output_device,*) "THE CAVITY CONFIGURATIONS IS TAKEN FROM FIX_PAR.93 FILE"

     if(smoothing.and.slab_calc) &
          call calc_images_slab_species_solvation(nn12)
     if(smoothing.and.lattice_calc) &
          call calc_images_species_solvation(nn12)
!!$print*,n_images_solv,nn12

      N_spheres = 0
      do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         N_spheres = N_spheres + 1
      enddo
      N_atom_spheres=N_spheres
      if(smoothing.and.(slab_calc.or.lattice_calc)) then
         do i=1,n_species
            typ=im_coor_solv(i)%type
            if(atoms(typ)%no_sphere) cycle
            N_spheres = N_spheres + n_images_solv
         end do
      end if

      allocate(at_index(N_spheres),stat=status)
      if(status/=0) call error_handler("MolMech:calc cavity_93: allocation at_index is failed")
      at_index(:)=0

      !definition of initial spheres
      k=0
      n_uni: do i=1,n_species
         typ=atoms_cart(i)%type
         if(atoms(typ)%no_sphere) cycle
         k=k+1
         if(smoothing) then
            s_f=atoms(typ)%scalefac_smooth
         else
            s_f=atoms(typ)%scalefac
         end if
         r_sphere_buf(k)=atoms(typ)%r_vdw*s_f!*a2b
         xyz_sphere_buf(k,:)=atoms_cart(i)%r!*a2b
         parents_buf(k,:)=0
         at_index(k)=i
      enddo n_uni
      if(smoothing.and.(slab_calc.or.lattice_calc)) then
         K_S=k
         do i=1,n_species
            typ=im_coor_solv(i)%type
            if(atoms(typ)%no_sphere) cycle
            s_f=atoms(typ)%scalefac_smooth
            do j=1,n_images_solv
               do l=1,K_S
                  d_min=r_sphere_buf(l)+atoms(typ)%r_vdw*s_f+nn12
                  r_current=xyz_sphere_buf(l,:)-im_coor_solv(i)%r(:,j)
                  d_current2=dot_product(r_current,r_current)
                  if(d_current2 <= d_min*d_min) goto 111  
               end do
               N_spheres=N_spheres-1
               cycle 
111            k=k+1
               r_sphere_buf(k)=atoms(typ)%r_vdw*s_f
               xyz_sphere_buf(k,:)=im_coor_solv(i)%r(:,j)
               at_index(k)=i
            end do
         end do
      end if
      ASSERT(k==N_spheres)

      do j=1,N_spheres
         do i=1,3
            if(abs(xyz_sphere_buf(j,i)) <= 1.0E-10_r8_kind) xyz_sphere_buf(j,i)=0.0_r8_kind
         enddo
      enddo

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
               if(use_stored ) then
                    if(fix_par(n_par_fix,1)==i .and. &
                       fix_par(n_par_fix,2)==j .and. &
                       fix_par(n_par_fix,3)==1) then
                       n_par_fix=n_par_fix+1
                     else
                        cycle second
                    endif
               endif

               !calculate the distance between the spheres
               delta_xyz12(:) = xyz_sphere_buf(j,:) - xyz_sphere_buf(i,:)
               distance12 = sqrt(dot_product(delta_xyz12,delta_xyz12))

               if(.not.use_stored) then
                  !can the solvent pass between spheres? if Yes choose new second sphere
                  if (distance12 >= (r_sphere_buf(i)+r_sphere_buf(j)+2.0_r8_kind*r_solv)) &
                       cycle second

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
                  if(distance12 <= r_max+r_min*ofac) cycle second

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
                     dist_3new = sqrt(dot_product(delta_xyz13,delta_xyz13))

                     !if no overlapping choose new the third sphere
                     if(dist_3new >= r_sphere_buf(k)+r_new) cycle third
                     !if overlapping
!!$                     if(r_new > r_sphere_buf(k)) cycle third
                     if(r_new >= r_sphere_buf(k)) cycle third         !one set of bulk spheres
                     !if overlapping too strong cancel the process and choose new second sphere
                     if(dist_3new <= r_sphere_buf(k)+r_new*ofac) cycle second
                  end do third

                  !Determine if new sphere contact with solvent
                  !(has non zero solvent accessible surface)
                  !if Yes choose new second sphere (define new pair)
                  surf_points:do k=1,N_centers_on_sphere
                     xyz_surf=surf_elem%xyz_centers(k,:)*(r_new+r_solv)/radius+xyz_new

                     atomss:do m=1,N_atom_spheres
!!$                        delta_xyz13(:) = xyz_new(:) - xyz_sphere_buf(m,:)
!!$                        dist_sk = sqrt(dot_product(delta_xyz13,delta_xyz13))
!!$                        if(dist_sk >= r_new+r_sphere_buf(k)+2.0_r8_kind*r_solv) cycle atoms

                        delta_xyz13(:) = xyz_surf(:) - xyz_sphere_buf(m,:)
                        dist_sk = sqrt(dot_product(delta_xyz13,delta_xyz13))

                        if(dist_sk < r_sphere_buf(m)+r_solv) cycle surf_points
                     end do atomss
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

               one_more_sphere=one_more_sphere+1
               xyz_sphere_buf(N_spheres+one_more_sphere,:)=xyz_new
               r_sphere_buf(N_spheres+one_more_sphere)=r_new
               parents_buf(N_spheres+one_more_sphere,1)=i
               parents_buf(N_spheres+one_more_sphere,2)=j
               parents_buf(N_spheres+one_more_sphere,3)=1
            end do second
         end do first

         N_spheres_old=N_spheres
         N_spheres=N_spheres+one_more_sphere

         exit bulk     !one set of bulk spheres
!!$         if(one_more_sphere == 0) exit bulk

      end do bulk

      do j=1,N_spheres
         do i=1,3
            if(abs(xyz_sphere_buf(j,i)) <= 1.0E-10_r8_kind) xyz_sphere_buf(j,i)=0.0_r8_kind
         enddo
      enddo

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
               if(use_stored ) then
                    if(fix_par(n_par_fix,1)==i .and. &
                       fix_par(n_par_fix,2)==j .and. &
                       fix_par(n_par_fix,3)==2) then
                       n_par_fix=n_par_fix+1
                    else
                        cycle second1
                    endif
               endif

               !calculate the distance between the spheres
               delta_xyz12(:) = xyz_sphere_buf(j,:) - xyz_sphere_buf(i,:)
               distance12 = sqrt(dot_product(delta_xyz12,delta_xyz12))

               if(.not.use_stored) then
                  !can the solvent pass between spheres? if Yes choose new second sphere
                  if (distance12 >= (r_sphere_buf(i)+r_sphere_buf(j)+2.0_r8_kind*r_solv)) &
                       cycle second1

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
                  if(distance12 < r_max+r_min*ofac) cycle second1
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
                  if(r_new <= r_mini) cycle second1
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
                     if(r_new <= r_mini) cycle second1
                     xyz_new=0.5_r8_kind*(xyz_sphere_buf(max_i,:)+xyz_sphere_buf(min_i,:))- &
                          (r_max-r_min)*(xyz_sphere_buf(max_i,:)-xyz_sphere_buf(min_i,:))/ &
                          (2.0_r8_kind*distance12)
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
                     dist_3new = sqrt(dot_product(delta_xyz13,delta_xyz13))

                     !if no overlapping choose new the third sphere
                     if(dist_3new >= r_sphere_buf(k)+r_new) cycle third1
                     !if overlapping
                     if(r_new >= r_sphere_buf(k)) cycle third1      !!! >
                     !if overlapping too strong cancel the process and choose new second sphere
                     if(dist_3new <= r_sphere_buf(k)+r_new*ofac) cycle second1
                  end do third1
               end if

               one_more_sphere=one_more_sphere+1
               xyz_sphere_buf(N_spheres+one_more_sphere,:)=xyz_new
               r_sphere_buf(N_spheres+one_more_sphere)=r_new
               parents_buf(N_spheres+one_more_sphere,1)=i
               parents_buf(N_spheres+one_more_sphere,2)=j
               parents_buf(N_spheres+one_more_sphere,3)=2
               if(N_spheres+one_more_sphere == N_max_spheres) exit first1
            end do second1
         end do first1

         N_spheres_old=N_spheres
         N_spheres=N_spheres+one_more_sphere
         if(N_spheres+one_more_sphere == N_max_spheres) exit crea

         if(one_more_sphere == 0) exit crea

      end do crea

1     allocate(r_sphere(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of r_sphere is failed (93)")
      allocate(zero_area(N_spheres), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of zero_area is failed (93)")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of xyz_sphere is failed (93)")
      allocate(parents(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of parents is failed (93)")

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
               dist_sk = sqrt(dot_product(delta_xyz13,delta_xyz13))

               if(dist_sk <= r_sphere_buf(j)) cycle three
            end do two
            zero_area(i)=.false.
            cycle one
         enddo three
      enddo one

      deallocate(r_sphere_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: deallocation of r_sphere_buf is failed (93)")
      deallocate(xyz_sphere_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: deallocation of xyz_sphere_buf is failed (93)")
      deallocate(parents_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: deallocation of parents_buf is failed (93)")

      if(fix_number_add .and. .not. parents_file_exists) then
         call write_fixed_parents(N_spheres,N_atom_spheres,parents)
      else if(fix_number_add) then
         deallocate(fix_par,stat=status)
         if( status /= 0) call error_handler( &
              "MolMech:points_on_cavity_surface: deallocation of fix_par is failed (93)")
      end if

      if(output_level >= 1) then
         N_S=N_spheres
         if(slab_calc.or.lattice_calc) N_S=N_atom_spheres
         write (output_device,'(a23,i4,a19)') 'The cavity consists of ',N_S,' overlaping spheres'
         write (output_device,*) 'number       radius(angs)                 coordinates(angs)'
         do i=1,N_S
            if(zero_area(i)) then
               write (output_device,'(1x,i4,6x,f11.8,6x,3(f13.9,1x),a)') i,r_sphere(i),xyz_sphere(i,:),' 0'
            else
               write (output_device,'(1x,i4,6x,f11.8,6x,3(f13.9,1x))')   i,r_sphere(i),xyz_sphere(i,:)
            end if
         enddo
         write (output_device,*) '---------------------------------------------------'
      endif

    end subroutine calc_cavity_93
    !------------------------------------------------------------
!!!MF - gradients >>>>
    !------------------------------------------------------------
    subroutine calc_cavity_gradient(cavity_flag)
      ! formulas see Cossi 1996
      !** End of interface *****************************************
      logical, intent(in) :: cavity_flag ! spheres with parents?

      integer(kind=i4_kind) :: i,j,k,l,l1,l2
      real(kind=r8_kind) :: r1s,r2s,ras,d12,vbuf(3),help
      real(kind=r8_kind) :: drdr1,drdr2,drda(3),dadr(3),dbdadiag1,dbdadiag2,&
           dbdanondiag,r_solv
      logical :: dist_check_j,dist_check_k
      real(kind=r8_kind), parameter :: small=1.e-8_r8_kind


      r_solv=solvent_radius
      call alloc_geom_deriv_part1(cagr)
      ! gradients of atomic centers and radia
      do i=1,N_atom_spheres
         k=at_index(i)
         do j=1,3
            cagr%dc(j,i)%xyz_grad(j,k,1)=1.0_r8_kind
         enddo
      enddo

      ! gradients of additional centers and radii
      ! spheres are ordered such that parents always before childs
      cav_f: if(cavity_flag) then
         do i=N_atom_spheres+1,N_spheres
            j=parents(i,1)
            k=parents(i,2)
            if(j>=i .or. k>=i) call error_handler("MolMech:calc_cavity_gradient: &
                 & wrong order of spheres")

            vbuf=xyz_sphere(j,:)-xyz_sphere(k,:)

            d12=sqrt(dot_product(vbuf,vbuf))
            r1s=r_sphere(j)+r_solv
            r2s=r_sphere(k)+r_solv
            ras=r_sphere(i)+r_solv

            vbuf=xyz_sphere(j,:)-r_sphere(j)/d12*(xyz_sphere(j,:)-xyz_sphere(k,:))&
                 -xyz_sphere(i,:)
            dist_check_j = dot_product(vbuf,vbuf)<small
            if(dist_check_j) then
               dist_check_k = .false.
            else
               vbuf=xyz_sphere(k,:)-r_sphere(k)/d12* &
                    (xyz_sphere(k,:)-xyz_sphere(j,:))&
                    -xyz_sphere(i,:)
               dist_check_k = dot_product(vbuf,vbuf)<small
            endif

            if(.not.dist_check_j .and. .not. dist_check_k) then
               ! case A/B formulas (18)-(21)
               drdr1=(3.0_r8_kind*r1s*(d12-r1s)-r2s*(d12-r2s)&
                    +2.0_r8_kind*r1s*r2s)/(4.0_r8_kind*d12*ras)
               drdr2=(3.0_r8_kind*r2s*(d12-r2s)-r1s*(d12-r1s)&
                    +2.0_r8_kind*r2s*r1s)/(4.0_r8_kind*d12*ras)
               drda(:)= (-d12**3+(r1s+r2s)*(r1s-r2s)**2)/&
                    (4.0_r8_kind*d12**3*ras) &
                    *(xyz_sphere(j,:)-xyz_sphere(k,:))
               dadr(:)=-(xyz_sphere(j,:)-xyz_sphere(k,:))/(2.0_r8_kind*d12)
               dbdadiag1= 0.5_r8_kind - (r1s-r2s)/(2.0_r8_kind*d12)
               dbdadiag2= 0.5_r8_kind + (r1s-r2s)/(2.0_r8_kind*d12)
               dbdanondiag = (r1s-r2s)/(2.0_r8_kind*d12**3)

               do l=1,n_species
                  do l1=1,3
                     help=xyz_sphere(j,l1)-xyz_sphere(k,l1) !beta_1-beta_2
                     do l2=1,3
                        cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,1) &
                             + dbdanondiag*help* &
                             (xyz_sphere(j,l2)-xyz_sphere(k,l2)) * &
                             cagr%dc(l2,j)%xyz_grad(:,l,1) &
                             + dbdanondiag*help* &
                             (xyz_sphere(k,l2)-xyz_sphere(j,l2)) * &
                             cagr%dc(l2,k)%xyz_grad(:,l,1)
                     enddo

                     cagr%dR(i)%xyz_grad(:,l,1) = cagr%dR(i)%xyz_grad(:,l,1) &
                          + drda(l1)*cagr%dc(l1,j)%xyz_grad(:,l,1) &
                          - drda(l1)*cagr%dc(l1,k)%xyz_grad(:,l,1)

                     cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                          cagr%dc(l1,i)%xyz_grad(:,l,1) &
                          + dadr(l1)*cagr%dR(j)%xyz_grad(:,l,1) &
                          - dadr(l1)*cagr%dR(k)%xyz_grad(:,l,1)

                     cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                          cagr%dc(l1,i)%xyz_grad(:,l,1) &
                          + dbdadiag1*cagr%dc(l1,j)%xyz_grad(:,l,1) &
                          + dbdadiag2*cagr%dc(l1,k)%xyz_grad(:,l,1)

                  enddo
                  cagr%dR(i)%xyz_grad(:,l,1) = cagr%dR(i)%xyz_grad(:,l,1) + &
                          drdr1*cagr%dR(j)%xyz_grad(:,l,1) + &
                          drdr2*cagr%dR(k)%xyz_grad(:,l,1)
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

               drdr1=(2.0_r8_kind*(d12*r1s+d12*r_sphere(j)- &
                    r1s*r_sphere(j)) - r1s**2 -r2s**2 +d12**2)/ &
                    (2.0_r8_kind*d12*ras)
               drdr2= (r2s*r_sphere(j))/(d12*r_sphere(i)**2)
               drda(:) = r_sphere(j)*(-d12**2+r1s**2-r2s**2)/ &
                    (2.0_r8_kind*d12**2 * ras**2) * &
                    (xyz_sphere(j,:)-xyz_sphere(k,:))
               dadr(:)=-(xyz_sphere(j,:)-xyz_sphere(k,:))/d12
               dbdadiag2=  r_sphere(j)/d12
               dbdadiag1= 1.0_r8_kind - dbdadiag2
               dbdanondiag = r_sphere(j)/d12**3

               do l=1,n_species
                  do l1=1,3
                     help=xyz_sphere(j,l1)-xyz_sphere(k,l1) !beta_1-beta_2
                     do l2=1,3
                        cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                             cagr%dc(l1,i)%xyz_grad(:,l,1) &
                             + dbdanondiag*help* &
                             (xyz_sphere(j,l2)-xyz_sphere(k,l2)) * &
                             cagr%dc(l2,j)%xyz_grad(:,l,1) &
                             + dbdanondiag*help* &
                             (xyz_sphere(k,l2)-xyz_sphere(j,l2)) * &
                             cagr%dc(l2,k)%xyz_grad(:,l,1)
                     enddo

                     cagr%dR(i)%xyz_grad(:,l,1) = cagr%dR(i)%xyz_grad(:,l,1) &
                          + drda(l1)*cagr%dc(l1,j)%xyz_grad(:,l,1) &
                          - drda(l1)*cagr%dc(l1,k)%xyz_grad(:,l,1)

                     cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                          cagr%dc(l1,i)%xyz_grad(:,l,1) &
                          + dadr(l1)*cagr%dR(j)%xyz_grad(:,l,1)

                     cagr%dc(l1,i)%xyz_grad(:,l,1) = &
                          cagr%dc(l1,i)%xyz_grad(:,l,1) &
                          + dbdadiag1*cagr%dc(l1,j)%xyz_grad(:,l,1) &
                          + dbdadiag2*cagr%dc(l1,k)%xyz_grad(:,l,1)

                  enddo
                  cagr%dR(i)%xyz_grad(:,l,1) = cagr%dR(i)%xyz_grad(:,l,1) + &
                       drdr1*cagr%dR(j)%xyz_grad(:,l,1) + &
                       drdr2*cagr%dR(k)%xyz_grad(:,l,1)
               enddo
            endif
         enddo
      endif cav_f

    end subroutine calc_cavity_gradient
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine alloc_geom_deriv_part1(ca)
      !allocation of dR,dc gradients of geom_deriv ca
      type(geom_deriv) :: ca
      !** End of interface *****************************************

      integer(kind=i4_kind) :: i,alloc_stat,j

      allocate(ca%dR(N_spheres),stat=alloc_stat)
      if ( alloc_stat /= 0) call error_handler( &
           "MolMech:alloc_geom_deriv_part: allocation of grad_atomic_center type  array is failed")
      allocate(ca%dc(3,N_spheres),stat=alloc_stat)
      if ( alloc_stat /= 0) call error_handler( &
           "MolMech:alloc_geom_deriv_part: allocation of grad_atomic_center type  array is failed")

      do i=1,N_spheres
         allocate(ca%dR(i)%xyz_grad(3,n_species,1),stat=alloc_stat)
         if ( alloc_stat /= 0) call error_handler( &
              "MolMech:alloc_geom_deriv_part: allocation of xyz_grad is failed")
         ca%dR(i)%xyz_grad(:,:,:)=0.0_r8_kind
         do j=1,3
            allocate(ca%dc(j,i)%xyz_grad(3,n_species,1),&
                 stat=alloc_stat)
            if ( alloc_stat /= 0) call error_handler( "MolMech:alloc_geom_&
                 &deriv_part: allocation of xyz_grad is failed")
            ca%dc(j,i)%xyz_grad(:,:,:)=0.0_r8_kind
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

      integer(kind=i4_kind) :: i,alloc_stat,j

      allocate(ca%darea(3,n_species,1),stat=alloc_stat)
      allocate(ca%dcenter(3,n_species,1),stat=alloc_stat)

      do i=1,n_species
         do j=1,3
            allocate(ca%darea(j,i,1)%m(n_max),stat=alloc_stat)
            allocate(ca%dcenter(j,i,1)%m(3,n_max),stat=alloc_stat)
            ca%darea(j,i,1)%m=0.0_r8_kind
            ca%dcenter(j,i,1)%m=0.0_r8_kind
         enddo
      enddo

    end subroutine alloc_geom_deriv_part2
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine generate_octahedron()
      !generate set of surface points on a sphere starting from a cube
      !and do subdivisions
      !** End of interface *****************************************

      real(kind=r8_kind),pointer :: xyz(:,:),xyz_cent(:,:)
      integer(kind=i4_kind), pointer :: indexx(:,:)
      real(kind=r8_kind) :: a,d


      real(kind=r8_kind) :: xyz_buf(3)
      integer(kind=i4_kind) :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
      integer(kind=i4_kind) :: N_dim_xyz_next,N_dim_cent_next
      integer(kind=i4_kind) :: i,j,l1,n1,i1,status1
      integer(kind=i4_kind) :: local_point_factor


      a=1.0_r8_kind
      radius=a
      d=sqrt(2.0_r8_kind)*a

      N_dim_cent=8
      N_dim_cent_st=8
      N_dim_xyz=6
      N_dim_xyz_st=6

      local_point_factor=2
      if(N_atom_spheres >= MAX_ATOMS1) local_point_factor=1

      i=1
      do
         if(i>local_point_factor-1) exit
         N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
         N_dim_cent=N_dim_cent*4
         if(i==1) then
            N_dim_xyz_next=N_dim_xyz
            N_dim_cent_next=N_dim_cent
         endif
         i=i+1
      enddo

      allocate(surf_elem%xyz(N_dim_xyz,3),stat=status1)
      if ( status1 /= 0) call error_handler &
           ("MolMech:generate_octaahedron: allocation surf_elem%xyz is failed")

      !definition of initial set of points
      xyz=>surf_elem%xyz
      xyz=0.0_r8_kind

      xyz(1,1)=a
      xyz(2,1)=-a
      xyz(3,2)=a
      xyz(4,2)=-a
      xyz(5,3)=a
      xyz(6,3)=-a

      !definition of triangles

      allocate(surf_elem%index(N_dim_cent,3),stat=status1)
      if ( status1 /= 0) call error_handler( &
           "MolMech:generate_octahedron: allocation surf_elem%index is failed")
      indexx=>surf_elem%index

      indexx(1,1)=1
      indexx(1,2)=3
      indexx(1,3)=5

      indexx(2,1)=1
      indexx(2,2)=4
      indexx(2,3)=5

      indexx(3,1)=1
      indexx(3,2)=3
      indexx(3,3)=6

      indexx(4,1)=1
      indexx(4,2)=4
      indexx(4,3)=6

      indexx(5,1)=2
      indexx(5,2)=3
      indexx(5,3)=5

      indexx(6,1)=2
      indexx(6,2)=4
      indexx(6,3)=5

      indexx(7,1)=2
      indexx(7,2)=3
      indexx(7,3)=6

      indexx(8,1)=2
      indexx(8,2)=4
      indexx(8,3)=6

      do j=1,local_point_factor-1
#ifdef _LINUX1
         call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
         call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
         N_dim_xyz_st=N_dim_xyz_next
         N_dim_cent_st=N_dim_cent_next
         N_dim_xyz_next= N_dim_xyz_st+(N_dim_cent_st*3)/2
         N_dim_cent_next= N_dim_cent_st*4
      enddo

      ! definition of triangle centers
      allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status1)
      if ( status1 /= 0) call error_handler( &
           "MolMech:generate_dodecahedron: allocation surf_elem%xyz_centers is failed")

      xyz_cent=>surf_elem%xyz_centers
      do i=1,N_dim_cent_st
         i1=indexx(i,1)
         l1=indexx(i,2)
         n1=indexx(i,3)
         xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
         xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
      enddo
      N_points_of_triangles=N_dim_xyz_st
      N_centers_on_sphere=N_dim_cent_st

      if(output_level >= 1) &
           write (output_device,'(a42,i5)') &
           'the number of triangles on each sphere is ', N_centers_on_sphere

    end subroutine generate_octahedron
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine generate_dodecahedron(gp93)
      !generate points due to a pentakis-dodecahedron
      !rotate in right symmetry
      !subdivide if needed
      !** End of interface *****************************************

      integer(kind=i4_kind), intent(in) :: gp93
      real(kind=r8_kind),pointer :: xyz(:,:),xyz_cent(:,:)
      integer(kind=i4_kind), pointer :: indexx(:,:)
      real(kind=r8_kind) :: t,s,h,r,b,c1,c2,h1,h2,rr,cosT,sinT

      real(kind=r8_kind) :: r_buf,xyz_buf(3)
      integer(kind=i4_kind) :: tang_buf(5),local_point_factor
      integer(kind=i4_kind) :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
      integer(kind=i4_kind) :: N_dim_xyz_next,N_dim_cent_next
      real(kind=r8_kind), parameter :: cos12=0.8660254037844_r8_kind
      real(kind=r8_kind), parameter :: sin12=0.5_r8_kind
      real(kind=r8_kind), parameter :: cos20=0.9510565162952_r8_kind
      real(kind=r8_kind), parameter :: sin20=0.3090169943749_r8_kind
      integer(kind=i4_kind) :: i,j,k,l,m,n,l1,n1,ind,i1,status2

      !definition a pentakis dodecahedron
      t=(sqrt(5.0_r8_kind)+1.0_r8_kind)/2.0_r8_kind
      s=sqrt(3.0_r8_kind-t)/2.0_r8_kind
      h=(2.0_r8_kind*s+(t+1.0_r8_kind)*sqrt(t+2.0_r8_kind))/4.0_r8_kind
      radius=sqrt(h**2+s**2)
      b=radius*sqrt(3.0_r8_kind)/3.0_r8_kind
      c1=(2.0_r8_kind*b+h)/5.0_r8_kind
      c2=(2.0_r8_kind*h+2.0_r8_kind*b+s)/5.0_r8_kind
      rr=sqrt(c1**2+c2**2)
      h1=c1*radius/rr
      h2=c2*radius/rr
      cosT=h2/sqrt(h1**2+h2**2)
      sinT=sqrt(1.0_r8_kind-cosT**2)
      r=sqrt(radius**2-1.0_r8_kind)

      N_dim_cent=60
      N_dim_cent_st=60
      N_dim_xyz=32
      N_dim_xyz_st=32

      local_point_factor=1

      if(gp93==1) then
         local_point_factor=NDIV
         i=1
         do
            if(i >local_point_factor-1) exit
            N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
            N_dim_cent=N_dim_cent*4
            if(i==1) then
               N_dim_xyz_next=N_dim_xyz
               N_dim_cent_next=N_dim_cent
            endif
            i=i+1
         end do
      else
         i=1
         do
            if(i>local_point_factor-1) exit
            N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
            N_dim_cent=N_dim_cent*4
            if(i==1) then
               N_dim_xyz_next=N_dim_xyz
               N_dim_cent_next=N_dim_cent
            endif
            i=i+1
         enddo
      end if

      allocate(surf_elem%xyz(N_dim_xyz,3),stat=status2)
      if ( status2 /= 0) call error_handler &
           ("MolMech:generate_dodecahedron: allocation surf_elem%xyz is failed")

      !definition of initial set of points
      xyz=>surf_elem%xyz
      xyz=0.0_r8_kind

      xyz(1:4,1)=0.0_r8_kind
      xyz(5,1)=s
      xyz(7,1)=s
      xyz(6,1)=-s
      xyz(8,1)=-s
      xyz(9:10,1)=h
      xyz(11:12,1)=-h
      xyz(13,1)=b
      xyz(15,1)=b
      xyz(17,1)=b
      xyz(19,1)=b
      xyz(14,1)=-b
      xyz(16,1)=-b
      xyz(18,1)=-b
      xyz(20,1)=-b
      xyz(21,1)=h1
      xyz(22,1)=-h1
      xyz(23,1)=h1
      xyz(24,1)=-h1
      xyz(25:29,1)=0.0_r8_kind
      xyz(29:30,1)=h2
      xyz(31:32,1)=-h2

      xyz(1,2)=-s
      xyz(3,2)=-s
      xyz(2,2)=s
      xyz(4,2)=s
      xyz(5:6,2)=h
      xyz(7:8,2)=-h
      xyz(9:12,2)=0.0_r8_kind
      xyz(13:14,2)=b
      xyz(17:18,2)=b
      xyz(15:16,2)=-b
      xyz(19:20,2)=-b
      xyz(21:24,2)=0.0_r8_kind
      xyz(25,2)=-h2
      xyz(26:27,2)=h2
      xyz(28,2)=-h2
      xyz(29,2)=-h1
      xyz(30,2)=h1
      xyz(31,2)=-h1
      xyz(32,2)=h1

      xyz(1:2,3)=h
      xyz(3:4,3)=-h
      xyz(5:8,3)=0.0_r8_kind
      xyz(9,3)=s
      xyz(11,3)=s
      xyz(10,3)=-s
      xyz(12,3)=-s
      xyz(13:16,3)=b
      xyz(17:20,3)=-b
      xyz(21:22,3)=h2
      xyz(23:24,3)=-h2
      xyz(25:26,3)=h1
      xyz(27:28,3)=-h1
      xyz(29:32,3)=0.0_r8_kind

      !definition of triangles

      allocate(surf_elem%index(N_dim_cent,3),stat=status2)
      if ( status2 /= 0) call error_handler( &
           "MolMech:generate_dodecahedron: allocation surf_elem%index is failed")
      indexx=>surf_elem%index
      ind=1
      do i=1,12
         k=i+20
         m=5*(i-1)+1
         indexx(m:m+4,1)=k

         n=1
         do j=1,20
            xyz_buf=xyz(k,:)-xyz(j,:)
            r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
            if (r_buf <= 1.06_r8_kind) then
               tang_buf(n)=j
               n=n+1
               if (n > 5) exit
            endif

         enddo

         do l=1,4
            l1=tang_buf(l)
            do n=l+1,5
               n1=tang_buf(n)
               xyz_buf=xyz(n1,:)-xyz(l1,:)
               r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
               if (r_buf <= 2.0_r8_kind*s+0.001_r8_kind) then
                  indexx(ind,2)=l1
                  indexx(ind,3)=n1
                  ind=ind+1
               endif
            enddo
         enddo
      enddo

      do j=1,local_point_factor-1
#ifdef _LINUX1
         call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
         call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
         N_dim_xyz_st=N_dim_xyz_next
         N_dim_cent_st=N_dim_cent_next
         N_dim_xyz_next= N_dim_xyz_st+(N_dim_cent_st*3)/2
         N_dim_cent_next= N_dim_cent_st*4
      enddo

      ! definition of triangle centers
      allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status2)
      if ( status2 /= 0) call error_handler( &
           "MolMech:generate_dodecahedron: allocation surf_elem%xyz_centers is failed")

      xyz_cent=>surf_elem%xyz_centers
      do i=1,N_dim_cent_st
         i1=indexx(i,1)
         l1=indexx(i,2)
         n1=indexx(i,3)
         xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
         xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
      enddo
      N_points_of_triangles=N_dim_xyz_st
      N_centers_on_sphere=N_dim_cent_st
      if(output_level >= 1 .and. gp93 /= 1) &
           write (output_device,'(a42,i5)') &
           'the number of triangles on each sphere is ', N_centers_on_sphere

    end subroutine generate_dodecahedron
    !------------------------------------------------------------

    !------------------------------------------------------------
#ifdef _LINUX1
    subroutine more_triangles(n_ind,n_xyz_st,n_ind_st)
#else
    subroutine more_triangles(xyz_t,ind_t,n_ind,n_xyz_st,n_ind_st)
#endif
      !subdivision routine, each triangle is devided in four by
      !halving the edges
      !** End of interface *****************************************

      real(kind=r8_kind), pointer :: xyz_t(:,:)
      integer(kind=i4_kind), pointer :: ind_t(:,:)
      integer(kind=i4_kind), intent(in) :: n_ind
      integer(kind=i4_kind), intent(inout) :: n_xyz_st,n_ind_st

      real(kind=r8_kind) :: xyz_t_buf(3),xyz_tt(3)
      integer(kind=i4_kind), allocatable :: ind_buf(:,:)
      integer(kind=i4_kind) :: new_numbers(6)
      integer(kind=i4_kind) :: status,next_i,neighbour,i1,i2
      integer(kind=i4_kind) :: i,j,k,l,m,m1,n,n1
      real(kind=r8_kind), parameter :: small=1.0e-11_r8_kind

#ifdef _LINUX1
      xyz_t=>surf_elem%xyz
      ind_t=>surf_elem%index
#endif

      allocate(ind_buf(n_ind,3),stat=status)
      if ( status /= 0) call error_handler( &
           "more_triangles: allocation ind_buf is failed")
      ind_buf=0


      next_i=0
      label_i: do i=1,n_ind_st

         m=0
         label_j: do j=1,3
            next_i=next_i+1
            neighbour=1
            i1=ind_t(i,j)
            ind_buf(next_i,neighbour)=i1

            label_k: do k=1,3
               if (k==j) cycle label_k
               i2=ind_t(i,k)
               xyz_t_buf=(xyz_t(i1,:)+xyz_t(i2,:))/2.0_r8_kind
               xyz_tt=(radius/sqrt(dot_product(xyz_t_buf,xyz_t_buf)))*xyz_t_buf

               label_l: do l=1,n_xyz_st

                  if(abs(xyz_tt(1)-xyz_t(l,1))<small .and. &
                       abs(xyz_tt(2)-xyz_t(l,2))<small .and. &
                       abs(xyz_tt(3)-xyz_t(l,3))<small) then
                     neighbour=neighbour+1
                     ind_buf(next_i,neighbour)=l
                     m=m+1
                     new_numbers(m)=l
                     cycle label_k
                  endif
               enddo label_l
               n_xyz_st=n_xyz_st+1
               xyz_t(n_xyz_st,:)=xyz_tt
               neighbour=neighbour+1
               ind_buf(next_i,neighbour)=n_xyz_st
               m=m+1
               new_numbers(m)=n_xyz_st
            enddo label_k
         enddo label_j
         next_i=next_i+1
         m1=1
         ind_buf(next_i,m1)=new_numbers(1)
         label_n: do n=2,m

            label_n1: do n1=1,n-1
               if(new_numbers(n) == new_numbers(n1)) cycle label_n
            enddo label_n1
            m1=m1+1
            ind_buf(next_i,m1)=new_numbers(n)
            if (m1==3) cycle label_i
         enddo label_n
      enddo label_i

      n_ind_st=next_i
      ind_t(1:n_ind,1:3)=ind_buf(1:n_ind,1:3)

      deallocate(ind_buf,stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:more_triangles: deallocation ind_buf is failed")
!     nullify(xyz_t,ind_t)
    end subroutine more_triangles
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine tesselation
      !creates points on each sphere
      !calls def_poligon_tesselation ("cutting routine")
      !stores resulting surface tesselation in data_tes
      !** End of interface *****************************************

      type(triangles), allocatable :: triangles_on_sphere(:)
      type(poligon) :: tes_temp
      real(kind=r8_kind) :: xyz_pol_c(3),area_pol,area_tot
      integer(kind=i4_kind) :: status,i,j,k,l,i1,j1,k1,l1
      integer(kind=i4_kind), allocatable :: i_sorted_spheres(:),N_cent_loc(:)
      integer(kind=i4_kind) :: ihelp
      real(kind=r8_kind) :: xyz_buf(3)
      integer(i4_kind) :: N_S

      N_S=N_spheres
      if(slab_calc.or.lattice_calc) N_S=N_atom_spheres

      allocate(i_sorted_spheres(N_s), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of i_sorted_spheres is failed")
      allocate(N_cent_loc(N_s), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of N_cent_loc is failed")

      i_sorted_spheres(1)=1
      do i=2,N_s
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

      allocate(triangles_on_sphere(N_s), stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:tesselation: allocation of triangles_on_sphere is failed")

      !projecting set of points of arbitrary sphere on each sphere forming the cavity
      labisorted: do ihelp=1,N_s
         i=i_sorted_spheres(ihelp)

         N_cent_loc(i)=N_centers_on_sphere
         do k=1,N_centers_on_sphere
            i1=surf_elem%index(k,1)
            l1=surf_elem%index(k,2)
            k1=surf_elem%index(k,3)
            xyz_buf=(surf_elem%xyz(i1,:)+surf_elem%xyz(l1,:)+surf_elem%xyz(k1,:))/3.0_r8_kind
            surf_elem%xyz_centers(k,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
         enddo
         allocate(triangles_on_sphere(i)%xyz(N_points_of_triangles,3), &
              triangles_on_sphere(i)%index(N_centers_on_sphere,3), &
              triangles_on_sphere(i)%xyz_centers(N_centers_on_sphere,3), &
              stat=status)
         if ( status /= 0) then
            call error_handler( &
                 "MolMech:tesselation: allocation of triangles_on_sphere(1) is failed")
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
         if(smoothing .and. do_electr) then
            triangles_on_sphere(i)%area=(4.0_r8_kind*pi*r_sphere(i)**2)/N_centers_on_sphere
         else
            triangles_on_sphere(i)%area=0.0_r8_kind
         end if
      enddo labisorted

      deallocate(i_sorted_spheres,surf_elem%xyz_centers,surf_elem%xyz,surf_elem%index, stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:tesselation: deallocation of i_sorted_spheres is failed")

      allocate(xyz_tes_c(N_s*N_centers_on_sphere,3), &
           area_tes(N_s*N_centers_on_sphere), &
           r_tes(N_s*N_centers_on_sphere), &
           sphere(N_s*N_centers_on_sphere), &
           cuttt(N_s*N_centers_on_sphere), &
           data_tes(N_s*N_centers_on_sphere),stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:tesselation: allocation of xyz_tes_c, area_tes, r_tes are failed")

      if(calc_gradients) then
         call alloc_geom_deriv_part2(cagr,N_s*N_centers_on_sphere)
      endif

      !definition of the cavity elements
      area_tot=0.0_r8_kind
      N_total=0
      lab1: do i=1,N_s
         if(do_electr) then
            if(gepol == 93) then
               if(zero_area(i)) cycle lab1
            end if
         end if
         lab2: do j=1,N_cent_loc(i)
           j1=triangles_on_sphere(i)%index(j,1)
           k1=triangles_on_sphere(i)%index(j,2)
           l1=triangles_on_sphere(i)%index(j,3)
!!$print*,i,j
!!$print*,triangles_on_sphere(i)%xyz(j1,:)
!!$print*,triangles_on_sphere(i)%xyz(k1,:)
!!$print*,triangles_on_sphere(i)%xyz(l1,:)
!!$           if(smoothing .and. not_cav_dis_rep) then
           if(smoothing .and. do_electr) then
              xyz_pol_c=triangles_on_sphere(i)%xyz_centers(j,:)
              area_pol=triangles_on_sphere(i)%area
              call smooth_tess_area1(i,j,xyz_pol_c,area_pol,tes_temp)
           elseif(.not.smoothing .or. .not.do_electr) then
             call def_poligon_tes(i, &
                   triangles_on_sphere(i)%xyz(j1,:), &
                   triangles_on_sphere(i)%xyz(k1,:), &
                   triangles_on_sphere(i)%xyz(l1,:), &
                   xyz_pol_c,area_pol,tes_temp)
           end if

           if(area_pol /= 0.0_r8_kind) then
              N_total=N_total+1
              xyz_tes_c(N_total,:)=xyz_pol_c
              area_tes(N_total)=area_pol
              area_tot=area_tot+area_tes(N_total)
              r_tes(N_total)=triangles_on_sphere(i)%radius
              sphere(N_total)=i

              cuttt(N_total)=.false.
              do l=1,tes_temp%n_vertises
                 if(tes_temp%n_sphere(l,1)/=0) cuttt(N_total) = .true.
              enddo

              data_tes(N_total)=tes_temp
           endif
        enddo lab2
      enddo lab1

      if(output_level >= 1) then
         write(output_device,*) '-------------------------------------------------------------'
         write(output_device,'(i5,a27)') N_total,' points have been generated'
         write(output_device,'(f18.9,a27)') area_tot,' ang^2 total surface area  '
         if(output_level >= 2) then
            write(output_device,*)'number                   coordinates(angs)               area   n_vert sphere'
            do i=1,N_total
               write(output_device,'(1x,i4,3x,3(1x,f13.9),1x,f14.9,1x,i3,1x,i3)') i,xyz_tes_c(i,:), &
                 area_tes(i),data_tes(i)%n_vertises,data_tes(i)%sphere
!!$            do j=1,data_tes(i)%n_vertises
!!$               write(output_device,*) data_tes(i)%xyz_vertex(j,:)
!!$            end do
            enddo
         endif
         write(output_device,*) '-------------------------------------------------------------'
      endif

      do i=1,N_s
         deallocate(triangles_on_sphere(i)%xyz, &
              triangles_on_sphere(i)%index, &
              triangles_on_sphere(i)%xyz_centers, &
              stat=status)
         if ( status /= 0) then
            call error_handler( &
                 "MolMech:tesselation: deallocation of triangles_on_sphere(1) is failed")
         endif

      enddo
      deallocate(triangles_on_sphere,stat=status)
      if ( status /= 0) call error_handler( &
         "MolMech:tesselation: deallocation of triangles_on_sphere(2) is failed")

      if(allocated(cut_rad_sort)) then
         deallocate(cut_rad_sort, stat=status)
         ASSERT(status==0)
      end if
      deallocate(N_cent_loc, stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: deallocation of N_cent_loc is failed")

      if(calc_gradients .and. .not. do_disp_rep) then
         call dealloc_geom_deriv_part1(cagr)
      endif

    end subroutine tesselation
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
         dr2r3=sqrt(dot_product(r2r3,r2r3))
         r1r3=r1-r3
         dr1r3=sqrt(dot_product(r1r3,r1r3))

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
         nn=dot_product(r1r5,r1r5)

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

      if(calc_gradients) then
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
    subroutine init_cut_rad_sort()
      implicit none
      ! *** end interface ***

      integer(kind=i4_kind) :: i,j,status,k,l
      real (kind=r8_kind),allocatable :: rcut(:)
      real (kind=r8_kind) :: help,vd(3)

      allocate(cut_rad_sort(N_spheres,N_spheres),rcut(N_spheres),stat=status)
      if(status/=0) call error_handler("MolMech:alloc of cut_rad_sort failed")
      cut_rad_sort=0

      do i=1,N_spheres
!!$     print*,xyz_sphere(i,:),i,'xyz_sphere'
         rcut=-1.0_r8_kind
         labj : do j=1,N_spheres
            if(i==j) cycle labj
            vd(:)=xyz_sphere(j,:)-xyz_sphere(i,:)
            help=sqrt(dot_product(vd,vd))
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
         xyz_triangl_1,xyz_triangl_2,xyz_triangl_3, &
         xyz_poligon_centre,poligon_area,poli)
      ! This is very important procedure calculating area and coordinates
      ! of geometrical center of tessarea of polygonal shape.
      ! The polygonal type of tessarea is formed as result of intersection
      ! between two or some spheres.
      !** End of interface *****************************************
      integer(kind=i4_kind), intent(in) :: sphere_number
      real(kind=r8_kind), intent(in) :: xyz_triangl_1(3), xyz_triangl_2(3), xyz_triangl_3(3)
      real(kind=r8_kind), intent(out) :: xyz_poligon_centre(3), poligon_area
      type(poligon), intent(out) :: poli

      real(kind=r8_kind) :: new_vertex(6)
      real(kind=r8_kind) :: alpha, beta, cos_gamma,v1(3),v2(3),v_buf(3)
      logical :: combination(MAX_POL_VER,MAX_POL_VER)
      logical :: inside_flag(MAX_POL_VER),cutted
      real(kind=r8_kind) :: d_rad,d_vv
      real(kind=r8_kind), parameter :: small_d=1.0e-8_r8_kind

      integer(kind=i4_kind) :: ncut,k_new,l1,l2,n_extra
      integer(kind=i4_kind) :: i,i1,ii,j,k,jh
      real(kind=r8_kind) :: tn(MAX_POL_VER,3),tns(MAX_POL_VER,3),vn(MAX_POL_VER,3),vns(MAX_POL_VER,3)
      logical :: extra_point(MAX_POL_VER),h_case(MAX_POL_VER,2),h_case_sphere
      real(kind=r8_kind), dimension(MAX_POL_VER) :: costhn,phin
      real(kind=r8_kind) :: cut_center(3),cut_radius,help,center_factor
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

      poli%sphere=sphere_number
      counter_n=counter_n+1

      poli%n_vertises=3

      poli%xyz_vertex(1,:)=xyz_triangl_1
      poli%xyz_vertex(2,:)=xyz_triangl_2
      poli%xyz_vertex(3,:)=xyz_triangl_3

      ! ensure anticlockwise sort
      v_buf=vector_product(xyz_triangl_1-xyz_sphere(sphere_number,:),&
           xyz_triangl_2-xyz_sphere(sphere_number,:))

      if(dot_product(v_buf,xyz_triangl_3-xyz_sphere(sphere_number,:))&
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
         if(do_electr) then
            if(gepol == 93) then
               if(zero_area(j)) cycle j_spheres
            end if
         end if
         v1=xyz_sphere(j,:)-xyz_sphere(sphere_number,:)
         d_rad=sqrt(dot_product(v1,v1))
         if (r_sphere(sphere_number)+r_sphere(j) < d_rad) cycle j_spheres

         d_vv=(d_rad**2+r_sphere(sphere_number)**2-r_sphere(j)**2)/&
              (2.0_r8_kind*d_rad)

         cut_center=xyz_sphere(sphere_number,:)+v1/d_rad*d_vv
         help=r_sphere(sphere_number)**2-d_vv**2
         if(help<0.0_r8_kind) help=0.0_r8_kind
         cut_radius=sqrt(help)

         !h_case: sphere more than half in an other, typically
         !for H-Atoms
         if(dot_product(v1,cut_center(:)-xyz_sphere(sphere_number,:))<0.0_r8_kind) then
            h_case_sphere=.true.
         endif

         ! definition of the vertises which are inside of J-th sphere

         do i=1,poli%n_vertises
            d_rad=sqrt(dot_product(poli%xyz_vertex(i,:)-xyz_sphere(j,:), &
                 poli%xyz_vertex(i,:)-xyz_sphere(j,:)))
            ! is vertex inside of sphere?
            if (d_rad < r_sphere(j) - small_d) then
               inside_flag(i)= .true.
               cutted=.true.
            endif
         enddo

         k=poli%n_vertises+1
         i_vertex: do i=1,poli%n_vertises
            i1=poli%bounds(i,1)

!!$print*,i,sphere_number,j,inside_flag(i),inside_flag(i1),'!!!!!!!!!!!!!'
            call line_cut(poli%xyz_vertex(i,:),poli%xyz_vertex(i1,:),&
                 poli%xyz_bound(i,1:3),poli%r_bound(i,1),&
                 xyz_sphere(j,:),r_sphere(j),&
                 inside_flag(i),inside_flag(i1), &
                 xyz_sphere(sphere_number,:),h_case(i,1),&
                 ncut,new_vertex)
!!$print*,i,sphere_number,j,inside_flag(i),inside_flag(i1),ncut,'!!!!!!!!!!!!!'

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
                     call error_handler("MolMech:def_poligon_tes: this case (1) &
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
                     call error_handler("MolMech:def_poligon_tes: this case (2) &
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
                     call error_handler("MolMech:this case (0) &
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
                  if(sqrt(dot_product(v_buf,v_buf))<small_d) then
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
            ! if(i.eq.0) print*,'i eq 0'
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
         d_vv= dot_product(v_buf,poli%xyz_vertex(i,:)-poli%xyz_bound(i,1:3))

         if(d_vv <small_d .and. -d_vv<small_d) then
           v_buf=v_buf/sqrt(dot_product(v_buf,v_buf))*poli%r_bound(i,1)+&
                poli%xyz_bound(i,1:3)
         else if(dot_product(v_buf,0.5_r8_kind*(poli%xyz_vertex(i,:)+&
              poli%xyz_vertex(i1,:))-poli%xyz_bound(i,1:3))< small_d ) then
            v_buf= -0.5_r8_kind*poli%xyz_vertex(i,:) &
                 -0.5_r8_kind*poli%xyz_vertex(i1,:) &
                 +poli%xyz_bound(i,1:3)
            v_buf= v_buf/sqrt(dot_product(v_buf,v_buf))*poli%r_bound(i,1)&
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
      do i=1,poli%n_vertises
         if(.not. extra_point(i)) then
            xyz_poligon_centre=xyz_poligon_centre+ & !(case 1)
                 (poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:))
         endif
      enddo
      center_factor=sqrt(dot_product(xyz_poligon_centre,xyz_poligon_centre))
      xyz_poligon_centre=xyz_poligon_centre*(r_sphere(sphere_number))/ &
           center_factor+xyz_sphere(sphere_number,:)
      center_factor=center_factor/r_sphere(sphere_number)

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

      first_vertex: do i=1,poli%n_vertises
         j=1
         i1=poli%bounds(i,j)

         if(poli%n_sphere(i,j)/=0) then
            v_buf=xyz_sphere(poli%n_sphere(i,j),:)-xyz_sphere(sphere_number,:)
            d_rad=sqrt(dot_product(v_buf,v_buf))
            cos_gamma=dot_product(poli%xyz_vertex(i,:)-xyz_sphere(sphere_number,:),v_buf)/&
                 (r_sphere(sphere_number)*d_rad)

            alpha=acosnum(dot_product(poli%xyz_vertex(i,:)-poli%xyz_bound(i,3*j-2:3*j), &
                 poli%xyz_vertex(i1,:)-poli%xyz_bound(i,3*j-2:3*j))/poli%r_bound(i,j)**2)

            ! because of the new sorting of vertices, always with neigbour
            !   number 1
            !if(j==1) then
            costhn(i)=cos_gamma
            phin(i)=alpha
            !else
            !   costhn(i1)=cos_gamma
            !   phin(i1)=alpha
            !endif
            poligon_area=poligon_area+alpha*cos_gamma
         endif
      enddo first_vertex

      do i=1,poli%n_vertises
         i1=poli%bounds(i,1)
         vn(i,:)=poli%xyz_vertex(i,:)-poli%xyz_bound(i,1:3)
         v_buf=vector_product(vn(i,:), &
              poli%xyz_vertex(i1,:)-poli%xyz_bound(i,1:3))

         v1=vector_product(vn(i,:),v_buf)

         v1=v1/sqrt(dot_product(v1,v1))
         i1=poli%bounds(i,2)
         vns(i,:)=poli%xyz_vertex(i,:)-poli%xyz_bound(i,4:6)

         v_buf=vector_product(vns(i,:), &
              poli%xyz_vertex(i1,:)-poli%xyz_bound(i,4:6))

         v2=vector_product(vns(i,:),v_buf)

         v2=v2/sqrt(dot_product(v2,v2))
         if(.not. extra_point(i)) then
            beta=pi-acosnum(dot_product(v1,v2))
            poligon_area=poligon_area-beta
         endif

         !v1,v2 are the tn,tn+1 in Cossi96
         if (calc_gradients) then
            tn(i,:)=v1
            tns(i,:)=v2
         endif

      enddo

      poligon_area=(poligon_area+2.0_r8_kind*pi*max_tes_part)*r_sphere(sphere_number)**2

      if (poligon_area < min_area) then
        DPRINT 'area small/negative (',poligon_area,') , canceled'
        DPRINT 'at center',xyz_poligon_centre
        poligon_area=0.0_r8_kind
      endif

      if(calc_gradients .and. poligon_area>0.0_r8_kind ) then
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
         call error_handler("MolMech:acos of value >1 or <-1 demanded?!")
      else if(cos_t>1.0_r8_kind) then
!!$      if(cos_t>1.0_r8_kind) then
         acosnum=acos(1.0_r8_kind)
      else if(cos_t<-1.0_r8_kind) then
         acosnum=acos(-1.0_r8_kind)
      else
         acosnum=acos(cos_t)
      end if

    end function acosnum
    !------------------------------------------------------------

!---------------------------------------------------------------------------
    subroutine line_cut(p1,p2,cent1,r1,cent2,r2,ins1i,ins2i,cent0,h_case,ncut,pcut)
      ! idea: any point of the edge considered (center at origin)
      ! can be written as:
      ! r=(v1+lam (v2-v1))/|(v1+lam (v2-v1))| *r_bound
      ! with lam in [0,1]
      ! look for solutions of |r-d|=R, where d is the vector to the
      ! cutting sphere center and R is its radius
      !** End of interface *****************************************

      real(kind=r8_kind),intent(in),dimension(3) :: p1,p2,cent1,cent2,cent0
      real(kind=r8_kind),intent(in) :: r1,r2
      logical, intent(in) :: ins1i, ins2i, h_case
      real(kind=r8_kind),intent(out),optional :: pcut(6)
      integer (kind=i4_kind), intent(out) :: ncut
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
      if(dot_product(vd,vd)<small .and. .not. dot_product(v12,v12) < small)&
                colin_case=.true. !vectors parallel!

      !test if edge has angle larger then pi -> only for cutting edges
      if(.not.colin_case) then
        vd=cent1-cent0
        if(dot_product(vd,vd)>small .and. .not. h_case ) then
           vd=vector_product(p1-cent0,p2-cent0)
!          larger pi case for h_case
           if(h_case) vd=-vd
           if(dot_product(vd,0.5_r8_kind*(p1-cent1)+0.5_r8_kind*(p2-cent1))&
                <0.0_r8_kind ) then
              smaller_pi=.false.
              vb= -0.5_r8_kind*p1 -0.5_r8_kind*p2 +cent1
              vb= vb/sqrt(dot_product(vb,vb)) *r1
              v12= vb -v1
              vd=vb-(cent2-cent1)
              if(dot_product(vd,vd)<r2**2) then
                 ins2=.true.
              else
                 ins2=.false.
              endif
           endif
        endif
      endif

      if(colin_case) then
        v12=vector_product(p1-cent0,p2-cent0)
        v12=v12/sqrt(dot_product(v12,v12)) * r1
        ! opposite for h_case (convex bond)
!       if (h_case) v12=-v12
        lower_bound=-1.0_r8_kind
      endif

      vd=cent2-cent1
      ncut=0

      d=sqrt(dot_product(vd,vd))
      del=r1**2+d**2-r2**2

1000  v1v12=dot_product(v1,v12)
      v1vd=dot_product(v1,vd)
      v12vd=dot_product(v12,vd)
      v12v12=dot_product(v12,v12)
      v1v1=dot_product(v1,v1)
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
           vb=vb/(sqrt(dot_product(vb,vb)))
         else
           if(lam(i)**2<1.0_r8_kind) then
              vb=lam(i)*v1+sqrt(1.0_r8_kind-lam(i)**2)*v12
           else
              vb=lam(i)*v1
           endif
           vb=vb/r1
         endif
         ! test if not solution of - vd*pcut = d**2+r1**2-r2**2
!!$         if(abs(-2.0_r8_kind*dot_product(vd,vb)*r1-del)<sqrt(small) .and. &
         if(abs(2.0_r8_kind*dot_product(vd,vb)*r1-del)>=sqrt(small) .and. &
              abs(dot_product(vd,vb))>sqrt(small)) then
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
             vb=vb/(sqrt(dot_product(vb,vb)))
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
           v1=v1/sqrt(dot_product(v1,v1))*r1
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
!        vb= vb/dot_product(vb,vb) *r1
!        pcut(1:3)=cent1+vb
!        ncut=11
!     endif

    end subroutine line_cut
    !------------------------------------------------------------

    !------------------------------------------------------------
    function vector_product(vector1,vector2)
      real(kind=r8_kind), intent(in) :: vector1(3), vector2(3)
      real(kind=r8_kind) :: vector_product(3)

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

      real(kind=r8_kind), dimension(3) :: vc,vd,tnsg,tng,&
           ep,no,nos,vbuf,dno,dnos,dik,diks
      real(kind=r8_kind) :: dp1,dpnoep,dpnosep,dpnotns,dpnostn,mu,mus,  &
           dist1,dist2,help,help1
      real(kind=r8_kind), dimension(3) :: dP,dTc,dTcs
      integer(kind=i4_kind) :: i,j,l,i1,i2,k,status
      logical :: b_cut,bs_cut
      real (kind=r8_kind), parameter :: small=1.0e-11_r8_kind

      type(arrmat2), allocatable :: cgdP(:,:,:),cgdv(:,:,:),cgdvs(:,:,:)
      type(arrmat1), allocatable :: cgdomeg(:,:,:),cgdcosth(:,:,:),cgdphi(:,:,:)

      !bound is the edge and vertex the corner of the poligon
      !all derivatives relative to sphere center
      !dP: vertex relative to sphere center
      !dv: vertex relative to bound center
      !dTc: bound center
      !tn.. bound tangents at vertices
      !phin: bound arc angle

      allocate(cgdP(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang0")
      allocate(cgdv(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang1")
      allocate(cgdvs(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang2")
      allocate(cgdomeg(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang3")
      allocate(cgdcosth(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang4")
      allocate(cgdphi(3,n_species,1),stat=status)
      if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang5")

      do i=1,n_species
         do j=1,3
            allocate(cgdP(j,i,1)%m(3,poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang6")
            allocate(cgdv(j,i,1)%m(3,poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang7")
            allocate(cgdvs(j,i,1)%m(3,poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang8")
            allocate(cgdomeg(j,i,1)%m(poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang9")
            allocate(cgdcosth(j,i,1)%m(poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang10")
            allocate(cgdphi(j,i,1)%m(poli%n_vertises),stat=status)
            if(status/=0) call error_handler("MolMech:alloc fail geom_grad_cutted_triang11")
            cgdP(j,i,1)%m=0.0_r8_kind
            cgdv(j,i,1)%m=0.0_r8_kind
            cgdvs(j,i,1)%m=0.0_r8_kind
            cgdcosth(j,i,1)%m=0.0_r8_kind
            cgdphi(j,i,1)%m=0.0_r8_kind
            cgdomeg(j,i,1)%m=0.0_r8_kind
         enddo
      enddo

      do i=1,poli%n_vertises
         !P here relative vertex vector
         ep=poli%xyz_vertex(i,:)-xyz_sphere(nsp,:)
         ep=ep/sqrt(dot_product(ep,ep))

         i1=poli%bounds(i,1)
         i2=poli%bounds(i,2)

         vbuf=vector_product(vn(i,:),vns(i1,:))
         no=vbuf/sqrt(dot_product(vbuf,vbuf))

         vbuf=vector_product(vns(i,:),vn(i2,:))
         nos=vbuf/sqrt(dot_product(vbuf,vbuf))

         b_cut=.false.
         if(poli%n_sphere(i,1)/=0) b_cut=.true.
         bs_cut=.false.
         if(poli%n_sphere(i,2)/=0) bs_cut=.true.

         if(b_cut) then
            dik=xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:)
            dist1=sqrt(dot_product(dik,dik))
            if(dot_product(no,dik)<0.0_r8_kind) no=-no
         endif
         if(bs_cut) then
            diks=xyz_sphere(poli%n_sphere(i,2),:)-xyz_sphere(nsp,:)
            dist2=sqrt(dot_product(diks,diks))
            if(dot_product(nos,diks)<0.0_r8_kind) nos=-nos
         endif

         ! compute products after prospective sign change!!!
         dpnotns=dot_product(no,tns(i,:))
         dpnostn=dot_product(nos,tn(i,:))
         dpnoep =dot_product(no ,ep)
         dpnosep=dot_product(nos,ep)

         do l=1,n_species
            do j=1,3
               dTc(:)=0.0_r8_kind
               dTcs(:)=0.0_r8_kind
               dno=0.0_r8_kind
               dnos=0.0_r8_kind

               !calculate dno and dnos
               !calculate dTc and dTcs
               if(b_cut) then
                  do k=1,3
                     vbuf(k)= cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,1) -&
                          cagr%dc(k,nsp)%xyz_grad(j,l,1)
                  enddo

                  help1=dot_product(dik,vbuf)
                  dno=vbuf/dist1-help1* dik/dist1**3

                  help=(r_sphere(nsp)**2- &
                       r_sphere(poli%n_sphere(i,1))**2)/dist1**2
                  dTc(:)=0.5_r8_kind*(1.0_r8_kind+help)*vbuf + &
                       dik*(r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,1)-  &
                       r_sphere(poli%n_sphere(i,1)) *          &
                       cagr%dR(poli%n_sphere(i,1))%xyz_grad(j,l,1)  &
                       -help1*help)/dist1**2
               endif
               if(bs_cut) then
                  do k=1,3
                     vbuf(k)= cagr%dc(k,poli%n_sphere(i,2))%xyz_grad(j,l,1)-&
                          cagr%dc(k,nsp)%xyz_grad(j,l,1)
                  enddo

                  help1=dot_product(diks,vbuf)
                  dnos=vbuf/dist2-help1* diks/dist2**3

                  help=(r_sphere(nsp)**2- &
                       r_sphere(poli%n_sphere(i,2))**2)/dist2**2
                  dTcs(:)=0.5_r8_kind*(1.0_r8_kind+help)*vbuf + &
                       diks*(r_sphere(nsp)*cagr%dR(nsp)%xyz_grad(j,l,1)-  &
                       r_sphere(poli%n_sphere(i,2))*               &
                       cagr%dR(poli%n_sphere(i,2))%xyz_grad(j,l,1)  -   &
                       help1*help)/dist2**2
               endif

               if(.not.extra_point(i)) then
                  mus=(dot_product(no,dTc(:))-                          &
                       dot_product(dno,vn(i,:)) -                       &
                       dpnoep*cagr%dR(nsp)%xyz_grad(j,l,1))/dpnotns
                  mu=(dot_product(nos,dTcs(:))-                         &
                       dot_product(dnos,vns(i,:)) -                     &
                       dpnosep*cagr%dR(nsp)%xyz_grad(j,l,1))/dpnostn

                  dP(:)= ep * cagr%dR(nsp)%xyz_grad(j,l,1)+ &
                       tn(i,:) * mu + tns(i,:) *mus
               else
                  vbuf=vector_product(ep,tn(i,:))
                  mu=(dot_product(no,dTc(:))-&
                       dot_product(dno,vn(i,:)) - &
                       dpnoep*cagr%dR(nsp)%xyz_grad(j,l,1))/&
                       dot_product(no,vbuf)
                  dP(:)= ep * cagr%dR(nsp)%xyz_grad(j,l,1) + &
                       vbuf(:) * mu
               endif

               do k=1,3
                  cgdP(j,l,1)%m(k,i)=dP(k)
                  if(.not.extra_point(i)) then
                     cagr%dcenter(j,l,1)%m(k,N_total+1)=dP(k) +&
                          cagr%dcenter(j,l,1)%m(k,N_total+1)
                  endif

                  cgdv(j,l,1)%m(k,i)=dP(k)-dTc(k) !! use for recalculating dTc
                  cgdvs(j,l,1)%m(k,i)=dP(k)-dTcs(k)

               enddo
               if(costhn(i)*costhn(i)>small) then
                  do k=1,3
                     vd(k)= cagr%dc(k,poli%n_sphere(i,1))%xyz_grad(j,l,1)-&
                          cagr%dc(k,nsp)%xyz_grad(j,l,1)
                  enddo
                  cgdcosth(j,l,1)%m(i)= &
                       deriv_cos(poli%xyz_vertex(i,:)-xyz_sphere(nsp,:),&
                       xyz_sphere(poli%n_sphere(i,1),:)-xyz_sphere(nsp,:),&
                       dP,vd,.false.,.false.)
               else
                  cgdcosth(j,l,1)%m(i)= zero
               endif

            enddo
         enddo
      enddo

      do i=1,poli%n_vertises
         i1=poli%bounds(i,1)
         i2=poli%bounds(i,2)

         do l=1,n_species
            do j=1,3

               vc=cgdv(j,l,1)%m(:,i)
               vd=cgdvs(j,l,1)%m(:,i1)

               cgdphi(j,l,1)%m(i)=  deriv_cos(vn(i,:),vns(i1,:),vc,vd,.true.,.false.)

               ! expansion of the double vector product terms ax(bxc)=b(ac)-c(ab)
               ! and setting (vv) = r_bound**2
               ! derivative of not normalized tangent
               dp1=dot_product(vn(i,:),vns(i1,:))
               tng =vn(i,:)*(dot_product(vn(i,:),vd)+dot_product(vc,vns(i1,:))) &
                    + vc *  dp1                                   &
                    - vd * poli%r_bound(i,1)**2                  &
                    - 2.0_r8_kind * vns(i1,:) * dot_product(vc,vn(i,:))

               ! norm of double vector product vx(vxv*)
               help=poli%r_bound(i,1)**2 * &
                    (poli%r_bound(i,1)**4 - dp1**2)
               if(help<0.0_r8_kind) help=0.0_r8_kind
               help=sqrt(help)

               tng=tng/help
               tng=tng-tn(i,:)*dot_product(tn(i,:),tng)

               vc=cgdvs(j,l,1)%m(:,i)
               vd=cgdv(j,l,1)%m(:,i2)

               dp1=dot_product(vns(i,:),vn(i2,:))

               tnsg=vns(i,:)*(dot_product(vns(i,:),vd)+dot_product(vc,vn(i2,:))) &
                    + vc *  dp1                                  &
                    - vd * poli%r_bound(i,2)**2                  &
                    - 2.0_r8_kind *vn(i2,:)* dot_product(vc,vns(i,:))

               help=poli%r_bound(i,2)**2 * &
                    (poli%r_bound(i,2)**4 - dp1**2)
               if(help<0.0_r8_kind) help=0.0_r8_kind
               help=sqrt(help)
               tnsg=tnsg/help
               tnsg=tnsg-tns(i,:)*dot_product(tns(i,:),tnsg)

               cgdomeg(j,l,1)%m(i)=-deriv_cos(tn(i,:),tns(i,:),tng,tnsg,.true.,extra_point(i))

               if(i==1) cagr%darea(j,l,1)%m(N_total+1) = 2.0_r8_kind * &
                    cagr%dR(nsp)%xyz_grad(j,l,1) *area / r_sphere(nsp)
               cagr%darea(j,l,1)%m(N_total+1)= &
                    cagr%darea(j,l,1)%m(N_total+1) +&
                    (phin(i)* cgdcosth(j,l,1)%m(i)    + &
                    costhn(i)* cgdphi(j,l,1)%m(i)     - &
                    cgdomeg(j,l,1)%m(i)) *r_sphere(nsp)**2
            enddo
         enddo
      enddo

      !representative point before scaling to sphere surface
      ep=(xyz_cent-xyz_sphere(nsp,:))*cent_f

      help=sqrt(dot_product(ep,ep))

      do l=1,n_species
         do j=1,3
            do k=1,3
               vd(k)=cagr%dcenter(j,l,1)%m(k,N_total+1)
            enddo
            vc=vd/help-ep*dot_product(ep,vd)/help**3
            do k=1,3
               cagr%dcenter(j,l,1)%m(k,N_total+1)= vc(k)*r_sphere(nsp) + &
                    cagr%dc(k,nsp)%xyz_grad(j,l,1) + &
                    ep(k)/help*cagr%dR(nsp)%xyz_grad(j,l,1)
            enddo
         enddo
      enddo

      do i=1,n_species
         do j=1,3
            deallocate(cgdP(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang1")
            deallocate(cgdv(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang2")
            deallocate(cgdvs(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang3")
            deallocate(cgdomeg(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang4")
            deallocate(cgdcosth(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang5")
            deallocate(cgdphi(j,i,1)%m,stat=status)
            if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang6")
         enddo
      enddo

      deallocate(cgdP,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang7")
      deallocate(cgdv,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang8")
      deallocate(cgdvs,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang9")
      deallocate(cgdomeg,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang10")
      deallocate(cgdcosth,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang11")
      deallocate(cgdphi,stat=status)
      if(status/=0) call error_handler("MolMech:dealloc fail geom_grad_cutted_triang12")

    end subroutine geom_grad_cutted_triang
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine geom_grad_orig_triang(nsp,area,xyz_cent)
      !gradients of non cutted surface triangles
      !** End of interface *****************************************
      integer(kind=i4_kind) :: l,j,k
      integer(kind=i4_kind) , intent(in):: nsp
      real(kind=r8_kind) , intent(in):: area,xyz_cent(3)

      do l=1,n_species
         do j=1,3
            cagr%darea(j,l,1)%m(N_total+1) = 2.0_r8_kind * area * &
                 cagr%dR(nsp)%xyz_grad(j,l,1) / r_sphere(nsp)
            do k=1,3
               !also valid for weight_cent_sin and weight_cent_mass
               !movement only by sphere movement
               cagr%dcenter(j,l,1)%m(k,N_total+1) = &
                    cagr%dc(k,nsp)%xyz_grad(j,l,1) + &
                    (xyz_cent(k) - xyz_sphere(nsp,k))* &
                    cagr%dR(nsp)%xyz_grad(j,l,1)/r_sphere(nsp)
            enddo
         enddo
      enddo

    end subroutine geom_grad_orig_triang
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
      integer(i4_kind) :: l,j,k
!      integer(i4_kind) :: m,i
      integer(i4_kind) :: n,na,n1,n2!,n3,n4
      real(r8_kind) :: e_cent(3)
      real(r8_kind) :: f1,f2,df1dx(3),df2dx(3)!,df1dy(3),df2dy(3),d2f1dxdy(3,3),d2f2dxdy(3,3)
      real(r8_kind) :: s1(3),rj(3),ri(3),r5(3),r15(3),nn(3),nn2,mm2
      real(r8_kind) :: rji(3),drji,r1i(3),dr1i,drjidx(3,3)!,drjidy(3,3)
      real(r8_kind) :: dr1idx(3,3)!,dr1idy(3,3)
      real(r8_kind) :: Rsj,Rsi,cosa,sina,cosb,sinb,aji,bji,cji,dji,eji,fji
      real(r8_kind) :: dnn2dx(3),dmm2dx(3),ds1dx(3,3),dridx(3,3),drjdx(3,3)
!      real(r8_kind) :: dnn2dy(3),dmm2dy(3),ds1dy(3,3),dridy(3,3),drjdy(3,3)
!      real(r8_kind) :: d2nn2dxdy(3,3),d2mm2dxdy(3,3),d2s1dxdy(3,3,3)
!      real(r8_kind) :: d2ridxdy(3,3,3),d2rjdxdy(3,3,3),d2rjidxdy(3,3,3),d2r1idxdy(3,3,3)
      real(r8_kind) :: dr5dx(3,3),dRsidx(3),dRsjdx(3),rbuf!,rbuf1,rbuf2
!      real(r8_kind) :: dr5dy(3,3),dRsidy(3),dRsjdy(3)
      real(r8_kind) :: dr15dx(3,3)!,dr15dy(3,3),d2r15dxdy(3,3,3)
!      real(r8_kind) :: d2r5dxdy(3,3,3),d2Rsidxdy(3,3),d2Rsjdxdy(3,3)
      real(r8_kind) :: dajidx(3),dbjidx(3),dcjidx(3),ddjidx(3),dejidx(3),dfjidx(3)
!      real(r8_kind) :: dajidy(3),dbjidy(3),dcjidy(3),ddjidy(3),dejidy(3),dfjidy(3)
!      real(r8_kind) :: d2ajidxdy(3,3),d2bjidxdy(3,3),d2cjidxdy(3,3),d2djidxdy(3,3),d2ejidxdy(3,3),d2fjidxdy(3,3)
      real(r8_kind) :: nn21,nn22,mm21,mm22
      real(r8_kind) :: sigma,sigma_tot,dsigmadx(3)!,sigma_tot1,dsigmady(3),d2sigmadxdy(3,3),sigma_buf
      !------------ Executable code --------------------------------
      e_cent=(xyz_cent(:)-xyz_sphere(nsp,:))/r_sphere(nsp)

      mm21=mm11*mm11; mm22=mm12*mm12; nn21=nn11*nn11; nn22=nn12*nn12

      s1=xyz_cent
      rj=xyz_sphere(nsp,:)
      Rsj=r_sphere(nsp)

      do l=1,n_species
         do j=1,3
            do k=1,3
               cagr%dcenter(j,l,1)%m(k,N_total+1) = &
                    cagr%dc(k,nsp)%xyz_grad(j,l,1) + &
                    e_cent(k)*cagr%dR(nsp)%xyz_grad(j,l,1)
            enddo

            cagr%darea(j,l,1)%m(N_total+1)=0.0_r8_kind
            d_area:do n=1,N_spheres
!                  if(zero_area(n)) cycle d_area
               if(n == nsp) cycle d_area
               na=n
               if((slab_calc.or.lattice_calc) .and. n > N_atom_spheres) na=at_index(n)

               dRsjdx(j)=cagr%dR(nsp)%xyz_grad(j,l,1)

               Rsi=r_sphere(n)
               dRsidx(j)=cagr%dR(na)%xyz_grad(j,l,1)

               ri=xyz_sphere(n,:)
               rji=rj-ri
               drji=sqrt(dot_product(rji,rji))

               do n1=1,3
                  ds1dx(j,n1)=cagr%dcenter(j,l,1)%m(n1,N_total+1)
                  drjdx(j,n1)=cagr%dc(n1,nsp)%xyz_grad(j,l,1)
                  dridx(j,n1)=cagr%dc(n1,na)%xyz_grad(j,l,1)
               end do
               drjidx(j,:)=drjdx(j,:)-dridx(j,:)

               r5=ri+rji*Rsi/drji

               rbuf=dot_product(rji,drjidx(j,:))
               do n1=1,3
                  dr5dx(j,n1)=dridx(j,n1)+drjidx(j,n1)*Rsi/drji+ &
                       rji(j)*dRsidx(n1)/drji-rji(n1)*rbuf*Rsi/drji**3
               end do
               dr15dx(j,:)=ds1dx(j,:)-dr5dx(j,:)

               nn=s1-r5; r15=s1-r5
               nn2=dot_product(nn,nn)

               r1i=s1-ri
               dr1i=sqrt(dot_product(r1i,r1i))
               dr1idx(j,:)=ds1dx(j,:)-dridx(j,:)

               aji=Rsj*Rsj+drji*drji-dr1i*dr1i
               bji=2.0_r8_kind*Rsj*drji
               cji=Rsj*Rsj+drji*drji-Rsi*Rsi
               dji=sqrt(bji*bji-aji*aji)
               eji=sqrt(bji*bji-cji*cji)
               fji=1.0_r8_kind/(bji*bji)

               dajidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+dot_product(rji,drjidx(j,:))- &
                    dot_product(r1i,dr1idx(j,:)))
               dbjidx(j)=2.0_r8_kind*(dRsjdx(j)*drji+ &
                    Rsj*dot_product(rji,drjidx(j,:))/drji)
               dcjidx(j)=2.0_r8_kind*(Rsj*dRsjdx(j)+dot_product(rji,drjidx(j,:))- &
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

                  dnn2dx(j)=2.0_r8_kind*(dot_product(r15,dr15dx(j,:)))
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
                  drji=sqrt(dot_product(rji,rji))

                  r5=ri+rji*Rsi/drji

                  nn=s1-r5
                  nn2=dot_product(nn,nn)

                  r1i=s1-ri
                  dr1i=sqrt(dot_product(r1i,r1i))

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
               cagr%darea(j,l,1)%m(N_total+1)= &
                    cagr%darea(j,l,1)%m(N_total+1)+dsigmadx(j)*sigma_tot
            end do d_area
            cagr%darea(j,l,1)%m(N_total+1)=cagr%darea(j,l,1)%m(N_total+1)*area

         end do
      enddo

    end subroutine geom_grad_reduced_triang1
    !------------------------------------------------------------

    !------------------------------------------------------------
    function deriv_cos(a,b,da,db,do_acos,extra_p)

     real(kind=r8_kind), dimension(3), intent(in) :: a,b,da,db
     real(kind=r8_kind) :: deriv_cos
     real(kind=r8_kind) :: costh,absa,absb,vbuf(3),sinth
     logical, intent(in) :: do_acos, extra_p
     real(kind=r8_kind), parameter :: small=1.e-9
     ! if do_acos == .true. the derivative of the angle is computed
     ! else the derivative of the cosine of the angle is computed

     ! see formula (31) in Cossi 96

     absa= sqrt(dot_product(a,a))
     absb= sqrt(dot_product(b,b))
     costh= dot_product(a,b)
     if(absa>small .and. absb>small) then
        costh=costh/(absa*absb)
        if(do_acos) then
                if(costh**2<1.0_r8_kind) then
                   sinth = sqrt(1.0_r8_kind-costh**2)
                else
                   sinth = 0.0_r8_kind !sign(max(small,abs(sinth)),sinth) !0.0_r8_kind
                endif
                if(sinth<small .and. .not. extra_p) then
                   call error_handler("deriv_cos: angle to small")
                else if (sinth<small) then
                   deriv_cos=0.0_r8_kind
                   return
                endif
        endif
        vbuf(:)=a(:)-costh*absa/absb* b(:)
        deriv_cos=dot_product(vbuf,db)
        vbuf(:)=b(:)-costh*absb/absa* a(:)
        deriv_cos=deriv_cos+dot_product(vbuf,da)
        deriv_cos=deriv_cos/(absa*absb)
        if(do_acos) deriv_cos = - deriv_cos/sinth
     else
        deriv_cos=0.0_r8_kind
     endif

    end function deriv_cos
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine symm_sorted_centers()

      integer(kind=i4_kind) :: i,status


      n_size=n_total
      allocate(tessarea(n_size),stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation TESSAREA is failed")

      do i=1,n_size
         tessarea(i) % n_equal = 1
         tessarea(i)%area=area_tes(i)
         tessarea(i)%r_tes=r_tes(i)
         tessarea(i)%cut_off=cuttt(i)

         allocate(tessarea(i)%xyz(1,3), &
              tessarea(i)%sphere(1), stat=status)

         tessarea(i)%xyz(1,:) = xyz_tes_c(i,:)
         tessarea(i)%sphere(1) = sphere(i)
      enddo

!!$      if(output_cavity_data) then
!!$         write (output_unit,*) '--------------------------------------------------------------'
!!$         write (output_unit,*) &
!!$              'The symmetrization of surface points has been done succesfully, number of unique points:', &
!!$              n_size
!!$         if(output_cavity_long) then
!!$           write(output_unit,*)'number n_equal                  coordinates(a.u.)               area'
!!$           do i=1,n_size
!!$            write(output_unit,'(1x,i4,2x,i4,3x,3(f13.9,1x),f14.9)') i,tessarea(i)%n_equal, &
!!$                 tessarea(i)%xyz(1,:),tessarea(i)%area
!!$           enddo
!!$          endif
!!$         write(output_unit,*) '---------------------------------------------------------------'
!!$      endif

      do i=1,n_size
         if(tessarea(i)%area <= 0.0_r8_kind)  then
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
           "MolMech:points_on_cavity_surface: deallocation xyz_tes_c, area_tes, r_tes are failed")

    end subroutine symm_sorted_centers
    !--------------------------------------------------------

    real(r8_kind) function norm(x)
      implicit none
      real(r8_kind), intent(in) :: x(3)
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

    if(.not.associated(cagr%dR) .and. .not.associated(cagr%dc)) return
    do i=1,N_spheres !size(cagr%dR)
       deallocate(cagr%dR(i)%xyz_grad,&
            stat=alloc_stat)
       if ( alloc_stat /= 0) call error_handler( &
            "MolMech:dealloc_geom_deriv_part: deallocation of xyz_grad is failed")
       do j=1,3
          deallocate(cagr%dc(j,i)%xyz_grad,&
               stat=alloc_stat)
          if ( alloc_stat /= 0) call error_handler( &
               "MolMech:dealloc_geom_deriv_part: deallocation of xyz_grad is failed")
       enddo
    enddo
    deallocate(cagr%dR,stat=alloc_stat)
    if ( alloc_stat /= 0) call error_handler( &
         "MolMech:dealloc_geom_deriv_part: deallocation of grad_atomic_center failed")
    deallocate(cagr%dc,stat=alloc_stat)
    if ( alloc_stat /= 0) call error_handler( &
         "MolMech:dealloc_geom_deriv_part: deallocation of grad_atomic_center failed")

  end subroutine dealloc_geom_deriv_part1
  !************************************************************

  !******************************************************
  subroutine dealloc_geom_deriv_part2
   !dealloc darea,dcenter of cagr
   !** End of interface *****************************************

    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: i,j

    do i=1,n_species
       do j=1,3
          deallocate(cagr%dcenter(j,i,1)%m, &
               cagr%darea(j,i,1)%m,stat=status)
          if ( status /= 0) call error_handler( &
               "MolMech:dealloc_geom_deriv_part2: deallocation of cagr%dcenter%m is failed")
       enddo
    enddo

    deallocate(cagr%dcenter,cagr%darea,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:dealloc_geom_deriv_part2: deallocation of cagr%dcenter is failed")

  end subroutine dealloc_geom_deriv_part2
  !*****************************************************

  !******************************************************
  subroutine dealloc_cavity_mm
   !deallocate tesserea
   !** End of interface *****************************************

    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: i

    do i=1,size(tessarea)
       deallocate(tessarea(i)%xyz, &
            tessarea(i)%sphere, stat=status)
       if (status .ne. 0 ) call error_handler( &
            "MolMech:dealloc_cavity: deallocation of TESSAREA%XYZ  failed")
    enddo
    deallocate(tessarea,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:dealloc_cavity: deallocation TESSAREA is failed")

    deallocate(r_sphere, stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:dealloc_cavity: deallocation r_sphere is failed")
!!$    if(allocated(iuniq))  then
!!$     deallocate(iuniq,stat=status)
!!$       if ( status /= 0) call error_handler( &
!!$            "MolMech:dealloc_cavity: deallocation of iuniq is failed")
!!$    endif

    deallocate(xyz_sphere, stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:dealloc_cavity: deallocation xyz_sphere is failed")

    if(allocated(parents)) then
       deallocate(parents, stat=status)
       if ( status /= 0) call error_handler( &
            "MolMech:dealloc_cavity: deallocation of parents is failed")
    endif

    if(allocated(at_index)) then
       deallocate(at_index,stat=status)
       if(status/=0) call error_handler( &
            "MolMech:dealloc_cavity: deallocation of at_index is failed")
    end if

    if(allocated(zero_area)) then
       deallocate(zero_area, stat=status)
       if ( status /= 0) call error_handler( &
            "MolMech:dealloc_cavity: deallocation of zero_area is failed")
    endif

  end subroutine dealloc_cavity_mm
  !******************************************************
end module cavity_module


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
module opt_data_module
  !---------------------------------------------------------------
  !
  !  Purpose:This module consists of main_data_module,
  !          coordinates_data_module and read_input subroutine.
  !          The reason of this merging is the fatal problems
  !          of compilation of  read_input.f90 file on DEC.
  !  References: ...
  !  Author: AS
  !  Date: 12/98
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use coortype_module
  use atom_data_module, only :dummy_charge, nuc_mass
  use iounitadmin_module
  use allocopt_module
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use efp_module, only: n_efp,efp_fixed,qm_fixed
#endif
#ifdef NEW_EPE
  use ewaldpc_module, only: epe_relaxation, qm_ref_run, opt_and_relax
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of constants and variables ------------
  integer(i4_kind),    public :: OPT_STDOUT = -1
#ifndef MAX_PATH
# define MAX_PATH 100
#endif
  type ,public :: epe_shells
     real(kind=r8_kind),dimension(3)::r
     real(kind=r8_kind),dimension(3)::s
     integer(kind=i4_kind)::k
     real(kind=r8_kind)::ant
  end type epe_shells

  integer(kind=i4_kind) :: num_int
  integer, parameter,  public :: filename_namelengthmax = MAX_PATH
  character(len=filename_namelengthmax), public :: opt_data_dir,opt_dir

  logical,public                           :: analitic_hessian_calculated=.false.

  !
  ! See a lengthy discussion at around the only place this variable is
  ! used. In short,  it should remain false by default  if we want the
  ! settings in optimizer.input to  have higher priority for backwards
  ! compatibility.   May/will   be  changed  after   reading  namelist
  ! operations from optimizer.input:
  !
  logical, private :: task_priority = .false.

!  use allocopt_module,only:analitic_hessian_calculated
  logical,public                           :: ts_scan = .false.
  logical,public                           :: tsscan_sphere = .false.
  logical,public                           :: select_sphere_vars = .false.
  logical,public                           :: qst_step = .false.
  logical,public                           :: tsscan_mix = .false.
  real(kind=r8_kind),public                :: rpmix=0.2_r8_kind
  logical,public                           :: exist_product = .false.
  logical,public                           :: linear_transit = .false.
  
  logical,public                           :: mass_center_print=.false.
  logical,public                           :: print_updated_geometry=.false.
  logical,public                           :: epe_interfaced_mode = .false. ! spec of epe impurities addeed
  logical,public                           :: stop_after_epe_relaxation=.true.
  logical,public                           :: epe_forces =.false.  ! short-range from EPE
  logical,public                           :: epeparameters =.false.  ! short-range from EPE
  logical,public                           :: save_epe_r =.false.
  logical,public:: list_epepar
  logical,public                           :: cart_format=.false.  ! Cartesian coordinates are used for optimization
  logical,public                           :: zmat_format=.true.   ! this means Vladimirs
  !                                                                ! gxfile-Format
  logical,public                           :: free_format=.false.  ! free format where internals
  !                                                                ! are specified in namelists
  logical,public                           :: valence_format=.false.
  logical,public                           :: keep_valence=.true.
  !
  logical,public                           :: atomic_units = .true. ! if set to .false. all 
  !                                                                 ! units will be taken as Angstrom
  !                                                                 ! and Degrees
  real(kind=r8_kind),public                :: alpha_valence=1.2_r8_kind
  ! 'zmat_format' and 'free_format' direct the way the connectivities are are defined.
  ! For 'zmat_format' only non-redundant internal coordinates can be defined, whereas
  ! for 'free_format' an arbirtrary set on redundant internal coordinates can be defined.
  ! If neither 'zmat_format' nor 'free_format' is set, the connectivities are taken from
  ! the list of covalent bondlengths in the valence_module, resulting usually in a 
  ! redundant set of (primitive) internals.
  logical, public                          :: line_search=.false., &
                                              estimate_grad=.false.,&
                                              estimate_hessian=.false.,line_search2=.false., &
                                              estimate_hessian_m = .false.,&
                                              force_directed=.false., &
                                              method_qn=.true.,&
                                              method_ah=.false.,&
                                              method_rfo=.false.,&
                                              method_gdiis=.false.,&
                                              quart=.false.,&
                                              dynamic=.true.
  logical,public                           :: step_restrict=.true.
  logical,public                           :: step_reset=.false.
  logical,public                           :: step_to_ts=.false.
  real(kind=r8_kind),public                :: step_max=0.3_r8_kind

  logical, public                          :: print_internals=.true.,&
                                              print_xmol=.false., &
                                              print_bmat=.false.,&
                                              print_gmat=.false.,&
                                              print_hesse=.false.,&
                                              print_debug=.false.,&
                                              convert_internals=.false., &
                                              calc_hessian=.false.,&
                                              calc_epeff_hessian=.false.,&
                                              optimization=.true., &
                                              ts_search=.false., &
                                              frequency_calculation=.false., &
                                              gx_test=.false., & ! switch to test gxfile before the actual
                                              ! start of the calculation
                                              calc_cart_hess=.false., & !!!!!!!!!AS
                                              calc_cart_grad=.false.    !!!!!!!!!AS
  !*** input switches for thermodynamic properties module ********************************
  integer(kind=i4_kind),public             :: symmetry_index=1              &
                                            , N_temperature_steps=1
  real(kind=r8_kind),public                :: temperature=273.15_r8_kind    &
                                            , delta_temperature=0.0_r8_kind &
                                            , pressure=100000.0_r8_kind

  ! input swtiches for the method of the Hessian update
  logical, public                          :: update_bfgs=.true.,&
                                              update_dfp=.true.,&
                                              update_fromcartessian=.false., &
                                              update_powell=.false.,&
                                              update_bofill=.false.,&
                                              update_direct=.true.,&
                                              update_CG    =.false.,  &
                                              estimate=.false.,&
                                              step_technique=.false.,&
                                              single_step=.false. ! make only one step per coordinate to
                                                                  ! calculate hessian
  ! integer input switches for method of Hessian handling
  real(kind=r8_kind),public                 :: step_size=0.01_r8_kind
  integer(kind=i4_kind),public              :: n_recalc=1000      ! set to a high value so that
  !                                                              ! in normal cases it will be calculated
  !                                                              ! only once.

  logical, public                          :: delocalized_coordinates=.false., &
                                              update_delocalized=.false.,&
                                              zmat_coordinates = .true., &
  ! delocalized_coordinates=t: use delocalized_coordinates internal coordinates after J.Baker
  ! zmat=t: use Vladimirs Z-Matrix coordinates
                                              cart_coordinates = .false.
  real(kind=r8_kind),public:: crossboundary_3b=0.0_r8_kind

  ! convergence criteria
  real(kind=r8_kind),public                :: rms_step = 1.0e-5_r8_kind, &
                                              rms_grad = 1.0e-5_r8_kind,&
                                              max_comp_step = 1.0e-5_r8_kind,&
                                              max_comp_grad = 1.0e-5_r8_kind,&
                                              max_dEdR_sphere = 1.0_r8_kind
  ! ts control parameters
  logical,public                           :: eigenmode_follow = .false.,&
                                              rfo_step=.false.,&
                                              wales=.true.,&
                                              minimization=.false.
  ! fixed orientation switches
  logical,public                           :: fixed_orientation=.false.
  integer(kind=i4_kind),public             :: fixed_atom_1=0_i4_kind,&
                                              fixed_atom_2=0_i4_kind,&
                                              fixed_atom_3=0_i4_kind

  logical,public :: logic_coor_map(500,3)

  namelist /coordinates/ delocalized_coordinates,update_delocalized,zmat_coordinates,cart_coordinates
  namelist /input_format/ cart_format,zmat_format,free_format,valence_format, &
#ifdef NEW_EPE
       atomic_units,alpha_valence,keep_valence, epe_interfaced_mode
#else
       atomic_units,alpha_valence,keep_valence, &
       epe_interfaced_mode
#endif
  namelist /method/ line_search, estimate_grad, &
       line_search2,method_qn,method_ah,method_rfo,method_gdiis,quart,dynamic, &
       force_directed
  namelist /hesse_method/ update_bfgs, update_dfp, update_direct, update_cg,&
       update_powell, update_bofill, single_step, estimate, n_recalc,&
       step_size,step_technique, update_fromcartessian
  namelist /step_parameters/ step_restrict, step_max, step_reset,step_to_ts
  namelist /operations/ print_internals,print_xmol,print_bmat,&
                        print_debug,&
                        print_gmat,&
                        print_hesse,&
                        convert_internals,&
                        calc_hessian,&
!#ifndef NEW_EPE
                        calc_epeff_hessian,&
!#endif
                        optimization,estimate_hessian,&
                        estimate_hessian_m,&
                        ts_search,&
                        frequency_calculation, &
                        gx_test, &
!#ifndef NEW_EPE
                        epe_forces, save_epe_r, &
                        epeparameters,list_epepar, &
                        stop_after_epe_relaxation, &
!#endif
                        calc_cart_hess,& !!!!!!!!!!AS
                        calc_cart_grad,& !!!!!!!!!!AS
                        analitic_hessian_calculated, &
                        mass_center_print,print_updated_geometry, &
                        task_priority,ts_scan,tsscan_sphere,tsscan_mix, qst_step, &
                        rpmix, linear_transit, &
                        symmetry_index,temperature,pressure, &
                        delta_temperature,N_temperature_steps
  namelist /convergence/ rms_step, max_comp_step, rms_grad, max_comp_grad, max_dEdR_sphere
  namelist /ts_parameters/ eigenmode_follow, rfo_step, wales, minimization
  namelist /orientation/ fixed_orientation,fixed_atom_1,fixed_atom_2,fixed_atom_3  

  ! namelists dont need to be public:
  private :: input_format,method,coordinates,operations,&
       hesse_method,step_parameters,convergence,ts_parameters,&
       orientation

  !------------ public functions and subroutines ------------------

  public :: filename_setup_opt
  public :: opt_read_input

  !===============================================================================
  integer(kind=i4_kind),public             :: io_flepo
  integer(kind=i4_kind),public             :: n_atoms,&      ! number of atoms
                                              n_internal,&   ! number of all internals
                                              n_primitive,&  ! primitive (Zmat) internals
                                              n_dummy, &     ! number of dummy atoms
                                              n_tot_atoms    ! number of atoms + dummies
  ! these are help variables to read the gxfile input -----------------------------
  integer(kind=i4_kind),parameter,public   :: max_atoms=999,max_equal=100
  integer(kind=i4_kind),public             :: zmat(3,max_atoms),numx(3,max_atoms), &
                                              index_unique(max_atoms), &
                                              index_eq(max_atoms), &
                                              impu(max_atoms)
  ! -------------------------------------------------------------------------------

  ! cartesian coordinates ---------------------------------------------------------
  type(atom_type),pointer,public       :: atom(:)
  type(atom_type),pointer,public       :: atom_reactant(:),atom_product(:)

  real(kind=r8_kind),public                :: &
       & x(max_atoms),y(max_atoms),z(max_atoms), & ! MUSTDIE!
       & xyz(3,max_atoms), &
       & xyz_reactant(3,max_atoms), xyz_product(3,max_atoms), &
!       & xyz_point_on_mep(3,max_atoms), &
       & charge(max_atoms)

#ifdef WITH_EFP
  real(r8_kind), allocatable, public :: xyz_torque(:,:)
  integer(i4_kind),public :: i_tq !file to keep torque displacements (not the real coordinates)
#endif
  logical,public :: dummy_list(max_atoms)
  ! -------------------------------------------------------------------------------

  ! variables for reading the input for frequency calculation
  integer(kind=i4_kind), public            :: n_defined_masses=0, & ! number of atoms 
       ! with user defined mases
                                              n_infinite_masses=0   ! number of atoms 
       ! with a infenitely large mass
  real(kind=r8_kind),allocatable, public   :: defined_masses(:) ! user defined masses
  integer(kind=i4_kind),allocatable,public :: defined_masses_index(:) ! indices of the
  ! atoms with user defined masses
  integer(kind=i4_kind),allocatable,public :: infinite_masses_index(:) ! indices of the 
  ! atoms with infinite masses
  logical, public                          :: calculate_intensities=.false. ! calculate
  ! intensities of  the vibrational  normal modes set  it to  false by
  ! default to allow task="freq_analyt" without optimizer.input FIXME:
  ! make  it  dependent  on   DIPOLE=true/false  or  on  existence  of
  ! "dipole.dat"

  namelist /frequency_parameters/ n_defined_masses, n_infinite_masses, calculate_intensities
  ! internal coordinates ----------------------------------------------------------
  integer(kind=i4_kind),public :: sphere_dependent_var
  real(kind=r8_kind),public :: distance_to_reactant,distance_to_product,tsscan_rp_var
  real(kind=r8_kind),public :: distance_to_ts
  type(int_coor),pointer,public :: s(:),s_prim(:)
  type(int_coor),pointer,public :: s_reactant(:),s_prim_reactant(:)
  type(int_coor),pointer,public :: s_pointonmep(:),s_prim_pointonmep(:)
  type(int_coor),pointer,public :: s_product(:),s_prim_product(:)
  real(kind=r8_kind),pointer,public        :: q(:),q_prim(:),q_prim_pointonmep(:)
  real(kind=r8_kind),pointer,public        :: q_prim_reactant(:),q_prim_product(:)
  real(kind=r8_kind),allocatable,public    :: new_internal(:)
  integer(kind=i4_kind),allocatable,public :: sym_int(:)
  ! -------------------------------------------------------------------------------


  ! variables to read the free format ---------------------------------------------
  integer(kind=i4_kind),public              :: n_stretch,n_bend,n_torsion
  integer(kind=i4_kind),allocatable,public  :: stretch(:,:),bend(:,:),torsion(:,:)
  logical,allocatable,public                :: var_stretch(:),var_bend(:),var_tors(:)
  ! stretch(i_stretch,1) = partner1, stretch(i_stretch,2) = partner2
  ! bend(i_bend,1) = partner1, bend(i_bend,2) = partner2, bend(i_bend,3) = apex
  ! torsion(i_tor,1) = partner1, torsion(i_tor,2) = partner2, 
  ! torsion(i_tor,3) = base1, torsion(i_tor,4) = base2
  ! -------------------------------------------------------------------------------

  ! constraint variables for the valence_format -----------------------------------
  integer(kind=i4_kind),allocatable,public  :: const_stretch(:,:),const_bend(:,:),&
                                               const_torsion(:,:)
  ! const_stretch(i_stretch,1) = partner1, const_stretch(i_stretch,2)=partner2
  ! where 'i_stretch' is the index of the bond length to be fixed and 
  ! partner1 and partner2 the atoms that make up the bond.
  ! dito for const_bend and const_torsion.
  !--------------------------------------------------------------------------------
  ! n_internal     : number of internal coordinates as specified 
  !                  in the input
  ! q              : values of actual internal coordinates, This is a pointer
  !                  to the variable 's(:)%value'.
  ! dummy_list     : true if corredponding atom is a DUMMY atom,
  !                  false otherwise
  ! index_unique   : array containing a number which corresponds
  !                  to the type of unique atom
  ! index_eq       : index of atom
  ! charge         : charge(i) contains the charge of atom i
  ! s              : data structure containing the types, involved atoms
  !                  etc. for the internal coordinates ( values
  !                  should be the same as in variable 'q')
  ! atom           : data structure containing the actual positions
  !                  of the atoms in cartesian coordinates.


  !------------ public functions and subroutines ------------------
  public coordinates_read,coordinates_write,read_internals, frequency_parameters
  !================================================================
  ! End of public interface of module
  !================================================================
  

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !---------------------------------------------------------------
  !
  !  Purpose: Contains Hesse-Matrix  (in internal coordinates)
  !           Routines for input and output of the Hesse-Matrix,
  !           Transformation Routines ...
  !
  !  Module called by: main_opt ...
  !
  !  References: ...
  !  Author: FN
  !  Date: 6/97
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------

  subroutine filename_setup_opt(optonly)
    !
    ! Purpose: reads environment variables TTFSINPUTDIR, TTFSOUTPUTDIR
    !
    use filename_module, only: env=>filename_env, filename_set_input_dir
    implicit none
    logical, intent(in) :: optonly
    ! *** end of interface ***

    ! ask environment variable "TTFSINPUTDIR" that should describe
    ! directory where input data are originally stored
    if ( .not. env("TTFSINPUTDIR", opt_data_dir) ) then
       opt_data_dir = "."
    endif

    ! ask environment variable "TTFSOUTPUTDIR" that should describe
    ! directory where output data will finally be copied to
    if ( .not. env("TTFSOUTPUTDIR", opt_dir) ) then
       opt_dir = opt_data_dir
    endif

    if ( optonly ) then
       call filename_set_input_dir(opt_dir)
    endif
  end subroutine filename_setup_opt


  subroutine coordinates_read(io_gx,geo_loop,xyz,step_counter)
    !  Purpose: preliminary routine to read the gx.file
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind),intent(inout) :: io_gx
    real(kind=r8_kind),intent(inout)    :: geo_loop
    integer(kind=i4_kind),intent(in)    :: step_counter
    real(kind=r8_kind),intent(out):: xyz(:,:)

    integer(kind=i4_kind) :: io_gx_reactant

    !------------ Executable code --------------------------------

    if (zmat_format) then

       if(tsscan_sphere.or.tsscan_mix) then
         io_gx_reactant=openget_iounit(status='old',form='formatted', &
          file=trim(opt_data_dir)//'/gx.reactant') 
         call zmatformat_read(xyz_reactant,io_gx_reactant,geo_loop)
         call returnclose_iounit(io_gx_reactant)
         inquire(EXIST=exist_product,FILE=trim(opt_data_dir)//'/gx.product') 
       endif

#if 0
       if(tsscan_sphere) then
         io_gx_reactant=openget_iounit(status='old',form='formatted', &
          file=trim(opt_data_dir)//'/gx.point_on_mep') 
         call zmatformat_read(xyz_point_on_mep,io_gx_reactant,geo_loop)
         call returnclose_iounit(io_gx_reactant)
       endif
#endif

       if(tsscan_mix.or.(tsscan_sphere.and.exist_product)) then
         io_gx_reactant=openget_iounit(status='old',form='formatted', &
          file=trim(opt_data_dir)//'/gx.product') 
        call zmatformat_read(xyz_product,io_gx_reactant,geo_loop)
        call returnclose_iounit(io_gx_reactant)
       endif

!     if(linear_transit) then
!      write(char_step_counter,'(i1)') step_counter+1
!      DPRINT char_step_counter
!      io_gx=openget_iounit(status='old',form='formatted', &
!          file=trim(opt_data_dir)//'/gx.'//char_step_counter)
!     endif

       call zmatformat_read(xyz,io_gx,geo_loop)   ! (1)

    elseif (free_format) then
       call freeformat_read(io_gx,geo_loop)
    elseif (cart_format) then
       call cartformat_read(io_gx,geo_loop)
    else
       call valenceformat_read(io_gx,geo_loop)
    endif

    ! set N_INTERNAL here for the case that internal coordinates
    ! are generated in the 'blanket-approach', n_primitive will then
    ! be set by the routine 'val_types'.
    if (.not.zmat_format.and..not.free_format.and..not.cart_format) then
       if ( (n_atoms+n_dummy) > 2_i4_kind ) then
          n_internal = 3_i4_kind*(n_atoms+n_dummy) - 6_i4_kind
       else
          n_internal = 3_i4_kind*(n_atoms+n_dummy) - 5_i4_kind
       endif
    endif
    n_tot_atoms=n_atoms+n_dummy
!    call set_masses()

  end subroutine coordinates_read

  subroutine set_masses()
    !  Purpose: set the masses for every atom
    !           Up to now only the masses from atomic_data_module
    !------------ Modules used ----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i
    !------------ Executable code --------------------------------
    print *,'n_tot_atoms:'
    do i=1,n_tot_atoms
       print *,'i',i
       print *,'charge:',charge(i)
       print *,'mass:',nuc_mass(int(charge(i),i4_kind))
       atom(i)%mass=nuc_mass(int(charge(i),i4_kind))
    end do

  end subroutine set_masses

  subroutine coordinates_write(io_gx,geo_loop,rp_mix)
    !  Purpose: preliminary routine to write the updated
    !           geometry to the gx.file
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    logical, optional, intent(in) :: rp_mix
    integer(kind=i4_kind)            :: io_gx
    real(kind=r8_kind)               :: geo_loop
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------


    if (zmat_format) then
    if(present(rp_mix)) then
       call zmatformat_write(io_gx,geo_loop,rp_mix)
    else
       call zmatformat_write(io_gx,geo_loop,.false.)
    endif
    elseif ( free_format ) then
       call freeformat_write(io_gx,geo_loop)
    elseif ( cart_format ) then
       call cartformat_write(io_gx,geo_loop)
    else
       call valenceformat_write(io_gx,geo_loop)
    endif
  
  end subroutine coordinates_write


  subroutine zmatformat_write(io_gx,geo_loop,rp_mix)
    ! purpose: write the gxfile in 'simol'format for updating
    !          a geometry step.
    ! ------------------------------------------------------------
    logical, intent(in):: rp_mix
    integer(kind=i4_kind),intent(in) :: io_gx
    real(kind=r8_kind),intent(in)    :: geo_loop
    integer(kind=i4_kind)            :: i,j,k,num_x(3)
    ! ------------------------------------------------------------
    real(kind=r8_kind):: fmix=0.2_r8_kind
    integer(kind=i4_kind):: ind
    fmix=rpmix


    ind=0
    
    output_loop: do i=1,n_atoms+n_dummy

        if(ts_scan.and.convert_internals) then
         num_x(:)=numx(:,i)
         do j=1,num_int
          do k=1,3
           if(sym_int(j).eq.numx(k,i)) num_x(k)=0
          enddo
         enddo
        endif

     if(rp_mix) then
      x(i)=xyz_reactant(1,i)*fmix+(1-fmix)*xyz_product(1,i)
      y(i)=xyz_reactant(2,i)*fmix+(1-fmix)*xyz_product(2,i)
      z(i)=xyz_reactant(3,i)*fmix+(1-fmix)*xyz_product(3,i)
     endif

!    always use proper orderning
     if(index_unique(i).ne.0)  then
      ind=ind+1
      index_unique(i)=ind
     endif
     index_eq(i)=i

       if(epe_interfaced_mode) then
         if(ts_scan.and.convert_internals) then
            write(io_gx,1001)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               num_x(1),num_x(2),num_x(3), impu(i)
         else
            write(io_gx,1001)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               numx(1,i),numx(2,i),numx(3,i),&
               impu(i)
         endif
       else
         if(ts_scan.and.convert_internals) then
            write(io_gx,1001)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               num_x(1),num_x(2),num_x(3)
         else
          write(io_gx,1000)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               numx(1,i),numx(2,i),numx(3,i)
         endif
       end if
    enddo output_loop

    if(epe_interfaced_mode)then
       write(io_gx,'(f6.1,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)') &
              REAL(geo_loop,KIND=r8_kind) , &
              0.0,0.0,0.0, 0,0,0,0,0,0,0,0,0
    else
       write(io_gx,'(f6.1,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)') &
              REAL(geo_loop,KIND=r8_kind) , &
              0.0,0.0,0.0, 0,0,0,0,0,0,0,0,0
    endif
    ! ------------------------------------
1000 format((f5.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4)) ! simol format
1001 format((f5.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)) ! simol format

  end subroutine zmatformat_write

  subroutine zmatformat_read(xyz,io_gx,geo_loop)
    !  Purpose: reads coordinates, z-matrix and definition of
    !           variable internal coordinates (numx) in
    ! ------------------------------------------------------------
    use gxfile, only: gxfile_read
    implicit none
    integer(kind=i4_kind),intent(inout) :: io_gx
    real(kind=r8_kind),intent(inout) :: geo_loop
    real(kind=r8_kind), intent(out):: xyz(:,:)

    ! --- Declaration of localc variables ------------------------
    integer(kind=i4_kind)   :: n_centers
    integer(kind=i4_kind)   :: io_ts_scan
    real(kind=r8_kind)::    xyz_dummy(3,max_atoms),geo_loop_dummy

    !------------ Executable code --------------------------------
    ! first initialize the input variables
    DPRINT "zmatformat_read: entered",io_gx
    x=0.0_r4_kind
    y=0.0_r4_kind
    z=0.0_r4_kind
    charge=0.0_r8_kind
    zmat=0_i4_kind
    numx=0_i4_kind
    index_unique=0_i4_kind
    index_eq=0_i4_kind
    n_atoms=0_i4_kind
    n_dummy=0_i4_kind
    dummy_list = .false.

    if(epe_interfaced_mode) then

     if(ts_scan.and.convert_internals) then
     call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
         & index_unique, index_eq, zmat, numx,iepe=impu)

     io_ts_scan = openget_iounit(status='old',form='formatted', &
                                 file=trim(opt_data_dir)//'/gx.ts_scan')
     call gxfile_read( io_ts_scan, n_centers, geo_loop_dummy, charge, xyz_dummy, &
         & index_unique, index_eq, zmat, numx,iepe=impu)
     call returnclose_iounit(io_ts_scan)
   else
     call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
         & index_unique, index_eq, zmat, numx,iepe=impu,calc_epeff_hessian=calc_epeff_hessian)
   endif
     

    else
     if(ts_scan.and.convert_internals) then
     call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
         & index_unique, index_eq, zmat, numx)

     io_ts_scan = openget_iounit(status='old',form='formatted', &
                                 file=trim(opt_data_dir)//'/gx.ts_scan')
     call gxfile_read( io_ts_scan, n_centers, geo_loop_dummy, charge, xyz_dummy, &
         & index_unique, index_eq, zmat, numx)
     call returnclose_iounit(io_ts_scan)
   else
     call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
         & index_unique, index_eq, zmat, numx,calc_epeff_hessian=calc_epeff_hessian)
   endif
  endif

    ! MUSTDIE:
    x(:n_centers) = xyz(1,:n_centers)
    y(:n_centers) = xyz(2,:n_centers)
    z(:n_centers) = xyz(3,:n_centers)

    dummy_list(:n_centers) = (charge(:n_centers).eq.dummy_charge)
    n_dummy = count(dummy_list(:n_centers))

    n_atoms = n_centers - n_dummy

    if (n_atoms == 0) then
       stop 'zmatformat_read: the gxfile was empty'
    endif
    
    ! now set the total number  of internal coordinates (n_internal)
    if (zmat_coordinates) then
       if ((n_atoms+n_dummy) == 2 ) then
          n_internal = 1
          n_primitive = n_internal
       else
          n_internal = maxval(reshape(abs(numx),(/3*max_atoms/)))
          n_primitive = 3_i4_kind*(n_atoms+n_dummy) - 6_i4_kind
       endif

    elseif( delocalized_coordinates) then
       if (n_atoms == 2 ) then
          n_internal = 1
       else
          n_internal = 3_i4_kind*(n_atoms+n_dummy) - 6_i4_kind
       endif
       n_primitive=n_internal
    endif

  end subroutine zmatformat_read

  subroutine cartformat_read(io_gx,geo_loop)
    !  Purpose: reads coordinates
    ! ------------------------------------------------------------
    use gxfile, only: gxfile_read
#ifdef WITH_EFP
    use iounitadmin_module, only: openget_iounit
    use filename_module, only: inpfile
#endif
    implicit none
    integer(kind=i4_kind),intent(inout) :: io_gx
    real(kind=r8_kind),intent(inout) :: geo_loop

    ! --- Declaration of localc variables ------------------------
    integer(kind=i4_kind)   :: n_centers
    logical :: read_cart
#ifdef WITH_EFP
    integer(kind=i4_kind)   :: status,i
    logical :: trq
#endif

    !------------ Executable code --------------------------------
    ! first initialize the input variables
    DPRINT "cartformat_read: entered",io_gx
    x=0.0_r4_kind
    y=0.0_r4_kind
    z=0.0_r4_kind
    charge=0.0_r8_kind
    zmat=0_i4_kind
    numx=0_i4_kind
    index_unique=0_i4_kind
    index_eq=0_i4_kind
    n_atoms=0_i4_kind
    n_dummy=0_i4_kind
    dummy_list = .false.
    read_cart=.true.

    if(epe_interfaced_mode) then
       call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
            & index_unique, index_eq, zmat, numx,iepe=impu,cart=read_cart)
    else
       call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
         & index_unique, index_eq, zmat, numx,cart=read_cart)
    endif

    ! MUSTDIE:
    x(:n_centers) = xyz(1,:n_centers)
    y(:n_centers) = xyz(2,:n_centers)
    z(:n_centers) = xyz(3,:n_centers)

    dummy_list(:n_centers) = (charge(:n_centers).eq.dummy_charge)
    n_dummy = count(dummy_list(:n_centers))
    ASSERT(n_dummy==0)

    n_atoms = n_centers - n_dummy
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) n_atoms=n_atoms-3*n_efp
#endif
    ASSERT(n_atoms >= 0)

    n_primitive=3_i4_kind*n_atoms
    n_internal=n_primitive

#ifdef WITH_EFP
    if(efp .and. n_efp > 0) then
       n_primitive=n_primitive+6*n_efp
       if(.not. efp_fixed) then
          if(.not.allocated(xyz_torque)) then
             allocate(xyz_torque(n_efp,3),stat=status)
             ASSERT(status==0)
          end if
          inquire(file=trim(inpfile('/torque_coor')),exist=trq)
          i_tq=openget_iounit(trim(inpfile('/torque_coor')), &
               form='formatted', status='unknown')
          if(abs(geo_loop)==1 .or. .not.trq) then
             xyz_torque=0.0_r8_kind
          else if(abs(geo_loop) > 1 .and. trq) then
             do i=1,n_efp
                read(i_tq,*) xyz_torque(i,:)
             end do
             rewind(i_tq)
          end if
       end if
       if(.not. efp_fixed) then
          if(qm_fixed) then
             n_internal=6*n_efp
          else
             n_internal=n_primitive
          end if
       end if
    end if
#endif

    ASSERT(n_primitive >= 0)

  end subroutine cartformat_read

  subroutine cartformat_write(io_gx,geo_loop)
    ! purpose: write the gxfile in 'simol'format for updating
    !          a geometry step.
    ! ------------------------------------------------------------
    integer(kind=i4_kind),intent(in) :: io_gx
    real(kind=r8_kind),intent(in)    :: geo_loop
    integer(kind=i4_kind)            :: i
    integer(kind=i4_kind)            :: n_particles
    ! ------------------------------------------------------------

    DPRINT 'cartformat_write: epe_interfaced_mode=',epe_interfaced_mode

    n_particles=n_atoms
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) n_particles=n_particles+3*n_efp
#endif
    output_loop: do i=1,n_particles+n_dummy

       if(epe_interfaced_mode) then
          write(io_gx,1001)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               numx(1,i),numx(2,i),numx(3,i),&
               impu(i)
       else
          write(io_gx,1000)charge(i),x(i),y(i),z(i),&
               index_unique(i),index_eq(i),&
               zmat(1,i),zmat(2,i),zmat(3,i),&
               numx(1,i),numx(2,i),numx(3,i)
       end if
    enddo output_loop

    if(epe_interfaced_mode)then
       write(io_gx,'(f6.1,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)') &
              REAL(geo_loop,KIND=r8_kind) , &
              0.0,0.0,0.0, 0,0,0,0,0,0,0,0,0
    else
!       write(io_gx,'(f5.1)') REAL(geo_loop,KIND=r8_kind)
       write(io_gx,'(f6.1,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)') &
              REAL(geo_loop,KIND=r8_kind) , &
              0.0,0.0,0.0, 0,0,0,0,0,0,0,0,0
    endif
    ! ------------------------------------
1000 format((f5.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4)) ! simol format
1001 format((f5.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4,i5)) ! simol format

  end subroutine cartformat_write

  subroutine freeformat_read(io_gx,geo_loop)
    ! Purpose: read in 
    !          I. list of cartesian coordinates
    !         II. a namlist containing the following items
    !             - total number of internal coordinates
    !             - number of stretches
    !                         bends
    !                         torsions
    !       III. n_stretches name lists 'STRETCH'
    !        IV. n_bends namelist 'BEND'
    !         V. n_torsions namelists 'TORSION'
    !
    ! ----------------------------------------------------------
    integer(kind=i4_kind),intent(inout)    :: io_gx
    real(kind=r8_kind),intent(inout)       :: geo_loop
    ! --- declaration of local variables -----------------------
    integer(kind=i4_kind)  :: partner1,partner2,apex,base1,base2
    integer(kind=i4_kind)  :: alloc_stat,stat,i,n_int
    logical                :: keep_constant 

    namelist /internals/ n_int, n_stretch, n_bend, n_torsion
    namelist /stretch_list/ partner1, partner2, keep_constant
    namelist /bend_list/ partner1,partner2,apex, keep_constant 
    namelist /torsion_list/ partner1,partner2,base1,base2, keep_constant 

    call read_coord(io_gx,geo_loop)

    read(io_gx,nml=internals)
    if (stat.gt.0) &
         stop 'freeformat_read: reading of namelist internals failed'
    ! error checking
    if(n_int/=(n_stretch+n_bend+n_torsion)) &
         stop ' freeformat_read: number of internals wrong'
    

    print*,'Number of internals = ',n_int
    print*,'Number of atoms =     ',n_atoms
    if (n_int < 3*(n_atoms+n_dummy)-6 ) then
       write(OPT_STDOUT,*)" You do not have specified enough internals to describe the molecule"
       write(OPT_STDOUT,*)" At least ",3*(n_atoms+n_dummy)-6," internal coordinates are required  "
       stop 1
    endif

    ! preliminary: this includes ONLY unconstrained optimization
    n_internal = 3_i4_kind*(n_atoms+n_dummy) - 6_i4_kind
    if ((n_atoms+n_dummy) <= 2_i4_kind) then
       n_internal = 3_i4_kind*(n_atoms+n_dummy) - 5_i4_kind
    endif
    n_primitive = n_int
    allocate(stretch(n_stretch,2),bend(n_bend,3),torsion(n_torsion,4),&
         var_stretch(n_stretch),var_bend(n_bend),var_tors(n_torsion),&
         STAT=alloc_stat)
    if(alloc_stat/=0) &
         stop ' freeformat_read: allocation (1) failed'
    
    stretch = 0_i4_kind
    bend = 0_i4_kind
    torsion = 0_i4_kind
    var_bend = .true.
    var_stretch=.true.
    var_tors=.true.

    do i=1,n_stretch
       keep_constant=.false.
       read(io_gx,nml=stretch_list,iostat=stat)
       if (stat.gt.0) then
          write(OPT_STDOUT,*) "freeformat_read: reading of namelist stretch_list failed",i
          stop
       endif
       stretch(i,1) = partner1
       stretch(i,2) = partner2
       if (keep_constant) var_stretch(i)=.false.
    enddo
    
    do i=1,n_bend
       keep_constant=.false.
       read(io_gx,nml=bend_list,iostat=stat)
       if (stat.gt.0) then
          write(OPT_STDOUT,*)"freeformat_read: reading of namelist bend_list failed",i
          stop 
       endif
       bend(i,1) = partner1
       bend(i,2) = partner2 
       bend(i,3) = apex
       if (keep_constant) var_bend(i)=.false.
    enddo

    do i=1,n_torsion
       keep_constant=.false.
       read(io_gx,nml=torsion_list,iostat=stat)
       if (stat.gt.0) then
            write(OPT_STDOUT,*)"freeformat_read: reading of namelist torsion_list failed",i
            stop
         endif
       torsion(i,1) = partner1
       torsion(i,2) = partner2 
       torsion(i,3) = base1
       torsion(i,4) = base2
       if (keep_constant) var_tors(i)=.false.
    enddo
    
  end subroutine freeformat_read
    
   !*************************************************************
  subroutine freeformat_write(io_gx,geo_loop)
    ! Purpose: write the gxfile in free format for updating a
    !          geometry step
    ! 
    ! Note: since I didnt want the namelists below to be global
    ! objects in this  module (because definitions of partner1,
    ! partner2 etc. do not make sense as global items) these namelists
    ! have to re-declraed for this routine.
    ! -------------------------------------------------------
    integer(kind=i4_kind),intent(in)   :: io_gx
    real(kind=r8_kind),intent(in)      :: geo_loop
    ! --- Declaration of local variables --------------------
    integer(kind=i4_kind)  :: partner1,partner2,apex,base1,base2
    integer(kind=i4_kind)  :: stat,i,n_int
    logical                :: keep_constant

    namelist /internals/ n_int, n_stretch, n_bend,n_torsion
    namelist /stretch_list/ partner1, partner2,  keep_constant
    namelist /bend_list/ partner1,partner2,apex, keep_constant 
    namelist /torsion_list/ partner1,partner2,base1,base2,  keep_constant

    output_loop_free: do i=1,n_atoms+n_dummy
       write(io_gx,1000)charge(i),x(i),&
                                  y(i),&
                                  z(i)
    enddo output_loop_free
1000 format(f5.2,3(2x,f21.12))
    ! this is only a preliminary measure
    print*,'Writing IWORK  = ',geo_loop,'to file'
    write(io_gx,*)geo_loop

    
    n_int = n_internal
    
    write(io_gx,nml=internals,iostat=stat)
    if (stat/=0) &
         stop 'writing of namelist internals falied'
    
    do i=1,n_stretch
       partner1 = stretch(i,1)
       partner2 = stretch(i,2)
       keep_constant = .not.var_stretch(i)
       write(io_gx,nml=stretch_list,iostat=stat)
       if (stat/=0) &
            stop 'freeformat_read: writing of namelist stretch_list failed'
    enddo

    do i=1,n_bend
       partner1 = bend(i,1)
       partner2 = bend(i,2)
       apex = bend(i,3)
       keep_constant = .not.var_bend(i)
       write(io_gx,nml=bend_list,iostat=stat)
       if (stat/=0) &
            stop 'freeformat_read: writing of namelist bend_list failed'
    enddo

    do i=1,n_torsion
       partner1 = torsion(i,1) 
       partner2 = torsion(i,2) 
       base1 = torsion(i,3)
       base2 = torsion(i,4)
       keep_constant = .not.var_tors(i)
       write(io_gx,nml=torsion_list,iostat=stat)
       if (stat/=0) &
            stop 'freeformat_read: writing of namelist torsion_list failed'
    enddo
  end subroutine freeformat_write

  !********************************************************************

  subroutine valenceformat_read(io_gx,geo_loop)
    ! Purpose: read in 
    !         I. the cartesian coordinates
    !        II. a namelist containing the following items
    !            - total number of constraints
    !            - number of bond-length constraints
    !            - number of bond-angle constraints
    !            - number of torsion constraints
    !      III.  n_stretch namelists 'STRETCH'
    !       IV.  n_bend namelists 'BEND'
    !        V.  n_torsion namelists 'TORSION'
    !      The contents of these namelists will be assigned to
    !      the variables 'const_stretch', 'const_angle' and
    !      'const_torsion' for further usage in val_setup.
    ! ----------------------------------------------------------
    integer(kind=i4_kind)    :: io_gx
    real(kind=r8_kind)       :: geo_loop
    ! --- declaration of local variables -----------------------
    integer(kind=i4_kind)  :: partner1,partner2,apex,base1,base2
    integer(kind=i4_kind)  :: alloc_stat,stat,i,n_int
    logical                :: keep_constant
    namelist /internals/ n_int, n_stretch, n_bend, n_torsion
    namelist /stretch_list/ partner1, partner2,keep_constant
    namelist /bend_list/ partner1,partner2,apex, keep_constant 
    namelist /torsion_list/ partner1,partner2,base1,base2, keep_constant  
   

    call read_coord(io_gx,geo_loop)
    
    n_int=0_i4_kind
    n_stretch=0_i4_kind
    n_bend=0_i4_kind
    n_torsion=0_i4_kind
    keep_constant=.true.

    read(io_gx,nml=internals,END=991,ERR=992)

    ! error checking
    if(n_int/=(n_stretch+n_bend+n_torsion)) call error_handler &
         (" freeformat_read: n_int /= n_stretch+n_bend+n_torsion ")
    
    if (n_int>=3_i4_kind*(n_atoms+n_dummy)-6_i4_kind) &
         call error_handler &
         ("valenceformat_read: you specified more constraints the degrees of&
         &freedom are available")
    allocate(const_stretch(n_stretch,2),STAT=alloc_stat)    
    if (alloc_stat/=0) call error_handler&
         ("valenceformat_read: allocation (2) failed")
    const_stretch=0_i4_kind
    allocate(const_bend(n_bend,3),STAT=alloc_stat)    
    if (alloc_stat/=0) call error_handler&
         ("valenceformat_read: allocation (2) failed")
    const_bend=0_i4_kind
    allocate(const_torsion(n_torsion,4),STAT=alloc_stat)    
    if (alloc_stat/=0) call error_handler&
         ("valenceformat_read: allocation (3) failed")
    const_torsion=0_i4_kind
    
    do i=1,n_stretch
       read(io_gx,nml=stretch_list,iostat=stat)
       if (stat.gt.0) then
          write(OPT_STDOUT,*) "valenceformat_read: reading of namelist stretch_list failed",i
          stop 1
       endif
       keep_constant=.true.
       const_stretch(i,1) = partner1
       const_stretch(i,2) = partner2
    enddo
    do i=1,n_bend
       read(io_gx,nml=bend_list,iostat=stat)
       if (stat.gt.0) then
          write(OPT_STDOUT,*)"valenceformat_read: reading of namelist bend_list failed",i
          stop 1
       endif
       keep_constant=.true.
       const_bend(i,1) = partner1
       const_bend(i,2) = partner2 
       const_bend(i,3) = apex
    enddo

    do i=1,n_torsion
       read(io_gx,nml=torsion_list,iostat=stat)
       if (stat.gt.0) then
          write(OPT_STDOUT,*)"valenceformat_read: reading of namelist torsion_list failed",i
          stop 1
       endif
       keep_constant=.true.
       const_torsion(i,1) = partner1
       const_torsion(i,2) = partner2 
       const_torsion(i,3) = base1
       const_torsion(i,4) = base2
    enddo
    return
991 call error_handler("valenceformat_read: end-of-file was reached")
992 call error_handler("valenceformat_read: error reading gxfile")
  end subroutine valenceformat_read

  !********************************************************************

  subroutine valenceformat_write(io_gx,geo_loop)
    ! Purpose: re-write the gxfile in the valenceformat. This is 
    !          the same as the freeformat. Only in order to provide
    !          greater flexibility for later changes has the been
    !          made a separate subroutine.
    !          For comments see freeformat_write
    ! ---------------------------------------------------------------
    integer(kind=i4_kind)   :: io_gx
    real(kind=r8_kind)      :: geo_loop
    ! --- Declaration of local variables --------------------
    integer(kind=i4_kind)  :: partner1,partner2,apex,base1,base2
    integer(kind=i4_kind)  :: stat,i,n_int
    logical                :: keep_constant   
    
    namelist /internals/ n_int, n_stretch, n_bend,n_torsion
    namelist /stretch_list/ partner1, partner2,  keep_constant
    namelist /bend_list/ partner1,partner2,apex, keep_constant 
    namelist /torsion_list/ partner1,partner2,base1,base2,  keep_constant  

    call write_coord(io_gx,geo_loop)
   
    n_int = n_stretch+n_bend+n_torsion
    write(io_gx,nml=internals,iostat=stat)
    if (stat/=0) call error_handler &
         ("writing of namelist internals failed")
    do i=1,n_stretch
       partner1 = const_stretch(i,1)
       partner2 = const_stretch(i,2)
       keep_constant = .true.
       write(io_gx,nml=stretch_list,iostat=stat)
       if (stat/=0) call error_handler &
            ("freeformat_read: writing of namelist stretch_list failed")
    enddo
    
    do i=1,n_bend
       partner1 = const_bend(i,1)
       partner2 = const_bend(i,2)
       apex = const_bend(i,3)
       keep_constant = .true.
       write(io_gx,nml=bend_list,iostat=stat)
       if (stat/=0) call error_handler &
            ("freeformat_read: writing of namelist bend_list failed")
    enddo

    do i=1,n_torsion
       partner1 = const_torsion(i,1) 
       partner2 = const_torsion(i,2) 
       base1 = const_torsion(i,3)
       base2 = const_torsion(i,4)
       keep_constant = .true.
       write(io_gx,nml=torsion_list,iostat=stat)
       if (stat/=0) call error_handler &
            ("freeformat_read: writing of namelist torsion_list failed")
    enddo
  end subroutine valenceformat_write

  !********************************************************************
  
  subroutine read_internals(io_unit)
    ! preliminary: this reads
    ! - the number of internal variables that have to be updated
    ! - the new list of internal variables
    !-----------------------------------------------
    integer(kind=i4_kind),intent(in)  :: io_unit
    !-----------------------------------------------
!    integer(kind=i4_kind) :: num_int, alloc_stat,i
    integer(kind=i4_kind) :: alloc_stat,i
    character(len=80) :: buffer
    integer(i4_kind) :: coor_type,ios
    
    read(io_unit,*) num_int
    allocate(new_internal(num_int),sym_int(num_int),STAT=alloc_stat)
    if (alloc_stat/=0) then
       write(OPT_STDOUT,*)"read_internals: allocation failed"
    endif
    DPRINT 'zmat_coordinates' ,zmat_coordinates
    if (zmat_coordinates) then
       do i=1,num_int
          read(io_unit,'(a80)') buffer
          read(buffer,*,iostat=ios)sym_int(i),new_internal(i),coor_type
          if(ios /= 0) then
             read(buffer,*,iostat=ios)sym_int(i),new_internal(i)
          else
             if(coor_type == 1) then
                new_internal(i)=new_internal(i)*1.889725988579_r8_kind
             elseif(coor_type == 2 .or. coor_type == 3) then
                new_internal(i)=new_internal(i)*3.141592653589793_r8_kind/180.0_r8_kind
             else
                call error_handler("Optimizer: wrong reading internals")
             end if
          end if
       enddo
    else
       do i=1,num_int
          read(io_unit,*)new_internal(i)
       enddo
    endif
       print*,'read_internals: ',num_int
    print*,new_internal
  end subroutine read_internals

  !*************************************************************

  subroutine  read_coord(io_gx,geo_loop)
    ! Purpose: read coordinates only
    ! --------------------------------------------------------
    use gxfile, only: gxfile_read
    implicit none
    integer(kind=i4_kind),intent(in)    :: io_gx
    real(kind=r8_kind),intent(out)      :: geo_loop
    ! --- declaration of local variables -----------------------

    integer(i4_kind) :: n_centers

    call gxfile_read( io_gx, n_centers, geo_loop, charge, xyz)

    ! MUSTDIE:
    x(:n_centers) = xyz(1,:n_centers)
    y(:n_centers) = xyz(2,:n_centers)
    z(:n_centers) = xyz(3,:n_centers)

    dummy_list(:n_centers) = (charge(:n_centers).eq.dummy_charge)
    n_dummy = count(dummy_list(:n_centers))

    n_atoms = n_centers - n_dummy

    if (n_atoms == 0) then
       stop 'read_coord: the gxfile was empty'
    endif

    ! check for fixed orientation switches
    if (fixed_orientation) then
       if (n_atoms+n_dummy < 3 ) then
          write(OPT_STDOUT,*)" read_coord: you specified a fixed orientation "
          write(OPT_STDOUT,*)"             with less than three atoms "
          call error_handler(" Rubbish!" )
       elseif (n_atoms+n_dummy == 3 ) then
          write(OPT_STDOUT,*)" read_coord : Dont  you think that fixing "
          write(OPT_STDOUT,*)" the orientation of three atoms is a bit too much?"
          write(OPT_STDOUT,*)" Anyway, you may go ahead and do this ... "
       endif
       if (maxval((/fixed_atom_1,fixed_atom_2,fixed_atom_3/)).gt. &
            n_atoms+n_dummy) then
          write(OPT_STDOUT,*)" read_coord: The atoms you want to fix are "
          write(OPT_STDOUT,*)"             not contained in your input   "
          call error_handler(" Rubbish!" ) 
       endif
    endif
  end subroutine read_coord
  !*************************************************************

  !*************************************************************
  subroutine  write_coord(io_gx,geo_loop)
    ! Purpose: write coordinates only
    ! --------------------------------------------------------
    integer(kind=i4_kind),intent(in)    :: io_gx
    real(kind=r8_kind),intent(inout)      :: geo_loop
    ! --- declaration of local variables -----------------------
    integer(kind=i4_kind)  :: i
    input_loop_free: do i=1,n_atoms+n_dummy
       write(io_gx,1000)charge(i),x(i),y(i),z(i)
    enddo input_loop_free
1000 format(f5.2,3(2X,f21.12)) ! simol format
    write(io_gx,*) geo_loop
  end subroutine write_coord
  !*************************************************************

subroutine opt_read_input(task, convert_internal)
  !----------------------------------------------------------------
  !  Purpose: read the input of the optimizer-program
  !
  !  Subroutine called by: main_opt
  !  Author: FN
  !  Date: 2/98
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
#ifndef FPP_OPTIMIZER
  use operations_module, only: operations_task, namelist_tasks_used, operations_qm_epe
  use solv_cavity_module, only: VTN
  use operations_module, only: operations_solvation_effect
#endif
!  use allocopt_module,only:analitic_hessian_calculated
  implicit none
  character(len=*), intent(in) :: task
  logical, optional, intent(in)::convert_internal 
  ! *** end of interface ***

  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------

  integer(kind=i4_kind)   :: status,io_optinput, i, alloc_stat, n,j
  logical                 :: exist,buffer(3)
  !------------ Executable code --------------------------------
  analitic_hessian_calculated=.false.

  DPRINT 'BEFORE call task_optimizer(',trim(task),')'
  DWRITE(*,nml=operations)

  ! set the internal controls according to the input argument task:
  call task_optimizer(task)

  DPRINT 'AFTER call task_optimizer(',trim(task),')'
  DWRITE(*,nml=operations)

#ifdef WITH_EFP
  if(efp) then
     delocalized_coordinates=.false.
     zmat_coordinates=.false.
     cart_coordinates=.true.
     cart_format=.true.
     zmat_format=.false.
     free_format=.false.
     calc_hessian=.false.
     calc_epeff_hessian=.false.
     estimate_hessian=.false.
     ts_search=.false.
     frequency_calculation=.false.
     gx_test=.false.
     print_internals=.false.
     convert_internals=.false.
     method_qn=.true.
     method_ah=.false.
     method_rfo=.false.
     line_search=.false.
     line_search2=.false.
     method_gdiis=.false.
     dynamic=.true.
#ifndef FPP_OPTIMIZER
     if(operations_solvation_effect) dynamic=.false.
#endif
     update_direct=.false.
     fixed_orientation=.false.
     rms_step = 1.0e-4_r8_kind
     rms_grad = 1.0e-4_r8_kind
     max_comp_step = 1.0e-4_r8_kind
     max_comp_grad = 1.0e-4_r8_kind
#ifndef FPP_OPTIMIZER
     if(operations_solvation_effect) then
        rms_step = 0.01_r8_kind
        max_comp_step = 0.1_r8_kind
     endif
#endif
  end if
#endif
  ! ..  then  check if  optimizer.input  exists  and  read it.   Thus,
  ! options set  in task_optimizer(task) will be  overwritten by those
  ! in  optimizer.input.   This should  be  default for  compatibility
  ! reasons.   To  ensure   migration  use  optimizer.input  only  for
  ! advanced options!
  io_optinput=30
  inquire(EXIST=exist,FILE=trim(opt_data_dir)//'/optimizer.input')
  if (.not.exist) then
     write(OPT_STDOUT,*) " read_input: no input file present - use defaults"
     return
  else
     DPRINT 'reading ',trim(opt_data_dir)//'/optimizer.input'
  endif

  io_optinput=openget_iounit(status='old',form='formatted',&
       file=trim(opt_data_dir)//'/optimizer.input')

  read(io_optinput,nml=coordinates,iostat=status)
#ifdef WITH_EFP
  if(efp) then
     delocalized_coordinates=.false.
     zmat_coordinates=.false.
     cart_coordinates=.true.
  end if
#endif
  if (status/=0) call error_handler("read_input: namelist coordinates failed")
  if (delocalized_coordinates.and.zmat_coordinates.and.cart_coordinates) then
     call error_handler&
          (" read_input: Choose either DELOCALIZED or ZMAT or CART coordinates")
  elseif (.not.delocalized_coordinates.and..not.zmat_coordinates.and..not.cart_coordinates) then
     call error_handler&
          (" read_input: Choose one of DELOCALIZED or ZMAT or CART coordinates")
  elseif ( zmat_coordinates .and. update_delocalized ) then
     write(OPT_STDOUT,*)" You have chosen ZMAT-COORDINATES and you set UPDATE_DELOCALIZED to true"
     write(OPT_STDOUT,*)" The UPDATE_DELOCALIZED swirch will be ignored"
  elseif(delocalized_coordinates.and.zmat_coordinates) then
     call error_handler&
          (" read_input: Choose either DELOCALIZED or ZMAT coordinates")
  elseif(zmat_coordinates.and.cart_coordinates) then
     call error_handler&
          (" read_input: Choose either ZMAT or CART coordinates")
  elseif(delocalized_coordinates.and.cart_coordinates) then
     call error_handler&
          (" read_input: Choose either DELOCALIZED or CART coordinates")
  endif

  read(io_optinput,nml=input_format,iostat=status)
  if (status/=0) call error_handler (" read_input: namelist input_format failed")
  if(cart_coordinates) then
     cart_format=.true.
     zmat_format=.false.
     free_format=.false.
  end if
 
  if ( count((/zmat_format,free_format,valence_format,cart_format/),1) /= 1) call error_handler&
       (" read_input: Please specify ALL variables in namelist INPUT_FORMAT correctly")

  read(io_optinput,nml=operations,iostat=status)
  if (status/=0) call error_handler (" read_input: namelist operations failed")
  DPRINT 'AFTER  call read(io_optinput,nml=operations,...)'
  DWRITE(*,nml=operations)
  DPRINT 'analitic_hessian_calculated', analitic_hessian_calculated

#ifndef FPP_OPTIMIZER
  !
  ! Do   not   call   task_optimizer(operations_task)  after   reading
  ! optimizer.input (see  above), the latter as a  legacy input method
  ! should   always  have  the   last  word.    You  may   still  call
  ! task_optimizer() with some subtask keyword.
  !
  ! This  was  reverted  by  the  patch  named  "NV:  TRANSIT  options
  ! [nasluzov] (3989)"  in 2009 which introduced  task_priority set to
  ! true  by default. With this change default was reset to false.
  !
  ! Calling task_optimizer() with operations_task from PG input breaks
  ! backwards  compatibility.   People  are  used to  specify  task  =
  ! "GeoOpt"  in  the PG  input,  but  set  more specific  options  in
  ! optimizer.input such  as for computing  numerical frequencies. For
  ! those users who believe they know  what they are doing and make an
  ! effort  to prepare  the  optimizer.input we  should respect  their
  ! settings.  Look above, we are calling task_optimizer() already.
  !
  ! Of course  for common tasks one  should encourage use  of the task
  ! keyword, such as  task = "freq_num".  That is  why we already pass
  ! the task keyword as an argument to optimizer and respect it here.
  !
  ! if(NEVER) call task_optimizer(operations_task)
  DPRINT 'namelist_tasks_used, task_priority =', namelist_tasks_used, task_priority
  if(namelist_tasks_used.and.task_priority) call task_optimizer(operations_task)
  !
#endif

  if(update_fromcartessian) then
   update_dfp=.false.
   update_bfgs=.false.
   update_bofill=.false.
   update_cg=.false.
  endif

  if(update_bofill) then
   update_dfp=.false.
   update_bfgs=.false.
   update_cg=.false.
  endif
  
  if(tsscan_mix) then
!   epe_forces=.false.
!   epeparameters=.false.
   ts_search=.false.
   calc_hessian=.false.
   ts_scan=.false.
  endif

  if (calc_cart_hess .or. calc_cart_grad) then

  ! calculation of the cartesian Hessian or gradients.
  ! Each atom is displaced in three directions ...
     if(calc_cart_hess .and. calc_cart_grad) &
          call error_handler("read_input: either calc_cart_hess or calc_cart_grad, not both")

     if (calc_cart_hess) write(OPT_STDOUT,*)"Calculation of Cartesian Hessian"
     calc_hessian=.false.
     calc_epeff_hessian=.false.
     optimization=.false.
     estimate_hessian=.false.
     ts_search=.false.
     frequency_calculation=.false.
     gx_test=.false.
     estimate_hessian_m=.false.
     save_epe_r=.false.
     zmat_format=.true.
     print_internals=.false.
     convert_internals=.false.
!VVVVVVVVVVVVVV
     if(calc_cart_grad) then
        logic_coor_map=.false.
        read(io_optinput,*,err=9999) n
        if(n < 0) n=0
        if(n/=0) then
           do i=1,n
              read(io_optinput,*,err=9999) j,buffer(:)
              logic_coor_map(j,:)=buffer
           end do
        goto 9998
        end if
9999    logic_coor_map=.true.
9998    continue
     end if
!^^^^^^^^^^^^^^^

     read(io_optinput,nml=hesse_method,iostat=status) !(1)
     goto 999 ! clean up and exit
     ! do not return, close io_optinput first!
  end if

  if (cart_format) then
     calc_hessian=.false.
     calc_epeff_hessian=.false.
     estimate_hessian=.false.
     ts_search=.false.
     frequency_calculation=.false.
     gx_test=.false.
!!!     print_internals=.false.
  end if

#ifdef NEW_EPE
  if(operations_qm_epe .and. (epe_relaxation .or. qm_ref_run)) optimization=.false.
  if(operations_qm_epe .and. opt_and_relax) optimization=.true.
#endif

  if (cart_format) convert_internals=.false.

  if (convert_internals.and.optimization ) call error_handler&
       (" read_input: convert_internals AND optimization not possible")
  if (optimization.and.ts_search ) then
     write(OPT_STDOUT,*)" read_input: Optimization AND transition state search not possible"
     write(OPT_STDOUT,*)"             If you want to minimize using the ts_search algorithm"
     write(OPT_STDOUT,*)"             please use the MINIMIZATION switch in the TS_PARAMETERS namelist"
     stop 1
  endif
  if (ts_search.and.convert_internals) call error_handler&
       (" read_input: convert_internals AND ts_search not possible")

  write(OPT_STDOUT,*)"    ---  Coordinates --- "
  write(OPT_STDOUT,*)" "
  if (delocalized_coordinates) then
     if (zmat_format) then
        write(OPT_STDOUT,*)" Using DELOCALIZED COORDINATES, connectivities specified in Z-Matrix Format"
     elseif(free_format) then
        write(OPT_STDOUT,*)" Using DELOCALIZED COORDINATES, connectivities specified in free format"
     else
        write(OPT_STDOUT,*)" Using DELOCALIZED COORDINATES, connectivities calculated based on colvalent radii"
        if(keep_valence) then
           write(OPT_STDOUT,*)" The primitive coordinates underlying the delocalized coordinates"
           write(OPT_STDOUT,*)" will be established in the first cycle only and subsequently kept"
           write(OPT_STDOUT,*)" in the file VAL_COOR.DAT"
        endif
     endif
  else if(zmat_coordinates) then
     write(OPT_STDOUT,*)" Using Z-MATRIX COORDINATES "
  else if(cart_coordinates) then
     write(OPT_STDOUT,*)" Using CARTESIAN COORDINATES "
  endif
  write(OPT_STDOUT,*)" "
  
  if(present(convert_internal)) then
   convert_internals=(convert_internal.and.ts_scan).or.convert_internals
   if(cart_format) convert_internals=.false.
    if(convert_internals) optimization=.false.
  endif
  if(convert_internals)  call read_internals(io_optinput)

  read(io_optinput,nml=method,iostat=status)
  if (status/=0) call error_handler(" read_input: reading of namelist METHOD failed")
  if(tsscan_mix) then
   update_fromcartessian=.true.
  endif
  if (cart_coordinates) then
     method_qn=.true.
     method_ah=.false.
     method_rfo=.false.
     line_search=.false.
     line_search2=.false.
     method_gdiis=.false.
     dynamic=.true.
#ifndef FPP_OPTIMIZER
     if(VTN .and. operations_solvation_effect) dynamic=.false.
#ifdef WITH_EFP
     if(efp .and. operations_solvation_effect) dynamic=.false.
#endif
#endif
  end if
  if (.not.line_search) then
     estimate_grad = .false. ! for consistency
  endif
  if (line_search2) then         !!!!!!!!!
     line_search = .false.       !!!!!!!!!
     estimate_grad = .false.     !!!!!!!!!
  end if                         !!!!!!!!!
  write(OPT_STDOUT,*)"    --- Optimization method --- "
  write(OPT_STDOUT,*)" "
  if ( line_search ) then
     if (estimate_grad) then
        write(OPT_STDOUT,*)" Linesearch (cubic_fit), gradient will be estimated for line search step"
     else
        write(OPT_STDOUT,*)" Linesearch (cubic_fit), gradient will be NOT estimated for line search step"
     endif
  elseif(line_search2) then
     write(OPT_STDOUT,*)" Linesearch (square_fit)"
  else
     write(OPT_STDOUT,*)" NO Linesearch "
  endif
  write(OPT_STDOUT,*)" "

  read(io_optinput,nml=hesse_method,iostat=status) !(2)
  if (status/=0) call error_handler(" read_input: namelist HESSE_METHOD failed")

  if(cart_coordinates) then
     update_direct=.false.
  end if

  if(frequency_calculation) then
   update_cg=.false.
  endif

  ! first check if input makes sense
   if (update_cg) then
      if( update_dfp.or.update_bfgs ) then
         write(*,*) "WARNING : CG and DFP/BFGS specified  : WARNING"
       write(*,*) "WARNING : Switching off DFP/BFGS     : WARNING"
       write(*,*) "WARNING : --- CG algorithm is ON --- : WARNING"
      end if
      update_dfp         = .false.
      update_bfgs        = .false.
      estimate_hessian   = .false.
      estimate_hessian_m = .false. ! (in principle "estimates" can be used)
      calc_hessian       = .false.
   end if

  if (.not.update_dfp.and.update_bfgs) call error_handler&
       (" read_input: choosing BFGS update without DFP not allowed")
  if (calc_hessian.and. .not.zmat_format) call error_handler&
       (" read_input: The hessian can be caclulated ONLY if zmat_format is chosen")
  if (n_recalc<0) call error_handler&
       (" read_input: N_RECALC must be positive")
  if (step_technique ) then
     if(.not.zmat_format) call error_handler&
          ("read_input: estimate of Hessian diagonals only with ZMAT_FORMAT possible")
     if(calc_hessian) call error_handler("read_input: either CALC_HESSIAN or STEP_TECHNIQUE")
     if (estimate_hessian) call error_handler("read_input: either ESTIMATE_HESSIAN or STEP_TECHNIQUE")
  endif
   if (estimate_hessian .and. estimate_hessian_m) &
      call error_handler("read_input: either ESTIMATE_HESSIAN or ESTIMATE_HESSIAN_M")
  if (calc_hessian ) then
     write(OPT_STDOUT,*)"   --- Initial calculation of Hessian --- "
     if (single_step) then
        write(OPT_STDOUT,*)"          single steps, i.e. number of steps required = number of internal variables"
     else

        write(OPT_STDOUT,*)"          double steps, i.e. number of steps required = 2 * number of internal variables"
     endif
     write(OPT_STDOUT,*)" "
  endif
  write(OPT_STDOUT,*) "    --- Update of  Hessian --- "
  write(OPT_STDOUT,*) " "
  if (.not.update_dfp.and..not.update_bfgs.and. &
      .not.update_powell.and..not.update_bofill .and. .not.update_cg ) then
     write(OPT_STDOUT,*)" The Hessian will not be updated. (Steepest descend method)"
  else
     if (update_direct .and. .not.update_CG) write(OPT_STDOUT,*)"  The Hessian will be updated directly "
     if (.not.update_direct .and. .not.update_CG)   write(OPT_STDOUT,*)"  The inverse Hessian will be updated"
     if (update_dfp .and.      update_bfgs)       write(OPT_STDOUT,*)"  using the BFGS-Algorithm"
     if (update_dfp .and. .not.update_bfgs)       write(OPT_STDOUT,*)"  using the DFP-Algorithm"

     if((.not.update_direct).and.update_powell) & 
          call error_handler("read_input: sorry, for Powell update only a direct update is implemented")
     if((.not.update_direct).and.update_bofill) &
          call error_handler("read_input: sorry, for Bofill update only a direct update is implemented")
     if(update_powell.and.(update_dfp.or.update_bfgs)) &
          call error_handler("read_input: sorry, using Powell update togehter with an other update does not make sense")

     if (update_powell.and..not.update_bofill) write(OPT_STDOUT,*)"  using the Powell-Algorithm"
     if(update_bofill) write(OPT_STDOUT,*)"  using the Bofill-Algorithm"
     if(update_bofill.and.(update_dfp.or.update_bfgs)) &
          call error_handler("read_input: sorry, using Bofill update togehter with an other update does not make sense")
  endif
  if (calc_hessian.and.n_recalc<1000_i4_kind) then
     write(OPT_STDOUT,*)"   --- Re-Calculation of Hessian --- "
     if (single_step) then
        write(OPT_STDOUT,*)"    Hessian will be re-calculated after",n_recalc," cycles"
        write(OPT_STDOUT,*)"    using the single-step strategy "
     else
        write(OPT_STDOUT,*)"    Hessian will be re-calculated after",n_recalc," cycles"
        write(OPT_STDOUT,*)"    using the double-step strategy"
     endif
  endif
  write(OPT_STDOUT,*)" "

  DPRINT 'read step_parameters'
  read(io_optinput,nml=step_parameters,iostat=status)
  if (status/=0) call error_handler ("read_input: namelist STEP_PARAMETERS failed")
  if (step_max<0.0_r8_kind) call error_handler("read_input: STEP_MAX negative")
  if (.not.optimization.and..not.ts_search) then
     write(OPT_STDOUT,*)"read_input: without OPTIMIZATION or TS_SEARCH switched on"
     write(OPT_STDOUT,*)"            step_parameters will be ignored"
  else
     write(OPT_STDOUT,*)"     --- Step Parameters --- "
     if (step_restrict) then
        write(OPT_STDOUT,'(" Step length will be limited to ",f9.5)')step_max
     else
        write(OPT_STDOUT,*)" Steps will not be limited "
     endif
     write(OPT_STDOUT,*)" "
  endif
  
  DPRINT 'nml=convergence'
  read(io_optinput,nml=convergence,iostat=status)
  if (status/=0) call error_handler ("read_input: namelist CONVERGENCE failed")
  if(cart_coordinates) then
     rms_step = 1.0e-4_r8_kind
     rms_grad = 1.0e-4_r8_kind
     max_comp_step = 1.0e-4_r8_kind
     max_comp_grad = 1.0e-4_r8_kind
#ifndef FPP_OPTIMIZER
     if(VTN .and. operations_solvation_effect) then
        rms_step = 0.001_r8_kind
        max_comp_step = 0.001_r8_kind
     endif
#ifdef WITH_EFP
     if(efp .and. operations_solvation_effect) then
        rms_step = 0.01_r8_kind
        max_comp_step = 0.01_r8_kind
     endif
#endif
#endif
  endif
  if (rms_step<0.0_r8_kind) call error_handler("read_input: rms_step must be positive")
  if (rms_grad<0.0_r8_kind) call error_handler("read_input: rms_grad must be positive")
  if (max_comp_step<0.0_r8_kind) call error_handler("read_input: max_comp_step must be positive")
  if (max_comp_grad<0.0_r8_kind) call error_handler("read_input: max_comp_grad must be positive")
  write(OPT_STDOUT,*)"     --- Convergence criteria ---"
  write(OPT_STDOUT,'("          rms_step: ",es12.5)')rms_step
  write(OPT_STDOUT,'("          rms_grad: ",es12.5)')rms_grad
  write(OPT_STDOUT,'("     max_comp_step: ",es12.5)')max_comp_step
  write(OPT_STDOUT,'("     max_comp_grad: ",es12.5)')max_comp_grad
  if(tsscan_sphere) write(OPT_STDOUT,'("     max_dEdR_sphere: ",es12.5)')max_dEdR_sphere

  read(io_optinput,nml=ts_parameters,iostat=status)
  if (status/=0) call error_handler("read_input: namelist TS_PARAMETERS failed")
  if (ts_search ) then
     write(OPT_STDOUT,*)"      --- TS-state search parameters ---"
     if (eigenmode_follow) then
        write(OPT_STDOUT,*)" Eigenmode following switched on. The lowest eigenmode"
        write(OPT_STDOUT,*)" of the initial Hessian will be followed, even if the update of "
        write(OPT_STDOUT,*)" produces other even lower lying modes. "
     else
        write(OPT_STDOUT,*)" Eigenmode following switched off. The lowest eigenmode"
        write(OPT_STDOUT,*)" of the actual Hessian will be followed even if the update"
        write(OPT_STDOUT,*)" has produced an artifial very low lying mode."
     endif
     if (rfo_step) then
        write(OPT_STDOUT,*)" If the Hessian has already one negative eigenvalue, the shift"
        write(OPT_STDOUT,*)" parameters for the RFO-Step will nevertheless be imposed"
     else
        write(OPT_STDOUT,*)" If the Hessian has already one negative eigenvalue, "
        write(OPT_STDOUT,*)" a regular Newton-Step instead of the RFO-step will"
        write(OPT_STDOUT,*)" be taken."
     endif
  endif
  ! now we read input for frequency calculation
  read(io_optinput,nml=frequency_parameters,iostat=status)
  if (status/=0) call error_handler("read_input: namelist FREQUENCY_PARAMETERS failed")  
  if(frequency_calculation) then
     write(OPT_STDOUT,*)"      --- FREQUENCY calculation parameters ---"
     write(OPT_STDOUT,*)"Number of atoms with user defined masses",n_defined_masses
     write(OPT_STDOUT,*)"Number of atoms with inifinite mass",n_infinite_masses
     if(n_defined_masses<0) call  error_handler("read_input: n_defined_masses negative")
     if(n_infinite_masses<0) call  error_handler("read_input: n_infinite_masses negative")
     if(n_defined_masses>0) then
        allocate(defined_masses(n_defined_masses),defined_masses_index(n_defined_masses),&
             stat=alloc_stat)
        if(alloc_stat/=0) call error_handler &
           ('read_input: allocation n_defined_masses  failed')
        do i=1,n_defined_masses
           read(io_optinput,*,iostat=status) defined_masses_index(i), &
                defined_masses(i)
           if (status/=0) then
              print *,'reading masses:',i
              call error_handler('read_input: user defined masses')
           end if
        end do
     end if
     if(n_infinite_masses>0) then 
        allocate(infinite_masses_index(n_infinite_masses),&
             stat=alloc_stat)        
        if(alloc_stat/=0) call error_handler &
           ('read_input: allocation infinite_masses failed')        
        do i=1,n_infinite_masses
           read(io_optinput,*,iostat=status) infinite_masses_index(i)
           if (status/=0) then
              print *,'reading atoms with infinite masses:',i
              call error_handler('read_input: infinite masses')
           end if
        end do
     end if
  end if
  
  read(io_optinput,nml=orientation,iostat=status)
#ifndef FPP_OPTIMIZER
  if(operations_qm_epe.and.optimization) fixed_orientation=.true.
#endif
  if (status/=0) call error_handler("read_input: namelist ORIENTATION failed") 
  if(cart_format) fixed_orientation=.false.
  if (fixed_orientation) then
     if (sum((/fixed_atom_1,fixed_atom_2,fixed_atom_3/))==0 ) then
        write(OPT_STDOUT,*)" read_input: If choosing a fixed orientation do not forget"
        write(OPT_STDOUT,*)"             to specify the atoms fixed_atom_1,fixed_atom_2 "
        write(OPT_STDOUT,*)"             and fixed_atom_3"
        call error_handler(" ")
     endif
     if (fixed_atom_1==fixed_atom_2) then
        call error_handler(" read_input: FIXED_ATOM_1 = FIXED_ATOM_2 is rubbish" )
     elseif (fixed_atom_1 == fixed_atom_3) then
        call error_handler(" read_input: FIXED_ATOM_1 = FIXED_ATOM_3 is rubbish" )
     elseif (fixed_atom_2 == fixed_atom_3 ) then
        call error_handler(" read_input: FIXED_ATOM_2 = FIXED_ATOM_3 is rubbish" )
     endif
     write(OPT_STDOUT,*)" "
     write(OPT_STDOUT,*)"        ---  Fixed Orientation Mode --- "
     write(OPT_STDOUT,*)" Atom ",fixed_atom_1," will remain in its original position "
     write(OPT_STDOUT,*)" Atom ",fixed_atom_2," will be moved along the connection line "
     write(OPT_STDOUT,*)" of atoms ",fixed_atom_1," and ",fixed_atom_2
     write(OPT_STDOUT,*)" Atom ",fixed_atom_3," will be moved within the plane spanned by "
     write(OPT_STDOUT,*)" atoms ",fixed_atom_1,", ",fixed_atom_2," and ",fixed_atom_3
     write(OPT_STDOUT,*)" "
  endif

999 CONTINUE ! close units clean up and exit

  DPRINT 'opt_read_input: returnclose_iounit(io_optinput)'
  call returnclose_iounit(io_optinput)
  DPRINT 'opt_read_input: done'

end subroutine opt_read_input


  recursive subroutine task_optimizer(task)
    use strings, only: tolower
    implicit none
    character(len=*), intent(in)  :: task
    ! *** end of interface ***

    character(len=11) :: newtask
    character(len=11) :: subtask
    integer(i4_kind)  :: I


    newtask = tolower(task)
    I = INDEX( newtask,',')
    if( I /= 0 )then
       subtask = newtask(1:I-1)
       call task_optimizer(subtask)
       subtask = newtask(I+1:11)
       call task_optimizer(subtask)
       RETURN
    endif

    select case (newtask)

    case('optimizer')

#ifndef NEW_EPE
    case('relax_epe')
    optimization=.false.
    calc_cart_hess=.false.
    calc_cart_grad=.false.
    calc_hessian=.false.
    frequency_calculation=.false.
    analitic_hessian_calculated=.false.
    convert_internals=.false.
#endif

    case('numcarthess')
    calc_cart_hess=.true.
    calc_cart_grad=.false.
    calc_hessian=.false.
    optimization=.false.
    frequency_calculation=.false.
    analitic_hessian_calculated=.false.
    convert_internals=.false.

    case ('hess')
    analitic_hessian_calculated=.false.
    calc_hessian=.false.
    convert_internals=.false.

    case ('freq_num   ')
    optimization=.false.
    calc_hessian=.true.
    frequency_calculation=.true.
    analitic_hessian_calculated=.false.
    estimate_hessian=.false.
    convert_internals=.false.

    case ('freq_analyt','freq_analit')

    ! FIXME: this will be set to true later if
    ! mod(step_counter, update_hessian_iteration) == 0
    ! But why not setting it right away?
    analitic_hessian_calculated=.false.

    optimization=.false.
    calc_hessian=.false.
    frequency_calculation=.true.
    estimate_hessian=.false.
    convert_internals=.false.
    ts_search=.false.

    case('ts')
    optimization=.false.
    ts_search=.true.
    convert_internals=.false.

    case ('nohess')
    analitic_hessian_calculated=.false.

    case ('geo','opt','geoopt','GeoOpt','transit')
    if(.not.ts_search) optimization=.true.
    calc_hessian=.false.
    frequency_calculation=.false.
    convert_internals=.false.

    case ('newinternal')
    optimization=.false.
    print_internals=.true.
    convert_internals=.true.
    analitic_hessian_calculated=.false.
    calc_hessian=.false.
    frequency_calculation=.false.
    estimate_hessian=.false.
    ts_search=.false.

    case default
      print *,'optimizer: dont know task=',task
      ABORT('no such task, see tty')
    end select
  end subroutine task_optimizer

end module opt_data_module

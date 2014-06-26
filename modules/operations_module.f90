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
module operations_module
!---------------------------------------------------------------
!
!  Purpose: Database to hold operations options that say which 
!           operations are wanted for present run.
!
!  Author: TB
!  Date: 10/95
!
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: KN
! Date:   26/10/99
! Description: allows g-tensor calculations
!
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   02/2001
! Description: allows potential calculations
!
!----------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
implicit none
private         ! by default, all names are private
save
!== Interrupt end of public interface of module =================

!------------- Declaration of public input parameters -----------
logical, public :: operations_symm             , &
                   operations_scf              , &
                   operations_integral         , &
                   operations_write_input      , &
                   operations_write_input_slave, &
                   operations_get_input_out    , & !!!!!!!!!!!!AS
                   operations_post_scf         , &
                   operations_gradients        , &
                   operations_geo_opt          , &
                   operations_transit          , &
                   operations_qm_mm            , &
#ifdef WITH_MOLMECH
                   operations_qm_mm_new        , & !!!!!!!!!!!!!AS qmmm
#endif
                   operations_symadapt_gx      , &
                   operations_make_gx          , &
                   operations_read_gx          , &
                   operations_dipole           , &
                   operations_old_input        , &
                   operations_gx_highprec      , &
                   operations_gx_test          , &
                   operations_properties       , &
                   operations_response         , &
                   operations_gx_epeformat     , &
                   operations_qm_epe           , & !!!!!!!!!!!!AS
                   operations_potential        , & !!!!!!!!!!!!AS
                   operations_pseudobonds   , &
                   operations_core_density     , &
                   operations_fitbasis_opt     , &
                   operations_gtensor          , &
                   operations_hfc              , &
                   operations_epe_lattice      , &
                   operations_solvation_effect , &
#ifdef WITH_MOLMECH
                   operations_mol_mech         , & !!!!!!!!!AS_mm
#endif
                   operations_optimizer_only   , &
                   ! not supported any more but kept for compatibility
                   operations_debug_orbitals   , &
                   operations_debug_grid       , &
                   operations_quadrupel_ob     , &
                   operations_quadrupel_ff

!------------- Declaration of public input switches -------------
integer(kind=i4_kind), public :: operations_echo_input_level

!------------ public functions and subroutines ------------------
public :: operations_read, operations_write

public :: operations_bcast!()


!================================================================
! End of public interface of module
!================================================================

!------------- Definition of private default input values -------
logical, private :: df_operations_symm              = .true. , &
                    df_operations_scf               = .true. , &
                    df_operations_integral          = .true. , &
                    df_operations_write_input       = .true. , &
                    df_operations_write_full_input  = .false., &
                    df_operations_write_short_input = .false., &
                    df_operations_write_input_slave = .false., &
                    df_operations_get_input_out     = .false., & !!!!!!!!!!!!AS
                    df_operations_post_scf          = .true. , &
                    df_operations_gradients         = .true. , &
                    df_operations_geo_opt           = .false., &
                    df_operations_transit           = .false., &
                    df_operations_qm_mm             = .false., &
#ifdef WITH_MOLMECH
                    df_operations_qm_mm_new         = .false., & !!!!!!!!!!!!!!AS qmmm
#endif
                    df_operations_symadapt_gx       = .false., &
                    df_operations_make_gx           = .false., &
                    df_operations_read_gx           = .false., &
                    df_operations_dipole            = .false., &
                    df_operations_old_input         = .false., &
                    df_operations_gx_highprec       = .true., & ! AM changed in V2.2.0.11
                    df_operations_gx_test           = .false., &
                    df_operations_properties        = .false., &
                    df_operations_response          = .false., &
                    df_operations_qm_epe            = .false., & !!!!!!!!!!!AS
                    df_operations_potential         = .false., & !!!!!!!!!!!AS
!!! MF change default (compatibility?!)
!                   df_operations_gx_epeformat      = .true.,  &
                    df_operations_gx_epeformat      = .false.,  &
                    df_operations_pseudobonds    = .false., &
                    df_operations_core_density      = .false., &
                    df_operations_gtensor           = .false., &
                    df_operations_hfc               = .false., &
                    df_operations_epe_lattice       = .false., &
                    df_operations_fitbasis_opt      = .false., &
                    df_operations_solvation_effect  = .false., &
#ifdef WITH_MOLMECH
                    df_operations_mol_mech          = .false., & !!!!!!AS_mm
#endif
                    df_operations_optimizer_only    = .false., & 
                    ! not supported any more but kept for compatibility
                    df_operations_debug_orbitals    = .false., &
                    df_operations_debug_grid        = .false., &
                    df_operations_quadrupel_ob      = .false., &
                    df_operations_quadrupel_ff      = .false.

character(len=11), private :: df_task = "SinglePoint"
character(len=7),  private :: df_echo_input_level = "default"


!------------- Declaration of private input parameters -----------
logical, private :: operations_write_full_input, &
                    operations_write_short_input
character(len=11), public, protected :: operations_task
character(len=11), private  :: task
character(len=7),  private  :: echo_input_level
logical, private :: dipole, write_input_slave, gx_highprec, read_gx, properties, &
                    response, gtensor, gx_epeformat, qm_epe, solvation, hfc, &
#ifdef WITH_MOLMECH
                    qm_mm_new, & !!!!!!!!!!!!!!AS qmmm
#endif
                    get_input_out !!!!!!!!!!!!AS
!------------- Declaration of private variables -----------------
!logical, private :: namelist_tasks_used
logical, public, protected :: namelist_tasks_used

!------------- Definition of namelist ---------------------------
namelist /operations/ operations_symm             , &
                      operations_scf              , &
                      operations_integral         , &
                      operations_write_input      , &
                      operations_write_full_input , &
                      operations_write_short_input, &
                      operations_write_input_slave, &
                      operations_get_input_out    , & !!!!!!!!!!!!AS
                      operations_post_scf         , &
                      operations_gradients        , &
                      operations_geo_opt          , &
                      operations_qm_mm            , &
#ifdef WITH_MOLMECH
                      operations_qm_mm_new        , & !!!!!!!!!!!!!!!!AS qmmm
#endif
                      operations_symadapt_gx      , &
                      operations_make_gx          , &
                      operations_read_gx          , &
                      operations_qm_epe           , & !!!!!!!!!!!!!AS
                      operations_potential        , & !!!!!!!!!!!!!AS
#ifndef NEW_EPE
                      operations_epe_lattice      , &
                      operations_gx_epeformat     , &
#endif
                      operations_dipole           , &
                      operations_old_input        , &
                      operations_gx_highprec      , &
                      operations_gx_test          , &
                      operations_properties       , &
                      operations_response         , &
                      operations_core_density     , &
                      operations_fitbasis_opt     , &
                      operations_gtensor          , &
                      operations_hfc              , &
                      operations_solvation_effect , &
#ifdef WITH_MOLMECH
                      operations_mol_mech         , & !!!!!!!!!AS_mm
#endif
                      operations_optimizer_only   , & 
                      operations_debug_orbitals   , &
                      operations_debug_grid       , &
                      operations_quadrupel_ob     , &
                      operations_quadrupel_ff     , &
                      operations_fitbasis_opt      

namelist /tasks/      task                        , &
                      echo_input_level            , &
                      get_input_out               , & !!!!!!!!!!!!!AS
                      dipole                      , &
                      gtensor                     , &
                      hfc                         , &
                      solvation                   , &
                      write_input_slave           , &
                      gx_highprec                 , &
                      properties                  , &
                      read_gx                     , &
#ifdef WITH_MOLMECH
                      qm_mm_new                   , & !!!!!!!!!!!!!!AS qmmm
#endif
#ifndef NEW_EPE
                      gx_epeformat                , &
#endif
                      qm_epe                       !!!!!!!!!!!!!!AS

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine operations_read()
    ! purpose: reads data from input file
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use input_module
    use echo_input_module, only: echo_level_full, echo_level_short, &
         echo_level_default
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: status, unit
    character(len=32) :: namelist_name
    logical           :: ok
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    ! read input
    operations_qm_epe            = df_operations_qm_epe !!!!!!!!!!!!!!!!!AS
    operations_potential         = df_operations_potential !!!!!!!!!!!!!!!!AS
    operations_get_input_out     = df_operations_get_input_out !!!!!!!!!!!!AS
    operations_pseudobonds       = df_operations_pseudobonds
    operations_gx_epeformat      = df_operations_gx_epeformat 
    operations_symm              = df_operations_symm
    operations_scf               = df_operations_scf
    operations_integral          = df_operations_integral
    operations_write_input       = df_operations_write_input
    operations_write_full_input  = df_operations_write_full_input
    operations_write_short_input = df_operations_write_short_input
    operations_write_input_slave = df_operations_write_input_slave
    operations_post_scf          = df_operations_post_scf
    operations_gradients         = df_operations_gradients
    operations_geo_opt           = df_operations_geo_opt
#ifdef WITH_MOLMECH
    operations_qm_mm_new         = df_operations_qm_mm_new !!!!!!!!!!!!!!AS qmmm
#endif
    operations_qm_mm             = df_operations_qm_mm
    operations_symadapt_gx       = df_operations_symadapt_gx
    operations_make_gx           = df_operations_make_gx
    operations_read_gx           = df_operations_read_gx
    operations_dipole            = df_operations_dipole
    operations_gtensor           = df_operations_gtensor
    operations_hfc               = df_operations_hfc
    operations_epe_lattice       = df_operations_epe_lattice
    operations_response          = df_operations_response
    operations_core_density      = df_operations_core_density
    operations_fitbasis_opt      = df_operations_fitbasis_opt
    operations_old_input         = df_operations_old_input
    operations_gx_highprec       = df_operations_gx_highprec
    operations_gx_test           = df_operations_gx_test
    operations_properties        = df_operations_properties
    operations_response          = df_operations_response
    operations_solvation_effect  = df_operations_solvation_effect
#ifdef WITH_MOLMECH
    operations_mol_mech          = df_operations_mol_mech !!!!!!!AS_mm
#endif
    operations_optimizer_only    = df_operations_optimizer_only

    ! input switches which are not fully supported any more
    operations_debug_orbitals    = df_operations_debug_orbitals
    operations_debug_grid        = df_operations_debug_grid
    operations_quadrupel_ob      = df_operations_quadrupel_ob
    operations_quadrupel_ff      = df_operations_quadrupel_ff


    task                         = df_task
    echo_input_level             = df_echo_input_level
    get_input_out                = df_operations_get_input_out !!!!!!!!!!!!!AS
    dipole                       = df_operations_dipole
    write_input_slave            = df_operations_write_input_slave
    gx_highprec                  = df_operations_gx_highprec
    properties                   = df_operations_properties
    response                     = df_operations_response
    read_gx                      = df_operations_read_gx
    gtensor                      = df_operations_gtensor
    hfc                          = df_operations_hfc
    qm_epe                       = df_operations_qm_epe !!!!!!!!!!!!!!AS
#ifdef WITH_MOLMECH
    qm_mm_new                    = df_operations_qm_mm_new !!!!!!!!!!!!!!AS qmmm
#endif
    gx_epeformat                 = df_operations_gx_epeformat
    solvation                    = df_operations_solvation_effect

    ! get the namelist_name and return ok:
    ok = input_which_namelist(namelist_name)
    if( .not. ok ) call input_error("operations_read: not a valid input line")

    DPRINT "operations_read: namelist_name=",namelist_name

    SELECT CASE(namelist_name)
    CASE("tasks")

       namelist_tasks_used = .true.

       call input_read_to_intermediate()
       unit = input_intermediate_unit()
       read(unit, nml=tasks, iostat=status)
       if (status /= 0)then
          print *,'operations_read: iostat=',status,'while reading nml=tasks'
          call input_error("operations_read: namelist tasks")
       endif

       operations_task=task
       call task_operations('___reset___')
!      print*,'operations task:'//task
       call task_operations(task)
!      print*,'operations task:'//task

       select case (echo_input_level)
       case ("default")
          operations_write_input       = .true.
          operations_write_full_input  = .false.
          operations_write_short_input = .false.
          operations_echo_input_level  = echo_level_default
       case ("full   ")
          operations_write_input       = .true.
          operations_write_full_input  = .true.
          operations_write_short_input = .false.
          operations_echo_input_level  = echo_level_full
       case ("short  ")
          operations_write_input       = .true.
          operations_write_full_input  = .false.
          operations_write_short_input = .true.
          operations_echo_input_level  = echo_level_short
       case ("none  ")
          operations_write_input       = .false.
          operations_write_full_input  = .false.
          operations_write_short_input = .false.
          operations_echo_input_level  = echo_level_short
       case default
          call input_error( &
               "operations_read: invalid echo_input_level " // echo_input_level // &
               ". Allowed values are: default full short none")
       end select

       operations_dipole = dipole

       ! for response calculations dipole moments are required !
       if(operations_response) operations_dipole = .true.
       operations_gtensor= gtensor
       operations_hfc= hfc
       operations_solvation_effect = solvation
       operations_write_input_slave = write_input_slave
       operations_gx_highprec = gx_highprec
       operations_read_gx = read_gx .or. operations_symadapt_gx
       operations_gx_epeformat = gx_epeformat
       operations_qm_epe = qm_epe  !!!!!!!!!!!!!!!!AS
#ifdef WITH_MOLMECH
       operations_qm_mm_new = qm_mm_new !!!!!!!!!!!!AS qmmm
#endif
       operations_get_input_out = get_input_out !!!!!!!!!!!!!!!!AS
       if(operations_fitbasis_opt .or. &
            operations_core_density .or. &
            operations_qm_mm ) operations_qm_epe = .false. !!!!!!!!!!!!!AS

       if ( operations_qm_epe ) operations_gx_epeformat = .true.

#ifdef WITH_MOLMECH
       if(operations_qm_mm .and. operations_qm_mm_new) call input_error( &
            "operations_read: namelist task: please choose&
            &either QM_MM or QM_MM_NEW for your run") !!!!!!!!!!!!!AS
#endif

#ifdef WITH_MOLMECH
       !if(operations_qm_mm_new) operations_gradients = .true.
#endif

    CASE("operations")

       namelist_tasks_used = .false.

       call input_read_to_intermediate()
       unit = input_intermediate_unit()
       read(unit, nml=operations, iostat=status)
       if (status /= 0)then
          print *,'operations_read: iostat=',status,'while reading nml=operations'
          call input_error("operations_read: namelist operations")
       endif

#ifndef NEW_EPE
       if( operations_epe_lattice) then
          operations_symm              = .false.
          operations_scf               = .false.
          operations_integral          = .false.
          operations_post_scf          = .false.
          operations_dipole            = .false.
          operations_gradients         = .false.
          operations_geo_opt           = .false.
#ifdef WITH_MOLMECH
          operations_qm_mm_new         = .false.
#endif
          operations_qm_mm             = .false.
          operations_symadapt_gx       = .false.
          operations_make_gx           = .false.
          operations_properties        = .false.
          operations_response          = .false.
          operations_core_density      = .false.
          operations_fitbasis_opt      = .false.
          operations_qm_epe            = .false.
          operations_solvation_effect  = .false.
          operations_gtensor           = .false.
          operations_potential         = .false.
          operations_gx_test           = .false.
          operations_epe_lattice       = .true.
#ifdef WITH_MOLMECH
          operations_mol_mech          = .false. !!!!!!!!AS_mm
#endif
       endif
#endif

#ifdef WITH_MOLMECH
       if( operations_mol_mech) then     !!!!!!!!!!AS_mm
          operations_symm              = .false.
          operations_scf               = .false.
          operations_integral          = .false.
          operations_post_scf          = .false.
          operations_dipole            = .false.
          operations_gradients         = .false.
          operations_geo_opt           = .false.
          operations_qm_mm_new         = .false.
          operations_qm_mm             = .false.
          operations_symadapt_gx       = .false.
          operations_make_gx           = .false.
          operations_properties        = .false.
          operations_response          = .false.
          operations_core_density      = .false.
          operations_fitbasis_opt      = .false.
          operations_qm_epe            = .false.
          operations_solvation_effect  = .false.
          operations_gtensor           = .false.
          operations_potential         = .false.
          operations_gx_test           = .false.
          operations_epe_lattice       = .false.
       endif
#endif
       if( operations_optimizer_only) then    
#if optimizer_only_flags
          operations_scf               = .false.
          operations_integral          = .false.
          operations_post_scf          = .false.
          operations_dipole            = .false.
          operations_gradients         = .false.
          operations_geo_opt           = .false.
          operations_qm_mm_new         = .false.
          operations_qm_mm             = .false.
          operations_symadapt_gx       = .false.
          operations_make_gx           = .false.
          operations_properties        = .false.
          operations_response          = .false.
          operations_core_density      = .false.
          operations_fitbasis_opt      = .false.
          operations_qm_epe            = .false.
          operations_solvation_effect  = .false.
          operations_gtensor           = .false.
          operations_potential         = .false.
          operations_gx_test           = .false.
          operations_epe_lattice       = .false.
          operations_mol_mech          = .false. !!!!!!!!AS_mm
#endif
       endif

       if ( operations_potential) then  !!!!!!!!!!!!!AS
          operations_symm              = .true.
          operations_scf               = .true.
          operations_integral          = .true.
          operations_post_scf          = .false.
          operations_gradients         = .false.
          operations_geo_opt           = .false.
          operations_qm_mm             = .false.
          operations_symadapt_gx       = .false.
          operations_make_gx           = .false.
          operations_core_density      = .false.
          operations_fitbasis_opt      = .false.
#ifdef WITH_MOLMECH
          operations_mol_mech          = .false. !!!!!!!!AS_mm
#endif
          operations_epe_lattice       = .false. !!!!!!!!!!!!!!AS
       endif

       if ( operations_symadapt_gx ) then
            operations_read_gx    =  .true.
            ! to facilitate switching between "QM-MM" and "SymAdaptGX" mode
            operations_scf        = .false.
            operations_integral   = .false.
            operations_post_scf   = .false.
            operations_gradients  = .false.
#ifdef WITH_MOLMECH
            operations_qm_mm_new  = .false. !!!!!!!!!!!!!AS qmmm
#endif
            operations_qm_mm      = .false. 
            operations_properties = .false.
#ifdef WITH_MOLMECH
            operations_mol_mech     = .false. !!!!!!!!AS_mm
#endif
            operations_epe_lattice  = .false. !!!!!!!!!!!!!!AS
            if (operations_make_gx) call input_error &
               ("operations_read: namelist operations: &
               &Make GX together with SymAdapt GX is not possible")
            if (operations_geo_opt) call input_error &
               ("operations_read: namelist operations: &
               &SymAdapt GX in combination with geo_opt is not possible")
       endif
       ! for response calculations dipole moments are required !
       if(operations_response) operations_dipole = .true.

       ! load operations_echo_input_level according to the write_input options
       if (operations_write_full_input .and. operations_write_short_input) &
            call input_error("operations_read: namelist operations: "// &
            "write_full_input not allowed together with write_short_input" )
       operations_echo_input_level = echo_level_default
       if (operations_write_full_input) then
          operations_write_input = .true.
          operations_echo_input_level = echo_level_full
       endif
       if (operations_write_short_input) then
          operations_write_input = .true.
          operations_echo_input_level = echo_level_short
       endif

       ! check the remaining input parameters
       if (operations_quadrupel_ob) call input_error( &
            "operations_read: namelist operations: operations_quadrupel_ob no longer supported")
       if (operations_quadrupel_ff) call input_error( &
            "operations_read: namelist operations: operations_quadrupel_ff no longer supported")
       if (operations_debug_orbitals) call input_error( &
            "operations_read: namelist operations: operations_debug_orbitals no longer supported")
       if (operations_debug_grid) call input_error( &
            "operations_read: namelist operations: operations_debug_grid no longer supported")
       if ( operations_symm .and. operations_scf .and. &
            .not. operations_integral ) call input_error( &
            "operations_read: namelist operations: for running&
            & scf part together with new symmetry part, integral part is required as well")
       if ( operations_post_scf .and. .not. operations_scf ) call input_error( &
            "operations_read: namelist operations: for running&
            & post scf part, scf part is required as well")
       if ( operations_dipole .and. .not. operations_scf ) call input_error( &
            "operations_read: namelist operations: for running&
            & dipole part, scf part is required as well")
       if ( operations_gtensor .and. .not. operations_scf ) call input_error( &
            "operations_read: namelist operations: for running&
            & gtensor part, scf part is required as well")
       if ( operations_hfc .and. .not. operations_scf ) call input_error( &
            "operations_read: namelist operations: for running&
            & hyperfine coupling constants part, scf part is required as well")
       if ( operations_gradients .and. (.not.operations_post_scf) ) call input_error( &
            "operations_read: namelist operations: &
            &Sorry, You can not calculate gradients without post_scf part")
       if ( operations_geo_opt.and. (.not.operations_gradients).and. (.not.operations_gx_test)) &
            call input_error &
            ("operations_read: namelist operations: &
            &For geometry optimization set also operations_gradients=t")
       if ( operations_geo_opt .and. operations_make_gx) call input_error &
            ("operations_read: namelist operations: &
            &Either geometry optimization or make_gx - both together do not work")
       if ( operations_fitbasis_opt .and. ( operations_post_scf .or. operations_gradients ) ) &
            call input_error("operations_read: namelist operations: &
            &Neiter Post SCF or Gradient Calculations allowed for Fit Basis Optimization")
       if (operations_qm_mm) then
           if ((.not.operations_gradients).and.(.not.operations_gx_test)) &
               call input_error &
               ("operations_read: namelist operations: &
               &A QM/MM run requires operations_gradients = true")
            if (operations_make_gx) call input_error &
               ("operations_read: namelist operations: &
               &A QM/MM run together with the make_gx options is not possible")
            if (operations_geo_opt) call input_error &
               ("operations_read: namelist operations: &
               &A QM/MM run together with the geo_opt options is not possible")
       endif
       if ( operations_qm_epe ) operations_gx_epeformat = .true.

#ifdef WITH_MOLMECH
       if(operations_qm_mm .and. operations_qm_mm_new) call input_error( &
            "operations_read: namelist operations: please choose&
            &either QM_MM or QM_MM_NEW for your run") !!!!!!!!!!!!!AS

       !if(operations_qm_mm_new) operations_gradients = .true.
#endif

    CASE DEFAULT

       call input_error("operations_read: either namelist tasks or namelist operations is required")

    END SELECT

    if(operations_potential .and. operations_solvation_effect) call input_error ( &
       "operations_read: electrostatic potential and solvation_effect cannot be calculated simulteniously")

  end subroutine operations_read

  recursive subroutine task_operations(task)
    use input_module, only: input_error
    use strings, only: tolower
    implicit none
    character(len=11), intent(in)  :: task
    ! *** end of interface ***

    character(len=11) :: newtask
    character(len=11) :: subtask
    integer(i4_kind)  :: I

!!$    print *,'task_operations(>'//task//'<)'

    newtask = tolower(task)
    I = INDEX( newtask,',')
    if( I /= 0 )then
!!$       print *,'found comma at pos ',I
       subtask = newtask(1:I-1)
       call task_operations(subtask)
       subtask = newtask(I+1:11)
       call task_operations(subtask)
       RETURN
    endif

    select case (newtask)
    case ('___reset___')
       operations_core_density      = .false.
       operations_dipole            = .false.
       operations_epe_lattice       = .false.
       operations_fitbasis_opt      = .false.
       operations_geo_opt           = .false.
       operations_gradients         = .false.
       operations_gtensor           = .false.
       operations_gx_test           = .false.
       operations_integral          = .false.
       operations_make_gx           = .false.
       operations_optimizer_only    = .false.
       operations_post_scf          = .false.
       operations_potential         = .false.
       operations_properties        = .false.
       operations_qm_epe            = .false.
       operations_qm_mm             = .false.
       operations_response          = .false.
       operations_scf               = .false.
       operations_solvation_effect  = .false.
       operations_symadapt_gx       = .false.
       operations_symm              = .false.
#ifdef WITH_MOLMECH
       operations_mol_mech          = .false.
       operations_qm_mm_new         = .false.
#endif
    case ('optimizer','newinternal')
      operations_optimizer_only     = .true.
    case ('scf')
       operations_symm              = .true.
       operations_integral          = .true.
       operations_scf               = .true.
    case ('singlepoint','sing','singl','single')
       subtask = 'scf'
       call task_operations(subtask)
       operations_post_scf          = .true.
    case ('potential','pot')
       subtask = 'scf'
       call task_operations(subtask)
       operations_potential         = .true.
    case ('gradients  ','gr','gra','grad','grads')
       subtask = 'scf'
       call task_operations(subtask)
       operations_post_scf          = .true.
       operations_gradients         = .true.
    case('hess','nohess')
       subtask = 'gradients'
       call task_operations(subtask)
    case ('geoopt     ','geo','opt')
       subtask = 'gradients'
       call task_operations(subtask)
       operations_geo_opt           = .true.
    case ('freq_analyt','freq_analit','freq_num','numcarthess','cart_hess', &
          'use_eperef','make_eperef','relax_epe','fixed_epe')
       subtask = 'geo'
       call task_operations(subtask)
    case ('ts_search','ts','transit')
       operations_transit=.true.
       subtask = 'geo'
       call task_operations(subtask)
    case ('qm-mm      ','qmmm')
       subtask = 'gradients'
       call task_operations(subtask)
       operations_qm_mm             = .true.
    case ('symadaptgx ','sygx')
       operations_symm              = .true.
       operations_symadapt_gx       = .true.
       operations_read_gx           = .true.
    case ('gxtest     ','gxt')
       operations_symm              = .true.
       operations_gx_test           = .true.
       operations_read_gx           = .true.
    case ('checkinput ','chki')
       operations_symm              = .true.
    case ('makegx     ','mkgx')
       operations_symm              = .true.
       operations_make_gx           = .true.
#ifndef NEW_EPE
    case ('epe_lattice','epe','make_regref','save_regref')
       operations_epe_lattice       = .true.
#endif
    case ('properties ','prop')
       operations_symm              = .true.
       operations_integral          = .true. ! will need overlap, I guess
       operations_properties        = .true.
    case ('response   ','resp')
       subtask = 'singlepoint'
       call task_operations(subtask)
       operations_response          = .true.
    case ('coredensity')
       subtask = 'singlepoint'
       call task_operations(subtask)
       operations_core_density      = .true.
    case ('fitbasisopt')
       subtask = 'scf'
       call task_operations(subtask)
       operations_fitbasis_opt      = .true.
#ifdef WITH_MOLMECH
    case ('mol_mech   ','mm','molmech')
       operations_mol_mech          = .true.
#endif
    case ('xxx        ')
       call input_error( &
            'operations_read: You must specify a Task in Namelist Tasks.&
            &Allowed values are: SinglePoint GeoOpt Gradients CheckInput &
            &MakeGX Properties Response CoreDensity FitBasisOpt QM-MM &
            &SymAdaptGX Potential Epe_lattice Mol_mech')
    case default
       call input_error( &
            'operations_read: invalid task ' // task // &
            '. Allowed values are: SinglePoint GeoOpt Gradients CheckInput &
            &MakeGX Properties Response CoreDensity FitBasisOpt QM-MM &
            &SymAdaptGX Potential Epe_lattice Mol_mech newinternal use_eperef transit')
    end select

    if(properties)then
       call task_operations('properties ')
    endif
  end subroutine task_operations
  !*************************************************************
#if 0
  recursive subroutine task_optimizer(task)
    use input_module, only: input_error
    use strings, only: tolower
    implicit none
    character(len=11), intent(in)  :: task
    ! *** end of interface ***

    character(len=11) :: newtask
    character(len=11) :: subtask
    integer(i4_kind)  :: I

!!$    print *,'task_operations(>'//task//'<)'

    newtask = tolower(task)
    I = INDEX( newtask,',')
    if( I /= 0 )then
!!$       print *,'found comma at pos ',I
       subtask = newtask(1:I-1)
       call task_operations(subtask)
       subtask = newtask(I+1:11)
       call task_operations(subtask)
       RETURN
    endif

    select case (newtask)
    case ('___reset___')
       operations_core_density      = .false.
       operations_dipole            = .false.
       operations_epe_lattice       = .false.
       operations_fitbasis_opt      = .false.
       operations_geo_opt           = .false.
       operations_gradients         = .false.
       operations_gtensor           = .false.
       operations_gx_test           = .false.
       operations_integral          = .false.
       operations_make_gx           = .false.
       operations_optimizer_only    = .false.
       operations_post_scf          = .false.
       operations_potential         = .false.
       operations_properties        = .false.
       operations_qm_epe            = .false.
       operations_qm_mm             = .false.
       operations_response          = .false.
       operations_scf               = .false.
       operations_solvation_effect  = .false.
       operations_symadapt_gx       = .false.
       operations_symm              = .false.
#ifdef WITH_MOLMECH
       operations_mol_mech          = .false.
       operations_qm_mm_new         = .false.
#endif
    case ('optimizer')
      operations_optimizer_only     = .true.
    case ('scf')
       operations_symm              = .true.
       operations_integral          = .true.
       operations_scf               = .true.
    case ('singlepoint','sing','singl','single')
       subtask = 'scf'
       call task_operations(subtask)
       operations_post_scf          = .true.
    case ('potential','pot')
       subtask = 'scf'
       call task_operations(subtask)
       operations_potential         = .true.
    case ('gradients  ','gr','gra','grad','grads')
       subtask = 'scf'
       call task_operations(subtask)
       operations_post_scf          = .true.
       operations_gradients         = .true.
    case ('geoopt     ','geo','opt')
       subtask = 'gradients'
       call task_operations(subtask)
       operations_geo_opt           = .true.
    case ('freq_analyt','freq_analit')
       subtask = 'geo'
       call task_operations(subtask)
    case ('qm-mm      ','qmmm')
       subtask = 'gradients'
       call task_operations(subtask)
       operations_qm_mm             = .true.
    case ('symadaptgx ','sygx')
       operations_symm              = .true.
       operations_symadapt_gx       = .true.
       operations_read_gx           = .true.
    case ('gxtest     ','gxt')
       operations_symm              = .true.
       operations_gx_test           = .true.
       operations_read_gx           = .true.
    case ('checkinput ','chki')
       operations_symm              = .true.
    case ('makegx     ','mkgx')
       operations_symm              = .true.
       operations_make_gx           = .true.
#ifndef NEW_EPE
    case ('epe_lattice','epe')
       operations_epe_lattice       = .true.
#endif
    case ('properties ','prop')
       operations_symm              = .true.
       operations_properties        = .true.
    case ('response   ','resp')
       subtask = 'singlepoint'
       call task_operations(subtask)
       operations_response          = .true.
    case ('coredensity')
       subtask = 'singlepoint'
       call task_operations(subtask)
       operations_core_density      = .true.
    case ('fitbasisopt')
       subtask = 'scf'
       call task_operations(subtask)
       operations_fitbasis_opt      = .true.
#ifdef WITH_MOLMECH
    case ('mol_mech   ','mm','molmech')
       operations_mol_mech          = .true.
#endif
    case ('xxx        ')
       call input_error( &
            'operations_read: You must specify a Task in Namelist Tasks.&
            &Allowed values are: SinglePoint GeoOpt Gradients CheckInput &
            &MakeGX Properties Response CoreDensity FitBasisOpt QM-MM &
            &SymAdaptGX Potential Epe_lattice Mol_mech')
    case default
       call input_error( &
            'operations_read: invalid task ' // task // &
            '. Allowed values are: SinglePoint GeoOpt Gradients CheckInput &
            &MakeGX Properties Response CoreDensity FitBasisOpt QM-MM &
            &SymAdaptGX Potential Epe_lattice Mol_mech')
    end select

    if(properties)then
       call task_operations('properties ')
    endif
  end subroutine task_optimizer
#endif


  !*************************************************************
  subroutine operations_write(iounit,full)
    !
    ! Purpose: writes data in input format
    !
    use echo_input_module, only: start, real, intg, flag, word, stop, &
         echo_level_full, word_format
    implicit none
    integer, intent(in)           :: iounit
    logical, intent(in), optional :: full
    !** End of interface ***************************************

    logical :: original_write_input

    !   << system dependent namelist format >>
    !   !------------ Declaration of local variables
    !   integer :: status
    !   write(iounit, fmt=*, iostat=status)
    !   if (status .gt. 0) call error_handler("operations_write")
    !   write(iounit, nml=operations, iostat=status)
    !   if (status .gt. 0) call error_handler( &
    !        "operations_write: write failed at namelist operations")
    !   write(iounit, fmt=*, iostat=status)
    !   if (status .gt. 0) call error_handler("operations_write")
    !

    ! set write_input option as originally specified in the input
    original_write_input = operations_write_input .and. &
         .not.operations_write_full_input .and. &
         .not.operations_write_short_input

    if ( namelist_tasks_used ) then

       word_format = '("    ",a," = ",a13:" # ",a)' ! including quotes

       if (present(full)) then
          call start("TASKS","OPERATIONS_WRITE",iounit,echo_level_full)
       else
          call start("TASKS","OPERATIONS_WRITE",iounit,operations_echo_input_level)
       endif
       call word("TASK             ",&
            task                         ,df_task                        )
#ifdef WITH_MOLMECH
       write (iounit,'(A100)') &
            "   # Allowed values are: SinglePoint GeoOpt Gradients &
                  &CheckInput MakeGX Properties", &
            "   # Response CoreDensity FitBasisOpt &
                  &QM-MM SymAdaptGX Potential Epe_lattice Mol_mech"
#else
       write (iounit,'(A100)') &
            "   # Allowed values are: SinglePoint GeoOpt Gradients &
                  &CheckInput MakeGX Properties", &
            "   # Response CoreDensity FitBasisOpt &
                  &QM-MM SymAdaptGX Potential Epe_lattice"
#endif
       call word("ECHO_INPUT_LEVEL ",&
            echo_input_level             ,df_echo_input_level            )
       if ( present(full) .or. operations_echo_input_level == echo_level_full .or. &
            echo_input_level .ne. df_echo_input_level ) write (iounit,*) &
            "   # Allowed values are: default full short none"
       call flag("DIPOLE           ",&
            operations_dipole            ,df_operations_dipole           )
       call flag("GTENSOR          ",&
            operations_gtensor           ,df_operations_gtensor          )
       call flag("HFC          ",&
            operations_hfc               ,df_operations_hfc              )
       call flag("SOLVATION        ",&
            operations_solvation_effect  ,df_operations_solvation_effect )
       call flag("GX_HIGHPREC      ",&
            operations_gx_highprec       ,df_operations_gx_highprec      )       
       call flag("PROPERTIES      ",&
            properties                   ,df_operations_properties      )
       call flag("WRITE_INPUT_SLAVE",&
            operations_write_input_slave ,df_operations_write_input_slave)
       call flag("READ_GX          ",&
            operations_read_gx           ,df_operations_read_gx          )
       call flag("QM_EPE          ",&
            operations_qm_epe            ,df_operations_qm_epe           ) !!!!!!!!!!!!!!!AS
!!$       call flag("pseudobonds          ",&
!!$            operations_pseudobonds,df_operations_pseudobonds          )
#ifndef NEW_EPE
       call flag("GX_EPEFORMAT            ",&
            operations_gx_epeformat,df_operations_gx_epeformat         )
#endif
       call flag("GET_INPUT_OUT    ",&
            operations_get_input_out     ,df_operations_get_input_out    ) !!!!!!!!!!!!!AS
#ifdef WITH_MOLMECH
       call flag("QM_MM_NEW        ",&
            operations_qm_mm_new         ,df_operations_qm_mm_new        ) !!!!!!!!!!!!!AS qmmm
#endif
       call stop()

    else

       if (present(full)) then
          if ( full ) then
             call start("OPERATIONS","OPERATIONS_WRITE",iounit,echo_level_full)
          else
             call start("OPERATIONS","OPERATIONS_WRITE",iounit,operations_echo_input_level)
          endif
       else
          call start("OPERATIONS","OPERATIONS_WRITE",iounit,operations_echo_input_level)
       endif
#ifndef NEW_EPE
       call flag("OPERATIONS_EPE_LATTICE      ",&
            operations_epe_lattice       ,df_operations_epe_lattice      )
#endif
       call flag("OPERATIONS_SYMM             ",&
            operations_symm              ,df_operations_symm             )
       call flag("OPERATIONS_SCF              ",&
            operations_scf               ,df_operations_scf              )
       call flag("OPERATIONS_INTEGRAL         ",&
            operations_integral          ,df_operations_integral         )
       call flag("OPERATIONS_WRITE_INPUT      ",&
            original_write_input         ,df_operations_write_input      )
       call flag("OPERATIONS_WRITE_FULL_INPUT ",&
            operations_write_full_input  ,df_operations_write_full_input )
       call flag("OPERATIONS_GET_INPUT_OUT    ",&
            operations_get_input_out    ,df_operations_get_input_out     )
       call flag("OPERATIONS_WRITE_SHORT_INPUT",&
            operations_write_short_input ,df_operations_write_short_input) !!!!!!!!!!!AS
       call flag("OPERATIONS_WRITE_INPUT_SLAVE",&
            operations_write_input_slave ,df_operations_write_input_slave)
       call flag("OPERATIONS_POST_SCF         ",&
            operations_post_scf          ,df_operations_post_scf         )
       call flag("OPERATIONS_GRADIENTS        ",&
            operations_gradients         ,df_operations_gradients        )
       call flag("OPERATIONS_GEO_OPT          ",&
            operations_geo_opt           ,df_operations_geo_opt          )
#ifdef WITH_MOLMECH
       call flag("OPERATIONS_QM_MM_NEW        ",&
            operations_qm_mm_new         ,df_operations_qm_mm_new        ) !!!!!!!!!!!!!AS qmmm
#endif
       call flag("OPERATIONS_QM_MM            ",&
            operations_qm_mm             ,df_operations_qm_mm            )
       call flag("OPERATIONS_SYMADAPT_GX      ",&
            operations_symadapt_gx       ,df_operations_symadapt_gx      )
       call flag("OPERATIONS_MAKE_GX          ",&
            operations_make_gx           ,df_operations_make_gx          )
       call flag("OPERATIONS_READ_GX          ",&
            operations_read_gx           ,df_operations_read_gx          )
       call flag("OPERATIONS_QM_EPE           ",&
            operations_qm_epe            ,df_operations_qm_epe           ) !!!!!!!!!!!!AS
       call flag("OPERATIONS_POTENTIAL        ",&
            operations_potential         ,df_operations_potential        ) !!!!!!!!!!!!AS
#ifndef NEW_EPE
       call flag("OPERATIONS_GX_EPEFORMAT            ", &
            operations_gx_epeformat             ,df_operations_gx_epeformat )
#endif
       call flag("OPERATIONS_DIPOLE           ",&
            operations_dipole            ,df_operations_dipole           )
       call flag("OPERATIONS_GTENSOR          ",&
            operations_gtensor           ,df_operations_gtensor          )
       call flag("OPERATIONS_HFC          ",&
            operations_hfc               ,df_operations_hfc              )
       call flag("OPERATIONS_OLD_INPUT        ",&
            operations_old_input         ,df_operations_old_input        )
       call flag("OPERATIONS_GX_HIGHPREC      ",&
            operations_gx_highprec       ,df_operations_gx_highprec      )
       call flag("OPERATIONS_GX_TEST          ",&
            operations_gx_test           ,df_operations_gx_test          )
       call flag("OPERATIONS_PROPERTIES       ",&
            operations_properties        ,df_operations_properties      )
       call flag("OPERATIONS_RESPONSE         ",&
            operations_response          ,df_operations_response         )
       call flag("OPERATIONS_CORE_DENSITY     ",&
            operations_core_density      ,df_operations_core_density     )
       call flag("OPERATIONS_FITBASIS_OPT     ",&
            operations_fitbasis_opt      ,df_operations_fitbasis_opt     )
       call flag("OPERATIONS_SOLVATION_EFFECT ",&
            operations_solvation_effect      ,df_operations_solvation_effect )
#ifdef WITH_MOLMECH
       call flag("OPERATIONS_MOL_MECH         ",&
            operations_mol_mech          ,df_operations_mol_mech ) !!!!!!!AS_mm
#endif
       call flag("OPERATIONS_OPTIMIZER_ONLY         ",&
            operations_optimizer_only          ,df_operations_optimizer_only)
       call stop()

    endif

  end subroutine operations_write
  !*************************************************************

  subroutine operations_bcast()
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call operations_pack()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call operations_unpack()
    endif
  end subroutine operations_bcast

  !*************************************************************
  subroutine operations_pack()
    ! purpose: packs data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: commpack
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call commpack(operations_scf,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_scf")
    call commpack(operations_integral,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_integral")
    call commpack(operations_write_input,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_write_input")
    call commpack(operations_echo_input_level,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_echo_input_level")
    call commpack(operations_write_input_slave,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_write_input_slave")
    call commpack(operations_post_scf,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_post_scf")
    call commpack(operations_symm,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_post_scf")
    call commpack(operations_gradients,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_gradients")
    call commpack(operations_dipole,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_dipole")
    call commpack(operations_response,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_response")
    call commpack(operations_core_density,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_core_density")
    call commpack(operations_geo_opt,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_geo_opt")
#ifdef WITH_MOLMECH
    call commpack(operations_qm_mm_new,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_qm_mm_new")
#endif
    call commpack(operations_qm_mm,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_qm_mm")
    call commpack(operations_symadapt_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_symadapt_gx")
    call commpack(operations_fitbasis_opt,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_fitbasis_opt")
    call commpack(operations_make_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_make_gx")
    call commpack(operations_read_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_read_gx")
    call commpack(operations_old_input,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_old_input")
    call commpack(operations_gx_highprec,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_gx_highprec")
    call commpack(operations_gx_test,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_gx_test")
    ! input switches which are not fully supported any more
    call commpack(operations_debug_orbitals,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_debug_orbitals")
    call commpack(operations_debug_grid,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_debug_grid")
    call commpack(operations_quadrupel_ob,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_quadrupel_ob")
    call commpack(operations_quadrupel_ff,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_quadrupel_ff")
    call commpack(operations_properties,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_properties")
    call commpack(operations_qm_epe,info) !!!!!!!!!!!!!!!!!AS
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_qm_epe")
    call commpack(operations_potential,info) !!!!!!!!!!!!!!!!!AS
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_potential")
    call commpack(operations_pseudobonds,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_pseudobonds")
    call commpack(operations_gtensor,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_gtensor")
    call commpack(operations_hfc,info)
    if (info .ne. 0) call error_handler( &
        "operations_pack: operations_hfc")
    call commpack(operations_solvation_effect,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_solvation_effect")
#ifdef WITH_MOLMECH
    call commpack(operations_mol_mech,info) !!!!!!!!!!AS_mm
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_mol_mech")
#endif
    call commpack(operations_optimizer_only,info) 
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_optimizer_only")
  end subroutine operations_pack
  !*************************************************************


  !*************************************************************
  subroutine operations_unpack()
    ! purpose: unpacks data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: communpack
    use echo_input_module, only : echo_level_full, echo_level_short
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call communpack(operations_scf,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_scf")
    call communpack(operations_integral,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_integral")
    call communpack(operations_write_input,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_write_input")
    call communpack(operations_echo_input_level,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_echo_input_level")
    call communpack(operations_write_input_slave,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_write_input_slave")
    call communpack(operations_post_scf,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_post_scf")
    call communpack(operations_symm,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_post_scf")
    call communpack(operations_gradients,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_gradients")
    call communpack(operations_dipole,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_dipole")
    call communpack(operations_response,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_response")
    call communpack(operations_core_density,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_core_density")
    call communpack(operations_geo_opt,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_geo_opt")
#ifdef WITH_MOLMECH
    call communpack(operations_qm_mm_new,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_qm_mm_new")
#endif
    call communpack(operations_qm_mm,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_qm_mm")
    call communpack(operations_symadapt_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_symadapt_gx")
    call communpack(operations_fitbasis_opt,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_fitbasis_opt")
    call communpack(operations_make_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_make_gx")
    call communpack(operations_read_gx,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_read_gx")
    call communpack(operations_old_input,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_old_input")
    call communpack(operations_gx_highprec,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_gx_highprec")
    call communpack(operations_gx_test,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_gx_test")
    ! input switches which are not fully supported any more
    call communpack(operations_debug_orbitals,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_debug_orbitals")
    call communpack(operations_debug_grid,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_debug_grid")
    call communpack(operations_quadrupel_ob,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_quadrupel_ob")
    call communpack(operations_quadrupel_ff,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_quadrupel_ff")
    call communpack(operations_properties,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_properties")
    call communpack(operations_qm_epe,info) !!!!!!!!!!!!!!!AS
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_qm_epe")
    call communpack(operations_potential,info) !!!!!!!!!!!!!!!AS
    if (info .ne. 0) call error_handler( &
         "operations_unpack: potential")
    call communpack(operations_pseudobonds,info)
    if (info .ne. 0) call error_handler( &
         "operations_pack: operations_pseudobonds")
    call communpack(operations_gtensor,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_gtensor")
    call communpack(operations_hfc,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_hfc")
    call communpack(operations_solvation_effect,info)
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_solvation_effect")
#ifdef WITH_MOLMECH
    call communpack(operations_mol_mech,info) !!!!!!!!!!!!AS_mm
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_mol_mech")
#endif
    call communpack(operations_optimizer_only,info) 
    if (info .ne. 0) call error_handler( &
         "operations_unpack: operations_optimizer_only")


    ! restore unsent input parameter
    operations_write_full_input =  &
         operations_echo_input_level == echo_level_full
    operations_write_short_input = &
         operations_echo_input_level == echo_level_short
    
  end subroutine operations_unpack
  !*************************************************************


!--------------- End of module ----------------------------------
end module operations_module

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
module  output_module
!---------------------------------------------------------------
!
!  Purpose: Database to hold output options that say how much
!           output is wanted.
!           The amount of output can be controlled in two ways:
!            1. Setting the output level. When the input is read,
!               the output level is evaluated first setting the
!               logical variables for individual output options
!               to default settings depending on the level
!            2. Setting individual options by logical variables
!            All cards are optional and can be omitted. The output
!            level (namelist output) must be given before the other
!            cards but apart from theat theit order is arbitrary.
!
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
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
#include "def.h"
use type_module ! type specification parameters
use iounitadmin_module ! iounitadmin_use_xxx are read in here
implicit none
save ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ Declaration of public input variables -------------

integer(kind=i4_kind ), public :: output_level ! default == 1

! namelist /output_trace/
logical, public :: output_slaveoperations = .false., &
                   output_main_master     = .true. , &
                   output_main_symm       = .false., &
                   output_main_scf        = .false., &
                   output_properties      = .true.,  &
                   output_read_input      = .false., &
                   output_write_input     = .false., &
                   output_pvm_msgtags     = .false., &
                   output_main_dipole     = .false., &
                   output_efield_calc     = .false.

! namelist /output_config/
logical, public :: output_pvmconfig  = .true. , &
                   output_operations = .true.

! namelist /output_unique_atoms/
logical, public :: output_atoms              = .true. , &
                   output_symadapt           = .false., &
                   output_renorm             = .false., &
                   output_orbitalprojections = .false.

! namelist /output_scf/
logical, public :: output_scfloops       = .true. , &
                   output_reoccup        = .false., &
                   output_fermi          = .false., &
                   output_fermi_newton   = .false., &
                   output_chargefit      = .false., &
                   output_condition      = .false., &
                   output_hamiltonian    = .false., &
                   output_eigendata      = .false., &
                   output_each_eigendata = .false., &
                   output_eigen_strategy = .false., &
                   output_densmat        = .false., &
                   output_grid           = .true. , &
                   output_overlap        = .false., &
                   output_data_saved     = .false., &
                   output_data_read      = .false., &
                   output_spectrum       = .false.

! namelist /output_scf/, unused entries. Will be accepted in the input
! but ignored:
logical, private :: output_occlevels = .false.

integer(kind=i4_kind), public :: output_n_density_dev, &
                                 output_n_coeff_dev

! namelist /output_sym/
logical, public :: output_symmetry        = .true. , &
                   output_nonsymequivvecs = .false.

! namelist /output_integral/
logical, public :: output_int_fitcontract      = .false., &
                   output_int_2c_fit           = .false., &
                   output_int_solhrules        = .false., &
                   output_int_parameters       = .true. , &
                   output_int_progress         = .false., &
                   output_int_detailedprogress = .false., &
                   output_int_taskdistribution = .false., &
                   output_int_quadrupelstore   = .false., &
                   output_int_loops            = .false., &
                   output_int_deeploops        = .false., &
                   output_int_data             = .false.

! namelist /output_timing/
logical, public :: output_timing_summary           = .true. , &
                   output_timing_detailedsummary   = .false., &
                   output_timing_integrals         = .true. , &
                   output_timing_detailedintegrals = .false., &
                   output_timing_scfloops          = .false., &
                   output_timing_scf               = .true. , &
                   output_timing_detailedscf       = .false., &
                   output_timing_post_scf          = .true. , &
                   output_timing_detailedpostscf   = .false., &
                   output_timing_slaves            = .false., &
                   output_timing_interrupts        = .false.

! following parameter renamed to ..scf, old values kept for consistency
! with old input files, will give an warning message
logical, private :: output_timing_posthoc           = .true. , &
                   output_timing_detailedposthoc   = .false., &
                   output_post_hoc_main = .false.

! namelist /output_post_scf/
logical, public :: output_post_scf_main = .false., &
                   output_main_gradient = .false.

! namelist /output_dipole/
logical, public :: output_dipole_detailed       = .false., &
                   output_dipole_transitionm_f  = .false., & ! formatted
                   output_dipole_transitionm_uf = .false., & ! unformatted
                   output_dipole_simol          = .false., &
                   output_dipole_optimizer      = .false.

! namelist /output_response/
logical, public :: output_response_general  = .false., &
                   output_response_detailed = .false., &
                   output_response_debug    = .false.

! namelist /output_solvation/
logical, public :: output_cavity_data, &
                   output_cavity_long,&
                   output_solv_grads

!------------ public functions and subroutines ------------------
public :: output_read, output_write

public :: output_bcast!()


!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of private input variables ------------

!------------ Definition of default variables -------------------

! namelist /output/
integer(kind=i4_kind) :: df_output_level = 1

! namelist /output_trace/
logical :: df_output_slaveoperations, &
           df_output_main_master    , &
           df_output_main_symm      , &
           df_output_main_scf       , &
           df_output_properties     , &
           df_output_read_input     , &
           df_output_write_input    , &
           df_output_pvm_msgtags    , &
           df_output_main_dipole    , &
           df_output_efield_calc

! namelist /output_config/
logical :: df_output_pvmconfig , &
           df_output_operations

! namelist /output_unique_atoms/
logical :: df_output_atoms             , &
           df_output_symadapt          , &
           df_output_renorm            , &
           df_output_orbitalprojections

! namelist /output_scf/
logical :: df_output_scfloops      , &
           df_output_reoccup       , &
           df_output_fermi         , &
           df_output_fermi_newton  , &
           df_output_chargefit     , &
           df_output_condition     , &
           df_output_hamiltonian   , &
           df_output_eigendata     , &
           df_output_each_eigendata, &
           df_output_eigen_strategy, &
           df_output_densmat       , &
           df_output_grid          , &
           df_output_overlap       , &
           df_output_data_saved    , &
           df_output_data_read     , &
           df_output_spectrum
integer(kind=i4_kind) :: df_output_n_density_dev, &
                         df_output_n_coeff_dev

! namelist /output_sym/
logical :: df_output_symmetry       , &
           df_output_nonsymequivvecs

! namelist /output_integral/
logical :: df_output_int_fitcontract     , &
           df_output_int_2c_fit          , &
           df_output_int_solhrules       , &
           df_output_int_parameters      , &
           df_output_int_progress        , &
           df_output_int_detailedprogress, &
           df_output_int_taskdistribution, &
           df_output_int_quadrupelstore  , &
           df_output_int_loops           , &
           df_output_int_deeploops       , &
           df_output_int_data

! namelist /output_timing/
logical :: df_output_timing_summary    , &
           df_output_timing_detailedsum, &
           df_output_timing_integrals  , &
           df_output_timing_detailedint, &
           df_output_timing_scfloops   , &
           df_output_timing_scf        , &
           df_output_timing_detailedscf, &
           df_output_timing_post_scf   , &
           df_output_timing_detailedps , &
           df_output_timing_slaves     , &
           df_output_timing_interrupts

! namelist /output_post_scf/
logical :: df_output_post_scf_main, &
           df_output_main_gradient

! namelist /output_dipole/
logical :: df_output_dipole_detailed      , &
           df_output_dipole_transitionm_f , &
           df_output_dipole_transitionm_uf, &
           df_output_dipole_simol, &
           df_output_dipole_optimizer

! namelist /output_response/
logical, public :: df_output_response_general , &
                   df_output_response_detailed , &
                   df_output_response_debug

! namelist /output_solvation/
logical, public :: df_output_cavity_data = .false. ,&
                   df_output_cavity_long = .false., &
                   df_output_solv_grads  = .false.

!------------ Definition of input namelist --------------

namelist /output/                output_level

namelist /output_trace/          output_slaveoperations, &
                                 output_main_master    , &
                                 output_main_symm      , &
                                 output_main_scf       , &
                                 output_properties     , &
                                 output_read_input     , &
                                 output_write_input    , &
                                 output_pvm_msgtags    , &
                                 output_main_dipole    , &
                                 output_efield_calc


namelist /output_config/         output_pvmconfig , &
                                 output_operations

namelist /output_unique_atoms/   output_atoms             , &
                                 output_symadapt          , &
                                 output_renorm            , &
                                 output_orbitalprojections

namelist /output_scf/            output_scfloops      , &
                                 output_reoccup       , &
                                 output_fermi         , &
                                 output_fermi_newton  , &
                                 output_chargefit     , &
                                 output_condition     , &
                                 output_hamiltonian   , &
                                 output_eigendata     , &
                                 output_each_eigendata, &
                                 output_eigen_strategy, &
                                 output_densmat       , &
                                 output_grid          , &
                                 output_occlevels     , & ! deprecated, ignored
                                 output_overlap       , &
                                 output_data_saved    , &
                                 output_data_read     , &
                                 output_n_density_dev , &
                                 output_n_coeff_dev,    &
                                 output_spectrum

namelist /output_sym/            output_symmetry       , &
                                 output_nonsymequivvecs

namelist /output_integral/       output_int_fitcontract     , &
                                 output_int_2c_fit          , &
                                 output_int_solhrules       , &
                                 output_int_parameters      , &
                                 output_int_progress        , &
                                 output_int_detailedprogress, &
                                 output_int_taskdistribution, &
                                 output_int_quadrupelstore  , &
                                 output_int_loops           , &
                                 output_int_deeploops       , &
                                 output_int_data

namelist /output_timing/         output_timing_summary          , &
                                 output_timing_detailedsummary  , &
                                 output_timing_integrals        , &
                                 output_timing_detailedintegrals, &
                                 output_timing_scfloops         , &
                                 output_timing_scf              , &
                                 output_timing_detailedscf      , &
                                 output_timing_posthoc          , &
                                 output_timing_post_scf         , &
                                 output_timing_detailedposthoc  , &
                                 output_timing_detailedpostscf  , &
                                 output_timing_slaves           , &
                                 output_timing_interrupts

namelist /output_post_scf/       output_post_hoc_main, &
                                 output_post_scf_main, &
                                 output_main_gradient

! doubles list namelist output_post_scf, is here for
! consistency with old input
namelist /output_post_hoc/       output_post_hoc_main, &
                                 output_post_scf_main, &
                                 output_main_gradient

namelist /output_dipole/         output_dipole_detailed      , &
                                 output_dipole_transitionm_f , &
                                 output_dipole_transitionm_uf, &
                                 output_dipole_simol,          &
                                 output_dipole_optimizer

namelist /output_response/       output_response_general, &
                                 output_response_detailed, &
                                 output_response_debug

namelist /output_solvation/      output_cavity_data, &
                                 output_cavity_long, &
                                 output_solv_grads

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine output_set_defaults
    ! purpose: set the defaults according to the output level
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use input_module, only: input_error
    !------------ Executable code --------------------------------

    df_output_n_density_dev = 5
    df_output_n_coeff_dev   = 5
    select case(output_level)
    case(0)
       ! namelist /output_trace/
       df_output_slaveoperations = .false.
       df_output_main_master     = .false.
       df_output_main_symm       = .false.
       df_output_main_scf        = .false.
       df_output_properties      =  .true.
       df_output_read_input      = .false.
       df_output_write_input     = .false.
       df_output_pvm_msgtags     = .false.
       df_output_main_dipole     = .false.
       df_output_efield_calc     = .false.

       ! namelist /output_config/
       df_output_pvmconfig  = .false.
       df_output_operations = .false.

       ! namelist /output_unique_atoms/
       df_output_atoms              = .false.
       df_output_symadapt           = .false.
       df_output_renorm             = .false.
       df_output_orbitalprojections = .false.

       ! namelist /output_scf/
       df_output_scfloops       = .false.
       df_output_reoccup        = .false.
       df_output_fermi          = .false.
       df_output_fermi_newton   = .false.
       df_output_chargefit      = .false.
       df_output_condition      = .false.
       df_output_hamiltonian    = .false.
       df_output_eigendata      = .false.
       df_output_each_eigendata = .false.
       df_output_eigen_strategy = .false.
       df_output_densmat        = .false.
       df_output_grid           = .false.
       df_output_overlap        = .false.
       df_output_data_saved     = .false.
       df_output_data_read      = .false.
       df_output_spectrum       = .false.

       ! namelist /output_sym/
       df_output_symmetry        = .false.
       df_output_nonsymequivvecs = .false.

       ! namelist /output_integral/
       df_output_int_fitcontract      = .false.
       df_output_int_2c_fit           = .false.
       df_output_int_solhrules        = .false.
       df_output_int_parameters       = .false.
       df_output_int_progress         = .false.
       df_output_int_detailedprogress = .false.
       df_output_int_taskdistribution = .false.
       df_output_int_quadrupelstore   = .false.
       df_output_int_loops            = .false.
       df_output_int_deeploops        = .false.
       df_output_int_data             = .false.

       ! namelist /output_timing/
       df_output_timing_summary     = .false.
       df_output_timing_detailedsum = .false.
       df_output_timing_integrals   = .false.
       df_output_timing_detailedint = .false.
       df_output_timing_scfloops    = .false.
       df_output_timing_scf         = .false.
       df_output_timing_detailedscf = .false.
       df_output_timing_post_scf     = .false.
       df_output_timing_detailedps  = .false.
       df_output_timing_slaves      = .false.
       df_output_timing_interrupts  = .false.

       ! namelist /output_post_scf/
       df_output_post_scf_main = .false.
       df_output_main_gradient = .false.

       ! namelist /output_dipole/
       df_output_dipole_detailed       = .false.
       df_output_dipole_transitionm_f  = .false.
       df_output_dipole_transitionm_uf = .false.
       df_output_dipole_simol          = .false.
       df_output_dipole_optimizer      = .false.

       ! namelist /output_response/
       df_output_response_general  = .false.
       df_output_response_detailed = .false.
       df_output_response_debug    = .false.

    case(1)
       ! namelist /output_trace/
       df_output_slaveoperations = .false.
       df_output_main_master     = .false.
       df_output_main_symm       = .false.
       df_output_main_scf        = .false.
       df_output_properties      = .true.
       df_output_read_input      = .false.
       df_output_write_input     = .false.
       df_output_pvm_msgtags     = .false.
       df_output_main_dipole     = .false.
       df_output_efield_calc     = .false.

       ! namelist /output_config/
       df_output_pvmconfig  = .true.
       df_output_operations = .true.

       ! namelist /output_unique_atoms/
       df_output_atoms              = .true.
       df_output_symadapt           = .false.
       df_output_renorm             = .false.
       df_output_orbitalprojections = .false.

       ! namelist /output_scf/
       df_output_scfloops       = .true.
       df_output_reoccup        = .false.
       df_output_fermi          = .true.
       df_output_fermi_newton   = .false.
       df_output_chargefit      = .false.
       df_output_condition      = .true.
       df_output_hamiltonian    = .false.
       df_output_eigendata      = .false.
       df_output_each_eigendata = .false.
       df_output_eigen_strategy = .false.
       df_output_densmat        = .false.
       df_output_grid           = .true.
       df_output_overlap        = .false.
       df_output_data_saved     = .false.
       df_output_data_read      = .false.
       df_output_spectrum       = .false.

       ! namelist /output_sym/
       df_output_symmetry        = .true.
       df_output_nonsymequivvecs = .false.

       ! namelist /output_integral/
       df_output_int_fitcontract      = .false.
       df_output_int_2c_fit           = .false.
       df_output_int_solhrules        = .false.
       df_output_int_parameters       = .true.
       df_output_int_progress         = .false.
       df_output_int_detailedprogress = .false.
       df_output_int_taskdistribution = .false.
       df_output_int_quadrupelstore   = .false.
       df_output_int_loops            = .false.
       df_output_int_deeploops        = .false.
       df_output_int_data             = .false.

       ! namelist /output_timing/
       df_output_timing_summary     = .true.
       df_output_timing_detailedsum = .false.
       df_output_timing_integrals   = .true.
       df_output_timing_detailedint = .false.
       df_output_timing_scfloops    = .false.
       df_output_timing_scf         = .true.
       df_output_timing_detailedscf = .false.
       df_output_timing_post_scf     = .true.
       df_output_timing_detailedps  = .false.
       df_output_timing_slaves      = .false.
       df_output_timing_interrupts  = .false.

       ! namelist /output_post_scf/
       df_output_post_scf_main = .false.
       df_output_main_gradient = .false.

       ! namelist /output_dipole/
       df_output_dipole_detailed       = .true.
       df_output_dipole_transitionm_f  = .false.
       df_output_dipole_transitionm_uf = .false.
       df_output_dipole_simol          = .false.
       df_output_dipole_optimizer      = .false.

       ! namelist /output_response/
       df_output_response_general  = .true.
       df_output_response_detailed = .false.
       df_output_response_debug    = .false.

    case(5)
       ! namelist /output_trace/
       df_output_slaveoperations = .true.
       df_output_main_master     = .true.
       df_output_main_symm       = .true.
       df_output_main_scf        = .true.
       df_output_properties      = .true.
       df_output_read_input      = .true.
       df_output_write_input     = .true.
       df_output_pvm_msgtags     = .false.
       df_output_main_dipole     = .true.
       df_output_efield_calc     = .true.

       ! namelist /output_config/
       df_output_pvmconfig  = .true.
       df_output_operations = .true.

       ! namelist /output_unique_atoms/
       df_output_atoms              = .true.
       df_output_symadapt           = .false.
       df_output_renorm             = .false.
       df_output_orbitalprojections = .true.

       ! namelist /output_scf/
       df_output_scfloops       = .true.
       df_output_reoccup        = .true.
       df_output_fermi          = .true.
       df_output_fermi_newton   = .false.
       df_output_chargefit      = .false.
       df_output_condition      = .true.
       df_output_hamiltonian    = .false.
       df_output_eigendata      = .false.
       df_output_each_eigendata = .false.
       df_output_eigen_strategy = .true.
       df_output_densmat        = .false.
       df_output_grid           = .true.
       df_output_overlap        = .false.
       df_output_data_saved     = .false.
       df_output_data_read      = .false.
       df_output_spectrum       = .true.

       ! namelist /output_sym/
       df_output_symmetry        = .true.
       df_output_nonsymequivvecs = .true.

       ! namelist /output_integral/
       df_output_int_fitcontract      = .false.
       df_output_int_2c_fit           = .true.
       df_output_int_solhrules        = .false.
       df_output_int_parameters       = .true.
       df_output_int_progress         = .true.
       df_output_int_detailedprogress = .true.
       df_output_int_taskdistribution = .true.
       df_output_int_quadrupelstore   = .false.
       df_output_int_loops            = .true.
       df_output_int_deeploops        = .false.
       df_output_int_data             = .false.

       ! namelist /output_timing/
       df_output_timing_summary  = .true.
       df_output_timing_detailedsum = .true.
       df_output_timing_integrals   = .true.
       df_output_timing_detailedint = .true.
       df_output_timing_scfloops    = .true.
       df_output_timing_scf         = .true.
       df_output_timing_detailedscf = .true.
       df_output_timing_post_scf     = .true.
       df_output_timing_detailedps  = .true.
       df_output_timing_slaves      = .true.
       df_output_timing_interrupts  = .false.

       ! namelist /output_post_scf/
       df_output_post_scf_main = .true.
       df_output_main_gradient = .true.

       ! namelist /output_dipole/
       df_output_dipole_detailed       = .true.
       df_output_dipole_transitionm_f  = .true.
       df_output_dipole_transitionm_uf = .false.
       df_output_dipole_simol          = .false.
       df_output_dipole_optimizer      = .false.

       ! namelist /output_response/
       df_output_response_general  = .true.
       df_output_response_detailed = .true.
       df_output_response_debug    = .false.

    case(10)
       ! namelist /output_trace/
       df_output_slaveoperations = .true.
       df_output_main_master     = .true.
       df_output_main_symm       = .true.
       df_output_main_scf        = .true.
       df_output_properties      = .true.
       df_output_read_input      = .true.
       df_output_write_input     = .true.
       df_output_pvm_msgtags     = .true.
       df_output_main_dipole     = .true.
       df_output_efield_calc     = .true.

       ! namelist /output_config/
       df_output_pvmconfig  = .true.
       df_output_operations = .true.

       ! namelist /output_unique_atoms/
       df_output_atoms              = .true.
       df_output_symadapt           = .true.
       df_output_renorm             = .true.
       df_output_orbitalprojections = .true.

       ! namelist /output_scf/
       df_output_scfloops       = .true.
       df_output_reoccup        = .true.
       df_output_fermi          = .true.
       df_output_fermi_newton   = .true.
       df_output_chargefit      = .true.
       df_output_condition      = .true.
       df_output_hamiltonian    = .true.
       df_output_eigendata      = .true.
       df_output_each_eigendata = .true.
       df_output_eigen_strategy = .true.
       df_output_densmat        = .true.
       df_output_grid           = .true.
       df_output_overlap        = .true.
       df_output_data_saved     = .true.
       df_output_data_read      = .true.
       df_output_spectrum       = .true.

       ! namelist /output_sym/
       df_output_symmetry        = .true.
       df_output_nonsymequivvecs = .true.

       ! namelist /output_integral/
       df_output_int_fitcontract      = .true.
       df_output_int_2c_fit           = .true.
       df_output_int_solhrules        = .true.
       df_output_int_parameters       = .true.
       df_output_int_progress         = .true.
       df_output_int_detailedprogress = .true.
       df_output_int_taskdistribution = .true.
       df_output_int_quadrupelstore   = .true.
       df_output_int_loops            = .true.
       df_output_int_deeploops        = .true.
       df_output_int_data             = .true.

       ! namelist /output_timing/
       df_output_timing_summary     = .true.
       df_output_timing_detailedsum = .true.
       df_output_timing_integrals   = .true.
       df_output_timing_detailedint = .true.
       df_output_timing_scfloops    = .true.
       df_output_timing_scf         = .true.
       df_output_timing_detailedscf = .true.
       df_output_timing_post_scf     = .true.
       df_output_timing_detailedps  = .true.
       df_output_timing_slaves      = .true.
       df_output_timing_interrupts  = .true.

       ! namelist /output_post_scf/
       df_output_post_scf_main = .true.
       df_output_main_gradient = .true.

       ! namelist /output_dipole/
       df_output_dipole_detailed       = .true.
       df_output_dipole_transitionm_f  = .true.
       df_output_dipole_transitionm_uf = .false.
       df_output_dipole_simol          = .false.
       df_output_dipole_optimizer      = .false.

       ! namelist /output_response/
       df_output_response_general  = .true.
       df_output_response_detailed = .true.
       df_output_response_debug    = .true.

    case default
       call input_error("output_read: illegal output_level")
    end select

  end subroutine output_set_defaults
  !*************************************************************


  !*************************************************************
  subroutine output_read
    ! purpose: reads data from input file
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use input_module
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: status, unit
    character(len=32) :: namelist_name
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    unit = input_intermediate_unit()
    ! read level
    output_level = df_output_level
    if ( input_line_is_namelist("output") ) then
       call input_read_to_intermediate
       read(unit, nml=output, iostat=status)
       if (status .gt. 0) call input_error( &
            "output_read: namelist output")
    endif

    ! set defaults according to the output_level specified
    call output_set_defaults

    ! now, initialize the input variables
    ! namelist /output_trace/
    output_slaveoperations = df_output_slaveoperations
    output_main_master     = df_output_main_master
    output_main_symm       = df_output_main_symm
    output_main_scf        = df_output_main_scf
    output_properties      = df_output_properties
    output_read_input      = df_output_read_input
    output_write_input     = df_output_write_input
    output_pvm_msgtags     = df_output_pvm_msgtags
    output_main_dipole     = df_output_main_dipole
    output_efield_calc     = df_output_efield_calc

    ! namelist /output_config/
    output_pvmconfig  = df_output_pvmconfig
    output_operations = df_output_operations

    ! namelist /output_unique_atoms/
    output_atoms              = df_output_atoms
    output_symadapt           = df_output_symadapt
    output_renorm             = df_output_renorm
    output_orbitalprojections = df_output_orbitalprojections

    ! namelist /output_scf/
    output_scfloops       = df_output_scfloops
    output_reoccup        = df_output_reoccup
    output_fermi          = df_output_fermi
    output_fermi_newton   = df_output_fermi_newton
    output_chargefit      = df_output_chargefit
    output_condition      = df_output_condition
    output_hamiltonian    = df_output_hamiltonian
    output_eigendata      = df_output_eigendata
    output_each_eigendata = df_output_each_eigendata
    output_eigen_strategy = df_output_eigen_strategy
    output_densmat        = df_output_densmat
    output_grid           = df_output_grid
    output_overlap        = df_output_overlap
    output_data_saved     = df_output_data_saved
    output_data_read      = df_output_data_read
    output_n_density_dev  = df_output_n_density_dev
    output_n_coeff_dev    = df_output_n_coeff_dev

    ! namelist /output_sym/
    output_symmetry        = df_output_symmetry
    output_nonsymequivvecs = df_output_nonsymequivvecs

    ! namelist /output_integral/
    output_int_fitcontract      = df_output_int_fitcontract
    output_int_2c_fit           = df_output_int_2c_fit
    output_int_solhrules        = df_output_int_solhrules
    output_int_parameters       = df_output_int_parameters
    output_int_progress         = df_output_int_progress
    output_int_detailedprogress = df_output_int_detailedprogress
    output_int_taskdistribution = df_output_int_taskdistribution
    output_int_quadrupelstore   = df_output_int_quadrupelstore
    output_int_loops            = df_output_int_loops
    output_int_deeploops        = df_output_int_deeploops
    output_int_data             = df_output_int_data

    ! namelist /output_timing/
    output_timing_summary           = df_output_timing_summary
    output_timing_detailedsummary   = df_output_timing_detailedsum
    output_timing_integrals         = df_output_timing_integrals
    output_timing_detailedintegrals = df_output_timing_detailedint
    output_timing_scfloops          = df_output_timing_scfloops
    output_timing_scf               = df_output_timing_scf
    output_timing_detailedscf       = df_output_timing_detailedscf
    output_timing_post_scf          = df_output_timing_post_scf
    output_timing_detailedpostscf   = df_output_timing_detailedps
    output_timing_slaves            = df_output_timing_slaves
    output_timing_interrupts        = df_output_timing_interrupts
    ! for consistency (will be redirected later on)
    output_timing_posthoc           = df_output_timing_post_scf
    output_timing_detailedposthoc   = df_output_timing_detailedps

    ! namelist /output_post_scf/
    output_post_hoc_main = df_output_post_scf_main ! input consistency
    output_post_scf_main = df_output_post_scf_main
    output_main_gradient = df_output_main_gradient

    ! namelist /output_dipole/
    output_dipole_detailed       = df_output_dipole_detailed
    output_dipole_transitionm_f  = df_output_dipole_transitionm_f
    output_dipole_transitionm_uf = df_output_dipole_transitionm_uf

    ! namelist /output_response/
    output_response_general  = df_output_response_general
    output_response_detailed = df_output_response_detailed
    output_response_debug    = df_output_response_debug
    output_dipole_simol      = df_output_dipole_simol
    output_dipole_optimizer  = df_output_dipole_optimizer

    ! namelist /output_solvation/
    output_cavity_data = df_output_cavity_data
    output_cavity_long = df_output_cavity_long
    output_solv_grads  = df_output_solv_grads

    ! set individual options
    do while ( input_which_namelist(namelist_name) )
       !
       ! Exit when namelist is not supposed to be handled here.
       ! By now all namelists that start with output_ are handled
       ! here:
       !
       select case ( namelist_name )

       case ("output_trace")
          call input_read_to_intermediate
          read(unit, nml=output_trace, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_trace")
       case ("output_config")
          call input_read_to_intermediate
          read(unit, nml=output_config, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_config")
       case ("output_unique_atoms")
          call input_read_to_intermediate
          read(unit, nml=output_unique_atoms, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_unique_atoms")
       case ("output_scf")
          call input_read_to_intermediate
          read(unit, nml=output_scf, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_scf")
       case ("output_sym")
          call input_read_to_intermediate
          read(unit, nml=output_sym, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_sym")
       case ("output_integral")
          call input_read_to_intermediate
          read(unit, nml=output_integral, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_integral")
       case ("output_timing")
          call input_read_to_intermediate
          read(unit, nml=output_timing, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_timing")
       case ("output_post_hoc")
          ! this case should be removed lateron, stays here
          ! only for consistency
          WARN('output_post_hoc has been renamed')
          print *, " name of this list should be output_post_scf"
          print *, " list will be redirected"
          print *, " Please conside updating your input file"
          WARN('Make sure not to use old and new one')
          call input_read_to_intermediate
          read(unit, nml=output_post_hoc, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_post_hoc")
       case ("output_post_scf")
          call input_read_to_intermediate
          read(unit, nml=output_post_scf, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_post_scf")
       case ("output_dipole")
          call input_read_to_intermediate
          read(unit, nml=output_dipole, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_dipole")
       case ("output_response")
          call input_read_to_intermediate
          read(unit, nml=output_response, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_response")
       case ("output_solvation")
          call input_read_to_intermediate
          read(unit, nml=output_solvation, iostat=status)
          if (status .gt. 0) call input_error( &
               "output_read: namelist output_solvation")
       case default
           exit ! do while loop
       end select

       if( output_cavity_long) output_cavity_data = .true.
    enddo


    if (output_post_hoc_main .neqv. df_output_post_scf_main) then
      WARN('output_post_hoc_main has been renamed')
      print *, " output_post_hoc_main is now output_post_scf_main"
      print *, " Now setting output_post_scf_main to output_post_hoc_main value"
      print *, " Please conside updating your input file"
      output_post_scf_main = output_post_hoc_main
    endif

    if (output_timing_posthoc .neqv. df_output_timing_post_scf) then
      WARN('output_timing_posthoc has been renamed')
      print *, " output_timing_posthoc is now output_timing_post_scf"
      print *, " Now setting output_timing_post_scf to output_timing_posthoc value"
      print *, " Please conside updating your input file"
      output_timing_post_scf = output_timing_posthoc
    endif

    if (output_timing_detailedposthoc .neqv. df_output_timing_detailedps) then
      WARN('output_timing_detailedposthoc has been renamed')
      print *, " output_timing_detailedposthoc is now output_timing_detailedpostscf"
      print *, " Now setting output_timing_detailedpostscf to output_timing_detailedposthoc value"
      print *, " Please conside updating your input file"
      output_timing_detailedpostscf = output_timing_detailedposthoc
    endif
  end subroutine output_read
  !*************************************************************


  !*************************************************************
  subroutine output_write(iounit)
    !
    ! Purpose: reads data from input file
    !
    use echo_input_module, only: start, real, intg, flag, stop, &
         echo_level_full
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    !** End of interface ***************************************

    !
    ! FIXME: use default format for integers/reals/flags ...
    !

    call start("OUTPUT","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call intg("OUTPUT_LEVEL",output_level,df_output_level)
    if ( operations_echo_input_level == echo_level_full .or. &
         output_level .ne. df_output_level ) &
         write(iounit,'(a34)') "    # legal values are 0, 1, 5, 10" !!!!!!!!!!!!!!!!!!
    call stop()

    call start("OUTPUT_TRACE","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_SLAVEOPERATIONS", &
               output_slaveoperations , df_output_slaveoperations)
    call flag("OUTPUT_MAIN_MASTER    ", &
               output_main_master     , df_output_main_master    )
    call flag("OUTPUT_MAIN_SYMM      ", &
               output_main_symm       , df_output_main_symm      )
    call flag("OUTPUT_MAIN_SCF       ", &
               output_main_scf        , df_output_main_scf       )
    call flag("OUTPUT_PROPERTIES     ", &
               output_properties      , df_output_properties     )
    call flag("OUTPUT_READ_INPUT     ", &
               output_read_input      , df_output_read_input     )
    call flag("OUTPUT_WRITE_INPUT    ", &
               output_write_input     , df_output_write_input    )
    call flag("OUTPUT_PVM_MSGTAGS    ", &
               output_pvm_msgtags     , df_output_pvm_msgtags    )
    call flag("OUTPUT_MAIN_DIPOLE    ", &
               output_main_dipole     , df_output_main_dipole    )
    call flag("OUTPUT_EFIELD_CALC    ", &
               output_efield_calc     , df_output_efield_calc    )
    call stop()

    call start("OUTPUT_CONFIG","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_PVMCONFIG ", &
               output_pvmconfig  , df_output_pvmconfig )
    call flag("OUTPUT_OPERATIONS", &
               output_operations , df_output_operations)
    call stop()

    call start("OUTPUT_UNIQUE_ATOMS","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_ATOMS             ", &
               output_atoms              , df_output_atoms             )
    call flag("OUTPUT_SYMADAPT          ", &
               output_symadapt           , df_output_symadapt          )
    call flag("OUTPUT_RENORM            ", &
               output_renorm             , df_output_renorm            )
    call flag("OUTPUT_ORBITALPROJECTIONS", &
               output_orbitalprojections , df_output_orbitalprojections)
    call stop()

    call start("OUTPUT_SCF","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_SCFLOOPS      ", &
               output_scfloops       , df_output_scfloops      )
    call flag("OUTPUT_REOCCUP       ", &
               output_reoccup        , df_output_reoccup       )
    call flag("OUTPUT_FERMI         ", &
               output_fermi          , df_output_fermi         )
    call flag("OUTPUT_FERMI_NEWTON  ", &
               output_fermi_newton   , df_output_fermi_newton  )
    call flag("OUTPUT_CHARGEFIT     ", &
               output_chargefit      , df_output_chargefit     )
    call flag("OUTPUT_CONDITION     ", &
               output_condition      , df_output_condition     )
    call flag("OUTPUT_HAMILTONIAN   ", &
               output_hamiltonian    , df_output_hamiltonian   )
    call flag("OUTPUT_EIGENDATA     ", &
               output_eigendata      , df_output_eigendata     )
    call flag("OUTPUT_EACH_EIGENDATA", &
               output_each_eigendata , df_output_each_eigendata)
    call flag("OUTPUT_EIGEN_STRATEGY", &
               output_eigen_strategy , df_output_eigen_strategy)
    call flag("OUTPUT_DENSMAT       ", &
               output_densmat        , df_output_densmat       )
    call flag("OUTPUT_GRID          ", &
               output_grid           , df_output_grid          )
    call flag("OUTPUT_OVERLAP       ", &
               output_overlap        , df_output_overlap       )
    call flag("OUTPUT_DATA_SAVED    ", &
               output_data_saved     , df_output_data_saved    )
    call flag("OUTPUT_DATA_READ     ", &
               output_data_read      , df_output_data_read     )
    call flag("OUTPUT_SPECTRUM      ", &
               output_spectrum       , df_output_spectrum             )
    call intg("OUTPUT_N_DENSITY_DEV ", &
               output_n_density_dev  , df_output_n_density_dev )
    call intg("OUTPUT_N_COEFF_DEV   ", &
               output_n_coeff_dev    , df_output_n_coeff_dev   )
    call stop()

    call start("OUTPUT_SYM","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_SYMMETRY       ", &
               output_symmetry        , df_output_symmetry       )
    call flag("OUTPUT_NONSYMEQUIVVECS", &
               output_nonsymequivvecs , df_output_nonsymequivvecs)
    call stop()

    call start("OUTPUT_INTEGRAL","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_INT_FITCONTRACT     ", &
               output_int_fitcontract      , df_output_int_fitcontract     )
    call flag("OUTPUT_INT_2C_FIT          ", &
               output_int_2c_fit           , df_output_int_2c_fit          )
    call flag("OUTPUT_INT_SOLHRULES       ", &
               output_int_solhrules        , df_output_int_solhrules       )
    call flag("OUTPUT_INT_PARAMETERS      ", &
               output_int_parameters       , df_output_int_parameters      )
    call flag("OUTPUT_INT_PROGRESS        ", &
               output_int_progress         , df_output_int_progress        )
    call flag("OUTPUT_INT_DETAILEDPROGRESS", &
               output_int_detailedprogress , df_output_int_detailedprogress)
    call flag("OUTPUT_INT_TASKDISTRIBUTION", &
               output_int_taskdistribution , df_output_int_taskdistribution)
    call flag("OUTPUT_INT_QUADRUPELSTORE  ", &
               output_int_quadrupelstore   , df_output_int_quadrupelstore  )
    call flag("OUTPUT_INT_LOOPS           ", &
               output_int_loops            , df_output_int_loops           )
    call flag("OUTPUT_INT_DEEPLOOPS       ", &
               output_int_deeploops        , df_output_int_deeploops       )
    call flag("OUTPUT_INT_DATA            ", &
               output_int_data             , df_output_int_data            )
    call stop()

    call start("OUTPUT_TIMING","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_TIMING_SUMMARY          ", &
               output_timing_summary           , df_output_timing_summary    )
    call flag("OUTPUT_TIMING_DETAILEDSUMMARY  ", &
               output_timing_detailedsummary   , df_output_timing_detailedsum)
    call flag("OUTPUT_TIMING_INTEGRALS        ", &
               output_timing_integrals         , df_output_timing_integrals  )
    call flag("OUTPUT_TIMING_DETAILEDINTEGRALS", &
               output_timing_detailedintegrals , df_output_timing_detailedint)
    call flag("OUTPUT_TIMING_SCFLOOPS         ", &
               output_timing_scfloops          , df_output_timing_scfloops   )
    call flag("OUTPUT_TIMING_SCF              ", &
               output_timing_scf               , df_output_timing_scf        )
    call flag("OUTPUT_TIMING_DETAILEDSCF      ", &
               output_timing_detailedscf       , df_output_timing_detailedscf)
    call flag("OUTPUT_TIMING_POST_SCF          ", &
               output_timing_post_scf          , df_output_timing_post_scf    )
    call flag("OUTPUT_TIMING_DETAILEDPOSTSCF  ", &
               output_timing_detailedpostscf   , df_output_timing_detailedps )
    call flag("OUTPUT_TIMING_SLAVES           ", &
               output_timing_slaves            , df_output_timing_slaves     )
    call flag("OUTPUT_TIMING_INTERRUPTS       ", &
               output_timing_interrupts        , df_output_timing_interrupts )
    call stop()

    call start("OUTPUT_POST_SCF","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_POST_SCF_MAIN", &
               output_post_scf_main , df_output_post_scf_main)
    call flag("OUTPUT_MAIN_GRADIENT", &
               output_main_gradient , df_output_main_gradient)
    call stop()

    call start("OUTPUT_DIPOLE","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_DIPOLE_DETAILED      ", &
               output_dipole_detailed       , df_output_dipole_detailed      )
    call flag("OUTPUT_DIPOLE_TRANSITIONM_F ", &
               output_dipole_transitionm_f  , df_output_dipole_transitionm_f )
    call flag("OUTPUT_DIPOLE_TRANSITIONM_UF", &
               output_dipole_transitionm_uf , df_output_dipole_transitionm_uf)
    call flag("OUTPUT_DIPOLE_SIMOL         ", &
               output_dipole_simol          , df_output_dipole_simol         )
    call flag("OUTPUT_DIPOLE_OPTIMIZER     ", &
               output_dipole_optimizer      , df_output_dipole_optimizer     )

    call stop()

    call start("OUTPUT_RESPONSE","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_RESPONSE_GENERAL ", &
               output_response_general  , df_output_response_general )
    call flag("OUTPUT_RESPONSE_DETAILED", &
               output_response_detailed , df_output_response_detailed)
    call flag("OUTPUT_RESPONSE_DEBUG   ", &
               output_response_debug    , df_output_response_debug)
    call stop()

    call start("OUTPUT_SOLVATION","OUTPUT_WRITE", &
         iounit,operations_echo_input_level)
    call flag("OUTPUT_CAVITY_DATA    ", &
               output_cavity_data     , df_output_cavity_data )
    call flag("OUTPUT_CAVITY_LONG    ", &
               output_cavity_long     , df_output_cavity_long )
    call flag("OUTPUT_SOLV_GRADS     ", &
               output_solv_grads      , df_output_solv_grads  )
    call stop()

  end subroutine output_write
  !*************************************************************

  subroutine output_bcast()
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call output_pack()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call output_unpack()
    endif
  end subroutine output_bcast

  !*************************************************************
  subroutine output_pack()
    ! purpose: packs data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: commpack
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call commpack(output_level,info)
    if (info .ne. 0) call error_handler("output_pack: output_level")
    call commpack(output_scfloops,info)
    if (info .ne. 0) call error_handler("output_pack: output_scfloops")
    call commpack(output_slaveoperations,info)
    if (info .ne. 0) call error_handler("output_pack: output_slaveoperations")
    call commpack(output_main_master,info)
    if (info .ne. 0) call error_handler("output_pack: output_main_master")
    call commpack(output_main_symm,info)
    if (info .ne. 0) call error_handler("output_pack: output_main_symm")
    call commpack(output_main_scf,info)
    if (info .ne. 0) call error_handler("output_pack: output_main_scf")
    call commpack(output_properties,info)
    if (info .ne. 0) call error_handler("output_pack: output_properties")
    call commpack(output_read_input,info)
    if (info .ne. 0) call error_handler("output_pack: output_read_input")
    call commpack(output_write_input,info)
    if (info .ne. 0) call error_handler("output_pack: output_write_input")
    call commpack(output_main_dipole,info)
    if (info .ne. 0) call error_handler("output_pack: output_main_dipole")
    call commpack(output_efield_calc,info)
    if (info .ne. 0) call error_handler("output_pack: output_efield_calc")
    call commpack(output_pvmconfig,info)
    if (info .ne. 0) call error_handler("output_pack: output_pvmconfig")
    call commpack(output_symmetry,info)
    if (info .ne. 0) call error_handler("output_pack: output_symmetry")
    call commpack(output_nonsymequivvecs,info)
    if (info .ne. 0) call error_handler("output_pack: output_nonsymequivvecs")
    call commpack(output_operations,info)
    if (info .ne. 0) call error_handler("output_pack: output_operations")
    call commpack(output_atoms,info)
    if (info .ne. 0) call error_handler("output_pack: output_atoms")
    call commpack(output_symadapt,info)
    if (info .ne. 0) call error_handler("output_pack: output_symadapt")
    call commpack(output_renorm,info)
    if (info .ne. 0) call error_handler("output_pack: output_renorm")
    call commpack(output_orbitalprojections,info)
    if (info .ne. 0) call error_handler("output_pack: output_orbitalprojections")
    call commpack(output_reoccup,info)
    if (info .ne. 0) call error_handler("output_pack: output_reoccup")
    call commpack(output_fermi,info)
    if (info.ne.0) call error_handler("output_pack: output_fermi")
    call commpack(output_fermi_newton,info)
    if (info.ne.0) call error_handler("output_pack: output_fermi_newton")
    call commpack(output_chargefit,info)
    if (info .ne. 0) call error_handler("output_pack: output_chargefit")
    call commpack(output_condition,info)
    if (info .ne. 0) call error_handler("output_pack: output_condition")
    call commpack(output_hamiltonian,info)
    if (info .ne. 0) call error_handler("output_pack: output_hamiltonian")
    call commpack(output_eigendata,info)
    if (info .ne. 0) call error_handler("output_pack: output_eigendata")
    call commpack(output_each_eigendata,info)
    if (info .ne. 0) call error_handler("output_pack: output_each_eigendata")
    call commpack(output_densmat,info)
    if (info .ne. 0) call error_handler("output_pack: output_densmat")
    call commpack(output_grid,info)
    if (info .ne. 0) call error_handler("output_pack: output_grid")
    call commpack(output_overlap,info)
    if (info .ne. 0) call error_handler("output_pack: output_overlap")
    call commpack(output_data_saved,info)
    if (info .ne. 0) call error_handler("output_pack: output_data_saved")
    call commpack(output_data_read,info)
    if (info .ne. 0) call error_handler("output_pack: output_data_read")
    call commpack(output_spectrum,info)
    if (info .ne. 0) call error_handler("output_pack: output_spectrum")
    call commpack(output_n_density_dev,info)
    if (info .ne. 0) call error_handler("output_pack: output_n_density_dev")
    call commpack(output_n_coeff_dev,info)
    if (info .ne. 0) call error_handler("output_pack: output_n_coeff_dev")
    call commpack(output_int_fitcontract,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_fitcontract")
    call commpack(output_int_2c_fit,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_2c_fit")
    call commpack(output_int_solhrules,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_solhrules")
    call commpack(output_int_parameters,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_parameters")
    call commpack(output_int_progress,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_progress")
    call commpack(output_int_detailedprogress,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_detailedprogress")
    call commpack(output_int_taskdistribution,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_taskdistribution")
    call commpack(output_int_quadrupelstore,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_quadrupelstore")
    call commpack(output_int_loops,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_loops")
    call commpack(output_int_deeploops,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_deeploops")
    call commpack(output_int_data,info)
    if (info .ne. 0) call error_handler("output_pack: output_int_data")
    call commpack(output_timing_summary,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_summary")
    call commpack(output_timing_detailedsummary,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_detailedsummary")
    call commpack(output_timing_integrals,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_integrals")
    call commpack(output_timing_detailedintegrals,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_detailedintegrals")
    call commpack(output_timing_scfloops,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_scfloops")
    call commpack(output_timing_scf,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_scf")
    call commpack(output_timing_detailedscf,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_detailedscf")
    call commpack(output_timing_post_scf,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_post_scf")
    call commpack(output_timing_detailedpostscf,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_detailedpostscf")
    call commpack(output_timing_slaves,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_slaves")
    call commpack(output_timing_interrupts,info)
    if (info .ne. 0) call error_handler("output_pack: output_timing_interrupts")
    call commpack(output_post_scf_main,info)
    if (info .ne. 0) call error_handler("output_pack: output_post_scf_main")
    call commpack(output_main_gradient,info)
    if (info .ne. 0) call error_handler("output_pack: output_post_main_gradient")
    call commpack(output_dipole_detailed,info)
    if (info .ne. 0) call error_handler("output_pack: output_dipole_detailed")
    call commpack(output_dipole_transitionm_f,info)
    if (info .ne. 0) call error_handler("output_pack: output_dipole_transitionm_f")
    call commpack(output_dipole_transitionm_uf,info)
    if (info .ne. 0) call error_handler("output_pack: output_dipole_transitionm_uf")
    call commpack(output_dipole_simol,info)
    if (info .ne. 0) call error_handler("output_pack: output_dipole_simol")
    call commpack(output_dipole_optimizer,info)
    if (info .ne. 0) call error_handler("output_pack: output_dipole_optimizer")
    call commpack(output_response_general,info)
    if (info .ne. 0) call error_handler("output_pack: output_response_general")
    call commpack(output_response_detailed,info)
    if (info .ne. 0) call error_handler("output_pack: output_response_detailed")
    call commpack(output_response_debug,info)
    if (info .ne. 0) call error_handler("output_pack: output_response_debug")
    call commpack(output_cavity_data,info)
    if (info .ne. 0) call error_handler("output_pack: output_cavity_data")
    call commpack(output_cavity_long,info)
    if (info .ne. 0) call error_handler("output_pack: output_cavity_long")
    call commpack(output_solv_grads,info)
    if (info .ne. 0) call error_handler("output_pack: output_solv_grads")
  end subroutine output_pack
  !*************************************************************


  !*************************************************************
  subroutine output_unpack()
    ! purpose: unpacks data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: communpack
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call communpack(output_level,info)
    if (info .ne. 0) call error_handler("output_unpack: output_level")
    ! restore the output leve dependent default values
    call output_set_defaults

    call communpack(output_scfloops,info)
    if (info .ne. 0) call error_handler("output_unpack: output_scfloops")
    call communpack(output_slaveoperations,info)
    if (info .ne. 0) call error_handler("output_unpack: output_slaveoperations")
    call communpack(output_main_master,info)
    if (info .ne. 0) call error_handler("output_unpack: output_main_master")
    call communpack(output_main_symm,info)
    if (info .ne. 0) call error_handler("output_unpack: output_main_symm")
    call communpack(output_main_scf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_main_scf")
    call communpack(output_properties,info)
    if (info .ne. 0) call error_handler("output_unpack: output_properties")
    call communpack(output_read_input,info)
    if (info .ne. 0) call error_handler("output_unpack: output_read_input")
    call communpack(output_write_input,info)
    if (info .ne. 0) call error_handler("output_unpack: output_write_input")
    call communpack(output_main_dipole,info)
    if (info .ne. 0) call error_handler("output_unpack: output_main_dipole")
    call communpack(output_efield_calc,info)
    if (info .ne. 0) call error_handler("output_unpack: output_efield_calc")
    call communpack(output_pvmconfig,info)
    if (info .ne. 0) call error_handler("output_unpack: output_pvmconfig")
    call communpack(output_symmetry,info)
    if (info .ne. 0) call error_handler("output_unpack: output_symmetry")
    call communpack(output_nonsymequivvecs,info)
    if (info .ne. 0) call error_handler("output_unpack: output_nonsymequivvecs")
    call communpack(output_operations,info)
    if (info .ne. 0) call error_handler("output_unpack: output_operations")
    call communpack(output_atoms,info)
    if (info .ne. 0) call error_handler("output_unpack: output_atoms")
    call communpack(output_symadapt,info)
    if (info .ne. 0) call error_handler("output_unpack: output_symadapt")
    call communpack(output_renorm,info)
    if (info .ne. 0) call error_handler("output_unpack: output_renorm")
    call communpack(output_orbitalprojections,info)
    if (info .ne. 0) call error_handler("output_unpack: output_orbitalprojections")
    call communpack(output_reoccup,info)
    if (info .ne. 0) call error_handler("output_unpack : output_reoccup")
    call communpack(output_fermi,info)
    if(info.ne.0) call error_handler("output_unpack : output_fermi")
    call communpack(output_fermi_newton,info)
    if(info.ne.0) call error_handler("output_unpack : output_fermi_newton")
    call communpack(output_chargefit,info)
    if (info .ne. 0) call error_handler("output_unpack : output_chargefit")
    call communpack(output_condition,info)
    if (info .ne. 0) call error_handler("output_unpack: output_condition")
    call communpack(output_hamiltonian,info)
    if (info .ne. 0) call error_handler("output_unpack: output_hamiltonian")
    call communpack(output_eigendata,info)
    if (info .ne. 0) call error_handler("output_unpack: output_eigendata")
    call communpack(output_each_eigendata,info)
    if (info .ne. 0) call error_handler("output_unpack: output_each_eigendata")
    call communpack(output_densmat,info)
    if (info .ne. 0) call error_handler("output_unpack: output_densmat")
    call communpack(output_grid,info)
    if (info .ne. 0) call error_handler("output_unpack: output_grid")
    call communpack(output_overlap,info)
    if (info .ne. 0) call error_handler("output_unpack: output_overlap")
    call communpack(output_data_saved,info)
    if (info .ne. 0) call error_handler("output_unpack: output_data_saved")
    call communpack(output_data_read,info)
    if (info .ne. 0) call error_handler("output_unpack: output_data_read")
    call communpack(output_spectrum,info)
    if (info .ne. 0) call error_handler("output_unpack: output_spectrum")
    call communpack(output_n_density_dev,info)
    if (info .ne. 0) call error_handler("output_unpack: output_n_density_dev")
    call communpack(output_n_coeff_dev,info)
    if (info .ne. 0) call error_handler("output_unpack: output_n_coeff_dev")
    call communpack(output_int_fitcontract,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_fitcontract")
    call communpack(output_int_2c_fit,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_2c_fit")
    call communpack(output_int_solhrules,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_solhrules")
    call communpack(output_int_parameters,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_parameters")
    call communpack(output_int_progress,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_progress")
    call communpack(output_int_detailedprogress,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_detailedprogress")
    call communpack(output_int_taskdistribution,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_taskdistribution")
    call communpack(output_int_quadrupelstore,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_quadrupelstore")
    call communpack(output_int_loops,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_loops")
    call communpack(output_int_deeploops,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_deeploops")
    call communpack(output_int_data,info)
    if (info .ne. 0) call error_handler("output_unpack: output_int_data")
    call communpack(output_timing_summary,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_summary")
    call communpack(output_timing_detailedsummary,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_detailedsummary")
    call communpack(output_timing_integrals,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_integrals")
    call communpack(output_timing_detailedintegrals,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_detailedintegrals")
    call communpack(output_timing_scfloops,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_scfloops")
    call communpack(output_timing_scf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_scf")
    call communpack(output_timing_detailedscf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_detailedscf")
    call communpack(output_timing_post_scf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_post_scf")
    call communpack(output_timing_detailedpostscf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_detailedpostscf")
    call communpack(output_timing_slaves,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_slaves")
    call communpack(output_timing_interrupts,info)
    if (info .ne. 0) call error_handler("output_unpack: output_timing_interrupts")
    call communpack(output_post_scf_main,info)
    if (info .ne. 0) call error_handler("output_unpack: output_post_scf_main")
    call communpack(output_main_gradient,info)
    if (info .ne. 0) call error_handler("output_unpack: output_main_gradient")
    call communpack(output_dipole_detailed,info)
    if (info .ne. 0) call error_handler("output_unpack: output_dipole_detailed")
    call communpack(output_dipole_transitionm_f,info)
    if (info .ne. 0) call error_handler("output_unpack: output_dipole_transitionm_f")
    call communpack(output_dipole_transitionm_uf,info)
    if (info .ne. 0) call error_handler("output_unpack: output_dipole_transitionm_uf")
    call communpack(output_dipole_simol,info)
    if (info .ne. 0) call error_handler("output_unpack: output_dipole_simol")
    call communpack(output_dipole_optimizer,info)
    if (info .ne. 0) call error_handler("output_unpack: output_dipole_optimizer")
    call communpack(output_response_general,info)
    if (info .ne. 0) call error_handler("output_unpack: output_response_general")
    call communpack(output_response_detailed,info)
    if (info .ne. 0) call error_handler("output_unpack: output_response_detailed")
    call communpack(output_response_debug,info)
    if (info .ne. 0) call error_handler("output_unpack: output_response_debug")
    call communpack(output_cavity_data,info)
    if (info .ne. 0) call error_handler("output_unpack: output_cavity_data")
    call communpack(output_cavity_long,info)
    if (info .ne. 0) call error_handler("output_unpack: output_cavity_long")
    call communpack(output_solv_grads,info)
    if (info .ne. 0) call error_handler("output_unpack: output_solv_grads")
  end subroutine output_unpack
  !*************************************************************


!--------------- End of module ----------------------------------
end module output_module

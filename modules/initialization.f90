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
!==============================================================================+
! Public interface of module
!==============================================================================+
module initialization
!
!  Purpose: different initialization purposes that need to run
!           in a parallel context
!
!           * substitute the old send/receive_initialisation
!             constructs by comm calls
!
!
!  Subroutine called by: main_master/ main_slave
!
!
!  Author: TMS
!  Date:  9/11
!
!==============================================================================+
! End of public interface of module
!==============================================================================+
!
!------------------------------------------------------------------------------+
! Modifications
!------------------------------------------------------------------------------+
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!------------------------------------------------------------------------------+

!------------ Modules used ----------------------------------------------------+
# include "def.h"
use type_module, only: IK => i4_kind, RK => r8_kind

implicit none
private

!------------ public functions and subroutines --------------------------------+

public :: initialize_environment!(), executed once before the input is read
public :: initialize_with_input!(), executed right after the input is read
public :: finalize_geometry!(), executed before leaving (or repeating) geometry loop
public :: finalize!(), executed before finializing the communication layer and exiting

!------------ Subroutines -----------------------------------------------------+
contains

  subroutine initialize_environment()
    !
    ! Executed  by all  workers before  the data  from the  input file
    ! becomes available. Only  the data that can be  obtained from the
    ! process environment and command line arguments may be used here.
    !
    use comm, only: comm_rank
    use time_module,  only: time_setup  ! timing routines
    use timer_module, only: timer_setup ! timing database
    use filename_module, only: filename_setup
    use iounitadmin_module, only: open_special_units
    use error_module, only: init_error_module
    implicit none
    ! *** end of interface ***

    integer(IK) :: rank

    rank = comm_rank()

    !
    ! Tell my index to set MyID prefix to many debug outputs:
    !
    call init_error_module(rank + 1)

    !
    ! Start the clock:
    !
    call time_setup()
    call timer_setup()

    !
    ! Read  setup from  environment variables,  this may  not  work on
    ! remote nodes,  so let the  master process do this  and broadcast
    ! results to all other workers.
    !
    ! OTOH, it the buisiness of filename_module how this is achieved:
    !
    DPRINT 'initialize_environment: call filename_setup()'
    call filename_setup()
    DPRINT 'initialize_environment: .'

    !
    ! Only master is supposed to use output_unit/trace_unit:
    !
    DPRINT 'initialize_environment: call open_special_units()'
    call open_special_units (rank)
    DPRINT 'initialize_environment: .'
  end subroutine initialize_environment

  subroutine initialize_with_input()
    !
    ! Executed in parallel context by  all workers soon after the data
    ! from input  file was  loaded (by master).  This is a  good entry
    ! point  for distributing initial  options. Indeed  a bulk  of the
    ! work done  here is plain broadcasting the  structures created on
    ! master and filled with input data.
    !
    ! Note  that just as  the reading  of the  input is  repeated evey
    ! geometry iteration, so is this the execution of this sub.
    !
    use comm, only: comm_rank
    use operations_module, only: operations_write_input_slave, operations_bcast, &
        operations_post_scf, operations_response
#ifdef WITH_MOLMECH
    use operations_module, only: operations_mol_mech
#endif
#ifdef WITH_EPE
    use operations_module, only: operations_epe_lattice
#endif
    use output_module, only: output_bcast, output_int_solhrules
    use iounitadmin_module, only: output_unit
    use machineparameters_module, only: machineparameters_bcast
    use symmetry_data_module, only: symmetry_data_setup, &
        symmetry_data_write_formatted
    use fit_coeff_module, only: fit_coeff_setup
    use options_module, only: options_data_bcast, options_spin_orbit
    use unique_atom_module, only: unique_atom_lmax_all
    use unique_atom_methods, only: unique_atom_setup
#ifdef WITH_EFP
    use efp_module,               only: efp_bcast
#endif
    use pointcharge_module,       only: pointcharge_bcast                      &
                                      , unique_timp_symadapt_bcast             &
                                      , N_moving_unique_timps
    use point_dqo_module,         only: external_centers_bcast
    use induced_dipoles_module,   only: pol_centers_bcast
#ifdef WITH_EPE
    use ewaldpc_module,           only: ewpc_bcast
#endif
    use orbitalprojection_module, only: orbitalprojection_setup
    use xc_cntrl,                 only: xc_input_bcast
    use spin_orbit_module,        only: spin_orbit_input_bcast
    use xcmda_hamiltonian,        only: mda_options_bcast
    use grid_module,              only: grid_bcast
    use ph_cntrl,                 only: post_scf_input_bcast
#ifdef WITH_RESPONSE
    use response_module,          only: response_input_bcast
#endif
    use orbital_plot_module,      only: orbital_plot_bcast
#ifdef WITH_ERI4C
    use eri4c_options,            only: eri4c_options_bcast
#endif
#ifdef WITH_DFTPU
    use dft_plus_u_module,        only: dft_plus_u_proj_bcast
#endif
    use solhrules_module, only: solhrules_setup, solhrules_print
    use interfaces, only: IMAST
    use occupation_module, only: get_n_elec
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc, &
        fit_coeff_dimensions_calc
    use efield_module, only: efield_print
    use symm_adapt_module, only: symm_adapt_module_setup
    implicit none
    ! *** end of interface ***

    real(RK) :: n_elec
    integer(IK) :: rank

#ifdef WITH_MOLMECH
    if(operations_mol_mech) return
#endif
#ifdef WITH_EPE
    if(operations_epe_lattice) return
#endif
    rank = comm_rank()

    !
    ! FIXME: for some reason broadcasting was sometimes disabled.
    !        Does it still need to be? In any case, abort until we get
    !        a test case. See also read_input() for the logic.
    !
#ifdef WITH_EPE
    ASSERT(.not.operations_epe_lattice)
#endif
#ifdef WITH_MOLMECH
    ASSERT(.not.operations_mol_mech)
#endif

    !
    ! Some of those *_bcast() subroutines below send/receive
    !
    !   msgtag_packed_message
    !
    ! and invoke symmetric combinaitons of pack/unpack.
    !
    ! TODO: transform rest of the pack/unpack subroutines into bcast
    ! versions
    !

    !
    ! OPERATIONS_MODULE does not use anything
    !
    call operations_bcast()

    !
    ! MACHINEPARAMETERS_MODULE uses
    !   OPERATIONS_MODULE
    !
    call machineparameters_bcast()

    !
    ! OUTPUT_MODULE uses
    !   OPERATIONS_MODULE
    !
    call output_bcast()

    !
    ! OPTIONS_MODULE uses
    !   OPERATIONS_MODULE
    !
    call options_data_bcast()

    !
    ! SYMMETRY_DATA_MODULE cannot use anything
    !
    ! set up mapping of irrep and partner indizes to combined index as
    ! used for dipoles
    !
    call symmetry_data_setup(options_spin_orbit)
    !
    ! SPIN_ORBIT_MODULE uses
    !   OPERATIONS_MODULE
    !   OPTIONS_MODULE
    !
    call spin_orbit_input_bcast()

    !
    ! UNIQUE_ATOM_MODULE uses
    !   OPTIONS_MODULE
    !   SYMMETRY_DATA_MODULE
    !   SPIN_ORBIT_MODULE
    !
    ! prepare data in unique_atom_module
    ! calculate dimensions of orbital basis for all irreps
    !   and charge + exchange fitfunction basis
    ! elliminate irreps with  basis dimension 0,
    ! also initiates setup of the modules it USEs:
    !
    call unique_atom_setup()

    !
    ! SYMM_ADAPT_MODULE uses
    !   OPTIONS_MODULE
    !   UNIQUE_ATOM_MODULE
    !   SPIN_ORBIT_MODULE
    !
    call symm_adapt_module_setup()

    !
    ! ORBITALPROJECTION_MODULE uses
    !   OUTPUT_MODULE
    !
    ! calculate indices describing projection of symmetryadapted basis
    ! functions to metaindex used in scf_part
    !
    call orbitalprojection_setup()

    !
    ! FIT_COEFF_MODULE uses
    !   OPERATIONS_MODULE
    !   unique_atom_module (indirectly?)
    !
    call fit_coeff_setup()

#ifdef WITH_EFP
    call efp_bcast()
#endif
    !
    ! bcast point charges
    call pointcharge_bcast()
    if( N_moving_unique_timps /= 0 ) call unique_timp_symadapt_bcast()
    !
    ! bcast external points
    call external_centers_bcast()
    !
    ! bcast external polarizable point data
    call pol_centers_bcast()
    !
#ifdef WITH_EPE
    ! bcast ewald point charges
    call ewpc_bcast()
#endif

    ! bcast xc_control input
    call xc_input_bcast()
    !
    ! bcast mda options and parameter
    call mda_options_bcast()
    !
    ! bcast grid input
    call grid_bcast( (operations_post_scf .or. operations_response) )
    !
    ! bcast post scf input
    call post_scf_input_bcast()
    !
#ifdef WITH_RESPONSE
    ! bcast response module input
    call response_input_bcast()
#endif
    ! bcast orbital_plot_input
    call orbital_plot_bcast()
    !
#ifdef WITH_ERI4C
    call eri4c_options_bcast()
#endif
    !
#ifdef WITH_DFTPU
    call dft_plus_u_proj_bcast()
#endif

    !
    ! Below this point it is not about broadcasting the data anymore
    ! ...
    !
    ! Based on the input data some immediate preparations are
    ! meaningfull at this point.
    !
    ! FIXME: write and share only information if fitted coulomb is
    ! employed

    !
    ! Was done before in main_slave, right after call to
    ! initialize_with_input():
    !
    if ( rank /= 0 ) then
        ! write input file if desired
        if ( operations_write_input_slave ) then
            call write_input()
        endif

    endif

    !
    ! Initialize product and differential rules for solid harmonics:
    !
    call solhrules_setup(unique_atom_lmax_all)

    !
    ! Write some setup informations of general interest:
    !
    if ( rank == 0 ) then
        call symmetry_data_write_formatted(output_unit)

        write(output_unit,*)
        write(output_unit,*) "Dimension of charge   basis: ", fit_coeff_n_ch()
        write(output_unit,*) "Dimension of exchange basis: ", fit_coeff_n_xc()
        write(output_unit,*)
        call get_n_elec(n_elec)
        write(output_unit,*) "Total number of electrons: ", n_elec
        write(output_unit,*)

        !
        ! Output of applied external electrival field
        !
        call efield_print(output_unit)

        if ( output_int_solhrules ) &
            call solhrules_print(output_unit)
    endif
    !
  end subroutine initialize_with_input

  subroutine finalize_geometry()
    !
    ! Executed in  parallel context by  all workers after the  most of
    ! the work has beed done  for the current geometry, but before the
    ! task  loop   proceeds  to  another  geometry   (or  leaves  upon
    ! convergence).
    !
    ! For  most  things  there  is  no  point  to  allocate/deallocate
    ! anything  between  geometry  loops. However,  historically  many
    ! things  are freed  just to  be  allocated again.   See (not  yet
    ! existent) initalize_geometry().
    !
    ! Note that this sub is called every geometry iteration.
    !
    use comm, only: comm_rank
    use operations_module, only: operations_scf, operations_integral, &
        operations_gradients, operations_symm, &
        operations_solvation_effect, operations_gx_test
    use occupation_module, only: occupation_shutdown, dealloc_occ_num
    use occupied_levels_module, only: eigvec_occ_dealloc
    use fit_coeff_module, only: fit_coeff_shutdown
    use mat_charge_module, only: free_mat_charge
    use eigen_data_module, only: eigen_data_free
    use orbitalprojection_module, only: orbitalprojection_close
    use symmetry, only: main_symm, MAIN_SYMM_DONE
    use symmetry_data_module, only: symmetry_data_close
    use solv_electrostat_module, only: shutdown_solvation
    use virtual_levels_module,only:  eigvec_vir_dealloc, print_viralloc
    use grid_module, only: print_alloc_grid
    use cpksdervs_matrices, only: cpks_gmat
    use integralpar_module, only: integralpar_cpksdervs
    use interfaces, only: IPARA
    use integralstore_module, only: integralstore_deallocate
    use pointcharge_module, only: pointcharge_close
    use point_dqo_module, only: dealloc_dqo, dealloc_center_inform, &
        moving_X_centers, moving_R_centers
    use unique_atom_methods, only: unique_atom_close
#ifdef WITH_EPE
    use ewaldpc_module, only: ewpc_close
    use epe_module, only: print_epemod_alloc
#endif
#ifdef WITH_EFP
    use efp_rep_module, only: dealloc_repf, dealloc_repf_inform
#endif
    use induced_dipoles_module, only: dealloc_id, dealloc_Pol_center_inform
#ifdef WITH_EFP
    use efp_rep_module, only: dealloc_repf, dealloc_repf_inform
#endif


    use orbital_module, only: print_orbmod_alloc
    use density_calc_module, only: print_alloc_density_calc
    use population_module, only: population_close
    use symm_adapt_module, only: symm_adapt_module_close
    implicit none
    ! *** end of interface ***

    integer(IK) :: rank, memstat

    rank = comm_rank()

    !
    ! Some of this code vas moved from main_gradient():
    !
    if (operations_gradients) then

        !
        ! New place for dervs related deallocations
        !
        if (integralpar_cpksdervs) then

            call eigvec_vir_dealloc (IPARA)

            ! FIXME: find matching for master and do it there, yet
            !        better don duplicate functionality of
            !        occupation_module
            if (rank /= 0) call dealloc_occ_num()

            ! FIXME: so far does not harm if not allocated:
            call integralstore_deallocate()

            !
            ! FIXME: make a module sub that does all of finalizaton
            !        for CPKS, why manimulating module varibale from
            !        everywhere?
            if (allocated (cpks_gmat)) then
                deallocate (cpks_gmat, stat=memstat)
                ASSERT(memstat==0)
            endif
        endif

#ifdef WITH_EPE
        ! both were always called after unique_atom_close in main_slave
        call ewpc_close()
        call print_epemod_alloc()
#endif

        !
        ! FIXME: These may be related to solvation, so why not making
        !        a "call solvation_shutdown()"?
        !
        call dealloc_dqo()
        if ( moving_X_centers .or. moving_R_centers ) &
            call dealloc_center_inform()
        call dealloc_id()
        call dealloc_Pol_center_inform()

#ifdef WITH_EFP
        if ( rank == 0 ) then
            call dealloc_repf()
            if(moving_R_centers) call dealloc_repf_inform()
        end if
#endif

        call eigvec_occ_dealloc() ! also eigenvalues
    endif ! operations_gradients

    !
    ! Clean up population_module (uses quite a lot):
    !
    call population_close()

    !
    ! Various shutdown and deallocation work:
    ! FIXME: why only if operations_scf?
    !
    if ( operations_scf ) then
        call occupation_shutdown()
    endif

    !
    ! All workers call fit_coeff_shutdown()
    !
    if (operations_integral.or.operations_scf) then
        ! no comm:
        call fit_coeff_shutdown()
    endif

    !
    ! no comm, idempotent:
    !
    call free_mat_charge()

    if (operations_scf) then
       ! no comm:
       call eigen_data_free()
    endif

    !
    ! orbitalprojection_module uses
    !   unique_atom_module
    !
    call orbitalprojection_close()

    !
    ! SYMM_ADAPT_MODULE uses
    !   OPTIONS_MODULE
    !   UNIQUE_ATOM_MODULE
    !   SPIN_ORBIT_MODULE
    !
    call symm_adapt_module_close()

    !
    ! Clean up unique_atom_module/unique_atom_methods, also initiates
    ! cleanup of the modules it USEs:
    !
    call unique_atom_close()

    !
    ! Clean up pointcharge_module:
    !
    call pointcharge_close()

    !
    ! main_symm() is not a parallel code, though we call it everywhere
    ! instead of broadcasting symmetrization coefficients:
    !
    if (operations_symm) then
       ! no comm:
       call main_symm (MAIN_SYMM_DONE)
    endif

    ! FIXME: why only then?
    if ( operations_gx_test ) then
        ! no comm:
        call symmetry_data_close()
    end if

    if (operations_solvation_effect) then
        ! no comm:
        call shutdown_solvation()
    endif

    ! debug, probably:
    call print_orbmod_alloc()
    call print_alloc_density_calc()

    ! no comm, debug only:
    call print_viralloc()

    ! no comm, debug only:
    call print_alloc_grid()
  end subroutine finalize_geometry

  subroutine finalize()
    !
    ! Executed in  parallel context  by all workers  after all  of the
    ! work has beed done, just before exiting.  One should assume that
    ! we  either  terminate  or  proceed with  a  completly  different
    ! calculation --- with different input, number of atoms, bases and
    ! methods. So  do not assume we  can re-use any  of the structures
    ! from the last run --- clean everything.
    !
    ! Note that this sub is usually called just once.
    !
    use comm, only: comm_rank
    use shgi_adkh, only: shgi_adkh_close
#ifdef WITH_SCHEDEIG
    use se_eigen_module, only: se_eigen_close
#endif
    use bounds_module, only: bounds_free
    use iounitadmin_module, only: close_special_units
    use solhrules_module, only: solhrules_free
    use orbital_module, only: orbital_shutdown
    use grid_module, only: grid_close
    implicit none
    ! *** end of interface ***

#ifdef WITH_SCHEDEIG
    !
    ! Clean up the parallel eigensolver subsystem:
    !
    call se_eigen_close()
#endif

    !
    ! Deallocate cached atomic transformatin matrices:
    !
    call shgi_adkh_close()

    !
    ! Free space used for calculating orbitals:
    !
    call orbital_shutdown()

    !
    ! Close   output    file.   See   the    call   to   corresponding
    ! open_special_units() above:
    !
    call close_special_units()

    !
    ! Free structs of solhrules_module:
    !
    call solhrules_free()

    !
    ! Free structs of bounds_module:
    !
    call bounds_free()

    !
    ! Clean up grid_module. "False" tells not to keep anything allocated:
    !
    call grid_close(.false.)
  end subroutine finalize

!*******************************************************************************
end module initialization

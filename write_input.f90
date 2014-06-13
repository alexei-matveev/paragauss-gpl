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
subroutine write_input
!----------------------------------------------------------------
!
!  Purpose: Writing an input file that can be interpreted by 
!           read_input() to file "$TTFSOUT/input.out".
!           The write routines of various modules are called in 
!           the same order as the read routines in read_input.
!
!  Subroutine called by: main_master
!
!  Author: TB
!  Date: 10/95
!
!================================================================
! End of public interface of module
!================================================================
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
  use type_module ! type specification parameters
  use unique_atom_methods
  use symmetry_data_module
  use iounitadmin_module
  use filename_module, only: outfile
  use operations_module
  use comm_module, only: comm_i_am_master
  use options_module, only: options_write_input
  use mixing_module, only: mixing_write_input
  use diis_fock_module, only: diis_write_input
  use fermi_module, only: fermi_write_input
  use convergence_module, only: convergence_write
  use output_module
  use machineparameters_module, only: machineparameters_write
  use grid_module, only: grid_write_ph, grid_write
  use ph_cntrl, only: post_scf_write_input
  use xc_cntrl,    only: xc_write_input
  use xcmda_hamiltonian, only: mda_options_write
  use occupation_module, only: occupation_write
  use properties_module, only: properties_write
  use frag_orb_analysis_module, only: frag_orb_analysis_write
  use orbital_plot_module, only: orbital_plot_write
  use population_module, only: population_write
  use pointcharge_module, only: pointcharge_write
  use dipole_module, only: dipole_write
#ifdef WITH_GTENSOR
  use hfc_module, only: magnetic_write
#endif
  use efield_module, only: efield_write
#ifdef WITH_RESPONSE
  use response_module, only: response_write_input
#endif
  use spin_orbit_module, only: spin_orbit_write_input
#ifdef WITH_EPE
  use ewaldpc_module, only: epe_write  !!!!!!!!!!!!!!!!AS
#endif
  use potential_calc_module, only: poten_calc_write !!!!!!!!!!!!!!!AS
  use solv_cavity_module, only: solvation_write
  use disp_rep_module, only: disp_rep_write
  use solv_charge_mixing_module, only: write_solv_mix
#ifdef WITH_MOLMECH
  use qmmm_interface_module, only: qmmm_write !!!!!!!!!!AS
#endif
#ifdef WITH_EFP
  use efp_module, only: efp_write_input
#endif
#if WITH_ERI4C == 1
  use eri4c_options, only: eri4c_write_input
#endif

  implicit none

  integer      :: unit

  ! FIXME: is it OK so?:
#ifndef NEW_EPE
  if(operations_epe_lattice) RETURN
#endif
#ifdef WITH_MOLMECH
  if( operations_mol_mech ) RETURN
#endif

  if (output_write_input) call write_to_output_units("write input: ")

  ! open file
  unit = openget_iounit(file=trim(outfile("input.out")))


  if (output_write_input) call write_to_output_units("write input: operations")
  call operations_write(unit)

  if (output_write_input) call write_to_output_units("write input: output")
  call output_write(unit)

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: population")
     call population_write(unit)
  else
     write(unit,'(a/)')" # Namelist POPULATION_LEVEL : not passed to the slaves"
  endif

#ifdef WITH_MOLMECH
  if (output_write_input) call write_to_output_units("write input: qmmm_write")
  call qmmm_write(unit)
#endif


  if (output_write_input) call write_to_output_units("write input: main_options")
  call options_write_input(unit)

  if ( operations_old_input ) then
     call error_handler("write_input: operations_old_input outdated")
  endif
     
     if (output_write_input) call write_to_output_units("write input: symmetry_group")
     call symmetry_data_group_write(unit)

     if (output_write_input) call write_to_output_units("write input: unique_atom_unique")
     call unique_atom_write(unit)
!!$  else
!!$     if (output_write_input) call write_to_output_units("write input: symmetry_data")
!!$     call symmetry_data_write(unit)
!!$
!!$     if (output_write_input) call write_to_output_units("write input: unique_atom")
!!$     call unique_atom_write(unit)
!!$  endif

#ifdef WITH_MOLMECH
  if(.not. operations_qm_mm_new) then   
#endif  
     if (output_write_input) call write_to_output_units("write input: pointcharge")
     call pointcharge_write(unit)
#ifdef WITH_MOLMECH
  end if
#endif  

#ifdef WITH_EFP
  if (output_write_input) call write_to_output_units("write input: efp_data")
  call efp_write_input(unit)
#endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: mixing")
     call mixing_write_input(unit)
  else
     write(unit,'(a/)')" # Namelist MIXING : not passed to the slaves"
  endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: diis")
     call diis_write_input(unit)
  else
     write(unit,'(a/)')" # Namelist DIIS : not passed to the slaves"
  endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: fermi")
     call fermi_write_input(unit)
  else
     write(unit,'(a/)')" # Namelist FERMI : not passed to the slaves"
  endif
  
  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: convergence")
     call convergence_write(unit)
  else
     write(unit,'(a/)')" # Namelist CONVERGENCE_LIST : not passed to the slaves"
  endif

  if (output_write_input) call write_to_output_units("write input: xc")
  call xc_write_input(unit)

#if WITH_ERI4C == 1
  if (output_write_input) call write_to_output_units("write input: ERI4C")
  call eri4c_write_input(unit)
#endif

  if (output_write_input) call write_to_output_units("write input: PostSCF")
  call post_scf_write_input(unit)

#ifdef WITH_RESPONSE
  if (output_write_input) call write_to_output_units("write input: Response")
  call response_write_input(unit)
#endif

  if (output_write_input) call write_to_output_units("write input: spin_orbit")
  call spin_orbit_write_input(unit)

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: occupation")
     call occupation_write(unit)
  else
     write(unit,'(a/)')" # Namelist OCCUPATION : not passed to the slaves"
  endif

  if (output_write_input) call write_to_output_units("write input: mda_options")
  call mda_options_write(unit)

#ifdef WITH_EPE
  if (output_write_input) call write_to_output_units("write input: epe_data")
  call epe_write(unit)
#endif

  if (output_write_input) call write_to_output_units("write input: poten_calc_write")
  call poten_calc_write(unit)
  if(operations_solvation_effect) then
     if (output_write_input) call write_to_output_units("write input: solvation")
     call solvation_write(unit)
     if (output_write_input) call write_to_output_units("write input:write_solv_mix")
     call write_solv_mix(unit)
     if (output_write_input) call write_to_output_units("write input:disper_repal")
     call disp_rep_write(unit)
  endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: properties")
     call properties_write(unit)
  else
     write(unit,'(a/)')" # Namelist PROPERTIES : not passed to the slaves"
  endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units("write input: frag_orb_analysis")
     call frag_orb_analysis_write(unit)
  else
     write(unit,'(a/)')" # Namelist FRAG_ORB_ANALYSIS : not passed to the slaves"
  endif

  if (output_write_input) call write_to_output_units("write input: orbital_plot")
  call orbital_plot_write(unit)

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units( &
          "write input: dipole transition moments")
     call dipole_write(unit)
  else
     write(unit,'(a/)')" # Namelist DIPOLE_TRANSITIONMOMENTS : "// &
           "not passed to the slaves"
  endif

#ifdef WITH_GTENSOR
  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units( &
          "write input: magnetic_properties")
     call magnetic_write(unit)
  else
     write(unit,'(a/)')" # Namelist MAGNETIC_PROPERTIES : "// &
          "not passed to the slaves"
  endif
#endif

  if (comm_i_am_master()) then
     if (output_write_input) call write_to_output_units( &
          "write input: external electrical field")
     call efield_write(unit)
  else
     write(unit,'(a/)')" # Namelist EFIELD : "// &
           "not passed to the slaves"
  endif

  if (output_write_input) call write_to_output_units("write input: machineparameters")
  call machineparameters_write(unit)

  if (output_write_input) call write_to_output_units("write input: grid")
#ifndef NEW_EPE
  if(.not.operations_epe_lattice) call grid_write(unit)
#else
  call grid_write(unit)
#endif

  if (output_write_input) call write_to_output_units("write input: unique_atom_basis")
  call unique_atom_write_basis(unit)

  if ( operations_old_input ) then
     if (output_write_input) call write_to_output_units("write input: symadapt")
     call unique_atom_symadapt_write(unit)

     if (output_write_input) call write_to_output_units("write input: symeqiv")
     call unique_atom_symequiv_write(unit)
  endif

  if ( (.not. operations_integral) .and. operations_scf ) then
     call error_handler("write_input: outdated branch")
!!$     if (output_write_input) call write_to_output_units("write input: renorm")
!!$     call unique_atom_renorm_write(unit)
  endif

  ! to really skip potential backspaced lines
  endfile(unit)

  ! close file
  call returnclose_iounit(unit)


  if (output_write_input) call write_to_output_units("write input: done")


end subroutine write_input

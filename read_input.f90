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
subroutine read_input(loop)
!
!  Reading  in  input file  and  performing  consistency checks.   The
!  input_module is initialised and the  input file is opened, then the
!  read  routines of  various  modules  are called  to  read the  data
!  beloging  to their  modules.  After  that routines  for consistency
!  checks are  called.  Then the  input file is closed.   All routines
!  reading the input file should use input_module: input_read(line) to
!  obtain string to read  from.  inut_error(message) to handle reading
!  errors An output can be automatically produced by output_write that
!  contains call to write  routines corresponding to the read routines
!  called here.
!
!  Subroutine called by: main_master
!
!
!  Author: TB
!  Date: 10/95
!
!================================================================
! End of public interface of module
!================================================================
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

!------------ Modules used --------------------------------------
#include "def.h"
use type_module ! type specification parameters
use input_module
use filename_module, only: inpfile, input_name
use iounitadmin_module
use options_module
use symmetry_data_module, only: symmetry_data_group_read
use unique_atom_module, only: N_unique_atoms, unique_atoms
use unique_atom_methods, only: unique_atom_read, unique_atom_read_basis
#ifdef WITH_GUILE
use unique_atom_methods, only: unique_atom_find_basis
#endif
use output_module
use operations_module
use mixing_module, only: mixing_read
use diis_fock_module, only: diis_read_input
use fermi_module, only: fermi_read, check_occ_fermi
use convergence_module, only: convergence_read
use grid_module, only: grid_read, grid_read_ph, grid_copy_to_ph
use machineparameters_module, only: machineparameters_read
use xc_cntrl, only: xc_read_input
#ifdef WITH_DFTPU
use dft_plus_u_module, only : dft_plus_u_read_input                            &
                            , dft_plus_u_mo_read_input
#endif
use xcmda_hamiltonian, only: mda_options_read
use occupation_module, only: occupation_read, count_electrons
use population_module, only: population_read
use pointcharge_module, only: pointcharge_read
use dipole_module, only: dipole_read
#ifdef WITH_GTENSOR
use hfc_module, only: magnetic_read
#endif
use efield_module, only: efield_read, efield_applied
use properties_module, only: properties_read
use frag_orb_analysis_module, only: frag_orb_analysis_read
use orbital_plot_module, only: orbital_plot_read
#ifdef WITH_RESPONSE
use response_module, only: response_read_input
#endif
use spin_orbit_module, only: spin_orbit_read_input

#ifdef WITH_EPE
use ewaldpc_module, only: epe_read
#endif
use solv_cavity_module, only: solvation_read, check_read
use disp_rep_module, only: disp_rep_read, check_read_dr
use solv_charge_mixing_module, only: read_solv_mix
use potential_calc_module, only: poten_calc_read !!!!!!!!!!!!!!!AS
#ifdef WITH_MOLMECH
use qmmm_interface_module, only: qmmm_read !!!!!!!!!!!AS
#endif
#ifdef WITH_EFP
use efp_module, only: efp_read_input
#endif
#ifdef WITH_ERI4C
use eri4c_options, only: eri4c_read_input
#endif

implicit none
integer(kind=i4_kind), intent(in):: loop
! *** end of interface ***

integer(i4_kind) :: i
character(len=32) :: namelist_name
logical :: out_read, ua_read, basis_read, gr_read, ph_gr_read, pc_read, &
     symadapt_read, sym_read, mix_read, ferm_read, conv_read, mp_read, opt_read, &
     rec_read, xc_read, occ_read, mda_read, pop_read, dip_read, ef_read, &
     prop_read, orb_plot_read, frag_orb_read, resp_read, spin_orbit_read, solvat_read, &
     grid_ph_only, dis_re_read,epe_rd, pot_calc_read, gten_read, empiricalmethods_read
logical :: solv_mix_read, diis_read
#ifdef WITH_EFP
logical :: efp_read
#endif
#ifdef WITH_MOLMECH
logical :: qm_mm_read  !!!!!!!!!!!AS
#endif
logical :: eri4c_read
logical :: dft_plus_u_read




!------------ Executable code -----------------------------------

out_read = .false.
ua_read = .false.
basis_read = .false.
gr_read = .false.
ph_gr_read = .false.
sym_read = .false.
pc_read = .false.
symadapt_read = .false.
mix_read = .false.
diis_read = .false.
ferm_read = .false.
conv_read = .false.
mp_read = .false.
opt_read = .false.
rec_read = .false.
xc_read = .false.
dft_plus_u_read = .false.
occ_read = .false.
prop_read = .false.
frag_orb_read = .false.
orb_plot_read = .false.
mda_read = .false.
pop_read = .false.
dip_read = .false.
gten_read = .false.
ef_read = .false.
empiricalmethods_read = .false.
resp_read = .false.
spin_orbit_read = .false.
epe_rd = .false.  !!!!!!!!!!!!!!!AS
pot_calc_read = .false. !!!!!!!!!!!!!!AS
solvat_read = .false.
dis_re_read = .false.
solv_mix_read= .false.
#ifdef WITH_MOLMECH
qm_mm_read = .false. !!!!!!!!!!!AS
#endif
#ifdef WITH_EFP
efp_read=.false.
#endif
eri4c_read = .false.

!
! Opening input file. The  startup script was historically copying the
! input file  (not necessarily named  "input") from user  directory to
! scratch.  Now  we write all temp  files into the  same directory and
! also  avoid  copying  of  the  input completely.  We  expect  it  in
! input_dir/input_name as returned  by the function inpfile() instead.
! This  requires  modifications  of  the  startup  script  ./bin/runpg
! recorded together  with an  earlier change.  Make  sure you  are not
! using the old script if you get errors here.
!
call input_open(inpfile(input_name))

! read operations or tasks namelist that is required at beginning
call operations_read() ! to load at least the default values

! short-circuit reading input:
if(operations_optimizer_only) goto 999 ! close input file and exit!
if(operations_epe_lattice   ) goto 999 ! close input file and exit!
#ifdef WITH_MOLMECH
if(operations_mol_mech      ) goto 999 ! close input file and exit!
#endif

! preset some defaults
call options_set_defaults()

! read namelists until end of input is reached
do while ( input_which_namelist(namelist_name) )

   DPRINT 'namelist_name=', namelist_name
   select case(namelist_name)

   case ("output", "output_trace", "output_config", "output_unique_atoms",  &
        "output_scf", "output_sym", "output_integral", "output_timing", &
        "output_post_scf", "output_iounitadmin", "output_dipole", "output_solvation" )
      if (out_read) call input_error( &
           "READ_INPUT: output options must all be given together" )
      ! read output options
      call output_read() ! to load at least the default values
      out_read = .true.

   case("machineparameters")
      if (mp_read) call input_error("READ_INPUT: namelist machineparameters appears twice")
      ! read machine dependent parameters like veclen
      if (output_read_input) call write_to_output_units("read input: machineparameters")
      call machineparameters_read() ! to load at least the default values
      mp_read = .true.

   case("main_options")
      if (opt_read) call input_error("READ_INPUT: namelist main_options appears twice")
      if (output_read_input) call write_to_output_units("READ_INPUT: main_options")
      call options_read_input()
      opt_read = .true.

   case("recover_options")
      if (rec_read) call input_error("READ_INPUT: namelist recover_options appears twice")
      if (output_read_input) call write_to_output_units("READ_INPUT: recover_options")
      call options_read_input()
      rec_read = .true.

   case("symmetry_group")
      if ( sym_read ) call input_error( "READ_INPUT: symmetry information appears twice" )
      if (output_read_input) call write_to_output_units("READ_INPUT: symmetry_group")
      call symmetry_data_group_read()
      sym_read = .true.

   case("unique_atom_number")
      if ( ua_read ) call input_error( "READ_INPUT: namelist unique_atom_number appears twice" )

      ! read in number and names of unique atoms and first equal atom
      if (output_read_input) call write_to_output_units("read input: unique_atom_read")
      call unique_atom_read(loop) ! was unique_atom_unique_read()
      ua_read = .true.

   case("unique_atom_basisset")
      if ( basis_read ) call input_error( &
           "READ_INPUT: namelist unique_atom_basisset appears to often" )
      if ( .not. ua_read ) call input_error( &
           "READ_INPUT: namelist unique_atom_number must preceed namelist unique_atom_basisset" )
      ! read information about basis sets unique atoms
      if (output_read_input) call write_to_output_units("read input: unique_atom_basisset")
      call unique_atom_read_basis()
      basis_read = .true.

   case("pointcharge_number")
      if (pc_read) call input_error( "READ_INPUT: namelist pointcharge_number appears twice" )
      if (output_read_input) call write_to_output_units("read input: pointcharge")
      call pointcharge_read()
      pc_read = .true.

   case("mixing")
      if (mix_read) call input_error("READ_INPUT: namelist mixing appears twice" )
      if (output_read_input) call write_to_output_units("read input: mixing")
      call mixing_read()
      mix_read = .true.
   case("diis")
      if (diis_read) call input_error("READ_INPUT: namelist diis appears twice" )
      if (output_read_input) call write_to_output_units("read input: diis")
      call diis_read_input()
      diis_read = .true.
   case("fermi")
      if (ferm_read) call input_error( "READ_INPUT: namelist fermi appears twice" )
      if (output_read_input) call write_to_output_units("read input: fermi")
      call fermi_read()
      ferm_read = .true.

   case("convergence_list")
      if (conv_read) call input_error( "READ_INPUT: namelist convergence appears twice" )
      if(output_read_input) call write_to_output_units("read input: convergence_list")
      call convergence_read()
      conv_read = .true.

   case("grid","gridatom")
      if (gr_read) call input_error( "READ_INPUT: namelist grid appears twice" )
      if (ph_gr_read) call input_error( "READ_INPUT: namelist grid must appear&
                                        &  b e f o r e  namelist grid_ph")
      if ( .not. ua_read ) call input_error( &
           "READ_INPUT: namelist unique_atom_number must preceed namelist grid" )
      if (output_read_input) call write_to_output_units("read input: grid")
      call grid_read()
      gr_read = .true.

   case("grid_ph","gridatom_ph")
      if (ph_gr_read) call input_error( "READ_INPUT: namelist grid_ph appears twice" )
      if ( .not. ua_read ) call input_error( &
           "READ_INPUT: namelist unique_atom_number must preceed namelist grid_ph" )
      if (output_read_input) call write_to_output_units("read input: grid_ph")
      call grid_read_ph(grid_ph_only)
      ph_gr_read = .true.

   case("xc_control")
      if (xc_read) call input_error( "READ_INPUT: namelist xc_control appears twice" )
      if (output_read_input) call write_to_output_units("read input: xc_control")
      call xc_read_input()
      xc_read = .true.

#ifdef WITH_DFTPU
   case("dft_plus_u")
      if (dft_plus_u_read) call input_error( "READ_INPUT: namelist dft_plus_u appears twice" )
      if (output_read_input) call write_to_output_units("read input: dft_plus_u")
      call dft_plus_u_read_input()
      dft_plus_u_read = .true.

   case("dft_plus_u_fragment")
      call dft_plus_u_mo_read_input()
#endif

   case("spin_orbit")
      if (spin_orbit_read) call input_error( "READ_INPUT: namelist spin_orbit appears twice" )
      if (output_read_input) call write_to_output_units("read input: spin_orbit")
      call spin_orbit_read_input()
      spin_orbit_read = .true.

   case("occupation")
      if (occ_read) call input_error( "READ_INPUT: namelist occupation appears twice" )
      if ( .not. opt_read ) call input_error( &
           "READ_INPUT: namelist main_options must preceed namelist occupation" )
      if ( .not. ua_read ) call input_error( &
           "READ_INPUT: namelist unique_atom_number must preceed namelist occupatio" )
      if (output_read_input) call write_to_output_units("read input: occupation")
      call occupation_read()
      occ_read = .true.

   case("properties")
      if(prop_read) call input_error( "READ_INPUT: namelist properties occurs twice" )
      if(output_read_input) call write_to_output_units &
           ("read input: Reading properties information")
      call properties_read()
      prop_read= .true.

   case("frag_orb_analysis")
      if(frag_orb_read) call input_error( "READ_INPUT: namelist frag_orb_analysis occurs twice" )
      if(output_read_input) call write_to_output_units &
           ("read input: Reading frag_orb_analysis information")
      call frag_orb_analysis_read()
      frag_orb_read = .true.

   case("orbital_plot")
      if(orb_plot_read) call input_error( "READ_INPUT: namelist orb_plot occurs twice" )
      if(output_read_input) call write_to_output_units &
           ("read input: Reading orbital_plot information")
      call orbital_plot_read()
      orb_plot_read=.true.

   case("mda_options")
      if (mda_read) call input_error( "READ_INPUT: namelist mda_options appears twice" )
      if (output_read_input) call write_to_output_units("read input: mda_options")
      call mda_options_read()
      mda_read = .true.

   case("population")
      if (pop_read) call input_error( "READ_INPUT: namelist population appears twice" )
      if (output_read_input) call write_to_output_units("read input: population")
      call population_read()
      pop_read = .true.

   case("dipole_transitionmoments")
      if (dip_read) call input_error( "READ_INPUT: namelist dipole_transitionmoment appears twice" )
      if (output_read_input) call write_to_output_units("read input: dipole_transitionmoment")
      call dipole_read()
      dip_read = .true.

   case("magnetic_properties")
      if (gten_read) call input_error( "READ_INPUT: namelist gtensor_orbitals appears twice" )
      if (output_read_input) call write_to_output_units("read input: magnetic_properties")
#ifdef WITH_GTENSOR
      call magnetic_read()
      gten_read = .true.
#else
      ABORT('recompile -DWITH_GTENSOR')
#endif

   case("efield")
      if (ef_read) call input_error( "READ_INPUT: namelist efield appears twice" )
      if (output_read_input) call write_to_output_units("read input: efield")
      call efield_read()
      ef_read = .true.

   case("empirical_methods")
      ! TODO: remove in later versions
      if (empiricalmethods_read) call input_error( "READ_INPUT: dft-d obsolete"&
                             //"remove unsupported namelist empirical_methods" )
      if (output_read_input) call write_to_output_units("read input: empirical_methods")
!     call empirical_methods_read()
      empiricalmethods_read = .true.

   case("solvation")
      if (solvat_read) call input_error( "READ_INPUT: namelist solvation appears twice" )
      if (output_read_input) call write_to_output_units("read input: solvation")
      call solvation_read()
      solvat_read = .true.

   case("disper_repal")
      if (dis_re_read) call input_error( "READ_INPUT: namelist disper_repal appears twice" )
      if (output_read_input) call write_to_output_units("read input: disper_repal")
      call disp_rep_read()
      dis_re_read = .true.

   case("solv_charge_mixing")
      if (solv_mix_read) call input_error( "READ_INPUT: namelist solv_charge_mixing appears twice" )
      if (output_read_input) call write_to_output_units("read input: solv_charge_mixing")
      call read_solv_mix()
      solv_mix_read = .true.

   case("response_control")
      if(output_read_input) call write_to_output_units &
           ("read input: Reading response_control information")
#ifdef WITH_RESPONSE
      call response_read_input()
      resp_read = .true.
#else
      ABORT('recompile -DWITH_RESPONSE')
#endif

   case("unique_atom","unique_atom_basis","unique_atom_glob_con", &
        "unique_atom_symadapt","pointcharge","popcolum","fragment","orbital_list")
      call input_error("READ_INPUT: namelist "//trim(namelist_name)//" in wrong position")

   case("epe_data")
#ifdef WITH_EPE
      if(epe_rd) call input_error( "READ_INPUT: namelist epe_data appears twice" )
      if (output_read_input) call write_to_output_units("read input: epe_data")
      call epe_read(N_unique_atoms)
      epe_rd = .true.
#else
      ABORT('recompile -DWITH_EPE')
#endif

#ifdef WITH_EFP
   case("efp_data")
      if(efp_read) call input_error( "READ_INPUT: namelist efp_data appears twice" )
      if(.not. ua_read) call input_error( "READ_INPUT: namelist efp_data has to be read in after definition of QM atoms" )
      if (output_read_input) call write_to_output_units("read input: efp_data")
      call efp_read_input()
      efp_read=.true.
#endif

   case("potential")
      if(pot_calc_read) call input_error( "READ_INPUT: namelist potential appears twice" )
      if (output_read_input) call write_to_output_units("read input: poten_calc_read")
      call poten_calc_read()
      pot_calc_read = .true.

#ifdef WITH_MOLMECH
   case("qmmm")
      if(qm_mm_read) call input_error( "READ_INPUT: namelist QMMM appears twice" )
      if (output_read_input) call write_to_output_units("read input: qmmm_read")
      call qmmm_read()
      qm_mm_read = .true.
#endif

#ifdef WITH_ERI4C
   case("eri4c")
      if (eri4c_read) call input_error( "READ_INPUT: namelist eri4c appears twice" )
      if (basis_read) call input_error( "READ_INPUT: namelist eri4c "//        &
                                  "must preceed namelist unique_atom_basisset" )
      if (output_read_input) call write_to_output_units("read input: eri4c")
      call eri4c_read_input()
      eri4c_read = .true.
#endif

   case default
      call input_error( "READ_INPUT: unknown namelist "//trim(namelist_name) )

   end select

enddo

!
! This is supposed to detect lines that are not recognized as a namelist,
! and complain loud:
!
if ( .not. input_end_of_file() ) then
   call input_error("READ_INPUT: invalid input line")
endif

#ifdef WITH_GUILE
if (ua_read .and. .not. basis_read) then
   do i = 1, size(unique_atoms)
      call unique_atom_find_basis(unique_atoms(i))
   end do

   ! This no more means unique_atom_read_basis() has been called:
   basis_read = .true.
endif
#endif

#ifdef NEW_EPE
if (.not. operations_epe_lattice) then
#endif
! check if all obligatory information was given
if ( .not. ua_read ) call input_error("READ_INPUT: no atoms given")
if ( .not. basis_read .and. N_unique_atoms /= 0) call input_error("READ_INPUT: no basissets given")

! initialize the remaining modules
! the defaults will be set and nothing will really be read
if ( .not. out_read ) call output_read()
if ( .not. mp_read ) call machineparameters_read()
if ( (.not. opt_read) .and. (.not. rec_read) ) call options_read_input()
if ( .not. sym_read ) call symmetry_data_group_read()
if ( .not. pc_read ) call pointcharge_read()
if ( .not. mix_read ) call mixing_read()
if ( .not. diis_read ) call diis_read_input()
if ( .not. ferm_read ) call fermi_read()
if ( .not. conv_read ) call convergence_read()
if ( .not. gr_read ) call grid_read()
if ( .not. ph_gr_read ) call grid_copy_to_ph()      ! not grid_read_ph()
if ( .not. xc_read ) call xc_read_input()
if ( .not. spin_orbit_read ) call spin_orbit_read_input()
if ( .not. occ_read ) call occupation_read()
if ( .not. prop_read ) call properties_read()
if ( .not. frag_orb_read ) call frag_orb_analysis_read()
if ( .not. orb_plot_read ) call orbital_plot_read()
if ( .not. mda_read ) call mda_options_read()
if ( .not. pop_read ) call population_read()
if ( .not. dip_read ) call dipole_read()
if ( .not. solvat_read ) call check_read()
if ( .not. dis_re_read) call check_read_dr()
if ( .not. solv_mix_read) call read_solv_mix()
#ifdef WITH_GTENSOR
if ( .not. gten_read ) call magnetic_read()
#endif
if ( .not. ef_read ) call efield_read()
#ifdef WITH_RESPONSE
if ( .not. resp_read ) call response_read_input()
#endif
#ifdef WITH_EPE
if ( .not. epe_rd .and. N_unique_atoms > 0) call epe_read(N_unique_atoms)
#endif
if ( .not. pot_calc_read) call poten_calc_read()
#ifdef WITH_MOLMECH
if(.not. qm_mm_read) call qmmm_read()
#endif
#ifdef WITH_EFP
if (.not. efp_read) call efp_read_input()
#endif

! performing global post input operations
if (ph_gr_read .and. grid_ph_only) call grid_copy_to_ph(gridatom_only=.true.)
call count_electrons()
call check_occ_fermi() !TS: fermi_fix_up_and_down is replaced
                       !    by fixed_spin_diff
                       !    fermi_unpaired is replaced
                       !    by magnetic_moment
#ifdef NEW_EPE
endif
#endif

! write read data as desired
if (output_read_input) call write_to_output_units("read input: writing")
if (output_operations) then
   if (output_unit > 0) then
      write (output_unit, '()')
      call operations_write (output_unit, full=.true.)
   endif
endif

! consistency checks
if ( operations_scf .and. &
     ( options_integrals_on_file() .and. .not.  operations_integral ) ) &
     call error_handler( &
     "read_input: mainoptions integrals_on_file requires operations_integral")

!f ( operations_gradients .and. efield_applied() ) &
!    call error_handler( &
!    "read_input: sorry, gradients of external electrical field are not yet implemented")

if( options_spin_orbit .and. operations_gradients )then
   WARN('SO: gradients not implemented, ignoring')
   operations_gradients = .false.
endif

999 CONTINUE ! close input file and exit
    ! closing input file
    call input_close()

    if (output_read_input) call write_to_output_units("read input: done")

end subroutine read_input

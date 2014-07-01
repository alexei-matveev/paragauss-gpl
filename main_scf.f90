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
subroutine main_scf()
  !
  !  Purpose: Main Level of the SCF part. Executed on all workers.
  !
  !  Subroutine called by: main_master
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !===============================================================
  ! End of Public interface of module
  !===============================================================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Author: UB
  ! Date:   27/5/97
  ! Description:  Introduction of the "USE_MODEL_DENSITY" option
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description:   ...
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description:   pvm -> comm
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description:   ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module, only: i4_kind, r8_kind ! type specification parameters
  use interfaces, only: chargefit
  use paragauss, only: toggle_legacy_mode
  use prescf_module, only: prescf_init, prescf_finalize
  use time_module, only: start_timer, stop_timer ! timing routines
  use timer_module, only: timer_grid_setup, timer_scf, &
       timer_scf_cycle, timer_scf_ham, timer_scf_preparations, &
       timer_scf_reoc, timer_scf_xc, timer_grid_small_to_large, &
       timer_print_scfcycle, timer_print_scf
  use occupation_module, only: fixed_hole, n_occo, reoccup,&
       print_occ_num, occupation_2d_correct, &
       occupation_print_spectrum, occupation_get_holes
  use fermi_module, only: fermi_reoccup, fermi_level_broad, fermi_get_entropy, &
       fermi_write_input, fermi_read_scfcontrol
  use occupied_levels_module, only: sndrcv_eigvec_occ
  use options_module, only: xcmode_model_density, xcmode_extended_mda, &
       xcmode_numeric_exch, recover_fitcoeff, recover_scfstate, recover_ksmatrix, &
       xcmode_exchange_fit, recover_eigenvec, recover_fragment, options_save_ksmatrix, &
       options_save_scfstate, options_save_interval, options_save_eigenvec_all, &
       options_save_eigenvec, options_save_densmatrix, options_save_fitcoeff, &
       options_xcmode, options_recover, lvshift_mixing, options_spin_orbit, &
       recover_nothing, options_save_as_fragment, options_reset_scfcycle
  use output_module, only: output_data_read, output_data_saved, output_timing_scfloops, &
       output_scfloops, output_timing_scf, output_timing_detailedscf, &
       output_spectrum, output_reoccup, output_each_eigendata, output_main_scf, &
       output_eigendata, output_densmat, output_chargefit
  use iounitadmin_module, only: get_iounit, return_iounit, output_unit, stdout_unit, &
       write_to_trace_unit, write_to_output_units, openget_iounit, returnclose_iounit
  use grid_module, only: grid_main, grid_close
  use xc_cntrl, only: xc_is_on=>is_on, xc_ANY
  use xc_hamiltonian, only: xc_setup, xc_build, &
       xc_hamiltonian_store, xc_hamiltonian_recover, xc_close, s_average
  use xcfit_hamiltonian, only: xcfit_setup, xcfit_build, xcfit_close
  use xcmda_hamiltonian, only: xcmda_setup, xcmda_coeff_store, &
       xcmda_build, xcmda_coeff_recover, xcmda_close, &
       rho_exact, comp_exact
  use convergence_module, only: convergence_abort_calculation, convergence, &
       convergence_check_density_dev, convergence_max_iterations, &
       convergence_state_store, convergence_check_coulomb_dev, &
       convergence_check_coeff_dev, convergence_state_recover, &
       convergence_clear_buffers, convergence_setup, convergence_resize_buffers, &
       convergence_put_energy, convergence_put_density_dev, convergence_put_coeff_dev, &
       convergence_put_coulomb_dev, convergence_shutdown
  use eigen_data_module, only: print_eigendata, eigvec_write, eigen_data_dump
  use fit_coeff_module, only: print_coeff_charge, &
       fit_coeff_store, fit_coeff_recover, fit_coeff_send, fit_coeff_sndrcv, &
       fit_coeff_normalize, ICHFIT
#ifdef WITH_CORE_DENS
  use fit_coeff_module, only: write_coeff_core
#endif
  use density_data_module, only: print_densmat, gendensmat_occ, &
       save_densmat, densmat
  use overlap_module, only: overlap
  use filename_module, only: outfile, recfile
  use hamiltonian_module, only: ham_tot, hamiltonian_store, hamiltonian_recover
  use eigen_data_module, only: eigen_data_solve, eigen_hole_shutdown, &
       build_lvsft_hamiltonian, make_level_shift
  use diis_fock_module, only: diis_fock_matrix, fixed_fock_matrix, &
       diis_fock_module_close, diis_on
  use energy_calc_module, only: write_energies, get_energy
  use operations_module, only: operations_geo_opt, &
       operations_core_density, operations_solvation_effect, operations_qm_mm
  use symmetry_data_module, only: symmetry_data_n_spin
  use ham_calc_module, only: ham_calc_main
  use mixing_module, only: mixing_state_store, mixing_state_recover, &
       mixing_close, mixing_write_input, mixing_reset_buffers, &
       mixing_beta_ch, mixing_beta_sp, mixing_beta_xc, level_shift, &
       start_after_cycle
#ifdef WITH_SCFCONTROL
  use mixing_module, only: mixing_read_scfcontrol
#endif
  use population_module, only: population_mulliken, population_spor_mulliken
  use comm, only: comm_parallel, comm_same, comm_bcast, comm_rank
  use pairs_module, only: n_pairs
  use readwriteblocked_module, only: readwriteblocked_tapehandle, &
       readwriteblocked_stopread, readwriteblocked_startread, &
       readwriteblocked_startwrite, readwriteblocked_read, readwriteblocked_stopwrite, &
       readwriteblocked_write
  use solv_electrostat_module, only: charge_mix_wrapper, build_solv_ham, &
       dealloc_ham_solv, solv_energy_el, calc_q_e
  use solv_cavity_module, only: sol_start_cycle, n_Q_update
  use solv_charge_mixing_module, only: mix_charges, Qs_mix, Q_dealloc
  use induced_dipoles_module, only: calc_Pol_centers, free_field_arrays, n_update
  use calc_id_module, only: build_Pol_ham, calc_id_energy, dealloc_Pol_ham
#ifdef WITH_EFP
  use induced_dipoles_module, only: do_pol_pcm
  use efp_module, only: n_efp, qm_fixed, print_id
  use efp_polar_module, only: deallocate_Efield
  use efp_solv_module, only: allocate_V_and_Q_id, calc_V_and_Q_id, allocate_E_cav, do_solv
  use efp_solv_module, only: deallocate_E_cav, deallocate_V_id
#endif
  use s2_expect, only: s2_calc
  implicit none
  ! *** end of interface ***

  !------------ Declaration of local variables --------------------
  integer (i4_kind) :: xcmode, recover_mode, loop
  integer (i4_kind) :: first_loop, trace_width(14)
  real (r8_kind) :: tot_en, spin_diff, time_for_last_cycle, entropy
  real (r8_kind) :: coeff_dev, coulomb_dev, density_dev
  character (len=16) :: trace_format(15), trace_line*101

  integer, parameter :: NCHAR = 18 ! string length
  character (len=NCHAR), parameter :: fit_file = 'saved_fitcoeff.dat'
  character (len=NCHAR), parameter :: scf_file = 'saved_scfstate.dat'
  character (len=NCHAR), parameter :: ham_file = 'saved_ksmatrix.dat'
  character (len=NCHAR), parameter :: eig_file = 'saved_eigenvec.dat'
  character (len=NCHAR), parameter :: frg_file = 'saved_fragment.dat'

  ! One of the above, including .dat extension:
  character (len=NCHAR) :: rec_file

  logical :: scf_conv ! signifies SCF convergence, broadcasted to slaves
  logical :: check_density_dev
  logical :: store_now, recover, etot_recovered
  logical :: new_trace_header
  logical :: update_scfcontrol = .FALSE.
  integer (i4_kind), save :: scf_count = 0 ! counts calls to main_scf, i.e. geo_loops
!ifdef WITH_EFP
  logical :: do_update
!endif

  !
  ! As of  now main_scf()  runs on all  workers.  If you  are thinking
  ! about adding code  to be executed on all  workers, consider adding
  ! code/calls  outside of  the legacy  mode blocks  or to  other subs
  ! executed in parallel:
  !
  !   - prescf_init()
  !   - grid_main()
  !   - ham_calc_main()
  !   - prescf_finalize()
  !   - ...
  !
  call say ("start")

  ! We are called once again, remember that:
  scf_count = scf_count + 1

  ! Write to trace unit to show progress of calculation
  call write_to_trace_unit ("Entering SCF part")

  call start_timer (timer_scf)
  call start_timer (timer_scf_preparations)

  ! Check  status  of  recover   files  (before  entering  prescf  and
  ! xc_setup_main)
  call init_recover (recover, recover_mode, rec_file) ! all intent (out)

  ! Preparatory tasks for scf-cycle
  call say ("prescf_init")
  !
  ! This sub is executed on all workers, put the code that is common
  ! to all workers inside (or rather call from there):
  !
  call prescf_init()

  xcmode = options_xcmode()

  ! Prepare calculation of exchange-correlation hamiltonian:
  if (xc_is_on (xc_ANY)) then
     ! Set up grid
     call say ("grid setup")
     call start_timer (timer_grid_setup)
     call say ("call grid_main()")
     call grid_main (.false.)
     call say ("done grid_main()")
     call stop_timer (timer_grid_setup)

     ! Call xc setup:
     call say ("call xc_setup()")
     select case (xcmode)
     case (xcmode_model_density, xcmode_extended_mda)
        call xcmda_setup ()
     case (xcmode_numeric_exch)
        call xc_setup ()
     case (xcmode_exchange_fit)
        call xcfit_setup ()
     end select
     call say ("done xc_setup()")
  endif

  if (output_chargefit) then
     call say ("coeff_charge is written to file")
     call print_coeff_charge (0)
  endif

  call say ("other preparations")

  ! Reset  convergence test  variables and  allocated  the convergence
  ! buffers
  call convergence_setup()

#ifdef WITH_SCFCONTROL
  ! Prepare file with steering parameters read at each cycle
  if (scf_count .eq. 1) then
      call write_scfcontrol()
  else
      call read_scfcontrol()
  endif
#endif /* ifdef WITH_SCFCONTROL */

  ! It was previousely done every SCF iteration in
  ! convergence_read_scfcontrol() called from read/write_scfcontrol
  call convergence_resize_buffers()

  call stop_timer (timer_scf_preparations)
  call say ("preparations done")
  call write_to_trace_unit ("SCF preparations done")

  ! Here the SCF-Cycle starts
  call write_trace_header()

  time_for_last_cycle = 0.0_r8_kind

  ! FIXME: what about tot_en?
  loop = 0
  etot_recovered = .false.
  ASSERT(comm_same(recover))
  if (recover) then
     !
     ! Loads   loop,   etot_recovered   and  tot_en.    When   calling
     ! hamiltonian_recover() the  Fock matrix in  ham_tot is allocated
     ! and filled with data. Master  seems to send the orders and data
     ! to slaves  from there, the correspondig  receives are initiated
     ! from main_slave():
     !
     do while (toggle_legacy_mode ())
        call do_recover (recover_mode, n_pairs, loop, tot_en, eof=etot_recovered)
     enddo
     call comm_bcast (loop)
     call comm_bcast (tot_en)
     call comm_bcast (etot_recovered)

     ! loop is incremented at the top of the scf_cycle:
     loop = loop - 1
  endif

  ! This is compared against base-1 loop:
  first_loop = loop + 1

  ASSERT(comm_same(fixed_hole))
  if (fixed_hole) then
     call say ("set up hole information in eigen_data_module")
     call occupation_get_holes()
  endif

  scf_cycle: do
     loop = loop + 1
     call start_timer (timer_scf_cycle)

     ! Saving state every  iteration may require much IO,  do it maybe
     ! every few iterations:
     store_now = store_now_p (loop, options_save_interval())

     if (output_timing_scfloops .or. output_scfloops) then
        call write_to_output_units (" ")
        call write_to_output_units ("####### Loop: ", loop)
        call write_to_output_units (" ")
     endif

     ! Check recover options and fork accordingly
     ASSERT(comm_same(loop))
     ASSERT(comm_same(first_loop))
     if (loop == first_loop .and. recover) then
        ASSERT(comm_same(recover_mode))
        select case (recover_mode)
        case (recover_fitcoeff); goto 1000
        case (recover_scfstate); goto 1000
        case (recover_ksmatrix); goto 2000
        case (recover_eigenvec); goto 3000
        case (recover_fragment); goto 3000
        end select
     endif

     ! Calculation of exchange-correlation hamiltonian:
     if (xc_is_on (xc_ANY)) then
        ASSERT(comm_same(xcmode))
        select case (xcmode)

        case (xcmode_model_density, xcmode_extended_mda)
           ABORT("verify SPMD")
           call say ("fit_coeff_sndrcv")
           call fit_coeff_sndrcv (ICHFIT)
           if ((rho_exact .or. comp_exact) .and. loop /= first_loop) then
              call sndrcv_eigvec_occ()
           endif
           call start_timer (timer_scf_xc)
           call say ("xcmda_build_ham")
           call xcmda_build ((loop == first_loop))
           call stop_timer (timer_scf_xc)

        case (xcmode_numeric_exch)
           if (loop /= first_loop) then
              call start_timer (timer_scf_xc)
              call say ("xc_build")
              ! Call with relative cycle number:
              call xc_build (loop - first_loop + 1)
              call stop_timer (timer_scf_xc)
           endif

        case (xcmode_exchange_fit)
           if (loop /= first_loop) then
              call start_timer (timer_scf_xc)
              call say ("xcfit_build")
              call xcfit_build (loop)
              call stop_timer (timer_scf_xc)
           endif
        end select
     endif

     ! The code inside  this DO block is executed  on master only. The
     ! slaves  are spinning  in  main_slave() and  waiting for  orders
     ! during that time:
     do while (toggle_legacy_mode ())

     if (calc_Pol_centers() .and. loop >= first_loop+1) then
        if (loop == first_loop+1) then
           do_update = .true.
        else
           do_update = (mod (loop - (first_loop + 1), n_update) == 0)
        end if
#ifdef WITH_EFP
        if (operations_solvation_effect) then
           if (do_pol_pcm .and. n_efp > 0) then
              call allocate_V_and_Q_id()
              call allocate_E_cav()
           end if
           do_solv = (loop > first_loop+sol_start_cycle)
        end if
        call build_Pol_ham (print_id, do_update)
#else
        call build_Pol_ham (do_update)
#endif
     end if

     ! Calculating solvation hamiltonian
     if (operations_solvation_effect .and. &
          loop >= first_loop+sol_start_cycle) then  !!!!!!!!!!!!!
        if (mod (loop - (first_loop + sol_start_cycle), n_Q_update) == 0) then
           call calc_Q_e()
        end if
#ifdef WITH_EFP
        if (do_pol_pcm .and. n_efp > 0) then
           call calc_V_and_Q_id()
        end if
#endif
        call charge_mix_wrapper (loop, first_loop + sol_start_cycle)
        call build_solv_ham()
     endif                                          !!!!!!!!!!!!!

     enddo                      ! while (toggle_legacy_mode())
     ! Parallel context from here on ...

     if (output_densmat) then
        call say ("densmat is written to file")
        call print_densmat (loop)
     endif

     ! Entry point for "read_fitcoeff" and "read_scfstate":
1000 continue

     ! Save fit coefficients and current SCF state if required
     if (options_save_fitcoeff()) then
        call do_fitcoeff_store (store_now) ! reads loop
     endif
     ! Save current SCF state if required
     if (options_save_scfstate()) then
        call do_scfstate_store (store_now) ! reads loop
     endif

     ! Calculation  of  other  (not  exchange  correlation)  parts  of
     ! hamiltonian and summation  to total hamiltonian and calculation
     ! of energies
     call say ("ham_calc_main")
     call start_timer (timer_scf_ham)

     !
     ! Build Coulomb  and Fock  matrices.  This also  call reset_ham()
     ! which allocates  and initializes  Fock matrix in  ham_tot. Some
     ! recover options jump over this call.
     !
     call ham_calc_main (loop - first_loop + 1) ! relative cycle number

     call stop_timer (timer_scf_ham)

     do while (toggle_legacy_mode ())

     if (operations_solvation_effect .and. &
          loop >= first_loop + sol_start_cycle) then
        call dealloc_ham_solv()
        call solv_energy_el ()
     endif

     if (calc_Pol_centers() .and. loop >= first_loop+1) then
        call dealloc_Pol_ham()
        call calc_id_energy()
     end if

     enddo                      ! while (toggle_legacy_mode())
     ! Parallel context from here on ...

     ! Write calculated  energies. Slaves  could always print  that to
     ! stdout, however the values are partial there.
     call say ("write_energies")
     if (output_scfloops .or. output_main_scf) then
        call write_energies (output_unit)
        if (comm_rank() == 0) then
           call write_energies (stdout_unit)
        endif
     endif
     if (etot_recovered) then
        if (output_scfloops .or. output_main_scf) then
500        format(2X, 'e_sum  =  ', F25.12, '  (recovered)')
           write (output_unit, 500) tot_en
           if (comm_rank() == 0) then
              write (stdout_unit, 500) tot_en
           endif
        endif
     else
        call get_energy (tot=tot_en)
     endif
     ! And pass energy to the convergence check module
     call convergence_put_energy (tot_en)

     if (store_now) then
        ! Write tot_en to a file:
        call store_total_energy()
     endif

     ! Entry point for "read_ksmatrix":
2000 continue

     ! Save Kohn-Sham matrix if required. NOOP on slaves:
     if (options_save_ksmatrix()) then
        call do_ksmatrix_store (store_now) ! reads loop
     endif

#ifdef WITH_SCFCONTROL
#warning "WITH_SCFCONTROL=1 is not useful for batch runs"
     ! Read file with steering  parameters first and modify the length
     ! of  convergence buffers  (if required)  keeping a  copy  of the
     ! original   data.   Caution:   reloads  options_save_interval(),
     ! options_save_fitcoeff().  If  you  want  every worker  to  know
     ! updated parameters, you either  have to let everyone read them,
     ! or broadcast them.
     call read_scfcontrol()
#endif /* ifdef WITH_SCFCONTROL */

     ! Check convergence and exit loop in case of convergence
     call say ("convergence is checked")

     ! Check SCF convergence.  FIXME: maybe make sure all workers have
     ! the data to determine convergence on their own?
     scf_conv = (convergence() .or. convergence_max_iterations (loop)) &
          .and. loop /= first_loop
     call comm_bcast (scf_conv)

     if (scf_conv) then
        call stop_timer (timer_scf_cycle)
        call timer_grid_small_to_large()

        ! Does nothing when the slaves have no output:
        if (output_timing_scfloops) then
           call timer_print_scfcycle()
        endif

        ! Uses  global variables loop,  tot_en and  old values  of the
        ! save_data options. NOOP on slaves.
        call do_final_store (n_vir=n_pairs)

        call say ("exit scf_cycle")

        ! Let  everyone exit  the SCF  block.  This  is the  only exit
        ! point:
        exit scf_cycle
     endif

     ! Now  update the save_data  options, deallocate  the temporarily
     ! kept  convergence  buffers (as  far  as  necessary), reset  the
     ! mixing buffers  (if required)  and update the  scf_control file
     ! accordingly
     store_now = store_now_p (loop, options_save_interval())
     call convergence_clear_buffers()

     call mixing_reset_buffers (loop, update_scfcontrol)
     ! Note that update_scfcontrol is also set by read_scfcontrol

#ifdef WITH_SCFCONTROL
     if (update_scfcontrol) call write_scfcontrol()
#endif /* ifdef WITH_SCFCONTROL */

     !
     ! In diis_fock_matrix() the total hamiltonian could be updated in
     ! a DIIS routine to get a convergence acceleration.  Only ham_tot
     ! will be  changed by  the routine.  It  could be switched  on by
     ! diis_on = true in the  DIIS input namelist.  Parameters need to
     ! be adapted; their best choice depend on other mixing strategies
     ! like Perturbation Theory and Fermi broadening.
     !
     if (diis_on) then
       ASSERT (.not.options_spin_orbit)
       ASSERT (allocated(overlap))
       ASSERT (allocated(ham_tot))
       ASSERT (allocated(densmat))
       call diis_fock_matrix (ham_tot, overlap, densmat, loop)

       call fixed_fock_matrix (ham_tot, loop)
     endif

     ! Eigensolver to determine orbital coefficients
     call say ("eigen_data_solve")

     !
     ! Perform Hamiltonian diagonalization, see eigen_data_module.
     ! FIXME: pass arguments explicitly.
     !
     call eigen_data_solve()

     ! Write eigendata to temporary file "eigenmist" for debugging
     if (output_each_eigendata) then
        call say ("print_eigendata")
        call print_eigendata()
     endif

     ! Reoccupation of orbitals
     call start_timer (timer_scf_reoc)

     if (fermi_level_broad) then
        call say ("fermi_reoccup")
        ! FIXME:  The input  differes  on master  and  slaves so  that
        ! Newton iteration fails:
        if (comm_rank() == 0) then
           call fermi_reoccup()
        endif
     else
        call say ("reoccup")
        call reoccup()
     endif

     call say ("occupation_2d_correct")
     if (.not. options_spin_orbit) then
        ! Does  no  communication.  Rewrites  occ_num(),  n_occo()  in
        ! place:
        call occupation_2d_correct()
     endif

     call stop_timer (timer_scf_reoc)

     ! Call on  master after eigen_data_solve()  and reocupp(). FIXME:
     ! The procedure  does no  communication, but fiddles  with global
     ! variables that may not be available on slaves. Verify this next
     ! time the ABORT() fires.
     if (lvshift_mixing) then
        ABORT("verify SPMD")
        call build_lvsft_hamiltonian &
             (n_occo, level_shift, loop .gt. START_AFTER_CYCLE)
     endif

     if (output_spectrum) then
        call say ("print out spectrum")
        call occupation_print_spectrum()
     endif

     if (output_reoccup .or. output_each_eigendata) then
        call say ("print_occ_num")
        call print_occ_num (loop)
     endif

     ! Entry point for "read_eigenvec":
3000 continue

     ! Save eigenstates and occupation if required. NOOP on slaves.
     if (options_save_eigenvec()) then
        call do_eigenvec_store (store_now, n_vir=n_pairs) ! reads loop
     endif

     ! Setup  of  density  matrix  (on  each  slave)  and  fit  charge
     ! fitfunction coefficients
     call say ("send_eigvec_occ")

     call sndrcv_eigvec_occ()

     call say ("gendensmat_occ")

     ! The function here seems to  just return the input parameter and
     ! thus should give the same value on all workers, modulo bugs, as
     ! always:
     check_density_dev = convergence_check_density_dev() &
          .and. .not. (recover .and. loop == first_loop)

     ASSERT(comm_same(check_density_dev))
     if (check_density_dev) then
        call gendensmat_occ (density_dev)
     else
        call gendensmat_occ()
        if (recover .and. loop == first_loop) then
           density_dev = 0.0_r8_kind
        else
           density_dev = huge (0.0_r8_kind)
        endif
     endif
     call convergence_put_density_dev (density_dev)

     !
     ! Moved from gendensmat_occ() here:
     !
     if (output_densmat) then
        call print_densmat (loop)
     endif

     call say ("chargefit")
     !
     ! One of  the first  thigs it  does is sending  a message  to the
     ! slaves  telling  them  to  also  call  chargefit()  with  dummy
     ! arguments that are assumed to be unused on slaves:
     !
     call chargefit (loop, coeff_dev, coulomb_dev)

     call convergence_put_coeff_dev (coeff_dev)
     call convergence_put_coulomb_dev (coulomb_dev)

     if (output_chargefit) then
        call say ("coeff_charge is written to file")
        call print_coeff_charge (loop)
     endif

     ! Calculate and print timing summary for SCF cycle
     call say ("timing")
     call stop_timer (timer_scf_cycle)
     time_for_last_cycle = timer_scf_cycle%real_timediff
     call timer_grid_small_to_large()
     if (output_timing_scfloops) call timer_print_scfcycle()

     ! Write to trace unit to show progress of calculation
     if (new_trace_header) then
        call write_trace_header()
     endif
     call write_trace_item()

     !
     ! In the next loop the energy is computed, not recovered:
     !
     etot_recovered = .false.

     call say ("cycle done")

  enddo scf_cycle

  ! This  checks convergence  of  energy AND  charge coefficients  AND
  ! density  matrix AND Coulomb  (?).  FIXME:  This function  seems to
  ! complete without crash but returns junk on slave workers:
  scf_conv = convergence()
  call comm_bcast (scf_conv)

  if (scf_conv) then
     call write_to_output_units ("")
     call write_to_output_units ("")
     call write_to_output_units ("####### convergence was reached in cycle ", inte=loop)
     call write_to_output_units ("")
     call write_to_output_units ("")
     call write_to_trace_unit ("MAIN_SCF: convergence was reached in cycle ", inte=loop)
  else
     call write_to_output_units ("")
     call write_to_output_units ("####### aborting at maximal number of cycles: ", inte=loop)
     call write_to_output_units ("")
     call write_to_output_units ("")
     call write_to_trace_unit ("MAIN_SCF: aborting at maximal number of cycles: ", inte=loop)
  endif

  ! Compute expectation  values of S2. Does some  output. FIXME: looks
  ! like we have a serial O(N^3) step here:
  ASSERT(allocated(n_occo))
  if (size (n_occo, 1) == 2) then ! unrestricted
     call s2_calc() ! from s2_expect.f90
  endif

#ifdef WITH_CORE_DENS
  if (operations_core_density) then
     ABORT("verify SPMD")
     call write_coeff_core (loop)
  end if
#endif

  !
  ! The diis_fock_matrix()  routine will not  be used anymore  in this
  ! SCF,    so   we    deallocate   the    storage    matrices.    The
  ! diis_charge_coeff() routine  is treated  the same way.   Make sure
  ! slaves do  no harm when executing  it!  FIXME: maybe  move it into
  ! prescf_finalize()?
  !
  call diis_fock_module_close()


  ! Clean up the XC part:
  if (xc_is_on (xc_ANY)) then
     call say ("call grid_close()")
     call grid_close (.true.)
     call say ("done grid_close()")

     call say ("call xc_close()")
     select case (xcmode)
     case (xcmode_model_density, xcmode_extended_mda)
        call xcmda_close ()
     case (xcmode_numeric_exch)
        call xc_close ()
     case (xcmode_exchange_fit)
        call xcfit_close ()
     end select
     call say ("done xc_close()")
  endif

  ! Write final calculated energies
  call say ("write final calculated energies")
  if (operations_solvation_effect) then
     do while (toggle_legacy_mode ())
        call calc_Q_e ()
        call solv_energy_el ()
        call Q_dealloc ()
     enddo
  endif

  if (calc_Pol_centers()) then
     call calc_id_energy()
#ifdef WITH_EFP
     if (.not. qm_fixed) call free_field_arrays()
#else
     call free_field_arrays()
#endif
#ifdef WITH_EFP
     if (n_efp > 1 .and. .not. qm_fixed) then
        call deallocate_Efield()
     end if
     if (n_efp > 0 .and. operations_solvation_effect .and. do_pol_pcm) then
        call deallocate_V_id()
        call deallocate_E_cav()
     end if
#endif
  end if

  call write_energies (output_unit)
  if (comm_rank() == 0) then
     if (output_main_scf) call write_energies (stdout_unit)
  endif
  call get_energy (tot=tot_en)

  ! Write to trace unit to show progress of calculation
  call write_to_trace_unit ("MAIN_SCF: final total energy: ", real=tot_en)

  call say ("population analysis in progress ...")
  if (options_spin_orbit) then
     call population_spor_mulliken()
  else
     call population_mulliken()
  endif

  ! Write eigendata to temporary file "eigvec.dat" for debugging
  if (output_eigendata) then
     call say ("print_eigendata")
     call print_eigendata()

     !
     ! Dumps all of hamiltonian, overlap, eigenvectors and eigenvalues
     ! into the named file in the output directory:
     !
     call eigen_data_dump ("eigendata.txt")

     if (options_spin_orbit) then
        DPRINT 'main_scf: call eigvec_write()'
        call eigvec_write()
        DPRINT 'main_scf: .'
     endif
  endif

  if (lvshift_mixing) then
     make_level_shift=.false.
  endif

  if (options_save_densmatrix()) call save_densmat !!!!!!!!!!!!AS

  call say ("call prescf_finalize()")
  call prescf_finalize()
  call say ("done prescf_finalize()")

  call say ("call mixing_close()")
  call mixing_close (loop, xcmode == xcmode_exchange_fit)
  call say ("done mixing_close()")


  if (fixed_hole) then
     call say ("shutdown hole information in eigen_data_module")
     call eigen_hole_shutdown()
  endif

  ! Print occupation of orbitals:
  call say ("call occupation_print_spectrum()")
  call occupation_print_spectrum ()
  call say ("done occupation_print_spectrum()")

  call stop_timer (timer_scf)

  ! Print timing summary for scf part:
  if (output_timing_scf .or.  output_timing_detailedscf) then
     call say ("timer_print_scf()")
     call timer_print_scf()
  endif

  if (convergence_abort_calculation()) &
       call error_handler (&
       "MAIN_SCF: aborting calculation because no convergence was reached")

  ! Deallocates the convergence test buffers:
  call say ("convergence_shutdown")
  call convergence_shutdown()

  call say ("end")

contains

  function store_now_p (loop, save_interval) result (yes)
    implicit none
    integer, intent (in) :: loop, save_interval
    logical :: yes
    ! *** end of interface ***

    if (save_interval /= 0) then
       yes = mod (loop, save_interval) == 0
    else
       yes = .false.
    endif
  end function store_now_p

#ifdef WITH_SCFCONTROL
  subroutine write_scfcontrol()
    !
    ! Purpose: write  variables that can  be changed at each  cycle to
    ! file.   This  function  when  executed  in  parallel,  does  not
    ! necessarily try writing to  the same file. The outfile() returns
    ! a path  to a  scratch directory for  slaves, different  for each
    ! worker.
    !
    ! To avoid  warning messages due to spurious  modifications of the
    ! scf control  data caused by the finite  accuracy achievable with
    ! formatted  write and  re-read,  the scf  control parametesr  are
    ! re-read immediatelly after writing (UB, 8/97)
    !
    use iounitadmin_module, only: openget_iounit, returnclose_iounit
    implicit none
    !** End of interface ***************************************

    integer :: unit

    ! not really an output file, but happens to reside in output directory:
    unit = openget_iounit (trim (outfile ("scfcontrol")), status='replace')
    call convergence_write (unit, full=.true.)
    call options_write_input (unit, scfcontrol=.true.)
    call mixing_write_input (unit, scfcontrol=.true.)
    call fermi_write_input (unit, scfcontrol=.true.)
    call returnclose_iounit (unit, status='keep')

    call read_scfcontrol (warning=.false.)
    update_scfcontrol = .false.
  end subroutine write_scfcontrol

  subroutine read_scfcontrol (warning)
    !
    ! Purpose: read variables  that can be changed at  each cycle from
    ! file. FIXME:  as is it looks  like each worker  would be reading
    ! different files as  outfile() points to scratch for  all but the
    ! rank-0 worker.
    !
    use input_module, only: input_open, input_close
    logical, optional, intent (in) :: warning
    !** End of interface ***************************************

    ! not really an output file, but happens to reside in output directory:
    call input_open (trim (outfile ("scfcontrol")))

    call convergence_read_scfcontrol (warning, new_trace_header)
    call options_read_scfcontrol (warning)
    call mixing_read_scfcontrol (warning)
    call fermi_read_scfcontrol (warning, new_trace_header, update_scfcontrol)

    call input_close()
  end subroutine read_scfcontrol

#endif /* ifdef WITH_SCFCONTROL */

  subroutine init_recover (recover, recover_mode, rec_file)
    !
    ! Check  recover  files   and  set  recover  options  accordingly.
    ! Occasionally does a  minumum of IO by INQUIREing  existance of a
    ! file.
    !
    ! INTENT (IN)
    ! character (len=*) :: fit_file, scf_file, ham_file, eig_file
    !
    implicit none
    logical, intent (out) :: recover
    integer, intent (out) :: recover_mode
    character (len=NCHAR), intent (out) :: rec_file
    !** End of interface ***************************************

    character (len=13) :: rec_option

    recover      = .false.
    recover_mode = options_recover()
    select case (recover_mode)
    case (recover_fitcoeff)
       recover    = .true.
       rec_file   = fit_file
       rec_option = 'read_fitcoeff'
    case (recover_scfstate)
       recover    = .true.
       rec_file   = scf_file
       rec_option = 'read_scfstate'
    case (recover_ksmatrix)
       recover    = .true.
       rec_file   = ham_file
       rec_option = 'read_ksmatrix'
    case (recover_eigenvec)
       recover    = .true.
       rec_file   = eig_file
       rec_option = 'read_eigenvec'
    case (recover_fragment)
       recover    = .true.
       rec_file   = frg_file
       rec_option = 'read_fragment'
    end select
    if (recover) then
       inquire (file=recfile (rec_file), exist=recover)
       if (.not. recover) then
          if (operations_geo_opt .or. operations_qm_mm) then
             recover_mode = recover_nothing
             call say ('Warning: No Recover File "' // rec_file // '"')
             call say ('         "'//trim (rec_option)//'" temporarily turned off !')
          else
             call error_handler &
                  ('Recover File "' // rec_file // '" not found !')
          endif
       endif
    endif
  end subroutine init_recover


  subroutine do_recover (recover_mode, n_vir, loop, tot_en, eof)
    !
    ! Purpose: performs  the recovering at the beginning  of the first
    ! SCF cycle and resets the SCF loop counter accordingly.
    !
    ! INTENT (IN)
    ! character (len=*) :: rec_file
    ! INTENT (OUT)
    ! integer (i4_kind) :: loop
    ! real   (r8_kind) :: tot_en
    ! logical          :: eof (== etot_recovered)
    use occupation_module, only: eigenstates_recover
    implicit none
    integer (i4_kind), intent (in) :: recover_mode
    integer (i4_kind), intent (in) :: n_vir
    integer (i4_kind), intent (out) :: loop
    real (r8_kind), intent (out) :: tot_en
    logical, intent (out) :: eof
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th
    real (r8_kind)                :: energy(1)

    ! recovery may  fail, to indicate  that the values of  output args
    ! are somehwat arbitrary this will be set to true:
    eof = .false.

    if (recover_mode == recover_fragment) then
       call readwriteblocked_startread (recfile (frg_file), th, variable_length=.true.)
       call eigenstates_recover (n_vir, th)

       ! no data to recover tot_en from:
       eof = .true.

       !
       ! Close the file, set tot_en  to dummy value and return. FIXME:
       ! what about setting loop in this case?
       !
       goto 999 ! clean up and return
    endif

    call readwriteblocked_startread (recfile (rec_file), th, variable_length=.true.)
    call fit_coeff_recover (th)
    call fit_coeff_normalize (spin_coeff=options_reset_scfcycle())

    call convergence_state_recover (th)
    call mixing_state_recover (loop, th)

    select case (recover_mode)

    case (recover_fitcoeff)
       if (xcmode == xcmode_model_density .or. &
            xcmode == xcmode_extended_mda) then
        if (comm_parallel()) then
          call say ("fit_coeff_send during recover")
          call fit_coeff_send()
        endif
          call say ("XC(MDA) coefficients re-computed during recover")
          call start_timer (timer_scf_xc)
          ABORT("SPMD?")
          call xcmda_build (.false.)
          call stop_timer (timer_scf_xc)
       endif

       ! tot_en will  be copied from  energy(1), but only if  not eof,
       ! otherwise a zero is returned:
       call readwriteblocked_read (energy, th, eof)

    case (recover_scfstate)
       if (xcmode == xcmode_model_density .or. &
            xcmode == xcmode_extended_mda) then
        if (comm_parallel()) then
          call say ("fit_coeff_send during recover")
          call fit_coeff_send()
        endif
       endif
       call xcmda_coeff_recover (th)
       call xc_hamiltonian_recover (th)

       ! tot_en will  be copied from  energy(1), but only if  not eof,
       ! otherwise a zero is returned:
       call readwriteblocked_read (energy, th, eof)

    case (recover_ksmatrix)
       call hamiltonian_recover (th)
       ! no data to recover tot_en from:
       eof = .true.

    case (recover_eigenvec)
       call eigenstates_recover (n_vir, th)
       ! no data to recover tot_en from:
       eof = .true.
    end select

999 CONTINUE
    call readwriteblocked_stopread (th, status="KEEP")

    if (eof) then
       tot_en = 0.0_r8_kind ! dummy value for trace printing
    else
       tot_en = energy(1)
       if (output_data_read) then
          write (output_unit, '(/ a     )') 'Recovered total energy :'
          write (output_unit, '(  a     )') 'tot_en'
          write (output_unit, '(4es20.13)') energy(1)
       endif
    endif
  end subroutine do_recover


  subroutine do_fitcoeff_store (store_now, mode, tot_en)
    !
    ! Purpose:   saves  the   information  for   the  "saved_fitcoeff"
    ! file. NOOP on slaves.
    !
    ! INTENT (IN)
    ! integer (i4_kind) :: loop
    ! character (len=*) :: fit_file
    use comm, only: comm_rank
    implicit none
    logical, intent (in) :: store_now
    integer (i4_kind), intent (in), optional :: mode
    real (r8_kind), intent (in), optional :: tot_en
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th

    if (comm_rank() /= 0) return

    if (store_now) then
       call readwriteblocked_startwrite (recfile (fit_file), th, variable_length=.true.)
       call fit_coeff_store (th)
       ! Present when called from do_final_store():
       call convergence_state_store (th, mode=mode)
       call mixing_state_store (loop, th)
       ! Present when called from do_final_store():
       if (present (tot_en)) then
          call readwriteblocked_write ([tot_en], th)
       endif
       call readwriteblocked_stopwrite (th)
    else
       ! Temporily saves  all data which might be  modified before the
       ! end of the current SCF cycle is reached
       call convergence_state_store (mode=recover_fitcoeff)
    endif
  end subroutine do_fitcoeff_store


  subroutine do_scfstate_store (store_now, mode, tot_en)
    !
    ! Purpose:   saves  the   information  for   the  "saved_scfstate"
    ! file. NOOP on slaves.
    !
    ! INTENT (IN)
    ! integer (i4_kind) :: loop
    ! character (len=*) :: fit_file
    use comm, only: comm_rank
    implicit none
    logical, intent (in) :: store_now
    integer (i4_kind), intent (in), optional :: mode
    real (r8_kind), intent (in), optional :: tot_en
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th

    if (comm_rank() /= 0) return

    if (store_now) then
       call readwriteblocked_startwrite (recfile (scf_file), th, variable_length=.true.)
       call fit_coeff_store (th)
       ! Present when called from do_final_store():
       call convergence_state_store (th, mode=mode)
       call mixing_state_store (loop, th)
       call xcmda_coeff_store (th)
       call xc_hamiltonian_store (th)
       ! Present when called from do_final_store():
       if (present (tot_en)) then
          call readwriteblocked_write ([tot_en], th)
       endif
       call readwriteblocked_stopwrite (th)
    else
       ! Temporily saves  all data which might be  modified before the
       ! end of the current SCF cycle is reached
       call convergence_state_store (mode=recover_scfstate)
    endif
  end subroutine do_scfstate_store


  subroutine do_ksmatrix_store (store_now, mode)
    !
    ! Purpose:   saves  the   information  for   the  "saved_ksmatrix"
    ! file. Does nothing on the slave workers.
    !
    ! INTENT (IN)
    ! integer (i4_kind) :: loop
    ! character (len=*) :: fit_file
    use comm, only: comm_rank
    implicit none
    logical, intent (in) :: store_now
    integer (i4_kind), intent (in), optional :: mode
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th

    if (comm_rank() /= 0) return

    if (store_now) then
       call readwriteblocked_startwrite (recfile (ham_file), th, variable_length=.true.)
       call fit_coeff_store (th)
       ! Present when called from do_final_store():
       call convergence_state_store (th, mode=mode)
       call mixing_state_store (loop, th)
       call hamiltonian_store (th)
       call readwriteblocked_stopwrite (th)
    endif
  end subroutine do_ksmatrix_store

  subroutine do_eigenvec_store (store_now, n_vir, mode)
    !
    ! Purpose:   saves  the   information  for   the  "saved_eigenvec"
    ! file. Does nothing on slaves.
    !
    ! INTENT (IN)
    ! integer (i4_kind) :: loop
    ! character (len=*) :: fit_file
    use comm, only: comm_rank
    use occupation_module, only: eigenstates_store
    implicit none
    logical, intent (in) :: store_now
    integer (i4_kind), intent (in) :: n_vir
    integer (i4_kind), intent (in), optional :: mode
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th

    if (comm_rank() /= 0) return

    ! Mode present when called from do_final_store().
    if (store_now) then
       call readwriteblocked_startwrite (recfile (eig_file), th, variable_length=.true.)
       call fit_coeff_store (th, mode=mode)
       call convergence_state_store (th, mode=mode)
       call mixing_state_store (loop, th, mode=mode)
       if (options_save_eigenvec_all()) then
          call eigenstates_store (th=th, mode=mode)
       else
          call eigenstates_store (n_vir, th=th, mode=mode)
       end if
       call readwriteblocked_stopwrite (th)
    else
       ! Temporily saves  all data which might be  modified before the
       ! end of the current SCF cycle is reached
       call fit_coeff_store (mode=recover_eigenvec)
       call convergence_state_store (mode=recover_eigenvec)
       call mixing_state_store (loop, mode=recover_eigenvec)
       if (options_save_eigenvec_all()) then
          call eigenstates_store (mode=recover_eigenvec)
       else
          call eigenstates_store (n_vir, mode=recover_eigenvec)
       end if
    endif
  end subroutine do_eigenvec_store


  subroutine store_total_energy
    !
    ! Purpose: adds the current total energy to the saved_fitcoeff.dat
    ! and save_scfstate.dat file.
    !
    ! INTENT (IN)
    ! real   (r8_kind) :: tot_en
    ! character (len=*) :: fit_file, scf_file
    use comm, only: comm_rank
    implicit none
    !** End of interface ***************************************

    type (readwriteblocked_tapehandle) :: th

    if (comm_rank() /= 0) return

    ! FIXME: this might be a dead code, is it ever used?
    if (options_save_fitcoeff()) then
       call readwriteblocked_startwrite (recfile (fit_file), th, &
            variable_length=.true., append=.true.)
       call readwriteblocked_write ([tot_en], th)
       call readwriteblocked_stopwrite (th)
       if (output_data_saved) then
          write (output_unit, '(/ a     )') 'Stored total energy :'
          write (output_unit, '(  a     )') 'tot_en'
          write (output_unit, '(4es20.13)') tot_en
       endif
    endif

    ! This one is used regualrly:
    if (options_save_scfstate()) then
       call readwriteblocked_startwrite (recfile (scf_file), th, &
            variable_length=.true., append=.true.)
       call readwriteblocked_write ([tot_en], th)
       call readwriteblocked_stopwrite (th)
       if (output_data_saved) then
          write (output_unit, '(/ a     )') 'Stored total energy :'
          write (output_unit, '(  a     )') 'tot_en'
          write (output_unit, '(4es20.13)') tot_en
       endif
    endif
  end subroutine store_total_energy


  subroutine do_final_store (n_vir)
    !
    ! Purpose: stores  the information  for "save_<type>" on  file (if
    ! not already done by do_<type>_store).
    !
    ! INTENT (IN)
    ! integer (i4_kind) :: loop
    ! real   (r8_kind) :: tot_en
    ! character (len=*) :: fit_file, scf_file, ham_file, eig_file
    use comm, only: comm_rank
    use occupation_module, only: eigenstates_store
    implicit none
    integer (i4_kind), intent (in) :: n_vir
    !** End of interface ***************************************

    integer :: save_interval
    logical :: stored_curr, stored_prev

    if (comm_rank() /= 0) return

    save_interval = options_save_interval()

    stored_curr = store_now_p (loop, save_interval)
    stored_prev = store_now_p (loop - 1, save_interval)

    ! Calling convergence_state_store with mode=recover_nothings means
    ! that   no    data   has    been   saved   previously    by   any
    ! convergence_state_store command but  the kept convergence buffer
    ! from the  read_scfcontrol command are  to be considered  for the
    ! storing operation.
    !
    ! Older versions  used to print a warining  that convergence state
    ! may be incomplete  if the program noticed that  the user changed
    ! any options by way of scfcontrol file. If you change any options
    ! during SCF you must know what you are doing.

    if ((options_save_fitcoeff() .and. .not. stored_curr)) then
       call do_fitcoeff_store (store_now=.true., &
            mode=recover_fitcoeff, tot_en=tot_en)
       if (output_data_saved) then
          write (output_unit, '(/ a     )') 'Stored total energy :'
          write (output_unit, '(  a     )') 'tot_en'
          write (output_unit, '(4es20.13)') tot_en
       endif
    endif

    if ((options_save_scfstate() .and. .not. stored_curr)) then
       call do_scfstate_store (store_now=.true., &
            mode=recover_scfstate, tot_en=tot_en)
       if (output_data_saved) then
          write (output_unit, '(/ a     )') 'Stored total energy :'
          write (output_unit, '(  a     )') 'tot_en'
          write (output_unit, '(4es20.13)') tot_en
       endif
    endif

    if ((options_save_ksmatrix() .and. .not. stored_curr)) then
       call do_ksmatrix_store (store_now=.true., mode=recover_nothing)
    endif

    if ((options_save_eigenvec() .and. .not. stored_prev)) then
       call do_eigenvec_store (store_now=.true., n_vir=n_vir, &
            mode=recover_eigenvec)
    endif

    if (options_save_as_fragment()) then
       block
          type (readwriteblocked_tapehandle) :: th

          call say ("do_save_as_fragment")
          call readwriteblocked_startwrite (recfile (frg_file), th, variable_length=.true.)
          call eigenstates_store (th=th)
          call readwriteblocked_stopwrite (th)
       end block
    endif
    ! Deallocate  the temporarily kept  convergence buffers  and reset
    ! the mixing state buffers (if required)
    call convergence_clear_buffers
    call mixing_reset_buffers (loop)

  end subroutine do_final_store

  subroutine write_trace_header
    use comm, only: comm_rank
    implicit none
    ! *** end of interface ***

    integer (i4_kind) :: trace_bgn, trace_end

    ! FIXME: some of the inquiry functions, notably
    ! occupation_spindiff() fail on slaves:
    if (comm_rank() /= 0) return

    trace_line = " "
    trace_end = 0
    ! loop
    trace_width(1) = 6 ! sum (trace_width) = 6
    trace_format(1) = '(I6, 1X)'
    trace_bgn = trace_end + 1
    trace_end = trace_end + trace_width(1)
    trace_line(trace_bgn:trace_end) = "  Loop"
    ! energy
    trace_width(2) = 21 ! sum (trace_width) = 27
    trace_format(2) = '(F20.11, 1X)'
    trace_bgn = trace_end + 1
    trace_end = trace_end + trace_width(2)
    trace_line(trace_bgn:trace_end) = "       Energy        "
    ! entropy
    if (fermi_level_broad) then
       trace_width(3) = 12 ! sum (trace_width) = 39
       trace_format(3) = '(F11.6, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(3)
       trace_line(trace_bgn:trace_end) = " Entr.cont. "
    endif
    ! deviation of charge fit coefficients
    if (convergence_check_coeff_dev()) then
       trace_width(4) = 8 ! sum (trace_width) =  47
       trace_format(4) = '(ES7.1, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(4)
       trace_line(trace_bgn:trace_end) = " |Da_k| "
    endif
    ! deviation of the fitted charge density
    if (convergence_check_coulomb_dev()) then
       trace_width(5) = 8 ! sum (trace_width) = 55
       trace_format(5) = '(ES7.1, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(5)
       trace_line(trace_bgn:trace_end) = " |Drho| "
    endif
    ! deviation of the density matrix
    if (convergence_check_density_dev()) then
       trace_width(6) = 8 ! sum (trace_width) = 63
       trace_format(6) = '(ES7.1, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(6)
       trace_line(trace_bgn:trace_end) = " |DPij| "
    endif
    ! spin difference
    if (symmetry_data_n_spin() == 2) then
       trace_width(7) = 9 ! sum (trace_width) = 72
       trace_format(7) = '(F8.5, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(7)
       trace_line(trace_bgn:trace_end) = "  Q_spin "
    endif
    ! spin orbit
    if (options_spin_orbit) then
       ! spin difference
       trace_width(8) = 9 ! sum (trace_width) = 72
       trace_format(8) = '(F8.5, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(8)
       trace_line(trace_bgn:trace_end) = "  Q_spin "
       ! total Jz
       trace_width(9) = 9 ! sum (trace_width) = 81
       trace_format(9) = '(F8.5, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(9)
       trace_line(trace_bgn:trace_end) = "  Jz     "
    endif
    ! mixing_ch
    trace_width(10) = 6 ! sum (trace_width) = 78
    trace_format(10) = '(F5.2, 1X)'
    trace_bgn = trace_end + 1
    trace_end = trace_end + trace_width(10)
    trace_line(trace_bgn:trace_end) = "chmix "
    ! mixing_sp
    if ((options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda) .and. &
         symmetry_data_n_spin() == 2) then
       trace_width(11) = 6 ! sum (trace_width) = 84
       trace_format(11) = '(F5.2, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(11)
       trace_line(trace_bgn:trace_end) = "spmix "
    endif
    ! mixing_xc
    if (options_xcmode() == xcmode_exchange_fit) then
       trace_width(12) = 6 ! sum (trace_width) = 90
       trace_format(12) = '(F5.2, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(12)
       trace_line(trace_bgn:trace_end) = "xcmix "
    endif
    ! Q solvation mixing
    if (operations_solvation_effect .and. mix_charges()) then
       trace_width(13) = 6 ! sum (trace_width) = 96
       trace_format(13) = '(F5.2, 1X)'
       trace_bgn = trace_end + 1
       trace_end = trace_end + trace_width(13)
       trace_line(trace_bgn:trace_end) = "Qs_mix"
    end if
    ! timing
    trace_width(14) = 7 ! sum (trace_width) = 103
    trace_format(14) = '(F7.1)'
    trace_bgn = trace_end + 1
    trace_end = trace_end + trace_width(14)
    trace_line(trace_bgn:trace_end) = "   Time"
    ! final writing
    call write_to_trace_unit (repeat ("-", trace_end))
    call write_to_trace_unit (trace_line(1:trace_end))
    call write_to_trace_unit (repeat ("-", trace_end))

    new_trace_header = .false.

  end subroutine write_trace_header

  subroutine write_trace_item ()
    use comm, only: comm_rank
    use occupation_module, only: occupation_spindiff, occupation_jz
    implicit none
    ! *** end of interface ***

    integer (i4_kind) :: trace_bgn, trace_end

    ! FIXME: some of the inquiry functions, notably
    ! occupation_spindiff() fail on slaves:
    if (comm_rank() /= 0) return

     trace_end = 0
     ! loop
     trace_bgn = trace_end + 1
     trace_end = trace_end + trace_width(1)
     write (trace_line(trace_bgn:trace_end), trace_format(1)) loop
     ! energy
     trace_bgn = trace_end + 1
     trace_end = trace_end + trace_width(2)
     write (trace_line(trace_bgn:trace_end), trace_format(2)) tot_en
     if (fermi_level_broad) then
        ! entropy
        entropy = fermi_get_entropy()
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(3)
        write (trace_line(trace_bgn:trace_end), trace_format(3)) entropy
     endif
     if (convergence_check_coeff_dev()) then
        ! deviation of charge fit coefficients
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(4)
        write (trace_line(trace_bgn:trace_end), trace_format(4)) coeff_dev
     endif
     if (convergence_check_coulomb_dev()) then
        ! deviation of the fitted charge density
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(5)
        write (trace_line(trace_bgn:trace_end), trace_format(5)) coulomb_dev
     endif
     if (convergence_check_density_dev()) then
        ! deviation of the density matrix
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(6)
        write (trace_line(trace_bgn:trace_end), trace_format(6)) density_dev
     endif
     if (symmetry_data_n_spin() == 2) then
        ! spin difference
        spin_diff = occupation_spindiff() ! FIXME: fails on slaves
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(7)
        write (trace_line(trace_bgn:trace_end), trace_format(7)) spin_diff
     endif
     if (options_spin_orbit) then
        ! spin difference
        spin_diff = s_average
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(8)
        write (trace_line(trace_bgn:trace_end), trace_format(8)) spin_diff
        ! total Jz
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(9)
        write (trace_line(trace_bgn:trace_end), trace_format(9)) occupation_jz()
     endif
     !mixing_ch
     trace_bgn = trace_end + 1
     trace_end = trace_end + trace_width(10)
     write (trace_line(trace_bgn:trace_end), trace_format(10)) mixing_beta_ch()
     !mixing_sp
     if ((options_xcmode() == xcmode_model_density .or. &
            options_xcmode() == xcmode_extended_mda) .and. &
          symmetry_data_n_spin() == 2) then
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(11)
        write (trace_line(trace_bgn:trace_end), trace_format(11)) mixing_beta_sp()
     endif
     !mixing_xc
     if (options_xcmode() == xcmode_exchange_fit) then
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(12)
        write (trace_line(trace_bgn:trace_end), trace_format(12)) mixing_beta_xc()
     endif
     ! Q solvation mixing
     if (operations_solvation_effect .and. mix_charges()) then
        trace_bgn = trace_end + 1
        trace_end = trace_end + trace_width(13)
        write (trace_line(trace_bgn:trace_end), trace_format(13)) Qs_mix
     end if
     ! timing
     trace_bgn = trace_end + 1
     trace_end = trace_end + trace_width(14)
     write (trace_line(trace_bgn:trace_end), trace_format(14)) time_for_last_cycle

     ! final writing
     call write_to_trace_unit (trace_line(1:trace_end))

  end subroutine write_trace_item

  subroutine say (phrase)
    use output_module, only: output_main_scf
    use iounitadmin_module, only: write_to_output_units
    implicit none
    character (len=*), intent (in) :: phrase
    ! *** end of interface ***

    if (output_main_scf) then
        call write_to_output_units ("MAIN_SCF: "//phrase)
    endif
  end subroutine say

end subroutine main_scf


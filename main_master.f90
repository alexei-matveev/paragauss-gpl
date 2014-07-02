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
subroutine main_master()
!
!  This routine encodes the MASTER PLAN for all workers.  Historically
!  it was executed only by  the rank-0 worker (usually called "master"
!  processor), hence  the name.  By now  it is executed  in a parallel
!  context.
!
!  The routines which divide the LCGTO into its major parts are called
!  from this level:
!
!  (1) main_symm() -> symmetry_part
!
!  (2) main_integral()
!
!  (3) main_scf() -> does the SCF cycle including Hamiltonian,
!      eigenvalue solution, Reoccupation of levels ...
!
!  (4) post_scf_main() -> does the calculation of the final total
!      xc-energy on the integration grid if operations_gradients is
!      true, also the xc-contributions to the energy gradient will be
!      calculated
!
!  (5) main_gradient() -> calculation of the energy gradients (except
!      xc-part)
!
!  Subroutine called by: main()
!
!
!  Author: TB, FN
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
! Author: TB
! Date:   12/95
! Description: added operations options for steering,
!              debug output to be switched on and off
!              and call to orbital_test.
!              moved call to read_start to read_input.
!
! Modification (Please copy before editing)
! Author: TB
! Date:   5/96
! Description: added call to integral part
!
! Modification (Please copy before editing)
! Author: MS
! Date:   3/97
! Description: added call to main_gradient
!
! Modification (Please copy before editing)
! Author: HH
! Date:   10/97
! Description: added call to response_main()
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: pvm -> comm
!
! Modification (Please copy before editing)
! Author: AS
! Date:   11-12/99
! Description: added calls to calculating solvent effect
!
! Modification (Please copy before editing)
! Author: KN
! Date:   26/7/99
! Description: added call to main_gtensor
!
! Modification
! Author: TS
! Date:   August 09
! Description: added call to empirical_methods (DFT-D)
!
! Modification (Please copy before editing)
! Author:      ...
! Date:        ...
! Description: ...
!----------------------------------------------------------------

# include "def.h"
  use type_module, only: i4_kind, r8_kind
  use operations_module ! defines which operations are to be performed
  use paragauss, only: toggle_legacy_mode
  use filename_module, only: filesystem_is_parallel
  use iounitadmin_module, only: output_unit, stdout_unit, &
       write_to_output_units, write_to_trace_unit
  use options_module, only: options_directaccess_integrals, &
       options_integrals_on_file, update_hessian_iteration
  use time_module, only: start_timer, stop_timer
  use timer_module, only: timer_initialisation, timer_print_summary, &
       timer_print_slavetiming
  use output_module, only: output_timing_detailedsummary, &
       output_timing_slaves, output_timing_summary
  use potential_calc_module, only: charge_constr, esp_map, pdc, &
       use_saved_densmatrix, V_electronic, calc_plane_grid, &
       grid2space_2d, get_poten_and_shutdown_2d, calc_shell_grid, &
       collect_poten_3d, calc_poten_derive_charges
  use integralpar_module, only: integralpar_set, integralpar_cpksdervs, &
       integralpar_int_part_name
  use integralstore_module, only: integralstore_deallocate, &
       integralstore_deallocate_pcm
  use initialization, only: initialize_with_input, finalize_geometry
  use xc_cntrl, only: xc_is_on=>is_on, xc_ANY
  use post_scf_module
  use energy_calc_module, only: write_energies, get_energy
#ifdef FPP_DEBUG
  use error_module, only: MyID
#endif
#ifdef WITH_RESPONSE
  use response_module, only: response_main
#endif
  use efield_module, only: efield_calculate_integrals, efield_applied
  use efield_module, only: efield_intensity, efield_change
  use unique_atom_module, only: unique_atom_iwork
  use unique_atom_methods, only: unique_atom_make_gx
  use occupation_module, only: occupation_symmetry_check
  use convergence_module, only: convergence_max_geo_iteration
  use properties_module, only: properties_main
#ifdef WITH_DFTPU
  use dft_plus_u_module, only: dft_plus_u_output
#endif
#ifdef WITH_EPE
  use epecom_module, only: get_epe_energies, &
       epe_side_optimized_energy_prev, epe_side_optimized_energy
#ifdef NEW_EPE
  use ewaldpc_module, only: epe_relaxation, get_qm_references, qm_ref_run
#else
  use ewaldpc_module, only: epe_relaxation
#endif
#endif
#ifdef WITH_EFP
  use efp_module, only: n_efp, read_gx_qm, def_efp_arrays, calc_X_points, calc_efield_points
  use efp_module, only: print_id, qm_fixed
  use efp_efp_module, only: efp_efp_energy
  use efp_only_opt_module, only: geom_converged
#endif
  use density_data_module, only: open_densmat
  use solv_cavity_module, only: stop_solv
  use potential_module, only: send_recv_space_point
  use elec_static_field_module
  use symmetry, only: main_symm
  use interfaces, only: main_integral
  use interfaces, only: potential_calculate
  use interfaces, only: main_molmech
#ifdef WITH_MOLMECH
  use qmmm_interface_module  !!!!!!!!!!!!!AS
  use qmmm1_interface_module !!!!!!!!!!!!AS
#endif
  use calc3c_switches, only: print_epe
#ifdef WITH_OPTIMIZER
  use opt_data_module, only: filename_setup_opt
#endif
  use induced_dipoles_module, only: calc_Pol_centers, dealloc_pol_center_inform
  USE_MEMLOG
! --------------------------------------------------
  implicit none
  character (len=128) :: version = FPP_PARAGAUSS_VERS
  integer (i4_kind) :: tasks=0  ! Number of tasks LEFT to be done
  integer (i4_kind) :: loop, max_geo_loop ! loop -- Laufvariable fuer Runs
  logical :: geometry_converged = .false. ! determines if Simol is converged or not
#ifdef WITH_EPE
  real (r8_kind) :: energy
  real (r8_kind) :: energy2, epe_latt_energy, cluster_regI
  real (r8_kind) :: eshort
  logical :: epe_side_energy_converged
#endif
  logical :: use_dens_mat


  DPRINT 'main_master: entered'

  !
  ! If  you are  thinking  about adding  code  to be  executed on  all
  ! workers, consider  adding it outside  of the "legacy"  mode blocks
  ! enclosed by
  !
  !   do while (toggle_legacy_mode())
  !      ...
  !   enddo
  !
  ! Alternatively add  it to main.f90 or  pick one of the  subs in the
  ! body such as
  !
  !     - initialize_with_input()
  !     - main_gradient()
  !     - properties_main()
  !     - finalize_geometry()
  !     - ...
  !
  ! that are executed by all workers and augment them.
  !
  ! FIXME: finish  converting to SPMD (single  program, multiple data)
  ! so that all  workers execute the same code.  There  is still a few
  ! cases of  a master-slave mode ---  search for toggle_legacy_mode()
  ! blocks.

  !
  ! Print the  version info and  machine config, uses  output_unit, so
  ! call it after initialize_environment():
  !
  call legal (version)

  call write_to_trace_unit (" ")
  call write_to_trace_unit (" -------------------------------------------")
  call write_to_trace_unit (" ")
  call write_to_trace_unit (" executing program : mainscf_" // trim (version))
  call write_to_trace_unit (" ")
  call write_to_trace_unit (" -------------------------------------------")
  call write_to_trace_unit (" ")

#ifdef WITH_EPE
  epe_side_optimized_energy_prev=0.0_r8_kind
#endif

#ifdef WITH_MOLMECH
  qm_mm_1_task=0  !!!!!!!!!!!AS
#endif

  ![[=== MAIN LOOP OVER TASKS/GEOMETRIES ===================================
2001 CONTINUE ! an "entry point" for the task-loop

  tasks        = 1 ! will be eventualy incremented after read_input!
  max_geo_loop = 1 ! will be eventualy reset to a higher value after read_input!
  loop         = 0
  geometry_converged = .false.

  geometry_loop: do while (tasks > 0)

     loop = loop + 1

     if (loop > max_geo_loop) then
       ! at loop 1 is not true anyway,
       ! after read_input max_geo_loop gets the  proper value ...
       call say ("maximum number of geo interations exceeded")
       tasks = tasks - 1
       cycle geometry_loop
     endif

     ! reset memory couters, if using MEMLOG
     MEMSET (0)

     call start_timer (timer_initialisation)

     ! Otherwise there will be multiple copies printed:
     if (output_unit > 0 .and. stdout_unit > 0) then ! yes, AND!
        call write_to_output_units (" ------------------------------------")
        call write_to_output_units (" -                                  -")
        call write_to_output_units (" - main_master: Run No. ", inte=loop)
        call write_to_output_units (" -                                  -")
        call write_to_output_units (" ------------------------------------")
     endif


     ! TODO: move call read_input() and
     !       subsequent control manupulations out of geometry loop!
     DPRINT 'main_master: call read_input()'
     call read_input (loop)
     DPRINT 'main_master: .'

     ! 7 runs for intensity calculation by finite difference of forces
     ! in presence of electric field.
     if (efield_intensity()) then
       max_geo_loop = 7
       tasks = tasks + 1
       ! change the direction of electric field in all 6 direction i.e +/-X, +/-Y and +/- Z
       call efield_change (loop)
     end if

#ifdef WITH_OPTIMIZER
     call filename_setup_opt (optonly=.false.)
#endif

#ifdef WITH_MOLMECH
     if (operations_qm_mm_new .and. qm_mm_1) then    !!!!!!!!!!!AS
        call def_qm_mm_1_tasks()                    !!!!!!!!!!!AS
     end if
#endif

     DPRINT 'main_master: options_directaccesa_integrals=', options_directaccess_integrals()
     if (options_directaccess_integrals()) then
        ASSERT (.not.filesystem_is_parallel)
     endif

     ! ... and set the variable 'max_geo_loop' according to the input
     if (operations_geo_opt .or. operations_optimizer_only) then
         max_geo_loop = convergence_max_geo_iteration()
     endif

     ![[=== decide wheather or not to do the hessian ===
     select case (update_hessian_iteration)

     case (0)
       ! do nothing, the default

     case (1:)
       ! geometry optimization with regular updates of the hessian
       if (MOD (loop-1, update_hessian_iteration) == 0) then
         ! at loop 1, update_hessian_iteration+1, 2*update_hessian_iteration+1, ...
         call integralpar_set ('SecondDervs')
         ! sets module vars:
         ! integralpar_2dervs    = .true.
         ! integralpar_cpksdervs = .true.
         print *,'main_master: do second derivatives at loop', loop &
                ,'(every', update_hessian_iteration,'iterations)'
       else
         call integralpar_set ('NoSecondDervs')
         ! integralpar_2dervs    = .false.
         ! integralpar_cpksdervs = .false.
       endif

     case (-1)
       ! frequencies after geometry optimization
       if (loop == 1) then
         ! add a task of computing second derivatives after geometry optimization:
         tasks = tasks + 1
       endif
       if (geometry_converged .or. loop == max_geo_loop) then
         ! enable second derivatives only after geometry is converged:
         update_hessian_iteration = 1
         call integralpar_set ('SecondDervs')
         ! sets module vars:
         ! integralpar_2dervs    = .true.
         ! integralpar_cpksdervs = .true.
         ! nothing below seems to depend on it:
         operations_geo_opt    = .false.
         print *,'main_master: do second derivatives at loop', loop &
                ,'(after geometry is converged)'
         ! mark the task completed (a bit early, of couse):
         tasks = tasks - 1
       endif

     case default
       print*,'ERROR: no such update_hessian_iteration=', update_hessian_iteration
       ABORT ('no such case yet!')
     end select
     !]]================================================

#ifdef WITH_EPE
     ! EPE lattice calculations
     if (operations_epe_lattice) then
        call say ("calling EPE-lattice optimization")
        call epe_lattice_optimization()
        call stop_timer (timer_initialisation)
        ! DONT exit geometry_loop
        tasks = tasks - 1
        cycle geometry_loop
     end if
#endif

      if (operations_optimizer_only) then
         call say ("calling opimizer")
         geometry_converged = optimizer_step (1)
         call stop_timer (timer_initialisation)
         ! DONT exit geometry_loop
         tasks = tasks - 1
         cycle geometry_loop
      endif

#ifdef WITH_MOLMECH
     if (operations_mol_mech) then                                !!!!!!!!!!!AS
        call say ("calling molecular mechanics program") !!!!!!!!!!!AS
        call main_molmech (mm_run, unique_atom_iwork, max_geo_loop) !!!!!!!!!!!AS
        if (unique_atom_iwork > 0) then                           !!!!!!!!!!!AS
           operations_geo_opt=.true.                             !!!!!!!!!!!AS
        end if                                                   !!!!!!!!!!!AS
     end if                                                      !!!!!!!!!!!AS
     if (operations_qm_mm_new .and. (qm_mm .or. qm_mm_1)) then    !!!!!!!!!!!AS
        call say ("distributing data for QM+MM run")
        call read_gx_qmmm()
        call main_molmech (qmmm_read_input)
        call qmmm2pc()
     end if                                                      !!!!!!!!!!AS
#endif

#ifdef WITH_EFP
     if (efp) then
        call say ("EFP - reading QM atoms from GX")
        call read_gx_qm()
     end if
#endif

     ! call symmetry part
     if (operations_symm .and. .not. operations_epe_lattice) then
        call say ("calling symmetry part")
        print_epe = loop .eq. 1
        call main_symm()
     endif

#ifdef WITH_EFP
     !generation of arrays of external points for EFP calculation
     if (efp) then
        call say ("EFP - generation of arrays of external points")
        call def_efp_arrays()
        call calc_X_points()
        if (calc_Pol_centers()) call calc_efield_points()
     end if
#endif

     ! write input file if desired
     if (operations_write_input .or. operations_get_input_out) then
        call say ("write_input ")
        call write_input()
        call say ("write_input done")
     endif

     if (operations_get_input_out) then
        call stop_timer (timer_initialisation)
        ! DONT exit geometry_loop
        tasks = tasks - 1
        cycle geometry_loop
     end if

     ! generating suface charge distribution (solvation effect)
     if (operations_solvation_effect) then
        call say ("call build_mol_surfaces()")
        do while (toggle_legacy_mode())
           call build_mol_surfaces()
        enddo
        if (stop_solv) then
           call stop_timer (timer_initialisation)
           ! DONT exit geometry_loop
           tasks = tasks - 1
           cycle geometry_loop
        end if
     endif

#ifdef WITH_EFP
     !calculating interfragment interactions (EFP)
     if (efp .and. n_efp > 0) then
        call say ("EFP - calculating interfragment interactions")
        call efp_efp_energy (print_id)
     end if
#endif

     ! write gx file
     if (operations_make_gx) then
        call say ("unique_atoms_make_gx")
#ifndef WITH_EPE
        ! FIXME: is it upto date?
        call unique_atom_make_gx (iloop= 1)
#else
        call unique_atom_make_epegx()
#endif
     end if

     !
     ! Initialize   various   modules    in   part   by   broadcasting
     ! initialisation  information taken  from  input file.   It is  a
     ! convenient place to put the  code to be executed on all workers
     ! ...
     !
     call say ("initialize_with_input")
     call initialize_with_input()

     call stop_timer (timer_initialisation)

     !
     ! At this point all information should be in place that is needed
     ! to check the fixed occupation numbers.
     !
     ! Check the 'n_nonempty_irreps'-specification in
     ! the namelist OCCUPATION with the variable 'symmetry_data_n_irreps'
     ! after the symmetry part has been run:
     !
     call occupation_symmetry_check()

     use_dens_mat = .false.
     if (operations_potential .and. esp_map) then
        if (.not. V_electronic .or. use_saved_densmatrix) goto 1111
     elseif (operations_potential .and. pdc) then
        if (charge_constr .and. use_saved_densmatrix) then
           use_dens_mat = .true.
           goto 1111
        endif
     endif

     ! do integral part
     if (operations_integral) then
        call say ("Starting the Integral Part")
        !
        ! Not all integrals are needed in a property run without prior
        ! SCF.  On the other hand if SCF is performed anyway, no extra
        ! integrals are needed for properties (that is the theory):
        !
        if (operations_scf) then
           call say ("call integralpar_set (Normal)")
           call integralpar_set ('Normal')
        elseif (operations_properties) then
           call say ("call integralpar_set (Properties)")
           call integralpar_set ('Properties')
        else
           ABORT ("ever happens?")
        endif
        call say ("done")

        call main_integral ()
        call say ("done with Integral Part")
     end if

     ! calculate integrals of external electrical field if one is applied
     if (efield_applied()) then
        call say ("Calculating integrals of external electrical field")

        call efield_calculate_integrals ()
     endif


     ! calculate integrals of electrostatic potential
     if (operations_solvation_effect .and. operations_integral) then
        call say ("Calculate integrals of electrostatic potential ...")
        call send_recv_space_point()
        call potential_calculate ('Solvation')
        call say ("Done with integrals of electrostatic potential.")
     endif

#ifdef WITH_EFP
     if (calc_Pol_centers() .and. operations_integral) then
        call field_calculate ()
     end if
#endif

     ! do scf part
     if (operations_scf) then
        MEMSET (0)
        call say ("Starting the main SCF routine ...")
        call main_scf()
        call say ("Done with the main SCF routine.")
        MEMSET (0)
     endif

     ! Calculate and print dipole moments:
     if (operations_dipole .or. operations_gtensor .or. operations_hfc) then
        call say ("call main_dipole()")
        call main_dipole()
        call say ("done main_dipole()")
     endif

     ! Potential derived charges here (not the same as solvation, even
     ! though  the same  buildign blocks  are used).  FIXME:  too much
     ! logic/branches for  a master plan.   Leave here only  the entry
     ! points, move logic to the respective modules:
1111 if (operations_potential) then
        call say ("Starting the potential routines ...")
        if (esp_map) then
           do while (toggle_legacy_mode())
              call calc_plane_grid()
           enddo
           call grid2space_2d()
           if (V_electronic) then
              call say ("call potential_calculate (Vel)")
              if (use_saved_densmatrix) then
                 call open_densmat() ! no comm, reads disk
              endif
              call potential_calculate ('Potential')
           endif
           call get_poten_and_shutdown_2d()
           if (.not. V_electronic) then
             tasks = tasks - 1
             cycle geometry_loop
           endif
        elseif (pdc) then
           call say ("call potential_calculate (PDC)")
           call calc_shell_grid()
           if (use_dens_mat) then
              call open_densmat() ! no comm, reads disk
           endif
           call potential_calculate ('Potential')
           call collect_poten_3d()
           do while (toggle_legacy_mode())
              call calc_poten_derive_charges()
           enddo
        endif
        call say ("Done with the potential routines.")
     endif

#ifdef WITH_DFTPU
     !
     ! DFT+U output, does nothing if not in use:
     !
     call dft_plus_u_output ()
#endif

     ! Do  properties   before  deallocating  the   integral  storage.
     ! Otherwise a  read_overlap() called from  properties_main() will
     ! fail. Slaves execute properties_main() too:
     if (operations_properties) then
        call say ("call properties_main()")
        call properties_main()
        call say ("done properties_main()")
     end if

     ! Kin and Nuc are now always allocated, dellaocate:
     call integralstore_deallocate (deallocate_kin=.true., deallocate_nuc=.true.)

     if (.not. options_integrals_on_file() .and. .not. integralpar_cpksdervs) then
        ! if integralpar_cpksdervs -- deallocate after
        ! cpks_g4constructs()

        ! Deallocate integral storage ...
        call integralstore_deallocate()

        ! Deallocate PCM-integral storage (if been allocated) ...
        call integralstore_deallocate_pcm()
     end if


     ! FIXME: please specify the goal of GOTO by words as well,
     !        not just by the number!
     ! GO TO: .... because ...
     if (operations_potential .and. use_saved_densmatrix) goto 1112

     ! do Post Scf calculation of Exchange Energy
     if (operations_post_scf) then

        ! only if XC /= off:
        if (xc_is_on (xc_ANY)) then
           call say ("call post_scf_main()")
           call post_scf_main()
           call say ("done post_scf_main()")
        endif

        call say ("write_energies")
        call write_energies (output_unit, post_scf=.true.)
     end if

#ifdef WITH_MOLMECH
     !Calculations electrostatic field produced by QM cluster
     !at MM atoms
     if (operations_qm_mm_new .and. qm_mm) then !!!!!!!!!!!!AS
        call QMfield_at_mm_points()             !!!!!!!!!!!!AS
     end if                                     !!!!!!!!!!!!AS

#endif

     ! do PostSCF calculation of matrix elements needed for
     ! response calculations with tdfrt response program
     if (operations_response) then
#ifdef WITH_RESPONSE
        call say ("call response_main()")
        call response_main ()
        call say ("done response_main()")
#else
        ABORT ('recompile -DWITH_RESPONSE')
#endif
     endif

     MEMSET (0)
#ifndef WITH_EPE
     if (operations_gradients) then
        call say ("call main_gradient()")
        call main_gradient (loop) ! (1)
        call say ("done main_gradient()")
     endif
#else
     ! do EPE calculations, FIXME: why inlining so much staff into this high
     ! level sub?
     if (operations_qm_epe .and. epe_relaxation) then
        call say ("main_epe_block")
        call main_epe_block()

        call get_energy (tot=energy)
        energy2 = energy
        call get_epe_energies (lattice_energy=epe_latt_energy, &
             epg_cluster_reg_I=cluster_regI, eshort_coupling_au=eshort)
        energy = energy + epe_latt_energy
        print*,'energy, energy2, eshort', energy, energy2, eshort
        print*,'epe_side_optimized_energy', energy+eshort
        epe_side_optimized_energy = energy + eshort
        print*,'cluster_regI', cluster_regI
        print*,'epe_latt_energy', epe_latt_energy

        call write_to_trace_unit ('epe_convergence_check')
        call epe_convergence_check (epe_side_energy_converged, loop)
     endif

     ! Calculate gradients
     if (operations_gradients .and. .not. epe_relaxation .or. &
          epe_relaxation .and. epe_side_energy_converged) then
        !
        ! Regular branch here!
        !
        call say ("Starting main_gradient() ...")
        if (epe_relaxation .and. epe_side_energy_converged) then
           call write_to_trace_unit ('epe_relaxation .and. epe_side_energy_converged')
        endif

        call main_gradient (loop) ! (2)
        call say ("Done with the integral part for gradients routine.")
     elseif (operations_gradients) then
        !
        ! FIXME: clean up is the task of finalize_geometry()
        !        that is called anyway. Why doing it here?
        !
        ABORT ('please adapt')
     endif
#endif

#ifdef WITH_MOLMECH
     if (operations_qm_mm_new) then  !!!!!!!!!!!!!!AS
        if (imomm) then
           call say ("calling molecular mechanics module to perform IMOMM job")
           call main_molmech (imomm_mm_small)
           call main_molmech (imomm_mm_large)
           call say ("QMMM job: summing up different gradients and writing gxfile")
           call sum_up_grads_and_write_gx()
        elseif (qm_mm .or. (qm_mm_1 .and. qm_mm_1_task ==0)) then
           call say ("calling molecular mechanics module to perform QM_MM (1) job")
           call main_molmech (qm_mm_run)
           call write_gx_qmmm()
        end if
     endif                           !!!!!!!!!!!!!!AS
#endif

        MEMSET (0)

     !
     ! Various shutdown and deallocation work:
     !
     call finalize_geometry()

1112 continue

     ! Otherwise there will be multiple copies printed:
     if (output_unit > 0 .and. stdout_unit > 0) then ! yes, AND!
        call write_to_output_units (" ------------------------------------")
        call write_to_output_units (" -                                  -")
        call write_to_output_units (" - main_master: End of Run No. ", inte=loop)
        call write_to_output_units (" -                                  -")
        call write_to_output_units (" ------------------------------------")
     endif

     !
     ! If max_geo_iteration was set to zero in the input,
     ! dont even try to run optimizer. This will not overwrite
     ! gxfile with the new geometry updated by optimizer algorithm.
     ! May be usefull for work with alternative external optimizers:
     !
     if (operations_geo_opt .and. loop > max_geo_loop) then
       WARN ('loop leaped beyond max_geo_loop')
       tasks = tasks - 1
       exit geometry_loop
       ! there is not point to cycle geometry_loop so far as the first
       ! condition it checks after entry is again "loop > max_geo_loop"
     endif

     if (operations_geo_opt .or. operations_gx_test) then
#ifdef WITH_EFP
        if (efp .and. qm_fixed) then
           if (geom_converged (loop)) then
              tasks = tasks - 1
           end if
           cycle geometry_loop
        end if
#endif
        call say ("call optimizer_step()")
        geometry_converged = optimizer_step (unique_atom_iwork)

        if (geometry_converged) then
           call say ("Geometry converged")
           ! DONT exit geometry_loop
           tasks = tasks - 1
           ! nothing below seems to depend on it:
           ! operations_geo_opt = .false.
           ! operations_gx_test = .false.
           ! But so far it doesnt work as the read_input()
           ! resets them to true again!
        endif
#ifdef NEW_EPE
        if (operations_qm_epe .and. geometry_converged) then
           if (get_qm_references .and. .not. qm_ref_run) then
              tasks = tasks + 1
              qm_ref_run=.true.
              !one more cycle to calculate and save
              !epe_reference and pg_epe_reference
           end if
        endif
#endif
     else
        ! a default calculation (energy/gradients) requires one run:
        tasks = tasks - 1
     endif

  enddo geometry_loop

#ifdef WITH_MOLMECH
  if (operations_qm_mm_new .and. qm_mm_1) then
     qm_mm_1_task=qm_mm_1_task+1
     ! FIXME: I need to repeat all once again because ...
     if (qm_mm_1_task == 1) goto 2001 ! enter the task loop again
     ! QUESTION: is this logic only to enter geometry_loop twice?
  end if
#endif
  !]]=== eof MAIN LOOP OVER TASKS/GEOMETRIES ================================

  ! print timing
  if (output_timing_summary .or. output_timing_slaves .or. output_timing_detailedsummary) then
    call say ("printing timing")
    if (output_timing_summary .or. output_timing_detailedsummary) &
         call timer_print_summary (integralpar_int_part_name)
    if (output_timing_slaves) &
         call timer_print_slavetiming (integralpar_int_part_name)
  endif

  call say ("done")

contains

#ifdef WITH_EPE
  subroutine epe_convergence_check (epe_side_energy_converged, i_iter)
    use epecom_module, only: epe_rel_converged, &
         epe_side_optimized_energy, epeside_energy_limit, &
         epe_side_optimized_energy_prev, epe_basic_action=>basic_action
    use filename_module, only: data_dir
    implicit none
    logical, intent (out) :: epe_side_energy_converged
    integer (i4_kind), intent (in) :: i_iter
    ! *** end of interface ***

    if (epe_relaxation .and. .not. epe_basic_action .eq. 0) then
       inquire (file=trim (data_dir) // "/epe_rel_unconverged", &
            exist=epe_rel_converged)
       epe_rel_converged = .not. epe_rel_converged
       print *, 'epe_convergence_check: epe_relaxation, epe_rel_converged', &
            epe_relaxation, epe_rel_converged
       epe_side_energy_converged = abs (epe_side_optimized_energy_prev - &
            epe_side_optimized_energy) .lt. 0.00002 !?  embed_convergence_limit

       if (.not. epe_side_energy_converged) &
            epe_side_energy_converged = abs (epe_side_optimized_energy_prev - &
            epe_side_optimized_energy) .lt. 0.0004 .and. i_iter .gt. 10

       print *, 'epe_side_energy_jump' , &
            epe_side_optimized_energy_prev - epe_side_optimized_energy, &
            epeside_energy_limit
       epe_side_optimized_energy_prev = epe_side_optimized_energy
       if (.not. epe_rel_converged .or. .not. epe_side_energy_converged) then
          print *, 'optimizer: epe relaxation is still not converged', &
               epe_rel_converged, epe_side_energy_converged, &
               epe_side_optimized_energy_prev - &
               epe_side_optimized_energy
          call write_to_trace_unit ('.not. epe_side_energy_converged')
       else
          call write_to_trace_unit ('epe_side_energy_converged')
       endif
    else
       epe_side_energy_converged = epe_basic_action .eq. 0
    endif
  end subroutine epe_convergence_check
#endif

  function optimizer_step (geo_loop) result (converged)
    !
    ! Executed  by  all workers,  though  most  of  the work  is  done
    ! serially.  The input "geo_loop"  schould be the same everywhere.
    ! The output will be the same on all workers.
    !
    use comm, only: comm_rank, comm_bcast, comm_barrier
    implicit none
    integer (i4_kind), intent (in) :: geo_loop
    logical :: converged
    ! *** end of interface ***

    call comm_barrier()         ! paranoya

    ! One of us does the work that includes some file-system mangling:
    if (comm_rank() == 0) then
       call do_optimizer_step (geo_loop, converged)
    endif

    ! Broadcast the result:
    call comm_bcast (converged)
    call comm_barrier()         ! paranoya
  end function optimizer_step

  subroutine do_optimizer_step (geo_loop, geo_conv)
    !
    ! FIXME: down  the call chain  there is some  file-system mangling
    ! such as  (re)writing gx- and  hessian files.  This code  may not
    ! work as intended if executed by more than one worker.
    !
#ifdef WITH_EPE
    use epecom_module, only: cross_boundary_3b, epe_rel_converged, &
         epe_basic_action => basic_action
#endif
#ifdef WITH_OPTIMIZER
    use optimizer, only: main_opt
    use operations_module, only: operations_task, namelist_tasks_used
#endif
    implicit none
    integer (i4_kind), intent (in) :: geo_loop
    logical, intent (out) :: geo_conv
    ! *** end of interface ***

    logical :: conv, stop_after_eperelaxation, convert_internal
    character (len=32) :: optimizer_task='GeoOpt'

    geo_conv = .false.
#ifdef WITH_EPE
    if (operations_qm_epe) then
       stop_after_eperelaxation=.true.
       if (epe_relaxation .and. .not. epe_basic_action .eq. 0) then
          if (.not. epe_rel_converged .or. .not. epe_side_energy_converged) then
             print*,'optimizer: epe relaxation is still not converged'
             call write_to_trace_unit ('optimizer epe relaxation is still not converged')
             return
          endif
       endif
    end if
#endif
#ifdef WITH_OPTIMIZER
    if (namelist_tasks_used) then
       optimizer_task=operations_task
    else
       optimizer_task='GeoOpt'
       WARN ('assuming GeoOpt in optimizer')
    endif
    print*,'optimizer: call main_opt(', optimizer_task,')'
    convert_internal=.false.
#ifdef WITH_EPE
    if (epe_relaxation) then
    call main_opt (task=optimizer_task , converged=conv       &
                 , stop_after_eperelaxation=stop_after_eperelaxation &
                 , convert_internal=convert_internal &
                 , cross_boundary_3b=cross_boundary_3b)
    else
#endif
    call main_opt (task=optimizer_task , converged=conv       &
                 , stop_after_eperelaxation=stop_after_eperelaxation &
                 , convert_internal=convert_internal)
#ifdef WITH_EPE
    endif
#endif
    if (convert_internal) &
    call main_opt (task=optimizer_task , converged=conv       &
                 , stop_after_eperelaxation=stop_after_eperelaxation &
                 , convert_internal=convert_internal)
    ! returns converged=true if converged
#else
    DPRINT 'optimizer: DONT call main_opt()'
    ABORT ('recompile with -DWITH_OPTIMIZER')
#endif
#ifdef WITH_EPE
    if (operations_qm_epe) then
       if (epe_relaxation .and. stop_after_eperelaxation) then
          conv=.true.
       endif
    end if
#endif

    ! set output flag:
    geo_conv = conv

    if (.not. conv) then
       call write_to_output_units&
            ("optimizer: geometry  not yet converged in loop", inte =geo_loop)
       call write_to_trace_unit&
            ("optimizer: geometry  not yet converged in loop", inte =geo_loop)
    endif
  end subroutine do_optimizer_step

  subroutine legal (version)
    use comm_module, only: comm_print_conf
    implicit none
    character (len=*), intent (in) :: version
    ! *** end of interface ***

    ! Slaves  may need  to have  output_unit <  0,  otherwise multiple
    ! copies of the header will be printed:
    if (output_unit <= 0) return

    !
    ! This prints machine config, both to output and to tty:
    !
    call comm_print_conf (output_unit)
    ASSERT(stdout_unit>0)
    call comm_print_conf (stdout_unit)

    !
    ! Legal header:
    !
    call write_to_output_units ('########################################################################')
    call write_to_output_units (' ')
    call write_to_output_units ('ParaGauss GPL ' // trim (version))
    call write_to_output_units (' ')
    call write_to_output_units ('features:')
#ifdef WITH_LIBDFTAUTO
    call write_to_output_units ('* [!] xc functionals from http://www.cse.clrc.ac.uk/qcg')
#endif
#ifdef NEW_INTEGRALS
    call write_to_output_units ('* [+] NEW INTEGRALS compiled in')
#else
    call write_to_output_units ('* [-] NEW INTEGRALS disabled')
#endif
#ifdef WITH_OPTIMIZER
    call write_to_output_units ('* [+] OPTIMIZER compiled in')
#else
    call write_to_output_units ('* [-] OPTIMIZER disabled')
#endif
#ifdef WITH_GTENSOR
    call write_to_output_units ('* [+] GTEN/HFCC support compiled in')
#else
    call write_to_output_units ('* [-] GTEN/HFCC disabled')
#endif
#if WITH_RELFIT || WITH_SHGI
    call write_to_output_units ('* [+] RELFIT support compiled in')
#else
    call write_to_output_units ('* [-] RELFIT disabled')
#endif
#ifdef WITH_EPE
    call write_to_output_units ('* [+] EPE support compiled in')
#else
    call write_to_output_units ('* [-] EPE disabled')
#endif
#ifdef WITH_EFP
    call write_to_output_units ('* [+] EFP support compiled in')
#else
    call write_to_output_units ('* [-] EFP disabled')
#endif
#ifdef WITH_MOLMECH
    call write_to_output_units ('* [+] MOLMECH support compiled in')
#else
    call write_to_output_units ('* [-] MOLMECH disabled')
#endif
#ifdef WITH_RESPONSE
    call write_to_output_units ('* [+] RESPONSE support compiled in')
#else
    call write_to_output_units ('* [-] RESPONSE disabled')
#endif
#ifdef WITH_ERI4C
    call write_to_output_units ('* [+] ERI4C support compiled in')
#else
    call write_to_output_units ('* [-] ERI4C disabled')
#endif
    call write_to_output_units ('########################################################################')
  end subroutine legal

  subroutine say (phrase)
    use output_module, only: output_main_master
    use iounitadmin_module, only: write_to_output_units
    implicit none
    character (len=*), intent (in) :: phrase
    ! *** end of interface ***

    if (output_main_master) then
        call write_to_output_units ("main_master: "//phrase)
    endif
  end subroutine say

end subroutine main_master

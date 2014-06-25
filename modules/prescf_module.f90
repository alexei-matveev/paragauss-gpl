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
module  prescf_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: contains all routines that perform preparatory
  !           tasks for the SCF-cyclces such as resetting
  !           I/O files, reading the kinetic and nuclear hamiltonian,
  !           setting up the initial density matrix etc.
  !
  !  Module called by: main_scf
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !  Contents:
  !     Routines:
  !  (i)  PUBLIC    reset_IO_files  -> do a simple open- and
  !                                    close statement on certain
  !                                    files to ensure they exist
  !                                    and are empty
  !                     prescf -> call all routines necessary for the
  !                               preparation on SCF-cycles
  !
  !  (ii) PRIVATE   pre_dens_master -> perform all necessary routines
  !                                    for setting up an initial density
  !                                    matrix. Options are
  !                                    * arbitrary densmat
  !                                    * inout from old lcgto
  !                                    The decicion is currently set in the
  !                                    DOIT- Routines (conversion from old to
  !                                    new lcgto)
  !                     read_dens ->   read in initial eigenvalues (c_dens) either
  !                                    from file or arbitrary values. This is hard-
  !                                    coded in this routine. No decision in DOIT

  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: pvm -> comm
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   12/99
  ! Description: calling "const_part_of_ham_and_energy" was added (solv. effect)
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters

  implicit none
  private
  save

  !------------ public functions and subroutines ------------------

  public :: prescf_init!(), is called by all workers
  public :: prescf_finalize!(), is called by all workers

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine reset_IO_files
    !  Purpose: reset the I/O-files that store information according
    !           to the ouput level.
    !           All files are to be found on $TTFSOUT.
    !           Since the files will be used with the 'append'-position
    !           they have to be opened and closed with the
    !           'replace'-status in the beginning to ensure they are
    !           existing and empty.
    !           This routine may be varied according to the actual needs
    !           concerning the output.
    !  List of Files (last modified :9.96)
    !  - eigvec.dat
    !  - occupation.out
    !  - hamiltonian
    !  - coeff_charge
    !  - densmat.new
    !** End of interface *****************************************
    ! ------------ Modules ---------------------------------------
    use filename_module, only: outfile
    use iounitadmin_module
    use output_module
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: io_u

    !------------ Executable code --------------------------------

    if (output_reoccup .or. output_eigendata) then
       io_u = get_iounit()
       open(io_u,status='replace',position='append',&
            file=trim(outfile('occupation.out')))
       call return_iounit(io_u)
       close(io_u)
    endif

    if (output_eigendata) then
       io_u = get_iounit()
       open(io_u,status='replace',file= &
            trim(outfile('eigvec.dat')))
       call return_iounit(io_u)
       close(io_u)
    endif

    if (output_hamiltonian) then
       io_u = get_iounit()
       open(io_u,status='replace',file= &
            trim(outfile('hamiltonian')))
       call return_iounit(io_u)
       close(io_u)
    endif

    if (output_chargefit) then
       io_u = get_iounit()
       open(io_u,status='replace',file= &
            trim(outfile('coeff_charge')))
       call return_iounit(io_u)
       close(io_u)
    endif

    if (output_densmat) then
       io_u = get_iounit()
       open(io_u,status='replace',file= &
            trim(outfile('densmat.new')))
       call return_iounit(io_u)
       close(io_u)
    endif

  end subroutine reset_IO_files

  !*************************************************************

  subroutine prescf_init()
    !
    ! Purpose: do preparations for the SCF cycles
    !
    ! (1) calculate the electrostatic interaction energy
    !     between all nuclei (later: also point charges)
    ! (2) calculate one center charge fit integrals if not
    !     already done by integral part
    ! (2a) read in charge overlap matrix, and
    !      if perturbation theory is not turned on or if the model
    !      density approach is used, decompose F and compute F^-1 Q
    ! (3) if a calculation with the ParaGAU integral part is
    !     to be done, the fitcoefficients are initialized
    !     appropriately. Otherwise they are read in from the
    !     old lcgto.
    ! (4) read the kinetic and nuclear part of Hamiltonian
    !     as it comes from the old lcgto, or from the
    !     new integral part. HAM_KIN and HAM_NUC
    !     are then 'public' in the ham_kin_nuc_module
    ! (5) setup bounds and distribute them among the slaves
    !
    ! Executed in parallel context.
    !
    ! author  : Folke Noertemann
    ! date    : 10/95
    !
    !------------ Modules used in this routine  ------------------
    use comm, only: comm_rank
    use energy_calc_module, only: core_interaction_calc
    use options_module, only: options_integrals_on_file
    use fit_coeff_module, only: fit_coeff_normalize, fit_coeff_initialize, &
         fit_coeff_calc_chargefitnorm, fit, get_fit
    use operations_module, only: operations_integral, &
         operations_solvation_effect
    use symmetry_data_module, only: ssym
    use hamiltonian_module, only: hamiltonian_setup
    use mat_charge_module, only: read_mat_charge
    use eigen_data_module, only: eigen_data_alloc
    use bounds_module, only: bounds_calc
    use solv_electrostat_module, only: const_part_of_ham_and_energy
    use induced_dipoles_module, only: calc_Pol_centers
    use elec_static_field_module, only: get_field_nuc
#ifdef WITH_EFP
    use efp_module, only: n_efp
    use efp_polar_module, only: allocate_Efield,calc_E_mp
#endif
    use integralstore_module, only: integralstore_kin_and_nuc_to_mem
    use density_data_module, only: density_data_alloc
    use occupation_module, only: alloc_occ_num
    use overlap_module, only: read_overlap
    implicit none
    ! *** end of interface ***

    type(fit) :: n_fit

    integer(i4_kind) :: rank

    rank = comm_rank()

    !
    ! Moved from main_scf():
    !
    call trace("PRESCF: reset_IO_files")
    !
    ! This includes files for debugging which can be set according to the
    ! users needs
    !
    if (rank == 0) &
         call reset_IO_files()  ! no comm!

    ! (1) calculate the electrostatic interaction energy
    !     between all nuclei (later: also point charges)
    if (rank == 0) &
         call core_interaction_calc() ! no comm!

    if(operations_solvation_effect) then
       call trace("PRESCF: const_contr_to_ham_and_energy_calc")
       !
       ! This broadcasts some data:
       !
       call const_part_of_ham_and_energy() ! does comm!
    endif

    if(calc_Pol_centers() .and. rank == 0) then
       call trace("PRESCF: get_field_nuc")
       call get_field_nuc() ! no comm!
#ifdef WITH_EFP
       if(n_efp > 1) then
          call allocate_Efield() ! no comm!
          call calc_E_mp() ! no comm?
       end if
#endif
    end if

    ! (2)  calculate one center  charge fit  integrals if  not already
    ! done by integral part
    if (.not. operations_integral) then
       call trace("PRESCF: fit_coeff_calc_chargefitnorm")
       if (rank == 0) &
            call fit_coeff_calc_chargefitnorm() ! no comm!
    endif

    ! (2a) read  in charge overlap matrix, and  if perturbation theory
    ! is  not turned  on or  if the  model density  approach  is used,
    ! decompose F and compute F^-1 Q
    call trace("PRESCF: read_mat_charge")
    if (rank == 0) &
         call read_mat_charge() ! no comm!

    ! (3) if  a calculation  with the ParaGAU  integral part is  to be
    ! done,   the  fitcoefficients   are   initialized  appropriately.
    ! Otherwise they are read in from the old lcgto.
    call trace("PRESCF: fitcoefficients are initialized")
    if (rank == 0) &
         call fit_coeff_initialize() ! no comm!

    call trace("PRESCF: fitcoefficients normalized")
    if (rank == 0) &
         call fit_coeff_normalize(spin_coeff=.true.) ! no comm!

    !
    ! (4a)  allocate  Fock   matrix  before  SCF.   See  corresponding
    ! hamiltonian_shutdown() in prescf_finalize():
    !
    call hamiltonian_setup(ssym)

    !
    ! (4b) if integrals are stored  on hard disk, transfer the kinetic
    ! and nuclear part of Hamiltonian to memory
    if (options_integrals_on_file()) then
       if (rank == 0) &
            call integralstore_kin_and_nuc_to_mem() ! no comm!
    endif

    !
    ! Read  overlap  from  disk  once,  before SCF  cycle.   This  was
    ! previousely done  every SCF cycle in  eigen_data_module. Then in
    ! main_scf  from a  master-only context,  now  it is  done from  a
    ! parallel context:
    !
    call read_overlap()

    !
    ! Allocate eigenvectors on all workers:
    !
    call eigen_data_alloc() ! no comm!

    ! Allocate,(a) density matrix, (b,master) occupation numbers
    call trace("PRESCF: allocate densmat and occupations")

    call density_data_alloc() ! no comm!

    if ( rank == 0 ) then
      ! occ_num seems to be present only on master?
      call alloc_occ_num() ! no comm!
    endif

    !
    ! (5) setup bounds on all workers:
    !
    call get_fit(n_fit) ! no comm!

    !
    ! This one is idempotent:
    !
    call bounds_calc(n_fit) ! no comm!
  end subroutine prescf_init

  subroutine prescf_finalize()
    !
    ! Purpose:  (conditionally)  release  densmat array  allocated  in
    ! prescf_init().  Pre-SCF  is a misnomer actually.   The prefix is
    ! to    keep    the    symmetry    between    prescf_init()    and
    ! prescf_finalize().   In  fact,  this  particular  subroutine  is
    ! called after SCF has converged.
    !
    ! Executed in parallel context.
    !
    ! author : Uwe Birkenheuer
    ! date   : 7/98
    !---------------------------------------------------------------
    ! Modifications
    !---------------------------------------------------------------
    ! Modification (Please copy before editing)
    ! Author:
    ! Date:
    ! Description:
    !--------------------------------------------------------------
    ! MODULES
    !--------------------------------------------------------------
    use density_data_module, only: density_data_free
    use operations_module, only: operations_potential
#ifdef WITH_MOLMECH
    use operations_module, only: operations_qm_mm_new
    use qmmm_interface_module, only: qm_mm
#endif
    use output_module, only: output_main_scf
    use iounitadmin_module, only: write_to_output_units
    use hamiltonian_module, only: hamiltonian_shutdown
#ifdef WITH_DFTPU
    use dft_plus_u_module, only: dft_plus_u_in_use, dft_plus_u_mo_in_use
#endif
    use overlap_module, only: dealloc_overlap
#ifdef WITH_BGY3D_NON_GPL
    use bgy3d, only: bgy3d_finalize
#endif
    implicit none
    ! *** end of interface ***

    logical :: do_post_dens

#ifdef WITH_DFTPU
    !
    ! Deallocate  overlap.   This  was  previousely   done  every  SCF
    ! iteration  after  diagonalization  in eigen_data_module.   Later
    ! from  main_scf() on  master. Slaves  always did  that  every SCF
    ! iteration from eigen_data_module.
    !
    if (.not. (dft_plus_u_in_use .or. dft_plus_u_mo_in_use)) then
       !
#endif
       call dealloc_overlap()
#ifdef WITH_DFTPU
    else
       ! In this case it  is done in dft_plus_u_grad_finalize() called
       ! from main_gradient()  after DFT+U contribution  to the forces
       ! have  been  computed. Deallocation  is  important  as it  may
       ! happen that the overlap will not be changed from one geometry
       ! iteration  into   another.   One   reason  for  it   is  that
       ! read_overlap() does nothing if it thinks overlap is valid. At
       ! some point,  the latest is  after completion of  SCF, current
       ! overlap should be invalidated.
    endif
#endif

    !
    ! Deallocate Hamiltonian (after output_eigendata). Deallocation of
    ! 'ham_tot'    and     'ham_lsft'    has    been     moved    into
    ! hamiltonian_shutdown().   See  corresponding hamiltonian_setup()
    ! in prescf_init().
    !
    if(output_main_scf) call write_to_output_units &
         ("MAIN_SCF: hamiltonian_shutdown")
    call hamiltonian_shutdown()

    do_post_dens = .true. ! do deallocate density matrix
    do_post_dens = do_post_dens .and. .not. operations_potential
#ifdef WITH_MOLMECH
    do_post_dens = do_post_dens .and. .not.(operations_qm_mm_new .and. qm_mm)
#endif
    if (do_post_dens) then
      ! deallocate density matrix
      call density_data_free()
    endif

#ifdef WITH_BGY3D_NON_GPL
    call bgy3d_finalize ()
#endif
  end subroutine prescf_finalize

  subroutine trace(msg)
    use output_module, only: output_main_scf
    use iounitadmin_module, only: write_to_output_units
    use comm, only: comm_rank
    implicit none
    character(len=*), intent(in) :: msg
    ! *** end of interface ***

    if(output_main_scf .and. comm_rank() == 0) &
      call write_to_output_units(msg)
  end subroutine trace
  !--------------- End of module ----------------------------------
end module prescf_module

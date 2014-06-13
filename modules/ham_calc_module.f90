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
module ham_calc_module
  !
  !  Purpose: Contains routines for calculating the
  !           'three-center'-parts of the hamiltonian
  !           as well as routines to add already existing
  !           integrals to the hamiltonian. Also the corresponding
  !           energies are calculated by different calls
  !           The Exchange part of the hamiltonian can
  !           either be done via the Fitfunctions or
  !           numerically. Only in the former case
  !           the contribution to and the energy are calculated here.
  !           The numerical calculation of the XC part takes place in
  !           'xc_hamiltonian'.
  !           In this module are routines for the master
  !           as well as for the slave. Each routine includes
  !           two lines to check they are called in the right
  !           way.
  !
  !
  !  Module called by: mainscf
  !
  !  References: still none
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification
  ! Author: TG
  ! Date:   9-10/96
  ! Description: 1. Subroutine build_hamiltonian of file
  !                 ../build_hamiltonian.f90 was integrated
  !                 in this module.
  !              2. Subroutine ham_calc_slave was superseded
  !                 by density_fit_coulomb_interaction
  !              3. Routines mentioned in point 1 and two
  !                 now sum up the hamiltonian part in a
  !                 cascadic algorithm.
  !              4. Spin dependence of ham_kin, ham_nuc and
  !                 ham_coul was removed.
  !
  ! Modification (Please copy before editing)
  ! Author: Uwe Birkenheuer
  ! Date:   14/7/97
  ! Description: 1. Control variable OPTIONS_XCMODE() introduced
  !              2. Control variable OPTIONS_ETOTMODE() introduced
  !              3. Model density Hamiltonian introduced
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   8/97
  ! Description: Calculation of bounds moved to prescf
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   12/99
  ! Description: extension to solvent effect
  !
  ! Modification
  ! Author: TMS
  ! Date:   07/11
  ! Description: larger restructurations, unification of the three
  !              OPTIONS_ETOTMODE paths and extension to exact exchange
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
  use datatype    ! user defined types
  use readwriteblocked_module, only: readwriteblocked_tapehandle
  use options_module, only: xcmode_numeric_exch, xcmode_exchange_fit           &
                          , xcmode_extended_mda, xcmode_model_density          &
                          , recover_scfstate                 &
                          , options_spin_orbit, options_spin_restricted        &
                          , options_integrals_on_file, options_xcmode          &
                          , options_recover
  implicit none
  private
  save

  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------

  public :: ham_calc_main!(loop)

  !------------ public formal parameters types to subroutines ---
  public :: arrmat3, arrmat2, readwriteblocked_tapehandle

!================================================================
! End of public interface of module
!================================================================

  !------------ Declaration of constants  -----------------------

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  subroutine ham_calc_main (iloop)
    !
    ! Entry point for building hamiltonian (Fock) matrix.  Executed in
    ! parallel  context by  all  workers. Called  from main_scf()  and
    ! man_slave() on master and slaves, respectively.
    !
    use comm, only: comm_rank, comm_bcast
    use comm_module, only: comm_init_send, comm_send, &
         comm_all_other_hosts                     ! MUSTDIE
    use msgtag_module, only: msgtag_ham_calc_main ! MUSTDIE
    use symmetry_data_module, only: ssym
    use unique_atom_module, only: unique_atoms
    use hamiltonian_module, only: reset_ham
    use energy_calc_module, only: init_energy
    use bounds_module, only: bounds_ch, bounds_xc
    use hamiltonian_module, only: ham_tot, ham_tot_real, ham_tot_imag
    use density_data_module, only: densmat, densmat_real, densmat_imag
#ifdef WITH_BGY3D
    use bgy3d, only: bgy3d_term
#endif
    implicit none
    integer (i4_kind), value :: iloop ! meaningfull on master only
    ! *** end of interface ***

    integer(i4_kind) :: rank
    integer(i4_kind) :: loop    ! will be valid on all workers

#ifdef WITH_BGY3D
    real(r8_kind) :: e_bgy      ! electrostatic energy
#endif

    rank = comm_rank()

    !
    ! Tell slaves to call ham_calc_main. FIXME: yet better call from a
    ! parallel context ...
    !
    if (rank == 0) then
      ! FIXME: legacy communication primitives here:
      call comm_init_send (comm_all_other_hosts, msgtag_ham_calc_main)
      ! ... the corresponding tag receive in main_slave
      call comm_send()
    endif

    ! Broadcast local copy to slaves, iloop is intent(in):
    loop = iloop
    call comm_bcast (loop)

    ! Dont use iloop below!
    iloop = -1

    ! RESET_HAM allocates and initializes the necessary parts
    ! of the hamiltonian:
    !    ham_tot  (spin polarized if required)
    call reset_ham(ssym)

    call init_energy() ! initializes the energies

    !
    ! For historical reasons much of the code is executed by master only (see
    ! build_1e_hamiltonian). The Coulomb and exchange terms are computed in
    ! parallel, however:
    !
    if (.not. options_spin_orbit) then
      !
      ! STANDARD SCF
      !
      call build_2e_hamiltonian (loop, bounds_ch, bounds_xc                    &
                               , h_matrix       = ham_tot                      &
                               , d_matrix       = densmat                      )
    else
      !
      ! SPIN ORBIT
      !
      call build_2e_hamiltonian (loop, bounds_ch, bounds_xc                    &
                               , h_matrix_real = ham_tot_real                  &
                               , h_matrix_imag = ham_tot_imag                  &
                               , d_matrix_real = densmat_real                  &
                               , d_matrix_imag = densmat_imag                  )
    endif

    if ( rank == 0 ) then
      if (.not. options_spin_orbit) then
        !
        ! STANDARD SCF
        !
        call build_1e_hamiltonian( loop, ham_tot      = ham_tot                &
                                       , densmat      = densmat                )
        !
      else
        !
        ! SPIN ORBIT
        !
        call build_1e_hamiltonian( loop, ham_tot_real = ham_tot_real           &
                                       , ham_tot_imag = ham_tot_imag           &
                                       , densmat_real = densmat_real           &
                                       , densmat_imag = densmat_imag           )
      endif
    endif

#ifdef WITH_BGY3D
    if (.not. options_spin_orbit) then
        !
        ! STANDARD SCF
        !
       ASSERT(allocated(ham_tot))
       ASSERT(allocated(densmat))

       ! Does nothing  for negative arguments, runs  pure solvent when
       ! argument  is zero  and updates  the solvent  distribution for
       ! positive arguments.   FIXME: if SCF converges  early the call
       ! may have no effect at all. The instance of the hamiltonian on
       ! master is incremented by  the solvation term. Other instances
       ! are not  affected. Energy  is effectively all-reduced  and is
       ! valid an all workers.
       call bgy3d_term (loop - 10 - 1, unique_atoms, densmat, ham_tot, e_bgy)
    else
       !
       ! SPIN ORBIT
       !
       ! Do   nothing.   FIXME:   this  will   silently   ignores  BGY
       ! option. This is  a workaround to make PG  regresison tests of
       ! SO pass for builds WITH_BGY3D.
    endif
#endif

    !
    ! Deallocation of hamiltonian that slaves did at this place before
    ! is  now done  together with  master  from hamiltonian_shutdown()
    ! which is itself called from prescf_finalize().
    !
    ! NOTE: the  storage in ham_tot(:) %  m(:, :, :) may  be used e.g.
    !       for receiving the irrep blocks to diagonalize.
    !
  end subroutine ham_calc_main

  !*****************************************************************************
  subroutine build_1e_hamiltonian( loop                                        &
                                 , ham_tot, ham_tot_real, ham_tot_imag         &
                                 , densmat, densmat_real, densmat_imag         )
    !
    !  Purpose: This routine assembles the Hamiltonian
    !           for the scf-cycles. Nuclear and kinetic part of
    !           the hamiltonian are provided by the routine
    !           "prescf" located in the prescf_module.
    !           Here only the coulomb- and  exchange-
    !           part are calculated.
    !           Furtheron, the exchange part is calculated ONLY
    !           if we are working with fitcoefficients for this
    !           part. If the appropriate flag NUM_EXCH is set
    !           the variable xc-term will be calcuated numerically
    !           by the routine 'xc_hamiltonian'.
    !
    !           The previously used 'options_energies_always' has
    !           been subsituted by 'options_etotmode', whose three
    !           different paths have now unified to a single one.
    !           From now on, 'options_etotmode' returns always
    !           'etotmode_etot_only' and all matrix terms are computed
    !           (or retrieved from memory / file). The energy
    !           contribution of these terms is computed, dumped inside
    !           energy_calc_module. The matrices are added directly to
    !           ham_tot and forgotten afterwards
    !
    !           The next step is the call to the eigensolving
    !           routines (eigen_data_module)
    !
    !  Subroutine called by: main_scf
    !
    !  Author: Folke Noertemann
    !  Date: 7/95
    !
    use print_module   !interface for outroutines
    use init_module,              only: init
    use symmetry_data_module,     only: ssym
    use output_module,            only: output_hamiltonian
    use xc_cntrl,                 only: xc_is_on=>is_on                        &
                                      , xc_ANY
    use xc_hamiltonian,           only: ham_xc_arr                             &
                                      , mat_initialized                        &
                                      , ham_xc_arr_real                        &
                                      , ham_xc_arr_imag
    use hamiltonian_module,       only: print_hamiltonian
    use iounitadmin_module,       only: output_unit                            &
                                      , stdout_unit
    use energy_calc_module,       only: set_energy                             &
                                      , write_energy
    use occupation_module,        only: get_diagonal_offset
    use efield_module,            only: efield_applied
    use solv_electrostat_module,  only: ham_solv_el
    use solv_cavity_module,       only: sol_start_cycle
    use operations_module,        only: operations_solvation_effect
    use overlap_module,           only: overlap
#ifdef WITH_DFTPU
    use dft_plus_u_module,        only: dft_plus_u_in_use                      &
                                      , dft_plus_u_term
#endif
    use induced_dipoles_module,   only: calc_Pol_centers                       &
                                      , ham_Pol
    use integralstore_module,     only: integralstore_2cob_kin                 &
                                      , integralstore_2cob_nuc
    use fit_coeff_module,         only: sp_initialized                         &
                                      , xc_initialized

    implicit none
    !------------ Declaration of formal parameters ----------------
    !
    integer(i4_kind),        intent(in)    :: loop ! actual SCF-cycle number
                                                   ! within the current SCF run
    !
    type(arrmat3), optional, intent(inout) :: ham_tot(:)
    type(arrmat2), optional, intent(inout) :: ham_tot_real(:)
    type(arrmat2), optional, intent(inout) :: ham_tot_imag(:)
    !
    type(arrmat3), optional, intent(in)    :: densmat(:)
    type(arrmat2), optional, intent(in)    :: densmat_real(:)
    type(arrmat2), optional, intent(in)    :: densmat_imag(:)
    !
    !** End of interface *****************************************

    !================================================================
    ! End of public interface of module
    !================================================================

    !------------ Declaration of subroutines ------------------------
    external error_handler
    !------------ Declaration of local variables --------------------
    integer(i4_kind)   :: counter, i, j, k, i_gamma, alloc_stat, n
    real(r8_kind) :: diagonal_offset
    logical :: numeric_exch, exchange_fit, model_density
    logical :: add_constant_diagonal

    integer(i4_kind)                     :: i_spin
    integer(i4_kind)                     :: n_irrep
    integer(i4_kind), allocatable        :: dim_irrep(:)
    ! n_irrep   : number of irreps
    ! dim_irrep : number of independent functions in irrep

    ! energies
    real(r8_kind) :: e_kin
    real(r8_kind) :: e_nuc
    real(r8_kind) :: e_dft_plus_u ! correction to the total energy due to DFT+U
    real(r8_kind) :: e_efield ! electric field contribution
    !------------------ Executable code ---------------------------
    !
    ! Check input variables
    if (.not. options_spin_orbit) then
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       ASSERT(present(ham_tot))
       ASSERT(present(densmat))
       !
    else ! options_spin_orbit
       !
       ! SPIN ORBIT
       !
       ASSERT(present(ham_tot_real))
       ASSERT(present(ham_tot_imag))
       ASSERT(present(densmat_real))
       ASSERT(present(densmat_imag))
       !
    endif ! options_spin_orbit
    !
    ! retrieve some initial information from modules
    numeric_exch  = options_xcmode() == xcmode_numeric_exch
    exchange_fit  = options_xcmode() == xcmode_exchange_fit
    model_density = options_xcmode() == xcmode_model_density .or. &
                    options_xcmode() == xcmode_extended_mda

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif ! options_spin_orbit

    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! KINETIC ENERGY TERM
    !
    ! Computing:          __
    !                    \  '       ^
    !             E   =   >   ( k | T | l )  P
    !              T     /__,                 kl
    !                    k,l
    !
    ! And incrementing the Hamiltonian matrix by
    !                          ^
    !             T    = ( k | T | l )
    !              kl
    ! Whereas:
    !
    !          k, l         :  basis functions
    !
    !          P            :  density matrix
    !           kl
    !               ^
    !         ( k | T | l ) :  two center kinetic energy integral
    !
    if (.not. options_spin_orbit) then
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call hamiltonian_term( integrals = integralstore_2cob_kin               &
                            , factor    = + 1.0_r8_kind                        &
                            , d_matrix  = densmat                              &
                            , h_matrix  = ham_tot                              &
                            , trace     = e_kin                                )
    else
       !
       ! SPIN ORBIT
       !
       call hamiltonian_term( integrals     = integralstore_2cob_kin           &
                            , factor        = + 1.0_r8_kind                    &
                            , d_matrix_real = densmat_real                     &
                            , d_matrix_imag = densmat_imag                     &
                            , h_matrix_real = ham_tot_real                     &
                            , h_matrix_imag = ham_tot_imag                     &
                            , trace         = e_kin                            )
    endif
    !
    call set_energy(kin = e_kin)
    !
    call write_energy(output_unit, print_kin = .TRUE.)
    call write_energy(stdout_unit, print_kin = .TRUE.)
    !
    ! end of KINETIC ENERGY TERM
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! NUCLEAR ATTRACTION TERM
    !
    ! Computing:            __
    !                      \  '       ^
    !             E   = -   >   ( k | Z | l )  P
    !              Z       /__,                 kl
    !                      k,l
    !
    ! And incrementing the Hamiltonian matrix by
    !                          ^
    !             Z    = ( k | Z | l )
    !              kl
    ! Whereas:
    !
    !          k, l         :  basis functions
    !
    !          P            :  density matrix
    !           kl
    !               ^
    !         ( k | Z | l ) :  two center nuclear attraction integral
    !
    if (.not. options_spin_orbit) then
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call hamiltonian_term( integrals = integralstore_2cob_nuc               &
                            , factor    = - 1.0_r8_kind                        &
                            , d_matrix  = densmat                              &
                            , h_matrix  = ham_tot                              &
                            , trace     = e_nuc                                )
    else
       !
       ! SPIN ORBIT
       !
       call hamiltonian_term( integrals     = integralstore_2cob_nuc           &
                            , factor        = - 1.0_r8_kind                    &
                            , d_matrix_real = densmat_real                     &
                            , d_matrix_imag = densmat_imag                     &
                            , h_matrix_real = ham_tot_real                     &
                            , h_matrix_imag = ham_tot_imag                     &
                            , trace         = e_nuc                            )
    endif
    !
    call set_energy(nuc = e_nuc)
    !
    call write_energy(output_unit, print_nuc = .TRUE.)
    call write_energy(stdout_unit, print_nuc = .TRUE.)
    !
    ! end of NUCLEAR ATTRACTION TERM
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    !
    ! Sum up contributions to hamiltonian:
    !
    irrep_loop : do i_gamma = 1, n_irrep
       if(operations_solvation_effect.and.loop >= sol_start_cycle+1) then
          ! NOTE: For UKS fill both axes directly
          do i_spin = 1, ssym%n_spin
             ham_tot(i_gamma)%m(:, :, i_spin) =                                &
                                            ham_tot(i_gamma)%m(:, :, i_spin)   &
                                          + ham_solv_el(i_gamma)%m
          enddo
       endif

       if(calc_Pol_centers() .and. loop > 1) then
          ! NOTE: For UKS fill both axes directly
          do i_spin = 1, ssym%n_spin
              ham_tot(i_gamma)%m(:, :, i_spin) =                               &
                                            ham_tot(i_gamma)%m(:, :, i_spin)   &
                                          + ham_Pol(i_gamma)%m
          enddo
       end if
    enddo irrep_loop

    if (numeric_exch) then
       if (loop /= 1 .or. options_recover() == recover_scfstate) then

          ![[=== Add XC contribution to the total hamiltonian ===
          if( xc_is_on(xc_ANY) )then

          ! get ham_xc_arr from hamiltonian_module as it
          ! comes from the xc_hamiltonian routine (lower triangle only!!)
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             counter = 0
             do i = 1, n_irrep
                do j = 1, dim_irrep(i)
                   do k = 1, j-1
                      counter = counter + 1
                      ! only H_ij with i <= j loaded here, real and imag parts
                      ham_tot_real(i)%m(k, j) = ham_tot_real(i)%m(k, j)        &
                                              + ham_xc_arr_real(counter)
                      ham_tot_real(i)%m(j, k) = ham_tot_real(i)%m(j, k)        &
                                              + ham_xc_arr_real(counter)
                      ham_tot_imag(i)%m(k, j) = ham_tot_imag(i)%m(k, j)        &
                                              + ham_xc_arr_imag(counter)
                      ham_tot_imag(i)%m(j, k) = ham_tot_imag(i)%m(j, k)        &
                                              - ham_xc_arr_imag(counter)
                   enddo!loop over nu
                   counter = counter + 1
                   !
                   ham_tot_real(i)%m(j, j) = ham_tot_real(i)%m(j, j)           &
                                           + ham_xc_arr_real(counter)
                   ham_tot_imag(i)%m(j, j) = ham_tot_imag(i)%m(j, j)           &
                                           + ham_xc_arr_imag(counter)
                enddo!loop over mu
             enddo! loop over irreps
          else ! options_spin_orbit
             !
             ! STANDARD SCF (NO SPIN ORBIT)
             !
             counter = 0
             do i = 1, n_irrep
                do j = 1, dim_irrep(i)
                   do k = 1, j-1
                      do i_spin = 1, ssym%n_spin
                         counter = counter + 1
                         ham_tot(i)%m(k,j,i_spin) = ham_tot(i)%m(k,j,i_spin)   &
                                                  + ham_xc_arr(counter)
                         ham_tot(i)%m(j,k,i_spin) = ham_tot(i)%m(j,k,i_spin)   &
                                                  + ham_xc_arr(counter)
                      enddo ! loop over spin
                   enddo!loop over nu
                   do i_spin = 1, ssym%n_spin
                      counter = counter + 1
                      ham_tot(i)%m(j,j,i_spin) = ham_tot(i)%m(j,j,i_spin)      &
                                               + ham_xc_arr(counter)
                   enddo ! loop over spin
                enddo!loop over mu
             enddo! loop over irreps
          endif ! options_spin_orbit
          endif
          !]]=== eof Add XC contribution to the total hamiltonian ===

       endif
    endif ! numeric_exch

    ! ----------------------------------------------------

    ! Add part of hamiltonian from efield directly to ham tot
    ! this is done here to allow the efield also to be applied
    ! in reverse directions for opposit spins
    ! the energy contribution of efield is also calculated
    if ( efield_applied() ) then
       !
       ! FIXME: maybe remove "integrals_on_file" case for electric field
       !        and use generic hamiltonian_term(...) instead:
       !
       call electric_field_term(densmat, ham_tot, e_efield)

       call set_energy(efield=e_efield)
    endif

#ifdef WITH_DFTPU
    ! ----------------------------------------------------
    ! DFT+U routine is called here. This should not called
    ! at an earlier point because the HAM(alpha) = HAM(beta)
    ! are kept same. See above.
    !------------------------------------------------------
    if( dft_plus_u_in_use )then
      !
      ! Add DFT+U part of the Hamiltonian,
      ! compute the energy correction:
      ASSERT(.not.options_spin_orbit)

      ! load overlap for checking and some DFT+U expressions:
      ! MOVED BEFORE SCF LOOP: call read_overlap()
      ASSERT(allocated(overlap))

      call dft_plus_u_term(overlap, densmat, ham_tot, e_dft_plus_u)

      ! copy the value into energy_calc_module:
      call set_energy(dft_plus_u=e_dft_plus_u)
    endif
#endif

    if (ssym%n_spin > 1) then
       ! if this is a spin-unrestricted calculation we may
       ! have to provide an initial spin difference here.
       call get_diagonal_offset(diagonal_offset)
       if (diagonal_offset /= 0.0_r8_kind) then
          ! The necessity of providing an initial! spin difference
          ! is now controlled by the <xx>_initialized flags (UB 7/97)
          if (numeric_exch) then
             add_constant_diagonal = mat_initialized
          elseif (exchange_fit) then
             add_constant_diagonal = xc_initialized
          elseif (model_density) then
             add_constant_diagonal = sp_initialized
          else
             ! make  compiler happy,  define  add_constant_diagonal in
             ! all cases:
             add_constant_diagonal = .false.
             ABORT("no such case")
          endif
          if (add_constant_diagonal) then
             do k=1,dim_irrep(1)
                ham_tot(1)%m(k,k,2)=ham_tot(1)%m(k,k,2)+diagonal_offset
             enddo!loop over nu
          endif
          if(model_density) sp_initialized = .false.
       endif
    endif
    ! ----------------------------------------------------

    if (numeric_exch) then
       ! Expand upper triangle, now ham_tot is accumulated completely
       ! check if needed!!!!!!
       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          do i=1,n_irrep
             do k=2,dim_irrep(i)
                ! Now H_ij with i > j is loaded
                ham_tot_real(i)%m(k,:k-1) = ham_tot_real(i)%m(:k-1,k)
                ham_tot_imag(i)%m(k,:k-1) = - ham_tot_imag(i)%m(:k-1,k)
             enddo
          enddo
       else ! options_spin_orbit
          !
          ! STANDARD SCF (NO SPIN ORBIT)
          !
          do i=1,n_irrep
             do i_spin = 1, ssym%n_spin
                do k=2,dim_irrep(i)
                   ! Now H_ij with i > j is loaded
                   ham_tot(i)%m(k,:k-1,i_spin) = ham_tot(i)%m(:k-1,k,i_spin)
                enddo
             enddo
          enddo
       endif ! options_spin_orbit
    endif

    ! if 'options_xcmode() /= xcmode_split_ham', print_hamiltonian simply does
    ! not print the unknown parts. (ham_xc_arr must be passed as optional
    ! parameter here because explicit reference by a use statement within
    ! print_hamiltonian would result in recursive module dependencies.)
    if (output_hamiltonian) then
       if (numeric_exch .and. &
            (loop /= 1 .or. options_recover() == recover_scfstate))then

          if(options_spin_orbit)then
             call print_hamiltonian("full","xc","tot", &
                  ham_xc_arr     =ham_xc_arr_real,&
                  ham_xc_arr_imag=ham_xc_arr_imag)
          else
             call print_hamiltonian("full","xc","tot", &
                  ham_xc_arr=ham_xc_arr)
          endif
       else
          call print_hamiltonian("full","xc","tot")
       endif
    endif

    deallocate(dim_irrep,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
  end subroutine build_1e_hamiltonian
  !*****************************************************************************

  !*****************************************************************************
  subroutine build_2e_hamiltonian (loop, bounds_ch, bounds_xc                  &
                                 , d_matrix                                    &
                                 , d_matrix_real, d_matrix_imag                &
                                 , h_matrix                                    &
                                 , h_matrix_real, h_matrix_imag                )
    !
    ! Purpose: calculation of the Coulomb matrix in parallel context
    !
    ! Input-Parameters
    !
    !    bounds_ch(:)            - number of chargefit coefficients per proc
    !    bounds_xc(:)            - number of exch.-fit coefficients per proc
    !    d_matrix(:)%m(:,:,:)    - density matrix standard scf
    !    d_matrix_real(:)%m(:,:) - density matrix spin orbit
    !    d_matrix_imag(:)%m(:,:) - density matrix spin orbit
    !    h_matrix(:)%m(:,:,:)    - hamiltonian matrix standard scf
    !    h_matrix_real(:)%m(:,:) - hamiltonian matrix spin orbit
    !    h_matrix_imag(:)%m(:,:) - hamiltonian matrix spin orbit
    !    trace                   - expectation value
    !
    !------------ Modules ----------------------------------------
    use comm, only: comm_reduce, comm_rank, comm_bcast, comm_same
    use iounitadmin_module,   only: output_unit                                &
                                  , stdout_unit
    use spin_orbit_module,    only: whatis                                     &
                                  , op_BackTrafo
    use dimensions,           only: IrrBasDim                                  &
                                  , IrrBasDimSpor
    use symmetry_data_module, only: symmetry_data_n_spin
    use fit_coeff_module, only: fit_coeff_n_ch, coeff_charge, coeff_charge_veff
    use energy_calc_module,   only: energ_coul_2z                              &
                                  , set_energy                                 &
                                  , write_energy                               &
                                  , direct_2c_energy_calc_and_add
    use eri4c_options, only: J_exact
#if WITH_ERI4C == 1
    use eri4c_options, only: K_exact
#endif
    use options_module,       only: options_spin_orbit
    use unique_atom_module,   only: unique_atoms
    use overlap_module,       only: overlap
#ifdef WITH_DFTPU
    use dft_plus_u_module,    only: dft_plus_u_mo_in_use                       &
                                  , dft_plus_u_proj
#endif
    !
    implicit none
    !
    !------------ Declaration of formal parameters ---------------
    !
    integer (i4_kind), intent(in)          :: loop
    integer(i4_kind), intent(in)           :: bounds_ch(:)
    integer(i4_kind), intent(in)           :: bounds_xc(:)
    !
    type(arrmat3), optional, intent(in)    :: d_matrix(:)
    type(arrmat2), optional, intent(in)    :: d_matrix_real(:)
    type(arrmat2), optional, intent(in)    :: d_matrix_imag(:)
    !
    type(arrmat3), optional, intent(inout) :: h_matrix(:)
    type(arrmat2), optional, intent(inout) :: h_matrix_real(:)
    type(arrmat2), optional, intent(inout) :: h_matrix_imag(:)
    !
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    logical                                :: numeric_exch                     &
                                            , exchange_fit                     &
                                            , model_density
    !
    integer(i4_kind)                       :: pid
    !
    integer(i4_kind)                       :: i_gamma
    integer(i4_kind)                       :: i_spin
    integer(i4_kind)                       :: i_start
    integer(i4_kind)                       :: i_stop
    integer(i4_kind)                       :: i_start_m
    integer(i4_kind)                       :: i_stop_m
    integer(i4_kind)                       :: n_triag
    !
    integer(i4_kind)                       :: n_spin
    integer(i4_kind)                       :: n_irrep
    integer(i4_kind)                       :: n_fit
    !
    integer(i4_kind)                       :: dim_irrep
    integer(i4_kind)                       :: dim_2c
    !
    ! energies
    real(r8_kind)                          :: e_coul
    real(r8_kind)                          :: e_2z_coul
    !
    real(r8_kind)                          :: e_exex
    real(r8_kind)                          :: e_dft_plus_u
    real(r8_kind)                          :: e_fit_xc
    !
    ! temporary arrays for all nodes...
    real(r8_kind), allocatable             :: j_array(:)
    real(r8_kind), allocatable             :: m_array(:)
    real(r8_kind), allocatable             :: f_array(:)
    real(r8_kind), allocatable             :: x_array(:)
    !
    !------------ Executable code --------------------------------
    !
    ! Check input variables
    if (.not. options_spin_orbit) then
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       ASSERT(present(h_matrix))
       ASSERT(present(d_matrix))
       !
    else ! options_spin_orbit
       !
       ! SPIN ORBIT
       !
       ASSERT(present(h_matrix_real))
       ASSERT(present(h_matrix_imag))
       ASSERT(present(d_matrix_real))
       ASSERT(present(d_matrix_imag))
       !
    endif ! options_spin_orbit
    !
    ! base-1 process ID:
    pid = 1 + comm_rank()

    ! See   ham_calc_main(),  loop   has   the  same   value  on   all
    ! workers. FIXME: rm when confident!
    ASSERT(comm_same(loop))

    e_coul    = 0.0_r8_kind
    e_2z_coul = 0.0_r8_kind
    e_exex    = 0.0_r8_kind
    e_fit_xc  = 0.0_r8_kind
    !
    numeric_exch  = options_xcmode() == xcmode_numeric_exch
    exchange_fit  = options_xcmode() == xcmode_exchange_fit
    model_density = options_xcmode() == xcmode_model_density .or. &
                    options_xcmode() == xcmode_extended_mda
    !
    n_spin = symmetry_data_n_spin()
    !
    if (.not. options_spin_orbit) then
       n_irrep = size(IrrBasDim)
       dim_2c  = sum(IrrBasDim * (IrrBasDim + 1) / 2)
    else
       n_irrep = size(IrrBasDimSpor)
       dim_2c  = sum(IrrBasDimspor * (IrrBasDimspor + 1) / 2)
    endif
    !
    ! Allocate tmp linear arrays
    if (.not. options_spin_orbit) then
       allocate(j_array( dim_2c ))
       j_array = 0.0_r8_kind
       if (model_density) then
          allocate(m_array( dim_2c * n_spin ))
          m_array = 0.0_r8_kind
        endif
       if (exchange_fit) then
          allocate(f_array( dim_2c )         )
          f_array = 0.0_r8_kind
          allocate(x_array( dim_2c * n_spin ))
          x_array = 0.0_r8_kind
        endif
    else
       allocate(j_array( dim_2c * 2 ))
    endif
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
#ifdef WITH_ERI4C
    if (.not. J_exact ) then
#endif
      !
      ! DENSITY FIT COULOMB TERM
      !
      ! Formation of the Coulomb matrix from 3-center integrals (exact Coulomb).
      ! Computing:                __
      !                          \  '
      !             J    =        >   c    ( i j | f  )
      !              ij          /__,  k            k
      !                           k
      !
      ! as well as the energies   __
      !                          \  '
      !             E    =        >   P    J
      !              J           /__,  ij   ij
      !                          i j
      !
      ! and the addional two center component resulting from density fit:
      !                           __
      !                       1  \  '
      !             E    = - ---  >   c   ( f   | f  )  c
      !              2c       2  /__,  k     k     l     l
      !                          k,l
      ! Whereas:
      !          i            :  KS orbitals
      !
      !          f            :  auxiliary fit basis functions
      !           k
      !
      !          c            :  charge fit coefficients
      !           k
      !
      !          P            :  density matrix elements
      !           ij
      !
      !        ( i j | f  )   :  three-center electron repulsion integral
      !                 k
      !
      !        ( f  | f  )    :  two center electron repulsion integral
      !           k    l
      !
      if ( comm_rank() == 0 ) then
        n_fit = size(coeff_charge)
ASSERT(n_fit==fit_coeff_n_ch())
        if(whatis(op_BackTrafo).eq.2)then
           call energ_coul_2z(coeff_charge + coeff_charge_veff, e_2z_coul)
        else
           call energ_coul_2z(coeff_charge, e_2z_coul) ! E_coul(fit only)
        endif
        !
        call set_energy(c_2z = e_2z_coul)
        !
      endif
      !
      ! initialize a working array on all procs
      if (.not. options_spin_orbit) then
        !
        if( .not. model_density .and. .not. exchange_fit) then
          !
          ! STANDARD SCF
          !
          call density_fit_coulomb_interaction( bounds_ch = bounds_ch          &
                                              , bounds_xc = bounds_xc          &
                                              , j_array   = j_array            )
          !
          call comm_reduce(j_array)
          !
          if ( comm_rank() == 0 ) then
            !
            i_start = 1
            i_stop  = 0
            do i_gamma = 1, n_irrep
              ! Address logic
              dim_irrep = IrrBasDim(i_gamma)
              n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
              i_stop    = i_stop + n_triag
              do i_spin = 1, n_spin
                call direct_2c_energy_calc_and_add(                            &
                                   f_t_arr    = 1.0_r8_kind                    &
                                 , f_trace    = 1.0_r8_kind                    &
                                 , t_arr      = j_array(i_start:i_stop)        &
                                 , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)&
                                 , e_2c       = e_coul                         &
                                 , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin))
              enddo
              i_start = i_start + n_triag
            enddo
          endif
          !
        elseif(model_density) then
          !
          ! MODEL DENSIY APPROACH
          !
          write(*,*) 'WARNING! model density approach is not verified'
          !
          call density_fit_coulomb_interaction( bounds_ch    = bounds_ch       &
                                              , bounds_xc    = bounds_xc       &
                                              , j_array      = j_array         &
                                              , m_array      = m_array         )
          !
          call comm_reduce(j_array)
          call comm_reduce(m_array)
          !
          if ( comm_rank() == 0 ) then
            !
            i_start   = 1
            i_stop    = 0
            i_start_m = 1
            i_stop_m  = 0
            do i_gamma = 1, n_irrep
              ! Address logic
              dim_irrep = IrrBasDim(i_gamma)
              n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
              i_stop    = i_stop   + n_triag
              i_stop_m  = i_stop_m + n_triag * n_spin
              do i_spin = 1, n_spin
                !
                call direct_2c_energy_calc_and_add(                            &
                                   f_t_arr    = 1.0_r8_kind                    &
                                 , f_trace    = 1.0_r8_kind                    &
                                 , t_arr      = j_array(i_start:i_stop)        &
                                 , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)&
                                 , e_2c       = e_coul                         &
                                 , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin))
                ! no contribution to energy
                call direct_2c_energy_calc_and_add(                            &
                               f_t_arr    = 1.0_r8_kind                        &
                             , f_trace    = 0.0_r8_kind                        &
                             , t_arr      = m_array(i_start_m:i_stop_m:n_spin) &
                             , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)    &
                             , e_2c       = e_coul                             &
                             , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin)    )
              enddo
              i_start   = i_start   + n_triag
              i_start_m = i_start_m + n_triag * n_spin
            enddo
          endif
          !
        elseif(exchange_fit) then
          !
          ! EXCHANGE FIT
          !
          write(*,*) 'WARNING! exchange fit path is not verified'
          !
          call density_fit_coulomb_interaction( bounds_ch    = bounds_ch       &
                                              , bounds_xc    = bounds_xc       &
                                              , j_array      = j_array         &
                                              , f_array      = f_array         &
                                              , x_array      = x_array         )
          !
          call comm_reduce(j_array)
          call comm_reduce(f_array)
          call comm_reduce(x_array)
          !
          if ( comm_rank() == 0 ) then
            !
            i_start   = 1
            i_stop    = 0
            i_start_m = 1
            i_stop_m  = 0
            do i_gamma = 1, n_irrep
              ! Address logic
              dim_irrep = IrrBasDim(i_gamma)
              n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
              i_stop    = i_stop   + n_triag
              i_stop_m  = i_stop_m + n_triag * n_spin
              do i_spin = 1, n_spin
                !
                call direct_2c_energy_calc_and_add(                            &
                                   f_t_arr    = 1.0_r8_kind                    &
                                 , f_trace    = 1.0_r8_kind                    &
                                 , t_arr      = j_array(i_start:i_stop)        &
                                 , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)&
                                 , e_2c       = e_coul                         &
                                 , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin))
                ! no contribution to energy
                call direct_2c_energy_calc_and_add(                            &
                               f_t_arr    = - 1.0_r8_kind                      &
                             , f_trace    = 0.0_r8_kind                        &
                             , t_arr      = x_array(i_start_m:i_stop_m:n_spin) &
                             , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)    &
                             , e_2c       = e_fit_xc                           &
                             , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin)    )
                ! only for energy
                call direct_2c_energy_calc_and_add(                            &
                               f_t_arr    = 0.0_r8_kind                        &
                             , f_trace    = - 1.0_r8_kind                      &
                             , t_arr      = f_array(i_start:i_stop)            &
                             , d_mat_real = d_matrix(i_gamma)%m(:,:,i_spin)    &
                             , e_2c       = e_fit_xc                           &
                             , h_mat_real = h_matrix(i_gamma)%m(:,:,i_spin)    )
              enddo
              i_start   = i_start   + n_triag
              i_start_m = i_start_m + n_triag * n_spin
            enddo
            !
            call set_energy(xc = e_fit_xc)
            !
          endif
          !
        endif
        !
      else
        !
        ! SPIN ORBIT
        !
        ASSERT(.not. model_density)
        if ( comm_rank() == 0 ) then
          ASSERT(present(h_matrix_real))
          ASSERT(present(h_matrix_imag))
          ASSERT(present(d_matrix_real))
          ASSERT(present(d_matrix_imag))
        endif
        !
        call density_fit_coulomb_interaction( bounds_ch     = bounds_ch        &
                                            , bounds_xc     = bounds_xc        &
                                            , j_array       = j_array          )
        !
        call comm_reduce(j_array)
        !
        if ( comm_rank() == 0 ) then
          !
          i_start = 1
          i_stop  = 0
          do i_gamma = 1, n_irrep
            ! Address logic
            dim_irrep = IrrBasDimSpor(i_gamma)
            n_triag   = dim_irrep * ( dim_irrep + 1 )
            i_stop    = i_stop + n_triag
            call direct_2c_energy_calc_and_add(                                &
                               f_t_arr    = 1.0_r8_kind                        &
                             , f_trace    = 1.0_r8_kind                        &
                             , t_arr      = j_array(i_start:i_stop)            &
                             , d_mat_real = d_matrix_real(i_gamma)%m(:,:)      &
                             , d_mat_imag = d_matrix_imag(i_gamma)%m(:,:)      &
                             , e_2c       = e_coul                             &
                             , h_mat_real = h_matrix_real(i_gamma)%m(:,:)      &
                             , h_mat_imag = h_matrix_imag(i_gamma)%m(:,:)      )
            i_start = i_start + n_triag
          enddo
        endif
        !
      endif
      !
      ! master sets and writes coulomb energy
      if (pid == 1) then
        !
        call set_energy(coul = e_coul)
        !
        call write_energy(output_unit, print_cou = .TRUE.)
        call write_energy(stdout_unit, print_cou = .TRUE.)
        !
      endif
      !
#ifdef WITH_ERI4C
    endif !.not. J_exact
#endif
    !
    ! end of DENSITY FIT COULOMB TERM
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! EXACT 2-ELECTRON INTERACTION TERMS
    !
    ! Formation of the Coulomb matrix from 4-center integrals.
    ! Computing:                __
    !                          \  '
    !             J    =        >   P    ( i j | k l )
    !              ij          /__,  kl
    !                          k l
    !
    ! as well as the energies:
    !                           __
    !                       1  \  '
    !             E    =   ---  >   P    J
    !              J        2  /__,  ij   ij
    !                          i j
    !
    !
    ! Formation of the exact exchange matrix from 4-center integrals
    ! Computing:                __
    !                          \  '
    !             K_il =    -   >   P_jk ( i j | k l )
    !                          /__,
    !                          j k
    !
    ! whereas an additional factor of 1/2 enters in the UKS case
    ! as well as the energies:  __
    !                       1  \  '
    !             E_K  = - ---  >   P_ij K_ij
    !                       2  /__,
    !                          i j
    !
    !
    ! Whereas:
    !          i, j, k, l   :  KS orbitals
    !
    !
    !          P            :  density matrix elements
    !           ij
    !
    !        ( i j | k l )  :  four-center electron repulsion integral
    !
    !
#if WITH_ERI4C == 1
    if ( J_exact .or. K_exact ) then ! includes also J_exact = K_exact = T
      if (options_spin_orbit) then
         call error_handler("coulomb_term: spin orbit in combination with "    &
                          //"exact exchange / coulomb not yet implemented")
      endif
      !
      call exact_2e_interaction (loop, d_matrix, h_matrix, e_coul, e_exex)
      !
    endif ! J_exact .or. K_exact
    !
    !
    ! master sets and writes coulomb energy
    if (pid == 1) then
      !
      if (J_exact) then
        !
        call set_energy(coul = e_coul)
        !
        call write_energy(output_unit, print_cou = .TRUE.)
        call write_energy(stdout_unit, print_cou = .TRUE.)
        !
      endif
      !
      if (K_exact) then
        !
        call set_energy(exex = e_exex)
        !
        call write_energy(output_unit, print_exx = .TRUE.)
        call write_energy(stdout_unit, print_exx = .TRUE.)
        !
      endif
    endif
#endif
    !
    ! end of EXACT 2-ELECTRON INTERACTION TERMS
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
#ifdef WITH_DFTPU
    if( dft_plus_u_mo_in_use )then
      !
      call dft_plus_u_proj( loop, unique_atoms, d_matrix, overlap, h_matrix, e_dft_plus_u )
      !
      ! copy the value into energy_calc_module:
      call set_energy(dft_plus_u=e_dft_plus_u)
      !
    endif
#endif
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! Allocate tmp linear arrays
    deallocate(j_array )
    if (model_density) then
       deallocate(m_array)
    endif
    if (exchange_fit) then
       deallocate(f_array)
       deallocate(x_array)
    endif
    !
  end subroutine build_2e_hamiltonian
  !*****************************************************************************

#if WITH_ERI4C == 1
  !*****************************************************************************
  subroutine exact_2e_interaction( loop                                        &
                                 , d_matrix                                    &
                                 , h_matrix                                    &
                                 , trace_J                                     &
                                 , trace_K                                     )
    !        The calculations are done in respective subroutines whereas
    !        this routine does the coordination of receiving/sending
    !         and allocating the data.
    !
    ! Input-Parameters
    !
    !
    !------------ Modules ----------------------------------------
    use comm,                 only : comm_rank                                 &
                                   , comm_bcast
    use eri_main,             only : eri_blk_main
    use unique_atom_module,   only : unique_atoms
    use energy_calc_module,   only : direct_mat_energy_calc_and_add
#ifdef WITH_ERI4C
    use eri4c_options,        only : J_exact                                   &
                                   , K_exact                                   &
                                   , QQP_tol                                   &
                                   , EXX_factor                                &
                                   , DeltaSCF
#endif
    use symm_positions,       only : symm_mapping_shells
    use group_module,         only : group_num_el
    use symmetry_data_module, only : symmetry_data_n_spin                      &
                                   , symmetry_data_n_irreps
    use direct_scf,           only : get_diff_P                                &
                                   , get_full_K                                &
                                   , get_full_J                                &
                                   , allocate_as
    implicit none
    !------------ Declaration of formal parameters ---------------
    !
    integer(i4_kind),        intent(in)    :: loop
    !
    type(arrmat3),           intent(in)    :: d_matrix(:)
    !
    type(arrmat3),           intent(inout) :: h_matrix(:) ! will be incremented
    !
    real(r8_kind),           intent(out)   :: trace_J
    real(r8_kind),           intent(out)   :: trace_K
    !
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !
    real(r8_kind)                          :: factor_K_mat
    real(r8_kind)                          :: factor_K_eny
    real(r8_kind)                          :: factor_J_mat
    real(r8_kind)                          :: factor_J_eny
    !
    type(arrmat2), allocatable             :: work_j_matrix(:)
    type(arrmat3), allocatable             :: work_k_matrix(:)
    type(arrmat3), allocatable             :: diff_d_matrix(:)
    !
    integer(i4_kind)                       :: i_irr
    integer(i4_kind)                       :: n_irr
    integer(i4_kind)                       :: i_spin
    integer(i4_kind)                       :: n_spin
    integer(i4_kind)                       :: n_atm, n_sym
    integer(i4_kind)                       :: n_shl
    integer(i4_kind), allocatable          :: symm_map(:,:)
    !
    !------------ Executable code ---------------------------------------------+
    !
    ! SET GENERAL SIZES and PREFACTORS
    n_atm  = sum( unique_atoms(:)%n_equal_atoms )
    n_shl  = sum( unique_atoms(:)%n_equal_atoms                                &
                * ( unique_atoms(:)%lmax_ob + 1 ) )
    n_irr  = symmetry_data_n_irreps()
    n_spin = symmetry_data_n_spin()
    !
    ! EXCHANGE MATRIX PREFACTOR IS -1 FOR UKS AND -1/2 FOR RKS
    ! due to the absence of opposite spin exchange interactions
    ! COULOMB MATRIX PREFACTOR IS +1
    ! ENERGY PREFACTORS ARE 1/2
    ! they result from the pairwise interaction of electrons
    ! i.e. from the formula E_ee = 1/2 * tr( P M_ee )
    factor_K_mat =-0.50_r8_kind * n_spin * EXX_factor
    factor_K_eny = 0.50_r8_kind * factor_K_mat
    factor_J_mat = 1.00_r8_kind
    factor_J_eny = 0.50_r8_kind * factor_J_mat
    !
    ! INITIALIZE EXPECTATION VALUES
    trace_J = 0.0_r8_kind
    trace_K = 0.0_r8_kind
    !
    ! GET SYMMETRY MAPPING INDICES FOR INDIVIDUAL SHELLS
    n_sym = group_num_el
    call comm_bcast( n_sym )
    IF ( comm_rank() == 0 ) THEN
      call symm_mapping_shells( unique_atoms, n_atm,                  symm_map )
    ELSE
      ALLOCATE( symm_map(n_sym, n_atm) )
    ENDIF
    call comm_bcast( symm_map )
    !
    ! TODO: move to right place
    IF ( J_exact ) stop 'J_exact not supported currently'
    !
    ! TODO: find a way to avoid excessive usage of intermediates
    call allocate_as( d_matrix, work_k_matrix )
    !
    IF ( DeltaSCF ) THEN
      !
      ! BUILD MATRICES WITH DELTA SCF
      call get_diff_P( d_matrix, diff_d_matrix )
      !
      IF     ( K_exact .and. .not. J_exact ) THEN
        !
        ! BUILD EXACT EXCHANGE MATRIX
        call eri_blk_main( QQP_tol, n_sym, n_atm                               &
                         , n_spin, n_shl, n_irr, loop                          &
                         , diff_d_matrix,                        work_k_matrix &
                         , symm_map, unique_atoms                              )
        !
      ELSEIF ( K_exact .and.       J_exact ) THEN
        !
        ! BUILD EXACT COULOMB AND EXCHANGE MATRICES
!       call eri4c_to_PG( loop, unique_atoms, QQP_tol, diff_d_matrix           &
!                       , jmat=work_j_matrix, kmat=work_k_matrix               )
        !
      ELSEIF ( J_exact .and. .not. K_exact ) THEN
        !
        ! BUILD EXACT COULOMB MATRIX
!       call eri4c_to_PG( loop, unique_atoms, QQP_tol, diff_d_matrix           &
!                       , jmat=work_j_matrix                                   )
        !
      ENDIF
      !
      ! Now convert the differential matrices to full ones
      !
      ! TODO: find a way to avoid excessive usage of intermediates
      IF ( K_exact ) THEN
        call get_full_K( work_k_matrix )
      ENDIF
      IF ( J_exact ) THEN
        call get_full_J( work_j_matrix )
      ENDIF
      !
    ELSE
      !
      ! BUILD MATRICES WITH STANDARD DIRECT SCF
      IF     ( K_exact .and. .not. J_exact ) THEN
        !
        ! BUILD EXACT EXCHANGE MATRIX
        call eri_blk_main( QQP_tol, n_sym, n_atm                               &
                         , n_spin, n_shl, n_irr, loop                          &
                         , d_matrix,                             work_k_matrix &
                         , symm_map, unique_atoms                              )
        !
      ELSEIF ( K_exact .and.       J_exact ) THEN
        !
        ! BUILD EXACT COULOMB AND EXCHANGE MATRICES
!       call eri4c_to_PG( loop, unique_atoms, QQP_tol, d_matrix                &
!                       , jmat=work_j_matrix, kmat=work_k_matrix               )
        !
      ELSEIF ( J_exact .and. .not. K_exact ) THEN
        !
        ! BUILD EXACT COULOMB MATRIX
!       call eri4c_to_PG( loop, unique_atoms, QQP_tol, d_matrix                &
!                       , jmat=work_j_matrix                                   )
        !
      ENDIF
    ENDIF
    !
    !
    IF ( comm_rank() == 0 ) THEN
      !
      ! CALCULATE TRACES AND ADD TO HAMILTONIAN (MASTER ONLY)
      ! For reduction of matrices see subroutine eri4c_to_PG
      IF (J_exact) THEN
        DO i_irr = 1, n_irr
          DO i_spin = 1, n_spin
            !
            ! SUBSECTION OF J
            !
            call direct_mat_energy_calc_and_add(                               &
                        f_trace    = factor_J_eny                              &
                      , f_t_mat    = factor_J_mat                              &
                      , t_mat_real = work_j_matrix(i_irr)%m                    &
                      , d_mat_real = d_matrix(i_irr)%m(:,:,i_spin)             &
                      , e_tr       = trace_J                                   &
                      , h_mat_real = h_matrix(i_irr)%m(:,:,i_spin)             )
          ENDDO
        ENDDO
      ENDIF
      !
      IF (K_exact) THEN
        DO i_irr = 1, n_irr
          DO i_spin = 1, n_spin
            !
            ! SUBSECTION OF K
            !
            call direct_mat_energy_calc_and_add(                               &
                        f_trace    = factor_K_eny                              &
                      , f_t_mat    = factor_K_mat                              &
                      , t_mat_real = work_k_matrix(i_irr)%m(:,:,i_spin)        &
                      , d_mat_real = d_matrix(i_irr)%m(:,:,i_spin)             &
                      , e_tr       = trace_K                                   &
                      , h_mat_real = h_matrix(i_irr)%m(:,:,i_spin)             )
          ENDDO
        ENDDO
      ENDIF
      !
    ENDIF
    !
  end subroutine exact_2e_interaction
#endif
  !*****************************************************************************

  !*****************************************************************************
  subroutine density_fit_coulomb_interaction( bounds_ch                        &
                                            , bounds_xc                        &
                                            , j_array                          &
                                            , m_array                          &
                                            , f_array                          &
                                            , x_array                          )
    !
    !        The calculations are done in respective subroutines whereas
    !        this routine does the coordination of receiving/sending
    !         and allocating the data.
    !
    ! Input-Parameters
    !
    !    bounds_ch(:) - number of chargefit coefficients per proc
    !    bounds_xc(:) - number of exch.-fit coefficients per proc
    !
    !
    !------------ Modules ----------------------------------------
    use comm, only: comm_reduce, comm_rank, comm_size, comm_scatter
    use msgtag_module
    use symmetry_data_module
    use dimensions,         only: IrrBasDim                                    &
                                , IrrBasDimSpor
    use strings, only: itoa
    use prepare_integralfiles_module
    use fit_coeff_module, only: coeff_charge, coeff_xcmda, coeff_xcmda_ts &
        , coeff_xc, coeff_en_xc, coeff_spin
    use readwriteblocked_module, only: readwriteblocked_tapehandle
    use iounitadmin_module
    use energy_calc_module, only: energy_calc_reduce

    implicit none
    !------------ Declaration of formal parameters ---------------
    !
    integer(i4_kind), intent(in)           :: bounds_ch(:)      ! (NPROCS)
    integer(i4_kind), intent(in)           :: bounds_xc(:)      ! (NPROCS)
    !
    ! Lower triangles of coulomb matrices
    real(r8_kind)          , intent(inout) :: j_array(:)
    ! Lower triangles of mda xc matrices
    real(r8_kind), optional, intent(inout) :: m_array(:)
    ! Lower triangles of xc-fit energy density matrices
    real(r8_kind), optional, intent(inout) :: f_array(:)
    ! Lower triangles of xc-fit potential matrices
    real(r8_kind), optional, intent(inout) :: x_array(:)
    !
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    type(readwriteblocked_tapehandle) :: th_ch,th_xc,th_ch_imag,th_xc_imag
    integer(kind=i4_kind) :: n_xc, n_ch, info, i_gamma, i_spin, &
         i_meta_ch, i_meta_xc
    integer(i4_kind) :: pid ! process ID == MPI rank + 1
    logical :: exchange_fit, numeric_exch, model_density

    integer(i4_kind)              :: i_start
    integer(i4_kind)              :: i_stop
    integer(i4_kind)              :: i_start_m
    integer(i4_kind)              :: i_stop_m
    integer(i4_kind)              :: n_triag

    integer(i4_kind)              :: ilower, iupper
    integer(i4_kind)              :: n_irrep
    integer(i4_kind)              :: n_spin
    integer(i4_kind)              :: dim_irrep

    !
    ! These names are intentionally kept similar to the names of
    ! variables declared in fit_coeff_module:
    !
    real(r8_kind), allocatable :: koeff_charge(:)
    real(r8_kind), allocatable :: koeff_xcmda(:, :), koeff_xcmda_ts(:, :)
    real(r8_kind), allocatable :: koeff_xc(:, :), koeff_en_xc(:)
    real(r8_kind) :: coeff_unused(0)

    !------------ Executable code --------------------------------

    exchange_fit  = options_xcmode() == xcmode_exchange_fit
    numeric_exch  = options_xcmode() == xcmode_numeric_exch
    model_density = options_xcmode() == xcmode_model_density .or. &
                    options_xcmode() == xcmode_extended_mda

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       n_spin  = 1
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym%n_irrep
       n_spin  = ssym%n_spin
    endif ! options_spin_orbit

    ! base-1 process ID:
    pid = 1 + comm_rank()

    ! number of charge- and xc fit functions assigned to me:
    n_ch = bounds_ch(pid)
    n_xc = bounds_xc(pid)
    !
    ! NOTE: number of xc fit functions supposed to
    !       be zero if xc fit is not in use!
    !

    ! -----------------------------------------------------------
    ! bounds_<xx>(i) : number of fit functions on the i-th CPU
    !                  with the master being the first CPU
    ! -----------------------------------------------------------

    !
    ! Distribute fit coefficients:
    !
    if (model_density) then
        ! NOTE: charge_coeff has already been sent to each slave!

        allocate( koeff_xcmda(n_ch, n_spin)                                    &
                , koeff_xcmda_ts(n_ch, n_spin), stat=info)
        ASSERT(info==0)

        ! FIXME: see default branch if this fires:
        ASSERT(allocated(coeff_xcmda))
        ASSERT(allocated(coeff_xcmda_ts))
        do i_spin = 1, n_spin
            call comm_scatter(coeff_xcmda(:, i_spin), koeff_xcmda(:, i_spin), bounds_ch)
            call comm_scatter(coeff_xcmda_ts(:, i_spin), koeff_xcmda_ts(:, i_spin), bounds_ch)
        enddo
    else
        !
        ! Default branch:
        !
        allocate(koeff_charge(n_ch), stat=info)
        ASSERT(info==0)

        !
        ! FIXME: coeff_charge is not associated on slaves, it is illegal to
        !        pass unassociated pointers where a plain array is expected.
        !        Maybe allocate(coeff_charge(0)) on slaves instead?
        !

        !
        ! First argument of comm_scatter() is unused on receiver side:
        !
        if ( pid == 1 ) then
          call comm_scatter(coeff_charge, koeff_charge, bounds_ch)
        else
          call comm_scatter(coeff_unused, koeff_charge, bounds_ch)
        endif
    endif

    if (exchange_fit) then
        allocate(koeff_xc(n_xc, n_spin), koeff_en_xc(n_xc), stat=info)
        ASSERT(info==0)

        ! FIXME: see default branch if this fires:
        ASSERT(allocated(coeff_xc))
        do i_spin = 1, n_spin
            call comm_scatter(coeff_xc(:, i_spin), koeff_xc(:, i_spin), bounds_xc)
        enddo
        ASSERT(allocated(coeff_en_xc))
        call comm_scatter(coeff_en_xc, koeff_en_xc, bounds_xc)
    endif

    ! --- ready with distributing data among workers ---

    if ( options_integrals_on_file() ) then
       ! open the integral files ----------------------
       if (options_spin_orbit) then
          if (n_ch.ne.0) &
               call prepare_integral_open('coul', th_ch, th_ch_imag)
          if (n_xc.ne.0 .and. exchange_fit) &
               call prepare_integral_open('exch', th_xc, th_xc_imag)
       else
          if (n_ch.ne.0) &
               call prepare_integral_open('coul', th_ch)
          if (n_xc.ne.0 .and. exchange_fit) &
               call prepare_integral_open('exch', th_xc)
       endif
    else
       i_meta_ch = 1
       i_meta_xc = 1
    endif

    !
    ! BUILD THE COULOMB MATRIX USING DENSITY EXPANSION COEFFICIENTS:
    !

    ! obligate calculation first
    ! decide how to call 'ham_calc' ------------
    if (n_ch.ne.0) then
       if (.not. model_density) then
          !
          ! NOTE: Master holds both a "long" coeff_charge, and
          ! a "short" koeff_charge, but slaves allocate and receive only a
          ! shorter section of that.
          !
          if (.not. options_spin_orbit) then
             !
             ! STANDARD SCF
             !
             i_start = 1
             i_stop  = 0
             do i_gamma=1,n_irrep
                dim_irrep = IrrBasDim(i_gamma)
                n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
                i_stop    = i_stop + n_triag
                !
                call ham_calc_ch( ssym                                         &
                                , i_gamma                                      &
                                , th_ch                                        &
                                , i_meta_ch                                    &
                                , koeff_charge                                 &
                                , j_subsection = j_array(i_start:i_stop)       )
                !
                i_start = i_start + n_triag
                !
             enddo ! loop over irreps
             !
          else
             !
             ! SPIN ORBIT
             !
             i_start = 1
             i_stop  = 0
             do i_gamma=1,n_irrep
                dim_irrep = IrrBasDimSpor(i_gamma)
                n_triag   = dim_irrep * ( dim_irrep + 1 )
                i_stop    = i_stop + n_triag
                !
                call ham_calc_ch( ssym                                         &
                                , i_gamma                                      &
                                , th_ch                                        &
                                , i_meta_ch                                    &
                                , koeff_charge                                 &
                                , th_ch_imag    = th_ch_imag                   &
                                , j_subsection  = j_array(i_start:i_stop)      )
                !
                i_start = i_start + n_triag
                !
             enddo ! loop over irreps
             !
          endif
          !
       else ! that is model_density
          !
          ! FIXME: In this branch all workers appear to hold
          ! the "long" version of coeff_charge. But see below.
          !
          ! My range of fit functions (ilower:iupper), apparently
          ! used only with model_density:
          !
          i_start   = 1
          i_stop    = 0
          i_start_m = 1
          i_stop_m  = 0
          do i_gamma=1,n_irrep
             dim_irrep = IrrBasDim(i_gamma)
             n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
             i_stop    = i_stop   + n_triag
             i_stop_m  = i_stop_m + n_triag * n_spin
             ilower = sum(bounds_ch(1:pid-1)) + 1
             iupper = sum(bounds_ch(1:pid  )) ! == ilower + n_ch - 1
             !
             call ham_calc_ch( ssym                                            &
                             , i_gamma                                         &
                             , th_ch                                           &
                             , i_meta_ch                                       &
                             , coeff_charge(ilower:iupper)                     &
                             , koeff_xcmda                                     &
                             , j_subsection = j_array(i_start:i_stop)          &
                             , m_subsection = m_array(i_start_m:i_stop_m)      )
             !
             i_start   = i_start   + n_triag
             i_start_m = i_start_m + n_triag * n_spin
             !
          enddo ! loop over irreps
          !
       endif
    endif

    if (n_xc.ne.0 .and. exchange_fit) then
       !
       i_start   = 1
       i_stop    = 0
       i_start_m = 1
       i_stop_m  = 0
       do i_gamma=1,n_irrep
          dim_irrep = IrrBasDim(i_gamma)
          n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
          i_stop    = i_stop   + n_triag
          i_stop_m  = i_stop_m + n_triag * ssym%n_spin
          !
          call ham_calc_xc( ssym                                               &
                          , i_gamma                                            &
                          , th_xc                                              &
                          , i_meta_xc                                          &
                          , koeff_xc                                           &
                          , koeff_en_xc                                        &
                          , f_subsection = f_array(i_start:i_stop)             &
                          , x_subsection = x_array(i_start_m:i_stop_m)         )
          !
          i_start   = i_start   + n_triag
          i_start_m = i_start_m + n_triag * n_spin
          !
       enddo
    endif

    !
    ! There are some not irrep specific things to be summed up
    ! over workers, like energy contributions of various types.
    !
    ! We delegate reduction of partial energy contributions
    ! stored in module private variables to the energy_calc_module:
    !
    ! TODO: Check if call below necessary. Currently energies should be
    ! calculated on master only.
    call energy_calc_reduce()

    if ( options_integrals_on_file() ) then
       ! close the integral files ------------------
       if (options_spin_orbit) then
          if (n_ch.ne.0) &
               call prepare_integral_close(th_ch=th_ch,th_ch_imag=th_ch_imag)
          if (n_xc.ne.0 .and. exchange_fit) &
               call prepare_integral_close(th_xc=th_xc,th_xc_imag=th_xc_imag)
       else
          if (n_ch.ne.0) &
               call prepare_integral_close(th_ch=th_ch)
          if (n_xc.ne.0 .and. exchange_fit) &
               call prepare_integral_close(th_xc=th_xc)
       endif
    endif

    if( pid /= 1 ) then
       !
       ! FIXME: not used anywhere in this module:
       !
       if(allocated(coeff_spin)) then
          deallocate( coeff_spin,STAT=info)
          ASSERT(info==0)
       endif
    endif

    if(allocated(koeff_charge))then
        deallocate(koeff_charge, STAT=info)
        ASSERT(info==0)
    endif

    if(allocated(koeff_xc))then
        deallocate(coeff_xc, STAT=info)
        ASSERT(info==0)
    endif

    if(allocated(koeff_en_xc))then
        deallocate(coeff_en_xc, STAT=info)
        ASSERT(info==0)
    endif

    if(allocated(koeff_xcmda))then
        deallocate(coeff_xcmda, STAT=info)
        ASSERT(info==0)
    endif

    if(allocated(koeff_xcmda_ts))then
        deallocate(coeff_xcmda_ts, STAT=info)
        ASSERT(info==0)
    endif
  end subroutine density_fit_coulomb_interaction
  !*****************************************************************************

  !*****************************************************************************
  subroutine ham_calc_ch( ssym                                                 &
                        , i_gamma                                              &
                        , th_ch                                                &
                        , i_meta                                               &
                        , coeff_charge                                         &
                        , coeff_xcmda                                          &
                        , th_ch_imag                                           &
                        , j_subsection                                         &
                        , m_subsection                                         )
    !
    !  Purpose: In this routine the following parts of
    !           the Hamiltonian are constructed:
    !     1. Coulomb-part: j_subsection (lower triangle)
    !        j_array = sum_k { a_k [mu nu | f_k ]}
    !        [mu nu | f_k ] := 3-center Integral / read in from file
    !        ak             := fit-coefficient for charge density
    !
    !     2. MDA xc-part: m_subsection (lower triangle)
    !        m_array = sum_k { b"_k [mu nu | f_k ]}
    !        [mu nu | f_k ] := 3-center Integral / read in from file
    !        b"_k           := Coulomb-type fit-coefficient for exchange
    !                          potential within the MDA approach
    !
    !  The calcualation of the exch.-corr. potential is OPTIONAL here.
    !  It will only be needed for the case that the MDA approach is used
    !
    !  Input parameter (not modified on output):
    !  Name:              Description/Range:
    !  ssym               symmetry information
    !  i_gamma            current IRREP
    !  coeff_charge       partial array of charge fit coeffs
    !  coeff_xcmda           "     "       exch. fit coeffs (Coulomb-type)
    !
    ! Output:
    !  the variables j_subsection and m_subsection (optional)
    !  contain the lower triangular part of the matrices resulting from
    !  the Coulomb term and the MDA XC part
    !
    ! Subroutine called by: density_fit_coulomb_interaction
    !
    ! References : see documentation on old lcgto
    !              (harrharr)
    ! Author: Folke Noertemann
    ! Date  : 7/95 , revised for modularization 10/95
    !
    ! Modification:
    ! Author: TG
    ! Date: 10/96
    ! Description: Removed spin dependency and
    !              very superfluent calculations
    !
    ! Modification:
    ! Author: UB
    ! Date: 7/97
    ! Description: use_model_density option introduced
    !              direct_energy_calc options introduced
    !
    ! Modification:
    ! Author: MM
    ! Date: 10/97
    ! Description: extension to spin orbit
    !
    ! Modification:
    ! Author: TS
    ! Date: 07/11
    ! Description: reorganizations. avoid intermediate quantities
    !
    ! Modification (Please copy before editing)
    ! Author: ...
    ! Date:   ...
    ! Description: ...
    !------------ Modules ---------------------------------------
    use filename_module
    use symmetry_data_module, only : sym
    use readwriteblocked_module
    use integralstore_module, only: integralstore_3c_co
    use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    !
    type(sym),                                   intent(in)    :: ssym
    integer(i4_kind),                            intent(in)    :: i_gamma
    type(readwriteblocked_tapehandle),           intent(inout) :: th_ch
    integer(i4_kind),                            intent(inout) :: i_meta
    real(r8_kind), dimension(:),                 intent(in)    :: coeff_charge
    real(r8_kind), dimension(:,:),     optional, intent(in)    :: coeff_xcmda
    type(readwriteblocked_tapehandle), optional, intent(inout) :: th_ch_imag
    real(r8_kind), dimension(:),                 intent(inout) :: j_subsection
    real(r8_kind), dimension(:),       optional, intent(inout) :: m_subsection
    !
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !
    integer(i4_kind)           :: item_arr_ch
    integer(i4_kind)           :: i_last
    integer(i4_kind)           :: i_arr
    integer(i4_kind)           :: i_arr_s
    integer(i4_kind)           :: dim_ir
    integer(i4_kind)           :: m
    integer(i4_kind)           :: n
    integer(i4_kind)           :: k
    integer(i4_kind)           :: alloc_stat
    integer(i4_kind)           :: spin
    integer(i4_kind)           :: n_spin
    real(r8_kind), pointer     :: coul_int(:)
    real(r8_kind), pointer     :: coul_int_real(:)
    real(r8_kind), pointer     :: coul_int_imag(:)
    logical                    :: model_density
    logical                    :: spread_ham_tot
    logical                    :: spin_polarized
    logical                    :: integrals_on_file

    ! ---------------------------------------------------------------

    model_density = options_xcmode() == xcmode_model_density .or. &
                    options_xcmode() == xcmode_extended_mda
    n_spin = ssym%n_spin
    spin_polarized = n_spin > 1
    spread_ham_tot = spin_polarized                                            &
                   .and. options_xcmode() == xcmode_exchange_fit
    integrals_on_file = options_integrals_on_file()

    item_arr_ch = size(coeff_charge,1)

    if (options_spin_orbit) then
       dim_ir = ssym%dim_proj(i_gamma)
    else
       dim_ir = ssym%dim(i_gamma)
    endif

    if (integrals_on_file) then
       if (options_spin_orbit) then
          allocate( coul_int_real(item_arr_ch),coul_int_imag(item_arr_ch),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       else
          allocate( coul_int(item_arr_ch),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       j_subsection = 0.0_r8_kind
       !
       i_arr = 1
       do m=1,dim_ir
          do n=1,m
             ! read in 3-Z integrals
             if (integrals_on_file) then
                call readwriteblocked_read(coul_int_real,th_ch)
                call readwriteblocked_read(coul_int_imag,th_ch_imag)
             else
                i_last = i_meta + 2*item_arr_ch - 1
                coul_int_real => integralstore_3c_co(i_meta:i_last-1:2)
                coul_int_imag => integralstore_3c_co(i_meta+1:i_last:2)
                i_meta = i_last + 1
             endif
             do k=1,item_arr_ch
                j_subsection(i_arr)     = j_subsection(i_arr)                  &
                                        + coeff_charge(k) * coul_int_real(k)
                ! sign needed to make convention of coul_int_imag compatible
                ! see description in energy_calc_module
                j_subsection(i_arr + 1) = j_subsection(i_arr + 1)              &
                                        - coeff_charge(k) * coul_int_imag(k)
             enddo! k-loop
             i_arr = i_arr + 2
          enddo! n-loop
       enddo!m-loop
       if (integrals_on_file) then
          deallocate(coul_int_real,coul_int_imag,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       else
          nullify(coul_int_real,coul_int_imag)
       endif
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       if ( .not. model_density ) then
          !
          i_arr = 1
          do m=1,dim_ir
             do n=1,m
                ! read in 3-Z integrals
                call get_coul_int() ! piece of code
                do k=1,item_arr_ch
                   j_subsection(i_arr)     = j_subsection(i_arr)               &
                                           + coeff_charge(k) * coul_int(k)
                enddo! k-loop
                i_arr = i_arr + 1
             enddo! n-loop
          enddo!m-loop
       else
          !
          ! i.e. model_density == same as above + more
          ! (absoft doesnt like ifs with optional args in loops?)
          i_arr   = 1
          i_arr_s = 1
          do m=1,dim_ir
             do n=1,m
                ! read in 3-Z integrals
                call get_coul_int() ! piece of code
                do k=1,item_arr_ch
                   j_subsection(i_arr)     = j_subsection(i_arr)               &
                                           + coeff_charge(k) * coul_int(k)
                enddo! k-loop
                i_arr = i_arr + 1
                !
                ! model_density only:
                do spin=1,n_spin
                   do k=1,item_arr_ch
                      ! storage order: alternating alpha & beta contributions
                      m_subsection(i_arr_s) = m_subsection(i_arr_s)            &
                                            + coeff_xcmda(k,spin) * coul_int(k)
                   enddo! k-loop
                   i_arr_s = i_arr_s + 1
                   !
                enddo
             enddo! n-loop
          enddo!m-loop
       endif

      if (integrals_on_file) then
         deallocate(coul_int,STAT=alloc_stat)
         ASSERT(alloc_stat==0)
      else
         nullify(coul_int)
      endif
    endif !options_spin_orbit

  contains
    subroutine get_coul_int()
      use readwriteblocked_module ! routines for blocked I/O
      implicit none
      ! *** end of interface ***

      ! read in 3-Z integrals
      if (integrals_on_file) then
         call readwriteblocked_read(coul_int,th_ch)
      else
         i_last = i_meta + item_arr_ch - 1
         coul_int => integralstore_3c_co(i_meta:i_last)
         i_meta = i_last + 1
      endif
    end subroutine get_coul_int

  end subroutine ham_calc_ch
  !*************************************************************

  !*************************************************************
  subroutine ham_calc_xc( ssym                                                 &
                        , i_gamma                                              &
                        , th_xc                                                &
                        , i_meta                                               &
                        , coeff_xc                                             &
                        , coeff_en_xc                                          &
                        , f_subsection                                         &
                        , x_subsection                                         )
    !
    ! TODO: THIS SUBROUTINE SHOULD BE UNIFIED WITH ham_calc_ch !
    !
    !  Purpose: In this routine the following parts of
    !           the Hamiltonian are constructed:
    !
    !    1.  XC potential ( rho^1/3 - spin dependent): x_subsection
    !        x_array = sum_l{ b_l <mu | g_k | nu> }
    !        <mu | g_l | nu>  := 3-center integral /read in from file
    !        b_l              := fit-coefficients for exchange-potential
    !
    !    2.  XC energy density (rho^4/3 - is a total): f_subsection
    !        f_array = sum_l{ b`_l <mu | g_l | nu> }
    !        b`_l           := fit-coefficients for exchange energy density
    !
    !  The calcualation of the exch.-corr. potential and
    !  exch.-corr. energy density is OPTIONAL here.
    !  It will only be needed for the case that these two
    !  variables will not be calculated numerically.
    !
    !  Input parameter (not modified on output):
    !  Name:              Description/Range:
    !  ssym               symmetry information
    !  i_gamma            current IRREP
    !  coeff_xc              "     "       exch. fit coeffs
    !  coeff_en_xc        dito for exch. energy density
    !
    ! Output:
    !  the variable x_subsection contains the lower triangular part of
    !  the matrix resulting from fitted XC term
    !
    ! Subroutine called by: density_fit_coulomb_interaction
    !
    ! References : see documentation on old lcgto
    !              (harrharr)
    ! Author: Folke Noertemann
    ! Date  : 7/95 , revised for modularization 10/95
    !
    !------------ Modules ---------------------------------------
    use filename_module
    use symmetry_data_module, only : sym
    use readwriteblocked_module
    use integralstore_module, only: integralstore_3c_xc
    use energy_calc_module, only: set_energy
    implicit none
    !------------ Declaration of formal parameters ---------------
    !
    type(sym),                                   intent(in)    :: ssym
    integer(i4_kind),                            intent(in)    :: i_gamma
    type(readwriteblocked_tapehandle),           intent(inout) :: th_xc
    integer(i4_kind),                            intent(inout) :: i_meta
    real(r8_kind), dimension(:,:),               intent(in)    :: coeff_xc
    real(r8_kind), dimension(:),                 intent(in)    :: coeff_en_xc
    real(r8_kind), dimension(:),                 intent(inout) :: f_subsection
    real(r8_kind), dimension(:),                 intent(inout) :: x_subsection
    !
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind)           :: alloc_stat
    integer(i4_kind)           :: item_arr_xc
    integer(i4_kind)           :: i_last
    integer(i4_kind)           :: dim_ir
    integer(i4_kind)           :: m
    integer(i4_kind)           :: n
    integer(i4_kind)           :: k
    integer(i4_kind)           :: spin
    integer(i4_kind)           :: n_spin
    integer(i4_kind)           :: i_arr
    integer(i4_kind)           :: i_arr_s
    real(r8_kind), pointer     :: exch_int(:)
    logical                    :: integrals_on_file
    ! ---------------------------------------------------------------
    integrals_on_file = options_integrals_on_file()

    item_arr_xc = size(coeff_xc,1)

    dim_ir = ssym%dim(i_gamma)
    n_spin = ssym%n_spin

    ! Note that coeff_en_xc actually holds (-1)*d_k and that coeff_xc
    ! actually holds (-1)*c_k here (see build_xcfit). Hence
    ! sum_xc_en = (-1)*Sum(k) d_k < xi_i xi_j | g_k > and
    ! sum_xc    = (-1)*Sum(k) c_k < xi_i xi_j | g_k > !

    if (integrals_on_file) then
       allocate(exch_int(item_arr_xc),STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    endif
    !
    do m = 1, dim_ir
       do n = 1, m
          if (integrals_on_file) then
             call readwriteblocked_read(exch_int,th_xc)
          else
             i_last = i_meta + item_arr_xc - 1
             exch_int => integralstore_3c_xc(i_meta:i_last)
             i_meta = i_last + 1
          endif
          !
          i_arr   = 1
          i_arr_s = 1
          !
          do k=1,item_arr_xc
             ! sign convention: see density_fit_coulomb_interaction
             f_subsection(i_arr)     = f_subsection(i_arr)                     &
                                     + coeff_en_xc(k) * exch_int(k)
          enddo ! k-loop
          i_arr = i_arr + 1
          !
          do spin = 1, n_spin
             do k = 1, item_arr_xc
                ! sign convention: see density_fit_coulomb_interaction
                ! storage order: alternating alpha & beta contributions
                x_subsection(i_arr_s) = x_subsection(i_arr_s)                  &
                                      + coeff_xc(k,spin) * exch_int(k)
             enddo ! k-loop
             i_arr_s = i_arr_s + 1
             !
          enddo!spin - loop
       enddo!n - loop
    enddo!m - loop
    if (integrals_on_file) then
       deallocate( exch_int,STAT=alloc_stat )
       ASSERT(alloc_stat==0)
    else
       nullify(exch_int)
    endif

  end subroutine ham_calc_xc
  !*************************************************************

  !*************************************************************
  subroutine electric_field_term(densmat, ham_tot, e_efield)
    !
    ! Purpose: adds contribution of efield to ham_tot and
    ! calculates contribution of efield to energy. TB 12/97
    !
    ! FIXME: maybe remove "integrals_on_file" case for electric field
    !        and use generic hamiltonian_term(...) instead.
    !        This may not work for "reverse_for_spins" case, though.
    !
    !        Maybe first "slurp" file into memory then call "in memory"
    !        case?
    !
    use symmetry_data_module
    use filename_module, only: tmpfile
    use readwriteblocked_module
    use integralstore_module, only: integralstore_2cob_efield
    use efield_module, only: efield_reverse_for_spins, efield_field
    use dipole_module, only: dipole_nuclear_calculate
    implicit none
    type(arrmat3), intent(in)    :: densmat(:) ! (n_irreps)
    type(arrmat3), intent(inout) :: ham_tot(:) ! (n_irreps)
    real(r8_kind), intent(out)   :: e_efield
    ! *** end of interface ***

    !------------ Declaration of local variables -----------------
    type(readwriteblocked_tapehandle)    :: th_efield
    integer(kind=i4_kind)                :: i_gamma,m,n_dim,i_meta,i_last
    real(kind=r8_kind),pointer           :: buffer_den(:),buffer_efield(:)
    logical :: spin_polarized, integrals_on_file, reverse_for_spins
    integer(i4_kind)                     :: sig ! sign

    spin_polarized = .not. options_spin_restricted()
    integrals_on_file = options_integrals_on_file()
    reverse_for_spins = efield_reverse_for_spins()

    if( reverse_for_spins ) then
      ! relative sign of electric fields acting on two spin directions:
      sig = -1
    else
      sig = +1
    endif

    if (integrals_on_file) then
       call readwriteblocked_startread(trim(tmpfile('efield.dat')), th_efield)
    else
       i_meta = 1
    endif
    e_efield = 0.0_r8_kind
    do i_gamma = 1, symmetry_data_n_irreps()
       n_dim =  symmetry_data_dimension(i_gamma)
       allocate(buffer_den(n_dim))
       if (integrals_on_file) allocate(buffer_efield(n_dim))
       do m = 1,n_dim
          ! only H_ij and P_ij with i <= j is used here
          if (integrals_on_file) then
             call readwriteblocked_read(buffer_efield(1:m),th_efield)
          else
             i_last = i_meta + m - 1
             buffer_efield => integralstore_2cob_efield(i_meta:i_last)
             i_meta = i_last + 1
          endif
          buffer_den(1:m) = densmat(i_gamma)%m(1:m,m,1)
          ! double the off-diagonal elements of the density matrix:
          buffer_den(1:m-1) = buffer_den(1:m-1) * 2

          ! increment the trace with the density matrix:
          e_efield = e_efield + sum(buffer_efield(1:m)*buffer_den(1:m))

          ! update the upper-right triangle of the hamiltonian with diagonal...
          ham_tot(i_gamma)%m(1:m,m,1) = ham_tot(i_gamma)%m(1:m,m,1) + buffer_efield(1:m)
          ! ... and lower-left triangle without diagonal:
          ham_tot(i_gamma)%m(m,1:m-1,1) = ham_tot(i_gamma)%m(m,1:m-1,1) + buffer_efield(1:m-1)

          if (spin_polarized) then
             buffer_den(1:m) = densmat(i_gamma)%m(1:m,m,2)
             ! double the off-diagonal elements of the density matrix:
             buffer_den(1:m-1) = buffer_den(1:m-1) * 2

             ! increment the trace with the density matrix:
             e_efield = e_efield +  sum(buffer_efield(1:m)*buffer_den(1:m)) * sig

             ! update the upper-right triangle of the hamiltonian...
             ham_tot(i_gamma)%m(1:m,m,2) = ham_tot(i_gamma)%m(1:m,m,2) + buffer_efield(1:m) * sig ! NOTE sig factor!
             ! ... and lower-left triangle:
             ham_tot(i_gamma)%m(m,1:m-1,2) = ham_tot(i_gamma)%m(m,1:m-1,2) + buffer_efield(1:m-1) * sig ! NOTE sig factor!
          endif
       enddo
       deallocate( buffer_den )
       if (integrals_on_file) deallocate(buffer_efield)
    enddo
    DPRINT  'electric_field_term: electron interacton energy with electric field =',e_efield

    ! add interaction energy of nuclei with external electric field.
    e_efield = e_efield - dot_product( efield_field(), dipole_nuclear_calculate() )
    DPRINT  'electric_field_term: total interacton energy with electric field =',e_efield

    if (integrals_on_file) call readwriteblocked_stopread(th_efield)
  end subroutine electric_field_term
  !*************************************************************

  !*****************************************************************************
  subroutine hamiltonian_term( integrals, factor                               &
                             , d_matrix                                        &
                             , d_matrix_real, d_matrix_imag                    &
                             , h_matrix                                        &
                             , h_matrix_real, h_matrix_imag                    &
                             , trace                                           )
    !  Purpose: Process integrals (kin, nuc, ...) stored in a linear array
    !           Substitutes all the different paths that were previously
    !           present in the read_ham_kin_nuc subroutine
    !
    !
    !  author: Thomas Soini
    !  date  : 6/11
    !
    !------------ Modules used -------------------------------------------------
    !
    use datatype,                       only: arrmat2                          &
                                            , arrmat3
    use energy_calc_module,             only: direct_2c_energy_calc_and_add
    !
    !------------ Declaration of formal parameters -----------------------------
    !
    implicit none
    real(r8_kind), intent(in)              :: integrals(:)
    real(r8_kind), intent(in)              :: factor
    !
    type(arrmat3), optional, intent(in)    :: d_matrix(:)
    type(arrmat2), optional, intent(in)    :: d_matrix_real(:)
    type(arrmat2), optional, intent(in)    :: d_matrix_imag(:)
    !
    type(arrmat3), optional, intent(inout) :: h_matrix(:)
    type(arrmat2), optional, intent(inout) :: h_matrix_real(:)
    type(arrmat2), optional, intent(inout) :: h_matrix_imag(:)
    !
    real(r8_kind), intent(out)             :: trace
    ! *** end of interface ***

    !------------ Declaration of local variables -------------------------------
    !
    integer(i4_kind)                       :: dim_irrep
    integer(i4_kind)                       :: n_irrep
    integer(i4_kind)                       :: i_irrep
    !
    integer(i4_kind)                       :: n_spin
    integer(i4_kind)                       :: i_spin
    !
    integer(i4_kind)                       :: i_start
    integer(i4_kind)                       :: i_stop
    !
    integer(i4_kind)                       :: n_triag
    !
    !------------ Executable code ----------------------------------------------
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    n_spin  = 1
    n_irrep = 1
    if     (      present(d_matrix)                                  ) then
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       ASSERT(present(h_matrix))
       ASSERT(.not.options_spin_orbit)
       !
       n_spin  = size(h_matrix(1)%m, 3)
       n_irrep = size(h_matrix)
       !
    elseif (      present(d_matrix_real)                             ) then
       !
       ! SPIN ORBIT
       !
       ASSERT(present(d_matrix_imag))
       ASSERT(present(h_matrix_real))
       ASSERT(present(h_matrix_imag))
       ASSERT(options_spin_orbit)
       !
       n_irrep = size(h_matrix_real)
       !
    else
       !
       call error_handler("hamiltonian_term: input variables not matching "//  &
                          "neither spin orbit nor standard scf")
       !
    endif ! options_spin_orbit
    !
    !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    trace    = 0.0_r8_kind
    !
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       i_start = 1
       i_stop  = 0
       do i_irrep = 1, n_irrep
          dim_irrep = size(h_matrix_real(i_irrep)%m, 1)
          ! factor 2 due to complex variable:
          n_triag   = dim_irrep * ( dim_irrep + 1 )
          i_stop    = i_stop + n_triag
          !
          call direct_2c_energy_calc_and_add(                                  &
                  f_t_arr = factor                                             &
                , f_trace = factor                                             &
                , t_arr      = integrals(i_start:i_stop)                       &
                , d_mat_real = d_matrix_real(i_irrep)%m(:,:)                   &
                , d_mat_imag = d_matrix_imag(i_irrep)%m(:,:)                   &
                , e_2c       = trace                                           &
                , h_mat_real = h_matrix_real(i_irrep)%m(:,:)                   &
                , h_mat_imag = h_matrix_imag(i_irrep)%m(:,:)                   )
          !
          i_start = i_start + n_triag
          !
       enddo
       !
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       i_start = 1
       i_stop  = 0
       do i_irrep = 1, n_irrep
          dim_irrep = size(h_matrix(i_irrep)%m, 1)
          n_triag   = dim_irrep * ( dim_irrep + 1 ) / 2
          i_stop    = i_stop + n_triag
          !
          do i_spin = 1, n_spin
             !
             call direct_2c_energy_calc_and_add(                               &
                     f_t_arr = factor                                          &
                   , f_trace = factor                                          &
                   , t_arr      = integrals(i_start:i_stop)                    &
                   , d_mat_real = d_matrix(i_irrep)%m(:,:,i_spin)              &
                   , e_2c       = trace                                        &
                   , h_mat_real = h_matrix(i_irrep)%m(:,:,i_spin)              )
             !
          enddo
          !
          i_start = i_start + n_triag
          !
       enddo
       !
    endif ! options_spin_orbit
    !
  end subroutine hamiltonian_term
  !*****************************************************************************


!--------------- End of module ----------------------------------
end module ham_calc_module

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
module  properties_module
  !---------------------------------------------------------------
  !
  !  Contains various routines for calculting properties up to now the
  !  following properties are available
  !
  !  - population analysis
  !  - splitting of populations
  !  - calculation of dipole_moments, transition_moments, xes-spectra,
  !    ...
  !  - plotting of orbitals, densities, spin-densities, ...
  !  - fragment orbital analysis
  !
  !  Before  entering the  properties  part a  minimal integral  part,
  !  which  simply calculates  the overlap  matrix, is  performed. The
  !  saved_eigenvec.dat  file must  be  available and  must have  been
  !  recorded with the save_eigenvec_all option.
  !
  !
  !  Module called by: main_master
  !
  !  References:
  !
  !
  !  Author: MS
  !  Date: 12/97
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
# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================



  !------------ public functions and subroutines ------------------
  public :: properties_read, properties_write, properties_main

!================================================================
! End of public interface of module
!================================================================

  !------------ Declaration of constants and variables ------------
  ! defaults for input switches
  logical :: df_mulliken          = .false., &
             df_frag_orb_analysis = .false., &
             df_plot_orbitals     = .false., &
             df_dipole            = .false.

  logical :: mulliken,          & ! do population analysis or not
             frag_orb_analysis, & ! do fragment orbital analysis
             plot_orbitals,     & ! plot orbitals
             dipole               ! do dipole moment integrals

  namelist /properties/ mulliken, frag_orb_analysis, plot_orbitals, dipole
!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  subroutine properties_read()
    !
    ! Purpose: read in input.
    !
    use input_module
    implicit none
    !** End of interface *****************************************

    integer :: unit, status

    mulliken = df_mulliken
    frag_orb_analysis = df_frag_orb_analysis
    plot_orbitals = df_plot_orbitals
    dipole = df_dipole

    if ( input_line_is_namelist("properties") ) then
       call input_read_to_intermediate
       unit = input_intermediate_unit()
       read(unit, nml=properties, iostat=status)
       if (status .gt. 0) call input_error( &
            "properties_read: namelist properties")
    endif
  end subroutine properties_read


  subroutine properties_write(iounit)
    !
    ! Write namelist population to iounit in input format.
    !
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer, intent(in) :: iounit
    !** End of interface *****************************************

    call start("PROPERTIES","PROPERTIES_WRITE", &
         iounit,operations_echo_input_level)
    call flag('MULLIKEN         ',mulliken         , df_mulliken         )
    call flag('PLOT_ORBITALS    ',plot_orbitals    , df_plot_orbitals    )
    call flag('FRAG_ORB_ANALYSIS',frag_orb_analysis, df_frag_orb_analysis)
    call flag('DIPOLE           ',dipole           , df_dipole           )
    call stop()

  end subroutine properties_write


  subroutine properties_main()
    !
    ! Purpose: main subroutine for  doing properties.  The first thing
    ! this  procedure does is  to tell  slaves to  call it.  From that
    ! point on it runs in a parallel context.
    !
    use filename_module, only: recover_dir, outfile
    use population_module, only: population_mulliken
    use eigen_data_module, only: eigen_data_alloc, eigen_data_free
    use overlap_module, only: read_overlap
    use readwriteblocked_module, only: readwriteblocked_tapehandle, &
         readwriteblocked_startread, readwriteblocked_startwrite, &
         readwriteblocked_stopread, readwriteblocked_stopwrite
    use mixing_module, only: mixing_state_recover
    use convergence_module, only: convergence_setup, &
         convergence_state_recover, convergence_shutdown
    use fit_coeff_module, only: fit_coeff_recover
    use occupied_levels_module, only: sndrcv_eigvec_occ
    use density_data_module, only: gendensmat_occ, density_data_free
    use iounitadmin_module, only: get_iounit, return_iounit
    use occupation_module, only: eigenstates_recover, eigenstates_store
    use eigen_data_module, only: eigen_data_bcast
    use frag_orb_analysis_module, only: frag_orb_analysis_main
    use orbital_plot_module, only: orbital_plot_main
    use options_module, only: options_save_as_fragment
    use operations_module, only: operations_scf
    use comm_module, only: comm_all_other_hosts, comm_init_send, comm_send
    use msgtag_module, only: msgtag_properties_main
    use comm, only: comm_rank, comm_bcast
    implicit none
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    type(readwriteblocked_tapehandle)   :: th
    integer(i4_kind) :: loop
    external error_handler
    !------------ Executable code --------------------------------

    ! first tell slaves to enter this sub
    if (comm_rank() == 0) then
       ! tell slaves to start plotting
       call comm_init_send(comm_all_other_hosts, msgtag_properties_main)
       call comm_send()
    end if

    if (.not. operations_scf) then
       ABORT('needs more work')
       ! Otherwise  the  eigenvectors  are  still available  from  the
       ! scf-part. If not, allocate:
       call eigen_data_alloc()
       ! Now read eigenvectors, occupations, ... from file
       call say ('reading saved_eigenvec.dat')
       call readwriteblocked_startread(trim(recover_dir)//&
            trim('/saved_eigenvec.dat'), th, variable_length=.true.)
       call fit_coeff_recover(th)
       call convergence_setup()
       call convergence_state_recover(th)
       call convergence_shutdown()
       call mixing_state_recover(loop, th)
       call eigenstates_recover(th=th)
       call readwriteblocked_stopread(th)
    end if

    ! Tell slaves to allocate eigenstaff and broadcast it:
    call eigen_data_bcast()

    call sndrcv_eigvec_occ()

    ! Generate density  matrix. Will allocate density matrix  if it is
    ! not     already     allocated.      See    the     corresponding
    ! density_data_free() towards the end.
    call gendensmat_occ()

    call say ('read_overlap')
    call read_overlap()

    ! Now  everything is  read in  and  initialized and  we can  start
    ! actual  properties.  First  population analysis.   FIXME: serial
    ! O(N^3) code here:
    if (comm_rank() == 0) then
       if (mulliken) then
          call say ('population_mulliken')
          call population_mulliken() ! does no comm
       endif
    endif

    call comm_bcast (plot_orbitals) ! FIXME: paranoya
    if (plot_orbitals) then
       call say ('orbital_plot_main')
       call orbital_plot_main()
    endif

    if (frag_orb_analysis) then
       ABORT("needs more work")

       call say ('frag_orb_analysis_main')
       call frag_orb_analysis_main()
    endif

    if (dipole) then
       ABORT("needs more work")

       call say ('main_dipole')
       call main_dipole()
    endif

    if (options_save_as_fragment().and..not.operations_scf) then
       call say ('save_as_fragment')

       call readwriteblocked_startwrite(trim(outfile("saved_fragment.dat")), &
            th,variable_length=.true.)
       call eigenstates_store(th=th)
       call readwriteblocked_stopwrite(th)
    end if

    ! At least not when operations_scf == true to avoid crashing after
    ! properties_main().   See  finalize_geometry()  for  another  use
    ! later in the program flow:
    if (.not. operations_scf) then
       call eigen_data_free()
    endif

    ! Because we  called gendensmat_occ() unconditionally,  see above.
    ! FIXME: but what if it was already allocated before?
    call density_data_free()

  contains

    subroutine say (phrase)
      use output_module, only: output_properties
      use iounitadmin_module, only: write_to_output_units
      use comm, only: comm_rank
      implicit none
      character(len=*), intent(in) :: phrase
      ! *** end of interface ***

      if (output_properties .and. comm_rank() == 0) then
         call write_to_output_units ("PROPERTIES_MAIN: "//phrase)
      endif
    end subroutine say

  end subroutine properties_main

  !--------------- End of module ----------------------------------
end module properties_module

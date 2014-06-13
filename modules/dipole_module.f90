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
Module dipole_module
  !---------------------------------------------------------------
  !
  !  Purpose: database to hold dipole matrix elements
  !           as well as calculated dipole moments.
  !           and tools to calculate diplols in orbital basis
  !           and to calculate overall dipole moments
  !           and output routines.
  !
  !  The data in this module only exist on the master !!!
  !
  !  Module called by: ...
  !
  !  References: ...
  !
  !  Author: TB
  !  Date: 7/97
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !
  !   Some important troutines:
  !         dipoleg_allocate  dipoleg_free
  !         dipoleg_calculate L973
  !         dipole_transitionmoment_f L1474
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Author: HH
  ! Date:   10/97
  ! Description: Add subroutine "dipol_trans_response".
  !              This subroutine writes all transition dip.mom.
  !              in an temporary, unformatted file in tmp_dir.
  !----------------------------------------------------------------
  !
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   10/99
  ! Description: Adaption to spin-orbit
  !
  !----------------------------------------------------------------
  !
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
# include "def.h"
  Use type_module ! type specification parameters
  Use symmetry_data_module
  Use occupation_module
  Use eigen_data_module
  Use options_module, Only: options_spin_orbit
  Use iounitadmin_module
  Implicit None
  Save                     ! save all variables defined in this module
  Private                       ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !
  !------------ Declaration of constants --------------------------
  Integer, Parameter, Public :: &
       DIPOLE_UNUSED      = 0, &
       DIPOLE_DIAGONAL    = 1, &
       DIPOLE_OFFDIAGONAL = 2
  ! parameters used within dipole_type
  !
  !------------ Declaration of types ------------------------------

  type, public :: dipole_integral_type
     integer :: use = DIPOLE_UNUSED
     ! can be: DIPOLE_UNUSED, DIPOLE_DIAGONAL or DIPOLE_OFFDIAGONAL
     !  for DIPOLE_UNUSED: neither diagonal(:), offdiagonal(:,:) used
     !  for DIPOLE_DIAGONAL: only diagonal(:) used
     !  for DIPOLE_OFFDIAGONAL: only offdiagonal(:,:) used
     real (kind=r8_kind), allocatable :: diagonal(:), offdiagonal (:,:)
     ! diagonal(dim_ir*(dim_ir+1)/2), offdiagonal(dim_ir2,dim_ir1)
     ! integrals are stored in symmetrical matrix fashion for
     ! matrix belonging to diagonal Irrep & Partner - block
     ! triangle storage format used:
     !    ij = 1
     !    do j = 1, n
     !       do i = 1, j
     !          at(ij) = a(i,j)
     !          ij = ij + 1
     !       enddo
     !    enddo
     real (kind=r8_kind), allocatable :: diagonal_imag(:), offdiagonal_imag(:,:)
  end type dipole_integral_type

  type, public :: dipole_type
     real (kind=r8_kind), allocatable :: occupied (:, :), real (:, :)
     ! xxx(i_xyz,i_orbital)
     ! occupied:  dipole moment of orbital with assumed occupation 1
     ! real:      dipole moment of orbital taking into account
     !             its real occupation
     real (kind=r8_kind) :: total (3)
  end type dipole_type

  type (dipole_integral_type), allocatable, target, public :: dipole_integrals(:,:,:)
  ! (3, n_ip, n_ip)
  ! 3    x-, y-, and z- component of dipole integrals
  ! n_ip is a combined index for Irrep and Partner
  !      as defined in symmetry_data_module
  ! only the upper triangle of the matrix is used
  !   i.e. : the secound index is only used up to the value of the third
  ! 3    x-, y-, and z- component of dipole integrals
  ! n_ip is a combined index for Irrep and Partner
  !      as defined in symmetry_data_module
  ! only the upper triangle of the matrix is used
  !   i.e. : the secound index is only used up to the value of the third

  Logical, Public :: dipole_offdiagonals_calculated
  Type (dipole_type), Allocatable, Target, Public :: dipole_moments(:, :) ! (n_spin, n_ip)
  ! n_ip is a combined index for Irrep and Partner
  !      as defined in symmetry_data_module
  ! n_ip is a combined index for Irrep and Partner
  !      as defined in symmetry_data_module
  ! (i_xyz)
  ! (i_xyz,i_spin)

  ! In  dipole_nuclear, both real  atoms and  point charges  are taken
  ! into account:
  real (r8_kind), public, protected :: dipole_total (3) ! x, y, and z
  real (r8_kind), public, protected :: dipole_nuclear (3) ! x, y, and z
  real (r8_kind), public, protected :: dipole_total_spin (3, 2) ! x, y, and z

  real (r8_kind), private :: dipole_surface_charge (3) ! x, y, and z
  real (r8_kind), private :: surface_charge

  real (r8_kind), allocatable, public, protected :: &
       dipole_total_irrep (:, :) ! (i_xyz ,i_ir)

  !------------ public functions and subroutines ------------------
  Public dipole_allocate, dipole_free, dipole_calculate, &
       & dipole_print, dipole_transitionmoment_f, &
       & dipole_transitionmoment_uf, dipole_read, dipole_write, &
       & dipole_trans_response, dipole_make_xes_spectra, &
       & dipole_xes_spectra, dipole_write_simol, dipole_write_optimizer

  public :: dipole_nuclear_calculate !() -> real(3)


  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of input variables -------------------

  ! note that switch in output namelist controls if any dipole
  ! transitionmoments are printed
  Logical, Private :: show_valence_excitations = .True.
  ! if this switch is set, dipole transition moments between occupied
  ! orbital and empty orbital are printed, if the corresponding
  ! excitation energy is below max_valence_excitation_energy
  Logical, Private :: show_core_excitations = .False.
  ! if this switch is set, core excitaition from core levels in the
  ! energy range between core_excitation_from_e_min and
  ! core_excitation_from_e_max to levels in the energy range between
  ! core_excitation_to_e_min and core_excitation_to_e_max are printed
  Real (Kind=r8_kind), Private :: max_valence_excitation_energy = &
       & 10.0_r8_kind, core_excitation_from_e_min, &
       & core_excitation_from_e_max, core_excitation_to_e_min, &
       & core_excitation_to_e_max

  ! explications are given above
  Logical :: make_xes_spectra ! calculate a xes spectrum

  ! the core_hole is described by the nect 4 variables number of the
  ! hole in the irrep
  Integer (Kind=i4_kind) :: core_hole_irrep, core_hole_spin, &
       & i_core_hole, core_hole_partner
  Real (Kind=r8_kind) :: x, y, z ! define the direction of the E_field
                                 ! of the emitted light


  ! Default valuse of the input variables
  Logical :: df_show_valence_excitations = .True., &
       & df_show_core_excitations = .False., df_make_xes_spectra = &
       & .False.

  Real (Kind=r8_kind) :: df_max_val_excitation_energy = &
       & 10.0_r8_kind, df_core_excitation_from_e_min = 0.0_r8_kind, &
       & df_core_excitation_from_e_max = 0.0_r8_kind, &
       & df_core_excitation_to_e_min = 0.0_r8_kind, &
       & df_core_excitation_to_e_max = 0.0_r8_kind, df_x = 0.0_r8_kind, &
       & df_y = 0.0_r8_kind, df_z = 0.0_r8_kind
  Integer (Kind=i4_kind) :: df_core_hole_irrep = 0, &
       & df_core_hole_spin = 0, df_core_hole_partner = 0, df_i_core_hole &
       & = 0

  Namelist / dipole_transitionmoments / show_valence_excitations, &
       & show_core_excitations, max_valence_excitation_energy, &
       & core_excitation_from_e_min, core_excitation_from_e_max, &
       & core_excitation_to_e_min, core_excitation_to_e_max, &
       & make_xes_spectra, core_hole_irrep, core_hole_spin, i_core_hole, &
       & core_hole_partner, x, y, z

  ! Stores a  unit vector along  the direction specified in  the input
  ! file:
  real (r8_kind), private :: x_vec (3) ! store x,y,z; see below

  !------------ Subroutines ---------------------------------------
Contains
!
  !*************************************************************
      Subroutine dipole_read ()
    ! purpose: read namelist dipole_transitionmoments that
    ! controls formatted output of dipole transition moments
    ! routine called by: read_input
    !** End of interface *****************************************
         Use input_module
    !------------ Declaration of local variables -----------------
         Integer (Kind=i4_kind) :: unit, status
         Real (Kind=r8_kind) :: norm
         External error_handler
    !------------ Executable code --------------------------------
!
         show_valence_excitations = df_show_valence_excitations
         show_core_excitations = df_show_core_excitations
         max_valence_excitation_energy = df_max_val_excitation_energy
         core_excitation_from_e_min = df_core_excitation_from_e_min
         core_excitation_from_e_max = df_core_excitation_from_e_max
         core_excitation_to_e_min = df_core_excitation_to_e_min
         core_excitation_to_e_max = df_core_excitation_to_e_max
         make_xes_spectra = df_make_xes_spectra
         core_hole_irrep = df_core_hole_irrep
         core_hole_spin = df_core_hole_spin
         core_hole_partner = df_core_hole_partner
         i_core_hole = df_i_core_hole
         x = df_x
         y = df_y
         z = df_z
!
!
         If (input_line_is_namelist("DIPOLE_TRANSITIONMOMENTS")) Then
            Call input_read_to_intermediate
            unit = input_intermediate_unit ()
            Read (Unit, Nml=dipole_transitionmoments, IoStat=Status)
            If (status .Gt. 0) Call input_error ("DIPOLE_READ_INPUT: na&
           &melist dipole_transitionmoments")
         End If
         If (show_core_excitations .And. (core_excitation_from_e_min &
        & .Eq. 0.0_r8_kind .Or. core_excitation_from_e_max .Eq. &
        & 0.0_r8_kind .Or. core_excitation_to_e_min .Eq. 0.0_r8_kind &
        & .Or. core_excitation_to_e_max .Eq. 0.0_r8_kind)) Call &
        & input_error&
        & ("DIPOLE_READ_INPUT: "//"You set show_core_excitations but give no energy ranges")
!
         If (make_xes_spectra) Then
            norm = x * x + y * y + z * z
            If (norm == 0.0_r8_kind) Call input_error &
            & ('DIPOLE_READ_INPUT: x, y, z must be set')
            x_vec (1) = x / Sqrt (norm)
            x_vec (2) = y / Sqrt (norm)
            x_vec (3) = z / Sqrt (norm)
            If (core_hole_irrep == 0 .Or. core_hole_spin == 0 .Or. core_hole_partner == 0 .Or. i_core_hole == 0) Call input_error &
           & ('DIPOLE_READ_INPUT: core_hole is not defined')
         End If
      End Subroutine dipole_read
  !*************************************************************
      Subroutine dipole_write (iounit)
        !
        ! Purpose: writes namelist dipole_transitionmoments
        !
        use echo_input_module, only: start, real, flag, intg, strng, stop, &
             echo_level_full
        Use operations_module, Only: operations_echo_input_level
        implicit none
        Integer, Intent (In) :: iounit
        !** End of interface ***************************************

         Call start ("DIPOLE_TRANSITIONMOMENTS", "DIPOLE_WRITE", &
        & iounit, operations_echo_input_level)
         Call flag ("SHOW_VALENCE_EXCITATIONS     ", &
        & show_valence_excitations, df_show_valence_excitations)
         Call flag ("SHOW_CORE_EXCITATIONS        ", &
        & show_core_excitations, df_show_core_excitations)
         Call real ("MAX_VALENCE_EXCITATION_ENERGY", &
        & max_valence_excitation_energy, df_max_val_excitation_energy)
         Call real ("CORE_EXCITATION_FROM_E_MIN   ", &
        & core_excitation_from_e_min, df_core_excitation_from_e_min)
         Call real ("CORE_EXCITATION_FROM_E_MAX   ", &
        & core_excitation_from_e_max, df_core_excitation_from_e_max)
         Call real ("CORE_EXCITATION_TO_E_MIN     ", &
        & core_excitation_to_e_min, df_core_excitation_to_e_min)
         Call real ("CORE_EXCITATION_TO_E_MAX     ", &
        & core_excitation_to_e_max, df_core_excitation_to_e_max)
         Call flag ("MAKE_XES_SPECTRA             ", make_xes_spectra, &
        & df_make_xes_spectra)
         Call intg ("CORE_HOLE_IRREP              ", core_hole_irrep, &
        & df_core_hole_irrep)
         Call intg ("CORE_HOLE_SPIN               ", core_hole_spin, &
        & df_core_hole_spin)
         Call intg ("CORE_HOLE_PARTNER            ", core_hole_partner, &
        & df_core_hole_partner)
         Call intg ("I_CORE_HOLE                  ", i_core_hole, &
        & df_i_core_hole)
         Call real ("X                            ", x, df_x)
         Call real ("Y                            ", y, df_y)
         Call real ("Z                            ", z, df_z)
         Call stop ()
      End Subroutine dipole_write


      Subroutine dipole_allocate (offdiagonals)
         !
         ! Purpose:  allocates   dipole_integrals(:,:,:).   Called  by
         ! int_send_dipole_setup().
         !
         Implicit None
         Logical, Intent (In) :: offdiagonals
         !** End of interface *****************************************

         Integer (Kind=i4_kind) :: i_ip1, i_ip2, i_ir1, i_ir2, i_pa1, &
        & i_pa2, n_ip, n_ir, i_xyz, status, dim_ir1, dim_ir2

         dipole_offdiagonals_calculated = offdiagonals

         If (options_spin_orbit) Then
            !
            ! SPIN ORBIT
            !
            n_ip = symmetry_data_n_ip_proj ()
            n_ir = symmetry_data_n_proj_irreps ()
         Else
            n_ip = symmetry_data_n_ip ()
            n_ir = symmetry_data_n_irreps ()
         End If

         ! Default initialization  for the %use component  of the type
         ! makes sure all of these are "unused" as of now:
         Allocate (dipole_integrals(3, n_ip, n_ip), Stat=status)
         ASSERT(status==0)

         If (options_spin_orbit) Then
            ! allocate diagonal elements
            i_ip1 = 1
            irrep_so: Do i_ir1 = 1, n_ir
               dim_ir1 = symmetry_data_dimension_proj (i_ir1)
               partner_so: Do i_pa1 = 1, symmetry_data_n_partners_proj (i_ir1)
                  Do i_xyz = 1, 3
                     Allocate (dipole_integrals(i_xyz, i_ip1, i_ip1) &
                          % diagonal(dim_ir1 * (dim_ir1 + 1) / 2), &
                          dipole_integrals(i_xyz, i_ip1, i_ip1) &
                          % diagonal_imag(dim_ir1 * (dim_ir1 + 1) / 2), Stat=status)
                     ASSERT(status==0)
                     dipole_integrals(i_xyz, i_ip1, i_ip1) % use = DIPOLE_DIAGONAL
                  End Do
                  i_ip1 = i_ip1 + 1
               End Do partner_so
            End Do irrep_so
         Else
            ! allocate diagonal elements
            i_ip1 = 0
            irrep: Do i_ir1 = 1, n_ir
               dim_ir1 = symmetry_data_dimension (i_ir1)
               partner: Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
                  i_ip1 = i_ip1 + 1
                  Do i_xyz = 1, 3
                     Allocate (dipole_integrals(i_xyz, i_ip1, i_ip1) &
                          % diagonal(dim_ir1 * (dim_ir1 + 1) / 2), Stat=status)
                     ASSERT(status==0)
                     dipole_integrals(i_xyz, i_ip1, i_ip1) % use = DIPOLE_DIAGONAL
                  End Do
               End Do partner
            End Do irrep
         End If


         If (options_spin_orbit) Then
            !
            ! SPIN ORBIT
            !
            ! if required, allocate the appropriate offdiagonal blocks
            off: If (offdiagonals) Then
               i_ip1 = 1
               irrep1_so: Do i_ir1 = 1, n_ir
                  dim_ir1 = symmetry_data_dimension_proj (i_ir1)
                  partner1_so: Do i_pa1 = 1, &
                 & symmetry_data_n_partners_proj (i_ir1)
                     i_ip2 = 1
                     irrep2_so: Do i_ir2 = 1, n_ir
                        dim_ir2 = symmetry_data_dimension_proj (i_ir2)
                        partner2_so: Do i_pa2 = 1, &
                       & symmetry_data_n_partners_proj (i_ir2)
                           xyz_so: Do i_xyz = 1, 3
!!$                         if(i_ir1.eq.i_ir2 .and. i_ip2 .lt. i_ip1 ) then
                              If (symmetry_data_pdipoles_exist(i_ir2, &
                             & i_ir1, i_xyz) .And. i_ip2 .Lt. i_ip1) &
                             & Then

                                 Allocate (dipole_integrals(i_xyz, i_ip2, i_ip1) &
                                      % offdiagonal(dim_ir2, dim_ir1), &
                                      dipole_integrals(i_xyz, i_ip2, i_ip1) &
                                      % offdiagonal_imag(dim_ir2, dim_ir1), Stat=status)
                                 ASSERT(status==0)
                                 dipole_integrals(i_xyz, i_ip2, i_ip1) % use = DIPOLE_OFFDIAGONAL
                              End If
                           End Do xyz_so
                           i_ip2 = i_ip2 + 1
                        End Do partner2_so
                     End Do irrep2_so
                     i_ip1 = i_ip1 + 1
                  End Do partner1_so
               End Do irrep1_so
            End If off
         Else ! not_option_spinorbit
            ! if required, allocate the appropriate offdiagonal blocks
            If (offdiagonals) Then
               i_ip1 = 0
               irrep1: Do i_ir1 = 1, n_ir
                  dim_ir1 = symmetry_data_dimension (i_ir1)
                  partner1: Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
                     i_ip1 = i_ip1 + 1
                     i_ip2 = 0
                     irrep2: Do i_ir2 = 1, n_ir
                        dim_ir2 = symmetry_data_dimension (i_ir2)
                        partner2: Do i_pa2 = 1, symmetry_data_n_partners (i_ir2)
                           i_ip2 = i_ip2 + 1
                           xyz: Do i_xyz = 1, 3
                              If (symmetry_data_dipoles_exist(i_ir2, i_ir1, i_xyz) &
                                   .And. i_ip2 .Lt. i_ip1) Then
                                 Allocate (dipole_integrals(i_xyz, i_ip2, i_ip1) &
                                      % offdiagonal(dim_ir2, dim_ir1), Stat=status)
                                 ASSERT(status==0)
                                 dipole_integrals(i_xyz, i_ip2, i_ip1) % use = DIPOLE_OFFDIAGONAL
                              End If
                           End Do xyz
                        End Do partner2
                     End Do irrep2
                  End Do partner1
               End Do irrep1
            End If
         End If
      End Subroutine dipole_allocate


  !*************************************************************
      Logical Function dipole_xes_spectra ()
    !  Purpose: give acces to public variable make_xes_spectra
    !** End of interface *****************************************
         dipole_xes_spectra = make_xes_spectra
      End Function dipole_xes_spectra
  !*************************************************************
      Subroutine dipole_free
         !  Purpose: deallocates dipole_integrals(:,:,:)
         !** End of interface *****************************************
         Implicit None
         !------------ Declaration of local variables -----------------
         Integer (Kind=i4_kind) :: status
         !------------ Executable code --------------------------------

         If (allocated(dipole_integrals)) Then
            !
            ! Recursive deallocation of struct with allocatable components:
            !
            Deallocate (dipole_integrals, Stat=status)
            ASSERT(status==0)
         End If

         If (allocated(dipole_moments)) Then
            !
            ! Recursive deallocation of struct with allocatable components:
            !
            Deallocate (dipole_moments, dipole_total_irrep, Stat=status)
            ASSERT(status==0)
         End If
      End Subroutine dipole_free
  !*************************************************************


      Subroutine dipole_calculate ()
        !
        ! Calculates the dipole moments of all orbitals plus the total
        ! dipole moment and the dipole  moments for each Irrep and the
        ! nuclear dipole moment.  The  values are all stored in public
        ! variables of this module,
        !
        ! Runs  on  all  workers,  but  check the  body,  some  module
        ! variables are allocated/set only on master.
        !
        use comm, only: comm_rank, comm_bcast
        use solv_electrostat_module, only: surface_charge_moments
        Implicit None
        !** End of interface *****************************************

         Integer (Kind=i4_kind) :: i_ir, n_ir, i_pa, n_pa, i_ip, n_ip, &
        & i_spin, n_spin, i_orb, n_orb, i_bas1, i_bas2, i_bas12, i_xyz, &
        & status
         Real (Kind=r8_kind), Pointer :: int_x (:), int_y (:), int_z &
        & (:), eigenvector (:, :), occupation (:), int_x_imag (:), &
        & int_y_imag (:), int_z_imag (:)
         Real (Kind=r8_kind), Pointer :: eigenvector_real (:, :), &
        & eigenvector_imag (:, :)
         Real (Kind=r8_kind) :: coeff, coeff_real, coeff_imag
         Logical :: allocate_necessary
         Type (dipole_type), Pointer :: dm


         allocate_necessary = .Not. allocated (dipole_moments)
         n_spin = symmetry_data_n_spin ()

         If (options_spin_orbit) Then
            n_ip = symmetry_data_n_ip_proj ()
            n_ir = symmetry_data_n_proj_irreps ()
         Else
            n_ip = symmetry_data_n_ip ()
            n_ir = symmetry_data_n_irreps ()
         End If

         ! Set global module variables:
         dipole_nuclear = dipole_nuclear_calculate()

         ! PCM contribution:
         call surface_charge_moments (surface_charge, dipole_surface_charge)

         ! Global module variable:
         If (.not. allocated (dipole_total_irrep)) Then
            Allocate (dipole_total_irrep(3, n_ir), Stat=status)
            ASSERT (status==0)
         End If

         !
         ! Do no work ---  get results, lazy slaves. Seriousely, there
         ! is SIGSEGV lurking somewhere ...
         !
         if (comm_rank() /= 0) goto 999

         ! This will only be available on rank-0:
         If (allocate_necessary) Then
            Allocate (dipole_moments(n_spin, n_ip), Stat=status)
            ASSERT (status==0)
         End If

         dipole_total_irrep = 0.0
         dipole_total_spin = 0.0

         i_ip = 1
         irrep: Do i_ir = 1, n_ir
            If (options_spin_orbit) Then
               n_orb = symmetry_data_dimension_proj (i_ir)
               n_pa = symmetry_data_n_partners_proj (i_ir)
            Else
               n_orb = symmetry_data_dimension (i_ir)
               n_pa = symmetry_data_n_partners (i_ir)
            End If
            partner: Do i_pa = 1, n_pa

               int_x => dipole_integrals(1, i_ip, i_ip)%diagonal
               int_y => dipole_integrals(2, i_ip, i_ip)%diagonal
               int_z => dipole_integrals(3, i_ip, i_ip)%diagonal

               If (options_spin_orbit) Then
                  int_x_imag => dipole_integrals(1, i_ip, &
                 & i_ip)%diagonal_imag
                  int_y_imag => dipole_integrals(2, i_ip, &
                 & i_ip)%diagonal_imag
                  int_z_imag => dipole_integrals(3, i_ip, &
                 & i_ip)%diagonal_imag
               End If
               spin: Do i_spin = 1, n_spin

                  If (options_spin_orbit) Then
                     dm => dipole_moments (1, i_ip)
                     eigenvector_real => eigvec_real(i_ir)%m(:, :)
                     eigenvector_imag => eigvec_imag(i_ir)%m(:, :)
                     occupation => occ_num(i_ir)%m(:, 1)
                  Else
                     dm => dipole_moments (i_spin, i_ip)
                     eigenvector => eigvec(i_ir)%m(:, :, i_spin)
                     occupation => occ_num(i_ir)%m(:, i_spin)
                  End If

                  If (allocate_necessary) Then
                     Allocate (dm%occupied(3, n_orb), dm%real(3, &
                    & n_orb), Stat=status)
                     If (status .Ne. 0) Call error_handler ("dipole_cal&
                    &culate: allocate of dipole_moments%dipoles_xxx fai&
                    &led")
                  End If

                  dm%occupied = 0.0_r8_kind

                  If (options_spin_orbit) Then
                     !
                     ! SPIN ORBIT
                     !
                     orbitals_so: Do i_orb = 1, n_orb

                        i_bas12 = 1
                        Do i_bas1 = 1, n_orb
                           Do i_bas2 = 1, i_bas1 - 1
                              coeff_real = (eigenvector_real(i_bas1, &
                             & i_orb)*eigenvector_real(i_bas2, &
                             & i_orb)+eigenvector_imag(i_bas1, &
                             & i_orb)*eigenvector_imag(i_bas2, i_orb)) &
                             & * 2.0_r8_kind
                              coeff_imag = (-eigenvector_real(i_bas1, &
                             & i_orb)*eigenvector_imag(i_bas2, &
                             & i_orb)+eigenvector_imag(i_bas1, &
                             & i_orb)*eigenvector_real(i_bas2, i_orb)) &
                             & * 2.0_r8_kind
                              dm%occupied (1, i_orb) = dm%occupied(1, &
                             & i_orb) + int_x (i_bas12) * coeff_real - &
                             & int_x_imag (i_bas12) * coeff_imag
                              dm%occupied (2, i_orb) = dm%occupied(2, &
                             & i_orb) + int_y (i_bas12) * coeff_real - &
                             & int_y_imag (i_bas12) * coeff_imag
                              dm%occupied (3, i_orb) = dm%occupied(3, &
                             & i_orb) + int_z (i_bas12) * coeff_real - &
                             & int_z_imag (i_bas12) * coeff_imag
                              i_bas12 = i_bas12 + 1
                           End Do
                           coeff_real = (eigenvector_real(i_bas1, &
                          & i_orb)*eigenvector_real(i_bas2, &
                          & i_orb)+eigenvector_imag(i_bas1, &
                          & i_orb)*eigenvector_imag(i_bas2, i_orb))
                           coeff_imag = (-eigenvector_real(i_bas1, &
                          & i_orb)*eigenvector_imag(i_bas2, &
                          & i_orb)+eigenvector_imag(i_bas1, &
                          & i_orb)*eigenvector_real(i_bas2, i_orb))
                           dm%occupied (1, i_orb) = dm%occupied(1, &
                          & i_orb) + int_x (i_bas12) * coeff_real - &
                          & int_x_imag (i_bas12) * coeff_imag
                           dm%occupied (2, i_orb) = dm%occupied(2, &
                          & i_orb) + int_y (i_bas12) * coeff_real - &
                          & int_y_imag (i_bas12) * coeff_imag
                           dm%occupied (3, i_orb) = dm%occupied(3, &
                          & i_orb) + int_z (i_bas12) * coeff_real - &
                          & int_z_imag (i_bas12) * coeff_imag
                           i_bas12 = i_bas12 + 1
                        End Do

                     End Do orbitals_so
                  Else
                     orbitals: Do i_orb = 1, n_orb

                        i_bas12 = 1
                        Do i_bas1 = 1, n_orb
                           Do i_bas2 = 1, i_bas1 - 1
                              coeff = eigenvector (i_bas1, i_orb) * &
                             & eigenvector (i_bas2, i_orb) * &
                             & 2.0_r8_kind
                              dm%occupied (1, i_orb) = dm%occupied(1, &
                             & i_orb) + int_x (i_bas12) * coeff
                              dm%occupied (2, i_orb) = dm%occupied(2, &
                             & i_orb) + int_y (i_bas12) * coeff
                              dm%occupied (3, i_orb) = dm%occupied(3, &
                             & i_orb) + int_z (i_bas12) * coeff
                              i_bas12 = i_bas12 + 1
                           End Do
                           coeff = eigenvector (i_bas1, i_orb) * &
                          & eigenvector (i_bas1, i_orb)
                           dm%occupied (1, i_orb) = dm%occupied(1, &
                          & i_orb) + int_x (i_bas12) * coeff
                           dm%occupied (2, i_orb) = dm%occupied(2, &
                          & i_orb) + int_y (i_bas12) * coeff
                           dm%occupied (3, i_orb) = dm%occupied(3, &
                          & i_orb) + int_z (i_bas12) * coeff
                           i_bas12 = i_bas12 + 1
                        End Do

                     End Do orbitals
                  End If

                  Do i_xyz = 1, 3
                     dm%real (i_xyz, :) = dm%occupied(i_xyz, :) * &
                    & occupation / n_pa
                     dm%total (i_xyz) = sum (dm%real(i_xyz, :))
                     dipole_total_irrep (i_xyz, i_ir) = &
                    & dipole_total_irrep (i_xyz, i_ir) + &
                    & dm%total(i_xyz)
                     dipole_total_spin (i_xyz, i_spin) = &
                    & dipole_total_spin (i_xyz, i_spin) + &
                    & dm%total(i_xyz)
                  End Do

               End Do spin

               i_ip = i_ip + 1
            End Do partner

         End Do irrep

999      continue
         !
         ! Make these results known everywhere:
         !
         call comm_bcast (dipole_total_spin)
         call comm_bcast (dipole_total_irrep)

         ! Another module global var:
         dipole_total = dipole_nuclear - &
              (dipole_total_spin (:, 1) + dipole_total_spin (:, 2))
      End Subroutine dipole_calculate


      function dipole_nuclear_calculate() result(dipole)
        !
        ! Purpose: calculates  nuclear dipole moment  stored in public
        ! variable dipole_nuclear.  Both real atoms  and point charges
        ! are taken into account.
        !
        use unique_atom_module, only: unique_atom_type, unique_atoms
        use unique_atom_methods, only: core_charge => unique_atom_core_charge
        use pointcharge_module, only: pointcharge_array, pointcharge_N
        use datatype, only: pointcharge_type
        implicit none
        real(r8_kind) :: dipole(3) ! result
        !** End of interface *****************************************

        integer (i4_kind) :: i_ua, i_ea, i_pc, i_ec
        type (unique_atom_type), pointer :: ua
        type (pointcharge_type), pointer :: pc

        dipole = 0.0
        do i_ua = 1, size (unique_atoms)
           ua => unique_atoms (i_ua)
           do i_ea = 1, ua % N_equal_atoms
              dipole = dipole + core_charge (ua) * ua % position(:, i_ea)
           enddo
        enddo

        do i_pc = 1, pointcharge_N
           pc => pointcharge_array (i_pc)
           do i_ec = 1, pc % N_equal_charges
              dipole = dipole + pc % z * pc % position(:, i_ec)
           enddo
        enddo
      end function dipole_nuclear_calculate


      Subroutine dipole_print (iounit, detailed)
        !
        ! Prints  the dipole moments  of all  orbitals plus  the total
        ! dipole  moment and  the  dipole moments  for  each Irrep  to
        ! iounit.  dipole_calculate must be called before
        !
        ! This should run on master only. See dipole_calculate() which
        ! allocates  and  computes  orbital-resolved  "dipole_moments"
        ! only on rank-0.
        !
        use constants, only: angstrom
        use operations_module, only: operations_solvation_effect
        Implicit None
        Integer (i4_kind), Intent (In) :: iounit
        Logical, Intent (In) :: detailed
        !** End of interface *****************************************

         Integer (i4_kind) :: i_ir, i_pa, i_ip, i_spin, n_spin, i_orb, n_orb
         Logical :: open_shell
         Real (r8_kind), Pointer :: eigenvalue (:), occupation (:)
         Type (dipole_type), Pointer :: dm
         Integer (i4_kind) :: n_irreps, i
         Integer (i4_kind), Allocatable :: n_partners (:), dim_irrep (:)
         Character (len=12), Allocatable :: irrepname (:)

         ! Slaves will complain right here:
         If (.Not. allocated (dipole_moments)) Call error_handler ("dip&
        &ole_print: dipole moments have not been calculated before")

         n_spin = symmetry_data_n_spin ()
         If (options_spin_orbit) Then
            n_irreps = symmetry_data_n_proj_irreps ()
            Allocate (n_partners(n_irreps), dim_irrep(n_irreps), &
           & irrepname(n_irreps))
            Do i = 1, n_irreps
               n_partners (i) = symmetry_data_n_partners_proj (i)
               dim_irrep (i) = symmetry_data_dimension_proj (i)
               irrepname (i) = symmetry_data_irrepname_proj (i)
            End Do
         Else
            n_irreps = symmetry_data_n_irreps ()
            Allocate (n_partners(n_irreps), dim_irrep(n_irreps), &
           & irrepname(n_irreps))
            Do i = 1, n_irreps
               n_partners (i) = symmetry_data_n_partners (i)
               dim_irrep (i) = symmetry_data_dimension (i)
               irrepname (i) = symmetry_data_irrepname (i)
            End Do
         End If
         open_shell = n_spin .Eq. 2

         ! Haeder
         Write (iounit, Fmt=*)
         Write (iounit, Fmt=*)
         Write (iounit, Fmt=*)
         Write (iounit, Fmt='(30X,"#####################")')
         Write (iounit, Fmt='(30X,"##  DIPOLE MOMENTS  ##")')
         Write (iounit, Fmt='(30X,"#####################")')
         Write (iounit, Fmt=*)
         Write (iounit, Fmt='(49X,"X",16X,"Y",16X,"Z")')

1001     format (A5, " Dipole Moment (A.U.)", 11X, 3F17.7)
1002     format (A5, " Dipole Moment (C*ANGSTROEM)", 4X, 3F17.7)
1003     format (A5, " Dipole Moment (C*M*10**-12)", 4X, 3F17.7)

         ! Total
         write (iounit, 1001) "Total", dipole_total
         write (iounit, 1002) "Total", dipole_total / angstrom
         write (iounit, 1003) "Total", dipole_total * 100 / angstrom

         ! Total using ES-format instead of F-format
         If (detailed) Then
            Write (iounit, Fmt=*)
            Write (iounit, Fmt=*)
            Write (iounit, Fmt='("Total Dipole Moment (A.U.)",11X,3ES17.10)') dipole_total
            Write (iounit, Fmt='("Total Dipole Moment (C*ANGSTROEM)",4X,3ES17.10)')&
           & dipole_total / angstrom
            Write (iounit, Fmt='("Total Dipole Moment (C*M*10**-12)",4X,3ES17.10)') &
           & dipole_total * 100 / angstrom
         End If

         ! Nuclear
         Write (iounit, Fmt=*)
         Write (iounit, Fmt='("Nuclear Dipole Moment",16X,3F17.7)') dipole_nuclear

         ! Total electronic
         Write (iounit, Fmt=*)
         Write (iounit, Fmt='("Total Electronic Dipole Moment",7X,3F17.7)') dipole_total_spin (:, 1) + &
        & dipole_total_spin (:, 2)

         ! Apparent surface charge. Dont print zeros if no solvation:
         if (operations_solvation_effect) then
            write (iounit, Fmt=*)
            write (iounit, 1001) "ASC", dipole_surface_charge
            write (iounit, 1002) "ASC", dipole_surface_charge / angstrom
            write (iounit, 1003) "ASC", dipole_surface_charge * 100 / angstrom
            write (iounit, Fmt=*)
            write (iounit, *) " ASC total charge", surface_charge
         endif

         ! Total Spin
         If (open_shell) Then
            Write (iounit, Fmt=*)
            Do i_spin = 1, 2
               Write (iounit, Fmt='("Total (all Irreps) for Spin",I3,7X,3F17.7)')&
              & i_spin, dipole_total_spin (:, &
              & i_spin)
            End Do
         End If

         ! For single Irreps
         i_ip = 1
         Do i_ir = 1, n_irreps
            Write (iounit, Fmt=*)
            Write (iounit, Fmt='("IRREP No",I3,A4,"  Total",15X,3F17.7)')&
           & i_ir, trim (irrepname(i_ir)), &
           & dipole_total_irrep (:, i_ir)
            Do i_pa = 1, n_partners (i_ir)
               Do i_spin = 1, n_spin
                  If (n_partners(i_ir) .Gt. 1) Then
                     If (open_shell) Then
                        Write (iounit, Fmt='(14X,"Partner",I3,"  Spin",I3,4X,3F17.7)')&
                       & i_pa, i_spin, &
                       & dipole_moments(i_spin, i_ip)%total
                     Else
                        Write (iounit, Fmt='(14X,"Partner",I3,13X,3F17.7)') i_pa, &
                       & dipole_moments(i_spin, i_ip)%total
                     End If
                  Else
                     If (open_shell) Then
                        Write (iounit, Fmt='(14X,"Spin",I3,16X,3F17.7)') i_spin, dipole_moments(i_spin, &
                       & i_ip)%total
                     End If
                  End If
               End Do
               i_ip = i_ip + 1
            End Do
         End Do

         Write (iounit, Fmt=*)
         Write (iounit, Fmt=*)

         ! for single orbitals
         orbital_output: If (detailed) Then

            Write (iounit, Fmt=*)
            Write (iounit, Fmt=*)
            Write (iounit, Fmt='(35X,"## DIPOLE MOMENTS OF SINGLE ORBITALS  ##")')

            i_ip = 1
            irrep: Do i_ir = 1, n_irreps
               n_orb = dim_irrep (i_ir)
               partner: Do i_pa = 1, n_partners (i_ir)
                  spin: Do i_spin = 1, n_spin

                     dm => dipole_moments (i_spin, i_ip)
                     eigenvalue => eigval(i_ir)%m(:, i_spin)
                     occupation => occ_num(i_ir)%m(:, i_spin)

                     Write (iounit, Fmt=*)
                     Write (iounit, Fmt=*)
                     If (open_shell) Then
                        Write (iounit, Fmt='("IRREP No",I3,A4,"  Partner",I3,"  Spin",I3)') i_ir, trim &
                       & (irrepname(i_ir)), i_pa, i_spin
                     Else
                        Write (iounit, Fmt='("IRREP No",I3,A4,"  Partner",I3)') i_ir, trim &
                       & (irrepname(i_ir)), i_pa
                     End If
                     Write (iounit, Fmt=*)
                     Write (iounit,&
                    & Fmt='(10X,"(A.U.)",32X,"With assumed Occupation 1.0",20X,"With real Occupation")')
                     Write (iounit,&
                    &Fmt='(" Orbital   Eigenvalue      Occupation",9X,"X",14X,"Y",14X,"Z",14X,"X",14X,"Y",14X,"Z")')

                     orbitals: Do i_orb = 1, n_orb
                        Write (iounit, Fmt='(I5,3X,F13.5,7F15.5)') i_orb, eigenvalue (i_orb) * &
                       & 27.211652_r8_kind, occupation (i_orb), &
                       & dm%occupied(:, i_orb), dm%real(:, i_orb)
                     End Do orbitals

                  End Do spin
                  i_ip = i_ip + 1
               End Do partner
            End Do irrep

            Write (iounit, Fmt=*)
            Write (iounit, Fmt=*)
         End If orbital_output

         Deallocate (dim_irrep, n_partners, irrepname)
      End Subroutine dipole_print


      Subroutine dipole_write_simol ()
    !  Purpose: writes total dipole moment in a.u. in formatted way
    !  to file "dipmom.sav" used by simole to calculate fequencies.
    !** End of interface *****************************************
         Use iounitadmin_module
         Use filename_module
         Integer (Kind=i4_kind) :: iostat, iounit
         If ( .Not. allocated(dipole_moments)) Call error_handler ("dip&
        &ole_write_simol: dipole moments have not been calculated befor&
        &e")
         iounit = openget_iounit (file=trim(data_dir)//'/dipmom.sav', &
        & form='FORMATTED', action='WRITE', position='APPEND',&
        &  status='UNKNOWN')
         Write (iounit, IoStat=IoStat, Fmt='(3f15.8)') dipole_total
         If (iostat .Ne. 0) Call error_handler ("dipole_write_simol: wr&
        &ite failed")
         Call returnclose_iounit (iounit, status='KEEP')
      End Subroutine dipole_write_simol

      subroutine dipole_write_optimizer(iounit)
         !  Purpose: appends total dipole moment in a.u. in formatted way
         !  to gxfile used by optimizer to calculate fequencies.
         integer(i4_kind), intent(in) :: iounit
         ! *** end of interface ***

         integer(i4_kind) :: iostat

         write (iounit, iostat=iostat, fmt='(3f15.8)') dipole_total
         if (iostat/=0) &
          call error_handler ("dipole_write_optimizer: write failed")
      end subroutine dipole_write_optimizer

      Subroutine dipole_transitionmoment_f (iounit)
    !  Purpose: calculates and prints all transition dipole moments
    !  in a formatted way to iounit. Includes spin-orbit option
         Implicit None
    !------------ Declaration of formal parameters ---------------
         Integer (Kind=i4_kind), Intent (In) :: iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
         Integer (Kind=i4_kind) :: n_ir, i_ir1, i_pa1, i_ip1, i_ir2, &
        & i_pa2, i_ip2, min_orb2, i_spin, n_spin, i_orb1, n_orb1, &
        & i_orb2, n_orb2, iostat, n_out, i
         Logical :: open_shell
         Real (Kind=r8_kind), Pointer :: eigenvector1 (:), eigenvector2 &
        & (:)
         Real (Kind=r8_kind), Pointer :: eigenvector1_real (:), &
        & eigenvector1_imag (:)
         Real (Kind=r8_kind), Pointer :: eigenvector2_real (:), &
        & eigenvector2_imag (:)
         Real (Kind=r8_kind) :: e1, e2
         Integer (Kind=i4_kind), Allocatable :: n_partners (:), &
        & dim_irrep (:)
         Character(len=12), Allocatable :: irrepname (:)
    !------------ Executable code --------------------------------
!
         n_spin = symmetry_data_n_spin ()
         If (options_spin_orbit) Then
            n_ir = symmetry_data_n_proj_irreps ()
            Allocate (n_partners(n_ir), dim_irrep(n_ir), &
           & irrepname(n_ir))
            Do i = 1, n_ir
               n_partners (i) = symmetry_data_n_partners_proj (i)
               dim_irrep (i) = symmetry_data_dimension_proj (i)
               irrepname (i) = symmetry_data_irrepname_proj (i)
            End Do
         Else
            n_ir = symmetry_data_n_irreps ()
            Allocate (n_partners(n_ir), dim_irrep(n_ir), &
           & irrepname(n_ir))
            Do i = 1, n_ir
               n_partners (i) = symmetry_data_n_partners (i)
               dim_irrep (i) = symmetry_data_dimension (i)
               irrepname (i) = symmetry_data_irrepname (i)
            End Do
         End If
         open_shell = n_spin .Eq. 2
!
    ! Haeder
         Write (iounit, IoStat=IoStat, Fmt=*)
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of blank failed")
         Write (iounit, IoStat=IoStat, Fmt=*)
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of blank failed")
         Write (iounit, IoStat=IoStat, Fmt=*)
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of blank failed")
         Write (iounit, IoStat=IoStat,&
        & Fmt='(30X,"################################")')
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of haeder failed")
         Write (iounit, IoStat=IoStat,&
        &  Fmt='(30X,"##  DIPOLE TRANSITION MOMENTS  ##")')
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of haeder failed")
         Write (iounit, IoStat=IoStat,&
        &  Fmt='(30X,"################################")')
         If (iostat .Ne. 0) Call error_handler (" dipole_transitionmome&
        &nt_f: write of haeder failed")
         Write (iounit, IoStat=IoStat, Fmt=*)
!
!
         valence_excitations: If (show_valence_excitations) Then
!
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
            Write (iounit, IoStat=IoStat,&
           & Fmt='(32X,"##  VALENCE EXCITATIONS  ##")')
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of haeder failed")
            Write (iounit, IoStat=IoStat,&
           &  Fmt='(25X,"Excitation Energies up to ",F10.3," eV")')&
           & max_valence_excitation_energy
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of haeder failed")
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
!
            i_ip1 = 0
            irrep1: Do i_ir1 = 1, n_ir
               n_orb1 = dim_irrep (i_ir1)
               partner1: Do i_pa1 = 1, n_partners (i_ir1)
                  i_ip1 = i_ip1 + 1
                  spin: Do i_spin = 1, n_spin
!
                     n_out = 0
                     Do i_orb1 = 1, n_occo (i_spin, i_ir1)
                        e1 = eigval(i_ir1)%m(i_orb1, i_spin) * &
                       & 27.211652_r8_kind
                        Do i_ir2 = 1, n_ir
                           n_orb2 = dim_irrep (i_ir2)
                           min_orb2 = n_occo (i_spin, i_ir2) + 1
                           Do i_pa2 = 1, n_partners (i_ir2)
                              Do i_orb2 = min_orb2, n_orb2
                                 e2 = eigval(i_ir2)%m(i_orb2, i_spin) * &
                                & 27.211652_r8_kind
                                 If (max_valence_excitation_energy .Ge. &
                                & e2-e1) n_out = n_out + 1
                              End Do
                           End Do
                        End Do
                     End Do
!
                     If (n_out .Eq. 0) Cycle spin
!
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     If (open_shell) Then
                        Write (iounit, IoStat=IoStat,&
                       & Fmt='("DIPOLE TRANSITION MOMENTS from IRREP No",I3,A4,"  Partner",I3,"  for Spin",I3,"  (A.U.)")')&
                       & i_ir1, trim (irrepname(i_ir1)), i_pa1, i_spin
                     Else
                        Write (iounit, IoStat=IoStat,&
                       & Fmt='("DIPOLE TRANSITION MOMENTS from IRREP No",I3,A4,"  Partner",I3,"  (A.U.)")')&
                       i_ir1, trim &
                       & (irrepname(i_ir1)), i_pa1
                     End If
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital dipole moments &
                    &haeder failed")
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     Write (iounit, IoStat=IoStat,&
                    & Fmt='("From",17X,"To",26X,"Excitation Energy               Dipole Transition Moment (A.U.)")')
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital legend failed")
                     Write (iounit, IoStat=IoStat,&
                    &  Fmt='("Orbital Eigenvalue   Irrep Partner Orbital Eigenvalue   delta E",10X,"ABS",15X,"X",17X,"Y",17X,"Z")')
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital legend failed")
!
                     orbitals1: Do i_orb1 = 1, n_occo (i_spin, i_ir1)
                        e1 = eigval(i_ir1)%m(i_orb1, i_spin) * &
                       & 27.211652_r8_kind
                        If (options_spin_orbit) Then
                           eigenvector1_real => eigvec_real(i_ir1)%m(:, &
                          & i_orb1)
                           eigenvector1_imag => eigvec_imag(i_ir1)%m(:, &
                          & i_orb1)
                        Else
                           eigenvector1 => eigvec(i_ir1)%m(:, i_orb1, &
                          & i_spin)
                        End If
!
                        i_ip2 = 0
                        irrep2: Do i_ir2 = 1, n_ir
                           n_orb2 = dim_irrep (i_ir2)
                           min_orb2 = n_occo (i_spin, i_ir2) + 1
                           partner2: Do i_pa2 = 1, n_partners (i_ir2)
                              i_ip2 = i_ip2 + 1
!
                              orbitals2: Do i_orb2 = min_orb2, n_orb2
                                 e2 = eigval(i_ir2)%m(i_orb2, i_spin) * &
                                & 27.211652_r8_kind
                                 If (options_spin_orbit) Then
                                    eigenvector2_real => &
                                   & eigvec_real(i_ir2)%m(:, i_orb2)
                                    eigenvector2_imag => &
                                   & eigvec_imag(i_ir2)%m(:, i_orb2)
                                 Else
                                    eigenvector2 => eigvec(i_ir2)%m(:, &
                                   & i_orb2, i_spin)
                                 End If
!
                                 If (max_valence_excitation_energy .Ge. &
                                & e2-e1) Call calc_write_tm ()
!
                              End Do orbitals2
                           End Do partner2
!
                        End Do irrep2
!
                     End Do orbitals1
                  End Do spin
               End Do partner1
            End Do irrep1
!
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
!
         End If valence_excitations
!
         core_excitations: If (show_core_excitations) Then
!
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
            Write (iounit, IoStat=IoStat,&
           &Fmt='(32X,"##  CORE EXCITATIONS  ##")')
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of haeder failed")
            Write (iounit, IoStat=IoStat,&
           & Fmt='(25X,"Low  Levels between ",F10.3," to ",F10.3," eV")')&
           & core_excitation_from_e_min, &
           & core_excitation_from_e_max
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of haeder failed")
            Write (iounit, IoStat=IoStat,&
           & Fmt='(25X,"High Levels between ",F10.3," to ",F10.3," eV")')&
           & core_excitation_to_e_min, &
           & core_excitation_to_e_max
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of haeder failed")
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
!
            i_ip1 = 0
            core_irrep1: Do i_ir1 = 1, n_ir
               n_orb1 = dim_irrep (i_ir1)
               core_partner1: Do i_pa1 = 1, n_partners (i_ir1)
                  i_ip1 = i_ip1 + 1
                  core_spin: Do i_spin = 1, n_spin
!
                     n_out = 0
                     Do i_orb1 = 1, n_orb1
                        e1 = eigval(i_ir1)%m(i_orb1, i_spin) * &
                       & 27.211652_r8_kind
                        If (e1 .Lt. core_excitation_from_e_min) Cycle
                        If (e1 .Gt. core_excitation_from_e_max) Exit
                        Do i_ir2 = 1, n_ir
                           n_orb2 = dim_irrep (i_ir2)
                           Do i_pa2 = 1, n_partners (i_ir2)
                              Do i_orb2 = 1, n_orb2
                                 e2 = eigval(i_ir2)%m(i_orb2, i_spin) * &
                                & 27.211652_r8_kind
                                 If (e2 .Lt. core_excitation_to_e_min) &
                                & Cycle
                                 If (e2 .Gt. core_excitation_to_e_max) &
                                & Exit
                                 n_out = n_out + 1
                              End Do
                           End Do
                        End Do
                     End Do
!
                     If (n_out .Eq. 0) Cycle core_spin
!
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     If (open_shell) Then
                        Write (iounit, IoStat=IoStat,&
                       & Fmt='("DIPOLE TRANSITION MOMENTS from IRREP No",I3,A4,"  Partner",I3,"  for Spin",I3,"  (A.U.)")')&
                       & i_ir1, trim (irrepname(i_ir1)), i_pa1, i_spin
                     Else
                        Write (iounit, IoStat=IoStat,&
                       & Fmt='("DIPOLE TRANSITION MOMENTS from IRREP No",I3,A4,"  Partner",I3,"  (A.U.)")')&
                       & i_ir1, trim &
                       & (irrepname(i_ir1)), i_pa1
                     End If
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital dipole moments &
                    &haeder failed")
                     Write (iounit, IoStat=IoStat, Fmt=*)
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of blank failed")
                     Write (iounit, IoStat=IoStat,&
                    & Fmt='("From",17X,"To",26X,"Excitation Energy             Dipole Transition Moment (A.U.)")')
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital legend failed")
                     Write (iounit, IoStat=IoStat,&
                    & Fmt='("Orbital Eigen value   Irrep Partner Orbital Eigenvalue   delta E",10X,"ABS",15X,"X",17X,"Y",17X,"Z")')
                     If (iostat .Ne. 0) Call error_handler (" dipole_tr&
                    &ansitionmoment_f: write of orbital legend failed")
!
                     core_orbitals1: Do i_orb1 = 1, n_orb1
                        e1 = eigval(i_ir1)%m(i_orb1, i_spin) * &
                       & 27.211652_r8_kind
!
                        If (e1 .Lt. core_excitation_from_e_min) Cycle &
                       & core_orbitals1
                        If (e1 .Gt. core_excitation_from_e_max) Exit &
                       & core_orbitals1
!
                        If (options_spin_orbit) Then
                           eigenvector1_real => eigvec_real(i_ir1)%m(:, &
                          & i_orb1)
                           eigenvector1_imag => eigvec_imag(i_ir1)%m(:, &
                          & i_orb1)
                        Else
                           eigenvector1 => eigvec(i_ir1)%m(:, i_orb1, &
                          & i_spin)
                        End If
!
                        i_ip2 = 0
                        core_irrep2: Do i_ir2 = 1, n_ir
                           n_orb2 = dim_irrep (i_ir2)
                           core_partner2: Do i_pa2 = 1, n_partners &
                          & (i_ir2)
                              i_ip2 = i_ip2 + 1
!
                              core_orbitals2: Do i_orb2 = 1, n_orb2
                                 e2 = eigval(i_ir2)%m(i_orb2, i_spin) * &
                                & 27.211652_r8_kind
!
                                 If (e2 .Lt. core_excitation_to_e_min) &
                                & Cycle core_orbitals2
                                 If (e2 .Gt. core_excitation_to_e_max) &
                                & Exit core_orbitals2
!
                                 If (options_spin_orbit) Then
                                    eigenvector2_real => &
                                   & eigvec_real(i_ir2)%m(:, i_orb2)
                                    eigenvector2_imag => &
                                   & eigvec_imag(i_ir2)%m(:, i_orb2)
                                 Else
                                    eigenvector2 => eigvec(i_ir2)%m(:, &
                                   & i_orb2, i_spin)
                                 End If
!
                                 Call calc_write_tm ()
!
                              End Do core_orbitals2
                           End Do core_partner2
                        End Do core_irrep2
!
                     End Do core_orbitals1
                  End Do core_spin
               End Do core_partner1
            End Do core_irrep1
!
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
            Write (iounit, IoStat=IoStat, Fmt=*)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_f: write of blank failed")
!
         End If core_excitations
!
      Contains
!
         Subroutine calc_write_tm ()!!!!
      ! Porpose: calculation and output of one transition moment between two orbitals
            Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12, i_xyz
            Real (Kind=r8_kind) :: transitionmoments (3)
            Real (Kind=r8_kind) :: transitionmoments_imag (3)
            Real (Kind=r8_kind) :: coeff_real_a, coeff_real_b, &
           & coeff_imag_a, coeff_imag_b
            Type (dipole_integral_type), Pointer :: di
            Do i_xyz = 1, 3
               If (i_ip2 .Gt. i_ip1) Then
                  di => dipole_integrals (i_xyz, i_ip1, i_ip2)
               Else
                  di => dipole_integrals (i_xyz, i_ip2, i_ip1)
               End If
               integralstorage: Select Case (di % use)
               Case (DIPOLE_DIAGONAL) integralstorage
            !triangular storage
                  i_bas12 = 1
                  transitionmoments (i_xyz) = 0.0_r8_kind
                  If (options_spin_orbit) Then
                     transitionmoments_imag (i_xyz) = 0.0_r8_kind
                     bas1_diagonal_so: Do i_bas1 = 1, n_orb1
                        bas2_diagonal_so: Do i_bas2 = 1, i_bas1 - 1
                           coeff_real_a = (eigenvector2_real(i_bas1)*&
                          & eigenvector1_real(i_bas2)+&
                          & eigenvector2_imag(i_bas1)*&
                          & eigenvector1_imag(i_bas2))
                           coeff_imag_a = (-eigenvector2_real(i_bas1)*&
                          & eigenvector1_imag(i_bas2)+&
                          & eigenvector2_imag(i_bas1)*&
                          & eigenvector1_real(i_bas2))
                           coeff_real_b = (eigenvector2_real(i_bas2)*&
                          & eigenvector1_real(i_bas1)+&
                          & eigenvector2_imag(i_bas2)*&
                          & eigenvector1_imag(i_bas1))
                           coeff_imag_b = (-eigenvector2_real(i_bas2)*&
                          & eigenvector1_imag(i_bas1)+&
                          & eigenvector2_imag(i_bas2)*&
                          & eigenvector1_real(i_bas1))
                           transitionmoments (i_xyz) = &
                          & transitionmoments (i_xyz) + &
                          & di%diagonal(i_bas12) * &
                          & (coeff_real_a+coeff_real_b)
                           transitionmoments_imag (i_xyz) = &
                          & transitionmoments_imag (i_xyz) + &
                          & di%diagonal(i_bas12) * &
                          & (coeff_imag_a+coeff_imag_b)
                           transitionmoments (i_xyz) = &
                          & transitionmoments (i_xyz) + &
                          & di%diagonal_imag(i_bas12) * &
                          & (coeff_imag_b-coeff_imag_a)
                           transitionmoments_imag (i_xyz) = &
                          & transitionmoments_imag (i_xyz) + &
                          & di%diagonal_imag(i_bas12) * &
                          & (coeff_real_a-coeff_real_b)
                           i_bas12 = i_bas12 + 1
                        End Do bas2_diagonal_so
                        coeff_real_a = (eigenvector2_real(i_bas1)*&
                       & eigenvector1_real(i_bas1)+&
                       & eigenvector2_imag(i_bas1)*&
                       & eigenvector1_imag(i_bas1))
                        coeff_imag_a = (-eigenvector2_real(i_bas1)*&
                       & eigenvector1_imag(i_bas1)+&
                       & eigenvector2_imag(i_bas1)*&
                       & eigenvector1_real(i_bas1))
                        transitionmoments (i_xyz) = transitionmoments &
                       & (i_xyz) + di%diagonal(i_bas12) * coeff_real_a
                        transitionmoments_imag (i_xyz) = &
                       & transitionmoments_imag (i_xyz) + &
                       & di%diagonal(i_bas12) * coeff_imag_a
                        i_bas12 = i_bas12 + 1
                     End Do bas1_diagonal_so
                  Else
                     bas1_diagonal: Do i_bas1 = 1, n_orb1
                        bas2_diagonal: Do i_bas2 = 1, i_bas1 - 1
                           transitionmoments (i_xyz) = &
                          & transitionmoments (i_xyz) + &
                          & (eigenvector1(i_bas1)*eigenvector2(i_bas2)+&
                          & eigenvector1(i_bas2)*eigenvector2(i_bas1)) &
                          & * di%diagonal(i_bas12)
                           i_bas12 = i_bas12 + 1
                        End Do bas2_diagonal
                        transitionmoments (i_xyz) = transitionmoments &
                       & (i_xyz) + eigenvector1 (i_bas1) * eigenvector2 &
                       & (i_bas1) * di%diagonal(i_bas12)
                        i_bas12 = i_bas12 + 1
                     End Do bas1_diagonal
                  End If
               Case (DIPOLE_OFFDIAGONAL) integralstorage
            ! full storage
                  transitionmoments (i_xyz) = 0.0_r8_kind
                  If (options_spin_orbit) Then
                     transitionmoments_imag (i_xyz) = 0.0_r8_kind
                     If (i_ip2 .Gt. i_ip1) Then
                        Do i_bas2 = 1, n_orb2
                           Do i_bas1 = 1, n_orb1
                              coeff_real_a = (eigenvector2_real(i_bas2)*&
                             & eigenvector1_real(i_bas1)+&
                             & eigenvector2_imag(i_bas2)*&
                             & eigenvector1_imag(i_bas1))
                              coeff_imag_a = - (-&
                             & eigenvector2_real(i_bas2)*&
                             & eigenvector1_imag(i_bas1)+&
                             & eigenvector2_imag(i_bas2)*&
                             & eigenvector1_real(i_bas1))
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) + &
                             & di%offdiagonal(i_bas1, i_bas2) * &
                             & coeff_real_a
                              transitionmoments_imag (i_xyz) = &
                             & transitionmoments_imag (i_xyz) + &
                             & di%offdiagonal(i_bas1, i_bas2) * &
                             & coeff_imag_a
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) - &
                             & di%offdiagonal_imag(i_bas1, i_bas2) * &
                             & coeff_imag_a
                              transitionmoments_imag (i_xyz) = &
                             & transitionmoments_imag (i_xyz) + &
                             & di%offdiagonal_imag(i_bas1, i_bas2) * &
                             & coeff_real_a
!
                           End Do
                        End Do
                     Else
                        Do i_bas1 = 1, n_orb1
                           Do i_bas2 = 1, n_orb2
                              coeff_real_a = (eigenvector2_real(i_bas2)*&
                             & eigenvector1_real(i_bas1)+&
                             & eigenvector2_imag(i_bas1)*&
                             & eigenvector1_imag(i_bas2))
                              coeff_imag_a = (-&
                             & eigenvector2_real(i_bas2)*&
                             & eigenvector1_imag(i_bas1)+&
                             & eigenvector2_imag(i_bas2)*&
                             & eigenvector1_real(i_bas1))
!
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) + &
                             & di%offdiagonal(i_bas2, i_bas1) * &
                             & coeff_real_a
                              transitionmoments_imag (i_xyz) = &
                             & transitionmoments_imag (i_xyz) + &
                             & di%offdiagonal(i_bas2, i_bas1) * &
                             & coeff_imag_a
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) + &
                             & di%offdiagonal_imag(i_bas2, i_bas1) * &
                             & coeff_imag_a
                              transitionmoments_imag (i_xyz) = &
                             & transitionmoments_imag (i_xyz) - &
                             & di%offdiagonal_imag(i_bas2, i_bas1) * &
                             & coeff_real_a
!
                           End Do
                        End Do
                     End If
                  Else
                     If (i_ip2 .Gt. i_ip1) Then
                        Do i_bas2 = 1, n_orb2
                           Do i_bas1 = 1, n_orb1
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) + eigenvector1 &
                             & (i_bas1) * eigenvector2 (i_bas2) * &
                             & di%offdiagonal(i_bas1, i_bas2)
                           End Do
                        End Do
                     Else
                        Do i_bas1 = 1, n_orb1
                           Do i_bas2 = 1, n_orb2
                              transitionmoments (i_xyz) = &
                             & transitionmoments (i_xyz) + eigenvector1 &
                             & (i_bas1) * eigenvector2 (i_bas2) * &
                             & di%offdiagonal(i_bas2, i_bas1)
                           End Do
                        End Do
                     End If
                  End If
               Case (DIPOLE_UNUSED) integralstorage
                  transitionmoments (i_xyz) = 0.0_r8_kind
               Case Default integralstorage
                  Call error_handler ("dipole_transitionmoment_f: forbi&
                 &dden case")
               End Select integralstorage
            End Do
!
            If (options_spin_orbit) Then
               Write (iounit, IoStat=IoStat,&
              & Fmt='(I5,F12.5,I6,A4,I4,I7,2F12.5,4F18.5)') i_orb1, e1, i_ir2, trim &
              & (irrepname(i_ir2)), i_pa2, i_orb2, e2, e2 - e1, Sqrt &
              & (sum(transitionmoments(:)**2+&
              & transitionmoments_imag(:)**2)), &
              & (transitionmoments_imag(i_xyz), i_xyz=1, 3)
               If (iostat .Ne. 0) Call error_handler (" dipole_transiti&
              &onmoment_f: write of orbital dipole moments failed")
            Else
               Write (iounit, IoStat=IoStat,&
              & Fmt='(I5,F12.5,I6,A4,I4,I7,2F12.5,4F18.5)') i_orb1, e1, i_ir2, trim &
              & (irrepname(i_ir2)), i_pa2, i_orb2, e2, e2 - e1, Sqrt &
              & (sum(transitionmoments(:)**2)), &
              & (transitionmoments(i_xyz), i_xyz=1, 3)
               If (iostat .Ne. 0) Call error_handler (" dipole_transiti&
              &onmoment_f: write of orbital dipole moments failed")
            End If
         End Subroutine calc_write_tm !
!
      End Subroutine dipole_transitionmoment_f
  !*************************************************************
      Subroutine dipole_transitionmoment_uf (iounit)
    !  Purpose: calculates and prints all transition dipole moments
    !  in a unformatted way to iounit.
         Implicit None
    !------------ Declaration of formal parameters ---------------
         Integer (Kind=i4_kind), Intent (In) :: iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
         Integer (Kind=i4_kind) :: n_ir, i_ir1, i_pa1, i_ip1, i_ir2, &
        & i_pa2, n_pa2, i_ip2, i_spin, n_spin, i_orb1, n_orb1, i_orb2, &
        & n_orb2, iostat
         Logical :: open_shell
         Real (Kind=r8_kind), Pointer :: eigenvector1 (:), eigenvector2 &
        & (:)
         Real (Kind=r8_kind) :: e1, e2
    !------------ Executable code --------------------------------
!
         n_spin = symmetry_data_n_spin ()
         open_shell = n_spin .Eq. 2
!
         n_ir = symmetry_data_n_irreps ()
!
         i_ip1 = 0
         irrep1: Do i_ir1 = 1, n_ir
            n_orb1 = symmetry_data_dimension (i_ir1)
            partner1: Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
               i_ip1 = i_ip1 + 1
!
               i_ip2 = 0
               irrep2: Do i_ir2 = 1, i_ir1
                  n_orb2 = symmetry_data_dimension (i_ir2)
                  If (i_ir1 .Eq. i_ir2) Then
                     n_pa2 = i_pa1
                  Else
                     n_pa2 = symmetry_data_n_partners (i_ir2)
                  End If
                  partner2: Do i_pa2 = 1, n_pa2
                     i_ip2 = i_ip2 + 1
!
                     spin: Do i_spin = 1, n_spin
!
                        orbitals1: Do i_orb1 = 1, n_orb1
                           e1 = eigval(i_ir1)%m(i_orb1, i_spin) * &
                          & 27.211652_r8_kind
                           eigenvector1 => eigvec(i_ir1)%m(:, i_orb1, &
                          & i_spin)
                           If (i_ip1 .Eq. i_ip2) n_orb2 = i_orb1
                           orbitals2: Do i_orb2 = 1, n_orb2
                              e2 = eigval(i_ir2)%m(i_orb2, i_spin) * &
                             & 27.211652_r8_kind
                              eigenvector2 => eigvec(i_ir2)%m(:, &
                             & i_orb2, i_spin)
!
                              Call calc_write_tm ()
!
                           End Do orbitals2
                        End Do orbitals1
!
                     End Do spin
!
                  End Do partner2
               End Do irrep2
            End Do partner1
         End Do irrep1
!
      Contains
!
         Subroutine calc_write_tm ()
      ! Porpose: calculation and output of one transition moment between two orbitals
            Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12, i_xyz
            Real (Kind=r8_kind) :: transitionmoments (3)
            Type (dipole_integral_type), Pointer :: di
            Do i_xyz = 1, 3
               If (i_ip2 .Gt. i_ip1) Then
                  di => dipole_integrals (i_xyz, i_ip1, i_ip2)
               Else
                  di => dipole_integrals (i_xyz, i_ip2, i_ip1)
               End If
               integralstorage: Select Case (di % use)
               Case (DIPOLE_DIAGONAL) integralstorage
            !triangular storage
                  i_bas12 = 1
                  transitionmoments (i_xyz) = 0.0_r8_kind
                  bas1_diagonal: Do i_bas1 = 1, n_orb1
                     bas2_diagonal: Do i_bas2 = 1, i_bas1 - 1
                        transitionmoments (i_xyz) = transitionmoments &
                       & (i_xyz) + &
                       & (eigenvector1(i_bas1)*eigenvector2(i_bas2)+&
                       & eigenvector1(i_bas2)*eigenvector2(i_bas1)) * &
                       & di%diagonal(i_bas12)
                        i_bas12 = i_bas12 + 1
                     End Do bas2_diagonal
                     transitionmoments (i_xyz) = transitionmoments &
                    & (i_xyz) + eigenvector1 (i_bas1) * eigenvector2 &
                    & (i_bas1) * di%diagonal(i_bas12)
                     i_bas12 = i_bas12 + 1
                  End Do bas1_diagonal
               Case (DIPOLE_OFFDIAGONAL) integralstorage
            ! full storage
                  transitionmoments (i_xyz) = 0.0_r8_kind
                  If (i_ip2 .Gt. i_ip1) Then
                     Do i_bas2 = 1, n_orb2
                        Do i_bas1 = 1, n_orb1
                           transitionmoments (i_xyz) = &
                          & transitionmoments (i_xyz) + eigenvector1 &
                          & (i_bas1) * eigenvector2 (i_bas2) * &
                          & di%offdiagonal(i_bas1, i_bas2)
                        End Do
                     End Do
                  Else
                     Do i_bas1 = 1, n_orb1
                        Do i_bas2 = 1, n_orb2
                           transitionmoments (i_xyz) = &
                          & transitionmoments (i_xyz) + eigenvector1 &
                          & (i_bas1) * eigenvector2 (i_bas2) * &
                          & di%offdiagonal(i_bas2, i_bas1)
                        End Do
                     End Do
                  End If
               Case (DIPOLE_UNUSED) integralstorage
                  transitionmoments (i_xyz) = 0.0_r8_kind
               Case Default integralstorage
                  Call error_handler ("dipole_transitionmoment_f: forbi&
                 &dden case")
               End Select integralstorage
            End Do
!
            Write (iounit, IoStat=IoStat) i_orb1, e1, i_ir2, i_pa2, &
           & i_orb2, e2, e2 - e1, (transitionmoments(i_xyz), i_xyz=1, &
           & 3)
            If (iostat .Ne. 0) Call error_handler (" dipole_transitionm&
           &oment_uf: write of orbital dipole moments failed")
         End Subroutine calc_write_tm
!
      End Subroutine dipole_transitionmoment_uf
  !*************************************************************
      Subroutine dipole_make_xes_spectra (iounit)
    !  Purpose: calculates the intensities for a xes spectrum
    !           the intensity is given by the square of the dipoltransitionmoment
    !           between a core orbital specified in the input and an other
    !           orbital times the Energiedifferenz to the power 4
    !           I= (Ei-Ef)^4*|<i|z|f>|^2
         Implicit None
    !------------ Declaration of formal parameters ---------------
         Integer :: iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
         Integer (Kind=i4_kind) :: n_ir, i_ir1, i_pa1, i_ip1, i_ir2, &
        & i_pa2, i_ip2, i_spin, n_spin, n_orb1, n_orb2, alloc_stat, ij, &
        & ij_2, i_xyz, i_bas1, i_bas2
         Real (Kind=r8_kind), Pointer :: eigenvector1 (:), eigenvector2 &
        & (:), mat_p1 (:), mat_p2 (:, :)
         Real (Kind=r8_kind) :: tm, eig_dif, intensity
         Real (Kind=r8_kind), Allocatable :: help_vec (:)
         Type (dipole_integral_type), Pointer :: di
    !------------ Executable code --------------------------------
!
         n_spin = symmetry_data_n_spin ()
!
         n_ir = symmetry_data_n_irreps ()
         If (core_hole_irrep > n_ir .Or. core_hole_irrep < 1) Call &
        & error_handler ('dipole_make_xes_spectra: wrong value for core_hole_irrep')
         If (core_hole_spin > n_spin .Or. core_hole_spin < 1) Call &
        & error_handler ('dipole_make_xes_spectra:  wrong value for core_hole_spin')
         If (core_hole_partner > symmetry_data_n_partners(core_hole_irrep)&
        & .Or. core_hole_spin < 1) Call error_handler&
        & ('dipole_ma e_xes_spectra:  wrong value for core_hole_partner')
         If (i_core_hole > symmetry_data_dimension(core_hole_irrep) &
        & .Or. i_core_hole < 1)&
        & Call error_handler ('dipole_make_xes_spectra:  wrong value for i_core_hole')
         i_ip1 = 1
         n_orb1 = symmetry_data_dimension (core_hole_irrep)
         Do i_ir1 = 1, core_hole_irrep - 1
            Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
               i_ip1 = i_ip1 + 1
            End Do
         End Do
         i_ip1 = i_ip1 + core_hole_partner - 1
         i_ip2 = 0
         eigenvector1 => eigvec(core_hole_irrep)%m(:, i_core_hole, &
        & core_hole_spin)
         Do i_ir2 = 1, n_ir
            n_orb2 = symmetry_data_dimension (i_ir2)
            Allocate (help_vec(n_orb2), Stat=alloc_stat)
            If (alloc_stat /= 0) Call error_handler ('dipole_transition moment_uf: allocate help_vec')
            Do i_pa2 = 1, symmetry_data_n_partners (i_ir2)
               i_ip2 = i_ip2 + 1
               Do i_xyz = 1, 3
                  If (i_ip2 .Gt. i_ip1) Then
                     di => dipole_integrals (i_xyz, i_ip1, i_ip2)
                  Else
                     di => dipole_integrals (i_xyz, i_ip2, i_ip1)
                  End If

                  If ( .Not. (di % use == DIPOLE_UNUSED .Or. x_vec(i_xyz) &
                 & == 0.0_r8_kind)) Then
                     help_vec = 0.0_r8_kind
                     If (di % use == DIPOLE_DIAGONAL) Then
                        ij = 1
                        mat_p1 => di%diagonal
                        Do i_bas2 = 1, n_orb1
                           Do i_bas1 = 1, i_bas2
                              help_vec (i_bas2) = help_vec (i_bas2) + &
                             & mat_p1 (ij) * eigenvector1 (i_bas1)
                              ij = ij + 1
                           End Do
                           ij_2 = ij + i_bas2 - 1
                           Do i_bas1 = i_bas2 + 1, n_orb1
                              help_vec (i_bas2) = help_vec (i_bas2) + &
                             & mat_p1 (ij_2) * eigenvector1 (i_bas1)
                              ij_2 = ij_2 + i_bas1
                           End Do
                        End Do
                     Else
                        mat_p2 => di%offdiagonal
                        If (i_ip2 .Gt. i_ip1) Then
                           Do i_bas2 = 1, n_orb2
                              Do i_bas1 = 1, n_orb1
                                 help_vec (i_bas2) = help_vec (i_bas2) &
                                & + mat_p2 (i_bas1, i_bas2) * &
                                & eigenvector1 (i_bas1)
                              End Do
                           End Do
                        Else
                           Do i_bas2 = 1, n_orb2
                              Do i_bas1 = 1, n_orb1
                                 help_vec (i_bas2) = help_vec (i_bas2) &
                                & + mat_p2 (i_bas2, i_bas1) * &
                                & eigenvector1 (i_bas1)
                              End Do
                           End Do
                        End If
                     End If
                     Do i_spin = 1, symmetry_data_n_spin ()
                        Do i_bas2 = 1, n_occo (i_ir2, i_spin)
                           eigenvector2 => eigvec(i_ir2)%m(:, i_bas2, &
                          & core_hole_spin)
                           eig_dif = &
                          & eigval(core_hole_irrep)%m(i_core_hole, &
                          & core_hole_spin) - eigval(i_ir2)%m(i_bas2, &
                          & i_spin)
                           tm = 0.0_r8_kind
                           Do i_bas1 = 1, n_orb2
                              tm = tm + help_vec (i_bas1) * &
                             & eigenvector2 (i_bas1)
                           End Do
                           intensity = tm * tm * eig_dif * eig_dif * &
                          & eig_dif * eig_dif * x_vec (i_xyz) * &
                          & occ_num(i_ir2)%m(i_bas2, i_spin)
                           Write (iounit, '(2F25.15)') &
                          & eigval(i_ir2)%m(i_bas2, i_spin), intensity
                        End Do ! i_bas2
                     End Do ! i_spin
                  End If
               End Do ! i_xyz
            End Do ! i_pa2
            Deallocate (help_vec, Stat=alloc_stat)
            If (alloc_stat /= 0) Call error_handler &
           & ('dipole_transition moment_uf: deallocate help_vec')
         End Do ! i_ir2
      End Subroutine dipole_make_xes_spectra
  !*************************************************************
      Subroutine dipole_trans_response(iounit)
        !  Purpose: calculates and prints all transition dipol moments
        !  in a unformatted way to iounit.
        !
        ! NOTE: At the moment this subroutine is almost identical
        !       to the subroutine "dipole_transitionmoment_uf" above !
        !       Only difference: no MO eigenvalues are considered
        !
        USE comm_module
        Implicit None
        !------------ Declaration of formal parameters ---------------
        Integer (Kind=i4_kind), Intent (In) :: iounit
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        Integer (Kind=i4_kind) :: n_ir, i_ir1, i_pa1, i_ip1, i_ir2, &
             & i_pa2, n_pa2, i_ip2, i_spin, n_spin, i_orb1, n_orb1, i_orb2, &
             & n_orb2, iostat !!$ m
        Logical :: open_shell
        Real (Kind=r8_kind), Pointer :: eigenvector1 (:), eigenvector2 &
             & (:)
        !------------ Executable code --------------------------------
        !
        n_spin = symmetry_data_n_spin ()
        open_shell = n_spin .Eq. 2
        !
        n_ir = symmetry_data_n_irreps ()
        !
        i_ip1 = 0
        irrep1: Do i_ir1 = 1, n_ir
           n_orb1 = symmetry_data_dimension (i_ir1)
           partner1: Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
              i_ip1 = i_ip1 + 1
              !
              i_ip2 = 0
              irrep2: Do i_ir2 = 1, i_ir1
                 n_orb2 = symmetry_data_dimension (i_ir2)
                 If (i_ir1 .Eq. i_ir2) Then
                    n_pa2 = i_pa1
                 Else
                    n_pa2 = symmetry_data_n_partners (i_ir2)
                 End If
                 partner2: Do i_pa2 = 1, n_pa2
                    i_ip2 = i_ip2 + 1
                    !
                    spin: Do i_spin = 1, n_spin
                       !
                       orbitals1: Do i_orb1 = 1, n_orb1
                          eigenvector1 => eigvec(i_ir1)%m(:, i_orb1, &
                               & i_spin)
                          If (i_ip1 .Eq. i_ip2) n_orb2 = i_orb1
                          orbitals2: Do i_orb2 = 1, n_orb2
                             eigenvector2 => eigvec(i_ir2)%m(:, &
                                  & i_orb2, i_spin)
                             !
                             Call calc_write_tm ()
                             !
                          End Do orbitals2
                       End Do orbitals1
                       !
                    End Do spin
                    !
                 End Do partner2
              End Do irrep2
           End Do partner1
        End Do irrep1

        !! FIXME: SB:Who need this staff??
        ! at the end also write total dipole moment
!!$        write(iounit,iostat=iostat) (dipole_total(m), m=1,3)
!!$        if ( iostat .ne. 0 ) call error_handler( &
!!$             " dipole_trans_response: write tot. dip. moment failed !")
        !
      Contains
        !
        Subroutine calc_write_tm ()
          use debug
          ! Porpose: calculation and output of one transition moment between two orbitals
          Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12, i_xyz
          Type (dipole_integral_type), Pointer :: di
          Real (Kind=r8_kind) :: transitionmoments (3)
          ASSERT(i_ip2 <= i_ip1)
          Do i_xyz = 1, 3
             If (i_ip2 .Gt. i_ip1) Then
                di => dipole_integrals (i_xyz, i_ip1, i_ip2)
             Else
                di => dipole_integrals (i_xyz, i_ip2, i_ip1)
             End If
             integralstorage: Select Case (di % use)
             Case (DIPOLE_DIAGONAL) integralstorage
                !triangular storage
                i_bas12 = 1
                transitionmoments (i_xyz) = 0.0_r8_kind
                bas1_diagonal: Do i_bas1 = 1, n_orb1
                   bas2_diagonal: Do i_bas2 = 1, i_bas1 - 1
                      transitionmoments (i_xyz) = transitionmoments &
                           & (i_xyz) + &
                           & (eigenvector1(i_bas1)*eigenvector2(i_bas2)+&
                           & eigenvector1(i_bas2)*eigenvector2(i_bas1)) * &
                           & di%diagonal(i_bas12)
                      i_bas12 = i_bas12 + 1
                   End Do bas2_diagonal
                   transitionmoments (i_xyz) = transitionmoments &
                        & (i_xyz) + eigenvector1 (i_bas1) * eigenvector2 &
                        & (i_bas1) * di%diagonal(i_bas12)
                   i_bas12 = i_bas12 + 1
                End Do bas1_diagonal
             Case (DIPOLE_OFFDIAGONAL) integralstorage
                ! full storage
                transitionmoments (i_xyz) = 0.0_r8_kind
                If (i_ip2 .Gt. i_ip1) Then
                   Do i_bas2 = 1, n_orb2
                      Do i_bas1 = 1, n_orb1
                         transitionmoments (i_xyz) = &
                              & transitionmoments (i_xyz) + eigenvector1 &
                              & (i_bas1) * eigenvector2 (i_bas2) * &
                              & di%offdiagonal(i_bas1, i_bas2)
                      End Do
                   End Do
                Else
                   Do i_bas1 = 1, n_orb1
                      Do i_bas2 = 1, n_orb2
                         transitionmoments (i_xyz) = &
                              & transitionmoments (i_xyz) + eigenvector1 &
                              & (i_bas1) * eigenvector2 (i_bas2) * &
                              & di%offdiagonal(i_bas2, i_bas1)
                      End Do
                   End Do
                End If
             Case (DIPOLE_UNUSED) integralstorage
                transitionmoments (i_xyz) = 0.0_r8_kind
             Case Default integralstorage
                Call error_handler ("dipole_transitionmoment_f: forbi&
                     &dden case")
             End Select integralstorage
          End Do
          !
          !!SB: i_spin, i_ir1, i_pa1 added
          Write (iounit, IoStat=IoStat) i_spin, i_ir1, i_pa1, i_orb1, i_ir2, i_pa2, i_orb2, &
               & (transitionmoments(i_xyz), i_xyz=1, 3)
          If (iostat .Ne. 0) Call error_handler&
               &(" dipole_transitionmoment_uf: write of orbital dipole moments failed")
        End Subroutine calc_write_tm
        !
      End Subroutine dipole_trans_response

      Subroutine dipole_trans_response_v2(iounit)
        !  Purpose: calculates and prints all transition dipol moments
        !  in a unformatted way to iounit.
        !
        ! NOTE: At the moment this subroutine is almost identical
        !       to the subroutine "dipole_transitionmoment_uf" above !
        !       Only difference: no MO eigenvalues are considered
        !
        USE comm_module
        Implicit None
        !------------ Declaration of formal parameters ---------------
        Integer (Kind=i4_kind), Intent (In) :: iounit
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        Integer(Kind=i4_kind) :: n_ir, i_ir1, i_pa1, i_orb1, n_orb1,&
             &                          i_ir2, i_pa2, i_orb2, n_orb2
        integer(kind=i4_kind) :: i_spin, n_spin, iostat
        Real (Kind=r8_kind), Pointer :: eigenvector1(:), eigenvector2(:)
        !------------ Executable code --------------------------------
        !
        n_spin = symmetry_data_n_spin ()
        n_ir   = symmetry_data_n_irreps ()

        spin: Do i_spin = 1, n_spin
           irrep1: Do i_ir1 = 1, n_ir
              n_orb1 = symmetry_data_dimension (i_ir1)
              irrep2: Do i_ir2 = 1, n_ir

                 print *," i_ir1 = ",i_ir1," i_ir2 = ",i_ir2

                 partner1: Do i_pa1 = 1, symmetry_data_n_partners (i_ir1)
                    partner2: Do i_pa2 = 1, symmetry_data_n_partners (i_ir2)
                       print *," i_pa1 = ",i_pa1," i_pa2 = ",i_pa2

                       orbitals1: Do i_orb1 = 1, n_orb1
                          eigenvector1 => eigvec(i_ir1)%m(:,i_orb1,i_spin)
                          orbitals2: Do i_orb2 = 1, n_orb2

                             eigenvector2 => eigvec(i_ir2)%m(:,i_orb2,i_spin)
                             call calc_write_tm (iounit,i_spin,i_ir1,i_pa1,i_orb1,&
                                  &                            i_ir2,i_pa2,i_orb2,&
                                  &                            eigenvector1, eigenvector2)
                          End Do orbitals2
                       End Do orbitals1
                    End Do partner2
                 End Do partner1
              End Do irrep2
           End Do irrep1
        End Do spin

      Contains

        Subroutine calc_write_tm (iounit,i_spin,i_ir1,i_pa1,i_orb1,&
             &                                  i_ir2,i_pa2,i_orb2,&
!!$             &                                  i_ip1, i_ip2, &
             &                                  eigenvector1, eigenvector2)
          USE debug
          integer(i4_kind), intent(IN) :: iounit,i_spin,i_ir1,i_pa1,i_orb1,i_ir2,i_pa2,i_orb2
!!$          integer(i4_kind), intent(IN) :: i_ip1, i_ip2
          real(r8_kind),    intent(IN) :: eigenvector1(:), eigenvector2(:)
          ! Porpose: calculation and output of one transition moment between two orbitals
          Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12, i_xyz
          Integer (Kind=i4_kind) :: i_ip1, i_ip2
          Type (dipole_integral_type), Pointer :: di
          Real (Kind=r8_kind) :: transitionmoments (3)
          Do i_xyz = 1, 3
             i_ip1 = symmetry_data_i_ip(i_ir1,i_pa1)
             i_ip2 = symmetry_data_i_ip(i_ir2,i_pa2)
!!$             di => dipole_integrals (i_xyz, i_ip1, i_ip2)
             print *,"i_ip1 = ",i_ip1," i_ip2 = ",i_ip2
             If (i_ip2 .Gt. i_ip1) Then
                di => dipole_integrals (i_xyz, i_ip1, i_ip2)
             Else
                di => dipole_integrals (i_xyz, i_ip2, i_ip1)
             End If

             integralstorage: Select Case (di % use)
             Case (DIPOLE_DIAGONAL) integralstorage
                !triangular storage
                i_bas12 = 1
                transitionmoments (i_xyz) = 0.0_r8_kind
                bas1_diagonal: Do i_bas1 = 1, n_orb1
                   bas2_diagonal: Do i_bas2 = 1, i_bas1 - 1
                      transitionmoments(i_xyz) = transitionmoments(i_xyz)  &
                           & + (eigenvector1(i_bas1)*eigenvector2(i_bas2)  &
                           & +  eigenvector1(i_bas2)*eigenvector2(i_bas1)) &
                           & *  di%diagonal (i_bas12)
                      i_bas12 = i_bas12 + 1
                   End Do bas2_diagonal
                   transitionmoments(i_xyz) = transitionmoments(i_xyz) &
                        & + eigenvector1(i_bas1) &
                        & * eigenvector2(i_bas1) &
                        & * di%diagonal (i_bas12)
                   i_bas12 = i_bas12 + 1
                End Do bas1_diagonal
             Case (DIPOLE_OFFDIAGONAL) integralstorage
                ! full storage
                print *,"shape(di%offdiagonal) = ",shape(di%offdiagonal)
                print *," n_orb1 = ",n_orb1," n_orb2 = ",n_orb2
                If (i_ip2 .Gt. i_ip1) Then  !!WHY THIS???
                   Do i_bas2 = 1, n_orb2
                      Do i_bas1 = 1, n_orb1
                         transitionmoments(i_xyz) = transitionmoments(i_xyz) &
                              & + eigenvector1  (i_bas1) &
                              & * eigenvector2  (i_bas2) &
                              & * di%offdiagonal(i_bas1,i_bas2)
                      End Do
                   End Do
                Else
                   Do i_bas1 = 1, n_orb1
                      Do i_bas2 = 1, n_orb2
                         transitionmoments(i_xyz) = transitionmoments(i_xyz) &
                              & + eigenvector1  (i_bas1) &
                              & * eigenvector2  (i_bas2) &
                              & * di%offdiagonal(i_bas2, i_bas1)
                      End Do
                   End Do
                End If
             Case (DIPOLE_UNUSED) integralstorage
                transitionmoments (i_xyz) = 0.0_r8_kind
             Case Default integralstorage
                Call error_handler ("dipole_trans_response: forbidden case")
             End Select integralstorage
          End Do
          !
          !!SB: i_spin, i_ir1, i_pa1 added
!!$          print *,"i_spin = ",i_spin," i_ir1 = ",i_ir1," i_pa1 = ",i_pa1," i_orb1 = ",i_orb1
!!$          print *,                   " i_ir2 = ",i_ir2," i_pa2 = ",i_pa2," i_orb2 = ",i_orb2
!!$          print *,"transitionmoments(1:3) = ",transitionmoments(:)
          Write (iounit, iostat=iostat) i_spin, i_ir1, i_pa1, i_orb1, i_ir2, i_pa2, i_orb2, &
               & (transitionmoments(i_xyz), i_xyz=1, 3)
          ASSERT(iostat==0)
        End Subroutine calc_write_tm
        !
      End Subroutine dipole_trans_response_v2


  !*************************************************************
End Module dipole_module

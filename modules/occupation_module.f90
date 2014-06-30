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
module occupation_module
  !
  !  Purpose: Contains  routines and information connected
  !           with the occupation of levels
  !  Contents:
  !    Variables (PUBLIC)   occ_num     occupation number
  !                                     for level (i,m,is)
  !                                     where i = IRREP,
  !                                     m = Orbital within i and
  !                                     is = spin index
  !                                     type(arrmat2) (:),i.e.
  !                                     occ_num(i)%m(m,is)
  !
  !                         n_occo      number of occupied levels
  !                                     par IRREP and spin
  !
  !                        occ_num_core occupation number
  !                                     for core level (i,m,is)
  !
  !                         n_occo_core number of occupied core levels
  !                                     par IRREP
  !
  !
  !   Routines
  !     (i)  PUBLIC         occupation_read read charge
  !                                     from input file
  !
  !                         occupation_write write in input format
  !
  !                         get_n_elec  access to private variable n_elec
  !                                     (total number of electrons)
  !
  !                     get_n_core_elec access to private variable
  !                                     n_core_elec
  !
  !                      get_spin_diff  access to priv. var. magn_moment
  !                                     (total magnetic moment of the elec.)
  !
  !                get_fix_spin_switch  access to priv. var. fixed_spin_diff
  !
  !                  get_min_spin_diff  access to priv. par. min_spin_diff
  !                                     (threshold for spin dens. normaliz.)
  !
  !                get_diagonal_offset  access to priv. par. diagonal_offset
  !                                     (amount of diagonal presetting)
  !
  !                         reoccup     Do the re-occupation of levels
  !                                     after the diagonalization of the
  !                                     Hamiltonian
  !
  !                          read_ocup  Read the variable OCUP written
  !                                     by the old lcgto from file and
  !                                     pack it into occ_num
  !
  !                      print_occ_num  Print out the variable occ_num
  !                                     in a pretty format well-suited for
  !                                     comparison to the old lcgto
  !
  !              occupation_level_sort  Sort the energy levels eigval
  !                                     in ascending order, producing
  !                                     a list which contains the index
  !                                     number of the corresponding
  !                                     energy level
  !
  !                      alloc_occ_num  allocate occ_num (and occ_num_core if required)
  !
  !                      alloc_n_occo   allocate n_occo (and n_occo_core if required)
  !
  !                occupation_spindiff   calculate the spin difference
  !
  !          occupation_print_spectrum  print out the occupation spectrum in nearly
  !                                     the same format as in the old LCGTO
  !
  !       occupation_print_popspectrum  print spectrum including a user defined
  !                                     population analysis
  !
  !   (ii)  PRIVATE  alloc_occ_num_occ  allocate occ_num_core
  !
  !
  !  Module called by: density_data_module->pre_dens_master
  !                    scf, chargefit,build_hamiltonian
  !
  !
  ! Fixed Hole Occupation:
  ! concept of variables: input variables are hole_list(:),hole_irrep(:) and n_holes.
  ! hole_irrep(i_hole)   Index of Irrep if hole i_hole
  ! hole_list(i_hole)    Index of Orbital within hole_irrep(i_hole) corresponding to i_hole
  ! hole_spin(i_hole)    Spin of i_hole
  ! These variables are mapped to the somewhat more convenient variables
  ! n_holes_per_irrep(:) and holes_per_irrep(:) during the routine
  ! 'eigen_hole_setup' which is called indirectly by 'do_recover' at the
  ! beginning of the SCF-cycle.
  ! These variables are stored  in the 'eigen_data_module'.
  ! Since it is not clear that the symmetry-part also being run (it is well possible
  ! that somebody only tried to read the input) we have to test if these arrays
  ! have been allocated in 'occupation_symmetry_check'.
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification: added occ_num_occ,send_occ_num_occ,
  !              recv_occ_num_occ,alloc_occ_num_occ,
  !              seek_occupied_levels
  ! Author: F.Noertemann
  ! Date:   12/95
  ! Description: see above documentation
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modification
  ! Author: Uwe Birkenheuer
  ! Date:   6/97
  ! Description: Introduction of private variables "magn_moment"
  !              to control the initialization (and progress) of
  !              a spin-polarized MDA calculation
  !              Introduction of private variables "diagonal_offset"
  !              to control the initialization of the xc matrix
  !
  ! Modification
  ! Author: Markus Staufer
  ! Date:   8/98
  ! Description: Additional subroutine occupation_2d_correct
  !              was added. See the subroutine header for details
  !
  ! Modification
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !
  ! Modification
  ! Author: Thomas Seemueller
  ! Date:   12/99
  ! Description: Subroutine reoccup was extended for the case of
  !              user defined number of unpaired electrons
  !              Introduction of private variable "fixed_spin_diff"
  !              Additional subroutine get_fix_spin_switch
  !              Additional subroutine check_occ_fermi (see
  !              subroutine header for detailed information)
  !
  ! Modification
  ! Author: AM
  ! Date: 01/2002
  ! Description: subroutine reoccup extended for another
  !              fixed occupation mode which fixes all occupation
  !              numbers, within irrep too,
  !              e.g.: 1a1^2, 2a1^0, 3a1^2
  !                    1a2^0
  !                    1b1^2, 2b1^1 ...
  !              (search for FORCE_SPECTRUM switch (lcs))
  !
  ! Modification (Please copy before crediting)
  ! Author:
  ! Date:
  ! Description:
  !
  !----------------------------------------------------------------
#include "def.h"
  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use datatype    ! user defined datatypes
  use symmetry_data_module  ! symmetry information
  use init_module    ! init
  use iounitadmin_module
  use operations_module, only: operations_echo_input_level, &
                               operations_core_density
#ifdef FPP_DEBUG
  use error_module, only: MyID
#endif
  implicit none
  private
  save
  !== Interrupt end of public interface of module =================

  !------------ Declaration of public constants and variables -----
  type, private :: irrep_t
     real(r8_kind), allocatable :: occ(:)
  end type irrep_t
  type(irrep_t), allocatable :: occupation_numbers(:,:) ! (n_spin,n_irr)

  public :: arrmat3,arrmat2,arrmat2int,sym
  integer(kind=i4_kind),allocatable,public :: n_occo(:,:)
  ! n_occo(i_spin,i_irrep)  number of occupied levels
  integer(kind=i4_kind),allocatable,public :: n_occo_core(:,:)
  ! n_occo_core(i_spin,i_irrep)  number of occupied core levels
  type(arrmat2), allocatable, target, public :: occ_num(:)
  ! occ_num(i_irrep)%m(i_orbital,i_spin) occupation number
  type(arrmat2),allocatable,public         :: occ_num_core(:)
  ! occ_num_core(i_irrep)%m(i_orbital,i_spin) core level occupation number
  real(kind=r8_kind),allocatable,public    :: noc_start(:,:),noc_fixed(:,:)
  ! noc_start(i_spin,n_irrep) number of occupied levels per irreps and spin
  ! for simple fixed occupation.
  ! noc_fixed(i_spin,i_irrep) is needed as help variable durig reoccupation
  ! with fixed occupation, since noc_start has to be left unchanged.
  ! variables needed for HOLES -------------------------------------------
  integer(kind=i4_kind),allocatable,public :: hole_list(:),hole_irrep(:),&
       hole_spin(:)
  real(kind=r8_kind),allocatable,public    :: hole_occnum(:)
  type(arrmat1),allocatable,public             :: occnum_hole(:)
  ! hole_list(i_hole)      index of MO with hole within irrep
  ! hole_irrep(i_hole)     irrep of hole
  ! hole_spin(i_hole)      spin of hole
  ! hole_occnum(i_hole)    remaining occupation number of hole
  ! occnum_hole(i_irrep)%m(i_hole)  occupation number of hole <i_hole> in
  !                                 Irrep i_irrep.
  ! ----------------------------------------------------------------------
  type(arrmat2int),allocatable,public      :: num_list(:)
  integer(kind=i4_kind),allocatable,public :: n_rot(:,:)
  type(arrmat3),allocatable,public         :: eigvec_kept(:)
  logical,public                           :: eigen_kept = .false.
  logical,public                           :: keep_eigen = .false.
  logical,public                           :: fixed,fixed_hole
  logical                                  :: hole_localization
  logical                                  :: hole_update
  logical                                  :: fixed_spin_diff
  !------------ public functions and subroutines ------------------
  public :: occupation_read,  occupation_write, &
       get_n_elec, get_n_core_elec, reoccup, read_ocup, print_occ_num, &
       occupation_level_sort, occupation_2d_correct, &
       alloc_occ_num, alloc_n_occo, occupation_spindiff, &
       occupation_print_spectrum,occupation_shutdown, &
       get_spin_diff, get_df_spin_diff, put_magn_moment,&
       get_fix_spin_switch, get_df_fix_spsw,&
       put_fix_spin_switch, get_min_spin_diff, &
       get_diagonal_offset, eigenstates_store, eigenstates_recover, &
       occupation_print_popspectrum, occupation_symmetry_check,&
       occupation_fixed,occupation_get_holes, occupation_fixed_hole,&
       occupation_jz, count_electrons, get_charge

  !
  ! FIXME: target is here only to make things like
  !
  !     v => eigvec_occ_real(irrep)%m
  !     n => occ_num_occ(irrep)%m
  !
  ! compile also when arrmat2 is declared with
  ! allocatable component instead of pointer.
  !

  !================================================================
  ! End of public interface of module
  !================================================================
  real(kind=r8_kind),parameter :: convert1 = 27.211652_r8_kind, &
                                  convert2 = 1.0_r8_kind / convert1

  !------------ Default values for input parameters ---------------
  real(kind=r8_kind) :: df_charge          = 0.0_r8_kind   , &
                        df_magn_moment     = 0.0_r8_kind   , &
                        df_min_spin_diff   = 1.0E-3_r8_kind, &
                        df_diag_offset_ev  = 1.0_r8_kind   , &
                        df_diagonal_offset = convert2  ! kept for compatibility
  integer(kind=i4_kind) :: df_n_nonempty_irreps = 0_i4_kind, &
                           df_n_holes           = 0_i4_kind
  logical ::&
       & df_fixed = .false.,&
       & df_hole_localization=.false.,&
       & df_hole_update=.false.,&
       & df_fixed_spin_diff = .false.,&
       & df_force_spectrum = .false.

  !------------ occupation input parameters -----------------------
  ! none of these variables is currently needed on the slaves, thus
  ! there is no occupation_pack and occupation_unpack routine so far (UB)
  real(kind=r8_kind) :: charge         , &
                        magn_moment    , &
                        min_spin_diff  , &
                        diag_offset_ev , &
                        diagonal_offset  ! kept for compatibility
  integer(kind=i4_kind) :: n_nonempty_irreps,n_holes

  !------------ Declaration of private constants and variables ----
  real(kind=r8_kind) :: n_elec = 0.0_r8_kind, n_core_elec = 0.0_r8_kind
  logical            :: check_num_list=.false.

  logical            :: force_spectrum

  namelist /occupation/ charge, fixed, magn_moment, min_spin_diff, &
                        diag_offset_ev, diagonal_offset, n_nonempty_irreps, &
                        hole_localization, hole_update, n_holes, &
                        fixed_spin_diff, force_spectrum


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine alloc_n_occo(ssym)
    !  Purpose: allocate the appropriate space for n_occo
    !** End of interface *****************************************
    use options_module, only: options_spin_orbit
    !------------ Declaration of local variables -----------------
    type(sym),intent(in)       :: ssym
    integer(kind=i4_kind)       :: alloc_stat
    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind)                :: n_irrep
    ! n_irrep    : number of irreps
    !------------ Executable code --------------------------------

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       n_irrep = ssym%n_proj_irrep
    else
       n_irrep = ssym%n_irrep
    endif

    allocate(n_occo(ssym%n_spin,n_irrep),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler&
         ("ALLOC_N_OCCO : allocation failed")
    if(operations_core_density)then
       allocate(n_occo_core(ssym%n_spin,n_irrep),STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler&
            ("ALLOC_N_OCCO : allocation (core) failed")
    endif

  end subroutine alloc_n_occo
  !*************************************************************

  subroutine alloc_occ_num()
    !
    ! Allocate the appropriate space for occ_num and initialize with
    ! zero
    !
    use options_module, only: options_spin_orbit
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: i_gamma, alloc_stat, n

    ! Dimensions of irreps (in order to account for SPIN ORBIT more
    ! easily)
    integer (i4_kind) :: n_irrep
    integer (i4_kind), allocatable :: dim_irrep(:)
    ! n_irrep: number of irreps
    ! dim_irrep: number of independent functions in irrep

    DPRINT MyID//'alloc_occ_num'

    ! Set appropriate  dimensions of irreps (use  projective irreps in
    ! case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym % n_proj_irrep
       allocate (dim_irrep(n_irrep))
       do n = 1, n_irrep
          dim_irrep(n) = ssym % dim_proj(n)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym % n_irrep
       allocate (dim_irrep(n_irrep))
       do n = 1, n_irrep
          dim_irrep(n) = ssym % dim(n)
       enddo
    endif ! options_spin_orbit

    allocate (occ_num(n_irrep), STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    do i_gamma = 1, n_irrep
       allocate (occ_num(i_gamma) % m(dim_irrep(i_gamma), ssym%n_spin), STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    enddo
    call init (occ_num)

    if(operations_core_density)then
       allocate (occ_num_core(n_irrep),STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler&
            ("ALLOC_OCC_NUM : allocation (core) failed")
       do i_gamma=1,n_irrep
          allocate( occ_num_core(i_gamma)%m(dim_irrep(i_gamma),&
               ssym%n_spin),STAT=alloc_stat)
          if(alloc_stat.ne.0) call error_handler&
               ("ALLOC_OCC_NUM : allocation (core) failed")
       enddo
       call init(occ_num_core)
    endif

    ! deallocate appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    deallocate(dim_irrep)

  end subroutine alloc_occ_num

  subroutine occupation_shutdown()
    !
    ! Deallocate the  variables n_occo, n_occo_core, occ_num, occupied,
    ! num_list at the end of the SCF-Cycles
    !
    ! Subroutine called by : main_scf()
    !
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: alloc_stat, i, j

    check_num_list = .false.

    if (allocated (n_occo)) then
       deallocate(n_occo, STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if (allocated (occ_num)) then
       do i = 1, size (occ_num)
          deallocate (occ_num(i) % m, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
       deallocate(occ_num, STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if (allocated (num_list)) then
        do i = 1, size(num_list)
            deallocate(num_list(i)%m, STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        enddo
        deallocate(num_list, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (n_occo_core)) then
        deallocate(n_occo_core, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (occ_num_core)) then
        do i = 1, size(occ_num_core)
            deallocate(occ_num_core(i)%m, STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        enddo
        deallocate(occ_num_core, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    endif

    if (allocated (noc_start)) then
        deallocate(noc_start, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (noc_fixed)) then
        deallocate(noc_fixed, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (occupation_numbers)) then
        do i = 1, size(occupation_numbers,1)
            do j = 1, size(occupation_numbers,2)
                deallocate(occupation_numbers(i, j)%occ, stat=alloc_stat)
                ASSERT(alloc_stat==0)
            enddo
        enddo
        deallocate(occupation_numbers, stat=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (hole_irrep)) then
        deallocate(hole_irrep, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (hole_list)) then
        deallocate(hole_list, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (hole_spin)) then
        deallocate(hole_spin, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    if (allocated (hole_occnum)) then
        deallocate(hole_occnum, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    ! this variable will only be allocated if the symmetry_part is run
    if (allocated (occnum_hole)) then
        do i = 1, size(occnum_hole)
            deallocate(occnum_hole(i)%m, STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        enddo
        deallocate(occnum_hole, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    end if

    DPRINT MyID//'occupation_shutdown done'
  end subroutine occupation_shutdown

  !*************************************************************
  subroutine occupation_read()
    !
    ! Read in  occupation control data from  input file.  informations
    ! in options_module are used and must exist before
    !
    use unique_atom_module
    use options_module, only: options_xcmode, &
                              xcmode_model_density, xcmode_extended_mda, &
                              options_n_spin, options_spin_restricted,&
                              options_recover,recover_eigenvec,&
                              & options_spin_orbit
    use input_module
    use iounitadmin_module, only: output_unit
    implicit none
    !** End of interface *****************************************

    integer :: i_ua, unit, status,i,is,alloc_stat,n_spin
    integer(i4_kind) :: n
    real(r8_kind)    :: real_buf(999), real_number
    ! increase for more occupation numbers per line


    n_elec      = 0.0_r8_kind
    n_core_elec = 0.0_r8_kind
    do i_ua = 1,N_unique_atoms
       n_elec      = n_elec      + unique_atoms(i_ua)%N_equal_atoms &
                                 * unique_atoms(i_ua)%Z
       n_core_elec = n_core_elec + unique_atoms(i_ua)%N_equal_atoms &
                                 * unique_atoms(i_ua)%Zc
    enddo

    if (output_unit > 0) then
       write (output_unit, *) ' number of electrons ', n_elec
    endif
    !------ Re-Definition of input dependent default valuse ----
    if ( ( options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda ) .and. &
         options_n_spin() > 1 ) then
       df_magn_moment = 1.0_r8_kind
       df_diag_offset_ev = 0.0_r8_kind
       df_diagonal_offset = 0.0_r8_kind ! kept for compatibility
    endif

    charge          = df_charge
    magn_moment     = df_magn_moment
    min_spin_diff   = df_min_spin_diff
    diag_offset_ev  = df_diag_offset_ev
    diagonal_offset = df_diagonal_offset ! kept for compatibility

    fixed           = df_fixed
    force_spectrum  = df_force_spectrum
    n_nonempty_irreps = df_n_nonempty_irreps

    fixed_spin_diff = df_fixed_spin_diff   !new
    hole_localization = df_hole_localization
    hole_update     = df_hole_update
    n_holes         = df_n_holes

    if ( input_line_is_namelist("occupation") ) then
       call input_read_to_intermediate
       unit = input_intermediate_unit()
       read(unit, nml=occupation, iostat=status)
       if (status .gt. 0) call input_error( &
            "occupation_read: namelist occupation.")
    endif

    fixed_spin_diff= (.not.options_spin_restricted()).and.fixed_spin_diff

    ! to ensure compatibility
    if (diagonal_offset /= df_diagonal_offset .and. &
         diag_offset_ev == df_diag_offset_ev) then
       diag_offset_ev = diagonal_offset * convert1
    endif
    ! charge is the electrostatic charge while n_elec counts electrons, thus
    n_elec = n_elec - charge
    if(n_core_elec > n_elec.and.operations_core_density)call input_error( &
         "occupation_read: number of core electrons exceeds number of valence electrons")
    if (n_core_elec > 0.0_r8_kind .and. .not.operations_core_density )then
       ! in case of a pseudopot calc n_elec holds the number of valence electrons
       n_elec = n_elec - n_core_elec
    endif
    if(operations_core_density .and. n_core_elec == 0.0_r8_kind) call input_error( &
         "occupation_read: Zc must be positive to create a core density.")

    if (options_spin_restricted() ) then
       n_spin=1
    else
       n_spin=2
    endif

    if (fixed) then
       ! check if number of non-empty irreps was specified
       if (n_nonempty_irreps ==df_n_nonempty_irreps) then
          call input_error("occupation_read: you forgot to specify the number of non-empty irreps")
       endif

       ! allocate the noc_start array
       allocate(noc_start(n_spin,n_nonempty_irreps),STAT=alloc_stat)
       if (alloc_stat/=0) call input_error &
            ("occupation_read: allocation of noc_start failed")
       allocate(noc_fixed(n_spin,n_nonempty_irreps),STAT=alloc_stat)
       if (alloc_stat/=0) call input_error &
            ("occupation_read: allocation of noc_fixed failed")

       noc_start=99.0_r8_kind
       noc_fixed=99.0_r8_kind
       ! read in the occupation numbers for IRREPS
       do i=1,n_nonempty_irreps
          call input_read_to_intermediate
          unit = input_intermediate_unit()
          read(unit, fmt=*, iostat=status) (noc_start(is,i),is=1,n_spin)
          if (status>0) call input_error("occupation_read: noc_start")
       enddo
       if (abs(sum(noc_start) - n_elec) > 1.0E-10_r8_kind) call input_error&
            ("occupation_read: your fixed occupation is inconsistent with the number of electrons")
       if (maxval(noc_start)==99.0_r8_kind) call input_error&
            ("occupation_read: check your list if occupation numbers.")
       noc_fixed = noc_start
    endif!if (fixed)

    if(force_spectrum)then
       ! check if number of non-empty irreps was specified
       if (n_nonempty_irreps ==df_n_nonempty_irreps) then
          call input_error("occupation_read: you forgot to specify the number of non-empty irreps")
       endif

       allocate(occupation_numbers(n_spin,n_nonempty_irreps),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       real_number = 0.0_r8_kind
       do i=1,n_nonempty_irreps
       do is=1,n_spin
          call input_read_to_intermediate() ! why that complicated????
          unit = input_intermediate_unit()
          call read_line(unit,n,real_buf)
          ASSERT(n.ge.0)
          allocate(occupation_numbers(is,i)%occ(n),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          occupation_numbers(is,i)%occ(:) = real_buf(:n)
          DPRINT 'irr=',i,'is=',is,'occ=',real_buf(:n)
          real_number = real_number + sum(occupation_numbers(is,i)%occ(:))
       enddo
       enddo
       if(abs(n_elec - real_number)>1.0E-7)then
          ABORT("provided spectrum conflicts with n_elec")
       endif
    endif!if(force_spectrum)

    ! set the variable 'fixed_hole'
    if (n_holes<0_i4_kind) then
       call input_error("occupation_read: n_holes should be positive")
    elseif (n_holes>0_i4_kind) then
       fixed_hole=.true.
    endif

    if (fixed_hole) then
       if(hole_localization) then
          if (options_recover() /= recover_eigenvec &
               & .and. .not. options_spin_orbit) call input_error&
               ("occupation_read: fixed holes ONLY with read_eigenvec=.TRUE.")
       else
          if (hole_update) then
             call write_to_output_units&
                  (" occupation_read: hole_update without hole_localization&
                  &does not make sense. Ignoring this option")
             hole_update=.false.
          endif
       endif
       if (fixed) call input_error&
            ("occupation_read: fixed occupation with fixed hole not possible yet")

       allocate(hole_list(n_holes),STAT=alloc_stat)
       if (alloc_stat/=0) call input_error&
            ("occupation_read: allocation of hole_list failed")
       allocate(hole_irrep(n_holes),STAT=alloc_stat)
       if (alloc_stat/=0) call input_error&
            ("occupation_read: allocation of hole_irrep failed")
       if (.not.options_spin_restricted()) then
          allocate(hole_spin(n_holes),STAT=alloc_stat)
          if (alloc_stat/=0) call input_error&
               ("occupation_read: allocation of hole_spin failed")
          hole_spin=0_i4_kind
       endif
       allocate(hole_occnum(n_holes),STAT=alloc_stat)
       if (alloc_stat/=0) call input_error&
            ("occupation_read: allocation of hole_occnum failed")

       hole_list=-99_i4_kind
       hole_irrep=-99_i4_kind
       hole_occnum=-99.0_r8_kind
       do i=1,n_holes
          call input_read_to_intermediate
          unit = input_intermediate_unit()
          if (options_spin_restricted()) then
             read(unit, fmt=*, iostat=status) hole_list(i),hole_irrep(i),hole_occnum(i)
             if (status>0) call input_error("occupation_read: hole_list")
          else
             read(unit,fmt=*, iostat=status) hole_list(i),hole_irrep(i),hole_spin(i),hole_occnum(i)
             if (status>0) call input_error("occupation_read : hole_list")
          endif
       enddo
       if(minval(hole_list)<=0_i4_kind) call input_error&
            ("occupation_read: check your hole_list")
       if(minval(hole_irrep)<=0_i4_kind) call input_error&
            ("occupation_read: check your hole_irrep")
       if (.not.options_spin_restricted()) then
          if(minval(hole_spin)<=0_i4_kind) call input_error&
               ("occupation_read: check your hole_spin list")
       endif
       if(minval(hole_occnum)<0.0_r8_kind) call input_error&
            ("occupation_read: check occupation number for holes")
    endif

  contains
    subroutine read_line(unit,n,arr)
      use io, only: readline, IO_EOR
      implicit none
      integer(i4_kind), intent(in)    :: unit
      integer(i4_kind), intent(out)   :: n
      real(r8_kind),    intent(inout) :: arr(:)
      ! *** end of interface ***

      real(r8_kind)    :: tmp(size(arr))
      integer(i4_kind) :: siz,m,i,status

      integer(i4_kind)   :: strlen,stat
      character(len=256) :: line

      strlen = LEN(line)
      call readline(unit,line,strlen,stat)
      DPRINT 'read_line: strlen=',strlen,' stat=',stat
      DPRINT '>'//line(:strlen)//'<'
      if ( stat /= IO_EOR ) then
         WARN('input line is too long!')
         ! GET RID OF input_intermediate!
         ! or increase 256:
      endif

      n = 0
      siz = size(arr)
      try: do m=1,siz
         ! unit is a one-line internal file:
         read(line(:strlen),fmt=*, iostat=status) (tmp(i),i=1,m)
         if(status.eq.0)then
            arr(:m) = tmp(:m)
            !REWIND(unit)
         else
            n = m - 1
            exit
         endif
      enddo try
      DPRINT 'om::read_line: read ',n,' numbers:',arr(:n)
    end subroutine read_line
  end subroutine occupation_read
  !*************************************************************

  !*************************************************************
  subroutine occupation_write(iounit)
    !
    ! Write namelist occupation to iounit in input format.
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop,&
         real_format1, real_format2, echo_level_short
    use options_module, only: options_spin_restricted, options_spin_orbit
    implicit none
    integer, intent(in) :: iounit
    !** End of interface *****************************************

    integer(kind=i4_kind)                 :: i,n_spin,alloc_stat,&
                                             n_irreps
    character(len=12),allocatable              :: irrepname(:)

    real_format1 = '("    ",a," = ", f10.6:" # ",a)'
    real_format2 = '("    ",a," = ",es10.3:" # ",a)'

    allocate(irrepname(n_nonempty_irreps),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("occupation_write: allocation of irrepname failed")
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_spin=1_i4_kind
       n_irreps = symmetry_data_n_proj_irreps()
       do i=1,n_nonempty_irreps
          irrepname(i) = symmetry_data_irrepname_proj(i)
       enddo
    else
       if (options_spin_restricted() ) then
          n_spin=1_i4_kind
       else
          n_spin=2_i4_kind
       endif
       n_irreps = symmetry_data_n_irreps()
       do i=1,n_nonempty_irreps
          irrepname(i) = symmetry_data_irrepname(i)
       enddo
    endif


!   hidden input parameter:
!   diagonal_offset

    call start("OCCUPATION","OCCUPATION_WRITE", &
         iounit,operations_echo_input_level)
    call real("CHARGE            ",charge         ,df_charge        ,1)
    call flag("FIXED_SPIN_DIFF   ",fixed_spin_diff,df_fixed_spin_diff)
    call real("MAGN_MOMENT       ",magn_moment    ,df_magn_moment   ,1)
    call real("MIN_SPIN_DIFF     ",min_spin_diff  ,df_min_spin_diff ,2)
    call real("DIAG_OFFSET_EV    ",diag_offset_ev ,df_diag_offset_ev,1)
    call flag("FIXED             ",fixed          ,df_fixed           )
    call intg("N_NONEMPTY_IRREPS ",n_nonempty_irreps,df_n_nonempty_irreps)
    call flag("HOLE_LOCALIZATION ",hole_localization,df_hole_localization )
    call flag("HOLE_UPDATE       ",hole_update,df_hole_update )
    call intg("N_HOLES           ",n_holes,df_n_holes                )
    call stop(empty_line=.false.)

    1000 format(" #    Fixed Occupation Numbers      #")
    1001 format( F8.3,3x,"  #",2x,a12,2x,"Irrep No. ",i3)
    1002 format(2F8.3,3x,"  #",2x,a12,2x,"Irrep No. ",i3)
    1011 format(" #", F8.3,3x,"  #",2x,a12,2x,"Irrep No. ",i3)
    1012 format(" #",2F8.3,3x,"  #",2x,a12,2x,"Irrep No. ",i3)
    if (fixed) then
       write(iounit,1000)
       if (n_spin==2_i4_kind) then
          do i=1,n_nonempty_irreps
             write(iounit,1002) noc_start(1,i),noc_start(2,i), &
                                irrepname(i),i
          enddo
       else
          do i=1,n_nonempty_irreps
             write(iounit,1001)noc_start(1,i),irrepname(i),i
          enddo
       endif
       write(iounit,'()')
    elseif (operations_echo_input_level /= echo_level_short) then
       write(iounit,1000)
       if (n_spin==2_i4_kind) then
          do i=1,n_nonempty_irreps
             write(iounit,1012)0.0,0.0,irrepname(i),i
          enddo
       else
          do i=1,n_nonempty_irreps
             write(iounit,1011)0.0,irrepname(i),i
          enddo
       endif
       write(iounit,'()')
    endif

    2000 format(" #    List of Holes                 #")
    2101 format(" # index   irrep   occ  ")
    2102 format(" # index   irrep   spin    occ ")
    2001 format(2(i5,3x),f6.3,"  # Hole No. ",i3)
    2002 format(3(i5,3x),f6.3,"  # Hole No. ",i3)
    2011 format(" # ",2(i5,3x),f6.3,"  # Hole No. ",i3)
    2012 format(" # ",3(i5,3x),f6.3,"  # Hole No. ",i3)
    if (fixed_hole) then
       write(iounit,2000)
       if (options_spin_restricted()) then
          write(iounit,2101)
       else
          write(iounit,2102)
       endif
       if(options_spin_restricted()) then
          do i=1,n_holes
             write(iounit,2001)hole_list(i),hole_irrep(i),hole_occnum(i),i
          enddo
       else
          do i=1,n_holes
             write(iounit,2002)hole_list(i),hole_irrep(i),hole_spin(i), &
                               hole_occnum(i),i
          enddo
       endif
       write(iounit,'()')
    elseif (operations_echo_input_level /= echo_level_short) then
       write(iounit,2000)
       if (options_spin_restricted()) then
          write(iounit,2101)
       else
          write(iounit,2102)
       endif
       if(options_spin_restricted()) then
          write(iounit,2011)1_i4_kind,1_i4_kind,0.0_r8_kind,1_i4_kind
       else
          write(iounit,2012)1_i4_kind,1_i4_kind,1_i4_kind,0.0_r8_kind,1_i4_kind
       endif
       write(iounit,'()')
    endif
    deallocate(irrepname,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("occupation_write: deallocation of irrepname failed")

  end subroutine occupation_write
  !*************************************************************

  !*************************************************************
  subroutine count_electrons()
    ! count the valence and core electrons and checks the number of
    ! occupied orbitals (in case of fixed occupations)
    ! informations in unique_atom_module are used and must exist before
    !** End of interface *****************************************
    use unique_atom_module, only: N_unique_atoms, unique_atoms
    use input_module
    integer :: i_ua
    !------------ Executable code --------------------------------

    n_elec      = 0.0_r8_kind
    n_core_elec = 0.0_r8_kind
    do i_ua = 1,N_unique_atoms
       n_elec      = n_elec      + real(unique_atoms(i_ua)%N_equal_atoms,r8_kind) &
                                 * unique_atoms(i_ua)%Z
       n_core_elec = n_core_elec + real(unique_atoms(i_ua)%N_equal_atoms &
                                 * unique_atoms(i_ua)%Zc,r8_kind)
    enddo

    ! charge is the electrostatic charge while n_elec counts electrons, thus
    n_elec = n_elec - charge

    if (n_core_elec > n_elec) call input_error( "occupation_read: &
         &number of core electrons exceeds total number of electrons" )
    if ( operations_core_density ) then
         if (n_core_elec == 0.0_r8_kind) call input_error( "occupation_read: &
              &Zc must be positive to create a core density." )
    else
         ! in case of a pseudopot calc n_elec holds number of valence electrons
         ! (otherwise n_core_elec is zero anyhow)
         n_elec = n_elec - n_core_elec
    endif

    if (fixed) then
       if ( abs(sum(noc_start)-n_elec) > 1.0e-10_r8_kind) call input_error( "occupation_read: &
            &your fixed occupation numbers don't sum up to the total &
            &number of electrons")
    endif

  end subroutine count_electrons
  !*************************************************************

  !*************************************************************
  subroutine get_n_elec(n_dummy)
    !  Purpose: make the private variable N_ELEC accessible
    !           very simple - very stupid.
    !           FN
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(out ) :: n_dummy
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    n_dummy = n_elec

  end subroutine get_n_elec
  !*************************************************************

  !*************************************************************
  subroutine get_n_core_elec(n_dummy)
    !  Purpose: make the private variable N_CORE_ELEC accessible
    !           very simple - very stupid.
    !           FN
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(out ) :: n_dummy
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    n_dummy = n_core_elec

  end subroutine get_n_core_elec
  !*************************************************************

  !*************************************************************
  subroutine get_spin_diff(spin_diff)
    !  Purpose: make private variable magn_moment accessible
    !           UB
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(out ) :: spin_diff
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    spin_diff = magn_moment

  end subroutine get_spin_diff
  !*************************************************************

  !*************************************************************
  subroutine get_df_spin_diff(df_spin_diff)
    !  Purpose: make private variable df_magn_moment accessible
    !           TS
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(out ) :: df_spin_diff
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    df_spin_diff = df_magn_moment

  end subroutine get_df_spin_diff
  !*************************************************************

  !*************************************************************
  subroutine put_magn_moment(fermi_unp)
    !  Purpose: set private variable magn_moment
    !           TS
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in) :: fermi_unp
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    magn_moment = fermi_unp

  end subroutine put_magn_moment
  !*************************************************************

  !*************************************************************
  subroutine get_fix_spin_switch(fixed_spin)
    !  Purpose: make private variable fixed_spin_diff accessible
    !           TS
    !------------ Declaration of formal parameters ---------------
    logical, intent(out ) :: fixed_spin
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    fixed_spin = fixed_spin_diff

  end subroutine get_fix_spin_switch
  !*************************************************************

  !*************************************************************
  subroutine get_df_fix_spsw(df_fixed_spin)
    !  Purpose: make private variable df_fixed_spin_diff accessible
    !           TS
    !------------ Declaration of formal parameters ---------------
    logical, intent(out ) :: df_fixed_spin
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    df_fixed_spin = df_fixed_spin_diff

  end subroutine get_df_fix_spsw
  !*************************************************************

  !*************************************************************
  subroutine put_fix_spin_switch(fermi_fix_ud)
    !  Purpose: set private variable fixed_spin_diff
    !           TS
    !------------ Declaration of formal parameters ---------------
    logical, intent(in) :: fermi_fix_ud
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    fixed_spin_diff = fermi_fix_ud

  end subroutine put_fix_spin_switch
  !*************************************************************

  !*************************************************************
  subroutine get_diagonal_offset(offset)
    !  Purpose: make private variable diagonal_offset accessible
    !           conversion from eV (input) to Ht (ham_calc) !
    !           UB
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(out ) :: offset
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    offset = diag_offset_ev / convert1

  end subroutine get_diagonal_offset
  !*************************************************************

  !*************************************************************
  real(kind=r8_kind) function get_min_spin_diff()
    !  Purpose: make private variable min_spin_diff accessible
    !           UB
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    get_min_spin_diff = min_spin_diff

  end function get_min_spin_diff
  !*************************************************************

  !*************************************************************
  real(kind=r8_kind) function get_charge()
    ! Purpose: make private variable charge accessible
    !          AS (04.02)
    !** End of interface *****************************************

    !------------ Executable code --------------------------------
    get_charge=charge

  end function get_charge
  !*************************************************************

  !*************************************************************
  subroutine reoccup ()
    !  Purpose: re-occupy the levels after the eigenvalues
    !           and eigenvectors are ready, i.e. after scf.
    !           (1) First a list with successive numbers of the
    !               eigenvalues is established:
    !               num_list(i)%m(m,is) is the number of
    !               orbital(i,m,is) with eigenvalue eigval(i)%m(m,is)
    !           (2) Then the electrons are (currently) filled in
    !               according to the AufbauPrinciple using the
    !               above list.
    !
    !  Warning:
    !  In this subroutine n_occo(i_g)%m(i_o,i_s) is loaded to hold the
    !  number of partially occupied orbitals. Actually, however, this
    !  variable should hold the highest orbital index of all occupied
    !  orbitals. This difference has to be considered as soon as core
    !  holes are introduced. (UB 7/97)
    !
    !  Subroutine called by: main_scf
    !
    !  Author: Folke Noertemann
    !  Date: 8/95
    !** End of interface *****************************************

    !------------ Modules ----------------------------------------
    use machineparameters_module, ONLY : machineparameters_DimCheck
    use print_module   ! nice interface for printout
    use eigen_data_module, only: eigval,n_holes_per_irrep,holes_per_irrep,&
         spin_of_hole
    use output_module, only: output_reoccup
    use options_module, only: options_spin_orbit
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: i,m,n,is,degen,num_orb,index, &
         n_elec_int,num_spin,i_hole, alloc_stat
    logical :: hole
    real(kind=r8_kind)  :: occ_hole,remmel,remmel_core
    real(kind=r8_kind),allocatable       :: remmel_m(:)
    ! remmel: remaining electrons for all variants of
    !         occupation except fixed_spin_diff
    ! remmel_m(:): rem. el. on either spin side in
    !              case of fixed_spin_diff
    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    integer(kind=i4_kind),allocatable    :: n_partners(:)
    ! n_irrep   : number of irreps
    ! dim_irrep : number of independent functions in irrep

    ! External an intrisic routines
    intrinsic max,nint
    external error_handler
    !------------ Executable code --------------------------------

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep),n_partners(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
          n_partners(n) = ssym%partner_proj(n)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep),n_partners(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
          n_partners(n) = ssym%partner(n)
       enddo
    endif ! options_spin_orbit

    if (.not. allocated (n_occo)) then
       call alloc_n_occo (ssym)
    endif
    if (.not. allocated (occ_num)) then
       call alloc_occ_num()
    endif
    call init(occ_num)
    call init(n_occo)
    if (operations_core_density) then
       call init(occ_num_core)
       call init(n_occo_core)
    endif

    ! first see if dimension are right
    call check()


    call occupation_level_sort(eigval)
    ! ---------------------------------------------------



    !---- now fill in the electrons according to  ----
    !---- Aufbau-Principle


    ! n_elec is of type REAL. This routine was originally
    ! designed to handle it as an integer.
    ! As long as this design is not changed we introduce
    ! n_elec_int as a help variable here
    n_elec_int = nint(n_elec,kind=i4_kind)

    num_orb=0_i4_kind
    do i=1,n_irrep
       num_orb=num_orb+dim_irrep(i)*n_partners(i)
    enddo
    num_orb=num_orb * 2_i4_kind ! total number of orbitals
    if (n_elec_int .gt. num_orb) call error_handler&
         ( "REOCCUP: more electrons than orbitals" )

    if (fixed_spin_diff) then

       allocate(remmel_m(2_i4_kind),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("reoccup: allocation of remmel_m failed")
       ! set counters of remaining electrons for each spin direction
       remmel_m(1) = (n_elec-magn_moment)/2.0_r8_kind ! Spin 1 = minority
       remmel_m(2) = n_elec-remmel_m(1)               ! Spin 2 = majority
    endif

    if (options_spin_orbit) then
       num_spin = 1
    else
       num_spin = 3_i4_kind - ssym%n_spin
    endif


    ! get the actual work done:
    if ( force_spectrum ) then
       call force(ssym%n_spin, n_occo, occ_num)
    else
       call aufbau()
    endif

    if (output_reoccup) then
       write(output_unit,*)"REOCCUP: Indices of highest occupied MOs per Irrep and Spin"
       write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
       do i=1,symmetry_data_n_irreps()
          write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo(is,i),is=1,ssym%n_spin)
       enddo
       write(output_unit,*)" "
       if(operations_core_density)then
          write(output_unit,*)"REOCCUP: Indices of highest occupied core orbital per Irrep and Spin"
          write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
          do i=1,symmetry_data_n_irreps()
             write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo_core(is,i),is=1,ssym%n_spin)
!:TST[ occ_num_core
!            do m=1,maxval(n_occo_core(:,i))
!               write(output_unit,*)(occ_num_core(i)%m(m,is),is=1,ssym%n_spin)
!            enddo
!:TST]
          enddo
          write(output_unit,*)" "
       endif
    endif

    ! deallocate appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    deallocate(dim_irrep,n_partners)
    if (fixed_spin_diff) then
       deallocate(remmel_m)
    endif
    return

  contains

    subroutine force(n_spin, n_occo, occ_num)
      implicit none
      integer(i4_kind), intent(in)  :: n_spin
      integer(i4_kind), intent(out) :: n_occo(:,:) ! (n_spin,:)
      type(arrmat2), intent(inout) :: occ_num(:) ! (n_irrep)%(dim_irrep, n_spin)
      ! *** end of interface ***

      integer(i4_kind) :: irr, is, degen
      integer(i4_kind) :: spin_degen
      real(r8_kind)    :: occ

      DPRINT 'reoccup.force: entered'

      select case (n_spin)
      case (1)
        spin_degen = 2
      case (2)
        spin_degen = 1
      case default
        ! Make compiler happy:
        spin_degen = -1
        ABORT('n_spin')
      end select

      n_occo(:,:) = 0
      ! copy occupation numbers:
      spin : do is = 1, n_spin
      irrep: do irr = 1, n_irrep
         degen = spin_degen * n_partners(irr)
         DPRINT 'spin=', is, 'irrep=', irr, 'degen=', degen, 'dim=',dim_irrep(irr)
         occ_num(irr)%m(:, is) = 0.0
         do m = 1, dim_irrep(irr)
            if (m > size(occupation_numbers(is,irr)%occ)) cycle !irrep

            occ = occupation_numbers(is,irr)%occ(m)
            n_occo(is, irr) = n_occo(is, irr) + 1
            DPRINT 'reoccup.force: put ', occ, ' to (n,irr,i)=', m, irr
            ASSERT(degen-occ>=-1.0E-7)
            occ_num(irr)%m(m, is) = occ
         enddo
         DPRINT 'N_OCCO(', is, irr,')=', n_occo(is, irr)
      enddo irrep
      enddo spin
    end subroutine force

    subroutine aufbau()
      implicit none
      ! *** end of interface ***
      if (fixed) then
         noc_fixed=noc_start
      endif
      !    remmel = n_elec_int   ! remmel = REMaining MELectrons
      remmel = n_elec
      if (operations_core_density) remmel_core = n_core_elec
    if(fixed_hole) then
    pre_mo : do index=1,num_orb
       pre_irrep: do i=1,n_irrep
          pre_spin: do is=1,ssym%n_spin
             pre_orbital: do m=1,dim_irrep(i)
                if(num_list(i)%m(m,is).eq.index) then
                   degen = num_spin*n_partners(i)

                      hole=.false.
                      pre_holes: do i_hole=1,n_holes_per_irrep(i)
                         if (holes_per_irrep(i)%m(i_hole).eq.m .and. &
                              spin_of_hole(i)%m(i_hole).eq.is) then
                            if (output_reoccup) then
                               write(output_unit,'("Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                                    m,i,eigval(i)%m(m,is)*convert1
                               write(stdout_unit,'("Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                                    m,i,eigval(i)%m(m,is)*convert1
                            endif
                            hole=.true.
                            occ_hole = occnum_hole(i)%m(i_hole)
                            if (remmel<occ_hole.and.occ_hole>0.0_r8_kind) then
                            call error_handler("reoccup: not enough electrons left to occupy hole")

                      else
                            occ_num(i)%m(m,is) = occ_hole
                            remmel = remmel - occ_hole
                         if(m.gt.n_occo(is,i)) n_occo(is,i) = m
                      endif
                      cycle pre_mo
                      endif
                      enddo pre_holes
               endif ! if (num_list .eq. index)
             enddo pre_orbital
          enddo pre_spin
       enddo pre_irrep
    enddo pre_mo
    endif
      mo : do index=1,num_orb
         irrep: do i=1,n_irrep
            spin: do is=1,ssym%n_spin
               orbital: do m=1,dim_irrep(i)
                  if(num_list(i)%m(m,is).eq.index) then
                     degen = num_spin*n_partners(i)
                     if (operations_core_density .and. &
                          remmel_core > 0.0_r8_kind) then
                        if (remmel_core > real(degen,kind=r8_kind)) then
                           occ_num_core(i)%m(m,is) = real(degen,kind=r8_kind)
                           n_occo_core(is,i) = n_occo_core(is,i) + 1_i4_kind
                           remmel_core = remmel_core - real(degen,kind=r8_kind)
                        else
                           occ_num_core(i)%m(m,is) = remmel_core
                           n_occo_core(is,i) = n_occo_core(is,i) + 1_i4_kind
                           remmel_core = 0.0_r8_kind
                        endif
                     endif
                     if (fixed) then
                        if (sum(noc_fixed)<=0.0_r8_kind) then
                           ! no irrep is left to be fllied any more
                           exit mo
                           elseif (noc_fixed(is,i)<=0.0_r8_kind) then
                              ! This irrep is filled completely according
                              ! to start occupation
                           cycle mo
                           elseif (noc_fixed(is,i)<=degen) then
                              ! There are electrons left to fill this irrep
                              ! but the actual mo can take more
                           occ_num(i)%m(m,is) = noc_fixed(is,i)
                           noc_fixed(is,i)=0.0_r8_kind
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           cycle mo
                        else
                           occ_num(i)%m(m,is) = real(degen,kind=r8_kind)
                           noc_fixed(is,i)=noc_fixed(is,i)-real(degen,kind=r8_kind)
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           cycle mo
                        endif! if (fixed)

                        elseif (fixed_hole) then

                           !find out if current orbital (i,m,is) was identified
                           !as a hole
                           ! In this implementation we will not notice if one
                           ! orbital is identified as hole TWICE.
                        hole=.false.
                        holes: do i_hole=1,n_holes_per_irrep(i)
                           if (holes_per_irrep(i)%m(i_hole).eq.m .and. &
                                spin_of_hole(i)%m(i_hole).eq.is) then
                              if (output_reoccup) then
                                 write(output_unit,'("Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                                      m,i,eigval(i)%m(m,is)*convert1
                                 write(stdout_unit,'("Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                                      m,i,eigval(i)%m(m,is)*convert1
                              endif
                              hole=.true.
                              occ_hole = occnum_hole(i)%m(i_hole)
                              exit holes
                           endif
                        enddo holes

                        if ( remmel<=0.0_r8_kind ) then
                           exit mo
                           elseif ( remmel<=real(degen,kind=r8_kind) ) then
                           if (.not.hole) then

                              occ_num(i)%m(m,is) = real(remmel,kind=r8_kind)
                              remmel = 0.0_r8_kind
                              if(m.gt.n_occo(is,i))   n_occo(is,i)=m
                           endif

                        else
                           if (.not.hole) then

                              occ_num(i)%m(m,is) = real(degen,kind=r8_kind)
                              remmel = max(remmel-real(degen,kind=r8_kind),0.0_r8_kind)
                              if(m.gt. n_occo(is,i))  n_occo(is,i) = m
                           endif

                        endif

                        elseif (fixed_spin_diff) then
                        if (remmel_m(1) .le. 0.0_r8_kind .and. remmel_m(2) .le. 0.0_r8_kind) then
                           exit mo
                        endif
                        if(remmel_m(is) .le. 0.0_r8_kind) then
                           cycle mo
                           elseif(remmel_m(is) .le. real(degen,kind=r8_kind)) then
                           occ_num(i)%m(m,is) = remmel_m(is)
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           remmel_m(is) = 0.0_r8_kind
                           cycle mo
                        else
                           occ_num(i)%m(m,is) = real(degen,kind=r8_kind)
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           remmel_m(is) = remmel_m(is)-real(degen,kind=r8_kind)
                        endif

                     else ! reoccupation according to Aufbau without holes
                        if(remmel.le.0.0_r8_kind) then  ! no more electrons left
                           exit mo
                           elseif(remmel.le.real(degen,kind=r8_kind)) then
                              ! all remaining electrons are
                              ! put into the current orbital
                           occ_num(i)%m(m,is) = remmel
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           remmel = 0.0_r8_kind
                           exit mo
                        else
                           ! regular case: there are more remaining
                           ! electrons than vacancies in the current orbital
                           occ_num(i)%m(m,is) = real(degen,kind=r8_kind)
                           n_occo(is,i) = n_occo(is,i) + 1_i4_kind
                           remmel = max((remmel-real(degen,kind=r8_kind)),0.0_r8_kind)
                        endif
                     endif! if (fixed)
                     cycle mo
                  endif! if (num_list .eq. index)
               enddo orbital
            enddo spin
         enddo irrep
      enddo mo
    end subroutine aufbau

    subroutine check()
      implicit none
      ! *** end of interface ***

      integer(i4_kind) :: i

      if(machineparameters_DimCheck) then
         if (n_irrep.ne.size(eigval,1)) &
              call error_handler( "reoccup : dim eigval wrong" )
         if (n_irrep.ne.size(occ_num,1)) &
              call error_handler( "reoccup : dim occ_num wrong")
         if(n_irrep.ne.size(n_occo,2))   &
              call error_handler( "reoccup : dim n_occo wrong")
         if(ssym%n_spin.ne.size(n_occo,1))   &
              call error_handler( "reoccup : dim n_occo wrong")
         if (operations_core_density) then
            if (n_irrep.ne.size(occ_num_core,1)) &
                 call error_handler( "reoccup : dim occ_num_core wrong")
            if(n_irrep.ne.size(n_occo_core,2))   &
                 call error_handler( "reoccup : dim n_occo_core wrong")
            if(ssym%n_spin.ne.size(n_occo_core,1))   &
                 call error_handler( "reoccup : dim n_occo_core wrong")
         endif

         do i=1,n_irrep
            if (dim_irrep(i).ne.size(eigval(i)%m,1)) &
                 call error_handler( "reoccup : dim (2) eigval wrong")
            if (ssym%n_spin.ne.size(eigval(i)%m,2)) &
                 call error_handler( "reoccup : dim (3) eigval wrong")
            if (dim_irrep(i).ne.size(occ_num(i)%m,1))  &
                 call error_handler( "reoccup : dim occ_num wrong")
            if (ssym%n_spin.ne.size(occ_num(i)%m,2))  &
                 call error_handler( "reoccup : dim occ_num wrong")
            if (operations_core_density) then
               if (dim_irrep(i).ne.size(occ_num_core(i)%m,1))  &
                    call error_handler( "reoccup : dim occ_num_core wrong")
               if (ssym%n_spin.ne.size(occ_num_core(i)%m,2))  &
                    call error_handler( "reoccup : dim occ_num_core wrong")
            endif
         enddo

      endif
    end subroutine check


  end subroutine reoccup
  !*************************************************************

  !*************************************************************
  subroutine read_ocup()
    ! purpose: read the variable OCUP produced by the old lcgto
    !          from the file $TTFSSTART/occ_num.dat and
    !          re-structure it into the variable occ_num.
    !          Just for testing purposes.
    !
    ! subroutine called by: main_scf
    !
    !  Author: Folke Noertemann
    !  Date: 8/95
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use filename_module, only: tmpfile
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)          :: io_u,num_orb, &
         icount,m,is,i
    real(kind=r8_kind),allocatable :: arr_help(:)

    ! first calculate number of orbitals
    num_orb=0
    do i=1,ssym%n_irrep
       num_orb=num_orb+ssym%dim(i)
    enddo
    num_orb=num_orb*ssym%n_spin

    ! allocate the help array
    allocate(arr_help(num_orb))

    ! now read arr_help from file
    io_u=get_iounit()
    open(io_u,form='formatted',status='old',&
         file=trim(tmpfile('occ_num.dat')))

    read(io_u,*)arr_help

    close(io_u)
    call return_iounit(io_u)


    ! take care of occ_nums allocation
    if(.not.allocated(occ_num)) &
         call alloc_occ_num()
    if(.not.allocated(n_occo)) &
         call alloc_n_occo(ssym)

    ! now remap to occ_num and set n_occo to the index of the HOMO
    icount=1
    do i=1,ssym%n_irrep
       do is=1,ssym%n_spin
          do m=1,ssym%dim(i)
             occ_num(i)%m(m,is)=arr_help(icount)
             icount=icount+1
             if (occ_num(i)%m(m,is) /= 0.0_r8_kind) n_occo(is,i)=m
          enddo
       enddo
    enddo

    deallocate(arr_help)
  end subroutine read_ocup
  !*************************************************************

  !*************************************************************
  subroutine occupation_level_sort(eigval,separate_spin)
    !  Purpose: produce the list num_list which contains
    !           the index of the corresponding eigenvalue
    !           in ascending order.
    !           If levels are to be sorted for each spin
    !           separately, the num_list(i)%m(m,1) contains
    !           indices for spin 1 and num_list(i_%m(m,2)
    !           contains indices for spin 2.
    !
    !  subroutine called by: 'reoccup', 'fermi'
    !  Author: FN
    !  Date: 4/96
    !------------ Modules ---------------------------------------
    use machineparameters_module, ONLY : machineparameters_DimCheck
    use options_module, only: options_spin_orbit
    implicit none
    type(arrmat2)         :: eigval(:)
    logical,optional      :: separate_spin
    !** End of interface *****************************************

    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep : number of independent functions in irrep
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)   :: i,is,m,n,mm,mstart,istart, &
         ii,iss,icount,alloc_stat
    logical                 :: separate
    external error_handler


    if (.not.present(separate_spin)) then
       separate = .false.
    endif  ! this is default value: no separate sorting

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

    ! see if dimensions are ok
    if(machineparameters_DimCheck) then
       if (n_irrep.ne.size(eigval,1)) &
            call error_handler &
            ( "occupation_level_sort : dim 1 eigval wrong" )
       do i=1,n_irrep
          if(dim_irrep(i).ne.size(eigval(i)%m,1)) &
               call error_handler &
               ( "occupation_level_sort : dim 2 eigval wrong" )
       enddo
    endif

    ! allocate the num_list which contains the indices of
    ! the sorted levels
    if ( .not.check_num_list) then
       allocate(num_list(n_irrep),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       do i = 1,n_irrep
          allocate(num_list(i)%m(dim_irrep(i),ssym%n_spin),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       enddo
       check_num_list =.true.
    endif

    ! first reset the variable num_list
    do i = 1, size(num_list) ! n_irrep
       num_list(i)%m = 0
    enddo

    ! --- fast and correct sorting algorithm ----
    if (.not.separate) then
       num_list(1)%m(1,1) = 0_i4_kind
 DPRINT 'occupations module: case not separate'

 do i = 1,n_irrep           ! Loop through all orbitals
    do is = 1,ssym%n_spin
       do m = 1,dim_irrep(i)


          istart = is
          icount = 1
          mstart = m

          do ii = i,1,-1              ! Loop just through those
             do iss = istart,1,-1     ! orbitals whose index is BELOW
                do mm = mstart,1,-1   ! orbitals whose index is BELOW

                   ! (1) examine how many of those orbitals are BELOW
                   !     the one picked by the first loop =
                   !     preliminary number num_list
                   ! (2) if eigval(orbital) larger or equal that one of
                   !     the orbital picked by the first loop,
                   !     num_list(orbital) must be larger than evaluated
                   !     by the last turn
                   if(eigval(ii)%m(mm,iss).lt.eigval(i)%m(m,is)) then
                      icount = icount + 1
                   else
                      num_list(ii)%m(mm,iss) = num_list(ii)%m(mm,iss) + 1_i4_kind
                   endif
                enddo
                mstart = dim_irrep(ii)
             enddo
             istart = ssym%n_spin
             if(ii.ne.1) mstart = dim_irrep(ii-1)
          enddo
          num_list(i)%m(m,is) = icount


       enddo
    enddo
 enddo


    else
       DPRINT 'occupations module: case separate'
       num_list(1)%m(1,1) = 0_i4_kind
       num_list(1)%m(1,2) = 0_i4_kind

       do is = 1,ssym%n_spin
          icount = 1
          do i = 1,n_irrep           ! Loop through all orbitals
             do m = 1,dim_irrep(i)

                mstart = m
                icount = 1


                do ii = i,1,-1              ! Loop just through those
                   do mm = mstart,1,-1      ! orbitals whose index is BELOW

                      ! (1) examine how many of those orbitals are BELOW
                      !     the one picked by the first loop =
                      !     preliminary number num_list
                      ! (2) if eigval(orbital) larger or equal that one of
                      !     the orbital picked by the first loop,
                      !     num_list(orbital) must be larger than evaluated
                      !     by the last turn
                      if(eigval(ii)%m(mm,is).lt.eigval(i)%m(m,is)) then
                         icount = icount + 1
                      else
                         num_list(ii)%m(mm,is) = num_list(ii)%m(mm,is) + 1_i4_kind
                      endif
                   enddo
                   if(ii.ne.1) mstart = dim_irrep(ii-1)
                enddo
                num_list(i)%m(m,is) = icount
             enddo
          enddo
       enddo

    endif

    deallocate(dim_irrep,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler( &
         "occupation_level_sort: deallocation of dim_irrep failed")
  end subroutine occupation_level_sort
  !*************************************************************

  !*************************************************************
  real(kind=r8_kind) function occupation_spindiff(signum)
  ! Purpose: calculate the spin difference.
  !
  ! Function called by: mainscf
  ! Author: FN
  ! Date  : 9/96
  !** End of interface *****************************************
  !------------ Declaration of local variables -----------------
  logical, intent(in), optional :: signum
  logical               :: sign_intern
  real(kind=r8_kind)    :: sum_up,sum_down
  integer(kind=i4_kind) :: i_gamma,i_m
  !------------ Executable code --------------------------------
  sign_intern=.false.
  if(present(signum)) sign_intern=signum

  if (ssym%n_spin.ne.2) then
     occupation_spindiff = 0.0_r8_kind
     return
  endif
  sum_up=0.0_r8_kind
  sum_down=0.0_r8_kind
  do i_gamma=1,ssym%n_irrep
     do i_m=1,n_occo(1,i_gamma)
        sum_up = sum_up + occ_num(i_gamma)%m(i_m,1)
     enddo
     do i_m=1,n_occo(2,i_gamma)
        sum_down = sum_down + occ_num(i_gamma)%m(i_m,2)
     enddo
  enddo
  occupation_spindiff = sum_up - sum_down
  if(.not. sign_intern) occupation_spindiff = abs(occupation_spindiff)
  end function occupation_spindiff
  !*************************************************************

  !*************************************************************
  subroutine occupation_2d_correct()
    !
    ! In  the case  of pseudo  2D  irreps occupy  each pseudo  partner
    ! equally.
    !
    ! Function called by: mainscf
    ! Author: MS
    ! Date  : 8/97
    !
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: i_gamma, i_m, i_s, counter

    do i_gamma = 1, ssym % n_irrep
       if (ssym % pseudo(i_gamma)) then
          do i_m = 1, ssym % dim(i_gamma), 2
             occ_num(i_gamma) % m(i_m, :) = 0.5_r8_kind * &
                  (occ_num(i_gamma) % m(i_m, :) + occ_num(i_gamma) % m(i_m + 1, :))
             occ_num(i_gamma) % m(i_m + 1, :) = occ_num(i_gamma) % m(i_m, :)
          enddo
          ! now correct n_occo
          do i_s = 1, ssym % n_spin
             counter = 0
             do i_m = 1, ssym % dim(i_gamma)
                if (occ_num(i_gamma) % m(i_m, i_s) > 0.0_r8_kind) counter = counter + 1
             end do
             n_occo(i_s, i_gamma) = counter
          end do
       end if
    enddo
  end subroutine occupation_2d_correct
  !*************************************************************

  !*************************************************************
  subroutine print_occ_num(loop)
    !  Purpose: print out the occupation number in a well-suited
    !           format to the file $TTFSOUT/occupation.out.
    !           The format is supposed to be comparable to
    !           the occupation numbers of the old lcgto
    !  subroutine called by: 'main_scf'
    !  Author: Folke Noertemann
    !  Date: 8/95
    !------------ Modules ---------------------------------------
    use filename_module, only: outfile
    use options_module, only: options_spin_orbit
    integer(kind=i4_kind),intent(in),optional :: loop
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)   :: io_u,i_gamma,m,is,i_dim
    ! number and dimensions of irreps (spin orbit and vector case)
    integer(kind=i4_kind)               :: n_irrep
    integer(kind=i4_kind),allocatable   :: dim_irrep(:)

    io_u=get_iounit()
    open(io_u,status='unknown',position='append',&
         form='formatted', &
         file=trim(outfile('occupation.out')))

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = symmetry_data_n_proj_irreps()
       allocate(dim_irrep(n_irrep))
       do m=1,n_irrep
          dim_irrep(m) = symmetry_data_dimension_proj(m)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = symmetry_data_n_irreps()
       allocate(dim_irrep(n_irrep))
       do m=1,n_irrep
          dim_irrep(m) = symmetry_data_dimension(m)
       enddo
    endif ! options_spin_orbit

    do i_gamma=1,n_irrep
       i_dim=dim_irrep(i_gamma)
       write(io_u,*)' '
       write(io_u,*)' '
       if(present(loop)) then
          write(io_u,*)' Loop  :',loop
          write(io_u,*)' '
       endif

       write(io_u,*)'---------- Irrep ',i_gamma,' -------------'

       do is=1,symmetry_data_n_spin()
          write(io_u,*)' Spin : ',is
          write(io_u,*)' '
          if (allocated(n_occo)) then
             write(io_u,*)'Index of the highest occupied orbital :', &
                           n_occo(is,i_gamma)
             write(io_u,*)' '
          endif

          write(io_u,'((A,F10.4))')&
               (' ',occ_num(i_gamma)%m(m,is),m=1,i_dim)

       enddo

    enddo

    ! deallocate dimensions of irreps
    deallocate(dim_irrep)

    close(io_u)
    call return_iounit(io_u)

  end subroutine print_occ_num
  !*************************************************************

  !*************************************************************
  subroutine occupation_print_spectrum()
    ! Purpose : provide a pretty output for the occupation
    !           spectrum. Should be in close analogy to old
    !           LCGTO.
    ! Author: FN
    ! Date: 10/96
    !
    use comm, only: comm_rank
    use eigen_data_module, only: eigval
    use options_module, only: options_spin_orbit
    !** End of interface *****************************************


    integer(i4_kind) :: i, i_count, m, is, counter_1, counter_2, &
         num_orb, num_print, num_occ
    intrinsic min

    ! FIXME: sigsegv?
    if (comm_rank() /= 0) return

    num_orb=0
    num_occ=0

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i=1,ssym%n_proj_irrep
          num_orb = num_orb + ssym%dim_proj(i)
          num_occ =num_occ + n_occo(1,i)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       do i=1,ssym%n_irrep
          num_orb = num_orb + ssym%dim(i)
          do is=1,ssym%n_spin
             num_occ =num_occ + n_occo(is,i)
          enddo
       enddo
    endif ! options_spin_orbit
    num_orb = num_orb*symmetry_data_n_spin()
    write(output_unit,*)'num_orb =',num_orb,'        num_occ=',num_occ
    num_print = min(num_occ+40,num_orb)

    if (symmetry_data_n_spin().eq.1) then
       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          call occupation_level_sort(eigval)
          write(output_unit,*)" "
          write(output_unit,*)" --------------- Spectrum ----------------"
          write(output_unit,*)" "
          counter_1 = num_print
          do i_count=num_print,1,-1
             do i=1,ssym%n_proj_irrep
                element : do m=1,ssym%dim_proj(i)
                   if (num_list(i)%m(m,1).eq.i_count) then
                      write(output_unit,1050)counter_1,m,ssym%name_proj(i),&
                           occ_num(i)%m(m,1),&
                           eigval(i)%m(m,1)*convert1,eigval(i)%m(m,1)
                      counter_1 = counter_1 - 1
                      cycle element
                   endif
                enddo element
             enddo
          enddo
       else ! options_spin_orbit
          !
          ! STANDARD SCF (NO SPIN ORBIT)
          !
          call occupation_level_sort(eigval)
          write(output_unit,*)" "
          write(output_unit,*)" --------------- Spectrum ----------------"
          write(output_unit,*)" "
          counter_1 = num_print
          do i_count=num_print,1,-1
             do i=1,ssym%n_irrep
                element2 : do m=1,ssym%dim(i)
                   if (num_list(i)%m(m,1).eq.i_count) then
                      write(output_unit,1000)counter_1,m,ssym%name(i),i,&
                           occ_num(i)%m(m,1),&
                           eigval(i)%m(m,1)*convert1
                      counter_1 = counter_1 - 1
                      cycle element2
                   endif
                enddo element2
             enddo
          enddo
       endif ! options_spin_orbit
       write(output_unit,*)" -----------------------------------------"
       write(output_unit,*)" "

    else
       ! sort the levels again since with n_spin=2,
       ! there is the possibility of fermi_fix_ups_and_down
       ! which causes the variable num_list to contain
       ! the wrong things
       call occupation_level_sort(eigval)

       write(output_unit,*)" "
       write(output_unit,*)" ------------------------------------- Spectrum ------------------------------------- "

!first let us calculate counter_1 and counter_2
       counter_1 = 0
       counter_2 = 0
       do i_count=1, num_print
          do i=1,ssym%n_irrep
             do m=1,ssym%dim(i)
                do is=1,symmetry_data_n_spin()

                   if (num_list(i)%m(m,is).eq.i_count) then
                      if (is.eq.1) then
                         counter_1 = counter_1 + 1
                      else
                         counter_2 = counter_2 + 1
                      endif
                   endif

                enddo
             enddo
          enddo

       enddo

       do i_count=num_print,1,-1
          do i=1,ssym%n_irrep
             do m=1,ssym%dim(i)
                do is=1,symmetry_data_n_spin()

                   if (num_list(i)%m(m,is).eq.i_count) then
                      if (is.eq.1) then
                         write(output_unit,1100)&
                              counter_1,m,ssym%name(i),i,is,&
                              occ_num(i)%m(m,is),&
                              eigval(i)%m(m,is)*convert1
                         counter_1 = counter_1 - 1
                      else
                         write(output_unit,1200)&
                              counter_2,m,ssym%name(i),i,is, &
                              occ_num(i)%m(m,is), &
                              eigval(i)%m(m,is)*convert1
                         counter_2 = counter_2 - 1
                      endif
                   endif

                enddo
             enddo
          enddo

       enddo
       write(output_unit,*)" "
       write(output_unit,*)" ------------------------------------------------------------------------------------ "
       write(output_unit,*)" "
    endif

#ifdef FPP_COMPAT
1000 format(I4,"   ",I3,"   ",A4,I3,"   ",F7.4,"  ",F13.6)
1050 format(I4,"   ",I3,"   ",A6,"   ",F7.4,"  ",F13.6,"  ",F13.6)
#else
1000 format(I4,"   ",I3,"   ",A4,I3,"   ",F7.4,"  ",F16.6)
1050 format(I4,"   ",I3,"   ",A6,"   ",F7.4,"  ",F16.6,"  ",F13.6)
#endif
1100 format(I4,"   ",I3,"   ",A4,2I3,"   ",F7.4,"  ",F16.6)
1200 format("                                             ",I4,"   ",I3,"   ",A4,2I3,"   ",F7.4,"  ",F16.6)
    return
  end subroutine occupation_print_spectrum
  !*************************************************************

  !*************************************************************
  subroutine occupation_print_popspectrum(pop_col,pop_store,eig_min,eig_max)
    ! Purpose : provide a pretty output for the occupation
    !           spectrum. Should be in close analogy to old
    !           LCGTO.
    ! Author: FN
    ! Date: 10/96
    !** End of interface *****************************************
    ! -----------------------------------------------------
    use eigen_data_module, only: eigval
    use iounitadmin_module
    type(arrmat3),optional :: pop_col(:) ! populations
    type(pop_store_type), optional :: pop_store(:,:,:) ! populations, coeffs,...
    real(kind=r8_kind),optional :: eig_min
    real(kind=r8_kind),optional :: eig_max
    ! --------- Declaration of local variables ------------
    integer(i4_kind) :: i, i_count, m, is, counter, num_orb, &
         num_print, num_occ, i_frag, i_print, alloc_stat, &
         n_fragments, n_cont
    real(kind=r8_kind) :: eig_min_loc, eig_max_loc
    integer(kind=i4_kind), allocatable :: counter_ir(:)
    integer(kind=i4_kind),pointer      :: pointer_index(:)
    real(kind=r8_kind), pointer        :: pointer_pop(:), pointer_coeff(:)
    logical                :: pop_col_loc, pop_store_loc
    intrinsic min
    ! --------- Executable code ---------------------------

    if(present(pop_store)) then
       pop_store_loc=.true.
       allocate(counter_ir(ssym%n_irrep), stat=alloc_stat)
       if(alloc_stat/=0) &
            call error_handler('occupation_print_popspectrum: allocation counter_ir failed')
       counter_ir=1
       n_fragments = size (pop_store, 3)
    else
       pop_store_loc=.false.
       n_fragments = 0
    end if

    if(present(pop_col)) then
       pop_col_loc=.true.
    else
       pop_col_loc=.false.
    end if

    if(present(eig_min)) then
       eig_min_loc=eig_min
       eig_max_loc=eig_max
    else
       eig_min_loc=minval(eigval(1)%m)-10.0_r8_kind
       eig_max_loc=maxval(eigval(1)%m)+10.0_r8_kind
    end if
    num_orb=0
    num_occ=0
    do i=1,ssym%n_irrep
       num_orb = num_orb + ssym%dim(i)
       do is=1,ssym%n_spin
          num_occ =num_occ + n_occo(is,i)
       enddo
    enddo
    num_orb = num_orb*symmetry_data_n_spin()
    write(output_unit,*)'num_orb =',num_orb,'        num_occ=',num_occ
    num_print = min(num_occ+40,num_orb)

    call occupation_level_sort(eigval)
    write(output_unit,*)" "
    write(output_unit,*)&
         " ------------------------------------- Spectrum ------------------------------------- "
    counter=1
    do i_count=1,num_print
       do i=1,ssym%n_irrep
          do m=1,ssym%dim(i)
             do is=1,symmetry_data_n_spin()
                if (num_list(i)%m(m,is).eq.i_count) then
                   if(eigval(i)%m(m,is)>eig_min_loc.and.&
                        eigval(i)%m(m,is)<eig_max_loc) then
                      if(pop_col_loc) &
                           write(output_unit,1000) &
                           counter,m,ssym%name(i),is, &
                           occ_num(i)%m(m,is),&
                           eigval(i)%m(m,is)*convert1,pop_col(i)%m(m,is,:)
                      if(pop_store_loc) then
                         write(output_unit,*)
                         write(output_unit,*)
                         write(output_unit,1000) counter,m,ssym%name(i),is, &
                           occ_num(i)%m(m,is),eigval(i)%m(m,is)*convert1
                         do i_frag = 1, n_fragments
                            pointer_pop=>pop_store(i,is,i_frag)%eig_cont(counter_ir(i))%pop
                            pointer_coeff=>pop_store(i,is,i_frag)%eig_cont(counter_ir(i))%coeff
                            pointer_index=>pop_store(i,is,i_frag)%eig_cont(counter_ir(i))%index
                            n_cont=pop_store(i,is,i_frag)%eig_cont(counter_ir(i))%n_cont
                            write(output_unit,'(A12,I5,A10,F10.7)') 'Fragment:',i_frag,'  Sum:',&
                                 pop_store(i,is,i_frag)%eig_cont(counter_ir(i))%sum
                            write(output_unit,'(5(I4,A,F7.4,A,F7.4))') ( &
                                 pointer_index(i_print), '/', &
                                 pointer_pop(i_print), '/', &
                                 pointer_coeff(i_print), &
                                 i_print=1,n_cont)
                         end do! loop over fragments
                         counter_ir(i)=counter_ir(i)+1
                      end if
                   end if
                   counter = counter+1
                endif
             enddo
          enddo
       enddo
    enddo
    write(output_unit,*)" "
    write(output_unit,*)" ------------------------------------------------------------------------------------ "
    write(output_unit,*)" "
    if(pop_store_loc) then
       deallocate(counter_ir, stat=alloc_stat)
       if(alloc_stat/=0) &
            call error_handler('occupation_print_popspectrum: deallocation counter_ir failed')
    end if
1000 format(I4," ",I3," ",A4," ",I3,F7.4," ",F13.6,20F8.5)

    return
  end subroutine occupation_print_popspectrum
  !*************************************************************


  !*************************************************************
  subroutine eigenstates_store(n_vir,th,mode)
    ! Purpose: stores the occupied orbitals (and their occupations)
    !          and the first few virtual ones on a readwriteblocked
    !          file in case the 'save_eigenvec' option is set.
    !
    ! << th >>     << mode >>   action
    ! PRESENT      NOT PRESENT  store present eigen states immediately
    ! NOT PRESENT  PRESENT      keep present eigen states for later storage
    ! PRESENT      PRESENT      store previously saved eigen states
    !
    ! routine called by: main_scf
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use readwriteblocked_module
    use options_module, only: options_perturbation_theory, options_spin_orbit
    use output_module     , only: output_main_scf, output_data_saved
    use eigen_data_module, only: eigval,eigvec
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), optional, intent(in) :: n_vir
    type(readwriteblocked_tapehandle), optional, intent(inout) :: th
    integer(kind=i4_kind), optional, intent(in) :: mode
    !------------ Declaration of local variables  ----------------
    integer(kind=i4_kind) :: i,s,n,k, n_spin,n_irrep,n_dim,n_orb,n_mod
    integer(kind=i4_kind), allocatable :: n_save(:,:)
    logical  :: reset
    integer  :: status
    external error_handler
    !------------ Executable code --------------------------------

    if ( options_spin_orbit ) then
      ABORT('SAVE_EIGENVEC not yet with SO')
    endif

    n_spin = symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()

    if (.not.present(th)) then
       ! keep modus
       ! take care: perturbation_theory flag may be changed during a SCF run
       ! only the master has the most recent value
       if (options_perturbation_theory()) then
          if (.not.eigen_kept) then
             allocate(n_rot(n_irrep,n_spin),stat=status)
             if (status /= 0) call error_handler &
                  ("EIGENSTATES_STORE: allocation of n_rot failed")
             allocate(eigvec_kept(n_irrep),stat=status)
             if (status /= 0) call error_handler &
                  ("EIGENSTATES_STORE: allocation of eigvec_kept failed")
             do i=1,n_irrep
                n_dim = symmetry_data_dimension(i)
                allocate(eigvec_kept(i)%m(n_dim,2,n_spin),stat=status)
                if (status /= 0) call error_handler &
                     ("EIGENSTATES_STORE: allocation of eigvec_kept(i)%m failed")
             enddo
             eigen_kept = .true.
          endif
          keep_eigen = .true. ! to envoke saving in subroutine chargefit
       endif
       return
    endif

    if (output_main_scf) call write_to_output_units &
         ("EIGENSTATES_STORE: saving eigenstates")

    if (present(mode)) then
       reset = options_perturbation_theory()
    else
       reset = .false.
    endif

    allocate(n_save(n_spin,n_irrep),stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_STORE : allocation of n_save failed")

    call readwriteblocked_write((/real(n_spin,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored eigenstates :'
       write(output_unit,'(  a     )')'n_spin'
       write(output_unit,'(4es20.13)')(/real(n_spin,r8_kind)/)
    endif
    do s=1,n_spin
       do i=1,n_irrep
          n_save(s,i) = symmetry_data_dimension(i)
       enddo
       if (present(n_vir)) then
          n_save(s,:) = min(n_save(s,:),n_occo(s,:)+n_vir)
       endif
       call readwriteblocked_write(real(n_occo(s,:),r8_kind),th)
       call readwriteblocked_write(real(n_save(s,:),r8_kind),th)
       if (output_data_saved) then
          write(output_unit,'( a,i1,a )')'n_occo(',s,',:)'
          write(output_unit,'(4es20.13)')real(n_occo(s,:),r8_kind)
          write(output_unit,'( a,i1,a )')'n_save(',s,',:)'
          write(output_unit,'(4es20.13)')real(n_save(s,:),r8_kind)
       endif
    enddo

    n_mod = -1
    do i=1,n_irrep
       n_dim = symmetry_data_dimension(i)
       do s=1,n_spin
          n_orb = n_save(s,i)
          if (n_orb > 0) then
             call readwriteblocked_write(eigval(i)%m(1:n_orb,s),th)
             call readwriteblocked_write(occ_num(i)%m(1:n_orb,s),th)
             if (output_data_saved) then
                1000 format(a,'(',i2,')%m(1:',i4,',',i1,')')
                write(output_unit, fmt = 1000 )'eigval',i,n_orb,s
                write(output_unit,'(4es20.13)')eigval(i)%m(1:n_orb,s)
                write(output_unit, fmt = 1000 )'occ_num',i,n_orb,s
                write(output_unit,'(4es20.13)')occ_num(i)%m(1:n_orb,s)
             endif
             if (reset) n_mod = n_rot(i,s)
             do n=1,n_orb
                k = n - n_mod + 1
                if (k==1 .or. k==2) then
                   call readwriteblocked_write(eigvec_kept(i)%m(:,k,s),th)
                   if (output_data_saved) then
                      write(output_unit, fmt = 2000 )i,n,s
                      write(output_unit,'(4es20.13)')eigvec_kept(i)%m(:,k,s)
                   endif
                else
                   call readwriteblocked_write(eigvec(i)%m(:,n,s),th)
                   if (output_data_saved) then
                      2000 format('eigvec(',i2,')%m(:,',i4,',',i1,')')
                      write(output_unit, fmt = 2000 )i,n,s
                      write(output_unit,'(4es20.13)')eigvec(i)%m(:,n,s)
                   endif
                endif
             enddo
          endif
       enddo
       if (reset) then
          deallocate(eigvec_kept(i)%m,stat=status)
          if (status /= 0) call error_handler &
               ("EIGENSTATES_STORE: deallocation of eigvec_kept(i)%m failed")
       endif
    enddo
    if (reset) then
       deallocate(eigvec_kept,stat=status)
       if (status /= 0) call error_handler &
            ("EIGENSTATES_STORE: deallocation of eigvec_kept failed")
       deallocate(n_rot,stat=status)
       if (status /= 0) call error_handler &
            ("EIGENSTATES_STORE: deallocation of n_rot failed")
       eigen_kept = .false.
    endif

    deallocate(n_save,stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_STORE : deallocation of n_save failed")

  end subroutine eigenstates_store

!*************************************************************
! record  1: n_spin
! record  2: n_occo(1,:)
! record  3: n_save(1,:)
! record  4: n_occo(2,:)                      [if n_spin > 1]
! record  5: n_save(2,:)                      [if n_spin > 1]
! FOREACH irrep i DO
! record  6: eigval (i)%m(  1:n_save(1,i),1)  [if n_save(1,i) > 0]
! record  7: occ_num(i)%m(  1:n_save(1,i),1)  [if n_save(1,i) > 0]
! record  8: eigvec (i)%m(:,1:n_save(1,i),1)  [if n_save(1,i) > 0]
! record  9: eigval (i)%m(  1:n_save(2,i),2)  [if n_spin > 1 & n_save(2,i)]
! record 10: occ_num(i)%m(  1:n_save(2,i),2)  [if n_spin > 1 & n_save(2,i)]
! record 11: eigvec (i)%m(:,1:n_save(2,i),2)  [if n_spin > 1 & n_save(2,i)]
! DONE
!*************************************************************

  subroutine eigenstates_recover(n_vir,th)
    ! Purpose: recovers the occupied orbitals (and their occupations)
    !          and the first few virtual ones from a readwriteblocked
    !          file in case the 'save_eigenvec' option is set.
    !
    ! routine called by: main_scf
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use readwriteblocked_module
    use output_module     , only: output_main_scf, output_data_read
    use eigen_data_module , only: eigval,eigvec
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), optional, intent(in) :: n_vir
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !------------ Declaration of local variables  ----------------
    allocatable           :: n_save, n_need, buffer
    integer(kind=i4_kind) :: i,s,n, n_spin, n_irrep, n_dim, n_orb, &
                             spin_stored, n_save(:,:), n_need(:,:)
    real(kind=r8_kind)    :: spin(1), buffer(:), half = 0.5_r8_kind
    integer               :: status
    !------------ Executable code --------------------------------
    external error_handler

    if (output_main_scf) call write_to_output_units &
         ("EIGENSTATES_RECOVER: reading eigenstates")
    if (.not.allocated(n_occo)) then
       call alloc_n_occo(ssym)
    endif
    if (.not.allocated(occ_num)) then
       call alloc_occ_num()
    endif

    n_spin = symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()
    call readwriteblocked_read(spin,th)
    spin_stored = int(spin(1),i4_kind)
    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered eigenstates :'
       write(output_unit,'(  a    )')'n_spin'
       write(output_unit,'(4es20.13)')spin(1)
    endif
    if (spin_stored > n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("EIGENSTATES_RECOVER: trying to convert spin-polarized data from")
       call write_to_output_units &
            ("                     tape into the spin-restricted data required")
       ! n_occo(i)        := max(s) n_occo(s,i) =: n_occo(smax(i),i)
       ! occ_num(i)%m(:)  := sum(s) occ_num(i)%m(:,s)
       ! eigval(i)%m(:)   := eigval(i)%m(:,smax(i))
       ! eigvec(i)%m(:,:) := eigval(i)%m(:,:,smax(i))
    endif
    if (spin_stored < n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("EIGENSTATES_RECOVER: trying to convert spin-restricted data from")
       call write_to_output_units &
            ("                     tape into the spin-polarized data required")
       ! n_occo(:,s)        := n_occo(:)
       ! occ_num(i)%m(:,s)  := occ_num(i)%m(:) / 2
       ! eigval(i)%m(:,s)   := eigval(i)%m(:)
       ! eigvec(i)%m(:,:,s) := eigvec(i)%m(:,:)
    endif

    allocate(n_save(spin_stored,n_irrep),buffer(n_irrep),stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_RECOVER : allocation (1) failed")
    do s=1,spin_stored
       call readwriteblocked_read(buffer,th)
       if (output_data_read) then
          write(output_unit,'( a,i1,a )')'n_occo(',s,',:)'
          write(output_unit,'(4es20.13)')buffer
       endif
       if (s <= n_spin) then
          n_occo(s,:) = int(buffer,i4_kind)
       else ! s = spin_stored > n_spin = 1
          n_occo(1,:) = max(n_occo(1,:),int(buffer,i4_kind))
       endif
       call readwriteblocked_read(buffer,th)
       if (output_data_read) then
          write(output_unit,'( a,i1,a )')'n_save(',s,',:)'
          write(output_unit,'(4es20.13)')buffer
       endif
       n_save(s,:) = int(buffer,i4_kind)
    enddo
    if (n_spin > 1 .and. spin_stored == 1) then
       n_occo(2,:) = n_occo(1,:)
    endif

    allocate(n_need(n_spin,n_irrep),stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_RECOVER : allocation (1') failed")
    do s=1,n_spin
       do i=1,n_irrep
          n_need(s,i) = symmetry_data_dimension(i)
       enddo
       if (present(n_vir)) n_need(s,:) = min(n_need(s,:),n_occo(s,:)+n_vir)

       if (spin_stored > 1 .and. n_spin == 1) then ! s = 1
          if (any( n_need(1,:) > max(n_save(1,:),n_save(2,:)) )) &
               call error_handler &
               ("EIGENSTATES_RECOVER : insufficient virtual orbitals found")
       elseif (s <= spin_stored) then
          if (any( n_need(s,:) > n_save(s,:) )) call error_handler &
               ("EIGENSTATES_RECOVER : insufficient virtual orbitals found")
       else ! s = n_spin > spin_stored = 1
          if (any( n_need(s,:) > n_save(1,:) )) call error_handler &
               ("EIGENSTATES_RECOVER : insufficient virtual orbitals found")
       endif
    enddo

    deallocate(n_need,buffer,stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_RECOVER : deallocation (1) failed")

    do i=1,n_irrep
       n_dim = symmetry_data_dimension(i)
       do s=1,spin_stored
          n_orb = n_save(s,i)
          if (n_orb > 0) then
             if (s <= n_spin) then
                call readwriteblocked_read(eigval(i)%m(1:n_orb,s),th)
                call readwriteblocked_read(occ_num(i)%m(1:n_orb,s),th)
                if (output_data_read) then
                    1000 format(a,'(',i2,')%m(1:',i4,',',i1,')')
                    1001 format(a,'(',i2,')%m(1:',i4,',',i1,') skipped')
                    write(output_unit, fmt = 1000 )'eigval',i,n_orb,s
                    write(output_unit,'(4es20.13)')eigval(i)%m(1:n_orb,s)
                    write(output_unit, fmt = 1000 )'occ_num',i,n_orb,s
                    write(output_unit,'(4es20.13)')occ_num(i)%m(1:n_orb,s)
                endif
                do n=1,n_orb
                   call readwriteblocked_read(eigvec(i)%m(:,n,s),th)
                   if (output_data_read) then
                      2000 format('eigvec(',i2,')%m(:,',i4,',',i1,')')
                      2001 format('eigvec(',i2,')%m(:,:,',i1,') skipped')
                      write(output_unit, fmt = 2000 )i,n,s
                      write(output_unit,'(4es20.13)')eigvec(i)%m(:,n,s)
                   endif
                enddo
             else ! s = spin_stored > n_spin = 1
                if (n_orb > n_save(1,i)) then
                   call readwriteblocked_read(eigval(i)%m(1:n_orb,1),th)
                   if (output_data_read) then
                      write(output_unit, fmt = 1000 )'eigval',i,n_orb,2
                      write(output_unit,'(4es20.13)')eigval(i)%m(1:n_orb,1)
                   endif
                else
                   call readwriteblocked_skipread(n_orb,th)
                   if (output_data_read) then
                      write(output_unit, fmt = 1001 )'eigval',i,n_orb,2
                   endif
                endif
                allocate(buffer(n_orb),stat=status)
                if(status.ne.0 ) call error_handler &
                     ("EIGENSTATES_RECOVER : allocation (2) failed")
                call readwriteblocked_read(buffer,th)
                if (output_data_read) then
                   write(output_unit, fmt = 1000 )'occ_num',i,n_orb,2
                   write(output_unit,'(4es20.13)')buffer
                endif
                occ_num(i)%m(1:n_orb,1) = occ_num(i)%m(1:n_orb,1) + buffer
                deallocate(buffer,stat=status)
                if(status.ne.0 ) call error_handler &
                     ("EIGENSTATES_RECOVER : deallocation (2) failed")
                if (n_orb > n_save(1,i)) then
                   do n=1,n_orb
                      call readwriteblocked_read(eigvec(i)%m(:,n,1),th)
                      if (output_data_read) then
                         write(output_unit, fmt = 2000 )i,n,2
                         write(output_unit,'(4es20.13)')eigvec(i)%m(:,n,1)
                      endif
                   enddo
                else
                   call readwriteblocked_skipread(n_dim*n_orb,th)
                   if (output_data_read) then
                      write(output_unit, fmt = 2001 )i,2
                   endif
                endif
             endif
          endif
       enddo
       if (n_spin > 1 .and. spin_stored == 1) then
          n_orb = n_save(1,i)
          occ_num(i)%m(  1:n_orb,1) = occ_num(i)%m(  1:n_orb,1) * half
          eigval (i)%m(  1:n_orb,2) = eigval (i)%m(  1:n_orb,1)
          occ_num(i)%m(  1:n_orb,2) = occ_num(i)%m(  1:n_orb,1)
          eigvec (i)%m(:,1:n_orb,2) = eigvec (i)%m(:,1:n_orb,1)
       endif
    enddo

    deallocate(n_save,stat=status)
    if(status.ne.0 ) call error_handler &
         ("EIGENSTATES_RECOVER : deallocation (1') failed")

  end subroutine eigenstates_recover
  !*************************************************************

  !*************************************************************
  subroutine occupation_symmetry_check()
    ! Purpose: check the 'n_nonempty_irreps'-specification in
    ! the namelist OCCUPATION with the variable 'symmetry_data_n_irreps'
    ! after the symmetry part has been run.
    !
    ! routine called by: main_master
    !--------------------------------------------------------
    use options_module, only: options_spin_restricted, options_spin_orbit
#ifdef WITH_MOLMECH
    use operations_module, only: operations_mol_mech
#endif
#ifdef WITH_EPE
    use operations_module, only: operations_epe_lattice
#endif
     !** End of interface *****************************************
    !------------ Declaration of local variables  ----------------
    integer(kind=i4_kind) :: n_spin,degen,i,num_spin,is,&
         alloc_stat,counter,i_hole,hole_count
    integer(kind=i4_kind) :: n_irreps
    integer(kind=i4_kind),allocatable :: dim_irrep(:),n_partners(:)
    !------------ Executable code --------------------------------

#ifdef WITH_MOLMECH
    if(operations_mol_mech) return
#endif
#ifdef WITH_EPE
    if(operations_epe_lattice) return
#endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !

       n_spin = 1_i4_kind
       num_spin = 1_i4_kind

       n_irreps = symmetry_data_n_proj_irreps()
       allocate(dim_irrep(n_irreps),n_partners(n_irreps),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("occupation_symmetry_check: allocation of dim_irrep and n_partners failed")
       dim_irrep = ssym%dim_proj
       n_partners = ssym%partner_proj
    else

       if (options_spin_restricted()) then
          n_spin = 1_i4_kind
       else
          n_spin = 2_i4_kind
       endif

       num_spin = 3_i4_kind-n_spin
       n_irreps = symmetry_data_n_irreps()
       allocate(dim_irrep(n_irreps),n_partners(n_irreps),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("occupation_symmetry_check: allocation of dim_irrep and n_partners failed")
       dim_irrep = ssym%dim
       n_partners = ssym%partner
    endif

    if (fixed) then
       if (n_nonempty_irreps/=n_irreps) then
          call error_handler &
               ("occupation_symmetry_check: n_nonempty_irreps is WRONG")
       endif
       do i=1,n_nonempty_irreps
          degen = n_partners(i)*num_spin
          do is=1,n_spin
             !f (noc_start(is,i)>=dim_irrep(i)*real(degen,kind=r8_kind)) then
             if (noc_start(is,i) > dim_irrep(i)*real(degen,kind=r8_kind)) then
                call error_handler&
                     ("occupation_symmetry_check: a fixed occupation number is too large")
             endif
          enddo
       enddo
    endif! if (fixed)

    if (fixed_hole) then
       ! check if the user-specified irreps where hole are located
       ! make sense in the current point group
       if (maxval(hole_irrep)>n_irreps) then
          call error_handler&
               ("occupation_symmetry_check: an irrep given in the list of holes is wrong")
       endif
       do i=1,n_holes
          degen=num_spin*n_partners(hole_irrep(i))
          if (hole_occnum(i)>real(degen,r8_kind)) then
             call error_handler&
                  ("occupation_symmetry_check: an occupation number for a hole is too large")
          endif
       enddo
       ! now initialize the variable 'occnum_hole'
       allocate(occnum_hole(n_irreps),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("occupation_symmetry_check: allocation (1) failed")
       do i=1,n_irreps
          counter=0
          do i_hole=1,n_holes
             if (hole_irrep(i_hole)==i) then
                counter=counter+1
             endif
          enddo
          hole_count=counter ! hole_count is now the number of hole is irrep i
          allocate(occnum_hole(i)%m(max(hole_count,1)),STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ("occupation_symmetry_check: allocation (2) failed")
          occnum_hole(i)%m = -99.0_r8_kind
          counter=0
          do i_hole=1,n_holes
             if (hole_irrep(i_hole)==i) then
                counter=counter+1
                if (counter>hole_count ) call error_handler &
                     ("occupation_symmetry_check: sth. went wrong when initializing holes")
                occnum_hole(i)%m(counter)=hole_occnum(i_hole)
             endif
          enddo
          !? print*,'List of occnums for holes in Irrep ',i,'     ',occnum_hole(i)%m
       enddo

    endif! if (fixed_hole)
    deallocate(dim_irrep,n_partners,STAT=alloc_stat)
    if(alloc_stat/=0) call error_handler&
         ("occupation_symmetry_check: deallocation of dim_irrep and n_partners failed")

  end subroutine occupation_symmetry_check
  !*************************************************************

  !*************************************************************
  function occupation_fixed()
    ! Purpose: gives acces to private variable 'fixed'
    !** End of interface *****************************************
    logical :: occupation_fixed
    if (fixed) then
       occupation_fixed=.true.
    else
       occupation_fixed=.false.
    endif
  end function occupation_fixed
  !*************************************************************

  !*************************************************************
  function occupation_fixed_hole()
    ! Purpose: gives acces to private variable 'fixed_hole'
    !** End of interface *****************************************
    logical :: occupation_fixed_hole
    if (fixed_hole) then
       occupation_fixed_hole=.true.
    else
       occupation_fixed_hole=.false.
    endif
  end function occupation_fixed_hole
  !*************************************************************

  !*************************************************************
  subroutine occupation_get_holes()
    !
    ! This   routine   is  actually   a   wrapper   for  the   routine
    ! eigen_hole_setup()    which    needs    variables    from    the
    ! occupation_module   but   which  cannot   import   them  via   a
    ! use-statement since  the occupation_module is  already using the
    ! eigen_data_module.   Of course,  since main_scf()  is  using the
    ! whole  occupation  module  with  all its  public  variables  the
    ! subsequent  call to  the routine  eigen_hole_setup()  could also
    ! have  been done  in main_scf(),  but would  have appeared  a bit
    ! awkward.
    !
    ! Subroutine called by: 'main_scf()'.
    !
    use options_module, only: options_spin_restricted
    use eigen_data_module, only: eigen_hole_setup
    use comm, only: comm_rank
    !** End of interface *****************************************

    logical :: force_hole

    ! The  staff below seems  to run  on master  only and  requires no
    ! communication but some IO:
    if (comm_rank() /= 0) return

    ! In other words, force_hole == .not. hole_localization
    if (hole_localization) then
       force_hole = .false.
    else
       force_hole = .true.
    endif

    if (options_spin_restricted()) then
       call eigen_hole_setup (force_hole, n_holes, hole_list, hole_irrep, &
            hole_update)
    else
       call eigen_hole_setup (force_hole, n_holes, hole_list, hole_irrep, &
            hole_update, hole_spin)
    endif
  end subroutine occupation_get_holes
  !*************************************************************

  !*************************************************************
  real(kind=r8_kind) function  occupation_jz()
    ! Purpose: determines total Jz (only makes sense for axial
    !          double groups
    !** End of interface *****************************************
    real(kind=r8_kind)       :: jz
    integer(kind=i4_kind)    :: i,m

    jz = 0
    do i=1,ssym%n_proj_irrep
       do m=1,ssym%dim_proj(i)
          jz = jz + occ_num(i)%m(m,1)*symmetry_data_jz(i)
       enddo
    enddo
    occupation_jz = jz
  end function occupation_jz
  !*************************************************************

  !--------------- End of module ----------------------------------
end module occupation_module

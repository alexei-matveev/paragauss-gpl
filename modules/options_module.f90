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
module options_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: Contains various steering variables to switch between
  !           different branches of program run.
  !
  !           Currentley the following information is contained:
  !              + spin restricted calculation ?
  !                (Note that this information was previously in
  !                symmetry_data_module and for sake of compability
  !                is copied there)
  !              + calculate exchange correlation numerically
  !                or via fitfunctions
  !              + Save SCF data for recovering at various stages:
  !                  save_fitcoeff: charge and XC fit coefficients as well
  !                                 as current convergence and mixing state
  !                  save_scfstate: as save_fitcoeff plus XC potential
  !                                 [the former save_xcks_matrix option]
  !                  save_ksmatrix: as save_fitcoeff plus Kohn-Sham matrix
  !                  save_eigenvec: as save_fitcoeff plus eigenstates
  !                                 and occupation numbers
  !              + Set up a continued SCF run using previously saved data
  !                  read_fitcoeff: charge and XC fit coefficients as well
  !                                 as actual convergence and mixing state
  !                  read_scfstate: as read_fitcoeff plus XC potential
  !                                 [the former read_xcks_matrix option]
  !                  read_ksmatrix: as read_fitcoeff plus Kohn-Sham matrix
  !                  read_eigenvec: as read_fitcoeff plus eigenstates
  !                                 and occupation numbers
  !              + Specify an SCF data save interval for regular
  !                saving during the calculation.
  !                [the former save_xcks_nbr options]
  !              + Use an approximated energy functional depending
  !                on a particular "model density" and the corres-
  !                ponding variational Hamiltonian instead of the
  !                standard Kohn-Sham procedure.
  !
  !           The values are taken from namelist main_options
  !           of the input file.
  !
  !           All informations are available on master and slaves.
  !
  !  Author: FN / TB
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   9/96
  ! Description: restructured, moved maximal number of itterations
  !              to convergence module, added spin_restricted
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   8/97
  ! Description: introduction of the "USE_MODEL_DENSITY" option
  !              introduction of the "PERTURBATION_THEORY" option
  !              introduction of the "DIRECT_ENERGY_CALC" option
  !              introduction of the "SAVE/READ_<SCF_DATA>" options
  !              introduction of the "SPLIT_GRADIENTS" option
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   8/97
  ! Description: default controlled formatted namelist output introduced
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: introduction of the "USE_EXTENDED_MDA" option
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description ...
  !
  ! Modification
  ! Author: am
  ! Date:   24.04.99
  ! Description: a key for debugging only:
  !                 function options_debug_key()
  !              and
  !                 (private?) variable debug_key

  !Modification (Please copy before editing)
  ! Author: DG
  ! Date: 30/11/2000
  ! Description: added options_kinematic_factors for g-tensor debugging:
  ! if (options_kinematic_factors == .false.) then don`t calculate kinematic factors
  ! default value is FALSE
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date: 03/02
  ! Description: save_densmatrix option was added
  !
  ! Modification: option_finite_nucleus is added
  ! Author: DG
  ! Date:   05/2003
  ! Description: option_finite_nucleus is added for finite nucleus approach
  !
  !
  ! Modification: perturbation_theory could be switched
  ! Author: AN
  ! Date:   06/2009
  ! Description: perturbation_theory could be switched outside this module
  !              ATTENTION: works only for the master
  !              the value of perturbation theory for the slaves
  !              has to be given via commpack from now on
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  use operations_module
  use symmetry_data_module
  use iounitadmin_module
  implicit none
  private
  save
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------

  public ::  option!(op) == true/false
  public :: ioption!(op) == integer value

  public :: options_spin_restricted, options_n_spin, &
       options_integrals_on_file, &
       options_read_input, options_read_scfcontrol, options_write_input, &
       options_save_ksmatrix, options_save_eigenvec, options_save_scfstate, &
       options_save_fitcoeff, options_save_interval, options_recover, &
       options_reset_scfcycle, options_xcmode, &
       options_perturbation_theory, options_split_gradients, &
       options_integral_expmax, options_save_eigenvec_all, &
       options_save_as_fragment, options_set_defaults, &
       options_gxfile_epeformat, options_debug_key, options_save_densmatrix, &
       options_directaccess_integrals,quadrupels_reclength, &
       options_directaccess_gradients, &
       options_orbitals_in_memory, &
       options_orbitals_on_file,&
       options_switch_pertth, &
       options_finite_nucleus

  public :: options_data_bcast!()

  !---- public parameters for enumeration type OPTIONS_XCMODE -----
  integer(kind=i4_kind), public, parameter :: xcmode_numeric_exch  = 0, &
                                              xcmode_exchange_fit  = 1, &
                                              xcmode_model_density = 2, &
                                              xcmode_extended_mda  = 3

  !---- public parameters for enumeration type OPTIONS_RECOVER -----
  integer(kind=i4_kind), public, parameter :: recover_nothing  = -1, &
                                              recover_ksmatrix =  1, &
                                              recover_eigenvec =  2, &
                                              recover_scfstate =  3, &
                                              recover_fitcoeff =  4, &
                                              recover_fragment =  5

  !---- public control variables -----------------------------------
  logical, public :: options_relativistic, options_spin_orbit, &
                     options_kinematic_factors

  integer(kind=i4_kind), public :: &
       df_max_geo_iteration =  1, &  ! kept for compatibility
       max_geo_iteration             ! kept for compatibility
  integer(kind=i4_kind), public ::  update_hessian_iteration, &
                                    df_update_hessian_iteration = 0

  !================================================================
  ! End of public interface of module
  !================================================================

  integer, parameter,    PRIVATE :: wlen=16
  character(wlen),       PRIVATE :: relativistic_keywords = 'false'

  integer(kind=i4_kind), PRIVATE :: adkh_inf  = 0
  integer(kind=i4_kind), PRIVATE :: adkh_fpfw = 0
  integer(kind=i4_kind), PRIVATE :: adkh_zrel = 0
  ! if not zero, enables DKH transformations of primitive
  ! integrals in SHGI-code so that electrons inside of
  ! heavy atoms become relativistic, but interatomic
  ! forces are not affected.
  ! Not to be used with options_relativistic!

  !------------- Definition of default input values ---------------
  character(wlen) :: &
             df_relativistic        = "false"
  logical :: df_numeric_exch        = .true. , &
             df_spin_restricted     = .true. , &
             df_save_xcks_matrix    = .false., & ! kept for compatibility
             df_read_xcks_matrix    = .false., & ! kept for compatibility
             df_energies_always     = .false., & ! kept for compatibility
             df_use_model_density   = .false., &
             df_use_extended_mda    = .false., &
             df_direct_energy_calc  = .false., & ! kept for compatibility
             df_perturbation_theory = .true. , &
             df_split_gradients     = .false., &
             df_spin_orbit          = .false., &
             df_kinematic_factors   = .false., &
             df_finite_nucleus      = .false., &
             df_lvshift_mixing      = .false., &
             df_orbitals_on_file    = .false., &
             df_orbitals_in_memory  = .false.
#ifdef FPP_BIG_MEMORY
  logical :: df_integrals_on_file = .false.
#else
  logical :: df_integrals_on_file = .true.
#endif
  logical :: df_directaccess_integrals = .false.
  integer(kind=i4_kind) :: df_quadrupels_reclength = 0
  integer(kind=i4_kind) :: df_save_xcks_nbr = 10 ! kept for compatibility
  logical :: df_gxfile_epeformat  = .true.
  logical :: df_save_ksmatrix     = .false., &
             df_read_ksmatrix     = .false., &
             df_save_eigenvec     = .false., &
             df_save_eigenvec_all = .false., &
             df_save_as_fragment  = .false., &
             df_read_eigenvec     = .false., &
             df_save_scfstate     = .false., &
             df_read_scfstate     = .false., &
             df_save_fitcoeff     = .false., &
             df_read_fitcoeff     = .false., &
             df_reset_scfcycle    = .true., &
             df_save_densmatrix    = .false.
  integer(kind=i4_kind) :: df_save_interval = 0
  real(kind=r8_kind) :: df_integral_expmax = 50.0_r8_kind


  !------------- Declaration of private variables -----------------
  character(len=wlen) :: &
             relativistic
  logical :: read_densmat       , & ! kept for compatibility
             numeric_exch       , &
             spin_restricted    , &
             save_xcks_matrix   , & ! kept for compatibility
             read_xcks_matrix   , & ! kept for compatibility
             energies_always    , & ! kept for compatibility
             use_model_density  , &
             use_extended_mda   , &
             direct_energy_calc , & ! kept for compatibility
             perturbation_theory, &
             split_gradients    , &
             spin_orbit         , &
             kinematic_factors  , &
             finite_nucleus     , &
             integrals_on_file  , &
             directaccess_integrals

  integer(kind=i4_kind) ::   quadrupels_reclength
  logical, public :: lvshift_mixing

  logical :: gxfile_epeformat, &
             orbitals_on_file   , &
             orbitals_in_memory , &
             options_finite_nucleus


  integer(kind=i4_kind) :: save_xcks_nbr ! kept for compatibility
  logical :: save_ksmatrix     , &
             read_ksmatrix     , &
             save_eigenvec     , &
             save_eigenvec_all , &
             save_as_fragment  , &
             read_eigenvec     , &
             save_scfstate     , &
             read_scfstate     , &
             save_fitcoeff     , &
             read_fitcoeff     , &
             read_oldlcgto     , & ! kept for compatibility
             reset_scfcycle    , &
             save_densmatrix
  integer(kind=i4_kind) :: recover_mode , &
                           save_interval
  real(kind=r8_kind) :: integral_expmax

  integer(kind=i4_kind) :: debug_key=0

  namelist /main_options/ read_densmat, numeric_exch, split_gradients, &
       spin_restricted, save_xcks_matrix, read_xcks_matrix, save_xcks_nbr, &
       energies_always, max_geo_iteration, relativistic, spin_orbit, kinematic_factors,finite_nucleus, &
       use_model_density, use_extended_mda, direct_energy_calc, &
       perturbation_theory, integral_expmax, integrals_on_file, &
       directaccess_integrals,quadrupels_reclength,&
       orbitals_on_file, orbitals_in_memory, &
       gxfile_epeformat, debug_key,update_hessian_iteration
  namelist /recover_options/ &
       save_ksmatrix, read_ksmatrix, save_eigenvec, save_eigenvec_all, read_eigenvec, &
       save_scfstate, read_scfstate, save_fitcoeff, read_fitcoeff, &
       save_interval, read_oldlcgto, reset_scfcycle, save_as_fragment, save_densmatrix

  !------------ Subroutines ---------------------------------------
contains

  recursive function option(op) result(val)
    implicit none
    character(len=*), intent(in) :: op!tion
    logical                      :: val ! true if set
    ! *** end of interface ***

    ! fallback to integer options:
    val = ( ioption(op) /= 0 )
  end function option

  recursive function ioption(op) result(val)
    ! takes a string as an option name, retruns its value
    implicit none
    character(len=*), intent(in) :: op!tion
    integer(i4_kind)             :: val!ue
    ! *** end of interface ***

    select case(op)
    case ('adkh')
      val = adkh_zrel
    case ('adkh.inf')
      val = adkh_inf
    case ('adkh.fpfw')
      val = adkh_fpfw
    case default
      val = -1                  ! make compiler happy
      ABORT('ioption: '//op)
    end select
  end function ioption


  subroutine options_set_defaults
    ! purpose: prepares reading read the namelists main_options
    !          and recover_options from input file by setting
    !          variables to defaults.
    ! operations_geo_opt is used to decide if reset_scfcycle,
    ! save_scfstate and read_scfstate are set.
    !** End of interface *****************************************
    !------ Re-definition of input dependent default valuse ----
    !am: df_reset_scfcycle = operations_geo_opt .or. operations_qm_mm
    !am: df_save_scfstate  = operations_geo_opt .or. operations_qm_mm
    !am: df_read_scfstate  = operations_geo_opt .or. operations_qm_mm
    ! to facilitate switching between "QM-MM" and "SymAdaptGX" mode
    if ( operations_symadapt_gx .or. operations_potential) then  !!!!!!!!!!!!!!AS
       df_save_ksmatrix     = .false.
       df_save_eigenvec     = .false.
       df_save_eigenvec_all = .false.
       df_save_scfstate     = .false.
       df_save_fitcoeff     = .false.
       df_read_ksmatrix     = .false.
       df_read_eigenvec     = .false.
       df_read_scfstate     = .false.
       df_read_fitcoeff     = .false.
    endif

    numeric_exch        = df_numeric_exch
    spin_restricted     = df_spin_restricted
    save_xcks_matrix    = df_save_xcks_matrix  ! kept for compatibility
    read_xcks_matrix    = df_read_xcks_matrix  ! kept for compatibility
    save_xcks_nbr       = df_save_xcks_nbr     ! kept for compatibility
    energies_always     = df_energies_always   ! kept for compatibility
    max_geo_iteration   = df_max_geo_iteration ! kept for compatibility
    update_hessian_iteration   = df_update_hessian_iteration
    use_model_density   = df_use_model_density
    use_extended_mda    = df_use_extended_mda
    direct_energy_calc  = df_direct_energy_calc! kept for compatibility
    perturbation_theory = df_perturbation_theory
    split_gradients     = df_split_gradients
    relativistic        = df_relativistic
    spin_orbit          = df_spin_orbit
    kinematic_factors   = df_kinematic_factors
    finite_nucleus      = df_finite_nucleus
    integrals_on_file   = df_integrals_on_file
    directaccess_integrals  = df_directaccess_integrals
    quadrupels_reclength  = df_quadrupels_reclength
    gxfile_epeformat    = df_gxfile_epeformat

    orbitals_on_file    = df_orbitals_on_file
    orbitals_in_memory  = df_orbitals_in_memory

    integral_expmax     = df_integral_expmax

    save_ksmatrix      = df_save_ksmatrix
    read_ksmatrix      = df_read_ksmatrix
    save_eigenvec      = df_save_eigenvec
    save_eigenvec_all  = df_save_eigenvec_all
    save_as_fragment   = df_save_as_fragment
    read_eigenvec      = df_read_eigenvec
    save_scfstate      = df_save_scfstate
    read_scfstate      = df_read_scfstate
    save_fitcoeff      = df_save_fitcoeff
    read_fitcoeff      = df_read_fitcoeff
    save_interval      = df_save_interval
    reset_scfcycle     = df_reset_scfcycle
    save_fitcoeff  = df_save_fitcoeff
    read_fitcoeff  = df_read_fitcoeff
    save_scfstate  = df_save_scfstate
    read_scfstate  = df_read_scfstate
    save_interval  = df_save_interval
    reset_scfcycle = df_reset_scfcycle
    save_densmatrix = df_save_densmatrix
    lvshift_mixing = df_lvshift_mixing
  end subroutine options_set_defaults

  !*************************************************************

  subroutine options_read_input()
    ! purpose: read the namelists main_options
    !          or recover_options from input file
    ! routine called by: read_input
    !** End of interface *****************************************
    use input_module
    use operations_module, only: operations_fitbasis_opt, &
         namelist_tasks_used, operations_task
    use iounitadmin_module, only: write_to_output_units,write_to_trace_unit
    use strings, only: tolower, takeword
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: unit,status
    external error_handler
    logical                :: force_recover
    character(len=wlen)    :: word
    integer(i4_kind)       :: pos
    !------------ Executable code --------------------------------

    if ( input_line_is_namelist("MAIN_OPTIONS") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=main_options, iostat=status)
       if (status .gt. 0) call input_error( &
            "OPTIONS_READ_INPUT: Namelist main_options")
       if(namelist_tasks_used) call task_options(operations_task)

    if(directaccess_integrals.and.quadrupels_reclength.eq.0) &
                                     quadrupels_reclength=1024
    elseif ( input_line_is_namelist("RECOVER_OPTIONS") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=recover_options, iostat=status)
       if (status .gt. 0) call input_error( &
            "OPTIONS_READ_INPUT: Namelist recover_options")
    endif

    if (orbitals_in_memory) then
      print *, "ERROR: orbitals_in_memory: This option is broken"
      print *, "      dynamic distribution of the grid jobs does"
      print *, "      not work together with storing the orbitals"
      ABORT("orbitals_in_memory has to be set to false")
      !Attention: code will also abort if reaching the point where
      !integral storage is created
    endif


    force_recover = operations_geo_opt .or. operations_qm_mm

    if(operations_potential) then  !!!!!!!!!!!!!!AS
       spin_orbit=.false.
!       integrals_on_file=.true.
    endif

    ![[=== parse main_options.relativistic ===
    adkh_zrel     = 0
    adkh_inf      = 0
    adkh_fpfw     = 0
    relativistic_keywords = tolower(relativistic)

    ! parse the legacy true/false strings (defaults to false on error):
    options_relativistic = truefalse(relativistic,stat=status)
    if( status /= 0 )then
      ! neither true, nor false:
      pos = 1
      do while( pos < len(relativistic_keywords) )
        status = 0
        ! take the next token from the list of keywords:
        word = takeword(relativistic_keywords,pos)
        print *,'options_read: word=>'//word//'<'

        select case( word )
        case ('atomic1')
          adkh_zrel =  1 ! H
        case ('atomic2')
          adkh_zrel =  3 ! He+1 = Li
        case ('atomic3')
          adkh_zrel = 11 ! Ne+1 = Na
        case ('atomic4')
          adkh_zrel = 19 ! Ar+1 = K
        case ('atomic5','atomic','adkh')
          adkh_zrel = 37 ! Kr+1 = Rb
        case ('inf')
          adkh_inf  = 1
        case ('fpfw')
          adkh_fpfw = 1
        case ('')
          exit ! do while loop
        case default
          status = 1
          exit ! do while loop
        end select
      enddo
    endif
    if( status /= 0 )then
      ! the default behaviour on error:
      print*,'ERROR: Invalid value of option relativistic=',relativistic
      print*,'       Valid values are:'
      print*,'       TRUE, FALSE (and synonymous)'
      print*,'       ADKH, ADKH:Inf, ADKH:fpFW, ...'
      ABORT('no such case, see tty')
    endif
    if( adkh_zrel /= 0 )then
        WARN('ADKH enabled!')
    endif
    if( adkh_fpfw /= 0 )then
        WARN('ADKH on interatomic ints only till fpFW!')
    endif
    !]]=======================================

    if (split_gradients .and. options_relativistic) call input_error( &
         "Sorry, option 'split_gradients' not yet implemented &
         &for relativistic forces")

    options_spin_orbit=spin_orbit
    options_kinematic_factors = kinematic_factors
    options_finite_nucleus = finite_nucleus
    if(options_spin_orbit)then
       if(.not.integrals_on_file)then
          WARN('SO requires INTEGRALS_ON_FILE, forcing')
          integrals_on_file = .true.
       endif
       if(save_eigenvec)then
          WARN('SO not yet with SAVE_EIGENVEC, revert')
          save_eigenvec     = .false.
          save_eigenvec_all = .false.
       endif
       if(.not.spin_restricted)then
          WARN('SO not yet with SPIN_RESTRICTED=F, revert')
          spin_restricted = .true.
       endif
       ASSERT(options_relativistic)
       ASSERT(numeric_exch)
    endif

    call symmetry_data_set(n_spin=options_n_spin())

    if (save_xcks_matrix .and. .not. save_scfstate) then
       save_scfstate = .true.
       if (save_xcks_nbr > 99999) then
          save_interval = 0
       else
          save_interval = save_xcks_nbr
       endif
       save_xcks_matrix = .false.
    endif
    if (read_xcks_matrix .and. .not. read_scfstate) then
       read_scfstate = .true.
       read_xcks_matrix = .false.
    endif

    ! to facilitate switching between "QM-MM" and "SymAdaptGX" mode
!   if ( operations_symadapt_gx .or. .not. operations_geo_opt .or. &
!        operations_potential ) then
!   MF bug fix
    if ( operations_symadapt_gx .or. operations_potential ) then
       save_ksmatrix     = .false.
       save_ksmatrix     = .false.
!       save_eigenvec     = .false.
!       save_eigenvec_all = .false.
       save_scfstate     = .false.
       save_fitcoeff     = .false.
       read_ksmatrix     = .false.
!       read_eigenvec     = .false.
       read_scfstate     = .false.
       read_fitcoeff     = .false.
    endif

    ! check and set recover_mode
    recover_mode = recover_nothing
    status = 0
    if (read_fitcoeff) then
       recover_mode = recover_fitcoeff ; status = status + 1
    endif
    if (read_scfstate) then
       recover_mode = recover_scfstate ; status = status + 1
       if ( .not. save_scfstate .and. force_recover ) then
          WARN("GeoOpt: forcing save_scfstate")
          save_scfstate = .true.
       endif
    endif
    if (read_ksmatrix) then
       recover_mode = recover_ksmatrix ; status = status + 1
       if ( .not. save_ksmatrix .and. force_recover ) then
           WARN("GeoOpt: forcing save_ksmatrix")
           save_ksmatrix = .true.
       endif
    endif
    if (read_eigenvec) then
       recover_mode = recover_eigenvec ; status = status + 1
    endif
    if (operations_fitbasis_opt) then
       recover_mode = recover_fragment ; status = status + 1
    endif
    if ( status .eq. 0 .and. force_recover ) then
       WARN("GeoOpt: forcing read/save_scfstate")
       read_scfstate = .true.
       save_scfstate = .true.
       recover_mode = recover_scfstate
       status = 1
    else if ( status > 1 ) then
       call error_handler('OPTIONS_READ_INPUT: '// &
         'At most one of the "READ_..." options must be used!')
    endif

  end subroutine options_read_input

  !*************************************************************

  subroutine options_read_scfcontrol(warning)
    ! purpose: read the save options in the namelist recover_options
    !          from the input file in $TTFSSTART/scfcontrol at each
    !          scf-cycle, giving a notice if parameters were changed.
    use input_module
    implicit none
    !------------ Declaration of formal parameters -------------
    logical, intent(in), optional :: warning ! if present and .false.
                                             ! warning messages are
                                             ! suppressed
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    logical :: save_ksmatrix_old     , &
               read_ksmatrix_old     , &
               save_eigenvec_old     , &
               save_eigenvec_all_old , &
               save_as_fragment_old  , &
               read_eigenvec_old     , &
               save_scfstate_old     , &
               read_scfstate_old     , &
               save_fitcoeff_old     , &
               read_fitcoeff_old     , &
               reset_scfcycle_old
    integer(kind=i4_kind) :: &
         recover_mode_old ,&
         save_interval_old, &
         unit, &
         status
    !------------ Executable code ------------------------------
    save_ksmatrix_old  = save_ksmatrix
    read_ksmatrix_old  = read_ksmatrix
    save_eigenvec_old  = save_eigenvec
    save_eigenvec_all_old = save_eigenvec_all
    save_as_fragment_old   = save_as_fragment
    read_eigenvec_old  = read_eigenvec
    save_scfstate_old  = save_scfstate
    read_scfstate_old  = read_scfstate
    save_fitcoeff_old  = save_fitcoeff
    read_fitcoeff_old  = read_fitcoeff
    reset_scfcycle_old = reset_scfcycle
    recover_mode_old   = recover_mode
    save_interval_old  = save_interval

!   read(inp_unit, nml=recover_options, iostat=status)           !!!!!!!!!!!!!!AS
!   if (status .gt. 0) call input_error( &
!        "OPTIONS_READ_SCFCONTROL: Namelist recover_options(1)")
    if ( input_line_is_namelist("RECOVER_OPTIONS") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=recover_options, iostat=status)
       if (status .gt. 0) call input_error( &
            "OPTIONS_READ_SCFCONTROL: Namelist recover_options")
    endif

    if (save_xcks_matrix .and. .not. save_scfstate) then
       save_scfstate = .true.
       if (save_xcks_nbr > 99999) then
          save_interval = 0
       else
          save_interval = save_xcks_nbr
       endif
       save_xcks_matrix = .false.
    endif
    if (read_xcks_matrix .and. .not. read_scfstate) then
       read_scfstate = .true.
       read_xcks_matrix = .false.
    endif

    ! check and set recover_mode
    recover_mode = recover_nothing
    status = 0
    if (read_fitcoeff) then
       recover_mode = recover_fitcoeff ; status = status + 1
    endif
    if (read_scfstate) then
       recover_mode = recover_scfstate ; status = status + 1
    endif
    if (read_ksmatrix) then
       recover_mode = recover_ksmatrix ; status = status + 1
    endif
    if (read_eigenvec) then
       recover_mode = recover_eigenvec ; status = status + 1
    endif
    if (operations_fitbasis_opt) then
       recover_mode = recover_fragment ; status = status + 1
    endif
    if (status > 1) call error_handler('OPTIONS_READ_SCFCONTROL: '// &
         'At most one of the "READ_..." options must be used!')

    read_ksmatrix  = read_ksmatrix_old
    read_eigenvec  = read_eigenvec_old
    read_scfstate  = read_scfstate_old
    read_fitcoeff  = read_fitcoeff_old
    reset_scfcycle = reset_scfcycle_old
    recover_mode   = recover_mode_old

    if (present(warning)) then
       if (.not.warning) return
    endif

    if (save_ksmatrix_old .neqv. save_ksmatrix) then
       call write_to_output_units(line("save_ksmatrix",save_ksmatrix))
       call write_to_trace_unit(line("save_ksmatrix",save_ksmatrix))
    endif
    if (save_eigenvec_old .neqv. save_eigenvec) then
       call write_to_output_units(line("save_eigenvec",save_eigenvec))
       call write_to_trace_unit(line("save_eigenvec",save_eigenvec))
    endif
    if (save_eigenvec_all_old .neqv. save_eigenvec_all) then
       call write_to_output_units(line("save_eigenvec_all",save_eigenvec_all))
       call write_to_trace_unit(line("save_eigenvec_old",save_eigenvec_old))
    endif
    if (save_as_fragment_old .neqv. save_as_fragment) then
       call write_to_output_units(line("save_as_fragment",save_as_fragment))
       call write_to_trace_unit(line("save_as_fragment",save_as_fragment))
    endif
    if (save_scfstate_old .neqv. save_scfstate) then
       call write_to_output_units(line("save_scfstate",save_scfstate))
       call write_to_trace_unit(line("save_scfstate",save_scfstate))
    endif
    if (save_fitcoeff_old .neqv. save_fitcoeff) then
       call write_to_output_units(line("save_fitcoeff",save_fitcoeff))
       call write_to_trace_unit(line("save_ksmatrix",save_ksmatrix))
    endif
    if (save_interval_old .ne. save_interval) then
       call write_to_output_units("options_read_scfcontrol: &
            &save_interval was altered. new value: ",inte=save_interval)
       call write_to_trace_unit("options_read_scfcontrol: &
            &save_interval was altered. new value: ",inte=save_interval)
    endif

    contains

    character(len=68) function line(name,value)
       character(len=13), intent(in) :: name
       logical          , intent(in) :: value
       line(1:63) = "options_read_scfcontrol: "//name// &
                    " was altered. New value: "
       if (value) then
          line(64:68) = " TRUE"
       else
          line(64:68) = "FALSE"
       endif
       return
    end function line

  end subroutine options_read_scfcontrol

  !*************************************************************

  subroutine options_write_input(iounit,scfcontrol)
    !
    ! Purpose:  writes  the  main_options  namelist with  its  default
    ! values to the input.out file.
    !
    use echo_input_module, only: start, real, intg, flag, strng, stop, &
         echo_level_full
    implicit none
    integer, intent(in) :: iounit
    logical, intent(in), optional :: scfcontrol
    !** End of interface ***************************************

    !
    ! FIXME: use default format for integers/reals/flags ...
    !

    if (present(scfcontrol)) then
       if (scfcontrol) goto 1000
    endif

!   hidden input parameter:
!   save_xcks_matrix, read_xcks_matrix, save_xcks_nbr, &
!   max_geo_iteration

    call start("MAIN_OPTIONS","OPTIONS_WRITE_INPUT", &
         iounit,operations_echo_input_level)
    call flag("NUMERIC_EXCH       ",numeric_exch       ,df_numeric_exch       )
    call flag("SPIN_RESTRICTED    ",spin_restricted    ,df_spin_restricted    )
 !  call flag("ENERGIES_ALWAYS    ",energies_always    ,df_energies_always    )
    call strng("RELATIVISTIC       ",relativistic       ,df_relativistic       )
    call flag("SPIN_ORBIT         ",spin_orbit         ,df_spin_orbit         )
    call flag("KINEMATIC_FACTORS  ",kinematic_factors  ,df_kinematic_factors  )
    call flag("FINITE_NUCLEUS     ",finite_nucleus     ,df_finite_nucleus     )
    call flag("INTEGRALS_ON_FILE  ",integrals_on_file  ,df_integrals_on_file  )
    call flag("DIRECTACCESS_INTEGRALS ", &
                          directaccess_integrals  ,df_directaccess_integrals  )
    call intg("QUADRUPELS_RECLENGTH ", quadrupels_reclength &
                                                 ,df_quadrupels_reclength )
    call intg("UPDATE_HESSIAN_ITERATION ", update_hessian_iteration &
                                                 ,df_update_hessian_iteration )
!   call flag("ORBITALS_ON_FILE   ",orbitals_on_file   ,df_orbitals_on_file   )
    call flag("ORBITALS_ON_FILE   ",orbitals_on_file   ,df_orbitals_on_file   )
    call flag("ORBITALS_IN_MEMORY ",orbitals_in_memory ,df_orbitals_in_memory )
    call flag("USE_MODEL_DENSITY  ",use_model_density  ,df_use_model_density  )
    call flag("USE_EXTENDED_MDA   ",use_extended_mda   ,df_use_extended_mda   )
!   call flag("DIRECT_ENERGY_CALC ",direct_energy_calc ,df_direct_energy_calc )
    call flag("PERTURBATION_THEORY",perturbation_theory,df_perturbation_theory)
    call real("INTEGRAL_EXPMAX    ",integral_expmax    ,df_integral_expmax    )
    call flag("SPLIT_GRADIENTS    ",split_gradients    ,df_split_gradients    )
!!$    call flag("GXFILE_EPEFORMAT   ",gxfile_epeformat   ,df_gxfile_epeformat   )
    call stop()

    call start("RECOVER_OPTIONS","OPTIONS_WRITE_INPUT", &
         iounit,operations_echo_input_level)
    call flag("SAVE_KSMATRIX     ",save_ksmatrix     ,df_save_ksmatrix     )
    call flag("READ_KSMATRIX     ",read_ksmatrix     ,df_read_ksmatrix     )
    call flag("SAVE_EIGENVEC     ",save_eigenvec     ,df_save_eigenvec     )
    call flag("SAVE_EIGENVEC_ALL ",save_eigenvec_all ,df_save_eigenvec_all )
    call flag("SAVE_AS_FRAGMENT  ",save_as_fragment  ,df_save_as_fragment  )
    call flag("READ_EIGENVEC     ",read_eigenvec     ,df_read_eigenvec     )
    call flag("SAVE_SCFSTATE     ",save_scfstate     ,df_save_scfstate     )
    call flag("READ_SCFSTATE     ",read_scfstate     ,df_read_scfstate     )
    call flag("SAVE_FITCOEFF     ",save_fitcoeff     ,df_save_fitcoeff     )
    call flag("READ_FITCOEFF     ",read_fitcoeff     ,df_read_fitcoeff     )
    call intg("SAVE_INTERVAL     ",save_interval     ,df_save_interval     )
!   call flag("READ_OLDLCGTO     ",read_oldlcgto     ,df_read_oldlcgto     )
    call flag("RESET_SCFCYCLE    ",reset_scfcycle    ,df_reset_scfcycle    )
    call flag("SAVE_DENSMATRIX   ",save_densmatrix   ,df_save_densmatrix   )
    call stop()

    return

    1000 continue ! entry point for the "scfcontrol" mode

    call start("RECOVER_OPTIONS","OPTIONS_WRITE_INPUT",iounit,echo_level_full)
    call flag("SAVE_KSMATRIX ",save_ksmatrix ,df_save_ksmatrix )
    call flag("SAVE_EIGENVEC ",save_eigenvec ,df_save_eigenvec )
    call flag("SAVE_EIGENVEC_ALL ",save_eigenvec_all ,df_save_eigenvec_all )
    call flag("SAVE_AS_FRAGMENT ",save_as_fragment ,df_save_as_fragment )
    call flag("SAVE_SCFSTATE ",save_scfstate ,df_save_scfstate )
    call flag("SAVE_FITCOEFF ",save_fitcoeff ,df_save_fitcoeff )
    call intg("SAVE_INTERVAL ",save_interval ,df_save_interval )
    call stop(empty_line=.false.) !!!!!!!!!!!!!AS

  end subroutine options_write_input

  !*************************************************************

  logical function options_integrals_on_file()
    options_integrals_on_file = integrals_on_file
  end function options_integrals_on_file

  logical function options_directaccess_integrals()
    options_directaccess_integrals = directaccess_integrals &
      .or. quadrupels_reclength.ne.0
  end function options_directaccess_integrals

  logical function options_directaccess_gradients()
!   options_directaccess_gradients = .false.
   options_directaccess_gradients = directaccess_integrals &
        .or. quadrupels_reclength.ne.0
  end function options_directaccess_gradients


  logical function options_gxfile_epeformat()
    ! TRUE  : read gxfile with specification of epe centers
    ! false : read gxfile using regular formats
    options_gxfile_epeformat= gxfile_epeformat
  end function options_gxfile_epeformat

  logical function options_orbitals_in_memory()
    ! TRUE  : keep orbitals on grid in memory during SCF
    ! FALSE : calculate or read orbitals from file
    !** End of interface *******************************
    options_orbitals_in_memory = orbitals_in_memory
  end function options_orbitals_in_memory

  logical function options_orbitals_on_file()
    ! TRUE  : keep calculated orbitals in file
    ! FALSE : calculate or read orbitals from memory
    !** End of interface *******************************
    options_orbitals_on_file = orbitals_on_file
  end function options_orbitals_on_file

  !*************************************************************

  logical function options_spin_restricted()
    ! TRUE : closed shell calculation
    ! FALSE: open shell calculation
    !** End of interface *****************************************
    options_spin_restricted = spin_restricted
  end function options_spin_restricted

  !*************************************************************

  integer function options_n_spin()
    ! 1: closed shell calculation
    ! 2: open shell calculation
    !** End of interface *****************************************
    if (spin_restricted) then
       options_n_spin = 1
    else
       options_n_spin = 2
    endif
  end function options_n_spin

  logical function options_perturbation_theory()
    ! TRUE : use perturbation theory during charge fit to improve
    !        the SCF convergence
    ! FALSE: perform the charge fit without any corrections
    ! ATTENTION: perturbation_theory flag may be changed during a SCF run
    ! only the master has the most recent value
    ! so make sure to use this function only on the master
    !** End of interface *****************************************
    options_perturbation_theory = perturbation_theory .and. &
                                  .not.operations_fitbasis_opt
  end function options_perturbation_theory

  !*************************************************************

  subroutine options_switch_pertth(new_value)
    ! Purpose: as the perturbation_theory flag is a private variable it could
    ! normally only changed in this module
    ! this routine changes it from outside on the master
    ! as only the value of the master is changed slaves have to
    ! get a copie of it from the master from now on before they use
    ! the perturbation_theory
    logical, intent(in)      :: new_value
    !** End of interface *****************************************
    perturbation_theory = new_value

  end subroutine options_switch_pertth

  !*************************************************************

  logical function options_split_gradients()
    ! TRUE : if the analytic total energy gradients should be
    !        split into their individual contributions on output
    ! FALSE: just compute the complete total energy gradients
    !** End of interface *****************************************
    options_split_gradients = split_gradients
  end function options_split_gradients

  !*************************************************************

  integer(kind=i4_kind) function options_save_interval()
    ! n: if n > 0 SCF data is to be saved each n-th SCF cycle
    !** End of interface *****************************************
    options_save_interval = save_interval
  end function options_save_interval

  !*************************************************************

  logical function options_save_ksmatrix()
    ! TRUE : if the Kohn Sham matrix is to be saved at the end
    !        of the SCF procedure (or each n-th SCF cycle if
    !        n = options_save_interval() > 0
    ! FALSE: The Kohn Sham matrix is saved
    !** End of interface *****************************************
    options_save_ksmatrix = save_ksmatrix
  end function options_save_ksmatrix

  !*************************************************************

  logical function options_save_eigenvec()
    ! TRUE : if the eigenvectors, eigenvalues and occupation
    !        numbers of all occupied and the first virtual states
    !        are to be saved at the end of the SCF procedure (or
    !        each n-th SCF cycle if n = options_save_interval() > 0
    ! FALSE: The eigenstates are not saved
    !** End of interface *****************************************
    options_save_eigenvec = save_eigenvec
  end function options_save_eigenvec

  !*************************************************************

  logical function options_save_eigenvec_all()
    ! TRUE : save all eigenvectors and not only the occupied and
    !        n_vir virtual
    !** End of interface *****************************************
    options_save_eigenvec_all = save_eigenvec_all
  end function options_save_eigenvec_all

  !*************************************************************

  logical function options_save_as_fragment()
    ! TRUE : save only eigenvec and not also coeffch ...
    !** End of interface *****************************************
    options_save_as_fragment = save_as_fragment
  end function options_save_as_fragment

  !*************************************************************

  logical function options_save_scfstate()
    ! TRUE : if the fit coefficients (and the XC matrix for
    !        numeric exchange) are to be saved at the end of
    !        the SCF procedure (or each n-th SCF cycle if
    !        n = options_save_interval() > 0
    ! FALSE: The fit coefficients (and XC matrix) are not saved
    !** End of interface *****************************************
    options_save_scfstate = save_scfstate
  end function options_save_scfstate

  !*************************************************************

  logical function options_save_fitcoeff()
    ! TRUE : if the fit coefficients only are to be saved at the
    !        end of the SCF procedure (or each n-th SCF cycle if
    !        n = options_save_interval() > 0
    ! FALSE: The fit coefficients (and XC matrix) are not saved
    !** End of interface *****************************************
    options_save_fitcoeff = save_fitcoeff
  end function options_save_fitcoeff

  !*************************************************************

  logical function options_reset_scfcycle()
    ! TRUE : if the SCF cycle number is reset to one during
    !        recovering due to one of the "read_..." options
    ! FALSE: the SCF cycle number is recoverd as well
    !** End of interface *****************************************
    options_reset_scfcycle = reset_scfcycle
  end function options_reset_scfcycle

  !*************************************************************

  logical function options_save_densmatrix()
    !** End of interface *****************************************
    options_save_densmatrix = save_densmatrix
  end function options_save_densmatrix

  !*************************************************************

  integer(kind=i4_kind) function options_xcmode()
    ! Return the current xc mode
    !** End of interface *****************************************
    if ( use_extended_mda ) then
      options_xcmode = xcmode_extended_mda
    elseif ( use_model_density ) then
      options_xcmode = xcmode_model_density
    elseif ( numeric_exch ) then
      options_xcmode = xcmode_numeric_exch
    else
      options_xcmode = xcmode_exchange_fit
    endif
  end function options_xcmode

  !*************************************************************

  integer(kind=i4_kind) function options_recover()
    ! Return the current recover mode
    !** End of interface *****************************************
    options_recover = recover_mode
  end function options_recover
  !*************************************************************

  real(kind=r8_kind) function options_integral_expmax()
    ! Return cutoff criterion for integral part
    !** End of interface *****************************************
    options_integral_expmax = integral_expmax
  end function options_integral_expmax

  subroutine options_data_bcast()
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call options_data_pack()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call options_data_unpack()
    endif
  end subroutine options_data_bcast

  subroutine options_data_pack()
    !  Purpose: pack the structure options into a comm buffer
    !** End of interface *****************************************
    use comm_module
    use xpack, only: pck
    integer(kind=i4_kind)   :: info
    external error_handler
    !------------ Executable code --------------------------------
    call commpack(numeric_exch,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(spin_restricted,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(energies_always,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(use_model_density,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(use_extended_mda,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(direct_energy_calc,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(perturbation_theory,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(split_gradients,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(save_interval,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(save_ksmatrix,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(save_eigenvec,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(save_scfstate,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(save_fitcoeff,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(recover_mode,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(reset_scfcycle,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(integrals_on_file,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(directaccess_integrals,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(quadrupels_reclength,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing failed")
    call commpack(orbitals_in_memory,info)
    if(info /= 0) call error_handler &
         ("OPTIONS_DATA_PACK : packing orbitals_in_memory failed")
    call commpack(orbitals_on_file, info)
    if(info /= 0) call error_handler &
         ("OPTIONS_DATA_PACK : packing orbitals_on_file failed")
    call commpack(options_relativistic,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing options_relativistic failed")
    call commpack(options_spin_orbit,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing options_spin_orbit failed")
    call commpack(options_kinematic_factors,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing options_kinematic_factors failed")
    call commpack(options_finite_nucleus,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing options_finite_nucleus failed")
    call commpack(integral_expmax,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing integral_expmax failed")
    call commpack(debug_key,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_PACK : packing debug_key failed")
    call commpack(lvshift_mixing,info)

    call pck(relativistic)
    call pck(adkh_zrel)
    call pck(adkh_inf)
  end subroutine options_data_pack

  !*************************************************************

  subroutine options_data_unpack()
    ! purpose: unpack a comm buffer into the structure options
    !** End of interface *****************************************
    use comm_module
    use xpack, only: upck
    integer(kind=i4_kind)    :: info
    external error_handler
    !------------ Executable code --------------------------------
    call communpack(numeric_exch,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(spin_restricted,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(energies_always,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(use_model_density,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(use_extended_mda,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(direct_energy_calc,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(perturbation_theory,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(split_gradients,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(save_interval,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(save_ksmatrix,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(save_eigenvec,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(save_scfstate,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(save_fitcoeff,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(recover_mode,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(reset_scfcycle,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(integrals_on_file,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(directaccess_integrals,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")
    call communpack(quadrupels_reclength,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking failed")

    call communpack(orbitals_in_memory,info)
    if(info /= 0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking orbitals_in_memory failed")

    call communpack(orbitals_on_file, info)
    if(info /= 0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking orbitals_on_file failed")

    call communpack(options_relativistic,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking options_relativistic failed")
    call communpack(options_spin_orbit,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking options_spin_orbit failed")
    call communpack(options_kinematic_factors,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking options_kinematic_factors  failed")
    call communpack(options_finite_nucleus,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking options_finite_nucleus  failed")
    call communpack(integral_expmax,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking integral_expmax failed")
    call communpack(debug_key,info)
    if(info.ne.0) call error_handler &
         ("OPTIONS_DATA_UNPACK : unpacking debug_key failed")
    call communpack(lvshift_mixing,info)

    call upck(relativistic)
    call upck(adkh_zrel)
    call upck(adkh_inf)

    ! restore unsent input parameter
    spin_orbit    = options_spin_orbit
    kinematic_factors = options_kinematic_factors
    finite_nucleus = options_finite_nucleus
    read_ksmatrix = recover_mode == recover_ksmatrix
    read_eigenvec = recover_mode == recover_eigenvec
    read_scfstate = recover_mode == recover_scfstate
    read_fitcoeff = recover_mode == recover_fitcoeff

    ! reload input dependent default values
    df_reset_scfcycle = operations_geo_opt .or. operations_qm_mm
    !am: df_save_scfstate  = operations_geo_opt .or. operations_qm_mm
    !am: df_read_scfstate  = operations_geo_opt .or. operations_qm_mm

  end subroutine options_data_unpack


  function options_debug_key(mask) result(key)
    !
    ! FOR DEBUG PURPOSES
    !
    implicit none
    integer(kind=i4_kind),optional :: mask   ! bitmask
    integer(kind=i4_kind)          :: key    !<<<result
    ! *** end of interface ***

    key = IAND( debug_key, mask )
  end function options_debug_key


  recursive subroutine task_options(task)
    use input_module, only: input_error
    use strings, only: tolower
    implicit none
    character(len=11), intent(in)  :: task
    ! *** end of interface ***

    character(len=11) :: newtask
    character(len=11) :: subtask
    integer(i4_kind)  :: I


    newtask = tolower(task)
    I = INDEX( newtask,',')
    if( I /= 0 )then
       subtask = newtask(1:I-1)
       call task_options(subtask)
       subtask = newtask(I+1:11)
       call task_options(subtask)
       RETURN
    endif

    select case (newtask)
    case ('freq_num   ','numcarthess')
      update_hessian_iteration=0
      if(max_geo_iteration.le.2) max_geo_iteration=100
    case ('freq_analyt','freq_analit')
      update_hessian_iteration=1
      max_geo_iteration=1
      NUMERIC_EXCH = .true.
    case ('ts_search')
      if(update_hessian_iteration.eq.0) update_hessian_iteration=10
      if(max_geo_iteration.eq.0.or.max_geo_iteration.eq.1) max_geo_iteration=19
    case ('hess')
      if(update_hessian_iteration.eq.0) update_hessian_iteration=1
    case ('nohess','transit')
      update_hessian_iteration=0
    case ('geo','opt','geoopt','GeoOpt')
      if(update_hessian_iteration.eq.1) update_hessian_iteration=15
      if(max_geo_iteration.eq.0) max_geo_iteration=29
    end select

  end subroutine task_options

!--------------- End of module ----------------------------------
end module options_module

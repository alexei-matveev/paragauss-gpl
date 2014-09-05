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
!=====================================================================
! Public interface of module
!=====================================================================
module convergence_module
  !-------------------------------------------------------------------
  !
  !  Purpose: contains a collection of routines to check
  !           convergence on set of variables.
  !  Currently, energy criterion and max_iteration are checked
  !
  !  Module called by: main_scf
  !
  !
  !  References: old lcgto
  !
  !
  !  Author: FN
  !  Date: 1/96
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------
  public convergence, convergence_max_iterations, convergence_setup, &
         convergence_shutdown, convergence_read, convergence_write, &
         convergence_clear_buffers, &
         convergence_state_store, convergence_state_recover, &
         convergence_put_energy, convergence_put_coeff_dev, &
         convergence_put_coulomb_dev, convergence_put_density_dev, &
         convergence_check_energy_dev, convergence_check_coeff_dev, &
         convergence_check_coulomb_dev, convergence_check_density_dev, &
         convergence_load_coeff_dev, convergence_load_metric_dev, &
         convergence_max_geo_iteration, convergence_abort_calculation

#ifdef WITH_SCFCONTROL
  public :: convergence_read_scfcontrol
#endif /* ifdef WITH_SCFCONTROL */

  public :: convergence_resize_buffers!()


  !===================================================================
  ! End of public interface of module
  !===================================================================

  ! ----------- Default values for input parameters -------------
  real(kind=r8_kind) :: df_energy_criterion  = 1.0E-5_r8_kind, &
                        df_coeff_criterion   = 1.0E-5_r8_kind, &
                        df_coulomb_criterion = 1.0E-5_r8_kind, &
                        df_density_criterion = 1.0E-5_r8_kind
  integer(kind=i4_kind) :: df_max_iteration       = 100, &
                           df_max_geo_iteration   =   1, &
                           df_energy_dev_checked  =   4, &
                           df_coeff_dev_checked   =   1, &
                           df_coulomb_dev_checked =   0, &
                           df_density_dev_checked =   1
  logical :: df_stop_if_not_converged = .false.

  ! ----------- convergence input parameter ----------------------
  real(kind=r8_kind) :: energy_criterion , &
                        coeff_criterion  , &
                        coulomb_criterion, &
                        density_criterion
  integer(kind=i4_kind) :: max_iteration      , &
                           max_geo_iteration  , &
                           energy_dev_checked , &
                           coeff_dev_checked  , &
                           coulomb_dev_checked, &
                           density_dev_checked
  logical :: stop_if_not_converged

  namelist /convergence_list/ energy_criterion     , &
                              coeff_criterion      , &
                              coulomb_criterion    , &
                              density_criterion    , &
                              max_iteration        , &
                              max_geo_iteration    , &
                              energy_dev_checked   , &
                              coeff_dev_checked    , &
                              coulomb_dev_checked  , &
                              density_dev_checked  , &
                              stop_if_not_converged

  ! ----------- local variables and arrays -----------------------
  integer(kind=i4_kind)       :: energy_conv_index , energy_kept_index , &
                                 coeff_conv_index  , coeff_kept_index  , &
                                 coulomb_conv_index, coulomb_kept_index, &
                                 density_conv_index, density_kept_index

  real(kind=r8_kind), pointer :: energy_conv_list(:) , energy_kept_list(:) , &
                                 coeff_conv_list(:)  , coeff_kept_list(:)  , &
                                 coulomb_conv_list(:), coulomb_kept_list(:), &
                                 density_conv_list(:), density_kept_list(:)
  ! auxiliary variable needed to keep change freom previous
  ! geometry loop
  integer(kind=i4_kind)       :: maxgeo_change = huge(1_i4_kind)
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************

  subroutine convergence_setup()
    ! purpose: resets private variables at start of scf part
    !          and allocates plus initializes the convergence buffers
    !** End of interface ***************************************
    integer :: alloc_stat

    allocate(energy_conv_list(energy_dev_checked+1), &
             coeff_conv_list(coeff_dev_checked)    , &
             coulomb_conv_list(coulomb_dev_checked), &
             density_conv_list(density_dev_checked), stat=alloc_stat)
    if (alloc_stat /= 0 ) call error_handler &
         ("convergence_setup: allocation failed")

    energy_conv_index = 1
    energy_conv_list = huge(0.0_r8_kind)
    coeff_conv_index = 1
    coeff_conv_list = huge(0.0_r8_kind)
    coulomb_conv_index = 1
    coulomb_conv_list = huge(0.0_r8_kind)
    density_conv_index = 1
    density_conv_list = huge(0.0_r8_kind)

  end subroutine convergence_setup

  !*************************************************************

  subroutine convergence_shutdown()
    ! purpose: deallocates the convergence buffers
    !** End of interface ***************************************
    integer :: alloc_stat

    deallocate(energy_conv_list , &
               coeff_conv_list  , &
               coulomb_conv_list, &
               density_conv_list, stat=alloc_stat)
    if (alloc_stat /= 0) call error_handler &
         ("convergence_convergence_shutdown: deallocation failed")

  end subroutine convergence_shutdown

  !*************************************************************

  subroutine convergence_read_scfcontrol(warning,new_trace_format)
    ! purpose: read in the namelist convergence from the
    !          input file  in $TTFSSTART/scfcontrol
    !          at each scf-cycle, giving a notice if
    !          parameters were changed.
    use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters -------------
    logical, intent(in   ), optional :: warning ! if present and .false. warning
                                                ! messages are suppressed
    logical, intent(inout), optional :: new_trace_format
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: max_iteration_old
    integer(kind=i4_kind) :: max_geo_iteration_old
    integer(kind=i4_kind) :: energy_dev_checked_old
    integer(kind=i4_kind) :: coeff_dev_checked_old
    integer(kind=i4_kind) :: coulomb_dev_checked_old
    integer(kind=i4_kind) :: density_dev_checked_old
    real(kind=r8_kind)    :: energy_criterion_old
    real(kind=r8_kind)    :: coeff_criterion_old
    real(kind=r8_kind)    :: coulomb_criterion_old
    real(kind=r8_kind)    :: density_criterion_old
    logical               :: stop_if_not_converged_old
    !------------ Executable code ------------------------------

    ! save input variables
    max_iteration_old         = max_iteration
    max_geo_iteration_old     = max_geo_iteration
    energy_dev_checked_old    = energy_dev_checked
    coeff_dev_checked_old     = coeff_dev_checked
    coulomb_dev_checked_old   = coulomb_dev_checked
    density_dev_checked_old   = density_dev_checked
    energy_criterion_old      = energy_criterion
    coeff_criterion_old       = coeff_criterion
    coulomb_criterion_old     = coulomb_criterion
    density_criterion_old     = density_criterion
    stop_if_not_converged_old = stop_if_not_converged


    call convergence_read(scfcontrol=.true.)  !!!!!!!!!!!!!!!!!!AS
    if (present(new_trace_format)) then
       new_trace_format = new_trace_format .or. &
            (coeff_dev_checked   > 0 .neqv. coeff_dev_checked_old   > 0) .or. &
            (coulomb_dev_checked > 0 .neqv. coulomb_dev_checked_old > 0) .or. &
            (density_dev_checked > 0 .neqv. density_dev_checked_old > 0)
    endif
    if (present(warning)) then
       if (.not.warning) return
    endif

    ! print warnings on those input variables which have changed
    call check_intg("MAX_ITERATION"        , &
                     max_iteration         , max_iteration_old        )
    call check_intg("MAX_GEO_ITERATION"    , &
                     max_geo_iteration     , max_geo_iteration_old    )
    call check_intg("ENERGY_DEV_CHECKED"   , &
                     energy_dev_checked    , energy_dev_checked_old   )
    call check_intg("COEFF_DEV_CHECKED"    , &
                     coeff_dev_checked     , coeff_dev_checked_old    )
    call check_intg("COULOMB_DEV_CHECKED"  , &
                     coulomb_dev_checked   , coulomb_dev_checked_old  )
    call check_intg("DENSITY_DEV_CHECKED"  , &
                     density_dev_checked   , density_dev_checked_old  )
    call check_real("ENERGY_CRITERION"     , &
                     energy_criterion      , energy_criterion_old     )
    call check_real("COEFF_CRITERION"      , &
                     coeff_criterion       , coeff_criterion_old      )
    call check_real("COULOMB_CRITERION"    , &
                     coulomb_criterion     , coulomb_criterion_old    )
    call check_real("DENSITY_CRITERION"    , &
                     density_criterion     , density_criterion_old    )
    call check_flag("STOP_IF_NOT_CONVERGED", &
                     stop_if_not_converged , stop_if_not_converged_old)

    ! resize convergence buffers without deallocating their original
    ! contents which are accessible via the <data>_kept_list pointers
    !
    ! NO: call convergence_resize_buffers()
    !
    ! This is more than just reading scfcontrol and was extracted
    ! into a separate sub that is called every scf iteration directly
    ! from main_scf()

  contains

    subroutine check_intg(name,val,val_old)
       character(len=*)     , intent(in) :: name
       integer(kind=i4_kind), intent(in) :: val, val_old
       if ( val .ne. val_old ) then
          call write_to_output_units("convergence_read_scfcontrol: "// &
               name//" was altered. New value: ",inte=val)
          call write_to_trace_unit("convergence_read_scfcontrol: "// &
               name//" was altered. New value: ",inte=val)
          ! for max_geo_iteration a change-variable has to be set
          ! since the input is read in for each geometry looop,
          ! which in turn would override the change made
          ! in this routine. Sorry for messing up this routine, Uwe.
          if (name=='MAX_GEO_ITERATION') then
             maxgeo_change=max_geo_iteration
          endif
       endif
    end subroutine check_intg

    subroutine check_flag(name,val,val_old)
       character(len=*)     , intent(in) :: name
       logical, intent(in) :: val, val_old
       if ( val .neqv. val_old ) then
          if (val) then
             call write_to_output_units("convergence_read_scfcontrol: Flag "// &
                  name//" was switched on.")
             call write_to_trace_unit("convergence_read_scfcontrol: Flag "// &
                  name//" was switched on.")
          else
             call write_to_output_units("convergence_read_scfcontrol: Flag "// &
                  name//" was switched off.")
             call write_to_trace_unit("convergence_read_scfcontrol: Flag "// &
                  name//" was switched off.")
          endif
       endif
     end subroutine check_flag

    subroutine check_real(name,val,val_old)
       character(len=*)  , intent(in) :: name
       real(kind=r8_kind), intent(in) :: val, val_old
       if ( val .ne. val_old ) then
          call write_to_output_units("convergence_read_scfcontrol: "// &
               name//" was altered. New value: ",re=val)
          call write_to_trace_unit("convergence_read_scfcontrol: "// &
               name//" was altered. New value: ",real=val)
       endif
    end subroutine check_real

  end subroutine convergence_read_scfcontrol

  !*************************************************************

  subroutine convergence_resize_buffers()
    ! resize convergence buffers without deallocating their original
    ! contents which are accessible via the <data>_kept_list pointers
    implicit none
    ! *** end of interface ***

    energy_kept_index = energy_conv_index
    energy_kept_list => energy_conv_list
    call resize("ENERGY_CONV_LIST", &
         energy_conv_list,energy_conv_index,energy_dev_checked+1)
    coeff_kept_index = coeff_conv_index
    coeff_kept_list => coeff_conv_list
    call resize("COEFF_CONV_LIST", &
         coeff_conv_list,coeff_conv_index,coeff_dev_checked)
    coulomb_kept_index = coulomb_conv_index
    coulomb_kept_list => coulomb_conv_list
    call resize("COULOMB_CONV_LIST", &
         coulomb_conv_list,coulomb_conv_index,coulomb_dev_checked)
    density_kept_index = density_conv_index
    density_kept_list => density_conv_list
    call resize("DENSITY_CONV_LIST", &
         density_conv_list,density_conv_index,density_dev_checked)
  end subroutine convergence_resize_buffers

  !*************************************************************

  subroutine resize(name,list,index,new_size)
     character(len=*)     , intent(in)    :: name
     real(kind=r8_kind)   , pointer       :: list(:)
     integer(kind=i4_kind), intent(inout) :: index
     integer(kind=i4_kind), intent(in)    :: new_size

     integer(kind=i4_kind)       :: i, offset, alloc_stat
     real(kind=r8_kind), pointer :: orig(:)

     offset = new_size - size(list)
     if (offset == 0) return
     orig => list
     allocate(list(new_size),stat=alloc_stat)
     if (alloc_stat /= 0) call error_handler &
          ("convergence_read_scfcontrol: allocation "//name//" failed")
     do i=new_size,max(offset+1,1),-1
        index = index - 1
        if (index == 0) index = size(orig)
        list(i) = orig(index)
     enddo
     list(:offset) = huge(0.0_r8_kind)
     index = 1
     nullify(orig)
  end subroutine resize

  !*************************************************************

  subroutine convergence_clear_buffers()
    ! purpose: releases the original data of the convergence buffers
    !** End of interface ***************************************
    !------------ Declaration of local variables   ---------------
    integer :: alloc_stat
    !------------ Executable code ------------------------------

    ! Caution: ASSOCIATED(PTR1,PTR2) yields .FALSE. for zero-sized
    !          pointers even if PTR1 and PTR2 are associated to each other

    if (size(energy_kept_list) /= size(energy_conv_list)) then
       deallocate(energy_kept_list,stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler &
            ("convergence_clear_buffers: deallocation (1) failed")
    endif

    if (size(coeff_kept_list) /= size(coeff_conv_list)) then
       deallocate(coeff_kept_list,stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler &
            ("convergence_clear_buffers: deallocation (2) failed")
    endif

    if (size(coulomb_kept_list) /= size(coulomb_conv_list)) then
       deallocate(coulomb_kept_list,stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler &
            ("convergence_clear_buffers: deallocation (3) failed")
    endif

    if (size(density_kept_list) /= size(density_conv_list)) then
       deallocate(density_kept_list,stat=alloc_stat)
       if (alloc_stat /= 0) call error_handler &
            ("convergence_clear_buffers: deallocation (4) failed")
    endif

  end subroutine convergence_clear_buffers

  !*************************************************************

  subroutine convergence_read(scfcontrol)
    ! purpose: read in the namelist convergence from the
    !          input file  in $TTFSSTART/input.
    !** End of interface ***************************************
    use input_module
    use options_module, only : options_max_geo_iter => max_geo_iteration, &
                               options_df_max_geo_iter => df_max_geo_iteration
    use operations_module, only : operations_fitbasis_opt, operations_qm_mm, &
                                  operations_symadapt_gx, operations_potential, &
                                  operations_geo_opt
    !------------ Declaration of formal parameters -------------
    logical, intent(in), optional :: scfcontrol      !!!!!!!!!!!!!!AS
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: unit,status
    !------------ Executable code ------------------------------

    energy_criterion      = df_energy_criterion
    coulomb_criterion     = df_coulomb_criterion
    coeff_criterion       = df_coeff_criterion
    density_criterion     = df_density_criterion
    max_iteration         = df_max_iteration
    max_geo_iteration     = df_max_geo_iteration
    energy_dev_checked    = df_energy_dev_checked
    coeff_dev_checked     = df_coeff_dev_checked
    coulomb_dev_checked   = df_coulomb_dev_checked
    density_dev_checked   = df_density_dev_checked
    stop_if_not_converged = df_stop_if_not_converged

    if ( input_line_is_namelist("convergence_list") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=convergence_list, iostat=status)
       if (status .gt. 0) call input_error( &
            "convergence_read: namelist convergence")
    endif

    ! to ensure compatibility
    if (max_geo_iteration == df_max_geo_iteration .and. &
        options_max_geo_iter /= options_df_max_geo_iter) then
       max_geo_iteration = options_max_geo_iter
    endif
    ! includes changes made in the previous geometry loop
    if (maxgeo_change/=huge(1_i4_kind)) then
       max_geo_iteration=maxgeo_change
    endif
    if (operations_fitbasis_opt) then
       max_iteration     = 1
       max_geo_iteration = 1
    endif
    if (operations_potential) max_geo_iteration = 1 !!!!!!!!!!!!!!!!AS
    if (.not.operations_geo_opt) max_geo_iteration = 1 !!!!!!!!!!!!!!!!AS
    if (operations_qm_mm .or. operations_symadapt_gx) max_geo_iteration = 1

    ! Only check if max_iteration is positive...
    ! Any upper restriction seems to be contraproductive today
    if (max_iteration .lt. 1) then
       call input_error &
            ("CONVERGENCE_READ: max_iteration too small")
    endif
    if (max_geo_iteration .ge.1000) then
       call input_error &
            ("CONVERGENCE_READ: max_geo_iteration too large")
    elseif (max_geo_iteration < 0 ) then
       ! use max_geo_iteration == 0 if you want gradients to be written
       ! to gxfile but dont want to run optimizer at all.
       call input_error &
            ("CONVERGENCE_READ: max_geo_iteration too small")
    endif
    if (energy_criterion .ge. 1.0_r8_kind) then
       call input_error &
            ("CONVERGENCE_READ: energy_criterion too large")
    elseif (energy_criterion .lt. 1.0E-16_r8_kind) then
       call input_error &
            ("CONVERGENCE_READ: energy_criterion too small")
    endif
    if (coeff_criterion .gt. 1.0_r8_kind) then
       call input_error &
            ("CONVERGENCE_READ: charge_coeff criterion too large")
    endif
    if (coulomb_criterion .gt. 1.0_r8_kind) then
       call input_error &
            ("CONVERGENCE_READ: charge_coulomb criterion too large")
    endif
    if (density_criterion .gt. 1.0_r8_kind) then
       call input_error &
            ("CONVERGENCE_READ: density criterion too large")
    endif
    if (energy_dev_checked .lt. 0) then
       call input_error &
            ("CONVERGENCE_READ: negative number of energy_devs checked")
    endif
    if (coeff_dev_checked .lt. 0) then
       call input_error &
            ("CONVERGENCE_READ: negative number of coeff_devs checked")
    endif
    if (coulomb_dev_checked .lt. 0) then
       call input_error &
            ("CONVERGENCE_READ: negative number of coulomb_devs checked")
    endif
    if (density_dev_checked .lt. 0) then
       call input_error &
            ("CONVERGENCE_READ: negative number of density_devs checked")
    endif
    if ( .not.convergence_check_energy_dev() .and. &
         .not.convergence_check_coeff_dev() .and. &
         .not.convergence_check_coulomb_dev() .and. &
         .not.convergence_check_density_dev() ) then
       call input_error &
            ("CONVERGENCE_READ: at least one convergence criterion be used")
    endif
  end subroutine convergence_read

  !*************************************************************

  subroutine convergence_write(iounit, full, scfcontrol)
    !
    ! Purpose: write the namelist convergence.
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in)           :: iounit
    logical, intent(in), optional :: full
    logical, intent(in), optional :: scfcontrol      !!!!!!!!!!!!!!AS
    !** End of interface ***************************************

    if (present(full)) then
       call start("CONVERGENCE_LIST","CONVERGENCE_WRITE", &
            iounit,echo_level_full)
    else
       call start("CONVERGENCE_LIST","CONVERGENCE_WRITE", &
            iounit,operations_echo_input_level)

    endif
    call intg("MAX_ITERATION        ", &
               max_iteration         , df_max_iteration        )
    call real("ENERGY_CRITERION     ", &
               energy_criterion      , df_energy_criterion     )
    call intg("ENERGY_DEV_CHECKED   ", &
               energy_dev_checked    , df_energy_dev_checked   )
    call real("COEFF_CRITERION      ", &
               coeff_criterion       , df_coeff_criterion      )
    call intg("COEFF_DEV_CHECKED    ", &
               coeff_dev_checked     , df_coeff_dev_checked    )
    call real("COULOMB_CRITERION    ", &
               coulomb_criterion     , df_coulomb_criterion    )
    call intg("COULOMB_DEV_CHECKED  ", &
               coulomb_dev_checked   , df_coulomb_dev_checked  )
    call real("DENSITY_CRITERION    ", &
               density_criterion     , df_density_criterion    )
    call intg("DENSITY_DEV_CHECKED  ", &
               density_dev_checked   , df_density_dev_checked  )
    call intg("MAX_GEO_ITERATION    ", &
               max_geo_iteration     , df_max_geo_iteration    )
    call flag("STOP_IF_NOT_CONVERGED", &
               stop_if_not_converged , df_stop_if_not_converged)

    if (present(scfcontrol)) then
       if (scfcontrol) then
          call stop(empty_line=.false.) !!!!!!!!!!!!!!!AS
       else
          call stop()
       end if
    else
       call stop()
    endif

  end subroutine convergence_write

  !*************************************************************

  function convergence_load_coeff_dev(coeff_old, coeff_new, n_output, text)
    !
    ! Purpose: gives the current deviation of the fit coefficients.
    !
    ! Input:
    !
    ! coeff_old - fit coefficients from last cycle
    ! coeff_new - fit coefficients from present cycle
    ! n_output  - number of deviations to be printed
    ! text      - header for the output
    !
    use iounitadmin_module, only: output_unit
    use fit_coeff_module, only: ff_map => fit_coeff_ff_map
    implicit none
    real(kind=r8_kind)   , intent(in) :: coeff_old(:),coeff_new(:)
    integer(kind=i4_kind), intent(in) :: n_output
    character(len=*)     , intent(in) :: text
    real(kind=r8_kind)                :: convergence_load_coeff_dev
    !** End of interface ***************************************

    integer(kind=i4_kind) :: i_cf, i_dev, j_dev, n_dev, alloc_stat
    real(kind=r8_kind)    :: dev0
    real(kind=r8_kind)   , allocatable :: dev(:)
    integer(kind=i4_kind), allocatable :: ind(:)
    integer(kind=i4_kind) :: n_ff


    n_ff = size(ff_map)
    ASSERT(n_ff==size(coeff_old))
    ASSERT(n_ff==size(coeff_new))

    n_dev = max(n_output,1)
    allocate(dev(n_dev),ind(n_dev),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("convergence_load_coeff_dev: allocation failed")

    !
    ! When the fit  basis is empty, the loop  in the else-branch below
    ! is  not  even  entered,  so  that dev  remains  at  its  initial
    ! value. The  fit is not used  in HF/Hybrid runs. Zero  might be a
    ! more meaningful values in this case:
    !
    ! dev = -1.0_r8_kind
    dev = 0.0
    ind = 0

    if (n_dev == 1) then
       !
       ! FIXME: maxval()  is not  well defined for  zero-length arrays
       ! and will  probably return some very negative  number when the
       ! fit basis is empty:
       !
       dev = maxval( &
            abs( &
              (coeff_old * sqrt(ff_map(:)%COUL_NORM)) &
            - (coeff_new * sqrt(ff_map(:)%COUL_NORM)) &
            ))
       ! if %renorm is on, %COUL_NORM should be anyway 1.0
       ! if %renorm is off, adjust the diff.
       ! parenthes necessary, because of numeric diffs
    else
       do i_cf=1,size(coeff_old)
          dev0 = abs( &
                  (coeff_old(i_cf) * sqrt(ff_map(i_cf)%COUL_NORM)) &
                - (coeff_new(i_cf) * sqrt(ff_map(i_cf)%COUL_NORM)) &
                )
          if (dev0 > dev(n_dev)) then
             ! find position of dev0 in dev(1) >= ... >= dev(n_dev)
             i_dev = 1
             do while (dev0 <= dev(i_dev))
                i_dev = i_dev + 1
             enddo
             do j_dev=n_dev,i_dev+1,-1
                dev(j_dev) = dev(j_dev-1)
                ind(j_dev) = ind(j_dev-1)
             enddo
             dev(i_dev) = dev0
             ind(i_dev) = i_cf
          endif
       enddo
    endif

    convergence_load_coeff_dev = dev(1)

    if (n_output > 0) then
       write(output_unit,*)" "
       write(output_unit,*)text//"   (delta/index)"
       write(output_unit,1000)(dev(i_dev),ind(i_dev),i_dev=1,n_dev)
       write(output_unit,*)" "
       1000 format((4(ES13.3,' /',I4)))
    endif

    deallocate(dev,ind,stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("convergence_load_coeff_dev: deallocation failed")

  end function convergence_load_coeff_dev

  !*************************************************************

  function convergence_load_metric_dev(coeff_old,coeff_new,metric,text)
    ! purpose: gives the current deviation of the fit coefficients
    !          within a given metric tensor (fit function overlap matrix)
    ! input  : coeff_old - fit coefficients from last cycle
    !          coeff_new - fit coefficients from present cycle
    !          metric    - the fit function overlap matrix
    !          text      - header for the output
    !----- Modules used --------------------------------------
    use iounitadmin_module, only : output_unit
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), intent(in) :: coeff_old(:),coeff_new(:)
    real(kind=r8_kind), intent(in) :: metric(:)
    character(len=*)  , intent(in) :: text
    real(kind=r8_kind)             :: convergence_load_metric_dev
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: i_cf, j_cf, n_cf, i_met, alloc_stat
    real(kind=r8_kind)    :: dev
    real(kind=r8_kind), allocatable :: coeff_diff(:)
    !------------ Executable code ------------------------------


    n_cf = ubound(coeff_old,1)
    allocate(coeff_diff(n_cf),stat=alloc_stat)
    if (alloc_stat /= 0) call error_handler &
         ("convergence_load_metric_dev: allocation failed")
    coeff_diff = coeff_new - coeff_old

    dev = 0.0_r8_kind
    i_met = 0
    do i_cf=1,n_cf
       do j_cf=1,i_cf-1 ! strict upper triangle only
          i_met = i_met + 1
          dev = dev + ( coeff_diff(i_cf)*coeff_diff(j_cf) + &
                        coeff_diff(j_cf)*coeff_diff(j_cf) ) * metric(i_met)
       enddo
       i_met = i_met + 1
       dev = dev + coeff_diff(i_cf)*coeff_diff(i_cf) * metric(i_met)
    enddo
    dev = sqrt(dev)

    convergence_load_metric_dev = dev

    write(output_unit,1000)text//" = ",dev
    1000 format(/1X,A,ES13.3/)

    deallocate(coeff_diff,stat=alloc_stat)
    if (alloc_stat /= 0) call error_handler &
         ("convergence_load_metric_dev: deallocation failed")

  end function convergence_load_metric_dev

  !*************************************************************

  integer(kind=i4_kind) function convergence_max_geo_iteration()
  ! number of geometry loops that 'main_master' will perform
  !** End of interface *****************************************
    convergence_max_geo_iteration = max_geo_iteration
  end function convergence_max_geo_iteration

  !*************************************************************

  logical function convergence()
    !
    ! Returns true  if convergence is reached.   Subroutine called by:
    ! main_scf(). FIXME: currently seems to return different values on
    ! master  and slaves.  Is  executed in  parallel context  at least
    ! once.
    !
    !** End of interface ***************************************

    convergence = convergence_energy()  .and. &
                  convergence_coeff()   .and. &
                  convergence_coulomb() .and. &
                  convergence_density()
  end function convergence

  !*************************************************************

  logical function convergence_max_iterations(loop)
    ! purpose: returns if maximal number of itterations is reached.
    ! subroutine called by: main_scf
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind) :: loop
    !** End of interface ***************************************
    convergence_max_iterations = loop .ge. max_iteration
  end function convergence_max_iterations

  !*************************************************************

  logical function convergence_abort_calculation()
    ! purpose: returns if calculation should be aborted because
    ! flag stop_if_not_converged is set and maximal number
    ! of scf cycles was exceeded.
    ! subroutine called by main_scf after end of scf cycles
    ! and output of spectrum.
    !** End of interface ***************************************
    convergence_abort_calculation = stop_if_not_converged .and. &
         .not. convergence()
  end function convergence_abort_calculation

  !*************************************************************

  logical function convergence_energy()
    ! purpose: check convergence on the total energy
    ! UB 8/97
    !** End of interface ***************************************
    !------------ Declaration of local variables  --------------
    real(kind=r8_kind) :: e_min, e_max
    !------------ Executable code ------------------------------
    if (convergence_check_energy_dev()) then
       e_min = minval(energy_conv_list)
       e_max = maxval(energy_conv_list)

       if (e_max == huge(0.0_r8_kind)) then
          convergence_energy = .false.
       else
          convergence_energy = e_max - e_min <= energy_criterion
       endif
    else
       convergence_energy = .true.
    endif
  end function convergence_energy

  !*************************************************************

  logical function convergence_coeff()
    ! purpose: check convergence on the absolute differences of
    !          the charge fit coefficients
    ! UB 8/97
    !** End of interface ***************************************
    !------------ Executable code ------------------------------
    if (convergence_check_coeff_dev()) then
       convergence_coeff = maxval(coeff_conv_list) <= coeff_criterion
    else
       convergence_coeff = .true.
    endif
  end function convergence_coeff

  !*************************************************************

  logical function convergence_coulomb()
    ! purpose: check convergence on the self-interaction of the
    !          fitted charge density differences
    ! UB 8/97
    !** End of interface ***************************************
    !------------ Executable code ------------------------------
    if (convergence_check_coulomb_dev()) then
       convergence_coulomb = maxval(coulomb_conv_list) <= coulomb_criterion
    else
       convergence_coulomb = .true.
    endif
  end function convergence_coulomb

  !*************************************************************

  logical function convergence_density()
    ! purpose: check convergence on the absolute differences of
    !          the density matrix elements
    ! UB 8/97
    !** End of interface ***************************************
    !------------ Executable code ------------------------------
    if (convergence_check_density_dev()) then
       convergence_density = maxval(density_conv_list) <= density_criterion
    else
       convergence_density = .true.
    endif
  end function convergence_density

  !*************************************************************

  logical function convergence_check_energy_dev()
    ! purpose: .TRUE. if total energy difference are checked for convergence
    !** End of interface ***************************************
    convergence_check_energy_dev = energy_dev_checked > 0
  end function convergence_check_energy_dev

  !*************************************************************

  logical function convergence_check_coeff_dev()
    ! purpose: .TRUE. if the total energy difference should be
    !                 checked to determine SCF convergence
    !** End of interface ***************************************
    convergence_check_coeff_dev = coeff_dev_checked > 0
  end function convergence_check_coeff_dev

  !*************************************************************

  logical function convergence_check_coulomb_dev()
    ! purpose: .TRUE. if absolute differences of charge fit coefficients
    !                 should be checked to determine SCF convergence
    !** End of interface ***************************************
    convergence_check_coulomb_dev = coulomb_dev_checked > 0
  end function convergence_check_coulomb_dev

  !*************************************************************

  logical function convergence_check_density_dev()
    ! purpose: .TRUE. if self-interaction of fitted charge density difference
    !                 should be checked to determine SCF convergence
    !** End of interface ***************************************
    convergence_check_density_dev = density_dev_checked > 0
  end function convergence_check_density_dev

  !*************************************************************

  subroutine convergence_put_energy(energy)
    ! purpose: puts energy into the total energy convergence buffer
    !          (called by main_scf)
    ! UB 8/97
    !------------ Declaration of formal parameter --------------
    real(kind=r8_kind), intent(in) :: energy
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: energy_num_conv
    !------------ Executable code ------------------------------
    energy_num_conv = size(energy_conv_list)
    if (energy_num_conv > 0) then
       energy_conv_list(energy_conv_index) = energy
       energy_conv_index = mod(energy_conv_index,energy_num_conv) + 1
    endif
  end subroutine convergence_put_energy

  !*************************************************************

  subroutine convergence_put_coeff_dev(coeff_dev)
    ! purpose: puts coeff_dev into the chargefit coefficient deviation buffer
    !          (called by main_scf)
    ! UB 8/97
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), intent(in) :: coeff_dev
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: coeff_num_conv
    !------------ Executable code ------------------------------
    coeff_num_conv = size(coeff_conv_list)
    if (coeff_num_conv > 0) then
       coeff_conv_list(coeff_conv_index) = coeff_dev
       coeff_conv_index = mod(coeff_conv_index,coeff_num_conv) + 1
    endif
  end subroutine convergence_put_coeff_dev

  !*************************************************************

  subroutine convergence_put_coulomb_dev(coulomb_dev)
    ! purpose: puts coulomb_dev into the fitted density difference buffer
    !          (called by main_scf)
    ! UB 8/97
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), intent(in) :: coulomb_dev
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: coulomb_num_conv
    !------------ Executable code ------------------------------
    coulomb_num_conv = size(coulomb_conv_list)
    if (coulomb_num_conv > 0) then
       coulomb_conv_list(coulomb_conv_index) = coulomb_dev
       coulomb_conv_index = mod(coulomb_conv_index,coulomb_num_conv) + 1
    endif
  end subroutine convergence_put_coulomb_dev

  !*************************************************************

  subroutine convergence_put_density_dev(density_dev)
    ! purpose: puts density_dev into the density matrix convergence buffer
    !          (called by main_scf)
    ! UB 8/97
    !** End of interface ***************************************
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), intent(in) :: density_dev
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: density_num_conv
    !------------ Executable code ------------------------------
    density_num_conv = size(density_conv_list)
    if (density_num_conv > 0) then
       density_conv_list(density_conv_index) = density_dev
       density_conv_index = mod(density_conv_index,density_num_conv) + 1
    endif
  end subroutine convergence_put_density_dev

  !*************************************************************

  subroutine convergence_state_store(th,mode)
    ! Purpose: stores the current state of the internal buffer of
    !          the convergence_module in case any of the options
    !          "save_scfstate", "save_ksmatrix" of "save_eigenvec"
    !          is used.
    !
    ! << th >>     << mode >>    action
    ! PRESENT      NOT PRESENT   store present state immediatelly
    ! NOT PRESENT  PRESENT       keep present state for later storage
    ! PRESENT      PRESENT       stores the previously saved state
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB, 8/97
    !------------ Modules ----------------------------------------
    use options_module, only: recover_ksmatrix, recover_eigenvec, &
                              recover_nothing ! except for the buffers
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_saved
    use readwriteblocked_module
    !------------ Declaration of formal paramaters ---------------
    type(readwriteblocked_tapehandle), optional, intent(inout) :: th
    integer(kind=i4_kind), optional, intent(in) :: mode
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)          :: energy, coeff, coulomb, density
    real(kind=r8_kind), save    :: energy_kept, coeff_kept, &
                                   coulomb_kept, density_kept
    integer(kind=i4_kind)       :: energy_work_index, coeff_work_index, &
                                   coulomb_work_index, density_work_index
    real(kind=r8_kind), pointer :: energy_work_list(:), coeff_work_list(:), &
                                   coulomb_work_list(:), density_work_list(:)
    integer(kind=i4_kind)       :: energy_num_conv, coeff_num_conv, &
                                   coulomb_num_conv, density_num_conv, &
                                   index
    !------------ Executable code ------------------------------------
    if (.not.present(th)) then
       ! keep modus
       if (mode /= recover_ksmatrix) then
          if ( size(energy_conv_list) > 0 ) &
               energy_kept = energy_conv_list(energy_conv_index)
       endif
       if (mode == recover_eigenvec) then
          if ( size(coeff_conv_list) > 0 ) &
               coeff_kept = coeff_conv_list(coeff_conv_index)
          if ( size(coulomb_conv_list) > 0 ) &
               coulomb_kept = coulomb_conv_list(coulomb_conv_index)
          if ( size(density_conv_list) > 0 ) &
             density_kept = density_conv_list(density_conv_index)
       endif
       return
    endif

    if (output_main_scf) call write_to_output_units &
         ("CONVERGENCE_STATE_STORE: saving convergence state")

    if (present(mode)) then
       ! reset modus processes the kept buffers
       energy_work_index  = energy_kept_index
       energy_work_list  => energy_kept_list
       coeff_work_index   = coeff_kept_index
       coeff_work_list   => coeff_kept_list
       coulomb_work_index = coulomb_kept_index
       coulomb_work_list => coulomb_kept_list
       density_work_index = density_kept_index
       density_work_list => density_kept_list
    else
       energy_work_index  = energy_conv_index
       energy_work_list  => energy_conv_list
       coeff_work_index   = coeff_conv_index
       coeff_work_list   => coeff_conv_list
       coulomb_work_index = coulomb_conv_index
       coulomb_work_list => coulomb_conv_list
       density_work_index = density_conv_index
       density_work_list => density_conv_list
    endif
    energy_num_conv  = size(energy_work_list)
    coeff_num_conv   = size(coeff_work_list)
    coulomb_num_conv = size(coulomb_work_list)
    density_num_conv = size(density_work_list)

    if (present(mode)) then
       ! reset modus
       if (mode /= recover_ksmatrix .and. mode /= recover_nothing) then
          if (energy_num_conv > 0) then
             index = energy_work_index - 1
             if (index == 0) index = energy_num_conv
             energy = energy_work_list(index)
             energy_work_list(index) = energy_kept
             energy_work_index       = index
          endif
       endif
       if (mode == recover_eigenvec) then
          if (coeff_num_conv > 0) then
             index = coeff_work_index - 1
             if (index == 0) index = coeff_num_conv
             coeff = coeff_work_list(index)
             coeff_work_list(index) = coeff_kept
             coeff_work_index       = index
          endif
          if (coulomb_num_conv > 0) then
             index = coulomb_work_index - 1
             if (index == 0) index = coulomb_num_conv
             coulomb = coulomb_work_list(index)
             coulomb_work_list(index) = coulomb_kept
             coulomb_work_index       = index
          endif
          if (density_num_conv > 0) then
             index = density_work_index - 1
             if (index == 0) index = density_num_conv
             density = density_work_list(index)
             density_work_list(index) = density_kept
             density_work_index       = index
          endif
       endif
    endif

    ! total energy buffer
    call readwriteblocked_write((/real(energy_num_conv,r8_kind), &
                                  real(energy_work_index,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored convergence state :'
       write(output_unit,'(  a     )')'energy_num_conv,energy_conv_index'
       write(output_unit,'(4es20.13)')(/real(energy_num_conv,r8_kind), &
                                        real(energy_work_index,r8_kind)/)
    endif
    if (energy_num_conv > 0) then
       call readwriteblocked_write(energy_work_list,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'energy_conv_list'
          write(output_unit,'(4es20.13)')energy_work_list
       endif
    endif
    ! charge fit coefficient deviation buffer
    call readwriteblocked_write((/real(coeff_num_conv,r8_kind), &
                                  real(coeff_work_index,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(  a     )')'coeff_num_conv,coeff_conv_index'
       write(output_unit,'(4es20.13)')(/real(coeff_num_conv,r8_kind), &
                                        real(coeff_work_index,r8_kind)/)
    endif
    if (coeff_num_conv > 0) then
       call readwriteblocked_write(coeff_work_list,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coeff_conv_list'
          write(output_unit,'(4es20.13)')coeff_work_list
       endif
    endif
    ! fitted charge density deviation buffer
    call readwriteblocked_write((/real(coulomb_num_conv,r8_kind), &
                                  real(coulomb_work_index,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(  a     )')'coulomb_num_conv,coulomb_conv_index'
       write(output_unit,'(4es20.13)')(/real(coulomb_num_conv,r8_kind), &
                                        real(coulomb_work_index,r8_kind)/)
    endif
    if (coulomb_num_conv > 0) then
       call readwriteblocked_write(coulomb_work_list,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coulomb_conv_list'
          write(output_unit,'(4es20.13)')coulomb_work_list
       endif
    endif
    ! density matrix deviation buffer
    call readwriteblocked_write((/real(density_num_conv,r8_kind), &
                                  real(density_work_index,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(  a     )')'density_num_conv,density_conv_index'
       write(output_unit,'(4es20.13)')(/real(density_num_conv,r8_kind), &
                                        real(density_work_index,r8_kind)/)
    endif
    if (density_num_conv > 0) then
       call readwriteblocked_write(density_work_list,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'density_conv_list'
          write(output_unit,'(4es20.13)')density_work_list
       endif
    endif

    ! restore the state of the working buffers because
    ! <data>_kept_list be may associated with <data>_conv_list
    if (present(mode)) then
       if (mode /= recover_ksmatrix .and. mode /= recover_nothing) then
          if ( energy_num_conv > 0 ) &
               energy_work_list(energy_work_index) = energy
       endif
       if (mode == recover_eigenvec) then
          if ( coeff_num_conv > 0 ) &
               coeff_work_list(coeff_work_index) = coeff
          if ( coulomb_num_conv > 0 ) &
               coulomb_work_list(coulomb_work_index) = coulomb
          if ( density_num_conv > 0 ) &
               density_work_list(density_work_index) = density
       endif
    endif

  end subroutine convergence_state_store

!*************************************************************
! record  1: energy_num_conv,energy_conv_index
! record  2: energy_conv_list                     [if energy_num_conv > 0]
! record  3: coeff_num_conv,coeff_conv_index
! record  4: coeff_conv_list                      [if coeff_num_conv > 0]
! record  5: coulomb_num_conv,coulomb_conv_index
! record  6: coulomb_conv_list                    [if coulomb_num_conv > 0]
! record  7: density_num_conv,density_conv_index
! record  8: density_conv_list                    [if density_num_conv > 0]
!*************************************************************

  subroutine convergence_state_recover(th)
    ! Purpose: recovers the current state from the internal buffer of
    !          the convergence_module in case any of the options
    !          "save_scfstate", "save_ksmatrix" of "save_eigenvec"
    !          is used.
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB 8/97
    !------------ Modules ----------------------------------------
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_read
    use readwriteblocked_module
    use options_module,     only: options_reset_scfcycle
    !------------ Declaration of formal_parameters ---------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    !------------ Declaration of local variables   ---------------
    real(kind=r8_kind)              :: dummy(2)
    real(kind=r8_kind), allocatable :: buffer(:)
    integer(kind=i4_kind)           :: energy_num_conv, coeff_num_conv, &
                                       coulomb_num_conv, density_num_conv
    integer(kind=i4_kind)           :: energy_load_size, coeff_load_size, &
                                       coulomb_load_size, density_load_size
    integer(kind=i4_kind)           :: alloc_stat
    !------------ Executable code ------------------------------------

    if (output_main_scf) call write_to_output_units &
         ("CONVERGENCE_STATE_RECOVER: reading convergence state")

    if (options_reset_scfcycle())then
       energy_load_size = 0
       coeff_load_size = 0
       coulomb_load_size = 0
       density_load_size = 0
    else
       energy_load_size = size(energy_conv_list)
       coeff_load_size = size(coeff_conv_list)
       coulomb_load_size = size(coulomb_conv_list)
       density_load_size = size(density_conv_list)
    endif

    call readwriteblocked_read(dummy,th)
    energy_num_conv   = int(dummy(1),i4_kind)
    energy_conv_index = int(dummy(2),i4_kind)
    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered convergence state :'
       write(output_unit,'(  a     )')'energy_num_conv,energy_conv_index'
       write(output_unit,'(4es20.13)')dummy
    endif
    if (energy_num_conv > 0) then
       if (energy_num_conv == energy_load_size) then
          call readwriteblocked_read(energy_conv_list,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'energy_conv_list'
             write(output_unit,'(4es20.13)')energy_conv_list
          endif
       elseif (energy_load_size > 0) then
          allocate(buffer(energy_num_conv),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: allocation (1) failed")
          call readwriteblocked_read(buffer,th)
          call reload("energy",energy_conv_list,energy_conv_index)
          deallocate(buffer,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: deallocation (1) failed")
       else
          call readwriteblocked_skipread(energy_num_conv,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'energy_conv_list skipped'
          endif
       endif
    else
       if (energy_load_size > 0) then
          if (output_data_read) then
             write(output_unit,'(  a     )')'energy_conv_list as initialized'
             write(output_unit,'(4es20.13)')energy_conv_list
          endif
       endif
    endif

    call readwriteblocked_read(dummy,th)
    coeff_num_conv   = int(dummy(1),i4_kind)
    coeff_conv_index = int(dummy(2),i4_kind)
    if (output_data_read) then
       write(output_unit,'(  a     )')'coeff_num_conv,coeff_conv_index'
       write(output_unit,'(4es20.13)')dummy
    endif
    if (coeff_num_conv > 0) then
       if (coeff_num_conv == coeff_load_size) then
          call readwriteblocked_read(coeff_conv_list,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_conv_list'
             write(output_unit,'(4es20.13)')coeff_conv_list
          endif
       elseif (coeff_load_size > 0) then
          allocate(buffer(coeff_num_conv),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: allocation (2) failed")
          call readwriteblocked_read(buffer,th)
          call reload("coeff",coeff_conv_list,coeff_conv_index)
          deallocate(buffer,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: deallocation (2) failed")
       else
          call readwriteblocked_skipread(coeff_num_conv,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_conv_list skipped'
          endif
       endif
    else
       if (coeff_load_size > 0) then
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_conv_list as initialized'
             write(output_unit,'(4es20.13)')coeff_conv_list
          endif
       endif
    endif

    call readwriteblocked_read(dummy,th)
    coulomb_num_conv   = int(dummy(1),i4_kind)
    coulomb_conv_index = int(dummy(2),i4_kind)
    if (output_data_read) then
       write(output_unit,'(  a     )')'coulomb_num_conv,coulomb_conv_index'
       write(output_unit,'(4es20.13)')dummy
    endif
    if (coulomb_num_conv > 0) then
       if (coulomb_num_conv == coulomb_load_size) then
          call readwriteblocked_read(coulomb_conv_list,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coulomb_conv_list'
             write(output_unit,'(4es20.13)')coulomb_conv_list
          endif
       elseif (coulomb_load_size > 0) then
          allocate(buffer(coulomb_num_conv),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: allocation (3) failed")
          call readwriteblocked_read(buffer,th)
          call reload("coulomb",coulomb_conv_list,coulomb_conv_index)
          deallocate(buffer,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: deallocation (3) failed")
       else
          call readwriteblocked_skipread(coulomb_num_conv,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coulomb_conv_list skipped'
          endif
       endif
    else
       if (coulomb_load_size > 0) then
          if (output_data_read) then
             write(output_unit,'(  a     )')'coulomb_conv_list as initialized'
             write(output_unit,'(4es20.13)')coulomb_conv_list
          endif
       endif
    endif

    call readwriteblocked_read(dummy,th)
    density_num_conv   = int(dummy(1),i4_kind)
    density_conv_index = int(dummy(2),i4_kind)
    if (output_data_read) then
       write(output_unit,'(  a     )')'density_num_conv,density_conv_index'
       write(output_unit,'(4es20.13)')dummy
    endif
    if (density_num_conv > 0) then
       if (density_num_conv == density_load_size) then
          call readwriteblocked_read(density_conv_list,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'density_conv_list'
             write(output_unit,'(4es20.13)')density_conv_list
          endif
       elseif (density_load_size > 0) then
          allocate(buffer(density_num_conv),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: allocation (4) failed")
          call readwriteblocked_read(buffer,th)
          call reload("density",density_conv_list,density_conv_index)
          deallocate(buffer,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("convergence_recover: deallocation (4) failed")
       else
          call readwriteblocked_skipread(density_num_conv,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'density_conv_list skipped'
          endif
       endif
    else
       if (density_load_size > 0) then
          if (output_data_read) then
             write(output_unit,'(  a     )')'density_conv_list as initialized'
             write(output_unit,'(4es20.13)')density_conv_list
          endif
       endif
    endif

  contains

    subroutine reload(name,list,index)
       character(len=*)     , intent(in   ) :: name
       real(kind=r8_kind)   , intent(inout) :: list(:)
       integer(kind=i4_kind), intent(inout) :: index

       integer(kind=i4_kind) :: i, offset

       offset = size(list) - size(buffer)
       do i=size(list),max(offset+1,1),-1
          index = index - 1
          if (index == 0) index = size(buffer)
          list(i) = buffer(index)
       enddo
       list(:offset) = huge(0.0_r8_kind)
       index = 1

       if (output_data_read) then
          write(output_unit,'(  a     )')name//'_conv_list'
          write(output_unit,'(4es20.13)')buffer
          write(output_unit,'(  a     )')'transformed into'
          write(output_unit,'(4es20.13)')list
          write(output_unit,'(  a     )')'with '//name//'_conv_index = 1'
       endif

    end subroutine reload

  end subroutine convergence_state_recover

  !*************************************************************

  !--------------- End of module --------------------------------
end module convergence_module

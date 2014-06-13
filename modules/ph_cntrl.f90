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
module ph_cntrl
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
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

# include "def.h"
  use type_module, only: i4_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !
  ! This flag is used in post_scf_module to indicate a GGA run:
  !
  logical, public, protected :: nl_calc_ph
  !
  ! FIXME: This flag and compatibility reasons (accepting older inputs)
  !        are the only reason for existence of this module.
  !        Get rid of it, whoever has the courage.
  !

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: post_scf_read_input
  public :: post_scf_write_input
  public :: post_scf_set_defaults
  public :: post_scf_input_bcast
  public :: ph_cntrl_setup

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine ph_cntrl_setup()
    ! called at the beginning of PH
!!$    use operations_module, only: operations_gradients
    implicit none
    ! *** end of interface ***

    ! does nothing currently, should import
    ! XC settings from PHXC into XC_CONTROL
    ! FIXME here, if you want different XC in PH and SCF
  end subroutine ph_cntrl_setup

  subroutine post_scf_read_input()
    !
    ! Kept only for compatibility reasons.
    !
    ! Purpose : Reading the input concerning the phxc_control
    !           The input consits of logicals which determine if
    !           a special functional has to be calculated or not.
    !           if a nonlocal functional is selected the logical
    !           nl_clac is set to true.
    use input_module
    use strings, strlen => stringlength_string
    implicit none
    !** End of interface **************************************
    integer(kind=i4_kind) :: unit,status
    character(len=strlen) :: XC = "undef"

    ! ----------- Post SCF XC control parameters ---------------------------
    logical ::&
         & becke88_ph,&
         & rbecke88_ph,&
         & perdew_ph ,&
         & perdewwang91c_ph ,&
         & rperdewwang91c_ph ,&
         & perdewwang91x_ph ,&
         & xalpha_ph ,&
         & rxalpha_ph ,&
         & vwn_ph,&
         & rvwn_ph,&
         & pbex_ph,&
         & pbec_ph,&
         & revPBEx_ph,&
         & PBENx_ph,&
         & revPW91c_ph,&
         & pwldac_ph, &
         & hcth_x_ph, &
         & hcth_c_ph, &
         & ecmv92_ph, &
         & nrecmv92_ph

    namelist /phxc_control/&
         XC,&
         becke88_ph,        &
         rbecke88_ph,       &
         perdew_ph ,        &
         perdewwang91c_ph , &
         rperdewwang91c_ph, &
         perdewwang91x_ph , &
         xalpha_ph ,        &
         rxalpha_ph ,       &
         vwn_ph,            &
         rvwn_ph,           &
         pbex_ph,           &
         pbec_ph,           &
         revPBEx_ph,        &
         PBENx_ph,          &
         revPW91c_ph,       &
         pwldac_ph,         &
         hcth_x_ph,         &
         hcth_c_ph,         &
         ecmv92_ph,         &
         nrecmv92_ph

    WARN('use of PHXC_CONTROL is deprecated!')

    write(*, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    write(*, *) "XXX                                              XXX"
    write(*, *) "XXX               W A R N I N G !                XXX"
    write(*, *) "XXX                                              XXX"
    write(*, *) "XXX Use of PHXC_CONTROL namelist is deprecated ! XXX"
    write(*, *) "XXX Specify XC functional in XC_CONTROL namelist XXX"
    write(*, *) "XXX                 instead.                     XXX"
    write(*, *) "XXX                                              XXX"
    write(*, *) "XXX   Settings in PHXC_CONTROL are IGNORED !     XXX"
    write(*, *) "XXX                                              XXX"
    write(*, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

    XC = "undef"
    if ( input_line_is_namelist("phxc_control") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=phxc_control, iostat=status)
       if (status .gt. 0) call input_error( &
            "phxc_read: namelist xc_control")
    endif

    ! ignore input, use default values:
    call post_scf_set_defaults()
  end subroutine post_scf_read_input

  subroutine post_scf_set_defaults()
    use xc_cntrl, SCF_ON => IS_ON
    use operations_module, only: operations_gradients
    implicit none
    !** End of interface **************************************

    !
    ! By default we do all GGAs in single-point post-scf run:
    !
    nl_calc_ph = .true.

    !
    ! Energy should be consistent with forces, in this case:
    !
    if ( operations_gradients ) then
        nl_calc_ph = SCF_ON(xc_gga)
    endif
  end subroutine post_scf_set_defaults

  subroutine post_scf_write_input(iounit)
    ! purpose : Writing the input concerning phxc_control to the output
    implicit none
    integer(i4_kind), intent(in) :: iounit
    !** End of interface *****************************************

    WARN('skip writing PHXC_CONTROL!')
  end subroutine post_scf_write_input

  subroutine post_scf_input_bcast()
    use comm, only: comm_bcast
    implicit none
    ! Purpose : broadcasting the phxc_input
    !** End of interface *****************************************

    call comm_bcast(nl_calc_ph)

  end subroutine post_scf_input_bcast

  !--------------- End of module ----------------------------------
end module ph_cntrl

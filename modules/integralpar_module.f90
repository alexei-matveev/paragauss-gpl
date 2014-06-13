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
module  integralpar_module
!---------------------------------------------------------------
!
!  Purpose:  This module contains  steering parameters  that determine
!           which  integrals  are  to  be calculated  and  with  witch
!           options.
!
!           Relativistic calculation  can be  turned on and  off here.
!           Posthos functionality  can be switched on  by changing the
!           parameters contained in this  module.  This module must be
!           peresnt both and master and slaves.
!
!           Routines  for  packing  and  unpacking  are  offered.   If
!           operations_integraltest    from    operations_module    is
!           .true. the  parameters are directly read in  as a namelist
!           from     input    by     integralpar_read().     Otherwise
!           integralpar_set_scf() or  integralpar_set_scf(0 are called
!           before the integral part is started.
!
!           This module it  is USEd at quite a  few places. Any import
!           of  another module here  has a  huge potential  for cyclic
!           dpendencies.  Dont do that (unnecessarily).

!
!  Prerequesits: integralpar_setup must be called at beginning of each
!           integral part to set derived variables.
!
!
!  Module called by: Routines all over and down to the bottom of
!                    the integral part
!
!
!  Author: TB
!  Date:   5/96
!
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification
! Author: MS
! Date:   3/97
! Description: New parameters have been introduced
!              in connection with imlementation of energy gradients
!              integralpar_send_2c
!              integralpar_send_3c
!              integralpar_2cob_ol_grad
!              integralpar_3cob_grad
!
! Modification (Please copy before editing)
! Author: TB
! Date:   7/97
! Description: New parameters for dipole integral run introduced
!
! Modification (Please copy before editing)
! Author: MS
! Date:   7/97
! Description: New parameters for relativistic and relativistic
!              energy gradients have been added
!              Also new subroutines to set parameters for relativistic
!              and relativistic gradients have been added.
!
! Modification (Please copy before editing)
! Author:      Uwe Birkenheuer
! Date:        8/98
! Description: new integral parameter for mixed 2-center overlap
!              integrals (required for the extended MDA) introduced.
!
! Modification (Please copy before editing)
! Author: MM
! Date:   5/98
! Description: New parameters for treatment of spin orbit added
!
! Modification
! Author: AM
! Date:   first: Tue Mar 2 1999
! Description: ...
!
! Modification (Please copy before editing)
! Author: AS
! Date:   11/99
! Description: new parameters for calculating electrostatic potential
!              of molecule have been added
!
! Modification (Please copy before editing)
! Author: AS
! Date:   03/00
! Description: new parameters for calculating electrostatic field
!              of molecule have been added
!
! Modification (Please copy before editing)
! Author: AS
! Date:   05/00
! Description: integral parameters of the solvation electrostatic
!              term has been added
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!----------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
implicit none

save ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ Declaration of constants and variables ------------

! default values for public variables:
logical, parameter :: &
           df_integralpar_2cob_kin         = .true. , &
           df_integralpar_2cob_nuc         = .true. , &
           df_integralpar_2cob_pvsp        = .false., &
           df_integralpar_2cob_pvxp        = .false., &
           df_integralpar_2cob_pvec        = .false., &
           df_integralpar_2cob_ol          = .true. , &
           df_integralpar_2cob_dipole      = .false., &
           df_integralpar_2cch_no          = .true. , &
           df_integralpar_2cch_pre         = .true. , &
           df_integralpar_2cxc_no          = .false., &
           df_integralpar_2c_mixed         = .false., &
           df_integralpar_1cch_no          = .true. , &
           df_integralpar_3c_xc            = .true. , &
           df_integralpar_3c_co            = .true. , &
           df_integralpar_3c_co_resp       = .false., &
           df_integralpar_3c_rcoul_pvsp    = .false., &
           df_integralpar_3c_rcoul_pvxp    = .false., &
           df_integralpar_3c_r2_pvsp       = .false., &
           df_integralpar_3c_r2_pvxp       = .false., &
           df_integralpar_relativistic     = .false., &
           df_integralpar_pseudo           = .false., & !AS
           df_integralpar_spor             = .false., &
           df_integralpar_gradients        = .false., &
           df_integralpar_2cob_ol_grad     = .false., &
           df_integralpar_3cob_grad        = .false., &
           df_integralpar_send_3c          = .true. , &
           df_integralpar_send_2c          = .true. , &
           df_integralpar_offdiag_dipoles  = .false., &
           df_integralpar_gten_kinematic   = .true. , &
           df_integralpar_rel_gradients    = .false., &
           df_integralpar_2cob_kin_grad    = .false., &
           df_integralpar_2cob_ol_rel_grad = .false., &
           df_integralpar_2cob_nuc_grad    = .false., &
           df_integralpar_2cob_pvsp_grad   = .false., &
           df_integralpar_solv_grad        = .false., &
           df_integralpar_Q_solv_grad      = .false., &
           df_integralpar_cpksdervs        = .false., &
           df_integralpar_2dervs           = .false., &
           df_integralpar_dervs            = .false., &
           df_integralpar_cpks_contribs    = .false., &
           df_integralpar_2cob_potential   = .false., &
           df_integralpar_pot_for_secderiv = .false., &
           df_integralpar_2cob_field       = .false., &
#ifdef WITH_EFP
           df_integralpar_2cob_field_efp   = .false., &
           df_integralpar_2cob_pc          = .false., &
           df_integralpar_2cob_X           = .false., &
           df_integralpar_efp_gradients    = .false., &
#endif
           df_integralpar_2cff             = .false., &
           df_integralpar_2cob3c           = .false., &
           df_integralpar_totalsymmetric   = .true.,  &
           df_integralpar_2cob_pc_grad     = .false., &
           df_integralpar_2cob_X_grad      = .false., &
           df_integralpar_2cob_ipd_grad    = .false.

integer(i4_kind), parameter :: df_integralpar_i_int_part = 1
  !----------------------------------------------------------------

! to check that nobody modifies these controls from outside of this
! module declare them ``protected'':
logical, PROTECTED, public  :: &
     integralpar_2cob_kin         = df_integralpar_2cob_kin, &
                                    ! 2-center orbital integral of kinetic energy
     integralpar_2cob_nuc         = df_integralpar_2cob_nuc, &
                                    ! 2-center orbital integral of nuclear attraction
     integralpar_2cob_pvsp        = df_integralpar_2cob_pvsp, &
                                    ! 2-center orbital integral of pvsp ,
                                    !  necessary for relativistic
     integralpar_2cob_pvxp        = df_integralpar_2cob_pvxp, &
                                    ! 2-center orbital integral of pvxp ,
                                    !  necessary for relativistic
     integralpar_2cob_pvec        = df_integralpar_2cob_pvec, &
                                    ! 2-center orbital integral of pvec ,
                                    !  necessary for raltivistic
     integralpar_2cob_ol          = df_integralpar_2cob_ol, &
                                    ! 2-center orbital integral of overlap
     integralpar_2cob_dipole      = df_integralpar_2cob_dipole, &
                                    ! 2-center orbital integral of dipole operator
     integralpar_2cob_potential   = df_integralpar_2cob_potential, &
     integralpar_pot_for_secderiv = df_integralpar_pot_for_secderiv, &
                                    ! 2-center orbital integral of electrostatic potential
     integralpar_2cob_field       = df_integralpar_2cob_field, &
#ifdef WITH_EFP
     integralpar_2cob_field_efp   = df_integralpar_2cob_field_efp, &
     integralpar_2cob_pc          = df_integralpar_2cob_pc, &
     integralpar_2cob_X           = df_integralpar_2cob_X , &
     integralpar_efp_gradients    = df_integralpar_efp_gradients, &
#endif
                                    ! 2-center orbital integral of electrostatic field
     integralpar_2cch_no          = df_integralpar_2cch_no, &
                                    ! 2-center charge norm integral of coulomb self-interaction
     integralpar_2cch_pre         = df_integralpar_2cch_pre, &
                                    ! estimate 2-center charge overlap integrals,
                                    !  used in response part for prescreening
     integralpar_2cxc_no          = df_integralpar_2cxc_no, &
                                    ! 2-center exchange norm integral
     integralpar_2c_mixed         = df_integralpar_2c_mixed, &
                                    ! mixed 2-center charge fit/exchange fit integral
     integralpar_1cch_no          = df_integralpar_1cch_no, &
                                    ! 1-center charge norm integral necessary for chargefit
     integralpar_3c_xc            = df_integralpar_3c_xc, &
                                    ! 3-center exchange integral
     integralpar_3c_co            = df_integralpar_3c_co, &
                                    ! 3-center coulomb integral for common use (totallysymmetric)
     integralpar_3c_co_resp       = df_integralpar_3c_co_resp, &
                                    ! 3-center coulomb integral for response (all symmetry)
     integralpar_3c_rcoul_pvsp    = df_integralpar_3c_rcoul_pvsp, &
                                    ! 3-center scalar relativistic coulomb integral
                                    !  with s-type fit function
     integralpar_3c_rcoul_pvxp    = df_integralpar_3c_rcoul_pvxp, &
                                    ! 3-center vector relativistic coulomb integral
                                    !  with s-type fit function
     integralpar_3c_r2_pvsp       = df_integralpar_3c_r2_pvsp, &
                                    ! 3-center scalar relativistic coulomb integral
                                    !  with r2-type fit function
     integralpar_3c_r2_pvxp       = df_integralpar_3c_r2_pvxp, &
                                    ! 3-center vector relativistic coulomb integral
                                    !  with r2-type fit function
     integralpar_relativistic     = df_integralpar_relativistic, &
                                    ! perform relativistic transformations
     integralpar_pseudo           = df_integralpar_pseudo, & !AS
     integralpar_spor             = df_integralpar_spor, &
                                    ! perform relativistic transformations
     integralpar_gradients        = df_integralpar_gradients, &
                                    ! calculate gradients of matrix elements
     integralpar_2cob_ol_grad     = df_integralpar_2cob_ol_grad, &
                                    ! calculate gradient of overlap matrix
     integralpar_3cob_grad        = df_integralpar_3cob_grad, &
                                    ! calculate derivatives of three center matrix elements
     integralpar_send_3c          = df_integralpar_send_3c, &
                                    ! send the three center integrals via comm
                                    !  has to be set to false for gradient run!
     integralpar_send_2c          = df_integralpar_send_2c, &
                                    ! send the two center integrals via comm
                                    !  has to be set to false for gradient run!
     integralpar_offdiag_dipoles  = df_integralpar_offdiag_dipoles, &
                                    ! are dipole matrix elements between
                                    !  diffrent Irreps required ?
     integralpar_gten_kinematic   = df_integralpar_gten_kinematic, &
                                    ! are kinematic factors have to be
                                    !  accounted for for g-tensor calc
     integralpar_2cff             = df_integralpar_2cff, &
                                    !
                                    ! calculate any 2-center fitfunction integrals
     integralpar_2cob3c           = df_integralpar_2cob3c, &
                                    ! calculate any 2-center orbital or 3 center integrals
     integralpar_totalsymmetric   = df_integralpar_totalsymmetric, &
                                    ! calculate totalsymmeteric integrals
     integralpar_rel_gradients    = df_integralpar_rel_gradients, &
                                    ! do a calculation of relativistic gradients
     integralpar_2cob_kin_grad    = df_integralpar_2cob_kin_grad, &
                                    ! gradients of kinetic energy
     integralpar_2cob_ol_rel_grad = df_integralpar_2cob_ol_rel_grad, &
                                    ! gradients of overlap matrix
     integralpar_2cob_nuc_grad    = df_integralpar_2cob_nuc_grad, &
                                    ! gradients of nuclear attraction
     integralpar_solv_grad        = df_integralpar_solv_grad, &
                                    ! calculate gradient of the electronic part
                                    !  of the electrostatic term to the solvation effect
     integralpar_Q_solv_grad      = df_integralpar_Q_solv_grad, &
                                    ! calculate surface charge contribution to
                                    ! gradient of the electronic part
                                    ! of the electrostatic term to the solvation effect
     integralpar_2cob_pvsp_grad   = df_integralpar_2cob_pvsp_grad, &
                                    ! gradients of pvsp integrals
     integralpar_cpksdervs        = df_integralpar_cpksdervs, &
     integralpar_2dervs           = df_integralpar_2dervs, &
     integralpar_dervs            = df_integralpar_dervs, &
     integralpar_cpks_contribs    = df_integralpar_cpks_contribs, &
                                    ! identifies second gradient run after CPKS
     integralpar_2cob_pc_grad     = df_integralpar_2cob_pc_grad, &
                                    ! gradients on PC
     integralpar_2cob_X_grad      = df_integralpar_2cob_X_grad, &
                                    ! gradients on eXternal centers (PD, PQ, PO)
     integralpar_2cob_ipd_grad    = df_integralpar_2cob_ipd_grad
                                    ! gradients on polarizable centers (IPD)

!
! Index of intergal part performed right now (for timing):
!
integer, public, protected :: integralpar_i_int_part = 1

! Total number  of integral parts (for timing).   FIXME: when changing
! this, adpat  a copy of this  number in timer_module.   There is some
! cyclic   dependency   if   one   imports   integralpar_module   into
! timer_module:
integer, public, parameter :: integralpar_n_int_parts = 4

!
! Names of integral parts for corresponding index:
!
character(len=8), public, parameter :: &
     integralpar_int_part_name(integralpar_n_int_parts) = &
     (/"Normal  ","PostSCF ","Gradient","Dipole  "/)

character(len=16), private, parameter :: df_integralpar_runtype="undef"
character(len=16), public, protected :: integralpar_runtype=df_integralpar_runtype ! one of the above

!------------ public functions and subroutines ------------------
public :: integralpar_set!(Task)
public :: integralpar_send_receive!()
public :: integralpar_setup!()

!================================================================
! End of public interface of module
!================================================================

!................................................................
! integralpar       SCF       PROPER   POSTHOC  GRADS    DIPOLE
! ----------------------------------------------------------------
! offdiag_dipoles   false     false    false    false    <OFFDIAG>
! gradients         false     false    false    true     false
! relativistic      <RELAT>   false    false    <RELAT>  false
! rel_gradients     false     false    false    <RELAT>  false
! send_3c           true      true     false    <RELAT>  false
! send_2c           true      true??   false    false    false
! 2cob_kin          true      false    false    false    false
! 2cob_nuc          true      false    false    false    false
! 2cob_potential    true      false    false    false    false
! 2cob_field        true      false    false    false    false
! 2cob_pvsp         <RELAT>   false    false    false    false
! 2cob_ol           true      true     false    false    false
! 2cob_dipole       false     false    false    false    true
! 2cch_no           true      false    false    true     false
! 2cxc_no           <EX_MDA>  false    false    false    false
! 2c_mixed          <EX-MDA>  false    false    not yet  false
! 1cch_no           true      false    false    false    false
! 3c_xc             <ANY XC>  false    false    false    false
! 3c_co             true      false    false    false    false
! 3c_co_resp        true      false    false    false    false
! renorm_ob         true      true     false    false    false
! renorm_ch         true      false    false    false    false
! renorm_xc         <ANY XC>  false    false    false    false
! 2cob_ol_grad      false     false    false    true     false
! 3cob_grad         false     false    false    true     false
! 2cob_nuc_grad     false     false    false    <RELAT>  false
! 2cob_pvsp_grad    false     false    false    <RELAT>  false
! 2cob_kin_grad     false     false    false    <RELAT>  false
! 2cob_ol_rel_grad  false     false    false    <RELAT>  false
! i_int_part        1         1        2        3        4
!
! <ANY XC> := <XC-FIT> or <EX-MDA>
!................................................................

! The following  variable is set  by the response module  and contains
! the value  of the input parameter  "num_2index_prescreening".  It is
! used by all "integralpar_set_*" subroutines  to set the value of the
! module  public variable  "integralpar_2cch_pre" defined  above. More
! elegant approaches I tried all ran into "cyclic module dependencies"
! HH 7/98. Descriptions appears outdated.
logical, parameter :: integralpar_resp_prescreen = .false.

  !------------ Subroutines ---------------------------------------
contains

  recursive subroutine integralpar_set(task)
    !
    ! To   affect   execution  mode   of   the   integral  part   call
    ! integralpar_set(task) prior to main_integral().
    !
    ! FIXME:  dont import  many things  here, integralpar  is  used in
    ! quite  a few modules.   Any additionl  import here  introduces a
    ! potential  cyclic  dependency.  integralpar_set(taks) should  be
    ! kept  explicit, that  is any  dependence  of its  action on  the
    ! "weather outside" is not welcome.
    !
    use pointcharge_module, only: moving_pc
    use point_dqo_module, only: moving_X_centers, moving_R_centers
    use induced_dipoles_module, only: moving_Pol_centers
    implicit none
    character(len=*), intent(in) :: task
    ! *** end of interface ***

    DPRINT 'integralpar_set(',task,'): entered'
    select case(task)
    case ('defaults')

      call integralpar_set_defaults()

    case ('Normal')

      call integralpar_set('defaults')
      call integralpar_set_scf()

    case ('Properties')

      call integralpar_set('defaults')
      call integralpar_set_properties()

    case ('Gradients')

      call integralpar_set('defaults')
      call integralpar_set_gradient()

      ! FIXME:  not really  a good  idea  to use  anything else  here,
      ! "task" as an argument is pretty explicit an should suffice:
      if(moving_pc) integralpar_2cob_pc_grad=.true.
      if(moving_X_centers .or. moving_R_centers) integralpar_2cob_X_grad=.true.
      if(moving_Pol_centers) integralpar_2cob_ipd_grad=.true.

    case ('SecondDervs')

      ! does not affect other flags, needs to be set in
      ! combination with (before?/after?) the integral runtype
      integralpar_2dervs    = .true.
      integralpar_cpksdervs = .true.

    case ('NoSecondDervs')

      ! does not affect other flags, needs to be set in
      ! combination with (before?/after?) the integral runtype
      integralpar_2dervs    = .false.
      integralpar_cpksdervs = .false.

    case ('Gradients2')

      call integralpar_set('Gradients')
      ! and invoke gradient contraction with CPKS results:
      integralpar_cpks_contribs = .true.
      integralpar_dervs         = .false.

    case ('RelGrads')

      call integralpar_set('Gradients')

      integralpar_relativistic = .true.
      integralpar_rel_gradients=.true.

      ! FIXME:  not really  a good  idea  to use  anything else  here,
      ! "task" as an argument is pretty explicit an should suffice:
      if(moving_pc.or. &
         moving_X_centers .or. &
         moving_R_centers .or. &
         moving_Pol_centers)  integralpar_pseudo=.true. !AS

      ! integrals have to be sent
      integralpar_send_3c=.true. ! FIXME: why?

      ! kin,nuc have to be separated:
      integralpar_2cob_ol_rel_grad=.true.
      integralpar_2cob_kin_grad=.true.
      integralpar_2cob_nuc_grad=.true.
      integralpar_2cob_pvsp_grad=.true.

    ! PLEASE CLARIFY WHO IS WHO HERE, I AM MISSING OVERVIEW:
    ! WHAT IS THE RELATION BETWEEN Potential, Field and Solvation?
    case ('Field')

      call integralpar_set('defaults')
      call integralpar_set_field()

#ifdef WITH_EFP
    case ('Field_at_EFP')

      call integralpar_set('defaults')
      call integralpar_set_field()
      integralpar_send_3c = .false.
      integralpar_2cob_field_efp = .true.

    case ('QM_EFP_energy')

      call integralpar_set('defaults')
      call integralpar_set_potential()
      integralpar_send_3c = .false.
      integralpar_2cob_potential = .false.
      integralpar_2cob_pc = .true.
      integralpar_2cob_X = .true.
      integralpar_2cob_nuc = .true.

   case ('EFP_gradients')
      call integralpar_set('defaults')
      call integralpar_set_field()
      integralpar_send_3c = .false.
      integralpar_2cob_field = .false.
      integralpar_efp_gradients = .true.

      ! FIXME:  not really  a good  idea  to use  anything else  here,
      ! "task" as an argument is pretty explicit an should suffice:
      if(moving_pc) integralpar_2cob_pc_grad=.true.
      if(moving_X_centers .or. moving_R_centers) integralpar_2cob_X_grad=.true.
      if(moving_Pol_centers) integralpar_2cob_ipd_grad=.true.

#endif

    case ('Potential') ! is it the same as Solvation?

      call integralpar_set('defaults')
      call integralpar_set_potential()

    case ('SolvGrads')

      call integralpar_set('defaults')
      call integralpar_set_grad_solv()

    case ('Q_SolvGrads')

      call integralpar_set('defaults')
      call integralpar_set_grad_solv()
      integralpar_solv_grad = .false.
      integralpar_Q_solv_grad = .true.

    case ('SolvGrads2')

      call integralpar_set('defaults')
      call integralpar_set_grad_solv()
      integralpar_cpks_contribs = .true.
      integralpar_dervs         = .false.

    case ('SolvDervs')

      call integralpar_set('defaults')
      call integralpar_set('Potential')
      !  and invoke gradient contraction with CPKS results:
      integralpar_pot_for_secderiv = .true.

    case ('SolvDervs2')

      call integralpar_set('defaults')
      call integralpar_set('Potential')
      !  and invoke gradient contraction with CPKS results:
      integralpar_pot_for_secderiv = .true.
      integralpar_cpks_contribs = .true.
      integralpar_dervs         = .false.

    case ('Dipole')
      call integralpar_set_dipole(offdiag_dipoles=.false.)

    case ('DipoleOff')
      call integralpar_set_dipole(offdiag_dipoles=.true.)

    case default
      print *,'integralpar_set(',task,'): ERROR: no such task!'
      ABORT('no such task, see tty')
    end select

    ! just for references purposes, e.g. for (debug) output:
    integralpar_runtype = task

    ! set the dependent/derived flags (MUSTDIE):
    call integralpar_setup()

    ! see the diff from defaults:
    DCALL integralpar_diff(task)
    DPRINT 'integralpar_set(',task,'): exit'
  end subroutine integralpar_set


  !*************************************************************
  subroutine integralpar_set_defaults()
    !  Purpose:  sets defaults
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of subroutines used ----------------
    !------------ Executable code --------------------------------
    ! integer:
    integralpar_i_int_part       = df_integralpar_i_int_part
    ! logical:
    integralpar_2cob_kin         = df_integralpar_2cob_kin
    integralpar_2cob_nuc         = df_integralpar_2cob_nuc
    integralpar_2cob_pvsp        = df_integralpar_2cob_pvsp
    integralpar_2cob_pvxp        = df_integralpar_2cob_pvxp
    integralpar_2cob_pvec        = df_integralpar_2cob_pvec
    integralpar_2cob_ol          = df_integralpar_2cob_ol
    integralpar_2cob_dipole      = df_integralpar_2cob_dipole
    integralpar_2cob_potential   = df_integralpar_2cob_potential
    integralpar_pot_for_secderiv = df_integralpar_pot_for_secderiv
    integralpar_2cob_field       = df_integralpar_2cob_field
#ifdef WITH_EFP
    integralpar_2cob_field_efp   = df_integralpar_2cob_field_efp
    integralpar_2cob_pc          = df_integralpar_2cob_pc
    integralpar_2cob_X           = df_integralpar_2cob_X
    integralpar_efp_gradients    = df_integralpar_efp_gradients
#endif
    integralpar_2cch_no          = df_integralpar_2cch_no
    integralpar_2cch_pre         = df_integralpar_2cch_pre
    integralpar_2cxc_no          = df_integralpar_2cxc_no
    integralpar_2c_mixed         = df_integralpar_2c_mixed
    integralpar_1cch_no          = df_integralpar_1cch_no
    integralpar_3c_xc            = df_integralpar_3c_xc
    integralpar_3c_co            = df_integralpar_3c_co
    integralpar_3c_co_resp       = df_integralpar_3c_co_resp
    integralpar_3c_rcoul_pvsp    = df_integralpar_3c_rcoul_pvsp
    integralpar_3c_rcoul_pvxp    = df_integralpar_3c_rcoul_pvxp
    integralpar_3c_r2_pvsp       = df_integralpar_3c_r2_pvsp
    integralpar_3c_r2_pvxp       = df_integralpar_3c_r2_pvxp
    integralpar_relativistic     = df_integralpar_relativistic
    integralpar_pseudo           = df_integralpar_pseudo !AS
    integralpar_spor             = df_integralpar_spor
    integralpar_gradients        = df_integralpar_gradients
    integralpar_2cob_ol_grad     = df_integralpar_2cob_ol_grad
    integralpar_3cob_grad        = df_integralpar_3cob_grad
    integralpar_send_3c          = df_integralpar_send_3c
    integralpar_send_2c          = df_integralpar_send_2c
    integralpar_offdiag_dipoles  = df_integralpar_offdiag_dipoles
    integralpar_gten_kinematic   = df_integralpar_gten_kinematic
    integralpar_2cff             = df_integralpar_2cff
    integralpar_2cob3c           = df_integralpar_2cob3c
    integralpar_totalsymmetric   = df_integralpar_totalsymmetric
    integralpar_rel_gradients    = df_integralpar_rel_gradients
    integralpar_2cob_kin_grad    = df_integralpar_2cob_kin_grad
    integralpar_2cob_ol_rel_grad = df_integralpar_2cob_ol_rel_grad
    integralpar_2cob_nuc_grad    = df_integralpar_2cob_nuc_grad
    integralpar_solv_grad        = df_integralpar_solv_grad
    integralpar_Q_solv_grad      = df_integralpar_Q_solv_grad
    integralpar_2cob_pvsp_grad   = df_integralpar_2cob_pvsp_grad
    ! these two should not be reset, because they are switched on/off
    ! from main_master independently of integral "runtype":
!   integralpar_cpksdervs        = df_integralpar_cpksdervs
!   integralpar_2dervs           = df_integralpar_2dervs
    integralpar_dervs            = df_integralpar_dervs
    integralpar_cpks_contribs    = df_integralpar_cpks_contribs
    integralpar_2cob_pc_grad     = df_integralpar_2cob_pc_grad
    integralpar_2cob_X_grad      = df_integralpar_2cob_X_grad
    integralpar_2cob_ipd_grad    = df_integralpar_2cob_ipd_grad
  end subroutine integralpar_set_defaults
  !*************************************************************

  subroutine integralpar_diff(runtype)
    !  Purpose:  sets defaults
    implicit none
    character(len=*), intent(in) :: runtype
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of subroutines used ----------------
    !------------ Executable code --------------------------------
#   define DIFFI(x,y) call diffi(x,y,__STRING(x))
#   define DIFFL(x,y) call diffl(x,y,__STRING(x))
    ! integer:
    DIFFI(integralpar_i_int_part       , df_integralpar_i_int_part )
    ! logical:
    DIFFL(integralpar_2cob_kin         , df_integralpar_2cob_kin )
    DIFFL(integralpar_2cob_nuc         , df_integralpar_2cob_nuc )
    DIFFL(integralpar_2cob_pvsp        , df_integralpar_2cob_pvsp )
    DIFFL(integralpar_2cob_pvxp        , df_integralpar_2cob_pvxp )
    DIFFL(integralpar_2cob_pvec        , df_integralpar_2cob_pvec )
    DIFFL(integralpar_2cob_ol          , df_integralpar_2cob_ol )
    DIFFL(integralpar_2cob_dipole      , df_integralpar_2cob_dipole )
    DIFFL(integralpar_2cob_potential   , df_integralpar_2cob_potential )
    DIFFL(integralpar_pot_for_secderiv , df_integralpar_pot_for_secderiv )
    DIFFL(integralpar_2cob_field       , df_integralpar_2cob_field )
#ifdef WITH_EFP
    DIFFL(integralpar_2cob_field_efp   , df_integralpar_2cob_field_efp )
    DIFFL(integralpar_2cob_pc          , df_integralpar_2cob_pc )
    DIFFL(integralpar_2cob_X           , df_integralpar_2cob_X )
    DIFFL(integralpar_efp_gradients    , df_integralpar_efp_gradients )
#endif
    DIFFL(integralpar_2cch_no          , df_integralpar_2cch_no )
    DIFFL(integralpar_2cch_pre         , df_integralpar_2cch_pre )
    DIFFL(integralpar_2cxc_no          , df_integralpar_2cxc_no )
    DIFFL(integralpar_2c_mixed         , df_integralpar_2c_mixed )
    DIFFL(integralpar_1cch_no          , df_integralpar_1cch_no )
    DIFFL(integralpar_3c_xc            , df_integralpar_3c_xc )
    DIFFL(integralpar_3c_co            , df_integralpar_3c_co )
    DIFFL(integralpar_3c_co_resp       , df_integralpar_3c_co_resp )
    DIFFL(integralpar_3c_rcoul_pvsp    , df_integralpar_3c_rcoul_pvsp )
    DIFFL(integralpar_3c_rcoul_pvxp    , df_integralpar_3c_rcoul_pvxp )
    DIFFL(integralpar_3c_r2_pvsp       , df_integralpar_3c_r2_pvsp )
    DIFFL(integralpar_3c_r2_pvxp       , df_integralpar_3c_r2_pvxp )
    DIFFL(integralpar_relativistic     , df_integralpar_relativistic )
    DIFFL(integralpar_pseudo           , df_integralpar_pseudo )
    DIFFL(integralpar_spor             , df_integralpar_spor )
    DIFFL(integralpar_gradients        , df_integralpar_gradients )
    DIFFL(integralpar_2cob_ol_grad     , df_integralpar_2cob_ol_grad )
    DIFFL(integralpar_3cob_grad        , df_integralpar_3cob_grad )
    DIFFL(integralpar_send_3c          , df_integralpar_send_3c )
    DIFFL(integralpar_send_2c          , df_integralpar_send_2c )
    DIFFL(integralpar_offdiag_dipoles  , df_integralpar_offdiag_dipoles )
    DIFFL(integralpar_gten_kinematic   , df_integralpar_gten_kinematic )
    DIFFL(integralpar_2cff             , df_integralpar_2cff )
    DIFFL(integralpar_2cob3c           , df_integralpar_2cob3c )
    DIFFL(integralpar_totalsymmetric   , df_integralpar_totalsymmetric )
    DIFFL(integralpar_rel_gradients    , df_integralpar_rel_gradients )
    DIFFL(integralpar_2cob_kin_grad    , df_integralpar_2cob_kin_grad )
    DIFFL(integralpar_2cob_ol_rel_grad , df_integralpar_2cob_ol_rel_grad )
    DIFFL(integralpar_2cob_nuc_grad    , df_integralpar_2cob_nuc_grad )
    DIFFL(integralpar_solv_grad        , df_integralpar_solv_grad )
    DIFFL(integralpar_Q_solv_grad      , df_integralpar_Q_solv_grad )
    DIFFL(integralpar_2cob_pvsp_grad   , df_integralpar_2cob_pvsp_grad )
    DIFFL(integralpar_cpksdervs        , df_integralpar_cpksdervs )
    DIFFL(integralpar_2dervs           , df_integralpar_2dervs )
    DIFFL(integralpar_dervs            , df_integralpar_dervs )
    DIFFL(integralpar_cpks_contribs    , df_integralpar_cpks_contribs )
    DIFFL(integralpar_2cob_pc_grad     , integralpar_2cob_pc_grad)
    DIFFL(integralpar_2cob_X_grad      , integralpar_2cob_X_grad)
    DIFFL(integralpar_2cob_ipd_grad    , integralpar_2cob_ipd_grad)
    contains

      subroutine diffl(op,df,nm)
        use error_module, only: MyID
        implicit none
        logical         , intent(in) :: op,df
        character(len=*), intent(in) :: nm
        ! *** end dof interface ***
        if( op .neqv. df )then
          print*,MyID//'integralpar_diff(',runtype,'): ',nm,'changed',df,'->',op
        endif
      end subroutine diffl

      subroutine diffi(op,df,nm)
        use error_module, only: MyID
        implicit none
        integer(i4_kind), intent(in) :: op,df
        character(len=*), intent(in) :: nm
        ! *** end dof interface ***
        if( op .ne. df )then
          print*,MyID//'integralpar_diff(',runtype,'): ',nm,'changed',df,'->',op
        endif
      end subroutine diffi
  end subroutine integralpar_diff

  !*************************************************************
  subroutine integralpar_setup()
    !  Purpose:  calculate derived variables depending upon values
    !  of other variables.
    !  called by: integralpar_set(task)
    implicit none
    !** End of interface *****************************************
    !------------ Executable code --------------------------------

    integralpar_2cff = integralpar_2cch_no .or. &
         integralpar_2cxc_no .or. &
         integralpar_2c_mixed .or. &
         integralpar_gradients
    integralpar_2cob3c = integralpar_2cob_kin .or. &
         integralpar_2cob_nuc .or. &
         integralpar_2cob_potential .or. &
         integralpar_2cob_field .or. &
#ifdef WITH_EFP
         integralpar_2cob_pc .or. &
         integralpar_2cob_X .or. &
         integralpar_efp_gradients .or. &
#endif
         integralpar_2cob_ol .or. &
         integralpar_3c_xc .or. &
         integralpar_3c_co .or. &
         integralpar_3c_co_resp .or. &
         integralpar_3c_rcoul_pvsp .or. &
         integralpar_3c_rcoul_pvxp .or. &
         integralpar_3c_r2_pvsp .or. &
         integralpar_3c_r2_pvxp .or. &
         integralpar_gradients

    integralpar_totalsymmetric = .not. (integralpar_3c_co_resp .or. integralpar_2cob_dipole)

    !!SB: PREVIOUS integralpar_totalsymmetric = .not. integralpar_2cob_dipole

    ! moved from main_integral (directly after call to integralpar_setup):
    if(integralpar_2cob_potential .or. &
       integralpar_2cob_field .or. &
#ifdef WITH_EFP
       integralpar_2cob_pc .or. &
       integralpar_2cob_X .or. &
#endif
       integralpar_solv_grad .or. &
       integralpar_Q_solv_grad) integralpar_2cff=.false.

       ! moved here from main_integral():
       if (integralpar_cpks_contribs) then
         integralpar_dervs = .false.
       else
         integralpar_dervs=integralpar_2dervs
       endif

       ! moved here from integral_setup():
       integralpar_dervs=integralpar_dervs.and.integralpar_gradients
  end subroutine integralpar_setup
  !*************************************************************

  !*************************************************************
  subroutine integralpar_set_field()
    !------------ Executable code --------------------------------
    integralpar_2cob_potential = .false.
    integralpar_pot_for_secderiv = .false.
    integralpar_2cob_field = .true.
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    ! 2-center orbital integral of electrostatic field
    integralpar_offdiag_dipoles = .false.
    ! no dipole matrix elements between diffrent Irreps
    integralpar_gradients = .false.
    ! no gradients
    integralpar_send_3c = .true.
    integralpar_send_2c = .false.
    ! integrals have to be sent
    integralpar_2cob_ol_grad = .false.
    ! gradients of overlap matrix
    integralpar_3cob_grad = .false.
    ! gradients of three center integrals
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.   ! options_relativistic
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.  !options_spin_orbit
    ! 2-center orbital integral of pvxp
#ifdef NEW_INTEGRALS
    integralpar_2cob_ol = .false.
#else
    integralpar_2cob_ol = .true.
#endif
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .false.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .false.
    ! 2-center charge norm integral of coulomb self-interaction
    !    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = .false.  !integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.  !options_relativistic
    ! perform relativistic transformations
    integralpar_pseudo = .false. !AS
    integralpar_spor = .false.  !options_spin_orbit
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do no calculation of relativistic gradients
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad=.false.
    integralpar_2cob_pvsp_grad=.false.
    ! gradients of pvsp
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 1
  end subroutine integralpar_set_field
  !*************************************************************



  !*************************************************************
  subroutine integralpar_set_potential()
    !------------ Executable code --------------------------------
    integralpar_2cob_potential = .true.
    ! 2-center orbital integral of potential
    integralpar_2cob_field = .false.
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    integralpar_offdiag_dipoles = .false.
    ! no dipole matrix elements between diffrent Irreps
    integralpar_gradients = .false.
    ! no gradients
    integralpar_send_3c = .true.
    integralpar_send_2c = .false.
    ! integrals have to be sent
    integralpar_2cob_ol_grad = .false.
    ! gradients of overlap matrix
    integralpar_3cob_grad = .false.
    ! gradients of three center integrals
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.   ! options_relativistic
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.  !options_spin_orbit
    ! 2-center orbital integral of pvxp
#ifdef NEW_INTEGRALS
    integralpar_2cob_ol = .false.
#else
    integralpar_2cob_ol = .true.
#endif
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .false.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .false.
    ! 2-center charge norm integral of coulomb self-interaction
    !    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = .false.  !integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.  !options_relativistic
    ! perform relativistic transformations
    integralpar_pseudo = .false. !AS
    integralpar_spor = .false.  !options_spin_orbit
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do no calculation of relativistic gradients
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad=.false.
    integralpar_2cob_pvsp_grad=.false.
    ! gradients of pvsp
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 1
  end subroutine integralpar_set_potential
  !*************************************************************


  !*************************************************************
  subroutine integralpar_set_scf()
    !  Purpose:  sets all variables contained in this module for
    !  run of integral part before scf part.
    !  called by main_master
    use operations_module, only: operations_response
    use spin_orbit_module, only: is_on, op_RelFit
    use pointcharge_module, only: present_X_centers !AS
    use options_module, only: options_relativistic, options_spin_orbit, &
         xcmode_extended_mda, xcmode_exchange_fit, options_xcmode
    !** End of interface *****************************************

    logical :: exchange_fit, extended_mda
    exchange_fit = options_xcmode() == xcmode_exchange_fit
    extended_mda = options_xcmode() == xcmode_extended_mda
    !------------ Executable code --------------------------------
    integralpar_2cob_pvsp = options_relativistic
    integralpar_2cob_pvxp = options_spin_orbit
    integralpar_2cob_pvec =  options_spin_orbit
    integralpar_2cch_pre = integralpar_resp_prescreen
    integralpar_2cxc_no = extended_mda
    integralpar_2c_mixed = extended_mda
    integralpar_3c_xc = exchange_fit .or. extended_mda

    integralpar_3c_co_resp = operations_response

    if( is_on(op_RelFit) ) then
       ! 3-center coulomb integral
       integralpar_3c_rcoul_pvsp = .true.
       ! 3-center scalar relativistic coulomb integral
       integralpar_3c_rcoul_pvxp = .true.
       ! 3-center vector relativistic coulomb integral
       integralpar_3c_r2_pvsp = .true.
       ! 3-center scalar relativistic coulomb integral
       integralpar_3c_r2_pvxp = .true.
       ! 3-center vector relativistic coulomb integral
    endif
    integralpar_relativistic = options_relativistic
    if(integralpar_relativistic .and. present_X_centers > 0) integralpar_pseudo=.true. !AS
    integralpar_spor = options_spin_orbit

    integralpar_i_int_part = 1
  end subroutine integralpar_set_scf
  !*************************************************************


  !*************************************************************
  subroutine integralpar_set_properties()
    !  Purpose:  sets all variables contained in this module for
    !  run of integral part before calculating properties.
    !  This integralpart is only run if normal; integralpart is not run.
    !  Only overlap and renorm coefficients are to be calculated
    !  called by main_master
    !** End of interface *****************************************
    !------------ Executable code --------------------------------
    integralpar_offdiag_dipoles = .false.
    ! no dipole matrix elements between diffrent Irreps
    integralpar_gradients = .false.
    ! no gradients
    integralpar_send_3c = .true.
    integralpar_send_2c = .true. ! why ??? (UB, 8/98)
    ! integrals have to be sent
    integralpar_2cob_ol_grad = .false.
    ! gradients of overlap matrix
    integralpar_3cob_grad = .false.
    ! gradients of three center integrals
    integralpar_2cob_potential = .false.
    integralpar_pot_for_secderiv = .false.
    ! 2-center orbital integral of potential
    integralpar_2cob_field = .false.
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.
    ! 2-center orbital integral of pvxp
    integralpar_2cob_pvec = .false. !<<<debug
    ! 2-center orbital integral of pvec
    integralpar_2cob_ol = .true.
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .false.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .false.
    ! 2-center charge norm integral of coulomb self-interaction
!    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.
    ! perform relativistic transformations
    integralpar_pseudo = .false. !AS
    integralpar_spor  = .false.
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do no calculation of relativistic gradients
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad=.false.
    integralpar_2cob_pvsp_grad=.false.
    ! gradients of pvsp
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 1
  end subroutine integralpar_set_properties
  !*************************************************************


  !*************************************************************
  subroutine integralpar_set_grad_solv()

    !------------ Executable code --------------------------------
    integralpar_offdiag_dipoles = .false.
    ! no dipole matrix elements between diffrent Irreps
    integralpar_gradients=.true.
    ! calculate gradients
    integralpar_send_3c=.false.
    integralpar_send_2c=.false.
    ! no integrals have to be sent
#ifdef NEW_INTEGRALS
    integralpar_2cob_ol_grad=.false.
#else
    integralpar_2cob_ol_grad=.true.
#endif
    ! gradients of overlap matrix
    integralpar_3cob_grad=.false.
    ! gradients of three center integrals
    integralpar_2cob_potential = .false.  !!!!!!!!!!!!!!!!!!!
    integralpar_pot_for_secderiv = .false.
    ! 2-center orbital integral of potential
    integralpar_2cob_field = .false.      !!!!!!!!!!!!!!!!!!!!!!
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.
    ! 2-center orbital integral of pvxp
    integralpar_2cob_ol = .false.
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .false.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .false.
    ! 2-center charge norm integral of coulomb self-interaction
!    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.
    ! perform relativistic transformations
    integralpar_spor = .false.
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do a calculation of relativistic gradients
    integralpar_pseudo = .false. !AS
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad= .true.
    integralpar_2cob_pvsp_grad=.false.
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 3
  end subroutine integralpar_set_grad_solv
  !*************************************************************


  !*************************************************************
  subroutine integralpar_set_gradient()
    !  Purpose:  sets all variables contained in this module for
    !  run of integral part in gradient part.
    !  called by main_master
    !** End of interface *****************************************
    !------------ Executable code --------------------------------
    integralpar_offdiag_dipoles = .false.
    ! no dipole matrix elements between diffrent Irreps
    integralpar_gradients=.true.
    ! calculate gradients
    integralpar_send_3c=.false.
    integralpar_send_2c=.false.
    ! no integrals have to be sent
    integralpar_2cob_ol_grad=.true.
    ! gradients of overlap matrix
    integralpar_3cob_grad=.true.
    ! gradients of three center integrals
    integralpar_2cob_potential = .false.
    integralpar_pot_for_secderiv = .false.
    ! 2-center orbital integral of potential
    integralpar_2cob_field = .false.
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.
    ! 2-center orbital integral of pvxp
    integralpar_2cob_pvec = .false.
    ! 2-center orbital integral of pvec
    integralpar_2cob_ol = .false.
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .false.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .true.
    ! 2-center charge norm integral of coulomb self-interaction
!    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.
    ! perform relativistic transformations
    integralpar_spor = .false.
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do a calculation of relativistic gradients
    integralpar_pseudo = .false. !AS
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad= .false.
    integralpar_2cob_pvsp_grad=.false.
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 3
  end subroutine integralpar_set_gradient
  !*************************************************************

  !*************************************************************
  subroutine integralpar_set_dipole(offdiag_dipoles)
    !  Purpose:  sets all variables contained in this module for
    !  run of dipole integral part.
    !  called by main_master
    !------------ Formal parameters ------------------------------
    logical, intent(in) :: offdiag_dipoles
    ! are dipole matrix elements between diffrent Irreps required ?
    !** End of interface *****************************************
    !------------ Executable code --------------------------------

    call integralpar_set('defaults')

    integralpar_offdiag_dipoles = offdiag_dipoles
    integralpar_gradients=.false.
    ! calculate gradients
    integralpar_send_3c=.false.
    integralpar_send_2c=.false.
    ! no integrals have to be sent
    integralpar_2cob_ol_grad=.false.
    ! gradients of overlap matrix
    integralpar_3cob_grad=.false.
    ! gradients of three center integrals
    integralpar_2cob_potential = .false.
    integralpar_pot_for_secderiv = .false.
    ! 2-center orbital integral of potential
    integralpar_2cob_field = .false.
#ifdef WITH_EFP
    integralpar_2cob_field_efp = .false.
    integralpar_2cob_pc = .false.
    integralpar_2cob_X = .false.
    integralpar_efp_gradients = .false.
#endif
    integralpar_2cob_kin = .false.
    ! 2-center orbital integral of kinetic energy
    integralpar_2cob_nuc = .false.
    ! 2-center orbital integral of nuclear attraction
    integralpar_2cob_pvsp = .false.
    ! 2-center orbital integral of pvsp
    integralpar_2cob_pvxp = .false.
    ! 2-center orbital integral of pvxp
    integralpar_2cob_pvec = .false.
    ! 2-center orbital integral of pvec
    integralpar_2cob_ol = .false.
    ! 2-center orbital integral of overlap
    integralpar_2cob_dipole = .true.
    ! 2-center orbital integral of dipole operator
    integralpar_2cch_no = .false.
    ! 2-center charge norm integral of coulomb self-interaction
!    integralpar_2cch_pre = integralpar_2cch_no
    integralpar_2cch_pre = integralpar_resp_prescreen
    ! estimate of 2-center charge overlap integral
    integralpar_2cxc_no = .false.
    ! 2-center exchange norm integral
    integralpar_2c_mixed = .false.
    ! mixed 2-center charge fit/exchange fit integral
    integralpar_1cch_no = .false.
    ! 1-center charge norm integral necessary for chargefit
    integralpar_3c_xc = .false.
    ! 3-center exchange integral
    integralpar_3c_co = .false.
    ! 3-center coulomb integral
    integralpar_3c_co_resp = .false.
    ! 3-center coulomb integral
    integralpar_3c_rcoul_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_rcoul_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_3c_r2_pvsp = .false.
    ! 3-center scalar relativistic coulomb integral
    integralpar_3c_r2_pvxp = .false.
    ! 3-center vector relativistic coulomb integral
    integralpar_relativistic = .false.
    ! perform relativistic transformations
    integralpar_spor = .false.
    ! perform spin orbit
    integralpar_rel_gradients=.false.
    ! do a calculation of relativistic gradients
    integralpar_pseudo = .false. !AS
    integralpar_2cob_nuc_grad=.false.
    ! gradients of nuclear attraction
    integralpar_solv_grad= .false.
    integralpar_2cob_pvsp_grad=.false.
    integralpar_2cob_kin_grad=.false.
    ! gradients of kinetic energy
    integralpar_2cob_ol_rel_grad=.false.

    integralpar_2cob_pc_grad     = .false.
    integralpar_2cob_X_grad      = .false.
    integralpar_2cob_ipd_grad    = .false.

    integralpar_i_int_part = 4
  end subroutine integralpar_set_dipole
  !*************************************************************


  subroutine integralpar_send_receive()
    use comm_module
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if(.not.comm_parallel() ) RETURN

    if(comm_i_am_master())then
       call comm_init_send(comm_all_other_hosts,msgtag_packed_message)
       call integralpar_pack()
       call comm_send()
    else
       call comm_save_recv(comm_master_host,msgtag_packed_message)
       call integralpar_unpack()
    endif
  end subroutine integralpar_send_receive


  subroutine integralpar_pack()
    ! purpose: packs data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: commpack
    use calc3c_switches
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call commpack(integralpar_2cob_potential,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_potential")
    call commpack(integralpar_pot_for_secderiv,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_pot_for_secderiv")
    call commpack(integralpar_2cob_field,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_field")
#ifdef WITH_EFP
    call commpack(integralpar_2cob_field_efp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_field_efp")
    call commpack(integralpar_2cob_pc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pc")
    call commpack(integralpar_2cob_X,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_X")
    call commpack(integralpar_efp_gradients,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_efp_gradients")
#endif
    call commpack(integralpar_2cob_kin,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_kin")
    call commpack(integralpar_2cob_nuc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_nuc")
    call commpack(integralpar_2cob_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pvsp")
    call commpack(integralpar_2cob_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pvxp")
    call commpack(integralpar_2cob_pvec,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pvec")
    call commpack(integralpar_2cob_ol,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_ol")
    call commpack(integralpar_2cob_dipole,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_dipole")
    call commpack(integralpar_2cch_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cch_no")
    call commpack(integralpar_2cch_pre,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cch_pre")
    call commpack(integralpar_2cxc_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cxc_no")
    call commpack(integralpar_2c_mixed,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2c_mixed")
    call commpack(integralpar_1cch_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_1cch_no")
    call commpack(integralpar_3c_xc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_xc")
    call commpack(integralpar_3c_co,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_co")
    call commpack(integralpar_3c_co_resp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_co_resp")
    call commpack(integralpar_3c_rcoul_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_rcoul_pvsp")
    call commpack(integralpar_3c_rcoul_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_rcoul_pvxp")
    call commpack(integralpar_3c_r2_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_r2_pvsp")
    call commpack(integralpar_3c_r2_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3c_r2_pvxp")
    call commpack(integralpar_relativistic,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_relativistic")
    call commpack(integralpar_pseudo,info)         !AS
    if (info .ne. 0) call error_handler( &         !AS
         "integralpar_pack: integralpar_pseudo")   !AS
    call commpack(integralpar_spor,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_spor")
    call commpack(integralpar_totalsymmetric,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_totalsymmetric")
    call commpack(integralpar_2cff,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cff")
    call commpack(integralpar_2cob3c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob3c")
    call commpack(integralpar_i_int_part,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_i_int_part")
    call commpack(integralpar_gradients,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_2cob_ol_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_ol_grad")
    call commpack(integralpar_3cob_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_3cob_grad")
    call commpack(integralpar_send_3c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_send_3c")
    call commpack(integralpar_send_2c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_send_2c")
    call commpack(integralpar_offdiag_dipoles,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_offdiag_dipoles")
    call commpack(integralpar_gten_kinematic,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_gten_kinematic")
    call commpack(integralpar_rel_gradients,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_rel_gradients")
    call commpack(integralpar_2cob_nuc_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_nuc_grad")
    call commpack(integralpar_2cob_pvsp_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pvsp_grad")
    call commpack(integralpar_2cob_kin_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_kin_grad")
    call commpack(integralpar_2cob_ol_rel_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_ol_rel_grad")
    call commpack(integralpar_solv_grad,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_Q_solv_grad,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_2cob_pc_grad,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_2cob_X_grad,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_2cob_ipd_grad,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_cpksdervs,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_2dervs,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_dervs,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_cpks_contribs,info)
    ASSERT(info.eq.0)
    call commpack(integralpar_runtype,info)
    ASSERT(info.eq.0)
  end subroutine integralpar_pack



  subroutine integralpar_unpack()
    ! purpose: unpacks data from comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: communpack
    use calc3c_switches
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    call communpack(integralpar_2cob_potential,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_potential")
    call communpack(integralpar_pot_for_secderiv,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_pot_for_secderiv")
    call communpack(integralpar_2cob_field,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_field")
#ifdef WITH_EFP
    call communpack(integralpar_2cob_field_efp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_field_efp")
    call communpack(integralpar_2cob_pc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_pc")
    call communpack(integralpar_2cob_X,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_2cob_X")
    call communpack(integralpar_efp_gradients,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_pack: integralpar_efp_gradients")
#endif
    call communpack(integralpar_2cob_kin,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_kin")
    call communpack(integralpar_2cob_nuc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_nuc")
    call communpack(integralpar_2cob_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_pvsp")
    call communpack(integralpar_2cob_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_pvxp")
    call communpack(integralpar_2cob_pvec,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_pvec")
    call communpack(integralpar_2cob_ol,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_ol")
    call communpack(integralpar_2cob_dipole,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_dipole")
    call communpack(integralpar_2cch_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cch_no")
    call communpack(integralpar_2cch_pre,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cch_pre")
    call communpack(integralpar_2cxc_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cxc_no")
    call communpack(integralpar_2c_mixed,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2c_mixed")
    call communpack(integralpar_1cch_no,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_1cch_no")
    call communpack(integralpar_3c_xc,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_xc")
    call communpack(integralpar_3c_co,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_co")
    call communpack(integralpar_3c_co_resp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_co_resp")
    call communpack(integralpar_3c_rcoul_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_rcoul_pvsp")
    call communpack(integralpar_3c_rcoul_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_rcoul_pvxp")
    call communpack(integralpar_3c_r2_pvsp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_r2_pvsp")
    call communpack(integralpar_3c_r2_pvxp,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3c_r2_pvxp")
    call communpack(integralpar_relativistic,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_relativistic")
    call communpack(integralpar_pseudo,info)     !AS
    if (info .ne. 0) call error_handler( &       !AS
         "integralpar_pack: integralpar_pseudo") !AS
    call communpack(integralpar_spor,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_spor")
    call communpack(integralpar_totalsymmetric,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_totalsymmetric")
    call communpack(integralpar_2cff,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cff")
    call communpack(integralpar_2cob3c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob3c")
    call communpack(integralpar_i_int_part,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_i_int_part")
    call communpack(integralpar_gradients,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_gradients")
    call communpack(integralpar_2cob_ol_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_ol_grad")
    call communpack(integralpar_3cob_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_3cob_grad")
    call communpack(integralpar_send_3c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_send_3c")
    call communpack(integralpar_send_2c,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_send_2c")
    call communpack(integralpar_offdiag_dipoles,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_offdiag_dipoles")
    call communpack(integralpar_gten_kinematic,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_gten_kinematic")
    call communpack(integralpar_rel_gradients,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_rel_gradients")
    call communpack(integralpar_2cob_nuc_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_nuc_grad")
    call communpack(integralpar_2cob_pvsp_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_pvsp_grad")
    call communpack(integralpar_2cob_kin_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_kin_grad")
    call communpack(integralpar_2cob_ol_rel_grad,info)
    if (info .ne. 0) call error_handler( &
         "integralpar_unpack: integralpar_2cob_ol_rel_grad")
    call communpack(integralpar_solv_grad,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_Q_solv_grad,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_2cob_pc_grad,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_2cob_X_grad,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_2cob_ipd_grad,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_cpksdervs,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_2dervs,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_dervs,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_cpks_contribs,info)
    ASSERT(info.eq.0)
    call communpack(integralpar_runtype,info)
    ASSERT(info.eq.0)
  end subroutine integralpar_unpack


  !--------------- End of module ----------------------------------
end module integralpar_module

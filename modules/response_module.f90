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
module response_module
  !-------------------------------------------------------------------
  !
  !  Purpose: This module generates the ground state data which
  !           is needed for calculation of dynamic linear response
  !           properties like
  !            - dynamical polarizabilities
  !            - excitation energies
  !           The data is written to several interface files which
  !           can be used with the "RESTD" response program.
  !
  !
  !  Module called by: main_master, main_slaves
  !
  !
  !  References: Notes A.G., H.H.
  !
  !
  !  Author: HH
  !  Date: 10/97
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification
  ! Author: HH
  ! Date:   11/98
  ! Description:
  ! - now correct 2-index-integral prescreening
  ! - reorganization of treatment of fractional occupation numbers
  ! - additional output of MO indices belonging to eigenvalue
  !   differences -> new version  response_write_eigenval_occ()
  ! - reorganization of the transformation of Coulomb integrals
  !   to MO basis in response_trans_Clb_integrals()
  !-------------------------------------------------------------------
  ! Modification
  ! Author: HH
  ! Date:   3/99
  ! Description:
  ! - Add additional debug output in 3-center-integrals
  !   now dump of rho,Vxc,fxc is possible -> only for spin=1
  ! - Add additional debug output in 2-/3-center-integrals
  !   now dump of <rho*fxc> and <vxc> possible.
  !   Note: rho & Vxc may be read in from file
  !         -> allows finite difference test on matrix elements
  !-------------------------------------------------------------------
  ! Modification
  ! Author: HH
  ! Date:   7/99
  ! Description: Adapted for use with comm_module of AS
  !              (AS did a similar thing in 7/98 on an outdated
  !               version of this module; however when the
  !               ttfs_hh_V20alpha-version was merged into V2.1
  !               it was easier to completely overwrite that version
  !               and to redo the changes)
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: HH
  ! Date:   10/99
  ! Description: Add treatment of hyperpolarizabilities
  !
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use operations_module, only: operations_dipole, operations_response
  use iounitadmin_module
  use readwriteblocked_module
  use filename_module
  use comm_module
  use msgtag_module
  use time_module    ! basic routines and datatype "timer_type" for timing
  use timer_module   ! variables of type(timer_type) for this program
  use input_module
  use output_module
  use machineparameters_module
  use eigen_data_module, only : eigval, eigen_data_bcast
  use occupied_levels_module
  use occupation_module, only: occ_num,n_occo
  use grid_module, only: more_grid_atom, atomicweight, grid_main, &
       grid_close, atomicweight_and_grad, weight_grads
  use grid_module, only: grid_loop_setup_atom
  use orbitalstore_module
  use orbital_module
  use density_data_module
  use density_calc_module, only: density_calc_setup,density_calc_nl
  use unique_atom_module
  use symmetry_data_module
  USE init_tddft_module, only: init_tddft_start
  USE tddft_diag, only: diag_init
  USE phys_param_module
  use resp_util_module
  use exchange
  use debug

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ public functions and subroutines ---------------------
  public :: response_read_input
  public :: response_write_input
  public :: response_input_bcast
  public :: response_main
  public :: response_do_dipole
  public :: response_setup
  public :: response_calc_2index_int_v2

  type, public :: two_index_type
     ! to contain calculated 2 index integrals for all Irreps
     real(kind=r8_kind), pointer :: o(:,:)
     ! o(dim_charge_matrix,n_spin*dim_factor)
     ! dim_charge_matrix = dimension_of_fit_ch(i_ir)*(dimension_of_fit_ch(i_ir)+1)/2
     ! n_spin*dim_factor = 2 or 4 or 8
  end type two_index_type

  integer(i4_kind), parameter :: & !! VERY IMPORTANT
       & X_NONE   = 0, &
       & C_NONE   = 0, &
       & C_VWN    = 1, &
       & C_PWlda  = 2, &
       & C_PBE    = 3, &
       & C_PW91   = 4, &
       & C_PERDEW = 5

  real(kind=r8_kind),      public :: rho_cutoff
  integer(kind = i4_kind), public :: X_KIND, C_KIND



  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------- Physical constants etc. ------------------------
  ! conversion factor from eV in hartree atomic units from NST
!  real(kind=r8_kind), parameter ::  ev2hartree =  0.036749309_r8_kind,&
!       & hartree2ev = 27.2113961_r8_kind


  ! namelist for ParaGau input file
  ! steering parameters for control of response calculation
  ! ----------- Default values for input parameters -------------
  logical               :: &
       & df_calc_all                   =.false.,&
       & df_lanczos                    =.false.,&
       & df_saved_XC                   =.false.,&
       & df_saved_2c_Q                 =.false.,&
       & df_saved_3c_Q                 =.false.,&
       & df_noRI                       =.false.,&
       & df_S_App                      =.false.,&
       & df_xalpha_resp                = .true.,&
       & df_vwn_resp                   =.false.,&

       & df_beckex_resp              =.false.,&       !!! SB B88    exchange
       & df_perdewc_resp               =.false.,&       !!! SB P      correlations

       & df_pw91x_resp                 =.false.,&       !!! SB PW91    exchange
       & df_pw91c_resp                 =.false.,&       !!! SB PW91    correlations
       & df_pbex_resp                  =.false.,&       !!! SB PBE     exchange
       & df_pbec_resp                  =.false.,&       !!! SB PBE     correlations

       & df_revpbex_resp               =.false.,&       !!! SB PBE     correlations
       & df_pbenx_resp                 =.false.,&       !!! SB PBE     correlations

       & df_pw_ldac_resp               =.false.,&       !!! SB PW_LDAC correlations
       & df_calc_osc_strength          =.false.,&
       & df_limit_unoccupied_levels    =.false.,&
       & df_NTO                        =.false.         !!!!MH NTO Calculations

  integer(kind=i4_kind) :: &
       & df_max_level_index            = -1_i4_kind,&
       & df_num_spectrum_levels        =  300_i4_kind, &
       & df_calc_n_low                 =  30_i4_kind,  &
       & df_max_iter                   =  30_i4_kind,  &
       & df_max_sp_trans               =  10_i4_kind

  real(kind=r8_kind)    :: &
       & df_LOWESTeV                   = 0.0_r8_kind, &
       & df_unoccupied_level_criterion = 1.0E-10_r8_kind,&
       & df_min_level_energy           =-1.0E+06_r8_kind,&
       & df_min_unocc_level_energy     =-1.0E+06_r8_kind,&
       & df_max_level_energy           = 1.0E+06_r8_kind,&
       & df_max_unocc_level_energy     = 1.0E+06_r8_kind,&
       & df_rho_cutoff                 = 1.0E-16_r8_kind,&
       & df_eigensolve_criterion       = 1.0E-16_r8_kind

  character(len=255)    :: df_hyperdatapath              ="<not used>"
  character(len=4)      :: df_target                     ="SS"
  ! ----------- Input parameter limits -------------------------
  real(kind=r8_kind) :: &
       & min_unoccupied_level_criterion = 1.0E-15_r8_kind,&
       & max_unoccupied_level_criterion = 1.0E-03_r8_kind
  ! ----------- Response control parameters --------------------
  logical               :: &
       & calc_all, &
       & lanczos,  &
       & saved_XC, &
       & saved_2c_Q, &
       & saved_3c_Q, &
       & noRI, &
       & S_App, &
       & xalpha_resp,&
       & vwn_resp,&
       & beckex_resp, &
       & perdewc_resp, &
       & pw91c_resp,&
       & pw91x_resp,&
       & pbec_resp, &
       & pbex_resp, &
       & revpbex_resp, &
       & pbenx_resp, &
       & pw_ldac_resp, &
       & calc_osc_strength,&
       & limit_unoccupied_levels, &
       & NTO                     !!!!MH NTO Calculations

  character(len=255)    :: hyperdatapath
  character(len=4)      :: target

  integer(kind=i4_kind) :: &
       & calc_n_low, &
       & max_sp_trans, &
       & max_iter, &
       & max_level_index,&
       & num_spectrum_levels

  real(kind=r8_kind)    :: &
       & LOWESTeV,&
       & min_level_energy,&
       & max_level_energy,&
       & min_unocc_level_energy,&
       & max_unocc_level_energy,&
       & eigensolve_criterion, &
       & unoccupied_level_criterion

  namelist /response_control/ &
       & target, &
       & xalpha_resp,&
       & vwn_resp,&
       & beckex_resp, &
       & perdewc_resp, &
       & pw91x_resp,&
       & pw91c_resp,&
       & pbec_resp, &
       & pbex_resp, &
       & revpbex_resp, &
       & pbenx_resp, &
       & pw_ldac_resp, &
       & calc_all, &
       & lanczos,  &
       & saved_XC, &
       & saved_2c_Q, &
       & saved_3c_Q, &
       & noRI, &
       & S_App, &
       & LOWESTeV, &
       & calc_n_low, &
       & max_sp_trans, &
       & max_iter, &
       & calc_osc_strength,&
       & rho_cutoff,&
       & unoccupied_level_criterion,&
       & limit_unoccupied_levels,&
       & max_level_index,&
       & min_level_energy,&
       & max_level_energy,&
       & min_unocc_level_energy, &
       & max_unocc_level_energy, &
       & num_spectrum_levels,&
       & eigensolve_criterion, &
       & NTO                           !!!!! MH NTO calculations

  ! names of the interface files for the response program which
  ! are written in this module:
  character(len=30),parameter ::&
       & header_file         = 'resp_header.out'      ,&
       & r2index_file        = 'xc_2c_'      ,&
       & r3index_file        = 'resp_3index.dat'      ,&
       & Clb_integral_file   = 'resp_coulomb.dat'     ,&
       & eigenoccind_file    = 'resp_eig_occ_ind.dat' ,&
       & dipole_file         = 'resp_dip'

  ! miscellaneous (set in response_setup):
  ! - vector length for this machine,
  ! - factor for dimension of response matrix elements
  integer(kind=i4_kind) :: vec_length, dim_factor

  ! num-of-:chargefit-functions,spins,irreps
  integer(kind=i4_kind) :: n_spin, n_irrep

  integer(kind=i4_kind), ALLOCATABLE :: n_chargefit(:) !! depend on irreps

  ! charge overlap matrix for prescreening(packed storage) as calculated
  ! in modules - integral_2c_fit_ch_module  (calculation)
  !            - int_send_2cff_module       (saving to tape)
  ! and its dimension, which will also be used for 2-index-integrals
  ! see module "mat_charge_module.f90" for description of storage

  real(kind=r8_kind),allocatable,dimension(:) :: charge_prematrix
  integer(kind=i4_kind),allocatable :: dim_charge_matrix(:) !! i_irrep

!!$  ! electronic density rho, energy density fxc and derivative of
!!$  ! V_xc with respect to rho
!!$  real(kind=r8_kind),allocatable  :: rho(:,:), fxc(:), dvdrho(:,:)

  ! orbitals and chargefit functions
  !  -> orbital_spin_type%o(:,:,:,:) ! (only for the response part needed)
  !     to contain calculated sym. orbitals of one irrep irreps(i)
  !     o(vec_len,symmetry_data_dimension(i),symmetry_data_n_partners(i),n_spin)
  type(orbital_spin_type),pointer :: phi_ob(:)   ! for 3 index integrals only
  type(orbital_type),pointer      :: orbs_ob(:)  ! for 2/3 index integrals
  type(orbital_gradient_type),pointer :: orbs_grads(:)
! type(orbital_type)              :: orbs_ch     ! chargefit function for all integrals
! real(kind=r8_kind),pointer      :: orbs_ch_p(:,:) ! pointer on component of orbs_ch

  ! gridpoints and weigths
  real(kind=r8_kind),pointer      :: grdpts(:,:),grdwts(:)

  ! num. 2 index integrals <f_k|dvdrho(:,:)|f_l> in packed storage mode
  real(kind=r8_kind),allocatable  :: r2_index_matrix(:,:)

  ! SB: second version for all symmetry
  ! SB: num. 2 index integrals <f_k|dvdrho(:,:)|f_l> in packed storage mode
  ! SB: added one more index for irreps
  type(two_index_type),pointer      :: r2_index_matrix_v2(:)

  ! num. 3 index integrals = array of arrays
  ! super matrix "r3_index_matrix(n_irrep,n_spin,dim_factor)%m"
  ! sub  matrix "m(n_is,max_n_chargefit_per_processor)"
  ! note:
  !   r3_index_matrix(a,sigma,tau)%m(as,k) = <phi_a phi_s| fxc(sigma,tau)|f_k>
  !   where
  !        fxc(sigma,tau) = d V_xc(sigma) / d rho(tau)
  !        sigma,tau      = alpha/beta-spin components
  !
  type(arrmat2),pointer   :: r3_index_matrix(:,:,:)

  real(r8_kind) :: max_level_energy_au, min_level_energy_au ! save!
  real(r8_kind) :: max_unocc_level_energy_au, min_unocc_level_energy_au ! save!

  ! contains analytical Coulomb integrals C(is,k)=<phi_i phi_s|f_k>
  ! for fixed spin and irrep and chargefit function k
  real(kind=r8_kind),pointer :: coulomb_matrix(:)

  !*********************************************************************
  !      ----**** for HYPERPOLARIZABILITY calculation ****----

  ! vector of expansion coefficients:
  !         rho1phi(1:2,1:xyz=3,1:Nas) into phi_a(r)*phi_s(r)
  !         rho1fit(1:2,1:xyz=3,1:Nk)  into g_k(r)
  ! for three components xyz of external perturbation (dipole moment operator)
  ! and two different frequencies omega1 and omega2 of ext. perturbation
  !
  ! the ...fit...  variables contain the expansion coefficients into the
  ! fit basis. They are needed for the Coulomb-contribution to the
  ! matrix F below [ rho1fitvec(1:N_fit,1:xyz=3,1:2) ]
  real(kind=r8_kind),allocatable :: rho1phi(:,:,:),rho1fit(:,:,:)

  ! the perturbing frequencies omega1,omega2
  real(kind=r8_kind)             :: omega(2)

  ! the differences of eigenvalues and the transition dipole moments
  !   ediff(1:N,1:N)
  !   transdip(1:3=x,y,z,1:N,1:N)
  ! where N = # of MOs,
  real(kind=r8_kind),allocatable :: resp_ediff(:,:),resp_transdip(:,:,:)

  ! total dipole moment in au
  real(kind=r8_kind) ::  resp_tot_dipmom(3)



  ! the matrices F(1:N,1:N,1:3,1:2) and G(1:3=xyz,1:3=x'y'z',1:Nas)
  ! where N = # of MOs,
  ! and 1:2 refers two the two different incident frequencies
  ! for performance reasons oemga 1 belongs to m1
  ! and omega2 to m2, and we use them as in G(m2,m1,:)
  real(kind=r8_kind),allocatable :: resp_Gmat(:,:,:),resp_Fmat(:,:,:,:)

  ! the vector for the constant contribution
  !  \int dr r_k*h(\omega1,\omega2,r_i,r_j;r)
  ! to the hyperpolarizability tensor beta_ijk
  ! resp_bsum(1:3,1:3,1:3,1:2) -> beta(i,j,k,+/-)
  !
  ! bsum(1:2,m3,m2,m1)  -> m1 belongs to omega1
  !                        m2 belongs to omega2
  !                        m3 belongs to component of observable=dipole operator
  !                        1:2 BECAUSE OF OMEGA1+OMEGA2 (=1), OMEGA1-OMEGA2 (=2)
  real(kind=r8_kind) ::  resp_bsum(2,3,3,3)

  ! the numerical 3-index integrals and the 3-index Coulomb integrals
  ! r3indall_matrix(1:N,1:N,1:Nfit)
  ! r3clball_matrix(1:N,1:N,1:Nfit)
  real(kind=r8_kind),allocatable :: r3indall_matrix(:,:,:),&
       & r3clball_matrix(:,:,:)


  ! the final RHS h(1:3,1:3,1:Nfit) as calculated in the last step
  real(kind=r8_kind),allocatable :: resp_hpluss(:,:,:),resp_hminus(:,:,:)


  !!FOR main_response
  FPP_TIMER_DECL(RESPONSE_ALL)
  FPP_TIMER_DECL(RESPONSE_SETUP)
  FPP_TIMER_DECL(RESPONSE_DIPOL)
  FPP_TIMER_DECL(COULOMB_2C)
  FPP_TIMER_DECL(XC_2C)
  FPP_TIMER_DECL(COULOMB_3C)
  FPP_TIMER_DECL(DIAG_PLUS_DVDSON)
  FPP_TIMER_DECL(RESPONSE_CLOSE)

  !! FOR local XC_2C
  FPP_TIMER_DECL(rho_calc)
  FPP_TIMER_DECL(XC_calc)
  FPP_TIMER_DECL(khi_calc)
  FPP_TIMER_DECL(loop_calc)
  FPP_TIMER_DECL(all)
  !!#endif

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  !**************************************************************
  subroutine response_read_input()
    ! Purpose : Reading the input namelist "response_control"
    !           which contains logicals which are used to determine
    !           what kind of interface data for the response
    !           program shall be calculated.
    !** End of interface **************************************
    !------------ Declaration of local variables ---------------------
    USE global_module, ONLY: gl_eig_crite
    integer(kind=i4_kind) :: unit,status !!$,iounit
!!$    character(len=255)    :: pathname
    !------------ Executable code ------------------------------------

    ! default values
    if (n_spin==2) then
       target                  = "open"
    else
       target                  = trim(adjustl(df_target))
    end if

    X_KIND = X_NONE
    C_KIND = C_NONE

    xalpha_resp                = df_xalpha_resp
    vwn_resp                   = df_vwn_resp

    beckex_resp              = df_beckex_resp     !!!! SB implementation of B88     exchange
    perdewc_resp               = df_perdewc_resp      !!!! SB implementation of PERDEW  correlation

    pw91x_resp                 = df_pw91x_resp         !!!! SB implementation of PW91    exchange
    pw91c_resp                 = df_pw91c_resp         !!!! SB implementation of PW91    correlation
    pbex_resp                  = df_pbex_resp          !!!! SB implementation of PBE     exchange
    pbec_resp                  = df_pbec_resp          !!!! SB implementation of PBE     correlation

    revpbex_resp               = df_revpbex_resp
    pbenx_resp                 = df_pbenx_resp

    pw_ldac_resp               = df_pw_ldac_resp       !!!! SB implementation of PW_LDAC correlation

    calc_all                   = df_calc_all
    lanczos                    = df_lanczos
    saved_XC                   = df_saved_XC
    saved_2c_Q                 = df_saved_2c_Q
    saved_3c_Q                 = df_saved_3c_Q
    noRI                       = df_noRI
    S_App                      = df_S_App
    calc_n_low                 = df_calc_n_low
    LOWESTeV                   = df_LOWESTeV
    max_sp_trans               = df_max_sp_trans
    max_iter                   = df_max_iter

    calc_osc_strength          = df_calc_osc_strength
    rho_cutoff                 = df_rho_cutoff
    unoccupied_level_criterion = df_unoccupied_level_criterion
    limit_unoccupied_levels    = df_limit_unoccupied_levels
    max_level_index            = df_max_level_index
    min_level_energy           = df_min_level_energy
    max_level_energy           = df_max_level_energy
    min_unocc_level_energy     = df_min_unocc_level_energy
    max_unocc_level_energy     = df_max_unocc_level_energy
    num_spectrum_levels        = df_num_spectrum_levels
    eigensolve_criterion       = df_eigensolve_criterion

    NTO                        = df_NTO                     !!!!!! MH NTO calculations

    if ( input_line_is_namelist("response_control") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=response_control, iostat=status)
       if (status .gt. 0) call input_error( &
            "response_read_input: namelist response_control")
    endif

    if ( calc_osc_strength .and. operations_response .and.&
         & (.not. operations_dipole) )&
         & call input_error( &
         "response_read_input: namelist response_control: for response&
         & part with dipole_moments = .true. it is necessary that &
         & calc_osc_strength = .true. as well")

    ! check if cut off criterion for the 2/3-index-integrals makes sense
    if( eigensolve_criterion < 0.0_r8_kind) then
       call write_to_output_units("response_read_input: warning: &
            &eigensolve_criterion < 0 ")
       eigensolve_criterion = abs(eigensolve_criterion)
       call write_to_output_units("                     &
            &eigensolve_criterion set to ",re=eigensolve_criterion)
    end if

    if( eigensolve_criterion < 1.0E-50_r8_kind) then
       call write_to_output_units("response_read_input: warning: &
            &eigensolve_criterion <  1.0E-50")
       eigensolve_criterion = 1.0E-50_r8_kind
       call write_to_output_units("                     &
            &eigensolve_criterion set to 1.0E-50")
    end if

    gl_eig_crite = eigensolve_criterion


    ! check if cut off criterion for the 2/3-index-integrals makes sense
    if( rho_cutoff < 0.0_r8_kind) then
       call write_to_output_units("response_read_input: warning: &
            &rho_cutoff < 0 ")
       rho_cutoff = abs(rho_cutoff)
       call write_to_output_units("                     &
            &rho_cutoff set to ",re=rho_cutoff)
    end if

    if( rho_cutoff < 1.0E-50_r8_kind) then
       call write_to_output_units("response_read_input: warning: &
            &rho_cutoff <  1.0E-50")
       rho_cutoff = 1.0E-50_r8_kind
       call write_to_output_units("                     &
            &rho_cutoff set to 1.0E-50")
    end if

    if( unoccupied_level_criterion<min_unoccupied_level_criterion) then
       call write_to_output_units("response_read_input: warning: &
            &unoccupied_level_criterion<",re=min_unoccupied_level_criterion)
       call write_to_output_units("                     &
            &unoccupied_level_criterion set to ",re=min_unoccupied_level_criterion)
       unoccupied_level_criterion=min_unoccupied_level_criterion
    end if
    if( unoccupied_level_criterion>max_unoccupied_level_criterion) then
       call write_to_output_units("response_read_input: warning: &
            &unoccupied_level_criterion>",re=max_unoccupied_level_criterion)
       call write_to_output_units("                     &
            &unoccupied_level_criterion set to ",re=max_unoccupied_level_criterion)
       unoccupied_level_criterion=max_unoccupied_level_criterion
    end if

    ! if using a cutoff for unoccupied orbitals is requested, i.e.
    ! "limit_unoccupied_levels=.true." then check which of the two
    ! mutually exclusive criteria were specified
    if(limit_unoccupied_levels) then
       if(max_level_index /= -1_i4_kind) then

          if(max_level_index <2_i4_kind) call input_error( &
               & "response_read_input: namelist response_control: for &
               & limit_unoccupied_levels=.true. a value of&
               & max_level_index > 1 is required !")

       else

          max_level_energy_au = max_level_energy * ev2hartree
          min_level_energy_au = min_level_energy * ev2hartree
          max_unocc_level_energy_au = max_unocc_level_energy * ev2hartree
          min_unocc_level_energy_au = min_unocc_level_energy * ev2hartree
       end if
    end if

    ! check the number of MO levels to be printed in the spectrum
    if(num_spectrum_levels<0) then
       call input_error( &
            "response_read_input: namelist response_control: num_spectrum_levels>0 &
            &is required !")
    end if

  end subroutine response_read_input
  !**************************************************************


  !**************************************************************
  subroutine response_write_input(iounit)
    ! purpose : Writing the whole namelist "response_control" to
    !           the output
    !
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(i4_kind),intent(in) :: iounit
    !** End of interface *****************************************

    word_format = '("    ",a," = ",a  :" # ",a)'

    call start("RESPONSE_CONTROL","RESPONSE_WRITE_INPUT", &
         iounit,operations_echo_input_level)

    call word("TARGET                     ", target                     ,&
         & df_target )
    call flag("XALPHA_RESP                ", xalpha_resp                ,&
         & df_xalpha_resp )
    call flag("VWN_RESP                   ", vwn_resp                   ,&
         & df_vwn_resp )
    call flag("BECKEX_RESP                ", beckex_resp                ,&
         & df_beckex_resp )
    call flag("PERDEWC_RESP               ", perdewc_resp               ,&
         & df_perdewc_resp )
    call flag("PW91X_RESP                 ", pw91x_resp                 ,&
         & df_pw91x_resp )
    call flag("PW91C_RESP                 ", pw91c_resp                 ,&
         & df_pw91c_resp )
    call flag("PBEX_RESP                  ", pbex_resp                  ,&
         & df_pbex_resp )
    call flag("PBEC_RESP                  ", pbec_resp                  ,&
         & df_pbec_resp )
    call flag("REVPBEX_RESP               ", revpbex_resp               ,&
         & df_revpbex_resp )
    call flag("PBENX_RESP                 ", pbenx_resp                 ,&
         & df_pbenx_resp )
    call flag("PW_LDAC_RESP               ", pw_ldac_resp               ,&
         & df_pw_ldac_resp )
    call flag("CALC_ALL                   ", calc_all                   ,&
         & df_calc_all )

    call flag("LANCZOS                    ", lanczos                    ,&
         & df_lanczos )

    call flag("SAVED_XC                   ", saved_XC                   ,&
         & df_saved_XC )
    call flag("SAVED_2C_Q                 ", saved_2c_Q                 ,&
         & df_saved_2c_Q )
    call flag("SAVED_3C_Q                 ", saved_3c_Q                 ,&
         & df_saved_3c_Q )

    call flag("noRI                       ", noRI                       ,&
         & df_noRI )
    call flag("S_APP                      ", S_App                      ,&
         & df_S_App )
    call intg("CALC_N_LOW                 ", calc_n_low                 ,&
         & df_calc_n_low)
    call real("LOWESTeV                   ", LOWESTeV                   ,&
         & df_LOWESTeV)
    call intg("MAX_SP_TRANS               ", max_sp_trans               ,&
         & df_max_sp_trans)
    call intg("MAX_ITER                   ", max_iter                   ,&
         & df_max_iter)
    call flag("CALC_OSC_STRENGTH          ", calc_osc_strength          ,&
         & df_calc_osc_strength )
    call real("RHO_CUTOFF                 ", rho_cutoff                 ,&
         & df_rho_cutoff )
    call real("EIGENSOLVE_CRITERION       ", eigensolve_criterion       ,&
         & df_eigensolve_criterion)
    call real("UNOCCUPIED_LEVEL_CRITERION ", unoccupied_level_criterion ,&
         & df_unoccupied_level_criterion)
    call flag("LIMIT_UNOCCUPIED_LEVELS    ", limit_unoccupied_levels    ,&
         & df_limit_unoccupied_levels)
    call intg("MAX_LEVEL_INDEX            ", max_level_index            ,&
         & df_max_level_index)
    call real("MIN_LEVEL_ENERGY           ", min_level_energy           ,&
         & df_min_level_energy)
    call real("MAX_LEVEL_ENERGY           ", max_level_energy           ,&
         & df_max_level_energy)
    call real("MIN_UNOCC_LEVEL_ENERGY     ", min_unocc_level_energy     ,&
         & df_min_unocc_level_energy)
    call real("MAX_UNOCC_LEVEL_ENERGY     ", max_unocc_level_energy     ,&
         & df_max_unocc_level_energy)

    call intg("NUM_SPECTRUM_LEVELS        ", num_spectrum_levels        ,&
         & df_num_spectrum_levels)
    call flag("NTO                        ", NTO                        ,&      !!!!!MH NTO calculations
         & df_NTO)

    call stop()

  end subroutine response_write_input


  !*************************************************************
  subroutine response_main()
    !
    ! Main routine of this module. Executed by all workers. FIXME: may
    ! need  some  work,  as  it  was previousely  executed  by  master
    ! only. Controls post-SCF calculation  of data needed for response
    ! calculations.
    !
    USE int_send_2c_resp, only: int_send_2c_resp_rewrite
    USE int_resp_module, only: int_resp_Clb_3c
    USE resp_dipole_module
    USE noRI_module, only: noRI_2c
    use comm, only: comm_rank
    implicit none
    !** End of interface *****************************************

#ifdef FPP_TIMERS
    FPP_TIMER_ZERO(RESPONSE_ALL)
    FPP_TIMER_ZERO(RESPONSE_SETUP)
    FPP_TIMER_ZERO(RESPONSE_DIPOL)
    FPP_TIMER_ZERO(COULOMB_2C)
    FPP_TIMER_ZERO(XC_2C)
    FPP_TIMER_ZERO(COULOMB_3C)
    FPP_TIMER_ZERO(DIAG_PLUS_DVDSON)
    FPP_TIMER_ZERO(RESPONSE_CLOSE)
#endif
    FPP_TIMER_START(RESPONSE_ALL)
    FPP_TIMER_START(RESPONSE_SETUP)

    call response_setup()

    if (.true.) then

       call write_to_trace_unit("Entering response part")

       call start_timer(timer_resp)

       ! *** Setup and Preparations ***
       call start_timer(timer_resp_preparations)

       call stop_timer(timer_resp_preparations)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% do calculate the integrals for LINEAR response only if %%%
       !%%% NO hyperpolarizability calculation is requested        %%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$       if(.NOT. (beta_A2 .or. beta_C) ) then

       ! write header-interface-file with information for response program
       call write_to_output_units("response_main: response_write_header")
       call start_timer(timer_resp_preparations)
       call start_timer(timer_resp_header)
       if (comm_rank() == 0) then
          call response_write_header()
       endif
       call stop_timer(timer_resp_header)

       ! write MO eigenvalues and occupation numbers to a tape
       call write_to_output_units("response_main: response_write_eigenval_occupation")
       if (comm_rank() == 0) then
          call response_write_eigenval_occ()
       endif
       call stop_timer(timer_resp_preparations)

       FPP_TIMER_STOP (RESPONSE_SETUP)
       FPP_TIMER_START(RESPONSE_DIPOL)

       ! *** Dipoles ***
       ! do we have to produce the tape with transition dipole moments ?
       if (calc_osc_strength) then
          call write_to_trace_unit  ("Rewrite tapes with transition dipole moments")
          call write_to_output_units("response_main: response_rewrite_dipoletape")
          call start_timer(timer_resp_dipole)
          if (comm_rank() == 0) then
             call resp_dipole_rewrite() ! serial
          endif
          call stop_timer(timer_resp_dipole)
       end if

       FPP_TIMER_STOP(RESPONSE_DIPOL)

       ! *** numerical 2-index-integrals ***
       ! user input: do we have to calculate the numerical 2 index integrals ?

       FPP_TIMER_START(RESPONSE_SETUP)

       call write_to_trace_unit("Calculate numerical 2 index integrals")
       call start_timer(timer_resp_2index)

       ! setup the grid and divide into parts and send those to slaves
       call write_to_output_units("response_main: grid_main")
       call grid_main (post_scf=.true.)

       FPP_TIMER_STOP (RESPONSE_SETUP)
       FPP_TIMER_START(COULOMB_2C)

       ! read <f_k|1/(r-r')|f_l> chargefit overlap integrals,
       ! and rewrite them to response interface file
       if (.not. saved_2c_Q) then
          call write_to_output_units("response_main: int_send_2c_resp_rewrite")
          call int_send_2c_resp_rewrite()
       end if

       FPP_TIMER_STOP (COULOMB_2C)
       FPP_TIMER_START(XC_2C)

       if (.not. saved_XC) then !! .and. (.not. noRI)) then
          call write_to_output_units("response_main: response_calc_2index_int_v2")
          call response_calc_2index_int_v2()
!!          call noRI_2c(X_KIND,C_KIND)
       end if

       call stop_timer(timer_resp_2index)

       ! Clean up stuff that was needed for numerical 2/3-index integrals
       call write_to_output_units("response_main: shutdown orbital_module")
       ! shutdown orbital_module
       call orbital_free(orbs_ob)
       call orbital_shutdown()

       FPP_TIMER_STOP (XC_2C)
       FPP_TIMER_START(COULOMB_3C)

       ! *** Transformation of Coulomb integrals ***
       ! user input: do we have to calculate the analytical 3-index integrals
       !             <phi_i phi_s|V_Clb|f_k> ?
       if (.not. saved_3c_Q) then
          call write_to_trace_unit("Transform 3 index Coulomb integrals")
          call start_timer(timer_resp_Coulomb)

          call write_to_output_units("response_main: int_resp_Clb_3c")

          call int_resp_Clb_3c()
          call stop_timer(timer_resp_Coulomb)
       end if

       ! SB: THE MAIN CALCULATION IS HERE
       call write_to_output_units("response_main: main response calculations (START)")

       FPP_TIMER_STOP (COULOMB_3C)
       FPP_TIMER_START(DIAG_PLUS_DVDSON)

       call init_tddft_start()

       call diag_init()

       call write_to_output_units &
            ("response_main: main response calculations (FINISH)")


       FPP_TIMER_STOP (DIAG_PLUS_DVDSON)
       FPP_TIMER_START(RESPONSE_CLOSE)

       ! clean up: deallocate leftovers (module private objects etc...)
       call write_to_output_units("response_main: shutdown response_module")
       call response_shutdown()

       ! do the timing (we always do timing since the response part will
       ! be the most time consuming part of the whole calculation)
       call write_to_output_units("response_main: call timer_print_response")
       call stop_timer(timer_resp)

       ! Shutdown the grid which was distributed among processors
       call write_to_output_units("response_main: grid_close")
       call grid_close (.false.)

       ! now we are finished
       call write_to_trace_unit  ("Exit response part")
       call write_to_output_units("response_main: done")
       FPP_TIMER_STOP(RESPONSE_CLOSE)
       FPP_TIMER_STOP(RESPONSE_ALL)

#ifdef FPP_TIMERS
       ! TIMER OUTPUT
       block
          real (r8_kind) :: tt
          tt = FPP_TIMER_VALUE(RESPONSE_SETUP)
          WRITE (*,*) "[m] RESPONSE TIMER "
          WRITE (*,*) "   * SETUP      ", tt
          tt = FPP_TIMER_VALUE(RESPONSE_DIPOL)
          WRITE (*,*) "   * DIPOLE     ", tt
          tt = FPP_TIMER_VALUE(COULOMB_2C)
          WRITE (*,*) "   * COULOMB 2C ", tt
          tt = FPP_TIMER_VALUE(XC_2C)
          WRITE (*,*) "   * XC 2C      ", tt
          tt = FPP_TIMER_VALUE(COULOMB_3C)
          WRITE (*,*) "   * COULOMB 3C ", tt
          tt = FPP_TIMER_VALUE(DIAG_PLUS_DVDSON)
          WRITE (*,*) "   * DIAG/DVDS  ", tt
          tt = FPP_TIMER_VALUE(RESPONSE_CLOSE)
          WRITE (*,*) "   * CLOSE      ", tt
          tt = FPP_TIMER_VALUE(RESPONSE_ALL)
          WRITE (*,*) "   * SUMMARY    ", tt
       end block
#endif

    end if !!i_am_master_

  end subroutine response_main
  !*************************************************************

  !*****************************************************************************
  subroutine response_input_bcast()
    ! Purpose : broadcasting the response_input
    ! called by subroutine send_recv_init_options
    use comm, only: comm_bcast
    implicit none
    !** End of interface *******************************************************

    call comm_bcast(target)
    call comm_bcast(vwn_resp)
    call comm_bcast(xalpha_resp)
    call comm_bcast(beckex_resp)
    call comm_bcast(perdewc_resp)
    call comm_bcast(pw91x_resp)
    call comm_bcast(pw91c_resp)
    call comm_bcast(pbex_resp)
    call comm_bcast(pbec_resp)
    call comm_bcast(pbenx_resp)
    call comm_bcast(revpbex_resp)
    call comm_bcast(pw_ldac_resp)
    call comm_bcast(calc_all)
    call comm_bcast(lanczos)
    call comm_bcast(saved_XC)
    call comm_bcast(saved_2c_Q)
    call comm_bcast(saved_3c_Q)
    call comm_bcast(calc_n_low)
    call comm_bcast(LOWESTeV)
    call comm_bcast(max_sp_trans)
    call comm_bcast(max_iter)
    call comm_bcast(calc_osc_strength)
    call comm_bcast(rho_cutoff)
    call comm_bcast(eigensolve_criterion)
    call comm_bcast(unoccupied_level_criterion)
    call comm_bcast(limit_unoccupied_levels)
    call comm_bcast(max_level_index)
    call comm_bcast(min_level_energy)
    call comm_bcast(max_level_energy)
    call comm_bcast(min_unocc_level_energy)
    call comm_bcast(max_unocc_level_energy)
  end subroutine response_input_bcast
  !*****************************************************************************

  !*************************************************************
  subroutine response_setup()
    !
    ! For  master  and every  slave:  Setup  module private  variables
    ! needed in several subroutines.
    !
    use ch_response_module, only: dimension_of_fit_ch
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: i_ir, status

    ! Get suitable vectorlength for this machine
    vec_length = machineparameters_veclen

    ! get symmetry information (from symmetry_data_module)
    n_spin  = ssym%n_spin      ! number of spins
    n_irrep = ssym%n_irrep     ! number of irreps

    ALLOCATE(n_chargefit(n_irrep),STAT=status)
    ASSERT(status==0)

    ! get the number of chargefit functions (from fit_coeff_module)
    do i_ir = 1, n_irrep
       n_chargefit(i_ir) = dimension_of_fit_ch(i_ir)
    end do

    ! user input: do also calculate integrals wrto abs-norms ?
    ! then we need larger dimensions in the arrays for the
    ! matrix elements etc.
    dim_factor=2

    if (comm_i_am_master()) then
       if(output_response_detailed) then
          call write_to_output_units(&
               & "response_setup: set local variables:" )
          call write_to_output_units(&
               & "    vec_length        = ",vec_length)
          call write_to_output_units(&
               & "    n_spin            = ",n_spin)
          call write_to_output_units(&
               & "    n_irrep           = ",n_irrep)

          do i_ir = 1, n_irrep
             call write_to_output_units(&
                  & "    i_irrep       = ", i_ir)
             call write_to_output_units(&
                  & "       n_chargefit       = ", n_chargefit(i_ir))
          end do

          call write_to_output_units(&
               & "    dim_factor        = ",dim_factor)
       end if
    end if

    ! user input: only if numerical 2 index integrals are to be
    ! calculated we need this information
    ! allocate memory for the chargefit overlap integrals
    ! overlap matrix is symmetric
    allocate (dim_charge_matrix(n_irrep),STAT=status)
    ASSERT(status==0)
    do i_ir = 1, n_irrep
       dim_charge_matrix(i_ir) = (n_chargefit(i_ir)*(n_chargefit(i_ir)+1))/2
    end do

    ! for numerical integrals we need the density and orbitals
    ! set up density calculation
    call write_to_output_units("response_setup: set up density calculation ->&
         & call density_calc_setup()")
    call density_calc_setup()

    ! allocation of orbitals
    call write_to_output_units(&
         & "response_setup: set up orbital calculations ->&
         & call orbital_setup(), orbital_allocate()")
    call orbital_setup   (vec_length,.true.,.true.)
    call orbital_allocate (orbs_ob, orbs_grads, phi_ob=phi_ob) !,orb_ch=orbs_ch)


    ! broadcast eigenvectos to slaves:
    call eigen_data_bcast()

    if (comm_i_am_master())then
       call resp_util_set_level_index(&
            & unoccupied_level_criterion,&
            & limit_unoccupied_levels,&
            & max_level_index,&
            & min_level_energy_au,&
            & max_level_energy_au, &
            & min_unocc_level_energy_au,&
            & max_unocc_level_energy_au, &
            & num_spectrum_levels)
    else
       call resp_util_upck_level_index()
    end if
  end subroutine response_setup
  !*************************************************************

  !*************************************************************
  subroutine response_write_header()
    !
    ! Write  header  interface  file  with  information  for  response
    ! program.  Called by master through response_main().
    !
    use echo_input_module
    use ch_response_module,   only: dimension_of_fit_ch
    use symmetry_data_module, only: symmetry_data_n_irreps
    use clebsch_gordan,       only: cg=>cg_eliminated, prod_bas
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: io_unit, io_stat
    integer(kind=i4_kind) :: i_spin, i_irrep, i_ir_a, i_ir_b, i_ir_c
    integer(kind=i4_kind) :: occ_times_unocc_dim, nmult, nc, ou_dim
    integer(kind=i4_kind) :: i_mlt
    character(6)          :: Xw,Cw

    io_unit = openget_iounit (status='unknown', &
         file= trim (resp_dir) // '/' // trim (header_file))

    n_irrep = symmetry_data_n_irreps()

    ! First write  number of spins, irreps and  chargefit functions as
    ! namelist  "system_data".  Set  echo level  to  suppress printing
    ! defaults.  Most values have no meaningful defaults even.
    call start ("SYSTEM_DATA","RESPONSE_WRITE_HEADER", io_unit, &
         echo_level_full)
    call intg ("NUM_SPINS         ", n_spin, -1)
    call intg ("NUM_IRREPS        ", n_irrep, -1)
    ! FIXME: this  is not a  standard namelist.  You cannot  have more
    ! than one entry  with the same name there.  This precludes use of
    ! standard namelist IO on the reader side!
    do i_irrep = 1, n_irrep
       nc = dimension_of_fit_ch(i_irrep)
       call intg ("NUM_CHARGEFITFCTS ", nc, -1)
    end do
    call stop()

    ! Now  write the  input options  for the  response  calculation as
    ! namelist
    call start ("RESPONSE_CONTROL","RESPONSE_WRITE_HEADER", io_unit, &
         echo_level_full)

    if (n_spin == 2) then
       call word ("TARGET            ", "SS", df_target)
    else
       call word ("TARGET            ", trim (adjustl (target)), df_target)
    end if

    if (xalpha_resp ) X_KIND = X_XALPHA
    if (vwn_resp    ) C_KIND = C_VWN
    if (beckex_resp ) X_KIND = X_BECKE88
    if (perdewc_resp  ) C_KIND = C_PERDEW
    if (pw91x_resp  ) X_KIND = X_PW91
    if (pw91c_resp  ) C_KIND = C_PW91
    if (pbex_resp   ) X_KIND = X_PBE
    if (pbec_resp   ) C_KIND = C_PBE
    if (revpbex_resp   ) X_KIND = X_REVPBE
    if (pbenx_resp  ) X_KIND = X_PBEN
    if (pw_ldac_resp) C_KIND = C_PWlda

    select case (X_KIND)
    case (X_XALPHA)
       Xw = "XA"
    case (X_BECKE88)
       Xw = "BECKE"
    case (X_PBE)
       Xw = "PBE"
    case (X_REVPBE)
       Xw = "REVPBE"
    case (X_PBEN)
       Xw = "PBEN"
    case (X_PW91)
       Xw = "PW91"
    case default
       Xw = "NONE"
    end select

    select case (C_KIND)
    case (C_VWN)
       Cw = "VWN"
    case (C_PERDEW)
       Cw = "PERDEW"
    case (C_PBE)
       Cw = "PBE"
    case (C_PW91)
       Cw = "PW91"
    case (C_PWlda)
       Cw = "PWlda"
    case default
       Cw = "NONE"
    end select

    call strng("EXCHANGE          ", TRIM (Xw), "NONE")
    call strng("CORRELATION       ", TRIM (Cw), "NONE")

    call intg ("MAX_ITER          ",  max_iter, df_max_iter)
    call flag ("CALC_ALL          ",  calc_all, df_calc_all)
    call flag ("LANCZOS           ",  lanczos, df_lanczos)
    call flag ("noRI              ",  noRI, df_noRI)
    call flag ("S_APP             ",  S_App, df_S_App)
    call intg ("CALC_N_LOW        ",  calc_n_low, df_calc_n_low)
    call real ("LOWESTeV          ",  LOWESTeV, df_LOWESTeV)
    call intg ("MAX_SP_TRANS      ",  max_sp_trans, df_max_sp_trans)

    ! MH: NTO calculations
    call flag ("NTO               ",  NTO, df_NTO)

    if (calc_osc_strength) then
       call flag("CALC_OSC_STR      ", .true., .false.)
    else
       call flag("CALC_OSC_STR      ", .false., .false.)
    end if
    call stop()

    ! Write the dimension (occupied * unoccupied) for each irrep/spin
    write (io_unit, '("#   IR C  SPIN  NUMBER OF TRANSITIONS")', iostat=io_stat)
    ASSERT(io_stat==0)
    write (output_unit, '("#   IR C  SPIN  NUMBER OF TRANSITIONS")', iostat=io_stat)
    ASSERT(io_stat==0)

    i_spin_: do i_spin=1,n_spin
       i_ir_c_: do i_ir_c = 1, n_irrep
          ou_dim=0_i4_kind
          i_ir_a_: do i_ir_a = 1,n_irrep
             i_ir_b_: do i_ir_b = 1,n_irrep
                nmult = cg(i_ir_c,i_ir_a,i_ir_b)%mult
                i_mlt_: do i_mlt = 1, nmult
                   call resp_util_calc_transitions_v2(i_ir_a, i_ir_b, i_ir_c, &
                        i_spin, occ_times_unocc_dim)
                   ou_dim = ou_dim + occ_times_unocc_dim
                end do i_mlt_
             end do i_ir_b_
          end do i_ir_a_
          write (io_unit, '(4X,I2,5X,I2,5X,I5)', iostat=io_stat) &
               i_ir_c, i_spin, ou_dim
          ASSERT(io_stat==0)
          write (output_unit, '(4X,I2,5X,I2,5X,I5)', iostat=io_stat) &
               i_ir_c, i_spin, ou_dim
          ASSERT(io_stat==0)
       end do i_ir_c_
    end do i_spin_

    call return_iounit(io_unit)
  end subroutine response_write_header
  !*************************************************************

  !*************************************************************
  subroutine response_2index_receive_v2 (n_irrep, n_spin)
    !  Purpose: receive the results from the slaves
    !  Values of 2-index-integrals on the grid portion of
    !  accessable to each slave are collected from the slaves
    !  and added up to yield the final matrix
    !------------ Modules used ------------------- ---------------
    use ch_response_module, only: dimension_of_fit_ch
    use xpack, only: upck
    implicit none
    integer(i4_kind) :: n_irrep, n_spin
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: alloc_stat,n_procs,i_proc
    integer(kind=i4_kind) :: n_dim,i_ir,dcm
    real(kind=r8_kind),allocatable :: help_mat(:,:)
    !------------ Executable code ------------------------------------

    n_procs=comm_get_n_processors()
    do i_proc=2,n_procs

       do i_ir=1,n_irrep

          n_dim = dimension_of_fit_ch(i_ir)
          dcm = (n_dim*(n_dim+1)/2)
          allocate(help_mat(dcm,n_spin*2),stat=alloc_stat)
          ASSERT(alloc_stat==0)
          help_mat  = 0.0_r8_kind

          call comm_save_recv(i_proc,msgtag_response_2in_send)
          if(comm_msgtag()/=msgtag_response_2in_send) &
               call error_handler('Wrong msgtag in response_2index_receive')

          call upck(help_mat)

          r2_index_matrix_v2(i_ir)%o(:,:)=&
               & r2_index_matrix_v2(i_ir)%o + help_mat

          deallocate(help_mat,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do

    end do

  end subroutine response_2index_receive_v2
  !*************************************************************

  !*************************************************************
  subroutine response_2index_tomaster_v2 (n_irrep)
    !  Purpose: send results from the slave to the master
    !  Values of 2-index-integrals on the grid portion of
    !  accessable to each slave are send to the master
    !------------ Modules used ------------------- ---------------
    use ch_response_module, only: dimension_of_fit_ch
    use xpack, only: pck
    implicit none
    integer(i4_kind), intent(in) :: n_irrep
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i_irrep
    !------------ Executable code ------------------------------------

    do i_irrep = 1,n_irrep
       call comm_init_send(comm_master_host,msgtag_response_2in_send)
       call pck(r2_index_matrix_v2(i_irrep)%o(:,:))
       call comm_send()
    end do

  end subroutine response_2index_tomaster_v2
  !*************************************************************

  !*************************************************************
  subroutine response_write_2index_totape_v2(n_irrep,n_spin,dim_factor)
    !  Purpose:
    !  Write 2index matrix row by row to a linear tape.
    !------------ Modules used ------------------- ---------------
!   use iounitadmin_module,   only: openget_iounit,returnclose_iounit
    use io, only: write_buffer
    implicit none
    integer(kind=i4_kind), intent(in) :: n_irrep,n_spin,dim_factor
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)             :: i_spin,i_dim,i_ir,idx
!   integer(kind=i4_kind)             :: io_unit
    !------------ Executable code ------------------------------------

    do i_ir = 1, n_irrep
       do i_spin=1, n_spin
          do i_dim=1, dim_factor
             idx = (i_spin-1)*dim_factor+i_dim
!            io_unit = openget_iounit(trim(resp_dir)//'/'&
!                 //resp_util_fname('xc_2c',i_ir,idx),form='unformatted',status='unknown')
             call write_buffer( trim(resp_dir)//'/'//resp_util_fname('xc_2c',i_ir,idx) &
                              , r2_index_matrix_v2(i_ir)%o(:,idx)                      &
                              )
!            call returnclose_iounit(io_unit)
          end do
       end do
    end do

  end subroutine response_write_2index_totape_v2
  !*************************************************************

  !*************************************************************
  subroutine response_write_eigenval_occ()
    use symmetry_data_module, only: symmetry_data_n_spin, symmetry_data_n_irreps
    use clebsch_gordan, only: cg=>cg_eliminated
    use resp_util_module
    ! Purpose:
    ! Write occupation numbers and eigenvalues to an interface
    ! tape for the response program.
    ! Note that in the response program only DIFFERENCES
    !    "(epsilon_i - epsilon_s)"  and  "(f_i - f_s)"
    ! of eigenvalues "epsilon_i,s" and (possibly fractional)
    ! occupation numbers "f_i,s" between occupied (index "i")
    ! and unoccupied (index "s") levels are needed,
    ! like e.g. in the
    !
    !                     2*(epsilon_i - epsilon_s)*(f_i - f_s)
    ! lambda(is;omega) = ---------------------------------------
    !                      (epsilon_i - epsilon_s)^2 - omega^2
    !
    ! quantity, where "omega" is the frequency.
    ! However, the ordering of the difference-meta-index "is"
    ! is crucial: it has to be the same as for the numerical
    ! and analytical 3 index integrals calculated elsewhere in
    ! this module !
    !------------ Modules used ------------------- ---------------
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)   :: io_unit, io_stat
    integer(kind=i4_kind)   :: i_ir_a, i_ir_b, i_spin, i_occ, i_unocc, i_ir_c
    integer(kind=i4_kind)   :: is_counter
    integer(kind=i4_kind)   :: occ_times_unocc_dim=0_i4_kind
    integer(kind=i4_kind)   :: occs, occe, unoccs, unocce
    real(kind=r8_kind)      :: eig_diff, occ_diff
    integer(kind=i4_kind)   :: na, nb, npaa, npab, nmult, i_mlt
    !------------ Executable code ------------------------------------

    n_spin  = symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()

    i_spin_: do i_spin=1, n_spin

       i_ir_c_: do i_ir_c = 1, n_irrep

          ! open the interfacefiles
          io_unit = openget_iounit(trim(resp_dir)//'/'//resp_util_fname('eg_oc',i_ir_c,i_spin),&
               & form='unformatted',status='unknown')

          is_counter = 0_i4_kind   ! combined metaindex a->b

          i_ir_a_: do i_ir_a = 1, n_irrep
             npaa = symmetry_data_n_partners(i_ir_a)
             na   = ssym%dim(i_ir_a)
             i_ir_b_: do i_ir_b = 1, n_irrep
                npab = symmetry_data_n_partners(i_ir_b)
                nb   = ssym%dim(i_ir_b)
                nmult = cg(i_ir_c,i_ir_a,i_ir_b)%mult
                i_mlt_: do i_mlt = 1, nmult

                   ! get the number of possible transitions
                   call resp_util_calc_transitions_v2(i_ir_a, i_ir_b, i_ir_c, &
                        i_spin, occ_times_unocc_dim)

                   if (occ_times_unocc_dim == 0) then
                      DPRINT 'No transitions... skip'
                      cycle !! check it out
                   end if

                   call resp_util_borders(i_ir_a,i_ir_b,i_spin,occs,occe,unoccs,unocce)

                   do i_occ= occs, occe
                      do i_unocc=unoccs, unocce
                         if (abs(eigval(i_ir_a)%m(i_occ,i_spin)-eigval(i_ir_b)%m(i_unocc,i_spin))<min_diff) cycle

                         ! compute differences
                         eig_diff =  eigval(i_ir_a)%m(i_occ,i_spin)      -  eigval(i_ir_b)%m(i_unocc,i_spin)
                         occ_diff = occ_num(i_ir_a)%m(i_occ,i_spin)/npaa - occ_num(i_ir_b)%m(i_unocc,i_spin)/npab

                         ! write computed differences to the future interface file
                         is_counter=is_counter + 1
                         write(unit=io_unit,iostat=io_stat) is_counter, i_ir_a, i_ir_b, i_occ, i_unocc, eig_diff, occ_diff
                         ASSERT(io_stat==0)
                      end do
                   end do
                end do i_mlt_
             end do i_ir_b_
          end do i_ir_a_
          ! close the interface file
          call returnclose_iounit(io_unit)
       end do i_ir_c_
    end do i_spin_

  end subroutine response_write_eigenval_occ
  !*************************************************************

  !*************************************************************
  subroutine response_shutdown
    !  Purpose: Clean up (deallocate private objects etc.)
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: alloc_stat
    !------------ Executable code ------------------------------------


    ! todo:
    ! call on every slave
    ! free_eigvec -> eigen_data module falls num. 2/3-index integrale

    ! deallocate other stuff in other modules ? grid after 3index ?
    ! bounds ?

    ! information about fractionally occupied levels
    deallocate(MO_status, begin_index, end_index,stat=alloc_stat)
    if(alloc_stat/=0) call error_handler(&
         & "response_shutdown: deallocation 1 failed")

  end subroutine response_shutdown
  !*************************************************************

  !*************************************************************
  logical function response_do_dipole()
    ! Purpose: returns value of the module private logical variable
    ! "dipole_moments", which is from the input namelist
    ! "response_control".
    implicit none
    !** End of interface ***************************************

    response_do_dipole = calc_osc_strength

  end function response_do_dipole
  !*************************************************************

  !*************************************************************
  subroutine response_calc_2index_int_v3()
    !  Purpose:
    ! calculate 2 index integrals <f_k| h_v |f_k>
    ! + parallelization over grid (like in post_scf_module)
    ! + integral prescreening using <f_k|f_l>
    !------------ Modules used ------------------- ---------------
    use vwnc
    use pw_ldac_module
    use gga_response_module
    use ch_response_module
    use exchange
    use constants
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: vla, &
         & i_atom,i_fit,i_fit2,&
         & i_vec, counter, i_dim, i_spin
    integer(i4_kind) :: alloc_stat, i_ir, n_pa, i_pa, n_ir, n_dim

!!! GGA variables

    real(kind=r8_kind), dimension(vec_length,2)  :: rho !!$, fxc(:), dvdrho(:,:) !!n_spin
    real(kind=r8_kind), dimension(vec_length,n_spin) :: nn

    real(kind=r8_kind), dimension(vec_length)   :: &
         Fx, Fc
    real(kind=r8_kind), dimension(vec_length,2) :: &
         dFxdn, dFcdn
    real(kind=r8_kind), dimension(vec_length,3) :: &
         dFxdg, dFcdg, dFxdndn, dFcdndn
    real(kind=r8_kind), dimension(vec_length,4) :: &
         dFdndn, dFdg
    real(kind=r8_kind), dimension(vec_length,6) :: &
         dFxdndg, dFcdndg
    real(kind=r8_kind), dimension(vec_length,6)   :: &
         dFxdgdg, dFcdgdg

    real(kind=r8_kind) :: &
         grarho(vec_length,3,2),gamma(vec_length,3)

    type(orbital_type),pointer             :: fcts_orbs_ch(:)
    type(orbital_gradient_type),pointer    :: grads_ch(:)

    real (r8_kind), pointer :: orbs_ch_p(:,:) ! pointer on component of orbs_ch
    real (r8_kind), pointer :: grad_orbs_ch_p(:,:,:)

    integer(i4_kind), parameter ::&
         A = 1, &
         B = 2

    integer(i4_kind), parameter ::&
         & XX = 1, &
         & YY = 2, &
         & ZZ = 3

    integer(i4_kind), parameter ::&
         & AA = 1, &
         & AB = 2, &
         & BB = 3, &
         & BA = 4

    integer(i4_kind), parameter ::&
         & AAA = 1, &
         & BBB = 2, &
         & AAB = 3, &
         & ABB = 4, &
         & BAB = 5, &
         & BAA = 6

    integer(i4_kind), parameter ::&
         & AAAA = 1, &
         & BBBB = 2, &
         & ABAB = 3, &
         & AABB = 4, &
         & AAAB = 5, &
         & BBAB = 6

    real(r8_kind) :: int

    !! For debug
    real(kind=r8_kind) :: &
         Flda(vec_length), dFldadn(vec_length,2), dFldadndn(vec_length,3)

    integer(i4_kind) :: CC

    !------------ Executable code ------------------------------------

    if (( .not. xalpha_resp ) .and. &
         ( .not. vwn_resp    ) .and. &
         ( .not. beckex_resp  ) .and. &
         ( .not. perdewc_resp  ) .and. &
         ( .not. pw91x_resp  ) .and. &
         ( .not. pw91c_resp  ) .and. &
         ( .not. pbex_resp   ) .and. &
         ( .not. revpbex_resp   ) .and. &
         ( .not. pbenx_resp   ) .and. &
         ( .not. pbec_resp   ) .and. &
         ( .not. pw_ldac_resp)) return

    print *,"WARNING: NO BLAS"
    print *," rho_cutoff = ",rho_cutoff

    n_ir = ssym%n_irrep

    ! allocate memory for 2-index-integrals
    allocate(r2_index_matrix_v2(n_ir),&
         & stat=alloc_stat)
    if(alloc_stat/=0) call error_handler(&
         & "response_calc_2index_int_v2: allocation r2_index_matrix_v2 failed")

    do i_ir = 1,n_ir
       n_dim = dimension_of_fit_ch(i_ir)
!!$       dim_charge_matrix = (n_dim*(n_dim+1)/2)
       allocate(r2_index_matrix_v2(i_ir)%o(dim_charge_matrix(i_ir),n_spin*dim_factor),&
            & stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            & "response_calc_2index_int_v2: allocation r2_index_matrix_v2%o failed")
       r2_index_matrix_v2(i_ir)%o(:,:) = zero
    enddo


    call fit_fct_allocate_response(fcts_ch = fcts_orbs_ch, grads = grads_ch)

    vla = machineparameters_veclen

    call orbital_setup_response(vla,.true.) !! do Gradients !!

    ! now do integration over the grid
    call write_to_trace_unit(&
         & "response_calc_2index_int_v2: calculating data for")

    !
    ! Loop over grid points:
    !
    call grid_loop_setup_atom()

    grid_points_: do while( more_grid_atom(vec_length, i_atom, grdpts, grdwts) )
       ASSERT(i_atom/=0)

       vla=size(grdpts,1)

       ! calculate gridweights
       grdwts(1:vla)=grdwts(1:vla)*&
            atomicweight(i_atom, grdpts(1:vla, :))*&
            unique_atoms(i_atom)%n_equal_atoms

       ! calculate orbitals and charge fitfunctions on the grid
       ! (orbitals are needed for the density below)
       call orbital_calculate(grdpts(1:vla,1:3),&
            vla,orbs_ob,orbs_grads)

       rho    = zero
       gamma  = zero
       grarho = zero
       call density_calc_nl(vla,nn,gamma,grarho,orbs_ob,orbs_grads)

       if (n_spin==1) then
          rho(:,1) = nn(:,1)/TWO
          rho(:,2) = nn(:,1)/TWO

          gamma(:,2) = gamma(:,1)/FOUR
          gamma(:,3) = gamma(:,1)/FOUR
          gamma(:,1) = gamma(:,1)/FOUR
       else
          rho(:,1) = nn(:,1)
          rho(:,2) = nn(:,2)
       end if


       !! CUT OFF
       where ( abs(rho) < rho_cutoff )
          rho = rho_cutoff
       endwhere

       !! EXCHANGE PART
       Flda      = zero
       dFldadn   = zero
       dFldadndn = zero

       Fx      = zero
       dFxdn   = zero
       dFxdg   = zero
       dFxdndn = zero
       dFxdndg = zero
       dFxdgdg = zero
       if ( xalpha_resp )  &
            call exchange_lda(X_XALPHA,vla,2,rho,Flda,dFldadn,dFxdndn)
       if ( pbex_resp   )  then
          call exchange_lda(X_XALPHA,vla,2,rho,Flda,dFldadn,dFldadndn)
          call exchange_gga(X_PBE   ,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
          Fx      = Fx      + Flda
          dFxdn   = dFxdn   + dFldadn
          dFxdndn = dFxdndn + dFldadndn
       end if

       !! CORRELATION PART
       Flda      = zero
       dFldadn   = zero
       dFldadndn = zero

       Fc      = zero
       dFcdn   = zero
       dFcdg   = zero
       dFcdndn = zero
       dFcdndg = zero
       dFcdgdg = zero
       if ( vwn_resp    )  &
            call vwn_resp_calc(rho,n_spin,vla,rho_cutoff,dFcdndn) !! FIXME: SHOULD BE REWRITTEN!!
       if ( pw_ldac_resp )  &
            call pw_ldac   (vla,4,rho,Fc,dFcdn,dFcdndn) !! 4 = unrestr case
       if ( pbec_resp   )  then
          call pw_ldac   (vla,4,rho,Flda,dFldadn,dFldadndn)
          call gga_correlation(C_PBE,2,2,rho,gamma,vla,&
               Fc,Flda,&
               dFcdn,dFcdg,dFldadn,&
               dFcdndn,dFcdndg,dFcdgdg,dFldadndn)
          Fc      = Fc      + Flda
          dFcdn   = dFcdn   + dFldadn
          dFcdndn = dFcdndn + dFldadndn
       end if

       if (n_spin == 1) then
          CC = 2
       else
          CC = 1
       end if

       dFxdg   = (dFxdg  +dFcdg  ) / CC
       dFxdndn = (dFxdndn+dFcdndn) / CC
       dFxdndg = (dFxdndg+dFcdndg) / CC
       dFxdgdg = (dFxdgdg+dFcdgdg) / CC

#ifdef WITH_ISNAN
       if(countNaN(dFxdg).ne.0)then
          print *,'response: ', countNaN(dFxdg),' NaNs in dFxdg'
          stop
       end if
       if(countNaN(dFxdndn).ne.0)then
          print *,'response: ', countNaN(dFxdndn),' NaNs in dFxdndn'
          stop
       end if
       if(countNaN(dFxdndg).ne.0)then
          print *,'response: ', countNaN(dFxdndg),' NaNs in dFxdndg'
          stop
       end if
       if(countNaN(dFxdgdg).ne.0)then
          print *,'response: ', countNaN(dFxdgdg),' NaNs in dFxdgdg'
          stop
       end if
#endif

       dFdg   = zero
       dFdndn = zero

       select case (n_spin)
       case (1)
          dFdg  (:,AA) = dFxdg(:,1)
          dFdg  (:,AB) = dFxdg(:,3) / TWO !! will be x2 later

          dFdndn(:,AA) = dFxdndn(:,AA)
          dFdndn(:,AB) = dFxdndn(:,AB)
       case (2)
          dFdg  (:,AA) = dFxdg(:,1)
          dFdg  (:,AB) = dFxdg(:,3) / TWO !! will be x2 later
          dFdg  (:,BB) = dFxdg(:,2)
          dFdg  (:,BA) = dFxdg(:,3) / TWO !! will be x2 later

          dFdndn(:,AA) = dFxdndn(:,AA)
          dFdndn(:,AB) = dFxdndn(:,AB)
          dFdndn(:,BB) = dFxdndn(:,BB)
          dFdndn(:,BA) = dFxdndn(:,AB)
       end select

       i_ir1_: do i_ir = 1,n_ir
          fcts_orbs_ch(i_ir)%o = zero
          grads_ch(i_ir)%o     = zero
       enddo i_ir1_

       call fit_fct_calculate_response(grdpts(1:vla,1:3),vla, &
            fcts_orbs_ch, grads_ch )

       call start_timer(timer_resp_2index_integration)

       i_ir_: do i_ir = 1,n_ir !!symmetry_data_n_irreps()
          n_pa = symmetry_data_n_partners(i_ir)
          i_pa_: do i_pa = 1, n_pa !!symmetry_data_n_irreps()

             orbs_ch_p          => fcts_orbs_ch(i_ir)%o(:,:,i_pa)
             grad_orbs_ch_p     => grads_ch(i_ir)%o(:,:,:,i_pa)

             i_spin_: DO i_spin = 1,n_spin

                i_dim_:  do i_dim=1,dim_factor

                   ! now calculate matrix elements <f_i|dvdrho(:,j)|f_k>, j=1..8
                   counter=0   ! combined index for packed storage of charge_prematrix
                   fit: do i_fit=1, dimension_of_fit_ch(i_ir) ! dimension of ch fit fct
                      fit2: do i_fit2=1,i_fit
                         counter = counter + 1
                         int = zero
                         i_vec_: do i_vec = 1, vla

                            !! dfxdndn (LDA)
                            int = int + &
                                 grdwts(i_vec) &
                                 * orbs_ch_p(i_vec,i_fit ) &
                                 * orbs_ch_p(i_vec,i_fit2) &
                                 * dFdndn(i_vec,(i_spin-1)*dim_factor+i_dim)

                            !! dfx_dgm
                            int = int + &
                                 grdwts(i_vec) * &
                                 ( &
                                 grad_orbs_ch_p(i_vec,1,i_fit) * grad_orbs_ch_p(i_vec,1,i_fit2) + &
                                 grad_orbs_ch_p(i_vec,2,i_fit) * grad_orbs_ch_p(i_vec,2,i_fit2) + &
                                 grad_orbs_ch_p(i_vec,3,i_fit) * grad_orbs_ch_p(i_vec,3,i_fit2) &
                                 ) &
                                 * TWO*dFdg(i_vec,(i_spin-1)*dim_factor+i_dim)

                            select case (n_spin)
                            case (1) !! closed shell
                               select case ((i_spin-1)*dim_factor+i_dim)
                               case (AA)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * (dFxdndg(i_vec,AAA)+dFxdndg(i_vec,AAB) / TWO)

                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * &
                                       (   dFxdgdg(i_vec,AAAA) &
                                       & + dFxdgdg(i_vec,AAAB) &
                                       & + dFxdgdg(i_vec,ABAB) / FOUR &
                                       ) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )
                               case (AB)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * (dFxdndg(i_vec,ABB)+dFxdndg(i_vec,AAB) / TWO)

                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * &
                                       (   dFxdgdg(i_vec,AABB) &
                                       & + dFxdgdg(i_vec,AAAB) &
                                       & + dFxdgdg(i_vec,ABAB) / FOUR &
                                       ) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )
                               end select
                            case (2) !! open shell
                               select case ((i_spin-1)*dim_factor+i_dim)
                               case (AA)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * TWO * dFxdndg(i_vec,AAA)

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * dFxdndg(i_vec,AAB)

                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * FOUR * dFxdgdg(i_vec,AAAA) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * TWO * dFxdgdg(i_vec,AAAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * TWO * dFxdgdg(i_vec,AAAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * dFxdgdg(i_vec,ABAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                               case (AB)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * TWO * dFxdndg(i_vec,ABB)

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * dFxdndg(i_vec,AAB)

                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * FOUR * dFxdgdg(i_vec,AABB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * TWO * dFxdgdg(i_vec,AAAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * TWO * dFxdgdg(i_vec,BBAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * dFxdgdg(i_vec,ABAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                               case (BB)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * TWO * dFxdndg(i_vec,BBB)

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * dFxdndg(i_vec,BAB)

                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * FOUR * dFxdgdg(i_vec,BBBB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * TWO * dFxdgdg(i_vec,BBAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * TWO * dFxdgdg(i_vec,BBAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * dFxdgdg(i_vec,ABAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )
                               case (BA)
                                  !! dFxdndg
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * TWO * dFxdndg(i_vec,ABB)

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit2) &
                                       + &
                                       (   &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       ) * orbs_ch_p(i_vec,i_fit) &
                                       ) * dFxdndg(i_vec,AAB)
                                  !! dfx_dgammadgamma
                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * FOUR * dFxdgdg(i_vec,AABB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * TWO * dFxdgdg(i_vec,AAAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,A) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,B) &
                                       ) * TWO * dFxdgdg(i_vec,BBAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )

                                  int = int + &
                                       grdwts(i_vec) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit) * grarho(i_vec,1,A) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit) * grarho(i_vec,2,A) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit) * grarho(i_vec,3,A) &
                                       ) * dFxdgdg(i_vec,ABAB) * &
                                       ( &
                                       grad_orbs_ch_p(i_vec,1,i_fit2) * grarho(i_vec,1,B) + &
                                       grad_orbs_ch_p(i_vec,2,i_fit2) * grarho(i_vec,2,B) + &
                                       grad_orbs_ch_p(i_vec,3,i_fit2) * grarho(i_vec,3,B) &
                                       )
                               end select
                            end select

                         enddo i_vec_

                         r2_index_matrix_v2(i_ir)%o(counter, (i_spin-1)*dim_factor+i_dim) = &
                              r2_index_matrix_v2(i_ir)%o(counter, (i_spin-1)*dim_factor+i_dim) &
                              + int/n_pa
                      enddo fit2
                   enddo fit
                enddo i_dim_
             enddo i_spin_
          enddo i_pa_ ! loop over partners
       enddo i_ir_ ! loop over irreps

       call stop_timer(timer_resp_2index_integration)

    enddo grid_points_ ! loop over gridpoints

    ! collect all results on the master
    if(comm_i_am_master()) then

       ! receiving results from the slaves
       if(output_response_detailed) call write_to_output_units(&
            & 'response_calc_2index_int_v2: receiving results from the slaves')
       call response_2index_receive_v2 (n_ir, n_spin)

       ! write final result to tape
       if(output_response_detailed) call write_to_output_units(&
            & 'response_calc_2index_int_v2: write results to tape')
       call response_write_2index_totape_v2(n_ir,n_spin,2)
    else
       if(output_response_detailed) call write_to_output_units(&
            & 'response_calc_2index_int_v2: sending the result to the master')
       ! sending the result to the master
       call response_2index_tomaster_v2 (n_ir)

    end if

    ! deallocate remaining things & clean up
    call fit_fct_free_response(fcts = fcts_orbs_ch, grads = grads_ch)

    do i_ir =1, n_ir
       deallocate(r2_index_matrix_v2(i_ir)%o, stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            & "response_calc_2index_int_v2: deallocation of r2_index_matrix_v2 failed")
    enddo

  end subroutine response_calc_2index_int_v3
  !*************************************************************

  subroutine response_calc_2index_int_v2()
    ! Purpose:
    ! calculate 2 index integrals <f_k| h_v |f_k>
    ! + parallelization over grid (like in post_scf_module)
    ! + integral prescreening using <f_k|f_l>
    !------------ Modules used ------------------- ---------------
    use vwnc
    use pw_ldac_module
    use gga_response_module
    use ch_response_module
    use linalg_module
    use exchange
    use constants
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    logical :: lda_case
    integer(kind=i4_kind) :: vla, &
         & i_atom,i_fit,i_fit2,&
         & counter, i_dim, i_spin !!, i_vec
    integer(i4_kind) :: alloc_stat, i_ir, n_pa, i_pa, n_ir, n_dim

!!! GGA variables
    real(kind=r8_kind), dimension(vec_length,2)  :: rho
    real(kind=r8_kind), dimension(vec_length,n_spin) :: nn

    real(kind=r8_kind), dimension(vec_length)   :: &
         Fx, Fc, eps
    real(kind=r8_kind), dimension(vec_length,2) :: &
         dFxdn, dFcdn
    real(kind=r8_kind), dimension(vec_length,3) :: &
         dFxdg, dFcdg, dFxdndn, dFcdndn
    real(kind=r8_kind), dimension(vec_length,4) :: &
         dFdndn, dFdg
    real(kind=r8_kind), dimension(vec_length,6) :: &
         dFxdndg, dFcdndg
    real(kind=r8_kind), dimension(vec_length,6)   :: &
         dFxdgdg, dFcdgdg

    real(kind=r8_kind) :: &
         grarho(vec_length,3,2),gamma(vec_length,3)

    type(orbital_type),pointer             :: fcts_orbs_ch(:)
    type(orbital_gradient_type),pointer    :: grads_ch(:)
    !!    type(orbital_sec_der_type),pointer     :: sec_ders_ch(:)

    real (r8_kind), pointer :: orbs_ch_p(:,:) ! pointer on component of orbs_ch
    real (r8_kind), pointer :: grad_orbs_ch_p(:,:,:)

    real(r8_kind), parameter :: FRTH = 0.25_r8_kind

    integer(i4_kind), parameter ::&
         & A = 1, &
         & B = 2

    integer(i4_kind), parameter ::&
         & XX = 1, &
         & YY = 2, &
         & ZZ = 3

    integer(i4_kind) :: ndc, idx, i_idx, i

    real(r8_kind)    :: CG2(2),CG3(4)
    integer(i4_kind) :: DG2(2),GG2(2),DG3(4),GG3(4),GG4(4)

    real(r8_kind),allocatable :: res_mat(:,:),gcr(:,:,:)
    real(r8_kind),allocatable :: LDA0(:,:),GGA1(:,:,:),GGA2(:,:,:),GGA3(:,:,:)

    integer(i4_kind), parameter ::&
         & AA = 1, &
         & AB = 2, &
         & BB = 3, &
         & BA = 4

    integer(i4_kind), parameter ::&
         & AAA = 1, &
         & BBB = 2, &
         & AAB = 3, &
         & ABB = 4, &
         & BAB = 5, &
         & BAA = 6

    integer(i4_kind), parameter ::&
         & AAAA = 1, &
         & BBBB = 2, &
         & ABAB = 3, &
         & AABB = 4, &
         & AAAB = 5, &
         & BBAB = 6

    !! For debug
    real(kind=r8_kind) :: &
         Flda(vec_length), dFldadn(vec_length,2), dFldadndn(vec_length,3),gw(vec_length)

    integer(i4_kind) :: CC

    integer(i4_kind) :: dim

    !------------ Executable code ------------------------------------

    if ( ( .not. xalpha_resp ) .and. &
         ( .not. vwn_resp    ) .and. &
         ( .not. beckex_resp  ) .and. &
         ( .not. perdewc_resp  ) .and. &
         ( .not. pw91x_resp  ) .and. &
         ( .not. pw91c_resp  ) .and. &
         ( .not. pbex_resp   ) .and. &
         ( .not. pbec_resp   ) .and. &
         ( .not. revpbex_resp   ) .and. &
         ( .not. pbenx_resp   ) .and. &
         ( .not. pw_ldac_resp)) return

    lda_case = .false.
    if (  xalpha_resp ) lda_case = .true.
    if (     vwn_resp ) lda_case = .true.
    if ( pw_ldac_resp ) lda_case = .true.

    FPP_TIMER_START(all)

    n_ir       = ssym%n_irrep
    n_spin     = ssym%n_spin
    dim_factor = 2

    ! allocate memory for 2-index-integrals
    allocate(r2_index_matrix_v2(n_ir), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    do i_ir = 1,n_ir
       n_dim = dimension_of_fit_ch(i_ir)
       dim = (n_dim*(n_dim+1))/2
!!       print *,"dim = ",dim
       allocate(r2_index_matrix_v2(i_ir)%o(dim,n_spin*dim_factor),&
            & STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       r2_index_matrix_v2(i_ir)%o = zero
    enddo

    call fit_fct_allocate_response(fcts_ch = fcts_orbs_ch, grads = grads_ch)

    vla = machineparameters_veclen

    call orbital_setup_response(vla,.true.,.true.)

    !
    ! Loop over grid points:
    !
    call grid_loop_setup_atom()

    grid_points_: do while( more_grid_atom(vec_length, i_atom, grdpts, grdwts) )
       ASSERT(i_atom/=0)

       vla = size(grdpts,1)


       ! calculate gridweights
       gw(1:vla) = grdwts(1:vla)&
            &    * atomicweight(i_atom, grdpts(1:vla, :)) &
            &    * unique_atoms(i_atom)%n_equal_atoms

       ! calculate orbitals and charge fitfunctions on the grid
       ! (orbitals are needed for the density below)
       FPP_TIMER_START(rho_calc)

       call orbital_calculate(grdpts(1:vla,1:3),vla,orbs_ob,orbs_grads)

       nn     = zero
       rho    = zero
       gamma  = zero
       grarho = zero

       call density_calc_nl(vla,nn,gamma,grarho,orbs_ob,orbs_grads)

       if (n_spin==1) then
          rho(:,1) = nn(:,1)/TWO
          rho(:,2) = nn(:,1)/TWO

          gamma(:,2) = gamma(:,1)/FOUR
          gamma(:,3) = gamma(:,1)/FOUR
          gamma(:,1) = gamma(:,1)/FOUR
       else
          rho(:,1) = nn(:,1)
          rho(:,2) = nn(:,2)
       end if

       !! CUT OFF
       where ( abs(rho) < rho_cutoff )
          rho   = rho_cutoff
       endwhere

       FPP_TIMER_STOP(rho_calc)

       ! calculate contributions from the functionals
       FPP_TIMER_START(XC_calc)

       !! EXCHANGE PART
       Flda      = zero
       dFldadn   = zero
       dFldadndn = zero

       Fx        = zero
       dFxdn     = zero
       dFxdg     = zero
       dFxdndn   = zero
       dFxdndg   = zero
       dFxdgdg   = zero
       if ( xalpha_resp )  &
            call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
       if ( pbex_resp   )  then
          call exchange_lda(X_XALPHA,vla,2,rho,Flda,dFldadn,dFldadndn)
          call exchange_gga(X_PBE   ,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
          Fx      = Fx      + Flda
          dFxdn   = dFxdn   + dFldadn
          dFxdndn = dFxdndn + dFldadndn
       end if

       !! CORRELATION PART
       Flda      = zero
       dFldadn   = zero
       dFldadndn = zero

       Fc        = zero
       dFcdn     = zero
       dFcdg     = zero
       dFcdndn   = zero
       dFcdndg   = zero
       dFcdgdg   = zero

       if ( vwn_resp    )  then
          eps(1:vla) = 1.0E-32_r8_kind
          call vwn_calcMDA(rho,dFcdn,2,Fc,vla,eps,dFcdndn)
          dFldadndn(:,3) = dFcdndn(:,3)
          dFcdndn(:,3)   = dFcdndn(:,2)
          dFcdndn(:,2)   = dFldadndn(:,3)
          dFldadndn      = 0.0_r8_kind
       end if

       if ( pw_ldac_resp )  &
            call pw_ldac   (vla,4,rho,Fc,dFcdn,dFcdndn) !! 4 = unrestr case
       if ( pbec_resp   )  then
          call pw_ldac   (vla,4,rho,Flda,dFldadn,dFldadndn)
          call gga_correlation(C_PBE,2,2,rho,gamma,vla,&
               Fc,Flda,&
               dFcdn,dFcdg,dFldadn,&
               dFcdndn,dFcdndg,dFcdgdg,dFldadndn)
          Fc      = Fc      + Flda
          dFcdn   = dFcdn   + dFldadn
          dFcdndn = dFcdndn + dFldadndn
       end if

       if (n_spin == 1) then
          CC = 2
       else
          CC = 1
       end if

       dFxdg   = (dFxdg  +dFcdg  ) / CC
       dFxdndn = (dFxdndn+dFcdndn) / CC
       dFxdndg = (dFxdndg+dFcdndg) / CC
       dFxdgdg = (dFxdgdg+dFcdgdg) / CC

       dFdg   = zero
       dFdndn = zero

       select case (n_spin)
       case (1)
          dFdg  (1:vla,AA) = gw(1:vla) * dFxdg(1:vla,1)
          dFdg  (1:vla,AB) = gw(1:vla) * dFxdg(1:vla,3) / TWO !! will be x2 later

          dFdndn(1:vla,AA) = gw(1:vla) * dFxdndn(1:vla,AA)
          dFdndn(1:vla,AB) = gw(1:vla) * dFxdndn(1:vla,AB)
       case (2)
          dFdg  (1:vla,AA) = gw(1:vla) * dFxdg(1:vla,1)
          dFdg  (1:vla,AB) = gw(1:vla) * dFxdg(1:vla,3) / TWO !! will be x2 later
          dFdg  (1:vla,BB) = gw(1:vla) * dFxdg(1:vla,2)
          dFdg  (1:vla,BA) = gw(1:vla) * dFxdg(1:vla,3) / TWO !! will be x2 later

          dFdndn(1:vla,AA) = gw(1:vla) * dFxdndn(1:vla,AA)
          dFdndn(1:vla,AB) = gw(1:vla) * dFxdndn(1:vla,AB)
          dFdndn(1:vla,BB) = gw(1:vla) * dFxdndn(1:vla,BB)
          dFdndn(1:vla,BA) = gw(1:vla) * dFxdndn(1:vla,AB)
       end select

       do i = 1,size(dFxdndg,2)
          dFxdndg(1:vla,i)  = gw(1:vla) * dFxdndg(1:vla,i)
          dFxdgdg(1:vla,i)  = gw(1:vla) * dFxdgdg(1:vla,i)
       end do

       FPP_TIMER_STOP(XC_calc)

       i_ir1_: do i_ir = 1,n_ir
          fcts_orbs_ch(i_ir)%o = zero
          grads_ch(i_ir)%o     = zero
       enddo i_ir1_

       FPP_TIMER_START(khi_calc)
       call fit_fct_calculate_response(grdpts(1:vla,1:3),vla, &
            fcts_orbs_ch, grads_ch )
       FPP_TIMER_STOP(khi_calc)

       call start_timer(timer_resp_2index_integration)

       FPP_TIMER_START(loop_calc)
       i_ir_: do i_ir = 1,n_ir
          n_pa = symmetry_data_n_partners(i_ir)
          ndc  = dimension_of_fit_ch(i_ir)
          if (ndc==0) exit
          allocate(res_mat(ndc,ndc), LDA0(vla,ndc), stat = alloc_stat)
          ASSERT(alloc_stat == 0)

          LDA0    = zero

          if (.not. lda_case) then
             allocate(GGA1(vla,ndc,3), &
                  &   GGA2(vla,ndc,2), &
                  &   GGA3(vla,ndc,4), &
                  &   gcr  (vla,ndc,n_spin), stat = alloc_stat)
             ASSERT(alloc_stat == 0)
          end if

          i_pa_: do i_pa = 1, n_pa

             orbs_ch_p      => fcts_orbs_ch(i_ir)%o(:,:,i_pa)
             grad_orbs_ch_p => grads_ch(i_ir)%o(:,:,:,i_pa)

             i_spin_: DO i_spin = 1,n_spin

                i_dim_:  do i_dim=1,dim_factor

                   idx = (i_spin-1)*dim_factor+i_dim

                   select case (n_spin)
                   case (1) !! closed shell
                      select case (idx)
                      case (AA)
                         CG2 = (/ONE, HALF/)
                         DG2 = (/AAA,  AAB/)
                         GG2 = (/  A,    A/)

                         CG3 = (/ONE,  HALF, HALF, FRTH/)
                         DG3 = (/AAAA, AAAB, AAAB, ABAB/)
                         GG3 = (/   A,    A,    A,    A/)
                         GG4 = (/   A,    A,    A,    A/)
                      case (AB)
                         CG2 = (/ONE, HALF/)
                         DG2 = (/ABB,  AAB/)
                         GG2 = (/  A,    A/)

                         CG3 = (/ONE,  HALF, HALF, FRTH/)
                         DG3 = (/AABB, AAAB, AAAB, ABAB/)
                         GG3 = (/   A,    A,    A,    A/)
                         GG4 = (/   A,    A,    A,    A/)
                      end select
                   case (2) !! open shell
                      select case (idx)
                      case (AA)
                         CG2 = (/TWO, ONE/)
                         DG2 = (/AAA, AAB/)
                         GG2 = (/  A,   B/)

                         CG3 = (/FOUR,  TWO,  TWO,  ONE/)
                         DG3 = (/AAAA, AAAB, AAAB, ABAB/)
                         GG3 = (/   A,    A,    B,    B/)
                         GG4 = (/   A,    B,    A,    B/)
                      case (AB)
                         CG2 = (/TWO, ONE/)
                         DG2 = (/ABB, AAB/)
                         GG2 = (/  B,   A/)

                         CG3 = (/FOUR,  TWO,  TWO,  ONE/)
                         DG3 = (/AABB, AAAB, BBAB, ABAB/)
                         GG3 = (/   A,    A,    B,    A/)
                         GG4 = (/   B,    A,    B,    B/)
                      case (BB)
                         CG2 = (/TWO, ONE/)
                         DG2 = (/BBB, BAB/)
                         GG2 = (/  B,   A/)

                         CG3 = (/FOUR,  TWO,  TWO,  ONE/)
                         DG3 = (/BBBB, BBAB, BBAB, ABAB/)
                         GG3 = (/   B,    B,    A,    A/)
                         GG4 = (/   B,    A,    B,    A/)
                      case (BA)
                         CG2 = (/TWO, ONE/)
!!                         DG2 = (/ABB, AAB/)
!!                         GG2 = (/  B,   A/)
                         DG2 = (/BAA, BAB/)
                         GG2 = (/  A,   B/)

                         CG3 = (/FOUR,  TWO,  TWO,  ONE/)
                         DG3 = (/AABB, AAAB, BBAB, ABAB/)
                         GG3 = (/   A,    A,    B,    A/)
                         GG4 = (/   B,    A,    B,    B/)
                      end select
                   end select

                   if (.not. lda_case) gcr  = zero

                   i_idx_: do i_idx = 1, vla
                      LDA0(i_idx,1:ndc) = orbs_ch_p(i_idx,1:ndc) * dFdndn(i_idx,idx)
                      if (.not. lda_case) then
                         i_: do i = XX,ZZ
                            GGA1(i_idx,1:ndc,i) = &
                                 &                TWO &
                                 &              * dFdg(i_idx,idx) &
                                 &              * grad_orbs_ch_p(i_idx,i,1:ndc)
                            gcr(i_idx,1:ndc,A) = gcr(i_idx,1:ndc,A) + grad_orbs_ch_p(i_idx,i,1:ndc) * grarho(i_idx,i,A)
                            if (n_spin==2) &
                                 gcr(i_idx,1:ndc,B) = gcr(i_idx,1:ndc,B) + grad_orbs_ch_p(i_idx,i,1:ndc) * grarho(i_idx,i,B)
                         end do i_
                         ! dFdndg
                         do i = 1,2
                            GGA2(i_idx,1:ndc,i) = &
                                 &                CG2(i) * dFxdndg(i_idx,DG2(i)) &
                                 &              * gcr(i_idx,1:ndc,GG2(i))
                         enddo
                         ! dFdgdg
                         do i = 1,4
                            GGA3(i_idx,1:ndc,i) = &
                                 &                CG3(i) * dFxdgdg(i_idx,DG3(i)) &
                                 &              * gcr(i_idx,1:ndc,GG3(i))
                         end do
                      end if
                   end do i_idx_

                   res_mat = zero

                   ! dFdndn
                   call matmatmul(orbs_ch_p(1:vla,1:ndc),LDA0,res_mat,'T','N')

                   if (.not. lda_case) then
                      ! dFdg
                      do i = XX,ZZ
                         call matmatmul(grad_orbs_ch_p(1:vla,i,1:ndc),GGA1(:,:,i),res_mat,'T','N',&
                              ONE, ONE)
                      end do
                      ! dFdndg
                      do i = 1,2
                         call matmatmul(orbs_ch_p(1:vla,1:ndc),       GGA2(:,:,i),res_mat,'T','N',&
                              TWO, ONE)
                      end do
                      ! dFdgdg
                      do i = 1,4
                         call matmatmul(gcr(1:vla,1:ndc,GG4(i)),      GGA3(:,:,i),res_mat,'T','N',&
                              ONE, ONE)
                      end do
                   end if

                   counter=0   ! combined index for packed storage of charge_prematrix
                   fit: do i_fit=1, ndc ! dimension of ch fit fct
                      fit2: do i_fit2=1,i_fit
                         counter = counter + 1
                         r2_index_matrix_v2(i_ir)%o(counter,idx) = r2_index_matrix_v2(i_ir)%o(counter,idx) + &
                              (res_mat(i_fit2,i_fit)+res_mat(i_fit,i_fit2))/n_pa/TWO
                      end do fit2
                   end do fit
                enddo i_dim_
             enddo i_spin_
          enddo i_pa_

          deallocate(LDA0, res_mat, stat = alloc_stat)
          ASSERT(alloc_stat==0)
          if (.not. lda_case) then
             deallocate(GGA1,GGA2,GGA3,gcr,stat = alloc_stat)
             ASSERT(alloc_stat==0)
          end if

       enddo i_ir_
       FPP_TIMER_STOP(loop_calc)

       call stop_timer(timer_resp_2index_integration)

    enddo grid_points_ ! loop over gridpoints

    ! collect all results on the master
    if(comm_i_am_master()) then

       ! receiving results from the slaves
       call write_to_output_units('response_calc_2index_int_v2: receiving results from the slaves')
       call response_2index_receive_v2 (n_ir, n_spin)

       ! write final result to tape
       call write_to_output_units('response_calc_2index_int_v2: write results to tape')
       call response_write_2index_totape_v2(n_ir,n_spin,2)
    else
       call write_to_output_units('response_calc_2index_int_v2: sending the result to the master')
       ! sending the result to the master
       call response_2index_tomaster_v2 (n_ir)
    end if

    ! deallocate remaining things & clean up
    call fit_fct_free_response(fcts = fcts_orbs_ch, grads = grads_ch)

    do i_ir =1, n_ir
       deallocate(r2_index_matrix_v2(i_ir)%o, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    enddo

    FPP_TIMER_STOP(all)
#ifdef FPP_TIMERS
    print *,'2c_resp TIMER: rho_calc  = ',FPP_TIMER_VALUE(rho_calc)
    print *,'2c_resp TIMER: XC_calc   = ',FPP_TIMER_VALUE(XC_calc)
    print *,'2c_resp TIMER: khi_calc  = ',FPP_TIMER_VALUE(khi_calc)
    print *,'2c_resp TIMER: loop_calc = ',FPP_TIMER_VALUE(loop_calc)
    print *,'2c_resp TIMER: all       = ',FPP_TIMER_VALUE(all)
#endif
  end subroutine response_calc_2index_int_v2


  !--------------- End of module -------------------------------------
end module response_module

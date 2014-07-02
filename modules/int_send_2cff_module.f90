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
module  int_send_2cff_module
!---------------------------------------------------------------
!
!  Purpose: Sending, receiving and storing of 
!           2 center fitfunction integrals.
!           The results of one quadrupel stored in
!           int_data_2cff_module are distributed.
!           When packing the integrals, the mapping to the scf-
!           metaindex is also done.
!           The integrals are send to the master.
!           
!
!
!  Module called by: integral_main_2cff, 
!    integral_calc_quad_2cff
!
!  References: Publisher Document: Concepts of Integral Part
! 
!
!  Author: TB
!  Date: 7/95
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!  Modification (Please copy before editing)
!  Author:      Uwe Birkenheuer
!  Date:        27.5.98
!  Description: algorithm for glob_glob_int summation corrected
! 
!  Modification (Please copy before editing)
!  Author:      Uwe Birkenheuer
!  Date:        8/98
!  Description: treatment of mixed exchange overlap integrals introduced
!               to this end, a new and more flexible algorithm for 
!               dealing with the meta storage index implemented
!
! Modification Master/Slave concept to DLB
! Author: AN
! Date:   4/11
! Description: for scheduling the 2cff integrals DLB is used
!              it replaces the master/slave concept
!              This removes all messages related to do or done of jobs
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
use integralpar_module
use output_module
use iounitadmin_module
use filename_module
use comm_module
use msgtag_module
use symmetry_data_module, only: get_totalsymmetric_irrep
use unique_atom_module
use int_data_2cff_module
use quadrupel_module
use timer_module
use time_module

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
public int_send_2cff_setup, int_send_2cff_shutdown, &
       int_send_2cff_send, int_send_2cff_receive


  !===================================================================
  ! End of public interface of module
  !===================================================================

!................................................................
! triangle storage mode
! ~~~~~~~~~~~~~~~~~~~~~
! 
!      |   H   |   O   |   For each combination of contraction types
!      | s p d | s p d |   loc_loc, loc_glob, glob_loc, and glob_glob,
!  ----+-------+-------+   only the quadruples contributions shown in
!    s | X     |       |   Fig. 1 are actually evaluated by the various
!  H p | X X   |       |   calls of the subroutine FITCONTRACT_2C.
!    d | X X X |       |   
!  ----+-------+-------+   These individual contributions have to be
!    s | X X X | X     |   gathered properly in the final 2-center
!  O p | X X X | X X   |   integral matrix F_ij.
!    d | X X X | X X X |   
!  ----+-------+-------+ FIG. 1
!
!
!         |  local fcts   | glob cont |    The order of the matrix 
!         |   H   |   O   |           |    elements in the F-matrix 
!         | s p d | s p d |  H  |  O  |    is shown in Fig. 2. Only 
!  -------+-------+-------+-----+-----+    the entries in the lower
!       s | X     |       |     |     |    triangle are stored. To this
!  l  H p | X X   |       |     |     |    end linear packed storage mode 
!  o    d | X X X |       |     |     |    according to A_ij -> A_k with
!  c  ----+-------+-------+-----+-----+    k = i*(i-1)/2 + j is applied.
!  a    s | X X X | X     |     |     |    
!  l  O p | X X X | X X   |     |     |    In addition the summation of 
!       d | X X X | X X X |     |     |    the individual l-contributions
!  -------+-------+-------+-----+-----+    to the global contractions have
!  g      |       |       |     |     |    to be carried out properly.
!  l    H | X X X | X X X |  X  |     |    
!  o    --+-------+-------+-----+-----+    Details on how this is achieved
!  b      |       |       |     |     |    are given below.
!  .    O | X X X | X X X |  X  |  X  |
!  -------+-------+-------+-----+-----+ FIG. 2
!
!
!  a) The lower triangle contributions [ l_i | l_j ] to the local-local
!     sub block of the F-matrix (Fig. 2) can directly be taken from the 
!     intermediates loc_loc(i,j) (see Fig. 1)
! 
!  b) The set of contributions [ g_n | l_j ] = Sum(l) [ g_n,l | l_j ] to
!     the entirely stored global-local sub block of the F-matrix (Fig. 2)
!     are split into two data sets (see Fig. 3)
!
!           |   H   |   O   |           |   H   |   O   |
!           | s p d | s p d |           | s p d | s p d |
!       ----+-------+-------+       ----+-------+-------+
!         s | X     |       |         s |   X X | X X X |
!       H p | X X   |       |       H p |     X | X X X |
!         d | X X X |       |         d |       | X X X |
!       ----+-------+-------+       ----+-------+-------+
!         s | X X X | X     |         s |       |   X X |
!       O p | X X X | X X   |       O p |       |     X |
!         d | X X X | X X X |         d |       |       |
!       ----+-------+-------+       ----+-------+-------+ FIG. 3
!
!    The first data set (left-hand side) can directly be taken from the
!    lower triangle intermediates glob_loc(n,j). The second data set
!    (right-hand side) is not availible directly. However, using the
!    relation [ g_n,l | l_j ] = [ l_j | g_n,l ] is can be mapped to the
!    lower triangle quadruple intermediates loc_glob(j,n) which are 
!    provided by FITCONTRACT_2C (see Fig. 1). The diagonal quadruples
!    are not involved in this mapping.
!
! c) The quadruple contributions required to load the lower triangle
!    entries [ g_n | g_m ] = Sum(l,l') [ g_n,l | g_m,l' ] of the 
!    F-matrix (Fig. 2) are displayed in Fig. 4.
!
!           |   H   |   O   |                |   H   |   O   |
!           | s p d | s p d |                | s p d | s p d |
!       ----+-------+-------+            ----+-------+-------+
!         s | X X X |       |              s |   X X |       |
!       H p | X X X |       |            H p |     X |       |
!         d | X X X |       |              d |       |       |
!       ----+-------+-------+            ----+-------+-------+
!         s | X X X | X X X |              s |       |   X X |
!       O p | X X X | X X X |            O p |       |     X |
!         d | X X X | X X X |              d |       |       |
!       ----+-------+-------+ FIG. 4     ----+-------+-------+ FIG. 5
!
!    The lower triangle contributions are directly availible via
!    glob_glob_(n,m). The missing quadruple contributions are shown
!    in Fig. 5. As in the case of the mixed matrix elements the 
!    symmetry relation [ g_n,l | g_m,l' ] = [ g_m,l' | g_n,l ] can
!    be used to map these missing terms onto the lower triangle
!    intermediates glob_glob(m,n). This has only to be done for quadrupels
!    with equal unique atoms but different l-quantum numbers l` > l.
!..............................................................................

!------------ Declaration of private types ----

! fitfct_index_type substitutes the old type index_ua1_type
type, private :: loc_index_type
   integer(kind=i4_kind)          :: N_lc
   integer(kind=i4_kind), pointer :: index(:,:) ! (i_exp,i_ind)
end type loc_index_type

type, private :: fitfct_index_type
   type(loc_index_type), pointer  :: loc_index(:)  ! (i_bas)
   integer(kind=i4_kind)          :: N_gc
   integer(kind=i4_kind), pointer :: glob_index(:) ! (i_gc)
end type fitfct_index_type
!------------ Declaration of constants and variables ----

! only on master:
integer(kind=i4_kind) :: n_missing_quadrupels, &
     n_tot_records_ch, n_tot_records_xc, n_tot_records_mx

real(kind=r8_kind), allocatable :: integral_2cch_no(:), integral_2cxc_no(:), &
                                   integral_2c_mixed(:)
! variables xx_fitfct_index replaces the old variables metaindex_of_..._xx
type(fitfct_index_type), allocatable, target :: ch_fitfct_index(:), & ! (i_ua)
                                                xc_fitfct_index(:)    ! (i_ua)
integer(kind=i4_kind), private               :: n_ch, n_xc
real(kind=r8_kind), allocatable :: &
     integral_2cch_pre(:)

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine fitfct_index_setup
    !  Purpose: 
    !     calculates the help arrays for mapping integrals to the
    !     1-dimensional data storage of results
    ! called by: int_send_2cff_setup (formerly done by metaindex_setup)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    implicit none
    integer(kind=i4_kind)                   :: status, pos, i_ir, i_ua, &
                                               i_bas, i_l, i_exp, i_ind, i_gc
    type(unique_atom_type)        , pointer :: ua
    type(unique_atom_basis_type)  , pointer :: uab
    type(unique_atom_partner_type), pointer :: uap
    type(loc_index_type)          , pointer :: loc_index(:)
    integer(kind=i4_kind)         , pointer :: glob_index(:)
    integer(kind=i4_kind)         , pointer :: index(:,:)
    integer(kind=i4_kind)         , pointer :: n_rec
    !------------ Executable code ------------------------------------
 
    i_ir = get_totalsymmetric_irrep()
 
    if (integralpar_2cch_no .or. integralpar_2c_mixed) then
       ! << load charge fit function index >>
       allocate( ch_fitfct_index(N_unique_atoms), stat=status )
       if (status /= 0) call error_handler( "fitfct_index_setup: &
            &allocate of ch_fitfct_index failed" )
       pos = 0
       ! << first all local contractions >>
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          allocate( loc_index(-1:ua%lmax_ch), stat=status )
          if (status /= 0) call error_handler( "fitfct_index_setup: &
               &allocate of ch_fitfct_index(i_ua)%loc_index failed" )
          ch_fitfct_index(i_ua)%loc_index => loc_index
          do i_bas=-1,ua%lmax_ch
             select case (i_bas)
             case (-1) ! s-type fit functions
                i_l =  0
                uab => ua%l_ch(0)
             case (0) ! r^2-type fit functions
                i_l =  0
                uab => ua%r2_ch
             case default ! l-type fit functions
                i_l =  i_bas
                uab => ua%l_ch(i_l)
             end select
             uap => ua%symadapt_partner(i_ir,i_l)
             allocate( index(uab%N_uncontracted_fcts + uab%N_contracted_fcts, &
                             uap%N_independent_fcts), stat=status )
             if (status /= 0) call error_handler( "fitfct_index_setup: &
                  &allocate of ch_fitfct_index(i_ua)%loc_index(i_bas)%index & 
                  &failed" )
             loc_index(i_bas)%index => index
             n_rec => loc_index(i_bas)%N_lc
             n_rec = pos
             do i_exp=1,uab%N_uncontracted_fcts + uab%N_contracted_fcts
                do i_ind=1,uap%N_independent_fcts
                   pos = pos + 1
                   index(i_exp,i_ind) = pos
                end do
             end do
             n_rec = pos - n_rec
          end do ! i_bas
       end do ! i_ua
       ! << now all global contractions >>
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          allocate( glob_index(ua%N_glob_cons_ch), stat=status )
          if (status /= 0) call error_handler( "fitfct_index_setup: &
               &allocate of ch_fitfct_index(i_ua)%glob_index failed" )
          ch_fitfct_index(i_ua)%glob_index => glob_index
          n_rec => ch_fitfct_index(i_ua)%N_gc
          n_rec = pos
          do i_gc=1,ua%N_glob_cons_ch
             pos = pos + 1
             glob_index(i_gc) = pos
          end do
          n_rec = pos - n_rec
       end do ! i_ua
       n_ch = pos
    end if
 
    if (integralpar_2cxc_no .or. integralpar_2c_mixed) then
       ! << load exchange fit function index >>
       allocate( xc_fitfct_index(N_unique_atoms), stat=status )
       if (status /= 0) call error_handler( "fitfct_index_setup: &
            &allocate of xc_fitfct_index failed" )
       pos = 0
       ! << first all local contractions >>
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          allocate( loc_index(-1:ua%lmax_xc), stat=status )
          if (status /= 0) call error_handler( "fitfct_index_setup: &
               &allocate of xc_fitfct_index(i_ua)%loc_index failed" )
          xc_fitfct_index(i_ua)%loc_index => loc_index
          do i_bas=-1,ua%lmax_xc
             select case (i_bas)
             case (-1) ! s-type fit functions
                i_l =  0
                uab => ua%l_xc(0)
             case (0) ! r^2-type fit functions
                i_l =  0
                uab => ua%r2_xc
             case default ! l-type fit functions
                i_l =  i_bas
                uab => ua%l_xc(i_l)
             end select
             uap => ua%symadapt_partner(i_ir,i_l)
             allocate( index(uab%N_uncontracted_fcts + uab%N_contracted_fcts, &
                             uap%N_independent_fcts), stat=status )
             if (status /= 0) call error_handler( "fitfct_index_setup: &
                  &allocate of xc_fitfct_index(i_ua)%loc_index(i_bas)%index & 
                  &failed" )
             loc_index(i_bas)%index => index
             n_rec => loc_index(i_bas)%N_lc
             n_rec = pos
             do i_exp=1,uab%N_uncontracted_fcts + uab%N_contracted_fcts
                do i_ind=1,uap%N_independent_fcts
                   pos = pos + 1
                   index(i_exp,i_ind) = pos
                end do
             end do
             n_rec = pos - n_rec
          end do ! i_bas
       end do ! i_ua
       ! << now all global contractions >>
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          allocate( glob_index(ua%N_glob_cons_xc), stat=status )
          if (status /= 0) call error_handler( "fitfct_index_setup: &
               &allocate of xc_fitfct_index(i_ua)%glob_index failed" )
          xc_fitfct_index(i_ua)%glob_index => glob_index
          n_rec => xc_fitfct_index(i_ua)%N_gc
          n_rec = pos
          do i_gc=1,ua%N_glob_cons_xc
             pos = pos + 1
             glob_index(i_gc) = pos
          end do
          n_rec = pos - n_rec
       end do ! i_ua
       n_xc = pos
    end if
 
    if (integralpar_2cch_no) then
       n_tot_records_ch = (n_ch*(n_ch+1))/2 ! triangle storage mode
       if (output_int_loops) call write_to_output_units("fitfct_index_setup:& 
            & total number of ch records: ",n_tot_records_ch)
    end if
    if (integralpar_2cxc_no) then
       n_tot_records_xc = (n_xc*(n_xc+1))/2 ! triangle storage mode
       if (output_int_loops) call write_to_output_units("fitfct_index_setup:& 
            & total number of xc records: ",n_tot_records_xc)
    end if
    if (integralpar_2c_mixed) then
       n_tot_records_mx = n_ch * n_xc
       if (output_int_loops) call write_to_output_units("fitfct_index_setup:& 
            & total number of mixed records: ",n_tot_records_mx)
    end if
 
  end subroutine fitfct_index_setup
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cff_setup(n_quads)
    !  Purpose: 
    !   + calculates help variables number of quadrupels
    !     to be received and total number of record in file
    !   + calculates help array for mapping integrals to
    !     1 dimensional data storage of results
    !   + Allocating storage for results
    ! called by: integral_main_2cff (only by master)
    !** End of interface *****************************************
    implicit none
    integer(kind=i4_kind), intent(in) :: n_quads
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: status

    ! calculate help arrays for mapping the integrals
    call fitfct_index_setup
 
    ! number of quadrupels needed to be send to master is number of
    ! all quadrupels to store
    n_missing_quadrupels = n_quads

    if ( output_int_progress ) then
       call write_to_output_units("")
       call write_to_output_units( &
            "total number of 2 center fitfunction quadrupels: ",n_quads)
       call write_to_output_units("")
    endif


    ! allocating storage for results
    if (integralpar_2cch_no) then
       allocate( integral_2cch_no(n_tot_records_ch), stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_setup: allocate of integral_2cch_no failed" )
       integral_2cch_no = 0.0_r8_kind
    endif
    if (integralpar_2cch_pre) then
       allocate( integral_2cch_pre(n_tot_records_ch), stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_setup: allocate of integral_2cch_pre failed" )
       integral_2cch_pre = 0.0_r8_kind
    endif
    if (integralpar_2cxc_no) then
       allocate( integral_2cxc_no(n_tot_records_xc), stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_setup: allocate of integral_2cxc_no failed" )
       integral_2cxc_no = 0.0_r8_kind
    endif
    if (integralpar_2c_mixed) then
       allocate( integral_2c_mixed(n_tot_records_mx), stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_setup: allocate of integral_2c_mixed failed" )
       integral_2c_mixed = 0.0_r8_kind
    endif

  end subroutine int_send_2cff_setup
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cff_shutdown
    !  Purpose:
    !    + waiting for missing integrals
    !    + deallocation of help arrays
    !    + writing results as sequential files as used in the scf part
    !    + waiting for messages of slaves that they finished shutdown
    ! called by: integral_main_2cff (only on master)
    !** End of interface *****************************************
    use operations_module
    use fit_coeff_module, only: fit_coeff_store_norm
    use mat_charge_module, only: mat_charge_set
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: status
    logical :: restart_timer
    !------------ Executable code ------------------------------------

    ! waiting for missing integrals
    if ( operations_integral ) then
       do while (n_missing_quadrupels .gt. 0 )
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_shutdown: waiting for quadrupels: ",n_missing_quadrupels)
          call comm_save_recv(comm_all_other_hosts, msgtag_int_2cff_result)
          if ( output_slaveoperations ) &
               call write_to_output_units("int_send_2cff_shutdown: integral_2cob_result")
          call int_send_2cff_receive()
       enddo
    endif

    
    if ( timer_int_idle_2cff(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_idle_2cff(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_rewrite_2cff(integralpar_i_int_part))

    ! writing results as sequential files as used in the scf part
    if (operations_integral) then

       if (  integralpar_2cch_no) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_shutdown: setting Coulomb metric")

          !
          ! This picks only the diagonal elements:
          !
          call fit_coeff_store_norm(integral_2cch_no)

          !
          ! This data was previousely written to file
          ! "mat_charge.dat":
          !
          call mat_charge_set(integral_2cch_no)
       endif

       if (  integralpar_2cch_pre) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_shutdown: writing prescreen")
          call write_integralfile( "ch", &
               trim(tmpfile("pre_charge.dat")), integral_2cch_pre)
       endif
       if ( integralpar_2cxc_no) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_shutdown: writing anal. exchange")
          call write_integralfile( "xc", &
               trim(tmpfile("mat_exchange.dat")), integral_2cxc_no)
       endif
       if ( integralpar_2c_mixed) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_shutdown: writing mixed overlap")
          call write_integralfile( "mx", &
               trim(tmpfile("mat_mixed.dat")), integral_2c_mixed)
       endif
    endif

    ! deallocation of help arrays
    call fitfct_index_shutdown

    ! deallocating storage for results
    if (integralpar_2cch_no) then
       deallocate( integral_2cch_no, stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_shutdown: deallocate of integral_2cch_no failed" )
    endif
    if (integralpar_2cch_pre) then
       deallocate( integral_2cch_pre, stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_shutdown: deallocate of integral_2cch_pre failed" )
    endif
    if (integralpar_2cxc_no) then
       deallocate( integral_2cxc_no, stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_shutdown: deallocate of integral_2cxc_no failed" )
    endif
    if (integralpar_2c_mixed) then
       deallocate( integral_2c_mixed, stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cff_shutdown: deallocate of integral_2c_mixed failed" )
    endif


    call stop_timer(timer_int_rewrite_2cff(integralpar_i_int_part))
    if (restart_timer) call start_timer(timer_int_idle_2cff(integralpar_i_int_part))

  end subroutine int_send_2cff_shutdown
  !*************************************************************
 

  !*************************************************************
  subroutine fitfct_index_shutdown
    ! Purpuse:
    !   deallocation of the help array for the intgral mapping
    ! called by: int_send_2cff_shutdown (formerly done by metaindex_setup)
    !** End of interface *************************************
    !------------ Declaration of local variables -------------
    implicit none
    integer(kind=i4_kind)           :: status, i_ua, i_bas
    type(unique_atom_type), pointer :: ua
    type(loc_index_type)  , pointer :: loc_index(:)

    if (integralpar_2cch_no .or. integralpar_2c_mixed) then
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          loc_index => ch_fitfct_index(i_ua)%loc_index
          do i_bas=-1,ua%lmax_ch
             deallocate( loc_index(i_bas)%index, stat=status )
             if (status /= 0) call error_handler("fitfct_index_shutdown: &
                  &deallocate of ch_fitfct_index(i_ua)%loc_index(i_bas)%index &
                  &failed" )
          end do
          deallocate( ch_fitfct_index(i_ua)%loc_index, stat=status )
          if (status /= 0) call error_handler("fitfct_index_shutdown: &
               &deallocate of ch_fitfct_index(i_ua)%loc_index failed" )
          deallocate( ch_fitfct_index(i_ua)%glob_index, stat=status )
          if (status /= 0) call error_handler("fitfct_index_shutdown: &
               &deallocate of ch_fitfct_index(i_ua)%glob_index failed" )
       end do
       deallocate( ch_fitfct_index, stat=status )
       if (status /= 0) call error_handler("fitfct_index_shutdown: &
            &deallocate of ch_fitfct_index failed" )
    end if

    if (integralpar_2cxc_no .or. integralpar_2c_mixed) then
       do i_ua=1,N_unique_atoms
          ua => unique_atoms(i_ua)
          loc_index => xc_fitfct_index(i_ua)%loc_index
          do i_bas=-1,ua%lmax_xc
             deallocate( loc_index(i_bas)%index, stat=status )
             if (status /= 0) call error_handler("fitfct_index_shutdown: &
                  &deallocate of xc_fitfct_index(i_ua)%loc_index(i_bas)%index &
                  &failed" )
          end do
          deallocate( xc_fitfct_index(i_ua)%loc_index, stat=status )
          if (status /= 0) call error_handler("fitfct_index_shutdown: &
               &deallocate of xc_fitfct_index(i_ua)%loc_index failed" )
          deallocate( xc_fitfct_index(i_ua)%glob_index, stat=status )
          if (status /= 0) call error_handler("fitfct_index_shutdown: &
               &deallocate of xc_fitfct_index(i_ua)%glob_index failed" )
       end do
       deallocate( xc_fitfct_index, stat=status )
       if (status /= 0) call error_handler("int_send_2cff_shutdown: &
            &deallocate of xc_fitfct_index failed" )
    end if

  end subroutine fitfct_index_shutdown
  !*************************************************************
 

#if 0
  !*************************************************************
    subroutine print_integralfiles(flag,integral,filename)
      implicit none
      !------------ Declaration of formal parameters -----------
      character(LEN=*), intent(in) :: flag ! "ch" or "xc" or "mx"
      real(kind=r8_kind), intent(in) :: integral(:)
      character(len=*), intent(in), optional :: filename
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_rec, in_unit,i_ua1,i_l1, &
           i_ind1,i_exp1,i_ua2,i_l2,i_ind2,i_exp2,i_gc1,i_gc2, &
           n_tot_records,i,imax,u_g,i_meta
      real(kind=r8_kind) :: org_integral(ubound(integral,1))
      integer(kind=i4_kind),parameter :: IO_len = 1024 ! rec_len = IO_len + 1
      character(len=5) :: storage ! "ch", "xc", "ch_xc", or "xc_ch_xc"
      !------------ Executable code ----------------------------
      write(debug_unit,*)
      write(debug_unit,*)
      write(debug_unit,*)
      select case (flag)
      case ("ch")
         write(debug_unit,*) '#################### 2 center coulomb integrals #####################'
         n_tot_records = n_tot_records_ch
         storage = "ch"
      case ("xc")
         write(debug_unit,*) '#################### 2 center exchange integrals ####################'
         n_tot_records = n_tot_records_xc
         storage = "xc"
         write(debug_unit,*) '################# mixed 2 center overlap integrals ##################'
         n_tot_records = n_tot_records_mx
         storage = "ch_xc"
      end select

      if ( present(filename) ) then
         in_unit = openget_iounit(filename, &
              status='OLD', form='unformatted')
         imax = n_tot_records/(IO_len+1)
         u_g=1
         do i=1,imax
            read(in_unit,ERR=98) org_integral(u_g:u_g+IO_len)
            u_g=u_g+IO_len+1
         enddo
         read(in_unit,ERR=98) org_integral(u_g:n_tot_records)
         call returnclose_iounit(in_unit,status='keep')
      else
         org_integral = 0.0_r8_kind
      endif

      write(debug_unit,*) 'rec  meta u1 l1 e1 i1 u2 l2 e2 i2 g1 g2 new integral   org integral'
      do i_rec = 1, n_tot_records
         call expand_integral_index(storage,i_rec, &
                                    i_ua1,i_l1,i_exp1,i_ind1,i_gc1, &
                                    i_ua2,i_l2,i_exp2,i_ind2,i_gc2  )
!        call quadrupel_set(quadrupel,i_ua1,i_l1,i_ua2,i_l2) <<< removed
         if ( i_ind1 .gt. 0 ) then
            call loc_loc_index(storage,i_ua1,i_l1,i_exp1,i_ind1, &
                                       i_ua2,i_l2,i_meta)
         else
            if ( i_ind2 .gt. 0 ) then
               call glob_loc_index(storage,i_ua1,i_gc1,i_ua2,i_l2,i_meta)
            else
               call glob_glob_index(storage,i_ua1,i_gc1,i_ua2,i_meta)
            endif
         endif
         write(debug_unit,'(2I5,10I3,2F15.10)',ERR=88) i_rec, i_meta, &
              i_ua1,i_l1,i_exp1,i_ind1,i_ua2,i_l2,i_exp2,i_ind2,i_gc1,i_gc2, &
              integral(i_rec),org_integral(i_rec)
      enddo
      select case (flag)
      case ("ch")
         write(debug_unit,*) '#################### 2 center coulomb integrals end #####################'
      case ("xc")
         write(debug_unit,*) '#################### 2 center exchange integrals end ####################'
         write(debug_unit,*) '################# mixed 2 center overlap integrals end ##################'
      end select
      write(debug_unit,*)
      write(debug_unit,*)
      return
88    call write_to_output_units( "print_integralfiles: writing failed at record",i_rec)
      call error_handler("print_integralfiles: writing failed")
98    call write_to_output_units( "print_integralfiles: reading of original integral failed at record",i_rec)
      call error_handler("print_integralfiles: reading failed")
    end subroutine print_integralfiles
#endif

    subroutine write_integralfile(flag,filename,integral)
      implicit none
      !------------ Declaration of formal parameters -----------
      character(LEN=*), intent(in)   :: flag ! "ch" or "xc" or "mx"
      character(len=*), intent(in) :: filename
      real(kind=r8_kind), intent(in) :: integral(:)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_rec, out_unit, n_tot_records, &
           i, imax, u_g
      integer(kind=i4_kind),parameter :: IO_len = 1024 ! rec_len = IO_len + 1
      !------------ Executable code ----------------------------
      select case (flag)
      case ("ch")
         n_tot_records = n_tot_records_ch
      case ("xc")
         n_tot_records = n_tot_records_xc
      case ("mx")
         n_tot_records = n_tot_records_mx
      end select
      out_unit = openget_iounit(filename, &
           status='REPLACE', form='unformatted')
!     << all three kinds of files are now treated in the same way >>
         imax = n_tot_records/(IO_len+1)
         u_g=1
         do i=1,imax
            write(out_unit,ERR=88) integral(u_g:u_g+IO_len)
            u_g=u_g+IO_len+1
         enddo
         write(out_unit,ERR=88) integral(u_g:n_tot_records)
      call returnclose_iounit(out_unit,status='keep')
      return
88    call write_to_output_units( " write_integralfile: writing failed at record",i_rec)
      call error_handler(" write_integralfile: writing failed")
    end subroutine write_integralfile

#if 0
    subroutine testread_integralfile(flag,filename,integral)
      implicit none
      !------------ Declaration of formal parameters -----------
      character(LEN=*), intent(in)   :: flag ! "ch" or "xc" or "mx"
      character(len=*), intent(in) :: filename
      real(kind=r8_kind), intent(in) :: integral(:)
      !** End of interface *************************************
      !------------ Declaration of external functions ----------
      logical, external :: precision_check
       !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_rec,i_ua1,i_l1,i_ind1,i_exp1, &
           i_ua2,i_l2,i_ind2,i_exp2,in_unit,i_gc1,i_gc2, &
           n_tot_records,i,imax,u_g
      real(kind=r8_kind) :: org_integral(ubound(integral,1))
      integer(kind=i4_kind),parameter :: IO_len = 1024 ! rec_len = IO_len + 1
      character(len=5) :: storage ! "ch", "xc", "ch_xc", or "xc_ch"
      !------------ Executable code ----------------------------
      select case (flag)
      case ("ch")
         n_tot_records = n_tot_records_ch
         storage = "ch"
      case ("xc")
         n_tot_records = n_tot_records_xc
         storage = "ch"
      case ("mx")
         n_tot_records = n_tot_records_mx
         storage = "ch_xc"
      end select
      in_unit = openget_iounit(filename, &
           status='OLD', form='unformatted')
      imax = n_tot_records/(IO_len+1)
      u_g=1
      do i=1,imax
         read(in_unit,ERR=98) org_integral(u_g:u_g+IO_len)
         u_g=u_g+IO_len+1
      enddo
      read(in_unit,ERR=98) org_integral(u_g:n_tot_records)
      call returnclose_iounit(in_unit,status='delete')
      do i_rec = 1, n_tot_records
         if ( precision_check(org_integral(i_rec),integral(i_rec)) ) then
            call expand_integral_index(storage,i_rec, &
                                       i_ua1,i_l1,i_exp1,i_ind1,i_gc1, &
                                       i_ua2,i_l2,i_exp2,i_ind2,i_gc2  )
            call write_to_output_units( "testread_integralfile:&
                 & incorrect integral detected:")
            call write_to_output_units( "  i_rec  ",i_rec  )
            call write_to_output_units( "  i_ua1  ",i_ua1  )
            call write_to_output_units( "  i_l1   ",i_l1   )
            call write_to_output_units( "  i_exp1 ",i_exp1 )
            call write_to_output_units( "  i_ind1 ",i_ind1 )
            call write_to_output_units( "  i_ua2  ",i_ua2  )
            call write_to_output_units( "  i_l2   ",i_l2   )
            call write_to_output_units( "  i_exp2 ",i_exp2 )
            call write_to_output_units( "  i_ind2 ",i_ind2 )
            call write_to_output_units( "  i_gc1  ",i_gc1  )
            call write_to_output_units( "  i_gc2  ",i_gc2  )
            call write_to_output_units( "  original integral ", re=org_integral(i_rec))
            call write_to_output_units( "  new      integral ", re=integral(i_rec))
            call error_handler( "testread_integralfile:&
                 & incorrect integral detected:")
         endif
      enddo
      return
98    call write_to_output_units( "testread_integralfile: &
           &reading org integrals failed at record",i_rec)
      call error_handler("int_send_2cff_shutdown: reading failed")
    end subroutine testread_integralfile
  !*************************************************************
#endif
 

  !*************************************************************
  subroutine int_send_2cff_send
    !  Purpose: sending of contracted and symmetry adapted
    !    integrals from int_data_2cff_module.
    !    When packing the integrals, the mapping to the scf-
    !    metaindex is also done.
    !    Two center integrals are send to the master and
    !    three center integrals are distributed over all hosts
    !    splitting them over the fitfunction metaindex.
    !  called by: integral_calc_quad_2cff
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    implicit none
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    call start_timer(timer_int_send_2cff(integralpar_i_int_part))

    if ( comm_i_am_master() ) then
       if ( integralpar_2cch_no .and. need_ch_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: writing coulomb")
          call store_integral( "ch", quadrupel, loc_loc_int_ch, &
               loc_glob_int_ch, glob_loc_int_ch, glob_glob_int_ch, &
               loc_glob_contrib_ch, glob_loc_contrib_ch, integral_2cch_no )
       endif
       if ( integralpar_2cch_pre .and. need_ch_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: writing pre_coulomb")
          call store_integral( "ch", quadrupel, loc_loc_int_pre, &
               loc_glob_int_pre, glob_loc_int_pre, glob_glob_int_pre, &
               loc_glob_contrib_pre, glob_loc_contrib_pre, integral_2cch_pre )
       endif
       if ( integralpar_2cxc_no .and. need_xc_xc ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: writing exchange")
          call store_integral( "xc", quadrupel, loc_loc_int_xc, &
               loc_glob_int_xc, glob_loc_int_xc, glob_glob_int_xc, &
               loc_glob_contrib_xc, glob_loc_contrib_xc, integral_2cxc_no )
       endif
       if ( integralpar_2c_mixed .and. need_ch_xc ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: writing exchange")
          call store_integral( "ch_xc", quadrupel, loc_loc_int_ch_xc, &
               loc_glob_int_ch_xc, glob_loc_int_ch_xc, glob_glob_int_ch_xc, &
               loc_glob_contrib_ch_xc, glob_loc_contrib_ch_xc, &
               integral_2c_mixed )
       endif
       if ( integralpar_2c_mixed .and. need_xc_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: writing exchange")
          call store_integral( "xc_ch", quadrupel, loc_loc_int_xc_ch, &
               loc_glob_int_xc_ch, glob_loc_int_xc_ch, glob_glob_int_xc_ch, &
               loc_glob_contrib_xc_ch, glob_loc_contrib_xc_ch, &
               integral_2c_mixed )
       endif
       n_missing_quadrupels = n_missing_quadrupels - 1
    else
       call comm_init_send(comm_master_host,msgtag_int_2cff_result)
       call quadrupel_pack(quadrupel)
       if ( integralpar_2cch_no .and. need_ch_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: packing coulomb")
          call pack_submatrix( loc_loc_int_ch, &
               loc_glob_int_ch, glob_loc_int_ch, glob_glob_int_ch, &
               loc_glob_contrib_ch, glob_loc_contrib_ch )
       endif
       if ( integralpar_2cch_pre .and. need_ch_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: packing coulomb")
          call pack_submatrix( loc_loc_int_pre, &
               loc_glob_int_pre, glob_loc_int_pre, glob_glob_int_pre, &
               loc_glob_contrib_pre, glob_loc_contrib_pre )
       endif
       if ( integralpar_2cxc_no .and. need_xc_xc ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: packing exchange")
          call pack_submatrix( loc_loc_int_xc, &
               loc_glob_int_xc, glob_loc_int_xc, glob_glob_int_xc, &
               loc_glob_contrib_xc, glob_loc_contrib_xc )
       endif
       ! pack both the < f_k | g_l > and the < g_k | f_l > contributions
       if ( integralpar_2c_mixed .and. need_ch_xc ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: packing exchange")
          call pack_submatrix( loc_loc_int_ch_xc, &
               loc_glob_int_ch_xc, glob_loc_int_ch_xc, glob_glob_int_ch_xc, &
               loc_glob_contrib_ch_xc, glob_loc_contrib_ch_xc )
       endif
       if ( integralpar_2c_mixed .and. need_xc_ch ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cff_send: packing exchange")
          call pack_submatrix( loc_loc_int_xc_ch, &
               loc_glob_int_xc_ch, glob_loc_int_xc_ch, glob_glob_int_xc_ch, &
               loc_glob_contrib_xc_ch, glob_loc_contrib_xc_ch )
       endif
       call comm_send()
    endif

    call stop_timer(timer_int_send_2cff(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cff_send: done")

  end subroutine int_send_2cff_send


    subroutine store_integral( flag, quadrupel, loc_loc_int, loc_glob_int, &
         glob_loc_int, glob_glob_int, loc_glob_contrib, &
         glob_loc_contrib, integral, n_dim )
      implicit none
      !------------ Declaration of formal parameters -----------
      character(LEN=*), intent(in)     :: flag ! "ch", "xc", "ch_xc", or "xc_ch"
      type(quadrupel_type), intent(in) :: quadrupel
      ! for the description of the sense of the following
      ! parameters, see int_data_2cff_module
      real(kind=r8_kind), intent(in) :: loc_loc_int(:,:,:,:), &
           loc_glob_int(:,:,:), glob_loc_int(:,:,:), glob_glob_int(:,:)
      integer(kind=i4_kind), intent(in) :: loc_glob_contrib(:), glob_loc_contrib(:)
      real(kind=r8_kind), dimension(:), intent(out) :: integral
      integer(kind=i4_kind), optional, intent(in) :: n_dim(:)
      ! if provided n_dim holds (/n_c1,n_c2,n_if1,n_if2,n_gc1,n_gc2/)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      logical               :: triangle_mode, transposed
      character(LEN=5)      :: storage
      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,i_c2,nn_if2,nn_c2, &
           i_meta,n_c1,n_c2,i_gc1,i_gc2,n_gc1,n_gc2,nn_gc2,nn_gc1,n_if1, n_if2
      logical :: diagonal
      !------------ Executable code ----------------------------
      diagonal = quadrupel%ua1 == quadrupel%ua2 .and. &
                 quadrupel%l1  == quadrupel%l2

      select case (flag)
      case ("ch")
         if (.not.present(n_dim)) then
         n_c1 = n_ch_c1
         n_c2 = n_ch_c2
         n_gc1 = n_ch_gc1
         n_gc2 = n_ch_gc2
         endif
         triangle_mode = .true.
         transposed = .false.
         storage = "ch"
      case ("xc")
         if (.not.present(n_dim)) then
         n_c1 = n_xc_c1
         n_c2 = n_xc_c2
         n_gc1 = n_xc_gc1
         n_gc2 = n_xc_gc2
         endif
         triangle_mode = .true.
         transposed = .false.
         storage = "xc"
      case ("ch_xc")
         if (.not.present(n_dim)) then
            n_c1 = n_ch_c1
            n_c2 = n_xc_c2
            n_gc1 = n_ch_gc1
            n_gc2 = n_xc_gc2
         endif
         triangle_mode = .false.
         transposed = .false.
         storage = "ch_xc"
      case ("xc_ch")
         if (.not.present(n_dim)) then
            n_c1 = n_xc_c1
            n_c2 = n_ch_c2
            n_gc1 = n_xc_gc1
            n_gc2 = n_ch_gc2
         endif
         triangle_mode = .false.
         transposed = .true.
         storage = "ch_xc" ! because transposed = .true.
      case default
         call error_handler("store_integral: invalid flag")
      end select
      if (.not.present(n_dim)) then
         n_if1 = n_indep1
         n_if2 = n_indep2
      else
         n_c1  = n_dim(1)
         n_c2  = n_dim(2)
         n_if1 = n_dim(3)
         n_if2 = n_dim(4)
         n_gc1 = n_dim(5)
         n_gc2 = n_dim(6)
      endif
 
      loc_loc_block: if ( triangle_mode .or. .not.transposed ) then
      ! << store loc_loc_int(i_c2,i_c1,i_if2,i_if1) >>
      do i_c1 = 1, n_c1
         if ( triangle_mode .and. diagonal ) then
            nn_c2 = i_c1
         else
            nn_c2 = n_c2
         endif
         
         do i_if1 = 1, n_if1
            if(n_if2>0) then
               call loc_loc_index(storage,quadrupel%ua1,quadrupel%l1,i_c1,i_if1, &
                    quadrupel%ua2,quadrupel%l2,i_meta)
               do i_c2 = 1, nn_c2
                  if ( triangle_mode .and. diagonal .and. i_c2 == i_c1 ) then
                     nn_if2 = i_if1
                  else
                     nn_if2 = n_if2
                  endif
                  integral(i_meta:i_meta+nn_if2-1) = &
                       loc_loc_int(i_c2,i_c1,1:nn_if2,i_if1)
                  i_meta = i_meta + nn_if2
               end do
            end if
         enddo
      enddo
      elseif ( .not. diagonal ) then loc_loc_block
      ! << store transposed loc_loc_int(i_c2,i_c1,i_if2,i_if1) >>
      do i_c2 = 1, n_c2
         do i_if2 = 1, n_if2
            if(n_if1>0) then
               call loc_loc_index(storage,quadrupel%ua2,quadrupel%l2,i_c2,i_if2, &
                    quadrupel%ua1,quadrupel%l1,i_meta)
               do i_c1 = 1, n_c1
                  integral(i_meta:i_meta+n_if1-1) = &
                       loc_loc_int(i_c2,i_c1,i_if2,1:n_if1)
                  i_meta = i_meta + n_if1
               enddo
            end if
         enddo
      enddo
      endif loc_loc_block
      glob_loc_block: if ( triangle_mode .or. .not.transposed ) then
      ! << store glob_loc_int(i_c2,i_if2,i_gc1) for all i_gc1 which >>
      ! << contribute to the current value of quadrupel%l1          >>
      do i_gc1 = 1,n_gc1
         if (glob_loc_contrib(i_gc1) .gt. 0) then
            call glob_loc_index(storage,quadrupel%ua1,i_gc1, &
                                        quadrupel%ua2,quadrupel%l2,i_meta)
            do i_c2 = 1, n_c2
               integral(i_meta:i_meta+n_if2-1) = &
                    integral(i_meta:i_meta+n_if2-1) + &
                    glob_loc_int(i_c2,:,i_gc1)
               i_meta = i_meta + n_if2
            enddo
         endif
      enddo
      elseif ( .not. diagonal ) then glob_loc_block
      ! << store transposed glob_loc_int(i_c2,i_if2,i_gc1) for all i_gc1 >>
      ! << which contribute to the current value of quadrupel%l1         >>
      do i_c2 = 1, n_c2
         do i_if2 = 1, n_if2
            if ( any( glob_loc_contrib(1:n_gc1) > 0 ) ) then
               call loc_glob_index(storage,quadrupel%ua2,quadrupel%l2,i_c2, &
                                     i_if2,quadrupel%ua1,i_meta)
               integral(i_meta:i_meta+n_gc1-1) = &
                   integral(i_meta:i_meta+n_gc1-1) + &
                   glob_loc_int(i_c2,i_if2,:)
            endif
         enddo
      enddo
      endif glob_loc_block
      loc_glob_block: if ( .not. triangle_mode .and. .not. transposed ) then
         ! << store loc_glob_int(i_gc2,i_c1,i_if1) for all i_gc2 which >> 
         ! << contribute to the current value of quadrupel%l2          >> 
         do i_c1 = 1, n_c1
            do i_if1 = 1, n_if1
               if ( any( loc_glob_contrib(1:n_gc2) > 0 ) ) then
                  call loc_glob_index(storage,quadrupel%ua1,quadrupel%l1,i_c1, &
                                        i_if1,quadrupel%ua2,i_meta)
                  integral(i_meta:i_meta+n_gc2-1) = & 
                       integral(i_meta:i_meta+n_gc2-1) + & 
                       loc_glob_int(:,i_c1,i_if1)
               endif
            enddo
         enddo
      elseif ( .not. diagonal ) then loc_glob_block
         ! << store transposed loc_glob_int(i_gc2,i_c1,i_if1) for all i_gc2 >>
         ! << which contribute to the current value of quadrupel%l2         >>
         do i_gc2 = 1,n_gc2
            if (loc_glob_contrib(i_gc2) .gt. 0) then
               call glob_loc_index(storage,quadrupel%ua2,i_gc2, &
                                           quadrupel%ua1,quadrupel%l1,i_meta)
               do i_c1 = 1, n_c1
                  integral(i_meta:i_meta+n_if1-1) = &
                       integral(i_meta:i_meta+n_if1-1) + &
                       loc_glob_int(i_gc2,i_c1,:)
                  i_meta = i_meta + n_if1
               enddo
            endif
         enddo
      endif loc_glob_block

      glob_glob_block_1: if ( triangle_mode .or. .not.transposed ) then
      ! << store glob_glob_int(i_gc2,i_gc1) for all (i_gc1,i_gc2) which   >>
      ! << contribute to the current value of (quadrupel%l1,quadrupel%l2) >>
      do i_gc1 = 1,n_gc1
         if (glob_loc_contrib(i_gc1) .gt. 0) then
            if ( triangle_mode .and. quadrupel%ua1 .eq. quadrupel%ua2) then
               nn_gc2 = i_gc1
            else 
               nn_gc2 = n_gc2
            endif
            if ( any( loc_glob_contrib(1:nn_gc2) > 0 ) ) then
               call glob_glob_index(storage,quadrupel%ua1,i_gc1, &
                                            quadrupel%ua2,i_meta)
               integral(i_meta:i_meta+nn_gc2-1) = &
                    integral(i_meta:i_meta+nn_gc2-1) + &
                    glob_glob_int(1:nn_gc2,i_gc1)
            end if
         end if
      enddo
      endif glob_glob_block_1

      ! note that triangle_mode and transposed are never set simultaneously !
      glob_glob_block_2: if ( .not. diagonal .and. ( transposed .or. &
           ( triangle_mode .and. quadrupel%ua1 .eq. quadrupel%ua2 ) ) ) then
         ! << store transposed glob_glob_int(i_gc2,i_gc1) as well >>
         do i_gc2=1,n_gc2
            if ( triangle_mode ) then ! ua1 == ua2
               nn_gc1 = i_gc2
            else 
               nn_gc1 = n_gc1
            endif
            if (loc_glob_contrib(i_gc2) .gt. 0) then
               if ( any( glob_loc_contrib(1:nn_gc1) > 0 ) ) then
                  call glob_glob_index(storage,quadrupel%ua2,i_gc2, &
                                               quadrupel%ua1,i_meta)
                  integral(i_meta:i_meta+nn_gc1-1) = &
                       integral(i_meta:i_meta+nn_gc1-1) + &
                       glob_glob_int(i_gc2,1:nn_gc1)
               end if
            end if
         end do
      endif glob_glob_block_2

    end subroutine store_integral
    !*************************************************************
 
 
    !*************************************************************
    subroutine pack_submatrix( loc_loc_int,  loc_glob_int , &
                               glob_loc_int, glob_glob_int, &
                               loc_glob_contrib, glob_loc_contrib )
      implicit none
      !------------ Declaration of formal parameters -----------
      real(kind=r8_kind), intent(in) :: loc_loc_int(:,:,:,:), &
                                        loc_glob_int(:,:,:), &
                                        glob_loc_int(:,:,:), &
                                        glob_glob_int(:,:)
      integer(kind=i4_kind), intent(in) :: glob_loc_contrib(:), &
                                           loc_glob_contrib(:)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: n_c1,n_c2,n_if1,n_if2,n_gc1,n_gc2,n_size,info
      !------------ Executable code ----------------------------
 
      ! << pack necessary dimensions >>
      n_c2  = size(loc_loc_int,1)
      n_c1  = size(loc_loc_int,2)
      n_if2 = size(loc_loc_int,3)
      n_if1 = size(loc_loc_int,4)
      n_gc2 = size(glob_glob_int,1)
      n_gc1 = size(glob_glob_int,2)
      call commpack((/n_c1,n_c2,n_if1,n_if2,n_gc1,n_gc2/),6,1,info)
      if (info /= 0) call error_handler( "pack_submatrix: &
           &packing of (n_c1,n_c2,n_if1,n_if2,n_gc1,n_gc2) failed")

      ! << pack global contraction selection arrays >>
      n_size = n_gc2
      if (n_size > 0) then
         call commpack(loc_glob_contrib,n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of loc_glob_contrib failed")
      endif
      n_size = n_gc1
      if (n_size > 0) then
         call commpack(glob_loc_contrib,n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of glob_loc_contrib failed")
      endif
 
      ! << pack integral contributions >>
      n_size = n_c2*n_c1*n_if2*n_if1
      if (n_size > 0) then
         call commpack(loc_loc_int(1,1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of loc_loc_int failed")
      endif
      n_size = n_c2*n_if2*n_gc1
      if (n_size > 0) then
         call commpack(glob_loc_int(1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of glob_loc_int failed")
      endif
      n_size = n_gc2*n_c1*n_if1
      if (n_size > 0) then
         call commpack(loc_glob_int(1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of loc_glob_int failed")
      endif
      n_size = n_gc2*n_gc1
      if (n_size > 0) then
         call commpack(glob_glob_int(1,1),n_size,1,info)
         if (info /= 0) call error_handler( "pack_submatrix: &
              &packing of glob_glob_int failed")
      endif
 
    end subroutine pack_submatrix
  !*************************************************************
 

  !*************************************************************
  subroutine int_send_2cff_receive
    !  Purpose:
    !    + receiving of contracted and symmetry adapted
    !      integrals for one quadrupel
    !    + Storing them in direct access files in scf symmetric
    !      storage mode with the record number as metindex
    !  called by: main_slave, integral_main_2cff, int_send_2cff_shutdown
    !** End of interface *****************************************
    implicit none
    !------------ Declaration of local variables ---------------------
    type(quadrupel_type)            :: quadrupel
    logical                         :: restart_timer,need_ch_ch,&
                                       need_ch_xc,need_xc_ch,need_xc_xc
    !------------ Executable code ------------------------------------

    if ( timer_int_idle_2cff(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_idle_2cff(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_receive_2cff(integralpar_i_int_part))

    ! get the quadrupel associated with the sent data
    call quadrupel_unpack(quadrupel)

    if ( output_int_loops ) then
       write(output_unit,*) &
            "int_send_2cff_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
       write(stdout_unit,*) &
            "int_send_2cff_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
    endif

    need_ch_ch = quadrupel%l1 <= unique_atoms(quadrupel%ua1)%lmax_ch .and. &
                 quadrupel%l2 <= unique_atoms(quadrupel%ua2)%lmax_ch
    need_ch_xc = quadrupel%l1 <= unique_atoms(quadrupel%ua1)%lmax_ch .and. &
                 quadrupel%l2 <= unique_atoms(quadrupel%ua2)%lmax_xc
    need_xc_ch = quadrupel%l1 <= unique_atoms(quadrupel%ua1)%lmax_xc .and. &
                 quadrupel%l2 <= unique_atoms(quadrupel%ua2)%lmax_ch
    need_xc_xc = quadrupel%l1 <= unique_atoms(quadrupel%ua1)%lmax_xc .and. &
                 quadrupel%l2 <= unique_atoms(quadrupel%ua2)%lmax_xc

    if ( integralpar_2cch_no .and. need_ch_ch ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cff_receive: unpacking coulomb")
       call unpack_submatrix( "ch", quadrupel, integral_2cch_no )
    endif
    if ( integralpar_2cch_pre .and. need_ch_ch ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cff_receive: unpacking prescreening")
       call unpack_submatrix( "ch", quadrupel, integral_2cch_pre )
    endif

    if ( integralpar_2cxc_no .and. need_xc_xc ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cff_receive: unpacking exchange")
       call unpack_submatrix( "xc", quadrupel, integral_2cxc_no )
    endif
    ! both < f_k | g_l > and < g_k | f_l > contributions are unpacked 
    if ( integralpar_2c_mixed .and. need_ch_xc ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cff_receive: unpacking charge/exchange")
       call unpack_submatrix( "ch_xc", quadrupel, integral_2c_mixed )
    endif
    if ( integralpar_2c_mixed .and. need_xc_ch ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cff_receive: unpacking exchange/charge")
       call unpack_submatrix( "xc_ch", quadrupel, integral_2c_mixed )
    endif


    n_missing_quadrupels = n_missing_quadrupels - 1


    call stop_timer(timer_int_receive_2cff(integralpar_i_int_part))
    if (restart_timer) call start_timer(timer_int_idle_2cff(integralpar_i_int_part))
   
    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cff_receive: done")

  end subroutine int_send_2cff_receive
  !*************************************************************
 
 
    !*************************************************************
    subroutine unpack_submatrix( flag, quadrupel, integral )
      implicit none
      !------------ Declaration of formal parameters -----------
      character(LEN=*), intent(in) :: flag ! "ch", "xc", "ch_xc", or "xc_ch"
      type(quadrupel_type), intent(in) :: quadrupel
      real(kind=r8_kind), dimension(:), intent(out) :: integral
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), allocatable :: loc_loc_int(:,:,:,:), &
                                         glob_loc_int(:,:,:), &
                                         loc_glob_int(:,:,:), &
                                         glob_glob_int(:,:)
      integer(kind=i4_kind), allocatable :: loc_glob_contrib(:), &
                                            glob_loc_contrib(:)
      integer(kind=i4_kind) :: n_c1,n_c2,n_if1,n_if2,n_gc1,n_gc2,n_dim(6), &
                               n_size,info,status
      !------------ Executable code ----------------------------
 
      ! << unpack required dimensions >>
      call communpack(n_dim,6,1,info)
      if (info /= 0) call error_handler( "unpack_submatrix: &
           &unpacking of (n_c1,n_c2,n_gc1,n_gc2) failed")
      n_c1  = n_dim(1)
      n_c2  = n_dim(2)
      n_if1 = n_dim(3)
      n_if2 = n_dim(4)
      n_gc1 = n_dim(5)
      n_gc2 = n_dim(6)
 
      ! << unpack global contraction selection arrays >>
      allocate(loc_glob_contrib(n_gc2),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocation of loc_glob_contrib failed")
      n_size = n_gc2
      if (n_size > 0) then
         call communpack(loc_glob_contrib,n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking of loc_glob_contrib failed")
      endif
 
      allocate(glob_loc_contrib(n_gc1),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocaton of glob_loc_contrib failed")
      n_size = n_gc1
      if (n_size > 0) then
         call communpack(glob_loc_contrib,n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking glob_loc_contrib failed")
      endif

      ! << unpack integral contributions >>
      allocate(loc_loc_int(n_c2,n_c1,n_if2,n_if1),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocaton of loc_loc_int failed")
      n_size = n_c2*n_c1*n_if2*n_if1
      if (n_size > 0) then
         call communpack(loc_loc_int(1,1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking loc_loc_int failed")
      endif
      allocate(glob_loc_int(n_c2,n_if2,n_gc1),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocaton of glob_loc_int failed")
      n_size = n_c2*n_if2*n_gc1
      if (n_size > 0) then
         call communpack(glob_loc_int(1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking glob_loc_int failed")
      endif
      allocate(loc_glob_int(n_gc2,n_c1,n_if1),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocaton of loc_glob_int failed")
      n_size = n_gc2*n_c1*n_if1
      if (n_size > 0) then
         call communpack(loc_glob_int(1,1,1),n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking loc_glob_int failed")
      endif
      allocate(glob_glob_int(n_gc2,n_gc1),stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &allocaton of glob_glob_int failed")
      n_size = n_gc2*n_gc1
      if (n_size > 0) then
         call communpack(glob_glob_int(1,1),n_size,1,info)
         if (info /= 0) call error_handler( "unpack_submatrix: &
              &unpacking glob_glob_int failed")
      endif

      ! << store integral contributions in the final array >>
      call store_integral( flag, quadrupel, loc_loc_int, loc_glob_int, &
           glob_loc_int, glob_glob_int, loc_glob_contrib, glob_loc_contrib, &
           integral, n_dim )
 
      ! << deallocate local arrays >>
      deallocate(loc_loc_int ,glob_loc_int , &
                 loc_glob_int,glob_glob_int, stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &deallocaton of xxx_yyy_int failed")
      deallocate(loc_glob_contrib,glob_loc_contrib,stat=status)
      if (status /= 0) call error_handler( "unpack_submatrix: &
           &deallocation of xxx_yyy_contrib failed")
 
    end subroutine unpack_submatrix
  !*************************************************************


  !*************************************************************
  subroutine loc_loc_index(flag,i_ua1,l1,i_exp1,i_ind1,i_ua2,l2,i_rec,n_rec)
    ! Purpose: returns the position of the first element of 
    !          loc_loc_int(i_ua1,l1,i_exp1,i_ind1,i_ua2,l2,:,:)
    !          within the linear storage mode of the 2-center integrals
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*), intent(in) :: flag ! "ch", "xc", "xc_ch", or "ch_xc"
    integer(kind=i4_kind), intent(in)  :: i_ua1,l1,i_ind1,i_exp1,i_ua2,l2
    integer(kind=i4_kind), intent(out) :: i_rec
    integer(kind=i4_kind), intent(out), optional :: n_rec
    !** End of interface ***************************************
    integer(kind=i4_kind) :: i, j, n
    
    ! first get the fitfct indices
    select case (flag)
    case ("ch")
!!$       print*,shape(ch_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1))
!!$       print*,i_exp1,i_ind1
       ! there is error in processing of the fit on DEC when r**2 functions are absent
       ! as seen with prints above.

       i = ch_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       if(size(ch_fitfct_index(i_ua2)%loc_index(l2)%index,1) == 0) then
          j=0
       else
          j = ch_fitfct_index(i_ua2)%loc_index(l2)%index(     1,     1)
       end if
       n = ch_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("xc")
       i = xc_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = xc_fitfct_index(i_ua2)%loc_index(l2)%index(     1,     1)
       n = xc_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("xc_ch")
       i = xc_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = ch_fitfct_index(i_ua2)%loc_index(l2)%index(     1,     1)
       n = ch_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("ch_xc")
       i = ch_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = xc_fitfct_index(i_ua2)%loc_index(l2)%index(     1,     1)
       n = xc_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    end select
 
    ! then load the integral index
    select case (flag)
    case ("ch","xc")
       i_rec = (i*(i-1))/2 + j
    case ("xc_ch")
       i_rec = (i-1)*n_ch + j
    case ("ch_xc")
       i_rec = (i-1)*n_xc + j
    end select
    if (present(n_rec)) n_rec = n
 
  end subroutine loc_loc_index
  !*************************************************************


  !*************************************************************
  subroutine glob_loc_index(flag,i_ua1,i_gc1,i_ua2,l2,i_rec,n_rec)
    ! Purpose: returns the position of the first element of 
    !          glob_loc_int(i_ua1,i_gc1,i_ua2,l2,:,:)
    !          within the linear storage mode of the 2-center integrals
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*), intent(in) :: flag ! "ch", "xc", "xc_ch", or "ch_xc"
    integer(kind=i4_kind), intent(in)  :: i_ua1,i_gc1,i_ua2,l2
    integer(kind=i4_kind), intent(out) :: i_rec
    integer(kind=i4_kind), intent(out), optional :: n_rec
    !** End of interface ***************************************
    integer(kind=i4_kind) :: i, j, n
    
    ! first get the fitfct indices
    select case (flag)
    case ("ch")
       i = ch_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = ch_fitfct_index(i_ua2)%loc_index(l2)%index(1,1)
       n = ch_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("xc")
       i = xc_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = xc_fitfct_index(i_ua2)%loc_index(l2)%index(1,1)
       n = xc_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("xc_ch")
       i = xc_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = ch_fitfct_index(i_ua2)%loc_index(l2)%index(1,1)
       n = ch_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    case ("ch_xc")
       i = ch_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = xc_fitfct_index(i_ua2)%loc_index(l2)%index(1,1)
       n = xc_fitfct_index(i_ua2)%loc_index(l2)%N_lc
    end select
 
    ! then load the integral index
    select case (flag)
    case ("ch","xc")
       i_rec = (i*(i-1))/2 + j
    case ("xc_ch")
       i_rec = (i-1)*n_ch + j
    case ("ch_xc")
       i_rec = (i-1)*n_xc + j
    end select
    if (present(n_rec)) n_rec = n
 
  end subroutine glob_loc_index
  !*************************************************************


  !*************************************************************
  subroutine loc_glob_index(flag,i_ua1,l1,i_exp1,i_ind1,i_ua2,i_rec,n_rec)
    ! Purpose: returns the position of the first element of 
    !          loc_glob_int(i_ua1,l1,i_exp1,i_ind1,i_ua2,:) transposed
    !          within the linear storage mode of the 2-center integrals
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*), intent(in) :: flag ! "ch", "xc", "xc_ch", or "ch_xc"
    integer(kind=i4_kind), intent(in)  :: i_ua1,l1,i_exp1,i_ind1,i_ua2
    integer(kind=i4_kind), intent(out) :: i_rec
    integer(kind=i4_kind), intent(out), optional :: n_rec
    !** End of interface ***************************************
    integer(kind=i4_kind) :: i, j, n
    
    ! first get the fitfct indices
    select case (flag)
    case ("ch")
       i = ch_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = ch_fitfct_index(i_ua2)%glob_index(1)
       n = ch_fitfct_index(i_ua2)%N_gc
    case ("xc")
       i = xc_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = xc_fitfct_index(i_ua2)%glob_index(1)
       n = xc_fitfct_index(i_ua2)%N_gc
    case ("xc_ch")
       i = ch_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = xc_fitfct_index(i_ua2)%glob_index(1)
       n = ch_fitfct_index(i_ua2)%N_gc
    case ("ch_xc")
       i = xc_fitfct_index(i_ua1)%loc_index(l1)%index(i_exp1,i_ind1)
       j = ch_fitfct_index(i_ua2)%glob_index(1)
       n = xc_fitfct_index(i_ua2)%N_gc
    end select
 
    ! then load the integral index
    select case (flag)
    case ("ch","xc")
       i_rec = (i*(i-1))/2 + j
    case ("xc_ch")
       i_rec = (i-1)*n_ch + j
    case ("ch_xc")
       i_rec = (i-1)*n_xc + j
    end select
    if (present(n_rec)) n_rec = n
 
  end subroutine loc_glob_index
  !*************************************************************


  !*************************************************************
  subroutine glob_glob_index(flag,i_ua1,i_gc1,i_ua2,i_rec,n_rec)
    ! Purpose: returns the position of the first element of 
    !          glob_glob_int(i_ua1,i_gc1,i_ua2,:)
    !          within the linear storage mode of the 2-center integrals
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*), intent(in) :: flag ! "ch", "xc", "xc_ch", or "ch_xc"
    integer(kind=i4_kind), intent(in)  :: i_ua1,i_gc1,i_ua2
    integer(kind=i4_kind), intent(out) :: i_rec
    integer(kind=i4_kind), intent(out), optional :: n_rec
    !** End of interface ***************************************
    integer(kind=i4_kind) :: i, j, n
    
    ! first get the fitfct indeces
    select case (flag)
    case ("ch")
       i = ch_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = ch_fitfct_index(i_ua2)%glob_index(1)
       n = ch_fitfct_index(i_ua2)%N_gc
    case ("xc")
       i = xc_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = xc_fitfct_index(i_ua2)%glob_index(1)
       n = xc_fitfct_index(i_ua2)%N_gc
    case ("xc_ch")
       i = xc_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = ch_fitfct_index(i_ua2)%glob_index(1)
       n = ch_fitfct_index(i_ua2)%N_gc
    case ("ch_xc")
       i = ch_fitfct_index(i_ua1)%glob_index(i_gc1)
       j = xc_fitfct_index(i_ua2)%glob_index(1)
       n = xc_fitfct_index(i_ua2)%N_gc
    end select

    ! then load the integral index
    select case (flag)
    case ("ch","xc")
       i_rec = (i*(i-1))/2 + j
    case ("xc_ch")
       i_rec = (i-1)*n_ch + j
    case ("ch_xc")
       i_rec = (i-1)*n_xc + j
    end select
    if (present(n_rec)) n_rec = n

  end subroutine glob_glob_index
  !*************************************************************


#if 0
  !*************************************************************
  subroutine expand_integral_index(flag,i_rec, &
                                   i_ua1,i_bas1,i_exp1,i_ind1,i_gc1, &
                                   i_ua2,i_bas2,i_exp2,i_ind2,i_gc2)
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*)     , intent(in)  :: flag ! "ch", "xc", "ch_xc", "xc_ch"
    integer(kind=i4_kind), intent(in)  :: i_rec
    integer(kind=i4_kind), intent(out) :: i_ua1,i_bas1,i_exp1,i_ind1,i_gc1, &
                                          i_ua2,i_bas2,i_exp2,i_ind2,i_gc2
    !** End of interface ***************************************
    integer(kind=i4_kind) :: i, j
    character(len=2)      :: flag1, flag2
    !------------ Executable code ------------------------------

    ! first find the fitfct indices
    if ( i_rec <= 0 ) call error_handler( &
         "expand_integral_index: non-positiv integral index provided" )
    select case (flag)
    case ("ch","xc") ! i_rec = (i*(i-1))/2 + j
       ! estimate first fitfct index i by (approximative) real arithmetic
       i = int(sqrt(real(2*i_rec,r8_kind)-1.75_r8_kind)+0.5_r8_kind,i4_kind)
       j = i_rec - (i*(i-1))/2
       ! correct first fitfct index, if necessary
       do while ( j <= 0 )
          i = i - 1
          j = j + i
       end do
       do while ( j > i )
         j = j - i
         i = i + 1
       end do
       flag1 = flag
       flag2 = flag
    case ("xc_ch") ! i_rec = (i-1)*n_ch + j
       i = (i_rec-1)/n_ch + 1
       j = i_rec - (i-1)*n_ch
       flag1 = "xc"
       flag2 = "ch"
    case ("ch_xc") ! i_rec = (i-1)*n_xc + j
       i = (i_rec-1)/n_xc + 1
       j = i_rec - (i-1)*n_xc
       flag1 = "ch"
       flag2 = "xc"
    end select

    ! then expand the individual fitfct indices
    call expand_fitfct_index(flag1,i,i_ua1,i_bas1,i_exp1,i_ind1,i_gc1)
    call expand_fitfct_index(flag2,j,i_ua2,i_bas2,i_exp2,i_ind2,i_gc2)

  contains

  subroutine expand_fitfct_index(type,index,i_ua,i_bas,i_exp,i_ind,i_gc)
    implicit none
    !------------ Declaration of formal parameters -------------
    character(LEN=*)     , intent(in)  :: type ! "ch", "xc"
    integer(kind=i4_kind), intent(in)  :: index
    integer(kind=i4_kind), intent(out) :: i_ua,i_bas,i_exp,i_ind,i_gc
    !** End of interface ***************************************
    integer(kind=i4_kind)            :: i_rec, n_rec, n_ind, n_lc, n_gc
    type(fitfct_index_type), pointer :: fitfct_index(:)
    type(loc_index_type)   , pointer :: loc_index(:)
    !------------ Executable code ------------------------------

    select case(type)
    case ("ch")
       fitfct_index => ch_fitfct_index
    case ("xc")
       fitfct_index => xc_fitfct_index
    end select
    n_rec = 0
    do i_ua=1,N_unique_atoms
       loc_index => fitfct_index(i_ua)%loc_index
       do i_bas=-1,size(loc_index)-2
          n_lc = loc_index(i_bas)%N_lc
          n_rec = n_rec + n_lc
          if (index <= n_rec) then
             ! fitfct index referes to one of the current local contractions
             i_rec = index - ( n_rec - n_lc )
             n_ind = size(loc_index(i_bas)%index,2)
             i_exp = (i_rec-1)/n_ind + 1
             i_ind = i_rec - (i_exp-1)*n_ind
             i_gc  = 0
             return
          end if
       end do
    end do

    ! fitfct index referes to a global contraction
    do i_ua=1,N_unique_atoms
       n_gc = fitfct_index(i_ua)%N_gc
       n_rec = n_rec + n_gc
       if (index <= n_rec) then
          i_bas = -2
          i_exp =  0 
          i_ind =  0
          i_gc  = index - ( n_rec - n_gc )
          return
       end if
    end do

    call error_handler("expand_fitfct_index: Wrong fitfct index provided")

  end subroutine expand_fitfct_index

  end subroutine expand_integral_index
  !*************************************************************
#endif

!--------------- End of module ----------------------------------
end module int_send_2cff_module

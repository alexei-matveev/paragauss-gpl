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
module  int_send_dipole_module
!---------------------------------------------------------------
!
!  Purpose: Sending, receiving and storing of dipole integrals.
!           The results of one quadrupel stored in
!           int_data_dipole_module are send from a slave to the
!           master and there stored in dipole_module.
!
!  Module called by: main_slave, integral_main_dipole,
!    integral_calc_quad_dipole
!
!  References: Publisher Document: Concepts of Integral Part
!
!
!  Author: TB
!  Date: 7/97
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   6/98
! Description: ...
!
! Modification put if(options_kinematic_factors) in some places
! Author: DG
! Date:   21/03/2001
! Description: I need to pass matrix elements in uncontracted basis
! in case of kinematic factors becouse they should be multiplyed to kf later
! and after that the contractation must be done.
! Dipole moments presently are contracted before commpack/uncommpack procedure
!
! Modification Master/Slave concept to DLB
! Author: AN
! Date:   4/11
! Description: for scheduling the dipole integrals DLB is used
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
use comm_module
use msgtag_module
use symmetry_data_module
use unique_atom_module
use dipole_module, only: dipole_integral_type, dipole_integrals, dipole_unused
#ifdef WITH_GTENSOR
use gtensor_module
use hfc_module
#endif
use quadrupel_module
use int_data_dipole_module
use orbitalprojection_module
use timer_module
use time_module
use options_module

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
public int_send_dipole_setup, int_send_dipole_shutdown, &
       int_send_dipole_send, int_send_dipole_receive
  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Declaration of private types ----

type, private :: index_ind1_exp1_ua2_type
   integer(kind=i4_kind), pointer :: l2(:) ! (i_l2)
end type index_ind1_exp1_ua2_type

type, private :: index_l1_type
   type(index_ind1_exp1_ua2_type), pointer :: ind1_exp1_ua2(:,:,:)
   ! (i_ind1,i_exp1,i_ua2)
end type index_l1_type

type, private :: index_ua1_type
   type(index_l1_type), pointer :: l1(:) ! (i_l1)
end type index_ua1_type

type, private :: index_ir_type
   type(index_ua1_type), pointer :: ua1(:) ! (i_ua1)
end type index_ir_type

! data type for spin-orbit

type, private :: index_ind1_exp1_ua2_type_p
   integer(kind=i4_kind), pointer :: l2(:) ! (i_l2)
end type index_ind1_exp1_ua2_type_p

type, private :: index_l1_type_p
   type(index_ind1_exp1_ua2_type_p), pointer :: ind1_exp1_ua2(:,:,:)
   ! (i_ind1,i_exp1,i_ua2)
end type index_l1_type_p

type, private :: index_ua1_type_p
   type(index_l1_type_p), pointer :: l1(:) ! (i_l1)
end type index_ua1_type_p

type, private :: index_ir_type_p
   type(index_ua1_type_p), pointer :: ua1(:) ! (i_ua1)
end type index_ir_type_p

!------------ Declaration of constants and variables ----
integer(kind=i4_kind) :: n_missing_quadrupels, &
     first_host, last_host, my_hostindex, n_file, &
     n_hosts_to_report_back, blocklength
type(index_ir_type), allocatable :: &
     metaindex_of_first_integral(:) ! (i_ir)
type(index_ir_type_p), allocatable :: &
     metaindex_of_first_integral_p(:),& ! (i_ir)
     metaindex_of_first_integral_k(:)
     ! to look up metaindex of first integral and number
     ! n_records of following records for a given vector
     ! (i_ir,i_ua1,i_l1,i_ind1,i_exp2,i_ua2,i_l2)
!------------ Subroutines ---------------------------------------
contains
  !*************************************************************
  subroutine int_send_dipole_setup(n_quads)
    !  Purpose:
    !   + calculates help variables number of quadrupels
    !     to be received and total number of record in file
    !   + calculates help array for mapping integrals to
    !     scf storage format.
    ! called by: integral_setup_dipole (only on master)
    !** End of interface *****************************************
    implicit none
    integer(kind=i4_kind) :: n_quads
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: status, i_ua1, i_ua2, &
         i_l1, n_l2, i_meta, i_ir, i_ind1, i_exp1, i_l2, &
         n_ind2,i_ind2, n_exp2, n_records
    type(unique_atom_type),             pointer :: ua1,ua2
    type(unique_atom_basis_type),       pointer :: uab1,uab2
    type(unique_atom_partner_type),     pointer :: uap1,uap2
    type(index_ind1_exp1_ua2_type),     pointer :: ind1_exp1_ua2(:,:,:)
    type(index_ind1_exp1_ua2_type_p),   pointer :: ind1_exp1_ua2_p(:,:,:)
    integer,                            pointer :: metaindex(:)
    integer,                            pointer :: metaindex_p(:)
    logical           :: diagonal,diagonal_unique
    !------------ Executable code ------------------------------------
    DPRINT '###int_send_dipole_setup:ENTERED'
    ! fetch my hostindex from comm_module
    my_hostindex = comm_myindex()
    ! calculate n_quadrupels
    DPRINT "integralpar_offdiag_dipoles",integralpar_offdiag_dipoles
    n_missing_quadrupels = n_quads
    if ( output_int_progress ) then
       call write_to_output_units("")
       call write_to_output_units( &
            "total number of dipole quadrupels: ",n_quads)
       call write_to_output_units("")
    endif
    ! during shutdown procedure, all slaves have to report back
    n_hosts_to_report_back = comm_get_n_processors() - 1

    if (options_spin_orbit .and. (.not. options_kinematic_factors)) then
       !
       ! SPIN ORBIT
       !
       ! calculate index_of_ua_l
       allocate( metaindex_of_first_integral_p(symmetry_data_n_proj_irreps()), &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral failed" )
       do i_ir = 1,symmetry_data_n_proj_irreps()
          allocate( metaindex_of_first_integral_p(i_ir)%ua1(N_unique_atoms), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ua1 failed" )
          i_meta = 1
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             allocate( metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l1 failed" )
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)
                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                allocate( ind1_exp1_ua2_p(uap1%N_independent_fcts, &
                     uab1%N_uncontracted_fcts + uab1%N_contracted_fcts, i_ua1), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2_p
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                      do i_ua2 = 1, i_ua1
                         ua2 => unique_atoms(i_ua2)
                         if ( i_ua1 .eq. i_ua2 ) then
                            n_l2 = i_l1
                            diagonal_unique = .true.
                         else
                            n_l2 = ua2%lmax_ob
                            diagonal_unique = .false.
                         endif
                         allocate(metaindex_p(0:n_l2), &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l2 failed" )
                         ind1_exp1_ua2_p(i_ind1,i_exp1,i_ua2)%l2 => metaindex_p
                         do i_l2 = 0, n_l2
                            diagonal = (diagonal_unique.and.(i_l1.eq.i_l2))
                            uab2 => ua2%l_ob(i_l2)
                            uap2 => ua2%symadapt_spor_partner(i_ir,i_l2)
                            if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2  ) then
                               n_ind2 = i_ind1
                            else
                               n_ind2 = uap2%N_independent_fcts
                            endif
                            n_records = 0
                            do i_ind2 = 1, n_ind2
                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                    .and. i_ind1 .eq. i_ind2  ) then
                                  n_exp2 = i_exp1
                               else
                                  n_exp2 = uab2%N_uncontracted_fcts + &
                                       uab2%N_contracted_fcts
                               endif
                               n_records = n_records + n_exp2
                               ! factor 2 accounts for real and imaginary parts
                            enddo
                            metaindex_p(i_l2) = i_meta
                            i_meta = i_meta + n_records
                         enddo! i_l2
                      enddo! i_ua2
                   enddo! i_exp1
                enddo! i_ind1
             enddo! i_l1
          enddo! i_ua1
       enddo! i_irrep
    end if
       !---------------------------------------------------------------------------------------------------
       kinemat: if (options_kinematic_factors .and. options_spin_orbit) then
          !
          ! SPIN ORBIT
          ! CASE FOR KINEMATIC FACTORS = ON
          ! calculate index_of_ua_l
          allocate( metaindex_of_first_integral_p(symmetry_data_n_proj_irreps()), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral failed" )
          do i_ir = 1,symmetry_data_n_proj_irreps()
             allocate( metaindex_of_first_integral_p(i_ir)%ua1(N_unique_atoms), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ua1 failed" )
             i_meta = 1
             do i_ua1 = 1, N_unique_atoms
                ua1 => unique_atoms(i_ua1)
                allocate( metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l1 failed" )
                do i_l1 = 0, ua1%lmax_ob
                   uab1 => ua1%l_ob(i_l1)
                   uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                   allocate( ind1_exp1_ua2_p(uap1%N_independent_fcts, &
                        uab1%N_exponents, i_ua1), &
                        stat=status )
                   if ( status .ne. 0 ) call error_handler( &
                        "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                   metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2_p
                   do i_ind1 = 1, uap1%N_independent_fcts
                      do i_exp1 = 1, uab1%N_exponents
                         do i_ua2 = 1, i_ua1
                            ua2 => unique_atoms(i_ua2)
                            if ( i_ua1 .eq. i_ua2 ) then
                               n_l2 = i_l1
                               diagonal_unique = .true.
                            else
                               n_l2 = ua2%lmax_ob
                               diagonal_unique = .false.
                            endif
                            allocate( metaindex_p(0:n_l2), &
                                 stat=status )
                            if ( status .ne. 0 ) call error_handler( &
                                 "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l2 failed" )
                            ind1_exp1_ua2_p(i_ind1,i_exp1,i_ua2)%l2 => metaindex_p
                            do i_l2 = 0, n_l2
                               diagonal = (diagonal_unique.and.(i_l1.eq.i_l2))
                               uab2 => ua2%l_ob(i_l2)
                               uap2 => ua2%symadapt_spor_partner(i_ir,i_l2)
                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2  ) then
                                  n_ind2 = i_ind1
                               else
                                  n_ind2 = uap2%N_independent_fcts
                               endif
                               n_records = 0
                               do i_ind2 = 1, n_ind2
                                  if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                       .and. i_ind1 .eq. i_ind2 ) then
                                     n_exp2 = i_exp1
                                  else
                                     n_exp2 = uab2%N_exponents
                                  endif
                                  n_records = n_records + n_exp2
                                  ! factor 2 accounts for real and imaginary parts
                               enddo
                               metaindex_p(i_l2) = i_meta
                               i_meta = i_meta + n_records
                            enddo! i_l2
                         enddo! i_ua2
                      enddo! i_exp1
                   enddo! i_ind1
                enddo! i_l1
             enddo! i_ua1
          enddo! i_irrep
       endif kinemat
!-------------------------------------------
        kfnoso: if (options_kinematic_factors .and. (.not. options_spin_orbit)) then
          !
          ! NO SPIN ORBIT
          ! CASE FOR KINEMATIC FACTORS = ON
          ! calculate index_of_ua_l
          allocate( metaindex_of_first_integral(symmetry_data_n_irreps()), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_setup: allocate of metaindex_of_first_integral failed" )
          do i_ir = 1,symmetry_data_n_irreps()
             allocate( metaindex_of_first_integral(i_ir)%ua1(N_unique_atoms), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_setup: allocate of metaindex_of_first_integral ua1 failed" )
             i_meta = 1
             do i_ua1 = 1, N_unique_atoms
                ua1 => unique_atoms(i_ua1)
                allocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_setup: allocate of metaindex_of_first_integral l1 failed" )
                do i_l1 = 0, ua1%lmax_ob
                   uab1 => ua1%l_ob(i_l1)
                   uap1 => ua1%symadapt_partner(i_ir,i_l1)
                   allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                        uab1%N_exponents, i_ua1), &
                        stat=status )
                   if ( status .ne. 0 ) call error_handler( &
                        "int_send_2cob3c_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                   metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2
                   do i_ind1 = 1, uap1%N_independent_fcts
                      do i_exp1 = 1, uab1%N_exponents
                         do i_ua2 = 1, i_ua1
                            ua2 => unique_atoms(i_ua2)
                            if ( i_ua1 .eq. i_ua2 ) then
                               n_l2 = i_l1
                               diagonal_unique = .true.
                            else
                               n_l2 = ua2%lmax_ob
                               diagonal_unique = .false.
                            endif
                            allocate( metaindex(0:n_l2), &
                                 stat=status )
                            if ( status .ne. 0 ) call error_handler( &
                                 "int_send_2cob3c_setup: allocate of metaindex_of_first_integral l2 failed" )
                            ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2 => metaindex
                            do i_l2 = 0, n_l2
                               diagonal = (diagonal_unique.and.(i_l1.eq.i_l2))
                               uab2 => ua2%l_ob(i_l2)
                               uap2 => ua2%symadapt_partner(i_ir,i_l2)
                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2  ) then
                                  n_ind2 = i_ind1
                               else
                                  n_ind2 = uap2%N_independent_fcts
                               endif
                               n_records = 0
                               do i_ind2 = 1, n_ind2
                                  if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                       .and. i_ind1 .eq. i_ind2 ) then
                                     n_exp2 = i_exp1
                                  else
                                     n_exp2 = uab2%N_exponents
                                  endif
                                  n_records = n_records + n_exp2
                                  ! factor 2 accounts for real and imaginary parts
                               enddo
                               metaindex(i_l2) = i_meta
                               i_meta = i_meta + n_records
                            enddo! i_l2
                         enddo! i_ua2
                      enddo! i_exp1
                   enddo! i_ind1
                enddo! i_l1
             enddo! i_ua1
          enddo! i_irrep
       endif kfnoso
       !---------------------------------------------------------------------------------------------------------------------------------
      nosonokf: if ((.not. options_spin_orbit) .and. (.not. options_kinematic_factors) ) then
       ! calculate index_of_ua_l
       allocate( metaindex_of_first_integral(symmetry_data_n_irreps()), &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_dipole_setup: allocate of metaindex_of_first_integral failed" )
       do i_ir = 1,symmetry_data_n_irreps()
          allocate( metaindex_of_first_integral(i_ir)%ua1(N_unique_atoms), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_dipole_setup: allocate of metaindex_of_first_integral ua1 failed" )
          i_meta = 1
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             allocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_dipole_setup: allocate of metaindex_of_first_integral l1 failed" )
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)
                uap1 => ua1%symadapt_partner(i_ir,i_l1)
                allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                     uab1%N_uncontracted_fcts + uab1%N_contracted_fcts, i_ua1), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_dipole_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                      do i_ua2 = 1, i_ua1
                         ua2 => unique_atoms(i_ua2)
                         if ( i_ua1 .eq. i_ua2 ) then
                            n_l2 = i_l1
                         else
                            n_l2 = ua2%lmax_ob
                         endif
                         allocate( metaindex(0:n_l2), &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_dipole_setup: allocate of metaindex_of_first_integral l2 failed" )
                         ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2 => metaindex
                         do i_l2 = 0, n_l2
                            uab2 => ua2%l_ob(i_l2)
                            uap2 => ua2%symadapt_partner(i_ir,i_l2)
                            if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 ) then
                               n_ind2 = i_ind1
                            else
                               n_ind2 = uap2%N_independent_fcts
                            endif
                            n_records = 0
                            do i_ind2 = 1, n_ind2
                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                    .and. i_ind1 .eq. i_ind2 ) then
                                  n_exp2 = i_exp1
                               else
                                  n_exp2 = uab2%N_uncontracted_fcts + &
                                       uab2%N_contracted_fcts
                               endif
                               n_records = n_records + n_exp2
                            enddo
                            metaindex(i_l2) = i_meta
                            i_meta = i_meta + n_records
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif nosonokf
    DPRINT '###int_send_dipole_setup:DONE'
  end subroutine int_send_dipole_setup
  !*************************************************************
  subroutine int_send_dipole_shutdown
    !  Purpose:
    !    + waiting for missing integrals
    !    + deallocation of help arrays
    ! called by: integral_main_dipole (only on master)
    !** End of interface *****************************************
    implicit none
    integer(kind=i4_kind) :: status, i_ir, i_ua1, i_ua2, i_l1, &
         i_ind1, i_exp1
    type(unique_atom_type),             pointer :: ua1
    type(unique_atom_basis_type),       pointer :: uab1
    type(unique_atom_partner_type),     pointer :: uap1
    type(index_ind1_exp1_ua2_type),     pointer :: ind1_exp1_ua2(:,:,:)
    type(index_ind1_exp1_ua2_type_p),   pointer :: ind1_exp1_ua2_p(:,:,:)
    !------------ Executable code ------------------------------------

    ! waiting for missing integrals

    do while (n_missing_quadrupels .gt. 0 )
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_dipole_shutdown: waiting for quadrupels: ",n_missing_quadrupels)
       call comm_save_recv(comm_all_other_hosts, msgtag_int_dipole_result)
          if ( output_slaveoperations ) &
               call write_to_output_units("int_send_dipole_shutdown: integral_dipole_result")
          call int_send_dipole_receive()
          !:VN receives both dipole and dipoleg matrix elements
    enddo

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i_ir = 1,symmetry_data_n_proj_irreps()
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)
                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                ind1_exp1_ua2_p => metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                      do i_ua2 = 1, i_ua1
                         deallocate( ind1_exp1_ua2_p(i_ind1,i_exp1,i_ua2)%l2, &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_dipole_shutdown: options_spin_orbit:ind1_exp1_ua2_p:&
                              & deallocate of metaindex_of_first_integral l2 failed" )
                      enddo
                   enddo
                enddo
                deallocate( ind1_exp1_ua2_p, stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_dipole_shutdown:options_spin_orbit:&
                     & deallocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
             enddo! i_l1
             deallocate( metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1, stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_dipole_shutdown:options_spin_orbit: deallocate of metaindex_of_first_integral l1 failed" )
          enddo! i_ua1
          deallocate( metaindex_of_first_integral_p(i_ir)%ua1, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_dipole_shutdown:options_spin_orbit: deallocate of metaindex_of_first_integral ua1 failed" )
       enddo
       deallocate( metaindex_of_first_integral_p, &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_dipole_shutdown: options_spin_orbit:deallocate of metaindex_of_first_integral failed" )
       !-----------------------------------------------------------------
       ! SPIN ORBIT
       ! CASE KINEMATIC FACTORS = ON
!!$       kinematic: if (options_spin_orbit .and. .false. ) then
!!$          do i_ir = 1,symmetry_data_n_proj_irreps()
!!$             do i_ua1 = 1, N_unique_atoms
!!$                ua1 => unique_atoms(i_ua1)
!!$                do i_l1 = 0, ua1%lmax_ob
!!$                   uab1 => ua1%l_ob(i_l1)
!!$                   uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
!!$                   ind1_exp1_ua2_p => metaindex_of_first_integral_p(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
!!$                   do i_ind1 = 1, uap1%N_independent_fcts
!!$                      do i_exp1 = 1, uab1%N_exponents
!!$                         do i_ua2 = 1, i_ua1
!!$                            deallocate( ind1_exp1_ua2_p(i_ind1,i_exp1,i_ua2)%l2, &
!!$                                 stat=status )
!!$                            if ( status .ne. 0 ) call error_handler( &
!!$                                 "int_send_dipole_shutdown:options_kinematic_factors:&
!!$                                 & deallocate of metaindex_of_first_integral_k l2 failed" )
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                   deallocate( ind1_exp1_ua2_k, stat=status )
!!$                   if ( status .ne. 0 ) call error_handler( &
!!$                        "int_send_dipole_shutdown:options_kinematic_factors:&
!!$                        & deallocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
!!$                enddo! i_l1
!!$                deallocate( metaindex_of_first_integral_k(i_ir)%ua1(i_ua1)%l1, stat=status )
!!$                if ( status .ne. 0 ) call error_handler( &
!!$                     "int_send_dipole_shutdown:options_kinematic_factors: allocate of metaindex_of_first_integral l1 failed" )
!!$             enddo! i_ua1
!!$             deallocate( metaindex_of_first_integral_k(i_ir)%ua1, &
!!$                  stat=status )
!!$             if ( status .ne. 0 ) call error_handler( &
!!$                  "int_send_dipole_shutdown:options_kinematic_factors: deallocate of metaindex_of_first_integral ua1 failed" )
!!$          enddo
!!$          deallocate( metaindex_of_first_integral_k, &
!!$               stat=status )
!!$          if ( status .ne. 0 ) call error_handler( &
!!$               "int_send_dipole_shutdown:options_kinematic_factors: deallocate of metaindex_of_first_integral failed" )
!!$       end if kinematic
       !----------------------------------------------------------------
    else
       ! deallocation of help arrays
       do i_ir = 1,symmetry_data_n_irreps()
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)
                uap1 => ua1%symadapt_partner(i_ir,i_l1)
                ind1_exp1_ua2 => metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                      do i_ua2 = 1, i_ua1
                         deallocate( ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2, &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_dipole_shutdown:else: deallocate of metaindex_of_first_integral l2 failed" )
                      enddo
                   enddo
                enddo
                deallocate( ind1_exp1_ua2, stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_dipole_shutdown:else: deallocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
             enddo
             deallocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1, stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_dipole_shutdown:else: allocate of metaindex_of_first_integral l1 failed" )
          enddo
          deallocate( metaindex_of_first_integral(i_ir)%ua1, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_dipole_shutdown:else: deallocate of metaindex_of_first_integral ua1 failed" )
       enddo
       deallocate( metaindex_of_first_integral, &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_dipole_shutdown:else: deallocate of metaindex_of_first_integral failed" )
    endif

    if ( output_int_loops ) call write_to_output_units("int_send_dipole_shutdown: done")

    DPRINT 'int_send_dipole_shutdown:DONE'

  end subroutine int_send_dipole_shutdown
  !*************************************************************
  subroutine int_send_dipole_send
    !  Purpose: sending of contracted and symmetry adapted
    !    integrals from int_data_dipole_module to the master.
    !    When packing the integrals, the mapping to the scf-
    !    metaindex is also done. If subroutine is called by
    !    master, integrals are directley stored in dipole_module.
    !  called by: integral_calc_quad_dipole
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use operations_module, only: operations_gtensor, operations_dipole, operations_hfc
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: info
    !------------ Executable code ------------------------------------
    DPRINT 'entered int_send_dipole_send ...'
    call start_timer(timer_int_send_2cob3c(integralpar_i_int_part))
   master: if ( comm_i_am_master() ) then
       if ( integralpar_2cob_dipole ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_dipole_send: storing dipole")
          call store_integralfile_dipole()
          !:VN this subroutines stores presently  both dipole and dipoleg matrix elements
       endif
       n_missing_quadrupels = n_missing_quadrupels - 1
    else ! I am slave

       call start_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
       call comm_init_send(comm_master_host,msgtag_int_dipole_result)
       call quadrupel_pack(quadrupel)

       if(options_kinematic_factors) then !DG I should pass uncontracted dimensions in this case
          call commpack(n_exp1,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: packing n_exp1  failed")
          call commpack(n_exp2,info)
          if (info .ne. 0) call error_handler( &
               "int_send_dipole_send: packing n_exp2 failed")

       end if

          call commpack(n_c1,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: packing n_c1 failed")
          call commpack(n_c2,info)
          if (info .ne. 0) call error_handler( &
               "int_send_dipole_send: packing n_c2 failed")


       if ( integralpar_2cob_dipole ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_dipole_send: packing dipole")
          call pack_integralfile_dipole()
          !VN: this subroutines packs presently both dipole and dipoleg matrix elements
       endif
       call stop_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
       call start_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
       call comm_send()
       call stop_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
    endif master

    call stop_timer(timer_int_send_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_dipole_send: done")
    DPRINT 'int_send_dipole_send: done'
  contains

    subroutine store_integralfile_dipole()

      implicit none
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
           saint_2cob_dipole,saint_2cob_dipole_real,saint_2cob_dipole_imag
      real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
           saint_2cob_dipoleg_real,saint_2cob_dipoleg_imag
       real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
           saint_2cob_hfc_real,saint_2cob_hfc_imag, saint_2cob_hfc
      integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
           i_ip1, i_ip2, n_if1, n_if2, i_if1, i_if2, nn_c2, i_c1, &
           i_meta, i_meta1, i_meta2, diag_offdiag, n_dim1, n_dim2 , i_ua
      type(dipole_integral_type), pointer :: dipole
      type(dipole_integral_type), pointer :: dipoleg
      type(dipole_integral_type), pointer :: hfc
      logical                         :: diagonal_spin
      !------------ Executable code ----------------------------
      DPRINT 'entered store_integralfile_dipole...'

      SO: if (options_spin_orbit) then
         !
         ! SPIN ORBIT
         !
        op_dp: if (operations_dipole) then
         xyz_so: do i_xyz = 1, 3
            diagonal_spin = diagonal
            i_ip1 = 1
            irrep1_so: do i_ir1 = 1, symmetry_data_n_proj_irreps()
               n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner1_so: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                  i_ip2 = 1
                  irrep2_so: do i_ir2 = 1, i_ir1
                     partner2_so: do i_pa2 = 1, symmetry_data_n_partners_proj(i_ir2)

                        needed_so: if ( associated( &
                             symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int ) ) then
                           saint_2cob_dipole_real => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int
                           saint_2cob_dipole_imag => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int_imag

                           dipole => dipole_integrals(i_xyz,i_ip2,i_ip1)

                           diagonal_block_so: if (i_ip1 .eq. i_ip2) then
                              !
                              ! diagonal block with triangular storage

                              qudrupel_required_so: if ( qudrupel_required_for_diagonal(quadrupel) ) then

                                 ind_fct_1_so: do i_if1 = 1, n_if1
                                    if ( diagonal_spin ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1_so: do i_c1 = 1, n_c1
                                       i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1)

                                       ind_fct_2_so: do i_if2 = 1, n_if2
                                          if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_c2
                                          endif

                                          if ((i_meta+nn_c2-1).gt.size(dipole%diagonal)) then
                                             DPRINT 'int_send_dipole_module:store_integralfile_dipole():'
                                             DPRINT '     :overflow found'
                                             call abort()
                                          endif
                                          dipole%diagonal(i_meta:i_meta+nn_c2-1) = &
                                               saint_2cob_dipole_real(1:nn_c2,i_c1,i_if2,i_if1)
                                          dipole%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                               saint_2cob_dipole_imag(1:nn_c2,i_c1,i_if2,i_if1)

                                          i_meta = i_meta + nn_c2
                                       enddo ind_fct_2_so
                                    enddo contraction1_so
                                 enddo ind_fct_1_so
                              endif qudrupel_required_so

                           else
                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                              ind_fct_1_od_so: do i_if1 = 1, n_if1
                                 contraction1_od_so: do i_c1 = 1, n_c1
                                    i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                    ind_fct_2_od_so: do i_if2 = 1, n_if2

                                       dipole%offdiagonal(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                            saint_2cob_dipole_real(1:n_c2,i_c1,i_if2,i_if1)

                                       dipole%offdiagonal_imag(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                            saint_2cob_dipole_imag(1:n_c2,i_c1,i_if2,i_if1)

                                       i_meta2 = i_meta2 + n_c2
                                    enddo ind_fct_2_od_so
                                    i_meta1 = i_meta1 + 1
                                 enddo contraction1_od_so
                              enddo ind_fct_1_od_so
                           endif diagonal_block_so
                        endif needed_so
                        i_ip2 = i_ip2 + 1
                     enddo partner2_so
                  enddo irrep2_so
                  i_ip1 = i_ip1 + 1
               enddo partner1_so
            enddo irrep1_so
         enddo xyz_so
      end if op_dp

#ifdef WITH_GTENSOR
         gten: if(operations_gtensor) then
            if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
              n_dim1 = n_exp1
              n_dim2 = n_exp2
              else
              n_dim1 = n_c1
              n_dim2 = n_c2
           end if
            xyz_sog: do i_xyz = 1, 10 !7
               diagonal_spin = diagonal
               diag_offdiag_cycle:  do diag_offdiag = 1,2
               i_ip1 = 1
               irrep1_sog: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                  n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                  partner1_sog: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)

                     i_ip2 = i_ip1
                     i_ir2 = i_ir1
                     i_pa2 = i_pa1
                           needed_sog: if ( associated( &
                                symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int ) ) then

                              saint_2cob_dipoleg_real => symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int
                              saint_2cob_dipoleg_imag => symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int_imag
                              !dipoleg => dipoleg_integrals(i_xyz,diag_offdiag,i_ip1)
                              if(diag_offdiag==1) dipoleg => gten_integrals_diag%integrals(i_xyz,i_ip1,1)
                              if(diag_offdiag==2) dipoleg => gten_integrals_offdiag%integrals(i_xyz,i_ip1,1)
                              diagonal_block_sog: if (diag_offdiag == 1) then

                                 ! diagonal block with triangular storage
                                 qudrupel_required_sog: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                    ind_fct_1_sog: do i_if1 = 1, n_if1
                                       if ( diagonal_spin ) then
                                          n_if2 = i_if1
                                       else
                                          n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                               N_independent_fcts
                                       endif
                                       contraction1_sog: do i_c1 = 1, n_dim1 !F1

                                                                                        !metaindexes (diff contract
                                             i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1) !and uncontr size)


                                          ind_fct_2_sog: do i_if2 = 1, n_if2
                                             if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                nn_c2 = i_c1
                                             else
                                                nn_c2 = n_dim2
                                             endif

                                             if ((i_meta+nn_c2-1).gt.size(dipoleg%diagonal)) then
                                                DPRINT 'overflow found'
                                                call abort()
                                             endif
                                             dipoleg%diagonal(i_meta:i_meta+nn_c2-1) = &
                                                  saint_2cob_dipoleg_real(1:nn_c2,i_c1,i_if2,i_if1)
                                             dipoleg%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                                  saint_2cob_dipoleg_imag(1:nn_c2,i_c1,i_if2,i_if1)
                                             i_meta = i_meta + nn_c2
                                          enddo ind_fct_2_sog
                                       enddo contraction1_sog
                                    enddo ind_fct_1_sog
                                 endif qudrupel_required_sog
                              else
                                 ! offdiagonal block with full storage
                                 n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                 if(options_kinematic_factors) then
                                         i_meta1 = orbitalprojection_spor_ob_k(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      else
                                         i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      end if
                                 ind_fct_1_od_sog: do i_if1 = 1, n_if1
                                    contraction1_od_sog: do i_c1 = 1, n_dim1

                                       if(options_kinematic_factors) then
                                          i_meta2 = orbitalprojection_spor_ob_k(i_ir2,quadrupel%l2,quadrupel%ua2)
                                       else
                                          i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                       end if
                                       ind_fct_2_od_sog: do i_if2 = 1, n_if2

                                          dipoleg%offdiagonal(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                               saint_2cob_dipoleg_real(1:n_dim2,i_c1,i_if2,i_if1)

                                          dipoleg%offdiagonal_imag(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                               saint_2cob_dipoleg_imag(1:n_dim2,i_c1,i_if2,i_if1)

                                          i_meta2 = i_meta2 + n_dim2
                                       enddo ind_fct_2_od_sog
                                       i_meta1 = i_meta1 + 1
                                    enddo contraction1_od_sog
                                 enddo ind_fct_1_od_sog
                              endif diagonal_block_sog
                           endif needed_sog
                     i_ip1 = i_ip1 + 1
                  enddo partner1_sog
               enddo irrep1_sog
            end do diag_offdiag_cycle
         enddo xyz_sog
      endif gten

       hfcop: if(operations_hfc) then
            if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
              n_dim1 = n_exp1
              n_dim2 = n_exp2
              else
              n_dim1 = n_c1
              n_dim2 = n_c2
           end if
           ualoop1: do i_ua=1, n_unique_atoms
            xyz_hfc1: do i_xyz = 1,9! dipole_hyperfine_xyz, L_xyz/r^3, f(r)* sigma_xyz
               diagonal_spin = diagonal
               diag_offdiag_cycle1:  do diag_offdiag = 1,2
               i_ip1 = 1
               irrep1_hfc1: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                  n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                  partner1_hfc1: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)

                     i_ip2 = i_ip1
                     i_ir2 = i_ir1
                     i_pa2 = i_pa1
                           needed_hfc1: if ( associated( &
                                symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int ) ) then

                              saint_2cob_hfc_real => symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int
                              saint_2cob_hfc_imag => symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int_imag
                             ! hfc => hfc_integrals(i_xyz,diag_offdiag,i_ip1,i_ua)
                              if(diag_offdiag == 1) hfc => hfc_integrals_diag%integrals(i_xyz,i_ip1,i_ua)
                              if(diag_offdiag == 2) hfc => hfc_integrals_offdiag%integrals(i_xyz,i_ip1,i_ua)
                              diagonal_block_hfc1: if (diag_offdiag == 1) then

                                 ! diagonal block with triangular storage
                                 qudrupel_required_hfc1: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                    ind_fct_1_hfc1: do i_if1 = 1, n_if1
                                       if ( diagonal_spin ) then
                                          n_if2 = i_if1
                                       else
                                          n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                               N_independent_fcts
                                       endif
                                       contraction1_hfc1: do i_c1 = 1, n_dim1 !F1


                                             i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1) !and uncontr size)


                                          ind_fct_2_hfc1: do i_if2 = 1, n_if2
                                             if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                nn_c2 = i_c1
                                             else
                                                nn_c2 = n_dim2
                                             endif

                                             if ((i_meta+nn_c2-1).gt.size(hfc%diagonal)) then
                                                DPRINT "overflow found"
                                                call abort() ! mpi_abort or mpi_abort_?
                                             endif
                                             hfc%diagonal(i_meta:i_meta+nn_c2-1) = &
                                                  saint_2cob_hfc_real(1:nn_c2,i_c1,i_if2,i_if1)
                                             hfc%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                                  saint_2cob_hfc_imag(1:nn_c2,i_c1,i_if2,i_if1)
                                             i_meta = i_meta + nn_c2
                                          enddo ind_fct_2_hfc1
                                       enddo contraction1_hfc1
                                    enddo ind_fct_1_hfc1
                                 endif qudrupel_required_hfc1
                              else
                                 ! offdiagonal block with full storage
                                 n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                 if(options_kinematic_factors) then
                                         i_meta1 = orbitalprojection_spor_ob_k(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      else
                                         i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      end if
                                 ind_fct_1_od_hfc1: do i_if1 = 1, n_if1
                                    contraction1_od_hfc1: do i_c1 = 1, n_dim1

                                       if(options_kinematic_factors) then
                                          i_meta2 = orbitalprojection_spor_ob_k(i_ir2,quadrupel%l2,quadrupel%ua2)
                                       else
                                          i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                       end if
                                       ind_fct_2_od_hfc1: do i_if2 = 1, n_if2

                                          hfc%offdiagonal(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                               saint_2cob_hfc_real(1:n_dim2,i_c1,i_if2,i_if1)

                                          hfc%offdiagonal_imag(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                               saint_2cob_hfc_imag(1:n_dim2,i_c1,i_if2,i_if1)

                                          i_meta2 = i_meta2 + n_dim2
                                       enddo ind_fct_2_od_hfc1
                                       i_meta1 = i_meta1 + 1
                                    enddo contraction1_od_hfc1
                                 enddo ind_fct_1_od_hfc1
                              endif diagonal_block_hfc1
                           endif needed_hfc1
                     i_ip1 = i_ip1 + 1
                  enddo partner1_hfc1
               enddo irrep1_hfc1
            end do diag_offdiag_cycle1
         enddo xyz_hfc1
      end do ualoop1
   endif hfcop
#endif

   else ! not options_spin_orbit
      if (operations_dipole) then
         xyz: do i_xyz = 1, 3
            i_ip1 = 1
            irrep1: do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                  i_ip2 = 1
                  irrep2: do i_ir2 = 1, i_ir1
                     partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                        needed: if ( associated( &
                             symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int ) ) then
                           saint_2cob_dipole => symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int
                           dipole => dipole_integrals(i_xyz,i_ip2,i_ip1)
                           diagonal_block: if (i_ip1 .eq. i_ip2) then
                              ! diagonal block with triangular storage
                              qudrupel_required: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_1: do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1: do i_c1 = 1, n_c1
                                       i_meta = metaindex(quadrupel,i_ir1,i_if1,i_c1)
                                       ind_fct_2: do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_c2
                                          endif
                                          dipole%diagonal(i_meta:i_meta+nn_c2-1) = &
                                               saint_2cob_dipole(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_meta = i_meta + nn_c2
                                       enddo ind_fct_2
                                    enddo contraction1
                                 enddo ind_fct_1
                              endif qudrupel_required
                           else
                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              i_meta1 = orbitalprojection_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                              ind_fct_1_od: do i_if1 = 1, n_if1
                                 contraction1_od: do i_c1 = 1, n_c1
                                    i_meta2 = orbitalprojection_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                    ind_fct_2_od: do i_if2 = 1, n_if2
                                       dipole%offdiagonal(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                            saint_2cob_dipole(1:n_c2,i_c1,i_if2,i_if1)
                                       i_meta2 = i_meta2 + n_c2
                                    enddo ind_fct_2_od
                                    i_meta1 = i_meta1 + 1
                                 enddo contraction1_od
                              enddo ind_fct_1_od
                           endif diagonal_block
                        endif needed
                        i_ip2 = i_ip2 + 1
                     enddo partner2
                  enddo irrep2
                  i_ip1 = i_ip1 + 1
               enddo partner1
            enddo irrep1
         enddo xyz

      end if

#ifdef WITH_GTENSOR
    hfcnr:if (operations_hfc) then
         if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
            n_dim1 = n_exp1
            n_dim2 = n_exp2
         else
            n_dim1 = n_c1
            n_dim2 = n_c2
         end if
         do i_ua = 1, n_unique_atoms
         xyznr: do i_xyz = 1, 7
            i_ip1 = 1
            irrep1nr: do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner1nr: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                     i_ip2 = i_ip1
                     i_ir2 = i_ir1
                     i_pa2 = i_pa1
                        needednr: if ( associated( &
                             symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int ) ) then
                           saint_2cob_hfc => symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int
                           hfc => hfc_integrals_diag%integrals(i_xyz,i_ip1,i_ua) !indexes!!!

                              ! diagonal block with triangular storage
                               if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                  do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    do i_c1 = 1, n_dim1
                                       i_meta = metaindex(quadrupel,i_ir1,i_if1,i_c1)
                                        do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_dim2
                                          endif
                                          hfc%diagonal(i_meta:i_meta+nn_c2-1) = &
                                               saint_2cob_hfc(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_meta = i_meta + nn_c2
                                       enddo
                                    enddo
                                 enddo
                              endif
                           end if needednr
                           i_ip1 = i_ip1 + 1
                        end do partner1nr
                     end do irrep1nr
                  end do xyznr

               end do

            end if hfcnr
#endif

         endif SO
      DPRINT 'store_integralfile_dipole: done'
    end subroutine store_integralfile_dipole

    subroutine pack_integralfile_dipole()

      implicit none
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
           saint_2cob_dipole, saint_2cob_dipole_real, saint_2cob_dipole_imag, saint_2cob_hfc
      real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
           saint_2cob_dipoleg_real, saint_2cob_dipoleg_imag,&
           saint_2cob_hfc_real, saint_2cob_hfc_imag
      integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
           i_ip1, i_ip2, n_if1, n_if2, i_if1, i_if2, &
           nn_c2, i_c1, status, info, i_buf_lower, i_buf_upper,&
           i_meta1, i_meta2, n_dim1, n_dim2, buf_size, diag_offdiag, i_ua
      integer(kind=i4_kind) :: n_pack
      real(kind=r8_kind), allocatable :: buffer(:)
      logical                         :: diagonal_spin
      !------------ Executable code ----------------------------

      if (options_spin_orbit) then
         !
         ! SPIN ORBIT
         !
         nopgt: if(operations_dipole) then
            pack_xyz_so: do i_xyz = 1, 3
               diagonal_spin = diagonal
               i_ip1 = 1
               irrep1_so: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                  n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                  partner1_so: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                     i_ip2 = 1
                     irrep2_so: do i_ir2 = 1, i_ir1
                        partner2_so: do i_pa2 = 1, symmetry_data_n_partners_proj(i_ir2)

                           needed_so: if ( associated( &
                                symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int ) ) then
                              saint_2cob_dipole_real => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int
                              saint_2cob_dipole_imag => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int_imag
                              allocate(buffer(2*size(saint_2cob_dipole_real)),stat=status)

                              if (status .ne. 0) call error_handler( &
                                   "int_send_dipole_send: allocating at packing failed")
                              i_buf_lower = 1
                              i_buf_upper = 0

                              diagonal_block_so: if (i_ip1 .eq. i_ip2) then

                                 ! diagonal block with triangular storage
                                 qudrupel_required_so: if ( qudrupel_required_for_diagonal(quadrupel) ) then

                                    ind_fct_1_so: do i_if1 = 1, n_if1
                                       if ( diagonal_spin ) then
                                          n_if2 = i_if1
                                       else
                                          n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                               N_independent_fcts
                                       endif
                                       contraction1_so: do i_c1 = 1, n_c1
                                          ind_fct_2_so: do i_if2 = 1, n_if2
                                             if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                nn_c2 = i_c1
                                             else
                                                nn_c2 = n_c2
                                             endif
                                             i_buf_upper = i_buf_upper + nn_c2

                                             if (i_buf_upper.gt.size(buffer))then
                                                DPRINT ' buffer overflow'

                                             endif
                                             buffer(i_buf_lower:i_buf_upper) = &
                                                  saint_2cob_dipole_real(1:nn_c2,i_c1,i_if2,i_if1)
                                             i_buf_lower = i_buf_lower + nn_c2
                                             i_buf_upper = i_buf_upper + nn_c2
                                             buffer(i_buf_lower:i_buf_upper) = &
                                                  saint_2cob_dipole_imag(1:nn_c2,i_c1,i_if2,i_if1)
                                             i_buf_lower = i_buf_lower + nn_c2
                                          enddo ind_fct_2_so
                                       enddo contraction1_so
                                    enddo ind_fct_1_so
                                 endif qudrupel_required_so
                              else
                                 ! offdiagonal block with full storage
                                 n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                 i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                 ind_fct_1_od_so: do i_if1 = 1, n_if1
                                    contraction1_od_so: do i_c1 = 1, n_c1
                                       i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                       ind_fct_2_od_so: do i_if2 = 1, n_if2
                                          i_buf_upper = i_buf_upper + n_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_dipole_real(1:n_c2,i_c1,i_if2,i_if1)

                                          i_buf_lower = i_buf_lower + n_c2
                                          i_buf_upper = i_buf_upper + n_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_dipole_imag(1:n_c2,i_c1,i_if2,i_if1)

                                          i_buf_lower = i_buf_lower + n_c2
                                          i_meta2 = i_meta2 + n_c2
                                       enddo ind_fct_2_od_so
                                       i_meta1 = i_meta1 + 1
                                    enddo contraction1_od_so
                                 enddo ind_fct_1_od_so
                              endif diagonal_block_so

                              call commpack(i_buf_upper,info)
                              if (info .ne. 0) call error_handler( &
                                   "int_send_dipole_send: packing buffer size failed")
                              call commpack(buffer,i_buf_upper,1,info)
                              if (info .ne. 0) call error_handler( &
                                   "int_send_dipole_send: packing buffer failed")
                              if (size(buffer).ge.10) then

                              endif

                              deallocate(buffer,stat=status)
                              if (status .ne. 0) call error_handler( &
                                   "int_send_dipole_send: deallocating at packing failed")

                           endif needed_so
                           i_ip2 = i_ip2 + 1
                        enddo partner2_so
                     enddo irrep2_so
                     i_ip1 = i_ip1 + 1
                  enddo partner1_so
               enddo irrep1_so
            enddo pack_xyz_so
         endif nopgt

#ifdef WITH_GTENSOR
         gten: if(operations_gtensor) then
            n_pack=0
            if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
               n_dim1 = n_exp1
               n_dim2 = n_exp2
            else
               n_dim1 = n_c1
               n_dim2 = n_c2
            end if

            xyz_sog: do i_xyz = 1,10! 7 !3->4
               diagonal_spin = diagonal
               diag_offdiag_pack: do diag_offdiag = 1,2 !Test
                  i_ip1 = 1
                  irrep1_sog: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                     n_if1 = ua1%symadapt_spor_partner&
                          (i_ir1,quadrupel%l1)%N_independent_fcts
                     partner1_sog: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                        i_ip2 = i_ip1
                        i_ir2 = i_ir1
                        i_pa2 = i_pa1
                        needed_sog: if ( associated( &
                             symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int ) ) then
                           saint_2cob_dipoleg_real => &
                                symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int
                           saint_2cob_dipoleg_imag => &
                                symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int_imag
                           buf_size = 2*size(saint_2cob_dipoleg_real)
                           allocate(buffer(buf_size),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: allocating at packing failed")
                           i_buf_lower = 1
                           i_buf_upper = 0
                           diagonal_block_sog: if ( diag_offdiag == 1) then
                              ! diagonal block with triangular storage
                              qudrupel_required_sog: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_1_sog: do i_if1 = 1, n_if1
                                    if ( diagonal_spin ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1_sog: do i_c1 = 1, n_dim1
                                       ind_fct_2_sog: do i_if2 = 1, n_if2
                                          if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_dim2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2

                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_dipoleg_real(1:nn_c2,i_c1,i_if2,i_if1)

                                          i_buf_lower = i_buf_lower + nn_c2
                                          i_buf_upper = i_buf_upper + nn_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_dipoleg_imag(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_buf_lower = i_buf_lower + nn_c2

                                       enddo ind_fct_2_sog
                                    enddo contraction1_sog
                                 enddo ind_fct_1_sog
                              endif qudrupel_required_sog

                           else

                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                              ind_fct_1_od_sog: do i_if1 = 1, n_if1
                                 contraction1_od_sog: do i_c1 = 1, n_dim1
                                    i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)

                                    ind_fct_2_od_sog: do i_if2 = 1, n_if2
                                       i_buf_upper = i_buf_upper + n_dim2
                                       buffer(i_buf_lower:i_buf_upper) = &
                                            saint_2cob_dipoleg_real(1:n_dim2,i_c1,i_if2,i_if1)

                                       i_buf_lower = i_buf_lower + n_dim2
                                       i_buf_upper = i_buf_upper + n_dim2
                                       buffer(i_buf_lower:i_buf_upper) = &
                                            saint_2cob_dipoleg_imag(1:n_dim2,i_c1,i_if2,i_if1)

                                       i_buf_lower = i_buf_lower + n_dim2
                                       i_meta2 = i_meta2 + n_dim2
                                    enddo ind_fct_2_od_sog
                                    i_meta1 = i_meta1 + 1
                                 enddo contraction1_od_sog
                              enddo ind_fct_1_od_sog
                           endif diagonal_block_sog
                           n_pack=n_pack+1
                           if (i_buf_upper > buf_size) then
                              DPRINT 'pack_integralfile_dipole :wrong buf_size !!!'
                           end if
                           call commpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer size failed")
                           call commpack(buffer,i_buf_upper,1,info)

                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer failed")

                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: deallocating at packing failed")

                        endif needed_sog
                        i_ip1 = i_ip1 + 1
                     enddo partner1_sog
                  enddo irrep1_sog
               end do diag_offdiag_pack
            enddo xyz_sog
         endif gten!

      hfcop1: if(operations_hfc) then
            n_pack=0
            if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
               n_dim1 = n_exp1
               n_dim2 = n_exp2
            else
               n_dim1 = n_c1
               n_dim2 = n_c2
            end if
            ua_loop:Do i_ua = 1, n_unique_atoms
            xyz_hfc: do i_xyz = 1, 9
               diagonal_spin = diagonal
               diag_offdiag_pack1: do diag_offdiag = 1,2
                  i_ip1 = 1
                  irrep1_hfc: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                     n_if1 = ua1%symadapt_spor_partner&
                          (i_ir1,quadrupel%l1)%N_independent_fcts
                     partner1_hfc: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                        i_ip2 = i_ip1
                        i_ir2 = i_ir1
                        i_pa2 = i_pa1
                        needed_hfc: if ( associated( &
                             symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int ) ) then
                           saint_2cob_hfc_real => &
                                symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int
                           saint_2cob_hfc_imag => &
                                symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int_imag
                           buf_size = 2*size(saint_2cob_hfc_real)
                           allocate(buffer(buf_size),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: allocating at packing failed")
                           i_buf_lower = 1
                           i_buf_upper = 0
                           diagonal_block_hfc: if ( diag_offdiag == 1) then
                              ! diagonal block with triangular storage
                              qudrupel_required_hfc: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_1_hfc: do i_if1 = 1, n_if1
                                    if ( diagonal_spin ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1_hfc: do i_c1 = 1, n_dim1
                                       ind_fct_2_hfc: do i_if2 = 1, n_if2
                                          if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_dim2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2

                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_hfc_real(1:nn_c2,i_c1,i_if2,i_if1)

                                          i_buf_lower = i_buf_lower + nn_c2
                                          i_buf_upper = i_buf_upper + nn_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_hfc_imag(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_buf_lower = i_buf_lower + nn_c2

                                       enddo ind_fct_2_hfc
                                    enddo contraction1_hfc
                                 enddo ind_fct_1_hfc
                              endif qudrupel_required_hfc

                           else

                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                              ind_fct_1_od_hfc: do i_if1 = 1, n_if1
                                 contraction1_od_hfc: do i_c1 = 1, n_dim1
                                    i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)

                                    ind_fct_2_od_hfc: do i_if2 = 1, n_if2
                                       i_buf_upper = i_buf_upper + n_dim2
                                       buffer(i_buf_lower:i_buf_upper) = &
                                            saint_2cob_hfc_real(1:n_dim2,i_c1,i_if2,i_if1)

                                       i_buf_lower = i_buf_lower + n_dim2
                                       i_buf_upper = i_buf_upper + n_dim2
                                       buffer(i_buf_lower:i_buf_upper) = &
                                            saint_2cob_hfc_imag(1:n_dim2,i_c1,i_if2,i_if1)

                                       i_buf_lower = i_buf_lower + n_dim2
                                       i_meta2 = i_meta2 + n_dim2
                                    enddo ind_fct_2_od_hfc
                                    i_meta1 = i_meta1 + 1
                                 enddo contraction1_od_hfc
                              enddo ind_fct_1_od_hfc
                           endif diagonal_block_hfc
                           n_pack=n_pack+1
                           if (i_buf_upper > buf_size) then
                              DPRINT 'pack_integralfile_dipole :wrong buf_size !!!'
                           end if
                           call commpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer size failed")
                           call commpack(buffer,i_buf_upper,1,info)

                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer failed")

                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: deallocating at packing failed")

                        endif needed_hfc
                        i_ip1 = i_ip1 + 1
                     enddo partner1_hfc
                  enddo irrep1_hfc
               end do diag_offdiag_pack1
            enddo xyz_hfc
         end Do ua_loop
         endif hfcop1!
#endif

      else ! not option_spin_orbit
         if (operations_dipole) then
         xyz: do i_xyz = 1, 3
            i_ip1 = 1
            irrep1: do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                  i_ip2 = 1
                  irrep2: do i_ir2 = 1, i_ir1
                     partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                        needed: if ( associated( &
                             symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int ) ) then
                           saint_2cob_dipole => symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int
                           allocate(buffer(size(saint_2cob_dipole)),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: allocating at packing failed")
                           i_buf_lower = 1
                           i_buf_upper = 0
                           diagonal_block: if (i_ip1 .eq. i_ip2) then
                              ! diagonal block with triangular storage
                              qudrupel_required: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_1: do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1: do i_c1 = 1, n_c1
                                       ind_fct_2: do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_c2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_dipole(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_buf_lower = i_buf_lower + nn_c2
                                       enddo ind_fct_2
                                    enddo contraction1
                                 enddo ind_fct_1
                              endif qudrupel_required
                           else
                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              ind_fct_1_od: do i_if1 = 1, n_if1
                                 contraction1_od: do i_c1 = 1, n_c1
                                    ind_fct_2_od: do i_if2 = 1, n_if2
                                       i_buf_upper = i_buf_upper + n_c2
                                       buffer(i_buf_lower:i_buf_upper) = &
                                            saint_2cob_dipole(1:n_c2,i_c1,i_if2,i_if1)
                                       i_buf_lower = i_buf_lower + n_c2
                                    enddo ind_fct_2_od
                                 enddo contraction1_od
                              enddo ind_fct_1_od
                           endif diagonal_block
                           call commpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer size failed")
                           call commpack(buffer,i_buf_upper,1,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer failed")
                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: deallocating at packing failed")
                        endif needed
                        i_ip2 = i_ip2 + 1
                     enddo partner2
                  enddo irrep2
                  i_ip1 = i_ip1 + 1
               enddo partner1
            enddo irrep1
         enddo xyz
      endif

#ifdef WITH_GTENSOR
      hfc_1:if (operations_hfc) then
         if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
            n_dim1 = n_exp1
            n_dim2 = n_exp2
         else
            n_dim1 = n_c1
            n_dim2 = n_c2
         end if

        ua_hfc1 : do i_ua = 1, n_unique_atoms
            xyz_hfc1: do i_xyz = 1, 7
                i_ip1 = 1
                irrep_hfc1: do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner_hfc1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                  i_ip2 = i_ip1
                  i_ir2 = i_ir1
                  i_pa2 = i_pa1
                        needed_hfc1: if ( associated( &
                             symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int ) ) then
                           saint_2cob_hfc => symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int
                           buf_size = size(saint_2cob_hfc)

                           allocate(buffer(buf_size),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: allocating at packing failed")
                           i_buf_lower = 1
                           i_buf_upper = 0

                              ! diagonal block with triangular storage
                              qudrupel_required_hfc1: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_hfc1: do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction_hfc1: do i_c1 = 1, n_dim1
                                       ind_fct_hfc2: do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_dim2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2
                                          buffer(i_buf_lower:i_buf_upper) = &
                                               saint_2cob_hfc(1:nn_c2,i_c1,i_if2,i_if1)
                                          i_buf_lower = i_buf_lower + nn_c2

                                       enddo ind_fct_hfc2

                                    enddo contraction_hfc1

                                 enddo ind_fct_hfc1

                              endif qudrupel_required_hfc1


                           call commpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer size failed")
                           call commpack(buffer,i_buf_upper,1,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_send: packing buffer failed")
                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_send: deallocating at packing failed")
                        end if needed_hfc1
                          i_ip1 = i_ip1 + 1
                     enddo partner_hfc1
                  enddo irrep_hfc1
               enddo xyz_hfc1
            end do ua_hfc1
         end if hfc_1
#endif

         endif
      DPRINT 'done pack_integralfile_dipole'
    end subroutine pack_integralfile_dipole

  end subroutine int_send_dipole_send
  !*************************************************************
  subroutine int_send_dipole_receive
    !  Purpose:
    !    + receiving of contracted and symmetry adapted
    !      integrals for one quadrupel on master
    !  called by: integral_main_dipole, int_send_dipole_shutdown
    !** End of interface *****************************************
    use operations_module, only: operations_gtensor, operations_dipole, operations_hfc
    implicit none
    !------------ Declaration of local variables ---------------------
    logical :: diagonal, restart_timer
    integer(kind=i4_kind) :: n_c1, n_c2,n_exp1,n_exp2, info
    type(quadrupel_type) :: quadrupel
    type(unique_atom_type), pointer :: ua1, ua2
    !------------ Executable code ------------------------------------
    DPRINT 'int_send_dipole_receive: ENTERED'
    if ( timer_int_idle_2cob3c(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_receive_2cob3c(integralpar_i_int_part))

    call quadrupel_unpack(quadrupel)


       if(options_kinematic_factors) then
          call communpack(n_exp1,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: unpacking n_exp1 failed")
          call communpack(n_exp2,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: unpacking n_exp2 failed")
       end if

          call communpack(n_c1,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: unpacking n_c1 failed")
          call communpack(n_c2,info)
          if (info .ne. 0) call error_handler( &
            "int_send_dipole_send: unpacking n_c2 failed")



    if ( output_int_loops ) then
       write(output_unit,*) &
            "int_send_dipole_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
       write(stdout_unit,*) &
            "int_send_dipole_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
    endif

    diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
               (quadrupel%l1 == quadrupel%l2)
    ua1 => unique_atoms(quadrupel%ua1)
    ua2 => unique_atoms(quadrupel%ua2)

    if ( integralpar_2cob_dipole ) then
       if ( output_int_deeploops ) call write_to_output_units( &
            "int_send_dipole_receive: unpacking dipole")
       call unpack_integralfile_dipole()
    endif

    n_missing_quadrupels = n_missing_quadrupels - 1

    call stop_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
    if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_dipole_receive: done")
  contains

    subroutine unpack_integralfile_dipole()

      implicit none
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
           i_ip1, i_ip2, n_if1, n_if2, i_if1, i_if2, &
           nn_c2, i_c1, i_meta, i_meta1, i_meta2, status, &
           info, i_buf_lower,i_buf_upper, n_dim1, n_dim2, diag_offdiag, i_ua
      real(kind=r8_kind), allocatable :: buffer(:)
      type(dipole_integral_type), pointer :: dipole
      type(dipole_integral_type), pointer :: dipoleg
      type(dipole_integral_type), pointer :: hfc
      logical                             :: diagonal_spin
      !------------ Executable code ----------------------------

      if (options_spin_orbit) then
         !
         ! SPIN ORBIT
         !
        nopgt:  if(operations_dipole) then
         unpack_xyz_so: do i_xyz = 1, 3
                  diagonal_spin = diagonal
                  i_ip1 = 1
                  irrep1_so: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                     n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                     partner1_so: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                        i_ip2 = 1
                        irrep2_so: do i_ir2 = 1, i_ir1
                           partner2_so: do i_pa2 = 1, symmetry_data_n_partners_proj(i_ir2)

                              dipole => dipole_integrals(i_xyz,i_ip2,i_ip1)
                              needed_so: if ( dipole%use .ne. dipole_unused) then

                                 call communpack(i_buf_upper,info)
                                 if (info .ne. 0) call error_handler( &
                                      "int_send_dipole_receive: unpacking buffersize failed")

                                 allocate(buffer(i_buf_upper),stat=status)
                                 if (status .ne. 0) call error_handler( &
                                      "int_send_dipole_receive so(1): allocating at unpacking failed")
                                 call communpack(buffer,i_buf_upper,1,info)
                                 if (info .ne. 0) call error_handler( &
                                      "int_send_dipole_receive: unpacking buffer failed")

                                 i_buf_lower = 1
                                 i_buf_upper = 0
                                 diagonal_block_so: if (i_ip1 .eq. i_ip2) then
                                    ! diagonal block with triangular storage
                                    qudrupel_required_so: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                       ind_fct_1_so: do i_if1 = 1, n_if1
                                          if ( diagonal_spin ) then
                                             n_if2 = i_if1
                                          else
                                             n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                                  N_independent_fcts
                                          endif
                                          contraction1_so: do i_c1 = 1, n_c1
                                             i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1)
                                             ind_fct_2_so: do i_if2 = 1, n_if2
                                                if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                   nn_c2 = i_c1
                                                else
                                                   nn_c2 = n_c2
                                                endif

                                                if ((i_meta+nn_c2-1).gt.size(dipole%diagonal)) then
                                                   DPRINT 'stopped because of overflow in dipole%diagonal'
                                                endif

                                                if (i_buf_upper.gt.size(buffer)) then
                                                   DPRINT 'stopped because of overflow in buffer'
                                                endif

                                                i_buf_upper = i_buf_upper + nn_c2
                                                dipole%diagonal(i_meta:i_meta+nn_c2-1) = &
                                                     buffer(i_buf_lower:i_buf_upper)

                                                i_buf_lower = i_buf_lower + nn_c2
                                                i_buf_upper = i_buf_upper + nn_c2
                                                if (i_buf_upper.gt.size(buffer)) then
                                                   DPRINT 'stopped because of overflow in buffer imag'
                                                endif
                                                dipole%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                                     buffer(i_buf_lower:i_buf_upper)

                                                i_buf_lower = i_buf_lower + nn_c2
                                                i_meta = i_meta + nn_c2
                                             enddo ind_fct_2_so
                                          enddo contraction1_so
                                       enddo ind_fct_1_so
                                    endif qudrupel_required_so
                                 else
                                    ! offdiagonal block with full storage
                                    n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                    i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                    ind_fct_1_od_so: do i_if1 = 1, n_if1
                                       contraction1_od_so: do i_c1 = 1, n_c1
                                          i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                          ind_fct_2_od_so: do i_if2 = 1, n_if2
                                             i_buf_upper = i_buf_upper + n_c2
                                             dipole%offdiagonal(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                                  buffer(i_buf_lower:i_buf_upper)

                                             i_buf_lower = i_buf_lower + n_c2
                                             i_buf_upper = i_buf_upper + n_c2
                                             dipole%offdiagonal_imag(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                                  buffer(i_buf_lower:i_buf_upper)

                                             i_buf_lower = i_buf_lower + n_c2
                                             i_meta2 = i_meta2 + n_c2
                                          enddo ind_fct_2_od_so
                                          i_meta1 = i_meta1 + 1
                                       enddo contraction1_od_so
                                    enddo ind_fct_1_od_so
                                 endif diagonal_block_so
                                 deallocate(buffer,stat=status)
                                 if (status .ne. 0) call error_handler( &
                                      "int_send_dipole_receive so(2): deallocating at unpacking failed")
                             !    endif test
                              endif needed_so
                             i_ip2 = i_ip2 + 1
                          enddo partner2_so
                      enddo irrep2_so
                      i_ip1 = i_ip1 + 1
                   enddo partner1_so
                enddo irrep1_so
       enddo unpack_xyz_so
        endif nopgt

#ifdef WITH_GTENSOR
        gten: if(operations_gtensor) then !this is unpack
           if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
              n_dim1 = n_exp1
              n_dim2 = n_exp2
              else
              n_dim1 = n_c1
              n_dim2 = n_c2
           end if

           xyzs_sog: do i_xyz = 1, 10 !DG 7+3 !3->4
                    diagonal_spin = diagonal
                    diag_offdiag_unpack:  do diag_offdiag = 1,2 !Test
                    i_ip1 = 1
                    irrep1_sog: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                       n_if1 = ua1%symadapt_spor_partner&
                            (i_ir1,quadrupel%l1)%N_independent_fcts
                       partner1_sog: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                          i_ip2 = i_ip1
                          i_ir2 = i_ir1
                          i_pa2 = i_pa1
                                   !dipoleg => dipoleg_integrals(i_xyz,diag_offdiag,i_ip1)
                          if(diag_offdiag==1) dipoleg => gten_integrals_diag%integrals(i_xyz,i_ip1,1)
                          if(diag_offdiag==2) dipoleg => gten_integrals_offdiag%integrals(i_xyz,i_ip1,1)
                                   needed_sog: if (dipoleg%use .ne. dipole_unused) then

                                   call communpack(i_buf_upper,info)

                                   if (info .ne. 0) call error_handler( &
                                        "int_send_dipole_receive: unpacking buffersize failed")

                                   if (i_buf_upper < 0) then
                                     DPRINT 'unpack_integralfile_dipole :wrong buffer size'
                                   end if
                                   allocate(buffer(i_buf_upper),stat=status)
                                   if (status .ne. 0) call error_handler( &
                                        "int_send_dipole_receivesog(3): allocating at unpacking failed")

                                   call communpack(buffer,i_buf_upper,1,info)

                                   if (info .ne. 0) call error_handler( &
                                        "int_send_dipole_receive: unpacking buffer failed")

                                   i_buf_lower = 1
                                   i_buf_upper = 0

                                   diagonal_block_sog: if ( diag_offdiag==1 ) then
                                      ! diagonal block with triangular storage
                                      qudrupel_required_sog:if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                      ind_fct_1_sog: do i_if1 = 1, n_if1
                                         if ( diagonal_spin ) then
                                            n_if2 = i_if1
                                         else
                                            n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                                 N_independent_fcts
                                         endif
                                         contraction1_sog: do i_c1 = 1, n_dim1
                                                                                            !metaindexes (diff contract
                                               i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1) !and uncontr size)

                                            ind_fct_2_sog: do i_if2 = 1, n_if2
                                               if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                  nn_c2 = i_c1
                                               else
                                                  nn_c2 = n_dim2
                                               endif

                                               i_buf_upper = i_buf_upper + nn_c2

                                               dipoleg%diagonal(i_meta:i_meta+nn_c2-1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + nn_c2
                                               i_buf_upper = i_buf_upper + nn_c2

                                               dipoleg%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + nn_c2
                                               i_meta = i_meta + nn_c2
                                            enddo ind_fct_2_sog
                                         enddo contraction1_sog
                                      enddo ind_fct_1_sog
                                   end if qudrupel_required_sog

                                   else diagonal_block_sog
                                      ! offdiagonal block with full storage
                                      n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                      if(options_kinematic_factors) then
                                         i_meta1 = orbitalprojection_spor_ob_k(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      else
                                         i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      end if
                                      ind_fct_1_od_sog: do i_if1 = 1, n_if1
                                         contraction1_od_sog: do i_c1 = 1, n_dim1
                                            if(options_kinematic_factors) then
                                               i_meta2 = orbitalprojection_spor_ob_k(i_ir2,quadrupel%l2,quadrupel%ua2)
                                            else
                                               i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                            end if
                                            ind_fct_2_od_sog: do i_if2 = 1, n_if2
                                               i_buf_upper = i_buf_upper + n_dim2

                                               dipoleg%offdiagonal(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                                    buffer(i_buf_lower:i_buf_upper)


                                               i_buf_lower = i_buf_lower + n_dim2
                                               i_buf_upper = i_buf_upper + n_dim2
                                               dipoleg%offdiagonal_imag(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + n_dim2
                                               i_meta2 = i_meta2 + n_dim2
                                            enddo ind_fct_2_od_sog
                                            i_meta1 = i_meta1 + 1
                                         enddo contraction1_od_sog
                                      enddo ind_fct_1_od_sog
                                   endif diagonal_block_sog
                                   deallocate(buffer,stat=status)
                                   if (status .ne. 0) call error_handler( &
                                        "int_send_dipole_receive sog off(4): deallocating at unpacking failed")
                                endif needed_sog
                          i_ip1 = i_ip1 + 1
                       enddo partner1_sog
                    enddo irrep1_sog
                 end do diag_offdiag_unpack
           enddo xyzs_sog
     endif gten ! operations_gten

    hfcop: if(operations_hfc) then !this is unpack
           if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
              n_dim1 = n_exp1
              n_dim2 = n_exp2
              else
              n_dim1 = n_c1
              n_dim2 = n_c2
           end if
        ualoop:  do i_ua = 1, n_unique_atoms
           xyzs_hfc: do i_xyz = 1, 9
                    diagonal_spin = diagonal
                    diag_offdiag_unpack_hfc:  do diag_offdiag = 1,2
                    i_ip1 = 1
                    irrep1_hfc: do i_ir1 = 1, symmetry_data_n_proj_irreps()
                       n_if1 = ua1%symadapt_spor_partner&
                            (i_ir1,quadrupel%l1)%N_independent_fcts
                       partner1_hfc: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                          i_ip2 = i_ip1
                          i_ir2 = i_ir1
                          i_pa2 = i_pa1

                                  if(diag_offdiag==1) hfc => hfc_integrals_diag%integrals(i_xyz,i_ip1,i_ua)
                                  if(diag_offdiag==2) hfc => hfc_integrals_offdiag%integrals(i_xyz,i_ip1,i_ua)
                                   needed_hfc: if (hfc%use .ne. dipole_unused) then

                                   call communpack(i_buf_upper,info)

                                   if (info .ne. 0) call error_handler( &
                                        "int_send_dipole_receive: unpacking buffersize failed")

                                   if (i_buf_upper < 0) then
                                     DPRINT 'unpack_integralfile_dipole :wrong buffer size'
                                   end if
                                   allocate(buffer(i_buf_upper),stat=status)
                                   if (status .ne. 0) call error_handler( &
                                        "int_send_dipole_receivesog(3): allocating at unpacking failed")

                                   call communpack(buffer,i_buf_upper,1,info)

                                   if (info .ne. 0) call error_handler( &
                                        "int_send_dipole_receive: unpacking buffer failed")

                                   i_buf_lower = 1
                                   i_buf_upper = 0

                                   diagonal_block_hfc: if ( diag_offdiag == 1 ) then
                                      ! diagonal block with triangular storage
                                      qudrupel_required_hfc:if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                      ind_fct_1_hfc: do i_if1 = 1, n_if1
                                         if ( diagonal_spin ) then
                                            n_if2 = i_if1
                                         else
                                            n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)% &
                                                 N_independent_fcts
                                         endif
                                         contraction1_hfc: do i_c1 = 1, n_dim1
                                                                                                !metaindexes (diff contract
                                               i_meta = metaindex_p(quadrupel,i_ir1,i_if1,i_c1) !and uncontr size)

                                            ind_fct_2_hfc: do i_if2 = 1, n_if2
                                               if ( diagonal_spin .and. i_if1 .eq. i_if2 ) then
                                                  nn_c2 = i_c1
                                               else
                                                  nn_c2 = n_dim2
                                               endif

                                               i_buf_upper = i_buf_upper + nn_c2

                                               hfc%diagonal(i_meta:i_meta+nn_c2-1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + nn_c2
                                               i_buf_upper = i_buf_upper + nn_c2

                                               hfc%diagonal_imag(i_meta:i_meta+nn_c2-1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + nn_c2
                                               i_meta = i_meta + nn_c2
                                            enddo ind_fct_2_hfc
                                         enddo contraction1_hfc
                                      enddo ind_fct_1_hfc
                                   end if qudrupel_required_hfc

                                else
                                      ! offdiagonal block with full storage
                                      n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                                      if(options_kinematic_factors) then
                                         i_meta1 = orbitalprojection_spor_ob_k(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      else
                                         i_meta1 = orbitalprojection_spor_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                                      end if
                                      ind_fct_1_od_hfc: do i_if1 = 1, n_if1
                                         contraction1_od_hfc: do i_c1 = 1, n_dim1
                                            if(options_kinematic_factors) then
                                               i_meta2 = orbitalprojection_spor_ob_k(i_ir2,quadrupel%l2,quadrupel%ua2)
                                            else
                                               i_meta2 = orbitalprojection_spor_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                            end if
                                            ind_fct_2_od_hfc: do i_if2 = 1, n_if2
                                               i_buf_upper = i_buf_upper + n_dim2

                                               hfc%offdiagonal(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                                    buffer(i_buf_lower:i_buf_upper)


                                               i_buf_lower = i_buf_lower + n_dim2
                                               i_buf_upper = i_buf_upper + n_dim2
                                               hfc%offdiagonal_imag(i_meta2:i_meta2+n_dim2-1,i_meta1) = &
                                                    buffer(i_buf_lower:i_buf_upper)

                                               i_buf_lower = i_buf_lower + n_dim2
                                               i_meta2 = i_meta2 + n_dim2
                                            enddo ind_fct_2_od_hfc
                                            i_meta1 = i_meta1 + 1
                                         enddo contraction1_od_hfc
                                      enddo ind_fct_1_od_hfc
                                   endif diagonal_block_hfc
                                   deallocate(buffer,stat=status)
                                   if (status .ne. 0) call error_handler( &
                                        "int_send_dipole_receive sog off(4): deallocating at unpacking failed")
                                endif needed_hfc
                          i_ip1 = i_ip1 + 1
                       enddo partner1_hfc
                    enddo irrep1_hfc
                 end do diag_offdiag_unpack_hfc
              enddo xyzs_hfc
           end do ualoop
        endif hfcop ! operations_hfc
#endif

      else ! not option_spin_orbit
         if (operations_dipole) then
         xyz: do i_xyz = 1, 3
            i_ip1 = 0
            irrep1: do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
               partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                  i_ip1 = i_ip1 + 1
                  i_ip2 = 0
                  irrep2: do i_ir2 = 1, i_ir1
                     partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                        i_ip2 = i_ip2 + 1
                        dipole => dipole_integrals(i_xyz, i_ip2, i_ip1)
                        needed: if (dipole % use .ne. dipole_unused) then
                           call communpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_receive: unpacking buffersize failed")
                           allocate(buffer(i_buf_upper),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_receive: allocating at unpacking failed")
                           call communpack(buffer,i_buf_upper,1,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_receive: unpacking buffer failed")
                           i_buf_lower = 1
                           i_buf_upper = 0
                           diagonal_block: if (i_ip1 .eq. i_ip2) then
                              ! diagonal block with triangular storage
                              qudrupel_required: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 ind_fct_1: do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction1: do i_c1 = 1, n_c1
                                       i_meta = metaindex(quadrupel,i_ir1,i_if1,i_c1)
                                       ind_fct_2: do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_c2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2
                                          dipole%diagonal(i_meta:i_meta+nn_c2-1) = &
                                               buffer(i_buf_lower:i_buf_upper)
                                          i_buf_lower = i_buf_lower + nn_c2
                                          i_meta = i_meta + nn_c2
                                       enddo ind_fct_2
                                    enddo contraction1
                                 enddo ind_fct_1
                              endif qudrupel_required
                           else
                              ! offdiagonal block with full storage
                              n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)%N_independent_fcts
                              i_meta1 = orbitalprojection_ob(i_ir1,quadrupel%l1,quadrupel%ua1)
                              ind_fct_1_od: do i_if1 = 1, n_if1
                                 contraction1_od: do i_c1 = 1, n_c1
                                    i_meta2 = orbitalprojection_ob(i_ir2,quadrupel%l2,quadrupel%ua2)
                                    ind_fct_2_od: do i_if2 = 1, n_if2
                                       i_buf_upper = i_buf_upper + n_c2
                                       dipole%offdiagonal(i_meta2:i_meta2+n_c2-1,i_meta1) = &
                                            buffer(i_buf_lower:i_buf_upper)
                                       i_buf_lower = i_buf_lower + n_c2
                                       i_meta2 = i_meta2 + n_c2
                                    enddo ind_fct_2_od
                                    i_meta1 = i_meta1 + 1
                                 enddo contraction1_od
                              enddo ind_fct_1_od
                           endif diagonal_block
                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_receive: deallocating at unpacking failed")
                       endif needed
                    enddo partner2
                 enddo irrep2
              enddo partner1
           enddo irrep1
        enddo xyz
     end if

#ifdef WITH_GTENSOR
     if (operations_hfc) then
        if(options_kinematic_factors) then ! I should use uncontracted dimensions in this case
            n_dim1 = n_exp1
            n_dim2 = n_exp2
         else
            n_dim1 = n_c1
            n_dim2 = n_c2
         end if
         do i_ua = 1, n_unique_atoms

          do i_xyz = 1, 7
            i_ip1 = 1
             do i_ir1 = 1, symmetry_data_n_irreps()
               n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                  i_ip2 = i_ip1
                  i_ir2 = i_ir1
                  i_pa2 = i_pa1
                        hfc => hfc_integrals_diag%integrals(i_xyz,i_ip1,i_ua)

                           call communpack(i_buf_upper,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_receive: unpacking buffersize failed")
                           allocate(buffer(i_buf_upper),stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_receive: allocating at unpacking failed")
                           call communpack(buffer,i_buf_upper,1,info)
                           if (info .ne. 0) call error_handler( &
                                "int_send_dipole_receive: unpacking buffer failed")
                           i_buf_lower = 1
                           i_buf_upper = 0

                              ! diagonal block with triangular storage
                              qudrupel_required_hfc1: if ( qudrupel_required_for_diagonal(quadrupel) ) then
                                 do i_if1 = 1, n_if1
                                    if ( diagonal ) then
                                       n_if2 = i_if1
                                    else
                                       n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)% &
                                            N_independent_fcts
                                    endif
                                    contraction_hfc1: do i_c1 = 1, n_dim1
                                       i_meta = metaindex(quadrupel,i_ir1,i_if1,i_c1)
                                       ind_fct_hfc2: do i_if2 = 1, n_if2
                                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                             nn_c2 = i_c1
                                          else
                                             nn_c2 = n_dim2
                                          endif
                                          i_buf_upper = i_buf_upper + nn_c2
                                          hfc%diagonal(i_meta:i_meta+nn_c2-1) = &
                                               buffer(i_buf_lower:i_buf_upper)
                                          i_buf_lower = i_buf_lower + nn_c2
                                          i_meta = i_meta + nn_c2
                                       enddo ind_fct_hfc2
                                    enddo contraction_hfc1
                                 enddo
                              endif qudrupel_required_hfc1

                           deallocate(buffer,stat=status)
                           if (status .ne. 0) call error_handler( &
                                "int_send_dipole_receive: deallocating at unpacking failed")

                 i_ip1 = i_ip1 + 1
              enddo
           enddo
        enddo
     end do
     end if
#endif

     endif
     DPRINT ' unpack_integralfile_dipole ... done'
    end subroutine unpack_integralfile_dipole

  end subroutine int_send_dipole_receive

integer(kind=i4_kind) function metaindex(quadrupel,i_ir,i_ind1,i_exp1)
  !------------ Declaration of formal parameters -------------
  type(quadrupel_type),  intent(in)  :: quadrupel
  integer(kind=i4_kind), intent(in)  :: i_ir,i_ind1,i_exp1
  !** End of interface ***************************************
  metaindex = &
       metaindex_of_first_integral(i_ir)%ua1(quadrupel%ua1)% &
       l1(quadrupel%l1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
       l2(quadrupel%l2)
end function metaindex

integer(kind=i4_kind) function metaindex_p(quadrupel,i_ir,i_ind1,i_exp1)
!------------ Declaration of formal parameters -------------
type(quadrupel_type),  intent(in)  :: quadrupel
integer(kind=i4_kind), intent(in)  :: i_ir,i_ind1,i_exp1
!** End of interface ***************************************

metaindex_p = &
     metaindex_of_first_integral_p(i_ir)%ua1(quadrupel%ua1)% &
     l1(quadrupel%l1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
     l2(quadrupel%l2)

end function metaindex_p



logical function qudrupel_required_for_diagonal(quadrupel)
  !------------ Declaration of formal parameters -------------
type(quadrupel_type),  intent(in)  :: quadrupel
!** End of interface ***************************************
if ( quadrupel%ua1 .gt. quadrupel%ua2 ) then
 qudrupel_required_for_diagonal = .true.
 elseif( quadrupel%ua1 .eq. quadrupel%ua2 ) then
 qudrupel_required_for_diagonal = &
      quadrupel%l1 .ge. quadrupel%l2
else
 qudrupel_required_for_diagonal = .false.
endif
end function qudrupel_required_for_diagonal

!--------------- End of module ----------------------------------
end module int_send_dipole_module

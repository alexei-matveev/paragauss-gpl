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
module  int_send_2cob3c_spor_module
!---------------------------------------------------------------
!
!  Purpose: Sending, receiving and storing of
!           2 center orbital and three center integrals
!           for spin orbit relativistics.
!           The results of one quadrupel stored in
!           int_data_2cob3c_module are distributed.
!           They are either stored in intermediate
!           quadrupel files or directly stored in memory
!           in integralstore_module. In case files are used,
!           when executing int_send_2cob3c_spor_shutdown(), the
!           results are rewritten to tapes in the form used
!           in the SCF part. The mapping to the scf-
!           metaindex is also done.
!           Two center integrals are send to the master and
!           three center integrals are distributed over all hosts
!           splitting them over the fitfunction metaindex.
!
!
!
!  Module called by: main_slave, integral_main_2cob3c,
!    integral_calc_quad_2cob3c, integral_interrupt_2cob3c,
!    integral_setup_2cob3c, integral_shutdown_2cob3c
!
!  References: Publisher Document: Concepts of Integral Part
!
!  This is an intermediate version.
!  Changes have only be done to part for options_integral_on_file.
!  Direct memory mode is claerly nonsense !!!
!
!
!  Author: TB
!  Date: 7/95
!
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: TB
! Date:   1/97
! Description: Moved intermediate storage of integrals from
!              large direct access files to seperate sequential
!              files for each IRREP and quadrupel that are
!              handled by readwriteblocked-module.
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
! Modification (Please copy before editing)
! Author: SM
! Date:   20/01/2004
! Description: ...
!
!----------------------------------------------------------------
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
use integralpar_module,      only: integralpar_2cob_kin                        &
                                 , integralpar_2cob_nuc                        &
                                 , integralpar_2cob_ol                         &
                                 , integralpar_2cob_pvec                       &
                                 , integralpar_2cob_pvsp                       &
                                 , integralpar_2cob_pvxp                       &
                                 , integralpar_3c_co                           &
                                 , integralpar_3c_xc                           &
                                 , integralpar_3c_r2_pvsp                      &
                                 , integralpar_3c_r2_pvxp                      &
                                 , integralpar_3c_rcoul_pvsp                   &
                                 , integralpar_3c_rcoul_pvxp                   &
                                 , integralpar_i_int_part                      &
                                 , integralpar_rel_gradients                   &
                                 , integralpar_relativistic
use output_module,           only: output_int_loops                            &
                                 , output_int_deeploops                        &
                                 , output_slaveoperations                      &
                                 , output_int_progress
use iounitadmin_module,      only: output_unit                                 &
                                 , stdout_unit                                 &
                                 , openget_iounit                              &
                                 , get_nbr_free_iounits                        &
                                 , write_to_output_units                       &
                                 , write_to_trace_unit                         &
                                 , returnclose_iounit
use filename_module,         only: tmpfile                                     &
                                 , filename_tmpdir                             &
                                 , filesystem_is_parallel
use comm_module,             only: comm_myindex                                &
                                 , comm_i_am_master                            &
                                 , comm_get_n_processors                       &
                                 , comm_all_other_hosts                        &
                                 , comm_init_send                              &
                                 , comm_send                                   &
                                 , commpack                                    &
                                 , communpack                                  &
                                 , comm_save_recv
use msgtag_module,           only: msgtag_int_2cob3c_result                    &
                                 , msgtag_int_2cob3c_filesize
use symmetry_data_module,    only: symmetry_data_n_irreps                      &
                                 , symmetry_data_n_proj_irreps
use unique_atom_module,      only: unique_atom_type                            &
                                 , unique_atom_basis_type                      &
                                 , unique_atom_partner_type                    &
                                 , unique_atoms                                &
                                 , N_unique_atoms
use quadrupel_module,        only: quadrupel_type                              &
                                 , quadrupel_pack                              &
                                 , quadrupel_unpack
use timer_module,            only: timer_int_pack_2cob3c                       &
                                 , timer_int_commsend_2cob3c                   &
                                 , timer_int_receive_2cob3c                    &
                                 , timer_int_write_2cob3c                      &
                                 , timer_int_rewrite_2cob3c                    &
                                 , timer_int_idle_2cob3c
use time_module,             only: start_timer                                 &
                                 , stop_timer
use readwriteblocked_module, only: readwriteblocked_tapehandle                 &
                                 , readwriteblocked_blocklength                &
                                 , readwriteblocked_startread                  &
                                 , readwriteblocked_read                       &
                                 , readwriteblocked_stopread                   &
                                 , readwriteblocked_startwrite                 &
                                 , readwriteblocked_write                      &
                                 , readwriteblocked_stopwrite
use options_module,          only: options_integrals_on_file
!use integralstore_module, only: integralstore_allocate

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
public int_send_2cob3c_spor_setup, int_send_2cob3c_spor_shutdown, &
       int_send_2cob3c_spor_send, int_send_2cob3c_spor_receive, &
       int_send_2cob3c_s_rec_filesizes


  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Declaration of private types ----

type, private ::  tapehandle_arr
   ! this type is for rel. gradients   ! for every irrep and every degree of freedom we
   ! need one tapehandle
   ! Note: the number of degrees of freedom is not the
   !       same for every irrep
   type(readwriteblocked_tapehandle),pointer :: th_p(:,:)
end type tapehandle_arr

type, private :: index_ind1_exp1_ua2_type
!!$   integer(kind=i4_kind), pointer :: l2(:,:) ! (i_l2)
   integer(kind=i4_kind), pointer :: l2(:)
end type index_ind1_exp1_ua2_type

type, private :: index_l1_type
   type(index_ind1_exp1_ua2_type), pointer :: ind1_exp1_ua2(:,:,:)
   ! (i_ind1,i_exp1,i_ua2)
end type index_l1_type

type, private :: index_ua1_type
!!$   type(index_l1_type), pointer :: l1(:,:) ! (i_l1)
   type(index_l1_type), pointer :: l1(:)
end type index_ua1_type

type, private :: index_ir_type
   type(index_ua1_type), pointer :: ua1(:) ! (i_ua1)
end type index_ir_type


!------------ Declaration of constants and variables ----
integer(kind=i4_kind) :: n_missing_quadrupels, &
     first_host, last_host, my_hostindex, n_file, n_tot_records, &
     n_tot_records_rel, &
     n_hosts_to_report_back, blocklength, n_missing_filesizes
integer(kind=i4_kind), dimension(:,:), allocatable :: &
     & borders_ch,  &
     & borders_s, &
     & borders_r2, &
     & borders_xc ! borders(3,comm_get_n_processors())
     ! 1. dim: 1: lower border, 2: upper border, 3: number of ff
integer(kind=i4_kind), allocatable :: ua_fileindex(:) ! (i_ua)
     ! to look up file index for l = 0
integer(kind=i4_kind), allocatable :: quadrupelfile_length(:,:,:), &
     quadrupelfile_length_allhosts(:,:,:,:)
type(index_ir_type), allocatable :: &
     metaindex_of_first_integral(:),metaindex_of_first_integral_rel(:) ! (i_ir)
     ! to look up metaindex of first integral and number
     ! n_records of following records for a given vector
     ! (i_ir,i_ua1,i_l1,i_ind1,i_exp2,i_ua2,i_l2)
type(tapehandle_arr),allocatable :: th_arr(:)
     ! th_arr(n_irrep): only used for rel. gradients; see tapehandle_arr


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine int_send_2cob3c_spor_setup(n_quads)
    !  Purpose:
    !   + calculates fitfunction index boundaries for
    !     splitting 3 center integrals
    !   + calculates help variables number of quadrupels
    !     to be received and total number of record in file
    !   + calculates help array for mapping integrals to
    !     1 dimensional direct access file.
    !   + Opening of direct access files.
    ! called by: integral_setup_2cob3c (on all hosts)
    !------------ Modules used ---------------------------------
    use fit_coeff_module,     only: fit_coeff_n_ch                             &
                                  , fit_coeff_n_s                              &
                                  , fit_coeff_n_r2                             &
                                  , fit_coeff_n_xc
    use int_send_2cob3c_spor, only: IX_CHFIT                                   &
                                  , IX_XCFIT                                   &
                                  , IX_CHFIT_S                                 &
                                  , IX_CHFIT_R2                                &
                                  , fitborders                                 &
                                  , setsplitting                               &
                                  , syncborders
    implicit none
    integer(kind=i4_kind) :: n_quads
    !** End of interface *****************************************

    integer (i4_kind) :: n_hosts, status, n_ch, n_s, n_r2, n_xc, i_ua1, &
         i_file, allocsize, basesize, n_free_iounits
    real (r8_kind), allocatable :: dummy(:)
#ifdef SO_INTEGRALS_IN_MEM
    ! FIXME: rm redundant decls:
    logical :: diagonal, diagonal_unique
    integer (i4_kind) :: i_ir, n_hosts, status, n_hosts_larger, &
         i_border, i_host, n_ch,&
         & n_s, n_r2,&
         & n_xc, i_ua1, i_ua2, &
         i_file, i_l1, allocsize, basesize, n_free_iounits
    integer (i4_kind) :: i_meta, i_ind1, i_exp1, i_l2, n_ind2, i_ind2, n_exp2
    type (unique_atom_type), pointer :: ua1, ua2
    type (unique_atom_basis_type), pointer :: uab1, uab2
    type (unique_atom_partner_type), pointer :: uap1, uap2
    type (readwriteblocked_tapehandle), pointer :: th_pointer
    integer (i4_kind), pointer :: metaindex(:)
    type (index_ind1_exp1_ua2_type), pointer :: ind1_exp1_ua2(:,:,:)
#endif
    !------------ Executable code ------------------------------------

    ! calculate borders
    first_host = 1
    last_host = comm_get_n_processors()
    n_hosts = last_host - first_host + 1
    if ( output_unit > 0 ) then
       write(output_unit,*)'int_send_2cob3c_spor_setup:'
       write(output_unit,*) " first_host:",first_host
       write(output_unit,*) " last_host:",last_host
    endif
    if ( n_hosts .lt. 1 ) call error_handler( &
         "int_send_2cob3c_spor_setup: wrong n_hosts" )

    allocate(&
        borders_ch(3,first_host:last_host), &
        borders_s(3,first_host:last_host), &
        borders_r2(3,first_host:last_host), &
        borders_xc(3,first_host:last_host), &
        stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "int_send_2cob3c_spor_setup: allocate of borders failed" )

    if (integralpar_3c_co) then
       n_ch = fit_coeff_n_ch()

       call FitBorders(n_ch,borders_ch)

       call SetSplitting(first_host,last_host,borders_ch(3,:),IX_CHFIT)
    else
       borders_ch = 0
    endif

    ! new for S integrals:
    if (integralpar_3c_rcoul_pvsp) then
       ASSERT(integralpar_3c_co)
       n_s = fit_coeff_n_s()

       call SyncBorders(borders_ch,borders_s,IX_CHFIT_S)
       n_s = n_s - sum(borders_s(3,:))
       ASSERT(n_s==0)

       call SetSplitting(first_host,last_host,borders_s(3,:),IX_CHFIT_S)

       if (.not.integralpar_3c_rcoul_pvxp) then
          ! dimensions same as borders_s
          ABORT("fix here: 1")
       endif
    else
       borders_s = 0
    endif

    ! new for R2 integrals:
    if (integralpar_3c_r2_pvsp) then
       ASSERT(integralpar_3c_co)
       n_r2 = fit_coeff_n_r2()

       call SyncBorders(borders_ch,borders_r2,IX_CHFIT_R2)
       n_r2 = n_r2 - sum(borders_r2(3,:))
       ASSERT(n_r2==0)

       call SetSplitting(first_host,last_host,borders_r2(3,:),IX_CHFIT_R2)

       if (.not.integralpar_3c_r2_pvxp) then
          ABORT("fix here: 2")
       endif
    else
       borders_r2 = 0
    endif

    if (integralpar_3c_xc) then
       n_xc = fit_coeff_n_xc()

       call FitBorders(n_xc,borders_xc)

       call SetSplitting(first_host,last_host,borders_xc(3,:),IX_XCFIT)
    else
       borders_xc = 0
    endif


    ! fetch my hostindex from comm_module
    my_hostindex = comm_myindex()



    if ( options_integrals_on_file() ) then


       !? print*,">>> setup: allocate and calculate ua_fileindex"
       ! allocate and calculate ua_fileindex
       allocate( ua_fileindex(N_unique_atoms), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: allocate of ua_fileindex failed" )
       i_file = 1
       do i_ua1 = 1, N_unique_atoms
          ua_fileindex(i_ua1) = i_file
          i_file = i_file + unique_atoms(i_ua1)%lmax_ob + 1
       enddo
       n_file = i_file - 1


       ! allocate quadrupelfile_length
       allocate( quadrupelfile_length(n_file,n_file,symmetry_data_n_proj_irreps()), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: allocate of quadrupelfile_length failed" )
       quadrupelfile_length = 0
       if ( filesystem_is_parallel ) then
          allocate( quadrupelfile_length_allhosts(n_file, n_file, &
               symmetry_data_n_proj_irreps(), comm_get_n_processors() ), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_setup: allocate of quadrupelfile_length_allhosts failed" )
          quadrupelfile_length = 0
          quadrupelfile_length_allhosts = 0
          n_missing_filesizes = comm_get_n_processors()
       endif


       ! blocklength for intermediate files to store quadrupels
       ! Test if sufficient memory is available, if not reduce blocksize
       ! take into account that swap (on a typical system, roughly
       ! two times hardware memory) will normally be used for allocation
       ! as well. This is all rough estimates only. The memory test will not
       ! this way on a shared memory architecture.
       blocklength = readwriteblocked_blocklength()
       n_free_iounits = get_nbr_free_iounits() * 3 ! swap
       basesize = 10 * blocklength
       allocsize = basesize + n_free_iounits * blocklength
       allocate( dummy(allocsize), stat=status )
       do while (status .ne. 0)
          blocklength = blocklength / 2
          if ( blocklength .lt. 10 ) call error_handler( &
               "int_send_2cob3c_spor_setup: insufficient memory for rewriting" )
          allocsize = basesize + n_free_iounits * blocklength
          allocate( dummy(allocsize), stat=status )
       enddo
       deallocate( dummy, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: deallocate of dummy failed" )

       if ( output_int_loops .and. output_unit > 0 ) then
          write (output_unit,*)
          write (output_unit,*) "##### using blocklength ", blocklength, " for quadrupel files #####"
          write (output_unit,*)
          write (stdout_unit,*)
          write (stdout_unit,*) "##### using blocklength ", blocklength, " for quadrupel files #####"
          write (stdout_unit,*)
       endif

    else ! .not.  options_integrals_on_file()
       ABORT('FIXME: ints in mem')
       ! do not compile:
#ifdef SO_INTEGRALS_IN_MEM
       ! calculate index_of_ua_l
       allocate( metaindex_of_first_integral(symmetry_data_n_proj_irreps()), &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral failed" )
       i_meta = 1
       do i_ir = 1,symmetry_data_n_proj_irreps()
          allocate( metaindex_of_first_integral(i_ir)%ua1(N_unique_atoms), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ua1 failed" )
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)

             allocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l1 failed" )
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)
                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                     uab1%N_uncontracted_fcts + uab1%N_contracted_fcts, i_ua1), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2
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

                         allocate( metaindex(0:n_l2), &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l2 failed" )
                         ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2 => metaindex
                         i_l2_:do i_l2 = 0, n_l2
                            diagonal = (diagonal_unique.and.(i_l1.eq.i_l2))
                            uab2 => ua2%l_ob(i_l2)

                            uap2 => ua2%symadapt_spor_partner(i_ir,i_l2)

                            if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2) then
                               n_ind2 = i_ind1
                            else
                               n_ind2 = uap2%N_independent_fcts
                            endif
                            n_records = 0
                            do i_ind2 = 1, n_ind2

                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 .and. i_ind1 .eq. i_ind2) then
                                  n_exp2 = i_exp1
                               else
                                  n_exp2 = uab2%N_uncontracted_fcts + &
                                       uab2%N_contracted_fcts
                               endif
                               n_records = n_records + n_exp2*2
                               ! factor 2 accounts for real and imaginary parts
                            enddo

                            metaindex(i_l2) = i_meta
                            i_meta = i_meta + n_records

                         enddo i_l2_! i_l2
                      enddo! i_ua2
                   enddo! i_exp1
                enddo! i_ind1

             enddo! i_l1
          enddo! i_ua1
       enddo! i_irrep
       n_tot_records = i_meta - 1
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_setup: total number of records: ",n_tot_records)

       ! calculate index_of_ua_l
       allocate( metaindex_of_first_integral_rel(symmetry_data_n_proj_irreps()), &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral failed" )
       i_meta = 1
       do i_ir = 1,symmetry_data_n_proj_irreps()
          allocate( metaindex_of_first_integral_rel(i_ir)%ua1(N_unique_atoms), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ua1 failed" )
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)

             allocate( metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                  stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l1 failed" )
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)

                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                     uab1%N_exponents, i_ua1), &
                     stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2
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
                              "int_send_2cob3c_spor_setup: allocate of metaindex_of_first_integral l2 failed" )
                         ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2 => metaindex
                         do i_l2 = 0, n_l2
                            diagonal = (diagonal_unique.and.(i_l1.eq.i_l2))
                            uab2 => ua2%l_ob(i_l2)

                            uap2 => ua2%symadapt_spor_partner(i_ir,i_l2)

                            if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2) then
                               n_ind2 = i_ind1
                            else
                               n_ind2 = uap2%N_independent_fcts
                            endif
                            n_records = 0
                            do i_ind2 = 1, n_ind2

                               if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 .and. i_ind1 .eq. i_ind2 ) then
                                  n_exp2 = i_exp1
                               else
                                  n_exp2 = uab2%N_exponents
                               endif
                               n_records = n_records + n_exp2*2
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
       n_tot_records_rel = i_meta - 1
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_setup: total number of records: ",n_tot_records_rel)

       ! for the case of cases:
       call integralstore_deallocate()

       ! allocating storage for integrals
       if ( integralstore_allocate ( &
            dim_2c = n_tot_records, &
            dim_2c_rel = n_tot_records_rel, &
            dim_3c_co = n_tot_records * borders_ch(3,my_hostindex), &
            dim_3c_xc = n_tot_records * borders_xc(3,my_hostindex), &
            need_2cob_kin = integralpar_2cob_kin .and. comm_i_am_master(), &
            need_2cob_nuc = integralpar_2cob_nuc .and. comm_i_am_master(), &
            need_2cob_kin_rel = integralpar_2cob_kin .and. comm_i_am_master() &
            .and. integralpar_relativistic , &
            need_2cob_nuc_rel = integralpar_2cob_nuc .and. comm_i_am_master() &
            .and. integralpar_relativistic, &
            need_2cob_pvsp = integralpar_2cob_pvsp .and. comm_i_am_master().and. &
            integralpar_relativistic, &
            need_2cob_pseudo = pseudopot_present .and. comm_i_am_master().and. &
            integralpar_relativistic, &
            need_2cob_pvxp = integralpar_2cob_pvxp .and. comm_i_am_master() &
            & .and. integralpar_relativistic, &
            need_2cob_pvec =& !<<<mdf
            & integralpar_2cob_pvec&
            & .and. comm_i_am_master()&
            & .and. integralpar_relativistic, &
            need_2cob_ol = integralpar_2cob_ol, &
            need_2cob_ol_rel = integralpar_2cob_ol .and. &
            integralpar_relativistic , &
            need_3c_xc = integralpar_3c_xc , &
            need_3c_co = integralpar_3c_co &
            )  ) &
            call error_handler( &
            "int_send_2cob3c_spor_setup: allocation of integral storage failed. set options_integrals_on_file .true.")
       !? print*," now call memstat"
#endif
    endif !  options_integrals_on_file()

    n_missing_quadrupels = n_quads
    if ( output_int_progress ) then
       call write_to_output_units("")
       call write_to_output_units( &
            "total number of 2 center orbital and 3 center quadrupels: ",n_quads)
       call write_to_output_units("")
    endif


    ! during shutdown procedure, all slaves have to report back
    n_hosts_to_report_back = comm_get_n_processors() - 1
  !  print*,"int_send_2cob3c_spor_setup: EXIT"


  end subroutine int_send_2cob3c_spor_setup
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cob3c_spor_shutdown
    !  Purpose:
    !    + waiting for missing integrals
    !    + deallocation of help arrays
    !    + reading of quadrupel files and writing
    !      them as sequential files as used in the scf part
    !    + for operations_integraltest, comparing integral files
    !      with the files produced by old lcgto
    !    + removal of quadrupel files
    ! called by: integral_shutdown_2cob3c (on all hosts)
    !** End of interface *****************************************
    use error_module, only: error
    use operations_module, only: operations_integral
    use spin_orbit_module, only: is_on                                         &
                               , op_FitTrafo !<<< options only
    implicit none

    !------------ Declaration of local types ---------------------
    type quadrupel_file_handle_type
       type(readwriteblocked_tapehandle) :: th
       integer(kind=i4_kind) :: n_if_c, total_length
       integer(kind=i4_kind), pointer ::  n_floats_if_c(:)
       ! the next elements are only used in the relativistic case
       integer(kind=i4_kind) :: n_if_c_uc
       integer(kind=i4_kind), pointer ::  n_floats_if_c_uc(:)
    end type quadrupel_file_handle_type

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: status, i_ir, i_unique,i_l, &
         i_file1, i_file2, n_file2, n_integrals, n_free_iounits, &
         i_file2_new, n_file2_new, n_concat, n_tot_records, n_tot_records_rel,&
         info, i_integral, i_if_c, n_file2_conc, n_files_left, n_if_c, n_irrep, &
         i_host, &
         unit_kin_nuc_real, unit_ol_real, unit_kin_nuc_imag, unit_ol_imag !!$,i_spin1,i_spin2
    integer(kind=i4_kind), allocatable :: blocksizes(:), &
         ! blocksizes(n_integrals)
    dimensions(:) ! dimensions of irreps, only used for relativistic

    type(readwriteblocked_tapehandle), allocatable :: &
         th_kin_real_arr(:),  th_nuc_real_arr(:),&
         th_pvsp_real_arr(:), th_pvxp_real_arr(:),& !mdf>>>
         th_ol_real_arr(:),   th_sigp_real_arr(:)   ! th_..(n_irrep), only used for relativistic
    type(readwriteblocked_tapehandle), allocatable :: &
         th_kin_imag_arr(:),  th_nuc_imag_arr(:),&
         th_pvsp_imag_arr(:), th_pvxp_imag_arr(:),& !mdf>>>
         th_ol_imag_arr(:),   th_sigp_imag_arr(:)   ! th_..(n_irrep), only used for relativistic

    character(len=5) :: ch_i_ir, ch_i_file1, ch_i_file2
    type(quadrupel_file_handle_type), allocatable :: qfh(:)
    type(readwriteblocked_tapehandle) ::  &
!        th_xc_real, th_co_real, th_ol_real, th_kin_real, th_nuc_real, th_pvsp_real,&
!        th_pvxp_real, th_sigp_real, & !<<<mdf
!        th_xc_imag, th_co_imag, th_ol_imag, th_kin_imag, th_nuc_imag, th_pvsp_imag,&
!        th_pvxp_imag, th_sigp_imag  !<<<mdf
         th_xc_real, &
         th_co_real, &
         th_rcoul_pvsp_real,&
         th_rcoul_pvxp_real,&
         th_r2_pvsp_real,&
         th_r2_pvxp_real, &
         th_ol_real, th_kin_real, th_nuc_real, th_pvsp_real,&
         th_pvxp_real, th_sigp_real, &
         th_xc_imag, th_co_imag, &
         th_rcoul_pvsp_imag,&
         th_rcoul_pvxp_imag,&
         th_r2_pvsp_imag,&
         th_r2_pvxp_imag, &
         th_ol_imag, th_kin_imag, th_nuc_imag, th_pvsp_imag,&
         th_pvxp_imag, th_sigp_imag  !<<<mdf

    logical :: restart_timer, rel
    logical, allocatable :: contract_flag(:)
    character(len=3) :: inter
    ! internal file

    logical               :: contract_3c_co !<<< mdf
    integer(kind=i4_kind) :: chklength
    integer(kind=i4_kind) :: n_if_ChFit

    !------------ Executable code ------------------------------------

    if(.not.options_integrals_on_file())then
       call error_handler(&
            & "int_send_2cob3c_spor_shutdown: I am afraid you need options_integrals_on_file")
    endif

    if ( operations_integral ) then
       if ( filesystem_is_parallel ) then
          ! send information about file sizes to other hosts
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_spor_shutdown: sending file sizes: ")
          do i_host = 1, comm_get_n_processors()
             if ( i_host .ne. my_hostindex ) then
                call comm_init_send(i_host,msgtag_int_2cob3c_filesize)
                do i_ir = 1, ubound(quadrupelfile_length,3)
                   do i_file2 = 1, ubound(quadrupelfile_length,2)
                      call commpack( &
                           quadrupelfile_length_allhosts(:,i_file2,i_ir,i_host), &
                           ubound(quadrupelfile_length,1), 1, info )
                      if (info .ne. 0) call error_handler( &
                           "int_send_2cob3c_spor_shutdown: packing of file sizes failed")
                   enddo
                enddo
                call comm_send()
             endif
          enddo

          ! receive information about file sizes from other hosts
          do while (n_missing_filesizes .gt. 1 )
             if ( output_int_loops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_shutdown: waiting for file sizes: ",inte=n_missing_filesizes)
             call comm_save_recv(comm_all_other_hosts, msgtag_int_2cob3c_filesize)

             if ( output_slaveoperations ) &
                  call write_to_output_units("int_send_2cob3c_spor_shutdown: integral_2cob_filesize")

             call int_send_2cob3c_s_rec_filesizes()
          enddo

          ! add own contrbution
          n_missing_filesizes = n_missing_filesizes - 1
          quadrupelfile_length = &
               quadrupelfile_length + quadrupelfile_length_allhosts(:,:,:,my_hostindex)

          ! free memory
          deallocate( quadrupelfile_length_allhosts, stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_shutdown: deallocate of quadrupelfile_length_allhosts failed" )

       else ! .not. ( filesystem_is_parallel )
          ! waiting for missing integrals
          do while (n_missing_quadrupels .gt. 0 )
             !? print*,"number of quadrupels missing: ",n_missing_quadrupels
             if ( output_int_loops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_shutdown: waiting for quadrupels: ",n_missing_quadrupels)
             call comm_save_recv(comm_all_other_hosts, msgtag_int_2cob3c_result)

             if ( output_slaveoperations ) &
                  call write_to_output_units("int_send_2cob3c_spor_shutdown: integral_2cob_result")

             call int_send_2cob3c_spor_receive()
          enddo

       endif

    endif

    integrals_on_file: if ( options_integrals_on_file() ) then


       if ( timer_int_idle_2cob3c(integralpar_i_int_part)%running ) then
          call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
          restart_timer = .true.
       else
          restart_timer = .false.
       endif
       call start_timer(timer_int_rewrite_2cob3c(integralpar_i_int_part))


       ! opening files, for operations integral those used in SCF-part
       ! and for operations_integraltest the same ones with "new" appended
       ! to be used for comparision with original files later on
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_shutdown: opening output files")
       if ( comm_i_am_master() ) then
          if (integralpar_2cob_kin) then
             call readwriteblocked_startwrite(trim(tmpfile("kin_real.dat")), th_kin_real)
             call readwriteblocked_startwrite(trim(tmpfile("kin_imag.dat")), th_kin_imag)
          endif
          if ( integralpar_2cob_nuc ) then
             call readwriteblocked_startwrite(trim(tmpfile("nuc_real.dat")), th_nuc_real)
             call readwriteblocked_startwrite(trim(tmpfile("nuc_imag.dat")), th_nuc_imag)
          endif
          if (integralpar_2cob_pvsp) then
             call readwriteblocked_startwrite(trim(tmpfile("pvsp_real.dat")), th_pvsp_real)
             call readwriteblocked_startwrite(trim(tmpfile("pvsp_imag.dat")), th_pvsp_imag)
          endif
          if (integralpar_2cob_pvxp) then
             call readwriteblocked_startwrite(trim(tmpfile("pvxp_real.dat")), th_pvxp_real)
             call readwriteblocked_startwrite(trim(tmpfile("pvxp_imag.dat")), th_pvxp_imag)
          endif
          !mdf>>>
          if (integralpar_2cob_pvec) then
             call readwriteblocked_startwrite(trim(tmpfile("sigp_real.dat")), th_sigp_real)
             call readwriteblocked_startwrite(trim(tmpfile("sigp_imag.dat")), th_sigp_imag)
          endif
       endif
       if (integralpar_2cob_ol.and.integralpar_relativistic) then
          call readwriteblocked_startwrite(trim(tmpfile("overlap_real.dat")), th_ol_real)
          call readwriteblocked_startwrite(trim(tmpfile("overlap_imag.dat")), th_ol_imag)
       endif
       if ( .not.integralpar_relativistic.and.integralpar_2cob_ol ) then
          unit_ol_real = openget_iounit(trim(tmpfile("overlap_real.dat")), &
               status='REPLACE', form='UNFORMATTED', action='WRITE')
          unit_ol_imag = openget_iounit(trim(tmpfile("overlap_imag.dat")), &
               status='REPLACE', form='UNFORMATTED', action='WRITE')
       endif
       if ( integralpar_3c_xc ) then
          call readwriteblocked_startwrite(trim(tmpfile("exch_real.dat")), th_xc_real)
          call readwriteblocked_startwrite(trim(tmpfile("exch_imag.dat")), th_xc_imag)
       endif
       if ( integralpar_3c_co ) then
          call readwriteblocked_startwrite(trim(tmpfile("coul_real.dat")), th_co_real)
          call readwriteblocked_startwrite(trim(tmpfile("coul_imag.dat")), th_co_imag)
       endif
       if ( integralpar_3c_rcoul_pvsp ) then
          call readwriteblocked_startwrite(trim(tmpfile("rcoul_pvsp_real.dat")), th_rcoul_pvsp_real)
          call readwriteblocked_startwrite(trim(tmpfile("rcoul_pvsp_imag.dat")), th_rcoul_pvsp_imag)
       endif
       if ( integralpar_3c_rcoul_pvxp ) then
          call readwriteblocked_startwrite(trim(tmpfile("rcoul_pvxp_real.dat")), th_rcoul_pvxp_real)
          call readwriteblocked_startwrite(trim(tmpfile("rcoul_pvxp_imag.dat")), th_rcoul_pvxp_imag)
       endif
       if ( integralpar_3c_r2_pvsp ) then
          call readwriteblocked_startwrite(trim(tmpfile("r2_pvsp_real.dat")), th_r2_pvsp_real)
          call readwriteblocked_startwrite(trim(tmpfile("r2_pvsp_imag.dat")), th_r2_pvsp_imag)
       endif
       if ( integralpar_3c_r2_pvxp ) then
          call readwriteblocked_startwrite(trim(tmpfile("r2_pvxp_real.dat")), th_r2_pvxp_real)
          call readwriteblocked_startwrite(trim(tmpfile("r2_pvxp_imag.dat")), th_r2_pvxp_imag)
       endif

       ! calculate the number of integrals n_integrals and blocksizes(n_integrals)
       ! used in concat_files(...)
       n_integrals = 0
       if ( comm_i_am_master()) then
          if (integralpar_2cob_kin) n_integrals = n_integrals + 2
          if (integralpar_2cob_nuc) n_integrals = n_integrals + 2
          if (integralpar_2cob_pvsp) n_integrals = n_integrals + 2
          if (integralpar_2cob_pvxp) n_integrals = n_integrals + 2
          if (integralpar_2cob_pvec) n_integrals = n_integrals + 2
       endif
       if ( integralpar_2cob_ol ) n_integrals = n_integrals + 2
       if ( integralpar_3c_xc ) n_integrals = n_integrals + 2
       if ( integralpar_3c_co ) n_integrals = n_integrals + 2
       if ( integralpar_3c_rcoul_pvsp ) n_integrals = n_integrals + 2
       if ( integralpar_3c_rcoul_pvxp ) n_integrals = n_integrals + 2
       if ( integralpar_3c_r2_pvsp ) n_integrals = n_integrals + 2
       if ( integralpar_3c_r2_pvxp ) n_integrals = n_integrals + 2
       allocate( blocksizes(n_integrals),contract_flag(n_integrals), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: allocate of blocksizes failed" )
       i_integral=1
       if ( comm_i_am_master()) then
          if (integralpar_2cob_kin) then
             if(integralpar_relativistic) contract_flag(i_integral)=.true.
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
          endif
          if (integralpar_2cob_nuc) then
             if(integralpar_relativistic) contract_flag(i_integral)=.true.
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
          endif
          if (integralpar_2cob_pvsp) then
             if(integralpar_relativistic) contract_flag(i_integral)=.true.
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
          endif
          if (integralpar_2cob_pvxp) then
             if(integralpar_relativistic) contract_flag(i_integral)=.true.
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
          endif
          !mdf>>>
          if (integralpar_2cob_pvec) then
             if(integralpar_relativistic) contract_flag(i_integral)=.true.
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
             blocksizes(i_integral) = 1
             i_integral = i_integral + 1
          endif
       endif
       if ( integralpar_2cob_ol ) then
          if(integralpar_relativistic) contract_flag(i_integral)=.true.
          blocksizes(i_integral) = 1
          i_integral = i_integral + 1
          blocksizes(i_integral) = 1
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_xc ) then
          blocksizes(i_integral) = borders_xc(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_xc(3,my_hostindex)
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_co ) then
          blocksizes(i_integral) = borders_ch(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_ch(3,my_hostindex)
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_rcoul_pvsp ) then
          blocksizes(i_integral) = borders_s(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_s(3,my_hostindex)
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_rcoul_pvxp ) then
          blocksizes(i_integral) = borders_s(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_s(3,my_hostindex)
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_r2_pvsp ) then
          blocksizes(i_integral) = borders_r2(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_r2(3,my_hostindex)
          i_integral = i_integral + 1
       endif
       if ( integralpar_3c_r2_pvxp ) then
          blocksizes(i_integral) = borders_r2(3,my_hostindex)
          i_integral = i_integral + 1
          blocksizes(i_integral) = borders_r2(3,my_hostindex)
          i_integral = i_integral + 1
       endif


       ! + reading of quadrupel integral files and writing
       !   them as sequential files as used in the scf part
       ! + removal of quadrupel integral files
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_shutdown: rewriting files")
       n_free_iounits = get_nbr_free_iounits()
       allocate( qfh(n_free_iounits), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "inotecordst_send_2cob3c_spor_shutdown: allocate of qfh failed" )
       n_tot_records = 0
       n_tot_records_rel = 0
       n_files_left = (n_file + 1) * n_file * symmetry_data_n_proj_irreps() / 2
       call write_to_trace_unit("int_send_2cob3c_spor_shutdown: &
            &total number of files to rewrite: ",inte=n_files_left)
       do i_ir = 1, symmetry_data_n_proj_irreps()
          write(ch_i_ir,'(i5)') i_ir
          ch_i_ir = adjustl(ch_i_ir)

          do i_file1 = 1, n_file
             write(ch_i_file1,'(i5)') i_file1
             ch_i_file1 = adjustl(ch_i_file1)

             n_file2 = i_file1

             ! in case there are more quadrupel files than free tapehandles,
             ! rewrite the last tapes to make up one tape
             i_file2_new = n_file2
             n_file2_conc = 0
             do while ( n_file2 + n_file2_conc .gt. n_free_iounits )
                !debug>>>
                call error("int_send_2cob3c_spor_shutdown: is concat supported? see cotract_flag")
                !<<<debug
                i_file2_new = i_file2_new + 1
                n_concat = n_file2 - n_free_iounits + 1 + n_file2_conc
                if (n_concat .gt. n_free_iounits-1) n_concat = n_free_iounits - 1
                call concat_files(n_file2-n_concat+1, n_file2, i_file2_new)
                n_file2_conc =  n_file2_conc + 1
                if (n_file2_conc .gt. n_free_iounits) call error_handler( &
                     "int_send_2cob3c_spor_shutdown: concating files failed" )
                n_file2 = n_file2 - n_concat
             enddo
             n_file2_new = i_file2_new

             ! opening the qu_file tapehandles
             do i_file2 = 1, n_file2
                call trace_message()
                write(ch_i_file2,'(i5)') i_file2

                ch_i_file2 = adjustl(ch_i_file2)
                qfh(i_file2)%total_length = quadrupelfile_length(i_file1,i_file2,i_ir)
                call readwriteblocked_startread( &
                     trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                     trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")),  &
                     qfh(i_file2)%th, blocklength=blocklength, &
                     total_length=qfh(i_file2)%total_length )
                call read_dimensions(qfh(i_file2))
             enddo
             do i_file2_new = i_file1+1, n_file2_new
                write(ch_i_file2,*) i_file2_new
                ch_i_file2 = adjustl(ch_i_file2)
                n_file2 = n_file2 + 1
                call readwriteblocked_startread( &
                     trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                     trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")),  &
                     qfh(n_file2)%th, blocklength=blocklength, &
                     total_length=qfh(i_file2)%total_length )
                call read_dimensions(qfh(n_file2))
             enddo

             ! sum up number of records
             do i_file2 = 1, n_file2
                n_tot_records = n_tot_records + sum(qfh(i_file2)%n_floats_if_c)
             enddo
             n_if_c=qfh(1)%n_if_c
             rel=.false.
             if(integralpar_relativistic) then
                do i_file2 = 1, n_file2
                   n_tot_records_rel = n_tot_records_rel + sum(qfh(i_file2)%n_floats_if_c_uc)
                enddo
                n_if_c=qfh(1)%n_if_c_uc
                rel=.true.
             end if

             ! do the rewriting
             if ( comm_i_am_master()) then
                if (integralpar_2cob_kin) then
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_kin_real,1, &
                              contracted=.false.)
                      enddo
                   enddo
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_kin_imag,1, &
                              contracted=.false.)
                      enddo
                   enddo
                endif
                if (integralpar_2cob_nuc) then
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_nuc_real,1, &
                              contracted=.false.)
                      enddo
                   enddo
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_nuc_imag,1, &
                              contracted=.false.)
                      enddo
                   enddo
                endif
                if (integralpar_2cob_pvsp) then
                    do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_pvsp_real,1, &
                              contracted=.false.)
                      enddo
                   enddo
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_pvsp_imag,1, &
                              contracted=.false.)
                      enddo
                   enddo
                endif
                if (integralpar_2cob_pvxp) then
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_pvxp_real,1, &
                              contracted=.false.)
                      enddo
                   enddo
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_pvxp_imag,1, &
                              contracted=.false.)
                      enddo
                   enddo
                endif
                !mdf>>>
                if (integralpar_2cob_pvec) then
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_sigp_real,1, &
                              contracted=.false.)
                      enddo
                   enddo
                   do i_if_c = 1,qfh(1)%n_if_c_uc
                      do i_file2 = 1, n_file2
                         call rewrite_blocked_integralfile(qfh(i_file2),th_sigp_imag,1, &
                              contracted=.false.)
                      enddo
                   enddo
                endif
                !<<<mdf
             endif
             if ( integralpar_2cob_ol ) then
                do i_if_c = 1,qfh(1)%n_if_c_uc
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile(qfh(i_file2),th_ol_real,1, &
                           contracted=.false.)
                   end do
                end do
                do i_if_c = 1,qfh(1)%n_if_c_uc
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile(qfh(i_file2),th_ol_imag,1, &
                           contracted=.false.)
                   end do
                end do
             endif
             if ( integralpar_3c_xc ) then
                do i_if_c = 1, qfh(1)%n_if_c
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_xc_real, &
                           borders_xc(3,my_hostindex), contracted=.true. )
                   enddo
                enddo
                do i_if_c = 1, qfh(1)%n_if_c
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_xc_imag, &
                           borders_xc(3,my_hostindex), contracted=.true. )
                   enddo
                enddo
             endif

             if(is_on(op_FitTrafo))then
                ! 3cInts uncontracted:
                contract_3c_co = .false.
                n_if_ChFit     =  qfh(1)%n_if_c_uc
             else
                ! 3cInts   contracted:
                contract_3c_co = .true.
                n_if_ChFit     =  qfh(1)%n_if_c
             endif

             if ( integralpar_3c_co ) then
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_co_real, &
                           borders_ch(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_co_imag, &
                           borders_ch(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
             endif
             if ( integralpar_3c_rcoul_pvsp ) then
                do i_if_c = 1,n_if_ChFit ! FIXME: what is the real number?
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_rcoul_pvsp_real, &
                           borders_s(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_rcoul_pvsp_imag, &
                           borders_s(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
             endif
             if ( integralpar_3c_rcoul_pvxp ) then
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_rcoul_pvxp_real, &
                           borders_s(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_rcoul_pvxp_imag, &
                           borders_s(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
             endif
             if ( integralpar_3c_r2_pvsp ) then
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_r2_pvsp_real, &
                           borders_r2(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_r2_pvsp_imag, &
                           borders_r2(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
             endif
             if ( integralpar_3c_r2_pvxp ) then
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_r2_pvxp_real, &
                           borders_r2(3,my_hostindex), contracted=contract_3c_co )
                   enddo
                enddo
                do i_if_c = 1,n_if_ChFit
                   do i_file2 = 1, n_file2
                      call rewrite_blocked_integralfile( qfh(i_file2),th_r2_pvxp_imag, &
                           borders_r2(3,my_hostindex), contracted=contract_3c_co ) !<<<mdf
                   enddo
                enddo
             endif
             ! closing the qu_file tapehandles
             do i_file2 = 1, n_file2
                call readwriteblocked_stopread(qfh(i_file2)%th,status='DELETE')
                call free_dimensions(qfh(i_file2))
             enddo

          enddo
       enddo! loop over irreps

       deallocate( qfh, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of qfh failed" )
       deallocate( blocksizes, contract_flag, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of blocksizes failed" )


       ! closing of files
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_shutdown: closing output files")
       if ( comm_i_am_master()) then
          if (integralpar_2cob_kin) then
             call readwriteblocked_stopwrite(th_kin_real)
             call readwriteblocked_stopwrite(th_kin_imag)
          endif
          if (integralpar_2cob_nuc) then
             call readwriteblocked_stopwrite(th_nuc_real)
             call readwriteblocked_stopwrite(th_nuc_imag)
          endif
          if (integralpar_2cob_pvsp) then
             call readwriteblocked_stopwrite(th_pvsp_real)
             call readwriteblocked_stopwrite(th_pvsp_imag)
          endif
          if (integralpar_2cob_pvxp) then
             call readwriteblocked_stopwrite(th_pvxp_real)
             call readwriteblocked_stopwrite(th_pvxp_imag)
          endif
          !mdf>>>
          if (integralpar_2cob_pvec) then
             call readwriteblocked_stopwrite(th_sigp_real)
             call readwriteblocked_stopwrite(th_sigp_imag)
          endif
       endif
       if (integralpar_relativistic.and. integralpar_2cob_ol ) then
          call readwriteblocked_stopwrite(th_ol_real)
          call readwriteblocked_stopwrite(th_ol_imag)
       endif
       if (.not.integralpar_relativistic.and. integralpar_2cob_ol ) then
          call returnclose_iounit(unit_ol_real,status='KEEP')
          call returnclose_iounit(unit_ol_imag,status='KEEP')
       endif
       if ( integralpar_3c_xc ) then
          call readwriteblocked_stopwrite(th_xc_real)
          call readwriteblocked_stopwrite(th_xc_imag)
       endif
       if ( integralpar_3c_co ) then
          call readwriteblocked_stopwrite(th_co_real,total_length=chklength)
          call readwriteblocked_stopwrite(th_co_imag,total_length=chklength)
       endif
       if ( integralpar_3c_rcoul_pvsp ) then
          call readwriteblocked_stopwrite(th_rcoul_pvsp_real,total_length=chklength)
          call readwriteblocked_stopwrite(th_rcoul_pvsp_imag,total_length=chklength)
       endif
       if ( integralpar_3c_rcoul_pvxp ) then
          call readwriteblocked_stopwrite(th_rcoul_pvxp_real,total_length=chklength)
          call readwriteblocked_stopwrite(th_rcoul_pvxp_imag,total_length=chklength)
       endif
       if ( integralpar_3c_r2_pvsp ) then
          call readwriteblocked_stopwrite(th_r2_pvsp_real,total_length=chklength)
          call readwriteblocked_stopwrite(th_r2_pvsp_imag,total_length=chklength)
       endif
       if ( integralpar_3c_r2_pvxp ) then
          call readwriteblocked_stopwrite(th_r2_pvxp_real,total_length=chklength)
          call readwriteblocked_stopwrite(th_r2_pvxp_imag,total_length=chklength)
       endif


!######### Am Besten die ganze Umschreiberei komplett Loeschen und nur geblockte Files verwenden !!!
       ! 'komplett loeschen' is good, Klasse !!! Weiteres ist noch zu ueberlegen

       if(integralpar_relativistic) then ! do some preparations for rewrite
          n_irrep=symmetry_data_n_proj_irreps()
          allocate(&
               & th_kin_real_arr(n_irrep),&
               & th_nuc_real_arr(n_irrep),&
               & th_ol_real_arr(n_irrep),&
               & th_pvsp_real_arr(n_irrep),&
               & th_pvxp_real_arr(n_irrep),&
               & th_sigp_real_arr(n_irrep),& !<<<mdf
               & dimensions(n_irrep),stat=status)
          if(status/=0) call error_handler(&
               "int_send_2cob3c_spor_shutdown: allocation unit_kin failed")
          allocate(&
               & th_kin_imag_arr(n_irrep),&
               & th_nuc_imag_arr(n_irrep),&
               & th_ol_imag_arr(n_irrep),&
               & th_pvsp_imag_arr(n_irrep),&
               & th_pvxp_imag_arr(n_irrep),&
               & th_sigp_imag_arr(n_irrep),&
               & stat=status)
          if(status/=0) call error_handler(&
               "int_send_2cob3c_spor_shutdown: allocation unit_kin failed")
          ! now calculate dimensions of uncontracted matrixes for every irrep
          do i_ir = 1, n_irrep
             dimensions(i_ir)=0
             do i_unique=1,n_unique_atoms
                do i_l=0,unique_atoms(i_unique)%lmax_ob
!!$                   do spin1 = 1,2
!!$                      if ((i_l.eq.0).and.(spin1.eq.1)) then
!!$                         cycle
!!$                      end if
                      dimensions(i_ir)=dimensions(i_ir)+unique_atoms(i_unique)%&
                           symadapt_spor_partner(i_ir,i_l)%n_independent_fcts*&
                           unique_atoms(i_unique)%l_ob(i_l)%n_exponents
!!$                   end do
                end do
             end do
             dimensions(i_ir)=dimensions(i_ir)*(dimensions(i_ir)+1)/2
          end do
       endif
       ! combining kin and nuc files
       if ( comm_i_am_master() .and. &
            integralpar_2cob_kin .and. integralpar_2cob_nuc .and.&
            (.not.integralpar_rel_gradients)) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_spor_shutdown: combining kin and nuc files")
          if (operations_integral) then
             ! cleaned...
          else
             call error_handler("int_send_2cob3c_spor_shutdown: wrong operations")
          endif
          call readwriteblocked_startread(trim(tmpfile("kin_real.dat")), th_kin_real)
          call readwriteblocked_startread(trim(tmpfile("kin_imag.dat")), th_kin_imag)
          call readwriteblocked_startread(trim(tmpfile("nuc_real.dat")), th_nuc_real)
          call readwriteblocked_startread(trim(tmpfile("nuc_imag.dat")), th_nuc_imag)
          if (integralpar_2cob_pvsp) then
             call readwriteblocked_startread(trim(tmpfile("pvsp_real.dat")), th_pvsp_real)
             call readwriteblocked_startread(trim(tmpfile("pvsp_imag.dat")), th_pvsp_imag)
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("pvsp_real.dat"//trim(inter))), th_pvsp_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("pvsp_imag.dat"//trim(inter))), th_pvsp_imag_arr(i_ir))
             end do
             call rewrite_one_integralfile(th_pvsp_real,th_pvsp_real_arr)
             call rewrite_one_integralfile(th_pvsp_imag,th_pvsp_imag_arr)
             call readwriteblocked_stopread(th_pvsp_real,status='DELETE')
             call readwriteblocked_stopread(th_pvsp_imag,status='DELETE')
          end if
          if (integralpar_2cob_pvxp) then
             call readwriteblocked_startread(trim(tmpfile("pvxp_real.dat")), th_pvxp_real)
             call readwriteblocked_startread(trim(tmpfile("pvxp_imag.dat")), th_pvxp_imag)
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("pvxp_real.dat"//trim(inter))), th_pvxp_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("pvxp_imag.dat"//trim(inter))), th_pvxp_imag_arr(i_ir))
             end do
             call rewrite_one_integralfile(th_pvxp_real,th_pvxp_real_arr)
             call rewrite_one_integralfile(th_pvxp_imag,th_pvxp_imag_arr)
             !mdf>>>
!!$             call readwriteblocked_stopread(th_pvxp_real,status='DELETE')
!!$             call readwriteblocked_stopread(th_pvxp_imag,status='DELETE')
             call readwriteblocked_stopread(th_pvxp_real,status='KEEP') !<<<debug
             call readwriteblocked_stopread(th_pvxp_imag,status='KEEP')
          end if
          !mdf>>>
          if (integralpar_2cob_pvec) then
             call readwriteblocked_startread(trim(tmpfile("sigp_real.dat")), th_sigp_real)
             call readwriteblocked_startread(trim(tmpfile("sigp_imag.dat")), th_sigp_imag)
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("sigp_real.dat"//trim(inter))), th_sigp_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("sigp_imag.dat"//trim(inter))), th_sigp_imag_arr(i_ir))
             end do
             call rewrite_one_integralfile(th_sigp_real,th_sigp_real_arr)
             call rewrite_one_integralfile(th_sigp_imag,th_sigp_imag_arr)
             call readwriteblocked_stopread(th_sigp_real,status='KEEP') !<<<debug
             call readwriteblocked_stopread(th_sigp_imag,status='KEEP')
          end if
          if (integralpar_2cob_ol.and.integralpar_relativistic) then
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("overlap_rel_real.dat"//trim(inter))), th_ol_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("overlap_rel_imag.dat"//trim(inter))), th_ol_imag_arr(i_ir))
             end do
             call readwriteblocked_startread(trim(tmpfile("overlap_real.dat")), th_ol_real)
             call readwriteblocked_startread(trim(tmpfile("overlap_imag.dat")), th_ol_imag)
             call rewrite_one_integralfile(th_ol_real,th_ol_real_arr)
             call rewrite_one_integralfile(th_ol_imag,th_ol_imag_arr)
             call readwriteblocked_stopread(th_ol_real,status='DELETE')
             call readwriteblocked_stopread(th_ol_imag,status='DELETE')
          end if
          if(integralpar_relativistic) then
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("kin_real.dat"//trim(inter))), th_kin_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("kin_imag.dat"//trim(inter))), th_kin_imag_arr(i_ir))
             end do
             call rewrite_one_integralfile(th_kin_real,th_kin_real_arr)
             call rewrite_one_integralfile(th_kin_imag,th_kin_imag_arr)
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("nuc_real.dat"//trim(inter))), th_nuc_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("nuc_imag.dat"//trim(inter))), th_nuc_imag_arr(i_ir))
             end do
             call rewrite_one_integralfile(th_nuc_real,th_nuc_real_arr)
             call rewrite_one_integralfile(th_nuc_imag,th_nuc_imag_arr)
          else
             do i_ir=1,n_irrep
                write(inter,'(i3)') i_ir
                inter=adjustl(inter)
                call readwriteblocked_startwrite&
                     (trim(tmpfile("kin_real.dat"//trim(inter))), th_kin_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("kin_imag.dat"//trim(inter))), th_kin_imag_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("nuc_real.dat"//trim(inter))), th_nuc_real_arr(i_ir))
                call readwriteblocked_startwrite&
                     (trim(tmpfile("nuc_imag.dat"//trim(inter))), th_nuc_imag_arr(i_ir))
             end do
             call combine_two_integralfiles(th_kin_real,th_nuc_real,unit_kin_nuc_real)
             call combine_two_integralfiles(th_kin_imag,th_nuc_imag,unit_kin_nuc_imag)
             call returnclose_iounit(unit_kin_nuc_real,status='KEEP')
             call returnclose_iounit(unit_kin_nuc_imag,status='KEEP')
          end if

          call readwriteblocked_stopread(th_kin_real,status='DELETE')
          call readwriteblocked_stopread(th_kin_imag,status='DELETE')
          call readwriteblocked_stopread(th_nuc_real,status='DELETE')
          call readwriteblocked_stopread(th_nuc_imag,status='DELETE')
       endif
       if (integralpar_relativistic) then
          deallocate(&
               & th_kin_real_arr,&
               & th_nuc_real_arr,&
               & th_ol_real_arr,&
               & th_pvsp_real_arr,&
               & th_pvxp_real_arr,&
               & th_sigp_real_arr,& !<<<mdf
               & dimensions,stat=status)
          if(status/=0) call error_handler(&
               "int_send_2cob3c_spor_shutdown: deallocating unit_kin_arr failed")
          deallocate(&
               & th_kin_imag_arr,&
               & th_nuc_imag_arr,&
               & th_ol_imag_arr,&
               & th_pvsp_imag_arr,&
               & th_pvxp_imag_arr,&
               & th_sigp_imag_arr,&
               & stat=status)
          if(status/=0) call error_handler(&
               "int_send_2cob3c_spor_shutdown: deallocating unit_kin_arr failed")
       endif


       ! deallocation of help arrays
       deallocate( ua_fileindex, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of ua_fileindex failed" )
       deallocate( quadrupelfile_length, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of quadrupelfile_length failed" )


       call stop_timer(timer_int_rewrite_2cob3c(integralpar_i_int_part))
       if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))



    else ! .not. options_integrals_on_file()
       ABORT('FIXME: ints in mem')
#ifdef SO_INTEGRALS_ON_FILE
       ! the integrals are already in their final location

       do i_ir = 1,symmetry_data_n_proj_irreps()
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)

                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)

                ind1_exp1_ua2 => metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                      do i_ua2 = 1, i_ua1
                         deallocate( ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2, &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral l2 failed" )
                      enddo
                   enddo
                enddo
                deallocate( ind1_exp1_ua2, stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )

             enddo! i_l1
             deallocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1, stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_shutdown: allocate of metaindex_of_first_integral l1 failed" )
          enddo! i_ua1
          deallocate( metaindex_of_first_integral(i_ir)%ua1, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral ua1 failed" )
       enddo
       deallocate( metaindex_of_first_integral, &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral failed" )


       do i_ir = 1,symmetry_data_n_proj_irreps()
          do i_ua1 = 1, N_unique_atoms
             ua1 => unique_atoms(i_ua1)
             do i_l1 = 0, ua1%lmax_ob
                uab1 => ua1%l_ob(i_l1)

                uap1 => ua1%symadapt_spor_partner(i_ir,i_l1)
                ind1_exp1_ua2 => metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
                do i_ind1 = 1, uap1%N_independent_fcts
                   do i_exp1 = 1, uab1%N_exponents
                      do i_ua2 = 1, i_ua1
                         deallocate( ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2, &
                              stat=status )
                         if ( status .ne. 0 ) call error_handler( &
                              "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral_rel l2 failed" )
                      enddo
                   enddo
                enddo
                deallocate( ind1_exp1_ua2, stat=status )
                if ( status .ne. 0 ) call error_handler( &
                     "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral_rel ind1_exp1_ua2 failed" )

             enddo! i_l1
             deallocate( metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1, stat=status )
             if ( status .ne. 0 ) call error_handler( &
                  "int_send_2cob3c_spor_shutdown: allocate of metaindex_of_first_integral_rel l1 failed" )
          enddo! i_ua1
          deallocate( metaindex_of_first_integral_rel(i_ir)%ua1, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral_rel ua1 failed" )
       enddo
       deallocate( metaindex_of_first_integral_rel, &
            stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_spor_shutdown: deallocate of metaindex_of_first_integral_rel failed" )

#endif
    endif integrals_on_file


    ! deallocation of help arrays
    deallocate( borders_ch, borders_s, borders_r2, borders_xc, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "int_send_2cob3c_spor_shutdown: deallocate of borders failed" )

    if ( output_int_loops ) call write_to_output_units("int_send_2cob3c_spor_shutdown: done")


  contains


    subroutine write_dimensions(qfh)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(quadrupel_file_handle_type), intent(inout)  :: qfh
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      real(kind=r8_kind) :: real_n_if_c(1)
      real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
      !------------ Executable code ----------------------------
      allocate( real_n_floats_if_c(qfh%n_if_c), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(qfh%n_if_c,r8_kind)
      real_n_floats_if_c = real(qfh%n_floats_if_c,r8_kind)
      call readwriteblocked_write(real_n_if_c, qfh%th)
      call readwriteblocked_write(real_n_floats_if_c,qfh%th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: deallocating at write_dimensions() failed")
      allocate( real_n_floats_if_c(qfh%n_if_c_uc), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(qfh%n_if_c,r8_kind)
      real_n_floats_if_c = real(qfh%n_floats_if_c_uc,r8_kind)
      call readwriteblocked_write(real_n_if_c, qfh%th)
      call readwriteblocked_write(real_n_floats_if_c,qfh%th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: deallocating at write_dimensions() failed")
    end subroutine write_dimensions


    subroutine read_dimensions(qfh)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(quadrupel_file_handle_type), intent(inout)  :: qfh
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      real(kind=r8_kind) :: real_n_if_c(1)
      real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
      !------------ Executable code ----------------------------
      call readwriteblocked_read(real_n_if_c,qfh%th)
      qfh%n_if_c = nint(real_n_if_c(1),i4_kind)
      allocate( real_n_floats_if_c(qfh%n_if_c), &
           qfh%n_floats_if_c(qfh%n_if_c), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at read_dimensions() failed")
      call readwriteblocked_read(real_n_floats_if_c,qfh%th)
      qfh%n_floats_if_c = nint(real_n_floats_if_c,i4_kind)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: deallocating at read_dimensions() failed")
      call readwriteblocked_read(real_n_if_c,qfh%th)
      qfh%n_if_c_uc = nint(real_n_if_c(1),i4_kind)
      allocate( real_n_floats_if_c(qfh%n_if_c_uc), &
           qfh%n_floats_if_c_uc(qfh%n_if_c_uc), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at read_dimensions() failed")
      call readwriteblocked_read(real_n_floats_if_c,qfh%th)
      qfh%n_floats_if_c_uc = nint(real_n_floats_if_c,i4_kind)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: deallocating at read_dimensions() failed")
    end subroutine read_dimensions


    subroutine sum_dimensions(qfh_in,qfh_sum)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(quadrupel_file_handle_type), intent(in)  :: qfh_in(:)
      type(quadrupel_file_handle_type), intent(inout)  :: qfh_sum
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status,i_f
      !------------ Executable code ----------------------------
      qfh_sum%n_if_c = qfh_in(1)%n_if_c
      allocate( qfh_sum%n_floats_if_c(qfh_sum%n_if_c), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at sum_dimensions() failed")
      qfh_sum%n_floats_if_c = 0
      do i_f = 1, ubound(qfh_in,1)
         qfh_sum%n_floats_if_c = qfh_sum%n_floats_if_c + qfh(i_f)%n_floats_if_c
      enddo
      qfh_sum%n_if_c_uc = qfh_in(1)%n_if_c_uc
      allocate( qfh_sum%n_floats_if_c_uc(qfh_sum%n_if_c_uc), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: allocating at sum_dimensions() failed")
      qfh_sum%n_floats_if_c_uc = 0
      do i_f = 1, ubound(qfh_in,1)
         qfh_sum%n_floats_if_c_uc = qfh_sum%n_floats_if_c_uc + qfh(i_f)%n_floats_if_c_uc
      enddo
    end subroutine sum_dimensions


    subroutine free_dimensions(qfh)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(quadrupel_file_handle_type), intent(inout) :: qfh
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      !------------ Executable code ----------------------------
      deallocate( qfh%n_floats_if_c,stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at free_dimensions() failed")
      deallocate( qfh%n_floats_if_c_uc,stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at free_dimensions() failed")
    end subroutine free_dimensions


    subroutine concat_files(i_file2_from,i_file2_to,i_file2_new)
      use error_module !<<< debug
      implicit none
      !------------ Declaration of formal parameters -----------
      integer(kind=i4_kind), intent(in) :: i_file2_from,i_file2_to,i_file2_new
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_integral,i_file2,i_f,n_f, n_if_c
      type(quadrupel_file_handle_type) :: qfh_new
      !------------ Executable code ----------------------------

      call warn("is23m/concat_files: ENTERED, ENTERED, ENTERED") !<<< debug

      ! opening the qu_file tapehandles to read from
      i_f = 0
      do i_file2 = i_file2_from, i_file2_to
         call trace_message()
         i_f = i_f + 1
         write(ch_i_file2,*) i_file2
         ch_i_file2 = adjustl(ch_i_file2)
         qfh(i_f)%total_length = quadrupelfile_length(i_file1,i_file2,i_ir)
         call readwriteblocked_startread( &
              trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
              trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")),  &
              qfh(i_f)%th, blocklength=blocklength, &
              total_length=qfh(i_f)%total_length )
         call read_dimensions(qfh(i_f))
      enddo
      n_f = i_f
      ! opening the tapehandle to write to
      write(ch_i_file2,*) i_file2_new
      ch_i_file2 = adjustl(ch_i_file2)
      call readwriteblocked_startwrite( &
           trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
           trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")),  &
           qfh_new%th, blocklength=blocklength )
      ! calculate and write size information for new file
      call sum_dimensions(qfh(1:n_f),qfh_new)
      call write_dimensions(qfh_new)
      ! do the rewriting
      do i_integral = 1, n_integrals
         if(integralpar_relativistic) then
            n_if_c=qfh_new%n_if_c_uc
         else
            n_if_c=qfh_new%n_if_c
         end if
         ! real
         do i_if_c = 1, n_if_c
            do i_f = 1, n_f
               call rewrite_blocked_integralfile(qfh(i_f), &
                    qfh_new%th,blocksizes(i_integral),contracted=contract_flag(i_integral))
            enddo
         enddo
         ! immaginary
         do i_if_c = 1, n_if_c
            do i_f = 1, n_f
               call rewrite_blocked_integralfile(qfh(i_f), &
                    qfh_new%th,blocksizes(i_integral),contracted=contract_flag(i_integral))
            enddo
         enddo
      enddo
      ! closing the qu_file tapehandles and removing the files
      do i_f = 1, n_f
         call readwriteblocked_stopread(qfh(i_f)%th,status='DELETE')
         call free_dimensions(qfh(i_f))
      enddo
      ! closing tapehandle that has been written to
      call readwriteblocked_stopwrite(qfh_new%th,total_length=qfh_new%total_length)
      call free_dimensions(qfh_new)
    end subroutine concat_files



    subroutine rewrite_blocked_integralfile(qfh_in,th_out,bordersize, contracted)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(quadrupel_file_handle_type), intent(inout) :: qfh_in
      type(readwriteblocked_tapehandle), intent(inout) :: th_out
      integer(kind=i4_kind), intent(in) :: bordersize
      logical :: contracted
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: n_floats,status
      real(kind=r8_kind), allocatable :: buffer(:)

      if(contracted) then
         n_floats = qfh_in%n_floats_if_c(i_if_c) * bordersize
      else
         n_floats = qfh_in%n_floats_if_c_uc(i_if_c) * bordersize
      end if
      allocate(buffer(n_floats), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: rewrite_blocked_integralfile: allocating failed")
      call readwriteblocked_read(buffer,qfh_in%th)
      call readwriteblocked_write(buffer,th_out)
      deallocate(buffer, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_shutdown: rewrite_blocked_integralfile: deallocating failed")
    end subroutine rewrite_blocked_integralfile


    subroutine combine_two_integralfiles(th1_in,th2_in,unit_out)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(readwriteblocked_tapehandle), intent(inout) :: th1_in,th2_in
      integer(kind=i4_kind), intent(in) :: unit_out
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_rec, n_rec, status, length, n_rec_left
      real(kind=r8_kind), allocatable :: buffer(:,:)
      !------------ Executable code ----------------------------
      length = readwriteblocked_blocklength()
      allocate(buffer(2,length), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_shutdown: combine_two_integralfiles: allocating failed")
      n_rec_left = n_tot_records_rel
      do while (n_rec_left .gt. 0)
         if (n_rec_left .gt. length) then
            n_rec_left = n_rec_left - length
            n_rec = length
         else
            n_rec = n_rec_left
            n_rec_left = 0
         endif
         call readwriteblocked_read(buffer(1,1:n_rec),th1_in)
         call readwriteblocked_read(buffer(2,1:n_rec),th2_in)
         do i_rec = 1, n_rec
            write(unit_out,ERR=99) buffer(:,i_rec)
         enddo
      enddo
      deallocate(buffer, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_shutdown: combine_two_integralfiles: deallocating failed")
      return
99    call error_handler("int_send_2cob3c_shutdown: combine_two_integralfiles: writing failed")
    end subroutine combine_two_integralfiles

    subroutine rewrite_one_integralfile(th_in,th_out)
      implicit none
      ! note: in the relativistic case we have one file for every irrep
      !------------ Declaration of formal parameters -----------
      type(readwriteblocked_tapehandle), intent(inout) :: th_in
      type(readwriteblocked_tapehandle), intent(inout) :: th_out(n_irrep)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_ir, status, length
      real(kind=r8_kind), allocatable :: buffer(:)
      !------------ Executable code ----------------------------
      do i_ir=1,n_irrep
         allocate(buffer(dimensions(i_ir)), stat=status)
         if (status .ne. 0) call error_handler( &
              "int_send_2cob3c_shutdown: rewrite_one_integralfile: allocating failed")
         call readwriteblocked_read(buffer(1:dimensions(i_ir)),th_in)
         call readwriteblocked_write(buffer(1:dimensions(i_ir)),th_out(i_ir))
         call readwriteblocked_stopwrite(th_out(i_ir),total_length=length)
         deallocate(buffer, stat=status)
         if (status .ne. 0) call error_handler( &
              "int_send_2cob3c_shutdown: rewrite_one_integralfile: deallocating failed")
      end do
    end subroutine rewrite_one_integralfile


    subroutine trace_message()
      n_files_left = n_files_left - 1
      if (mod(n_files_left,50) .eq. 0) &
           call write_to_trace_unit("int_send_2cob3c_spor_shutdown: &
           &number of files left: ",inte=n_files_left)
    end subroutine trace_message


  end subroutine int_send_2cob3c_spor_shutdown
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cob3c_spor_send
    !  Purpose: sending of contracted and symmetry adapted
    !    integrals from int_data_2cob3c_module.
    !    When packing the integrals, the mapping to the scf-
    !    metaindex is also done.
    !    Two center integrals are send to the master and
    !    three center integrals are distributed over all hosts
    !    splitting them over the fitfunction metaindex.
    !    Integrals to be hold on present hosts are directly
    !    written to file or stored in memory in integralstore_module.
    !    In case the integral scratch is parallel and files are used,
    !    everything is directley written to files.
    !  called by: integral_calc_quad_2cob3c
    use type_module, RK=>r8_kind
    use int_data_2cob3c_module
    use spin_orbit_module, only: is_on,op_FitTrafo
    use int_send_aux_module
    use quadrupel_module, only: quadrupel_transpose, operator(.eq.)
    use symm_adapt_int,   only: hconjug
    !** End of interface ***************************************

    type(quadrupel_type) :: tquadrupel

    if ( options_integrals_on_file() ) then

       call int_send_2cob3c_spor_send_file()

       if(is_on(op_FitTrafo))then
          !
          ! store data from int_data_2cob3c_module
          !
          call int_write_sa_int(sa_2c_ints(:, int_sa_sigp), quadrupel)

          !---------------------------------
          tquadrupel = quadrupel_transpose(quadrupel)

          if(quadrupel.eq.tquadrupel)then
             !
             ! diagonal quadrupel:
             !
             call int_write_sa_int(sa_alphap(:, 1), quadrupel, EXT="bLS")
          else
             !
             ! not a diagonal quadrupel:
             !
             call int_write_sa_int(sa_alphap(:, 1), quadrupel, EXT="bLS")
             call hconjug(sa_alphap(:, 2))
             call int_write_sa_int(sa_alphap(:, 2), tquadrupel, EXT="bLS")
          endif
          !---------------------------------
       endif
    else
       call int_send_2cob3c_spor_send_mem()
    endif
  end subroutine int_send_2cob3c_spor_send
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cob3c_spor_send_file
    !  Purpose: sending of contracted and symmetry adapted
    !    integrals from int_data_2cob3c_module.
    !    When packing the integrals, the mapping to the scf-
    !    metaindex is also done.
    !    Two center integrals are send to the master and
    !    three center integrals are distributed over all hosts
    !    splitting them over the fitfunction metaindex.
    !  called by: integral_calc_quad_2cob3c
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use int_data_2cob3c_module, only: ua1                                      &
                                    , ua2                                      &
                                    , n_c1                                     &
                                    , n_c2                                     &
                                    , n_exp1                                   &
                                    , n_exp2                                   &
                                    , quadrupel                                &
                                    , diagonal                                 &
                                    , int_sa_co                                &
                                    , int_sa_kin                               &
                                    , int_sa_nuc                               &
                                    , int_sa_ol                                &
                                    , int_sa_pvsp                              &
                                    , int_sa_pvxp                              &
                                    , int_sa_r2_pvsp                           &
                                    , int_sa_r2_pvxp                           &
                                    , int_sa_rcoul_pvsp                        &
                                    , int_sa_rcoul_pvxp                        &
                                    , int_sa_sigp                              &
                                    , int_sa_xc                                &
                                    , sa_2c_ints                               &
                                    , sa_3c_ints

    use symm_adapt_int,         only: Write2cIntBlock                          &
                                    , Write3cIntBlock                          &
                                    , Pack2cIntBlock                           &
                                    , Pack3cIntBlock
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i_host, i_ir, n_ir, n_if_c_max, &
         n_if_c_uc_max
    type(readwriteblocked_tapehandle),target :: th
    type(readwriteblocked_tapehandle),pointer :: th_pointer
    integer(kind=i4_kind), allocatable :: &
         n_floats_if_c(:,:), n_floats_if_c_max(:), n_if_c(:),&
         n_floats_if_c_uc(:,:), n_floats_if_c_uc_max(:), n_if_c_uc(:)
    ! n_floats_if_c(i_if_c,i_ir), n_floats_if_c_max(i_ir), n_if_c(i_ir)
    ! n_floats_if_c(i_if_c_uc,i_ir), n_floats_if_c_uc_max(i_ir), n_if_c_uc(i_ir)

    n_ir = symmetry_data_n_proj_irreps()


!!$    if ( quadrupel%l1 == 0 ) then
!!$       n_spin1 = 1
!!$    else
!!$       n_spin1 = 2
!!$    endif
!!$    if ( quadrupel%l2 == 0 ) then
!!$       n_spin2 = 1
!!$    else
!!$       n_spin2 = 2
!!$    endif

    call calc_dimensions()

    hosts: do i_host = first_host, last_host

       parallel_filesystem: if ( filesystem_is_parallel ) then
          call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))
          do i_ir = 1, n_ir

             call readwriteblocked_startwrite( &
                  trim(filename_tmpdir(i_host)) // "/" // trim(quadrupel_filename(i_ir,quadrupel)), &
                  th, blocklength=blocklength )
             th_pointer=> th
             call write_dimensions(th,i_ir)

             if ( i_host .eq. 1) then
                if ( integralpar_2cob_kin ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing kin")
                   call Write2cIntBlock(&
                        & sa_2c_ints(i_ir, int_sa_kin),&
                        & diagonal,&
                        & th_pointer)
                endif

                if ( integralpar_2cob_nuc ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing nuc")
                   call Write2cIntBlock(&
                        & sa_2c_ints(i_ir, int_sa_nuc),&
                        & diagonal,&
                        & th_pointer)
                endif
                if ( integralpar_2cob_pvsp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing pvsp")
                   call Write2cIntBlock(&
                        & sa_2c_ints(i_ir, int_sa_pvsp),&
                        & diagonal,&
                        & th_pointer)
                endif
                if ( integralpar_2cob_pvxp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing pvxp")
                   call Write2cIntBlock(&
                        & sa_2c_ints(i_ir, int_sa_pvxp),&
                        & diagonal,&
                        & th_pointer)
                endif
                if ( integralpar_2cob_pvec ) then
                   if ( output_int_deeploops ) call write_to_output_units( & !<<<debug
                        "int_send_2cob3c_spor_send: writing sigp")
                   call Write2cIntBlock(&
                        & sa_2c_ints(i_ir, int_sa_sigp),&
                        & diagonal,&
                        & th_pointer)
                endif
             endif
             if ( integralpar_2cob_ol ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing overlap")
                call Write2cIntBlock(&
                     & sa_2c_ints(i_ir, int_sa_ol),&
                     & diagonal,&
                     & th_pointer)
             endif
             if ( integralpar_3c_xc ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing exchange")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_xc),&
                     & borders_xc(:,i_host),&
                     & diagonal,&
                     & th_pointer)
             endif
#ifdef FPP_DEBUG
             print *, 'isend: q=', quadrupel%ua1,quadrupel%l1,quadrupel%ua2,quadrupel%l2
#endif
             if ( integralpar_3c_co ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing coulomb")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_co),&
                     & borders_ch(:,i_host),&
                     & diagonal,&
                     & th_pointer)
#ifdef FPP_DEBUG
             print *, 'isend: coul(',i_ir,')=',&
                  & sum(sa_3c_ints(i_ir, int_sa_co)%re),&
                  & sum(sa_3c_ints(i_ir, int_sa_co)%im)
#endif
             endif
             if ( integralpar_3c_rcoul_pvsp ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing rcoul_pvsp")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_rcoul_pvsp),&
                     & borders_s(:,i_host),&
                     & diagonal,&
                     & th_pointer)
#ifdef FPP_DEBUG
                print *, 'isend: rcoul_pvsp(',i_ir,')=',&
                     & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvsp)%re),&
                     & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvsp)%im)
#endif
             endif
             if ( integralpar_3c_rcoul_pvxp ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing rcoul_pvxp")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_rcoul_pvxp),&
                     & borders_s(:,i_host),&
                     & diagonal,&
                     & th_pointer)
#ifdef FPP_DEBUG
                print *, 'isend: rcoul_pvxp(',i_ir,')=',&
                     & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvxp)%re),&
                     & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvxp)%im)
#endif
             endif
             if ( integralpar_3c_r2_pvsp ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing r2_pvsp")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_r2_pvsp),&
                     & borders_r2(:,i_host),&
                     & diagonal,&
                     & th_pointer)
#ifdef FPP_DEBUG
               print *, 'isend: r2_pvsp(',i_ir,')=',&
                    & sum(sa_3c_ints(i_ir, int_sa_r2_pvsp)%re),&
                    & sum(sa_3c_ints(i_ir, int_sa_r2_pvsp)%im)
#endif
             endif
             if ( integralpar_3c_r2_pvxp ) then
                if ( output_int_deeploops ) call write_to_output_units( &
                     "int_send_2cob3c_spor_send: writing r2_pvsp")
                call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_r2_pvxp),&
                     & borders_r2(:,i_host),&
                     & diagonal,&
                     & th_pointer)
#ifdef FPP_DEBUG
                print *, 'isend: r2_pvxp(',i_ir,')=',&
                     & sum(sa_3c_ints(i_ir, int_sa_r2_pvxp)%re),&
                     & sum(sa_3c_ints(i_ir, int_sa_r2_pvxp)%im)
#endif
             endif

             ! quadrupelfile_length is defined HERE:
             call readwriteblocked_stopwrite(th, total_length= &
                  quadrupelfile_length_allhosts( &
                  index_ua_l(quadrupel%ua1,quadrupel%l1), &
                  index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir, i_host ) )

          enddo
          call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))

       else ! not parallel_filesystem

          write_own_part: if ( i_host .eq. my_hostindex) then
             call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))
             do i_ir = 1, n_ir
                !? print*,">>> int_send_2cob3c_spor_send_file: i_ir",i_ir
                call readwriteblocked_startwrite( &
                     trim(tmpfile(quadrupel_filename(i_ir,quadrupel))), &
                     th, blocklength=blocklength )
                th_pointer=> th
                call write_dimensions(th,i_ir)
                if ( i_host .eq. 1 ) then
                   if ( integralpar_2cob_kin ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: writing kin")
!!$                      call write_integralfile_2c(symadapt_int_2cob_kin_p(i_ir,:,:), &
!!$                           contracted=.false. )
                      call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_kin),diagonal,th_pointer)
                   endif
                   if ( integralpar_2cob_nuc ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: writing nuc")
!!$                      call write_integralfile_2c(symadapt_int_2cob_nuc_p(i_ir,:,:), &
!!$                           contracted=.false. )
                      call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_nuc),diagonal,th_pointer)
                   endif
                   if ( integralpar_2cob_pvsp ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: writing pvsp")
!!$                      call write_integralfile_2c(symadapt_int_2cob_pvsp_p(i_ir,:,:), &
!!$                           contracted=.false. )
                      call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_pvsp),diagonal,th_pointer)
                   endif
                   if ( integralpar_2cob_pvxp ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: writing pvxp")
!!$                      call write_integralfile_2c(symadapt_int_2cob_pvxp_p(i_ir,:,:), &
!!$                           contracted=.false. )
                      call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_pvxp),diagonal,th_pointer)
                   endif
                   !mdf>>>
                   if ( integralpar_2cob_pvec ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: writing sigp")
!!$                      call write_integralfile_2c(symadapt_int_2cob_sigp_p(i_ir,:,:), &
!!$                           contracted=.false. )
                      call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_sigp),diagonal,th_pointer)
                   endif
                end if
                if ( integralpar_2cob_ol ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing overlap")
!!$                   call write_integralfile_2c(symadapt_int_2cob_ol_p(i_ir,:,:), &
!!$                        contracted=.false. )
                   call Write2cIntBlock(sa_2c_ints(i_ir, int_sa_ol),diagonal,th_pointer)
                endif
                if ( integralpar_3c_xc ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing coulomb")
!!$                   call write_integralfile_3c( symadapt_int_3c_xc_p(i_ir,:,:), &
!!$                        borders_xc(:,i_host) )
                   call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_xc),&
                     & borders_xc(:,i_host),&
                     & diagonal,&
                     & th_pointer)
                endif
                if ( integralpar_3c_co ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing exchange")
!!$                   call write_integralfile_3c( symadapt_int_3c_co_p(i_ir,:,:), &
!!$                        borders_ch(:,i_host) )
                   call Write3cIntBlock(&
                     & sa_3c_ints(i_ir, int_sa_co),&
                     & borders_ch(:,i_host),&
                     & diagonal,&
                     & th_pointer)
                endif
#ifdef FPP_DEBUG
                print *, 'isend: q=', quadrupel%ua1,quadrupel%l1,quadrupel%ua2,quadrupel%l2
#endif
                if ( integralpar_3c_rcoul_pvsp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing rcoul_pvsp")
                   call Write3cIntBlock(&
                        & sa_3c_ints(i_ir, int_sa_rcoul_pvsp),&
                        & borders_s(:,i_host),&
                        & diagonal,&
                        & th_pointer)
#ifdef FPP_DEBUG
                   print *, 'isend: rcoul_pvsp(',i_ir,')=',&
                        & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvsp)%re),&
                        & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvsp)%im)
#endif
                endif
                if ( integralpar_3c_rcoul_pvxp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing rcoul_pvxp")
                   call Write3cIntBlock(&
                        & sa_3c_ints(i_ir, int_sa_rcoul_pvxp),&
                        & borders_s(:,i_host),&
                        & diagonal,&
                        & th_pointer)
#ifdef FPP_DEBUG
                   print *, 'isend: rcoul_pvxp(',i_ir,')=',&
                        & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvxp)%re),&
                        & sum(sa_3c_ints(i_ir, int_sa_rcoul_pvxp)%im)
#endif
                endif
                if ( integralpar_3c_r2_pvsp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing r2_pvsp")
                   call Write3cIntBlock(&
                        & sa_3c_ints(i_ir, int_sa_r2_pvsp),&
                        & borders_r2(:,i_host),&
                        & diagonal,&
                        & th_pointer)
#ifdef FPP_DEBUG
                   print *, 'isend: r2_pvsp(',i_ir,')=',&
                        & sum(sa_3c_ints(i_ir, int_sa_r2_pvsp)%re),&
                        & sum(sa_3c_ints(i_ir, int_sa_r2_pvsp)%im)
#endif
                endif
                if ( integralpar_3c_r2_pvxp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: writing r2_pvsp")
                   call Write3cIntBlock(&
                        & sa_3c_ints(i_ir, int_sa_r2_pvxp),&
                        & borders_r2(:,i_host),&
                        & diagonal,&
                        & th_pointer)
#ifdef FPP_DEBUG
                   print *, 'isend: r2_pvxp(',i_ir,')=',&
                        & sum(sa_3c_ints(i_ir, int_sa_r2_pvxp)%re),&
                        & sum(sa_3c_ints(i_ir, int_sa_r2_pvxp)%im)
#endif
                endif
                call readwriteblocked_stopwrite(th, total_length= &
                     quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                     index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ) )

             enddo
             call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
             n_missing_quadrupels = n_missing_quadrupels - 1
          else ! not write_own_part
             call start_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
             call comm_init_send(i_host,msgtag_int_2cob3c_result)
             call quadrupel_pack(quadrupel)
             call pack_dimensions()
             do i_ir = 1, n_ir
                if ( i_host .eq. 1) then
                   if ( integralpar_2cob_kin ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: packing kin")
                      call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_kin),diagonal)
                   endif
                   if ( integralpar_2cob_nuc ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: packing nuc")
                      call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_nuc),diagonal)
                   endif
                   if ( integralpar_2cob_pvsp ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: packing pvsp")
                      call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_pvsp),diagonal)
                   endif
                   if ( integralpar_2cob_pvxp ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: packing pvxp")
                      call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_pvxp),diagonal)
                   endif
                   if ( integralpar_2cob_pvec ) then
                      if ( output_int_deeploops ) call write_to_output_units( &
                           "int_send_2cob3c_spor_send: packing sigp")
                      call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_sigp),diagonal)
                   endif
                endif
                if ( integralpar_2cob_ol ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing overlap")
                   call Pack2cIntBlock(sa_2c_ints(i_ir, int_sa_ol),diagonal)
                endif
                if ( integralpar_3c_xc ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing coulomb")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_xc),borders_xc(:,i_host),diagonal)
                endif
                if ( integralpar_3c_co ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing exchange")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_co),borders_ch(:,i_host),diagonal)
                endif
                if ( integralpar_3c_rcoul_pvsp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing rcoul_pvsp")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_rcoul_pvsp),borders_s(:,i_host),diagonal)
                endif
                if ( integralpar_3c_rcoul_pvxp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing rcoul_pvxp")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_rcoul_pvxp),borders_s(:,i_host),diagonal)
                endif
                if ( integralpar_3c_r2_pvsp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing r2_pvsp")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_r2_pvsp),borders_r2(:,i_host),diagonal)
                endif
                if ( integralpar_3c_r2_pvxp ) then
                   if ( output_int_deeploops ) call write_to_output_units( &
                        "int_send_2cob3c_spor_send: packing r2_pvxp")
                   call Pack3cIntBlock(sa_3c_ints(i_ir, int_sa_r2_pvxp),borders_r2(:,i_host),diagonal)
                endif
             end do
             call stop_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
             call start_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
             call comm_send()
             call stop_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
          end if write_own_part

       end if parallel_filesystem

    end do hosts

    call free_dimensions()

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cob3_send: done")


  contains


    subroutine calc_dimensions()
      !** End of interface *************************************
      implicit none
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,n_if1_max, &
           n_if1,n_if2,n_floats,i_if_c,status !!$,i_spin1,i_spin2,nn_spin2
      !------------ Executable code ----------------------------


      ! We need information both for integrals that have
      ! contracted exponents and uncontracted exponents

      ! Contracted exponents:
!!$      n_if1_max = maxval(ua1%symadapt_spor_partner(:,quadrupel%l1,3-n_spin1:2)%N_independent_fcts)
      n_if1_max = maxval(ua1%symadapt_spor_partner(:,quadrupel%l1)%N_independent_fcts)
      n_if_c_max = n_if1_max * n_c1 !!$* n_spin1
!!$      print*,">>> calc_dimensions: n_if1_max",n_if1_max
!!$      print*,">>> calc_dimensions: n_c1",n_c1
!!$      print*,">>> calc_dimensions: n_spin1",n_spin1
      !allocate( n_floats_if_c(n_if_c_max,n_ir), n_if_c(n_ir), &
      !     n_floats_if_c_max(n_ir), stat=status)
      !if (status .ne. 0) call error_handler( &
      !     "int_send_2cob3c_spor_send: allocating at calc_dimensions() failed")
!!$      print*,">>> calc_dimensions:  n_if_c_max",n_if_c_max," n_ir ",n_ir
      allocate( n_floats_if_c(n_if_c_max,n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating of n_floats_if_c at calc_dimensions() failed")
      allocate(n_if_c(n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating of n_if_c at calc_dimensions() failed")
      allocate( n_floats_if_c_max(n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating of n_floats_if_c_max at calc_dimensions() failed")
      do i_ir = 1,n_ir
         i_if_c = 1
!!$         do i_spin1 = 1, 2
!!$            if ((n_spin1.eq.1).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1)%N_independent_fcts
!!$            if (diagonal) then
!!$               nn_spin2 = i_spin1
!!$            else
!!$               nn_spin2 = 2
!!$            endif
         do i_if1 = 1, n_if1
            do i_c1 = 1, n_c1
               n_floats = 0
!!$                  do i_spin2 = 1, nn_spin2
!!$                     if ((n_spin2.eq.1).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
               if ( diagonal ) then
                  n_if2 = i_if1
               else
!!$                  n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                       N_independent_fcts
                  n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2)% &
                       N_independent_fcts
               endif
               do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
                  if ( diagonal .and. i_if1 .eq. i_if2 ) then
                     n_floats = n_floats + i_c1
                  else
                     n_floats = n_floats + n_c2
                  endif
               enddo
!!$                  enddo
               n_floats_if_c(i_if_c,i_ir) = n_floats
               i_if_c = i_if_c + 1
            enddo
!!$         enddo
         enddo
      n_if_c(i_ir) = i_if_c - 1
      if ( n_if_c(i_ir) .gt. 0 ) then
         n_floats_if_c_max(i_ir) = maxval(n_floats_if_c(1:n_if_c(i_ir),i_ir))
      else
         n_floats_if_c_max(i_ir) = 0
      endif
   enddo

      ! Uncontracted exponents:
!!$      n_if1_max = maxval(ua1%symadapt_spor_partner(:,quadrupel%l1,3-n_spin1:2)%N_independent_fcts)
      n_if1_max = maxval(ua1%symadapt_spor_partner(:,quadrupel%l1)%N_independent_fcts)
      n_if_c_uc_max = n_if1_max * n_exp1 !!$ * n_spin1
      allocate( n_floats_if_c_uc(n_if_c_uc_max,n_ir), n_if_c_uc(n_ir), &
           n_floats_if_c_uc_max(n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating at calc_dimensions() failed")
      do i_ir = 1,n_ir
         i_if_c = 1
!!$         do i_spin1 = 1, 2
!!$            if ((n_spin1.eq.1).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1)%N_independent_fcts
!!$            if (diagonal) then
!!$               nn_spin2 = i_spin1
!!$            else
!!$               nn_spin2 = 2
!!$            endif
         do i_if1 = 1, n_if1
            do i_c1 = 1, n_exp1
               n_floats = 0
!!$                  do i_spin2 = 1, nn_spin2
!!$                     if ((n_spin2.eq.1).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
               if ( diagonal ) then
                  n_if2 = i_if1
               else
!!$                  n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                       N_independent_fcts
                  n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2)% &
                       N_independent_fcts
               endif
               do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
                  if ( diagonal .and. i_if1 .eq. i_if2 ) then
                     n_floats = n_floats + i_c1
                  else
                     n_floats = n_floats + n_exp2
                  endif
               enddo
!!$                  enddo
               n_floats_if_c_uc(i_if_c,i_ir) = n_floats
               i_if_c = i_if_c + 1
            enddo
!!$            enddo
         enddo
         n_if_c_uc(i_ir) = i_if_c - 1
         if ( n_if_c_uc(i_ir) .gt. 0 ) then
            n_floats_if_c_uc_max(i_ir) = maxval(n_floats_if_c_uc(1:n_if_c_uc(i_ir),i_ir))
         else
            n_floats_if_c_uc_max(i_ir) = 0
         endif
      enddo
    end subroutine calc_dimensions


    subroutine write_dimensions(th,i_ir)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(readwriteblocked_tapehandle), intent(inout) :: th
      integer(kind=i4_kind), intent(in) :: i_ir
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      real(kind=r8_kind) :: real_n_if_c(1)
      real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
      !------------ Executable code ----------------------------
      !contracted
      allocate( real_n_floats_if_c(n_if_c(i_ir)), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(n_if_c(i_ir),r8_kind)
      real_n_floats_if_c = real(n_floats_if_c(1:n_if_c(i_ir),i_ir),r8_kind)
      call readwriteblocked_write(real_n_if_c,th)
      call readwriteblocked_write(real_n_floats_if_c,th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at write_dimensions() failed")
      !uncontracted
      allocate( real_n_floats_if_c(n_if_c_uc(i_ir)), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(n_if_c_uc(i_ir),r8_kind)
      real_n_floats_if_c = real(n_floats_if_c_uc(1:n_if_c_uc(i_ir),i_ir),r8_kind)
      call readwriteblocked_write(real_n_if_c,th)
      call readwriteblocked_write(real_n_floats_if_c,th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at write_dimensions() failed")
    end subroutine write_dimensions


    subroutine free_dimensions()
      !** End of interface *************************************
      implicit none
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      !------------ Executable code ----------------------------
      deallocate( n_floats_if_c, n_if_c, n_floats_if_c_max, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at free_dimensions() failed")
      deallocate( n_floats_if_c_uc, n_if_c_uc, n_floats_if_c_uc_max, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: deallocating at free_dimensions() failed")
    end subroutine free_dimensions


    subroutine pack_dimensions()
      !** End of interface *************************************
      implicit none
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: info, i_ir
      !------------ Executable code ----------------------------
      !contracted
      call commpack(n_if_c_max,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      call commpack(n_if_c,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      do i_ir = 1, n_ir
         call commpack(n_floats_if_c(:,i_ir),n_if_c_max,1,info)
         if (info .ne. 0) call error_handler( &
              "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      enddo
      call commpack(n_floats_if_c_max,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      !uncontracted
      call commpack(n_if_c_uc_max,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      call commpack(n_if_c_uc,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      do i_ir = 1, n_ir
         call commpack(n_floats_if_c_uc(:,i_ir),n_if_c_uc_max,1,info)
         if (info .ne. 0) call error_handler( &
              "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
      enddo
      call commpack(n_floats_if_c_uc_max,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_send: packing at pack_dimensions()  failed")
    end subroutine pack_dimensions


!!$    subroutine write_integralfile_3c(sa_int_3c,borders)
!!$      implicit none
!!$      !------------ Declaration of formal parameters -----------
!!$      !type(symadapt_totsym_3c_int_type), pointer :: sa_int_3c(:,:)
!!$      ! (spin1,spin2)
!!$      type(symadapt_totsym_3c_int_type), intent(in) :: sa_int_3c(:,:)
!!$      ! (spin1,spin2)
!!$      integer(kind=i4_kind), intent(in) :: borders(3)
!!$      !** End of interface *************************************
!!$      !------------ Declaration of local variables -------------
!!$      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,i_c2,n_if1,n_if2, &
!!$           nn_c2,n_floats,i_buf_lower,i_buf_upper,status, &
!!$           i_spin1,i_spin2,nn_spin2
!!$      real(kind=r8_kind), allocatable :: buffer(:)
!!$      !real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
!!$      !------------ Executable code ----------------------------
!!$      n_floats = borders(3) * n_floats_if_c_max(i_ir)
!!$      allocate( buffer(n_floats),stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: allocating at writing 3 center failed")
!!$
!!$      ! real part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_c1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_3c(i_spin1,i_spin2)%int
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_c2
!!$                     endif
!!$                     do i_c2 = 1, nn_c2
!!$                        i_buf_upper = i_buf_upper + borders(3)
!!$                        !buffer(i_buf_lower:i_buf_upper) = &
!!$                        !     integral(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                        buffer(i_buf_lower:i_buf_upper) = &
!!$                             sa_int_3c(i_spin1,i_spin2)%int(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                        i_buf_lower = i_buf_lower + borders(3)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$               call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      ! immaginary part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_c1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_3c(i_spin1,i_spin2)%int_imag
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_c2
!!$                     endif
!!$                     do i_c2 = 1, nn_c2
!!$                        i_buf_upper = i_buf_upper + borders(3)
!!$                        buffer(i_buf_lower:i_buf_upper) = &
!!$                             sa_int_3c(i_spin1,i_spin2)%int_imag(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                        i_buf_lower = i_buf_lower + borders(3)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$               call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      deallocate(buffer,stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: deallocating at writing 3 center failed")
!!$    end subroutine write_integralfile_3c

!!$
!!$    subroutine write_integralfile_2c(sa_int_2c,contracted)
!!$      implicit none
!!$      !------------ Declaration of formal parameters -----------
!!$      !type(symadapt_totsym_2c_int_type), pointer :: sa_int_2c(:,:)
!!$      type(symadapt_totsym_2c_int_type), intent(in) :: sa_int_2c(:,:)
!!$      ! (spin1,spin2)
!!$      logical, intent(in) :: contracted
!!$      !** End of interface *************************************
!!$      !------------ Declaration of local variables -------------
!!$      real(kind=r8_kind), allocatable :: buffer(:)
!!$      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,n_if1,n_if2, &
!!$           nn_c2,n_floats,i_buf_lower,i_buf_upper,status, &
!!$           n_dim1,n_dim2,i_spin1,i_spin2,nn_spin2
!!$      real(kind=r8_kind), pointer :: integral(:,:,:,:)
!!$      !------------ Executable code ----------------------------
!!$      if (contracted) then
!!$         n_floats = n_floats_if_c_max(i_ir)
!!$         n_dim1=n_c1
!!$         n_dim2=n_c2
!!$      else
!!$         n_floats = n_floats_if_c_uc_max(i_ir)
!!$         n_dim1=n_exp1
!!$         n_dim2=n_exp2
!!$      end if
!!$
!!$      allocate( buffer(n_floats),stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: allocating at writing 2 center failed")
!!$
!!$      ! real part
!!$      do i_spin1 = 1, 2
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_dim1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_2c(i_spin1,i_spin2)%int
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_dim2
!!$                     endif
!!$                     i_buf_upper = i_buf_upper + nn_c2
!!$                     !buffer(i_buf_lower:i_buf_upper) = &
!!$                     !     integral(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     buffer(i_buf_lower:i_buf_upper) = &
!!$                          sa_int_2c(i_spin1,i_spin2)%int(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     i_buf_lower = i_buf_lower + nn_c2
!!$                  enddo
!!$               enddo
!!$               call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      ! immaginary part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_dim1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_2c(i_spin1,i_spin2)%int_imag
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_dim2
!!$                     endif
!!$                     i_buf_upper = i_buf_upper + nn_c2
!!$                     !buffer(i_buf_lower:i_buf_upper) = &
!!$                     !     integral(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     buffer(i_buf_lower:i_buf_upper) = &
!!$                          sa_int_2c(i_spin1,i_spin2)%int_imag(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     i_buf_lower = i_buf_lower + nn_c2
!!$                  enddo
!!$               enddo
!!$               call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      deallocate(buffer,stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: deallocating at writing 2 center failed")
!!$
!!$    end subroutine write_integralfile_2c


!!$    subroutine pack_integralfile_3c(sa_int_3c,borders)
!!$      implicit none
!!$      !------------ Declaration of formal parameters -----------
!!$      !type(symadapt_totsym_3c_int_type), pointer :: sa_int_3c(:,:)
!!$      ! (spin1,spin2)
!!$      type(symadapt_totsym_3c_int_type), intent(in) :: sa_int_3c(:,:)
!!$      ! (spin1,spin2)
!!$      integer(kind=i4_kind), intent(in) :: borders(3)
!!$      !** End of interface *************************************
!!$      !------------ Declaration of local variables -------------
!!$      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,i_c2,n_if1,n_if2, &
!!$           nn_c2,n_floats,i_buf_lower,i_buf_upper,status,info, &
!!$           i_spin1,i_spin2,nn_spin2
!!$      real(kind=r8_kind), allocatable :: buffer(:)
!!$      real(kind=r8_kind) :: checksum
!!$      real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
!!$      !------------ Executable code ----------------------------
!!$      n_floats = borders(3) * n_floats_if_c_max(i_ir)
!!$      allocate( buffer(n_floats),stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: allocating at packing 3 center failed")
!!$
!!$      ! real part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_c1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_3c(i_spin1,i_spin2)%int
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_c2
!!$                     endif
!!$                     do i_c2 = 1, nn_c2
!!$                        i_buf_upper = i_buf_upper + borders(3)
!!$                        buffer(i_buf_lower:i_buf_upper) = &
!!$                             sa_int_3c(i_spin1,i_spin2)%int(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                        i_buf_lower = i_buf_lower + borders(3)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$               call commpack(buffer,i_buf_upper,1,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 3 center failed")
!!$               checksum = sum(buffer(1:i_buf_upper))
!!$               call commpack(checksum,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 3 center checksum failed")
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      ! immaginary part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_c1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                   end if
!!$                  !integral => sa_int_3c(i_spin1,i_spin2)%int_imag
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_c2
!!$                     endif
!!$                     do i_c2 = 1, nn_c2
!!$                        i_buf_upper = i_buf_upper + borders(3)
!!$                        buffer(i_buf_lower:i_buf_upper) = &
!!$                             sa_int_3c(i_spin1,i_spin2)%int_imag(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                        i_buf_lower = i_buf_lower + borders(3)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$               call commpack(buffer,i_buf_upper,1,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 3 center failed")
!!$               checksum = sum(buffer(1:i_buf_upper))
!!$               call commpack(checksum,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 3 center checksum failed")
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      deallocate(buffer,stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: deallocating at packing 3 center failed")
!!$    end subroutine pack_integralfile_3c
!!$
!!$
!!$    subroutine pack_integralfile_2c(sa_int_2c,contracted)
!!$      implicit none
!!$      !------------ Declaration of formal parameters -----------
!!$      !type(symadapt_totsym_2c_int_type), pointer :: sa_int_2c(:,:)
!!$      ! (spin1,spin2)
!!$      type(symadapt_totsym_2c_int_type), intent(in) :: sa_int_2c(:,:)
!!$      ! (spin1,spin2)
!!$      logical, intent(in) :: contracted
!!$      !** End of interface *************************************
!!$      !------------ Declaration of local variables -------------
!!$      real(kind=r8_kind), allocatable :: buffer(:)
!!$      integer(kind=i4_kind) :: i_if1,i_if2,i_c1,info,n_if1,n_if2, &
!!$           nn_c2,n_floats,i_buf_lower,i_buf_upper,status, &
!!$           n_dim1,n_dim2,i_spin1,i_spin2,nn_spin2
!!$      real(kind=r8_kind) :: checksum
!!$      real(kind=r8_kind), pointer :: integral(:,:,:,:)
!!$      !------------ Executable code ----------------------------
!!$      if (contracted) then
!!$         n_floats = n_floats_if_c_max(i_ir)
!!$         n_dim1=n_c1
!!$         n_dim2=n_c2
!!$      else
!!$         n_floats = n_floats_if_c_uc_max(i_ir)
!!$         n_dim1=n_exp1
!!$         n_dim2=n_exp2
!!$      endif
!!$
!!$      allocate(buffer(n_floats),stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: allocating at packing 2 center failed")
!!$
!!$      ! real part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_dim1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_2c(i_spin1,i_spin2)%int
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_dim2
!!$                     endif
!!$                     i_buf_upper = i_buf_upper + nn_c2
!!$                     buffer(i_buf_lower:i_buf_upper) = &
!!$                          sa_int_2c(i_spin1,i_spin2)%int(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     i_buf_lower = i_buf_lower + nn_c2
!!$                  enddo
!!$               enddo
!!$               call commpack(buffer,i_buf_upper,1,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 2 center failed")
!!$               checksum = sum(buffer(1:i_buf_upper))
!!$               call commpack(checksum,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 2 center checksum failed")
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      ! immaginary part
!!$      do i_spin1 = 1, 2
!!$         if ((i_spin1.eq.1).and.(n_spin1.eq.1)) then
!!$            cycle
!!$         end if
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         if (diagonal) then
!!$            nn_spin2 = i_spin1
!!$         else
!!$            nn_spin2 = 2
!!$         endif
!!$         do i_if1 = 1, n_if1
!!$            do i_c1 = 1, n_dim1
!!$               i_buf_lower = 1
!!$               i_buf_upper = 0
!!$               do i_spin2 = 1, nn_spin2
!!$                  if ((i_spin2.eq.1).and.(n_spin2.eq.1)) then
!!$                     cycle
!!$                  end if
!!$                  !integral => sa_int_2c(i_spin1,i_spin2)%int_imag
!!$                  if ( diagonal .and. i_spin2 .eq. i_spin1 ) then
!!$                     n_if2 = i_if1
!!$                  else
!!$                     n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)% &
!!$                          N_independent_fcts
!!$                  endif
!!$                  do i_if2 = 1, n_if2
!!$                     if ( diagonal .and. i_spin2 .eq. i_spin1 .and. i_if1 .eq. i_if2 ) then
!!$                        nn_c2 = i_c1
!!$                     else
!!$                        nn_c2 = n_dim2
!!$                     endif
!!$                     i_buf_upper = i_buf_upper + nn_c2
!!$                     buffer(i_buf_lower:i_buf_upper) = &
!!$                          sa_int_2c(i_spin1,i_spin2)%int_imag(1:nn_c2,i_c1,i_if2,i_if1)
!!$                     i_buf_lower = i_buf_lower + nn_c2
!!$                  enddo
!!$               enddo
!!$               call commpack(buffer,i_buf_upper,1,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 2 center failed")
!!$               checksum = sum(buffer(1:i_buf_upper))
!!$               call commpack(checksum,info)
!!$               if (info .ne. 0) call error_handler( &
!!$                    "int_send_2cob3c_spor_send: packing 2 center checksum failed")
!!$            enddo
!!$         enddo
!!$      enddo
!!$
!!$      deallocate(buffer,stat=status)
!!$      if (status .ne. 0) call error_handler( &
!!$           "int_send_2cob3c_spor_send: deallocating at packing 2 center failed")
!!$    end subroutine pack_integralfile_2c

  end subroutine int_send_2cob3c_spor_send_file
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cob3c_spor_send_mem
    !  Purpose: sending of contracted and symmetry adapted
    !    integrals from int_data_2cob3c_module.
    !    When packing the integrals, the mapping to the scf-
    !    metaindex is also done.
    !    Two center integrals are send to the master and
    !    three center integrals are distributed over all hosts
    !    splitting them over the fitfunction metaindex.
    !    integrals to be hold on present hosts are directly
    !    stored in memory in integralstore_module
    !  called by: integral_calc_quad_2cob3c
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    implicit none

    ABORT('FIXME: ints in mem')
#ifdef SO_INTEGRALS_IN_MEM
    n_ir=symmetry_data_n_proj_irreps()
    do i_host = first_host, last_host
       if ( i_host .eq. my_hostindex ) then
          call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))
          if ( i_host .eq. 1 ) then
             if ( integralpar_2cob_kin ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing kin")
                if(.not.integralpar_relativistic) then
                   call store_integral_2c(symadapt_int_2cob_kin_p,integralstore_2cob_kin)
                else
                   call store_integral_2c(symadapt_int_2cob_kin_p,integralstore_2cob_kin_rel)
                end if
             endif
             if ( integralpar_2cob_nuc ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing nuc")
                if(.not.integralpar_relativistic) then
                   call store_integral_2c(symadapt_int_2cob_nuc_p,integralstore_2cob_nuc)
                else
                   call store_integral_2c(symadapt_int_2cob_nuc_p,integralstore_2cob_nuc_rel)
                end if
             endif
             if ( integralpar_2cob_pvsp ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing nuc")
                call store_integral_2c(symadapt_int_2cob_pvsp_p,integralstore_2cob_pvsp)
             endif
             if ( integralpar_2cob_pvxp ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing pvxp")
                call store_integral_2c(sa_2c_ints(:, int_sa_pvxp),integralstore_2cob_pvxp)
             endif

             if ( integralpar_2cob_pvec ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing sigp")
                call store_integral_2c(symadapt_int_2cob_sigp_p,integralstore_2cob_sigp)
             endif
          endif! if ( i_host .eq. 1 )
             if ( integralpar_2cob_ol ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing overlap")
                if(.not.integralpar_relativistic) then
                   call store_integral_2c(symadapt_int_2cob_ol_p,integralstore_2cob_ol)
                else
                   call store_integral_2c(symadapt_int_2cob_ol_p,integralstore_2cob_ol_rel)
                end if
             endif
             if ( integralpar_3c_xc ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing coulomb")
                call store_integral_3c( symadapt_int_3c_xc_p, &
                     integralstore_3c_xc, borders_xc(:,my_hostindex) )
             endif
             if ( integralpar_3c_co ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: writing exchange")
                call store_integral_3c( symadapt_int_3c_co_p, &
                     integralstore_3c_co, borders_ch(:,my_hostindex) )
             endif
          call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
          n_missing_quadrupels = n_missing_quadrupels - 1
       else ! packing
          call start_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
          call comm_init_send(i_host,msgtag_int_2cob3c_result)
          call quadrupel_pack(quadrupel)
          call commpack(n_c1,info)
          if (info .ne. 0) call error_handler( &
               "int_send_2cob3c_send: packing n_c1 failed")
          call commpack(n_c2,info)
          if (info .ne. 0) call error_handler( &
               "int_send_2cob3c_send: packing n_c2 failed")
          if(integralpar_relativistic) then
             call commpack(n_exp1,info)
             if (info .ne. 0) call error_handler( &
                  "int_send_2cob3c_send: packing n_exp1 failed")
             call commpack(n_exp2,info)
             if (info .ne. 0) call error_handler( &
                  "int_send_2cob3c_send: packing n_exp2 failed")
          end if
          if ( i_host .eq. 1 ) then
             if ( integralpar_2cob_kin ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: packing kin")
                call pack_integralfile_2c(symadapt_int_2cob_kin_p)
             endif
             if ( integralpar_2cob_nuc ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: packing nuc")
                call pack_integralfile_2c(symadapt_int_2cob_nuc_p)
             endif
             if ( integralpar_2cob_pvsp ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: packing pvsp")
                call pack_integralfile_2c(symadapt_int_2cob_pvsp_p)
             endif
             if ( integralpar_2cob_pvxp ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: packing pvxp")
                call pack_integralfile_2c(symadapt_int_2cob_pvxp_p)
             endif

             if ( integralpar_2cob_pvec ) then
                if ( output_int_loops ) call write_to_output_units( &
                     "int_send_2cob3c_send: packing sigp")
                call pack_integralfile_2c(symadapt_int_2cob_sigp_p)
             endif
          endif
          if ( integralpar_2cob_ol ) then
             if ( output_int_loops ) call write_to_output_units( &
                  "int_send_2cob3c_send: packing overlap")
             call pack_integralfile_2c(symadapt_int_2cob_ol_p)
          endif
          if ( integralpar_3c_xc ) then
             if ( output_int_loops ) call write_to_output_units( &
                  "int_send_2cob3c_send: packing coulomb")
             call pack_integralfile_3c( symadapt_int_3c_xc_p, &
                  borders_xc(:,i_host) )
          endif
          if ( integralpar_3c_co ) then
             if ( output_int_loops ) call write_to_output_units( &
                  "int_send_2cob3c_send: packing exchange")
             call pack_integralfile_3c( symadapt_int_3c_co_p, &
                  borders_ch(:,i_host) )
          endif
          call stop_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
          call start_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
          call comm_send()
          call stop_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
       end if! i_host.eq.my_hostindex
    enddo! hosts
    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cob3_send: done")

  contains


    subroutine store_integral_3c(symadapt_int,integralstore,border)
      !------------ Declaration of formal parameters -----------
      type(symadapt_totsym_3c_int_type), pointer :: symadapt_int(:)
!!$      real(kind=r8_kind), intent(out) :: integralstore(:)
!!$      integer(kind=i4_kind), intent(in) :: border(3)
      real(kind=r8_kind) :: integralstore(:) ! avoid warnings
      integer(kind=i4_kind) :: border(3)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer :: integral_real(:,:,:,:,:),integral_imag(:,:,:,:,:)
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_c2,i_meta, &
           n_if1,n_if2,nn_c2,i_spin1,i_spin2
      !------------ Executable code ----------------------------
      call error_handler("store_integral_3c: need an update")
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$            do i_if1 = 1, n_if1
!!$               do i_c1 = 1, n_c1
!!$                  do i_spin2 = 1,2
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     integral_real => symadapt_int(i_ir,i_spin1,i_spin2)%int
!!$                     integral_imag => symadapt_int(i_ir,i_spin1,i_spin2)%int_imag
!!$                     i_meta = (metaindex(quadrupel,i_ir,i_spin1,i_spin2,i_if1,i_c1) - 1) * border(3) + 1
!!$                     if (diagonal.and.(i_spin1.eq.i_spin2)) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                     endif
!!$                     do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and.(i_spin1.eq.i_spin2) ) then
!!$                           nn_c2 = i_c1
!!$                        else
!!$                           nn_c2 = n_c2
!!$                        endif
!!$                        do i_c2 = 1, nn_c2
!!$                           integralstore(i_meta:i_meta+2*border(3)-2:2) = &
!!$                                integral_real(i_c2,i_c1,border(1):border(2),i_if2,i_if1)
!!$                           integralstore(i_meta+1:i_meta+2*border(3)-1:2) = &
!!$                                integral_imag(i_c2,i_c1,border(1):border(2),i_if2,i_if1)
!!$                           i_meta = i_meta + 2*border(3)
!!$                        enddo! i_c2
!!$                     enddo! i_if2
!!$                  enddo ! i_spin2
!!$               enddo! i_c1
!!$            enddo! i_if1
!!$         enddo ! i_spin1
!!$      enddo ! i_ir
    end subroutine store_integral_3c


    subroutine store_integral_2c(symadapt_int,integralstore)
      !------------ Declaration of formal parameters -----------
      type(symadapt_totsym_2c_int_type), pointer :: symadapt_int(:)
!!$      real(kind=r8_kind), intent(out) :: integralstore(:)
      real(kind=r8_kind) :: integralstore(:) ! avoid warnings
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer :: integral_real(:,:,:,:),integral_imag(:,:,:,:)
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_c2,i_meta, &
           n_if1,n_if2,nn_c2,n_dim1,n_dim2,counter !!$,i_spin1,i_spin2
      logical :: rel
      !------------ Executable code ----------------------------
      call error_handler("store_integral_2c: need an update")
!!$      if(.not.integralpar_relativistic) then
!!$         rel=.false.
!!$         n_dim1=n_c1
!!$         n_dim2=n_c2
!!$      else
!!$         rel=.true.
!!$         n_dim1=n_exp1
!!$         n_dim2=n_exp2
!!$      end if
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1)%N_independent_fcts
!!$            do i_if1 = 1, n_if1
!!$               do i_c1 = 1, n_dim1
!!$                  do i_spin2 = 1,2
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     integral_real => symadapt_int(i_ir,i_spin1,i_spin2)%int
!!$                     integral_imag => symadapt_int(i_ir,i_spin1,i_spin2)%int_imag
!!$                     i_meta = metaindex(quadrupel,i_ir,i_spin1,i_spin2,i_if1,i_c1,rel)
!!$                  integral_real => symadapt_int(i_ir)%int
!!$                  integral_imag => symadapt_int(i_ir)%int_imag
!!$                  i_meta = metaindex(quadrupel,i_ir,i_if1,i_c1,rel)
!!$                     counter = i_meta
!!$                     if (diagonal.and.(i_spin1.eq.i_spin2)) then
!!$                      if (diagonal) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2)%N_independent_fcts
!!$                     endif
!!$                     do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and.(i_spin1.eq.i_spin2) ) then
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2) ) then
!!$                           nn_c2 = i_c1
!!$                        else
!!$                           nn_c2 = n_dim2
!!$                        endif
!!$                        do i_c2 = 1,nn_c2
!!$                           ! real part
!!$                           integralstore(counter) = &
!!$                                integral_real(i_c2,i_c1,i_if2,i_if1)
!!$                           counter = counter + 1
!!$                           ! imaginary part
!!$                           integralstore(counter) = &
!!$                                integral_imag(i_c2,i_c1,i_if2,i_if1)
!!$                           counter = counter + 1
!!$                        enddo ! i_c2
!!$                     enddo! i_if2
!!$                  enddo ! i_spin2
!!$               enddo ! i_c1
!!$            enddo! i_if1
!!$         enddo! i_spin1
!!$      enddo ! i_ir
!!$      print*,"store_integral_2c: LEFT"
      !? print*,"store_integral_2c: LEFT"
    end subroutine store_integral_2c

   subroutine pack_integralfile_2c(symadapt_int)
      !------------ Declaration of formal parameters -----------
     type(symadapt_totsym_2c_int_type), pointer :: symadapt_int(:)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer :: integral_real(:,:,:,:),integral_imag(:,:,:,:)
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
           n_if1,n_if2,nn_c2,i_buf_lower,i_buf_upper,status,&
           n_dim1,n_dim2,i_c2,counter !!$,i_spin1,i_spin2
      real(kind=r8_kind), allocatable :: buffer(:)
      integer(kind=i4_kind)           :: length,length_new
      ! for debugging >>>
      integer(kind=i4_kind)  :: number_packed
      ! for debugging <<<
       !------------ Executable code ----------------------------
      call error_handler("pack_integralfile_2c: need an update")
!!$      if(.not.integralpar_relativistic) then
!!$         n_dim1=n_c1
!!$         n_dim2=n_c2
!!$      else
!!$         n_dim1=n_exp1
!!$         n_dim2=n_exp2
!!$      end if
!!$      number_packed = 0
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         length = 0
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            do i_spin2 = 1,2
!!$               if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                  cycle
!!$               endif
!!$               if (associated(symadapt_int(i_ir,i_spin1,i_spin2)%int)) then
!!$                if (associated(symadapt_int(i_ir)%int)) then
!!$                  length_new = size(symadapt_int(i_ir,i_spin1,i_spin2)%int)
!!$                  length_new = size(symadapt_int(i_ir)%int)
!!$                  if (length_new.gt.length) then
!!$                     length = length_new
!!$                  endif
!!$               endif
!!$            enddo
!!$         enddo
!!$         allocate( buffer(24*length) ,stat=status)
!!$         ! factor 4 since: one factor two from complex and another factor two
!!$         ! from two spins
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_send: allocating at packing 2 center failed")
!!$         !i_buf_lower = 1
!!$         !i_buf_upper = 0
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            counter = 0
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1)%N_independent_fcts
!!$            print*,">> n_if1: ",n_if1
!!$            call buggy(">> n_if1: ",inte=n_if1)
!!$            do i_if1 = 1, n_if1
!!$               print*," i_if1",i_if1
!!$               do i_c1 = 1, n_dim1
!!$                  print*," i_c1",i_c1
!!$                  do i_spin2 = 1,2
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     integral_real => symadapt_int(i_ir,i_spin1,i_spin2)%int
!!$                     integral_imag => symadapt_int(i_ir,i_spin1,i_spin2)%int_imag
!!$                     if (diagonal.and.(i_spin1.eq.i_spin2)) then
!!$                     integral_real => symadapt_int(i_ir)%int
!!$                     integral_imag => symadapt_int(i_ir)%int_imag
!!$                     if (diagonal) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2)%N_independent_fcts
!!$                     endif
!!$                     do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and.(i_spin1.eq.i_spin2) ) then
!!$                           nn_c2 = i_c1
!!$                        else
!!$                           nn_c2 = n_dim2
!!$                        endif
!!$                        do i_c2 = 1,nn_c2
!!$                           !print*,"   i_c2",i_c2
!!$                           ! real part
!!$                           counter = counter + 1
!!$                           buffer(counter) = &
!!$                                integral_real(i_c2,i_c1,i_if2,i_if1)
!!$                           !print*,"   counter ",counter
!!$                           ! imaginary part
!!$                           counter = counter + 1
!!$                           buffer(counter) = &
!!$                                integral_imag(i_c2,i_c1,i_if2,i_if1)
!!$                           !print*,"   counter ",counter
!!$                        enddo ! i_c2
!!$                     enddo ! i_if2
!!$                  enddo ! i_spin2
!!$               enddo ! i_c1
!!$            enddo ! i_if1
!!$            call commpack(buffer(1:counter),counter,1,info)
!!$            number_packed = number_packed + counter
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_send: packing 2 center failed")
!!$         enddo! i_spin1
!!$         deallocate(buffer,stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_send: deallocating at packing 2 center failed")
!!$      enddo! i_ir
    end subroutine pack_integralfile_2c


    subroutine pack_integralfile_3c(symadapt_int,borders)
      !------------ Declaration of formal parameters -----------
      type(symadapt_totsym_3c_int_type), pointer :: symadapt_int(:)
      integer(kind=i4_kind), intent(in) :: borders(3)
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      real(kind=r8_kind), pointer :: integral_real(:,:,:,:,:),integral_imag(:,:,:,:,:)
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_c2,info, &
           n_if1,n_if2,nn_c2,i_buf_lower,i_buf_upper,status !!$,i_spin1,i_spin2
      real(kind=r8_kind), allocatable :: buffer(:)
      integer(kind=i4_kind)           :: length,length_new
      !------------ Executable code ----------------------------

      call error_handler("pack_integralfile_3c: need modification")
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         length = 0
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            do i_spin2 = 1,2
!!$               if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                  cycle
!!$               endif
!!$               if (associated(symadapt_int(i_ir,i_spin1,i_spin2)%int)) then
!!$                  length_new = size(symadapt_int(i_ir,i_spin1,i_spin2)%int)
!!$                  if (length_new.gt.length) then
!!$                     length = length_new
!!$                  endif
!!$               endif
!!$            enddo
!!$         enddo
!!$         if (length.eq.0) then
!!$            cycle
!!$         endif
!!$         allocate (buffer(4*borders(3)*length),stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_send: allocating of buffer at packing 3 center failed")
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$            i_buf_lower = 1
!!$            i_buf_upper = 0
!!$            do i_if1 = 1, n_if1
!!$               if ( diagonal ) n_if2 = i_if1
!!$               do i_c1 = 1, n_c1
!!$                  do i_spin2 = 1,2
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     integral_real => symadapt_int(i_ir,i_spin1,i_spin2)%int
!!$                     integral_imag => symadapt_int(i_ir,i_spin1,i_spin2)%int_imag
!!$                     if (diagonal.and.(i_spin1.eq.i_spin2)) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                     endif
!!$                     do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and.(i_spin1.eq.i_spin2) ) then
!!$                           nn_c2 = i_c1
!!$                        else
!!$                           nn_c2 = n_c2
!!$                        endif
!!$                        do i_c2 = 1, nn_c2
!!$                           ! real part
!!$                           i_buf_upper = i_buf_upper + 2*borders(3)
!!$                           buffer(i_buf_lower:i_buf_upper-1:2) = &
!!$                                integral_real(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                           ! imaginary part
!!$                           buffer(i_buf_lower+1:i_buf_upper:2) = &
!!$                                integral_imag(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
!!$                           i_buf_lower = i_buf_lower + 2*borders(3)
!!$                        enddo ! i_c2
!!$                     enddo ! i_if2
!!$                  enddo ! i_spin2
!!$               enddo ! i_c1
!!$            enddo ! i_if1
!!$            call commpack(buffer,i_buf_upper,1,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_send: packing 3 center failed")
!!$         enddo ! i_spin1
!!$         deallocate(buffer,stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_send: deallocating at packing 3 center failed")
!!$      enddo ! i_ir
    end subroutine pack_integralfile_3c
#endif
  end subroutine int_send_2cob3c_spor_send_mem
  !*************************************************************



  !*************************************************************
  subroutine int_send_2cob3c_spor_receive
    !  Purpose:
    !    + receiving of contracted and symmetry adapted
    !      integrals for one quadrupel
    !    + Storing them either in memory or on file
    !  called by: main_slave, integral_main_2cob3c,
    !    integral_interrupt_2cob3c, int_send_2cob3c_spor_shutdown
    !** End of interface ***************************************
    if ( options_integrals_on_file() ) then
       call int_send_2cob3c_sp_receive_file()
    else
       call int_send_2cob3c_sp_receive_mem()
    endif
  end subroutine int_send_2cob3c_spor_receive
  !*************************************************************



  !*************************************************************
  subroutine int_send_2cob3c_sp_receive_file
    !  Purpose:
    !    + receiving of contracted and symmetry adapted
    !      integrals for one quadrupel
    !    + Storing them in direct access files in scf symmetric
    !      storage mode with the record number as metindex
    !  called by: main_slave, integral_main_2cob3c,
    !    integral_interrupt_2cob3c, int_send_2cob3c_spor_shutdown
    !** End of interface *****************************************
    use symm_adapt_int, only: UnPackIntBlock
    implicit none
    !------------ Declaration of local variables ---------------------
    type(readwriteblocked_tapehandle) :: th
    type(quadrupel_type) :: quadrupel
    integer(kind=i4_kind) :: i_ir, n_ir, n_if_c_max, n_if_c_uc_max
    logical :: diagonal, restart_timer
    integer(kind=i4_kind), allocatable :: &
         n_floats_if_c(:,:), n_floats_if_c_max(:), n_if_c(:), &
         n_floats_if_c_uc(:,:), n_floats_if_c_uc_max(:), n_if_c_uc(:)
    ! n_floats_if_c(i_if_c,i_ir), n_floats_if_c_max(i_ir), n_if_c(i_ir)
    ! n_floats_if_c_uc(i_if_c,i_ir), n_floats_if_c_uc_max(i_ir), n_if_c_uc(i_ir)
    !------------ Executable code ------------------------------------

    if ( timer_int_idle_2cob3c(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))

    call quadrupel_unpack(quadrupel)

    if ( output_int_loops .and. output_unit > 0 ) then
       write(output_unit,*) &
            "int_send_2cob3c_spor_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
       write(stdout_unit,*) &
            "int_send_2cob3c_spor_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
    endif

    diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
         (quadrupel%l1 == quadrupel%l2)

    n_ir = symmetry_data_n_proj_irreps()

    call unpack_dimensions()

    do i_ir = 1, n_ir
       call readwriteblocked_startwrite( &
            trim(tmpfile(quadrupel_filename(i_ir,quadrupel))), &
            th, blocklength=blocklength )
       call write_dimensions(th,i_ir)
       if ( comm_i_am_master()) then
          if ( integralpar_2cob_kin ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking kin")
             call UnPackIntBlock(th)
          endif
          if ( integralpar_2cob_nuc ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking nuc")
             call UnPackIntBlock(th)
          endif
          if ( integralpar_2cob_pvsp ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking pvsp")
             call UnPackIntBlock(th)
          endif
          if ( integralpar_2cob_pvxp ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking pvxp")
             call UnPackIntBlock(th)
          endif
          if ( integralpar_2cob_pvec ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking sigp")
             call UnPackIntBlock(th)
          endif
       endif
       if ( integralpar_2cob_ol ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_2cob3c_spor_receive: unpacking overlap")
          call UnPackIntBlock(th)
       endif
       if ( integralpar_3c_xc ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_2cob3c_spor_receive: unpacking exchange")
          call UnPackIntBlock(th)
       endif
       if ( integralpar_3c_co ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_2cob3c_spor_receive: unpacking coulomb")
          call UnPackIntBlock(th)
       endif
       if ( integralpar_3c_rcoul_pvsp ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_2cob3c_spor_receive: unpacking rcoul_pvsp")
          call UnPackIntBlock(th)
       endif
       if ( integralpar_3c_rcoul_pvxp ) then
          if ( output_int_deeploops ) call write_to_output_units( &
               "int_send_2cob3c_spor_receive: unpacking rcoul_pvxp")
          ! UnPackIntBlock : calling the integrals and storing in 'th'
          call UnPackIntBlock(th)
       endif
          if ( integralpar_3c_r2_pvsp ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking r2_pvsp")
             call UnPackIntBlock(th)
          endif
          if ( integralpar_3c_r2_pvxp ) then
             if ( output_int_deeploops ) call write_to_output_units( &
                  "int_send_2cob3c_spor_receive: unpacking r2_pvxp")
             call UnPackIntBlock(th)
          endif
       call readwriteblocked_stopwrite(th, total_length= &
            quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
            index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ) )
    enddo! loop over irreps

    n_missing_quadrupels = n_missing_quadrupels - 1

    call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
    call stop_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
    if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cob3c_spor_receive: done")

  contains


    subroutine write_dimensions(th,i_ir)
      implicit none
      !------------ Declaration of formal parameters -----------
      type(readwriteblocked_tapehandle), intent(inout) :: th
      integer(kind=i4_kind), intent(in) :: i_ir
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      real(kind=r8_kind) :: real_n_if_c(1)
      real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
      !------------ Executable code ----------------------------
      !contracted
      allocate( real_n_floats_if_c(n_if_c(i_ir)), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(n_if_c(i_ir),r8_kind)
      real_n_floats_if_c = real(n_floats_if_c(1:n_if_c(i_ir),i_ir),r8_kind)
      call readwriteblocked_write(real_n_if_c,th)
      call readwriteblocked_write(real_n_floats_if_c,th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: deallocating at write_dimensions() failed")
      !uncontracted
      allocate( real_n_floats_if_c(n_if_c_uc(i_ir)), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: allocating at write_dimensions() failed")
      real_n_if_c(1) = real(n_if_c_uc(i_ir),r8_kind)
      real_n_floats_if_c = real(n_floats_if_c_uc(1:n_if_c_uc(i_ir),i_ir),r8_kind)
      call readwriteblocked_write(real_n_if_c,th)
      call readwriteblocked_write(real_n_floats_if_c,th)
      deallocate( real_n_floats_if_c, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: deallocating at write_dimensions() failed")
    end subroutine write_dimensions


    subroutine free_dimensions()
      !** End of interface *************************************
      implicit none
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: status
      !------------ Executable code ----------------------------
      deallocate( n_floats_if_c, n_if_c, n_floats_if_c_max, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: deallocating at free_dimensions() failed")
      deallocate( n_floats_if_c_uc, n_if_c_uc, n_floats_if_c_uc_max, stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: deallocating at free_dimensions() failed")
    end subroutine free_dimensions


    subroutine unpack_dimensions()
      !** End of interface *************************************
      implicit none
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: info, status, i_ir
      !------------ Executable code ----------------------------
      !contracted
      call communpack(n_if_c_max,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(1)  failed")
      allocate( n_floats_if_c(n_if_c_max,n_ir), &
           n_if_c(n_ir), n_floats_if_c_max(n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: allocating at unpack_dimensions(2) failed")
      call communpack(n_if_c,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(3)  failed")
      do i_ir = 1, n_ir
         call communpack(n_floats_if_c(:,i_ir),n_if_c_max,1,info)
         if (info .ne. 0) call error_handler( &
              "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(4)  failed")
      enddo
      call communpack(n_floats_if_c_max,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(5)  failed")
      !uncontracted
      call communpack(n_if_c_uc_max,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(6)  failed")
      allocate( n_floats_if_c_uc(n_if_c_uc_max,n_ir), &
           n_if_c_uc(n_ir), n_floats_if_c_uc_max(n_ir), stat=status)
      if (status .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: allocating at unpack_dimensions(7) failed")
      call communpack(n_if_c_uc,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(8)  failed")
      do i_ir = 1, n_ir
         call communpack(n_floats_if_c_uc(:,i_ir),n_if_c_uc_max,1,info)
         if (info .ne. 0) call error_handler( &
              "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(9)  failed")
      enddo
      call communpack(n_floats_if_c_uc_max,n_ir,1,info)
      if (info .ne. 0) call error_handler( &
           "int_send_2cob3c_spor_receive: unpacking at unpack_dimensions(10)  failed")
    end subroutine unpack_dimensions


!!$    subroutine unpack_integralfile(bordersize,contracted)
!!$      implicit none
!!$      !------------ Declaration of formal parameters -----------
!!$      integer(kind=i4_kind), intent(in) :: bordersize
!!$      logical :: contracted
!!$      !** End of interface *************************************
!!$      !------------ Declaration of local variables -------------
!!$      integer(kind=i4_kind) :: info,n_floats,status,i_if_c
!!$      real(kind=r8_kind), allocatable :: buffer(:)
!!$      real(kind=r8_kind) :: checksum_received, checksum_calculated
!!$      character(len=200) :: string
!!$      !------------ Executable code ----------------------------
!!$
!!$      if ( contracted ) then
!!$
!!$         n_floats = n_floats_if_c_max(i_ir) * bordersize
!!$         allocate(buffer(n_floats),stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_spor_receive: allocating at unpacking integrals failed")
!!$
!!$         ! real
!!$         do i_if_c = 1, n_if_c(i_ir)
!!$            n_floats = n_floats_if_c(i_if_c,i_ir) * bordersize
!!$            call communpack(buffer,n_floats,1,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking integrals failed")
!!$            checksum_calculated = sum(buffer(1:n_floats))
!!$            call communpack(checksum_received,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking checksum failed")
!!$            if ( checksum_calculated .ne. checksum_received) then
!!$               write(string,*) "int_send_2cob3c_spor_receive: checksum error: received ", &
!!$                    checksum_received, " and calculated ", checksum_calculated
!!$               call error_handler(trim(string))
!!$            endif
!!$            call readwriteblocked_write(buffer(1:n_floats),th)
!!$         enddo
!!$
!!$         ! immaginary
!!$         do i_if_c = 1, n_if_c(i_ir)
!!$            n_floats = n_floats_if_c(i_if_c,i_ir) * bordersize
!!$            call communpack(buffer,n_floats,1,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking integrals failed")
!!$            checksum_calculated = sum(buffer(1:n_floats))
!!$            call communpack(checksum_received,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking checksum failed")
!!$            if ( checksum_calculated .ne. checksum_received) then
!!$               write(string,*) "int_send_2cob3c_spor_receive: checksum error: received ", &
!!$                    checksum_received, " and calculated ", checksum_calculated
!!$               call error_handler(trim(string))
!!$            endif
!!$            call readwriteblocked_write(buffer(1:n_floats),th)
!!$         enddo
!!$
!!$         deallocate(buffer, stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_spor_receive: deallocating at unpacking failed")
!!$
!!$      else ! not contracted
!!$
!!$         n_floats = n_floats_if_c_uc_max(i_ir) * bordersize
!!$         allocate(buffer(n_floats),stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_spor_receive: allocating at unpacking integrals failed")
!!$
!!$         ! real
!!$         do i_if_c = 1, n_if_c_uc(i_ir)
!!$            n_floats = n_floats_if_c_uc(i_if_c,i_ir) * bordersize
!!$            call communpack(buffer,n_floats,1,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking integrals failed")
!!$            checksum_calculated = sum(buffer(1:n_floats))
!!$            call communpack(checksum_received,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking checksum failed")
!!$            if ( checksum_calculated .ne. checksum_received) then
!!$               write(string,*) "int_send_2cob3c_spor_receive: checksum error: received ", &
!!$                    checksum_received, " and calculated ", checksum_calculated
!!$               call error_handler(trim(string))
!!$            endif
!!$            if(.not.integralpar_rel_gradients) then
!!$               call readwriteblocked_write(buffer(1:n_floats),th)
!!$            else
!!$               call readwriteblocked_write(buffer(1:n_floats),th_pointer)
!!$            end if
!!$         enddo
!!$
!!$         ! immaginary
!!$         do i_if_c = 1, n_if_c_uc(i_ir)
!!$            n_floats = n_floats_if_c_uc(i_if_c,i_ir) * bordersize
!!$            call communpack(buffer,n_floats,1,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking integrals failed")
!!$            checksum_calculated = sum(buffer(1:n_floats))
!!$            call communpack(checksum_received,info)
!!$            if (info .ne. 0) call error_handler( &
!!$                 "int_send_2cob3c_spor_receive: unpacking checksum failed")
!!$            if ( checksum_calculated .ne. checksum_received) then
!!$               write(string,*) "int_send_2cob3c_spor_receive: checksum error: received ", &
!!$                    checksum_received, " and calculated ", checksum_calculated
!!$               call error_handler(trim(string))
!!$            endif
!!$            if(.not.integralpar_rel_gradients) then
!!$               call readwriteblocked_write(buffer(1:n_floats),th)
!!$            else
!!$               call readwriteblocked_write(buffer(1:n_floats),th_pointer)
!!$            end if
!!$         enddo
!!$
!!$         deallocate(buffer, stat=status)
!!$         if (status .ne. 0) call error_handler( &
!!$              "int_send_2cob3c_spor_receive: deallocating at unpacking failed")
!!$      end if
!!$    end subroutine unpack_integralfile

  end subroutine int_send_2cob3c_sp_receive_file
  !*************************************************************




  !*************************************************************
  subroutine int_send_2cob3c_sp_receive_mem
    !  Purpose:
    !    + receiving of contracted and symmetry adapted
    !      integrals for one quadrupel
    !    + Storing them in memory in integralstore_module
    !  called by: main_slave, integral_main_2cob3c,
    !    integral_interrupt_2cob3c, int_send_2cob3c_spor_shutdown
    implicit none
    !** End of interface *****************************************

    ABORT('FIXME: ints in mem')
#ifdef SO_INEGRALS_IN_MEM
    type(quadrupel_type)            :: quadrupel
    integer(kind=i4_kind)           :: n_c1, n_c2, n_exp1, n_exp2, info
    type(unique_atom_type), pointer :: ua1,ua2
    logical                         :: diagonal, restart_timer

    if ( timer_int_idle_2cob3c(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))

    call quadrupel_unpack(quadrupel)

    if ( output_int_loops .and. output_unit > 0 ) then
       write(output_unit,*) &
            "int_send_2cob3c_spor_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
       write(stdout_unit,*) &
            "int_send_2cob3c_spor_receive: start with quadrupel ", &
            quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
    endif

    ua1 => unique_atoms(quadrupel%ua1)
    ua2 => unique_atoms(quadrupel%ua2)
    diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
               (quadrupel%l1 == quadrupel%l2)
    call communpack(n_c1,info)
    if (info .ne. 0) call error_handler( &
         "int_send_2cob3c_spor_receive: unpacking n_c1 failed")
    call communpack(n_c2,info)
    if (info .ne. 0) call error_handler( &
         "int_send_2cob3c_spor_receive: unpacking n_c2 failed")

    if(integralpar_relativistic) then
       call communpack(n_exp1,info)
       if (info .ne. 0) call error_handler( &
            "int_send_2cob3c_receive: unpacking n_exp1 failed")
       call communpack(n_exp2,info)
       if (info .ne. 0) call error_handler( &
         "int_send_2cob3c_receive: unpacking n_exp2 failed")
    end if

    if ( comm_i_am_master()) then
       if ( integralpar_2cob_kin ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_receive: unpacking kin")
          if(.not.integralpar_relativistic) then
             call unpack_integralfile_2c(integralstore_2cob_kin)
          else
             call unpack_integralfile_2c(integralstore_2cob_kin_rel)
          end if
       endif
       if ( integralpar_2cob_nuc ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_receive: unpacking nuc")
          if(.not.integralpar_relativistic) then
             call unpack_integralfile_2c(integralstore_2cob_nuc)
          else
             call unpack_integralfile_2c(integralstore_2cob_nuc_rel)
          end if
       endif
       if ( integralpar_2cob_pvsp ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_receive: unpacking pvsp")
          call unpack_integralfile_2c(integralstore_2cob_pvsp)
       endif
       if ( integralpar_2cob_pvxp ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_receive: unpacking pvxp")
          call unpack_integralfile_2c(integralstore_2cob_pvxp)
       endif

       if ( integralpar_2cob_pvec ) then
          if ( output_int_loops ) call write_to_output_units( &
               "int_send_2cob3c_receive: unpacking sigp")
          call unpack_integralfile_2c(integralstore_2cob_sigp)
       endif
    endif
    if ( integralpar_2cob_ol ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_receive: unpacking overlap")
       if(.not.integralpar_relativistic) then
          call unpack_integralfile_2c(integralstore_2cob_ol)
       else
          call unpack_integralfile_2c(integralstore_2cob_ol_rel)
       end if
    endif
    if ( integralpar_3c_xc ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_receive: unpacking exchange")
       call unpack_integralfile_3c(integralstore_3c_xc, &
         borders_xc(3,my_hostindex))
    endif
    if ( integralpar_3c_co ) then
       if ( output_int_loops ) call write_to_output_units( &
            "int_send_2cob3c_spor_receive: unpacking coulomb")
       call unpack_integralfile_3c(integralstore_3c_co, &
         borders_ch(3,my_hostindex))
    endif


    n_missing_quadrupels = n_missing_quadrupels - 1

    call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
    call stop_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
    if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "int_send_2cob3c_spor_receive: done")

  contains


    subroutine unpack_integralfile_2c(integralstore)
      !------------ Declaration of formal parameters -----------
!!$      real(kind=r8_kind), intent(out) :: integralstore(:)
      real(kind=r8_kind) :: integralstore(:) ! avoid warnings
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
           n_if1,n_if2,nn_c2,i_meta,n_dim1,n_dim2 !!$,i_spin1,i_spin2,
      ! for debugging >>>
      integer(kind=i4_kind)  :: number_read
      ! for debugging <<<
      logical :: rel
      !------------ Executable code ----------------------------
      call error_handler("unpack_integralfile_2c: need an update")
!!$      if(.not.integralpar_relativistic) then
!!$         rel=.false.
!!$         n_dim1=n_c1
!!$         n_dim2=n_c2
!!$      else
!!$         rel=.true.
!!$         n_dim1=n_exp1
!!$         n_dim2=n_exp2
!!$      end if
!!$      number_read = 0
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$         n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1)%N_independent_fcts
!!$            do i_if1 = 1, n_if1
!!$               do i_c1 = 1, n_dim1
!!$                  do i_spin2 = 1,2
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!$                      if ( diagonal.and.(i_spin1.eq.i_spin2) ) then
!!$                  if ( diagonal ) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2)%N_independent_fcts
!!$                     endif
!!$                     i_meta = metaindex(quadrupel,i_ir,i_spin1,i_spin2,i_if1,i_c1,rel)
!!$                     i_meta = metaindex(quadrupel,i_ir,i_if1,i_c1,rel)
!!$                     do i_if2 = 1, n_if2
!!$                        ! the factor two in the following accounts for real and imaginary part
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and. (i_spin1.eq.i_spin2)) then
!!$                           nn_c2 = 2*i_c1
!!$                        else
!!$                           nn_c2 = 2*n_dim2
!!$                        endif
!!$                        call communpack(integralstore(i_meta:i_meta+nn_c2-1),nn_c2,1,info)
!!$                        if (info .ne. 0) call error_handler( &
!!$                             "int_send_2cob3c_spor_receive: unpacking 2 center failed")
!!$                        i_meta = i_meta + nn_c2
!!$                        number_read = number_read + nn_c2
!!$                     enddo ! i_if2
!!$                  enddo ! i_spin2
!!$               enddo ! i_c1
!!$            enddo ! i_if1
!!$         enddo ! i_spin1
!!$      enddo ! i_ir
    end subroutine unpack_integralfile_2c


    subroutine unpack_integralfile_3c(integralstore,n_receive_ff)
      !------------ Declaration of formal parameters -----------
      integer(kind=i4_kind), intent(in) :: n_receive_ff
!!$      real(kind=r8_kind), intent(out) :: integralstore(:)
      real(kind=r8_kind) :: integralstore(:) ! avoid warnings
      !** End of interface *************************************
      !------------ Declaration of local variables -------------
      integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
           n_if1,n_if2,nn_c2,i_meta,i_c2,i_spin1,i_spin2
      !------------ Executable code ----------------------------

      call error_handler("unpack_integralfile_3c: need an update")
!!$      do i_ir = 1,symmetry_data_n_proj_irreps()
!!$         do i_spin1 = 1,2
!!$            if ((quadrupel%l1.eq.0).and.(i_spin1.eq.1)) then
!!$               cycle
!!$            endif
!!$            n_if1 = ua1%symadapt_spor_partner(i_ir,quadrupel%l1,i_spin1)%N_independent_fcts
!!$            do i_if1 = 1, n_if1
!!$               do i_c1 = 1, n_c1
!!$                  do i_spin2 = 1,2
!!$                     if ((quadrupel%l2.eq.0).and.(i_spin2.eq.1)) then
!!$                        cycle
!!$                     endif
!!$                     if (diagonal.and.(i_spin2.gt.i_spin1)) then
!!$                        cycle
!!$                     endif
!!!                     i_meta = (metaindex(quadrupel,i_ir,i_spin1,i_spin2,i_if1,i_c1) - 1) * n_receive_ff + 1
!!$                     i_meta = (metaindex(quadrupel,i_ir,i_if1,i_c1) - 1) * n_receive_ff + 1
!!$                     if ( diagonal.and.(i_spin1.eq.i_spin2) ) then
!!$                        n_if2 = i_if1
!!$                     else
!!$                        n_if2 = ua2%symadapt_spor_partner(i_ir,quadrupel%l2,i_spin2)%N_independent_fcts
!!$                     endif
!!$                     do i_if2 = 1, n_if2
!!$                        if ( diagonal .and. (i_if1 .eq. i_if2).and.(i_spin1.eq.i_spin2) ) then
!!$                           nn_c2 = i_c1
!!$                        else
!!$                           nn_c2 = n_c2
!!$                        endif
!!$                        do i_c2 = 1, nn_c2
!!$                           call communpack( integralstore(i_meta:i_meta+2*n_receive_ff-1), &
!!$                                2*n_receive_ff, 1, info)
!!$                           if (info .ne. 0) call error_handler( &
!!$                                "int_send_2cob3c_spor_receive: unpacking 3 center failed")
!!$                           i_meta = i_meta + 2*n_receive_ff
!!$                        enddo! i_c2
!!$                     enddo! i_if2
!!$                  enddo! i_spin2
!!$               enddo! i_c1
!!$            enddo! i_if1
!!$         enddo! i_spin1
!!$      enddo ! i_ir

    end subroutine unpack_integralfile_3c
#endif
  end subroutine int_send_2cob3c_sp_receive_mem
  !*************************************************************



  !*************************************************************
  character(len=5) function ua_l_fileindex(i_ua,i_l)
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in)  :: i_ua,i_l
    !** End of interface ***************************************
    write(ua_l_fileindex,'(i5)') ua_fileindex(i_ua) + i_l
    ua_l_fileindex = adjustl(ua_l_fileindex)
  end function ua_l_fileindex
  !*************************************************************



  !*************************************************************
  character(len=28) function quadrupel_filename(i_ir,q)
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in) :: i_ir
    type(quadrupel_type), intent(in)  :: q
    !** End of interface ***************************************
    character(len=4) :: charnbr
    write (charnbr,'(i4)') i_ir !!!!!!!!!!!!!!!!!! no free format for DEC f90
    charnbr = adjustl(charnbr)
    quadrupel_filename = repeat(" ",28)
    quadrupel_filename = "quadrupel_" // trim(charnbr) // &
         "_" // trim(ua_l_fileindex(q%ua1,q%l1)) // "_" // &
         trim(ua_l_fileindex(q%ua2,q%l2)) // ".dat"
  end function quadrupel_filename
  !*************************************************************



  !*************************************************************
  integer(kind=i4_kind) function index_ua_l(i_ua,i_l)
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in)  :: i_ua,i_l
    !** End of interface ***************************************
    index_ua_l = ua_fileindex(i_ua) + i_l
  end function index_ua_l
  !*************************************************************



  !*************************************************************
!!$  integer(kind=i4_kind) function metaindex(quadrupel,i_ir,i_spin1,i_spin2,i_ind1,i_exp1,rel_dummy)
  integer(kind=i4_kind) function metaindex(quadrupel,i_ir,i_ind1,i_exp1,rel_dummy)
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type),  intent(in)  :: quadrupel
    integer(kind=i4_kind), intent(in)  :: i_ir,i_ind1,i_exp1 !!$,i_spin1,i_spin2
    logical,optional :: rel_dummy
    !** End of interface ***************************************
    logical :: rel
    rel=.false.
    if(present(rel_dummy)) rel=rel_dummy
!!$    if(.not.rel) then
!!$       metaindex = &
!!$            metaindex_of_first_integral(i_ir)%ua1(quadrupel%ua1)% &
!!$            l1(quadrupel%l1,i_spin1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
!!$            l2(quadrupel%l2,i_spin2)
!!$    else
!!$       metaindex = &
!!$            metaindex_of_first_integral_rel(i_ir)%ua1(quadrupel%ua1)% &
!!$            l1(quadrupel%l1,i_spin1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
!!$            l2(quadrupel%l2,i_spin2)
!!$    endif
    if(.not.rel) then
       metaindex = &
            metaindex_of_first_integral(i_ir)%ua1(quadrupel%ua1)% &
            l1(quadrupel%l1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
            l2(quadrupel%l2)
    else
       metaindex = &
            metaindex_of_first_integral_rel(i_ir)%ua1(quadrupel%ua1)% &
            l1(quadrupel%l1)%ind1_exp1_ua2(i_ind1,i_exp1,quadrupel%ua2)% &
            l2(quadrupel%l2)
    endif
  end function metaindex
  !*************************************************************



  !*************************************************************
  subroutine indices(i_meta,i_ir,i_ua1,i_l1,i_ind1,i_exp1, &
       i_ua2,i_l2,i_ind2,i_exp2)
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in)  :: i_meta
    integer(kind=i4_kind)::  i_ir,i_ua1,i_l1, &
         i_ind1,i_exp1,i_ua2,i_l2,i_ind2,i_exp2
    !** End of interface ***************************************
    integer(kind=i4_kind)  :: n_ind2, n_l2, n_exp2, i_m
    type(unique_atom_type),             pointer :: ua1,ua2
    type(unique_atom_basis_type),       pointer :: uab1,uab2
    type(unique_atom_partner_type),     pointer :: uap1,uap2
    !------------ Executable code ------------------------------
    i_m = 1
    do i_ir = 1,symmetry_data_n_irreps()
       do i_ua1 = 1, N_unique_atoms
          ua1 => unique_atoms(i_ua1)
          do i_l1 = 0, ua1%lmax_ob
             uab1 => ua1%l_ob(i_l1)
             uap1 => ua1%symadapt_partner(i_ir,i_l1)
             do i_ind1 = 1, uap1%N_independent_fcts
                do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                   do i_ua2 = 1, i_ua1
                      ua2 => unique_atoms(i_ua2)
                      if ( i_ua1 .eq. i_ua2 ) then
                         n_l2 = i_l1
                      else
                         n_l2 = ua2%lmax_ob
                      endif
                      do i_l2 = 0, n_l2
                         uab2 => ua2%l_ob(i_l2)
                         uap2 => ua2%symadapt_partner(i_ir,i_l2)
                         if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 ) then
                            n_ind2 = i_ind1
                         else
                            n_ind2 = uap2%N_independent_fcts
                         endif
                         do i_ind2 = 1, n_ind2
                            if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                 .and. i_ind1 .eq. i_ind2 ) then
                               n_exp2 = i_exp1
                            else
                               n_exp2 = uab2%N_uncontracted_fcts + &
                                    uab2%N_contracted_fcts
                            endif
                            do i_exp2 = 1, n_exp2
                               if ( i_m .eq. i_meta ) return
                               i_m = i_m + 1
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine indices
  !*************************************************************


  !*************************************************************
  subroutine int_send_2cob3c_s_rec_filesizes()
    !  Purpose: During shutdown procedure of 2cob3c integral part,
    !    if the file system is parallel all host have to
    !    send to each other information about file sizes they
    !    wrote by sending msgtag_int_2cob3c_filesize.
    !    This messagetag is received by host either while
    !    beeing in subroutine integral_shutdown_2cob3c()
    !    ( => This subroutine is called ) or while executing
    !    int_send_2cob3c_shutdown().
    !  called by: integral_shutdown_2cob3c,
    !** End of interface *****************************************
    integer(kind=i4_kind),allocatable  :: quadrupelfile_length_loc(:,:,:)
    integer(kind=i4_kind) :: info, i_ir, i_file2, alloc_stat
    allocate(quadrupelfile_length_loc( &
         size(quadrupelfile_length,1), &
         size(quadrupelfile_length,2), &
         size(quadrupelfile_length,3)),stat=alloc_stat )
        if(alloc_stat/=0) call error_handler(&
        "int_send_2cob3c_rec_filesizes: allocation failed")
    do i_ir = 1, size(quadrupelfile_length,3)
       do i_file2 = 1, size(quadrupelfile_length,2)
          call communpack( &
               quadrupelfile_length_loc(:,i_file2,i_ir), &
               size(quadrupelfile_length,1), 1, info )
          if (info .ne. 0) call error_handler( &
               "int_send_2cob3c_rec_filesizes: unpacking failed")
       enddo
    enddo
    quadrupelfile_length = &
         quadrupelfile_length + quadrupelfile_length_loc
    n_missing_filesizes = n_missing_filesizes - 1
    deallocate(quadrupelfile_length_loc,stat=alloc_stat )
        if(alloc_stat/=0) call error_handler(&
        "int_send_2cob3c_rec_filesizes: deallocation failed")
      end subroutine int_send_2cob3c_s_rec_filesizes
  !*************************************************************

!--------------- End of module ----------------------------------
end module int_send_2cob3c_spor_module

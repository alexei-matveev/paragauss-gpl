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
module  int_send_2cob3c_module
!---------------------------------------------------------------
!
!  Purpose: Sending, receiving and storing of
!           2 center orbital and three center integrals.
!           The results of one quadrupel stored in
!           int_data_2cob3c_module are distributed.
!           They are either stored in intermediate
!           quadrupel files or directly stored in memory
!           in integralstore_module. In case files are used,
!           when executing int_send_2cob3c_shutdown(), the
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
! Modification (Please copy before editing)
! Author: Uwe Birkenheuer
! Date:   June 1998
! Description: Moving_Unique_Atom concept introduced
!              Split_Gradients concept introduced
!              Gradients for Model_Density_Approach introduced
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: pvm -> comm
!
! Modification (Please copy before editing)
! Author: AS
! Date:   11-12/99
! Description: integrals of electrostatic potential are added
!
! Modification Master/Slave concept to DLB
! Author: AN
! Date:   4/11
! Description: for scheduling the 2cob3c integrals DLB is used
!              it replaces the master/slave concept
!              This removes all messages related to do or done of jobs
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!define FPP_NO_CHECKSUM
#include "def.h"
use type_module ! type specification parameters
#ifdef _COMPAC_FORTRAN
use datatype
#endif
use integralpar_module,   only: integralpar_2cob_kin                           &
                              , integralpar_2cob_nuc                           &
                              , integralpar_2cob_field                         &
                              , integralpar_2cob_ol                            &
                              , integralpar_2cob_potential                     &
                              , integralpar_2cob_pvsp                          &
                              , integralpar_3c_co                              &
                              , integralpar_3c_xc                              &
                              , integralpar_relativistic                       &
                              , integralpar_gradients                          &
                              , integralpar_i_int_part
use output_module,        only: output_int_loops                               &
                              , output_int_deeploops                           &
                              , output_int_progress                            &
                              , output_slaveoperations
use iounitadmin_module
use filename_module
use comm_module
use msgtag_module,        only: msgtag_int_2cob3c_filesize                     &
                              , msgtag_int_2cob3c_result
use symmetry_data_module, only: symmetry_data_n_irreps
use unique_atom_module ! FIXME: only?
use quadrupel_module,     only: quadrupel_type                                 &
                              , quadrupel_pack                                 &
                              , quadrupel_unpack
use timer_module
use time_module
use readwriteblocked_module
use options_module,       only: options_integrals_on_file
use integralstore_module, only: integralstore_allocate                         &
                              , integralstore_deallocate                       &
                              , integralstore_allocate_pcm                     &
                              , integralstore_3c_xc                            &
                              , integralstore_3c_co                            &
                              , integralstore_3c_poten                         &
                              , integralstore_3c_field                         &
                              , integralstore_2cob_pvsp                        &
                              , integralstore_2cob_ol                          &
                              , integralstore_2cob_ol_rel                      &
                              , integralstore_2cob_nuc                         &
                              , integralstore_2cob_nuc_rel                     &
                              , integralstore_2cob_kin                         &
                              , integralstore_2cob_kin_rel
use gradient_data_module
use potential_module
use elec_static_field_module
use error_module

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
public int_send_2cob3c_setup, int_send_2cob3c_shutdown, &
       int_send_2cob3c_send,                          &
                                 int_send_2cob3c_rec_filesizes

public :: int_send_2cob3c_receive
public :: int_send_2cob3c_receive_all

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


!------------ Declaration of constants and variables ----
integer(kind=i4_kind) :: n_missing_quadrupels, &
     first_host, last_host, my_hostindex, n_file, n_tot_records, &
     n_tot_records_rel, &
     n_hosts_to_report_back, blocklength, n_missing_filesizes
integer(kind=i4_kind), dimension(:,:), allocatable :: &
     borders_ch, borders_xc, borders_poten, borders_field ! borders(3,comm_get_n_hosts())
     ! 1. dim: 1: lower border, 2: upper border, 3: number of ff
integer(kind=i4_kind), allocatable :: ua_fileindex(:) ! (i_ua)
     ! to look up file index for l = 0
integer(kind=i4_kind), allocatable :: quadrupelfile_length(:,:,:), &
     quadrupelfile_length_allhosts(:,:,:,:)
type(index_ir_type), allocatable :: &
     metaindex_of_first_integral(:) ! (i_ir)
     ! to look up metaindex of first integral and number
     ! n_records of following records for a given vector
     ! (i_ir,i_ua1,i_l1,i_ind1,i_exp2,i_ua2,i_l2)
type(index_ir_type), allocatable :: &
     metaindex_of_first_integral_rel(:) ! (i_ir)
     ! to look up metaindex of first integral in
     ! a uncontracted basis
     ! (i_ir,i_ua1,i_l1,i_ind1,i_exp2,i_ua2,i_l2)

#ifdef _SR8000
character(len=6), parameter :: DELETE='KEEP'
#else
character(len=6), parameter :: DELETE='DELETE'
#endif
integer(kind=i4_kind) :: da_rec_no
integer(kind=i4_kind) :: dag_rec_no(5)
integer(kind=i4_kind) :: dag_unit(5)
integer(kind=i4_kind), allocatable, dimension(:,:,:) :: darec_no

type(readwriteblocked_tapehandle),target :: th_da
!type(readwriteblocked_tapehandle),pointer :: th_da
type(readwriteblocked_tapehandle),pointer :: thdaf_pointer

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine int_send_2cob3c_setup(n_quad)
    !  Purpose:
    !   + calculates fitfunction index boundaries for
    !     splitting 3 center integrals
    !   + calculates help variables number of quadrupels
    !     to be received and total number of record in file
    !   + calculates help array for mapping integrals to
    !     1 dimensional direct access file.
    !   + Opening of direct access files.
    ! called by: integral_setup_2cob3c (on all hosts)
    !** End of interface *****************************************
    !------------ Modules used ---------------------------------
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc
    use options_module, only: options_directaccess_integrals,quadrupels_reclength
    implicit none
    integer(kind=i4_kind), intent(in)  :: n_quad
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: i_ir, n_hosts, status, n_hosts_larger, &
         n_per_host, i_border, i_host, n_ch, n_xc, n_poten, n_field, &
         i_ua1, i_ua2, &
         i_file, i_l1, n_l2, allocsize, basesize, n_free_iounits, &
         i_meta, i_ind1, i_exp1, i_l2, n_ind2, n_records, &
         i_ind2, n_exp2
    real(kind=r8_kind), allocatable :: dummy(:)
    type(unique_atom_type),             pointer :: ua1,ua2
    type(unique_atom_basis_type),       pointer :: uab1,uab2
    type(unique_atom_partner_type),     pointer :: uap1,uap2
    type(index_ind1_exp1_ua2_type),     pointer :: ind1_exp1_ua2(:,:,:)
    integer,                            pointer :: metaindex(:)
!! norge
    !------------ Declaration of formal parameters for function integralstore_allocate-------------
    integer(kind=i4_kind)  :: &
           dim_2c, dim_2c_rel, dim_3c_co, dim_3c_xc, dim_3c_poten, dim_3c_field
    logical  :: &
           need_2cob_kin, need_2cob_nuc, need_2cob_kin_rel, need_2cob_nuc_rel, &
           need_2cob_pvsp, need_2cob_pvxp, need_2cob_ol, &
           need_2cob_ol_rel, need_2cob_pvec, need_3c_xc, need_3c_co, &
           need_3c_poten, need_3c_field
!! norge end
    !------------ Executable code ------------------------------------


    ! Nothing to be communicated for gradeints
    if( integralpar_gradients )then
       n_missing_quadrupels = 0
       RETURN
    endif

    n_missing_quadrupels = n_quad

    if ( output_int_progress ) then
       call write_to_output_units("")
       call write_to_output_units( &
            "total number of 2 center orbital and 3 center quadrupels: ",n_quad)
       call write_to_output_units("")
    endif

    ! calculate borders
    first_host = 1
    last_host = comm_get_n_processors()
    n_hosts = last_host - first_host + 1
    if ( n_hosts .lt. 1 ) call error_handler( &
         "int_send_2cob3c_setup: wrong n_hosts" )

       allocate( borders_ch(3,first_host:last_host), &
            borders_xc(3,first_host:last_host),  stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_setup: allocate of borders failed" )

       allocate( borders_poten(3,first_host:last_host), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "int_send_2cob3c_setup: allocate of poten borders failed" )

       allocate( borders_field(3,first_host:last_host), stat=status )
       if ( status .ne. 0 ) call error_handler( & 
            "int_send_2cob3c_setup: allocate of field borders failed" )

    if (integralpar_3c_co) then
       n_ch = fit_coeff_n_ch()
       n_hosts_larger = mod(n_ch,n_hosts)
       n_per_host = n_ch / n_hosts
       i_border = 1
       do i_host = first_host, first_host + n_hosts_larger - 1
          borders_ch(1,i_host) = i_border
          i_border = i_border + n_per_host
          borders_ch(2,i_host) = i_border
          i_border = i_border + 1
          borders_ch(3,i_host) = n_per_host + 1
       enddo
       do i_host = first_host + n_hosts_larger, last_host
          borders_ch(1,i_host) = i_border
          i_border = i_border + n_per_host - 1
          borders_ch(2,i_host) = i_border
          i_border = i_border + 1
          borders_ch(3,i_host) = n_per_host
       enddo
       if ( output_int_loops .and. output_unit > 0 ) then
          write (output_unit,*)
          write (output_unit,*) "##### borders for distribution of ch integrals to hosts #####"
          write (output_unit,*) "i_host  from    to     n"
          do i_host = first_host, last_host
             write (output_unit,'(4I6)') i_host, borders_ch(1,i_host), &
                  borders_ch(2,i_host), borders_ch(3,i_host)
          enddo
          write (output_unit,*) "### borders for distribution of ch integrals to hosts end ###"
          write (output_unit,*)
       endif
    else
       borders_ch = 0
    endif

    if (integralpar_3c_xc) then
       n_xc = fit_coeff_n_xc()
       n_hosts_larger = mod(n_xc,n_hosts)
       n_per_host = n_xc / n_hosts
       i_border = 1
       do i_host = first_host, first_host + n_hosts_larger - 1
          borders_xc(1,i_host) = i_border
          i_border = i_border + n_per_host
          borders_xc(2,i_host) = i_border
          i_border = i_border + 1
          borders_xc(3,i_host) = n_per_host + 1
       enddo
       do i_host = first_host + n_hosts_larger, last_host
          borders_xc(1,i_host) = i_border
          i_border = i_border + n_per_host - 1
          borders_xc(2,i_host) = i_border
          i_border = i_border + 1
          borders_xc(3,i_host) = n_per_host
       enddo
       if ( output_int_loops .and. output_unit > 0 ) then
          write (output_unit,*)
          write (output_unit,*) "##### borders for distribution of xc integrals to hosts #####"
          write (output_unit,*) "i_host  from    to     n"
          do i_host = first_host, last_host
             write (output_unit,'(4I6)') i_host, borders_xc(1,i_host), &
                  borders_xc(2,i_host), borders_xc(3,i_host)
          enddo
          write (output_unit,*) "### borders for distribution of xc integrals to hosts end ###"
          write (output_unit,*)
       endif
    else
       borders_xc = 0
    endif

    if (integralpar_2cob_potential) then
       n_poten = N_points
       n_hosts_larger = mod(n_poten,n_hosts)
       n_per_host = n_poten / n_hosts
       i_border = 1
       do i_host = first_host, first_host + n_hosts_larger - 1
          borders_poten(1,i_host) = i_border
          i_border = i_border + n_per_host
          borders_poten(2,i_host) = i_border
          i_border = i_border + 1
          borders_poten(3,i_host) = n_per_host + 1
       enddo
       do i_host = first_host + n_hosts_larger, last_host
          borders_poten(1,i_host) = i_border
          i_border = i_border + n_per_host - 1
          borders_poten(2,i_host) = i_border
          i_border = i_border + 1
          borders_poten(3,i_host) = n_per_host
       enddo
       if ( output_int_loops .and. output_unit > 0 ) then
          write (output_unit,*)
          write (output_unit,*) "##### borders for distribution of poten integrals to hosts #####"
          write (output_unit,*) "i_host  from    to     n"
          do i_host = first_host, last_host
             write (output_unit,'(4I6)') i_host, borders_poten(1,i_host), &
                  borders_poten(2,i_host), borders_poten(3,i_host)
          enddo
          write (output_unit,*) "### borders for distribution of poten integrals to hosts end ###"
          write (output_unit,*)
       endif
    else
       borders_poten = 0
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (integralpar_2cob_field) then
       n_field = N_surface_points
       n_hosts_larger = mod(n_field,n_hosts)
       n_per_host = n_field / n_hosts
       i_border = 1
       do i_host = first_host, first_host + n_hosts_larger - 1
          borders_field(1,i_host) = surf_points_grad_index(i_border)
          i_border = i_border + n_per_host
          borders_field(2,i_host) = surf_points_grad_index(i_border+1)-1
          i_border = i_border + 1
          borders_field(3,i_host) = borders_field(2,i_host)-borders_field(1,i_host)+1
       end do
       do i_host = first_host + n_hosts_larger, last_host
          borders_field(1,i_host) = surf_points_grad_index(i_border)
          i_border = i_border + n_per_host-1
          borders_field(2,i_host) = surf_points_grad_index(i_border+1)-1
          i_border = i_border + 1
          borders_field(3,i_host) = borders_field(2,i_host)-borders_field(1,i_host)+1
       end do
       if ( output_int_loops .and. output_unit > 0 ) then
          write (output_unit,*)
          write (output_unit,*) "##### borders for distribution of field integrals to hosts #####"
          write (output_unit,*) "i_host  from    to     n"
          do i_host = first_host, last_host
             write (output_unit,'(4I6)') i_host, borders_field(1,i_host), &
                  borders_field(2,i_host), borders_field(3,i_host)
          enddo
          write (output_unit,*) "### borders for distribution of field integrals to hosts end ###"
          write (output_unit,*)
       endif
    else
       borders_field = 0
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! fetch my hostindex from comm_module
    my_hostindex = comm_myindex()



    if ( options_integrals_on_file() ) then



               ! allocate and calculate ua_fileindex
               allocate( ua_fileindex(N_unique_atoms), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_send_2cob3c_setup: allocate of ua_fileindex failed" )
               i_file = 1
               do i_ua1 = 1, N_unique_atoms
                  ua_fileindex(i_ua1) = i_file
                  i_file = i_file + unique_atoms(i_ua1)%lmax_ob + 1
               enddo
               n_file = i_file - 1


               ! allocate quadrupelfile_length
               allocate( quadrupelfile_length(n_file,n_file,symmetry_data_n_irreps()), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_send_2cob3c_setup: allocate of quadrupelfile_length failed" )
               quadrupelfile_length = 0
               if ( filesystem_is_parallel ) then
                  allocate( quadrupelfile_length_allhosts(n_file, n_file, &
                       symmetry_data_n_irreps(), comm_get_n_processors() ), &
                       stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "int_send_2cob3c_setup: allocate of quadrupelfile_length_allhosts failed" )
                  quadrupelfile_length = 0
                  quadrupelfile_length_allhosts = 0
                  n_missing_filesizes = comm_get_n_processors()
               endif


               ! blocklength for intermediate files to store quadrupels
               ! Test if sufficient memory is available, if not reduce blocksize
               ! take into account that swap (on a typical system, roughly
               ! two times hardware memory) will normally be used for allocation
               ! as well. This is all rough estimates only. The memory test will not
               ! work this way on a shared memory architecture.
               blocklength = readwriteblocked_blocklength()
               n_free_iounits = get_nbr_free_iounits() * 3 ! swap
               basesize = 10 * blocklength
               allocsize = basesize + n_free_iounits * blocklength
!!!!!AG               allocate( dummy(allocsize), stat=status ) ! the attempt to avoid large memory jumps
               status=0
               do while (status .ne. 0)
                  blocklength = blocklength / 2
                  if ( blocklength .lt. 10 ) call error_handler( &
                       "int_send_2cob3c_setup: insufficient memory for rewriting" )
                  allocsize = basesize + n_free_iounits * blocklength
                  allocate( dummy(allocsize), stat=status )
               enddo
!!!!!AG        deallocate( dummy, stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_send_2cob3c_setup: deallocate of dummy failed" )
               if ( output_int_loops .and. output_unit > 0 ) then
                  write (output_unit,*)
                  write (output_unit,*) "using blocklength ",blocklength," for quadrupel files"
                  write (output_unit,*)
               endif

            else ! .not.  options_integrals_on_file()

                  ! calculate index_of_ua_l
                  allocate( metaindex_of_first_integral(symmetry_data_n_irreps()), &
                       stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "send_2cob3c_setup: allocate of metaindex_of_first_integral failed" )
                  i_meta = 1
                  do i_ir = 1,symmetry_data_n_irreps()
                     allocate( metaindex_of_first_integral(i_ir)%ua1(N_unique_atoms), &
                          stat=status )
                     if ( status .ne. 0 ) call error_handler( &
                          "send_2cob3c_setup: allocate of metaindex_of_first_integral ua1 failed" )
                     do i_ua1 = 1, N_unique_atoms
                        ua1 => unique_atoms(i_ua1)
                        allocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                             stat=status )
                        if ( status .ne. 0 ) call error_handler( &
                             "send_2cob3c_setup: allocate of metaindex_of_first_integral l1 failed" )
                        do i_l1 = 0, ua1%lmax_ob
                           uab1 => ua1%l_ob(i_l1)
                           uap1 => ua1%symadapt_partner(i_ir,i_l1)
                           allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                                uab1%N_uncontracted_fcts + uab1%N_contracted_fcts, i_ua1), &
                                stat=status )
                           if ( status .ne. 0 ) call error_handler( &
                                "send_2cob3c_setup: allocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
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
                                         "send_2cob3c_setup: allocate of metaindex_of_first_integral l2 failed" )
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
                  n_tot_records = i_meta - 1
                  if ( output_int_loops ) call write_to_output_units( &
                       "send_2cob3c_setup: total number of records: ",n_tot_records)

               if(integralpar_relativistic) then
                  ! calculate index_of_ua_l
                  allocate( metaindex_of_first_integral_rel(symmetry_data_n_irreps()), &
                       stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "send_2cob3c_setup: allocate of metaindex_of_first_integral_rel failed" )
                  i_meta = 1
                  do i_ir = 1,symmetry_data_n_irreps()
                     ! we want to have every irrep in a seperate array
                     allocate( metaindex_of_first_integral_rel(i_ir)%ua1(N_unique_atoms), &
                          stat=status )
                     if ( status .ne. 0 ) call error_handler( &
                          "send_2cob3c_setup: allocate of metaindex_of_first_integral_rel ua1 failed" )
                     do i_ua1 = 1, N_unique_atoms
                        ua1 => unique_atoms(i_ua1)
                        allocate( metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(0:ua1%lmax_ob), &
                             stat=status )
                        if ( status .ne. 0 ) call error_handler( &
                             "send_2cob3c_setup: allocate of metaindex_of_first_integral_rel l1 failed" )
                        do i_l1 = 0, ua1%lmax_ob
                           uab1 => ua1%l_ob(i_l1)
                           uap1 => ua1%symadapt_partner(i_ir,i_l1)
                           allocate( ind1_exp1_ua2(uap1%N_independent_fcts, &
                                uab1%N_exponents, i_ua1), &
                                stat=status )
                           if ( status .ne. 0 ) call error_handler( &
                                "send_2cob3c_setup: allocate of metaindex_of_first_integral_rel ind1_exp1_ua2 failed" )
                           metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2 => ind1_exp1_ua2
                           do i_ind1 = 1, uap1%N_independent_fcts
                              do i_exp1 = 1, uab1%N_exponents
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
                                         "send_2cob3c_setup: allocate of metaindex_of_first_integral_rel l2 failed" )
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
                                             n_exp2 = uab2%N_exponents
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
                     enddo! i_ua1
                  enddo! irreps
               end if! integralpar_relativistic
               n_tot_records_rel = i_meta - 1
               if ( output_int_loops ) call write_to_output_units( &
                    "send_2cob3c_setup: total number of records_rel: ",n_tot_records_rel)
               ! allocating storage for integrals
               if(.not.integralpar_relativistic) n_tot_records_rel=n_tot_records

               if( .not. integralpar_2cob_potential &
                    .and. .not. integralpar_2cob_field) then
!! norge
                       dim_2c = n_tot_records
                       dim_2c_rel = n_tot_records_rel
                       dim_3c_co = n_tot_records * borders_ch(3,my_hostindex)
                       dim_3c_xc = n_tot_records * borders_xc(3,my_hostindex)
                       need_2cob_kin = integralpar_2cob_kin .and. comm_i_am_master()
                       need_2cob_nuc = integralpar_2cob_nuc .and. comm_i_am_master()
                       need_2cob_kin_rel = integralpar_2cob_kin .and. comm_i_am_master() &
                       .and. integralpar_relativistic
                       need_2cob_nuc_rel = integralpar_2cob_nuc .and. comm_i_am_master() &
                       .and. integralpar_relativistic
                       need_2cob_pvsp = integralpar_2cob_pvsp .and. comm_i_am_master().and. &
                       integralpar_relativistic
!!!!!!!!               need_2cob_pseudo = pseudopot_present .and. comm_i_am_master().and. &
!!!!!!!!               integralpar_relativistic
                       need_2cob_pvxp = .false.
                       need_2cob_ol = integralpar_2cob_ol
                       need_2cob_ol_rel = integralpar_2cob_ol .and. &
                       integralpar_relativistic
                       need_2cob_pvec = .false.    !<<<mdf AM???
                       need_3c_xc = integralpar_3c_xc
                       need_3c_co = integralpar_3c_co

                  ! for the case of cases:
                  call integralstore_deallocate()

                  if ( integralstore_allocate ( &
                       dim_2c = dim_2c , &
                       dim_2c_rel=dim_2c_rel , &
                       dim_3c_co=dim_3c_co , &
                       dim_3c_xc=dim_3c_xc , &
                       need_2cob_kin=need_2cob_kin , &
                       need_2cob_nuc=need_2cob_nuc , &
                       need_2cob_kin_rel=need_2cob_kin_rel , &
                       need_2cob_nuc_rel=need_2cob_nuc_rel , &
                       need_2cob_pvsp=need_2cob_pvsp , &
!               need_2cob_pseudo=need_2cob_pseudo , &
                       need_2cob_pvxp=need_2cob_pvxp , &
                       need_2cob_ol=need_2cob_ol , &
                       need_2cob_ol_rel=need_2cob_ol_rel , &
                       need_2cob_pvec=need_2cob_pvec , &
                       need_3c_xc=need_3c_xc , &
                       need_3c_co=need_3c_co &
                       )  ) &
                       call error_handler( &
                       "send_2cob3c_setup: allocation of integral storage failed. set options_integrals_on_file .true.")
!! norge end
               else 
!! norge 
                      dim_2c = n_tot_records
                      dim_3c_poten=n_tot_records*borders_poten(3,my_hostindex)
                      dim_3c_field=n_tot_records*borders_field(3,my_hostindex)
                      need_3c_poten = &
                             integralpar_2cob_potential
                      need_3c_field = &
                             integralpar_2cob_field
                  if ( integralstore_allocate_pcm(&
                      dim_2c = dim_2c , &
                      dim_3c_poten=dim_3c_poten,&
                      dim_3c_field=dim_3c_field,&
                      need_3c_poten=need_3c_poten , &
                      need_3c_field=need_3c_field ) ) &

                       call error_handler( &
                       "send_2cob3c_setup: allocation of integral storage pcm failed. &
                       &set options_integrals_on_file .true.")
               end if
!! norge end
            endif!  options_integrals_on_file()/else


            ! during shutdown procedure, all slaves have to report back
            n_hosts_to_report_back = comm_get_n_processors() - 1

            if(options_directaccess_integrals() ) then
            print*,'check where da_buff is deallocated'
            da_rec_no=0
            call readwriteblocked_startwrite(trim(tmpfile("da")), th_da, &
                                    blocklength=quadrupels_reclength,rec=da_rec_no )
            thdaf_pointer=> th_da
            allocate(darec_no(symmetry_data_n_irreps(),n_file,n_file)) !dealloc in shutdown
            endif


  end subroutine int_send_2cob3c_setup
  !*************************************************************


          !*************************************************************
          subroutine int_send_2cob3c_shutdown
            !  Purpose:
            !    + waiting for missing integrals
            !    + deallocation of help arrays
            !    + reading of quadrupel files and writing
            !      them as sequential files as used in the scf part
            !      with the files produced by old lcgto
            !    + removal of quadrupel files
            ! called by: integral_shutdown_2cob3c (on all hosts)
            !** End of interface *****************************************
            use operations_module
            use options_module, only: options_directaccess_integrals,quadrupels_reclength
            implicit none

            !------------ Declaration of local types ---------------------
            type quadrupel_file_handle_type
               type(readwriteblocked_tapehandle) :: th
               integer(kind=i4_kind) :: n_if_c, total_length
               integer(kind=i4_kind), pointer ::  n_floats_if_c(:)
               ! the next elements are only used in the relativistic case
               integer(kind=i4_kind) :: n_if_rel
               integer(kind=i4_kind), pointer ::  n_floats_if_rel(:)
            end type quadrupel_file_handle_type

            !------------ Declaration of local variables -----------------
            integer(kind=i4_kind) :: status, i_ir, unit_kin_nuc, unit_ol, i_unique,i_l,&
                 i_file1, i_file2, n_file2, n_integrals, n_free_iounits, &
                 i_file2_new, n_file2_new, n_concat, n_tot_records, n_tot_records_rel, &
                 info, i_integral, i_if_c, n_file2_conc, n_files_left, n_if_c, n_irrep, &
                 i_host, i_ua1, i_l1, i_ind1, i_exp1, i_ua2
!!$            integer(kind=i4_kind) :: i
            integer(kind=i4_kind) :: da_unit
!!$            logical:: da_status
            integer(kind=i4_kind), allocatable :: blocksizes(:), &
                 ! blocksizes(n_integrals)
            dimensions(:) ! dimensions of irreps, only used for relativistic

            type(readwriteblocked_tapehandle), allocatable :: &
                 th_kin_arr(:), th_nuc_arr(:), &
!!!!!!!!         th_nuc_pseudo_arr(:),&
                 th_pvsp_arr(:), th_ol_arr(:) ! th_..(n_irrep), only used for relativistic

            character(len=5) :: ch_i_ir, ch_i_file1, ch_i_file2
            type(quadrupel_file_handle_type), allocatable :: qfh(:)
            type(quadrupel_file_handle_type) :: qfh_da
            type(readwriteblocked_tapehandle) :: th_xc, th_co, th_ol, th_kin, th_nuc,  th_pvsp
            type(readwriteblocked_tapehandle) :: th_poten
            type(readwriteblocked_tapehandle) :: th_field
!!!!!!!!    type(readwriteblocked_tapehandle) :: th_nuc_pseudo !split nuc atraction m. el.

            logical :: restart_timer, rel
            logical, allocatable :: rel_flag(:)
            character(len=3)                    :: inter
            ! internal file
!!$            character(len=4)                    :: inter1
            ! internal file

            type(unique_atom_type),             pointer :: ua1
            type(unique_atom_basis_type),       pointer :: uab1
            type(unique_atom_partner_type),     pointer :: uap1
            type(index_ind1_exp1_ua2_type),     pointer :: ind1_exp1_ua2(:,:,:)
            !------------ Executable code --------------------------------

            ! Nothing to be done for gradeints
            if( integralpar_gradients )then
               RETURN
            endif

            if( n_missing_quadrupels /= 0 )then
               ABORT('call _receive_all first!')
            endif

            do_integrals: if ( (operations_integral .or. operations_properties) ) then

               filesizes: if ( filesystem_is_parallel &
                    .and. options_integrals_on_file() ) then

                  ![[=== send/recv information about file sizes to other hosts ===
                  if( comm_parallel() )then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_shutdown: sending file sizes: ")
                  do i_host = 1, comm_get_n_processors()
                     ! THIS LOOP IS EXECUTED n_proc - 1 TIMES:
                     if( i_host == my_hostindex ) cycle

                        DPRINT  MyID,'int_send_2cob3c_shutdown: SEND filesizes to i_host=',i_host
                        ! SENDING TO i_host:
                        call comm_init_send(i_host,msgtag_int_2cob3c_filesize)
                        do i_ir = 1, ubound(quadrupelfile_length,3)
                           do i_file2 = 1, size(quadrupelfile_length,2)
                              call commpack( &
                                   quadrupelfile_length_allhosts(:,i_file2,i_ir,i_host), &
                                   size(quadrupelfile_length,1), 1, info )
                              if (info .ne. 0) call error_handler( &
                                   "int_send_2cob3c_shutdown: packing of file sizes failed")
                           enddo
                        enddo
                        call comm_send() 

                        ! RECEIVING FROM i_host (probably). ARRIVAL ORDER NOT GUARANTEED!
                        if ( output_int_loops ) call write_to_output_units( &
                           "int_send_2cob3c_shutdown: waiting for file sizes: ",inte=n_missing_filesizes)

                        call comm_save_recv(comm_all_other_hosts,msgtag_int_2cob3c_filesize)

                        if ( output_slaveoperations ) &
                             call write_to_output_units("int_send_2cob3c_shutdown: integral_2cob_filesize")

                        ! FIXME: name is unfortunate, this one UNPACKS i.e. receivs:
                        call int_send_2cob3c_rec_filesizes()
                        DPRINT  MyID,'int_send_2cob3c_shutdown: RECV filesizes indexed by',i_host
                  enddo
                  endif ! comm_parallel()
                  !]]=============================================================

                  DPRINT  MyID,'int_send_2cob3c_shutdown: n_missing_filesizes=',n_missing_filesizes
                  ! add own contrbution
                  n_missing_filesizes = n_missing_filesizes - 1
                  DPRINT  MyID,'int_send_2cob3c_shutdown: DECREMENT n_missing_filesizes to',n_missing_filesizes
                  quadrupelfile_length = &
                       quadrupelfile_length + quadrupelfile_length_allhosts(:,:,:,my_hostindex)

                  ! free memory
                  deallocate( quadrupelfile_length_allhosts, stat=status )
                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: deallocate of quadrupelfile_length_allhosts failed" )
               endif filesizes

            endif do_integrals

           integrals_on_file:  if ( options_integrals_on_file() ) then


               if ( timer_int_idle_2cob3c(integralpar_i_int_part)%running ) then
                  call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
                  restart_timer = .true.
               else
                  restart_timer = .false.
               endif
               call start_timer(timer_int_rewrite_2cob3c(integralpar_i_int_part))

                  ! opening files, for operations integral those used in SCF-part
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_shutdown: opening output files")

                  if ( comm_i_am_master() ) then
!    shutdown
                     if (integralpar_2cob_kin) &
                          call readwriteblocked_startwrite(trim(tmpfile("kin.dat")), th_kin)
                     if ( integralpar_2cob_nuc ) then
                        call readwriteblocked_startwrite(trim(tmpfile("nuc.dat")), th_nuc)
                     endif
                     if (integralpar_2cob_pvsp) &
                          call readwriteblocked_startwrite(trim(tmpfile("pvsp.dat")), th_pvsp)

                  endif! comm_i_am_master

                  if (integralpar_2cob_ol.and.integralpar_relativistic) &
                       call readwriteblocked_startwrite(trim(tmpfile("overlap.dat")), th_ol)

                  if (operations_integral.or.operations_properties) then
                     if ( .not.integralpar_relativistic.and.integralpar_2cob_ol ) &
                          unit_ol = openget_iounit(trim(tmpfile("overlap.dat")), &
                          status='REPLACE', form='UNFORMATTED', action='WRITE')
                     if ( integralpar_3c_xc ) &
                          call readwriteblocked_startwrite(trim(tmpfile("exch.dat")), th_xc)
                     if ( integralpar_3c_co ) &
                          call readwriteblocked_startwrite(trim(tmpfile("coul.dat")), th_co)
                     if ( integralpar_2cob_potential ) &
                          call readwriteblocked_startwrite(trim(tmpfile("poten.dat")), th_poten)
                     if ( integralpar_2cob_field ) &
                          call readwriteblocked_startwrite(trim(tmpfile("field.dat")), th_field)
                  else
                     call error_handler("int_send_2cob3c_shutdown: wrong operations")
                  endif

                  ! calculate the number of integrals n_integrals and blocksizes(n_integrals)
                  ! used in concat_files(...)
                  n_integrals = 0

                  if ( comm_i_am_master()) then
                     if (integralpar_2cob_kin) n_integrals = n_integrals + 1
                     if (integralpar_2cob_nuc) then 
                        n_integrals = n_integrals + 1
                     endif

                     if (integralpar_2cob_pvsp) n_integrals = n_integrals + 1
!!!!!!!!             !not relativistic gradients
!!!!!!!!             if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                n_integrals = n_integrals + 1
!!!!!!!!             endif
                  endif

                  if ( integralpar_2cob_ol ) n_integrals = n_integrals + 1

                  if ( integralpar_3c_xc ) n_integrals = n_integrals + 1
                  if ( integralpar_3c_co ) n_integrals = n_integrals + 1
                  if (integralpar_2cob_potential) n_integrals = n_integrals + 1
                  if (integralpar_2cob_field) n_integrals = n_integrals + 1

                  allocate( blocksizes(n_integrals),rel_flag(n_integrals), stat=status )
                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: allocate of blocksizes failed" )
                  i_integral=1
                  if ( comm_i_am_master()) then   
                     if (integralpar_2cob_kin) then
                        if(integralpar_relativistic) rel_flag(i_integral)=.true.
                        blocksizes(i_integral) = 1
                        i_integral = i_integral + 1
                     endif
                     if (integralpar_2cob_nuc) then
                        if(integralpar_relativistic) rel_flag(i_integral)=.true. 
                        blocksizes(i_integral) = 1
                        i_integral = i_integral + 1
                     endif
                     if (integralpar_2cob_pvsp) then
                        if(integralpar_relativistic) rel_flag(i_integral)=.true. 
                        blocksizes(i_integral) = 1
                        i_integral = i_integral + 1
                     endif
!!!!!!!!             if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                rel_flag(i_integral)=.true.
!!!!!!!!                blocksizes(i_integral) = 1
!!!!!!!!                i_integral = i_integral + 1
!!!!!!!!             endif
                  endif! comm_i_am_master()

                  if ( integralpar_2cob_ol ) then
                     if(integralpar_relativistic) rel_flag(i_integral)=.true. 
                     blocksizes(i_integral) = 1
                     i_integral = i_integral + 1
                  endif
                  if ( integralpar_3c_xc ) then
                     blocksizes(i_integral) = borders_xc(3,my_hostindex)
                     i_integral = i_integral + 1
                  endif
                  if ( integralpar_3c_co ) then
                     blocksizes(i_integral) = borders_ch(3,my_hostindex)
                     i_integral = i_integral + 1
                  endif
                  if (integralpar_2cob_potential) then
                     blocksizes(i_integral) = borders_poten(3,my_hostindex)
                     i_integral = i_integral + 1
                  endif
                  if (integralpar_2cob_field) then
                     blocksizes(i_integral) = borders_field(3,my_hostindex)
                     i_integral = i_integral + 1
                  endif

                  ! + reading of quadrupel integral files and writing
                  !   them as sequential files as used in the scf part
                  ! + removal of quadrupel integral files
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_shutdown: rewriting files")
                  n_free_iounits = get_nbr_free_iounits()

                  if(options_directaccess_integrals()) then
                  allocate( qfh(n_file), stat=status )
                  else
                  allocate( qfh(n_free_iounits), stat=status )
                  endif

                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: allocate of qfh failed" )
                  n_tot_records = 0
                  n_tot_records_rel = 0
                  n_files_left = (n_file + 1) * n_file * symmetry_data_n_irreps() / 2

                  !
                  ! Trace file is not always availabel for slaves:
                  !
                  if ( comm_i_am_master() ) then
                     call write_to_trace_unit("int_send_2cob3c_shutdown: &
                          &total number of files to rewrite: ",inte=n_files_left)
                  endif

                if(options_directaccess_integrals()) then
                 call readwriteblocked_returnclose(th_da)
                           da_rec_no = 0
                endif

                  do i_ir = 1, symmetry_data_n_irreps()
                     write(ch_i_ir,'(i5)') i_ir                !!!!!!!!!!!!!!!!!!!!
                     ch_i_ir = adjustl(ch_i_ir)

                     do i_file1 = 1, n_file
                        write(ch_i_file1,'(i5)') i_file1        !!!!!!!!!!!!!!
                        ch_i_file1 = adjustl(ch_i_file1)

                        n_file2 = i_file1

                        i_file2_new = n_file2
                        n_file2_conc = 0

                        if(.not.options_directaccess_integrals()) then
                        ! in case there are more quadrupel files than free tapehandles,
                        ! rewrite the last tapes to make up one tape          
                        do while ( n_file2 + n_file2_conc .gt. n_free_iounits )
                           i_file2_new = i_file2_new + 1
                           n_concat = n_file2 - n_free_iounits + 1 + n_file2_conc
                           if (n_concat .gt. n_free_iounits-1) n_concat = n_free_iounits - 1
                           call concat_files(n_file2-n_concat+1, n_file2, i_file2_new)
                           n_file2_conc =  n_file2_conc + 1
                           if (n_file2_conc .gt. n_free_iounits) call error_handler( & 
                                "int_send_2cob3c_shutdown: concating files failed" )
                           n_file2 = n_file2 - n_concat
                        enddo
                        endif

                        n_file2_new = i_file2_new

                        ! opening the qu_file tapehandles
                        do i_file2 = 1, n_file2
                           call trace_message()
                           write(ch_i_file2,'(i5)') i_file2        !!!!!!!!!!!!!!!
                           ch_i_file2 = adjustl(ch_i_file2)
                           qfh(i_file2)%total_length = quadrupelfile_length(i_file1,i_file2,i_ir)

                         if(options_directaccess_integrals() ) then
                           qfh_da%total_length = quadrupelfile_length(i_file1,i_file2,i_ir)
                           call readwriteblocked_startread(  trim(tmpfile("da")), &
                                qfh(i_file2)%th, blocklength=quadrupels_reclength, &
                                total_length=qfh(i_file2)%total_length,rec=da_rec_no, unit=da_unit )
!                            print*,'rea dim',darec_no(i_ir,i_file1,i_file2)
                            call read_dimensions(qfh(i_file2),rec=darec_no(i_ir,i_file1,i_file2) )
                         else 
                           call readwriteblocked_startread( &
                                trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                                trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")), &
                                qfh(i_file2)%th, blocklength=blocklength, &
                                total_length=qfh(i_file2)%total_length )
                            call read_dimensions(qfh(i_file2))
                         endif
                        enddo

                        do i_file2_new = i_file1+1, n_file2_new
                           write(ch_i_file2,'(i5)') i_file2_new                !!!!!!!!!!!!
                           ch_i_file2 = adjustl(ch_i_file2)
                           n_file2 = n_file2 + 1
                           if(options_directaccess_integrals()) then
                           print*,'startread in i_file2_new cycle'
                           call readwriteblocked_startread( trim(tmpfile("da")), &
                                qfh(n_file2)%th, blocklength=quadrupels_reclength, &
                                total_length=qfh(i_file2)%total_length,rec=darec_no(i_ir,i_file1,i_file2) )
                           call read_dimensions(qfh(n_file2),rec=darec_no(i_ir,i_file1,i_file2))
                           else
                           call readwriteblocked_startread( &
                                trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                                trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")), &
                                qfh(n_file2)%th, blocklength=blocklength, &
                                total_length=qfh(i_file2)%total_length )
                           call read_dimensions(qfh(n_file2))
                           endif
                        enddo

                        ! sum up number of records
                        do i_file2 = 1, n_file2
                           n_tot_records = n_tot_records + sum(qfh(i_file2)%n_floats_if_c)
                        enddo
                        n_if_c=qfh(1)%n_if_c
                        rel=.false.

                        if(integralpar_relativistic) then
                           do i_file2 = 1, n_file2
                              n_tot_records_rel = n_tot_records_rel + sum(qfh(i_file2)%n_floats_if_rel)
                           enddo
                           n_if_c=qfh(1)%n_if_rel
                           rel=.true.
                        end if

                        ! do the rewriting
                        if ( comm_i_am_master()) then
                           if (integralpar_2cob_kin) then
                              do i_if_c = 1, n_if_c
                                 do i_file2 = 1, n_file2
                                if(options_directaccess_integrals()) then
!                                print*,'rew kin',darec_no(i_ir,i_file1,i_file2)
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_kin,1,rel, &
                                                                rec=darec_no(i_ir,i_file1,i_file2))
                                else
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_kin,1,rel)
                                  endif
                                 enddo
                              enddo
                           endif

                           if (integralpar_2cob_nuc) then
                              do i_if_c = 1,n_if_c
                                 do i_file2 = 1, n_file2
                                    if(options_directaccess_integrals()) then
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_nuc,1,rel &
                                                             ,rec=darec_no(i_ir,i_file1,i_file2))
                                    else
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_nuc,1,rel)
                                    endif
                                 enddo
                              enddo
                           endif
                           if (integralpar_2cob_pvsp) then
                              do i_if_c = 1,n_if_c
                                 do i_file2 = 1, n_file2
                                    if(options_directaccess_integrals()) then
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_pvsp,1,rel &
                                                              ,rec=darec_no(i_ir,i_file1,i_file2))
                                    else
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_pvsp,1,rel)
                                    endif
                                 enddo
                              enddo
                           endif
!!!!!!!!                    if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                        do i_if_c = 1,n_if_c
!!!!!!!!                         do i_file2 = 1, n_file2
!!!!!!!!                          if(options_directaccess_integrals()) then
!!!!!!!!                           call rewrite_blocked_integralfile(qfh(i_file2),th_nuc_pseudo,1,rel &
!!!!!!!!                                                           ,rec=darec_no(i_ir,i_file1,i_file2))
!!!!!!!!                          else
!!!!!!!!                           call rewrite_blocked_integralfile(qfh(i_file2),th_nuc_pseudo,1,rel)
!!!!!!!!                          endif
!!!!!!!!                         enddo
!!!!!!!!                       enddo
!!!!!!!!                    endif
                         endif! comm_i_am_master()

                        if ( integralpar_2cob_ol ) then
                        
                           if(.not.integralpar_relativistic) then
                              do i_if_c = 1, qfh(1)%n_if_c
                                 do i_file2 = 1, n_file2
                                 if(options_directaccess_integrals()) then
                                    call rewrite_integralfile(qfh(i_file2),unit_ol,1 &
                                                 ,rec=darec_no(i_ir,i_file1,i_file2) )
                                 else
                                    call rewrite_integralfile(qfh(i_file2),unit_ol,1)
                                 endif
                                 enddo
                              enddo
                           else
                              do i_if_c = 1, n_if_c
                                 do i_file2 = 1, n_file2               
                                if(options_directaccess_integrals()) then
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_ol,1,rel,&
                                                rec=darec_no(i_ir,i_file1,i_file2))
                                else
                                    call rewrite_blocked_integralfile(qfh(i_file2),th_ol,1,rel)
                                endif
                                 end do
                              end do
                           end if
                        endif
                        if ( integralpar_3c_xc ) then
                           do i_if_c = 1, qfh(1)%n_if_c
                              do i_file2 = 1, n_file2
                             if(options_directaccess_integrals()) then
                                 call rewrite_blocked_integralfile( qfh(i_file2),th_xc, &
                                      borders_xc(3,my_hostindex),rec=darec_no(i_ir,i_file1,i_file2) )
                             else
                                 call rewrite_blocked_integralfile( qfh(i_file2),th_xc, &
                                      borders_xc(3,my_hostindex) )
                             endif
                              enddo
                           enddo
                        endif

                        if ( integralpar_3c_co ) then
                           do i_if_c = 1, qfh(1)%n_if_c
                              do i_file2 = 1, n_file2
                                if(options_directaccess_integrals()) then
                                 call rewrite_blocked_integralfile( qfh(i_file2),th_co, &
                                      borders_ch(3,my_hostindex),rec=darec_no(i_ir,i_file1,i_file2) )
                                else
                                 call rewrite_blocked_integralfile( qfh(i_file2),th_co, &
                                      borders_ch(3,my_hostindex) )
                                endif
                              enddo
                           enddo
                        endif

                        if ( integralpar_2cob_potential ) then
                           do i_if_c = 1, qfh(1)%n_if_c
                              do i_file2 = 1, n_file2
                                 if(options_directaccess_integrals()) then
                                    call rewrite_blocked_integralfile( qfh(i_file2),th_poten, &
                                         borders_poten(3,my_hostindex),rec=darec_no(i_ir,i_file1,i_file2) )
                                 else
                                    call rewrite_blocked_integralfile( qfh(i_file2),th_poten, &
                                         borders_poten(3,my_hostindex) )
                                 end if
                              enddo
                           enddo
                        endif
                        if ( integralpar_2cob_field ) then
                           do i_if_c = 1, qfh(1)%n_if_c
                              do i_file2 = 1, n_file2
                                 if(options_directaccess_integrals()) then
                                    call rewrite_blocked_integralfile( qfh(i_file2),th_field, &
                                         borders_field(3,my_hostindex),rec=darec_no(i_ir,i_file1,i_file2) )
                                 else
                                    call rewrite_blocked_integralfile( qfh(i_file2),th_field, &
                                         borders_field(3,my_hostindex) )
                                 endif
                              enddo
                           enddo
                        endif
!        print*,         ' closing the qu_file tapehandles with stopread'
                        do i_file2 = 1, n_file2
                           if(options_directaccess_integrals()) then
                           call readwriteblocked_stopread(qfh(i_file2)%th, &
                                         rec=darec_no(i_ir,i_file1,i_file2),status='KEEP')
                           else
                           call readwriteblocked_stopread(qfh(i_file2)%th,status=DELETE)
                           endif
                           call free_dimensions(qfh(i_file2))
                        enddo
                     enddo
                  enddo! loop over irreps
                  if(options_directaccess_integrals()) then
!                   print*,'shutdown returnclose da_unit',da_unit
!                    call returnclose_iounit(da_unit,status='DELETE')
                    call returnclose_iounit(da_unit,status='KEEP')
                  endif
                  deallocate( qfh, stat=status )
                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: deallocate of qfh failed" )
                     deallocate( blocksizes, rel_flag, stat=status )
                     if ( status .ne. 0 ) call error_handler( & 
                          "int_send_2cob3c_shutdown: deallocate of blocksizes failed" )

                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_shutdown: total number of records:",inte=n_tot_records)

                  ! closing of files
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_shutdown: closing output files")
                  if ( comm_i_am_master()) then
                     if (integralpar_2cob_kin) &
                          call readwriteblocked_stopwrite(th_kin)
                     if (integralpar_2cob_nuc) &
                          call readwriteblocked_stopwrite(th_nuc)
                     if (integralpar_2cob_pvsp) &
                          call readwriteblocked_stopwrite(th_pvsp)
                     !if (integralpar_relativistic.and. integralpar_2cob_ol ) &
                     !     call readwriteblocked_stopwrite(th_ol)
!!!!!!!!             if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                call readwriteblocked_stopwrite(th_nuc_pseudo)
!!!!!!!!             endif! 
                  endif! comm_i_am_master()

                  if (integralpar_relativistic.and. integralpar_2cob_ol ) &
                       call readwriteblocked_stopwrite(th_ol)
                  if (.not.integralpar_relativistic.and. integralpar_2cob_ol ) &
                       call returnclose_iounit(unit_ol,status='KEEP')
                  if ( integralpar_3c_xc ) &
                       call readwriteblocked_stopwrite(th_xc)
                  if ( integralpar_3c_co ) &
                       call readwriteblocked_stopwrite(th_co)
                  if ( integralpar_2cob_potential ) &
                       call readwriteblocked_stopwrite(th_poten)
                  if ( integralpar_2cob_field ) &
                       call readwriteblocked_stopwrite(th_field)
                  if(integralpar_relativistic) then ! do some preparations for rewrite
                     n_irrep=symmetry_data_n_irreps()
                     allocate(th_kin_arr(n_irrep),th_nuc_arr(n_irrep),th_ol_arr(n_irrep),&
                          th_pvsp_arr(n_irrep),dimensions(n_irrep),stat=status)
                     if(status/=0) call error_handler(&
                          "int_send_2cob3c_shutdown: allocation unit_kin failed")
                     ! now calculate dimensions of uncontracted matrixes for every irrep
!!!!!!!!             if(pseudopot_present) then
!!!!!!!!                allocate(th_nuc_pseudo_arr(n_irrep),stat=status)
!!!!!!!!                if(status/=0) call error_handler(&
!!!!!!!!                     "int_send_2cob3c_shutdown: allocation unit_nuc_pseudo failed")
!!!!!!!!             endif
                        
          
                     do i_ir = 1, n_irrep
                        dimensions(i_ir)=0
                        do i_unique=1,n_unique_atoms
                           do i_l=0,unique_atoms(i_unique)%lmax_ob
                              dimensions(i_ir)=dimensions(i_ir)+unique_atoms(i_unique)%&
                                   symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                                   unique_atoms(i_unique)%l_ob(i_l)%n_exponents
                           end do
                        end do
                        dimensions(i_ir)=dimensions(i_ir)*(dimensions(i_ir)+1)/2
                     end do
                  endif

                  !shutdown:  combining kin and nuc files
                  if ( comm_i_am_master() .and. &
                       integralpar_2cob_kin .and. integralpar_2cob_nuc ) then

                     if ( output_int_loops ) call write_to_output_units( &
                          "int_send_2cob3c_shutdown: combining kin and nuc files")
                     if (operations_integral) then
                        if(integralpar_relativistic) then
                           ! for every irrep we use a seperate file
                           do i_ir=1,n_irrep
                              write(inter,'(i3)') i_ir
                              inter=adjustl(inter)
!                              print*,'startwrite kin.dat 1 '
                              call readwriteblocked_startwrite& ! kin.dat
                                   (trim(tmpfile("kin.dat"//trim(inter))), th_kin_arr(i_ir))
                              call readwriteblocked_startwrite&  ! nuc.dat
                                   (trim(tmpfile("nuc.dat"//trim(inter))), th_nuc_arr(i_ir))
                              call readwriteblocked_startwrite& ! overlap_rel.dat
                                   (trim(tmpfile("overlap_rel.dat"//trim(inter))), th_ol_arr(i_ir))
                           end do
                        else ! not integralpar_relativistic
                           unit_kin_nuc = openget_iounit(trim(tmpfile("ham_kin_nuc.dat")), &
                                status='REPLACE', form='UNFORMATTED', action='WRITE')
                        end if
                        if (integralpar_2cob_pvsp) then
                           do i_ir=1,n_irrep
                              write(inter,'(i3)') i_ir
                              inter=adjustl(inter)
                              call readwriteblocked_startwrite& !pvsp.dat
                                   (trim(tmpfile("pvsp.dat"//trim(inter))), th_pvsp_arr(i_ir))
                           end do
                        end if

                     else
                        call error_handler("int_send_2cob3c_shutdown: wrong operations")
                     endif

                     call readwriteblocked_startread(trim(tmpfile("kin.dat")), th_kin)
                     call readwriteblocked_startread(trim(tmpfile("nuc.dat")), th_nuc)
                     if (integralpar_2cob_pvsp) then
                        call readwriteblocked_startread(trim(tmpfile("pvsp.dat")), th_pvsp)
                        call rewrite_one_integralfile(th_pvsp,th_pvsp_arr)
                        call readwriteblocked_stopread(th_pvsp,status=DELETE)
                     end if

                     if (integralpar_2cob_ol.and.integralpar_relativistic) then
                        call readwriteblocked_startread(trim(tmpfile("overlap.dat")), th_ol)
                        call rewrite_one_integralfile(th_ol,th_ol_arr)
                        call readwriteblocked_stopread(th_ol,status=DELETE)
                     end if

                     if(integralpar_relativistic) then
                        call rewrite_one_integralfile(th_kin,th_kin_arr)
                        call rewrite_one_integralfile(th_nuc,th_nuc_arr)
!!!!!!!!                if(pseudopot_present)then
!!!!!!!!                   call rewrite_one_integralfile(th_nuc_pseudo, th_nuc_pseudo_arr)
!!!!!!!!                endif
                     else
                        call combine_two_integralfiles(th_kin,th_nuc,unit_kin_nuc)
                        call returnclose_iounit(unit_kin_nuc,status='KEEP')
                     end if
                     call readwriteblocked_stopread(th_kin,status=DELETE)
                     call readwriteblocked_stopread(th_nuc,status=DELETE)
!!!!!!!!             if(pseudopot_present.and.integralpar_relativistic.and.integralpar_2cob_nuc)then
!!!!!!!!SR8000
!!!!!!!!                call readwriteblocked_stopread(th_nuc_pseudo,status='KEEP')
!!!!!
!!!!!!!!                call readwriteblocked_stopread(th_nuc_pseudo,status='DELETE')
!!!!!!
!!!!!!!!             endif
                  endif! comm_i_am_master() .and. ...
!                       files combined

                  if(integralpar_relativistic) then                        
                     deallocate(th_kin_arr,th_nuc_arr,th_ol_arr,&
                          th_pvsp_arr,dimensions,stat=status)
                     if(status/=0) call error_handler(&
                          "int_send_2cob3c_shutdown: deallocating unit_kin_arr failed")
!!!!!!!!             if(pseudopot_present) then
!!!!!!!!                deallocate(th_nuc_pseudo_arr,stat=status)
!!!!!!!!                if(status/=0) call error_handler(&
!!!!!!!!                     "int_send_2cob3c_shutdown: deallocating unit th_nuc_pseudo_arr  failed")
!!!!!!!!             endif
                  endif

               ! deallocation of help arrays
               deallocate( ua_fileindex, stat=status )
               if ( status .ne. 0 ) call error_handler( & 
                    "int_send_2cob3c_shutdown: deallocate of ua_fileindex failed" )
               deallocate( quadrupelfile_length, stat=status )
               if ( status .ne. 0 ) call error_handler( & 
                    "int_send_2cob3c_shutdown: deallocate of quadrupelfile_length failed" )
               
               call stop_timer(timer_int_rewrite_2cob3c(integralpar_i_int_part))
               if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

           else! .not. options_integrals_on_file()

               ! the integrals are already in their final location

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
                                         "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral l2 failed" )
                                 enddo
                              enddo
                           enddo
                           deallocate( ind1_exp1_ua2, stat=status )
                           if ( status .ne. 0 ) call error_handler( & 
                                "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral ind1_exp1_ua2 failed" )
                        enddo
                        deallocate( metaindex_of_first_integral(i_ir)%ua1(i_ua1)%l1, stat=status )
                        if ( status .ne. 0 ) call error_handler( & 
                             "int_send_2cob3c_shutdown: allocate of metaindex_of_first_integral l1 failed" )
                     enddo
                     deallocate( metaindex_of_first_integral(i_ir)%ua1, &
                          stat=status )
                     if ( status .ne. 0 ) call error_handler( & 
                          "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral ua1 failed" )
                  enddo
                  deallocate( metaindex_of_first_integral, &
                       stat=status )
                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral failed" )

               if(integralpar_relativistic) then
                  do i_ir = 1,symmetry_data_n_irreps()
                     do i_ua1 = 1, N_unique_atoms
                        ua1 => unique_atoms(i_ua1)
                        do i_l1 = 0, ua1%lmax_ob
                           uab1 => ua1%l_ob(i_l1)
                           uap1 => ua1%symadapt_partner(i_ir,i_l1)
                           ind1_exp1_ua2 => metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1(i_l1)%ind1_exp1_ua2
                           do i_ind1 = 1, uap1%N_independent_fcts
                              do i_exp1 = 1, uab1%N_exponents
                                 do i_ua2 = 1, i_ua1
                                    deallocate( ind1_exp1_ua2(i_ind1,i_exp1,i_ua2)%l2, &
                                         stat=status )
                                    if ( status .ne. 0 ) call error_handler( & 
                                         "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral_rel l2 failed" )
                                 enddo
                              enddo
                           enddo
                           deallocate( ind1_exp1_ua2, stat=status )
                           if ( status .ne. 0 ) call error_handler( & 
                                "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral_rel ind1_exp1_ua2 failed" )
                        enddo
                        deallocate( metaindex_of_first_integral_rel(i_ir)%ua1(i_ua1)%l1, stat=status )
                        if ( status .ne. 0 ) call error_handler( & 
                             "int_send_2cob3c_shutdown: allocate of metaindex_of_first_integral_rel l1 failed" )
                     enddo
                     deallocate( metaindex_of_first_integral_rel(i_ir)%ua1, &
                          stat=status )
                     if ( status .ne. 0 ) call error_handler( & 
                          "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral_rel ua1 failed" )
                  enddo
                  deallocate( metaindex_of_first_integral_rel, &
                       stat=status )
                  if ( status .ne. 0 ) call error_handler( & 
                       "int_send_2cob3c_shutdown: deallocate of metaindex_of_first_integral_rel failed" )
               end if! integralpar_relativistic


               
            endif integrals_on_file!/else

            ! deallocation of help arrays
            deallocate( borders_ch, borders_xc, borders_poten, borders_field, stat=status )
            
            if ( status .ne. 0 ) call error_handler( & 
                    "int_send_2cob3c_shutdown: deallocate of borders failed" )

            if(options_directaccess_integrals() ) deallocate(darec_no,stat=status )
        
            if ( status .ne. 0 ) call error_handler( &
                    "int_send_2cob3c_shutdown: deallocate darec_no failed" )

            if ( output_int_loops ) call write_to_output_units("int_send_2cob3c_shutdown: done")

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
                   "int_send_2cob3c_shutdown: allocating at write_dimensions() failed")
              real_n_if_c(1) = real(qfh%n_if_c,r8_kind)
              real_n_floats_if_c = real(qfh%n_floats_if_c,r8_kind)
              call readwriteblocked_write(real_n_if_c, qfh%th)
              call readwriteblocked_write(real_n_floats_if_c,qfh%th)
              deallocate( real_n_floats_if_c, stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_shutdown: deallocating at write_dimensions() failed")
              if(integralpar_relativistic) then
                 allocate( real_n_floats_if_c(qfh%n_if_rel), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: allocating at write_dimensions() failed")
                 real_n_if_c(1) = real(qfh%n_if_c,r8_kind)
                 real_n_floats_if_c = real(qfh%n_floats_if_rel,r8_kind)
                 call readwriteblocked_write(real_n_if_c, qfh%th)
                 call readwriteblocked_write(real_n_floats_if_c,qfh%th)
                 deallocate( real_n_floats_if_c, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: deallocating at write_dimensions() failed")
              end if
            end subroutine write_dimensions


            subroutine read_dimensions(qfh,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              type(quadrupel_file_handle_type), intent(inout)  :: qfh
              integer(kind=i4_kind), intent(inout),optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: status
              real(kind=r8_kind) :: real_n_if_c(1)
              real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
              !------------ Executable code ----------------------------
              if(present(rec)) then
               call readwriteblocked_read(real_n_if_c,qfh%th,rec=rec)
              else
               call readwriteblocked_read(real_n_if_c,qfh%th)
              endif

              qfh%n_if_c = nint(real_n_if_c(1),i4_kind)
              allocate( real_n_floats_if_c(qfh%n_if_c), &
                   qfh%n_floats_if_c(qfh%n_if_c), stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_shutdown: allocating at read_dimensions() failed")

              if(present(rec)) then
               call readwriteblocked_read(real_n_floats_if_c,qfh%th,rec=rec)
              else
               call readwriteblocked_read(real_n_floats_if_c,qfh%th)
              endif

              qfh%n_floats_if_c = nint(real_n_floats_if_c,i4_kind)
              deallocate( real_n_floats_if_c, stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_shutdown: deallocating at read_dimensions() failed")

              if(integralpar_relativistic) then
              if(present(rec)) then
                 call readwriteblocked_read(real_n_if_c,qfh%th,rec=rec)
              else
                 call readwriteblocked_read(real_n_if_c,qfh%th)
              endif 

                 qfh%n_if_rel = nint(real_n_if_c(1),i4_kind)
                 allocate( real_n_floats_if_c(qfh%n_if_rel), &
                      qfh%n_floats_if_rel(qfh%n_if_rel), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: allocating at read_dimensions() failed")
              if(present(rec)) then
                 call readwriteblocked_read(real_n_floats_if_c,qfh%th,rec=rec)
              else        
                 call readwriteblocked_read(real_n_floats_if_c,qfh%th)
              endif
                 qfh%n_floats_if_rel = nint(real_n_floats_if_c,i4_kind)
                 deallocate( real_n_floats_if_c, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: deallocating at read_dimensions() failed")
              end if

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
                   "int_send_2cob3c_shutdown: allocating at sum_dimensions() failed")
              qfh_sum%n_floats_if_c = 0
              do i_f = 1, size(qfh_in,1)
                 qfh_sum%n_floats_if_c = qfh_sum%n_floats_if_c + qfh(i_f)%n_floats_if_c
              enddo
              if(integralpar_relativistic) then
                 qfh_sum%n_if_rel = qfh_in(1)%n_if_rel
                 allocate( qfh_sum%n_floats_if_rel(qfh_sum%n_if_rel), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: allocating at sum_dimensions() failed")
                 qfh_sum%n_floats_if_rel = 0
                 do i_f = 1, size(qfh_in,1)
                    qfh_sum%n_floats_if_rel = qfh_sum%n_floats_if_rel + qfh(i_f)%n_floats_if_rel
                 enddo
              end if

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
                   "send_2cob3c_send: deallocating at free_dimensions() failed")
              if(integralpar_relativistic) then
                 deallocate( qfh%n_floats_if_rel,stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: deallocating at free_dimensions() failed")
              end if
            end subroutine free_dimensions


            subroutine concat_files(i_file2_from,i_file2_to,i_file2_new,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              integer(kind=i4_kind), intent(in) :: i_file2_from,i_file2_to,i_file2_new
              integer(kind=i4_kind), intent(inout),optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_integral,i_file2,i_f,n_f, n_if_c
              type(quadrupel_file_handle_type) :: qfh_new
              !------------ Executable code ----------------------------
              ! opening the qu_file tapehandles to read from
              i_f = 0
              do i_file2 = i_file2_from, i_file2_to
                 call trace_message()
                 i_f = i_f + 1
                 write(ch_i_file2,'(i5)') i_file2        !!!!!!!!!!!!
                 ch_i_file2 = adjustl(ch_i_file2)
                 qfh(i_f)%total_length = quadrupelfile_length(i_file1,i_file2,i_ir)
                if(present(rec)) then
                print*,'concating files is not requiered for direct_access integrals'
                else
                 call readwriteblocked_startread( &
                      trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                      trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")), &
                      qfh(i_f)%th, blocklength=blocklength, &
                      total_length=qfh(i_f)%total_length )
                endif
                 call read_dimensions(qfh(i_f))
              enddo
              n_f = i_f 
              ! opening the tapehandle to write to
              write(ch_i_file2,'(i5)') i_file2_new        !!!!!!!!!!!!!!!!!!
              ch_i_file2 = adjustl(ch_i_file2)
              call readwriteblocked_startwrite( & ! quadrupel - concat_files
                   trim(tmpfile("quadrupel_" // trim(ch_i_ir) // "_" // &
                   trim(ch_i_file1) // "_" // trim(ch_i_file2) // ".dat")), &
                   qfh_new%th, blocklength=blocklength )
              ! calculate and write size information for new file
              call sum_dimensions(qfh(1:n_f),qfh_new)
              call write_dimensions(qfh_new)
              ! do the rewriting
              do i_integral = 1, n_integrals
                 if(integralpar_relativistic) then
                    n_if_c=qfh_new%n_if_rel
                 else
                    n_if_c=qfh_new%n_if_c
                 end if
                 do i_if_c = 1, n_if_c
                    do i_f = 1, n_f
                       call rewrite_blocked_integralfile(qfh(i_f), &
                            qfh_new%th,blocksizes(i_integral),rel_dummy=rel_flag(i_integral))
                    enddo
                 enddo
              enddo
              ! closing the qu_file tapehandles and removing the files
              do i_f = 1, n_f
                 call readwriteblocked_stopread(qfh(i_f)%th,status=DELETE)
                 call free_dimensions(qfh(i_f))
              enddo
              ! closing tapehandle that has been written to
              call readwriteblocked_stopwrite(qfh_new%th,total_length=qfh_new%total_length)
              call free_dimensions(qfh_new)
            end subroutine concat_files


            subroutine rewrite_blocked_integralfile(qfh_in,th_out,bordersize, rel_dummy,test,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              type(quadrupel_file_handle_type), intent(inout) :: qfh_in
              type(readwriteblocked_tapehandle), intent(inout) :: th_out
              integer(kind=i4_kind), intent(in) :: bordersize
              logical,optional :: rel_dummy,test
              integer(kind=i4_kind), intent(inout), optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: n_floats,status
              real(kind=r8_kind), allocatable :: buffer(:)
              logical :: rel
              !------------ Executable code ----------------------------
              rel=.false.
              if(present(rel_dummy)) rel=rel_dummy
              if(.not.rel) then
                 n_floats = qfh_in%n_floats_if_c(i_if_c) * bordersize
              else
                 n_floats = qfh_in%n_floats_if_rel(i_if_c) * bordersize
              end if
              allocate(buffer(n_floats), stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_shutdown: rewrite_blocked_integralfile: allocating failed")
                if(present(test)) print*,'rewrite_blocked_integralfile: read n_floats ',n_floats
              if(present(rec)) then
               call readwriteblocked_read(buffer,qfh_in%th,rec=rec)
              else
               call readwriteblocked_read(buffer,qfh_in%th)
              endif
                if(present(test)) then
                 print*,buffer
                endif
                  call readwriteblocked_write(buffer,th_out)
              deallocate(buffer, stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_shutdown: rewrite_blocked_integralfile: deallocating failed")
            end subroutine rewrite_blocked_integralfile


            subroutine rewrite_integralfile(qfh_in,unit_out,bordersize,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              type(quadrupel_file_handle_type), intent(inout) :: qfh_in
              integer(kind=i4_kind), intent(in) :: unit_out,bordersize
              integer(kind=i4_kind),optional, intent(inout) :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: n_floats,status,i_float
              real(kind=r8_kind), allocatable :: buffer(:)
              !------------ Executable code ----------------------------
              if(.not.integralpar_relativistic) then
                 n_floats = qfh_in%n_floats_if_c(i_if_c) * bordersize
              else
                 n_floats = qfh_in%n_floats_if_rel(i_if_c) * bordersize
              end if
              allocate(buffer(n_floats), stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_shutdown: rewrite_integralfile: allocating failed")
              if(present(rec)) then
              call readwriteblocked_read(buffer,qfh_in%th,rec=rec)
              else
              call readwriteblocked_read(buffer,qfh_in%th)
              endif
              do i_float = 1, n_floats
                 write(unit_out,ERR=99) buffer(i_float)
              enddo
              deallocate(buffer, stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_shutdown: rewrite_integralfile: deallocating failed")
              return
        99    call error_handler("send_2cob3c_shutdown: rewrite_integralfile: writing failed")
            end subroutine rewrite_integralfile


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
                   "send_2cob3c_shutdown: combine_two_integralfiles: allocating failed")
              if(.not.integralpar_relativistic) then
                 n_rec_left = n_tot_records
              else
                 n_rec_left = n_tot_records_rel
              end if
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

            subroutine rewrite_one_integralfile(th_in,th_out,test)
              implicit none
              ! note: in the relativistic case we have one file for every irrep
              !------------ Declaration of formal parameters -----------
              type(readwriteblocked_tapehandle), intent(inout) :: th_in
              type(readwriteblocked_tapehandle), intent(inout) :: th_out(n_irrep)
              logical, optional :: test
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
                if(present(test).and.i_ir.eq.1) then
        !? print*,buffer(1:dimensions(i_ir))
                endif
                 call readwriteblocked_write(buffer(1:dimensions(i_ir)),th_out(i_ir))         
                 call readwriteblocked_stopwrite(th_out(i_ir),total_length=length)
                 deallocate(buffer, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_shutdown: rewrite_one_integralfile: deallocating failed")
              end do
            end subroutine rewrite_one_integralfile


            subroutine trace_message()
              n_files_left = n_files_left - 1

              !
              ! Trace file is not always availabel for slaves:
              !
              if ( comm_i_am_master() ) then
                 if (mod(n_files_left,50) .eq. 0) &
                      call write_to_trace_unit("int_send_2cob3c_shutdown: &
                      &number of files left: ",inte=n_files_left)
              endif
            end subroutine trace_message


          end subroutine int_send_2cob3c_shutdown
          !*************************************************************
         

          !*************************************************************
          subroutine int_send_2cob3c_send
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
            !    everything is directly written to files.
            !  called by: integral_calc_quad_2cob3c
            !** End of interface ***************************************
            implicit none
            !------------ Declaration of local variables -----------------
            ! Nothing to be done for gradients:
            if( integralpar_gradients ) RETURN

            !***********************************************************
            ! Some two center integrals (kin, nuc, ...)
            ! are send directly to memory, since the integrals_on_file
            ! is sentenced to death for these hamiltonian terms
            !***********************************************************

            !***********************************************************
            ! For the other integrals the usual path is followed ...
            !***********************************************************
            if ( options_integrals_on_file() ) then
               call int_send_2cob3c_send_file()
            else
               call int_send_2cob3c_send_mem()
            endif
            ! MOVED from inside of above two subs:
            n_missing_quadrupels = n_missing_quadrupels - 1
          end subroutine int_send_2cob3c_send
          !*************************************************************
         
!DG
          !*************************************************************
          subroutine int_send_2cob3c_send_file
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
            use int_data_2cob3c_module
            use options_module, only: options_directaccess_integrals,quadrupels_reclength
            implicit none
            !------------ Declaration of local variables -----------------
            integer(kind=i4_kind) :: i_host, i_ir, n_ir, n_if_c_max, n_if_rel_max
            type(readwriteblocked_tapehandle),target :: th
            type(readwriteblocked_tapehandle),pointer :: th_pointer
            integer(kind=i4_kind), allocatable :: &
                 n_floats_if_c(:,:), n_floats_if_c_max(:), n_if_c(:),&
                 n_floats_if_rel(:,:), n_floats_if_rel_max(:), n_if_rel(:)    
            ! n_floats_if_c(i_if_c,i_ir), n_floats_if_c_max(i_ir), n_if_c(i_ir)
            ! n_floats_if_c(i_if_rel,i_ir), n_floats_if_rel_max(i_ir), n_if_rel(i_ir)
            ! write some informations about quadrupel to file
            !------------ Executable code --------------------------------

            n_ir = symmetry_data_n_irreps()

            call calc_dimensions()

            hosts: do i_host = first_host, last_host
            

               parallel_filesystem: if ( filesystem_is_parallel ) then

                  call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))
                  
                  do i_ir = 1, n_ir

                        if(options_directaccess_integrals())  then
                          call readwriteblocked_startwrite(trim(tmpfile("da")), th_da, &
                                                  blocklength=quadrupels_reclength,rec=da_rec_no)
                            thdaf_pointer=> th_da !used in subr
                            darec_no(i_ir,index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                         index_ua_l(quadrupel%ua2,quadrupel%l2))=da_rec_no
                        print*,' filesystem_is_parallel readwriteblocked_startwrite access_direct',&
                                index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                index_ua_l(quadrupel%ua2,quadrupel%l2),da_rec_no,comm_myindex()
                            call write_dimensions(th_da,i_ir,rec=da_rec_no)
                else

                     call readwriteblocked_startwrite( & !quadrupel_filename parallel_filesys
                          trim(filename_tmpdir(i_host)) // "/" // trim(quadrupel_filename(i_ir,quadrupel)), &
                          th, blocklength=blocklength )
                     th_pointer=> th
                     call write_dimensions(th,i_ir)
                endif

                     if ( i_host .eq. 1) then
                        if ( integralpar_2cob_kin ) then
                           if ( output_int_deeploops ) call write_to_output_units( &
                                "send_2cob3c_send: writing kin")
                         if(options_directaccess_integrals())  then
                           call write_integrals_2c_file(symadapt_int_2cob_kin(i_ir)%int,rec=da_rec_no)
                         else
                           call write_integrals_2c_file(symadapt_int_2cob_kin(i_ir)%int)
                         endif
                        endif
                        if ( integralpar_2cob_nuc ) then
                           if ( output_int_deeploops ) call write_to_output_units( &
                                "send_2cob3c_send: writing nuc")
                         if(options_directaccess_integrals())  then
                           call write_integrals_2c_file(symadapt_int_2cob_nuc(i_ir)%int,rec=da_rec_no)
                         else
                           call write_integrals_2c_file(symadapt_int_2cob_nuc(i_ir)%int)
                         endif
                        endif
                        if ( integralpar_2cob_pvsp ) then
                           if ( output_int_deeploops ) call write_to_output_units( &
                                "send_2cob3c_send: writing pvsp")
                         if(options_directaccess_integrals())  then
                           call write_integrals_2c_file(symadapt_int_2cob_pvsp(i_ir)%int,rec=da_rec_no)
                         else
                           call write_integrals_2c_file(symadapt_int_2cob_pvsp(i_ir)%int)
                         endif
                        endif
!!!!!!!!                if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                   if ( output_int_deeploops ) call write_to_output_units( &
!!!!!!!!                        "send_2cob3c_send: writing pseudo")

!!!!!!!!                 if(options_directaccess_integrals())  then
!!!!!!!!                   call write_integrals_2c_file(symadapt_int_2cob_nuc_pseudo(i_ir)%int,rec=da_rec_no)
!!!!!!!!                 else
!!!!!!!!                   call write_integrals_2c_file(symadapt_int_2cob_nuc_pseudo(i_ir)%int)
!!!!!!!!                 endif

!!!!!!!!                endif
                     endif! i_host .eq. 1

                     if ( integralpar_2cob_ol ) then
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_send: writing overlap")
                         if(options_directaccess_integrals())  then
                          call write_integrals_2c_file(symadapt_int_2cob_ol(i_ir)%int,rec=da_rec_no)
                         else
                          call write_integrals_2c_file(symadapt_int_2cob_ol(i_ir)%int)
                         endif
                     endif
                     if ( integralpar_3c_xc ) then
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_send: writing coulomb")
                      if(options_directaccess_integrals())  then
                        call write_integrals_3c_file( symadapt_int_3c_xc(i_ir)%int, &
                             borders_xc(:,i_host),rec=da_rec_no )
                      else
                        call write_integrals_3c_file( symadapt_int_3c_xc(i_ir)%int, &
                             borders_xc(:,i_host) )
                      endif
                     endif
                     if ( integralpar_3c_co ) then
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_send: writing exchange")
                      if(options_directaccess_integrals())  then
                        call write_integrals_3c_file( symadapt_int_3c_co(i_ir)%int, &
                             borders_ch(:,i_host),rec=da_rec_no )
                      else
                        call write_integrals_3c_file( symadapt_int_3c_co(i_ir)%int, &
                             borders_ch(:,i_host) )
                      endif
                     endif
                     if ( integralpar_2cob_potential ) then
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_send: writing potential 1")
                         if(options_directaccess_integrals())  then
                            call write_integrals_3c_file( symadapt_int_2cob_poten(i_ir)%int, &
                                 borders_poten(:,i_host),rec=da_rec_no )
                         else
                            call write_integrals_3c_file( symadapt_int_2cob_poten(i_ir)%int, &
                                 borders_poten(:,i_host) )
                         endif
                     endif
                     if ( integralpar_2cob_field ) then
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_send: writing field")
                         if(options_directaccess_integrals())  then
                            call write_integrals_3c_file( symadapt_int_2cob_field(i_ir)%int, &
                                 borders_field(:,i_host),rec=da_rec_no )
                         else
                            call write_integrals_3c_file( symadapt_int_2cob_field(i_ir)%int, &
                                 borders_field(:,i_host) )
                         endif
                     endif
                         if(options_directaccess_integrals())  then
                            call readwriteblocked_stopwrite(th_da, total_length= &
                                quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ),rec=da_rec_no)
                         else
                     call readwriteblocked_stopwrite(th, total_length= &
                          quadrupelfile_length_allhosts( &
                          index_ua_l(quadrupel%ua1,quadrupel%l1), &
                          index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir, i_host ) )
                         endif
                  enddo
                  call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))

               else parallel_filesystem

                  myhost: if ( i_host .eq. my_hostindex) then
!!!!!!!!!!!!!!!!     DPRINT MyID, ' pseudopot_present = ', pseudopot_present
                     call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))

                     irreps: do i_ir = 1, n_ir

                         if(options_directaccess_integrals())  then
                          call readwriteblocked_startwrite(trim(tmpfile("da")), th_da, &
                                                  blocklength=quadrupels_reclength,rec=da_rec_no)
                          thdaf_pointer=> th_da
                            darec_no(i_ir,index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                         index_ua_l(quadrupel%ua2,quadrupel%l2))=da_rec_no
!                            print*, 'wr dim',da_rec_no
                            call write_dimensions(th_da,i_ir,rec=da_rec_no)
                         else

                           call readwriteblocked_startwrite( & 
                                trim(tmpfile(quadrupel_filename(i_ir,quadrupel))), &
                                th, blocklength=blocklength )
                           th_pointer=> th
                           call write_dimensions(th,i_ir)
                         endif


                           if ( i_host .eq. 1 ) then
                              if ( integralpar_2cob_kin ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: writing kin")
                                if(options_directaccess_integrals()) then
!                                print*, 'wr kin',da_rec_no
                                 call write_integrals_2c_file(symadapt_int_2cob_kin(i_ir)%int,rec=da_rec_no )
                                else
                                 call write_integrals_2c_file(symadapt_int_2cob_kin(i_ir)%int)
                                endif
                              endif

                              if ( integralpar_2cob_nuc ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: writing nuc")
                                if(options_directaccess_integrals()) then
                                 call write_integrals_2c_file(symadapt_int_2cob_nuc(i_ir)%int,rec=da_rec_no)
                                else
                                 call write_integrals_2c_file(symadapt_int_2cob_nuc(i_ir)%int)
                                endif
                              endif

                              if ( integralpar_2cob_pvsp ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: writing pvsp")
                                if(options_directaccess_integrals()) then
                                 call write_integrals_2c_file(symadapt_int_2cob_pvsp(i_ir)%int,rec=da_rec_no)
                                else
                                 call write_integrals_2c_file(symadapt_int_2cob_pvsp(i_ir)%int)
                                endif
                              endif
!!!!!!!!                      if(integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                         if ( output_int_deeploops ) call write_to_output_units( &
!!!!!!!!                              "send_2cob3c_send: writing pseudo")
!!!!!!!!                        if(options_directaccess_integrals()) then
!!!!!!!!                         call write_integralfile_2c(symadapt_int_2cob_nuc_pseudo(i_ir)%int,rec=da_rec_no)
!!!!!!!!                        else
!!!!!!!!                         call write_integralfile_2c(symadapt_int_2cob_nuc_pseudo(i_ir)%int)
!!!!!!!!                        endif
!!!!!!!!                      endif
                           end if! i_host .eq. 1

                           if ( integralpar_2cob_ol ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: writing overlap")
                            if(options_directaccess_integrals()) then
                              call write_integrals_2c_file(symadapt_int_2cob_ol(i_ir)%int,rec=da_rec_no)
                            else
                              call write_integrals_2c_file(symadapt_int_2cob_ol(i_ir)%int)
                            endif        
                           endif


                           if ( integralpar_3c_xc ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: writing coulomb")
                              if(options_directaccess_integrals()) then
                              call write_integrals_3c_file( symadapt_int_3c_xc(i_ir)%int, &
                                   borders_xc(:,i_host), rec=da_rec_no)
                              else
                              call write_integrals_3c_file( symadapt_int_3c_xc(i_ir)%int, &
                                   borders_xc(:,i_host) )
                              endif
                           endif
                           if ( integralpar_3c_co ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: writing exchange")
                              if(options_directaccess_integrals()) then
                              call write_integrals_3c_file( symadapt_int_3c_co(i_ir)%int, &
                                   borders_ch(:,i_host),rec=da_rec_no)
                              else
                              call write_integrals_3c_file( symadapt_int_3c_co(i_ir)%int, &
                                   borders_ch(:,i_host) )
                              endif

                           endif

                           if ( integralpar_2cob_potential ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: writing potential 2")
                              if(options_directaccess_integrals()) then
                                 call write_integrals_3c_file( symadapt_int_2cob_poten(i_ir)%int, &
                                      borders_poten(:,i_host),rec=da_rec_no )
                              else
                                 call write_integrals_3c_file( symadapt_int_2cob_poten(i_ir)%int, &
                                      borders_poten(:,i_host) )
                               endif
                           endif
                           if ( integralpar_2cob_field ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: writing field")
                              if(options_directaccess_integrals()) then
                                 call write_integrals_3c_file( symadapt_int_2cob_field(i_ir)%int, &
                                      borders_field(:,i_host),rec=da_rec_no )
                              else
                                 call write_integrals_3c_file( symadapt_int_2cob_field(i_ir)%int, &
                                      borders_field(:,i_host) )
                              endif
                           endif

                         if(options_directaccess_integrals()) then
                           call readwriteblocked_stopwrite(th_da, total_length= &
                                quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ),rec=da_rec_no) 
                         else
                           call readwriteblocked_stopwrite(th, total_length= &
                                quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ) )
                         endif

                     enddo irreps

                     call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
                     !MOVE TO int_send_2cob3c_send: n_missing_quadrupels = n_missing_quadrupels - 1

                  else myhost

                     call start_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
                     call comm_init_send(i_host,msgtag_int_2cob3c_result)  
                     call quadrupel_pack(quadrupel)
                     call pack_dimensions()

                     irreps2: do i_ir = 1, n_ir

                           if ( i_host .eq. 1) then
                              if ( integralpar_2cob_kin ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: packing kin")
                                 DPRINT MyID, ' pack integralpar_2cob_kin '
                                 call pack_integrals_2c_file(symadapt_int_2cob_kin(i_ir)%int)
                              endif
                              if ( integralpar_2cob_nuc ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: packing nuc")
                                   DPRINT MyID, ' pack integralpar_2cob_nuc '
                                 call pack_integrals_2c_file(symadapt_int_2cob_nuc(i_ir)%int)
                              endif
                              if ( integralpar_2cob_pvsp ) then
                                 if ( output_int_deeploops ) call write_to_output_units( &
                                      "send_2cob3c_send: packing pvsp")
                                 DPRINT MyID, ' pack  integralpar_2cob_pvsp '
                                 call pack_integrals_2c_file(symadapt_int_2cob_pvsp(i_ir)%int)
                              endif
!!!!!!!!                      if(pseudopot_present.and.integralpar_relativistic.and.integralpar_2cob_nuc)then
!!!!!!!!                         if ( output_int_deeploops ) call write_to_output_units( &
!!!!!!!!                              "send_2cob3c_send: packing pseudo")
!!!!!!!!                       
!!!!!!!!                         call pack_integrals_2c_file(symadapt_int_2cob_nuc_pseudo(i_ir)%int)
!!!!!!!!                      endif
                           endif! i_host .eq. 1

                           if ( integralpar_2cob_ol ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: packing overlap")
                              call pack_integrals_2c_file(symadapt_int_2cob_ol(i_ir)%int)
                               DPRINT MyID, ' integralpar_2cob_ol '
                           endif
                           if ( integralpar_3c_xc ) then
                              DPRINT MyID, '  send integralpar_3c_xc'  !DG
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: packing exchange")
                              call pack_integrals_3c_file( symadapt_int_3c_xc(i_ir)%int, &
                                   borders_xc(:,i_host) )
                           endif
                           if ( integralpar_3c_co ) then
                               DPRINT MyID, '  send integralpar_3c_co for i_ir = ', i_ir
                               DPRINT MyID, 'sizeof(buffer) = ', size(symadapt_int_3c_co(i_ir)%int)
                               DPRINT MyID, ' pack buffer = '
                               DWRITE(*,'(10F14.7)') symadapt_int_3c_co(i_ir)%int
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: packing coulomb")
                              if(borders_ch(3,i_host).ne.0) &
                                call pack_integrals_3c_file( symadapt_int_3c_co(i_ir)%int, &
                                   borders_ch(:,i_host) )

                           endif
                           if ( integralpar_2cob_potential ) then
                              DPRINT MyID, '  send integralpar_2cob_potential'
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: packing potential")
                              call pack_integrals_3c_file( symadapt_int_2cob_poten(i_ir)%int, &
                                   borders_poten(:,i_host) )
                           endif
                           if ( integralpar_2cob_field ) then
                              if ( output_int_deeploops ) call write_to_output_units( &
                                   "send_2cob3c_send: packing field")
                              call pack_integrals_3c_file( symadapt_int_2cob_field(i_ir)%int, &
                                   borders_field(:,i_host) )
                           endif

                     end do irreps2

                     call stop_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
                     call start_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
                     call comm_send()   
                     call stop_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))

                  end if myhost

               end if parallel_filesystem! /else


            end do hosts

            call free_dimensions()

            if ( output_int_loops ) call write_to_output_units( &
                 "int_send_2cob3_send: done")

          contains !send


            subroutine calc_dimensions()
              !** End of interface *************************************
              implicit none
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_if1,i_if2,i_c1,n_if1_max, &
                   n_if1,n_if2,n_floats,i_if_c,status
              !------------ Executable code ----------------------------
              n_if1_max = maxval(ua1%symadapt_partner(:,quadrupel%l1)%N_independent_fcts)
              n_if_c_max = n_if1_max * n_c1
              allocate( n_floats_if_c(n_if_c_max,n_ir), n_if_c(n_ir), &
                   n_floats_if_c_max(n_ir), stat=status)
              if (status .ne. 0) call error_handler( &
                   " allocating at calc_dimensions() failed")
              do i_ir = 1,n_ir
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 n_if_c(i_ir) = n_if1 * n_c1
                 if ( n_if_c(i_ir) .gt. 0 ) then
                    i_if_c = 1
                    do i_if1 = 1, n_if1
                       if ( diagonal ) n_if2 = i_if1
                       do i_c1 = 1, n_c1
                          n_floats = 0
                          do i_if2 = 1, n_if2
                             if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                n_floats = n_floats + i_c1
                             else
                                n_floats = n_floats + n_c2
                             endif
                          enddo
                          n_floats_if_c(i_if_c,i_ir) = n_floats
                          i_if_c = i_if_c + 1
                       enddo
                    enddo
                    n_floats_if_c_max(i_ir) = maxval(n_floats_if_c(1:n_if_c(i_ir),i_ir))
                    DPRINT MyID, ' calc dimensions : i_ir =  ',i_ir ,' n_floats = ', n_floats_if_c_max(i_ir) 
                 else
                    n_floats_if_c_max(i_ir) = 0
                 endif
              enddo
              ! for the relativistic case we have to consider, that the two center 
              ! matrices are uncontracted
              if(integralpar_relativistic) then
                 n_if1_max = maxval(ua1%symadapt_partner(:,quadrupel%l1)%N_independent_fcts)
                 n_if_rel_max = n_if1_max * n_exp1
                 allocate( n_floats_if_rel(n_if_rel_max,n_ir), n_if_rel(n_ir), &
                      n_floats_if_rel_max(n_ir), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: allocating at calc_dimensions() failed")
                 do i_ir = 1,n_ir
                    n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                    n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                    n_if_rel(i_ir) = n_if1 * n_exp1
                    if ( n_if_rel(i_ir) .gt. 0 ) then
                       i_if_c = 1
                       do i_if1 = 1, n_if1
                          if ( diagonal ) n_if2 = i_if1
                          do i_c1 = 1, n_exp1
                             n_floats = 0
                             do i_if2 = 1, n_if2
                                if ( diagonal .and. i_if1 .eq. i_if2 ) then
                                   n_floats = n_floats + i_c1
                                else
                                   n_floats = n_floats + n_exp2
                                endif
                             enddo
                             n_floats_if_rel(i_if_c,i_ir) = n_floats
                             i_if_c = i_if_c + 1
                          enddo
                       enddo
                       n_floats_if_rel_max(i_ir) = maxval(n_floats_if_rel(1:n_if_rel(i_ir),i_ir))
                      
                    else
                       n_floats_if_rel_max(i_ir) = 0
                    endif
                 enddo
              end if
            end subroutine calc_dimensions


            subroutine write_dimensions(th,i_ir,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              type(readwriteblocked_tapehandle), intent(inout) :: th
              integer(kind=i4_kind), intent(in) :: i_ir
              integer(kind=i4_kind), intent(inout), optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: status
              real(kind=r8_kind) :: real_n_if_c(1)
              real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
              !------------ Executable code ----------------------------
              allocate( real_n_floats_if_c(n_if_c(i_ir)), stat=status)
              if (status .ne. 0) call error_handler( &
                   " allocating at write_dimensions() failed")
              real_n_if_c(1) = real(n_if_c(i_ir),r8_kind)
              real_n_floats_if_c = real(n_floats_if_c(1:n_if_c(i_ir),i_ir),r8_kind)
              if(present(rec)) then
              call readwriteblocked_write(real_n_if_c,th,rec=rec)
              call readwriteblocked_setrec(th,rec)
              call readwriteblocked_write(real_n_floats_if_c,th,rec=rec)
              call readwriteblocked_setrec(th,rec)
              else
              call readwriteblocked_write(real_n_if_c,th)
              call readwriteblocked_write(real_n_floats_if_c,th)
              endif
              deallocate( real_n_floats_if_c, stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at write_dimensions() failed")
              ! for relativistic case write also data for uncontracted basis sizes
              if(integralpar_relativistic) then
                 allocate( real_n_floats_if_c(n_if_rel(i_ir)), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: allocating at write_dimensions() failed")
                 real_n_if_c(1) = real(n_if_rel(i_ir),r8_kind)
                 real_n_floats_if_c = real(n_floats_if_rel(1:n_if_rel(i_ir),i_ir),r8_kind)
                if(present(rec)) then
                 call readwriteblocked_write(real_n_if_c,th,rec=rec)
                 call readwriteblocked_setrec(th,rec)
                 call readwriteblocked_write(real_n_floats_if_c,th,rec=rec)
                 call readwriteblocked_setrec(th,rec)
                else
                 call readwriteblocked_write(real_n_if_c,th)
                 call readwriteblocked_write(real_n_floats_if_c,th)
                endif
                 deallocate( real_n_floats_if_c, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: deallocating at write_dimensions() failed")
              end if
            end subroutine write_dimensions


            subroutine free_dimensions()
              !** End of interface *************************************
              implicit none
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: status
              !------------ Executable code ----------------------------
              deallocate( n_floats_if_c, n_if_c, n_floats_if_c_max, stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at free_dimensions() failed")
              if(integralpar_relativistic) then
                 deallocate( n_floats_if_rel, n_if_rel, n_floats_if_rel_max, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: deallocating at free_dimensions() failed")
              end if
            end subroutine free_dimensions


            subroutine pack_dimensions()
              !** End of interface *************************************
              implicit none
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: info, i_ir
              !------------ Executable code ----------------------------
              call commpack(n_if_c_max,info)
              if (info .ne. 0) call error_handler( &
                   "send_2cob3c_send: packing at pack_dimensions()  failed")
              call commpack(n_if_c,n_ir,1,info)
              if (info .ne. 0) call error_handler( &
                   "send_2cob3c_send: packing at pack_dimensions()  failed")
              do i_ir = 1, n_ir
                 call commpack(n_floats_if_c(:,i_ir),n_if_c_max,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing at pack_dimensions()  failed")
              enddo
              call commpack(n_floats_if_c_max,n_ir,1,info)
              if (info .ne. 0) call error_handler( &
                   "send_2cob3c_send: packing at pack_dimensions()  failed")

              if(integralpar_relativistic) then
                 call commpack(n_if_rel_max,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing at pack_dimensions()  failed")
                 call commpack(n_if_rel,n_ir,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing at pack_dimensions()  failed")
                 do i_ir = 1, n_ir
                    call commpack(n_floats_if_rel(:,i_ir),n_if_rel_max,1,info)
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_send: packing at pack_dimensions()  failed")
                 enddo
                 call commpack(n_floats_if_rel_max,n_ir,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing at pack_dimensions()  failed")
              end if
            end subroutine pack_dimensions


            subroutine write_integrals_3c_file(integral,borders,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
              integer(kind=i4_kind), intent(in) :: borders(3)
              integer(kind=i4_kind), intent(inout),optional :: rec 
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_if1,i_if2,i_c1,i_c2,n_if1,n_if2, &
                   nn_c2,n_floats,i_buf_lower,i_buf_upper,status
              real(kind=r8_kind), allocatable :: buffer(:)
              !------------ Executable code ----------------------------
              n_floats = borders(3) * n_floats_if_c_max(i_ir)
              allocate( buffer(n_floats),stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: allocating at writing 3 center failed")
              n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
              n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
              do i_if1 = 1, n_if1
                 if ( diagonal ) n_if2 = i_if1
                 do i_c1 = 1, n_c1
                    i_buf_lower = 1
                    i_buf_upper = 0
                    do i_if2 = 1, n_if2
                       if ( diagonal .and. i_if1 .eq. i_if2 ) then
                          nn_c2 = i_c1
                       else
                          nn_c2 = n_c2
                       endif
                       do i_c2 = 1, nn_c2
                          i_buf_upper = i_buf_upper + borders(3)
                          buffer(i_buf_lower:i_buf_upper) = &
                               integral(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
                          i_buf_lower = i_buf_lower + borders(3)
                       enddo
                    enddo
                    if(present(rec)) then
                     call readwriteblocked_write(buffer(1:i_buf_upper), thdaf_pointer, rec=rec )
                    else
                    call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
                    endif
                 enddo
              enddo
              deallocate(buffer,stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at writing 3 center failed")
            end subroutine write_integrals_3c_file


            subroutine write_integrals_2c_file(integral,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              real(kind=r8_kind), pointer :: integral(:,:,:,:)
              integer(kind=i4_kind), intent(inout), optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), allocatable :: buffer(:)
              integer(kind=i4_kind) :: i_if1,i_if2,i_c1,n_if1,n_if2, &
                   nn_c2,n_floats,i_buf_lower,i_buf_upper,status,n_dim1,n_dim2
              !------------ Executable code ----------------------------
              if(.not.integralpar_relativistic) then
                 n_floats = n_floats_if_c_max(i_ir)
                 n_dim1=n_c1
                 n_dim2=n_c2
              else
                 n_floats = n_floats_if_rel_max(i_ir)
                 n_dim1=n_exp1
                 n_dim2=n_exp2
              end if
              allocate( buffer(n_floats),stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: allocating at writing 2 center failed")
              n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
              n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
              do i_if1 = 1, n_if1
                 if ( diagonal ) n_if2 = i_if1
                 do i_c1 = 1, n_dim1
                    i_buf_lower = 1
                    i_buf_upper = 0
                    do i_if2 = 1, n_if2
                       if ( diagonal .and. i_if1 .eq. i_if2 ) then
                          nn_c2 = i_c1
                       else
                          nn_c2 = n_dim2
                       endif
                       i_buf_upper = i_buf_upper + nn_c2
                       buffer(i_buf_lower:i_buf_upper) = &
                            integral(1:nn_c2,i_c1,i_if2,i_if1)
                       i_buf_lower = i_buf_lower + nn_c2
                    enddo
                     if(present(rec)) then
!                     print*,i_buf_upper,rec,'size rec'
                        call readwriteblocked_write(buffer(1:i_buf_upper), thdaf_pointer, rec=rec)
                        call readwriteblocked_setrec(thdaf_pointer,rec)
                     else        
                       call readwriteblocked_write(buffer(1:i_buf_upper), th_pointer )
                     endif
                  enddo
              enddo
              deallocate(buffer,stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at writing 2 center failed")

            end subroutine write_integrals_2c_file


            subroutine pack_integrals_3c_file(integral,borders)
              implicit none
              !------------ Declaration of formal parameters -----------
              real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
              integer(kind=i4_kind), intent(in) :: borders(3)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_if1,i_if2,i_c1,i_c2,n_if1,n_if2, &
                   nn_c2,n_floats,i_buf_lower,i_buf_upper,status,info
              real(kind=r8_kind), allocatable :: buffer(:)
#ifndef FPP_NO_CHECKSUM
              real(kind=r8_kind) :: checksum
#endif
              !------------ Executable code ----------------------------
              DPRINT MyID,' running  pack_integralfile_3c with borders = ', borders
              n_floats = borders(3) * n_floats_if_c_max(i_ir)
              DPRINT MyID,' n_float = ', n_floats
              allocate( buffer(n_floats),stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: allocating at packing 3 center failed")
              n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
              n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
              DPRINT MyID, ' n_if1 = ', n_if1, ' n_if2 = ', n_if2
              DPRINT MyID, ' n_c1 = ', n_c1, ' n_c2 = ',n_c2
              DPRINT MyID, ' shape = ', shape (integral)
              do i_if1 = 1, n_if1
                 if ( diagonal ) n_if2 = i_if1
                 do i_c1 = 1, n_c1
                    i_buf_lower = 1
                    i_buf_upper = 0
                    do i_if2 = 1, n_if2
                       if ( diagonal .and. i_if1 .eq. i_if2 ) then
                          nn_c2 = i_c1
                       else
                          nn_c2 = n_c2
                       endif
                       do i_c2 = 1, nn_c2
                          i_buf_upper = i_buf_upper + borders(3)
                          buffer(i_buf_lower:i_buf_upper) = &
                               integral(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
                          i_buf_lower = i_buf_lower + borders(3)
                       enddo
                    enddo
                   if(n_floats.ne.0) then
                    call commpack(buffer,i_buf_upper,1,info)
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_send: packing 3 center failed")
#ifndef FPP_NO_CHECKSUM
                    checksum = sum(buffer(1:i_buf_upper))
                    call commpack(checksum,info)                                  !DG1
                    DPRINT MyID,'SEND_3c: checksum=',checksum,'n_floats= ',i_buf_upper
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_send: packing 3 center checksum failed")
#endif
                   endif
                 enddo
              enddo
              deallocate(buffer,stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at packing 3 center failed")
            end subroutine pack_integrals_3c_file


            subroutine pack_integrals_2c_file(integral)
              implicit none
              !------------ Declaration of formal parameters -----------
              real(kind=r8_kind), pointer :: integral(:,:,:,:)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), allocatable :: buffer(:)
              integer(kind=i4_kind) :: i_if1,i_if2,i_c1,info,n_if1,n_if2, &
                   nn_c2,n_floats,i_buf_lower,i_buf_upper,status,n_dim1,n_dim2
#ifndef FPP_NO_CHECKSUM
              real(kind=r8_kind) :: checksum
#endif
              !------------ Executable code ----------------------------
              DPRINT MyID, ' pack_integralfile_2c runns ...'
              DPRINT ' integralpar_relativistic = ', integralpar_relativistic
              if(.not.integralpar_relativistic) then
                 n_floats = n_floats_if_c_max(i_ir)
                 n_dim1=n_c1
                 n_dim2=n_c2
              else
                 n_floats = n_floats_if_rel_max(i_ir)
                 n_dim1=n_exp1
                 n_dim2=n_exp2
              end if
              DPRINT MyID, '  n_floats = ',  n_floats, ' n_dim1= ', n_dim1,  ' n_dim2= ', n_dim2
              allocate(buffer(n_floats),stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: allocating at packing 2 center failed")
              n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
              n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
              DPRINT MyID,' n_if1 = ', n_if1, ' n_if2 = ', n_if2
              do i_if1 = 1, n_if1
                 if ( diagonal ) n_if2 = i_if1
                 do i_c1 = 1, n_dim1
                    i_buf_lower = 1
                    i_buf_upper = 0
                    do i_if2 = 1, n_if2
                       if ( diagonal .and. i_if1 .eq. i_if2 ) then
                          nn_c2 = i_c1
                       else
                          nn_c2 = n_dim2
                       endif
                       i_buf_upper = i_buf_upper + nn_c2
                       buffer(i_buf_lower:i_buf_upper) = &
                            integral(1:nn_c2,i_c1,i_if2,i_if1)
                       i_buf_lower = i_buf_lower + nn_c2
                    enddo
                   if(n_floats.ne.0) then 
                    call commpack(buffer,i_buf_upper,1,info)
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_send: packing 2 center failed")
#ifndef FPP_NO_CHECKSUM
                    checksum = sum(buffer(1:i_buf_upper))
                    call commpack(checksum,info)
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_send: packing 3 center checksum failed")
                    DPRINT MyID,'SEND_2c: checksum=',checksum,'n_floats= ', i_buf_upper
#endif
                   endif
                 enddo
              enddo
              deallocate(buffer,stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_send: deallocating at packing 2 center failed")
            end subroutine pack_integrals_2c_file

          end subroutine int_send_2cob3c_send_file
          !*************************************************************



          !*************************************************************
          subroutine int_send_2cob3c_send_mem
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
            use int_data_2cob3c_module, only: symadapt_totsym_2c_int_type      &
                                            , symadapt_totsym_3c_int_type      &
                                            , symadapt_int_2cob_ol             &
                                            , symadapt_int_2cob_kin            &
                                            , symadapt_int_2cob_nuc            &
                                            , symadapt_int_2cob_field          &
                                            , symadapt_int_2cob_poten          &
                                            , symadapt_int_2cob_pvsp           &
                                            , symadapt_int_3c_co               &
                                            , symadapt_int_3c_xc               &
                                            , diagonal                         &
                                            , quadrupel                        &
                                            , n_c1                             &
                                            , n_c2                             &
                                            , n_exp1                           &
                                            , n_exp2                           &
                                            , ua1                              &
                                            , ua2 ! ua1 & ua2 really from here??
            implicit none
            !------------ Declaration of local variables -----------------
            integer(kind=i4_kind)                :: i_host, info, n_ir
!!$            integer(kind=i4_kind) :: i
        !!$    type(symadapt_totsym_2c_int_type), pointer :: symadap_ps(:)
            !------------ Executable code --------------------------------
            n_ir=symmetry_data_n_irreps()
            do i_host = first_host, last_host
               myhost: if ( i_host .eq. my_hostindex ) then
                  call start_timer(timer_int_write_2cob3c(integralpar_i_int_part))

                     if ( i_host .eq. 1 ) then
                        if ( integralpar_2cob_kin ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: writing kin")
                           if(.not.integralpar_relativistic) then
                              call store_integrals_2c_mem(symadapt_int_2cob_kin,integralstore_2cob_kin)
                           else
                              call store_integrals_2c_mem(symadapt_int_2cob_kin,integralstore_2cob_kin_rel)
                           end if
                        endif
                        if ( integralpar_2cob_nuc ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: writing nuc")
                           if(.not.integralpar_relativistic) then
                              call store_integrals_2c_mem(symadapt_int_2cob_nuc,integralstore_2cob_nuc)
                           else
                              call store_integrals_2c_mem(symadapt_int_2cob_nuc,integralstore_2cob_nuc_rel)
                           end if
                        endif
                        if ( integralpar_2cob_pvsp ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: writing nuc")
                           call store_integrals_2c_mem(symadapt_int_2cob_pvsp,integralstore_2cob_pvsp)
                        endif
!!!!!!!!                if (pseudopot_present.and.integralpar_relativistic) then
!!!!!!!!                   if ( output_int_loops ) call write_to_output_units( &
!!!!!!!!                        "send_2cob3c_send: writing pseudo")
!!!!!!!!                   call store_integrals_2c(symadapt_int_2cob_nuc_pseudo,integralstore_2cob_pseudo)
!!!!!!!!                end if
                     endif! i_host .eq. 1

                     if ( integralpar_2cob_ol ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: writing overlap")
                        if(.not.integralpar_relativistic) then
                           call store_integrals_2c_mem(symadapt_int_2cob_ol,integralstore_2cob_ol)
                        else
                           call store_integrals_2c_mem(symadapt_int_2cob_ol,integralstore_2cob_ol_rel)
                        end if
                     endif
                     if ( integralpar_3c_xc ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: writing coulomb")
                        call store_integrals_3c_mem( symadapt_int_3c_xc, &
                             integralstore_3c_xc, borders_xc(:,my_hostindex) )
                     endif
                     if ( integralpar_3c_co ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: writing exchange")
                        call store_integrals_3c_mem( symadapt_int_3c_co, &
                             integralstore_3c_co, borders_ch(:,my_hostindex) )
                     endif
                     if ( integralpar_2cob_potential ) then
                        if ( output_int_loops) call write_to_output_units(&
                                "send_2cob3c_send: writing potential 3")
                        call store_integrals_3c_mem(symadapt_int_2cob_poten,&
                                integralstore_3c_poten,&
                                borders_poten(:,my_hostindex) )
                     endif
                     if ( integralpar_2cob_field ) then
                        if ( output_int_loops) call write_to_output_units(&
                                "send_2cob3c_send: writing field")
                        call store_integrals_3c_mem(symadapt_int_2cob_field,&
                                integralstore_3c_field,&
                                borders_field(:,my_hostindex) )
                     endif

                  call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
                  !MOVE TO int_send_2cob3c_send: n_missing_quadrupels = n_missing_quadrupels - 1

               else myhost ! packing
                  call start_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
                  call comm_init_send(i_host,msgtag_int_2cob3c_result)
                  call quadrupel_pack(quadrupel)
                  call commpack(n_c1,info)
                  if (info .ne. 0) call error_handler( &
                       "send_2cob3c_send: packing n_c1 failed")
                  call commpack(n_c2,info)
                  if (info .ne. 0) call error_handler( &
                       "send_2cob3c_send: packing n_c2 failed")
                  if(integralpar_relativistic) then
                     call commpack(n_exp1,info)
                     if (info .ne. 0) call error_handler( &
                          "send_2cob3c_send: packing n_exp1 failed")
                     call commpack(n_exp2,info)
                     if (info .ne. 0) call error_handler( &
                          "send_2cob3c_send: packing n_exp2 failed")
                  end if

                     if ( i_host .eq. 1 ) then
                        if ( integralpar_2cob_kin ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: packing kin")
                           call pack_integrals_2c_mem(symadapt_int_2cob_kin)
                        endif
                        if ( integralpar_2cob_nuc ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: packing nuc")
                           call pack_integrals_2c_mem(symadapt_int_2cob_nuc)
                        endif
                        if ( integralpar_2cob_pvsp ) then
                           if ( output_int_loops ) call write_to_output_units( &
                                "send_2cob3c_send: packing pvsp")
                           call pack_integrals_2c_mem(symadapt_int_2cob_pvsp)
                        endif
!!!!!!!!                if ( pseudopot_present.and.integralpar_relativistic.and.integralpar_2cob_nuc)then
!!!!!!!!                   if ( output_int_loops ) call write_to_output_units( &
!!!!!!!!                        "send_2cob3c_send: packing pseudo")
!!!!!!!!                   call pack_integrals_2c_mem(symadapt_int_2cob_nuc_pseudo)
!!!!!!!!                end if
                     endif! i_host .eq. 1

                     if ( integralpar_2cob_ol ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: packing overlap")
                        call pack_integrals_2c_mem(symadapt_int_2cob_ol)
                     endif
                     if ( integralpar_3c_xc ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: packing coulomb")
                        call pack_integrals_3c_mem( symadapt_int_3c_xc, &
                             borders_xc(:,i_host) )
                     endif
                     if ( integralpar_3c_co ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: packing exchange")
                        call pack_integrals_3c_mem( symadapt_int_3c_co, &
                             borders_ch(:,i_host) )
                     endif
                     if ( integralpar_2cob_potential ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: packing poten")
                        call pack_integrals_3c_mem( symadapt_int_2cob_poten, &
                             borders_poten(:,i_host) )
                     endif
                     if ( integralpar_2cob_field ) then
                        if ( output_int_loops ) call write_to_output_units( &
                             "send_2cob3c_send: packing field")
                        call pack_integrals_3c_mem( symadapt_int_2cob_field, &
                             borders_field(:,i_host) )
                     endif

                  call stop_timer(timer_int_pack_2cob3c(integralpar_i_int_part))
                  call start_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
                  call comm_send()
                  call stop_timer(timer_int_commsend_2cob3c(integralpar_i_int_part))
               end if myhost ! i_host.eq.my_hostindex
            enddo! hosts

            if ( output_int_loops ) call write_to_output_units( &
                 "int_send_2cob3_send: done")

          contains


            subroutine store_integrals_3c_mem(symadapt_int, integralstore, border)
              !------------ Declaration of formal parameters -----------
              type(symadapt_totsym_3c_int_type), intent(in) :: symadapt_int(:)
              real(kind=r8_kind), intent(out) :: integralstore(:)
              integer(kind=i4_kind), intent(in) :: border(3)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_c2,i_meta, &
                   n_if1,n_if2,nn_c2
              !------------ Executable code ----------------------------
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 if ( .not. diagonal ) n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 integral => symadapt_int(i_ir)%int
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_c1
                       i_meta = (metaindex(quadrupel,i_ir,i_if1,i_c1) - 1) * border(3) + 1
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_c2
                          endif
                          do i_c2 = 1, nn_c2
                             integralstore(i_meta:i_meta+border(3)-1) = &
                                  integral(i_c2,i_c1,border(1):border(2),i_if2,i_if1)
                             i_meta = i_meta + border(3)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
            end subroutine store_integrals_3c_mem


            subroutine store_integrals_2c_mem(symadapt_int, integralstore)
              !------------ Declaration of formal parameters -----------
              type(symadapt_totsym_2c_int_type), intent(in) :: symadapt_int(:)
              real(kind=r8_kind), intent(out) :: integralstore(:)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), pointer :: integral(:,:,:,:)
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_meta, &
                   n_if1,n_if2,nn_c2,n_dim1,n_dim2
              logical :: rel
              !------------ Executable code ----------------------------
              if(.not.integralpar_relativistic) then
                 rel=.false.
                 n_dim1=n_c1
                 n_dim2=n_c2
              else
                 rel=.true.
                 n_dim1=n_exp1
                 n_dim2=n_exp2
              end if
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 if ( .not. diagonal ) n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 integral => symadapt_int(i_ir)%int
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_dim1
                       i_meta = metaindex(quadrupel,i_ir,i_if1,i_c1,rel)
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_dim2
                          endif
                          integralstore(i_meta:i_meta+nn_c2-1) = &
                               integral(1:nn_c2,i_c1,i_if2,i_if1)
                          i_meta = i_meta + nn_c2
                       enddo
                    enddo
                 enddo
              enddo
            end subroutine store_integrals_2c_mem

            subroutine pack_integrals_2c_mem(symadapt_int)
              !------------ Declaration of formal parameters -----------
              type(symadapt_totsym_2c_int_type), intent(in) :: symadapt_int(:)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), pointer :: integral(:,:,:,:)
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
                   n_if1,n_if2,nn_c2,i_buf_lower,i_buf_upper,status,&
                   n_dim1,n_dim2
              real(kind=r8_kind), allocatable :: buffer(:)
              !------------ Executable code ----------------------------
              if(.not.integralpar_relativistic) then
                 n_dim1=n_c1
                 n_dim2=n_c2
              else
                 n_dim1=n_exp1
                 n_dim2=n_exp2
              end if
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 if ( .not. diagonal ) n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 integral => symadapt_int(i_ir)%int
                 if(size(integral).eq.0) cycle
                 allocate( buffer( size(integral) ),stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: allocating at packing 2 center failed")
                 i_buf_lower = 1
                 i_buf_upper = 0
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_dim1
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_dim2
                          endif
                          i_buf_upper = i_buf_upper + nn_c2
                          buffer(i_buf_lower:i_buf_upper) = &
                               integral(1:nn_c2,i_c1,i_if2,i_if1)
                          i_buf_lower = i_buf_lower + nn_c2
                       enddo
                    enddo
                 enddo
                 call commpack(buffer,i_buf_upper,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing 2 center failed")
                 deallocate(buffer,stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: deallocating at packing 2 center failed")
              enddo
            end subroutine pack_integrals_2c_mem

            subroutine pack_integrals_3c_mem(symadapt_int, borders)
              !------------ Declaration of formal parameters -----------
              type(symadapt_totsym_3c_int_type), intent(in) :: symadapt_int(:)
              integer(kind=i4_kind), intent(in) :: borders(3)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              real(kind=r8_kind), pointer :: integral(:,:,:,:,:)
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,i_c2,info, &
                   n_if1,n_if2,nn_c2,i_buf_lower,i_buf_upper,status,alloc_size
              real(kind=r8_kind), allocatable :: buffer(:)
              !------------ Executable code ----------------------------
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 integral => symadapt_int(i_ir)%int
                 alloc_size=  size(integral,1) * size(integral,2) * size(integral,4) *  &
                              size(integral,5) * borders(3)
                 if(alloc_size.eq.0) cycle
                 allocate( buffer( size(integral,1) * size(integral,2) * size(integral,4) *  &
                      size(integral,5) * borders(3) ), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: allocating at packing 3 center failed")
                 i_buf_lower = 1
                 i_buf_upper = 0
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_c1
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_c2
                          endif
                          do i_c2 = 1, nn_c2
                             i_buf_upper = i_buf_upper + borders(3)
                             buffer(i_buf_lower:i_buf_upper) = &
                                  integral(i_c2,i_c1,borders(1):borders(2),i_if2,i_if1)
                             i_buf_lower = i_buf_lower + borders(3)
                          enddo
                       enddo
                    enddo
                 enddo

                 call commpack(buffer,i_buf_upper,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_send: packing 3 center failed")
                 deallocate(buffer,stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_send: deallocating at packing 3 center failed")
              enddo
            end subroutine pack_integrals_3c_mem

          end subroutine int_send_2cob3c_send_mem
          !*************************************************************

          !*************************************************************
          subroutine int_send_2cob3c_receive_all
            !  Purpose:
            !    + receiving of contracted and symmetry adapted
            !      integrals for THE REST OF quadrupels
            use comm, only: comm_barrier
            implicit none
            ! *** end of interface ***

            DPRINT  MyID,'int_send_2cob3c_receive_all: entered'
            DPRINT  MyID,'int_send_2cob3c_receive_all: n_missing_quadrupels=',n_missing_quadrupels

            if( options_integrals_on_file() .and. filesystem_is_parallel )then
              ! the other hosts didnt need to send any data (quadrupels),
              ! but have written it directly to my disk instead --
              ! so the number of missing quadrupels doesnt matter.
              n_missing_quadrupels = 0

              ! Just wait until the other slaves finsh their work
              ! (writing data to my disks and come to this barrier):
              call comm_barrier()
              DPRINT  MyID,'int_send_2cob3c_receive_all: filesystem_is_parallel, return!'
            else
              ! wait untill all data (quadrupels) arrive and
              ! either write them to local discs or into memory:
              do while (n_missing_quadrupels .gt. 0 )

                 if ( output_int_loops ) call write_to_output_units( &
                      "int_send_2cob3c_receive_all: waiting for quadrupels: ",n_missing_quadrupels)

                 call comm_save_recv(comm_all_other_hosts,msgtag_int_2cob3c_result)

                 if ( output_slaveoperations ) &
                      call write_to_output_units("int_send_2cob3c_receive_all: integral_2cob_result")

                 call int_send_2cob3c_receive() ! decrements n_missing_quadrupels by 1
              enddo
            endif
            DPRINT  MyID,'int_send_2cob3c_receive_all: exit'
          end subroutine int_send_2cob3c_receive_all
          !*************************************************************

          !*************************************************************
          subroutine int_send_2cob3c_receive
            !  Purpose:
            !    + receiving of contracted and symmetry adapted
            !      integrals for one quadrupel
            !    + Storing them either in memory or on file
            !  called by: main_slave, integral_main_2cob3c,
            !    integral_interrupt_2cob3c, int_send_2cob3c_shutdown
            ! used if .not. filesystem_is_parallel in integral part
            !** End of interface ***************************************
             DPRINT MyID, ' int_send_2cob3c_receive run...'
            if ( options_integrals_on_file() ) then

               call int_send_2cob3c_receive_file()
            else
               call int_send_2cob3c_receive_mem()
            endif
            ! MOVED HERE from inside of the above two:
            n_missing_quadrupels = n_missing_quadrupels - 1
          end subroutine int_send_2cob3c_receive
          !*************************************************************

          !*************************************************************
          subroutine int_send_2cob3c_receive_file
            !  Purpose:
            !    + receiving of contracted and symmetry adapted
            !      integrals for one quadrupel
            !    + Storing them in direct access files in scf symmetric
            !      storage mode with the record number as metindex
            !  called by: main_slave, integral_main_2cob3c,
            !    integral_interrupt_2cob3c, int_send_2cob3c_shutdown
            !** End of interface *****************************************
          use options_module, only: options_directaccess_integrals,quadrupels_reclength
            implicit none
            !------------ Declaration of local variables -----------------
            type(readwriteblocked_tapehandle) :: th
            type(quadrupel_type) :: quadrupel
            integer(kind=i4_kind) :: i_ir, n_ir, n_if_c_max, n_if_rel_max
            logical :: diagonal, restart_timer
            integer(kind=i4_kind), allocatable :: &
                 n_floats_if_c(:,:), n_floats_if_c_max(:), n_if_c(:), &
                 n_floats_if_rel(:,:), n_floats_if_rel_max(:), n_if_rel(:)
            ! n_floats_if_c(i_if_c,i_ir), n_floats_if_c_max(i_ir), n_if_c(i_ir)
            ! n_floats_if_rel(i_if_c,i_ir), n_floats_if_rel_max(i_ir), n_if_rel(i_ir)
            !------------ Executable code --------------------------------
             DPRINT  MyID, ' run int_send_2cob3c_receive_file'
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
                    "int_send_2cob3c_receive: start with quadrupel ", &
                    quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
               write(stdout_unit,*) &
                    "int_send_2cob3c_receive: start with quadrupel ", &
                    quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
            endif

            diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
                       (quadrupel%l1 == quadrupel%l2)

            n_ir = symmetry_data_n_irreps()

            call unpack_dimensions()

            do i_ir = 1, n_ir

                if(options_directaccess_integrals()) then
                  call readwriteblocked_startwrite( &
                       trim(tmpfile("da")), &
                       th_da, blocklength=quadrupels_reclength,rec=da_rec_no )
                            darec_no(i_ir,index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                         index_ua_l(quadrupel%ua2,quadrupel%l2))=da_rec_no
                  call write_dimensions(th_da,i_ir,rec=da_rec_no)
                else
                  call readwriteblocked_startwrite( & ! quadrupel_filename send_2cob3c_receive_file
                       trim(tmpfile(quadrupel_filename(i_ir,quadrupel))), &
                       th, blocklength=blocklength )
                  call write_dimensions(th,i_ir)
                endif
                 master: if ( comm_i_am_master()) then
                     DPRINT MyID, ' I am master'
                     if ( integralpar_2cob_kin ) then
                        DPRINT  MyID, ' Receive integralpar_2cob_kin'
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "send_2cob3c_receive: unpacking kin")
                        if(.not.integralpar_relativistic) then
                           call unpack_integralfile(1) !18
                        else
                           call unpack_integralfile(1,rel_dummy=.true.) !17
                        end if
                     endif
                     if ( integralpar_2cob_nuc ) then
                         DPRINT  MyID, ' Receive integralpar_2cob_nuc'
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "int_send_2cob3c_receive: unpacking nuc")
                        if(.not.integralpar_relativistic) then
                           call unpack_integralfile(1) !16
                        else
                           call unpack_integralfile(1,rel_dummy=.true.) !15
                        end if
                     endif
                     if ( integralpar_2cob_pvsp ) then
                         DPRINT  MyID, ' Receive integralpar_2cob_pvsp'
                        if ( output_int_deeploops ) call write_to_output_units( &
                             "int_send_2cob3c_receive: unpacking pvsp")
                        if(.not.integralpar_relativistic) then
                           call unpack_integralfile(1) !14
                        else
                           call unpack_integralfile(1,rel_dummy=.true.)!13
                        end if
                     endif
                     !mdf>>> ! moved inside if
!!!!!!!!             if (integralpar_2cob_nuc.and.pseudopot_present.and.integralpar_relativistic)then
!!!!!!!!                 DPRINT  MyID, ' Receive pseudo'
!!!!!!!!                if ( output_int_deeploops ) call write_to_output_units( &
!!!!!!!!                     "int_send_2cob3c_receive: unpacking pseudo")
!!!!!!!!                call unpack_integralfile(1,rel_dummy=.true.) !12
!!!!!!!!             endif
                  endif master! comm_i_am_master()

                  if ( integralpar_2cob_ol ) then
                     if ( output_int_deeploops ) call write_to_output_units( &
                          "int_send_2cob3c_receive: unpacking overlap")
                      DPRINT  MyID, ' Receive overlap'
!                      print*,"unpacking overlap",comm_myindex()
                     if(.not.integralpar_relativistic) then
                        call unpack_integralfile(1) !11
                     else
                        call unpack_integralfile(1,rel_dummy=.true.) !10
                     end if
                  endif

                  if ( integralpar_3c_xc ) then
                      DPRINT  MyID, 'integralpar_3c_xc '
                     if ( output_int_deeploops ) call write_to_output_units( &
                          "send_2cob3c_receive: unpacking exchange")
                     !print*,"unpacking exchange",comm_myindex()
                     call unpack_integralfile(borders_xc(3,my_hostindex)) !9
                  endif
                  if ( integralpar_3c_co ) then
                     if ( output_int_deeploops ) call write_to_output_units( &
                          "send_2cob3c_receive: unpacking coulomb")
                     call unpack_integralfile(borders_ch(3,my_hostindex))    !8
                       DPRINT  MyID, 'integralpar_3c_co '
                  endif
                  if ( integralpar_2cob_potential ) then
                     DPRINT  MyID, ' integralpar_2cob_potential '
                     if ( output_int_deeploops ) call write_to_output_units( &
                          "send_2cob3c_receive: unpacking poten")
                     call unpack_integralfile(borders_poten(3,my_hostindex))   !7
                  endif
                  if ( integralpar_2cob_field ) then
                        DPRINT  MyID, ' integralpar_2cob_field '
                     if ( output_int_deeploops ) call write_to_output_units( &
                          "send_2cob3c_receive: unpacking field")
                     call unpack_integralfile(borders_field(3,my_hostindex))    !6
                  endif

                  if(options_directaccess_integrals()) then
                  call readwriteblocked_stopwrite(th_da, total_length= &
                       quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                                      index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ),rec=da_rec_no)
                  else
                  call readwriteblocked_stopwrite(th, total_length= &
                       quadrupelfile_length( index_ua_l(quadrupel%ua1,quadrupel%l1), &
                          index_ua_l(quadrupel%ua2,quadrupel%l2), i_ir ) )
                  endif

            enddo! loop over irreps

            call free_dimensions()

            !MOVED TO _receive: n_missing_quadrupels = n_missing_quadrupels - 1

            call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
            call stop_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
            if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

            if ( output_int_loops ) call write_to_output_units( &
                 "int_send_2cob3c_receive: done")


          contains


            subroutine write_dimensions(th,i_ir,rec)
              implicit none
              !------------ Declaration of formal parameters -----------
              type(readwriteblocked_tapehandle), intent(inout) :: th
              integer(kind=i4_kind), intent(in) :: i_ir
              integer(kind=i4_kind), intent(inout),optional :: rec
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: status
              real(kind=r8_kind) :: real_n_if_c(1)
              real(kind=r8_kind), allocatable :: real_n_floats_if_c(:)
              !------------ Executable code ----------------------------
              allocate( real_n_floats_if_c(n_if_c(i_ir)), stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_receive: allocating at write_dimensions() failed")
              real_n_if_c(1) = real(n_if_c(i_ir),r8_kind)
              real_n_floats_if_c = real(n_floats_if_c(1:n_if_c(i_ir),i_ir),r8_kind)
              if(present(rec)) then
              call readwriteblocked_write(real_n_if_c,th,rec=rec)
              call readwriteblocked_write(real_n_floats_if_c,th,rec=rec)
              else
              call readwriteblocked_write(real_n_if_c,th)
              call readwriteblocked_write(real_n_floats_if_c,th)
              endif
              deallocate( real_n_floats_if_c, stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_receive: deallocating at write_dimensions() failed")
              if(integralpar_relativistic) then
                 allocate( real_n_floats_if_c(n_if_rel(i_ir)), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: allocating at write_dimensions() failed")
                 real_n_if_c(1) = real(n_if_rel(i_ir),r8_kind)
                 real_n_floats_if_c = real(n_floats_if_rel(1:n_if_rel(i_ir),i_ir),r8_kind)
                 if(present(rec)) then
                 call readwriteblocked_write(real_n_if_c,th,rec=rec)
                 call readwriteblocked_write(real_n_floats_if_c,th,rec=rec)
                 else
                 call readwriteblocked_write(real_n_if_c,th)
                 call readwriteblocked_write(real_n_floats_if_c,th)
                 endif
                 deallocate( real_n_floats_if_c, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: deallocating at write_dimensions() failed")
              end if
            end subroutine write_dimensions


            subroutine free_dimensions()
              !** End of interface *************************************
              implicit none
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: status
              !------------ Executable code ----------------------------
              deallocate( n_floats_if_c, n_if_c, n_floats_if_c_max, stat=status)
              if (status .ne. 0) call error_handler( &
                   "int_send_2cob3c_receive: deallocating at free_dimensions() failed")
              if(integralpar_relativistic) then
                 deallocate( n_floats_if_rel, n_if_rel, n_floats_if_rel_max, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: deallocating at free_dimensions() failed")
              end if
            end subroutine free_dimensions


            subroutine unpack_dimensions()
              !** End of interface *************************************
              implicit none
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: info, status, i_ir
              !------------ Executable code ----------------------------
              call communpack(n_if_c_max,info)
              if (info .ne. 0) call error_handler( &
                   "send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
              allocate( n_floats_if_c(n_if_c_max,n_ir), &
                   n_if_c(n_ir), n_floats_if_c_max(n_ir), stat=status)
              if (status .ne. 0) call error_handler( &
                   "send_2cob3c_receive: allocating at unpack_dimensions() failed")
              call communpack(n_if_c,n_ir,1,info)
              if (info .ne. 0) call error_handler( &
                   "send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
              do i_ir = 1, n_ir
                 call communpack(n_floats_if_c(:,i_ir),n_if_c_max,1,info)
                 if (info .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
              enddo
              call communpack(n_floats_if_c_max,n_ir,1,info)
              if (info .ne. 0) call error_handler( &
                   "int_send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
              if(integralpar_relativistic) then
                 call communpack(n_if_rel_max,info)
                 if (info .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
                 allocate( n_floats_if_rel(n_if_rel_max,n_ir), &
                   n_if_rel(n_ir), n_floats_if_rel_max(n_ir), stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: allocating at unpack_dimensions() failed")
                 call communpack(n_if_rel,n_ir,1,info)
                 if (info .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
                 do i_ir = 1, n_ir
                    call communpack(n_floats_if_rel(:,i_ir),n_if_rel_max,1,info)
                    if (info .ne. 0) call error_handler( &
                         "send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
                 enddo
                 call communpack(n_floats_if_rel_max,n_ir,1,info)
                 if (info .ne. 0) call error_handler( &
                      "send_2cob3c_receive: unpacking at unpack_dimensions()  failed")
              end if
            end subroutine unpack_dimensions


            subroutine unpack_integralfile(bordersize,rel_dummy)
              implicit none
              !------------ Declaration of formal parameters -----------
              integer(kind=i4_kind), intent(in) :: bordersize
              logical,optional :: rel_dummy
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: info,n_floats,status,i_if_c
              real(kind=r8_kind), allocatable :: buffer(:)
#ifndef FPP_NO_CHECKSUM
              real(kind=r8_kind) :: checksum_received, checksum_calculated
#endif
              character(len=200) :: string
              logical :: rel
              !------------ Executable code ----------------------------
              DPRINT MyID, ' run unpack_integralfile... for i_ir =  ', i_ir, ' borgersize = ', bordersize
              rel=.false.
              if(present(rel_dummy)) then
                 rel=rel_dummy
                 DPRINT MyID, ' present rel_dummy = ',rel_dummy
              end if
              DPRINT 'rel = ', rel
              not_rel: if(.not.rel) then
                 n_floats = n_floats_if_c_max(i_ir) * bordersize


                 allocate(buffer(n_floats),stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_receive: allocating at unpacking integrals failed")
                 do i_if_c = 1, n_if_c(i_ir)
                    n_floats = n_floats_if_c(i_if_c,i_ir) * bordersize
                    if(n_floats.ne.0)  then
                     call communpack(buffer,n_floats,1,info)
                     if (info .ne. 0) call error_handler( &
                         "send_2cob3c_receive: unpacking integrals failed")
#ifndef FPP_NO_CHECKSUM
                     checksum_calculated = sum(buffer(1:n_floats))
                     DPRINT MyID,'RECV_nr: checksum=',checksum_calculated,'n_floats=',n_floats
                     call communpack(checksum_received,info)
                     if (info .ne. 0) call error_handler( &
                         "send_2cob3c_receive: unpacking checksum failed")
                     if ( checksum_calculated .ne. checksum_received) then
                       write(string,*) "send_2cob3c_receive: checksum error 1: received ", &
                            checksum_received, " and calculated ", checksum_calculated   &
                            ,"diff=",checksum_received-checksum_calculated
                       ! call error_handler(trim(string))
                       WARN(trim(string))
                     endif
#endif
                    endif
                    if(options_directaccess_integrals()) then
!                    print*,comm_myindex(),'size->',n_floats
                    call readwriteblocked_write(buffer(1:n_floats),th_da,rec=da_rec_no)
                    else
                    call readwriteblocked_write(buffer(1:n_floats),th)
                    endif
                 enddo
                 deallocate(buffer, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: deallocating at unpacking failed")
              else not_rel
                 n_floats = n_floats_if_rel_max(i_ir) * bordersize
                 allocate(buffer(n_floats),stat=status)
                 if (status .ne. 0) call error_handler( &
                      "send_2cob3c_receive: allocating at unpacking integrals failed")
                 do i_if_c = 1, n_if_rel(i_ir)
                    n_floats = n_floats_if_rel(i_if_c,i_ir) * bordersize
                    if(n_floats.ne.0)  then
                     call communpack(buffer,n_floats,1,info)
                     if (info .ne. 0) call error_handler( &
                         "int_send_2cob3c_receive: unpacking integrals failed")
#ifndef FPP_NO_CHECKSUM
                     checksum_calculated = sum(buffer(1:n_floats))
                     DPRINT MyID,'RECV_rl: checksum=',checksum_calculated,'n_floats=',n_floats
                     call communpack(checksum_received,info)
                     if (info .ne. 0) call error_handler( &
                         "send_2cob3c_receive: unpacking checksum failed")
                     if ( checksum_calculated .ne. checksum_received) then
                        write(string,*) "int_send_2cob3c_receive: checksum error 2: received ", &
                             checksum_received, " and calculated ", checksum_calculated        &
                            ,"diff=",checksum_received-checksum_calculated
                       ! call error_handler(trim(string))
                       WARN(trim(string))
                     endif
#endif
                    endif
                       if(options_directaccess_integrals()) then
!                       print*,comm_myindex(),da_rec_no,n_floats
                       call readwriteblocked_write(buffer(1:n_floats),th_da,rec=da_rec_no)
                       else
                       call readwriteblocked_write(buffer(1:n_floats),th)
                       endif
                 enddo

                 deallocate(buffer, stat=status)
                 if (status .ne. 0) call error_handler( &
                      "int_send_2cob3c_receive: deallocating at unpacking failed")
              endif not_rel
            end subroutine unpack_integralfile

          end subroutine int_send_2cob3c_receive_file
          !*************************************************************




          !*************************************************************
          subroutine int_send_2cob3c_receive_mem
            !  Purpose:
            !    + receiving of contracted and symmetry adapted
            !      integrals for one quadrupel
            !    + Storing them in memory in integralstore_module
            !  called by: main_slave, integral_main_2cob3c,
            !    integral_interrupt_2cob3c, int_send_2cob3c_shutdown
            !** End of interface *****************************************
            implicit none
            !------------ Declaration of local variables -----------------
            type(quadrupel_type)            :: quadrupel
            integer(kind=i4_kind)           :: n_c1, n_c2, n_exp1, n_exp2, &
                                               info
!!$            integer(kind=i4_kind)           :: i
            type(unique_atom_type), pointer :: ua1,ua2
            logical                         :: diagonal, restart_timer
            !------------ Executable code --------------------------------

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
                    "int_send_2cob3c_receive: start with quadrupel ", &
                    quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
               write(stdout_unit,*) &
                    "int_send_2cob3c_receive: start with quadrupel ", &
                    quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
            endif

            ua1 => unique_atoms(quadrupel%ua1)
            ua2 => unique_atoms(quadrupel%ua2)
            diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
                       (quadrupel%l1 == quadrupel%l2)
            call communpack(n_c1,info)
            if (info .ne. 0) call error_handler( &
                 "int_send_2cob3c_receive: unpacking n_c1 failed")
            call communpack(n_c2,info)
            if (info .ne. 0) call error_handler( &
                 "int_send_2cob3c_receive: unpacking n_c2 failed")

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
                        call unpack_integrals_2c_mem(integralstore_2cob_kin)
                     else
                        call unpack_integrals_2c_mem(integralstore_2cob_kin_rel)
                     end if
                  endif
                  !** this is subroutine int_send_2cob3c_receive_mem
                  if ( integralpar_2cob_nuc ) then
                     if ( output_int_loops ) call write_to_output_units( &
                          "int_send_2cob3c_receive: unpacking nuc")
                     if(.not.integralpar_relativistic) then
                        call unpack_integrals_2c_mem(integralstore_2cob_nuc)
                     else
                        call unpack_integrals_2c_mem(integralstore_2cob_nuc_rel)
                     end if
                  endif
                  if ( integralpar_2cob_pvsp ) then
                     if ( output_int_loops ) call write_to_output_units( &
                          "int_send_2cob3c_receive: unpacking pvsp")
                     call unpack_integrals_2c_mem(integralstore_2cob_pvsp)
                  endif
!!!!!!!!          if(pseudopot_present.and.integralpar_relativistic.and.integralpar_2cob_nuc)then
!!!!!!!!             if ( output_int_loops ) call write_to_output_units( &
!!!!!!!!                  "int_send_2cob3c_receive: unpacking pseudo")
!!!!!!!!             call unpack_integrals_2c_mem(integralstore_2cob_pseudo)
!!!!!!!!          end if
               endif! comm_i_am_master()

               if ( integralpar_2cob_ol ) then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_receive: unpacking overlap")
                  if(.not.integralpar_relativistic) then
                     call unpack_integrals_2c_mem(integralstore_2cob_ol)
                  else
                     call unpack_integrals_2c_mem(integralstore_2cob_ol_rel)
                  end if
               endif
               if ( integralpar_3c_xc ) then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_receive: unpacking exchange")
                  call unpack_integrals_3c_mem(integralstore_3c_xc, &
                       borders_xc(3,my_hostindex))
               endif
               if ( integralpar_3c_co ) then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_receive: unpacking coulomb")
                  if(borders_ch(3,my_hostindex).ne.0) &
                     call unpack_integrals_3c_mem(integralstore_3c_co, &
                       borders_ch(3,my_hostindex))
               endif
               if ( integralpar_2cob_potential ) then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_receive: unpacking poten")
                  call unpack_integrals_3c_mem(integralstore_3c_poten, &
                       borders_poten(3,my_hostindex))
               endif
               if ( integralpar_2cob_field ) then
                  if ( output_int_loops ) call write_to_output_units( &
                       "int_send_2cob3c_receive: unpacking field")
                  call unpack_integrals_3c_mem(integralstore_3c_field, &
                       borders_field(3,my_hostindex))
               endif

            !MOVED TO _receive: n_missing_quadrupels = n_missing_quadrupels - 1
            call stop_timer(timer_int_write_2cob3c(integralpar_i_int_part))
            call stop_timer(timer_int_receive_2cob3c(integralpar_i_int_part))
            if (restart_timer) call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

            if ( output_int_loops ) call write_to_output_units( &
                 "int_send_2cob3c_receive: done")

          contains


            subroutine unpack_integrals_2c_mem(integralstore)
              !------------ Declaration of formal parameters -----------
              real(kind=r8_kind), intent(out) :: integralstore(:)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
                   n_if1,n_if2,nn_c2,i_meta,n_dim1,n_dim2
              logical :: rel
              !------------ Executable code ----------------------------
              if(.not.integralpar_relativistic) then
                 rel=.false.
                 n_dim1=n_c1
                 n_dim2=n_c2
              else
                 rel=.true.
                 n_dim1=n_exp1
                 n_dim2=n_exp2
              end if
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 if(n_if1.le.0) cycle
                 if ( .not. diagonal ) then
                 n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                 if(n_if2.eq.0) cycle
                 endif
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_dim1
                       i_meta = metaindex(quadrupel,i_ir,i_if1,i_c1,rel)
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_dim2
                          endif
                          call communpack(integralstore(i_meta:i_meta+nn_c2-1),nn_c2,1,info)
                          if (info .ne. 0) call error_handler( &
                               "int_send_2cob3c_receive: unpacking 3 center failed")
                          i_meta = i_meta + nn_c2
                       enddo
                    enddo
                 enddo
              enddo
            end subroutine unpack_integrals_2c_mem

            subroutine unpack_integrals_3c_mem(integralstore,n_receive_ff)
              !------------ Declaration of formal parameters -----------
              integer(kind=i4_kind), intent(in) :: n_receive_ff
              real(kind=r8_kind), intent(out) :: integralstore(:)
              !** End of interface *************************************
              !------------ Declaration of local variables -------------
              integer(kind=i4_kind) :: i_ir,i_if1,i_if2,i_c1,info, &
                   n_if1,n_if2,nn_c2,i_meta,i_c2
              !------------ Executable code ----------------------------
              do i_ir = 1,symmetry_data_n_irreps()
                 n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
                 if(n_if1.eq.0) cycle
                 if ( .not. diagonal ) then
                  n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
                  if(n_if2.eq.0) cycle
                 endif
                 do i_if1 = 1, n_if1
                    if ( diagonal ) n_if2 = i_if1
                    do i_c1 = 1, n_c1
                       i_meta = (metaindex(quadrupel,i_ir,i_if1,i_c1) - 1) * n_receive_ff + 1
                       do i_if2 = 1, n_if2
                          if ( diagonal .and. i_if1 .eq. i_if2 ) then
                             nn_c2 = i_c1
                          else
                             nn_c2 = n_c2
                          endif
                          do i_c2 = 1, nn_c2
                             call communpack( integralstore(i_meta:i_meta+n_receive_ff-1), &
                                  n_receive_ff, 1, info)
                             if (info .ne. 0) call error_handler( &
                                  "int_send_2cob3c_receive: unpacking 2 center failed")
                             i_meta = i_meta + n_receive_ff
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
            end subroutine unpack_integrals_3c_mem

          end subroutine int_send_2cob3c_receive_mem
          !*************************************************************



          !*************************************************************
          character(len=5) function ua_l_fileindex(i_ua,i_l)
            implicit none
            !------------ Declaration of formal parameters -------------
            integer(kind=i4_kind), intent(in)  :: i_ua,i_l
            !** End of interface ***************************************
            write(ua_l_fileindex,'(i5)') ua_fileindex(i_ua) + i_l        !!!!!!!!!!!
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
            write (charnbr,'(i4)') i_ir        !!!!!!!!!!!!!!!!!!!!!!!!
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
          integer(kind=i4_kind) function metaindex(quadrupel,i_ir,i_ind1,i_exp1,rel_dummy)
            !------------ Declaration of formal parameters -------------
            type(quadrupel_type),  intent(in)  :: quadrupel
            integer(kind=i4_kind), intent(in)  :: i_ir,i_ind1,i_exp1
            logical,optional :: rel_dummy
            !** End of interface ***************************************
            logical :: rel
            rel=.false.
            if(present(rel_dummy)) rel=rel_dummy
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
            end if
          end function metaindex
          !*************************************************************



          !*************************************************************
          subroutine indices(i_meta,i_ir,i_ua1o,i_l1o,i_ind1o,i_exp1o, &
               i_ua2o,i_l2o,i_ind2o,i_exp2o)
            implicit none
            !------------ Declaration of formal parameters -------------
            integer(kind=i4_kind), intent(in)  :: i_meta
            integer(kind=i4_kind), intent(out) :: i_ir,i_ua1o,i_l1o, &
                 i_ind1o,i_exp1o,i_ua2o,i_l2o,i_ind2o,i_exp2o
            !** End of interface ***************************************
            integer(kind=i4_kind)  :: n_ind2, n_l2, n_exp2, i_m,i_ua1,i_ua2,i_l1,i_l2
            integer(kind=i4_kind)  :: i_ind1,i_ind2,i_exp1,i_exp2
            type(unique_atom_type),             pointer :: ua1,ua2
            type(unique_atom_basis_type),       pointer :: uab1,uab2
            type(unique_atom_partner_type),     pointer :: uap1,uap2
            !------------ Executable code ------------------------------
            i_m = 1
            do i_ir = 1,symmetry_data_n_irreps()
               do i_ua1 = 1, N_unique_atoms
                  i_ua1o = i_ua1
                  ua1 => unique_atoms(i_ua1)
                  do i_l1 = 0, ua1%lmax_ob
                     i_l1o=i_l1
                     uab1 => ua1%l_ob(i_l1)
                     uap1 => ua1%symadapt_partner(i_ir,i_l1)
                     do i_ind1 = 1, uap1%N_independent_fcts
                        i_ind1o = i_ind1
                        do i_exp1 = 1, uab1%N_uncontracted_fcts + uab1%N_contracted_fcts
                           i_exp1o = i_exp1
                           do i_ua2 = 1, i_ua1
                              i_ua2o = i_ua2
                              ua2 => unique_atoms(i_ua2)
                              if (i_ua1.eq.i_ua2) then
                                 n_l2 = i_l1
                              else
                                 n_l2 = ua2%lmax_ob
                              endif
                              do i_l2 = 0, n_l2
                                 i_l2o = i_l2
                                 uab2 => ua2%l_ob(i_l2)
                                 uap2 => ua2%symadapt_partner(i_ir,i_l2)
                                 if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 ) then
                                    n_ind2 = i_ind1
                                 else
                                    n_ind2 = uap2%N_independent_fcts
                                 endif
                                 do i_ind2 = 1, n_ind2
                                    i_ind2o=i_ind2
                                    if ( i_ua1 .eq. i_ua2 .and. i_l1 .eq. i_l2 &
                                         .and. i_ind1 .eq. i_ind2 ) then
                                       n_exp2 = i_exp1
                                    else
                                       n_exp2 = uab2%N_uncontracted_fcts + &
                                            uab2%N_contracted_fcts
                                    endif
                                    do i_exp2 = 1, n_exp2
                                       i_exp2o=i_exp2
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
          subroutine int_send_2cob3c_rec_filesizes()
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
          end subroutine int_send_2cob3c_rec_filesizes
          !*************************************************************

          !*************************************************************
          function max_fitindex()
            !  Purpose: returns max. number of fitfcts on this proc.
            !** End of interface ***************************************
            use fit_coeff_module, only: fit_coeff_n_ch
            implicit none
            !------------  Declaration of return value  ---------------
            integer(kind=i4_kind)             :: max_fitindex
            !------------ Declaration of local variables ---------------
            integer(kind=i4_kind) :: first_host,last_host,n_hosts,n_ch,&
                 & n_hosts_larger, n_per_host
            !------------ Executable code ------------------------------

            ! calculate borders
            first_host = 1

            last_host = comm_get_n_processors()
            n_hosts   = last_host - first_host + 1
            if ( n_hosts < 1 ) call error_handler( &
                 "send_2cob3c_module : max_fitindex : wrong n_hosts" )

            n_ch = fit_coeff_n_ch()
            n_hosts_larger = mod(n_ch,n_hosts)
            n_per_host = n_ch / n_hosts

            ! Note:
            ! If number of chargefit functions is not exactly divisible by
            ! the number of proc, then the REMAINDER will be distributed
            ! evenly among the FIRST few procs.
            ! Since the REMAINDER is less then the number of procs., each
            ! of the first procs. will get ONE additional chargefit function.
            if(comm_myindex() > first_host + n_hosts_larger - 1 ) then
               ! I get nothing of the remainder
               max_fitindex = n_per_host
            else
               ! I get something of the remainder
               max_fitindex = n_per_host + 1
            end if

          end function max_fitindex
          !*************************************************************

!--------------- End of module ----------------------------------
end module int_send_2cob3c_module

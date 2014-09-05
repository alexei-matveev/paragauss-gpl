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
module machineparameters_module
!---------------------------------------------------------------
!
!  Purpose: Database to hold machine dependent parameters
!
!  Author: TB
!  Date: 10/95
!
!
!---------------------------------------------------------------------
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!---------------------------------------------------------------------
!
!
! Modification (Please copy before editing)
! Author: AG
! Date:   01
! Description: Variable size Read/Write buffers introduced
!
!---------------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
implicit none
private         ! by default, all names are private
save

!------------ Declaration of Default Parameters -----------------

integer(i4_kind), parameter :: df_machineparameters_veclen    = 128 ! was 2048 on VPP
integer(i4_kind), parameter :: df_machineparameters_rw_buffer = 128 ! 8192 (previous)
logical, parameter :: df_machineparameters_DimCheck = .TRUE.
logical, parameter :: df_machineparameters_PEIGVAL = .FALSE.
real(r8_kind), parameter :: df_machineparameters_checkpreci = 1.0E-6_r8_kind

!== Interrupt end of public interface of module =================

!------------ Declaration of constants and variables ------------

integer(i4_kind), public, protected :: machineparameters_veclen &
     = df_machineparameters_veclen

integer(i4_kind), public, protected :: machineparameters_rw_buffer &
     = df_machineparameters_rw_buffer

integer(i4_kind), public, protected :: lenblk

logical, public, protected :: machineparameters_DimCheck   = &
     df_machineparameters_DimCheck

logical, public, protected :: machineparameters_PEIGVAL   = &
     df_machineparameters_PEIGVAL

real(r8_kind), public, protected :: machineparameters_checkpreci = &
     df_machineparameters_checkpreci
! to switch on internal Checking of dimensions at critical points

!------------ public functions and subroutines ------------------
public :: machineparameters_read, machineparameters_write

public :: machineparameters_bcast!()


  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Declaration of constants and variables ----

integer(i4_kind), private :: machineparameters_ph_veclen = -1 ! UNUSED

namelist /machineparameters/ machineparameters_veclen    , &
                             machineparameters_ph_veclen , &
                             machineparameters_DimCheck  , &
                             machineparameters_PEIGVAL, &
                             machineparameters_checkpreci, &
                             machineparameters_rw_buffer

!---------------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine machineparameters_read
    ! purpose: reads data from input file
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use input_module
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: status, unit
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------
    ! read input
    machineparameters_veclen     = df_machineparameters_veclen
    machineparameters_ph_veclen  = -1
    machineparameters_DimCheck   = df_machineparameters_DimCheck
    machineparameters_PEIGVAL = df_machineparameters_PEIGVAL
    machineparameters_checkpreci = df_machineparameters_checkpreci
    machineparameters_rw_buffer  = df_machineparameters_rw_buffer
    if (input_line_is_namelist("machineparameters")) then

       call input_read_to_intermediate()
       unit = input_intermediate_unit()
       read(unit, nml=machineparameters, iostat=status)
       if (status .gt. 0) call input_error( &
            "machineparameters_read: namelist machineparameters")
    endif

    if ( machineparameters_ph_veclen /= -1 ) then
        WARN('machineparameters_ph_veclen is ignored')
    endif

    if ( machineparameters_veclen .le. 0 .or. &
         machineparameters_veclen .gt. 10000 ) &
       call input_error( &
           "machineparameters_read: machineparameters_veclen out of range")

    if(machineparameters_rw_buffer .le. 1) &
       call input_error( &
           "machineparameters_read: machineparameters_rw_buffer too small")
    lenblk=machineparameters_rw_buffer

  end subroutine machineparameters_read
  !*************************************************************


  !*************************************************************
  subroutine machineparameters_write(iounit)
    !
    ! Purpose: reads data from input file.
    !
    use echo_input_module, only: start, real, intg, flag, stop
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    !** End of interface ***************************************

!   << system dependent namelist format >>
!   integer(kind=i4_kind) :: status
!   write(iounit, nml=machineparameters, iostat=status)
!   if (status .gt. 0) call error_handler( &
!        "machineparameters_write: write failed at namelist machineparameters")
!

     call start("MACHINEPARAMETERS","MACHINEPARAMETERS_WRITE", &
          iounit,operations_echo_input_level)
     call intg("MACHINEPARAMETERS_VECLEN    ", &
                machineparameters_veclen     , df_machineparameters_veclen    )
     call flag("MACHINEPARAMETERS_DIMCHECK  ", &
                machineparameters_DimCheck   , df_machineparameters_DimCheck  )
     call flag("MACHINEPARAMETERS_PEIGVAL  ", &
                machineparameters_PEIGVAL   , df_machineparameters_PEIGVAL  )
     call real("MACHINEPARAMETERS_CHECKPRECI", &
                machineparameters_checkpreci , df_machineparameters_checkpreci)
     call intg("MACHINEPARAMETERS_RW_BUFFER" , &
                machineparameters_rw_buffer  , df_machineparameters_rw_buffer)
     call stop()

  end subroutine machineparameters_write
  !*************************************************************

  subroutine machineparameters_bcast()
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call machineparameters_pack()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call machineparameters_unpack()
    endif
  end subroutine machineparameters_bcast

  !*************************************************************
  subroutine machineparameters_pack()
    ! purpose: packs data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: commpack
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------
    call commpack(machineparameters_veclen,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_pack: machineparameters_veclen")
    call commpack(machineparameters_DimCheck,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_pack: machineparameters_DimCheck")
    call commpack(machineparameters_PEIGVAL,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_pack: machineparameters_PEIGVAL")
    call commpack(machineparameters_checkpreci,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_pack: machineparameters_checkpreci")
    call commpack(lenblk,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_pack: machineparameters_rw_buffer => lenblk)")
  end subroutine machineparameters_pack
  !*************************************************************


  !*************************************************************
  subroutine machineparameters_unpack()
    ! purpose: unpacks data into comm buffer
    !** End of interface ***************************************
    !------------ Modules used -----------------------------------
    use comm_module, only: communpack
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: info
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------
    call communpack(machineparameters_veclen,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_unpack: machineparameters_veclen")
    call communpack(machineparameters_DimCheck,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_unpack: machineparameters_DimCheck")
    call communpack(machineparameters_PEIGVAL,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_unpack: machineparameters_PEIGVAL")
    call communpack(machineparameters_checkpreci,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_unpack: machineparameters_checkpreci")
    call communpack(lenblk,info)
    if (info .ne. 0) call error_handler( &
         "machineparameters_unpack: machineparameters_rw_buffer => lenblk")
   end subroutine machineparameters_unpack
  !*************************************************************


!--------------- End of module ----------------------------------
end module machineparameters_module

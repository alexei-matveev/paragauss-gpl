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
module  iounitadmin_module
  !---------------------------------------------------------------
  !
  !  Purpose: administers unit numbers for read and write to files.
  !           All unit numbers should be demanded with the
  !           get_iounit() function before being used and
  !           returned with the return_iounit(unit) subroutine
  !           when the corresponding file is beeing closed.
  !           The get_iounit() routine returns 0 if no more units
  !           are available.
  !           The number of free unit can be inquired with the
  !           get_nbr_free_iounits() function.
  !
  !           Units for special purpose are defined. Those special
  !           purpose units that should stay open all the program
  !           run are opened and closed to the the
  !           open_special_units() and close_special_units()
  !           subroutines
  !
  !
  !  Module called by: everything opening or closing files
  !
  !
  !  Author: TB
  !  Date: 16.06.95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !
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
  use type_module
  use filename_module, only: max_length => filename_namelengthmax
  use iso_fortran_env, only: stdin_unit => input_unit, &
       stdout_unit => output_unit
  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module =================

  !------------ Declaration of public constants and variables -----

  !
  ! Units reserved for special purpose
  !
  ! The output_unit will  be set to 1 and trace_unit will  be set to 3
  ! after  call  open_special_units(). Both  will  remain negative  on
  ! slaves that do not do not call open_special_units():
  !
  integer, public, protected :: output_unit = -1 ! == 1 after opening
  integer :: trace_unit = -3                     ! == 3 after opening
  public :: stdin_unit                ! == 5
  public :: stdout_unit               ! == 6

#ifdef WITH_EFP
  logical, public :: no_trace_output=.false.
  logical, public :: no_output_unit_output=.false.
#else
  logical, parameter, public :: no_output_unit_output=.false.
#endif
  !------------ public functions and subroutines ------------------
  public get_iounit, return_iounit, get_nbr_free_iounits, &
       openget_iounit, returnclose_iounit
#ifndef FPP_OPTIMIZER
  public :: open_special_units, close_special_units, &
            write_to_output_units, write_to_trace_unit
#endif


  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----

  character(len=12), parameter :: output_filename="output"
  ! intended for regular output
  character(len=12), parameter :: trace_filename="trace_output"
  ! intended for tracing output that can be used to monitor
  ! state of calculation

  !
  ! This determines if subroutine write_to_output_units() writes to
  ! STDOUT in addition to output_unit:
  !
  logical, parameter :: iounitadmin_use_stdout = .true.

  integer, parameter, private :: maxunits=199  ! maximal unit nbr
!!!  integer, parameter, private :: maxunits=99  ! maximal unit nbr
  integer, parameter, private :: startunit=7   ! first free unit
  integer, parameter, private :: maxfreeunits=maxunits-startunit+1
  logical,            private :: unit_free(maxunits)
  logical,            private :: special_units_open = .false.
  integer,            private :: nbr_free_units, first_free_unit
  character(len=max_length), private :: trace_complete_filename

#ifndef WITH_EFP
  logical, parameter :: no_trace_output=.false.
#endif

  !----------- Initialising private variables --------------------
  data unit_free /maxunits * .true./
  data nbr_free_units /maxfreeunits/
  data first_free_unit /startunit/

#ifdef _VPP
#define FPP_HAVE_FLUSH
#endif

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  integer function get_iounit()
    !  Purpose: returns free unit number and reserves it
    !  returns 0 if no more free units
    !** End of interface *****************************************
    !------------ Declaration of variables -----------------------
    integer :: i
    !------------ Executable code --------------------------------
    get_iounit = first_free_unit
    if ( first_free_unit .ne. 0 ) then
       unit_free(first_free_unit) = .false.
       nbr_free_units = nbr_free_units - 1
       do i  = first_free_unit+1,maxunits
          if ( unit_free(i) ) then
             first_free_unit = i
             return
          endif
       enddo
    endif
    first_free_unit = 0

  end function get_iounit
  !*************************************************************


  !*************************************************************
  subroutine return_iounit(unit)
    !  Purpose: frees unit number
    !** End of interface *****************************************
    !------------ Declaration of formal parameters ---------------
    integer, intent(in ) :: unit
    !------------ Executable code --------------------------------
    if ( unit .ge. startunit .and. &
         .not. unit_free(unit)      ) then
       unit_free(unit) = .true.
       nbr_free_units = nbr_free_units + 1
       if ( unit .lt. first_free_unit .or. &
            first_free_unit .eq. 0         ) then
          first_free_unit = unit
       endif
    endif
  end subroutine return_iounit
  !*************************************************************


  !*************************************************************
  integer function get_nbr_free_iounits()
    !  Purpose: returns number of free units
    !** End of interface *****************************************
    !------------ Executable code --------------------------------
    get_nbr_free_iounits = nbr_free_units
    return
  end function get_nbr_free_iounits
  !*************************************************************


  !*************************************************************
  integer function openget_iounit(file,status,access, &
       form,recl,blank,position,action,delim,pad)
    !  Purpose: returns free unit number and reserves it,
    !  and opens file, performing error handling. All arguments
    !  to the standard f90 open commands are possible in the
    !  standard syntax.
    !------------ Declaration of formal parameters ---------------
    character(len=*),           intent(in)  :: file
    character(len=*), optional, intent(in)  :: status
    character(len=*), optional, intent(in)  :: access
    character(len=*), optional, intent(in)  :: form
    integer,          optional, intent(in)  :: recl
    character(len=*), optional, intent(in)  :: blank
    character(len=*), optional, intent(in)  :: position
    character(len=*), optional, intent(in)  :: action
    character(len=*), optional, intent(in)  :: delim
    character(len=*), optional, intent(in)  :: pad
    !** End of interface *****************************************
    character(len=20) :: istatus, iaccess, iform, iblank, &
         & iposition,iaction,idelim,ipad
    integer irecl
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    DPRINT 'ioua::openget_iounit: entered, file=',file,'<<'
    if(present(status)) then
       istatus = trim(status)
    else
       istatus = 'unknown'
    endif
    if(present(access)) then
       iaccess = trim(access)
    else
       iaccess = 'sequential'
    endif
    if(present(form)) then
       iform = trim(form)
    else
       if (iaccess .eq. 'direct' .or. iaccess .eq. 'DIRECT' ) then
          iform =       'unformatted'
       else
          iform =       'formatted'
       endif
    endif
    if(present(recl)) then
       irecl = recl
    else
       irecl = 65536
    endif
    if(present(blank)) then
       iblank = blank
    else
       iblank = 'null'
    endif
    if(present(position)) then
       iposition = trim(position)
    else
       iposition = 'asis'
    endif
    if(present(action)) then
       iaction = trim(action)
    else
       iaction = 'readwrite'
    endif
    if(present(delim)) then
       idelim = trim(delim)
    else
       idelim = 'none'
    endif
    if(present(pad)) then
       ipad = trim(pad)
    else
       ipad = 'yes'
    endif

    openget_iounit =  get_iounit()
    if ( openget_iounit .eq. 0 ) &
         call error_handler("openget_iounit: no more units")

    if ( iaccess .eq. 'direct' .or. iaccess .eq. 'DIRECT' ) then
       if ( iform .eq. 'formatted' .or. iform .eq. 'FORMATTED' ) then
          open(openget_iounit,err=888,file=trim(file),status=trim(istatus), &
               access=trim(iaccess),form=trim(iform),recl=irecl,blank=trim(iblank), &
               action=trim(iaction),delim=trim(idelim),pad=trim(ipad))
       else
          open(openget_iounit,err=888,file=trim(file),status=trim(istatus), &
               access=trim(iaccess),form=trim(iform),recl=irecl,action=trim(iaction))
       endif
    else
       if ( iform .eq. 'formatted' .or. iform .eq. 'FORMATTED' ) then
          open(openget_iounit,err=888,file=trim(file),status=trim(istatus), &
!AG            access=trim(iaccess),form=trim(iform),recl=irecl,blank=trim(iblank), &
               access=trim(iaccess),form=trim(iform),           blank=trim(iblank), &
               position=trim(iposition),action=trim(iaction), &
               delim=trim(idelim),pad=trim(ipad))
       else
#ifndef _DEC
          open(openget_iounit,err=888,file=trim(file),status=trim(istatus), &
               access=trim(iaccess),form=trim(iform),recl=irecl, &
               position=trim(iposition),action=trim(iaction))
#else
          open(openget_iounit,file=trim(file),status=trim(istatus), &
               access=trim(iaccess),form=trim(iform),recl=irecl, &
               position=trim(iposition),action=trim(iaction))
#endif
       endif
    endif
    DPRINT 'ioua::openget_iounit: unit(',openget_iounit,')=',file,'<<'
    return
888 call error_handler("openget_iounit: failed for file "//file//" .")
  end function openget_iounit
  !*************************************************************


  !*************************************************************
  subroutine returnclose_iounit(unit,status)
    !  Purpose:  frees unit number and closes unit,
    !  performing error handling
    !------------ Declaration of formal parameters ---------------
    integer,                    intent(in) :: unit
    character(len=*), optional, intent(in) :: status
    !** End of interface *****************************************
    character(len=20) :: istatus
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------
    DPRINT 'ioua::returnclose_iounit: close unit(',unit,')'
    istatus='keep'
    if(present(status)) istatus=status
    call return_iounit(unit)
    close(unit,err=999,status=trim(istatus))
    return
999 call error_handler( "returnclose_iounit: closing unit failed")
  end subroutine returnclose_iounit
  !*************************************************************


#ifndef FPP_OPTIMIZER
  !*************************************************************
  subroutine open_special_units (rank)
    !  Purpose: special purpose units that should stay open all
    !           the program run are opened
    use filename_module, only: outfile
    implicit none
    integer, intent (in) :: rank
    !** End of interface *****************************************

    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Declaration of variables -----------------------
    integer      :: iostat
    character(len=120) :: message
    logical, parameter :: append_local = .false.
    !------------ Executable code --------------------------------

    ! We  used  to have  slaves  do some  output  too.   At many,  but
    ! possibly not all  places IO to output_unit will  only be done if
    ! the unit  number is positive.   An early return here  makes sure
    ! the  unit numbers on  all slave  workers remain  negative.  This
    ! helps enforcing  the master-only IO.  Comment this,  if you want
    ! all slave to do IO too:
    !
    if (rank /= 0) return

    ! Output and trace unit numbers remains negative on slaves that do
    ! not call this subroutine:
    output_unit = 1
    trace_unit = 3

#ifdef FPP_DEBUG
# define MyID "[??]"
#endif
    DPRINT MyID,'open_special_units: entered; output_unit=',output_unit
    DPRINT MyID,'open_special_units: output_filename=',trim(output_filename),"<"

    if (append_local) then
       open(output_unit,file=trim(outfile(output_filename)), &
            iostat=iostat,position="append",status="old")
    else
       open(output_unit,file=trim(outfile(output_filename)), &
            iostat=iostat,status = "replace")
    endif

    if (iostat.ne.0) then
       message = "open_special_units: opening failed for file "// &
            trim(outfile(output_filename))
       call error_handler(message)
    endif

    special_units_open = .true.

    if (append_local) return

    !
    ! Set private global var, the path may be constructed differently for
    ! master and slaves:
    !
    trace_complete_filename = outfile(trace_filename)

    open(trace_unit,file=trim(trace_complete_filename), &
         status = "replace", iostat=iostat)
    if (iostat.ne.0) call error_handler( &
         "open_special_units: opening failed for " // &
         trim(trace_complete_filename) )

#ifndef FPP_HAVE_FLUSH
    close(trace_unit, iostat=iostat)
    if (iostat.ne.0) call error_handler( &
         "open_special_units: closing failed for " // &
         trim(adjustl(trace_complete_filename)) )
#endif
  end subroutine open_special_units
  !*************************************************************


  !*************************************************************
  subroutine close_special_units ()
    !
    ! Purpose:  special purpose units  that should  stay open  all the
    ! program run are closed.
    !
    implicit none
    !** End of interface *****************************************

    external error_handler

    integer :: iostat

    if (.not. special_units_open) return

    close(output_unit,iostat=iostat)
    if (iostat.ne.0) call error_handler( &
         "close_special_units: closing failed for " // &
         trim(adjustl(output_filename)) )

#ifdef FPP_HAVE_FLUSH
    close(trace_unit, iostat=iostat)
    if (iostat.ne.0) call error_handler( &
         "close_special_units: closing failed for file " // &
         trim(adjustl(trace_complete_filename)) )
#endif

    output_unit = -1
    trace_unit = -3
    special_units_open = .false.
  end subroutine close_special_units
  !*************************************************************


  !*************************************************************
  subroutine write_to_output_units(message,inte,re)
    !  Purpose: writes message to output, debug op and trace op
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: message
    integer(kind=i4_kind),intent(in), optional :: inte
    real(kind=r8_kind),intent(in), optional    :: re
    !** End of interface *****************************************
    if (no_output_unit_output) goto 100
    if ( special_units_open ) then
       if(present(inte).and..not.present(re)) then
          write(output_unit,*) message,inte
       elseif(present(re).and..not.present(inte)) then
          write(output_unit,*) message,re
       elseif(present(inte).and.present(re)) then
          write(output_unit,*) message,inte,re
       else
          write(output_unit,*) message
       endif
    endif

100 continue
    if (iounitadmin_use_stdout) then
       if (present(inte).and..not.present(re)) then
          print *, message,inte
       elseif(present(re).and..not.present(inte)) then
          print *, message,re
       elseif(present(re).and.present(inte)) then
          print *, message,inte,re
       else
          print *, message
       endif
    endif
  end subroutine write_to_output_units
  !*************************************************************


  !*************************************************************
  subroutine write_to_trace_unit(msg, inte, real)
    !  Purpose: writes message to end of trace file, opening and
    !  and closing the file
    use comm, only: comm_rank
    use time_module, only: clktime
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: msg
    integer(kind=i4_kind), intent(in), optional :: inte
    real(kind=r8_kind), intent(in), optional :: real
    !** End of interface *****************************************

    logical, save :: warned = .false.
    integer(i4_kind)           :: stat
    real(r8_kind)              :: time
    character(len=11)          :: prefix
    character(len=12+len(msg)) :: message

    !
    ! This  is  used in  a  single  file, efp_only_opt_module.f90,  to
    ! temporarily disable trace output. Dont ask me why.
    !
    if (no_trace_output) return

    !
    ! Some platforms have problems  with many workers writing to files
    ! simultaneousely, so trace output may be only enabled for master:
    !
    if ( trace_unit <= 0 ) then
       if ( .not. warned ) then
          print *, "write_to_trace_unit(", msg, "... ) on rank", comm_rank()
          print *, "WARNING: writing to trace file is disabled for workers"
          print *, "         that did not call open_special_units()."
          warned = .true.
       endif
       RETURN ! *** RETURN POINT ***
    endif
    ASSERT(trace_unit>0)

    ! prepend current time to each line of the trace_output:
    time = clktime()
    write(prefix,'("[",F9.2,"]")') time

    message = prefix // " " // msg
    ! will look like this:
    !
    ! [    13.99] Entering SCF part
    !

#ifndef FPP_HAVE_FLUSH
    open(trace_unit, file=trim(trace_complete_filename), &
         position="append", action="write", iostat=stat)
    if( stat /= 0 ) then
        WARN("Error opening trace_unit")
        goto 9999 ! abort
    endif
#endif

    !
    ! F20.11 is also used in main_scf() for SCF line entries:
    !
    if (present(inte).and..not.present(real)) then

       write(trace_unit, '(A, I8)', iostat=stat) message, inte

    elseif(present(real).and..not.present(inte)) then

       write(trace_unit, '(A, F20.11)', iostat=stat) message, real

    elseif(present(real).and.present(inte)) then

       write(trace_unit, '(A, I8, F20.11)', iostat=stat) message, inte, real

    else
       write(trace_unit, '(A)', iostat=stat) message
    endif
    if( stat /= 0 ) goto 9999 ! abort

#ifndef FPP_HAVE_FLUSH
    close(trace_unit, iostat=stat)
    if( stat /= 0 ) then
        WARN("Error closing trace_unit")
        goto 9999 ! abort
    endif
#else
#ifdef FPP_HAVE_F2003_FLUSH
    flush(trace_unit, iostat=stat)
    if( stat /= 0 ) then
        WARN("Error flushing trace_unit")
        goto 9999 ! abort
    endif
#else
    call flush(trace_unit)
#endif
#endif

    return

9999 CONTINUE
    print *, "iostat=", stat
    print *, "msg=", message
    if( present(inte) ) print *, "inte=", inte
    if( present(real) ) print *, "real=", real
    ABORT('Error: write_to_trace_unit: failed')
  end subroutine write_to_trace_unit
#endif
  !*************************************************************

  !--------------- End of module ----------------------------------
end module iounitadmin_module

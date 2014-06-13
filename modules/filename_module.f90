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
module  filename_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: defining of directoties for input and output:
  !
  !  Defines:
  !
  !     inpfile(name), outfile(name), tmpfile(name),
  !     data_dir, sharedfs_dir,
  !     input_name
  !
  ! The names are read in from environment variables
  ! TTFSTMP, TTFSOUT, TTFSDATADIR, TTFSSTART, TTFSOUTPUTDIR,
  ! TTFSINPUT,TTFSSHAREDFS, TTFSRECOVERDIR
  !
  ! call filename_setup() to read in dir names
  ! and use filename_bcast for sending
  ! information to slaves
  !
  ! in case a parralel file system is used, TTFSSHAREDFS must be set and
  ! contain the name of a directory that contains the directories tmp_dir
  ! for the seperate processors named according to their processor ids
  ! (for instance master: setenv TTFSTMP $TTFSSHAREDFS/1 )
  !
  !  Author: TB
  !  Date: 9/95
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   11/96
  ! Description: added TTFSOUTPUTDIR stored in output_dir
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   2/97
  ! Description: added TTFSSHAREDFS stored in sharedfs_dir
  !              to pupport parallel file systems
  !
  ! Modification
  ! Author: HH
  ! Date:   5/99
  ! Description: Added two new environment variables
  ! TTFSRESPDIR    -> Final output from response module, i.e.
  !                   binary interface-integral-files for response
  !                   program package RESTDD
  !
  !                   UPDATE: use TTFSTMP instead:
  !                   For analytical 3-index Coulomb-integrals in
  !                   primitive basis that are created in the integral part
  !                   (int_send_2cob3c_module.f90p) and read in the
  !                   response module
  !                   Note: these files may become very large (~30GB)
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
#include "def.h"
#ifndef MAX_PATH
# define MAX_PATH 100
#endif
  implicit none
  private         ! by default, all names are private
  save
  !== Interrupt end of public interface of module =================

  integer, parameter, private :: max_path = MAX_PATH

  !------------ Declaration of public variables -------------------
  integer, parameter, public :: filename_namelengthmax = max_path

  character (len=max_path), public, protected :: data_dir, &
       sharedfs_dir, resp_dir

  character (len=max_path), public, protected :: input_name

  logical, public, protected :: filesystem_is_parallel = .FALSE.

  !------------ public functions and subroutines ------------------

  interface filename_env ! name for Public Relations
     module procedure env
  end interface
  public :: filename_env

  public :: filename_setup
  public :: filename_shutdown
  public :: filename_set_input_dir

  public :: filename_tmpdir

  !
  ! Your  should assume  that tmpfile(name)  will return  a path  to a
  ! node-local directory  potentially different for  every worker. You
  ! can assume that the outfile(name)  on the master will be delivered
  ! to the user but the very  same funciton on other wokers may return
  ! a path to some scratch file.
  !
  public :: inpfile ! (name) -> path, returns full path to an input file
  public :: outfile ! (name) -> path, returns full path to an output file
  public :: tmpfile ! (name) -> path, returns full path to a temp file
  public :: recfile ! (name) -> path, returns full path to a recover file


  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private variables ------------------

  integer, private :: my_hostindex = -1 ! initialized in filename_setup()

  character (len=max_path), private :: tmp_base ! rank non-specific
  character (len=max_path), private :: tmp_dir ! rank-specific
  character (len=max_path), private :: output_dir ! rank-specific
  character (len=max_path), private :: input_dir ! shared by all workers
  character (len=max_path), private :: recover_dir ! shared by all workers


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function tmpfile (name) result (path)
    !
    !  Returns full path for within tmp_dir unique for each worker
    !
    !  NOTE: the path is constructed as
    !
    !           path == TMP_BASE/xxxx-name
    !
    !         with xxxx being the zero-padded base-1 worker ID
    !
    implicit none
    character(len=*), intent(in)          :: name
    character(len=max_path) :: path ! result
    ! *** end of interface ***

    character(len=4) :: xxxx

    ! call filename_setup() first:
    ASSERT(my_hostindex>=0)

    ! FIXME: needs wider char buf:
    ASSERT(my_hostindex<10000)

    write(xxxx, '(I4.4)') my_hostindex

    !
    ! tmp_dir is supposed to be unique for each rank:
    !
    path = trim(tmp_base) // '/' // xxxx // "-" // trim(adjustl(name))
  end function tmpfile

  function outfile (name) result (path)
    !
    ! Returns filename with output dir (of this process) prepended.
    !
    ! NOTE:  historically slaves  produce output  files with  the same
    ! names  as  the  master  assuming their  output  directories  are
    ! private.   On the  other  hand  the output  of  slaves is  never
    ! (rarely) examined.   This sub redirects the output  of slaves to
    ! their rank-specific  files in tmp_dir  and the output  of master
    ! goes to output_dir.
    !
    implicit none
    character(len=*), intent(in)          :: name
    character(len=max_path) :: path ! result
    ! *** end of interface ***

    ! call filename_setup() first:
    ASSERT(my_hostindex>=0)

    if ( my_hostindex == 1 ) then
        !
        ! Master writes to final output dir:
        !
        path = trim (output_dir) // '/' // trim (adjustl (name))
    else
        !
        ! Rank-specific output of slaves goes to temp dir:
        !
        path = tmpfile (name)
    endif
  end function outfile

  function inpfile (name) result (path)
    !
    !  Returns full path for within input_dir
    !
    implicit none
    character(len=*), intent(in)          :: name
    character(len=max_path) :: path ! result
    ! *** end of interface ***

    !
    ! input_dir is supposed to be shared (if slaves ever need it)
    !
    path = trim (input_dir) // '/' // trim (adjustl (name))
  end function inpfile


  function recfile (name) result (path)
    implicit none
    character (len=*), intent(in) :: name
    character (len=max_path) :: path ! result
    ! *** end of interface ***

    !
    ! Recover files  may need to be  saved between runs. Let  it be in
    ! the user directory  and point all workers to  the same files. By
    ! default recover_dir is the same as input_dir:
    !
    path = trim (recover_dir) // '/' // trim (adjustl (name))
  end function recfile


  subroutine filename_set_input_dir(dir)
    implicit none
    character(len=*), intent(in) :: dir
    ! *** end of interface ***

    input_dir = dir
  end subroutine filename_set_input_dir

  subroutine filename_setup()
    !
    ! Runs in parallel context:
    !
    use comm, only: comm_rank
    implicit none
    ! *** end of interface ***

    integer :: rank

    rank = comm_rank()

    !
    ! Do not abuse global module state:
    !
    my_hostindex = rank + 1

    !
    ! Remote nodes may not have proper environment,
    ! so let the master process do it ...
    !
    if ( rank == 0 ) then
        call guess_environment()
    endif

    !
    ! ... and broadcast results to all workers:
    !
    call bcast()

    !
    ! This sets worker specific tmp_dir, traditionally
    ! using base-1 indices here. FIXME: should we?
    !
    tmp_dir = append_hostindex(tmp_base, rank + 1)
  end subroutine filename_setup

  !*************************************************************
  subroutine guess_environment()
    !  Purpose: reads environment variables TTFSTMP, TTFSOUT ,
    !  TTFSDATADIR TTFSSTART and TTFSOUTPUTDIR to determine filenames
    implicit none
    !** End of interface *****************************************

    !
    ! Ask   environment  variable   "TTFSTMP"  that   should  describe
    ! directory  for scratch  data. If  "TTFSSHAREDFS" is  set, assume
    ! tmp_dir is shared among all nodes:
    !
    filesystem_is_parallel = .false.
    sharedfs_dir = "/do/not/use"

    if (env ("TTFSTMP", tmp_base)) then
       ! fine, nothing to be done
    else if (env ("TTFSSHAREDFS", sharedfs_dir)) then
       filesystem_is_parallel = .true.
       tmp_base = sharedfs_dir
    else
       WARN ("TTFSTMP not set, using ." )
       tmp_base = "."
    endif

    ! Ask  environment variable  "TTFSINPUTDIR"  that should  describe
    ! directory where input data are originally stored
    if (.not. env ("TTFSINPUTDIR", input_dir)) then
      ! WARN ("TTFSINPUTDIR not set, using .")
      input_dir = "."
    endif

    ! Ask  environment variable  "TTFSOUTPUTDIR" that  should describe
    ! directory where output data will finally be copied to
    if (.not. env ("TTFSOUTPUTDIR", output_dir)) then
      ! WARN ("TTFSOUTPUTDIR not set, using ." )
      output_dir = "."
    endif

    ! Ask  environment  variable  "TTFSDATADIR" that  should  describe
    ! directory for intermediate files living longer than the job
    if (.not. env ("TTFSDATADIR", data_dir)) then
      ! WARN ("TTFSDATADIR not set, using "//trim(tmp_base))
      data_dir = tmp_base ! FIXME: or rather output_dir?
    endif

    ! Ask environment  variable "TTFSRECOVERDIR" that  should describe
    ! directory where scf recover files are stored
    if (.not. env ("TTFSRECOVERDIR", recover_dir)) then
      ! WARN ("TTFSRECOVERDIR not set, using "//trim(input_dir))
      recover_dir = input_dir
    endif

    ! Ask environment variable "TTFSINPUT" that is name of input file
    if (.not. env ("TTFSINPUT", input_name)) then
      ! WARN ("TTFSINPUT not set, using 'input'")
      input_name = "input"
    endif

    ! Ask  environment  variable "TTFSRESPDIR"  ->  Final output  from
    ! response  module,   i.e.   binary  interface-integral-files  for
    ! response program package RESTDD
    if (.not. env ("TTFSRESPDIR", resp_dir)) then
      ! WARN ("TTFSRESPDIR not set, using "//trim(tmp_base))
      resp_dir = tmp_base
    endif
  end subroutine guess_environment
  !*************************************************************


  !*************************************************************
  function append_hostindex(i_name,hostindex) result(o_name)
    !  Purpose: returns name associated with processor
    !  with index hostindex
    !------------ Declaration of formal parameters ---------------
    character(len=*),intent(in)           :: i_name
    integer, intent(in)                   :: hostindex
    character(len=max_path) :: o_name ! <<< result
    !** End of interface *****************************************

    character(len=4) :: char

    write(char,'(I4)') hostindex
    char = adjustl(char)

    o_name = repeat(" ",max_path)
    o_name = trim(i_name) //"/"//trim(char) ! // "/"
  end function append_hostindex
  !*************************************************************

#ifndef NO_COMM

  subroutine bcast()
    use comm, only: comm_bcast
    implicit none
    ! *** end of interface ***

    call comm_bcast(filesystem_is_parallel)
    call comm_bcast(sharedfs_dir)
    call comm_bcast(tmp_base)
    call comm_bcast(data_dir)
    call comm_bcast(output_dir)
    call comm_bcast(input_dir)
    call comm_bcast(input_name)
    call comm_bcast(resp_dir)
  end subroutine bcast

#endif

  character(len=max_path) function filename_tmpdir(hostindex)
    !  Purpose: returns name of tmpdir associated with processor
    !  with index hostindex
    !------------ Declaration of formal parameters ---------------
    integer, intent(in) :: hostindex
    !** End of interface *****************************************

    ! FIXME: one should not even try to access directories of
    !        other workers as it is not always possible.
    !        Maybe we should implement tmpfile(name, rank) instead?
    ASSERT(filesystem_is_parallel)

    filename_tmpdir = append_hostindex(tmp_base, hostindex)
  end function filename_tmpdir

  subroutine filename_shutdown()
    ! Purpose: deallocates the private variables
    !** End of interface *****************************************

    ! FIXME: should we do this?
    ! my_hostindex = -1
   end subroutine filename_shutdown
  !*************************************************************

  recursive function env (var, val) result (def)
    !
    !
    implicit none
    character (len=*), intent (in) :: var
    character (len=*), intent (inout) :: val
    logical :: def              ! result
    ! *** end of interface **

    integer :: stat, length

#ifdef _ITANIUM_NSK
#warning "recursive call to env(var, val), really?"
    !
    ! FIXME: even in this case it is assumed that getting
    !        $PWD from environment should work (see below).
    !        Why should't it work for one (or a few) more
    !        variables like $TTFSTMP?
    !
    if ( var /= "PWD" ) then
      ! read TTFS-specific env. vars from file
      def = get_ttfs_env(var,val)
      ! better GET RID of TTFS-specific vars completely!
      return
    endif
#endif

    !
    ! F2003 feature, may not be available with older compilers:
    !
    call GET_ENVIRONMENT_VARIABLE(var, val, length, stat)

    ! make compiler happy, always set return value:
    def = .false.

    if ( stat == 0 ) then
      ! defined and fits into val(:)
      def = .true.
      if (length > 0) then
         if ( val(length:length) .eq. '/' ) then
            WARN('DONT USE SLASH AT THE END OF DIRS!')
         endif
      endif
    else if ( stat == 1 ) then
      ! variable undefined, val(:) should be filled with blanks,
      ! nothing to do, return value already set to false ...
    else if ( stat == -1 ) then
      ! value to long:
      ABORT("too long: $"//var//"="//val)
    else if ( stat == 2 ) then
      ! platform does not support environment
      ABORT("no environment")
    else
      ABORT("should not occur")
    endif
  end function env

#ifdef _ITANIUM_NSK
  function get_ttfs_env(var,val) result(def)
    use type_module
    implicit none
    character(LEN=*),intent(in)    :: var
    character(LEN=*),intent(inout) :: val
    logical                        :: def !<<<result
    ! *** end of interface **
    character(len=max_path) :: work_dir
    integer :: namelength,i,nn
    integer(kind=i4_kind) :: ttfs_env
    character(len=200) :: input_line

    if(.not. env("PWD",work_dir)) &
        call error_handler( &
               "filename_setup: Error: environment variable PWD must be set" )

    open(ttfs_env,file=trim(work_dir)//'/TTFS_ENV', &
        form='formatted', status='unknown')
    do
       read(ttfs_env,'(a200)',end=200,err=200) input_line
       if(index(input_line,trim(var)) == 1) exit
    enddo
    close(ttfs_env)

    nn=index(input_line,"=")
    val=input_line(nn+1:)
    namelength = LEN(val)
    do i=1,LEN(val)
       if (val(i:i) .eq. ' ') then
          namelength = i-1
          exit
       endif
    enddo
    if ( namelength .eq. 0 ) then
       def = .false.
    else
       def = .true.
       if ( val(namelength:namelength) .eq. '/' ) then
          WARN('DONT USE SLASH AT THE END OF DIRS!')
       endif
    endif
    return
200 def = .false.
    close(ttfs_env)
  end function get_ttfs_env
#endif

  !--------------- End of module ----------------------------------
end module filename_module

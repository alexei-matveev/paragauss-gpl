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
subroutine epe_lattice_optimization()
  !----------------------------------------------------------------
  !
  !  Purpose: to make EPE lattice optimization
  !
  !
  !  Subroutine called by: main_master
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
#include "def.h"
  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use iounitadmin_module
  use comm_module
  use epecom_module ,  only: periodic_optimization, output_epe, &
                             epedata_dir, pg_cfield_calculated, &
                             epe_input_dir, do_print, qm_interfaced_mode, &
                             n_iterations,read_configuration, epeit_count, n_centres_of_generation
  use main_epe_module, only: main_epe,start_epe,finish_epe, &
                             epe_send_init_to_slave, &
                             epe_send_finish_to_slave
  implicit none
  !------------ Declaration of subroutines ------------------------
  external error_handler
  !------------ Declaration of local variables --------------------
  integer (kind=i4_kind):: conv_unit, namelength,i, i_slash
  character(len=500)    :: file_name, sys_command
#ifndef NEW_EPE
  integer (kind=i4_kind):: input_ep, local_file
  character(len=500)    :: first_line
#endif
  logical, parameter :: dealloc=.true.
  logical, parameter :: no_dealloc=.false.
  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------

   call write_to_output_units("=======================================")
   call write_to_output_units("=                                     =")
   call write_to_output_units("=           Start EPE program         =")
   call write_to_output_units("=                                     =")
   call write_to_output_units("=======================================")

  pg_cfield_calculated=.false.
  !!! rm files from TTFSDATADIR else they will substitued 
  !!! those from epedata_dir
  if(.not.senv("TTFSDATADIR",epedata_dir,slash=.true.)) &
   call error_handler("epe_lattice_optimization:  TTFSDATADIR not defined")   
  
  print*,'PG epedata_dir----------',trim(epedata_dir)
!!$  call system("ls -lrt "//trim(epedata_dir))
!  call system("rm "//trim(epedata_dir)//"/epe.pcs")
!  call system("rm "//trim(epedata_dir)//"/epe.pcc")
!  print *,' SEE that TTFSDATADIR do not contains EPE files'
!  call system("ls -lrt "//trim(epedata_dir))

!AG           make epe_script, copying epe-files => /scratch/...
  if(.not.senv("TTFSINPUTDIR", epe_input_dir,slash=.false.)) &
       call error_handler("epe_lattice_optimization: TTFSINPUTDIR not defined")
        print*,'EPE_INPUT_DIR ',trim(epe_input_dir) 
  namelength =100
  do i=1,100
     if (epe_input_dir(i:i).eq.' ') then
        namelength = i-1
        exit
     endif
  enddo
   do i=1,namelength
     if(epe_input_dir(namelength+1-i:namelength+1-i) .eq. '/') then
      i_slash=namelength+1-i
      exit
     end if
   end do

   file_name = epe_input_dir(i_slash+1:namelength)
   print*, 'epe_input_file ', trim(epe_input_dir)//'/'//file_name

  namelength =100
  do i=1,100
     if (epedata_dir(i:i).eq.' ') then
        namelength = i-1
        exit
     endif
  enddo
  if ( epedata_dir(namelength:namelength) .ne. '/' ) then
     namelength = namelength + 1
     epedata_dir(namelength:namelength) = '/'
  endif
!       call make_out2in()

  if(periodic_optimization) then 
     conv_unit=openget_iounit (trim(epedata_dir)//'conv',status='unknown')
     call returnclose_iounit (conv_unit)
  endif ! periodic_optimization

#ifdef NEW_EPE
  call read_epeinput()
#else
  print*,' open epe.input'
  input_epe=openget_iounit(trim(epe_input_dir)//'/epe.input', &
       form='formatted', status='old')
  read(input_epe,'(a)') first_line
  call returnclose_iounit(input_epe)
  print*,' done'

  if (index(first_line,"work_options") /= 0 .or. &
       index(first_line,"WORK_OPTIONS") /= 0) then
     call read_epe_input
  else if (index(first_line,"tasks") /= 0 .or. &
       index(first_line,"TASKS") /= 0) then
     call read_epe_input_new
  else
     call error_handler("Check your EPE input file. I think it is wrong")
  endif
#endif

#ifdef NEW_EPE
  do_print=.true.
  write(output_epe,*) "****************** OPTIMIZATION OF EPE ENVIRONMENT *******************"
  write(output_epe,*) ""
#endif
  call start_epe
    if( comm_parallel() ) then 
       call write_to_output_units("lattice_opt => epe_send_init_to_slave")
          call epe_send_init_to_slave()
       call write_to_output_units("lattice_opt <= epe_send_init_to_slave")
    end if
  call main_epe
 !print*, ' main_epe done 1'
#ifdef NEW_EPE
  call finish_epe(no_dealloc)
#else
  call finish_epe(dealloc)
#endif
  if(comm_parallel() ) call epe_send_finish_to_slave()

#ifdef NEW_EPE
  write(output_epe,*) "*********** EPE.PCR, EPE.PCS, EPE.PCC and REG_REFERENCE **************"
  write(output_epe,*) ""
  read_configuration=.true.
  do_print=.false.
  qm_interfaced_mode=.true.
  n_iterations=0
  epeit_count=0
  n_centres_of_generation=0
  call start_epe
  if( comm_parallel() ) then 
     call epe_send_init_to_slave()
  end if
  call main_epe
  call finish_epe(dealloc)
  if(comm_parallel() ) call epe_send_finish_to_slave()
#endif

  call returnclose_iounit(output_epe)
 print*,'output_epe returned'
!!$  call returnclose_iounit(conv_unit)
  conv_unit=openget_iounit (trim(epedata_dir)//'conv',status='unknown')
  call returnclose_iounit(conv_unit,status='delete')

!                      backward copying epe-files
  sys_command=trim(epe_input_dir)//'/./epe_script ' &
            //trim(epe_input_dir)//' '//trim(file_name)//' finish'
  print*, 'command  ', trim(sys_command)
!  call system( trim(sys_command) )

  call write_to_output_units("=======================================")
  call write_to_output_units("=                                     =")
  call write_to_output_units("=   Finish EPE lattice optimization   =")
  call write_to_output_units("=                                     =")
  call write_to_output_units("=======================================")

contains

  recursive function senv(var,val,slash) result(def)
    implicit none
    character(LEN=*),intent(in)    :: var
    character(LEN=*),intent(inout) :: val
    logical,optional,intent(in)    :: slash
    logical                        :: def !<<<result
    ! *** end of interface **

    integer :: namelength,i
    logical :: slash_

    slash_ = .false.
    if(present(slash)) slash_ = slash

#ifdef _ITANIUM_NSK
    if ( var /= "PWD" ) then
      ! read TTFS-specific env. vars from file
      def = get_ttfs_env(var,val)
      ! better GET RID of TTFS-specific vars completely!
      return
    endif
#endif
    call GETENV(var,val)

!!$    print *,'env: var=',var,'<'
!!$    print *,'env: val=',val,'<'

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
       if ( slash_ .and. val(namelength:namelength) .ne. '/' ) then
          namelength = namelength + 1
          val(namelength:namelength) = '/'
       endif
    endif
  end function senv

#ifdef _ITANIUM_NSK
  function get_ttfs_env(var,val) result(def)
    use type_module
    use filename_module, only: filename_namelengthmax
    implicit none
    character(LEN=*),intent(in)    :: var
    character(LEN=*),intent(inout) :: val
    logical                        :: def !<<<result
    ! *** end of interface **
    character(len=filename_namelengthmax) :: work_dir
    integer :: namelength,i,nn
    integer(kind=i4_kind) :: ttfs_env
    character(len=200) :: input_line

    if(.not. senv("PWD",work_dir)) &
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
end subroutine epe_lattice_optimization

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
subroutine epe_driver
!
  !===================================================================
  ! End of public interface of module
  !===================================================================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
! aaaaaaa
!---------------------------------------------------------------------

use type_module
use comm_module
use epecom_module, only: epedata_dir, n_iterations, ex_pgdata, &
     epe_input_dir
use iounitadmin_module
use filename_module,only: data_dir, inpfile
use main_epe_module, only: main_epe,start_epe, &
                           epe_send_init_to_slave

implicit none
 integer (kind=i4_kind):: conv_unit, namelength,i
#ifndef NEW_EPE
 integer (kind=i4_kind):: input_epe
 character(len=500)    :: first_line
#endif

        call getenv("TTFSDATADIR",epedata_dir)
        print*,'epe_driver: epedata_dir= ',epedata_dir
        namelength =100
        do i=1,100
           if (epedata_dir(i:i).eq.' ') then
              namelength = i-1
              exit
           endif
        enddo
        epedata_dir=data_dir
  if ( epedata_dir(namelength:namelength) .ne. '/' ) then
       namelength = namelength + 1
       epedata_dir(namelength:namelength) = '/'
    endif

        ! FIXME: quick hack, inpfile() with an empty file name returns directory
        !        with slash appended, alternatively, make input_dir public again.
        ! epe_input_dir = input_dir
        epe_input_dir = inpfile('')
!!$     print*,'epe_driver: epe_input_dir ',epe_input_dir

#ifdef NEW_EPE
  call read_epeinput()
#else
  print*,'open '//trim(epe_input_dir)//'/epe.input',get_nbr_free_iounits()
  input_epe=openget_iounit(trim(epe_input_dir)//'/epe.input', &
       form='formatted', status='old')
  print*,'done',input_epe
  read(input_epe,'(a)') first_line
  call returnclose_iounit(input_epe)

  if (index(first_line,"work_options") /= 0 .or. &
       index(first_line,"WORK_OPTIONS") /= 0) then
     call read_epe_input
  else if (index(first_line,"tasks") /= 0 .or. &
       index(first_line,"TASKS") /= 0) then
     print*,'read_epe_input_new'
     call read_epe_input_new
  else
     call error_handler("Check your EPE input file. I think it is wrong")
  endif
#endif

  n_iterations=0
  ex_pgdata=.false.
!  call mrel_out2in()
!!$     print*, 'eped_driver: out2in'
        ! if no epein then use epeout instead
!  call system(trim(epedata_dir)//'out2in')

  call start_epe

!AG[
    if( comm_parallel() ) then 
!!$       call write_to_output_units("=> epe_send_init_to_slave")
          call epe_send_init_to_slave()
!!$       call write_to_output_units("<= epe_send_init_to_slave")
    end if
!AG]
  call main_epe   
    ex_pgdata=.true. 

  conv_unit=openget_iounit (trim(epedata_dir)//'conv',status='unknown')
  call returnclose_iounit(conv_unit,status='delete')

!  contains
  
!    subroutine mrel_out2in()
!    !  Purpose: makes "out2in" to cp files in make_reg_reference mode
!    !------------ Modules used ----------------------------------
!    implicit none
!   !------------ Declaration of local variables -----------------
!   integer             :: to_epe_script_file
!   character(len=100)  :: sys_command
!   !------------ Executable code --------------------------------
!
!   to_epe_script_file = openget_iounit (trim(epedata_dir)//'out2in')
!!!$        print*,trim(epedata_dir)//'out2in'
!   write(to_epe_script_file, '("#!/bin/csh ")' )
!   write(to_epe_script_file,'(a)')  &
!       "if (-e "//trim(epe_input_dir)//"/epeout && ! -e "//trim(epedata_dir)//"epein  ) then"
!   write(to_epe_script_file,'(a)') &
! "cp "//trim(epe_input_dir)//"/epeout "//trim(epedata_dir)//"epein "
!   write(to_epe_script_file,'(a)') &
!        "echo cp "//trim(epe_input_dir)//"/epeout "//trim(epedata_dir)//"epein "
!   write(to_epe_script_file, '("endif " )' )
!   call returnclose_iounit(to_epe_script_file)
!   call system('chmod 777 '//trim(epedata_dir)//'out2in')
!
!  end subroutine mrel_out2in


end subroutine epe_driver

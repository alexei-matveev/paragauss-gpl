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
program epe
!
!================================================================
! End of public interface of module
!================================================================
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

use type_module
use epecom_module , only :periodic_optimization,output_epe,epedata_dir, &
			  pg_cfield_calculated
use iounitadmin_module
use filename_module,only:data_dir
use main_epe_module, only: main_epe,start_epe,finish_epe

implicit none
 integer (kind=i4_kind):: conv_unit, namelength,i

	pg_cfield_calculated=.false.
	call getenv("TTFSDATADIR",epedata_dir)
	print*,'epedata_dir',epedata_dir
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

if(periodic_optimization) then 
  conv_unit=openget_iounit (trim(epedata_dir)//'conv',status='unknown')
  call returnclose_iounit (conv_unit)
endif ! periodic_optimization

call read_epe_input
call start_epe
!	do i=1,20
call main_epe
!	enddo ! k=1,20
call finish_epe()
 
 call returnclose_iounit(output_epe)
 call returnclose_iounit(conv_unit)
  conv_unit=openget_iounit (trim(epedata_dir)//'conv',status='unknown')
  call returnclose_iounit(conv_unit,status='delete')

end program epe

subroutine error_handler(message)
  ! write message to Standard error and stop
  implicit none
  character(len=*) :: message
  write(0,*) message
  stop 1
end subroutine error_handler

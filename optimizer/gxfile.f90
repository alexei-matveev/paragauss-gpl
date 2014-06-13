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
module gxfile
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
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

# include "def.h"
  use type_module, only:&
       & IK=>i4_kind, &
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------
!!$  type, public ::  gxfile_
!!$  end type gxfile_

  !------------ Declaration of constants and variables ------------
!!$  integer(kind=IK), parameter, public  :: gxfile_
!!$  real(kind=RK),    parameter, public  :: gxfile_
!!$  logical,               parameter, public  :: gxfile_
!!$  character,             parameter, public  :: gxfile_
!!$  integer(kind=IK),            public  :: gxfile_
!!$  real(kind=RK),               public  :: gxfile_
!!$  logical,                          public  :: gxfile_
!!$  character,                        public  :: gxfile_


  !------------ Interface statements ------------------------------
!!$  interface gxfile_
!!$  end interfacegxfile_
!!$  public gxfile_

  !------------ public functions and subroutines ------------------
  public gxfile_read

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------
!!$  type
!!$  end type

  !------------ Declaration of constants and variables ----
!!$  integer(IK), parameter :: max_atoms = 999
!!$  real(kind=RK),    parameter ::
!!$  logical,               parameter ::
!!$  character,             parameter ::
!!$  integer(kind=IK),           ::
!!$  real(kind=RK),              ::
!!$  logical,                         ::
!!$  character,                       ::



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine gxfile_read( io_gx, n_centers, geo_loop, charge, xyz, &
       & unique, equiv, zmat, numx, iepe, cart, calc_epeff_hessian)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use io
!    use opt_data_module, only: calc_epeff_hessian
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in)  :: io_gx
    real(RK),    intent(out) :: geo_loop
    integer(IK), intent(out) :: n_centers
    real(RK),    intent(inout) :: charge(:) ! (max_atoms)
    real(RK),    intent(out) :: xyz(:,:)  ! (3,max_atoms)
    integer(IK), intent(out) :: zmat(:,:) ! (3,max_atoms)
    integer(IK), intent(out) :: unique(:) ! (max_atoms)
    integer(IK), intent(out) :: equiv(:)  ! (max_atoms)
    integer(IK), intent(out) :: numx(:,:) ! (3,max_atoms)
    integer(IK), intent(out) :: iepe(:)   ! (max_atoms)
    logical,     intent(in)  :: cart,calc_epeff_hessian
    logical :: epeff_hessian
    optional :: unique, equiv, zmat, numx, iepe, cart,calc_epeff_hessian
    ! the first one that is missing makes
    ! ALL the subsequent to be ignored
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(IK)  :: new_numx(size(charge))
    integer(IK)  :: i,j, max_atoms
    integer(IK)  :: stat
    integer(IK)  :: num
    character    :: char
#if FPP_NO_F90_IO
    character (len=128) :: line
#endif
    !------------ Executable code --------------------------------

    DPRINT 'gxfile_read: enetered'
    if(present(calc_epeff_hessian)) then
    epeff_hessian=calc_epeff_hessian
    else
    epeff_hessian=.false.
    endif
    max_atoms = size(xyz,2)
    if (present(zmat)) then
       ASSERT(all(shape(xyz)==shape(zmat)))
    endif
    if (present(numx)) then
       ASSERT(all(shape(xyz)==shape(numx)))
    endif
    if (present(unique)) then
       ASSERT(size(unique)==max_atoms)
    endif
    if (present(equiv)) then
       ASSERT(size(equiv)==max_atoms)
    endif
    if( present(iepe) ) then
       ASSERT(size(iepe)==max_atoms)
    endif

    n_centers = 0
#if FPP_NO_F90_IO
    input_loop: do i=1,max_atoms
       DPRINT "gxfile_read: i=",i

       ! read one line:
       read(io_gx,'(A128)',iostat=stat) line
       !print *,'L=>'//line//'<, iostat=',stat
       ASSERT(stat==0)

       if( .not. present(iepe) )then
         read(line,*,iostat=stat)          &
                       charge(i)           &
                     , (xyz(j,i),j=1,3)    &
                     , unique(i), equiv(i) &
                     , (zmat(j,i),j=1,3)   &
                     , (numx(j,i),j=1,3)
       else
         ! one more field on the right:
         read(line,*,iostat=stat)          &
                       charge(i)           &
                     , (xyz(j,i),j=1,3)    &
                     , unique(i), equiv(i) &
                     , (zmat(j,i),j=1,3)   &
                     , (numx(j,i),j=1,3)   &
                     , iepe(i)
if(iepe(i).ne.0) then
ASSERT(charge(i)-aint(charge(i)).gt.0.0001)
endif
       endif
       !print *,'z=',charge(i),'iostat=',stat

       ! may fail here for legacy gxfiles
       ! dont ASSERT(stat==0)

       ! Instead, try reading only the first (negative!) charge:
       if( stat/=0 )then
         WARN('iostat/=0 while reading gxfile')
         read(line,*,iostat=stat)          &
                       charge(i)
         !print *,'g=',charge(i),'iostat=',stat
         ASSERT(stat==0)
         ASSERT(charge(i)<0)
       endif

       if (charge(i)<=0.5_RK) then
          geo_loop = charge(i)
          exit input_loop
       endif
       n_centers = n_centers + 1
    enddo input_loop
#else /* use nonadvancing IO of f90 */
    input_loop: do i=1,max_atoms
       DPRINT "gxfile_read: i=",i

       call readline(io_gx,charge(i),stat)
       DPRINT 'gxfile_read: charge<<',charge(i)
       ! may fail here for hand-made gxfiles
       !ASSERT(stat==IO_OK)

       if (charge(i)<=0.5_RK) then
          geo_loop = charge(i)
          exit input_loop
       endif
       ! but not here!
       ASSERT(stat==IO_OK)

       n_centers = n_centers + 1

       do j=1,3
          call readline(io_gx,xyz(j,i),stat)
          DPRINT 'gxfile_read: xyz(',j,')<<',xyz(j,i)
          ASSERT(stat==IO_OK)
       enddo

       if(.not.present(unique)) goto 888 ! cycle io line

       call readline(io_gx,unique(i),stat)
       DPRINT 'gxfile_read: ua<<',unique(i)
       ASSERT(stat==IO_OK)

       if(.not.present(equiv)) goto 888 ! cycle io line

       call readline(io_gx,equiv(i),stat)
       DPRINT 'gxfile_read: ea<<',equiv(i)
       ASSERT(stat==IO_OK)

       if(.not.present(zmat)) goto 888 ! cycle io line

       do j=1,3
          call readline(io_gx,zmat(j,i),stat)
          DPRINT 'gxfile_read: zmat(',j,')<<',zmat(j,i)
          ASSERT(stat==IO_OK)
       enddo

       if(.not.present(numx)) goto 888 ! cycle io line

       do j=1,3
          call readline(io_gx,numx(j,i),stat)
          DPRINT 'gxfile_read: numx(',j,')<<',numx(j,i)
          ASSERT(stat==IO_OK.or.j==3) !  stat may be IO_EOR here
       enddo

       if(.not.present(iepe)) goto 888 ! cycle io line

       call readline(io_gx,iepe(i),stat)
       DPRINT 'gxfile_read: iepe(',i,')<<',iepe(i)
       ! stat may be IO_EOR here
       if(iepe(i).ne.0) then

if(.not.epeff_hessian) then
ASSERT(charge(i)-aint(charge(i)).gt.0.0001)
endif
       endif

888    continue ! seek end of line:
       do while(stat==IO_OK)
          call readline(io_gx,char,stat)
       enddo
       ASSERT(stat==IO_EOR)

    enddo input_loop

    ! seek end of line:
    do while(stat==IO_OK)
       call readline(io_gx,char,stat)
    enddo
    ASSERT(stat==IO_EOR)
#endif

    DPRINT 'gxfile_read: loop=',geo_loop
#if 1
 if(.not.present(cart)) then
 num=0
 do j=1,3
 new_numx=0
  do i=2,n_centers
   if(numx(j,i).eq.0) cycle
   if(numx(j,i).ne.numx(j,i-1)) num=num+1
   new_numx(i)=num
  enddo
  numx(j,1:n_centers)=new_numx(1:n_centers)
 enddo
 do i=1,n_centers
  equiv(i)=i
 enddo
 end if
#endif
  end subroutine gxfile_read
  !*************************************************************


  !--------------- End of module ----------------------------------
end module gxfile

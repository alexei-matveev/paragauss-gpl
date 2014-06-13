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
module punchfile
  !---------------------------------------------------------------
  !
  !  Creates  GAMESS-style  punchfiles.   Files  in  this  format  are
  !  understood by CCP1 GUI, which uses VTK for visualizing molecules,
  !  densities, and orbitals.
  !
  !
  !  Module called by: orbital_plot_module.f90
  !
  !
  !  References: GAMESS-UK user guide and reference manual,
  !              http://www.cfs.dl.ac.uk/docs/html/part11/node3.html
  !
  !
  !  Author: AM
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

  use type_module, only: &
       IK=>i4_kind, &
       RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------
  interface pun_block_header
     module procedure pun_block_header_0
     module procedure pun_block_header_ind
     module procedure pun_block_header_att ! ( io, name, records, att, val )
  end interface

  !------------ public functions and subroutines ------------------
  public :: pun_block_header
  public :: pun_title
  public :: pun_coordinates
  public :: pun_connectivity
  public :: pun_grid_title
  public :: pun_grid_axes
  public :: pun_grid_mapping
  public :: pun_grid_data


  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine pun_block_header_0( io, name, records )
    !block=grid_title records= 1
    implicit none
    integer(IK)     , intent(in) :: io
    character(len=*), intent(in) :: name
    integer(IK)     , intent(in) :: records
    ! *** end of interface ***

    write(io,111) trim(name),records
111 FORMAT('block=',A,' records=',I6)
  end subroutine pun_block_header_0

  subroutine pun_block_header_ind( io, name, records, index )
    !block=grid_title records= 1 index =   1
    implicit none
    integer(IK)     , intent(in) :: io
    character(len=*), intent(in) :: name
    integer(IK)     , intent(in) :: records
    integer(IK)     , intent(in) :: index
    ! *** end of interface ***

    write(io,111) trim(name),records,index
111 FORMAT('block=',A,' records=',I6,' index=',I6)
  end subroutine pun_block_header_ind

  subroutine pun_block_header_att( io, name, records, att, val )
    !block=grid_data records=    2500 index=  1 elements =   1
    implicit none
    integer(IK)     , intent(in) :: io
    character(len=*), intent(in) :: name
    integer(IK)     , intent(in) :: records
    character(len=*), intent(in) :: att
    integer(IK)     , intent(in) :: val
    ! *** end of interface ***

    write(io,111) trim(name),records, att,val
111 FORMAT('block=',A,' records=',I6,1X,A,'=',I6)
  end subroutine pun_block_header_att

  subroutine pun_title( io, title )
    !block=grid_title records= 1 index =   1
    !density
    implicit none
    integer(IK)     , intent(in) :: io
    character(len=*), intent(in) :: title
    ! *** end of interface ***

    call pun_block_header(io,'title',1)
    write(io,111) title
111 FORMAT(A)
  end subroutine pun_title

  subroutine pun_coordinates( io, SYM, CRD )
    !block = coordinates records =    12 unit = au
    !  c             2.6541359     -0.0000157      0.0000000
    !  c             1.3269813      2.2986467      0.0000000
    ! ...
    !  h             4.7392596     -0.0000157      0.0000000
    !  h             2.3697321      4.1044689      0.0000000
    implicit none
    integer(IK)     , intent(in) :: io
    character*12    , intent(in) :: SYM(:)   ! (NA)
    real(RK)        , intent(in) :: CRD(:,:) ! (3,NA)
    ! *** end of interface ***

    integer(IK) :: NA,a

    NA = size(CRD,2)
    call pun_block_header(io,'coordinates',NA)
    do a=1,NA
       write(io,111) SYM(a), CRD(:3,a)
    enddo
111  FORMAT(2X,A8,3(3X,F12.7))
  end subroutine pun_coordinates

  subroutine pun_connectivity( io, BMAT )
    !block = connectivity records =   12
    !    1    2
    !    1    6
    !    1    7
    !    2    3
    implicit none
    integer(IK)     , intent(in) :: io
    integer(IK)     , intent(in) :: BMAT(:,:)
    ! *** end of interface ***

    integer(IK) :: I,J,N,NB

    N = min(size(BMAT,1),size(BMAT,2))
    NB = 0
    do J=1,N
       do I=1,J-1
          if( BMAT(I,J)/=0 ) NB = NB + 1
       enddo
    enddo
    call pun_block_header(io,'connectivity',NB)
!!$    do J=1,N
!!$       write(io,*) BMAT(J,:)
!!$    enddo
    do J=1,N
       do I=1,J-1
          if( BMAT(I,J)/=0 )then
             write(io,111) I,J
          endif
       enddo
    enddo
111 FORMAT(2I5)
  end subroutine pun_connectivity

  subroutine pun_grid_title( io, title )
    !block=grid_title records= 1 index =   1
    !density
    implicit none
    integer(IK)     , intent(in) :: io
    character(len=*), intent(in) :: title
    ! *** end of interface ***

    call pun_block_header(io,'grid_title',1)
    write(io,111) title
111 FORMAT(A)
  end subroutine pun_grid_title

  subroutine pun_grid_axes( io, N, D )
    !block=grid_axes records=   3 index=  1
    !  51   0.000000  16.156798 0  au xaxis
    !  51   0.000000  13.758351 0  au yaxis
    !  51   0.000000  15.626145 0  au zaxis
    implicit none
    integer(IK)     , intent(in) :: io
    integer(IK)     , intent(in) :: N(3)
    real(RK)        , intent(in) :: D(3)
    ! *** end of interface ***

    real(RK) , parameter :: O=0.0_rk

    call pun_block_header(io,'grid_axes',3)
    write(io,111) N(1), O,D(1), 0, 'au', 'xaxis'
    write(io,111) N(2), O,D(2), 0, 'au', 'yaxis'
    write(io,111) N(3), O,D(3), 0, 'au', 'zaxis'
111 FORMAT(1X,I3,2F11.6,I2,A3,1X,A5)
  end subroutine pun_grid_axes

  subroutine pun_grid_mapping( io, X0, X1, X2, X3 )
    !block=grid_mapping records=   3 index=  1
    ! -8.537405 -6.359682 -7.368608  7.619393 -6.359682 -7.368608
    ! -8.537405 -6.359682 -7.368608 -8.537405  7.398670 -7.368608
    ! -8.537405 -6.359682 -7.368608 -8.537405 -6.359682  8.257537
    implicit none
    integer(IK)     , intent(in) :: io
    real(RK)        , intent(in) :: X0(3),X1(3),X2(3),X3(3)
    ! *** end of interface ***

    call pun_block_header(io,'grid_mapping',3)
    write(io,111) X0, X1
    write(io,111) X0, X2
    write(io,111) X0, X3
111 FORMAT(6F10.6)
  end subroutine pun_grid_mapping

  subroutine pun_grid_data( io, data, records, index )
    !block=grid_mapping records=   3 index=  1
    ! -8.537405 -6.359682 -7.368608  7.619393 -6.359682 -7.368608
    ! -8.537405 -6.359682 -7.368608 -8.537405  7.398670 -7.368608
    ! -8.537405 -6.359682 -7.368608 -8.537405 -6.359682  8.257537
    implicit none
    integer(IK)     , intent(in) :: io
    real(RK)        , intent(in) :: data(:)
    integer(IK)     , intent(in) :: records, index
    optional :: records, index
    ! *** end of interface ***

    integer(IK) :: i

    if( present(records) )then
       ! first the header
       call pun_block_header(io,'grid_data',records, index)
    endif
    do i=1,size(data)
       write(io,111) data(i)
    enddo
111 FORMAT(F11.6)
  end subroutine pun_grid_data

  !--------------- End of module ----------------------------------
end module punchfile

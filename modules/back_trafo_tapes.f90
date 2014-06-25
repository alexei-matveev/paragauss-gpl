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
!===============================================================
! Public interface of module
!===============================================================
module back_trafo_tapes
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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
  use type_module, only: IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use filename_module, only: strlen=>filename_namelengthmax
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  interface put_matrix
     module procedure put_matrix_c
     module procedure put_matrix_d
  end interface

  interface get_matrix
     module procedure get_matrix_c
     module procedure get_matrix_d
  end interface

  !------------ public functions and subroutines ------------------

  public :: get_matrix!(name, irrep, matrix)
  public :: put_matrix!(name, irrep, matrix)

  public :: get_hmatrix!(name, irrep, matrix), legacy on-disk format?

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine get_hmatrix(fn, irr, MAT)
    use matrix_module, only: cmatrix, square
    use comm, only: comm_rank, comm_bcast
    implicit none
    character(len=*), intent(in) :: fn
    integer(IK), intent(in) :: irr
    type(cmatrix), intent(inout) :: MAT
    ! *** end of interface ***

    character(LEN=10)     :: buf
    character(len=strlen) :: fname

    if(irr.eq.-1)return !<<< NO_IRREP from back_trafo_module

    ASSERT(square(MAT))

    write(buf,'(I3)') irr
    buf = adjustl(buf)

    if (comm_rank() == 0) then
       fname = trim(fn) // '_real.dat' // trim(buf)
       call read_smatrix(fname, MAT%re)
       fname = trim(fn) // '_imag.dat' // trim(buf)
       call read_amatrix(fname, MAT%im)
    endif
    call comm_bcast(MAT%re)
    call comm_bcast(MAT%im)
  end subroutine get_hmatrix

  subroutine read_smatrix(fn,MM)
    ! readwriteblocked Symmetric matrix
    use readwriteblocked_module
    use filename_module, only: tmpfile ! environment
    implicit none
    character(len=*),intent(in)  :: fn
    real(RK)        ,intent(out) :: MM(:,:)
    ! *** end of interface ***

    character(len=strlen) :: fname
    type(readwriteblocked_tapehandle) :: tph
    integer(IK) :: m,MDIM,NDIM,TRIDIM


    MDIM = size(MM,1)
    NDIM = size(MM,2)
    ASSERT(MDIM==NDIM)
    ! ^ symmetric == square

    TRIDIM = MDIM*(MDIM+1)/2

    fname = trim(adjustl(tmpfile(fn)))
    DPRINT 'read_smatrix: <',trim(fname)
    call readwriteblocked_startread( fname, tph, TOTAL_LENGTH=TRIDIM )

    do m=1,MDIM
       call readwriteblocked_read( MM(1:m,m), tph )
       MM(m,1:m) = MM(1:m,m)
    enddo
    call readwriteblocked_stopread( tph, STATUS='keep' )
  end subroutine read_smatrix

  subroutine read_amatrix(fn,MM)
    ! readwriteblocked Asymmetric matrix
    use readwriteblocked_module
    use filename_module, only: tmpfile ! environment
    implicit none
    character(len=*),intent(in)  :: fn
    real(RK)        ,intent(out) :: MM(:,:)
    ! *** end of interface ***

    character(len=strlen) :: fname
    type(readwriteblocked_tapehandle) :: tph
    integer(IK) :: m,MDIM,NDIM,TRIDIM


    MDIM = size(MM,1)
    NDIM = size(MM,2)
    ASSERT(MDIM==NDIM)
    ! ^ symmetric == square

    TRIDIM = MDIM*(MDIM+1)/2

    fname = trim(adjustl(tmpfile(fn)))
    DPRINT 'read_amatrix: <',trim(fname)
    call readwriteblocked_startread( fname, tph, TOTAL_LENGTH=TRIDIM )

    do m=1,MDIM
       call readwriteblocked_read( MM(1:m,m), tph )
       MM(m,1:m) = - MM(1:m,m)
    enddo
    call readwriteblocked_stopread( tph, STATUS='keep' )
  end subroutine read_amatrix

  subroutine put_matrix_c(fn, irr, MAT)
    !
    ! See get_matrix_c(). Runs in parallel context.
    !
    use matrix_module, only: cmatrix
    use quadrupel_fname, only: qfilename
    use io, only: file_write
    implicit none
    character(len=*), intent(in)   , optional :: fn
    integer(IK)     , intent(in)   , optional :: irr
    type(cmatrix)   , intent(inout), optional :: MAT
    ! *** end of interface ***

    real(RK), allocatable :: buf(:, :, :) ! FIXME: temp storage

    ASSERT(irr/=-1)
    !<<< NO_IRREP from back_trafo_module

    allocate(buf(size(MAT%re,1), size(MAT%re,2), 2))
    buf(:, :, 1) = MAT%re
    buf(:, :, 2) = MAT%im

    call file_write(qfilename(trim(fn), irr, EXT='CMX'), buf)
  end subroutine put_matrix_c

  subroutine get_matrix_c(fn, irr, MAT)
    !
    ! See put_matrix_c().
    !
    use matrix_module, only: cmatrix
    use quadrupel_fname, only: qfilename
    use io, only: file_read
    implicit none
    character(len=*), intent(in)    :: fn
    integer(IK)     , intent(in)    :: irr
    type(cmatrix)   , intent(inout) :: MAT ! needs to be allocated
    ! *** end of interface ***

    real(RK) :: buf(size(MAT%re,1), size(MAT%re,2), 2) ! FIXME: temp storage

    if(irr.eq.-1)return !<<< NO_IRREP from back_trafo_module

    call file_read(qfilename(trim(fn), irr, 'CMX'), buf)
    MAT%re = buf(:, :, 1)
    MAT%im = buf(:, :, 2)
  end subroutine get_matrix_c

  subroutine put_matrix_d(fn, irr, MAT)
    !
    ! Saves a matrix for  future retrival by get_matrix_d().  Runs
    ! in parallel context.
    !
    use matrix_module, only: rdmatrix
    use io, only: file_write
    use quadrupel_fname, only: qfilename
    implicit none
    character(len=*), intent(in) :: fn
    integer(IK), intent(in) :: irr
    type(rdmatrix), intent(in) :: MAT
    ! *** end of interface ***

    ASSERT(irr/=-1)
    !<<< NO_IRREP from back_trafo_module

    call file_write(qfilename(trim(fn), irr, 'RDX'), MAT%d)
  end subroutine put_matrix_d

  subroutine get_matrix_d(fn,irr,MAT)
    !
    ! Counterpart to put_matrix_d().
    !
    use matrix_module, only: rdmatrix
    use quadrupel_fname, only: qfilename
    use io, only: file_read
    implicit none
    character(len=*), intent(in)    :: fn
    integer(IK)     , intent(in)    :: irr
    type(rdmatrix)  , intent(inout) :: MAT
    ! *** end of interface ***

    if (irr.eq.-1) return !<<< NO_IRREP from back_trafo_module

    call file_read(qfilename(trim(fn), irr, 'RDX'), MAT%d)
  end subroutine get_matrix_d

  !--------------- End of module ----------------------------------
end module back_trafo_tapes

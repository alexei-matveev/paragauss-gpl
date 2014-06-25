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
!================================================================
! Public interface of module
!================================================================
module prepare_integralfiles_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: opens and closes the integral files which were produced in
  !           the integral part of the (up to now OLD) lcgto.
  !           It is necessary to perform this in extra routines,since
  !           the buffered I/O requires a 'startread' and a 'stopread'.
  !           (For details see: readwrite_blocked_module). It is intended
  !           to keep the integral files open during the loop over IRREPs
  !           in the routine 'build_hamiltonian', thus keeping the 
  !           tapehandle`s io_unit fixed.
  !           If NUM_EXCH is set, the routines should leave out the
  !           exchange files.
  !
  !  Contents:
  !      Routines:        prepare_integral_open -> determine which file to open,
  !                                     open it and call the routine
  !                                     'startread' which yields a 
  !                                     tapehandle. The tapehandle will be
  !                                     given to the routine build_hamiltonian
  !                                     which in turn hands it to
  !                                     'ham_calc_master'.
  !                                     If NUM_EXCH is set, an empty tapehandle
  !                                     th_xc is returned to build_ham.
  !                                     This is because derived types cannot
  !                                     be optional arguments.
  !
  !                       prepare_integral_close -> performs the 'stopread' 
  !                                      routine,
  !                                      returns the io-units and closes
  !                                      the files. Empty tapehandles should
  !                                      be handled correctly.
  !
  !  Module called by: build_hamiltonian
  !
  !
  !  References: yet none
  ! 
  !
  !  Author: Folke Noertemann
  !  Date: 12/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Adaption to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use iounitadmin_module              ! provides I/O units
  use readwriteblocked_module,only : readwriteblocked_tapehandle, &
       readwriteblocked_startread,&
       readwriteblocked_stopread
  use options_module, only: options_spin_orbit
  implicit none
  private
  save


  !------------ public functions and subroutines ------------------
  public prepare_integral_open,prepare_integral_close


!================================================================
! End of public interface of module
!================================================================

  public readwriteblocked_tapehandle


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine prepare_integral_open(basename, th_real, th_imag)
    !
    ! Purpose: open files, perform 'startread'
    !
    ! Modules:
    use dimensions,   only: IrrDim=>IrrBasDimSpor !<<< contracted dims
    use bounds_module, only: bounds_ch, bounds_xc
    use comm,  only: comm_rank
    use filename_module, only: tmpfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)                             :: basename ! 'coul' or 'exch'
    type(readwriteblocked_tapehandle), intent(out)           :: th_real
    type(readwriteblocked_tapehandle), intent(out), optional :: th_imag ! in case of spin orbit
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: rank
    integer(i4_kind) :: n_ij, n_k, n_ijk
    
    if ( .not. options_spin_orbit ) then
       call readwriteblocked_startread(trim(tmpfile(basename // '.dat')), th_real)
    else
       !
       ! SPIN ORBIT
       !
       ASSERT(present(th_imag))
       ASSERT(allocated(IrrDim))

       rank = comm_rank()

       !
       ! Now find out the bounds_xx`s: --------------
       ! they are needed to decide which files are to be opened
       !
       ! We refer to both bounds_ch and bounds_xc whether
       ! num_exch=.true. or not, because we want to make use
       ! of the emtpy bounds_xc array in the case num_exch=.true.
       !
       ASSERT(allocated(bounds_ch))
       ASSERT(allocated(bounds_xc))
       select case ( basename )
       case ( 'coul' )
           n_k = bounds_ch(rank + 1)
       case ( 'exch' )
           n_k = bounds_xc(rank + 1)
           WARN('verify total_length!')
       case default
           n_k = -1
           ABORT('no such case')
       end select

       if ( n_k <= 0 ) then
           WARN('opening empty file?')
       endif

       n_ij = sum( IrrDim(:) * (IrrDim(:) + 1) / 2 )
       n_ijk = n_ij * n_k

       ! FIXME: total_length for exchange?
       call readwriteblocked_startread( trim(tmpfile(basename // '_real.dat')) &
                                      , th_real, total_length=n_ijk)
       call readwriteblocked_startread( trim(tmpfile(basename // '_imag.dat')) &
                                      , th_imag, total_length=n_ijk)
    endif
  end subroutine prepare_integral_open
  !*************************************************************


  !*************************************************************
  subroutine prepare_integral_close(th_ch,th_xc,th_ch_imag,th_xc_imag)
    !  Purpose: performs a 'stopread' on the files contained
    !           in the specified tapehandle,
    !           returns the corresponding io_unit and
    !           closes the files
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle),intent(inout),optional :: th_ch
    type(readwriteblocked_tapehandle),intent(inout),optional :: th_xc
    type(readwriteblocked_tapehandle),intent(inout),optional :: th_ch_imag
    type(readwriteblocked_tapehandle),intent(inout),optional :: th_xc_imag
    !** End of interface *****************************************      
    if (present(th_ch)) then
       call readwriteblocked_stopread(th_ch)
       if (options_spin_orbit) then
          call readwriteblocked_stopread(th_ch_imag)
       endif
     endif
    if (present(th_xc)) then
       call readwriteblocked_stopread(th_xc)
       if (options_spin_orbit) then
          call readwriteblocked_stopread(th_xc_imag)
       endif
    endif
  end subroutine prepare_integral_close
  !*************************************************************


!--------------- End of module ----------------------------------
end module prepare_integralfiles_module

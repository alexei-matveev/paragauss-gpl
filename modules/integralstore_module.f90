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
module  integralstore_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Hold precalculated integrals in case they are not
  !           written to file
  !
  !  Module called by: int_send_2cob3c_module,
  !
  !
  !  Author: TB
  !  Date: 6/97
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   9/97
  ! Description: New arrays for relativistic matrixelements have been
  !              added
  !
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  use options_module, only: options_spin_orbit
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of constants and variables ---------------
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_kin
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_nuc
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_efield
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_kin_rel
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_nuc_rel
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_pvsp
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_pvxp
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_sigp
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_ol
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_2cob_ol_rel
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_3c_xc
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_3c_co
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_3c_poten
  real(r8_kind), allocatable, target, dimension(:), public :: integralstore_3c_field

  ! datasructures for relativistiv gradients
  !------------ public functions and subroutines ---------------------
  public :: integralstore_allocate
  public :: integralstore_deallocate
  public :: integralstore_allocate_efield
  public :: integralstore_deallocate_efield
  public :: integralstore_allocate_pcm
  public :: integralstore_deallocate_pcm
  public :: integralstore_kin_and_nuc_to_mem

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  function integralstore_allocate ( dim_2c, dim_2c_rel, dim_3c_co, &
       dim_3c_xc, need_2cob_kin, need_2cob_nuc, need_2cob_kin_rel, &
       need_2cob_nuc_rel, need_2cob_pvsp,&
       need_2cob_pvxp,&
       need_2cob_pvec,&
       need_2cob_ol, need_2cob_ol_rel, &
       need_3c_xc, need_3c_co ) result(error)
    !  Purpose: allocates the arrays specified by the logical
    !           variables to the given dimensions.
    !  returns .true. if allocation failed
    implicit none
    !------------ Declaration of formal parameters -------------
    logical :: error
    ! dimensions
    integer(kind=i4_kind), intent(in), optional :: dim_2c                      &
                                                 , dim_2c_rel                  &
                                                 , dim_3c_co                   &
                                                 , dim_3c_xc
    ! flags, denoting required arrays
    logical, intent(in), optional :: need_2cob_kin                             &
                                   , need_2cob_nuc                             &
                                   , need_2cob_kin_rel                         &
                                   , need_2cob_nuc_rel                         &
                                   , need_2cob_pvsp                            &
                                   , need_2cob_pvxp                            &
                                   , need_2cob_pvec                            &
                                   , need_2cob_ol                              &
                                   , need_2cob_ol_rel                          &
                                   , need_3c_xc                                &
                                   , need_3c_co
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: status,n_complex
    !------------ Executable code ------------------------------

    n_complex = 1
    if(options_spin_orbit) then
       n_complex = 2 ! complex variables require double space
    endif
    error = .true.
    if (present(need_2cob_kin)) then
      if ( need_2cob_kin ) then
        ASSERT(present(dim_2c))
        allocate ( integralstore_2cob_kin(n_complex*dim_2c), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_ol)) then
      if ( need_2cob_ol ) then
        ASSERT(present(dim_2c))
        allocate ( integralstore_2cob_ol(n_complex*dim_2c), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_nuc)) then
      if ( need_2cob_nuc ) then
        ASSERT(present(dim_2c))
        allocate ( integralstore_2cob_nuc(n_complex*dim_2c), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_kin_rel)) then
      if ( need_2cob_kin_rel ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_kin_rel(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_ol_rel)) then
      if ( need_2cob_ol_rel ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_ol_rel(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_nuc_rel)) then
      if ( need_2cob_nuc_rel ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_nuc_rel(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_pvsp)) then
      if ( need_2cob_pvsp ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_pvsp(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_pvxp)) then
      if ( need_2cob_pvxp ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_pvxp(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_2cob_pvec)) then
      if ( need_2cob_pvec ) then
        ASSERT(present(dim_2c_rel))
        allocate ( integralstore_2cob_sigp(n_complex*dim_2c_rel), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_3c_xc)) then
      if ( need_3c_xc ) then
        ASSERT(present(dim_3c_xc))
        allocate ( integralstore_3c_xc(n_complex*dim_3c_xc), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    if (present(need_3c_co)) then
      if ( need_3c_co ) then
        ASSERT(present(dim_3c_co))
        allocate ( integralstore_3c_co(n_complex*dim_3c_co), stat=status )
        if ( status .ne. 0 ) return
      endif
    endif
    error = .false.
  end function integralstore_allocate
  !*************************************************************

  !*************************************************************
  subroutine integralstore_deallocate( deallocate_kin                          &
                                     , deallocate_kin_rel                      &
                                     , deallocate_ol                           &
                                     , deallocate_ol_rel                       &
                                     , deallocate_nuc                          &
                                     , deallocate_nuc_rel                      &
                                     , deallocate_pvsp                         &
                                     , deallocate_pvxp                         &
                                     , deallocate_sigp                         &
                                     , deallocate_efield                       &
                                     , deallocate_3c_xc                        &
                                     , deallocate_3c_co                        )
    !  Purpose: deallocates all allocated integral arrays
    !** End of interface ***************************************
    logical, intent(in), optional :: deallocate_kin                            &
                                   , deallocate_kin_rel                        &
                                   , deallocate_ol                             &
                                   , deallocate_ol_rel                         &
                                   , deallocate_nuc                            &
                                   , deallocate_nuc_rel                        &
                                   , deallocate_pvsp                           &
                                   , deallocate_pvxp                           &
                                   , deallocate_sigp                           &
                                   , deallocate_efield                         &
                                   , deallocate_3c_xc                          &
                                   , deallocate_3c_co
    !------------ Declaration of local variables ---------------
    logical                       :: dealloc_all                               &
                                   , dealloc_kin                               &
                                   , dealloc_kin_rel                           &
                                   , dealloc_ol                                &
                                   , dealloc_ol_rel                            &
                                   , dealloc_nuc                               &
                                   , dealloc_nuc_rel                           &
                                   , dealloc_pvsp                              &
                                   , dealloc_pvxp                              &
                                   , dealloc_sigp                              &
                                   , dealloc_efield                            &
                                   , dealloc_3c_xc                             &
                                   , dealloc_3c_co
    integer(kind=i4_kind)         :: status
    !------------ Executable code ------------------------------
    !
    dealloc_all     = .FALSE.
    dealloc_kin     = .FALSE.
    dealloc_kin_rel = .FALSE.
    dealloc_ol      = .FALSE.
    dealloc_ol_rel  = .FALSE.
    dealloc_nuc     = .FALSE.
    dealloc_nuc_rel = .FALSE.
    dealloc_pvsp    = .FALSE.
    dealloc_pvxp    = .FALSE.
    dealloc_sigp    = .FALSE.
    dealloc_efield  = .FALSE.
    dealloc_3c_xc   = .FALSE.
    dealloc_3c_co   = .FALSE.
    !
    if (.not. present(deallocate_kin    ) .and.                                &
        .not. present(deallocate_kin_rel) .and.                                &
        .not. present(deallocate_ol     ) .and.                                &
        .not. present(deallocate_ol_rel ) .and.                                &
        .not. present(deallocate_nuc    ) .and.                                &
        .not. present(deallocate_nuc_rel) .and.                                &
        .not. present(deallocate_pvsp   ) .and.                                &
        .not. present(deallocate_pvxp   ) .and.                                &
        .not. present(deallocate_sigp   ) .and.                                &
        .not. present(deallocate_efield ) .and.                                &
        .not. present(deallocate_3c_xc  ) .and.                                &
        .not. present(deallocate_3c_co  )       ) then
       !
       ! No argument present. Deallocate all allocated arrays (default)
       !
       dealloc_all = .TRUE.
       !
    else
       !
       ! Some arguments present. Deallocate all only some specific arrays
       !
       dealloc_all = .FALSE.
       !
       if (present(deallocate_kin    )) dealloc_kin     = deallocate_kin
       if (present(deallocate_kin_rel)) dealloc_kin_rel = deallocate_kin_rel
       if (present(deallocate_ol     )) dealloc_ol      = deallocate_ol
       if (present(deallocate_ol_rel )) dealloc_ol_rel  = deallocate_ol_rel
       if (present(deallocate_nuc    )) dealloc_nuc     = deallocate_nuc
       if (present(deallocate_nuc_rel)) dealloc_nuc_rel = deallocate_nuc_rel
       if (present(deallocate_pvsp   )) dealloc_pvsp    = deallocate_pvsp
       if (present(deallocate_pvxp   )) dealloc_pvxp    = deallocate_pvxp
       if (present(deallocate_sigp   )) dealloc_sigp    = deallocate_sigp
       if (present(deallocate_efield )) dealloc_efield  = deallocate_efield
       if (present(deallocate_3c_xc  )) dealloc_3c_xc   = deallocate_3c_xc
       if (present(deallocate_3c_co  )) dealloc_3c_co   = deallocate_3c_co
       !
    endif
    !
    if (dealloc_all .or. dealloc_kin    ) then
       if ( allocated(integralstore_2cob_kin) ) then
          deallocate( integralstore_2cob_kin, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_kin failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_kin_rel) then
       if ( allocated(integralstore_2cob_kin_rel) ) then
          deallocate( integralstore_2cob_kin_rel, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_kin_rel failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_ol     ) then
       if ( allocated(integralstore_2cob_ol) ) then
          deallocate( integralstore_2cob_ol, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_ol failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_ol_rel ) then
       if ( allocated(integralstore_2cob_ol_rel) ) then
          deallocate( integralstore_2cob_ol_rel, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_ol_rel failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_nuc    ) then
       if ( allocated(integralstore_2cob_nuc) ) then
          deallocate( integralstore_2cob_nuc, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_nuc failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_nuc_rel) then
       if ( allocated(integralstore_2cob_nuc_rel) ) then
          deallocate( integralstore_2cob_nuc_rel, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_nuc_rel failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_pvsp   ) then
       if ( allocated(integralstore_2cob_pvsp) ) then
          deallocate( integralstore_2cob_pvsp, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_pvsp failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_pvxp   ) then
       if ( allocated(integralstore_2cob_pvxp) ) then
          deallocate( integralstore_2cob_pvxp, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_pvxp failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_sigp   ) then
       if ( allocated(integralstore_2cob_sigp) ) then
          deallocate( integralstore_2cob_sigp, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_sigp failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_efield ) then
       if ( allocated(integralstore_2cob_efield) ) then
          deallocate( integralstore_2cob_efield, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_2cob_efield failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_3c_xc  ) then
       if( allocated(integralstore_3c_xc) ) then
          deallocate( integralstore_3c_xc, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_3c_xc failed" )
       endif
    endif
    !
    if (dealloc_all .or. dealloc_3c_co  ) then
       if ( allocated(integralstore_3c_co) ) then
          deallocate( integralstore_3c_co, stat=status )
          if ( status .ne. 0) call error_handler( &
               "integralstore_deallocate: of integralstore_ failed" )
       endif
    endif
    !
  end subroutine integralstore_deallocate
  !*************************************************************

  !*************************************************************
  subroutine integralstore_allocate_efield(dim_2c)
    !  Purpose: allocates integralstore_2cob_efield
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in) :: dim_2c
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: status
    !------------ Executable code ------------------------------
    allocate( integralstore_2cob_efield(dim_2c), stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integralstore_allocate_efield: allocate failed" )
  end subroutine integralstore_allocate_efield
  !*************************************************************

  !*************************************************************
  subroutine integralstore_deallocate_efield()
    !  Purpose: allocates integralstore_2cob_efield
    !** End of interface ***************************************
    implicit none
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: status
    !------------ Executable code ------------------------------
    deallocate( integralstore_2cob_efield, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integralstore_deallocate_efield: deallocate failed" )
  end subroutine integralstore_deallocate_efield
  !*************************************************************

  !*************************************************************
  function integralstore_allocate_pcm(dim_2c,dim_3c_poten,dim_3c_field,&
        need_3c_poten,need_3c_field) &
     result(error)
    !  Purpose: allocates integralstore_2cob_poten
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in) :: dim_3c_poten,dim_3c_field,dim_2c
    logical              , intent(in) :: need_3c_poten,need_3c_field
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    logical :: error
    integer(kind=i4_kind) :: status
    !------------ Executable code ------------------------------
    error=.true.
    if(need_3c_poten .or. need_3c_field) then                        !!!!!!!!!!!!!
       if(.not. allocated(integralstore_2cob_ol)) then
          allocate ( integralstore_2cob_ol(dim_2c), stat=status )       !!!!!!!!!!!!!
          if ( status .ne. 0 ) return                                   !!!!!!!!!!!!!
       end if
    end if                                                           !!!!!!!!!!!!!
    if(need_3c_poten) then
      allocate( integralstore_3c_poten(dim_3c_poten), stat=status )
      if ( status .ne. 0 ) return
    endif
    if(need_3c_field) then
      allocate( integralstore_3c_field(dim_3c_field), stat=status )
      if ( status .ne. 0 ) return
    endif
    error=.false.
  end function integralstore_allocate_pcm
  !*************************************************************

  !*************************************************************               !!TODO merge with integralstore_deallocate
  subroutine integralstore_deallocate_pcm()                                    !
    !  Purpose: allocates integralstore_3c_poten                               !
    !** End of interface ***************************************               !
    implicit none                                                              !
    !------------ Declaration of local variables ---------------               !
    integer(kind=i4_kind) :: status                                            !
    !------------ Executable code ------------------------------               !
    if ( allocated(integralstore_2cob_ol) ) then                               !
       deallocate( integralstore_2cob_ol, stat=status )                        !
       if ( status .ne. 0) call error_handler( &                               !
            "integralstore_deallocate: of integralstore_2cob_ol failed(1)" )   !
    endif                                                                      !
    if(allocated( integralstore_3c_poten)) then                                !
        deallocate( integralstore_3c_poten, stat=status )                      !
        if ( status .ne. 0 ) call error_handler( &                             !
         "integralstore_deallocate_pcm: deallocate _3c_poten failed" )         !
    endif                                                                      !
    if(allocated( integralstore_3c_field)) then                                !
        deallocate( integralstore_3c_field, stat=status )                      !
        if ( status .ne. 0 ) call error_handler( &                             !
         "integralstore_deallocate_pcm: deallocate _3c_field failed" )         !
    endif                                                                      !
  end subroutine integralstore_deallocate_pcm                                  !
  !*************************************************************               !!TODO merge with integralstore_deallocate

  !*************************************************************
  subroutine integralstore_kin_and_nuc_to_mem()
    ! Purpose: Transfer integrals for kin and nuc stored in files
    !          to memory
    use options_module, only      : options_integrals_on_file                  &
                                  , options_spin_orbit
    use dimensions, only          : IrrBasDim                                  &
                                  , IrrBasDimSpor
    use iounitadmin_module,   only: get_iounit                                 &
                                  , return_iounit
    use filename_module,      only: tmpfile
    !
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind)                :: dim_2c
    !
    integer(kind=i4_kind)                :: io_u
    integer(kind=i4_kind)                :: io_u_real, io_u_imag
    !
    integer(kind=i4_kind)                :: i_meta
    !------------ Executable code ------------------------------
    !
    ! Skip if integrals already in memory
    if ( .not. options_integrals_on_file()) return
    !
    ! Set dimensions
    if (options_spin_orbit) then
      !
      !  SPIN ORBIT
      !
      dim_2c = sum(IrrBasDimspor * (IrrBasDimspor + 1) / 2)
      !
    else ! options_spin_orbit
      !
      ! STANDARD SCF (NO SPIN ORBIT)
      !
      dim_2c =  sum(IrrBasDim * (IrrBasDim + 1) / 2)
      !
    endif
    !
    ! allocation call
    if (integralstore_allocate( dim_2c        = dim_2c                         &
                              , need_2cob_kin = .TRUE.                         &
                              , need_2cob_nuc = .TRUE.))                       &
    call error_handler("integralstore_file_to_mem: allocation of integral storage failed")
    !
    !
    ! open file(s)
    if (options_spin_orbit) then
      !
      !  SPIN ORBIT
      !
      io_u_real = get_iounit()
      open(io_u_real,form='unformatted',status='unknown',                      &
           file=trim(tmpfile('ham_kin_nuc_real.dat')))
      io_u_imag = get_iounit()
      open(io_u_imag,form='unformatted',status='unknown',                      &
           file=trim(tmpfile('ham_kin_nuc_imag.dat')))
      !
      ! transfer data
      do i_meta = 1, dim_2c * 2, 2
        !
        read(io_u_real) integralstore_2cob_kin(i_meta)                         &
                      , integralstore_2cob_nuc(i_meta)
        !
        read(io_u_imag) integralstore_2cob_kin(i_meta + 1)                     &
                      , integralstore_2cob_nuc(i_meta + 1)
        !
      enddo
      !
      call return_iounit(io_u_real)
      close (io_u_real)
      call return_iounit(io_u_imag)
      close (io_u_imag)
      !
    else ! options_spin_orbit
      !
      ! STANDARD SCF (NO SPIN ORBIT)
      !
      ! open file
      io_u      = get_iounit()
      open(io_u,form='unformatted',status='unknown',                           &
           file=trim(tmpfile('ham_kin_nuc.dat')))
      !
      ! transfer data
      do i_meta = 1, dim_2c
        !
        read(io_u) integralstore_2cob_kin(i_meta)                              &
                 , integralstore_2cob_nuc(i_meta)
        !
      enddo
      !
      call return_iounit(io_u)
      close (io_u)
      !
    endif
    !
  end subroutine integralstore_kin_and_nuc_to_mem
  !*************************************************************


  !--------------- End of module -------------------------------
end module integralstore_module

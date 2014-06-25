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
module quadrupel_fname
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
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  public :: init!()
  public :: done!()
  public :: qfilename!(name, int, ext) -> "name_int.ext"
  public :: qfilename3!(name, int, int, int, ext) -> "name_int_int_int.ext"
  public :: bipel_index
  public :: max_bipel_index

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of constants and variables ----

  logical,private             :: initialized = .false.
  integer(IK),private         :: refcount = 0

  integer(IK),private         ::  MyIndex     = -1

  integer(IK),allocatable,private :: ua_fileindex(:)
  integer(IK),private             :: max_bipel_index_ = -1

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init()
    use comm, only: comm_rank
    use filename_module, only: filesystem_is_parallel
    use unique_atom_module
    implicit none
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: a,ix

    refcount = refcount + 1

    if(initialized)then
       WARN("qn/init: already called")
       return
    endif

    ! set global var:
    MyIndex = comm_rank() + 1

    if ( filesystem_is_parallel ) then
       !
       ! FIXME: workers may share files ...
       !
       ABORT("needs update")
    endif

    allocate(ua_fileindex(n_unique_atoms),STAT=memstat)
    ASSERT(memstat==0)

    ix = 0
    do a = 1,n_unique_atoms
       ix = ix +1
       ua_fileindex(a) = ix
       ix = ix + unique_atoms(a)%lmax_ob
    enddo
    max_bipel_index_ = ix

    initialized = .true.
  end subroutine init

  subroutine done()
    implicit none
    ! *** end of interface ***

    integer :: memstat

    refcount = refcount - 1

    if(.not.initialized)then
       WARN("qn/done: not initialized")
       return
    endif

    if(refcount > 0)then
       WARN("qn/done: refcount nonzero")
       return
    endif

    deallocate(ua_fileindex,STAT=memstat)
    ASSERT(memstat==0)

    max_bipel_index_ = -1

    initialized = .false.
  end subroutine done

  function max_bipel_index() result(ix)
    implicit none
    integer(IK) :: ix !<<< ressult
    ! *** end of interface ***

    ASSERT(initialized)

    ix = max_bipel_index_
  end function max_bipel_index

  function bipel_index(i_ua,i_l) result(ix)
    implicit none
    integer(IK), intent(in) :: i_ua,i_l
    integer(IK)             :: ix !<<<result  
    ! *** end of interface ***

    ASSERT(initialized)

    ix = ua_fileindex(i_ua) + i_l
  end function bipel_index

  function qfilename3(name, irr, ix1, ix2, ext) result(path)
    use filename_module, only: strlen=>filename_namelengthmax, tmpfile
    implicit none
    character(len=*), intent(in) :: name
    integer(IK),      intent(in) :: irr, ix1, ix2
    character(len=*), intent(in) :: ext
    character(len=strlen)        :: path ! result
    ! *** end of interface ***

    character(len=9) :: charnbr

    path = repeat(" ", strlen)
    
    ! prefix:
    path = trim(name)

    ! irrep number:
    ASSERT(irr<10**9)
    write (charnbr,'(i9)') irr
    charnbr = adjustl(charnbr)

    path = trim(path) // "_" // trim(charnbr)

    ! bipel 1:
    ASSERT(ix1<10**9)
    write (charnbr,'(i9)') ix1
    charnbr = adjustl(charnbr)

    path = trim(path) // "_" // trim(charnbr)

    ! bipel 2:
    ASSERT(ix2<10**9)
    write (charnbr, '(i9)') ix2
    charnbr = adjustl(charnbr)

    path = trim(path) // "_" // trim(charnbr)

    ! extension:
    path = trim(path) // "." // ext

    ! full path:
    path = tmpfile(path)
  end function qfilename3

  function qfilename(name, irr, ext, FullPath) result(path)
    use filename_module, only: strlen=>filename_namelengthmax, tmpfile
    implicit none
    character(len=*),    intent(in) :: name
    integer(IK),         intent(in) :: irr
    character(len=*),    intent(in) :: ext
    logical, optional,   intent(in) :: FullPath ! FIXME: must die!
    character(len=strlen)           :: path ! result
    ! *** end of interface ***

    character(len=9) :: charnbr
    logical :: absolute

    path = repeat(" ", strlen)

    ASSERT(irr<10**9)
    write (charnbr,'(i9)') irr
    charnbr = adjustl(charnbr)

    path = trim(name) // "_" // trim(charnbr) // "." // trim(ext)

    absolute = .true.
    if ( present(FullPath) ) then
        absolute = FullPath
    endif

    if ( absolute ) path = tmpfile(path)
  end function qfilename

  !--------------- End of module ----------------------------------
end module quadrupel_fname

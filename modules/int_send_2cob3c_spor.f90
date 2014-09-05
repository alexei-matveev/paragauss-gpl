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
module int_send_2cob3c_spor
  !-------------------------------------------------------------------
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
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind! type specification parameters
  !use error_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  interface MaxFitIndex
     module procedure MaxFitIndex_Which
     module procedure MaxFitIndex_ChFit
  end interface

  interface free
     module procedure free_SplitDims
  end interface

  !------------ public functions and subroutines ---------------------

  public :: &
       FitBorders, &
       SetSplitting, &
       MaxFitIndex, &
       MyFitIndex, &
       SyncBorders

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  type,private :: SplitDims
     integer(IK)          :: first,last
     integer(IK),pointer  :: dims(:)
     !...(first:last)
  end type SplitDims

  !------------ Declaration of constants and variables ---------------

  logical,parameter :: YES=.true.,NO=.false.
  integer(IK)       :: i_

  integer(IK),parameter,public  :: &
       & IX_CHFIT     = 1, &
       & IX_XCFIT     = 2,&
       & IX_CHFIT_S   = 3, &
       & IX_CHFIT_R2  = 4
  integer(IK),parameter,private :: &
       & NN           = 4 ! max of the above

  type(SplitDims),private       :: List(NN)
  logical,private               :: Set(NN) = (/(NO,i_=1,NN)/)


  logical,private :: initialized=NO  
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine init()
    implicit none
    ! *** end of interface ***

    initialized = YES
  end subroutine init

  subroutine done()
    implicit none
    ! *** end of interface ***

    integer(IK) :: ix

    do ix=1,NN
       if(Set(ix))call free(List(ix))
    enddo

    initialized = NO
  end subroutine done

  function MyFitIndex(typ,ff) result(yes)
    ! does the function ff
    ! of type typ (r2,s,..) is assigned to me?
    use comm_module
    implicit none
    integer(IK),intent(in) :: typ,ff
    logical                :: yes !<<< result
    ! *** end of interface ***

    integer(IK)     :: id

    id = comm_myindex()

    yes = (id .eq. WhoseFitIndex(typ,ff))
  end function MyFitIndex

  function WhoseFitIndex(typ,ff) result(proc)
    ! whom (as a processor) the function ff
    ! of type typ (r2,s,..) is assigned to 
    implicit none
    integer(IK),intent(in) :: typ,ff
    integer(IK)            :: proc !<<< result
    ! *** end of interface ***

    type(SplitDims) :: sd
    integer(IK)     :: tot, p

    ASSERT(Set(typ))
    ASSERT(ff>=0)

    sd = List(typ)

    proc = -1
    tot = 0
    do p = sd%first, sd%last
       tot = tot + sd%dims(p)
       if( tot >= ff )then
          proc = p
          exit ! do loop
       endif
    enddo
    ASSERT(proc>0)
  end function WhoseFitIndex

  function MaxFitIndex_Which(typ) result(indx)
    use error_module
    use comm_module
    implicit none
    integer(IK),intent(in) :: typ
    integer(IK)            :: indx !<<< result
    ! *** end of interface ***

    integer(IK)     :: id
    type(SplitDims) :: sd

    ASSERT(Set(typ))

    id = comm_myindex()

    sd = List(typ)

    if(id<sd%first.or.id>sd%last)&
         & call error("is23s/MaxFitIndex_Which: nothing for you")

    indx = sd%dims(id)
  end function MaxFitIndex_Which

  function MaxFitIndex_ChFit() result(indx)
    use error_module
    use comm_module
    implicit none
    integer(IK)            :: indx !<<< result
    ! *** end of interface ***

    integer(IK)     :: id
    type(SplitDims) :: sd

    ASSERT(Set(ix_ChFit))

    id = comm_myindex()

    sd = List(ix_ChFit)

    if(id<sd%first.or.id>sd%last)&
         & call error("is23s/MaxFitIndex_ChFit: nothing for you")

    indx = sd%dims(id)
  end function MaxFitIndex_ChFit

  subroutine SyncBorders(borders_ch,borders,typ)
    ! 
    use fit_coeff_module, only: ff_map=>fit_coeff_ff_map
    implicit none
    integer(IK),intent(in)    :: borders_ch(:,:)
    integer(IK),intent(inout) :: borders(:,:) !...(3,n_hosts)
    ! 1. dim: 1: lower border, 2: upper border, 3: number of ff
    integer(IK),intent(in)    :: typ
    ! *** end of interface ***

    integer(IK) :: np,p,start,stopp,last_stopp,nff,LL

    np = size(borders_ch,2)
    ASSERT(np==size(borders,2))

    nff = sum(borders_ch(3,:))
    ASSERT(nff==size(ff_map))

    last_stopp = 0
    do p=1,np
       start = borders_ch(1,p)
       stopp = borders_ch(2,p)
       nff   = stopp - start + 1
       ASSERT(nff==borders_ch(3,p))
       ! FIXME: can it be different?
       ASSERT(start==last_stopp+1)
       last_stopp = stopp

       select case(typ)
       case(IX_CHFIT_S)
          LL = -1 ! s-type
       case(IX_CHFIT_R2)
          LL =  0 ! r2-type
       case default
          ABORT('no such case yet')
       end select

       ! count how many LL-type functions fit
       ! into the range of ALL assigned to processor
       ! p.
       nff = count(ff_map(start:stopp)%L==LL)

       borders(3,p) = nff
    enddo

    last_stopp = 0
    do p=1,np
       borders(1,p) = last_stopp + 1
       borders(2,p) = last_stopp + borders(3,p)
       last_stopp   = last_stopp + borders(3,p)
    enddo
  end subroutine SyncBorders

  subroutine FitBorders(n_ff,borders)
    use error_module
    implicit none
    integer(IK),intent(in)    :: n_ff
    integer(IK),intent(inout) :: borders(:,:) !...(3,n_hosts)
    ! 1. dim: 1: lower border, 2: upper border, 3: number of ff
    ! *** end of interface ***

    integer(IK) :: first_host,last_host,n_hosts,n_per_host,rest_ff,&
         &         i_host,i_border

    if(size(borders,1)/=3)call error("is23s/FitBorders: borders shape?")

    first_host = lbound(borders,2)
    last_host  = ubound(borders,2)

    n_hosts    = last_host - first_host +1

    n_per_host = n_ff/n_hosts

    rest_ff = mod(n_ff,n_hosts)

    i_border = 1
    do i_host = first_host, first_host + rest_ff - 1
       borders(1,i_host) = i_border
       i_border = i_border + n_per_host
       borders(2,i_host) = i_border
       i_border = i_border + 1
       borders(3,i_host) = n_per_host + 1
    enddo
    do i_host = first_host + rest_ff, last_host
       borders(1,i_host) = i_border
       i_border = i_border + n_per_host - 1
       borders(2,i_host) = i_border
       i_border = i_border + 1
       borders(3,i_host) = n_per_host
    enddo
  end subroutine FitBorders

  subroutine SetSplitting(first,last,dims,ix)
    use error_module
    implicit none
    integer(IK),intent(in) :: first,last
    integer(IK),intent(in) :: dims(first:last)
    integer(IK),intent(in) :: ix ! ix_ChFit,...
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: n_hosts

    n_hosts = last - first + 1
    call error(n_hosts/=size(dims),"is23s/SetSplitting: n_hosts ?")

    List(ix)%first = first
    List(ix)%last  = last
    allocate(List(ix)%dims(first:last),STAT=memstat)
    call error(memstat,"is23s/SetSplitting: alloc failed")

    List(ix)%dims(first:last) = dims(first:last)
    Set(ix) = .true.
  end subroutine SetSplitting

  subroutine free_SplitDims(d)
    use error_module
    implicit none
    type(SplitDims),intent(inout) :: d
    ! *** end of interfaace ***

    integer :: memstat

    d%first = -1
    d%last  = -1
    deallocate(d%dims,STAT=memstat)
    call error(memstat,"is23s/free_SplitDims: dealloc failed")
  end subroutine free_SplitDims
    
  !--------------- End of module -------------------------------------
end module int_send_2cob3c_spor

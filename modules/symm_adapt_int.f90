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
module symm_adapt_int
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
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  use error_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !=======================================================
  !
  !  DECLARATION OF THESE TYPES MOVED FROM INT_DATA_2COB3C_MODULE ...
  !

!!$  type, public :: symadapt_totsym_2c_int_type
!!$     ! to hold symmetry-adapted total symmetric 
!!$     ! two center integrals of one IRREP
!!$     real(RK), pointer :: int(:,:,:,:)
!!$     ! int(n_c2, n_c1, n_if2, n_if1)
!!$     ! for relativistic integrals, n_c = n_exp
!!$     real(RK), pointer :: int_imag(:,:,:,:)
!!$     ! int(n_c2, n_c1, n_if2, n_if1)
!!$     ! for relativistic integrals, n_c = n_exp
!!$     ! imaginary part of integrals in case of spin orbit
!!$  end type symadapt_totsym_2c_int_type

  type, public :: symadapt_nottotsym_2c_int_type
     ! to hold symmetry-adapted not total symmetric
     ! two center integrals of one IRREP
     real(RK), pointer :: int(:,:,:,:,:) => NULL()
     ! int(n_c2, n_c1, n_if2, n_if1, n_partner)
     ! for relativistic integrals, n_c = n_exp
  end type symadapt_nottotsym_2c_int_type

!!$  type, public :: symadapt_totsym_3c_int_type
!!$     ! to hold symmetry-adapted total symmetric 
!!$     ! three center integrals of one IRREP
!!$     real(RK), pointer :: int(:,:,:,:,:)
!!$     ! int(n_c2, n_c1, n_ff, n_if2, n_if1)
!!$     ! for relativistic integrals, n_c = n_exp
!!$     real(RK), pointer :: int_imag(:,:,:,:,:)
!!$     ! int(n_c2, n_c1, n_ff, n_if2, n_if1)
!!$     ! for relativistic integrals, n_c = n_exp
!!$     ! imaginary part of integrals in case of spin orbit
!!$  end type symadapt_totsym_3c_int_type

  type, public :: symadapt_nottotsym_3cint_type
     ! to hold symmetry-adapted not total symmetric
     ! three center integrals of one IRREP
     real(RK), pointer :: int(:,:,:,:,:,:) => NULL()
     ! int(n_c2, n_c1, n_ff, n_if2, n_if1, n_partner)
     ! for relativistic integrals, n_c = n_exp
  end type symadapt_nottotsym_3cint_type

  !=======================================================
  !
  ! SYMMETRY ADAPTED INTEGRALS >>>
  !
  ! two-center:
  !
  type,public ::  sa_int_block
     real(RK), pointer :: re(:,:,:,:) => NULL()
     real(RK), pointer :: int(:,:,:,:) => NULL()
     real(RK), pointer :: im(:,:,:,:) => NULL()
     real(RK), pointer :: int_imag(:,:,:,:) => NULL()
     ! re,im(n_c2, n_c1, n_if2, n_if1)
     ! re => int, im => int_imag
  end type sa_int_block
  !---------------------------------

  !---------------------------------
  !
  ! three center:
  !
  type,public ::  sa_3c_int_block
     real(RK), pointer :: re(:,:,:,:,:) => NULL()
     real(RK), pointer :: int(:,:,:,:,:) => NULL()
     real(RK), pointer :: im(:,:,:,:,:) => NULL()
     real(RK), pointer :: int_imag(:,:,:,:,:) => NULL()
     ! re,im(n_c2, n_c1, n_fit, n_if2, n_if1)
     ! re => int, im => int_imag
  end type sa_3c_int_block
  !---------------------------------

  !
  ! prim_*_int declarations and methods
  ! moved to a separate
  !           prim_int_store.f90
  ! module
  !

  !------------ Declaration of constants and variables ------------

  integer(IK),parameter,public ::&
       & DEFAULT_INT_BLOCK = 0,&
       & LL_INT_BLOCK      = 1,&
       & LS_INT_BLOCK      = 2,&
       & SL_INT_BLOCK      = 3,&
       & SS_INT_BLOCK      = 4


  !=======================================================
  !------------ Interface statements ---------------------

  interface typealloc ! public alloc
     module procedure ralloc_sa_int    ! use data from dimensions,
     module procedure ralloc_sa_3c_int ! only: SymDimSpor
  end interface

  interface alloc     ! private alloc
     module procedure alloc_sa_int_block
     module procedure alloc_sa_3c_int_block
  end interface

  interface typedealloc ! public dealloc
     module procedure dealloc_sa_int_block
     module procedure dealloc_sa_3c_int_block
  end interface

  interface dealloc     ! private dealloc
     module procedure dealloc_sa_int_block
     module procedure dealloc_sa_3c_int_block
  end interface

  !------------ public functions and subroutines ------------------
  public &
       & typealloc,&
       & typedealloc,&
       & Write2cIntBlock, Write3cIntBlock,&
       & Pack2cIntBlock,  Pack3cIntBlock,&
       & UnPackIntBlock

  public :: hconjug ! conjugation (hermitean)

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine Write2cIntBlock(blck,diag,th)
    use error_module
    use readwriteblocked_module,only:&
         & readwriteblocked_tapehandle,&
         & readwriteblocked_write
    implicit none
    type(sa_int_block), intent(in)    :: blck
    logical,            intent(in)    :: diag
    type(readwriteblocked_tapehandle),intent(inout) :: th
    ! *** end of interface ***

    real(RK),dimension(1:size(blck%re),1:2) :: buf
    integer(IK)                         :: n_int,i

    if(size(blck%re).eq.0) return !<<< compiler bug in the next line?
    
    call pck_2c(blck%re, blck%im, diag, buf, n_int)

    do i=1,2 ! re,im
       call readwriteblocked_write(buf(1:n_int,i),th)
    end do
  end subroutine Write2cIntBlock

  subroutine Write3cIntBlock(blck,bord,diag,th)
    use error_module
    use readwriteblocked_module,only:&
         & readwriteblocked_tapehandle,&
         & readwriteblocked_write
    implicit none
    type(sa_3c_int_block),intent(in) :: blck
    logical,intent(in)               :: diag
    integer(IK),intent(in)           :: bord(3)
    type(readwriteblocked_tapehandle),intent(inout) :: th
    ! *** end of interface ***

    real(RK),dimension(size(blck%re),2) :: buf
    integer(IK)                         :: n_int,i

    call pck_3c(blck%re, blck%im, bord, diag, buf, n_int)
    
    do i=1,2 ! re,im
       call readwriteblocked_write(buf(1:n_int,i),th)
    enddo
  end subroutine Write3cIntBlock

  subroutine Pack2cIntBlock(blck,diag)
    use xpack
    implicit none
    type(sa_int_block), intent(in)    :: blck
    logical,            intent(in)    :: diag
    ! *** end of interface ***

    real(RK),dimension(size(blck%re),2) :: buf
    integer(IK)                         :: n_int,i
#ifndef FPP_GFORTRAN_BUGS
    real(RK) :: chk
#endif

    call pck_2c(blck%re, blck%im, diag, buf, n_int)

    call pck(n_int)

    do i = 1, 2 ! re, im
       call pck(buf(1:n_int, i))

#ifndef FPP_GFORTRAN_BUGS       /* with 4.6 At -O1 checksum fails */
       chk = sum(buf(1:n_int, i))
       call pck(chk)
#endif
    enddo
  end subroutine Pack2cIntBlock

  subroutine Pack3cIntBlock(blck,bord,diag)
    use xpack
    implicit none
    type(sa_3c_int_block), intent(in) :: blck
    integer(IK),intent(in)            :: bord(3)
    logical,            intent(in)    :: diag
    ! *** end of interface ***

    real(RK),dimension(size(blck%re),2) :: buf
    integer(IK)                         :: n_int,i
#ifndef FPP_GFORTRAN_BUGS
    real(RK) :: chk
#endif

    call pck_3c(blck%re, blck%im, bord, diag, buf, n_int)

    call pck(n_int)

    do i = 1, 2 ! re, im
       call pck(buf(1:n_int, i))

#ifndef FPP_GFORTRAN_BUGS       /* with 4.6 At -O1 checksum fails */
       chk = sum(buf(1:n_int, i))
       call pck(chk)
#endif
    enddo
  end subroutine Pack3cIntBlock

  subroutine UnPackIntBlock(th)
    use error_module
    use xpack
    use readwriteblocked_module,only:&
         & readwriteblocked_tapehandle,&
         & readwriteblocked_write
    implicit none
    type(readwriteblocked_tapehandle),intent(inout) :: th
    ! *** end of interface ***

    real(RK),pointer :: buf(:)
    integer(IK)      :: n_int,memstat,i
#ifndef FPP_GFORTRAN_BUGS
    real(RK) :: chk1, chk2
#endif

    call upck(n_int)

    allocate(buf(n_int),STAT=memstat)
    ASSERT(memstat==0)

    do i = 1, 2 ! re, im
       call upck(buf)

#ifndef FPP_GFORTRAN_BUGS       /* with 4.6 At -O1 checksum fails */
       !
       ! FIXME: the difference is  small, 1e-16 .. 1e-13. I attributed
       ! it to the artifacts of "extended" precision in floating point
       ! arithmetics:
       !
       call upck(chk1)
       chk2 = sum(buf)          ! FIXME: force this number to 64 bits!
       if (chk1 /= chk2) then
          print *, MyID, "chk1 = ", chk1, "chk2 = ", chk2, "diff =", chk2 - chk1
          ABORT("checksum failed")
       endif
#endif

       call readwriteblocked_write(buf,th)
    enddo

    deallocate(buf,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine UnPackIntBlock

  subroutine pck_2c(re, im, diag, buf, n_int)
    use error_module
    implicit none
    real(RK), intent(in) :: re(:,:,:,:), im(:,:,:,:)
    logical,            intent(in)    :: diag
    real(RK),intent(inout)            :: buf(:,:)
    integer(IK),intent(out)           :: n_int
    ! *** end of interface ***

    integer(IK),parameter               :: perm(4) = (/2,1,4,3/)

    integer(IK) :: n_rad,n_ind,n_fun,f1,f2,ij
    integer(IK) :: rad1,rad2,ind1,ind2

    if(diag)then
       n_rad = size(re, 1) ! radial (contraction)
       n_ind = size(re, 3) ! independent

       n_fun = n_rad * n_ind

       n_int = n_fun*(n_fun+1)/2

       ij = 0
       do f1=0,n_fun-1
          do f2=0,f1

             rad2 = mod(f2,n_rad) + 1
             rad1 = mod(f1,n_rad) + 1
             ind2 = f2/n_rad      + 1
             ind1 = f1/n_rad      + 1

             ij = ij + 1
             buf(ij,1) = re(rad2, rad1, ind2, ind1)
             buf(ij,2) = im(rad2, rad1, ind2, ind1)
          enddo
       enddo
    else
       n_int = size(re)

       buf(:,1) = reshp(re)
       buf(:,2) = reshp(im)
    endif
    
  contains
    function reshp(a4d) result(a1d)
      real(RK),intent(in) :: a4d(:,:,:,:)
      real(RK)            :: a1d(size(a4d)) !<<< result
      ! *** end of interface ***

      integer(IK),parameter :: ord(4) = (/1,3,2,4/)
      integer(IK)           :: sh(4),shord(4)

      sh  = shape(a4d)
        shord= sh(ord)
!      a1d = pack(reshape(a4d,sh(ord),ORDER=ord),.true.)
      a1d = pack(reshape(a4d,shord   ,ORDER=ord),.true.)
    end function reshp
  end subroutine pck_2c

  subroutine pck_3c(re, im, bord, diag, buf, n_int)
    use error_module
    use readwriteblocked_module,only:&
         & readwriteblocked_tapehandle,&
         & readwriteblocked_write
    implicit none
    real(RK), intent(in) :: re(:,:,:,:,:), im(:,:,:,:,:)
    logical,intent(in)               :: diag
    integer(IK),intent(in)           :: bord(3)
    real(RK),intent(inout)           :: buf(:,:) ! (1:n_int+,2)
    integer(IK),intent(out)          :: n_int
    ! *** end of interface ***

    integer(IK),parameter               :: perm(5) = (/2,1,3,5,4/)

    integer(IK) :: n_rad,n_ind,n_fun,f1,f2,ijk_s,ijk_e
    integer(IK) :: rad1,rad2,ind1,ind2,n_ff,ff_s,ff_e


    if(bord(2)-bord(1)+1/=bord(3))call error("sai/Write3cIntBlock: borders?")

    n_ff = bord(3)
    ff_s = bord(1)
    ff_e = bord(2)

    if(diag)then
       n_rad = size(re, 1) ! radial (contraction)
       n_ind = size(re, 4) ! independent

       n_fun = n_rad * n_ind

       n_int = n_ff * n_fun*(n_fun+1)/2

       ijk_s = 1
       ijk_e = n_ff
       do f1=0,n_fun-1
          do f2=0,f1

             rad2 = mod(f2,n_rad) + 1
             rad1 = mod(f1,n_rad) + 1
             ind2 = f2/n_rad      + 1
             ind1 = f1/n_rad      + 1

             buf(ijk_s:ijk_e, 1) = re(rad2, rad1, ff_s:ff_e, ind2, ind1)
             buf(ijk_s:ijk_e, 2) = im(rad2, rad1, ff_s:ff_e, ind2, ind1)

             ijk_s = ijk_s + n_ff
             ijk_e = ijk_e + n_ff
          enddo
       enddo
    else
       n_int = size(re) / size(re, 3) * n_ff

       buf(1:n_int, 1) = reshp(re(:, :, ff_s:ff_e, :, :))
       buf(1:n_int, 2) = reshp(im(:, :, ff_s:ff_e, :, :))
    endif
    
  contains
    function reshp(a5d) result(a1d)
      real(RK),intent(in) :: a5d(:,:,:,:,:)
      real(RK)            :: a1d(size(a5d)) !<<< result
      ! *** end of interface ***

      integer(IK),parameter :: prm(5) = (/3,1,4,2,5/)
      integer(IK),parameter :: ord(5) = (/2,4,1,3,5/)
      integer(IK)           :: shp(5),shpprm(5)

      shp = shape(a5d)

!!!      a1d = pack(reshape(a5d,shp(prm),ORDER=ord),.true.)
        shpprm=shp(prm)
      a1d = pack(reshape(a5d,shpprm,ORDER=ord),.true.)
    end function reshp
  end subroutine pck_3c

  subroutine transp(sa)
    implicit none
    type(sa_int_block), intent(inout) :: sa(:) ! (n_irr)
    ! *** end of interface ***

    type(sa_int_block) :: bt

    integer(IK)           :: irr
    integer(IK)           :: sh(4)
    integer(IK),parameter :: ord(4)=(/2,1,4,3/)

    do irr = 1, size(sa)
       sh = shape(sa(irr)%re)

       call alloc(sh(2), sh(1), sh(4), sh(3), bt)
    
       if(product(sh)>0)then
          bt%re(:, :, :, :) = reshape(sa(irr)%re, shape(bt%re), ORDER=ord)
          bt%im(:, :, :, :) = reshape(sa(irr)%im, shape(bt%im), ORDER=ord)
       end if

       call dealloc(sa(irr))

       sa(irr) = bt
    enddo
  end subroutine transp

  subroutine conjug(sa)
    implicit none
    type(sa_int_block), intent(inout) :: sa(:) ! (n_irr)
    ! *** end of interface ***

    integer(IK) :: irr

    do irr = 1, size(sa)
       sa(irr)%im = - sa(irr)%im
    enddo
  end subroutine conjug

  subroutine hconjug(sa)
    implicit none
    type(sa_int_block), intent(inout) :: sa(:) ! (n_irr)
    ! *** end of interface ***

    call transp(sa)
    call conjug(sa)
  end subroutine hconjug

  !
  ! MEMORY MANAGEMENT FOR SYMM-ADAPT INT:
  !

  subroutine ralloc_sa_int(a2, a1, L2, L1, ne2, ne1, sa, what)
    use dimensions, only: SubSpaceDim, &
        SymDimSpor, SymDimSporL, SymDimSporS
    implicit none
    integer(IK), intent(in)          :: a2, a1, L2, L1, ne2, ne1
    type(sa_int_block), intent(inout) :: sa(:) ! (n_irr)
    integer(IK), intent(in), optional :: what
    ! *** end of interface ***

    integer(IK) :: n2, n1, n_irr, irr

    integer(IK) :: iwhat
    type(SubSpaceDim) :: d1, d2

    iwhat = DEFAULT_INT_BLOCK
    if(present(what)) iwhat = what

    select case(iwhat)
    case (DEFAULT_INT_BLOCK)
       !
       ! FIXME: are these copies or just refs?
       !
       d2 = SymDimSpor(a2)
       d1 = SymDimSpor(a1)
    case (LL_INT_BLOCK)
       d2 = SymDimSporL(a2)
       d1 = SymDimSporL(a1)
    case (LS_INT_BLOCK)  
       d2 = SymDimSporL(a2)
       d1 = SymDimSporS(a1)
    case (SL_INT_BLOCK)
       d2 = SymDimSporS(a2)
       d1 = SymDimSporL(a1)
    case (SS_INT_BLOCK)
       d2 = SymDimSporS(a2)
       d1 = SymDimSporS(a1)
    case default
       ABORT("no such case")
    end select

    ASSERT(d2%n_irr==d1%n_irr)

    n_irr = d2%n_irr
    ASSERT(n_irr==size(sa))

    do irr = 1, size(sa) ! n_irr
       n2 = d2%lm(L2, irr) ! n_indep 2
       n1 = d1%lm(L1, irr) ! n_indep 1
       call alloc(ne2, ne1, n2, n1, sa(irr))
    enddo
  end subroutine ralloc_sa_int

  subroutine ralloc_sa_3c_int(a2, a1, L2, L1, ne2, ne1, nff, sa)
    use dimensions, only: SymDimSpor
    implicit none
    integer(IK), intent(in) :: a2, a1, L2, L1, ne2, ne1, nff
    type(sa_3c_int_block), intent(inout) :: sa(:) ! (n_irr)
    ! *** end of interface ***

    integer(IK) :: n2, n1, n_irr, irr

    call error(SymDimSpor(a2)%n_irr /= SymDimSpor(a1)%n_irr,&
         & "sai/ralloc_sa_3c_int: n_irr ???")

    n_irr = SymDimSpor(a2)%n_irr
    ASSERT(n_irr==size(sa))

    do irr = 1, size(sa) ! n_irr
       n2 = SymDimSpor(a2)%lm(L2,irr) ! n_indep 2
       n1 = SymDimSpor(a1)%lm(L1,irr) ! n_indep 1
       call alloc(ne2, ne1, nff, n2, n1, sa(irr))
    enddo
  end subroutine ralloc_sa_3c_int

  subroutine alloc_sa_int_block(ne2,ne1,n2,n1,b)
    implicit none
    integer(IK),intent(in)           :: ne2,ne1,n2,n1
    type(sa_int_block),intent(inout) :: b
    ! *** end of interface ***

    integer            :: memstat
    real(RK),parameter :: zero = 0.0_rk

    allocate(b%int(ne2,ne1,n2,n1),b%int_imag(ne2,ne1,n2,n1),STAT=memstat)
    ASSERT(memstat==0)
    b%re => b%int
    b%im => b%int_imag
    b%re = zero
    b%im = zero
  end subroutine alloc_sa_int_block

  subroutine alloc_sa_3c_int_block(ne2,ne1,nff,n2,n1,b)
    implicit none
    integer(IK),intent(in)              :: ne2,ne1,nff,n2,n1
    type(sa_3c_int_block),intent(inout) :: b
    ! *** end of interface ***

    integer            :: memstat
    real(RK),parameter :: zero = 0.0_rk

    allocate(b%int(ne2,ne1,nff,n2,n1),b%int_imag(ne2,ne1,nff,n2,n1),STAT=memstat)
    ASSERT(memstat==0)
    b%re => b%int
    b%im => b%int_imag
    b%re = zero
    b%im = zero
  end subroutine alloc_sa_3c_int_block

  subroutine dealloc_sa_int_block(b)
    implicit none
    type(sa_int_block),intent(inout) :: b
    ! *** end of interface ***

    integer            :: memstat

    deallocate(b%int,b%int_imag,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine dealloc_sa_int_block

  subroutine dealloc_sa_3c_int_block(b)
    implicit none
    type(sa_3c_int_block),intent(inout) :: b
    ! *** end of interface ***

    integer            :: memstat

    deallocate(b%int,b%int_imag,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine dealloc_sa_3c_int_block

  !--------------- End of module ----------------------------------
end module symm_adapt_int

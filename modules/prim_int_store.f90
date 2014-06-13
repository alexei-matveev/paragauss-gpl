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
module prim_int_store
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
  USE_MEMLOG
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !
  ! types and methods originally coming from
  !          symm_adapt_int.f90
  ! module
  !

  !=======================================================
  !
  ! PRIMITIVE,[CONTRACTED] INTEGRALS >>>
  !
  ! kin,nuc,overlap ...:
  !
  type, public :: prim_int

     real(RK),pointer     :: prim(:,:,:,:)
     ! prim(n_exp,n_exp,2l+1,2l+1)

     !---------------------------------
     logical,pointer :: mask(:,:,:,:)
     ! non zero elements

     integer(IK)     :: shp(4) 
     ! shp = shape(prim)

     integer(IK)     :: n_rad
     ! number of radial (independent) functions
     ! n_rad = product(shp(1:2))
  end type prim_int
  !---------------------------------
  !
  ! P_vec * V x P_vec, P_vec:
  !
  type, public :: prim_vec_int

     real(RK),pointer     :: prim(:,:,:,:,:)
     ! prim(n_exp2,n_exp1,2l2+1,2l1+1,3)

     logical,pointer :: mask(:,:,:,:,:)
     ! non zero elements

     integer(IK)     :: shp(5) 
     ! shp = shape(prim)

     integer(IK)     :: n_rad
     ! number of radial (independent) functions
     ! n_rad = product(shp(1:2))
  end type prim_vec_int
  !---------------------------------
  !
  ! charge, exchange fit integrals:
  !
  type, public :: prim_3c_int

     real(RK),pointer        :: prim(:,:,:,:,:)
     ! prim(n_exp2,n_exp1,n_ff,2L2+1,2L1+1)

     !---------------------------------
     logical,pointer  :: mask(:,:,:,:,:)
     ! non zero elements

     integer(IK)      :: shp(5) 
     ! shp = shape(prim)

     integer(IK)      :: n_rad
     ! number of radial (independent) functions
     ! n_rad = product(shp(1:3))

     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     logical,pointer  :: rmask(:,:,:) ! =any(any(mask,4),4)

     integer(IK)      :: n_red        ! =count(rmask)

     real(RK),pointer :: pck(:,:,:) ! pck(n_red,nf2,nf1)
  end type prim_3c_int
  !---------------------------------
  !
  ! charge, exchange fit integrals (vector type):
  !
  type, public :: prim_3cv_int

     real(RK),pointer        :: prim(:,:,:,:,:,:)
     ! prim(n_exp2,n_exp1,n_ff,2L2+1,2L1+1,3)

     !---------------------------------
     logical,pointer  :: mask(:,:,:,:,:,:)
     ! non zero elements

     integer(IK)      :: shp(6)
     ! shp = shape(prim)

     integer(IK)      :: n_rad
     ! number of radial (independent) functions
     ! n_rad = product(shp(1:3))

  end type prim_3cv_int
  !---------------------------------


  !------------ Declaration of constants and variables ------------

  !=======================================================
  !------------ Interface statements ---------------------

  interface set_mask
     module procedure set_mask_prim_3c
     module procedure set_mask_prim_2c
     module procedure set_mask_prim_vec
  end interface

  interface tpack
     module procedure pack_prim_3c_int
  end interface

  interface free_mask
     module procedure free_mask_prim_3c
     module procedure free_mask_prim_2c
     module procedure free_mask_prim_vec
  end interface

  interface alloc
     module procedure alloc_prim_int
     module procedure alloc_prim_vec_int
     module procedure alloc_prim_3c_int
     module procedure alloc_prim_3cv_int
  end interface

  interface dealloc
     module procedure dealloc_prim_int
     module procedure dealloc_prim_vec_int
     module procedure dealloc_prim_3c_int
     module procedure dealloc_prim_3cv_int
  end interface 

  !------------ public functions and subroutines ------------------
  public &
       set_mask, &
       free_mask, &
       tpack, &
       alloc, &
       dealloc
       
  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine alloc_prim_int(ne2,ne1,nlm2,nlm1,p2c)
    implicit none
    integer(IK),   intent(in)    :: ne2,ne1,nlm2,nlm1
    type(prim_int),intent(inout) :: p2c
    ! *** end of interface ***

    integer(IK) :: memstat

    allocate(p2c%prim(ne2,ne1,nlm2,nlm1),stat=memstat)
    ASSERT(memstat==0)
    MEMLOG(+size(p2c%prim))
    p2c%shp = shape(p2c%prim)
  end subroutine alloc_prim_int

  subroutine dealloc_prim_int(p2c)
    implicit none
    type(prim_int),intent(inout) :: p2c
    ! *** end of interface ***

    integer(IK) :: memstat

    MEMLOG(-size(p2c%prim))
    deallocate(p2c%prim,stat=memstat)
    ASSERT(memstat==0)

    p2c%shp = -1
  end subroutine dealloc_prim_int

  subroutine alloc_prim_vec_int(ne2,ne1,nlm2,nlm1,p2c)
    implicit none
    integer(IK),      intent(in)     :: ne2,ne1,nlm2,nlm1
    type(prim_vec_int),intent(inout) :: p2c
    ! *** end of interface ***

    integer(IK) :: memstat

    allocate(p2c%prim(ne2,ne1,nlm2,nlm1,3),stat=memstat)
    ASSERT(memstat==0)
    MEMLOG(+size(p2c%prim))
    p2c%shp = shape(p2c%prim)
  end subroutine alloc_prim_vec_int

  subroutine dealloc_prim_vec_int(p2c)
    implicit none
    type(prim_vec_int),intent(inout) :: p2c
    ! *** end of interface ***

    integer(IK) :: memstat

    MEMLOG(-size(p2c%prim))
    deallocate(p2c%prim,stat=memstat)
    ASSERT(memstat==0)

    p2c%shp = -1
  end subroutine dealloc_prim_vec_int

  subroutine alloc_prim_3c_int(ne2,ne1,n3c,nlm2,nlm1,p3c)
    implicit none
    integer(IK),      intent(in)    :: ne2,ne1,n3c,nlm2,nlm1
    type(prim_3c_int),intent(inout) :: p3c
    ! *** end of interface ***

    integer(IK) :: memstat

    allocate(p3c%prim(ne2,ne1,n3c,nlm2,nlm1),stat=memstat)
    ASSERT(memstat==0)
    MEMLOG(+size(p3c%prim))
    p3c%shp = shape(p3c%prim)
  end subroutine alloc_prim_3c_int

  subroutine dealloc_prim_3c_int(p3c)
    implicit none
    type(prim_3c_int),intent(inout) :: p3c
    ! *** end of interface ***

    integer(IK) :: memstat

    MEMLOG(-size(p3c%prim))
    deallocate(p3c%prim,stat=memstat)
    ASSERT(memstat==0)

    p3c%shp = -1
  end subroutine dealloc_prim_3c_int

  subroutine alloc_prim_3cv_int(ne2,ne1,n3c,nlm2,nlm1,p3c)
    implicit none
    integer(IK),      intent(in)     :: ne2,ne1,n3c,nlm2,nlm1
    type(prim_3cv_int),intent(inout) :: p3c
    ! *** end of interface ***

    integer(IK) :: memstat

    allocate(p3c%prim(ne2,ne1,n3c,nlm2,nlm1,3),stat=memstat)
    ASSERT(memstat==0)
    MEMLOG(+size(p3c%prim))
    p3c%shp = shape(p3c%prim)
  end subroutine alloc_prim_3cv_int

  subroutine dealloc_prim_3cv_int(p3c)
    implicit none
    type(prim_3cv_int),intent(inout) :: p3c
    ! *** end of interface ***

    integer(IK) :: memstat

    MEMLOG(-size(p3c%prim))
    deallocate(p3c%prim,stat=memstat)
    ASSERT(memstat==0)

    p3c%shp = -1
  end subroutine dealloc_prim_3cv_int

  !=======================================================
  !
  ! MASKING NON-ZERO ELEMENTS >>>
  !

  !---------------------------------
  !
  ! 3-CENTER INTEGRALS
  !

  subroutine set_mask_prim_3c(pint)
    use options_module, only: options_integral_expmax
    implicit none
    type(prim_3c_int),intent(inout) :: pint
    ! *** end of interface ***

    integer :: memstat
    integer(IK)         :: shp(5)
    real(RK),parameter  :: ZERO = 0.0_rk

    call error(.not.associated(pint%prim),"set_mask_prim_3c: 1")

    pint%shp   = shape(pint%prim)

    pint%n_rad = product(pint%shp(1:3))

    shp = pint%shp
    allocate(&
         & pint%mask(shp(1),shp(2),shp(3),shp(4),shp(5)), &
         & pint%rmask(shp(1),shp(2),shp(3)),              &
         & STAT=memstat)
    call error(memstat,"set_mask_prim_3c: 2")

    !debug>>>
!!$    ! it`s easier to find "zeros" :
!!$    !
!!$    pint%mask = (pint%prim .eq. ZERO)
!!$
!!$    pint%mask = pint%mask .or. &
!!$         & (2*exponent(pint%prim)/3 < -INT(options_integral_expmax()))
!!$    !
!!$    ! since ln(2) ~= 3/2
!!$    !
!!$
!!$    pint%mask = .not. pint%mask
    !<<<debug

    pint%mask = (pint%prim .ne. ZERO)

!!$    print *,'sam/set_mask_prim_3c: of ',size(pint%mask),&
!!$         & ' (  all ) zeros make ',count(.NOT.pint%mask),&
!!$         & ' ratio=',100*count(.NOT.pint%mask)/size(pint%mask)
    
    pint%rmask = any(any(pint%mask,DIM=4),DIM=4)

    pint%n_red = count(pint%rmask)

!!$    print *,'sam/set_mask_prim_3c: of ',size(pint%rmask),&
!!$         & ' (radial) zeros make ',count(.NOT.pint%rmask),&
!!$         & 'ratio=',100*count(.NOT.pint%rmask)/size(pint%rmask)

  end subroutine set_mask_prim_3c

  subroutine pack_prim_3c_int(p)
    implicit none
    type(prim_3c_int),intent(inout) :: p
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: m2,m1
  
 !  print*,'I am here in pack_prim_3c_int 1'

    call error(any(p%shp/=shape(p%prim)),"sam/pack_prim_3c_int: shape wrong")

    allocate(p%pck(p%n_red,p%shp(4),p%shp(5)),STAT=memstat)
    call error(memstat,"sam/pack_prim_3c_int: alloc failed")

    if(p%n_red==0) return !<<< for HP

    do m2=1,p%shp(4)
       do m1=1,p%shp(5)
          p%pck(:,m2,m1) = pack(p%prim(:,:,:,m2,m1),p%rmask)
       enddo
    enddo
  end subroutine pack_prim_3c_int

  ! ^^^ 3-c integrals
  !---------------------------------

  subroutine set_mask_prim_2c(pint)
    implicit none
    type(prim_int),intent(inout) :: pint
    ! *** end of interface ***

    integer :: memstat
    integer(IK)         :: shp(4)
    real(RK),parameter  :: ZERO = 0.0_rk

    call error(.not.associated(pint%prim),"set_mask_prim_2c: 1")

    pint%shp   = shape(pint%prim)

    pint%n_rad = product(pint%shp(1:2))

    shp = pint%shp
    allocate(&
         & pint%mask(shp(1),shp(2),shp(3),shp(4)),&
         & STAT=memstat)
    call error(memstat,"set_mask_prim_2c: 2")

    pint%mask = (pint%prim /= ZERO)

!!$    print *,'set_mask_prim_2c: of ',size(pint%mask),&
!!$         & ' zeros make ',count(.NOT.pint%mask)
  end subroutine set_mask_prim_2c

  subroutine set_mask_prim_vec(pint)
    implicit none
    type(prim_vec_int),intent(inout) :: pint
    ! *** end of interface ***

    integer :: memstat
    integer(IK)         :: shp(5)
    real(RK),parameter  :: ZERO = 0.0_rk

    call error(.not.associated(pint%prim),"set_mask_prim_vec: 1")

    pint%shp   = shape(pint%prim)

    pint%n_rad = product(pint%shp(1:2))

    shp = pint%shp
    allocate(&
         & pint%mask(shp(1),shp(2),shp(3),shp(4),shp(5)),&
         & STAT=memstat)
    call error(memstat,"set_mask_prim_vec: 2")

    pint%mask = (pint%prim /= ZERO)

!!$    print *,'set_mask_prim_vec: of ',size(pint%mask),&
!!$         & ' zeroes make ',count(.NOT.pint%mask)
  end subroutine set_mask_prim_vec

  subroutine free_mask_prim_3c(pint,PCK)
    implicit none
    type(prim_3c_int),intent(inout) :: pint
    logical,optional,intent(in)     :: PCK
    ! *** end of interface ***

    integer :: memstat
    logical :: pck_

    pck_ = .false.
    if(present(PCK)) pck_ = PCK

    deallocate(pint%mask,STAT=memstat)
    call error(memstat,"free_mask_prim_3c: 1")

    if(pck_)then
       deallocate(pint%pck,STAT=memstat)
       call error(memstat,"free_mask_prim_3c: 2")
    endif
  end subroutine free_mask_prim_3c

  subroutine free_mask_prim_2c(pint)
    implicit none
    type(prim_int),intent(inout) :: pint
    ! *** end of interface ***

    integer :: memstat

    deallocate(pint%mask,STAT=memstat)
    call error(memstat,"free_mask_prim_2c")
  end subroutine free_mask_prim_2c

  subroutine free_mask_prim_vec(pint)
    implicit none
    type(prim_vec_int),intent(inout) :: pint
    ! *** end of interface ***
    
    integer :: memstat

    deallocate(pint%mask,STAT=memstat)
    call error(memstat,"free_mask_prim_vec")
  end subroutine free_mask_prim_vec

  !--------------- End of module ----------------------------------
end module prim_int_store

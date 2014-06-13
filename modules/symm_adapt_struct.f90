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
module symm_adapt_struct
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
  use dimensions, only: SubSpaceDim
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !=======================================================
  !
  ! The third, I think, struct for symmetry 
  ! adaption in THE PROGRAM :
  !

  type, public :: symm_adapt_ul   ! a subspace {UA,L}, equiv to that above

     integer(IK)                     :: n_ea,lmax,n_irr
     type(symm_adapt_irrbas),pointer :: aos(:,:,:)
     logical,pointer                 :: exists(:,:)
     !
     ! A_tomic O_rbitals :: aos(n_equiv_atoms,0:Lmax,n_irr)
     !                   :: exists(0:Lmax,n_irr)
  end type symm_adapt_ul
  !---------------------------------

  type,public :: symm_adapt_irrbas ! Irr_ep BAS_is
     !
     ! contributions to SALCAOs centered on an atom
     ! for each partner and independent function
     !
     integer(IK)                     :: dim_irr,n_indep
     type(symm_adapt_spinor),pointer :: bas(:,:)
     !
     ! Irr_ep BAS_is, bas(dim_irr,n_indep)
  end type symm_adapt_irrbas
  !---------------------------------

  type,public :: symm_adapt_ao ! A_tomic O_rbital
     !
     ! (complex) combination of |lm>
     !
     integer(IK)          :: n
     integer(IK),pointer  :: m(:)  ! m(n)
     real(RK),pointer     :: re(:) ! re(n)
     real(RK),pointer     :: im(:) ! im(n)
     !---------------------------------
     real(RK),pointer     :: c(:,:) ! c(:,2)
     ! re => c(:,1), im => c(:,2)
  end type symm_adapt_ao
  !---------------------------------

  type,public :: symm_adapt_spinor
     !
     ! 2-spinor
     !
     type(symm_adapt_ao) :: psi(2) 
  end type symm_adapt_spinor
  !---------------------------------

  type(symm_adapt_ul),public,pointer :: LSymAdp(:),LSymAdpL(:),LSymAdpS(:) ! LSymAdp(n_ua)

  !=======================================================
  !
  ! SIGMA-MATRIX REPRESENTATION:
  !
  type,public :: symm_adapt_4x4mat

     !--------- natural: --------------      !     | 1  0 |
     integer(IK)              :: m(4,4)      ! 1 = |      |
     !--------- packed: ---------------      !     | 0  1 |
     integer(IK),dimension(4) :: alpha,beta  !
     integer(IK),dimension(4) :: a,b         !     | 0 -1 |
     integer(IK),dimension(4) :: c           ! i = |      |
  end type symm_adapt_4x4mat                 !     | 1  0 |


  !=======================================================
  !------------ Declaration of constants and variables ---

  integer(IK),private            :: E(4,4),IE(4,4),SIGMA(4,4,3),i__,j__
  type(symm_adapt_4x4mat),public :: tE(2), tSIGMA(2,3)

  data E(:,:)         /&
       &  1, 0,  0, 0, &
       &  0, 1,  0, 0, &
       &  0, 0,  1, 0, &
       &  0, 0,  0, 1  /  

  data ((IE(i__,j__),j__=1,4),i__=1,4) /&
       &  0,-1,  0, 0,&
       &  1, 0,  0, 0,&
       &  0, 0,  0,-1,&
       &  0, 0,  1, 0 /

  data SIGMA(:,:,:)   /&
       &  0, 0,  1, 0, &
       &  0, 0,  0, 1, &
       &  1, 0,  0, 0, &
       &  0, 1,  0, 0  &
       &             , &
       &  0, 0,  0,-1, &
       &  0, 0,  1, 0, &
       &  0, 1,  0, 0, &
       & -1, 0,  0, 0  &
       &             , &
       &  1, 0,  0, 0, &
       &  0, 1,  0, 0, &
       &  0, 0, -1, 0, &
       &  0, 0,  0,-1  /
  !---------------------------------

  !=======================================================
  !------------ Interface statements ---------------------

  interface alloc     ! private alloc
     module procedure alloc_symm_adapt_ul
     module procedure alloc_symm_adapt_irrbas
     module procedure alloc_symm_adapt_ao
  end interface

  interface dealloc     ! private dealloc
     module procedure dealloc_symm_adapt_ul
     module procedure dealloc_symm_adapt_irrbas
     module procedure dealloc_symm_adapt_ao
  end interface

  interface free       ! recursive dealloc
     module procedure free_symm_adapt_ul
     module procedure free_arr_symm_adapt_ul
     module procedure free_symm_adapt_irrbas
  end interface

  interface tpack
     module procedure pack_symm_adapt_4x4mat
  end interface

  interface tprint
     module procedure print_4x4mat
     module procedure print_symm_adapt_4x4mat
  end interface

  interface typeprint
     module procedure print_4x4mat
     module procedure print_symm_adapt_4x4mat
  end interface

  !------------ public functions and subroutines ------------------

  public init_symm_adapt_struct,&
       & done_symm_adapt_struct,&
       & alloc,free,&
       & typeprint !<<< debug only

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init_symm_adapt_struct()
    implicit none
    ! *** end of interface ***

    call pack_Sigma(E,IE,SIGMA,tE,tSIGMA)
  end subroutine init_symm_adapt_struct

  subroutine done_symm_adapt_struct()
    implicit none
    ! *** end of interface ***

  end subroutine done_symm_adapt_struct

  subroutine pack_Sigma(E,I,S,tE,tS)
    implicit none
    integer(IK),intent(in)                :: E(4,4),I(4,4),S(4,4,3)
    type(symm_adapt_4x4mat),intent(inout) :: tE(2),tS(2,3)
    ! *** end of interface ***

    integer(IK) :: k

    call tpack( E,tE(1))
    call tpack(-I,tE(2))

!!$    if(debug)then
!!$       print *,'---------- Product (Re) ---------'
!!$       call tprint(tE(1))
!!$       print *,'---------- Product (Im) ---------'
!!$       call tprint(tE(2))
!!$    endif

    do k=1,3

       call tpack(matmul( E,S(:,:,k)),tS(1,k))
       call tpack(matmul(-I,S(:,:,k)),tS(2,k))

!!$       if(debug)then
!!$          print *,'----- SIGMA(',k,'), (Re) ------'
!!$          call tprint(tS(1,k))
!!$
!!$          print *,'----- SIGMA(',k,'), (Im) ------'
!!$          call tprint(tS(2,k))
!!$       endif
    enddo
  end subroutine pack_Sigma

  subroutine pack_symm_adapt_4x4mat(m, tm)
    implicit none
    integer(IK),intent(in)                :: m(:, :) ! (4, 4)
    type(symm_adapt_4x4mat),intent(inout) :: tm
    ! *** end of interface ***

    logical :: msk(4,4)

    ASSERT(size(m,1)==4)
    ASSERT(size(m,2)==4)

    tm%m = m

    msk = (m /= 0)

    tm%c = pack(m,msk)

    tm%alpha = pack(spread((/1,1,2,2 /),2,4),msk)
    tm%beta  = pack(spread((/1,1,2,2 /),1,4),msk)
    tm%a     = pack(spread((/1,2,1,2 /),2,4),msk)
    tm%b     = pack(spread((/1,2,1,2 /),1,4),msk)
  end subroutine pack_symm_adapt_4x4mat

  subroutine print_4x4mat(M)
    implicit none
    integer(IK),intent(in) :: M(4,4)
    ! *** end of interface ***

    integer(IK) :: i,j

    write(*,*)
    write(*,'(4I3)') ((M(i,j),j=1,4),i=1,4)
    write(*,*)
  end subroutine print_4x4mat

  subroutine print_symm_adapt_4x4mat(tm)
    implicit none
    type(symm_adapt_4x4mat),intent(in) :: tm
    ! *** end of interface ***

    write(*,'(A)') 'unpacked >>>'

    call tprint(tm%m)
    
    write(*,'(A)') 'packed   >>>'

    write(*,*)
    write(*,'("alpha | ",4I3)') tm%alpha
    write(*,'("beta  | ",4I3)') tm%beta
    write(*,'("a     | ",4I3)') tm%a
    write(*,'("b     | ",4I3)') tm%b
    write(*,'("------+-------------")')
    write(*,'("c     | ",4I3)') tm%c
    write(*,*)   
  end subroutine print_symm_adapt_4x4mat
    

  !-------------------------------------------------------
  !
  ! FOR TYPE(SYMM_ADAPT_XXX) >>>
  !

  subroutine free_symm_adapt_ul(ul)
    implicit none
    type(symm_adapt_ul),intent(inout) :: ul
    ! *** end of interface ***
    
    integer(IK) :: ea,l,irr

    l_:do l=1,ul%lmax
       irr_:do irr=1,ul%n_irr

          if(ul%exists(l,irr))then
             ea_:do ea=1,ul%n_ea
                call free(ul%aos(ea,l,irr))
             enddo ea_
          endif

       enddo irr_
    enddo l_
    call dealloc(ul)
  end subroutine free_symm_adapt_ul

  subroutine free_arr_symm_adapt_ul(ul)
    implicit none
    type(symm_adapt_ul),intent(inout) :: ul(:)
    ! *** end of interface ***
    
    integer(IK) :: i

    do i=1,size(ul)
       call free(ul(i))
    enddo
  end subroutine free_arr_symm_adapt_ul

  subroutine free_symm_adapt_irrbas(b)
    implicit none
    type(symm_adapt_irrbas),intent(inout) :: b
    ! *** end of interface ***

    integer(IK) :: p,f,c

    p_: do p=1,b%dim_irr
       f_: do f=1,b%n_indep
          c_: do c=1,2
             call dealloc(b%bas(p,f)%psi(c))
          enddo c_
       enddo f_
    enddo p_
    call dealloc(b)
  end subroutine free_symm_adapt_irrbas
    
  !-------------------------------------------------------

  subroutine alloc_symm_adapt_ul(n_ea,lmax,n_irr,ua)
    use error_module
    implicit none
    integer(IK),intent(in)            :: n_ea,lmax,n_irr
    type(symm_adapt_ul),intent(inout) :: ua
    ! *** end of interface ***

    integer :: memstat

    ua%n_ea  = n_ea
    ua%lmax  = lmax
    ua%n_irr = n_irr

    allocate(&
         & ua%aos(n_ea,0:lmax,n_irr),&
         & ua%exists(0:lmax,n_irr),&
         & STAT=memstat)
    call error(memstat,"sas/alloc_symm_adapt_ul: alloc failed")

    ua%exists = .false.
  end subroutine alloc_symm_adapt_ul

  subroutine dealloc_symm_adapt_ul(ua)
    use error_module
    implicit none
    type(symm_adapt_ul),intent(inout) :: ua
    ! *** end of interface ***

    integer :: memstat

    ua%n_ea  = -1
    ua%lmax  = -1
    ua%n_irr = -1

    deallocate(ua%aos,ua%exists, STAT=memstat)
    call error(memstat,"sas/dealloc_symm_adapt_ul: alloc failed")
  end subroutine dealloc_symm_adapt_ul

  subroutine alloc_symm_adapt_ao(n,o)
    use error_module
    implicit none
    integer(IK),intent(in)            :: n
    type(symm_adapt_ao),intent(inout) :: o
    ! *** end of interface ***

    integer :: memstat

    o%n = n

    allocate(o%m(n), o%c(n,2), STAT=memstat)
    call error(memstat,"sas/alloc_symm_adapt_ao: alloc failed")

    o%re => o%c(:,1)
    o%im => o%c(:,2)
  end subroutine alloc_symm_adapt_ao

  subroutine dealloc_symm_adapt_ao(o)
    use error_module
    implicit none
    type(symm_adapt_ao),intent(inout) :: o
    ! *** end of interface ***

    integer :: memstat

    o%n = -1

    deallocate(o%m, o%c, STAT=memstat)
    call error(memstat,"sas/dealloc_symm_adapt_ao: alloc failed")

    nullify(o%re,o%im)
  end subroutine dealloc_symm_adapt_ao

  subroutine alloc_symm_adapt_irrbas(dim_irr,n_indep,b)
    use error_module
    implicit none
    integer(IK),intent(in)                :: dim_irr,n_indep
    type(symm_adapt_irrbas),intent(inout) :: b
    ! *** end of interface ***

    integer :: memstat

    b%dim_irr = dim_irr
    b%n_indep = n_indep

    allocate(b%bas(dim_irr,n_indep), STAT=memstat)
    call error(memstat,"sas/alloc_symm_adapt_bas: alloc failed")
  end subroutine alloc_symm_adapt_irrbas

  subroutine dealloc_symm_adapt_irrbas(b)
    use error_module
    implicit none
    type(symm_adapt_irrbas),intent(inout) :: b
    ! *** end of interface ***

    integer :: memstat

    b%dim_irr = -1
    b%n_indep = -1

    deallocate(b%bas,STAT=memstat)
    call error(memstat,"sas/dealloc_symm_adapt_irrbas: dealloc failed")
  end subroutine dealloc_symm_adapt_irrbas

  !--------------- End of module ----------------------------------
end module symm_adapt_struct

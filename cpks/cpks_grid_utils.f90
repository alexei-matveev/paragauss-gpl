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
module cpks_grid_utils
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
  use type_module, only:   &
                   i4_kind &
                 , r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: gridInt
  public :: gridInt2
  public :: rho1
  public :: rho2
  public :: eigMO
  public :: symWtsNGra
  public :: cartNGra
  public :: cartNDer

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function rho1( vl, o, phi, u1 ) result(r1)
    !
    ! Compute density response r1 due to orbitals shifts U1(i,j)
    !
    !   r1 = SUM[i,a] { 2 * phi(i) * phi(a) *     U1(i,a)       }
    !
    ! note also the use for occ x occ response:
    !
    !   r1 = SUM[i,j] { 2 * phi(i) * phi(j) * ( - S1(i,j) / 2 ) }
    !
    implicit none
    integer(i4_kind), intent(in) :: vl, o
    real(r8_kind)   , intent(in) :: phi(:,:)   ! (:vl,:n)
    real(r8_kind)   , intent(in) :: u1(1:,o+1:) ! (ni,ni) or (ni,ni+1:ni+na)
    real(r8_kind)                :: r1(vl)     ! result
    ! *** end of interface ***

    r1 = rho2( vl, o, phi, phi, u1 )
  end function rho1

  function rho2( vl, o, phi, chi, u1 ) result(r1)
    !
    ! Compute response of the density gradient r1
    ! due to orbitals shifts U1(i,j)
    !
    !   r1 = SUM[i,a] { 2 * phi(i) * chi(a) *     U1(i,a)       }
    !
    ! note also the use for occ x occ response:
    !
    !   r1 = SUM[i,j] { 2 * phi(i) * chi(j) * ( - SX(i,j) / 2 ) }
    !
    ! with chi = d/dY phi and SX = d/dX S
    !
    implicit none
    integer(i4_kind), intent(in) :: vl, o
    real(r8_kind)   , intent(in) :: phi(:,:)    ! (:vl,:n)
    real(r8_kind)   , intent(in) :: chi(:,:)    ! (:vl,:n)
    real(r8_kind)   , intent(in) :: u1(1:,o+1:) ! (ni,ni) or (ni,ni+1:ni+na)
    real(r8_kind)                :: r1(vl)      ! result
    ! *** end of interface ***

    integer(i4_kind) :: n,ni,na
#ifdef FPP_NOBLAS
    integer(i4_kind) :: i,a
    real(r8_kind)    :: pre(vl)
#else
    integer(i4_kind) :: i
    real(r8_kind)    :: mm(vl,size(u1,1)) ! (vl,ni)
#endif

    n  = size(phi,2)
    ni = size(u1,1)
    na = size(u1,2)

#if SB_COMMENT
    if( o /= 0 )then
      ASSERT(o==ni)
      ASSERT(n==ni+na)
      ASSERT(n==size(chi,2))
    else
      ASSERT(na==ni)
      ASSERT(ni<=n)
      ASSERT(ni<=size(chi,2))
    endif
#endif

    ASSERT(vl<=size(phi,1))
    ASSERT(vl<=size(chi,1))

!     do i=1,ni
!     do a=o+1,o+na
!     do p=1,vl
!       r1(p) += 2 * phi(p,i) * u1(i,a) * chi(p,a)

#ifdef FPP_NOBLAS
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning
    r1(:) = 0.0_r8_kind
    do i=1,ni
      pre(:) = 0.0_r8_kind
      do a=o+1,o+na
        pre(:) = pre(:) + u1(i,a) * chi(:vl,a)
      enddo
      r1(:) = r1(:) + 2 * phi(:vl,i) * pre(:)
    enddo
#else
    ! mm(p,i) := SUM[a]{ 2 * chi(p,a) * u1(i,a) }:
    !         o+1 <= a <= o+na
    call dgemm( 'n','t', vl, ni, na,     2.0_r8_kind &
              , chi(1,o+1), size(chi,1)              &
              , u1        , size(u1 ,1), 0.0_r8_kind &
              , mm        , size(mm ,1)              &
              )
    r1(:) = 0.0_r8_kind
    do i=1,ni
      r1(:) = r1(:) + phi(:vl,i) * mm(:,i)
    enddo
#endif
  end function rho2

  subroutine eigMO(vl,o,ob,ev,mo)
    !
    ! mo(:vl,m) += SUM[k] ob(:vl,k) * ev(k,m)
    !
    ! sum over subrange of basis functions:
    !
    !         o + 1 <= k <= o + nk
    !
!   use f77_blas, only: dgemm
    implicit none
    integer(i4_kind), intent(in)  :: vl, o
    real(r8_kind)   , intent(in)  :: ob(:,o+1:) ! (:vl,o+1:o+nk)
    real(r8_kind)   , intent(in)  :: ev(:,:)    ! (  n,=nmo)
    real(r8_kind)   , intent(out) :: mo(:,:)    ! (:vl,=nmo)
    ! *** end of interface ***

    integer(i4_kind) :: nmo,nk
#ifdef FPP_NOBLAS
    integer(i4_kind) :: m,i,k
#endif

    nmo = size(mo,2)
    nk  = size(ob,2)

    ASSERT(o+nk<=size(ev,1))
    ASSERT(nmo ==size(ev,2))

#ifdef FPP_NOBLAS
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning
    ! DONT: mo = 0.0_r8_kind
!   if( o /= 0 )then
      do m=1,nmo    ! molecular orbitlas
      do k=o+1,o+nk ! range of basis functions (on single ua)
      do i=1,vl     ! grid points
        mo(i,m) = mo(i,m) + ob(i,k) * ev(k,m)
      enddo
      enddo
      enddo
!   else
!     call dgemm( 'n','n', vl, nmo, nk, 1.0_r8_kind &
!               , ob, size(ob,1)                    &
!               , ev, size(ev,1),       1.0_r8_kind &
!               , mo, size(mo,1)                    &
!               )
!   endif
#else
!error "First make the BLAS version work for o /= 0 !"
    ! DONT: mo = 0.0_r8_kind, note beta==1
!   ! FIXME: copy-in-out:
!   call dgemm( 'n','n', vl, nmo, nk,   1.0_r8_kind &
!             , ob        , size(ob,1)              &
!             , ev(o+1:o+nk,:), nk    , 1.0_r8_kind &
!             , mo        , size(mo,1)              &
!             )
    call dgemm( 'n','n', vl, nmo, nk,  1.0_r8_kind &
              , ob       , size(ob,1)              &
              , ev(o+1,1), size(ev,1), 1.0_r8_kind &
              , mo       , size(mo,1)              &
              )
#endif
  end subroutine eigMO

  subroutine gridInt(vl,o,vx,phi,q)
    !
    ! Integrates over grid to obtain matrix elements
    !
    !   q(i,a) += < i | vx | a > =
    !
    !           = SUM[r]{ phi(r,i) * vx(r) * phi(r,a) }
    !
    implicit none
    integer(i4_kind), intent(in)    :: vl, o
    real(r8_kind)   , intent(in)    :: vx(:)      ! (:vl)
    real(r8_kind)   , intent(in)    :: phi(:,:)   ! (:vl,ni+na)
    real(r8_kind)   , intent(inout) :: q(1:,o+1:) ! (1:ni,o+1:o+na)
    ! *** end of interface ***

    call gridInt2(vl,o,vx,phi,phi,q)
  end subroutine gridInt

  subroutine gridInt2(vl,o,vx,phi,chi,q)
    !
    ! Integrates over grid to obtain matrix elements
    !
    !   q(i,a) += < i | vx | a > =
    !
    !           = SUM[r]{ phi(r,i) * vx(r) * chi(r,a) }
    !
!   use f77_blas, only: dgemm
    implicit none
    integer(i4_kind), intent(in)    :: vl, o
    real(r8_kind)   , intent(in)    :: vx(:)      ! (:vl)
    real(r8_kind)   , intent(in)    :: phi(:,:)   ! (:vl,ni+na)
    real(r8_kind)   , intent(in)    :: chi(:,:)   ! (:vl,ni+na)
    real(r8_kind)   , intent(inout) :: q(1:,o+1:) ! (ni,o+1:o+na)
    ! *** end of interface ***

    integer(i4_kind) :: n,ni,na
!define FPP_NOBLAS
#ifdef FPP_NOBLAS
    integer(i4_kind) :: i,a,r
#else
    integer(i4_kind) :: i
    real(r8_kind)    :: phiv(vl,size(q,1)) ! (vl,ni)
#endif

    n  = size(phi,2)
    ni = size(q,1)
    na = size(q,2)

#if SB_COMMENT
    if( o /= 0 )then
      ASSERT(o==ni)
      ASSERT(n==ni+na)
    else
      ASSERT(na==ni)
      ASSERT(ni<=n)
    endif
#endif

    ASSERT(vl<=size(phi,1))
    ASSERT(vl<=size(chi,1))

#if SB_COMMENT
    ASSERT(n ==size(chi,2))
#endif

#ifdef FPP_NOBLAS
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning
    ! F77 VERSION:
    do i=1,ni
    do a=o+1,o+na
      ! DONT: q(i,a) = 0
      do r=1,vl
        q(i,a) = q(i,a) + phi(r,i) * vx(r) * chi(r,a)
      enddo
    enddo
    enddo
#else
    ! BLAS VERSION WITH TEMP MATRIX:
    do i=1,ni
      phiv(:vl,i) = phi(:vl,i) * vx(:vl)
    enddo

!   q = q + matmul(transpose(phiv),chi(:vl,o+1:o+na))

    ! NOTE: beta == 1 -- add the contribution:
    call dgemm( 't', 'n', ni, na, vl,    1.0_r8_kind &
              , phiv(1,1) , size(phiv,1)             &
              , chi(1,o+1), size(chi,1), 1.0_r8_kind &
              , q(1,o+1)  , size(q,1)                &
              )
#endif
  end subroutine gridInt2

  subroutine symWtsNGra(vl,wx,wq)
    !
    ! Transform cartesian gradients of the weights
    ! from a legacy storage ``rx'' to gradients wrt symmetric modes
    !
    !               wq(:vl,i_grad)
    !
    use unique_atom_module, only: unique_atom_grad_info
    use datatype, only: arrmat3
    implicit none
    integer(i4_kind), intent(in) :: vl
    type(arrmat3), intent(in)    :: wx(:)   ! (n_ua)%m(:vl,k,i_ea)
    real(r8_kind), intent(out)   :: wq(:,:) ! (:vl,i_grad)
    ! *** end of interface ***

    integer(i4_kind) :: i_ua,i_ea,k
    integer(i4_kind) :: off,ng
    real(r8_kind)    :: w3(vl,3)
    real(r8_kind), pointer :: rot(:,:,:) ! (ng,3,n_ea)

    wq(:vl,:) = 0.0_r8_kind
    off = 0
    do i_ua=1,size(unique_atom_grad_info)

      rot => unique_atom_grad_info(i_ua)%m(:,:,:)
      ng  = size(rot,1)

      do i_ea=1,size(rot,3)
        if( is_unity(rot(:,:,i_ea)) )then
           ! contribution to tot. sym. grad due to this (ua,ea):
           call add(off,ng,wx(i_ua)%m(:vl,:,i_ea),wq)
        else
           do k=1,ng
             w3(:vl,k) = wx(i_ua)%m(:vl,1,i_ea) * rot(k,1,i_ea) &
                       + wx(i_ua)%m(:vl,2,i_ea) * rot(k,2,i_ea) &
                       + wx(i_ua)%m(:vl,3,i_ea) * rot(k,3,i_ea)
           enddo
           ! contribution to tot. sym. grad due to this (ua,ea):
           call add(off,ng,w3,wq)
        endif
      enddo

      ! advance symm mode offset to the next ua:
      off = off + ng
    enddo
  contains
    subroutine add(off,ng,v3,vx)
      implicit none
      integer(i4_kind), intent(in)    :: off,ng
      real(r8_kind)   , intent(in)    :: v3(:,:) ! (vl,<=3)
      real(r8_kind)   , intent(inout) :: vx(:,:) ! (vl,n_modes)
      ! *** end of interface ***

      integer(i4_kind) :: k

      do k=1,ng
        vx(:,off+k) = vx(:,off+k) + v3(:,k)
      enddo
    end subroutine add
  end subroutine symWtsNGra

  function is_unity(rot) result(yes)
    USE DEBUG, only: NaN
    implicit none
    real(r8_kind), intent(in) :: rot(:,:) ! (3?,3?)
    logical                   :: yes ! result
    ! *** end of interface ***

    real(r8_kind), parameter :: unity(3,3) &
       = reshape( (/ &
       1.0_r8_kind, 0.0_r8_kind, 0.0_r8_kind, &
       0.0_r8_kind, 1.0_r8_kind, 0.0_r8_kind, &
       0.0_r8_kind, 0.0_r8_kind, 1.0_r8_kind  /), &
       (/3,3/))

    yes = .false.

    if( size(rot,1) /= size(rot,2) ) RETURN

    yes = sum( (rot - unity)**2 ) < 1.0e-7_r8_kind
  end function is_unity

 subroutine cartNGra(vl,symm,cart)
   !
   ! transforms gradients 
   !
   !           symm(:vl,m,spin)
   !
   ! from symmetric mode into cartesian gradients
   ! to store them in a legacy structure
   !
   !           cart(ua)%m(:vl,3,n_ea,spin)
   !
   use unique_atom_module, only: unique_atom_grad_info
   use datatype, only: arrmat4
   implicit none
   integer(i4_kind), intent(in)    :: vl
   real(r8_kind)   , intent(in)    :: symm(:,:,:) ! (:vl,n_modes,ispin)
   type(arrmat4)   , intent(inout) :: cart(:) ! cartesian grads array
                                    ! (ua)%m(:vl,3,n_ea,ispin)
   ! *** end of interface ***

   real(r8_kind), pointer :: rot(:,:,:) ! => (ng,3,n_ea)
   integer(i4_kind) :: ng,k
   integer(i4_kind) :: nu,u
   integer(i4_kind) :: ne,e
   integer(i4_kind) :: m
   integer(i4_kind) :: z
   integer(i4_kind) :: spin

   nu = size(unique_atom_grad_info)

   ![[ loop over modes:
   m = 0
   do u=1,nu
   rot => unique_atom_grad_info(u)%m!(ng,1:3,n_ea)
   ng  =  size(rot,1)
   ne  =  size(rot,3)
   cart(u)%m(:vl,:,:,:) = 0.0_r8_kind
   do k=1,ng ! symm modes
   m = m + 1

     do spin=1,size(symm,3) !(( ispin

       do e=1,ne !(( equal atoms
       do z=1,3   ! cartesian x,y,z

         cart(u)%m(:vl,z,e,spin) =            &
         cart(u)%m(:vl,z,e,spin) + rot(k,z,e) &
                                 * symm(:vl,m,spin)
       enddo
       enddo !)) equal atoms
     enddo !)) ispin
   enddo
   enddo !]] loop over modes
   ASSERT(m==size(symm,2))
 end subroutine cartNGra

 subroutine cartNDer(symm,cart)
   !
   ! transforms second derivatives 
   !
   !           symm(m1,m2,spin)
   !
   ! from symmetric modes m1,m2 into cartesian derivatives
   ! to store them in a legacy structure
   !
   !           cart(ua1,ua2)%m(3,n_ea1,3,n_ea2,spin)
   !
   use unique_atom_module, only: unique_atom_grad_info
   use datatype, only: arrmat5
   implicit none
   real(r8_kind)   , intent(in)    :: symm(:,:,:) ! (n_modes,n_modes,ispin)
   type(arrmat5)   , intent(inout) :: cart(:,:) ! cartesian dervs array
                                    ! (u1,u2)%m(3,n_ea1,3,n_ea2,ispin)
   ! *** end of interface ***

   real(r8_kind), pointer :: rot1(:,:,:) ! => (ng,3,n_ea)
   real(r8_kind), pointer :: rot2(:,:,:) ! => (ng,3,n_ea)
   integer(i4_kind) :: ng1,k1,ng2,k2
   integer(i4_kind) :: nu,u1,u2
   integer(i4_kind) :: ne1,e1,ne2,e2
   integer(i4_kind) :: m1,m2
   integer(i4_kind) :: z1,z2
   integer(i4_kind) :: spin

   nu = size(unique_atom_grad_info)


   do u2=1,nu
   do u1=1,nu
     cart(u1,u2)%m = 0.0_r8_kind
   enddo
   enddo

   ![[ loop over modes 2:
   m2 = 0
   do u2=1,nu
   rot2 => unique_atom_grad_info(u2)%m!(ng,1:3,n_ea)
   ng2  =  size(rot2,1)
   ne2  =  size(rot2,3)
   do k2=1,ng2 ! symm modes 2
   m2 = m2 + 1

   ![[ loop over modes 1:
   m1 = 0
   do u1=1,nu
   rot1 => unique_atom_grad_info(u1)%m!(ng,1:3,n_ea)
   ng1  =  size(rot1,1)
   ne1  =  size(rot1,3)
   do k1=1,ng1 ! symm modes 1
   m1 = m1 + 1

     do spin=1,size(symm,3) !(( ispin

       do e2=1,ne2 !(( equal atoms 2
       do z2=1,3   ! cartesian x,y,z

       do e1=1,ne1 !(( equal atoms 1
       do z1=1,3   ! cartesian x,y,z

         cart(u1,u2)%m(z1,e1,z2,e2,spin) =                &
         cart(u1,u2)%m(z1,e1,z2,e2,spin) + rot1(k1,z1,e1) &
                                         * rot2(k2,z2,e2) &
                                         * symm(m1,m2,spin)
       enddo
       enddo !)) equal atoms 1

       enddo
       enddo !)) equal atoms 2
     enddo !)) ispin
   enddo
   enddo !]] loop over modes 1
   ASSERT(m1==size(symm,1))

   enddo
   enddo !]] loop over modes 2
   ASSERT(m2==size(symm,2))
   ASSERT(m2==m1)
 end subroutine cartNDer

  !--------------- End of module ----------------------------------
end module cpks_grid_utils

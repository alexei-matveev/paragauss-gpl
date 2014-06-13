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
module cpks_utils
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
  use type_module, only: &
      IK=>i4_kind,       &
      RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: cpks_h1_store
  public :: cpks_bcast_hr1
  public :: cpks_free_hr1

  public :: cpks_calc_p1w1

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine cpks_calc_p1w1(iorder,p0,w0,p1,w1)
    !
    !  Purpose: Computes (Energy) Density Matirx Gradients
    !
    !------------ Modules used ------------------- ---------------
    use datatype, only: arrmat2, arrmat3
    use cpksdervs_matrices, only: cpks ! %h1, %s1, %HBH
    use eigen_data_module, only: eigval, eigvec
    use matrix_matmult, only: mm=>matmult
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK)  , intent(in)              :: iorder  ! 0 or 1
    type(arrmat3), intent(inout), optional :: p0(:)   ! (n_irr)%m(dim_irr,dim_irr,n_spn)
    type(arrmat2), intent(inout), optional :: w0(:)   ! (n_irr)%m(dim_irr,dim_irr)
    type(arrmat3), intent(inout), optional :: p1(:,:) ! (n_irr,n_gra)%m(dim_irr,dim_irr,n_spn)
    type(arrmat2), intent(inout), optional :: w1(:,:) ! (n_irr,n_gra)%m(dim_irr,dim_irr)
    ! *** end of interface ***

    integer(IK)              :: n,ni,na,ns,nir,ng
    integer(IK)              :: i,j,a,p,q
    integer(IK)              :: ir,s,x
    real(RK), allocatable    :: D(:,:) ! (n ,n )
    real(RK), allocatable    :: o(:,:) ! (ni,ni)
    real(RK), allocatable    :: v(:,:) ! (ni,na)

    FPP_TIMER_DECL(tot)

    FPP_TIMER_START(tot)
    DPRINT  'cpks_calc_p1w1() entered'

    select case(iorder)
    case(0)
      ASSERT(present(p0))
      ASSERT(present(w0))
    case(1)
      ASSERT(present(p1))
      ASSERT(present(w1))
    case default
      ABORT('no such case')
    end select

    if( .not. allocated(eigval) )then
      ABORT('eigval not alloc 1')
    endif
    if( .not. allocated(eigval(1)%m) )then
      ABORT('eigval not alloc 2')
    endif
    if( .not. allocated(eigvec) )then
      ABORT('eigvec not alloc 1')
    endif
    if( .not. allocated(eigvec(1)%m) )then
      ABORT('eigvec not alloc 2')
    endif
    if( .not. allocated(cpks) )then
      ABORT('cpks not alloc')
    endif
    if( .not. associated(cpks(1,1,1)%h1) )then
      ABORT('%h1 not alloc')
    endif
    if( .not. associated(cpks(1,1,1)%HBH) )then
      ABORT('%HBH not alloc')
    endif

    ng  = size(cpks,1)
    nir = size(cpks,2)
    ns  = size(cpks,3)

    ASSERT(ns==1)
    do s=1,ns ! spin
    do ir=1,nir
      n =size(eigvec(ir)%m,1)
      ni=size(cpks(1,ir,s)%HBH,1)
      na=size(cpks(1,ir,s)%HBH,2)
      ASSERT(ni+na==n)

      allocate(D(n ,n ))
      allocate(o(ni,ni))
      allocate(v(ni,na))

      ! (occupied) Eigenvalues:

      ![[=== Density Matrix ==========================
      if( present(p0) )then
         D = 0.0_rk
         do i=1,ni
          do p=1,n
          do q=1,n
             D(p,q) = D(p,q)              &
                    + occ(ir,i,s)         &
                    * eigvec(ir)%m(p,i,s) &
                    * eigvec(ir)%m(q,i,s)
          enddo
          enddo
         enddo
         ASSERT(allocated(p0(ir)%m))
         p0(ir)%m(:,:,s) = D(:,:)
      endif
      !]]=============================================

      ![[=== Energy-Density Matrix ===================
      if( present(w0) )then
         D = 0.0_rk
         do i=1,ni
          do p=1,n
          do q=1,n
             D(p,q) = D(p,q)              &
                    + occ(ir,i,s)         &
                    * eigval(ir)%m(i,s)   &
                    * eigvec(ir)%m(p,i,s) &
                    * eigvec(ir)%m(q,i,s)
          enddo
          enddo
         enddo
         ASSERT(allocated(w0(ir)%m))
         w0(ir)%m(:,:) = D(:,:)
      endif
      !]]=============================================

      do x=1,ng

        if( present(p1) )then
        ![[=== Gradient of Density Matrix ===================
           D = 0.0_rk

           ! P1, Occupied x Virtual:  ni * X1(i,a)
           do a=1,na
           do i=1,ni
              v(i,a) = cpks(x,ir,s)%HBH(i,a) &
                     * occ(ir,i,s)
           enddo
           enddo

#if _SLOW_REF
           do a=1,na
           do i=1,ni
              do q=1,n
              do p=1,n
                 D(p,q) = D(p,q)                                         &
                        + v(i,a)                                         &
                        * ( eigvec(ir)%m(p,i,s) * eigvec(ir)%m(q,ni+a,s) )
              enddo
              enddo
           enddo
           enddo
#else
           call mm2(ni,eigvec(ir)%m(:,:,s),v,D)
#endif
           D = D + transpose(D)

           ! P1, Occupied x Occupied: ( ni + nj ) * ( - S1(i,j) / 2 )
           do i=1,ni
           do j=1,ni
              o(i,j) = cpks(x,ir,s)%s1(i,j) &
                     * ( - ( occ(ir,i,s) + occ(ir,j,s) ) / 2 )
           enddo
           enddo

#if _SLOW_REF
           do j=1,ni
           do i=1,ni
              do q=1,n
              do p=1,n
                 D(p,q) = D(p,q)              &
                        + o(i,j)              &
                        * eigvec(ir)%m(p,i,s) &
                        * eigvec(ir)%m(q,j,s)
              enddo
              enddo
           enddo
           enddo
#else
           call mm2(0,eigvec(ir)%m(:,:,s),o,D)
#endif

           ASSERT(allocated(p1(ir,x)%m))
           p1(ir,x)%m(:,:,s) = D(:,:)
        endif
        !]]==================================================

        ![[=== Gradient of Energy-Density Matrix ============
        if( present(w1) )then
           D = 0.0_rk

           ! W1, Occupied x Virtual:  ni * e(i) * X1(i,a)
           do a=1,na
           do i=1,ni
              v(i,a) = cpks(x,ir,s)%HBH(i,a) &
                     * ( eigval(ir)%m(i,s) * occ(ir,i,s) )
           enddo
           enddo

#if _SLOW_REF
           do a=1,na
           do i=1,ni
              do q=1,n
              do p=1,n
                 D(p,q) = D(p,q)                                         &
                        + v(i,a)                                         &
                        * ( eigvec(ir)%m(p,i,s) * eigvec(ir)%m(q,ni+a,s) )
              enddo
              enddo
           enddo
           enddo
#else
           call mm2(ni,eigvec(ir)%m(:,:,s),v,D)
#endif
           D = D + transpose(D)

           ! === Two similar terms: ===
           ! W1.1, Occupied x Occupied: F1(i,j) + ( ei + ej ) * ( - S1(i,j) / 2 )
           ! NOTE: diagonal elements of are the derivatives of eigenvalues
           !                  e1(i) = f1(i,i)
           !
           ! W1.2, Occupied x Occupied: ( ni * e(i) + nj * e(j) ) * ( - S1(i,j) / 2 )

           do j=1,ni
           do i=1,ni
              o(i,j) = ( cpks(x,ir,s)%h1(i,j)                        & ! W1.1
                       - cpks(x,ir,s)%s1(i,j) / 2                    &
                       * (  eigval(ir)%m(i,s) + eigval(ir)%m(j,s)  ) &
                       ) * occ(ir,i,s)                               &
                       - cpks(x,ir,s)%s1(i,j) / 2                    & ! W1.2
                       * ( eigval(ir)%m(i,s) * occ(ir,i,s)           &
                         + eigval(ir)%m(j,s) * occ(ir,j,s) )
           enddo
           enddo

#if _SLOW_REF
           ! W1, Diag of Occupied x Occupied: ni * e1(i)
           do j=1,ni
           do i=1,ni
              do q=1,n
              do p=1,n
                 D(p,q) = D(p,q)               &
                        + o(i,j)               &
                        * eigvec(ir)%m(p,i,s)  &
                        * eigvec(ir)%m(q,j,s)
              enddo
              enddo
           enddo
           enddo
#else
           call mm2(0,eigvec(ir)%m(:,:,s),o,D)
#endif

           ASSERT(allocated(w1(ir,x)%m))
           w1(ir,x)%m(:,:) = D(:,:)
        endif
        !]]==================================================

      enddo ! x

      deallocate(D)
      deallocate(o,v)
    enddo   ! irrep
    enddo   ! spin

    FPP_TIMER_STOP(tot)
    DPRINT  'cpks_calc_p1w1() exit'
#ifdef FPP_TIMERS
    print *,'cpks_calc_p1w1() time req=',FPP_TIMER_VALUE(tot)
#endif
  contains

     function occ(ir,i,s) result(n)
       use occupied_levels_module, only: occ_num_occ
       implicit none
       integer(IK), intent(in) :: ir,i,s
       real(RK)                :: n
       ! *** end of interface ***

       n = occ_num_occ(ir)%m(i,s)
      !n = 1.0_rk
     end function occ

#ifndef _SLOW_REF
     subroutine mm2( l, e, v, d)
       ! increments d(p,q) by
       !
       !   d(p,q) += SUM(i,a) e(p,i) * v(i,a) * e(q,a)
       !
       ! in two matrix multiplications
       ! DONT FORGET TO CLEAR d(p,q)!
       implicit none
       integer(IK), intent(in)    :: l ! 0 or ni
       real(RK)   , intent(in)    :: e(:,:)
       real(RK)   , intent(in)    :: v(:,l+1:)
       real(RK)   , intent(inout) :: d(:,:)
       ! *** end of interface ***

       integer(IK) :: n,ni,na
#ifdef FPP_NOBLAS
       integer(IK) :: p,q,i,a
#endif

       real(RK) :: w(size(e,1),size(v,1)) ! w(n,ni)

       n  = size(e,1)
       ni = size(v,1)
       na = size(v,2)

#ifdef FPP_NOBLAS
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning
       do i=1,ni
         do p=1,n
           w(p,i) = 0.0_rk
           do a=l+1,l+na
             w(p,i) = w(p,i) + v(i,a) * e(p,a)
           enddo
         enddo
       enddo
       do q=1,n
         do p=1,n
           do i=1,ni
             d(p,q) = d(p,q) + w(p,i) * e(q,i)
           enddo
         enddo
       enddo
#else
      if( ni==0 .or. na==0 ) RETURN
      call dgemm( 'n', 't', n,ni,na,   1.0D0 &
                , e(1,l+1), size(e,1)        &
                , v(1,l+1), size(v,1), 0.0D0 &
                , w(1,1)  , size(w,1)        &
                )
      ! NOTE: beta == 1 -- add the contribution:
      call dgemm( 'n', 't', n,n,ni,    1.0D0 &
                , e(1,1)  , size(e,1)        &
                , w(1,1)  , size(w,1), 1.0D0 &
                , d(1,1)  , size(d,1)        &
                )
#endif
     end subroutine mm2
#endif
  end subroutine cpks_calc_p1w1

  subroutine cpks_h1_store(irr,igr,h,iop)
    !
    !  Purpose: adds another contribution to the first
    !           derivative of hamiltonian H1.
    !           H1(i,j) and H(i,a) belong to CPKS input.
    !
    !           So far is called only in relativistic case
    !           to add the 1st derivative of the
    !           *relativistically transformed* Hamiltonian
    !                        Hrel=Trel+Vrel
    !           So far is called once for each irrep and gradient
    !
    !------------ Modules used ------------------- ---------------
    use cpksdervs_matrices, only: cpks
    use eigen_data_module, only: eigvec
    use matrix_module, only: sim
    USE_MEMLOG
    USE DEBUG, only: show
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: irr, igr
    real(RK)   , intent(in) :: h(:,:)   ! (dim_irr,dim_irr)
    integer(IK), intent(in) :: iop      ! 0,+/-1
    ! *** end of interface ***

    integer(IK) :: n,ni,na,ns
    integer(IK) :: s
    real(RK)    :: h1(size(h,1),size(h,2)) ! (n,n)
    integer(IK) :: memstat

    DPRINT  'cpks_h1_store(',irr,igr,'iop=',iop,')'

    if( .not. allocated(eigvec) )then
      ABORT('eigvec not alloc 1')
    endif
    if( .not. allocated(eigvec(irr)%m) )then
      ABORT('eigvec not alloc 2')
    endif
    if( .not. allocated(cpks) )then
      ABORT('cpks not alloc')
    endif
    if( .not. associated(cpks(igr,irr,1)%h1) )then
      ABORT('%h1 not alloc')
    endif
    if( .not. associated(cpks(igr,irr,1)%h1ai) )then
      ABORT('%h1ai not alloc')
    endif

    n = size(h,1)
    ASSERT(n==size(h,2))
    ASSERT(n==size(eigvec(irr)%m,1))
    ASSERT(n==size(eigvec(irr)%m,2))

    ![[=== Store the rel. ham. grad. till after CPKS here
    !      (in the symmetrized basis)
    if( .not. associated(cpks(igr,irr,1)%hr1) )then
      ! n == dimension of irrep basis:
      allocate( cpks(igr,irr,1)%hr1(n,n), stat=memstat )
      ASSERT(memstat==0)
      MEMLOG(size(cpks(igr,irr,1)%hr1))
    else
      ASSERT(n==size(cpks(igr,irr,1)%hr1,1))
      ASSERT(n==size(cpks(igr,irr,1)%hr1,2))
      WARN('%hr1 already alloced?')
    endif

    DPRINT  'cpks_h1_store(',igr,irr,'): saving rel ham'
    cpks(igr,irr,1)%hr1 = h
    !]]============================================

    ![[=== Add the rel. ham. to the Fock mat. grad.
    !      (in the eigenvector basis)
    ns = size(eigvec(irr)%m,3)
    do s=1,ns ! spin
      ni=size(cpks(igr,irr,s)%h1  ,1)
      na=size(cpks(igr,irr,s)%h1ai,2)
      ASSERT(ni==size(cpks(igr,irr,s)%h1,2))
      ASSERT(ni==size(cpks(igr,irr,s)%h1ai,1))
      if ( ni + na /= n ) then
         WARN('case with holes!')
      endif

      DPRINT  'RGR: n_occ=',ni,' n_vir=',na,' n=',n
      DPRINT  'RGR: shape(h1  )=',shape(cpks(igr,irr,s)%h1)
      DPRINT  'RGR: shape(h1ai)=',shape(cpks(igr,irr,s)%h1ai)
      DPRINT  'RGR: h1  <-',1,':',ni,'x',1,':',ni
      DPRINT  'RGR: h1ai<-',1,':',ni,'x',1+n-na,':',n

      ! FIXME: vir x vir block is actually unneded:
      ! sim(h,U) = U^T * h * U
      h1 = sim( h, eigvec(irr)%m(:,:,s) )

      ! FIXME: if occupied are not the the first ones?
      !        and  virtual are not the the last ones?
      ! increment occ x occ and occ x vir gradients:
      select case (iop)
      case (  0 )
        cpks(igr,irr,s)%h1   =   h1(:ni,:ni)
        cpks(igr,irr,s)%h1ai =   h1(:ni,1+n-na:)
      case ( +1 )
        cpks(igr,irr,s)%h1   = + h1(:ni,:ni)     + cpks(igr,irr,s)%h1
        cpks(igr,irr,s)%h1ai = + h1(:ni,1+n-na:) + cpks(igr,irr,s)%h1ai
      case ( -1 )
        cpks(igr,irr,s)%h1   = - h1(:ni,:ni)     + cpks(igr,irr,s)%h1
        cpks(igr,irr,s)%h1ai = - h1(:ni,1+n-na:) + cpks(igr,irr,s)%h1ai
      case default
        ABORT('no such case')
      end select
    enddo ! spin
    !]]============================================
  end subroutine cpks_h1_store

  subroutine cpks_bcast_hr1()
    !
    !  Purpose: Copies cpks%hr1 to others
    !
    !------------ Modules used ------------------- ---------------
    use cpksdervs_matrices, only: cpks
    use comm
    implicit none
    !------------ Declaration of formal parameters ---------------
    ! *** end of interface ***

    integer(IK) :: n, irr, igr
    integer(IK) :: root
    integer(IK) :: memstat

    if( .not. comm_parallel() ) RETURN

    do irr = 1, size(cpks, 2)
    do igr = 1, size(cpks, 1)

      ! FIXME: is there a rule by which this data is distributed?
      if( associated(cpks(igr, irr, 1)%hr1) )then
        ! sender:
        root = comm_rank()
        n = size(cpks(igr, irr, 1)%hr1,1)
      else
        ! receivers:
        root = 0 ! zero for reduction by addition
      endif

      ! FIXME: a better way to negotiate the root role?
      call comm_reduce(root) ! on master
      call comm_bcast(root) ! from master

      ! tell receivers the dimension:
      call comm_bcast(n, root=root)

      if( .not. associated(cpks(igr, irr, 1)%hr1) )then
        ! receiver:
        allocate(cpks(igr, irr, 1)%hr1(n, n), stat=memstat)
        ASSERT(memstat==0)
      endif

      ! now actual data:
      call comm_bcast(cpks(igr, irr, 1)%hr1, root=root)
    enddo
    enddo
  end subroutine cpks_bcast_hr1

  subroutine cpks_free_hr1()
    !
    !  Purpose: Copies cpks%hr1 to others
    !
    !------------ Modules used ------------------- ---------------
    use cpksdervs_matrices, only: cpks
    implicit none
    !------------ Declaration of formal parameters ---------------
    ! *** end of interface ***

    integer(IK) :: irr,igr
    integer(IK) :: memstat

    ASSERT(allocated(cpks))
    do irr=1,size(cpks,2)
    do igr=1,size(cpks,1)
      ASSERT(associated(cpks(igr,irr,1)%hr1))
      deallocate( cpks(igr,irr,1)%hr1, STAT=memstat )
      ASSERT(memstat==0)
    enddo
    enddo
  end subroutine cpks_free_hr1

  !--------------- End of module ----------------------------------
end module cpks_utils

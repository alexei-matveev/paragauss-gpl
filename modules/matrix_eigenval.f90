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
module matrix_eigenval
  !-------------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  ! Copyright (c) Raghunathan Ramakrishnan
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

#include "def.h"

  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind,&
       & CK=>c16_kind ! type specification parameters
  use matrix_types
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------
  interface eigs
     !
     ! Eigensolver (Simple Eigenvalue Problem)
     !
     module procedure eigs_r_rd_r
     module procedure eigs_c_rd_c

     module procedure eigs_dsyev
     module procedure eigs_zheev
  end interface

  interface geigs
     !
     ! Eigensolver (General Eigenvalue Problem)
     !
     module procedure geigs_r_r_rd_r
     module procedure geigs_c_c_rd_c

     module procedure geigs_dsygv
     module procedure geigs_zhegv
  end interface

  interface svd
     !
     ! SVD decomposition
     !
     module procedure svd_c_c_rd_c

     module procedure svd_plain
     module procedure zgesvd90
  end interface

  !------------ public functions and subroutines ---------------------
  public :: eigs, geigs, svd
  public :: jacobi

  public :: dgesvd90!(A, U, s, VT)
  public :: dsygv90!(A, B, e)

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine geigs_dsygv(H,S,eigval,EVEC)
    implicit none
    real(RK),dimension(:,:),intent(in)  :: H,S
    real(RK),dimension(:),  intent(out) :: eigval
    real(RK),dimension(:,:),intent(out) :: EVEC
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: N,i
    real(RK),allocatable    :: A(:,:),B(:,:)

    N = size(eigval)

    if ( N .eq. 0 ) return ! for the case of cases

    ASSERT(N.eq.size(H,1))
    ASSERT(N.eq.size(H,2))
    ASSERT(N.eq.size(S,1))
    ASSERT(N.eq.size(S,2))
    ASSERT(N.eq.size(EVEC,1))
    ASSERT(N.eq.size(EVEC,2))

    allocate(A(N,N),B(N,N),STAT=memstat)
    ASSERT(memstat==0)

    A = H
    B = S

    call dsygv90(A,B,eigval)

    do i=1,N
       EVEC(:,i) =  A(:,i)
    enddo

    deallocate(A,B,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine geigs_dsygv

  subroutine geigs_zhegv(H_re,H_im,S_re,S_im,eigval,EVEC_re,EVEC_im)
    implicit none
    real(RK),dimension(:,:),intent(in)  :: H_re,H_im,S_re,S_im
    real(RK),dimension(:),  intent(out) :: eigval
    real(RK),dimension(:,:),intent(out) :: EVEC_re,EVEC_im
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: N,i!!$,j
    complex(CK),allocatable :: A(:,:),B(:,:) !!$,s(:)

    N = size(eigval)

    if ( N .eq. 0 ) return ! for the case of cases

    ASSERT(N.eq.size(H_re,1))
    ASSERT(N.eq.size(H_re,2))
    ASSERT(N.eq.size(H_im,1))
    ASSERT(N.eq.size(H_im,2))
    ASSERT(N.eq.size(S_re,1))
    ASSERT(N.eq.size(S_re,2))
    ASSERT(N.eq.size(S_im,1))
    ASSERT(N.eq.size(S_im,2))
    ASSERT(N.eq.size(EVEC_re,1))
    ASSERT(N.eq.size(EVEC_re,2))
    ASSERT(N.eq.size(EVEC_im,1))
    ASSERT(N.eq.size(EVEC_im,2))

    allocate(A(N,N),B(N,N),STAT=memstat)
    ASSERT(memstat==0)

    A = cmplx(H_re,H_im,RK)
    B = cmplx(S_re,S_im,RK)

    call zhegv90(A,B,eigval)

    do i=1,N
       EVEC_re(:,i) =  real(A(:,i))
       EVEC_im(:,i) = aimag(A(:,i))
    enddo

    deallocate(A,B,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine geigs_zhegv

  subroutine eigs_dsyev(H,eigval,EVEC)
    implicit none
    real(RK),dimension(:,:),intent(in)  :: H
    real(RK),dimension(:),  intent(out) :: eigval
    real(RK),dimension(:,:),intent(out) :: EVEC
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: N,i
    real(RK),allocatable    :: A(:,:)

    N = size(eigval)

    if ( N .eq. 0 ) return ! for the case of cases

    ASSERT(N.eq.size(H,1))
    ASSERT(N.eq.size(H,2))
    ASSERT(N.eq.size(EVEC,1))
    ASSERT(N.eq.size(EVEC,2))

    allocate(A(N,N),STAT=memstat)
    ASSERT(memstat.eq.0)

    A = H

    call dsyev90(A,eigval)

    do i=1,N
       EVEC(:,i) =  A(:,i)
    enddo

    deallocate(A,STAT=memstat)
    ASSERT(memstat.eq.0)
  end subroutine eigs_dsyev

  subroutine eigs_zheev(H_re,H_im,eigval,EVEC_re,EVEC_im)
    implicit none
    real(RK),dimension(:,:),intent(in)  :: H_re,H_im
    real(RK),dimension(:),  intent(out) :: eigval
    real(RK),dimension(:,:),intent(out) :: EVEC_re,EVEC_im
    ! *** end of interface ***

    real(RK),parameter      :: one = 1.0_rk 
    integer(IK)             :: memstat
    integer(IK)             :: N,i!!$,j
    complex(CK),allocatable :: A(:,:)!!$,B(:,:),s(:)

    N = size(eigval)

    if ( N .eq. 0 ) return ! for the case of cases

    ASSERT(N.eq.size(H_re,1))
    ASSERT(N.eq.size(H_re,2))
    ASSERT(N.eq.size(H_im,1))
    ASSERT(N.eq.size(H_im,2))
    ASSERT(N.eq.size(EVEC_re,1))
    ASSERT(N.eq.size(EVEC_re,2))
    ASSERT(N.eq.size(EVEC_im,1))
    ASSERT(N.eq.size(EVEC_im,2))

    allocate(A(N,N),STAT=memstat)
    ASSERT(memstat.eq.0)

    A = cmplx(H_re,H_im,RK)

    call zheev90(A,eigval)

    do i=1,N
       EVEC_re(:,i) =  real(A(:,i))
       EVEC_im(:,i) = aimag(A(:,i))
    enddo

    deallocate(A,STAT=memstat)
    ASSERT(memstat.eq.0)
  end subroutine eigs_zheev

  subroutine dsygv90(A, B, EIGVAL)
    !
    ! General Eigenvalue Problem Solver for
    ! real symmetric matrices (LAPACK/DSYGV)
    !
    ! A * Psi = E * B * Psi
    !
    ! A(inout)   : in : A, upper triangle
    !              out: eigenvectors
    ! B(inout)   : in : B, upper triangle(?)
    !              out: Cholesky factorization matrix
    ! EIGVAL(out): eigenvalues 
    !
    use f77_lapack, only: DSYGV, ILAENV
!   use debug, only: show
    implicit none
    real(RK), intent(inout) :: A(:,:), B(:,:) ! LAPACK overwrites them ...
    real(RK), intent(out) :: EIGVAL(:)
    ! *** end of interface ***

    integer(IK)             :: N, NB
    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    real(RK),allocatable    :: WORK(:)


    N = size(EIGVAL)
    DPRINT 'mev/dsygv90: entered, N=',N

    if(N.eq.0)return
    ASSERT(N.eq.size(A,1))
    ASSERT(N.eq.size(A,2))
    ASSERT(N.eq.size(B,1))
    ASSERT(N.eq.size(B,2))

    NB = ILAENV(1,'DSYTRD','U',N,-1,-1,-1)
    DPRINT 'mev/dsygv90: ILAENV returned ',NB

    LWORK = (NB+2)*N
    ASSERT(LWORK>=max(1,3*N-1))
    allocate(WORK(LWORK),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/dsygv90: call DSYGV(',N,')'
    call DSYGV(1, 'V', 'U', N, A, N, B, N, EIGVAL, WORK, LWORK, INFO)
    DPRINT 'mev/dsygv90: done DSYGV; INFO=',INFO,' LWORK mismatch by ',LWORK - NINT(WORK(1))
    if(INFO.ne.0)then
      print *,'mev/dsygv90: DSYGV(1,V,U,',N,'A,',N,'B,',N,'EIGVAL,WORK,LWORK,INFO) failed!'
      print *,'mev/dsygv90: DSYGV returned INFO=',INFO,' LWORK mismatch by ',LWORK - NINT(WORK(1))
!     call show('A',A)
!     call show('B',B)
      ABORT('DSYGV returned INFO/=0, see tty')
    endif

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/dsygv90: exit'
  end subroutine dsygv90

  subroutine zhegv90(A,B,EIGVAL)
    !
    ! General Eigenvalue Problem Solver for
    ! complex hermitean matrices (LAPACK/ZHEGV)
    !
    ! A * Psi = E * B * Psi
    !
    ! A(inout)   : in : A, upper triangle
    !              out: eigenvectors
    ! B(inout)   : in : B, upper triangle(?)
    !              out: Cholesky factorization matrix
    ! EIGVAL(out): eigenvalues 
    !
    use f77_lapack, only: ZHEGV, ILAENV
    implicit none
    complex(RK),intent(inout) :: A(:,:),B(:,:)
    real(RK),intent(out)      :: EIGVAL(:)
    ! *** end of interface ***

    integer(IK)             :: N, NB
    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    complex(RK),allocatable :: WORK(:)
    real(RK),allocatable    :: RWORK(:)


    N = size(EIGVAL)
    DPRINT 'mev/zhegv90: entered, N=',N

    if(N.eq.0)return
    ASSERT(N.eq.size(A,1))
    ASSERT(N.eq.size(A,2))
    ASSERT(N.eq.size(B,1))
    ASSERT(N.eq.size(B,2))

    NB = ILAENV(1,'ZHETRD','U',N,-1,-1,-1)
    DPRINT 'mev/zhegv90: ILAENV returned ',NB

    LWORK = (NB+1)*N
    ASSERT(LWORK>=max(1,2*N-1))
    allocate(WORK(LWORK),RWORK(3*N-2),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/zhegv90: call ZHEGV(',N,')'
    call ZHEGV(1,'V','U',N,A,N,B,N,EIGVAL,WORK,LWORK,RWORK,INFO)
    DPRINT 'mev/zhegv90: done ZHEGV; INFO=',INFO,' LWORK mismatch by ',LWORK - NINT(REAL(WORK(1)))
    ASSERT(INFO.eq.0)

    deallocate(WORK,RWORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/zhegv90: exit'
  end subroutine zhegv90

  subroutine dsyev90(A,EIGVAL)
    !
    ! Symmetric Eigenvalue Problem Solver for
    ! real symmetric matrices (LAPACK/DSYEV)
    !
    ! A * Psi = E * Psi
    !
    ! A(inout)   : in : A, upper triangle
    !              out: eigenvectors
    ! EIGVAL(out): eigenvalues 
    !
    use f77_lapack, only: DSYEV !,ILAENV
    implicit none
    real(RK),intent(inout)    :: A(:,:)
    real(RK),intent(out)      :: EIGVAL(:)
    ! *** end of interface ***

    integer(IK)             :: N
    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    real(RK),allocatable    :: WORK(:)

    N = size(EIGVAL)
    DPRINT 'mev/dsyev90: entered, N=',N

    if(N.eq.0)return
    ASSERT(N.eq.size(A,1))
    ASSERT(N.eq.size(A,2))

    allocate(WORK(1),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/dsyev90: inquire DSYEV(',N,')'
    call DSYEV('V','U',N,A,N,EIGVAL,WORK,-1,INFO)
    ASSERT(INFO.eq.0)

    LWORK = NINT(WORK(1))
    DPRINT 'mev/dsyev90: optimal LWORK=',LWORK

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    allocate(WORK(LWORK),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/dsyev90: call DSYEV(',N,')'
    call DSYEV('V','U',N,A,N,EIGVAL,WORK,LWORK,INFO)
    DPRINT 'mev/dsyev90: done DSYEV, INFO=',INFO
    ASSERT(INFO.eq.0)

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/dsyev90: exit'
  end subroutine dsyev90

  subroutine zheev90(A,EIGVAL)
    !
    ! Hermitean Eigenvalue Problem Solver for
    ! complex hermitean matrices (LAPACK/ZHEEV)
    !
    ! A * Psi = E * Psi
    !
    ! A(inout)   : in : A, upper triangle
    !              out: eigenvectors
    ! EIGVAL(out): eigenvalues 
    !
    use f77_lapack, only: ZHEEV !,ILAENV
    implicit none
    complex(RK),intent(inout)    :: A(:,:)
    real(RK),intent(out)         :: EIGVAL(:)
    ! *** end of interface ***

    integer(IK)             :: N
    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    complex(RK),allocatable :: WORK(:)
    real(RK),allocatable    :: RWORK(:)


    N = size(EIGVAL)
    DPRINT 'mev/zheev90: entered, N=',N

    if(N.eq.0)return
    ASSERT(N.eq.size(A,1))
    ASSERT(N.eq.size(A,2))

    allocate(WORK(1),RWORK(3*N-2),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/zheev90: inquire ZHEEV(',N,')'
    call ZHEEV('V','U',N,A,N,EIGVAL,WORK,-1,RWORK,INFO)
    ASSERT(INFO.eq.0)

    LWORK = NINT(REAL(WORK(1)))
    DPRINT 'mev/zheev90: optimal LWORK=',LWORK

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    allocate(WORK(LWORK),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/zheev90: call ZHEEV(',N,')'
    call ZHEEV('V','U',N,A,N,EIGVAL,WORK,LWORK,RWORK,INFO)
    DPRINT 'mev/zheev90: done ZHEEV, INFO=',INFO
    ASSERT(INFO.eq.0)

    deallocate(WORK,RWORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/zheev90: exit'
  end subroutine zheev90

  subroutine geigs_r_r_rd_r(H, S, e, V)
    use matrix_types
    use matrix_check, only: square, conform
    use matrix_methods, only: alloc
    implicit none
    type(rmatrix), intent(in) :: H, S
    type(rdmatrix),intent(out) :: e
    type(rmatrix), intent(out) :: V
    ! *** end of interface ***

    if(verbose)print *, 'mm/geigs_r_r_rd_r: entered'

    ASSERT(square(H))
    ASSERT(conform(H,S))

    call alloc(size(H%m,1), size(H%m,2), V)
    call alloc(size(H%m,1), e)

    call geigs(H%m, S%m, e%d, V%m)
  end subroutine geigs_r_r_rd_r

  subroutine geigs_c_c_rd_c(H, S, e, V)
    use matrix_types
    use matrix_check, only: square, conform, mconform
    use matrix_methods, only: alloc
    implicit none
    type(cmatrix), intent(in)  :: H, S
    type(rdmatrix),intent(out) :: e
    type(cmatrix), intent(out) :: V
    ! *** end of interface ***

    if(verbose)print *, 'mm/geigs_c_c_rd_c: entered'

    ASSERT(square(H))
    ASSERT(conform(H,S))

    call alloc(size(H%re,1), size(H%re,2), V)
    call alloc(size(H%re,1), e)

    call geigs(H%re, H%im, S%re, S%im, e%d, V%re, V%im)
  end subroutine geigs_c_c_rd_c

  subroutine eigs_r_rd_r(H, e, V)
    use matrix_types
    use matrix_check, only: square
    use matrix_methods, only: alloc
    implicit none
    type(rmatrix), intent(in) :: H
    type(rdmatrix), intent(out) :: e
    type(rmatrix), intent(out) :: V
    ! *** end of interface ***

    if(verbose)print *, 'mm/eigs_r_rd_r: entered'

    ASSERT(square(H))

    call alloc(size(H%m,1), size(H%m,2), V)
    call alloc(size(H%m,1), e)

    call eigs(H%m, e%d, V%m)
  end subroutine eigs_r_rd_r

  subroutine eigs_c_rd_c(H, e, V)
    use matrix_types, only: cmatrix, rdmatrix
    use matrix_check, only: square, conform, mconform
    use matrix_methods, only: alloc
    implicit none
    type(cmatrix), intent(in) :: H
    type(rdmatrix), intent(out) :: e
    type(cmatrix), intent(out) :: V
    ! *** end of interface ***

    if(verbose)print *, 'mm/eigs_c_rd_c: entered'

    ASSERT(square(H))

    call alloc(size(H%re,1), size(H%re,2), V)
    call alloc(size(H%re,1), e)

    call eigs(H%re, H%im, e%d, V%re, V%im)
  end subroutine eigs_c_rd_c

  subroutine svd_c_c_rd_c(A, U, s, VT)
    !
    ! SVD (singular value decomposition)
    ! of complex MxN matrices (LAPACK/ZGESVD)
    !
    ! A = U * s * V^T
    !
    implicit none
    type(cmatrix), intent(in)   :: A
    type(cmatrix), intent(out)  :: U, VT
    type(rdmatrix), intent(inout) :: s ! FIXME: assumed to be allocated
    ! *** end of interface ***
    
    call svd(A%re, A%im, U%re, U%im, s%d, VT%re, VT%im)
  end subroutine svd_c_c_rd_c

  subroutine svd_plain(A_re,A_im,U_re,U_im,s,VT_re,VT_im)
    !
    ! SVD (singular value decomposition)
    ! of complex MxN matrices (LAPACK/ZGESVD)
    !
    ! A = U * s * V^T
    !
    implicit none
    real(RK), intent(in), dimension(:,:) ::&
         & A_re, A_im
    real(RK), intent(out), dimension(:,:) ::&
         & U_re,U_im,VT_re,VT_im
    real(RK), intent(out) ::&
         & s(:)
    ! *** end of interface ***
    
    integer(IK) :: M,N,memstat
    complex(RK), allocatable, dimension(:,:) ::&
         & A,U,VT

    M = size(A_re,1)
    N = size(A_re,2)

    DPRINT 'mev/svd_plain: entered, M=',M,' N=',N

    if(N.eq.0.or.M.eq.0) return

    ASSERT(M.eq.size(U_re,1))
    ASSERT(M.eq.size(U_re,2))
    ASSERT(N.eq.size(VT_re,1))
    ASSERT(N.eq.size(VT_re,2))
    ASSERT(M.eq.size(U_im,1))
    ASSERT(M.eq.size(U_im,2))
    ASSERT(N.eq.size(VT_im,1))
    ASSERT(N.eq.size(VT_im,2))
    ASSERT(size(s).eq.min(M,N))

    allocate(A(M,N),U(M,M),VT(N,N),STAT=memstat)
    ASSERT(memstat.eq.0)

    A = cmplx(A_re,A_im,RK)

    call zgesvd90(A,U,s,VT)

    U_re = Real(U)
    U_im = aImag(U)
    VT_re = Real(VT)
    VT_im = aImag(VT)

    deallocate(A,U,VT,STAT=memstat)
    ASSERT(memstat.eq.0)
  end subroutine svd_plain

  subroutine zgesvd90(A,U,s,VT)
    !
    ! SVD (singular value decomposition)
    ! of complex MxN matrices (LAPACK/ZGESVD)
    !
    ! A = U * s * V^H
    !
    use f77_lapack, only: ZGESVD
    implicit none
    complex(RK),intent(in)    :: A(:,:) ! MxN
    complex(RK),intent(out)   :: U(:,:) ! MxM
    complex(RK),intent(out)   :: VT(:,:) ! NxN
    real(RK),intent(out)      :: s(:) ! min(M,N)
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    complex(RK),allocatable :: WORK(:)
    real(RK),allocatable    :: RWORK(:)
    integer(IK)             :: M,N

    M = size(A,1)
    N = size(A,2)

    DPRINT 'mev/zgesvd90: entered, M=',M,' N=',N

    if(N.eq.0.or.M.eq.0) return

    ASSERT(M.eq.size(U,1))
    ASSERT(M.eq.size(U,2))
    ASSERT(N.eq.size(VT,1))
    ASSERT(N.eq.size(VT,2))
    ASSERT(size(s).eq.min(M,N))

    ! inquire optimal WORK size:
    allocate(WORK(1),RWORK(5*min(M,N)),STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/zgesvd90: inquire ZGESVD(',M,',',N,')'
    call ZGESVD('A','A',M,N,A,M,s,U,M,VT,N,WORK,-1,RWORK,INFO)
    ASSERT(INFO.eq.0)

    LWORK = NINT(REAL(WORK(1)))
    DPRINT 'mev/zgesvd90: optimal LWORK=',LWORK

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    allocate(WORK(LWORK),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/zgesvd90: call ZGESVD(',M,',',N,')'
    call ZGESVD('A','A',M,N,A,M,s,U,M,VT,N,WORK,LWORK,RWORK,INFO)

    DPRINT 'mev/zgesvd90: done ZGESVD; INFO=',INFO
    ASSERT(INFO.eq.0)

    deallocate(WORK,RWORK,STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/zgesvd90: exit'
  end subroutine zgesvd90

  subroutine dgesvd90(A,U,s,VT)
    !
    ! SVD (singular value decomposition)
    ! of real MxN matrices (LAPACK/DGESVD)
    !
    ! A = U * s * V^T
    !
    implicit none
    real(RK),intent(in)    :: A(:,:) ! MxN
    real(RK),intent(out)   :: U(:,:) ! MxM
    real(RK),intent(out)   :: VT(:,:) ! NxN
    real(RK),intent(out)   :: s(:) ! min(M,N)
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: LWORK,INFO
    real(RK),allocatable    :: WORK(:)
    integer(IK)             :: M,N

    M = size(A,1)
    N = size(A,2)

    DPRINT 'mev/dgesvd90: entered, M=',M,' N=',N

    if(N.eq.0.or.M.eq.0) return

    ASSERT(M.eq.size(U,1))
    ASSERT(M.eq.size(U,2))
    ASSERT(N.eq.size(VT,1))
    ASSERT(N.eq.size(VT,2))
    ASSERT(size(s).eq.min(M,N))

    ! inquire optimal WORK size:
    allocate(WORK(1),STAT=memstat)
    ASSERT(memstat.eq.0)
    DPRINT 'mev/dgesvd90: inquire DGESVD(',M,',',N,')'
    call DGESVD('S','S',M,N,A,M,s,U,M,VT,N,WORK,-1,INFO)
    ASSERT(INFO.eq.0)

    LWORK = NINT(REAL(WORK(1)))
    DPRINT 'mev/dgesvd90: optimal LWORK=',LWORK

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)
    allocate(WORK(LWORK),STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/dgesvd90: call DGESVD(',M,',',N,')'
    call DGESVD('A','A',M,N,A,M,s,U,M,VT,N,WORK,LWORK,INFO)

    DPRINT 'mev/dgesvd90: done DGESVD; INFO=',INFO
    ASSERT(INFO.eq.0)

    deallocate(WORK,STAT=memstat)
    ASSERT(memstat.eq.0)

    DPRINT 'mev/dgesvd90: exit'
  end subroutine dgesvd90


!!$  subroutine eigvalues_typed(H,h_diag)
!!$    use error_module
!!$    use matrix_types
!!$    use matrix_check, only: square
!!$    implicit none
!!$    type(cmatrix),intent(in)     :: H
!!$    type(rdmatrix),intent(inout) :: h_diag
!!$    ! *** end of interface ***
!!$
!!$    if(verbose)print *, 'mm/eigvalues_typed: entered'
!!$
!!$    if(.not.(square(H).and.(H%n.eq.h_diag%n)))&
!!$         & call error("mm/eigvalues_typed: shape conflict")
!!$
!!$    call eigs(H%n,H%re,H%im,h_diag%d)
!!$  end subroutine eigvalues_typed
!!$
!!$  subroutine eigvalues_c_d(H,h_diag)
!!$    use error_module
!!$    use matrix_check
!!$    implicit none
!!$    type(cmatrix),intent(in) :: H
!!$    real(RK),intent(inout)   :: h_diag(:)
!!$    ! *** end of interface ***
!!$
!!$    if(verbose)print *, 'mm/eigvalues_c_d: entered'
!!$
!!$    if(.not.(square(H).and.(H%n.eq.size(h_diag))))&
!!$         & call error("mm/eigvalues_c_d: shape conflict")
!!$
!!$    call eigs(H%n,H%re,H%im,h_diag)
!!$  end subroutine eigvalues_c_d
!!$ 
!!$  subroutine eigvalues_plain_eispack(n,M_real,M_imag,m_diag)
!!$    !
!!$    ! HX=EX type of problems; EISPACK/CH based
!!$    !
!!$    use error_module
!!$    implicit none
!!$    external ch ! <<< eispack
!!$    integer(IK),intent(in)             :: n
!!$    real(RK),dimension(n,n),intent(in) :: M_real,M_imag
!!$    real(RK),dimension(n),intent(out)  :: m_diag
!!$    ! *** end of interface ***
!!$
!!$    real(RK),dimension(n,n) :: M_re,M_im,U_re,U_im
!!$    real(RK)                :: wrk1(n),wrk2(n),wrk3(2,n)
!!$    integer(IK)             :: ierr
!!$
!!$    if(n.eq.0)return !<<< special case, a BIG question
!!$
!!$    M_re = M_real !<<< to secure intent(in) arguments
!!$    M_im = M_imag
!!$
!!$    call ch(n,n,&
!!$         & M_re,M_im,&
!!$         & m_diag,0,&
!!$         & U_re,U_im,&
!!$         & wrk1,wrk2,wrk3,&
!!$         & ierr)
!!$    ASSERT(ierr.eq.0)
!!$  end subroutine eigvalues_plain_eispack
  
  subroutine jacobi(n,AA,lda,e)
    !Diagonalize AA using Jacobi overwriting a with
    !the eigenvectors and returning in e the eigenvalues
    implicit none
    integer(IK), intent(in) :: n,lda
    real(RK), intent(inout) :: AA(lda,n)
    real(RK), intent(out)   :: e(n)
    !
    real(RK), allocatable :: B(:)
    integer(IK) :: status,ij,i,j,info

    allocate(B(n*(n+1)/2),stat=status)
    ASSERT(status==0)

    !Lower triangle of a into B
    ij = 0
    do i = 1, n
       do j = 1, i
          B(1+ij) = AA(j,i)
          ij = ij + 1
       enddo
    enddo

    !unit matrix into AA
    AA=0.0_RK
    do i=1,n
       AA(i,i)=1.0_RK
    end do

    call sjacobi(lda,n,B,e,AA,info)
    ASSERT(info==0)

    deallocate(B,stat=status)
    ASSERT(status==0)

  contains

    subroutine sjacobi(nm,n,a,w,eivr,ierr)
      !use jacobi method to emulate eispack rs
      !this routine uses a variable threshold jacobi method
      !it gives very good eigenvalues and eigenvectors
      !the routine is much faster than the old hdiag routine written
      !at m.i.t. that uses the jacobi method but not the variable
      !threshold technique that is applied here   

      integer(IK) :: nm,n,ierr
      real(RK) :: a(n*(n+1)/2),w(n),eivr(nm,n)

      real(RK) :: sq2inv, t1, t2, avgf, dstop, d, thrsh, aij, daij, aii
      real(RK) :: ajj,s,ds,c,t,u,temp,atop
      integer(IK) :: i,j,k,ij,iflag,jcol,jcoltr,jcol1,irow,irowtr
      integer(IK) :: indxa, i2, itr, jtr, idiag, isort

      !-----parameters---------------
      sq2inv = 1.0_RK/sqrt(2.d0)
      t1     = 1.0e-12_RK
      t2     = 1.0e-12_RK
      !------------------------------
      ierr   = 0
      avgf   = real(n*(n-1))*0.55_RK
      if(n.lt.1) then
         goto 160
      else if(n.lt.2) then
         eivr(1,1)=1.0_RK
         goto 160
      end if
      do j=1,n
         do i=1,n
            eivr(i,j)=0.0_RK
         end do
         eivr(j,j)=1.0_RK
      end do

      !Find the absolutely largest element of a.
      atop=0.0_RK
      ij = 0
      do i=1,n
         do j=1,i
            ij = ij + 1
            if(atop < abs(a(ij))) atop =abs(a(ij))
         end do
      end do
      if(atop <= 0.0_RK) then
         ierr = 1
         goto 160
      end if

      !Calculate the stopping criterion -- dstop.
      d = 0.0d0
      ij = 0
      do i=1,n
         do j=1,i-1
            ij = ij + 1
            d = d + a(ij)**2
         end do
         ij = ij + 1
      end do
      dstop=t1*d

      !Calculate the threshold, thrsh.
      thrsh= sqrt(d/avgf)

      !Start a sweep.
70    continue
      iflag=0
      do jcol=2,n
         jcoltr = jcol*(jcol-1)/2
         jcol1=jcol-1
         i140: do irow=1,jcol1
            indxa = jcoltr + irow
            aij=a(indxa)

            !Compare the off-diagonal element with thrsh.
            daij = abs(aij)
            if(daij.le.thrsh) cycle i140
            irowtr = irow*(irow-1)/2
            indxa = irowtr + irow
            aii=a(indxa)
            indxa = jcoltr + jcol
            ajj=a(indxa)
            s=ajj-aii
            ds = abs(s)

            !The chosen rotation is less than the rounding.
            !Do not rotate.
            if (daij.lt.t2*ds) cycle i140
            iflag=1
            if(t2*daij.ge.ds)then

               !Rotation is very close to 45 degrees,
               !sin and cos = 1/(root 2).

               s = sq2inv
               c = s
            else

               !Rotation is not very close to 45 degrees.
               t = aij/s
               u = 0.25_RK/sqrt(0.25_RK+t*t)
               c = sqrt(0.5_RK+u)
               s = 2.0_RK*t*u/c
            end if

            !Calculate new elements of matrix a.
            do i=1,irow
               t         = a(irowtr + i)
               u         = a(jcoltr + i)
               a(irowtr + i) = c*t-s*u
               a(jcoltr + i) = s*t+c*u
            end do
            i2 = irow+2
            if (i2.le.jcol) then
               do i=i2,jcol
                  itr = (i-1)*(i-2)/2
                  t           = a(jcoltr + i-1)
                  u           = a(itr + irow)
                  a(jcoltr + i-1) = s*u+c*t
                  a(itr + irow) = c*u-s*t
               end do
            end if
            a(jcoltr + jcol) = s*aij+c*ajj
            a(irowtr + irow) = c*a(irowtr + irow)-s*(c*aij-s*ajj)
            do j=jcol,n
               jtr = j*(j-1)/2
               t         = a(jtr + irow)
               u         = a(jtr + jcol)
               a(jtr + irow) = c*t-s*u
               a(jtr + jcol) = s*t+c*u
            end do

            !Rotation completed. see if eigenvectors are wanted by
            !user.
            do i=1,n
               t=eivr(i,irow)
               eivr(i,irow)=c*t-eivr(i,jcol)*s
               eivr(i,jcol)=s*t+eivr(i,jcol)*c
            end do

            !Calculate the new norm d and compare with dstop.
            s=aij
            d=d-s*s
            if(d.lt.dstop) then
               
               !Recalculate dstop and thrsh to discard rounding errors.
               d=0.0_RK
               ij = 0
               do i=1,n
                  do j=1,i-1
                     ij = ij + 1
                     d = d + a(ij)**2
                  enddo
                  ij = ij + 1
               enddo
               dstop=t1*d
            end if
            thrsh=sqrt(d/avgf)
         enddo i140
      enddo
      if(iflag.ne.0) goto 70
160   continue

      !Fill eigenvalue vector.
      idiag = 0
      do i=1,n
         idiag = idiag + i
         w(i) = a(idiag)
      end do
      !if(sorting) then
      if(.true.) then

         !Arrange eigenvalues & vectors in ascending order.
         isort=1
180      continue
         if(isort.eq.1)then
            isort=0
            do i=1,n-1
               if(w(i).gt.w(i+1))then
                  temp=w(i)
                  w(i)=w(i+1)
                  w(i+1)=temp
                  do k = 1,n
                     temp=eivr(k,i)
                     eivr(k,i)=eivr(k,i+1)
                     eivr(k,i+1)=temp
                  end do
                  isort=1
               endif
            end do
            goto 180
         endif
      endif
    end subroutine sjacobi
  end subroutine jacobi

  !--------------- End of module -------------------------------------
end module matrix_eigenval

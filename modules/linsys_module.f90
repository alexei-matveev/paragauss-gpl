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
!================================================================
! Public interface of module
!================================================================
module linsys_module 
  !---------------------------------------------------------------
  !-------------- Module specification ---------------------------
  !---------------------------------------------------------------
  !
  !  Purpose: Provides subroutines for solving the linear equation
  !           systems occuring during the fitting procedures
  !
  !           a) unconstrained fit
  !
  !           Equation to solve:
  !           Sum(j) F_ij a_j = b_i
  !
  !           Subroutine calls required:
  !           CALL decompose_linsys(n,F) -- computes the decomposition of F
  !           CALL solve_linsys(n,F,b)   -- computes the solution a
  !
  !           b) constrained fit
  !
  !           Equations to solve:    
  !
  !           Basic algorithm:
  !           x_j = Sum(i) (F^-1)_ji b_i           <=> Sum(j) F_ij x_j = b_i
  !           c_j = Sum(i) (F^-1)_ji Q_i           <=> Sum(j) F_ij c_j = Q_i
  !           sum = Sum(i) c_i Q_i                 
  !           lam = ( N - Sum(i) c_i b_i ) / sum
  !               = ( N - Sum(j) Q_j x_j ) / sum
  !           a_j = x_j + lam * c_j
  !            
  !           Subroutine calls required:
  !           CALL decompose_linsys(n,F,Q,c) -- computes the decomposition of F
  !                                             the dual contrain vector c and
  !                                             and its F-norm
  !           CALL solve_linsys(n,F,b,c,N)   -- computes the solution a and
  !                                             the Lagrange multiplier lam
  !                                             using Sum(i) c_i b_i
  !       or  CALL solve_linsys(n,F,b,c,N,Q) -- computes the solution a and
  !                                             the Lagrange multiplier lam
  !                                             using Sum(j) Q_j x_j
  !
  !           If the unconstrained solution x = F^-1 b is already known
  !           the call of solve_linsys can be replaced by 
  !
  !           CALL constrain_linsys(n,x,b,c,N) -- computes the solution a and
  !                                               the Lagrange multiplier lam
  !                                               based on the unconstrained
  !                                               solution x = F^-1 b using 
  !                                               using Sum(i) c_i b_i
  !       or  CALL constrain_linsys(n,x,c,N,Q) -- computes the solution a and
  !                                               the Lagrange multiplier lam
  !                                               based on the unconstrained
  !                                               solution x = F^-1 b using 
  !                                               using Sum(i) Q_j x_j
  ! 
  !  Author: Uwe Birkenheuer (extracted from lin_solve by FN)
  !  Date: 7/7/97
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
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
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module
  implicit none
  private
  save
  !== Interrupt end of public interface of module =================
  !------------ public functions and subroutines ------------------
  public decompose_linsys, solve_linsys, constrain_linsys_cb, constrain_linsys_Qx 
  public gmat_linsys

  !------------ Interfaces ----------------------------------------
  interface constrain_linsys
     module procedure constrain_linsys_cb
     module procedure constrain_linsys_Qx
  end interface
  public constrain_linsys

!================================================================
! End of public interface of module
!================================================================

  !----------------------------------------------------------------
  ! LAPACK Routines used 
  !----------------------------------------------------------------
  ! SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
  ! CHARACTER          UPLO
  ! INTEGER            INFO, N
  ! DOUBLE PRECISION   AP( * )
  ! Purpose:
  ! DPPTRF computes the Cholesky factorization of a real symmetric
  ! positive definite matrix A stored in packed format.
  ! The factorization has the form
  !    A = U**T * U,  if UPLO = 'U', or
  !    A = L  * L**T,  if UPLO = 'L'
  ! where U is an upper triangular matrix and L is lower triangular.
  ! Arguments:
  ! UPLO    (input) CHARACTER*1
  !         = 'U':  Upper triangle of A is stored
  !         = 'L':  Lower triangle of A is stored.
  ! N       (input) INTEGER
  !         The order of the matrix A.  N >= 0.
  ! AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
  !         On entry, the upper or lower triangle of the symmetric matrix
  !         packed columnwise in a linear array.  The j-th column of A
  !         is stored in the array AP as follows:
  !         if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j
  !         if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
  !         On exit, if INFO = 0, the triangular factor U or L from the
  !         Cholesky factorization A = U**T*U or A = L*L**T, in the same
  !         storage format as A.
  ! INFO    (output) INTEGER
  !         = 0:  successful exit
  !         < 0:  if INFO = -i, the i-th argument had an illegal value
  !         > 0:  if INFO = i, the leading minor of order i is not
  !               positive definite, and the factorization could not be
  !               completed.
  !----------------------------------------------------------------
  ! SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
  ! CHARACTER          UPLO
  ! INTEGER            INFO, LDB, N, NRHS
  ! DOUBLE PRECISION   AP( * ), B( LDB, * )
  ! Purpose:
  ! DPPTRS solves a system of linear equations A*X = B with a symmetric
  ! positive definite matrix A in packed storage using the Cholesky
  ! factorization A = U**T*U or A = L*L**T computed by DPPTRF.
  ! Arguments:
  ! UPLO    (input) CHARACTER*1
  !         = 'U':  Upper triangle of A is stored
  !         = 'L':  Lower triangle of A is stored.
  ! N       (input) INTEGER
  ! The order of the matrix A.  N >= 0.
  ! NRHS    (input) INTEGER
  !         The number of right hand sides, i.e., the number of columns
  !         of the matrix B.  NRHS >= 0.
  ! AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
  !         The triangular factor U or L from the Cholesky factorization
  !         A = U**T*U or A = L*L**T, packed columnwise in a linear
  !         array.  The j-th column of U or L is stored in the array AP
  !         as follows:
  !         if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j
  !         if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
  ! B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
  !         On entry, the right hand side matrix B.
  !         On exit, the solution matrix X.
  ! LDB     (input) INTEGER
  !         The leading dimension of the array B.  LDB >= max(1,N).
  ! INFO    (output) INTEGER
  !         = 0:  successful exit
  !         < 0:  if INFO = -i, the i-th argument had an illegal value
  !----------------------------------------------------------------

  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************

  subroutine decompose_linsys(n_dim,f_mat,q_vec,c_vec)
    ! Purpose: Decompose the n x n fit overlap matrix F
    !          and prepare the constrain correction - lam * c
    ! Input  : n_dim      -- the dimension n of the linear equation system
    !          f_mat(:)   -- the overlap matrix F_ij in packed storage format
    !          q_vec(:)   -- the constrain contributions Q_i
    ! Output : f_mat(:)   -- the decomposed F matrix in packed storage format
    !          c_vec(1:n) -- the transformed constrain contributions c = F^-1 Q
    !          c_vec(n+1) -- the normalization Sum(i) c_i Q_i
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind)          , intent(in   ) :: n_dim
    real   (kind=r8_kind)          , intent(inout) :: f_mat(:)
    real   (kind=r8_kind), optional, intent(in   ) :: q_vec(:)
    real   (kind=r8_kind), optional, intent(  out) :: c_vec(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: info, n_sum !,i,j
!   real   (kind=r8_kind)::summ
    !------------ Executable code --------------------------------

!  summ=0.0_r8_kind
!  do j=1,n_dim
!   do i=1,j
!   summ=summ+f_mat(i + (j-1)*j/2)
!   if(j.ne.i) summ=summ+f_mat(i + (j-1)*j/2)
!   enddo
!  enddo
!   print*,'decompose linsys sum f_mat', summ,sum(f_mat)

    
    call dpptrf('U',n_dim,f_mat,info)
    if (info.lt.0) then
       call error_handler &
            ("DECOMPOSE_LINSYS: dpptrf detected illegal argument")
    else if (info.gt.0) then
       call error_handler &
            ("DECOMPOSE_LINSYS: dpptrf detected non positive definite matrix")
    endif
   
    if (present(q_vec) .and. present(c_vec)) then
       n_sum = n_dim+1
       c_vec(1:n_dim) = q_vec

       ! SOLVE LINEAR SYSTEM: f_mat * X = c_vec
       ! and put X -> c_vec
       call dpptrs('U',n_dim,1,f_mat,c_vec,n_sum,info)

       if (info.lt.0) then
          call error_handler &
               ("DECOMPOSE_LINSYS: dpptrs detected illegal argument")
       else if (info.gt.0) then
          call error_handler &
               ("DECOMPOSE_LINSYS: unknown error situation after dpptrs")
       endif
!      print*,sum(c_vec(1:n_dim)),'decompose sum of a", solutions of auxilary system'
       c_vec(n_sum) = sum(c_vec(1:n_dim)*q_vec)
    endif
   
  end subroutine decompose_linsys

  subroutine gmat_linsys(n_dim, f_inv, q_vec, c_vec, g_mat, zmat, matcharge)
    !
    ! Purpose:  construct cpks  equiations gmat  from  the constrained
    ! charge fit equations using the decomposed overlap matrix and the
    ! transformed constrain  contributions computed by  the subroutine
    ! DECOMPOSE_LINSYS
    !
    ! Input  : n_dim      -- the dimension n of the linear equation system
    !          f_inv(:)   -- the decomposed F matrix in packed storage format
    !          q_vec(:)   -- the original constrain contributions Q_i
    ! Output : zmat       -- inverted matcharge matrix
    use cpksdervs_matrices, only: cpksalloc, cpks_gmat_calculated
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind)          , intent(in   ) :: n_dim
    real   (kind=r8_kind)          , intent(in   ) :: f_inv(:)
    real   (kind=r8_kind), optional, intent(in   ) :: matcharge(:)
    real(r8_kind), intent(in) :: q_vec(:)
    real   (kind=r8_kind), optional, intent(in   ) :: c_vec(:)
    real   (kind=r8_kind), intent(out) :: g_mat(n_dim,n_dim)
    real(kind=r8_kind),intent(out) :: zmat(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: info,i,j
    real(kind=r8_kind),allocatable :: fmat(:,:),matcha(:,:)
    !------------ Executable code --------------------------------
  if(present(matcharge) ) then
   allocate(matcha(n_dim,n_dim))
   do j=1,n_dim
    do i=1,j
    matcha(i,j)=matcharge(i + (j-1)*j/2)
    matcha(j,i)=matcharge(i + (j-1)*j/2)
    enddo
   enddo
   endif

    call dpptri('U',n_dim,f_inv,info)
    if (info.lt.0) then
       call error_handler &
            ("GMAT_LINSYS: dpptri detected illegal argument")
    else if (info.gt.0) then
       call error_handler &
            ("GMAT_LINSYS: unknown error situation after dpptri")
    endif
   allocate(fmat(n_dim,n_dim),stat=cpksalloc(18))
   ASSERT(cpksalloc(18).eq.0)
   
!  print*,sum(c_vec(:n_dim)),'linsys gmat : sum a"'
   do j=1,n_dim
    do i=1,j
             !   a"       Q        sum Qxa"
    fmat(i,j)=-c_vec(i)*q_vec(j)/c_vec(n_dim+1)
    fmat(j,i)=-c_vec(j)*q_vec(i)/c_vec(n_dim+1)
    enddo
    fmat(j,j)=fmat(j,j)+1.0_r8_kind
   enddo
   do j=1,n_dim
    do i=1,j
    zmat(i,j)=f_inv(i + (j-1)*j/2)
    zmat(j,i)=f_inv(i + (j-1)*j/2)
    enddo
   enddo

   g_mat=matmul(fmat,zmat)
   ! print*,'g_mat fmat zmat', sum(g_mat),sum(fmat),sum(zmat)
   deallocate(fmat,stat=cpksalloc(18))
   ASSERT(cpksalloc(18).eq.0)
   cpksalloc(18)=1
   cpks_gmat_calculated=.true.

   if(present(matcharge) ) then
    ! print*,'zmat matcha check sum',sum(matmul(matcha,zmat))
    deallocate(matcha)
   endif
  end subroutine gmat_linsys

  subroutine solve_linsys(n_dim,f_inv,b_vec,c_vec,n_val,q_vec)
    ! Purpose: Solve the constrained charge fit equations using the decomposed
    !          overlap matrix and the transformed constrain contributions
    !          computed by the subroutine DECOMPOSE_LINSYS
    ! Input  : n_dim      -- the dimension n of the linear equation system
    !          f_inv(:)   -- the decomposed F matrix in packed storage format
    !          b_vec(:)   -- the right hand side b to the constrained fit
    !          c_vec(1:n) -- the transformed constrain contributions c = F^-1 Q
    !          c_vec(n+1) -- the normalization Sum(i) c_i Q_i
    !          n_val      -- the constrain target N
    !          q_vec(:)   -- the original constrain contributions Q_i
    ! Output : b_vec(:)   -- the solutions a of the constrained fit
    !          n_val      -- the Lagrange multiplier lambda
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind)          , intent(in   ) :: n_dim
    real   (kind=r8_kind)          , intent(in   ) :: f_inv(:)
    real   (kind=r8_kind)          , intent(inout) :: b_vec(:)
    real   (kind=r8_kind), optional, intent(in   ) :: c_vec(:)
    real   (kind=r8_kind), optional, intent(inout) :: n_val
    real   (kind=r8_kind), optional, intent(in   ) :: q_vec(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: info, n_sum
    !------------ Executable code --------------------------------

    if (present(c_vec) .and. present(n_val) .and..not.present(q_vec)) then
      n_sum = n_dim + 1
      n_val = ( n_val - sum(c_vec(1:n_dim)*b_vec) ) / c_vec(n_sum)
    endif
   
    !
    ! Compute b_vec (charge fit coefficients) by an unconstrained density fit
    ! procedure:
    !
    if ( n_dim > 0 ) then
        ! DPPTRS does not cope with zero-sized arrays ...
        call dpptrs('U', n_dim, 1, f_inv, b_vec, n_dim, info)
        ASSERT(info==0)
    else
        b_vec(:n_dim) = -1.0 ! NOOP, as lhs is zero sized.
    endif

    if (present(c_vec) .and. present(n_val)) then
      if (present(q_vec)) then
         n_sum = n_dim + 1
         n_val = ( n_val - sum(q_vec*b_vec) ) / c_vec(n_sum)
      endif
      b_vec = b_vec + n_val * c_vec(1:n_dim)
    endif
   
  end subroutine solve_linsys

  

  subroutine constrain_linsys_cb(n_dim,x_vec,b_vec,c_vec,n_val)
    ! Purpose: Solve the constrained charge fit equations using the solution 
    !          from the corresponding unconstrained fit the transf. constrain
    !          contributions computed by the subroutine DECOMPOSE_LINSYS
    ! Input  : n_dim      -- the dimension n of the linear equation system
    !          x_vec(:)   -- the unconstrained solution x = F^-1 b
    !          b_vec(:)   -- the right hand side b to the constrained fit
    !          c_vec(1:n) -- the transformed constrain contributions c = F^-1 Q
    !          c_vec(n+1) -- the normalization Sum(i) c_i Q_i
    !          n_val      -- the constrain target N

    ! Output : b_vec(:)   -- the solutions a of the constrained fit
    !          n_val      -- the Lagrange multiplier lambda
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in   ) :: n_dim
    real   (kind=r8_kind), intent(in   ) :: x_vec(:)
    real   (kind=r8_kind), intent(inout) :: b_vec(:)
    real   (kind=r8_kind), intent(in   ) :: c_vec(:)
    real   (kind=r8_kind), intent(inout) :: n_val
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: n_sum
    !------------ Executable code --------------------------------

    n_sum = n_dim + 1
    n_val = ( n_val - sum(c_vec(1:n_dim)*b_vec) ) / c_vec(n_sum)

    b_vec = x_vec + n_val * c_vec(1:n_dim)
   
  end subroutine constrain_linsys_cb


  subroutine constrain_linsys_Qx(n_dim,x_vec,c_vec,n_val,q_vec)
    ! Purpose: Solve the constrained charge fit equations using the solution 
    !          from the corresponding unconstrained fit the transf. constrain
    !          contributions computed by the subroutine DECOMPOSE_LINSYS
    ! Input  : n_dim      -- the dimension n of the linear equation system
    !          x_vec(:)   -- the unconstrained solution x = F^-1 b
    !          c_vec(1:n) -- the transformed constrain contributions c = F^-1 Q
    !          c_vec(n+1) -- the normalization Sum(i) c_i Q_i
    !          n_val      -- the constrain target N
    !          q_vec(:)   -- the original constrain contributions Q_i
    !
    ! Output : x_vec(:)   -- the solutions a of the constrained fit
    !          n_val      -- the Lagrange multiplier lambda
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in   ) :: n_dim
    real   (kind=r8_kind), intent(inout) :: x_vec(:)
    real   (kind=r8_kind), intent(in   ) :: c_vec(:)
    real   (kind=r8_kind), intent(inout) :: n_val
    real   (kind=r8_kind), intent(in   ) :: q_vec(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: n_sum
    !------------ Executable code --------------------------------

    n_sum = n_dim + 1
    n_val = ( n_val - sum(q_vec*x_vec) ) / c_vec(n_sum)

    x_vec = x_vec + n_val * c_vec(1:n_dim)
   
  end subroutine constrain_linsys_Qx


end module linsys_module

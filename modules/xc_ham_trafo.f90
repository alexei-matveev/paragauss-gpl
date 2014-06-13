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
module xc_ham_trafo
  !---------------------------------------------------------------
  !
  !  Purpose: Transform vector representaion into
  !           the spinor representation.
  !
  ! The vector/projective irrep-based ordering of the
  ! basis functions emphasizes the block structure of
  ! spin-free operators (and facilitates bulk operation on matrices)
  !
  !   each projective basis: pIrr=1,n_p_irr
  !     consists of vector irrep origins: vIrr=1,n_v_irr
  !       (such that: vIrr x SU(2) = pIrr)
  !         from each vector basis function: vf=1,dim(vIrr)
  !           obtain new projective basis function: pf++
  ! graphically:
  ! _________________________
  ! |      p1    |     p2   |  (1)*
  ! | v1   | v2  |  v1 | v2 |
  !
  ! *) note grouping of (vi)pj basis functions
  !
  !
  ! However PG uses different order of the basis functions
  !
  !   each irrep basis: Irr=1,n_irr
  !     is obtained by reduction of ua x L shells: uaL=1,uaL_max
  !       for each instanse of this reduction: mult=1,mult(uaL->Irr)
  !         for each radial contraction: e=1,n_e(uaL)
  !            obtain new basis function: f++
  !
  ! This last scheme is employed for both vector and projective irreps.
  ! Graphically it corresponds to
  ! _________________________
  ! |     v1    |     v2    |  (2)
  ! | uaL1 uaL2 | uaL1 uaL2 |
  !
  ! and
  ! _________________________
  ! |     p1    |     p2    |  (3)
  ! | uaL1 uaL2 | uaL1 uaL2 |
  !
  ! Even though projective irreps are generated via the
  ! product of SU(2) and reduced to vector irreps uaL space
  ! the clear block structure of (1) is hidden:
  ! ________________________________
  ! |       p1      |       p2      |
  ! |  uaL1 |  uaL2 |  uaL1 | uaL2  | (4)*
  ! | v1 v2 | v1 v2 | v1 v2 | v1 v2 |
  !
  ! *) note fragmentation of (vi)pj functions: there are
  !    two ranges of (v1)p1 that differ by uaL
  !
  ! This module operates on blocks assuming (1) and then
  ! by calling permute_matrix() brings it to (4)
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
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: so_trafo

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine so_trafo(pham, vham)
    use dimensions, only: vdims => IrrBasDim, pdims => IrrBasDimSpor
    use clebsch_gordan, only: cg => vsu2cg_eliminated
    use datatype, only: arrmat2
    implicit none
    real(RK), intent(out) :: pham(:) ! xc_length_proj
    real(RK), intent(in) :: vham(:) ! xc_length_vec
    ! *** end of interface ***

    integer(IK) :: i
    type(arrmat2) :: ham_vec(size(vdims)) ! (n_virr)
    type(arrmat2) :: ham_proj(size(pdims)) ! (n_pirr)

    !
    ! Unpack the block diagonal orbital representation of the
    ! hamiltonian matrix:
    !
    call unpack_block_diagonal(vdims, vham, ham_vec)

    !
    ! Build projective matrices:
    !
    do i = 1, size(ham_proj)
       ham_proj(i) = xsu2(cg(i, :)%mult, ham_vec(:))
       ASSERT(dim(ham_proj(i))==pdims(i))
    enddo

    !
    ! This reshuffles elements, pass multiplicities through:
    !
    call reshuffle(cg(:, :)%mult, ham_proj)

    !
    ! Pack blocks into a linear array:
    !
    call pack_block_diagonal(ham_proj, pham)
  end subroutine so_trafo

  function dim(m) result(n)
    !
    ! Dimension of a square arrmat2
    !
    use datatype, only: arrmat2
    implicit none
    type(arrmat2), intent(in) :: m
    integer(IK) :: n
    ! *** end of interface ***

    n = size(m%m, 1)
    ASSERT(n==size(m%m,2))
  end function dim

  function xsu2(cg, ham_vec) result(mat_proj)
    !
    ! Distributes matrix elements over projective matrix.
    !
    ! Array cg(:) is supposed to hold multiplicities of the direct
    ! product:
    !
    !       vector irrep * SU(2)
    !               -> projective irrep cg(vector irrep) times 
    !
    use datatype, only: arrmat2
    implicit none
    integer(IK), intent(in) :: cg(:) ! (n_virr)
    type(arrmat2), intent(in) :: ham_vec(:) ! (n_virr)
    type(arrmat2) :: mat_proj
    ! *** end of interface ***

    integer(IK) :: i, k, n, m, off

    ASSERT(size(cg)==size(ham_vec))

    !
    ! First compute the dimension of mat_proj:
    !
    n = 0
    do i = 1, size(ham_vec) ! n_virr
       n = n + dim(ham_vec(i)) * cg(i)
    enddo

    !
    ! Allocate output matrix:
    !
    allocate(mat_proj%m(n, n))

    !
    ! The off diagonal blocks will not be set, zero them:
    !
    mat_proj%m = 0.0

    !
    ! The vector representaiton "ham_vec" is block-diagonal,
    ! so run only over the diagonal pairs of vector irreps (i, i):
    !
    off = 0
    do i = 1, size(ham_vec)
       m = dim(ham_vec(i))
       !
       ! Here we assume that blocks off-diagonal in
       ! irrep instance indices (k1, k2) vanish too.
       !
       ! FIXME: cg(i) is only > 1 in C1, we should somehow make
       !        obvious that this fact holds for C1
       !
       do k = 1, cg(i)
          mat_proj%m(off+1:off+m, off+1:off+m) = ham_vec(i)%m

          !
          ! Advance the offset counter:
          !
          off = off + m
       enddo
    enddo
    ASSERT(off==dim(mat_proj))
  end function xsu2

  subroutine unpack_block_diagonal(dims, packed, unpacked)
    use datatype, only: arrmat2
    implicit none
    integer(IK), intent(in) :: dims(:) ! (nirr)
    real(RK), intent(in) :: packed(:) ! (sum(dims*(dims+1)/2))
    type(arrmat2), intent(out) :: unpacked(:) ! (nirr)
    ! *** end of interface ***

    integer(IK) :: i, a, b, ab, n

    ab = 0
    do i = 1, size(dims)
       n = dims(i)

       allocate(unpacked(i)%m(n, n))

       do b = 1, n
          do a = 1, b
             ab = ab + 1
             unpacked(i)%m(a, b) = packed(ab)
             unpacked(i)%m(b, a) = packed(ab)
          enddo
       enddo
    enddo
    ASSERT(ab==size(packed))
  end subroutine unpack_block_diagonal

  subroutine pack_block_diagonal(unpacked, packed)
    use datatype, only: arrmat2
    implicit none
    type(arrmat2), intent(in) :: unpacked(:) ! (nirr)
    real(RK), intent(out) :: packed(:) ! (sum(dims*(dims+1)/2))
    ! *** end of interface ***

    integer(IK) :: i, a, b, ab

    ab = 0
    do i = 1, size(unpacked)
       do b = 1, dim(unpacked(i))
          do a = 1, b
             ab = ab + 1
             packed(ab) = unpacked(i)%m(a, b)
          enddo
       enddo
    enddo
    ASSERT(ab==size(packed))
  end subroutine pack_block_diagonal

  subroutine reshuffle(cg, ham_proj)
    use datatype, only: arrmat2
    implicit none
    integer(IK), intent(in) :: cg(:, :) ! (n_pirr, n_virr)
    type(arrmat2), intent(inout) :: ham_proj(:) ! (n_pirr)
    ! *** end of interface ***

    integer(IK) :: i
    integer(IK), allocatable :: order(:) ! (uaL_max * sum(cg(i, :)))
    integer(IK), allocatable :: sizes(:) ! (uaL_max * sum(cg(i, :)))


    do i = 1, size(ham_proj)
       call compute_permutation(cg(i, :), order, sizes)
       call permute_matrix(ham_proj(i)%m, order, sizes)
    enddo
  end subroutine reshuffle

  subroutine compute_permutation(cg, order, sizes)
    use dimensions, only: uaL_max, uaL_vec_mult, uaL_nrad
    implicit none
    integer(IK), intent(in) :: cg(:) ! (n_virr)
    integer(IK), allocatable, intent(out) :: order(:) ! (uaL_max*??)
    integer(IK), allocatable, intent(out) :: sizes(:) ! (uaL_max*??)
    ! *** end of interface ***

    integer(IK) :: shell, irrep, i, k, block
    integer(IK) :: off(size(cg)) ! (n_virr)

!!$ integer(IK), allocatable :: seq1(:, :) ! (4, blocks)
    integer(IK), allocatable :: seq2(:, :) ! (4, blocks)

    !
    ! Number of spinors "cg(irrep)" that can be obtained from an orbital
    ! of irrep "irrep" varies. These "offsets" will help us addressing
    ! the shells:
    !
    off(1) = 0
    do irrep = 2, size(cg)
       off(irrep) = off(irrep - 1) + cg(irrep - 1)
    enddo

    !
    !
    ! We will be permuting in blocks. A block corresponds to a tuple
    !
    !       (shell, irrep, i, k)
    !
    ! with
    !
    ! shell -- cumulative shell index (ua, L)
    ! irrep -- orbital symmetry
    ! i     -- instance of that irrep, in case there are multile
    ! k     -- instance of the projective obtained by product with SU(2)
    !
    block = 0
    do shell = 1, uaL_max
       do irrep = 1, size(cg)
          do i = 1, uaL_vec_mult(shell, irrep)
             do k = 1, cg(irrep)
                block = block + 1
             enddo
          enddo
       enddo
    enddo
    allocate(order(block), sizes(block))

!!$ allocate(seq1(4, block))
!!$ !
!!$ ! Seq1: output sequence
!!$ !
!!$ block = 0
!!$ do shell = 1, uaL_max
!!$ do irrep = 1, size(cg)
!!$ do i = 1, uaL_vec_mult(shell, irrep)
!!$ do k = 1, cg(irrep)
!!$     block = block + 1
!!$     seq1(1, block) = k
!!$     seq1(2, block) = i
!!$     seq1(3, block) = irrep
!!$     seq1(4, block) = shell
!!$ enddo
!!$ enddo
!!$ enddo
!!$ enddo
!!$ ASSERT(block==size(order))

    !
    ! Seq2: input sequence
    !
    allocate(seq2(4, block))
    block = 0
    do irrep = 1, size(cg)
       do k = 1, cg(irrep)
          do shell = 1, uaL_max
             do i = 1, uaL_vec_mult(shell, irrep)
                block = block + 1
                seq2(1, block) = k
                seq2(2, block) = i
                seq2(3, block) = irrep
                seq2(4, block) = shell
             enddo
          enddo
       enddo
    enddo
    ASSERT(block==size(order))

    !
    ! This loop structure is intended to reproduce the
    ! ordering of basis spinors used by the rest of PG:
    !

    !
    ! For orbitals originating from each atomic shell (ua, L):
    !
    block = 0
    do shell = 1, uaL_max
       !
       ! ... classified by their vector irrep lables:
       !
       do irrep = 1, size(cg) ! n_virr
          !
          ! ... with zero or more possible orbital combinations of that irrep:
          !
          do i = 1, uaL_vec_mult(shell, irrep)
             !
             ! ... there are zero or more possibilities to construct spinors:
             !
             do k = 1, cg(irrep)
                !
                ! ... from each of uaL_nrad(shell, irrep) radial shapes:
                !
                block = block + 1

                !
                ! Source block from position src(shell, irrep, k) goes to destination
                ! block at position "block" == dst(shell, irrep, k):
                !
                order(block) = src(shell, irrep, i, k)

                !
                ! Size of the block:
                !
                sizes(order(block)) = uaL_nrad(shell)
             enddo
          enddo
       enddo
    enddo

  contains

    function src(shell, irrep, i, k) result(seq)
      !
      ! Source of the block (shell, irrep, i, k)
      !
      implicit none
      integer(IK), intent(in) :: shell, irrep, i, k
      integer(IK) :: seq
      ! *** end of interface ***

      integer(IK) :: j

      seq = -1
      do j = 1, size(seq2, 2)
         if (k     == seq2(1, j) .and. &
              i     == seq2(2, j) .and. &
              irrep == seq2(3, j) .and. &
              shell == seq2(4, j)) then
            seq = j
            exit ! loop
         endif
      enddo
      ASSERT(seq/=-1)
    end function src

  end subroutine compute_permutation

    subroutine permute_matrix(mat, ord, siz)
      !
      ! In blocks of size siz(i) x siz(j) do:
      !
      !         block(i, j) <- block(ord(i), ord(j))
      !
      ! Example: if order = [2, 3, 1] and siz = [1, 1, 1]
      ! then reshuffle the elements as shown below:
      !
      !         mat(in)        mat(out)
      !         22 23 21       11 12 13
      !         32 33 31  ->   21 22 23
      !         12 13 11       31 32 33
      !
      implicit none
      real(RK), intent(inout) :: mat(:,:) ! (n, n) 
      integer(IK), intent(in) :: ord(:) ,siz(:) ! (m)
      ! *** end of interface ***

      integer(IK) :: i
      integer(IK) :: counter
      integer(IK) :: ii,ai,ei,aii,eii
      real(RK)    :: wrk(size(mat,1),size(mat,2))
      integer(IK) :: off(size(ord))
      integer(IK) :: n,m

      n = size(mat,1)
      ASSERT(n==size(mat,2))

      m = size(ord)
      ASSERT(m==size(siz))
      ASSERT(n==sum(siz))

      ! calculate offsets:
      counter = 0
      do i=1,m
         off(i)  = counter
         counter = counter + siz(i)
      enddo
      ASSERT(counter==n)

      ! permute the second index
      ! into a temp storage:
      counter = 0
      do i=1,m
         ii  = ord(i)
         aii = off(ii) + 1
         eii = off(ii) + siz(ii)
         ai  = counter + 1
         ei  = counter + siz(ii)
         DWRITE(*,'("xcht/xc_so_trafo: ",I3,A4,I3,A2,I3,A1)') ii,' -> ',i, ' (',siz(ii),')'

         wrk(:,ai:ei) = mat(:,aii:eii)

         counter = counter + siz(ii)
      enddo
      ASSERT(counter==n)

      ! permute the first index
      ! back into original storage:
      counter = 0
      do i=1,m
         ii  = ord(i)
         aii = off(ii) + 1
         eii = off(ii) + siz(ii)
         ai  = counter + 1
         ei  = counter + siz(ii)

         mat(ai:ei,:) = wrk(aii:eii,:)

         counter = counter + siz(ii)
      enddo
      ASSERT(counter==n)
    end subroutine permute_matrix

  !--------------- End of module ----------------------------------
end module xc_ham_trafo

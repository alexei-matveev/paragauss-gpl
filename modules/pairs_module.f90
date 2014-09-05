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
  !===================================================================
! Public interface of module
  !===================================================================
module pairs_module
!-------------- Module specification ---------------------------
!
!  Purpose: Contains necessary data for the pairs of
!           of the perturbation theory
!
!  Contents:
!    - Variables:  ptpair1
!                  ptpair2
!                  ptweit
!                  n_pairs
!                  ptlim
!                  npair0
!    - Subroutines:
!                  ptpairs
!                  deallocate_pairs
!
!  Author: TG
!  Date: 11/10/96
!
!
!---------------------------------------------------------------------
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: UB
! Date:   7/97
! Description: Perturbation theory is now based on eigval_occ and
!              occ_num_occ from the occupied_levels_module and on
!              eigval_vir from the virtual_levels_module.
!
! Modification (Please copy before editing)
! Author: MM
! Date:   6/98
! Description: Extension to Spin Orbit
!
! Modification
! Author: AM
! Date:   03.05.99
! Description: deallocate_pairs sub modified
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!---------------------------------------------------------------------
#include "def.h"
use type_module, only: i4_kind, r8_kind ! type specification parameters
use datatype, only: arrmat2, arrmat2int
use options_module, only: options_spin_orbit

implicit none
private         ! by default, all names are private
save

!== Interrupt end of public interface of module =================

type(arrmat2int), allocatable, public, protected :: ptpair1(:), ptpair2(:)
type(arrmat2), allocatable, public, protected :: ptweit(:) ! PT-weight, a misnomer
integer(i4_kind), allocatable, public, protected :: n_pairs0(:,:)
integer(kind=i4_kind), public, parameter :: n_pairs = 10
real(kind=r8_kind), public, parameter :: ptlim = 5.0_r8_kind


!------------ public functions and subroutines ------------------

public :: ptpairs!()
public :: deallocate_pairs!()


  !===================================================================
  ! End of public interface of module
  !===================================================================

!---------------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


   !*************************************************************
  subroutine allocate_pairs()
    ! Purpose: Short allocation routine for pairs stuff
    use symmetry_data_module, only: ssym
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i_i
    ! dimensions
    integer(kind=i4_kind)             :: n_irrep,n_spin
    integer(kind=i4_kind),allocatable :: dim_irrep(:)
    !------------ Executable code ------------------------------------

    ! set appropriate dimensions of irreps
    n_spin = ssym%n_spin
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do i_i=1,n_irrep
          dim_irrep(i_i) = ssym%dim_proj(i_i)
       enddo
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do i_i=1,n_irrep
          dim_irrep(i_i)  = ssym%dim(i_i)
       enddo
    endif ! options_spin_orbit

    allocate(ptweit(n_irrep))
    allocate(ptpair1(n_irrep))
    allocate(ptpair2(n_irrep))
    allocate(n_pairs0(n_spin,n_irrep))
    do i_i=1,n_irrep
       allocate(ptpair1(i_i)%m(n_pairs,n_spin))
       allocate(ptpair2(i_i)%m(n_pairs,n_spin))
       allocate(ptweit(i_i)%m(n_pairs,n_spin))
    enddo

    ! deallocate dimensions of irreps again
    deallocate(dim_irrep)
  end subroutine allocate_pairs

  subroutine ptpairs()
    !  Purpose: Calculates the pairs of the perturbation theory
    !  and the associated weights of one IRREP
    !
    !  This subrotine runs in a parallel context.
    !
    !  Subroutine called by: pert_coeff_module
    !
    !  Author: T. Grauschopf
    !
    !  Date: 7/95
    !--------------------------------------------------------------
    !
    !  Documentation
    !
    !  The perturbative prefactor of a pair of eigenstates (i,j) is given by
    !
    !     w(i,j) := ( n_i - n_j ) / MAX( | eps_i - eps_j | / 2 , 0.01 ) .
    !
    !  Based on this weight those  N_pair pairs (i,j) with the largest absoule
    !  weights |w(i,j)| are selected.
    !
    !  If the eigenvalues are sorted  increasingly, i.g. eps_i+1 >= eps_i, and
    !  N_occo and N_virt  are the highest (partially) occupied  and the lowest
    !  not fully occupied orbital, respectively, then
    !
    !     w(i<=N_occo,j>N_occo) = 2 * n_i / MAX( eps_j - eps_i , 0.02 )
    !
    !  is monotonically decreasing with increasing j, and thus for any given i
    !  <= N_occo  the heighest virtual orbital  which has to  be considered in
    !  the search is the (N_occo + N_pair)-th one. Similary,
    !
    !     w(i<N_virt,j>=N_virt) = 2 * ( 1 - n_j ) / MAX( eps_j - eps_i , 0.02 )
    !
    !  is monotonically decreasing with decreasing i, and thus for any given j
    !  >= N_virt the  lowest occupied orbital which has to  be included in the
    !  search is  the (N_virt  - N_PAIR)-th one.  Finally if  both, i and  j >
    !  N_occo  or <  N_virt  w(i,j)  = 0  since  the corresponding  occupation
    !  numbers n_i and n_j are equal.
    !
    !  Therefore  it is  sufficient  to  consider only  pairs  (i,j) from  the
    !  following simple and and rather limited index set:
    !
    !  { (i,j) | N_virt - N_pair <= i < j <= N_occo + N_pair }
    !
    !  Author: Uwe Birkenheuer
    !  Date:   7/97
    !--------------------------------------------------------------
    use comm, only: comm_rank
    use symmetry_data_module
    use datatype
    use type_module
    use output_module, only: output_chargefit
    use iounitadmin_module, only: output_unit
    use occupation_module, ONLY : n_occo
    use occupied_levels_module, ONLY : eigval_occ, occ_num_occ
    use virtual_levels_module, only: virtual_levels_bcast, eigval_vir
    implicit none
    !** End of interface *****************************************

    !------------ Declaration of subroutines ------------------------
    intrinsic max,abs

    !------------ Declaration of local constants --------------------
    real(kind=r8_kind), parameter :: eps=4.0D-4

    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) i_i,i_s,i_p,i_o,j_o
    real(kind=r8_kind) dfctrs,ptmin,deigis,eigval
    integer(i4_kind) :: i, j, k
    ! dimensions
    integer(kind=i4_kind)             :: n_irrep,n_spin
    integer(kind=i4_kind),allocatable :: dim_irrep(:)
    !------------ Executable code ------------------------------------

    !
    ! Load eigvec_vir and eigval_vir on each slave.
    !
    call virtual_levels_bcast(n_pairs)

    !
    ! Was done from outside before:
    !
    call allocate_pairs()

    !
    ! FIXME: why dont we let slaves compute them too?
    !        Verify if pre-requisites are available.
    !
    if ( comm_rank() /= 0 ) goto 999 ! receive at broadcast and exit

    ! set appropriate dimensions of irreps
    n_spin = ssym%n_spin
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do i_i=1,n_irrep
          dim_irrep(i_i) = ssym%dim_proj(i_i)
       enddo
    else
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do i_i=1,n_irrep
          dim_irrep(i_i)  = ssym%dim(i_i)
       enddo
    endif

    !  search the n_spin*n_pairs pairs of basis function indices
    !  whose associated weights are biggest.
    irrep: do i_i=1,n_irrep
       ptweit(i_i)%m(:, :) = 0.0_r8_kind
       ptpair1(i_i)%m(:, :) = 0 ! integers
       ptpair2(i_i)%m(:, :) = 0 ! integers
       a: do i_s=1,n_spin
          b: do i_o = 1,min(dim_irrep(i_i)-1,n_occo(i_s,i_i))

             !  the upper bound is a little bit tricky. this is
             !  due to fact, that the calculated weight differs from
             !  zero only in the case that one orbital is occupied and
             !  the second not. The pairs are ordered so that
             !  the associated weights are ascendent:
             !  ptweit(n_pairs,i_s)->max!

             c: do j_o=i_o+1,min(dim_irrep(i_i),n_occo(i_s,i_i)+n_pairs)
                if (j_o > n_occo(i_s,i_i)) then ! occ_num(i_i)%m(j_o,i_s) = 0
                   dfctrs=occ_num_occ(i_i)%m(i_o,i_s)
                   eigval=eigval_vir(i_i)%m(j_o-n_occo(i_s,i_i),i_s)
                else
                   dfctrs=occ_num_occ(i_i)%m(i_o,i_s) - &
                          occ_num_occ(i_i)%m(j_o,i_s)
                   eigval=eigval_occ(i_i)%m(j_o,i_s)
                endif
                if (abs(dfctrs) <= eps) cycle ! screening

                deigis=max(1.0e-2_r8_kind,0.5_r8_kind*abs(eigval - &
                     eigval_occ(i_i)%m(i_o,i_s)))
                dfctrs=dfctrs/deigis

                !
                ! If the weight is larger  or equal to the smalles one from
                ! the list, update the list:
                !
                if (dfctrs < ptweit(i_i)%m(1, i_s)) cycle ! no need to update

                ! take pair (i_o,j_o) and sort it in accord. to w(i_o,j_o)
                do i_p = 2, n_pairs
                   if (dfctrs < ptweit(i_i)%m(i_p, i_s)) exit ! the shift loop

                   !
                   ! Shift to the left:
                   !
                   ptweit(i_i)%m(i_p - 1, i_s) = ptweit(i_i)%m(i_p, i_s)
                   ptpair1(i_i)%m(i_p - 1, i_s) = ptpair1(i_i)%m(i_p, i_s)
                   ptpair2(i_i)%m(i_p - 1, i_s) = ptpair2(i_i)%m(i_p, i_s)
                enddo

                !
                ! Insert new entry, eventually into the rightmost position:
                !
                ptweit(i_i)%m(i_p - 1, i_s) = dfctrs
                ptpair1(i_i)%m(i_p - 1, i_s) = i_o
                ptpair2(i_i)%m(i_p - 1, i_s) = j_o
             enddo c
          enddo b

          !
          ! Take only those pairs, namely in the range from n_pairs0(i_s, i_i)
          ! to n_pairs (which  may be an empty list,  of course) whose weights
          ! satisfy an empiric criterion:
          !
          ptmin = max(ptlim, (ptweit(i_i)%m(n_pairs, i_s) * 1.0e-1_r8_kind))
          n_pairs0(i_s, i_i) = 1
          do i_p = 1, n_pairs
             if (ptweit(i_i)%m(i_p, i_s) .lt. ptmin) n_pairs0(i_s, i_i) = i_p + 1
          enddo

!!$          !
!!$          ! FIXME: This  looks like  a workaround for  the code that  does not
!!$          ! fully respect the range of valid pairs computed above:
!!$          !
!!$          if (n_pairs0(i_s, i_i) .eq. n_pairs + 1) then
!!$             ptpair1(i_i)%m(n_pairs, i_s) = 0
!!$          endif

        enddo a
      enddo irrep

      ! deallocate dimensions of irreps again
      deallocate(dim_irrep)

      !
      ! This was called in chargefit code right after "call ptpairs()"
      !
      if( output_chargefit ) then
          do i = 1, size(ptweit)
              do j = 1, n_spin
                  write(output_unit,*) 'Paare IRREP ',i,' Spin ',j
                  write(output_unit,*) 'PTPAIR1'
                  write(output_unit,*) (ptpair1(i)%m(k, j), k=1,n_pairs)
                  write(output_unit,*) 'PTPAIR2'
                  write(output_unit,*) (ptpair2(i)%m(k, j), k=1,n_pairs)
                  write(output_unit,*) 'PTWEITS'
                  write(output_unit,*) (ptweit(i)%m(k, j), k=1,n_pairs)
              enddo
          enddo
      endif

999   continue
      !
      ! Slaves dont compute them (historically or for a reason) but
      ! just receive them from the master:
      !
      call pairs_bcast()
    end subroutine ptpairs

  subroutine deallocate_pairs()
    ! Purpose: Free allocated space for pairs stuff
    implicit none
    external error_handler
    !** End of interface *****************************************

    integer :: memstat
    integer(kind=i4_kind) :: i_ir,n_irr

    n_irr = size(ptpair1)

    do i_ir=1,n_irr
       deallocate( ptpair1(i_ir)%m, ptpair2(i_ir)%m,STAT=memstat)
       if(memstat/=0)call error_handler("pm/deallocate_pairs: dealloc 1 failed")

       deallocate( ptweit(i_ir)%m,STAT=memstat)
       if(memstat/=0)call error_handler("pm/deallocate_pairs: dealloc 2 failed")
    enddo

    deallocate(ptpair1, ptpair2, n_pairs0,STAT=memstat)
    if(memstat/=0)call error_handler("pm/deallocate_pairs: dealloc 3 failed")

    deallocate(ptweit,STAT=memstat)
    if(memstat/=0)call error_handler("pm/deallocate_pairs: dealloc 4 failed")
  end subroutine deallocate_pairs

  subroutine pairs_bcast()
    !
    ! Runs in parallel context.
    !
    use comm, only: comm_rank, comm_bcast
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: irrep

    ASSERT(allocated(ptpair1))
    ASSERT(allocated(ptpair2))
    ASSERT(allocated(ptweit))
    ASSERT(size(ptweit)==size(ptpair1))
    ASSERT(size(ptweit)==size(ptpair2))

    do irrep = 1, size(ptweit)
        call comm_bcast(ptpair1(irrep)%m)
        call comm_bcast(ptpair2(irrep)%m)

        ! FIXME: weights? Slaves dont seem to use them:
        call comm_bcast(ptweit(irrep)%m)
    enddo
    call comm_bcast(n_pairs0)
  end subroutine pairs_bcast


!--------------- End of module ----------------------------------
end module pairs_module

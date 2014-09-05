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
MODULE F77_SCALAPACK
  !
  ! Interface declarations of (selected) ScaLAPACK subroutines.
  !
  ! Everybody who needs other than listed here is encouraged to extend
  ! this interface.
  !
  ! Usage:
  !
  !  subroutine sub (...)
  !   use f77_scalapack, only: pdgemm
  !   implicit none
  !
  !   call pdgemm (...)
  !  end subroutine sub
  !
  ! AM, 11/2012

  IMPLICIT NONE

  INTERFACE

     pure integer function indxl2g(indxloc, nb, iproc, isrcproc, nprocs)
       implicit none
       integer, intent(in) :: indxloc, iproc, isrcproc, nb, nprocs
     end function indxl2g

      pure INTEGER FUNCTION NUMROC (N, NB, IPROC, ISRCPROC, NPROCS)
        implicit none
        INTEGER, intent(in) :: IPROC, ISRCPROC, N, NB, NPROCS
      END FUNCTION NUMROC

     integer function blacs_pnum(icontxt, prow, pcol)
       implicit none
       integer, intent(in) :: icontxt, prow, pcol
     end function blacs_pnum

     subroutine blacs_gridmap(icontxt, usermap, ldumap, nprow, npcol)
       implicit none
       integer, intent(inout) :: icontxt
       integer, intent(in) :: usermap(ldumap, npcol)
       integer, intent(in) :: ldumap, nprow, npcol
     end subroutine blacs_gridmap

     subroutine blacs_gridexit (icontxt)
       implicit none
       integer, intent(in) :: icontxt
     end subroutine blacs_gridexit

     subroutine blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)
       implicit none
       integer, intent(in) :: icontxt
       integer, intent(out) :: nprow, npcol, myprow, mypcol
     end subroutine blacs_gridinfo

      subroutine descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info)
        integer :: icsrc, ictxt, irsrc, lld, m, mb, n, nb
        integer, intent(out) :: desc(*)
        integer, intent(out) :: info
      end subroutine descinit

     subroutine pdgemm(transa, transb, m, n, k, &
          alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
       implicit none
       character :: transa, transb
       integer :: m, n, k, ia, ja, ib, jb, ic, jc
       integer :: desca(*), descb(*), descc(*)
       double precision :: alpha, beta
       double precision :: a(*), b(*), c(*)
     end subroutine pdgemm

     subroutine  pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
       implicit none
       integer :: ia, ic, ja, jc, m, n
       double precision :: alpha, beta
       integer :: desca(*), descc(*)
       double precision :: a(*), c(*)
     end subroutine pdtran

     subroutine pdsygvx(ibtype, jobz, range, uplo, n, a, ia, ja, &
          desca, b, ib, jb, descb, vl, vu, il, iu, &
          abstol, m, nz, w, orfac, z, iz, jz, descz, &
          work, lwork, iwork, liwork, ifail, iclustr, &
          gap, info)
       implicit none
       character :: jobz, range, uplo
       integer :: ia, ib, ibtype, il, info, iu, iz, ja, jb, jz, liwork, lwork, m, n, nz
       double precision :: abstol, orfac, vl, vu
       integer :: desca(*), descb(*), descz(*), iclustr(*), ifail(*), iwork(*)
       double precision :: a(*), b(*), gap(*), w(*), work(*), z(*)
       intent(inout) :: a, b
     end subroutine pdsygvx

  END INTERFACE

END MODULE F77_SCALAPACK

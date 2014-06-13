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
module commpack_module
!-----------------------------------------------------------
!--------------------------------------------------------------

!------------ Modules used ------------------------------------
implicit none
save
private

interface commpack

  subroutine mpix_pkdouble_vec( dp, nitem, stride, info) bind(C)
    use iso_c_binding
    real(c_double), dimension(*)  :: dp
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_pkdouble_vec

  subroutine mpix_pkint_vec( np, nitem, stride, info) bind(C)
    use iso_c_binding
    integer(c_int), dimension(*)  :: np
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_pkint_vec

  subroutine mpix_pkdouble_vecsc( dp, nitem, stride, info) bind(C)
    use iso_c_binding
    real(c_double)                :: dp
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_pkdouble_vecsc

  subroutine mpix_pkint_vecsc( np, nitem, stride, info) bind(C)
    use iso_c_binding
    integer(c_int)                :: np
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_pkint_vecsc

  subroutine mpix_pkdouble_scalar( p, info ) bind(C)
    use iso_c_binding
    real(c_double)                :: p
    integer(c_int)                :: info
  end subroutine mpix_pkdouble_scalar

  subroutine mpix_pkint_scalar( p, info ) bind(C)
    use iso_c_binding
    integer(c_int)                :: p
    integer(c_int)                :: info
  end subroutine mpix_pkint_scalar

  module procedure mpix_packstring

  module procedure mpix_packlogical

  module procedure mpix_packlogicalvec

end interface


interface communpack

  subroutine mpix_upkdouble_vec( dp, nitem, stride, info) bind(C)
    use iso_c_binding
    real(c_double), dimension(*)  :: dp
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_upkdouble_vec

  subroutine mpix_upkint_vec( np, nitem, stride, info) bind(C)
    use iso_c_binding
    integer(c_int), dimension(*)  :: np
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_upkint_vec

  subroutine mpix_upkdouble_vecsc( dp, nitem, stride, info) bind(C)
    use iso_c_binding
    real(c_double)                :: dp
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_upkdouble_vecsc

  subroutine mpix_upkint_vecsc( np, nitem, stride, info) bind(C)
    use iso_c_binding
    integer(c_int)                :: np
    integer(c_int)                :: nitem, stride, info
  end subroutine mpix_upkint_vecsc

  subroutine mpix_upkdouble_scalar( p, info ) bind(C)
    use iso_c_binding
    real(c_double)                :: p
    integer(c_int)                :: info
  end subroutine mpix_upkdouble_scalar

  subroutine mpix_upkint_scalar( p, info ) bind(C)
    use iso_c_binding
    integer(c_int)                :: p
    integer(c_int)                :: info
  end subroutine mpix_upkint_scalar

  module procedure mpix_unpackstring

  module procedure mpix_unpacklogical

  module procedure mpix_unpacklogicalvec

end interface

public :: commpack
public :: communpack

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

   !*************************************************************
   subroutine mpix_packlogical(p,info)
   use iso_c_binding
   implicit none
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical,               intent(in)  :: p
   integer(c_int), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer(c_int)              :: i

   if ( p ) then
      i = 1
   else
      i = 0
   endif
   call commpack(i, info)
   end subroutine mpix_packlogical
   !*************************************************************

   !*************************************************************
   subroutine mpix_unpacklogical(p,info)
   use iso_c_binding
   implicit none
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical,               intent(out) :: p
   integer(c_int), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer(c_int)              :: i

   call communpack(i, info)
   if ( info .ne. 0 ) return
   p =  ( i .eq. 1 )
   end subroutine mpix_unpacklogical
   !*************************************************************

   !*************************************************************
   subroutine mpix_packlogicalvec(p,nitem, stride,info)
   use iso_c_binding
   implicit none
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical, dimension(*), intent(in)  :: p
   integer(c_int), intent(out) :: info
   integer(c_int), intent(in)  :: nitem, stride
   !------------ Declaration of local variables -----------------
   integer(c_int)              :: i(nitem)
   integer(c_int)              :: j

   do j=1,nitem
      if ( p(j) ) then
         i(j) = 1
      else
         i(j) = 0
      endif
   enddo
   call commpack(i, nitem, stride, info)
   end subroutine mpix_packlogicalvec
   !*************************************************************

   !*************************************************************
   subroutine mpix_unpacklogicalvec(p,nitem, stride,info)
   use iso_c_binding
   implicit none
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical, dimension(*), intent(out) :: p
   integer(c_int), intent(out) :: info
   integer(c_int), intent(in)  :: nitem, stride
   !------------ Declaration of local variables -----------------
   integer(c_int)              :: i(nitem), j

   call communpack(i, nitem, stride, info)
   if ( info .ne. 0 ) return
   do j=1,nitem
      p(j) =  ( i(j) .eq. 1 )
   enddo
   end subroutine mpix_unpacklogicalvec
   !*************************************************************

   !*************************************************************
   subroutine mpix_packstring(p,info)
   use iso_c_binding
   implicit none
   ! Purpose: packs string by packing first length, then contents
   !------------ Declaration of formal parameters ---------------
   character(len=*), intent(in) :: p
   integer(c_int), intent(out)  :: info
   ! *** end of interface ***

   integer(c_int) :: i, array(len(p))

   do i = 1, len(p)
     array(i) = iachar(p(i:i))
   enddo

   call commpack(len(p), info)
   if ( info .ne. 0 ) return

   call commpack(array, len(p), 1, info)
   end subroutine mpix_packstring
   !*************************************************************

   !*************************************************************
   subroutine mpix_unpackstring(p,info)
   use iso_c_binding
   implicit none
   ! Purpose: unpacks string by unpacking first length, then contents
   !------------ Declaration of formal parameters ---------------
   character(len=*), intent(out) :: p
   integer(c_int), intent(out)   :: info
   ! *** end of interface ***

   integer(c_int) :: i, length, array(len(p))

   call communpack(length, info)
   if ( info .ne. 0 ) return
   if ( len(p) .lt. length ) info = -100
   if ( info .ne. 0 ) return

   call communpack(array, length, 1, info)

   p = repeat(" ", len(p))
   do i = 1, length
     p(i:i) = achar(array(i))
   enddo
   end subroutine mpix_unpackstring
   !*************************************************************

!--------------- End of module ----------------------------------
end module commpack_module


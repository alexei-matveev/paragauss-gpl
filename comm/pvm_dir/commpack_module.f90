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
!===============================================================
! Public interface of module
!===============================================================
module  commpack_module
!-----------------------------------------------------------
!
!  Purpose: Substituting the standard pvmfpack and pvmfunpack
!           subroutines that do not work with calls for 
!           different data types because of more rigorose
!           type checking in f90 than in f77.
!           The Specific c-routines contained in pvmpack.c 
!           are called. Since calling syntax for c-routines 
!           might differ for different compilers and machines
!           this file might need to be replaced in such cases.
!           Version with underscore appended to subroutine
!           (pvmpack_.c) and without (pvmpack.c) already exist.
!           However, machine dependency is kept in one place.
!           commpack.c shoulb be compiled seperately and included
!           when linking.
!           A more portable solution is possible for the price of
!           lower performance by replacing the c-routines by 
!           fortran subroutines that call the standard pvmfpack
!           routine and that are compiled in seperately. Some of
!           the necessary functions already exist.
!
!           The following generic subroutines are defined:
!
!           Packing:
!
!           subroutine commpack( p, nitem, stride, info)
!             <different standard data types>, intention(in) :: p(*)
!             integer(kind=i4_kind) ::  nitem, stride, info
!
!             special form for scalar arguments and strings:
!
!           subroutine commpack( p, info )
!             <different standard data types>, intention(in) :: p
!             integer(kind=i4_kind) ::  info
!
!           Unpacking:
!
!           subroutine communpack( p, nitem, stride, info)
!             <different standard data types>, intention(out) :: p(*)
!             integer(kind=i4_kind) ::  nitem, stride, info
!
!             special form for scalar arguments and strings:
!
!           subroutine communpack( p, info )
!             <different standard data types>, intention(out) :: p
!             integer(kind=i4_kind) ::  info
!
!           Only those generic subroutines should be used !
!
!           Attention for packing strings and logicals !!!
!             Intermediate variables generrated by wrapper
!             routines implemented here are packed. Thus,
!             PvmDataInPlace does not work for these types !!!
!
!================================================================
! End of public interface of module
!================================================================
!
!           For strings, the the following subroutines are defined
!
!           subroutine commpackstring( p, info )
!             character(*), intention(in) :: p
!             integer(kind=i4_kind) ::  info
!
!           subroutine communpackstring( p, info )
!             character(*), intention(out) :: p
!             integer(kind=i4_kind) ::  info
!
!           these subroutines can also be accessed via
!           the generic interface:
!
!           subroutine commpack( p, info )
!             character(*), intention(in) :: p
!             integer(kind=i4_kind) ::  info
!
!           subroutine communpack( p, info )
!             character(*), intention(out) :: p
!             integer(kind=i4_kind) ::  info
!
!
!
!  Module called by: Various subroutines that pack and unpack PVM data
!
!
!  References: PVM-Manual; man pvm_pack; man pvm_unpack
! 
!
!  Author: TB
!  Date: 20.06.95
!
!
!---------------------------------------------------------------
! Modifications
!---------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!--------------------------------------------------------------

!------------ Modules used ------------------------------------
use type_module ! type specification parameters

implicit none

save



!------------ Interface statements ---------------------------

interface commpack

  subroutine pvm_pkbyte_vec( cp, nitem, stride, info)
    use type_module
    character, dimension(*)              :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkbyte_vec

  subroutine pvm_pkcplx_vec( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind), dimension(*)  :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkcplx_vec

  subroutine pvm_pkdcplx_vec( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind), dimension(*) :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkdcplx_vec

  subroutine pvm_pkdouble_vec( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind), dimension(*)     :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkdouble_vec

  subroutine pvm_pkfloat_vec( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind), dimension(*)     :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkfloat_vec

  subroutine pvm_pkint_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkint_vec

  subroutine pvm_pkshort_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkshort_vec

  subroutine pvm_pkbyte_vecsc( cp, nitem, stride, info)
    use type_module
    character                            :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkbyte_vecsc

  subroutine pvm_pkcplx_vecsc( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind)                :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkcplx_vecsc

  subroutine pvm_pkdcplx_vecsc( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind)               :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkdcplx_vecsc

  subroutine pvm_pkdouble_vecsc( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind)                   :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkdouble_vecsc

  subroutine pvm_pkfloat_vecsc( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind)                   :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkfloat_vecsc

  subroutine pvm_pkint_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkint_vecsc

  subroutine pvm_pkshort_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_pkshort_vecsc

  subroutine pvm_pkcplx_scalar( p, info )
    use type_module
    complex(kind=c8_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkcplx_scalar

  subroutine pvm_pkdcplx_scalar( p, info )
    use type_module
    complex(kind=c16_kind)               :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkdcplx_scalar

  subroutine pvm_pkdouble_scalar( p, info )
    use type_module
    real(kind=r8_kind)                   :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkdouble_scalar

  subroutine pvm_pkfloat_scalar( p, info )
    use type_module
    real(kind=r4_kind)                   :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkfloat_scalar

  subroutine pvm_pkint_scalar( p, info )
    use type_module
    integer(kind=i4_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkint_scalar

  subroutine pvm_pkshort_scalar( p, info )
    use type_module
    integer(kind=i2_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_pkshort_scalar

  module procedure commpackstring

  module procedure commpacklogical

  module procedure commpacklogicalvec

end interface


interface communpack

  subroutine pvm_upkbyte_vec( cp, nitem, stride, info)
    use type_module
    character, dimension(*)              :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkbyte_vec

  subroutine pvm_upkcplx_vec( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind), dimension(*)  :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkcplx_vec

  subroutine pvm_upkdcplx_vec( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind), dimension(*) :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkdcplx_vec

  subroutine pvm_upkdouble_vec( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind), dimension(*)     :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkdouble_vec

  subroutine pvm_upkfloat_vec( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind), dimension(*)     :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkfloat_vec

  subroutine pvm_upkint_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkint_vec

  subroutine pvm_upkshort_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkshort_vec

  subroutine pvm_upkbyte_vecsc( cp, nitem, stride, info)
    use type_module
    character                            :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkbyte_vecsc

  subroutine pvm_upkcplx_vecsc( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind)                :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkcplx_vecsc

  subroutine pvm_upkdcplx_vecsc( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind)               :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkdcplx_vecsc

  subroutine pvm_upkdouble_vecsc( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind)                   :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkdouble_vecsc

  subroutine pvm_upkfloat_vecsc( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind)                   :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkfloat_vecsc

  subroutine pvm_upkint_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkint_vecsc

  subroutine pvm_upkshort_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
  end subroutine pvm_upkshort_vecsc

  subroutine pvm_upkcplx_scalar( p, info )
    use type_module
    complex(kind=c8_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkcplx_scalar

  subroutine pvm_upkdcplx_scalar( p, info )
    use type_module
    complex(kind=c16_kind)               :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkdcplx_scalar

  subroutine pvm_upkdouble_scalar( p, info )
    use type_module
    real(kind=r8_kind)                   :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkdouble_scalar

  subroutine pvm_upkfloat_scalar( p, info )
    use type_module
    real(kind=r4_kind)                   :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkfloat_scalar

  subroutine pvm_upkint_scalar( p, info )
    use type_module
    integer(kind=i4_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkint_scalar

  subroutine pvm_upkshort_scalar( p, info )
    use type_module
    integer(kind=i2_kind)                :: p
    integer(kind=i4_kind)                :: info
  end subroutine pvm_upkshort_scalar

  module procedure communpackstring

  module procedure communpacklogical

  module procedure communpacklogicalvec

end interface




!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

   !*************************************************************
   subroutine commpacklogical(p,info)
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical,               intent(in)  :: p
   integer(kind=i4_kind), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind)              :: i
   !------------ External subroutine used -----------------------
   !------------ Executable code --------------------------------
   if ( p ) then
      i = 1
   else
      i = 0
   endif
   call pvm_pkint_scalar(i,info)
   end subroutine commpacklogical
   !*************************************************************

   !*************************************************************
   subroutine communpacklogical(p,info)
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical,               intent(out) :: p
   integer(kind=i4_kind), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind)              :: i
   !------------ External subroutine used -----------------------
   !------------ Executable code --------------------------------
   call pvm_upkint_scalar(i,info)
   if ( info .ne. 0 ) return
   p =  ( i .eq. 1 )
   end subroutine communpacklogical
   !*************************************************************

   !*************************************************************
   subroutine commpacklogicalvec(p,nitem,stride,info)
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical, dimension(*), intent(in)  :: p
   integer(kind=i4_kind), intent(out) :: info
   integer(kind=i4_kind), intent(in)  :: nitem,stride
   !------------ Declaration of local variables -----------------
   integer(kind=i2_kind)              :: i(nitem)
   integer(kind=i4_kind)              :: j
   !------------ External subroutine used -----------------------
   !------------ Executable code --------------------------------
   do j=1,nitem,stride
      if ( p(j) ) then
         i(j) = 1
      else
         i(j) = 0
      endif
   enddo
   call commpack(i,nitem,stride,info)
   end subroutine commpacklogicalvec
   !*************************************************************

   !*************************************************************
   subroutine communpacklogicalvec(p,nitem,stride,info)
   ! Purpose: packs logical as short int
   !------------ Declaration of formal parameters ---------------
   logical, dimension(*), intent(out) :: p
   integer(kind=i4_kind), intent(out) :: info
   integer(kind=i4_kind), intent(in)  :: nitem,stride
   !------------ Declaration of local variables -----------------
   integer(kind=i2_kind)              :: i(nitem)
   integer(kind=i4_kind)              :: j
   !------------ External subroutine used -----------------------
   !------------ Executable code --------------------------------
   call communpack(i,nitem,stride,info)
   if ( info .ne. 0 ) return
   do j=1,nitem,stride
      p(j) =  ( i(j) .eq. 1 )
   enddo
   end subroutine communpacklogicalvec
   !*************************************************************

   !*************************************************************
   subroutine commpackstring(p,info)
   ! Purpose: packs string by packing first length, then contents
   !------------ Declaration of formal parameters ---------------
   character*(*),         intent(in)  :: p
   integer(kind=i4_kind), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer                              :: length
   !------------ External subroutine used -----------------------
   !external pvm_pkstring
   !------------ Executable code --------------------------------
   length=len(p)
   call commpack(length,info)
   if ( info .ne. 0 ) return
   call pvm_pkstring(p,length,info)
   end subroutine commpackstring
   !*************************************************************

   !*************************************************************
   subroutine communpackstring(p,info)
   ! Purpose: unpacks string by unpacking first length, then contents
   !------------ Declaration of formal parameters ---------------
   character*(*),         intent(out) :: p
   integer(kind=i4_kind), intent(out) :: info
   !------------ Declaration of local variables -----------------
   integer                              :: length
   !------------ External subroutine used -----------------------
   !external pvm_upkstring
   !------------ Executable code --------------------------------
   call communpack(length,info)
   if ( len(p) .lt. length ) info = -100_i4_kind
   if ( info .ne. 0 ) return
   p=repeat(" ",len(p))
   call pvm_upkstring(p,length,info)
   end subroutine communpackstring
   !*************************************************************


!--------------- End of module ----------------------------------
end module commpack_module


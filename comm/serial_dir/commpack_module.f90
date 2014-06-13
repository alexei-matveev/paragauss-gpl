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
#include <def.h>
!===============================================================
! Public interface of module
!===============================================================
module  commpack_module
!-----------------------------------------------------------
  !  Version for serial work
  !
  !  This module contains dummy versions of COMMPACK and COMMUNPACK
  !  subroutines. This module is only needed for correct compilation
  !  of Paragaus`s serial version.
  !================================================================
  ! End of public interface of module
  !================================================================
  !
  !  Author: AS
  !  Date: 7/98
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

  module procedure pvm_pkbyte_vec

  module procedure pvm_pkcplx_vec

  module procedure pvm_pkdcplx_vec

  module procedure pvm_pkdouble_vec

  module procedure pvm_pkfloat_vec

  module procedure pvm_pkint_vec

  module procedure pvm_pkshort_vec

  module procedure pvm_pkbyte_vecsc

  module procedure pvm_pkcplx_vecsc

  module procedure pvm_pkdcplx_vecsc

  module procedure pvm_pkdouble_vecsc

  module procedure pvm_pkfloat_vecsc

  module procedure pvm_pkint_vecsc

  module procedure pvm_pkshort_vecsc

  module procedure pvm_pkcplx_scalar

  module procedure pvm_pkdcplx_scalar

  module procedure pvm_pkdouble_scalar

  module procedure pvm_pkfloat_scalar

  module procedure pvm_pkint_scalar

  module procedure pvm_pkshort_scalar

  module procedure commpackstring

  module procedure commpacklogical

  module procedure commpacklogicalvec

end interface


interface communpack

  module procedure pvm_upkbyte_vec

  module procedure pvm_upkcplx_vec

  module procedure pvm_upkdcplx_vec

  module procedure pvm_upkdouble_vec

  module procedure pvm_upkfloat_vec

  module procedure pvm_upkint_vec

  module procedure pvm_upkshort_vec

  module procedure pvm_upkbyte_vecsc

  module procedure pvm_upkcplx_vecsc

  module procedure pvm_upkdcplx_vecsc

  module procedure pvm_upkdouble_vecsc

  module procedure pvm_upkfloat_vecsc

  module procedure pvm_upkint_vecsc

  module procedure pvm_upkshort_vecsc

  module procedure pvm_upkcplx_scalar

  module procedure pvm_upkdcplx_scalar

  module procedure pvm_upkdouble_scalar

  module procedure pvm_upkfloat_scalar

  module procedure pvm_upkint_scalar

  module procedure pvm_upkshort_scalar

  module procedure communpackstring

  module procedure communpacklogical

  module procedure communpacklogicalvec

end interface




!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains
  subroutine pvm_pkbyte_vec( cp, nitem, stride, info)
    use type_module
    character, dimension(*)              :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkbyte_vec

  subroutine pvm_pkcplx_vec( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind), dimension(*)  :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkcplx_vec

  subroutine pvm_pkdcplx_vec( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind), dimension(*) :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdcplx_vec

  subroutine pvm_pkdouble_vec( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind), dimension(*)     :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdouble_vec

  subroutine pvm_pkfloat_vec( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind), dimension(*)     :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkfloat_vec

  subroutine pvm_pkint_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkint_vec

  subroutine pvm_pkshort_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkshort_vec

  subroutine pvm_pkbyte_vecsc( cp, nitem, stride, info)
    use type_module
    character                            :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkbyte_vecsc

  subroutine pvm_pkcplx_vecsc( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind)                :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkcplx_vecsc

  subroutine pvm_pkdcplx_vecsc( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind)               :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdcplx_vecsc

  subroutine pvm_pkdouble_vecsc( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind)                   :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdouble_vecsc

  subroutine pvm_pkfloat_vecsc( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind)                   :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkfloat_vecsc

  subroutine pvm_pkint_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkint_vecsc

  subroutine pvm_pkshort_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkshort_vecsc

  subroutine pvm_pkcplx_scalar( p, info )
    use type_module
    complex(kind=c8_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkcplx_scalar

  subroutine pvm_pkdcplx_scalar( p, info )
    use type_module
    complex(kind=c16_kind)               :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdcplx_scalar

  subroutine pvm_pkdouble_scalar( p, info )
    use type_module
    real(kind=r8_kind)                   :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkdouble_scalar

  subroutine pvm_pkfloat_scalar( p, info )
    use type_module
    real(kind=r4_kind)                   :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkfloat_scalar

  subroutine pvm_pkint_scalar( p, info )
    use type_module
    integer(kind=i4_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkint_scalar

  subroutine pvm_pkshort_scalar( p, info )
    use type_module
    integer(kind=i2_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_pkshort_scalar




  subroutine pvm_upkbyte_vec( cp, nitem, stride, info)
    use type_module
    character, dimension(*)              :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkbyte_vec

  subroutine pvm_upkcplx_vec( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind), dimension(*)  :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkcplx_vec

  subroutine pvm_upkdcplx_vec( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind), dimension(*) :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdcplx_vec

  subroutine pvm_upkdouble_vec( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind), dimension(*)     :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdouble_vec

  subroutine pvm_upkfloat_vec( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind), dimension(*)     :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkfloat_vec

  subroutine pvm_upkint_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkint_vec

  subroutine pvm_upkshort_vec( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind), dimension(*)  :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkshort_vec

  subroutine pvm_upkbyte_vecsc( cp, nitem, stride, info)
    use type_module
    character                            :: cp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkbyte_vecsc

  subroutine pvm_upkcplx_vecsc( xp, nitem, stride, info)
    use type_module
    complex(kind=c8_kind)                :: xp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkcplx_vecsc

  subroutine pvm_upkdcplx_vecsc( zp, nitem, stride, info)
    use type_module
    complex(kind=c16_kind)               :: zp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdcplx_vecsc

  subroutine pvm_upkdouble_vecsc( dp, nitem, stride, info)
    use type_module
    real(kind=r8_kind)                   :: dp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdouble_vecsc

  subroutine pvm_upkfloat_vecsc( fp, nitem, stride, info)
    use type_module
    real(kind=r4_kind)                   :: fp
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkfloat_vecsc

  subroutine pvm_upkint_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i4_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkint_vecsc

  subroutine pvm_upkshort_vecsc( np, nitem, stride, info)
    use type_module
    integer(kind=i2_kind)                :: np
    integer(kind=i4_kind)                :: nitem, stride, info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkshort_vecsc

  subroutine pvm_upkcplx_scalar( p, info )
    use type_module
    complex(kind=c8_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkcplx_scalar

  subroutine pvm_upkdcplx_scalar( p, info )
    use type_module
    complex(kind=c16_kind)               :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdcplx_scalar

  subroutine pvm_upkdouble_scalar( p, info )
    use type_module
    real(kind=r8_kind)                   :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkdouble_scalar

  subroutine pvm_upkfloat_scalar( p, info )
    use type_module
    real(kind=r4_kind)                   :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkfloat_scalar

  subroutine pvm_upkint_scalar( p, info )
    use type_module
    integer(kind=i4_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkint_scalar

  subroutine pvm_upkshort_scalar( p, info )
    use type_module
    integer(kind=i2_kind)                :: p
    integer(kind=i4_kind)                :: info
    info = 0
    ABORT('serial')
  end subroutine pvm_upkshort_scalar

  subroutine commpacklogical(p,info)
    use type_module
   logical,               intent(in)  :: p
   integer(kind=i4_kind), intent(out) :: info
    info = 0
    ABORT('serial')
   end subroutine commpacklogical

   subroutine communpacklogical(p,info)
     use type_module
     logical,               intent(out) :: p
     integer(kind=i4_kind), intent(out) :: info
    info = 0
    ABORT('serial')
   end subroutine communpacklogical

   subroutine commpacklogicalvec(p,nitem,stride,info)
     use type_module
     logical, dimension(*), intent(in)  :: p
     integer(kind=i4_kind), intent(out) :: info
     integer(kind=i4_kind), intent(in)  :: nitem,stride
    info = 0
    ABORT('serial')
   end subroutine commpacklogicalvec

  subroutine communpacklogicalvec(p,nitem,stride,info)
    use type_module
    logical, dimension(*), intent(out) :: p
    integer(kind=i4_kind), intent(out) :: info
    integer(kind=i4_kind), intent(in)  :: nitem,stride
    info = 0
    ABORT('serial')
  end subroutine communpacklogicalvec

   subroutine commpackstring(p,info)
     use type_module
     character*(*),         intent(in)  :: p
     integer(kind=i4_kind), intent(out) :: info
    info = 0
    ABORT('serial')
   end subroutine commpackstring

   subroutine communpackstring(p,info)
     use type_module
     character*(*)                      :: p
     integer(kind=i4_kind), intent(out) :: info
    p = ''
    info = 0
    ABORT('serial')
   end subroutine communpackstring


!--------------- End of module ----------------------------------
end module commpack_module


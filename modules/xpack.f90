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
module xpack
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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

  interface pck
     module procedure pck_int
     module procedure pck_int1D
     module procedure pck_int2D
     module procedure pck_int3D
     !----------------------
     module procedure pck_double
     module procedure pck_double1D
     module procedure pck_double2D
     module procedure pck_double3D
     module procedure pck_double4D
     module procedure pck_double_arr
     !----------------------
     module procedure pck_logic
     module procedure pck_logic1D
     module procedure pck_logic2D
     module procedure pck_logic3D
     !----------------------
     module procedure pck_string
  end interface

  interface upck
     module procedure upck_int
     module procedure upck_int1D
     module procedure upck_int2D
     module procedure upck_int3D
     !----------------------
     module procedure upck_double
     module procedure upck_double1D
     module procedure upck_double2D
     module procedure upck_double3D
     module procedure upck_double4D
     module procedure upck_double_arr
     !----------------------
     module procedure upck_logic
     module procedure upck_logic1D
     module procedure upck_logic2D
     module procedure upck_logic3D
     !----------------------
     module procedure upck_string
  end interface

!!$  interface pck_arr ! private
!!$     module procedure pck_double_arr
!!$  end interface
!!$
!!$  interface upck_arr ! private
!!$     module procedure upck_double_arr
!!$  end interface

  !------------ public functions and subroutines ------------------

  public pck, upck

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine pck_int(i)
    use error_module, only: error
    use comm_module, only: commpack
    implicit none
    integer(IK),intent(in) :: i
    ! *** end of interface ***

    integer(IK) :: info

    call commpack(i,info)
    if(info/=0)call error("xp/pck_int: pck failed")

!!$    print  *,MyID,'pck_int:',i
  end subroutine pck_int

  subroutine upck_int(i)
    use error_module, only: error
    use comm_module, only: communpack
    implicit none
    integer(IK),intent(out) :: i
    ! *** end of interface ***

    integer(IK) :: info

    call communpack(i,info)
    if(info/=0)call error("xp/upck_int: upck failed")

!!$    print  *,MyID,'upck_int:',i
  end subroutine upck_int

  subroutine pck_int_arr(v,n)
    use error_module, only: error
    use comm_module, only: commpack
    implicit none
    integer(IK),intent(in) :: v( * )
    integer(IK),intent(in) :: n
    ! *** end of interface ***

    integer(IK) :: info

    call commpack(v,n,1,info)
    if(info/=0)call error("xp/pck_int_arr: pck failed")
  end subroutine pck_int_arr

  subroutine upck_int_arr(v,n)
    use error_module, only: error
    use comm_module, only: communpack
    implicit none
    integer(IK),intent(out) :: v( * )
    integer(IK),intent(in)  :: n
    ! *** end of interface ***

    integer(IK) :: info

    call communpack(v,n,1,info)
    if(info/=0)call error("xp/upck_int_arr: pck failed")
  end subroutine upck_int_arr

  subroutine pck_int1D(v)
    implicit none
    integer(IK),intent(in) :: v(:)
    ! *** end of interface ***

    call pck_int_arr(v,size(v))
  end subroutine pck_int1D

  subroutine upck_int1D(v)
    implicit none
    integer(IK),intent(out) :: v(:)
    ! *** end of interface ***

    call upck_int_arr(v,size(v))
  end subroutine upck_int1D

  subroutine pck_int2D(v)
    implicit none
    integer(IK),intent(in) :: v(:,:)
    ! *** end of interface ***

    call pck_int_arr(v,size(v))
  end subroutine pck_int2D

  subroutine upck_int2D(v)
    implicit none
    integer(IK),intent(out) :: v(:,:)
    ! *** end of interface ***

    call upck_int_arr(v,size(v))
  end subroutine upck_int2D

  subroutine pck_int3D(v)
    implicit none
    integer(IK),intent(in) :: v(:,:,:)
    ! *** end of interface ***

    call pck_int_arr(v,size(v))
  end subroutine pck_int3D

  subroutine upck_int3D(v)
    implicit none
    integer(IK),intent(out) :: v(:,:,:)
    ! *** end of interface ***

    call upck_int_arr(v,size(v))
  end subroutine upck_int3D

  subroutine pck_double(v)
    use error_module, only: error
    use comm_module, only: commpack
    implicit none
    real(RK),intent(in) :: v
    ! *** end of interface ***

    integer(IK) :: info

    call commpack(v,info)
    if(info/=0)call error("xp/pck_double: pck failed")
  end subroutine pck_double

  subroutine upck_double(v)
    use error_module, only: error
    use comm_module, only: communpack
    implicit none
    real(RK),intent(out) :: v
    ! *** end of interface ***

    integer(IK) :: info

    call communpack(v,info)
    if(info/=0)call error("xp/pck_double: upck failed")
  end subroutine upck_double

  subroutine pck_double_arr(v,n)
    use error_module, only: error
    use comm_module, only: commpack
    implicit none
    real(RK),intent(in)    :: v( * )
    integer(IK),intent(in) :: n
    ! *** end of interface ***

    integer(IK) :: info

    call commpack(v,n,1,info)
    if(info/=0)call error("xp/pck_double_arr: pck failed")
  end subroutine pck_double_arr

  subroutine upck_double_arr(v,n)
    use error_module, only: error
    use comm_module, only: communpack
    implicit none
    real(RK),intent(out)   :: v( * )
    integer(IK),intent(in) :: n
    ! *** end of interface ***

    integer(IK) :: info

    call communpack(v,n,1,info)
    if(info/=0)call error("xp/upck_double_arr: upck failed")
  end subroutine upck_double_arr

  subroutine pck_double1D(v)
    implicit none
    real(RK),intent(in) :: v(:)
    ! *** end of interface ***

    call pck_double_arr(v,size(v))
  end subroutine pck_double1D

  subroutine upck_double1D(v)
    implicit none
    real(RK),intent(out) :: v(:)
    ! *** end of interface ***

    call upck_double_arr(v,size(v))
  end subroutine upck_double1D

  subroutine pck_double2D(v)
    implicit none
    real(RK),intent(in) :: v(:,:)
    ! *** end of interface ***

    call pck_double_arr(v,size(v))
  end subroutine pck_double2D

  subroutine upck_double2D(v)
    implicit none
    real(RK),intent(out) :: v(:,:)
    ! *** end of interface ***

    call upck_double_arr(v,size(v))
  end subroutine upck_double2D

  subroutine pck_double3D(v)
    implicit none
    real(RK),intent(in) :: v(:,:,:)
    ! *** end of interface ***

    call pck_double_arr(v,size(v))
  end subroutine pck_double3D

   subroutine upck_double3D(v)
    implicit none
    real(RK),intent(out) :: v(:,:,:)
    ! *** end of interface ***

    call upck_double_arr(v,size(v))
  end subroutine upck_double3D

  subroutine pck_double4D(v)
    implicit none
    real(RK),intent(in) :: v(:,:,:,:)
    ! *** end of interface ***

    call pck_double_arr(v,size(v))
  end subroutine pck_double4D

   subroutine upck_double4D(v)
    implicit none
    real(RK),intent(out) :: v(:,:,:,:)
    ! *** end of interface ***

    call upck_double_arr(v,size(v))
  end subroutine upck_double4D

  subroutine pck_logic(v)
    implicit none
    logical,intent(in) :: v
    ! *** end of interface ***

    integer(IK) :: iv

    iv = merge(1,0,v)

    call pck(iv)
  end subroutine pck_logic

  subroutine upck_logic(v)
    implicit none
    logical,intent(out) :: v
    ! *** end of interface ***

    integer(IK) :: iv

    call upck(iv)

    v = (iv/=0)
  end subroutine upck_logic

  subroutine pck_logic1D(v)
    implicit none
    logical,intent(in) :: v(:)
    ! *** end of interface ***

    integer(IK),dimension(size(v)) :: ibuf

    ibuf = merge(1,0,v)

    call pck(ibuf)
  end subroutine pck_logic1D

  subroutine upck_logic1D(v)
    implicit none
    logical,intent(out) :: v(:)
    ! *** end of interface ***

    integer(IK),dimension(size(v)) :: ibuf

    call upck(ibuf)

    v = (ibuf/=0)
  end subroutine upck_logic1D

  subroutine pck_logic2D(m)
    implicit none
    logical,intent(in) :: m(:,:)
    ! *** end of interface ***

    logical,    dimension(size(m)) ::  buf
    integer(IK),dimension(size(m)) :: ibuf

    buf  = reshape(m,(/size(m)/))

    ibuf = merge(1,0,buf)

    call pck(ibuf)
  end subroutine pck_logic2D

  subroutine upck_logic2D(m)
    implicit none
    logical,intent(out) :: m(:,:)
    ! *** end of interface ***

    logical,    dimension(size(m)) ::  buf
    integer(IK),dimension(size(m)) :: ibuf

    call upck(ibuf)

    buf = (ibuf/=0)

    m = reshape(buf,shape(m))
  end subroutine upck_logic2D

  subroutine pck_logic3D(m)
    implicit none
    logical,intent(in) :: m(:,:,:)
    ! *** end of interface ***

    logical,    dimension(size(m)) ::  buf
    integer(IK),dimension(size(m)) :: ibuf

    buf  = reshape(m,(/size(m)/))

    ibuf = merge(1,0,buf)

    call pck(ibuf)
  end subroutine pck_logic3D

  subroutine upck_logic3D(m)
    implicit none
    logical,intent(out) :: m(:,:,:)
    ! *** end of interface ***

    logical,    dimension(size(m)) ::  buf
    integer(IK),dimension(size(m)) :: ibuf

    call upck(ibuf)

    buf = (ibuf/=0)

    m = reshape(buf,shape(m))
  end subroutine upck_logic3D

  subroutine pck_string(s)
    use error_module
    use comm_module, only: commpack
    implicit none
    character(len=*),intent(in) :: s
    ! *** end of interface ***

    integer(IK) :: info

    call commpack(s,info)
    if(info/=0)call error("xp/pck_string: pck failed")
  end subroutine pck_string

  subroutine upck_string(s)
    use error_module
    use comm_module, only: communpack
    implicit none
    character(len=*),intent(inout) :: s
    ! *** end of interface ***

    integer(IK) :: info

    call communpack(s,info)
    if(info/=0)call error("xp/upck_string: upck failed")
  end subroutine upck_string

  !--------------- End of module ----------------------------------
end module xpack

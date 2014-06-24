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
module debug
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
       IK=>i4_kind, &
       RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------
!!$  type, public ::  debug_
!!$  end type debug_

  !------------ Declaration of constants and variables ------------
!!$  integer(kind=IK), parameter, public  :: debug_
!!$  real(kind=RK),    parameter, public  :: debug_
!!$  logical,               parameter, public  :: debug_
!!$  character,             parameter, public  :: debug_
!!$  integer(kind=IK),            public  :: debug_
!!$  real(kind=RK),               public  :: debug_
!!$  logical,                          public  :: debug_
!!$  character,                        public  :: debug_


  !------------ Interface statements ------------------------------
  interface show
     module procedure show_real_matrix
     ! show_real_matrix(name,real(:,:))
     module procedure show_vec
     ! show_vec(int,vec(int))
     module procedure show_2_vec
     ! show_2_vec(int,vec1(int),vec2(int))
     module procedure show_3_vec
     ! show_3_vec(int,vec1(int),vec2(int),vec3(int))
     module procedure show_4_vec
     module procedure show_5_vec
     module procedure show_cmplx
     ! show_cmplx(int,M_real(int,int),M_imag(int,int))
     module procedure show_cmtrx
     ! show_cmtrx(M)  type(cmatrix):: M
  end interface

  interface zhow
     module procedure zhow_real_matrix
  end interface

  interface octave
     module procedure octave_real_matrix
     module procedure octave_real_vector
     module procedure octave_real_3D
     module procedure octave_real_4D
  end interface

  interface show_diff
     module procedure show_diff_real_1D
  end interface

#ifdef WITH_ISNAN
  interface isNaN
     module procedure isNaN_s
     module procedure isNaN_1D
  end interface

  interface countNaN
     module procedure count_NaN_scalar
     module procedure count_NaN_1D
     module procedure count_NaN_2D
     module procedure count_NaN_3D
     module procedure count_NaN_4D
     module procedure count_NaN_5D
     module procedure count_NaN_buf
  end interface

  public :: isNaN
  public :: countNaN
#endif

#ifdef WITH_ISINF
  interface countInf
     module procedure count_Inf_scalar
     module procedure count_Inf_1D
     module procedure count_Inf_2D
     module procedure count_Inf_3D
     module procedure count_Inf_4D
     module procedure count_Inf_5D
     module procedure count_Inf_buf
  end interface

  interface countNInf
     module procedure count_NInf_scalar
     module procedure count_NInf_1D
     module procedure count_NInf_2D
     module procedure count_NInf_3D
     module procedure count_NInf_4D
     module procedure count_NInf_5D
  end interface

  !------------ public functions and subroutines ------------------
  public :: isInf
  public :: countInf
#endif

#if defined WITH_ISNAN && defined WITH_ISINF
  public :: countNInf ! countNaN + countInf
#endif

#ifdef FPP_GFORTRAN_BUGS /* uses GETPID and SLEEP intrinsics avail with -std=gnu */
  public :: print_pids_and_sleep!(sec)
#endif
  public :: NaN
  public :: show
  public :: show_diff_real_buf
  public :: octave
  public :: disp
  public :: display_oct_sqr_mat

  interface dump
     module procedure dump_real_buf
     module procedure dump_real_1D
     module procedure dump_real_2D
     module procedure dump_real_3D
     module procedure dump_int_buf
     module procedure dump_int
  end interface

  public :: dump

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------
!!$  type 
!!$  end type 

  !------------ Declaration of constants and variables ----

  real(RK), parameter :: zero = 0.0_RK 

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

#if 0 /* seems to be unused */
  subroutine wait()
    use error_module, only:MyID
    implicit none
    ! *** end of interface ***

   integer :: rate,cur,old,tim
   integer :: long
   integer :: ticks,tick, tikk
   real    :: slice = 0.0

   call SYSTEM_CLOCK(old,rate)
!!$   print *,MyID,'wait: entered, old=',old,'rate=',rate

   long  = 1*rate
   ticks = 3
   tim = 0
   tikk = 0
   do while ( tim < long ) ! 1 sec
      call SYSTEM_CLOCK(cur)
      tim = cur - old
      tick = (ticks*tim)/long
      if( tick /= tikk )then
         print *,MyID,'wait  ...',tick
         tikk = tick
      end if
   end do
  end subroutine wait
#endif

#ifdef FPP_GFORTRAN_BUGS /* uses GETPID and SLEEP intrinsics avail with -std=gnu */
  subroutine print_pids_and_sleep(sec)
    !
    ! Print PIDs and sleep for some time to allow
    ! the developer to attach the debugger to any
    ! of the processors, e.g. by
    !
    !   gdb -p PID
    !
    use comm, only: comm_rank, comm_size, comm_reduce ! be sure not to induce circular deps!
    implicit none
    integer(ik), intent(in) :: sec ! number of seconds to sleep
    ! *** end of interface ***

    integer(ik)              :: p, np
    integer(ik), allocatable :: pids(:)

    ! print PIDs of all processors e.g. to attach a debugger:
    np = comm_size()

    allocate(pids(np))

    p = comm_rank() + 1

    pids = 0
    pids(p) = getpid()
    call comm_reduce(pids)

    if( p == 1 )then
      do p=1,np
        print *,'Process ',p,' has PID',pids(p)
      enddo
      print *,'going to sleep for', sec,' seconds ...'
    endif
    deallocate(pids)

    call sleep(sec)
    print *,'Procsess ',p,' awake!'
  end subroutine print_pids_and_sleep
#endif

  function NaN()
    implicit none
    real(8) :: NaN
    ! *** end of interface ***

    real(8) :: z

    z = -1.0_8
    NaN = sqrt(z)
  end function NaN

#ifdef WITH_ISNAN
  function count_NaN_scalar(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v
    integer(ik)             :: NaNs
    ! *** end of interface ***

    integer(ik) :: cnt

    cnt=0
    if(isNaN(v)) cnt=cnt+1
    NaNs = cnt
  end function count_NaN_scalar

  function count_NaN_1D(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v(:)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    NaNs = count_NaN_buf(size(v),v)
  end function count_NaN_1D

  function count_NaN_2D(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v(:,:)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    NaNs = count_NaN_buf(size(v),v)
  end function count_NaN_2D

  function count_NaN_3D(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    NaNs = count_NaN_buf(size(v),v)
  end function count_NaN_3D

  function count_NaN_4D(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    NaNs = count_NaN_buf(size(v),v)
  end function count_NaN_4D

  function count_NaN_5D(v) result(NaNs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:,:)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    NaNs = count_NaN_buf(size(v),v)
  end function count_NaN_5D

  function count_NaN_buf(vl,v) result(NaNs)
    implicit none
    integer(ik), intent(in) :: vl
    real(rk), intent(in)    :: v(*)
    integer(ik)             :: NaNs
    ! *** end of interface ***

    integer(ik) :: i,cnt

    cnt=0
    do i=1,vl
       if(isNaN(v(i))) cnt=cnt+1
    enddo
    NaNs = cnt
  end function count_NaN_buf

  function isNaN_s(v) result(yes)
    implicit none
    real(rk), intent(in)    :: v
    logical                 :: yes
    ! *** end of interface ***

#ifndef INTRINSIC_ISNAN
    logical, external :: isNaN ! external func, provided by Intel
#else
    intrinsic :: isNaN ! intrinsic func, provided by Gfortran
#endif

    yes = isNaN(v)
  end function isNaN_s

  function isNaN_1D(v) result(mask)
    implicit none
    real(rk), intent(in)    :: v(:)
    logical                 :: mask(size(v))
    ! *** end of interface ***

    integer(ik) :: i

    do i=1,size(v)
       mask(i) = isNaN(v(i))
    enddo
  end function isNaN_1D
#endif

#ifdef WITH_ISINF
  function count_Inf_scalar(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v
    integer(ik)             :: Infs
    ! *** end of interface ***

    integer(ik) :: cnt

    cnt=0
    if(isInf(v)) cnt=cnt+1
    Infs = cnt
  end function count_Inf_scalar

  function count_Inf_1D(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v(:)
    integer(ik)             :: Infs
    ! *** end of interface ***

    Infs = count_Inf_buf(size(v),v)
  end function count_Inf_1D

  function count_Inf_2D(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v(:,:)
    integer(ik)             :: Infs
    ! *** end of interface ***

    Infs = count_Inf_buf(size(v),v)
  end function count_Inf_2D

  function count_Inf_3D(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:)
    integer(ik)             :: Infs
    ! *** end of interface ***

    Infs = count_Inf_buf(size(v),v)
  end function count_Inf_3D

  function count_Inf_4D(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:)
    integer(ik)             :: Infs
    ! *** end of interface ***

    Infs = count_Inf_buf(size(v),v)
  end function count_Inf_4D

  function count_Inf_5D(v) result(Infs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:,:)
    integer(ik)             :: Infs
    ! *** end of interface ***

    Infs = count_Inf_buf(size(v),v)
  end function count_Inf_5D

  function count_Inf_buf(vl,v) result(Infs)
    implicit none
    integer(ik), intent(in) :: vl
    real(rk), intent(in)    :: v(*)
    integer(ik)             :: Infs
    ! *** end of interface ***

    integer(ik) :: i,cnt

    cnt=0
    do i=1,vl
       if(isInf(v(i))) cnt=cnt+1
    enddo
    Infs = cnt
  end function count_Inf_buf

  function isInf(v) result(yes)
    implicit none
    real(rk), intent(in)    :: v
    logical                 :: yes
    ! *** end of interface ***

    yes = abs(v) > huge(v)
  end function isInf
#endif

#if defined WITH_ISNAN && defined WITH_ISINF
  function count_NInf_scalar(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_scalar

  function count_NInf_1D(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v(:)
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_1D

  function count_NInf_2D(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v(:,:)
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_2D

  function count_NInf_3D(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:)
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_3D

  function count_NInf_4D(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:)
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_4D

  function count_NInf_5D(v) result(NInfs)
    implicit none
    real(rk), intent(in)    :: v(:,:,:,:,:)
    integer(ik)             :: NInfs
    ! *** end of interface ***

    NInfs = countNaN(v) + countInf(v)
  end function count_NInf_5D
#endif

  !*************************************************************
  subroutine show_real_matrix(name,MM)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)       :: name
    real(RK),dimension(:,:),intent(in) :: MM
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(IK),parameter       :: zmax=5
    integer(IK)           :: i,j,n,m,zn,zm
    integer(IK)           ::         dn,dm
    !------------ Executable code --------------------------------

    n  = size(MM,1)
    m  = size(MM,2)
    zn = min(n-1,zmax)
    zm = min(m-1,zmax)

    ! for 1 x 1 matrices, do not divide by zero:
    dn = max(zn,1)
    dm = max(zm,1)

    ! 1) the (part of) the matrix:
    write(*,'(A)') '===================='

    print *, '<<< ',trim(name),' (',n,'x',m,') >>>'
    write(*,'(10I21)') (1+j*(m-1)/dm, j=0,zm)
    if( size(MM) > 0 )then
    do i=0,zn
       write(*,'(I3," ",10(1PE20.12," "))')&
            1+i*(n-1)/dn, (MM(1+i*(n-1)/dn, 1+j*(m-1)/dm), j=0,zm)
    enddo
    endif

    ! 2) traces:
    if( m == n )then
      print *, 'Trace of the matrix =',trace(MM)
      if( m < 100 ) &
      print *, 'Trace of its square =',trace(matmul(MM,MM))
    endif

    ! 3) show the pattern of non-zero matrix elements:
    call zhow(MM,cut=1.0E-3_rk)

    ! 4) dump the content in matlab/Octave format:
!   call octave(name,MM)
    contains

      function trace(M) result(tr)
        implicit none
        real(RK), intent(in) :: M(:,:)
        real(RK)             :: tr
        ! *** end of interface ***

        integer(IK) :: i

        ASSERT(size(M,1)==size(M,2))

        tr = 0.0
        do i=1,size(M,1)
          tr = tr + M(i,i)
        enddo
      end function trace

  end subroutine show_real_matrix
  !*************************************************************
  
  !*************************************************************
  subroutine zhow_real_matrix(mat,cut)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(RK),dimension(:,:),intent(in) :: mat
    real(RK),optional      ,intent(in) :: cut
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(IK)           :: i,j,n,m
    integer(IK),parameter :: many=128
    integer(IK)           :: d_n,d_m
    real(RK)              :: cut_
    !------------ Executable code --------------------------------

    if( size(mat) == 0 ) RETURN

    n = size(mat,1)
    m = size(mat,2)

    d_n = 1
    d_m = 1
    if(n>many) d_n = max(n/many,2)
    if(m>many) d_m = max(m/many,2)

    cut_ = 1.0E-7_rk
    if(present(cut)) cut_ = cut

    print *,'Pattern (',n,'x',m,') of non-zero elements with cutoff ',cut_

    if( m <= 99 )then
      if(m>9 ) write(*,'(4X  ,100I2)') (mod(i/10 ,10), i=1,m)
               write(*,'(4X  ,100I2)') (mod(i    ,10), i=1,m)
      do i=1,n
               write(*,'(I3,X,100A2)') i, ( zchar(mat(i,j)), j=1,m)
      enddo
    else
               write(*,'(4X  ,200I1)') (    i/100    , i=1,m,d_m)
               write(*,'(4X  ,200I1)') (mod(i/10 ,10), i=1,m,d_m)
               write(*,'(4X  ,200I1)') (mod(i    ,10), i=1,m,d_m)
      do i=1,n,d_n
               write(*,'(I3,X,200A1)') i, ( zchar(mat(i,j)), j=1,m,d_m)
      enddo
    endif
  contains
    function zchar(x) result(c)
      implicit none
      real(RK), intent(in) :: x
      character(len=1)     :: c ! result
      ! *** end of interface ***

      if ( abs(x) < cut_ ) then
        c = ' '
      else
        if( x > 0.0_rk )then
          c = '+'
        else
          c = '-'
        endif
      endif
#ifdef WITH_ISNAN
      if (isNaN (x)) then
        c = 'N'
      endif
#endif
    end function zchar
  end subroutine zhow_real_matrix
  !*************************************************************

  !*************************************************************
  subroutine octave_real_matrix(name,MM)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)       :: name
    real(RK),dimension(:,:),intent(in) :: MM
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(IK)           :: i,j,n,m
    !------------ Executable code --------------------------------

    n = size(MM,1)
    m = size(MM,2)

    write(*,'(A)') '<octave>'
    write(*,*)     '# ',trim(var(name)),' (',n,'x',m,')'
    write(*,'((A)," = [")',ADVANCE='yes') trim(var(name))
    if( size(MM) > 0 )then
    do i=1,n
    do j=1,m
       write(*,'(1PE24.16E3,A2)',ADVANCE='no')  &
          MM(i,j), term(j,m)
       if( mod(j,5)==0 ) write(*,'()',ADVANCE='yes')
    enddo
       if( mod(m,5)/=0 ) write(*,'()',ADVANCE='yes')
    enddo
    endif
    write(*,'("];")',ADVANCE='yes')
    !rite(*,*)     '#',trim(var(name))//'_sum    = ',sum(MM)     , ';'
    !rite(*,*)     '#',trim(var(name))//'_sumabs = ',sum(abs(MM)), ';'
    write(*,'(A)') '</octave>'
  end subroutine octave_real_matrix

  subroutine octave_real_vector(name,MM)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)       :: name
    real(RK),dimension(:),intent(in)   :: MM
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(IK)           :: i,n
    !------------ Executable code --------------------------------

    n = size(MM)

    write(*,'(A)') '<octave>'
    write(*,*)     '# ',trim(var(name)),' (',n,')'
    write(*,'((A)," = [")',ADVANCE='yes') trim(var(name))
    write(*,'(5(1PE24.16E3,"; "))')  ( MM(i), i=1,n )
    write(*,'("];")',ADVANCE='yes')
    write(*,'(A)') '</octave>'
  end subroutine octave_real_vector

  subroutine octave_real_3D(name,MM)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)       :: name
    real(RK),dimension(:,:,:),intent(in)   :: MM
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    write(*,'(A)') '<octave>'
    write(*,*)     '# ',trim(var(name)),' (',shape(MM),') in fortran order'
    write(*,'((A)," = [")',ADVANCE='yes') trim(var(name))
    write(*,'(5(1PE24.16E3,"; "))')  MM
    write(*,'("];")',ADVANCE='yes')
    write(*,'(A)') '</octave>'
  end subroutine octave_real_3D

  subroutine octave_real_4D(name,MM)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)       :: name
    real(RK),dimension(:,:,:,:),intent(in)   :: MM
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    write(*,'(A)') '<octave>'
    write(*,*)     '# ',trim(var(name)),' (',shape(MM),') in fortran order'
    write(*,'((A)," = [")',ADVANCE='yes') trim(var(name))
    write(*,'(5(1PE24.16E3,"; "))')  MM
    write(*,'("];")',ADVANCE='yes')
    write(*,'(A)') '</octave>'
  end subroutine octave_real_4D

  function term(i,n) result(comma)
    integer(IK), intent(in) :: i,n
    character(len=2)        :: comma
    ! *** end of interface ***

    comma = ', '
    if( mod(i,5)==0 ) comma = ',\'
    if(     i   ==n ) comma = '; '
  end function term

  function var(name) result(str)
    character(len=*)         :: name
    character(len=LEN(name)) :: str
    ! *** end of interface ***

    integer(IK) :: i

    str = name
    do i=1,LEN(trim(str))
      if( str(i:i) == ' ' ) str(i:i) = '_'
    enddo
  end function var
  !*************************************************************

  subroutine disp(name,F)
    implicit none
    character(len=*), intent(in) :: name
    real(RK)        , intent(in) :: F(:,:)
    ! *** end of interface ***

    integer(IK) :: ii,i
    integer(IK)            :: n
    integer(IK), parameter :: nel=6
    integer(IK)            :: js(nel)

    n = size(F,1)
    ASSERT(n==size(F,2))

    print *,'disp: '//name//', matrix size=',n,'x',n

    js( 1)     =  1
    js( 2)     =  2
    js( 3)     =  3
    js( 4)     =  n-2
    js( 5)     =  n-1
    js( 6)     =  n

    where(js>n) js = n

    print    '(4X,(7I14),/: ("... ",4X,7I14))', js
    do ii=1,nel
       i = js(ii)
       print '(I4,(7E14.6),/: ("... ",4X,7E14.6))', i,&
            F(i,js)
    enddo
  end subroutine disp

  !*************************************************************
  subroutine show_diff_real_buf(name,buf1,buf2,siz)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: name
    real(RK), intent(in)         :: buf1(*), buf2(*)
    integer(IK), intent(in)      :: siz
    !** End of interface *****************************************

    call show_diff_real_1D(name,buf1(1:siz),buf2(1:siz))
  end subroutine show_diff_real_buf
  !*************************************************************

  !*************************************************************
  subroutine show_diff_real_1D(name,buf1,buf2)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: name
    real(RK), intent(in)         :: buf1(:), buf2(:)
    !** End of interface *****************************************

    integer(IK), parameter :: N_PREC=3
    real(RK), parameter :: &
         eps(N_PREC) = (/ 1.0E-4_rk, 1.0E-6_rk, 1.0E-8_rk /)

    integer(IK) :: arrloc(1),loc,n,i
    real(RK)    :: rms1,rms2,rmsdiff,maxdiff,val1,val2,abs1,abs2

    n = size(buf1)
    ASSERT(n==size(buf2))

    rms1    = sqrt(sum(buf1**2)/n)
    rms2    = sqrt(sum(buf2**2)/n)
    rmsdiff = sqrt(sum((buf1-buf2)**2)/n)

    print *,name,' rms(1)   =',rms1
    print *,name,' rms(2)   =',rms2

    arrloc = maxloc(abs(buf1-buf2))
    loc = arrloc(1)
    val1 = buf1(loc)
    val2 = buf2(loc)
    maxdiff = val1 - val2
    abs1 = abs(val1)
    abs2 = abs(val2)

    if( rms1+rms2 /= zero )then
       print *,name,' rms(1-2) =',rmsdiff, 2*rmsdiff/(rms1+rms2),'*100%'
    else
       print *,name,' rms(1-2) =',rmsdiff, '(zero arrays?)'
    endif

    if( abs1+abs2 /= zero )then
       print *,name,' max(1-2) =',maxdiff, 2*maxdiff/(abs1+abs2),'*100% (',val1,val2,')'
    else
       print *,name,' max(1-2) =',maxdiff, '(crossing zero?) (',val1,val2,')'
    endif

    do i=1,N_PREC
       print *,name,' { threshold ',eps(i),'}'
       if(rms1/=zero)then
          print *,name,' count(1)=',count(abs(buf1/rms1)>eps(i)),'/',n
       endif
       if(rms2/=zero)then
          print *,name,' count(2)=',count(abs(buf2/rms2)>eps(i)),'/',n
       endif
       if( rms1+rms2 /= zero )then
          print *,name,' count(1-2)=',count(2*abs(buf1-buf2)/(rms1+rms2)>eps(i)),'/',n
       endif
    enddo
  end subroutine show_diff_real_1D
  !*************************************************************

  subroutine show_vec(name,v)
    character(len=*), intent(in) :: name
    real(RK)        ,intent(in)  :: v(:)
    ! *** end of interface ***

    integer(IK) :: i,n,nn

    n = size(v)
    nn=merge(n,24,n<24)

    write(*,'(A)') '===================='

    print *, '<<< ',trim(name),' (',n,') >>>'
    do i=1,nn
       write(*,'(I3," ",F25.10)') i,v(i)
    enddo
    if(nn<n)write(*,'(A)') '... the rest skipped ...'
    call octave(name,v)
  end subroutine show_vec

  subroutine show_2_vec(name,v1,v2)
    character(len=*), intent(in) :: name
    real(RK),intent(in)          :: v1(:),v2(:)

    integer :: i,n,nn

    n = size(v1)
    nn=merge(n,24,n<24)

    write(*,'(A)') '===================='

    print *, '<<< ',trim(name),' (',n,') >>>'
    do i=1,nn
       write(*,'(I3," ",F25.10," ",F25.10)') i,v1(i),v2(i)
    enddo
    if(nn<n)write(*,'(A)') '... the rest skipped ...'
  end subroutine show_2_vec

  subroutine show_3_vec(n,v1,v2,v3)
    integer(IK),intent(in) :: n
    real(RK),intent(in)    :: v1(n),v2(n),v3(n)

    integer :: i,nn

    nn=merge(n,24,n<24)

    do i=1,nn
       write(*,'(I3," ",F25.10," ",F25.10," ",F25.10)') i,v1(i),v2(i),v3(i)
    enddo
    if(nn<n)write(*,'(A)') '... the rest skipped ...'
  end subroutine show_3_vec

  subroutine show_4_vec(n,v1,v2,v3,v4)
    integer(IK),intent(in) :: n
    real(RK),intent(in)    ::&
         & v1(n),v2(n),v3(n),v4(n)

    integer :: i,nn

    nn=merge(n,24,n<24)

    do i=1,nn
       write(*,'(I3," ",F18.5," ",F18.5," ",F18.5," ",F18.5)')&
            & i,v1(i),v2(i),v3(i),v4(i)
    enddo
    if(nn<n)write(*,'(A)') '... the rest skipped ...'
  end subroutine show_4_vec

  subroutine show_5_vec(n,v1,v2,v3,v4,v5)
    integer(IK),intent(in) :: n
    real(RK),intent(in)    ::&
         & v1(n),v2(n),v3(n),v4(n),v5(n)

    integer :: i,nn

    nn=merge(n,24,n<24)

    do i=1,nn
       write(*,'(I3," ",F12.5," ",F12.5," ",F12.5," ",F12.5," ",F12.5," ")')&
            & i,v1(i),v2(i),v3(i),v4(i),v5(i)
    enddo
    if(nn<n)write(*,'(A)') '... the rest skipped ...'
  end subroutine show_5_vec

  subroutine show_cmtrx(M)
    use matrix_types, only: cmatrix
    implicit none
    type(cmatrix),intent(in) :: M
    ! *** end of interface ***

    call show(M%re,M%im)
  end subroutine show_cmtrx

  subroutine show_cmplx(M_real,M_imag)
    use error_module
    real(RK),dimension(:,:),intent(in) :: M_real,M_imag
    ! *** end of interface ***

    integer(IK),parameter       :: many=7
    real(RK),dimension(4) :: norms
    integer(IK)           :: i,j,n,m,d_n,d_m
    
    call error(any(shape(M_real)/=shape(M_imag)),"isam/show_cmplx: shapes ???")

    n = size(M_real,1)
    m = size(M_real,2)

    d_n = 1
    d_m = 1
    if(n>many) d_n = max(n/many,2)
    if(m>many) d_m = max(m/many,2)

    DPRINT'isam/show_cmplx: shape=',n,'x',m

    ! print a block>>>
    write(*,'(A)') '\real>>>'
    write(*,'(8I16)') (i,i=1,m,d_m)
    do i=1,n,d_n
       write(*,'(I3," ",8(E15.8," "))')&
            & i, (M_real(i,j), j=1,m,d_m)
    enddo
    
    write(*,'(A)') '\imag>>>'
    write(*,'(8I16)') (i,i=1,m,d_m)
    do i=1,n,d_n
       write(*,'(I3," ",8(E15.8," "))')&
            & i, (M_imag(i,j), j=1,m,d_m)
    enddo

    if(n.eq.m)then
       norms = examine_cmplx(M_real,M_imag)
       write(*,'(A)') 'Summary:'
       write(*,'("real part: diag= ",E20.10," offdiag= ",E20.10 )')&
            & norms(1:2)
       write(*,'("imag part: diag= ",E20.10," offdiag= ",E20.10 )')&
            & norms(3:4)
    endif
  end subroutine show_cmplx

  subroutine display_oct_sqr_mat(mat_name, A)
    !-------------------------------------------------------------------------------------
    !  Purpose: To display a square matrix in octave format. Matrix will be deisplayed
    !           by a width of 5 columns.
    !
    !  Author: RR
    !  Date: February 2009
    !
    !------------ Modules used -----------------------------------------------------------
    use type_module
    implicit none

    !------------ Declaration of formal parameters ---------------------------------------
    character(len = *)           :: mat_name
    real(r8_kind), intent(in)    :: A(:,:)

    !******* End of interface ************************************************************

    !------------ Declaration of local variables -----------------------------------------
    integer(i4_kind)             :: i, j, l, memstat
    integer(i4_kind)             :: N, M, M_mod_nc
    real(r8_kind), allocatable   :: temp_A(:)
    !------------ Executable code --------------------------------------------------------
    !------------ Size of the matrix -----------------------------------------------------
    M = size(A, 1)
    if (size(A, 1) .ne. size(A, 2)) then
      stop "ERROR: In subroutine display_sqr_mat  argument 'A' is not a square matrix"
    endif

    if (M .gt. 999) write(*, '(/,4x, a)') "Warning(DISPLAY): Size the matrix exceeds 999 "

    allocate( temp_A(M**2), stat = memstat )
    if(memstat /= 0) stop "ERROR: Allocation of 'temp_A' failed in display_module"

    write(*,'(4x, a)') '=============================='

    temp_A = reshape(transpose(A), (/M**2/))

    write (*, '( 4x, 3a, i3, a, i3, a )')"<<< ", trim(mat_name)," (",M,"x",M,") >>>"

    N = M / 5
    !******* N is the number of blocks of M x 5 matrices to be printed *******************

    M_mod_nc = mod(M,5)
    !**** The mod(M,5) of the matrix decides the specific formatting of the text displayed
select case(M_mod_nc)

      case(0)
      !******* If the number of columns is an integer multiple of 5 **********************

        do l = 0, N-1
          write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*l+1," through ",5*l+5
          do i = 0, M-1
            write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*l+1, M*i + 5*l+5)
          enddo
        enddo

      case(1)
      !******* If the number of columns is an integer multiple of 5 plus 1 ***************

        do l = 0, N-1
          write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*l+1," through ",5*l+5
          do i = 0, M-1
            write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*l+1, M*i + 5*l+5)
          enddo
        enddo

        write (*, '(/,4x, a, i4 /)')"Column ",5*N+1
        do i = 0, M-1
          write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*N+1, M*i + 5*N+mod(M,5))
        enddo

      case(2)
      !******* If the number of columns is an integer multiple of 5 plus 2 ***************

        do l = 0, N-1
          write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*l+1," through ",5*l+5
          do i = 0, M-1
            write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*l+1, M*i + 5*l+5)
          enddo
        enddo

        write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*N+1," and ",5*N+ mod(M,5)
        do i = 0, M-1
          write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*N+1, M*i + 5*N+mod(M,5))
        enddo

      case default
      !******* If the number of columns is an integer multiple of 5 plus 3 or 4 *******

        do l = 0, N-1
          write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*l+1," through ",5*l+5
          do i = 0, M-1
            write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*l+1, M*i + 5*l+5)
          enddo
        enddo

        write (*, '(/,4x, a, i4, a, i4,/)')"Columns ",5*N+1," through ",5*N+ mod(M,5)
        do i = 0, M-1
          write (*, '(4x, 5e17.7)') (temp_A(j), j = M*i + 5*N+1, M*i + 5*N+mod(M,5))
        enddo

    end select

    deallocate( temp_A, stat = memstat )
    if(memstat /= 0) stop "ERROR: Deallocation of 'temp_A' failed in display_module"

    !--------------- End of subroutine ---------------------------------------------------
    return

  end subroutine display_oct_sqr_mat

  !***************************************************************************************
  
  function examine_cmplx(M_real,M_imag) result(res)
    use error_module
    implicit none
    real(RK),dimension(:,:),intent(in) :: M_real,M_imag
    real(RK),dimension(4) :: res  ! <<<result
    ! *** end of interface ***

    integer(IK) :: n

    real(RK) ::&
         & diag_norm_real, diag_norm_imag,&
         & offdiag_norm_real, offdiag_norm_imag

    call error(any(shape(M_real)/=shape(M_imag)),"isam/examine_cmplx: shapes ???")
    call error(size(M_real,1)/=size(M_real,2),"isam/examine_cmplx: non square !!!")

    n = size(M_real,1)

    diag_norm_real =   sum(pack(M_real,diag_mask(n))**2)
    diag_norm_imag =   sum(pack(M_imag,diag_mask(n))**2)
    offdiag_norm_real= sum(pack(M_real,offdiag_mask(n))**2)
    offdiag_norm_imag= sum(pack(M_imag,offdiag_mask(n))**2)
    
    res=(/ diag_norm_real,offdiag_norm_real,&
         & diag_norm_imag,offdiag_norm_imag /)
    return
  end function examine_cmplx

  function diag_mask(n) result(m)
    integer(IK),intent(in) :: n
    logical,dimension(n,n) :: m  !<<<result
    ! *** end of interface ***
    
    integer(IK) :: i

    m = .false.
    do i=1,n
       m(i,i) = .true.
    enddo
  end function diag_mask
  
  function offdiag_mask(n) result(m)
    integer(IK),intent(in) :: n
    logical,dimension(n,n) :: m  !<<<result
    ! *** end of interface ***

    integer(IK) :: i

    m = .true.
    do i=1,n
       m(i,i) = .false.
    enddo
  end function offdiag_mask

  !
  ! Dump the data to unit, here in text form.
  ! Save with enough info on structure of the arrays.
  !

  subroutine dump_int_buf(iou, buf, siz)
    implicit none
    integer(IK), intent(in) :: iou
    integer(IK), intent(in) :: buf(*)
    integer(IK), intent(in) :: siz
    ! *** end of interface ***

    write(iou,'(12I10)') buf(:siz)
  end subroutine dump_int_buf

  subroutine dump_int(iou, buf)
    implicit none
    integer(IK), intent(in) :: iou
    integer(IK), intent(in) :: buf
    ! *** end of interface ***

    call dump_int_buf(iou, (/ buf /), 1)
  end subroutine dump_int

  subroutine dump_real_buf(iou, buf, siz)
    implicit none
    integer(IK), intent(in) :: iou
    real(RK)   , intent(in) :: buf(*)
    integer(IK), intent(in) :: siz
    ! *** end of interface ***

    write(iou,'(5E25.16)') buf(:siz)
  end subroutine dump_real_buf

  subroutine dump_real_1D(iou, buf)
    implicit none
    integer(IK), intent(in) :: iou
    real(RK)   , intent(in) :: buf(:)
    ! *** end of interface ***

    call dump_int_buf(iou, shape(buf), 1)
    call dump_real_buf(iou, buf, size(buf))
  end subroutine dump_real_1D

  subroutine dump_real_2D(iou, buf)
    implicit none
    integer(IK), intent(in) :: iou
    real(RK)   , intent(in) :: buf(:,:)
    ! *** end of interface ***

    call dump_int_buf(iou, shape(buf), 2)
    call dump_real_buf(iou, buf, size(buf))
  end subroutine dump_real_2D

  subroutine dump_real_3D(iou, buf)
    implicit none
    integer(IK), intent(in) :: iou
    real(RK)   , intent(in) :: buf(:,:,:)
    ! *** end of interface ***

    call dump_int_buf(iou, shape(buf), 3)
    call dump_real_buf(iou, buf, size(buf))
  end subroutine dump_real_3D

  !--------------- End of module ----------------------------------
end module debug

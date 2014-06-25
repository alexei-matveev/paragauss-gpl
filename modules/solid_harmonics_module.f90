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
module solid_harmonics_module
  !-----------------------------------------------------------
  !
  !  Purpose: Calculation of solid harmonics.
  !           Defines solid_harmonics_type that holds all
  !           solid harmonics up to lmax for a vector of
  !           vector_length gridpoints. Implements subroutines
  !           for memory allocation and freeing and calculation
  !           of solid harmonics.
  !
  !           This module also contains a routine 
  !           'solid_harmonics_calc', which is needed to evaluate
  !           the solid harmonics on a single point in contrast
  !           to 'solid_harmonics_calculate' which calculates
  !           them on a grid. Please do not mistake them
  !
  !  Subroutine called by: orbital_module,ss_calculate
  !
  !  References: ...
  ! 
  !
  !  Author: TB
  !  Date: 26.06.95
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !---------------------------------------------------------------
  ! Modifications
  !---------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   6/96
  ! Description: routine 'solid_harmonics_calc' added
  !              In this routine l and m are treated as one
  !              metaindex
  !              routine 'solid_harmonics_scalar' added
  !              Can be used for only one argument, and not for a 
  !              vector
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   10/96
  ! Description: gridpoints have been transposed in 
  !              solid_harmonics_calculate
  !
  ! Modification (Please copy before editing)
  ! Author: 
  ! Date:  
  ! Description: 
  !--------------------------------------------------------------

# include "def.h"
  use type_module
  USE_MEMLOG

  implicit none
  private
  save
  !== Interrupt end of public interface of module =================


  type solid_harmonics_ltype
     real(kind=r8_kind), pointer    :: m(:,:) ! (vector_length,2*l+1)
  end type solid_harmonics_ltype

  type solid_harmonics_grads_ltype
     real(kind=r8_kind), pointer    :: m(:,:,:) ! (vector_length,3,2*l+1)
     ! secound dimension: x,y,z
  end type solid_harmonics_grads_ltype

  type solid_harmonics_sec_der_ltype
     real(kind=r8_kind), pointer    :: m(:,:,:) ! (vector_length,6,2*l+1)
     ! the thirs dimension is (a symmetric) double derivative:
     ! xx(1)
     ! xy(2) yy(4)
     ! xz(3) yz(5) zz(6)
  end type solid_harmonics_sec_der_ltype

  type solid_harmonics_type
     integer                            :: vector_length
     integer                            :: lmax
     type(solid_harmonics_ltype), pointer   :: l(:) ! (0:lmax)
     real(kind=r8_kind), pointer    :: r2(:) ! (vector_length)  square of distance center to point
  end type solid_harmonics_type

  type solid_harmonics_grads_type
     integer                            :: vector_length
     integer                            :: lmax
     type(solid_harmonics_grads_ltype), pointer   :: l(:) ! (0:lmax)
  end type solid_harmonics_grads_type

  type solid_harmonics_sec_der_type
     integer                            :: vector_length
     integer                            :: lmax
     type(solid_harmonics_sec_der_ltype), pointer   :: l(:) ! (0:lmax)
  end type solid_harmonics_sec_der_type


  !------------ public functions and subroutines ------------------
  public :: solid_harmonics_ltype
  public :: solid_harmonics_type
  public :: solid_harmonics_grads_ltype
  public :: solid_harmonics_grads_type
  public :: solid_harmonics_sec_der_ltype
  public :: solid_harmonics_sec_der_type

  public solid_harmonics_setup, solid_harmonics_shutdown, &
       solid_harmonics_allocate, solid_harmonics_free, &
       solid_harmonics_calculate, &
       solid_harmonics_calculate_grads, &
       solid_harmonics_calc_sec_der, &
       solid_harmonics_calc, solid_harmonics_scalar

  public :: solid_harmonics_calc_3rd_der


  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of private constants and variables ----
  real(kind=r8_kind), allocatable, dimension(:), private :: a,b,x,y,z
  ! a(vector_length), b(vector_length) intermediates
  integer, private :: vec_len = 0
  logical, private :: allocated = .false.

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains 

  !******************************************************************
  subroutine solid_harmonics_setup( vector_length )
    ! allocates memory for private variables
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer, intent(in) :: vector_length
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer status
    !-----------  executable part ----------------------------------
    if ( vector_length .ne. vec_len .and. allocated ) &
         call solid_harmonics_shutdown()
    if ( .not. allocated ) then
       vec_len = vector_length
       allocate( a(vector_length), b(vector_length), x(vector_length), &
            y(vector_length), z(vector_length), stat=status )
       if ( status .ne. 0 ) &
            call error_handler("solid_harmonics_setup: allocate failed")
       MEMLOG(5*vector_length)
       allocated = .true.
    endif
  end subroutine solid_harmonics_setup


  subroutine solid_harmonics_shutdown()
    ! allocates memory for private variables
    !** End of interface ***************************************
    implicit none
    !-----------  private variables --------------------------------
    integer status
    !-----------  executable part ----------------------------------
    if ( allocated ) then
       MEMLOG(-5*vec_len)
       deallocate( a, b, x, y, z, stat=status )
       if ( status .ne. 0 ) &
            call error_handler("solid_harmonics_shutdown: deallocate failed")
       allocated = .false.
       vec_len = 0
    endif
  end subroutine solid_harmonics_shutdown
  !******************************************************************


  subroutine solid_harmonics_allocate( vector_length, lmax, sh, shg, shs, sht )
    ! allocates memory
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer, intent(in) :: vector_length, lmax
    type(solid_harmonics_type), intent (inout), optional :: sh
    type(solid_harmonics_grads_type), intent (inout), optional :: shg
    type(solid_harmonics_sec_der_type), intent (inout), optional :: shs
    type(solid_harmonics_sec_der_type), intent (inout), optional :: sht
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer I_l, status, statussum
    !-----------  executable part ----------------------------------
    if ( present(sh) ) then
       sh%vector_length = vector_length
       sh%lmax = lmax
       allocate( sh%l(0:lmax), stat=status )
       ASSERT(status.eq.0)
       statussum = 0
       do I_l = 0, lmax
          allocate( sh%l(I_l)%m(vector_length,2*I_l+1), stat=status )
          statussum = statussum + status
          MEMLOG(vector_length*(2*I_l+1))
       enddo
       allocate( sh%r2(vector_length), stat=status )
       MEMLOG(vector_length)
       statussum = statussum + status
       ASSERT(statussum.eq.0)
    endif
    if ( present(shg) ) then
       shg%vector_length = vector_length
       shg%lmax = lmax
       allocate( shg%l(0:lmax), stat=status )
       ASSERT(status.eq.0)
       statussum = 0
       do I_l = 0, lmax
          allocate( shg%l(I_l)%m(vector_length,3,2*I_l+1), stat=status )
          statussum = statussum + status
          MEMLOG(vector_length*3*(2*I_l+1))
       enddo
       ASSERT(statussum.eq.0)
    endif

    if ( present(shs) ) then
       shs%vector_length = vector_length
       shs%lmax = lmax
       allocate( shs%l(0:lmax), stat=status )
       ASSERT(status.eq.0)
       statussum = 0
       do I_l = 0, lmax
          allocate( shs%l(I_l)%m(vector_length,6,2*I_l+1), stat=status )
          statussum = statussum + status
          MEMLOG(vector_length*6*(2*I_l+1))
       enddo
       ASSERT(statussum.eq.0)
    endif

    if ( present(sht) ) then
       sht%vector_length = vector_length
       sht%lmax = lmax
       allocate( sht%l(0:lmax), stat=status )
       ASSERT(status.eq.0)
       statussum = 0
       do I_l = 0, lmax
          allocate( sht%l(I_l)%m(vector_length,10,2*I_l+1), stat=status )
          statussum = statussum + status
          MEMLOG(vector_length*10*(2*I_l+1))
       enddo
       ASSERT(statussum.eq.0)
    endif

  end subroutine solid_harmonics_allocate



  subroutine solid_harmonics_free( sh, shg, shs, sht )
    ! frees memory
    implicit none
    !------------ Declaration of formal parameters -----------------
    type(solid_harmonics_type), intent (inout), optional :: sh
    type(solid_harmonics_grads_type), intent (inout), optional :: shg
    type(solid_harmonics_sec_der_type), intent (inout), optional :: shs
    type(solid_harmonics_sec_der_type), intent (inout), optional :: sht
    !** End of interface ***************************************
    !------------ Declaration of subroutines -----------------------
    intrinsic  associated, size
    !-----------  private variables --------------------------------
    integer I_l, status, statussum
    !-----------  executable part ----------------------------------
    if ( present(sh) ) then
       statussum = 0
       if ( associated(sh%l) ) then
          do I_l = 0, size(sh%l)-1
             if ( associated( sh%l(I_l)%m ) ) then
                MEMLOG(-size(sh%l(I_l)%m))
                deallocate( sh%l(I_l)%m, stat=status )
                statussum = statussum + status
             endif
          enddo
          deallocate( sh%l, stat=status )
          statussum = statussum + status
       endif
       if ( associated(sh%r2) ) then
          MEMLOG(-size(sh%r2))
          deallocate( sh%r2, stat=status )
          statussum = statussum + status
       endif
       ASSERT(statussum.eq.0)
       sh%lmax = 0
    endif
    if ( present(shg) ) then
       statussum = 0
       if ( associated(shg%l) ) then
          do I_l = 0, size(shg%l)-1
             if ( associated( shg%l(I_l)%m ) ) then
                MEMLOG(-size(shg%l(I_l)%m))
                deallocate( shg%l(I_l)%m, stat=status )
                statussum = statussum + status
             endif
          enddo
          deallocate( shg%l, stat=status )
          statussum = statussum + status
       endif
       ASSERT(statussum.eq.0)
       shg%lmax = 0
    endif

    if ( present(shs) ) then
       statussum = 0
       if ( associated(shs%l) ) then
          do I_l = 0, size(shs%l)-1
             if ( associated( shs%l(I_l)%m ) ) then
                MEMLOG(-size(shs%l(I_l)%m))
                deallocate( shs%l(I_l)%m, stat=status )
                statussum = statussum + status
             endif
          enddo
          deallocate( shs%l, stat=status )
          statussum = statussum + status
       endif
       ASSERT(statussum.eq.0)
       shs%lmax = 0
    endif

    if ( present(sht) ) then
       statussum = 0
       if ( associated(sht%l) ) then
          do I_l = 0, size(sht%l)-1
             if ( associated( sht%l(I_l)%m ) ) then
                MEMLOG(-size(sht%l(I_l)%m))
                deallocate( sht%l(I_l)%m, stat=status )
                statussum = statussum + status
             endif
          enddo
          deallocate( sht%l, stat=status )
          statussum = statussum + status
       endif
       ASSERT(statussum.eq.0)
       shs%lmax = 0
    endif

  end subroutine solid_harmonics_free



  subroutine solid_harmonics_calculate( sh, grid_points, center_point, N_points )
    ! calculates solid harmonics
    implicit none
    !------------ Declaration of formal parameters -----------------
    type(solid_harmonics_type), intent (inout) :: sh
    real(kind=r8_kind), intent(in) :: grid_points(:,:), &
         center_point(1:3)
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer :: I_l, I_vec, gridlength
    integer :: I_m         ! complex magnetical quantum number
    integer :: m, m1       ! real magnetical quantum number
    integer :: ll1, ll2    ! intermediates calculated from I_l
    real(kind=r8_kind)  :: lr, l1r, ll1r   ! intermediates calculated from I_l
    real(kind=r8_kind)  :: rlm, rl1m, sim, rll    ! intermediates
    real(kind=r8_kind), pointer, dimension(:,:)  :: mp, mp1, mp2
    !-----------  executable part ----------------------------------

    if ( present(N_points) ) then
       gridlength = N_points
    else
       gridlength = sh%vector_length
    endif


    ! Calculate sh for l=0 and r2 and intermediates x, y, z
    mp => sh%l(0)%m
    x(1:gridlength) = grid_points(1:gridlength,1) - center_point(1)
    y(1:gridlength) = grid_points(1:gridlength,2) - center_point(2)
    z(1:gridlength) = grid_points(1:gridlength,3) - center_point(3)
    sh%r2(1:gridlength) = x(1:gridlength) * x(1:gridlength) &
         + y(1:gridlength) * y(1:gridlength) +&
         z(1:gridlength) * z(1:gridlength)
    mp(1:gridlength,1) = 1.0_r8_kind

    if ( sh%lmax == 0 ) return

    ! Calculate sh for l=1
    mp => sh%l(1)%m
    do  I_vec = 1, gridlength
       mp(I_vec,1) = z(I_vec)
       mp(I_vec,2) = x(I_vec)
       mp(I_vec,3) = y(I_vec)
    enddo

    if ( sh%lmax == 1 ) return


    do I_l = 2, sh%lmax
       mp => sh%l(I_l)%m
       mp1 => sh%l(I_l-1)%m
       mp2 => sh%l(I_l-2)%m

       ! intermediates
       lr = real( I_l, r8_kind )
       l1r = real( I_l-1, r8_kind )
       ll1r = real( 2*I_l-1, r8_kind )

       ! m = 1
       do  I_vec = 1, gridlength
          mp(I_vec,1) = ( ll1r * z(I_vec) * mp1(I_vec,1) - &
               l1r * sh%r2(I_vec) * mp2(I_vec,1) ) / lr
       enddo

       ! m = 2, .. , 2 * I_l - 3
       do I_m = 1, I_l - 2

          ! intermediates
          sim = real( I_m * I_m , r8_kind )
          rlm = sqrt( lr * lr - sim )
          rlm = 1 / rlm
          rl1m = sqrt( l1r * l1r - sim )
          m = 2 * I_m
          m1 = m + 1
          do  I_vec = 1, gridlength
             a(I_vec) = ll1r * z(I_vec) * rlm
             b(I_vec) = rl1m * sh%r2(I_vec) * rlm
             mp(I_vec,m)  = a(I_vec) * mp1(I_vec,m) - b(I_vec) * mp2(I_vec,m) 
             mp(I_vec,m1) = a(I_vec) * mp1(I_vec,m1) - b(I_vec) * mp2(I_vec,m1)
          enddo

       enddo

       ! m = 2 * I_l - 2,  2 * I_l - 1
       I_m = I_l - 1
       sim = real( I_m * I_m , r8_kind )
       rlm = sqrt( lr * lr - sim )
       rlm = 1 / rlm
       m = 2 * I_m
       m1 = m + 1
       do  I_vec = 1, gridlength
          a(I_vec) = ll1r * z(I_vec) * rlm
          mp(I_vec,m)  = a(I_vec) * mp1(I_vec,m)
          mp(I_vec,m1) = a(I_vec) * mp1(I_vec,m1)
       enddo

       ! m = 2 * I_l,  2 * I_l + 1
       rll = sqrt ( ll1r / ( 2 * lr ) )
       m = 2 * I_l
       m1 = m + 1
       ll1 =  2 * I_l - 1
       ll2 =  2 * I_l - 2
       do  I_vec = 1, gridlength
          mp(I_vec,m) =  ( x(I_vec) * mp1(I_vec,ll2) - &
               y(I_vec) * mp1(I_vec,ll1) ) * rll
          mp(I_vec,m1) = ( y(I_vec) * mp1(I_vec,ll2) + &
               x(I_vec) * mp1(I_vec,ll1) ) * rll
       enddo

    enddo

  end subroutine solid_harmonics_calculate
  !******************************************************************



  !******************************************************************
  subroutine solid_harmonics_calculate_grads( sh, shg, N_points )
    ! calculates garadients of solid harmonics
    use solhrules_module
    implicit none
    !------------ Declaration of formal parameters -----------------
    type(solid_harmonics_type), intent (in) :: sh
    type(solid_harmonics_grads_type), intent (inout) :: shg
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer :: gridlength,lm,l,m,lm2,l2,m2,i_sum,i_vec
    type(solhrules_differential_type), pointer    :: shrd
    real(kind=r8_kind), pointer, dimension(:,:)   :: sh_m
    real(kind=r8_kind), pointer, dimension(:,:,:) :: shg_m
    real(kind=r8_kind)                            :: coef
    !-----------  executable part ----------------------------------

    if ( present(N_points) ) then
       gridlength = N_points
    else
       gridlength = sh%vector_length
    endif

    do l = 0,sh%lmax
       shg%l(l)%m = 0.0_r8_kind
    enddo

    ! Note how the solid harmonics for l=1 look like:
    ! sh(l=1,m=1)=z  sh(l=1,m=2)=x  sh(l=1,m=3)=y
    ! Thus, we can write the nable operator as the solid
    ! harmonic differential operator for l=1, only
    ! exchanging x,z,z

    lm2=0
    do l2 =  0,sh%lmax
       do m2 = 1, 2*l2+1
          lm2 = lm2+1
          shg_m => shg%l(l2)%m

          shrd => solhrules_differential(3,lm2) ! differential rule for x
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             sh_m => sh%l(l)%m
             coef = shrd%coef(i_sum)
             do i_vec = 1, gridlength
                shg_m(i_vec,1,m2) = shg_m(i_vec,1,m2) + coef * sh_m(i_vec,m)
             enddo
          enddo

          shrd => solhrules_differential(4,lm2) ! differential rule for y
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             sh_m => sh%l(l)%m
             coef = shrd%coef(i_sum)
             do i_vec = 1, gridlength
                shg_m(i_vec,2,m2) = shg_m(i_vec,2,m2) + coef * sh_m(i_vec,m)
             enddo
          enddo

          shrd => solhrules_differential(2,lm2) ! differential rule for z
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             sh_m => sh%l(l)%m
             coef = shrd%coef(i_sum)
             do i_vec = 1, gridlength
                shg_m(i_vec,3,m2) = shg_m(i_vec,3,m2) + coef * sh_m(i_vec,m)
             enddo
          enddo

       enddo
    enddo


  end subroutine solid_harmonics_calculate_grads


  subroutine solid_harmonics_calc_sec_der( shg, shs, N_points )
    ! calculates secound derivatives of solid harmonics
    use solhrules_module
    implicit none
    !------------ Declaration of formal parameters -----------------
    type(solid_harmonics_grads_type), intent (in) :: shg
    type(solid_harmonics_sec_der_type), intent (inout) :: shs
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer :: gridlength,lm,l,m,lm2,l2,m2,i_sum,i_vec,i_xyz
    type(solhrules_differential_type), pointer      :: shrd
    real(kind=r8_kind), pointer, dimension(:,:,:)   :: shg_m
    real(kind=r8_kind), pointer, dimension(:,:,:)   :: shs_m
    real(kind=r8_kind)                              :: coef
    integer(i4_kind) :: D12
    !-----------  executable part ----------------------------------

    if ( present(N_points) ) then
       gridlength = N_points
    else
       gridlength = shg%vector_length
    endif

    do l = 0,shg%lmax
       shs%l(l)%m = 0.0_r8_kind
    enddo

    ! Note how the solid harmonics for l=1 look like:
    ! sh(l=1,m=1)=z  sh(l=1,m=2)=x  sh(l=1,m=3)=y
    ! Thus, we can write the nable operator as the solid
    ! harmonic differential operator for l=1, only
    ! exchanging x,z,z

    !AM: store second derivatives in trianlge mode:
    ! D12 = 0
    ! do D1=1,3
    !    do D2=D1,3
    !      D12 = D12 + 1
    !-------------------
    ! i.e.: column-wise:
    ! xx(1)
    ! xy(2) yy(4)
    ! xz(3) yz(5) zz(6)
    ! saves 3/9 of memory

    lm2=0
    do l2 =  0,shg%lmax
       do m2 = 1, 2*l2+1
          lm2 = lm2+1
          shs_m => shs%l(l2)%m

          shrd => solhrules_differential(3,lm2) ! differential rule for x
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shg_m => shg%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 0
             do i_xyz = 1,3
                D12 = D12 + 1
                do i_vec = 1, gridlength
                   shs_m(i_vec,D12,m2) = shs_m(i_vec,D12,m2) + coef * shg_m(i_vec,i_xyz,m)
                enddo
             enddo
          enddo

          shrd => solhrules_differential(4,lm2) ! differential rule for y
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shg_m => shg%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 3 ! xz -- one before yy
             do i_xyz = 2,3 ! start from y
                D12 = D12 + 1
                do i_vec = 1, gridlength
                   shs_m(i_vec,D12,m2) = shs_m(i_vec,D12,m2) + coef * shg_m(i_vec,i_xyz,m)
                enddo
             enddo
          enddo

          shrd => solhrules_differential(2,lm2) ! differential rule for z
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shg_m => shg%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 5 ! yz -- one before zz
             do i_xyz = 3,3 ! start from z
                D12 = D12 + 1
                do i_vec = 1, gridlength
                   shs_m(i_vec,D12,m2) = shs_m(i_vec,D12,m2) + coef * shg_m(i_vec,i_xyz,m)
                enddo
             enddo
             ASSERT(D12==6)
          enddo

!!$          ! debug >>>
!!$          do i_sum=1,3
!!$             do i_xyz=1,3
!!$                print *,'sh_sd(',L2,M2,'): trans=', &
!!$                     & maxval(abs(shs_m(:,i_sum,i_xyz,M2)-shs_m(:,i_xyz,i_sum,M2)))
!!$             enddo
!!$          enddo
!!$          ! <<< debug

       enddo
    enddo


  end subroutine solid_harmonics_calc_sec_der

  subroutine solid_harmonics_calc_3rd_der( shs, sht, N_points )
    ! calculates secound derivatives of solid harmonics
    use solhrules_module
    implicit none
    !------------ Declaration of formal parameters -----------------
    type(solid_harmonics_sec_der_type), intent (in) :: shs
    type(solid_harmonics_sec_der_type), intent (inout) :: sht
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    !** End of interface ***************************************
    !-----------  private variables --------------------------------
    integer :: gridlength,lm,l,m,lm2,l2,m2,i_sum,i_vec,i_xyz
    type(solhrules_differential_type), pointer      :: shrd
    real(kind=r8_kind), pointer, dimension(:,:,:)   :: sht_m
    real(kind=r8_kind), pointer, dimension(:,:,:)   :: shs_m
    real(kind=r8_kind)                              :: coef
    integer(i4_kind) :: D12
    !-----------  executable part ----------------------------------

    if ( present(N_points) ) then
       gridlength = N_points
    else
       gridlength = shs%vector_length
    endif

    do l = 0,shs%lmax
       sht%l(l)%m = 0.0_r8_kind
    enddo

    ! Note how the solid harmonics for l=1 look like:
    ! sh(l=1,m=1)=z  sh(l=1,m=2)=x  sh(l=1,m=3)=y
    ! Thus, we can write the nable operator as the solid
    ! harmonic differential operator for l=1, only
    ! exchanging x,z,z

    !AM: store second derivatives in trianlge mode:
    ! D12 = 0
    ! do D1=1,3
    !    do D2=D1,3
    !      D12 = D12 + 1
    !-------------------
    ! i.e.: column-wise:
    ! xx(1)
    ! xy(2) yy(4)
    ! xz(3) yz(5) zz(6)
    ! saves 3/9 of memory

    lm2=0
    do l2 =  0,shs%lmax
       do m2 = 1, 2*l2+1
          lm2 = lm2+1
          sht_m => sht%l(l2)%m

          shrd => solhrules_differential(3,lm2) ! differential rule for x
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shs_m => shs%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 0
             do i_xyz = 1,6
                D12 = D12 + 1
                   !!! xxx, xyx, xzx terms
                   sht_m(:gridlength,D12,m2) = sht_m(:gridlength,D12,m2)  &
                   + coef * shs_m(:gridlength,i_xyz,m)
             enddo
          enddo

          shrd => solhrules_differential(4,lm2) ! differential rule for y
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shs_m => shs%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 6 ! xz -- one before yy
             do i_xyz = 4,6 ! start from y
                D12 = D12 + 1
                   !!! yyy yzy zzy
                   sht_m(:gridlength,D12,m2) = sht_m(:gridlength,D12,m2)  &
                 + coef * shs_m(:gridlength,i_xyz,m)
             enddo
          enddo

          shrd => solhrules_differential(2,lm2) ! differential rule for z
          do i_sum=1,shrd%n_summands
             lm = shrd%lm_sh(i_sum)
             l = solhrules_l_and_m_of_lm(1,lm)
             m = solhrules_l_and_m_of_lm(2,lm)
             shs_m => shs%l(l)%m
             coef = shrd%coef(i_sum)
             D12 = 9 ! yz -- one before zz
             do i_xyz = 6,6 ! start from z
                D12 = D12 + 1
                do i_vec = 1, gridlength
                   sht_m(i_vec,D12,m2) = sht_m(i_vec,D12,m2) + coef * shs_m(i_vec,i_xyz,m)
                enddo
             enddo
             ASSERT(D12==10)
          enddo

!!$          ! debug >>>
!!$          do i_sum=1,3
!!$             do i_xyz=1,3
!!$                print *,'sh_sd(',L2,M2,'): trans=', &
!!$                     & maxval(abs(shs_m(:,i_sum,i_xyz,M2)-shs_m(:,i_xyz,i_sum,M2)))
!!$             enddo
!!$          enddo
!!$          ! <<< debug

       enddo
    enddo


  end subroutine solid_harmonics_calc_3rd_der


  function solid_harmonics_calc(lmax,args) result(yl)
    ! this function calculates the solid harmonics for a vector of
    ! arguments for all angular momentums up to lmax
    ! In this routine l and m are treated as one metaindex
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer(kind=i4_kind) :: lmax ! maximum angular momentum
    real(kind=r8_kind),dimension(:,:),intent(in) :: args ! array of arguments
    ! args(number of args, 3 )
    real(kind=r8_kind)  :: yl(size(args,1),(lmax+1)**2) ! results
    ! first index  : over alpha, beta pairs
    ! second index : metaindex over lm
    !** End of interface ***************************************
    real(kind=r8_kind)  :: sim,rlm,rl1m
    real(kind=r8_kind),dimension(size(args,1)) :: a,b,dist
    real(kind=r8_kind) :: lr,l1r,ll1r,rll
    integer(kind=i4_kind) :: counter_old2,counter,i_l,i_m,m1,m, &
         counter_old1,ll1,ll2

    ! square of the vector argument:
    dist = args(:, 1)**2 + args(:, 2)**2 + args(:, 3)**2

    yl(:,1)=1.0_r8_kind
    if ( lmax == 0 ) return
    ! Calculate sh for l=1
    yl(:,2)=args(:,3)
    yl(:,3)=args(:,1)
    yl(:,4)=args(:,2)

    counter_old1=1
    counter=2
    do I_l = 2, lmax
       counter_old2=counter_old1 ! begin of yl(i_l-2)
       counter_old1=counter  ! begin 
       counter=counter+2*i_l-1


       ! intermediates
       lr = real( I_l, r8_kind )
       l1r = real( I_l-1, r8_kind )
       ll1r = real( 2*I_l-1, r8_kind )

       ! m = 1
       yl(:,counter) = ( ll1r * args(:,3) * yl(:,counter_old1) - &
            l1r * dist(:) * yl(:,counter_old2) ) / lr

       ! m = 2, .. , 2 * I_l - 3
       do I_m = 1, I_l - 2

          ! intermediates
          sim = real( I_m ** 2 , r8_kind )
          rlm = sqrt( lr ** 2 - sim )
          rl1m = sqrt( l1r ** 2 - sim )
          m = 2 * I_m-1
          m1 = m+1

          a(:) = ll1r * args(:,3) / rlm
          b(:) = rl1m * dist(:) / rlm
          yl(:,counter+m)  = a(:) * yl(:,counter_old1+m) - b(:) *&
               yl(:,counter_old2+m) 
          yl(:,counter+m1) = a(:) * yl(:,counter_old1+m1) - b(:) * &
               yl(:,counter_old2+m1)


       enddo

       ! m = 2 * I_l - 2,  2 * I_l - 1
       I_m = I_l - 1
       sim = real( I_m ** 2 , r8_kind )
       rlm = sqrt( lr ** 2 - sim )
       rl1m = sqrt( l1r ** 2 - sim )
       m = 2 * I_m-1
       m1 = m+1 

       a(:) = ll1r * args(:,3)/ rlm

       yl(:,counter+m)  = a(:) * yl(:,counter_old1+m)
       yl(:,counter+m1) = a(:) * yl(:,counter_old1+m1)


       ! m = 2 * I_l,  2 * I_l + 1
       rll = sqrt ( real(ll1r / ( 2 * lr ),r8_kind) )
       m = 2 * I_l-1
       m1 = m + 1
       ll1 =  2 * I_l - 2
       ll2 =  2 * I_l - 3

       yl(:,counter+m) =  ( args(:,1) * yl(:,counter_old1+ll2) - &
            args(:,2) * yl(:,counter_old1+ll1) ) * rll
       yl(:,counter+m1) = ( args(:,2) * yl(:,counter_old1+ll2) + &
            args(:,1) * yl(:,counter_old1+ll1) ) * rll

    enddo
  end function solid_harmonics_calc
  !******************************************************************

  !******************************************************************
  function solid_harmonics_scalar(lmax,arg) result(yl)
    ! this function calculates the solid harmonics for one
    ! arguments for all angular momentums up to lmax
    ! In this routine l and m are treated as one metaindex
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer(kind=i4_kind) :: lmax ! maximum angular momentum
    real(kind=r8_kind),dimension(3),intent(in) :: arg ! argument
    real(kind=r8_kind)  :: yl((lmax+1)**2) ! results
    !** End of interface ***************************************
    real(kind=r8_kind)  :: sim,rlm,rl1m
    real(kind=r8_kind) :: a,b,dist
    real(kind=r8_kind) :: lr,l1r,ll1r,rll
    integer(kind=i4_kind) :: counter_old2,counter,i_l,i_m,m1,m, &
         counter_old1,ll1,ll2


    ! Calculate sh for l=0 and r2 and intermediates x, y, z

    dist=sum(arg*arg)
    yl(1)=1.0_r8_kind
    if ( lmax == 0 ) return
    ! Calculate sh for l=1
    yl(2)=arg(3)
    yl(3)=arg(1)
    yl(4)=arg(2)

    counter_old1=1
    counter=2
    do I_l = 2, lmax
       counter_old2=counter_old1 ! begin of yl(i_l-2)
       counter_old1=counter  ! begin 
       counter=counter+2*i_l-1


       ! intermediates
       lr = real( I_l, r8_kind )
       l1r = real( I_l-1, r8_kind )
       ll1r = real( 2*I_l-1, r8_kind )

       ! m = 1
       yl(counter) = ( ll1r * arg(3) * yl(counter_old1) - &
            l1r * dist * yl(counter_old2) ) / lr

       ! m = 2, .. , 2 * I_l - 3
       do I_m = 1, I_l - 2

          ! intermediates
          sim = real( I_m ** 2 , r8_kind )
          rlm = sqrt( lr ** 2 - sim )
          rl1m = sqrt( l1r ** 2 - sim )
          m = 2 * I_m-1
          m1 = m+1

          a = ll1r * arg(3) / rlm
          b = rl1m * dist / rlm
          yl(counter+m)  = a * yl(counter_old1+m) - b *&
               yl(counter_old2+m) 
          yl(counter+m1) = a * yl(counter_old1+m1) - b * &
               yl(counter_old2+m1)


       enddo

       ! m = 2 * I_l - 2,  2 * I_l - 1
       I_m = I_l - 1
       sim = real( I_m ** 2 , r8_kind )
       rlm = sqrt( lr ** 2 - sim )
       rl1m = sqrt( l1r ** 2 - sim )
       m = 2 * I_m-1
       m1 = m+1 

       a = ll1r * arg(3) / rlm

       yl(counter+m)  = a * yl(counter_old1+m)  
       yl(counter+m1) = a * yl(counter_old1+m1) 


       ! m = 2 * I_l,  2 * I_l + 1
       rll = sqrt ( real(ll1r / ( 2 * lr ),r8_kind) )
       m = 2 * I_l-1
       m1 = m + 1
       ll1 =  2 * I_l - 2
       ll2 =  2 * I_l - 3

       yl(counter+m) =  ( arg(1) * yl(counter_old1+ll2) - &
            arg(2) * yl(counter_old1+ll1) ) * rll
       yl(counter+m1) = ( arg(2) * yl(counter_old1+ll2) + &
            arg(1) * yl(counter_old1+ll1) ) * rll

    enddo
    
  end function solid_harmonics_scalar
  !******************************************************************




end module solid_harmonics_module

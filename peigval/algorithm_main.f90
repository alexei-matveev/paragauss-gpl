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


  module algorithm_main

  !---------------------------------------------------------------
  !
  !  Purpose: Distribute a matrices between processors in optimal way
  !
  !  Module called by: peigval module
  !
  !  Author: MG
  !  Date: ?
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: SB
  ! Date:   2.3.2004
  ! Description: More convenient to use and optimal to understand
  !
  !----------------------------------------------------------------


  use type_module, only: IK=>i4_kind, RK=>r8_kind
  use linux_part

!  use mpi_headers_module

  implicit none

  private

  public :: algorithm1

  contains

!============================================================================
!============================================================================
     subroutine algorithm1( n, p, dim, distribution, count_in)

    implicit none

!   .. Parameters ..
    integer(kind=IK), intent(in) :: n, p, count_in   ! p - # of processors, n - # of matrices
    integer(kind=IK), dimension(:), intent(in) :: dim ! set of dimensionalities
    integer(kind=IK), dimension(:,:), intent(out) :: distribution ! how i must calculate

!   .. Local variables ..
    integer(kind=IK) :: summ, r, t, steps
    integer(kind=IK), dimension(n) :: unsolved_matrices

!   Executable code:

!   k_opt - corresponding optimal number of processors
!   summ - number of processors for optimal calculations

    summ = k_opt( dim(1) )
    do r = 2, n
      summ = summ + k_opt( dim(r) )
    end do

    do r = 1, n
      do t = 1, n
        distribution(r,t) = 0
      end do
    end do
!    print *,n,p,summ


    if( p >= summ ) then !case then # of processors >= # of processors that we need
       do r = 1, n
         distribution(1,r) = k_opt( dim(r) )
       end do

     else if( dim(1) < 300 .and. p < n .or. p == 2 ) then

!    else if( dim(1) >= 700 .and. p < 3 .or. dim(1) >= 500 .and. &   ! case if not
!             dim(1) < 700 .and. p < 4 ) then


       if( mod(n,p) /= 0 ) then
          steps = n / p + 1
       else
          steps = n / p
       end if

       do r = 1, p
         distribution(1,r) = 1 ! fill the columns (all processorses must solve one matrix)
       end do                                                                         !
                                                                                      !
                                                                                      !
       do t = 2, steps                                                                !
         do r = (t - 1) * p + 1, t * p                                                !
           if( r <= n ) then                                                          !
             distribution(t,r) = 1 ! fill other steps ---------------------------------
           end if
         end do
       end do
    else ! if we have matrices and a few processors

       do r = 1, n
         unsolved_matrices(r) = r
       end do

       call algorithm2( n, p, unsolved_matrices, dim, distribution, count_in )

    end if


    end subroutine algorithm1

!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================

    subroutine algorithm2( n, p, unsolved_matrices, dim, distribution, &
                           count_in )

    implicit none

!   .. Parameters ..
    integer(kind=IK), intent(in) :: n, p
    integer(kind=IK), intent(in), optional :: count_in
    integer(kind=IK), dimension(:), intent(in) :: unsolved_matrices, dim
    integer(kind=IK), dimension(:,:), intent(inout) :: distribution

!   .. Local variables ..
    integer(kind=IK) :: step, count, min_location, r, s, t, k, extent, b
    real(kind=RK) :: a
    logical :: loop, steps = .true. !Mircea's case = .false. !???????
    integer(kind=IK), dimension(n) :: unsolved_matrices_temp
    real(kind=RK), allocatable, dimension(:,:) :: time
    integer(kind=IK), allocatable, dimension(:,:,:) :: distribution_temp
    real(kind=RK), allocatable, dimension(:) :: time_sum
    integer(kind=IK), allocatable, dimension(:) :: count_forbidden
    integer(kind=IK), dimension(n) :: unsolved_matrices_out, &
                                      distribution_alg3
! My variables
    logical :: flag
    integer(IK) :: st


!   Executable code:
!    print *,dim(unsolved_matrices(1))"
    extent = k_opt( dim(unsolved_matrices(1)) )

!    print *,"extent = ",extent
!    print *,"dim(unsolved_matrices(1))=",dim(unsolved_matrices(1))

    allocate( time(n,extent), distribution_temp(n,extent,n), &
              time_sum(extent), count_forbidden(extent) )

!   set up variables Mircea's case
!    do r = 1, n
!      unsolved_matrices_temp(r) = 0   ! not really needed
!      unsolved_matrices_temp(r) = unsolved_matrices(r)
!      do s = 1, extent
!        time(r,s) = 0.
!        time_sum(s) = 99999.
!        count_forbidden(s) = 0
!        do t = 1, n
!          distribution_temp(r,s,t) = 0
!        end do
!      end do
!    end do

      unsolved_matrices_temp = 0   ! not really needed
      unsolved_matrices_temp = unsolved_matrices
      time = 0.
      time_sum = 99999.
      count_forbidden = 0
      distribution_temp = 0


    if( present(count_in) ) then
       count = count_in
      if( count <= 1 .or. count > extent ) then
        print *,"ALGORITHM: STOP /Can't find solution, choose simple case/ "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maybe it's not necessary, but
! SB implementation into Mircea's algorithm, if alg can't find solution
! we must choose a simple case all matrix solved on 1 processorces
! Why Mircea's alg can't find solution? I don't know.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if( mod(n,p) /= 0 ) then
          st = n / p + 1
       else
          st = n / p
       end if

       do r = 1, p
         distribution(1,r) = 1 ! fill the columns (all processorses must solve one matrix)
       end do                                                                         !
                                                                                      !
                                                                                      !
       do t = 2, st                                                                   !
         do r = (t - 1) * p + 1, t * p                                                !
           if( r <= n ) then                                                          !
             distribution(t,r) = 1 ! fill other steps ---------------------------------
           end if
         end do
       end do
       go to 99

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of SB implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MC:   stop

      end if
      do r = 1, count - 1
        count_forbidden(r) = r
      end do
    else
      count = 1
    end if

    b = count
    step = 1
20  loop = .false.

    if( p >= k_opt( dim(unsolved_matrices_temp(1)) ) - f( count, steps ) + 1 ) then
       k = k_opt( dim(unsolved_matrices_temp(1)) ) - f( count, steps ) + 1
    else
       k = p - f( count, steps ) + 1
    end if

    if( .not. permitted( k ) ) then
       if( steps ) then ! take the first permitted value
          do s = k - 1, 1,-1
            t = s
            if( permitted( s ) ) then
!              print *,"if( permitted( s ) ) exit"
              exit
            end if
          end do
          k = t
       else ! keep the unpermitted count for comparison
          count_forbidden(b) = count
          b = b + 1
          go to 50
       end if
    end if

    if( permitted( k ) ) then ! not really needed: a double check
!       Mircea's case
!       do r = 1, n
!         unsolved_matrices_out(r) = 0
!         distribution_alg3(r) = 0
!       end do

       unsolved_matrices_out = 0
       distribution_alg3 = 0

       call algorithm3( k, p, n, dim, unsolved_matrices_temp, &
                        unsolved_matrices_out, distribution_alg3 )
!      Mircea's case
!       do r = 1, n
!         distribution_temp(step,count,r) = distribution_alg3(r)
!       end do

       distribution_temp(step,count,:) = distribution_alg3
!       print *,"distribution_alg3"
!       print *,distribution_alg3
!       print *,"+++++++++++++++++++++++++++++++++++++++++++++++++++"
!       print *,"distribution_temp"
!       print *,distribution_temp
!       print *,"+++++++++++++++++++++++++++++++++++++++++++++++++++"


       time(step,count) = time_fit( k, dim(unsolved_matrices_temp(1)) )
       if( maxval(unsolved_matrices_out) /= 0 ) then
             step = step + 1
             do r = 1, n
               unsolved_matrices_temp(r) = 0 ! needed here
               unsolved_matrices_temp(r) = unsolved_matrices_out(r)
             end do
             loop = .true.
             steps = .true.
       end if
    end if

    if( loop ) go to 20

50  if( count + 1 <= extent ) then
       count = count + 1
       steps = .false.
       do r = 1, n
         unsolved_matrices_temp(r) = 0 ! not needed here
         unsolved_matrices_temp(r) = unsolved_matrices(r)
       end do
       go to 20
    end if

    s = 1
30  a = 0. ! build up the array "time_sum"

    do t = 1, extent
      if( s == count_forbidden(t) ) go to 60
    end do

    do r = 1, n
      a = a + time(r,s)
    end do

    time_sum(s) = a
60  s = s + 1
    if( s <= extent ) go to 30 ! end of time_sum
    t = 1

40  flag = .true.
    s = 1
     do while ( flag )
! Mircea's case:    do s = 1, extent-1 ! it doesn't have to be permitted; find its minimum loc
      if( time_sum(s) < time_sum(t) ) then
         min_location = s
! Mircea's case:        exit
          flag = .false.
      else
         min_location = t
      end if
      s = s + 1
      if ( s == extent-1 ) flag = .false.
    end do

    if( min_location /= t ) then
       t = min_location
!         print *,"GOTO 40"
         go to 40
    end if ! end of minloc

!  Mircea's case
!    do r = 1, n
!      do t = 1, n
!        distribution(r,t) = distribution_temp(r,min_location,t)
!      end do
!    end do

    distribution(:,:) = distribution_temp(:,min_location,:)

!    print *,distribution

!    print *,"START DEALLOCATE"
99    deallocate( time, distribution_temp, time_sum, count_forbidden )
!    print *,"FINISHED DEALLOCATE"


    end subroutine algorithm2
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
     subroutine algorithm3( k, p, n, dim, unsolved_matrices_in, &
                            unsolved_matrices_out, distribution )
     implicit none

!    .. Parameters ..
     integer(kind=IK), intent(in) :: n, p, k
     integer(kind=IK), dimension(:), intent(in) :: unsolved_matrices_in, dim
     integer(kind=IK), dimension(:), intent(inout) :: distribution, &
                                                    unsolved_matrices_out

!    .. Local variables ..
     integer(kind=IK) :: m, i, q, procs, s, n_temp, int_temp1, int_temp2, t

!    .. SB variables ..
     logical :: flag

!    Executable code:

     i = 1
     procs = p - k                                ! available procs for a
     distribution(unsolved_matrices_in(1)) = k    ! parallel solution
!     print *,"distribution(unsolved_matrices_in(1)) =",k
     q = 2
     if( minval(unsolved_matrices_in) /= 0 ) then ! the first positions are
        n_temp = n                                ! occupied: determine here
     else                                         ! how many
        do m = 1, n
          n_temp = m - 1
          if( unsolved_matrices_in(m) == 0 ) then
!            print *,"if( unsolved_matrices_in(m) == 0 ) exit"
            exit
          end if
        end do
     end if

10   if( time_fit( 1, dim(unsolved_matrices_in(q)) ) <= &
         time_fit( k, dim(unsolved_matrices_in(1)) ) ) then
        if( procs >= n_temp - q + 1 ) then
           do m = q, n_temp
             distribution(unsolved_matrices_in(m)) = 1
           end do
        else
           do m = i, i + n_temp - q - procs
             unsolved_matrices_out(m) = &
             unsolved_matrices_in(m + q + procs - i)
           end do
           do m = q, q + procs - 1
             distribution(unsolved_matrices_in(m)) = 1
           end do
        end if
     else
        int_temp1 = k_opt( dim(unsolved_matrices_in(q)) )
        int_temp2 = k0( dim(unsolved_matrices_in(q)) ) ! 1 to 3 instead of
        do m = int_temp1, int_temp2, -1                ! int_temp2 should
          t = m + 1                                    ! do just fine
          if( permitted( m ) .and. time_fit( m, dim(unsolved_matrices_in(q)) &
              ) > time_fit( k, dim(unsolved_matrices_in(1)) ) ) then
!            print *,"if( ... ) exit"
            exit
          end if
        end do
        do m = t, int_temp1 ! select first permitted value
          s = m
          if( permitted( m ) ) then
!            print *,"if( permitted( m ) ) exit"
            exit
          end if
        end do
        if( procs >= s ) then
           distribution(unsolved_matrices_in(q)) = s

           q = q + 1
           procs = procs - s
           if( q <= n_temp .and. procs > 0 ) go to 10
        else
           unsolved_matrices_out(i) = unsolved_matrices_in(q)
           i = i + 1
           q = q + 1
           if( q <= n_temp ) go to 10
        end if
     end if

     end subroutine algorithm3
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end module algorithm_main


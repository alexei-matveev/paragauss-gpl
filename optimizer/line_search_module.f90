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
!=====================================================================
! Public interface of module
!=====================================================================
module  line_search_module
  !-------------------------------------------------------------------
  !  Purpose: Contains all necessary routines for the line-search
  !           subproblem. Currently only a cubic fit to the
  !           Energy along the line q_last-q is implemented.
  !
  !  Module called by: step_module
  !
  !  References: Schlegel, dipl.-Arb, Gaussian...
  !  Author: FN
  !  Date: 7/97
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  use type_module ! type specification parameters
  use coortype_module
  use math_module
  use opt_data_module !!!!!!!!!!
  use gradient_module, only: grad_intern,energy,energy_ph
  use iounitadmin_module
  use filename_module, only: inpfile
  
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of constants and variables ---------------
  integer(kind=i4_kind):: io_linmin = 20
  real(kind=r8_kind),allocatable        :: q_last(:),g_last(:),line_step(:)
  real(kind=r8_kind)                    :: e_last
  logical                               :: last_line_search
  real(kind=r8_kind)                    :: g_last_proj,g_proj
  real(kind=r8_kind)                    :: a3,a2,a1,a0
  ! io_linmin : io-unit for the file 'line_search.dat'
  ! q_last,g_last : coordinates and gradient of the previous point
  ! e_last        : energy of the last point
  ! g_last_proj   : g_last*(q-q_last) = component of the last gradient
  !                 along the line q-q_last
  ! g_proj        : grad_intern*(q-q_last)
  ! a3 - a0       : coefficients for the cubic polynomial fit
  !                 this can be easily extended to a quartic fit
  !                 once I know how this is done.
  
  public line_search_main, line_search_main2
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  !*************************************************************
  subroutine line_search_read()
    !  Purpose: read in the file line_search.dat. This file
    !           contains the following information:
    !           -last coordinate vector q_last
    !           -last gradient vector g_last (internal cordinates)
    !           -last Energy e_last
    !           -last_line_search : TRUE if last point was
    !            the result of a line search
    !                               FALSE if not.
    !------------ Modules used ---------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)      :: alloc_stat
    !------------ Executable code -------------------------------
    allocate(q_last(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_read : allocation (1) failed'
    allocate(g_last(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_read : allocation (2) failed'
    allocate(line_step(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_read : allocation (3) failed'
    line_step=zero
    
    io_linmin=openget_iounit(status='unknown',form='formatted',&
         file=trim(inpfile('line_search.dat')))
!         file=trim(opt_data_dir)//'/line_search.dat')

    read(io_linmin,*)q_last
    read(io_linmin,*)g_last
    read(io_linmin,*)e_last,last_line_search
    call returnclose_iounit(io_linmin)
  end subroutine line_search_read
   !*************************************************************
  subroutine line_search_write(dealloc)
    !  Purpose: write in the file line_search.dat. This file
    !           contains the following information:
    !           -last coordinate vector q_last
    !           -last gradient vector g_last (internal cordinates)
    !           -last Energy e_last
    !           -last_line_search : TRUE if last point was
    !            the result of a line search
    !                               FALSE if not.
    !------------ Modules used ---------------------------------
    !------------ Declaration of formal parameters ---------------
    logical,optional    :: dealloc
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)      :: alloc_stat
    !------------ Executable code -------------------------------
    
    io_linmin=openget_iounit(status='unknown',form='formatted',&
         file=trim(inpfile('line_search.dat')))
!         file=trim(opt_data_dir)//'/line_search.dat')

    write(io_linmin,*)q
    write(io_linmin,*)grad_intern
    write(io_linmin,*)energy,last_line_search
    call returnclose_iounit(io_linmin)

    if (present(dealloc).and.dealloc.or. .not.present(dealloc) ) then
       deallocate(q_last,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop 'line_search_write : deallocation (1) failed'
       deallocate(g_last,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop 'line_search_write : deallocation (2) failed'
       deallocate(line_step,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop 'line_search_write : deallocation (3) failed'
    endif
  end subroutine line_search_write
   !*************************************************************
  subroutine line_search_main(do_newton,geo_loop)
    ! Purpose: main routine for the line search.
    !          1. Test if line search makes sense at all:
    !          g_last*(q-q_last) = g_last_proj < 0
    !                   AND
    !          grad_intern*(q-q_last) = g_proj > 0
    !          i.e. along the line connecting q and q_last lies
    !          at least one local minimum.
    !          If test is positive, do a line search by cubic fit
    !          and finding the minima of the first derivative of this
    !          polynomial.
    !
    !          The newstep is found by:
    !          q_new = 1/2*(q+q_last) + (alpha-1/2)*(q-q_last)
    !
    !          2. If the above test is negative, make another newton step
    !          with fixed limitation. Lateron methods will be implemented
    !          that checks the history of steps...
    !   
    !          3. In case of successfull line search, two options are
    !          possible:
    !          - produce a step only. This 'step' is sym-checked and
    !          the updated coordinate 'q' is given to the routine
    !          'internal_to_cart'.
    !          - produce a step (sym-check it!) and estimate the gradient
    !          at the new position 'q'. Take this as new input variables
    !          'q' and grad_intern and do a newton step thus saving one
    !          point of DFT-calculation.
    !------------ Modules used ---------------------------------
    use coordinates_module,only: sym_check,q_old
    use hesse_module
    !------------ Declaration of formal parameters ---------------
    logical,intent(out)                :: do_newton
    integer(kind=i4_kind),intent(in)   :: geo_loop
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable     :: delta_q(:),help_grad(:)
    integer(kind=i4_kind)              :: alloc_stat,i,info_fit
    real(kind=r8_kind),parameter       :: small_alpha = 1.0e-4_r8_kind
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)   :: alpha
    
    allocate(delta_q(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("line_search_main: allocation (1) failed")
    delta_q = zero
    do_newton=.true.
    a3=zero
    a2=zero
    a1=zero
    a0=zero
    alpha=zero
    if (geo_loop==1) then
       call line_search_write(dealloc=.false.)
       return
    end if
    call line_search_read()
    do i=1,n_internal
       delta_q(i) = q(i)-q_last(i)
       if (s(i)%typ == d_angle ) then
          if (delta_q(i) > pi ) delta_q(i) = -two*pi+delta_q(i)
          if (delta_q(i) < - pi ) delta_q(i) =  two*pi+delta_q(i)
       endif
    enddo
    
    !test
    g_last_proj = dot_product(g_last,delta_q) / abs_value(delta_q) ! AG ! norm added
    g_proj = dot_product(grad_intern,delta_q) / abs_value(delta_q) ! AG ! norm added
    !AG : line-coordinates redefined
    write(OPT_STDOUT,"(' line_search_main : q1 q2 =',2(f14.6,4x))") 0.0_r8_kind, abs_value(q - q_last)
    !-------------------------------
    
    if ( g_last_proj*g_proj<zero ) then
       write(io_flepo,*)" line_search_main:  Line search will be attempted  "
       write(io_flepo,'(" Energies: e_last, e_actual    :",2(2x,f20.6))')e_last,energy
       write(io_flepo,'(" g_last_proj                   :",1ES12.5)')g_last_proj
       write(io_flepo,'(" g_proj                        :",1ES12.5)')g_proj
       call cubic_fit(info_fit)
       write(OPT_STDOUT,*) 'info_fit',info_fit
       if (info_fit/=0) then
          do_newton=.true.
       else
          call find_zero(alpha)
          if (alpha/=zero) then
             if (abs(alpha) <= small_alpha) then
                write(OPT_STDOUT,*)" line_search_main :  Line search successfull, but alpha below threshold"
                write(OPT_STDOUT,*)"                     Do an Newton step"
             else
                write(OPT_STDOUT,*)" line_search_main : Line search succesfull - no Newton step required"
                do_newton = .false.
                ! put the newly found q in the step format such that its symmetry can be
                ! be checked by the same routine as the newton step:
                line_step = zero
                where (s%var)
                   line_step = (alpha-half)*delta_q
                end where
                write(OPT_STDOUT,*) 's%var'
                write(OPT_STDOUT,*) s%var
                write(OPT_STDOUT,*) line_step 
        
                write(OPT_STDOUT,*) q_last
                write(OPT_STDOUT,*) q 
                q = half*(q_last+q) + line_step
                write(OPT_STDOUT,*) q
                if (estimate_grad) then
                   !This routine will estimate the gradient at the new point 'q'
                   !and update the variable 'grad_intern'.
                   allocate(help_grad(n_internal),STAT=alloc_stat)
                   if (alloc_stat/=0) call error_handler &
                        (" line_search_main: allocation (2) failed")
                   help_grad=grad_intern
                   call line_search_grad(delta_q,alpha)
                   do i=1,n_internal
                      delta_q(i) = q(i) - q_old(i)
                      if (s(i)%typ == d_angle ) then
                         if (delta_q(i) > pi ) delta_q(i) = -two*pi+delta_q(i)
                         if (delta_q(i) < pi ) delta_q(i) =  two*pi+delta_q(i)
                      endif
                   enddo
                   deallocate(help_grad,STAT=alloc_stat)
                   if(alloc_stat/=0) call error_handler &
                        (" line_search_main: deallocation (2) failed")
                   do_newton=.false.
                endif
             endif
          endif
       endif
    else
       write(io_flepo,*)" No Line search possible "
       write(io_flepo,'(" Energies: e_last, e_actual    :",2(2x,F20.6))')e_last,energy
       write(io_flepo,'(" g_last_proj                   :",1ES12.5)')g_last_proj
       write(io_flepo,'(" g_proj                        :",1ES12.5)')g_proj
       do_newton=.true.
    endif

    deallocate(delta_q,STAT=alloc_stat)
    if (alloc_stat/=0 ) call error_handler &
         ("line_search_main: deallocation (1) failed")

    call line_search_write(dealloc=.true.)
  end subroutine line_search_main
   !*************************************************************
  subroutine cubic_fit(info_fit)
    ! Purpose: Does a one-dimensional cubic polynomial fit to the
    !          the energy along the line connecting
    !          q_last and q, using the information
    !          - E1 = e_last(q_last)
    !          - E2 = energy(q)   => Why not energy_ph???
    !          - g1 = g_last_proj (see line_search_main)
    !          - g2 = g_proj (see line_search_main) :
    !
    !          e_fit = a3*q**3 + a2*q**2 + a1*q +a0
    !    
    !          Formal parameters:
    !          info_fit = 0   cubic fit successfull
    !          info_fit /=0   cubic fit failed.
    !           -1 : interpolation points (absolute values of coordinate
    !                vectors in internal ccords. are the same within
    !                assumed accuracy.
    !           -2 : LU factorization (dgetrf) failed
    !           -3 : solution of LES (dgetrs) failed
    !           -4 : test1 and test2 /= e_last,energy,see below.
    !
    !------------ Modules used ---------------------------------
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)   :: info_fit
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------  
    real(kind=r8_kind)      :: a_mat(4,4),b_vec(4),q1,q2,test1,test2,qi
    integer(kind=i4_kind)   :: ipiv(4),info,i,ii
    real(kind=r8_kind),parameter  :: small_delta=1.0e-8_r8_kind
    real(kind=r8_kind),parameter  :: small_edelta=5.0e-4_r8_kind

    external dgetrs,dgetrf

    q1 = abs_value(q_last)
    q2 = abs_value(q)

!AG : line-coordinates redefined
    q1 = 0.0_r8_kind
    q2 = abs_value(q - q_last)
    write(OPT_STDOUT,"(' line_search_main/cubic_fit : q1 q2 =',2(f14.6,4x))") q1,q2
!-------------------------------

        info_fit = 0

    if (abs(q1-q2)<=small_delta) then
       write(OPT_STDOUT,"(1x,'q1  q2 ',2(f14.6,4x))") q1, q2
       write(OPT_STDOUT,*) ' line_search_main/cubic_fit : interpolation points are too close'
       write(OPT_STDOUT,*) '                             cubic fit failed - do a newton step'
       info_fit = -1
       return
    endif

    ! fill in b_vec
    b_vec(1) = e_last
    b_vec(2) = energy
    b_vec(3) = g_last_proj
    b_vec(4) = g_proj

    ! fill in a_mat
    a_mat(1,1) = q1**3
    a_mat(1,2) = q1**2
    a_mat(1,3) = q1
    a_mat(1,4) = one

    a_mat(2,1) = q2**3
    a_mat(2,2) = q2**2
    a_mat(2,3) = q2
    a_mat(2,4) = one

    a_mat(3,1) = three*q1**2
    a_mat(3,2) = two*q1
    a_mat(3,3) = one
    a_mat(3,4) = zero

    a_mat(4,1) = three*q2**2
    a_mat(4,2) = two*q2
    a_mat(4,3) = one
    a_mat(4,4) = zero

    ! To do the cubic fit, solve the linear system of equations:
    ! a_mat * x = b_vec. This is done by 'dgetrs' which uses the
    ! LU - factorization of 'dgetrf'
    call dgetrf(4,4,a_mat,4,ipiv,info)
    if (info/=0) then
       write(OPT_STDOUT,*)' line_search/cubic_fit : LU factorization (dgetrf) failed'
       write(OPT_STDOUT,*)'                         cubic fit failed - do a newton step'
       info_fit = -2
       return
    endif

    call dgetrs('N',4,1,a_mat,4,ipiv,b_vec,4,info)
    if (info/=0) then
       write(OPT_STDOUT,*)' line_search/cubic_fit : solution of LES (dgetrs) failed'
       write(OPT_STDOUT,*)'                         cubic fit failed - do a newton step'
       info_fit = -3
       return
    endif

    ! test the fit
    !AG: E-profile (start printing)
    write(io_flepo,*) "********* Cubic fit E-profile *********"
    do ii=1,11
       qi = q1 + (ii-1)*(q2-q1)/10
       test1 =  b_vec(1)*qi**3 + b_vec(2)*qi**2 + b_vec(3)*qi + b_vec(4)
       write(io_flepo,"(3x, 'q E(q) ',2(f18.10,3x))") qi,test1 
    end do
    !AG: E-profile (end printing)
    test1 = b_vec(1)*q1**3 + b_vec(2)*q1**2 + b_vec(3)*q1 + b_vec(4)
    test2 = b_vec(1)*q2**3 + b_vec(2)*q2**2 + b_vec(3)*q2 + b_vec(4)

    if (print_debug) then
       write(OPT_STDOUT,1000)(b_vec(i),i=1,4)
       write(OPT_STDOUT,1020)test1,e_last
       write(OPT_STDOUT,1030)test2,energy
    endif
    if (abs(test1-e_last) >= small_edelta .or. abs(test2-energy) >= small_edelta ) then
       write(OPT_STDOUT,*)' line_search/cubic_fit : fit does not reproduce energy values'
       write(OPT_STDOUT,*)'                         cubic fit failed - do a newton step'
       write(OPT_STDOUT,*) abs(test1-e_last),small_edelta
       write(OPT_STDOUT,*) abs(test2-energy),small_edelta
       info_fit=-4
       return
    endif
    
    ! fill in the coefficients:
    if(q2 < q1) b_vec = -b_vec ! AG
    a3=b_vec(1)
    a2=b_vec(2)
    a1=b_vec(3)
    a0=b_vec(4)

1000 format("coeffs of cubic polynomial  :",4(1ES12.5,3x))
1020 format("cubic_fit : test1 , e_last  :",2(f14.6,3x))
1030 format("cubic_fit : test2 , energy  :",2(f14.6,3x))

  end subroutine cubic_fit
   !*************************************************************
  subroutine find_zero(alpha)
    ! Purpose: find the zeros of the first derivative of the
    !          fitted polynomial. For a cubic fit this is a
    !          simple p-q-formula. For a quartic fit either
    !          a regula-falsi or Newton-procedure has to be done.
    ! 
    ! Subroutine called by: line_search_main
    !
    !------------ Modules used ---------------------------------
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(out)   :: alpha
    ! alpha is the zero of the first derivative of the fitted
    ! polynimial, i.e. the size of the new step along the line
    ! q-q_last.
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)    :: x1,x2,d,q1,q2,test
    integer(kind=i4_kind) :: counter
    logical               :: is1,is2
    !--- executable code -----------------------------------------
    d = (a2/(three*a3))**2 - a1/(three*a3)
    if (d < zero) then
       write(OPT_STDOUT,*)" line_search/find_zero : no real solution to the linmin problem could be found"
       alpha = zero
       return
    endif

    is1=.false.
    is2=.false.
    x1 = -(a2/(three*a3)) + sqrt(d)
    x2 = -(a2/(three*a3)) - sqrt(d)
    q1 = abs_value(q_last)
    q2 = abs_value(q)
!AG : line-coordinates redefined
    q1 = 0.0_r8_kind
    q2 = abs_value(q - q_last)
!-------------------------------
    write(OPT_STDOUT,*)" line_search/find_zero    :"
    write(OPT_STDOUT,1000)q1,q2,x1,x2
1000 format("  interpolation pts: q1,q2 ",2(f11.6,2x),3x," zeros of cubic polynomial  x1,x2 ",2(f11.6,2x))
    if (x1==q1.or.x1==q2) then
       stop ' line_search/find_zero : the gradient must have been zero at beginning of fit'
    elseif (x2==q1.or.x2==q2) then
       stop ' line_search/find_zero : the gradient must have been zero at beginning of fit'
    endif

    counter=0
    is1=.false.
    is2=.false.
    if ( (x1>q1.and.x1<q2) .or. (x1>q2.and.x1<q1) ) then
       test = 6.0_r8_kind*a3*x1 + two*a2
       if (test < zero ) then
          write(OPT_STDOUT,*)" line_search/find_zero: x1 is a maximum."
       else
          counter=counter+1
          is1=.true.
       endif
    else
       write(OPT_STDOUT,*)" line_search/find_zero: x1 is beyond the interpolation points "
    endif
    if( (x2>q1.and.x2<q2) .or. (x2>q2.and.x2<q1)) then
       test = 6.0_r8_kind*a3*x2 + two*a2  
       if (test < zero ) then
          write(OPT_STDOUT,*)" line_search/find_zero: x2 is a maximum."
       else
          counter=counter+1
          is2=.true.
       endif
    else
       write(OPT_STDOUT,*)" line_search/find_zero: x2 is beyond the interpolation points "
    endif
       

    if (counter == 0 ) then
       write(OPT_STDOUT,*)" line_search/find_zero: line_search failed. No minimum found"
       alpha=zero
    elseif(counter == 1 ) then
       if (is1) then
          alpha = (q1-x1)/(q1-q2)
       else
          alpha = (q1-x2)/(q1-q2)
       endif
    elseif(counter == 2 ) then
       write(OPT_STDOUT,*)" line_search/find_zero: both points are minima. Proceed with x1"
       alpha = (q1-x1)/(q1-q2)
    else
       write(OPT_STDOUT,*)" rubbish in routine FIND_ZERO"
       alpha=zero
    endif
    write(OPT_STDOUT,*)" line_search/find_zero : alpha is     = ",alpha
    !AG: scale "alpha to make q-update consistent with line-coordinate"
    !alpha = alpha*(q2-q1)/abs_value(q-q_last)
    !write(OPT_STDOUT,*) " line_search/find_zero : scaled alpha =", alpha
  end subroutine find_zero
  
   !*************************************************************
  subroutine line_search_grad(delta_q,alpha)
    ! Purpose: estimate the gradient at the new position 'q' after
    !          succesfull line_search.
    ! Input: delta_q = q - q_last, referring to the not yet updated
    !        variable q.
    !        alpha : line_search parameter resulting from 'find_zero'.
    !        Note that in the notoaion of F. Eckerts Dipl.Arb. this
    !        is used as xm = alpha - 1/2.
    ! ----------------------------------------------------------
    real(kind=r8_kind),intent(in)   :: delta_q(:)
    real(kind=r8_kind),intent(in)   :: alpha
    ! ----- Declaration of local variables ---------------------
    real(kind=r8_kind)   :: xt
    real(kind=r8_kind),allocatable :: help(:)
    integer(kind=i4_kind) :: alloc_stat,i

    write(OPT_STDOUT,*)" "
    allocate(help(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         (" line_search_grad : allocation (1) failed")

    help = half*(grad_intern+g_last) + &
         (alpha-half)*(grad_intern-g_last)

    xt = dot_product(delta_q,help)/dot_product(delta_q,delta_q)
    
    grad_intern = zero
    where (s%var)
       grad_intern = help - xt*delta_q
    end where
    write(OPT_STDOUT,*)"  line_search_grad : estimated gradient ="
    do i=1,n_internal
       write(OPT_STDOUT,1000)i,grad_intern(i)
    enddo
    write(OPT_STDOUT,*)" "
1000 format(3x,i3,f13.10)
    
    deallocate(help,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         (" line_search_grad: deallocation (1) failed")

  end subroutine line_search_grad

   !*************************************************************
  subroutine line_search_main2(do_newton,geo_loop)
    !------------ Declaration of formal parameters ---------------
    logical,intent(out)                :: do_newton
    integer(kind=i4_kind),intent(in)   :: geo_loop
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable     :: delta_q(:)
    integer(kind=i4_kind)              :: alloc_stat,i
    real(kind=r8_kind),parameter       :: small_alpha = 1.0e-4_r8_kind
    real(kind=r8_kind)   :: alpha
    !** End of interface *****************************************
    
    if (mod(geo_loop,2) /= 0) then
       do_newton=.true.
       call line_search_write(dealloc=.false.)
       return
    end if
    do_newton=.false.
       
    allocate(delta_q(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("line_search_main: allocation (1) failed")
    delta_q = zero

    call line_search_read()
    do i=1,n_internal
       delta_q(i) = q(i)-q_last(i)
       if (s(i)%typ == d_angle ) then
          if (delta_q(i) > pi ) delta_q(i) = -two*pi+delta_q(i)
          if (delta_q(i) < -pi ) delta_q(i) =  two*pi+delta_q(i)
       endif
    enddo
    
    !test
    g_last_proj = dot_product(g_last,delta_q)

    call square_fit()

    line_step = zero
    where (s%var)
       line_step=(alpha-half)*delta_q
    end where
    q=half*(q_last+q) + line_step

    deallocate(delta_q,STAT=alloc_stat)
    if (alloc_stat/=0 ) call error_handler &
         ("line_search_main2: deallocation (1) failed")
    deallocate(q_last,STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_write2 : deallocation (2) failed'
    deallocate(g_last,STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_write2 : deallocation (3) failed'
    deallocate(line_step,STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop 'line_search_write2 : deallocation (4) failed'



  contains
    subroutine square_fit()
      !------------ Modules used ---------------------------------
      !------------ Declaration of formal parameters -------------
      !** End of interface ***************************************
      !------------ Declaration of local variables ---------------
      real(kind=r8_kind) :: q1,q2,aa0,aa1,aa2,e_m,q_m,beta
      real(kind=r8_kind),parameter  :: small_edelta=1.0e-10_r8_kind

      q1 = abs_value(q_last)
      q2 = abs_value(q)

      beta=1.0_r8_kind
      aa0=e_last
      aa1=g_last_proj
      aa2=(energy-e_last-g_last_proj*beta)/beta**2

      if(abs(energy-e_last) <= small_edelta) then
         q_m=beta
      elseif(aa2 > 0.0_r8_kind) then
         q_m=-aa1/(2.0_r8_kind*aa2)
         if((energy-e_last) < 0.0_r8_kind) then
            if(q_m >= 5.0_r8_kind*beta) q_m=5.0_r8_kind*beta
         else
            if(q_m > beta) q_m=beta/2.0_r8_kind
         end if
      else
         if((energy-e_last) < 0.0_r8_kind) then
            q_m=2.0_r8_kind*beta
         else
            q_m=beta
         end if
      end if

      e_m=aa2*q_m**2+aa1*q_m+aa0
      alpha=q_m/beta

      write(OPT_STDOUT,*)" line_search2: square fit"
      write(OPT_STDOUT,'(3(a5,f19.11))')" x0 :",0.0    ," x1 :",beta    ," xm :",q_m
      write(OPT_STDOUT,'(3(a5,f19.11),/)')" E0 :",e_last," E1 : ",energy," Em :",e_m
      
    end subroutine square_fit
  end subroutine line_search_main2


!--------------- End of module ----------------------------------
end module line_search_module

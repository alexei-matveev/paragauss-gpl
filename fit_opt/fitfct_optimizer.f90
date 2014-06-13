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
!
!____________ Optimizer program______________________________________________
!
! Author: Valentin Igoshin
! begin : 19.05.99
!
! ---------------------------------------------------------------------------
! Multidimensional MINIMIZATION OF THE FUNCTION funk(x) where x(1:ndim)
! is a vector in ndim dimensions, by the downhill simplex method of Nelder
! and Mear. The matrix p(1:ndim+1,1:ndim) is input. Its ndim+1 rows
! are ndim-dimensional vectors which are the vertices of the starting simplex
! Also input is the vector y(1:ndim+1), whose components must be the
! pre-initialized to the values of funk evaluated at the ndim+1
! vertices (rows) of p; and ftol the fractional convergence tolerance to be 
! achived in the function value. On output, p and y will have been reset to 
! ndim+1 new points all within ftol of a minimum function value, and 
! iter gives the number of function evaluations taken.
!
! "Numerical recipies in Fortran" William H. Press and others
!  
!  Program was tested on the multidimenzion surfaces. Now as a function we
!  use ParaGauss program and Coulomb integral is the result of the function. 
!----------------------------------------------------------------------------
!
PROGRAM fitfct_optimizer 
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER                 :: r8_kind = SELECTED_REAL_KIND(15)
  !
  ! NMAX - maximal number of dimensions;  ITMAX - maximal number of iterations
  INTEGER, PARAMETER                 :: NMAX=40, ITMAX=300 
  !
  ! In this type we have the input(start) value in the opt_inp file
  ! One modification: during the program all exponents substitut by its logarithms
  TYPE var
     CHARACTER(len=1)                :: type
     INTEGER                         :: lm, n_exp
     REAL(r8_kind)                   :: alpha_min, alpha_max, q
  END TYPE var
  !
  ! This type we use when printing the fit exponents in the ch_fit.op file
  ! In this type we use real exponents ( not the logarithms )
  TYPE fit
     CHARACTER(len=1)                :: type
     INTEGER                         :: lm, n_exp
     REAL(r8_kind)                   :: q
     REAL(r8_kind), POINTER          :: expn(:)
  END TYPE fit
  !
  LOGICAL                            :: new_iter
  !
  INTEGER                            :: i, j, k, mp, ndim, np, iter, miny, l, loc(1)
  INTEGER                            :: output_level, l_max, error, status, n 
  !
  REAL(r8_kind)                      :: delta, ftol, rtol, b, Ecoul, s
  REAL(r8_kind)                      :: e_target, delta_e
  REAL(r8_kind), ALLOCATABLE         :: p(:,:), y(:), x(:), result(:)
  !
  TYPE(var), ALLOCATABLE, TARGET     :: var_set(:)
  TYPE(var), POINTER                 :: v
  TYPE(fit), ALLOCATABLE, TARGET     :: fit_set(:)
  TYPE(fit), POINTER                 :: f
  !
  namelist/control/output_level, ftol, l_max, e_target, delta_e
  !
  ! output_level = 0 - only start value of the function and result value of this function after
  !                    miminization in the file opt_output
  !
  ! output_level = 1 - also value of function on each step in the file opt_output
  !
  ! output_level = 2 - also input arrays p and coordinats on each step in the file opt_output
  !
  ! Open file opt_output. In this file we'll have output of the optimizer
  !
  open(unit=6, file="opt_output",IOSTAT=status)
  if( status /= 0 ) print *,"open file opt_output error in main_opt  subroutine"
  write (unit=6,fmt='(A)')"___________ Output of the optimizer____________ "
  !
  ftol=0.0001 ! convergence criterion, default
  new_iter = .false.
  n = 0       
  !
  ! Read input of the optimizer. This subroutine give us
  ! the start point for optimization, ndim, ftol and output_level
3 print *," "
  call read_input(ftol)
  print *," Convergency creiterion ftol=", ftol
  print *," ndim= ",ndim
  print *," e_target =", e_target
  print *," delta_e  =", delta_e
  print *," ndim     =", ndim
  !
  !ndim=(l_max+2)*3 ! number of dimensions in our case
  !
  mp = ndim+1
  np = ndim
  ALLOCATE (p (mp,np))
  ALLOCATE (y (mp))
  delta = 0.1 ! use when create p array
  iter=0   ! number of iterations
  ! create an input for amoeba subroutine
  ! make p array - include coordinates of all ndim+1 points
  !
  p(1,:) = x(:)
  if (output_level == 2) then
     write (unit=6,fmt=*)"Start point x :"
     do i=1, ndim
        write (unit=6,fmt=*)"x(",i,")=",x(i)
     enddo
  endif
  ! setting up the ndim+1 starting points
  print *," Corrections of the elements in P array"
  if (output_level == 2) write (unit=6,fmt=*) " Corrections of the elements in P array"
  do i=2, ndim+1
     do j=1,ndim
        p(i,j)=x(j)
     enddo
     p(i,i-1) = p(i,i-1) + delta*p(i,i-1)
     print *,"i=",i," p(",i,",",i-1,") ",p(i,i-1)
     if (output_level == 2) then
        write (unit=6,fmt=*)"i=",i," p(",i,",",i-1,") ",p(i,i-1) 
     endif
  enddo
  if (output_level == 2) then
     write (unit=6,fmt=*) "BEGIN Start array"
     do i=1,ndim+1
        write (unit=6,fmt=*) "Start array"
        do j=1,ndim
           write (unit=6,fmt=*) "p(",i,",",j,")=",p(i,j) 
        enddo
     enddo
     write (unit=6,fmt=*) "END Start array"
  endif
  ! Calculate the start values of function in ndim+1 point
  write (unit=6,fmt=*) " Convergency creiterion ftol=", ftol
  print *," Convergency creiterion ftol=", ftol
  k=0
  y=0
  do while (k<=ndim)
     k=k+1
     print *,"Start point number ********>",k
     y(k)=funk(p(k,:))
     write (unit=6,fmt=*) "start point number ",k,"  y=",y(k)  
  enddo
  ! End calculate of the start point
  !
  ! Start optimization algorithm
  call amoeba(p,y,mp,np,ndim,ftol,iter)
  ! End miminization algorithm
  !
  if (iter >= ITMAX) then
     write (unit=6,fmt=*) "OPTIMIZER: convergence is NOT reached  " 
     write (unit=6,fmt=*) "OPTIMIZER: convergence criterion: rtol =",rtol
     print *,"OPTIMIZER: convergence is NOT reached  "
     print *,"OPTIMIZER: convergence criterion: rtol =",rtol
  else
     write (unit=6,fmt=*) "OPTIMIZER: convergence reached in ",iter,"iterations"
     write (unit=6,fmt=*) "OPTIMIZER: convergence criterion: rtol =",rtol
     print *,"OPTIMIZER: convergence reached in ",iter,"iterations"
     print *,"OPTIMIZER: convergence criterion: rtol =",rtol
  endif
  k=0
  b=1.0E+20  ! border of values of function funk
  !
  ! find the best point in p array
  ! P.S. In principle all point in p
  ! have a good accuracy 
  loc = MINLOC(y)
  miny=loc(1)
  print *,"OPTIMIZER: optimal value of function is y(",miny,")=",y(miny)
  print *,"OPTIMIZER: optimal value of variables"
  do i=1, ndim
     print *,"OPTIMIZER: x(",i,")=",p(miny,i)
     if(abs(p(miny,i))> b ) k=k+1
  enddo
  if (k>0) then
     print *,"WARNING: Absolute value of ",k," variables is mach than",b
     print *,"WARNING: Your function is infinite, minimum is NOT found"
  endif
  ! Call function "funk" for optimized value
  ! then we will have optimized charge fit in file "ch_fit.op"
  s = funk(p(miny,:))
  !
  if ( delta_e /= 0) then
     if ( (e_target+s) > delta_e ) then
        if (new_iter) then
           print *," n_exp =", n
           stop
        else
           n = n+1
        endif
     else if ( (e_target+s) < delta_e) then
        ! the fit is good 
        ! then we keep our date in file "ch_fit.store".
        call store_fit
        if ( (e_target+s) > (delta_e - delta_e*0.1) ) then
           print *,"number of exponents in the optimized case =", n
           stop
        endif
        n = n-1
        new_iter = .true.
     endif
  endif
  !
  write (unit=6,fmt=*) " ------------------------------------------------- "
  write (unit=6,fmt=*) " "
  write (unit=6,fmt=*) " FUNCTION AFTER MIMINIZATION =", y(miny)
  write (unit=6,fmt=*) " "
  call from_x_to_var_result_print(p(miny,:))
  !
  ! close opt_output file
  close(unit=6)
  ! Deallocating 
  DEALLOCATE (p,y)
  DEALLOCATE (var_set)
  DEALLOCATE (fit_set)
  DEALLOCATE (x)
  if ( delta_e /= 0) goto 3
  !  
  !
CONTAINS
    !
    SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,iter)
    !
    ! Optimization algorithm (simplex method)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT)           :: iter
    INTEGER, INTENT(IN)              :: mp, ndim, np
    REAL(r8_kind), INTENT(INOUT)     :: p(1:mp,1:np), y(1:mp)
    REAL(r8_kind), INTENT(IN)        :: ftol
    ! *** end of interface ***
    INTEGER                          :: i, ihi, ilo, inhi, j, m, n
    REAL(r8_kind)                    :: summ, swap, ysave, ytry, psum(NMAX) 
    !
    print *,"ameoba enter "
    do i=1, mp
       print *,"x=",p(i,:),"y=",y(i)
    enddo
    iter=0
    ! Enter here when starting or have just overal contracted.
1   do n=1, ndim   
       summ = 0  !Recompute psum.
       do m=1, ndim+1
          summ = summ+p(m,n)
       enddo
       psum(n) = summ
    enddo
2   ilo=1   !Enter here when have just changed a single point
    ! Determine which point is the highest (worst), next-highest,
    ! and lowest (best)
    do i=1, mp
       print *,"iter",iter,"x=",p(i,:),"y=",y(i)
       if (output_level > 0) then 
          write (unit=6,fmt=*)"iter",iter,"  y=",y(i)
       endif
       if (output_level > 0) write (unit=6,fmt=*)" "
       if (output_level == 2) call  from_x_to_var(p(i,:))
    enddo
    print *,"------------ new iteration -----------------------  "
    if (output_level > 0) then 
       write (unit=6,fmt=*) "------------ new iteration -----------------------  "
       write (unit=6,fmt=*) " "
    endif
    if (y(1) > y(2)) then
       ihi = 1
       inhi = 2
    else
       ihi = 2
        inhi = 1
    endif
    do i=1, ndim+1  ! by looping over the points in the simplex
       if (y(i) <=  y(ilo)) ilo = i
       if (y(i) >   y(ihi)) then
          inhi = ihi
          ihi = i
       else if(y(i) > y(inhi)) then
          if (i /= ihi) inhi = i
       endif
    enddo
    rtol = 2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
    print *,"------- rtol = ",rtol
    ! compute the fractional range from highest to lowest and return if 
    ! satisfactory.
    if (rtol < ftol) then ! If returning, put best point and value in slot 1.
       swap = y(1)
       y(i) = y(ilo)
       y(ilo) = swap
       do n=1, ndim
          swap = p(1,n)
          p(1,n) = p(ilo,n)
          p(ilo,n) = swap
       enddo
       !print *,"amoeba: convergence reached in ",iter,"iterations"
       !print *,"convergence criterion: rtol =",rtol
       return
    endif
    !if (iter >= ITMAX) pause 'ITMAX exceeded in amoeba'
    if (iter >= ITMAX) then
       print *," subroutine amoeba: Aborting program at the maximal number of iterations"
       print *," subroutine amoeba: Maximal number of iterations is ITMAX =", ITMAX       
       write (unit=6,fmt=*) " subroutine amoeba: Aborting program at the maximal number of iterations"
       write (unit=6,fmt=*) " subroutine amoeba: Maximal number of iterations is ITMAX =", ITMAX       
       return
    endif
    iter = iter + 2
    ! Begin a new iterration. First extrapolate by a factor -1
    ! through the face of the simplex across from the high point,
    ! i.e. reflect the simplex from the high point.
    ytry = amotry(p, y, psum, mp, np, ndim, ihi, -1.0)
    if (ytry <= y(ilo)) then
       ! Gives a result better than the best point, so try an additional
       ! extrapolation by a factor 2.
       ytry=amotry(p, y, psum, mp, np, ndim, ihi, 2.0)
    else if(ytry >= y(inhi)) then
       ! the reflected point is worse than the second-highest, so
       ! look for an intermediate lower point, i.e., do a one-dimensional
       ! contraction.
       ysave = y(ihi)
       ytry=amotry(p, y, psum, mp, np, ndim, ihi, 0.5)
       ! Can't seem to get rid of that high point. Better contact around the
       ! lowest (best) point.
       if (ytry >= ysave) then 
          do i=1, ndim
             if (i /= ilo) then
                do j=1, ndim
                   psum(j) = 0.5*(p(i,j) + p(ilo,j))
                   p(i,j) = psum(j)
                enddo
                y(i) = funk(psum)
             endif
          enddo
          iter = iter + ndim
          !Keep trac of function evaluations.
          goto 1
       endif
    else
       iter = iter - 1 ! correct the evaluation count
    endif
    goto 2
  END SUBROUTINE amoeba
  !
  FUNCTION amotry (p, y, psum, mp, np, ndim, ihi, fac) RESULT(amo)
    !
    ! Extrapolates by a factor fac through the face of the simplex across from
    ! the high point, tries it, and replaces the high point if the new point
    ! is better
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN)           :: ihi, mp, ndim, np
    REAL,INTENT(IN)              :: fac
    REAL(r8_kind),INTENT(INOUT)  :: p(1:mp,1:np), psum(1:np), y(1:mp)
    REAL(r8_kind)                :: amo ! <<<result
    !*** end of interface ***
    INTEGER                      :: j
    REAL                         :: fac1, fac2
    REAL(r8_kind)                :: ytry, ptry(NMAX)
    !
    fac1 = (1.0 - fac)/ndim
    fac2 = fac1 - fac
    do j=1, ndim
       ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
    enddo
    ytry = funk(ptry)
    ! Evaluate the function at the trial point.
    if (ytry < y(ihi)) then ! If it's better then the highest,then replace it
       y(ihi) = ytry
       do j=1, ndim
          psum(j) = psum(j) - p(ihi,j) + ptry(j)
          p(ihi,j) = ptry(j)
       enddo
    endif
    amo = ytry
    return
  END FUNCTION amotry
  !
  FUNCTION funk(x) RESULT(f)
    !
    ! This is the function wich we would like to optimized
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(IN) :: x(1:ndim)
    REAL(r8_kind)             :: f !<<<result
    ! *** end of interface ***
    ! 
    !   if ( MINVAL(x) > 0 ) then
       call make_fit_exp(x,fit_set)    ! create fit exponents for input file
       call put_fit_in_input(fit_set)  ! write fit exponents in the input file
       call call_script                ! start the script of the optimizer
       call get_E(Ecoul, error)        ! get the result
       f = -Ecoul
       if (error == 1) then
          print *,"OPTIMIZER :: (function funk) error in ParaGauss execution !!!"
          write (unit=6,fmt=*) "OPTIMIZER :: (function funk) Error in ParaGauss execution !!!"
          stop
       end if
       ! else
       ! Optimizer algorithm give unphysical value ( smaller then zero )
       ! it make sense if we optimized exponents (not the logarithm of it!)
       !f = 10
       !endif
   !
  END FUNCTION funk
    !  
  SUBROUTINE read_input(ftol)
    !
    ! Read from opt_inp file input variables for q, n_exp, alpha_min, alpha_max for each type
    ! functions r^2, s, p, d ... And also read output_level and l_max
    !
    ! IMPORTANT !!!
    ! In this program as in the ParaGauss
    ! r^2 function -   lm = 0
    ! s   function -   lm = -1
    ! p   function -   lm = 1
    ! ...................
    ! REM: In some modules of the ParaGauss we have the reverse order 
    !      for r**2 and s functions :-) 
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(OUT) :: ftol
    ! *** end of interfece ***
    !
    CHARACTER(len=1)           :: type
    INTEGER                    :: i, lm, n_exp, k, status
    REAL(r8_kind)              :: q, alpha_min, alpha_max  
    !
    namelist/set_parameters/type, n_exp, q, alpha_min, alpha_max
    ! 
    open(unit=1, file="opt_inp",IOSTAT=status)
    if( status /= 0 ) print *,"open file error in read_input subroutine"
    read(unit=1,nml=control,IOSTAT=status)
    if( status /= 0 ) print *,"read error in read_input subroutine, namelist control" 
    ALLOCATE (var_set(-1:l_max))
    do i=-1, l_max
       read(unit=1,nml=set_parameters,IOSTAT=status)
       if( status /= 0 ) print *,"read error in read_input subroutine, namelist set_parameters"
       v => var_set(i)
       v%type=type
       select case (v%type)
       case ("r") 
          v%lm = 0
       case ("s")
          v%lm = -1
       case ("p")
          v%lm = 1
       case ("d")
          v%lm = 2
       case default
          v%lm = ICHAR(v%type(1:1)) - ICHAR("c")
       end select
       v%n_exp = n_exp 
!-10.09       
       if ( l_max == 0 .and. v%lm == -1 ) then
          if ( new_iter ) then
             v%n_exp = n
          else 
             v%n_exp = n_exp
          endif
       endif
!-10.09
       v%q         = q
       print *,"q=",q
       ! Now we use logarithm coordinate
       v%alpha_min = LOG(alpha_min)
       v%alpha_max = LOG(alpha_max)
    enddo
    close(unit=1)
    ! Now made an input vector x
    ! this vector iclude logarithm of the exponents
    v => var_set(0)
    if (v%n_exp /= 0) then
       ndim=(l_max+2)*3 ! number of dimensions with r**2 functions
    else
       ndim=(l_max+1)*3 ! number of dimensions without r**2 functions
    endif
    ALLOCATE (x(1:ndim))
    print *,"read_input >>>"
    k=0
    do i=-1, l_max
       v => var_set(i)
       if (v%n_exp /= 0) then
          x(i+2+k)=v%q
          x(i+3+k)=v%alpha_min
          x(i+4+k)=v%alpha_max
          print *,"x(",i+2+k,")", x(i+2+k)
          print *,"x(",i+3+k,")", EXP(x(i+3+k))
          print *,"x(",i+4+k,")", EXP(x(i+4+k))
          k=k+2
       endif
    enddo
    ALLOCATE (fit_set(-1:l_max))
    do i=-1, l_max
       v => var_set(i)
       f => fit_set(i)
       f%lm = v%lm
       f%type=v%type
       f%n_exp = v%n_exp
       if (f%n_exp /= 0) then
          f%q  = v%q
          ALLOCATE(f%expn(1:f%n_exp))
       endif
    enddo
    ! now we have the start vector x !!!!
    !
  END SUBROUTINE read_input
  !
  SUBROUTINE make_fit_exp( x, fit_set)
    !
    ! The subroutine make charge fit exponents for PG input file
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(IN)              :: x(1:ndim)
    TYPE(fit), TARGET, INTENT(OUT)         :: fit_set(-1:l_max)
    !*** end of interfece ***
    INTEGER                                :: status, output_level, k, j, N_max
    REAL(r8_kind)                          :: a_min, a_max
    !    
    TYPE(var), POINTER                     :: vv
    TYPE(fit), POINTER                     :: ff
    !
    ! trasfer from x in to var_set
    k=0
    do i=1, ndim, 3
       vv => var_set(i-2-k)
       vv%q         = x(i)
       vv%alpha_min = x(i+1)
       vv%alpha_max = x(i+2)
       k=k+2
    enddo
    print *,"make_fit_exp ----- >>>"
    do i=1, ndim, 3
       print *,"x(",i,")=",x(i)
       print *,"x(",i+1,")=",EXP(x(i+1))
       print *,"x(",i+2,")=",EXP(x(i+2))
    enddo
    do k=-1, l_max
       vv => var_set(k)
       ff => fit_set(k)
       !ff%q  = vv%q
       ff%lm = vv%lm
       ff%n_exp = vv%n_exp
       ff%type  = vv%type
       if (ff%n_exp == 2) then ! We have 2 exponents
          ff%q  = vv%q ! just dummy argument in this case
          ff%expn(1) = EXP(vv%alpha_min)
          ff%expn(2) = EXP(vv%alpha_max)
       else if (ff%n_exp == 1) then ! We have 1 exponent
          ff%expn(1) = EXP(vv%alpha_min)
       else if (ff%n_exp == 0) then ! We have no exponents
          ! Nothing, no exponents, no "q" factor.
          ! This is impliment for case without r2-functions
       else
          ! We have more than 2 exponents and we create a charge basis for PG input
          ff%q  = vv%q
          N_max = ff%n_exp-1
          do j=0, N_max
             a_min = EXP(vv%alpha_min)
             a_max = EXP(vv%alpha_max)
             !
             ff%expn(j+1)=(a_min*(a_max/a_min)**(real(j,kind=r8_kind)/real(N_max,kind=r8_kind))) &
                  *(vv%q**(real(j*(N_max-j),kind=r8_kind)/(real(2*N_max-2,kind=r8_kind))))
             !
             if (ff%expn(j+1) < 0) then
                print *,"OPTIMIZER  ::  OPTIMIZATION STOP !"
                print *,"Subroutine make_fit_exp: Error fit exponent smaller then zero"
                write (unit=6,fmt=*)"Subroutine make_fit_exp: Error fit exponent smaller then zero"
                stop
             endif
          enddo
       endif
    enddo
    !
    print *,"<<< ----- make_fit_exp" 
    !
  END SUBROUTINE make_fit_exp
  !
  SUBROUTINE put_fit_in_input(fit_set)
    !
    ! The subroutine write exponents in the file, then this file is used by PG
    ! input. In this file we have charge fit exponents
    !
    IMPLICIT NONE
    !
    TYPE(fit), INTENT(IN), TARGET          :: fit_set(-1:l_max)
    !
    !*** end of interfece ***
    INTEGER                                :: status, i
    REAL(r8_kind)                          :: a(5)
    TYPE(var), POINTER                     :: v
    !
    open(unit=2, file="ch_fit.op",IOSTAT=status)
    if( status /= 0 ) print *,"open file error in put_x subroutine"
    !r^2 - exponents
    f => fit_set(0)
    write(unit=2,fmt='(A   )') " &UNIQUE_ATOM_BASIS # unique atom  ? , charge fit basis set : r**2"
    write(unit=2,fmt='(A,I6)') "    N_EXPONENTS         =",f%n_exp
    write(unit=2,fmt='(A,I6)') "    N_UNCONTRACTED_FCTS =",f%n_exp
    write(unit=2,fmt='(A,I6)') "    N_CONTRACTED_FCTS   =     0 # (the default)"
    write(unit=2,fmt='(A,I6)') "    AUTOMATIC           = FALSE # (the default)"
    write(unit=2,fmt='(A,I6)') " /UNIQUE_ATOM_BASIS"
    if (f%n_exp /= 0) then
       write(unit=2,fmt=*) "  # exponents   q =",f%q
       write(unit=2,fmt='(3ES25.15:"  %")') (f%expn(j),j=1,f%n_exp)
    endif
    ! s - exponents
    f => fit_set(-1)
    write(unit=2,fmt='(    )')
    write(unit=2,fmt='(A   )') " &UNIQUE_ATOM_BASIS # unique atom  ? , charge fit basis set : l = 0"
    write(unit=2,fmt='(A,I6)') "    N_EXPONENTS         =",f%n_exp
    write(unit=2,fmt='(A,I6)') "    N_UNCONTRACTED_FCTS =",f%n_exp
    write(unit=2,fmt='(A,I6)') "    N_CONTRACTED_FCTS   =     0 # (the default)"
    write(unit=2,fmt='(A,I6)') "    AUTOMATIC           = FALSE # (the default)"
    write(unit=2,fmt='(A,I6)') " /UNIQUE_ATOM_BASIS"
    write(unit=2,fmt=*) "  # exponents   q =",f%q 
    write(unit=2,fmt='(3ES25.15:"  %")') (f%expn(j),j=1,f%n_exp)
    !p,d,f .. exponents
    do i=1, l_max
       f => fit_set(i)
       write(unit=2,fmt='(    )')
       write(unit=2,fmt='(A,I2)') " &UNIQUE_ATOM_BASIS # unique atom  ? , charge fit basis set : l = ",f%lm
       write(unit=2,fmt='(A,I6)') "    N_EXPONENTS         =",f%n_exp
       write(unit=2,fmt='(A,I6)') "    N_UNCONTRACTED_FCTS =",f%n_exp
       write(unit=2,fmt='(A,I6)') "    N_CONTRACTED_FCTS   =     0 # (the default)"
       write(unit=2,fmt='(A,I6)') "    AUTOMATIC           = FALSE # (the default)"
       write(unit=2,fmt='(A,I6)') " /UNIQUE_ATOM_BASIS"
       write(unit=2,fmt=*) "  # exponents   q =",f%q
       write(unit=2,fmt='(3ES25.15:"  %")') (f%expn(j),j=1,f%n_exp)
    enddo
    if (l_max == 0) then
       ! Standart fashion
       do i=1, 5
          a(i) = 0.1_r8_kind*(2.5_r8_kind**(i-1))
       end do
       write(unit=2,fmt='(    )')
       WRITE(unit=2,fmt='(A)')' &UNIQUE_ATOM_BASIS # unique atom  ? , charge fit basis set : l = 1'
       WRITE(unit=2,fmt='(A)')'    N_EXPONENTS         =     5'
       WRITE(unit=2,fmt='(A)')'    N_UNCONTRACTED_FCTS =     5'
       WRITE(unit=2,fmt='(A)')'    N_CONTRACTED_FCTS   =     0 # (the default)'
       WRITE(unit=2,fmt='(A)')'    AUTOMATIC           = FALSE # (the default)'
       WRITE(unit=2,fmt='(A)')' /UNIQUE_ATOM_BASIS'
       WRITE(unit=2,fmt='(A)')'  # exponents'
       WRITE(unit=2,'(3ES25.15:"  %")') (a(j),j=1,5)
       do i=1, 5
          a(i) = (0.1_r8_kind)*(2_r8_kind)*(2.5_r8_kind**(i-1))
          !d%alpha(i) = 0.1_r8_kind*(Real(d%lm,r8_kind))*(2.5_r8_kind**(i-1))
       end do
       write(unit=2,fmt='(    )')
       WRITE(unit=2,fmt='(A)')' &UNIQUE_ATOM_BASIS # unique atom  ? , charge fit basis set : l = 2'
       WRITE(unit=2,fmt='(A)')'    N_EXPONENTS         =     5'
       WRITE(unit=2,fmt='(A)')'    N_UNCONTRACTED_FCTS =     5'
       WRITE(unit=2,fmt='(A)')'    N_CONTRACTED_FCTS   =     0 # (the default)'
       WRITE(unit=2,fmt='(A)')'    AUTOMATIC           = FALSE # (the default)'
       WRITE(unit=2,fmt='(A)')' /UNIQUE_ATOM_BASIS'
       WRITE(unit=2,fmt='(A)')'  # exponents'
       WRITE(unit=2,'(3ES25.15:"  %")') (a(j),j=1,5)
    endif
    close(unit=2,IOSTAT=status)
    if( status /= 0 ) print *,"close file (ch_fit.op) error in put_x subroutine" 
    print *," put_fit_in_input << ---"
  END SUBROUTINE put_fit_in_input
  !
  SUBROUTINE call_script
    !
    ! The subroutine call the script for optimizer.
    ! This script start the ParaGauss program.
    ! If you would like to change the version of executable file
    ! of the PG program, you can make it in the script file "opt_script"
    !
    print *," ------ Script IN -------"
    call system("./opt_script")
    print *," ------ Script OUT ------"
    !
  END SUBROUTINE call_script
  !
  SUBROUTINE  store_fit
    call system("cp ch_fit.op ch_fit.store")
  END SUBROUTINE store_fit
  !
  SUBROUTINE get_E(Ecoul, error)
    !
    ! The subroutine read the result of PG calculation
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(OUT)       :: Ecoul
    INTEGER, INTENT(OUT)             :: error
    ! *** end of interfece ***
    CHARACTER(LEN=70)                :: input_line
    INTEGER                          :: status
    !
    open(unit=3, file="opt_helpE",IOSTAT=status)
    if( status /= 0 ) print *,"open file error in get_x subroutine"
    read(unit=3,fmt='(10X,F30.15)') Ecoul
    print *,"E=",Ecoul
    close(unit=3)
    open(unit=4, file="opt_error",IOSTAT=status)
    if( status /= 0 ) print *,"open file error in get_x subroutine"
    read(unit=4,fmt='(A)',IOSTAT=status) input_line
    if( status < 0 ) then
       error=0
    else
       if (input_line(1:20)=="                    ") then 
          error = 0
       else
          error = 1
       endif
    endif
    close(unit=4)
    print *,"error =", error
    print *,Ecoul,"<<<< Ecoul <<<<<"
  END SUBROUTINE get_E
  !
  SUBROUTINE from_x_to_var(x)
    !
    ! The subroutine converte the data from x-type (used by optimizer algorithm)
    ! to var-type (used by input and output)
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(IN)              :: x(1:ndim)
    ! *** end of interfece ***
!1/07/99 ->
    !TYPE(var), TARGET                      :: var_set(-1:l_max)
!1/07/99 <-
    ! 
    INTEGER                                :: i, k
    TYPE(var), POINTER                     :: vv
    ! trasfer from x in to var_set
    k=0
    do i=1, ndim, 3
       vv => var_set(i-2-k)
       vv%q         = x(i)
       vv%alpha_min = x(i+1)
       vv%alpha_max = x(i+2)
       k=k+2
    enddo
    ! Print the result of optimization
    do i=-1, l_max
       vv => var_set(i)
       vv%lm = i
       select case (vv%lm)
       case (-1) 
          vv%type = "s"
       case (0)
          vv%type = "r"
       case (1)
          vv%type = "p"
       case (2)
          vv%type = "d"
       case default
          vv%type = CHAR(ICHAR("c")+i)
       end select
       !       
       if (output_level > 0) then 
          write (unit=6,fmt=*) " ",vv%type," type functions"
          write (unit=6,fmt=*) "n_exp=", vv%n_exp
          if (vv%n_exp /= 0) then
             write (unit=6,fmt=*) "q=",vv%q
             write (unit=6,fmt=*) "alpha_min =", EXP(vv%alpha_min)
             write (unit=6,fmt=*) "alpha_max =", EXP(vv%alpha_max)
          endif
          write (unit=6,fmt=*) " "
       endif
    enddo
    !
  END SUBROUTINE from_x_to_var
  !
  SUBROUTINE from_x_to_var_result_print (x)
    !
    ! The subroutine write the result of optimization in "opt_output"
    !
    IMPLICIT NONE
    !
    REAL(r8_kind), INTENT(IN)              :: x(1:ndim)
    ! *** end of interfece ***
    !
    write (unit=6,fmt=*) " "
    write (unit=6,fmt=*) "       THE PARAMETERS AFTER OPTIMIZATION "
    write (unit=6,fmt=*) "  (optimised charge fit is in the ch_fit.op file)"
    write (unit=6,fmt=*) " "
    !
    output_level=1
    call from_x_to_var (x)
    !
    write (unit=6,fmt=*) " ------------------------------------------------- "
  END SUBROUTINE from_x_to_var_result_print
  !
 
  !
END PROGRAM fitfct_optimizer 



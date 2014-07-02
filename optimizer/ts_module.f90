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
module  ts_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Contains all routines necessary to locate
  !           transition states.
  !
  !  References: Jon Baker, An algorithm for the Location of 
  !              Transition States, J.Comp.Chem. 7, 385 (1985)
  !
  !  Author: FN
  !  Date: 3/98
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
# include "def.h"
  use type_module ! type specification parameters
  use iounitadmin_module
  use opt_data_module, only: OPT_STDOUT
  use hesse_module, only: step_counter
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of constants and variables ---------------
  !VVP
  real(kind=r8_kind),allocatable       :: step(:)
  real(kind=r8_kind),allocatable       :: ahesse(:,:),scl(:,:),help(:,:)
  real(kind=r8_kind),allocatable       :: ahesse_eigvec(:,:),ahesse_eigval(:)
  real(kind=r8_kind)                   :: r_curr, step_length
  real(kind=r8_kind)                   :: e_prev                 ! SAVE by default
  real(kind=r8_kind),allocatable       :: q_save(:),grad_save(:) ! SAVE by default
  logical                              :: less_t_r_curr=.false.
  !END OF MY VARIABLES
  real(kind=r8_kind),allocatable  :: f_vec(:),lambda_vec(:)
  real(kind=r8_kind)              :: lambda_p,lambda_n,alpha
  real(kind=r8_kind),save          :: de_mod,de_obj,de_mod_prev
  integer(kind=i4_kind):: io_eigmod=40
!  integer(kind=i4_kind)           :: step_counter=0                   ! SAVE by default
  logical                         :: accepted
  !------------ public functions and subroutines ---------------------
  public ts_main, ts_module_persistent_state
  
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !------------ Subroutines ------------------------------------------
contains
  
  
  subroutine ts_main(geo_loop)
    !  Purpose: main routine to control location of transition 
    !           states or minima.
    use hesse_module
    use math_module
    use gradient_module, only: grad_intern,energy
    use opt_data_module
    use coordinates_module
    use matrix_eigenval
    use step_module, only: step_max_comp, step_mean_square
    USE_DEBUG
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in) :: geo_loop
    !------------ Declaration of local variables ---------------------
    !VVP

    integer(kind=i4_kind)            :: n_negative,j
    integer(kind=i4_kind)            :: i,mode,tv_mode,num_n_eig
    character(len=80)                :: method_w
    logical                          :: recalc,err,desired
    real(kind=r8_kind)               :: a0_rfo,denum,a0_ah,e_cur
    real(kind=r8_kind)               :: lambda_rs_i_rfo,lambda_ah,lambda_rs_rfo,tht,ad
    real(kind=r8_kind),parameter     :: scale1=0.6_r8_kind,small_eig=0.0001_r8_kind
    real(kind=r8_kind),parameter     :: eps1=0.00005,eps2=0.0000005
    real(kind=r8_kind),save          :: e_save
    real(kind=r8_kind)               :: lambda_ahrs,a0_ahrs,mu
    real(kind=r8_kind)               :: step_internal(n_internal)
    !------------ Executable code ------------------------------------
    e_cur=energy
    ! it is usefull alwais to find energy in flepo output 
    write(io_flepo,*) step_counter,e_cur,e_prev, 'step_counter, e_cur,e_prev 1' 

    err=.false.
    ! first determine local characteristics of Hessian
    if(wales) then 
      method_w="wales"
    elseif(rfo_step) then
      method_w="baker"
    elseif(eigenmode_follow) then
      method_w="EF" 

    elseif(dynamic.and.method_rfo.and.ts_search) then
      method_w="rs_i_rfo"
    elseif(dynamic.and.method_ah.and.ts_search) then
      method_w="ahrs"
    else
      call error_handler("Please check settings in optimizer input")
    endif

    call step_setup(method_w)
    
    meth: if((method_w.ne."rs_i_rfo").and.(method_w.ne."ahrs")) then

      f_vec=zero
      n_negative=0_i4_kind
      do i=1,n_internal
         if(hesse_eigval(i)<zero)  n_negative=n_negative+1
         f_vec(i) = dot_product(hesse_eigvec(:,i),grad_intern)
      enddo
      alpha=two ! Wales-Setting: although he uses this as a parameter,
                ! it is always set to two.
      if(eigenmode_follow) then
        call eigenmode_identify(mode,err)
        if(err) then
          mode = max(1,n_negative)
        endif
      else
        if(n_negative>1) then
          mode = 2
       else
          mode = 1
       endif
      endif

      lambda_vec = zero
      if((.not.minimization.and.n_negative==0) .or.&
        (n_negative>0.and.rfo_step)) then
        lambda_p = half*hesse_eigval(mode) + half*sqrt(hesse_eigval(mode)**2 + &
                   alpha**2*f_vec(mode)**2)
        call newton(lambda_vec,mode)   ! (1)
        lambda_n = minval(lambda_vec)
! play safe: not only must minval(lambda_vec) be nagetive,but it also must
! be lower than hesse_eigval(2) (if the lowest eigenmode is being followed)
! or lower than hesse_eigval(1) ( if some other mode is being followed):
        recalc=.false.
        if(mode == 1) then
          if(lambda_n>hesse_eigval(2)) recalc=.true.
        else
          if(lambda_n>hesse_eigval(1)) recalc=.true.
        endif
        if(recalc) then
          write(OPT_STDOUT,*)" ts_main: calculation of shift parameters using "
          write(OPT_STDOUT,*)"          J. Bakers method (J.Comp.Chem. 7,385-396,(1986))"
          write(OPT_STDOUT,*)"          failed. Re-calculate them using the method of"
          write(OPT_STDOUT,*)"          Wales (J.Chem.Soc.FARADAY TRANS.,1990,86, 3505-3517)"
          lambda_n = half*hesse_eigval(mode) - half*sqrt(hesse_eigval(mode)**2 + &
          alpha**2*f_vec(mode)**2)  
        endif
        lambda_vec = lambda_n     
        lambda_vec(mode) = lambda_p

      elseif(wales) then
        do i=1,n_internal
           lambda_n = half*hesse_eigval(i) - half*sqrt(hesse_eigval(i)**2 + &
           alpha**2*f_vec(i)**2) 
           lambda_vec(i) = lambda_n
        enddo
        lambda_p = half*hesse_eigval(mode) + half*sqrt(hesse_eigval(mode)**2 + &
        alpha**2*f_vec(mode)**2)
        lambda_vec(mode) = lambda_p
      endif
       
    
      if(n_negative==0_i4_kind) then
        write(OPT_STDOUT,*)"ts_main : we are still in the vicinity of a"
        write(OPT_STDOUT,*)"          local minimum. Follow the lowest"
        write(OPT_STDOUT,*)"          eigenmode with eigenvalue"
        write(OPT_STDOUT,'("         ",f10.6)')hesse_eigval(mode)
      endif

      call ts_step(mode) !(1)
      step_length=dsqrt(dot_product(step,step))
      if(step_counter.ne.0) de_obj=e_cur-e_prev

      if(eigenmode_follow) then
        call eigenmode_store(mode)
      endif
      print*,'call eigenmode_store passed'

     write(io_flepo,*) 'Total energy:                            ', e_cur,e_prev,de_obj
     write(io_flepo,*) 'Step length:                             ', step_length 
     write(io_flepo,*) 'Predictable change of energy in Pade is: ', de_mod,de_mod_prev
     if(step_counter.ne.0) write(io_flepo,*) 'Ratio of actual and predicted energy change:', de_obj/de_mod_prev
     e_prev=e_cur

    elseif(method_w.eq."rs_i_rfo") then
!Working in RS-I-RFO approach

!Inicialization block. Define all local variables and set to initial values

     if(step_counter.eq.0) then
      de_mod_prev=0.0_r8_kind

     else
      de_obj=e_cur-e_prev
      de_mod_prev=de_mod
     endif

     de_mod=0.0_r8_kind
     desired=.false.
     num_n_eig=0
     tv_mode=1
     err=.false.
!TR model
     write(io_flepo,*) 'de_mod_prev', de_mod_prev

     call tr(geo_loop,de_mod_prev,less_t_r_curr,r_curr,e_cur,accepted)

     nacc: if(.not.accepted) then
     !Restore geometry and gradient from prevuois step
      write(io_flepo,'(A,7x,I10)')"Previous step was unacceptable. Using geometry from step:",geo_loop-1
      q=q_save
      grad_intern=grad_save
      if(less_t_r_curr.and.(r_curr>=step_length*scale1)) then
        r_curr=step_length*scale1
        write(io_flepo,'(A)')"Rescale Trust Radius[multiply by factor 0.6]"
      endif
     endif nacc

!Finding the index of the transition vector:
    write(io_flepo,'(A)')"Rational Function Optimization method Image Function method with Restricted Step "
     if(step_counter==0) then
      call eigenmode_store(tv_mode)
     else 
      call eigenmode_identify(tv_mode,err)
      if(err) then
         call error_handler("Transition vector was missed")
      else 
         call eigenmode_store(tv_mode)
      endif   
      write(io_flepo,'(A)')"Storing transition vector"
    endif
!End of looking for transition vector

!Usually, tv_mode==1, imho.
    write(io_flepo,'(A,4X,I4)')"The index of transition vector is:",tv_mode
    f_vec=zero
    do i=1,n_internal
       f_vec(i) = dot_product(hesse_eigvec(:,i),grad_intern)
    enddo
!At first pass we construct normal Augmented Hessian
    ahesse=zero
    ahesse(1,2:)=f_vec(:)
    ahesse(2:,1)=f_vec(:)
    do i=1,n_internal
       ahesse(i+1,i+1)=hesse_eigval(i)
    enddo
    ahesse(1,tv_mode+1)=-ahesse(1,tv_mode+1)
    ahesse(tv_mode+1,1)=-ahesse(tv_mode+1,1)
    ahesse(tv_mode+1,tv_mode+1)=-ahesse(tv_mode+1,tv_mode+1)
!End construct of Augmented Hessian

!Check inertia of current Hessian
    do i=1,n_internal 
       if(hesse_eigval(i)<=zero) num_n_eig=num_n_eig+1
    enddo
    if(num_n_eig==1) desired=.true.
    write(io_flepo,'(A,4X,I4)')"Number of negative eigenvalues:",num_n_eig

!!!!!!!!MAIN BLOCK!!!!!!!!    
!Solve secular equation for AH
    call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
    DCALL show("Augmented Hessian",ahesse(:,:))
    DCALL show("Eigenvectors of Augmented Hessian",ahesse_eigvec(:,:))
    DCALL show("Eigenvalues of Augmented Hessian",ahesse_eigval(:))
    a0_ah=ahesse_eigvec(1,1)
    lambda_ah=ahesse_eigval(1)
    step=ahesse_eigvec(2:,1)/a0_ah
    step=matmul(hesse_eigvec,step)
!Eliminate small components of step
    call eliminate_step(step)
    step_length=dsqrt(dot_product(step,step))
    step=matmul(tmat,step)
    write(io_flepo,*)""
    write(io_flepo,'(27X,A)')"--- Augmented Hessian ---"
    write (io_flepo,'(A,F15.8)')"a0_ah      =",abs(a0_ah)
    write (io_flepo,'(A,F15.8)')"lambda_ah  =",lambda_ah

    des: if(desired.and.(step_length<=r_curr)) then
      less_t_r_curr=.true.
      write (io_flepo,'(A)')"Making AH step"
      lambda_ah=lambda_ah+(two*f_vec(tv_mode)*step(tv_mode) &
       +hesse_eigval(tv_mode)*step(tv_mode)**2)/(one+one*dot_product(step,step))
      write(io_flepo,*) 'de_mod case 1'
#if 0
      de_mod=lambda_ah/two*a0_ah*a0_ah
#else
      de_mod=sum(step_internal*f_vec)  
#endif
      tht=ang(step)
      call write_step("AH",step,f_vec,hesse_eigval,zero,tht)
      q_save = q
      grad_save = grad_intern
    else des
!Step[2]
      write (io_flepo,*) 'At first pass we construct New Augmented Hessian'
      ahesse=zero
      ahesse_eigvec=zero
      ahesse_eigval=zero
      ahesse(1,2:)=f_vec(:)
      ahesse(2:,1)=f_vec(:)
      do i=1,n_internal
         ahesse(i+1,i+1)=hesse_eigval(i)
      enddo
      ahesse(1,tv_mode+1)=-ahesse(1,tv_mode+1)
      ahesse(tv_mode+1,1)=-ahesse(tv_mode+1,1)
      ahesse(tv_mode+1,tv_mode+1)=-ahesse(tv_mode+1,tv_mode+1)
      less_t_r_curr=.false.
      lambda_rs_i_rfo=lambda_ah

      write(io_flepo,'(A,"calculated length=",f20.8,", maximum value=",f15.8)')&
      "AH method: step length exceeds maximum value: ",step_length,r_curr

      write(io_flepo,'(27X,A)')"--- Solve RS-RFO equations ---"
!Start of microiterative process
      DPRINT "===================="
      DPRINT "Compute     Alpha"

      alpha=one
      j=1
      DPRINT "Iter	   Alpha",j,alpha
      do while((abs(step_length-r_curr)>eps1).or.(abs(ad)>eps2))
!Calculate d(step_length)/dalpha
              j=j+1
              ad=delta(lambda_rs_i_rfo,alpha)
              alpha=alpha+ad

 DPRINT 'Solve general secular equation',j,denum,r_curr,step_length,ad,alpha
              call new_step(alpha,step_length,lambda_rs_i_rfo)
!              ad=delta(lambda_rs_i_rfo,alpha)
!              alpha=alpha+0.01
!              call new_step(alpha,step_length,lambda_rs_i_rfo)
!stop
      enddo
     

      DCALL show("New Augmented Hessian[RS-I-RFO]",ahesse(:,:))
      DCALL show("New Eigenvectors of Augmented Hessian[RS-I-RFO]",ahesse_eigvec(:,:))
      DCALL show("New Eigenvalues of Augmented Hessian[RS-I-RFO]",ahesse_eigval(:))
!End of microiterative process

!So, we have solution of RS-I-RFO equations
!Recalculation lambda:
      lambda_rs_rfo=lambda_rs_i_rfo+(two*f_vec(tv_mode)*step(tv_mode) &
      +hesse_eigval(tv_mode)*step(tv_mode)**2)/(one+alpha*dot_product(step,step))
!End of recalculation lambda

!Bakwards transformation of step.
      step_internal=step
      step=matmul(hesse_eigvec,step)
      step=matmul(tmat,step)
!Eliminate small components of step
      call eliminate_step(step)
      step_length=dsqrt(dot_product(step,step))
      write(io_flepo,*) 'de_mod case 2'
      de_mod=lambda_rs_rfo*(one+alpha*step_length**2)/two
#if 1
      de_mod=sum(step_internal*f_vec)
#else
      de_mod=lambda_rs_rfo*(one+alpha*step_length**2)/two
#endif
      write (io_flepo,'(A,F15.8)')"a0_rfo     =",abs(a0_rfo)
      write (io_flepo,'(A,F15.8)')"lambda_rs_rfo =",lambda_rs_rfo
      write (io_flepo,'(A,F15.8)')"alpha      =",alpha
      write (io_flepo,'(A)')"Making RS-I-RFO step"
      tht=ang(step)
      call write_step("RS-I-RFO",step,f_vec,hesse_eigval,zero,tht,step_internal)
      grad_save=grad_intern
      q_save=q
      e_save=e_cur
      endif  des

   else meth
   write(io_flepo,*) 'use augmented hessian method'
    
!#
    call tr(geo_loop,de_mod,less_t_r_curr,r_curr,e_cur,accepted)
    nacc2: if(.not.accepted) then
!Restore geometry and gradient from prevuois step
      write(io_flepo,'(A,7x,I10)')"Prevuoisly step was unacceptable.Using geometry from step:",geo_loop-1
      q=q_save
      grad_intern=grad_save
      if(less_t_r_curr.and.(r_curr>=step_length*scale1)) then
        r_curr=step_length*scale1
        write(io_flepo,'(A)')"Rescale Trust Radius[multiply by factor 0.6]"
      endif
    endif   nacc2
     f_vec=zero
     do i=1,n_internal
        f_vec(i) = dot_product(hesse_eigvec(:,i),grad_intern)
     enddo               
     write(io_flepo,'(A)')"Using Augmented hessian technique"
!Construct Augmented Hessian
     !call constr_ah(zero,zero,zero,hesse_eigval,f_vec,ahesse)
     ahesse=zero
     ahesse(1,2:)=f_vec(:)
     ahesse(2:,1)=f_vec(:)
     do i=1,n_internal
       ahesse(i+1,i+1)=hesse_eigval(i)
     enddo
!Solve secular equation for AH
     call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
!     call eigensolver(ahesse,n_internal+1,ahesse_eigval,ahesse_eigvec)
     DCALL show("Augmented Hessian",ahesse(:,:))
     DCALL show("Eigenvectors of Augmented Hessian",ahesse_eigvec(:,:))
     DCALL show("Eigenvalues of Augmented Hessian",ahesse_eigval(:))
!Construct unitary transformation matrix
     help=zero
     help(1,1)=one
     help(2:,2:)=hesse_eigvec
!     call show("Unitary transformation",help(:,:))
!     call show("hesse eigvec",hesse_eigvec(:,:))
     ahesse_eigvec=matmul(transpose(help),ahesse_eigvec)
!Compute a0
     write (io_flepo,'(A)') "Augmented Hessian"
     a0_ah=ahesse_eigvec(1,2)
     step=ahesse_eigvec(2:,2)/a0_ah
     step_internal=step
     DCALL show("second row",ahesse_eigvec(:,2))
     DCALL show("first row",ahesse_eigvec(:,1))
     print *,"a0 from second row",ahesse_eigvec(1,2)
!Bakwards transformation of step.
     step=matmul(hesse_eigvec,step)
     step=matmul(tmat,step)
!Eliminate small components of step
     call eliminate_step(step)
     step_length=dsqrt(dot_product(step,step))
     lambda_ah=ahesse_eigval(2)
     write (io_flepo,'(A,F15.8)') "a0_ah=",abs(a0_ah)
     do i=1,n_internal 
       if(hesse_eigval(i)<=zero) num_n_eig=num_n_eig+1
     enddo
     if(num_n_eig==1) desired=.true.
!call desrd("AH",n_negative,desired)
     if(abs(a0_ah)<=0.75_r8_kind) then
     write(io_flepo,*) "a0<0.75. So current point is far from solution"
     write(io_flepo,*)"Using method suggested by Khait"
     write(io_flepo,*)"Follow the lowest"
        write(io_flepo,*)"           eigenmode with eigenvalue"
        write(io_flepo,'("         ",f10.6)')hesse_eigval(1)
        step=hesse_eigvec(:,1)*f_vec(1)/hesse_eigval(1)
        step_internal=0.0_r8_kind
        step_internal(1)=f_vec(1)/hesse_eigval(1)
        if(step_counter.ne.1) de_mod_prev=de_mod
!        de_mod=-half*f_vec(1)**2/hesse_eigval(1)
        de_mod=half*f_vec(1)**2/hesse_eigval(1)
        if(step_counter.ne.0) then 
         de_obj=e_cur-e_prev
         write(io_flepo,*) "e_cur-e_prev",e_cur,e_prev
        endif
        step = matmul(tmat,step)
        step_length=dsqrt(dot_product(step,step))
        write(io_flepo,*)"Perform scale of calculated step without change the direction"
        step_length=dsqrt(dot_product(step,step))
        tht=ang(step)
        call write_step("AH-Khait",step,f_vec,hesse_eigval,zero,tht,step_internal)
        step = step/step_length*step_max
        step_internal=step_internal/step_length*step_max
        step_length=dsqrt(dot_product(step,step))
     else
!Step[1]
      write(io_flepo,'(A)')"a0>=0.75.Current point is sufficiently close to the solution"
      write(io_flepo,'(A)')"Using Augmented Hessian technique"
      write(io_flepo,'(A,F15.8)')"lambda:",lambda_ah
      do i=1,n_internal 
         if(hesse_eigval(i)-lambda_ah<=zero) num_n_eig=num_n_eig+1
      enddo
      if(num_n_eig==1) desired=.true.
!      call inertia(zero,lambda_ah,n_negative)
!      call desrd("AH",n_negative,desired)
      des3: if(desired.and.(step_length<=r_curr)) then
        less_t_r_curr=.true.
        write (io_flepo,'(A)')"Making AH step"
        write (io_flepo,'(A,F15.8)')"lambda_ah=",lambda_ah
        de_mod_prev=de_mod
        de_mod=lambda_ah/2*a0_ah*a0_ah
        tht=ang(step)
        call write_step("AH",step,f_vec,hesse_eigval,zero,tht,step_internal)
      else des3
!Step[2]	
          less_t_r_curr=.false.
          ahesse_eigvec=zero
          ahesse_eigval=zero
          write(io_flepo,*)"AH method: step length exceeds maximum value "
          write(io_flepo,'("         calculated length ",f20.5," maximum value ",f9.5)')step_length,&
          r_curr
          write (io_flepo,'(A)')"Making AHRS step"
          write(io_flepo,'(A)')"Compute mu. Solve  Augmented Hessian with Restricted Step equation [AHRS]"
          call hebden(mu,zero,r_curr)
          write(io_flepo,*)"mu=",mu
!          call constr_ah(mu,lambda_ah,r_curr,hesse_eigval,f_vec,ahesse)
          ahesse=zero
          do i=1,n_internal
             ahesse(i+1,i+1)=hesse_eigval(i)
          enddo
          ahesse(1,2:)=f_vec(:)
          ahesse(2:,1)=f_vec(:)
          ahesse(1,1)=-mu*r_curr**2
          call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
          DCALL show("New Augmented Hessian[AHRS]",ahesse(:,:))
          DCALL show("New Eigenvectors of Augmented Hessian[AHRS]",ahesse_eigvec(:,:))
          DCALL show("New Eigenvalues of Augmented Hessian[AHRS]",ahesse_eigval(:))
          a0_ahrs=ahesse_eigvec(1,2)
          write(io_flepo,*)"a0_ahrs=",abs(a0_ahrs)
          step=ahesse_eigvec(2:,2)/a0_ahrs
          lambda_ahrs=ahesse_eigval(2)
!Bakwards transformation of step.
          step_internal=step
          step=matmul(hesse_eigvec,step)
          step=matmul(tmat,step)
!Eliminate small components of step
          call eliminate_step(step)
          step=matmul(tmat,step)
          de_mod_prev=de_mod
          de_mod=lambda_ahrs/2*a0_ahrs*a0_ahrs
          step_length=dsqrt(dot_product(step,step))
          write (io_flepo,*) "lambda_ahrs=",lambda_ahrs
          tht=ang(step)
          call write_step("AHRS",step,f_vec,hesse_eigval,zero,tht,step_internal)
    endif des3
    endif
    q_save = q
    grad_save = grad_intern
   endif meth

     q = q + step ! update q
     step_max_comp = maxval(step)
     step_mean_square = sqrt(sum(step**2)/n_internal)

   step_counter=step_counter+1
   call step_shutdown(method_w)
  contains 
              subroutine new_step(alpha,step_length,lambda_rs_i_rfo)
              real(r8_kind), intent(in):: alpha
              real(r8_kind), intent(out):: step_length,lambda_rs_i_rfo
              scl=zero
              scl(1,1)=one
              do i=1,n_internal
                scl(i+1,i+1)=alpha
              enddo
!    CALL show("Augmented Hessian",ahesse(:,:))
!    CALL show("scl",scl(:,:))
              call geigs(ahesse,scl,ahesse_eigval,ahesse_eigvec)
!    DPRINT 'done geigs'
              a0_rfo=ahesse_eigvec(1,1)
              lambda_rs_i_rfo=ahesse_eigval(1)
              step=ahesse_eigvec(2:,1)/a0_rfo
              step_length=dsqrt(dot_product(step,step))
              write (*,'(A,F15.8)')"current step_length=",step_length
              DCALL show("Eigenvalues of AH",ahesse_eigval(:))
        end subroutine new_step


              function delta(lambda_rs_i_rfo,alpha)
              real(kind=r8_kind):: lambda_rs_i_rfo,delta,alpha
              denum=zero
              do i=1,n_internal
                 denum=denum+f_vec(i)**2/(hesse_eigval(i)-lambda_rs_i_rfo*alpha)**3
              enddo
              denum=denum*lambda_rs_i_rfo/(one+alpha*step_length**2)
              alpha=alpha+0.01_r8_kind
              call new_step(alpha,step_length,lambda_rs_i_rfo)
              denum=step_length**2
              alpha=alpha-0.01_r8_kind
              call new_step(alpha,step_length,lambda_rs_i_rfo)
              denum=50.0_r8_kind*(denum-step_length**2)
!End of calculate of d(step_length)/dalpha
              delta=(r_curr*step_length-step_length**2)/denum
              end function delta

  end subroutine ts_main
  !*************************************************************

  subroutine ts_step(mode_index)
    ! Purpse :preliminary, make a step that maximizes in the direction
    !         of mode mode_index and minimizes in all other 
    !         directions
    ! --------------------------------------------------------
    use hesse_module, only: hesse_eigvec,hesse_eigval
    use coordinates_module, only: tmat
    use math_module, only: zero,one,abs_value
    use opt_data_module, only: step_max, n_internal
    implicit none
    integer(kind=i4_kind),intent(in)  :: mode_index
    ! -------------------------------------------------------
    integer(kind=i4_kind)             :: i
    real(kind=r8_kind)                :: lambda,scaling_fac
    
    step = zero
    lambda=zero
    if(step_counter.ne.0) de_mod_prev=de_mod
    de_mod=zero
    do i=1,n_internal
       if (i==mode_index) then
          lambda=lambda_p
       else
          lambda = lambda_vec(i)
       endif
       step = step + (-one)*f_vec(i)*hesse_eigvec(:,i)/(hesse_eigval(i) - lambda_vec(i))
       de_mod=de_mod-0.5_r8_kind*f_vec(i)**2/(hesse_eigval(i) - lambda_vec(i))
    enddo
    step = matmul(tmat,step)
    write(OPT_STDOUT,*)"step",step
    ! where (.not.s(:)%var ) step=zero
    if (abs_value(step)>step_max) then
       write(OPT_STDOUT,*)"ts_step: step length exceeds maximum value "
       write(OPT_STDOUT,'("         calculated length ",f19.5," maximum value ",f9.5)') abs_value(step),&
            step_max 
       write(OPT_STDOUT,*)"         scaling step down to maximum length"      
       scaling_fac=step_max/abs_value(step)
       de_mod=de_mod*scaling_fac*(2.0_r8_kind-scaling_fac)
       step = step*scaling_fac
    endif
    write(OPT_STDOUT,*)"                          --- TS_STEP --- "
    write(OPT_STDOUT,*)"   mode    grad. comp.       shift parameter     eigenvalue     coord       STEP "
    do i=1,n_internal
       write(OPT_STDOUT,'(3x,i3,a3,es12.5,10x,es12.5,6x,es12.5,4x,i3,6x,f9.5,2f13.7)')&
            i,'m  ',f_vec(i),lambda_vec(i),hesse_eigval(i),i,step(i), &
            -f_vec(i)/(hesse_eigval(i) - lambda_vec(i))*scaling_fac, &
            -0.5_r8_kind*f_vec(i)**2/(hesse_eigval(i) - lambda_vec(i))*scaling_fac 
    enddo
    write(OPT_STDOUT,*)"                          --------------- "
  end subroutine ts_step
  !*************************************************************
  subroutine eigenmode_store(index_mode)
    ! Purpose: stores the eigenmode followed uphill (with or without
    !          shift parameter is not important here) on the
    !          file 'eigenmode.dat'
    ! --------------------------------------------------------
    use hesse_module, only: hesse_eigvec
    use opt_data_module, only: print_debug,opt_data_dir,n_internal
    implicit none
    integer(kind=i4_kind)    :: index_mode
    ! --------------------------------------------------------
    io_eigmod=openget_iounit(form='formatted',status='replace',&
         file=trim(opt_data_dir)//'/eigenmode.dat')
    if (print_debug) then
       write(OPT_STDOUT,*)" eigenmode_store: writing eigenmode No.",index_mode,&
            "to file eigenmode.dat"
    endif
    write(io_eigmod,*)index_mode,n_internal
    write(io_eigmod,*)hesse_eigvec(:,index_mode)
    call returnclose_iounit(io_eigmod)
  end subroutine eigenmode_store
  !*************************************************************
  subroutine eigenmode_identify(index_mode,err)
    ! Purpose:test if the file 'eigenmode.dat' is present. If not,
    !         issue a warning but continue to follow lowest eigenmode,
    !         i. e. set 'index_mode' to one on output.
    !         If the file is present, read from it the index of the mode
    !         followed in the previous cycle, the dimension of the 
    !         eigenvector (n_internal) and the eigenvector itself.
    !         calculate then all dot-products of this vector
    !         with all eigenvectors of the actual Hessian and see,
    !         which one corresponds to the largest overlap.
    ! --------------------------------------------------------
    use hesse_module, only: hesse_eigvec
    use math_module, only: zero
    use opt_data_module, only: opt_data_dir,n_internal
    implicit none
    integer(kind=i4_kind)      :: index_mode
    logical                    :: err
    ! --- declaration of local variables ---------------------
    logical :: exist
    real(kind=r8_kind),allocatable  :: last_mode(:),prod(:)
    integer(kind=i4_kind)           :: n_int, index_old,&
         alloc_stat,index_largest,i
    real(kind=r8_kind)              :: max_val
    ! --- executable code ------------------------------------
    exist = .false.
    err=.false.
    inquire(EXIST=exist,file=trim(opt_data_dir)//'/eigenmode.dat')
    if (.not.exist) then
       write(OPT_STDOUT,*)" "
       write(OPT_STDOUT,*)"    --- Eigenmode_Identify: WARNING ---"
       write(OPT_STDOUT,*)" The file eigenmode.dat does not exist "
       write(OPT_STDOUT,*)" This should only occur in the first geometry cycle"
       write(OPT_STDOUT,*)" "
       err=.true.
       return
    endif
    io_eigmod=openget_iounit(form='formatted',status='old',&
         file=trim(opt_data_dir)//'/eigenmode.dat')
    read(io_eigmod,*)index_old,n_int
    if (n_int /= n_internal ) call error_handler&
         ("eigenmode_identify: dimensions in file wrong")

    allocate(last_mode(n_internal),&
         prod(n_internal),&
         STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" eigenmode_identify: allocation (1) failed")
    prod = zero
    read(io_eigmod,*)last_mode(:)

    index_largest = 0
    max_val = zero
    do i=1,n_internal
       prod(i) = abs(dot_product(last_mode,hesse_eigvec(:,i)))
       if (prod(i)>max_val) then
          max_val = prod(i)
          index_largest = i
       endif
    enddo
       
    write(OPT_STDOUT,*)" eigenmode_identify: list of overlaps with eigenvector ",index_old,&
         " from last cycle"
    do i=1,n_internal
       write(OPT_STDOUT,'(3x,f9.5)')prod(i)
    enddo
    write(OPT_STDOUT,*)" "
    write(OPT_STDOUT,*)" eigenmode_identify: -------------------------"
    write(OPT_STDOUT,*)"   index of mode followed in previous cycle:",index_old
    write(OPT_STDOUT,*)"   index of mode with largest overlap      :",index_largest
    write(OPT_STDOUT,*)"   Eigenmode being followed now:"
    do i=1,n_internal
       write(OPT_STDOUT,'(3x,f9.5)')hesse_eigvec(i,index_largest)
    enddo
    write(OPT_STDOUT,*)" "
    index_mode = index_largest
    call returnclose_iounit(io_eigmod)
  end subroutine eigenmode_identify

  subroutine step_setup(meth)
    !-------------------------------------------------------------
    ! Purpose: depending on the chosen method (Quasi-Newton,
    !          rationale function/ augmented hessian etc.) 
    !          all required variablea are allocated and initialized
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use opt_data_module, only: n_internal
    !------------ Declaration of formal parameters ---------------
    character(len=*),    intent(in)  :: meth
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)            :: alloc_stat
    !------------ Executable code --------------------------------  

    allocate(step(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    if(.not.allocated(q_save)) then
    allocate(q_save(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    endif     
    if(.not.allocated(grad_save)) then
    allocate(grad_save(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    endif 
    allocate(f_vec(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    allocate(help(n_internal+1,n_internal+1),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    if((meth.ne."rs_i_rfo").and.(meth.ne."ahrs")) then
      allocate (lambda_vec(n_internal),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
    else   
      allocate(ahesse_eigval(n_internal+1),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      allocate(ahesse_eigvec(n_internal+1,n_internal+1),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      allocate(ahesse(n_internal+1,n_internal+1),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      if(meth.eq."rs_i_rfo") then
        allocate(scl(n_internal+1,n_internal+1),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
      endif
    endif
 
 end subroutine step_setup

 
 
 subroutine step_shutdown(mth)
     !-------------------------------------------------------------
     ! Purpose: all used array are deallocating
     !-------------------------------------------------------------
     !------------ Declaration of formal parameters ---------------
     character(len=*),    intent(in)      :: mth
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind)   :: alloc_stat
     !------------ Executable code --------------------------------
     deallocate(step,STAT=alloc_stat)
     ASSERT(alloc_stat==0)
     deallocate(f_vec,STAT=alloc_stat)
     ASSERT(alloc_stat==0)
     deallocate(help,STAT=alloc_stat)
     ASSERT(alloc_stat==0)
     if((mth.eq."rs_i_rfo").or.(mth.eq."ahrs")) then
       deallocate(ahesse_eigval,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(ahesse_eigvec,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       deallocate(ahesse,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       if(mth.eq."rs_i_rfo") then
       deallocate(scl,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       endif
     else
     deallocate(lambda_vec,STAT=alloc_stat)
     ASSERT(alloc_stat==0)
     endif
 end subroutine step_shutdown

  subroutine ts_module_persistent_state(act,iounit)
    !
    ! Save/Restore module vars that need to persist
    ! across optimizer calls if optimizer is compiled
    ! as a standalone program or for "optimizer_only" runs
    !
    implicit none
    character(len=*), intent(in) :: act
    integer(i4_kind), intent(in) :: iounit
    ! *** end of interface ***

    integer(i4_kind) :: n_coords, allocstat

    ! global module vars:
    namelist /ts_module_state/ de_mod      &
                               , e_prev & 
                               , step_length &
                               , less_t_r_curr

    !
    ! Gfortran adds an empty line after write(iou,nml=...)
    ! so that it is difficult to put comments like section
    ! headers into the "optimizer.state" file
    !
    select case(act)

    case('restore')

      ! 4) random vars:
      read(iounit,nml=ts_module_state)
    !
    ! Gfortran adds an empty line after write(iou,nml=...)
    ! so that it is difficult to put comments like section
    ! headers into the "optimizer.state" file
    !
      ! 5) q_save, grad_save (these are used to back off to
      !    the previous geometry, if the new is too bad):
      read(iounit,*) n_coords
      if( n_coords >= 0 )then
        if( .not.allocated(q_save) )then
          allocate(q_save(n_coords),grad_save(n_coords), stat = allocstat)
          ASSERT(allocstat==0)
        endif
        read(iounit,*) q_save(:)
        read(iounit,*) grad_save(:)
      endif

    case('save')

      ! 4) random vars:
      write(iounit,nml=ts_module_state)

      ! 5) q_save, grad_save:
      n_coords = -1

     if( allocated(q_save) ) n_coords = size(q_save)

      write(iounit,*) n_coords
      if( n_coords >= 0 )then
        write(iounit,*) q_save(:)
        write(iounit,*) grad_save(:)
      endif

    case default
      ABORT('no such action')
    end select
  end subroutine ts_module_persistent_state

 subroutine tr(geo_loop,de_mod,less_t_r_curr,r_curr,e_cur,accepted)
    !-------------------------------------------------------------
    !  Purpose: Calculate Trust Radius at each iteration.
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use math_module,only:two,zero
    use gradient_module,only:energy,grad_mean_square,grad_max_comp
    use opt_data_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in)  :: geo_loop
    real(kind=r8_kind),   intent(in)  :: de_mod
    logical,              intent(in)  :: less_t_r_curr
    real(kind=r8_kind),   intent(out) :: r_curr,e_cur
    logical,              intent(out) :: accepted
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)                :: r,de_obj,e_curr
    real(kind=r8_kind)                :: r1e,r1i,rue,rui
    real(kind=r8_kind),parameter      :: r_min=0.75_r8_kind,r_good=0.85_r8_kind
    real(kind=r8_kind),parameter      :: r_0=0.3_r8_kind,sf=two,gr_small=0.001_r8_kind
    real(kind=r8_kind),parameter      :: en_small=0.000015_r8_kind
    !------------ Executable code --------------------------------  
    write(io_flepo,*)"--- Trust Radius --- for step", step_counter
    if((grad_max_comp<=gr_small).and.(grad_mean_square<=gr_small).and.(step_counter/=0)) then
      write(io_flepo,2200)
      r_curr=r_curr
      e_cur=energy
      write(io_flepo,2100) r_curr
        de_obj=e_curr-e_prev
        if(abs(de_mod).gt.0.00000001_r8_kind)r=de_obj/de_mod
      write(io_flepo,1800) r
      accepted=.true.
    else
!r1e=0.75;r1i=1.25;rue=1.15;rui=0.85    
      r1e=r_min;r1i=2-r_min;rue=2-r_good;rui=r_good
      if((geo_loop==1).or.(step_counter==0)) then
!Using R0 as value of initial TRUST Raduis    
        e_prev=energy
        e_curr=energy
        accepted=.true.
        r_curr=step_max
        write(io_flepo,2000) r_curr
      else
        e_curr=energy
        r_curr=r_curr
!Compute change energy beetwen previuos and current iterations
!at ab initio (DFT) level of theory.
        de_obj=e_curr-e_prev
        if((abs(de_obj)<=en_small.and.abs(de_mod)<=en_small).and.(step_counter/=0)) then 
        if(abs(de_mod).gt.0.0000001_r8_kind)then
          r=de_obj/de_mod
          write(io_flepo,1800) r
        endif
          write(io_flepo,2200)
          r_curr=r_curr
          e_cur=energy
          write(io_flepo,2100) r_curr
          accepted=.true.
          return
        endif
!Compute parameter r, defining quality of approximating the 
!objective function by model(quadratic) function.
        r=de_obj/de_mod
!Save energy value
          if(less_t_r_curr) then
            r_curr=r_curr
          else
!Analyzing of r-value
            if((r<=r1e).or.(r>=r1i)) then
              r_curr=r_curr/sf
            elseif((r<=rue).and.(r>=rui)) then
              r_curr=r_curr*dsqrt(sf)
            else 
              r_curr=r_curr
            endif
          endif
    if((r_curr<=0.005_r8_kind).and.(step_counter/=0_i4_kind))  r_curr=0.005_r8_kind
      write (io_flepo,1600) de_obj
      write (io_flepo,1700) de_mod
      write (io_flepo,1800) r
      write (io_flepo,1900) r_curr

!    if((r<=zero).or.(r>=1000_r8_kind)) then
    if(abs(r)>=1000_r8_kind) then
      accepted=.false.
      e_curr=e_prev
      e_prev=e_prev
    else    
      accepted=.true.
      e_prev=energy
    endif
    endif
    e_cur=e_curr
    endif  
    DPRINT "step_counter",step_counter
    write(io_flepo,'(27X,A)') "--------------------"
1600 format("Actual change energy is      [DFT] :",F15.8)
1700 format("Predictable change energy is [LQA] :",F15.8)
1800 format("Ratio=[DFT]/[LQA]                  :",F15.8)
1900 format("Current value of TR is             :",F15.8)
2000 format("Initial value of TR is             :",F15.8)
2100 format("Frozen value of TR is              :",F15.8)
2200 format(27X,"--- Is switch off ---")
 end subroutine tr
  
    
 function ang(step)
  !-------------------------------------------------------------
  !Compute angle(in degree) between step and force vectors
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
    use math_module, only : convert1,zero
    use gradient_module,only:grad_intern
    USE_DEBUG! only: isNaN
    implicit none
  !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in) :: step(:)
    real(kind=r8_kind)            :: ang

    ang=zero
    ang=acos(-dot_product(grad_intern,step)/(sqrt(dot_product(step,step))*&
       sqrt(dot_product(grad_intern,grad_intern))))*convert1
#ifdef WITH_ISNAN
    if(isNaN(ang)) ang=zero
#else
    if( ang /= ang ) ang = zero ! FIXME: should holds for NaNs
#endif
 end function ang


 subroutine eliminate_step(step)
  !-------------------------------------------------------------
  !This subroutine write step information in flepo file.
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
  use opt_data_module
  use math_module,only:zero
  implicit none
  !------------ Declaration of formal parameters ---------------
  real(kind=r8_kind),intent(inout)  :: step(:)
  !------------ Declaration of local variables -----------------
  integer(kind=i4_kind)             :: i
  real(kind=r8_kind),parameter      :: sm_comp=0.000001
  !------------ Executable code --------------------------------  
!Eliminate small components of step
  do i=1,n_internal
     if(abs(step(i))<=sm_comp) then
        step(i)=zero
        write(io_flepo,'(A,4X,I4)')"Eliminating component of step:",i
     endif 
  enddo
 end subroutine eliminate_step


 subroutine write_step(meth,step,f_vec,hesse_eigval,mu,tht,step_internal)
  !-------------------------------------------------------------
  !This subroutine write step information in flepo file.
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
  use opt_data_module
  use math_module,only:convert1
  use gradient_module,only:grad_intern
  implicit none
  !------------ Declaration of formal parameters ---------------
  character(len=*),intent(in)       :: meth
  real(kind=r8_kind),intent(in)     :: step(:),f_vec(:),hesse_eigval(:)
  real(kind=r8_kind),intent(in),optional:: step_internal(:)
  real(kind=r8_kind),intent(in)     :: mu,tht
  !------------ Declaration of local variables -----------------
  character(len=4)                  :: std="--- "
  character(len=80)                 :: title,en
  integer(kind=i4_kind)             :: ln,i
  !------------ Executable code --------------------------------  
      title=std//meth//'_STEP'//adjustr(std)
      ln=len_trim(title)
      en=repeat('-',ln)
      write(io_flepo,'(27X,A)') title
      if(meth.eq."GDIIS".or.meth.eq."Line_Search") then
        write(io_flepo,'(27X,A)')"  type         step"
        do i=1,n_internal
           write(io_flepo,'(29x,i3,6x,es12.5)')i,step(i)
        enddo
      else
        write(io_flepo,*)"   mode    grad. intern.     grad. comp.       eigenvalue     type       step"
        do i=1,n_internal
       if(present(step_internal)) then
           write(io_flepo,'(3x,i3,6x,es12.5,5x,es12.5,6x,es12.5,4x,i3,5x,2es12.5)')&
           i,grad_intern(i),f_vec(i),hesse_eigval(i),i,step(i),step_internal(i)
       else
           write(io_flepo,'(3x,i3,6x,es12.5,5x,es12.5,6x,es12.5,4x,i3,5x,es12.5)')&
           i,grad_intern(i),f_vec(i),hesse_eigval(i),i,step(i)
       endif
        enddo
      endif
      write(io_flepo,'(27X,A)') en
      write(io_flepo,'(A,1X,A,F15.8)')"Theta(Degree) between step and forces ",":",tht
 end subroutine write_step


 
                                                                                                                          
  !*************************************************************
  subroutine newton(lambda_val,follow_mode)
    ! Purpose: solves the equation
    !          sum_i(Fi**2/(lambda-eigval(i)) - lambda) = 0
    !          iteratively via a Newton-Procedure.
    !          eigval(i)  = Eigenvalues of the Hessian (n_internal)
    !          lambda_val(j) = Eigenvalues of the augmented Hessian
    !                          (actually n_internal+1, but n_internal in
    !                           practice)
    !          f_vec(i) = Component of the gradients along the 
    !                     eigenmodes of the Hessian.
    !          follow_mode: index of the mode being followed. The
    !                       eigenvalue belonging to this index
    !                       has to be left out of the procedure, thus
    !                       leaving n_internal solutions lambda_val
    !
    ! Reference: J.Baker: An Algorithm for the Location of 
    !                      Transition States, J.Comp.Chem. 7,
    !                      385-395,(1986)
    !        
    ! ---------------------------------------------------------
    use hesse_module, only: hesse_eigval
    use math_module, only: small,zero 
    implicit none
    real(kind=r8_kind),intent(out)   :: lambda_val(:)
    integer(kind=i4_kind)            :: follow_mode
    ! ---------------------------------------------------------
    integer(kind=i4_kind)            :: i,newt_counter
    real(kind=r8_kind)               :: lam_start
    logical                          :: conv
    integer(kind=i4_kind),parameter  :: max_newt = 50_i4_kind
    real(kind=r8_kind),parameter     :: delta = 0.00001_r8_kind
    ! ---------------------------------------------------------
    lambda_val = zero
    lambdas: do i=1,ubound(hesse_eigval,1)
       if (i==follow_mode) cycle lambdas
       lam_start = hesse_eigval(i)-delta
       
       conv = .false.
       newt_counter=0
       newt: do
          if (abs(func(follow_mode,lam_start))<=small) then
             conv=.true.
             exit newt
          else
             conv=.false.
          endif
          newt_counter=newt_counter+1
          if (newt_counter>max_newt) exit newt
          lam_start = lam_start - &
               func(follow_mode,lam_start)/func_prim(follow_mode,lam_start)
       enddo newt
       if (conv) then
          lambda_val(i) = lam_start
       endif
    enddo lambdas
    write(OPT_STDOUT,*)" newton: The calculated shift parameters for the eigenvalues are:"
    write(OPT_STDOUT,*)"         eigval         shift par."
    do i=1,ubound(hesse_eigval,1)
       write(OPT_STDOUT,'(10x,f9.5,9x,f9.5)')hesse_eigval(i),lambda_vec(i)
    enddo
    write(OPT_STDOUT,*)" "
    

  contains

    function func(follow_mode,lam)
      use math_module, only : zero
      implicit none
      real(kind=r8_kind)     :: func
      real(kind=r8_kind)     :: lam
      integer(kind=i4_kind)  :: follow_mode,i
      
      func=zero
      do i=1,ubound(hesse_eigval,1)
         if (i/=follow_mode) then
            func = func + f_vec(i)**2/(lam-hesse_eigval(i))
         endif
      enddo
      func = func - lam
         
    end function func

    function func_prim(follow_mode,lam)
      use math_module, only : zero,one
      implicit none
      real(kind=r8_kind)     :: func_prim
      real(kind=r8_kind)     :: lam
      integer(kind=i4_kind)  :: follow_mode,i
      
      func_prim=zero
      do i=1,ubound(hesse_eigval,1)
         if( i/=follow_mode) then
            func_prim = func_prim - f_vec(i)**2/(lam-hesse_eigval(i))**2
         endif
      enddo
      func_prim = func_prim - one
    end function func_prim
  end subroutine newton

subroutine hebden(mu,lambda,r_curr)
    !-------------------------------------------------------------
    !Purpose: solves the equation
    ! 	  one-sum_i[(f_vec(i)/(hesse_eigval(i)-lambda+mu))]^2/r_curr=zero
    !         iteratively via a Newton-Method
    !	  hesse_eigval(i) = Eigenvalues of current Hessian (n_internal)
    !         f_vec(i)	  = Component of gradients along eigenmodes of the
    !			    current Hessian	
    !         lambda,mu       = shift parameter's. lambda<=zero;mu>=zero
    !
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use type_module
    use opt_data_module
    use math_module,only:one,zero,small
    use hesse_module,only:hesse_eigval
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in)   :: lambda,r_curr
    real(kind=r8_kind),intent(out)  :: mu
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)           :: i
    integer(kind=i4_kind),parameter :: maxiter=1000
    real(kind=r8_kind)              :: ad,mu_prev
    real(kind=r8_kind)              :: delta = 0.005_r8_kind
    !------------ Executable code ------------------------------------
    mu=-hesse_eigval(1)+delta 
    DPRINT "===================="
    DPRINT "Compute mu"
    DPRINT "Iter		mu"
hebd: do i=1,maxiter
           ad=-func(mu,lambda)/dfunc(mu,lambda)
!Using Secant method
           if(dfunc(mu,lambda)<=small) then 
             DPRINT "Using Secant method"
             ad=-func(mu,lambda)*(mu-mu_prev)/(func(mu,lambda)-func(mu_prev,lambda))
           endif
           if((abs(ad))<small) exit hebd
           mu=mu+ad
           mu_prev=mu-ad
           DPRINT i,mu
        end do hebd
    DPRINT "End compute mu"
    contains
    
    function func(mu,lambda)
    use math_module, only : zero,one
    implicit none
    real(kind=r8_kind),intent(in)     :: mu,lambda
    real(kind=r8_kind)                :: func
    integer(kind=i4_kind)             :: i
    func=zero
    do i=1,n_internal
       func=func+(f_vec(i)/(hesse_eigval(i)-lambda+mu))**2
    enddo
    func=dsqrt(func)
    func=1-func/r_curr
    end function func
        
    function dfunc(mu,lambda)
    use math_module, only : zero,one 
    implicit none
    real(kind=r8_kind),intent(in)     :: mu,lambda
    real(kind=r8_kind)                :: num,denum
    real(kind=r8_kind)                :: dfunc
    integer(kind=i4_kind)             :: i
    dfunc=zero;num=zero;denum=zero
    do i=1,n_internal
       denum=denum+(f_vec(i)/(hesse_eigval(i)-lambda+mu))**2
       num=num+(f_vec(i))**2/(hesse_eigval(i)-lambda+mu)**3
    enddo
    denum=r_curr*dsqrt(denum)
    dfunc=num/denum
    end function dfunc        
  end subroutine hebden
  !*************************************************************
end module ts_module

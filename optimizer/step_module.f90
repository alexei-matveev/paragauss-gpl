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
module step_module
  !-------------------------------------------------------------------
  !  Purpose: PRELIMINARY. Contains all routines and data
  !           to calculate the next step. Could also
  !           be called newton_module or update_module
  !
  !  Module called by: main_opt
  !  References: ...
  ! 
  !
  !  Author: FN
  !  Date: 7/97
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: VVP
  ! Date:   11/05
  ! Description: Added:
  !                        1) AHRS,RQN,RS-RFO methods;
  !                        2) Spectral Representation of hesse matrix on each 
  !                           step and analyse of inertia;
  !                        3) Trust  Radius algorithm[TR];
  !                        4) Khait algorithm;
  !                        5) Quartic inexact line search;
  !                        5) GDIIS in various modifications;
  !                        etc.[Please see articles]
  ! References: 
  !   AH-a0:    1. Khait, Yu.G., Panin, A.I., Averyanov, A.S. 
  ! [  Khait       Search for Stationary Points of Arbitrary Index by Augmented Hessian Method 
  ! algortihm]     // Int.J.Quant.Chem. - 1995. - V. 54. - N. 6 - P. 329-336
  !             2. Khait, Yu. G., Panin, Yu. G. Search for stationary points on potential energy surfaces 
  !                // Theochem. - 1997. - N. 398-399. - P.101-109, 103.
  !   AHRS:         3. Anglanda, J.M., Bofill, J.M. On the restricted step method coupled with the augmented Hessian 
  !                        for the search of stationary points of any continuous function // Int.J.Quant.Chem. - 1997. 
  !                - V. 62. - N. 2 - P. 153-165
  !   Hebden:   4. R. Flecther. Practical methods of optimization.Second edition. 1988. JOHN WILEY & SONS
  !                [Wiley-Interscience publication]. P.104-106
  !   TR:           5. Culot, P., Dive, G., Nguyen, V.H., Ghuysen, J.M. A quasi-newton algorithm
  !                        for first-order saddle-point location  // Theor.Chim.Acta. - 1992. - V. 82. - P.189-205, 198.  
  !                and [3].
  !  RS-RFO:    6. Besalu, E., Bofill, J.M. On the automatic restricted-step ratsional-function-optimization
  !                        method // Theo. Chem. Acc. - 1998. - V. 100. - P. 263-274
  !  [Wales     7. Wales, D.J.  Structural and topological consequences of anisotropic interactions in 
  !  shift]        clusters // J. Chem. SOC., Faraday Trans. - 1990. - V. 86. - N. 21 - P. 3505-3517
  !  Greenstadt 8. Schlegel, H.B. Optimization of equilibrium geometries and transition structures //
  !                J. Comp. Chem. - 1982. - V. 3. - N. 2 - P. 214-218 
  !  GDISS      9. Farkas, O., Schlegel, H.B. Methods for optimizing large molecules. Part III. An improved 
  !                algortihm for geometry optimization using direct inversion in the iterative subspace (GDIIS) // 
  !                Phys.Chem.Chem.Phys. - 2001. - V. 14. - P. 11-15.
  !            10. Eckert, F., Pulay P., Werner, H.J. Ab initio optimization for large molecules // J.Comp.Chem. -
  !                V. 18. - N 12. - P. 1473-1483
  !   All      11. Vysotskiy, V. Modern methods of optimization of geometry and search of transition states. // 
  ! together       2005. - TUM. [/home/vysotskiy/DOCS/Report_MEMO/MEMO_OPTI.pdf]
  !Note: GDIIS used in the end of optimization process for convergence acceleration. Accepted by Prof. Dr. J.M.Bofill.  
  !
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use gradient_module, only: grad_intern,energy!,energy_ph
  use math_module, only: zero,small,abs_value
  use coordinates_module
  use coortype_module
  use hesse_module, only:hesse,hesse_eigvec,hesse_eigval,step_counter,mep_point
  use opt_data_module !!!!!!!!!!!!!
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of constants and variables ---------------
  real(kind=r8_kind),       public     :: step_mean_square,step_max_comp

  ! these variables are made public only to allow saving/restoring persistent state:
  real(kind=r8_kind),       public     :: r_curr ! Trust radius
  real(kind=r8_kind),       public     :: e_prev ! Energy from the previous geo-iteration

  real(kind=r8_kind)                   :: lambda,mu
  real(kind=r8_kind)                   :: de_mod      ! (persistent) Energy change dE as predicted by the current model
  real(kind=r8_kind)                   :: e_save      ! (persistent)
  real(kind=r8_kind)                   :: step_length ! (persistent)
  logical                              :: less_t_r_curr != .true.  ! (persistent), initial value ok?
  real(kind=r8_kind),allocatable       :: grad_save(:),q_save(:)  ! (persistent), used to step back from bad trial geometries

  integer(kind=i4_kind)                :: n_negative
  real(kind=r8_kind),allocatable       :: step(:)
  real(kind=r8_kind),allocatable       :: ahesse(:,:),scl(:,:)
  real(kind=r8_kind),allocatable       :: ahesse_eigvec(:,:),ahesse_eigval(:)
  real(kind=r8_kind),allocatable       :: f_vec(:)
  logical,public:: noline,new_mep_point
  integer(kind=i4_kind),public:: number_nolines=0
  !------------ public functions and subroutines ---------------------
  public :: newton_step
  public :: step_module_persistent_state ! saves/restores gdiis_work(:num_stpt) ...
  
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------
  
  ! This type is very helpful for GDIIS  
  type gdiis_point
   integer(kind=i4_kind)             :: num_pt = -1
   real(kind=r8_kind),pointer        :: q_gdiis(:) => NULL()
   real(kind=r8_kind),pointer        :: g_gdiis(:) => NULL()
   real(kind=r8_kind)                :: ener = 0.0
  end type gdiis_point

  integer(kind=i4_kind),parameter    :: max_point = 4
  integer(kind=i4_kind),parameter    :: ns_point = 4
  integer(kind=i4_kind)              :: num_stpt = 0          ! (persistent)
  type(gdiis_point)                  :: gdiis_work(max_point) ! (persistent) Stores info about previous "num_stpt" points
  logical                            :: gdiis_storage_initialized = .false.

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine step_setup(mth)
    !-------------------------------------------------------------
    ! Purpose: depending on the chosen method (Quasi-Newton,
    !          rationale function/ augmented hessian etc.) 
    !          all required variablea are allocated and initialized
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use opt_data_module, only: n_internal
    !------------ Declaration of formal parameters ---------------
    character(len=*),    intent(in)  :: mth
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)            :: alloc_stat,ah_dim
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

    if(tsscan_sphere) then
    allocate(f_vec(n_internal-1),STAT=alloc_stat)
    else
    allocate(f_vec(n_internal),STAT=alloc_stat)
    endif
    ASSERT(alloc_stat==0)

    if(mth.ne."qn") then
      ah_dim=n_internal+1
      if(tsscan_sphere) ah_dim=n_internal
      allocate(ahesse_eigval(ah_dim),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      allocate(ahesse_eigvec(ah_dim,ah_dim),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      allocate(ahesse(ah_dim,ah_dim),STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      if(mth.eq."rfo") then
        allocate(scl(ah_dim,ah_dim),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
      endif
    endif

 end subroutine step_setup
  
 subroutine constr_ah(mu,lambda,r_curr,hesse_eigval,f_vec,ahesse)
    !-------------------------------------------------------------
    ! Purpose: Construct augmented hessian matrix
    !------------ Modules used -----------------------------------
    use math_module,only:zero
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),    intent(in)    :: hesse_eigval(:),f_vec(:)
    real(kind=r8_kind),    intent(in)    :: mu,lambda,r_curr
    real(kind=r8_kind),    intent(inout) :: ahesse(:,:)
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: i
    !------------ Executable code ------------------------------------
    ahesse=zero
    ahesse(1,1)=-mu*r_curr*r_curr
    ahesse(1,2:)=f_vec(:)
    ahesse(2:,1)=f_vec(:)
    do i=1,size(hesse_eigval)
       ahesse(i+1,i+1)=hesse_eigval(i)+mu-lambda
    end do
 end subroutine constr_ah
  
 subroutine step_shutdown(mth)
    !-------------------------------------------------------------
    ! Purpose: all used array are deallocating
    !-------------------------------------------------------------
    !------------ Declaration of formal parameters ---------------
    character(len=*),    intent(in)      :: mth
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)   :: alloc_stat
    !------------ Executable code ------------------------------------
    deallocate(step,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    deallocate(f_vec,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    if(mth.ne."qn") then
      deallocate(ahesse_eigval,STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      deallocate(ahesse_eigvec,STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      deallocate(ahesse,STAT=alloc_stat)
      ASSERT(alloc_stat==0)
      if(mth.eq."rfo") then
        deallocate(scl,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
      endif
    endif    
 end subroutine step_shutdown
 
 function inertia(mu,lambda) result(n_negative)
    !-------------------------------------------------------------
    !  Purpose: compute inertia of hessian
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use type_module
    use hesse_module,only:hesse_eigval
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in)        :: mu,lambda
    integer(kind=i4_kind)                :: n_negative
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: i
    !------------ Executable code ------------------------------------
    n_negative=0_i4_kind
    do i=1,size(hesse_eigval)
       if (hesse_eigval(i)+mu-lambda<zero) then
       n_negative=n_negative+1
       endif
    enddo
 end function inertia

! MARK ts_serach2 useful routine?    
 subroutine newton_step(geo_loop)
    !-------------------------------------------------------------
    !Compute step in various approximation to PES
    !-------------------------------------------------------------
    !  Purpose: does the quasi Newton-Step:
    !           s_vec = H**(-1)*grad_vec
    !  VVP 11/05: Solve more general problem:
    !        [1] Perform diagonalization of Hessian;
    !        [2] Then perform analysis of its spectra;
    !        [3] Choose method depend on inertia and TR model;
    !        [4] Calculate step. 
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use type_module
    use math_module, only: zero,small,abs_value,print_matrix,one,half
    use opt_data_module
    use coordinates_module
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in) :: geo_loop
    !------------ Declaration of local variables ---------------------
    character(len=4)                     :: mth

    !------------ Executable code ------------------------------------

    ! Define method:
    quart=.false.
    if(method_qn)  mth="qn"
    if(method_ah)  mth="ah"
    if(method_rfo) mth="rfo"
    if(method_qn)  quart=.true.

    ! allocates module global vars step(:) and more:
    call step_setup(mth)

    ! Main block
    ! fill module global vars step(:) with data:
    call compute_step(mth,geo_loop, step) !(1)
    if(.not.new_mep_point) step_counter=step_counter+1

    ! deallocates module global vars step(:) and more:
    call step_shutdown(mth)
  end subroutine newton_step
  
 subroutine compute_step(meth,geo_loop, step)
    !-------------------------------------------------------------
    !  Purpose: Solve RQN or AHRS or RS-RFO equations,
    !           compute step in this approach.
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use type_module
    use opt_data_module,only:n_internal,dynamic,step_max,method_gdiis, make_qst_step=>qst_step
    use gradient_module, only: grad_intern,energy,grad_max_sphere,grad_mean_sphere
    use math_module
    use matrix_eigenval
    use coortype_module
    use hesse_module, only: hesse_eigvec,hesse_eigval
    use hesse_module,only:hess_and_grad_project,write_old_hess
    USE_DEBUG
     use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*),intent(in)          :: meth
    integer(kind=i4_kind),intent(in)     :: geo_loop
    real(r8_kind), intent(inout)         :: step(:)

    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)                :: n_negative,i,j,alloc_stat,io_pointonmep
    real(kind=r8_kind)                   :: mu,lambda_ahrs,lambda_ah,a0_ahrs
    real(kind=r8_kind)                   :: lambda_rfo,a0_rfo,denum,alpha,e_new,beta,a0_ah
!   real(kind=r8_kind),save              :: de_mod ! local declaration eclipses the global declaration
!   real(kind=r8_kind),save              :: e_save ! made global
!   real(kind=r8_kind),save              :: step_length ! made global
    real(kind=r8_kind),allocatable       :: q_new(:),g_new(:)
    real(kind=r8_kind)                   :: e_cur
    real(kind=r8_kind)                   :: tht,ad
    real(kind=r8_kind),parameter         :: scale1=0.6_r8_kind,small_eig=0.0001_r8_kind
    real(kind=r8_kind),parameter         :: eps1=0.00005,eps2=0.0000005
    logical                              :: desired,accepted,do_gdiis
!   logical,save                         :: less_t_r_curr ! made global
    !------------ Executable code ------------------------------------

   if(step_reset) r_curr=step_max

   if(tsscan_sphere) then
    new_mep_point=new_mep_point.and.grad_max_sphere.lt.max_comp_grad*2.and.grad_mean_sphere.lt.rms_grad*2
    new_mep_point=grad_max_sphere.lt.max_comp_grad.and.grad_mean_sphere.lt.rms_grad
!    new_mep_point=new_mep_point.and.step_counter.ne.0
   else
    new_mep_point=.false.
   endif

   if(tsscan_sphere) then

    if(mep_point.eq.0) then
     write(io_flepo,*) ' point on mep is set to reactant'
     s_pointonmep(:)%value=s_reactant(:)%value
     io_pointonmep=openget_iounit(status='unknown',form='formatted', &
          file=trim(opt_data_dir)//'/point_on_mep.dat')
     write(io_pointonmep,*) s_pointonmep(:)%value
     mep_point=mep_point+1
     call returnclose_iounit(io_pointonmep) 
    endif

    if(new_mep_point) step=step_to_new_mep_point(q_old,de_mod) 

   endif

    dynamic1: if(dynamic.and..not.new_mep_point) then
      call tr(geo_loop,de_mod,less_t_r_curr,r_curr,e_cur,accepted)
      if(.not.accepted) then
    !Bifurcation: at first pass we try to use inexact/partial line search by quartic interpolation
          noline=.true.  !(1)
          allocate(q_new(n_internal),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(g_new(n_internal),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          if (quart) call quartic(q_save,q_old,grad_save,grad_intern,e_save,energy,q_new,e_new,noline)
    !If line search is unsuccescfull used prevuoisly geometry and gradient.
          if(noline) then
            write(io_flepo,'(A,7x,I10)')"Prevuois step was unacceptable. Using geometry from step:",geo_loop-1
            q_old=q_save
            grad_intern=grad_save
            if(step_reset) then
            r_curr=step_max
            elseif(less_t_r_curr.and.(r_curr>=step_length*scale1)) then
               r_curr=step_length*scale1
               write(io_flepo,'(A)')"Rescale Trust Radius[multiply by factor 0.6]"
            endif
          else
          number_nolines=0
     !If line search is succescfull used results estimated geometry. My practical experience 
     !showed, that estimaed gradient by line search method in unappropriate. .           
            if(.not. cart_coordinates) then
               e_prev=energy
            else
               e_prev=e_save
            end if
            q=q_new
            if(.not. cart_coordinates) then
               step=q_new-q_old
            else
               step=q_new-q_save
            end if
            if(.not. cart_coordinates) call correct_step(step)
            if(.not. cart_coordinates) then
               tht=ang(step,grad_intern)
            else
               tht=ang(step,grad_save)
            end if
            f_vec=zero
            call write_step("Line_Search",step,f_vec,hesse_eigval,tht)
!!$            if(cart_coordinates .and. tht > 89.0_r8_kind) reinit_hess=.true.
            if(.not. cart_coordinates) then
               de_mod=e_new-energy
            else
               de_mod=e_new-e_save
            end if
            step_mean_square = dsqrt((sum(step**2))/n_internal)
            step_max_comp = maxval(abs(step))
            step_length=dsqrt(dot_product(step,step))
            write(io_flepo,1100) energy
            write(io_flepo,1200) de_mod
            write(io_flepo,1400) step_length
            deallocate(q_new,STAT=alloc_stat)
            ASSERT(alloc_stat==0)
            deallocate(g_new,STAT=alloc_stat)
            ASSERT(alloc_stat==0)
            if(cart_coordinates) then
               call write_old_hess()
            end if
            return
          endif
        deallocate(q_new,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        deallocate(g_new,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
      endif
    else
    e_cur=energy
    endif dynamic1

  if(.not.new_mep_point) then

    number_nolines=0

     !Method can be:qn or ah, rfo. 
     !Default:meth.eq.'qn'
     !Quasi-Newton method

    if(.not.cart_coordinates) then
       f_vec=eigvec_grads(grad_intern,hesse_eigvec,tsscan_sphere)
    else
       call hess_and_grad_project(f_vec,grad_intern)
    end if

    select case(meth)

    case("qn")

    n_negative= inertia(zero,zero)
    call desrd("QN",n_negative,desired) 
    if(.not.cart_coordinates) then
       step=qn_step(sphere_dependent_var,f_vec)
    else
       !Project the gradient onto the Hessian eigenvectors
       f_vec=matmul(transpose(hesse_eigvec),f_vec)
!!$print*,'GRADS in special form'
!!$print*,f_vec
       do i=1,n_internal
          step(i)=-f_vec(i)/hesse_eigval(i)
       end do
!!$print*,'STEP in special form'
!!$print*,step
       !Transform back to optimization space
       step=matmul(hesse_eigvec,step)
!!$print*,'STEP normal'
!!$print*,step
       step = matmul(tmat,step)
!       call eliminate_step(step)
    end if

    step_length=dsqrt(dot_product(step,step))
    write(io_flepo,'(27X,A)')"--- Quasi-Newton ---"

    if(desired.and.(step_length<=r_curr).and.dynamic.and.accepted) then

      !GDIIS block: start
      less_t_r_curr =.true.

      if(method_gdiis) then
        allocate(q_new(n_internal),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        allocate(g_new(n_internal),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        call gdiis(geo_loop,q_old,zero,meth,q_new,g_new,do_gdiis)

        if(do_gdiis) then
          write(io_flepo,'(A)')"Using geometry and gradient getting from GDIIS"
          grad_save=grad_intern
          q_save=q_old
          e_save=e_cur
          step=q_new-q_old
          de_mod=dot_product(grad_intern,step)+half*dot_product(step,matmul(hesse,step)) !(1)
          q=q_new
          f_vec=zero;
          tht=ang(step,grad_intern)
          call eliminate_step(step) !!!!!
          step_length=abs_value(step)
          call write_step("GDIIS",step,f_vec,hesse_eigval,tht)
          step_mean_square = dsqrt((sum(step**2))/n_internal)
          step_max_comp = maxval(abs(step))
          write(io_flepo,1100) e_cur
          write(io_flepo,1400) step_length
          write(io_flepo,1200) de_mod
          deallocate(q_new,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(g_new,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          return
!GDIIS block: end       
        else
          de_mod=zero
          do i=1,size(hesse_eigval)
             de_mod=de_mod-half*f_vec(i)**2/hesse_eigval(i)
          enddo
          tht=ang(step,grad_intern)
          write (io_flepo,'(A)')"Making QN step (GDIIS)"
          call write_step("QN",step,f_vec,hesse_eigval,tht)
        endif

      else
      write (io_flepo,*) '.not.method_gdiis'
        if(tsscan_sphere) then
         write(io_flepo,*) 'sphere_noreturn  executed'
         do  i=1,9
         if(.not.sphere_noreturn(sphere_dependent_var,step)) exit
         step_max=step_max*0.7
          write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
          "RQN method: step length exceeds maximum value: ",step_length,step_max
          call restrict(step_max, step)
         enddo
         step_length=sqrt(sum(step(1:n_internal-1)**2))
        endif
       de_mod=zero
       do i=1,size(hesse_eigval)
          de_mod=de_mod-half*f_vec(i)**2/hesse_eigval(i)
       enddo
       tht=ang(step,grad_intern)
       write (io_flepo,'(A)')"Making QN step (default)"
       call write_step("QN",step,f_vec,hesse_eigval,tht)
      endif 

    else ! i.e. not all is fine

       if(dynamic) then
       write(io_flepo,*) 'dynamic'
         less_t_r_curr=.false.
         if(step_length > r_curr) then
            write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
                 "QN method: step length exceeds maximum value: ",step_length,r_curr
         end if
         write(io_flepo,'(27X,A)')"--- Solve RQN equations ---" !(1)
         call hebden(mu,zero,r_curr)
         if(.not.desired.and.mu.lt.-hesse_eigval(1) ) mu=-hesse_eigval(1)*1.03
         write (io_flepo,'(A,F15.8)')"mu     =",mu
         write(io_flepo,'(A)')"Making RQN-step"

         if(.not.cart_coordinates) then
            step=rqn_step(sphere_dependent_var,mu)
         else
            do i=1,n_internal
               step(i)=-f_vec(i)/(hesse_eigval(i)+mu)
            end do
            step=matmul(hesse_eigvec,step)
            step = matmul(tmat,step)
         end if

        if(tsscan_sphere) then
         write(io_flepo,*) 'sphere_noreturn  executed'
         do  i=1,19
         if(.not.sphere_noreturn(sphere_dependent_var,step)) exit
         step_max=step_max*0.7
          write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
          "RQN method: step length exceeds maximum value: ",step_length,step_max
          call restrict(step_max, step)
         enddo
        endif

         call eliminate_step(step)
         step_length=dsqrt(dot_product(step,step))
         de_mod=zero
         do i=1,size(hesse_eigval)
            de_mod=de_mod-(half*hesse_eigval(i)+mu)*f_vec(i)**2/((hesse_eigval(i)+mu)**2)
         enddo
         tht=ang(step,grad_intern)
         call write_step("RQN",step,f_vec,hesse_eigval,tht)

       else ! i.e. not dynamic
         if(.not.desired) then
           write(io_flepo,'(A)')"Use Greenstadt method"
           step = Greenstadt_step(sphere_dependent_var)
           step_length=dsqrt(dot_product(step,step))
       endif

         if(step_length>=step_max) then
           !First stage: maybe current hessian have some small eigenvalues ?
           step=Wales_step(small_eig,sphere_dependent_var)
           write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
           "QN method: step length exceeds maximum value: ",step_length,step_max
           call restrict(step_max, step)
           step_length=dsqrt(dot_product(step,step))
         endif
         tht=ang(step,grad_intern)
         call write_step("QN",step,f_vec,hesse_eigval,tht)
       endif
    endif

!Augmented Hessian method:
    case("ah")
    write(io_flepo,'(A)')"Using Augmented hessian technique"
!Construct Augmented Hessian
    call constr_ah(zero,zero,zero,hesse_eigval,f_vec,ahesse)
!Solve secular equation for AH
    call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
      DCALL show("Augmented Hessian",ahesse(:,:))
      DCALL show("Eigenvectors of Augmented Hessian",ahesse_eigvec(:,:))
      DCALL show("Eigenvalues of Augmented Hessian",ahesse_eigval(:))
!Compute a0
    write (io_flepo,'(A)') "Augmented Hessian"
    step=ah_step(a0_ah,sphere_dependent_var)
    call eliminate_step(step)

    step_length=dsqrt(dot_product(step,step))
    lambda_ah=ahesse_eigval(1)
    write (io_flepo,'(A,F15.8)') "a0_ah=",abs(a0_ah)

    if(abs(a0_ah)<=0.75_r8_kind) then
    write(io_flepo,*) "a0<0.75. So current point is far from solution"
    write(io_flepo,*)"Using method suggested by Khait"
    write(io_flepo,*)"Follow the lowest"
        write(io_flepo,*)"           eigenmode with eigenvalue"
        write(io_flepo,'("         ",f10.6)')hesse_eigval(1)

        de_mod=-half*abs(hesse_eigval(1))*f_vec(1)**2/((hesse_eigval(1))**2)
        step=Khait_step(sphere_dependent_var)

        step_length=dsqrt(dot_product(step,step))

        write(io_flepo,*)"Perform scale of calculated step without change the direction"
        call restrict(step_max, step)
        step_length=dsqrt(dot_product(step,step))
        tht=ang(step,grad_intern)
        call write_step("AH-Khait",step,f_vec,hesse_eigval,tht)
        step = step/step_length*step_max
        de_mod=de_mod/step_length*step_max

    else
!Step[1]
      write(io_flepo,'(A)')"a0>=0.75.Current point is sufficiently close to the solution"
      write(io_flepo,'(A)')"Using Augmented Hessian technique"
      write(io_flepo,'(A,F15.8)')"lambda:",lambda_ah
      n_negative= inertia(zero,lambda_ah)
      call desrd("AH",n_negative,desired)
      if(desired.and.(step_length<=r_curr).and.dynamic) then
        less_t_r_curr=.true.
        write (io_flepo,'(A)')"Making AH step 3"
        write (io_flepo,'(A,F15.8)')"lambda_ah=",lambda_ah
        de_mod=lambda_ah/2*a0_ah*a0_ah
        tht=ang(step,grad_intern)
        call write_step("AH",step,f_vec,hesse_eigval,tht)
      else
!Step[2]
        if(dynamic) then
          less_t_r_curr=.false.
          ahesse_eigvec=zero
          ahesse_eigval=zero
          write(io_flepo,*)"AH method: step length exceeds maximum value "
          write(io_flepo,'("         calculated length ",f9.5," maximum value ",f9.5)')step_length,&
          r_curr
          write (io_flepo,'(A)')"Making AHRS step"
          write(io_flepo,'(A)')"Compute mu. Solve  Augmented Hessian with Restricted Step equation [AHRS]"
          call hebden(mu,zero,r_curr)
          write(io_flepo,*)"mu=",mu
          call constr_ah(mu,lambda_ah,r_curr,hesse_eigval,f_vec,ahesse)
          call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
          DCALL show("New Augmented Hessian[AHRS]",ahesse(:,:))
          DCALL show("New Eigenvectors of Augmented Hessian[AHRS]",ahesse_eigvec(:,:))
          DCALL show("New Eigenvalues of Augmented Hessian[AHRS]",ahesse_eigval(:))

          step=ah_step(a0_ahrs,sphere_dependent_var)
          write(io_flepo,*)"a0_ahrs=",abs(a0_ahrs)
          lambda_ahrs=ahesse_eigval(1)

!Eliminate small components of step
          call eliminate_step(step)
          step=matmul(tmat,step)
          de_mod=lambda_ahrs/2*a0_ahrs*a0_ahrs
          step_length=dsqrt(dot_product(step,step))
          write (io_flepo,*) "lambda_ahrs=",lambda_ahrs
          tht=ang(step,grad_intern)
          call write_step("AHRS",step,f_vec,hesse_eigval,tht)
        else
            if(step_length>=step_max) then
              write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
              "AH method: step length exceeds maximum value: ",step_length,step_max
               call restrict(step_max, step)
               step_length=dsqrt(dot_product(step,step))
            endif
            tht=ang(step,grad_intern)
            call write_step("AH",step,f_vec,hesse_eigval,tht)
        endif
    endif
    endif

    case ("rfo")
    write (io_flepo,'(A)')"use rfo technique"
    call constr_ah(zero,zero,zero,hesse_eigval,f_vec,ahesse)
!Solve secular equation for AH
    call eigs(ahesse,ahesse_eigval,ahesse_eigvec)
    DCALL show("Augmented Hessian",ahesse(:,:))
    DCALL show("Eigenvectors of Augmented Hessian",ahesse_eigvec(:,:))
    DCALL show("Eigenvalues of Augmented Hessian",ahesse_eigval(:))
    lambda_ah=ahesse_eigval(1)
    write (io_flepo,'(A)')"calculate rfo ah_step"
    step=ah_step(a0_ah,sphere_dependent_var)
    call eliminate_step(step)
    step_length=dsqrt(dot_product(step,step))

    write(io_flepo,*)""
    write(io_flepo,'(27X,A)')"--- Augmented Hessian ---"
    write (io_flepo,'(A,F15.8)')"a0_ah      =",abs(a0_ah)
    write (io_flepo,'(A,F15.8)')"lambda_ah  =",lambda_ah

    n_negative= inertia(zero,zero)
    call desrd("AH",n_negative,desired) 

    des: if(desired.and.(step_length<=r_curr).and.dynamic.and.accepted) then
      less_t_r_curr=.true.
      if(method_gdiis) then
        allocate(q_new(n_internal),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        allocate(g_new(n_internal),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        call gdiis(geo_loop,q_old,a0_ah,meth,q_new,g_new,do_gdiis)
        if(do_gdiis) then
          write(io_flepo,'(A)')"Using geometry and gradient getting from GDIIS"
          grad_save=grad_intern
          q_save=q_old
          e_save=e_cur
          step=q_new-q_old
          de_mod=dot_product(grad_intern,step)+half*dot_product(step,matmul(hesse,step))
          q=q_new
          f_vec=zero;
          tht=ang(step,grad_intern)
          call eliminate_step(step)
          step_length=abs_value(step)
          call write_step("GDIIS",step,f_vec,hesse_eigval,tht)
          step_mean_square = dsqrt((sum(step**2))/n_internal)
          step_max_comp = maxval(abs(step))
          write(io_flepo,1100) e_cur
          write(io_flepo,1400) step_length
          write(io_flepo,1200) de_mod
          return
        else
          write (io_flepo,'(A)')"Making AH step 1"
          de_mod=lambda_ah/two*a0_ah*a0_ah
          tht=ang(step,grad_intern)
          call write_step("AH",step,f_vec,hesse_eigval,tht)
        endif
          deallocate(q_new,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(g_new,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      else
        write (io_flepo,'(A)')"Making AH step 2"
        de_mod=lambda_ah/two*a0_ah*a0_ah
        tht=ang(step,grad_intern)
        call write_step("AH",step,f_vec,hesse_eigval,tht)
      endif

    else des
!Step[2]
      DPRINT 'desired=',desired,'dynamic=', dynamic,'accepted=',accepted
      DPRINT 'step_length<=r_curr is',step_length<=r_curr,'as',step_length,'cmp',r_curr
      if(dynamic.and.step_length>r_curr) then
        less_t_r_curr=.false.
        ahesse_eigvec=zero
        ahesse_eigval=zero
        lambda_rfo=lambda_ah
        write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
        "AH method: step length exceeds maximum value: ",step_length,r_curr
        write(io_flepo,'(27X,A)')"--- Solve RS-RFO equations ---"

!Start of microiterative process
        alpha=one;
        j=1_i4_kind
        DPRINT "===================="
        DPRINT "Compute     Alpha"
        DPRINT "Iter       Alpha"

        call constr_ah(zero,zero,zero,hesse_eigval,f_vec,ahesse)
        step=rsrfo_step(alpha,a0_rfo,lambda_rfo,sphere_dependent_var)
        call eliminate_step(step)

         de_mod=lambda_rfo*(one+alpha*step_length**2)/two
         step_length=dsqrt(dot_product(step,step))
         write (io_flepo,'(A,F15.8)')"a0_rfo     =",abs(a0_rfo)
         write (io_flepo,'(A,F15.8)')"lambda_rfo =",lambda_rfo
         write (io_flepo,'(A,F15.8)')"alpha      =",alpha
         write (io_flepo,'(A)')"Making RS-RFO step"
         tht=ang(step,grad_intern)
         call write_step("RS-RFO",step,f_vec,hesse_eigval,tht)
      else 

        if(step_length>=step_max) then
          write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
          "AH method: step length exceeds maximum value: ",step_length,step_max
          call restrict(step_max, step)
        endif

        if(tsscan_sphere) then
         do  i=1,5
         if(.not.sphere_noreturn(sphere_dependent_var,step)) exit
         step_max=step_max*0.7
          write(io_flepo,'(A,"calculated length=",f15.8,", maximum value=",f15.8)')&
          "AH method: step length exceeds maximum value: ",step_length,step_max
          call restrict(step_max, step)
         enddo
        endif

        tht=ang(step,grad_intern)
        call write_step("AH",step,f_vec,hesse_eigval,tht)
      endif
    endif des

    case("cg")
     step=cg_step()

    case default
    write(io_flepo,'(A)') "Something is wrong"
    call error_handler ("Compute_step: Unknown method ")
   end select
  endif ! .not.new_mep_point

!Write in flepo file results:
    if(dynamic) then
      if(new_mep_point) then
         write(io_flepo,*) "Total energy                           :", e_cur,.not.new_mep_point
      else
         write(io_flepo,*) "Total energy                           :", e_cur,.not.new_mep_point
      end if
      if(meth.eq."qn")  write(io_flepo,1200) de_mod
      if(meth.ne."qn")  write(io_flepo,1300) de_mod
      grad_save=grad_intern
      q_save=q_old
      e_save=e_cur
    else 
      write(io_flepo,*) "Total energy                           :", energy,.not.new_mep_point
    endif

    write(io_flepo,1400) step_length
    step_mean_square = dsqrt((sum(step**2))/n_internal)
    step_max_comp = maxval(abs(step))

!___________________________
    q=q_old+step ! update q
    if(tsscan_sphere.and..not.new_mep_point) then
     if(exist_product) then
      q(sphere_dependent_var)=rp_dependent_var(sphere_dependent_var)

!      s(sphere_dependent_var)%value=q(sphere_dependent_var)
!      call dkdt(sphere_dependent_var)
!      tsscan_rp_var=tsscan_rp_var+0.001_r8_kind
!      q(sphere_dependent_var)=rp_dependent_var(sphere_dependent_var)
!      s(sphere_dependent_var)%value=q(sphere_dependent_var)
!      call dkdt(sphere_dependent_var)

!      call dkdi(sphere_dependent_var)
!      q(1)=q(1)+0.001_r8_kind
!      q(sphere_dependent_var)=rp_dependent_var(sphere_dependent_var)
!      call dkdi(sphere_dependent_var)
      
     else
      q(sphere_dependent_var)=dependent_var(sphere_dependent_var)
     endif
    endif

1100 format("Total energy                           :",F20.8)
1200 format("Predictable change of energy in LQA  is:",F20.8)
1300 format("Predictable change of energy in Pade is:",F20.8)
1400 format("Step length                            :",F20.8)

  contains

  subroutine dkdt(kk)
  integer(kind=i4_kind):: kk
  real(kind=r8_kind):: P2,R2,ssrk,sspk,ssk,rp_var,dRdT
   
   P2=sum((s(:)%value-s_product(:)%value)**2 )
   R2=sum((s(:)%value-s_reactant(:)%value)**2)
   distance_to_reactant=sqrt(R2)
   distance_to_product=sqrt(P2)
   rp_var=distance_to_product/(distance_to_product+distance_to_reactant)
    ssrk=s(kk)%value-s_reactant(kk)%value
    sspk=s(kk)%value-s_product(kk)%value
    ssk=rp_var**2*ssrk-(1-rp_var)**2*sspk
    dRdT=-(P2*(1-rp_var)+R2*rp_var)/ssk
    write(io_flepo,*) s(kk)%value, dRdT, &
     ( P2-R2 -4*dRdT*(rp_var*ssrk +(1-rp_var)*sspk ) -(2*rp_var-1)*dRdT**2 )/ssk,'dRdT d2R/dt2'

  
  end subroutine dkdt
  subroutine dkdi(k)
  integer(kind=i4_kind):: k
  real(kind=r8_kind) :: t,d1,d2
  t=tsscan_rp_var
  d1=-(t**2*(q(1)-s_reactant(1)%value)-(1-t)**2*(q(1)-s_product(1)%value))&
     /(t**2*(q(k)-s_reactant(k)%value)-(1-t)**2*(q(k)-s_product(k)%value))
  d2=-(2*t-1)*(t**2*(q(1)-s_reactant(1)%value)-(1-t)**2*(q(1)-s_product(1)%value))**2 &
    /(t**2*(q(k)-s_reactant(k)%value)-(1-t)**2*(q(k)-s_product(k)%value))**3 &
    -(2*t-1)/(t**2*(q(k)-s_reactant(k)%value)-(1-t)**2*(q(k)-s_product(k)%value))
  write(io_flepo,*) q(1),d1,d2 ,'q(1),d1,d2'
  
  end subroutine dkdi

 function qst_step(s,s_r,s_p) result(step)
  type(int_coor), intent(in):: s(:),s_r(:), s_p(:)
  real(kind=r8_kind):: step(n_internal),gamma(n_internal),d_r,d_p, p, f
  d_r=sqrt(sum( (s(:)%value-s_r(:)%value)**2) )
  d_p=sqrt(sum( (s(:)%value-s_p(:)%value)**2 ) )
  p=d_r/(d_r+d_p)
  write(io_flepo,*) 'current value of reaction coordinate',p
  gamma=( s(:)%value-(1-p)*s_r(:)%value-p*s_p(:)%value )/p/(1-p)
  f=p-step_max
  step=(1-f)*s_r(:)%value+f*s_p(:)%value+f*(1-f)*gamma-s(:)%value
 end function qst_step

  function step_to_new_mep_point(q_old,de_mod) result(step)
    real(kind=r8_kind), intent(out):: q_old(:),de_mod
    real(kind=r8_kind):: step(n_internal),new_distance_to_reactant,new_distance_to_product,new_rp_var
    logical:: mep_point_is_reactant,over_step
     io_pointonmep=openget_iounit(status='old',form='formatted', &
          file=trim(opt_data_dir)//'/point_on_mep.dat')
     read(io_pointonmep,*) s_pointonmep(:)%value
     write(io_flepo,*) s_pointonmep(1)%value,' s_pointonmep read'
     call returnclose_iounit(io_pointonmep) 
     step_counter=0

      step(:)=s(:)%value-s_pointonmep(:)%value
      q_old=s(:)%value
     write(io_flepo,*) '               NEW MEP POINT reached with step of',sqrt(sum(step**2))

     mep_point_is_reactant=sum((s_reactant(:)%value-s_pointonmep(:)%value)**2).lt.0.00001_r8_kind

     if(mep_point_is_reactant) then
     write(io_flepo,*) 'reverse step direction'
       step=-step
       step=grad_intern
     endif

    if(make_qst_step.and..not.mep_point_is_reactant) then
     write(io_flepo,*) 'making qst_step'
     step=qst_step(s,s_reactant,s_pointonmep)
    else
     step_length=sqrt(sum(step(:)**2))
      write(io_flepo,*) 'unrestricted step',sqrt(sum(step(:)**2))
      if(step_max.lt.distance_to_reactant) call scale_step(step,step_max)

      write(io_flepo,*) '  restricted step step_max',sqrt(sum(step(:)**2)),step_max
      call eliminate_step(step)
      step_length=sqrt(sum(step(:)**2))
     endif

     new_distance_to_reactant=new_distance(step)


     write(io_flepo,*) 'new_distance_to_reactant old_distance_to_reactant', &
                                new_distance_to_reactant,distance_to_reactant
     if(exist_product) then
      new_distance_to_product=distance_product(step)
      new_rp_var=new_distance_to_product /(new_distance_to_reactant+new_distance_to_product)
      write(io_flepo,*) new_rp_var,tsscan_rp_var,new_rp_var-tsscan_rp_var,distance_to_ts, &
                      ' new_rp_var(1) old diff dist_to_TS'
      over_step=abs(new_rp_var-tsscan_rp_var).gt.abs(distance_to_ts)
      if( step_to_ts.and.over_step) then
       step_max=0.9_r8_kind*step_length*abs(distance_to_ts/(new_rp_var-tsscan_rp_var))
       call scale_step(step,step_max)
       call eliminate_step(step)
       step_length=sqrt(sum(step(:)**2))
       write(io_flepo,*) 'step scaled to 0.9 of distance to TS ', step_length
       new_distance_to_reactant=new_distance(step)
       new_distance_to_product=distance_product(step)
       new_rp_var=new_distance_to_product /(new_distance_to_reactant+new_distance_to_product)
      write(io_flepo,*) new_rp_var,tsscan_rp_var,new_rp_var-tsscan_rp_var,distance_to_ts, &
                      ' new_rp_var(2) old diff dist_to_TS'
       
      endif
      if( new_rp_var.lt.tsscan_rp_var) then
       write(io_flepo,*) 'step inverted'
       step(:)=-step(:)
       new_distance_to_reactant=new_distance(step)
       new_distance_to_product=distance_product(step)
       new_rp_var=new_distance_to_product /(new_distance_to_reactant+new_distance_to_product)
       write(io_flepo,*) new_rp_var,tsscan_rp_var,new_rp_var-tsscan_rp_var,distance_to_ts, &
                      ' new_rp_var(3) old diff dist_to_TS'
      endif

     else
      if(new_distance_to_reactant.ge. distance_to_reactant) step=-step
      step_length=sqrt(sum(step(:)**2))
      over_step=abs(distance_to_reactant-new_distance_to_reactant).gt.abs(distance_to_ts)
      if( step_to_ts.and.over_step) then
       step_max=step_length*distance_to_ts/(distance_to_reactant-new_distance_to_reactant)
       call scale_step(step,step_max)
       call eliminate_step(step)
       step_length=sqrt(sum(step(:)**2))
       write(io_flepo,*) 'step scaled to distance to TS ', step_length
      endif
     endif

      de_mod=sum(step*grad_intern)

    if(grad_max_sphere.lt.max_comp_grad.and.grad_mean_sphere.lt.rms_grad) then
     io_pointonmep=openget_iounit(status='unknown',form='formatted', &
          file=trim(opt_data_dir)//'/point_on_mep.dat')
     write(io_pointonmep,*) s(:)%value
     mep_point=mep_point+1
     call returnclose_iounit(io_pointonmep) 
    endif

   end  function step_to_new_mep_point

  function new_distance(step)
     real(kind=r8_kind), intent(in):: step(:)
     real(kind=r8_kind):: new_distance,distance_to_reactant,contib
     integer(kind=i4_kind):: i
     distance_to_reactant=0.0_r8_kind
   do i=1,size(step)
    if(select_sphere_vars.and..not.s(i)%sphere) cycle
    contib=(s(i)%value+step(i)-s_reactant(i)%value)**2
    distance_to_reactant=distance_to_reactant+contib
     write(io_flepo,*) i,contib,s(i)%value+step(i),s_reactant(i)%value,'contib s s_r',s(i)%value,step(i)
   enddo
   ASSERT(distance_to_reactant.gt.small)
   new_distance=sqrt(distance_to_reactant)
  end function new_distance

  function distance_product(step)
     real(kind=r8_kind), intent(in):: step(:)
     real(kind=r8_kind):: distance_product,distance_to_product,contib
     integer(kind=i4_kind):: i
     distance_to_product=0.0_r8_kind
   do i=1,size(step)
    if(select_sphere_vars.and..not.s(i)%sphere) cycle
    contib=(s(i)%value+step(i)-s_product(i)%value)**2
    distance_to_product=distance_to_product+contib
     write(io_flepo,*) i,contib,s(i)%value+step(i),s_product(i)%value,'contib s s_p',s(i)%value,step(i)
   enddo
   ASSERT(distance_to_product.gt.small)
   distance_product=sqrt(distance_to_product)
  end function distance_product

    function cg_step() result(step1)
!Using Polak-Ribiere algorithm

    real(kind=r8_kind):: step1(n_internal)
    if(.not.dynamic) then
      if(step_counter==0) then
!Inicializaton of CG
        step1=-grad_intern
        grad_save=grad_intern
      else
        beta=dot_product(grad_intern,(grad_intern-grad_save))&
        /dot_product(grad_save,grad_save)
        step1=-grad_intern+beta*step
        grad_save=grad_intern
      endif
    endif  
   end function cg_step

    function sphere_noreturn(kk,step1) result(no_return)
    logical:: no_return
    real(kind=r8_kind),intent(in):: step1(:)
    real(kind=r8_kind):: q_dep2,sqt,sqtm,aa,bb,sr,sp,cc,det
    integer(kind=i4_kind),intent(in):: kk
    integer(kind=i4_kind):: i
    q=q_old+step1
    
  if(exist_product) then
    sqt=tsscan_rp_var**2
    sqtm=(1-tsscan_rp_var)**2

    aa=2*tsscan_rp_var-1
    bb=-2.0_r8_kind*(sqt*s_reactant(kk)%value-sqtm*s_product(kk)%value)

    sr=0.0_r8_kind
    sp=0.0_r8_kind
    do i=1,n_internal
     if(i.eq.kk.or.(select_sphere_vars.and..not.s(i)%sphere)) cycle
     sr=sr+(s_reactant(i)%value-q(i))**2
     sp=sp+(s_product(i)%value-q(i))**2
    enddo
    cc=sqt*(sr+s_reactant(kk)%value**2)-sqtm*(sp+s_product(kk)%value**2)
    det=bb**2-4*aa*cc
    no_return=.not.(det.ge.0.0_r8_kind)
    if(no_return) write(io_flepo,*) det,' no returnt to t_var_surface'
   else
    
    
    q_dep2=distance_to_reactant**2
    do i=1,n_internal
     if(i.eq.kk.or.(select_sphere_vars.and..not.s(i)%sphere)) cycle
     q_dep2=q_dep2-(s_reactant(i)%value-s(i)%value-step1(i))**2
    enddo
    no_return=.not.(q_dep2.ge.0.0_r8_kind)
    if(no_return) write(io_flepo,*) q_dep2,' no returnt to sphere'
   endif
    end function sphere_noreturn

    function dependent_var(kk) result(q_dep_var)
    real(kind=r8_kind):: q_dep_var,q_dep2,dq1,dq2
    integer(kind=i4_kind),intent(in):: kk
    integer(kind=i4_kind):: i

    q_dep2=distance_to_reactant**2
    do i=1,n_internal
     if(i.eq.kk.or.(select_sphere_vars.and..not.s(i)%sphere)) cycle
     q_dep2=q_dep2-(s_reactant(i)%value-q(i))**2
    enddo
    ASSERT(q_dep2.ge.0.0_r8_kind)
    dq1=abs(s(kk)%value-s_reactant(kk)%value-sqrt(q_dep2))
    dq2=abs(s(kk)%value-s_reactant(kk)%value+sqrt(q_dep2))
    if(dq1.lt.dq2) then
    q_dep_var=s_reactant(kk)%value+sqrt(q_dep2)
    write(io_flepo,*)'using dq1 sphere ',dq1, ' distance_to_reactant ',distance_to_reactant
    else
    q_dep_var=s_reactant(kk)%value-sqrt(q_dep2)
    write(io_flepo,*)'using dq2 sphere ',dq2, ' distance_to_reactant ',distance_to_reactant
    endif
    end function dependent_var

    function rp_dependent_var(kk) result(q_dep_var)
    real(kind=r8_kind):: q_dep_var,sr,sp,aa,bb,cc,sqt,sqtm,pp,rr
    real(kind=r8_kind):: q_dep_var2
    integer(kind=i4_kind),intent(in):: kk
    integer(kind=i4_kind):: i

    sqt=tsscan_rp_var**2
    sqtm=(1-tsscan_rp_var)**2

    aa=2*tsscan_rp_var-1
    bb=-2.0_r8_kind*(sqt*s_reactant(kk)%value-sqtm*s_product(kk)%value)

    sr=0.0_r8_kind
    sp=0.0_r8_kind
    do i=1,n_internal
     if(i.eq.kk.or.(select_sphere_vars.and..not.s(i)%sphere)) cycle
     sr=sr+(s_reactant(i)%value-q(i))**2
     sp=sp+(s_product(i)%value-q(i))**2
    enddo
    cc=sqt*(sr+s_reactant(kk)%value**2)-sqtm*(sp+s_product(kk)%value**2)
    ASSERT(bb**2-4*aa*cc.ge.0.0_r8_kind)
    
    q_dep_var=(-sqrt(bb**2-4*aa*cc)-bb)/(2*aa)
    q_dep_var2=(sqrt(bb**2-4*aa*cc)-bb)/(2*aa)
    if(abs(s(kk)%value-q_dep_var).gt.abs(s(kk)%value-q_dep_var2)) q_dep_var=q_dep_var2
    
    rr=sqrt(sr+(q_dep_var-s_reactant(kk)%value)**2)
    pp=sqrt(sp+(q_dep_var-s_product(kk)%value)**2)
    write(io_flepo,*)'tscan_rp_var dep_var ',  pp/(rr+pp), q_dep_var,kk,aa*q_dep_var**2+bb*q_dep_var+cc
    end function rp_dependent_var

    function rsrfo_step(alpha,a0_rfo,lambda_rfo,kk) result(step1)

    real(kind=r8_kind):: step1(n_internal),one=1.0_r8_kind
    integer(kind=i4_kind),intent(in):: kk
    real(kind=r8_kind),intent(out):: a0_rfo
    real(kind=r8_kind),intent(inout):: alpha,lambda_rfo
    integer(kind=i4_kind):: i,ifin,it
    logical :: no_return
    ifin=size(hesse_eigval)
        
        no_return=tsscan_sphere
        do  it=1,5
        do while((abs(step_length-r_curr)>eps1).or.(abs(ad)>eps2))
!Calculate d(step_length)/dalpha
                j=j+1
                denum=zero
                do i=1,ifin
!                  denum=denum+f_vec(i)**2/(hesse_eigval(i)-lambda_rs_i_rfo*alpha)**3
                   denum=denum+f_vec(i)**2/(hesse_eigval(i)-lambda_rfo*alpha)**3
                enddo
                write(io_flepo,*) 'step denum: ', denum,step_length,r_curr
                denum=denum*lambda_rfo/(one+alpha*step_length**2)
!End of calculate of d(step_length)/dalpha
                ad=(r_curr*step_length-step_length**2)/denum
                alpha=alpha+ad
!Solve general secular equation
                DPRINT j,alpha
                scl=zero
                scl(1,1)=one
                do i=1,ifin
                   scl(i+1,i+1)=alpha
                enddo
                call geigs(ahesse,scl,ahesse_eigval,ahesse_eigvec)
                a0_rfo=ahesse_eigvec(1,1)
                lambda_rfo=ahesse_eigval(1)
                step1(:ifin)=ahesse_eigvec(2:,1)/a0_rfo
                step_length=dsqrt(dot_product(step1(:ifin),step1(:ifin)))
                write (*,'(A,F15.8)')"current step_length=",step_length
         enddo

         DCALL show("New Augmented Hessian[RS-RFO]",ahesse(:,:))
         DCALL show("New Eigenvectors of Augmented Hessian[RS-RFO]",ahesse_eigvec(:,:))
         DCALL show("New Eigenvalues of Augmented Hessian[RS-RFO]",ahesse_eigval(:))
!End of microiterative process
!So, we have solution of RS-RFO equations
!Bakwards transformation of step.
         step1(:ifin)=matmul(hesse_eigvec,step1(:ifin))
           if(tsscan_sphere) then
            do i=ifin,kk,-1
             step1(i+1)=step1(i)
            enddo
            step1(kk)=0.0_r8_kind
           endif
          step1=matmul(tmat,step1)
          if(tsscan_sphere)  no_return=sphere_noreturn(sphere_dependent_var,step1)
          if(.not.no_return) exit
          if(.not.step_reset) r_curr=r_curr*0.7_r8_kind
         enddo
         end function rsrfo_step

    function qn_step(kk,f_vec) result(step1)
    real(kind=r8_kind):: step1(n_internal),one=1.0_r8_kind
    integer(kind=i4_kind),intent(in):: kk
    real(kind=r8_kind),intent(in):: f_vec(:)
    integer(kind=i4_kind):: i,ifin

    ifin=size(hesse_eigval)
    step1=zero
    do i=1,ifin
    step1(1:ifin) = step1(1:ifin) + (-one)*f_vec(i)*hesse_eigvec(:,i) / hesse_eigval(i)
    enddo
    if(tsscan_sphere) then
     do i=ifin,kk,-1
      step1(i+1)=step1(i)
     enddo
     step1(kk)=0.0_r8_kind
    endif
    step1 = matmul(tmat,step1)
    end function qn_step

        function Khait_step(kk) result(step1)
        real(kind=r8_kind):: step1(n_internal)
        integer(kind=i4_kind),intent(in):: kk
         integer(kind=i4_kind):: i,ifin

        ifin=size(hesse_eigval)
        step1(:ifin)=-hesse_eigvec(:,1)*f_vec(1)/abs(hesse_eigval(1))
        if(tsscan_sphere) then
         do i=ifin,kk,-1
          step1(i+1)=step1(i)
         enddo
         step1(kk)=0.0_r8_kind
        endif

        step1 = matmul(tmat,step1)
        end function Khait_step

    function ah_step(a0_ah,kk) result(step1)
    real(kind=r8_kind):: step1(n_internal)
    real(kind=r8_kind),intent(out) :: a0_ah
    integer(kind=i4_kind),intent(in):: kk
    integer(kind=i4_kind):: i,ifin
    
    a0_ah=ahesse_eigvec(1,1)
    step1=ahesse_eigvec(2:,1)/a0_ah

    !Bakwards transformation of step.
    ifin=size(hesse_eigval)
    step1(1:ifin)=matmul(hesse_eigvec,step1(:ifin))
        if(tsscan_sphere) then
     do i=ifin,kk,-1
      step1(i+1)=step1(i)
     enddo
     step1(kk)=0.0_r8_kind
    endif
    step1=matmul(tmat,step1)
    end function ah_step

           function Wales_step(small_eig,kk) result(step1)

           !Using Wales scheme for eliminate of shift parameters

            real(kind=r8_kind):: step1(n_internal),one=1.0_r8_kind
            real(kind=r8_kind), intent(in):: small_eig
            real(kind=r8_kind):: shift 
            integer(kind=i4_kind),intent(in):: kk
            integer(kind=i4_kind):: i,ifin
            ifin=size(hesse_eigval)

           step1=zero
           do i=1,ifin
              if(abs(hesse_eigval(i))<=small_eig) then ! 0.0001
                shift=zero
                write(io_flepo,'(A,4X,I4)') "Small eigenvalues:",i
                shift = (hesse_eigval(i) - sqrt(hesse_eigval(i)**2+2**2*f_vec(i)**2))/2
                write(io_flepo,'(A,4X,F12.8)') "Wales shift for this mode:",shift
                step1 = step1+(-one)*f_vec(i)*hesse_eigvec(:,i)/(hesse_eigval(i)-shift)
               else
                 step1 = step1+(-one)*f_vec(i)*hesse_eigvec(:,i)/abs((hesse_eigval(i)))
               endif
           enddo

           if(tsscan_sphere) then
             do i=ifin,kk,-1
              step1(i+1)=step1(i)
             enddo
             step1(kk)=0.0_r8_kind
            endif

           end function Wales_step

           function Greenstadt_step(kk) result(step1)
           real(kind=r8_kind):: step1(n_internal),one=1.0_r8_kind
           integer(kind=i4_kind):: kk,i,ifin

           ifin=size(hesse_eigval)
           step1=zero
           do i=1,ifin
              step1(1:ifin) = step1(1:ifin) +  (-one)*f_vec(i)*hesse_eigvec(:,i)/abs(hesse_eigval(i))
           enddo
           if(tsscan_sphere) then
             do i=ifin,kk,-1
              step1(i+1)=step1(i)
             enddo
             step1(kk)=0.0_r8_kind
            endif
           end function Greenstadt_step

         function rqn_step(kk,mu) result(step1)
         real(kind=r8_kind):: step1(n_internal),one=1.0_r8_kind
         real(kind=r8_kind),intent(in):: mu
         integer(kind=i4_kind),intent(in):: kk
         integer(kind=i4_kind):: i,ifin
         step1=zero
         ifin=size(hesse_eigval)
         do i=1,ifin
            step1 = step1+(-one)*f_vec(i)*hesse_eigvec(:,i)/(hesse_eigval(i)+mu)
         enddo

         if(tsscan_sphere) then
          do i=ifin,kk,-1
           step1(i+1)=step1(i)
          enddo
          step1(kk)=0.0_r8_kind
         endif

         step1 = matmul(tmat,step1)
         end function rqn_step

   function eigvec_grads(grad_intern,hesse_eigvec,tsscan_sphere) result(f_vec1)
    !Projecting gradient on eigenvectors of current hessian    
    use gradient_module, only: grad_sphere
    real(kind=r8_kind), intent(in):: grad_intern(:),hesse_eigvec(:,:)
    real(kind=r8_kind):: f_vec1(size(grad_intern))
    logical, intent(in):: tsscan_sphere
    integer(kind=i4_kind):: i
    f_vec1=zero
    if(tsscan_sphere) then
    do i=1,size(grad_sphere)
       f_vec1(i) = dot_product(hesse_eigvec(:,i),grad_sphere)
    enddo
    else
    do i=1,size(grad_intern)
       f_vec1(i) = dot_product(hesse_eigvec(:,i),grad_intern)
    enddo
    endif

   end function eigvec_grads

  end subroutine compute_step

  subroutine hebden(mu,lambda,r_curr)
    !-------------------------------------------------------------
    !Purpose: solves the equation
    !     one-sum_i[(f_vec(i)/(hesse_eigval(i)-lambda+mu))]^2/r_curr=zero
    !         iteratively via a Newton-Method
    !     hesse_eigval(i) = Eigenvalues of current Hessian (n_internal)
    !         f_vec(i)    = Component of gradients along eigenmodes of the
    !                       current Hessian
    !         lambda,mu       = shift parameter's. lambda<=zero;mu>=zero
    !
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use type_module
    use math_module,only:one,zero,small
    use hesse_module,only:hesse_eigval
    implicit none
    !------------ Declaration of formal parameters ---------------
        real(kind=r8_kind),intent(in)   :: lambda,r_curr
        real(kind=r8_kind),intent(out)  :: mu
        !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)           :: i
    integer(kind=i4_kind),parameter :: maxiter=1000
        real(kind=r8_kind)              :: ad,mu_prev
    real(kind=r8_kind)              :: delta = 0.0002_r8_kind
    !------------ Executable code ------------------------------------
!    mu=max(-hesse_eigval(1),0.0_r8_kind)+delta 
    mu=-hesse_eigval(1)+delta 
    mu_prev=mu+delta
    DPRINT "===================="
    DPRINT "Compute mu"
    DPRINT "Iter                mu"
  hebd: do i=1,maxiter
           if(abs(dfunc(mu,lambda))<=small) then 
             !Using Secant method
             DPRINT  "Using Secant method",i,mu,mu_prev,func(mu_prev,lambda)
             ad=-func(mu,lambda)*(mu-mu_prev)/(func(mu,lambda)-func(mu_prev,lambda))
            else
             ad=-func(mu,lambda)/dfunc(mu,lambda) 
             DPRINT 'func dfunc',i,func(mu,lambda),dfunc(mu,lambda)
           endif
!           print*, i,mu,func(mu,lambda),ad,small,mu-mu_prev, &
!                    func(mu,lambda)-func(mu_prev,lambda)
           if((abs(ad))<small) exit hebd
           mu_prev=mu
           mu=mu+ad
        end do hebd
    DPRINT  "End compute mu",func(mu,lambda),mu,lambda
    contains
    
    function func(mu,lambda)
    use math_module, only : zero,one
    implicit none
    real(kind=r8_kind),intent(in)     :: mu,lambda
    real(kind=r8_kind)                :: func
    integer(kind=i4_kind)             :: i
    func=zero
    do i=1,size(hesse_eigval)
       func=func+(f_vec(i)/(hesse_eigval(i)-lambda+mu))**2
    enddo
    func=dsqrt(func)
!    print*, func,r_curr
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
    do i=1,size(hesse_eigval)
       denum=denum+(f_vec(i)/(hesse_eigval(i)-lambda+mu))**2
       num=num+(f_vec(i))**2/(hesse_eigval(i)-lambda+mu)**3
    enddo
    denum=r_curr*dsqrt(denum)
    dfunc=num/denum
    end function dfunc        

  end subroutine hebden

  subroutine desrd(method,n_negative,desired)
    !---------------------------------------------------------------
    ! Purpose: check has hessian desired structure
    !---------------------------------------------------------------
    !------------ Modules used -------------------------------------
    use type_module
    use opt_data_module, only: io_flepo
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer(kind=i4_kind),intent(in)  :: n_negative
    character(len=*),intent(in)       :: method
    logical,intent(out)               :: desired
    !------------ Executable code ----------------------------------        
    if(n_negative>0_i4_kind) then
      write(io_flepo,'(A,A)')"Method:",method 
      write(io_flepo,'(A)')"Current hessian has UNDESIRED structure:has a superfluous number of negative eigenvalues"
      write(io_flepo,'(A,I10)')"Number of negative eigenvalues is:",n_negative
      desired=.false.
    else
      write(io_flepo,'(A,A)')"Method:",method 
      write(io_flepo,'(A)')"Current hessian has DESIRED structure:all eigenvalues of hesse matrix is positive"
      desired=.true.
    endif
 end subroutine desrd
  
 subroutine scale_step(step,step_max) 
    !---------------------------------------------------------------
    ! Purpose: do a simple step limitation to 0.5 au/0.3 rad.
    !          Many authors suggest that a constant step restriction
    !          rather than a calculated trust radius is sufficient.
    !          See e.g. J.Baker, B.H.Schlegel ...
    !          12/05 VVP:
    !          For Augmented Hessian method please see [8].
    !---------------------------------------------------------------
    real(kind=r8_kind), intent(in):: step_max
    real(kind=r8_kind), intent(inout):: step(:)
    real(kind=r8_kind):: step_length
       step_length=sqrt( sum( step**2) )
       step = step/step_length*step_max
       write(io_flepo,*)"          step is scaled to maximum step length",sqrt(sum( step**2) )
 end subroutine scale_step

 subroutine restrict(step_max, step)
    !---------------------------------------------------------------
    ! Purpose: do a simple step limitation to 0.5 au/0.3 rad.
    !          Many authors suggest that a constant step restriction
    !          rather than a calculated trust radius is sufficient.
    !          See e.g. J.Baker, B.H.Schlegel ...
    !          12/05 VVP:
    !          For Augmented Hessian method please see [8].
    !---------------------------------------------------------------
    use math_module, only: abs_value
    !    use opt_data_module, only: step_max
    implicit none
    real(kind=r8_kind), intent(in)    :: step_max
    real(kind=r8_kind), intent(inout) :: step(:)
    ! *** end of interface ***

    if (abs_value(step)>=step_max.or.new_mep_point) then
       write(io_flepo,*)"          step will be scaled to maximum step length"
       step = step/abs_value(step)*step_max
    endif
 end subroutine restrict


 subroutine tr(geo_loop,de_mod,less_t_r_curr,r_curr,e_cur,accepted)
    !-------------------------------------------------------------
    !  Purpose: Calculate Trust Radius at each iteration.
    !-------------------------------------------------------------
    !------------ Modules used -----------------------------------
    use math_module,only:two
    use gradient_module,only:energy,grad_mean_square,grad_max_comp
    use opt_data_module,only:step_max
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in)  :: geo_loop
    real(kind=r8_kind),   intent(in)  :: de_mod
    logical,              intent(in)  :: less_t_r_curr
    real(kind=r8_kind),   intent(out) :: r_curr,e_cur
    logical,              intent(out) :: accepted
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)                :: r,de_obj,e_curr
!    real(kind=r8_kind),save           :: e_prev
    real(kind=r8_kind)                :: r1e,r1i,rue,rui
    real(kind=r8_kind),parameter      :: r_min=0.75_r8_kind,r_good=0.85_r8_kind
    real(kind=r8_kind),parameter      :: r_0=0.3_r8_kind,sf=two,gr_small=0.0001_r8_kind
    real(kind=r8_kind),parameter      :: en_small=0.00000015_r8_kind
    !------------ Executable code --------------------------------  
    write(io_flepo,'(27X,A)')"--- Trust Radius ---"
    if((grad_max_comp<=gr_small).and.(grad_mean_square<=gr_small).and.(step_counter/=0)) then
      write(io_flepo,2200)
      r_curr=r_curr
      e_cur=energy
      write(io_flepo,2100) r_curr
      accepted=.true.
    else
!r1e=0.75;r1i=1.45;rue=1.30;rui=0.85    
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
          write(io_flepo,2200)
          r_curr=r_curr
          e_cur=energy
          write (io_flepo,1600) de_obj
          write (io_flepo,1700) de_mod
          write (io_flepo,2300) en_small
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
          elseif(.not.step_reset) then
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

    if((r<=zero).or.(r>=1000_r8_kind)) then
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
    DPRINT "Counter",step_counter
    write(io_flepo,'(27X,A)') "--------------------"
1600 format("Actual change energy is      [DFT] :",F15.8)
1700 format("Predictable change energy is [LQA] :",F15.8)
1800 format("Ratio=[DFT]/[LQA]                  :",F15.8)
1900 format("Current value of TR is             :",F15.8)
2000 format("Initial value of TR is             :",F15.8)
2100 format("Frozen value of TR is              :",F15.8)
2200 format(27X,"--- Is switch off ---")
2300 format("Minimal energy deviation           :",F15.8)
 end subroutine tr
  
    
 function ang(step,grad)
  !-------------------------------------------------------------
  ! Compute angle(in degree) between step and force vectors
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
    use math_module, only : convert1,zero
    USE_DEBUG! only: isNaN(x)
    implicit none
  !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in) :: step(:)
    real(kind=r8_kind),intent(in) :: grad(:)
    real(kind=r8_kind)            :: ang ! in degrees
    ! *** end of interface ***

    real(kind=r8_kind) :: steplen, gradlen

    steplen = sqrt(dot_product(step,step))
    gradlen = sqrt(dot_product(grad,grad))

    if(steplen*gradlen /= 0.0)then
       ! angle in degrees:
       ang=acos(-dot_product(grad,step)/(steplen*gradlen)) &
           * convert1
    else
       ang=zero
    endif
 end function ang

 subroutine write_step(meth,step,f_vec,hesse_eigval,tht)
  !-------------------------------------------------------------
  !This subroutine write step information in flepo file.
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
  use opt_data_module,only:n_internal
  use math_module,only:convert1
  use gradient_module,only:grad_intern,grad_sphere
  implicit none
  !------------ Declaration of formal parameters ---------------
  character(len=*),intent(in)       :: meth
  real(kind=r8_kind),intent(in)     :: step(:),f_vec(:),hesse_eigval(:)
  real(kind=r8_kind),intent(in)     :: tht
  !------------ Declaration of local variables -----------------
  character(len=4)                  :: std="--- "
  character(len=80)                 :: title,en
  integer(kind=i4_kind)             :: ln,i,ii
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
         if(.not.cart_coordinates) then
            write(io_flepo,*)"   mode    grad. intern.     grad. comp.       eigenvalue     type       step"
            if(tsscan_sphere) then
               do i=1,size(hesse_eigval)
                  ii=i
                  if(i.ge.sphere_dependent_var) ii=i+1
                  if(i.eq.sphere_dependent_var)write(io_flepo,'(A)') '___________________________________________'
                  write(io_flepo,'(3x,i3,a,5x,es12.5,5x,es12.5,6x,es12.5,4x,i3,5x,es12.5)')&
                       i,'m',grad_sphere(i),f_vec(i),hesse_eigval(i),i,step(ii)
               enddo
            else
               do i=1,size(hesse_eigval)
                  write(io_flepo,'(3x,i3,a,5x,es12.5,5x,es12.5,6x,es12.5,4x,i3,5x,es12.5)')&
                       i,'m',grad_intern(i),f_vec(i),hesse_eigval(i),i,step(i)
               enddo
            endif
         else
            write(io_flepo,*)"   mode    grad. intern       eigenvalue.     type       step"
            do i=1,n_internal
               write(io_flepo,'(3x,i3,a,5x,es12.5,5x,es12.5,4x,i3,5x,es12.5)')&
                    i,'m',grad_intern(i),hesse_eigval(i),i,step(i)
            enddo
         end if
      endif

      write(io_flepo,'(27X,A)') en
      write(io_flepo,'(A,1X,A,F15.8)')"Theta(Degree) between step and forces ",":",tht

 end subroutine write_step

 subroutine eliminate_step(step)
  !-------------------------------------------------------------
  ! This subroutine sets tiny components of step(:) to plain zero
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
  use opt_data_module,only:n_internal
  use math_module,only:zero
  implicit none
  !------------ Declaration of formal parameters ---------------
  real(kind=r8_kind),intent(inout)  :: step(:)
  !------------ Declaration of local variables -----------------
  integer(kind=i4_kind)             :: i
  real(kind=r8_kind),parameter      :: sm_comp=0.000001
  !------------ Executable code --------------------------------  
  ! Eliminate small components of step
  do i=1,n_internal
     if(abs(step(i))<=sm_comp) then
        step(i)=zero
        write(io_flepo,'(A,4X,I4)')"Eliminating component of step:",i
     endif 
  enddo
 end subroutine eliminate_step


  subroutine gdiis(geo_loop,geom,a0,meth,q_new,g_new,do_gdiis)
  !-------------------------------------------------------------
  !This subroutine perform GDIIS extrapolation/interpolation
  !-------------------------------------------------------------
  !------------ Modules used -----------------------------------
  use opt_data_module,only:n_internal
  use type_module
  use math_module
  use hesse_module, only: hesse_inv
  use gradient_module,only:grad_mean_square,grad_max_comp,energy,grad_intern
  USE_DEBUG
! !This type is very helpful for GDIIS  
! type gdiis_point
!  integer(kind=i4_kind)             :: num_pt
!  real(kind=r8_kind),pointer        :: q_gdiis(:)
!  real(kind=r8_kind),pointer        :: g_gdiis(:)
!  real(kind=r8_kind)                :: ener
! end type gdiis_point
  !------------ Declaration of formal parameters ---------------
  integer(kind=i4_kind),intent(in)  :: geo_loop
  real(kind=r8_kind),intent(inout)  :: geom(:)
  character(len=*),intent(in)       :: meth
  real(kind=r8_kind),intent(in)     :: a0
  real(kind=r8_kind),intent(out)    :: q_new(:),g_new(:)
  logical,intent(out)               :: do_gdiis

! see globals module vars:
! integer(kind=i4_kind),parameter    :: max_point=4 !initital value
! integer(kind=i4_kind),parameter    :: ns_point=4
! type(gdiis_point), save            :: gdiis_work(max_point)
! integer(kind=i4_kind),save         :: num_stpt=0
! integer(kind=i4_kind),save         :: calls=0

  integer(kind=i4_kind)              :: i,alloc_stat,j,rep_point
  real(kind=r8_kind)                 :: length,max_length,grad_sm,err
  real(kind=r8_kind),allocatable     :: gwm(:,:), cg(:),e_i(:),e_j(:),rs(:)
  real(kind=r8_kind),allocatable     :: delta_q(:)
  character(len=80)                  :: gdiis_method
  logical,save                       :: accpt
  
  accpt=.false.
  do_gdiis=.false.
  grad_sm = 0.0
!three various scheme for calculate of error vector. For QN.
  select case(meth)
         case("qn")
!Suggested by Pulay(original
         gdiis_method="energy"
         grad_sm=0.01_r8_kind
!Next two suggested by Werner&&Forkas. For (RFO|AH).
         case("rfo")
         gdiis_method="energy"
         grad_sm=0.001_r8_kind
         case("ah")
         gdiis_method="geometries"
         grad_sm=0.001_r8_kind
         case default
         write(io_flepo,'(A,A)') "This method does not exist", meth
  end select
  ASSERT (grad_sm /= 0.0)
  
  write(io_flepo,'(27X,A)')"--- GDIIS ---"
  if( .not.gdiis_storage_initialized )then
    write(io_flepo,'(A)')"Initialization of GDIIS"
    do i=1,max_point
    allocate(gdiis_work(i)%q_gdiis(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    allocate(gdiis_work(i)%g_gdiis(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    enddo
    gdiis_storage_initialized = .true.
  endif
  
  if(num_stpt<max_point) then
    if((grad_mean_square<=grad_sm).and.(grad_max_comp<=grad_sm).and.(meth.eq."qn")) accpt=.true.
    if((abs(a0)>=0.9995_r8_kind).and.(meth.eq."rfo")) accpt=.true.
    if(accpt) then
      write(io_flepo,'(A)')"Current geometry is acceptable for future GDIIS interpolation"
      num_stpt=num_stpt+1
      gdiis_work(num_stpt)%num_pt=geo_loop
      gdiis_work(num_stpt)%q_gdiis=geom
      gdiis_work(num_stpt)%g_gdiis=grad_intern
      gdiis_work(num_stpt)%ener=energy
    endif
    write (io_flepo,'(A,2X,I4)')"Numbers of stored points               :",num_stpt
!   calls=calls+1
  else
  !Compute index of replacing point
    allocate(delta_q(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    max_length=zero
    ! rep_point is supposed to be set in the following loop in any case
    ! initialisation for making the compiler happy
    rep_point = -num_stpt
    do i=1,num_stpt
          delta_q=gdiis_work(i)%q_gdiis-geom
          if(.not. cart_coordinates) call correct_step(delta_q)
          length=abs_value(delta_q)
          if(length>=max_length) then
             max_length=length
             rep_point=i
          endif
          write(io_flepo,'(A,I4,A,F15.8)')"Certain distance between current point and ",gdiis_work(i)%num_pt,&
          " stored geometries:",length
    enddo
    write(io_flepo,'(A,F15.8)') "Maximal distance:",max_length
    write(io_flepo,'(A,4X,I4,A,I4)')"Replacing point :",gdiis_work(rep_point)%num_pt,"->",geo_loop
    gdiis_work(rep_point)%q_gdiis=geom
    gdiis_work(rep_point)%g_gdiis=grad_intern
    gdiis_work(rep_point)%num_pt=geo_loop
    gdiis_work(rep_point)%ener=energy
    deallocate(delta_q,STAT=alloc_stat)
    ASSERT(alloc_stat==zero)
  endif
!Write information in flepo file
 if(num_stpt>=ns_point) then
   write(io_flepo,'(A)')"------------------------------------------------------------------"
   write(io_flepo,*)"geo_loop   max.comp.of.grad        rms.grad.       energy"
   write(io_flepo,'(A)')"------------------------------------------------------------------"
   do i=1,num_stpt
      write(io_flepo,'(3x,i3,8x,F15.8,5x,F15.8,2X,F15.8)')&
      gdiis_work(i)%num_pt, &
      maxval(abs(gdiis_work(i)%g_gdiis)), &
      abs_value(gdiis_work(i)%g_gdiis / sqrt(real(n_internal, r8_kind))), &
      gdiis_work(i)%ener
   enddo
   write(io_flepo,'(A)')"------------------------------------------------------------------"
 endif
 if(num_stpt>=ns_point) then
!Construct GDIIS works matrix and needed arraies 
   DCALL show("Current H^(-1)",hesse_inv(:,:))
   allocate(gwm(num_stpt+1,num_stpt+1),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   gwm=one
   allocate(cg(num_stpt+1),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   allocate(rs(num_stpt+1),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   rs=zero;rs(num_stpt+1)=one
   allocate(e_j(n_internal),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   allocate(e_i(n_internal),STAT=alloc_stat)
   ASSERT(alloc_stat==0)

   do i=1,num_stpt
     do j=1,num_stpt
     select case(gdiis_method)
     case("gradient")
     e_i=gdiis_work(i)%g_gdiis
     e_j=gdiis_work(j)%g_gdiis
     gwm(i,j)=dot_product(e_i,e_j)
     case("energy")
     e_i=-matmul(hesse_inv,gdiis_work(i)%g_gdiis)
     e_j=-matmul(hesse_inv,gdiis_work(j)%g_gdiis)
     gwm(i,j)=dot_product(e_i,e_j)
     case("geometries")
     e_i=gdiis_work(i)%g_gdiis
     e_j=gdiis_work(j)%g_gdiis
     gwm(i,j)=dot_product(e_i,matmul(hesse_inv,e_j))
     case default
       ABORT('no such method')
     end select
     enddo
   enddo

   gwm(num_stpt+1,num_stpt+1)=zero
!End of construct works matrix and needed arraies
   DCALL show("GDIIS works matrix",gwm(:,:))
   call gauss_solve(gwm,rs,cg)
   DCALL show("GDIIS coefficients",cg(1:num_stpt))
   DPRINT "Lambda",cg(num_stpt+1)
!Construct new geometry and new gradient
   q_new=zero;g_new=zero
   do i=1,num_stpt
     select case(gdiis_method)
     case("gradient")
     e_i=gdiis_work(i)%g_gdiis
     q_new=q_new+cg(i)*gdiis_work(i)%q_gdiis!+cg(i)*e_i
     g_new=g_new+cg(i)*gdiis_work(i)%g_gdiis!+cg(i)*e_i
     case("energy")
     e_i=-matmul(hesse_inv,gdiis_work(i)%g_gdiis)
     q_new=q_new+cg(i)*gdiis_work(i)%q_gdiis+cg(i)*e_i
     g_new=g_new+cg(i)*gdiis_work(i)%g_gdiis+cg(i)*e_i
     case("geometries")
     e_i=gdiis_work(i)%g_gdiis
     q_new=q_new+cg(i)*gdiis_work(i)%q_gdiis!+cg(i)*e_i
     g_new=g_new+cg(i)*gdiis_work(i)%g_gdiis!+cg(i)*e_i
     case default
       ABORT('no such method')
     end select
   enddo
   DCALL show("GDIIS geometry",q_new(:))
   DCALL show("GDIIS gradient",g_new(:))
   allocate(delta_q(n_internal),STAT=alloc_stat)
   ASSERT(alloc_stat==zero)
   delta_q=q_new-geom
   DCALL show("Delta Q",delta_q(:))
   if(.not. cart_coordinates) call correct_step(delta_q)

   ! Compute magnitude of the residual vector:
   err = 0.0
   do i=1,num_stpt
     do j=1,num_stpt
     err=err+cg(i)*gwm(i,j)*cg(j)
     enddo
   enddo
   do_gdiis=check_gdiis(cg,num_stpt,delta_q,err)
   deallocate(cg,STAT=alloc_stat);ASSERT(alloc_stat==0)
   deallocate(gwm,STAT=alloc_stat);ASSERT(alloc_stat==0)
   deallocate(rs,STAT=alloc_stat);ASSERT(alloc_stat==0)
   deallocate(e_i,STAT=alloc_stat);ASSERT(alloc_stat==0)
   deallocate(e_j,STAT=alloc_stat);ASSERT(alloc_stat==0)
   deallocate(delta_q,STAT=alloc_stat);ASSERT(alloc_stat==0)
 endif
!  calls=calls+1
 write(io_flepo,'(27X,A)')"-------------"
  
  
  
  contains  
  
  function gdiis_angle(cos_anbs,max_point)
  real(kind=r8_kind),intent(in)       :: cos_anbs
  integer(kind=i4_kind),intent(in)    :: max_point
  logical                             :: gdiis_angle
  
  gdiis_angle=.false.
  select case(max_point) 
         case(2_i4_kind) 
             if(cos_anbs>0.97_r8_kind) gdiis_angle=.true.
         case(3_i4_kind)
             if(cos_anbs>0.84_r8_kind) gdiis_angle=.true.
         case(4_i4_kind)
             if(cos_anbs>0.71_r8_kind) gdiis_angle=.true.
         case(5_i4_kind)
             if(cos_anbs>0.67_r8_kind) gdiis_angle=.true.
         case(6_i4_kind)
             if(cos_anbs>0.62_r8_kind) gdiis_angle=.true.
         case(7_i4_kind)
             if(cos_anbs>0.56_r8_kind) gdiis_angle=.true.
         case(8_i4_kind)
             if(cos_anbs>0.49_r8_kind) gdiis_angle=.true.
         case(9_i4_kind)
             if(cos_anbs>0.41_r8_kind) gdiis_angle=.true.
         case(10_i4_kind)
                                        gdiis_angle=.true.
  end select
  end function gdiis_angle

   function check_gdiis(cg,num_stpt,delta_q,err)
 
   use hesse_module, only: hesse,hesse_inv
   use gradient_module, only: grad_intern
   use opt_data_module
   use math_module
   
   real(kind=r8_kind),intent(in)    :: cg(:)
   integer(kind=i4_kind),intent(in) :: num_stpt
   real(kind=r8_kind),intent(in)    :: delta_q(:)
   real(kind=r8_kind),intent(in)    :: err
   logical                          :: check_gdiis
 
   real(kind=r8_kind)               :: n_sum_cg,p_sum_cg,sum,cos_anbs
   real(kind=r8_kind),allocatable   :: help(:)
   integer(kind=i4_kind)            :: i,alloc_stat
   
   allocate(help(n_internal),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   check_gdiis=.true.
   n_sum_cg=zero;p_sum_cg=zero
   do i=1,num_stpt
      if(cg(i)>zero) then
        p_sum_cg=p_sum_cg+cg(i)
      else
        n_sum_cg=n_sum_cg+cg(i)
      endif
   enddo
   if((p_sum_cg>15_r8_kind).or.(n_sum_cg<-15_r8_kind))            check_gdiis=.false.
   sum=p_sum_cg+n_sum_cg
   if(abs(sum-one)>=small)                      check_gdiis=.false.
   if(maxval(abs(cg/err))>=100000000_r8_kind)   check_gdiis=.false.
 
   de_mod=dot_product(grad_intern,delta_q)+half*dot_product(delta_q,matmul(hesse,delta_q))
   if(de_mod>=zero)                             check_gdiis=.false.
   help=-matmul(hesse_inv,grad_intern)
   if(abs_value(delta_q)/abs_value(help)>10_r8_kind) check_gdiis=.false.
   cos_anbs=dot_product(delta_q,help)/(abs_value(delta_q)*abs_value(help))
   if(.not.gdiis_angle(cos_anbs,num_stpt)) check_gdiis=.false.
   write(io_flepo,1000)"Sum of negative coefficients           :",n_sum_cg
   write(io_flepo,1000)"Sum of positive coefficients           :",p_sum_cg
   write(io_flepo,1000)"Sum of all GDIIS coefficients          :",sum
   write(io_flepo,1100)"Magnitude of residium vector [e^2]     :",err
   write(io_flepo,1000)"Relation between GDIIS and QN steps    :",abs_value(delta_q)/abs_value(help)
   write(io_flepo,1000)"Cos(theta) between GDIIS and QN steps  :",cos_anbs
   deallocate(help,STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   return
1000 format(A,2X,F15.8)
1100 format(A,2X,ES15.4)
   end function check_gdiis
 
 end subroutine gdiis
  
 
 

 subroutine correct_step(delta_q)
   !-------------------------------------------------------------
   !This subroutine correct step
   !-------------------------------------------------------------
   !------------ Modules used -----------------------------------
   use math_module
   use coortype_module
   use opt_data_module
   !------------ Declaration of formal parameters ---------------
   real(kind=r8_kind),intent(inout)   :: delta_q(:)
   !------------ Declaration of local variables -----------------
   integer(kind=i4_kind)              :: i
   !------------ Executable code --------------------------------  
   do i=1,n_internal
      if(s(i)%typ==d_angle) then
        if(delta_q(i)>pi)  delta_q(i)=-two*pi+delta_q(i)
        if(delta_q(i)<-pi) delta_q(i)= two*pi+delta_q(i)
      endif
   enddo     
 end subroutine correct_step

 subroutine quartic(q_prev,q_curr,g_prev,g_curr,e_prev,e_curr,q_new,e_new,noline)
   !-------------------------------------------------------------
   !This subroutine perfrom partial(inexact) line search by quartic polynom.
   !-------------------------------------------------------------
   !------------ Modules used -----------------------------------
   use type_module
   use math_module
   use opt_data_module
   USE_DEBUG
   !------------ Declaration of formal parameters ---------------
   implicit none
   real(kind=r8_kind),intent(in)  :: q_prev(:),q_curr(:),g_prev(:),g_curr(:)
   real(kind=r8_kind),intent(in)  :: e_prev,e_curr
   real(kind=r8_kind),intent(out) :: q_new(:)
   real(kind=r8_kind),intent(out) :: e_new
   logical,intent(inout)          :: noline
   !------------ Declaration of local variables -----------------
   real(kind=r8_kind),allocatable :: delta_q(:)
   real(kind=r8_kind),parameter   :: stepmx=one,tol=1d-9,eps=1d-14
   real(kind=r8_kind)             :: convl,e1a,e1b,a1,a2,a3,sss
   real(kind=r8_kind)             :: p,test,a4,xq,fq,xr,fr,fea,feb
   real(kind=r8_kind)             :: fe1a,fe1b,tt,ss,a0,dx
   real(kind=r8_kind)             :: x1,x2,pp,qq,dd
   integer(kind=i4_kind)          :: i,alloc_stat
   !------------ Executable code --------------------------------  
   write(io_flepo,'(27X,A)')"---- Line search ---"
   noline=.false.     
   convl = zero
!all in terms of vector from previously to current steps
   e1a=zero
   e1b=zero
   allocate (delta_q(n_internal),STAT=alloc_stat)
   ASSERT(alloc_stat==0)
   delta_q=q_curr-q_prev
   if(.not. cart_coordinates) call correct_step(delta_q)
   e1a=dot_product(delta_q,g_prev)
   e1b=dot_product(delta_q,g_curr)
   write(io_flepo,1000)"Previous   energy                    :",e_prev
   write(io_flepo,1000)"Current    energy                    :",e_curr
   write(io_flepo,1000)"Projected previuos   gradient on step:",e1a
   write(io_flepo,1000)"Projected    current gradient on step:",e1b
!when energies are equal, dont rely on quartic fit to do anything
!reasonable at all
   if(abs((e_curr-e_prev)/e_prev)<=eps) then
     q_new = (q_curr+q_prev)*.5d0
     e_new = (e_curr+e_prev)*.5d0
     write(io_flepo,'(A)')"Energies are equal"
     return
   endif
   a1 = .25_r8_kind*(6_r8_kind*e_curr-6_r8_kind*e_prev-e1b-e1a)
   a3 = -6_r8_kind *(2_r8_kind*e_curr-2_r8_kind*e_prev-e1b-e1a)

   sss = e1b-e1a
   if(sss<=zero) then
!gradient decreases, i.e. negative definite
     if(e_prev<=e_curr) then
       q_new = q_prev
     else
       q_new = q_curr
     end if
      deallocate(delta_q)
   return
   endif

   p = a3**2_r8_kind/48_r8_kind
   test = sss**2-4_r8_kind*p
!!$   if(test<=zero) then
   if(test<=1.0e-8_r8_kind) then !??
!must use quadratic only
#if 1
a1=e1a
a3=2*(e_prev-e_curr)+e1b+e1a
a2=e_curr-e_prev-e1a-a3
pp=2*a2/3.0_r8_kind/a3
qq=a1/3.0_r8_kind/a3
dd=pp**2/4-qq
if(dd.ge.0.0_r8_kind) then
x1=-pp/2+sqrt(dd)
x2=-pp/2-sqrt(dd)
if(x1.gt.0.0_r8_kind.and.x2.gt.0.0_r8_kind) then
xr=min(x1,x2)
e_new=e_prev+a1*xr+a2*xr**2+a3*xr**3
elseif(x1.gt.0.0_r8_kind.and.x1.lt.1.0_r8_kind) then
xr=x1
e_new=e_prev+a1*xr+a2*xr**2+a3*xr**3
elseif(x2.gt.0.0_r8_kind.and.x2.lt.1.0_r8_kind) then
xr=x2
e_new=e_prev+a1*xr+a2*xr**2+a3*xr**3
else
noline=.true.
e_new=min(e_prev,e_curr)
if(e_prev.gt.e_curr) then
xr=0.0_r8_kind
else
xr=1.0_r8_kind
endif
endif
endif
   write (io_flepo,*)"New energy(a1x+a2x**2+a3x**3), xr             :",e_new,xr
 do i=1,n_internal
  q_new(i) = q_prev(i)+xr*delta_q(i)
 enddo
      deallocate(delta_q)
      call quartic_shutdown(noline)
      return
#endif
     a3 = zero
     a4 = zero
     a2 = sss
     a1 = e_curr-e_prev
     noline = .true.    !(3)
   else
!use quartic
     test = dsqrt(test)
     a2 = .5_r8_kind*(sss+test)
     a4 = 12_r8_kind*(sss-test)
   endif
!find gradient=0 by secant method
   if(e_prev<=e_curr) then
     xq = -.5_r8_kind
     fq = e1a
   else
     xq = .5_r8_kind
     fq = e1b
   endif
   xr = -a2/a1
   fr = a1+xr*(a2+.5_r8_kind*xr*(a3+xr*a4/3_r8_kind))
100   dx = (xq-xr)*fr / (fr-fq)
   xq = xr
   fq = fr
   xr = xr+dx
   fr = a1+xr*(a2+.5_r8_kind*xr*(a3+xr*a4/3_r8_kind))
   if(dx>tol.and.abs(fr-fq)>tol.and.abs(xr)<stepmx) goto 100
   if(abs(xr)>stepmx) then
     noline = .true.     !(2)
     xr = sign(stepmx,xr)
   endif
   a0 = .5_r8_kind*(e_prev+e_curr)-a2/8_r8_kind-a4/(24_r8_kind*16_r8_kind)
   e_new = a0+xr*(a1+(xr/2_r8_kind)*(a2+(xr/3_r8_kind)*(a3+(xr/4_r8_kind)*a4)))
   write(io_flepo,'(A)')"Coefficients of quartic polynom:"
   write(io_flepo,'(A)')"------------------"
   write(io_flepo,1000)"a0=",a0
   write(io_flepo,1000)"a1=",a1
   write(io_flepo,1000)"a2=",a2
   write(io_flepo,1000)"a3=",a3
   write(io_flepo,1000)"a4=",a4
   write(io_flepo,'(A)')"------------------"
   write (io_flepo,1000)"New energy                          :",e_new
!construct new x,g
!g interpolated linearly for orthogonal directions
!g annihilated in search direction if search has been 'satisfactory'
   tt = zero
   ss = zero
   do i = 1, n_internal
      ss = ss + (delta_q(i))**2
      tt = tt + (delta_q(i))*(.5d0*(g_curr(i)+g_prev(i))+xr*(g_curr(i)-g_prev(i)))
   enddo
      tt = -tt/ss
      if(noline) tt = zero
      do i = 1, n_internal
         q_new(i) = .5d0*(q_curr(i)+q_prev(i))+xr*delta_q(i)
      enddo
      convl=sqrt(ss)*(xr-.5d0)
      if (xr.lt.0) convl=sqrt(ss)*(xr+.5d0)
      convl = abs(convl)
      xq = -.5d0
      fea = a0+xq*(a1+(xq/2)*(a2+(xq/3)*(a3+(xq/4)*a4)))
      fe1a = a1+xq*(a2+(xq/2)*(a3+(xq/3)*(a4)))
      xq = .5d0
      feb = a0+xq*(a1+(xq/2)*(a2+(xq/3)*(a3+(xq/4)*a4)))
      fe1b = a1+xq*(a2+(xq/2)*(a3+(xq/3)*(a4)))
      DCALL show("q_new",q_new(:))
      deallocate(delta_q)
      call quartic_shutdown(noline)
      return
1000  format(A,F20.8)

      contains 
       
       subroutine quartic_shutdown(noline)
       logical,intent(in)      :: noline
       if(noline) then
         write(io_flepo,'(A)')"Line search is unsuccessfull"
       else 
         write(io_flepo,'(A)')"Line search is successfull"
       endif  
       return
       end subroutine quartic_shutdown
      end subroutine quartic

  subroutine step_module_persistent_state(act,iounit)
    !
    ! Save/Restore module vars that need to persist
    ! across optimizer calls if optimizer is compiled
    ! as a standalone program or for "optimizer_only" runs
    !
    implicit none
    character(len=*), intent(in) :: act
    integer(i4_kind), intent(in) :: iounit
    ! *** end of interface ***

    integer(i4_kind) :: i,n_coords

    ! global module vars:
    namelist /step_module_state/ de_mod      &
                               , e_save      &
                               , step_length &
                               , less_t_r_curr

    !
    ! Gfortran adds an empty line after write(iou,nml=...)
    ! so that it is difficult to put comments like section
    ! headers into the "optimizer.state" file
    !
    select case(act)

    case('restore')

      ! 1) random vars:
      read(iounit,nml=step_module_state)

      ! 2) gdiis_work:
      read(iounit,*) num_stpt, n_coords

      ! need to allocate all, data only for num_stpt
      if( n_coords >= 0 )then
        do i=1,size(gdiis_work)
          if( associated(gdiis_work(i)%q_gdiis) ) CYCLE
          allocate(gdiis_work(i)%q_gdiis(n_coords),gdiis_work(i)%g_gdiis(n_coords))
        enddo
        gdiis_storage_initialized = .true.
      endif

      do i=1,num_stpt ! data only for num_stpt
        read(iounit,*) gdiis_work(i)%num_pt, gdiis_work(i)%ener
        if( n_coords < 0 ) CYCLE

        ! also read data:
        read(iounit,*) gdiis_work(i)%q_gdiis(:)
        read(iounit,*) gdiis_work(i)%g_gdiis(:)
      enddo

      ! 3) q_save, grad_save (these are used to back off to
      !    the previous geometry, if the new is too bad):
      read(iounit,*) n_coords
      if( n_coords >= 0 )then
        if( .not.allocated(q_save) )then
          allocate(q_save(n_coords),grad_save(n_coords))
        endif
        read(iounit,*) q_save(:)
        read(iounit,*) grad_save(:)
      endif
    
    case('save')

      ! 1) random vars:
      write(iounit,nml=step_module_state)

      ! 2) gdiis_work:
      n_coords = -1
      if( associated(gdiis_work(1)%q_gdiis) ) n_coords = size(gdiis_work(1)%q_gdiis)

      write(iounit,*) num_stpt, n_coords

      do i=1,num_stpt ! write only valid data
        write(iounit,*) gdiis_work(i)%num_pt, gdiis_work(i)%ener
        if( n_coords < 0 ) CYCLE

        write(iounit,*) gdiis_work(i)%q_gdiis(:)
        write(iounit,*) gdiis_work(i)%g_gdiis(:)
      enddo

      ! 3) q_save, grad_save:
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
  end subroutine step_module_persistent_state
 
end module step_module
  !--------------- End of module -------------------------------------

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
module  gradient_2c_fit_ch_module
  !---------------------------------------------------------------
  !  Purpose: contains all routines necessary to perform
  !           the two center fitfct gradient integrals:
  !           - [Da Fk|Fl] (Da (mat_charge))
  !  Module called by: integral_calc_quad_2cff.f90
  !
  !  References: ...
  !  Author: FN
  !  Date: 7/96
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: FN
  ! Date:   19.9.96
  ! Description: Discovered a sign error using the solid
  !              harmonics. Most corrections done in lr2_calc.
  !              One was done in ll_calc.
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
!define FPP_TIMERS 2
#include "def.h"
  use type_module ! type specification parameters
  use symmetry_data_module, only : get_totalsymmetric_irrep
  use unique_atom_module
  use solid_harmonics_module, only : solid_harmonics_scalar,&
       solid_harmonics_calc
  use solhrules_module, only : prod_rule, diff_rule,&
       solhrules_differential
  use gamma_module, only: gamma
  use int_data_2cff_module
  use fitcontract_2c_grad_module, only: fitcontract_2c_grad, &
                                        fitcontract_2c_dervs
  use output_module, only: output_int_fitcontract,output_int_2c_fit
  use iounitadmin_module, only: write_to_output_units
  use gradient_data_module, only: grad_fit_ch_cartesian ! for check only
  use calc3c_switches
  USE_MEMLOG
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !------------ public functions and subroutines ------------------
  public charge_overlap_grad


  !================================================================
  ! End of public interface of module
  !================================================================
  ! constants
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind

  real(kind=r8_kind),parameter    :: two=2.0_r8_kind,three=3.0_r8_kind, &
       one=1.0_r8_kind,four=4.0_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  ! variables needed by all subroutines: these are set in the
  ! routine get_exponent_data(na,la,nb,lb)
  integer(kind=i4_kind)           :: naexps, nbexps
  real(kind=r8_kind), allocatable :: fact0(:), fact1(:), fact2(:)
  logical,allocatable             :: cutoff(:, :)
  integer(kind=i4_kind)           :: num, i_ir
  logical                         :: diag_unique
  ! mapping to get the cartesian components of grad yl right
  integer(kind=i4_kind),dimension(3) :: yl_map
  data yl_map /3,4,2/
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  subroutine charge_overlap_grad()
    !  Purpose: calculate grad(mat_charge)
    !** End of interface *****************************************
    implicit none
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)     :: la,lb,na,nb
    !------------ Executable code --------------------------------
 
    select case (quadrupel%l1)
    case (-1)
       la = 0
    case (0)
       la = -1
    case default
       la=quadrupel%l1
    end select
    select case (quadrupel%l2)
    case (-1)
       lb = 0
    case (0)
       lb = -1
    case default
       lb=quadrupel%l2
    end select
    na=quadrupel%ua1
    nb=quadrupel%ua2

    diag_unique = (quadrupel%ua1.eq.quadrupel%ua2)
    
FPP_TIMER_START(t_2c_total)
    if ( (la.gt.0.and.lb.ge.0).or.(la.ge.0.and.lb.gt.0)) then
       call ll_calc_grad(na,la,nb,lb)
    elseif (la.eq.0.and.lb.eq.0) then
       call ss_calc_grad(na,nb)
    elseif( (la.eq.-1.and.lb.ge.0).or.(lb.eq.-1.and.la.ge.0)) then
       call lr2_calc_grad(na,la,nb,lb)
    elseif (la.eq.-1.and.lb.eq.-1) then
       call r2_calc_grad(na,nb)
    else
       call error_handler("CHARGE_OVERLAP_GRAD : something fishy")
    endif
FPP_TIMER_STOP(t_2c_total)

  end subroutine charge_overlap_grad


  subroutine ss_calc_grad(na,nb)
   use cpksdervs_matrices,only: cpksalloc
   use calc3c_switches
 
    !         
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************

    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: &
         intermediate_dervs(:,:,:),intermediate(:,:),gamma_help(:,:), &
         fit_int_grad_a(:,:,:,:,:),fit_int_grad_b(:,:,:,:,:), &
         fit_int_dervs_a(:,:,:,:,:,:,:),fit_int_dervs_b(:,:,:,:,:,:,:)


    integer(kind=i4_kind) :: a_eq,b_eq,alloc_stat,i,begin,k2dr
    integer(kind=i4_kind) :: ma, mb
    logical :: moving_a,moving_b

    

    ma = unique_atoms(na)%moving_atom
    mb = unique_atoms(nb)%moving_atom

    moving_a = ma > 0
    moving_b = mb > 0
    if (.not.moving_a .and. .not.moving_b) return
    call get_exponent_data(na,0,nb,0)

    if(integralpar_2dervs.and..not.integralpar_cpks_contribs) then
     allocate( gamma_help(num,3),intermediate(num,3), &
               intermediate_dervs(num,3,3),  STAT=cpksalloc(25))
    ASSERT(cpksalloc(25).eq.0)
    intermediate_dervs=zero
    MEMLOG(size(gamma_help)+size(intermediate)+size(intermediate_dervs))
    else
     allocate(gamma_help(num,2),intermediate(num,3),STAT=cpksalloc(57))
    ASSERT(cpksalloc(57).eq.0)
    endif


    gamma_help=zero
    intermediate=zero

    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! Loop over ALL equal atoms for symmetry adaption

    ! First calculate grad_a [Fk(a)|Fl(b)] for each equal atom of nb

    mova: if (moving_a) then

       allocate(fit_int_grad_a(naexps,nbexps,1,1,3),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_a=zero

    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
     allocate(fit_int_dervs_a(naexps,nbexps,1,1,3,3,unique_atoms(nb)%n_equal_atoms), &
                stat=cpksalloc(56))
       ASSERT(cpksalloc(56).eq.0)
       MEMLOG(size(fit_int_dervs_a))
       fit_int_dervs_a=zero
    endif

    a_eq=1
    xa=unique_atoms(na)%position(:,a_eq)

    begin=1
    equb: do b_eq=begin,unique_atoms(nb)%n_equal_atoms
       xb=unique_atoms(nb)%position(:,b_eq)
       
       arg=sum((xa-xb)**2)

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
        gamma_help(:,1:3) = gamma(3,fact1/fact0*arg)
       else
        gamma_help(:,1:2) = gamma(2,fact1/fact0*arg)
       endif

       ! contribs from all equalb summed up
       grad_ia: do i=1,3
          intermediate(:,i) = intermediate(:,i) + &
               two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
               two*fact1/fact0*(xa(i)-xb(i))*(-gamma_help(:,2))

      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

         intermediate_dervs(:,i,:)= &
               spread(two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0)) &
              *two*fact1/fact0*(xa(i)-xb(i)) &
              *two*fact1/fact0*gamma_help(:,3),2,3)* &
               spread((xa(:)-xb(:)),1,num)

         intermediate_dervs(:,i,i)=intermediate_dervs(:,i,i) &
        +two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
         two*fact1/fact0*(-gamma_help(:,2))
        
        do k2dr=1,3
         fit_int_dervs_a(:,:,1,1,i,k2dr,b_eq) = &
                 unpack(intermediate_dervs(:,i,k2dr),cutoff,zero)
         enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
      endif

     enddo grad_ia

    enddo equb

    do i=1,3
       fit_int_grad_a(:,:,1,1,i) = unpack(intermediate(:,i),cutoff,zero)
    enddo


       ! do contractions and store result

        call fitcontract_2c_grad(fit_int_grad_a,ma,1) !ss

       deallocate(fit_int_grad_a,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

       
       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then

!       print*,'fit_int_dervs_a ss',na,nb
!       if(size(fit_int_dervs_a,7).gt.1) then
!        do k2dr=1,3
!        print*,sum(fit_int_dervs_a(:,:,1,1,1,k2dr,1)), &
!               sum(fit_int_dervs_a(:,:,1,1,2,k2dr,1)), &
!               sum(fit_int_dervs_a(:,:,1,1,3,k2dr,1)), &
!               sum(fit_int_dervs_a(:,:,1,1,1,k2dr,2)), &
!               sum(fit_int_dervs_a(:,:,1,1,2,k2dr,2)), &
!               sum(fit_int_dervs_a(:,:,1,1,3,k2dr,2))
!        enddo
!       else
!        do k2dr=1,3
!        print*,sum(fit_int_dervs_a(:,:,1,1,k2dr,1,1)), &
!                  sum(fit_int_dervs_a(:,:,1,1,k2dr,2,1)), &
!                     sum(fit_int_dervs_a(:,:,1,1,k2dr,3,1))
!        enddo
!       endif

FPP_TIMER_START(t_contract_2c_dervs)
         call fitcontract_2c_dervs(fit_int_dervs_a,ma,mb)  !ss(1)
FPP_TIMER_STOP(t_contract_2c_dervs)

        MEMLOG(-size(fit_int_dervs_a))
        deallocate(fit_int_dervs_a,stat=cpksalloc(56))
        ASSERT(cpksalloc(56).eq.0)
cpksalloc(56)=1

       endif

    endif mova

    !if (na.ne.nb) then

       intermediate = zero
       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) &
                                       intermediate_dervs = zero

       ! Now calculate grad_b [Fk(a)|Fl(b)] for each equal atom of na

       if (moving_b) then
          allocate(fit_int_grad_b(naexps,nbexps,1,1,3),STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          fit_int_grad_b=zero
          
          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           allocate(fit_int_dervs_b(naexps,nbexps,1,1,3,3,unique_atoms(na)%n_equal_atoms), &
                    STAT=cpksalloc(66))
          ASSERT(cpksalloc(66).eq.0)
          MEMLOG(size(fit_int_dervs_b))
          fit_int_dervs_b=zero
          endif

       begin=1
       b_eq=1
       xb = unique_atoms(nb)%position(:,b_eq)
       equa: do a_eq=begin,unique_atoms(na)%n_equal_atoms
          xa=unique_atoms(na)%position(:,a_eq)
          
          arg=sum((xa-xb)**2)

         if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          gamma_help(:,1:3) = gamma(3,fact1/fact0*arg)
         else
          gamma_help(:,1:2) = gamma(2,fact1/fact0*arg)
         endif

          grad_ib: do i=1,3
             intermediate(:,i) = intermediate(:,i) + &
                  two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
                  two*fact1/fact0*(xa(i)-xb(i))*&
                  gamma_help(:,2)

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

             intermediate_dervs(:,i,:) = & !!! intermediate_dervs(:,i,:) - &
              -spread( two*(xa-xb),1,num)* &
               spread( two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
                  two*fact1/fact0*(xa(i)-xb(i))*&
                  (-gamma_help(:,3))*fact1/fact0,2,3)

             intermediate_dervs(:,i,i) = intermediate_dervs(:,i,i) &
                 -two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
                  two*fact1/fact0* gamma_help(:,2)

           do k2dr=1,3
            fit_int_dervs_b(:,:,1,1,i,k2dr,a_eq) =  &
                   unpack(intermediate_dervs(:,i,k2dr),cutoff,zero)
           enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
           endif
        
          enddo grad_ib
       enddo equa

       do i=1,3
          fit_int_grad_b(:,:,1,1,i) = unpack(intermediate(:,i),cutoff,zero)
       enddo

          ! do contractions and store result to grad_fit_ch_cartesian

!    print*,'fit_int_grad_b ss',na,nb
!    print*,sum(fit_int_grad_b(:,:,1,1,1)), &
!              sum(fit_int_grad_b(:,:,1,1,2)), &
!                 sum(fit_int_grad_b(:,:,1,1,3))

           call fitcontract_2c_grad(fit_int_grad_b,mb,2)   ! ss

          deallocate(fit_int_grad_b,STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
         
          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then


!       print*,'fit_int_dervs_b ss',mb,ma
!       if(size(fit_int_dervs_b,7).gt.1) then
!        do k2dr=1,3
!        print*,sum(fit_int_dervs_b(:,:,1,1,1,k2dr,1)), &
!               sum(fit_int_dervs_b(:,:,1,1,2,k2dr,1)), &
!               sum(fit_int_dervs_b(:,:,1,1,3,k2dr,1)), &
!               sum(fit_int_dervs_b(:,:,1,1,1,k2dr,2)), &
!               sum(fit_int_dervs_b(:,:,1,1,2,k2dr,2)), &
!               sum(fit_int_dervs_b(:,:,1,1,3,k2dr,2))
!        enddo
!       else
!        do k2dr=1,3
!        print*,sum(fit_int_dervs_b(:,:,1,1,k2dr,1,1)), &
!                  sum(fit_int_dervs_b(:,:,1,1,k2dr,2,1)), &
!                     sum(fit_int_dervs_b(:,:,1,1,k2dr,3,1))
!        enddo
!       endif

FPP_TIMER_START(t_contract_2c_dervs)
            call fitcontract_2c_dervs(fit_int_dervs_b,mb,ma) ! ss (2)
FPP_TIMER_STOP(t_contract_2c_dervs)
           MEMLOG(-size(fit_int_dervs_b))
           deallocate(fit_int_dervs_b,STAT=cpksalloc(66))
          ASSERT(cpksalloc(66).eq.0)
          cpksalloc(66)=1
          endif

       endif

    
    ! deallocation
    if(integralpar_2dervs.and..not.integralpar_cpks_contribs) then
    MEMLOG(-size(intermediate)-size(gamma_help)-size(intermediate_dervs))
     deallocate(intermediate,gamma_help,intermediate_dervs,STAT=cpksalloc(25))
     ASSERT(cpksalloc(25).eq.0)
     cpksalloc(25)=1
    else
     deallocate(intermediate,gamma_help,STAT=cpksalloc(57))
     ASSERT(cpksalloc(57).eq.0)
     cpksalloc(57)=1
    endif

    call dump_exponent_data()

  end subroutine ss_calc_grad


  subroutine ll_calc_grad(na,la,nb,lb)
    ! Purpose: calculate primitives for both angular momenta
    !          greater than zero.
    !          Symmetry adaption included.
    !------------ Modules used ----------------------------------
    !------------ Declaration of formal parameters ---------------

use calc3c_switches
use cpksdervs_matrices
    integer(kind=i4_kind),intent(inout)  :: la,lb
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help arrays for gamma-function
    real(kind=r8_kind),allocatable,dimension(:,:)   :: gamma_arg
    real(kind=r8_kind),allocatable,dimension(:,:)   :: gamma_help
    real(kind=r8_kind),allocatable,dimension(:,:,:) :: yl_arr
    real(kind=r8_kind),allocatable,dimension(:,:,:) :: yl_dervs
    real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: yl_grad

    real(kind=r8_kind),allocatable,dimension(:)     :: yl_arr_grad
    ! intermediate arrays
    real(kind=r8_kind),allocatable,dimension(:,:)   :: intermediate2,&
         inter1,inter2,term_1,term_2
    real(kind=r8_kind),allocatable,dimension(:,:) :: intermediate1
    real(kind=r8_kind),allocatable,dimension(:,:,:) :: inter1_dervs,term_1_dervs
    real(kind=r8_kind),allocatable,dimension(:,:,:) :: inter2_dervs,term_2_dervs
    real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: intermediate1_dervs
    real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: intermediate2_dervs
    ! uncontracted but symmetryadapted integrals:
    ! fit_int(num,i_ind1,i_ind2)
    real(kind=r8_kind),allocatable    :: &
         fit_int_grad_a(:,:,:,:,:),fit_int_dervs_a(:,:,:,:,:,:,:), &
         fit_int_grad_b(:,:,:,:,:),fit_int_dervs_b(:,:,:,:,:,:,:)
    ! number of independents and contributing fcts
    integer(kind=i4_kind)   :: n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind),allocatable  :: n_contributing_a(:), &
         n_contributing_b(:)
    real(kind=r8_kind),pointer       :: coeff_a(:),coeff_b(:)
    integer(kind=i4_kind),pointer,dimension(:)   :: magn_a,&
         magn_b,eq_atom_a,eq_atom_b
    real(kind=r8_kind),allocatable :: args(:,:)
    integer(kind=i4_kind)   :: counter,i_l,i_m, &
         alloc_stat,i_ind1,i_ind2,i_cont1,i_cont2,&
         l_max,eq_b,eq_a,l_prod,i,lm,i_sum

    integer(kind=i4_kind) :: ma, mb, begin, k2dr
    logical :: moving_a,moving_b
    intrinsic max
    ! ---------- Executable code --------------------------------

    ma = unique_atoms(na)%moving_atom
    mb = unique_atoms(nb)%moving_atom
    moving_a = ma > 0
    moving_b = mb > 0
    if (.not.moving_a .and. .not.moving_b) return
    call get_exponent_data(na,la,nb,lb)

    n_indep_a=unique_atoms(na)%symadapt_partner(i_ir,la)%n_independent_fcts
    n_indep_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%n_independent_fcts
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    l_max=max(la,lb)

    allocate( &
         n_contributing_a(n_indep_a), &
         n_contributing_b(n_indep_b), &
         yl_arr(n_equal_a,n_equal_b,(l_max+1)**2), &
         yl_grad(n_equal_a,n_equal_b,(l_max+1)**2,3), &
         gamma_arg(num,3), STAT=alloc_stat)
         ASSERT(alloc_stat.eq.0)

    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
     allocate(gamma_help(num,la+lb+3),STAT=alloc_stat)
    else
     allocate(gamma_help(num,la+lb+2),STAT=alloc_stat)
    endif
    ASSERT(alloc_stat.eq.0)

    yl_arr = zero
    yl_grad = zero
    gamma_help=zero
    gamma_arg=zero

    if (n_indep_a /= 0) then
       n_contributing_a=unique_atoms(na)%symadapt_partner(i_ir,la)%&
            &symadapt(:,1)%n_fcts
    endif

    if (n_indep_b/=0) then
       n_contributing_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
            &symadapt(:,1)%n_fcts
    endif

    ! To avoid doublication of code double DO-loops over all pairs of equal
    ! atoms are used, though only the following pairs are final processed
    !   a) if moving_a : the pairs (1,eq_b) and
    !   b) if moving_b : the pairs (eq_a,1) .

    ! do a precalculation of solid harmonics 
    ! gradients are calculated during the actual evaluation
    ! of primitive matrix elements ------------------------------------

    begin = 1
    if (moving_a) then
    allocate(args(n_equal_b,3),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    eq_a = 1
       args =  transpose (spread(unique_atoms(na)%position(:,eq_a),2,n_equal_b) - &
            unique_atoms(nb)%position)

       yl_arr(eq_a,:,:) = solid_harmonics_calc(l_max,args)

                      do eq_b=1,n_equal_b

                       do i=1,3
                       do lm=1,(l_max+1)**2
                        do i_sum=1,solhrules_differential(yl_map(i),lm)%n_summands
                         yl_grad(eq_a,eq_b,lm,i) = &
                              yl_grad(eq_a,eq_b,lm,i) + &
                              solhrules_differential(yl_map(i),lm)%coef(i_sum)* &
                              yl_arr(eq_a,eq_b,&
                              solhrules_differential(yl_map(i),lm)%&
                              lm_sh(i_sum))
                          enddo 
                         enddo
                       enddo 
                    enddo

    deallocate( args,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
      begin = 2
    endif

    if (moving_b) then
       allocate(args(1,3),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
          eq_b = 1
       do eq_a = begin,n_equal_a
          args(1,:) = unique_atoms(na)%position(:,eq_a) - &
                      unique_atoms(nb)%position(:,1)

          yl_arr(eq_a,1:1,:) = solid_harmonics_calc(l_max,args)

                       do i=1,3
                       do lm=1,(l_max+1)**2
                        do i_sum=1,solhrules_differential(yl_map(i),lm)%n_summands
                         yl_grad(eq_a,eq_b,lm,i) = &
                              yl_grad(eq_a,eq_b,lm,i) + &
                              solhrules_differential(yl_map(i),lm)%coef(i_sum)* &
                              yl_arr(eq_a,eq_b,&
                              solhrules_differential(yl_map(i),lm)%&
                              lm_sh(i_sum))
                          enddo 
                         enddo
                       enddo 

       enddo
       deallocate( args,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    endif

    ! ------------------------------------------------------------------
    if (moving_a) then
       allocate( fit_int_grad_a(naexps,nbexps,n_indep_a,n_indep_b,3), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_a=zero

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
       allocate( fit_int_dervs_a(naexps,nbexps,n_indep_a,n_indep_b,3,3,n_equal_b), &
                 STAT=cpksalloc(74))
       ASSERT(cpksalloc(74).eq.0)
       MEMLOG(size(fit_int_dervs_a))
       fit_int_dervs_a=zero
     endif

    endif

    if (moving_b) then
       allocate( fit_int_grad_b(naexps,nbexps,n_indep_a,n_indep_b,3), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_b=zero

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
       allocate( fit_int_dervs_b(naexps,nbexps,n_indep_a,n_indep_b,3,3,n_equal_a), &
                 STAT=cpksalloc(75))
       ASSERT(cpksalloc(75).eq.0)
       MEMLOG(size(fit_int_dervs_b))
       fit_int_dervs_b=zero
     endif
    endif

    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b

    if (la.gt.0.and.lb.eq.0) then ! la > 0 and lb = 0 ----------------
       !                           => ( -2 ) **la

       if (moving_a) then
          allocate( intermediate1(num,1), STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          allocate( intermediate1_dervs(num,1,3,3), STAT=cpksalloc(76))
          ASSERT(cpksalloc(76).eq.0)
          MEMLOG(size(intermediate1_dervs))
          intermediate1_dervs=zero
        endif
       endif

       if (moving_b) then
          allocate( intermediate2(num,1), STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)

        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          allocate( intermediate2_dervs(num,1,3,3), STAT=cpksalloc(78))
          ASSERT(cpksalloc(78).eq.0)
          MEMLOG(size(intermediate2_dervs))
          intermediate2_dervs=zero
        endif

       endif

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           allocate(yl_dervs((la+1)**2,3,3),stat=cpksalloc(77))
           MEMLOG(size(yl_dervs))
           ASSERT(cpksalloc(77).eq.0)
           yl_dervs=zero
          endif

    cartesian_ls :do i=1,3

          do i_ind1=1,n_indep_a
             if (moving_a)  then
               intermediate1 = zero
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) intermediate1_dervs=zero
             endif
             if (moving_b) then
               intermediate2 = zero
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) intermediate2_dervs=zero
             endif
                eq_atom_a => &
                     unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%I_equal_atom
                coeff_a => &
                     unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%c  
                magn_a => &
                     unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%m

             cont1 : do i_cont1=1,n_contributing_a(i_ind1)

                eq_a=eq_atom_a(i_cont1)
                equal_b: do eq_b=1,n_equal_b

                   if ( (eq_atom_a(i_cont1) /= 1 .or. .not.moving_a) .and. &
                        (eq_b /= 1 .or. .not.moving_b) ) cycle

                   lm = la**2 + magn_a(i_cont1)

                   xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
                   xb = unique_atoms(nb)%position(:,eq_b)

                   arg = sum((xa-xb)**2)

                   if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                   gamma_help(:,1:la+3) = gamma(la+3,fact1/fact0*arg)
                   ! now yl_dervs calculated
                   yl_dervs(lm,i,:)=zero
                   if (.not.((na.eq.nb).and.(eq_b.eq.eq_atom_a(i_cont1)))) then
                      do k2dr=1,3
                      do i_sum=1,solhrules_differential(yl_map(k2dr),lm)%n_summands
                         yl_dervs(lm,i,k2dr) =  yl_dervs(lm,i,k2dr) + &
                              solhrules_differential(yl_map(k2dr),lm)%coef(i_sum)* &
                              yl_grad(eq_atom_a(i_cont1),eq_b,&
                              solhrules_differential(yl_map(k2dr),lm)%lm_sh(i_sum),i)
                      enddo
                      enddo
                   endif
FPP_TIMER_STOP(t_calc_2c_dervs)
                   else
                    gamma_help(:,1:la+2) = gamma(la+2,fact1/fact0*arg)
                   endif
                   
                   ! First evaluate grad_a [ Fk(a) | Fl(b)]
                   if (moving_a .and. eq_atom_a(i_cont1) == 1) then
                      intermediate1(:num,1) = intermediate1(:,1)  &
                          -yl_arr(1,eq_b,lm) * &
                           gamma_help(:,la+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_a(i_cont1)
                      
                      intermediate1(:,1) = intermediate1(:,1) + &
                           gamma_help(:,la+1)*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_a(i_cont1)

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      intermediate1_dervs(:num,1,i,:3) =  &          !!!d/dr yl_arr
                       -spread( yl_grad(1,eq_b,lm,:3),1,num)*&
                        spread( gamma_help(:,la+2)*(xa(i)-xb(i))* &
                                two*fact1/fact0*coeff_a(i_cont1),2,3)

                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:) - &
                        spread( two*(xa-xb),1,num)*&
                        spread( yl_arr(1,eq_b,lm) * &
                           (-gamma_help(:,la+3))*(xa(i)-xb(i))* &     !!! d/dr gamma_help
                           two*fact1/fact0*coeff_a(i_cont1)*fact1/fact0,2,3)

                      intermediate1_dervs(:,1,i,i) = intermediate1_dervs(:,1,i,i) - &
                           yl_arr(1,eq_b,lm) * &
                           gamma_help(:,la+2)*two*fact1/fact0*coeff_a(i_cont1) !!! d/dr r

                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:)  &
                        +spread( two*(xa-xb),1,num)*&                !!! d/dr gamma_help
                         spread( (-gamma_help(:,la+2))*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_a(i_cont1)*fact1/fact0,2,3)

                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:)  &
                       +spread( yl_dervs(lm,i,:),1,num)* &
                        spread( gamma_help(:,la+1)*coeff_a(i_cont1),2,3)

                      do k2dr=1,3
                       fit_int_dervs_a(:,:,i_ind1,1,i,k2dr,eq_b) = &
                       fit_int_dervs_a(:,:,i_ind1,1,i,k2dr,eq_b) + &
                         unpack((-two*fact1/fact0)**la/(fact1*sqrt(fact0)) * &
                         intermediate1_dervs(:,1,i,k2dr),cutoff,zero)* &
                         two*pi*pi*sqrt(pi)
                      enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif
                   endif

                   ! Now calculate grad_b [Fk(a) | Fl(b) ] - since the fitfcts
                   ! are already symmetry-adapted, this is not necessarily
                   ! = - grad_a [Fk(a) | Fl(b) ]

                   if (moving_b .and. eq_b == 1) then
                      intermediate2(:,1) = intermediate2(:,1) + &
                           yl_arr(eq_atom_a(i_cont1),eq_b,lm) * &
                           gamma_help(:,la+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_a(i_cont1)
                      
                      intermediate2(:,1) = intermediate2(:,1) - &
                           gamma_help(:,la+1)*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_a(i_cont1)

                     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      intermediate2_dervs(:,1,i,:) =  &
                          -spread( yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num) * &
                           spread( gamma_help(:,la+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_a(i_cont1),2,3)

                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:)  &
                       -spread( two*(xa-xb),1,num)*&
                        spread( yl_arr(eq_atom_a(i_cont1),eq_b,lm) * &
                           (-gamma_help(:,la+3))*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_a(i_cont1)*fact1/fact0,2,3)

                      intermediate2_dervs(:,1,i,i) = intermediate2_dervs(:,1,i,i)  &
                          -yl_arr(eq_atom_a(i_cont1),eq_b,lm) * &
                           gamma_help(:,la+2)*two*fact1/fact0*coeff_a(i_cont1)

                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:)  &
                       +spread( two*(xa-xb),1,num)*&
                        spread( (-gamma_help(:,la+2))*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_a(i_cont1)*fact1/fact0,2,3)

                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:)  &
                       +spread( yl_dervs(lm,i,:),1,num)*&
                        spread( gamma_help(:,la+1)*coeff_a(i_cont1),2,3)

                      do k2dr=1,3
                       fit_int_dervs_b(:,:,i_ind1,1,i,k2dr,eq_a) =  &
                       fit_int_dervs_b(:,:,i_ind1,1,i,k2dr,eq_a)+ &
                            unpack((-two*fact1/fact0)**la/(fact1*sqrt(fact0)) * &
                            intermediate2_dervs(:,1,i,k2dr),cutoff,zero)* &
                            two*pi*pi*sqrt(pi)
                      enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                     endif
                   endif

                enddo equal_b
             enddo cont1

           map_ls_contrib: if(.true.) then
             if (moving_a) then
             fit_int_grad_a(:,:,i_ind1,1,i) = &
                  unpack((-two*fact1/fact0)**la/(fact1*sqrt(fact0)) * &
                  intermediate1(:,1),cutoff,zero)* &
                  two*pi*pi*sqrt(pi)
             endif
!             !if(na.ne.nb) then

                if (moving_b) then
                fit_int_grad_b(:,:,i_ind1,1,i) = &
                     unpack((-two*fact1/fact0)**la/(fact1*sqrt(fact0)) * &
                     intermediate2(:,1),cutoff,zero)* &
                     two*pi*pi*sqrt(pi)
                endif

             !endif
            endif map_ls_contrib
             
          enddo

       enddo cartesian_ls
          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then                           !!!!!!!
           MEMLOG(-size(yl_dervs))
           deallocate(yl_dervs,stat=cpksalloc(77))
           ASSERT(cpksalloc(77).eq.0)
           cpksalloc(77)=1
          endif

       if (moving_a) then
          deallocate( intermediate1, STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          MEMLOG(-size(intermediate1_dervs))
          deallocate( intermediate1_dervs, STAT=cpksalloc(76))
          ASSERT(cpksalloc(76).eq.0)
          cpksalloc(76)=1
        endif
       endif

       if (moving_b) then
          deallocate( intermediate2, STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          MEMLOG(-size(intermediate2_dervs))
          deallocate( intermediate2_dervs, STAT=cpksalloc(78))
          ASSERT(cpksalloc(78).eq.0)
          cpksalloc(78)=1
        endif
       endif

!    print*,'ls fit_int_grad_a fit_int_dervs',na,nb,xa
!    print*,sum(fit_int_grad_a(:,:,1,1,1)),sum(fit_int_dervs_a(:,:,1,1,1,1,1))


    elseif (la.eq.0.and.lb.gt.0) then ! la > 0 and lb = 0 ----------------
       !                           => ( -2 ) **la

       if (moving_a) then
          allocate( intermediate1(num,1), STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          allocate( intermediate1_dervs(num,1,3,3), STAT=cpksalloc(76))
          ASSERT(cpksalloc(76).eq.0)
          MEMLOG(size(intermediate1_dervs))
        endif
       endif
       if (moving_b) then
          allocate( intermediate2(num,1), STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)

        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          allocate( intermediate2_dervs(num,1,3,3), STAT=cpksalloc(78))
          ASSERT(cpksalloc(78).eq.0)
          MEMLOG(size(intermediate2_dervs))
        endif
       endif

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           allocate(yl_dervs((lb+1)**2,3,3),stat=cpksalloc(77))
           ASSERT(cpksalloc(77).eq.0)
           MEMLOG(size(yl_dervs))
           yl_dervs=zero
          endif

       cartesian_sl :do i=1,3
          do i_ind2=1,n_indep_b
             if (moving_a) then
              intermediate1 = zero
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) intermediate1_dervs=zero
             endif
             if (moving_b) then
              intermediate2 = zero
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) intermediate2_dervs=zero
             endif

                eq_atom_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                     &symadapt(i_ind2,i_ir)%I_equal_atom
                coeff_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                     &symadapt(i_ind2,i_ir)%c  
                magn_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                     &symadapt(i_ind2,i_ir)%m
             cont1_sl : do i_cont2=1,n_contributing_b(i_ind2)
                eq_b=eq_atom_b(i_cont2)
                equal_a: do eq_a=1,n_equal_a

                   if ( (eq_a /= 1 .or. .not.moving_a) .and. & 
                        (eq_atom_b(i_cont2) /= 1 .or. .not.moving_b) ) cycle

                   ! calculate gradient of solid harmonics
                   lm = lb**2 + magn_b(i_cont2)

                   xa = unique_atoms(na)%position(:,eq_a)
                   xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont2))

                   arg = sum((xa-xb)**2)
                  if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                   gamma_help(:,1:lb+3) = gamma(lb+3,fact1/fact0*arg)
                   ! now yl_dervs calculated
                   yl_dervs(lm,i,:) = zero
                   if (.not.((na.eq.nb).and.(eq_a.eq.eq_atom_b(i_cont2)))) then
                      do k2dr=1,3
                      do i_sum=1,solhrules_differential(yl_map(k2dr),lm)%n_summands
                         yl_dervs(lm,i,k2dr) = &
                              yl_dervs(lm,i,k2dr) + &
                              solhrules_differential(yl_map(k2dr),lm)%coef(i_sum)* &
                              yl_grad(eq_a,eq_atom_b(i_cont2),&
                              solhrules_differential(yl_map(k2dr),lm)%lm_sh(i_sum),i)
                      enddo
                      enddo
                   endif
FPP_TIMER_STOP(t_calc_2c_dervs)
                  else
                   gamma_help(:,1:lb+2) = gamma(lb+2,fact1/fact0*arg)
                  endif
                   
                   ! First evaluate grad_a [ Fk(a) | Fl(b)]
                   if (moving_a .and. eq_a == 1) then

                      intermediate1(:,1) = intermediate1(:,1) - &
                           yl_arr(1,eq_atom_b(i_cont2),lm) * &
                           gamma_help(:,lb+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2)
                      
                      intermediate1(:,1) = intermediate1(:,1) + &
                           gamma_help(:,lb+1)*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_b(i_cont2)

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      intermediate1_dervs(:,1,i,:) =  &
                          -spread(yl_grad(1,eq_atom_b(i_cont2),lm,:),1,num) * &
                           spread(gamma_help(:,lb+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2),2,3)

                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:) - &
                           spread(two*(xa-xb),1,num)*&
                           spread(yl_arr(1,eq_atom_b(i_cont2),lm) * &
                           (-gamma_help(:,lb+3))*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2)*fact1/fact0,2,3)

                      intermediate1_dervs(:,1,i,i) = intermediate1_dervs(:,1,i,i) - &
                           yl_arr(1,eq_atom_b(i_cont2),lm) * &
                           gamma_help(:,lb+2)*two*fact1/fact0*coeff_b(i_cont2)
                      
                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:) + &
                           spread(two*(xa-xb),1,num)*& 
                           spread((-gamma_help(:,lb+2))*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_b(i_cont2)*fact1/fact0,2,3)

                      intermediate1_dervs(:,1,i,:) = intermediate1_dervs(:,1,i,:) + &
                           spread(yl_dervs(lm,i,:),1,num)*&
                           spread(gamma_help(:,lb+1)*coeff_b(i_cont2),2,3)

                      do k2dr=1,3
                      fit_int_dervs_a(:,:,1,i_ind2,i,k2dr,eq_b) = &
                      fit_int_dervs_a(:,:,1,i_ind2,i,k2dr,eq_b) + &
                         unpack((two*fact1/fact0)**lb/(fact1*sqrt(fact0)) * &
                         intermediate1_dervs(:,1,i,k2dr),cutoff,zero)* &
                         two*pi*pi*sqrt(pi)
                      enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif
                   endif
                   ! Now calculate grad_b [Fk(a) | Fl(b) ] - since the fitfcts
                   ! are already symmetry-adapted, this is not necessarily
                   ! = - grad_a [Fk(a) | Fl(b) ]
                   if (moving_b .and. eq_atom_b(i_cont2) == 1) then

                      intermediate2(:,1) = intermediate2(:,1) + &
                           yl_arr(eq_a,eq_atom_b(i_cont2),lm) * &
                           gamma_help(:,lb+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2)
                      
                      intermediate2(:,1) = intermediate2(:,1) - &
                           gamma_help(:,lb+1)*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_b(i_cont2)

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                      intermediate2_dervs(:,1,i,:) = &
                          -spread(yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num) * &
                           spread(gamma_help(:,lb+2)*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2),2,3)

                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:)  &
                          -spread(two*(xa-xb),1,num)*&
                           spread(yl_arr(eq_a,eq_atom_b(i_cont2),lm) * &
                           (-gamma_help(:,lb+3))*(xa(i)-xb(i))* &
                           two*fact1/fact0*coeff_b(i_cont2)*fact1/fact0,2,3)

                      intermediate2_dervs(:,1,i,i) = intermediate2_dervs(:,1,i,i)  &
                          -yl_arr(eq_a,eq_atom_b(i_cont2),lm) * &
                           gamma_help(:,lb+2)*two*fact1/fact0*coeff_b(i_cont2)
                      
                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:) + &
                           spread(two*(xa-xb),1,num)*&
                           spread((-gamma_help(:,lb+2))*yl_grad(eq_a,eq_b,lm,i) * &
                           coeff_b(i_cont2)*fact1/fact0,2,3)

                      intermediate2_dervs(:,1,i,:) = intermediate2_dervs(:,1,i,:) + &
                           spread(yl_dervs(lm,i,:),1,num)*&
                           spread(gamma_help(:,lb+1)*coeff_b(i_cont2),2,3)

                      do k2dr=1,3
                       fit_int_dervs_b(:,:,1,i_ind2,i,k2dr,eq_a) = &
                       fit_int_dervs_b(:,:,1,i_ind2,i,k2dr,eq_a) + &
                            unpack((two*fact1/fact0)**lb/(fact1*sqrt(fact0)) * &
                            intermediate2_dervs(:,1,i,k2dr),cutoff,zero)* &
                            two*pi*pi*sqrt(pi)
                      enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif
                   endif
                      

                enddo equal_a
             enddo cont1_sl

         map_sl_contrib: if(.true.) then
             if (moving_a) then
             fit_int_grad_a(:,:,1,i_ind2,i) = &
                  unpack((two*fact1/fact0)**lb/(fact1*sqrt(fact0)) * &
                  intermediate1(:,1),cutoff,zero)* &
                  two*pi*pi*sqrt(pi)
             endif
             !if(na.ne.nb) then
                if (moving_b) then
                fit_int_grad_b(:,:,1,i_ind2,i) = &
                     unpack((two*fact1/fact0)**lb/(fact1*sqrt(fact0)) * &
                     intermediate2(:,1),cutoff,zero)* &
                     two*pi*pi*sqrt(pi)
                endif
             !endif
          endif map_sl_contrib
             
          enddo

       enddo cartesian_sl

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           MEMLOG(-size(yl_dervs))
           deallocate(yl_dervs,stat=cpksalloc(77))
           ASSERT(cpksalloc(77).eq.0)
           cpksalloc(77)=1
          endif

       if (moving_a) then
          deallocate( intermediate1, STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          MEMLOG(-size(intermediate1_dervs))
          deallocate( intermediate1_dervs, STAT=cpksalloc(76))
          ASSERT(cpksalloc(76).eq.0)
          cpksalloc(76)=1
        endif
       endif
       if (moving_b) then
          deallocate( intermediate2, STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
        if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
          MEMLOG(-size(intermediate2_dervs))
          deallocate( intermediate2_dervs, STAT=cpksalloc(78))
          ASSERT(cpksalloc(78).eq.0)
          cpksalloc(78)=1
        endif
       endif

!    print*,'sl fit_int_grad_a fit_int_dervs',na,nb,la,lb
!    print*,sum(fit_int_grad_a(:,:,1,1,1)),sum(fit_int_dervs_a(:,:,1,1,1,1,1))

!    print*,'sl fit_int_grad_b fit_int_dervs',na,nb,la,lb
!    print*,sum(fit_int_grad_b(:,:,1,1,1)),sum(fit_int_dervs_b(:,:,1,1,1,1,1))

    elseif ( la.gt.0.and.lb.gt.0) then

       allocate(  intermediate1(num,1:(lb+1)**2), &
            intermediate2(num,1:(la+1)**2), &
            inter1(num,1:(la+1)**2), &
            inter2(num,1:(la+1)**2), &
            term_1(num,1), &
            term_2(num,1), &
            yl_arr_grad((l_max+1)**2), &
            STAT=alloc_stat)
            ASSERT(alloc_stat.eq.0)
       intermediate1 = zero
       intermediate2 = zero
       inter1 = zero
       inter2 = zero
       term_1 = zero
       term_2 = zero

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
        allocate(intermediate1_dervs(num,1:(lb+1)**2,3,3), &
                 intermediate2_dervs(num,1:(la+1)**2,3,3), &
                 yl_dervs((l_max+1)**2,3,3), &
                 inter1_dervs(num,1:(la+1)**2,3), &
                 inter2_dervs(num,1:(la+1)**2,3), &
                 term_1_dervs(num,1,3),& 
                 term_2_dervs(num,1,3),& 
                 STAT=cpksalloc(79))
        MEMLOG(size(intermediate1_dervs)+size(intermediate2_dervs))
        MEMLOG(size(yl_dervs)+size(inter1_dervs)+size(inter2_dervs))
        MEMLOG(size(term_1_dervs)+size(term_2_dervs))
        ASSERT(cpksalloc(79).eq.0)
        intermediate1_dervs=zero
        intermediate2_dervs=zero
        inter1_dervs=zero
        inter2_dervs=zero
        term_1_dervs = zero
        term_2_dervs = zero
       endif

       yl_arr_grad = zero

       cartesian_ll: do i=1,3

          independent1: do i_ind1=1,n_indep_a
                coeff_a =>  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%c
                eq_atom_a =>unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%I_equal_atom
                magn_a =>unique_atoms(na)%symadapt_partner(i_ir,la)%&
                     &symadapt(i_ind1,i_ir)%m
                independent2: do i_ind2=1,n_indep_b
                      coeff_b => unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                           &symadapt(i_ind2,i_ir)%c
                      eq_atom_b => unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                           &symadapt(i_ind2,i_ir)%I_equal_atom
                      magn_b =>  unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                           &symadapt(i_ind2,i_ir)%m

             contributing1: do i_cont1=1,n_contributing_a(i_ind1)
                   contributing2: do i_cont2=1,n_contributing_b(i_ind2)

                      eq_a = eq_atom_a(i_cont1)
                      eq_b = eq_atom_b(i_cont2)

                      if ( (eq_a /= 1 .or. .not.moving_a) .and. &
                           (eq_b /= 1 .or. .not.moving_b) ) cycle

                   
                      xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
                      xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont2))
                      arg = sum((xa-xb)**2)

                     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                      gamma_help(:,1:la+lb+3) = gamma(la+lb+3,fact1/fact0*arg)
                      yl_dervs(:,i,:) = zero
                      if (.not.((na.eq.nb).and.(eq_a.eq.eq_b))) then
                         lm = (l_max+1)**2
                         do counter = 1,lm
                         do k2dr=1,3
                            do i_sum=1,solhrules_differential(yl_map(k2dr),counter)%n_summands
                               yl_dervs(counter,i,k2dr) =  yl_dervs(counter,i,k2dr) + &
                                    solhrules_differential(yl_map(k2dr),counter)%coef(i_sum)* &
                                    yl_grad(eq_a,eq_b,&
                                    solhrules_differential(yl_map(k2dr),counter)%lm_sh(i_sum),i)
                            enddo
                         enddo
                         enddo
                      endif
FPP_TIMER_STOP(t_calc_2c_dervs)
                     else
                      gamma_help(:,1:la+lb+2) = gamma(la+lb+2,fact1/fact0*arg)
                     endif



                      ! prepare the array containing the 'yl_arr' entering the diff_rule
                      do counter =1,(lb+1)**2
                         intermediate1(:,counter) = spread(yl_arr(eq_a,eq_b,counter),1,num)

                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                         intermediate1_dervs(:,counter,i,:)= &
                               spread(yl_grad(eq_a,eq_b,counter,:),1,num)
FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif

                      enddo

                      inter1(:,1:(la+1)**2)= &
                           diff_rule(intermediate1,1,(la+1)**2,lb**2+magn_b(i_cont2))
                      inter1 = inter1 * coeff_b(i_cont2)

                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                        do k2dr=1,3
                         inter1_dervs(:,1:(la+1)**2,k2dr)= &
                           diff_rule(intermediate1_dervs(:,:,i,k2dr),1,(la+1)**2,lb**2+magn_b(i_cont2))
                         inter1_dervs(:,:,k2dr) = inter1_dervs(:,:,k2dr) * coeff_b(i_cont2)
                        enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif
                      
                      ! prepare the array containing the 'yl_arr_grad' entering the diff_rule
                      intermediate1 = zero
                      do counter=1,(lb+1)**2
                         intermediate1(:,counter) = spread(yl_grad(eq_a,eq_b,counter,i),1,num)

                         if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                          intermediate1_dervs(:,counter,i,:)= spread(yl_dervs(counter,i,:),1,num)
FPP_TIMER_STOP(t_calc_2c_dervs)
                         endif

                      enddo

                      inter2(:,1:(la+1)**2)= &
                        diff_rule(intermediate1,1,(la+1)**2,lb**2+magn_b(i_cont2))
                      inter2 = inter2 * coeff_b(i_cont2)

                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                       do k2dr=1,3
                       inter2_dervs(:,1:(la+1)**2,k2dr)=diff_rule &
                           (intermediate1_dervs(:,:,i,k2dr),1,(la+1)**2,lb**2+magn_b(i_cont2))
                       enddo
                       inter2_dervs=inter2_dervs* coeff_b(i_cont2) !!! lost stat
FPP_TIMER_STOP(t_calc_2c_dervs)
                      endif
                      
                      l_prod=la**2+magn_a(i_cont1)
                      
                      ! now produce term_1 = prod_rule(inter2,(yl_arr(la-i_l)*gamma(la+lb-i_l)))
                      counter=1
                      do i_l=0,la
                         do i_m=1,2*i_l+1
                            intermediate2(:,counter) =  &
                                 yl_arr(eq_a,eq_b,counter) * &
                                 gamma_help(:,lb+i_l+1) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l

                            if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                            intermediate2_dervs(:,counter,i,:) =  &
                                 spread(yl_grad(eq_a,eq_b,counter,:),1,num) * &
                                 spread(gamma_help(:,lb+i_l+1) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l,2,3)

                            intermediate2_dervs(:,counter,i,:) =  &
                                              intermediate2_dervs(:,counter,i,:)+  &
                                spread(two*(xa-xb),1,num)*&
                                spread(yl_arr(eq_a,eq_b,counter)*(-gamma_help(:,lb+i_l+2)) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l*fact1/fact0,2,3)
FPP_TIMER_STOP(t_calc_2c_dervs)
                            endif

                            counter=counter+1
                         enddo
                      enddo
                      term_1(:,1:1) = prod_rule(inter2,intermediate2,l_prod,l_prod)
                      term_1 = term_1 * coeff_a(i_cont1)
                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                       do k2dr=1,3
                        term_1_dervs(:,1:1,k2dr)= &
                           prod_rule(inter2_dervs(:,:,k2dr),intermediate2,l_prod,l_prod) &
                          +prod_rule(inter2,intermediate2_dervs(:,:,i,k2dr),l_prod,l_prod)
                        term_1_dervs(:,1:1,k2dr)=term_1_dervs(:,1:1,k2dr)*coeff_a(i_cont1)
                       enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                      endif
                      
                      ! produce term_2 = prod_rule(inter1,(yl_arr_grad(la-i_l)*gamma(la+lb-i_l)))

                      counter=1
                      intermediate2 = zero
                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) &
                                                 intermediate2_dervs(:,:,i,:)=zero 
                      do i_l=0,la
                         do i_m=1,2*i_l+1
                            intermediate2(:,counter) = &
                                 yl_grad(eq_a,eq_b,counter,i)*gamma_help(:,lb+i_l+1) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l

                            if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                            intermediate2_dervs(:,counter,i,:) = &
                                 spread(yl_dervs(counter,i,:),1,num)* &
                                 spread(gamma_help(:,lb+i_l+1) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l,2,3)

                            intermediate2_dervs(:,counter,i,:) = &
                            intermediate2_dervs(:,counter,i,:) + &
                              spread(two*(xa-xb),1,num)*&
                              spread(yl_grad(eq_a,eq_b,counter,i)*(-gamma_help(:,lb+i_l+2)) * &
                                 (two*fact1/fact0)**(i_l+lb)*(-one)**i_l*fact1/fact0,2,3)
FPP_TIMER_STOP(t_calc_2c_dervs)
                            endif

                            counter=counter+1
                         enddo
                      enddo
                      
                      ! produce term_3=prod_rule(inter1,(yl_arr(la-i_l)*grad_a(gamma(la+lb-i_l))))
                      counter=1
                      do i_l=0,la
                         do i_m=1,2*i_l+1                  
                            intermediate2(:,counter) = &
                                 intermediate2(:,counter) + &
                                 yl_arr(eq_a,eq_b,counter)*&
                                 (-gamma_help(:,lb+i_l+2))*two*fact1/fact0*(xa(i)-xb(i))*&
                                 (two*fact1/fact0)**(lb+i_l)*(-one)**i_l

                           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                            intermediate2_dervs(:,counter,i,:) = &
                                 intermediate2_dervs(:,counter,i,:) + &
                                 spread(yl_grad(eq_a,eq_b,counter,:),1,num)*&
                                 spread((-gamma_help(:,lb+i_l+2))*two*fact1/fact0*(xa(i)-xb(i))*&
                                 (two*fact1/fact0)**(lb+i_l)*(-one)**i_l,2,3)

                            intermediate2_dervs(:,counter,i,:) = &
                                 intermediate2_dervs(:,counter,i,:) + &
                                 spread(two*(xa-xb),1,num)*&
                                 spread(yl_arr(eq_a,eq_b,counter)*&
                                 (gamma_help(:,lb+i_l+3))*two*fact1/fact0*(xa(i)-xb(i))*&
                                 (two*fact1/fact0)**(lb+i_l)*(-one)**i_l*fact1/fact0,2,3)

                            intermediate2_dervs(:,counter,i,i) = &
                                 intermediate2_dervs(:,counter,i,i) + &
                                 yl_arr(eq_a,eq_b,counter)*&
                                 (-gamma_help(:,lb+i_l+2))*two*fact1/fact0*&
                                 (two*fact1/fact0)**(lb+i_l)*(-one)**i_l
FPP_TIMER_STOP(t_calc_2c_dervs)
                           endif
                            counter=counter+1
                         enddo
                      enddo

                      term_2 = prod_rule(inter1,intermediate2,l_prod,l_prod)
                      term_2 = term_2 *coeff_a(i_cont1) 

                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                      do k2dr=1,3
                       term_2_dervs(:,:,k2dr) = prod_rule(inter1_dervs(:,:,k2dr),intermediate2,l_prod,l_prod) &
                                              + prod_rule(inter1,intermediate2_dervs(:,:,i,k2dr),l_prod,l_prod)
                       term_2_dervs(:,:,k2dr) = term_2_dervs(:,:,k2dr) *coeff_a(i_cont1) 
                      enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                      endif
                      
                    ll_contrib_map: if(.true.) then
                      if (moving_a .and. eq_a == 1) then
                         fit_int_grad_a(:,:,i_ind1,i_ind2,i) = &
                              fit_int_grad_a(:,:,i_ind1,i_ind2,i) + &
                              unpack(two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))*( &
                              term_1(:,1) + term_2(:,1)),cutoff,zero)

                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                         do k2dr=1,3
                         fit_int_dervs_a(:,:,i_ind1,i_ind2,i,k2dr,eq_b) = &
                              fit_int_dervs_a(:,:,i_ind1,i_ind2,i,k2dr,eq_b) + &
                              unpack(two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
                              (term_1_dervs(:,1,k2dr) + term_2_dervs(:,1,k2dr)),cutoff,zero)
                              !-----------             ------------
                          enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif
                      endif
                      !if (na.ne.nb) then
                         if (moving_b .and. eq_b == 1) then
                            fit_int_grad_b(:,:,i_ind1,i_ind2,i) = &
                                 fit_int_grad_b(:,:,i_ind1,i_ind2,i) - &
                                 unpack( two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))*( &
                                 term_1(:,1) + term_2(:,1)),cutoff,zero)
                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
!if(nb.eq.1.and.na.ne.nb.and.la.gt.1.and.lb.gt.1) then
!   print*,'term_1',sum(term_1(:,1)),sum(term_1_dervs(:,1,3)),eq_a,magn_a(i_cont1),magn_b(i_cont2)
!endif
                 
                         do k2dr=1,3
                            fit_int_dervs_b(:,:,i_ind1,i_ind2,i,k2dr,eq_a) = &
                                 fit_int_dervs_b(:,:,i_ind1,i_ind2,i,k2dr,eq_a) + &
                                 unpack( two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &
                                 (term_1_dervs(:,1,k2dr) + term_2_dervs(:,1,k2dr)),cutoff,zero)
                         enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif
                         endif
                     endif ll_contrib_map
                      !endif
                      
                   enddo contributing2
                enddo contributing1
             enddo independent2
          enddo independent1
       enddo cartesian_ll

       deallocate( yl_arr_grad, term_2, term_1, &
            intermediate1, intermediate2, inter1, inter2, &
            STAT=alloc_stat)
        ASSERT(alloc_stat.eq.0)

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
       MEMLOG(-size(intermediate1_dervs)-size(intermediate2_dervs))
       MEMLOG(-size(inter1_dervs)-size(inter2_dervs)-size(yl_dervs))
       MEMLOG(-size(term_1_dervs)-size(term_2_dervs))
        deallocate(intermediate1_dervs,intermediate2_dervs, &
                   inter1_dervs,inter2_dervs,yl_dervs, &
                   term_1_dervs,term_2_dervs,STAT=cpksalloc(79))
        ASSERT(cpksalloc(79).eq.0)
        cpksalloc(79)=1
       endif
       
    endif
!!if(na.eq.1) then
!if(na.eq.1.and.na.ne.nb.and.la.gt.1.and.lb.gt.1) then
! print*,'ll fit_int_grad_a fit_int_dervs_a',na,nb,la,lb
!!! print*,sum(fit_int_grad_a),sum(fit_int_dervs_a(:,:,:,:,:,3,1)),3
! print*,sum(fit_int_grad_a(:,:,:,:,1)),sum(fit_int_dervs_a(:,:,:,:,1,3,:)),3
! print*,sum(fit_int_grad_a(:,:,:,:,2)),sum(fit_int_dervs_a(:,:,:,:,2,3,:)),3
! print*,sum(fit_int_grad_a(:,:,:,:,3)),sum(fit_int_dervs_a(:,:,:,:,3,3,:)),3
!! !print*,sum(fit_int_grad_a(:,:,1,1,1)),sum(fit_int_dervs_a(:,:,1,1,1,1,1))
!endif
    
    ! Now do contractions
    if (moving_a) then

       call fitcontract_2c_grad(fit_int_grad_a,ma,1) !(ll) (1)
       MEMLOG(-size(fit_int_grad_a))
       deallocate( fit_int_grad_a, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!!!      if(la.gt.0.and.lb.gt.0) &
FPP_TIMER_START(t_contract_2c_dervs)
       call fitcontract_2c_dervs(fit_int_dervs_a,ma,mb) !(ll)(1)
FPP_TIMER_STOP(t_contract_2c_dervs)
       MEMLOG(-size(fit_int_dervs_a))
       deallocate( fit_int_dervs_a, STAT=cpksalloc(74))
       ASSERT(cpksalloc(74).eq.0)
       cpksalloc(74)=1
     endif

    endif

!!if(nb.eq.1) then
!if(nb.eq.1.and.na.ne.nb.and.la.gt.1.and.lb.gt.1) then
! print*,'ll fit_int_grad_b fit_int_dervs_b',na,nb,la,lb
!!! print*,sum(fit_int_grad_b),sum(fit_int_dervs_b(:,:,:,:,:,3,1))
! print*,sum(fit_int_grad_b(:,:,:,:,1)),sum(fit_int_dervs_b(:,:,:,:,1,3,:))
! print*,sum(fit_int_grad_b(:,:,:,:,2)),sum(fit_int_dervs_b(:,:,:,:,2,3,:))
! print*,sum(fit_int_grad_b(:,:,:,:,3)),sum(fit_int_dervs_b(:,:,:,:,3,3,:))
!!print*,sum(fit_int_grad_b(:,:,1,1,1)),sum(fit_int_dervs_b(:,:,1,1,1,1,1))
!endif

    if (moving_b) then
       call fitcontract_2c_grad(fit_int_grad_b,mb,2) !(ll) (2)
       MEMLOG(-size(fit_int_grad_b))
       deallocate( fit_int_grad_b, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)


     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!       print*,'fit_int_dervs_b'
!       do k2dr=1,3
!       print*,sum(fit_int_dervs_b(:,:,:,:,1,k2dr,1)), &
!                 sum(fit_int_dervs_b(:,:,:,:,2,k2dr,1)), &
!                    sum(fit_int_dervs_b(:,:,:,:,3,k2dr,1))
!       enddo

!!!      if(la.gt.0.and.lb.gt.0) &
FPP_TIMER_START(t_contract_2c_dervs)
       call fitcontract_2c_dervs(fit_int_dervs_b,mb,ma) !(ll) (2)
FPP_TIMER_STOP(t_contract_2c_dervs)
       MEMLOG(-size(fit_int_dervs_b))
       deallocate( fit_int_dervs_b, STAT=cpksalloc(75))
       ASSERT(cpksalloc(75).eq.0)
       cpksalloc(75)=1
     endif

    endif

    ! deallocation
    deallocate(  gamma_help,  gamma_arg,  yl_arr, yl_grad, &
         n_contributing_a,  n_contributing_b, &
         STAT=alloc_stat)
         ASSERT(alloc_stat.eq.0)
    call dump_exponent_data()

    
  end subroutine ll_calc_grad


  subroutine lr2_calc_grad(na,la,nb,lb)

    ! Purpose: calculate primitives for one angular momentum
    !          greater equal zero and the other being r2-type
    !          (l=-1).
    !          Symmetry adaption included.
    !------------ Declaration of formal parameters ---------------

   use cpksdervs_matrices,only: cpksalloc
   use calc3c_switches

    integer(kind=i4_kind),intent(inout)  :: na,nb
    integer(kind=i4_kind),intent(inout)  :: la,lb

    !** End of interface *****************************************
    ! atomic positions



    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    real(kind=r8_kind),allocatable   :: yl_dervs(:,:,:),yl_grad(:,:,:,:)
    ! help factors
    real(kind=r8_kind),allocatable  :: gamma_help(:,:)
    real(kind=r8_kind),allocatable  :: yl_arr(:,:,:)
    real(kind=r8_kind),allocatable  :: inter1(:,:),inter2(:,:),&
         inter1_dervs(:,:,:),inter2_dervs(:,:,:),&
         intermediate1(:),intermediate2(:), &
         intermediate1_dervs(:,:,:),intermediate2_dervs(:,:,:)
    real(kind=r8_kind),allocatable  :: args(:,:)

    integer(kind=i4_kind),allocatable  :: n_contributing_a(:)
    integer(kind=i4_kind),allocatable  :: n_contributing_b(:)
    integer(kind=i4_kind),pointer      :: magn_a(:),eq_atom_a(:)
    integer(kind=i4_kind),pointer      :: magn_b(:),eq_atom_b(:)
    real(kind=r8_kind),pointer         :: coeff_b(:)
    real(kind=r8_kind),pointer         :: coeff_a(:)

    ! fitintegral that enters contraction

    real(kind=r8_kind),allocatable  :: &
         fit_int_grad_a(:,:,:,:,:),fit_int_dervs_a(:,:,:,:,:,:,:), &
         fit_int_grad_b(:,:,:,:,:),fit_int_dervs_b(:,:,:,:,:,:,:) 

    integer(kind=i4_kind) :: eq_a,eq_b,n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b,lm,i_sum
    integer(kind=i4_kind) :: alloc_stat
    integer(kind=i4_kind) :: i_ind1,i_ind2,i_cont1,i_cont2,&
         l_max,max_order,i
    integer(kind=i4_kind) :: ma, mb, begin,llm,k2dr
    logical :: moving_a,moving_b
    intrinsic max
    !------------ Executable code --------------------------------

    ma = unique_atoms(na)%moving_atom
    mb = unique_atoms(nb)%moving_atom
    moving_a = ma > 0
    moving_b = mb > 0
    if (.not.moving_a .and. .not.moving_b) return
    call get_exponent_data(na,la,nb,lb)

    l_max=max(la,lb)

    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms

    if(lb.eq.-1) then
       n_indep_a = unique_atoms(na)%symadapt_partner(i_ir,la)%&
            &n_independent_fcts 
       n_indep_b= 1
       allocate(  n_contributing_a(n_indep_a),n_contributing_b(1), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       if (n_indep_a/=0) then
          n_contributing_a=unique_atoms(na)%symadapt_partner(i_ir,la)%&
               &symadapt(:,1)%n_fcts
       endif
       n_contributing_b(1)=n_equal_b
    else
       n_indep_b = unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
            &n_independent_fcts 
       n_indep_a= 1 
       allocate( n_contributing_b(n_indep_b),n_contributing_a(1), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

       if (n_indep_b/=0) then
          n_contributing_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
               &symadapt(:,1)%n_fcts
       endif
       n_contributing_a(1)=n_equal_a
    endif

    if (moving_a) then
       allocate( fit_int_grad_a(naexps,nbexps,n_indep_a,n_indep_b,3), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_a = zero

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
       allocate( fit_int_dervs_a(naexps,nbexps,n_indep_a,n_indep_b,3,3,n_equal_b), &
                 STAT=cpksalloc(60))
       ASSERT(cpksalloc(60).eq.0)
       MEMLOG(size(fit_int_dervs_a))
       fit_int_dervs_a = zero
     endif
    endif

    if (moving_b) then
       allocate( fit_int_grad_b(naexps,nbexps,n_indep_a,n_indep_b,3), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_b = zero

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
       allocate( fit_int_dervs_b(naexps,nbexps,n_indep_a,n_indep_b,3,3,n_equal_a), &
                 STAT=cpksalloc(67))
       ASSERT(cpksalloc(67).eq.0)
       MEMLOG(size(fit_int_dervs_b))
       fit_int_dervs_b = zero
     endif

    endif

    ! To avoid doublication of code double DO-loops over all pairs of equal
    ! atoms are used, though only the following pairs are final processed
    !   a) if moving_a : the pairs (1,eq_b) and
    !   b) if moving_b : the pairs (eq_a,1) .

    if (l_max.gt.0) then
       allocate(yl_arr(n_equal_a,n_equal_b,(l_max+1)**2), &
                yl_grad(n_equal_a,n_equal_b,(l_max+1)**2,3), &
                                         STAT=cpksalloc(65))
        ASSERT(cpksalloc(65).eq.0)

                      yl_grad = zero
  
       ! do a pre-calculation of solid_harmonics------------------------

       begin = 1
       if (moving_a) then
          allocate(args(n_equal_b,3),STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          eq_a = 1
          args =  transpose (spread(unique_atoms(na)%position(:,eq_a),2,n_equal_b) - &
               unique_atoms(nb)%position)

          yl_arr(eq_a,:,:) = solid_harmonics_calc(l_max,args)

                      do eq_b=1,n_equal_b

                       do i=1,3
                       do llm=1,(l_max+1)**2
                         do i_sum=1,solhrules_differential(yl_map(i),llm)%n_summands
                            yl_grad(eq_a,eq_b,llm,i) =  yl_grad(eq_a,eq_b,llm,i) + &
                            solhrules_differential(yl_map(i),llm)%coef(i_sum)* &
                                                   yl_arr(eq_a,eq_b,&
                                                   solhrules_differential(yl_map(i),llm)%lm_sh(i_sum))
                         enddo
                       enddo 
                    enddo
             enddo

       deallocate( args, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

          begin = 2
       endif

       if (moving_b) then
          allocate(args(1,3),STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
          do eq_a = begin,n_equal_a
             args(1,:) = unique_atoms(na)%position(:,eq_a) - &
                         unique_atoms(nb)%position(:,1)
             yl_arr(eq_a,1:1,:) = solid_harmonics_calc(l_max,args)

                       do i=1,3
                       do llm=1,(l_max+1)**2
                         do i_sum=1,solhrules_differential(yl_map(i),llm)%n_summands
                            yl_grad(eq_a,1,llm,i) =  yl_grad(eq_a,1,llm,i) + &
                            solhrules_differential(yl_map(i),llm)%coef(i_sum)* &
                                                   yl_arr(eq_a,1,&
                                                   solhrules_differential(yl_map(i),llm)%lm_sh(i_sum))
                         enddo
                       enddo 
                    enddo
          enddo
          deallocate( args, STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)
       endif
       ! ------------------------------------------------------------------
    endif


    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! fact2 = a/b  (if lb=-1)
    ! fact2 = b/a  (if la=-1)

    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
     max_order=max(4,l_max+4)
    else
     max_order=max(3,l_max+3)
    endif

    allocate( gamma_help(num,max_order), STAT=alloc_stat)
     ASSERT(alloc_stat.eq.0)
     gamma_help=zero


    if (lb.eq.-1) then     ! ========== lb = r2 ==============================
       if (la.eq.0) then   ! la = 0 and lb = r2

          if (moving_a) then
             allocate( inter1(num,3), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             inter1 = zero

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( inter1_dervs(num,3,3), STAT=cpksalloc(61))
             ASSERT(cpksalloc(61).eq.0)
             MEMLOG(size(inter1_dervs))
             inter1_dervs = zero
           endif

          endif

          if (moving_b) then
             allocate( inter2(num,3), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             inter2 = zero

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( inter2_dervs(num,3,3), STAT=cpksalloc(68))
             ASSERT(cpksalloc(68).eq.0)
             MEMLOG(size(inter2_dervs))
             inter2_dervs = zero
           endif

          endif


          a_eq: do eq_a=1,n_equal_a
             b_eq: do eq_b=1,n_equal_b
 
                if ( (eq_a /= 1 .or. .not.moving_a) .and. &
                     (eq_b /= 1 .or. .not.moving_b) ) cycle

                xa=unique_atoms(na)%position(:,eq_a)
                xb=unique_atoms(nb)%position(:,eq_b)

                if ((na.eq.nb).and.(eq_a.eq.eq_b)) then
                   arg = zero
                else
                   arg=sum((xa-xb)**2)
                endif

               if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
                gamma_help(:,1:4) = gamma(4,fact1/fact0*arg)
               else
                gamma_help(:,1:3) = gamma(3,fact1/fact0*arg)
               endif

                cart_grad: do i=1,3
                   if (moving_a .and. eq_a == 1) then

                      inter1(:,i) = inter1(:,i)  &
                          -(three*gamma_help(:,2) + &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3))* (xa(i)-xb(i))

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      inter1_dervs(:,i,:) = &! inter1_dervs(:,i,:)  &          ! d/da gamma_help(:,2)
                          -spread(two*(xa(i)-xb(i))*(xa-xb),1,num)* &
                            spread(three*(-gamma_help(:,3))*fact1/fact0,2,3)  
                            !             ^             ^
                      inter1_dervs(:,i,i)=inter1_dervs(:,i,i)  &
                        - three*gamma_help(:,2)

                      inter1_dervs(:,i,:) = inter1_dervs(:,i,:) - &           ! d/da arg
                           spread(four*(xa(i)-xb(i))*(xa-xb),1,num)* &
                            spread(fact1*fact2/fact0*gamma_help(:,3),2,3)

                      inter1_dervs(:,i,:) = inter1_dervs(:,i,:) - &           ! d/da gamma_help(:,3)
                           spread((xa(i)-xb(i))*four*(xa-xb),1,num)* &
                           spread(fact1/fact0* &
                                  fact1*fact2/fact0*arg*(-gamma_help(:,4)),2,3)
                                  !                                    ^

                      inter1_dervs(:,i,i) = inter1_dervs(:,i,i) - &
                           ( two*fact1*fact2/fact0*arg*gamma_help(:,3))

                       do k2dr=1,3
                         fit_int_dervs_a(:,:,1,1,i,k2dr,eq_b) = two*pi*pi*sqrt(pi)*&
                         unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))* &
                                                 inter1_dervs(:,i,k2dr),cutoff,zero)
                       enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif

                   endif

                   if (moving_b .and. eq_b == 1) then
                      inter2(:,i) = inter2(:,i) + &
                           (three*gamma_help(:,2) + &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3))*&
                           (xa(i)-xb(i))

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      inter2_dervs(:,i,:) = &!inter2_dervs(:,i,:) &
                      -spread(two*(xa-xb),1,num)* &
                       spread(three*(-gamma_help(:,3))*(xa(i)-xb(i))*fact1/fact0,2,3)

                      inter2_dervs(:,i,i) = inter2_dervs(:,i,i) - three*gamma_help(:,2)

                      inter2_dervs(:,i,:) = inter2_dervs(:,i,:) - &
                       spread(two*(xa-xb),1,num)* &
                       spread( (two*fact1*fact2/fact0*gamma_help(:,3))*(xa(i)-xb(i)),2,3)

                      inter2_dervs(:,i,:) = inter2_dervs(:,i,:) - &
                       spread( two*(xa-xb),1,num)* &
                       spread( (two*fact1*fact2/fact0*arg*(-gamma_help(:,4)))*(xa(i)-xb(i))* &
                            fact1/fact0,2,3)

                      inter2_dervs(:,i,i) = inter2_dervs(:,i,i) - &
                           ( two*fact1*fact2/fact0*arg*gamma_help(:,3))

                      do k2dr=1,3
                        fit_int_dervs_b(:,:,1,1,i,k2dr,eq_a) = two*pi*pi*sqrt(pi)*&
                          unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))* &
                                                 inter2_dervs(:,i,k2dr),cutoff,zero)
                      enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif

                   endif
                enddo cart_grad


             enddo b_eq
          enddo a_eq

          do i=1,3
             if (moving_a) then
              fit_int_grad_a(:,:,1,1,i) = two*pi*pi*sqrt(pi)*&
                   unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))*inter1(:,i),cutoff,zero)
            endif
             !if(na.ne.nb) then
                if (moving_b) then
                 fit_int_grad_b(:,:,1,1,i) = two*pi*pi*sqrt(pi)*&
                    unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))*inter2(:,i),cutoff,zero)


                endif
             !endif
          enddo

          if (moving_a) then
             deallocate( inter1, STAT=alloc_stat)
!       print*,' fit_int_grad_a lr2',na,nb,la,lb
!       print*,sum(fit_int_grad_a(:,:,1,1,1)),sum(fit_int_dervs_a(:,:,1,1,1,1,1))
           ASSERT(alloc_stat.eq.0)

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             MEMLOG(-size(inter1_dervs))
             deallocate( inter1_dervs, STAT=cpksalloc(61) )
!       print*,' fit_int_grad_b lr2',na,nb,la,lb
!       print*,sum(fit_int_grad_b(:,:,1,1,1)),sum(fit_int_dervs_b(:,:,1,1,1,1,1))
           ASSERT(cpksalloc(61).eq.0)
           cpksalloc(61)=1
           endif

          endif

          if (moving_b) then
             deallocate( inter2, STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             MEMLOG(-size(inter2_dervs))
             deallocate( inter2_dervs, STAT=cpksalloc(68))
             ASSERT(cpksalloc(68).eq.0)
             cpksalloc(68)=1
           endif
          endif


       else ! la > 0 AND lb=r2

          if (moving_a) then
             allocate( intermediate1(num), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             intermediate1 = zero

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( intermediate1_dervs(num,3,3), STAT=cpksalloc(62))
             ASSERT(cpksalloc(62).eq.0)
             MEMLOG(size(intermediate1_dervs))
             intermediate1_dervs = zero
           endif

          endif

          if (moving_b) then
             allocate( intermediate2(num), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             intermediate2 = zero

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( intermediate2_dervs(num,3,3), STAT=cpksalloc(69))
             ASSERT(cpksalloc(69).eq.0)
             MEMLOG(size(intermediate2_dervs))
             intermediate2_dervs = zero
           endif

          endif

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           allocate(yl_dervs((la+1)**2,3,3),stat=cpksalloc(63))
           ASSERT(cpksalloc(63).eq.0)
           MEMLOG(size(yl_dervs))
           yl_dervs=zero
          endif

          cartesian_lr2: do i=1,3
             
             indeps_a: do i_ind1=1,n_indep_a

                   eq_atom_a => &
                        unique_atoms(na)%symadapt_partner(i_ir,la)%&
                        &symadapt(i_ind1,i_ir)%I_equal_atom
                   coeff_a => &
                   unique_atoms(na)%symadapt_partner(i_ir,la)%symadapt(i_ind1,i_ir)%c  
                   magn_a => &
                   unique_atoms(na)%symadapt_partner(i_ir,la)%symadapt(i_ind1,i_ir)%m

                cont1_lr2 : do i_cont1=1,n_contributing_a(i_ind1) 
                   eq_a=eq_atom_a(i_cont1)
                   equal_lr2: do eq_b=1,unique_atoms(nb)%n_equal_atoms

                      if ( (eq_atom_a(i_cont1) /= 1 .or. .not.moving_a) .and. &
                           (eq_b /= 1 .or. .not.moving_b) ) cycle

                      ! calculate gradients of solid harmonics
                      lm = la**2 + magn_a(i_cont1)


                      xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
                      xb = unique_atoms(nb)%position(:,eq_b)
                      arg = sum((xa-xb)**2)
!                     print*,'arg =',arg,two*(xa(1)-xb(1)),two

                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                       yl_dervs(:,i,:)=zero
                      if (.not.((na.eq.nb).and.(eq_atom_a(i_cont1).eq.eq_b))) then
                       do llm=1,(la+1)**2
                       do k2dr=1,3
                         do i_sum=1,solhrules_differential(yl_map(k2dr),llm)%n_summands
                            yl_dervs(llm,i,k2dr) =  yl_dervs(llm,i,k2dr) + &
                            solhrules_differential(yl_map(k2dr),llm)%coef(i_sum)* &
                             yl_grad(eq_a,eq_b,solhrules_differential(yl_map(k2dr),llm)%lm_sh(i_sum),i)
                         enddo
                       enddo 
                       enddo 
                      endif
                       gamma_help(:,1:la+4) = gamma(la+4,fact1/fact0*arg)
FPP_TIMER_STOP(t_calc_2c_dervs)
                      else
                       gamma_help(:,1:la+3) = gamma(la+3,fact1/fact0*arg)
                      endif
                    
                      if (moving_a .and. eq_atom_a(i_cont1) == 1) then

                         intermediate1 =  intermediate1 + ( &                         ! (1)
                              (three/two-fact2*real(la-1,kind=r8_kind))*( &
                              yl_grad(eq_a,eq_b,lm,i)*gamma_help(:,la+1) &
                                      - yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                             gamma_help(:,la+2) *two*fact1/fact0*(xa(i)-xb(i)) &
                                ) )*coeff_a(i_cont1)

                         intermediate1 = intermediate1 + ( &                         ! (2)
                              fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                              gamma_help(:,la+2))* coeff_a(i_cont1)
                         
                         intermediate1 = intermediate1 + (&                          ! (3)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)*(&
                              two*(xa(i)-xb(i))*gamma_help(:,la+2) - &
                              arg*gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i))))*&
                              coeff_a(i_cont1)
                      
                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                         intermediate1_dervs(:,i,:) =                               &       ! (1.1)
                             +spread(yl_dervs(lm,i,:),1,num) * &
                                  spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                                     gamma_help(:,la+1) * coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &     ! (1.2)
                            spread( two*(xa-xb),1,num)*  &
                              spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                              yl_grad(eq_a,eq_b,lm,i)* &
                               (-gamma_help(:,la+2)) * &
                                 coeff_a(i_cont1) * fact1/fact0, 2,3)
                              
                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:)  &      ! (1.3)
                             -spread(yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num) * &
                              spread( (three/two-fact2*real(la-1,kind=r8_kind)) &
                                *gamma_help(:,la+2)*two*fact1/fact0*(xa(i)-xb(i)) &
                                      *coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:)   &
                            -spread( two*(xa-xb),1,num)*&                                  ! (1.4)
                              spread( (three/two-fact2*real(la-1,kind=r8_kind))*( &
                                yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                                (-gamma_help(:,la+3)) &
                                  *two*fact1/fact0*(xa(i)-xb(i)) &
                                    )*coeff_a(i_cont1)*fact1/fact0,2,3)
                                                                                    !-----------
                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) +  &     ! (1.5)
                              (three/two-fact2*real(la-1,kind=r8_kind))*( &
                                    - yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                              gamma_help(:,la+2)*two*fact1/fact0 )*coeff_a(i_cont1)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &      ! (2.1)
                              spread( yl_dervs(lm,i,:),1,num) * spread( &
                              fact1*fact2/fact0*arg*gamma_help(:,la+2)*coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &      ! (2.2)
                         spread(two*(xa-xb),1,num)* &
                         spread(fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*&
                                gamma_help(:,la+2)* coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) +  &      ! (2.3)
                              spread(two*(xa-xb),1,num)* &
                              spread( fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                              (-gamma_help(:,la+3))* coeff_a(i_cont1) * &
                                      fact1/fact0,2,3) 

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:)+ &        ! (3.1)
                           spread(yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num)* &      
                           spread(fact1*fact2/fact0*&
                            two*(xa(i)-xb(i))*gamma_help(:,la+2)*coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) + &       ! (3.2)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                              two*gamma_help(:,la+2)*coeff_a(i_cont1)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &       ! (3.3)
                          spread(two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                                  two*(xa(i)-xb(i))*(-gamma_help(:,la+3))*coeff_a(i_cont1)*&
                                  fact1/fact0,2,3)
                                  
                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &       ! (3.4)
                           spread(yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num)* &
                           spread( fact1*fact2/fact0* &
                                   (-arg*gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i)) )*&
                                   coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &        ! (3.5)
                           spread(two*(xa-xb),1,num)* &
                           spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                                   (-gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i)) )*&
                                   coeff_a(i_cont1),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &         ! (3.6)
                           spread(two*(xa-xb),1,num)* &
                           spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                                   (-arg*(-gamma_help(:,la+4))*two*fact1/fact0*(xa(i)-xb(i)) )*&
                                   coeff_a(i_cont1)*fact1/fact0,2,3)

                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) + &          ! (3.7)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                               (-arg*gamma_help(:,la+3)*two*fact1/fact0 )*coeff_a(i_cont1)

                         do k2dr=1,3
                          fit_int_dervs_a(:,:,i_ind1,1,i,k2dr,eq_b) = &
                             fit_int_dervs_a(:,:,i_ind1,1,i,k2dr,eq_b) + &
                             unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**la*&
                             one/(fact1*fact0*sqrt(fact0))*intermediate1_dervs(:,i,k2dr), cutoff,zero)
                          enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif

                      endif

                      if (moving_b .and. eq_b == 1) then
                         intermediate2 = intermediate2 + ( &                     ! (1)
                              (three/two-fact2*real(la-1,kind=r8_kind))*( &
                              -yl_grad(eq_a,eq_b,lm,i)*gamma_help(:,la+1)  &
                                      + yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                              gamma_help(:,la+2)*two*fact1/fact0*(xa(i)-xb(i))))*coeff_a(i_cont1)

                         intermediate2 = intermediate2 - ( &                     ! (2)
                              fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                              gamma_help(:,la+2))* coeff_a(i_cont1)

                         intermediate2 = intermediate2 + ( &                     ! (3)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)*( &
                              -two*(xa(i)-xb(i))*gamma_help(:,la+2) + &
                              arg*gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i))))*&
                              coeff_a(i_cont1)

                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                         intermediate2_dervs(:,i,:) =                                  &       ! (1.1)
                           spread( yl_dervs(lm,i,:),1,num)* &
                           spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                              gamma_help(:,la+1)*coeff_a(i_cont1),2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &       ! (1.2)
                           spread( two*(xa-xb),1,num)* &
                           spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                              yl_grad(eq_a,eq_b,lm,i)*(-gamma_help(:,la+2))*coeff_a(i_cont1)* &
                              fact1/fact0,2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &        ! (1.3)
                           spread(yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num)* &
                           spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                              gamma_help(:,la+2)*two*fact1/fact0*(xa(i)-xb(i))*coeff_a(i_cont1),2,3)



                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &        ! (1.4)
                           spread( two*(xa-xb),1,num)* &
                           spread( (three/two-fact2*real(la-1,kind=r8_kind))* &
                               yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                              (-gamma_help(:,la+3))*two*fact1/fact0*(xa(i)-xb(i))*coeff_a(i_cont1)* &
                               fact1/fact0,2,3)

                         intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) -  &         ! (1.5)
                              (three/two-fact2*real(la-1,kind=r8_kind))* &
                               yl_arr(eq_atom_a(i_cont1),eq_b,lm)*&
                              gamma_help(:,la+2)*two*fact1/fact0*coeff_a(i_cont1)


                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &         ! (2.1)
                           spread(yl_dervs(lm,i,:),1,num)* &
                           spread( fact1*fact2/fact0*arg*&
                                   gamma_help(:,la+2)* coeff_a(i_cont1),2,3)


                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &         ! (2.2)
                           spread( two*(xa-xb),1,num)*&
                           spread( fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*&
                              gamma_help(:,la+2)* coeff_a(i_cont1),2,3)


                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &         ! (2.3)
                           spread( two*(xa-xb),1,num)* &
                           spread( fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                                   (-gamma_help(:,la+3))* coeff_a(i_cont1)*fact1/fact0,2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &         ! (3.1)
                           spread( yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num)* &
                           spread( fact1*fact2/fact0* &
                              two*(xa(i)-xb(i))*gamma_help(:,la+2) * coeff_a(i_cont1),2,3)


                         intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) +  &         ! (3.2)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                              two*gamma_help(:,la+2) * coeff_a(i_cont1)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &         ! (3.3)
                           spread( two*(xa-xb),1,num)* &
                           spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                              two*(xa(i)-xb(i))*(-gamma_help(:,la+3))*coeff_a(i_cont1)* &
                              fact1/fact0,2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &         ! (3.4)
                           spread( yl_grad(eq_atom_a(i_cont1),eq_b,lm,:),1,num)* &
                           spread( fact1*fact2/fact0* &
                              arg*gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i))*&
                              coeff_a(i_cont1),2,3 )


                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &         ! (3.5)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                              gamma_help(:,la+3)*two*fact1/fact0*(xa(i)-xb(i))*&
                              coeff_a(i_cont1),2,3 )

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &          ! (3.6)
                           spread( two*(xa-xb),1,num)* &
                           spread( fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                              arg*(-gamma_help(:,la+4))*two*fact1/fact0*(xa(i)-xb(i))*&
                              coeff_a(i_cont1)*fact1/fact0,2,3 )

                         intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) -  &          ! (3.7)
                              fact1*fact2/fact0*yl_arr(eq_atom_a(i_cont1),eq_b,lm)* &
                              arg*gamma_help(:,la+3)*two*fact1/fact0*coeff_a(i_cont1)


                         do k2dr=1,3
                         fit_int_dervs_b(:,:,i_ind1,1,i,k2dr,eq_a) =  &
                           fit_int_dervs_b(:,:,i_ind1,1,i,k2dr,eq_a) + &
                            unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**la*&
                            one/(fact1*fact0*sqrt(fact0))*intermediate2_dervs(:,i,k2dr), cutoff,zero)
                         enddo
FPP_TIMER_STOP(t_calc_2c_dervs)

                       endif

                      endif
                         

                   enddo equal_lr2
                enddo cont1_lr2

                if (moving_a) then
                fit_int_grad_a(:,:,i_ind1,1,i) = fit_int_grad_a(:,:,i_ind1,1,i) + &
                     unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**la*&
                     one/(fact1*fact0*sqrt(fact0))*intermediate1,&
                     cutoff,zero)
                   intermediate1 = zero
                endif

                !if(na.ne.nb) then
                   if (moving_b) then
                   fit_int_grad_b(:,:,i_ind1,1,i) = fit_int_grad_b(:,:,i_ind1,1,i) + &
                        unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**la*&
                        one/(fact1*fact0*sqrt(fact0))*intermediate2,&
                        cutoff,zero) 
                      intermediate2 = zero
                   endif
                !endif

             enddo indeps_a
          enddo cartesian_lr2 !(i)

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then                           !!!!!!!
           MEMLOG(-size(yl_dervs))
           deallocate(yl_dervs,stat=cpksalloc(63))
           ASSERT(cpksalloc(63).eq.0)
           cpksalloc(63)=1
          endif

          if (moving_a) then

             deallocate( intermediate1, STAT=alloc_stat)
           ASSERT(alloc_stat.eq.0)

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!    print*,'lr2= fit_int_grad_a',na,nb
!   print*,sum(fit_int_grad_a(:,:,:,:,1)),sum(fit_int_dervs_a(:,:,:,:,1,1,1))
             MEMLOG(-size(intermediate1_dervs))
             deallocate( intermediate1_dervs, STAT=cpksalloc(62))
           ASSERT(cpksalloc(62).eq.0)
           cpksalloc(62)=1
           endif

          endif
          if (moving_b) then
             deallocate( intermediate2, STAT=alloc_stat)
           ASSERT(alloc_stat.eq.0)
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!    print*,'lr2= fit_int_grad_b',na,nb
!   print*,sum(fit_int_grad_b(:,:,:,:,1)),sum(fit_int_dervs_b(:,:,:,:,1,1,1))
             MEMLOG(-size(intermediate2_dervs))
             deallocate( intermediate2_dervs, STAT=cpksalloc(69))
             ASSERT(cpksalloc(69).eq.0)
             cpksalloc(69)=1
           endif
          endif
       endif

    else ! la=r2

       if (lb.eq.0) then

          if (moving_a) then
             allocate( inter1(num,3), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             inter1 = zero
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( inter1_dervs(num,3,3), STAT=cpksalloc(70))
             ASSERT(cpksalloc(70).eq.0)
             inter1_dervs=zero
           endif
          endif

          if (moving_b) then
             allocate( inter2(num,3), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             inter2 = zero
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( inter2_dervs(num,3,3), STAT=cpksalloc(71))
             ASSERT(cpksalloc(71).eq.0)
             inter2_dervs=zero
           endif
          endif
          
          a_eq_lar2: do eq_a=1,n_equal_a
             b_eq_lar2: do eq_b=1,n_equal_b

                if ( (eq_a /= 1 .or. .not.moving_a) .and. &
                     (eq_b /= 1 .or. .not.moving_b) ) cycle

                xa=unique_atoms(na)%position(:,eq_a)
                xb=unique_atoms(nb)%position(:,eq_b)
                arg=sum((xa-xb)**2)
              
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
                gamma_help(:,1:4) = gamma(4,fact1/fact0*arg)
              else
                gamma_help(:,1:3) = gamma(3,fact1/fact0*arg)
              endif

                grad_i_lar2: do i=1,3

                   if (moving_a .and. eq_a == 1) then
                      inter1(:,i) = inter1(:,i) - &
                           (three*gamma_help(:,2) + &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3))*&
                           (xa(i)-xb(i))

                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      inter1_dervs(:,i,:) =  &
                      -spread( two*(xa-xb),1,num)* &
                       spread( three*(-gamma_help(:,3))*(xa(i)-xb(i))*fact1/fact0,2,3)

                      inter1_dervs(:,i,i) = inter1_dervs(:,i,i)-three*gamma_help(:,2)

                      inter1_dervs(:,i,:) = inter1_dervs(:,i,:) - &
                       spread( two*(xa-xb),1,num)* &
                       spread( two*fact1*fact2/fact0*gamma_help(:,3)*&
                               (xa(i)-xb(i)),2,3)

                      inter1_dervs(:,i,:) = inter1_dervs(:,i,:) - &
                       spread( two*(xa-xb),1,num)* &
                       spread( two*fact1*fact2/fact0*arg*(-gamma_help(:,4))*&
                               (xa(i)-xb(i))*fact1/fact0,2,3)

                      inter1_dervs(:,i,i) = inter1_dervs(:,i,i) - &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3)

                      do k2dr=1,3
                      fit_int_dervs_a(:,:,1,1,i,k2dr,eq_b) = two*pi*pi*sqrt(pi)*&
                         unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))* &
                                inter1_dervs(:,i,k2dr),cutoff,zero)
                      enddo
FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif
                   endif

                   if (moving_b .and. eq_b == 1) then
                      inter2(:,i) = inter2(:,i) + &
                           (three*gamma_help(:,2) + &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3))*&
                           (xa(i)-xb(i))
                    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                      inter2_dervs(:,i,:) =  &
                       -spread( two*(xa-xb),1,num)* &
                        spread( three*(-gamma_help(:,3))*(xa(i)-xb(i))*fact1/fact0,2,3)

                      inter2_dervs(:,i,i) = inter2_dervs(:,i,i)-three*gamma_help(:,2)

                      inter2_dervs(:,i,:) = inter2_dervs(:,i,:) - &
                        spread( two*(xa-xb),1,num)* &
                        spread( two*fact1*fact2/fact0*gamma_help(:,3)*&
                           (xa(i)-xb(i)),2,3)

                      inter2_dervs(:,i,:) = inter2_dervs(:,i,:) - &
                        spread( two*(xa-xb),1,num)* &
                        spread( two*fact1*fact2/fact0*arg*(-gamma_help(:,4))*&
                           (xa(i)-xb(i))*fact1/fact0,2,3)

                      inter2_dervs(:,i,i) = inter2_dervs(:,i,i)- &
                           two*fact1*fact2/fact0*arg*gamma_help(:,3)

                      do k2dr=1,3
                      fit_int_dervs_b(:,:,1,1,i,k2dr,eq_a) = two*pi*pi*sqrt(pi)*&
                         unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))* &
                                inter2_dervs(:,i,k2dr),cutoff,zero)
                      enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
                    endif
                   endif
                enddo grad_i_lar2
!             print*,'inter1_dervs lar2'
!             do k2dr=1,3
!              print*, sum(inter1_dervs(:,1,k2dr)), &
!                          sum(inter1_dervs(:,2,k2dr)), &
!                             sum(inter1_dervs(:,3,k2dr))
!             enddo

             enddo b_eq_lar2
          enddo a_eq_lar2

          do i=1,3
             if (moving_a) then
             fit_int_grad_a(:,:,1,1,i) = two*pi*pi*sqrt(pi)*&
                  unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))*inter1(:,i),cutoff,zero)
             endif
             !if(na.ne.nb) then
                if (moving_b) then
                fit_int_grad_b(:,:,1,1,i) = two*pi*pi*sqrt(pi)*&
                     unpack(fact1/(fact0*fact1*fact0*sqrt(fact0))*inter2(:,i),cutoff,zero)
                endif
             !endif
          enddo

          if (moving_a) then
             deallocate (inter1,STAT=alloc_stat)
           ASSERT(alloc_stat.eq.0)

           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             deallocate (inter1_dervs,STAT=cpksalloc(70))
             ASSERT(cpksalloc(70).eq.0)
             cpksalloc(70)=1
           endif

          endif
          if (moving_b) then
             deallocate (inter2,STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             deallocate (inter2_dervs,STAT=cpksalloc(71))
             ASSERT(cpksalloc(71).eq.0)
             cpksalloc(71)=1
           endif
          endif
!      print*,'lr2x fit_int_grad_a',na,nb
!      print*,sum(fit_int_grad_a(:,:,:,:,1)),sum(fit_int_dervs_a(:,:,:,:,1,1,1))
!      print*,'lr2x fit_int_grad_b',na,nb
!      print*,sum(fit_int_grad_b(:,:,:,:,1)),sum(fit_int_dervs_b(:,:,:,:,1,1,1))


       elseif (lb.gt.0) then
          if (moving_a) then
             allocate( intermediate1(num), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             intermediate1 = zero
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( intermediate1_dervs(num,3,3), STAT=cpksalloc(72))
             ASSERT(cpksalloc(72).eq.0)
             MEMLOG(size(intermediate1_dervs))
            intermediate1_dervs=zero
           endif
          endif

          if (moving_b) then
             allocate( intermediate2(num), STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
             intermediate2 = zero
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             allocate( intermediate2_dervs(num,3,3), STAT=cpksalloc(73))
             ASSERT(cpksalloc(73).eq.0)
             MEMLOG(size(intermediate2_dervs))
            intermediate2_dervs=zero
           endif
          endif
          
          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           allocate(yl_dervs((lb+1)**2,3,3),stat=cpksalloc(63))
           ASSERT(cpksalloc(63).eq.0)
           MEMLOG(size(yl_dervs))
           yl_dervs=zero
          endif

          cartesian_r2l: do i=1,3
             
             do i_ind2=1,n_indep_b
                   eq_atom_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%I_equal_atom
                   coeff_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%c  
                   magn_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%m
                cont1_r2l : do i_cont2=1,n_contributing_b(i_ind2) 

                   eq_b=eq_atom_b(i_cont2)
                   equal_r2l: do eq_a=1,unique_atoms(na)%n_equal_atoms

                      if ( (eq_a /= 1 .or. .not.moving_a) .and. &
                           (eq_atom_b(i_cont2) /= 1 .or. .not.moving_b) ) cycle

                      ! calculate gradients of solid harmonics
                      lm = lb**2 + magn_b(i_cont2)

                      if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)
                       yl_dervs=zero
                       do k2dr=1,3
                       do llm=1,(lb+1)**2
                         do i_sum=1,solhrules_differential(yl_map(k2dr),llm)%n_summands
                            yl_dervs(llm,i,k2dr) =  yl_dervs(llm,i,k2dr) + &
                            solhrules_differential(yl_map(k2dr),llm)%coef(i_sum)* &
                             yl_grad(eq_a,eq_b,solhrules_differential(yl_map(k2dr),llm)%lm_sh(i_sum),i)
                         enddo
                       enddo 
                       enddo 
FPP_TIMER_STOP(t_calc_2c_dervs)
                      endif


                      xa = unique_atoms(na)%position(:,eq_a)
                      xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont2))
                      arg = sum((xa-xb)**2)

                     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
                      gamma_help(:,1:lb+4) = gamma(lb+4,fact1/fact0*arg)
                     else
                      gamma_help(:,1:lb+3) = gamma(lb+3,fact1/fact0*arg)
                     endif

                      if (moving_a .and. eq_a == 1) then
                         intermediate1 = intermediate1 + ( &                       ! (1)
                              (three/two-fact2*real(lb-1,kind=r8_kind))*( &
                              yl_grad(eq_a,eq_b,lm,i)*gamma_help(:,lb+1)  &
                                        - yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              gamma_help(:,lb+2)*two*fact1/fact0*(xa(i)-xb(i))))*coeff_b(i_cont2)
                         intermediate1 = intermediate1 + ( &                       ! (2)
                              fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                              gamma_help(:,lb+2))* coeff_b(i_cont2)
                         intermediate1 = intermediate1 + (&                        ! (3)
                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*(&
                              two*(xa(i)-xb(i))*gamma_help(:,lb+2) - &
                              arg*gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))))*&
                              coeff_b(i_cont2)

                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                         intermediate1_dervs(:,i,:) =  &  ! (1.1)
                          spread(yl_dervs(lm,i,:),1,num)* &
                          spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                                  gamma_help(:,lb+1)*coeff_b(i_cont2), 2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) +  &  ! (1,2)
                              spread( two*(xa-xb),1,num)* &
                              spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              yl_grad(eq_a,eq_b,lm,i)*(-gamma_help(:,lb+2))*coeff_b(i_cont2)* &
                              fact1/fact0,2,3 )

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) -  &  ! (1.3)
                          spread( yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
                          spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              gamma_help(:,lb+2)*two*fact1/fact0*(xa(i)-xb(i))*coeff_b(i_cont2),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:)-  &   ! (1.4)
                          spread(two*(xa-xb),1,num)* &
                          spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                                  yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                                  (-gamma_help(:,lb+3))*two*fact1/fact0*(xa(i)-xb(i))*coeff_b(i_cont2)* &
                                  fact1/fact0,2,3 )

                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) -  &  ! (1.5)
                              (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              gamma_help(:,lb+2)*two*fact1/fact0*coeff_b(i_cont2)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) +  &  ! (2.1)
                          spread(yl_dervs(lm,i,:),1,num)* &
                          spread( fact1*fact2/fact0*arg*&
                                  gamma_help(:,lb+2)* coeff_b(i_cont2),2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) +  &  ! (2.2)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*&
                                  gamma_help(:,lb+2)* coeff_b(i_cont2),2,3 )

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) +  &  ! (2.3)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                                  (-gamma_help(:,lb+3))* coeff_b(i_cont2)* &
                                  fact1/fact0,2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &   ! (3.1)
                          spread( yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
                          spread( fact1*fact2/fact0*two*(xa(i)-xb(i))* &
                                  gamma_help(:,lb+2)*coeff_b(i_cont2),2,3 )

                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) + &   ! (3.2)
                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              two*gamma_help(:,lb+2)*coeff_b(i_cont2)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) + &   ! (3.3)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                                  two*(xa(i)-xb(i))*(-gamma_help(:,lb+3))*coeff_b(i_cont2)*&
                                  fact1/fact0,2,3)

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.4) yl
                          spread(yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
                          spread( fact1*fact2/fact0*&
                                  arg*gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))*&
                                  coeff_b(i_cont2),2,3 )

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.5)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                                  gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))*&
                                  coeff_b(i_cont2),2,3 )

                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.6)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                                  arg*(-gamma_help(:,lb+4))*two*fact1/fact0*(xa(i)-xb(i))*&
                                  coeff_b(i_cont2)*fact1/fact0,2,3 )

                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) - &   ! (3.7)
                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              arg*gamma_help(:,lb+3)*two*fact1/fact0*&
                              coeff_b(i_cont2)

                        do k2dr=1,3
                        fit_int_dervs_a(:,:,1,i_ind2,i,k2dr,eq_b) = &
                         fit_int_dervs_a(:,:,1,i_ind2,i,k2dr,eq_b) + (-one)**lb* &
                            unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**lb*&
                            one/(fact1*fact0*sqrt(fact0))*intermediate1_dervs(:,i,k2dr),&
                            cutoff,zero)
                        enddo

                           intermediate1_dervs = zero

FPP_TIMER_STOP(t_calc_2c_dervs)
                       endif

                      endif
                      if (moving_b .and. eq_atom_b(i_cont2) == 1) then
                         intermediate2 = intermediate2 + ( &                       ! (1)
                              (three/two-fact2*real(lb-1,kind=r8_kind))*( &
                              -yl_grad(eq_a,eq_b,lm,i)*gamma_help(:,lb+1) &
                              + yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              gamma_help(:,lb+2)*two*fact1/fact0*(xa(i)-xb(i))))*coeff_b(i_cont2)
                         intermediate2 = intermediate2 - ( &                       ! (2)
                              fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                              gamma_help(:,lb+2))* coeff_b(i_cont2)
                         intermediate2 = intermediate2 + ( &                       ! (3)
                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*( &
                              -two*(xa(i)-xb(i))*gamma_help(:,lb+2) + &
                              arg*gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))))*&
                              coeff_b(i_cont2)
                       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                         intermediate2_dervs(:,i,:) =  &      ! (1.1)
                         +spread( yl_dervs(lm,i,:),1,num)*&
                          spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              gamma_help(:,lb+1)*coeff_b(i_cont2),2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) +  &      ! (1.2)
                           spread( two*(xa-xb),1,num)*&
                           spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              yl_grad(eq_a,eq_b,lm,i)*(-gamma_help(:,lb+2))*coeff_b(i_cont2)*&
                              fact1/fact0,2,3)


                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &      ! (1.3)
                           spread( yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)*&
                           spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                              gamma_help(:,lb+2)*two*fact1/fact0*(xa(i)-xb(i))*coeff_b(i_cont2),2,3)

                         intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:) -  &      ! (1.4)
                           spread( two*(xa-xb),1,num)*&
                           spread( (three/two-fact2*real(lb-1,kind=r8_kind))* &
                               yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              (-gamma_help(:,lb+3))*two*fact1/fact0*(xa(i)-xb(i))*coeff_b(i_cont2)*&
                               fact1/fact0,2,3)

                         intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) -  &      ! (1.5)
                              (three/two-fact2*real(lb-1,kind=r8_kind))* &
                               yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
                              gamma_help(:,lb+2)*two*fact1/fact0*coeff_b(i_cont2)

!!                         intermediate2 = intermediate2 - ( &                       ! (2)
!!                              fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
!!                              gamma_help(:,lb+2))* coeff_b(i_cont2)
                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)+&   !(2.1)
                           spread( yl_dervs(lm,i,:),1,num)* &
                           spread( fact1*fact2/fact0*arg*gamma_help(:,lb+2)* &
                           coeff_b(i_cont2),2,3)

                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)+&   !(2.2)
                          spread( two*(xa-xb),1,num)*&
                          spread(fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)* &
                                 gamma_help(:,lb+2)* coeff_b(i_cont2),2,3 )
                          
                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)+&   !(2.3)
                          spread( two*(xa-xb),1,num)*&
                          spread(fact1*fact2/fact0*yl_grad(eq_a,eq_b,lm,i)*arg*&
                                 (-gamma_help(:,lb+3))* coeff_b(i_cont2)*fact1/fact0,2,3)

!!                         intermediate2 = intermediate2 + ( &                       ! (3)
!!                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*( &
!!                              -two*(xa(i)-xb(i))*gamma_help(:,lb+2) + &
!!                              arg*gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i)) ) )*&
!!                              coeff_b(i_cont2)

                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)+&   !(3.1)
                           spread( yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
                           spread(fact1*fact2/fact0*two*(xa(i)-xb(i))* &
                                  gamma_help(:,lb+2)*coeff_b(i_cont2),2,3)
                          
                          intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i)+&   !(3.2)
                          fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)* &
                          two*gamma_help(:,lb+2)*coeff_b(i_cont2)
                          
                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)+&   !(3.3)
                          spread( two*(xa-xb),1,num)*  &
                          spread(fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)* &
                          two*(xa(i)-xb(i))*(-gamma_help(:,lb+3))* &
                                  coeff_b(i_cont2)*fact1/fact0,2,3)

     
                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)-&   !(3.4)
                          spread( yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
                          spread( fact1*fact2/fact0*arg*gamma_help(:,lb+3)* &
                                  two*fact1/fact0*(xa(i)-xb(i))*coeff_b(i_cont2),2,3)

                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)-&   !(3.5)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)* &
                                  gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))* &
                                  coeff_b(i_cont2),2,3)

                          intermediate2_dervs(:,i,:) = intermediate2_dervs(:,i,:)-&   !(3.6)
                          spread( two*(xa-xb),1,num)* &
                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)* &
                                  arg*(-gamma_help(:,lb+4))*two*fact1/fact0*(xa(i)-xb(i))* &
                                  coeff_b(i_cont2)*fact1/fact0,2,3)

                          intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i)-&   !(3.7)
                          fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)* &
                          arg*gamma_help(:,lb+3)*two*fact1/fact0*coeff_b(i_cont2)




!                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.4) yl
!                          spread(yl_grad(eq_a,eq_atom_b(i_cont2),lm,:),1,num)* &
!                          spread( fact1*fact2/fact0*&
!                                  arg*gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))*&
!                                  coeff_b(i_cont2),2,3 )

!                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.5)
!                          spread( two*(xa-xb),1,num)* &
!                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
!                                  gamma_help(:,lb+3)*two*fact1/fact0*(xa(i)-xb(i))*&
!                                  coeff_b(i_cont2),2,3 )

!                         intermediate1_dervs(:,i,:) = intermediate1_dervs(:,i,:) - &   ! (3.6)
!                          spread( two*(xa-xb),1,num)* &
!                          spread( fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
!                                  arg*(-gamma_help(:,lb+4))*two*fact1/fact0*(xa(i)-xb(i))*&
!                                  coeff_b(i_cont2)*fact1/fact0,2,3 )

!                         intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) - &   ! (3.7)
!                              fact1*fact2/fact0*yl_arr(eq_a,eq_atom_b(i_cont2),lm)*&
!                              arg*gamma_help(:,lb+3)*two*fact1/fact0*&
!                              coeff_b(i_cont2)

                         do k2dr=1,3
                          fit_int_dervs_b(:,:,1,i_ind2,i,k2dr,eq_a) = &
                              fit_int_dervs_b(:,:,1,i_ind2,i,k2dr,eq_a) + (-one)**lb* &
                               unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**lb*&
                               one/(fact1*fact0*sqrt(fact0))*intermediate2_dervs(:,i,k2dr),&
                               cutoff,zero) 
                         enddo

                         intermediate2_dervs=zero

FPP_TIMER_STOP(t_calc_2c_dervs)

                       endif
                      endif

                   enddo equal_r2l
                enddo cont1_r2l
                if (moving_a) then

                fit_int_grad_a(:,:,1,i_ind2,i) = fit_int_grad_a(:,:,1,i_ind2,i) + (-one)**lb* &
                     unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**lb*&
                     one/(fact1*fact0*sqrt(fact0))*intermediate1,&
                     cutoff,zero)
                   intermediate1 = zero
                endif
                !if(na.ne.nb) then
                   if (moving_b) then
                   fit_int_grad_b(:,:,1,i_ind2,i) = fit_int_grad_b(:,:,1,i_ind2,i) + (-one)**lb* &
                        unpack(two*pi*pi*sqrt(pi)*(-two*fact1/fact0)**lb*&
                        one/(fact1*fact0*sqrt(fact0))*intermediate2,&
                        cutoff,zero) 
                   intermediate2 = zero
                endif
                !endif
             enddo
                  
          enddo cartesian_r2l
!         print*,'lr2l fit_int_grad_a',na,nb
!         print*,sum(fit_int_grad_a(:,:,:,:,1)),sum(fit_int_dervs_a(:,:,:,:,1,1,1))
!         print*,'lr2l fit_int_grad_b',na,nb
!         print*,sum(fit_int_grad_b(:,:,:,:,1)),sum(fit_int_dervs_b(:,:,:,:,1,1,1))
!          print*,'intermediate1_dervs r2l'
!          do k2dr=1,3
!           print*,sum(intermediate1_dervs(:,1,k2dr)), &
!                    sum(intermediate1_dervs(:,2,k2dr)), & 
!                       sum(intermediate1_dervs(:,3,k2dr))
!          enddo

          if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
           MEMLOG(-size(yl_dervs))
           deallocate(yl_dervs,stat=cpksalloc(63))
           ASSERT(cpksalloc(63).eq.0)
           cpksalloc(63)=1
          endif

          if (moving_a) then
             deallocate(intermediate1,STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             MEMLOG(-size(intermediate1_dervs))
             deallocate(intermediate1_dervs,STAT=cpksalloc(72))
             ASSERT(cpksalloc(72).eq.0)
             cpksalloc(72)=1
           endif
          endif

          if (moving_b) then
             deallocate(intermediate2,STAT=alloc_stat)
             ASSERT(alloc_stat.eq.0)
           if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
             MEMLOG(-size(intermediate2_dervs))
             deallocate(intermediate2_dervs,STAT=cpksalloc(73))
             ASSERT(cpksalloc(73).eq.0)
             cpksalloc(73)=1
           endif
          endif
          
       endif
    endif

    if (moving_a) then

!   if(.not.integralpar_cpks_contribs) then
!   print*,'fit_int_dervs_a lr2',ma,shape(fit_int_grad_a),'shape fit_int_grad_a'
!   print*,sum(fit_int_grad_a(:,:,1,:,1)),sum(fit_int_dervs_a(:,:,1,:,1,1,1))
!   endif
!
        call fitcontract_2c_grad(fit_int_grad_a,ma,1) !lr2
       MEMLOG(-size(fit_int_grad_a))
       deallocate( fit_int_grad_a, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!    print*,'fit_int_dervs_a lr2'
!     do k2dr=1,3
!     print*,sum( fit_int_dervs_a(:,:,1,:,1,k2dr,1) ), &
!               sum( fit_int_dervs_a(:,:,1,:,2,k2dr,1) ), &
!                   sum( fit_int_dervs_a(:,:,1,:,3,k2dr,1) )
!     enddo

FPP_TIMER_START(t_contract_2c_dervs)
         call fitcontract_2c_dervs(fit_int_dervs_a,ma,mb) !lr2(1)
FPP_TIMER_STOP(t_contract_2c_dervs)
        MEMLOG(-size(fit_int_dervs_a))
        deallocate( fit_int_dervs_a, STAT=cpksalloc(60))
        ASSERT(cpksalloc(60).eq.0)
        cpksalloc(60)=1
       endif

    endif

    if (moving_b) then
!   if(.not.integralpar_cpks_contribs) then
!   print*,'fit_int_dervs_b lr2'
!   print*,sum(fit_int_grad_b(:,:,1,1,1)),sum(fit_int_dervs_b(:,:,1,1,1,1,1))
!   endif
        call fitcontract_2c_grad(fit_int_grad_b,mb,2) !lr2
       MEMLOG(-size(fit_int_grad_b))
       deallocate( fit_int_grad_b, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then

!     print*,'fit_int_dervs_b lr2'
!     do k2dr=1,3
!     print*,sum( fit_int_dervs_b(:,:,1,:,1,k2dr,1) ), &
!               sum( fit_int_dervs_b(:,:,1,:,2,k2dr,1) ), &
!                   sum( fit_int_dervs_b(:,:,1,:,3,k2dr,1) )
!     enddo

FPP_TIMER_START(t_contract_2c_dervs)
        call fitcontract_2c_dervs(fit_int_dervs_b,mb,ma) !lr2 (1)
FPP_TIMER_STOP(t_contract_2c_dervs)
       MEMLOG(-size(fit_int_dervs_b))
       deallocate( fit_int_dervs_b, STAT=cpksalloc(67))
       ASSERT(cpksalloc(67).eq.0)
       cpksalloc(67)=1
     endif
    endif

    if (l_max.gt.0) then
       deallocate(yl_arr,yl_grad,STAT=cpksalloc(65))
       ASSERT(cpksalloc(65).eq.0)
       cpksalloc(65)=1
    endif

    deallocate(  n_contributing_a, n_contributing_b, gamma_help, &
         STAT=alloc_stat)

    ASSERT(alloc_stat.eq.0)

    
    call dump_exponent_data()

  end subroutine lr2_calc_grad


  subroutine r2_calc_grad(na,nb)
   use cpksdervs_matrices
   use calc3c_switches

    ! Purpose: calculate primitives for both angular momenta
    !          of r2-type (la=lb=-1)
    !          Symmetry adaption included.
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help arrays
    real(kind=r8_kind),allocatable   :: &
         intermediate1(:),intermediate1_dervs(:,:,:), &
         intermediate2(:),intermediate2_dervs(:,:,:), &
         gamma_help(:,:),a1(:),a2(:),a3(:)
    ! symmetry adapted integrals
    real(kind=r8_kind),allocatable   :: &
     fit_int_grad_a(:,:,:,:,:), fit_int_dervs_a(:,:,:,:,:,:,:),&
     fit_int_grad_b(:,:,:,:,:), fit_int_dervs_b(:,:,:,:,:,:,:)
    integer(kind=i4_kind)            :: eq_a,eq_b,n_equal_a,n_equal_b, &
         alloc_stat,i,k2dr
    integer(kind=i4_kind) :: ma, mb
    logical :: moving_a,moving_b
    !------------ Executable code --------------------------------

    ma = unique_atoms(na)%moving_atom
    mb = unique_atoms(nb)%moving_atom
    moving_a = ma > 0
    moving_b = mb > 0
    if (.not.moving_a .and. .not.moving_b) return
    call get_exponent_data(na,-1,nb,-1)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms

    ! To avoid doublication of code double DO-loops over all pairs of equal
    ! atoms are used, though only the following pairs are final processed

    !   a) if moving_a : the pairs (1,eq_b) and
    !   b) if moving_b : the pairs (eq_a,1) .

    if (moving_a) then
       allocate(fit_int_grad_a(naexps,nbexps,1,1,3),intermediate1(num), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_a= zero
       intermediate1 = zero

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
         allocate(fit_int_dervs_a(naexps,nbexps,1,1,3,3,n_equal_b), &
                  intermediate1_dervs(num,3,3), STAT=cpksalloc(58)) 
         ASSERT(cpksalloc(58).eq.0)
         MEMLOG(size(fit_int_dervs_a)+size(intermediate1_dervs))
       fit_int_dervs_a=zero
       intermediate1_dervs=zero
       endif
    endif

    if (moving_b) then
       allocate(fit_int_grad_b(naexps,nbexps,1,1,3),intermediate2(num), &
            STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       fit_int_grad_b= zero
       intermediate2 = zero

       if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
         allocate(fit_int_dervs_b(naexps,nbexps,1,1,3,3,n_equal_a), &
                  intermediate2_dervs(num,3,3), STAT=cpksalloc(64)) 
         ASSERT(cpksalloc(64).eq.0)
       MEMLOG(size(fit_int_dervs_b)+size(intermediate2_dervs))
       fit_int_dervs_b=zero
       intermediate2_dervs=zero
       endif
    endif

    if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
    allocate(gamma_help(num,5),a1(num),a2(num),a3(num), &
         STAT=cpksalloc(59))
    else
    allocate(gamma_help(num,4),a1(num),a2(num),a3(num), &
         STAT=cpksalloc(59))
    endif
     ASSERT(cpksalloc(59).eq.0)

    gamma_help = zero
    a1 = zero
    a2 = zero
    a3 = zero


    cartesian: do i=1,3
       do eq_a = 1,n_equal_a
          do eq_b = 1,n_equal_b

             if ( (eq_a /= 1 .or. .not.moving_a) .and. &
                  (eq_b /= 1 .or. .not.moving_b) ) cycle

             xa=unique_atoms(na)%position(:,eq_a)
             xb=unique_atoms(nb)%position(:,eq_b)

             arg=sum((xa-xb)**2)

             if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
              gamma_help(:,1:5) = gamma(5,fact1/fact0*arg)
             else
              gamma_help(:,1:4) = gamma(4,fact1/fact0*arg)
             endif

             a1 = three/(two*fact1)+three/(four*fact0**2)
             a2 = three/(two*fact0)-three*fact1/(fact0*fact0**2)
             a3 = fact1**2/fact0**4

             if (moving_a .and. eq_a == 1) then

                intermediate1 =  intermediate1 + &
                 two*(xa(i)-xb(i))*gamma_help(:,2)*(a2-fact1/fact0*a1) 

                intermediate1 = intermediate1 + &
                 two*arg*(xa(i)-xb(i))*gamma_help(:,3)*(two*a3-fact1/fact0*a2)

                intermediate1 = intermediate1 - &
                 a3*arg**2*gamma_help(:,4)*two*fact1/fact0*(xa(i)-xb(i))

             if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)


              intermediate1_dervs(:,:,i)= & !!! intermediate1_dervs(:,:,i)  &      ! (1,2)
                 +spread(two*(xa(i)-xb(i))* &
                    (-gamma_help(:,3))*(a2-fact1/fact0*a1)*fact1/fact0,2,3)* &
                         spread(two*(xa-xb),1,num)

              intermediate1_dervs(:,i,i)=intermediate1_dervs(:,i,i) &       ! (1.1)
             +two*gamma_help(:,2)*(a2-fact1/fact0*a1)

                intermediate1_dervs(:,:,i) = intermediate1_dervs(:,:,i) + &           ! (2.1)
                 spread(four*(xa-xb),1,num)* &
                 spread((xa(i)-xb(i))*gamma_help(:,3)*(two*a3-fact1/fact0*a2),2,3)

                intermediate1_dervs(:,:,i) = intermediate1_dervs(:,:,i) + &           ! (2.2)
                 spread(two*arg*(xa(i)-xb(i))*&
                        (-gamma_help(:,4))*(two*a3-fact1/fact0*a2)* &
                        fact1/fact0,2,3) * spread(two*(xa-xb),1,num)

                intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) + & ! (2.3)
                 two*arg*gamma_help(:,3)*(two*a3-fact1/fact0*a2)

                intermediate1_dervs(:,:,i) = intermediate1_dervs(:,:,i) - &
                 spread(two*(xa-xb),1,num)* &
                  spread(a3*two*arg*gamma_help(:,4)*two*fact1/fact0*(xa(i)-xb(i)),2,3)   ! (3.1)

                intermediate1_dervs(:,:,i) = intermediate1_dervs(:,:,i) - &
                 spread(two*(xa-xb),1,num)*spread( &
                 a3*arg**2*(-gamma_help(:,5))*two*fact1/fact0*(xa(i)-xb(i))* &
                 fact1/fact0,2,3)

                intermediate1_dervs(:,i,i) = intermediate1_dervs(:,i,i) - &
                 a3*arg**2*gamma_help(:,4)*two*fact1/fact0

                 do k2dr=1,3
                 fit_int_dervs_a(:,:,1,1,i,k2dr,eq_b) = &
                   two*pi*pi*sqrt(pi)* &
                   unpack(intermediate1_dervs(:,k2dr,i)/(fact1*sqrt(fact0)),cutoff,zero)
                 enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
            endif

             endif

             if (moving_b .and. eq_b == 1) then

                intermediate2 =  intermediate2 - &
                     two*(xa(i)-xb(i))*gamma_help(:,2)*(a2-fact1/fact0*a1) 

                intermediate2 = intermediate2 - &
                     two*arg*(xa(i)-xb(i))*gamma_help(:,3)*&
                     (two*a3-fact1/fact0*a2)

                intermediate2 = intermediate2 + &
                     a3*arg**2*gamma_help(:,4)*two*fact1/fact0*&
                     (xa(i)-xb(i))
              if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
FPP_TIMER_START(t_calc_2c_dervs)

                intermediate2_dervs(:,:,i) = & !!!  intermediate2_dervs(:,:,i) + &
                 spread(two*(xa-xb),1,num)* &
                 spread( two*(xa(i)-xb(i))*(-gamma_help(:,3))*(a2-fact1/fact0*a1) * &
                         fact1/fact0,2,3)

                intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) + &
                     two*gamma_help(:,2)*(a2-fact1/fact0*a1) 


                intermediate2_dervs(:,:,i) = intermediate2_dervs(:,:,i) + &
                   spread(two*(xa-xb),1,num)* &
                   spread( two*(xa(i)-xb(i))*gamma_help(:,3)*&
                     (two*a3-fact1/fact0*a2),2,3 )

                intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) + &
                     two*arg*gamma_help(:,3)*(two*a3-fact1/fact0*a2)

                intermediate2_dervs(:,:,i) = intermediate2_dervs(:,:,i) + &
                  spread(two*(xa-xb),1,num)* &
                  spread( two*arg*(xa(i)-xb(i))*(-gamma_help(:,4))*&
                          (two*a3-fact1/fact0*a2)*fact1/fact0,2,3 )

                intermediate2_dervs(:,:,i) = intermediate2_dervs(:,:,i) - &
                  spread( two*(xa-xb),1,num)* &
                  spread( a3*arg*two*gamma_help(:,4)*two*fact1/fact0*&
                     (xa(i)-xb(i)),2,3 )

                intermediate2_dervs(:,:,i) = intermediate2_dervs(:,:,i) - &
                  spread( two*(xa-xb),1,num)* &
                  spread( a3*arg**2*(-gamma_help(:,5))*two*fact1/fact0*&
                     (xa(i)-xb(i))*fact1/fact0,2,3 )

                intermediate2_dervs(:,i,i) = intermediate2_dervs(:,i,i) - &
                     a3*arg**2*gamma_help(:,4)*two*fact1/fact0

                do k2dr=1,3
                fit_int_dervs_b(:,:,1,1,i,k2dr,eq_a) = &
                  two*pi*pi*sqrt(pi)* &
                   unpack(intermediate2_dervs(:,k2dr,i)/(fact1*sqrt(fact0)),cutoff,zero)
                enddo

FPP_TIMER_STOP(t_calc_2c_dervs)
              endif
             endif

          enddo
       enddo

       if (moving_a) then
       fit_int_grad_a(:,:,1,1,i) = &
           two*pi*pi*sqrt(pi)*unpack(intermediate1/(fact1*sqrt(fact0)),cutoff,zero)
          intermediate1 = zero
       endif


       !if (na.ne.nb) then
       if (moving_b) then
       fit_int_grad_b(:,:,1,1,i) =  &
            two*pi*pi*sqrt(pi)*unpack(intermediate2/(fact1*sqrt(fact0)),cutoff,zero)
          intermediate2 = zero
       endif
       !endif

    enddo cartesian

    if (moving_a) then

        call fitcontract_2c_grad(fit_int_grad_a,ma,1) !r2
       MEMLOG(-size(fit_int_grad_a)-size(intermediate1))
       deallocate(fit_int_grad_a,intermediate1,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!    print*,'fit_int_dervs_a  r2'
!    if(size(fit_int_dervs_a,7).gt.1) then
!    do k2dr=1,3
!    print*,sum(fit_int_dervs_a(:,:,1,1,1,k2dr,1)), &
!             sum(fit_int_dervs_a(:,:,1,1,2,k2dr,1)),  &
!              sum(fit_int_dervs_a(:,:,1,1,3,k2dr,1)), &
!               sum(fit_int_dervs_a(:,:,1,1,1,k2dr,2)),  &
!                 sum(fit_int_dervs_a(:,:,1,1,2,k2dr,2)), &
!                   sum(fit_int_dervs_a(:,:,1,1,3,k2dr,2))
!    enddo
!    else
!    do k2dr=1,3
!    print*,sum(fit_int_dervs_a(:,:,1,1,1,k2dr,1)),  &
!              sum(fit_int_dervs_a(:,:,1,1,2,k2dr,1)),  &
!                 sum(fit_int_dervs_a(:,:,1,1,3,k2dr,1))
!    enddo
!    endif

FPP_TIMER_START(t_contract_2c_dervs)
        call fitcontract_2c_dervs(fit_int_dervs_a,ma,mb) ! r2(1)
FPP_TIMER_STOP(t_contract_2c_dervs)
      MEMLOG(-size(fit_int_dervs_a)-size(intermediate1_dervs))
      deallocate(fit_int_dervs_a,intermediate1_dervs,STAT=cpksalloc(58))
      ASSERT(cpksalloc(58).eq.0)
      cpksalloc(58)=1
     endif

    endif

    if (moving_b) then
        call fitcontract_2c_grad(fit_int_grad_b,mb,2)  !r2
       MEMLOG(-size(fit_int_grad_b)-size(intermediate2))
       deallocate(fit_int_grad_b,intermediate2,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)

     if(.not.integralpar_cpks_contribs.and.integralpar_2dervs) then
!    print*,'fit_int_dervs_b  r2'
!   if(size(fit_int_dervs_b,7).gt.1) then
!    do k2dr=1,3
!    print*,sum(fit_int_dervs_b(:,:,1,1,1,k2dr,1)),sum(fit_int_dervs_b(:,:,1,1,2,k2dr,1)),  &
!              sum(fit_int_dervs_b(:,:,1,1,3,k2dr,1)),sum(fit_int_dervs_b(:,:,1,1,1,k2dr,2)),  &
!                 sum(fit_int_dervs_b(:,:,1,1,2,k2dr,2)),sum(fit_int_dervs_b(:,:,1,1,3,k2dr,2))
!    enddo
!   else
!    do k2dr=1,3
!    print*,sum(fit_int_dervs_b(:,:,1,1,1,k2dr,1)),  &
!              sum(fit_int_dervs_b(:,:,1,1,2,k2dr,1)),  &
!                 sum(fit_int_dervs_b(:,:,1,1,3,k2dr,1))
!    enddo
!   endif

FPP_TIMER_START(t_contract_2c_dervs)
        call fitcontract_2c_dervs(fit_int_dervs_b,mb,ma) !r2 (1)
FPP_TIMER_STOP(t_contract_2c_dervs)
      MEMLOG(-size(fit_int_dervs_b)-size(intermediate2_dervs))
      deallocate(fit_int_dervs_b,intermediate2_dervs,STAT=cpksalloc(64))
      ASSERT(cpksalloc(64).eq.0)
      cpksalloc(64)=1
     endif
    endif

    deallocate(gamma_help,a1,a2,a3,STAT=cpksalloc(59))
    ASSERT(cpksalloc(59).eq.0)
    cpksalloc(59)=1

    call dump_exponent_data()

  end subroutine r2_calc_grad
    

  subroutine get_exponent_data(na, la, nb, lb)
    !
    !  Purpose: allocate space for the exponent arrays, the help arrays
    !           and get the numbers from the unique_atom_module.
    !
    !           Global module variables set in this routine:
    !
    !           naexps, nbexps
    !           fact0, fact1, fact2
    !           cutoff, num, i_ir
    !
    !  Subroutine called by: all xx_calc routines in this module
    !  Author: FN
    !  Date: 7/96
    !------------ Modules used --------------------------------------
    use type_module ! type specification parameters
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(kind=i4_kind), intent(in) :: na,la,nb,lb
    !** End of interface *****************************************
    !------------ Declaration of local variables --------------------
    real(kind=r8_kind),allocatable  :: fact0_arr(:,:),fact1_arr(:,:),&
         fact2_arr(:,:)
    real(kind=r8_kind), allocatable :: aexps(:), bexps(:)
    integer(kind=i4_kind)   :: alloc_stat

    if ( la .eq. -1 ) then
       naexps = unique_atoms(na)%r2_ch%n_exponents ! global

       allocate(aexps(naexps),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       aexps = unique_atoms(na)%r2_ch%exponents(:)
    else
       naexps = unique_atoms(na)%l_ch(la)%n_exponents ! global

       allocate(aexps(naexps),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       aexps = unique_atoms(na)%l_ch(la)%exponents(:)
    endif

    if ( lb .eq. -1 ) then
       nbexps = unique_atoms(nb)%r2_ch%n_exponents ! global

       allocate(bexps(nbexps),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       bexps = unique_atoms(nb)%r2_ch%exponents(:)
    else
       nbexps = unique_atoms(nb)%l_ch(lb)%n_exponents ! global

       allocate(bexps(nbexps),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       bexps = unique_atoms(nb)%l_ch(lb)%exponents(:)
    endif

    i_ir = get_totalsymmetric_irrep() ! global

    allocate(fact0_arr(naexps,nbexps),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(fact1_arr(naexps,nbexps),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(cutoff(naexps,nbexps),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    ! There is no cutoff criterion for the charge overlap
    ! matrix but in order to be able to use the 'prod_rule',
    ! 'diff_rule' etc. we pack th exponents to a linear array
    cutoff = .true.
    num=count(cutoff)

    allocate (fact0(num), fact1(num), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    fact0_arr=(spread(aexps,2,nbexps)+spread(bexps,1,naexps))
    fact1_arr=(spread(aexps,2,nbexps)*spread(bexps,1,naexps))

    if ( la .eq. -1 .or. lb .eq. -1 ) then

       allocate(fact2_arr(naexps,nbexps),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       allocate (fact2(num),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       if(lb.eq.-1) then
          fact2_arr=(spread(aexps,2,nbexps)/spread(bexps,1,naexps))
       else
          fact2_arr=(spread(bexps,1,naexps)/spread(aexps,2,nbexps))
       endif

       ! FIXME: segfaults here for naxeps x nbexps being zero:
       if ( naexps * nbexps > 0 ) then
          !
          ! This is global:
          !
          fact2 = pack(fact2_arr, cutoff)
       endif

       deallocate(fact2_arr,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    endif

    if ( naexps * nbexps > 0 ) then
       !
       ! These are global module variables:
       !
       fact0 = pack(fact0_arr, cutoff)
       fact1 = pack(fact1_arr, cutoff)
    endif

    deallocate(fact0_arr,fact1_arr,aexps,bexps,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
  end subroutine get_exponent_data


  subroutine dump_exponent_data()
    ! Purpose : deallocate all allocatable data set in the routine
    !           'get_exponent_data'.
    ! Routine called by: all xx_calc  routines in this module
    !** End of interface *****************************************
    integer(kind=i4_kind)   :: alloc_stat
    ! --------------- Executable code --------------------------
    naexps = 0
    nbexps = 0
    i_ir = 0
    num = 0

    deallocate(fact0,fact1,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (1) failed")
    if(allocated(fact2)) then
       deallocate(fact2,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("dump_exponent_data: deallocate (3) failed")
    endif
    deallocate(cutoff,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (4) failed")
  end subroutine dump_exponent_data

  !*************************************************************

end module gradient_2c_fit_ch_module

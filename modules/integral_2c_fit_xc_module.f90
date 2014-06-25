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
module  integral_2c_fit_xc_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains all routines necessary to perform
  !           the two center fitfct integrals:
  !           - <fk|fl>
  !
  !  Module called by: integral_calc_quad_2cff.f90
  !
  !
  !  References: ...
  !  Author: FN
  !  Date: 8/96
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing) 
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: subroutine extended to the treatment of mixed
  !              overlap integrals <f_k|g_l> and <g_k|f_l>. To
  !              this end an additional argument "INT_TYPE" has been
  !              introduced for the subroutine XC_OVERLAP
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  use type_module ! type specification parameters
  use symmetry_data_module, only : get_totalsymmetric_irrep
  use unique_atom_module
  use solid_harmonics_module, only : solid_harmonics_scalar
  use solhrules_module, only: prod_rule, diff_rule
  use fitcontract_2c_module, only: fitcontract_2c
  use output_module, only: output_int_fitcontract,output_int_2c_fit
  use iounitadmin_module, only: write_to_output_units 
  use time_module, only: stop_timer, start_timer
  use timer_module, only: timer_int_prim_2cff, timer_int_cont_2cff
  use integralpar_module, only: integralpar_i_int_part

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !------------ public functions and subroutines ------------------
  public xc_overlap


  !================================================================
  ! End of public interface of module
  !================================================================
  ! constants
  real(kind=r8_kind),parameter    :: expmax=50.0_r8_kind
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: two=2.0_r8_kind,three=3.0_r8_kind, &
       one=1.0_r8_kind,four=4.0_r8_kind,five=5.0_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind

  ! variables needed by all subroutines: these are set in the
  ! routine get_exponent_data(na,la,nb,lb)
  integer(kind=i4_kind)           :: naexps,nbexps
  real(kind=r8_kind),allocatable  :: aexps(:),bexps(:)
  real(kind=r8_kind),allocatable  :: fact0(:),fact1(:),fact2(:)
  logical,allocatable             :: cutoff(:,:)
  integer(kind=i4_kind)           :: num,i_ir
  character(len=5)                :: flag ! "xc_xc", "xc_ch", or "ch_xc"

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************

  subroutine xc_overlap(int_type)
    !  Purpose: calculate <g_k|g_l>, <g_k|f_l>, or <f_k|g_l>
    !** End of interface *****************************************
    !------------ Modules used --------------------------------
    use int_data_2cff_module
    implicit none
    !------------ Declaration of formal parameter ----------------
    character(len=5) :: int_type ! "xc_xc", "xc_ch", or "ch_xc"
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
    flag=int_type

    if ( (la.gt.0.and.lb.ge.0).or.(la.ge.0.and.lb.gt.0)) then
       if (output_int_2c_fit) call write_to_output_units &
            ("XC_OVERLAP: calling ll_calc")
       call ll_calc(na,la,nb,lb)
    elseif (la.eq.0.and.lb.eq.0) then
       if (output_int_2c_fit) call write_to_output_units &
            ("XC_OVERLAP: calling ss_calc")
       call ss_calc(na,nb)
    elseif( (la.eq.-1.and.lb.ge.0).or.(lb.eq.-1.and.la.ge.0)) then
       if (output_int_2c_fit) call write_to_output_units &
            ("XC_OVERLAP: calling lr2_calc")
       call lr2_calc(na,la,nb,lb)
    elseif (la.eq.-1.and.lb.eq.-1) then
       if (output_int_2c_fit) call write_to_output_units &
            ("XC_OVERLAP: calling r2_calc")
       call r2_calc(na,nb)
    else
       call error_handler &
            ("XC_OVERLAP : something fishy")
    endif
  end subroutine xc_overlap

  !*************************************************************

  subroutine ss_calc(na,nb)
    ! Purpose: calculate primitives for both angular momenta
    !          equal zero.
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: intermediate(:), &
         fit_int(:,:,:,:)
    integer(kind=i4_kind) :: a_eq,b_eq,alloc_stat
    ! ------------- Executable code ---------------------------
    
    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,0,nb,0)
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: allocation (1) failed")
    allocate(intermediate(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: allocation (2) failed") 
    intermediate = zero
    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b 
    ! Loop over ALL equal atoms for symmetry adaption
    a_eq=1
    xa=unique_atoms(na)%position(:,a_eq)
    do b_eq=1,unique_atoms(nb)%n_equal_atoms
       xb=unique_atoms(nb)%position(:,b_eq)
       arg=sum((xa-xb)**2)
       intermediate = intermediate + exp(-fact1/fact0*arg)
    enddo

    intermediate = intermediate*pi*sqrt(pi)/(fact0*sqrt(fact0))* &
         unique_atoms(na)%n_equal_atoms
    fit_int(:,:,1,1) = unpack(intermediate,cutoff,zero)
    call dump_exponent_data()
    deallocate(intermediate,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("SS_CALC: deallocation (1) failed")

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do  LOCAL contraction
    if (output_int_fitcontract) call write_to_output_units &
         ("SS_CALC / XC : calling fitcontract_2c")   
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,flag)
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))
    
    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("SS_CALC: deallocation (2) failed")    
    
  end subroutine ss_calc

  !*************************************************************

  subroutine ll_calc(na,la,nb,lb)
    ! Purpose: calculate primitives for both angular momenta
    !          greater than zero.
    !          Symmetry adaption included.
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: la,lb
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help arrays for gamma-function
    real(kind=r8_kind),allocatable,dimension(:)   :: yl_arr
    ! intermediate arrays
    real(kind=r8_kind),allocatable,dimension(:,:)   :: intermediate1, &
         intermediate2,intermediate3,inter1
    real(kind=r8_kind),allocatable    :: fact(:)
    ! uncontracted but symmetryadapted integrals:
    real(kind=r8_kind),allocatable    :: fit_int(:,:,:,:)

    ! number of independents and contributing fcts
    integer(kind=i4_kind)   :: n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind),allocatable  :: n_contributing_a(:), &
         n_contributing_b(:)
    real(kind=r8_kind),pointer       :: coeff_a(:),coeff_b(:)
    integer(kind=i4_kind),pointer,dimension(:)   :: magn_a,&
         magn_b,eq_atom_a,eq_atom_b
    
    integer(kind=i4_kind)   :: counter,i_l,i_m, &
         alloc_stat,alloc_sum,i_ind1,i_ind2,i_cont1,i_cont2,&
         l_max,eq_b,eq_a,l_prod

    intrinsic max
    ! ---------- Executable code --------------------------------

    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,la,nb,lb)
    n_indep_a=unique_atoms(na)%symadapt_partner(i_ir,la)%n_independent_fcts
    n_indep_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%n_independent_fcts
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    l_max=max(la,lb)
    alloc_sum=0
    allocate(n_contributing_a(n_indep_a),STAT=alloc_stat)
    alloc_sum=alloc_sum+alloc_stat
    allocate(n_contributing_b(n_indep_b),STAT=alloc_stat)
    alloc_sum=alloc_sum+alloc_stat
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation (1) failed")

    allocate(yl_arr((l_max+1)**2),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC : allocation (2) failed")
    yl_arr = zero  
    
    if (n_indep_a/=0) then
       n_contributing_a=unique_atoms(na)%symadapt_partner(i_ir,la)%&
            &symadapt(:,1)%n_fcts
    endif
    if (n_indep_b/=0) then
       n_contributing_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
            &symadapt(:,1)%n_fcts
    endif

    allocate(fact(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation (3) failed")
    fact=pi/fact0
       
    allocate(fit_int(naexps,nbexps,n_indep_a,n_indep_b),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation of (4) failed")
    fit_int=zero

    if (la.gt.0.and.lb.eq.0) then ! la > 0 and lb = 0 ----------------
       allocate(intermediate1(num,1),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LL_CALC: allocation (5) failed")
       do i_ind1=1,n_indep_a
          intermediate1=zero
          do i_cont1=1,n_contributing_a(i_ind1)
             coeff_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%c  
             eq_atom_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%I_equal_atom
             magn_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%m
             eq_b=1
             xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
             xb = unique_atoms(nb)%position(:,eq_b)
             arg = sum((xa-xb)**2)
             yl_arr = solid_harmonics_scalar(la,xa-xb)
             intermediate1(:,1) = intermediate1(:,1) + &
                  yl_arr(la**2+magn_a(i_cont1))* &
                  coeff_a(i_cont1)*exp(-fact1/fact0*arg)
          enddo

          intermediate1(:,1) = intermediate1(:,1)* &
               pi*sqrt(pi)/(fact0*sqrt(fact0))* &
               (-two*fact1/fact0)**la*          &
               unique_atoms(nb)%n_equal_atoms
          fit_int(:,:,i_ind1,1) = unpack(intermediate1(:,1), &
               cutoff,zero)
       enddo
       deallocate(intermediate1,fact,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("LL_CALC: deallocation (1) failed")

    elseif (la.eq.0.and.lb.gt.0) then ! la=0 and lb > 0
       allocate(intermediate1(num,1),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LL_CALC: allocation (6) failed")
       do i_ind2=1,n_indep_b
          intermediate1=zero
          do i_cont2=1,n_contributing_b(i_ind2)
             coeff_b => &
                  unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                  &symadapt(i_ind2,i_ir)%c
             eq_atom_b => &
                  unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                  &symadapt(i_ind2,i_ir)%I_equal_atom
             magn_b => &
                  unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                  &symadapt(i_ind2,i_ir)%m

             eq_a=1
             xa=unique_atoms(na)%position(:,eq_a)
             xb=unique_atoms(nb)%position(:,eq_atom_b(i_cont2))
             arg = sum((xa-xb)**2)
             yl_arr = solid_harmonics_scalar(lb,xa-xb)
             
             intermediate1(:,1) = intermediate1(:,1) + &
                  yl_arr(lb**2+magn_b(i_cont2))*coeff_b(i_cont2)*&
                  exp(-fact1/fact0*arg)
          enddo
          !ATTENTION: this branch is not tested yet.
          fit_int(:,:,1,i_ind2) = fit_int(:,:,1,i_ind2)+&
               unpack(pi*sqrt(pi)/(fact0*sqrt(fact0))*&
               (two*fact1/fact0)**lb*intermediate1(:,1),&
               cutoff,zero)
          ! There might still be sign errrors, FN 4.8.96
       enddo
       deallocate(intermediate1,fact,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("LL_CALC: deallocation (2) failed")
          
    else  ! la > 0 and lb > 0
       alloc_sum=0
       allocate(intermediate1(num,(lb+1)**2),STAT=alloc_stat)
       alloc_sum=alloc_sum+1
       intermediate1=zero
       allocate(intermediate2(num,(la+1)**2),STAT=alloc_stat)
       alloc_sum=alloc_sum+1
       intermediate2=zero
       allocate(intermediate3(num,1),STAT=alloc_stat)
       alloc_sum=alloc_sum+1
       intermediate3=zero
       allocate(inter1(num,(la+1)**2),STAT=alloc_stat)
       alloc_sum=alloc_sum+1
       inter1=zero
       if(alloc_stat.ne.0) call error_handler &
            ("LL_CALC: allocation (7) failed")

       independent1: do i_ind1=1,n_indep_a
          contributing1: do i_cont1=1,n_contributing_a(i_ind1)
             coeff_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%c
             eq_atom_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%I_equal_atom
             magn_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,la)%&
                  &symadapt(i_ind1,i_ir)%m
             independent2: do i_ind2=1,n_indep_b
                contributing2: do i_cont2=1,n_contributing_b(i_ind2)
                   coeff_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%c
                   eq_atom_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%I_equal_atom
                   magn_b => &
                        unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
                        &symadapt(i_ind2,i_ir)%m

                   ! build intermediate2 = intermediate2(num,1:la(la+2))
                   ! this constitutes the two factors that enter the product
                   ! rule, the third one is the diff_rule below.
                   xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
                   xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont2))
                    
                   arg = sum((xa-xb)**2)
                   yl_arr = solid_harmonics_scalar(l_max,xa-xb)
                   ! this contains all solid harmonics up to la
                   ! and it includes of course all appropriate magnetic
                   ! quantum numbers ---------------------------     

                   counter=1
                   do i_l=0,la
                      do i_m=1,2*i_l+1
                         intermediate2(:,counter) = yl_arr(counter) * &
                              (-one)**i_l*(two*fact1/fact0)**(i_l+lb)
                         counter=counter+1
                      enddo
                   enddo
                   intermediate2 = intermediate2*coeff_a(i_cont1)* &
                        spread(exp(-fact1/fact0*arg),2,(la+1)**2)

                   ! build up intermediate1(num,1:lb(lb+2))
                   ! this variable is used to enter the diff_rule, the
                   ! result of which enters the prod_rule together
                   ! with intermediate2
                   counter=1
                   do i_l=0,lb
                      do i_m=1,2*i_l+1
                         intermediate1(:, counter) = &
                              spread(yl_arr(counter),1,num)
                         counter=counter+1
                      enddo
                   enddo
                   intermediate1 = intermediate1*coeff_b(i_cont2)

                   !now make an intermediate3 out of the two former ones
                   inter1(:,1:(la+1)**2) = diff_rule &
                        (intermediate1,1,(la+1)**2,lb**2+magn_b(i_cont2))

                   l_prod=la**2+magn_a(i_cont1)
                   intermediate3(:,1:1) = &
                        prod_rule(inter1,intermediate2,l_prod,l_prod)

                   intermediate3(:,1) = intermediate3(:,1)* &
                        pi*sqrt(pi)/(fact0*sqrt(fact0))

                   fit_int(:,:,i_ind1,i_ind2) = fit_int(:,:,i_ind1,i_ind2)+&
                        unpack(intermediate3(:,1),cutoff,zero)
                   
                enddo contributing2
             enddo independent2

          enddo contributing1
       enddo independent1
                   
       deallocate(intermediate1,intermediate2,intermediate3,&
            inter1,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("LL_CALC: deallocation (3) failed")
                   
    endif


    deallocate(n_contributing_a,n_contributing_b,yl_arr,&
         STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation of fit_int (1) failed")

    call dump_exponent_data()

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do  LOCAL contraction
    if (output_int_fitcontract) call write_to_output_units &
         ("LL_CALC / XC : calling fitcontract_2c")     
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,flag)
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))

    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation (2) failed")    

  end subroutine ll_calc

  !****************************************************************

  subroutine lr2_calc(na,la,nb,lb)
    ! Purpose: calculates uncontracted, symmetryadapted integral
    !          < fl,fk >.
    !----------------Declaration of formal parameters -----------
    integer(kind=i4_kind),intent(in)   :: na,nb
    integer(kind=i4_kind),intent(in)   :: la,lb
    !** End of interface ****************************************
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: yl_arr(:),exp_fact(:)
    real(kind=r8_kind),allocatable  :: intermediate1(:), &
         inter2(:), inter3(:)
    integer(kind=i4_kind),allocatable  :: n_contributing_a(:)
    integer(kind=i4_kind),allocatable  :: n_contributing_b(:)
    integer(kind=i4_kind),pointer      :: magn_a(:),eq_atom_a(:)
    integer(kind=i4_kind),pointer      :: magn_b(:),eq_atom_b(:)
    real(kind=r8_kind),pointer         :: coeff_b(:)
    real(kind=r8_kind),pointer         :: coeff_a(:)

    ! fitintegral that enters contraction
    real(kind=r8_kind),allocatable  :: fit_int(:,:,:,:) 

    integer(kind=i4_kind) :: eq_a,eq_b,n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind) :: alloc_stat
    integer(kind=i4_kind) :: i_ind,i_cont,l_max
    intrinsic max
    !------------ Executable code -------------------------------- 

    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,la,nb,lb)
    l_max=max(la,lb)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    if(lb.eq.-1) then
       n_indep_a = unique_atoms(na)%symadapt_partner(i_ir,la)%&
            &n_independent_fcts 
       n_indep_b= 1
       allocate(n_contributing_a(n_indep_a),&
            n_contributing_b(1),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("LR2_CALC : allocation (1) failed")
       if (n_indep_a/=0) then
          n_contributing_a=unique_atoms(na)%symadapt_partner(i_ir,la)%&
               &symadapt(:,1)%n_fcts
       endif
       n_contributing_b(1)=n_equal_b
       allocate(fit_int(naexps,nbexps,n_indep_a,1),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (2) failed")
    else
       n_indep_b = unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
            &n_independent_fcts 
       n_indep_a= 1 
       allocate(n_contributing_b(n_indep_b),&
            n_contributing_a(1),STAT=alloc_stat) 
       if( alloc_stat.ne.0) call error_handler &
            ("LR2_CALC : allocation (3) failed")
       if (n_indep_b/=0) then
          n_contributing_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
               &symadapt(:,1)%n_fcts
       endif
       n_contributing_a(1)=n_equal_a
       allocate(fit_int(naexps,nbexps,1,n_indep_b),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (4) failed")
    endif
    fit_int = zero
    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! fact2 = a**2 (if lb=-1)   
    ! fact2 = b**2  (if la=-1)   

    if ((lb.eq.-1.and.la.eq.0).or.(la.eq.-1.and.lb.eq.0)) then
       allocate(intermediate1(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (5) failed")
       intermediate1 = zero

       eq_a=1
       xa=unique_atoms(na)%position(:,eq_a)
       equal_b_1: do eq_b=1,n_equal_b
          xb=unique_atoms(nb)%position(:,eq_b)
          arg = sum((xa-xb)**2)
          intermediate1 = intermediate1 + &
               (three/two + fact2/fact0*arg) * &
               exp(-fact1/fact0*arg)
       enddo equal_b_1

       intermediate1 = n_equal_a*intermediate1*&
            pi*sqrt(pi)/(fact0*fact0*sqrt(fact0))
       fit_int(:,:,1,1) = unpack(intermediate1,cutoff,zero)
       deallocate(intermediate1,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: deallocation (1) failed")
    endif

    if (lb.eq.-1.and.la.gt.0) then  ! =========== lb=-1and la>0 ==========
       allocate(yl_arr((la+1)**2),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (6) failed")
       allocate(exp_fact(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (7) failed")      
       allocate(inter2(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (8) failed") 
       allocate(inter3(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (9) failed") 
       do i_ind = 1,n_indep_a
          inter2=zero
          inter3=zero
          coeff_a => &
               unique_atoms(na)%symadapt_partner(i_ir,la)%&
               &symadapt(i_ind,i_ir)%c
          eq_atom_a => &
               unique_atoms(na)%symadapt_partner(i_ir,la)%&
               &symadapt(i_ind,i_ir)%I_equal_atom
          magn_a => &
               unique_atoms(na)%symadapt_partner(i_ir,la)%&
               &symadapt(i_ind,i_ir)%m
          do i_cont=1,n_contributing_a(i_ind)
             do eq_b=1,n_equal_b
                xa = unique_atoms(na)%position(:,eq_atom_a(i_cont))
                xb = unique_atoms(nb)%position(:,eq_b)
                arg = sum((xa-xb)**2)               
                yl_arr=solid_harmonics_scalar(la,xa-xb)
                exp_fact = exp(-fact1/fact0*arg)
                if (arg.eq.zero) then
                   inter3=inter3+zero
                else
                   inter3 = inter3 + &
                        fact2/fact0*arg*exp_fact*&
                        yl_arr(la**2+magn_a(i_cont))*coeff_a(i_cont)
                endif
                
                inter2 = inter2 + &
                     (three/two-fact2/fact1*real(la,kind=r8_kind))*&
                     exp_fact*coeff_a(i_cont)* &
                     yl_arr(la**2+magn_a(i_cont))
             enddo
          enddo
          fit_int(:,:,i_ind,1) = fit_int(:,:,i_ind,1) + &
               pi*sqrt(pi)* &
               unpack( ( (-two*fact1/fact0)**la*(inter2+inter3) / &
               (fact0*fact0*sqrt(fact0)) ),cutoff,zero)
       enddo
       deallocate(inter2,inter3,yl_arr,exp_fact,&
            STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: deallocation (2) failed")      
    else if (la.eq.-1.and.lb.gt.0) then !======== la=-1 and lb>0 =======

       allocate(yl_arr((lb+1)**2),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (10) failed")
       allocate(exp_fact(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (11) failed")
       allocate(inter2(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (12) failed") 
       allocate(inter3(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: allocation (13) failed") 
       do i_ind =1,n_indep_b
          inter2=zero
          inter3=zero
          coeff_b => &
               unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
               &symadapt(i_ind,i_ir)%c
          eq_atom_b => &
               unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
               &symadapt(i_ind,i_ir)%I_equal_atom
          magn_b => &
               unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
               &symadapt(i_ind,i_ir)%m
          do i_cont = 1,n_contributing_b(i_ind)
             do eq_a = 1,n_equal_a
                xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont))
                xa = unique_atoms(na)%position(:,eq_a)
                arg = sum((xa-xb)**2)
                yl_arr = solid_harmonics_scalar(lb,xb-xa)
                exp_fact = exp(-fact1/fact0*arg)
                inter2 = inter2 + &
                     (three/two-fact2/fact1*real(lb,kind=r8_kind))*&
                     exp_fact*yl_arr(lb**2+magn_b(i_cont))*&
                     coeff_b(i_cont)

                if (arg.eq.zero) then
                   inter3 = inter3 + zero
                else
                   inter3 = inter3 + &
                        fact2/fact0*arg*exp_fact* &
                        yl_arr(lb**2+magn_b(i_cont))*coeff_b(i_cont)
                endif
             enddo
          enddo
          fit_int(:,:,1,i_ind) = fit_int(:,:,1,i_ind) +&
               pi*sqrt(pi)* &
               unpack ( ( (-two*fact1/fact0)**lb*(inter2+inter3)/ &
               (fact0*fact0*sqrt(fact0)) ),cutoff,zero)
       enddo
       deallocate(inter2,inter3,yl_arr,exp_fact,&
            STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC: deallocation (3) failed") 
    endif

    deallocate(n_contributing_a,n_contributing_b,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC: deallocation (4) failed") 
    call dump_exponent_data()

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do  LOCAL contraction
    if (output_int_fitcontract) call write_to_output_units &
         ("LR2_CALC / XC : calling fitcontract_2c")      
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,flag)
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))
    
    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LR2_CALC: deallocation (5) failed") 
    return
  end subroutine lr2_calc

  !*************************************************************

  subroutine r2_calc(na,nb)
    ! Purpose : calculate overlap integrals <fk| fl> for
    !           both angular momenta = -1.
    !
    ! Subroutine called by: hm ?!
    ! Author: FN  8/96
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help arrays
    real(kind=r8_kind),allocatable   :: intermediate3(:)
    real(kind=r8_kind),allocatable   :: intermediate1(:),&
         intermediate2(:),help(:)
    ! symmetry adapted integrals
    real(kind=r8_kind),allocatable   :: fit_int(:,:,:,:)
    integer(kind=i4_kind)            :: eq_a,eq_b,n_equal_a,n_equal_b, &
         alloc_stat
    !------------ Executable code --------------------------------

    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,-1,nb,-1)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    allocate(intermediate3(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (1) failed")
    intermediate3 = zero
    allocate(intermediate1(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (1) failed")
    intermediate1 = zero
    allocate(intermediate2(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (1) failed")
    intermediate2 = zero
    allocate(help(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (1) failed")
    help = zero
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (2) failed")
    fit_int=zero

    eq_a = 1
    xa=unique_atoms(na)%position(:,eq_a)
    do eq_b=1,n_equal_b
       xb=unique_atoms(nb)%position(:,eq_b)
       arg=sum((xa-xb)**2)
       help = exp(-fact1/fact0*arg)
       intermediate1 = intermediate1 + help
       intermediate2 = intermediate2 + arg*help
       intermediate3 = intermediate3 + arg**2*help
    enddo
    help = n_equal_a*pi*sqrt(pi)/(fact0*sqrt(fact0))* &
         ( 3.75_r8_kind/fact0**2*intermediate1 + &
         (three/(two*fact0) - five*fact1/(fact0**3))*intermediate2 + &
         fact1**2/(fact0**4)*intermediate3)

    fit_int(:,:,1,1) = unpack(help,cutoff,zero)
    deallocate(intermediate1,intermediate2,&
         intermediate3,help,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: deallocation (1) failed")        

    call dump_exponent_data()

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do  LOCAL contraction
    if (output_int_fitcontract) call write_to_output_units &
         ("R2_CALC / XC : calling fitcontract_2c")     
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,flag)
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))

    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("R2_CALC: deallocation (2) failed")
    return
  end subroutine r2_calc

  !*************************************************************

  subroutine get_exponent_data(na,la,nb,lb)
    !----------------------------------------------------------------
    !  Purpose: allocate space for the exponent arrys, the help arrays
    !           and get the numbers from the unique_atom_module.
    !           Variables set in this routine:
    !           naexps,nbexps
    !           fact0,fact1,fact2
    !           cutoff,num,i_ir
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
    integer(kind=i4_kind)   :: alloc_stat,alloc_sum
    !------------ Executable code -----------------------------------
    if(la.eq.-1) then
       if (flag(1:2) == "xc") then
       naexps = unique_atoms(na)%r2_xc%n_exponents
       allocate(aexps(naexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (1) failed")
       aexps = unique_atoms(na)%r2_xc%exponents(:)
       else ! flag(1:2) == "ch" 
          naexps = unique_atoms(na)%r2_ch%n_exponents
          allocate(aexps(naexps),STAT=alloc_stat)
          if( alloc_stat.ne.0) call error_handler &
               ("get_exponent_data : allocation (1') failed")
          aexps = unique_atoms(na)%r2_ch%exponents(:)
       endif
    else
       if (flag(1:2) == "xc") then
       naexps = unique_atoms(na)%l_xc(la)%n_exponents
       allocate(aexps(naexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (2) failed")
       aexps = unique_atoms(na)%l_xc(la)%exponents(:)
       else ! flag(1:2) == "ch" 
          naexps = unique_atoms(na)%l_ch(la)%n_exponents
          allocate(aexps(naexps),STAT=alloc_stat)
          if( alloc_stat.ne.0) call error_handler &
               ("get_exponent_data : allocation (2') failed")
          aexps = unique_atoms(na)%l_ch(la)%exponents(:)
       endif
    endif
    if(lb.eq.-1) then
       if (flag(4:5) == "xc") then
       nbexps = unique_atoms(nb)%r2_xc%n_exponents
       allocate(bexps(nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (3) failed")
       bexps = unique_atoms(nb)%r2_xc%exponents(:)
       else ! flag(4:5) == "ch" 
          nbexps = unique_atoms(nb)%r2_ch%n_exponents
          allocate(bexps(nbexps),STAT=alloc_stat)
          if( alloc_stat.ne.0) call error_handler &
               ("get_exponent_data : allocation (3') failed")
          bexps = unique_atoms(nb)%r2_ch%exponents(:)
       endif
    else
       if (flag(4:5) == "xc") then
       nbexps = unique_atoms(nb)%l_xc(lb)%n_exponents
       allocate(bexps(nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (4) failed")
       bexps = unique_atoms(nb)%l_xc(lb)%exponents(:)
       else ! flag(4:5) == "ch" 
          nbexps = unique_atoms(nb)%l_ch(lb)%n_exponents
          allocate(bexps(nbexps),STAT=alloc_stat)
          if( alloc_stat.ne.0) call error_handler &
               ("get_exponent_data : allocation (4') failed")
          bexps = unique_atoms(nb)%l_ch(lb)%exponents(:)
       endif
    endif
    i_ir = get_totalsymmetric_irrep()
    allocate(fact0_arr(naexps,nbexps),STAT=alloc_stat)
    if( alloc_stat.ne.0) call error_handler &
         ("get_exponent_data : allocation (5) failed")
    allocate(fact1_arr(naexps,nbexps),STAT=alloc_stat)
    if( alloc_stat.ne.0) call error_handler &
         ("get_exponent_data : allocation (6) failed")
    allocate(cutoff(naexps,nbexps),STAT=alloc_stat)
    if( alloc_stat.ne.0) call error_handler &
         ("get_exponent_data : allocation (7) failed")
    ! There is no cutoff criterion for the charge overlap
    ! matrix but in order to be able to use the 'prod_rule',
    ! 'diff_rule' etc. we pack th exponents to a linear array
    cutoff = .true.
    num=count(cutoff)
    alloc_sum=0
    allocate (fact0(num),STAT=alloc_stat)
    alloc_sum=alloc_sum+alloc_stat
    allocate (fact1(num),STAT=alloc_stat)
    alloc_sum = alloc_sum+alloc_stat
    if (alloc_sum.ne.0) call error_handler &
         ("get_exponent_data: allocation (8) failed")
    fact0_arr=(spread(aexps,2,nbexps)+spread(bexps,1,naexps))
    fact1_arr=(spread(aexps,2,nbexps)*spread(bexps,1,naexps))
    if(la.eq.-1.or.lb.eq.-1) then
       allocate(fact2_arr(naexps,nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (9) failed")
       allocate (fact2(num),STAT=alloc_stat)
       if(lb.eq.-1) then
          fact2_arr=(spread(aexps**2,2,nbexps))
       else
          fact2_arr=(spread(bexps**2,1,naexps))
       endif
       fact2=pack(fact2_arr,cutoff)
       deallocate(fact2_arr,STAT=alloc_stat)
       if (alloc_stat.ne.0 ) call error_handler &
            ("get_exponent_data: deallocate (1) failed")
    endif
    fact0=pack(fact0_arr,cutoff)
    fact1=pack(fact1_arr,cutoff)

    deallocate(fact0_arr,fact1_arr,STAT=alloc_stat)
    if (alloc_stat.ne.0 ) call error_handler &
         ("get_exponent_data: deallocate (2) failed")
    return
  end subroutine get_exponent_data

  !*************************************************************

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
    deallocate(aexps,bexps,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (1) failed")
    deallocate(fact0,fact1,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (2) failed")
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

end module integral_2c_fit_xc_module

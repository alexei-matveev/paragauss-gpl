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
module  integral_2c_fit_ch_module
  !-------------------------------------------------------------------
  !  Purpose: contains all routines necessary to perform
  !           the two center fitfct integrals:
  !           - [Fk|Fl] (mat_charge)
  !           - norms for fit functions <1|fk>
  !  Module called by: integral_calc_quad_2cff.f90
  !
  !  References: ...
  !  Author: FN
  !  Date: 7/96
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: FN
  ! Date:   19.9.96
  ! Description: Discovered a sign error using the solid
  !              harmonics. Most corrections done in lr2_calc.
  !              One was done in ll_calc.
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: Subroutine adapted to the new meaning of the
  !              argument FLAG of subroutine FITCONTRACT_2C
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use symmetry_data_module, only : get_totalsymmetric_irrep
  use unique_atom_module
  use solid_harmonics_module, only : solid_harmonics_scalar
  use solhrules_module, only: prod_rule, diff_rule
  use gamma_module, only: gamma
  use int_data_2cff_module
  use fitcontract_2c_module, only: fitcontract_2c
  use output_module, only: output_int_fitcontract,output_int_2c_fit
  use iounitadmin_module, only: write_to_output_units
  use time_module, only: stop_timer, start_timer
  use timer_module, only: timer_int_prim_2cff, timer_int_cont_2cff
  use integralpar_module, only: integralpar_i_int_part, integralpar_2cch_pre

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================
  !------------ public functions and subroutines ---------------------
  public charge_overlap


  !===================================================================
  ! End of public interface of module
  !===================================================================
  ! constants
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: two=2.0_r8_kind,three=3.0_r8_kind, &
       one=1.0_r8_kind,four=4.0_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  ! variables needed by all subroutines: these are set in the
  ! routine get_exponent_data(na,la,nb,lb)
  integer(kind=i4_kind)           :: naexps,nbexps
  real(kind=r8_kind),allocatable  :: fact0(:),fact1(:),fact2(:)
  logical,allocatable             :: cutoff(:,:)
  integer(kind=i4_kind)           :: num,i_ir
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************

  subroutine charge_overlap()
    use comm_module, only: comm_i_am_master
    !  Purpose: calculate mat_charge
    !** End of interface *****************************************
    !------------ Modules used --------------------------------
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)     :: la,lb,na,nb
    !------------ Executable code ------------------------------------
 
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

#if OLD_CODE
    if ( (la.gt.0.and.lb.ge.0).or.(la.ge.0.and.lb.gt.0)) then
       if (output_int_2c_fit) call write_to_output_units &
            ("CHARGE_OVERLAP: calling ll_calc")      
       call ll_calc(na,la,nb,lb)
       if(integralpar_2cch_pre) then
          if (output_int_2c_fit) call write_to_output_units &
               ("CHARGE_OVERLAP: calling ll_calc_pre")
          call ll_calc_pre(na,la,nb,lb)          
       end if
    elseif (la.eq.0.and.lb.eq.0) then
       if (output_int_2c_fit) call write_to_output_units &
            ("CHARGE_OVERLAP: calling ss_calc")      
       call ss_calc(na,nb)
       if(integralpar_2cch_pre) then
          if (output_int_2c_fit) call write_to_output_units &
               ("CHARGE_OVERLAP: calling ss_calc_pre")
          call ss_calc_pre(na,nb)          
       end if
    elseif( (la.eq.-1.and.lb.ge.0).or.(lb.eq.-1.and.la.ge.0)) then
       if (output_int_2c_fit) call write_to_output_units &
            ("CHARGE_OVERLAP: calling lr2_calc")      
       call lr2_calc(na,la,nb,lb)
       if(integralpar_2cch_pre) then
          if (output_int_2c_fit) call write_to_output_units &
               ("CHARGE_OVERLAP: calling lr2_calc_pre")
          call lr2_calc_pre(na,la,nb,lb)          
       end if       
    elseif (la.eq.-1.and.lb.eq.-1) then
       if (output_int_2c_fit) call write_to_output_units &
            ("CHARGE_OVERLAP: calling r2_calc")      
       call r2_calc(na,nb)
       if(integralpar_2cch_pre) then
          if (output_int_2c_fit) call write_to_output_units &
               ("CHARGE_OVERLAP: calling r2_calc_pre")
          call r2_calc_pre(na,nb)          
       end if
    else
       call error_handler &
            ("CHARGE_OVERLAP : something fishy")
    endif

#else
!!$    print *,'call ll_calc_v2(',na,la,nb,lb,')'
    call ll_calc_v2(na,la,nb,lb)
#endif
  end subroutine charge_overlap

  !*************************************************************

#if OLD_CODE
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
    real(kind=r8_kind),allocatable  :: intermediate(:),gamma_help(:,:), &
         fit_int(:,:,:,:)
    integer(kind=i4_kind) :: a_eq,b_eq,alloc_stat
    ! ------------- Executable code ---------------------------

    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,0,nb,0)
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: allocation (1) failed")
    fit_int=zero
    allocate(gamma_help(num,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: allocation (2) failed")
    gamma_help=zero
    allocate(intermediate(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: allocation (3) failed")
    intermediate=zero

    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! Loop over ALL equal atoms for symmetry adaption
    a_eq=1
    xa=unique_atoms(na)%position(:,a_eq)
    do b_eq=1,unique_atoms(nb)%n_equal_atoms
       xb=unique_atoms(nb)%position(:,b_eq)
       
       arg=sum((xa-xb)**2)
       gamma_help(:,1:1) = gamma(1,fact1/fact0*arg)
       
       intermediate = intermediate + &
            two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))*gamma_help(:,1)
    enddo
    intermediate = unique_atoms(na)%n_equal_atoms* &
         intermediate
    fit_int(:,:,1,1) = unpack(intermediate,cutoff,zero)

    ! deallocation
    call dump_exponent_data()
    deallocate(intermediate,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("SS_CALC: deallocation (1) failed")
    deallocate(gamma_help,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: deallocation (2) failed")

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do contractions and store result
    if (output_int_fitcontract) call write_to_output_units &
         ("SS_CALC / CH : calling fitcontract_2c")
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,'ch_ch')
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))


    deallocate(fit_int,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC: deallocation (3) failed")

  end subroutine ss_calc

  !*************************************************************

  subroutine ss_calc_pre(na,nb)
    ! Purpose: calculate primitives for both angular momenta
    !          equal zero.
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: intermediate(:) ,&
         fit_int(:,:,:,:)
    integer(kind=i4_kind) :: a_eq,b_eq,alloc_stat
    ! ------------- Executable code ---------------------------


    call get_exponent_data(na,0,nb,0)
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC_PRE: allocation (1) failed")
    fit_int=zero
    allocate(intermediate(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC_PRE: allocation (2) failed")
    intermediate=zero

    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! Loop over ALL equal atoms for symmetry adaption
    a_eq=1
    xa=unique_atoms(na)%position(:,a_eq)
    do b_eq=1,unique_atoms(nb)%n_equal_atoms
       xb=unique_atoms(nb)%position(:,b_eq)
       
       arg=sum((xa-xb)**2)
       intermediate = intermediate + &
            exp(-fact1/fact0*arg)
    enddo
    intermediate = unique_atoms(na)%n_equal_atoms* &
         intermediate*pi*sqrt(pi)/(fact0*sqrt(fact0))
    fit_int(:,:,1,1) = unpack(intermediate,cutoff,zero)
    
    ! deallocation
    call dump_exponent_data()
    deallocate(intermediate,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("SS_CALC_PRE: deallocation (1) failed")

    ! now do contractions and store result
    if (output_int_fitcontract) call write_to_output_units &
         ("SS_CALC_PRE / CH : calling fitcontract_2c")

    call fitcontract_2c(fit_int,'pre')

    deallocate(fit_int,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SS_CALC_PRE: deallocation (2) failed")

  end subroutine ss_calc_pre

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
    real(kind=r8_kind),allocatable,dimension(:,:)  :: gamma_help
    real(kind=r8_kind),allocatable,dimension(:)   :: yl_arr
    ! intermediate arrays
    real(kind=r8_kind),allocatable,dimension(:,:)   :: intermediate1, &
         intermediate2,intermediate3,inter1
    ! uncontracted but symmetryadapted integrals:
    ! fit_int(num,i_ind1,i_ind2)
    real(kind=r8_kind),allocatable    :: fit_int(:,:,:,:)
    ! number of independents and contributing fcts
    integer(kind=i4_kind)   :: n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind),allocatable  :: n_contributing_a(:), &
         n_contributing_b(:)
    real(kind=r8_kind),pointer       :: coeff_a(:),coeff_b(:)
    integer(kind=i4_kind),pointer,dimension(:)   :: magn_a,&
         magn_b,eq_atom_a,eq_atom_b
    real(kind=r8_kind),allocatable :: help_arr(:,:),help_intermediate2(:,:)
    integer(kind=i4_kind)   :: counter,i_l,i_m, &
         alloc_stat,i_ind1,i_ind2,i_cont1,i_cont2,&
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
    allocate(n_contributing_a(n_indep_a), &
         n_contributing_b(n_indep_b), STAT=alloc_stat)
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


    allocate(gamma_help(num,la+lb+1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation of gamma_help failed")
    gamma_help=zero

    allocate(fit_int(naexps,nbexps,n_indep_a,n_indep_b),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation (4) failed")
    fit_int=zero
    allocate(intermediate1(num,1:(lb+1)**2),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC: allocation (5) failed")

    intermediate1=zero

    allocate(help_arr(naexps,nbexps))
    allocate(help_intermediate2(num,(la+1)**2))
    help_arr=zero
    help_intermediate2=zero

    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    if (la.gt.0.and.lb.eq.0) then ! la > 0 and lb = 0 ----------------
       !                           => ( -2 ) **la
       do i_ind1=1,n_indep_a
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
             equal_b: do eq_b=1,unique_atoms(nb)%n_equal_atoms

                xa = unique_atoms(na)%position(:,eq_atom_a(i_cont1))
                xb = unique_atoms(nb)%position(:,eq_b)
                arg = sum((xa-xb)**2)
                gamma_help(:,1:la+1) = gamma(la+1,fact1/fact0*arg)
                yl_arr = solid_harmonics_scalar(la,xa-xb)

                intermediate1(:,1) = coeff_a(i_cont1) * &
                     yl_arr(la**2+magn_a(i_cont1)) *    &
                     gamma_help(:,la+1)

                fit_int(:,:,i_ind1,1) = fit_int(:,:,i_ind1,1) + &
                     unpack(intermediate1(:,1),cutoff,zero)

             enddo equal_b
          enddo
          fit_int(:,:,i_ind1,1) = fit_int(:,:,i_ind1,1) * &
               unpack((two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))* &   
               (-two*fact1/fact0)**la),cutoff,zero)
       enddo
       elseif (la.eq.0.and.lb.gt.0) then ! la=0 and lb > 0
          !                                => (+2)**lb
       do i_ind2=1,n_indep_b
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
             equal_a: do eq_a=1,unique_atoms(na)%n_equal_atoms
                xa = unique_atoms(na)%position(:,eq_a) 
                xb = unique_atoms(nb)%position(:,eq_atom_b(i_cont2))
                arg = sum((xa-xb)**2)
                gamma_help(:,1:lb+1) = gamma(lb+1,fact1/fact0*arg)
                yl_arr = solid_harmonics_scalar(lb,xb-xa)

                intermediate1(:,1) = coeff_b(i_cont2) * &
                     yl_arr(lb**2+magn_b(i_cont2)) * &
                     gamma_help(:,lb+1)

                fit_int(:,:,1,i_ind2) = fit_int(:,:,1,i_ind2) + &
                     unpack(intermediate1(:,1),cutoff,zero)
             enddo equal_a

          enddo
          fit_int(:,:,1,i_ind2) = fit_int(:,:,1,i_ind2) * &
               unpack((two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0)) * &
               (-two*fact1/fact0)**lb),cutoff,zero)
       enddo
    else
       ! allocate the intermediates here ... why should we waste
       ! the storage if we don`t need it ?!
       allocate(intermediate2(num,1:(la+1)**2),intermediate3(num,1), &
            inter1(num,1:(la+1)**2),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LL_CALC: allocation (6) failed")
       intermediate2 = zero
       intermediate3 = zero
       inter1 = zero

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
                   gamma_help(:,1:la+lb+1) = gamma(la+lb+1,fact1/fact0*arg)
                   ! -------------------------------------------
                   yl_arr = solid_harmonics_scalar(l_max,xa-xb)
                   ! this contains all solid harmonics up to la
                   ! and it includes of course all appropriate magnetic
                   ! quantum numbers ---------------------------
                   counter=1
                   ! prepare a whole array the enters the prod_rule
                   do i_l=0,la
                      do i_m=1,2*i_l+1
                         intermediate2(:,counter) =  &
                              yl_arr(counter) * gamma_help(:,lb+i_l+1) * &
                              (two*fact1/fact0)**(i_l+lb)*(-one)**i_l
                         counter=counter+1
                      enddo
                   enddo

                   ! build up intermediate1(num,1:lb(lb+2))
                   ! this variable is used to enter the diff_rule, the
                   ! result of which enters the prod_rule together with
                   ! intermediate2
                   counter=1
                   do i_l=0,lb
                      do i_m=1,2*i_l+1
                         intermediate1(:,counter) = spread(yl_arr(counter),&
                              1,num)
                         counter=counter+1
                      enddo
                   enddo

                   ! now make an intermediate3 out of the former ones
                   ! inter1(num,la(la+2))
                   ! intermediate(num,1))
                   inter1(:,1:(la+1)**2)=diff_rule &
                        (intermediate1,1,(la+1)**2,lb**2+magn_b(i_cont2))

                   inter1=inter1*coeff_b(i_cont2)

                   l_prod=la**2+magn_a(i_cont1)

                   intermediate3(:,1:1) = &
                        prod_rule(inter1,intermediate2,l_prod,l_prod)
                   intermediate3(:,1) = intermediate3(:,1)*two*pi*pi*sqrt(pi) / &
                        (fact1*sqrt(fact0))*coeff_a(i_cont1)

                   fit_int(:,:,i_ind1,i_ind2)=fit_int(:,:,i_ind1,i_ind2) + &
                        unpack(intermediate3(:,1),cutoff,zero)

                enddo contributing2
             enddo independent2

          enddo contributing1
       enddo independent1
       deallocate(intermediate2,intermediate3,inter1,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("LL_CALC: deallocation(1) failed")
    endif


    ! deallocation
    call dump_exponent_data()
    deallocate(intermediate1,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation(2) failed")   
    deallocate(n_contributing_a,n_contributing_b, &
         yl_arr,help_arr,help_intermediate2,gamma_help) 

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! Now do contractions
    if (output_int_fitcontract) call write_to_output_units &
         ("LL_CALC / CH : calling fitcontract_2c")    
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,'ch_ch')
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))

    
    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation(3) failed")    

  end subroutine ll_calc
#endif /* of if OLD_CODE */

  subroutine ll_calc_v2(na,la,nb,lb)
    ! Purpose: calculate primitives for both angular momenta
    !          greater than zero.
    !          Symmetry adaption included.
    !------------ Declaration of formal parameters ---------------
    use shgi, only: shgi_drv_coul2
    use symmetry_data_module, only : symmetry_data_n_irreps, &
         symmetry_data_n_partners 
#ifdef WITH_RESPONSE
    use int_send_2c_resp, only: int_send_2c_resp_save
    use operations_module, only: operations_response
#endif
    use comm_module
    use xpack, only: pck, upck
    use msgtag_module, only: msgtag_resp_2center
    implicit none
    integer(kind=i4_kind),intent(inout)  :: la,lb
    integer(kind=i4_kind),intent(inout)  :: na,nb
    !** End of interface *****************************************

    real(kind=r8_kind),allocatable    :: fit_int(:,:,:,:)
    ! number of independents and contributing fcts
    integer(kind=i4_kind)   :: n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind)   :: n_contr_a, n_contr_b
    real(kind=r8_kind),pointer       :: coeff_a(:),coeff_b(:)
    integer(kind=i4_kind),pointer,dimension(:)   :: magn_a,&
         magn_b,eq_atom_a,eq_atom_b
    integer(kind=i4_kind)   :: alloc_stat,i_ind1,i_ind2,i_cont1,i_cont2,l_max
    real(r8_kind)    :: CFA,CFB,CFAB,NORMAB
    integer(i4_kind) :: EQA,EQB,MA,MB
    integer(i4_kind) :: LLA,LLB
    real(r8_kind),allocatable :: CO(:,:,:,:,:,:)
    integer(i4_kind) :: nir, npa, i_ir, i_pa, la2, lb2
    ! ---------- Executable code --------------------------------

    LLA = LA
    if(LA==-1) then 
       LLA = 0
       naexps = unique_atoms(na)%r2_ch%n_exponents
    else
       naexps = unique_atoms(na)%l_ch(LLA)%n_exponents 
    end if

    LLB = LB
    if(LB==-1) then 
       LLB = 0
       nbexps = unique_atoms(nb)%r2_ch%n_exponents
    else
       nbexps = unique_atoms(nb)%l_ch(LLB)%n_exponents
    end if

    ! The density expands over totally symmetric functions
    ! Therefore do only one irrep A1:
    nir = 1 ! but see below the case of response calculation ...

#ifdef WITH_RESPONSE
    if (operations_response) then
      ! In response calculations one works with the
      ! density response that is expanded also over
      ! not-totally symmetric fit-functions.
      ! Therefore, we evaluate the self-repulsion integrals
      ! for fit-funcitons of all irreps:
      nir = symmetry_data_n_irreps()
    end if
#endif

    i_ir_loop_: do i_ir = 1, nir 

       n_indep_a=unique_atoms(na)%symadapt_partner(i_ir,LLA)%n_independent_fcts
       if (n_indep_a==0) cycle
       n_indep_b=unique_atoms(nb)%symadapt_partner(i_ir,LLB)%n_independent_fcts
       if (n_indep_b==0) cycle

       call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

       n_equal_a=unique_atoms(na)%n_equal_atoms
       n_equal_b=unique_atoms(nb)%n_equal_atoms

       ! FIXME: norm of s- and r2-fit funcs is different
       NORMAB = 1.0_r8_kind
       if(LLA==0) NORMAB = NORMAB * SQRT(REAL(n_equal_a,r8_kind)) !! for r^2
       if(LLB==0) NORMAB = NORMAB * SQRT(REAL(n_equal_b,r8_kind)) !! for r^2

       l_max=max(LLA,LLB)

       ! Hope that the rest of the code handles zero-sized arrays
       ! gracefully:
       ASSERT(naexps   >=0)
       ASSERT(nbexps   >=0)
       ASSERT(n_indep_a>=0)
       ASSERT(n_indep_b>=0)
       allocate(fit_int(naexps,nbexps,n_indep_a,n_indep_b),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       fit_int=zero

       ! FIXME:
       allocate(CO(naexps,nbexps,2*LLA+1,2*LLB+1,n_equal_a,n_equal_b),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       do EQA=1,n_equal_a
          do EQB=1,n_equal_b
             call shgi_drv_coul2(NA,EQA,LA,NB,EQB,LB,CO(:,:,:,:,EQA,EQB))
          enddo
       enddo

!!$       print *,"2c: NA = ",NA, " LA = ",LA,"NB = ",NB,"LB = ",LB, " IRR= ",i_ir
!!$       print *,'2c: NINDA=',n_indep_a,'NINDB=',n_indep_b

       npa = symmetry_data_n_partners(i_ir)

       i_pa_: do i_pa = 1, npa

          independent1: do i_ind1=1,n_indep_a
             coeff_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,LLA)%&
                  &symadapt(i_ind1,i_pa)%c
             eq_atom_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,LLA)%&
                  &symadapt(i_ind1,i_pa)%I_equal_atom
             magn_a => &
                  unique_atoms(na)%symadapt_partner(i_ir,LLA)%&
                  &symadapt(i_ind1,i_pa)%m
!!$             print *,'2c: INDA=',i_ind1,'M=',magn_a,'E=',eq_atom_a
             independent2: do i_ind2=1,n_indep_b
                coeff_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,LLB)%&
                     &symadapt(i_ind2,i_pa)%c
                eq_atom_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,LLB)%&
                     &symadapt(i_ind2,i_pa)%I_equal_atom
                magn_b => &
                     unique_atoms(nb)%symadapt_partner(i_ir,LLB)%&
                     &symadapt(i_ind2,i_pa)%m
!!$                print *,'2c: INDB=',i_ind1,'M=',magn_b,'E=',eq_atom_b

                n_contr_a = unique_atoms(na)%symadapt_partner(i_ir,LLA)%symadapt(i_ind1,i_pa)%n_fcts 
                n_contr_b = unique_atoms(nb)%symadapt_partner(i_ir,LLB)%symadapt(i_ind2,i_pa)%n_fcts
                contributing1: do i_cont1=1,n_contr_a
                   contributing2: do i_cont2=1,n_contr_b

                      EQA = eq_atom_a(i_cont1)
                      EQB = eq_atom_b(i_cont2)
                      MA  = magn_a(i_cont1)
!!$                      ASSERT(MA==0)
                      MB  = magn_b(i_cont2)
!!$                      ASSERT(MB==0)
                      CFA = coeff_a(i_cont1)
                      CFB = coeff_b(i_cont2)
                      CFAB = CFA * CFB * NORMAB

!!$                      print *, 'adding ',CFA,'x',CFB,' of ',MA,MB,EQA,EQB, ' and ',NORMAB

                      fit_int(:,:,i_ind1,i_ind2)=fit_int(:,:,i_ind1,i_ind2) &
                           + CFAB * CO(:,:,MA,MB,EQA,EQB)
                   enddo contributing2
                enddo contributing1

             enddo independent2
          enddo independent1

       end do i_pa_

!!$       print *,"2c: NA = ",NA, " LA = ",LA,"NB = ",NB,"LB = ",LB, "fit_int = ",SUM(fit_int)

       call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

       ! Now do contractions
       if (output_int_fitcontract) call write_to_output_units &
            ("LL_CALC / CH : calling fitcontract_2c")    

       call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
       fit_int=fit_int/npa

       if( i_ir == 1 )then ! Totally symmetric A1 irrep ...
         ! Here the two-center self repulsion integrals over totally
         ! symmetric fit functions are saved in fitcontract_2c_module:
         call fitcontract_2c(fit_int,'ch_ch')
       endif

       call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))

#ifdef WITH_RESPONSE
       !! to tape and forget !!
       if (operations_response) then

          if (la == 0) then
             la2 = -1
          elseif (la == -1) then
             la2 = 0
          else
             la2 = la
          end if

          if (lb == 0) then
             lb2 = -1
          elseif (lb == -1) then
             lb2 = 0
          else
             lb2 = lb
          end if

          call int_send_2c_resp_save(na,nb,la2,lb2,fit_int,i_ir)

       end if
#endif

       deallocate(fit_int,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       deallocate(CO,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

    end do i_ir_loop_

  end subroutine ll_calc_v2

  !*************************************************************

#if OLD_CODE
  subroutine ll_calc_pre(na,la,nb,lb)
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
    ! intermediate arrays
    real(kind=r8_kind),allocatable,dimension(:)   :: intermediate, &
         alpha, alpha2, alpha3, alpha4, beta, beta2, beta3, beta4, &
         help1, help2
    ! uncontracted but symmetryadapted integrals:
    ! fit_int(num,i_ind1,i_ind2)
    real(kind=r8_kind),allocatable    :: fit_int(:,:,:,:)
    ! number of independents and contributing fcts
    integer(kind=i4_kind)   :: n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b, lalpha, lbeta
    integer(kind=i4_kind)   :: a_eq, b_eq, &
         alloc_stat,i_ind_a,i_ind_b
    intrinsic max
    ! ---------- Executable code --------------------------------


    call get_exponent_data(na,la,nb,lb)
    n_indep_a=unique_atoms(na)%symadapt_partner(i_ir,la)%n_independent_fcts
    n_indep_b=unique_atoms(nb)%symadapt_partner(i_ir,lb)%n_independent_fcts
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
 


    allocate(fit_int(naexps,nbexps,n_indep_a,n_indep_b),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC_PRE: allocation (4) failed")
    fit_int=zero
    allocate(intermediate(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC_PRE: allocation (5) failed")

    intermediate=zero
    if(na==nb) then
       ! if the basis funcions are on one center, we do not
       ! calculate, but make only a crude estimation, e. g. these
       ! matrixelements will be never thrown away in the response part.
       fit_int=1.0_r8_kind
       goto 1111
    end if
    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    a_eq=1
    xa=unique_atoms(na)%position(:,a_eq)
    if(la<lb) then
       lalpha=lb
       lbeta=la
    else
       lalpha=la
       lbeta=lb       
    end if
       
    if(lbeta==0) then
       if(lalpha<=2) then
          if(lb<la) then
             allocate(beta2(num),&
                  help1(naexps),help2(nbexps),stat=alloc_stat)             
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (6a) failed")   
             help2=unique_atoms(nb)%l_ch(lb)%exponents(:)**2
             beta2=pack(spread(help2,1,naexps),cutoff)  
          else
             allocate(beta2(num),&
                  help1(nbexps),help2(naexps),stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (6b) failed")   
             help2=unique_atoms(na)%l_ch(la)%exponents(:)**2
             beta2=pack(spread(help2,2,nbexps),cutoff)              
          end if
          do b_eq=1,n_equal_b
             xb=unique_atoms(nb)%position(:,b_eq)
             arg=sum((xa-xb)**2)
             intermediate=intermediate+1.5_r8_kind+arg*beta2/fact0                
          end do
          intermediate=intermediate*exp(-fact1/fact0*arg)/fact0**2.5_r8_kind
          deallocate(beta2,&
               help1,help2,stat=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LL_CALC_PRE: deallocation (6) failed")             
       elseif(lalpha<=4) then
          if(lb<la) then
             allocate(beta2(num),beta4(num),&
                  help1(naexps),help2(nbexps),stat=alloc_stat)             
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (7a) failed")   
             help2=unique_atoms(nb)%l_ch(lb)%exponents(:)**2
             beta2=pack(spread(help2,1,naexps),cutoff)
             help2=help2**2
             beta4=pack(spread(help2,1,naexps),cutoff)             
          else
             allocate(beta2(num),&
                  help1(nbexps),help2(naexps),stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (7b) failed")   
             help2=unique_atoms(na)%l_ch(la)%exponents(:)**2
             beta2=pack(spread(help2,2,nbexps),cutoff)
             help2=help2**2
             beta4=pack(spread(help2,2,nbexps),cutoff)
          end if
          do b_eq=1,n_equal_b
             xb=unique_atoms(nb)%position(:,b_eq)
             arg=sum((xa-xb)**2)
             intermediate=intermediate+(3.75_r8_kind+5.0_r8_kind*arg*beta2/fact0+&
                  arg**2*beta4/fact0/fact0)
          end do
          intermediate=intermediate*exp(-fact1/fact0*arg)/fact0**3.5_r8_kind
          deallocate(beta2,beta4,&
               help1,help2,stat=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LL_CALC_PRE: deallocation (7) failed")             
       else
          call error_handler(&
               'll_calc_pre: Sorry, the estimation of the overlap of the'//&
               ' charge fitting functions is only possible for l up to 4')
       end if
    elseif(lbeta<=2) then
       if(lalpha<=2) then
          allocate(alpha(num),beta(num),alpha2(num),beta2(num),&
               help1(naexps),help2(nbexps),stat=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LL_CALC_PRE: allocation (8) failed")   
          help1=unique_atoms(na)%l_ch(la)%exponents(:)
          alpha=pack(spread(help1,2,nbexps),cutoff)
          help1=help1**2
          alpha2=pack(spread(help1,2,nbexps),cutoff)
          help2=unique_atoms(nb)%l_ch(lb)%exponents(:)
          beta=pack(spread(help2,1,naexps),cutoff)
          help2=help2**2
          beta2=pack(spread(help2,1,naexps),cutoff)   
          do b_eq=1,n_equal_b
             xb=unique_atoms(nb)%position(:,b_eq)
             arg=sum((xa-xb)**2)
             intermediate=intermediate+3.75_r8_kind+arg/fact0*&
                  (1.5_r8_kind*(alpha2+beta2)-2.0_r8_kind*alpha*beta)&
                  +arg**2*alpha2*beta2/fact0/fact0
          end do
          intermediate=intermediate*exp(-fact1/fact0*arg)/fact0**3.5_r8_kind
          deallocate(alpha,beta,alpha2,beta2,&
               help1,help2,stat=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LL_CALC_PRE: deallocation (8) failed")           
       elseif(lalpha<=4) then
          if(lb<la) then
             allocate(alpha(num),beta(num),alpha2(num),beta2(num),beta3(num),beta4(num),&
                  help1(naexps),help2(nbexps),stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (9a) failed")             
             help1=unique_atoms(na)%l_ch(la)%exponents(:)
             alpha=pack(spread(help1,2,nbexps),cutoff)
             help1=help1**2
             alpha2=pack(spread(help1,2,nbexps),cutoff)
             help2=unique_atoms(nb)%l_ch(lb)%exponents(:)
             beta=pack(spread(help2,1,naexps),cutoff)
             help2=help2**2
             beta2=pack(spread(help2,1,naexps),cutoff)
             help2=help2*unique_atoms(nb)%l_ch(lb)%exponents(:)
             beta3=pack(spread(help2,1,naexps),cutoff)
             help2=help2*unique_atoms(nb)%l_ch(lb)%exponents(:)
             beta4=pack(spread(help2,1,naexps),cutoff)
          else
             allocate(alpha(num),beta(num),alpha2(num),beta2(num),beta3(num),beta4(num),&
                  help1(nbexps),help2(naexps),stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (9b) failed")             
             help1=unique_atoms(nb)%l_ch(lb)%exponents(:)
             alpha=pack(spread(help1,1,naexps),cutoff)
             help1=help1**2
             alpha2=pack(spread(help1,1,naexps),cutoff)
             help2=unique_atoms(na)%l_ch(la)%exponents(:)
             beta=pack(spread(help2,2,nbexps),cutoff)
             help2=help2**2
             beta2=pack(spread(help2,2,nbexps),cutoff)
             help2=help2*unique_atoms(na)%l_ch(la)%exponents(:)
             beta3=pack(spread(help2,2,nbexps),cutoff)
             help2=help2*unique_atoms(na)%l_ch(la)%exponents(:)
             beta4=pack(spread(help2,2,nbexps),cutoff)             
          end if
          do b_eq=1,n_equal_b
             xb=unique_atoms(nb)%position(:,b_eq)
             arg=sum((xa-xb)**2)
             intermediate=intermediate+13.125_r8_kind+5.0_r8_kind*arg/fact0*&
                  (0.75_r8_kind*alpha2+2.5_r8_kind*beta2-2.0_r8_kind*alpha*beta)+&
                  arg**2/fact0/fact0*&
                  (5.0_r8_kind*alpha2*beta2-4.0_r8_kind*alpha*beta3+1.5_r8_kind*beta4)+&
                  arg**3*alpha2*beta4/fact0/fact0/fact0
          end do
          intermediate=intermediate*exp(-fact1/fact0*arg)/fact0**4.5_r8_kind
          deallocate(alpha,beta,alpha2,beta2,beta3,beta4,help1,help2,stat=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("LL_CALC_PRE: allocation (9b) failed")           
       else
          call error_handler(&
               'll_calc_pre: Sorry, the estimation of the overlap of the'//&
               ' charge fitting functions is only possible for l up to 4')
       end if
    elseif(lbeta<=4) then
       allocate(alpha(num),beta(num),alpha2(num),beta2(num),alpha3(num),&
            beta3(num),alpha4(num),beta4(num),help1(naexps),help2(nbexps), &
            stat=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
         ("LL_CALC_PRE: allocation (10) failed")   
       help1=unique_atoms(na)%l_ch(la)%exponents(:)
       alpha=pack(spread(help1,2,nbexps),cutoff)
       help1=help1**2
       alpha2=pack(spread(help1,2,nbexps),cutoff)
       help1=help1*unique_atoms(na)%l_ch(la)%exponents(:)
       alpha3=pack(spread(help1,2,nbexps),cutoff)
       help1=help1*unique_atoms(na)%l_ch(la)%exponents(:)
       alpha4=pack(spread(help1,2,nbexps),cutoff)
       help2=unique_atoms(nb)%l_ch(lb)%exponents(:)
       beta=pack(spread(help2,1,naexps),cutoff)
       help2=help2**2
       beta2=pack(spread(help2,1,naexps),cutoff)
       help2=help2*unique_atoms(nb)%l_ch(lb)%exponents(:)
       beta3=pack(spread(help2,1,naexps),cutoff)
       help2=help2*unique_atoms(nb)%l_ch(lb)%exponents(:)
       beta4=pack(spread(help2,1,naexps),cutoff)       
    
       do b_eq=1,n_equal_b
          xb=unique_atoms(nb)%position(:,b_eq)
          arg=sum((xa-xb)**2)
          intermediate=intermediate+59.0625_r8_kind+5.0_r8_kind*arg/fact0*&
               (8.75_r8_kind*(alpha2+beta2)-14.0_r8_kind*alpha*beta)+&
               (arg/fact0)**2*(3.75*(alpha4+beta4)-20.0_r8_kind*&
               (alpha3*beta+alpha*beta3)+47.0_r8_kind*alpha2*beta2)+&
               (arg/fact0)**3*&
               (5.0_r8_kind*(alpha4*beta2+alpha2*beta4)-8.0_r8_kind*alpha3*beta3)+&
               (arg/fact0)**4*alpha4*beta4
       end do
       intermediate=intermediate*exp(-fact1/fact0*arg)/fact0**5.5_r8_kind
       deallocate(alpha,beta,alpha2,beta2,alpha3,beta3,alpha4,beta4,&
            stat=alloc_stat)       
       if (alloc_stat.ne.0) call error_handler &
            ("LL_CALC_PRE: deallocation (10) failed")         
    else
       call error_handler(&
            'll_calc_pre: Sorry, the estimation of the overlap of the'//&
            ' charge fitting functions is only possible for l up to 4')       
    end if

    intermediate = n_equal_a*intermediate*pi*sqrt(pi)
    do i_ind_b=1,n_indep_b
       do i_ind_a=1,n_indep_a
          fit_int(:,:,i_ind_a,i_ind_b) = unpack(intermediate,cutoff,zero)
       end do
    end do
    
    ! deallocation
1111 call dump_exponent_data()
    deallocate(intermediate,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation(2) failed")   
    

    ! Now do contractions
    if (output_int_fitcontract) call write_to_output_units &
         ("LL_CALC_PRE / CH : calling fitcontract_2c")    

    call fitcontract_2c(fit_int,'pre')

    
    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LL_CALC: deallocation(3) failed")    

  end subroutine ll_calc_pre

  !*************************************************************

  subroutine lr2_calc(na,la,nb,lb)
    ! Purpose: calculate primitives for one angular momentum
    !          greater equal zero and the other being r2-type
    !          (l=-1).
    !          Symmetry adaption included.
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(inout)  :: na,nb
    integer(kind=i4_kind),intent(inout)  :: la,lb
    !** End of interface *****************************************
    ! atomic positions
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: gamma_help(:,:)
    real(kind=r8_kind),allocatable  :: yl_arr(:)
    real(kind=r8_kind),allocatable  :: inter2(:),inter3(:)

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
    integer(kind=i4_kind) :: i_ind,i_cont,l_max,max_order
!!$    integer(kind=i4_kind) :: n3, n4
    intrinsic max
    !------------ Executable code ------------------------------------

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
       if(n_indep_b/=0) then
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
    ! fact2 = a/b  (if lb=-1)
    ! fact2 = b/a  (if la=-1)
    max_order=max(2,l_max+2)
    allocate(gamma_help(num,max_order),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LR2_CALC: allocation (5) failed")
    gamma_help=zero

    if (lb.eq.-1) then     ! ========== lb = r2 ==============================
       if (la.eq.0) then   ! la = 0 and lb = r2

          equal_a_1: do eq_a=1,n_equal_a
             equal_b_1: do eq_b=1,n_equal_b
                xa=unique_atoms(na)%position(:,eq_a)
                xb=unique_atoms(nb)%position(:,eq_b)
                arg = sum((xa-xb)**2)
                gamma_help(:,1:2) = gamma(2,fact1/fact0*arg)

                fit_int(:,:,1,1) = fit_int(:,:,1,1) + &
                     unpack(((three/two+fact2)*gamma_help(:,1) + &
                     fact1*fact2/fact0 * arg * gamma_help(:,2)),cutoff,zero)
             enddo equal_b_1
          enddo equal_a_1
          fit_int(:,:,1,1) = fit_int(:,:,1,1) * &
               two*pi*pi*sqrt(pi)/          &
               unpack((fact1*fact0*sqrt(fact0)),cutoff,zero)
       else   ! la > 0 and lb=r2
          allocate(yl_arr((la+1)**2),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (6) failed")
          yl_arr = zero
          allocate(inter2(num),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (7) failed")
          inter2=zero
          allocate(inter3(num),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (8) failed")
          inter3=zero

          do i_ind =1,n_indep_a
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
             do i_cont = 1,n_contributing_a(i_ind)
                do eq_b = 1,n_equal_b
                   xa = unique_atoms(na)%position(:,eq_atom_a(i_cont))
                   xb = unique_atoms(nb)%position(:,eq_b)
                   arg = sum((xa-xb)**2)
                   yl_arr = solid_harmonics_scalar(la,xa-xb)
                   gamma_help(:,1:la+2) = gamma(la+2,fact1/fact0*arg)

                   inter2 = inter2 + &
                        (three/two-fact2*real(la-1,kind=r8_kind)) * &
                        gamma_help(:,la+1)*yl_arr(la**2+magn_a(i_cont))* &
                        coeff_a(i_cont)
                   if (arg.eq.zero) then
                      inter3 = inter3 + zero
                   else
                      inter3 = inter3 + &
                           fact2*fact1/fact0*arg*gamma_help(:,la+2) * &
                           yl_arr(la**2+magn_a(i_cont)) * &
                           coeff_a(i_cont)
                   endif
  
                enddo
             enddo
             fit_int(:,:,i_ind,1) = fit_int(:,:,i_ind,1) + &
                  two*pi*pi*sqrt(pi) *&
                  unpack( ( (-two*fact1/fact0)**la * (inter2+inter3) /&
                  (fact1*fact0*sqrt(fact0)) ),cutoff,zero)
             
          enddo
          ! deallocate
          deallocate(yl_arr,inter2,inter3,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: deallocation (1) failed")
       endif

    else if (la.eq.-1) then  ! ============== la = r2 =====================
 
       if (lb.eq.0) then ! la=r2 and lb=0
          eq_a=1
          xa=unique_atoms(na)%position(:,eq_a)
          equal_b_2: do eq_b=1,n_equal_b
             xb=unique_atoms(nb)%position(:,eq_b)
             arg = sum((xa-xb)**2)
             gamma_help(:,1:2) = gamma(2,fact1/fact0*arg)
             fit_int(:,:,1,1) = fit_int(:,:,1,1) + &
                  unpack(((three/two+fact2)*gamma_help(:,1) + &
                  fact1*fact2/fact0 * arg * gamma_help(:,2)),cutoff,zero)
          enddo equal_b_2
          fit_int(:,:,1,1) = fit_int(:,:,1,1) * &
               n_equal_a*two*pi*pi*sqrt(pi)/          &
               unpack((fact1*fact0*sqrt(fact0)),cutoff,zero)

       else   ! la=r2 and lb > 0

          allocate(yl_arr((lb+1)**2),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (9) failed")
          yl_arr = zero
          allocate(inter2(num),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (10) failed")
          inter2 = zero
          allocate(inter3(num),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: allocation (11) failed")
          inter3 = zero

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
                   gamma_help(:,1:lb+2) = gamma(lb+2,fact1/fact0*arg)
                   
                   inter2 = inter2 + &
                        (three/two-fact2*real(lb-1,kind=r8_kind)) * &
                        gamma_help(:,lb+1)*yl_arr(lb**2+magn_b(i_cont)) * &
                        coeff_b(i_cont)
                   if (arg.eq.zero) then
                      inter3 = inter3 + zero
                   else
                      inter3 = inter3 + &
                           fact2*fact1/fact0*arg*gamma_help(:,lb+2) * &
                           yl_arr(lb**2+magn_b(i_cont)) * &
                           coeff_b(i_cont)
                   endif

                enddo
             enddo
             fit_int(:,:,1,i_ind) = fit_int(:,:,1,i_ind) + &
                  two*pi*pi*sqrt(pi) * &
                  unpack( ( (-two*fact1/fact0)**lb*(inter2+inter3)/ &
                  (fact1*fact0*sqrt(fact0)) ),cutoff,zero)
          enddo
          ! deallocate
          deallocate(yl_arr,inter2,inter3,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("LR2_CALC: deallocation (2) failed")
       endif
    endif
    deallocate(gamma_help,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC: deallocation (3) failed")
    deallocate(n_contributing_a,n_contributing_b,STAT=alloc_stat)
    if ( alloc_stat.ne.0) call error_handler &
         ("LR2_CALC : deallocation (4) failed")

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do contractions
    if (output_int_fitcontract) call write_to_output_units &
         ("LR2_CALC / CH : calling fitcontract_2c")    
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,'ch_ch')
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))


    deallocate(fit_int,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC: deallocation (5) failed")
    call dump_exponent_data()

    return
  end subroutine lr2_calc

  !*************************************************************

  subroutine lr2_calc_pre(na,la,nb,lb)
    ! Purpose: calculates uncontracted, symmetryadapted integral
    !          < fl,fk >.
    !----------------Declaration of formal parameters -----------
    integer(kind=i4_kind),intent(in)   :: na,nb
    integer(kind=i4_kind),intent(in)   :: la,lb
    !** End of interface ****************************************
    real(kind=r8_kind),dimension(3)  :: xa,xb
    real(kind=r8_kind)               :: arg
    ! help factors
    real(kind=r8_kind),allocatable  :: intermediate1(:), fact3(:)

    ! fitintegral that enters contraction
    real(kind=r8_kind),allocatable  :: fit_int(:,:,:,:), fact2_arr(:,:), &
         fact3_arr(:,:)

    integer(kind=i4_kind) :: eq_a,eq_b,n_indep_a,n_indep_b, &
         n_equal_a,n_equal_b
    integer(kind=i4_kind) :: alloc_stat
    integer(kind=i4_kind) :: i_ind_a, i_ind_b
    intrinsic max
    !---begin debugging
    integer :: was_here
    !---end debugging
    !------------ Executable code -------------------------------- 

    !---begin debugging
    was_here = 0
    !---end debugging

    call get_exponent_data(na,la,nb,lb)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    allocate(fact2_arr(naexps,nbexps),fact3_arr(naexps,nbexps), &
         fact3(num), STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC_PRE: allocation (1) failed")    
    if(lb.eq.-1) then
       n_indep_a = unique_atoms(na)%symadapt_partner(i_ir,la)%&
            &n_independent_fcts 
       n_indep_b= 1
       allocate(fit_int(naexps,nbexps,n_indep_a,1),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC_PRE: allocation (2) failed")
       fact3_arr=spread(unique_atoms(na)%l_ch(la)%exponents(:),2,nbexps)
       fact2_arr=spread(unique_atoms(nb)%l_ch(lb)%exponents(:),1,naexps)
    else
       n_indep_b = unique_atoms(nb)%symadapt_partner(i_ir,lb)%&
            &n_independent_fcts 
       n_indep_a= 1 
       allocate(fit_int(naexps,nbexps,1,n_indep_b),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("LR2_CALC_PRE: allocation (4) failed")
       fact3_arr=spread(unique_atoms(nb)%l_ch(lb)%exponents(:),1,naexps)
       fact2_arr=spread(unique_atoms(na)%l_ch(la)%exponents(:),2,nbexps)
    endif
    fact2=pack(fact2_arr,cutoff)
    fact3=pack(fact3_arr,cutoff)
    fit_int = zero
    ! List of *facts* at the beginning
    ! fact0 = a + b
    ! fact1 = a * b
    ! fact2 = a (if la=-1)   
    ! fact3 = b (if la=-1)       
    ! fact2 = b (if lb=-1)   
    ! fact3 = a (if lb=-1)   

    allocate(intermediate1(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC_PRE: allocation (5) failed")
    intermediate1 = zero

    if((lb==-1.and.la==0).or.(lb==0.and.la==-1)) then
       eq_a=1
       xa=unique_atoms(na)%position(:,eq_a)
       equal_b_1: do eq_b=1,n_equal_b
          xb=unique_atoms(nb)%position(:,eq_b)
          arg = sum((xa-xb)**2)
          intermediate1 = intermediate1 + &
               (three/two + fact3*fact3/fact0*arg) * &
               exp(-fact1/fact0*arg)
       enddo equal_b_1

       intermediate1 = intermediate1/( fact0*fact0*sqrt(fact0) )

    elseif((lb==-1.and.la<=2).or.(la==-1.and.lb<=2)) then
       eq_a=1
       xa=unique_atoms(na)%position(:,eq_a)
       equal_b_2: do eq_b=1,n_equal_b
          xb=unique_atoms(nb)%position(:,eq_b)
          arg = sum((xa-xb)**2)
          intermediate1 = intermediate1 + (&
               & 3.75_r8_kind &
               & + arg*1.5_r8_kind*( fact0**2 - (10.0_r8_kind/3.0_r8_kind)*fact1)/fact0&
               & + (arg*fact1/fact0)**2&
               & )*exp(-arg*fact1/fact0)
       enddo equal_b_2
       
       intermediate1=intermediate1/fact0**3.5 
    else
       eq_a=1
       xa=unique_atoms(na)%position(:,eq_a)
       equal_b_3: do eq_b=1,n_equal_b
          xb=unique_atoms(nb)%position(:,eq_b)
          arg = sum((xa-xb)**2)
          intermediate1 = intermediate1 + (&
               & 13.125_r8_kind&
               & + (2.5_r8_kind*fact2*fact2-2.0_r8_kind*fact1+0.75_r8_kind*fact3*fact3)&
               &   *5.0_r8_kind*(arg/fact0)&
               & + (1.5_r8_kind*fact2**4 - 4.0_r8_kind*fact1*fact2*fact2&
               &    +5.0_r8_kind*fact1*fact1) * (arg/fact0)**2&
               & + fact1*fact1*fact2*fact2*(arg/fact0)**3&
               & ) * exp(-arg*fact1/fact0)          
       enddo equal_b_3
       intermediate1=intermediate1/fact0**4.5 
    end if
    intermediate1=intermediate1*n_equal_a*pi*sqrt(pi)
    do i_ind_b=1,n_indep_b
       do i_ind_a=1,n_indep_a
          fit_int(:,:,i_ind_a,i_ind_b) = unpack(intermediate1,cutoff,zero)
       end do
    end do
    
    deallocate(intermediate1,fact2_arr, fact3_arr, fact3, &
         STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("LR2_CALC_PRE: deallocation (1) failed")
    
    call dump_exponent_data()

    ! now do  LOCAL contraction 

    call fitcontract_2c(fit_int,'pre')
    
    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("LR2_CALC_PRE: deallocation (5) failed") 
    return
  end subroutine lr2_calc_pre

  !*************************************************************

  subroutine r2_calc(na,nb)
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
    real(kind=r8_kind),allocatable   :: intermediate(:), &
         gamma_help(:,:)
    ! symmetry adapted integrals
    real(kind=r8_kind),allocatable   :: fit_int(:,:,:,:)
    integer(kind=i4_kind)            :: eq_a,eq_b,n_equal_a,n_equal_b, &
         alloc_stat
!!$    integer(kind=i4_kind)            :: n1, n2
    !------------ Executable code ------------------------------------

    call start_timer(timer_int_prim_2cff(integralpar_i_int_part))

    call get_exponent_data(na,-1,nb,-1)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms
    allocate(gamma_help(num,3),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (1) failed")
    allocate(intermediate(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (2) failed")
    intermediate = zero
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC: allocation (3) failed")
    fit_int=zero

    eq_a=1
    xa=unique_atoms(na)%position(:,eq_a)       
    equal_b: do eq_b=1,n_equal_b
       xb=unique_atoms(nb)%position(:,eq_b)
       arg = sum((xa-xb)**2)
       gamma_help(:,1:3) = gamma(3,fact1/fact0*arg)           
       intermediate = intermediate + &
            ((three/(two*fact1) + three/(four*fact0**2)) * &
            gamma_help(:,1) + &
            (three/(two*fact0) - three*fact1/fact0**3) * arg * &
            gamma_help(:,2) + &
            (fact1/fact0**2 * arg) **2 * gamma_help(:,3)) * &
            two*pi*pi*sqrt(pi)/(fact1*sqrt(fact0))
    enddo equal_b
    fit_int(:,:,1,1) = fit_int(:,:,1,1) + &
         unpack(n_equal_a*intermediate,cutoff,zero)
    call dump_exponent_data()
    deallocate(intermediate,gamma_help,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("R2_CALC: deallocation (1) failed")

    call stop_timer(timer_int_prim_2cff(integralpar_i_int_part))

    ! now do LOCAL contractions
    if (output_int_fitcontract) call write_to_output_units &
         ("R2_CALC / CH : calling fitcontract_2c")    
    call start_timer(timer_int_cont_2cff(integralpar_i_int_part))
    call fitcontract_2c(fit_int,'ch_ch')
    call stop_timer(timer_int_cont_2cff(integralpar_i_int_part))

    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("R2_CALC: deallocation (2) failed")    
    return
  end subroutine r2_calc

  !*************************************************************

  subroutine r2_calc_pre(na,nb)
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
    real(kind=r8_kind),allocatable   :: intermediate1(:), &
         intermediate2(:), intermediate3(:), help(:)
    ! symmetry adapted integrals
    real(kind=r8_kind),allocatable   :: fit_int(:,:,:,:)
    integer(kind=i4_kind)            :: eq_a,eq_b,n_equal_a,n_equal_b, &
         alloc_stat
    !------------ Executable code ------------------------------------


    call get_exponent_data(na,-1,nb,-1)
    n_equal_a=unique_atoms(na)%n_equal_atoms
    n_equal_b=unique_atoms(nb)%n_equal_atoms

    allocate(intermediate1(num), intermediate2(num), &
         intermediate3(num), help(num), STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC_PRE: allocation (2) failed")
    intermediate1 = zero
    intermediate2 = zero
    intermediate3 = zero
    allocate(fit_int(naexps,nbexps,1,1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC_PRE: allocation (3) failed")
    fit_int=zero

    eq_a=1
    xa=unique_atoms(na)%position(:,eq_a)       
    equal_b: do eq_b=1,n_equal_b
       xb=unique_atoms(nb)%position(:,eq_b)
       arg = sum((xa-xb)**2)
       help = exp(-arg*fact1/fact0)           
       intermediate1 = intermediate1 + help
       intermediate2 = intermediate2 + help*arg
       intermediate3 = intermediate3 + help*(arg**2)
    enddo equal_b
    help = ( intermediate1*3.75_r8_kind/fact0**2 &
         & + intermediate2*1.5_r8_kind*( (fact0**2) - (10.0_r8_kind/3.0_r8_kind)*fact1 )/(fact0**3)&
         & + fact1**2/(fact0**4)*intermediate3)
    help = help * n_equal_a*pi*sqrt(pi)/(fact0*sqrt(fact0))

    fit_int(:,:,1,1) = unpack(help,cutoff,zero)
    deallocate(intermediate1,intermediate2,&
         intermediate3,help,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("R2_CALC_PRE: deallocation (1) failed")

    call dump_exponent_data()

    if(alloc_stat.ne.0) call error_handler &
         ("R2_CALC_PRE: deallocation (1) failed")

    ! now do LOCAL contractions
    if (output_int_fitcontract) call write_to_output_units &
         ("R2_CALC_PRE / CH : calling fitcontract_2c")    

    call fitcontract_2c(fit_int,'pre')


    deallocate(fit_int,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("R2_CALC_PRE: deallocation (2) failed")    
    return
  end subroutine r2_calc_pre

  !*************************************************************

  subroutine get_exponent_data(na,la,nb,lb)
    !----------------------------------------------------------------
    !  Purpose: allocate space for the exponent arrys, the help arrays
    !           and get the numbers from the unique_atom_module.
    !           Variables set in this routine:
    !           naexps,nbexps
    !           aexps,bexps
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
    real(kind=r8_kind),allocatable  :: aexps(:),bexps(:)
    integer(kind=i4_kind)   :: alloc_stat
    !------------ Executable code -----------------------------------
    if(la.eq.-1) then
       naexps = unique_atoms(na)%r2_ch%n_exponents
       allocate(aexps(naexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (1) failed")
       aexps = unique_atoms(na)%r2_ch%exponents(:)
    else
       naexps = unique_atoms(na)%l_ch(la)%n_exponents
       allocate(aexps(naexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (2) failed")
       aexps = unique_atoms(na)%l_ch(la)%exponents(:)
    endif
    if(lb.eq.-1) then
       nbexps = unique_atoms(nb)%r2_ch%n_exponents
       allocate(bexps(nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (3) failed")
       bexps = unique_atoms(nb)%r2_ch%exponents(:)
    else
       nbexps = unique_atoms(nb)%l_ch(lb)%n_exponents
       allocate(bexps(nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (4) failed")
       bexps = unique_atoms(nb)%l_ch(lb)%exponents(:)
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
    allocate (fact0(num),fact1(num),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("get_exponent_data: allocation (8) failed")
    fact0_arr=(spread(aexps,2,nbexps)+spread(bexps,1,naexps))
    fact1_arr=(spread(aexps,2,nbexps)*spread(bexps,1,naexps))
    if(la.eq.-1.or.lb.eq.-1) then
       allocate(fact2_arr(naexps,nbexps),STAT=alloc_stat)
       if( alloc_stat.ne.0) call error_handler &
            ("get_exponent_data : allocation (9) failed")
       allocate (fact2(num),STAT=alloc_stat)
       if(lb.eq.-1) then
          fact2_arr=(spread(aexps,2,nbexps)/spread(bexps,1,naexps))
       else
          fact2_arr=(spread(bexps,1,naexps)/spread(aexps,2,nbexps))
       endif
       fact2=pack(fact2_arr,cutoff)
       deallocate(fact2_arr,STAT=alloc_stat)
       if (alloc_stat.ne.0 ) call error_handler &
            ("get_exponent_data: deallocate (1) failed")
    endif
    fact0=pack(fact0_arr,cutoff)
    fact1=pack(fact1_arr,cutoff)

    deallocate(fact0_arr,fact1_arr,aexps,bexps,STAT=alloc_stat)
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

    if (allocated(fact0) .and. allocated(fact1)) then
    deallocate(fact0,fact1,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (1) failed")
    end if
    if(allocated(fact2)) then
       deallocate(fact2,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("dump_exponent_data: deallocate (3) failed")
    endif
    if (allocated(cutoff)) then
    deallocate(cutoff,STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("dump_exponent_data: deallocate (4) failed")
    end if
  end subroutine dump_exponent_data

  !*************************************************************
#endif /* of if OLD_CODE */

end module integral_2c_fit_ch_module

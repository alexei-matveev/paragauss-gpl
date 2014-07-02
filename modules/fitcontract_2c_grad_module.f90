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
!define FPP_TIMERS 2
module fitcontract_2c_grad_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Perform fitcontractions + renormaliation +
  !           multiplication with fitcoefficients +
  !           summation over fitfunctions for the
  !           two-center fitfunctions integrals in the
  !           gradient part:
  !           F_fit = sum_{k,l}{ ak al Dx [F_k | F_l] }
  !
  !
  !  Module called by: integral_2c_fit_ch_grad
  !
  !  Author: FN
  !  Date: 12/96
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author:
  ! Date:
  ! Description:
  !
  !-------------------------------------------------------------------
#include "def.h"
  use type_module ! type specification parameters
#ifdef _COMPAC_FORTRAN
  use datatype
#endif
  use unique_atom_module
  use int_data_2cff_module
  use output_module, only: output_int_fitcontract
  use fit_coeff_module, only: coeff_charge, coeff_spin, coeff_xcmda,coeff_xcmda_ts
  use gradient_data_module, only: grad_fit_ch_cartesian, cpks_grad_fit_ch_cartesian, &
                                  grad_mda_rhofit_cartesian, &
                                  grad_mda_xcpot_cartesian
USE_MEMLOG


  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public fitcontract_2c_grad,fitcontract_2c_dervs


  !===================================================================
  ! End of public interface of module
  !===================================================================

!..............................................................................
! << OUTPUT ARRAYS >>
! ===================
! F1 = [ rhofit | d/dR rhofit ]
!    = Sum(k,l) a_k,tot a_l,tot [ f_k | d/dR f_l ]
!
! options_xcmode == xcmode_model_density
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! F2 = Sum(s) < rhofit_s | Sum(l) b_l,s V_H[ d/dR f_l ] >
!    = Sum(k,l,s) a_k,s b_l,s [ f_k | d/dR f_l ]
!    = Sum(k,l,t) a_k,t b_l,t [ f_k | d/dR f_l ]
!      with t = tot,spin
! F3 = < Sum(k) a_k,s d/dR f_k | V_X,s >
!    = Sum(k,l,s) a_k,s b_l,s [ d/dR f_k | f_l ]
!    = Sum(k,l,t) a_k,t b_l,t [ d/dR f_k | f_l ]
!      with t = tot,spin
!
! Storage (negative d/dR E_tot contributions)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! options_split_gradients()             | FALSE  FALSE     TRUE   TRUE
! options_xcmode() == xcmode_model_dens.| FALSE  TRUE      FALSE  TRUE
! --------------------------------------+-----------------------------
! grad_fit_ch_cartesian    (1:3,1:n_ua) | F1     F1+F2+F3  F1     F1
! grad_mda_xcpot_cartesian (1:3,1:n_ua) | --     --        --     F2
! grad_mda_rhofit_cartesian(1:3,1:n_ua) | --     --        --     F3
!
! << WORKING ARRAYS >>
! ====================
! grad_fit(1:3) -- Sum(k,l) a_k,tot [ d/dR f_k |      f_l ] a_l,tot   and
!                  Sum(k,l) a_k,tot [      f_k | d/dR f_l ] a_l,tot
! grad_vxc(1:3) -- Sum(k,l) b_k,t   [ d/dR f_k |      f_l ] a_l,t     and
!                  Sum(k,l) a_k,t   [      f_k | d/dR f_l ] b_l,t
! grad_rho(1:3) -- Sum(k,l) a_k,t   [ d/dR f_k |      f_l ] b_l,t     and
!                  Sum(k,l) b_k,t   [      f_k | d/dR f_l ] a_l,t
!..............................................................................

  type glob_con_simple_type
     ! to hold contribution of one independent function
     !  and one l to global contraction
     integer(kind=i4_kind) :: N_contributing_fcts
     integer(kind=i4_kind), pointer :: index_exp(:)
     real(kind=r8_kind), pointer :: coefs(:)
  end type glob_con_simple_type

  !------------ Declaration of constants and variables ---------------
  real(kind=r8_kind), parameter :: zero=0.0_r8_kind, one=1.0_r8_kind, &
                                   half=0.5_r8_kind, two=2.0_r8_kind
  type(unique_atom_basis_type), pointer :: basis1,basis2
! type(unique_atom_renorm_indfct_type), pointer :: &
!      renormalization1, renormalization2
  logical :: do_glob_cons1, do_glob_cons2, norm_only
! logical :: calc_renormalization_coeffs
  integer(kind=i4_kind) :: n_gc1,n_gc2,n_c1,n_c2,n_exp1,n_exp2
  type(unique_atom_glob_con_type),pointer :: glob_con1(:), glob_con2(:)
  real(kind=r8_kind),pointer :: loc_loc_int(:,:,:,:), &
       loc_glob_int(:,:,:), glob_loc_int(:,:,:), glob_glob_int(:,:)

!  integer(kind=i4_kind),pointer :: glob_loc_contrib(:), loc_glob_contrib(:)

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  subroutine fitcontract_2c_grad(fit_int_grad, index_moving_unique, i_derivative)

    !  Purpose: main routine for performing the contractions
    !           and renormaliations of the two center
    !           fit integrals.
    !           Integrals are multiplied with the approporiate
    !           fitcoefficients and summed over fit-indices,
    !           yielding a (three dimensional) force on
    !           every unique.
    !------------ Modules used ----------------------------------
    use type_module
    use iounitadmin_module, only: write_to_output_units
    implicit none
    real(kind=r8_kind), intent(inout) :: fit_int_grad(:,:,:,:,:)
    integer(kind=i4_kind),intent(in) :: index_moving_unique
    integer(kind=i4_kind),intent(in) :: i_derivative
    ! *** end of interface ***

    ! fit_int_grad:
    ! 1. dim: naexps
    ! 2. dim: nbexps
    ! 3. dim: n_indep_a
    ! 4. dim: n_indep_b
    ! 5. dim: cartesian coordinate x,y,z
    ! attention: fit_int_grad is spoiled afterwards
    ! index_moving_unique:
    ! index of moving unique atom (within all moving unique atoms)
    ! i_derivative:
    !  = 1 : [ d/dR f_k | f_l ] are passed in fit_int_grad(:,:,:,:,1:3)
    !  = 2 : [ f_k | d/dR f_l ] are passed in fit_int_grad(:,:,:,:,1:3)
    !** End of interface *****************************************

    !------------ Executable code ------------------------------------

!   renormalization1 => renormalization1_ch
!   renormalization2 => renormalization2_ch
    basis1 => ua1_basis_ch
    basis2 => ua2_basis_ch
    do_glob_cons1 = n_ch_gc1 .gt. 0
    do_glob_cons2 = n_ch_gc2 .gt. 0
    n_gc1 = n_ch_gc1
    n_gc2 = n_ch_gc2
    n_c1 = n_ch_c1
    n_c2 = n_ch_c2
    glob_con1 => ua1%glob_con_ch
    glob_con2 => ua2%glob_con_ch

    if (output_int_fitcontract) &
         call write_to_output_units &
         ("FITCONTRACT_2C_GRAD: processing integral charge overlap integral")

    n_exp1 = basis1%n_exponents
    n_exp2 = basis2%n_exponents

    if (output_int_fitcontract)  call write_to_output_units &
         ("FITCONTRACT_2C_GRAD: calling contract_2c_grad ")

    call contract_2c_grad(fit_int_grad,index_moving_unique,i_derivative)

  end subroutine fitcontract_2c_grad

  subroutine fitcontract_2c_dervs(fit_int_dervs,index_unique_a,index_unique_b)

    !  Purpose: main routine for performing the contractions
    !           and renormaliations of the two center
    !           fit integrals.
    !           Integrals are multiplied with the approporiate
    !           fitcoefficients and summed over fit-indices,
    !           yielding a (three dimensional) force on
    !           every unique.
    !------------ Modules used ----------------------------------
    use type_module
    use iounitadmin_module
    implicit none

    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(inout) :: fit_int_dervs(:,:,:,:,:,:,:)
    integer(kind=i4_kind),intent(in) :: index_unique_a,index_unique_b

    ! fit_int_grad:
    ! 1. dim: naexps
    ! 2. dim: nbexps
    ! 3. dim: n_indep_a
    ! 4. dim: n_indep_b
    ! 5. dim: cartesian coordinate x,y,z
    ! attention: fit_int_grad is spoiled afterwards
    ! index_moving_unique:
    ! index of moving unique atom (within all moving unique atoms)
    ! i_derivative:
    !  = 1 : [ d/dR f_k | f_l ] are passed in fit_int_grad(:,:,:,:,1:3)
    !  = 2 : [ f_k | d/dR f_l ] are passed in fit_int_grad(:,:,:,:,1:3)
    !** End of interface *****************************************

    !------------ Executable code ------------------------------------

!   renormalization1 => renormalization1_ch
!   renormalization2 => renormalization2_ch
    basis1 => ua1_basis_ch
    basis2 => ua2_basis_ch
    do_glob_cons1 = n_ch_gc1 .gt. 0
    do_glob_cons2 = n_ch_gc2 .gt. 0
    n_gc1 = n_ch_gc1
    n_gc2 = n_ch_gc2
    n_c1 = n_ch_c1
    n_c2 = n_ch_c2
    glob_con1 => ua1%glob_con_ch
    glob_con2 => ua2%glob_con_ch

    n_exp1 = basis1%n_exponents
    n_exp2 = basis2%n_exponents

    call contract_2c_dervs(fit_int_dervs,index_unique_a,index_unique_b)

  end subroutine fitcontract_2c_dervs





  subroutine  contract_2c_grad(fit_int_grad,index_moving_unique,i_derivative)
  use gradient_data_module, only: cpks_add_fitmat_ch_grads

    !  Purpose: perform contractions for the two
    !           center fit gradient integrals.
    !           Multiplication with fitcoefficents
    !           and summation over fit indices included.
    !           This routine yields a variable
    !           grad_fit...
    !
    !  Attention: this routines destroys the contents of
    !             the input variable 'fit_int_grad'.

    !  Subroutine called by: fitcontract_2c_grad
    !
    !  Author: FN
    !  Date: 12/96
    !

    ! output: grad_fit_ch_cartesian  to be tranformed to  gradient_totalsym

    !------------ Modules used --------------------------------------
    use type_module ! type specification parameters
    use orbitalprojection_module, only: orbitalprojection_ch, &
         orbitalprojection_globcontr_ch
    use options_module, only: options_xcmode, xcmode_model_density, &
                              options_split_gradients, options_spin_restricted
    use  calc3c_switches
#ifdef cpks_grad_fitmat
    use  cpksdervs_matrices,only: cpksalloc
#endif

    implicit none

    !------------ Declaration of formal parameters ------------------

    real(kind=r8_kind),intent(inout)        :: fit_int_grad(:,:,:,:,:)
    ! fit_int_grad(naexps,nbexps,n_indep_a,n_indep_b,3)
    ! the cartesian gradients of the primitives:
    ! either [d/dR_a f_k(a)|f_l(b)] or [f_k(a)|d/dR_b f_l(b)]

    integer(kind=i4_kind),intent(in)        :: index_moving_unique
    ! the index of the moving unique atom those gradient data is passed:
    ! unique_atoms(na)%moving_atom if d/dR_a gradients are passed and
    ! unique_atoms(nb)%moving_atom if d/dR_b gradients are passed

    integer(kind=i4_kind),intent(in)        :: i_derivative
    ! i_derivative = 1 if d/dR_a gradients are passed and
    ! i_derivative = 2 if d/dR_b gradients are passed
    !

    ! The output data is accumulated in the global arrays
    ! grad_fit_ch_cartesian(:,index_moving_unique)     [gradient_data_module]
    ! grad_mda_rhofit_cartesian(:,index_moving_unique) [gradient_data_module]
    ! grad_mda_xcpot_cartesian(:,index_moving_unique)  [gradient_data_module]



    !** End of interface *****************************************

    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: status, &
         n_uncont1,n_uncont2,i_cont1,i_cont2, &
         n_cont1,n_cont2,first_cont1,first_cont2, &
         i_ind1,i_ind2,i_exp1,i_exp2,i_c1,i_c2,i_gc1,i_gc2,&
         k_start,l_start,k_start_0,l_start_0,&
         i_gr, &
         gc1,gc2,n_independent_1,n_independent_2
    integer(kind=i4_kind)  :: l_glob_start,k_glob_start
    real(kind=r8_kind), allocatable :: intermediate(:,:), &
         intermediate_gc2(:,:), intermediate2_gc2(:,:)
    real(kind=r8_kind), pointer :: &
         contract1(:,:),contract2(:,:)
!   real(kind=r8_kind), pointer :: renorm_coeff_contr1(:), renorm_coeff_contr2(:), &
!        renorm_coeff_uncontr1(:), renorm_coeff_uncontr2(:)
    type(glob_con_simple_type), allocatable :: &
         glob_con_simple1(:), glob_con_simple2(:)
    real(kind=r8_kind), allocatable :: help_grad(:,:)
    real(kind=r8_kind) :: factor,grad_fit(3),grad_vxc(3),grad_rho(3)
#ifdef cpks_grad_fitmat
    real(kind=r8_kind),allocatable :: cpks_grad_fitmat(:,:,:)
#endif
    logical :: dra_passed, drb_passed, spin_polarized
!   integer(kind=i4_kind)::kk
    !------------ Executable code -----------------------------------

    n_uncont1 = basis1%n_uncontracted_fcts
    n_cont1 = basis1%n_contracted_fcts
    first_cont1 = n_uncont1 + 1
    n_uncont2 = basis2%n_uncontracted_fcts
    n_cont2 = basis2%n_contracted_fcts
    first_cont2 = n_uncont2 + 1
    contract1 => basis1%contractions
    contract2 => basis2%contractions
    grad_fit = zero

#ifdef cpks_grad_fitmat
    if(integralpar_cpksdervs) then
     allocate(cpks_grad_fitmat(size(coeff_charge),size(coeff_charge),3), & !!! in contract_2c_grad
                                                   stat=cpksalloc(157))
     MEMLOG(size(cpks_grad_fitmat))
     cpks_grad_fitmat=0.0_r8_kind
    endif
#endif

    mda1: if (options_xcmode() == xcmode_model_density) then
       dra_passed = i_derivative == 1
       drb_passed = i_derivative == 2
       spin_polarized = .not.options_spin_restricted()
       if (spin_polarized) then
          ! temporarily transform XCMDA coeffients b_k,up/dn into b_k,tot/spin
!         coeff_xcmda(:,1)=half*(coeff_xcmda(:,1)+coeff_xcmda(:,2)) ! b_k,tot
!         coeff_xcmda(:,2)=      coeff_xcmda(:,1)-coeff_xcmda(:,2)  ! b_k,spin
!should be already initializid, but for secure:
          coeff_xcmda_ts(:,1)=half*(coeff_xcmda(:,1)+coeff_xcmda(:,2))
          coeff_xcmda_ts(:,2)=half*(coeff_xcmda(:,1)-coeff_xcmda(:,2))
       endif

       grad_vxc = zero
       grad_rho = zero

    else mda1 ! i.e. nonmda
       dra_passed = .false.
       drb_passed = .false.
    endif mda1

    ! n_c2 -> int_data_2cff_module
    ! n_c1 ...

    allocate ( intermediate(n_c2,n_exp1),  intermediate_gc2(n_gc2,n_exp1), &
         intermediate2_gc2(n_gc2,n_c1),  stat = status )
         ASSERT(status.eq.0)

    allocate ( glob_con_simple1(n_gc1), glob_con_simple2(n_gc2), &
         stat = status )
    ASSERT(status.eq.0)

! -------------------------------------------------------------------------
! help_grad(ic1,ic2) -- working array to hold one gradient component of one
!                       pair of independent functions of a quadrupel for
!                       all basis functions (or their contributions)
!                       ic1: uncontracted primitives of ua1 +
!                            locally contracted functions of ua1 +
!                            l1-contributons to global contractions of ua1
!                       ic2: uncontracted primitives of ua2 +
!                            locally contracted functions of ua2 +
!                            l2-contributions to global contractions of ua2
! -------------------------------------------------------------------------

    allocate ( help_grad((n_c2+n_gc2),(n_c1+n_gc1)), STAT=status)
     ASSERT(status.eq.0)

    k_start_0 = orbitalprojection_ch(quadrupel%l1,quadrupel%ua1)
    l_start_0 = orbitalprojection_ch(quadrupel%l2,quadrupel%ua2)
    k_glob_start = orbitalprojection_globcontr_ch(quadrupel%ua1)
    l_glob_start = orbitalprojection_globcontr_ch(quadrupel%ua2)

    n_independent_1=ubound(fit_int_grad,3)
    n_independent_2=ubound(fit_int_grad,4)

    cartesian: do i_gr=1,3

       indep1: do i_ind1=1,n_independent_1
!         renorm_coeff_contr1 => renormalization1%renorm(i_ind1)%c_contr
!         renorm_coeff_uncontr1 => renormalization1%renorm(i_ind1)%c_exp
          k_start = k_start_0+(i_ind1-1)
          if ( do_glob_cons1 ) then
             ! the following routine maps the data for the global
             ! contraction on ONE unique from the more complicated
             ! data structures of the unique atom module into
             ! better suited ones.
             call calculate_glob_con_simple(glob_con_simple1, &
                  glob_con1, quadrupel%l1, i_ind1 )
          endif

          indep2: do i_ind2=1,n_independent_2
             l_start = l_start_0+(i_ind2-1)
!            renorm_coeff_contr2 => renormalization2%renorm(i_ind2)%c_contr
!            renorm_coeff_uncontr2 => renormalization2%renorm(i_ind2)%c_exp
             if ( do_glob_cons2 ) then
                call calculate_glob_con_simple(glob_con_simple2, &
                     glob_con2, quadrupel%l2, i_ind2 )
             endif

             ! First dimension 2 is handled, i.e. contractions are
             ! done for unique2 while the exponents of unique1
             ! are left untouched in the variables intermediate
             ! and intermediate_gc2. ---------------------------------


             ! local contractions, to be renormalized later
             ! inter(n_u2+c2,e1) := Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2)
             i_c2=first_cont2
             do i_cont2 = 1, n_cont2
                intermediate(i_c2,:) = zero
                do i_exp2 = 1, n_exp2
                   intermediate(i_c2,:) = intermediate(i_c2,:) + &
                        contract2(i_exp2,i_cont2) * &
                        fit_int_grad(:,i_exp2,i_ind1,i_ind2,i_gr)
                enddo
                i_c2=i_c2+1
             enddo

!            ! inter(u2,e1) := u_norm2(u2) * prim(e1,u2,i1,i2)
!            ! if do_glob_cons2 then
!            ! prim`(e1,e2,i1,i2) := u_norm2(e2) * prim(e1,e2,i1,i2)
!            if (do_glob_cons2) then
!               do i_exp2 = 1, n_exp2
!                  fit_int_grad(:,i_exp2,i_ind1,i_ind2,i_gr) = &
!                       fit_int_grad(:,i_exp2,i_ind1,i_ind2,i_gr) * &
!                       renorm_coeff_uncontr2(i_exp2)
!               enddo
                do i_exp2 = 1, n_uncont2
                   intermediate(i_exp2,:) = &
                        fit_int_grad(:,i_exp2,i_ind1,i_ind2,i_gr)
!               enddo
!            else
!               do i_exp2=1,n_uncont2
!                  intermediate(i_exp2,:) =  &
!                       fit_int_grad(:,i_exp2,i_ind1,i_ind2,i_gr) !!!* &
!                       renorm_coeff_uncontr2(i_exp2)
                enddo
!            endif

             ! inter(n_u2+c2,e1) := c_norm2(c2) * ...
             !                      Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2
             i_c2=first_cont2
             do i_cont2 = 1,n_cont2
!               intermediate(i_c2,:) = intermediate(i_c2,:)*renorm_coeff_contr2(i_cont2)
                i_c2=i_c2+1
             enddo

             ! inter(g2,e1) := Sum*(e2) coeff2(e2,g2) prim`(e1,e2,i1,i2)
             if (do_glob_cons2) then
                do i_gc2 = 1, n_gc2
                   intermediate_gc2(i_gc2,:) = zero
                   do i_cont2 = 1, glob_con_simple2(i_gc2)%N_contributing_fcts
                      intermediate_gc2(i_gc2,:) = &
                           intermediate_gc2(i_gc2,:) + &
                           glob_con_simple2(i_gc2)%coefs(i_cont2) * &
                           fit_int_grad(:,glob_con_simple2(i_gc2)% &
                           index_exp(i_cont2),i_ind1,i_ind2,i_gr)
                   enddo
                enddo
             endif

             ! at this point:
             ! intermediate(ic2,i_exp1) -
             ! ic2 : local contractions of unique 2, renorm_coeffs included for 2
             ! i_exp1 : ALL uncontracted exponents of unique 1
             !
             ! intermediate_gc2(igc2,i_exp1) -
             ! igc2 : global contractions of unique 2, renorm coeffs included

             ! now handle dimension 1

             ! help(u&c2,n_u1+c1) := Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
             i_c1 = first_cont1
             do i_cont1=1,n_cont1
                help_grad(1:n_c2,i_c1) = zero
                do i_exp1=1,n_exp1
                   help_grad(1:n_c2,i_c1) = &
                        help_grad(1:n_c2,i_c1) + &
                        contract1(i_exp1,i_cont1)*intermediate(:,i_exp1)
                enddo
                i_c1=i_c1+1
             enddo

             ! temp(g2,n_u1+c1) := Sum(e1) coeff1(e1,c1) inter(g2,e1)
             ! temp(g2,u1) := inter(g2,u1)
             if (do_glob_cons2) then
                i_c1=first_cont1
                do i_cont1 = 1, n_cont1
                   intermediate2_gc2(:,i_c1) = zero
                   do i_exp1 = 1, n_exp1
                      intermediate2_gc2(:,i_c1) = intermediate2_gc2(:,i_c1) + &
                           contract1(i_exp1,i_cont1)*intermediate_gc2(:,i_exp1)
                   enddo
                   i_c1=i_c1+1
                enddo

                do i_exp1=1,n_uncont1
                   intermediate2_gc2(:,i_exp1) = intermediate_gc2(:,i_exp1)
                enddo

             endif

!            ! help(u&c2,u1) := u_norm1(u1) * inter(u&c2,e1)
!            ! if do_glob_cons1 then
!            ! inter`(u&c2,e1) := u_norm1(e1) * inter(u&c2,e1)
!            if (do_glob_cons1) then
!               do i_exp1 = 1, n_exp1
!                  intermediate(:,i_exp1) = &
!                       intermediate(:,i_exp1) * &
!                       renorm_coeff_uncontr1(i_exp1)
!               enddo
                do i_exp1 = 1, n_uncont1
                   help_grad(1:n_c2,i_exp1) = &
                        intermediate(1:n_c2,i_exp1)
                enddo
!            else
!               do i_exp1=1,n_uncont1
!                  help_grad(1:n_c2,i_exp1) =  &
!                       intermediate(1:n_c2,i_exp1)* &
!                       renorm_coeff_uncontr1(i_exp1)
!               enddo
!            endif

             ! help(u&c2,n_u1+c1) := c_norm1(c1) * ...
             !                       Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
             i_c1=first_cont1
             do i_cont1 = 1, n_cont1
!               help_grad(1:n_c2,i_c1) = &
!                    help_grad(1:n_c2,i_c1)*renorm_coeff_contr1(i_cont1)
                i_c1=i_c1+1
             enddo


             ! help(n_u&c2+g2,u1) := u_norm1(u1) * inter(g2,u1)
             ! if do_glob_cons1 then
             ! inter`(g2,e1) := u_norm1(e1) * inter(g2,e1)
             if (do_glob_cons2) then
!               if (do_glob_cons1) then
!                  do i_exp1 =1, n_exp1
!                     intermediate_gc2(:,i_exp1) =  &
!                          intermediate_gc2(:,i_exp1) * &
!                          renorm_coeff_uncontr1(i_exp1)
!                  enddo
                   do i_exp1 = 1, n_uncont1
                      help_grad((n_c2+1):(n_c2+n_gc2),i_exp1) = &
                           intermediate_gc2(:,i_exp1)
                   enddo
!               else
!                  do i_exp1 = 1, n_uncont1
!                     help_grad((n_c2+1):(n_c2+n_gc2),i_exp1) = &
!                          intermediate2_gc2(:,i_exp1)* &
!                          renorm_coeff_uncontr1(i_exp1)
!                  enddo
!               endif

                ! help(n_u&c2+g2,c1) := c_norm1(c1) * ...
                !                       Sum(e1) coeff1(e1,c1) inter(g2,e1)
                i_c1=first_cont1
                do i_cont1=1,n_cont1
                   help_grad((n_c2+1):(n_c2+n_gc2),i_c1) = &
                        intermediate2_gc2(:,i_c1) !!!* &
!                       renorm_coeff_contr1(i_cont1)
                   i_c1=i_c1+1
                enddo
             endif


             if (do_glob_cons1) then
                ! help(u&c2,n_u&c1+g1) := Sum*(e1) coeff1(e1,g1) inter`u&g2,e1)
                gc1 = n_c1+1
                do i_gc1=1, n_gc1
                   help_grad(1:n_c2,gc1) = zero
                   do i_cont1 = 1, glob_con_simple1(i_gc1)%N_contributing_fcts
                      help_grad(1:n_c2,gc1) = &
                           help_grad(1:n_c2,gc1) + &
                           glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                           intermediate(:,glob_con_simple1(i_gc1)% &
                           index_exp(i_cont1))
                   enddo
                   gc1=gc1+1
                enddo

                ! help(n_u&c2+g2,n_u&c1+g1) := Sum*(e1) coeff1(e1,g1) ...
                !                              inter`(g2,e1)
                if (do_glob_cons2) then
                   gc1=n_c1+1
                   gc2=n_c2+1
                   do i_gc1=1,n_gc1
                      help_grad(gc2:(n_c2+n_gc2),gc1) = zero
                      do i_cont1 = 1,glob_con_simple1(i_gc1)%N_contributing_fcts
                         help_grad(gc2:(n_c2+n_gc2),gc1) = &
                              help_grad(gc2:(n_c2+n_gc2),gc1) + &
                              glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                              intermediate_gc2(:,glob_con_simple1(i_gc1)% &
                              index_exp(i_cont1))
                      enddo
                      gc1=gc1+1
                   enddo
                endif
             endif

             ! now multiply with coeff_charge/coeff_xcmda and sum
             ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ! Sum(k,l) a_k,tot [ d/dR f_k | f_l ] a_l,tot +
             !          a_k,tot [ f_k | d/dR f_l ] a_l,tot

             grad_fit(i_gr)=grad_fit(i_gr)+sum_up(coeff_charge,coeff_charge)

             cpks1: if(integralpar_cpksdervs) then

FPP_TIMER_START(t_cpks_grad_fitmat)
#if 1
             call cpksmat_sum_up(size(coeff_charge))
#else
             if(k_start.eq.l_start) then
               cpks_grad_fitmat(:,:,i_gr)=cpks_grad_fitmat(:,:,i_gr)+ &
                         2.0_r8_kind*cpksmat_sum_up(size(coeff_charge))
             else
               cpks_grad_fitmat(:,:,i_gr)=cpks_grad_fitmat(:,:,i_gr)+ &
                                       cpksmat_sum_up(size(coeff_charge))
             endif
#endif
FPP_TIMER_STOP(t_cpks_grad_fitmat)


             endif cpks1

             mda2: if (dra_passed) then ! model_density is used
                ! Sum(k,l,t) b_k,t [ d/dR f_k | f_l ] a_l,t

                grad_vxc(i_gr)=grad_vxc(i_gr)+sum_up(coeff_xcmda_ts(:,1),coeff_charge)
                if (spin_polarized) then
                   grad_vxc(i_gr)=grad_vxc(i_gr)+sum_up(coeff_xcmda_ts(:,2),coeff_spin)
                endif

                ! Sum(k,l,t) a_k,t [ d/dR f_k | f_l ] b_l,t

                grad_rho(i_gr)=grad_rho(i_gr)+sum_up(coeff_charge,coeff_xcmda_ts(:,1))
                if (spin_polarized) then
                   grad_rho(i_gr)=grad_rho(i_gr)+sum_up(coeff_spin,coeff_xcmda_ts(:,2))
                endif
             endif mda2

             mda3: if (drb_passed) then ! model_density is used
                ! Sum(k,l,t) b_k,t [ f_k | d/dR f_l ] a_l,t
                grad_rho(i_gr)=grad_rho(i_gr)+sum_up(coeff_xcmda_ts(:,1),coeff_charge)
                if (spin_polarized) then
                   grad_rho(i_gr)=grad_rho(i_gr)+sum_up(coeff_xcmda_ts(:,2),coeff_spin)
                endif
                ! Sum(k,l,t) a_k,t [ f_k | d/dR f_l ] b_l,t
                grad_vxc(i_gr)=grad_vxc(i_gr)+sum_up(coeff_charge,coeff_xcmda_ts(:,1))
                if (spin_polarized) then
                   grad_vxc(i_gr)=grad_vxc(i_gr)+sum_up(coeff_spin,coeff_xcmda_ts(:,2))
                endif
             endif mda3

             call free_glob_con_simple(glob_con_simple2)
          enddo indep2
          call free_glob_con_simple(glob_con_simple1)
       enddo indep1


    enddo cartesian

    deallocate(help_grad,STAT=status)
    ASSERT(status.eq.0)

    deallocate( glob_con_simple1,glob_con_simple2, STAT=status)
    ASSERT(status.eq.0)

    deallocate(intermediate,intermediate_gc2,intermediate2_gc2, &
         STAT=status)
    ASSERT(status.eq.0)

    ! add this to the total gradient in cartesian coordinates
    if (diagonal) then

       factor=0.5_r8_kind

       ![f_k|f_k] is independent of
       !the position of its nuclear center, concerning only
       !derivatives with other equivalent points? [f_k|f_k~eq]

    else
       factor=1.0_r8_kind
    endif

    nonmda: if (options_xcmode() /= xcmode_model_density) then


    grad_fit_ch_cartesian(:,index_moving_unique) = &
            grad_fit_ch_cartesian(:,index_moving_unique) + factor*grad_fit

    cpks2: if(integralpar_cpksdervs.and..not.integralpar_cpks_contribs) then
#if 0
FPP_TIMER_START(t_cpks_grad_fitmat)
      do kk=1,size(cpks_grad_fitmat,3)
      cpks_grad_fit_ch_cartesian(:,kk,index_moving_unique)= &
       cpks_grad_fit_ch_cartesian(:,kk,index_moving_unique)+ &
           factor*matmul(cpks_grad_fitmat(:,:,kk),coeff_charge)
      enddo
FPP_TIMER_STOP(t_cpks_grad_fitmat)
#endif
    endif cpks2

    else nonmda ! i.e. MDA
       if (spin_polarized) then
          ! transform XCMDA coeffients b_k,tot/spin back to b_k,up,dn
!         coeff_xcmda(:,1)=coeff_xcmda(:,1)+    coeff_xcmda(:,2) ! b_k,up
!         coeff_xcmda(:,2)=coeff_xcmda(:,1)-two*coeff_xcmda(:,2) ! b_k,dn
       endif

       if (options_split_gradients()) then
          grad_fit_ch_cartesian(:,index_moving_unique) = &
               grad_fit_ch_cartesian(:,index_moving_unique) + &
               factor*grad_fit
          grad_mda_rhofit_cartesian(:,index_moving_unique) = &
               grad_mda_rhofit_cartesian(:,index_moving_unique) + &
               factor*grad_rho
          grad_mda_xcpot_cartesian(:,index_moving_unique) = &
               grad_mda_xcpot_cartesian(:,index_moving_unique) + &
               factor*grad_vxc
       else
          grad_fit_ch_cartesian(:,index_moving_unique) = &
               grad_fit_ch_cartesian(:,index_moving_unique) + &
               factor*( grad_fit + grad_rho + grad_vxc )
       endif
    endif nonmda

#ifdef cpks_grad_fitmat
    if(integralpar_cpksdervs) then
     MEMLOG(-size(cpks_grad_fitmat))
     deallocate(cpks_grad_fitmat, stat=cpksalloc(157))
     ASSERT(cpksalloc(157).eq.0)
     cpksalloc(157)=1
    endif
#endif

  contains !of contract_2c_grad

    real(kind=r8_kind) function sum_up(coeff1,coeff2)
      !----------------------------------------------------------------
      !  Purpose: Returns sums of the form Sum(k,l) a_k F_kl b_l
      !           with F_kl being any matrix in the charge fit basis
      !           and a_k and b_l being any type of charge fit coefficients
      !           UB (9/97)
      !  Input:
      !  n_c1, n_gc1 : dimension of the summation indices split into local
      !  n_c2, n_gc2   contributions (1:n_c1 and 1:n_c2) and global contri-
      !                butions (n_c1+1:n_c1+n_gc1 and n_c2+1:n_c2+n_gc2)
      !  help_grad(:,:) : the matrix elements F_kl stored in transposed fashion
      !                   F_kl = help_grad(l,k)
      !  k_start, k_glob_start : pointer on the first element of the local
      !  l_start, l_glob_start   (k_start and l_start) and the global
      !                          (k_glob_start and l_glob_start) contributions
      !                          of the coefficient vectors a_k and b_l
      !  n_independent_1 : storage increment of the local contributions
      !  n_independent_2   to a_k: a_k = a(k_start + (k-1)*n_independent_1) and
      !                    to b_l: b_l = b(l_start + (l-1)*n_independent_2)
      !----------------------------------------------------------------
      implicit none
      !------------ Declaration of formal parameters ------------------
      real(kind=r8_kind), intent(in) :: coeff1(:), coeff2(:)
      ! the coefficient vectors a_k and b_l in linear packed storage mode
      !** End of interface *****************************************
      !------------ Declaration of local variables --------------------
      integer(kind=i4_kind) :: i_a,i_b,k,l!,nt
      real   (kind=r8_kind) :: temp!,temp1 ! to hold Sum(l) F_kl b_l
      !------------ Executable code -----------------------------------
!      nt=0
!      temp1 = 0.0_r8_kind
      sum_up = 0.0_r8_kind
      ! < a | F | b >
      k = k_start
      do i_a=1,n_c1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
!            nt=nt+1
!            temp1 = temp1 + help_grad(i_b,i_a) * coeff2(l) * coeff1(k)
!            print*,help_grad(i_b,i_a),coeff2(l),coeff1(k),nt,help_grad(i_b,i_a)*coeff2(l)*coeff1(k),temp1,l,k
            l = l + n_independent_2
         enddo
         sum_up = sum_up + coeff1(k) * temp
         k = k + n_independent_1
      enddo

!      if(abs(sum_up).gt.1.e-10) &
!                print*,k_start,l_start,sum_up

      ! < ag | F | b >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + n_independent_2
         end do
         sum_up = sum_up + coeff1(k) * temp
         k = k + 1
      end do
      if (n_gc2 == 0) return
      ! < a | F | bg >
      k = k_start
      do i_a=1,n_c1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up = sum_up + coeff1(k) * temp
         k = k + n_independent_1
      end do
      ! < ag | F | bg >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up = sum_up + coeff1(k) * temp
         k = k + 1
      end do
    end function sum_up

#if 1
     subroutine cpksmat_sum_up(n_ch)
      !----------------------------------------------------------------
      !  Purpose: Returns sums of the form Sum(k,l) a_k F_kl b_l
      !           with F_kl being any matrix in the charge fit basis
      !           and a_k and b_l being any type of charge fit coefficients
      !           UB (9/97)
      !  Input:
      !  n_c1, n_gc1 : dimension of the summation indices split into local
      !  n_c2, n_gc2   contributions (1:n_c1 and 1:n_c2) and global contri-
      !                butions (n_c1+1:n_c1+n_gc1 and n_c2+1:n_c2+n_gc2)
      !  help_grad(:,:) : the matrix elements F_kl stored in transposed fashion
      !                   F_kl = help_grad(l,k)
      !  k_start, k_glob_start : pointer on the first element of the local
      !  l_start, l_glob_start   (k_start and l_start) and the global
      !                          (k_glob_start and l_glob_start) contributions
      !                          of the coefficient vectors a_k and b_l
      !  n_independent_1 : storage increment of the local contributions
      !  n_independent_2   to a_k: a_k = a(k_start + (k-1)*n_independent_1) and
      !                    to b_l: b_l = b(l_start + (l-1)*n_independent_2)
      !----------------------------------------------------------------
      implicit none
      !------------ Declaration of formal parameters ------------------
      integer(kind=i4_kind),intent(in):: n_ch
!     real(kind=r8_kind) :: sum_up(n_ch,n_ch)
      ! the coefficient vectors a_k and b_l in linear packed storage mode
      !** End of interface *****************************************
      !------------ Declaration of local variables --------------------
      integer(kind=i4_kind) :: i_a,i_b,k,l
      real   (kind=r8_kind) :: fact ! to hold Sum(l) F_kl b_l
      !------------ Executable code -----------------------------------
!      sum_up=0.0_r8_kind
      if(diagonal) then
      fact=0.5_r8_kind
      else
      fact=1.0_r8_kind
      endif
      k = k_start
      do i_a=1,n_c1
         l = l_start
         do i_b=1,n_c2
#if 0
            sum_up(l,k) =  help_grad(i_b,i_a)
            sum_up(k,l) =  help_grad(i_b,i_a)
            cpks_grad_fitmat(k,l,i_gr)=cpks_grad_fitmat(k,l,i_gr)+help_grad(i_b,i_a)
            cpks_grad_fitmat(l,k,i_gr)=cpks_grad_fitmat(l,k,i_gr)+help_grad(i_b,i_a)
#else
      cpks_grad_fit_ch_cartesian(k,i_gr,index_moving_unique)= &
      cpks_grad_fit_ch_cartesian(k,i_gr,index_moving_unique)+ &
           fact*coeff_charge(l)*help_grad(i_b,i_a)
      cpks_grad_fit_ch_cartesian(l,i_gr,index_moving_unique)= &
      cpks_grad_fit_ch_cartesian(l,i_gr,index_moving_unique)+ &
           fact*coeff_charge(k)*help_grad(i_b,i_a)
#endif
            l = l + n_independent_2
         end do
         k = k + n_independent_1
      end do


#if 0
      print*,' < a | F | b >', &
       l_start,':',l_start+n_independent_2*(n_c2-1), &
       k_start,':',k_start+n_independent_1*(n_c1-1), &
       sum(abs(sum_up)), sum(abs(sum_up(k_start:k_start+n_independent_1*(n_c1-1), &
                                        l_start:l_start+n_independent_2*(n_c2-1)) ))+ &
                         sum(abs(sum_up(l_start:l_start+n_independent_2*(n_c2-1), &
                                        k_start:k_start+n_independent_1*(n_c1-1)) ))
#endif
ASSERT(n_gc1.eq.0)
ASSERT(n_gc2.eq.0)
    end subroutine cpksmat_sum_up

#else
     function cpksmat_sum_up(n_ch) result(sum_up)
      !----------------------------------------------------------------
      !  Purpose: Returns sums of the form Sum(k,l) a_k F_kl b_l
      !           with F_kl being any matrix in the charge fit basis
      !           and a_k and b_l being any type of charge fit coefficients
      !           UB (9/97)
      !  Input:
      !  n_c1, n_gc1 : dimension of the summation indices split into local
      !  n_c2, n_gc2   contributions (1:n_c1 and 1:n_c2) and global contri-
      !                butions (n_c1+1:n_c1+n_gc1 and n_c2+1:n_c2+n_gc2)
      !  help_grad(:,:) : the matrix elements F_kl stored in transposed fashion
      !                   F_kl = help_grad(l,k)
      !  k_start, k_glob_start : pointer on the first element of the local
      !  l_start, l_glob_start   (k_start and l_start) and the global
      !                          (k_glob_start and l_glob_start) contributions
      !                          of the coefficient vectors a_k and b_l
      !  n_independent_1 : storage increment of the local contributions
      !  n_independent_2   to a_k: a_k = a(k_start + (k-1)*n_independent_1) and
      !                    to b_l: b_l = b(l_start + (l-1)*n_independent_2)
      !----------------------------------------------------------------
      implicit none
      !------------ Declaration of formal parameters ------------------
      integer(kind=i4_kind),intent(in):: n_ch
      real(kind=r8_kind) :: sum_up(n_ch,n_ch)
      ! the coefficient vectors a_k and b_l in linear packed storage mode
      !** End of interface *****************************************
      !------------ Declaration of local variables --------------------
      integer(kind=i4_kind) :: i_a,i_b,k,l
      real   (kind=r8_kind) :: temp ! to hold Sum(l) F_kl b_l
      !------------ Executable code -----------------------------------
      sum_up(:,:) = 0.0_r8_kind
      k = k_start
      do i_a=1,n_c1
         l = l_start
         do i_b=1,n_c2
            sum_up(l,k) =  help_grad(i_b,i_a)
            sum_up(k,l) =  help_grad(i_b,i_a)
            l = l + n_independent_2
         end do
         k = k + n_independent_1
      end do
      print*,' < a | F | b >', &
       l_start,':',l_start+n_independent_2*(n_c2-1), &
       k_start,':',k_start+n_independent_1*(n_c1-1), &
       sum(abs(sum_up))

ASSERT(n_gc1.eq.0)
ASSERT(n_gc2.eq.0)
#if 0 /* no_glob_contractions */
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_start
         do i_b=1,n_c2
            sum_up(l,k) = sum_up(l,k) + help_grad(i_b,i_a)
            l = l + n_independent_2
         end do
         k = k + 1
      end do


      if (n_gc2 == 0) return
      ! < a | F | bg >
      k = k_start
      do i_a=1,n_c1
         l = l_glob_start
         do i_b=n_c2+1,n_c2+n_gc2
            sum_up(l,k) = sum_up(l,k) + help_grad(i_b,i_a)
            l = l + 1
         end do
         k = k + n_independent_1
      end do
      ! < ag | F | bg >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_glob_start
         do i_b=n_c2+1,n_c2+n_gc2
            sum_up(l,k) = sum_up(l,k) + help_grad(i_b,i_a)
            l = l + 1
         end do
         k = k + 1
      end do
#endif
    end function cpksmat_sum_up
#endif

#if 0
     function cpks_sum_up(coeff2,n_ch) result(sum_up)
      !----------------------------------------------------------------
      !  Purpose: Returns sums of the form Sum(k,l) a_k F_kl b_l
      !           with F_kl being any matrix in the charge fit basis
      !           and a_k and b_l being any type of charge fit coefficients
      !           UB (9/97)
      !  Input:
      !  n_c1, n_gc1 : dimension of the summation indices split into local
      !  n_c2, n_gc2   contributions (1:n_c1 and 1:n_c2) and global contri-
      !                butions (n_c1+1:n_c1+n_gc1 and n_c2+1:n_c2+n_gc2)
      !  help_grad(:,:) : the matrix elements F_kl stored in transposed fashion
      !                   F_kl = help_grad(l,k)
      !  k_start, k_glob_start : pointer on the first element of the local
      !  l_start, l_glob_start   (k_start and l_start) and the global
      !                          (k_glob_start and l_glob_start) contributions
      !                          of the coefficient vectors a_k and b_l
      !  n_independent_1 : storage increment of the local contributions
      !  n_independent_2   to a_k: a_k = a(k_start + (k-1)*n_independent_1) and
      !                    to b_l: b_l = b(l_start + (l-1)*n_independent_2)
      !----------------------------------------------------------------
      implicit none
      !------------ Declaration of formal parameters ------------------
      real(kind=r8_kind), intent(in) ::  coeff2(:)
      integer(kind=i4_kind),intent(in):: n_ch
      real(kind=r8_kind) :: sum_up(n_ch)
      ! the coefficient vectors a_k and b_l in linear packed storage mode
      !** End of interface *****************************************
      !------------ Declaration of local variables --------------------
      integer(kind=i4_kind) :: i_a,i_b,k,l
      real   (kind=r8_kind) :: temp ! to hold Sum(l) F_kl b_l
      !------------ Executable code -----------------------------------
      sum_up(:) = 0.0_r8_kind
      ! < a | F | b >
      k = k_start
      do i_a=1,n_c1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + n_independent_2
         end do
         sum_up(k) = sum_up(k) + temp
         k = k + n_independent_1
      end do
      ! < ag | F | b >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + n_independent_2
         end do
         sum_up(k) = sum_up(k) + temp
         k = k + 1
      end do
      if (n_gc2 == 0) return
      ! < a | F | bg >
      k = k_start
      do i_a=1,n_c1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up(k) = sum_up(k) + temp
         k = k + n_independent_1
      end do
      ! < ag | F | bg >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up(k) = sum_up(k) + temp
         k = k + 1
      end do
    end function cpks_sum_up
#endif

  end subroutine contract_2c_grad

  subroutine  contract_2c_dervs(fit_int_dervs,unique_1,unique_2)

    !  Purpose: perform contractions for the two
    !           center fit gradient integrals.
    !           Multiplication with fitcoefficents
    !           and summation over fit indices included.
    !           This routine yields a variable
    !           grad_fit...
    !
    !  Attention: this routines destroys the contents of
    !             the input variable 'fit_int_dervs'.

    !  Subroutine called by: fitcontract_2c_grad
    !
    !------------ Modules used --------------------------------------
    use type_module ! type specification parameters
    use orbitalprojection_module, only: orbitalprojection_ch, &
         orbitalprojection_globcontr_ch
    use options_module, only: options_xcmode, xcmode_model_density, &
                              options_split_gradients, options_spin_restricted
    use  calc3c_switches
    use  gradient_data_module, only: dervs_totalsym,gradient_index

    implicit none

    !------------ Declaration of formal parameters ------------------

    real(kind=r8_kind),intent(inout)        ::   fit_int_dervs(:,:,:,:,:,:,:)
    integer(kind=i4_kind),intent(in)  :: unique_1,unique_2

    !

    ! The output data is accumulated in the global arrays
    ! grad_fit_ch_cartesian(:,index_moving_unique)     [gradient_data_module]
    ! grad_mda_rhofit_cartesian(:,index_moving_unique) [gradient_data_module]
    ! grad_mda_xcpot_cartesian(:,index_moving_unique)  [gradient_data_module]

    !** End of interface *****************************************

    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: n_equal_atoms
    integer(kind=i4_kind) :: status, &
         n_uncont1,n_uncont2,i_cont1,i_cont2, &
         n_cont1,n_cont2,first_cont1,first_cont2, &
         i_ind1,i_ind2,i_exp1,i_exp2,i_c1,i_c2,i_gc1,i_gc2,&
         k_start,l_start,k_start_0,l_start_0,&
         i_gr,k2dr, &
         gc1,gc2,n_independent_1,n_independent_2
    integer(kind=i4_kind)  :: l_glob_start,k_glob_start
    integer(kind=i4_kind)  :: equal_atom
    real(kind=r8_kind), allocatable :: intermediate(:,:), &
         intermediate_gc2(:,:), intermediate2_gc2(:,:)
    real(kind=r8_kind), pointer :: &
         contract1(:,:),contract2(:,:)
    type(glob_con_simple_type), allocatable :: &
         glob_con_simple1(:), glob_con_simple2(:)
    real(kind=r8_kind), allocatable :: help_grad(:,:)
    real(kind=r8_kind) :: factor
    real(kind=r8_kind),allocatable :: dervs_fit(:,:,:)
    !------------ Executable code -----------------------------------

!    print*,quadrupel%ua1,quadrupel%l1,'///',quadrupel%ua2,quadrupel%l2

    n_uncont1 = basis1%n_uncontracted_fcts
    n_cont1 = basis1%n_contracted_fcts
    first_cont1 = n_uncont1 + 1
    n_uncont2 = basis2%n_uncontracted_fcts
    n_cont2 = basis2%n_contracted_fcts
    first_cont2 = n_uncont2 + 1
    contract1 => basis1%contractions
    contract2 => basis2%contractions



    ! n_c2 -> int_data_2cff_module
    ! n_c1 ...

    allocate ( intermediate(n_c2,n_exp1), &
         intermediate_gc2(n_gc2,n_exp1), &
         intermediate2_gc2(n_gc2,n_c1), &
         stat = status )
         ASSERT(status.eq.0)

    allocate ( glob_con_simple1(n_gc1), glob_con_simple2(n_gc2), &
         stat = status )
    ASSERT(status.eq.0)

! -------------------------------------------------------------------------
! help_grad(ic1,ic2) -- working array to hold one gradient component of one
!                       pair of independent functions of a quadrupel for
!                       all basis functions (or their contributions)
!                       ic1: uncontracted primitives of ua1 +
!                            locally contracted functions of ua1 +
!                            l1-contributons to global contractions of ua1
!                       ic2: uncontracted primitives of ua2 +
!                            locally contracted functions of ua2 +
!                            l2-contributions to global contractions of ua2
! -------------------------------------------------------------------------

    allocate ( help_grad((n_c2+n_gc2),(n_c1+n_gc1)), STAT=status)
     ASSERT(status.eq.0)

    k_start_0 = orbitalprojection_ch(quadrupel%l1,quadrupel%ua1)
    l_start_0 = orbitalprojection_ch(quadrupel%l2,quadrupel%ua2)
    k_glob_start = orbitalprojection_globcontr_ch(quadrupel%ua1)
    l_glob_start = orbitalprojection_globcontr_ch(quadrupel%ua2)

    n_independent_1=ubound(fit_int_dervs,3)
    n_independent_2=ubound(fit_int_dervs,4)
    n_equal_atoms=size(fit_int_dervs,7)
    allocate(dervs_fit(3,3,n_equal_atoms),stat=status)
    ASSERT(status.eq.0)
    dervs_fit = zero

    cartesian_gr: do i_gr=1,3
               equa: do equal_atom=1,n_equal_atoms

               cartesian_2dr:do k2dr=1,3

       indep1: do i_ind1=1,n_independent_1
!          renorm_coeff_contr1 => renormalization1%renorm(i_ind1)%c_contr
!          renorm_coeff_uncontr1 => renormalization1%renorm(i_ind1)%c_exp
          k_start = k_start_0+(i_ind1-1)
          if ( do_glob_cons1 ) then
             ! the following routine maps the data for the global
             ! contraction on ONE unique from the more complicated
             ! data structures of the unique atom module into
             ! better suited ones.
             call calculate_glob_con_simple(glob_con_simple1, &
                  glob_con1, quadrupel%l1, i_ind1 )
          endif

          indep2: do i_ind2=1,n_independent_2
             l_start = l_start_0+(i_ind2-1)
!             renorm_coeff_contr2 => renormalization2%renorm(i_ind2)%c_contr
!             renorm_coeff_uncontr2 => renormalization2%renorm(i_ind2)%c_exp
             if ( do_glob_cons2 ) then
                call calculate_glob_con_simple(glob_con_simple2, &
                                               glob_con2, quadrupel%l2, i_ind2 )
             endif

             ! First dimension 2 is handled, i.e. contractions are
             ! done for unique2 while the exponents of unique1
             ! are left untouched in the variables intermediate
             ! and intermediate_gc2. ---------------------------------


             ! local contractions, to be renormalized later
             ! inter(n_u2+c2,e1) := Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2)
             i_c2=first_cont2
             do i_cont2 = 1, n_cont2
                intermediate(i_c2,:) = zero
                do i_exp2 = 1, n_exp2
                   intermediate(i_c2,:) = intermediate(i_c2,:) + &
                        contract2(i_exp2,i_cont2) * &
                        fit_int_dervs(:,i_exp2,i_ind1,i_ind2,i_gr,k2dr,equal_atom)
                enddo
                i_c2=i_c2+1
             enddo

!             ! inter(u2,e1) := u_norm2(u2) * prim(e1,u2,i1,i2)
!             ! if do_glob_cons2 then
!             ! prim`(e1,e2,i1,i2) := u_norm2(e2) * prim(e1,e2,i1,i2)
!             if (do_glob_cons2) then
!                do i_exp2 = 1, n_exp2
!                   fit_int_dervs(:,i_exp2,i_ind1,i_ind2,i_gr,k2dr,equal_atom) = &
!                        fit_int_dervs(:,i_exp2,i_ind1,i_ind2,i_gr,k2dr,equal_atom) * &
!                        renorm_coeff_uncontr2(i_exp2)
!                enddo
                do i_exp2 = 1, n_uncont2
                   intermediate(i_exp2,:) = &
                        fit_int_dervs(:,i_exp2,i_ind1,i_ind2,i_gr,k2dr,equal_atom)
                enddo
!             else
!                do i_exp2=1,n_uncont2
!                   intermediate(i_exp2,:) =  &
!                        fit_int_dervs(:,i_exp2,i_ind1,i_ind2,i_gr,k2dr,equal_atom) * &
!                        renorm_coeff_uncontr2(i_exp2)
!                enddo
!             endif

             ! inter(n_u2+c2,e1) := c_norm2(c2) * ...
             !                      Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2
             i_c2=first_cont2
             do i_cont2 = 1,n_cont2
!                intermediate(i_c2,:) = intermediate(i_c2,:)*renorm_coeff_contr2(i_cont2)
                i_c2=i_c2+1
             enddo

             ! inter(g2,e1) := Sum*(e2) coeff2(e2,g2) prim`(e1,e2,i1,i2)
             if (do_glob_cons2) then
                do i_gc2 = 1, n_gc2
                   intermediate_gc2(i_gc2,:) = zero
                   do i_cont2 = 1, glob_con_simple2(i_gc2)%N_contributing_fcts
                      intermediate_gc2(i_gc2,:) = &
                           intermediate_gc2(i_gc2,:) + &
                           glob_con_simple2(i_gc2)%coefs(i_cont2) * &
                       fit_int_dervs( :,glob_con_simple2(i_gc2)%index_exp(i_cont2), &
                                      i_ind1,i_ind2,i_gr,k2dr,equal_atom)
                   enddo
                enddo
             endif

             ! at this point:
             ! intermediate(ic2,i_exp1) -
             ! ic2 : local contractions of unique 2, renorm_coeffs included for 2
             ! i_exp1 : ALL uncontracted exponents of unique 1
             !
             ! intermediate_gc2(igc2,i_exp1) -
             ! igc2 : global contractions of unique 2, renorm coeffs included

             ! now handle dimension 1

             ! help(u&c2,n_u1+c1) := Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
             i_c1 = first_cont1
             do i_cont1=1,n_cont1
                help_grad(1:n_c2,i_c1) = zero
                do i_exp1=1,n_exp1
                   help_grad(1:n_c2,i_c1) = &
                        help_grad(1:n_c2,i_c1) + &
                        contract1(i_exp1,i_cont1)*intermediate(:,i_exp1)
                enddo
                i_c1=i_c1+1
             enddo

             ! temp(g2,n_u1+c1) := Sum(e1) coeff1(e1,c1) inter(g2,e1)
             ! temp(g2,u1) := inter(g2,u1)
             if (do_glob_cons2) then
                i_c1=first_cont1
                do i_cont1 = 1, n_cont1
                   intermediate2_gc2(:,i_c1) = zero
                   do i_exp1 = 1, n_exp1
                      intermediate2_gc2(:,i_c1) = intermediate2_gc2(:,i_c1) + &
                           contract1(i_exp1,i_cont1)*intermediate_gc2(:,i_exp1)
                   enddo
                   i_c1=i_c1+1
                enddo

                do i_exp1=1,n_uncont1
                   intermediate2_gc2(:,i_exp1) = intermediate_gc2(:,i_exp1)
                enddo

             endif

!             ! help(u&c2,u1) := u_norm1(u1) * inter(u&c2,e1)
!             ! if do_glob_cons1 then
!             ! inter`(u&c2,e1) := u_norm1(e1) * inter(u&c2,e1)
!             if (do_glob_cons1) then
!                do i_exp1 = 1, n_exp1
!                   intermediate(:,i_exp1) = &
!                        intermediate(:,i_exp1) * &
!                        renorm_coeff_uncontr1(i_exp1)
!                enddo
                do i_exp1 = 1, n_uncont1
                   help_grad(1:n_c2,i_exp1) = &
                        intermediate(1:n_c2,i_exp1)
                enddo
!             else
!                do i_exp1=1,n_uncont1
!                   help_grad(1:n_c2,i_exp1) =  &
!                        intermediate(1:n_c2,i_exp1)* &
!                        renorm_coeff_uncontr1(i_exp1)
!                enddo
!             endif

             ! help(u&c2,n_u1+c1) := c_norm1(c1) * ...
             !                       Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
             i_c1=first_cont1
             do i_cont1 = 1, n_cont1
!                help_grad(1:n_c2,i_c1) = &
!                     help_grad(1:n_c2,i_c1)*renorm_coeff_contr1(i_cont1)
                i_c1=i_c1+1
             enddo


             ! help(n_u&c2+g2,u1) := u_norm1(u1) * inter(g2,u1)
             ! if do_glob_cons1 then
             ! inter`(g2,e1) := u_norm1(e1) * inter(g2,e1)
             if (do_glob_cons2) then
!                if (do_glob_cons1) then
!                   do i_exp1 =1, n_exp1
!                      intermediate_gc2(:,i_exp1) =  &
!                           intermediate_gc2(:,i_exp1) * &
!                           renorm_coeff_uncontr1(i_exp1)
!                   enddo
                   do i_exp1 = 1, n_uncont1
                      help_grad((n_c2+1):(n_c2+n_gc2),i_exp1) = &
                           intermediate_gc2(:,i_exp1)
                   enddo
!                else
!                   do i_exp1 = 1, n_uncont1
!                      help_grad((n_c2+1):(n_c2+n_gc2),i_exp1) = &
!                           intermediate2_gc2(:,i_exp1)* &
!                           renorm_coeff_uncontr1(i_exp1)
!                   enddo
!                endif

                ! help(n_u&c2+g2,c1) := c_norm1(c1) * ...
                !                       Sum(e1) coeff1(e1,c1) inter(g2,e1)
                i_c1=first_cont1
                do i_cont1=1,n_cont1
                   help_grad((n_c2+1):(n_c2+n_gc2),i_c1) = &
                        intermediate2_gc2(:,i_c1)!!!*renorm_coeff_contr1(i_cont1)
                   i_c1=i_c1+1
                enddo
             endif


             if (do_glob_cons1) then
                ! help(u&c2,n_u&c1+g1) := Sum*(e1) coeff1(e1,g1) inter`u&g2,e1)
                gc1 = n_c1+1
                do i_gc1=1, n_gc1
                   help_grad(1:n_c2,gc1) = zero
                   do i_cont1 = 1, glob_con_simple1(i_gc1)%N_contributing_fcts
                      help_grad(1:n_c2,gc1) = &
                           help_grad(1:n_c2,gc1) + &
                           glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                           intermediate(:,glob_con_simple1(i_gc1)% &
                           index_exp(i_cont1))
                   enddo
                   gc1=gc1+1
                enddo

                ! help(n_u&c2+g2,n_u&c1+g1) := Sum*(e1) coeff1(e1,g1) ...
                !                              inter`(g2,e1)
                if (do_glob_cons2) then
                   gc1=n_c1+1
                   gc2=n_c2+1
                   do i_gc1=1,n_gc1
                      help_grad(gc2:(n_c2+n_gc2),gc1) = zero
                      do i_cont1 = 1,glob_con_simple1(i_gc1)%N_contributing_fcts
                         help_grad(gc2:(n_c2+n_gc2),gc1) = &
                              help_grad(gc2:(n_c2+n_gc2),gc1) + &
                              glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                              intermediate_gc2(:,glob_con_simple1(i_gc1)% &
                              index_exp(i_cont1))
                      enddo
                      gc1=gc1+1
                   enddo
                endif
             endif

             ! now multiply with coeff_charge/coeff_xcmda and sum
             ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ! Sum(k,l) a_k,tot [ d/dR f_k | f_l ] a_l,tot +
             !          a_k,tot [ f_k | d/dR f_l ] a_l,tot

             dervs_fit(i_gr,k2dr,equal_atom)=dervs_fit(i_gr,k2dr,equal_atom) &
                                      +sum_up_dervs(coeff_charge,coeff_charge)


             call free_glob_con_simple(glob_con_simple2)
          enddo indep2
          call free_glob_con_simple(glob_con_simple1)
       enddo indep1


    enddo cartesian_2dr
    enddo equa
    enddo cartesian_gr


    deallocate(help_grad,STAT=status)
    ASSERT(status.eq.0)

    deallocate( glob_con_simple1,glob_con_simple2, STAT=status)
    ASSERT(status.eq.0)

    deallocate(intermediate,intermediate_gc2,intermediate2_gc2, &
         STAT=status)
    ASSERT(status.eq.0)

    ! add this to the total gradient in cartesian coordinates
    if (diagonal) then

       factor=0.25_r8_kind

       ![f_k|f_k] is independent of
       !the position of its nuclear center, concerning only
       !derivatives with other equivalent points? [f_k|f_k~eq]

    else
       factor=0.5_r8_kind
    endif

    dervs_fit =  factor*dervs_fit

!     do equal_atom=2,n_equal_atoms
!      dervs_fit(:,:,equal_atom)=0.0_r8_kind !!! dervs_fit(:,:,equal_atom)*2.0_r8_kind
!     enddo
!    if(unique_1.ne.unique_2.or.unique_1.eq.2) then
!     do equal_atom=1,n_equal_atoms
!      dervs_fit(:,:,equal_atom)=0.0_r8_kind
!     enddo
!    endif



!   do  equal_atom=1,size(
!  if(unique_1.eq.1) &
!  print*, sum(dervs_fit(:,3,:)),'sum dervs_fit',unique_1,unique_2,size(dervs_fit,3)
#define IMPL_CONTRIBS_TO_DERVS
#ifdef IMPL_CONTRIBS_TO_DERVS
   call addto_totsym(dervs_totalsym,dervs_fit,unique_1,unique_2) !!!(1)
#endif
!  print*, 'now fit contrib to dervs_totalsym on'

!  call addto_totsym(dervs_totalsym_fit,dervs_fit,unique_1,unique_2) ! for test only
!   this is only to see dervs_totalsym_fit contribs separate
!   to be printed in main_grad

   deallocate(dervs_fit,stat=status)
   ASSERT(status.eq.0)


  contains

   subroutine addto_totsym(dervs_totalsym,dervs_fit,unique_1,unique_2)

   real(kind=r8_kind), intent(inout):: dervs_totalsym(:,:)
   real(kind=r8_kind), intent(in):: dervs_fit(:,:,:)
   integer(kind=i4_kind), intent(in)::unique_1,unique_2

   real(kind=r8_kind),allocatable :: dervs_totsym_temp1(:,:,:), &
                                     dervs_totsym_temp2(:,:,:)
   real(kind=r8_kind)    :: weight

   integer(kind=i4_kind):: center_1,index_1,grad_dim_1, &
                           center_2,index_2,grad_dim_2, &
                           i,i_gr,alloc_stat,equal_atom,n_equal_atoms

   ! each ab contrib results in four 3x3 blocks of deriv matrix
!  print*
!  print*,'dervs_fit addto_totsym unique_1,unique_2',unique_1,unique_2
!  if(size(dervs_fit,3).gt.1) then
!   print*,dervs_fit(1,:,1),dervs_fit(1,:,2)
!   print*,dervs_fit(2,:,1),dervs_fit(2,:,2)
!   print*,dervs_fit(3,:,1),dervs_fit(3,:,2)
!  else
!   print*,dervs_fit(1,:,1)
!   print*,dervs_fit(2,:,1)
!   print*,dervs_fit(3,:,1)
!  endif



    center_1=moving_unique_atom_index(unique_1)
    weight=unique_atoms(center_1)%n_equal_atoms
    index_1=gradient_index(unique_1)-1
    grad_dim_1=gradient_index(unique_1+1)-gradient_index(unique_1)

    center_2=moving_unique_atom_index(unique_2)
    n_equal_atoms=unique_atoms(center_2)%n_equal_atoms
    index_2=gradient_index(unique_2)-1
    grad_dim_2=gradient_index(unique_2+1)-gradient_index(unique_2)

    allocate( dervs_totsym_temp1(grad_dim_1,3,size(dervs_fit,3)), &
              stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    dervs_totsym_temp1=0.0_r8_kind

    do i=1,grad_dim_1
     do i_gr=1,3
     dervs_totsym_temp1(i,:,:) = dervs_totsym_temp1(i,:,:) + &
     weight * unique_atom_grad_info(unique_1)%m(i,i_gr,1)* &
              dervs_fit(i_gr,:,:)
     enddo
    enddo

    do i=1,grad_dim_1
    do equal_atom=1,n_equal_atoms
     do k2dr=1,3
     dervs_totalsym(index_1+1:index_1+grad_dim_1,index_1+i) = &
     dervs_totalsym(index_1+1:index_1+grad_dim_1,index_1+i) - &  !(1)
     unique_atom_grad_info(unique_1)%m(i,k2dr,1)* &
                        dervs_totsym_temp1(:,k2dr,equal_atom)
     enddo
    enddo
    enddo

!  print*,'dervs_totalsym 1x1',index_1+1,index_1+grad_dim_1,index_1+1,index_1+grad_dim_1
!   do i=1,grad_dim_1
!    print*,dervs_totalsym(index_1+1:index_1+grad_dim_1,index_1+i)
!   enddo

    do i=1,grad_dim_2
    do equal_atom=1,n_equal_atoms
     do k2dr=1,3
     dervs_totalsym(index_1+1:index_1+grad_dim_1,index_2+i) = &
     dervs_totalsym(index_1+1:index_1+grad_dim_1,index_2+i) + &   !(2)
     unique_atom_grad_info(unique_2)%m(i,k2dr,equal_atom)* &
                        dervs_totsym_temp1(:,k2dr,equal_atom)
     enddo
    enddo
    enddo
! print*,'dervs_totalsym 1x2',index_1+1,index_1+grad_dim_1,index_2+1,index_2+grad_dim_2
! do i=1,grad_dim_2
!   print*,dervs_totalsym(index_1+1:index_1+grad_dim_1,index_2+i)
! enddo

   deallocate(dervs_totsym_temp1, stat=alloc_stat)
            ASSERT(alloc_stat.eq.0)
   !-------------------------------------------------------------






   !2nd ind   -->  1st proj
   !          -->  2nd proj

    allocate(dervs_totsym_temp2(grad_dim_2,3,n_equal_atoms),stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    dervs_totsym_temp2=0.0_r8_kind

    do i=1,grad_dim_2
    do equal_atom=1,n_equal_atoms
     do i_gr=1,3
     dervs_totsym_temp2(i,:,equal_atom) = dervs_totsym_temp2(i,:,equal_atom) + &
      unique_atom_grad_info(unique_2)%m(i,i_gr,equal_atom)* &
                dervs_fit(i_gr,:,equal_atom)*weight
     enddo
    enddo
    enddo

!  print*,'dervs_totsym_temp2 for check'
!  do i=1,grad_dim_2
!   print*,dervs_totsym_temp2(i,:,1)
!  enddo

    do i=1,grad_dim_1
    do equal_atom=1,n_equal_atoms
     do k2dr=1,3
     dervs_totalsym(index_2+1:index_2+grad_dim_2,index_1+i) = &
     dervs_totalsym(index_2+1:index_2+grad_dim_2,index_1+i) + &   !(3)
     unique_atom_grad_info(unique_1)%m(i,k2dr,1)* &
                        dervs_totsym_temp2(:,k2dr,equal_atom)
     enddo
     enddo
    enddo
! print*,'dervs_totalsym 2x1',index_2+1,index_2+grad_dim_2,index_1+1,index_1+grad_dim_1
! do i=1,grad_dim_1
!   print*,dervs_totalsym(index_2+1:index_2+grad_dim_2,index_1+i)
! enddo

    do i=1,grad_dim_2
     do equal_atom=1,n_equal_atoms
     do k2dr=1,3
     dervs_totalsym(index_2+1:index_2+grad_dim_2,index_2+i) = &
     dervs_totalsym(index_2+1:index_2+grad_dim_2,index_2+i) - &  !(4)
     unique_atom_grad_info(unique_2)%m(i,k2dr,equal_atom)* &
                        dervs_totsym_temp2(:,k2dr,equal_atom)
     enddo
     enddo
    enddo
!  print*,'dervs_totalsym 2x2',index_2+1,index_2+grad_dim_2,index_2+1,index_2+grad_dim_2, &
!                              unique_atom_grad_info(unique_2)%m(2,2,1)
!   do i=1,grad_dim_2
!    print*,dervs_totalsym(index_2+1:index_2+grad_dim_2,index_2+i)
!    enddo

   deallocate(dervs_totsym_temp2, stat=alloc_stat)
            ASSERT(alloc_stat.eq.0)

   end subroutine addto_totsym

    real(kind=r8_kind) function sum_up_dervs(coeff1,coeff2)
      !----------------------------------------------------------------
      !  Purpose: Returns sums of the form Sum(k,l) a_k F_kl b_l
      !           with F_kl being any matrix in the charge fit basis
      !           and a_k and b_l being any type of charge fit coefficients
      !           UB (9/97)
      !  Input:
      !  n_c1, n_gc1 : dimension of the summation indices split into local
      !  n_c2, n_gc2   contributions (1:n_c1 and 1:n_c2) and global contri-
      !                butions (n_c1+1:n_c1+n_gc1 and n_c2+1:n_c2+n_gc2)
      !  help_grad(:,:) : the matrix elements F_kl stored in transposed fashion
      !                   F_kl = help_grad(l,k)
      !  k_start, k_glob_start : pointer on the first element of the local
      !  l_start, l_glob_start   (k_start and l_start) and the global
      !                          (k_glob_start and l_glob_start) contributions
      !                          of the coefficient vectors a_k and b_l
      !  n_independent_1 : storage increment of the local contributions
      !  n_independent_2   to a_k: a_k = a(k_start + (k-1)*n_independent_1) and
      !                    to b_l: b_l = b(l_start + (l-1)*n_independent_2)
      !----------------------------------------------------------------
      implicit none
      !------------ Declaration of formal parameters ------------------
      real(kind=r8_kind), intent(in) :: coeff1(:), coeff2(:)
      ! the coefficient vectors a_k and b_l in linear packed storage mode
      !** End of interface *****************************************
      !------------ Declaration of local variables --------------------
      integer(kind=i4_kind) :: i_a,i_b,k,l
      real   (kind=r8_kind) :: temp ! to hold Sum(l) F_kl b_l
      !------------ Executable code -----------------------------------
      sum_up_dervs = 0.0_r8_kind
      ! < a | F | b >
      k = k_start
      do i_a=1,n_c1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + n_independent_2
         end do
         sum_up_dervs = sum_up_dervs + coeff1(k) * temp
         k = k + n_independent_1
      end do
      ! < ag | F | b >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_start
         temp = 0.0_r8_kind
         do i_b=1,n_c2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + n_independent_2
         end do
         sum_up_dervs = sum_up_dervs + coeff1(k) * temp
         k = k + 1
      end do
      if (n_gc2 == 0) return
      ! < a | F | bg >
      k = k_start
      do i_a=1,n_c1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up_dervs = sum_up_dervs + coeff1(k) * temp
         k = k + n_independent_1
      end do
      ! < ag | F | bg >
      k = k_glob_start
      do i_a=n_c1+1,n_c1+n_gc1
         l = l_glob_start
         temp = 0.0_r8_kind
         do i_b=n_c2+1,n_c2+n_gc2
            temp = temp + help_grad(i_b,i_a) * coeff2(l)
            l = l + 1
         end do
         sum_up_dervs = sum_up_dervs + coeff1(k) * temp
         k = k + 1
      end do
    end function sum_up_dervs


  end subroutine contract_2c_dervs


  subroutine calculate_glob_con_simple( glob_con_simple,glob_con,l,i_if)
    !----------------------------------------------------------------
    !  Purpose: allocation of arrays in and and calculation of data
    !  in glob_con_simple.
    !  The summands with spezified (i_if, l_quadrupel) are copied
    !  from glob_con to glob_con_simple to avoid unnecessary
    !  ifs in inner loops
    !----------------------------------------------------------------
    implicit none
    !------------ Declaration of formal parameters ------------------
    type(glob_con_simple_type), intent(out) :: glob_con_simple(:)
    type(unique_atom_glob_con_type), intent(in) :: glob_con(:)
    integer(kind=i4_kind), intent(in) :: i_if, l
    !** End of interface *****************************************
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: i_gc, i_cont, i_cont_simple, &
         n_cont_simple, status
    !------------ Executable code -----------------------------------
    do i_gc = 1, ubound(glob_con,1)
       n_cont_simple = 0
       do i_cont = 1, glob_con(i_gc)%N_contributing_fcts
          if ( l .eq. glob_con(i_gc)%l(i_cont) .and. &
               i_if .eq. glob_con(i_gc)%index_ind_fct(i_cont) ) then
             n_cont_simple = n_cont_simple + 1
          endif
       enddo
       glob_con_simple(i_gc)%N_contributing_fcts = n_cont_simple
       allocate( glob_con_simple(i_gc)%index_exp(n_cont_simple), &
            glob_con_simple(i_gc)%coefs(n_cont_simple), &
            stat = status )
       if ( status .ne. 0 ) call error_handler( &
            "calculate_glob_con_simple: allocate failed" )
       i_cont_simple = 0
       do i_cont = 1, glob_con(i_gc)%N_contributing_fcts
          if ( l .eq. glob_con(i_gc)%l(i_cont) .and. &
               i_if .eq. glob_con(i_gc)%index_ind_fct(i_cont) ) then
             i_cont_simple = i_cont_simple + 1
             glob_con_simple(i_gc)%index_exp(i_cont_simple) = &
                  glob_con(i_gc)%index_exp(i_cont)
             glob_con_simple(i_gc)%coefs(i_cont_simple) = &
                  glob_con(i_gc)%coefs(i_cont)
          endif
       enddo
    enddo
  end subroutine calculate_glob_con_simple


  subroutine free_glob_con_simple( glob_con_simple )
    !----------------------------------------------------------------
    !  Purpose: allocation of arrays in and and calculation of data
    !  in glob_con_simple.
    !  The summands with spezified (i_if, l_quadrupel) are copied
    !  from glob_con to glob_con_simple to avoid unnecessary
    !  ifs in inner loops
    !----------------------------------------------------------------
    implicit none
    !------------ Declaration of formal parameters ------------------
    type(glob_con_simple_type), intent(inout) :: glob_con_simple(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: i_gc, status
    !------------ Executable code -----------------------------------
    do i_gc = 1, ubound(glob_con_simple,1)
       deallocate( glob_con_simple(i_gc)%index_exp, &
            glob_con_simple(i_gc)%coefs, &
            stat = status )
       if ( status .ne. 0 ) call error_handler( &
            "calculate_glob_con_free: deallocate failed" )
    enddo
  end subroutine free_glob_con_simple

end module fitcontract_2c_grad_module


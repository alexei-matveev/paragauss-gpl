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
module fitcontract_2c_module
  !---------------------------------------------------------------
  !
  !  Purpose: Perform fitcontractions on two-center fit integrals
  !           Re-Normaliation included.
  !
  !  Module called by: integral_calc_2cff
  !
  !  References: Publisher Document: Concepts of Integral Part
  !
  !  Author: FN, TB
  !  Date: 8/96
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: contraction for mixed integrals such as <f_k|g_l>
  !              introduced (see new meaning of parameter FLAG)
  !
  ! Modification (Please copy before editing)
  ! Author:
  ! Date:
  ! Description:
  !
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use unique_atom_module
  use int_data_2cff_module
  use output_module, only: output_int_fitcontract
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public fitcontract_2c


  !================================================================
  ! End of public interface of module
  !================================================================

  type glob_con_simple_type
     ! to hold contribution of one independent function
     !  and one l to global contraction
     integer(kind=i4_kind) :: N_contributing_fcts
     integer(kind=i4_kind), pointer :: index_exp(:)
     real(kind=r8_kind), pointer :: coefs(:)
  end type glob_con_simple_type

  !------------ Declaration of constants and variables ----
  real(kind=r8_kind), parameter :: zero=0.0_r8_kind, one=1.0_r8_kind
  type(unique_atom_basis_type), pointer :: basis1,basis2
  logical :: do_glob_cons1, do_glob_cons2, calc_norm, calc_int
  integer(kind=i4_kind) :: n_gc1,n_gc2,n_c1,n_c2,n_exp1,n_exp2
  type(unique_atom_glob_con_type),pointer :: glob_con1(:), glob_con2(:)
  real(kind=r8_kind),pointer :: loc_loc_int(:,:,:,:), &
       loc_glob_int(:,:,:), glob_loc_int(:,:,:), glob_glob_int(:,:)
  integer(kind=i4_kind),pointer :: glob_loc_contrib(:), loc_glob_contrib(:)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  subroutine fitcontract_2c(fit_int,flag)
    !  Purpose: main routine for performing the contractions
    !           and renormaliations of the two center
    !           fit integrals
    !------------ Modules used ----------------------------------
    use type_module
    use iounitadmin_module, only: write_to_output_units
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(inout) :: fit_int(:,:,:,:)
    ! attention: fit_int is spoiled afterwards
    character(*), intent(in) :: flag ! "ch_ch", "xc_xc", "ch_xc", "xc_ch", or "pre"
    !** End of interface *****************************************
    !------------ Executable code --------------------------------
    ! << set pointer to the first type of basis functions >>
    select case ( flag )
    case ( 'ch_ch', 'ch_xc', 'pre' ) ! charge fit basis functions f_k
       basis1 => ua1_basis_ch
       do_glob_cons1 = n_ch_gc1 .gt. 0
       n_gc1 = n_ch_gc1
       n_c1 = n_ch_c1
       glob_con1 => ua1%glob_con_ch
    case ( 'xc_xc', 'xc_ch' ) ! exchange fit basis functions g_k
       basis1 => ua1_basis_xc
       do_glob_cons1 = n_xc_gc1 .gt. 0
       n_gc1 = n_xc_gc1
       n_c1 = n_xc_c1
       glob_con1 => ua1%glob_con_xc
    end select
    ! << set pointer to the second type of basis functions >>
    select case ( flag )
    case ( 'ch_ch', 'xc_ch', 'pre' ) ! charge fit basis functions f_l
       basis2 => ua2_basis_ch
       do_glob_cons2 = n_ch_gc2 .gt. 0
       n_gc2 = n_ch_gc2
       n_c2 = n_ch_c2
       glob_con2 => ua2%glob_con_ch
    case ( 'xc_xc', 'ch_xc' ) ! exchange fit basis functions g_l
       basis2 => ua2_basis_xc
       do_glob_cons2 = n_xc_gc2 .gt. 0
       n_gc2 = n_xc_gc2
       n_c2 = n_xc_c2
       glob_con2 => ua2%glob_con_xc
    end select
    ! << set pointer to the type of integral >>
    select case ( flag )
    case ( 'ch_ch' ) ! [ f_k | f_l ]
       calc_norm = .true.
       calc_int  = integralpar_2cch_no
       if (calc_int) then
          loc_loc_int => loc_loc_int_ch
          loc_glob_int => loc_glob_int_ch
          glob_loc_int => glob_loc_int_ch
          glob_glob_int => glob_glob_int_ch
          glob_loc_contrib => glob_loc_contrib_ch
          loc_glob_contrib => loc_glob_contrib_ch
       endif

       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: processing integral charge overlap integral")
    case ( 'xc_xc' ) ! < g_k | g_l >
       calc_norm = .true.
       calc_int  = integralpar_2cxc_no
       if (calc_int) then
          loc_loc_int => loc_loc_int_xc
          loc_glob_int => loc_glob_int_xc
          glob_loc_int => glob_loc_int_xc
          glob_glob_int => glob_glob_int_xc
          glob_loc_contrib => glob_loc_contrib_xc
          loc_glob_contrib => loc_glob_contrib_xc
       endif

       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: processing exchange overlap integral")
    case ( 'ch_xc' ) ! < f_k | g_l >
       calc_norm = .false.
       calc_int  = integralpar_2c_mixed
       if (calc_int) then
          loc_loc_int => loc_loc_int_ch_xc
          loc_glob_int => loc_glob_int_ch_xc
          glob_loc_int => glob_loc_int_ch_xc
          glob_glob_int => glob_glob_int_ch_xc
          glob_loc_contrib => glob_loc_contrib_ch_xc
          loc_glob_contrib => loc_glob_contrib_ch_xc
       endif

       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: processing <f_k|g_l> overlap integral")
    case ( 'xc_ch' ) ! < g_k | f_l >
       calc_norm = .false.
       calc_int  = integralpar_2c_mixed
       if (calc_int) then
          loc_loc_int => loc_loc_int_xc_ch
          loc_glob_int => loc_glob_int_xc_ch
          glob_loc_int => glob_loc_int_xc_ch
          glob_glob_int => glob_glob_int_xc_ch
          glob_loc_contrib => glob_loc_contrib_xc_ch
          loc_glob_contrib => loc_glob_contrib_xc_ch
       endif

       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: processing <g_k|f_l> overlap integral")
    case ( 'pre' )

       calc_norm = .false.
       calc_int  = .true.
       if ( calc_int ) then
          loc_loc_int => loc_loc_int_pre
          loc_glob_int => loc_glob_int_pre
          glob_loc_int => glob_loc_int_pre
          glob_glob_int => glob_glob_int_pre
          glob_loc_contrib => glob_loc_contrib_pre
          loc_glob_contrib => loc_glob_contrib_pre
       endif

       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: processing integral charge&
            & estimate overlap integral")

    end select

    n_exp1 = basis1%n_exponents
    n_exp2 = basis2%n_exponents

    if (calc_int) then
       if (output_int_fitcontract) &
            call write_to_output_units &
            ("CONTRACT_2C_FIT: calling contract_2c ")
       call contract_2c(fit_int)
       nullify(loc_loc_int)
       nullify(glob_loc_int)
       nullify(loc_glob_int)
       nullify(glob_glob_int)
       nullify(glob_loc_contrib)
       nullify(loc_glob_contrib)
    endif

  end subroutine fitcontract_2c

  !*************************************************************

  subroutine  contract_2c(fit_int)
    !  Purpose: perform contractions for the two
    !           center fit integrals.
    !
    !  Subroutine called by: fitcontract_2c
    !
    !  Author: FN
    !  Date: 8/96
    !
    !------------ Modules used --------------------------------------
    use type_module ! type specification parameters
    implicit none
    !------------ Declaration of formal parameters ------------------
    real(kind=r8_kind),intent(inout)        :: fit_int(:,:,:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: status, &
         n_uncont1,n_uncont2,i_cont1,i_cont2, &
         n_cont1,n_cont2,first_cont1,first_cont2, &
         i_if1,i_if2,i_exp1,i_exp2,i_c1,i_c2,i_gc1,i_gc2
    real(kind=r8_kind), allocatable :: intermediate(:,:), &
         intermediate_gc2(:,:), intermediate2_gc2(:,:)
    real(kind=r8_kind), pointer :: &
         contract1(:,:),contract2(:,:)
    type(glob_con_simple_type), allocatable :: &
         glob_con_simple1(:), glob_con_simple2(:)
    !------------ Executable code -----------------------------------

    allocate ( intermediate(n_c2,n_exp1), &
         intermediate_gc2(n_gc2,n_exp1), &
         intermediate2_gc2(n_gc2,n_c1), &
         stat = status )
    if ( status .ne. 0 ) call error_handler( &
         "contract_2c: allocate of intermediates failed")

    allocate ( glob_con_simple1(n_gc1), glob_con_simple2(n_gc2), &
         stat = status )
    if ( status .ne. 0 ) call error_handler( &
         "contract_2c: allocate of glob_con_simple failed")

    n_uncont1 = basis1%n_uncontracted_fcts
    n_cont1 = basis1%n_contracted_fcts
    first_cont1 = n_uncont1 + 1
    n_uncont2 = basis2%n_uncontracted_fcts
    n_cont2 = basis2%n_contracted_fcts
    first_cont2 = n_uncont2 + 1

    contract1 => basis1%contractions
    contract2 => basis2%contractions

    loc_loc_int   = zero
    loc_glob_int  = zero
    glob_loc_int  = zero
    glob_glob_int = zero
    glob_loc_contrib = 0
    loc_glob_contrib = 0

    indep1: do i_if1 = 1,ubound(fit_int,3)

       if ( do_glob_cons1 ) then
          call calculate_glob_con_simple(glob_con_simple1, &
               glob_con1, quadrupel%l1, i_if1 )
          glob_loc_contrib = glob_loc_contrib + &
               glob_con_simple1(:)%N_contributing_fcts
       endif
       indep2: do i_if2 = 1,ubound(fit_int,4)

          if ( do_glob_cons2 ) then
             call calculate_glob_con_simple(glob_con_simple2, &
                  glob_con2, quadrupel%l2, i_if2 )
             if ( i_if1 .eq. 1 ) &
                  loc_glob_contrib = loc_glob_contrib + &
                  glob_con_simple2(:)%N_contributing_fcts
          endif


          intermediate = zero
          intermediate_gc2 = zero


          !------ first handle dimension 2 ------
          ! take primitives from fit_int
          ! store results in intermediate(i_c2,i_exp1)
          ! and intermediate_gc2(i_gc2,i_exp1)


          ! local contractions with unrenormalized integrals
          ! inter(n_u2+c2,e1) := Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2)
          i_c2 = first_cont2
          do i_cont2 = 1, n_cont2
             do i_exp2 = 1, n_exp2
                intermediate(i_c2,:) = &
                     intermediate(i_c2,:) +  &
                     contract2(i_exp2,i_cont2) * &
                     fit_int(:,i_exp2,i_if1,i_if2)
             enddo
             i_c2 = i_c2 + 1
          enddo

             ! copy uncontracted to intermediate
             do i_exp2 = 1, n_uncont2
                intermediate(i_exp2,:) = &
                     fit_int(:,i_exp2,i_if1,i_if2)
             enddo

          ! inter(n_u2+c2,e1) := c_norm2(c2) * ...
          !                      Sum(e2) coeff2(e2,c2) prim(e1,e2,i1,i2
          i_c2 = first_cont2
          do i_cont2 = 1, n_cont2
             i_c2 = i_c2 + 1
          enddo

          ! global contractions
          ! inter(g2,e1) := Sum*(e2) coeff2(e2,g2) prim`(e1,e2,i1,i2)
          do i_gc2 = 1, n_gc2
             do i_cont2 = 1, glob_con_simple2(i_gc2)%N_contributing_fcts
                intermediate_gc2(i_gc2,:) = &
                     intermediate_gc2(i_gc2,:) + &
                     glob_con_simple2(i_gc2)%coefs(i_cont2) * &
                     fit_int(:,glob_con_simple2(i_gc2)% &
                     index_exp(i_cont2),i_if1,i_if2)
             enddo
          enddo



          !------ then handle dimension 1 ------
          ! take primitives from intermediate(i_c2,i_exp1)
          ! and intermediate_gc2(i_gc2,i_exp1).
          ! For local contractions store results
          ! in loc_loc_int(i_c2,i_c1,i_if2,i_if1).
          ! For global contractions, add contributions to
          ! glob_loc_int(i_c2,i_if2,i_gc1), loc_glob_int(i_gc2,i_c1,i_if1),
          ! and glob_glob_int(i_gc2,i_gc1)



          ! local contractions with unrenormalized integrals
          ! ll(u&c2,n_u1+c1,i1,i2) := Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
          i_c1 = first_cont1
          do i_cont1 = 1, n_cont1
             do i_exp1 = 1, n_exp1
                loc_loc_int(:,i_c1,i_if2,i_if1) = &
                     loc_loc_int(:,i_c1,i_if2,i_if1) +  &
                     contract1(i_exp1,i_cont1) * &
                     intermediate(:,i_exp1)
             enddo
             i_c1 = i_c1 + 1
          enddo
          ! temp(g2,n_u1+c1) := Sum(e1) coeff1(e1,c1) inter(g2,e1)
          if ( do_glob_cons2 ) then
             intermediate2_gc2 = zero
             i_c1 = first_cont1
             do i_cont1 = 1, n_cont1
                do i_exp1 = 1, n_exp1
                   intermediate2_gc2(:,i_c1) = &
                        intermediate2_gc2(:,i_c1) +  &
                        contract1(i_exp1,i_cont1) * &
                        intermediate_gc2(:,i_exp1)
                enddo
                i_c1 = i_c1 + 1
             enddo
          endif


             ! copy uncontracted to result
             do i_exp1 = 1, n_uncont1
                loc_loc_int(:,i_exp1,i_if2,i_if1) = &
                     intermediate(:,i_exp1)
             enddo

          ! ll(u&c2,n_u1+c1,i1,i2) := c_norm1(c1) * ...
          !                           Sum(e1) coeff1(e1,c1) inter(u&c2,e1)
          i_c1 = first_cont1
          do i_cont1 = 1, n_cont1
             i_c1 = i_c1 + 1
          enddo


          ! do renormalization of global contraction 2
          ! a secound intermediate array is required
          ! because we want to renormalize a summand only
          ! temp(g2,u1) := u_norm1(u1) * inter(g2,u1)
          ! if do_glob_cons1 then
          ! inter`(g2,e1) := u_norm1(e1) * inter(g2,e1)
          if ( do_glob_cons2 ) then
                ! copy uncontracted to result
                do i_exp1 = 1, n_uncont1
                   intermediate2_gc2(:,i_exp1) = &
                        intermediate_gc2(:,i_exp1)
                enddo

             ! temp(g2,n_u1+c1) := c_norm(c1) * ...
             !                     Sum(e1) coeff1(e1,c1) inter(g2,e1)
             i_c1 = first_cont1
             do i_cont1 = 1, n_cont1
                i_c1 = i_c1 + 1
             enddo

             ! lg(g2,u&c1,i1) := temp(g2,u&c1)
             loc_glob_int(:,:,i_if1) = &
                  loc_glob_int(:,:,i_if1) + intermediate2_gc2
          endif


          ! global contractions
          ! gl(u&c2,i2,g1) := Sum*(e1) coeff1(e1,g1) inter`(u&c2,e1)
          do i_gc1 = 1, n_gc1
             do i_cont1 = 1, glob_con_simple1(i_gc1)%N_contributing_fcts
                glob_loc_int(:,i_if2,i_gc1) = &
                     glob_loc_int(:,i_if2,i_gc1) + &
                     glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                     intermediate(:,glob_con_simple1(i_gc1)% &
                     index_exp(i_cont1))
             enddo
          enddo

          ! gg(g2,g1) := Sum*(e1) coeff1(e1,g1) inter`(g2,e1)
          if ( do_glob_cons2 ) then
             do i_gc1 = 1, n_gc1
                do i_cont1 = 1, glob_con_simple1(i_gc1)%N_contributing_fcts
                   glob_glob_int(:,i_gc1) = &
                        glob_glob_int(:,i_gc1) + &
                        glob_con_simple1(i_gc1)%coefs(i_cont1) * &
                        intermediate_gc2(:,glob_con_simple1(i_gc1)% &
                        index_exp(i_cont1))
                enddo
             enddo
          endif


          call free_glob_con_simple(glob_con_simple2)
       enddo indep2

       call free_glob_con_simple(glob_con_simple1)
    enddo indep1

    deallocate ( intermediate, &
         intermediate_gc2, &
         intermediate2_gc2, &
         stat = status )
    if ( status .ne. 0 ) call error_handler( &
         "contract_2c: deallocate of intermediates failed")

    deallocate ( glob_con_simple1, glob_con_simple2, &
         stat = status )
    if ( status .ne. 0 ) call error_handler( &
         "contract_2c: deallocate of glob_con_simple failed")

  end subroutine contract_2c

  !*************************************************************

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

  !*************************************************************

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



end module fitcontract_2c_module

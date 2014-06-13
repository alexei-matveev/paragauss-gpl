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
!===============================================================
! Public interface of module
!===============================================================
module integral_calc_quad_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
#if defined (_COMPAC_FORTRAN) || (_ITANIUM_NSK)
# define FPP_NO_BIG_AUTOMATIC_ARRAYS
#endif
# include "def.h"
  use type_module, only: IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use datatype, only: arrmat2
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  integer(IK), parameter, public :: IDONT_CONTRACT          = 1
  integer(IK), parameter, public :: IDONT_SUM_OVER_PARTNERS = 2


  ! these hold the sections of density matrix and energy-weighted
  ! density matrix of a single quadrupel:
  type(arrmat2), pointer, public :: quad_pmat(:,:) ! (n_irr,n_spin)
  type(arrmat2), pointer, public :: quad_wmat(:)   ! (n_irr)


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: integral_calc_quad_densmat ! (U1,L1,U2,L2,quad_pmat,quad_wmat,use_spin_dens)
  public :: integral_calc_quad_close   ! ()
  public :: contrsym2c                 ! (U1,E1,L1,U2,E2,L2,weight,prim,sym,imode)
  public :: contrsym3c                 ! (U1,E1,L1,U2,E2,L2,weight,prim,sym,imode)
  public :: contr3c                    ! (prim_int_3c_co)
  public :: quad_density_mat           ! (U1,E1,L1,U2,E2,L2,weight,prim_int_dens_mat)
  public :: unsymadp2c                 ! (U1,E1,L1,U2,E2,L2,weight,dmat,cnt)
  public :: symadp2c                   ! (U1,E1,L1,U2,E2,L2,weight,dmat,cnt)
  public :: symadapt_spin
  public :: unsymadapt_spin
  public :: desymmetrize_block
  public :: symmetrize_block
  public :: unsymadaptvec

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  subroutine contrsym2c(U1,E1,L1,U2,E2,L2,weight,prim,sym,imode)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !------------ Declaration of formal parameters -------------
    use int_data_2cob3c_module, only: symadapt_totsym_2c_int_type
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    use int_data_2cob3c_module, only: n_c1, n_c2, n_m1, n_m2
#endif
    implicit none
    integer(IK)                 , intent(in)    :: U1,E1,L1,U2,E2,L2
    real(RK)               , intent(in)    :: weight
    real(RK)               , intent(in)    :: prim(:,:,:,:)
    type(symadapt_totsym_2c_int_type), intent(inout) :: sym(:) ! (n_irreps)
    integer(IK), optional       , intent(in)    :: imode
    !------------ Declaration of local variables ---------------
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    real(RK), allocatable :: cont(:,:,:,:)
    integer(IK)   :: status
#endif
    integer(IK)   :: mode
    !------------ Executable code ------------------------------

    mode = 0
    if( present(imode) ) mode = imode

    if( IAND(mode,IDONT_CONTRACT) /= 0 ) then
      call symadp2c(U1,E1,L1,U2,E2,L2,weight,prim,sym,mode)
    else
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
      allocate(cont(n_c2,n_c1,n_m2,n_m1),stat=status)
      ASSERT(status==0)
      call contract_2center(prim,cont)
      call symadp2c(U1,E1,L1,U2,E2,L2,weight,cont,sym,mode)
      deallocate(cont,stat=status)
      ASSERT(status==0)
#else
      call symadp2c(U1,E1,L1,U2,E2,L2,weight,contr2c(prim),sym,mode)
#endif
    endif
  end subroutine contrsym2c

  !**************************************************************
  subroutine contrsym3c(U1,E1,L1,U2,E2,L2,weight,prim,sym,imode)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !------------ Declaration of formal parameters -------------
    use int_data_2cob3c_module, only: symadapt_totsym_3c_int_type
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    use int_data_2cob3c_module, only: n_c1, n_c2, n_m1, n_m2
#endif
    implicit none
    integer(IK)                 , intent(in)    :: U1,E1,L1,U2,E2,L2
    real(RK)               , intent(in)    :: weight
    real(RK)               , intent(in)    :: prim(:,:,:,:,:)
    type(symadapt_totsym_3c_int_type), intent(inout) :: sym(:) ! (n_irreps)
    integer(IK), optional       , intent(in)    :: imode
    !------------ Declaration of local variables ---------------
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    real(RK), allocatable :: cont(:,:,:,:,:)
    integer(IK)   :: status
#endif
    integer(IK)   :: mode
    !------------ Executable code ------------------------------

    mode = 0
    if( present(imode) ) mode = imode

    if( IAND(mode,IDONT_CONTRACT) /= 0 ) then
      call symadp3c(U1, E1, L1, U2, E2, L2, weight, prim, sym, mode)
    else
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
      allocate(cont(n_c2, n_c1, size(prim, 3), n_m2, n_m1), stat=status)
      ASSERT(status==0)
      call contract_3center(prim, cont)
      call symadp3c(U1, E1, L1, U2, E2, L2, weight, cont, sym, mode)
      deallocate(cont, stat=status)
      ASSERT(status==0)
#else
      call symadp3c(U1, E1, L1, U2, E2, L2, weight, contr3c(prim), sym, mode)
#endif
    endif
  end subroutine contrsym3c

  subroutine contract_2center(prim,cont)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !------------ Declaration of formal parameters -------------
    use int_data_2cob3c_module, only: n_uncontr1, n_uncontr2       &
                                    , n_m1, n_m2, n_exp1, n_exp2   &
                                    , n_c1, n_c2                   &
                                    , contractions1, contractions2
    implicit none
    real(RK), intent(in) , dimension(:,:,:,:) :: prim
    real(RK), intent(out), dimension(:,:,:,:) :: cont
    !------------ Declaration of local variables ---------------
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    real(RK), allocatable :: int(:)
    integer(IK) :: status
#else
    real(RK)    :: int(n_c2)
#endif
    real(RK)    :: coef
    integer(IK) :: i_m1, i_m2, i_exp1, i_exp2, i_c1, &
                             i_c2, i_co1, i_co2
    !------------ Executable code ------------------------------
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    allocate( int(n_c2), stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_2cob3c: contract_2center: allocate of int failed")
#endif
    cont(:,n_uncontr1+1:,:,:) = 0.0_RK
    m1: do i_m1 = 1, n_m1
       m2: do i_m2 = 1, n_m2
          exp1_uc: do i_exp1 = 1, n_uncontr1
             do i_c2 =1, n_uncontr2
                int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1)
             enddo
             i_co2 = 1
             do i_c2 = n_uncontr2 + 1, n_c2
                int(i_c2) = 0.0_RK
                do i_exp2 = 1, n_exp2
                   int(i_c2) = int(i_c2) + &
                        contractions2(i_exp2,i_co2) * &
                        prim(i_exp2,i_exp1,i_m2,i_m1)
                enddo
                i_co2 = i_co2 + 1
             enddo
             i_co1 = 1
             do i_c1= n_uncontr1 + 1, n_c1
                coef = contractions1(i_exp1,i_co1)
                do i_c2 = 1, n_c2
                   cont(i_c2,i_c1,i_m2,i_m1) = &
                        cont(i_c2,i_c1,i_m2,i_m1) + &
                        coef * int(i_c2)
                enddo
                i_co1 = i_co1 + 1
             enddo
             do i_c2 = 1, n_c2
                cont(i_c2,i_exp1,i_m2,i_m1) = int(i_c2)
             enddo
          enddo exp1_uc
          exp1_c: do i_exp1 = n_uncontr1 + 1, n_exp1
             do i_c2 =1, n_uncontr2
                int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1)
             enddo
             i_co2 = 1
             do i_c2 = n_uncontr2 + 1, n_c2
                int(i_c2) = 0.0_RK
                do i_exp2 = 1, n_exp2
                   int(i_c2) = int(i_c2) + &
                        contractions2(i_exp2,i_co2) * &
                        prim(i_exp2,i_exp1,i_m2,i_m1)
                enddo
                i_co2 = i_co2 + 1
             enddo
             i_co1 = 1
             do i_c1= n_uncontr1 + 1, n_c1
                coef = contractions1(i_exp1,i_co1)
                do i_c2 = 1, n_c2
                   cont(i_c2,i_c1,i_m2,i_m1) = &
                        cont(i_c2,i_c1,i_m2,i_m1) + &
                        coef * int(i_c2)
                enddo
                i_co1 = i_co1 + 1
             enddo
          enddo exp1_c
       enddo m2
    enddo m1
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    deallocate( int, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_2cob3c: contract_2center: deallocate of int failed")
#endif
  end subroutine contract_2center

  subroutine uncontract_2center(cont,prim)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !------------ Declaration of formal parameters -------------
    use int_data_2cob3c_module, only: n_uncontr1, n_uncontr2       &
                                    , n_m1, n_m2, n_exp1, n_exp2   &
                                    , n_c1, n_c2                   &
                                    , contractions1, contractions2
    implicit none
    real(RK), intent(in)  :: cont(:,:,:,:)
    real(RK), intent(out) :: prim(:,:,:,:)
    !------------ Declaration of local variables ---------------
    real(RK)    :: int(n_c2)
    real(RK)    :: coef
    integer(IK) :: i_m1, i_m2, i_exp1, i_exp2, i_c1, &
                             i_c2, i_co1, i_co2 !, status
    !------------ Executable code ------------------------------

    m1: do i_m1 = 1, n_m1
       m2: do i_m2 = 1, n_m2

          exp1_uc: do i_exp1 = 1, n_uncontr1

             prim(:,i_exp1,i_m2,i_m1) = 0.0_rk

             int=0.0_rk

             do i_c2 = 1, n_c2
                int(i_c2)=cont(i_c2,i_exp1,i_m2,i_m1)
             enddo

             i_co1 = 1
             do i_c1= n_uncontr1 + 1, n_c1
                coef = contractions1(i_exp1,i_co1)
                do i_c2 = 1, n_c2
                   int(i_c2)=int(i_c2)+cont(i_c2,i_c1,i_m2,i_m1)*coef
                enddo
                i_co1 = i_co1 + 1
             enddo

             i_co2 = 1
             do i_c2 = n_uncontr2 + 1, n_c2
                do i_exp2 = 1, n_exp2
                prim(i_exp2,i_exp1,i_m2,i_m1)=prim(i_exp2,i_exp1,i_m2,i_m1)+&
                     contractions2(i_exp2,i_co2) * int(i_c2)
                enddo
                i_co2 = i_co2 + 1
             enddo

             do i_c2 =1, n_uncontr2
                prim(i_c2,i_exp1,i_m2,i_m1)=prim(i_c2,i_exp1,i_m2,i_m1)+int(i_c2)
             enddo

          enddo exp1_uc

          exp1_c: do i_exp1 = n_uncontr1 + 1, n_exp1

             prim(:,i_exp1,i_m2,i_m1) = 0.0_rk

             int=0.0_rk

             i_co1 = 1
             do i_c1= n_uncontr1 + 1, n_c1
                coef = contractions1(i_exp1,i_co1)
                do i_c2 = 1, n_c2
                int(i_c2)=int(i_c2)+cont(i_c2,i_c1,i_m2,i_m1)*coef
                enddo
                i_co1 = i_co1 + 1
             enddo

             i_co2 = 1
             do i_c2 = n_uncontr2 + 1, n_c2
                do i_exp2 = 1, n_exp2
                 prim(i_exp2,i_exp1,i_m2,i_m1)=prim(i_exp2,i_exp1,i_m2,i_m1) &
                +contractions2(i_exp2,i_co2) * int(i_c2)
                enddo
                i_co2 = i_co2 + 1
             enddo

             do i_c2 =1, n_uncontr2
                prim(i_c2,i_exp1,i_m2,i_m1)=prim(i_c2,i_exp1,i_m2,i_m1)+int(i_c2)
             enddo

          enddo exp1_c

       enddo m2
    enddo m1

  end subroutine uncontract_2center
  !**************************************************************


  !**************************************************************
  subroutine contract_3center(prim,cont)
    !  Purpose: does contraction on three center integrals
    !  the result is in cont that is allocated by this subroutine
    use int_data_2cob3c_module, only: n_uncontr1, n_uncontr2       &
                                    , n_m1, n_m2, n_exp1, n_exp2   &
                                    , n_c1, n_c2                   &
                                    , contractions1, contractions2
    implicit none
    !------------ Declaration of formal parameters -------------
    real(RK), intent(in) , dimension(:,:,:,:,:) :: prim
    real(RK), intent(out), dimension(:,:,:,:,:) :: cont
    !------------ Declaration of local variables ---------------
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    real(RK), allocatable :: int(:,:,:,:), prim_help(:,:,:,:,:), &
         cont_help(:,:,:,:,:)
    integer(IK) :: status
#else
    real(RK)    :: int(      size(prim,3),n_m2,n_m1,n_c2)
    real(RK)    :: cont_help(size(prim,3),n_m2,n_m1,n_c2,n_c1)
    real(RK)    :: prim_help(size(prim,3),n_m2,n_m1,n_exp2,n_exp1)
#endif
    real(RK)    :: coef
    integer(IK) :: i_exp1, i_exp2, i_c1, &
         i_c2, i_co1, i_co2, n_f
    !------------ Executable code ------------------------------
    FPP_TIMER_START(cntr3)
    n_f = ubound(prim,3)
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    allocate( int(n_f,n_m2,n_m1,n_c2), cont_help(n_f,n_m2,n_m1,n_c2,n_c1),&
         prim_help(n_f,n_m2,n_m1,n_exp2,n_exp1),stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_2cob3c: contract_3center: allocate of int failed")
#endif
    cont_help(:,:,:,:,n_uncontr1+1:) = 0.0_RK
    do i_exp1 = 1,n_exp1
       do i_exp2 = 1,n_exp2
          prim_help(:,:,:,i_exp2,i_exp1)=prim(i_exp2,i_exp1,:,:,:)
       end do
    end do
    exp1_uc: do i_exp1 = 1, n_uncontr1
       int(:,:,:,1:n_uncontr2) = prim_help(:,:,:,1:n_uncontr2,i_exp1)
       i_co2 = 1
       do i_c2 = n_uncontr2 + 1, n_c2
          int(:,:,:,i_c2) = 0.0_RK
          do i_exp2 = 1, n_exp2
             int(:,:,:,i_c2) = int(:,:,:,i_c2) + &
                  contractions2(i_exp2,i_co2) * &
                  prim_help(:,:,:,i_exp2,i_exp1)
          enddo
          i_co2 = i_co2 + 1
       enddo
       i_co1 = 1
       do i_c1= n_uncontr1 + 1, n_c1
          coef = contractions1(i_exp1,i_co1)
          do i_c2 = 1, n_c2
             cont_help(:,:,:,i_c2,i_c1) = &
                  cont_help(:,:,:,i_c2,i_c1) + &
                  coef * int(:,:,:,i_c2)
          enddo
          i_co1 = i_co1 + 1
       enddo
       cont_help(:,:,:,:,i_exp1) = int(:,:,:,:)
    enddo exp1_uc
    exp1_c: do i_exp1 = n_uncontr1 + 1, n_exp1
       int(:,:,:,1:n_uncontr2) = prim_help(:,:,:,1:n_uncontr2,i_exp1)
       i_co2 = 1
       do i_c2 = n_uncontr2 + 1, n_c2
          int(:,:,:,i_c2) = 0.0_RK
          do i_exp2 = 1, n_exp2
             int(:,:,:,i_c2) = int(:,:,:,i_c2) + &
                  contractions2(i_exp2,i_co2) * &
                  prim_help(:,:,:,i_exp2,i_exp1)
          enddo
          i_co2 = i_co2 + 1
       enddo
       i_co1 = 1
       do i_c1= n_uncontr1 + 1, n_c1
          coef = contractions1(i_exp1,i_co1)
          do i_c2 = 1, n_c2
             cont_help(:,:,:,i_c2,i_c1) = &
                  cont_help(:,:,:,i_c2,i_c1) + &
                  coef * int(:,:,:,i_c2)
          enddo
          i_co1 = i_co1 + 1
       enddo
    enddo exp1_c
    do i_c1 = 1,n_c1
       do i_c2 = 1,n_c2
          cont(i_c2,i_c1,:,:,:)=cont_help(:,:,:,i_c2,i_c1)
       end do
    end do
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
    deallocate( int, prim_help, cont_help, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_2cob3c: contract_3center: deallocate of int failed")
#endif
    FPP_TIMER_STOP(cntr3)
  end subroutine contract_3center
  !**************************************************************

  !**************************************************************
  function contr2c(prim) result(cont)
    ! wrapper for contract_2center
    use int_data_2cob3c_module, only: n_c1, n_c2, n_m1, n_m2
    implicit none
    real(RK), intent(in)  :: prim(:,:,:,:)
    real(RK)              :: cont(n_c2,n_c1,n_m2,n_m1)
    ! *** end of interface ***

    call contract_2center(prim,cont)
  end function contr2c

  !**************************************************************
  function contr3c(prim) result(cont)
    ! wrapper for contract_3center
    use int_data_2cob3c_module, only: n_c1, n_c2, n_m1, n_m2
    implicit none
    real(RK), intent(in)  :: prim(:,:,:,:,:)
    real(RK)              :: cont(n_c2,n_c1,size(prim,3),n_m2,n_m1)
    ! *** end of interface ***

    call contract_3center(prim,cont)
  end function contr3c

  !**************************************************************
  subroutine symadp2c(U1,E1,L1,U2,E2,L2,weight,cnt,sym,imode)
    ! add contribution of one not symmetry equivalent pair
    ! of equal atom to symmetry adapted integrals
    ! the case of total symmetric integrals
    use int_data_2cob3c_module, only: symadapt_totsym_2c_int_type
    use unique_atom_module, only: unique_atom_partner_type &
                                , unique_atom_sa_int_type  &
                                , unique_atoms
    use symmetry_data_module, only: symmetry_data_n_irreps &
                                  , symmetry_data_n_partners
    implicit none
    integer(IK)                 , intent(in)    :: U1,E1,L1,U2,E2,L2
    real(RK)                    , intent(in)    :: weight
    real(RK)                    , intent(in)    :: cnt(:,:,:,:) ! (nc2,nc1,2*L2+1,2*L1+1)
    type(symadapt_totsym_2c_int_type), intent(inout) :: sym(:)   ! (n_irreps)
    integer(IK), optional       , intent(in)    :: imode
    !------------ Declaration of local variables ---------------
    real(RK), pointer, dimension(:,:,:,:) :: sai
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(IK) :: i_ir, i_pa, i_if1, i_if2, &
         i_cf1, i_cf2, m1, m2
    real(RK) :: coef1, coef, weight2
    integer(IK) :: n_pa
    integer(IK) :: mode
    !------------ Executable code ------------------------------

    mode = 0
    if( present(imode) ) mode = imode


    irreps: do i_ir = 1, symmetry_data_n_irreps()
       sap1 => unique_atoms(U1)%symadapt_partner(i_ir,L1)
       sap2 => unique_atoms(U2)%symadapt_partner(i_ir,L2)

       sai => sym(i_ir)%int

       if( IAND(mode,IDONT_SUM_OVER_PARTNERS)/=0 )then
         ! all totally symmetric integrals dont depend
         ! on partner index, compute only i_pa=1:
         n_pa    = 1
         weight2 = weight
       else
         ! only symmetry inequivalent distances are used,
         ! weighted sum over partners is required:
!        weight = symequivatoms%weight(i_symequiv)
         n_pa    = symmetry_data_n_partners(i_ir)
         weight2 = weight / n_pa
       endif

       partners: do i_pa = 1, n_pa
          ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
             sat1 => sap1%sa_int(E1,i_if1,i_pa)
             ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                sat2 => sap2%sa_int(E2,i_if2,i_pa)
                m_sum_1: do i_cf1 = 1, sat1%N_fcts
                   m1 = sat1%m(i_cf1)
                   coef1 = sat1%c(i_cf1) * weight2
                   m_sum_2: do i_cf2 = 1, sat2%N_fcts
                      m2 = sat2%m(i_cf2)
                      coef = sat2%c(i_cf2) * coef1

                         sai(:,:,i_if2,i_if1) =  sai(:,:,i_if2,i_if1) &
                                        + coef * cnt(:,:,m2,m1)

                   enddo m_sum_2
                enddo m_sum_1
             enddo ind_fct_2
          enddo ind_fct_1
       enddo partners
    enddo irreps
  end subroutine symadp2c

  subroutine unsymadp2c(U1, E1, L1, U2, E2, L2, weight, dmat, cnt)
    ! add contribution of one not symmetry equivalent pair
    ! of equal atom to symmetry adapted integrals
    ! the case of total symmetric integrals
    use unique_atom_module, only: unique_atom_partner_type &
                                , unique_atom_sa_int_type  &
                                , unique_atoms
    use symmetry_data_module, only: symmetry_data_n_partners
    implicit none
    integer(IK)                 , intent(in)    :: U1,E1,L1,U2,E2,L2
    real(RK)                    , intent(in)    :: weight
    type(arrmat2)               , intent(in)    :: dmat(:)      ! (n_irr)%(:,:)
    real(RK)                    , intent(out)   :: cnt(:,:,:,:) ! (n_c2,n_c1,2*L2+1,2*L1+1)
    !------------ Declaration of local variables ---------------
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(IK) :: i_ir,i_pa, i_if1, i_if2, &
                   i_cf1, i_cf2, m1, m2
    real(RK) :: coef1, coef, weight2
    integer(IK) :: n_pa
    integer(IK) :: i_exp1,i_exp2,mu,nu
    !------------ Executable code ------------------------------

     cnt(:,:,:,:) = 0.0

     irreps: do i_ir = 1, size(dmat) ! n_irr
         sap1 => unique_atoms(U1)%symadapt_partner(i_ir, L1)
         sap2 => unique_atoms(U2)%symadapt_partner(i_ir, L2)

         n_pa    = symmetry_data_n_partners(i_ir)
         weight2 = weight / n_pa

         partners: do i_pa = 1, n_pa
             !
             ! Loop over nu == (i_exp1, i_if1):
             !
             nu = 0
             ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
             sat1 => sap1%sa_int(E1, i_if1, i_pa)
             do i_exp1 = 1, size(cnt, 2)
                 nu = nu + 1

                 !
                 ! Loop over mu == (i_exp2, i_if2):
                 !
                 mu = 0
                 ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                 sat2 => sap2%sa_int(E2, i_if2, i_pa)
                 do i_exp2 = 1,size(cnt, 1)
                     mu = mu + 1

                     !
                     ! Each symmetry independent function |i_if| is a linear
                     ! combination of |N_fcts| harmonics with different magnetic
                     ! moments |m| and weight |coef|:
                     !
                     m_sum_1: do i_cf1 = 1, sat1%N_fcts
                         m1 = sat1%m(i_cf1)
                         coef1 = sat1%c(i_cf1) * weight2
                         m_sum_2: do i_cf2 = 1, sat2%N_fcts
                             m2 = sat2%m(i_cf2)
                             coef = sat2%c(i_cf2) * coef1

                             cnt(i_exp2, i_exp1, m2, m1) = cnt(i_exp2, i_exp1, m2, m1) + coef * dmat(i_ir)%m(mu, nu)

                         enddo m_sum_2
                     enddo m_sum_1
                 enddo
                 enddo ind_fct_2
             enddo
             enddo ind_fct_1
         enddo partners
     enddo irreps
  end subroutine unsymadp2c

  !**************************************************************
  subroutine desymmetrize_block( d1, d2, b, i, u1, e1, c1, l1, u2, e2, c2, l2  &
                               , block, blk_sa, uas                            )
    use type_module,               only : IK => i4_kind                        &
                                        , RK => r8_kind
    use unique_atom_module,        only : unique_atom_type
    use symmetry_data_module,      only : symmetry_data_n_partners
    use dimensions,                only : dimoff                               &
                                        , CNF => CONTRACTED
    implicit none
    integer(IK),                        intent(in)    :: d1, d2, b, i
    integer(IK),                        intent(in)    :: u2, e2, c2, l2
    integer(IK),                        intent(in)    :: u1, e1, c1, l1
    real(RK),                           intent(inout) :: block(d1, d2)
    real(RK),                           intent(in)    :: blk_sa(b, b)
    real(RK)                                          :: wgt, Z1, Z2
    integer(IK)                                       :: ipr, nprt
    integer(IK)                                       :: nif1, nif2, nf1, nf2
    integer(IK)                                       :: x1, x2, f1, f2, m1, m2
    integer(IK)                                       :: r1, r2, a1, a2, w1, w2
    integer(IK)                                       :: s1, s2
    integer(IK)                                       :: mo, no, m, n
    type(unique_atom_type),             intent(in)    :: uas(:)

    associate( sap2 => uas(u2)%symadapt_partner(i, l2)                         &
             , sap1 => uas(u1)%symadapt_partner(i, l1)                         )
      nif2 = sap2%N_independent_fcts
      nif1 = sap1%N_independent_fcts
      IF ( nif2 < 1 .or. nif1 < 1 ) RETURN
      x2   = dimoff(i, u2, l2, CNF)
      x1   = dimoff(i, u1, l1, CNF)
      !
      w2   = 2 * l2 + 1
      w1   = 2 * l1 + 1
      !
      nprt = symmetry_data_n_partners(i)
      wgt  = 1.0_RK / real( nprt, RK )
      !
      DO ipr = 1, nprt
        no = x2
        DO f2 = 1, nif2
          nf2  = sap2%sa_int(e2, f2, ipr)%N_fcts
          s2   = 0
          DO r2 = 1, c2
            n = no + r2
            DO a2 = 1, nf2
              m2 = sap2%sa_int(e2, f2, ipr)%m(a2) + s2
              Z2 = sap2%sa_int(e2, f2, ipr)%c(a2) * wgt
              mo = x1
              DO f1 = 1, nif1
                nf1  = sap1%sa_int(e1, f1, ipr)%N_fcts
                s1   = 0
                DO r1 = 1, c1
                  m = mo + r1
                  DO a1 = 1, nf1
                    m1 = sap1%sa_int(e1, f1, ipr)%m(a1) + s1
                    Z1 = sap1%sa_int(e1, f1, ipr)%c(a1) * Z2
                    block(m1, m2) = block(m1, m2) + Z1 * blk_sa(m, n)
                  ENDDO
                  s1 = s1 + w1
                ENDDO
                mo = mo + c1
              ENDDO
            ENDDO
            s2 = s2 + w2
          ENDDO
          no = no + c2
        ENDDO
      ENDDO
    END associate

  end subroutine desymmetrize_block

  !**************************************************************
  subroutine symmetrize_block( d1, d2, b, i, u1, e1, c1, l1, u2, e2, c2, l2    &
                               , block, blk_sa, uas                            )
    use type_module,               only : IK => i4_kind                        &
                                        , RK => r8_kind
    use unique_atom_module,        only : unique_atom_type
    use symmetry_data_module,      only : symmetry_data_n_partners
    use dimensions,                only : dimoff                               &
                                        , CNF => CONTRACTED
    implicit none
    integer(IK),                        intent(in)    :: d1, d2, b, i
    integer(IK),                        intent(in)    :: u2, e2, c2, l2
    integer(IK),                        intent(in)    :: u1, e1, c1, l1
    real(RK),                           intent(in)    :: block(d1, d2)
    real(RK),                           intent(inout) :: blk_sa(b, b)
    real(RK)                                          :: wgt, Z1, Z2
    integer(IK)                                       :: ipr, nprt
    integer(IK)                                       :: nif1, nif2, nf1, nf2
    integer(IK)                                       :: x1, x2, f1, f2, m1, m2
    integer(IK)                                       :: r1, r2, a1, a2, w1, w2
    integer(IK)                                       :: s1, s2
    integer(IK)                                       :: mo, no, m, n
    type(unique_atom_type),             intent(in)    :: uas(:)

    associate( sap2 => uas(u2)%symadapt_partner(i, l2)                         &
             , sap1 => uas(u1)%symadapt_partner(i, l1)                         )
      nif2 = sap2%N_independent_fcts
      nif1 = sap1%N_independent_fcts
      IF ( nif2 < 1 .or. nif1 < 1 ) RETURN
      x2   = dimoff(i, u2, l2, CNF)
      x1   = dimoff(i, u1, l1, CNF)
      !
      w2   = 2 * l2 + 1
      w1   = 2 * l1 + 1
      !
      nprt = symmetry_data_n_partners(i)
      wgt  = 1.0_RK / real( nprt, RK )
      !
      DO ipr = 1, nprt
        no = x2
        DO f2 = 1, nif2
          nf2  = sap2%sa_int(e2, f2, ipr)%N_fcts
          s2   = 0
          DO r2 = 1, c2
            n = no + r2
            DO a2 = 1, nf2
              m2 = sap2%sa_int(e2, f2, ipr)%m(a2) + s2
              Z2 = sap2%sa_int(e2, f2, ipr)%c(a2) * wgt
              mo = x1
              DO f1 = 1, nif1
                nf1  = sap1%sa_int(e1, f1, ipr)%N_fcts
                s1   = 0
                DO r1 = 1, c1
                  m = mo + r1
                  DO a1 = 1, nf1
                    m1 = sap1%sa_int(e1, f1, ipr)%m(a1) + s1
                    Z1 = sap1%sa_int(e1, f1, ipr)%c(a1) * Z2
                    blk_sa(m, n) = blk_sa(m, n) + Z1 * block(m1, m2)
                  ENDDO
                  s1 = s1 + w1
                ENDDO
                mo = mo + c1
              ENDDO
            ENDDO
            s2 = s2 + w2
          ENDDO
          no = no + c2
        ENDDO
      ENDDO
    END associate

  end subroutine symmetrize_block

  subroutine unsymadaptvec( uas, irr, u1, l1, e1, c1, vec_sy, a, b,     vec_ao )
    use type_module,               only : IK => i4_kind                        &
                                        , RK => r8_kind
    use datatype,                  only : arrmat3
    use unique_atom_module,        only : unique_atom_type                     &
                                        , unique_atom_partner_type             &
                                        , unique_atom_sa_int_type
    use symmetry_data_module,      only : symmetry_data_n_partners
    use dimensions,                only : dimoff                               &
                                        , CNF => CONTRACTED
    type(unique_atom_type),       target, intent(in)  :: uas(:)
    real(RK),                             intent(in)  :: vec_sy(:)
    integer(IK),                          intent(in)  :: irr, u1, e1, l1, c1
    integer(IK),                          intent(in)  :: a, b
    real(RK),                             intent(out) :: vec_ao(a,b)
    integer(IK)                                       :: ipr
    integer(IK)                                       :: m, mo
    integer(IK)                                       :: a1, f1, nprt, nif1, nf1
    integer(IK)                                       :: r1, s1, w1, x1, m1
    real(RK)                                          :: Z1

    vec_ao = 0.0_RK
    associate( sap1 => uas(u1)%symadapt_partner(irr, l1)                       )
      nif1 = sap1%N_independent_fcts
      IF ( nif1 < 1 ) RETURN
      x1   = dimoff(irr, u1, l1, CNF)
      !
      w1   = 2 * l1 + 1
      !
      nprt = symmetry_data_n_partners(irr)
      !
      DO ipr = 1, nprt
        mo = x1
        DO f1 = 1, nif1
          nf1  = sap1%sa_int(e1, f1, ipr)%N_fcts
          s1   = 0
          DO r1 = 1, c1
            m = mo + r1
            DO a1 = 1, nf1
              m1 = sap1%sa_int(e1, f1, ipr)%m(a1) + s1
              Z1 = sap1%sa_int(e1, f1, ipr)%c(a1)
              vec_ao(m1, ipr) = vec_ao(m1, ipr) + Z1 * vec_sy(m)
            ENDDO
            s1 = s1 + w1
          ENDDO
          mo = mo + c1
        ENDDO
      ENDDO
    END associate
  end subroutine unsymadaptvec

  !**************************************************************
  subroutine symadapt_spin( uas, mat_ao,                                mat_sy )
    use type_module,               only : IK => i4_kind                        &
                                        , RK => r8_kind
    use datatype,                  only : arrmat3
    use unique_atom_module,        only : unique_atom_type                     &
                                        , unique_atom_partner_type             &
                                        , unique_atom_sa_int_type
    use symmetry_data_module,      only : symmetry_data_n_partners             &
                                        , symmetry_data_n_irr                  &
                                        , symmetry_data_n_spin
    use dimensions,                only : dimoff                               &
                                        , dimens                               &
                                        , CNF => CONTRACTED
    type(unique_atom_type),       target, intent(in)  :: uas(:)
    type(arrmat3),                        intent(in)  :: mat_ao(:,:)
    type(arrmat3),           allocatable, intent(out) :: mat_sy(:)
    integer(IK)                                       :: nirr, nprt, spin
    integer(IK)                                       :: irr, ipr
    integer(IK)                                       :: m, n, s
    integer(IK)                                       :: u1, u2, e1, e2, l1, l2
    integer(IK)                                       :: a1, a2, f1, f2, d1
    integer(IK)                                       :: r1, r2, s1, s2, x1, x2
    integer(IK)                                       :: m1, m2, c1, c2
    real(RK)                                          :: wgt, fct
    type(unique_atom_partner_type),         pointer   :: sap1, sap2
    type(unique_atom_sa_int_type),          pointer   :: sat1, sat2

    nirr = symmetry_data_n_irr()
    spin = symmetry_data_n_spin()

    ALLOCATE( mat_sy(nirr) )
    DO irr = 1, nirr
      d1 = dimens(irr, CNF)
      ALLOCATE( mat_sy(irr)%m(d1, d1, spin) )
      mat_sy(irr)%m = 0.0_RK
    ENDDO

    s2 = 0
    DO u2 = 1, size(uas)
      DO l2 = 0, uas(u2)%lmax_ob
        DO e2 = 1, uas(u2)%n_equal_atoms
          c2 = uas(u2)%l_ob(l2)%N_contracted_fcts                              &
             + uas(u2)%l_ob(l2)%N_uncontracted_fcts
          s2 = s2 + 1
          s1 = 0
          DO u1 = 1, size(uas)
            DO l1 = 0, uas(u1)%lmax_ob
              DO e1 = 1, uas(u1)%n_equal_atoms
                c1 = uas(u1)%l_ob(l1)%N_contracted_fcts                        &
                   + uas(u1)%l_ob(l1)%N_uncontracted_fcts
                s1 = s1 + 1
                !
                DO s = 1, spin
                  DO irr = 1, nirr
                    sap2 => uas(u2)%symadapt_partner(irr, l2)
                    sap1 => uas(u1)%symadapt_partner(irr, l1)
                    x1 = dimoff(irr, u1, l1, CNF)
                    x2 = dimoff(irr, u2, l2, CNF)
                    !
                    nprt = symmetry_data_n_partners(irr)
                    wgt  = 1.0_RK / nprt
                    DO ipr = 1, nprt
                      n = x2
                      DO f2 = 1, sap2%N_independent_fcts
                        sat2 => sap2%sa_int(e2, f2, ipr)
                        DO r2 = 1, c2
                          n = n + 1
                          m = x1
                          DO f1 = 1, sap1%N_independent_fcts
                            sat1 => sap1%sa_int(e1, f1, ipr)
                            DO r1 = 1, c1
                              m = m + 1
                              DO a2 = 1, sat2%N_fcts
                                m2 = sat2%m(a2)
                                fct = sat2%c(a2) * wgt
                                DO a1 = 1, sat1%N_fcts
                                  m1 = sat1%m(a1)
                                  ! FIXME: optimize indexation, use long loops
                                  mat_sy(irr)%m(m, n, s) =                     &
                                      mat_sy(irr)%m(m, n, s)                   &
                                    + fct * sat1%c(a1)                         &
                                          * mat_ao(s1, s2)%m(r1+(m1-1)*c1&
                                                           , r2+(m2-1)*c2&
                                                           , s)
!                                         * mat_ao(s1, s2)%m(m1+(2*l1+1)*(r1-1)&
!                                                          , m2+(2*l2+1)*(r2-1)&
!                                                          , s)
                                ENDDO
                              ENDDO
                            ENDDO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  end subroutine symadapt_spin

  !**************************************************************
  subroutine unsymadapt_spin( uas, mat_sy,                              mat_ao )
    use type_module,               only : IK => i4_kind                        &
                                        , RK => r8_kind
    use datatype,                  only : arrmat3
    use unique_atom_module,        only : unique_atom_type                     &
                                        , unique_atom_partner_type             &
                                        , unique_atom_sa_int_type
    use symmetry_data_module,      only : symmetry_data_n_partners             &
                                        , symmetry_data_n_irr                  &
                                        , symmetry_data_n_spin
    use dimensions,                only : dimoff                               &
                                        , CNF => CONTRACTED
    type(unique_atom_type),       target, intent(in)  :: uas(:)
    type(arrmat3),                        intent(in)  :: mat_sy(:)
    type(arrmat3),           allocatable, intent(out) :: mat_ao(:,:)
    integer(IK)                                       :: nirr, nshl, nprt, spin
    integer(IK)                                       :: irr, ipr
    integer(IK)                                       :: m, n, s
    integer(IK)                                       :: u1, u2, e1, e2, l1, l2
    integer(IK)                                       :: a1, a2, f1, f2
    integer(IK)                                       :: r1, r2, s1, s2, x1, x2
    integer(IK)                                       :: m1, m2, c1, c2
    real(RK)                                          :: wgt, fct
    type(unique_atom_partner_type),         pointer   :: sap1, sap2
    type(unique_atom_sa_int_type),          pointer   :: sat1, sat2

    nirr = symmetry_data_n_irr()
    nshl = sum( uas(:)%n_equal_atoms * (uas(:)%lmax_ob + 1) )
    spin = symmetry_data_n_spin()

    ALLOCATE( mat_ao(nshl, nshl) )

    s2 = 0
    DO u2 = 1, size(uas)
      DO l2 = 0, uas(u2)%lmax_ob
        DO e2 = 1, uas(u2)%n_equal_atoms
          c2 = uas(u2)%l_ob(l2)%N_contracted_fcts                              &
             + uas(u2)%l_ob(l2)%N_uncontracted_fcts
          s2 = s2 + 1
          s1 = 0
          DO u1 = 1, size(uas)
            DO l1 = 0, uas(u1)%lmax_ob
              DO e1 = 1, uas(u1)%n_equal_atoms
                c1 = uas(u1)%l_ob(l1)%N_contracted_fcts                        &
                   + uas(u1)%l_ob(l1)%N_uncontracted_fcts
                s1 = s1 + 1
                ALLOCATE( mat_ao(s1, s2)%m((2*l1+1)*c1, (2*l2+1)*c2, spin) )
                mat_ao(s1, s2)%m = 0.0_RK
                !
                DO irr = 1, nirr
                  sap2 => uas(u2)%symadapt_partner(irr, l2)
                  sap1 => uas(u1)%symadapt_partner(irr, l1)
                  x1 = dimoff(irr, u1, l1, CNF)
                  x2 = dimoff(irr, u2, l2, CNF)
                  !
                  nprt = symmetry_data_n_partners(irr)
                  wgt  = 1.0_RK / nprt
                  !
                  DO s = 1, spin
                    DO ipr = 1, nprt
                      n = x2
                      DO f2 = 1, sap2%N_independent_fcts
                        sat2 => sap2%sa_int(e2, f2, ipr)
                        DO r2 = 1, c2
                          n = n + 1
                          m = x1
                          DO f1 = 1, sap1%N_independent_fcts
                            sat1 => sap1%sa_int(e1, f1, ipr)
                            DO r1 = 1, c1
                              m = m + 1
                              DO a2 = 1, sat2%N_fcts
                                m2 = sat2%m(a2)
                                fct = sat2%c(a2) * wgt
                                DO a1 = 1, sat1%N_fcts
                                  m1 = sat1%m(a1)
                                  mat_ao(s1, s2)%m(r1+(m1-1)*c1                &
                                                 , r2+(m2-1)*c2, s) =          &
                                    mat_ao(s1, s2)%m(r1+(m1-1)*c1              &
                                                   , r2+(m2-1)*c2, s)          &
                                    + sat1%c(a1) * fct * mat_sy(irr)%m(m, n, s)
                                ENDDO
                              ENDDO
                            ENDDO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  end subroutine unsymadapt_spin

  !**************************************************************
  subroutine symadp3c(U1, E1, L1, U2, E2, L2, weight, cnt, sym, imode)
    ! add contribution of one not symmetry equivalent pair
    ! of equal atom to symmetry adapted integrals
    ! the case of total symmetric integrals
    use int_data_2cob3c_module, only: symadapt_totsym_3c_int_type
    use unique_atom_module, only: unique_atom_partner_type &
                                , unique_atom_sa_int_type  &
                                , unique_atoms
    use symmetry_data_module, only: symmetry_data_n_irreps &
                                  , symmetry_data_n_partners
    implicit none
    integer(IK)                 , intent(in)    :: U1,E1,L1,U2,E2,L2
    real(RK)                    , intent(in)    :: weight
    real(RK)                    , intent(in)    :: cnt(:,:,:,:,:) ! (nc2,nc1,:,2*L2+1,2*L1+1)
    type(symadapt_totsym_3c_int_type), intent(inout) :: sym(:)   ! (n_irreps)
    integer(IK), optional       , intent(in)    :: imode
    !------------ Declaration of local variables ---------------
    real(RK), pointer, dimension(:,:,:,:,:) :: sai
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(IK) :: i_ir, i_pa, i_if1, i_if2, &
                   i_cf1, i_cf2, m1, m2
    real(RK)    :: coef1, coef, weight2
    integer(IK) :: n_pa
    integer(IK) :: mode
    !------------ Executable code ------------------------------

    FPP_TIMER_START(sym3c)
    mode = 0
    if( present(imode) ) mode = imode

    irreps: do i_ir = 1, symmetry_data_n_irreps()
       sap1 => unique_atoms(U1)%symadapt_partner(i_ir,L1)
       sap2 => unique_atoms(U2)%symadapt_partner(i_ir,L2)

       sai => sym(i_ir)%int

       if( IAND(mode,IDONT_SUM_OVER_PARTNERS)/=0 )then
         ! all totally symmetric integrals dont depend
         ! on partner index, compute only i_pa=1:
         n_pa    = 1
         weight2 = weight
       else
         ! only symmetry inequivalent distances are used,
         ! wheited sum over partners is required:
!        weight = symequivatoms%weight(i_symequiv)
         n_pa    = symmetry_data_n_partners(i_ir)
         weight2 = weight / n_pa
       endif

       partners: do i_pa = 1, n_pa
          ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
             sat1 => sap1%sa_int(E1,i_if1,i_pa)
             ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                sat2 => sap2%sa_int(E2,i_if2,i_pa)
                m_sum_1: do i_cf1 = 1, sat1%N_fcts
                   m1 = sat1%m(i_cf1)
                   coef1 = sat1%c(i_cf1) * weight2
                   m_sum_2: do i_cf2 = 1, sat2%N_fcts
                      m2 = sat2%m(i_cf2)
                      coef = sat2%c(i_cf2) * coef1

                         sai(:, :, :, i_if2, i_if1) =  sai(:, :, :, i_if2, i_if1) &
                                        + coef * cnt(:, :, :, m2, m1)

                   enddo m_sum_2
                enddo m_sum_1
             enddo ind_fct_2
          enddo ind_fct_1
       enddo partners
    enddo irreps
    FPP_TIMER_STOP(sym3c)
  end subroutine symadp3c

  !*************************************************************
  subroutine integral_calc_quad_densmat(U1,L1,U2,L2,P,W)
    !
    !  Purpose: prepare the sections of the densmat P and energy-weighted
    !           densmat W that correspond to quadrupel (U1,L1,U2,L2)
    !
    !  Note:    P(irr,1) -- is a total density matrix
    !           P(irr,2) -- is a spin-density matrix available only if n_spin==2
    !
    use datatype, only: arrmat2
    use dimensions, only: dimens, dimoff
    use occupied_levels_module, only: eigvec_occ, eigval_occ, occ_num_occ
    use symmetry_data_module, only: ssym
    implicit none
    integer(IK), intent(in) :: U1,L1,U2,L2
    type(arrmat2), pointer  :: P(:,:) ! (n_irr,n_spin)
    type(arrmat2), pointer  :: W(:)   ! (n_irr)
    ! *** end of interface ***

    integer(IK), parameter :: CONTRACTED=1
    integer(IK)            :: dim1,off1
    integer(IK)            :: dim2,off2
    integer(IK)            :: n_irr,n_spin
    integer(IK)            :: irr,s,nu,mu,i_occ
    integer(IK)            :: stat
    logical                :: diagonal
    real(RK)               :: sums
    real(RK)               :: pa,pb

    diagonal = (U1 == U2) .and. (L1 == L2)

    n_irr  = ssym%n_irrep
    n_spin = ssym%n_spin

    ASSERT(.not.associated(P))
    ASSERT(.not.associated(W))
    allocate(P(n_irr,n_spin),W(n_irr),stat=stat)
    ASSERT(stat==0)

    do irr=1,n_irr

      dim1 = dimens(irr,U1,L1,CONTRACTED)
      off1 = dimoff(irr,U1,L1,CONTRACTED)
      dim2 = dimens(irr,U2,L2,CONTRACTED)
      off2 = dimoff(irr,U2,L2,CONTRACTED)

ASSERT(.not.allocated(w(irr)%m))
      allocate(W(irr)%m(dim2,dim1),stat=stat)
      ASSERT(stat==0)
      W(irr)%m(:,:) = 0

      stat = size(eigvec_occ(irr)%m,3)
      ASSERT(stat==n_spin)

      do s = 1,n_spin
ASSERT(.not.allocated(P(irr,s)%m))
        allocate(P(irr,s)%m(dim2,dim1),stat=stat)
        ASSERT(stat==0)
        P(irr,s)%m(:,:) = 0

        ! eigenvector sections corresponding to
        ! (irr,U1,L1) and (irr,U2,L2):
        ! eigenvalues for energy-weighted densmat:
        ! corresponding occupation numbers:
        associate( eigv1 => eigvec_occ(irr)%m(1+off1:dim1+off1,:,s)            &
                 , eigv2 => eigvec_occ(irr)%m(1+off2:dim2+off2,:,s)            &
                 , eigval => eigval_occ(irr)%m(:,s)                            &
                 , occ => occ_num_occ(irr)%m(:,s) )

ASSERT(size(eigval)==size(occ))
ASSERT(size(eigval)==size(eigv1,2))
ASSERT(size(eigval)==size(eigv2,2))

        do i_occ=1,size(eigval) ! occupied orbitals
          do nu=1,dim1
            do mu=1,dim2
              sums = occ(i_occ) * eigv2(mu,i_occ) * eigv1(nu,i_occ)

              if(.not.diagonal) sums = sums * 2

              P(irr,s)%m(mu,nu) = P(irr,s)%m(mu,nu) + sums
              W(irr)%m(mu,nu)   = W(irr)%m(mu,nu)   - sums * eigval(i_occ)
            enddo
          enddo
        enddo ! loop over occupied orbitals
        end associate
      enddo ! spin s

      if(n_spin==2)then
        ! go from (alpha,beta) to (total,spin):
        do nu=1,dim1
          do mu=1,dim2
            pa = P(irr,1)%m(mu,nu)
            pb = P(irr,2)%m(mu,nu)
            P(irr,1)%m(mu,nu) = pa + pb ! total densmat
            P(irr,2)%m(mu,nu) = pa - pb ! spin densmat
          enddo
        enddo
      endif
    enddo ! irrep irr
  end subroutine integral_calc_quad_densmat
  !*************************************************************

  !*************************************************************
  subroutine quad_density_mat(U1,E1,L1,U2,E2,L2,weight,densmat)
   !
   ! Purpose:  to provide density matrix in terms of primitive GTO for coulomb gradient matrixes
   !
   use int_data_2cob3c_module, only: n_c2,n_c1,n_m2,n_m1
   implicit none
   integer(IK), intent(in)  :: U1,E1,L1,U2,E2,L2
   real(RK)   , intent(in)  :: weight
   real(RK)   , intent(out) :: densmat(:,:,:,:) ! (nexp2,nexp1,2*L2+1,2*L1+1)
   ! *** end of interface ***

   real(RK) :: cnt(n_c2,n_c1,n_m2,n_m1)

ASSERT(size(densmat,3)==2*L2+1)
ASSERT(size(densmat,4)==2*L1+1)

ASSERT(associated(quad_pmat))

   ! transform densmat in MO representation to (contracted) AO representation:
   call unsymadp2c(U1,E1,L1,U2,E2,L2,weight,quad_pmat(:,1),cnt)

   ! transform densmat in contracted AO rep to uncontracted AO rep:
   call uncontract_2center(cnt,densmat)
  end subroutine quad_density_mat
  !*************************************************************

  !*************************************************************
  subroutine integral_calc_quad_close()
    !
    !  Purpose: deallocate global vars before proceeding to
    !           the next quadrupel
    implicit none
    ! *** end of interface ***

    integer(IK) :: irr,s,stat


    ASSERT(associated(quad_pmat))
    ASSERT(associated(quad_wmat))

    do irr=1,size(quad_wmat)
      deallocate(quad_wmat(irr)%m,stat=stat)
      ASSERT(stat==0)
      do s=1,size(quad_pmat,2)
        deallocate(quad_pmat(irr,s)%m,stat=stat)
        ASSERT(stat==0)
      enddo
    enddo
    deallocate(quad_wmat,quad_pmat,stat=stat)
    ASSERT(stat==0)
  end subroutine integral_calc_quad_close
  !*************************************************************

  !--------------- End of module ----------------------------------
end module integral_calc_quad_module

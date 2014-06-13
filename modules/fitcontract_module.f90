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
module  fitcontract_module
  !---------------------------------------------------------------
  !  Purpose: contains all routines and data necessary to
  !           perform fitcontractions on the 3-center-primitive
  !           integrals.
  !           The renormalisation of the integrals concerning
  !           fitfunctions and the mapping of fitfunctions to
  !           metaindex described in orbitalprojection_module
  !           is also done.
  !           Stores the result in  data structures in
  !           'int_data_2cob3c_module'.
  !
  !  PUBLIC ROUTINES:
  !       - fitcontract(flag,num,i_ua,cutoff,integral)
  !         main routine for perfomrming the fitcontraction
  !         on a three center primitive integral.
  !
  !  Module used by: ss_calculate,ls_calculate,ll_calculate...
  !
  !  References: ...
  !
  !  Author: FN
  !  Date: 6/96
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   3/97
  ! Description: The module is no also able to contract
  !              the gradients of 3c-Integrals. For this purpose
  !              the flag grad has been added.
  !              The summation over fitting functions is done
  !              automatically
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype, only: three_center_l_v2, arrmat5
  use unique_atom_module
  use iounitadmin_module
  use fit_coeff_module, only: coeff_charge, coeff_spin, coeff_xcmda, coeff_xcmda_ts, &
                              coeff_charge_eperef
  use output_module, only: output_int_fitcontract
  use options_module, only: options_spin_restricted, &
                            options_xcmode, xcmode_model_density
  use symmetry_data_module,only: get_totalsymmetric_irrep

  implicit none
  private         ! by default, all names are private
  save            ! by default, all data are to be kept
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public :: fitcontract,print_fcalloc
#ifdef WITH_EPE
  public :: fitcontract_eperef
#endif
  !!SB: should be only for resp, but...
  public :: fitcontract_v2

  !================================================================
  ! End of public interface of module
  !================================================================

!..............................................................................
! << OUTPUT ARRAYS >>
! ===================
!
! flags = "ch"
! ~~~~~~~~~~~~
! I1 = [ prim_i | f_k | prim_j ]
!
! F1 : prim_int_3c_co(1:n_aexp,1:n_bexp,1:n_ch,1:n_m1,1:n_m2)
!
! flag = "xc"
! ~~~~~~~~~~~~~~~~~~~~~
! I1 = < prim_i | g_k | prim_j >
!
! F1 : prim_int_3c_xc(1:n_aexp,1:n_bexp,1:n_ch,1:n_m1,1:n_m2)
!
! flag = "grads"
! ~~~~~~~~~~~~~~
! F1 = Sum(k) a_k,tot  [ d/dRa prim_i |          f_k |       prim_j ]   or
!      Sum(k) a_k,tot  [       prim_i |          f_k | d/dRb prim_j ]   or
!      Sum(k) a_k,tot  [       prim_i | dsym/dRc f_k |       prim_j ]
!
! F2 = Sum(k) b_k,tot  [ d/dRa prim_i |          f_k |       prim_j ]   or
!      Sum(k) b_k,tot  [       prim_i |          f_k | d/dRb prim_j ]   or
!      Sum(k) b_k,tot  [       prim_i | dsym/dRc f_k |       prim_j ]
!
! F3 = Sum(k) b_k,spin [ d/dRa prim_i |          f_k |       prim_j ]   or
!      Sum(k) b_k,spin [       prim_i |          f_k | d/dRb prim_j ]   or
!      Sum(k) b_k,spin [       prim_i | dsym/dRc f_k |       prim_j ]
!
! Storage
! ~~~~~~~
! PRESENT(mda_xcpot_gradients)                   | FALSE  FALSE  TRUE   TRUE
! options_xcmode() == xcmode_model_density       | FALSE  TRUE   FALSE  TRUE
! -----------------------------------------------+--------------------------
! gradients          (1:n_a,1:n_b,1:n_m1,1:n_m2) | F1     F1+F2  F1     F1
! mda_xcpot_gradients(1:n_a,1:n_b,1:n_m1,1:n_m2) | --     --     --     F2
! mda_spin_gradients (1:n_a,1:n_b,1:n_m1,1:n_m2) | --     [F3]   --     [F3]
!
! n_a : n_aexp , n_b : n_bexp , [F3] : F3 if options_n_spin() > 1
!
! << WORKING ARRAYS >>
! ====================
! intermediate(1:num,1:n_m1,1:n_m2,1:n_dim) : [ d/dR prim_i | f_k | prim_j ] or
! intermediate(1:num                      )   [ prim_i | d/dR f_k | prim_j ] or
!                                             [ prim_i | f_k | d/dR prim_j ]
!..............................................................................

  !------------ Declaration of constants and variables ----
  integer(kind=i4_kind)            :: counter
  logical                          :: do_glob_cons
  real(kind=r8_kind),parameter     :: zero=0.0_r8_kind
  type(unique_atom_type),  pointer :: ua
  type(unique_atom_basis_type), pointer :: uabl(:),uabr
!  type(unique_atom_renorm_indfct_type), pointer :: renormaliation_partner(:)
  type(unique_atom_glob_con_type), pointer :: glob_con(:)
  integer, pointer :: orbitalprojection(:,:)
  integer, pointer :: orbitalprojection_globcontr(:)
!!$  real(kind=r8_kind), pointer   :: prim_int_3c(:,:,:,:,:)
  integer :: lmax
  integer(kind=i4_kind):: fcalloc(4)=0
  logical :: model_density, spin_polarized

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************

  subroutine fitcontract_v2(num,i_ua,cutoff,integral)
    !
    !  Neither global nor local contraction used for coulomb fit fct
    !  Procedure just resend integral into prim_int_3c_co
    !
    ! author: SB
    ! date: 02/05
    !------------ Modules used -----------------------------------
    use integralpar_module
    use datatype
    use orbitalprojection_module
    use int_data_2cob3c_module, only: prim_int_3c_co
    use symmetry_data_module,only: symmetry_data_n_irreps,  &
         symmetry_data_n_partners,&
         get_totalsymmetric_irrep
    use ch_response_module, only: dimension_of_fit_ch
    use operations_module, only: operations_response
    use debug, only:show
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind)             :: num
    ! metaindex for exponents of first and secound unique atom
    integer(kind=i4_kind),intent(in)  :: i_ua
    ! index of equal atom processed
    logical,dimension(:,:),intent(in) :: cutoff
    ! this is the logical mask used to  unpack from the linear
    ! index back to the two dimensional array of exponents
    ! alpha and beta
    type(three_center_l_v2)              :: integral(:)
    ! the fitfunction-symmetry-adapted three center integral
    ! integral(i_ir)%l(-1:lmax_ch)%m(num,ncexps,n_independent_fcts,n_m1,n_m2,n_pa)
    ! with num beeing the metaindex for (i_exp1,i_exp2)
    ! [ the input integrals, spoiled on output! ]

    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
!!$    type(arrmat1int), allocatable :: proj_index_ch(:,:)

!!$    integer(kind=i4_kind)            :: i_uncont, &
!!$         i_l,i_ind,i_ir,counter_old,counter_ind,m1,m2, &
!!$         dim_intermediate,n_independent,i_a,i_b,i_ab, &
!!$         n_m1,n_m2, i_pa
!!$    integer(kind=i4_kind)            :: alloc_stat, i_ua_alloc
!!$    real(kind=r8_kind),allocatable        :: intermediate(:,:,:,:)
!!$    type(unique_atom_basis_type), pointer :: uab
!!$    real(kind=r8_kind),pointer            :: pointer_int(:,:,:,:,:)
    integer(kind=i4_kind), save           :: irrep_ffc(20,20)
    integer(kind=i4_kind)                 :: ffc, nexp, nind, LL, i_c, i_pa, i_ir
    integer(kind=i4_kind)                 :: i_l, i_ind, m2, m1, npa, nirr, dim
    !------------ Executable code --------------------------------


    if (operations_response) then
       nirr = symmetry_data_n_irreps()
       ASSERT(nirr<=20)
    else
       nirr = 1_i4_kind
    end if

!! debug start
!!$    do i_ir = 1, nirr
!!$       npa = symmetry_data_n_partners(i_ir)
!!$       do i_pa = 1, npa
!!$          do LL = -1,size(integral(i_ir)%l)-2
!!$             print *," i_ir = ",i_ir," i_pa = ",i_pa," LL= ", LL
!!$             call show("int",integral(i_ir)%l(LL)%m(1,1,1,:,:,i_pa))
!!$          end do
!!$       end do
!!$    end do
!! debug end

    if (i_ua==1) then
       ffc = 0
       do i_ir = 1, nirr
          npa = symmetry_data_n_partners(i_ir)
          ASSERT(npa<=20)
          dim = dimension_of_fit_ch(i_ir)
          do i_pa = 1, npa
             irrep_ffc(i_ir,i_pa) = ffc
!!$             print *,"irrep_ffc(",i_ir,",",i_pa,") = ",ffc
             ffc = ffc+dim
          end do
       end do
    end if

    ua => unique_atoms(i_ua)

    uabr => ua%r2_ch
    uabl => ua%l_ch

    lmax =  ua%lmax_ch

!!$    print *,'fitcontract_v2: resp entry'

    i_ir_: do i_ir = 1, nirr

       i_pa_: do i_pa=1,symmetry_data_n_partners(i_ir)


          ffc = irrep_ffc(i_ir,i_pa)


!!$          print *,"i_ir = ",i_ir," i_pa= ",i_pa," counter = ", ffc

          i_l_: do i_l = -1,lmax
!!$             print *, " i_l = ", i_l, " counter = ", ffc
             LL = i_l
             select case (i_l)
             case(-1) !! i_l = s
                nexp  = uabl(0)%N_uncontracted_fcts
                nind = ua%symadapt_partner(i_ir,0)%N_independent_fcts
                LL = 0
             case (0) !! i_l = r2
                nexp  = uabr%N_uncontracted_fcts
                nind = ua%symadapt_partner(i_ir,0)%N_independent_fcts
                LL = -1
             case default
                nexp  = uabl(i_l)%N_uncontracted_fcts
                nind = ua%symadapt_partner(i_ir,i_l)%N_independent_fcts
             end select

             do i_c = 1, nexp
                do i_ind = 1, nind
                   ffc = ffc + 1
!! debug:
!!$                   print *,"i_ir = ",i_ir," i_pa = ",i_pa," i_ind = ",i_ind, " i_l = ",LL
!!$                   print *,"shape int ",shape(integral(i_ir)%l(LL)%m)
!!$                   call show(" int ",integral(i_ir)%l(LL)%m(:,:,i_ind,1,1,i_pa))
!! end debug
                   do m2 = 1, size(prim_int_3c_co,4)
                      do m1 = 1, size(prim_int_3c_co,5)
                         prim_int_3c_co(:,:,ffc,m2,m1) =  &
                              unpack(integral(i_ir)%l(LL)%m(:,i_c,i_ind,m2,m1,i_pa),cutoff,0.0_r8_kind)
                      end do
                   end do
                end do
             end do

          end do i_l_
          irrep_ffc(i_ir,i_pa) = ffc
       end do i_pa_

    end do i_ir_

!!$    print *," shape(prim_int_3c_co)=",shape(prim_int_3c_co)
!!$    do i_l = 1,size(prim_int_3c_co,3)
!!$       print *," ff= ",i_l
!!$       call show('prim_int_3c_co',prim_int_3c_co(:,:,i_l,1,1))
!!$    end do

  end subroutine fitcontract_v2

  subroutine fitcontract(flag,num,i_ua,cutoff,integral,gradients, &
             mda_spin_gradients,mda_xcpot_gradients,cpks_gradients)
    !  Purpose: main routine for performing fit_contractions:
    !
    !  (1) local contractions : contraction of exponents on
    !      one specific unique atom and for one specific angular
    !      momentum to a single exponent.
    !      These contractions are done first for all uniques.
    !
    !  (2) global contractions : contraction of exponents on
    !      diverse unique atoms and diverse angular momenta
    !      including different independent functions to a
    !      single exponent.
    !      These contractions are appended to the local ones in
    !      in the sense that we first loop over the uniques to
    !      do the local contraction and then loop again over
    !      the uniques to do the global contractions.
    !
    !  The renormalisation of the integrals concerning
    !  fitfunctions and the mapping of fitfunctions to
    !  metaindex described in orbitalprojection_module
    !  is also done
    !  Stores the result in  data structures in
    !  'int_data_2cob3c_module'.
    !
    ! author: FN
    ! date: 7/96
    !------------ Modules used -----------------------------------
    use orbitalprojection_module
    use int_data_2cob3c_module, only: prim_int_3c_co, prim_int_3c_xc
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(*),intent(in)           :: flag ! 'ch' or 'xc' or 'grad'
    integer(kind=i4_kind)             :: num
      ! metaindex for exponents of first and secound unique atom
    integer(kind=i4_kind),intent(in)  :: i_ua
      ! index of equal atom processed
    logical,dimension(:,:),intent(in) :: cutoff
      ! this is the logical mask used to  unpack from the linear
      ! index back to the two dimensional array of exponents
      ! alpha and beta
    type(arrmat5) , target               :: integral(-1:) ! (-1:lmax_xx)
      ! the fitfunction-symmetry-adapted three center integral
      ! integral(-1:lmax_xx)%m(num,ncexps,n_independent_fcts,n_m1,n_m2)
      ! with num beeing the metaindex for (i_exp1,i_exp2)
      ! [ the input integrals, spoiled on output! ]
    real(kind=r8_kind),optional :: gradients(:,:,:,:), &
         mda_spin_gradients(:,:,:,:), mda_xcpot_gradients(:,:,:,:)
#ifdef no_cpks_coul_grads
    real(kind=r8_kind),optional, intent(inout) :: cpks_gradients(:)
#else
    real(kind=r8_kind),optional,pointer :: cpks_gradients(:,:,:,:,:)
#endif
      ! only used for 'grad' case [ the output integrals ]
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                   :: i_glob,N_glob_cons
    real(kind=r8_kind), pointer   :: prim_int_3c(:,:,:,:,:)
    !------------ Executable code --------------------------------
    ua => unique_atoms(i_ua)
    select case ( flag )
    case ( 'ch' )
       uabr => ua%r2_ch
       uabl => ua%l_ch
       glob_con => ua%glob_con_ch
       orbitalprojection => orbitalprojection_ch
       orbitalprojection_globcontr => orbitalprojection_globcontr_ch
       prim_int_3c => prim_int_3c_co
       do_glob_cons = ua%N_glob_cons_ch.gt.0
       N_glob_cons = ua%N_glob_cons_ch
       lmax = ua%lmax_ch
       if (output_int_fitcontract) call write_to_output_units &
            ("FIT_CONTRACT : processing Coulomb integral")
    case ( 'xc' )
       uabr => ua%r2_xc
       uabl => ua%l_xc
       glob_con => ua%glob_con_xc
       orbitalprojection => orbitalprojection_xc
       orbitalprojection_globcontr => orbitalprojection_globcontr_xc
       prim_int_3c => prim_int_3c_xc
       do_glob_cons = ua%N_glob_cons_xc.gt.0
       N_glob_cons = ua%N_glob_cons_xc
       lmax = ua%lmax_xc
       if (output_int_fitcontract)  call write_to_output_units &
            ("FIT_CONTRACT : processing Exchange integral")

    case ( 'grad' )
       uabr => ua%r2_ch
       uabl => ua%l_ch
       glob_con => ua%glob_con_ch
       orbitalprojection => orbitalprojection_ch
       orbitalprojection_globcontr => orbitalprojection_globcontr_ch
       do_glob_cons = ua%N_glob_cons_ch.gt.0
       N_glob_cons = ua%N_glob_cons_ch
       lmax = ua%lmax_ch
       model_density = options_xcmode() == xcmode_model_density
       spin_polarized = .not. options_spin_restricted()
       if (output_int_fitcontract) &
            call write_to_output_units  ("FIT_CONTRACT : processing Coulomb integral")

       counter = orbitalprojection(-1,i_ua)
       if (output_int_fitcontract)  call write_to_output_units &
            ("FIT_CONTRACT : calling contract_loc")
       mda: if (model_density)then
          if (spin_polarized) then
             if (present(mda_xcpot_gradients)) then
                call contract_loc(integral,i_ua,num,cutoff, &
                     gradients=gradients, &
                     spin_grads=mda_spin_gradients, &
                     xcpot_grads=mda_xcpot_gradients)
             else
                call contract_loc(integral,i_ua,num,cutoff, &
                     gradients=gradients, &
                     spin_grads=mda_spin_gradients)
             endif
          else
             if (present(mda_xcpot_gradients)) then
                call contract_loc(integral,i_ua,num,cutoff, &
                     gradients=gradients, &
                     xcpot_grads=mda_xcpot_gradients)
             else
                call contract_loc(integral,i_ua,num,cutoff, &
                     gradients=gradients)
             endif
          endif
       else ! i.e.  normal nonmda case
          if(present(cpks_gradients)) then
             call contract_loc(integral,i_ua,num,cutoff, &
                  gradients=gradients, &
                  cpks_gradients=cpks_gradients)         !(1)
          else
             call contract_loc(integral,i_ua,num,cutoff, &
                  gradients=gradients)
          endif
       endif mda

       if (do_glob_cons) then
          counter = orbitalprojection_globcontr(i_ua)
          if (output_int_fitcontract) &
               write(*,*) &
               "FIT_CONTRACT : counter first global contraction", counter
          if (output_int_fitcontract) &
               call write_to_output_units("FIT_CONTRACT : calling contract_glob")
          if (model_density)then
             if (spin_polarized) then
                if (present(mda_xcpot_gradients)) then
                   do i_glob=1,N_glob_cons
                      call contract_glob(integral,i_glob,num,cutoff, &
                           gradients=gradients, &
                           spin_grads=mda_spin_gradients, &
                           xcpot_grads=mda_xcpot_gradients )
                   enddo
                else
                   do i_glob=1,N_glob_cons
                      call contract_glob(integral,i_glob,num,cutoff, &
                           gradients=gradients, &
                           spin_grads=mda_spin_gradients)
                   enddo
                endif
             else
                if (present(mda_xcpot_gradients)) then
                   do i_glob=1,N_glob_cons
                      call contract_glob(integral,i_glob,num,cutoff, &
                           gradients=gradients, &
                           xcpot_grads=mda_xcpot_gradients)
                   enddo
                else
                   do i_glob=1,N_glob_cons
                      call contract_glob(integral,i_glob,num,cutoff, &
                           gradients=gradients)
                   enddo
                endif
             endif
          else ! normal case
          do i_glob=1,N_glob_cons
             call contract_glob(integral,i_glob,num,cutoff, &
                  gradients=gradients)
          enddo
          endif
       endif
       return
    end select

    counter = orbitalprojection(-1,i_ua)
    if (output_int_fitcontract) &
         call write_to_output_units("FIT_CONTRACT : calling contract_loc")
    call contract_loc(integral,i_ua,num,cutoff, &
         prim_int_3c)

    if (do_glob_cons) then
       counter = orbitalprojection_globcontr(i_ua)
       if (output_int_fitcontract) &
            write(*,*) &
            "FIT_CONTRACT : counter first global contraction", counter
       if (output_int_fitcontract) &
            call write_to_output_units &
            ("FIT_CONTRACT : calling contract_glob")
       do i_glob=1,N_glob_cons
          call contract_glob(integral,i_glob,num,cutoff, &
               prim_int_3c)
       enddo
    endif

  end subroutine fitcontract


  subroutine contract_loc(integral,i_ua,num,cutoff,&
                         prim_int_3c, &
                         gradients, &
                         spin_grads, &
                         xcpot_grads, &
                         cpks_gradients)
    ! purpose: perform the local contraction only and
    !          assign the result directory to the
    !          GreatDataStructure 'prim_int_3c' by
    !          unpacking it using the logical mask 'cutoff'
    !
    !          As as security measure, the variable 'counter'
    !          is cross-checked with 'orbitalprojection'
    !          which contains the beginning index of the
    !          (scf-numbered) index for the fitfunctions
    !          for a given unique and angular momentum.
    !          Please pay attention to the dreadful order:
    !          In orbitalprojection_module l=-1 means s-type
    !          and l=0 means r2-type.
    !
    ! routine called by: fitcontract
    ! author: FN
    ! date: 6/96
    !------------ Modules used ------------------ ---------------
    use symmetry_data_module, only: get_totalsymmetric_irrep
    use machineparameters_module, only: machineparameters_DimCheck
    use int_data_2cob3c_module, only: prim_int_dens_mat
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(arrmat5), target            :: integral(-1:) ! (-1:lmax_xx)
    integer(kind=i4_kind),intent(in) :: i_ua
    integer(kind=i4_kind),intent(in) :: num
    logical,intent(in)               :: cutoff(:,:)
    real(kind=r8_kind),intent(inout) :: prim_int_3c(:,:,:,:,:)
    optional                         :: prim_int_3c ! not for the gradient part!
    real(kind=r8_kind),optional :: gradients(:,:,:,:), &
     spin_grads(:,:,:,:), xcpot_grads(:,:,:,:)
#ifdef no_cpks_coul_grads
 real(kind=r8_kind),optional, intent(inout):: cpks_gradients(:)
 integer(kind=i4_kind):: k_fit
#else
  real(kind=r8_kind),optional,pointer :: cpks_gradients(:,:,:,:,:)
#endif
    ! only present for gradient case
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)            :: i_uncont,i_cont,i_exp, &
         i_l, i_ind, i_ir, counter_old, counter_ind, m1, m2, &
         l_max,l_min,dim_intermediate,n_independent,i_a,i_b,i_ab, &
         n_m1,n_m2
    integer(kind=i4_kind)            :: i_fit
    real(kind=r8_kind),allocatable   :: intermediate(:,:,:,:)
    real(kind=r8_kind),allocatable   :: temp(:)
    real(kind=r8_kind),pointer       :: contract(:,:)
    type(unique_atom_basis_type),         pointer :: uab
    real(kind=r8_kind),pointer :: pointer_int(:,:,:,:,:)
    external error_handler
    !------------ Executable code --------------------------------

    ! first keep the counter
    counter_old = counter
    !allocate space for the intermediate variable
    i_ir = get_totalsymmetric_irrep()
    l_max=size(integral)
    l_min = 1
    n_m2 = size(integral(0)%m,5)
    n_m1 = size(integral(0)%m,4)

    ! calculate how much to allocate for intermediate
    dim_intermediate = uabl(0)%N_uncontracted_fcts + &
         uabl(0)%N_contracted_fcts +     &
         uabr%N_uncontracted_fcts +     &
         uabr%N_contracted_fcts

    do i_l=1,lmax
       dim_intermediate = dim_intermediate + &
            ua%symadapt_partner(i_ir,i_l)%N_independent_fcts * &
            ( uabl(i_l)%N_uncontracted_fcts +          &
            uabl(i_l)%N_contracted_fcts )
    enddo

    allocate( &
         intermediate( num, n_m1, n_m2, counter_old:counter_old+dim_intermediate ), &
         STAT=fcalloc(2)) !contract_loc
    ASSERT(fcalloc(2).eq.0)
    fcalloc(2)=1
    if (model_density) then
       allocate( temp(num), STAT=fcalloc(4))
       ASSERT(fcalloc(4).eq.0)
       fcalloc(4)=1
    end if

    counter = counter_old

    intermediate = 0.0_r8_kind

    ! first do s-type functions : this must be the order
    ! in which they do it in the old code ----------------------
    ! ATTENTION : in the variable orbitalprojection
    !             l=0  r2-type
    !             l=-1 s-type      (wuerg!)
    if (machineparameters_DimCheck) then
       if (counter.ne.orbitalprojection(-1,i_ua)) &
            call error_handler("contract_loc : something fishy (1)")
    endif

    uab => uabl(0)
    contract => uab%contractions
      pointer_int => integral(0)%m
      if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter l = 0 start", counter
    counter = counter + uab%N_uncontracted_fcts
    if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter l = 0 first contraction", counter
    do i_cont = 1,uab%N_contracted_fcts
       do i_exp=1,uab%N_exponents
          intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
               contract(i_exp,i_cont) * pointer_int(:,i_exp,1,:,:)
       enddo

       counter = counter + 1
    enddo

    counter = counter - uab%N_uncontracted_fcts - uab%N_contracted_fcts
    if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter l = 0 first primitive", counter
    do i_uncont = 1,uab%N_uncontracted_fcts
       intermediate(:,:,:,counter) = pointer_int(:,i_uncont,1,:,:)
       counter = counter + 1
    enddo
    counter = counter + uab%N_contracted_fcts
    ! ----------- end of s-type --------------------------------

    ! now do the r2-type (as soon as debugging is over we can
    ! drop the awkward order again) ----------------------------
    if (machineparameters_DimCheck) then
       if (counter.ne.orbitalprojection(0,i_ua)) &
            call error_handler &
            ("contract_loc : something fishy (2)")
    endif

    uab => uabr
    contract => uab%contractions
    pointer_int => integral(-1)%m
    if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter r**2 start", counter
    counter = counter + uab%N_uncontracted_fcts
    if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter r**2 first contracted", counter
    do i_cont = 1,uab%N_contracted_fcts
       do i_exp=1,uab%N_exponents
          intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
               contract(i_exp,i_cont) * pointer_int(:,i_exp,1,:,:)
       enddo

       counter = counter + 1
    enddo

    counter = counter - uab%N_uncontracted_fcts - uab%N_contracted_fcts
    if (output_int_fitcontract) &
         write(*,*) &
         "FIT_CONTRACT : counter r**2 first primitive", counter
    do i_uncont = 1,uab%N_uncontracted_fcts
       intermediate(:,:,:,counter) = pointer_int(:,i_uncont,1,:,:)
       counter = counter + 1
    enddo
    counter = counter + uab%N_contracted_fcts
    ! ---------- end of r2-type --------------------------------

    ! loop over angular momenta
    angular_mom: do i_l = 1,lmax
       pointer_int => integral(i_l)%m
       n_independent = &
            ua%symadapt_partner(i_ir,i_l)%n_independent_fcts

       if (machineparameters_DimCheck) then
          if (counter.ne.orbitalprojection(i_l,i_ua)) &
               call error_handler &
               ("contract_loc : something fishy (3)")
       endif

       uab => uabl(i_l)
       if (output_int_fitcontract) &
            write(*,*) &
            "FIT_CONTRACT : counter l = ",i_l," start", counter
       counter_ind = counter

       ! loop over independent fcts
       do i_ind = 1,n_independent
          counter = counter_ind + i_ind - 1
          if (output_int_fitcontract) &
               write(*,*) &
               "FIT_CONTRACT : counter l = ",i_l," i_ind = ",i_ind," start", counter

          ! local contract before renormalisation and renormalise result
          contract => uab%contractions
          counter = counter + uab%N_uncontracted_fcts * n_independent
          if (output_int_fitcontract) &
               write(*,*) &
               "FIT_CONTRACT : counter l = ",i_l," i_ind = ",i_ind," first contacted", counter
          do i_cont = 1,uab%N_contracted_fcts
             do i_exp=1,uab%N_exponents
                intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
                     contract(i_exp,i_cont) * pointer_int(:,i_exp,i_ind,:,:)
             enddo

             counter = counter + n_independent
          enddo

          counter = counter - n_independent * &
               (uab%N_uncontracted_fcts + uab%N_contracted_fcts)
          if (output_int_fitcontract) &
               write(*,*) &
               "FIT_CONTRACT : counter l = ",i_l," i_ind = ",i_ind," first primitive", counter
          ! pick those which remain uncontracted: these are
          ! assumed to be the FIRST N_uncontracted functions for
          ! each primitive fitbasis
          do i_uncont = 1,uab%N_uncontracted_fcts
             intermediate(:,:,:,counter) =  &
                  pointer_int(:,i_uncont,i_ind,:,:)
             counter = counter + n_independent
          enddo
          counter = counter + uab%N_contracted_fcts * n_independent
       enddo

       if (n_independent.ne.0) counter = counter - n_independent + 1

    enddo angular_mom
    if (output_int_fitcontract) &
         write(*,*) "FIT_CONTRACT : counter end", counter

    if(present(gradients)) then ! gradient case
       ! do summation over fitfunctions
!      m2_loop: do m2 = 1, n_m2
!         m1_loop: do m1 = 1, n_m1

       ! Intel compiler does not like
       !      if(present(optional))...
       ! inside of loops. Especially if they are .not.present!

             mda: if (model_density) then
                if (present(xcpot_grads)) then
                   ! ADD Sum(k) a_k,tot [ ... f_k ... ]
                  do m2 = 1, n_m2
                  do m1 = 1, n_m1
                   temp = zero
                   do i_fit=counter_old,counter-1
                      temp = temp + coeff_charge(i_fit) * &
                           intermediate(:,m1,m2,i_fit)
                   end do
                   gradients(:,:,m1,m2) = gradients(:,:,m1,m2) + &
                        unpack(temp,cutoff,zero)
                   ! ADD Sum(k) b_k,tot [ ... f_k ... ]
                   temp = zero
                   do i_fit=counter_old,counter-1
!                     temp = temp + coeff_xcmda(i_fit,1) * &
                      temp = temp + coeff_xcmda_ts(i_fit,1) * &
                           intermediate(:,m1,m2,i_fit)
                   end do
                   xcpot_grads(:,:,m1,m2) = xcpot_grads(:,:,m1,m2) + &
                        unpack(temp,cutoff,zero)
                  enddo
                  enddo
                else
                   ! ADD Sum(k) ( a_k,tot + b_k,tot ) [ ... f_k ... ]
                  do m2 = 1, n_m2
                  do m1 = 1, n_m1
                   temp = zero
                   do i_fit=counter_old,counter-1
                      temp = temp + ( coeff_charge(i_fit) + &
!                          coeff_xcmda(i_fit,1) ) * intermediate(:,m1,m2,i_fit)
                           coeff_xcmda_ts(i_fit,1) ) * intermediate(:,m1,m2,i_fit)
                   end do
                   gradients(:,:,m1,m2) = gradients(:,:,m1,m2) + &
                        unpack(temp,cutoff,zero)
                  enddo
                  enddo
                end if
                if (spin_polarized) then
                   ! ADD Sum(k) b_k,spin [ ... f_k ... ]
                  do m2 = 1, n_m2
                  do m1 = 1, n_m1
                   temp = zero
                   do i_fit=counter_old,counter-1
!                     temp = temp + coeff_xcmda(i_fit,2) * &
                      temp = temp + coeff_xcmda_ts(i_fit,2) * &
                           intermediate(:,m1,m2,i_fit)
                   end do
                   spin_grads(:,:,m1,m2) = spin_grads(:,:,m1,m2) + &
                        unpack(temp,cutoff,zero)
                  enddo
                  enddo
                endif

             else mda
             ! ADD Sum(k) a_k,tot [ ... f_k ... ]  (old fashion technique)
            do m2 = 1, n_m2
            do m1 = 1, n_m1
             i_ab = 1
             do i_b = 1, size(gradients,2)
                do i_a  = 1, size(gradients,1)
                   if ( cutoff(i_a,i_b) ) then
                      gradients(i_a,i_b,m1,m2) = gradients(i_a,i_b,m1,m2) + &
                           sum( coeff_charge(counter_old:counter-1) * &
                           intermediate(i_ab,m1,m2,counter_old:counter-1) )
                      i_ab = i_ab + 1
                   endif
                enddo
             enddo
            enddo
            enddo
             if(present(cpks_gradients)) then
            do m2 = 1, n_m2
            do m1 = 1, n_m1
             i_ab = 1
             do i_b = 1, size(gradients,2)
                do i_a  = 1, size(gradients,1)
                   if ( cutoff(i_a,i_b) ) then
#ifdef no_cpks_coul_grads
             do k_fit=counter_old,counter-1
              cpks_gradients(k_fit)=cpks_gradients(k_fit) &
             +intermediate(i_ab,m1,m2,k_fit)*prim_int_dens_mat(i_a,i_b,m1,m2)
             enddo
#else
                      cpks_gradients(i_a,i_b,counter_old:counter-1,m1,m2) =  &
                       cpks_gradients(i_a,i_b,counter_old:counter-1,m1,m2) + &
                                intermediate(i_ab,m1,m2,counter_old:counter-1)
#endif
                      i_ab = i_ab + 1
                   endif
                enddo
             enddo
            enddo
            enddo
             endif

             endif mda

!         enddo m1_loop
!      enddo m2_loop

    else !'normal' non gradient case
       ASSERT(present(prim_int_3c))
       i_ab = 1
       do i_b = 1, size(prim_int_3c,2)
          do i_a  = 1, size(prim_int_3c,1)
             if ( cutoff(i_a,i_b) ) then
                do m2 = 1, n_m2
                   do m1 = 1, n_m1
                      prim_int_3c(i_a,i_b,counter_old:counter-1,m1,m2) = &
                           intermediate(i_ab,m1,m2,counter_old:counter-1)
                   enddo
                enddo
                i_ab = i_ab + 1
             else
                prim_int_3c(i_a,i_b,counter_old:counter-1,:,:) = zero
             endif
          enddo
       enddo
    end if

    if (model_density) then
       deallocate( temp, STAT=fcalloc(4))
       ASSERT(fcalloc(4).eq.0)
    end if

    deallocate(intermediate,STAT=fcalloc(2))
    if (fcalloc(2).ne.0 ) call error_handler &
         ("contract_loc : deallocation (1) failed")

  end subroutine contract_loc

  !*************************************************************

  subroutine contract_glob(integral,i_glob,num,cutoff,&
                           prim_int_3c, &
                           gradients, &
                           spin_grads, &
                           xcpot_grads)
    ! purpose: perform the global contraction for a given
    !          unique atom and assign the result to the
    !          GreatDataStructure in int_data_2cob3c_module.
    !          The index of the fitfunction (as it is used
    !          in the SCF-part is taken from 'counter'.
    !
    ! routine called by: fit_contract
    ! author: FN
    ! date: 6/96
    !------------ Modules used ------------------ ---------------
    use int_data_2cob3c_module
    use symmetry_data_module, only: get_totalsymmetric_irrep
    use iounitadmin_module, only: write_to_output_units
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(arrmat5), intent(in)        :: integral(-1:) ! (-1:lmax_xx)
    integer(kind=i4_kind),intent(in) :: i_glob
    integer(kind=i4_kind),intent(in) :: num
    logical,intent(in)               :: cutoff(:,:)
    real(r8_kind), intent(inout)     :: prim_int_3c(:,:,:,:,:)
    optional                         :: prim_int_3c ! not for the gradient part
    real(kind=r8_kind),optional      :: gradients(:,:,:,:), &
         spin_grads(:,:,:,:), xcpot_grads(:,:,:,:)
    !** End of interface *****************************************
    !------------- Declaration of local variables ----------------
    integer(kind=i4_kind)            :: n_contributing, &
         i_cont, i_ind, i_exp, i_l ,m1, m2, &
         min_m1, max_m1, min_m2,max_m2,i_ir
    real(kind=r8_kind)               :: coef
    real(kind=r8_kind),allocatable   :: intermediate(:)
    external error_handler
    !------------ Executable code --------------------------------

    ! first get some things from the 'unique_atom_module'
    n_contributing = &
         glob_con(i_glob)%N_contributing_fcts

    i_ir = get_totalsymmetric_irrep()
    min_m1 = 1
    max_m1 = size(integral(0)%m,4)
    min_m2 = 1
    max_m2 = size(integral(0)%m,5)

    allocate( intermediate(num),STAT=fcalloc(3)) ! contract_glob
    ASSERT(fcalloc(3).eq.0)
    fcalloc(3)=1

    m2_loop: do m2 = min_m2, max_m2
       m1_loop: do m1 = min_m1, max_m1

          intermediate = zero

          contributing: do i_cont = 1,n_contributing

             ! ATTENTION i_l = -1 means s-type
             !           i_l = 0  means r2-type
             i_l = glob_con(i_glob)%l(i_cont)
             if (i_l .eq. -1) then
                i_l = 0
             elseif (i_l .eq. 0 ) then
                i_l = -1
             endif
             i_ind = glob_con(i_glob)%index_ind_fct(i_cont)
             i_exp = glob_con(i_glob)%index_exp(i_cont)
             coef = glob_con(i_glob)%coefs(i_cont)
             intermediate(:) = intermediate(:) + &
                  integral(i_l)%m(:,i_exp,i_ind,m1,m2) * coef
          enddo contributing

          ! now re-map
          if(present(gradients)) then ! gradient case
             if (model_density) then
                if (present(xcpot_grads)) then
                   ! ADD Sum(k) a_k,tot [ ... f_k ... ]
                   gradients(:,:,m1,m2) = gradients(:,:,m1,m2) + &
                        unpack( coeff_charge(counter)*intermediate(:), &
                        cutoff,zero)
                   ! ADD Sum(k) b_k,tot [ ... f_k ... ]
                   xcpot_grads(:,:,m1,m2) = xcpot_grads(:,:,m1,m2) + &
!                       unpack( coeff_xcmda(counter,1)*intermediate(:), &
                        unpack( coeff_xcmda_ts(counter,1)*intermediate(:), &
                        cutoff,zero)
                else
                   ! ADD Sum(k) ( a_k,tot + b_k,tot ) [ ... f_k ... ]
                   gradients(:,:,m1,m2) = gradients(:,:,m1,m2) + &
                        unpack( ( coeff_charge(counter) + &
!                       coeff_xcmda(counter,1) )*intermediate(:), &
                        coeff_xcmda_ts(counter,1) )*intermediate(:), &
                        cutoff,zero)
                end if
                if (spin_polarized) then
                   ! ADD Sum(k) b_k,spin [ ... f_k ... ]
                   spin_grads(:,:,m1,m2) = spin_grads(:,:,m1,m2) + &
!                       unpack( coeff_xcmda(counter,2)*intermediate(:), &
                        unpack( coeff_xcmda_ts(counter,2)*intermediate(:), &
                        cutoff,zero)
                endif
             else
             ! ADD Sum(k) a_k,tot [ ... f_k ... ]  (old fashion technique)

             ! do summation over fitfunctions
             gradients(:,:,m1,m2) = &
                  gradients(:,:,m1,m2)+coeff_charge(counter)*&
                  unpack(intermediate(:),cutoff,zero)
             endif
          else !'normal' case
             ASSERT(present(prim_int_3c))
             prim_int_3c(:,:,counter,m1,m2)= &
                  unpack(intermediate(:),cutoff,zero)
          end if
       enddo m1_loop
    enddo m2_loop

    counter = counter + 1

    deallocate(intermediate,STAT=fcalloc(3))
    if (fcalloc(3).ne.0 ) call error_handler &
         ("contract_glob : deallocation (1) failed")

  end subroutine contract_glob

#ifdef WITH_EPE
  subroutine fitcontract_eperef(flag,num,i_ua,cutoff,integral,gradients, &
                         mda_spin_gradients,mda_xcpot_gradients)
    !  Purpose: main routine for performing fit_contractions:
    !
    !  (1) local contractions : contraction of exponents on
    !      one specific unique atom and for one specific angular
    !      momentum to a single exponent.
    !      These contractions are done first for all uniques.
    !
    !  (2) global contractions : contraction of exponents on
    !      diverse unique atoms and diverse angular momenta
    !      including different independent functions to a
    !      single exponent.
    !      These contractions are appended to the local ones in
    !      in the sense that we first loop over the uniques to
    !      do the local contraction and then loop again over
    !      the uniques to do the global contractions.
    !
    !  The renormalisation of the integrals concerning
    !  fitfunctions and the mapping of fitfunctions to
    !  metaindex described in orbitalprojection_module
    !  is also done
    !  Stores the result in  data structures in
    !  'int_data_2cob3c_module'.
    !
    !  reshaped from fit_contract by VN
    !------------ Modules used -----------------------------------
    use orbitalprojection_module
    use int_data_2cob3c_module, only: prim_int_3c_co, prim_int_3c_xc
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(*),intent(in)           :: flag ! 'ch' or 'xc' or 'grad'
    integer(kind=i4_kind)             :: num
      ! metaindex for exponents of first and secound unique atom
    integer(kind=i4_kind),intent(in)  :: i_ua
      ! index of equal atom processed
    logical,dimension(:,:),intent(in) :: cutoff
      ! this is the logical mask used to  unpack from the linear
      ! index back to the two dimensional array of exponents
      ! alpha and beta
    type(arrmat5), target              :: integral(-1:) ! (-1:lmax_xx)
      ! the fitfunction-symmetry-adapted three center integral
      ! integral(-1:lmax_xx)%m(num,ncexps,n_independent_fcts,n_m1,n_m2)
      ! with num beeing the metaindex for (i_exp1,i_exp2)
      ! [ the input integrals, spoiled on output! ]
!    real(kind=r8_kind),poiner,dimension(:):: coeff_charge_current
    real(kind=r8_kind),optional :: gradients(:,:,:,:), &
         mda_spin_gradients(:,:,:,:), mda_xcpot_gradients(:,:,:,:)
      ! only used for 'grad' case [ the output integrals ]
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                   :: N_glob_cons
!!$    integer(kind=i4_kind)                   :: i_glob
    !------------ Executable code --------------------------------
    ua => unique_atoms_eperef(i_ua)
    orbitalprojection => orbitalprojection_ch_eperef
       uabr => ua%r2_ch
       uabl => ua%l_ch

       do_glob_cons = ua%N_glob_cons_ch.gt.0
       N_glob_cons = ua%N_glob_cons_ch
       if(do_glob_cons) glob_con => ua%glob_con_ch
       lmax = ua%lmax_ch
       model_density = options_xcmode() == xcmode_model_density
       spin_polarized = .not. options_spin_restricted()
       counter = orbitalprojection(-1,i_ua)
       call contract_loc_eperef(integral,i_ua,num,cutoff,gradients)
       if (do_glob_cons) then
        print*,'Global contraction of EPE relaxation do not implemented'
        stop
       endif
       return

    counter = orbitalprojection(-1,i_ua)
    call contract_loc(integral,i_ua,num,cutoff)
    if (do_glob_cons) then
        print*,'Global contraction of EPE relaxation do not implemented'
        stop
    endif

  end subroutine fitcontract_eperef

  subroutine contract_loc_eperef(integral,i_ua,num,cutoff,gradients, &
                          spin_grads,xcpot_grads)
    ! purpose: perform the local contraction only and
    !          assign the result directory to the
    !          GreatDataStructure 'prim_int_3c' by
    !          unpacking it using the logical mask 'cutoff'
    !
    !          As as security measure, the variable 'counter'
    !          is cross-checked with 'orbitalprojection'
    !          which contains the beginning index of the
    !          (scf-numbered) index for the fitfunctions
    !          for a given unique and angular momentum.
    !          Please pay attention to the dreadful order:
    !          In orbitalprojection_module l=-1 means s-type
    !          and l=0 means r2-type.
    !
    ! routine called by: fitcontract
    ! author: FN
    ! date: 6/96
    !------------ Modules used ------------------ ---------------
    use symmetry_data_module, only: get_totalsymmetric_irrep
    use machineparameters_module, only: machineparameters_DimCheck
    use epecom_module, only: i_ir_eperef
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(arrmat5), target            :: integral(-1:) ! (-1:lmax_ch)
    integer(kind=i4_kind),intent(in) :: i_ua
    integer(kind=i4_kind),intent(in) :: num
    logical,intent(in)               :: cutoff(:,:)
    real(kind=r8_kind),optional :: gradients(:,:,:,:), &
         spin_grads(:,:,:,:), xcpot_grads(:,:,:,:)
    ! only present for gradient case
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)            :: i_uncont,i_cont,i_exp, &
         i_l, i_ind, i_ir, counter_old, counter_ind, m1, m2, &
         l_max,l_min,dim_intermediate,n_independent,i_a,i_b,i_ab, &
         n_m1,n_m2
!!$ integer(kind=i4_kind)            :: i_fit
    real(kind=r8_kind),allocatable   :: intermediate(:,:,:,:)
!!$ real(kind=r8_kind),allocatable   :: temp(:)
    real(kind=r8_kind),pointer       :: contract(:,:)
    type(unique_atom_basis_type),         pointer :: uab
    real(kind=r8_kind),pointer :: pointer_int(:,:,:,:,:)
    external error_handler
    !------------ Executable code --------------------------------

    ! first keep the counter
    counter_old = counter
    !allocate space for the intermediate variable
!    i_ir = get_totalsymmetric_irrep()

    i_ir= i_ir_eperef
    l_max=size(integral) ! coul_int_c allocated for ref case
    l_min = 1
    n_m2 = size(integral(0)%m,5)
    n_m1 = size(integral(0)%m,4)

    ! calculate how much to allocate for intermediate
    ! ua => unique_atoms_eperef(i_ua)
    dim_intermediate = & ! unique_atoms_eperef(i_ua)$l_ch(0)
       uabl(0)%N_uncontracted_fcts+uabl(0)%N_contracted_fcts +     &
             uabr%N_uncontracted_fcts + uabr%N_contracted_fcts

    do i_l=1,lmax
       dim_intermediate = dim_intermediate + &
         ua%symadapt_partner(i_ir,i_l)%N_independent_fcts * &
          ( uabl(i_l)%N_uncontracted_fcts + uabl(i_l)%N_contracted_fcts )
    enddo

    allocate( intermediate &
              (num,n_m1,n_m2,counter_old:counter_old+dim_intermediate), &
                                                         stat=fcalloc(1))
    ASSERT(fcalloc(1).eq.0)
    fcalloc(1)=1
    counter = counter_old

    intermediate = 0.0_r8_kind

    ! first do s-type functions : this must be the order
    ! in which they do it in the old code ----------------------
    ! ATTENTION : in the variable orbitalprojection
    !             l=0  r2-type
    !             l=-1 s-type      (wuerg!)
    if (machineparameters_DimCheck) then
       if (counter.ne.orbitalprojection(-1,i_ua)) &
            call error_handler("contract_loc : something fishy (1)")
    endif

    uab => uabl(0)
    contract => uab%contractions
    pointer_int=>integral(0)%m
    counter = counter + uab%N_uncontracted_fcts
    do i_cont = 1,uab%N_contracted_fcts
       do i_exp=1,uab%N_exponents
          intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
               contract(i_exp,i_cont) * pointer_int(:,i_exp,1,:,:)
       enddo
       counter = counter + 1
    enddo

    counter = counter - uab%N_uncontracted_fcts - uab%N_contracted_fcts
    do i_uncont = 1,uab%N_uncontracted_fcts
       intermediate(:,:,:,counter) = pointer_int(:,i_uncont,1,:,:)
       counter = counter + 1
    enddo
    counter = counter + uab%N_contracted_fcts
    ! ----------- end of s-type --------------------------------

    ! now do the r2-type (as soon as debugging is over we can
    ! drop the awkward order again) ----------------------------
    if (machineparameters_DimCheck) then
       if (counter.ne.orbitalprojection(0,i_ua)) &
            call error_handler &
            ("contract_loc : something fishy (2)")
    endif

    uab => uabr
    contract => uab%contractions
    pointer_int=>integral(-1)%m
    counter = counter + uab%N_uncontracted_fcts
    do i_cont = 1,uab%N_contracted_fcts
       do i_exp=1,uab%N_exponents
          intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
               contract(i_exp,i_cont) * pointer_int(:,i_exp,1,:,:)
       enddo

       counter = counter + 1
    enddo

    counter = counter - uab%N_uncontracted_fcts - uab%N_contracted_fcts
    do i_uncont = 1,uab%N_uncontracted_fcts
       intermediate(:,:,:,counter) = pointer_int(:,i_uncont,1,:,:)
       counter = counter + 1
    enddo
    counter = counter + uab%N_contracted_fcts
    ! ---------- end of r2-type --------------------------------

    ! loop over angular momenta
    angular_mom: do i_l = 1,lmax
       pointer_int=>integral(i_l)%m
       n_independent = &
            ua%symadapt_partner(i_ir,i_l)%n_independent_fcts

       if (machineparameters_DimCheck) then
          if (counter.ne.orbitalprojection(i_l,i_ua)) &
               call error_handler &
               ("contract_loc : something fishy (3)")
       endif

       uab => uabl(i_l)
       counter_ind = counter

       ! loop over independent fcts
       do i_ind = 1,n_independent
          counter = counter_ind + i_ind - 1

          ! local contract before renormalisation and renormalise result
          contract => uab%contractions
          counter = counter + uab%N_uncontracted_fcts * n_independent
          do i_cont = 1,uab%N_contracted_fcts
             do i_exp=1,uab%N_exponents
                intermediate(:,:,:,counter) = intermediate(:,:,:,counter) + &
                     contract(i_exp,i_cont) * pointer_int(:,i_exp,i_ind,:,:)
             enddo

             counter = counter + n_independent
          enddo

          counter = counter - n_independent * &
               (uab%N_uncontracted_fcts + uab%N_contracted_fcts)
          ! pick those which remain uncontracted: these are
          ! assumed to be the FIRST N_uncontracted functions for
          ! each primitive fitbasis
          do i_uncont = 1,uab%N_uncontracted_fcts
             intermediate(:,:,:,counter) =  &
                  pointer_int(:,i_uncont,i_ind,:,:)
             counter = counter + n_independent
          enddo
          counter = counter + uab%N_contracted_fcts * n_independent
       enddo

       if (n_independent.ne.0) counter = counter - n_independent + 1

    enddo angular_mom

    if(present(gradients)) then ! gradient case
       ! do summation over fitfunctions
       m2_loop: do m2 = 1, n_m2
          m1_loop: do m1 = 1, n_m1
             ! ADD Sum(k) a_k,tot [ ... f_k ... ]  (old fashion technique)
             i_ab = 1
             do i_b = 1, size(gradients,2)
                do i_a  = 1, size(gradients,1)
                   if ( cutoff(i_a,i_b) ) then
                      gradients(i_a,i_b,m1,m2) = gradients(i_a,i_b,m1,m2) + &
                           sum( coeff_charge_eperef(counter_old:counter-1) * &
                           intermediate(i_ab,m1,m2,counter_old:counter-1) )
                      i_ab = i_ab + 1
                   endif
                enddo
             enddo

          enddo m1_loop
       enddo m2_loop
    endif

    deallocate(intermediate,STAT=fcalloc(1))
    ASSERT(fcalloc(1).eq.0)
  end subroutine contract_loc_eperef
#endif

  subroutine print_fcalloc()
  integer (kind=i4_kind):: i
   do i=1,size(fcalloc)
    if(fcalloc(i).ne.0) print*, i, 'fcalloc ne 0'
   enddo
  end subroutine print_fcalloc

end module fitcontract_module

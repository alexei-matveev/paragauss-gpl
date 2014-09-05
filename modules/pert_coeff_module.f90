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
  !===================================================================
! Public interface of module
  !===================================================================
module pert_coeff_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: Calculate coefficients of the perturbation theory.
  !
  !  Module called by: chargefit.f90 decide.f90
  !
  !  Public Data: pert_vec   containing [rho_spin|f_k]
  !               pert_vec_xc containing <rho_spin|g_l> (for EXT-MDA)
  !               pert_vec2  see below
  !               pert_mat   containing [phi_i,spin,irrep^2|f_k]
  !
  !  References: ...
  !
  !  Author: TG
  !  Date: 10/95
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   2/7/97
  ! Description: Introduction of the OPTIONS_PERTURBATION_THEORY switch
  !              Unused arrays like PERT_MAT_HELP removed
  !              Module Description Updated
  !              Subroutine UPDATE_CHARGE_OVERLAP_MATRIX from module
  !              LINSYS_CHARGEFIT_MODULE included here
  !              PERT_VEC made spin polarized (for the MDA approach)
  !
  ! Modification
  ! Author: MM
  ! Date:   10/97
  ! Description: Extension to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: Calculation of EXT-MDA projection <rho_spin|g_l> introduced
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  use datatype       ! contains user defined data types
  use options_module, only: options_xcmode, &
                            xcmode_model_density, xcmode_extended_mda, &
                            options_spin_orbit
  use operations_module, only: operations_core_density

  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module ====================


  !------------ Declaration of variables ------------
  real(kind=r8_kind), allocatable, public  :: pert_mat(:,:,:,:) ! (nk, ...) but (nff, *) on master
  ! pert_mat(i_ff,i_irrep,i_spin,i_pairs)
  ! with i_ff running over the part of fitfunctions stored on local host
  !
  ! in case of spin orbit pert_mat is complex
  real(kind=r8_kind), allocatable, public  :: pert_mat_real(:,:,:) ! (nk, n_irrep, n_pairs) but (nff, *) on master
  real(kind=r8_kind), allocatable, public  :: pert_mat_imag(:,:,:) ! (nk, n_irrep, n_pairs) but (nff, *) on master

  ! pert_mat(i_ff,i_irrep,i_spin,i_pairs)
  real(kind=r8_kind), allocatable, public  :: pert_vec(:,:) ! (nk, n_proj) but (nff, *) on master
  ! pert_vec(i_ff,i_p)
  ! with i_p from 1 to n_proj (if xcmode_model_density n_spin, else 1)
  ! with i_ff running over the part of fitfunctions stored on local host
  real(kind=r8_kind), allocatable, public  :: pert_vec_core(:) ! (nk) but (nff) on master
  real(kind=r8_kind), allocatable, public  :: pert_vec_xc(:,:) ! (i_ff,i_p)
  real(kind=r8_kind), allocatable, public  :: pert_vec2(:, :, :) ! (nk, n_irrep, n_spin) but (nff, *) on master
  !
  ! pert_vec2(:, irrep, spin) corresponds to the fit expansion of the
  ! density "variation" which would result if one removed an electron
  ! from an occupied level "e2" and put it into virtual level "e1".
  ! The corresponding density matrix is symmetric/hermitean and is given
  ! by the rank-2 expression:
  !
  !                H         H
  !     dP =  v * v  -  v * v
  !            2   2     1   1
  !
  ! Here H = T in regular (no spin-orbit) case.
  !
  ! In fact, the total density variation is a sum of such terms for
  ! each irrep because in general PT eventually selects an occupied/virtual
  ! pair for several irreps.
  !
  ! Thus this must be wrong:
  !
  !     "in case of spin orbit pert_vec2 is complex"
  !

  ! pert_vec2(i_ff,i_irrep,i_spin)
  ! with i_ff running over the part of fitfunctions stored on local host
  ! stores sum_ij (c_ih*c_jh-c_il*c_ij)*[ij|k]
  ! where h and l are the HOMO and LUMO
  ! for every irrep and every spin


  !------------ public functions and subroutines ---------------------
  public :: pcalc_pert_coeff
  public :: pert_coeff_free
  public :: update_charge_overlap_matrix
  public calc_cpks3c_ai


  !===================================================================
  ! End of public interface of module
  !===================================================================

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine calc_pert_coeff(nk, nj, perturbation_theory)
    !
    !  Purpose:  Computes  the righthand-side  for  the  density  fit and  the
    !           coefficients of  which the  corrections to density  matrix are
    !           assembled later
    !
    !  This subroutine is executed in parallel context by all workers.
    !
    !  Changed by subroutine:
    !  pert_vec         righthand-side for density fit
    !  pert_mat            coefficients
    !  pert_vec_xc      exchange fit projection for the extended MDA
    !
    !  Subroutine called by: pcorr_calc
    !
    !  Author: T. Grauschopf
    !  Date: 7/95
    !
    use comm, only: comm_rank
    use symmetry_data_module , only: symmetry_data_n_spin
    use bounds_module, only: bounds_ch
    use density_data_module, ONLY: densmat,densmat_real,densmat_imag,&
                                   core_densmat
    use occupation_module, ONLY: n_occo
    use occupied_levels_module, ONLY: eigvec_occ,eigvec_occ_real,eigvec_occ_imag
    use virtual_levels_module, ONLY: eigvec_vir, eigvec_vir_real, eigvec_vir_imag
    use integralstore_module, only: integralstore_3c_co
    use options_module, only: options_integrals_on_file
    implicit none
    !------------ Declaration of formal parameters -----------------
    integer(kind=i4_kind), intent(in) :: nk ! number of my (charge) fit functions
    integer(kind=i4_kind), intent(in) :: nj ! number of my (xc) fit functions
    logical, intent(in)               :: perturbation_theory
    !** End of interface *****************************************

    !
    !  * A few words of introduction to  a concept of "pairs" of orbitals used
    !    in this file. Much of  the wording below is the "reverse engeneering"
    !    not done by the original author. Beware of the misinterpretation.
    !
    !  n_pairs
    !
    !    Number of pairs of the  perturbation theory. A parameter, fixed to 10
    !    at the  moment (see  pairs_module for details).   However not  all of
    !    "pairs" are meaningful.  Also those  that are meaningfull are not the
    !    leading ones but the *trailing* ones.
    !
    !  n_pairs0(spin, irrep)
    !
    !    Least  index,  whose  associated  pair  satifies  a  certain  quality
    !    criterion (associated weight greater  than some threshold).  Only the
    !    trailing pairs  in the range starting from  n_pairs0(spin, irrep) and
    !    ending with  n_pairs (==10)  are meaningfull. I  am not yet  making a
    !    statment if the list of valid pairs is necessarily non-empty.
    !
    !  ptpair1(irrep)%m(pair, spin),
    !  ptpair2(irrep)%m(pair, spin), n_pairs0(spin, irrep) <= pair <= n_pairs
    !
    !    First and  second index  of an orbital  pair within  the perturbation
    !    theory. For a valid pair it is in the range [1, dim(irrep)]. It seems
    !    to be assumed that  if a pair that has an index  0 then it is invalid
    !    ("if, but  not only  if", see n_pairs0  above).  If you  consider the
    !    array axis indexed by "pair" only  zero or more (FIXME: one or more?)
    !    trailng pairs may be meanigful.
    !
    !  eigvec_occ(irrep)%m(:, :, spin)
    !
    !    Orbital coefficients of  the occupied orbitals. As some  of the pairs
    !    involve virtual orbitals this array is not always sufficient to build
    !    the matrix  elements of the  perturbation theory. The second  axis of
    !    the %m array is the orbital index.
    !
    !  eigvec_vir(irrep)%m(:, :, spin)
    !
    !    Orbital coefficients of the first  virtual orbitals. Not all of them,
    !    just  enough to get  PT theory  working. Therefore  it should  not be
    !    directly subscripted with the indices from ptpair1/2.
    !
    !  densmat            Density matrix
    !
    !--------------------------------------------------------------
    !  Modifications:
    !
    !  1. The double sum over orbitals has been changed considerably.  the sum
    !     now runs ONLY over the lower triangular matrix.  The reason for this
    !     change is, that  the 3-center integrals only exist  for unique pairs
    !     of  orbitals (because  the integrals  are symmetric!).  Thus  if the
    !     indices run over the full  square matrix, the integrals are not read
    !     in correctly.   Unfortunately this means, that  one if-statement has
    !     to be included. I have looked this up from the OLD lcgto - obviously
    !     its authors also didn`t find  a workaround for the if-statement.  If
    !     anybody has a better idea, please let me know.  FN 29.8.95
    !

    !  2. The perturbation_theory flag is  now given as an argument instead of
    !     the   old  variante   where  it   has  been   read  in   trough  the
    !     options_perturbation_theory function the main  reason for it is that
    !     now during a run the perturbation  theory flag can be changed, as by
    !     the  diis method so  now its  important that  all the  functions and
    !     processes always have the most recent variant of it AN 7.2009
    !
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind)  :: irrep, n_spin, n_proj
    logical :: integrals_on_file, extended_mda

    integer(i4_kind) :: rank ! base-0
    integer(i4_kind), allocatable :: dimensions(:) ! (n_irrep), irrep dimensions
    !
    ! Note that  the 3-center  integrals may be  viewed here as  a rectangular
    ! matrix of the shape
    !
    !   (nk, nij)
    !
    ! To  build the  RHS for  the density  fitting procedure  we  form "matrix
    ! product"  of these  fit coefficients  with the  elements of  the density
    ! matrix viewed as a vector of length nij.  Dont forget the factor two for
    ! the off-diagonal elements!
    !

    DPRINT 'calc_pert_coeff: entered'

    ! verify consistency of redundant data:
    rank = comm_rank()
    ASSERT(nk==bounds_ch(1+rank))

    n_spin = symmetry_data_n_spin()

    integrals_on_file = options_integrals_on_file()
    extended_mda = options_xcmode() == xcmode_extended_mda
    if (options_xcmode() == xcmode_model_density .or. extended_mda) then
      n_proj = n_spin
    else
      n_proj = 1
    endif

    pert_vec(:nk, :) = 0.0_r8_kind
    if (operations_core_density) pert_vec_core(:nk) = 0.0_r8_kind
    if (extended_mda) then
       pert_vec_xc=0.0_r8_kind
    endif

    pert_init: if (perturbation_theory) then
       pert_vec2(:nk, :, :) = 0.0_r8_kind

       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          pert_mat_real(:nk, :, :) = 0.0_r8_kind
          pert_mat_imag(:nk, :, :) = 0.0_r8_kind
       else ! options_spin_orbit
          !
          ! STANDARD SCF (NO SPIN ORBIT)
          !
          pert_mat(:nk, :, :, :) = 0.0_r8_kind
       endif ! options_spin_orbit
    endif pert_init

  if (options_spin_orbit) then
    !
    ! SPIN ORBIT
    !

    !
    ! Irrep dimensions:
    !
    allocate(dimensions(size(densmat_real)))
    do irrep = 1, size(densmat_real)
        dimensions(irrep) = size(densmat_real(irrep)%m, 1)
    enddo

    call spin_orbit(nk, dimensions) ! piece of code, see contains
  else ! options_spin_orbit
    !
    ! STANDARD SCF (NO SPIN ORBIT)
    !

    !
    ! Irrep dimensions:
    !
    allocate(dimensions(size(densmat)))
    do irrep = 1, size(densmat)
        dimensions(irrep) = size(densmat(irrep)%m, 1)
    enddo

    call standard(nk, dimensions) ! pice of code, see contains
  endif ! options_spin_orbit

    DPRINT 'calc_pert_coeff: exit'

contains

  subroutine standard(nk, dimensions)
    use istore, only: istore_t
    use istore, only: istore_new, istore_apply_and_advance
    implicit none
    integer(i4_kind), intent(in) :: nk
    integer(i4_kind), intent(in) :: dimensions(:) ! (n_irrep)
    ! *** end of interface ***

    integer(i4_kind) :: irrep
    integer(i4_kind) :: irhs, nrhs
    real(r8_kind), allocatable :: rhs(:, :) ! (n*(n+1)/2, nrhs)
    real(r8_kind), allocatable :: lhs(:, :) ! (nk, nrhs)
    integer(i4_kind) :: nij
    type(istore_t) :: coul

    ! FIXME: needs update:
    ASSERT(.not.extended_mda)

    !
    ! This is enough to hold all irreps of a symmetric matrix:
    !
    nij = sum( dimensions * (dimensions + 1) / 2 )

    !
    ! Integral store object:
    !
    coul = istore_new("coul", nk, nij, integralstore_3c_co)

    !
    ! Now loop over irreps:
    !
    do irrep = 1, size(dimensions)
        !
        ! Count number of rhs to be computed in one shot
        ! by convolution with 3-center integrals:
        !
        nrhs = 0
        call tasks(irrep, nrhs)

        ! Storage for rhs:
        allocate(rhs(dimensions(irrep)*(dimensions(irrep)+1)/2, nrhs))

        !
        ! Now fill the inputs for convolution procedure for all tasks:
        !
        irhs = 0
        call tasks(irrep, irhs, inputs=rhs)
        ASSERT(irhs==nrhs)

        !
        ! Sorage for convolution results of coulomb ints and rhs:
        !
        allocate(lhs(nk, nrhs))

        !
        ! Convolution of the 3c integrals and the density matrix,
        ! sum over subrange of integrals corresponding to one irrep
        ! here:
        !
        call istore_apply_and_advance(coul, rhs, lhs)

        !
        ! Now let the tasks get the outputs and do what they were
        ! intended for:
        !
        irhs = 0
        call tasks(irrep, irhs, outputs=lhs)
        ASSERT(irhs==nrhs)

        ! next irrep will have different dimensions:
        deallocate(rhs, lhs)
    enddo
  end subroutine standard

  subroutine spin_orbit(nk, dimensions)
    use istore, only: istore_t
    use istore, only: istore_new, istore_apply_and_advance
    implicit none
    integer(i4_kind), intent(in) :: nk
    integer(i4_kind), intent(in) :: dimensions(:) ! (n_irrep)
    ! *** end of interface ***

    integer(i4_kind) :: irrep
    integer(i4_kind) :: irhs, nrhs
    real(r8_kind), allocatable :: rhs(:,:,:) ! (2, n*(n+1)/2, nrhs)
    real(r8_kind), allocatable :: lhs(:, :) ! (nk, nrhs)
    integer(i4_kind) :: nij
    type(istore_t) :: coul

    if ( .not. integrals_on_file ) then
        WARN("untested")
    endif

    !
    ! This is enough to hold all irreps of a symmetric matrix:
    !
    nij = sum( dimensions * (dimensions + 1) / 2 )

    !
    ! Integral store object (complex integrals, double size):
    !
    coul = istore_new("coul", nk, 2*nij, integralstore_3c_co)

    !
    ! Now loop over irreps:
    !
    do irrep = 1, size(dimensions)

        !
        ! Count number of rhs to be computed in one shot:
        !
        nrhs = 0
        call spin_orbit_tasks(irrep, nrhs)

        !
        ! Storage for rhs, real and image:
        !
        allocate(rhs(2, dimensions(irrep)*(dimensions(irrep)+1)/2, nrhs))

        !
        ! Fill the inputs to be convoluted with coulomb integrals:
        !
        irhs = 0
        call spin_orbit_tasks(irrep, irhs, inputs=rhs)
        ASSERT(irhs==nrhs)

        ! storage for lhs:
        allocate(lhs(nk, nrhs))

        !
        ! Convolute 3c integrals and all of rhs, advance
        ! the storage pointer to the next irrep:
        !
        call istore_apply_and_advance(coul, rhs, lhs)

        !
        ! Provide the outputs to the code that knows what
        ! to do with them:
        !
        irhs = 0
        call spin_orbit_tasks(irrep, irhs, outputs=lhs)
        ASSERT(irhs==nrhs)

        ! next irrep will have different dimensions:
        deallocate(rhs, lhs)
    enddo
  end subroutine spin_orbit

  subroutine tasks(irrep, itask, inputs, outputs)
    !
    ! Logic of what needs to be convoluted with the coulomb integrals
    ! depending on the global settings and irrep indices is not quite
    ! trivial. This sub is responisble for this logic.
    !
    ! a) If only integer itask is present than it is responsible to increment
    ! it for each new subtask for this irrep
    !
    ! b) If inputs are present, this sub is responsible for filling the input array
    ! at the position specified by the running itask.
    !
    ! c) If outputs are present, take them from the respective positions and
    ! do whatever needs to be done with it.
    !
    use pairs_module, only: ptpair1, ptpair2, n_pairs0, n_pairs
    implicit none
    integer(i4_kind), intent(in) :: irrep
    integer(i4_kind), intent(inout) :: itask
    real(r8_kind), optional, intent(out) :: inputs(:, :) ! (n*(n+1)/2, nrhs)
    real(r8_kind), optional, intent(in) :: outputs(:, :) ! (nk, nrhs)
    ! *** end of interface ***

    integer(i4_kind) :: spin, pair
    integer(i4_kind) :: n1, n2

    !
    ! itask couns number of rhs to be computed in one shot
    ! by convolution with 3-center integrals for this particular
    ! irrep.
    !
    ASSERT(n_proj==1)

    !
    ! Convolution with density matrix is unconditional,
    ! we also need contributions of all irreps:
    !
    ! Sum of all irrep contributions goes to
    !
    !   pert_vec(:nk, 1)
    !
    itask = itask + 1
    if ( present(inputs) ) then
        inputs(:, itask) = tripack(sum(densmat(irrep)%m, 3))
    endif
    if ( present(outputs) ) then
        pert_vec(:nk, 1) = pert_vec(:nk, 1) + outputs(:, itask)
    endif

    !
    ! Handling of core_densmat to calculate pert_vec_core
    ! is quite similar:
    !
    ! Sum of all irrep contributions goes to
    !
    !   pert_vec_core(:nk)
    !
    if (operations_core_density) then
        itask = itask + 1
        if ( present(inputs) ) then
            inputs(:, itask) = tripack(sum(core_densmat(irrep)%m, 3))
        endif
        if ( present(outputs) ) then
            pert_vec_core(:nk) = pert_vec_core(:nk) + outputs(:, itask)
        endif
    endif

    !
    ! Perturbation theory gives rize to zero or more
    ! irrep blocks, the number differes from irrep to irrep:
    !
    if (perturbation_theory) then
        !
        ! If any irrep has a valid ptpair1/ptpair2, add one rhs,
        !
        ! Irrep/spin contributions go to
        !
        !       pert_vec2(:nk, irrep, spin)
        !
        do spin = 1, n_spin
           !
           ! Here one is apparently interested  in a single pair from the list
           ! with the largest weight. The  list of pairs correspond to a range
           ! [n_pairs0(spin, irrep), n_pairs], so do not confuse n_pairs0 with
           ! the number of pairs.:
           !
           if ( n_pairs0(spin, irrep) > n_pairs ) cycle ! to next spin

           !
           ! The last pair indexed by n_pairs with the largest weigth is in
           ! the list! Also note the order of axes in n_pairs0(:, :).
           !
           itask = itask + 1

           ! These two  indices, together with irrep and  spin identify a
           ! pair of orbitals:
           n1 = ptpair1(irrep)%m(n_pairs, spin)
           n2 = ptpair2(irrep)%m(n_pairs, spin)

           if ( present(inputs) ) then
              call dmat_pert(-1, irrep, spin, n1, n2, inputs(:, itask))
           endif
           if ( present(outputs) ) then
              pert_vec2(:nk, irrep, spin) = outputs(:, itask)
           endif
        enddo

        !
        ! For pert_mat one rhs per spin per pair (whatever that is):
        !
        ! Irrep/spin/pair contributions go to
        !
        !       pert_mat(:nk, irrep, spin, pair)
        !
        do spin = 1, n_spin
            do pair = n_pairs0(spin, irrep), n_pairs
                itask = itask + 1

                ! These two  indices, together with irrep and  spin identify a
                ! pair of orbitals:
                n1 = ptpair1(irrep)%m(pair, spin)
                n2 = ptpair2(irrep)%m(pair, spin)

                if ( present(inputs) ) then
                    call dmat_pert(+1, irrep, spin, n1, n2, inputs(:, itask))
                endif
                if ( present(outputs) ) then
                    pert_mat(:nk, irrep, spin, pair) = outputs(:, itask)
                endif
            enddo
        enddo
    endif
  end subroutine tasks

  subroutine spin_orbit_tasks(irrep, itask, inputs, outputs)
    !
    ! See tasks ...
    !
    use pairs_module, only: ptpair1, ptpair2, n_pairs0, n_pairs
    implicit none
    integer(i4_kind), intent(in) :: irrep
    integer(i4_kind), intent(inout) :: itask
    real(r8_kind), optional, intent(out) :: inputs(:, :, :) ! (2, n*(n+1)/2, nrhs)
    real(r8_kind), optional, intent(in) :: outputs(:, :) ! (nk, nrhs)
    ! *** end of interface ***

    integer(i4_kind) :: pair, n1, n2

    !
    ! itask couns number of rhs to be computed in one shot
    ! by convolution with 3-center integrals for this particular
    ! irrep.
    !

    !
    ! Convolution with density matrix is unconditional,
    ! we also need contributions of all irreps:
    !
    ! Sum of all irrep contributions goes to
    !
    !   pert_vec(:nk, 1)
    !
    itask = itask + 1

    if ( present(inputs) ) then
        !
        ! FIXME: istore given the expansion coefficients (x, y)
        !        over real and imaginary parts of the upper (lower?)
        !        triangle of the coulomb matrix
        !
        !               (coul_int_real, coul_int_imag)
        !
        !        computes the "dot product" not a "complex product"
        !
        !               result = x * coul_int_real + y * coul_int_imag
        !
        !        Apparently to get what we want here we need to flip
        !        the sign (or pack the lower triangle):
        !
        inputs(1, :, itask) =   tripack(densmat_real(irrep)%m)
        inputs(2, :, itask) = - tripack(densmat_imag(irrep)%m)
    endif

    if ( present(outputs) ) then
        pert_vec(:nk, 1) = pert_vec(:nk, 1) + outputs(:, itask)
    endif

    !
    ! Perturbation theory gives rize to zero or more
    ! irrep blocks, the number differes from irrep to irrep:
    !
    if (perturbation_theory) then
        !
        ! If any irrep has a valid ptpair1/ptpair2, add one rhs,
        !
        ! Irrep/spin contributions go to
        !
        !       pert_vec2(:nk, irrep, spin)
        !
        if ( n_pairs0(1, irrep) <= n_pairs ) then

            itask = itask + 1

            n1 = ptpair1(irrep)%m(n_pairs, 1)
            n2 = ptpair2(irrep)%m(n_pairs, 1)

            if ( present(inputs) ) then
                if (n2 > n_occo(1,irrep)) then
                    n2 = n2 - n_occo(1, irrep)

                    ! (homo-lumo) pair density diff:
                    call rank2h(eigvec_occ_real(irrep)%m(:, n1), eigvec_occ_imag(irrep)%m(:, n1), &
                                eigvec_vir_real(irrep)%m(:, n2), eigvec_vir_imag(irrep)%m(:, n2), &
                                inputs(1:2, :, itask))
                else
                    call rank2h(eigvec_occ_real(irrep)%m(:, n1), eigvec_occ_imag(irrep)%m(:, n1), &
                                eigvec_occ_real(irrep)%m(:, n2), eigvec_occ_imag(irrep)%m(:, n2), &
                                inputs(1:2, :, itask))
                endif
            endif

            if ( present(outputs) ) then
                pert_vec2(:nk, irrep, 1) = outputs(:, itask)
            endif
        endif

        !
        ! For pert_mat one rhs per spin per pair (whatever that is):
        !
        ! Irrep/pair contributions go to
        !
        !       pert_mat_real/imag(:nk, irrep, pair)
        !
        do pair = n_pairs0(1, irrep), n_pairs
            !
            ! Need real an imaginary parts of perturbation,
            ! carefull when indexing into inputs/outputs:
            !
            itask = itask + 2

            if ( present(inputs) ) then
                n1 = ptpair1(irrep)%m(pair, 1)
                n2 = ptpair2(irrep)%m(pair, 1)

                if ( n2 > n_occo(1, irrep) ) then
                    n2 = n2 - n_occo(1, irrep)

                    call pmat(eigvec_occ_real(irrep)%m(:, n1), eigvec_occ_imag(irrep)%m(:, n1), &
                              eigvec_vir_real(irrep)%m(:, n2), eigvec_vir_imag(irrep)%m(:, n2), &
                              inputs(1:2, :, itask-1:itask))
                else
                    call pmat(eigvec_occ_real(irrep)%m(:, n1), eigvec_occ_imag(irrep)%m(:, n1), &
                              eigvec_occ_real(irrep)%m(:, n2), eigvec_occ_imag(irrep)%m(:, n2), &
                              inputs(1:2, :, itask-1:itask))
                endif
            endif

            if ( present(outputs) ) then
                pert_mat_real(:nk, irrep, pair) = outputs(:, itask-1)
                pert_mat_imag(:nk, irrep, pair) = outputs(:, itask-0)
            endif
        enddo
    endif
  end subroutine spin_orbit_tasks

  function tripack(m) result(t)
    !
    ! Pack the upper triangle of symmetric,
    ! weight off-diagonal elements by factor two.
    !
    implicit none
    real(r8_kind), intent(in) :: m(:, :) ! (n, n)
    real(r8_kind) :: t(size(m,1)*(size(m,1)+1)/2) ! (n*(n+1)/2)
    ! *** end of interface ***

    integer(i4_kind) :: n, a, b, ab

    n = size(m, 1)
    ASSERT(n==size(m, 2))
    ASSERT(n*(n+1)/2==size(t))

    ab = 0
    do b = 1, n
        do a = 1, b - 1
            ab = ab + 1
            ! NOTE: factor two here:
            t(ab) = 2 * m(a, b)
        enddo
        ab = ab + 1
        ! NOTE: no factor two here:
        t(ab) = m(b, b)
    enddo
    ASSERT(ab==n*(n+1)/2)
  end function tripack

    subroutine dmat_pert(sgn, irrep, spin, n1, n2, dm)
      implicit none
      integer(i4_kind), intent(in) :: sgn
      integer(i4_kind), intent(in) :: irrep, spin, n1, n2
      real(r8_kind), intent(out)   :: dm(:) ! (n*(n+1)/2)
      ! *** end of interface ***

      integer :: a, b
      logical :: virtual

      ASSERT(n1>0)
      ASSERT(n2>0)

      !
      ! To  eventually access the  orbital coefficients  of a  virtual orbital
      ! from  the incomplete eigvec_vir(irrep)%m(:,  start:end, spin)  one may
      ! need to re-base the index:
      !
      virtual = (n2 > n_occo(spin, irrep))
      if ( virtual ) then
         a = n1
         b = n2 - n_occo(spin, irrep)
      else
         a = n1
         b = n2
      endif

      ASSERT(a>0)
      ASSERT(b>0)

      select case(sgn)
      case (+1)
         if( virtual )then
ASSERT(a<=size(eigvec_occ(irrep)%m, 2))
ASSERT(b<=size(eigvec_vir(irrep)%m, 2))
            call rank1a( eigvec_occ(irrep)%m(:, a, spin), &
                         eigvec_vir(irrep)%m(:, b, spin), &
                         dm )
         else
ASSERT(a<=size(eigvec_occ(irrep)%m, 2))
ASSERT(b<=size(eigvec_occ(irrep)%m, 2))
            call rank1a( eigvec_occ(irrep)%m(:, a, spin), &
                         eigvec_occ(irrep)%m(:, b, spin), &
                         dm )
         endif
      case (-1)
         if( virtual )then
ASSERT(a<=size(eigvec_occ(irrep)%m, 2))
ASSERT(b<=size(eigvec_vir(irrep)%m, 2))
            call rank2s( eigvec_occ(irrep)%m(:, a, spin), &
                         eigvec_vir(irrep)%m(:, b, spin), &
                        dm )
         else
ASSERT(a<=size(eigvec_occ(irrep)%m, 2))
ASSERT(b<=size(eigvec_occ(irrep)%m, 2))
            call rank2s( eigvec_occ(irrep)%m(:, a, spin), &
                         eigvec_occ(irrep)%m(:, b, spin), &
                         dm )
         endif
      case default
         dm = huge(dm)
         ABORT('no such case')
      end select
    end subroutine dmat_pert

    subroutine rank1a(e1, e2, f)
      !
      ! This is a non-symmetric rank-1 "perturbation"
      !
      !                   T
      !         f = e  * e
      !              1    2
      !
      ! For the purpose of forming a trace over (a, b)
      ! indices with symmetric coulomb integrals
      !
      !         [ab|k]
      !
      ! we return the sum of the off-diagonal positions
      !
      !         f(ab) := f(a, b) + f(b, a),  a /= b
      !
      ! Of course no such doubling for diagonal f(a, a).
      !
      implicit none
      real(r8_kind), intent(in)  :: e1(:), e2(:) ! (n)
      real(r8_kind), intent(out) :: f(:) ! (n*(n+1)/2)
      ! *** end of interface ***

      integer(i4_kind) :: n, a, b, ab

      n = size(e1)
      ASSERT(n==size(e2))
      ASSERT(n*(n+1)/2==size(f))

      ab = 0
      do b = 1, n
         do a = 1, b - 1
            ab = ab + 1
            f(ab) = e1(a) * e2(b) + e1(b) * e2(a)
         enddo
            ab = ab + 1
            f(ab) = e1(b) * e2(b)
      enddo
    end subroutine rank1a

    subroutine rank2s(e1, e2, f)
      !
      ! This is a symmetric rank-2 contribution:
      !
      !                   T           T
      !         f = e  * e   -   e * e
      !              1    1       2   2
      !
      ! For the purpose of forming a trace over (a, b)
      ! indices with symmetric coulomb integrals
      !
      !         [ab|k] = [ba|k]
      !
      ! we return the sum of the off-diagonal positions
      !
      !         f(ab) := f(a, b) + f(b, a),  a /= b
      !
      ! Of course no such doubling for diagonal f(a, a).
      !
      ! For this symmetric expression the sum amounts
      ! to a factor two for the off-diagonal positions.
      !
      implicit none
      real(r8_kind), intent(in)  :: e1(:), e2(:) ! (n)
      real(r8_kind), intent(out) :: f(:) ! (n*(n+1)/2)
      ! *** end of interface ***

      integer(i4_kind) :: n, a, b, ab

      n = size(e1)
      ASSERT(n==size(e2))
      ASSERT(n*(n+1)/2==size(f))

      ab = 0
      do b = 1, n
         do a = 1, b - 1
            ab = ab + 1
            f(ab) = 2 * (e1(a) * e1(b) - e2(a) * e2(b))
         enddo
            ab = ab + 1
            f(ab) = e1(b) * e1(b) - e2(b) * e2(b)
      enddo
    end subroutine rank2s

    subroutine rank2h(r1, i1, r2, i2, f)
      !
      ! This must be the symmetric rank-2 contribution:
      !
      !                   H           H
      !         f = v  * v   -   v * v
      !              1    1       2   2
      !
      ! For the purpose of forming a trace over (a, b)
      ! indices with hermitean coulomb integrals
      !
      !         [ab|k] = [ba|k]*
      !
      ! we return the sum of the off-diagonal position
      ! and the conjugate of the corresponding
      ! transposed position
      !
      !         f(ab) := f(a, b) + f(b, a)*,  a /= b
      !
      ! Of course no such doubling for diagonal f(a, a).
      !
      ! For this hermitean expression the sum amounts
      ! to a factor two for the off-diagonal positions.
      !
      implicit none
      real(r8_kind), intent(in)  :: r1(:), i1(:) ! (n), real/imag eigenvec v1
      real(r8_kind), intent(in)  :: r2(:), i2(:) ! (n), real/imag eigenvec v2
      real(r8_kind), intent(out) :: f(:, :) ! (2, n*(n+1)/2), real/imag matrix
      ! *** end of interface ***

      integer(i4_kind) :: n, a, b, ab

      n = size(r1)
      ASSERT(n==size(i1))
      ASSERT(n==size(r2))
      ASSERT(n==size(i2))
      ASSERT(2==size(f,1))
      ASSERT(n*(n+1)==size(f))

      ab = 0
      do b = 1, n
         do a = 1, b
            ab = ab + 1
            f(1, ab) = r1(a) * r1(b) - r2(a) * r2(b) + i1(a) * i1(b) - i2(a) * i2(b)
            f(2, ab) = i1(a) * r1(b) - r1(a) * i1(b) - i2(a) * r2(b) + r2(a) * i2(b)

            if ( a /= b ) then
                f(1, ab) = f(1, ab) * 2
                f(2, ab) = f(2, ab) * 2
            endif
         enddo
      enddo
    end subroutine rank2h

    subroutine pmat(r1, i1, r2, i2, p)
      !
      ! This must be the non-symmetric rank-1 "perturbation":
      !
      !                   H
      !         f = v  * v
      !              1    2
      !
      ! NOTE: the matrix is not symmetric, for a trace with an
      !       upper triangle of a hermitean matrix we need to sum
      !       off-diagonal elements
      !
      ! With c(1:2) being the real and imaginary parts of the coul ints
      ! the result p(1:2, 1:2) contains the coefficients to build the
      ! complex off-diagonal matrix element of the perturbation:
      !
      !       < v | f | v >
      !          2       1
      !
      ! with real and imaginary part given by linear combination:
      !
      !       f(i) = sum  c(k) * p(k, i)
      !                 k
      !
      ! here i = 1, 2 and k = 1, 2 corresponding to real and imaginary parts
      ! and the ab-axis is omitted.
      !
      implicit none
      real(r8_kind), intent(in)  :: r1(:), i1(:) ! (n), real/imag eigenvec v1
      real(r8_kind), intent(in)  :: r2(:), i2(:) ! (n), real/imag eigenvec v2
      real(r8_kind), intent(out) :: p(:,:,:) ! (2, n*(n+1)/2, 2)
      ! *** end of interface ***

      integer(i4_kind) :: n, a, b, ba
      real(r8_kind) :: x, y

      n = size(r1)
      ASSERT(n==size(i1))
      ASSERT(n==size(r2))
      ASSERT(n==size(i2))
      ASSERT(n*(n+1)/2==size(p,2))
      ASSERT(2==size(p,1))
      ASSERT(2==size(p,3))

      ba = 0
      do a = 1, n
          do b = 1, a
              ba = ba + 1

              x = r1(b) * r2(a) + i1(b) * i2(a)
              y = r1(b) * i2(a) - i1(b) * r2(a)

              p(1, ba, 1) = + x
              p(2, ba, 1) = - y
              p(1, ba, 2) = + y
              p(2, ba, 2) = + x

              if ( a /= b ) then
                  x = r1(a) * r2(b) + i1(a) * i2(b)
                  y = r1(a) * i2(b) - i1(a) * r2(b)

                  p(1, ba, 1) = p(1, ba, 1) + x
                  p(2, ba, 1) = p(2, ba, 1) + y
                  p(1, ba, 2) = p(1, ba, 2) + y
                  p(2, ba, 2) = p(2, ba, 2) - x
              endif
          enddo
      enddo
    end subroutine pmat

  end subroutine calc_pert_coeff

  subroutine calc_cpks3c_ai(ilower,iupper,th_mo_coul)
    !
    !   Purpose:
    !   calculates G_ai_kl*S_kl of Qai using  cpks_gvec intermediate
    !
    !  Changed by subroutine:
    !  pert_vec         righthand-side for density fit
    !  pert_mat            coefficients
    !  pert_vec_xc      exchange fit projection for the extended MDA
    !
    !  Subroutine called by: pcorr_calc
    !
    !
    use comm_module
    use msgtag_module
    use symmetry_data_module
    use iounitadmin_module ! provides IO-units
    use readwriteblocked_module    ! contains the routines for buffered I/O
    use filename_module      ! set I/O-Filenames
    use prepare_integralfiles_module
    use integralstore_module, only: integralstore_3c_co
    use options_module, only: options_integrals_on_file
    use cpksdervs_matrices
    use fit_coeff_module, only: get_fit, fit
    use eigen_data_module,only: eigvec
    implicit none
    !------------ Declaration of formal parameters -----------------
    !  ilower/iupper : Range of charge fit functions to work on
    !  jlower/jupper : Range of exchange fit functions to work on
    integer(kind=i4_kind),           intent(in) :: ilower,iupper
    type(readwriteblocked_tapehandle),intent(inout) :: th_mo_coul
    !** End of interface *****************************************


    !  Modifications:
    !  1. The double sum over orbitals has been changed considerably.
    !     the sum now runs ONLY over the lower triangular matrix.
    !     The reason for this change is, that the 3-center integrals
    !     only exist for unique pairs of orbitals (because the integrals
    !     are symmetric!). Thus if the indices run over the full square
    !     matrix, the integrals are not read in correctly.
    !     Unfortunately this means, that one if-statement has to be
    !     included. I have looked this up from the OLD lcgto - obviously
    !     its authors also didn`t find a workaround for the if-statement.
    !     If anybody has a better idea, please let me know.
    !     FN 29.8.95
    !

    !------------ Declaration of local variables --------------------
    type(fit)                   :: n_fit

    integer(kind=i4_kind)  :: i_i, i_s, i_a, i_b, n1, n2, &
         i_meta, i_last, i_step, n_spin, n_irrep
    integer(kind=i4_kind)  :: k_step, k_dim, klower, kupper, n_split

    integer(kind=i4_kind)  :: k, i_occ ! required to make dummy allocation
    integer(kind=i4_kind)  :: eig_dim,occ_dim,vir_dim

    real(kind=r8_kind),pointer :: coul_int(:), eigvec1(:)
    real(kind=r8_kind), allocatable:: co_i(:,:,:)
    real(kind=r8_kind), allocatable:: co_ai(:)
    type(readwriteblocked_tapehandle) :: th_ch
    logical :: integrals_on_file

    real(kind=r8_kind) :: summ
    real(kind=r8_kind) :: fact
    !------------ Executable code ------------------------------

    n_split=n_cpks_split
    n_spin =symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()
    integrals_on_file = options_integrals_on_file()

    call get_fit(n_fit)

    i_step = iupper - ilower + 1

#define no_co_ai
#ifdef no_co_ai
 call readwriteblocked_startwrite(trim(tmpfile('mo_coul.dat')), th_mo_coul)
#endif

    !
    ! These iterations over fit function indices define format
    ! of the on-disk file mo_coul.dat and therefore have to be
    ! done identically at the place where it is read, see
    ! cpks_g4constructs.f90, calc_cpksQai.f90.
    !
    ! k_step == 1 means processing one fit-function at a time:
    k_step = i_step / n_split + 1
    kupper = 0
    do while ( kupper < i_step )
        ! process in pieces the range of k from klower to kupper:
        klower = kupper + 1
        kupper = kupper + k_step

        ! last piece may be smaller:
        kupper = min(kupper, i_step)

        ASSERT(kupper-klower>=0)
!       print *, "WRITE: klower=", klower, "kupper=", kupper

        k_dim  = kupper - klower + 1 ! <= k_step

    spl_spins: do i_s=1,n_spin

    if (integrals_on_file) then
       allocate(coul_int(i_step),STAT=cpksalloc(6))
       ASSERT(cpksalloc(6).eq.0)
       call prepare_integral_open('coul', th_ch) ! incl.  start read
    else
       i_meta = 1
    endif

    do i_i=1,n_irrep

    eig_dim=size(eigvec(i_i)%m,1)
    occ_dim=cpks3c(i_i,i_s)%occ_dim
    vir_dim=cpks3c(i_i,i_s)%vir_dim

    allocate(co_i(k_dim, occ_dim, eig_dim), co_ai(vir_dim), stat=cpksalloc(161))
    ASSERT(cpksalloc(161).eq.0)
    co_i=0.0_r8_kind
       do i_a=1,symmetry_data_dimension(i_i)
          do i_b=1,i_a
             if (i_a.eq.i_b) then
                fact = 0.5_r8_kind
             else
               fact = 1.0_r8_kind
             endif

             if (integrals_on_file) then
                call readwriteblocked_read(coul_int,th_ch)
             else
                i_last = i_meta + i_step - 1
                coul_int => integralstore_3c_co(i_meta:i_last)
                i_meta = i_last + 1
             endif


                do n1=1,occ_dim
                 eigvec1 => eigvec(i_i)%m(:,n1,i_s)
                  co_i(:,n1,i_b)=co_i(:,n1,i_b)+fact*coul_int(klower:kupper)*eigvec1(i_a)
                  co_i(:,n1,i_a)=co_i(:,n1,i_a)+fact*coul_int(klower:kupper)*eigvec1(i_b)
                enddo

          enddo! i_b loop
       enddo!i_a loop

    do k=1,size(co_i,1)
      do i_occ=1,occ_dim
        do n2=1,vir_dim
          summ = 0.0_r8_kind
          do i_a=1,size(co_i,3)
            summ = summ + co_i(k, i_occ, i_a) * eigvec(i_i)%m(i_a, n2 + eig_dim - vir_dim, i_s)
          enddo
          co_ai(n2) = summ
        enddo
        call readwriteblocked_write(co_ai(:),th_mo_coul)
      enddo
    enddo

       deallocate(co_i,stat=cpksalloc(161))
       ASSERT(cpksalloc(161).eq.0)
       deallocate(co_ai,stat=cpksalloc(161))
       ASSERT(cpksalloc(161).eq.0)
       cpksalloc(161)=1
    enddo! i_ir

    if (integrals_on_file) then
       call prepare_integral_close(th_ch)
       deallocate(coul_int,STAT=cpksalloc(6))
       ASSERT(cpksalloc(6).eq.0)
       cpksalloc(6)=1
    else
      nullify(coul_int)
    endif
  enddo spl_spins
    enddo ! while

#ifdef no_co_ai
  call readwriteblocked_stopwrite(th_mo_coul)
#endif

  end subroutine calc_cpks3c_ai

  subroutine pcalc_pert_coeff()
    !
    ! Purpose: wrapper of assemblation of right hand side and
    !          perturbation corrections
    !
    ! To be executed in parallel context, called by cahrgefit()
    !
    !  Author: T. Grauschopf
    !  Date: 8/95
    !
    !
    !--------------------------------------------------------------
    ! Modifications:
    ! 1. The receiving has been deleted, since the slave-branch
    !    in 'mainscf' does the receiving whereas the slave only
    !    unpacks the most recently received buffer
    !
    ! 2. The case ilower > iupper, i.e. the local part of the
    !    fitfunctions is empty, is implemented.
    !    For this purpose an initialization module 'init_module'
    !    has been included, which sets those variables to zero
    !    that are to be computed by 'corr_calc'. If ilower>iupper
    !    'corr_calc' is NOT called and the variables remain zero.
    ! 3. AS 3/99
    !    I deleted communpack(n_pairs,info)
    ! 4.  the perturbation_theory flag is now recieved from master
    !     instead of the old variante where it has been read in trough
    !     the options_perturbation_theory function
    !     the main reason for it is that now during a run the perturbation
    !     theory flag can be changed, as by the diis method
    !     so now its important that all the functions and processes always
    !     have the most recent variant of it
    !     This is especially important in this case because before the slaves had
    !     their own copy of the flag which was independent from the one of master
    !     AN 7.2009
    !
    use comm, only: comm_rank, comm_bcast
    use symmetry_data_module, only: symmetry_data_n_irreps, symmetry_data_n_proj_irreps, &
        symmetry_data_n_spin
    use options_module, only: options_xcmode, xcmode_model_density, &
        xcmode_extended_mda, options_spin_orbit, options_perturbation_theory
    use operations_module, only: operations_core_density
    use bounds_module, only: bounds_ch, bounds_xc
    use pairs_module, only: n_pairs ! a const
    use pairs_module, only: ptpairs, deallocate_pairs
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: pid ! rank + 1
    integer(i4_kind) :: nk, nj, n_proj, n_irrep, n_spin
    logical :: perturbation_theory, extended_mda

    pid = comm_rank() + 1

    !
    ! NOTE: perturbation_theory flag may be changed during the SCF run
    ! only the master has the most recent value
    ! FIXME: make options_perturbation_theory() return consistent results
    !
    perturbation_theory = options_perturbation_theory()

    ! master knows better, tells it to slaves:
    call comm_bcast(perturbation_theory)

    if ( options_spin_orbit ) then
        n_irrep = symmetry_data_n_proj_irreps()
    else
        n_irrep = symmetry_data_n_irreps()
    endif

    n_spin = symmetry_data_n_spin()

    extended_mda = options_xcmode() == xcmode_extended_mda
    if (options_xcmode() == xcmode_model_density .or. extended_mda) then
        n_proj = n_spin
    else
        n_proj = 1
    endif

    if ( perturbation_theory ) then
      !
      ! Calculation of the pairs of the perturbation theory
      !
      call ptpairs()
    endif

    if ( pid == 1 ) then
        !
        ! NOTE: on master, size(pert_vec, 1) and co is the total range of fit funcitons.
        ! Slaves, on the other hand, allocate, compute and deliver only the
        ! smaller range assigned to them:
        !
        nk = sum(bounds_ch(:))
        nj = sum(bounds_xc(:))
    else
        !
        ! Slaves allocate data only for their subrange:
        !
        nk = bounds_ch(pid)
        nj = bounds_xc(pid)
    endif
    call pert_coeff_alloc(nk, nj, n_proj, n_irrep, n_spin)

    !
    ! However the hard work is shared evenly, so that even the master
    ! has to process only a range of fit-funcitons.
    !
    ! Number of charge- and xc-fit functions assigned to this worker:
    !
    nk = bounds_ch(pid)
    nj = bounds_xc(pid)

    !
    ! Pass the dimensions unconditionally. The perturbation_theory will be
    ! given as an argument to make sure that the the function uses the most
    ! recent value (may be changed during the scf).
    !
    ! This computes the charge density projections [rho|f_k]
    ! and (eventually) the extended MDA exchange fit projection <rho|g_l>
    !
    call calc_pert_coeff(nk, nj, perturbation_theory)

    !
    ! Gather results on master:
    !
    call gather()

    !
    ! Cleanup is performed at the end of chargefit() which
    ! is our caller.
    !
    contains

      subroutine gather()
        !
        !  Purpose: gather (short) results of each worker into (long)
        !  results on master.
        !
        ! FIXME: maybe one can reduce the number of calls to GATHERV?
        !
        use comm, only: comm_gather
        implicit none
        ! *** end of interface ***

        integer(i4_kind) :: irrep, spin, pair

        !
        ! Receive results (right hand side and coefficients from the sons).
        !
        ! Note that master "gathers" resutls from slaves by putting
        ! each contribution into its own array range by using an in-place
        ! version of GATHERV.
        !

        do spin = 1, size(pert_vec, 2)
            call comm_gather(pert_vec(:, spin), bounds_ch)
        enddo

        if (operations_core_density) then
            !
            ! The core fit functions are distributed precisely as the charge
            ! fit functions:
            !
            call comm_gather(pert_vec_core, bounds_ch)
        endif

        if (perturbation_theory) then
            do spin = 1, size(pert_vec2, 3)
                do irrep = 1, size(pert_vec2, 2)
                    call comm_gather(pert_vec2(:, irrep, spin), bounds_ch)
                enddo
            enddo

            ! FIXME: not all (irrep, pair) entries contain meanigfull data:
            if (options_spin_orbit) then
                !
                ! SPIN ORBIT
                !
                do pair = 1, size(pert_mat_real, 3)
                    do irrep = 1, size(pert_mat_real, 2)
                        call comm_gather(pert_mat_real(:, irrep, pair), bounds_ch)
                        call comm_gather(pert_mat_imag(:, irrep, pair), bounds_ch)
                    enddo
                enddo
            else
                !
                ! STANDARD SCF (NO SPIN ORBIT)
                !
                do pair = 1, size(pert_mat, 4)
                    do spin = 1, size(pert_mat, 3)
                        do irrep = 1, size(pert_mat, 2)
                            call comm_gather(pert_mat(:, irrep, spin, pair), bounds_ch)
                        enddo
                    enddo
                enddo
            endif ! options_spin_orbit
        endif ! perturbation_theory

        if (extended_mda) then
            do spin = 1, size(pert_vec_xc, 2)
                call comm_gather(pert_vec_xc(:, spin), bounds_xc)
            enddo
        endif
      end subroutine gather

  subroutine pert_coeff_alloc(nk, nj, n_proj, n_irrep, n_spin)
    !
    ! Kept here as an internal sub because of logical
    ! flags affecting execution.
    !
    implicit none
    integer(i4_kind), intent(in) :: nk, nj, n_proj, n_irrep, n_spin
    ! *** end of interface ***

    integer(i4_kind) :: alloc_stat

    allocate(pert_vec(nk, n_proj), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    if ( operations_core_density ) then
        allocate(pert_vec_core(nk), STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    endif

    if ( extended_mda ) then
        allocate(pert_vec_xc(nj, n_proj), STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    endif

    if ( perturbation_theory ) then
        allocate(pert_vec2(nk, n_irrep, n_spin), STAT=alloc_stat)
        ASSERT(alloc_stat==0)

        if ( options_spin_orbit ) then
            !
            ! SPIN ORBIT
            !
            allocate(pert_mat_real(nk, n_irrep, n_pairs), &
                     pert_mat_imag(nk, n_irrep, n_pairs), STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        else
            !
            ! STANDARD SCF (NO SPIN ORBIT)
            !
            allocate(pert_mat(nk, n_irrep, n_spin, n_pairs), STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        endif
    endif
  end subroutine pert_coeff_alloc

  end subroutine pcalc_pert_coeff

  subroutine  update_charge_overlap_matrix(f_mat)
    !  Purpose: assembles the corrections to the matrix elements
    !           of the density fit from the coefficients computed
    !           by the slaves
    !
    !  This code is called from chargefit.f90 and runs in master only
    !  context
    !
    !  Data used:
    !
    !  Name:              Description/Range:
    !
    !  ssym               symmetry flags  -> datatypes
    !  n_pairs            Number of pairs in perturbation theory
    !  n_pairs0           Least indices of pairs whose weights
    !                     satisfy quality criterion
    !  ptweit             Weights from perturbation theory
    !  nbch               Number of fit functions
    !  pert_mat            Coefficients from the slaves of which
    !                     we will assemble the fit matrix
    !
    !  Modified data:
    !
    !  Name:              Description/Range:
    !
    !  f_mat              The matrix which will be updated
    !                     packed storage mode introduced (UB: 7/97)
    !                     ij=0; DO i=1,n_ch; DO j=1,i; ij=ij+1; ij=>(i,j) ; &
    !                     ENDDO; ENDDO
    !
    !  Subroutine called by: chargefit
    !
    !  Author: T.Grauschopf
    !  Date: 7/95
    !
    !--------------------------------------------------------------

    use symmetry_data_module
    use pairs_module, ONLY : n_pairs, n_pairs0, ptweit
    use fit_coeff_module, ONLY : get_fit,fit
    use iounitadmin_module ! provides IO-units
    implicit none
    !** End of interface *****************************************

    type(fit)                   :: n_fit

    ! --- Declaration of formal parameters ------------------------
    real(kind=r8_kind),intent(inout) :: f_mat(:)
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) i_i,i_s,i_p,i,j, ij
    real(kind=r8_kind) f

    call get_fit(n_fit)


    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i_i=1,ssym%n_proj_irrep
             do i_p=n_pairs0(1,i_i),n_pairs
                f=ptweit(i_i)%m(i_p,1)
                ij = 0
                do j=1,n_fit%n_ch
                   do i=1,j
                      ij = ij + 1 ! (j,i)
                      f_mat(ij)=f_mat(ij)+f*(pert_mat_real(i,i_i,i_p)&
                           & *pert_mat_real(j,i_i,i_p) + &
                           pert_mat_imag(i,i_i,i_p)*pert_mat_imag(j,i_i,i_p))
                   enddo
                enddo
             enddo
       enddo
    else ! options_spin_orbit
       do i_i=1,ssym%n_irrep
          do i_s=1,ssym%n_spin
             do i_p=n_pairs0(i_s,i_i),n_pairs
                f=ptweit(i_i)%m(i_p,i_s)
                ij = 0
                do j=1,n_fit%n_ch
                   do i=1,j
                      ij = ij + 1 ! (j,i)
                      f_mat(ij)=f_mat(ij)+f*pert_mat(i,i_i,i_s,i_p)&
                           & *pert_mat(j,i_i,i_s,i_p)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif! options_spin_orbit
  end subroutine update_charge_overlap_matrix

  subroutine pert_coeff_free()
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: alloc_stat

    deallocate(pert_vec, STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    if ( allocated(pert_vec_core) ) then
        deallocate(pert_vec_core, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    endif

    if ( allocated(pert_vec_xc) ) then
        deallocate(pert_vec_xc, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
    endif

    if ( allocated(pert_vec2) ) then
        deallocate(pert_vec2, STAT=alloc_stat)
        ASSERT(alloc_stat==0)

        if ( allocated(pert_mat_real) ) then
            !
            ! SPIN ORBIT
            !
            deallocate(pert_mat_real, pert_mat_imag, STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        else
            !
            ! STANDARD SCF (NO SPIN ORBIT)
            !
            deallocate(pert_mat, STAT=alloc_stat)
            ASSERT(alloc_stat==0)
        endif
    endif
  end subroutine pert_coeff_free

end module pert_coeff_module

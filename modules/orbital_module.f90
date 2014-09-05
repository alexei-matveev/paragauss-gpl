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
module orbital_module
  !-------------------------------------------------------------------
  !
  !  Calcualation  of contracted  and symmetry  adapted  orbital basis
  !  functions and charge and  exchange fitfunction basis functions on
  !  an array of gridpoints.
  !
  !  The  results  are given  in  arrays  of  type orbital_type.   The
  !  subroutines orbital_allocate,  orbital_free and orbital_write can
  !  be used to handle theses arrays.
  !
  !  This     module    uses    data     from    symmetry_data_module,
  !  unique_atom_module, orbitalprojection_module.   All these modules
  !  must be initialised before.
  !
  !  Before  calculations  can  be  performed, orbital_setup  must  be
  !  called to perform  preparatory calculations and memory allocation
  !  for intermediates. The memory can be freed with orbital_shutdown.
  !
  !  The  calculation   of  orbitals  is   performed  with  subroutine
  !  orbital_calculate.
  !
  !
  !  Module called by: ...
  !
  !
  !  References: Publisher-document: orbital-calculation
  !
  !
  !  Author: TB
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
  ! Author: TB
  ! Date:   29.4.96
  ! Description: implementation of contractions changed
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   10/96
  ! Description: ordering of indices in grads variable changed
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date: 2/97
  ! Description: second derivatives of primitive orbitals added
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:  3/97
  ! Description: Additional subroutines, which calculate derivatives
  !              with respect to nuclear displacements and not with
  !              respect to the electronic coordinate
  !
  ! Modification
  ! Author: HH
  ! Date:   11/97
  ! Description: Merged changes introduced by MS in V11alpha version
  !              modified for use with the response program (view
  !              "ttfs_hh_V11alpha"):
  !              - new datatype "orbital_spin_type"
  !              - optional parameters "phi_ob" of this type in
  !                subroutines orbital_allocate(), orbital_free()
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date: 12/97
  ! Description: Derivatives of fitting functions introduced
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:  1/98
  ! Description: Addaption to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author:
  ! Date:
  ! Description:
  !
  !-------------------------------------------------------------------

  !------------ Modules used -----------------------------------------
#ifdef FPP_DEBUG
! define FPP_TIMERS 2
#endif
# include "def.h"
  use type_module ! type specification parameters
  use solid_harmonics_module, only: solid_harmonics_type, &
    solid_harmonics_grads_type, solid_harmonics_sec_der_type, &
    solid_harmonics_sec_der_type ! calculates solid harmonics
  use orbitalstore_module, only: orbital_type, orbital_gradient_type, &
    spinor_gradient_type, orbital_sec_der_type, orbital_nuclear_gradient_type, &
    orbital_nuclear_sec_der_type, orbital_spin_type, spinor_type
    ! declaration of grid-orbital-types and RW procedures
  USE_MEMLOG

  implicit none
  save ! save all variables defined in this module
  private         ! by default, all names are private

  !== Interrupt end of public interface of module ====================

  ! The various types of derivatives of symmetry adapted orbitals xi_i
  ! which are evaluated in this module are
  !                 xi_i(r) = Sum(n) c_ni           f_n(r-R_n)
  !           d/dr  xi_i(r) = Sum(n) c_ni      d/dr f_n(r-R_n)
  !         [-d/dR] xi_i(r) = Sum(n) c_ni      d/dr f_n(r-R_n) delta_R,R_n
  !   d/dr    d/dr  xi_i(r) = Sum(n) c_ni d/dr d/dr f_n(r-R_n)
  ! [-d/dR]   d/dr  xi_i(r) = Sum(n) c_ni d/dr d/dr f_n(r-R_n) delta_R,R_n
  ! [-d/dR] [-d/dS] xi_i(r) = delta_R,S * [-d/dR] d/dr xi_i(r)
  ! with f_n(r-R_n) being any primitive Gaussian centered at an atom R_n

  !------------ Declaration of public types -----------------------
  !
  ! Declarations of these types have been moved to
  ! orbitalstore_module:
  !
  ! xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,lin(ua,i),prt)
  ! type, public :: orbital_type
  ! type, public :: spinor_type
  ! type, public :: orbital_spin_type ! for the response part

  ! d/dr_n xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,lin(ua,i),prt)
  ! type, public :: orbital_gradient_type

  ! d/dr_n xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,lin(ua,i),prt)
  ! type, public :: spinor_gradient_type

  ! d/dr_n d/dr_m xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,m,lin(ua,i),prt)
  ! type, public :: orbital_sec_der_type

  ! [-d/dR_n] xi[ua,irr,prt]_i(r_a) = <var>(irr,ua)%o(a,n,R,i,prt)
  ! type, public :: orbital_nuclear_gradient_type

  ! d/dR_n d/dR_m xi[ua,irr,prt]_i(r_a) = <var>(irr,ua)%o(a,lin(n,m),R,i,prt)
  ! type, public :: orbital_nuclear_sec_der_type

  ! reduced variants of the orbital types used for core density fitting fcts.
  ! xi[ua]_i(r_a) = <var>%o(a,lin(ua,i))
  ! type, public :: core_orbital_type

  ! d/dr_n xi[ua]_i(r_a) = <var>%o(a,n,lin(ua,i))
  ! type, public :: core_orbital_gradient_type

  ! d/dr_n d/dr_m xi[ua]_i(r_a) = <var>%o(a,n,m,lin(ua,i))
  ! type, public :: core_orbital_sec_der_type

  ! [-d/dR_n] xi[ua]_i(r_a) = <var>(ua)%o(a,n,R,i)
  ! type, public :: core_orbital_nuc_gradient_type

  ! d/dR_n d/dR_m xi[ua]_i(r_a) = <var>(ua)%o(a,lin(n,m),R,i)
  ! type, public :: core_orbital_nuc_sec_der_type

  !------------ public functions and subroutines ---------------------
  public orbital_allocate, orbital_free, &
       orbital_setup, orbital_shutdown, &
       orbital_calculate, orbital_write, &
       fit_fct_allocate, fit_fct_free, &
       fit_fct_calculate, &
       coulomb_pot_of_bas,print_orbmod_alloc

#ifdef WITH_CORE_DENS
  public :: fit_core_calculate
  public :: core_orbital_allocate
  public :: core_orbital_free
#endif

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of private types ----------------------

  ! The various types of derivatives of solid harmonics used to
  ! evaluate the orbitals and their derivatives on the grid are
  !               Y_lm(r_a-R_ua,ea) = <var>(ua)%sh (ea)%l(l)%m(a,m)
  !        d/dr_i Y_lm(r_a-R_ua,ea) = <var>(ua)%shg(ea)%l(l)%m(a,i,m)
  ! d/dr_i d/dr_j Y_lm(r_a-R_ua,ea) = <var>(ua)%shs(ea)%l(l)%m(a,i,j,m)
  ! Note, that x-Rx = Y_12(r-R), y-Ry = Y_13(r-R), and z-Rz = Y_10(r-R)

  type solid_harmonics_array_type
     ! solid harmonics of all equal atoms of one unique atome
     type(solid_harmonics_type), pointer :: sh(:)
     type(solid_harmonics_grads_type), pointer :: shg(:)
     type(solid_harmonics_sec_der_type), pointer :: shs(:)
     type(solid_harmonics_sec_der_type), pointer :: sht(:)
  end type solid_harmonics_array_type

  type glob_con_intermediate_l_type
     ! to store intermediate values of global contraction of one unique atom:
     real(r8_kind), pointer :: oo  (:,:,:) ! kept for compatibility
     ! oo  (grid_length,N_exponents,N_independent_fcts)
     real(r8_kind), pointer :: o   (:,:,:)
     ! o   (grid_length,N_exponents,N_independent_fcts)
     real(r8_kind), pointer :: og  (:,:,:,:)
     ! og  (grid_length,3,N_exponents,N_independent_fcts)
     real(r8_kind), pointer :: osd (:,:,:,:)
     ! osd (grid_length,6,N_exponents,N_independent_fcts)
     real(r8_kind), pointer :: ong (:,:,:,:,:)
     ! ong (grid_length,3,N_equal_atoms,N_exponents,N_independent_fcts)
     real(r8_kind), pointer :: onsd(:,:,:,:,:)
     ! onsd(grid_length,6,N_equal_atoms,N_exponents,N_independent_fcts)
     ! 1: xx  2: xy  3: xz  4: yy  5: yz  6: zz
  end type glob_con_intermediate_l_type

  type glob_con_intermediate_type
     ! to store intermediate values of global contraction of one unique atom:
     type(glob_con_intermediate_l_type), pointer :: l(:) ! l(-1:lmax) s,r2,p,...
  end type glob_con_intermediate_type

  !------------ Declaration of private constants and variables ----
  real(r8_kind), parameter, private :: pi = 3.14159265358979324_r8_kind
  integer,                                 private :: vec_len = 0
  ! vector length: number of grid points processed together
  integer,                                 private :: N_max_exponents
  ! maximal number of exponents of primitive basis
  integer,                                 private :: N_max_equal_atoms
  ! maximal number of equal atoms for all unique atoms
  type(solid_harmonics_array_type), allocatable, target, private :: &
       solid_harmonics(:) ! solid_harmonics(N_unique_atoms)
  !
  ! FIXME: these  arrays are allocatad  once and used as  temp working
  !        arrays.  Why  again  are  we  not  using  local  subroutine
  !        varibales?
  !
  ! Intermediate for subroutine calculate_orbital, these are allocated
  ! in orbital_setup() for a  specific "vector_length" and system size
  ! (see   N_max_exponents,   N_max_equal_atoms   above).   They   are
  ! deallocated in orbital_shutdown().
  !
  real(r8_kind), allocatable, private :: &
       prim_orb_fac(:,:,:),  cont_orb_fac(:,:,:), &
       prim_grad_fac(:,:,:), cont_grad_fac(:,:,:), &
       prim_sd_fac(:,:,:),   cont_sd_fac(:,:,:)
  !
  ! These   working  arrays   above  are   global  because   they  are
  ! pre-allocated to the maximum size  that will be encountered in the
  ! loop over  atoms and the  re-used at every integral  batch. FIXME:
  ! make them local vars of respective subs.  Until then, they need to
  ! be threadprivate:
  !
  !$omp threadprivate(                 &
  !$omp& prim_orb_fac,  cont_orb_fac,  &
  !$omp& prim_grad_fac, cont_grad_fac, &
  !$omp& prim_sd_fac,   cont_sd_fac    )
  !
  real(r8_kind), allocatable, private :: prim_td_fac(:,:,:), &
       cont_td_fac(:,:,:)
  ! prim_xxx_fac(vec_len,N_max_exponents,N_max_equal_atom)
  ! cont_xxx_fac(vec_len,N_max_exponents,N_max_equal_atom)
  !
  ! Factors  of  primitive  orbitals   to  be  multiplied  with  solid
  ! hormonics intermediate for subroutine orbital_calculate.

  type(glob_con_intermediate_type), pointer :: uncont_orbs_ch(:), &
                                               uncont_orbs_xc(:)
  ! Uncontracted primitive basis functions needed as intermediates for
  ! subroutine  calculate_orbital  if global  contractions  are to  be
  ! done.

  logical :: do_glob_con_ch, do_glob_con_xc
  ! true if  global contractions  for charge or  exchange fitfunctions
  ! are to be done

  logical :: memory_allocated= .false., &
       grads_allocated = .false., sec_der_allocated = .false.
  logical :: third_der_allocated = .false.

  integer(i4_kind), parameter ::&
       & RE = 1, &
       & IM = 2, &
       & UP = 1, &
       & DN = 2
  integer(i4_kind), parameter :: &
       & VALUE = 0, &
       & DX    = 1, &
       & DY    = 2, &
       & DZ    = 3

  integer(i4_kind), parameter ::&
       & XX = 1, &
       & XY = 2, &
       & XZ = 3, &
       & YY = 4, &
       & YZ = 5, &
       & ZZ = 6

  real(r8_kind), parameter ::&
       & ZERO = 0.0_r8_kind,&
       & ONE  = 1.0_r8_kind,&
       & TWO  = 2.0_r8_kind

 integer(i4_kind) :: status(63) = 1



  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine print_orbmod_alloc()
    integer(kind=i4_kind):: i
    do i = 1, size(status)
       if (status(i) .eq. 0) print*, i, ' orbmod_alloc '
    enddo
    DPRINT ' orbital module deallocations checked'
  end subroutine print_orbmod_alloc

  subroutine orbital_allocate (orbs_ob, grads, sec_der, nuc_grads, &
       nuc_sec_der, orb_xc, phi_ob, phi_gr_nuc, num_gr_nuc, &
       orbs_spinor_ob, spinor_grads, nuc_3rd_der)
    !
    ! Allocate   orbs(N_irreps),  grads(N_irreps),  sec_der(N_irreps),
    ! nuc_grads, nuc_sec_der, orb_xc, phi_orb
    !
    use datatype, only: arrmat5
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc, fit_coeff_n_cd
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalstore_module, only: orbitalstore_initialized, orbitalstore_allocate
    use symmetry_data_module ! description of irreps
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(orbital_type), pointer,          optional :: orbs_ob(:)
    type(orbital_gradient_type), pointer, optional :: grads(:)
    type(spinor_gradient_type), pointer, optional :: spinor_grads(:)
    type(orbital_sec_der_type), pointer,  optional :: sec_der(:)
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_sec_der(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_3rd_der(:,:)
    type(orbital_type), intent(inout), optional :: orb_xc
    ! argument  kept  for  compatibility.  To  allocate  exchange  fit
    ! functions   on  the   grid  please   use  fit_fct_allocate("xc",
    ! fcts=orb_xc, ...).  But see xcmda_hamiltonian.f90
    type(orbital_spin_type),pointer,      optional :: phi_ob(:)
    type(arrmat5)          ,pointer,      optional :: phi_gr_nuc(:)
    integer(i4_kind)       ,intent(in),   optional :: num_gr_nuc
    type(spinor_type), pointer,           optional :: orbs_spinor_ob(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_ir, n_orbitals, i_unique, i_l
    integer(i4_kind) :: bas_dim,n_partn !,ic
    integer(i4_kind) :: n_v_irr, n_ea, n_spin, n_ff
    !------------ Executable code ------------------------------------

    if (present (orbs_ob)) then
       if (orbitalstore_initialized) then
          call orbitalstore_allocate (orbs_ob=orbs_ob)
       else
          n_v_irr = symmetry_data_n_irreps()
          allocate (orbs_ob(n_v_irr), stat=status(58))
          ASSERT(status(58).eq.0)
          do i_ir = 1, n_v_irr
             bas_dim = symmetry_data_dimension(i_ir)
             n_partn = symmetry_data_n_partners(i_ir)
             allocate (orbs_ob(i_ir)%o(vec_len, bas_dim, n_partn), stat=status(1))
             ASSERT(status(1).eq.0)
             MEMLOG(size(orbs_ob(i_ir)%o))
          enddo
       end if
    endif

    if ( present(orbs_spinor_ob) ) then
       !
       ! SPIN ORBIT
       !
       if( orbitalstore_initialized )then
          call orbitalstore_allocate(orbs_spinor_ob=orbs_spinor_ob)
       else
          allocate( orbs_spinor_ob(symmetry_data_n_proj_irreps()),stat=status(2) )
           ASSERT(status(2).eq.0)
          do i_ir = 1, symmetry_data_n_proj_irreps()
             bas_dim = symmetry_data_dimension_proj(i_ir)
             n_partn = symmetry_data_n_partners_proj(i_ir)

             allocate(&
                  & orbs_spinor_ob(i_ir)%o(vec_len,bas_dim,n_partn,2,2),&
                  & STAT=status(2))
             if(status(2).ne.0) call error_handler&
                  & ("allocate_orbital: alloc orbs_spinor_ob(i_ir)%o failed")
             MEMLOG(vec_len*bas_dim*n_partn*2*2)
!            do ic=UP,DN
!               orbs_spinor_ob(i_ir)%spinor(ic)%o      => orbs_spinor_ob(i_ir)%o(:,:,:,RE,ic)
!               orbs_spinor_ob(i_ir)%spinor(ic)%o_imag => orbs_spinor_ob(i_ir)%o(:,:,:,IM,ic)
!            enddo
          enddo
       endif
    endif ! present(orbs_spinor_ob)

    if ( present(spinor_grads) ) then
       if( orbitalstore_initialized )then
          call orbitalstore_allocate(spinor_grads=spinor_grads)
       else
          allocate( spinor_grads(symmetry_data_n_proj_irreps()), stat=status(3) )
          if ( status(3) .ne. 0 ) call error_handler( &
               "orbital_allocate: allocate of spinor_grads failed")
          do i_ir = 1, symmetry_data_n_proj_irreps()
             bas_dim = symmetry_data_dimension_proj(i_ir)
             n_partn = symmetry_data_n_partners_proj(i_ir)

             allocate( spinor_grads(i_ir)%o(vec_len, 3, bas_dim, n_partn, RE:IM, UP:DN), &
                  stat=status(3) )
             ASSERT(status(3)==0)
             MEMLOG(vec_len*3*bas_dim*n_partn*2*2)
!            do ic=UP,DN
!               spinor_grads(i_ir)%spinor(ic)%o      => spinor_grads(i_ir)%o(:,:,:,:,RE,ic)
!               spinor_grads(i_ir)%spinor(ic)%o_imag => spinor_grads(i_ir)%o(:,:,:,:,IM,ic)
!            enddo
          enddo
       endif ! else of orbitalstore_initialized
    endif ! spinor_grads

    if ( present(grads) ) then
       if( orbitalstore_initialized )then
          call orbitalstore_allocate(grads=grads)
       else
          allocate(grads(symmetry_data_n_irreps()), stat=status(59) )
          ASSERT(status(59).eq.0)

          do i_ir = 1, symmetry_data_n_irreps()
             bas_dim = symmetry_data_dimension(i_ir)
             n_partn = symmetry_data_n_partners(i_ir)
             allocate( grads(i_ir)%o( vec_len, 3, bas_dim, n_partn ), &
                  stat=status(4) )
             ASSERT(status(4).eq.0)
             MEMLOG(size(grads(i_ir)%o))
          enddo
       end if
    endif

    if ( present(sec_der) ) then
       allocate( sec_der(symmetry_data_n_irreps()), stat=status(60) )
       ASSERT(status(60).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          bas_dim = symmetry_data_dimension(i_ir)
          n_partn = symmetry_data_n_partners(i_ir)
          allocate( sec_der(i_ir)%o( vec_len, 6, bas_dim, n_partn ), &
                    stat=status(5) )
          ASSERT(status(5).eq.0)
          MEMLOG(size(sec_der(i_ir)%o))
       enddo
    endif

    if ( present(nuc_grads) ) then

       allocate( nuc_grads(n_unique_atoms,symmetry_data_n_irreps()), &
            stat=status(61) )
       ASSERT(status(61).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          do i_unique=1,n_unique_atoms
             n_ea = unique_atoms(i_unique)%n_equal_atoms
             n_partn = symmetry_data_n_partners(i_ir)

             n_orbitals=0
             do i_l=0,unique_atoms(i_unique)%lmax_ob
                n_orbitals=n_orbitals+unique_atoms(i_unique)%&
                     symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                     (unique_atoms(i_unique)%l_ob(i_l)%n_contracted_fcts+&
                     unique_atoms(i_unique)%l_ob(i_l)%n_uncontracted_fcts)
             end do

             allocate( nuc_grads(i_unique,i_ir)%o( vec_len, 3 , n_ea, n_orbitals, n_partn ), &
                  stat=status(6) )
             ASSERT(status(6).eq.0)
             MEMLOG(size(nuc_grads(i_unique,i_ir)%o))
!!$#ifdef FPP_DEBUG
!!$             n_orbitals = n_orbitals - sum(BasDim(i_unique)%LM(:,i_ir))
!!$             ASSERT(n_orbitals==0)
!!$#endif
          enddo
       end do
    endif

    if ( present(nuc_sec_der) ) then
       allocate( nuc_sec_der(n_unique_atoms,symmetry_data_n_irreps()), &
            stat=status(62) )
       ASSERT(status(62).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          do i_unique=1,n_unique_atoms
             n_ea = unique_atoms(i_unique)%n_equal_atoms
             n_partn = symmetry_data_n_partners(i_ir)

             n_orbitals=0
             do i_l=0,unique_atoms(i_unique)%lmax_ob
                n_orbitals=n_orbitals+unique_atoms(i_unique)%&
                     symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                     (unique_atoms(i_unique)%l_ob(i_l)%n_contracted_fcts+&
                     unique_atoms(i_unique)%l_ob(i_l)%n_uncontracted_fcts)
             end do

             allocate( nuc_sec_der(i_unique,i_ir)%o( vec_len, 6 , n_ea, n_orbitals, n_partn ), &
                  stat=status(7) )
             ASSERT(status(7).eq.0)
             MEMLOG(size(nuc_sec_der(i_unique,i_ir)%o))
          enddo
       end do
    endif

    if ( present(nuc_3rd_der) ) then
       allocate( nuc_3rd_der(n_unique_atoms,symmetry_data_n_irreps()), &
            stat=status(51) )
       ASSERT(status(51).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          do i_unique=1,n_unique_atoms
             n_ea = unique_atoms(i_unique)%n_equal_atoms
             n_partn = symmetry_data_n_partners(i_ir)

             n_orbitals=0
             do i_l=0,unique_atoms(i_unique)%lmax_ob
                n_orbitals=n_orbitals+unique_atoms(i_unique)%&
                     symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                     (unique_atoms(i_unique)%l_ob(i_l)%n_contracted_fcts+&
                     unique_atoms(i_unique)%l_ob(i_l)%n_uncontracted_fcts)
             enddo

             allocate( nuc_3rd_der(i_unique,i_ir)%o( vec_len, 10 , n_ea, n_orbitals, n_partn ), &
                  stat=status(52) )
             ASSERT(status(52).eq.0)
             MEMLOG(size(nuc_3rd_der(i_unique,i_ir)%o))
          enddo
       end do
    endif

    if ( present(orb_xc) ) then
       i_ir = get_totalsymmetric_irrep()
       n_ff = fit_coeff_n_xc()
       allocate( orb_xc%o( vec_len, n_ff, 1 ), &
            stat=status(9) )
       if ( status(9) .ne. 0 ) call error_handler( &
            "orbital_allocate: allocate of orbs%o failed")
       status(9)=1
       MEMLOG(vec_len*n_ff)
       ! allocation of uncont_orbs_xc moved to orbital_allocate
    endif

    if ( present(phi_ob) ) then
       allocate( phi_ob(symmetry_data_n_irreps()), stat=status(63) )
       ASSERT(status(63).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          bas_dim = symmetry_data_dimension(i_ir)
          n_partn = symmetry_data_n_partners(i_ir)
          n_spin  = symmetry_data_n_spin()

          allocate( phi_ob(i_ir)&
               & %o( vec_len, bas_dim, n_partn , n_spin ),&
               stat=status(10) )
          ASSERT(status(10).eq.0)
          MEMLOG(vec_len*bas_dim*n_partn*n_spin)
       enddo
    endif

    if ( present(phi_gr_nuc) ) then
       ASSERT(present(num_gr_nuc))
       allocate( phi_gr_nuc(symmetry_data_n_irreps()), stat=status(63) )
       ASSERT(status(63).eq.0)

       do i_ir = 1, symmetry_data_n_irreps()
          bas_dim = symmetry_data_dimension(i_ir)
          n_partn = symmetry_data_n_partners(i_ir)
          n_spin  = symmetry_data_n_spin()
          ! num_gr_nuc = gradient_data_n_gradients ! from gradient_data_module

          allocate( phi_gr_nuc(i_ir)%m( vec_len, bas_dim, n_partn , num_gr_nuc, n_spin ),&
               stat=status(10) )
          ASSERT(status(10).eq.0)
          MEMLOG(vec_len*bas_dim*n_partn*n_spin*num_gr_nuc)
       enddo
    endif
  end subroutine orbital_allocate


#ifdef WITH_CORE_DENS
  subroutine core_orbital_allocate(orbs_ob_core,grads_core,&
             sec_der_core,nuc_grads_core,nuc_sec_der_core)
    !  Purpose: allocate orbs_core, grads_core
    !           sec_der_core, nuc_grads_core,nuc_sec_der_core
    use orbitalstore_module, only: core_orbital_type, core_orbital_gradient_type, &
        core_orbital_sec_der_type, core_orbital_nuc_gradient_type,  &
        core_orbital_nuc_sec_der_type
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc, fit_coeff_n_cd
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(core_orbital_type),   optional       :: orbs_ob_core
    type(core_orbital_gradient_type),optional :: grads_core
    type(core_orbital_sec_der_type), optional :: sec_der_core
    type(core_orbital_nuc_gradient_type), pointer, optional :: &
                                                 nuc_grads_core(:)
    type(core_orbital_nuc_sec_der_type), pointer, optional :: &
                                                 nuc_sec_der_core(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer              :: n_orbitals, i_unique
    !------------ Executable code ------------------------------------
    if ( present(orbs_ob_core) ) then
          allocate( orbs_ob_core%o( vec_len, &
               fit_coeff_n_cd()), &
               stat=status(11) )
          if ( status(11) .ne. 0 ) call error_handler( &
               "core_orbital_allocate: allocate of orbs%o failed")
          status(11)=1
    endif

    if ( present(grads_core) ) then
          allocate( grads_core%o( vec_len,3, &
               fit_coeff_n_cd()), &
               stat=status(12) )
          if ( status(12) .ne. 0 ) call error_handler( &
               "core_orbital_allocate: allocate of grads%o failed")
          status(12)=1
    endif

    if ( present(sec_der_core) ) then
          allocate( sec_der_core%o( vec_len,6, &
               fit_coeff_n_cd()), &
               stat=status(13) )
          if ( status(13) .ne. 0 ) call error_handler( &
               "orbital_allocate: allocate of sec_der%o failed")
          status(13)=1
    endif

    if ( present(nuc_grads_core) ) then
          allocate( nuc_grads_core(n_unique_atoms),stat=status(14) )
          if ( status(14) .ne. 0 ) call error_handler( &
          "core_orbital_allocate: allocate of nuc_grads_core failed")
          do i_unique=1,n_unique_atoms
                n_orbitals=unique_atoms(i_unique)%s_core&
                     %n_exponents + unique_atoms(i_unique)%r2_core&
                     %n_exponents
             allocate( nuc_grads_core(i_unique)%o( vec_len, 3 , &
                  unique_atoms(i_unique)%n_equal_atoms,&
                  n_orbitals) , &
                  stat=status(14) )
             if ( status(14) .ne. 0 ) call error_handler( &
                  "core_orbital_allocate: allocate of nuc_grads%o failed")
          enddo
          status(14)=1
    endif

    if ( present(nuc_sec_der_core) ) then
          allocate( nuc_sec_der_core(n_unique_atoms),stat=status(15) )
          if ( status(15) .ne. 0 ) call error_handler( &
          "core_orbital_allocate: allocate of nuc_sec_der_core failed")
          do i_unique=1,n_unique_atoms
                n_orbitals=unique_atoms(i_unique)%s_core&
                     %n_exponents + unique_atoms(i_unique)%r2_core&
                     %n_exponents
             allocate( nuc_sec_der_core(i_unique)%o( vec_len, 6 , &
                  unique_atoms(i_unique)%n_equal_atoms,&
                  n_orbitals), &
                  stat=status(15) )
             if ( status(15) .ne. 0 ) call error_handler( &
                  "core_orbital_allocate: allocate of nuc_sec_der%o failed")
          enddo
          status(15)=1
    endif

  end subroutine core_orbital_allocate
#endif

  subroutine orbital_free (orbs_ob, grads, sec_der, nuc_grads, &
       nuc_sec_der, orb_xc, phi_ob, phi_gr_nuc, orbs_spinor_ob, &
       spinor_grads, nuc_3rd_der)
    !
    ! Deallocate orbs(N_irreps), and other.
    !
    use datatype, only: arrmat5
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalstore_module, only: orbitalstore_initialized
    implicit none
    !------------ Declaration of formal parameters -------------
    type(orbital_type), pointer,          optional :: orbs_ob(:)
    type(spinor_type), pointer,          optional :: orbs_spinor_ob(:)
    type(orbital_gradient_type), pointer, optional :: grads(:)
    type(spinor_gradient_type), pointer, optional :: spinor_grads(:)
    type(orbital_sec_der_type),  pointer, optional :: sec_der(:)
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_sec_der(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_3rd_der(:,:)
    type(orbital_type), intent(inout),    optional :: orb_xc
    ! argument  kept for  compatibility.  To  deallocate  exchange fit
    ! functions on the grid please use fit_fct_free("xc", fcts=orb_xc,
    ! ...).  But see xcmda_hamiltonian.f90.
    type(orbital_spin_type),pointer,      optional :: phi_ob(:)
    type(arrmat5)          ,pointer,      optional :: phi_gr_nuc(:)
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer              :: i_ir, i_unique
    !------------ Executable code ------------------------------

    if ( present(orbs_ob) ) then
       if(.NOT.orbitalstore_initialized) then
          do i_ir = 1, size(orbs_ob,1)
             MEMLOG(-size(orbs_ob(i_ir)%o))
             deallocate( orbs_ob(i_ir)%o, stat=status(1) )
             ASSERT(status(1).eq.0)
             status(1)=1
          enddo

          deallocate( orbs_ob, stat=status(58) )
          ASSERT(status(58).eq.0)
          status(58)=1
       endif
    end if

    if ( present(orbs_spinor_ob) ) then
       !
       ! SPIN ORBIT
       !
       if(.NOT.orbitalstore_initialized) then
          DPRINT 'of: free orbs_spinor_ob...'
          do i_ir = 1, ubound(orbs_spinor_ob,1)
             MEMLOG(-size(orbs_spinor_ob(i_ir)%o))
             deallocate( orbs_spinor_ob(i_ir)%o,STAT=status(2))
             if ( status(2) .ne. 0 ) call error_handler( &
                  "free_orbital: deallocate of orbs_spinor%o failed")
          enddo
          deallocate( orbs_spinor_ob, stat=status(2) )
          if ( status(2) .ne. 0 ) call error_handler( &
               "free_orbital: deallocate of orbs_spinor failed")
       endif
    endif ! present(orbs_spinor_ob)

    if ( present(phi_ob) ) then
       do i_ir = 1, ubound(orbs_ob,1)
          MEMLOG(-size(phi_ob(i_ir)%o))
          deallocate( phi_ob(i_ir)%o, stat=status(10) )
          ASSERT(status(10).eq.0)
          status(10)=1
       enddo
       deallocate( phi_ob, stat=status(63) )
       ASSERT(status(63).eq.0)
       status(63)=1
    endif

    if ( present(phi_gr_nuc) ) then
       do i_ir = 1, ubound(phi_gr_nuc,1)
          MEMLOG(-size(phi_gr_nuc(i_ir)%m))
          deallocate( phi_gr_nuc(i_ir)%m, stat=status(10) )
          ASSERT(status(10).eq.0)
          status(10)=1
       enddo
       deallocate( phi_gr_nuc, stat=status(63) )
       ASSERT(status(63).eq.0)
       status(63)=1
    endif

    if ( present(grads) ) then
       if(.NOT.orbitalstore_initialized) then
          do i_ir = 1, size(grads,1)
             MEMLOG(-size(grads(i_ir)%o))
             deallocate( grads(i_ir)%o, stat=status(4) )
             ASSERT(status(4).eq.0)
             status(4)=1
          enddo
          deallocate( grads, stat=status(59) )
          ASSERT(status(59).eq.0)
          status(59)=1
       end if
    endif

    if ( present(spinor_grads) ) then
       if(.NOT.orbitalstore_initialized) then
          do i_ir = 1, size(spinor_grads,1)
             MEMLOG(-size(spinor_grads(i_ir)%o))
             deallocate(spinor_grads(i_ir)%o, stat=status(3))
             ASSERT(status(3)==0)
          enddo
          deallocate( spinor_grads, stat=status(3) )
          if ( status(3) .ne. 0 ) call error_handler( &
               "orbital_free: deallocate of grads failed")
       endif
    endif ! spinor_grads

    if ( present(sec_der) ) then
       do i_ir = 1, size(sec_der,1)
          MEMLOG(-size(sec_der(i_ir)%o))
          deallocate( sec_der(i_ir)%o, stat=status(5) )
          ASSERT(status(5).eq.0)
          status(5)=1
       enddo

       deallocate(sec_der, stat=status(60) )
       ASSERT(status(60).eq.0)
       status(60)=1
    endif
    if ( present(nuc_grads) ) then
       do i_unique=1,n_unique_atoms
          do i_ir = 1, size(nuc_grads,2)
             MEMLOG(-size(nuc_grads(i_unique,i_ir)%o))
             deallocate( nuc_grads(i_unique,i_ir)%o, stat=status(6) )
           ASSERT(status(6).eq.0)
           status(6)=1
          enddo
       end do
       deallocate( nuc_grads, stat=status(61) )
       ASSERT(status(61).eq.0)
       status(61)=1
    end if
    if ( present(nuc_sec_der) ) then
       do i_unique=1,n_unique_atoms
          do i_ir = 1, size(nuc_sec_der,2)
             MEMLOG(-size(nuc_sec_der(i_unique,i_ir)%o))
             deallocate( nuc_sec_der(i_unique,i_ir)%o, stat=status(7) )
          ASSERT(status(7).eq.0)
          status(7)=1
          enddo
       end do
       deallocate( nuc_sec_der, stat=status(62) )
       ASSERT(status(62).eq.0)
       status(62)=1
    end if

    if ( present(nuc_3rd_der) ) then
       do i_unique=1,n_unique_atoms
          do i_ir = 1, size(nuc_3rd_der,2)
             MEMLOG(-size(nuc_3rd_der(i_unique,i_ir)%o))
             deallocate( nuc_3rd_der(i_unique,i_ir)%o, stat=status(52) )
          ASSERT(status(52).eq.0)
       status(52)=1
          enddo
       end do
       deallocate( nuc_3rd_der, stat=status(51) )
       ASSERT(status(51).eq.0)
       status(51)=1
    end if

    if ( present(orb_xc) ) then
       MEMLOG(-size(orb_xc%o))
       deallocate( orb_xc%o, stat=status(9) )
       if ( status(9) .ne. 0 ) call error_handler( &
            "orbital_free: deallocate of orb_xc%o failed")
       ! deallocation of uncont_orbs_xc moved to orbital_allocate
    endif
  end subroutine orbital_free

#ifdef WITH_CORE_DENS
  subroutine core_orbital_free(orbs_ob_core,grads_core,sec_der_core,&
                               nuc_grads_core,nuc_sec_der_core)
    !  Purpose: deallocate orbs_core
    use orbitalstore_module, only: core_orbital_type, core_orbital_gradient_type, &
        core_orbital_sec_der_type, core_orbital_nuc_gradient_type,  &
        core_orbital_nuc_sec_der_type
    use unique_atom_module ! desription of unique atoms, must be calculated before
    implicit none
    !------------ Declaration of formal parameters -------------
    type(core_orbital_type),              optional :: orbs_ob_core
    type(core_orbital_gradient_type),     optional :: grads_core
    type(core_orbital_sec_der_type),      optional :: sec_der_core
    type(core_orbital_nuc_gradient_type), pointer, optional :: &
                                                 nuc_grads_core(:)
    type(core_orbital_nuc_sec_der_type), pointer, optional ::  &
                                                 nuc_sec_der_core(:)
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer              :: i_unique
    !------------ Executable code ------------------------------
    if ( present(orbs_ob_core) ) then
          deallocate( orbs_ob_core%o, stat=status(11) )
          if ( status(11) .ne. 0 ) call error_handler( &
               "core_orbital_free: deallocate of orbs_ob%o failed")
    endif

    if ( present(grads_core) ) then
          deallocate( grads_core%o, stat=status(12) )
          if ( status(12) .ne. 0 ) call error_handler( &
               "core_orbital_free: deallocate of grads%o failed")
    endif
    if ( present(sec_der_core) ) then
          deallocate( sec_der_core%o, stat=status(13) )
          if ( status(13) .ne. 0 ) call error_handler( &
               "core_orbital_free: deallocate of sec_der%o failed")
    endif
    if ( present(nuc_grads_core) ) then
       do i_unique=1,n_unique_atoms
             deallocate( nuc_grads_core(i_unique)%o, stat=status(14) )
             if ( status(14) .ne. 0 ) call error_handler( &
                  "core_orbital_free: deallocate of grads%o failed")
       end do
       deallocate( nuc_grads_core, stat=status(14) )
       if ( status(14) .ne. 0 ) call error_handler( &
            "core_orbital_free: deallocate of nuc_grads failed")
    end if
    if ( present(nuc_sec_der_core) ) then
       do i_unique=1,n_unique_atoms
             deallocate( nuc_sec_der_core(i_unique)%o, stat=status(15) )
             if ( status(15) .ne. 0 ) call error_handler( &
                  "core_orbital_free: deallocate of nuc_sec_der%o failed")
       end do
       deallocate( nuc_sec_der_core, stat=status(15) )
       if ( status(15) .ne. 0 ) call error_handler( &
            "core_orbital_free: deallocate of nuc_sec_der failed")
    end if

  end subroutine core_orbital_free
#endif


  subroutine orbital_setup(vector_length, do_gradients, do_sec_der, do_3rd_der)
    !  Purpose: initialising and allocating module-privat variables
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use solid_harmonics_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: vector_length
    logical, optional, intent(in) :: do_gradients, do_sec_der,do_3rd_der
    ! default: false
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)       :: i_ua,i_ea,N_ea,i_l
    type(unique_atom_type), pointer   :: ua
    logical :: do_grad_local, do_sec_der_local, do_3rd_der_local
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    if ( present(do_gradients) ) then
       do_grad_local = do_gradients
    else
       do_grad_local = .false.
    endif

    if ( present(do_sec_der) ) then
       do_sec_der_local = do_sec_der
    else
       do_sec_der_local = .false.
    endif

    if ( present(do_3rd_der) ) then
       do_3rd_der_local = do_3rd_der
    else
       do_3rd_der_local = .false.
    endif

    ! we require gradients to calculate secound derivatives
    if (do_3rd_der_local) do_sec_der_local = .true.
    if (do_sec_der_local) do_grad_local = .true.

    !
    ! FIXME: honestly one would need to re-allocate global arrays
    !        in the case ANY of the
    !
    !                (vector_length, N_max_exponents, N_max_equal_atoms)
    !
    !        changing from the previous values.

    if  ( vec_len .ne. vector_length .and. memory_allocated) then
       !
       ! memory_allocated is set to .false. as a result:
       !
       call orbital_shutdown()
       ASSERT(.not.memory_allocated)
    endif

    !
    ! Store vector_length:
    !
    vec_len = vector_length

    !
    ! Determine N_max_exponents, N_equal_atoms, do_glob_con_xx and N_max_l_xx:
    !
!   if ( you dont mind spending a few CPU cycles ) then
       N_max_exponents = 0
       N_max_equal_atoms = 0
       do_glob_con_ch = .false.
       do_glob_con_xc = .false.
       do i_ua = 1, N_unique_atoms
          ua => unique_atoms(i_ua)
          do i_l = 0, ua%lmax_ob
             N_max_exponents = max( N_max_exponents, ua%l_ob(i_l)%N_exponents )
          enddo
          do i_l = 0, ua%lmax_ch
             N_max_exponents = max( N_max_exponents, ua%l_ch(i_l)%N_exponents )
          enddo
          N_max_exponents = max( N_max_exponents, ua%r2_ch%N_exponents )
          do i_l = 0, ua%lmax_xc
             N_max_exponents = max( N_max_exponents, ua%l_xc(i_l)%N_exponents )
          enddo
          N_max_exponents = max( N_max_exponents, ua%r2_xc%N_exponents )
          N_max_equal_atoms = max( N_max_equal_atoms, ua%N_equal_atoms )
          do_glob_con_ch = do_glob_con_ch .or. ua%N_glob_cons_ch .gt. 0
          do_glob_con_xc = do_glob_con_xc .or. ua%N_glob_cons_xc .gt. 0
       enddo
!   endif

    ! initialise solid_harmonics_module
    call solid_harmonics_setup( vector_length )


    if ( .not. memory_allocated ) then

       ! allocate solid_harmonics
       allocate( solid_harmonics(N_unique_atoms), stat=status(16) )
       ASSERT(status(16).eq.0)
       status(16)=1

       if ( do_3rd_der_local ) then
          do i_ua = 1, N_unique_atoms
             N_ea = unique_atoms(i_ua)%N_equal_atoms
             allocate( solid_harmonics(i_ua)%sh(N_ea), &
                  solid_harmonics(i_ua)%shg(N_ea), &
                  solid_harmonics(i_ua)%shs(N_ea), &
                  solid_harmonics(i_ua)%sht(N_ea), &
                  stat=status(53) )
             ASSERT(status(53).eq.0)
             do i_ea = 1,N_ea
                call solid_harmonics_allocate( vec_len, &
                     unique_atoms(i_ua)%lmax_all, &
                     solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea), &
                     solid_harmonics(i_ua)%shs(i_ea), &
                     solid_harmonics(i_ua)%sht(i_ea) )
             enddo
          enddo

       elseif ( do_sec_der_local ) then
          do i_ua = 1, N_unique_atoms
             N_ea = unique_atoms(i_ua)%N_equal_atoms
             allocate( solid_harmonics(i_ua)%sh(N_ea), &
                  solid_harmonics(i_ua)%shg(N_ea), &
                  solid_harmonics(i_ua)%shs(N_ea), &
                  stat=status(17) )
             ASSERT(status(17).eq.0)
             do i_ea = 1,N_ea
                call solid_harmonics_allocate( vec_len, &
                     unique_atoms(i_ua)%lmax_all, &
                     solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea), &
                     solid_harmonics(i_ua)%shs(i_ea) )
             enddo
          enddo
       status(17)=1
       elseif ( do_grad_local ) then
          do i_ua = 1, N_unique_atoms
             N_ea = unique_atoms(i_ua)%N_equal_atoms
             allocate( solid_harmonics(i_ua)%sh(N_ea), &
                  solid_harmonics(i_ua)%shg(N_ea), &
                  stat=status(18) )
        ASSERT(status(18).eq.0)
             do i_ea = 1,N_ea
                call solid_harmonics_allocate( vec_len, &
                     unique_atoms(i_ua)%lmax_all, &
                     solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea) )
             enddo
          enddo
       else
          do i_ua = 1, N_unique_atoms
             N_ea = unique_atoms(i_ua)%N_equal_atoms
             allocate( solid_harmonics(i_ua)%sh(N_ea), stat=status(19) )
             ASSERT(status(19).eq.0)
             do i_ea = 1,N_ea
                call solid_harmonics_allocate( vec_len, &
                     unique_atoms(i_ua)%lmax_all, &
                     solid_harmonics(i_ua)%sh(i_ea) )
             enddo
          enddo
       endif

       ! allocate prim_orb_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( prim_orb_fac(vec_len, N_max_exponents, N_max_equal_atoms), stat=status(20) )
       ASSERT(status(20).eq.0)
       MEMLOG(size(prim_orb_fac))

       ! allocate cont_orb_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( cont_orb_fac(vec_len, N_max_exponents, N_max_equal_atoms), stat=status(21) )
       ASSERT(status(21).eq.0)
       MEMLOG(size(cont_orb_fac))


       if ( do_grad_local ) then

          ! allocate prim_grad_fac(vec_len,N_max_exponents,N_max_equal_atoms)
          allocate( prim_grad_fac(vec_len, N_max_exponents, N_max_equal_atoms), stat=status(22) )
          ASSERT(status(22).eq.0)
          MEMLOG(size(prim_grad_fac))

          ! allocate cont_grad_fac(vec_len,N_max_exponents,N_max_equal_atoms)
          allocate( cont_grad_fac(vec_len, N_max_exponents, N_max_equal_atoms), stat=status(23) )
          ASSERT(status(23).eq.0)
          MEMLOG(size(cont_grad_fac))

          grads_allocated = .true.

       else
          ! dummy allocation to facilitate call of combine_rad_and_y(00/lm)
          allocate( prim_grad_fac(1,1,1), cont_grad_fac(1,1,1), stat=status(22) )
          MEMLOG(size(prim_grad_fac)+size(cont_grad_fac))
          ASSERT(status(22).eq.0)
       endif


       if ( do_sec_der_local ) then

          ! allocate prim_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms)
          allocate( prim_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms), stat=status(24) )
          ASSERT(status(24).eq.0)
          MEMLOG(size(prim_sd_fac))

          allocate( cont_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms), stat=status(25) )
          ASSERT(status(25).eq.0)
          MEMLOG(size(cont_sd_fac))

          sec_der_allocated = .true.

       else
          ! dummy allocation to facilitate call of combine_rad_and_y(00/lm)
          allocate( prim_sd_fac(1,1,1), cont_sd_fac(1,1,1), stat=status(24) )
          ASSERT(status(24).eq.0)
          MEMLOG(size(prim_sd_fac)+size(cont_sd_fac))
       endif

       if ( do_3rd_der_local ) then

          ! allocate prim_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms)
          allocate( prim_td_fac(vec_len,N_max_exponents,N_max_equal_atoms), stat=status(54) )
          ASSERT(status(54).eq.0)
          MEMLOG(size(prim_td_fac))

          allocate( cont_td_fac(vec_len,N_max_exponents,N_max_equal_atoms), stat=status(55) )
          ASSERT(status(55).eq.0)
          MEMLOG(size(cont_td_fac))

          third_der_allocated = .true.

       else
          ! dummy allocation to facilitate call of combine_rad_and_y(00/lm)
          allocate( prim_td_fac(1,1,1), cont_td_fac(1,1,1), stat=status(54) )
          ASSERT(status(54).eq.0)
          MEMLOG(size(prim_td_fac)+size(cont_td_fac))
       endif

       ! for global contractions
       if ( do_glob_con_ch ) call glob_con_setup(uncont_orbs_ch, &
                                                 unique_atoms(:)%lmax_ch)
       if ( do_glob_con_xc ) call glob_con_setup(uncont_orbs_xc, &
                                                 unique_atoms(:)%lmax_xc)

       memory_allocated = .true.

    endif

  end subroutine orbital_setup


  subroutine glob_con_setup(uncont_orbs,lmax)
    !  Purpose: allocates the root components for the intermediates which
    !  hold the uncontracted fit functions and their derivatives
    !
    !  called by orbital_setup
    use unique_atom_module ! desription of unique atoms, must be calculated before
    implicit none
    type(glob_con_intermediate_type), pointer    :: uncont_orbs(:)
    integer                         , intent(in) :: lmax(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer :: i_ua
    !------------ Executable code ------------------------------------
    allocate( uncont_orbs(N_unique_atoms), stat=status(57) )
    if ( status(57) /= 0 ) call error_handler("glob_con_setup:&
         & allocation of uncont_orbs failed")

    do i_ua=1,N_unique_atoms
       allocate( uncont_orbs(i_ua)%l(-1:lmax(i_ua)), stat=status(26) )
       if ( status(26) /= 0 ) call error_handler("glob_con_setup:&
            & allocation of uncont_orbs(i_ua)%l failed")
    end do! loop over unique atoms

  end subroutine glob_con_setup


  subroutine orbital_shutdown()
    !
    ! Purpose: deallocating module-privat variables
    !
    ! Sets global flag memory_allocated = false
    !
    use solid_harmonics_module, only: solid_harmonics_free, solid_harmonics_shutdown
    implicit none
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)             :: i_ua,i_ea
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    if ( .not. memory_allocated ) then
        WARN('why calling?')
        return
    endif

       ! free solid_harmonics
       if ( third_der_allocated ) then
          do i_ua = 1, size(solid_harmonics)
             do i_ea = 1, size(solid_harmonics(i_ua)%sh)
                call solid_harmonics_free( solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea), &
                     solid_harmonics(i_ua)%shs(i_ea), &
                     solid_harmonics(i_ua)%sht(i_ea)  )
             enddo
             deallocate( solid_harmonics(i_ua)%sh, &
                  solid_harmonics(i_ua)%shg, &
                  solid_harmonics(i_ua)%shs, &
                  solid_harmonics(i_ua)%sht, &
                  stat=status(53) )
          ASSERT(status(53).eq.0)
          status(53)=1
          enddo
       elseif ( sec_der_allocated ) then
          do i_ua = 1, size(solid_harmonics)
             do i_ea = 1, size(solid_harmonics(i_ua)%sh)
                call solid_harmonics_free( solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea), &
                     solid_harmonics(i_ua)%shs(i_ea) )
             enddo
             deallocate( solid_harmonics(i_ua)%sh, &
                  solid_harmonics(i_ua)%shg, &
                  solid_harmonics(i_ua)%shs, &
                  stat=status(17) )
          ASSERT(status(17).eq.0)
          status(17)=1
          enddo
       elseif ( grads_allocated ) then
          do i_ua = 1, size(solid_harmonics)
             do i_ea = 1, size(solid_harmonics(i_ua)%sh)
                call solid_harmonics_free( solid_harmonics(i_ua)%sh(i_ea), &
                     solid_harmonics(i_ua)%shg(i_ea) )
             enddo
             deallocate( solid_harmonics(i_ua)%sh,  solid_harmonics(i_ua)%shg, &
                  stat=status(18) )
          ASSERT(status(18).eq.0)
          status(18)=1
          enddo
       else
          do i_ua = 1, size(solid_harmonics)
             do i_ea = 1, size(solid_harmonics(i_ua)%sh)
                call solid_harmonics_free( solid_harmonics(i_ua)%sh(i_ea) )
             enddo
             deallocate( solid_harmonics(i_ua)%sh, stat=status(19) )
             ASSERT(status(19).eq.0)
             status(19)=1
          enddo
       endif

       deallocate( solid_harmonics, stat=status(16) )
       ASSERT(status(16).eq.0)
       status(16)=1

       ! shut down solid_harmonics module
       call solid_harmonics_shutdown()

       ! free prim_orb_fac
       MEMLOG(-size(prim_orb_fac))
       deallocate( prim_orb_fac, stat=status(20) )
       ASSERT(status(20).eq.0)
       status(20)=1

       ! free cont_orb_fac
       MEMLOG(-size(cont_orb_fac))
       deallocate( cont_orb_fac, stat=status(21) )
       ASSERT(status(21).eq.0)
       status(21)=1


       ! free prim_grad_fac
       MEMLOG(-size(prim_grad_fac))
       deallocate( prim_grad_fac, stat=status(22) )
       ASSERT(status(22).eq.0)
       status(22)=1

       ! free cont_grad_fac
       MEMLOG(-size(cont_grad_fac))
       deallocate( cont_grad_fac, stat=status(23) )
       ASSERT(status(23).eq.0)
       status(23)=1

       grads_allocated = .false.

       ! free prim_sd_fac
       MEMLOG(-size(prim_sd_fac))
       deallocate( prim_sd_fac, stat=status(24) )
       ASSERT(status(24).eq.0)
       status(24)=1

       ! free cont_sd_fac
       MEMLOG(-size(cont_sd_fac))
       deallocate( cont_sd_fac, stat=status(25) )
       ASSERT(status(25).eq.0)
       status(25)=1

       sec_der_allocated = .false.

       ! free prim_td_fac
       MEMLOG(-size(prim_td_fac))
       deallocate( prim_td_fac, stat=status(54) )
       ASSERT(status(54).eq.0)
       status(54)=1

       ! free cont_sd_fac
       MEMLOG(-size(cont_td_fac))
       deallocate( cont_td_fac, stat=status(55) )
       ASSERT(status(55).eq.0)
       status(55)=1

       third_der_allocated = .false.

       ! for global contractions
       if ( do_glob_con_ch ) call glob_con_shutdown(uncont_orbs_ch)
       if ( do_glob_con_xc ) call glob_con_shutdown(uncont_orbs_xc)

       memory_allocated = .false.
  end subroutine orbital_shutdown


  subroutine glob_con_shutdown(uncont_orbs)
    !  Purpose: deallocates the root components for the intermediates which
    !  have hold the uncontracted fit functions and their derivatives
    !
    !  called by orbital_shutdown
    implicit none
    type(glob_con_intermediate_type), pointer    :: uncont_orbs(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer :: i_ua
    !------------ Executable code ------------------------------------
    do i_ua = 1, size(uncont_orbs)
       deallocate( uncont_orbs(i_ua)%l, stat=status(26) )
       if ( status(26)/= 0 ) call error_handler("glob_con_shutdown:&
            & allocation of uncont_orbs(i_ua)%l failed")
    end do! loop over unique atoms

    deallocate( uncont_orbs, stat=status(26) )
    if ( status(26) /= 0 ) call error_handler("glob_con_shutdown:&
         & deallocation of uncont_orbs failed")

  end subroutine glob_con_shutdown
  !*************************************************************


  !*************************************************************
  subroutine uncont_fit_fct_allocate(type,&
             fcts,grads,sec_ders,nuc_grads,sec_nuc_ders,orbs)
    !  Purpose: allocates the intermediates holding the uncontracted
    !  fit functions and their derivatives on a set of grid points
    !  which are used to evaluate the global contractions
    !
    !  called by fit_fct_allocate and orbital_allocate
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use symmetry_data_module ! description of irreps
    implicit none
    character(len=2) , intent(in) :: type ! "ch" or "xc"
    logical, optional, intent(in) :: fcts, grads, sec_ders, &
                                     nuc_grads, sec_nuc_ders
    logical, optional, intent(in) :: orbs ! kept for compatibility
    !** End of interface *****************************************
!AG #ifdef FPP_NOMDA
!AG     call error_handler("om/uncont_fit_fct_allocate: recompile w/o FPP_NOMDA")
!AG #else
    !------------ Declaration of local variables ---------------------
    integer :: i_ir, i_ua, lmax, i_bas
    logical :: is_charge_fit, moving_atom, alloc_fcts, alloc_grads, &
               alloc_sec_ders, alloc_nuc_grads, alloc_sec_nuc_ders, &
               alloc_orbs
    type(unique_atom_type),             pointer :: ua
    type(unique_atom_basis_type),       pointer :: uab
    type(unique_atom_partner_type),     pointer :: uap
    type(glob_con_intermediate_l_type), pointer :: uncont_orb(:)
    type(glob_con_intermediate_l_type), pointer :: uncont_orb_l
    !------------ Executable code ------------------------------------
    select case (type)
    case ("ch")
       is_charge_fit = .true.
    case ("xc")
       is_charge_fit = .false.
    case default
       call error_handler( &
            "uncont_fit_fct_allocate: unknown type of fit functions")
    end select

    if (present(orbs)) then
       alloc_orbs = orbs
    else
       alloc_orbs = .false.
    endif
    if (present(fcts)) then
       alloc_fcts = fcts
    else
       alloc_fcts = .false.
    endif
    if (present(grads)) then
       alloc_grads = grads
    else
       alloc_grads = .false.
    endif
    if (present(sec_ders)) then
       alloc_sec_ders = sec_ders
    else
       alloc_sec_ders = .false.
    endif
    if (present(nuc_grads)) then
       alloc_nuc_grads = nuc_grads
    else
       alloc_nuc_grads = .false.
    endif
    if (present(sec_nuc_ders)) then
       alloc_sec_nuc_ders = sec_nuc_ders
    else
       alloc_sec_nuc_ders = .false.
    endif

    i_ir = get_totalsymmetric_irrep()

    ! loop over unique atoms
    do i_ua=1,N_unique_atoms
       ua => unique_atoms(i_ua)
       moving_atom = ua%moving_atom > 0
       if ( is_charge_fit ) then
          lmax = ua%lmax_ch
          uncont_orb => uncont_orbs_ch(i_ua)%l
       else
          lmax = ua%lmax_xc
          uncont_orb => uncont_orbs_xc(i_ua)%l
       endif

       ! loop over basis types
       do i_bas = -1,lmax
          select case (i_bas)
          case (-1) ! s-type
             if (is_charge_fit) then
                uab => ua%l_ch(0)
             else
                uab => ua%l_xc(0)
             endif
             uap => ua%symadapt_partner(i_ir,0)
          case (0) ! r2-type
             if (is_charge_fit) then
                uab => ua%r2_ch
             else
                uab => ua%r2_xc
             endif
             uap => ua%symadapt_partner(i_ir,0)
          case default
             if (is_charge_fit) then
                uab => ua%l_ch(i_bas)
             else
                uab => ua%l_xc(i_bas)
             endif
             uap => ua%symadapt_partner(i_ir,i_bas)
          end select
          uncont_orb_l => uncont_orb(i_bas)

          if (alloc_orbs) then
             allocate( uncont_orb_l%oo(vec_len, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(27) )
             if ( status(27) /= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%oo failed")
             status(27)=1
          endif
          if (alloc_fcts) then
             allocate( uncont_orb_l%o(vec_len, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(28) )
             if ( status(28) /= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%o failed")
             status(28)=1
          endif
          if (alloc_grads) then
             allocate( uncont_orb_l%og(vec_len,3, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(29) )
             if ( status(29)/= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%og failed")
             status(29)=1
          endif
          if (alloc_sec_ders) then
             allocate( uncont_orb_l%osd(vec_len,6, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(30) )
             if ( status(30) /= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%osd failed")
             status(30)=1
          endif
          if (alloc_nuc_grads .and. moving_atom) then
             allocate( uncont_orb_l%ong(vec_len,3,ua%N_equal_atoms, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(31) )
             if ( status(31) /= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%ong failed")
             status(31)=1
          endif
          if (alloc_sec_nuc_ders .and. moving_atom) then
             allocate( uncont_orb_l%onsd(vec_len,6,ua%N_equal_atoms, &
                       uab%N_exponents,uap%N_independent_fcts), stat=status(32) )
             if ( status(32) /= 0 ) call error_handler("uncont_fit_fct_allocate:&
                  & allocation of uncont_orbs(i_ua)%l(i_bas)%onsd failed")
             status(32)=1
          endif

       end do! loop over i_bas
    end do! loop over unique atoms
!AG#endif
  end subroutine uncont_fit_fct_allocate
  !*************************************************************


  !*************************************************************
  subroutine uncont_fit_fct_free(type,&
             fcts,grads,sec_ders,nuc_grads,sec_nuc_ders,orbs)
    !  Purpose: deallocates the intermediates holding the uncontracted
    !  fit functions and their derivatives on a set of grid points
    !  which have been used to evaluate the global contractions
    !
    !  called by fit_fct_free and orbital_free
    use unique_atom_module ! desription of unique atoms, must be calculated before
    implicit none
    character(len=2) , intent(in) :: type ! "ch" or "xc"
    logical, optional, intent(in) :: fcts, grads, sec_ders, &
                                     nuc_grads, sec_nuc_ders
    logical, optional, intent(in) :: orbs ! kept for compatibility
    !** End of interface *****************************************
#ifdef FPP_NOMDA
    call error_handler("om/uncont_fit_fct_free: recompile w/o FPP_NOMDA")
#else
    !------------ Declaration of local variables ---------------------
    integer :: i_ua, lmax, i_bas
    logical :: is_charge_fit, moving_atom, dealloc_fcts, dealloc_grads, &
               dealloc_sec_ders, dealloc_nuc_grads, dealloc_sec_nuc_ders, &
               dealloc_orbs
    type(unique_atom_type),             pointer :: ua
    type(glob_con_intermediate_l_type), pointer :: uncont_orb(:)
    type(glob_con_intermediate_l_type), pointer :: uncont_orb_l
    !------------ Executable code ------------------------------------
    select case (type)
    case ("ch")
       is_charge_fit = .true.
    case ("xc")
       is_charge_fit = .false.
    case default
       call error_handler( &
            "uncont_fit_fct_free: unknown type of fit functions")
    end select

    if (present(orbs)) then
       dealloc_orbs = orbs
    else
       dealloc_orbs = .false.
    endif
    if (present(fcts)) then
       dealloc_fcts = fcts
    else
       dealloc_fcts = .false.
    endif
    if (present(grads)) then
       dealloc_grads = grads
    else
       dealloc_grads = .false.
    endif
    if (present(sec_ders)) then
       dealloc_sec_ders = sec_ders
    else
       dealloc_sec_ders = .false.
    endif
    if (present(nuc_grads)) then
       dealloc_nuc_grads = nuc_grads
    else
       dealloc_nuc_grads = .false.
    endif
    if (present(sec_nuc_ders)) then
       dealloc_sec_nuc_ders = sec_nuc_ders
    else
       dealloc_sec_nuc_ders = .false.
    endif

    ! loop over unique atoms
    do i_ua=1,N_unique_atoms
       ua => unique_atoms(i_ua)
       moving_atom = ua%moving_atom > 0

       if ( is_charge_fit ) then
          lmax = ua%lmax_ch
          uncont_orb => uncont_orbs_ch(i_ua)%l
       else
          lmax = ua%lmax_xc
          uncont_orb => uncont_orbs_xc(i_ua)%l
       endif

       ! loop over basis types
       do i_bas = -1,lmax
          uncont_orb_l => uncont_orb(i_bas)

          if (dealloc_orbs) then
             deallocate( uncont_orb_l%oo, stat=status(27) )
             if ( status(27) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%oo failed")
          endif
          if (dealloc_fcts) then
             deallocate( uncont_orb_l%o, stat=status(28) )
             if ( status(28) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%o failed")
          endif
          if (dealloc_grads) then
             deallocate( uncont_orb_l%og, stat=status(29) )
             if ( status(29) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%og failed")
          endif
          if (dealloc_sec_ders) then
             deallocate( uncont_orb_l%osd, stat=status(30) )
             if ( status(30) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%osd failed")
          endif
          if (dealloc_nuc_grads .and. moving_atom) then
             deallocate( uncont_orb_l%ong, stat=status(31) )
             if ( status(31) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%ong failed")
          endif
          if (dealloc_sec_nuc_ders .and. moving_atom) then
             deallocate( uncont_orb_l%onsd, stat=status(32) )
             if ( status(32) /= 0 ) call error_handler("uncont_fit_fct_free:&
                  & deallocation of uncont_orbs(i_ua)%l(i_bas)%onsd failed")
             status(32)=1
          endif

       end do! loop over i_bas
    end do! loop over unique atoms
#endif
  end subroutine uncont_fit_fct_free
  !*************************************************************

  subroutine calc_spinors(vl,orb,sporb)
    use type_module,    only: IK=>i4_kind, RK=>r8_kind
    use error_module,   only: error
    use clebsch_gordan, only: cg => vsu2cg_eliminated,&
         &                    prod_bas ! a type
    use dimensions,     only: uaLdims => uaL_vec_dims,&
         &                    uaL_max
    implicit none
    integer(IK),intent(in)        :: vl ! vector length
    type(orbital_type),intent(in) ::   orb(:) ! (n_vec_irrep)
    type(spinor_type),intent(out) :: sporb(:) ! (n_proj_irrep)
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: n_v_irr,n_p_irr,ip,iv
    integer(IK)             :: pcount
    integer(IK),allocatable :: vcount(:)
    integer(IK)             :: uaL,mult,ic,pp,pv,ap,ep,av,ev,nf
    type(prod_bas),pointer  :: cf


    n_v_irr = size(orb)
    n_p_irr = size(sporb)

    allocate(vcount(n_v_irr),STAT=memstat)
    call error(memstat,"om/calc_spinors: alloc failed")

    ip_:do ip=1,n_p_irr
       pcount = 0
       vcount = 0
       uaL_:do uaL=1,uaL_max
          iv_:do iv=1,n_v_irr

             nf = uaLdims(uaL,iv)
             av = vcount(iv) +  1
             ev = vcount(iv) + nf

             mult_:do mult=1,cg(ip,iv)%mult
                ap = pcount +  1
                ep = pcount + nf
                cf => cg(ip,iv)%sub(mult)

                !
                ! now F(mu,gamma) = SUM_pi{ C(mu,gamma,pi) * F(pi) }
                ! or  F(ic,pp)    = SUM_pv{ C(ic,pp,   pv) * F(pv) }
                ! where mu/ic    -- spinor component
                !       gamma/pp -- projectve partner
                !       pi/pv    -- vector partner
                !       C(mu,gamma,pi)/C(ic,pp,pv) -- CG for Proj = Vec * SU(2)
                !
                ic_: do ic = 1, size(cf%z, 3) ! nc==2, a spinor component
                   pp_: do pp = 1, size(cf%z, 1) ! na, proj partner
                      pv_: do pv = 1, size(cf%z, 2) ! nb, vec partner

                         sporb(ip)%o(:vl, ap:ep, pp, RE, ic) = sporb(ip)%o(:vl, ap:ep, pp, RE, ic) &
                              + cf%re(pp, pv, ic) * orb(iv)%o(:vl, av:ev, pv)

                         sporb(ip)%o(:vl, ap:ep, pp, IM, ic) = sporb(ip)%o(:vl, ap:ep, pp, IM, ic) &
                              + cf%im(pp, pv, ic) * orb(iv)%o(:vl, av:ev, pv)

                      enddo pv_
                   enddo pp_
                enddo ic_
                pcount = pcount + nf
             enddo mult_

             vcount(iv) = vcount(iv) + nf
          enddo iv_
       enddo uaL_
    enddo ip_

    deallocate(vcount,STAT=memstat)
    call error(memstat,"om/calc_spinors: dealloc failed")

  end subroutine calc_spinors

  subroutine radial(vl, r2, exps, norms, n_unc, contrs, &
       & rad0, &
       & rad1, &
       & rad2, &
       & rad3 )
    ! computes radial parts
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(in)    :: r2(:)    ! (>=vl)
    real(r8_kind), intent(in)    :: exps(:)  ! (n_exps)
    real(r8_kind), intent(in)    :: norms(:) ! (n_exps)
    integer(i4_kind), intent(in) :: n_unc
    real(r8_kind), intent(in)    :: contrs(:,:) ! (n_exps,n_contrs)
    real(r8_kind), intent(out)   :: rad0(:,:)   ! (>=vl,>=n_exps)
    real(r8_kind), intent(out)   :: rad1(:,:)   ! (>=vl,>=n_exps)
    real(r8_kind), intent(out)   :: rad2(:,:)   ! (>=vl,>=n_exps)
    real(r8_kind), intent(out)   :: rad3(:,:)   ! (>=vl,>=n_exps)
    optional :: rad1
    optional :: rad2
    optional :: rad3
    ! *** end of interface ***

    integer(i4_kind) :: n_exps, n_contrs
    logical          :: D1,D2,D3
    integer(i4_kind) :: deriv
#ifdef NEVER
    integer(i4_kind) :: total = 0
    integer(i4_kind) :: zeros = 0
#endif

    deriv = 0
    D1 = present(rad1)
    if(D1) deriv = deriv + 1
    D2 = present(rad2)
    if(D2) deriv = deriv + 1
    ASSERT(D2.eqv.deriv.eq.2)
    D3 = present(rad3)
    if(D3) deriv = deriv + 1
    ASSERT(D3.eqv.deriv.eq.3)

    ASSERT(vl<=size(rad0,1))
    ASSERT(vl<=size(r2))

    n_exps   = size(contrs,1)
    n_contrs = size(contrs,2)

    ASSERT(n_exps==size(exps))
    ASSERT(n_exps==size(norms))
    ASSERT(n_exps<=size(rad0,2))

    select case(deriv)
    case (0)
       ! calculate primitive orbital factors for all exponents and gridpoints
       call radial0()
       call contract(n_unc, rad0, contrs)
    case (1)
       ! calculate primitive orbital factors for all exponents and gridpoints
       call radial1()
       call contract(n_unc, rad0, contrs)
       call contract(n_unc, rad1, contrs)
    case (2)
       ! calculate primitive orbital factors for all exponents and gridpoints
       call radial2()
       call contract(n_unc, rad0, contrs)
       call contract(n_unc, rad1, contrs)
       call contract(n_unc, rad2, contrs)
    case (3)
       ! calculate primitive orbital factors for all exponents and gridpoints
       call radial3()
       call contract(n_unc, rad0, contrs)
       call contract(n_unc, rad1, contrs)
       call contract(n_unc, rad2, contrs)
       call contract(n_unc, rad3, contrs)
    end select
#ifdef NEVER
      total = total +  size(rad0(:vl,:n_unc+n_contrs))
      zeros = zeros + count(rad0(:vl,:n_unc+n_contrs)==ZERO)
      DPRINT 'om: total=',total,' zeros=',zeros
#endif
  contains

    subroutine radial0()
      integer(i4_kind) :: i_exp
      real(r8_kind)    :: norm, alpha
      real(r8_kind) :: ar2(vl)
      do i_exp = 1, n_exps
         alpha = exps(i_exp)
         norm  = norms(i_exp)
         ar2 = alpha * r2(:vl)
         where( ar2 .lt. 50.0_r8_kind )
            rad0(:vl,i_exp) = norm * exp( - ar2(:vl) )
         elsewhere
            rad0(:vl,i_exp) = 0.0_r8_kind
         endwhere
      enddo
    end subroutine radial0

    subroutine radial1()
      integer(i4_kind) :: i_exp
      real(r8_kind)    :: norm, alpha
      real(r8_kind) :: ar2(vl)
      do i_exp = 1, n_exps
         alpha = exps(i_exp)
         norm  = norms(i_exp)
         ar2 = alpha * r2(:vl)
         where( ar2 .lt. 50.0_r8_kind )
            rad0(:vl,i_exp) = norm * exp( - ar2(:) )
            rad1(:vl,i_exp) = (-TWO * alpha) * rad0(:vl,i_exp)
         elsewhere
            rad0(:vl,i_exp) = 0.0_r8_kind
            rad1(:vl,i_exp) = 0.0_r8_kind
         endwhere
      enddo
    end subroutine radial1

    subroutine radial2()
      integer(i4_kind) :: i_exp
      real(r8_kind)    :: norm, alpha
      real(r8_kind) :: ar2(vl)
      do i_exp = 1, n_exps
         alpha = exps(i_exp)
         norm  = norms(i_exp)
         ar2 = alpha * r2(:vl)
         where( ar2 .lt. 50.0_r8_kind )
            rad0(:vl,i_exp) = norm * exp( - ar2(:) )
            rad1(:vl,i_exp) = (-TWO * alpha) * rad0(:vl,i_exp)
            rad2(:vl,i_exp) = (-TWO * alpha) * rad1(:vl,i_exp)
         elsewhere
            rad0(:vl,i_exp) = 0.0_r8_kind
            rad1(:vl,i_exp) = 0.0_r8_kind
            rad2(:vl,i_exp) = 0.0_r8_kind
         endwhere
      enddo
    end subroutine radial2

    subroutine radial3()
      integer(i4_kind) :: i_exp
      real(r8_kind)    :: norm, alpha
      real(r8_kind) :: ar2(vl)
      do i_exp = 1, n_exps
         alpha = exps(i_exp)
         norm  = norms(i_exp)
         ar2 = alpha * r2(:vl)
         where( ar2 .lt. 50.0_r8_kind )
            rad0(:vl,i_exp) = norm * exp( - ar2(:) )
            rad1(:vl,i_exp) = (-TWO * alpha) * rad0(:vl,i_exp)
            rad2(:vl,i_exp) = (-TWO * alpha) * rad1(:vl,i_exp)
            rad3(:vl,i_exp) = (-TWO * alpha) * rad2(:vl,i_exp)
         elsewhere
            rad0(:vl,i_exp) = 0.0_r8_kind
            rad1(:vl,i_exp) = 0.0_r8_kind
            rad2(:vl,i_exp) = 0.0_r8_kind
            rad3(:vl,i_exp) = 0.0_r8_kind
         endwhere
      enddo
    end subroutine radial3

    subroutine contract(n_unc, rad, contrs)
      ! Good!
      use matrix_module, only: matmult
      implicit none
      integer(i4_kind), intent(in) :: n_unc
      real(r8_kind), intent(inout) :: rad(:,:)
      real(r8_kind), intent(in)    :: contrs(:,:)
      ! *** end of interface ***

      integer(i4_kind) :: n_exps, n_contrs

      n_exps   = size(contrs,1)
      n_contrs = size(contrs,2)

      !
      ! Calculate contracted orbitals factors for all contractions and grid points
      !
      rad(:vl, n_unc+1:n_unc+n_contrs) = matmult(vl, n_contrs, n_exps, rad, contrs, "nn")

      !
      ! NOTE: occasionally "vl" is not equal to "size(rad, 1)" so that doing a
      !
      !         matmult(rad(:vl, :n_exps), contrs)
      !
      ! will cause creating a temporary array before passing the data to dgemm.
      !
      ! The in-place modification of "rad(:, :)" cannot be translated to a plain
      ! dgemm call that is realized as a subroutine.
      !
      ! FIXME: may be just use intrinsic MATMUL here?
      !
    end subroutine contract

  end subroutine radial

  subroutine orbital_calculate(grid_points,N_points,orbs_ob,grads, &
       sec_der,nuc_grads,nuc_sec_der,orb_xc,orb_ch,orbs_spinor_ob,spinor_grads, &
       nuc_3rd_der)
    !  Purpose: calculates orbitals, gradients, secound derivatives,
    !           charge and exchange fitfunctions (depending on which
    !           optional arguments are used) on a set of grid points
    use solhrules_module
    use xc_cntrl, only: whatis, xc_orbcalc_version, xc_nuc_sec_der_version
    use calc3c_switches
    use options_module, only : options_orbitals_in_memory, &
                               options_orbitals_on_file
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalprojection_module ! projection of indices to metaindex used in scf
    use orbitalstore_module, only: orbitalstore_initialized, orbitalstore_rw, &
         stored_orbs_ob, stored_grads, stored_orbs_spinor_ob, stored_spinor_grads
    use symmetry_data_module ! description of irreps
    use solid_harmonics_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in), optional :: grid_points(:,:)
    ! give this argument only at first call for given grid_points.
    ! Solid harmonics are evaluated when this argument is present and
    ! remain stored
    integer(i4_kind), intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    ! Old declaration returned to call "orbitalstore_rw"
    type(orbital_type), pointer, optional :: orbs_ob(:) ! intent(out)
    ! orbs_ob(N_irreps)
    type(spinor_type), pointer, optional :: orbs_spinor_ob(:) ! intent(out)
    ! orbs_ob(N_irreps)
    ! Old declaration returned to call "orbitalstore_rw"
    type(orbital_gradient_type), pointer, optional :: grads(:)
    ! intent(out) grad_ob(N_irreps)
    type(spinor_gradient_type), pointer, optional :: spinor_grads(:)
    ! intent(out) grad_ob(N_irreps)
    type(orbital_sec_der_type), intent(out), optional :: sec_der(:)
    ! intent(out) sec_der(N_irreps)
    type(orbital_nuclear_gradient_type), intent(out), optional :: nuc_grads(:,:)
    ! intent(out) nuc_grads(N_irreps,N_unique_atoms)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_sec_der(:,:)
    ! intent(out) nuc_sec_der(N_irreps,N_unique_atoms)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_3rd_der(:,:)
    type(orbital_type), intent(out), optional :: orb_ch ! intent(out)
    type(orbital_type), intent(out), optional :: orb_xc ! intent(out)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer             :: i_ua, i_ir, n_ir, i_ea, gridlength, i_l
    logical             :: last_atom
    ! loop indices
    type(unique_atom_type),             pointer :: ua
    type(solid_harmonics_array_type),   pointer :: sha
    logical                                     :: orbital_read
    FPP_TIMER_DECL(oc)
    FPP_TIMER_DECL(sh)
    FPP_TIMER_DECL(oco)
    FPP_TIMER_DECL(ini)
    !------------ Executable code ------------------------------------

    FPP_TIMER_START(oc)
    DLPRINT(9) 'om/orbital_calculate: entered'

    ASSERT(memory_allocated)

    if ( (present(grads) .or. present(sec_der)) ) then
        ASSERT(grads_allocated)
    endif

    if ( present(N_points) ) then
       ASSERT(N_points<=vec_len)
       gridlength = N_points
    else
       gridlength = vec_len
    endif

    ! initialising orbs and orbs_dim_indices
    n_ir = symmetry_data_n_irreps()
    orbital_read=.false.

    if ( present(orbs_ob) ) then
       if(orbitalstore_initialized) then
          call orbitalstore_rw(orbs_ob=orbs_ob)
          if(stored_orbs_ob) orbital_read=.true.
!!$       else
!!$          FPP_TIMER_START(ini)
!!$          do i_ir = 1,n_ir
!!$             ! FIXME: test init
!!$             orbs_ob(i_ir)%o = HUGE(1.0_r8_kind)
!!$          enddo
!!$          FPP_TIMER_STOP(ini)
          ! DO NOT INITIALIZE, TAKES A LOT OF TIME
          ! better fix the sub!
       end if
    end if

    if ( present(grads) ) then
       if(orbitalstore_initialized) then
          call orbitalstore_rw(grads=grads)
          if(stored_grads) orbital_read=.true.
!!$       else
!!$          FPP_TIMER_START(ini)
!!$          do i_ir = 1,n_ir
!!$             ! FIXME: test init
!!$             grads(i_ir)%o = HUGE(1.0_r8_kind)
!!$          enddo
!!$          FPP_TIMER_STOP(ini)
          ! DO NOT INITIALIZE, TAKES A LOT OF TIME
          ! better fix the sub!
       endif
    endif

    if ( present(orbs_spinor_ob) ) then
       !
       ! SPIN ORBIT
       !
       n_ir = symmetry_data_n_proj_irreps()

       if( orbitalstore_initialized )then
          call orbitalstore_rw(orbs_spinor_ob=orbs_spinor_ob)
          if(stored_orbs_spinor_ob) orbital_read=.true.
      endif
      ! DO NOT INITIALIZE, TAKES A LOT OF TIME
    endif ! present(orbs_spinor_ob)

    if ( present(spinor_grads) ) then
       !
       ! SPIN ORBIT
       !
       n_ir = symmetry_data_n_proj_irreps()
       if( orbitalstore_initialized ) then
          call orbitalstore_rw(spinor_grads=spinor_grads)
          if(stored_spinor_grads) orbital_read=.true.
      endif
      ! DO NOT INITIALIZE, TAKES A LOT OF TIME
    endif ! present(orbs_spinor_ob)

    if ( present(sec_der) ) then
       FPP_TIMER_START(ini)
       do i_ir = 1,n_ir
          sec_der(i_ir)%o = 0.0_r8_kind
       enddo
       FPP_TIMER_STOP(ini)
    endif

    if ( present(nuc_grads) ) then
       FPP_TIMER_START(ini)
       do i_ua=1,n_unique_atoms
          do i_ir = 1,n_ir
             nuc_grads(i_ua,i_ir)%o = 0.0_r8_kind
          enddo
       end do
       FPP_TIMER_STOP(ini)
    endif

    if ( present(nuc_sec_der) ) then
       FPP_TIMER_START(ini)
       do i_ua=1,n_unique_atoms
          do i_ir = 1,n_ir
             nuc_sec_der(i_ua,i_ir)%o = 0.0_r8_kind
          enddo
       end do
       FPP_TIMER_STOP(ini)
    endif

    if ( present(nuc_3rd_der) ) then
       FPP_TIMER_START(ini)
       do i_ua=1,n_unique_atoms
          do i_ir = 1,n_ir
             nuc_3rd_der(i_ua,i_ir)%o = 0.0_r8_kind
          enddo
       end do
       FPP_TIMER_STOP(ini)
    endif

    if ( present(orb_ch) ) then
       i_ir = get_totalsymmetric_irrep()
       FPP_TIMER_START(ini)
       orb_ch%o = 0.0_r8_kind
       FPP_TIMER_STOP(ini)
       if (do_glob_con_ch) call uncont_fit_fct_allocate("ch",orbs=.true.)
    endif

    if ( present(orb_xc) ) then
       i_ir = get_totalsymmetric_irrep()
       FPP_TIMER_START(ini)
       orb_xc%o = 0.0_r8_kind
       FPP_TIMER_STOP(ini)
       if (do_glob_con_xc) call uncont_fit_fct_allocate("xc",orbs=.true.)
    endif

    ! *** RETURN POINT *** for calculated orbitals etc.
    if(orbitalstore_initialized .and. orbital_read) RETURN
    ! *** RETURN POINT ***

    ! orbitals do have to be computed:
    ! loop over unique atoms
    i_ua_:do i_ua = 1, N_unique_atoms
       last_atom = i_ua==N_unique_atoms
       ua => unique_atoms(i_ua)
       sha => solid_harmonics(i_ua)

       FPP_TIMER_START(sh)
       if ( present(grid_points) ) then
          ! prepare solid harmonics for all equal atoms
          do i_ea = 1, ua%N_equal_atoms
             call solid_harmonics_calculate( sha%sh(i_ea), &
                  grid_points,ua%position(1:3,i_ea),N_points )
          enddo
       endif

       if ( present(grads) .or. present(sec_der).or. &
            present(nuc_sec_der).or.present(spinor_grads) ) then

          do i_ea = 1, ua%N_equal_atoms
             call solid_harmonics_calculate_grads( sha%sh(i_ea), &
                  sha%shg(i_ea), N_points )
          enddo
       endif

       if ( present(sec_der).or.present(nuc_sec_der) ) then
          do i_ea = 1, ua%N_equal_atoms
             call solid_harmonics_calc_sec_der( sha%shg(i_ea), &
                  sha%shs(i_ea), N_points )
          enddo
       endif

       if ( present(nuc_3rd_der) ) then
          do i_ea = 1, ua%N_equal_atoms
             call solid_harmonics_calc_3rd_der( sha%shs(i_ea), &
                  sha%sht(i_ea), N_points )
          enddo
       endif
       FPP_TIMER_STOP(sh)
       !
       ! Now calculate symmetrized basis on the grid:
       !
       FPP_TIMER_START(oco)
       if(.not.present(grads).and..not.present(spinor_grads))then
          !
          ! Gradients are not needed
          !
          if( present(orbs_ob) )then
             !
             ! STANDARD (NO SPIN-ORBIT) CASE
             !
             call calculate_orbitals(orbs_ob,orbitalprojection_ob,ua%lmax_ob, ua%l_ob)
             if(options_orbitals_on_file() .and. last_atom) call orbitalstore_rw(orbs_ob=orbs_ob)

          else if( present(orbs_spinor_ob) )then
             !
             ! SPIN-ORBIT CASE
             !
             call spinors_and_grads (ua%lmax_ob, ua%l_ob, orbitalprojection_spor_ob, &
                  orbs_spinor_ob)
          endif
          !
          ! see below the case when both are "present()"
          !
       endif

       if ( present(orbs_ob) .and. present(grads).and. .not.&
            present(nuc_grads) ) then
          if ( whatis(xc_orbcalc_version) .eq. 3 ) then

             ! Version 3
             call calculate_orbs_and_grads_v3(&
                  & gridlength, ua%n_equal_atoms,&
                  & orbs_ob,grads, orbitalprojection_ob, ua%lmax_ob, ua%l_ob)
          else
             ABORT('no more supported!')
          endif

          if(options_orbitals_on_file() .and. last_atom) call orbitalstore_rw(orbs_ob=orbs_ob,grads=grads)
       endif

       if ( present(orbs_spinor_ob) .and. present(spinor_grads).and. .not.&
            present(nuc_grads) ) then
          call spinors_and_grads (ua%lmax_ob, ua%l_ob, orbitalprojection_spor_ob, &
               orbs_spinor_ob, spinor_grads)
       endif

       if ( present(orbs_ob) .and. present(grads) .and. present(nuc_grads) )&
            then
          call calculate_orbs_and_nuc_grads( &
               orbs_ob,grads,nuc_grads,orbitalprojection_ob, &
               ua%lmax_ob, ua%l_ob)
       endif
       if ( .not. present(orbs_ob) .and. present(grads) ) then
          call calculate_gradients( &
               grads,orbitalprojection_ob, ua%lmax_ob, ua%l_ob)
       endif

       if ( present(sec_der) ) then
          call calculate_sec_der( &
               sec_der,orbitalprojection_ob, ua%lmax_ob, ua%l_ob)
       endif

       if ( present(nuc_sec_der) ) then
          if (whatis(xc_nuc_sec_der_version).eq.2) then
             call calculate_nuc_sec_der_v2( &
               nuc_sec_der,orbitalprojection_ob, ua%lmax_ob, ua%l_ob)
          else
             ABORT('no more supported!')
          endif
       endif

      if ( present(nuc_3rd_der) ) then
       call calculate_nuc_3rd_der_v2( &
               nuc_3rd_der,orbitalprojection_ob, ua%lmax_ob, ua%l_ob)
      endif

       if ( present(orb_ch) ) then
          if ( ua%N_glob_cons_ch .gt. 0 ) then
             call calculate_fitfunctions( &
                  orb_ch,orbitalprojection_ch, ua%lmax_ch, ua%l_ch, &
                  ua%r2_ch, &
                  ua%N_glob_cons_ch, ua%glob_con_ch, uncont_orbs_ch(i_ua), &
                  orbitalprojection_globcontr_ch(i_ua) )
          else
             call calculate_fitfunctions( &
                  orb_ch,orbitalprojection_ch, ua%lmax_ch, ua%l_ch, &
                  ua%r2_ch )
          endif
       endif

       if ( present(orb_xc) ) then
          if ( ua%N_glob_cons_xc .gt. 0 ) then
             call calculate_fitfunctions( &
                  orb_xc,orbitalprojection_xc, ua%lmax_xc, ua%l_xc, &
                  ua%r2_xc, &
                  ua%N_glob_cons_xc,  ua%glob_con_xc, uncont_orbs_xc(i_ua), &
                  orbitalprojection_globcontr_xc(i_ua) )
          else
             call calculate_fitfunctions( &
                  orb_xc,orbitalprojection_xc, ua%lmax_xc, ua%l_xc, &
                  ua%r2_xc )
          endif
       endif

       FPP_TIMER_STOP(oco)
    enddo i_ua_ ! unique atoms

    if(present(orbs_spinor_ob).AND.present(orbs_ob))then
       !
       ! for the case when both have to be calculated,
       ! calculate spinors using already evaluated orbitals:
       !
       call calc_spinors(N_points,orbs_ob,orbs_spinor_ob)
    endif

    if ( present(orb_ch) ) then
       if (do_glob_con_ch) call uncont_fit_fct_free("ch",orbs=.true.)
    endif
    if ( present(orb_xc) ) then
       if (do_glob_con_xc) call uncont_fit_fct_free("xc",orbs=.true.)
    endif

    DLPRINT(9) 'om/orbital_calculate: exit'
    FPP_TIMER_STOP(oc)

    FPP_TIMER_PRINT(sh)
    FPP_TIMER_PRINT(ini)
    FPP_TIMER_PRINT(oco)
    FPP_TIMER_PRINT(oc)
!    write(*,*) "<=== Orbital calculate"
  contains


    subroutine calculate_orbitals( &
         orbs,orbitalprojection,lmax,l_bas) !!$,rpa)
      !  Purpose: does calculation for a given unique atom
      !  for orbitals.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_type),intent(inout) :: orbs(:) ! orbs(N_irreps), intent(out)
!!$      type(orbital_type), pointer :: orbs(:)
      integer(i4_kind),intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
!!$      integer, pointer :: orbitalprojection(:,:,:)
      ! orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)  :: lmax
!!$      integer, intent(in)  :: lmax
      type(unique_atom_basis_type),intent(in) :: l_bas(0:) ! intent(in)
!!$      type(unique_atom_basis_type), pointer :: l_bas(:)
      target :: l_bas
      ! l_bas(lmax)
      !------------ Declaration of local variables -----------------
      integer(i4_kind) :: i_vec, i_ea, i_exp, &
           i_pa, i_fu, i_if, i_l, i_ir, &
           i_sa, i_sa_intermediate
      ! loop indices
      integer(i4_kind)             :: m
!!$      real(kind=r8_kind)  :: alpha, exp_orb, norm
      ! intermediates
      real(kind=r8_kind)  :: coef !!$, coef_renorm
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:)
      real(kind=r8_kind), pointer    :: so(:,:,:)
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      integer(i4_kind) :: n_unc, n_fcts
      integer(i4_kind) :: i_sa1
      !------------ Executable code --------------------------------


      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            r2 => sh%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea) )

         enddo! equal atoms


         ! symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.
         do i_ir = 1,symmetry_data_n_irreps()  ! irreps
            uap => ua%symadapt_partner(i_ir,i_l)
            so => orbs(i_ir)%o
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               i_sa =orbitalprojection(i_ir,i_l,i_ua)
               i_sa_intermediate = i_sa
               do i_if = 1, uap%N_independent_fcts ! independent fcts
                  uas => uap%symadapt(i_if,i_pa)

                  ! first nullify:
                  i_sa1 = i_sa ! make a copy, not to advance i_sa
                  do i_exp = 1, n_fcts
                     so(:gridlength,i_sa1,i_pa) = zero
                     i_sa1 = i_sa1 + 1
                  enddo

                  do i_fu = 1, uas%N_fcts ! contributing functions
                     i_sa = i_sa_intermediate
                     i_ea = uas%I_equal_atom(i_fu)
                     if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           do i_vec = 1, gridlength
                              so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                                   coef * prim_orb_fac(i_vec,i_exp,i_ea)
                           enddo
                           i_sa = i_sa + 1
                        enddo! uncontracted exponents
                     else ! i_l .gt. 0  take solid harmonics into account
                        sh => sha%sh(i_ea)
                        shm => sh%l(i_l)%m
                        m = uas%m(i_fu)
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           do i_vec = 1, gridlength
                              so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                                   shm(i_vec,m) * coef * &
                                   prim_orb_fac(i_vec,i_exp,i_ea)
                           enddo
                           i_sa = i_sa + 1
                        enddo!exponents
                     endif! i_l .eq. 0
                  enddo! contributing functions
                  i_sa_intermediate = i_sa
               enddo! independent fcts
            enddo! partners
         enddo! irreps
      enddo! l

    end subroutine calculate_orbitals

    subroutine calculate_orbs_and_grads_v3 (vl, n_ea, orbs, grads, orbitalprojection, lmax, l_bas)
      !
      ! Purpose: does calculation for a given unique atom for orbitals
      ! and   gradients.   performs   l-specific   contraction  before
      ! symmetry adaption.
      !
#ifdef FPP_F77_COPY_IN_OUT
      use f77_blas, only: dgemm
#endif
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind), intent(in)             :: vl, n_ea
      type(orbital_type), intent(inout)          :: orbs(:) ! orbs(N_irreps), intent(out)
      type(orbital_gradient_type), intent(inout) :: grads(:)  ! grads(N_irreps), intent(out)
      integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
      !orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)             :: lmax
      type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
      target :: l_bas
      ! l_bas(lmax)
      !------------ Declaration of local variables -----------------

      integer(i4_kind)               ::&
           & i_ea, i_pa, i_fu, i_c, i_if, i_l, i_ir,&
           & i_sa, m
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module

      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas

      integer(i4_kind)                      :: xyz
      integer(i4_kind)                      :: n_fcts
      integer(i4_kind)                      :: n_unc
      real(r8_kind), dimension(n_ea)        :: sh_fx_sa
      real(r8_kind), dimension(vl,0:3,n_ea) :: sh_dfdx_sa
      real(r8_kind), dimension(vl,1:3,n_ea) :: rvec
      real(r8_kind), dimension(vl)          :: help_vec

      FPP_TIMER_DECL(tt)
      FPP_TIMER_DECL(rd)
      FPP_TIMER_DECL(ag)
      FPP_TIMER_DECL(fn)
      !------------ Executable code --------------------------------

      FPP_TIMER_START(tt)
      DLPRINT(9) 'om/calculate_orbs_and_grads_v3: entered'

      do i_ea=1,n_ea
         rvec(:vl,1,i_ea) = sha%sh(i_ea)%L(1)%m(:vl,2)
         rvec(:vl,2,i_ea) = sha%sh(i_ea)%L(1)%m(:vl,3)
         rvec(:vl,3,i_ea) = sha%sh(i_ea)%L(1)%m(:vl,1)
      enddo

      ! loop over l
      do i_l = 0, lmax
         DLPRINT(9) 'om/calculate_orbs_and_grads_v3: L=',i_l
         uab => l_bas(i_l)

         DLPRINT(9) 'om/calculate_orbs_and_grads_v3: doing radial distributions ...'
         ! loop over equal atoms
         FPP_TIMER_START(rd)
         do i_ea = 1, n_ea !, ua%N_equal_atoms
            r2 => sha%sh(i_ea)%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(vl, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea) )

         enddo! i_ea
         FPP_TIMER_STOP(rd)
         DLPRINT(9) 'om/calculate_orbs_and_grads_v3: doing symmetry adaption ...'


         ! Symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! See header  of obitalprojection_module to  understand order
         ! of results in so calculated.
         do i_ir = 1, symmetry_data_n_irreps() ! irreps
            uap => ua % symadapt_partner(i_ir, i_l)
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               associate (so => orbs(i_ir) % o (:, :, i_pa), sog => grads(i_ir) % o(:, :, :, i_pa))
                 !
                 ! Here i_sa is an index of symmetry adapted orbitals,
                 ! a  cumulative index  for independent  functions and
                 ! radial shapes.  Incremented for every  radial shape
                 ! at the top of the loop, hence -1 here:
                 !
                 i_sa = orbitalprojection(i_ir, i_l, i_ua) - 1

                 do i_if = 1, uap % N_independent_fcts ! independent fcts
                    uas => uap % symadapt(i_if, i_pa)

                    if ( i_l .eq. 0 ) then
                       sh_fx_sa   = ZERO
                    else
                       sh_fx_sa   = ZERO
                       sh_dfdx_sa = ZERO
                    endif

                    FPP_TIMER_START(ag)
                    do i_fu = 1, uas%N_fcts ! contributing functions


                       ! contribution to Sym. Adapt. Orb.:
                       i_ea  = uas%I_equal_atom(i_fu)
                       m     = uas%m(i_fu)
                       !coef = uas%c(i_fu)

                       if ( i_l .eq. 0 ) then
                          ! solid harm. are 1.0. their gradients are  0.0:
                          sh_fx_sa(i_ea) = uas%c(i_fu)
                       else
                          ! compute the VALUE and the DERIVATVES (X,Y,Z) of the sym-adapt harmonics:
                          sh_fx_sa(i_ea) = ONE ! just to indicate that i_ea is needed!
                          sh_dfdx_sa(:vl,VALUE,i_ea) = sh_dfdx_sa(:vl,VALUE,i_ea)&
                               & + uas%c(i_fu)&
                               & * sha%sh(i_ea)%L(i_l)%m(:vl,m)

                          sh_dfdx_sa(:vl,DX:DZ,i_ea) = sh_dfdx_sa(:vl,DX:DZ,i_ea)&
                               & + uas%c(i_fu)&
                               & * sha%shg(i_ea)%L(i_l)%m(:vl,DX:DZ,m)
                       endif
                    enddo
                    FPP_TIMER_STOP(ag)
                    FPP_TIMER_START(fn)

                    if (i_l .eq. 0) then ! solid harmonics are 1.0, use sh_fx_sa:
                       do i_c = 1, n_fcts ! radial functions
                          i_sa = i_sa + 1

                          ! Sum over symmetry-equivalent centers:
                          so(:vl, i_sa) = 0.0
                          do i_ea = 1, n_ea
                             if (sh_fx_sa(i_ea) == 0.0) cycle
                             so(:vl, i_sa) = so(:vl, i_sa) &
                                  + sh_fx_sa(i_ea) * prim_orb_fac(:vl, i_c, i_ea)
                          enddo

                          ! Same for the gradients:
                          sog(:vl, 1:3, i_sa) = 0.0
                          do i_ea = 1, n_ea
                             if (sh_fx_sa(i_ea) == 0.0) cycle
                             help_vec(:vl) = sh_fx_sa(i_ea) * prim_grad_fac(:vl,i_c,i_ea)
                             do xyz = DX, DZ
                                sog(:vl, xyz, i_sa) = sog(:vl, xyz, i_sa) &
                                     + help_vec(:vl) * rvec(:vl, xyz, i_ea)
                             enddo
                          enddo
                       enddo
                    else
                       do i_c = 1, n_fcts ! radial functions
                          i_sa = i_sa + 1

                          ! Sum over symmetry-equivalent centers:
                          so(:vl, i_sa) = 0.0
                          do i_ea = 1, n_ea
                             if (sh_fx_sa(i_ea) .eq. zero) cycle
                             so(:vl, i_sa) = so(:vl, i_sa) + &
                                  sh_dfdx_sa(:vl, VALUE, i_ea) * prim_orb_fac(:vl, i_c, i_ea)
                          enddo

                          ! Same for the gradients:
                          sog(:vl, 1:3, i_sa) = 0.0
                          do i_ea = 1, n_ea
                             if (sh_fx_sa(i_ea) .eq. zero) cycle
                             help_vec(:vl) = sh_dfdx_sa(:vl, VALUE, i_ea) * prim_grad_fac(:vl, i_c, i_ea)
                             do xyz = DX, DZ
                                sog(:vl, xyz, i_sa) = sog(:vl, xyz, i_sa) &
                                     + rvec(:vl, xyz, i_ea) * help_vec(:vl) &
                                     + sh_dfdx_sa(:vl, xyz, i_ea) * prim_orb_fac(:vl, i_c, i_ea)
                             enddo
                          enddo
                       enddo
                    endif
                    FPP_TIMER_STOP(fn)
                 enddo          ! independent fcts
               end associate    ! so  => ..., sog => ...
            enddo    ! partners
         enddo       ! irreps
      enddo          ! l

      DLPRINT(9) 'om/calculate_orbs_and_grads_v3: exit'
      FPP_TIMER_STOP(tt)
      FPP_TIMER_PRINT(fn)
      FPP_TIMER_PRINT(ag)
      FPP_TIMER_PRINT(rd)
      FPP_TIMER_PRINT(tt)
    end subroutine calculate_orbs_and_grads_v3

    subroutine spinors_and_grads( lmax, l_bas, oproj, orbs, grads )
      !  Purpose: does calculation for a given unique atom
      !  for orbitals and gradients.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(i4_kind), intent(in)              :: lmax
      type(unique_atom_basis_type), intent(in)  :: l_bas(0:)      ! (0:lmax)
      target :: l_bas
      integer(i4_kind)          , intent(in)    :: oproj(:,0:,:)
      ! orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      type(spinor_type)         , intent(inout) :: orbs(:)   ! orbs(N_irreps)
      type(spinor_gradient_type), intent(inout) :: grads(:)  ! grads(N_irreps)
      optional :: grads
      !------------ Declaration of local variables -----------------
      integer(i4_kind)    :: i_ea, vl, i_ir, i_sa
      ! loop indices
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap

      integer(i4_kind) :: n_unc, n_fcts
      real(r8_kind)    :: beta
      FPP_TIMER_DECL(sog)

      !------------ Executable code --------------------------------

      FPP_TIMER_START(sog)
      vl=gridlength

      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         beta = zero ! for the first EA
         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts

            if(present(grads))then
               ! we do one EA at a time, so that only prim_orb_fac(:,:,1) is used
               call radial(vl, sha%sh(i_ea)%r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                    & prim_orb_fac(:,:,1), &
                    & prim_grad_fac(:,:,1) )
            else
               call radial(vl, sha%sh(i_ea)%r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                    & prim_orb_fac(:,:,1) )
            endif

            i_ir_:do i_ir = 1,symmetry_data_n_proj_irreps()  ! irreps
               uap => ua%symadapt_spor_partner(i_ir,i_l)

               i_sa =oproj(i_ir,i_l,i_ua)

               ! spinors values:
               call do_spinor(i_sa, n_fcts, prim_orb_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, &
                    vl, orbs(i_ir)%o(:,:,:,RE,UP), RE, uap%sa_spor_int(i_ea,UP,:,:), &
                    beta)

               call do_spinor(i_sa, n_fcts, prim_orb_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, &
                    vl, orbs(i_ir)%o(:,:,:,IM,UP), IM, uap%sa_spor_int(i_ea,UP,:,:), &
                    beta)

               call do_spinor(i_sa, n_fcts, prim_orb_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, &
                    vl, orbs(i_ir)%o(:,:,:,RE,DN), RE, uap%sa_spor_int(i_ea,DN,:,:), &
                    beta)

               call do_spinor(i_sa, n_fcts, prim_orb_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, &
                    vl, orbs(i_ir)%o(:,:,:,IM,DN), IM, uap%sa_spor_int(i_ea,DN,:,:), &
                    beta)

               if(.not.present(grads)) cycle ! irrep
               ! spinor gradients:
               call do_spinor_grad(i_sa, n_fcts, prim_orb_fac(:,:,1), prim_grad_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, sha%sh(i_ea)%L(1)%m, sha%shg(i_ea)%L(i_l)%m, &
                    vl, grads(i_ir)%o(:,:,:,:,RE,UP), RE, uap%sa_spor_int(i_ea,UP,:,:), &
                    beta)

               call do_spinor_grad(i_sa, n_fcts, prim_orb_fac(:,:,1), prim_grad_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, sha%sh(i_ea)%L(1)%m, sha%shg(i_ea)%L(i_l)%m, &
                    vl, grads(i_ir)%o(:,:,:,:,IM,UP), IM, uap%sa_spor_int(i_ea,UP,:,:), &
                    beta)

               call do_spinor_grad(i_sa, n_fcts, prim_orb_fac(:,:,1), prim_grad_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, sha%sh(i_ea)%L(1)%m, sha%shg(i_ea)%L(i_l)%m, &
                    vl, grads(i_ir)%o(:,:,:,:,RE,DN), RE, uap%sa_spor_int(i_ea,DN,:,:), &
                    beta)

               call do_spinor_grad(i_sa, n_fcts, prim_orb_fac(:,:,1), prim_grad_fac(:,:,1), &
                    i_l, sha%sh(i_ea)%L(i_l)%m, sha%sh(i_ea)%L(1)%m, sha%shg(i_ea)%L(i_l)%m, &
                    vl, grads(i_ir)%o(:,:,:,:,IM,DN), IM, uap%sa_spor_int(i_ea,DN,:,:), &
                    beta)

            enddo i_ir_! irreps
            ! for the next UA, add up:
            beta = one
         enddo ! EA
      enddo! l
      FPP_TIMER_STOP(sog)
      FPP_TIMER_PRINT(sog)
    end subroutine spinors_and_grads

    subroutine calculate_gradients( &
         grads,orbitalprojection,lmax,l_bas ) !!$,rpa)
      !  Purpose: does calculation for a given unique atom
      !  for gradients.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_gradient_type), intent(inout) :: grads(:)  ! grads(N_irreps), intent(out)
!!$      type(orbital_gradient_type), pointer :: grads(:)  ! grads(N_irreps), intent(out)
      integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
!!$      integer, pointer ::orbitalprojection(:,:,:) ! intent(in)
      !orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)  :: lmax
!!$      integer, intent(in)  :: lmax
      type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
!!$      type(unique_atom_basis_type), pointer :: l_bas(:) ! intent(in)
      target :: l_bas
      ! l_bas(lmax)
      !------------ Declaration of local variables -----------------
      integer(i4_kind)    :: i_ea, i_exp, &
           i_pa, i_fu, i_if, i_l, i_ir, &
           i_sa, i_sa_intermediate, vl
      ! loop indices
      integer(i4_kind)    :: m
!!$      real(kind=r8_kind)  :: alpha, exp_orb, norm, grad_norm
      real(kind=r8_kind)  :: sym_fac(gridlength)
      ! intermediates
      real(kind=r8_kind)  :: coef !!$, coef_renorm
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:), shm_1(:,:), shm_grad(:,:,:), &
           sog(:,:,:,:), x(:), y(:), z(:)
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      type(solid_harmonics_grads_type),     pointer :: shg
      integer(i4_kind) :: n_unc, n_fcts
      integer(i4_kind) :: i_sa1
      !------------ Executable code --------------------------------


      vl=gridlength
      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            shg => sha%shg(i_ea)
            r2 => sh%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea) )

         enddo! equal atoms


         ! symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.
         do i_ir = 1,symmetry_data_n_irreps()  ! irreps
            uap => ua%symadapt_partner(i_ir,i_l)
            sog => grads(i_ir)%o
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               i_sa =orbitalprojection(i_ir,i_l,i_ua)
               i_sa_intermediate = i_sa
               do i_if = 1, uap%N_independent_fcts ! independent fcts
                  uas => uap%symadapt(i_if,i_pa)

                  ! first nullify:
                  i_sa1 = i_sa ! make a copy, not to advance i_sa
                  do i_exp = 1, n_fcts
                     sog(:vl,1:3,i_sa1,i_pa) = zero
                     i_sa1 = i_sa1 + 1
                  enddo

                  do i_fu = 1, uas%N_fcts ! contributing functions
                     i_sa = i_sa_intermediate
                     i_ea = uas%I_equal_atom(i_fu)
                     sh => sha%sh(i_ea)
                     shm_1 => sh%l(1)%m
                     x=>shm_1(1:vl,2)
                     y=>shm_1(1:vl,3)
                     z=>shm_1(1:vl,1)
                     if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                        ! gradients of solid harmonics are 0.0
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) + &
                                coef * prim_grad_fac(1:vl,i_exp,i_ea) * x
                           sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) + &
                                coef * prim_grad_fac(1:vl,i_exp,i_ea) * y
                           sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) + &
                                coef * prim_grad_fac(1:vl,i_exp,i_ea) * z
                           i_sa = i_sa + 1
                        enddo! uncontracted exponents
                     else ! i_l .gt. 0  take solid harmonics into account
                        shm => sh%l(i_l)%m
                        sh => sha%sh(i_ea)
                        shg => sha%shg(i_ea)
                        shm_grad => shg%l(i_l)%m
                        m = uas%m(i_fu)
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           sym_fac(1:vl) = shm(1:vl,m) * coef
                           sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) &
                                + sym_fac * x * prim_grad_fac(1:vl,i_exp,i_ea) &
                                + coef * shm_grad(1:vl,1,m) * prim_orb_fac(1:vl,i_exp,i_ea)
                           sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) &
                                + sym_fac * y * prim_grad_fac(1:vl,i_exp,i_ea) &
                                + coef * shm_grad(1:vl,2,m) * prim_orb_fac(1:vl,i_exp,i_ea)
                           sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) &
                                + sym_fac * z * prim_grad_fac(1:vl,i_exp,i_ea) &
                                + coef * shm_grad(1:vl,3,m) * prim_orb_fac(1:vl,i_exp,i_ea)
                           i_sa = i_sa + 1
                        enddo!exponents
                     endif! i_l .eq. 0
                  enddo! contributing functions
                  i_sa_intermediate = i_sa
               enddo! independent fcts
            enddo! partners
         enddo! irreps
      enddo! l

    end subroutine calculate_gradients


      subroutine calculate_orbs_and_nuc_grads( &
           orbs,grads,nuc_grads,orbitalprojection,lmax,l_bas) !!$,rpa)
        !  Purpose: does calculation for a given unique atom
        !  for orbitals and gradients.
        !  performs l-specific contraction before symmetry adaption
        implicit none
        !------------ Declaration of formal parameters ---------------
        type(orbital_type), intent(inout) :: orbs(:) ! orbs(N_irreps), intent(out)
!!$        type(orbital_type), pointer :: orbs(:) ! orbs(N_irreps), intent(out)
        type(orbital_gradient_type), intent(inout) :: grads(:)  ! grads(N_irreps), intent(out)
!!$        type(orbital_gradient_type), pointer :: grads(:)  ! grads(N_irreps), intent(out)
        type(orbital_nuclear_gradient_type), intent(inout) :: nuc_grads(:,:)
!!$        type(orbital_nuclear_gradient_type), pointer :: nuc_grads(:,:)
             ! grads(n_unique_atoms,N_irreps), intent(out)
        integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
!!$        integer, pointer ::orbitalprojection(:,:,:) ! intent(in)
        !orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
        integer(i4_kind), intent(in)  :: lmax
!!$        integer, intent(in)  :: lmax
        type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
!!$        type(unique_atom_basis_type), pointer :: l_bas(:) ! intent(in)
        target :: l_bas
           ! l_bas(lmax)
        !------------ Declaration of local variables -----------------
        integer(i4_kind)    :: i_ea, i_exp, &
                               i_pa, i_fu, i_if, i_l, i_ir, &
                               i_sa, i_sa_intermediate, i_orbs, i_orbs_intermediate
                               ! loop indices
        integer(i4_kind)    :: m,vl
!!$        real(kind=r8_kind)  :: alpha, exp_orb, norm, grad_norm
        real(kind=r8_kind)  :: sym_fac(gridlength), help_vec(gridlength), &
             & help_vec2(gridlength), help_vec3(gridlength)
           ! intermediates
        real(kind=r8_kind)  :: coef !!$, coef_renorm
           ! contraction coefficient
        real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
           ! square of distance center to point calculated by
           ! solid_harmonics_module
        real(kind=r8_kind), pointer    :: shm(:,:), shm_1(:,:), shm_grad(:,:,:)
        real(kind=r8_kind), pointer    :: so(:,:,:)
        real(kind=r8_kind), pointer    :: sog(:,:,:,:)
        real(kind=r8_kind), pointer    :: song(:,:,:,:,:)
        type(unique_atom_basis_type),         pointer :: uab
        type(unique_atom_partner_type),       pointer :: uap
        type(unique_atom_symadapt_type),      pointer :: uas
        type(solid_harmonics_type),           pointer :: sh
        type(solid_harmonics_grads_type),     pointer :: shg
        integer(i4_kind) :: n_unc, n_fcts
        integer(i4_kind) :: i_sa1
        !------------ Executable code --------------------------------

        vl=gridlength
        ! loop over l
        do i_l = 0, lmax
           uab => l_bas(i_l)

           ! loop over equal atoms
           do i_ea = 1, ua%N_equal_atoms
              sh => sha%sh(i_ea)
              shg => sha%shg(i_ea)
              r2 => sh%r2

              n_unc  = uab%N_uncontracted_fcts
              n_fcts = n_unc + uab%N_contracted_fcts
              call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea) )

           enddo! equal atoms


           ! symmetry adaption:
           ! for all irreps
           !         partners,
           !         independent fcts,
           !         uncontracted exponents and contractions,
           !         gridpoints:
           !   sum over equal atoms and m
           !
           ! see haeder of obitalprojection_module to understand order of
           ! results in so calculated.
           do i_ir = 1,symmetry_data_n_irreps()  ! irreps
              uap => ua%symadapt_partner(i_ir,i_l)
              so => orbs(i_ir)%o
              sog => grads(i_ir)%o
              song => nuc_grads(i_ua,i_ir)%o
              do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
                 i_sa =orbitalprojection(i_ir,i_l,i_ua)
                 i_sa_intermediate = i_sa
                 i_orbs=orbitalprojection(i_ir,i_l,i_ua)+1-orbitalprojection(i_ir,0,i_ua)
                 i_orbs_intermediate=i_orbs
                 do i_if = 1, uap%N_independent_fcts ! independent fcts
                    uas => uap%symadapt(i_if,i_pa)

                    ! first nullify:
                    i_sa1 = i_sa ! make a copy, not to advance i_sa
                    do i_exp = 1, n_fcts
                       so (:vl    ,i_sa1,i_pa) = zero
                       sog(:vl,1:3,i_sa1,i_pa) = zero
                       i_sa1 = i_sa1 + 1
                    enddo
                    ! FIXME: you still need to initialize "song"

                    do i_fu = 1, uas%N_fcts ! contributing functions
                       i_sa = i_sa_intermediate
                       i_orbs = i_orbs_intermediate
                       i_ea = uas%I_equal_atom(i_fu)
                       sh => sha%sh(i_ea)
                       shm_1 => sh%l(1)%m
                       if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                                              ! gradients of solid harmonics are 0.0
                          coef = uas%c(i_fu)
                          do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                             help_vec(1:vl) = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                             so(1:vl,i_sa,i_pa) = so(1:vl,i_sa,i_pa) + &
                                  coef * prim_orb_fac(1:vl,i_exp,i_ea)
                             help_vec2(1:vl)=help_vec(1:vl) * shm_1(1:vl,2)
                             song(1:vl,1,i_ea,i_orbs,i_pa)=song(1:vl,1,i_ea,i_orbs,i_pa)+&
                                     help_vec2(1:vl)
                             sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) + &
                                  help_vec2(1:vl)
                             help_vec2(1:vl)=help_vec(1:vl) * shm_1(1:vl,3)
                             song(1:vl,2,i_ea,i_orbs,i_pa)=song(1:vl,2,i_ea,i_orbs,i_pa) +&
                                  help_vec2(1:vl)
                             sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) + &
                                  help_vec2(1:vl)
                             help_vec2(1:vl)=help_vec(1:vl) * shm_1(1:vl,1)
                             song(1:vl,3,i_ea,i_orbs,i_pa)=song(1:vl,3,i_ea,i_orbs,i_pa) +&
                                  help_vec2(1:vl)
                             sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) + &
                                  help_vec2(1:vl)
                             i_sa = i_sa + 1
                             i_orbs=i_orbs + 1
                          enddo! uncontracted exponents
                       else ! i_l .gt. 0  take solid harmonics into account
                          shm => sh%l(i_l)%m
                          sh => sha%sh(i_ea)
                          shg => sha%shg(i_ea)
                          shm_grad => shg%l(i_l)%m
                          m = uas%m(i_fu)
                          coef = uas%c(i_fu)
                          do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                             sym_fac(1:vl) = shm(1:vl,m) * coef
                             help_vec(1:vl) = sym_fac(1:vl)*prim_grad_fac(1:vl,i_exp,i_ea)
                             help_vec2 = coef * prim_orb_fac(1:vl,i_exp,i_ea)
                             so(1:vl,i_sa,i_pa) = so(1:vl,i_sa,i_pa) &
                                  + sym_fac * prim_orb_fac(1:vl,i_exp,i_ea)
                             help_vec3(1:vl)=  shm_1(1:vl,2) * help_vec(1:vl) &
                                  + shm_grad(1:vl,1,m) * help_vec2(1:vl)
                             sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) &
                                  + help_vec3(1:vl)
                             song(1:vl,1,i_ea,i_orbs,i_pa) = song(1:vl,1,i_ea,i_orbs,i_pa) &
                                  + help_vec3(1:vl)
                             help_vec3(1:vl)=  shm_1(1:vl,3) * help_vec(1:vl) &
                                  + shm_grad(1:vl,2,m) * help_vec2(1:vl)
                             sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) &
                                  +  help_vec3(1:vl)
                             song(1:vl,2,i_ea,i_orbs,i_pa) = song(1:vl,2,i_ea,i_orbs,i_pa) &
                                  +   help_vec3(1:vl)
                             help_vec3(1:vl)=  shm_1(1:vl,1) * help_vec(1:vl) &
                                  + shm_grad(1:vl,3,m) * help_vec2(1:vl)
                             sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) &
                                  + help_vec3(1:vl)
                             song(1:vl,3,i_ea,i_orbs,i_pa) = song(1:vl,3,i_ea,i_orbs,i_pa) &
                                  +  help_vec3(1:vl)
                             i_sa = i_sa + 1
                             i_orbs = i_orbs + 1
                          enddo!exponents
                       endif ! i_l .eq. 0
                    enddo ! contributing functions
                    i_sa_intermediate = i_sa
                    i_orbs_intermediate = i_orbs

                 enddo ! independent fcts
              enddo ! partners
           enddo ! irreps
        enddo! l

      end subroutine calculate_orbs_and_nuc_grads

    subroutine calculate_sec_der( &
         sec_der,orbitalprojection,lmax,l_bas) !!$,rpa)
      !  Purpose: does calculation for a given unique atom
      !  for secound derivatives.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_sec_der_type), intent(inout) :: sec_der(:)  ! sec_der(N_irreps), intent(out)
!!$      type(orbital_sec_der_type), pointer :: sec_der(:)  ! sec_der(N_irreps), intent(out)
      integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
!!$      integer, pointer ::orbitalprojection(:,:,:) ! intent(in)
      !orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)  :: lmax
!!$      integer, intent(in)  :: lmax
      type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
!!$      type(unique_atom_basis_type), pointer :: l_bas(:) ! intent(in)
      target :: l_bas
      ! l_bas(lmax)
      !------------ Declaration of local variables -----------------
      integer(i4_kind)    :: i_ea, i_exp, i_pa, i_fu, i_if, i_l, i_ir, &
           i_sa, i_sa_intermediate, vl, m
      ! loop indices
!!$      real(kind=r8_kind)  :: alpha, exp_orb, norm, grad_norm
      real(kind=r8_kind)  ::  sym_fac(gridlength), &
           help_vec(gridlength), help_vec2(gridlength), help_vec3(gridlength), &
           help_vec4(gridlength)
      ! intermediates
      real(kind=r8_kind)  :: coef !!$, coef_renorm
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:), shm_1(:,:), shm_grad(:,:,:), &
           x(:), y(:), z(:)
      real(kind=r8_kind), pointer    :: sosd(:,:,:,:)
      real(kind=r8_kind), pointer    :: shm_sd(:,:,:)
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      type(solid_harmonics_grads_type),     pointer :: shg
      type(solid_harmonics_sec_der_type),   pointer :: shsd
      integer(i4_kind) :: n_unc, n_fcts
      !------------ Executable code --------------------------------


      vl=gridlength
      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            shg => sha%shg(i_ea)
            shsd => sha%shs(i_ea)
            r2 => sh%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea), &
                 & prim_sd_fac(:,:,i_ea) )

         enddo! equal atoms


         ! symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.
         do i_ir = 1,symmetry_data_n_irreps()  ! irreps
            uap => ua%symadapt_partner(i_ir,i_l)
            sosd => sec_der(i_ir)%o
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               i_sa =orbitalprojection(i_ir,i_l,i_ua)
               i_sa_intermediate = i_sa
               do i_if = 1, uap%N_independent_fcts ! independent fcts
                  uas => uap%symadapt(i_if,i_pa)
                  do i_fu = 1, uas%N_fcts ! contributing functions
                     i_sa = i_sa_intermediate
                     i_ea = uas%I_equal_atom(i_fu)
                     sh => sha%sh(i_ea)
                     shm_1 => sh%l(1)%m
                     x=>shm_1(1:vl,2)
                     y=>shm_1(1:vl,3)
                     z=>shm_1(1:vl,1)
                     if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                        ! gradients of solid harmonics are 0.0
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           help_vec = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_sd_fac(1:vl,i_exp,i_ea)
                           !sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) + &
                           !     help_vec(1:vl) * x
                           sosd(1:vl,XX,i_sa,i_pa) = sosd(1:vl,XX,i_sa,i_pa) + &
                                help_vec + help_vec2 * x * x
                           sosd(1:vl,XY,i_sa,i_pa) = sosd(1:vl,XY,i_sa,i_pa) + &
                                help_vec2 * x * y
                           sosd(1:vl,XZ,i_sa,i_pa) = sosd(1:vl,XZ,i_sa,i_pa) + &
                                help_vec2 * x * z
                           !sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) + &
                           !     help_vec(1:vl) * y
                           sosd(1:vl,YY,i_sa,i_pa) = sosd(1:vl,YY,i_sa,i_pa) + &
                                help_vec + help_vec2 * y * y
                           sosd(1:vl,YZ,i_sa,i_pa) = sosd(1:vl,YZ,i_sa,i_pa) + &
                                help_vec2 * y * z
                           !sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) + &
                           !     help_vec(1:vl) * z
                           sosd(1:vl,ZZ,i_sa,i_pa) = sosd(1:vl,ZZ,i_sa,i_pa) + &
                                help_vec + help_vec2 * z * z
                           i_sa = i_sa + 1
                        enddo! uncontracted exponents
                     else ! i_l .gt. 0  take solid harmonics into account
                        shm => sh%l(i_l)%m
                        sh => sha%sh(i_ea)
                        shg => sha%shg(i_ea)
                        shm_grad => shg%l(i_l)%m
                        shm_sd => shsd%l(i_l)%m
                        m = uas%m(i_fu)
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           sym_fac(1:vl) = shm(1:vl,m) * coef
                           help_vec =  coef * shm(1:vl,m) * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_orb_fac(1:vl,i_exp,i_ea)
                           help_vec3 = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec4 = coef * shm(1:vl,m) * prim_sd_fac(1:vl,i_exp,i_ea)

                           !sog(1:vl,1,i_sa,i_pa) = sog(1:vl,1,i_sa,i_pa) &
                           !     + x * help_vec(1:vl) &
                           !     + shm_grad(1:vl,1,m) * help_vec2(1:vl)
                           sosd(1:vl,XX,i_sa,i_pa) = sosd(1:vl,XX,i_sa,i_pa) + &
                                help_vec + &
                                x * shm_grad(1:vl,1,m) * help_vec3 * 2.0_r8_kind + &
                                x * x * help_vec4 + &
                                shm_sd(1:vl,XX,m) * help_vec2(1:vl)
                           sosd(1:vl,XY,i_sa,i_pa) = sosd(1:vl,XY,i_sa,i_pa) + &
                                x * shm_grad(1:vl,2,m) * help_vec3 + &
                                x * y * help_vec4 + &
                                shm_sd(1:vl,XY,m) * help_vec2(1:vl) + &
                                shm_grad(1:vl,1,m) * y * help_vec3
                           sosd(1:vl,XZ,i_sa,i_pa) = sosd(1:vl,XZ,i_sa,i_pa) + &
                                x * shm_grad(1:vl,3,m) * help_vec3 + &
                                x * z * help_vec4 + &
                                shm_sd(1:vl,XZ,m) * help_vec2(1:vl) + &
                                shm_grad(1:vl,1,m) * z * help_vec3

                           !sog(1:vl,2,i_sa,i_pa) = sog(1:vl,2,i_sa,i_pa) &
                           !     + y *  help_vec(1:vl) &
                           !     + shm_grad(1:vl,2,m) *  help_vec2(1:vl)
                           sosd(1:vl,YY,i_sa,i_pa) = sosd(1:vl,YY,i_sa,i_pa) + &
                                help_vec + &
                                y * shm_grad(1:vl,2,m) * help_vec3 * 2.0_r8_kind + &
                                y * y * help_vec4 + &
                                shm_sd(1:vl,YY,m) * help_vec2(1:vl)
                           sosd(1:vl,YZ,i_sa,i_pa) = sosd(1:vl,YZ,i_sa,i_pa) + &
                                y * shm_grad(1:vl,3,m) * help_vec3 + &
                                y * z * help_vec4 + &
                                shm_sd(1:vl,YZ,m) * help_vec2(1:vl) + &
                                shm_grad(1:vl,2,m) * z * help_vec3

                           !sog(1:vl,3,i_sa,i_pa) = sog(1:vl,3,i_sa,i_pa) &
                           !     + z * help_vec(1:vl) &
                           !     + shm_grad(1:vl,3,m) * help_vec2(1:vl)
                           sosd(1:vl,ZZ,i_sa,i_pa) = sosd(1:vl,ZZ,i_sa,i_pa) + &
                                help_vec + &
                                z * shm_grad(1:vl,3,m) * help_vec3 * 2.0_r8_kind + &
                                z * z * help_vec4 + &
                                shm_sd(1:vl,ZZ,m) * help_vec2(1:vl)

                           i_sa = i_sa + 1
                        enddo!exponents
                     endif! i_l .eq. 0
                  enddo! contributing functions
!!$                  ! copy symmetric secound derivatives
!!$                  i_sa = i_sa_intermediate
!!$                  do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
!!$                     sosd(1:vl,2,1,i_sa,i_pa) = sosd(1:vl,1,2,i_sa,i_pa)
!!$                     sosd(1:vl,3,1,i_sa,i_pa) = sosd(1:vl,1,3,i_sa,i_pa)
!!$                     sosd(1:vl,3,2,i_sa,i_pa) = sosd(1:vl,2,3,i_sa,i_pa)
!!$                     i_sa = i_sa + 1
!!$                  enddo! uncontracted exponents
                  i_sa_intermediate = i_sa
               enddo! independent fcts
            enddo! partners
         enddo! irreps
      enddo! l

    end subroutine calculate_sec_der

    subroutine calculate_nuc_sec_der_v2 (nuc_sec_der, orbitalprojection, &
         lmax, l_bas)
      !  Purpose: does calculation for a given unique atom
      !  for secound derivatives.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_nuclear_sec_der_type), intent(inout) :: nuc_sec_der(:,:)
      ! sec_der(N_unique_atoms,N_irreps)
      integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
      ! orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)  :: lmax
      type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
      target :: l_bas           ! l_bas(lmax)
      ! *** end of interface ***

      integer(i4_kind)    :: i_ea, i_exp, i_pa, i_fu, i_if, i_l, i_ir, &
           i_sa, i_sa_intermediate, i_orbs,i_orbs_intermediate, vl, m
      ! loop indices
!!$      real(kind=r8_kind)  :: alpha, exp_orb, norm, grad_norm
      real(kind=r8_kind)  :: sym_fac(gridlength), &
           help_vec(gridlength), help_vec2(gridlength), help_vec3(gridlength), &
           help_vec4(gridlength)
      ! intermediates
      real(kind=r8_kind)  :: coef !!$, coef_renorm
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:), shm_1(:,:), shm_grad(:,:,:), &
           sosd(:,:,:,:,:)
      real(kind=r8_kind), pointer    :: shm_sd(:,:,:)
!!$      real(kind=r8_kind), pointer    :: x(:), y(:), z(:)
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      type(solid_harmonics_grads_type),     pointer :: shg
      type(solid_harmonics_sec_der_type),   pointer :: shsd

      real(r8_kind)    :: xyz(gridlength,3)
      integer(i4_kind) :: D1,D2,D12
      integer(i4_kind) :: n_unc, n_fcts
      !------------ Executable code --------------------------------

      vl=gridlength
      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            shg => sha%shg(i_ea)
            shsd => sha%shs(i_ea)
            r2 => sh%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea), &
                 & prim_sd_fac(:,:,i_ea) )

         enddo! equal atoms


         ! symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.
         do i_ir = 1,symmetry_data_n_irreps()  ! irreps
            uap => ua%symadapt_partner(i_ir,i_l)
            sosd => nuc_sec_der(i_ua,i_ir)%o
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               i_sa =orbitalprojection(i_ir,i_l,i_ua)
               i_sa_intermediate = i_sa
               i_orbs=orbitalprojection(i_ir,i_l,i_ua)+1-orbitalprojection(i_ir,0,i_ua)
               i_orbs_intermediate=i_orbs
               do i_if = 1, uap%N_independent_fcts ! independent fcts
                  uas => uap%symadapt(i_if,i_pa)
                  do i_fu = 1, uas%N_fcts ! contributing functions
                     i_sa = i_sa_intermediate
                     i_orbs = i_orbs_intermediate
                     i_ea = uas%I_equal_atom(i_fu)
                     sh => sha%sh(i_ea)
                     shm_1 => sh%l(1)%m

                     xyz(:,1) = shm_1(1:vl,2)
                     xyz(:,2) = shm_1(1:vl,3)
                     xyz(:,3) = shm_1(1:vl,1)
                     if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                        ! gradients of solid harmonics are 0.0
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           help_vec = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_sd_fac(1:vl,i_exp,i_ea)

                           D12 = 0
                           do D1=1,3
                              do D2=D1,3
                                 D12 = D12 + 1
                                 sosd(1:vl,D12,i_ea,i_orbs,i_pa) = sosd(1:vl,D12,i_ea,i_orbs,i_pa) &
                                      & + help_vec2 * xyz(:,D1) * xyz(:,D2) &
                                      & + merge(help_vec,zero,D1.eq.D2) ! * Delta(D1,D2)
                              enddo
                           enddo
                           i_sa = i_sa + 1
                           i_orbs=i_orbs + 1
                        enddo! uncontracted exponents
                     else ! i_l .gt. 0  take solid harmonics into account
                        shm => sh%l(i_l)%m
                        sh => sha%sh(i_ea)
                        shg => sha%shg(i_ea)
                        shsd => sha%shs(i_ea)
                        shm_grad => shg%l(i_l)%m
                        shm_sd => shsd%l(i_l)%m
                        m = uas%m(i_fu)
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           sym_fac(1:vl) = shm(1:vl,m) * coef
                           help_vec =  coef * shm(1:vl,m) * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_orb_fac(1:vl,i_exp,i_ea)
                           help_vec3 = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec4 = coef * shm(1:vl,m) * prim_sd_fac(1:vl,i_exp,i_ea)

                           D12 = 0
                           do D1=1,3
                              do D2=D1,3
                                 D12 = D12 + 1
                                 sosd(1:vl,D12,i_ea,i_orbs,i_pa) = sosd(1:vl,D12,i_ea,i_orbs,i_pa) &
                                      & +    xyz(:,D1) * xyz(:,D2)             * help_vec4 &
                                      & + (  xyz(:,D1) * shm_grad(1:vl,D2,m) &
                                      &    + xyz(:,D2) * shm_grad(1:vl,D1,m) ) * help_vec3 &
                                      & +    shm_sd(1:vl,D12,m)                * help_vec2 &
                                      & + merge(help_vec,zero,D1.eq.D2) ! * Delta(D1,D2)
!!!                                 + prim_sd_fac(1:vl,i_exp,i_ea)*coef
                              enddo
                           enddo
                           i_sa = i_sa + 1
                           i_orbs=i_orbs+1
                        enddo!exponents
                     endif! i_l .eq. 0
                  enddo! contributing functions
                  ! copy symmetric secound derivatives
                  i_sa_intermediate = i_sa
                  i_orbs_intermediate = i_orbs
               enddo! independent fcts
            enddo! partners
         enddo! irreps
      enddo! l

    end subroutine calculate_nuc_sec_der_v2

    subroutine calculate_nuc_3rd_der_v2( &
         nuc_3rd_der,orbitalprojection,lmax,l_bas)
      !  Purpose: does calculation for a given unique atom
      !  for secound derivatives.
      !  performs l-specific contraction before symmetry adaption
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_nuclear_sec_der_type), intent(inout) :: nuc_3rd_der(:,:)
      ! sec_der(N_unique_atoms,N_irreps)
      integer(i4_kind), intent(in) :: orbitalprojection(:,0:,:) ! intent(in)
      ! orbitalprojection(N_irreps,N_max_l,N_unique_atoms)
      integer(i4_kind), intent(in)  :: lmax
      type(unique_atom_basis_type), intent(in) :: l_bas(0:) ! intent(in)
      target :: l_bas
      ! l_bas(lmax)
      !------------ Declaration of local variables -----------------
      integer(i4_kind)    :: i_ea, i_exp, i_pa, i_fu, i_if, i_l, i_ir, &
           i_sa, i_sa_intermediate, i_orbs,i_orbs_intermediate, vl, m
      ! loop indices
!!$      real(kind=r8_kind)  :: alpha, exp_orb, norm, grad_norm
      real(kind=r8_kind)  :: sym_fac(gridlength), &
           help_vec(gridlength), help_vec2(gridlength), help_vec3(gridlength), &
           help_vec4(gridlength),DD(gridlength,3,3),help_vec5(gridlength), &
           help_vec6(gridlength)
      ! intermediates
      real(kind=r8_kind)  :: coef !!$, coef_renorm
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:), shm_1(:,:), shm_grad(:,:,:), &
           sotd(:,:,:,:,:)
      real(kind=r8_kind), pointer    :: shm_sd(:,:,:)
      real(kind=r8_kind), pointer    :: shm_td(:,:,:)
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      type(solid_harmonics_grads_type),     pointer :: shg
      type(solid_harmonics_sec_der_type),   pointer :: shsd
      type(solid_harmonics_sec_der_type),   pointer :: shtd

      real(r8_kind)    :: xyz(gridlength,3)
      integer(i4_kind) :: D1,D2,D12,D3
      integer(i4_kind) :: n_unc, n_fcts
      !------------ Executable code --------------------------------

      vl=gridlength
      ! loop over l
      do i_l = 0, lmax
         uab => l_bas(i_l)

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            shg => sha%shg(i_ea)
            shsd => sha%shs(i_ea)
            shtd => sha%sht(i_ea)
            r2 => sh%r2

            n_unc  = uab%N_uncontracted_fcts
            n_fcts = n_unc + uab%N_contracted_fcts
            call radial(gridlength, r2, uab%exponents, uab%norms, n_unc, uab%contractions, &
                 & prim_orb_fac(:,:,i_ea), &
                 & prim_grad_fac(:,:,i_ea), &
                 & prim_sd_fac(:,:,i_ea), &
                 & prim_td_fac(:,:,i_ea) )

         enddo! equal atoms


         ! symmetry adaption:
         ! for all irreps
         !         partners,
         !         independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.
         do i_ir = 1,symmetry_data_n_irreps()  ! irreps
            uap => ua%symadapt_partner(i_ir,i_l)
            sotd => nuc_3rd_der(i_ua,i_ir)%o
            do i_pa = 1, symmetry_data_n_partners(i_ir) ! partners
               i_sa =orbitalprojection(i_ir,i_l,i_ua)
               i_sa_intermediate = i_sa
               i_orbs=orbitalprojection(i_ir,i_l,i_ua)+1-orbitalprojection(i_ir,0,i_ua)
               i_orbs_intermediate=i_orbs
               do i_if = 1, uap%N_independent_fcts ! independent fcts
                  uas => uap%symadapt(i_if,i_pa)
                  do i_fu = 1, uas%N_fcts ! contributing functions
                     i_sa = i_sa_intermediate
                     i_orbs = i_orbs_intermediate
                     i_ea = uas%I_equal_atom(i_fu)

                     sh => sha%sh(i_ea)
                     shm_1 => sh%l(1)%m

                     xyz(:,1) = shm_1(1:vl,2)
                     xyz(:,2) = shm_1(1:vl,3)
                     xyz(:,3) = shm_1(1:vl,1)
                     if ( i_l .eq. 0 ) then ! solid harmonics are 1.0
                        ! gradients of solid harmonics are 0.0
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           help_vec = coef * prim_grad_fac(1:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_sd_fac(1:vl,i_exp,i_ea)
                           help_vec3 = coef * prim_td_fac(1:vl,i_exp,i_ea)

                           D12 = 0
                           do D3=1,3
                           do D1=D3,3
                              do D2=D1,3
                                 D12 = D12 + 1
                                 sotd(:vl,D12,i_ea,i_orbs,i_pa) = sotd(:vl,D12,i_ea,i_orbs,i_pa) &
                                      & + help_vec3 * xyz(:,D1) * xyz(:,D2)*xyz(:,D3) &
                                      & + xyz(:,D1)*merge(help_vec2,zero,D2.eq.D3)    &
                                      & + xyz(:,D2)*merge(help_vec2,zero,D1.eq.D3)    &
                                      & + xyz(:,D3)*merge(help_vec2,zero,D1.eq.D2)
                              enddo
                           enddo
                           enddo
                           i_sa = i_sa + 1
                           i_orbs=i_orbs + 1
                        enddo! uncontracted exponents

                     else ! i_l .gt. 0  take solid harmonics into account
                        shm => sh%l(i_l)%m
                        sh => sha%sh(i_ea)
                        shg => sha%shg(i_ea)
                        shsd => sha%shs(i_ea)
                        shtd => sha%sht(i_ea)

                        shm_grad => shg%l(i_l)%m
                        shm_sd => shsd%l(i_l)%m
                        shm_td => shtd%l(i_l)%m
                        m = uas%m(i_fu)
                        coef = uas%c(i_fu)
                        do i_exp = 1, n_fcts !uab%N_uncontracted_fcts ! uncontracted exponents
                           sym_fac(:vl) = shm(:vl,m) * coef
                           help_vec =  coef * shm(:vl,m) * prim_grad_fac(:vl,i_exp,i_ea)
                           help_vec2 = coef * prim_orb_fac(:vl,i_exp,i_ea)
                           help_vec3 = coef * prim_grad_fac(:vl,i_exp,i_ea)
                           help_vec4 = coef * shm(:vl,m) * prim_sd_fac(:vl,i_exp,i_ea)
                           help_vec5 = coef * shm(:vl,m) * prim_td_fac(:vl,i_exp,i_ea)
                           help_vec6 = coef * prim_sd_fac(:vl,i_exp,i_ea)

                           DD(:vl,1,1)=shm_sd(:vl,1,m)
                           DD(:vl,1,2)=shm_sd(:vl,2,m)
                           DD(:vl,2,1)=shm_sd(:vl,2,m)
                           DD(:vl,1,3)=shm_sd(:vl,3,m)
                           DD(:vl,3,1)=shm_sd(:vl,3,m)
                           DD(:vl,2,2)=shm_sd(:vl,4,m)
                           DD(:vl,2,3)=shm_sd(:vl,5,m)
                           DD(:vl,3,2)=shm_sd(:vl,5,m)
                           DD(:vl,3,3)=shm_sd(:vl,6,m)

                           D12 = 0
                           do D3=1,3
                           do D1=D3,3
                              do D2=D1,3
                                 D12 = D12 + 1
!                                 sosd(:vl,D12,i_ea,i_orbs,i_pa) = sosd(:vl,D12,i_ea,i_orbs,i_pa) &
!                                      & +    xyz(:,D1) * xyz(:,D2)             * help_vec4 &    ! (1)
!                                      & + (  xyz(:,D1) * shm_grad(1:vl,D2,m) &                  ! (2)
!                                      &    + xyz(:,D2) * shm_grad(1:vl,D1,m) ) * help_vec3 &    ! (3)
!                                      & +    shm_sd(1:vl,D12,m)                * help_vec2 &    ! (4)
!                                      & + merge(help_vec,zero,D1.eq.D2)                         ! (5)
                                 sotd(:vl,D12,i_ea,i_orbs,i_pa) = sotd(:vl,D12,i_ea,i_orbs,i_pa) &
                                      & +    xyz(:,D1) * xyz(:,D2) * xyz(:,D3) * help_vec5 &              !(1.1)[d3 prim_sd_fac]
                                      & +    xyz(:,D1) * xyz(:,D2) * shm_grad(:vl,D3,m)* help_vec6 &      !(1.4)[d3 shm]
                                      & +    xyz(:,D2) * merge(help_vec4,zero,D1.eq.D3) &                 !(1.2)[d3 xyz(:,D1)]
                                      & +    xyz(:,D1) * merge(help_vec4,zero,D2.eq.D3) &                 !(1.3)[d3 xyz(:,D2)]
                                      & +    xyz(:,D1) * shm_grad(:vl,D2,m) * xyz(:,D3) * help_vec6 &     !(2.1)[d3 help_vec3]
                                      & +    xyz(:,D1) * DD(:vl,D2,D3) * help_vec3 &                      !(2.2)[d3 shm_grad(D2)]
                                      & +    merge(shm_grad(:vl,D2,m),zero,D1.eq.D3) * help_vec3 &        !(2.3)[d3 xyz(:,D1)]
                                      & +    xyz(:,D2) * shm_grad(1:vl,D1,m) * xyz(:,D3) * help_vec6 &    !(3.1)[d3 help_vec3]
                                      & +    xyz(:,D2) * DD(:vl,D1,D3) * help_vec3 &                      !(3.2)[d3 shm_grad(D1)]
                                      & +    merge(shm_grad(:vl,D1,m),zero,D2.eq.D3)*help_vec3 &          !(3.3)[d3 xyz(:,D2)]
                                      & +    DD(:vl,D1,D2)*xyz(:,D3) * help_vec3 &                        !(4.1)[d3 help_vec2]
                                      & +    shm_td(:vl,D12,m)*help_vec2 &                                !(4.2)[d3 shm_sd]
                                      & +    merge(xyz(:,D3)*help_vec4,zero,D1.eq.D2) &                   !(5.1)[d3 prim_grad_fac]
                                      & +    merge(shm_grad(:vl,D3,m)*help_vec3,zero,D1.eq.D2)            !(5.2)[d3 shm]
!!!                                       +xyz(:,D3) * prim_td_fac(:vl,i_exp,i_ea)*coef
                              enddo
                           enddo
                           enddo
                           i_sa = i_sa + 1
                           i_orbs=i_orbs+1
                        enddo!exponents
                     endif! i_l .eq. 0
                  enddo! contributing functions
                  ! copy symmetric secound derivatives
                  i_sa_intermediate = i_sa
                  i_orbs_intermediate = i_orbs
               enddo! independent fcts
            enddo! partners
         enddo! irreps
      enddo! l

    end subroutine calculate_nuc_3rd_der_v2

    subroutine calculate_fitfunctions( &
         orb, orbitalprojection, lmax, l_bas, r2_bas, N_gc, gc, &
         uncont_orb, i_sa_globcontr)
      !  Purpose: does calculation for a given unique atom
      !  for fitfunctions of charge or exchange
      !  performs l-specific contraction before symmetry adaption
      !  and if optional arguments are given global contraction in the end
    use symmetry_data_module ! description of irreps
      implicit none
      !------------ Declaration of formal parameters ---------------
      type(orbital_type), intent(inout) :: orb
      integer, pointer :: orbitalprojection(:,:) ! intent(in)
      !orbitalprojection(N_max_l,N_unique_atoms)
      integer, intent(in)  :: lmax
      type(unique_atom_basis_type)          :: l_bas(0:) ! l_bas(lmax) intent(in)
      type(unique_atom_basis_type)          :: r2_bas ! intent(in)
      target :: l_bas, r2_bas
      ! the following arguments must either all be omitted (=> no global contraction)
      ! or given (=> global contraction)
      integer, optional, intent(in)  :: N_gc
      type(unique_atom_glob_con_type), optional :: gc(:) ! intent(in)
      type(glob_con_intermediate_type), optional :: uncont_orb
      integer, optional, intent(in) :: i_sa_globcontr
      ! Index where to store global contractions
      !------------ Declaration of local variables -----------------
      integer             :: i_vec, i_ea, i_exp, i_gc, i_gcf, &
           i_pa, i_fu, i_c, i_if, i_l, i_bas, i_ir, &
           i_sa, i_sa_intermediate, i_sa_step
      integer             :: m
      real(kind=r8_kind)  :: alpha, exp_orb, norm
      ! intermediates
      real(kind=r8_kind)  :: coef
      ! contraction coefficient
      real(kind=r8_kind), pointer    :: r2(:) ! r2(vec_len)
      ! square of distance center to point calculated by
      ! solid_harmonics_module
      real(kind=r8_kind), pointer    :: shm(:,:)
      real(kind=r8_kind), pointer    :: so(:,:,:)
      real(kind=r8_kind), pointer    :: uncont_orb_l(:,:,:)
      logical                        :: is_r2function
      logical                        :: do_glob_con
      type(unique_atom_basis_type),         pointer :: uab
      type(unique_atom_partner_type),       pointer :: uap
      type(unique_atom_symadapt_type),      pointer :: uas
      type(solid_harmonics_type),           pointer :: sh
      !------------ Executable code --------------------------------

      ! interpret optional arguments
      if ( present(N_gc) ) then
         if ( .not. present(gc) .or. &
              .not. present(uncont_orb) .or. &
              .not. present(i_sa_globcontr) ) then
            call error_handler( &
                 "orbital_calculate: contract_symadapt: invalid arguments")
         else
            do_glob_con = .true.
         endif
      else
         if ( present(gc) .or. present(i_sa_globcontr) .or. present(uncont_orb) ) then
            call error_handler( &
                 "orbital_calculate: contract_symadapt: invalid arguments")
         else
            do_glob_con = .false.
         endif
      endif
      i_ir = get_totalsymmetric_irrep()
      i_pa = 1
      so => orb%o


      ! loop over l
      do i_bas = -1, lmax
         ! interpret i_bas to set i_l and is_r2function
         if ( i_bas .eq. -1 ) then
            i_l = 0
            is_r2function = .false.
            uab => l_bas(0)
         elseif ( i_bas .eq. 0 ) then
            i_l = 0
            is_r2function = .true.
            uab => r2_bas
         else
            i_l = i_bas
            is_r2function = .false.
            uab => l_bas(i_l)
         endif

         if (do_glob_con ) uncont_orb_l => uncont_orb%l(i_bas)%oo

         ! loop over equal atoms
         do i_ea = 1, ua%N_equal_atoms
            sh => sha%sh(i_ea)
            r2 => sh%r2

            ! calculate primitive orbital factors for all exponents and gridpoints
            if ( is_r2function ) then
               do i_exp = 1, uab%N_exponents
                  alpha = uab%exponents(i_exp)
                  do i_vec = 1, gridlength
                     exp_orb = alpha * r2(i_vec)
                     prim_orb_fac(i_vec,i_exp,i_ea) = r2(i_vec) * exp( - exp_orb )
                  enddo
               enddo
            else
               do i_exp = 1, uab%N_exponents
                  alpha = uab%exponents(i_exp)
                  norm = uab%norms(i_exp)
                  do i_vec = 1, gridlength
                     exp_orb = alpha * r2(i_vec)
                     prim_orb_fac(i_vec,i_exp,i_ea) = norm * exp( - exp_orb )
                  enddo
               enddo
            endif


            ! calculate contracted orbitals factors for all contractions and gridpoints
            do i_c = 1, uab%N_contracted_fcts ! contractions
               do i_vec = 1, gridlength
                  cont_orb_fac(i_vec,i_c,i_ea) = 0.0_r8_kind
               enddo
               do i_exp = 1, uab%N_exponents ! contributing exponents
                  coef = uab%contractions(i_exp,i_c)
                  do i_vec = 1, gridlength
                     cont_orb_fac(i_vec,i_c,i_ea) = cont_orb_fac(i_vec,i_c,i_ea) + &
                          coef * prim_orb_fac(i_vec,i_exp,i_ea)
                  enddo
               enddo! contributing exponents
            enddo! contractions


         enddo! equal atoms



         ! symmetry adaption:
         ! for all independent fcts,
         !         uncontracted exponents and contractions,
         !         gridpoints:
         !   sum over equal atoms and m
         !
         ! see haeder of obitalprojection_module to understand order of
         ! results in so calculated.

         uap => ua%symadapt_partner(i_ir,i_l)
         i_sa =orbitalprojection(i_bas,i_ua)
         i_sa_intermediate = i_sa
         i_sa_step = uap%N_independent_fcts

         if ( i_l .eq. 0 ) then ! solid harmonics are 1.0 and symmetry adaption is
            ! performed as sum over all equal atoms with wight 1.0.
            ! there is only one independent function

            i_if = 1
            if ( do_glob_con ) then
               do i_exp = 1, uab%N_exponents ! all exponents
                  do i_vec = 1, gridlength
                     uncont_orb_l(i_vec,i_exp,i_if) =  0.0_r8_kind
                  enddo
               enddo! all exponents
            endif
            do i_ea = 1,ua%N_equal_atoms
               i_sa = i_sa_intermediate
               if ( do_glob_con ) then
                  do i_exp = 1, uab%N_exponents ! all exponents
                     do i_vec = 1, gridlength
                        uncont_orb_l(i_vec,i_exp,i_if) =  &
                             uncont_orb_l(i_vec,i_exp,i_if) + &
                             prim_orb_fac(i_vec,i_exp,i_ea)
                     enddo
                  enddo! all exponents
                  i_sa = i_sa + uab%N_uncontracted_fcts
               else ! .not. do_glob_con
                  do i_exp = 1, uab%N_uncontracted_fcts ! uncontracted exponents
                     do i_vec = 1, gridlength
                        so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                             prim_orb_fac(i_vec,i_exp,i_ea)
                     enddo
                     i_sa = i_sa + 1
                  enddo! uncontracted exponents
               endif! do_glob_con
               do i_c = 1, uab%N_contracted_fcts  ! contractions
                  do i_vec = 1, gridlength
                     so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                          cont_orb_fac(i_vec,i_c,i_ea)
                  enddo
                  i_sa = i_sa + 1
               enddo! contractions
            enddo! equal atoms
            i_sa = i_sa_intermediate
            if ( do_glob_con ) then
               i_sa = i_sa + uab%N_uncontracted_fcts
            else ! .not. do_glob_con
               do i_exp = 1, uab%N_uncontracted_fcts ! uncontracted exponents
                  i_sa = i_sa + 1
               enddo! uncontracted exponents
            endif! do_glob_con
            do i_c = 1, uab%N_contracted_fcts  ! contractions
               i_sa = i_sa + 1
            enddo! contractions
            ! in case of global contraction copy to result
            i_sa = i_sa_intermediate
            if ( do_glob_con ) then
               do i_exp = 1, uab%N_uncontracted_fcts ! uncontacted exponents
                  do i_vec = 1, gridlength
                     so(i_vec,i_sa,i_pa) = uncont_orb_l(i_vec,i_exp,i_if)
                  enddo
                  i_sa = i_sa + 1
               enddo! uncontacted exponents
            endif

         else ! i_l .gt. 0  take solid harmonics into account and allow several
            ! independent functions. Use description of symmetry adaption
            ! given in uas

            do i_if = 1, uap%N_independent_fcts ! independent fcts
               uas => uap%symadapt(i_if,i_pa)
               if ( do_glob_con ) then
                  do i_exp = 1, uab%N_exponents ! all exponents
                     do i_vec = 1, gridlength
                        uncont_orb_l(i_vec,i_exp,i_if) =  0.0_r8_kind
                     enddo
                  enddo! all exponents
               endif
               do i_fu = 1, uas%N_fcts ! contributing functions
                  i_sa = i_sa_intermediate
                  i_ea = uas%I_equal_atom(i_fu)
                  sh => sha%sh(i_ea)
                  shm => sh%l(i_l)%m
                  m = uas%m(i_fu)
                  coef = uas%c(i_fu)
                  if ( do_glob_con ) then
                     do i_exp = 1, uab%N_exponents ! all exponents
                        do i_vec = 1, gridlength
                           uncont_orb_l(i_vec,i_exp,i_if) =  &
                                uncont_orb_l(i_vec,i_exp,i_if) + &
                                shm(i_vec,m) * coef * &
                                prim_orb_fac(i_vec,i_exp,i_ea)
                        enddo
                     enddo! all exponents
                     i_sa = i_sa + uab%N_uncontracted_fcts * i_sa_step
                  else ! .not. do_glob_con
                     do i_exp = 1, uab%N_uncontracted_fcts ! uncontracted exponents
                        do i_vec = 1, gridlength
                           so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                                shm(i_vec,m) * coef * &
                                prim_orb_fac(i_vec,i_exp,i_ea)
                        enddo
                        i_sa = i_sa + i_sa_step
                     enddo! uncontracted exponents
                  endif! do_glob_con
                  do i_c = 1, uab%N_contracted_fcts  ! contractions
                     do i_vec = 1, gridlength
                        so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                             shm(i_vec,m) * coef * &
                             cont_orb_fac(i_vec,i_c,i_ea)
                     enddo
                     i_sa = i_sa + i_sa_step
                  enddo! contractions
               enddo! contributing functions
               i_sa_intermediate = i_sa_intermediate + 1
            enddo! independent fcts

            if ( do_glob_con ) then
               ! copy uncontracted orbitals to result
               uap => ua%symadapt_partner(i_ir,i_l)
               i_pa = 1 ! glob contractions only for fitfcts
               i_sa =orbitalprojection(i_bas,i_ua)
               do i_exp = 1, uab%N_uncontracted_fcts ! uncontracted exponents
                  do i_if = 1, uap%N_independent_fcts ! independent fcts
                     do i_vec = 1, gridlength
                        so(i_vec,i_sa,i_pa) = uncont_orb_l(i_vec,i_exp,i_if)
                     enddo
                     i_sa = i_sa + 1
                  enddo! independent fcts
               enddo! uncontracted exponents
            endif

         endif! i_l .gt. 0

      enddo! l


      if ( do_glob_con ) then
         ! calculate globally contracted functions
         ! and store them in the end of so
         i_sa = i_sa_globcontr
         i_pa = 1
         do i_gc = 1, N_gc ! global contractions
            do i_gcf = 1, gc(i_gc)%N_contributing_fcts ! contributing_fcts
               i_bas = gc(i_gc)%l(i_gcf)
               i_exp  = gc(i_gc)%index_exp(i_gcf)
               i_if = gc(i_gc)%index_ind_fct(i_gcf)
               coef = gc(i_gc)%coefs(i_gcf)
               uncont_orb_l => uncont_orb%l(i_bas)%oo
               do i_vec = 1, gridlength
                  so(i_vec,i_sa,i_pa) = so(i_vec,i_sa,i_pa) + &
                       coef * uncont_orb_l(i_vec,i_exp,i_if)
               enddo
            enddo! contributing_fcts
            i_sa = i_sa + 1
         enddo! global contractions
      endif! calculate globally contracted functions


    end subroutine calculate_fitfunctions


  end subroutine orbital_calculate

    subroutine do_angular(vl,Y,cs,ms,sh,shg)
      implicit none
      integer(i4_kind), intent(in)  :: vl
      real(r8_kind)   , intent(out) :: Y(:,0:)    ! (vl,0:3)
      real(r8_kind)   , intent(in)  :: cs(:)
      integer(i4_kind), intent(in)  :: ms(:)
      real(r8_kind)   , intent(in)  :: sh(:,:)    ! (vl  ,2L+1)
      real(r8_kind)   , intent(in)  :: shg(:,:,:) ! (vl,3,2L+1)
      optional :: shg
      ! *** end of interface ***

      integer(i4_kind) :: i_fu,m
      real(r8_kind)    :: c

      ASSERT(size(Y,1)>=vl)
      ASSERT(size(Y,2)>=1)

      Y(:vl,VALUE) = zero

      do i_fu = 1, size(cs)
         c = cs(i_fu)
         m = ms(i_fu)
         Y(:vl,VALUE) = Y(:vl,VALUE) + c * sh(:vl,m)
      enddo

      if(.not.present(shg)) return

      ASSERT(size(Y,2)==4)

      Y(:vl,DX:DZ) = zero

      do i_fu = 1, size(cs)
         c = cs(i_fu)
         m = ms(i_fu)
         Y(:vl,DX:DZ) = Y(:vl,DX:DZ) + c * shg(:vl,1:3,m)
      enddo
    end subroutine do_angular

    subroutine do_spinor(isa,NEX,rad,LL,sh,vl,o,REIM,uas,beta)
      use unique_atom_module ! desription of unique atoms, must be calculated before
      implicit none
      integer(i4_kind)             , intent(in)  :: isa
      integer(i4_kind)             , intent(in)  :: NEX
      real(r8_kind)                , intent(in)  :: rad(:,:) ! (vl,>=NEX)
      integer(i4_kind)             , intent(in)  :: LL
      real(r8_kind)                , intent(in)  :: sh(:,:)  ! (vl,m)
      integer(i4_kind)             , intent(in)  :: vl
      real(r8_kind)                , intent(out) :: o(:,:,:)
      ! o(vl,NBAS,NPA)
      integer(i4_kind)             , intent(in)  :: REIM
      type(unique_atom_sa_int_type), intent(in)  :: uas(:,:)
      real(r8_kind)                , intent(in)  :: beta
      ! uas(NIF,NPA)
      ! *** end of interface ***

      integer(i4_kind) :: NIF,NPA
      integer(i4_kind) :: i_sa_cur,i_sa,i_if,i_pa,i_exp
      integer(i4_kind) :: vll
      real(r8_kind)    :: Y(vl,0:0)


      NIF = size(uas,1)
      NPA = size(uas,2)
      ASSERT(NPA==size(o,3))

      if(LL==0)then
         ! for L=0 spherical harmonic is constant
         vll = 1
         ASSERT(vl>0)
      else
         vll = vl
      endif

      i_sa_cur = isa

      do i_if = 1, NIF
         do i_pa = 1, NPA

            ! nullify:
            if(beta==zero)then
               i_sa = i_sa_cur
               do i_exp = 1, NEX
                  o(1:vl,i_sa,i_pa) = zero
                  i_sa = i_sa + 1
               enddo
            endif

            select case(REIM)
            case (RE)
               call do_angular(vll, Y, uas(i_if,i_pa)%re, uas(i_if,i_pa)%m, sh)
            case (IM)
               call do_angular(vll, Y, uas(i_if,i_pa)%im, uas(i_if,i_pa)%m, sh)
            case default
               ABORT('no such REIM')
            end select

            i_sa = i_sa_cur
            if ( LL .eq. 0 ) then ! solid harmonics are 1.0
               ! gradients of solid harmonics are 0.0
               do i_exp = 1, NEX
                  o(1:vl,i_sa,i_pa) = o(1:vl,i_sa,i_pa) + &
                       Y(1,VALUE) * rad(1:vl,i_exp)
                  i_sa = i_sa + 1
               enddo
            else ! LL .gt. 0,  take solid harmonics into account
               do i_exp = 1, NEX
                  o(1:vl,i_sa,i_pa) = o(1:vl,i_sa,i_pa) &
                       + Y(:vl,VALUE) * rad(1:vl,i_exp)
                  i_sa = i_sa + 1
               enddo
            endif! i_l .eq. 0
         enddo ! i_pa
         i_sa_cur = i_sa_cur + NEX
      enddo
    end subroutine do_spinor

    subroutine do_spinor_grad(isa,NEX,rad0,rad1,LL,sh,sh1,shg,vl,o,REIM,uas,beta)
      use unique_atom_module ! desription of unique atoms, must be calculated before
      implicit none
      integer(i4_kind)              , intent(in)  :: isa
      integer(i4_kind)              , intent(in)  :: NEX
      real(r8_kind)                 , intent(in)  :: rad0(:,:) ! (vl,>=NEX)
      real(r8_kind)                 , intent(in)  :: rad1(:,:) ! (vl,>=NEX)
      integer(i4_kind)              , intent(in)  :: LL
      real(r8_kind)                 , intent(in)  :: sh(:,:)   ! (vl  ,2*L+1)
      real(r8_kind)                 , intent(in)  :: sh1(:,:)  ! (vl  ,2*1+1)
      real(r8_kind)                 , intent(in)  :: shg(:,:,:)! (vl,3,2*L+1)
      integer(i4_kind)              , intent(in)  :: vl
      real(r8_kind)                 , intent(out) :: o(:,:,:,:)
      ! o(vl,3,NBAS,NPA)
      integer(i4_kind)              , intent(in)  :: REIM
      type(unique_atom_sa_int_type) , intent(in)  :: uas(:,:)
      ! uas(NIF,NPA)
      real(r8_kind)                 , intent(in)  :: beta
      ! *** end of interface ***

      integer(i4_kind) :: NIF,NPA
      integer(i4_kind) :: i_sa_cur,i_sa,i_if,i_pa,i_exp
      integer(i4_kind) :: vll
      real(r8_kind)    :: Y(vl,0:3)
      real(r8_kind)    :: YDR(vl), RDY(vl)
      real(r8_kind)    :: rv(vl,3)

      NIF = size(uas,1)
      NPA = size(uas,2)
      ASSERT(NPA==size(o,4))

      if(LL==0)then
         ! for L=0 spherical harmonic is constant
         vll = 1
         ASSERT(vl>0)
      else
         vll = vl
      endif

      rv(:,DX) = sh1(:vl,2)
      rv(:,DY) = sh1(:vl,3)
      rv(:,DZ) = sh1(:vl,1)

      i_sa_cur = isa

      do i_if = 1, NIF
         do i_pa = 1, NPA

            ! nullify:
            if(beta==zero)then
               i_sa = i_sa_cur
               do i_exp = 1, NEX
                  o(1:vl,1:3,i_sa,i_pa) = zero
                  i_sa = i_sa + 1
               enddo
            endif

            select case(REIM)
            case (RE)
               call do_angular(vll, Y, uas(i_if,i_pa)%re, uas(i_if,i_pa)%m, sh, shg)
            case (IM)
               call do_angular(vll, Y, uas(i_if,i_pa)%im, uas(i_if,i_pa)%m, sh, shg)
            case default
               ABORT('no such REIM')
            end select

            i_sa = i_sa_cur
            if ( LL .eq. 0 ) then ! solid harmonics are 1.0
               ! gradients of solid harmonics are 0.0
               do i_exp = 1, NEX

                  YDR = Y(1,VALUE) * rad1(1:vl,i_exp) * rv(:,DX)
                  o(1:vl,DX,i_sa,i_pa) = o(1:vl,DX,i_sa,i_pa) + YDR

                  YDR = Y(1,VALUE) * rad1(1:vl,i_exp) * rv(:,DY)
                  o(1:vl,DY,i_sa,i_pa) = o(1:vl,DY,i_sa,i_pa) + YDR

                  YDR = Y(1,VALUE) * rad1(1:vl,i_exp) * rv(:,DZ)
                  o(1:vl,DZ,i_sa,i_pa) = o(1:vl,DZ,i_sa,i_pa) + YDR

                  i_sa = i_sa + 1
               enddo! uncontracted exponents
            else ! LL .gt. 0,  take solid harmonics into account
               do i_exp = 1, NEX

                  YDR = Y(:vl,VALUE) * rad1(1:vl,i_exp) * rv(:,DX)
                  RDY = Y(:vl,DX   ) * rad0(1:vl,i_exp)

                  o(1:vl,DX,i_sa,i_pa) = o(1:vl,DX,i_sa,i_pa) + YDR + RDY

                  YDR = Y(:vl,VALUE) * rad1(1:vl,i_exp) * rv(:,DY)
                  RDY = Y(:vl,DY   ) * rad0(1:vl,i_exp)

                  o(1:vl,DY,i_sa,i_pa) = o(1:vl,DY,i_sa,i_pa) + YDR + RDY

                  YDR = Y(:vl,VALUE) * rad1(1:vl,i_exp) * rv(:,DZ)
                  RDY = Y(:vl,DZ   ) * rad0(1:vl,i_exp)

                  o(1:vl,DZ,i_sa,i_pa) = o(1:vl,DZ,i_sa,i_pa) + YDR + RDY

                  i_sa = i_sa + 1
               enddo
            endif! i_l .eq. 0
         enddo ! i_pa
         i_sa_cur = i_sa_cur + NEX
      enddo
    end subroutine do_spinor_grad

  !*************************************************************
  subroutine orbital_write(outunit,orbs_ob,orb_ch,orb_xc,grads, &
       sec_der,complete)
    !  Purpose: writes debug output
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc, fit_coeff_n_cd
    use symmetry_data_module ! description of irreps
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer, intent(in) :: outunit
    type(orbital_type), pointer, optional :: orbs_ob(:) ! orbs_ob(N_irreps) intent(in)
    type(orbital_gradient_type), intent(in), optional :: grads(:) ! grads(N_irreps) intent(in)
    type(orbital_sec_der_type), intent(in), optional :: sec_der(:) ! sec_der(N_irreps) intent(in)
    type(orbital_type), intent(in), optional :: orb_ch
    type(orbital_type), intent(in), optional :: orb_xc
    logical, intent(in), optional :: complete ! write complete list for all
    ! gridpoints, otherwise only for point 10. default: .false.
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer  :: i_vec, i_sa, i_pa, i_ir, i, i_x
    ! loop indices
    real(kind=r8_kind), pointer :: so(:,:,:), sog(:,:,:,:)
    real(kind=r8_kind), pointer :: sosd(:,:,:,:)
    logical :: complete_local
    !------------ Executable code ------------------------------------

    if ( present( complete ) ) then
       complete_local = complete
    else
       complete_local = .false.
    endif

    if ( present(orbs_ob) ) then
       write (outunit,*)
       write (outunit,*) "Orbitals atomic basis"
       write (outunit,*)
       i = 0
       do i_ir = 1, symmetry_data_n_irreps()
          so => orbs_ob(i_ir)%o
          !write(outunit,*) "IRREP: ", i_ir
          do i_sa = 1, symmetry_data_dimension(i_ir)
             do i_pa = 1, symmetry_data_n_partners(i_ir)
                i = i + 1
                !write(outunit,'("irrep ",i2,"  partner ",i2," : ",2i4)')  i_ir, i_pa, i_sa, i
                if ( complete_local) then
                   do i_vec=1,vec_len
                      write(outunit,'(2I4,E25.14E3)')  i, i_vec, so(i_vec,i_sa,i_pa)
                   enddo
                else
                   i_vec = 10
                   write(outunit,'(I4,E25.14E3)')  i, so(i_vec,i_sa,i_pa)
                endif
             enddo
          enddo
       enddo
    endif

    if ( present(grads) ) then
       write (outunit,*)
       write (outunit,*) "Gradients atomic basis"
       write (outunit,*)
       i = 0
       do i_ir = 1, symmetry_data_n_irreps()
          sog => grads(i_ir)%o
          !write(outunit,*) "IRREP: ", i_ir
          do i_sa = 1, symmetry_data_dimension(i_ir)
             do i_pa = 1, symmetry_data_n_partners(i_ir)
                i = i + 1
                !write(outunit,'("irrep ",i2,"  partner ",i2," : ",2i4)')  i_ir, i_pa, i_sa, i
                if ( complete_local) then
                   do i_vec=1,vec_len
                      write(outunit,'(2I4,3E25.14E3)')  i, i_vec, (sog(i_vec,i_x,i_sa,i_pa), i_x=1,3)
                   enddo
                else
                   i_vec = 10
                   write(outunit,'(I4,3E25.14E3)')  i, (sog(i_vec,i_x,i_sa,i_pa), i_x=1,3)
                endif
             enddo
          enddo
       enddo
    endif

    if ( present(sec_der) ) then
       write (outunit,*)
       write (outunit,*) "secound derivatives of atomic basis"
       write (outunit,*)
       i = 0
       do i_ir = 1, symmetry_data_n_irreps()
          sosd => sec_der(i_ir)%o
          !write(outunit,*) "IRREP: ", i_ir
          do i_sa = 1, symmetry_data_dimension(i_ir)
             do i_pa = 1, symmetry_data_n_partners(i_ir)
                i = i + 1
                !write(outunit,'("irrep ",i2,"  partner ",i2," : ",2i4)')  i_ir, i_pa, i_sa, i
                if ( complete_local) then
                   do i_vec=1,vec_len
                      write(outunit,'(2I4,6E15.4E3)')  i, i_vec, &
                           sosd(i_vec,XX,i_sa,i_pa), sosd(i_vec,XY,i_sa,i_pa), &
                           sosd(i_vec,XZ,i_sa,i_pa), sosd(i_vec,YY,i_sa,i_pa), &
                           sosd(i_vec,YZ,i_sa,i_pa), sosd(i_vec,ZZ,i_sa,i_pa)
                   enddo
                else
                   i_vec=10
                   write(outunit,'(I4,6E15.4E3)')  i, &
                        sosd(i_vec,XX,i_sa,i_pa), sosd(i_vec,XY,i_sa,i_pa), &
                        sosd(i_vec,XZ,i_sa,i_pa), sosd(i_vec,YY,i_sa,i_pa), &
                        sosd(i_vec,YZ,i_sa,i_pa), sosd(i_vec,ZZ,i_sa,i_pa)
                endif
             enddo
          enddo
       enddo
    endif

    if ( present(orb_ch) ) then
       write (outunit,*)
       write (outunit,*) "Orbitals charge fitfunctions"
       write (outunit,*)
       i = 0
       i_ir = get_totalsymmetric_irrep()
       i_pa = 1
       so => orb_ch%o
       do i_sa = 1, fit_coeff_n_ch()
          i = i + 1
          if ( complete_local) then
             write(outunit,'("irrep ",i2,"  partner ",i2," : ",2i4)')  i_ir, i_pa, i_sa, i
             write(outunit,'(3E25.14E3)')  (so(i_vec,i_sa,i_pa), i_vec=1,vec_len)
          else
             write(outunit,'(3E25.14E3,4I4)')  (so(i_vec,i_sa,i_pa), i_vec=10,12), i_ir, i_pa, i_sa, i
          endif
       enddo
    endif

    if ( present(orb_xc) ) then
       write (outunit,*)
       write (outunit,*) "Orbitals exchange fitfunctions"
       write (outunit,*)
       i = 0
       i_ir = get_totalsymmetric_irrep()
       i_pa = 1
       so => orb_xc%o
       do i_sa = 1, fit_coeff_n_xc()
          i = i + 1
          if ( complete_local) then
             write(outunit,'("irrep ",i2,"  partner ",i2," : ",2i4)')  i_ir, i_pa, i_sa, i
             write(outunit,'(3E25.14E3)')  (so(i_vec,i_sa,i_pa), i_vec=1,vec_len)
          else
             write(outunit,'(3E25.14E3,4I4)')  (so(i_vec,i_sa,i_pa), i_vec=10,12), i_ir, i_pa, i_sa, i
          endif
       enddo
    endif

    write (outunit,*)
    write (outunit,*)

  end subroutine orbital_write
  !*************************************************************


  !*************************************************************
  subroutine fit_fct_allocate(type, fcts, grads, sec_ders, nuc_grads, &
       sec_nuc_ders)
    !  Purpose: allocates the arrays to hold fit functions
    !  and their various derivates on a set of grid points
    use options_module, only : options_orbitals_in_memory
    use fit_coeff_module, only: fit_coeff_n_ch, fit_coeff_n_xc, fit_coeff_n_cd
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalstore_module, only: orbitalstore_initialized, orbitalstore_allocate
    use symmetry_data_module ! description of irreps
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=2), intent(in)                           :: type
    type(orbital_type)                          , optional :: fcts
    type(orbital_gradient_type)                 , optional :: grads
    type(orbital_sec_der_type)                  , optional :: sec_ders
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:)
    type(orbital_nuclear_sec_der_type) , pointer, optional :: sec_nuc_ders(:)
    !** End of interface *****************************************
!AG #ifdef FPP_NOMDA
!AG     call error_handler("om/fit_fct_allocate: recompile w/o FPP_NOMDA")
!AG #else
    !------------ Declaration of local variables ---------------------
    integer :: n_dim, i_ir,  i_ua, i_ma, n_atoms, n_orbs
    logical :: is_charge_fit
    !------------ Declaration of local handles -------------------
    type(unique_atom_type), pointer :: ua
    !------------ Executable code ------------------------------------
    select case (type)
    case ("ch")
       is_charge_fit = .true.
       n_dim = fit_coeff_n_ch()
    case ("xc")
       is_charge_fit = .false.
       n_dim = fit_coeff_n_xc()
    case default
       call error_handler("fit_fct_allocate: unknown type of fit functions")
    end select
    i_ir = get_totalsymmetric_irrep()

    if (present(fcts)) then
      if(.NOT.orbitalstore_initialized .or. .NOT.options_orbitals_in_memory() ) then
       allocate( fcts%o(vec_len,n_dim,1), stat=status(33))
       if (status(33) /= 0) call error_handler( &
            "fit_fct_allocate: allocate of fcts%o failed")
      end if
        if(orbitalstore_initialized) call orbitalstore_allocate(fcts=fcts,fcts_dim=n_dim)
    endif
!   print*, ' alloc alloc',shape(fcts%o),vec_len,n_dim,type,orbitalstore_initialized,present(fcts)
    if (present(grads)) then
      if(.NOT.orbitalstore_initialized .or. .NOT.options_orbitals_in_memory() ) then
       allocate( grads%o(vec_len,3,n_dim,1), stat=status(34))
       if (status(34) /= 0) call error_handler( &
            "fit_fct_allocate: allocate of grads%o failed")
      endif
     if(orbitalstore_initialized) call orbitalstore_allocate(fcts_grads=grads,fcts_dim=n_dim)
    endif
    if (present(sec_ders)) then
       allocate( sec_ders%o(vec_len,6,n_dim,1), stat=status(35))
       if (status(35) /= 0) call error_handler( &
            "fit_fct_allocate: allocate of sec_ders%o failed")
    endif
    if (present(nuc_grads)) then
       allocate( nuc_grads(N_moving_unique_atoms), stat=status(36))
       if (status(36) /= 0) call error_handler( &
            "fit_fct_allocate: allocate of nuc_grads failed")
    endif
    if (present(sec_nuc_ders)) then
       allocate( sec_nuc_ders(N_moving_unique_atoms), stat=status(37))
       if (status(37) /= 0) call error_handler( &
            "fit_fct_allocate: allocate of sec_nuc_ders failed")
    endif

    if (present(nuc_grads) .or. present(sec_nuc_ders)) then
       do i_ma=1,N_moving_unique_atoms
          i_ua = moving_unique_atom_index(i_ma)
          ua => unique_atoms(i_ua)
          ! find dimensions
          n_atoms = ua%N_equal_atoms
          if (is_charge_fit) then
             n_orbs = number_of_fit_fcts(ua%lmax_ch,ua%r2_ch,ua%l_ch)
             n_orbs = n_orbs + ua%N_glob_cons_ch
          else
             n_orbs = number_of_fit_fcts(ua%lmax_xc,ua%r2_xc,ua%l_xc)
             n_orbs = n_orbs + ua%N_glob_cons_xc
          endif
          ! allocate components
          if (present(nuc_grads)) then
             allocate( nuc_grads(i_ma)%o(vec_len,3,n_atoms,n_orbs,1), &
                  stat=status(38) )
             if (status(38) /= 0) call error_handler( &
                  "fit_fct_allocate: allocate of nuc_grads(i_ua)%o failed")
          endif
          if (present(sec_nuc_ders)) then
             allocate( sec_nuc_ders(i_ma)%o(vec_len,6,n_atoms,n_orbs,1), &
                  stat=status(39) )
             if (status(39) /= 0) call error_handler( &
                  "fit_fct_allocate: allocate of sec_nuc_ders(i_ua)%o failed")
          endif
       enddo
    endif

  contains

    function number_of_fit_fcts(lmax,r2_bas,l_bas)
       !  Purpose: returns the number of fit_functions contained in the
       !  atomic fit function basis set passed via lmax, r2_bas, and l_bas
       implicit none
       !------------ Declaration of formal parameters ---------------
       integer(kind=i4_kind)       , intent(in) :: lmax
       type(unique_atom_basis_type)             :: r2_bas   ! intent(in)
       type(unique_atom_basis_type)             :: l_bas(0:) ! intent(in)
       target :: r2_bas, l_bas
       integer(kind=i4_kind)                    :: number_of_fit_fcts
       !** End of interface *****************************************
       !------------ Declaration of local variables -----------------
       integer(kind=i4_kind)                   :: i_bas
       type(unique_atom_basis_type)  , pointer :: uab
       type(unique_atom_partner_type), pointer :: uap
       !------------ Executable code --------------------------------
       number_of_fit_fcts = 0
       do i_bas=-1,lmax
         select case (i_bas)
         case (-1) ! s-type
           uab => l_bas(0)
           uap => ua%symadapt_partner(i_ir,0)
         case (0) ! r^2-type
           uab => r2_bas
           uap => ua%symadapt_partner(i_ir,0)
         case default ! l>0-type
           uab => l_bas(i_bas)
           uap => ua%symadapt_partner(i_ir,i_bas)
         end select
         number_of_fit_fcts = number_of_fit_fcts + uap%N_independent_fcts &
              * ( uab%N_uncontracted_fcts + uab%N_contracted_fcts )
       end do
    end function number_of_fit_fcts
!AG #endif
  end subroutine fit_fct_allocate
  !*************************************************************


  !*************************************************************
  subroutine fit_fct_free(type, fcts, grads, sec_ders, nuc_grads, sec_nuc_ders)
    !  Purpose: deallocates the arrays ment to hold fit functions
    !  and their various derivates on a set of grid points
    use options_module, only : options_orbitals_in_memory, &
                               options_orbitals_on_file
    use unique_atom_module ! desription of unique atoms, must be calculated before
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=2), intent(in)                           :: type
    type(orbital_type)                          , optional :: fcts
    type(orbital_gradient_type)                 , optional :: grads
    type(orbital_sec_der_type)                  , optional :: sec_ders
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:)
    type(orbital_nuclear_sec_der_type) , pointer, optional :: sec_nuc_ders(:)
    !** End of interface *****************************************
!AG #ifdef FPP_NOMDA
!AG     call error_handler("om/fit_fct_free: recompile w/o FPP_NOMDA")
!AG #else
    !------------ Declaration of local variables ---------------------
    integer :: i_ma
    logical :: is_charge_fit
    !------------ Executable code ------------------------------------
    select case (type)
    case ("ch")
       is_charge_fit = .true.
    case ("xc")
       is_charge_fit = .false.
    case default
       call error_handler("fit_fct_free: unknown type of fit functions")
    end select

    if (present(fcts)) then
       if(.NOT. options_orbitals_in_memory() ) then
       deallocate( fcts%o, stat=status(33))
       ASSERT(status(33).eq.0)
       status(33)=1
       end if
    endif
    if (present(grads)) then
       if(.NOT. options_orbitals_in_memory() ) then
       deallocate( grads%o, stat=status(34))
       ASSERT(status(34).eq.0)
       status(34)=1
       endif
    endif
    if (present(sec_ders)) then
       deallocate( sec_ders%o, stat=status(35))
       if (status(35) /= 0) call error_handler( &
            "fit_fct_free: deallocate of sec_ders%o failed")
       status(35)=1
    endif
    if (present(nuc_grads)) then
       do i_ma=1,N_moving_unique_atoms
          deallocate( nuc_grads(i_ma)%o, stat=status(38) )
          if (status(38) /= 0) call error_handler( &
               "fit_fct_free: deallocate of nuc_grads(i_ua)%o failed")
       status(38)=1
       enddo
       deallocate( nuc_grads, stat=status(36))
       if (status(36) /= 0) call error_handler( &
            "fit_fct_free: deallocate of nuc_grads failed")
       status(36)=1
    endif
    if (present(sec_nuc_ders)) then
       do i_ma=1,N_moving_unique_atoms
          MEMLOG(-size(sec_nuc_ders(i_ma)%o))
          deallocate( sec_nuc_ders(i_ma)%o, stat=status(39) )
          if (status(39) /= 0) call error_handler( &
               "fit_fct_free: deallocate of sec_nuc_ders(i_ua)%o failed")
          status(39)=1
       enddo
       deallocate( sec_nuc_ders, stat=status(37))
       if (status(37) /= 0) call error_handler( &
            "fit_fct_free: deallocate of sec_nuc_ders failed")
       status(37)=1
    endif
!AG#endif
  end subroutine fit_fct_free


  subroutine fit_fct_calculate(grid_points,N_points,type, &
       fcts, grads, sec_ders, nuc_grads, sec_nuc_ders)
    !  Purpose: calculates the values of the fit functions
    !  and their various derivates on a set of grid points
    use options_module, only : options_orbitals_in_memory, &
                               options_orbitals_on_file
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalprojection_module ! projection of indices to metaindex used in scf
    use orbitalstore_module, only: orbitalstore_initialized, stored_fcts, &
         stored_fcts_grads, orbitalstore_rw
    use symmetry_data_module ! description of irreps
    use solid_harmonics_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in), optional :: grid_points(:,:)
    ! give this argument only at first call for given grid_points.
    ! Solid harmonics are evaluated when this argument is present and
    ! remain stored
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    character(len=2), intent(in)                          :: type
    character(len=2)                                      :: typef
    ! either "ch" for charge fit functions or "xc" for exchange fit functions
    type(orbital_type)                 , target, optional :: fcts
    type(orbital_gradient_type)        , target, optional :: grads
    type(orbital_sec_der_type)         , target, optional :: sec_ders
    type(orbital_nuclear_gradient_type), target, optional :: nuc_grads(:)
    type(orbital_nuclear_sec_der_type) , target, optional :: sec_nuc_ders(:)
    ! the various types of fit function arrays to load (all intent(out))
    !** End of interface *****************************************
!AG #ifdef FPP_NOMDA
!AG     call error_handler("om/fit_fct_calculate: recompile w/o FPP_NOMDA")
!AG #else
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind), parameter :: zero = 0.0_r8_kind
    integer(kind=i4_kind) :: vl, i_ma, i_ua, i_ea, i_ir, &
                             i_bas, i_l, i_exp, i_c, i_if, i_fu, m, &
                             i_gc, i_gcf, i_step, i_sa, i_sa_of_i_if, &
                             i_so, i_so_of_i_if, i_so_local
    real(kind=r8_kind)    :: coef
    logical               :: do_sec_ders, do_grads, do_fcts, do_nuc_grads, &
                             do_glob_con, is_charge_fit, is_r2function, &
                             moving_atom, do_moving_sec_ders, do_moving_grads, &
                             do_moving_fcts, do_moving_nuc_grads, alloc_uncont,&
                             do_coulomb_fcts, orbital_read
    ! pointers to the output variables
    real(kind=r8_kind), pointer :: so   (:,:,:)     ! => fcts%o
    real(kind=r8_kind), pointer :: sog  (:,:,:,:)   ! => grads%o
    real(kind=r8_kind), pointer :: sosd (:,:,:,:)   ! => sec_ders%o
    real(kind=r8_kind), pointer :: song (:,:,:,:,:) ! => nuc_grads(i_ma)%o
    real(kind=r8_kind), pointer :: sonsd(:,:,:,:,:) ! => sec_nuc_ders(i_ma)%o
    ! pointers to the intermediates for global contractions
    type(glob_con_intermediate_l_type), pointer :: uncont_orb(:)
    type(glob_con_intermediate_l_type), pointer :: uncont_orb_l
    real(kind=r8_kind),pointer :: uso   (:,:,:)     ! => uncont_orb_l%o
    real(kind=r8_kind),pointer :: usog  (:,:,:,:)   ! => uncont_orb_l%og
    real(kind=r8_kind),pointer :: usosd (:,:,:,:)   ! => uncont_orb_l%osd
    real(kind=r8_kind),pointer :: usong (:,:,:,:,:) ! => uncont_orb_l%ong
    real(kind=r8_kind),pointer :: usonsd(:,:,:,:,:) ! => uncont_orb_l%onsd
    ! basis set information
    type(unique_atom_type)              , pointer :: ua
    integer(kind=i4_kind)                         :: lmax
    type(unique_atom_basis_type)        , pointer :: r2_bas
    type(unique_atom_basis_type)        , pointer :: l_bas(:)
    type(unique_atom_basis_type)        , pointer :: uab
    integer(kind=i4_kind)               , pointer :: orbproj(:,:)
    ! symmetrization information
    type(unique_atom_partner_type)      , pointer :: uap
    type(unique_atom_symadapt_type)     , pointer :: uas
    ! global contraction information
    integer(kind=i4_kind)                         :: N_glob_cons
    type(unique_atom_glob_con_type)     , pointer :: glob_con(:)
    type(unique_atom_glob_con_type)     , pointer :: gci
    integer(kind=i4_kind)                         :: globproj
    ! solid harmonics information
    type(solid_harmonics_array_type)    , pointer :: sha
    type(solid_harmonics_type)          , pointer :: sh
    type(solid_harmonics_grads_type)    , pointer :: shg
    type(solid_harmonics_sec_der_type)  , pointer :: shsd
    real(kind=r8_kind)                  , pointer :: r2(:)
    real(kind=r8_kind)                  , pointer :: x(:), y(:), z(:)
    real(kind=r8_kind)                  , pointer :: shm(:,:)
    real(kind=r8_kind)                  , pointer :: shm_grad(:,:,:)
    real(kind=r8_kind)                  , pointer :: shm_sd(:,:,:)
    ! working arrays
    real(kind=r8_kind), allocatable :: temp(:)
    real(kind=r8_kind), allocatable :: help(:), help_g(:), help_gy(:), &
                                       help_x(:), help_y(:), help_z(:), &
                                       help_xx(:), help_xy(:), help_xz(:), &
                                       help_yy(:), help_yz(:), help_zz(:)
    !------------ Declaration of external functions --------------
    external error_handler
    !------------ Executable code ------------------------------------
    do_sec_ders  = present(sec_ders) .or. present(sec_nuc_ders)
    do_grads     = do_sec_ders .or. present(grads) .or. present(nuc_grads)
    do_fcts      = do_grads .or. present(fcts)
    do_nuc_grads = present(nuc_grads) .or. present(sec_nuc_ders)
    do_coulomb_fcts=.false.

    ! check input parameters and set gridlength
    if (do_sec_ders .and. .not.sec_der_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for second derivatives")
    if (do_grads .and. .not.grads_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for gradients")
    if (do_fcts .and. .not.memory_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for fit functions")
    if (present(N_points)) then
       if (N_points > vec_len) call error_handler( &
            "fit_fct_calculate: N_points larger than vec_len")
       vl = N_points
    else
       vl = vec_len
    endif

    typef=type
    ! check type of fitting functions
    select case (type)
    case ("ch")
       is_charge_fit = .true.
       alloc_uncont = do_glob_con_ch
    case ("xc")
       is_charge_fit = .false.
       alloc_uncont = do_glob_con_xc
    case ("cl")
       is_charge_fit = .true.
       alloc_uncont = do_glob_con_ch
       do_coulomb_fcts=.true.
       typef="ch"
    case default
       call error_handler("fit_fct_calculate: unknown type of fit functions")
    end select

    ! set pointers on the output arrays (if present) else nullify
    ! first for the fit functions and their electronic derivatives
    orbital_read=.false.
    if (present(fcts)) then
       if( orbitalstore_initialized ) then
          call orbitalstore_rw(fcts=fcts)
          if(stored_fcts) orbital_read=.true.
       end if
       if(.not.orbital_read) then
          so => fcts%o
          so = zero
       end if
    else
       nullify(so)
    end if
    if (present(grads)) then
!      print*,'shape grads%o' , shape(grads%o),orbital_read
       if( orbitalstore_initialized ) then
          call orbitalstore_rw(fcts_grads=grads)
          if(stored_fcts_grads) orbital_read=.true.
       end if
       if(.not.orbital_read) then
          sog => grads%o
          sog = zero
       end if
    else
       nullify(sog)
    end if
    if (present(sec_ders)) then
       sosd => sec_ders%o
       sosd = zero
    else
       nullify(sosd)
    end if

    ! Return point for calculated orbitals
    if(orbitalstore_initialized .and. orbital_read) then
       return
    end if

    ! allocation of uncont_orbs_xx for global contractions
    if ( alloc_uncont ) call uncont_fit_fct_allocate(typef, &
         present(fcts),present(grads),present(sec_ders), &
         present(nuc_grads),present(sec_nuc_ders))

    ! allocate working arrays
    allocate( temp(N_max_exponents), stat=status(40) )
    if (status(40) /= 0) call error_handler( &
         "fit_fct_calculate: allocate temp failed")
    if (do_fcts) then
       allocate( help(vl), stat=status(41) )
       if (status(41) /= 0) call error_handler( &
            "fit_fct_calculate: allocate help failed")
    end if
    if (do_grads) then
       allocate( help_x(vl), help_y(vl), help_z(vl), &
                 help_g(vl), help_gy(vl), stat=status(42) )
       if (status(42) /= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{x,y,z} and help_g{,y} failed")
    end if
    if (do_sec_ders) then
       allocate( help_xx(vl), help_xy(vl), help_xz(vl), &
                 help_yy(vl), help_yz(vl), help_zz(vl), stat=status(43) )
       if (status(43) /= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{xx,xy,xz,yy,yz,zz} failed")
    end if

    ! loop over unique atoms
    do i_ua=1,N_unique_atoms
       ua => unique_atoms(i_ua)
       sha => solid_harmonics(i_ua)

       i_ma = ua%moving_atom
       moving_atom = i_ma > 0
       do_moving_sec_ders  = present(sec_ders) .or. &
                             ( moving_atom .and. present(sec_nuc_ders) )
       do_moving_grads     = do_moving_sec_ders .or. present(grads) .or. &
                             ( moving_atom .and. present(nuc_grads) )
       do_moving_fcts      = do_moving_grads .or. present(fcts)
       do_moving_nuc_grads = moving_atom .and. &
                             ( present(nuc_grads) .or. present(sec_nuc_ders) )

       ! now set pointers to the nuclear derivatives of the output arrays
       ! (if present) else nullify
       if (present(nuc_grads) .and. moving_atom) then
          song => nuc_grads(i_ma)%o
          song = zero
       else
          nullify(song)
       end if
       if (present(sec_nuc_ders) .and. moving_atom) then
          sonsd => sec_nuc_ders(i_ma)%o
          sonsd = zero
       else
          nullify(sonsd)
       end if

       ! load solid harmonics and their derivatives on the grid
       do i_ea=1,ua%N_equal_atoms
          ! load solid_harmonics even if do_moving_fcts = .false.
          if (present(grid_points)) call solid_harmonics_calculate( &
               sha%sh(i_ea), grid_points, ua%position(:,i_ea),N_points )
          if (do_moving_grads) call solid_harmonics_calculate_grads( &
               sha%sh(i_ea), sha%shg(i_ea), N_points )
          if (do_moving_sec_ders) call solid_harmonics_calc_sec_der( &
               sha%shg(i_ea), sha%shs(i_ea), N_points )
       end do

       ! set pointer on the basis set information
       if (is_charge_fit) then
          N_glob_cons = ua%N_glob_cons_ch
          lmax = ua%lmax_ch
          r2_bas => ua%r2_ch
          l_bas => ua%l_ch
          orbproj => orbitalprojection_ch
          do_glob_con = N_glob_cons > 0
          if (do_glob_con) then
             glob_con => ua%glob_con_ch
             globproj = orbitalprojection_globcontr_ch(i_ua)
             uncont_orb => uncont_orbs_ch(i_ua)%l
          end if
       else
          N_glob_cons = ua%N_glob_cons_xc
          do_glob_con = N_glob_cons > 0
          lmax = ua%lmax_xc
          r2_bas => ua%r2_xc
          l_bas => ua%l_xc
          orbproj => orbitalprojection_xc
          do_glob_con = N_glob_cons > 0
          if (do_glob_con) then
             glob_con => ua%glob_con_xc
             globproj = orbitalprojection_globcontr_xc(i_ua)
             uncont_orb => uncont_orbs_xc(i_ua)%l
          end if
       end if

       ! start evaluating the fit functions on the grid
       i_ir = get_totalsymmetric_irrep()

       ! loop over all types of basis functions (s,r2,p,d,...)
       do i_bas = -1,lmax
          ! check i_bas to set i_l and is_r2function
          select case (i_bas)
          case (-1) ! s-type
             i_l = 0
             is_r2function = .false.
             uab => l_bas(0)
          case (0) ! r2-type
             i_l = 0
             is_r2function = .true.
             uab => r2_bas
          case default
             i_l = i_bas
             is_r2function = .false.
             uab => l_bas(i_l)
          end select

          ! set pointer to global contraction intermediates
          if (do_glob_con) then
             uncont_orb_l => uncont_orb(i_bas)
             if (associated(so))then
                uso => uncont_orb_l%o
                uso = zero
             else
                nullify(uso)
             end if
             if (associated(sog)) then
                usog => uncont_orb_l%og
                usog = zero
             else
                nullify(usog)
             end if
             if (associated(sosd)) then
                usosd => uncont_orb_l%osd
                usosd = zero
             else
                nullify(usosd)
             end if
             if (associated(song)) then
                usong => uncont_orb_l%ong
                usong = zero
             else
                nullify(usong)
             end if
             if (associated(sonsd)) then
                usonsd => uncont_orb_l%onsd
                usonsd = zero
             else
                nullify(usonsd)
             end if
          end if

          ! first evaluate the radial part
          if (do_moving_grads) &
             temp(:uab%N_exponents) = - ( uab%exponents + uab%exponents )
          do i_ea = 1,ua%N_equal_atoms
             sh => sha%sh(i_ea)
             r2 => sh%r2

             ! calculate the radial part of the primitive orbitals
             if (is_r2function) then
                !              R(r) = r^2 exp(-ar^2)
                ! [1/r d/dr]   R(r) = ( 2 - 2ar^2 ) exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = ( 4a^2r^2 - 8a ) exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   if (do_moving_fcts) then
                     if (do_coulomb_fcts) then
                      call coulomb_pot_of_bas(r2,vl,i_l,.true.,&
                                uab%exponents(i_exp),help)
                      prim_orb_fac(:vl,i_exp,i_ea) =  help
                     else
                      help = exp( - uab%exponents(i_exp) * r2(:vl) )
                      prim_orb_fac(:vl,i_exp,i_ea) = r2(:vl) * help
                     endif
                   end if
                   if (do_moving_grads) then
                      help = help + help
                      prim_grad_fac(:vl,i_exp,i_ea) = help + &
                           temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   end if
                   if (do_moving_sec_ders) then
                      prim_sd_fac(:vl,i_exp,i_ea) = temp(i_exp) * ( help + &
                           prim_grad_fac(:vl,i_exp,i_ea) )
                   end if
                end do
             else
                !              R(r) = N exp(-ar^2)
                ! [1/r d/dr]   R(r) = -2a N exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = 4a^2 N exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   if (do_moving_fcts) then
                     if (do_coulomb_fcts) then
                      !norms already inculded in help
                      call coulomb_pot_of_bas(r2,vl,i_l,.false.,&
                                uab%exponents(i_exp),help)
                      prim_orb_fac(:vl,i_exp,i_ea) =  help
                     else
                      prim_orb_fac(:vl,i_exp,i_ea) = &
                           uab%norms(i_exp) * exp( - uab%exponents(i_exp) * r2(:vl) )
                     endif
                   end if
                   if (do_moving_grads) then
                      prim_grad_fac(:vl,i_exp,i_ea) = &
                           temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   end if
                   if (do_moving_sec_ders) then
                      prim_sd_fac(:vl,i_exp,i_ea) = &
                           temp(i_exp) * prim_grad_fac(:vl,i_exp,i_ea)
                   end if
                end do
             end if

             ! calculate the radial part of the local contractions
             ! [1/r d/dr]^n C_i(r) = Sum(a) coeff(a,i) [1/r d/dr]^n R_a(r)
             do i_c = 1,uab%N_contracted_fcts ! contractions
                if (do_moving_fcts) then
                   cont_orb_fac(:vl,i_c,i_ea) = zero
                   do i_exp = 1,uab%N_exponents ! contributing exponents
                      coef = uab%contractions(i_exp,i_c)
                      cont_orb_fac(:vl,i_c,i_ea) = cont_orb_fac(:vl,i_c,i_ea)+&
                           coef * prim_orb_fac(:vl,i_exp,i_ea)
                   end do ! contributing exponents
                end if
                if (do_moving_grads) then
                   cont_grad_fac(:vl,i_c,i_ea) = zero
                   do i_exp = 1,uab%N_exponents ! contributing exponents
                      coef = uab%contractions(i_exp,i_c)
                      cont_grad_fac(:vl,i_c,i_ea)=cont_grad_fac(:vl,i_c,i_ea)+&
                           coef * prim_grad_fac(:vl,i_exp,i_ea)
                   end do! contributing exponents
                end if
                if (do_moving_sec_ders) then
                   cont_sd_fac(:vl,i_c,i_ea) = zero
                   do i_exp = 1,uab%N_exponents ! contributing exponents
                      coef = uab%contractions(i_exp,i_c)
                      cont_sd_fac(:vl,i_c,i_ea) = cont_sd_fac(:vl,i_c,i_ea) + &
                           coef * prim_sd_fac(:vl,i_exp,i_ea)
                   end do! contributing exponents
                end if
             end do! contractions

          end do! equal atoms

          ! symmetry adaption:
          ! for all independent fcts,
          !         uncontracted exponents and contractions,
          !         gridpoints:
          !   sum over equal atoms and m
          !
          uap => ua%symadapt_partner(i_ir,i_l)
          i_step = uap%N_independent_fcts
          i_sa_of_i_if = orbproj(i_bas,i_ua)
          i_so_of_i_if = i_sa_of_i_if - orbproj(-1,i_ua) + 1
          do i_if = 1,uap%N_independent_fcts ! independent fcts
             uas => uap%symadapt(i_if,1)

             if (i_l == 0) then
                ! special case of s-type functions (just one independet fct)
                ! i_fu = i_ea, m = 0 , coef = 1.0, Y_lm(r) = 1.0
                do i_ea = 1,ua%N_equal_atoms ! contributing functions
                   if (do_moving_grads) then
                      x => sha%sh(i_ea)%l(1)%m(:,2)
                      y => sha%sh(i_ea)%l(1)%m(:,3)
                      z => sha%sh(i_ea)%l(1)%m(:,1)
                   end if
                   if (do_glob_con) then
                      ! accumulate primitives in uncont_orb
                      ! and skip primitives in fcts, grads, ...
                      i_sa = 1
                      i_so = 1
                      call combine_rad_and_y00(uab%N_exponents, 1.0_r8_kind, &
                           prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                           uso, usog, usosd, usong, usonsd )
                      i_sa = i_sa_of_i_if + uab%N_uncontracted_fcts
                      i_so = i_so_of_i_if + uab%N_uncontracted_fcts
                   else
                      ! accumulate primitives in fcts, grads, ...
                      i_sa = i_sa_of_i_if
                      i_so = i_so_of_i_if
                      call combine_rad_and_y00(uab%N_uncontracted_fcts, 1.0_r8_kind, &
                           prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                           so, sog, sosd, song, sonsd )
                   end if
                   ! accumulate contractions in fcts, grads, ...
                   call combine_rad_and_y00(uab%N_contracted_fcts, 1.0_r8_kind, &
                        cont_orb_fac, cont_grad_fac, cont_sd_fac, &
                        so, sog, sosd, song, sonsd )

                end do ! contributing functions

             else ! i_l > 0
                ! the standard case including the Y_lm contributions

                do i_fu = 1,uas%N_fcts ! contributing functions
                   i_ea = uas%I_equal_atom(i_fu)
                   m = uas%m(i_fu)
                   coef = uas%c(i_fu)
                   if (do_moving_fcts) then
                      sh => sha%sh(i_ea)
                      shm => sh%l(i_l)%m
                   end if
                   if (do_moving_grads) then
                      x => sh%l(1)%m(:,2)
                      y => sh%l(1)%m(:,3)
                      z => sh%l(1)%m(:,1)
                      shg => sha%shg(i_ea)
                      shm_grad => shg%l(i_l)%m
                   end if
                   if (do_moving_sec_ders) then
                      shsd => sha%shs(i_ea)
                      shm_sd => shsd%l(i_l)%m
                   end if

                   if (do_glob_con) then
                      ! accumulate primitives in uncont_orb
                      ! and skip primitives in fcts, grads, ...
                      i_sa = 1
                      i_so = 1
                      call combine_rad_and_ylm(uab%N_exponents, &
                           coef, 1, i_if, &
                           prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                           uso, usog, usosd, usong, usonsd )
                      i_sa = i_sa_of_i_if + uab%N_uncontracted_fcts * i_step
                      i_so = i_so_of_i_if + uab%N_uncontracted_fcts * i_step
                   else
                      ! accumulate primitives in fcts, grads, ...
                      i_sa = i_sa_of_i_if
                      i_so = i_so_of_i_if
                      call combine_rad_and_ylm(uab%N_uncontracted_fcts, &
                           coef, i_step, 1, &
                           prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                           so, sog, sosd, song, sonsd )
                   end if
                   ! accumulate contractions in fcts, grads, ...
                   call combine_rad_and_ylm(uab%N_contracted_fcts, &
                        coef, i_step, 1, &
                        cont_orb_fac, cont_grad_fac, cont_sd_fac, &
                        so, sog, sosd, song, sonsd )
!                  print*,'grads%o grads%o',sum(sog),shape(sog)

                end do ! contributing functions

             end if

             i_sa_of_i_if = i_sa_of_i_if + 1
             i_so_of_i_if = i_so_of_i_if + 1
          end do ! independent fcts

          i_so_local = i_so - i_step ! index of last local orbital of current
                                     ! basis type within actual unique atom

          if (do_glob_con) then
             ! copy uncontracted orbitals from intermediate into result
             i_sa = orbproj(i_bas,i_ua)
             i_so = i_sa - orbproj(-1,i_ua) + 1
             do i_exp = 1,uab%N_uncontracted_fcts ! uncontracted exponents
                do i_if = 1,uap%N_independent_fcts ! independent fcts
                   if (associated(so)) then
                      so(:vl,i_sa,1) = uso(:vl,i_exp,i_if)
                   end if
                   if (associated(sog)) then
                      sog(:vl,:,i_sa,1) = usog(:vl,:,i_exp,i_if)
                   end if
                   if (associated(sosd)) then
                      sosd(:vl,:,i_sa,1) = usosd(:vl,:,i_exp,i_if)
                   end if
                   if (associated(song)) then
                      song(:vl,:,:,i_so,1) = usong(:vl,:,:,i_exp,i_if)
                   end if
                   if (associated(sonsd)) then
                      sonsd(:vl,:,:,i_so,1) = usonsd(:vl,:,:,i_exp,i_if)
                   end if
                   i_sa = i_sa + 1
                   i_so = i_so + 1
                end do! independent fcts
             end do! uncontracted exponents
          end if

       end do ! i_bas
       ! now i_so_local holds the index of the last local orbital within
       ! entire set of orbitals of the current unique atom

       if (do_glob_con) then
          ! calculate globally contracted functions
          ! and store them at the end of the result arrays
          i_sa = globproj
          i_so = i_so_local + 1
          do i_gc = 1,N_glob_cons ! global contractions
             gci => glob_con(i_gc)
             do i_gcf = 1,gci%N_contributing_fcts ! contributing_fcts
                i_bas = gci%l(i_gcf)
                i_exp = gci%index_exp(i_gcf)
                i_if  = gci%index_ind_fct(i_gcf)
                coef  = gci%coefs(i_gcf)

                uncont_orb_l => uncont_orb(i_bas)
                if (associated(so)) then
                   so(:vl,i_sa,1) = so(:vl,i_sa,1) + &
                       coef * uncont_orb_l%o(:vl,i_exp,i_if)
                end if
                if (associated(sog)) then
                   sog(:vl,:,i_sa,1) = sog(:vl,:,i_sa,1) + &
                        coef * uncont_orb_l%og(:vl,:,i_exp,i_if)
                end if
                if (associated(sosd)) then
                   sosd(:vl,:,i_sa,1) = sosd(:vl,:,i_sa,1) + &
                        coef * uncont_orb_l%osd(:vl,:,i_exp,i_if)
                end if
                if (associated(song)) then
                    song(:vl,:,:,i_so,1) = song(:vl,:,:,i_so,1) + &
                        coef * uncont_orb_l%ong(:vl,:,:,i_exp,i_if)
                end if
                if (associated(sonsd)) then
                   sonsd(:vl,:,:,i_so,1) = sonsd(:vl,:,:,i_so,1) + &
                        coef * uncont_orb_l%onsd(:vl,:,:,i_exp,i_if)
                end if
             end do ! contributing_fcts
             i_sa = i_sa + 1
             i_so = i_so + 1
          end do ! global contractions
       end if

    end do ! unique atom

    if( options_orbitals_on_file() ) then
       if( present(fcts)  ) call orbitalstore_rw(fcts=fcts)
       if( present(grads) ) call orbitalstore_rw(fcts_grads=grads)
    end if

    ! deallocate working arrays
    if (do_sec_ders) then
       deallocate( help_xx, help_xy, help_xz, &
                   help_yy, help_yz, help_zz, stat=status(43) )
       if (status(43) /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help_{xx,xy,xz,yy,yz,zz} failed")
       status(43)=1
    end if
    if (do_grads) then
       deallocate( help_x, help_y, help_z, &
                   help_g, help_gy, stat=status(42) )
       if (status(42) /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help_{x,y,z} and help_g{,y} failed")
       status(42)=1
    end if
    if (do_fcts) then
       deallocate( help, stat=status(41) )
       if (status(41) /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help failed")
       status(41)=1
    end if
    deallocate( temp, stat=status(40) )
    if (status(40) /= 0) call error_handler( &
         "fit_fct_calculate: deallocate temp failed")
    status(40)=1

    ! deallocation of uncont_orbs_xx
    if ( alloc_uncont ) call uncont_fit_fct_free(typef, &
         present(fcts),present(grads),present(sec_ders), &
         present(nuc_grads),present(sec_nuc_ders))
  contains

    subroutine combine_rad_and_ylm(N_fcts, coeff, i_incr, i_sel, &
         rad, rad_g, rad_sd, o, og, osd, ong, onsd)
      ! Purpose : evaluation of the primitive Gaussians and their various
      !           derivatives on a set of grid points
      ! ----- Declaration of formal parameters -----------------------------
      ! alternating control parameter
      integer(kind=i4_kind), intent(in) :: N_fcts
      real   (kind=r8_kind), intent(in) :: coeff
      integer(kind=i4_kind), intent(in) :: i_incr
      integer(kind=i4_kind), intent(in) :: i_sel
      ! the radial functions to process
      real(kind=r8_kind), intent(in) :: rad   (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_g (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_sd(:,:,:)
      ! the result types (in arbitrary occurance)
      real(kind=r8_kind), optional, pointer :: o   (:,:,:)
      real(kind=r8_kind), optional, pointer :: og  (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: osd (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: ong (:,:,:,:,:)
      real(kind=r8_kind), optional, pointer :: onsd(:,:,:,:,:)
      ! ----- Declaration of local variables -------------------------------
      integer(kind=i4_kind) :: i_fct
      !---------------------------------------------------------------------
      !               phi(r) = C R(r) Y_lm(r)
      !        d/dr_i phi(r) = C [1/r d/dr] R(r) r_i Y_lm(r)
      !                      + C R(r) d/dr_i Y_lm(r)
      ! d/dr_i d/dr_j phi(r) = C [1/r d/dr]^2 R(r) r_j r_i Y_lm(r)
      !                      + C [1/r d/dr] R(r) delta_ij Y_lm(r)
      !                      + C [1/r d/dr] R(r) r_i d/dr_j Y_lm(r)
      !                      + C [1/r d/dr] R(r) r_j d/dr_i Y_lm(r)
      !                      + C R(r) d/dr_j d/dr_i Y_lm(r)
      !------------ Executable code --------------------------------
      do i_fct = 1,N_fcts

         if (do_fcts) then
            ! help(r) := C * R(r)
            help = coeff * rad(:vl,i_fct,i_ea)
            if (associated(o)) then
               o(:vl,i_sa,i_sel) = o(:vl,i_sa,i_sel) + &
                    help * shm(:vl,m)
            end if
         end if

         if (do_grads) then
            ! help_g(r)   := C * [1/r d/dr] R(r)
            ! help_gy(r)  := C * [1/r d/dr] R(r) * Y_lm(r)
            ! help_r_i(r) := d/dr_i phi(r)
            help_g  = coeff * rad_g(:vl,i_fct,i_ea)
            help_gy = help_g * shm(:vl,m)
            help_x  = help_gy * x(:vl) + &
                      help * shm_grad(:vl,1,m)
            help_y  = help_gy * y(:vl) + &
                      help * shm_grad(:vl,2,m)
            help_z  = help_gy * z(:vl) + &
                      help * shm_grad(:vl,3,m)
            if (associated(og)) then
               og(:vl,1,i_sa,i_sel) = og(:vl,1,i_sa,i_sel) + help_x
               og(:vl,2,i_sa,i_sel) = og(:vl,2,i_sa,i_sel) + help_y
               og(:vl,3,i_sa,i_sel) = og(:vl,3,i_sa,i_sel) + help_z
            end if
            if (associated(ong)) then
               ong(:vl,1,i_ea,i_so,i_sel) = ong(:vl,1,i_ea,i_so,i_sel) + help_x
               ong(:vl,2,i_ea,i_so,i_sel) = ong(:vl,2,i_ea,i_so,i_sel) + help_y
               ong(:vl,3,i_ea,i_so,i_sel) = ong(:vl,3,i_ea,i_so,i_sel) + help_z
            end if
         end if

         if (do_sec_ders) then
            ! help_r_ir_j(r) := d/dr_i d/dr_j phi(r)
            help_xx = help * shm_sd(:vl,XX,m) + help_gy
            help_xy = help * shm_sd(:vl,XY,m)
            help_xz = help * shm_sd(:vl,XZ,m)
            help_yy = help * shm_sd(:vl,YY,m) + help_gy
            help_yz = help * shm_sd(:vl,YZ,m)
            help_zz = help * shm_sd(:vl,ZZ,m) + help_gy
            ! help_r_i(r)    := C * [1/r d/dr] R(r) * r_i
            help_x  = help_g * x(:vl)
            help_y  = help_g * y(:vl)
            help_z  = help_g * z(:vl)
            help_xx = help_xx + help_x * shm_grad(:vl,1,m) * 2.0_r8_kind
            help_xy = help_xy + help_x * shm_grad(:vl,2,m) &
                              + help_y * shm_grad(:vl,1,m)
            help_xz = help_xz + help_x * shm_grad(:vl,3,m) &
                              + help_z * shm_grad(:vl,1,m)
            help_yy = help_yy + help_y * shm_grad(:vl,2,m) * 2.0_r8_kind
            help_yz = help_yz + help_y * shm_grad(:vl,3,m) &
                              + help_z * shm_grad(:vl,2,m)
            help_zz = help_zz + help_z * shm_grad(:vl,3,m) * 2.0_r8_kind
            ! help_gy(r)     := C * [1/r d/dr]^2 R(r) * Y_lm(r)
            ! help_r_i(r)    := C * [1/r d/dr]^2 R(r) * r_i * Y_lm(r)
            help_gy = coeff * rad_sd(:vl,i_fct,i_ea) * shm(:vl,m)
            help_x  = help_gy * x(:vl)
            help_y  = help_gy * y(:vl)
            help_z  = help_gy * z(:vl)
            help_xx = help_xx + help_x * x(:vl)
            help_xy = help_xy + help_x * y(:vl)
            help_xz = help_xz + help_x * z(:vl)
            help_yy = help_yy + help_y * y(:vl)
            help_yz = help_yz + help_y * z(:vl)
            help_zz = help_zz + help_z * z(:vl)
            if (associated(osd)) then
               osd(:vl,XX,i_sa,i_sel) = osd(:vl,XX,i_sa,i_sel) + help_xx
               osd(:vl,XY,i_sa,i_sel) = osd(:vl,XY,i_sa,i_sel) + help_xy
               osd(:vl,XZ,i_sa,i_sel) = osd(:vl,XZ,i_sa,i_sel) + help_xz
               osd(:vl,YY,i_sa,i_sel) = osd(:vl,YY,i_sa,i_sel) + help_yy
               osd(:vl,YZ,i_sa,i_sel) = osd(:vl,YZ,i_sa,i_sel) + help_yz
               osd(:vl,ZZ,i_sa,i_sel) = osd(:vl,ZZ,i_sa,i_sel) + help_zz
!!$               osd(:vl,2,1,i_sa,i_sel) = osd(:vl,2,1,i_sa,i_sel) + help_xy
!!$               osd(:vl,3,1,i_sa,i_sel) = osd(:vl,3,1,i_sa,i_sel) + help_xz
!!$               osd(:vl,3,2,i_sa,i_sel) = osd(:vl,3,2,i_sa,i_sel) + help_yz
            end if
            if (associated(onsd)) then
               onsd(:vl,1,i_ea,i_so,i_sel)=onsd(:vl,1,i_ea,i_so,i_sel)+help_xx
               onsd(:vl,2,i_ea,i_so,i_sel)=onsd(:vl,2,i_ea,i_so,i_sel)+help_xy
               onsd(:vl,3,i_ea,i_so,i_sel)=onsd(:vl,3,i_ea,i_so,i_sel)+help_xz
               onsd(:vl,4,i_ea,i_so,i_sel)=onsd(:vl,4,i_ea,i_so,i_sel)+help_yy
               onsd(:vl,5,i_ea,i_so,i_sel)=onsd(:vl,5,i_ea,i_so,i_sel)+help_yz
               onsd(:vl,6,i_ea,i_so,i_sel)=onsd(:vl,6,i_ea,i_so,i_sel)+help_zz
            end if
         end if

         i_sa = i_sa + i_incr
         i_so = i_so + i_incr
      end do

    end subroutine combine_rad_and_ylm

    subroutine combine_rad_and_y00(N_fcts, coeff, &
         rad, rad_g, rad_sd, o, og, osd, ong, onsd)
      ! Purpose : evaluation of the primitive Gaussians and their various
      !           derivatives on a set of grid points
      ! ----- Declaration of formal parameters -----------------------------
      ! alternating control parameter
      integer(kind=i4_kind), intent(in) :: N_fcts
      real   (kind=r8_kind), intent(in) :: coeff
      ! the radial functions to process
      real(kind=r8_kind), intent(in) :: rad   (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_g (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_sd(:,:,:)
      ! the result types (in arbitrary occurance)
      real(kind=r8_kind), optional, pointer :: o   (:,:,:)
      real(kind=r8_kind), optional, pointer :: og  (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: osd (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: ong (:,:,:,:,:)
      real(kind=r8_kind), optional, pointer :: onsd(:,:,:,:,:)
      ! ----- Declaration of local variables -------------------------------
      integer(kind=i4_kind) :: i_fct
      !---------------------------------------------------------------------
      !               phi(r) = C R(r)
      !        d/dr_i phi(r) = C [1/r d/dr] R(r) r_i
      ! d/dr_i d/dr_j phi(r) = C [1/r d/dr]^2 R(r) r_j r_i
      !                      + C [1/r d/dr] R(r) delta_ij
      !------------ Executable code --------------------------------
      do i_fct = 1,N_fcts

         if (do_fcts) then
            if (associated(o)) then
               o(:vl,i_sa,1) = o(:vl,i_sa,1) + &
                               coeff * rad(:vl,i_fct,i_ea)
            end if
         end if

         if (do_grads) then
            ! help_g(r)   := C * [1/r d/dr] R(r)
            ! help_r_i(r) := d/dr_i phi(r)
            help_g  = coeff * rad_g(:vl,i_fct,i_ea)
            help_x  = help_g * x(:vl)
            help_y  = help_g * y(:vl)
            help_z  = help_g * z(:vl)
            if (associated(og)) then
               og(:vl,1,i_sa,1) = og(:vl,1,i_sa,1) + help_x
               og(:vl,2,i_sa,1) = og(:vl,2,i_sa,1) + help_y
               og(:vl,3,i_sa,1) = og(:vl,3,i_sa,1) + help_z
            end if
            if (associated(ong)) then
               ong(:vl,1,i_ea,i_so,1) = ong(:vl,1,i_ea,i_so,1) + help_x
               ong(:vl,2,i_ea,i_so,1) = ong(:vl,2,i_ea,i_so,1) + help_y
               ong(:vl,3,i_ea,i_so,1) = ong(:vl,3,i_ea,i_so,1) + help_z
            end if
         end if

         if (do_sec_ders) then
            ! help(r)        := C * [1/r d/dr]^2 R(r)
            ! help_r_i(r)    := C * [1/r d/dr]^2 R(r) * r_i
            ! help_r_ir_j(r) := d/dr_i d/dr_j phi(r)
            help    = coeff * rad_sd(:vl,i_fct,i_ea)
            help_x  = help * x(:vl)
            help_y  = help * y(:vl)
            help_z  = help * z(:vl)
            help_xx = help_x * x(:vl) + help_g
            help_xy = help_x * y(:vl)
            help_xz = help_x * z(:vl)
            help_yy = help_y * y(:vl) + help_g
            help_yz = help_y * z(:vl)
            help_zz = help_z * z(:vl) + help_g
            if (associated(osd)) then
               osd(:vl,XX,i_sa,1) = osd(:vl,XX,i_sa,1) + help_xx
               osd(:vl,XY,i_sa,1) = osd(:vl,XY,i_sa,1) + help_xy
               osd(:vl,XZ,i_sa,1) = osd(:vl,XZ,i_sa,1) + help_xz
               osd(:vl,YY,i_sa,1) = osd(:vl,YY,i_sa,1) + help_yy
               osd(:vl,YZ,i_sa,1) = osd(:vl,YZ,i_sa,1) + help_yz
               osd(:vl,ZZ,i_sa,1) = osd(:vl,ZZ,i_sa,1) + help_zz
!!$               osd(:vl,2,1,i_sa,1) = osd(:vl,2,1,i_sa,1) + help_xy
!!$               osd(:vl,3,1,i_sa,1) = osd(:vl,3,1,i_sa,1) + help_xz
!!$               osd(:vl,3,2,i_sa,1) = osd(:vl,3,2,i_sa,1) + help_yz
            end if
            if (associated(onsd)) then
               onsd(:vl,1,i_ea,i_so,1) = onsd(:vl,1,i_ea,i_so,1) + help_xx
               onsd(:vl,2,i_ea,i_so,1) = onsd(:vl,2,i_ea,i_so,1) + help_xy
               onsd(:vl,3,i_ea,i_so,1) = onsd(:vl,3,i_ea,i_so,1) + help_xz
               onsd(:vl,4,i_ea,i_so,1) = onsd(:vl,4,i_ea,i_so,1) + help_yy
               onsd(:vl,5,i_ea,i_so,1) = onsd(:vl,5,i_ea,i_so,1) + help_yz
               onsd(:vl,6,i_ea,i_so,1) = onsd(:vl,6,i_ea,i_so,1) + help_zz
            end if
         end if

         i_sa = i_sa + 1
         i_so = i_so + 1
      end do

    end subroutine combine_rad_and_y00
!AG #endif
  end subroutine fit_fct_calculate
  !*************************************************************

#ifdef WITH_CORE_DENS
  ! the subroutines for core_orbital_calculation with fitting functions
    subroutine fit_core_calculate(grid_points,N_points,&
       fcts, grads, sec_ders, nuc_grads, sec_nuc_ders)
    !  Purpose: calculates the values of the fit core functions
    !  and their various derivates on a set of grid points
    use orbitalstore_module, only: core_orbital_type, core_orbital_gradient_type, &
        core_orbital_sec_der_type, core_orbital_nuc_gradient_type,  &
        core_orbital_nuc_sec_der_type
    use unique_atom_module ! desription of unique atoms, must be calculated before
    use orbitalprojection_module ! projection of indices to metaindex used in scf
    use solid_harmonics_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in), optional :: grid_points(:,:)
    ! give this argument only at first call for given grid_points.
    ! Solid harmonics are evaluated when this argument is present and
    ! remain stored
    integer, intent(in), optional :: N_points
    ! default: vector_length from orbital_setup
    ! either "ch" for charge fit functions or "xc" for exchange fit functions
    type(core_orbital_type)                 , target, optional :: fcts
    type(core_orbital_gradient_type)        , target, optional :: grads
    type(core_orbital_sec_der_type)         , target, optional :: sec_ders
    type(core_orbital_nuc_gradient_type)    , target, optional :: nuc_grads(:)
    type(core_orbital_nuc_sec_der_type)     , target, optional :: sec_nuc_ders(:)
    ! the various types of fit function arrays to load (all intent(out))
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind), parameter :: zero = 0.0_r8_kind
    integer(kind=i4_kind) :: vl, i_ma, i_ua, i_ea, &
                             i_bas, i_l, i_exp, &
                             i_step, i_sa, i_sa_of_i_if, &
                             i_so, i_so_of_i_if
    logical               :: do_sec_ders, do_grads, do_fcts, do_nuc_grads, &
                             is_r2function, &
                             moving_atom, do_moving_sec_ders, do_moving_grads, &
                             do_moving_fcts, do_moving_nuc_grads
    ! pointers to the output variables
    real(kind=r8_kind), pointer :: so   (:,:)     ! => fcts%o
    real(kind=r8_kind), pointer :: sog  (:,:,:)   ! => grads%o
    real(kind=r8_kind), pointer :: sosd (:,:,:)   ! => sec_ders%o
    real(kind=r8_kind), pointer :: song (:,:,:,:) ! => nuc_grads(i_ma)%o
    real(kind=r8_kind), pointer :: sonsd(:,:,:,:) ! => sec_nuc_ders(i_ma)%o
    ! basis set information
    type(unique_atom_type)              , pointer :: ua
    type(unique_atom_atomic_dens_type)        , pointer :: r2_bas
    type(unique_atom_atomic_dens_type)        , pointer :: s_bas
    type(unique_atom_atomic_dens_type)        , pointer :: uab
    integer(kind=i4_kind)               , pointer :: orbproj(:,:)
    ! solid harmonics information
    type(solid_harmonics_array_type)    , pointer :: sha
    type(solid_harmonics_type)          , pointer :: sh
!!$    type(solid_harmonics_grads_type)    , pointer :: shg
!!$    type(solid_harmonics_sec_der_type)  , pointer :: shsd
    real(kind=r8_kind)                  , pointer :: r2(:)
    real(kind=r8_kind)                  , pointer :: x(:), y(:), z(:)
!!$    real(kind=r8_kind)                  , pointer :: shm(:,:)
!!$    real(kind=r8_kind)                  , pointer :: shm_grad(:,:,:)
!!$    real(kind=r8_kind)                  , pointer :: shm_sd(:,:,:)
    ! working arrays
    real(kind=r8_kind), allocatable :: temp(:)
    real(kind=r8_kind), allocatable :: help(:), help_g(:), help_gy(:), &
                                       help_x(:), help_y(:), help_z(:), &
                                       help_xx(:), help_xy(:), help_xz(:), &
                                       help_yy(:), help_yz(:), help_zz(:)
    !------------ Declaration of external functions --------------
    external error_handler
    !------------ Executable code ------------------------------------
    do_sec_ders  = present(sec_ders) .or. present(sec_nuc_ders)
    do_grads     = do_sec_ders .or. present(grads) .or. present(nuc_grads)
    do_fcts      = do_grads .or. present(fcts)
    do_nuc_grads = present(nuc_grads) .or. present(sec_nuc_ders)

    ! check input parameters and set gridlength
    if (do_sec_ders .and. .not.sec_der_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for second derivatives")
    if (do_grads .and. .not.grads_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for gradients")
    if (do_fcts .and. .not.memory_allocated) call error_handler( &
         "fit_fct_calculate: no memory allocated for fit functions")
    if (present(N_points)) then
       if (N_points > vec_len) call error_handler( &
            "fit_fct_calculate: N_points larger than vec_len")
       vl = N_points
    else
       vl = vec_len
    endif

    ! allocate working arrays
    allocate( temp(N_max_exponents), stat=status(44) )
    if (status(44) /= 0) call error_handler( &
         "fit_fct_calculate: allocate temp failed")
    if (do_fcts) then
       allocate( help(vl), stat=status(45) )
       if (status(45) /= 0) call error_handler( &
            "fit_fct_calculate: allocate help failed")
    end if
    if (do_grads) then
       allocate( help_x(vl), help_y(vl), help_z(vl), &
                 help_g(vl), help_gy(vl), stat=status(46) )
       if (status(46)/= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{x,y,z} and help_g{,y} failed")
    end if
    if (do_sec_ders) then
       allocate( help_xx(vl), help_xy(vl), help_xz(vl), &
                 help_yy(vl), help_yz(vl), help_zz(vl), stat=status(47) )
       if (status(47) /= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{xx,xy,xz,yy,yz,zz} failed")
    end if

    ! set pointers on the output arrays (if present) else nullify
    ! first for the fit functions and their electronic derivatives
    if (present(fcts)) then
       so => fcts%o
       so = zero
    else
       nullify(so)
    end if
    if (present(grads)) then
       sog => grads%o
       sog = zero
    else
       nullify(sog)
    end if
    if (present(sec_ders)) then
       sosd => sec_ders%o
       sosd = zero
    else
       nullify(sosd)
    end if

    ! loop over unique atoms
    do i_ua=1,N_unique_atoms
       ua => unique_atoms(i_ua)
       sha => solid_harmonics(i_ua)

       i_ma = ua%moving_atom
       moving_atom = i_ma > 0
       do_moving_sec_ders  = present(sec_ders) .or. &
                             ( moving_atom .and. present(sec_nuc_ders) )
       do_moving_grads     = do_moving_sec_ders .or. present(grads) .or. &
                             ( moving_atom .and. present(nuc_grads) )
       do_moving_fcts      = do_moving_grads .or. present(fcts)
       do_moving_nuc_grads = moving_atom .and. &
                             ( present(nuc_grads) .or. present(sec_nuc_ders) )

       ! now set pointers to the nuclear derivatives of the output arrays
       ! (if present) else nullify
       if (present(nuc_grads) .and. moving_atom) then
          song => nuc_grads(i_ma)%o
          song = zero
       else
          nullify(song)
       end if
       if (present(sec_nuc_ders) .and. moving_atom) then
          sonsd => sec_nuc_ders(i_ma)%o
          sonsd = zero
       else
          nullify(sonsd)
       end if

       ! load solid harmonics and their derivatives on the grid
       do i_ea=1,ua%N_equal_atoms
          ! load solid_harmonics even if do_moving_fcts = .false.
          if (present(grid_points)) call solid_harmonics_calculate( &
               sha%sh(i_ea), grid_points, ua%position(:,i_ea),N_points )
          if (do_moving_grads) call solid_harmonics_calculate_grads( &
               sha%sh(i_ea), sha%shg(i_ea), N_points )
          if (do_moving_sec_ders) call solid_harmonics_calc_sec_der( &
               sha%shg(i_ea), sha%shs(i_ea), N_points )
       end do

       ! set pointer on the basis set information

       orbproj => orbitalprojection_cd
       r2_bas => ua%r2_core
       s_bas  => ua%s_core

       ! loop over all types of basis functions (s,r2)
       do i_bas = -1, 0
          ! check i_bas to set i_l and is_r2function
          select case (i_bas)
          case (-1) ! s-type
             i_l = 0
             is_r2function = .false.
             uab => s_bas
          case (0) ! r2-type
             i_l = 0
             is_r2function = .true.
             uab => r2_bas
          case default
          end select

          ! first evaluate the radial part
          if (do_moving_grads) &
             temp(:uab%N_exponents) = - ( uab%exponents + uab%exponents )
          do i_ea = 1,ua%N_equal_atoms
             sh => sha%sh(i_ea)
             r2 => sh%r2

             ! calculate the radial part of the primitive orbitals
             if (is_r2function) then
                !              R(r) = r^2 exp(-ar^2)
                ! [1/r d/dr]   R(r) = ( 2 - 2ar^2 ) exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = ( 4a^2r^2 - 8a ) exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   if (do_moving_fcts) then
                      help = exp( - uab%exponents(i_exp) * r2(:vl) )
                      prim_orb_fac(:vl,i_exp,i_ea) = r2(:vl) * help
                   end if
                   if (do_moving_grads) then
                      help = help + help
                      prim_grad_fac(:vl,i_exp,i_ea) = help + &
                           temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   end if
                   if (do_moving_sec_ders) then
                      prim_sd_fac(:vl,i_exp,i_ea) = temp(i_exp) * ( help + &
                           prim_grad_fac(:vl,i_exp,i_ea) )
                   end if
                end do
             else
                !              R(r) = N exp(-ar^2)
                ! [1/r d/dr]   R(r) = -2a N exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = 4a^2 N exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   if (do_moving_fcts) then
                      prim_orb_fac(:vl,i_exp,i_ea) = &
                                 exp( - uab%exponents(i_exp) * r2(:vl) )
                   end if
                   if (do_moving_grads) then
                      prim_grad_fac(:vl,i_exp,i_ea) = &
                           temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   end if
                   if (do_moving_sec_ders) then
                      prim_sd_fac(:vl,i_exp,i_ea) = &
                           temp(i_exp) * prim_grad_fac(:vl,i_exp,i_ea)
                   end if
                end do
             end if

          end do! equal atoms

          i_step = 1
          i_sa_of_i_if = orbproj(i_bas,i_ua)
          i_so_of_i_if = i_sa_of_i_if - orbproj(-1,i_ua) + 1

                ! special case of s-type functions (just one independet fct)
                ! i_fu = i_ea, m = 0 , coef = 1.0, Y_lm(r) = 1.0
                do i_ea = 1,ua%N_equal_atoms ! contributing functions
                   if (do_moving_grads) then
                      x => sha%sh(i_ea)%l(1)%m(:,2)
                      y => sha%sh(i_ea)%l(1)%m(:,3)
                      z => sha%sh(i_ea)%l(1)%m(:,1)
                   end if

                      ! accumulate primitives in fcts, grads, ...
                      i_sa = i_sa_of_i_if
                      i_so = i_so_of_i_if
                      call combine_rad_and_y00_core(uab%N_exponents,&
                           prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                           so, sog, sosd, song, sonsd )

                end do ! contributing functions


           !  i_sa_of_i_if = i_sa_of_i_if + 1
           !  i_so_of_i_if = i_so_of_i_if + 1

           !i_so_local = i_so - i_step ! index of last local orbital of current
           !                          ! basis type within actual unique atom

       end do ! i_bas
       ! now i_so_local holds the index of the last local orbital within
       ! entire set of orbitals of the current unique atom

    end do ! unique atom
    if (do_sec_ders) then
       MEMLOG(-size(help_xx)*6)
       deallocate( help_xx, help_xy, help_xz, &
                   help_yy, help_yz, help_zz, stat=status(47) )
       ASSERT(status(47).eq.0)
       status(47)=1
    end if
    if (do_grads) then
       MEMLOG(-size(help_x)*3-size(help_g)-size(help_gy))
       deallocate( help_x, help_y, help_z, &
                   help_g, help_gy, stat=status(46) )
       ASSERT(status(46).eq.0)
       status(46)=1
    end if
    if (do_fcts) then
       MEMLOG(-size(help))
       deallocate( help, stat=status(45) )
      ASSERT(status(45).eq.0)
      status(45)=1
    end if
    MEMLOG(-size(temp))
    deallocate( temp, stat=status(44) )
    ASSERT(status(44).eq.0)
    status(44)=1

  contains

    subroutine combine_rad_and_y00_core(N_fcts, &
         rad, rad_g, rad_sd, o, og, osd, ong, onsd)
      ! Purpose : evaluation of the primitive Gaussians and their various
      !           derivatives on a set of grid points
      ! ----- Declaration of formal parameters -----------------------------
      ! alternating control parameter
      integer(kind=i4_kind), intent(in) :: N_fcts
      ! the radial functions to process
      real(kind=r8_kind), intent(in) :: rad   (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_g (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_sd(:,:,:)
      ! the result types (in arbitrary occurance)
      real(kind=r8_kind), optional, pointer :: o   (:,:)
      real(kind=r8_kind), optional, pointer :: og  (:,:,:)
      real(kind=r8_kind), optional, pointer :: osd (:,:,:)
      real(kind=r8_kind), optional, pointer :: ong (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: onsd(:,:,:,:)
      ! ----- Declaration of local variables -------------------------------
      integer(kind=i4_kind) :: i_fct
      !---------------------------------------------------------------------
      !               phi(r) = C R(r)
      !        d/dr_i phi(r) = C [1/r d/dr] R(r) r_i
      ! d/dr_i d/dr_j phi(r) = C [1/r d/dr]^2 R(r) r_j r_i
      !                      + C [1/r d/dr] R(r) delta_ij
      !------------ Executable code --------------------------------
      do i_fct = 1,N_fcts

         if (do_fcts) then
            if (associated(o)) then
               o(:vl,i_sa) = o(:vl,i_sa) + &
                                rad(:vl,i_fct,i_ea)
            end if
         end if

         if (do_grads) then
            ! help_g(r)   := C * [1/r d/dr] R(r)
            ! help_r_i(r) := d/dr_i phi(r)
            help_g  =  rad_g(:vl,i_fct,i_ea)
            help_x  = help_g * x(:vl)
            help_y  = help_g * y(:vl)
            help_z  = help_g * z(:vl)
            if (associated(og)) then
               og(:vl,1,i_sa) = og(:vl,1,i_sa) + help_x
               og(:vl,2,i_sa) = og(:vl,2,i_sa) + help_y
               og(:vl,3,i_sa) = og(:vl,3,i_sa) + help_z
            end if
            if (associated(ong)) then
               ong(:vl,1,i_ea,i_so) = ong(:vl,1,i_ea,i_so) + help_x
               ong(:vl,2,i_ea,i_so) = ong(:vl,2,i_ea,i_so) + help_y
               ong(:vl,3,i_ea,i_so) = ong(:vl,3,i_ea,i_so) + help_z
            end if
         end if

         if (do_sec_ders) then
            ! help(r)        := C * [1/r d/dr]^2 R(r)
            ! help_r_i(r)    := C * [1/r d/dr]^2 R(r) * r_i
            ! help_r_ir_j(r) := d/dr_i d/dr_j phi(r)
            help    =  rad_sd(:vl,i_fct,i_ea)
            help_x  = help * x(:vl)
            help_y  = help * y(:vl)
            help_z  = help * z(:vl)
            help_xx = help_x * x(:vl) + help_g
            help_xy = help_x * y(:vl)
            help_xz = help_x * z(:vl)
            help_yy = help_y * y(:vl) + help_g
            help_yz = help_y * z(:vl)
            help_zz = help_z * z(:vl) + help_g
            if (associated(osd)) then
               osd(:vl,XX,i_sa) = osd(:vl,XX,i_sa) + help_xx
               osd(:vl,XY,i_sa) = osd(:vl,XY,i_sa) + help_xy
               osd(:vl,XZ,i_sa) = osd(:vl,XZ,i_sa) + help_xz
               osd(:vl,YY,i_sa) = osd(:vl,YY,i_sa) + help_yy
               osd(:vl,YZ,i_sa) = osd(:vl,YZ,i_sa) + help_yz
               osd(:vl,ZZ,i_sa) = osd(:vl,Zz,i_sa) + help_zz
!!$               osd(:vl,2,1,i_sa) = osd(:vl,2,1,i_sa) + help_xy
!!$               osd(:vl,3,1,i_sa) = osd(:vl,3,1,i_sa) + help_xz
!!$               osd(:vl,3,2,i_sa) = osd(:vl,3,2,i_sa) + help_yz
            end if
            if (associated(onsd)) then
               onsd(:vl,1,i_ea,i_so) = onsd(:vl,1,i_ea,i_so) + help_xx
               onsd(:vl,2,i_ea,i_so) = onsd(:vl,2,i_ea,i_so) + help_xy
               onsd(:vl,3,i_ea,i_so) = onsd(:vl,3,i_ea,i_so) + help_xz
               onsd(:vl,4,i_ea,i_so) = onsd(:vl,4,i_ea,i_so) + help_yy
               onsd(:vl,5,i_ea,i_so) = onsd(:vl,5,i_ea,i_so) + help_yz
               onsd(:vl,6,i_ea,i_so) = onsd(:vl,6,i_ea,i_so) + help_zz
            end if
         end if

         i_sa = i_sa + 1
         i_so = i_so + 1
      end do

    end subroutine combine_rad_and_y00_core

  end subroutine fit_core_calculate
#endif
  !*************************************************************
!end module orbital_module

!.............................................................................
!Evaluation of the fit functions
!===============================
!xx = ch or xc
!
!>>the unique atoms
!  type(unique_atom_type) :: unique_atoms( N_unique_atoms )
!  ua => unique_atoms(i_ua)
!  >>the equal atoms
!    real(kind=r8_kind) :: ua%position( 3 , ua%N_equal_atoms )
!  >>the fitting function basis sets
!    type(unique_atom_basis_type) :: ua%l_xx( 0:ua%lmax_xx ) [s,p,...]
!    type(unique_atom_basis_type) :: ua%r2_xx                [r^2]
!           ua%l_xx(0)     : i_bas = -1  (s-type)
!    uab => ua%r2_xx       : i_bas =  0  (r^2-type)
!           ua%l_xx(i_bas) : i_bas >= 1  (l>0-type)
!    >>the exponents of the primitives
!      real(kind=r8_kind) :: uab%exponents( uab%N_exponents )
!    >>the normalization factor of the primitives
!      real(kind=r8_kind) :: uab%norms( uab%N_exponents )
!      [ (2a)^l for l>0-type functions and 1 for s- and r^2-type functions ]
!    >>the local contraction coefficients
!      real(kind=r8_kind) :: uab%contractions( uab%N_exponents , &
!                                              uab%N_contracted_fcts )
!  >>the symmetry adaptation information for each l-type (excluding r^2)
!    i_ir = get_totalsymmetric_irrep()
!    type(unique_atom_partner_type) :: ua%symadapt_partner( i_ir:i_ir , &
!                                                           0:ua%lmax_all )
!    uap => ua%symadapt_partner(i_ir,i_l)  [s,p,...]
!    >>the l-specific symmetry adaptation
!      type(unique_atom_symadapt_type) :: uap%symadapt( &
!                                         uap%N_independent_fcts , 1 )
!      uas => uap%symadapt(i_if,1)  [only one-dim. total symm. IRREPs]
!      >>the contributing equal atoms
!        integer(kind=i4_kind) :: uas%I_equal_atom( uas%N_fcts )
!      >>the contributing m-quantum numbers
!        integer(kind=i4_kind) :: uas%m( uas%N_fcts )
!      >>the symmetry adaptation coefficients
!        real(kind=r8_kind) :: uas%c( uas%N_fcts )
!  >>the normalization factors of the fitting function [including r^2]
!    type(unique_atom_renorm_indfct_type) :: ua%renormaliation_partner_xx( &
!                                            -1:ua%lmax_xx )
!    rp => ua%renormaliation_partner_xx(i_bas)  [s,r2,p,...]
!    >>the basis-type specific normalization
!      type(unique_atom_renorm_type) :: rp%renorm( uap%N_independent_fcts )
!      r => rp%renorm(i_if)
!      >>normalization factor of the uncontracted and primitive orbitals
!        real(kind=r8_kind) :: r%c_exp( uab%N_exponents )
!      >>normalization factor of the contracted orbitals
!        real(kind=r8_kind) :: r%c_contr( uab%N_contracted_fcts )
!  >>the global contraction scheme
!    type(unique_atom_glob_con_type) :: ua%glob_con_xx( ua%N_glob_cons_xx )
!    gci => ua%glob_con_xx(i_gc)
!    >>the contributing orbital type [s,r^2,p,...]
!      integer(kind=i4_kind) :: gci%l( gci%N_contributing_fcts )
!    >>the contributing l-specififc symmetry adaptations
!      integer(kind=i4_kind) :: gci%index_ind_fct( gci%N_contributing_fcts )
!    >>the contributing primitive(exponent)
!      integer(kind=i4_kind) :: gci%index_exp( gci%N_contributing_fcts )
!    >>the global contraction coefficients
!      real(kind=r8_kind) :: gci%coefs( gci%N_contributing_fcts )
!>>the solid harmonics on the grid
!  type(solid_harmonics_array_type) :: solid_harmonics( N_unique_atoms )
!  sha => solid_harmonics(i_ua)
!  >>the individual solid harmonic Y_lm(r-R) on the grid
!    type(solid_harmonics_type) :: sha%sh( ua%N_equal_atoms )
!    sh => sha%sh(i_ea)
!    >>the values of (r-R)^2 on the grid
!      real(kind=r8_kind) :: sh%r2( grid_length )
!      r2 => sh%r2
!    >>the values of r-R on the grid
!      x => sh%l(1)%m(:,2)
!      y => sh%l(1)%m(:,3)
!      z => sh%l(1)%m(:,1)
!    >>the Y_lm terms of order l
!      type(solid_harmonics_ltype) :: sh%l( 0:l_dim ) ! l_dim >= ua%lmax_xx
!      shm => sh%l(i_l)%m
!      >>the values on the grid
!        real(kind=r8_kind) :: shm( grid_length , m_dim )
!  >>the d/dr Y_lm(r-R) contributions
!    type(solid_harmonics_grads_type) :: sha%shg( ua%N_equal_atoms )
!    shg => sha%shg(i_ea)
!    >>the d/dr Y_lm terms of order l
!      type(solid_harmonics_grad_ltype) :: shg%l( 0:l_dim )
!      shm_grad => shg%l(i_l)%m
!      >>the values on the grid
!        real(kind=r8_kind) :: shm_grad( grid_length , 3 , m_dim)
!  >>the d^2/dr^2 Y_lm(r-R) contributions
!    type(solid_harmonics_sec_der_type) :: sha%shs( ua%N_equal_atoms )
!    shsd => sha%shs(i_ea)
!    >>the d^2/dr^2 terms of order l
!      type(solid_harmonics_sec_der_ltype) :: shsd%l( 0:l_dim )
!      shm_sd => shsd%l(i_l)%m
!      >>the values on the grid
!        real(kind=r8_kind) :: shm_sd( grid_length , 3 , 3 , m_dim )
!>>the radial part R(r) of the primitives on the grid
!  real(kind=r8_kind) :: prim_orb_fac( vec_length , &
!                                      N_max_exponents , &
!                                      N_max_equal_atoms )
!>>the first derivative (1/r d/dr) R(r) of the radial part
!  real(kind=r8_kind) :: prim_grad_fac( vec_length , &
!                                       N_max_exponents , &
!                                       N_max_equal_atoms )
!>>the second derivative (1/r d/dr)^2 R(r) of the radial part
!  real(kind=r8_kind) :: prim_sd_fac( vec_length , &
!                                     N_max_exponents , &
!                                     N_max_equal_atoms )
!
!>>the radial part C(r) of the contractions on the grid
!  real(kind=r8_kind) :: cont_orb_fac( vec_length , &
!                                      N_max_exponents , &
!                                      N_max_equal_atoms )
!>>the first derivative (1/r d/dr) C(r) of the radial part
!  real(kind=r8_kind) :: cont_grad_fac( vec_length , &
!                                       N_max_exponents , &
!                                       N_max_equal_atoms )
!>>the second derivative (1/r d/dr)^2 C(r) of the radial part
!  real(kind=r8_kind) :: cont_sd_fac( vec_length , &
!                                     N_max_exponents , &
!                                     N_max_equal_atoms )
!>>the primitves on the grid as intermediates for the global contractions
!  type(glob_con_intermediate_type) :: uncont_orbs_xx( N_unique_atoms )
!  uncont_orb => uncont_orbs_xx(i_ua)   [in orbital_calculate]
!  uncont_orb => uncont_orbs_xx(i_ua)%l
!  >>the l-contributions to the global contractions
!    type(glob_con_intermediate_l_type) :: uncont_orb(-1:ua%lmax_xx)
!    uncont_orb_l => uncont_orb%l(i_bas)%oo  [in orbital calculate]
!    uncont_orb_l => uncont_orb(i_bas)       [s,r^2,p,...]
!    uso    => uncont_orb_l%o
!    usog   => uncont_orb_l%og
!    usosd  => uncont_orb_l%osd
!    usong  => uncont_orb_l%ong
!    usonsd => uncont_orb_l%onsd
!    >>the values on the grid
!      real(kind=r8_kind) :: uso( grid_length , &
!                                 uab%N_exponents , &
!                                 uap%N_independent_fcts )
!    >>the values of the gradients on the grid
!      real(kind=r8_kind) :: usog( grid_length , 1:3 , &
!                                  uab%N_exponents , &
!                                  uap%N_independent_fcts )
!    >>the values of the second derivatives on the grid
!      real(kind=r8_kind) :: usosd( grid_length , 1:3 , 1:3 , &
!                                   uab%N_exponents , &
!                                   uap%N_independent_fcts )
!    >>the values of the nuclear gradients on the grid
!      real(kind=r8_kind) :: usong( grid_length , 1:3 , &
!                                   ua%N_equal_atoms , &
!                                   uab%N_exponents , &
!                                   uap%N_independent_fcts )
!    >>the values of the second nuclear derivatives on the grid
!      real(kind=r8_kind) :: usonsd( grid_length , 1:6 , &
!                                    ua%N_equal_atoms , &
!                                    uab%N_exponents , &
!                                    uap%N_independent_fcts )
!>>the local fitting function index association
!  integer(kind=i4_kind) :: orbitalprojection_xx( -1:N_max_l_xx, & [s,r^2,p..]
!                                                 N_unique_atoms )
!>>the global contraction index association
!  integer(kind=i4_kind) :: orbitalprojection_globcontr_xx( N_unique_atoms )
!>>the internal order of the fitting functions
!  i_sa = orbitalprojection_xx(i_bas,i_ua)
!  DO i_unc = 1,uab%N_uncontracted_fcts
!     DO i_if = 1, uap%N_independent_fcts
!        pos(i_ua,i_bas,i_if,i_unc) = i_sa
!        i_sa = i_sa+1
!     END DO
!  END DO
!  DO i_con = 1,uab%N_contracted_fcts
!     DO i_if = 1, uap%N_independent_fcts
!        pos(i_ua,i_bas,i_if,i_con) = i_sa
!        i_sa = i_sa + 1
!     END DO
!  END DO
!  i_sa = orbitalprojection_globcontr_xx(i_ua)
!  DO i_gc = 1,ua%N_glob_cons_xx
!     pos(i_ua,i_gc) = i_sa
!     i_sa = i_sa + 1
!  END DO
!  >>the internal order of orbitalprojection_xx(i_bas,i_au)
!  >>and orbitalprojection_globcontr_xx(i_ua)
!    i_sa = 1
!    DO i_ua = 1,N_unique_atoms
!       DO i_bas = -1,ua%lmax_xx
!          orbitalprojection_xx(i_bas,i_ua) = i_sa
!          i_sa = i_sa + ( uab%N_uncontracted_fcts + uab%N_contracted_fcts ) &
!                      * uap%N_independent_fcts
!       END DO
!    END DO
!    DO i_ua = 1,N_unique_atoms
!       orbitalprojection_globcontr_xx(i_ua) = i_sa
!       i_sa = i_sa + ua%N_glob_cons_xx
!    END DO
!>>the symmetrized fitting functions f_k on the grid
!  type(orbital_type) :: fcts_xx
!  so => fcts_xx%o
!  >>the values on the grid
!    real(kind=r8_kind) :: so( grid_length , fit_coeff_n_xx() , 1:1 )
!>>the first drivative d/dr_i f_k on the grid
!  type(orbital_gradient_type) :: grads_xx
!  sog => grads_xx%o
!  >>the values on the grid
!    real(kind=r8_kind) :: sog( grid_length , 3 , fit_coeff_n_xx() , 1:1 )
!>>the second drivative d/dr_i d/dr_j f_k on the grid
!  type(orbital_sec_der_type) :: sec_ders_xx
!  sosd => sec_ders_xx%o
!  >>the values on the grid
!    real(kind=r8_kind) :: sosd( grid_length , 3 , 3 , fit_coeff_n_xx() , 1:1 )
!>>the nuclear drivatives [-d/dR_i] f_k on the grid
!  type(orbital_nuclear_gradient_type) :: nuc_grads_xx( &
!                                         N_moving_unique_atoms )
!  song => nuc_grads_xx%(i_ma)%o
!  >>the values on the grid
!    real(kind=r8_kind) :: song( grid_length , 3 , &
!                                ua%N_equal_atoms, N_fitfcts_on_ua , 1:1 )
!>>the second nuclear drivatives [-d/dR_i][-d/dR_j] f_k on the grid
!>>the mixed nuclear derivates   [-d/dR_i]  d/dr    f_k on the grid
!  type(orbital_nuclear_sec_der_type) :: sec_nuc_ders_xx( &
!                                        N_moving_unique_atoms )
!  sonsd => sec_nuc_ders_xx(i_ma)%o
!  >>the values on the grid
!    real(kind=r8_kind) :: sonsd( grid_length , 6 , &
!                                 ua%N_equal_atoms, N_fitfcts_on_ua , 1:1 )
!                                 1: xx  2: xy  3: xz  4: yy  5: yz  6: zz
!.............................................................................


!!!!! MF testing

  subroutine coulomb_pot_of_bas(r2,N_points,i_l,is_r2function,expo&
                ,prim_int_array)
    use gamma_module
    !  Purpose: calculates the coulomb potential of the primitive
    !charge fit functionson a set of grid points
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), pointer :: r2(:)
    real(kind=r8_kind), intent(in) :: expo
    integer(kind=i4_kind), intent(in) :: N_points,i_l
    real(kind=r8_kind), intent(inout) :: prim_int_array(:)
    logical, intent(in) :: is_r2function
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: allo_stat
    ! working arrays
    real(kind=r8_kind), allocatable :: gam_arg(:),gam_fct(:,:)
    real(kind=r8_kind):: pfac,pfac1
    real(kind=r8_kind),parameter :: pi=3.14159265358979324_r8_kind, two=2.0_r8_kind
    !------------ Declaration of external functions --------------
    external error_handler
    !------------ Executable code ------------------------------------

    prim_int_array(1:N_points) = 0.0_r8_kind

    allocate( gam_arg(N_points) ,stat=allo_stat)
           if(allo_stat/=0) call error_handler(&
                "coul_pot_of_bas: allocation of gam (1) failed")
    gam_arg(1:N_points)=expo*r2(1:N_points)

    if(is_r2function) then
       allocate( gam_fct(N_points,2) ,stat=allo_stat)
           if(allo_stat/=0) call error_handler(&
                "coul_pot_of_bas: allocation of gam (2) failed")
       pfac1=two*pi/expo
       pfac=pfac1/expo
       gam_fct(:,1:2) = gamma(2,gam_arg)
       prim_int_array(:)= pfac*gam_fct(:,1)+pfac1*r2(:)*gam_fct(:,2)
    else
       allocate( gam_fct(N_points,i_l+1) ,stat=allo_stat)
           if(allo_stat/=0) call error_handler(&
                "coul_pot_of_bas: allocation of gam (3) failed")
       pfac=two**(i_l+1) * pi * expo**(i_l-1)
       gam_fct(:,1:i_l+1) = gamma(i_l+1,gam_arg(1:N_points))
       prim_int_array(:)= pfac*gam_fct(1:N_points,i_l+1)
    endif

    deallocate(gam_arg,gam_fct,stat=allo_stat)
    if(allo_stat/=0) call error_handler("rspace_coul: dealloc of ints failed")


end subroutine coulomb_pot_of_bas

end module orbital_module

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
module unique_atom_module
  !
  !  SEE UNIQUE_ATOM_METHODS ...
  !
  !  Database   for  Coordinates   and   Coefficients  of   primitive,
  !  contracted and symmetry adapted basis functions.  Subroutines for
  !  reading and  writing this database  and packing and  unpacking it
  !  are  also defined.   Memory allocation  is automatically  done by
  !  these subroutines.  Subroutines  doing memory allocation are made
  !  public for data that will be calculated externally.
  !
  !  References: publisher document: orbital_calculation
  !
  !
  !  Author: TB
  !  Date: 18.07.95
  !
  !== Interrupt of public interface of module ========================
  !
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
  ! Author: FN
  ! Date:   8/96
  ! Description: Shift the calculation of the fitbasis
  !              dimensions to the fit_coeff_module due
  !              to cross-dependencies (use statements)
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   11/96
  ! Description:
  !   moving the reading of global contractions in input file:
  !   they should be read in together with the rest of basis set
  !   for a given unique atom.
  !   Deleted unique_atom_glob_cons_read/write.
  !   Reading/Writing of global contractions now in unique_atom_read/write.
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: Uwe Birkenheuer
  ! Date:   June 1998
  ! Description: Moving_Unique_Atom concept introduced
  !----
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
! Author: AH
! Date:   11/98
! Description: add pseudopotential basis set
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!---------------------------------------------------------------------

!------------ Modules used --------------------------------------
#include "def.h"
use type_module  ! type specification parameters
use datatype, only : arrmat3
use uatom_symmadapt, only:&
     & unique_atom_symadapt_type => symadapt_type,&
     & unique_atom_sa_int_type   => sa_int_type,&
     & unique_atom_partner_type  => partner_type

implicit none
private ! by default, everything is private
!== Interrupt end of public interface of module =================

!------------ Declaration of types ------------------------------

public unique_atom_symadapt_type

public unique_atom_sa_int_type

public unique_atom_partner_type

type, public ::  unique_atom_glob_con_type
   ! Describes the global contraction of primitive to contracted basis
   ! functions  of charge  fit basis  or exchange  fit  basis.  Theses
   ! contractions allow summation  of symmetry-adapted functions build
   ! from  primitive  functions  (one  exponent, no  local  l-specific
   ! contraction) over all independent functions  of all l and r2. One
   ! global-contracted basis function is  obtained by summing over all
   ! basis  functions  described  by  the indices  (l,  index_ind_fct,
   ! index_exp) given here.
   integer :: N_contributing_fcts
       ! Number of exponents per contracted basis functions
   integer, allocatable :: l(:) ! l(N_contributing_fcts)
       ! l of fcts contributing, -1 means l=0, 0 means r2 function
   integer, allocatable :: index_ind_fct(:) ! index_ind_fct(N_contributing_fcts)
       ! Indices of independent functions contributing
   integer, allocatable :: index_exp(:) ! index_exp(N_contributing_fcts)
       ! Indices of exponents contributing
   real(r8_kind), allocatable :: coefs(:) ! coefs(N_contributing_fcts)
       ! contraction coefficients
end type unique_atom_glob_con_type


type, public ::  unique_atom_basis_type
   ! Angular quantum  number dependent information  for either orbital
   ! basis or charge fit basis or exchange fit basis.  information for
   ! r**2 fitbasis of charge and exchange
   integer                           :: N_exponents
       ! Number of exponents
   integer                           :: N_uncontracted_fcts
       ! Number of uncontracted basis functions
   real(r8_kind), allocatable :: exponents(:) ! exponents(N_exponents)
       ! exponents of primitive basis functions
       ! it is assumed that the first N_uncontracted_exp are used directly
       ! without contractions
   real(r8_kind), allocatable :: norms(:) ! norms(N_exponents)
       ! norms of primitive basis functions
       ! ignored for r**2 type
   integer                           :: N_contracted_fcts
       ! Number of contracted basis functions
   real(r8_kind), allocatable :: contractions(:,:)
       ! contractions(N_exponents, N_contracted_fcts)
       ! contraction coefficients for the exponents

   !
   ! Linear combinations of basis functions, e.g. density expansion:
   !
   real(r8_kind), allocatable   :: functions(:, :) ! (nbas, nfcts)
end type unique_atom_basis_type

type, public :: unique_atom_pseudopot_type
   ! Angular   quantum  number   dependent   information  for   pseudo
   ! potentials
   integer                           :: N_exponents
       ! Number of exponents
   real(r8_kind), allocatable :: exponents(:) ! exponents(N_exponents)
       ! exponents of the Gaussian functions
   real(r8_kind), allocatable :: coefficients(:) ! coefficients(N_exponents)
       ! contraction coefficient of the Gaussian functions
   integer(i4_kind), allocatable :: powers(:) ! powers(N_exponents)
       ! power of the pre-factor r^n of the Gaussian functions
end type unique_atom_pseudopot_type

#ifdef WITH_CORE_DENS
type, public :: unique_atom_atomic_dens_type
   ! Information about fitting basis sets for atomic densities just s-
   ! and r^2-type functions (no charge fit renormalization!)
   integer                     :: N_exponents = -1
       ! Number of exponents
   real(r8_kind), allocatable :: exponents(:) ! exponents(N_exponents)
       ! exponents of primitive basis functions
   real(r8_kind), allocatable :: contractions(:) ! contractions(N_exponents)
       ! contraction coefficients for the atomic density
   integer(i4_kind), allocatable :: linked_fitfcts(:) ! linked_fitfcts(N_exponents)
       ! link to the corresponding charge fitting function
end type unique_atom_atomic_dens_type
#endif

type, public :: unique_atom_type
   !
   ! Describes   all   properties   of   primitive,   contracted   and
   ! symmetry-adapted basis functions of an unique atom.
   !
   character(len=12) :: name    ! name of unique atom

   real(r8_kind) :: Z           ! charge of nucleus

   ! The core charge.  Set when reading the basis set/pseudo potential
   ! description. Should remain zero if ECPs are not used.
   real (r8_kind) :: ZC = 0.0

   integer :: nuclear_spin ! the spin of the nucleus is nuclear_spin *
                           ! 1/2

   ! Nuclear magnetic moment in nuclear magnetons:
   real(r8_kind) :: nuclear_magnetic_moment

   real(r8_kind) :: nuclear_radius ! Nuclear radius in fm

   integer :: N_equal_atoms     ! Number of partners of unique atom

   integer :: moving_atom ! = 0 if the position of the unique atom is
   ! fixed => no gradients and the index of the unique atom among all
   ! moving unique atoms.

   integer :: lmax_ob     ! maximal angular momentum of orbital basis

   integer :: lmax_pseudo ! maximal angualr momentum of
                          ! pseudopotentials

   integer :: lmax_ch   ! maximal angular momentum of charge fit basis

   integer :: lmax_xc ! maximal angular momentum of exchange fit basis

   integer :: lmax_all          ! max(lmax_orbit, lmax_ch, lmax_xc)

   logical :: core_density_use ! some potentials with core density
                               ! dependent

   real(r8_kind), allocatable :: position(:, :)
   ! position (3, N_equal_atoms) positions of partners ("equal atoms")
   ! of unique atom

   real(r8_kind) :: position_first_ea(3) ! positions of first partner
   ! ("equal atom") of unique atom

   type(unique_atom_basis_type), allocatable :: l_ob(:)
   ! l_orbit (0:lmax_ob) angular quantum number dependent information
   ! for orbital basis

   type(unique_atom_basis_type) :: r2_ch ! r2 information for charge
                                         ! fit basis

   type(unique_atom_basis_type), allocatable :: l_ch(:)
   ! l_ch (0:lmax_ch) angular quantum number dependent information for
   ! charge fit basis

#ifdef WITH_CORE_DENS
   type(unique_atom_atomic_dens_type) :: s_core ! s-type information
   ! for the core charge fit basis

   type(unique_atom_atomic_dens_type) :: r2_core ! r2-type information
   ! for the core charge fit basis
#endif

   type(unique_atom_basis_type) :: r2_xc ! r2 information for exchange
                                         ! fit basis

   type(unique_atom_basis_type), allocatable :: l_xc(:) ! l_xc(0:lmax_xc)
   ! angular quantum number dependent information for exchange fit
   ! basis

   type(unique_atom_pseudopot_type), allocatable :: l_pseudopot(:)
   ! l_pseudopot(0:lmax_ob+1) for pseudopotentials

   type(unique_atom_partner_type), allocatable :: symadapt_partner(:, :)
   ! symadapt_partner(N_irreps, 0:lmax_sym) information for symmetry
   ! adaption of any basis

   type(unique_atom_partner_type), allocatable :: symadapt_spor_partner(:, :)
   ! symadapt_partner(N_irreps, 0:lmax_sym) information for symmetry
   ! adaption of any basis of the spin orbitals

   !
   ! These are rarely used in practical calculations:
   !
   integer :: N_glob_cons_ch ! number of global (over all l)
   ! contractions of charge fitfunctions if set to zero only local
   ! contractions are performed

   integer :: N_glob_cons_xc ! number of global (over all l)
   ! contractions of exchange fitfunctions if set to zero only local
   ! contractions are performed

   type(unique_atom_glob_con_type), allocatable :: glob_con_ch(:)
   ! glob_con_ch(N_glob_cons_ch) desribes global contraction (over all
   ! independent functions of all l and r2) of charge fitfunctions

   type(unique_atom_glob_con_type), allocatable :: glob_con_xc(:)
   ! glob_con_ch(N_glob_cons_xc) desribes global contraction (over all
   ! independent functions of all l and r2) of exchange fitfunctions
end type unique_atom_type


type,public :: unique_atom_symequivatoms_type
   ! contains necessary information to perform integrals
   ! according to A.Goerling , Dissertation, page 34,
   ! Eq. (2.73)
   integer(i4_kind)     :: n_equiv_dist
   ! number of different sym.equiv. distances for a
   ! unique atom
   integer(i4_kind), allocatable :: index_eq_atom2(:)
   ! indices of those atoms with a unique sym.equiv
   ! distance
   integer(i4_kind), allocatable :: weight(:)
   ! weights for the above atoms
end type unique_atom_symequivatoms_type


!------------ Declaration of public variables -------------------
integer(i4_kind), public           :: N_unique_atoms = 0
   ! Number of unique atoms
integer(i4_kind), public           :: N_moving_unique_atoms
   ! Number of non-fixed unique atoms
type(unique_atom_type), allocatable, target, public :: unique_atoms(:)
   ! unique_atom(N_unique_atoms)
type(unique_atom_type), allocatable, target, public :: unique_atoms_eperef(:)
   ! unique_atoms(N_unique_atoms_eperef)
integer(i4_kind), public  :: &
     unique_atom_lmax_all = -1, &
     unique_atom_lmax_ob = -1, &
     unique_atom_lmax_ch = -1, &
     unique_atom_lmax_xc = -1, &
     unique_atom_lmax_pseudo=-1

   ! maximal l values for all atoms
integer(i4_kind), public           :: unique_atom_iwork
   ! iwork as read from gxfile. This is a control parameter for the
   ! geometry optimizer

logical, public                         :: pseudopot_present  = .false. ! true if any atom is a pp-atom
logical, public                         :: core_density_setup = .false. ! true to enforce use of core density

type(unique_atom_symequivatoms_type), allocatable, target, public :: unique_atom_symequiv(:,:)
   ! unique_atom_symequiv(N_unique_atoms,N_unique_atoms)

type(arrmat3), allocatable, target, public :: unique_atom_grad_info(:)
   ! unique_atom_grad_info(N_moving_unique_atoms)%m(n_indep,3,n_equals)
   !  secound dimension for magnetical quantum numbers:  1: m=2, 2: m=3, 3: m=1
   ! Equivalent to variable
   !  unique_atom(i)%symapdapt_partner(irrep,l)%symadapt(n_indep,n_partner)%
   !    N_fcts
   !    I_equal_atom
   !    c ...
   ! with i_irrep = 1, l = 1
   ! but better suited for the gradient part where
   ! we loop over equal atoms. (-> better performance)
integer(i4_kind), allocatable, public :: moving_unique_atom_index(:)
   ! moving_unique_atom_index(N_moving_unique_atoms)
   ! the i-th moving unique atom is actually the
   ! MOVING_UNIQUE_ATOM_INDEX(i)-th unique atom

!
! These are constructors for unique_atom_basis_type:
!
public :: unique_atom_make_empty_basis!()
public :: unique_atom_make_basis!(L, exponents, contractions)

!
! These  two  encode legacy  normalization  factors  for the  gaussian
! functions:
!
public :: unique_atom_basis_norm_ob!(L, alpha)
public :: unique_atom_basis_norm_ch!(L, alpha)

  !===================================================================
  ! End of public interface of module
  !===================================================================

contains

  !
  ! Maybe in the future this module will also contain constructors for
  ! the types it declares ...
  !

  function unique_atom_make_empty_basis() result(basis)
    implicit none
    type(unique_atom_basis_type) :: basis
    ! *** end of interface ***

    real(r8_kind) :: empty1(0), empty2(0, 0)

    basis = unique_atom_basis_type(0, 0, empty1, empty1, 0, empty2, empty2)
  end function unique_atom_make_empty_basis

  function unique_atom_make_basis(L, exponents, contractions) result(basis)
    !
    ! Constructor   for  a   typical  basis,   that  does   not  have
    ! uncontracted  exponents (or rather  represents them  as trivial
    ! contractions,  if necessary).  The  legacy representation  of a
    ! shell basis assumes that  contraction coefficients are those of
    ! a  linear combination over  NORMALIZED gaussians.   These norms
    ! are precomputed  and stored separately,  this is also  the only
    ! reason we pass angular momentum to this function.
    !
    implicit none
    integer(i4_kind), intent(in) :: L
    real(r8_kind), intent(in) :: exponents(:), contractions(:, :)
    type(unique_atom_basis_type) :: basis
    ! *** end of interface ***

    integer(i4_kind) :: ne, nc
    real(r8_kind) :: empty2(0, 0)

    ne = size(exponents)
    nc = size(contractions, 2)
    ASSERT(ne==size(contractions, 1))

    basis = unique_atom_basis_type(ne, 0, exponents, &
         unique_atom_basis_norm_ob(L, exponents), nc, contractions, empty2)
  end function unique_atom_make_basis

  elemental function unique_atom_basis_norm_ob(L, alpha) result(norm)
    !
    ! Legacy normalization factor for orbital gaussians.
    !
    implicit none
    integer(i4_kind), intent(in) :: L
    real(r8_kind), intent(in) :: alpha
    real(r8_kind) :: norm
    ! *** end of interface ***

    ! FIXME: use constants:
    real(r8_kind), parameter :: pi = 3.14159265358979324_r8_kind
    real(r8_kind) :: df

    ! double factorial, convert to real:
    df = dfac(L)

    norm = (2 * alpha / pi)**0.75_r8_kind * (2 * sqrt(alpha))**L / sqrt(df)
  end function unique_atom_basis_norm_ob

  elemental function unique_atom_basis_norm_ch(L, alpha) result(norm)
    !
    ! Legacy normalization factor for fit gaussians.
    !
    integer, intent(in) :: L
    real(r8_kind), intent(in) :: alpha
    real(r8_kind) :: norm
    ! *** end of interface ***

    if (L >= 0) then
       ! Zeroth  power  of  a  number  is one,  this  holds  also  for
       ! s-exponents:
       norm = (alpha + alpha)**L
    else
       ! Special case for r2-exponents:
       norm = 1.0
    endif
  end function unique_atom_basis_norm_ch

  elemental function dfac(k) result(n)
    !
    ! dfac(k) = 1 * 3 * ... * (2k - 1)
    !
    ! k |   0   1   2   3
    ! --+----------------
    ! n |   1   1   3  15
    !
    ! FIXME: log2(dfac(10)) = 29.2, log2(dfac(17)) = 62.5
    !
    implicit none
    integer(i4_kind), intent(in) :: k
    integer(i4_kind) :: n
    ! *** end of interface ***

    integer :: i

    n = 1
    do i = 1, k
       n = n * (2 * i - 1)
    enddo
  end function dfac

end module unique_atom_module

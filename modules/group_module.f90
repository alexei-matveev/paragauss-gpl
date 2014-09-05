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
module group_module
  !-------------------------------------------------------------------
  !
  !  Purpose: performs operations between symmetry elements.
  !
  !  Author: MM
  !  Date: 07.07.1996
  !
  !-------------------------------------------------------------------
!== Interrupt of public interface of module =====================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author:
  ! Date:
  ! Description:
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables in the module until explicit deallocation
  private         ! by default, all names are private
!== Interrupt end of public interface of module =================


  !------------ Declaration of types ---------------------------------

  type, public ::  group_symm_el
     ! symmetry operation of the pointgroup, given in quaternion form
     real(kind=r8_kind)          :: quaternion(5)
     ! quaternion parameters
     character(len=8)            :: name
     ! name of symmetry operation
     integer(kind=i4_kind)       :: inverse
     ! inverse element
     integer(kind=i4_kind)       :: inverse_sign
     ! sign of the inverse element
     ! inverse_sign = -1 : inverse = -inverse for spinor rotations
     ! inverse_sign = +1 : inverse = +inverse for spinor rotations
     integer(kind=i4_kind)       :: klass
     ! class of element
  end type group_symm_el


  type, public :: group_class
     ! classes of the group
     integer(kind=i4_kind), allocatable :: conjug_el(:)
     ! group elements conjug_el(number)
     integer(kind=i4_kind)           :: number
     ! number of class elements
     integer(kind=i4_kind)           :: inverse
     ! inverse element
     logical                         :: ambivalent
     ! is class ambivalent?
     logical                         :: ambivalent_proj
     ! is class ambivalent in rep group?
     logical                         :: regular
     ! is class regular?
  end type group_class

  type, public :: group_irrep
     ! irreps of the group

     !
     ! variables which identify the irrep
     !
     character(len=8)                :: label
     ! irrep label
     real(kind=r8_kind)              :: eigen_value
     ! eigenvalue of CSCO of irrep
     integer(kind=i4_kind)           :: dimension
     ! dimension of the irrep
     integer(kind=i4_kind)           :: time_dimension
     ! time inversion invariant dimension of the irrep
     logical                         :: pseudo
     ! .true. if the irrep is a pseudo-2D irrep
     real(kind=r8_kind), allocatable :: characters(:)
     ! characters of the irrep: characters(group_num_cl)

     !
     ! the following variables are used to store the
     ! representation matrices of the irrep and to generate the
     ! representation basis
     !
     real(kind=r8_kind)              :: gen_eigen_value
     ! eigenvalue of the eigenvector, from which the basis of the
     ! representation is built
     logical                         :: build_projector
     ! build projector onto that the eigenvector of eigenvalue gen_eigen_value
     ! this must be done if the eigenvector was forced to be real
     real(kind=r8_kind), allocatable :: irrep_matrix(:,:,:)
     ! representation matrix: irrep_matrix(dimension,dimension,group_num_el)

     integer(kind=i4_kind), allocatable :: irrep_gen_oper(:)
     ! symmetry elements, which generate the basis of the irrep
     ! irrep_gen_oper(dimension)
     real(kind=r8_kind), allocatable :: irrep_gen_coeff(:,:)
     ! coefficients, which linearily combine the symmetry elements
     ! selected in irrep_gen_oper in order to directly give the
     ! basis functions:
     ! irrep_gen_coeff(operator,basisvector)
  end type group_irrep

  type, public :: group_proj_irrep
     ! projective irreps of the group

     !
     ! variables which identify the irrep
     !
     character(len=8)                :: label
     ! irrep label
     real(kind=r8_kind)              :: eigen_value
     ! eigenvalue of CSCO of irrep
     integer(kind=i4_kind)           :: dimension
     ! dimension of the irrep
     integer(kind=i4_kind)           :: time_dimension
     ! time inversion invariant dimension of the irrep
     logical                         :: pseudo
     ! .true. if the irrep is a pseudo-2D irrep
     real(kind=r8_kind), allocatable :: characters(:)
     ! characters of the irrep: characters(group_num_cl)
     real(kind=r8_kind)              :: sign_of_jz
     ! for axial groups: indicates if jz is positive or negative

     !
     ! the following variables are used to store the
     ! representation matrices of the irrep and to generate the
     ! representation basis
     !
     real(kind=r8_kind)              :: gen_eigen_value
     ! eigenvalue of the eigenvector, from which the basis of the
     ! representation is built
     logical                         :: build_projector
     ! build projector onto that the eigenvector of eigenvalue gen_eigen_value
     ! this must be done if the eigenvector was forced to be real
     complex(kind=c16_kind), allocatable :: irrep_matrix(:,:,:)
     ! representation matrix: irrep_matrix(dimension,dimension,group_num_el)

     integer(kind=i4_kind), allocatable :: irrep_gen_oper(:)
     ! symmetry elements, which generate the basis of the irrep
     ! irrep_gen_oper(dimension)
     complex(kind=c16_kind), allocatable :: irrep_gen_coeff(:,:)
     ! coefficients, which linearily combine the symmetry elements
     ! selected in irrep_gen_oper in order to directly give the
     ! basis functions:
     ! irrep_gen_coeff(operator,basisvector)
  end type group_proj_irrep

  type, public ::  sub_group
     ! subgroups of the group
     integer(kind=i4_kind)           :: id
     ! id number of subgroup
     character(len=4)                :: name
     ! name of the subgroup
     integer(kind=i4_kind)           :: num_el
     ! number of elements in the subgroup
     integer(kind=i4_kind), allocatable :: elements(:)
     ! ids of those elements of the whole group, which
     ! are elements of the subgroup as well (elements(num_el))
     integer(kind=i4_kind)           :: num_cl
     ! number of classes of the subgroup
     integer(kind=i4_kind)           :: num_re
     ! number of regular classes of the subgroup
     integer(kind=i4_kind)           :: num_ir
     ! number of (pseudo) irreps of the subgroup
     type(group_class), allocatable  :: klass(:)
     ! classes of the subgroup
     type(group_class), allocatable  :: pseudo_class(:)
     ! pseudo classes of the subgroup (The class and the class of the inverse are put
     ! into one pseudo class)
     integer(kind=i4_kind), allocatable :: super_class(:)
     ! ids of classes in the super group
     complex(kind=c16_kind), allocatable :: characters(:,:)
     ! character table for pseudo irreps characters(num_ir,num_ir)
     logical                         :: invariant
     ! flag for invariant subgroups
     integer(kind=i4_kind)           :: ncsco(3)
     ! classes used in csco of subgroup
     integer(kind=i4_kind)           :: fcsco(3)
     ! factor used for class in csco of subgroup
     integer(kind=i4_kind)           :: ncsco_proj(3)
     ! classes used in csco of subgroup for rep groups
     integer(kind=i4_kind)           :: fcsco_proj(3)
     ! factor used for class in csco of subgroup for rep groups
     integer(kind=i4_kind)           :: inversion,horizontal
     ! is inversion or horizontal plane present
     integer(kind=i4_kind)           :: primary_axis,secondary_axis,main_axis
     ! symmetry axes
     integer(kind=i4_kind)           :: order_primary_axis,order_secondary_axis
     ! orders of axes
     type(group_irrep), allocatable  :: irrep_can(:)
     ! canonically ordered irreps of the subgroup
     integer(kind=i4_kind), allocatable :: super_irrep(:)
     ! index of irreps to irreps in the super group (only possible
     ! for invariant subgroups)
  end type sub_group


  type, public :: group_coset
     ! cosets of the group for some subgroup
     integer(kind=i4_kind)           :: num_cs
     ! number of cosets
     integer(kind=i4_kind)           :: num_el
     ! number of elements in one coset
     integer(kind=i4_kind), allocatable :: elements(:,:)
     ! elements(num_el,num_cs)
  end type group_coset

  ! the following type contains the transformation matrices for some
  ! angular momentum

  type, public ::  symm_transformation
     real(kind=r8_kind), allocatable :: matrix(:,:,:)
     ! matrix(dim_trafo,dim_trafo,group_num_el)
     !                                 ^ group operations
  end type symm_transformation

  type, public ::  symm_transformation_c
     complex(kind=c16_kind), allocatable :: matrix(:,:,:)
     ! matrix(dim_trafo,dim_trafo,group_num_el)
     !                                 ^ group operations
  end type symm_transformation_c

  type, public ::  symm_basis_transformation_c
     complex(kind=c16_kind), allocatable :: matrix(:,:)
     ! matrix(dim_trafo,dim_trafo)
  end type symm_basis_transformation_c

  type, public ::  symm_transformation_int
     integer(kind=i4_kind), allocatable :: matrix(:,:,:)
     ! matrix(dim_trafo,dim_trafo,group_num_el)
     !                                 ^ group operations
  end type symm_transformation_int


  !------------ Declaration of constants and variables ---------------

  integer (i4_kind), parameter, public :: nptgrp=75
  ! number of available pointgroups
  integer (i4_kind), public, protected :: ncscop_proj(3, nptgrp)
  ! ncsco for projective groups
  integer (i4_kind), public, protected :: fcscop_proj(3, nptgrp)
  ! fcsco for projective groups
  character (len=4), public :: group_name ! unprotected
  ! name of symmetry group
  character (len=4), public, protected :: nonamb_group_name
  ! in case the symmetrygroup is nonambivalent this is the actual group name,
  ! while group_name is an ambivalent supergroup
  type (sub_group), public, protected :: nonamb_group
  ! in case the symmetrygroup is nonambivalent this is the actual group
  logical, public, protected :: nonambivalent
  ! point group is nonambivalent
  integer (i4_kind), public, protected :: igroup
  ! identification number of symmetry group
  integer (i4_kind), public, protected :: group_num_el
  ! number of group elements
  integer (i4_kind), public, protected :: group_num_cl
  ! number of classes
  integer (i4_kind), public, protected :: group_num_re
  ! number of regular classes
  integer (i4_kind), public, protected :: group_num_ir
  ! number of irreps
  integer (i4_kind), public, protected :: group_num_pir
  ! number of projective irreps
  integer (i4_kind), public, protected :: group_num_sg
  ! number of subgroups
  integer (i4_kind), public, protected :: ncsco(3)
  ! classes used in csco
  integer (i4_kind), public, protected :: fcsco(3)
  ! factor used for class in csco
  integer (i4_kind), public, protected :: ncsco_proj(3)
  ! classes used in csco of rep group
  integer (i4_kind), public, protected :: fcsco_proj(3)
  ! factor used for class in csco of rep group

  ! generator of group

  integer (i4_kind), public, protected :: inversion, horizontal
  ! inversion and horizontal reflection plane
  integer (i4_kind) :: primary_axis, secondary_axis, main_axis
  ! primary and secondary axis, the one of higher order is the main axis
  integer (i4_kind) :: order_primary_axis, order_secondary_axis
  ! orders of the two generator axis


  integer (i4_kind), allocatable, public, protected :: group_multab(:, :)
  ! multiplication table of symmetry elements
  ! R(i)*R(j) = group_multab(i, j)
  integer (i4_kind), allocatable, public, protected :: group_factab(:, :)
  ! table of the factor system [i, j] of the group
  integer (i4_kind), allocatable, public, protected :: class_mult(:, :, :)
  ! multiplication table of symmetry elements
  integer (i4_kind), allocatable, public, protected :: class_fact(:, :, :)
  ! factor system of class multiplication
  complex (c16_kind), allocatable, public, protected :: characters(:, :)
  ! character table: character(k, j) k=symm_el, j=irrep
  complex (c16_kind), allocatable, public, protected :: proj_characters(:, :)
  ! character table of proj irreps: proj_character(k, j) k=symm_el, j=irrep
  real (r8_kind), allocatable, public, protected :: eig_val(:)
  ! eigenvalues of the CSCO (irrep_labels)
  real (r8_kind), allocatable, public, protected :: proj_eig_val(:)
  ! eigenvalues of the CSCO of rep group (proj irrep_labels)
  type (group_irrep), allocatable, public :: irrep_can(:) ! unprotected
  ! irreps in canonical order: irrep_can(group_num_cl)
  type (group_proj_irrep), pointer, public :: proj_irrep_can(:)
  ! irreps in canonical order: irrep_can(group_num_cl), FIXME: make it
  ! allocatable, protected
  integer (i4_kind), allocatable, public, protected :: regular(:)
  ! vector subscript list of regular classes regular(group_num_re)
  integer (i4_kind), allocatable, public, protected :: regular_inv(:)
  ! inverted vector subscript list of regular classes regular(group_num_re)


  type (group_symm_el), allocatable, public, protected :: symm_el(:)
  ! symmetry elements of the group symm_el(group_num_el)
  type (group_class), allocatable, public, protected :: klass(:)
  ! classes of the group class(group_num_cl)
  type (sub_group), allocatable, public, protected :: subgroups(:)
  ! canonical subgroup chain of the group
  type (sub_group), allocatable, public, protected :: subgroups_proj(:)
  ! canonical subgroup chain of the group (for projective irreps)

  type (symm_transformation), allocatable, public :: ylm_trafos(:) ! unprotected
  ! transformation matrix of the spherical harmonics of the group
  ! efm_trafos(unique_atom_lmax_all)
  type (symm_transformation_c), allocatable, public :: ylm_proj_trafos(:, :) ! unprotected
  ! transformation matrix of the spherical harmonics of the group
  ! efm_trafos(unique_atom_lmax_all,2)
  type (symm_basis_transformation_c), allocatable, public :: comp_to_real(:, :, :) ! unprotected
  ! transformation matrix of the spherical harmonics of the group
  ! efm_trafos(unique_atom_lmax_all, 2)
  type (symm_transformation), public, protected :: E_irrep
  ! This is a real representation of the E Irrep of the group T, which
  ! was introduced in order to obtain real representation matrices

  !------------ Interface statements ---------------------------------

  interface operator(*)
     module procedure group_symmel_mult
  end interface

  interface operator(==)
     module procedure group_symmel_id
  end interface

  interface contains_inversion
     module procedure contains_inversion_el
  end interface

  !------------ public functions and subroutines ---------------------
  public &
        group_symmel_alloc, &
        group_symmel_read, &
        group_symmel_dealloc,&
        group_symmel_multab,&
        group_class_gen,&
        group_char_gen,&
        group_generator_gen,&
        group_irrep_label,&
        subgroup_chain_gen, &
        group_coset_decomp, &
        group_proj_char_gen, &
        group_close,  &
        group_proj_irrep_pre_label,&
        contains_inversion


  !===================================================================
  ! End of public interface of module
  !===================================================================



  ! name of pointgroup:
  character (len=4), public, protected :: ptgrps(nptgrp) = &
       &(/'C1  ','C2  ','C3  ','C4  ','C5  ','C6  ','C7  ','C8  ','C9  ',&
       &  'C10 ','Ci  ','CS  ','S4  ','S6  ','S8  ','S10 ','S12 ','S14 ',&
       &  'S16 ','S18 ','S20 ','D2  ','D3  ','D4  ','D5  ','D6  ','D7  ',&
       &  'D8  ','D9  ','D10 ','D2H ','D3H ','D4H ','D5H ','D6H ','D7H ',&
       &  'D8H ','D9H ','D10H','Dinh','D2D ','D3D ','D4D ','D5D ','D6D ',&
       &  'D7D ','D8D ','D9D ','D10D','C2V ','C3V ','C4V ','C5V ','C6V ',&
       &  'C7V ','C8V ','C9V ','C10V','Cinv','C2H ','C3H ','C4H ','C5H ',&
       &  'C6H ','C7H ','C8H ','C9H ','C10H','O   ','T   ','OH  ','TH  ',&
       &  'TD  ','I   ','IH  '/)

  ! read order of point groups
  ! order of pointgroup:
  integer (i4_kind), public, protected :: korder(nptgrp) = &
       &(/   1,     2,     3,     4,     5,     6,     7,     8,     9,&
       &    10,     2,     2,     4,     6,     8,    10,    12,    14,&
       &    16,    18,    20,     4,     6,     8,    10,    12,    14,&
       &    16,    18,    20,     8,    12,    16,    20,    24,    28,&
       &    32,    36,    40,    16,     8,    12,    16,    20,    24,&
       &    28,    32,    36,    40,     4,     6,     8,    10,    12,&
       &    14,    16,    18,    20,     8,     4,     6,     8,    10,&
       &    12,    14,    16,    18,    20,    24,    12,    48,    24,&
       &    24,    60,   120 /)

  ! read CSCO`s of projective representations
  integer (i4_kind), private :: k_, i_ ! counters
  data ((ncscop_proj(i_, k_), i_ = 1, 3), k_ = 1, 18) &
       & /  1, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0,&
       &    2, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0,  2, 0, 0, &
       &    4, 0, 0,  6, 0, 0,  8, 0, 0, 10, 0, 0, 12, 0, 0, 14, 0, 0 /
  data ((fcscop_proj(i_, k_), i_ = 1, 3), k_ = 1, 18) &
       & /  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,&
       &    1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0, &
       &    1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0 /

  data ((ncscop_proj(i_, k_), i_ = 1, 3), k_ = 19, 36) &
       & / 16, 0, 0, 18, 0, 0, 20, 0, 0,  1, 0, 0,  3, 0, 0,  2, 0, 0, &
       &    2, 4, 0,  2, 0, 0,  2, 5, 0,  2, 0, 0,  2, 6, 0,  2, 0, 0, &
       &    5, 0, 0,  5, 0, 0,  2, 6, 0,  2, 6, 0,  2, 7, 0,  2, 7, 0 /
  data ((fcscop_proj(i_, k_), i_ = 1, 3), k_ = 19, 36) &
       & /  1, 0, 0,  1, 0, 0,  1, 0, 0,  2, 0, 0,  1, 0, 0,  2, 0, 0,&
       &    2, 1, 0,  2, 0, 0,  2, 1, 0,  2, 0, 0,  2, 1, 0,  2, 0, 0, &
       &    1, 0, 0,  1, 0, 0,  4, 1, 0,  4, 1, 0,  4, 1, 0,  4, 1, 0 /

  data ((ncscop_proj(i_, k_), i_ = 1, 3), k_ = 37, 54) &
       & /  2, 8, 0,  2, 7, 0,  2, 9, 0,  2, 4, 6,  4, 0, 0,  4, 6, 0, &
       &    6, 0, 0,  7, 8, 0,  8, 0, 0,  9,10, 0, 10, 0, 0, 11,12, 0, &
       &   12, 0, 0,  1, 0, 0,  3, 0, 0,  2, 0, 0,  2, 4, 0,  2, 0, 0 /
  data ((fcscop_proj(i_, k_), i_ = 1, 3), k_ = 37, 54) &
       & /  2, 2, 0,  2, 1, 0,  4, 1, 0,  2, 1, 1,  2, 0, 0,  2, 1, 0,&
       &    3, 0, 0,  2, 1, 0,  4, 0, 0,  2, 1, 0,  2, 0, 0,  2, 1, 0, &
       &    2, 0, 0,  2, 0, 0,  1, 0, 0,  2, 0, 0,  2, 1, 0,  2, 0, 0 /

  data ((ncscop_proj(i_, k_), i_ = 1, 3), k_ = 55, 72) &
       & /  2, 5, 0,  2, 0, 0,  2, 6, 0,  2, 0, 0,  2, 4, 0,  2, 3, 0, &
       &    2, 4, 0,  2, 5, 0,  2,10, 0,  2, 7, 0,  2, 8, 0,  2, 9, 0, &
       &    2,10, 0,  2,11, 0,  4, 0, 0,  3, 0, 0,  4, 6, 0,  3, 5, 0 /
  data ((fcscop_proj(i_, k_), i_ = 1, 3), k_ = 55, 72) &
       & /  2, 1, 0,  2, 0, 0,  2, 1, 0,  2, 0, 0,  2, 1, 0,  2, 1, 0,&
       &    1, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0, &
       &    2, 5, 0,  2, 5, 0,  1, 0, 0,  1, 0, 0,  2, 1, 0,  1, 5, 0 /

  data ((ncscop_proj(i_, k_), i_ = 1, 3), k_ = 73, 75) &
       & /  4, 0, 0,  2, 0, 0,  2, 6, 0 /
  data ((fcscop_proj(i_, k_), i_ = 1, 3), k_ = 73, 75) &
       & /  1, 0, 0,  2, 0, 0,  2, 1, 0 /


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  function contains_inversion_el(el) result(yes)
    implicit none
    integer(i4_kind),intent(in) :: el
    logical                     :: yes !<<<result
    ! *** end of interface ***

    if( el.le.0.or.el.gt.group_num_el )&
         & call error_handler("grpm/contains_inversion_el: el out of range")

    yes = ( NINT(symm_el(el)%quaternion(5)).eq.-1 )
    ! i.e. contains inversion
  end function contains_inversion_el

  !*************************************************************
  subroutine group_symmel_alloc
    !  Purpose: allocates space for symmetry elements
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: k,error
    ! counter, error_code
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    ! identify redundant names

    if (trim(adjustl(group_name)).eq.'Dinh') then
       group_name = 'D4h'
    end if
    if (trim(adjustl(group_name)).eq.'Cinv') then
       group_name = 'C4v'
    end if

    ! determine group number igroup
    igroup = 0
    do k=1,nptgrp
       DPRINT 'gm::GROUP_SYMMEL_ALLOC: cmp',k,'>'//group_name//'<>'//ptgrps(k)//'<'
!!$       if (trim(adjustl(group_name)).eq.trim(adjustl(ptgrps(k)))) then
       if(group_name.eq.ptgrps(k))then
          igroup = k
          exit
       end if
    end do
    if (igroup.eq.0) then
       call error_handler("GROUP_SYMMEL_ALLOC: wrong point group")
    endif

    ! determine group size
    group_num_el = korder(k)

    ! allocate symmetry_elements

    allocate (symm_el(group_num_el),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of symm_el failed")

    ! allocate multiplication table

    allocate (group_multab(group_num_el,group_num_el),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of multab failed")

    ! allocate factor system

    allocate (group_factab(group_num_el,group_num_el),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of factab failed")


  end subroutine group_symmel_alloc
  !************************************************************


  !*************************************************************
  subroutine group_number_size(group_name,number,group_size)
    !  Purpose: determines the group_id number and the group size
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none
    character(len=4), intent(in)       :: group_name
    ! name of group to be investigated
    integer(kind=i4_kind),intent(out)  :: number
    ! group id (to be determined)
    integer(kind=i4_kind),intent(out)  :: group_size
    ! size of the group (to be determined)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: k
    ! counter
    !------------ Executable code ------------------------------------

    ! determine group number igroup
    number = 0
    do k=1,nptgrp
       if (group_name.eq.ptgrps(k)) then
          number = k
          exit
       end if
    end do

    if (number.eq.0) call error_handler( &
         "group_number_size: number.eq.0" )

    ! determine group size
    group_size = korder(number)



  end subroutine group_number_size
  !************************************************************


  !*************************************************************
  subroutine group_symmel_read
    !  Purpose: reads symmetry data
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable,dimension(:,:) :: elem
    ! work array for quaternions of elements of point group
    integer(kind=i4_kind)                :: k,l
    ! counter, error_code
    character(len=8), allocatable, dimension(:) :: namele
    ! work array for name of elements of point group
    !------------ Executable code ------------------------------------

    ! allocate work arrays
    allocate( elem(5,group_num_el),namele(group_num_el))

    ! generate group elements
    call gengrp(group_name,igroup,group_num_el,elem,namele,&
         ncsco,fcsco)

    ! do not use preliminary identification of ambivalent classes
    ncsco = abs(ncsco)

    ! set CSCO of the representation group
    ncsco_proj = ncscop_proj(:,igroup)
    fcsco_proj = fcscop_proj(:,igroup)

    ! load group elements to new data structure
    do k=1,group_num_el
       symm_el(k)%quaternion = elem(:,k)
       symm_el(k)%name = namele(k)
    end do

    ! write to output file
    if (output_unit > 0) then
       write (output_unit, 605) group_name, &
            (l,symm_el(l)%name, (symm_el(l)%quaternion(k),k=1,5),l=1,group_num_el)
    endif
605 format(/1x,'Elements of group ',a4,', quaternionic parameters:'//&
         &      (1x,i4,1x,a8,4f13.8,f5.0))

    ! deallocate work arrays
    deallocate(elem,namele)

  end subroutine group_symmel_read
  !*************************************************************


  !*************************************************************
  subroutine subgroup_symmel_read(subgroup_name,subgroup)
    !  Purpose: read symmetry elements of the subgroup,
    !           identifies them with those of main group
    !           Then the classes of the subgroup are determined
    !           by subroutine subgroup_class_gen
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    character(len=4), intent(in)         :: subgroup_name
    ! name of subgroup
    type(sub_group)                      :: subgroup
    ! points to newly created subgroup
    !** End of interface *****************************************
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable,dimension(:,:) :: elem
    ! work array for quaternions of elements of point group
    integer(kind=i4_kind)                :: k,l
    ! counter, error_code
    integer(kind=i4_kind)                :: group_size
    ! size of the subgroup
    type(group_symm_el)                  :: symmetry_el
    ! some symmetryelement
    character(len=8), allocatable, dimension(:) :: namele
    ! work array for name of elements of point group
    !------------ Executable code ------------------------------------

    ! determine id and size of the subgroup
    call group_number_size(subgroup_name,subgroup%id,subgroup%num_el)

    group_size = subgroup%num_el

    ! allocate list of elements of subgroup
    allocate(subgroup%elements(group_size))

    ! allocate work arrays
    allocate(elem(5,group_size),namele(group_size))

    ! generate group elements
    call gengrp(subgroup_name,subgroup%id,group_size,elem,namele,&
         subgroup%ncsco,subgroup%fcsco)

    ! give name to subgroup
    subgroup%name = subgroup_name

    ! do not use preliminary identification of ambivalent classes
    subgroup%ncsco = abs(subgroup%ncsco)

    ! set CSCO of the representation group
    subgroup%ncsco_proj = ncscop_proj(:,subgroup%id)
    subgroup%fcsco_proj = fcscop_proj(:,subgroup%id)

    ! check if subgroup is contained in the group (subgroup test)
    outest: do
       outer:     do k=1,group_size
          symmetry_el%quaternion = elem(:,k)
          do l=1,group_num_el
             if (symmetry_el==symm_el(l)) then
                ! symmetry element is contained
                subgroup%elements(k) = l
                cycle outer
             end if
          end do
          ! symmetry element was not found in group
          ! => subgroup test failed
          if (output_unit > 0) then
             write(output_unit,*) subgroup_name," is not a subgroup of ", group_name
          endif
          exit outest
       end do outer
       ! all symmetry elements of subgroup were found in the group
       ! => subgroup test passed
       if (output_unit > 0) then
          write(output_unit,*)
          write(output_unit,*)"----------------------------"
          write(output_unit,*) "Subgroup ",subgroup_name," of ", group_name
          write(output_unit,*)"----------------------------"
       endif
       ! now generate the classes of the subgroup
       call subgroup_class_gen(subgroup%num_el,subgroup%elements,subgroup,&
            &subgroup%num_cl)
       exit
    end do outest

    ! deallocate work arrays
    deallocate(elem,namele)

  end subroutine subgroup_symmel_read
  !************************************************************


  !************************************************************
  subroutine subgroup_chain_gen
    !  Purpose: generates a canonical subgroup chain
    !** End of interface *****************************************
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of local variables ---------------------
    !integer(kind=i4_kind)                :: k


    !------------ Executable code ------------------------------------

    ! select by group
    select case(group_name)

       ! polyhedral groups
    case('T   ','TH  ')
       group_num_sg = 1
       allocate(subgroups(group_num_sg))
       allocate(subgroups_proj(group_num_sg))
       call subgroup_symmel_read('D2  ',subgroups(1))
       call subgroup_symmel_read('D2  ',subgroups_proj(1))
    case('TD  ','O   ','OH  ','I   ','IH  ')
       group_num_sg = 2
       allocate(subgroups(group_num_sg))
       allocate(subgroups_proj(group_num_sg))
       call subgroup_symmel_read('T   ',subgroups(1))
       call subgroup_symmel_read('T   ',subgroups_proj(1))
       call subgroup_symmel_read('D2  ',subgroups(2))
       call subgroup_symmel_read('C2  ',subgroups_proj(2))
       !case('Td  ')
       !   group_num_sg = 1
       !   allocate(subgroups(group_num_sg))
       !   call subgroup_symmel_read('S4  ',subgroups(1))

       ! axial groups
    case default
       if (secondary_axis.ne.0) then
          ! build subgroup from the two elements secondary axis(C2-rotation
          ! or vertical reflection plane) and identity

          group_num_sg = 1
          allocate(subgroups(group_num_sg))
          allocate(subgroups_proj(group_num_sg))
          allocate(subgroups(1)%elements(2))

          ! some name for the subgroup
          subgroups(1)%name = 'SEC '
          ! the subgroup has two elements
          subgroups(1)%num_el = 2
          ! first element is the identity
          subgroups(1)%elements(1) = 1
          ! second element is the secondary axis
          subgroups(1)%elements(2) = klass(secondary_axis)%conjug_el(1)
          ! second element is CSCO of subgroup

          subgroups(1)%ncsco(1) = 2
          subgroups(1)%fcsco(1) = 1

          subgroups(1)%ncsco(2:3) = 0
          subgroups(1)%fcsco(2:3) = 0

          call subgroup_class_gen(subgroups(1)%num_el,subgroups(1)%elements,subgroups(1),&
               &subgroups(1)%num_cl)

          if (group_name.eq.'D4H ') then
             ! this will also allocate subgroups_proj(1)%elements(:) ...
             call subgroup_symmel_read('S4  ', subgroups_proj(1))
          else
             ! otherwise allocate the space here:
             allocate(subgroups_proj(1)%elements(2))

             ! some name for the subgroup
             subgroups_proj(1)%name = 'SEC '
             ! the subgroup has two elements
             subgroups_proj(1)%num_el = 2
             ! first element is the identity
             subgroups_proj(1)%elements(1) = 1
             ! second element is the secondary axis
             subgroups_proj(1)%elements(2) = klass(secondary_axis)%conjug_el(1)
             ! second element is CSCO of subgroup


             ! second element is CSCO of subgroup
             subgroups_proj(1)%ncsco_proj(1) = 2
             subgroups_proj(1)%fcsco_proj(1) = 1

             subgroups_proj(1)%ncsco_proj(2:3) = 0
             subgroups_proj(1)%fcsco_proj(2:3) = 0

             ! generate classes of the subgroup
             call subgroup_class_gen(subgroups_proj(1)%num_el,subgroups_proj(1)%elements,subgroups_proj(1),&
                  &subgroups_proj(1)%num_cl)
          endif
       else
          ! to avoid special cases, allocate anyway:
          group_num_sg = 0
          allocate(subgroups(group_num_sg))
          allocate(subgroups_proj(group_num_sg))
       end if

    end select

  end subroutine subgroup_chain_gen
  !************************************************************


  !************************************************************
  subroutine group_symmel_dealloc
    !  Purpose: dealocates space of symmetry elements
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: error


    !------------ Executable code ------------------------------------

    deallocate (symm_el,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of symm_el failed")

    deallocate (group_multab,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of group_multab failed")

    deallocate (group_factab,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of group_factab failed")

    deallocate (class_mult,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of class_mult failed")

    deallocate (class_fact,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of class_fact  failed")

    deallocate (characters,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of characters failed")

    deallocate (eig_val,stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_dealloc: deallocation of eig_val failed")


  end subroutine group_symmel_dealloc
  !************************************************************


  !************************************************************
  subroutine group_close()
    !
    !  Purpose: clean up module. This is called when finalizing a
    !           geometry calculation. That is in a geometry
    !           optimization it will be called many times. But,
    !           historically, for each geometry all phases are
    !           repeated so that the group data will be re-generated
    !           again. This sub is supposed to do its best to clean
    !           the state of the module. Irrespective if the next
    !           calculation is using the same or a different symmetry.
    !
    use type_module
    use iounitadmin_module
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: error, i, j, n

    ! deallocate 'proj_irrep_can' --------------------------
    if ( associated(proj_irrep_can) ) then
       !
       ! CASE OF SPIN-ORBIT
       !
       n = ubound(proj_irrep_can,1)
       do i=1,n
          deallocate(proj_irrep_can(i)%characters,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (1) failed")
          deallocate(proj_irrep_can(i)%irrep_matrix,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (2) failed")
          deallocate(proj_irrep_can(i)%irrep_gen_oper,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (3) failed")
          deallocate(proj_irrep_can(i)%irrep_gen_coeff,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (4) failed")
       enddo
       ! FIXME: with allocatable components this should suffice:
       deallocate(proj_irrep_can,stat=error)
       if (error/=0) call error_handler &
            ("group_close : deallocation (5) failed")

       deallocate(regular, stat=error)
       ASSERT(error==0)

       deallocate(regular_inv, stat=error)
       ASSERT(error==0)

       deallocate(proj_characters, stat=error)
       ASSERT(error==0)

       deallocate(proj_eig_val, stat=error)
       ASSERT(error==0)
    endif
    ! -----------------------------------------------------

    ! FIXME: what about irrep_can?

    ! deallocate 'klass' ----------------------------------
    n = ubound(klass,1)
    do i=1,n
       deallocate(klass(i)%conjug_el,stat=error)
       if (error/=0) call error_handler &
            ("group_close : deallocation (12) failed")
    enddo
    deallocate(klass,stat=error)
    if (error/=0) call error_handler &
         ("group_close : deallocation (13) failed")
    ! -----------------------------------------------------

    if( allocated(subgroups) ) then
       ! deallocate 'subgroups' ------------------------------
       do i = 1, size(subgroups)
          ! class and pseudo_class (type group_class)
          do j = 1, size(subgroups(i)%klass)
             deallocate(subgroups(i)%klass(j)%conjug_el,stat=error)
             if (error/=0) call error_handler &
                  ("group_close : deallocation of subgroups (14) failed")
          enddo
          do j = 1, size(subgroups(i)%pseudo_class)
             if (allocated(subgroups(i)%pseudo_class(j)%conjug_el)) then
                deallocate(subgroups(i)%pseudo_class(j)%conjug_el,stat=error)
                if (error/=0) call error_handler &
                     ("group_close : deallocation (15) failed")
             endif
          enddo
          deallocate(subgroups(i)%super_class,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (16) failed")
       enddo
       deallocate(subgroups,stat=error)
       if (error/=0) call error_handler &
            ("group_close : deallocation (23) failed")
    endif

    if( allocated(subgroups_proj) ) then
       ! deallocate 'subgroups_proj' ------------------------------
       do i = 1, size(subgroups_proj)
          ! class and pseudo_class (type group_class)
          do j = 1, size(subgroups_proj(i)%klass)
             deallocate(subgroups_proj(i)%klass(j)%conjug_el,stat=error)
             if (error/=0) call error_handler &
                  ("group_close : deallocation of subgroup_proj (14) failed")
          enddo
          do j = 1, size(subgroups_proj(i)%pseudo_class)
             if (allocated(subgroups_proj(i)%pseudo_class(j)%conjug_el)) then
                deallocate(subgroups_proj(i)%pseudo_class(j)%conjug_el,stat=error)
                if (error/=0) call error_handler &
                     ("group_close : deallocation (15) failed")
             endif
          enddo
          deallocate(subgroups_proj(i)%super_class,stat=error)
          if (error/=0) call error_handler &
               ("group_close : deallocation (16) failed")
       enddo
       deallocate(subgroups_proj,stat=error)
       if (error/=0) call error_handler &
            ("group_close : deallocation (23) failed")
    end if
    ! ---------------------------------------------------------
  end subroutine group_close
  !************************************************************


  !************************************************************
  type(group_symm_el) function group_symmel_mult(s1,s2)
    !  Purpose: multiplies two symmetry elements
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(group_symm_el),intent(in)   :: s1,s2
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)           :: x1(5),x2(5),y(5)

    !------------ Executable code ------------------------------------

    x1=s1%quaternion
    x2=s2%quaternion

    y(2) = x1(1)*x2(2) + x1(2)*x2(1) + x1(3)*x2(4) - x1(4)*x2(3)
    y(3) = x1(1)*x2(3) + x1(3)*x2(1) + x1(4)*x2(2) - x1(2)*x2(4)
    y(4) = x1(1)*x2(4) + x1(4)*x2(1) + x1(2)*x2(3) - x1(3)*x2(2)
    y(1) = x1(1)*x2(1) - x1(2)*x2(2) - x1(3)*x2(3) - x1(4)*x2(4)

    y(5) = x1(5)*x2(5)

    group_symmel_mult%quaternion = y

  end function group_symmel_mult
  !************************************************************


  !************************************************************
  logical function group_symmel_id(s1,s2)
    !
    !  Purpose: verifies if two symmetry elements are identical
    !
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none

    !------------ Declaration of formal parameters ---------------
    type(group_symm_el),intent(in)   :: s1,s2

    !** End of interface *****************************************
    !------------ Declaration of local constants  ----------------
    real(kind=r8_kind),parameter :: small = 1.e-10

    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)           :: x1(5),x2(5)
    real(kind=r8_kind)           :: diff1, diff2

    !------------ Executable code ------------------------------------

    x1=s1%quaternion
    x2=s2%quaternion


    diff1 = ABS(x1(1)-x2(1)) + ABS(x1(2)-x2(2))+ &
         ABS(x1(3)-x2(3)) + ABS(x1(4)-x2(4))+    &
         ABS(x1(5)-x2(5))

    diff2 = ABS(x1(1)+x2(1)) + ABS(x1(2)+x2(2))+ &
         ABS(x1(3)+x2(3)) + ABS(x1(4)+x2(4))+    &
         ABS(x1(5)-x2(5))


    if ((diff1.lt.small).or.(diff2.lt.small)) then
       group_symmel_id = .true.
    else
       group_symmel_id = .false.
    end if

  end function group_symmel_id
  !************************************************************


  !************************************************************
  integer function group_symmel_factor(s1, s2)
    !  Purpose: determines factor between two symmetry identical
    !           symmetry elements
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(group_symm_el), intent(in) :: s1, s2
    !** End of interface *****************************************

    external error_handler
    !------------ Declaration of local constants  ----------------
    real(kind=r8_kind),parameter :: small = 1.e-10

    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)           :: x1(5),x2(5)
    real(kind=r8_kind)           :: diff1, diff2

    !------------ Executable code ------------------------------------

    x1=s1%quaternion
    x2=s2%quaternion


    diff1 = ABS(x1(1)-x2(1)) + ABS(x1(2)-x2(2))+ &
         ABS(x1(3)-x2(3)) + ABS(x1(4)-x2(4))+    &
         ABS(x1(5)-x2(5))

    diff2 = ABS(x1(1)+x2(1)) + ABS(x1(2)+x2(2))+ &
         ABS(x1(3)+x2(3)) + ABS(x1(4)+x2(4))+    &
         ABS(x1(5)-x2(5))


    if ((diff1.lt.small).or.(diff2.lt.small)) then
       if (diff1.lt.small) then
          group_symmel_factor = 1
       else
          group_symmel_factor = -1
       end if
    else
       group_symmel_factor = 0 ! to make compiler happy
       call error_handler( &
            "group_symmel_factor: symmetry elements differ")
    end if
  end function group_symmel_factor
  !************************************************************


  !*************************************************************
  subroutine group_symmel_multab
    !  Purpose: generates group multiplication table
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: i,j,k,jmax
    ! counters
    type(group_symm_el)                  :: y
    ! temporary results of group multiplication

    !------------ Executable code ------------------------------------

    do i=1,group_num_el
       do j=1,group_num_el
          ! multiply two group elements
          y = symm_el(i) * symm_el(j)
          ! check for result of multiplication
          do k=1,group_num_el
             if (y==symm_el(k)) then
                group_multab(i,j) = k
                ! determine [i,j]
                group_factab(i,j) = group_symmel_factor(y,symm_el(k))
                if (k==1) then
                   symm_el(i)%inverse = j
                   ! determine sign of inverse in spinor space
                   if (group_factab(i,j).lt.0) then
                      symm_el(i)%inverse_sign = -1
                   else
                      symm_el(i)%inverse_sign = +1
                   endif
                end if
                exit
             end if
          end do
       end do
    end do

    ! show group multiplication table
    if (output_unit > 0) then
       write(output_unit,fmt=701) group_name
       do j = 1, group_num_el, 15
          jmax = min(j+14,group_num_el)
          write(output_unit,fmt=702) (k, k=j,jmax)
          write(output_unit,fmt=703) ('----', k=j,jmax)
          do i = 1,group_num_el
             write(output_unit,fmt=704) i, symm_el(i)%name, (group_multab(i,k),&
                  k=j,jmax)
          end do
       end do

       ! show factor system
       write(output_unit,fmt=705) group_name
       do j = 1, group_num_el, 15
          jmax = min(j+14,group_num_el)
          write(output_unit,fmt=702) (k, k=j,jmax)
          write(output_unit,fmt=703) ('----', k=j,jmax)
          do i = 1,group_num_el
             write(output_unit,fmt=704) i, symm_el(i)%name, (group_factab(i,k),&
                  k=j,jmax)
          end do
       end do
    endif
    !
701 format(/1x,'Group table of group ',a)
702 format(/15x,15i4)
703 format(15x,15a)
704 format(1x,i4,1x,a8,'|',15i4)
705 format(/1x,'factor system ',a)

  end subroutine group_symmel_multab
  !*************************************************************


  !*************************************************************
  subroutine group_class_gen
    !  Purpose: generates classes
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Definition of local types ----------------------

    type class_list
       type(group_class)            :: klass
       type(class_list), pointer    :: next => NULL()
    end type class_list

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: remain(group_num_el)
    ! contains not referenced group elements
    integer(kind=i4_kind)                :: counter(group_num_el)
    ! help array which determines, if element is in class
    integer(kind=i4_kind)                :: counter_sign(group_num_el)
    ! help array which determines, if element is in class
    character(len=3), allocatable        :: regular_char(:)
    ! labels for the output
    integer(kind=i4_kind)                :: i,j,k,l,m,iclass
    ! counters
    integer(kind=i4_kind)                :: conjug_sign
    ! sign of conjugated element in spinor space
    logical                              :: class_regular
    ! .true. if the class is regular
    type(class_list),pointer             :: top,bottom,previous
    ! pointer to head of class list

    !initialize remain
    remain  = 1
    group_num_cl = 0
    outer: do
       !initialize count
       counter   = 0
       class_regular = .true.
       ! search for remaining elements
       do i=1,group_num_el
             ! if remaining element found generate class
             if (remain(i) == 0) cycle
             group_num_cl = group_num_cl + 1
             ! generate conjugated elements
             do j=1,group_num_el
                k = group_multab(j,group_multab(i,symm_el(j)%inverse))
                ! determine sign of conjugated element in spinor space
                conjug_sign = symm_el(j)%inverse_sign*group_factab(i,symm_el(j)%inverse)*&
                     &group_factab(j,group_multab(i,symm_el(j)%inverse))
                ! if the sign is negative the class is irregular
                if (conjug_sign.lt.0) then
                   class_regular = .false.
                endif
                counter(k)=1
             end do
             ! determine number of elements in class
             iclass = 0
             do j=1,group_num_el
                iclass = iclass + counter(j)
             end do
             ! allocate new class element
             allocate(top)
             allocate(top%klass%conjug_el(iclass))
             ! populate new class
             k=1
             do j=1, group_num_el
                if (counter(j)== 1) then
                   top%klass%conjug_el(k) = j
                   remain(j) = 0
                   ! attach class label to symmetry element
                   symm_el(j)%klass = group_num_cl
                   k = k + 1
                end if
             end do
             top%klass%number  = iclass
             top%klass%regular = class_regular
             if (group_num_cl>1) then
                previous%next => top
             else
                bottom => top
             end if
             previous => top
             cycle outer
       end do
       ! no more remaining elements found
       exit
    end do outer


    ! allocate global class list
    allocate(klass(group_num_cl))
    allocate(regular_char(group_num_cl))

    ! read generated classes to global list and deallocate linked list:
    top => bottom
    do j = 1, group_num_cl ! while ( associated(top) )
       klass(j) =  top%klass

       ! to be deallocated:
       previous => top

       ! advance head:
       top => top%next

       deallocate(previous) ! and clean the previous entry
    end do

    ! determine inverse class and number of regular classes
    group_num_re = 0
    do j=1,group_num_cl
       klass(j)%inverse = symm_el(symm_el(klass(j)%conjug_el(1))%inverse)&
            &%klass
       if (klass(j)%inverse == j) then
          klass(j)%ambivalent = .true.
          ! check ambivalence of rep group
          if (symm_el(klass(j)%conjug_el(1))%inverse_sign.gt.0) then
             klass(j)%ambivalent_proj = .true.
          else
             klass(j)%ambivalent_proj = .false.
          endif
       else
          klass(j)%ambivalent = .false.
          klass(j)%ambivalent_proj = .false.
       end if

       ! determine output label of regularity
       if (klass(j)%regular) then
          group_num_re = group_num_re + 1
          regular_char(j) = 'reg'
       else
          regular_char(j) = 'irr'
       end if
    end do

    ! allocate class multiplication table
    allocate(class_mult(group_num_cl,group_num_cl,group_num_cl))
    allocate(class_fact(group_num_cl,group_num_cl,group_num_cl))

    ! generate class multiplication and factor table

    ! loop over classes
    do i=1,group_num_cl
       do j=1,group_num_cl
          counter = 0
          counter_sign = 0
          ! loop over class elements
          do k=1,klass(i)%number
             do l=1,klass(j)%number
                m = group_multab(klass(i)%conjug_el(k),klass(j)%conjug_el(l))
                conjug_sign = group_factab(klass(i)%conjug_el(k),klass(j)%conjug_el(l))
                counter_sign(m) = counter_sign(m) + conjug_sign
                counter(m) = counter(m) + 1
             end do
          end do
          do k=1,group_num_cl
             class_mult(k,i,j) = counter(klass(k)%conjug_el(1))
             class_fact(k,i,j) = sign(1,counter_sign(klass(k)%conjug_el(1)))
          end do
       end do
    end do

    ! print class list
    if (output_unit > 0) then
       write(output_unit,602) group_name, group_num_cl,group_num_re
       write(output_unit,603)
       do i = 1, group_num_cl
          write(output_unit,604) i,klass(i)%inverse,regular_char(i),&
               &       (symm_el(klass(i)%conjug_el(j))%name, j = 1,klass(i)%number)
       end do
    endif
602 format(/1x,'Group ',a,' has ',i2,' classes (',i2,' are regular).')
603 format(/1x,' their elements are:'   /)
604 format(1x,i3,' (',i2,') ',a,' : ',8a8:,3(/12x,8a8:))
!05 format(1x,i3,' : ',5i3)

  end subroutine group_class_gen
  !*************************************************************


  !*************************************************************
  subroutine subgroup_class_gen(num_el,elements,subgroup,num_cl)
    !  Purpose: generates classes of the subgroup
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    integer(kind=i4_kind), intent(in)       :: num_el
    ! number of elements of the subgroup
    integer(kind=i4_kind), intent(in)       :: elements(:)
    ! elements of the subgroup
    integer(kind=i4_kind), intent(out)      :: num_cl
    ! number of classes in the subgroup
    type(sub_group)                         :: subgroup
    ! contains all information of the subgroup
    !** End of interface *****************************************
    !------------ Definition of local types ----------------------

    type class_list
       type(group_class)            :: klass
       type(class_list), pointer    :: next => NULL()
    end type class_list

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: remain(group_num_el)
    ! contains not referenced group elements
    integer(kind=i4_kind)                :: counter(group_num_el)
    ! help array which determines, if element is in class
    character(len=3), allocatable        :: regular_char(:)
    ! labels for the output
    integer(kind=i4_kind)                :: elem_class(group_num_el)
    ! contains class of subgroup, in which element is contained
    integer(kind=i4_kind)                :: i,i_elem,j,j_elem,k,iclass
    ! counters
    integer(kind=i4_kind)                :: class_size
    ! number of elements in class
    integer(kind=i4_kind)                :: conjug_sign
    ! sign of conjugated element in spinor space
    integer(kind=i4_kind)                :: inv_sign
    ! sign of inverse element in spinor space
    logical                              :: class_regular
    ! .true. if class is regular
    integer(kind=i4_kind)                :: num_ir
    ! number of pseudo irreps
    type(class_list),pointer             :: top,bottom,previous
    ! pointer to head of class list

    !initialize remain

    elem_class = 0
    remain  = 1
    num_cl = 0
    outer: do
       !initialize count
       counter   = 0
       class_regular = .true.
       ! search for remaining elements
       do i=1,num_el
             i_elem = elements(i)
             ! if remaining element found generate class
             if (remain(i_elem) == 0) cycle
             num_cl = num_cl + 1
             ! generate conjugated elements
             do j=1,num_el
                j_elem = elements(j)
                k = group_multab(j_elem,group_multab(i_elem,symm_el(j_elem)%inverse))
                ! determine sign of conjugated element in spinor space
                conjug_sign = symm_el(j_elem)%inverse_sign*group_factab(i_elem,symm_el(j_elem)%inverse)*&
                     &group_factab(j_elem,group_multab(i_elem,symm_el(j_elem)%inverse))
                ! if the sign is negative the class is irregular
                if (conjug_sign.lt.0) then
                   class_regular = .false.
                endif

                counter(k)=1
             end do
             ! determine number of elements in class
             iclass = 0
             do j=1,num_el
                j_elem = elements(j)
                iclass = iclass + counter(j_elem)
             end do
             ! allocate new class element
             allocate(top)
             allocate(top%klass%conjug_el(iclass))
             ! populate new class
             k=1
             do j=1,group_num_el
                if (counter(j)== 1) then
                   remain(j) = 0
                   elem_class(j) = num_cl
                   top%klass%conjug_el(k) = j
                   ! attach class label to symmetry element
                   !symm_el(j)%klass = num_cl
                   k = k + 1
                end if
             end do
             top%klass%number  = iclass
             top%klass%regular = class_regular
             if (num_cl>1) then
                previous%next => top
             else
                bottom => top
             end if
             previous => top
             cycle outer
       end do
       ! no more remaining elements found
       exit
    end do outer


    ! allocate subgroup class list
    allocate(subgroup%klass(num_cl))
    allocate(regular_char(num_cl))

    ! read generated classes to subgroup class list and deallocate linked list:
    top => bottom
    do j = 1, num_cl ! while ( associated(top) )
       subgroup%klass(j) = top%klass

       ! to be deallocated:
       previous => top

       ! advance head:
       top => top%next

       deallocate(previous) ! and clean the previous entry
    end do

    ! determine inverse class and number of regular classes
    ! set number of (pseudo) irreps = number of pseudo classes
    subgroup%num_ir = num_cl
    subgroup%num_re = 0
    do j=1,num_cl
       subgroup%klass(j)%inverse = elem_class(symm_el(subgroup%klass(j)%conjug_el(1))&
            &%inverse)
       if (subgroup%klass(j)%inverse == j) then
          subgroup%klass(j)%ambivalent = .true.
          ! check ambivalence of rep group
          if (symm_el(subgroup%klass(j)%conjug_el(1))%inverse_sign.gt.0) then
             subgroup%klass(j)%ambivalent_proj = .true.
          else
             subgroup%klass(j)%ambivalent_proj = .false.
          endif
       else
          subgroup%klass(j)%ambivalent = .false.
          subgroup%klass(j)%ambivalent_proj = .false.
          ! two nonambivalent class give rise to one pseudo class
          if (subgroup%klass(j)%inverse.gt.j) then
             subgroup%num_ir = subgroup%num_ir - 1
          end if
       end if

       ! determine output label of regularity
       if (subgroup%klass(j)%regular) then
          subgroup%num_re = subgroup%num_re + 1
          regular_char(j) = 'reg'
       else
          regular_char(j) = 'irr'
       end if

    end do

    !
    ! build pseudo classes
    !

    num_ir = subgroup%num_ir
    ! allocate pseudo classes
    allocate(subgroup%pseudo_class(num_ir))
    ! allocate identification list of pseudo classes with classes of super group
    allocate(subgroup%super_class(num_ir))

    i=1
    do j=1,num_cl
       if (subgroup%klass(j)%ambivalent) then
          subgroup%pseudo_class(i)%number = subgroup%klass(j)%number
          class_size = subgroup%pseudo_class(i)%number
          allocate(subgroup%pseudo_class(i)%conjug_el(class_size))
          subgroup%pseudo_class(i)%conjug_el = subgroup%klass(j)%conjug_el
          i = i + 1
       else
          if (subgroup%klass(j)%inverse.lt.j) then
             cycle
          else
             subgroup%pseudo_class(i)%number = 2*subgroup%klass(j)%number
             class_size = subgroup%pseudo_class(i)%number
             allocate(subgroup%pseudo_class(i)%conjug_el(class_size))
             subgroup%pseudo_class(i)%conjug_el(1:class_size/2) = subgroup%klass(j)%conjug_el
             subgroup%pseudo_class(i)%conjug_el(class_size/2+1:class_size) = &
                  &subgroup%klass(subgroup%klass(j)%inverse)%conjug_el
             i=1+1
          end if
       end if
    end do

    ! determine if subgroup is invariant

    subgroup%invariant = .true.

    outest: do j=1,num_ir
       do i=1,group_num_cl
          ! check if subgroup%klass(j) = klass(i)
          if (equal_class(subgroup%pseudo_class(j),klass(i))) then
             subgroup%super_class(j) = i
             cycle outest
          end if
       end do
       subgroup%invariant = .false.
       exit outest
    end do outest


    ! print class list
    if (output_unit > 0) then
       write(output_unit,*)
       write(output_unit,*)"number of elements in subgroup",num_el
       write(output_unit,602) subgroup%name,num_cl,subgroup%num_re
       write(output_unit,603)
       do i = 1, num_cl
          inv_sign = symm_el(subgroup%klass(i)%conjug_el(1))%inverse_sign
          !if ((.not.subgroup%klass(j)%ambivalent_proj).and.(subgroup%klass(j)%ambivalent)) then
          !   inv_sign = -1
          !end if
          write(output_unit,604) i,inv_sign*subgroup%klass(i)%inverse,regular_char(i),&
               &       (symm_el(subgroup%klass(i)%conjug_el(j))%name, j = 1,subgroup%klass(i)%number)
          if (subgroup%klass(i)%ambivalent_proj) then
             write(output_unit,*) "ambivalent"
          else
             write(output_unit,*) "NON ambivalent"
          end if

       end do

       if (subgroup%invariant) then
          write(output_unit,*)
          write(output_unit,*) "   subgroup is invariant"
          write(output_unit,*) "   subgroup has ",num_ir," pseudo classes:"
          write(output_unit,*)
          do i=1,num_ir
             write(output_unit,*) "   pseudo_class (",i,") = group_class(",subgroup%super_class(i),")"
             write(output_unit,605) i,i,&
                  &(symm_el(subgroup%pseudo_class(i)%conjug_el(j))%name, j = 1,subgroup%pseudo_class(i)%number)
          end do
       end if
    endif

602 format(/1x,'subgroup ',a,' has ',i2,' classes (',i2,' are regular).')
603 format(/1x,' their elements are:'   /)
604 format(1x,i3,' (',i2,') ',a,' : ',8a8:,3(/12x,8a8:))
605 format(1x,i3,' (',i2,') : ',8a8:,3(/12x,8a8:))

  contains

    logical function equal_class(class1,class2)
      ! checks if two classes are identical
      ! (class from subgroup and group)
      !
      !------------ formal parameters -------------------------
      type(group_class),intent(in)       :: class1,class2
      ! classes to compare
      !------------ Declaration of local variables ------------
      integer(kind=i4_kind)              :: m,n
      ! counters
      integer(kind=i4_kind)              :: size1,size2
      ! sizes of classes 1 and 2

      ! initialize
      equal_class = .false.

      size1 = class1%number
      size2 = class2%number

      if (size1.eq.size2) then
         equal_class = .true.
         ! only if both classes have the same size they might be identical
         control: do m=1,size1
            do n=1,size1
               if (class1%conjug_el(m).eq.class2%conjug_el(n)) then
                  cycle control
               end if
            end do
            equal_class = .false.
            exit control
         end do control
      end if


    end function equal_class

  end subroutine subgroup_class_gen
  !*************************************************************


  !*************************************************************
  subroutine group_coset_decomp(num_cs, subgroup, coset, coset_trafo)
    !  Purpose: determines the coset decomposition of the group
    !           for a given subgroup
    !           determines the transformation matrix of the cosets
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of dummy variables -----------------
    integer(kind=i4_kind),intent(inout)        :: num_cs
    ! number of cosets
    type(sub_group),intent(in)                 :: subgroup
    ! contains all information of the subgroup
    type(group_coset),intent(out)              :: coset
    ! coset of the group
    integer(kind=i4_kind), intent(out), allocatable :: coset_trafo(:,:,:)
    ! how the cosets transform with respect to group operations
    ! coset_trafo(num_cs,num_cs,group_num_el)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                      :: num_el
    ! number of elements in the subgroup
    integer (i4_kind) :: i, j, k, l, count
    ! counters
    integer(kind=i4_kind),allocatable          :: used(:)
    ! operations which already were associated with a coset

    !------------ Declaration of subroutines used ----------------
    external error_handler

    !------------ Executable code ------------------------------------

    ! |subgroup| is an integer divider of |group|
    num_cs = group_num_el/subgroup%num_el
    num_el = subgroup%num_el
    coset%num_cs = num_cs
    coset%num_el = num_el

    ! allocate local variables
    allocate(used(group_num_el))

    ! allocate coset and coset_trafo
    allocate(coset%elements(num_el,num_cs))
    allocate(coset_trafo(num_cs,num_cs,group_num_el))

    !
    ! now determine cosets
    !
    used = 0
    do i=1,num_cs
       ! look for group operation, which was not associated yet
       do j=1,group_num_el
          if (used(j).eq.0) then
             ! apply operation to subgroup and
             ! fill the corresponding coset
             do l=1,num_el
                coset%elements(l,i) = group_multab(j,&
                     subgroup%elements(l))
                used(coset%elements(l,i)) = 1
             end do
             exit
          end if
       end do
    end do

    ! test if all symmetry elements reside in cosets
    do j=1,group_num_el
       if (used(j).eq.0) then
          write(output_unit,*) "coset decomposition failed for element :",j
          call error_handler("group_coset_decomp: coset decomposition failed")
       end if
    end do

    !
    ! now the transformation matrices are determined
    !

    coset_trafo = 0
    ! loop over all symmetry operations
    do i=1,group_num_el
       ! loop over all cosets
       do j=1,num_cs
          count = group_multab(i,coset%elements(1,j))
          ! in which coset does count reside?
          search: do k=1,num_cs
             do l=1,num_el
                if (count.eq.coset%elements(l,k)) then
                   coset_trafo(k,j,i) = 1
                   exit search
                end if
             end do
          end do search
       end do
    end do

    ! deallocate local variables
    deallocate(used)

  end subroutine group_coset_decomp
  !*************************************************************


  !*************************************************************
  subroutine group_char_gen(group_num_cl)
    !
    !  Purpose: generates character table
    !
    use type_module
    use iounitadmin_module
    implicit none
    integer(i4_kind), intent(in) :: group_num_cl
    !** End of interface *****************************************

    !------------ Declaration of local parameters ----------------
    integer(kind=i4_kind),parameter     :: matz=1
    !------------ Declaration of local variables ---------------------

    ! matrix of CSCO:
    complex(c16_kind) :: csco_mat(group_num_cl, group_num_cl)
    complex(c16_kind) :: eig_vec(group_num_cl, group_num_cl)

    ! real and imaginary part of matrix of CSCO:
    real(r8_kind) :: eig_vec_r(group_num_cl, group_num_cl)
    real(r8_kind) :: eig_vec_i(group_num_cl, group_num_cl)

    ! work array for diagonal matrix:
    real(r8_kind) :: diag(group_num_cl)

    ! work arrays for eigensolver:
    real(r8_kind) :: work1(group_num_cl)
    real(r8_kind) :: work2(group_num_cl)
    real(r8_kind) :: work3(2, group_num_cl)

    ! counters
    integer(i4_kind) :: i, j, k, ierr, stat

    ! matrix of CSCO
    real(r8_kind) :: factor


    ! allocate characters table
    allocate(characters(group_num_cl,group_num_cl), stat=stat)
    if ( stat .ne. 0 ) call error_handler("group_char_gen: allocate characters table failed")

    ! allocate eigenvalues of CSCO(irrep labels)
    allocate(eig_val(group_num_cl), stat=stat)
    if ( stat .ne. 0 ) call error_handler("group_char_gen: allocate eigenvalues of CSCO failed")

    ! build up matrix of CSCO C

    csco_mat = (0.0_r8_kind,0.0_r8_kind)
    diag = klass%number
    diag = sqrt(diag)

    do k=1,3
       if (ncsco(k).eq.0) then
          cycle
       end if

       eig_vec = class_mult(:,ncsco(k),:)

       ! transform C as:
       ! C` = G*C*(G^-1), where G = diag(sqrt(klass(i)%number))
       do j=1,group_num_cl
          do i=1,group_num_cl
             eig_vec(i,j) = -diag(i)/diag(j)*eig_vec(i,j)
          end do
       end do

       if (klass(ncsco(k))%ambivalent) then
          csco_mat = csco_mat + fcsco(k)*eig_vec
       else
          factor = 1/(2*sqrt(2.0_r8_kind))
          csco_mat = csco_mat + &
               & factor*3.0_r8_kind*fcsco(k)*(transpose(eig_vec) + eig_vec) + &
               & factor*fcsco(k)*cmplx(0.0_r8_kind,1.0_r8_kind,r8_kind)*(eig_vec - &
               & transpose(eig_vec))
       end if
    end do

    ! now solve equation:
    ! C`*y=lambda*y

    call ch(group_num_cl,group_num_cl,real(csco_mat),aimag(csco_mat),&
         eig_val,matz,eig_vec_r,eig_vec_i,&
         work1,work2,work3,ierr)

    eig_vec = cmplx(eig_vec_r,eig_vec_i,r8_kind)
    if (ierr.ne.0) then
       write(output_unit,*)"group_char_gen: stopping due to error in eigensolver"
       write(output_unit,*)"IERR = ",ierr
       call error_handler( "group_char_gen: error in eigensolver" )
    end if

    ! obtain characters as
    ! character = sqrt(group_num_cl)*q, where:
    ! q = (G^-1)*y

    factor=sqrt(dble(group_num_el))
    do i=1,group_num_cl
       characters(i,:) = factor*eig_vec(i,:)/diag(i)
    end do

    ! normalize, such that char(E) is real

    do j=1,group_num_cl
       characters(:,j) = characters(:,j)/characters(1,j)&
            *abs(characters(1,j))
    end do

    ! presentation
    if (output_unit > 0) then
       write(output_unit,*) "------------------------------------------------"
       write(output_unit,*) " Determination of Characters from class algebra "
       write(output_unit,*) "------------------------------------------------"
       write(output_unit,*)
       write(output_unit,*) " CSCO used for group ",group_name,": "
       do k=1,3
          if (ncsco(k).ne.0) then
             write(output_unit,*) ' ',fcsco(k),'* ',symm_el(klass(ncsco(k))%conjug_el(1))%name
          end if
       end do
    endif
  end subroutine group_char_gen
  !*************************************************************


  !*************************************************************
  subroutine group_proj_char_gen(group_num_re)
    !
    !  Purpose: generates character table of projective representations
    !
    use type_module
    use iounitadmin_module
    implicit none
    integer(i4_kind), intent(in) :: group_num_re
    !** End of interface *****************************************

    !------------ Declaration of local parameters ----------------
    integer(kind=i4_kind),parameter     :: matz=1
    !------------ Declaration of local variables ---------------------

    ! matrix of CSCO:
    complex(c16_kind) :: csco_mat(group_num_re, group_num_re)
    complex(c16_kind) :: eig_vec(group_num_re, group_num_re)

    ! real and imaginary part of matrix of CSCO:
    real(r8_kind) :: eig_vec_r(group_num_re, group_num_re)
    real(r8_kind) :: eig_vec_i(group_num_re, group_num_re)

    ! work array for diagonal matrix:
    real(r8_kind) :: diag(group_num_re)

    ! work arrays for eigensolver:
    real(r8_kind) :: work1(group_num_re)
    real(r8_kind) :: work2(group_num_re)
    real(r8_kind) :: work3(2, group_num_re)

    ! counters
    integer(i4_kind) :: i, j, k, ierr, error

    ! matrix of CSCO
    real(r8_kind) :: factor

    ! allocate vector subscripts of regular classes
    allocate(regular(group_num_re),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of regular failed")

    ! allocate vector subscripts of regular classes
    allocate(regular_inv(group_num_cl),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of regular_inv failed")

    ! allocate characters table
    allocate(proj_characters(group_num_re,group_num_re),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of proj_characters failed")

    ! allocate eigenvalues of CSCO(irrep labels)
    allocate(proj_eig_val(group_num_re),stat=error)
    if (error.ne.0) call error_handler( &
         "symmel_alloc: allocation of proj_eig_val failed")

    ! generate subscript list of regular classes
    i = 0
    ! preset inverse of regular to identity
    regular_inv = 1
    do k=1,group_num_cl
       if (klass(k)%regular) then
          i = i+1
          regular(i) = k
          regular_inv(k) = i
       end if
    end do

    ! build up matrix of CSCO C

    csco_mat = (0.0_r8_kind,0.0_r8_kind)
    diag = klass(regular)%number
    diag = sqrt(diag)

    do k=1,3
       if (ncsco_proj(k).eq.0) then
          cycle
       end if

       eig_vec = class_mult(regular,ncsco_proj(k),regular)*class_fact(regular,ncsco_proj(k),regular)

       ! transform C as:
       ! C` = G*C*(G^-1), where G = diag(sqrt(klass(i)%number))
       do j=1,group_num_re
          do i=1,group_num_re
             eig_vec(i,j) = -diag(i)/diag(j)*eig_vec(i,j)
          end do
       end do

       if (klass(ncsco_proj(k))%ambivalent_proj) then
          csco_mat = csco_mat + fcsco_proj(k)*eig_vec
       else
          factor = 1/(2*sqrt(2.0_r8_kind))
          csco_mat = csco_mat + &
               & factor*3.0_r8_kind*fcsco_proj(k)*(transpose(eig_vec) + eig_vec) + &
               & factor*fcsco_proj(k)*cmplx(0.0_r8_kind,1.0_r8_kind,r8_kind)*(eig_vec - &
               & transpose(eig_vec))
       end if
    end do

    ! now solve equation:
    ! C`*y=lambda*y

    call ch(group_num_re,group_num_re,real(csco_mat),aimag(csco_mat),&
         proj_eig_val,matz,eig_vec_r,eig_vec_i,&
         work1,work2,work3,ierr)

    eig_vec = cmplx(eig_vec_r,eig_vec_i,r8_kind)
    if (ierr.ne.0) then
       write(*,*) "group_proj_char_gen: stopped due to error in eigensolver, IERR = ",ierr
       call error_handler("group_proj_char_gen: stopped due to error in eigensolver")
    end if

    ! obtain characters as
    ! character = sqrt(group_num_cl)*q, where:
    ! q = (G^-1)*y

    factor=sqrt(dble(group_num_el))
    do i=1,group_num_re
       proj_characters(i,:) = factor*eig_vec(i,:)/diag(i)
    end do

    ! normalize, such that char(E) is real

    do j=1,group_num_re
       proj_characters(:,j) = proj_characters(:,j)/proj_characters(1,j)&
            *abs(proj_characters(1,j))
    end do

    ! presentation
    if (output_unit > 0) then
       write(output_unit,*) "------------------------------------------------"
       write(output_unit,*) " projective representations: "
       write(output_unit,*) " Determination of Characters from class algebra "
       write(output_unit,*) "------------------------------------------------"
       write(output_unit,*)
       write(output_unit,*) " CSCO used for group # ",igroup," : "
       do k=1,3
          if (ncsco_proj(k).ne.0) then
             write(output_unit,*) ncsco_proj(k),' ',fcsco_proj(k),'* ',symm_el(klass(ncsco_proj(k))%conjug_el(1))%name
          end if
       end do

       ! show results of calculation
       write(output_unit,601) group_name, (symm_el(klass(regular(j))%conjug_el(1))&
            &%name,&
            &          j = 1,group_num_re)
       do j = 1,group_num_re
          write(output_unit,602) proj_eig_val(j), (proj_characters(k,j),&
               &          k = 1,group_num_re)
       end do
    endif

601 format(/1x,'Group ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,f7.2,2x,4(f7.2,f7.2:)/4(10x,4(f7.2,f7.2:)/))
  end subroutine group_proj_char_gen
  !*************************************************************


  !*************************************************************
  subroutine subgroup_char_gen(subgroup)
    !  Purpose: generates character table of subgroup
    !           note, that the subgroup must be invariant, such that
    !           its classes are a subgroup of the groups classes!
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    type(sub_group)                         :: subgroup
    ! contains all information about the subgroup
    !** End of interface *****************************************
    !------------ Declaration of local constants -----------------
    real(kind=r8_kind),parameter  :: small = 1e-5

    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: num_ir
    ! number of pseudo irreps in subgroup
    integer(kind=i4_kind)                :: i,j,k,l
    ! counters
    real(kind=r8_kind),allocatable       :: work1(:)
    ! work array

    num_ir = subgroup%num_ir

    ! allocate work arrays
    allocate(work1(num_ir))

    ! allocate character table of the subgroup
    allocate(subgroup%characters(num_ir,num_ir))

    if (.not.subgroup%invariant) then
       call error_handler("subgroup_char_gen: subgroup is not invariant")
    end if

    ! pseudo irrep number l
    l = 1
    ! loop over all irreps of the (super) group
    outer:do j=1,group_num_cl
       do i=1,num_ir
          work1(i) = characters(subgroup%super_class(i),j)
       end do
       ! look if irrep already appeared
       do k=1,j-1
          if (maxval(abs(subgroup%characters(:,k)-work1)).lt.small) then
             ! if yes take next irrep
             cycle outer
          end if
       end do
       subgroup%characters(:,l) = work1(1:num_ir)
       l = l + 1
    end do outer

    ! show results of calculation
    if (output_unit > 0) then
       write(output_unit,*) "number of pseudo classes in subgroup ",num_ir
       write(output_unit,601) subgroup%name, (symm_el(subgroup%pseudo_class(j)%conjug_el(1))&
            &%name,&
            &          j = 1,num_ir)
       do j = 1,num_ir
          write(output_unit,602) (subgroup%characters(k,j),&
               &          k = 1,num_ir)
       end do
    endif

601 format(/1x,'subgroup ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,4(f17.12,f17.12:)/4(10x,4(f17.12,f17.12:)/))


    ! deallocate working arrays
    deallocate(work1)

  end subroutine subgroup_char_gen
  !*************************************************************


  !*************************************************************
  subroutine group_generator_gen
    !  Purpose: generates generator of group
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Definition of local types ----------------------
    integer(kind=i4_kind)  :: j
    ! counter
    character(len=8)       :: symm_op

    ! initialize symmetry elements
    inversion = 0
    horizontal = 0
    primary_axis = 0
    secondary_axis = 0
    main_axis = 0
    order_primary_axis = 0
    order_secondary_axis = 0

    !
    ! check presence of inversion
    !

    inversion = present('i       ')
    horizontal = present('sigh    ')

    !
    ! determine primary axis
    !

    if ((inversion.eq.0).and.(horizontal.eq.0)) then
       ! case: no inversion or horizontal reflection plane is present
       if (present('C5,1+   ').ne.0) then
          ! icoshedral groups
          primary_axis = present('C5,1+   ')
          order_primary_axis = 5
          elseif(present('C3,1+   ').ne.0) then
             ! octahedral groups
          primary_axis = present('C3,1+   ')
          order_primary_axis = 3
       else
          ! axial groups
          ! look for maximal improper axis
          do j=1,4
             write(symm_op,'(A,I1,A)') 'S',2*j,'+'
             if (present(symm_op).ne.0) then
                primary_axis = present(symm_op)
                order_primary_axis = 2*j
             end if
          end do
          do j=5,10
             write(symm_op,'(A,I2,A)') 'S',2*j,'+'
             if (present(symm_op).ne.0) then
                primary_axis = present(symm_op)
                order_primary_axis = 2*j
             end if
          end do

          if (primary_axis.eq.0) then
             ! if no improper axis was found then
             ! look for maximal proper axis
             if (present('C2      ').ne.0) then
                primary_axis = present('C2      ')
                order_primary_axis = 2
             end if
             do j=3,9
                write(symm_op,'(A,I1,A)') 'C',j,'+'
                if (present(symm_op).ne.0) then
                   primary_axis = present(symm_op)
                   order_primary_axis = j
                end if
             end do
             if (present('C10+    ').ne.0) then
                primary_axis = present('C10+    ')
                order_primary_axis = 10
             end if

          end if
       end if
    else
       ! case: inversion or horizontal reflection plane is present
       if (present('C5,1+   ').ne.0) then
          ! icosahedral groups
          primary_axis = present('C5,1+   ')
          order_primary_axis = 5
          elseif(present('C3,1+   ').ne.0) then
             ! octahedral groups
          primary_axis = present('C3,1+   ')
          order_primary_axis = 3
       else
          ! axial groups
          if (present('C2      ').ne.0) then
             primary_axis = present('C2      ')
             order_primary_axis = 2
          end if
          do j=3,9
             write(symm_op,'(A,I1,A)') 'C',j,'+'
             if (present(symm_op).ne.0) then
                primary_axis = present(symm_op)
                order_primary_axis = j
             end if
          end do
          if (present('C10+    ').ne.0) then
             primary_axis = present('C10+    ')
             order_primary_axis = 10
          end if
       end if
    end if

    !
    ! determine secondary axis
    !

    if ((inversion.eq.0).and.(horizontal.eq.0)) then
       ! case: no inversion or horizontal reflection plane is present
       if (present('S4x+    ').ne.0) then
          secondary_axis = present('S4x+    ')
          order_secondary_axis = 4
          elseif (present('C4x+    ').ne.0) then
          secondary_axis = present('C4x+    ')
          order_secondary_axis = 4
          elseif (present('C2x     ').ne.0) then
          secondary_axis = present('C2x     ')
          order_secondary_axis = 2
          elseif (present('C2,1''   ').ne.0) then
          secondary_axis = present('C2,1''   ')
          order_secondary_axis = 2
          elseif (present('sigv1   ').ne.0) then
          secondary_axis = present('sigv1   ')
          order_secondary_axis = 1
          elseif (present('sigv1''  ').ne.0) then
          secondary_axis = present('sigv1''  ')
          order_secondary_axis = 1
       end if
    else
       ! case: inversion or horizontal reflection plane is present
       if (present('C4x+    ').ne.0) then
          secondary_axis = present('C4x+    ')
          order_secondary_axis = 4
          elseif (present('C2x     ').ne.0) then
          secondary_axis = present('C2x     ')
          order_secondary_axis = 2
          elseif (present('C2,1''   ').ne.0) then
          secondary_axis = present('C2,1''   ')
          order_secondary_axis = 2
          elseif (present('sigv1   ').ne.0) then
          secondary_axis = present('sigv1   ')
          order_secondary_axis = 1
          elseif (present('sigv1''  ').ne.0) then
          secondary_axis = present('sigv1''  ')
          order_secondary_axis = 1
       end if
    end if

    ! special cases
    if (group_name.eq.'C2V ') then
       secondary_axis = present('sigx    ')
       order_secondary_axis = 1
    end if

    if (output_unit > 0) then
       write(output_unit,*)
       write(output_unit,*) "group ",group_name," has the generators: "
       write(output_unit,*)
    endif
    if (inversion.ne.0) then
       if (output_unit > 0) then
          write(output_unit,*)  "inversion ",symm_el(inversion)%name
       endif
       inversion  = symm_el(inversion)%klass
    else
       if (output_unit > 0) then
          write(output_unit,*)  "no inversion"
       endif
    endif

    if (horizontal.ne.0) then
       if (output_unit > 0) then
          write(output_unit,*)  "horizontal plane ",symm_el(horizontal)%name
       endif
       horizontal = symm_el(horizontal)%klass
    else
       if (output_unit > 0) then
          write(output_unit,*)  "no horizontal plane"
       endif
    endif

    if (primary_axis.ne.0) then
       if (output_unit > 0) then
          write(output_unit,*) "primary axis: ", symm_el(primary_axis)%name
          write(output_unit,*)  "order of primary axis ",order_primary_axis
       endif
       primary_axis = symm_el(primary_axis)%klass
    else
       if (output_unit > 0) then
          write(output_unit,*)  "no primary axis"
       endif
    endif

    if (secondary_axis.ne.0) then
       if (output_unit > 0) then
          write(output_unit,*)  "secondary axis: ",symm_el(secondary_axis)%name
          write(output_unit,*)  "order of secondary_axis ",order_secondary_axis
       endif
       secondary_axis = symm_el(secondary_axis)%klass
    else
       if (output_unit > 0) then
          write(output_unit,*)  "no secondary axis"
       endif
    endif

    if (order_secondary_axis.ge.order_primary_axis) then
       main_axis = secondary_axis
    else
       main_axis = primary_axis
    endif

    if (main_axis.ne.0) then
       if (output_unit > 0) then
          write(output_unit,*)  "main axis (axis of maximum order) ",symm_el(klass(main_axis)%conjug_el(1))%name
       endif
    endif



  contains

    integer function present(operation)
      !
      ! Purpose: checks if >>operation<< is a symmetry element
      ! of the group
      !
      !------------ Declaration of local variables ----------------
      character(len=8), intent(in) :: operation
      ! operation to look for
      integer(kind=i4_kind)  :: j
      ! counter

      present = 0
      do j=1,group_num_el
         if (operation.eq.symm_el(j)%name) then
            present = j
            exit
         end if
      end do

    end function present

  end subroutine group_generator_gen
  !*************************************************************


  !*************************************************************
  subroutine subgroup_generator_gen(num_el,elements,subgroup)
    !  Purpose: checks which generator elements of the group are left
    !           in the subgroup
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    implicit none
    integer(kind=i4_kind), intent(in)       :: num_el
    ! number of elements of the subgroup
    integer(kind=i4_kind), intent(in)       :: elements(:)
    ! elements of the subgroup
    type(sub_group)                         :: subgroup
    ! contains all information of the subgroup
    !** End of interface *****************************************
    !------------ Definition of local types ----------------------
    integer(kind=i4_kind)  :: j
    ! counter

    ! initialize symmetry elements
    subgroup%inversion = 0
    subgroup%horizontal = 0
    subgroup%primary_axis = 0
    subgroup%secondary_axis = 0
    subgroup%main_axis = 0
    subgroup%order_primary_axis = 0
    subgroup%order_secondary_axis = 0

    !
    ! check presence of inversion
    !

    if (inversion.ne.0) then
       if (present_int(klass(inversion)%conjug_el(1)).ne.0) then
          subgroup%inversion = inversion
          if (output_unit > 0) then
             write(output_unit,*) "Inversion still present"
             write(output_unit,*)  symm_el(inversion)%name
          endif
       endif
    endif

    !
    ! check presence of horizontal plane
    !

    if (horizontal.ne.0) then
       if (present_int(klass(horizontal)%conjug_el(1)).ne.0) then
          subgroup%horizontal = horizontal
          if (output_unit > 0) then
             write(output_unit,*) "Horizontal plane still present"
             write(output_unit,*)  symm_el(horizontal)%name
          endif
       endif
    endif

    !
    ! is primary axis still present?
    !

    if (primary_axis.ne.0) then
       do j=1,klass(primary_axis)%number
          if (present_int(klass(primary_axis)%conjug_el(j)).ne.0) then
             subgroup%primary_axis = primary_axis
             subgroup%order_primary_axis = order_primary_axis
             if (output_unit > 0) then
                write(output_unit,*) "Primary axis is still present"
                write(output_unit,*)  symm_el(klass(primary_axis)%conjug_el(1))%name
             endif
             exit
          endif
       end do
    endif

    !
    ! is secondary_axis axis still present?
    !

    if (secondary_axis.ne.0) then
       do j=1,klass(secondary_axis)%number
          if (present_int(klass(secondary_axis)%conjug_el(j)).ne.0) then
             subgroup%secondary_axis = secondary_axis
             if (output_unit > 0) then
                write(output_unit,*) "Secondary axis is still present"
                write(output_unit,*)  symm_el(klass(secondary_axis)%conjug_el(1))%name
             endif
             subgroup%order_secondary_axis = order_secondary_axis
             exit
             elseif (present('C2x     ').ne.0) then
             subgroup%secondary_axis = symm_el(present('C2x     '))%klass
             subgroup%order_secondary_axis = 2
             if (output_unit > 0) then
                write(output_unit,*) "Secondary axis was replaced"
                write(output_unit,*)  symm_el(klass(subgroup%secondary_axis)%conjug_el(1))%name
             endif
             exit
          endif
       end do
    endif

    !
    ! is main_axis axis still present?
    !

    if (main_axis.ne.0) then
       do j=1,klass(main_axis)%number
          if (present_int(klass(main_axis)%conjug_el(j)).ne.0) then
             subgroup%main_axis = main_axis
             if (output_unit > 0) then
                write(output_unit,*) "Main axis is still present"
                write(output_unit,*)  symm_el(klass(main_axis)%conjug_el(1))%name
             endif
             exit
          endif
       end do
    endif



  contains

    integer function present_int(operation)
      !
      ! Purpose: checks if >>operation<< is a symmetry element
      ! of the group
      !
      !------------ Declaration of local variables ----------------
      integer,intent(in)     :: operation
      ! operation to look for
      integer(kind=i4_kind)  :: j
      ! counter

      present_int = 0
      do j=1,num_el
         if (operation.eq.elements(j)) then
            present_int = elements(j)
            exit
         end if
      end do

    end function present_int

    integer function present(operation)
      !
      ! Purpose: checks if >>operation<< is a symmetry element
      ! of the group
      !
      !------------ Declaration of local variables ----------------
      character(len=8), intent(in) :: operation
      ! operation to look for
      integer(kind=i4_kind)  :: j
      ! counter

      present = 0
      do j=1,group_num_el
         if (operation.eq.symm_el(j)%name) then
            present = j
            exit
         end if
      end do

    end function present

  end subroutine subgroup_generator_gen
  !*************************************************************


  !*************************************************************
  subroutine group_irrep_label
    !  Purpose: labels irreps using their character table
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module

    implicit none

    !------------ Definition of local types ----------------------
    type irrep
       integer(kind=i4_kind)  :: irrep_number
       real(kind=r8_kind)     :: main_axis_value
       logical                :: flag
       type(irrep), pointer   :: next => NULL()
    end type irrep
    ! a special irrep

    type irrep_list
       type(irrep), pointer :: start => NULL()
       type(irrep), pointer :: end => NULL()
    end type irrep_list
    ! list of a group of irreps of the same size
    ! and inversion symmetry

    !------------ Definition of local constants  -----------------

    real(kind=r8_kind),parameter  :: small = 1e-5
    integer(kind=i4_kind),parameter,dimension(4) :: a_vec = (/1,1,1,1/)
    integer(kind=i4_kind),parameter,dimension(4) :: b1_vec = (/1,1,-1,-1/)
    integer(kind=i4_kind),parameter,dimension(4) :: b2_vec = (/1,-1,-1,1/)
    integer(kind=i4_kind),parameter,dimension(4) :: b3_vec = (/1,-1,1,-1/)

    !------------ Definition of local variables  -----------------

    integer(kind=i4_kind)       :: i,j,k,i_lab,inv_ind,dim_ind, error
    ! counter

    integer(kind=i4_kind)       :: inv_hor
    ! is inversion or horizontal plane present

    logical                     :: single
    ! is .true. if there is only one irrep of a special kind

    character              :: even_char, odd_char
    ! even and odd respective inversion or horizontal plane

    type(irrep), pointer       :: new_irrep,irr_point,previous
    ! newly generated irrep, pointer to some irrep

    type(irrep), pointer       :: irr_point2,irr_pointnew,irr_pointold
    ! newly generated irrep, pointer to some irrep

    type(irrep_list)           :: irrep_kinds(2,0:5)
    ! irrep_kinds points to irrep list of a special kind of
    ! irreps, which is specified by inversion symmetry and
    ! dimension

    !-------- Executable Code ------------------------------------

    ! check existence of inversion and horizontal plane
    if (inversion.ne.0) then
       inv_hor = inversion
       even_char = 'G'
       odd_char  = 'U'
    elseif (horizontal.ne.0) then
       inv_hor = horizontal
       even_char = ''''
       odd_char  = '"'
    else
       inv_hor = 0
       ! these are unused in this case:
       even_char = '?'
       odd_char = '?'
    endif

    ! initialize number of irreps
    group_num_ir = 2*group_num_cl

    ! loop over all irreps
    do i = 1, group_num_cl

       ! process first character of label
       i_lab = 1

       ! generate new irrep
       allocate(new_irrep)
       new_irrep%irrep_number = i
       new_irrep%flag = .false.

       ! determine inversion symmetry of irrep
       inv_ind = 1
       if (inv_hor.ne.0) then
          if (real(characters(inv_hor,i)).ge.0) then
             inv_ind = 1
          else
             inv_ind = 2
          end if
       end if

       ! determine dimension of irrep
       dim_ind = int(real(characters(1,i))+small)

       if ((group_name.eq.'D2  ').or.(group_name.eq.'D2H ')) then
          ! special case of groups D2 and D2H
          if (maxval(abs(characters(1:4,i)-a_vec)).lt.small) then
             dim_ind = 0
             new_irrep%main_axis_value = 1
          elseif (maxval(abs(characters(1:4,i)-b1_vec)).lt.small) then
             dim_ind = 1
             new_irrep%main_axis_value = 3
          elseif (maxval(abs(characters(1:4,i)-b2_vec)).lt.small) then
             dim_ind = 1
             new_irrep%main_axis_value = 2
          elseif (maxval(abs(characters(1:4,i)-b3_vec)).lt.small) then
             dim_ind = 1
             new_irrep%main_axis_value = 1
          endif

       else
          ! normal case
          ! is one dimensional irrep A or B?
          if (dim_ind.eq.1) then
             ! A irrep
             dim_ind = 0
             if (primary_axis.ne.0) then
                if (abs(aimag(characters(primary_axis,i))).gt.small) then
                   ! case of nonambivalent group, where conjugated
                   ! irrep pairs appear, which form a two dim. irrep
                   dim_ind = 2
                   ! decrement number of irreps
                   group_num_ir = group_num_ir - 1
                elseif (real(characters(primary_axis,i)).lt.0) then
                   ! B irrep
                   dim_ind = 1
                endif
             endif
          endif

          ! determine by which value irreps are ordered
          if (dim_ind.gt.1) then
             if (main_axis.ne.0) then
                new_irrep%main_axis_value = real(characters(main_axis,i))
             end if
          else
             ! treatment of A and B
             if (secondary_axis.ne.0) then
                new_irrep%main_axis_value = real(characters(secondary_axis,i))
             end if
          end if
       endif


       if (associated(irrep_kinds(inv_ind,dim_ind)%start)) then
          ! extend old list
          irrep_kinds(inv_ind,dim_ind)%end%next => new_irrep
          irrep_kinds(inv_ind,dim_ind)%end => new_irrep
       else
          ! make new list
          irrep_kinds(inv_ind,dim_ind)%start => new_irrep
          ! first irrep exists double, because interchange is simpler
          allocate(new_irrep)
          new_irrep = irrep_kinds(inv_ind,dim_ind)%start
          irrep_kinds(inv_ind,dim_ind)%start%main_axis_value = 20.0
          nullify(new_irrep%next)
#ifdef _COMPAC_FORTRAN
          new_irrep%main_axis_value=0.0_r8_kind !!!!!!!!!!!!!!!!!!!AS
#endif
          irrep_kinds(inv_ind,dim_ind)%start%next => new_irrep
          irrep_kinds(inv_ind,dim_ind)%end => new_irrep
       end if

    end do

    ! determine number of irreps
    group_num_ir = group_num_ir/2

    ! allocate canonically ordered irreps
    allocate(irrep_can(group_num_ir))

    ! initialize irrep_can, which contains all information
    ! about the irreps (in canonical order)
    do i=1,group_num_ir
       allocate(irrep_can(i)%characters(group_num_cl))
       irrep_can(i)%label = '        '
       irrep_can(i)%pseudo = .false.
    enddo

    ! irrep number
    i = 1
    ! process all kinds of irreps
    do inv_ind = 1,2
       do dim_ind = 0,5
             if (.not.associated(irrep_kinds(inv_ind,dim_ind)%start)) cycle
             ! irrep_kinds(inv_ind,dim_ind)%end%flag = .true.
             ! sort irreps of one kind
             do
                irr_pointold => irrep_kinds(inv_ind,dim_ind)%start
                irr_point => irr_pointold%next
                single = .false.
                if (.not.associated(irr_point%next)) then
                   single = .true.
                   exit
                   elseif (irr_point%next%flag) then
                   exit
                end if
                do
                   if (irr_point%main_axis_value.lt.irr_point%next%&
                        &main_axis_value) then
                      irr_point2        => irr_point%next
                      irr_pointnew      => irr_point%next%next

                      irr_point%next%next =>irr_point
                      irr_point%next => irr_pointnew
                      irr_pointold%next => irr_point2

                      irr_point         => irr_point2
                   end if
                   irr_pointold => irr_point
                   irr_point => irr_point%next
                   if (.not.associated(irr_point%next)) then
                      irr_point%flag = .true.
                      exit
                      elseif (irr_point%next%flag) then
                      irr_point%flag = .true.
                      exit
                   end if
                end do
             end do
             ! number j of multi dimensional irrep
             ! e.g. E1, E2
             j=1
             previous => irrep_kinds(inv_ind,dim_ind)%start
             irr_point => irrep_kinds(inv_ind,dim_ind)%start%next
             do while ( associated(irr_point) )
                   ! write character table and eigenvalues according to
                   ! the canonical order of the irreps

                   ! check for pseudo-2D irreps
                   if (abs(irr_point%main_axis_value-previous%&
                        &main_axis_value).lt.small) then
                      j = j - 1 ! decrement number of multi
                      i = i - 1 ! decrement number of irrep
                      irrep_can(i)%pseudo = .true.
                      ! set time inversion invariant dimension of irrep
                      irrep_can(i)%dimension = 2
                      irrep_can(i)%time_dimension = 1
                      ! check if only one pseudo-2D irrep occurs
                      if ((.not.associated(irr_point%next)).and.(j.eq.1)) then
                         single = .true.
                         irrep_can(i)%label(2:3) = '  '
                      endif
                   else
                      ! set dimensions
                      if (dim_ind.eq.0) then
                         irrep_can(i)%dimension = 1
                         irrep_can(i)%time_dimension = 1
                      else
                         irrep_can(i)%dimension = dim_ind
                         irrep_can(i)%time_dimension = dim_ind
                      endif
                   end if
                   irrep_can(i)%characters = characters(:,irr_point%irrep_number)
                   irrep_can(i)%eigen_value = eig_val(irr_point%irrep_number)
                   !char_tab(:,i)  = characters(:,irr_point%irrep_number)
                   !eig_val_can(i) = eig_val(irr_point%irrep_number)
                   select case(dim_ind)
                   case(0)
                      irrep_can(i)%label(1:1) = 'A'
                   case(1)
                      irrep_can(i)%label(1:1) = 'B'
                   case(2)
                      irrep_can(i)%label(1:1) = 'E'
                   case(3)
                      irrep_can(i)%label(1:1) = 'T'
                   case(4)
                      irrep_can(i)%label(1:1) = 'G'
                   case(5)
                      irrep_can(i)%label(1:1) = 'H'
                   end select
                   if (.not.single) then
                      write(irrep_can(i)%label(2:2),'(I1)') j
                   end if
                   if (inv_hor.ne.0) then
                      if (.not.single) then
                         if (real(characters(inv_hor,irr_point%irrep_number)).ge.0)&
                              & then
                            irrep_can(i)%label(3:3) = even_char
                         else
                            irrep_can(i)%label(3:3) = odd_char
                         endif
                      else
                         if (real(characters(inv_hor,irr_point%irrep_number)).ge.0)&
                              & then
                            irrep_can(i)%label(2:2) = even_char
                         else
                            irrep_can(i)%label(2:2) = odd_char
                         endif

                      end if
                   endif
                   j = j + 1 ! increment number of multi
                   i = i + 1 ! increment irrep number
                   previous  => irr_point
                   irr_point => irr_point%next
             end do
       end do! dimensions
    end do! inv_hor

    ! Cleaning list of all kinds of irreps
    DPRINT "Cleaning list of all kinds of irreps"
    !call mmon_get_mem_info("group_irrep_label : Cleaning list")
    do inv_ind = 1,2
       do dim_ind = 0,5
          ! deallocate list starting at irrep_kinds(inv_ind, dim_ind)%start ...

          ! start with the head:
          irr_point => irrep_kinds(inv_ind, dim_ind)%start

          do while ( associated(irr_point) )
             ! save the pointer to the tail:
             new_irrep => irr_point%next

             ! deallocate head:
             deallocate(irr_point, stat=error)
             if( error /= 0) call error_handler("group_irrep_label : Deallocation failed")

             ! advance head:
             irr_point => new_irrep
          enddo
       end do
    end do
    DPRINT "Cleaning list of all kinds of irreps, Done"
    !call mmon_get_mem_info("group_irrep_label : Cleaning list, Done")

    if (output_unit > 0) then
       write(output_unit,601) group_name, (symm_el(klass(j)%conjug_el(1))&
            &%name,&
            &          j = 1,group_num_cl)
       do j = 1,group_num_ir
          write(output_unit,602) irrep_can(j)%label, (irrep_can(j)%characters(k),&
               &          k = 1,group_num_cl)
       end do
    endif

601 format(/1x,'Group ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,A3,2x,4(f8.3,f8.3:)/4(10x,4(f8.3,f8.3:)/))


  end subroutine group_irrep_label
  !*************************************************************


  !*************************************************************
  subroutine group_proj_irrep_label
    !  Purpose: labels irreps using their character table
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module

    implicit none

    !------------ Definition of local types ----------------------
    type irrep
       integer(kind=i4_kind)  :: irrep_number
       real(kind=r8_kind)     :: main_axis_value
       integer(kind=i4_kind)  :: dimension
       logical                :: flag
       type(irrep),pointer    :: next => NULL()
    end type irrep
    ! a special irrep

    type irrep_list
       type(irrep),pointer  :: start => NULL()
       type(irrep),pointer  :: end => NULL()
    end type irrep_list
    ! list of a group of irreps of the same size
    ! and inversion symmetry

    !------------ Definition of local constants  -----------------

    real(kind=r8_kind),parameter  :: small = 1e-5

    !------------ Definition of local variables  -----------------

    integer(kind=i4_kind)       :: i,j,k,i_lab,inv_ind,dim_ind
    ! counter

    integer(kind=i4_kind)       :: inv_hor
    ! is inversion or horizontal plane present

    logical                     :: single
    ! is .true. if there is only one irrep of a special kind

    character              :: even_char, odd_char
    ! even and odd respective inversion or horizontal plane

    type(irrep), pointer       :: new_irrep,irr_point,previous
    ! newly generated irrep, pointer to some irrep

    type(irrep), pointer       :: irr_point2,irr_pointnew,irr_pointold
    ! newly generated irrep, pointer to some irrep

    type(irrep_list)           :: irrep_kinds(2)
    ! irrep_kinds points to irrep list of a special kind of
    ! irreps, which is specified by inversion symmetry and
    ! dimension

    !-------- Executable Code ------------------------------------

    ! check existence of inversion and horizontal plane
    if (inversion.ne.0) then
       inv_hor = inversion
       even_char = 'G'
       odd_char  = 'U'
    elseif (horizontal.ne.0) then
       inv_hor = horizontal
       even_char = ''''
       odd_char  = '"'
    else
       inv_hor = 0
       ! these are unused in this case:
       even_char = '?'
       odd_char = '?'
    endif

    ! initialize number of irreps
    group_num_pir = 2*group_num_re

    ! loop over all irreps
    do i = 1, group_num_re

       ! process first character of label
       i_lab = 1

       ! generate new irrep
       allocate(new_irrep)
       new_irrep%irrep_number = i
       new_irrep%flag = .false.

       ! determine inversion symmetry of irrep
       inv_ind = 1
       if (inv_hor.ne.0) then
          if (real(proj_characters(regular_inv(inv_hor),i)).ge.0) then
             inv_ind = 1
          else
             inv_ind = 2
          end if
       end if

       ! determine dimension of irrep
       dim_ind = int(real(proj_characters(1,i))+small)
       new_irrep%dimension = dim_ind

       ! check for conjugated irrep pair
       do j=1,group_num_re
          if (abs(aimag(proj_characters(j,i))).gt.small) then
             ! case of nonambivalent group, where conjugated
             ! irrep pairs appear, which form a two dim. irrep
             ! decrement number of irreps
             group_num_pir = group_num_pir - 1
             ! adjust dimension
             dim_ind = 2*dim_ind
             exit
          end if
       end do

       ! determine by which value irreps are ordered
       if (main_axis.ne.0) then
          new_irrep%main_axis_value = real(proj_characters(regular_inv(main_axis),i))/&
               new_irrep%dimension
       else if (inv_hor.ne.0) then
          new_irrep%main_axis_value = real(proj_characters(regular_inv(inv_hor),i))/&
               new_irrep%dimension
       end if


       if (associated(irrep_kinds(inv_ind)%start)) then
          ! extend old list
          irrep_kinds(inv_ind)%end%next => new_irrep
          irrep_kinds(inv_ind)%end => new_irrep
       else
          ! make new list
          irrep_kinds(inv_ind)%start => new_irrep
          ! first irrep exists double, because interchange is simpler
          allocate(new_irrep)
          new_irrep = irrep_kinds(inv_ind)%start
          irrep_kinds(inv_ind)%start%main_axis_value = 20.0
          nullify(new_irrep%next)
          irrep_kinds(inv_ind)%start%next => new_irrep
          irrep_kinds(inv_ind)%end => new_irrep
       end if

    end do

    ! determine number of irreps
    group_num_pir = group_num_pir/2

    ! allocate canonically ordered irreps
    allocate(proj_irrep_can(group_num_pir))

    ! initialize irrep_can, which contains all information
    ! about the irreps (in canonical order)
    do i=1,group_num_pir
       allocate(proj_irrep_can(i)%characters(group_num_re))
       proj_irrep_can(i)%label = '        '
       proj_irrep_can(i)%pseudo = .false.
    enddo

    ! irrep number
    i = 1
    ! process all kinds of irreps
    do inv_ind = 1,2
          ! number j of multi dimensional irrep
          ! e.g. E1, E2
          j = 1
          if (.not.associated(irrep_kinds(inv_ind)%start)) cycle
          ! irrep_kinds(inv_ind,dim_ind)%end%flag = .true.
          ! sort irreps of one kind
          do
             irr_pointold => irrep_kinds(inv_ind)%start
             irr_point => irr_pointold%next
             single = .false.
             if (.not.associated(irr_point%next)) then
                single = .true.
                exit
                elseif (irr_point%next%flag) then
                exit
             end if
             do
                if (irr_point%main_axis_value.lt.irr_point%next%&
                     &main_axis_value) then
                   irr_point2        => irr_point%next
                   irr_pointnew      => irr_point%next%next

                   irr_point%next%next =>irr_point
                   irr_point%next => irr_pointnew
                   irr_pointold%next => irr_point2

                   irr_point         => irr_point2
                end if
                irr_pointold => irr_point
                irr_point => irr_point%next
                if (.not.associated(irr_point%next)) then
                   irr_point%flag = .true.
                   exit
                   elseif (irr_point%next%flag) then
                   irr_point%flag = .true.
                   exit
                end if
             end do
          end do
          previous => irrep_kinds(inv_ind)%start
          irr_point => irrep_kinds(inv_ind)%start%next
          do while ( associated(irr_point) )
                dim_ind = irr_point%dimension
                ! write character table and eigenvalues according to
                ! the canonical order of the irreps

                ! check for pseudo-2D irreps
                if (abs(irr_point%main_axis_value-previous%&
                     &main_axis_value).lt.small) then
                   j = j - 1 ! decrement number of multi
                   i = i - 1 ! decrement number of irrep
                   proj_irrep_can(i)%pseudo = .true.
                   ! set time inversion invariant dimension of irrep
                   proj_irrep_can(i)%dimension = 2*dim_ind
                   proj_irrep_can(i)%time_dimension = dim_ind
                   dim_ind = 2*dim_ind
                   ! check if only one pseudo-2D irrep occurs
                else
                   ! set dimensions
                   proj_irrep_can(i)%dimension = dim_ind
                   proj_irrep_can(i)%time_dimension = dim_ind
                end if
                proj_irrep_can(i)%characters = proj_characters(:,irr_point%irrep_number)
                proj_irrep_can(i)%eigen_value = proj_eig_val(irr_point%irrep_number)
                !char_tab(:,i)  = characters(:,irr_point%irrep_number)
                !eig_val_can(i) = eig_val(irr_point%irrep_number)
                select case(dim_ind)
                case(1)
                   proj_irrep_can(i)%label(1:1) = 'A'
                case(2)
                   proj_irrep_can(i)%label(1:1) = 'E'
                case(3)
                   proj_irrep_can(i)%label(1:1) = 'T'
                case(4)
                   proj_irrep_can(i)%label(1:1) = 'G'
                case(5)
                   proj_irrep_can(i)%label(1:1) = 'H'
                case(6)
                   proj_irrep_can(i)%label(1:1) = 'I'
                end select
                write(proj_irrep_can(i)%label(2:4),'(I1,a2)') 2*j-1,'/2'
                if (inv_hor.ne.0) then
                   if (real(proj_characters(regular_inv(inv_hor),irr_point%irrep_number)).ge.0)&
                        & then
                      proj_irrep_can(i)%label(5:5) = even_char
                   else
                      proj_irrep_can(i)%label(5:5) = odd_char
                   endif
                endif
                j = j + 1 ! increment number of multi
                i = i + 1 ! increment irrep number
                previous  => irr_point
                irr_point => irr_point%next
          end do
    end do! inv_hor

    if (output_unit > 0) then
       write(output_unit,*) "found ",group_num_re," regular classes "
       write(output_unit,*) "and   ",group_num_pir,"projective irreps "
       write(output_unit,601) group_name, (symm_el(klass(regular(j))%conjug_el(1))&
            &%name,&
            &          j = 1,group_num_re)
       do j = 1,group_num_pir
          write(output_unit,602) proj_irrep_can(j)%label, (proj_irrep_can(j)%characters(k),&
               &          k = 1,group_num_re)
       end do
    endif

601 format(/1x,'Group ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,A5,2x,4(f8.3,f8.3:)/4(10x,4(f8.3,f8.3:)/))


  end subroutine group_proj_irrep_label
  !*************************************************************


  !*************************************************************
  subroutine group_proj_irrep_pre_label
    !  Purpose: labels irreps using their character table
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module

    implicit none

    !------------ Definition of local constants  -----------------

    real(kind=r8_kind),parameter  :: small = 1e-5

    !------------ Definition of local variables  -----------------

    integer(kind=i4_kind)       :: i, j, k, i_lab, dim_ind, dim_max
    ! counter
    integer(i4_kind) :: error
    !-------- Executable Code ------------------------------------

    DPRINT 'proj_label: entered'

    dim_max = 0

    ! initialize number of irreps
    group_num_pir = group_num_re

    ! allocate canonically ordered irreps
    allocate(proj_irrep_can(group_num_pir))

    ! loop over all irreps
    do i = 1, group_num_re

       ! process first character of label
       i_lab = 1

       allocate(proj_irrep_can(i)%characters(group_num_re),STAT=error)
       ASSERT(error.eq.0)
       proj_irrep_can(i)%label = 'Irr     '
       write(proj_irrep_can(i)%label(5:7),'(i3)') i
       proj_irrep_can(i)%characters = proj_characters(:,i)
       proj_irrep_can(i)%eigen_value = proj_eig_val(i)
       dim_ind = int(real(proj_characters(1,i))+small)
       proj_irrep_can(i)%dimension = dim_ind
       proj_irrep_can(i)%time_dimension = dim_ind
       DPRINT 'pre_label: main_axis=',main_axis
       DPRINT 'pre_label: size(regular_inv)=',size(regular_inv)
       DPRINT 'pre_label: regular_inv(main_axis)=',regular_inv(main_axis)
       DPRINT 'pre_label: shape(proj_characters)=',shape(proj_characters)
       if(main_axis.ne.0)then
          if (aimag(proj_characters(regular_inv(main_axis),i)).gt.small) then
             proj_irrep_can(i)%sign_of_jz = 1
             elseif(aimag(proj_characters(regular_inv(main_axis),i)).lt.-small) then
             proj_irrep_can(i)%sign_of_jz = -1
          else
             proj_irrep_can(i)%sign_of_jz = 0
          endif
       else
          WARN("sign of Jz may be determined wrong")
          if(group_name.ne."CS  ")then
             WARN("sign of Jz is probably wrong")
          endif
          if(group_name.ne."C1  ")then
          ! workaround:
          error = size(proj_characters,1) - 2
          ASSERT(error.ge.0)
          if(aimag(proj_characters(2,i)).gt.small)then
             proj_irrep_can(i)%sign_of_jz = +1
          else
             proj_irrep_can(i)%sign_of_jz = -1
          endif
          else
          WARN("jz for C1 not implemented")
          proj_irrep_can(i)%sign_of_jz = 0
          endif
       endif


       ! check for conjugated irrep pair
       proj_irrep_can(i)%pseudo = .false.
       do j=1,group_num_re
          if (abs(aimag(proj_characters(j,i))).gt.small) then
             ! case of nonambivalent group, where conjugated
             ! irrep pairs appear, which form a two dim. irrep
             ! decrement number of irreps
             proj_irrep_can(i)%pseudo = .true.
             exit
          end if
       end do

    end do

    if (output_unit > 0) then
       write(output_unit,*) ">>> Preliminary Labeling of Irreps <<<"
       write(output_unit,*) "found ",group_num_re," regular classes and Irreps"
       write(output_unit,601) group_name, (symm_el(klass(regular(j))%conjug_el(1))&
            &%name,&
            &          j = 1,group_num_re)
       do j = 1,group_num_pir
          write(output_unit,602) proj_irrep_can(j)%label, (proj_irrep_can(j)%characters(k),&
               &          k = 1,group_num_re)
       end do
    endif

601 format(/1x,'Group ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,A5,2x,4(f8.3,f8.3:)/4(10x,4(f8.3,f8.3:)/))


    DPRINT 'proj_label: exit'
  end subroutine group_proj_irrep_pre_label
  !*************************************************************



  !*************************************************************
  subroutine subgroup_irrep_label(subgroup)
    !  Purpose: determines labels for (invariant) subgroups
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    type(sub_group)                         :: subgroup
    ! contains all information about the subgroup
    !** End of interface *****************************************
    !------------ Definition of local types ----------------------
    type irrep
       integer(kind=i4_kind)  :: irrep_number
       real(kind=r8_kind)     :: main_axis_value
       logical                :: flag
       type(irrep),pointer    :: next => NULL()
    end type irrep
    ! a special irrep

    type irrep_list
       type(irrep),pointer  :: start => NULL()
       type(irrep),pointer  :: end => NULL()
    end type irrep_list
    ! list of a group of irreps of the same size
    ! and inversion symmetry

    !------------ Definition of local constants  -----------------

    real(kind=r8_kind),parameter  :: small = 1e-5

    !------------ Definition of local variables  -----------------

    integer(kind=i4_kind)       :: i,j,k,i_lab,inv_ind,dim_ind
    ! counter

    integer(kind=i4_kind)       :: inv_hor
    ! is inversion or horizontal plane present

    integer(kind=i4_kind)       :: num_cl,num_ir
    ! number of classes and irreps in subgroup

    integer(kind=i4_kind)       :: sg_primary_axis,sg_secondary_axis,sg_main_axis
    ! axes of the subgroup

    integer(kind=i4_kind)       :: sg_inversion,sg_horizontal
    ! inversion and horizontal plane of the subgroup

    logical                     :: single
    ! is .true. if there is only one irrep of a special kind

    character                  :: even_char, odd_char
    ! even and odd respective inversion or horizontal plane

    type(irrep), pointer       :: new_irrep,irr_point,previous
    ! newly generated irrep, pointer to some irrep

    type(irrep), pointer       :: irr_point2,irr_pointnew,irr_pointold
    ! newly generated irrep, pointer to some irrep

    type(irrep_list)           :: irrep_kinds(2,0:5)
    ! irrep_kinds points to irrep list of a special kind of
    ! irreps, which is specified by inversion symmetry and
    ! dimension

    !-------- Executable Code ------------------------------------

    ! set abreviations for derived types
    num_cl = subgroup%num_cl
    num_ir = subgroup%num_ir
    sg_primary_axis = subgroup%primary_axis
    sg_secondary_axis = subgroup%secondary_axis
    sg_main_axis = subgroup%main_axis
    sg_inversion = subgroup%inversion
    sg_horizontal = subgroup%horizontal

    ! allocate canonically ordered irreps
    allocate(subgroup%irrep_can(num_ir))

    ! initialize irrep_can, which contains all information
    ! about the irreps (in canonical order)
    do i=1,num_ir
       allocate(subgroup%irrep_can(i)%characters(num_ir))
       subgroup%irrep_can(i)%label = '        '
    enddo

    ! check existence of inversion and horizontal plane
    if (sg_inversion.ne.0) then
       inv_hor = sg_inversion
       even_char = 'G'
       odd_char  = 'U'
    elseif (sg_horizontal.ne.0) then
       inv_hor = sg_horizontal
       even_char = ''''
       odd_char  = '"'
    else
       inv_hor = 0
       ! these are unused in this case:
       even_char = '?'
       odd_char = '?'
    endif

    ! loop over all irreps
    do i = 1, num_ir

       ! process first character of label
       i_lab = 1

       ! generate new irrep
       allocate(new_irrep)
       new_irrep%irrep_number = i
       new_irrep%flag = .false.

       ! determine inversion symmetry of irrep
       inv_ind = 1
       if (inv_hor.ne.0) then
          if (real(subgroup%characters(inv_hor,i)).ge.0) then
             inv_ind = 1
          else
             inv_ind = 2
          end if
       end if

       ! determine dimension of irrep
       dim_ind = int(real(subgroup%characters(1,i))+small)

       ! is one dimensional irrep A or B?
       if (dim_ind.eq.1) then
          ! A irrep
          dim_ind = 0
          if (sg_primary_axis.ne.0) then
             if (abs(aimag(subgroup%characters(sg_primary_axis,i))).gt.small) then
                ! case of nonambivalent group, where conjugated
                ! irrep pairs appear, which form a two dim. irrep
                dim_ind = 2
                elseif (real(subgroup%characters(sg_primary_axis,i)).lt.0) then
                   ! B irrep
                dim_ind = 1
             endif
          endif
       endif

       ! determine by which value irreps are ordered
       if (dim_ind.gt.1) then
          if (sg_main_axis.ne.0) then
             new_irrep%main_axis_value = real(subgroup%characters(sg_main_axis,i))
          end if
       else
          if (sg_secondary_axis.ne.0) then
             new_irrep%main_axis_value = real(subgroup%characters(sg_secondary_axis,i))
          end if
       end if


       if (associated(irrep_kinds(inv_ind,dim_ind)%start)) then
          ! extend old list
          irrep_kinds(inv_ind,dim_ind)%end%next => new_irrep
          irrep_kinds(inv_ind,dim_ind)%end => new_irrep
       else
          ! make new list
          irrep_kinds(inv_ind,dim_ind)%start => new_irrep
          ! first irrep exists double, because interchange is simpler
          allocate(new_irrep)
          new_irrep = irrep_kinds(inv_ind,dim_ind)%start
          irrep_kinds(inv_ind,dim_ind)%start%main_axis_value = 20.0
          nullify(new_irrep%next)
          irrep_kinds(inv_ind,dim_ind)%start%next => new_irrep
          irrep_kinds(inv_ind,dim_ind)%end => new_irrep
       end if

    end do

    ! irrep number
    i = 1
    ! process all kinds of irreps
    do inv_ind = 1,2
       do dim_ind = 0,5
             if (.not.associated(irrep_kinds(inv_ind,dim_ind)%start)) cycle
             ! irrep_kinds(inv_ind,dim_ind)%end%flag = .true.
             ! sort irreps of one kind
             do
                irr_pointold => irrep_kinds(inv_ind,dim_ind)%start
                irr_point => irr_pointold%next
                single = .false.
                if (.not.associated(irr_point%next)) then
                   single = .true.
                   exit
                   elseif (irr_point%next%flag) then
                   exit
                end if
                do
                   if (irr_point%main_axis_value.lt.irr_point%next%&
                        &main_axis_value) then
                      irr_point2        => irr_point%next
                      irr_pointnew      => irr_point%next%next

                      irr_point%next%next =>irr_point
                      irr_point%next => irr_pointnew
                      irr_pointold%next => irr_point2

                      irr_point         => irr_point2
                   end if
                   irr_pointold => irr_point
                   irr_point => irr_point%next
                   if (.not.associated(irr_point%next)) then
                      irr_point%flag = .true.
                      exit
                      elseif (irr_point%next%flag) then
                      irr_point%flag = .true.
                      exit
                   end if
                end do
             end do
             ! number j of multi dimensional irrep
             ! e.g. E1, E2
             j=1
             previous => irrep_kinds(inv_ind,dim_ind)%start
             irr_point => irrep_kinds(inv_ind,dim_ind)%start%next
             do while ( associated(irr_point) )
                   ! write character table and eigenvalues according to
                   ! the canonical order of the irreps
                   subgroup%irrep_can(i)%characters = subgroup%characters(:,irr_point%irrep_number)
                   subgroup%irrep_can(i)%eigen_value = eig_val(irr_point%irrep_number)
                   if (dim_ind.eq.0) then
                      subgroup%irrep_can(i)%dimension = 1
                   else
                      subgroup%irrep_can(i)%dimension = dim_ind
                   endif
                   !char_tab(:,i)  = characters(:,irr_point%irrep_number)
                   !eig_val_can(i) = eig_val(irr_point%irrep_number)
                   select case(dim_ind)
                   case(0)
                      subgroup%irrep_can(i)%label(1:1) = 'A'
                   case(1)
                      subgroup%irrep_can(i)%label(1:1) = 'B'
                   case(2)
                      subgroup%irrep_can(i)%label(1:1) = 'E'
                   case(3)
                      subgroup%irrep_can(i)%label(1:1) = 'T'
                   case(4)
                      subgroup%irrep_can(i)%label(1:1) = 'G'
                   case(5)
                      subgroup%irrep_can(i)%label(1:1) = 'H'
                   end select
                   if (abs(irr_point%main_axis_value-previous%&
                        &main_axis_value).lt.small) then
                      j = j - 1 ! decrement number of multi
                   end if
                   if (.not.single) then
                      write(subgroup%irrep_can(i)%label(2:2),'(I1)') j
                   end if
                   if (inv_hor.ne.0) then
                      if (.not.single) then
                         if (real(subgroup%characters(inv_hor,irr_point%irrep_number)).ge.0)&
                              & then
                            subgroup%irrep_can(i)%label(3:3) = even_char
                         else
                            subgroup%irrep_can(i)%label(3:3) = odd_char
                         endif
                      else
                         if (real(subgroup%characters(inv_hor,irr_point%irrep_number)).ge.0)&
                              & then
                            subgroup%irrep_can(i)%label(2:2) = even_char
                         else
                            subgroup%irrep_can(i)%label(2:2) = odd_char
                         endif

                      end if
                   endif
                   j = j + 1 ! increment number of multi
                   i = i + 1 ! increment irrep number
                   previous  => irr_point
                   irr_point => irr_point%next
             end do
       end do! dimensions
    end do! inv_hor

    if (output_unit > 0) then
       write(output_unit,"(1x,'subgroup ',a,' eigenvalues of class operator ','and primitive characters :')")&
            & subgroup%name
       write(output_unit,"(5(10x,4(a8,8x :)))")&
            & (symm_el(subgroup%klass(j)%conjug_el(1))%name,j = 1,num_ir)
       do j = 1,num_ir
          write(output_unit,"(1x,A3,2x,4(f8.3,f8.3:)/4(10x,4(f8.3,f8.3:)/))")&
               &  subgroup%irrep_can(j)%label,&
               & (subgroup%irrep_can(j)%characters(k), k = 1,num_ir)
       end do
    endif
  end subroutine subgroup_irrep_label
  !*************************************************************


end module group_module

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
module datatype

!  These are the user-defined data types

  use type_module

  !  Array of real matrices of variable size
  type arrmat1
     real(kind=r8_kind), allocatable :: m(:)
  endtype arrmat1

  !  Array of complex matrices of variable size
  type arrmat1_c
     complex(kind=c16_kind), pointer :: m(:) => NULL()
  endtype arrmat1_c

  !  Array of real matrices of variable size
  type arrmat2
     real(kind=r8_kind), allocatable :: m(:, :)
  endtype arrmat2

  !  Array of complex matrices of variable size
  type arrmat2_c
     complex(kind=c16_kind), pointer :: m(:,:) => NULL()
  endtype arrmat2_c

  !  Array of integer matrices of variable size
  type arrmat1int
     integer(kind=i4_kind), allocatable :: m(:)
  endtype arrmat1int

  !  Array of integer matrices of variable size
  type arrmat2int
     integer(kind=i4_kind), allocatable :: m(:, :)
  endtype arrmat2int

  !  Array of real 3D-matrices of variable size
  type arrmat3
     real(kind=r8_kind), allocatable :: m(:,:,:)
  endtype arrmat3

  !  Array of real 3D-matrices of variable size
  type arrmat4
     real(kind=r8_kind), allocatable :: m(:,:,:,:)
  endtype arrmat4

  ! Array of real 5D-matrices of variables size
  type arrmat5
     real(kind=r8_kind), allocatable :: m(:,:,:,:,:)
  endtype arrmat5

  type arrmat6
     real(kind=r8_kind), allocatable :: m(:,:,:,:,:,:)
  endtype arrmat6

  type arrmat7
     real(kind=r8_kind), allocatable :: m(:,:,:,:,:,:,:)
  endtype arrmat7

  ! Array of logical 3D-matrices of variable size
  type arrmat3log
     logical,pointer :: m(:,:,:) => NULL()
  end type arrmat3log

  type three_center_l
     type(arrmat5), allocatable :: l(:)
  end type three_center_l

  ! SB: second version for Coulomb
  type three_center_l_v2
     type(arrmat6), allocatable :: l(:)
  end type three_center_l_v2

  type intmat1
     integer(kind=i4_kind), pointer :: m(:) => NULL()
  end type intmat1

  type intmat2
     integer(kind=i4_kind), pointer :: m(:,:) => NULL()
  end type intmat2

  type fragment_type ! used to describe fragments
     !in fragment orbital analysis
     integer(kind=i4_kind) :: n_atoms ! number of atoms in this fragment
     integer(kind=i4_kind), pointer :: atoms(:) => NULL() ! indices of
     ! atoms in the fragment
  end type fragment_type

  type eig_cont_type
     ! this type is used to store results of fragment orbital analysis
     integer(kind=i4_kind) :: n_cont ! number of cont. orbitals
     integer(kind=i4_kind), pointer :: index(:) => NULL() ! indices of cont orbitals
     real(kind=r8_kind), pointer    :: pop(:) => NULL() ! population
     real(kind=r8_kind), pointer    :: coeff(:) => NULL() ! orbital coefficient
     real(kind=r8_kind)             :: sum ! sum of all contributions of one fragment
  end type eig_cont_type

  type pop_store_type
     ! used to get a array of eig_cont_type
     type(eig_cont_type), pointer :: eig_cont(:) => NULL()
  end type pop_store_type

type, public ::  pointcharge_type
   character(len=12)                 :: name
       ! name of unique point charge
   real(kind=r8_kind)                :: Z
       ! charge of point charge
   integer                           :: N_equal_charges
       ! Number of partners of unique point charge
   real(kind=r8_kind), pointer       :: position(:,:) => NULL() ! position(3,N_equal_charges)
       ! positions of partners ("equal charges") of unique charge
   real(kind=r8_kind)                :: position_first_ec(3)
       ! positions of first partner ("equal charge") of unique charge
   real(kind=r8_kind)                :: c
   real(kind=r8_kind)                :: a
   real(kind=r8_kind)                :: cf
   real(kind=r8_kind)                :: af
       ! screening parameters (EFP model)
   integer(i4_kind), pointer :: group(:) => NULL()
end type pointcharge_type

end module datatype
!================================================================
! End of public interface of module
!================================================================

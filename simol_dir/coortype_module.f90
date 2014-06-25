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
module  coortype_module
  !---------------------------------------------------------------
  !
  !  Purpose: Contains the data types for coordinates:
  !           -bond_length, bond_angle, dihedral_angle
  !           -int_coor, atom_type
  !            Furthermore the parameters
  !           -b_length, b_angle, d_angle
  !  Author: FN
  !  Date: 2/98
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------

  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------
  type,public ::  bond_length
     integer(kind=i4_kind)  :: partner1,partner2
  end type bond_length
  type,public ::  bond_angle
     integer(kind=i4_kind)   :: partner1,partner2
     integer(kind=i4_kind)   :: apex
  end type bond_angle
  type, public ::  dihedral_angle
     integer(kind=i4_kind)   :: partner1,partner2
     integer(kind=i4_kind)   :: base1,base2
  end type dihedral_angle
  type, public ::  int_coor
     integer(kind=i4_kind):: typ  ! can be b_length,b_angle or d_angle
     type(bond_length)    :: length 
     type(bond_angle)     :: angle
     type(dihedral_angle) :: dihedral
     real(kind=r8_kind)   :: value
     logical              :: unique   ! serves to distinguish the internal coord.
     !                                ! from the NON symmetry equivalent coords.
     logical              :: var      ! if TRUE : this component is to be varied
     !                                ! if FALSE: this component is to be klept constant
     !                                !           (contrained optimization)
     integer(kind=i4_kind),pointer :: equal(:) ! contains the symmetry-equivalent 
     !                                         ! atoms if unique=.true.
     integer(kind=i4_kind)         :: n_equal
     integer(kind=i4_kind),pointer :: equal_sign(:) ! in case the symmetry-equivalent
     !                                ! variable has the same magnitude but opposite sign
     !                                ! with reference to the 'unique' coordinate, this
     !                                ! array contains a '-'. Otherwise it will be "+"
  end type int_coor
  type, public ::  atom_type
     real(kind=r8_kind)   :: x(3)
     real(kind=r8_kind)   :: x_old(3)
     real(kind=r8_kind)   :: charge
     real(kind=r8_kind)   :: mass
     logical              :: dummy ! true if atom is a dummy
     logical              :: heavy ! true if atom has infinite mass,
     ! e. g. is kept fixed. Only used for frequnecy calculation
  end type atom_type

!------------ Declaration of constants and variables ------------
  integer(kind=i4_kind),parameter,public :: b_length=1, &
                                            b_angle=2,&
                                            d_angle=3
                      
  character(len=20),public               :: bond_type(3)
  data bond_type /'bond_length','bond_angle','dihedral_angle'/

!================================================================
! End of public interface of module
!================================================================


!--------------- End of module ----------------------------------
end module coortype_module

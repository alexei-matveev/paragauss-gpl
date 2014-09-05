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
subroutine ss_calculate_dipole()
  !
  !  Purpose: calculation of primitive dipole integrals
  !       for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,l2,equal_atom1,equal_atom2)
  !       with l1 = l2 = 0.
  !
  !  Author: TB
  !  Date:   9/97
  !
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  use type_module
  use unique_atom_module
  use integralpar_module
  use int_data_dipole_module
  use options_module, only: options_integral_expmax

  implicit none

  integer(kind=r4_kind) :: naexps,nbexps
  real(kind=r8_kind),pointer :: aexps(:),bexps(:)

  ! constants
  real(kind=r8_kind),parameter    :: two=2.0_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind

  ! mapping of exponents to one dimension and cutoff of small integrals
  logical,allocatable   :: cutoff(:,:) ! (naexps,nbexps)
  integer(kind=i4_kind) :: num ! metaindex for (naexps,nbexps) > cutoff

  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: &
       fact0_arr, fact1_arr, fact2_arr ! (naexps,nbexps)
  real(kind=r8_kind),allocatable,dimension(:)  :: &
       fact0, fact1, fact2 ! (num) metaindex for (naexps,nbexps) > cutoff

  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:)  :: gamma_arg ! (num,3)

  ! help variables
  real(kind=r8_kind) :: arg
  real(kind=r8_kind),dimension(3)  :: xa,xb
  integer(kind=i4_kind)  :: i_xyz, alloc_stat

  ! integrals
  real(kind=r8_kind),allocatable ::  &
       overlap(:), & ! (num)
       dipole(:,:) ! (num,3)



  naexps = ua1_basis%n_exponents
  nbexps = ua2_basis%n_exponents
  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("ss_calculate_dipole : allocation (1) failed")

  xa = center1
  xb = center2

  aexps => ua1_basis%exponents(:)
  bexps => ua2_basis%exponents(:)


  arg=sum((xa-xb)**2)

  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  where(fact2_arr*arg>options_integral_expmax()) ! cutoff: where almost no overlap
     cutoff=.false.              ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)

  if(num==0) then ! all integrals are equal zero
     prim_int_2cob_dipole = 0.0_r8_kind
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("ss_calculate_dipole: deallocation (1/1) failed")
     return
  end if

  allocate ( &
       fact0(num), &
       fact1(num), &
       fact2(num), &
       gamma_arg(num,3), &
       overlap(num), &
       dipole(num,3), &
       STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ss_calculate_dipole: allocation (2) failed")

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)

  deallocate( &
       fact0_arr, &
       fact1_arr, &
       fact2_arr, &
       STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ss_calculate_dipole: deallocation (1/2a) failed")

  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  do i_xyz = 1, 3
     gamma_arg(:,i_xyz)=(pack(spread(aexps*xa(i_xyz),1,nbexps) + &
          spread(bexps*xb(i_xyz),2,naexps),cutoff))/fact0
  enddo

  ! overlap ----------------------
  overlap = (two*sqrt(fact1)/fact0)* &
       sqrt((two*sqrt(fact1)/fact0))*exp(-fact2*arg)

  ! dipole ----------------------
  do i_xyz = 1, 3
     dipole(:,i_xyz) = overlap * gamma_arg(:,i_xyz)
  enddo

  ! re-map them to the int_data_2cob3c_stuff
  do i_xyz = 1, 3
     prim_int_2cob_dipole(:,:,1,1,i_xyz) = unpack(dipole(:,i_xyz),cutoff,zero)
  enddo


  deallocate ( &
       fact0, &
       fact1, &
       fact2, &
       gamma_arg, &
       overlap, &
       dipole, &
       cutoff, &
       STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ss_calculate_dipole: deallocation (1/2b) (2) failed")

end subroutine ss_calculate_dipole

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
subroutine ll_calculate_dipole(na,nb,la,lb)
  !
  !  Purpose: calculation of all primitive dipole integrals
  !       for a given set of indizes
  !       (unique_atom1,unique_atom2,la,lb).
  !
  !  Author: MS
  !  Date:   8/96
  !
  use unique_atom_module
  use type_module
  use solid_harmonics_module, only : solid_harmonics_calc,solid_harmonics_scalar
  use int_data_dipole_module  
  use solhrules_module
  use integralpar_module
  use options_module, only: options_integral_expmax
  implicit none

  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
  integer(kind=i4_kind),intent(in) :: la ! angular momentum of unique atom a
  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
  !===================================================================
  ! End of public interface of module
  !===================================================================

  ! constants
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: gam=1.0_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
  real(kind=r8_kind),parameter,dimension(0:8) :: dfac= (/ &
       1.0_r8_kind, 1.0_r8_kind, 3.0_r8_kind, 15.0_r8_kind, 105.0_r8_kind, &
       945.0_r8_kind, 10395.0_r8_kind, 135135.0_r8_kind, 2027025.0_r8_kind /)

  ! variables
  integer(kind=i4_kind) :: num,m,alloc_stat,i_xyz,i_l,i_lma,i_lmb, &
       n_lma,n_lmb,naexps,nbexps,ma,mb,l
  logical,allocatable   :: cutoff(:,:)
  real(kind=r8_kind),pointer,dimension(:) :: &
       aexps,bexps
  real(kind=r8_kind),allocatable,dimension(:) :: &
       fact0,fact1,fact2,fact4,fact6,tau,aexp_arr,bexp_arr,clmamb_scalar
  real(kind=r8_kind),allocatable,dimension(:,:) :: &
       fact0_arr,fact1_arr,fact2_arr,clmamb,diff_arr0
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: &
       overlap
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: &
       dipole,diff_rule_result
  real(kind=r8_kind),dimension(3)  :: xa,xb,xd
  real(kind=r8_kind) :: arg



  naexps = unique_atoms(na)%l_ob(la)%n_exponents
  nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents


  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ll_calculate_dipole: allocation (1) failed") 


  xa = center1
  xb = center2
  xd =xa-xb
  aexps => unique_atoms(na)%l_ob(la)%exponents(:)
  bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)

  n_lma = (la+1)**2
  n_lmb = (lb+1)**2

  arg=sum(xd**2)

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
          ("ll_calculate_dipole: deallocation (1/1) failed")
     return
  end if

  allocate ( &
       fact0(num), &
       fact1(num), &
       fact2(num), &
       fact4(num), &
       fact6(num), &
       tau(num), &
       overlap(num,n_lma,n_lmb), &
       dipole(num,2*la+1,2*lb+1,3), &
       diff_rule_result(num,n_lma,n_lmb,3), &
       aexp_arr(num), &
       bexp_arr(num), &
       clmamb_scalar((max(la,lb)+1)**2), &
       clmamb(num,(la+1)**2), &
       diff_arr0((la+1)**2,(lb+1)**2), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls_calculate_dipole: allocation (2) failed")


  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)

  aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  bexp_arr=pack(spread(bexps,2,naexps),cutoff)


  deallocate( &
       fact0_arr, &
       fact1_arr, &
       fact2_arr, &
       stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("ll_calculate_dipole: deallocation  (1/2a) failed")


  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(max(la,lb),xd)
  fact4=1.0_r8_kind
  i_lma=1
  tau=fact2*arg ! a*b/(a+b)*(A-B)**2
  do l=0,la
     do m=1,2*l+1
        clmamb(:,i_lma)=clmamb_scalar(i_lma)*fact4
        i_lma=i_lma+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo



  ! calculate overlap 
  fact6= exp(-tau) * (4.0_r8_kind*fact2/fact0)**0.75_r8_kind  &
       / ( sqrt(aexp_arr**la*dfac(la)) * sqrt(bexp_arr**lb*dfac(lb)) )
  i_lmb=1
  do i_l=0,lb
     do mb=1,2*i_l+1
        diff_arr0(:,i_lmb)= &
             reshape( &
             diff_rule( spread(clmamb_scalar,1,1), 1, n_lma, i_lmb ), &
             (/n_lma/) &
             )
        i_lmb=i_lmb+1
     end do
  end do
  i_lmb=1
  do i_l=0,lb
     do mb=1,2*i_l+1
        overlap(:,:,i_lmb)= &
             spread( fact6*(2.0_r8_kind*fact2)**i_l, 2, n_lma ) * &
             prod_rule( &
             spread( diff_arr0(:,i_lmb), 1, num ), &
             clmamb(:,:), 1, n_lma &
             )
        i_lmb = i_lmb+1
     end do
  enddo


  ! calculate double differential rule on C(1,m) and multiply all scaling factors
  diff_rule_result = 0.0_r8_kind
  diff_rule_result(:,1,1,1) = &
       pack( &
       spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps), &
       cutoff &
       ) / fact0
  diff_rule_result(:,1,1,2) = &
       pack( &
       spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps), &
       cutoff &
       ) / fact0
  diff_rule_result(:,1,1,3) = &
       pack( &
       spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps), &
       cutoff &
       ) / fact0
  fact6 = aexp_arr / fact0
  do i_lma = 2, 4
     diff_rule_result(:,i_lma,1,i_lma-1) = fact6
  enddo
  fact6 = bexp_arr / fact0
  do i_lmb = 2, 4
     diff_rule_result(:,1,i_lmb,i_lmb-1) = fact6
  enddo


  ! double product rule with respect to a and b
  do i_xyz = 1, 3
     dipole(:,:,:,i_xyz) = &
          prod_rule_double( &
          diff_rule_result(:,:,:,i_xyz), &
          overlap, &
          la**2+1,la**2+2*la+1, &
          lb**2+1,lb**2+2*lb+1 &
          )
  enddo


  ! re-map to int_data_2cob_dipole
  ! Take into account special mapping of i_xyz to x,y,z
  ! determined by solidharmonics for l=1 and m=1,2,3
  do ma=1,2*la+1
     do mb=1,2*lb+1
        prim_int_2cob_dipole(:,:,mb,ma,1) = &
             unpack(dipole(:,ma,mb,2),cutoff,zero)
        prim_int_2cob_dipole(:,:,mb,ma,2) = &
             unpack(dipole(:,ma,mb,3),cutoff,zero)
        prim_int_2cob_dipole(:,:,mb,ma,3) = &
             unpack(dipole(:,ma,mb,1),cutoff,zero)
     end do
  end do


  deallocate ( &
       fact0, &
       fact1, &
       fact2, &
       fact4, &
       fact6, &
       tau, &
       overlap, &
       dipole, &
       diff_rule_result, &
       aexp_arr, &
       bexp_arr, &
       clmamb_scalar, &
       clmamb, &
       diff_arr0, &
       cutoff, &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls_calculate_dipole: deallocation (1/2b) (2) failed")

end subroutine ll_calculate_dipole

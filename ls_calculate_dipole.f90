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
subroutine ls_calculate_dipole(na,nb,la_in,lb)
  !
  !  Purpose: calculation of all primitive dipole integrals
  !           for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,equal_atom2).
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
  integer(kind=i4_kind),intent(in) :: la_in ! angular momentum of unique atom a
  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b

!================================================================
! End of public interface of module
!================================================================

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
  integer(kind=i4_kind) :: num,m,alloc_stat,i_xyz,i_l,i_m,i_lm,n_lm, &
       nm_la,nm_lb,la_org,lb_org,la,naexps,nbexps
  logical,allocatable   :: cutoff(:,:)
  real(kind=r8_kind),pointer,dimension(:) :: &
       aexps,bexps
  real(kind=r8_kind),allocatable,dimension(:,:) :: &
       fact0_arr,fact1_arr,fact2_arr,overlap,solh,gamma_arg
  real(kind=r8_kind),allocatable,dimension(:) :: &
       fact0,fact1,fact2,fact4,fact6,tau,aexp_arr,clmamb_scalar
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: &
       dipole,diff_rule_result
  real(kind=r8_kind),dimension(3)  :: xa,xb,xd
  real(kind=r8_kind) :: arg
  logical :: laltlb  ! flag to decide if la is lower then lb or not




  nm_la=2*la_in+1
  nm_lb=2*lb+1
  la_org=la_in
  lb_org=lb

  naexps = unique_atoms(na)%l_ob(la_in)%n_exponents
  nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents

  if(lb>la_in) then
     laltlb=.true.
     la=lb
  else
     laltlb=.false.
     la=la_in
  end if

  n_lm = (la+1)**2

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls_calculate_dipole: allocation (1) failed")

  xa = center1  ! from int_data_module
  xb = center2  ! from int_data_module
  xd =xa-xb

  aexps => unique_atoms(na)%l_ob(la_in)%exponents(:)
  bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)

  arg=sum(xd**2)

  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  ! cutoff: where almost no overlap is present calculation is not necessary
  where(fact2_arr*arg>=options_integral_expmax())
     cutoff=.false.
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)


  if(num==0) then    ! all integrals are equal zero
     prim_int_2cob_dipole = 0.0_r8_kind
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("ls_calculate_dipole: deallocation (1/1) failed")
     return
  end if

  allocate ( &
       fact0(num), &
       fact1(num), &
       fact2(num), &
       fact4(num), &
       fact6(num), &
       tau(num), &
       overlap(num,n_lm), &
       solh(num,n_lm), &
       dipole(num,2*la+1,3), &
       diff_rule_result(num,n_lm,3), &
       gamma_arg(num,3), &
       aexp_arr(num), &
       clmamb_scalar((la+1)**2),&
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


  if(.not.laltlb) then
     aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  else
     aexp_arr=pack(spread(bexps,2,naexps),cutoff)
  end if

  deallocate( &
       fact0_arr, &
       fact1_arr, &
       fact2_arr, &
       stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("ls_calculate_dipole: deallocation  (1/2a) failed")


  ! precalculation of solid harmonics
  if(laltlb) then
     xd=-xd
  end if
  clmamb_scalar=solid_harmonics_scalar(la,xd)


  ! calculate overlap * sqrt(aexp_arr**i_l*dfac(i_l))
  tau = fact2 * arg
  fact6 = exp(-tau) * (4.0_r8_kind*fact2/fact0)**0.75_r8_kind
  fact4=1.0_r8_kind
  i_lm=1
  do i_l=0,la
     do i_m=1,2*i_l+1
        overlap(:,i_lm) = fact6 * fact4 * clmamb_scalar(i_lm)
        i_lm=i_lm+1
     enddo
     fact4 = - fact4 * fact2 * 2.0_r8_kind
  enddo


  ! calculate differential rule on C(1,m) and multiply all scaling factors
  ! gamma_arg = (a*vec_a + b*vec_b) / a
  do i_xyz = 1, 3
     gamma_arg(:,i_xyz)= &
          pack( &
          spread( aexps*xa(i_xyz), 1, nbexps ) + &
          spread( bexps*xb(i_xyz), 2, naexps ), &
          cutoff &
          ) / aexp_arr
  enddo
  solh = solid_harmonics_calc(la,gamma_arg)
  ! multiply scaling factor already here to save work:
  ! solh = solh * a / (a+b) / sqrt(a**la * dfac(la))
  fact6 = aexp_arr / ( fact0 * sqrt(aexp_arr**la * dfac(la)) )
  do i_lm = 1, n_lm
     solh(:,i_lm) = &
          solh(:,i_lm) * fact6
  enddo
  do i_xyz = 1, 3
     diff_rule_result(:,:,i_xyz) = diff_rule( solh, 1, n_lm, i_xyz+1 )
  enddo


  ! product rule with respect to a
  do i_xyz = 1, 3
     dipole(:,:,i_xyz) = &
          prod_rule( &
          diff_rule_result(:,:,i_xyz), &
          overlap, &
          la**2+1,la**2+2*la+1 &
          )
  enddo



  ! re-map to int_data_2cob_dipole
  ! Take into account special mapping of i_xyz to x,y,z
  ! determined by solidharmonics for l=1 and m=1,2,3
  if(.not.laltlb) then
     magnetic_number_map: do m=1,2*la+1
        prim_int_2cob_dipole(:,:,1,m,1) = &
             unpack(dipole(:,m,2),cutoff,zero)
        prim_int_2cob_dipole(:,:,1,m,2) = &
             unpack(dipole(:,m,3),cutoff,zero)
        prim_int_2cob_dipole(:,:,1,m,3) = &
             unpack(dipole(:,m,1),cutoff,zero)
     end do magnetic_number_map
  else
     magnetic_number_map_2: do m=1,2*la+1
        prim_int_2cob_dipole(:,:,m,1,1) = &
             unpack(dipole(:,m,2),cutoff,zero)
        prim_int_2cob_dipole(:,:,m,1,2) = &
             unpack(dipole(:,m,3),cutoff,zero)
        prim_int_2cob_dipole(:,:,m,1,3) = &
             unpack(dipole(:,m,1),cutoff,zero)
     end do magnetic_number_map_2
  end if


  deallocate ( &
       fact0, &
       fact1, &
       fact2, &
       fact4, &
       fact6, &
       tau, &
       overlap, &
       solh, &
       dipole, &
       diff_rule_result, &
       gamma_arg, &
       aexp_arr, &
       clmamb_scalar,&
       cutoff, &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls_calculate_dipole: deallocation (1/2b) (2) failed")

end subroutine ls_calculate_dipole

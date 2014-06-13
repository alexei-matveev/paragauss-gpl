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
module dip_prim_module

  use type_module, only: &
       IK=>i4_kind, &
       RK=>r8_kind
  use unique_atom_module, only: &
       unique_atom_basis_type
  use solid_harmonics_module, only: &
       solid_harmonics_calc, &
       solid_harmonics_scalar
  use solhrules_module
  use options_module, only:&
       options_integral_expmax

private

public  ll_ovrl, ll_dip, ll_momnt

logical,parameter :: debug=.false.


contains

!!$function options_integral_expmax()
!!$  real(RK) :: options_integral_expmax
!!$  ! not to use options_module, only: options_integral_expmax
!!$  options_integral_expmax=50.0_RK
!!$end function options_integral_expmax

!************************************************
!***         XY_CALCULATE_ WHATEVER           ***
!************************************************
subroutine ll_ovrl(bexps, aexps, lb, la, center2, center1, ovrl)
  !
  !  Purpose: wrapper for three ss_, ls_ and ll_ subroutines
  !           
  implicit none
  intent(in)  aexps, bexps, la, lb, center1, center2
  intent(out) ovrl
  real(RK), dimension(:)        :: aexps, bexps
  integer(IK)                   :: la ! angular momentum of unique atom a
  integer(IK)                   :: lb ! angular momentum of unique atom b
  real(RK),dimension(3)         :: center1,center2
  real(RK),dimension(:,:,:,:)   :: ovrl

  if(la==0.and.lb==0)then
     call ss(bexps, aexps, center2, center1, OVRL_INT=ovrl)
     return
  endif
  if(la>0.and.lb>0)then
     call ll(bexps, aexps, lb, la, center2, center1, OVRL_INT=ovrl)
     return
  endif
     call ls(bexps, aexps, lb, la, center2, center1, OVRL_INT=ovrl)
     return
end subroutine ll_ovrl

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ll_dip(bexps, aexps, lb, la, center2, center1, dip)
  !
  !  Purpose: wrapper for three ss_, ls_ and ll_ subroutines
  !           
  implicit none
  intent(in)  aexps, bexps, la, lb, center1, center2
  intent(out) dip
  real(RK),  dimension(:)        :: aexps, bexps
  integer(IK)                   :: la ! angular momentum of unique atom a
  integer(IK)                   :: lb ! angular momentum of unique atom b
  real(RK),dimension(3)         :: center1,center2
  real(RK),dimension(:,:,:,:,:) :: dip

  if(la==0.and.lb==0)then
     call ss(bexps, aexps, center2, center1, DIP_INT=dip)
     return
  endif
  if(la>0.and.lb>0)then
     call ll(bexps, aexps, lb, la, center2, center1, DIP_INT=dip)
     return
  endif
     call ls(bexps, aexps, lb, la, center2, center1, DIP_INT=dip)
     return
end subroutine ll_dip

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ll_momnt(bexps, aexps, lb, la, center2, center1, momnt)
  !
  !  Purpose: wrapper for three ss_, ls_ and ll_ subroutines
  !           
  implicit none
  intent(in)  aexps, bexps,la,lb,center1,center2
  intent(out) momnt
  real(RK), dimension(:)        :: aexps, bexps
  integer(IK)                   :: la ! angular momentum of unique atom a
  integer(IK)                   :: lb ! angular momentum of unique atom b
  real(RK),dimension(3)         :: center1,center2
  real(RK),dimension(:,:,:,:,:) :: momnt

  if(la==0.and.lb==0)then
     call ss(bexps, aexps, center2, center1, MOMNT_INT=momnt)
     return
  endif
  if(la>0.and.lb>0)then
     call ll(bexps, aexps, lb, la, center2, center1, MOMNT_INT=momnt)
     return
  endif
     call ls(bexps, aexps, lb, la, center2, center1, MOMNT_INT=momnt)
     return
end subroutine ll_momnt

!************************************************
!***         SS_CALCULATE_(DIPOLE)            ***
!************************************************
subroutine ss(bexps, aexps, center2, center1, &
                               momnt_int,&
                               dip_int, &
                               ovrl_int)
  !
  !  Purpose: calculation of primitive dipole integrals
  !       for a given set of indizes
  !       with l1 = l2 = 0.

  implicit none
  intent(in)  aexps, bexps,&
              center1,center2
  intent(out) momnt_int,dip_int,ovrl_int
  optional    momnt_int,dip_int,ovrl_int
  real(RK), dimension(:)         :: aexps, bexps
  real(RK), dimension(3)         :: center1,center2
  real(RK), dimension(:,:,:,:,:) :: momnt_int,dip_int
  real(RK), dimension(:,:,:,:)   :: ovrl_int
  ! *** end of interface **

  integer(IK)      :: naexps,nbexps

  ! constants
  real(RK),parameter    :: two=2.0_RK
  real(RK),parameter    :: very_small=1.0e-100_RK
  real(RK),parameter    :: very_big=1.0e100_RK
  real(RK),parameter    :: zero=0.0_RK

  ! mapping of exponents to one dimension and cutoff of small integrals
  logical,allocatable   :: cutoff(:,:) ! (naexps,nbexps)
  integer(IK) :: num ! metaindex for (naexps,nbexps) > cutoff

  ! help factors
  real(RK),allocatable,dimension(:,:):: &
       fact0_arr, fact1_arr, fact2_arr ! (naexps,nbexps)
  real(RK),allocatable,dimension(:)  :: &
       fact0, fact1, fact2 ! (num) metaindex for (naexps,nbexps) > cutoff

  ! help arrays for gamma-function
  real(RK),allocatable,dimension(:,:)  :: gamma_arg ! (num,3)

  ! help variables
  real(RK) :: arg
  real(RK),dimension(3)  :: xa,xb
  integer(IK)  :: i_xyz, alloc_stat

  ! integrals
  real(RK),allocatable ::  &
       overlap(:), & ! (num)
       dipole(:,:) ! (num,3)
  intrinsic maxval



  if(debug)print *, 'dip_prim/ss: entered'

  naexps = size(aexps)
  nbexps = size(bexps)

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("ss : allocation (1) failed")

  xa = center1
  xb = center2

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
     if(present(momnt_int))momnt_int = 0.0_RK
     if(present(dip_int))  dip_int   = 0.0_RK
     if(present(ovrl_int)) ovrl_int  = 0.0_RK
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("ss: deallocation (1/1) failed")
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
       ("ss: allocation (2) failed")

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
       ("ss: deallocation (1/2a) failed")

  ! overlap ----------------------
  overlap = (two*sqrt(fact1)/fact0)* &
       sqrt((two*sqrt(fact1)/fact0))*exp(-fact2*arg)

  if(present(ovrl_int))then
     ovrl_int(:,:,1,1) = unpack(overlap,cutoff,zero)
  endif

  if(present(dip_int))then
     ! dipole ----------------------
     ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
     do i_xyz = 1, 3
        gamma_arg(:,i_xyz)=(pack(spread(aexps*xa(i_xyz),1,nbexps) + &
             spread(bexps*xb(i_xyz),2,naexps),cutoff))/fact0

        dipole(:,i_xyz) = overlap * gamma_arg(:,i_xyz)

        ! re-map them to the int_data_2cob3c_stuff
        dip_int(:,:,1,1,i_xyz) = unpack(dipole(:,i_xyz),cutoff,zero)
     enddo
  endif

  if(present(momnt_int))then
     ! analog to dipole:
     ! <a|d/dx|b>:
     ! gamma_arg= b*{(a*vec_a + b*vec_b)/(a + b) - vec_b} ==
     ! == -2ab/(a + b)*(vec_a - vec_b)
     do i_xyz = 1, 3
        gamma_arg(1:num,i_xyz)=fact2*(xa(i_xyz)-xb(i_xyz))

        dipole(:,i_xyz) = -two * overlap * gamma_arg(:,i_xyz)

        ! re-map them to the int_data_2cob3c_stuff
        momnt_int(:,:,1,1,i_xyz) = unpack(dipole(:,i_xyz),cutoff,zero)
     enddo
     if(debug)print *, 'dip_prim/ss: abs momnt=',sum(momnt_int**2)
  endif

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
       ("ss: deallocation (1/2b) (2) failed")

end subroutine ss

!***********************************************
!****        LS_CALCULATE_(DIPOLE)          ****
!***********************************************

subroutine ls(bexps, aexps, lb, la_in, center2, center1, &
                               momnt_int, &
                               dip_int, &
                               ovrl_int)
  !
  !  Purpose: calculation of all primitive dipole integrals
  !           for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,equal_atom2).
  implicit none

  intent(in)  aexps, bexps,la_in,lb,center1,center2
  intent(out) momnt_int,dip_int,ovrl_int
  optional    momnt_int,dip_int,ovrl_int
  real(RK), dimension(:)        :: aexps, bexps
  integer(IK)                   :: la_in ! angular momentum of unique atom a
  integer(IK)                   :: lb    ! angular momentum of unique atom b
  real(RK),dimension(3)         :: center1,center2
  real(RK),dimension(:,:,:,:,:) :: momnt_int,dip_int
  real(RK),dimension(:,:,:,:)   :: ovrl_int
  ! *** end of interface ***

  ! constants
  real(RK),parameter    :: two=2.0_RK
  real(RK),parameter    :: pi=3.14159265358979324_RK
  real(RK),parameter    :: gam=1.0_RK
  real(RK),parameter    :: very_small=1.0e-100_RK
  real(RK),parameter    :: very_big=1.0e100_RK
  real(RK),parameter    :: zero=0.0_RK
  real(RK),parameter,dimension(0:8) :: dfac= (/ &
       1.0_RK, 1.0_RK, 3.0_RK, 15.0_RK, 105.0_RK, &
       945.0_RK, 10395.0_RK, 135135.0_RK, 2027025.0_RK /)

  ! variables
  integer(IK) :: num,m,alloc_stat,i_xyz,i_l,i_m,i_lm,n_lm, &
       nm_la,nm_lb,la_org,lb_org,la,naexps,nbexps

  logical,allocatable   :: cutoff(:,:)
  real(RK),allocatable,dimension(:,:) :: &
       fact0_arr,fact1_arr,fact2_arr,overlap,solh,gamma_arg
  real(RK),allocatable,dimension(:) :: &
       fact0,fact1,fact2,fact4,fact6,tau,aexp_arr,bexp_arr,clmamb_scalar
  real(RK),allocatable,dimension(:,:,:) :: &
       dipole,diff_rule_result
  real(RK),dimension(3)  :: xa,xb,xd
  real(RK) :: arg
  logical :: laltlb  ! flag to decide if la is lower then lb or not

  if(debug)print *, 'dip_prim/ls: entered'

  if(lb==0.and.la_in==0)call error_handler&
       & ("ls: calculate ss integrals with ss_...")

  nm_la=2*la_in+1
  nm_lb=2*lb+1
  la_org=la_in
  lb_org=lb

  naexps = size(aexps)
  nbexps = size(bexps)

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
       ("ls: allocation (1) failed")

  xa = center1  ! from int_data_module
  xb = center2  ! from int_data_module
  xd =xa-xb

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
     if(present(momnt_int)) momnt_int = 0.0_RK
     if(present(dip_int))   dip_int   = 0.0_RK
     if(present(ovrl_int))  ovrl_int  = 0.0_RK
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("ls: deallocation (1/1) failed")
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
       solh(num,max(n_lm,4)), &
       dipole(num,2*la+1,3), &
       diff_rule_result(num,n_lm,3), &
       gamma_arg(num,3), &
       aexp_arr(num), &
       bexp_arr(num),&
       clmamb_scalar((la+1)**2),&
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls: allocation (2) failed")


  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)

  if(.not.laltlb) then
     aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
     bexp_arr=pack(spread(bexps,2,naexps),cutoff)
  else
     aexp_arr=pack(spread(bexps,2,naexps),cutoff)
     bexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  end if

  deallocate( &
       fact0_arr, &
       fact1_arr, &
       fact2_arr, &
       stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("ls: deallocation  (1/2a) failed")

  ! precalculation of solid harmonics
  if(laltlb) then
     xd=-xd
  end if
  clmamb_scalar=solid_harmonics_scalar(la,xd)

  ! calculate overlap * sqrt(aexp_arr**i_l*dfac(i_l))
  tau = fact2 * arg
  fact6 = exp(-tau) * (4.0_RK*fact2/fact0)**0.75_RK
  fact4=1.0_RK
  i_lm=1
  do i_l=0,la
     do i_m=1,2*i_l+1
        overlap(:,i_lm) = fact6 * fact4 * clmamb_scalar(i_lm)
        i_lm=i_lm+1
     enddo
     fact4 = - fact4 * fact2 * 2.0_RK
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

  solh = solid_harmonics_calc(max(la,1),gamma_arg)
  ! multiply scaling factor already here to save work:
  ! solh = solh * a / (a+b) / sqrt(a**la * dfac(la))
  fact6 = aexp_arr / ( fact0 * sqrt(aexp_arr**la * dfac(la)) )
  do i_lm = 1, max(n_lm,4)
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

  if(present(ovrl_int).or.present(momnt_int))then
     overlap(:,la*la+1:(la+1)**2)=&
          & overlap(:,la*la+1:(la+1)**2)&
          & /sqrt(spread(aexp_arr**la,2,2*la+1)*dfac(la))
  endif

  ! re-map to int_data_2cob_dipole
  ! Take into account special mapping of i_xyz to x,y,z
  ! determined by solidharmonics for l=1 and m=1,2,3
  if(present(dip_int))then
     if(.not.laltlb) then
        magnetic_number_map: do m=1,2*la+1
           dip_int(:,:,1,m,1) = &
                unpack(dipole(:,m,2),cutoff,zero)
           dip_int(:,:,1,m,2) = &
                unpack(dipole(:,m,3),cutoff,zero)
           dip_int(:,:,1,m,3) = &
                unpack(dipole(:,m,1),cutoff,zero)
        end do magnetic_number_map
     else
        magnetic_number_map_2: do m=1,2*la+1
           dip_int(:,:,m,1,1) = &
                unpack(dipole(:,m,2),cutoff,zero)
           dip_int(:,:,m,1,2) = &
                unpack(dipole(:,m,3),cutoff,zero)
           dip_int(:,:,m,1,3) = &
                unpack(dipole(:,m,1),cutoff,zero)
        end do magnetic_number_map_2
     end if
  endif

  if(present(momnt_int))then
     ! <2|d/dx|1>= -2b( <2|x_vec|1> - b_vec<2|1> )
     do m=1,2*la+1
        dipole(:,m,(/2,3,1/))=&    !(/2,3,1/) cause we add a vector
             & spread(bexp_arr,2,3)&
             & *(&
             &   -two*dipole(:,m,(/2,3,1/))&
             &   +spread(xa+xb-xd,1,num)&        !a_vec+b_vec-(+/-)(a_vec-b_vec),
             &   *spread(overlap(:,la*la+m),2,3)&!but for both cases: a<->b
             &  )
     enddo
     if(.not.laltlb)then
        do m=1,2*la+1
           momnt_int(:,:,1,m,1) = &
                unpack(dipole(:,m,2),cutoff,zero)
           momnt_int(:,:,1,m,2) = &
                unpack(dipole(:,m,3),cutoff,zero)
           momnt_int(:,:,1,m,3) = &
                unpack(dipole(:,m,1),cutoff,zero)
        enddo
     else
        ! the sign changes because
        ! <a|d/dx|b> == - <b|d/dx|a>
        do m=1,2*la+1
           momnt_int(:,:,m,1,1) = &
                unpack(-dipole(:,m,2),cutoff,zero)
           momnt_int(:,:,m,1,2) = &
                unpack(-dipole(:,m,3),cutoff,zero)
           momnt_int(:,:,m,1,3) = &
                unpack(-dipole(:,m,1),cutoff,zero)
        enddo
     end if
     if(debug)print *, 'dip_prim/ls: abs momnt=',sum(momnt_int**2)
  endif

  if(present(ovrl_int))then
     if(.not.laltlb)then
        i_lm=la*la
        do m=1,2*la+1
           i_lm=i_lm+1
           ovrl_int(:,:,1,m)=unpack(overlap(:,i_lm),cutoff,zero)
        enddo
     else
        i_lm=la*la
        do m=1,2*la+1
           i_lm=i_lm+1
           ovrl_int(:,:,m,1)=unpack(overlap(:,i_lm),cutoff,zero)
        enddo
     endif
  endif

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
       bexp_arr, &
       clmamb_scalar,&
       cutoff, &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ls: deallocation (1/2b) (2) failed")

end subroutine ls

!**************************************************
!***         LL_CALCULATE_(DIPOLE)              ***
!**************************************************

subroutine ll(bexps, aexps, lb, la, center2, center1, &
                               momnt_int,&
                               dip_int,&
                               ovrl_int)
  !
  !  Purpose: calculation of all primitive dipole integrals
  !       for a given set of primitives

  implicit none

  intent(in)  aexps, bexps,&
              la,lb,&
              center1,center2
  intent(out) momnt_int,dip_int,ovrl_int
  optional    momnt_int,dip_int,ovrl_int
  real(RK), dimension(:)       :: aexps, bexps
  integer(IK)                  :: la ! angular momentum of unique atom a
  integer(IK)                  :: lb ! angular momentum of unique atom b
  real(RK),dimension(3)        :: center1,center2 ! positions
  real(RK),dimension(:,:,:,:,:):: momnt_int,dip_int
  real(RK),dimension(:,:,:,:)  :: ovrl_int
  ! **** End of interface ***

  ! constants
  real(RK),parameter    :: two=2.0_RK
  real(RK),parameter    :: pi=3.14159265358979324_RK
  real(RK),parameter    :: gam=1.0_RK
  real(RK),parameter    :: very_small=1.0e-100_RK
  real(RK),parameter    :: very_big=1.0e100_RK
  real(RK),parameter    :: zero=0.0_RK
  real(RK),parameter,dimension(0:8) :: dfac= (/ &
       1.0_RK, 1.0_RK, 3.0_RK, 15.0_RK, 105.0_RK, &
       945.0_RK, 10395.0_RK, 135135.0_RK, 2027025.0_RK /)

  ! variables
  integer(IK) :: num,m,alloc_stat,i_xyz,i_l,i_lma,i_lmb, &
       n_lma,n_lmb,naexps,nbexps,ma,mb,l
  logical,allocatable   :: cutoff(:,:)
  real(RK),allocatable,dimension(:) :: &
       fact0,fact1,fact2,fact4,fact6,tau,aexp_arr,bexp_arr,clmamb_scalar
  real(RK),allocatable,dimension(:,:) :: &
       fact0_arr,fact1_arr,fact2_arr,clmamb,diff_arr0
  real(RK),allocatable,dimension(:,:,:) :: &
       overlap
  real(RK),allocatable,dimension(:,:,:,:) :: &
       dipole,diff_rule_result
  real(RK),dimension(3)  :: xa,xb,xd
  real(RK) :: arg

  if(debug)print *, 'dip_prim/ll: entered'

  if(la==0.or.lb==0)call error_handler&
       & ("ll: you`d better use ls_ or ss_ instead")

  naexps = size(aexps)
  nbexps = size(bexps)

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ll: allocation (1) failed")

  xa = center1
  xb = center2
  xd =xa-xb

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
     if(present(dip_int))   dip_int   = 0.0_RK
     if(present(momnt_int)) momnt_int = 0.0_RK
     if(present(ovrl_int))  ovrl_int = 0.0_RK
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("ll: deallocation (1/1) failed")
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
       diff_arr0((la+1)**2,(lb+1)**2),&
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("ll: allocation (2) failed")

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
       ("ll: deallocation  (1/2a) failed")

  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(max(la,lb),xd)
  fact4=1.0_RK
  i_lma=1
  tau=fact2*arg ! a*b/(a+b)*(A-B)**2
  do l=0,la
     do m=1,2*l+1
        clmamb(:,i_lma)=clmamb_scalar(i_lma)*fact4
        i_lma=i_lma+1
     enddo
     fact4=-fact4*fact2*2.0_RK
  enddo

  ! calculate overlap
  fact6= exp(-tau) * (4.0_RK*fact2/fact0)**0.75_RK  &
       / ( sqrt(aexp_arr**la*dfac(la)) * sqrt(bexp_arr**lb*dfac(lb)) )
  i_lmb=1
  do i_l=0,lb
     do mb=1,2*i_l+1
!!$        diff_arr0(:,i_lmb)= &
!!$             reshape( &
!!$             diff_rule( spread(clmamb_scalar,1,1), 1, n_lma, i_lmb ), &
!!$             (/n_lma/) &
!!$             )
        ! try scalar version:
        diff_arr0(:,i_lmb)=diff_rule(clmamb_scalar,1,n_lma,i_lmb)
        i_lmb=i_lmb+1
     end do
  end do
 
  i_lmb=1
  do i_l=0,lb
     do mb=1,2*i_l+1
        overlap(:,:,i_lmb)= &
             spread( fact6*(2.0_RK*fact2)**i_l, 2, n_lma ) * &
             prod_rule( &
             spread( diff_arr0(:,i_lmb), 1, num ), &
             clmamb(:,:), 1, n_lma &
             )
        i_lmb = i_lmb+1
     end do
  enddo

  if(present(ovrl_int))then
     ! return overlap:
     i_lma=la**2+1
     do ma=1,2*la+1
        i_lmb=lb**2+1
        do mb=1,2*lb+1
           ovrl_int(1:nbexps,1:naexps,mb,ma)=&
                & unpack(overlap(:,i_lma,i_lmb),cutoff,zero)
           i_lmb=i_lmb+1
        enddo
        i_lma=i_lma+1
     enddo
  endif
                  
  if(present(momnt_int))then
     ! complete and return <a|d/dx|b>:
     diff_rule_result = 0.0_RK
     diff_rule_result(:,1,1,(/2,3,1/)) =& ! (/x,y,z/) in solh terms
          & -two * spread(fact2,2,3) * spread(xa-xb,1,num)
     if(la>0)then
        do i_lma = 2, 4
           diff_rule_result(:,i_lma,1,i_lma-1) = -two * fact2
        enddo
     endif
     if(lb>0)then
        do i_lmb = 2, 4
           diff_rule_result(:,1,i_lmb,i_lmb-1) = +two * fact2
        enddo
     endif

     ! double product rule with respect to a and b
     ! dipole is not dipole any more, but <a|d/dx|b>
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
           momnt_int(:,:,mb,ma,1) = &
                unpack(dipole(:,ma,mb,2),cutoff,zero)
           momnt_int(:,:,mb,ma,2) = &
                unpack(dipole(:,ma,mb,3),cutoff,zero)
           momnt_int(:,:,mb,ma,3) = &
                unpack(dipole(:,ma,mb,1),cutoff,zero)
        end do
     end do
     if(debug)print *, 'dip_prim/ll: abs moment=',sum(momnt_int**2)
  endif

  if(present(dip_int))then
     ! complete and return dipole:
     ! calculate double differential rule on C(1,m) and multiply all scaling factors

     diff_rule_result = 0.0_RK
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
     if(la>0)then
        fact6 = aexp_arr / fact0
        do i_lma = 2, 4
           diff_rule_result(:,i_lma,1,i_lma-1) = fact6
        enddo
     endif
     if(lb>0)then
        fact6 = bexp_arr / fact0
        do i_lmb = 2, 4
           diff_rule_result(:,1,i_lmb,i_lmb-1) = fact6
        enddo
     endif

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
           dip_int(:,:,mb,ma,1) = &
                unpack(dipole(:,ma,mb,2),cutoff,zero)
           dip_int(:,:,mb,ma,2) = &
                unpack(dipole(:,ma,mb,3),cutoff,zero)
           dip_int(:,:,mb,ma,3) = &
                unpack(dipole(:,ma,mb,1),cutoff,zero)
        end do
     end do
  endif

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
       diff_arr0,&
       cutoff, &
       stat=alloc_stat)
!!$  write(*,*) 'll_calc...:dealloc . ok=',alloc_stat !debug
  if (alloc_stat.ne.0) call error_handler &
       ("ll: deallocation (1/2b) (2) failed")

end subroutine ll

end module dip_prim_module

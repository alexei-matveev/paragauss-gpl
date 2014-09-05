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
subroutine ll_calculate(na,nb,la,lb,imode,many_3c)
  !
  !  Purpose: calculation of all primitive 2 center orbital
  !           and 3 center integrals for a given set of indizes
  !       (unique_atom1,unique_atom2,la,lb).
  !       For three center integrals, contraction and symmetry-
  !       adaption concerning fitfunctions is also performed.
  !
  !
  !  Author: MS
  !  Date:   8/96
  !
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   11-12/99
  ! Description: integrals of electrostatic potential are added
  !
  ! Modification 
  ! Author: SB
  ! Date:   02/05
  ! Description: symmetrization for fit functions in 3c Coulomb integrals 
  !              was added
  !              NEW: l_fit_symmetry_adaption_v2 
  !              also coulomb calculations was rebuilded.
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  use unique_atom_module, noname=>pseudopot_present
  use gamma_module
  use type_module
  use datatype
  use solid_harmonics_module, only : solid_harmonics_calc,solid_harmonics_scalar
  use int_data_2cob3c_module  
  use solhrules_module
  use fitcontract_module
  use integralpar_module
  use pointcharge_module
  use options_module, only: options_integral_expmax
  use potential_module        
  use elec_static_field_module 
  use calc3c_switches,only: old_potential,old_3c_co,old_elfield
  use symmetry_data_module,only: symmetry_data_n_irreps,  &
       symmetry_data_n_partners,& 
       get_totalsymmetric_irrep 
!!$  use iounitadmin_module
  use shgi_cntrl, only: IPSEU
  implicit none

  !== Interrupt end of public interface of module ====================

  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
  integer(kind=i4_kind),intent(in) :: la ! angular momentum of unique atom a
  integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
  integer(kind=i8_kind),intent(in) :: imode ! for control
  real(r8_kind),optional,intent(out) :: many_3c(:,:,:,:,:)
  ! many_3c(nbexp,naexp,N_INTS*index3c,nlmb,nlma)
  ! Stored as illustrated by
  ! do ua=1,N_UA
  !   Int1[ua] stored at [3*(ua-1) + OFF_PVSP ] ! PVSP(ua)
  !   Int2[ua] stored at [3*(ua-1) + OFF_V    ] ! V(ua)
  !   Int3[ua] stored at [3*(ua-1) + OFF_VFIN ] ! V_{fin}(ua)
  !   ...
  ! enddo
  ! offsets OFF_* defined in int_data_2cob3c_module
  !===================================================================
  ! End of public interface of module
  !===================================================================

  integer(i4_kind), parameter :: &
       AxB = 1, &
       BxA = 2

  integer(kind=i4_kind) :: naexps,nbexps,ncexps
  real(kind=r8_kind),pointer :: aexps(:),bexps(:)
  real(kind=r8_kind),pointer     :: cexps(:)
  real(kind=r8_kind) :: z  ! charge
  real(kind=r8_kind) :: zc ! core charge
  integer(kind=i4_kind)          :: max_order

  ! constants
! real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
  real(kind=r8_kind),dimension(3,3) :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))
  real(kind=r8_kind),parameter :: &
       pi=3.14159265358979324_r8_kind, &
       very_small=1.0e-100_r8_kind, &
       very_big=1.0e100_r8_kind, &
       zero=0.0_r8_kind, &
       one=1.0_r8_kind, &
       two=2.0_r8_kind, &
       four=4.0_r8_kind, &
       six=6.0_r8_kind
! real(kind=r8_kind),parameter,dimension(0:8) :: dfac= &
  real(kind=r8_kind),dimension(0:8) :: dfac= &
       (/ 1.0_r8_kind, 1.0_r8_kind, 3.0_r8_kind, 15.0_r8_kind, 105.0_r8_kind, &
       945.0_r8_kind, 10395.0_r8_kind, 135135.0_r8_kind, 2027025.0_r8_kind /)
  integer(kind=i4_kind) :: one_i,zero_i
  integer(kind=i4_kind) :: num,counter,m,ma,mb,alloc_stat(40)=0
  integer(kind=i4_kind) :: memstat

  logical, allocatable   :: cutoff(:,:)

  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: &
       fact0_arr, fact1_arr, fact2_arr, fact10
  real(kind=r8_kind),allocatable,dimension(:)  :: &
       fact0, fact1, fact2, fact4, fact5, fact6, fact7, fact8, rcsabc, tau

  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:) :: gamma_arg, gamma_arg2, gamma_help

  ! help arrays for solid harmincs
  real(kind=r8_kind),allocatable  :: &
       yl_arr(:,:,:), yl_arr2(:,:,:), clmamb(:,:), clmamb2(:,:), clmamb_scalar(:)

  ! help arrays for product_rule and diff_rule
  real(kind=r8_kind),allocatable  :: &
       prod_arr(:,:,:,:,:), diff_arr(:,:,:), diff_arr0(:,:), &
       intermediate(:,:,:,:,:,:)

  real(kind=r8_kind) :: arg

  ! cartesian coordinates
  real(kind=r8_kind),dimension(3)  :: xa,xb,xc,xd

  integer(kind=i4_kind)  :: i,j,i_l,j_l,i_lb,k,i_ind,i_cnt,l !,l_cf,l_j
  integer(kind=i4_kind)  :: lmax_ch,lmax_xc,lmax_abs,ly_max
  integer(kind=i4_kind)  :: n_equals,n_independent_fcts,n_contributing_fcts

  integer(kind=i4_kind), pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind), pointer      :: coeff(:)

  real(kind=r8_kind),allocatable :: aexp_arr(:),bexp_arr(:)
  real(kind=r8_kind),allocatable :: nested2_fac1(:,:),nested2_fac2(:,:),nested2_fac12(:,:,:)
  
  real(kind=r8_kind)               :: expmax
  type(unique_atom_type), pointer  :: ua_pointer

  ! the calculated integrals
  real(kind=r8_kind),allocatable  :: potential(:,:,:,:), field(:,:,:,:), intermed_3c(:)  !!!!!!!!!!!!
  real(kind=r8_kind),allocatable  :: prod_arr_gr(:,:,:,:,:),help_vec(:),prod_arr_gr_vec(:,:,:), &
       help_mat(:,:,:)
  real(kind=r8_kind),allocatable :: diff_arr_xyz(:,:,:,:),yl_arr_xyz(:,:,:)

  real(kind=r8_kind),allocatable  :: &
       overlap(:,:,:), kinetic(:), &
       nuc(:,:,:), nuc_pseudo(:,:,:)
  logical                         :: pseudopot_present ! same name as in UA module
  real(kind=r8_kind),allocatable  :: nuc_pc_timps(:,:,:)
  type(three_center_l)            :: xc_int !! coul_int 

  !! SB: new store for coul_int
  !! coul_int(n_irreps)%l(-1:lmx_co)%m(num,ncexps,n_if,mb,ma,n_pa)
  type(three_center_l_v2),allocatable :: coul_int(:)
  !! end of new stores

  type nested2_opt
   integer(kind=i4_kind):: n
!???   type(nested2_vars),allocatable,dimension(:)::summands 
   type(nested2_vars),pointer,dimension(:)::summands 
  end type nested2_opt
  type(nested2_opt),allocatable,dimension(:):: nested2_summands
  integer(kind=i4_kind):: nested2_l1_max,nested2_l3_max,l1_max,l3_max
  logical:: opt_nested2=.true.

  logical :: split3c
  real(r8_kind),allocatable  :: this(:,:,:) ! (num,nlmA,nlmB)
!!$  real(r8_kind),allocatable  :: this5d(:,:,:,:,:) ! (num,1,1,nlmB,nlmA)
  real(r8_kind),allocatable  :: this6d(:,:,:,:,:,:) ! (num,1,1,nlmB,nlmA,1)
  ! should it be like that --- AxB, BxA ?! ...
  real(r8_kind) :: zexps(1) ! for finite nuc only
  logical :: with_timps
  integer(i4_kind) :: N_length
  FPP_TIMER_DECL(pll)
  integer(i4_kind) :: i_ir, i_pa, n_pa

  intrinsic max

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'll_calculate: PP=',pseudopot_present,' imode=',imode

! print*, 'in ll_calc'

  split3c = present(many_3c)

      if(opt_nested2) then
       allocate(nested2_summands((lb+1)**2),stat=alloc_stat(1))
       if(alloc_stat(1).ne.0) call error_handler('nested2_summands not allocated')
       alloc_stat(1)=1
       counter=0
       nested2_l1_max=0
       nested2_l3_max=0
       do i_l=0,lb
          do mb=1,2*i_l+1
             counter=counter+1
	  nested2_summands(counter)%n=nsum_prod_rule_nested2(la,i_l,mb)
!          print*,nested2_summands(counter)%n,'nested2_summands n'
	  allocate (nested2_summands(counter)%summands(nested2_summands(counter)%n), &
	            stat=alloc_stat(2))
          if(alloc_stat(2).ne.0) call error_handler('2 nested2_summands%summands allocate failed')
	  alloc_stat(2)=1
          call summands_prod_rule_nested2(la,i_l,mb,nested2_summands(counter)%summands,l1_max,l3_max)
          if(l1_max.gt.nested2_l1_max) nested2_l1_max=l1_max
          if(l3_max.gt.nested2_l3_max) nested2_l3_max=l3_max
          end do
       end do
      endif




  one_i=1_i4_kind
  zero_i=0_i4_kind
  naexps = unique_atoms(na)%l_ob(la)%n_exponents
  nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat(3))
  if (alloc_stat(3).ne.0) call error_handler &
       ("LL_CALCULATE: allocation (1) failed") 
	alloc_stat(3)=1
        alloc_stat(4)=1 !cutoff


  xa = center1
  xb = center2
  xd =xa-xb
  aexps => unique_atoms(na)%l_ob(la)%exponents(:)
  bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)

  arg=sum(xd**2)

  fact0_arr=(spread(aexps,1,nbexps)+spread(bexps,2,naexps))
  fact1_arr=(spread(aexps,1,nbexps)*spread(bexps,2,naexps))

  where(fact0_arr>=very_small) ! prevent division by zero
     fact2_arr=fact1_arr/fact0_arr
  elsewhere
     fact2_arr=very_big
  end where

  expmax = options_integral_expmax()
  where(fact2_arr*arg> expmax ) ! cutoff: where almost no overlap
     cutoff=.false.              ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where
  num=count(cutoff)
  if(num==0) then ! all integrals are equal zero
     if (integralpar_2cob_ol) then
        prim_int_2cob_ol = 0.0_r8_kind
     end if
     if (integralpar_2cob_kin) then
        prim_int_2cob_kin= 0.0_r8_kind
     end if
     if (integralpar_2cob_nuc) then
        prim_int_2cob_nuc(:,:,:,:)=0.0_r8_kind
     end if
     if (integralpar_2cob_potential) then
        prim_int_2cob_poten(:,:,:,:,:)=0.0_r8_kind      !!!!!!!!!!!!!!!
     end if 
     if (integralpar_2cob_field) then
        prim_int_2cob_field(:,:,:,:,:)=0.0_r8_kind      !!!!!!!!!!!!!!!
     end if
     if(integralpar_relativistic)then
        prim_int_2cob_pvsp = zero
     endif
     if(integralpar_3c_co) then
        prim_int_3c_co=0.0_r8_kind
     end if
     if(integralpar_3c_xc) then
        prim_int_3c_xc=0.0_r8_kind
     end if
     if( split3c )then
        ASSERT(integralpar_relativistic)
        many_3c = zero
     endif
     deallocate(fact0_arr,fact1_arr,&
          fact2_arr,cutoff,stat=alloc_stat(3))
     if (alloc_stat(3).ne.0) call error_handler &
          ("LL_CALCULATE: deallocation (1) failed")
	alloc_stat(4)=0
     !return
     goto 999 ! clean up and exit
  end if

  allocate (&
       fact0(num),fact1(num),fact2(num),fact4(num),fact5(num),&
       fact6(num),fact7(num),fact8(num),rcsabc(num),tau(num),& !tau 6
       gamma_arg(num,3),aexp_arr(num),bexp_arr(num),& 
       overlap(num,(la+1)**2,(lb+1)**2),& !5
       clmamb_scalar((max(la,lb)+1)**2),& !5
       clmamb(num,(la+1)**2),& !5
       clmamb2(num,(la+1)**2),& !6
       diff_arr(num,(la+1)**2,(lb+1)**2),& !5
       diff_arr0((la+1)**2,(lb+1)**2),& !5
       stat=alloc_stat(5))
  if (alloc_stat(5).ne.0) call error_handler &
       ("LL_CALCULATE: allocation (2) failed")
  alloc_stat(5)=1
  alloc_stat(6)=1 !tau clmamb2

  ! AxB
  allocate( this(num,2*la+1,2*lb+1), stat=memstat)
  ASSERT(memstat==0)
  if(split3c)then
     ! BxA
!!$     allocate( this5d(num,1,1,2*lb+1,2*la+1), stat=memstat)
!!$     ASSERT(memstat==0)
     allocate( this6d(num,1,1,2*lb+1,2*la+1,1), stat=memstat)
     ASSERT(memstat==0)
  endif

  if (integralpar_2cob_kin) then
     allocate(kinetic(num),stat=alloc_stat(7))
     if (alloc_stat(7).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (3) failed")
	alloc_stat(7)=1
  end if

  if (integralpar_2cob_nuc) then
     allocate(nuc(num,2*la+1,2*lb+1), stat=alloc_stat(8))
     if (alloc_stat(8).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (4) failed")
	alloc_stat(8)=1
     nuc=0.0_r8_kind
     if (pseudopot_present) then
        allocate(nuc_pseudo(num,2*la+1,2*lb+1), &
             stat=alloc_stat(9))
        if (alloc_stat(9).ne.0) call error_handler &
             ("LL_CALCULATE: allocation (4) failed")
	alloc_stat(9)=1
        nuc_pseudo = 0.0_r8_kind
        with_timps = pointcharge_N+n_timps .gt. 0
        if (with_timps .and. integralpar_relativistic)  then
	   allocate(nuc_pc_timps(num,2*la+1,2*lb+1), &
                stat=alloc_stat(10))
           if (alloc_stat(10).ne.0) call error_handler &
                ("LL_CALCULATE: allocation nuc_pc_timps failed")
	   alloc_stat(10)=1
           nuc_pc_timps = 0.0_r8_kind
        endif
     end if ! pseudopot_present
  end if ! integralpar_2cob_nuc

  if (integralpar_2cob_potential) then
     allocate(potential(N_points,num,2*la+1,2*lb+1), stat=alloc_stat(11))  
     if (alloc_stat(11).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (5) failed")
     potential=0.0_r8_kind
     allocate(intermed_3c(num), stat=alloc_stat(11))  
     if (alloc_stat(11).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (5a) failed")
	alloc_stat(11)=1
  end if
#if 0
  if (integralpar_2cob_field) then
     if(calc_normal) then
        N_length=N_surface_points
     else
        N_length=totsym_field_length
     end if
     allocate(field(N_length,num,2*lb+1,2*la+1), stat=alloc_stat(11)) 
     if (alloc_stat(11).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (6) failed")
     field=0.0_r8_kind
     allocate(intermed_3c(num), stat=alloc_stat(11)) 
     if (alloc_stat(11).ne.0) call error_handler &
          ("LL_CALCULATE: allocation (6a) failed")
	alloc_stat(11)=1
  end if
#endif

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  ! fact7= 1/sqrt(a**l*(2l-1)!!)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)

  aexp_arr=pack(spread(aexps,1,nbexps),cutoff)
  bexp_arr=pack(spread(bexps,2,naexps),cutoff)

  if(opt_nested2) then
   allocate(nested2_fac1(size(aexp_arr,1),0:nested2_l3_max), &
                nested2_fac2(size(bexp_arr,1),0:nested2_l1_max),stat=alloc_stat(12))
   if (alloc_stat(12)/=0) call error_handler("allocation nested2_fac failed")
   allocate(nested2_fac12(size(aexp_arr,1),0:nested2_l3_max,0:nested2_l1_max), &
                                                                 stat=alloc_stat(12))
   if (alloc_stat(12)/=0) call error_handler("allocation nested2_fac12 failed")
       alloc_stat(12)=1
	nested2_fac1(:,0)=1.0_r8_kind
	nested2_fac2(:,0)=1.0_r8_kind
	do i_l=1,nested2_l3_max
        nested2_fac1(:,i_l)=nested2_fac1(:,i_l-1)*(-2.0_r8_kind*aexp_arr)
        enddo
	do i_l=1,nested2_l1_max
        nested2_fac2(:,i_l)=nested2_fac2(:,i_l-1)*(-2.0_r8_kind*bexp_arr)
        enddo
	do i_l=0,nested2_l3_max
	do j_l=0,nested2_l1_max
         nested2_fac12(:,i_l,j_l)=nested2_fac1(:,i_l)*nested2_fac2(:,j_l)
	enddo
	enddo
  endif    
       

  deallocate(fact0_arr,fact1_arr,fact2_arr,stat=alloc_stat(3))
  if (alloc_stat(3)/=0) call error_handler &
       ("LL_CALCULATE: deallocation (2) failed")

  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  gamma_arg(:,1)=(pack(spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps),cutoff))/fact0

  gamma_arg(:,2)=(pack(spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps),cutoff))/fact0

  gamma_arg(:,3)=(pack(spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps),cutoff))/fact0


  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(max(la,lb),xd)
  fact4=1.0_r8_kind
  counter=1
  tau=fact2*arg ! a*b/(a+b)*(A-B)**2
  do l=0,la
     do m=1,2*l+1
        clmamb(:,counter)=clmamb_scalar(counter)*fact4
        clmamb2(:,counter)=clmamb_scalar(counter)&
             *fact4*(tau-real(l,kind=r8_kind))
        counter=counter+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo


  ! first calculating 2-center integrals----------------
  ! fact5=fact2*(3.0_r8_kind-2.0_r8_kind*tau+2.0_r8_kind*la)  
  ! a*b/(a+b)(3-2*tau+2*l)
  fact6=1.0_r8_kind/sqrt(aexp_arr**la*dfac(la))/&
       sqrt(bexp_arr**lb*dfac(lb))*exp(-tau)*&
       (4.0_r8_kind*fact2/fact0)**0.75_r8_kind  
  fact5=fact2*fact6
  fact7=(fact2*2.0_r8_kind)**lb
  counter=1
  do i_l=0,lb
     do mb=1,2*i_l+1
        diff_arr0(:,counter)=reshape(diff_rule(spread(clmamb_scalar,1,1),&
             1,(la+1)**2,counter),(/(la+1)**2/))
        counter=counter+1
     end do
  end do
  counter=1
  do i_l=0,lb
     magnetic_number_b: do mb=1,2*i_l+1
        ! overlap
        overlap(:,1:(la+1)**2,counter)=spread(fact6*(2.0_r8_kind*fact2)**i_l,&
             2,(la+1)**2)*&
             prod_rule(spread(diff_arr0&
             (:,counter),1,num),clmamb(:,:),1,&
             (la+1)**2)
        counter=counter+1
     end do magnetic_number_b
  enddo

  if (integralpar_2cob_ol) then
     do ma=1,2*la+1
        do mb=1,2*lb+1
           prim_int_2cob_ol(:,:,mb,ma) = unpack&
                (overlap(:,la**2+ma,lb**2+mb),cutoff,zero)
        end do
     end do
  end if

  if (integralpar_2cob_kin) then
     do ma=1,2*la+1
        do mb=1,2*lb+1
           ! kinetic energy
           kinetic=fact5*fact7*reshape(prod_rule(spread(diff_arr0&
                (:,lb**2+mb),1,num),(3.0_r8_kind+2.0_r8_kind*lb)&
                *clmamb(:,:)-2.0_r8_kind*clmamb2,(la)**2+ma,&
                (la)**2+ma),(/num/))
           ! re-map them to the int_data_2cob3c_stuff
           prim_int_2cob_kin(:,:,mb,ma)= unpack(kinetic,cutoff,zero)
        end do
     end do
     deallocate(kinetic,STAT=alloc_stat(7))
     if (alloc_stat(7).ne.0) call error_handler &
          ("LL_CACLULATE : deallocation (3) failed")
  endif

  call integral_interrupt_2cob3c()

  deallocate(clmamb2,tau,stat=alloc_stat(6))
  if (alloc_stat(6).ne.0) call error_handler &
       ("LL_CACLULATE : deallocation (4) failed")


  ! calculate integrals that involve third center,
  ! i.e. fit integrals, nuclear attraction and relativistic pv scalar p
  third_center_required: if ( integralpar_2cob_nuc .or. integralpar_3c_xc &
       .or. integralpar_3c_co .or. integralpar_relativistic &
!!! MF merge bug fix
!      .or. integralpar_2cob_potential) then
     ) then

     fact8=2.0_r8_kind*sqrt(fact0/pi)

     unique_atom_loop: do i=1,n_unique_atoms + n_timps  ! loop over third center
        if(i<=n_unique_atoms) then
           ua_pointer=>unique_atoms(i)
           lmax_ch= ua_pointer%lmax_ch      ! maximum l  for chargefit
           lmax_xc= ua_pointer%lmax_xc      ! maximum l  for xcfit  
           ! determine the maximal angular momentum

           ly_max=max(la,lb,lmax_ch,lmax_xc)
           if (.not.integralpar_3c_xc) then
              lmax_abs=lmax_ch
           else
              lmax_abs=max(lmax_ch,lmax_xc)
           endif
           ly_max=max(la,lb,lmax_ch,lmax_xc)  
           max_order=max(1+la+lb+lmax_abs,3+la+lb)
           z= ua_pointer%z                  ! charge 
           zc= ua_pointer%zc                ! core charge
           n_equals=ua_pointer%n_equal_atoms
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'll_calc: ua=',i,', zero its charge!'
        zc = zero
        z  = zero
           allocate ( &
                gamma_help(num,max_order), &
                gamma_arg2(num,n_equals), &
                stat=alloc_stat(13))
           if (alloc_stat(13)/=0) call error_handler &
                ("LL_CACLULATE : allocation (5) failed")
	   alloc_stat(13)=1


           ! --- further allocation ----------------------------------
           ! num : number of pairs(a,b) which are inside the cutoff
           ! for s-and r2-type there is only 1 indep. fct
           if(integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
              allocate (coul_int(symmetry_data_n_irreps()),stat=alloc_stat(14))

              i_ir_alloc_: DO i_ir=1,symmetry_data_n_irreps() !!allocation for coul_int

                 n_pa = symmetry_data_n_partners(i_ir)

                 allocate (coul_int(i_ir)%l(-1:lmax_ch),stat=alloc_stat(14))
                 if (alloc_stat(14)/=0) call error_handler &
                      ("LL_CACLULATE : allocation coul_int%l failed")
                 ncexps = ua_pointer%r2_ch%n_exponents

                 n_independent_fcts  = &
                         ua_pointer%symadapt_partner(i_ir,0)%n_independent_fcts 

                 allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,n_independent_fcts,&
                      2*lb+1,2*la+1,n_pa),stat=alloc_stat(14))
                 if (alloc_stat(14)/=0) call error_handler &
                      ("LL_CACLULATE : allocation coul_int%l(-1)%m failed")

                 ncexps = ua_pointer%l_ch(0)%n_exponents
                 allocate(coul_int(i_ir)%l(0)%m(num,ncexps,n_independent_fcts,&
                      2*lb+1,2*la+1,n_pa),stat=alloc_stat(14))
                 if (alloc_stat(14)/=0) call error_handler &
                      ("LL_CACLULATE : allocation  coul_int%l(0)%m failed")
                 alloc_stat(14)=1
                 do i_l=1,lmax_ch
                    ncexps = ua_pointer%l_ch(i_l)%n_exponents
                    n_independent_fcts  = &
                         ua_pointer%symadapt_partner(i_ir,i_l)%n_independent_fcts 
                    allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                         2*lb+1,2*la+1,n_pa),&
                         stat=alloc_stat(14))
                    if (alloc_stat(14)/=0) call error_handler &
                         ("LL_CACLULATE : allocation (8) failed")
                    alloc_stat(14)=1
                 end do

                 do i_l = -1,lmax_ch
                    coul_int(i_ir)%l(i_l)%m = 0.0_r8_kind
                 end do

              END DO i_ir_alloc_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
           elseif (integralpar_3c_co) then

              i_ir = get_totalsymmetric_irrep()
              allocate (coul_int(i_ir),stat=alloc_stat(14))
              n_pa = 1

              allocate (coul_int(i_ir)%l(-1:lmax_ch),stat=alloc_stat(14))
              if (alloc_stat(14)/=0) call error_handler &
                   ("LL_CACLULATE : allocation coul_int%l failed")
              ncexps = ua_pointer%r2_ch%n_exponents
              allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,1,2*lb+1,2*la+1,n_pa),stat=alloc_stat(14))
              if (alloc_stat(14)/=0) call error_handler &
                   ("LL_CACLULATE : allocation coul_int%l(-1)%m failed")
 
              ncexps = ua_pointer%l_ch(0)%n_exponents
              allocate(coul_int(i_ir)%l(0)%m(num,ncexps,1,2*lb+1,2*la+1,n_pa),stat=alloc_stat(14))
              if (alloc_stat(14)/=0) call error_handler &
                   ("LL_CACLULATE : allocation  coul_int%l(0)%m failed")
	      alloc_stat(14)=1
              do i_l=1,lmax_ch
                 ncexps = ua_pointer%l_ch(i_l)%n_exponents
                 n_independent_fcts  = &
                      ua_pointer%symadapt_partner(i_ir,i_l)%n_independent_fcts 
                 allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                      2*lb+1,2*la+1,n_pa),&
                      stat=alloc_stat(14))
                 if (alloc_stat(14)/=0) call error_handler &
                      ("LL_CACLULATE : allocation (8) failed")
                 alloc_stat(14)=1
              end do

              do i_l = -1,lmax_ch
                 coul_int(i_ir)%l(i_l)%m = 0.0_r8_kind
              end do

           end if

           if(integralpar_3c_xc) then
              allocate (xc_int%l(-1:lmax_xc),stat=alloc_stat(15)) 
              if (alloc_stat(15)/=0) call error_handler &
                   ("LL_CACLULATE : allocation xc_int%l failed")
              ncexps = ua_pointer%r2_xc%n_exponents
              allocate(xc_int%l(-1)%m(num,ncexps,1,2*lb+1,2*la+1),stat=alloc_stat(15))
              if (alloc_stat(15)/=0) call error_handler &
                   ("LL_CACLULATE : allocation xc_int%l(-1)%m failed")
              ncexps = ua_pointer%l_xc(0)%n_exponents
              allocate(xc_int%l(0)%m(num,ncexps,1,2*lb+1,2*la+1),stat=alloc_stat(15))
              if (alloc_stat(15)/=0) call error_handler &
                   ("LL_CACLULATE : allocation xc_int%l(0)%m failed")
		alloc_stat(15)=1
              xc_int%l(0)%m  = 0.0_r8_kind
              xc_int%l(-1)%m = 0.0_r8_kind
           endif

           ly_max=max(la,lb,lmax_ch,lmax_xc)       
           allocate( &
                yl_arr(num,(ly_max+1)**2,n_equals),&
                fact10(num,(ly_max+1)**2),&
                yl_arr2(num,(ly_max+1)**2,n_equals),&
                prod_arr(num,(la+1)**2,(lb+1)**2,0:la+lb,n_equals),&
                stat=alloc_stat(16))
           if (alloc_stat(16)/=0) call error_handler &
                ("LL_CACLULATE : allocation (6) failed")
	   alloc_stat(16)=1
           prod_arr=0.0_r8_kind
           ! do a precalculation of a factor needed for the
           ! product rule
           counter=1
           fact4=1.0_r8_kind
           do i_l=0,ly_max
              do ma=1,2*i_l+1
                 fact10(:,counter)=fact4
                 counter=counter+1
              enddo
              fact4=fact4*aexp_arr/(fact0)
           enddo

           ! precalculate prod_arr and calculate nuclear attraction
           call precalculate_and_nuc(this)
           if( integralpar_2cob_nuc )then
              nuc = nuc + this
           endif

           if( split3c )then
              call unpack_many_3c(this,i,OFF_V,form=AxB)
           end if

           !
           ! now calculating fit integrals
           !

           ! XC part was not changed, but separated

           ! s-type xc fit integrals
           if(integralpar_3c_xc) call s_xc()
           
           ! r2-type exchange fit integral
           if(integralpar_3c_xc) call r2_xc()

           ! l_type coloumb and exchange fit integral
           if(integralpar_3c_xc .or. integralpar_3c_co) then
              do i_l=1,lmax_abs
                 n_independent_fcts  = &
                      ua_pointer%symadapt_partner(1,i_l)%n_independent_fcts
              
                 allocate( &
                      intermediate(num,2*la+1,2*lb+1,n_independent_fcts,n_equals,0:la+lb) &
                      ,stat=alloc_stat(17))
                 if (alloc_stat(17)/=0) call error_handler('LL_CALCULATE: allocation (7) failed')
		     alloc_stat(17)=1

                 ! symmetry adaption for l-type xc fit integrals
                 ! results are stored in intermediate array
                 call l_fit_symmetry_adapt()

                 ! now the same for exchange
                 if(integralpar_3c_xc.and.lmax_xc>=i_l) call l_xc()

                 deallocate(intermediate,stat=alloc_stat(17))
                 if (alloc_stat(17)/=0) call error_handler &
                      ("LL_CACLULATE : deallocation (5) failed")  
              end do! loop over lc
              ! finished with l-type fit integrals
           endif

           ! END of XC part

           ! Coulomb part was rebuilded for symmetry

           if (integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
              i_ir_: do i_ir=1,symmetry_data_n_irreps()
                 i_pa_: do i_pa=1,symmetry_data_n_partners(i_ir)
                    i_l = 0 ! r2 and s i_l=0

                    n_independent_fcts  = &
                         unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts  

                    if (n_independent_fcts .ne. 0) then

                       ! s-type coulomb fit integrals
                       call s_coulomb( &
                            unique_atoms(i)%l_ch(0)%exponents(:), &
                            coul_int(i_ir)%l(0)%m &
                            )

                       if( split3c )then
                          zexps(1) = (3.0_r8_kind/2.0_r8_kind) / unique_atoms(i)%nuclear_radius**2
                          call s_coulomb( zexps, this6d )
                          call unpack_many_3c(this6d(:,1,1,:,:,1),i,OFF_VFIN,form=BxA)
                       end if

                       ! r2-type coloumb fit integral
                       call r2_coulomb()

                    end if

                    ! l_type coloumb fit integral
                    do i_l=1,lmax_ch

                       n_independent_fcts  = &
                            unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts 

                       if (n_independent_fcts .ne. 0) then
                          allocate( &
                               intermediate(num,2*la+1,2*lb+1,n_independent_fcts,n_equals,0:la+lb) &
                               ,stat=alloc_stat(17))
                          if (alloc_stat(17)/=0) call error_handler('LL_CALCULATE: allocation (7) failed')
                          alloc_stat(17)=1

                          ! symmetry adaption for l-type charge fit integrals
                          ! results are stored in intermediate array
                          call l_fit_symmetry_adapt_v2(i,i_l,i_ir,i_pa,intermediate)

                          ! coulomb integrals
                          call l_coulomb()

                          deallocate(intermediate,stat=alloc_stat(17))
                          if (alloc_stat(17)/=0) call error_handler &
                               ("LL_CACLULATE : deallocation (5) failed")
                       end if
                       
                    end do! loop over lc
                       ! finished with l-type fit integrals

                 end do i_pa_
              end do i_ir_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif

           elseif (integralpar_3c_co) then

              i_ir = get_totalsymmetric_irrep()
              i_pa = 1

              i_l = 0 ! r2 and s i_l=0

              n_independent_fcts  = &
                   unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts  

              if (n_independent_fcts .ne. 0) then
!!$                 allocate( &
!!$                      intermediate(num,2*la+1,2*lb+1,n_independent_fcts,n_equals,0:la+lb) &
!!$                      ,stat=alloc_stat(17))
!!$                 if (alloc_stat(17)/=0) call error_handler('LL_CALCULATE: allocation (7) failed')
!!$                 alloc_stat(17)=1

!!$                 call l_fit_symmetry_adapt_v2(i,i_l,i_ir,i_pa,intermediate)

                 ! s-type coulomb fit integrals
                 call s_coulomb( &
                      unique_atoms(i)%l_ch(0)%exponents(:), &
                      coul_int(i_ir)%l(0)%m &
                      )

                 if( split3c )then
                    zexps(1) = (3.0_r8_kind/2.0_r8_kind) / unique_atoms(i)%nuclear_radius**2
                    call s_coulomb( zexps, this6d )
                    call unpack_many_3c(this6d(:,1,1,:,:,1),i,OFF_VFIN,form=BxA)
                 end if

                 ! r2-type coloumb fit integral
                 call r2_coulomb()

!!$                 deallocate(intermediate,stat=alloc_stat(17))
!!$                 if (alloc_stat(17)/=0) call error_handler &
!!$                      ("LL_CACLULATE : deallocation (5) failed")
              end if

              ! l_type coloumb fit integral
              if(integralpar_3c_co .and. old_3c_co) then
                 do i_l=1,lmax_ch

                    n_independent_fcts  = &
                         unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts 

                    if (n_independent_fcts .ne. 0) then
                       allocate( &
                            intermediate(num,2*la+1,2*lb+1,n_independent_fcts,n_equals,0:la+lb) &
                            ,stat=alloc_stat(17))
                       if (alloc_stat(17)/=0) call error_handler('LL_CALCULATE: allocation (7) failed')
                       alloc_stat(17)=1

                       ! symmetry adaption for l-type charge fit integrals
                       ! results are stored in intermediate array
                       call l_fit_symmetry_adapt_v2(i,i_l,i_ir,i_pa,intermediate)

                       ! coulomb integrals
                       call l_coulomb()

                       deallocate(intermediate,stat=alloc_stat(17))
                       if (alloc_stat(17)/=0) call error_handler &
                            ("LL_CACLULATE : deallocation (5) failed")
                    end if

                 end do! loop over lc
                 ! finished with l-type fit integrals
              endif

           end if


           ! contract the fit integrals with respect to fit dimension
           ! and write them to their final location in int_data_2cob3c_module

           if(integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE              
              call fitcontract_v2(num,i,cutoff,coul_int)            
              do i_ir=1,symmetry_data_n_irreps()              
                 do i_l = -1, lmax_ch
                    deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat(14))
                    if(alloc_stat(14).ne.0) call error_handler &
                         ("LL_CALCULATE : deallocation coul_int%l%m failed")
                 enddo
                 deallocate (coul_int(i_ir)%l,STAT=alloc_stat(14))
                 if(alloc_stat(14).ne.0) call error_handler &
                      ("LL_CALCULATE : deallocation coul_int%l failed")
              end do
              deallocate(coul_int,STAT=alloc_stat(14))
              if(alloc_stat(14).ne.0) call error_handler &
                   ("LL_CALCULATE : deallocation coul_int%l failed")
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif

           elseif(integralpar_3c_co) then
              call fitcontract_v2(num,i,cutoff,coul_int)
              i_ir=get_totalsymmetric_irrep() 
              do i_l = -1, lmax_ch
                 deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat(14))
                 if(alloc_stat(14).ne.0) call error_handler &
                      ("LL_CALCULATE : deallocation coul_int%l%m failed")
              enddo
              deallocate (coul_int(i_ir)%l,STAT=alloc_stat(14))
              if(alloc_stat(14).ne.0) call error_handler &
                   ("LL_CALCULATE : deallocation coul_int%l failed")
              deallocate(coul_int,STAT=alloc_stat(14))
              if(alloc_stat(14).ne.0) call error_handler &
                   ("LL_CALCULATE : deallocation coul_int%l failed")
           endif
           if(integralpar_3c_xc) then
              call fitcontract('xc',num,i,cutoff,xc_int)
              do i_l = -1, lmax_xc
                 deallocate(xc_int%l(i_l)%m,STAT=alloc_stat(15)) 
                 if(alloc_stat(15).ne.0) call error_handler &
                      ("LL_CALCULATE : deallocation xc_int%l%m failed")
              enddo
              deallocate (xc_int%l,STAT=alloc_stat(15))
              if(alloc_stat(15).ne.0) call error_handler &
                   ("LL_CALCULATE : deallocation xc_int%l failed")
           end if

           deallocate(yl_arr,yl_arr2,gamma_help,gamma_arg2,fact10,&
                prod_arr,stat=alloc_stat(13)) !gamma_help,gamma_arg2
           if (alloc_stat(13)/=0) call error_handler &
                ("LL_CALCULATE : deallocation (7) failed")
		alloc_stat(16)=0 ! yl_arr,yl_arr2 fact10 prod_arr
           !  end do unique_atom_loop
        else !timp
           ua_pointer=>unique_timps(i-n_unique_atoms)
           z= ua_pointer%z                  ! charge 
           zc= ua_pointer%zc                ! core charge
           n_equals=ua_pointer%n_equal_atoms
        end if
        if(zc/=0.0_r8_kind .and. .not.integralpar_2cob_potential) then  ! pseudopotential contributions
           ABORT('not supported')
        endif ! end of pseudopotential contributions

     end do unique_atom_loop

     ! add contribution of point charges to nuclear attraction
     if ( (integralpar_2cob_nuc .or. integralpar_relativistic) &
          .and. pointcharge_N+n_timps .gt. 0) call add_pointcharges()         

  end if third_center_required
  DPRINT  'TIMER: ll_pseudo=',FPP_TIMER_VALUE(pll)

  if (integralpar_2cob_potential.and.old_potential) call calc_potential()             !!!!!!!!!!!!!!!!1
#if 0
  if (integralpar_2cob_field ) &
       call calc_field()                 !!!!!!!!!!!!!!!!!
#endif

  if (integralpar_2cob_nuc) then
     if (pseudopot_present &
          .and.(.not.integralpar_relativistic)) then
        do mb=1,2*lb+1
           do ma=1,2*la+1 
              prim_int_2cob_nuc(:,:,mb,ma)= &
                   unpack(nuc(:,ma,mb)+nuc_pseudo(:,ma,mb),cutoff,zero)
           enddo
        end do
     else ! i.e. relativistic
        do mb=1,2*lb+1
           do ma=1,2*la+1
              prim_int_2cob_nuc(:,:,mb,ma)=unpack(nuc(:,ma,mb),cutoff,zero)
           enddo
        enddo
        if(pseudopot_present) then
           do mb=1,2*lb+1
              do ma=1,2*la+1
                 prim_int_2cob_nuc_pseudo(:,:,mb,ma)=unpack(nuc_pseudo(:,ma,mb),cutoff,zero)
              enddo
           enddo
           with_timps = pointcharge_N+n_timps .gt. 0
           if(with_timps) then
              do mb=1,2*lb+1
                 do ma=1,2*la+1
                    prim_int_2cob_nuc_pseudo(:,:,mb,ma)=prim_int_2cob_nuc_pseudo(:,:,mb,ma)&
                         + unpack(nuc_pc_timps(:,ma,mb),cutoff,zero)
                 enddo
              enddo
           endif
        endif ! pseudopot_present
     endif ! else of PP and not rel

     deallocate(nuc,stat=alloc_stat(8))
     if(alloc_stat(8)/=0) call error_handler &
          ("LL_CALCULATE: deallocation (9) failed")

     if (pseudopot_present) then
        deallocate(nuc_pseudo,stat=alloc_stat(9))
        if(alloc_stat(9)/=0) call error_handler &
             ("LL_CALCULATE: deallocation (9) failed")

        with_timps = pointcharge_N+n_timps .gt. 0
        if (with_timps .and. integralpar_relativistic)  then
	   deallocate(nuc_pc_timps, stat=alloc_stat(10))
           if (alloc_stat(10).ne.0) call error_handler &
                ("LL_CALCULATE: deallocation nuc_pc_timps failed")
        endif ! 
     end if ! pseudopot_present
!!!     call print_nuclear()
  end if ! integralpar_2cob_nuc

!      print*,sum(prim_int_2cob_nuc),na,la,nb,lb,'sum prim_int_2cob_nuc'

  if (integralpar_2cob_potential) then
     do i=1,N_points
        do mb=1,2*lb+1
           do ma=1,2*la+1
              !restoring needed for VPP
              intermed_3c(:)=potential(i,:,ma,mb)
              prim_int_2cob_poten(:,:,i,mb,ma)=unpack(intermed_3c(:),cutoff,zero) 
           enddo
        enddo
     end do
     deallocate(potential,stat=alloc_stat(11))
     if(alloc_stat(11)/=0) call error_handler &
          ("LL_CALCULATE: deallocation (10) failed")
     deallocate(intermed_3c,stat=alloc_stat(11))
     if(alloc_stat(11)/=0) call error_handler &
          ("LL_CALCULATE: deallocation (10a) failed")
  end if
#if 0
  if (integralpar_2cob_field) then
     if(calc_normal) then
        N_length=N_surface_points
     else
        N_length=totsym_field_length
     end if
     do i=1,N_length
        do mb=1,2*lb+1
           do ma=1,2*la+1
              intermed_3c(:)=field(i,:,mb,ma)
              prim_int_2cob_field(:,:,i,mb,ma)=unpack(intermed_3c,cutoff,zero) 
           enddo
        enddo
     end do
     deallocate(field,stat=alloc_stat(11))
     if(alloc_stat(11)/=0) call error_handler &
          ("LL_CALCULATE: deallocation (10) failed")
     deallocate(intermed_3c,stat=alloc_stat(11))
     if(alloc_stat(11)/=0) call error_handler &
          ("LL_CALCULATE: deallocation (10a) failed")
  end if
#endif

  deallocate(this,stat=memstat)
  ASSERT(memstat==0)
  if(split3c)then
!!$     deallocate(this5d,stat=memstat)
!!$     ASSERT(memstat==0)
     deallocate(this6d,stat=memstat)
     ASSERT(memstat==0)
  endif

  deallocate(clmamb,clmamb_scalar,fact6,fact7,fact2,fact1,fact8,fact5,rcsabc,&
       fact0,fact4, overlap,gamma_arg,aexp_arr,bexp_arr, &
       diff_arr,diff_arr0,cutoff,stat=alloc_stat(5)) ! facts rcsabc
  if(alloc_stat(5)/=0) call error_handler &
       ("LL_CALCULATE: deallocation (10) failed")  
	alloc_stat(4)=0 !cutoff


999 CONTINUE ! cleanup and exit:

  if(opt_nested2) then
       counter=0
       do i_l=0,lb
          do mb=1,2*i_l+1
             counter=counter+1
          deallocate (nested2_summands(counter)%summands,stat=alloc_stat(2))
	  if(alloc_stat(2).ne.0) &
             call error_handler("2 deallocate nested2_summands%summands")
          end do
       end do
       deallocate(nested2_summands,stat=alloc_stat(1))
       if (alloc_stat(1).ne.0) call error_handler("1 deallocate nested2_summands failed")

  if( allocated(nested2_fac1) )then
  deallocate(nested2_fac1,nested2_fac2,nested2_fac12,stat=alloc_stat(12))
  if (alloc_stat(12).ne.0) call error_handler("12 deallocate nested2_facs failed")
  endif
  endif


  !write(*,'(A9,F17.10)') 'Result:',sum(prim_int_2cob_nuc)

contains

  !**************************************************************

  subroutine unpack_many_3c(this,i,offset,form)
    implicit none
    real(r8_kind), intent(in)    :: this(:,:,:)
    integer(i4_kind), intent(in) :: i,offset,form
    ! *** end of iterface ***

    integer(i4_kind) :: ma,mb, pos

    pos = OFF_STRIDE*(i-1)+offset
    ASSERT(pos<=size(many_3c,3))
    ! remap ints to the appropriate data structure
    select case(form)
    case (AxB)
       ASSERT(size(this,2)==2*la+1)
       ASSERT(size(this,3)==2*lb+1)
       do mb=1,2*lb+1
          do ma=1,2*la+1 
             many_3c(:,:,pos,mb,ma)=&
                  unpack(this(:,ma,mb),cutoff,zero)
          enddo
       end do
    case (BxA)
       ASSERT(size(this,2)==2*lb+1)
       ASSERT(size(this,3)==2*la+1)
       do mb=1,2*lb+1
          do ma=1,2*la+1 
             many_3c(:,:,pos,mb,ma)=&
                  unpack(this(:,mb,ma),cutoff,zero)
          enddo
       end do
    case default
       ABORT('no such case')
    end select
  end subroutine unpack_many_3c

  !**************************************************************

  subroutine  precalculate_and_nuc(nuc)
    ! precalculate prod_arr and calculate nuclear attraction
    real(r8_kind), intent(out) :: nuc(:,:,:) ! (num,2*la+1,2*lb+1)
    ! *** end of interface ***

    nuc = zero

    equal_atoms_nuc: do j=1,n_equals 
       call integral_interrupt_2cob3c()

       xc=unique_atoms(i)%position(:,j)
       yl_arr(:,:,j)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
            spread(xc,1,num))
       gamma_arg2(:,j)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
            (gamma_arg(:,3)-xc(3))**2)*fact0
       gamma_help(:,1:1+la+lb)=gamma(1+la+lb,gamma_arg2(:,j))
       fact8=2.0_r8_kind*sqrt(fact0/pi)
       yl_arr2(:,:,j)=yl_arr(:,:,j)/fact10(:,1:((ly_max+1)**2))

       counter=1
       do i_l=0,lb
          do mb=1,2*i_l+1
             diff_arr(:,:,counter)=spread((aexp_arr/fact0)**i_l,2,&
                  (la+1)**2)*diff_rule(yl_arr2(:,:,j)&
                  ,1,(la+1)**2,counter)
             counter=counter+1
          end do
       end do

    if(opt_nested2) then
       do i_l=0,lb
          do mb=1,2*i_l+1
             prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=opt_prod_rule_nested2&
                  (yl_arr(:,:,j),&
                  diff_arr(:,:,:),overlap(:,:,:),la,i_l,nested2_l1_max,nested2_l3_max, &
                  nested2_summands(i_l**2+mb)%summands,&
                  nested2_fac12)   
          end do
       end do 
     else
       do i_l=0,lb
          do mb=1,2*i_l+1
             prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=prod_rule_nested2&
                  (yl_arr(:,:,j),&
                  diff_arr(:,:,:),overlap(:,:,:),la,i_l,mb,&
                  (-2.0_r8_kind*aexp_arr),(-2.0_r8_kind*bexp_arr))   
          end do
       end do 
     endif

       ! nuclear attraction
       if (integralpar_2cob_nuc) then 
          do mb=1,2*lb+1
             do i_l =0,lb+la

		if(zc.ne.zero.and.integralpar_relativistic) then               
                   nuc_pseudo(:,:,mb)=nuc_pseudo(:,:,mb)+(z-zc)*prod_arr(:,la**2+1:&
                     (la+1)**2,lb**2+mb,i_l,j)*spread&
                     (fact8*gamma_help(:,i_l+1)&
                     ,2,2*la+1)
                else
                nuc(:,:,mb)=nuc(:,:,mb)+(z-zc)*prod_arr(:,la**2+1:&
                     (la+1)**2,lb**2+mb,i_l,j)*spread&
                     (fact8*gamma_help(:,i_l+1)&
                     ,2,2*la+1)
             end if
             
                
             end do
          end do
       end if

    end do equal_atoms_nuc
  end subroutine precalculate_and_nuc

  !**************************************************************

  subroutine s_coulomb(cexps,coul)
    ! s-type coulomb fit integrals
    real(kind=r8_kind),dimension(:),         intent(in)  :: cexps
    real(kind=r8_kind),dimension(:,:,:,:,:,:), intent(out) :: coul
    ! *** end of interface ***

    real(kind=r8_kind) :: help_vec(num)
    integer(i4_kind)   :: ncexps
    integer(i4_kind)   :: n_cf, n_if, i_ind, ic
    real(r8_kind)      :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if = unique_atoms(i)%symadapt_partner(i_ir,0)%n_independent_fcts

    ind_: do i_ind = 1, n_if

       n_cf    =  unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%n_fcts
       coeff   => unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%c
       eq_atom => unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%I_equal_atom

       sc = 0.0_r8_kind
       do ic = 1, n_cf
          sc(eq_atom(ic)) = sc(eq_atom(ic)) + coeff(ic)
       end do


    ncexps = size(cexps)

       equal_atoms_s_ch: do j=1,n_equals
          call integral_interrupt_2cob3c()
          do k=1,ncexps ! loop over fitexponents
             ! precalculation of two factors
             rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
             fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))
             gamma_help(:,1:1+la+lb)=gamma(1+la+lb,gamma_arg2(:,j)*rcsabc(:))
             do i_l=0,la+lb
                help_vec=fact8*gamma_help(:,i_l+1)
                do mb=1,2*lb+1
                   do ma=1,2*la+1            
                      coul(:,k,i_ind,mb,ma,i_pa)=coul(:,k,i_ind,mb,ma,i_pa)+&
                           ( prod_arr(:,la**2+ma,lb**2+mb,i_l,j) * & 
                           help_vec ) * NORM * sc(j)
                   enddo! loop over ma
                end do! mb
                fact8=fact8*rcsabc(:)
             end do!i_l
          end do! loop over fitfunctions
       end do equal_atoms_s_ch
    end do ind_
  end subroutine s_coulomb

  !**************************************************************

  subroutine s_xc()
    ! s-type xc fit integrals
    cexps => unique_atoms(i)%l_xc(0)%exponents(:)
    ncexps = unique_atoms(i)%l_xc(0)%n_exponents
    equal_atoms_s_xc: do j=1,n_equals           
       call integral_interrupt_2cob3c()
       do k=1,ncexps ! loop over fitexponents
          ! precalculation of two factors
          rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
          fact8=sqrt(fact0/(fact0+cexps(k)))&
               *(fact0/(fact0+cexps(k)))*exp(-rcsabc*gamma_arg2(:,j))
          ! now the acutal calculation starts
          do mb=1,2*lb+1
             do i_l=0,la+lb
                xc_int%l(0)%m(:,k,1,mb,:)=xc_int%l(0)%m&
                     (:,k,1,mb,:)+&
                     prod_arr(:,la**2+1:(la+1)**2,lb**2+mb,i_l,j)*&
                     spread(fact8*(rcsabc(:))**i_l,2,2*la+1)     
             end do! loop over l
          end do
       end do! loop over fitfunctions
    end do equal_atoms_s_xc
  end subroutine s_xc

  !**************************************************************

  subroutine r2_coulomb()
    real(kind=r8_kind),pointer,dimension(:,:,:,:,:,:) :: pointer_coul
    real(kind=r8_kind) :: help_vec(num),help_vec2(num)
    ! r2-type coulomb fit integrals

    integer(i4_kind)   :: n_cf, n_if, i_ind, ic
    real(r8_kind)      :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if = unique_atoms(i)%symadapt_partner(i_ir,0)%n_independent_fcts
    
    ind_: do i_ind = 1, n_if
       n_cf    =  unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%n_fcts
       coeff   => unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%c
       eq_atom => unique_atoms(i)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%I_equal_atom

       sc = 0.0_r8_kind
       do ic = 1, n_cf
          sc(eq_atom(ic)) = sc(eq_atom(ic)) + coeff(ic)
       end do

    ncexps = unique_atoms(i)%r2_ch%n_exponents
    cexps => unique_atoms(i)%r2_ch%exponents(:)
    pointer_coul=>coul_int(i_ir)%l(-1)%m

       equal_atoms_r2_ch: do j=1,n_equals
          call integral_interrupt_2cob3c()
          do k=1,ncexps ! loop over fitexponents
             ! precalculation of two factors
             rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
             fact8=2.0_r8_kind*pi/(cexps(k)**2)*sqrt(fact0/(fact0+cexps(k)))&
                  *(fact0/(fact0+cexps(k)))
             fact7=(fact0+1.5_r8_kind*cexps(k))/fact0
             gamma_help(:,1:2+la+lb)=gamma(2+la+lb,gamma_arg2(:,j)*&
                  rcsabc(:))
             help_vec2=gamma_arg2(:,j)*rcsabc(:)
             do i_l=0,la+lb
                help_vec=fact8*((fact7-i_l)&
                     *gamma_help(:,i_l+1)+help_vec2*gamma_help(:,i_l+2))
                do mb=1,2*lb+1
                   do ma=1,2*la+1         
                      pointer_coul(:,k,i_ind,mb,ma,i_pa) = &
                           pointer_coul(:,k,i_ind,mb,ma,i_pa)+&
                           ( prod_arr(:,la**2+ma,lb**2+mb,i_l,j)* &
                           help_vec ) * NORM * sc(j)
                   enddo! loop over ma
                end do! loop over mb
                fact8=fact8*rcsabc(:)
             end do!loop over l
          end do! loop over fitfunctions
       end do equal_atoms_r2_ch
    end do ind_
  end subroutine r2_coulomb

  !**************************************************************

  subroutine r2_xc()
    ! r2-type xc fit integrals
    ncexps = unique_atoms(i)%r2_xc%n_exponents
    cexps => unique_atoms(i)%r2_xc%exponents(:)
    equal_atoms_r2_xc: do j=1,n_equals
       call integral_interrupt_2cob3c()
       do k=1,ncexps ! loop over fitexponents
          ! precalculation of two factors
          rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
          fact8=sqrt(fact0/(fact0+cexps(k)))*fact0/(fact0+cexps(k))&
               /(fact0+cexps(k))*exp(-gamma_arg2(:,j)*rcsabc(:))
          do mb=1,2*lb+1
             do i_l=0,la+lb            
                xc_int%l(-1)%m(:,k,1,mb,:)=xc_int%l(-1)%m(:,k,1,mb,:)+&
                     prod_arr(:,la**2+1:(la+1)**2,lb**2+mb,i_l,j)*&
                     spread(fact8*rcsabc(:)**i_l*&
                     (1.5_r8_kind+fact0/cexps(k)*&
                     (gamma_arg2(:,j)*rcsabc(:)-i_l)),2,2*la+1)
             enddo! loop over l
          end do! loop over mb
       end do! loop over fitfunctions
    end do equal_atoms_r2_xc
  end subroutine r2_xc

  !**************************************************************

  subroutine l_fit_symmetry_adapt ! eq(3.63) M. Staufer
    ! symmetry adaption for l-type charge and xc fit integrals
    ! for a given l i_l
    ! results are stored in intermediate array
    intermediate=0.0_r8_kind
    independents: do i_ind=1,n_independent_fcts
       n_contributing_fcts = &
            unique_atoms(i)%symadapt_partner(1,i_l)%&
            symadapt(i_ind,1)%n_fcts
       coeff => unique_atoms(i)%symadapt_partner(1,i_l)%&
            symadapt(i_ind,1)%c
       magn => unique_atoms(i)%symadapt_partner(1,i_l)%&
            symadapt(i_ind,1)%m
       eq_atom => unique_atoms(i)%symadapt_partner(1,i_l)%&
            symadapt(i_ind,1)%&
            I_equal_atom 
       contributing: do i_cnt=1,n_contributing_fcts
          counter=1
          do i_lb=0,lb
             do mb=1,2*i_lb+1
                diff_arr(:,:,counter)=&
                     spread(&
                     (& 
                     aexp_arr/ fact0  )**i_l   &   !   (a/(a+b))^i_l
                     *(bexp_arr/aexp_arr)**i_lb, & ! * (b/a)^i_lb
                     2,        &     
                     (la+1)**2 &
                     ) &
                     * diff_rule_nested(&          ! res: p(lm_grad) = sh_grad(counter) (sh_grad(lm_grad) f(i_l^2+magn(i_cnt)))
                     yl_arr2(:,:,eq_atom(i_cnt)),& ! C(num,(ly_max+1)^2,eq_atom(i_cnt))
                     1,&
                     (la+1)**2,&
                     counter,&
                     i_l**2+magn(i_cnt))           ! eq(3.61) M.Staufer, without X
                counter=counter+1
             end do
          end do
          do mb=1,2*lb+1
             intermediate(:,:,mb,i_ind,eq_atom(i_cnt),:)=&
                  intermediate(:,:,mb,i_ind,eq_atom(i_cnt),:)&
                  + coeff(i_cnt)&
                  * prod_rule_nested3(diff_arr(:,:,:),prod_arr(:,:,:,:,eq_atom(i_cnt)),la,lb,mb,la+lb,i_l) ! eq (3.63) M.Staufer
          end do
       end do contributing
    end do independents
  end subroutine l_fit_symmetry_adapt
  !**************************************************************
  
  subroutine l_fit_symmetry_adapt_v2(i,i_l,i_ir,i_pa,ZL)
    implicit none
    integer(i4_kind), intent(in) :: i, i_l, i_ir, i_pa 
    real(r8_kind), intent(out)   :: ZL(:,:,:,:,:,0:)
    !*** end of interface ***
    integer(i4_kind)         :: n_if, n_cf, i_cnt, counter, i_lb, mb 
    real(r8_kind),pointer    :: coeff(:)
    integer(i4_kind),pointer :: magn(:), eq_atom(:)
    real(r8_kind)            :: NORM
    
    ! eq(3.63) from M. Staufer rebuild for symmetry
    ! symmetry adaption for l-type charge and xc fit integrals
    ! for a given l -> i_l, i_ir and i_pa,
    ! results are stored in intermediate array ZL(num,ma,mb,i_ind,n_equal,0:la+lb)

!!$    NORM  =  SQRT(REAL(n_equals,r8_kind))    
    NORM = 1.0D0
    ZL(:,:,:,:,:,:)=0.0D0
    n_if = unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

    independents: do i_ind=1,n_if

       ! i_ir and i_pa was added because of the symmetrization of fit_fct

       n_cf = &
                  unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%n_fcts

       coeff   => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%c

       magn    => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%m

       eq_atom => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%I_equal_atom 

       contributing: do i_cnt=1,n_cf
          counter=1
          do i_lb=0,lb
             do mb=1,2*i_lb+1
                diff_arr(:,:,counter)=&
                     spread(&
                     (& 
                     aexp_arr/ fact0  )**i_l   &   !   (a/(a+b))^i_l
                     *(bexp_arr/aexp_arr)**i_lb, & ! * (b/a)^i_lb
                     2,        &     
                     (la+1)**2 &
                     ) &
                     * diff_rule_nested(&          ! res: p(lm_grad) = sh_grad(counter) (sh_grad(lm_grad) f(i_l^2+magn(i_cnt)))
                     yl_arr2(:,:,eq_atom(i_cnt)),& ! C(num,(ly_max+1)^2,eq_atom(i_cnt))
                     1,&
                     (la+1)**2,&
                     counter,&
                     i_l**2+magn(i_cnt))           ! eq(3.61) M.Staufer, without X
                counter=counter+1
             end do
          end do
          do mb=1,2*lb+1

             ZL(:,:,mb,i_ind,eq_atom(i_cnt),:)=&  ! i_ir and i_pa was added for symmetrization
                  ZL(:,:,mb,i_ind,eq_atom(i_cnt),:)&
                  + NORM * coeff(i_cnt)&                  ! Here you can see the influence of symmetry (C from Eq. 3.63 MS)
                  * prod_rule_nested3(diff_arr(:,:,:),prod_arr(:,:,:,:,eq_atom(i_cnt)),la,lb,mb,la+lb,i_l) ! eq (3.63) M.Staufer
          end do
       end do contributing
    end do independents
  end subroutine l_fit_symmetry_adapt_v2

  !**************************************************************

  subroutine l_coulomb()  ! eq(3.64) M.Staufer 
    ! l-type coulomb fit integrals for a given i_l
    real(kind=r8_kind) :: help_vec(num),help_vec2(num)
    real(kind=r8_kind),pointer,dimension(:,:,:,:,:,:) :: pointer_coul
    ncexps =  unique_atoms(i)%l_ch(i_l)%n_exponents
    cexps => unique_atoms(i)%l_ch(i_l)%exponents(:)
    
    pointer_coul=>coul_int(i_ir)%l(i_l)%m

    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       equals_ch:  do j=1,n_equals
          rcsabc=cexps(k)/(fact0+cexps(k))                       ! gm/a+b+gm
          gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2&  ! I(lc+l)
               (:,j)*rcsabc(:))
          fact8=2.0_r8_kind*pi/cexps(k)*sqrt(fact0/&             ! Pc = 2Pi/gm*(a+b/a+b+gm)^(1/2)
               (fact0+cexps(k)))
          counter=1
          help_vec=(2.0_r8_kind*fact0*rcsabc)**i_l*fact8         ! Pc*(...)*I(lc+l)*(...) from Eq. 3.64 M.Staufer
          do l=0,la+lb
             help_vec2=help_vec*gamma_help(:,l+1+i_l)
             do i_ind=1,n_independent_fcts
                do mb=1,2*lb+1
                   do ma=1,2*la+1
                      pointer_coul(:,k,i_ind,mb,ma,i_pa)=&       ! Pc*(...)*I(lc+l)*(...)*Z Eq. 3.64 as it M.Staufer
                           pointer_coul(:,k,i_ind,mb,ma,i_pa)+&
                           intermediate(:,ma,mb,i_ind,j,l)*&
                           help_vec2
                   end do! loop over ma
                end do! loop over mb
             end do! loop over i_ind
             help_vec=help_vec*rcsabc(:)
          end do! loop over l
       end do equals_ch
    end do! loop over fitfunctions
  end subroutine l_coulomb

  !**************************************************************

  subroutine l_xc()
    ! l-type xc fit integrals for a given i_l
    cexps => unique_atoms(i)%l_xc(i_l)%exponents(:)
    ncexps = unique_atoms(i)%l_xc(i_l)%n_exponents
    allocate(xc_int%l(i_l)%m(num,ncexps,n_independent_fcts,&
         2*lb+1,2*la+1),&
         stat=alloc_stat(15))
    if(alloc_stat(15)/=0) call error_handler &
         ("LL_CALCULATE: allocation (9) failed")  
	alloc_stat(15)=1
    xc_int%l(i_l)%m=0.0_r8_kind
    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       equals_xc:  do j=1,n_equals
          rcsabc=cexps(k)/(fact0+cexps(k))
          fact8=sqrt(fact0/(fact0+cexps(k)))*fact0/&
               (fact0+cexps(k))
          gamma_help(:,1)=exp(-gamma_arg2(:,j)*rcsabc(:))
          do i_ind=1,n_independent_fcts
             do mb=1,2*lb+1
                do l=0,la+lb
                   xc_int%l(i_l)%m(:,k,i_ind,mb,:)=xc_int%&
                        l(i_l)%m(:,k,i_ind,mb,:)+&
                        intermediate(:,:,mb,i_ind,j,l)*spread(&
                        (2.0_r8_kind*fact0*rcsabc)**i_l*&
                        (rcsabc(:))**l*fact8*&
                        gamma_help(:,1),2,2*la+1)
                end do! loop over l
             end do! loop over mb
          end do! loop over i_ind   
       end do equals_xc! loop over equal atoms
    end do! loop over k 
  end subroutine l_xc


  subroutine add_pointcharges
    ! add contribution of point charges to nuclear attraction
    ly_max=max(la,lb)
    max_order = 3+la+lb
    unique_charge_loop: do i=1,pointcharge_N+n_timps  ! loop over third center
       if(i<=n_timps) then
          z= unique_timps(i)%z -  unique_timps(i)%zc   ! charge 
          n_equals=unique_timps(i)%n_equal_atoms          
       else
          cycle unique_charge_loop ! pointcharges go to SHGI !!!!!!!!!!!!!AS
          z= pointcharge_array(i-n_timps)%z                  ! charge 
          n_equals=pointcharge_array(i-n_timps)%n_equal_charges
       end if
       allocate ( gamma_help(num,max_order),&
            gamma_arg2(num,n_equals),&
            yl_arr(num,(ly_max+1)**2,n_equals),&
            fact10(num,(ly_max+1)**2),&
            yl_arr2(num,(ly_max+1)**2,n_equals),&
            prod_arr(num,(la+1)**2,(lb+1)**2,0:la+lb,n_equals),&
            stat=alloc_stat(27))
       if (alloc_stat(27)/=0) call error_handler &
            ("LL_CACLULATE : allocation (10) failed")
       alloc_stat(27)=1
       prod_arr=0.0_r8_kind
       ! do a precalculation of a factor needed for the
       ! product rule
       counter=1
       fact4=1.0_r8_kind
       do i_l=0,ly_max
          do ma=1,2*i_l+1
             fact10(:,counter)=fact4
             counter=counter+1
          enddo
          fact4=fact4*aexp_arr/(fact0)
       enddo

       equal_charge_nuc: do j=1,n_equals
          if(i<=n_timps) then
             xc=unique_timps(i)%position(:,j)
          else
             xc=pointcharge_array(i-n_timps)%position(:,j)
          end if
          yl_arr(:,:,j)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
               spread(xc,1,num))
          gamma_arg2(:,j)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
               (gamma_arg(:,3)-xc(3))**2)*fact0
          gamma_help(:,1:1+la+lb)=gamma(1+la+lb,gamma_arg2(:,j))
          fact8=2.0_r8_kind*sqrt(fact0/pi)
          counter=1
          yl_arr2(:,:,j)=yl_arr(:,:,j)/fact10(:,1:((ly_max+1)**2))

          do i_l=0,lb
             do mb=1,2*i_l+1
                diff_arr(:,:,counter)=spread((aexp_arr/fact0)**i_l,2,&
                     (la+1)**2)*diff_rule(yl_arr2(:,:,j)&
                     ,1,(la+1)**2,counter)
                counter=counter+1
             end do
          end do

 if(opt_nested2) then
          do i_l=0,lb
             do mb=1,2*i_l+1

             prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=opt_prod_rule_nested2&
                  (yl_arr(:,:,j),&
                  diff_arr(:,:,:),overlap(:,:,:),la,i_l,nested2_l1_max,nested2_l3_max, &
                  nested2_summands(i_l**2+mb)%summands,&
                  nested2_fac12)   
             end do
          end do
 else
          do i_l=0,lb
             do mb=1,2*i_l+1
                prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=prod_rule_nested2&
                     (yl_arr(:,:,j),&
                     diff_arr(:,:,:),overlap(:,:,:),la,i_l,mb,&
                     (-2.0_r8_kind*aexp_arr),(-2.0_r8_kind*bexp_arr))   
             end do
          end do
 endif

          do mb=1,2*lb+1
             do i_l =0,lb+la
                if(pseudopot_present.and.integralpar_relativistic) then	
                   nuc_pc_timps(:,:,mb)=nuc_pc_timps(:,:,mb) &
                        +z*prod_arr(:,la**2+1:&
                        (la+1)**2,lb**2+mb,i_l,j)*spread&
                        (fact8*gamma_help(:,i_l+1)&
                        ,2,2*la+1)
                else
                   nuc(:,:,mb)=nuc(:,:,mb)+z*prod_arr(:,la**2+1:&
                        (la+1)**2,lb**2+mb,i_l,j)*spread&
                        (fact8*gamma_help(:,i_l+1)&
                        ,2,2*la+1)
                endif
             end do
          end do

       end do equal_charge_nuc

       deallocate(yl_arr,yl_arr2,gamma_help,gamma_arg2,fact10,&
            prod_arr,stat=alloc_stat(27))
       if (alloc_stat(27)/=0) call error_handler &
            ("LL_CACLULATE : deallocation (8) failed")
    end do unique_charge_loop

  end subroutine add_pointcharges

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !**************************************************************

  !********************************************************************************
  subroutine calc_potential

    ! integrals of potential
    ly_max=max(la,lb)
    max_order = 3+la+lb

    do i=1,N_points
       n_equals=point_in_space(i)%n_equal_points
       allocate ( gamma_help(num,max_order),&
            gamma_arg2(num,n_equals),&
            yl_arr(num,(ly_max+1)**2,n_equals),&
            fact10(num,(ly_max+1)**2),&
            yl_arr2(num,(ly_max+1)**2,n_equals),&
            prod_arr(num,(la+1)**2,(lb+1)**2,0:la+lb,n_equals),&
            stat=alloc_stat(29))
       if (alloc_stat(29)/=0) call error_handler &
            ("LL_CACLULATE : allocation (10) failed")
	alloc_stat(29)=1
       prod_arr=0.0_r8_kind
       ! do a precalculation of a factor needed for the
       ! product rule
       counter=1
       fact4=1.0_r8_kind
       do i_l=0,ly_max
          do ma=1,2*i_l+1
             fact10(:,counter)=fact4
             counter=counter+1
          enddo
          fact4=fact4*aexp_arr/(fact0)
       enddo
       equal_points: do j=1,n_equals
          call integral_interrupt_2cob3c()
          xc=point_in_space(i)%position(:,j)
          yl_arr(:,:,j)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
               spread(xc,1,num))
          gamma_arg2(:,j)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
               (gamma_arg(:,3)-xc(3))**2)*fact0
          gamma_help(:,1:1+la+lb)=gamma(1+la+lb,gamma_arg2(:,j))
          fact8=2.0_r8_kind*sqrt(fact0/pi)
          yl_arr2(:,:,j)=yl_arr(:,:,j)/fact10(:,1:((ly_max+1)**2))

          counter=1
          do i_l=0,lb
             do mb=1,2*i_l+1
                diff_arr(:,:,counter)=spread((aexp_arr/fact0)**i_l,2,&
                     (la+1)**2)*diff_rule(yl_arr2(:,:,j)&
                     ,1,(la+1)**2,counter)
                counter=counter+1
             end do
          end do

 if(opt_nested2) then
          do i_l=0,lb
             do mb=1,2*i_l+1

             prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=opt_prod_rule_nested2&
                  (yl_arr(:,:,j),&
                  diff_arr(:,:,:),overlap(:,:,:),la,i_l,nested2_l1_max,nested2_l3_max, &
                  nested2_summands(i_l**2+mb)%summands,&
                  nested2_fac12)   
             end do
          end do
 else
          do i_l=0,lb
             do mb=1,2*i_l+1
                prod_arr(:,:,i_l**2+mb,0:la+i_l,j)=prod_rule_nested2&
                     (yl_arr(:,:,j),&
                     diff_arr(:,:,:),overlap(:,:,:),la,i_l,mb,&
                     (-2.0_r8_kind*aexp_arr),(-2.0_r8_kind*bexp_arr))   
             end do
          end do
 endif

          do mb=1,2*lb+1
             do i_l =0,lb+la
                potential(i,:,:,mb)=potential(i,:,:,mb)+prod_arr(:,la**2+1:&
                     (la+1)**2,lb**2+mb,i_l,j)*spread&
                     (fact8*gamma_help(:,i_l+1)&
                     ,2,2*la+1)
             end do
          end do

       end do equal_points

       deallocate(yl_arr,yl_arr2,gamma_help,gamma_arg2,fact10,&
            prod_arr,stat=alloc_stat(29))
       if (alloc_stat(29)/=0) call error_handler &
            ("LL_CACLULATE : deallocation (10) failed")
    enddo

  end subroutine calc_potential
  !********************************************************************************

  !***********************************************************
#if 0
  subroutine calc_field

    integer(kind=i4_kind) :: i,j,f_dim,i_l,counter,l,m,l2pm,k_gr,k,mb,i_la,i_lb,i_grad,alloc_stat1
    integer(kind=i4_kind) :: nm_la,nm_lb,lasq,lbsq,la_index,lb_index,xyz_map(3),n_equals
    real(kind=r8_kind),allocatable :: poten_grad(:,:,:,:),grad_poten(:,:,:,:)
    real(kind=r8_kind),allocatable :: help_arr(:,:), help_arr0(:,:), aqapb(:,:)
    real(kind=r8_kind),allocatable :: help_arr_gr(:,:,:,:)
    real(kind=r8_kind),pointer :: rotmat(:,:)
    logical :: do_rotation

    nm_la=2*la+1
    nm_lb=2*lb+1
    lasq=la**2
    lbsq=lb**2

    xyz_map(1)=3_i4_kind
    xyz_map(2)=4_i4_kind
    xyz_map(3)=2_i4_kind
    ly_max=max(la,lb)
    max_order=la+lb+4

    allocate (grad_poten(num,3,nm_lb,nm_la),stat=alloc_stat1)
    if (alloc_stat1/=0) call error_handler &
         ("LL_CALCULATE : allocation grad_poten failed")

    allocate ( gamma_help(num,max_order), gamma_arg2(num,1),  & 
         yl_arr(num,(ly_max+1)**2,1),yl_arr_xyz(num,(ly_max+1)**2,3),&
         help_arr(num,(ly_max+1)**2),help_arr0(num,(la+1)**2),&
         fact10(num,(ly_max+1)**2),yl_arr2(num,(ly_max+1)**2,1),&
         aqapb(num,0:ly_max),diff_arr_xyz(num,(la+1)**2,(lb+1)**2,3), &
         help_vec(num),help_arr_gr(num,3,nm_lb,nm_la), &
         prod_arr_gr(num,3,(la+1)**2,(lb+1)**2,0:la+lb+1), &
         help_mat(num,(la+1)**2,(lb+1)**2),prod_arr_gr_vec(num*6*(la+lb+2),(la+1)**2,(lb+1)**2), &
         stat=alloc_stat1)
    if (alloc_stat1/=0) call error_handler &
         ("LL_CALCULATE : allocation (1) failed")

    aqapb(:,0)=1.0_r8_kind
    aqapb(:,1)=aexp_arr/fact0

    Nsp: do i=1,N_surface_points
       n_equals=surface_points(i)%N_equal_points

       f_dim=surf_points_grad_index(i+1)-surf_points_grad_index(i)

       allocate ( poten_grad(num,nm_lb,nm_la,f_dim),stat=alloc_stat1)
       if (alloc_stat1/=0) call error_handler &
         ("LL_CALCULATE : allocation poten_grad failed")
       poten_grad=0.0_r8_kind

       do i_l=2,ly_max
          aqapb(:,i_l)=aqapb(:,i_l-1)*aqapb(:,1)
       end do

       ! do a precalculation of a factor needed for the diff rule
       counter=1
       fact4=1.0_r8_kind
       do i_l=0,ly_max
          do ma=1,2*i_l+1
             fact10(:,counter)=fact4
             counter=counter+1
          enddo
          fact4=fact4*aqapb(:,1)
       enddo

       n_equals_j: do j=1,n_equals

          if(f_dim == 3) then
             if(sum((surf_points_grad_info(i)%m(:,:,j)-unity_matrix)**2)<&
                  1.0e-7_r8_kind) then
                do_rotation=.false.
             else
                do_rotation=.true.
                rotmat=>surf_points_grad_info(i)%m(:,:,j)
             endif
          else
             do_rotation=.true.
             rotmat=>surf_points_grad_info(i)%m(:,:,j)
          endif

          xc=surface_points(i)%position(:,j)

          yl_arr_xyz=0.0_r8_kind

          call integral_interrupt_2cob3c()

          yl_arr(:,:,1)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
               spread(xc,1,num))
          gamma_arg2(:,1)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
               (gamma_arg(:,3)-xc(3))**2)*fact0
          gamma_help(:,1:max_order)=gamma(max_order,gamma_arg2(:,1))
          fact8=2.0_r8_kind*sqrt(fact0/pi)
          counter=1
          yl_arr2(:,:,1)=yl_arr(:,:,1)/fact10(:,1:((ly_max+1)**2))
          ! now calculation of derivatives of yl_arr with respect to -c
          do l=0,ly_max
             do m=1,2*l+1
                l2pm=l*l+m
                do k_gr=1,3
                   do k=1,solhrules_differential(xyz_map(k_gr),l2pm)%n_summands
                      yl_arr_xyz(:,l2pm,k_gr)=yl_arr_xyz(:,l2pm,k_gr)+&
                           solhrules_differential(xyz_map(k_gr),l2pm)%coef(k)*&
                           yl_arr(:,solhrules_differential(xyz_map(k_gr),l2pm)%lm_sh(k),1)
                   end do
                end do
             end do
          end do!loop over l

          do i_l=0,lb
             help_arr0=spread(aqapb(:,i_l),2,(la+1)**2)
             do mb=1,2*i_l+1
                diff_arr(:,:,counter)=&
                     help_arr0*diff_rule(yl_arr2(:,:,1),1,(la+1)**2,counter)
                ! now let`s differentiate diff_arr
                do k_gr=1,3
                   diff_arr_xyz(:,:,counter,k_gr)=&
                        help_arr0*diff_rule(yl_arr_xyz(:,:,k_gr)/fact10,1,(la+1)**2,counter)
                end do
                counter=counter+1
             end do
          end do
          do mb=1,(lb+1)**2
             do ma=1,(la+1)**2
                help_mat(:,ma,mb)= overlap(:,ma,mb)*aqapb(:,1)
             enddo
          enddo
          help_vec=-2.0_r8_kind*fact0

          call calc_prod_arr()
          
          help_arr=0.0_r8_kind
          help_arr_gr=0.0_r8_kind
          do i_l =0,lb+la+1
             help_vec=fact8*gamma_help(:,i_l+1)
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   la_index=lasq+i_la
                   lb_index=lbsq+i_lb
                   do k_gr=1,3
                      help_arr_gr(:,k_gr,i_lb,i_la)=help_arr_gr(:,k_gr,i_lb,i_la)&
                           +prod_arr_gr(:,k_gr,la_index,lb_index,i_l)*help_vec
                   end do! loop over k_gr
                enddo! loop over i_lb
             end do! loop over i_la
          end do! loop over i_l

          if(do_rotation) then
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   ! make gradient totalsymmetric before adding
                   do i_grad=1,f_dim 
                      poten_grad(:,i_lb,i_la,i_grad)=poten_grad(:,i_lb,i_la,i_grad)-&
                           rotmat(i_grad,1)*help_arr_gr(:,1,i_lb,i_la)-&
                           rotmat(i_grad,2)*help_arr_gr(:,2,i_lb,i_la)-&
                           rotmat(i_grad,3)*help_arr_gr(:,3,i_lb,i_la)
                   enddo
                end do
             enddo
          else
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   do i_grad=1,3
                      poten_grad(:,i_lb,i_la,i_grad)=poten_grad(:,i_lb,i_la,i_grad)-&
                           help_arr_gr(:,i_grad,i_lb,i_la)
                   enddo
                end do
             enddo
          end if
       enddo n_equals_j

       if(calc_normal) then
          grad_poten=0.0_r8_kind
          rotmat=>surf_points_grad_info(i)%m(:,:,1)
          do j=1,f_dim
             do k=1,3
                grad_poten(:,k,:,:)= grad_poten(:,k,:,:)+rotmat(j,k)*poten_grad(:,:,:,j)
             enddo
          enddo

          do k=1,3
             field(i,:,:,:)=field(i,:,:,:)+grad_poten(:,k,:,:)*surface_points(i)%out_normal(k)
          enddo
       else
          k=surf_points_grad_index(i)
          do j=1,f_dim
             field(k,:,:,:)=poten_grad(:,:,:,j)
             k=k+1
          end do
       end if

       deallocate (poten_grad,stat=alloc_stat1)
       if (alloc_stat1/=0) call error_handler &
            ("LL_CALCULATE : deallocation poten_grad failed")

    end do Nsp

!!$       print*,'========FIELD_LL========='
!!$       print*,field(1,1,1,1)
!!$       print*,'========================='

    deallocate ( grad_poten,stat=alloc_stat1)
    if (alloc_stat1/=0) call error_handler &
         ("LL_CALCULATE : deallocation grad_poten failed")
    deallocate ( gamma_help, gamma_arg2,yl_arr,yl_arr_xyz,&
         help_arr,help_arr0,fact10,yl_arr2,help_arr_gr,  &  !!!prod_arr, &
         aqapb,diff_arr_xyz,help_vec,prod_arr_gr,help_mat,prod_arr_gr_vec, stat=alloc_stat1)
    if (alloc_stat1/=0) call error_handler &
         ("LL_CALCULATE : allocation (1) failed")

  end subroutine calc_field
#endif
  !**********************************************************

  !**********************************************************
  subroutine calc_prod_arr
    real(kind=r8_kind) :: help_arr1(3*num,(ly_max+1)**2),help_arr2(3*num,0:la+lb,(la+1)**2),&
         coef,twoa(num),twob(num)
    integer(kind=i4_kind) :: l_sum,i3_p_j3,i3_p2_j3,i_p_j1,i2_p_j2,i0_p_j1,meta_min,meta_max
    real(kind=r8_kind),dimension(num) :: help_vec1,help_vec0,help_vec2,help_vec3
    integer(kind=i4_kind) :: i1,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j3, l_j1,k_gr,i_ma,i_i,lb_index
    integer(kind=i4_kind),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(kind=r8_kind),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    real(kind=r8_kind) :: overlap_help(num*6),help_spread(num*6)   !!!!!!!overlap_xyz_help(num*3)
    
    twoa=-2.0_r8_kind*aexp_arr
    twob=-2.0_r8_kind*bexp_arr

    prod_arr_gr_vec=0.0_r8_kind
    do l1=1,(ly_max+1)**2
       do k_gr=1,3
          help_arr1((k_gr-1)*num+1:k_gr*num,l1)=yl_arr(:,l1,1)*(gamma_arg(:,k_gr)-xc(k_gr)) 
       end do
    end do

    help_vec=-2.0_r8_kind*fact0

    do i_l=0,lb
       do mb=1,2*i_l+1
          help_arr2=0.0_r8_kind
          lb_index=i_l**2+mb
          lm_in=lb_index
          lm_min=1
          lm_max=(la+1)**2
          index_p=>solhrules_product(lm_in)%lm_sh1
          index0_p=>solhrules_product(lm_in)%lm_sh2
          coef1_p=>solhrules_product(lm_in)%coef
          n_sum1=solhrules_product(lm_in)%n_summands
          do i1=lm_min,lm_max
             index1_p=>solhrules_product(i1)%lm_sh1
             index2_p=>solhrules_product(i1)%lm_sh2
             n_sum2=solhrules_product(i1)%n_summands
             coef2_p=>solhrules_product(i1)%coef
             do j_1=1,n_sum1
                i_p_j1=index_p(j_1)
                i0_p_j1=index0_p(j_1)
                l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                help_vec2=twob**l1
                do j_2=1,n_sum2
                   i2_p_j2=index2_p(j_2)
                   index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
                   index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
                   n_sum3=solhrules_product(index1_p(j_2))%n_summands
                   coef3_p=>solhrules_product(index1_p(j_2))%coef
                   overlap_help(1:num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(num+1:2*num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(2*num+1:3*num) = help_mat(:,i2_p_j2,i0_p_j1)
                   overlap_help(3*num+1:4*num) = overlap(:,i2_p_j2,i0_p_j1)
                   overlap_help(4*num+1:5*num) = overlap(:,i2_p_j2,i0_p_j1)
                   overlap_help(5*num+1:6*num) = overlap(:,i2_p_j2,i0_p_j1)
                   do j_3=1,n_sum3
                      i3_p_j3=index3_p(j_3)
                      i3_p2_j3=index3_p2(j_3)
                      l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                      l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                      l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                      if(l_j3<=l_j1) then
                         coef=coef2_p(j_2)*coef3_p(j_3)*coef1_p(j_1)
                         help_vec1=coef*help_vec2*twoa**l3
                         l_sum=l1+l3
                         meta_min = l_sum * 6 * num + 1
                         meta_max = l_sum * 6 * num + 3 * num
                         help_vec0=help_vec1*&
                              yl_arr(:,i3_p_j3,1)*diff_arr(:,i3_p2_j3,i_p_j1)
                         meta_max = (l_sum+1) * 6 * num
                         if(l_j3<l_j1) then
                            help_vec3=help_vec1*yl_arr(:,i3_p_j3,1)
                            help_spread(1:num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,1)
                            help_spread(num+1:2*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,2)
                            help_spread(2*num+1:3*num) = help_vec3*diff_arr_xyz(:,i3_p2_j3,i_p_j1,3)
!!$!FPP:VPP!!OCL NOVREC(help_spread)
                            help_spread(3*num+1:6*num) = help_spread(1:3*num)
                            prod_arr_gr_vec(meta_min:meta_max,i1,lb_index) = &
                                 prod_arr_gr_vec(meta_min:meta_max,i1,lb_index) + &
                                 help_spread * overlap_help
                         endif
                         help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)
                         help_spread(1:num) = help_vec3*yl_arr_xyz(:,i3_p_j3,1)
                         help_spread(num+1:2*num) = help_vec3*yl_arr_xyz(:,i3_p_j3,2)
                         help_spread(2*num+1:3*num) = help_vec3*yl_arr_xyz(:,i3_p_j3,3)
!!$!FPP:VPP!!OCL NOVREC(help_spread)
                         help_spread(3*num+1:6*num) = help_spread(1:3*num)
                         prod_arr_gr_vec(meta_min:meta_max,i1,lb_index) = &
                              prod_arr_gr_vec(meta_min:meta_max,i1,lb_index) + &
                              help_spread * overlap_help
                         help_vec3=help_vec1*diff_arr(:,i3_p2_j3,i_p_j1)*&
                              overlap(:,i2_p_j2,i0_p_j1)
                         help_spread(1:num) = help_vec3
                         help_spread(num+1:2*num) = help_vec3
                         help_spread(2*num+1:3*num) = help_vec3
                         help_arr2(:,l_sum,i1) = help_arr2(:,l_sum,i1) + &
                              help_spread(1:3*num) * help_arr1(:,i3_p_j3)
                      endif
                   end do! j_3
                end do! j_2
             enddo! j_1
          end do! i1=1
          help_spread(3*num+1:4*num) = help_vec
          help_spread(4*num+1:5*num) = help_vec
          help_spread(5*num+1:6*num) = help_vec
          do i_ma=1,(la+1)**2
             do i_i=1,1+la+i_l
                prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) = &
                     prod_arr_gr_vec((i_i*6+3)*num+1:(i_i+1)*6*num,i_ma,lb_index) + &
                     help_spread(3*num+1:6*num) * help_arr2(:,i_i-1,i_ma)
             enddo
          end do
       end do! mb
    end do! i_l
    call map_prod_arr_gr()
    
  end subroutine calc_prod_arr
  !**********************************************************
  
  !**********************************************************
  subroutine map_prod_arr_gr
    ! Purpose: piece of the code 
    ! maps  prod_arr_gr_vec to  prod_arr_gr 
    integer(kind=i4_kind) :: i_ma, i_mb, k_gr, l, i_min, i_max
    do i_mb=1,(lb+1)**2
       do i_ma=1,(la+1)**2
          i_min = 1
          i_max = num
          do l = 0, la+lb+1
             do k_gr=1,6
                if(k_gr >= 4) then
                   prod_arr_gr(:,k_gr-3,i_ma,i_mb,l) = &
                        prod_arr_gr_vec(i_min:i_max,i_ma,i_mb)
                endif
                i_min = i_min + num
                i_max = i_max + num
             end do
          end do
       end do
    end do
  end subroutine map_prod_arr_gr
  !**********************************************************

#if 0
   !**************************************************************!
   ! Just to print nuclear integrals !
   subroutine print_nuclear()
     use iounitadmin_module
     integer(kind=i4_kind) :: i_exp1, i_exp2, nuc_int_unit
     logical               :: file_exists

     nuc_int_unit = get_iounit()
     inquire(file="NUC_INT", exist=file_exists)
     if(file_exists) then
        open(nuc_int_unit, file="NUC_INT",position="append")
     else
        open(nuc_int_unit, file="NUC_INT")
     end if
     write(nuc_int_unit,"('********* LL ********** [',2i2,' |',2i2,' ] ********* LL ***********')")  na, la, nb, lb
     do i_exp1 = 1,naexps
        do i_exp2 = 1,nbexps
           do mb=1,2*lb+1
              do ma=1,2*la+1 
                 write(nuc_int_unit,"(5x,'exp2 = ',i3,3x,'exp1 = ',i3,3x, 'm',2i3,4x,f14.8)") &
                      i_exp2, i_exp1, mb,ma, prim_int_2cob_nuc(i_exp2,i_exp1,mb,ma)
              end do
           end do
        end do
     end do
     close(nuc_int_unit, status="keep")    
     call return_iounit(nuc_int_unit)

   end subroutine print_nuclear
   !**************************************************************!
#endif


end subroutine ll_calculate

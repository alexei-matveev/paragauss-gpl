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
subroutine ss_calculate(na,nb,imode,many_3c)
  !
  !  Purpose: calculation of all primitive 2 center orbital
  !           and 3 center integrals for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,l2,equal_atom1,equal_atom2).
  !       For three center integrals, contraction and symmetry-
  !       adaption concerning fitfunctions is also performed.
  !
  !
  !  Author: FN
  !  Date:   7/96
  !
!================================================================
! End of public interface of module
!================================================================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! letzter debug check am 17.7.96, nur stichprobenartig
  ! kontrolliert: co,h2,cuh
  ! FN
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date: 30.7.96
  ! Description: Restructuring of subroutine.
  !              Lots of missing deallocations added
  !              Steering variables from integralpar_module added
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:  27.12.96
  ! Description: Relativistic pv scalar p matrix elements added
  !
  !
  ! Modification (Please copy before editing)
  ! Author: AH
  ! Date:   22.09.98
  ! Description: pseudopotential
  !
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   11-12/99
  ! Description: integrals of electrostatic potential are added
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   03/00
  ! Description: integrals of electrostatic field are added
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  use unique_atom_module, noname=>pseudopot_present
  use gamma_module
  use type_module
  use datatype
  use solid_harmonics_module, only : solid_harmonics_calc,&
                                     solid_harmonics_scalar
  use solhrules_module
  use fitcontract_module, only: fitcontract
  !!SB: should be only for resp, but...
  use fitcontract_module, only: fitcontract_v2

  use integralpar_module
  use iounitadmin_module
  use int_data_2cob3c_module, only : &
       prim_int_2cob_ol, &
       prim_int_2cob_kin, &
       prim_int_2cob_poten, &
       prim_int_2cob_nuc, &
       prim_int_2cob_pvsp, &
       prim_int_3c_co, &
       prim_int_3c_xc, &
       prim_int_2cob_nuc_pseudo, &
       prim_int_2cob_field,&
       center1, &
       center2, &
       n_m1, &
       n_m2, &
       ua1_basis, &
       ua2_basis, &
       OFF_STRIDE, OFF_V, OFF_PVSP, OFF_VFIN
  use pointcharge_module
  use options_module, only: options_integral_expmax
  use potential_module
  use elec_static_field_module
  use operations_module, only: operations_core_density
  use symmetry_data_module, only: symmetry_data_n_irreps, &
       symmetry_data_n_partners, &
       get_totalsymmetric_irrep
  use calc3c_switches,only: old_potential,old_3c_co,old_elfield
  use shgi_cntrl, only: IPSEU
  implicit none
  integer(kind=i4_kind), intent(in) :: na,nb
  integer(kind=i8_kind), intent(in) :: imode ! control bitmask
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
  ! *** end of interface ***

  integer(kind=i4_kind) :: m1,m2

! real(kind=r8_kind),dimension(3,3),parameter :: unity_matrix=reshape&
  real(kind=r8_kind),dimension(3,3) :: unity_matrix=reshape&
       ((/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind,&
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/),(/3,3/))

  integer(kind=i4_kind) :: naexps,nbexps,ncexps
  real(kind=r8_kind),pointer     :: aexps(:),bexps(:)
  real(kind=r8_kind),pointer     :: cexps(:)
  real(kind=r8_kind),pointer     :: rotmat(:,:)  !!!!!!!!!!!!!!!!
  real(kind=r8_kind) :: z  ! charge
  real(kind=r8_kind) :: zc ! core charge
  ! constants
  real(kind=r8_kind),parameter    :: pi=3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter    :: two=2.0_r8_kind,three=3.0_r8_kind
  real(kind=r8_kind),parameter    :: four=4.0_r8_kind,six=6.0_r8_kind
  real(kind=r8_kind),parameter    :: very_small=1.0e-100_r8_kind
  real(kind=r8_kind),parameter    :: very_big=1.0e100_r8_kind
  real(kind=r8_kind),parameter    :: zero=0.0_r8_kind
!!! parameter-bug on sgi
! real(kind=r8_kind),parameter,dimension(0:8) :: dfac= &
  real(kind=r8_kind), dimension(0:8) :: dfac= &
       (/ 1.0_r8_kind, 1.0_r8_kind, 3.0_r8_kind, 15.0_r8_kind, 105.0_r8_kind, &
       945.0_r8_kind, 10395.0_r8_kind, 135135.0_r8_kind, 2027025.0_r8_kind /)

  ! mapping of exponents to one dimension and cutoff of small integrals
  logical,allocatable   :: cutoff(:,:) ! (naexps,nbexps)
  logical :: do_rotation  !!!!!!!!!!!!!!!!!!!
  integer(kind=i4_kind) :: num ! metaindex for (naexps,nbexps) > cutoff

  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: fact0_arr, &
       fact1_arr,fact2_arr ! (naexps,nbexps)
  real(kind=r8_kind),allocatable,dimension(:)  :: fact0,fact1, &
       fact2,fact3 ! (num) metaindex for (naexps,nbexps) > cutoff

  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:)  :: gamma_arg ! (num,3)
  real(kind=r8_kind),allocatable,dimension(:)    :: gamma_arg2 ! (num)
  real(kind=r8_kind),allocatable,dimension(:,:)  :: gamma_help ! (num,max_order(l)/1)
  integer(kind=i4_kind)                          :: max_order,max_gamma

  ! help arrays for solid harmincs
  real(kind=r8_kind),allocatable  :: yl_arr(:,:,:) ! (num,(lmax_abs+1)**2,n_equals)
  real(kind=r8_kind),allocatable  :: yl_arg(:,:)   ! (num,3)

  ! help variables
  real(kind=r8_kind) :: arg,charge_c
  real(kind=r8_kind),dimension(3)  :: xa,xb,xc
  integer(kind=i4_kind)  :: i_l,k,i_ind,i_cont,i_exp,i_ua,i_ea,i,f_dim,j
  integer(kind=i4_kind)  :: lmax_ch,lmax_xc,lmax_abs,alloc_stat
  integer(kind=i4_kind)  :: memstat
  integer(kind=i4_kind)  :: n_equals,n_independent_fcts, &
       n_contributing_fcts,independent_max

  ! pointing to arrays in unique_atom_symadapt_type
  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coeff(:)

  real(kind=r8_kind),allocatable  :: sym_coeff(:,:,:,:)
    ! (num,n_equals,independent_max,lmax_abs)

  ! integrals
  real(kind=r8_kind),allocatable,dimension(:) :: overlap,kinetic,nuclear ! (num)
  real(kind=r8_kind),allocatable,dimension(:,:) :: potential  !!!!!!!!!!!!!!!!!!
  real(kind=r8_kind),allocatable,dimension(:,:) :: poten_grad !!!!!!!!!!!!!!!!!!
  real(kind=r8_kind),allocatable,dimension(:,:) :: grad_poten !!!!!!!!!!!!!!!!!!
#if 0
  real(kind=r8_kind),allocatable,dimension(:,:) :: field      !!!!!!!!!!!!!!!!!!
#endif
  real(kind=r8_kind),allocatable,dimension(:) :: intermed_3c  !!!!!!!!!!!!!!!!!!

  real(kind=r8_kind),allocatable,dimension(:) ::nuclear_pc_timps
  real(kind=r8_kind),allocatable,dimension(:) :: nuc_pseudo
  type(three_center_l)                        :: xc_int
     ! xx_int%l(-1:lmax_xx)%m(num,ncexps,n_independent_fcts,n_m1,n_m2)

  real(kind=r8_kind)                              :: expmax
  logical                                         :: pseudopot_present ! same name as in UA module
  ! end of the variables for pseudopotentials
  integer(kind=i4_kind) :: i4_zero = 0_i4_kind

  logical :: split3c
  real(r8_kind),allocatable  :: this(:) ! (num)
!!$  real(r8_kind),allocatable  :: this5d(:,:,:,:,:) ! (num,1,1,nlmB,nlmA)
  real(r8_kind),allocatable  :: this6d(:,:,:,:,:,:) ! (num,1,1,nlmB,nlmA,1)
  real(r8_kind) :: zexps(1) ! for finite nuc only
  logical :: with_timps
  type(three_center_l_v2), allocatable :: coul_int(:)
  integer(i4_kind) :: i_ir, i_pa, n_pa, n_ir
  FPP_TIMER_DECL(pss)

  intrinsic max,maxval
  external error_handler

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'ss_calculate: PP=',pseudopot_present,' imode=',imode

  split3c = present(many_3c)

  naexps = ua1_basis%n_exponents
  nbexps = ua2_basis%n_exponents
  allocate(fact0_arr(nbexps,naexps),STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : allocation (1) failed")
  allocate(fact1_arr(nbexps,naexps),STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : allocation (2) failed")
  allocate(fact2_arr(nbexps,naexps),STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : allocation (3) failed")
  allocate(cutoff(nbexps,naexps),STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : allocation (4) failed")


  xa = center1
  xb = center2
  m1=n_m1
  m2=n_m2
  ASSERT(m1==1)
  ASSERT(m2==1)

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

  expmax = options_integral_expmax()
  where(fact2_arr*arg>expmax) ! cutoff: where almost no overlap
     cutoff=.false.           ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)

  if(num==0) then ! all integrals are equal zero
     if (integralpar_2cob_ol) then
        prim_int_2cob_ol = 0.0_r8_kind
     end if
     if (integralpar_2cob_potential) then             !!!!!!!!!!!!
        prim_int_2cob_poten(:,:,:,:,:) = 0.0_r8_kind  !!!!!!!!!!!!
     end if
     if (integralpar_2cob_field) then  !!!!!!!!!!
        prim_int_2cob_field(:,:,:,:,:) = 0.0_r8_kind  !!!!!!!!!
     end if
     if (integralpar_2cob_kin) then
        prim_int_2cob_kin= 0.0_r8_kind
     end if
     if (integralpar_2cob_nuc) then
        prim_int_2cob_nuc(:,:,:,:)=0.0_r8_kind
     end if
     if (integralpar_relativistic) then
        prim_int_2cob_pvsp(:,:,:,:)=0.0_r8_kind
     end if
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
          fact2_arr,cutoff,stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: deallocation fact0_arr ... failed")
    return
  end if

  allocate (fact0(num),fact1(num),fact2(num),STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE: allocation of fact(i) failed")

  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  fact0=pack(fact0_arr,cutoff)
  fact1=pack(fact1_arr,cutoff)
  fact2=pack(fact2_arr,cutoff)

  deallocate(fact0_arr,fact1_arr,fact2_arr,STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE: deallocation (1) failed")


  if ( integralpar_2cob_nuc .or. &
       integralpar_3c_xc .or. &
       integralpar_3c_co .or. &
       integralpar_2cob_potential .or. &
       integralpar_2cob_field ) then
     allocate(gamma_arg(num,3),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of gamma_arg failed")
     allocate(gamma_arg2(num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of gamma_arg2 failed")

     ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
     gamma_arg(:,1)=(pack(spread(aexps*xa(1),1,nbexps) + &
          spread(bexps*xb(1),2,naexps),cutoff))/fact0

     gamma_arg(:,2)=(pack(spread(aexps*xa(2),1,nbexps) + &
          spread(bexps*xb(2),2,naexps),cutoff))/fact0

     gamma_arg(:,3)=(pack(spread(aexps*xa(3),1,nbexps) + &
          spread(bexps*xb(3),2,naexps),cutoff))/fact0
  endif


  ! first calculating 2-center integrals----------------

  allocate(this(num),stat=memstat)
  ASSERT(memstat==0)
  if( split3c )then
     allocate(this6d(num,1,1,1,1,1),stat=memstat)
     ASSERT(memstat==0)
  endif

  ! overlap ----------------------
  allocate(overlap(num),STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE: allocation overlap failed")
  overlap = (two*sqrt(fact1)/fact0)* &
       sqrt((two*sqrt(fact1)/fact0))*exp(-fact2*arg)
  ! -----------------------------

  ! kinetic energy --------------
  if (integralpar_2cob_kin) then
     allocate(kinetic(num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE: allocation of kinetic failed")
     kinetic=(three*fact2-two*fact2**two*arg)*overlap
  endif
  ! ----------------------------

  ! nuclear attraction ---------
  if (integralpar_2cob_nuc) then
     allocate(nuclear(num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of nuclear failed")
     allocate(gamma_help(num,1),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : allocation gamma_help (1) failed")
     if (pseudopot_present)then
       allocate(nuc_pseudo(num),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of nuc_pseudo failed")
       nuc_pseudo=0.0_r8_kind
     end if
     with_timps = n_timps+pointcharge_N .gt. 0
     if(with_timps .and. pseudopot_present &
          .and.integralpar_relativistic) then
	allocate(nuclear_pc_timps(num),STAT=alloc_stat)
	if (alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE: allocation of nuclear_pc_timps failed")
	nuclear_pc_timps=0.0_r8_kind
     endif ! pointcharge_N+n_timps.ne.0

     nuclear=0.0_r8_kind
     unique_at: do i_ua=1,n_unique_atoms
        charge_c=unique_atoms(i_ua)%Z
        zc      =unique_atoms(i_ua)%zc
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'ss_calc: ua=',i_ua,', zero its charge!'
        zc = zero
        charge_c  = zero
        this = zero
        do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
           xc=unique_atoms(i_ua)%position(:,i_ea) ! coordinates
           gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                (gamma_arg(:,3)-xc(3))**2)
           gamma_help(:,1:1)=gamma(1,fact0*gamma_arg2)

           this = this + (charge_c-zc) * gamma_help(:,1)
        enddo
        this = this * two*sqrt(fact0/pi)*overlap

        if(zc.ne.zero.and.integralpar_relativistic) then 
           nuc_pseudo = nuc_pseudo + this
        else
           nuclear = nuclear + this ! this is contrib from atoms
        endif

        if( split3c )then
           call unpack_many_3c(this,i_ua,OFF_V)
        endif
     enddo unique_at
    
     
     ! pseudopotential calculation

     if(pseudopot_present.and. &
          (.not.operations_core_density)) then
        ABORT('not supported')
     end if ! endif for pseudopot_present

     ! now add contributions from point charges and timps
     do i_ua=1,pointcharge_N+n_timps
        if(i_ua<=n_timps) then
           charge_c=unique_timps(i_ua)%Z - unique_timps(i_ua)%ZC
           n_equals=unique_timps(i_ua)%n_equal_atoms
        else
           cycle ! pointcharges go to SHGI  !!!!!!!!!!!!!!AS
           charge_c=pointcharge_array(i_ua-n_timps)%Z
           n_equals=pointcharge_array(i_ua-n_timps)%n_equal_charges
        end if
        this = zero
        do i_ea=1,n_equals
           if(i_ua<=n_timps) then
              xc=unique_timps(i_ua)%position(:,i_ea)
           else
              xc=pointcharge_array(i_ua-n_timps)%position(:,i_ea) ! coordinates
           end if
           gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                (gamma_arg(:,3)-xc(3))**2)
           gamma_help(:,1:1)=gamma(1,fact0*gamma_arg2)

           this = this + charge_c * gamma_help(:,1)
        enddo
        this = this * two*sqrt(fact0/pi)*overlap
	if(integralpar_relativistic &
	  .and.pseudopot_present) then
	   nuclear_pc_timps=nuclear_pc_timps+ this
	else
           nuclear = nuclear + this ! this is PC contrib
	endif ! integralpar_relativistic
     enddo

     deallocate(gamma_help,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation gamma_help (1) failed") 

     if (pseudopot_present.and.(.not.integralpar_relativistic)) then
        nuclear = nuclear + nuc_pseudo
     end if

  endif

  !calculate integrals of electrostatic potential
  if (integralpar_2cob_potential.and.old_potential) then
     allocate(potential(N_points,num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of potential failed")
     potential = 0.0_r8_kind
     allocate(intermed_3c(num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of intermed_3c failed")
     allocate(gamma_help(num,1),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE : allocation gamma_help (1) failed")
     do i=1,N_points
        n_equals=point_in_space(i)%N_equal_points
        do i_ea=1,n_equals
           xc=point_in_space(i)%position(:,i_ea) ! coordinates
           gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                (gamma_arg(:,3)-xc(3))**2)
           gamma_help(:,1:1)=gamma(1,fact0*gamma_arg2)
           potential(i,:) = potential(i,:) + gamma_help(:,1)
        enddo
     enddo
     deallocate(gamma_help,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation gamma_help (1) failed")
     do i=1,N_points
        potential(i,:) = potential(i,:)*two*sqrt(fact0/pi)*overlap
     enddo
!!$print*,'pot',na,nb,sum(potential(:,1))
  endif

#if 0
  !calculate integrals of electrostatic field (projection)
  if (integralpar_2cob_field .and. old_elfield) then
     if(calc_normal) then
        allocate(field(N_surface_points,num),stat=alloc_stat)
     else
        allocate(field(totsym_field_length,num),stat=alloc_stat)
     end if
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of field is failed")
     field = 0.0_r8_kind
     allocate(intermed_3c(num),stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of indermed_3c for field is failed")
     lmax_ch = maxval(unique_atoms(:)%lmax_ch)
     max_order = max(lmax_ch+2,3)
     max_gamma=2
!!$     if(integralpar_relativistic) then
!!$        max_gamma=4
!!$        max_order=max(max_order,4)
!!$     end if
     allocate( grad_poten(num,3),STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation grad_poten failed")
     grad_poten = zero
     allocate (gamma_help(num,max_order),STAT=alloc_stat)
     if(alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE : allocation gamma_help(2) failed ")

     i_nsp: do i=1,N_surface_points
        n_equals=surface_points(i)%N_equal_points
        f_dim=surf_points_grad_index(i+1)-surf_points_grad_index(i)
        allocate(poten_grad(num,f_dim),stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("SS_CACLULATE : allocation failed poten_grad")
        poten_grad=0.0_r8_kind
        
        n_equals_i: do i_ea=1,n_equals
           if(f_dim == 3) then
              if(sum((surf_points_grad_info(i)%m(:,:,i_ea)-unity_matrix)**2)<&
                   1.0e-7_r8_kind) then
                 do_rotation=.false.
              else
                 do_rotation=.true.
                 rotmat=>surf_points_grad_info(i)%m(:,:,i_ea)
              endif
           else
              do_rotation=.true.
              rotmat=>surf_points_grad_info(i)%m(:,:,i_ea)
           endif

           xc=surface_points(i)%position(:,i_ea)
           gamma_arg2 = ((gamma_arg(:,1)-xc(1))**2 + &
                (gamma_arg(:,2)-xc(2))**2 + &
                (gamma_arg(:,3)-xc(3))**2)
           gamma_help(:,1:max_gamma) = gamma(max_gamma,fact0*gamma_arg2)
           do j=1,3
              grad_poten(:,j) = four*sqrt(fact0/pi)*overlap*&
                   fact0*gamma_help(:,2)* &
                   (gamma_arg(:,j)-xc(j))
           enddo

           if(do_rotation) then
              do j=1,f_dim 
                 poten_grad(:,j)=poten_grad(:,j)+&
                      rotmat(j,1)*grad_poten(:,1)+&
                      rotmat(j,2)*grad_poten(:,2)+&
                      rotmat(j,3)*grad_poten(:,3)
              enddo
           else
              do j=1,3
                 poten_grad(:,j)=poten_grad(:,j)+grad_poten(:,j)
              enddo
           end if
        enddo n_equals_i

        if(calc_normal) then
           grad_poten=0.0_r8_kind
           rotmat=>surf_points_grad_info(i)%m(:,:,1)
           do j=1,f_dim
              do k=1,3
                 grad_poten(:,k)= grad_poten(:,k)+rotmat(j,k)*poten_grad(:,j)
              enddo
           enddo

           do k=1,3
              field(i,:)=field(i,:)+grad_poten(:,k)*surface_points(i)%out_normal(k)
           enddo
        else
           k=surf_points_grad_index(i)
           do j=1,f_dim
              field(k,:)=poten_grad(:,j)
              k=k+1
           end do
        end if

        deallocate(poten_grad,stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("SS_CACLULATE : deallocation failed poten_grad")
     enddo i_nsp

     deallocate( grad_poten,STAT=alloc_stat )
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: deallocation grad_poten failed")
     deallocate (gamma_help,STAT=alloc_stat)
     if(alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE : deallocation gamma_help(2) failed ")
  endif
#endif
  ! re-map them to the int_data_2cob3c_stuff

  if (integralpar_2cob_ol) then
     prim_int_2cob_ol(:,:,m1,m2) = unpack(overlap,cutoff,zero)
  endif
  if (integralpar_2cob_kin) then
     prim_int_2cob_kin(:,:,m1,m2)= unpack(kinetic,cutoff,zero)
     deallocate(kinetic,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation kinetic failed")
  endif

  if (integralpar_2cob_nuc) then
     prim_int_2cob_nuc(:,:,m1,m2)= unpack(nuclear,cutoff,zero)
     ! 1-true nuclear, 2-pseudo,3-PC,4-EWPC
     deallocate(nuclear,STAT=alloc_stat)
     ASSERT(alloc_stat.eq.0)

     if (pseudopot_present) then
	if(integralpar_relativistic) then
           with_timps = n_timps+pointcharge_N .gt. 0
	 if( with_timps ) then
	  prim_int_2cob_nuc_pseudo(:,:,m1,m2)= & !!! prim_int_2cob_nuc_pseudo(:,:,m1,m2)+ &
	    unpack(nuc_pseudo,cutoff,zero)+ &
              unpack(nuclear_pc_timps,cutoff,zero)
      else 
      prim_int_2cob_nuc_pseudo(:,:,m1,m2)= &
					    unpack(nuc_pseudo,cutoff,zero)
	 endif ! pointcharge_N+n_timps.ne.0
	endif ! integralpar_relativistic
        deallocate(nuc_pseudo,STAT=alloc_stat)
        if (alloc_stat.ne.0) call error_handler &
            ("SS_CACLULATE : deallocation nuc_pseudo failed")
     end if
!!!     call print_nuclear()
  endif
  if (integralpar_2cob_potential.and.old_potential) then
     do i=1,N_points
        intermed_3c(:)=potential(i,:)
        prim_int_2cob_poten(:,:,i,m1,m2)= unpack(intermed_3c(:),cutoff,zero)
     enddo
!!$print*,na,nb,sum(prim_int_2cob_poten(1,1,:,1,1)),m1,m2
     deallocate(potential,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation potential failed")
     deallocate(intermed_3c,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation intermed_3c failed")
  endif
#if 0
  if (integralpar_2cob_field .and. old_elfield) then  !!!!!!!!!!!!!!!!1
     if(calc_normal) then 
        do i=1,N_surface_points
           intermed_3c(:)=field(i,:)
           prim_int_2cob_field(:,:,i,m1,m2)= unpack(intermed_3c(:),cutoff,zero)
        enddo
     else
        do i=1,totsym_field_length
           intermed_3c(:)=field(i,:)
           prim_int_2cob_field(:,:,i,m1,m2)= unpack(intermed_3c(:),cutoff,zero)
        enddo
     end if
     deallocate(field,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation field failed")
     deallocate(intermed_3c,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CACLULATE : deallocation intermed_3c for field failed")
  endif
#endif
  if (integralpar_relativistic) then
     with_timps = n_timps+pointcharge_N .gt. 0
     if( with_timps .and.pseudopot_present) then
        deallocate(nuclear_pc_timps,STAT=alloc_stat)
        if (alloc_stat.ne.0) call error_handler &
             ("SS_CACLULATE : deallocation nuclear_pc_timps failed")
     endif
  endif

  ! ----------------------------------------------------

  call integral_interrupt_2cob3c()

  ! now calculation of the 3-center integrals-----------
  three_center_integrals: if (integralpar_3c_xc .or. integralpar_3c_co) then

     allocate (fact3(num),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: allocation of fact3 failed")

     unique_atom_loop: do i_ua = 1,n_unique_atoms  ! loop over third center
        lmax_ch= unique_atoms(i_ua)%lmax_ch      ! maximum l  for chargefit
        lmax_xc= unique_atoms(i_ua)%lmax_xc      ! maximum l  for xcfit  
        ! determine the maximal angular momentum
        if (integralpar_3c_co .and. .not. integralpar_3c_xc) then
           lmax_abs=lmax_ch
        elseif (.not. integralpar_3c_co .and. integralpar_3c_xc) then
           lmax_abs=lmax_xc
        else
           lmax_abs=max(lmax_ch,lmax_xc)
        endif
        z= unique_atoms(i_ua)%z                  ! charge 
        ! --- further allocation ----------------------------------
        ! num : number of pairs(a,b) which are inside the cutoff
        ! for s-and r2-type there is only 1 indep. fct
        max_order=max((lmax_abs+1),2)

        allocate(gamma_help(num,max_order),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation gamma_help (2) ")

        if (integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
           n_ir = symmetry_data_n_irreps() 
           allocate (coul_int(n_ir),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (2) failed")

           i_ir_alloc_: DO i_ir=1,n_ir

              n_pa = symmetry_data_n_partners(i_ir)

           ncexps = unique_atoms(i_ua)%r2_ch%n_exponents
           n_independent_fcts = &
                unique_atoms(i_ua)%symadapt_partner(i_ir,0)%n_independent_fcts

              allocate(coul_int(i_ir)%l(-1:lmax_ch),STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : allocation coul_int (3) failed")
              allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,n_independent_fcts,&
                   n_m1,n_m2,n_pa),STAT=alloc_stat)              
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (3) failed")
           ncexps = unique_atoms(i_ua)%l_ch(0)%n_exponents
              allocate(coul_int(i_ir)%l(0)%m(num,ncexps,n_independent_fcts,&
                   n_m1,n_m2,n_pa),STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : allocation coul_int (4) failed")

              do i_l=1,lmax_ch
                 ncexps = unique_atoms(i_ua)%l_ch(i_l)%n_exponents
                 n_independent_fcts = &
                      unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts
                 allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                      n_m1,n_m2,n_pa),STAT=alloc_stat)
                 if(alloc_stat.ne.0) call error_handler &
                      ("SS_CALCULATE : allocation coul_int (4) failed")
              end do

              do i_l=-1,lmax_ch
                 coul_int(i_ir)%l(i_l)%m=zero
              enddo

           end do i_ir_alloc_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
        elseif(integralpar_3c_co) then
           i_ir = get_totalsymmetric_irrep()
           n_pa = 1
           allocate (coul_int(i_ir),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (2) failed")
           allocate(coul_int(i_ir)%l(-1:lmax_ch),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (3) failed")
           ncexps = unique_atoms(i_ua)%r2_ch%n_exponents
           allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,1,n_m1,n_m2,n_pa),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (3) failed")
           ncexps = unique_atoms(i_ua)%l_ch(0)%n_exponents
           allocate(coul_int(i_ir)%l(0)%m(num,ncexps,1,n_m1,n_m2,n_pa),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation coul_int (4) failed")
           do i_l=-1,0
              coul_int(i_ir)%l(i_l)%m=zero
           enddo

           do i_l=1,lmax_ch
              ncexps = unique_atoms(i_ua)%l_ch(i_l)%n_exponents
              n_independent_fcts = &
                   unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts
              allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                   n_m1,n_m2,n_pa),STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : allocation coul_int (4) failed")
              coul_int(i_ir)%l(i_l)%m = zero
           end do

        endif

        if(integralpar_3c_xc) then
           allocate (xc_int%l(-1:lmax_xc),STAT=alloc_stat) 
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation xc_int (2) failed")
           ncexps = unique_atoms(i_ua)%r2_xc%n_exponents
           allocate(xc_int%l(-1)%m(num,ncexps,1,n_m1,n_m2),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation xc_int (3) failed")
           ncexps = unique_atoms(i_ua)%l_xc(0)%n_exponents
           allocate(xc_int%l(0)%m(num,ncexps,1,n_m1,n_m2),STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation xc_int (4) failed")
           do i_l=-1,0
              xc_int%l(i_l)%m=zero
           enddo
        endif

        n_equals=unique_atoms(i_ua)%n_equal_atoms

        allocate(yl_arr(num,(lmax_abs+1)**2,n_equals),STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : allocation yl_arr failed")

        if( split3c )then
           zexps(1) = (3.0_r8_kind/2.0_r8_kind) / unique_atoms(i_ua)%nuclear_radius**2
           this6d = zero
        endif

        equal_atoms: do i_ea=1,n_equals

           xc=unique_atoms(i_ua)%position(:,i_ea) ! coordinates
           ! gamma_arg2 = (gamma_arg**2 - vec_c**2)*(a+b)
           gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                (gamma_arg(:,3)-xc(3))**2)



           ! s-type and r2-type exchange fit integrals
           if(integralpar_3c_xc) call s_r2_xc()

           ! do a pre-calculation of the solid harmonics --------
           ! calculate them for ALL equal_atoms and for l_max_abs
           ! it will cost more storage but will be faster
           if (lmax_abs.gt.0) then
              allocate(yl_arg(num,3),STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : allocation yl_arg failed")
              xc=unique_atoms(i_ua)%position(:,i_ea)
              yl_arg(:,1) = gamma_arg(:,1) - xc(1)
              yl_arg(:,2) = gamma_arg(:,2) - xc(2)
              yl_arg(:,3) = gamma_arg(:,3) - xc(3)

              yl_arr(:,:,i_ea) = solid_harmonics_calc(lmax_abs,yl_arg)
              deallocate(yl_arg,STAT=alloc_stat) ! this is not needed anymore
              if (alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation yl_arg failed")
           endif

        enddo equal_atoms

        independent_max = maxval(unique_atoms(i_ua)%symadapt_partner(1,:)%n_independent_fcts)

        lmax_gt_zero: if (lmax_abs.gt.0) then 
        
           allocate(sym_coeff(num,n_equals,independent_max,lmax_abs),STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : allocation sym_coeff failed")

           ! precalculate symmetry adaption coefficient sym_coeff
           call l_fit_symmetry_adapt()

           ! l-type exchange fit integrals
           if(integralpar_3c_xc) call l_xc()

           deallocate(sym_coeff,STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE: deallocation sym_coeff failed ")

        endif lmax_gt_zero

        !
        ! - CO calculation
        !    

        if (integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE

           i_ir_: DO i_ir=1,symmetry_data_n_irreps()
              i_pa_: DO i_pa=1,symmetry_data_n_partners(i_ir)

                 ea_: do i_ea=1,n_equals

                    xc=unique_atoms(i_ua)%position(:,i_ea) ! coordinates
                    ! gamma_arg2 = (gamma_arg**2 - vec_c**2)*(a+b)
                    gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                         (gamma_arg(:,3)-xc(3))**2)

                    i_l=0

                    n_independent_fcts = &
                         unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts

                    if (n_independent_fcts .ne. 0) then

                       ! s-type and r2-type coulomb fit integrals
                       call s_coulomb( &
                            unique_atoms(i_ua)%l_ch(0)%exponents(:), &
                            coul_int(i_ir)%l(0)%m)
                       call r2_coulomb()

                       if( split3c )then
                          ! adds this EA only:
                          call s_coulomb( zexps, this6d )
                       endif

                       allocate(yl_arg(num,3),STAT=alloc_stat)
                       if(alloc_stat.ne.0) call error_handler &
                            ("SS_CALCULATE : allocation yl_arg failed")
                       xc=unique_atoms(i_ua)%position(:,i_ea)
                       yl_arg(:,1) = gamma_arg(:,1) - xc(1)
                       yl_arg(:,2) = gamma_arg(:,2) - xc(2)
                       yl_arg(:,3) = gamma_arg(:,3) - xc(3)
                       yl_arr(:,:,i_ea) = solid_harmonics_calc(lmax_abs,yl_arg)
                       deallocate(yl_arg,STAT=alloc_stat) ! this is not needed anymore
                       if (alloc_stat.ne.0) call error_handler &
                            ("SS_CALCULATE : deallocation yl_arg failed")

                    end if

                 end do ea_

                 do i_l=1,lmax_ch
                    n_independent_fcts = &
                         unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts
                    if (n_independent_fcts .ne. 0) then

                       allocate(sym_coeff(num,n_equals,n_independent_fcts,i_l),&
                            STAT=alloc_stat)
                       if (alloc_stat.ne.0) call error_handler &
                            ("SS_CALCULATE : allocation sym_coeff failed")

                       ! precalculate symmetry adaption coefficient sym_coeff
                       call l_fit_symmetry_adapt_v2(i_l,i_ir)

                       ! l-type coulomb fit integrals
                       call l_coulomb()

                       deallocate(sym_coeff,STAT=alloc_stat)
                       if (alloc_stat.ne.0) call error_handler &
                            ("SS_CALCULATE: deallocation sym_coeff failed ")

                    end if

                 end do

              end do i_pa_
           end do i_ir_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
        elseif (integralpar_3c_co) then
           i_ir = get_totalsymmetric_irrep()
           i_pa = 1
           ea2_: do i_ea=1,n_equals

              xc=unique_atoms(i_ua)%position(:,i_ea) ! coordinates
              ! gamma_arg2 = (gamma_arg**2 - vec_c**2)*(a+b)
              gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                   (gamma_arg(:,3)-xc(3))**2)

              i_l=0

              n_independent_fcts = &
                   unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts

              if (n_independent_fcts .ne. 0) then

                 ! s-type and r2-type coulomb fit integrals
                 call s_coulomb( &
                      unique_atoms(i_ua)%l_ch(0)%exponents(:), &
                      coul_int(i_ir)%l(0)%m )
                 call r2_coulomb()

                 if( split3c )then
                    ! adds this EA only:
                    call s_coulomb( zexps, this6d )
                 endif

                 allocate(yl_arg(num,3),STAT=alloc_stat)
                 if(alloc_stat.ne.0) call error_handler &
                      ("SS_CALCULATE : allocation yl_arg failed")
                 xc=unique_atoms(i_ua)%position(:,i_ea)
                 yl_arg(:,1) = gamma_arg(:,1) - xc(1)
                 yl_arg(:,2) = gamma_arg(:,2) - xc(2)
                 yl_arg(:,3) = gamma_arg(:,3) - xc(3)
                 yl_arr(:,:,i_ea) = solid_harmonics_calc(lmax_abs,yl_arg)
                 deallocate(yl_arg,STAT=alloc_stat) ! this is not needed anymore
                 if (alloc_stat.ne.0) call error_handler &
                      ("SS_CALCULATE : deallocation yl_arg failed")

              end if
           end do ea2_

           do i_l=1,lmax_ch
              n_independent_fcts = &
                   unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts
              if (n_independent_fcts .ne. 0) then

                 allocate(sym_coeff(num,n_equals,n_independent_fcts,i_l),&
                      STAT=alloc_stat)
                 if (alloc_stat.ne.0) call error_handler &
                      ("SS_CALCULATE : allocation sym_coeff failed")

                 ! precalculate symmetry adaption coefficient sym_coeff
                 call l_fit_symmetry_adapt_v2(i_l,i_ir)

                 ! l-type coulomb fit integrals
                 call l_coulomb()

                 deallocate(sym_coeff,STAT=alloc_stat)
                 if (alloc_stat.ne.0) call error_handler &
                      ("SS_CALCULATE: deallocation sym_coeff failed ")
              end if

           end do

        end if

        !
        ! - end of CO calculation
        !


        deallocate(yl_arr,STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : deallocation yl_arr failed")
        deallocate(gamma_help,STAT=alloc_stat)
        if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : deallocation gamma_help (2) failed")

        if(integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
           call fitcontract_v2(num,i_ua,cutoff,coul_int)
           do i_ir=1,symmetry_data_n_irreps()
           do i_l = -1, lmax_ch
                 deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation coul_int%l%m failed")
           enddo
              deallocate (coul_int(i_ir)%l,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation coul_int%l failed")
           end do
           deallocate (coul_int, STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : deallocation coul_int%l failed")
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
        elseif(integralpar_3c_co) then
           !!SB: shouldbe only for resp, but...
           call fitcontract_v2(num,i_ua,cutoff,coul_int)
           i_ir = get_totalsymmetric_irrep()
           do i_l = -1, lmax_ch
              deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation coul_int%l%m failed")
           end do
           deallocate (coul_int(i_ir)%l,STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
             ("SS_CALCULATE : deallocation coul_int%l failed")
           deallocate (coul_int, STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : deallocation coul_int%l failed")

        endif

        if(integralpar_3c_xc) then
           call fitcontract('xc',num,i_ua,cutoff,xc_int)
           do i_l = -1, lmax_ch
              deallocate(xc_int%l(i_l)%m,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("SS_CALCULATE : deallocation xc_int%l%m failed")
           enddo
           deallocate (xc_int%l,STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("SS_CALCULATE : deallocation xc_int%l failed")
        endif

     enddo unique_atom_loop
  
     deallocate (fact3,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: deallocation of fact3 failed")

  endif three_center_integrals


  deallocate (fact0,fact1,fact2,STAT=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE: deallocation of fact(i) failed")

  deallocate(overlap,STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : deallocation of overlap failed")

  deallocate(cutoff,STAT=alloc_stat)
  if( alloc_stat.ne.0) call error_handler &
       ("SS_CALCULATE : deallocation of cutoff failed")

  if ( integralpar_2cob_nuc .or. &
       integralpar_3c_xc .or. &
       integralpar_3c_co .or. &
       integralpar_2cob_potential .or. &
       integralpar_2cob_field ) then
     deallocate(gamma_arg,gamma_arg2,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("SS_CALCULATE: deallocation of gamma_arg failed")
  endif

  deallocate(this,stat=memstat)
  ASSERT(memstat==0)
  if( split3c )then
     deallocate(this6d,stat=memstat)
     ASSERT(memstat==0)
  endif

contains

  subroutine unpack_many_3c(this,i,offset)
    implicit none
    real(r8_kind), intent(in)    :: this(:)
    integer(i4_kind), intent(in) :: i,offset
    ! *** end of iterface ***

    integer(i4_kind) :: pos

    pos = OFF_STRIDE*(i-1)+offset
    ASSERT(pos<=size(many_3c,3))
    ! remap ints to the appropriate data structure
    many_3c(:,:,pos,1,1)=unpack(this(:),cutoff,zero)
  end subroutine unpack_many_3c

  !**************************************************************

  subroutine s_coulomb(cexps,coul)
    ! s-type coulomb fit integrals
    ! differs from LL and LS counterparts by not
    ! summing over EAs
    real(r8_kind), intent(in)    :: cexps(:)
    real(r8_kind), intent(inout) :: coul(:,:,:,:,:,:)
    ! *** end of iterface ***

    integer(i4_kind) :: ncexps

    integer(i4_kind)           :: n_cf, n_if, i_ind, ic
    real(r8_kind)              :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if  =  unique_atoms(i_ua)%symadapt_partner(i_ir,0)%n_independent_fcts

    do i_ind = 1, n_if

       n_cf    =  unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%n_fcts
       coeff   => unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%c
       eq_atom => unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%I_equal_atom

    sc = 0.0_r8_kind
    do ic = 1, n_cf
       sc(eq_atom(ic)) = sc(eq_atom(ic)) + coeff(ic)
    end do

       ncexps = size(cexps)

!!$    ncexps = unique_atoms(i_ua)%l_ch(0)%n_exponents
!!$    cexps => unique_atoms(i_ua)%l_ch(0)%exponents(:)

    do k=1,ncexps ! loop over exponents of third center, s-fit-fct.
       fact1 = fact0    / (fact0+cexps(k))      ! = (a+b)   / (a+b+c)
       fact3 = cexps(k) * fact1              ! = (a+b)*c / (a+b+c)
       gamma_help(:,1:1) = gamma(1,gamma_arg2*fact3)

          coul(:,k,i_ind,m1,m2,i_pa) = coul(:,k,i_ind,m1,m2,i_pa) + &
               ( two * pi/cexps(k)     * &
               sqrt(fact1) * overlap * &
               gamma_help(:,1) &
               ) * NORM * sc(i_ea) 
          
       enddo
    end do
    call integral_interrupt_2cob3c()
  end subroutine s_coulomb

  subroutine r2_coulomb()
    ! r2-type coulomb fit integrals
    ! differs from LL and LS counterparts by not
    ! summing over EAs

    integer(i4_kind)           :: n_cf, i_ind, n_if, ic
    real(r8_kind)              :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if  =  unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts

    do i_ind = 1, n_if
       n_cf    =  unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%n_fcts
       coeff   => unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%c
       eq_atom => unique_atoms(i_ua)%symadapt_partner(i_ir,0)%symadapt(i_ind,i_pa)%I_equal_atom

    sc = 0.0_r8_kind
    do ic = 1, n_cf
       sc(eq_atom(ic)) = sc(eq_atom(ic)) + coeff(ic)
    end do

       ncexps = unique_atoms(i_ua)%r2_ch%n_exponents
       cexps => unique_atoms(i_ua)%r2_ch%exponents(:)
       do k=1,ncexps
          fact1=fact0/(fact0+cexps(k))              ! = (a+b)/(a+b+c)
          fact2=(fact0+1.5_r8_kind*cexps(k))/fact0  ! = (a+b+3/2c)/(a+b)
          fact3=cexps(k)*fact1                      ! = (a+b)*c/(a+b+c)
          gamma_help(:,1:2)=gamma(2,gamma_arg2*fact3)
          
          coul_int(i_ir)%l(-1)%m(:,k,i_ind,m1,m2,i_pa) = & 
               coul_int(i_ir)%l(-1)%m(:,k,i_ind,m1,m2,i_pa) + &
               ( two*pi/(cexps(k)**2)*overlap*fact1*sqrt(fact1)* &
               ( fact2*gamma_help(:,1) + &
               fact3*gamma_arg2* &
               gamma_help(:,2) ) &  
               ) * NORM * sc(i_ea)
          
       enddo
    end do
    call integral_interrupt_2cob3c()
  end subroutine r2_coulomb

  !**************************************************************

  subroutine s_r2_xc()
    ! s-type and r2-type xc fit integral

    ncexps = unique_atoms(i_ua)%l_xc(0)%n_exponents
    cexps => unique_atoms(i_ua)%l_xc(0)%exponents(:)       
    do k=1,ncexps ! loop over s-type fitfcts.
       fact1=fact0/(fact0+cexps(k))      ! = (a+b)/(a+b+c)
       fact3=cexps(k)*fact1              ! = (a+b)*c/(a+b+c)
       xc_int%l(0)%m(:,k,1,m1,m2) = xc_int%l(0)%m(:,k,1,m1,m2) + &
            fact1*sqrt(fact1) * exp( -(fact3*gamma_arg2) ) &
            * overlap
    enddo
    call integral_interrupt_2cob3c()

    ncexps = unique_atoms(i_ua)%r2_xc%n_exponents
    cexps => unique_atoms(i_ua)%r2_xc%exponents(:)  
    do k=1,ncexps ! loop over r2-type fitfcts.
       fact1=fact0/(fact0+cexps(k))      ! = (a+b)/(a+b+c)
       fact2=(fact0+three/two*cexps(k))/fact0  ! = (a+b+3/2c)/(a+b)
       fact3=cexps(k)*fact1              ! = (a+b)*c/(a+b+c)

       xc_int%l(-1)%m(:,k,1,m1,m2) = xc_int%l(-1)%m(:,k,1,m1,m2) + &
            fact0/(fact0+cexps(k))*sqrt(fact0/(fact0+cexps(k))) &
            *overlap/(fact0+cexps(k))* &
            (three/two &
            + fact0**2/(fact0+cexps(k))*gamma_arg2 )* &
            exp(-(fact3*gamma_arg2))
    enddo
    call integral_interrupt_2cob3c()
  end subroutine s_r2_xc

  !**************************************************************

     subroutine l_fit_symmetry_adapt_v2(i_l,i_ir)
       implicit none
       integer(i4_kind), intent(in) :: i_l,i_ir
       ! *** end of interface ***

       ! eclipse global vars:
       integer(i4_kind) :: n_independent_fcts, i_ind
       integer(i4_kind) :: n_contributing_fcts, i_cont
       integer(i4_kind), pointer :: magn(:), eq_atom(:)
       real(r8_kind), pointer    :: coeff(:)

       ! precalculate symmetry adaption coefficient sym_coeff
       ! for l-type xc and coulomb fit integral

       sym_coeff = 0.0_r8_kind

       n_independent_fcts = &
            unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%n_independent_fcts

       independents :do i_ind = 1,n_independent_fcts
          n_contributing_fcts = &
               unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%n_fcts
          coeff => unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%c
          magn  => unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%m
          eq_atom => unique_atoms(i_ua)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%I_equal_atom

          ! construct a coefficient sym_coeff(:,i_equal,i_ind,i_l)
          ! sum over contributing functions
          do i_cont=1,n_contributing_fcts
             sym_coeff(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coeff(:,eq_atom(i_cont),i_ind,i_l) + &
                  coeff(i_cont)*yl_arr(:,i_l**2+magn(i_cont),eq_atom(i_cont))
          enddo
       enddo independents

     end subroutine l_fit_symmetry_adapt_v2

  subroutine l_fit_symmetry_adapt
    implicit none
    ! *** end of interface ***

    ! eclipse global vars:
    integer(i4_kind) :: i_l
    integer(i4_kind) :: n_independent_fcts, i_ind
    integer(i4_kind) :: n_contributing_fcts, i_cont
    integer(i4_kind), pointer :: magn(:), eq_atom(:)
    real(r8_kind), pointer    :: coeff(:)

    ! precalculate symmetry adaption coefficient sym_coeff
    ! for l-type xc and coulomb fit integral

    sym_coeff = 0.0_r8_kind

    ang_momentum : do i_l=1,lmax_abs ! now loop over higher ang. momenta
       ! first calculate the part which is independent of 'ch' or 'xc'
       ! i.e. the symm_coeffs which do not contain
       ! cexps (i.e. gamma)
       n_independent_fcts = &
            unique_atoms(i_ua)%symadapt_partner(1,i_l)%n_independent_fcts

       independents :do i_ind = 1,n_independent_fcts
          n_contributing_fcts = &
               unique_atoms(i_ua)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%n_fcts
          coeff => unique_atoms(i_ua)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c
          magn => unique_atoms(i_ua)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m
          eq_atom => unique_atoms(i_ua)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%I_equal_atom

          ! construct a coefficient sym_coeff(:,i_equal,i_ind,i_l)
          ! sum over contributing functions
          do i_cont=1,n_contributing_fcts
             sym_coeff(:,eq_atom(i_cont),i_ind,i_l) = &
                  sym_coeff(:,eq_atom(i_cont),i_ind,i_l) + &
                  coeff(i_cont)*yl_arr(:,i_l**2+magn(i_cont),eq_atom(i_cont))
          enddo
       enddo independents

    enddo ang_momentum
  end subroutine l_fit_symmetry_adapt

  !**************************************************************

  subroutine l_coulomb()
    ! l-type coulomb fit integrals

       ncexps = unique_atoms(i_ua)%l_ch(i_l)%n_exponents
       cexps => unique_atoms(i_ua)%l_ch(i_l)%exponents(:)

       ! now loop over exponents gamma
       c_exponents_ch: do i_exp=1,ncexps
          fact1=fact0/(fact0+cexps(i_exp))  ! = (a+b)/(a+b+c)
          fact3=cexps(i_exp)*fact1          ! = (a+b)*c/(a+b+c)

          ! loop now over ALL centers: it is then possible to
          ! use the loop over independents as the innermost loop
          ! and thus we do not have to calculate the gamma-function
          ! for each indep. fct.

          ! now sum over equals
          equals_ch :do i_ea=1,n_equals
             xc=unique_atoms(i_ua)%position(:,i_ea)
             gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                  (gamma_arg(:,3)-xc(3))**2)*fact3
             gamma_help(:,1:i_l+1)=gamma(i_l+1,gamma_arg2)

             do i_ind=1,n_independent_fcts
                coul_int(i_ir)%l(i_l)%m(:,i_exp,i_ind,m1,m2,i_pa) = &
                     coul_int(i_ir)%l(i_l)%m(:,i_exp,i_ind,m1,m2,i_pa)+ &
                     two*pi*overlap/cexps(i_exp)*sqrt(fact1) * &
                     (fact3+fact3)**i_l * &
                     sym_coeff(:,i_ea,i_ind,i_l)*gamma_help(:,i_l+1)
             enddo
          enddo equals_ch
       enddo c_exponents_ch
       call integral_interrupt_2cob3c()
  end subroutine l_coulomb

  !**************************************************************

  subroutine l_xc()
    ! l-type exchange fit integrals
    ang_momentum_xc: do i_l = 1,lmax_xc
       n_independent_fcts = &
            unique_atoms(i_ua)%symadapt_partner(1,i_l)%n_independent_fcts
       ncexps = unique_atoms(i_ua)%l_xc(i_l)%n_exponents
       cexps => unique_atoms(i_ua)%l_xc(i_l)%exponents(:)
       allocate(xc_int%l(i_l)%m(num,ncexps,n_independent_fcts,n_m1,n_m2),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("SS_CALCULATE: allocation xc_int (6) failed")
       xc_int%l(i_l)%m = 0.0_r8_kind

       ! now loop over exponents gamma
       c_exponents_xc: do i_exp=1,ncexps
          fact1=fact0/(fact0+cexps(i_exp))  ! = (a+b)/(a+b+c)
          fact3=cexps(i_exp)*fact1          ! = (a+b)*c/(a+b+c)

          equals_xc :do i_ea=1,n_equals
             xc=unique_atoms(i_ua)%position(:,i_ea)
             gamma_arg2=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
                  (gamma_arg(:,3)-xc(3))**2)

             do i_ind=1,n_independent_fcts
                xc_int%l(i_l)%m(:,i_exp,i_ind,m1,m2) = &
                     xc_int%l(i_l)%m(:,i_exp,i_ind,m1,m2)+ &
                     fact1*sqrt(fact1) * (two*fact3)**i_l * &
                     sym_coeff(:,i_ea,i_ind,i_l) * &
                     exp(-fact3*gamma_arg2) * overlap
             enddo
          enddo equals_xc
       enddo c_exponents_xc
       call integral_interrupt_2cob3c()
    enddo ang_momentum_xc
  end subroutine l_xc

  !**************************************************************


  !**************************************************************!
  ! Just to print nuclear integrals !
  subroutine print_nuclear()
    use iounitadmin_module
    integer(kind=i4_kind) :: i_exp1, i_exp2, nuc_int_unit, ma,mb
    logical               :: file_exists
    
    nuc_int_unit = get_iounit()
    inquire(file="NUC_INT", exist=file_exists)
    if(file_exists) then
       open(nuc_int_unit, file="NUC_INT",position="append")
    else
       open(nuc_int_unit, file="NUC_INT")
    end if    
    write(nuc_int_unit,"('********* SS ********** [',2i2,' |',2i2,' ] ********* SS ***********')")  na, 0, nb, 0
    do i_exp1 = 1,naexps
       do i_exp2 = 1,nbexps
          do mb=1,n_m1    !
             do ma=1,n_m2 ! 
                write(nuc_int_unit,"(5x,'exp2 = ',i3,3x,'exp1 = ',i3,3x, 'm',2i3,4x,f26.20)") &
                     i_exp2, i_exp1, mb, ma, prim_int_2cob_nuc(i_exp2,i_exp1,mb,ma)
             end do
          end do
       end do
    end do
    close(nuc_int_unit, status="keep")    
    call return_iounit(nuc_int_unit)

  end subroutine print_nuclear
  !**************************************************************!

end subroutine ss_calculate

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
subroutine ls_calculate(na,nb,la_in,lb,imode,many_3c)
  !
  !  Purpose: calculation of all primitive 2 center orbital
  !           and 3 center integrals for a given set of indizes
  !       (unique_atom1,unique_atom2,l1,equal_atom2).
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
  use symmetry_data_module,only: symmetry_data_n_irreps,  &
       symmetry_data_n_partners,& 
       get_totalsymmetric_irrep 
  use calc3c_switches,only: old_potential,old_3c_co,old_elfield
  use shgi_cntrl, only: IPSEU
  implicit none


  !== Interrupt end of public interface of module ====================

  integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
  integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
  integer(kind=i4_kind),intent(in) :: la_in ! angular momentum of unique atom a
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
       three_half=1.5_r8_kind

! real(kind=r8_kind),parameter,dimension(0:8) :: dfac= &
  real(kind=r8_kind),dimension(0:8) :: dfac= &
       (/ 1.0_r8_kind, 1.0_r8_kind, 3.0_r8_kind, 15.0_r8_kind, 105.0_r8_kind,&
       945.0_r8_kind, 10395.0_r8_kind, 135135.0_r8_kind, 2027025.0_r8_kind /)

  logical,allocatable   :: cutoff(:,:)

  integer(kind=i4_kind) :: num,counter,m,ma,alloc_stat
  integer(kind=i4_kind) :: memstat
  ! help factors
  real(kind=r8_kind),allocatable,dimension(:,:):: fact0_arr, &
       fact1_arr,fact2_arr,fact10
  real(kind=r8_kind),allocatable,dimension(:)  :: fact0,fact1, &
       fact2,fact4,fact5,fact6,fact7,fact8,rcsabc,tau

  ! help arrays for gamma-function
  real(kind=r8_kind),allocatable,dimension(:,:)     :: gamma_arg,gamma_arg2
  real(kind=r8_kind),allocatable,dimension(:,:,:)   :: gamma_help

  ! help arrays for solid harmincs
  real(kind=r8_kind),allocatable                 :: yl_arr(:,:,:)
  real(kind=r8_kind),allocatable                 :: clmamb(:,:),clmamb_scalar(:)

  ! help arrays for product_rule and diff_rule
  real(kind=r8_kind),allocatable                 :: prod_arr(:,:,:,:),&
       diff_arr(:,:),intermediate(:,:,:,:,:)

  real(kind=r8_kind) :: arg

  ! cartesian coordinates
  real(kind=r8_kind),dimension(3)  :: xa,xb,xc,xd

  logical :: laltlb  ! flag to decide if la is lower then lb or not
  integer(kind=i4_kind)  :: i,j,i_la,i_lb,i_l,k,i_ind,i_cnt,nm_la,nm_lb,&
       la_org,lb_org,la,l !,l_cg
  integer(kind=i4_kind)  :: lmax_ch,lmax_xc,lmax_abs,ly_max
  integer(kind=i4_kind)  :: n_equals,n_independent_fcts, &
       n_contributing_fcts

  integer(kind=i4_kind),pointer   :: eq_atom(:),magn(:)
  real(kind=r8_kind),pointer      :: coeff(:)

  real(kind=r8_kind),allocatable :: overlap(:,:), &
       kinetic(:),aexp_arr(:)
  real(kind=r8_kind),allocatable,dimension(:,:)    :: nuc 
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: potential 
#if 0
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: field  
#endif
  real(kind=r8_kind),allocatable,dimension(:) :: intermed_3c
  real(kind=r8_kind),allocatable  :: prod_arr_gr(:,:,:,:),help_arr_gr(:,:,:,:)
  real(kind=r8_kind),allocatable  :: diff_arr_grad(:,:,:),rsaba(:),overlap1(:,:) 
  real(kind=r8_kind),allocatable,dimension(:,:)    :: nuc_pc_timps

  real(kind=r8_kind),allocatable,dimension(:,:)    :: nuc_pseudo
  real(kind=r8_kind)                              :: expmax
  logical                                         :: pseudopot_present ! same name is in UA module
  type(three_center_l)     :: xc_int
  type(three_center_l_v2),allocatable :: coul_int(:)
  type(unique_atom_type), pointer  :: ua_pointer

  logical :: split3c
  real(r8_kind), allocatable :: this(:,:) ! (num,nlma)
!!$  real(r8_kind), allocatable :: this5d(:,:,:,:,:) ! (num,1,1,nlmB,nlmA)
  real(r8_kind), allocatable :: this6d(:,:,:,:,:,:) ! (num,1,1,nlmB,nlmA)
  real(r8_kind) :: zexps(1) ! for finite nuc only
  logical :: with_timps
  integer(i4_kind) :: N_length !!!!!!!!!!!!
  FPP_TIMER_DECL(pls)

  integer(i4_kind) :: i_ir, i_pa, n_pa

  intrinsic max

  pseudopot_present = IAND(imode,IPSEU) .ne. 0
  DPRINT 'ls_calculate: PP=',pseudopot_present,' imode=',imode

  split3c = present(many_3c)

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

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("LS_CALCULATE: allocation fact0_arr ... failed")   
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

  expmax = options_integral_expmax()
  where(fact2_arr*arg>= expmax ) ! cutoff: where almost no overlap
     cutoff=.false.              ! is present calculation is not necessary
  elsewhere
     cutoff=.true.
  end where

  num=count(cutoff)

  if(num==0) then    ! all integrals are equal zero
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
        prim_int_2cob_poten(:,:,:,:,:)=0.0_r8_kind   !!!!!!!!!!!!!!
     end if
     if (integralpar_2cob_field ) then
        prim_int_2cob_field(:,:,:,:,:)=0.0_r8_kind   !!!!!!!!!!!!!!
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
          fact2_arr,cutoff,stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE: deallocation fact0_arr ... failed")
     return
  end if

  allocate ( &
       fact0(num), fact1(num), fact2(num), fact4(num), &
       fact5(num), fact6(num), fact7(num), fact8(num), &
       aexp_arr(num), rcsabc(num), tau(num),&
       overlap(num,2*la+1), &
       gamma_arg(num,3),&
       clmamb_scalar((la+1)**2), &
       clmamb(num,(la+1)**2), &
       diff_arr(num,(la+1)**2), &
       stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("LS_CALCULATE: allocation fact0 ... failed")

  allocate(this(num,2*la+1),stat=memstat)
  ASSERT(memstat==0)
  if( split3c )then
     ! original (true) dimensions nm_la, nm_lb:
     allocate(this6d(num,1,1,nm_lb,nm_la,1),stat=memstat)
     ASSERT(memstat==0)
  endif

  if (integralpar_2cob_kin) then
     allocate(kinetic(num),stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE: allocation kinetic failed")
  end if

  if (integralpar_2cob_nuc) then
     allocate(nuc(num,2*la+1),stat=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LS_CALCULATE: allocation  failed")
     nuc=0.0_r8_kind
     if (pseudopot_present) then
         allocate(nuc_pseudo(num,2*la+1),&
                  stat=alloc_stat)
         if (alloc_stat.ne.0) call error_handler &
            ("LS_CALCULATE: allocation(nuc_pseudo) failed")
         nuc_pseudo=0.0_r8_kind
     end if
     with_timps = pointcharge_N +n_timps .gt. 0
     if ( integralpar_relativistic .and. pseudopot_present .and. &
	   with_timps ) then
	allocate(nuc_pc_timps(num,2*la+1),&
                  stat=alloc_stat)
	if (alloc_stat.ne.0) call error_handler &
            ("LS_CALCULATE: allocation(nuc_pc_timps) failed")
         nuc_pc_timps=0.0_r8_kind
     endif ! 
  end if ! integralpar_2cob_nuc

  if (integralpar_2cob_potential.and.old_potential) then
     allocate(potential(N_points,num,2*la+1),stat=alloc_stat)  
     if (alloc_stat.ne.0) call error_handler &                
          ("LS_CALCULATE: allocation of potential failed")   
     potential=0.0_r8_kind                                  
     allocate(intermed_3c(num),stat=alloc_stat)  
     if (alloc_stat.ne.0) call error_handler &                
          ("LS_CALCULATE: allocation of intermed_3c failed")   
  end if

#if 0
  if (integralpar_2cob_field.and.old_elfield) then
     if( calc_normal) then
        N_length=N_surface_points
     else
        N_length=totsym_field_length
     end if
     allocate(field(N_length,num,2*la+1),stat=alloc_stat)   
     if (alloc_stat.ne.0) call error_handler &                     
          ("LS_CALCULATE: allocation of field failed")
     field=0.0_r8_kind                                           
     allocate(intermed_3c(num),stat=alloc_stat)   
     if (alloc_stat.ne.0) call error_handler &                     
          ("LS_CALCULATE: allocation of intermed_3c for field failed")            
  end if
#endif

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

  deallocate(fact0_arr,fact1_arr,fact2_arr,stat=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("LS_CALCULATE: deallocation fact0_arr ... failed")

  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
  gamma_arg(:,1)=(pack(spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps),cutoff))/fact0

  gamma_arg(:,2)=(pack(spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps),cutoff))/fact0

  gamma_arg(:,3)=(pack(spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps),cutoff))/fact0

  if(laltlb) then
     xd=-xd
  end if


  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(la,xd)

  fact4=1.0_r8_kind
  counter=1
  do i=0,la
     do m=1,2*i+1
        clmamb(:,counter)=clmamb_scalar(counter)*fact4
        counter=counter+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo


  ! first calculating 2-center integrals----------------
  tau=fact2*arg
  fact5=fact2*(3.0_r8_kind-2.0_r8_kind*tau+2*la)  ! a*b/(a+b)(3-2*tau+2*l)
  fact6=1.0_r8_kind/sqrt(aexp_arr**la*dfac(la))*exp(-tau)*&
       (4.0_r8_kind*fact2/fact0)**0.75_r8_kind ! 
  fact2=2.0_r8_kind*fact6*sqrt(fact0/pi)

  if(.not.laltlb) then
     if (integralpar_2cob_ol) then
        magnetic_number: do m=1,2*la+1
           ! overlap
           overlap(:,m)=fact6*clmamb(:,(la)**2+m)
           ! re-map them to the int_data_2cob3c_stuff
           prim_int_2cob_ol(:,:,1,m) = unpack(overlap(:,m),cutoff,zero)    
        end do magnetic_number
     end if
     if (integralpar_2cob_kin) then
        magnetic_number_kin: do m=1,2*la+1
           ! kinetic energy
           kinetic=fact5*overlap(:,m)
           ! re-map them to the int_data_2cob3c_stuff
           prim_int_2cob_kin(:,:,1,m)=unpack(kinetic,cutoff,zero)
        end do magnetic_number_kin
     endif
  else
     if (integralpar_2cob_ol) then
        magnetic_number_2: do m=1,2*la+1
           ! overlap
           overlap(:,m)=fact6*clmamb(:,(la)**2+m)
           ! re-map them to the int_data_2cob3c_stuff
           prim_int_2cob_ol(:,:,m,1) = unpack(overlap(:,m),cutoff,zero)    
        end do magnetic_number_2
     end if

     if (integralpar_2cob_kin) then
        magnetic_number_kin_2: do m=1,2*la+1
           ! kinetic energy
           kinetic=fact5*overlap(:,m)
           ! re-map them to the int_data_2cob3c_stuff
           prim_int_2cob_kin(:,:,m,1)=unpack(kinetic,cutoff,zero)
        end do magnetic_number_kin_2
     endif
  end if

  if (integralpar_2cob_kin) then
     deallocate(kinetic,STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("LS_CACLULATE : deallocation kinetic failed")
  endif

!!!  deallocate(tau,fact5,overlap,stat=alloc_stat)
  deallocate(tau,fact5,stat=alloc_stat)
  if (alloc_stat.ne.0) call error_handler &
       ("LS_CACLULATE : deallocation tau ... failed")

  call integral_interrupt_2cob3c()


  ! calculate integrals that involve third center,
  ! i.e. fit integrals, nuclear attraction and relativistic pv scalar p
  third_center_required: if ( integralpar_2cob_nuc .or. integralpar_3c_xc &
       .or. integralpar_3c_co .or. integralpar_relativistic &
       ) then

     unique_atom_loop: do i=1,n_unique_atoms + n_timps  ! loop over third center
        if(i<=n_unique_atoms) then
           ua_pointer=>unique_atoms(i)
           lmax_ch= ua_pointer%lmax_ch      ! maximum l  for chargefit
           lmax_xc= ua_pointer%lmax_xc      ! maximum l  for xcfit  
           ! determine the maximal angular momentum
           if (integralpar_3c_xc) then
              lmax_abs=lmax_ch
           else
              lmax_abs=max(lmax_ch,lmax_xc)
           endif
           ly_max=max(la,lmax_ch,lmax_xc)
           max_order=max(1+la+lmax_abs,3+la)
           z= ua_pointer%z                  ! charge 
           zc = ua_pointer%zc               ! charge
           n_equals=unique_atoms(i)%n_equal_atoms
        ! NUC and PP is handled by SHGI, skip the NUC:
        DPRINT   'ls_calc: ua=',i,', zero its charge!'
        zc = zero
        z  = zero
        allocate ( &
             gamma_help(num,max_order,n_equals), &
             gamma_arg2(num,n_equals), &
             fact10(num,(ly_max+1)**2), &
             stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("LL_CACLULATE : allocation gamma_help ... failed")
        counter=1
        fact4=1.0_r8_kind
        do i_l=0,ly_max
           do ma=1,2*i_l+1
              fact10(:,counter)=fact4
              counter=counter+1
           enddo
           fact4=fact4*aexp_arr/(fact0)
        enddo

        allocate( &
             yl_arr(num,(ly_max+1)**2,n_equals), &
             prod_arr(num,(la+1)**2,0:la,n_equals), &
             stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("LS_CACLULATE : allocation yl_arr, prod_arr failed")

        ! calculate help arrays and nuclear attraction
        call precalculate_and_nuc(this)
        if( integralpar_2cob_nuc )then
           nuc = nuc + this
        endif

        if( split3c )then
           call unpack_many_3c(this,i,OFF_V)
        end if

        ! --- further allocation ----------------------------------
        ! num : number of pairs(a,b) which are inside the cutoff
        ! for s-and r2-type there is only 1 indep. fct

           if(integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
              allocate (coul_int(symmetry_data_n_irreps()),stat=alloc_stat)
              if (alloc_stat/=0) call error_handler &
                   ("LS_CACLULATE : allocation coul_int%l failed")

              i_ir_alloc_: DO i_ir=1,symmetry_data_n_irreps()

                 n_pa = symmetry_data_n_partners(i_ir)

                 allocate (coul_int(i_ir)%l(-1:lmax_ch),stat=alloc_stat)
                 if (alloc_stat/=0) call error_handler &
                      ("LS_CACLULATE : allocation coul_int%l failed")

                 ncexps = unique_atoms(i)%r2_ch%n_exponents

                 n_independent_fcts = &
                      ua_pointer%symadapt_partner(i_ir,0)%n_independent_fcts
                 
                 allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,n_independent_fcts,&
                      nm_lb,nm_la,n_pa),stat=alloc_stat)
                 if (alloc_stat/=0) call error_handler &
                      ("LS_CACLULATE : allocation coul_int%l(-1)%m failed")

                 ncexps = unique_atoms(i)%l_ch(0)%n_exponents
                 allocate(coul_int(i_ir)%l(0)%m(num,ncexps,n_independent_fcts,&
                      nm_lb,nm_la,n_pa),stat=alloc_stat)
                 if (alloc_stat/=0) call error_handler &
                      ("LS_CACLULATE : allocation coul_int%l(0)%m failed")

                 coul_int(i_ir)%l(0)%m  = 0.0_r8_kind
                 coul_int(i_ir)%l(-1)%m = 0.0_r8_kind

                 do i_l=1,lmax_ch
                    ncexps = unique_atoms(i)%l_ch(i_l)%n_exponents
                    n_independent_fcts = &
                         ua_pointer%symadapt_partner(i_ir,i_l)%n_independent_fcts
                    allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                         nm_lb,nm_la,n_pa),stat=alloc_stat)
                    if (alloc_stat/=0) call error_handler &
                         ("LS_CACLULATE : allocation coul_int%l(0)%m failed")

                    coul_int(i_ir)%l(i_l)%m = 0.0_r8_kind

                 end do

              end do i_ir_alloc_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
           elseif (integralpar_3c_co) then
              i_ir=get_totalsymmetric_irrep()
              allocate (coul_int(i_ir),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation coul_int%l failed")

                 n_pa = 1

                 allocate (coul_int(i_ir)%l(-1:lmax_ch),stat=alloc_stat)
                 if (alloc_stat/=0) call error_handler &
                      ("LS_CACLULATE : allocation coul_int%l failed")

           ncexps = unique_atoms(i)%r2_ch%n_exponents
                 allocate(coul_int(i_ir)%l(-1)%m(num,ncexps,1,nm_lb,nm_la,n_pa),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation coul_int%l(-1)%m failed")
           ncexps = unique_atoms(i)%l_ch(0)%n_exponents
                 allocate(coul_int(i_ir)%l(0)%m(num,ncexps,1,nm_lb,nm_la,n_pa),stat=alloc_stat)
                 if (alloc_stat/=0) call error_handler &
                      ("LS_CACLULATE : allocation coul_int%l(0)%m failed")
                 coul_int(i_ir)%l(0)%m=0.0_r8_kind
                 coul_int(i_ir)%l(-1)%m=0.0_r8_kind

                 do i_l=1,lmax_ch
                    ncexps = unique_atoms(i)%l_ch(i_l)%n_exponents
                    n_independent_fcts = &
                         ua_pointer%symadapt_partner(i_ir,i_l)%n_independent_fcts
                    allocate(coul_int(i_ir)%l(i_l)%m(num,ncexps,n_independent_fcts,&
                         nm_lb,nm_la,n_pa),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation coul_int%l(0)%m failed")
                    coul_int(i_ir)%l(i_l)%m = 0.0_r8_kind
                 end do

        end if
        if(integralpar_3c_xc) then
           allocate (xc_int%l(-1:lmax_xc),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation xc_int%l failed")
           ncexps = unique_atoms(i)%r2_xc%n_exponents
           allocate(xc_int%l(-1)%m(num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation xc_int%l(-1)%m failed")
           ncexps = unique_atoms(i)%l_xc(0)%n_exponents
           allocate(xc_int%l(0)%m(num,ncexps,1,nm_lb,nm_la),stat=alloc_stat)
           if (alloc_stat/=0) call error_handler &
                ("LS_CACLULATE : allocation xc_int%l(0)%m failed")
           xc_int%l(0)%m=0.0_r8_kind
           xc_int%l(-1)%m=0.0_r8_kind
        endif


        !
           ! now calculating XC fit integrals
           !

           ! s-type xc integrals
           if(integralpar_3c_xc) call s_xc()

           ! now treating r2-type exchange integral
           if(integralpar_3c_xc) call r2_xc()

           ! l-type fit integrals
           do i_l=1,lmax_abs
              n_independent_fcts  = &
                   unique_atoms(i)%symadapt_partner(1,i_l)%n_independent_fcts

              allocate(intermediate(num,2*la+1,n_independent_fcts,n_equals,0:la),&
                   stat=alloc_stat)
              if(alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE: allocation intermediate failed")

              ! precalculate symmetry adaption and store result in intermediate
              call l_fit_symmetry_adapt()

              ! l-type exchange integrals
              if(integralpar_3c_xc.and.lmax_xc>=i_l) call l_xc()

              deallocate(intermediate,stat=alloc_stat)
              if(alloc_stat/=0) call error_handler &
                   ("LS_CALCULATE: deallocation intermediate failed")

           end do! loop over lc
           ! finished with l-type fit integrals

           !
           ! - CO calculations
        !

           if (integralpar_3c_co_resp) then    
#ifdef WITH_RESPONSE
              i_ir_: DO i_ir=1, symmetry_data_n_irreps()
                 i_pa_: DO i_pa=1,symmetry_data_n_partners(i_ir)

                    i_l=0

                    n_independent_fcts = &
                         unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

                    if (n_independent_fcts .ne. 0) then

        ! s-type coulomb
           call s_coulomb( &
                unique_atoms(i)%l_ch(0)%exponents(:), &
                            coul_int(i_ir)%l(0)%m &
                )
        
        if( split3c )then
           zexps(1) = (3.0_r8_kind/2.0_r8_kind) / unique_atoms(i)%nuclear_radius**2
                          call s_coulomb( zexps, this6d )
           if( laltlb )then
                             call unpack_many_3c(this6d(:,1,1,:,1,1),i,OFF_VFIN)
           else
                             call unpack_many_3c(this6d(:,1,1,1,:,1),i,OFF_VFIN)
           endif
        end if

        ! r2-type coloumb integral
                       call r2_coulomb()

                    end if
        ! l-type fit integrals
                    do i_l=1,lmax_ch

           n_independent_fcts  = &
                            unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

                       if (n_independent_fcts .ne. 0) then

           allocate(intermediate(num,2*la+1,n_independent_fcts,n_equals,0:la),&
                stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ("LS_CALCULATE: allocation intermediate failed")
 
           ! precalculate symmetry adaption and store result in intermediate
                          call l_fit_symmetry_adapt_v2()

           ! l-type coulomb integrals
                          call l_coulomb

                          deallocate(intermediate,stat=alloc_stat)
                          if(alloc_stat/=0) call error_handler &
                               ("LS_CALCULATE: deallocation intermediate failed")

                       end if

                    end do! loop over lc

                 end do i_pa_
              end do i_ir_
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
           elseif(integralpar_3c_co) then 

              i_ir = get_totalsymmetric_irrep()
              i_pa = 1
              i_l=0

              n_independent_fcts = &
                   unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

              if (n_independent_fcts .ne. 0) then
                 ! s-type coulomb
                 call s_coulomb( &
                      unique_atoms(i)%l_ch(0)%exponents(:), &
                      coul_int(i_ir)%l(0)%m &
                      )

                 if( split3c )then
                    zexps(1) = (3.0_r8_kind/2.0_r8_kind) / unique_atoms(i)%nuclear_radius**2
                    call s_coulomb( zexps, this6d )
                    if( laltlb )then
                       call unpack_many_3c(this6d(:,1,1,:,1,1),i,OFF_VFIN)
                    else
                       call unpack_many_3c(this6d(:,1,1,1,:,1),i,OFF_VFIN)
                    endif
                 end if

                 ! r2-type coloumb integral
                 call r2_coulomb()

              end if
              ! l-type fit integrals
              do i_l=1,lmax_ch
                 n_independent_fcts  = &
                      unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

                 if (n_independent_fcts .ne. 0) then

                    allocate(intermediate(num,2*la+1,n_independent_fcts,n_equals,0:la),&
                         stat=alloc_stat)
                    if(alloc_stat/=0) call error_handler &
                         ("LS_CALCULATE: allocation intermediate failed")

                    ! precalculate symmetry adaption and store result in intermediate
                    call l_fit_symmetry_adapt_v2()

                    ! l-type coulomb integrals
                    call l_coulomb

           deallocate(intermediate,stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ("LS_CALCULATE: deallocation intermediate failed")

                 end if

        end do! loop over lc

           end if
             

           !
           ! - End of CO calculations
           !

        ! contract the fit integrals with respect to fit dimension
        ! and write them to their final location in int_data_2cob3c_module
           if(integralpar_3c_co_resp) then
#ifdef WITH_RESPONSE
              call fitcontract_v2(num,i,cutoff,coul_int)
              do i_ir=1,symmetry_data_n_irreps()
           do i_l = -1, lmax_ch
                    deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("LS_CALCULATE : deallocation coul_int%l%m failed")
                 end do
                 deallocate (coul_int(i_ir)%l,STAT=alloc_stat)
                 if(alloc_stat.ne.0) call error_handler &
                      ("LS_CALCULATE : deallocation coul_int%l failed")
              end do
              deallocate (coul_int,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("LS_CALCULATE : deallocation coul_int%l failed")
#else
      ABORT('recompile w/ -DWITH_RESPONSE')
#endif
           elseif(integralpar_3c_co) then
              call fitcontract_v2(num,i,cutoff,coul_int)
              i_ir=get_totalsymmetric_irrep()
              do i_l = -1, lmax_ch
                 deallocate(coul_int(i_ir)%l(i_l)%m,STAT=alloc_stat)
                 if(alloc_stat.ne.0) call error_handler &
                      ("LS_CALCULATE : deallocation coul_int%l%m failed")
              end do
              deallocate (coul_int(i_ir)%l,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("LS_CALCULATE : deallocation coul_int%l failed")
              
              deallocate (coul_int,STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("LS_CALCULATE : deallocation coul_int%l failed")
        endif

        if(integralpar_3c_xc) then
           call fitcontract('xc',num,i,cutoff,xc_int)
           do i_l = -1, lmax_xc
              deallocate(xc_int%l(i_l)%m,STAT=alloc_stat)
              if(alloc_stat.ne.0) call error_handler &
                   ("LS_CALCULATE : deallocation xc_int%l%m failed")
           enddo
           deallocate (xc_int%l,STAT=alloc_stat)
           if(alloc_stat.ne.0) call error_handler &
                ("LS_CALCULATE : deallocation xc_int%l failed")
        end if


        deallocate(yl_arr,gamma_help,gamma_arg2,fact10,&
             prod_arr,stat=alloc_stat)
        if(alloc_stat/=0) call error_handler &
             ("LS_CALCULATE: deallocation yl_arr ... failed")
     else
        ua_pointer=>unique_timps(i-n_unique_atoms)
        n_equals=ua_pointer%n_equal_atoms
        z= ua_pointer%z                  ! charge 
        zc = ua_pointer%zc               ! charge
     end if
     
     if(zc/=0.0_r8_kind .and. .not. integralpar_2cob_potential) then ! pseudopotential 
        ABORT('not supported')
  end if ! pseudopotential contributions

     end do unique_atom_loop
 
     ! add contribution of point charges to nuclear attraction
     if ( (integralpar_2cob_nuc .or. integralpar_relativistic) &
          .and. (pointcharge_N +n_timps).gt. 0) then
        call  add_pointcharges()
     endif

  end if third_center_required

  DPRINT  'TIMER: ls_pseudo=',FPP_TIMER_VALUE(pls)
  if (integralpar_2cob_potential.and.old_potential) &
       call calc_potential()       !!!!!!!!!!!!!!
  
#if 0
  if (integralpar_2cob_field .and. old_elfield) &
       call calc_field()           !!!!!!!!!!!!!!
#endif

  deallocate(overlap,stat=alloc_stat)   
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation overlap failed")

  ! remap nuc to the appropriate data structure
  if (integralpar_2cob_nuc) then
    if (pseudopot_present.and.(.not.integralpar_relativistic)) then
       if(.not.laltlb) then
          do m=1,2*la+1 
             prim_int_2cob_nuc(:,:,1,m)=unpack(nuc(:,m)+&
                     nuc_pseudo(:,m),cutoff,zero)
          enddo
       else
          do m=1,2*la+1 
             prim_int_2cob_nuc(:,:,m,1)=unpack(nuc(:,m)+&
                     nuc_pseudo(:,m),cutoff,zero)
          enddo
       end if
    else !  relativistic 
       if(.not.laltlb) then
          do m=1,2*la+1
             prim_int_2cob_nuc(:,:,1,m)=unpack(nuc(:,m),cutoff,zero)
          enddo
       else
          do m=1,2*la+1
             prim_int_2cob_nuc(:,:,m,1)=unpack(nuc(:,m),cutoff,zero)
          enddo
       endif
       if(pseudopot_present)  then
          if(.not.laltlb) then
             do m=1,2*la+1
                prim_int_2cob_nuc_pseudo(:,:,1,m)=unpack(nuc_pseudo(:,m),cutoff,zero)
             enddo
          else
             do m=1,2*la+1
                prim_int_2cob_nuc_pseudo(:,:,m,1)=unpack(nuc_pseudo(:,m),cutoff,zero)
             enddo
          endif

          with_timps = pointcharge_N +n_timps .gt. 0
          if( with_timps ) then
             ! add point charge and TIMP to PP contributions
             if(.not.laltlb) then
                do m=1,2*la+1
                   prim_int_2cob_nuc_pseudo(:,:,1,m)=prim_int_2cob_nuc_pseudo(:,:,1,m) &
                        + unpack(nuc_pc_timps(:,m),cutoff,zero)
                enddo
             else
                do m=1,2*la+1
                   prim_int_2cob_nuc_pseudo(:,:,m,1)=prim_int_2cob_nuc_pseudo(:,:,m,1) &
                        + unpack(nuc_pc_timps(:,m),cutoff,zero)
                enddo
             endif
          endif
       endif ! pseudopot_present
    end if ! else of PP and not rel

     deallocate(nuc,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation nuc failed")
    if (pseudopot_present) then
        deallocate(nuc_pseudo,stat=alloc_stat)
        if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation nuc_pseudo failed")
    end if
!!!    call print_nuclear()
  end if

  if (integralpar_2cob_potential.and.old_potential) then
     if(.not.laltlb) then              
        do i=1,N_points                
           do m=1,2*la+1               
              !restoring needed for VPP
              intermed_3c(:)=potential(i,:,m)
              prim_int_2cob_poten(:,:,i,1,m)=unpack(intermed_3c(:),cutoff,zero) 
           enddo
        enddo
     else
        do i=1,N_points
           do m=1,2*la+1
              intermed_3c(:)=potential(i,:,m)
              prim_int_2cob_poten(:,:,i,m,1)=unpack(intermed_3c(:),cutoff,zero) 
           enddo
        enddo
     end if
     deallocate(potential,stat=alloc_stat)        !!!!!!!!!!!!!!!
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation potential failed")
     deallocate(intermed_3c,stat=alloc_stat) 
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation intermed_3c failed")
  endif

#if 0
  if (integralpar_2cob_field .and. old_elfield) then
     if(calc_normal) then
        N_length=N_surface_points
     else
        N_length=totsym_field_length
     end if
     if(.not.laltlb) then              
        do i=1,N_length                
           do m=1,2*la+1               
              intermed_3c(:)=field(i,:,m)
              prim_int_2cob_field(:,:,i,1,m)=unpack(intermed_3c,cutoff,zero) 
           enddo
        enddo
     else
        do i=1,N_length
           do m=1,2*la+1
              intermed_3c(:)=field(i,:,m)
              prim_int_2cob_field(:,:,i,m,1)=unpack(intermed_3c,cutoff,zero)
           enddo
        enddo
     end if
     deallocate(field,stat=alloc_stat)           
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation field failed")
     deallocate(intermed_3c,stat=alloc_stat)           
     if(alloc_stat/=0) call error_handler &
          ("LS_CALCULATE: deallocation intermed_3c for field failed")
  endif
#endif
  ! remap pvscalarp to the appropriate data structure
  if (integralpar_relativistic) then
     with_timps = pointcharge_N +n_timps .gt. 0
     if( with_timps .and. pseudopot_present ) then
        deallocate(nuc_pc_timps, stat=alloc_stat)
        if (alloc_stat.ne.0) call error_handler &
             ("LS_CALCULATE: deallocation(nuc_pc_timps) failed")
     endif ! 
  end if ! integralpar_relativistic

  deallocate(this,stat=memstat)
  ASSERT(memstat==0)
  if( split3c )then
     deallocate(this6d,stat=memstat)
     ASSERT(memstat==0)
  endif

  deallocate(clmamb,clmamb_scalar,fact0,fact1,fact2,fact4,fact6,&
       fact7,fact8,rcsabc,gamma_arg,diff_arr,aexp_arr,cutoff,stat=alloc_stat)
  if(alloc_stat/=0) call error_handler &
       ("LS_CALCULATE: deallocation clmamb ... failed")

contains

  subroutine unpack_many_3c(this,i,offset)
    implicit none
    real(r8_kind), intent(in)    :: this(:,:)
    integer(i4_kind), intent(in) :: i,offset
    ! *** end of iterface ***

    integer(i4_kind) :: m, pos

    pos = OFF_STRIDE*(i-1)+offset
    ASSERT(pos<=size(many_3c,3))
    ASSERT(size(this,2)==2*la+1)
    ! remap ints to the appropriate data structure
    if(.not.laltlb) then
       do m=1,2*la+1 
          many_3c(:,:,pos,1,m)=unpack(this(:,m),cutoff,zero)
       enddo
    else
       do m=1,2*la+1 
          many_3c(:,:,pos,m,1)=unpack(this(:,m),cutoff,zero)
       enddo
    endif
  end subroutine unpack_many_3c

  !**************************************************************

  subroutine  precalculate_and_nuc(nuc)
    ! precalculate prod_arr and calculate nuclear attraction
    implicit none
    real(r8_kind), intent(out) :: nuc(:,:) ! (num,2*la+1)
    ! *** end of iterface ***

    nuc = zero

    equal_atoms_nuc: do j=1,n_equals 
       xc=unique_atoms(i)%position(:,j)
       ! xc(3)=xc(3)+(i_l-1)*0.001_r8_kind
       yl_arr(:,:,j)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
            spread(xc,1,num))
       gamma_arg2(:,j)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
            (gamma_arg(:,3)-xc(3))**2)*fact0
       ! first calculation of the nuclear attraction
       gamma_help(:,1:3+la,j)=gamma(3+la,gamma_arg2(:,j))
       prod_arr(:,:,:,j)=prod_rule2(yl_arr(:,:,j),clmamb,la)
       if (integralpar_2cob_nuc) then           
          counter=1
          do l=0,la
             if(zc.ne.zero.and.integralpar_relativistic) then
               nuc_pseudo(:,:)=nuc_pseudo(:,:) &
                    +(z-zc)*prod_arr(:,la**2+1:(la+1)**2,l,j)*&
                  spread(gamma_help(:,l+1,j)*&
                  (-2.0_r8_kind*aexp_arr)**l*fact2,2,2*la+1) 
               else
             nuc(:,:)=nuc(:,:)+(z-zc)*prod_arr(:,la**2+1:(la+1)**2,l,j)*&
                  spread(gamma_help(:,l+1,j)*&
                  (-2.0_r8_kind*aexp_arr)**l*fact2,2,2*la+1)
          end if
          
          enddo
          ! calculation of nuclear attraction finished
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
    integer(i4_kind)   :: n_cf, i_ind, n_if, ic
    real(r8_kind)      :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if  =  unique_atoms(i)%symadapt_partner(i_ir,0)%n_independent_fcts

    do i_ind = 1, n_if
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
             fact8=2.0_r8_kind*fact6*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))     
             gamma_help(:,1:1+la,1)=gamma(1+la,gamma_arg2(:,j)*rcsabc(:))
             do l=0,la
                help_vec=gamma_help(:,l+1,1)*fact8
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      coul(:,k,i_ind,i_lb,i_la,i_pa)=&
                           coul(:,k,i_ind,i_lb,i_la,i_pa)+&
                           ( prod_arr(:,la**2-1+i_la+i_lb,l,j)*&
                           help_vec ) * NORM * sc(j)
                   end do
                end do
                fact8=fact8*(-2.0_r8_kind*aexp_arr*rcsabc(:))
             enddo! loop over l
          end do! loop over fitfunctions
       end do equal_atoms_s_ch
    end do
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
          fact8=fact6*sqrt(fact0/(fact0+cexps(k)))&
               *(fact0/(fact0+cexps(k)))*exp(-rcsabc*gamma_arg2(:,j))
          ! now the acutal calculation starts
          do l=0,la
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   xc_int%l(0)%m(:,k,1,i_lb,i_la)=&
                        xc_int%l(0)%m(:,k,1,i_lb,i_la)+&
                        prod_arr(:,la**2-1+i_la+i_lb,l,j)*&
                        fact8*(-2.0_r8_kind*aexp_arr*rcsabc(:))**l
                end do
             end do
          end do
       end do! loop over fitfunctions
    end do equal_atoms_s_xc
  end subroutine s_xc

  !**************************************************************

  subroutine r2_coulomb()
    real(kind=r8_kind) :: help_vec(num),help_vec2(num)
    real(kind=r8_kind),pointer,dimension(:,:,:,:,:,:) :: pointer_coul
    ! r2-type coloumb fit integral
    integer(i4_kind)   :: n_cf, n_if, i_ind, ic
    real(r8_kind)      :: NORM, sc(n_equals)

    NORM  =  SQRT(REAL(n_equals,r8_kind))
    n_if  =  unique_atoms(i)%symadapt_partner(i_ir,i_l)%n_independent_fcts

    do i_ind = 1, n_if

       n_cf    =  unique_atoms(i)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%n_fcts     
       coeff   => unique_atoms(i)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%c
       eq_atom => unique_atoms(i)%symadapt_partner(i_ir,i_l)%symadapt(i_ind,i_pa)%I_equal_atom

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
             fact8=2.0_r8_kind*fact6*pi/(cexps(k)**2)*sqrt(fact0/(fact0+cexps(k)))&
                  *(fact0/(fact0+cexps(k)))
             fact7=(fact0+1.5_r8_kind*cexps(k))/fact0
             gamma_help(:,1:2+la,1)=gamma(2+la,gamma_arg2(:,j)*rcsabc(:))
             help_vec2=gamma_arg2(:,j)*rcsabc(:)
             do l=0,la
                help_vec=fact8*((fact7-real(l,r8_kind))*gamma_help(:,l+1,1)+&
                     help_vec2*gamma_help(:,l+2,1))
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      pointer_coul(:,k,i_ind,i_lb,i_la,i_pa)=&
                           pointer_coul(:,k,i_ind,i_lb,i_la,i_pa)+&
                           ( prod_arr(:,la**2-1+i_la+i_lb,l,j)*&
                           help_vec ) * NORM * sc(j)
                   end do
                end do
                fact8=fact8*(-2.0_r8_kind*aexp_arr*rcsabc(:))
             enddo! loop over l
          end do! loop over fitfunctions
       end do equal_atoms_r2_ch
    end do
  end subroutine r2_coulomb

  !**************************************************************

  subroutine r2_xc()
    ! r2-type xc fit integral
    ncexps = unique_atoms(i)%r2_xc%n_exponents
    cexps => unique_atoms(i)%r2_xc%exponents(:)
    equal_atoms_r2_xc: do j=1,n_equals
       call integral_interrupt_2cob3c()
       do k=1,ncexps ! loop over fitexponents
          ! precalculation of two factors
          rcsabc=cexps(k)/(fact0+cexps(k)) ! c/(a+b+c)
          fact8=fact6*sqrt(fact0/(fact0+cexps(k)))*fact0/(fact0+cexps(k))&
               /(fact0+cexps(k))*exp(-gamma_arg2(:,j)*rcsabc(:))
          do l=0,la            
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   xc_int%l(-1)%m(:,k,1,i_lb,i_la)=&
                        xc_int%l(-1)%m(:,k,1,i_lb,i_la)+&
                        prod_arr(:,la**2-1+i_la+i_lb,l,j)*&
                        fact8*&
                        (-2.0_r8_kind*aexp_arr*rcsabc(:))**l*&
                        (1.5_r8_kind+fact0/cexps(k)*&
                        (gamma_arg2(:,j)*rcsabc(:)-l))
                end do
             end do
          enddo! loop over l
       end do! loop over fitfunctions
    end do equal_atoms_r2_xc
  end subroutine r2_xc

  !**************************************************************

  subroutine l_fit_symmetry_adapt_v2() 
    ! eq(3.63) from M. Staufer rebuild for symmetry
    ! symmetry adaption for l-type charge and xc fit integrals
    ! for a given l -> i_l, i_ir and i_pa,
    ! results are stored in intermediate array intermediate(num,ma,mb,i_ind,n_equal,0:la+lb)
    real(r8_kind) :: NORM

    intermediate=0.0_r8_kind
    NORM = 1.0_r8_kind

    independents: do i_ind=1,n_independent_fcts
       ! i_ir and i_pa was added because of the symmetrization of fit_fct
       n_contributing_fcts = &
            unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%n_fcts
       coeff   => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%c
       magn    => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%m
       eq_atom => unique_atoms(i)%symadapt_partner(i_ir,i_l)%&
            symadapt(i_ind,i_pa)%I_equal_atom 
       contributing: do i_cnt=1,n_contributing_fcts
          diff_arr(:,:)=spread((aexp_arr/fact0)**i_l,2,(la+1)**2)*&
               diff_rule(yl_arr(:,:,eq_atom(i_cnt))/&
               fact10(:,1:((ly_max+1)**2)),1,(la+1)**2,i_l**2+&
               magn(i_cnt))
          do l=0,la
             intermediate(:,:,i_ind,eq_atom(i_cnt),l)=&
                  intermediate(:,:,i_ind,eq_atom(i_cnt),l)+&
                  NORM * coeff(i_cnt)*prod_rule(prod_arr(:,:,l,eq_atom(i_cnt)),&
                  diff_arr(:,:),la**2+1,(la+1)**2)
          end do
       end do contributing
    end do independents
  end subroutine l_fit_symmetry_adapt_v2


  subroutine l_fit_symmetry_adapt
    ! precalculate symmetry adaption and store result in intermediate
    ! for a given i_l of l-type xc and coulomb fit integral
    intermediate=0.0_r8_kind
    independents: do i_ind=1,n_independent_fcts
       n_contributing_fcts = &
            unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%n_fcts
       coeff => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%c
       magn => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%m
       eq_atom => unique_atoms(i)%symadapt_partner(1,i_l)%symadapt(i_ind,1)%&
            I_equal_atom 
       contributing: do i_cnt=1,n_contributing_fcts
          diff_arr(:,:)=spread((aexp_arr/fact0)**i_l,2,(la+1)**2)*&
               diff_rule(yl_arr(:,:,eq_atom(i_cnt))/&
               fact10(:,1:((ly_max+1)**2)),1,(la+1)**2,i_l**2+&
               magn(i_cnt))
          do l=0,la
             intermediate(:,:,i_ind,eq_atom(i_cnt),l)=&
                  intermediate(:,:,i_ind,eq_atom(i_cnt),l)+&
                  coeff(i_cnt)*prod_rule(prod_arr(:,:,l,eq_atom(i_cnt)),&
                  diff_arr(:,:),la**2+1,(la+1)**2)
          end do
       end do contributing
    end do independents
  end subroutine l_fit_symmetry_adapt

  !**************************************************************

  subroutine l_coulomb()
    real(kind=r8_kind) :: help_vec(num),help_vec2(num)
    real(kind=r8_kind),pointer,dimension(:,:,:,:,:,:) :: pointer_coul
    ! l-type coulomb fit integrals for given i_l
    ncexps=  unique_atoms(i)%l_ch(i_l)%n_exponents
    cexps => unique_atoms(i)%l_ch(i_l)%exponents(:)
    ! allocation for coul_int
    pointer_coul=>coul_int(i_ir)%l(i_l)%m

    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       equals_ch:  do j=1,n_equals
          rcsabc=cexps(k)/(fact0+cexps(k))
          gamma_help(:,1:max_order,1)=gamma(max_order,gamma_arg2(:,j)*rcsabc(:))
          fact8=2.0_r8_kind*fact6*pi/cexps(k)*sqrt(fact0/(fact0+cexps(k)))
          help_vec=(2.0_r8_kind*fact0*rcsabc)**i_l*fact8
          do l=0,la
             help_vec2=help_vec*gamma_help(:,l+1+i_l,1)
             do i_ind=1,n_independent_fcts
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      pointer_coul(:,k,i_ind,i_lb,i_la,i_pa)=&
                           pointer_coul(:,k,i_ind,i_lb,i_la,i_pa)+&
                           intermediate(:,i_la+i_lb-1,i_ind,j,l)*&
                           help_vec2
                   end do
                end do
             end do! loop over i_ind
             help_vec=help_vec*(-2.0_r8_kind*aexp_arr*rcsabc(:))
          end do! loop over l
       end do equals_ch
    end do! loop over k 
  end subroutine l_coulomb

  !**************************************************************

  subroutine l_xc()
    ! l-type exchange fit integrals for given i_l
    cexps => unique_atoms(i)%l_xc(i_l)%exponents(:)
    ncexps = unique_atoms(i)%l_xc(i_l)%n_exponents
    allocate(xc_int%l(i_l)%m(num,ncexps,n_independent_fcts,nm_lb,nm_la),&
         stat=alloc_stat)
    if(alloc_stat/=0) call error_handler &
         ("LS_CALCULATE: allocation xc_int%l(i_l)%m failed")  
    xc_int%l(i_l)%m=0.0_r8_kind
    call integral_interrupt_2cob3c()
    do k=1,ncexps ! loop over fitexponents
       equals_xc:  do j=1,n_equals
          rcsabc=cexps(k)/(fact0+cexps(k))
          fact8=fact6*sqrt(fact0/(fact0+cexps(k)))*fact0/&
               (fact0+cexps(k))
          gamma_help(:,1,1)=exp(-gamma_arg2(:,j)*rcsabc(:))
          counter=1
          do i_ind=1,n_independent_fcts
             do l=0,la
                do i_la=1,nm_la
                   do i_lb=1,nm_lb
                      xc_int%l(i_l)%m(:,k,i_ind,i_lb,i_la)=&
                           xc_int%l(i_l)%m&
                           (:,k,i_ind,i_lb,i_la)+&
                           intermediate(:,i_la+i_lb-1,i_ind,j,l)*&
                           (2.0_r8_kind*fact0*rcsabc)**i_l*&
                           (-2.0_r8_kind*aexp_arr*rcsabc(:))**l &
                           *fact8*gamma_help(:,1,1)
                   end do
                end do
             end do! loop over l
          end do! loop over i_ind   
       end do equals_xc
    end do! loop over k 
  end subroutine l_xc

  !**************************************************************

  subroutine add_pointcharges
    ! add contribution of point charges to nuclear attraction
    unique_charge_loop: do i=1,pointcharge_N+n_timps  ! loop over third center
       if(i<=n_timps) then
          z= unique_timps(i)%z - unique_timps(i)%zc
          n_equals= unique_timps(i)%n_equal_atoms
       else
          cycle unique_charge_loop ! pointcharges go to SHGI !!!!!!!!!!!!!!AS
          z= pointcharge_array(i-n_timps)%z                  ! charge 
          n_equals=pointcharge_array(i-n_timps)%n_equal_charges
       end if
       ly_max = la
       max_order = la + 3
       allocate ( &
            gamma_help(num,max_order,n_equals), &
            gamma_arg2(num,n_equals),&
            yl_arr(num,(ly_max+1)**2,n_equals), &
            prod_arr(num,(la+1)**2,0:la,n_equals),&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LS_CACLULATE : allocation add_pointcharges failed")
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
          gamma_help(:,1:3+la,j)=gamma(3+la,gamma_arg2(:,j))
          prod_arr(:,:,:,j)=prod_rule2(yl_arr(:,:,j),clmamb,la)
          counter=1
          do l=0,la
	if(integralpar_relativistic.and.pseudopot_present) then
	 nuc_pc_timps(:,:)=nuc_pc_timps(:,:)&
              +z*prod_arr(:,la**2+1:(la+1)**2,l,j)*&
                  spread(gamma_help(:,l+1,j)*&
                  (-2.0_r8_kind*aexp_arr)**l*fact2,2,2*la+1)
	else
             nuc(:,:)=nuc(:,:)+z*prod_arr(:,la**2+1:(la+1)**2,l,j)*&
                  spread(gamma_help(:,l+1,j)*&
                  (-2.0_r8_kind*aexp_arr)**l*fact2,2,2*la+1)
          endif! integralpar_relativistic.and.pseudopot_present/else
          enddo
       end do equal_charge_nuc

       deallocate(yl_arr,gamma_help,gamma_arg2,prod_arr,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler &
            ("LS_CALCULATE: deallocation add_pointcharges failed")
    end do unique_charge_loop
  end subroutine add_pointcharges
  !**************************************************************

  !**************************************************************
  subroutine calc_potential
    ! integrals of potential

    integer(kind=i4_kind) ::i_ps

    do i_ps=1,N_points
       call integral_interrupt_2cob3c()
       n_equals=point_in_space(i_ps)%N_equal_points
       ly_max = la
       max_order = la + 3
       allocate ( &
            gamma_help(num,max_order,n_equals), &
            gamma_arg2(num,n_equals),&
            yl_arr(num,(ly_max+1)**2,n_equals), &
            prod_arr(num,(la+1)**2,0:la,n_equals),&
            stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LS_CACLULATE : allocation calc_potential failed")
       equal_points: do j=1,n_equals
          xc=point_in_space(i_ps)%position(:,j)
          yl_arr(:,:,j)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
               spread(xc,1,num))
          gamma_arg2(:,j)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2+&
               (gamma_arg(:,3)-xc(3))**2)*fact0
          gamma_help(:,1:3+la,j)=gamma(3+la,gamma_arg2(:,j))
          prod_arr(:,:,:,j)=prod_rule2(yl_arr(:,:,j),clmamb,la)
          counter=1
          do l=0,la
             potential(i_ps,:,:)=potential(i_ps,:,:)+prod_arr(:,la**2+1:(la+1)**2,l,j)*&
                  spread(gamma_help(:,l+1,j)*&
                  (-2.0_r8_kind*aexp_arr)**l*fact2,2,2*la+1)
          enddo
       end do equal_points

       deallocate(yl_arr,gamma_help,gamma_arg2,prod_arr,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler &
            ("LS_CALCULATE: deallocation of calc_potential failed")
    enddo
  end subroutine calc_potential
  !**************************************************************

  !**************************************************************
#if 0
  subroutine calc_field

    integer(kind=i4_kind) ::i,l,k,i_la,i_lb,i_grad,counter,j,f_dim,m,n_equals
    real(kind=r8_kind),allocatable :: poten_grad(:,:,:),grad_poten(:,:,:) 
    real(kind=r8_kind),pointer :: rotmat(:,:)
    logical :: do_rotation

    allocate ( grad_poten(num,3,2*la+1), &
         stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("LS_CALCULATE : allocation grad_poten failed")

    ly_max=la
    max_order=4+la

    allocate(overlap1(num,(la+1)**2),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE : allocation overlap1 failed")
    overlap1=zero

    do l=0,la
       do m=1,2*l+1
          overlap1(:,l**2+m)=fact6*clmamb(:,l**2+m)
       enddo
    enddo

    allocate (yl_arr(num,(ly_max+1)**2,1), &
         gamma_help(num,max_order,1), &
         gamma_arg2(num,1),diff_arr_grad(num,(ly_max+1)**2,3),&
         prod_arr(num,(la+1)**2,0:la,1), prod_arr_gr(num,3,(la+1)**2,0:la+1), &
         help_arr_gr(num,3,nm_lb,nm_la),rsaba(num),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE : allocation (1) failed")

    rsaba=fact0/aexp_arr

    Nsp: do i=1,N_surface_points
       n_equals=surface_points(i)%N_equal_points
       f_dim=surf_points_grad_index(i+1)-surf_points_grad_index(i)

       allocate(poten_grad(num,2*la+1,f_dim),stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LS_CALCULATE : allocation poten_grad failed")
       poten_grad=0.0_r8_kind

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

          yl_arr(:,:,1)=solid_harmonics_calc(ly_max,gamma_arg(:,:)-&
               spread(xc,1,num))
          gamma_arg2(:,1)=((gamma_arg(:,1)-xc(1))**2+(gamma_arg(:,2)-xc(2))**2&
               +(gamma_arg(:,3)-xc(3))**2)*fact0
          diff_arr_grad=0.0_r8_kind
          do l=1,(ly_max+1)**2
             do k=1,solhrules_differential(3,l)%n_summands
                diff_arr_grad(:,l,1)=diff_arr_grad(:,l,1)+&
                     solhrules_differential(3,l)%coef(k)*&
                     yl_arr(:,solhrules_differential(3,l)%lm_sh(k),1)  
             end do
             do k=1,solhrules_differential(4,l)%n_summands
                diff_arr_grad(:,l,2)=diff_arr_grad(:,l,2)+&
                     solhrules_differential(4,l)%coef(k)*&
                     yl_arr(:,solhrules_differential(4,l)%lm_sh(k),1)
             end do
             do k=1,solhrules_differential(2,l)%n_summands
                diff_arr_grad(:,l,3)=diff_arr_grad(:,l,3)+&
                     solhrules_differential(2,l)%coef(k)*&
                     yl_arr(:,solhrules_differential(2,l)%lm_sh(k),1)
             end do
          end do

          call calculate_helpers()

          if(do_rotation) then
             counter=1
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   do i_grad=1,f_dim 
                      poten_grad(:,counter,i_grad)=poten_grad(:,counter,i_grad)-&
                           rotmat(i_grad,1)*help_arr_gr(:,1,i_lb,i_la)-&
                           rotmat(i_grad,2)*help_arr_gr(:,2,i_lb,i_la)-&
                           rotmat(i_grad,3)*help_arr_gr(:,3,i_lb,i_la)
                   enddo
                   counter=counter+1
                end do
             enddo
          else
             counter=1
             do i_la=1,nm_la
                do i_lb=1,nm_lb
                   do i_grad=1,3
                      poten_grad(:,counter,i_grad)=poten_grad(:,counter,i_grad)-&
                           help_arr_gr(:,i_grad,i_lb,i_la)
                   enddo
                   counter=counter+1
                end do
             enddo
          end if
       enddo n_equals_j

       if(calc_normal) then
          grad_poten=0.0_r8_kind
          rotmat=>surf_points_grad_info(i)%m(:,:,1)
          do j=1,f_dim
             do k=1,3
                grad_poten(:,k,:)= grad_poten(:,k,:)+rotmat(j,k)*poten_grad(:,:,j)
             enddo
          enddo

          do k=1,3
             field(i,:,:)=field(i,:,:)+grad_poten(:,k,:)*surface_points(i)%out_normal(k)
          enddo
       else
          k=surf_points_grad_index(i)
          do j=1,f_dim
             field(k,:,:)=poten_grad(:,:,j)
             k=k+1
          end do
       end if

       deallocate (poten_grad,stat=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("LS_CALCULATE : deallocationi poten_grad failed")

    enddo Nsp

!!$print*,'===========FIELD-LS========'
!!$print*,field(1,1,1)
!!$print*,'========================'

    deallocate(overlap1,stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
             ("LS_CALCULATE : deallocation overlap1 failed")

    deallocate ( yl_arr, gamma_help, &
         gamma_arg2,diff_arr_grad, prod_arr, prod_arr_gr, &
         help_arr_gr,rsaba,stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("LS_CALCULATE : deallocation (1) failed")
    
    deallocate (grad_poten,stat=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("LS_CALCULATE : deallocationi grad_poten failed")

  end subroutine calc_field
#endif
  !************************************************************

  !************************************************************
  subroutine calculate_helpers
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind) :: help_vec_arr(num,3),help_vec1(num),help_vec2(num),help_vec(num),fact21(num)
    integer(kind=i4_kind) :: ii,jj,l_i,l_index
    integer(kind=i4_kind) :: i1,i2,k_gr
    integer(kind=i4_kind), pointer :: index1_p(:),index2_p(:)
    real(kind=r8_kind) :: coeff
    real(kind=r8_kind), pointer :: coeff_p(:)
    !------------ Executable code ------------------------------

    gamma_help(:,1:4+la,1)=gamma(4+la,gamma_arg2(:,1))
    prod_arr=0.0_r8_kind
    prod_arr_gr=0.0_r8_kind

    do ii=1,(la+1)**2
       index1_p=>solhrules_product(ii)%lm_sh1
       index2_p=>solhrules_product(ii)%lm_sh2
       coeff_p=>solhrules_product(ii)%coef
       do jj=1,solhrules_product(ii)%n_summands
          i1=index1_p(jj)
          l_i=solhrules_l_and_m_of_lm(1,i1)
          i2=index2_p(jj)
          coeff=coeff_p(jj)
          help_vec1=coeff*yl_arr(:,i1,1)
          prod_arr(:,ii,l_i,1)=prod_arr(:,ii,l_i,1)+help_vec1*overlap1(:,i2)
          help_vec2=coeff*overlap1(:,i2)
          do k_gr=1,3
             prod_arr_gr(:,k_gr,ii,l_i)=prod_arr_gr(:,k_gr,ii,l_i) &  !!!!!!k_gr+3=>k_gr
                  +diff_arr_grad(:,i1,k_gr)*help_vec2
          end do
       end do
    enddo
    do l=0,la
       do k=1,(la+1)**2
          do k_gr=1,3
             prod_arr_gr(:,k_gr,k,l+1) = prod_arr_gr(:,k_gr,k,l+1)+&  !!!!!!k_gr+3=>k_gr
                  prod_arr(:,k,l,1)*(gamma_arg(:,k_gr)-xc(k_gr))*rsaba
          end do
       end do
    end do

    fact21=2.0_r8_kind*sqrt(fact0/pi)
    help_arr_gr=0.0_r8_kind
    do l=0,la+1
       help_vec=gamma_help(:,l+1,1)*&
            (-2.0_r8_kind*aexp_arr)**l*fact21
       do k_gr=1,3 !!!!!6
          help_vec_arr(:,k_gr)=help_vec
       end do
       do i_la=1,nm_la
          do i_lb=1,nm_lb
             l_index=la**2-1+i_la+i_lb 
             help_arr_gr(:,:,i_lb,i_la)=help_arr_gr(:,:,i_lb,i_la) &
                  +prod_arr_gr(:,:,l_index,l)*help_vec_arr
          enddo
       end do
    enddo
  end subroutine calculate_helpers
  !**************************************************************

#if 0
  !**************************************************************!
  ! Just to print nuclear integrals (duplicated) !
  subroutine print_nuclear()
    use iounitadmin_module
    integer(kind=i4_kind) :: i_exp1, i_exp2, nuc_int_unit, ne1,ne2
    logical               :: file_exists
    
    nuc_int_unit = get_iounit()
    inquire(file="NUC_INT", exist=file_exists)
    if(file_exists) then
       open(nuc_int_unit, file="NUC_INT",position="append")
    else
       open(nuc_int_unit, file="NUC_INT")
    end if  

    write(nuc_int_unit,"('********* LS ********** [',2i2,' |',2i2,' ] ********* LS ***********')")  na, la_in, nb, lb
    if(laltlb) then
       ne1=naexps
       ne2=nbexps
    else
       ne1=nbexps
       ne2=naexps
    end if
    do i_exp1 = 1,ne1
       do i_exp2 = 1,ne2
          do ma=1,2*la+1
             if(.not.laltlb) then
                write(nuc_int_unit,"(5x,'exp2 = ',i3,3x,'exp1 = ',i3,3x, 'm',2i3,4x,f14.8)") &
                     i_exp2, i_exp1, 1, ma, prim_int_2cob_nuc(i_exp2,i_exp1,1,ma)
             else
                write(nuc_int_unit,"(5x,'exp2 = ',i3,3x,'exp1 = ',i3,3x, 'm',2i3,4x,f14.8)") &
                     i_exp2, i_exp1, ma, 1, prim_int_2cob_nuc(i_exp2,i_exp1,ma,1)               
             end if
          end do
       end do
    end do
    close(nuc_int_unit, status="keep")
    call return_iounit(nuc_int_unit)

  end subroutine print_nuclear
  !**************************************************************!
#endif


end subroutine ls_calculate

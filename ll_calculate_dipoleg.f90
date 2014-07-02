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
!       
! Public interface of module
!=====================================================================
subroutine ll_calculate_dipoleg(na,nb,la,lb)
  !
  !  Purpose: calculation of all primitive dipole integrals
  !       for a given set of indizes
  !       (unique_atom1,unique_atom2,la,lb).
  !
  !  Author: MS
  !  Date:   8/96
  !
#include <def.h>

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
  integer(kind=i4_kind) :: num,m,alloc_stat,i_xyz,k,i_l,i_lma,i_lmb, &
       n_lma,n_lmb,naexps,nbexps,ma,mb,l,iota
  integer(kind=i4_kind), dimension(2:4) :: lm2xyz = (/3,1,2/)
! transforms lm-indices 2(z), 3(x), 4(y) into 1(x), 2(y), 3(z)

  logical,allocatable   :: cutoff(:,:)
  real(kind=r8_kind),pointer,dimension(:) :: &
       aexps,bexps
  real(kind=r8_kind),allocatable,dimension(:) :: &
       fact0,fact1,fact2,fact4,fact6,tau,aexp_arr,bexp_arr,clmamb_scalar
  real(kind=r8_kind),allocatable,dimension(:,:) :: &
       fact0_arr,fact1_arr,fact2_arr,clmamb,diff_arr0, &
       grad_clmamb_scalar(:,:)
  
  real(kind=r8_kind),allocatable,dimension(:,:,:) :: &
       overlap,diff_arr0_grad,clmamb_grad
  real(kind=r8_kind),allocatable,dimension(:,:,:,:) :: &
       diff_rule_result,overlap_grad,orbzeeman, corrzeeman
  real(kind=r8_kind),allocatable,dimension(:,:,:,:,:) :: &
        diff_rule_result_grad,rnabla_primitive
  real(kind=r8_kind),dimension(3)  :: xa,xb,xd
  real(kind=r8_kind) :: arg
  integer(kind=i4_kind),parameter:: i1=1,i2=2,i4=4

  naexps = unique_atoms(na)%l_ob(la)%n_exponents
  nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents

  allocate( &
       fact0_arr(nbexps,naexps), &
       fact1_arr(nbexps,naexps), &
       fact2_arr(nbexps,naexps), &
       cutoff(nbexps,naexps), &
       stat=alloc_stat)

  xa = center1
  xb = center2
!!$  xb(2) = xb(2)+0.0001
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
     prim_int_2cob_dipoleg = 0.0_r8_kind
     deallocate( &
          fact0_arr, &
          fact1_arr, &
          fact2_arr, &
          cutoff, &
          stat=alloc_stat)
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
       overlap_grad(num,n_lma,n_lmb,2:4), &
       orbzeeman(num,2*la+1,2*lb+1,3), &
       rnabla_primitive(num,2*la+1,2*lb+1,3,3), &
       corrzeeman(num,2*la+1,2*lb+1,3), &
       diff_rule_result(num,n_lma,n_lmb,3), &
       diff_rule_result_grad(num,n_lma,n_lmb,3,2:4), &
       aexp_arr(num), &
       bexp_arr(num), &
       clmamb_scalar((max(la,lb)+1)**2), &
       grad_clmamb_scalar(2:4,(la+1)**2), &
       clmamb(num,(la+1)**2), &
       clmamb_grad(num,(la+1)**2,2:4), &
       diff_arr0((la+1)**2,(lb+1)**2), &
       diff_arr0_grad((la+1)**2,(lb+1)**2,2:4), &
       stat=alloc_stat)


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

  ! precalculation of solid harmonics
  clmamb_scalar=solid_harmonics_scalar(max(la,lb),xd)
  tau=fact2*arg ! a*b/(a+b)*(A-B)**2

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
        do i_xyz=2,4
        diff_arr0_grad(:,i_lmb,i_xyz)=&
         reshape( &
          diff_rule_nested(spread(clmamb_scalar,1,1), 1, n_lma, i_xyz,i_lmb), &
                (/n_lma/) ) 
        enddo
!       do i_lma=1,n_lma
!        print*,i_lma,clmamb_scalar(i_lmb),diff_arr0(i_lma,i_lmb),&
!               diff_arr0_grad(i_lma,i_lmb,2)
!       enddo 
        i_lmb=i_lmb+1
     end do
  end do
   do i_lma=1,n_lma
     grad_clmamb_scalar(2:4,i_lma)=&
     reshape(diff_rule(spread(clmamb_scalar,1,1),i2,i4,i_lma),(/3/))
!print*,i_lma,clmamb_scalar(i_lma),grad_clmamb_scalar(2,i_lma)  
   enddo

  fact4=1.0_r8_kind
  i_lma=1
  do l=0,la
     do m=1,2*l+1
        clmamb(:,i_lma)=clmamb_scalar(i_lma)*fact4
        clmamb_grad(:,i_lma,2)=grad_clmamb_scalar(2,i_lma)*fact4
        clmamb_grad(:,i_lma,3)=grad_clmamb_scalar(3,i_lma)*fact4
        clmamb_grad(:,i_lma,4)=grad_clmamb_scalar(4,i_lma)*fact4
        i_lma=i_lma+1
     enddo
     fact4=-fact4*fact2*2.0_r8_kind
  enddo

  i_lmb=1
  do i_l=0,lb
     do mb=1,2*i_l+1
        overlap(:,:,i_lmb)= &
             spread( fact6*(2.0_r8_kind*fact2)**i_l, 2, n_lma ) * &
             prod_rule( &
             spread( diff_arr0(:,i_lmb), 1, num ), &
             clmamb(:,:), 1, n_lma &
        )

        do iota=2,4
        overlap_grad(:,:,i_lmb,iota)= &
         spread( -xd(lm2xyz(iota))*fact6*(2.0_r8_kind*fact2)**(i_l+1), 2, n_lma )  *  &
         prod_rule( &
             spread( diff_arr0(:,i_lmb), 1, num ), &
             clmamb(:,:), 1, n_lma &
          )+ &
         spread( fact6*(2.0_r8_kind*fact2)**i_l, 2, n_lma ) * ( &
         prod_rule( &
             spread( diff_arr0(:,i_lmb),1, num ), &
             clmamb_grad(:,:,iota), i1, n_lma &
             )+ &
        prod_rule( &
             spread(diff_arr0_grad(:,i_lmb,iota), 1, num ), &
             clmamb(:,:), i1, n_lma &
             ) )
        enddo
!       print*,overlap(1,1,i_lmb),overlap_grad(1,1,i_lmb,2)
        i_lmb = i_lmb+1
     end do
  enddo

 ! calculate double differential rule on C(1,m) and multiply all scaling factors
  diff_rule_result = 0.0_r8_kind
  diff_rule_result_grad = 0.0_r8_kind
  diff_rule_result(:,1,1,1) = & !z
       pack( &
       spread(aexps*xa(3),1,nbexps) + &
       spread(bexps*xb(3),2,naexps), &
       cutoff &
       ) / fact0
   diff_rule_result_grad(:,1,1,1,2) = bexp_arr / fact0
!   diff_rule_result_grad(:,1,1,1,2) = pack(spread(bexps,2,naexps),cutoff)&
!        / fact0
!   print*,'xb xa',xb(3),xa(3)

  diff_rule_result(:,1,1,2) = & !x
       pack( &
       spread(aexps*xa(1),1,nbexps) + &
       spread(bexps*xb(1),2,naexps), &
       cutoff &
       ) / fact0
   diff_rule_result_grad(:,1,1,2,3)=bexp_arr / fact0
  diff_rule_result(:,1,1,3) = & !y
       pack( &
       spread(aexps*xa(2),1,nbexps) + &
       spread(bexps*xb(2),2,naexps), &
       cutoff &
       ) / fact0
   diff_rule_result_grad(:,1,1,3,4)= bexp_arr / fact0
  fact6 = aexp_arr / fact0
  if(n_lma.ge.4) then
  do i_lma = 2, 4
     diff_rule_result(:,i_lma,1,i_lma-1) = fact6
  enddo
  end if

  fact6 = bexp_arr / fact0
  if(n_lmb.ge.4) then
  do i_lmb = 2, 4
     diff_rule_result(:,1,i_lmb,i_lmb-1) = fact6
  enddo
 end if

!!$print*,'diff_rule_result'
!       do i_lmb =1, n_lmb
!       do i_lma = 1,n_lma
!       print*,diff_rule_result(1,i_lma,i_lmb,1), &
!               diff_rule_result_grad(1,i_lma,i_lmb,1,2)
!       enddo
!       enddo

  ! double product rule with respect to a and b
 do i_xyz = 1, 3
    do iota=2,4
       rnabla_primitive(:,:,:,lm2xyz(i_xyz+1),lm2xyz(iota)) = &
            prod_rule_double( &
            diff_rule_result_grad(:,:,:,i_xyz,iota), &
            overlap, &
            la**2+1,la**2+2*la+1, &
            lb**2+1,lb**2+2*lb+1 &
            )- &
            prod_rule_double( &
            diff_rule_result(:,:,:,i_xyz), &
            overlap_grad(:,:,:,iota), &
            la**2+1,la**2+2*la+1, &
            lb**2+1,lb**2+2*lb+1 &
            )
    enddo
 enddo

if(.false. .and. la.eq.0.and.lb.eq.2)then
DPRINT 'orbzeeman ll=2,0'
    do ma=1,2*la+1
     do mb=1,2*lb+1
        print*,ma,mb
!!$        print*,rnabla_primitive(1,ma,mb,1,1),rnabla_primitive(1,ma,mb,2,1),&
!!$             rnabla_primitive(1,ma,mb,3,1)
        print*,rnabla_primitive(1,ma,mb,1,2),rnabla_primitive(1,ma,mb,2,2),&
             rnabla_primitive(1,ma,mb,3,2)
        
     end do
  end do
  
endif

orbzeeman(:,:,:,:)=&
     reshape(vec_prod_vector(reshape(rnabla_primitive, & ! the - sign deleted 5.12.99
     (/num*(2*la+1)*(2*lb+1),3,3/)))&
     ,(/num,(2*la+1),(2*lb+1),3/))

corrzeeman(:,:,:,:)=0.0_r8_kind

   do ma=1,2*la+1
      do mb=1,2*lb+1
         do k = 1, 13
            select case(k)
            case(4) ! Sigma
               prim_int_2cob_dipoleg(:,:,mb,ma,4) = &
                    unpack(overlap(:,la**2+ma,lb**2+mb),cutoff,zero)
            case (1:3)
               prim_int_2cob_dipoleg(:,:,mb,ma,k) = &
                    unpack(orbzeeman(:,ma,mb,k),cutoff,zero)

               !DG gauge correction----------------------------------------------------

!!$            corrzeeman(:,:,:,:)=&
!!$                 reshape(calc_correction_integrals(reshape(rnabla_primitive,&
!!$                 (/num*(2*la+1)*(2*lb+1),3,3/)),k)&
!!$                 ,(/num,(2*la+1),(2*lb+1),3/))
!!$            prim_int_2cob_dipoleg(:,:,mb,ma,4+(k-1)*3+1) = & !DG 5,8,11
!!$                 unpack(corrzeeman(:,ma,mb,1),cutoff,zero)
!!$            prim_int_2cob_dipoleg(:,:,mb,ma,4+(k-1)*3+2) = & !6,9,12
!!$                 unpack(corrzeeman(:,ma,mb,2),cutoff,zero)
!!$            prim_int_2cob_dipoleg(:,:,mb,ma,4+(k-1)*3+3) = & !7,10,13
!!$                 unpack(corrzeeman(:,ma,mb,3),cutoff,zero)
!!$            case (5:13)
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k) = &
!!$                    unpack(overlap(:,la**2+ma,lb**2+mb),cutoff,zero)
!!$            case default
!!$
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k) = &
!!$                    unpack(corrzeeman(:,ma,mb,3),cutoff,zero)
!!$               case (5,8,11)
!!$               corrzeeman(:,:,:,:)=&
!!$                 reshape(calc_correction_integrals(reshape(rnabla_primitive,&
!!$                 (/num*(2*la+1)*(2*lb+1),3,3/)),k)&
!!$                 ,(/num,(2*la+1),(2*lb+1),3/))   
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k) = &
!!$                 unpack(corrzeeman(:,ma,mb,1),cutoff,zero)
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k+1) = & 
!!$                 unpack(corrzeeman(:,ma,mb,2),cutoff,zero)
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k+2) = & 
!!$                 unpack(corrzeeman(:,ma,mb,3),cutoff,zero)
!!$            case (5,9,13)
!!$               prim_int_2cob_dipoleg(:,:,mb,ma,k) = &
!!$                    unpack(overlap(:,la**2+ma,lb**2+mb),cutoff,zero)   
              case default

               prim_int_2cob_dipoleg(:,:,mb,ma,k) = 0.0_r8_kind
                  
            end select
         end do
         if(.false..and.la.eq.1.and.lb.eq.1) then
            print*,ma,mb
            print*,prim_int_2cob_dipoleg(1,1,mb,ma,1), &
                 prim_int_2cob_dipoleg(1,1,mb,ma,2),&
                 prim_int_2cob_dipoleg(1,1,mb,ma,3)
         end if

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
       overlap_grad, &
       rnabla_primitive, &
        orbzeeman, &
       corrzeeman,&
       diff_rule_result, &
       diff_rule_result_grad, &
       aexp_arr, &
       bexp_arr, &
       clmamb_scalar, &
       grad_clmamb_scalar, &
       clmamb, &
       clmamb_grad, &
       diff_arr0, &
       diff_arr0_grad, &
       cutoff, &
       stat=alloc_stat)

  contains

    function vec_prod_vector(r_nabla_elements) result (p)
      !
      ! vector product constartcts
      !
      !----------------------------------------------------------------
      ! Modifications
      !----------------------------------------------------------------
      !
      ! Modification (Please copy before editing)
      ! Author: ...DG
      ! Date:   ...11/03/2001
      ! Description: ...Changed sign in Lz component in finction
      ! 
      !
      !----------------------------------------------------------------

      USE type_module
      IMPLICIT NONE
     real(kind=r8_kind), intent(in), dimension(:,:,:):: r_nabla_elements
     real(kind=r8_kind),dimension(size(r_nabla_elements,1),3)::p

      p(:,1)=(r_nabla_elements(:,2,3)-r_nabla_elements(:,3,2))
      p(:,2)=(r_nabla_elements(:,3,1)-r_nabla_elements(:,1,3))
      !FIXME: sign changed (undo whenever you find and correct the
      ! reason this change is needed for):
      p(:,3)=-(r_nabla_elements(:,1,2)-r_nabla_elements(:,2,1))


    end function vec_prod_vector

    function calc_correction_integrals (r_nabla_elements,k) result (p)
      
      USE type_module
      IMPLICIT NONE
     real(kind=r8_kind), intent(in), dimension(:,:,:):: r_nabla_elements
     integer(kind=i4_kind), intent (in) :: k
     real(kind=r8_kind),dimension(size(r_nabla_elements,1),3)::p

    select case (k)
     case(5)
        p(:,1) = r_nabla_elements(:,2,2)+r_nabla_elements(:,3,3)
        p(:,2) = -r_nabla_elements(:,2,1)
        p(:,3) = -r_nabla_elements(:,3,1)
     case(8)
        p(:,1) = -r_nabla_elements(:,1,2) 
        p(:,2) = r_nabla_elements(:,1,1)+r_nabla_elements(:,3,3)
        p(:,3) = -r_nabla_elements(:,3,2)
     case(11)
        p(:,1) =  -r_nabla_elements(:,1,3) 
        p(:,2) =  -r_nabla_elements(:,2,3)
        p(:,3) = r_nabla_elements(:,1,1)+r_nabla_elements(:,2,2) 
     end select

    end function calc_correction_integrals
end subroutine ll_calculate_dipoleg


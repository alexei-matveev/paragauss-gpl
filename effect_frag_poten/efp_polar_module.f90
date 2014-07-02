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
module efp_polar_module
  !-------------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
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

# include <def.h>
  use type_module ! type specification parameters
  use datatype    ! user defined types
  use common_data_module
  use induced_dipoles_module
  use elec_static_field_module
  use operations_module, only: operations_solvation_effect
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  logical, public :: calc_polar=.false.

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public calc_E_mp, calc_E_id, allocate_Efield, deallocate_Efield, calc_id
  public calc_id_FFenergy, calc_id_FFgrads, total_id, calc_id1

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  subroutine allocate_Efield
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code ------------------------------------

    if(.not.allocated(E_mp)) then
       allocate(E_mp(N_ipd),E_id(N_ipd),E_id1(N_ipd), stat=status)
       ASSERT(status==0)
    end if

    do i=1,N_ipd
       if(.not.allocated(E_mp(i)%m)) then
          allocate(E_mp(i)%m(3,n_eq),E_id(i)%m(3,n_eq),E_id1(i)%m(3,n_eq), stat=status)
          ASSERT(status==0)
       end if

       E_mp(i)%m=zero
       E_id(i)%m=zero
       E_id1(i)%m=zero
    end do

  end subroutine allocate_Efield
  !*************************************************************

  !*************************************************************
  subroutine deallocate_Efield
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code ------------------------------------

    do i=1,N_ipd
       n_eq=ipd_array(i)%N_equal_dipoles
       deallocate(E_mp(i)%m,E_id(i)%m,E_id1(i)%m, stat=status)
       ASSERT(status==0)
    end do
    deallocate(E_mp,E_id,E_id1, stat=status)
    ASSERT(status==0)

  end subroutine deallocate_Efield
  !*************************************************************

  !*************************************************************
  subroutine calc_E_mp()
    !  Purpose: calculate electrostatic field produced by different
    !           multipole centers at positions of induced dipoles
    !------------ Modules used -----------------------------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: po_array, N_po
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iu_id,ie_id,n_eq,iu_mp,ie_mp,i,j,k
    real(r8_kind) :: E(3),r_id(3),r_mp(3),rr(3),dr,E1(3)
    real(r8_kind) :: Z,D(3),Q(3,3)
!!$    real(r8_kind) ::V,V1
    integer(i4_kind) :: e_dim,index,id_grp,mp_grp
    real(r8_kind), pointer :: rotmat(:,:)
    !------------ Executable code ------------------------------------

    ASSERT(N_ipd==N_surface_points)

    totalsym_field=zero

    I_PD: do iu_id=1,N_surface_points

!!$       V=zero
       E=zero
       n_eq=surface_points(iu_id)%N_equal_points
       e_dim=surf_points_grad_index(iu_id+1)-surf_points_grad_index(iu_id)
       index=surf_points_grad_index(iu_id)-1

       do ie_id=1,n_eq

          id_grp=ipd_array(iu_id)%group(ie_id)
          r_id=surface_points(iu_id)%position(:,ie_id)
          rotmat=>surf_points_grad_info(iu_id)%m(:,:,ie_id)

          !POINT CHARGES
          PC: do iu_mp=1,pointcharge_N

             Z=pointcharge_array(iu_mp)%Z
          
             do ie_mp=1,pointcharge_array(iu_mp)%N_equal_charges

                mp_grp=pointcharge_array(iu_mp)%group(ie_mp)
                if(id_grp==mp_grp) cycle
                r_mp=pointcharge_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_id
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!!$                V=V + Z/dr
                E = E - Z*rr/dr**3
             enddo
          enddo PC

          !POINT DIPOLES
          PD: do iu_mp=1,N_pd
             do ie_mp=1,pd_array(iu_mp)%N_equal_dipoles

                mp_grp=pd_array(iu_mp)%group(ie_mp)
                if(id_grp==mp_grp) cycle
                D=pd_array(iu_mp)%dipole(:,ie_mp)
                r_mp=pd_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_id
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!!$                V=V + dot_product(D,rr)/dr**3
                E=E + three*dot_product(D,rr)*rr/dr**5-D/dr**3
             end do
          end do PD

          !POINT QUADRUPOLES
          PQ: do iu_mp=1,N_pq
             do ie_mp=1,pq_array(iu_mp)%N_equal_qpoles

                mp_grp=pq_array(iu_mp)%group(ie_mp)
                if(id_grp==mp_grp) cycle
                Q=pq_array(iu_mp)%quadrupole(:,:,ie_mp)
                r_mp=pq_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_id
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!!$                V1=zero
                E1=zero
                do i=1,3
                   do j=1,3
!!$                      V1=V1+three*rr(i)*rr(j)*Q(i,j)/dr**5
!!$                      if(i==j) V1=V1-Q(i,j)/dr**3

                      E1=E1-(fifteen*rr(i)*rr(j)*rr/dr**7)*Q(i,j)
                      do k=1,3
                         if(k==i) E1(k)=E1(k)+(three*rr(j)/dr**5)*Q(i,j)
                         if(k==j) E1(k)=E1(k)+(three*rr(i)/dr**5)*Q(i,j)
                      end do
!                      if(i==j) E1=E1+(three*rr/dr**5)*Q(i,j)
                   end do
                end do
!!$                V=V+V1/three
                E=E+E1/three
             end do
          end do PQ

#if 0
          !POINT OCTOPOLES
          PO: do iu_mp=1,N_po
             do ie_mp=1,po_array(iu_mp)%N_equal_opoles

                mp_grp=po_array(iu_mp)%group(ie_mp)
                if(id_grp==mp_grp) cycle
                O=po_array(iu_mp)%octopole(:,:,:,ie_mp)
                r_mp=po_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_id
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!!$                V1=zero
                E1=zero
                do i=1,3
                   do j=1,3
                      do k=1,3
!!$                         V1=V1+fifteen*rr(i)*rr(j)*rr(k)*O(i,j,k)/dr**7
!!$                         if(i==j) V1=V1-three*rr(k)*O(i,j,k)/dr**5
!!$                         if(i==k) V1=V1-three*rr(j)*O(i,j,k)/dr**5
!!$                         if(j==k) V1=V1-three*rr(i)*O(i,j,k)/dr**5

                         E1=E1-(105.0_r8_kind*rr(i)*rr(j)*rr(k)/dr**9)*rr*O(i,j,k)
                         do l=1,3
                            if(l==i) E1(l)=E1(l)+(fifteen*rr(j)*rr(k)/dr**7)*O(i,j,k)
                            if(l==j) E1(l)=E1(l)+(fifteen*rr(i)*rr(k)/dr**7)*O(i,j,k)
                            if(l==k) E1(l)=E1(l)+(fifteen*rr(i)*rr(j)/dr**7)*O(i,j,k)
                         end do
                         if(i==j) E1=E1+(fifteen*rr(k)*rr/dr**7)*O(i,j,k)
                         if(i==k) E1=E1+(fifteen*rr(j)*rr/dr**7)*O(i,j,k)
                         if(j==k) E1=E1+(fifteen*rr(i)*rr/dr**7)*O(i,j,k)
                         do l=1,3
                            if(i==j .and. k==l) E1(l)=E1(l)-(three*O(i,j,k)/dr**5)
                            if(i==k .and. j==l) E1(l)=E1(l)-(three*O(i,j,k)/dr**5)
                            if(j==k .and. i==l) E1(l)=E1(l)-(three*O(i,j,k)/dr**5)
                         end do
                      end do
                   end do
                end do
!!$                V=V + V1/fifteen
                E=E - E1/fifteen
             end do
          end do PO
#endif

          do i=1,e_dim
             totalsym_field(index+i)=totalsym_field(index+i)- & !+
                  sum(rotmat(i,:)*E(:))
          enddo

       end do
    end do I_PD

    call transform_to_cart_field(E_mp,totalsym_field)
!!$do i=1,N_surface_points
!!$print*,'E_mp',E_mp(i)%m(:,1)
!!$end do

    totalsym_field=0.0_r8_kind

  end subroutine calc_E_mp
  !*************************************************************

  !*************************************************************
  subroutine calc_E_id()
    !  Purpose: calculate electrostatic field produced by 
    !           induced dipoles at positions of other induced dipoles
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iu_id1,ie_id1,n_eq1,iu_id2,ie_id2,n_eq2,i
    real(r8_kind) :: E(3),r_id1(3),r_id2(3),r12(3),dr
    real(r8_kind) :: D2(3)
!    real(r8_kind) ::V
    integer(i4_kind) :: e_dim,index,id_grp1,id_grp2
    real(r8_kind), pointer :: rotmat(:,:)
    !------------ Executable code ------------------------------------

    ASSERT(N_ipd==N_surface_points)

    totalsym_field=0.0_r8_kind

    I_PD1: do iu_id1=1,N_ipd

!       V=zero
       E=zero
       n_eq1=surface_points(iu_id1)%N_equal_points
       e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
       index=surf_points_grad_index(iu_id1)-1

       do ie_id1=1,n_eq1

          id_grp1=ipd_array(iu_id1)%group(ie_id1)
          r_id1=surface_points(iu_id1)%position(:,ie_id1)
          rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

          I_PD2: do iu_id2=1,N_ipd
             n_eq2=ipd_array(iu_id2)%N_equal_dipoles

             do ie_id2=1,n_eq2

                id_grp2=ipd_array(iu_id2)%group(ie_id2)
                if(id_grp1==id_grp2) cycle
                D2=ipd_array(iu_id2)%idipole(:,ie_id2)
                r_id2=ipd_array(iu_id2)%position(:,ie_id2)
                r12=r_id2-r_id1
                dr=sqrt(dot_product(r12,r12))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!                V=V + dot_product(D2,r12)/dr**3
                E=E + three*dot_product(D2,r12)*r12/dr**5-D2/dr**3
             end do
          end do I_PD2
          do i=1,e_dim
             totalsym_field(index+i)=totalsym_field(index+i)- & !+
                  sum(rotmat(i,:)*E(:))
          end do
       end do
    end do I_PD1

    call transform_to_cart_field(E_id,totalsym_field)

    totalsym_field=0.0_r8_kind

    I_PD3: do iu_id1=1,N_ipd

!       V=zero
       E=zero
       n_eq1=surface_points(iu_id1)%N_equal_points
       e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
       index=surf_points_grad_index(iu_id1)-1

       do ie_id1=1,n_eq1

          id_grp1=ipd_array(iu_id1)%group(ie_id1)
          r_id1=surface_points(iu_id1)%position(:,ie_id1)
          rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

          I_PD4: do iu_id2=1,N_ipd
             n_eq2=ipd_array(iu_id2)%N_equal_dipoles

             do ie_id2=1,n_eq2

                id_grp2=ipd_array(iu_id2)%group(ie_id2)
                if(id_grp1==id_grp2) cycle
                D2=ipd_array(iu_id2)%idipole1(:,ie_id2)
                r_id2=ipd_array(iu_id2)%position(:,ie_id2)
                r12=r_id2-r_id1
                dr=sqrt(dot_product(r12,r12))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
!                V=V + dot_product(D2,r12)/dr**3
                E=E + three*dot_product(D2,r12)*r12/dr**5-D2/dr**3
             end do
          end do I_PD4
          do i=1,e_dim
             totalsym_field(index+i)=totalsym_field(index+i)- & !+
                  sum(rotmat(i,:)*E(:))
          end do
       end do
    end do I_PD3

    call transform_to_cart_field(E_id1,totalsym_field)

    totalsym_field=0.0_r8_kind

  end subroutine calc_E_id
  !*************************************************************

  !*************************************************************
  subroutine calc_E_id1(iu_id1,ie_id1)
    !  Purpose: calculate electrostatic field produced by 
    !           induced dipoles at positions of other induced dipoles
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind) :: iu_id1,ie_id1
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n_eq1,iu_id2,ie_id2,n_eq2
    real(r8_kind) :: E(3),r_id1(3),r_id2(3),r12(3),dr
    real(r8_kind) :: D2(3)
!    real(r8_kind) ::V
    integer(i4_kind) :: e_dim,index,id_grp1,id_grp2
    real(r8_kind), pointer :: rotmat(:,:)
    !------------ Executable code ------------------------------------

    ASSERT(N_ipd==N_surface_points)

    totalsym_field=0.0_r8_kind

!       V=zero
    E=zero
    n_eq1=surface_points(iu_id1)%N_equal_points
    e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
    index=surf_points_grad_index(iu_id1)-1

    id_grp1=ipd_array(iu_id1)%group(ie_id1)
    r_id1=surface_points(iu_id1)%position(:,ie_id1)
    rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

    I_PD2: do iu_id2=1,N_ipd
       n_eq2=ipd_array(iu_id2)%N_equal_dipoles

       do ie_id2=1,n_eq2

          id_grp2=ipd_array(iu_id2)%group(ie_id2)
          if(id_grp1==id_grp2) cycle
          D2=ipd_array(iu_id2)%idipole(:,ie_id2)
          r_id2=ipd_array(iu_id2)%position(:,ie_id2)
          r12=r_id2-r_id1
          dr=sqrt(dot_product(r12,r12))
          if(dr < 0.001_r8_kind) dr=0.001_r8_kind
          !V=V + dot_product(D2,r12)/dr**3
          E=E + three*dot_product(D2,r12)*r12/dr**5-D2/dr**3
       end do
    end do I_PD2
    E_id(iu_id1)%m(:,ie_id1)=E

  end subroutine calc_E_id1
  !*************************************************************

  !*************************************************************
  subroutine calc_E_id11(iu_id1,ie_id1)
    !  Purpose: calculate electrostatic field produced by 
    !           induced dipoles at positions of other induced dipoles
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind) :: iu_id1,ie_id1
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n_eq1,iu_id2,ie_id2,n_eq2
    real(r8_kind) :: E(3),r_id1(3),r_id2(3),r12(3),dr
    real(r8_kind) :: D2(3)
!    real(r8_kind) ::V
    integer(i4_kind) :: e_dim,index,id_grp1,id_grp2
    real(r8_kind), pointer :: rotmat(:,:)
    !------------ Executable code ------------------------------------

    ASSERT(N_ipd==N_surface_points)

    totalsym_field=0.0_r8_kind

!       V=zero
    E=zero
    n_eq1=surface_points(iu_id1)%N_equal_points
    e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
    index=surf_points_grad_index(iu_id1)-1

    id_grp1=ipd_array(iu_id1)%group(ie_id1)
    r_id1=surface_points(iu_id1)%position(:,ie_id1)
    rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

    I_PD2: do iu_id2=1,N_ipd
       n_eq2=ipd_array(iu_id2)%N_equal_dipoles

       do ie_id2=1,n_eq2

          id_grp2=ipd_array(iu_id2)%group(ie_id2)
          if(id_grp1==id_grp2) cycle
          D2=ipd_array(iu_id2)%idipole1(:,ie_id2)
          r_id2=ipd_array(iu_id2)%position(:,ie_id2)
          r12=r_id2-r_id1
          dr=sqrt(dot_product(r12,r12))
          if(dr < 0.001_r8_kind) dr=0.001_r8_kind
          !V=V + dot_product(D2,r12)/dr**3
          E=E + three*dot_product(D2,r12)*r12/dr**5-D2/dr**3
       end do
    end do I_PD2
    E_id1(iu_id1)%m(:,ie_id1)=E

  end subroutine calc_E_id11
  !*************************************************************

  !*************************************************************
  subroutine calc_id(cycle)
    ! Purpose: Calculate induced dipole moments
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: cycle
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i,j,n_eq
    real(r8_kind) :: alpha(3,3),E(3)
    !--- executable code-------------------------------------

    E=0.0_r8_kind

    do i=1,N_ipd
       n_eq=ipd_array(i)%N_equal_dipoles
       do j=1,n_eq
          alpha=ipd_array(i)%pol_tensor(:,:,j)
          E=E_mp(i)%m(:,j)
          if(cycle > 1) E=E+E_id(i)%m(:,j)
          ipd_array(i)%idipole(:,j)=matmul(alpha,E)
          E=E_mp(i)%m(:,j)
          if(cycle > 1) E=E+E_id1(i)%m(:,j)
          ipd_array(i)%idipole1(:,j)=matmul(transpose(alpha),E)
       end do
    end do

  end subroutine calc_id
  !*************************************************************

  !*************************************************************
  subroutine calc_id1()
    ! Purpose: Calculate induced dipole moments
    !------------ Modules used -----------------------------------
    use solv_electrostat_module, only: A_matrix_inv
    use solv_cavity_module, only: dielectric_constant
    use efp_solv_module, only: calc_V_and_Q_id, calc_E_cav1, calc_E_cav11
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i,j,n_eq
    real(r8_kind) :: alpha(3,3),E(3),E1(3)
    !--- executable code-------------------------------------

    E=0.0_r8_kind

    do i=1,N_ipd
       n_eq=ipd_array(i)%N_equal_dipoles
       do j=1,n_eq
          alpha=ipd_array(i)%pol_tensor(:,:,j)
          E1=E_mp(i)%m(:,j)
          if(calc_polar) E1=E1+E_ele(i)%m(:,j)+E_nuc(i)%m(:,j)
          call calc_E_id1(i,j)
          E=E1+E_id(i)%m(:,j)
          if(operations_solvation_effect .and. do_pol_pcm) then 
             call calc_V_and_Q_id()
             call calc_E_cav1(i,j)
             E=E+E_cav(i)%m(:,j)
          end if
          ipd_array(i)%idipole(:,j)=matmul(alpha,E)

          call calc_E_id11(i,j)
          E=E1+E_id1(i)%m(:,j)
          if(operations_solvation_effect .and. do_pol_pcm) then 
             call calc_E_cav11(i,j)
             E=E+E_cav(i)%m(:,j)
          end if
          ipd_array(i)%idipole1(:,j)=matmul(transpose(alpha),E)
       end do
    end do

  end subroutine calc_id1
  !*************************************************************

  !*************************************************************
 subroutine calc_id_FFenergy(en_ind_dipole)
    ! Purpose: Calculate polarizable efp-efp energy
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(out) :: en_ind_dipole
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_uc,i_ec
    integer(i4_kind) :: n_ec
    real(r8_kind),pointer  :: idipole(:,:)
    real(r8_kind) :: E_f(3)
    !--- executable code-------------------------------------

    en_ind_dipole=0.0_r8_kind

    if(N_ipd > 0) then
       do i_uc=1,N_ipd
          n_ec=ipd_array(i_uc)%N_equal_dipoles
          idipole=>ipd_array(i_uc)%idipole

          do i_ec=1,n_ec
             E_f=E_mp(i_uc)%m(:,i_ec)
             if(calc_polar) E_f=E_f+E_ele(i_uc)%m(:,i_ec)+E_nuc(i_uc)%m(:,i_ec)
             en_ind_dipole=en_ind_dipole- &
                  dot_product(E_f,idipole(:,i_ec))*half
          end do
       end do
!print*,en_ind_dipole,'en_ind_dipole'
    end if

  end subroutine calc_id_FFenergy
  !*************************************************************

  !*************************************************************
  subroutine calc_id_FFgrads
    !  Purpose: calculate gradients of polarization energy
    !  The carrent approach considers polarizability tensors
    !  as symmetric to simplify expression for gradient
    !  calculations. Therefore, the used approach is not
    !  exact but approximate
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use pointcharge_module, only: unique_pc_grad_info,unique_pc_index
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: dipoles_grad_index,dipoles_grad_info
    use point_dqo_module, only: gradient_dip_totalsym,torque_dip_totalsym
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: qpoles_grad_index,qpoles_grad_info
    use point_dqo_module, only: gradient_quad_totalsym,torque_quad_totalsym
#if 0
    use point_dqo_module, only: po_array, N_po
    use point_dqo_module, only: opoles_grad_index,opoles_grad_info
    use point_dqo_module, only: gradient_oct_totalsym
#endif
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: iu_id,ie_id,ne_id,iu_mp,ie_mp,ne_mp,i,j,k,i1,j1
    integer(i4_kind) :: iu_id2,ie_id2,ne_id2,ind2
    integer(i4_kind) :: grad_dim1,grad_dim2,index1,index2,i_grp1,i_grp2
    real(r8_kind) :: r1(3),r2(3),r12(3),dr,G(3)
    real(r8_kind) :: ID(3),ID1(3),ID2(3),Z2,D2(3),Q2(3,3),gradient1(3),gradient2(3)
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    real(r8_kind) :: torque1(3), torque2(3),E(3),E1(3),E2(3),Qr(3),EE(3)
    !------------ Executable code ------------------------------------

    do iu_id=1,N_ipd
       ne_id=ipd_array(iu_id)%N_equal_dipoles
       grad_dim1=ind_dip_grad_index(iu_id+1)-ind_dip_grad_index(iu_id)
       index1=ind_dip_grad_index(iu_id)-1
        
       do ie_id=1,ne_id
          rotmat1=>ind_dip_grad_info(iu_id)%m(:,:,ie_id)
          gradient1=zero
          torque1=zero
          i_grp1=ipd_array(iu_id)%group(ie_id)
          r1=ipd_array(iu_id)%position(:,ie_id)
          ID=ipd_array(iu_id)%idipole(:,ie_id)+ipd_array(iu_id)%idipole1(:,ie_id)
          ID1=ipd_array(iu_id)%idipole(:,ie_id)

          !ID - PC
          id_pc: do iu_mp=1,pointcharge_N
             grad_dim2=unique_pc_index(iu_mp+1)-unique_pc_index(iu_mp)
             index2=unique_pc_index(iu_mp)-1
             ne_mp=pointcharge_array(iu_mp)%N_equal_charges
             
             do ie_mp=1,ne_mp
                i_grp2=pointcharge_array(iu_mp)%group(ie_mp)
                if(i_grp2 == i_grp1) cycle

                rotmat2=>unique_pc_grad_info(iu_mp)%m(:,:,ie_mp)
                gradient2=zero
                r2=pointcharge_array(iu_mp)%position(:,ie_mp)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind
                Z2=pointcharge_array(iu_mp)%Z

                G=-three*Z2*dot_product(ID,r12)*r12/dr**5+Z2*ID/dr**3
                gradient1=gradient1-G
                gradient2=gradient2+G

                E=Z2*r12/dr**3
                torque1=torque1-vector_product(ID,E)

                do j1=1,grad_dim2
                   gradient_pc_totalsym(index2+j1)=gradient_pc_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*gradient2(:))*half
                enddo
             end do
          end do id_pc

          !ID - PD
          id_pd: do iu_mp=1,N_pd
             grad_dim2=dipoles_grad_index(iu_mp+1)-dipoles_grad_index(iu_mp)
             index2=dipoles_grad_index(iu_mp)-1
             ne_mp=pd_array(iu_mp)%N_equal_dipoles

             do ie_mp=1,ne_mp
                i_grp2=pd_array(iu_mp)%group(ie_mp)
                if(i_grp2 == i_grp1) cycle

                rotmat2=>dipoles_grad_info(iu_mp)%m(:,:,ie_mp)
                gradient2=zero
                torque2=zero
                r2=pd_array(iu_mp)%position(:,ie_mp)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind
                D2=pd_array(iu_mp)%dipole(:,ie_mp)

                G=-three*dot_product(ID,D2)*r12/dr**5- &
                     three*(ID*dot_product(D2,r12)+D2*dot_product(ID,r12))/dr**5+ &
                     fifteen*dot_product(ID,r12)*dot_product(D2,r12)*r12/dr**7
                gradient1=gradient1+G
                gradient2=gradient2-G

                E2=D2/dr**3-three*r12*dot_product(D2,r12)/dr**5
                E1=ID/dr**3-three*r12*dot_product(ID,r12)/dr**5
                torque1=torque1+vector_product(ID,E2)
                torque2=torque2+vector_product(D2,E1)

                do j1=1,grad_dim2
                   gradient_dip_totalsym(index2+j1)=gradient_dip_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*gradient2(:))*half
                   torque_dip_totalsym(index2+j1)=torque_dip_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*torque2(:))*half
                enddo
             end do
          end do id_pd

          !ID - PQ
          id_pq: do iu_mp=1,N_pq
             grad_dim2=qpoles_grad_index(iu_mp+1)-qpoles_grad_index(iu_mp)
             index2=qpoles_grad_index(iu_mp)-1
             ne_mp=pq_array(iu_mp)%N_equal_qpoles
 
             do ie_mp=1,ne_mp
                i_grp2=pq_array(iu_mp)%group(ie_mp)
                if(i_grp2 == i_grp1) cycle

                rotmat2=>qpoles_grad_info(iu_mp)%m(:,:,ie_mp)
                gradient2=zero
                torque2=zero
                r2=pq_array(iu_mp)%position(:,ie_mp)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind
                Q2=pq_array(iu_mp)%quadrupole(:,:,ie_mp)

                G=zero
                E2=zero
                EE=zero
                do i=1,3
                   do j=1,3
                      G=G-(fifteen*(ID(i)*r12(j)+ID(j)*r12(i))*r12/dr**7)*Q2(i,j)
                      do k=1,3
                         if(k==i) G(k)=G(k)+(three*ID(j)/dr**5)*Q2(i,j)
                         if(k==j) G(k)=G(k)+(three*ID(i)/dr**5)*Q2(i,j)
                      end do
                      G=G+(105.0_r8_kind*r12(i)*r12(j)*dot_product(r12,ID)*r12/dr**9)*Q2(i,j)
                      do k=1,3
                         if(k==i) G(k)=G(k)-(fifteen*r12(j)*dot_product(r12,ID)/dr**7)*Q2(i,j)
                         if(k==j) G(k)=G(k)-(fifteen*r12(i)*dot_product(r12,ID)/dr**7)*Q2(i,j)
                      end do
                      G=G-(fifteen*r12(i)*r12(j)*ID/dr**7)*Q2(i,j)

                      E2=E2-fifteen*r12(i)*r12(j)*r12*Q2(i,j)/dr**7
                      do k=1,3
                         if(i==k) E2(k)=E2(k)+three*r12(j)*Q2(i,j)/dr**5
                         if(j==k) E2(k)=E2(k)+three*r12(i)*Q2(i,j)/dr**5
                      end do
                   end do
                   EE(1)=EE(1)+(Q2(2,i)*ID(i)*r12(3)-Q2(3,i)*ID(i)*r12(2))/dr**5
                   EE(2)=EE(2)+(Q2(3,i)*ID(i)*r12(1)-Q2(1,i)*ID(i)*r12(3))/dr**5
                   EE(3)=EE(3)+(Q2(1,i)*ID(i)*r12(2)-Q2(2,i)*ID(i)*r12(1))/dr**5
                end do

                gradient1=gradient1+G/three
                gradient2=gradient2-G/three

                E2=E2/three
                torque1=torque1+vector_product(ID,E2)
                Qr=matmul(Q2,r12)
                torque2=torque2+two*EE
                EE=two*(ID/dr**5-five*dot_product(ID,r12)*r12/dr**7)
                torque2=torque2+vector_product(Qr,EE)

                do j1=1,grad_dim2
                   gradient_quad_totalsym(index2+j1)=gradient_quad_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*gradient2(:))*half
                   torque_quad_totalsym(index2+j1)=torque_quad_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*torque2(:))*half
                enddo
             end do
          end do id_pq

#if 0
          !ID - PO
          id_po: do iu_mp=1,N_po
             grad_dim2=opoles_grad_index(iu_mp+1)-opoles_grad_index(iu_mp)
             index2=opoles_grad_index(iu_mp)-1
             ne_mp=po_array(iu_mp)%N_equal_opoles
 
             do ie_mp=1,ne_mp
                i_grp2=po_array(iu_mp)%group(ie_mp)
                if(i_grp2 == i_grp1) cycle

                rotmat2=>opoles_grad_info(iu_mp)%m(:,:,ie_mp)
                gradient2=zero
                r2=po_array(iu_mp)%position(:,ie_mp)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind
                O2=po_array(iu_mp)%octopole(:,:,:,ie_mp)

                G=zero
                do i=1,3
                   do j=1,3
                      do k=1,3
                         G=G+945.0_r8_kind*r12(i)*r12(j)*r12(k)*dot_product(ID,r12)*r12*O2(i,j,k)/dr**11
                         do l=1,3
                            if(l==i) G(l)=G(l)- &
                                 105.0_r8_kind*r12(j)*r12(k)*dot_product(ID,r12)*O2(i,j,k)/dr**9- &
                                 105.0_r8_kind*r12(i)*r12(j)*r12(k)*ID(i)*O2(i,j,k)/dr**9
                            if(l==j) G(l)=G(l)- &
                                 105.0_r8_kind*r12(i)*r12(k)*dot_product(ID,r12)*O2(i,j,k)/dr**9- &
                                 105.0_r8_kind*r12(i)*r12(j)*r12(k)*ID(j)*O2(i,j,k)/dr**9
                            if(l==k) G(l)=G(l)- &
                                 105.0_r8_kind*r12(i)*r12(j)*dot_product(ID,r12)*O2(i,j,k)/dr**9- &
                                 105.0_r8_kind*r12(i)*r12(j)*r12(k)*ID(k)*O2(i,j,k)/dr**9
                         end do
                         G=G-105.0_r8_kind*(ID(i)*r12(j)*r12(k)+ID(j)*r12(i)*r12(k)+ID(k)*r12(i)*r12(j))*r12*O2(i,j,k)/dr**9
                         do l=1,3
                            if(l==i) G(l)=G(l)+fifteen*(ID(j)*r12(k)+ID(k)*r12(j))*O2(i,j,k)/dr**7
                            if(l==j) G(l)=G(l)+fifteen*(ID(i)*r12(k)+ID(k)*r12(i))*O2(i,j,k)/dr**7
                            if(l==k) G(l)=G(l)+fifteen*(ID(i)*r12(j)+ID(j)*r12(i))*O2(i,j,k)/dr**7
                         end do
                         if(i==j) G=G-105.0_r8_kind*r12(k)*dot_product(ID,r12)*r12*O2(i,j,k)/dr**9
                         if(i==k) G=G-105.0_r8_kind*r12(j)*dot_product(ID,r12)*r12*O2(i,j,k)/dr**9
                         if(j==k) G=G-105.0_r8_kind*r12(i)*dot_product(ID,r12)*r12*O2(i,j,k)/dr**9
                         do l=1,3
                            if(i==j .and. k==l) G(l)=G(l)+fifteen*dot_product(ID,r12)*O2(i,j,k)/dr**7
                            if(i==k .and. j==l) G(l)=G(l)+fifteen*dot_product(ID,r12)*O2(i,j,k)/dr**7
                            if(j==k .and. i==l) G(l)=G(l)+fifteen*dot_product(ID,r12)*O2(i,j,k)/dr**7
                         end do
                         do l=1,3
                            if(i==j) then
                               if(i==l) G(l)=G(l)+fifteen*r12(k)*ID(i)*O2(i,j,k)/dr**7
                               if(j==l) G(l)=G(l)+fifteen*r12(k)*ID(j)*O2(i,j,k)/dr**7
                               if(k==l) G(l)=G(l)+fifteen*r12(k)*ID(k)*O2(i,j,k)/dr**7
                            end if
                            if(i==k) then
                               if(i==l) G(l)=G(l)+fifteen*r12(j)*ID(i)*O2(i,j,k)/dr**7
                               if(j==l) G(l)=G(l)+fifteen*r12(j)*ID(j)*O2(i,j,k)/dr**7
                               if(k==l) G(l)=G(l)+fifteen*r12(j)*ID(k)*O2(i,j,k)/dr**7
                            end if
                            if(j==k) then
                               if(i==l) G(l)=G(l)+fifteen*r12(i)*ID(i)*O2(i,j,k)/dr**7
                               if(j==l) G(l)=G(l)+fifteen*r12(i)*ID(j)*O2(i,j,k)/dr**7
                               if(k==l) G(l)=G(l)+fifteen*r12(i)*ID(k)*O2(i,j,k)/dr**7
                            end if
                         end do
                         if(i==j) G=G+fifteen*ID(k)*O2(i,j,k)/dr**7
                         if(i==k) G=G+fifteen*ID(j)*O2(i,j,k)/dr**7
                         if(j==k) G=G+fifteen*ID(i)*O2(i,j,k)/dr**7
                     end do
                   end do
                end do

                gradient1=gradient1+G/fifteen
                gradient2=gradient2-G/fifteen
                do j1=1,grad_dim2
                   gradient_oct_totalsym(index2+j1)=gradient_oct_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*gradient2(:))*half
                enddo
             end do
          end do id_po
#endif

          !IPD - IPD
          ind2=0
          id_id1: do iu_id2=1,N_ipd
             ne_id2=ipd_array(iu_id2)%N_equal_dipoles
             grad_dim2=ind_dip_grad_index(iu_id2+1)-ind_dip_grad_index(iu_id2)
             index2=ind_dip_grad_index(iu_id2)-1
        
             do ie_id2=1,ne_id2
                i_grp2=ipd_array(iu_id2)%group(ie_id2)
                if(i_grp2 == i_grp1) cycle

                rotmat2=>ind_dip_grad_info(iu_id2)%m(:,:,ie_id2)
                gradient2=zero
                torque2=zero
                r2=ipd_array(iu_id2)%position(:,ie_id2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                ID2=ipd_array(iu_id2)%idipole1(:,ie_id2)

                G=-three*dot_product(ID1,ID2)*r12/dr**5- &
                     three*(ID1*dot_product(ID2,r12)+ID2*dot_product(ID1,r12))/dr**5+ &
                     fifteen*dot_product(ID1,r12)*dot_product(ID2,r12)*r12/dr**7
                gradient1=gradient1+G
                gradient2=gradient2-G

                E2=ID2/dr**3-three*r12*dot_product(ID2,r12)/dr**5
                E1=ID1/dr**3-three*r12*dot_product(ID1,r12)/dr**5
                torque1=torque1+vector_product(ID1,E2)
                torque2=torque2+vector_product(ID2,E1)

                do j1=1,grad_dim2
                   grad_idip_totalsym(index2+j1)=grad_idip_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*gradient2(:))*half
                   torque_idip_totalsym(index2+j1)=torque_idip_totalsym(index2+j1)+ &
                        sum(rotmat2(j1,:)*torque2(:))*half
                enddo
             end do
          end do id_id1
          
          !Gradients on polarizable points
          do i1=1,grad_dim1
             grad_idip_totalsym(index1+i1)=grad_idip_totalsym(index1+i1)+ &
                  sum(rotmat1(i1,:)*gradient1(:))*half
             torque_idip_totalsym(index1+i1)=torque_idip_totalsym(index1+i1)+ &
                  sum(rotmat1(i1,:)*torque1(:))*half
          enddo

       end do
    end do

  end subroutine calc_id_FFgrads
  !*************************************************************

  !*************************************************************
  function vector_product(v1,v2)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: vector_product(3)
    real(r8_kind) :: v1(3),v2(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------
    
    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !*************************************************************

  !*************************************************************
  function total_id()
    !  Purpose: Calculate total induced dipole moment of EFP system
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: total_id(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i,j,ne
    !------------ Executable code ------------------------------------

    total_id=zero

    do i=1,N_ipd
       ne=ipd_array(i)%N_equal_dipoles
       do j=1,ne
          total_id=total_id+ipd_array(i)%idipole(:,j)
       end do
    end do

  end function total_id
  !*************************************************************


  !--------------- End of module -------------------------------------
end module efp_polar_module

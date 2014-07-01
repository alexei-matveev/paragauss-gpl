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
module efp_solv_module
  !---------------------------------------------------------------
  !
  !  Purpose:
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: AS
  !  Date: 07.08
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
# include "def.h"

  use type_module ! type specification parameters
  use common_data_module
  use induced_dipoles_module
  use elec_static_field_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  logical, public :: do_solv

  !------------ Interface statements ------------------------------
  !------------ public functions and subroutines ------------------
  public calc_mp_potential, deallocate_V_mp, efp_mp_solv_energy
  public allocate_V_and_Q_id, allocate_E_cav, deallocate_E_cav, calc_E_cav1
  public efp_id_mp_solv_energy, deallocate_V_id, calc_V_and_Q_id, calc_E_cav11
  public calc_E_cav

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------
  !------------ Declaration of constants and variables ------------


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains
  !*************************************************************
  subroutine calc_mp_potential()
    !  Purpose: Calculate electrostatic potential produced by
    !           multipole centers at solute cavity points
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_mp, N_points, point_in_space
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: po_array, N_po1
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status, iu_sp,ie_sp,n_eq,iu_mp,ie_mp,i,j,k
    real(r8_kind) :: V,V1,r_sp(3),r_mp(3),rr(3),dr,dr2
    real(r8_kind) :: Z,D(3),Q(3,3),O(3,3,3),C,A
    !------------ Executable code --------------------------------

    allocate(V_pot_mp(N_points),stat=status)
    ASSERT(status==0)
    V_pot_mp=zero

    I_SP: do iu_sp=1,N_points

       V=zero
       n_eq=point_in_space(iu_sp)%N_equal_points

       do ie_sp=1,n_eq

          r_sp=point_in_space(iu_sp)%position(:,ie_sp)

          !POINT CHARGES
          PC: do iu_mp=1,pointcharge_N

             Z=pointcharge_array(iu_mp)%Z
             C=pointcharge_array(iu_mp)%C
             A=pointcharge_array(iu_mp)%A

             do ie_mp=1,pointcharge_array(iu_mp)%N_equal_charges

                r_mp=pointcharge_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_sp
                dr2=dot_product(rr,rr)
                dr=sqrt(dr2)
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                V=V + Z/dr
                if(C /= 0.0_r8_kind) then
                   V=V-C*exp(-A*dr2)*Z/dr
                end if
             enddo
          enddo PC

          !POINT DIPOLES
          PD: do iu_mp=1,N_pd
             do ie_mp=1,pd_array(iu_mp)%N_equal_dipoles

                D=pd_array(iu_mp)%dipole(:,ie_mp)
                r_mp=pd_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_sp
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                V=V - dot_product(D,rr)/dr**3
             end do
          end do PD

          !POINT QUADRUPOLES
          PQ: do iu_mp=1,N_pq
             do ie_mp=1,pq_array(iu_mp)%N_equal_qpoles

                Q=pq_array(iu_mp)%quadrupole(:,:,ie_mp)
                r_mp=pq_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_sp
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                V1=zero
                do i=1,3
                   do j=1,3
                      V1=V1+three*rr(i)*rr(j)*Q(i,j)/dr**5
                      if(i==j) V1=V1-Q(i,j)/dr**3
                   end do
                end do
                V=V+V1/three
             end do
          end do PQ

          !POINT OCTOPOLES
          PO: do iu_mp=1,N_po1
             do ie_mp=1,po_array(iu_mp)%N_equal_opoles

                O=po_array(iu_mp)%octopole(:,:,:,ie_mp)
                r_mp=po_array(iu_mp)%position(:,ie_mp)
                rr=r_mp-r_sp
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                V1=zero
                do i=1,3
                   do j=1,3
                      do k=1,3
                         V1=V1+fifteen*rr(i)*rr(j)*rr(k)*O(i,j,k)/dr**7
                         if(i==j) V1=V1-three*rr(k)*O(i,j,k)/dr**5
                         if(i==k) V1=V1-three*rr(j)*O(i,j,k)/dr**5
                         if(j==k) V1=V1-three*rr(i)*O(i,j,k)/dr**5
                      end do
                   end do
                end do
                V=V + V1/fifteen
             end do
          end do PO
          V_pot_mp(iu_sp)=V_pot_mp(iu_sp)+V
       end do
    end do I_SP

  end subroutine calc_mp_potential
  !*************************************************************

  !*************************************************************
  subroutine deallocate_V_mp()
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_mp
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status
    !------------ Executable code --------------------------------

    deallocate(V_pot_mp,stat=status)
    ASSERT(status==0)

  end subroutine deallocate_V_mp
  !*************************************************************

  !*************************************************************
  subroutine efp_mp_solv_energy(efp_efp_en)
    !  Purpose: Calculate constant part of solvation energy due to
    !           efp multipole centers
    !------------ Modules used -----------------------------------
    use solv_electrostat_module, only: A_matrix_inv
    use solv_cavity_module, only: dielectric_constant,Q_mp
    use potential_module, only: V_pot_mp, N_points, point_in_space
    use pointcharge_module, only : print_energy
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: efp_efp_en
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind), allocatable :: V_buf(:)
    real(r8_kind) :: efp_mpl_solv_energy
    integer(i4_kind) :: status,n_eql
    integer(i4_kind) :: i
    !------------ Executable code --------------------------------

    allocate(V_buf(N_points),stat=status)
    ASSERT(status==0)

    do i=1,N_points
       n_eql=point_in_space(i)%N_equal_points
       V_buf(i)=V_pot_mp(i)/real(n_eql,kind=r8_kind)
    enddo

    allocate(Q_mp(N_points),stat=status)
    ASSERT(status==0)

    Q_mp= -((dielectric_constant-one)/(dielectric_constant))* &
         MATMUL(A_matrix_inv,V_buf)

    deallocate(V_buf,stat=status)
    ASSERT(status==0)

    efp_mpl_solv_energy=0.0_r8_kind
    do i=1,N_points
       efp_mpl_solv_energy=efp_mpl_solv_energy+Q_mp(i)*V_pot_mp(i)*half
    enddo
    if(print_energy) print*,'SOLVATION ENERGY DUE TO EFP MULTIPOLES - 0.5*Q_mp*V_mp:', efp_mpl_solv_energy
    efp_efp_en=efp_efp_en+efp_mpl_solv_energy

  end subroutine efp_mp_solv_energy
  !*************************************************************

  !*************************************************************
  subroutine allocate_E_cav()
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code --------------------------------

    if(.not.allocated(E_cav)) then
       allocate(E_cav(N_ipd), stat=status)
       ASSERT(status==0)
       allocate(E_cav1(N_ipd), stat=status)
       ASSERT(status==0)
    end if

    do i=1,N_ipd
       n_eq=ipd_array(i)%N_equal_dipoles
       if(.not.allocated(E_cav(i)%m)) then
          allocate(E_cav(i)%m(3,n_eq), stat=status)
          ASSERT(status==0)
          allocate(E_cav1(i)%m(3,n_eq), stat=status)
          ASSERT(status==0)
       end if

       E_cav(i)%m=zero
       E_cav1(i)%m=zero
    end do

  end subroutine allocate_E_cav
  !*************************************************************

  !*************************************************************
  subroutine deallocate_E_cav()
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i
    !------------ Executable code --------------------------------

    do i=1,N_ipd
       deallocate(E_cav(i)%m, stat=status)
       ASSERT(status==0)
       deallocate(E_cav1(i)%m, stat=status)
       ASSERT(status==0)
    end do
    deallocate(E_cav, stat=status)
    ASSERT(status==0)
    deallocate(E_cav1, stat=status)
    ASSERT(status==0)

  end subroutine deallocate_E_cav
  !*************************************************************

  !*************************************************************
  subroutine allocate_V_and_Q_id()
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_id, V_pot_id1, N_points
    use solv_cavity_module, only: Q_id, Q_id1
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: status
    !--- executable code-------------------------------------

    if(.not. allocated(V_pot_id)) then
       allocate(V_pot_id(N_points),stat=status)
       ASSERT(status==0)
       V_pot_id=zero
    end if

    if(.not. allocated(V_pot_id1)) then
       allocate(V_pot_id1(N_points),stat=status)
       ASSERT(status==0)
       V_pot_id1=zero
    end if

    if(.not. allocated(Q_id)) then
       allocate(Q_id(N_points),stat=status)
       ASSERT(status==0)
       Q_id=zero
    end if

    if(.not. allocated(Q_id1)) then
       allocate(Q_id1(N_points),stat=status)
       ASSERT(status==0)
       Q_id1=zero
    end if

  end subroutine allocate_V_and_Q_id
  !*************************************************************

  !*************************************************************
  subroutine deallocate_V_id()
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_id, V_pot_id1
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: status
    !--- executable code-------------------------------------

    deallocate(V_pot_id,stat=status)
    ASSERT(status==0)

    deallocate(V_pot_id1,stat=status)
    ASSERT(status==0)

  end subroutine deallocate_V_id
  !*************************************************************

  !*************************************************************
  subroutine calc_V_and_Q_id()
    ! Purpose: Calculate electrostatic potential and point charges
    !          at solute cavity surface points due to efp induced
    !          dipoles
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_id, V_pot_id1, N_points, point_in_space
    use solv_electrostat_module, only: A_matrix_inv
    use solv_cavity_module, only: dielectric_constant, Q_id, Q_id1
    use unique_atom_module, only: N_unique_atoms
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    real(r8_kind) :: V,V1,r_sp(3),r_ip(3),D(3),D1(3),rr(3),dr
    integer(i4_kind) :: iu_sp,n_eq,ie_sp,iu_ip,ie_ip,n_eql,i
    real(r8_kind), allocatable :: V_buf(:)
    integer(i4_kind) :: status
    !--- executable code-------------------------------------

    V_pot_id=zero;  V_pot_id1=zero

    I_SP: do iu_sp=1,N_points

       V=zero; V1=zero
       n_eq=point_in_space(iu_sp)%N_equal_points

       do ie_sp=1,n_eq

          r_sp=point_in_space(iu_sp)%position(:,ie_sp)

          ID: do iu_ip=1,N_ipd
             do ie_ip=1,ipd_array(iu_ip)%N_equal_dipoles

                D=ipd_array(iu_ip)%idipole(:,ie_ip)
                D1=ipd_array(iu_ip)%idipole1(:,ie_ip)

                r_ip=ipd_array(iu_ip)%position(:,ie_ip)
                rr=r_sp-r_ip
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                V=V + dot_product(D,rr)/dr**3
                V1=V1 + dot_product(D1,rr)/dr**3
             end do
          end do ID

          if(N_unique_atoms > 0) then
             V_pot_id(iu_sp)=V_pot_id(iu_sp)+V
             V_pot_id1(iu_sp)=V_pot_id1(iu_sp)+V1
          else
             V_pot_id(iu_sp)=V_pot_id(iu_sp)+V
             V_pot_id1(iu_sp)=V_pot_id1(iu_sp)+V1
          end if
       end do
    end do I_SP

    allocate(V_buf(N_points),stat=status)
    ASSERT(status==0)

    do i=1,N_points
       n_eql=point_in_space(i)%N_equal_points
       V_buf(i)=V_pot_id(i)/real(n_eql,kind=r8_kind)
    enddo

    Q_id= -((dielectric_constant-one)/(dielectric_constant))* &
         MATMUL(A_matrix_inv,V_buf)

    do i=1,N_points
       n_eql=point_in_space(i)%N_equal_points
       V_buf(i)=V_pot_id1(i)/real(n_eql,kind=r8_kind)
    enddo

    Q_id1= -((dielectric_constant-one)/(dielectric_constant))* &
         MATMUL(A_matrix_inv,V_buf)

    deallocate(V_buf,stat=status)
    ASSERT(status==0)

  end subroutine calc_V_and_Q_id
  !*************************************************************

  !*************************************************************
  subroutine calc_E_cav()
    ! Purpose: Calculate electrostatic field produced by the solute
    !          cavity surface charges at positions of induced point
    !          dipoles
    !------------ Modules used -----------------------------------
    use potential_module, only: N_points, point_in_space
    use solv_cavity_module, only: Q_mp, Q_id, Q_id1, Q_n, Q_e
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: iu_id1,ie_id1,n_eq1,iu_sp,n_eq2,ie_sp,i
    real(r8_kind) :: E(3),Q,r_id(3),r_sp(3),rr(3),dr
    integer(i4_kind) :: e_dim,index
    real(r8_kind), pointer :: rotmat(:,:)
    !--- executable code-------------------------------------

    totalsym_field=0.0_r8_kind

    I_PD1: do iu_id1=1,N_ipd

       E=zero
       n_eq1=surface_points(iu_id1)%N_equal_points
       e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
       index=surf_points_grad_index(iu_id1)-1

       do ie_id1=1,n_eq1

          r_id=surface_points(iu_id1)%position(:,ie_id1)
          rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

          I_SP1: do iu_sp=1,N_points
             n_eq2=point_in_space(iu_sp)%N_equal_points

             Q=Q_e(iu_sp)+Q_n(iu_sp)+Q_id(iu_sp)

             do ie_sp=1,n_eq2
                r_sp=point_in_space(iu_sp)%position(:,ie_sp)
                rr=r_id-r_sp
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                E=E-Q*rr/dr**3
             end do
          end do I_SP1

          do i=1,e_dim
             totalsym_field(index+i)=totalsym_field(index+i)+ & !+
                  sum(rotmat(i,:)*E(:))
          end do
       end do
    end do I_PD1

    call transform_to_cart_field(E_cav,totalsym_field)

    totalsym_field=0.0_r8_kind

    I_PD2: do iu_id1=1,N_ipd

       E=zero
       n_eq1=surface_points(iu_id1)%N_equal_points
       e_dim=surf_points_grad_index(iu_id1+1)-surf_points_grad_index(iu_id1)
       index=surf_points_grad_index(iu_id1)-1

       do ie_id1=1,n_eq1

          r_id=surface_points(iu_id1)%position(:,ie_id1)
          rotmat=>surf_points_grad_info(iu_id1)%m(:,:,ie_id1)

          I_SP2: do iu_sp=1,N_points
             n_eq2=point_in_space(iu_sp)%N_equal_points

             Q=Q_e(iu_sp)+Q_n(iu_sp)+Q_id1(iu_sp)

             do ie_sp=1,n_eq2
                r_sp=point_in_space(iu_sp)%position(:,ie_sp)
                rr=r_id-r_sp
                dr=sqrt(dot_product(rr,rr))
                if(dr < 0.001_r8_kind) dr=0.001_r8_kind
                E=E-Q*rr/dr**3
             end do
          end do I_SP2

          do i=1,e_dim
             totalsym_field(index+i)=totalsym_field(index+i)+ & !+
                  sum(rotmat(i,:)*E(:))
          end do
       end do
    end do I_PD2

    call transform_to_cart_field(E_cav1,totalsym_field)

    totalsym_field=0.0_r8_kind

  end subroutine calc_E_cav
  !*************************************************************

  !*************************************************************
  subroutine calc_E_cav1(iu_id1,ie_id1)
    ! Purpose: Calculate electrostatic field produced by the solute
    !          cavity surface charges at positions of induced point
    !          dipoles
    !------------ Modules used -----------------------------------
    use potential_module, only: N_points, point_in_space
    use solv_cavity_module, only: Q_mp, Q_id
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind) :: iu_id1,ie_id1,n_eq1,iu_sp,n_eq2,ie_sp
    real(r8_kind) :: E(3),Q,r_id(3),r_sp(3),rr(3),dr
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    !--- executable code-------------------------------------

    E=zero
    n_eq1=surface_points(iu_id1)%N_equal_points
    r_id=surface_points(iu_id1)%position(:,ie_id1)

    I_SP: do iu_sp=1,N_points
       n_eq2=point_in_space(iu_sp)%N_equal_points

       Q=Q_mp(iu_sp)+Q_id(iu_sp)

       do ie_sp=1,n_eq2
          r_sp=point_in_space(iu_sp)%position(:,ie_sp)
          rr=r_id-r_sp
          dr=sqrt(dot_product(rr,rr))
          if(dr < 0.001_r8_kind) dr=0.001_r8_kind
          E=E+Q*rr/dr**3
       end do
    end do I_SP

    E_cav(iu_id1)%m(:,ie_id1)=E

  end subroutine calc_E_cav1
  !*************************************************************

  !*************************************************************
  subroutine calc_E_cav11(iu_id1,ie_id1)
    ! Purpose: Calculate electrostatic field produced by the solute
    !          cavity surface charges at positions of induced point
    !          dipoles
    !------------ Modules used -----------------------------------
    use potential_module, only: N_points, point_in_space
    use solv_cavity_module, only: Q_mp, Q_id1
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind) :: iu_id1,ie_id1,n_eq1,iu_sp,n_eq2,ie_sp
    real(r8_kind) :: E(3),Q,r_id(3),r_sp(3),rr(3),dr
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    !--- executable code-------------------------------------

    E=zero
    n_eq1=surface_points(iu_id1)%N_equal_points
    r_id=surface_points(iu_id1)%position(:,ie_id1)

    I_SP: do iu_sp=1,N_points
       n_eq2=point_in_space(iu_sp)%N_equal_points

       Q=Q_mp(iu_sp)+Q_id1(iu_sp)

       do ie_sp=1,n_eq2
          r_sp=point_in_space(iu_sp)%position(:,ie_sp)
          rr=r_id-r_sp
          dr=sqrt(dot_product(rr,rr))
          if(dr < 0.001_r8_kind) dr=0.001_r8_kind
          E=E+Q*rr/dr**3
       end do
    end do I_SP

    E_cav1(iu_id1)%m(:,ie_id1)=E

  end subroutine calc_E_cav11
  !*************************************************************

  !*************************************************************
  subroutine efp_id_mp_solv_energy(en_solv)
    ! Purpose: Calculate polarizable efp-efp energy
    !------------ Modules used -----------------------------------
    use potential_module, only: V_pot_id, V_pot_mp, N_points
    use solv_cavity_module, only: Q_mp, Q_id
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(out) :: en_solv
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i
    !--- executable code-------------------------------------

    en_solv=zero
    if(.not.do_pol_pcm) return
    do i=1,N_points
       en_solv=en_solv+Q_id(i)*V_pot_mp(i)*half
!!$       en_solv=en_solv+(Q_id(i)*V_pot_mp(i)+Q_mp(i)*V_pot_id(i)+Q_id(i)*V_pot_id(i))*half
    enddo

  end subroutine efp_id_mp_solv_energy
  !*************************************************************

  !--------------- End of module ----------------------------------
end module efp_solv_module

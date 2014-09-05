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
module efp_efp_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Calculate interfragment interactions
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AS
  !  Date: 11.07
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
#include <def.h>

  use type_module ! type specification parameters
  use common_data_module
  use pointcharge_module, only : print_energy
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  real(r8_kind), public :: efp_efp_en=zero
  real(r8_kind), public :: qm_efp_energy
  !------------ Interface statements ---------------------------------
  !------------ public functions and subroutines ---------------------
  public efp_efp_energy, efp_efp_gradients

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of types ---------------------------------
  !------------ Declaration of constants and variables ---------------
  integer(i4_kind), parameter :: energy=1
  integer(i4_kind), parameter :: gradient=2
  logical :: print_id_scc

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains
  !*************************************************************
  subroutine efp_efp_energy(print_id_conv)
    !  Purpose: Calculate energies of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N
    use point_dqo_module, only: N_pd,N_pq,N_po
    use induced_dipoles_module, only: N_ipd
    use efp_rep_module, only: N_rcf
    use unique_atom_module, only: N_unique_atoms
    use efp_polar_module, only: calc_polar
    use operations_module, only: operations_solvation_effect
    use efp_solv_module, only: calc_mp_potential,efp_mp_solv_energy,deallocate_V_mp
    !------------ Declaration of formal parameters ---------------
    logical, intent(in) :: print_id_conv
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    print_id_scc=print_id_conv

    efp_efp_en=zero

    if(pointcharge_N > 0) then
       !Charge-charge
       call efp_efp_pc_pc(energy)

       !Charge-dipole
       if(N_pd > 0) then
          call efp_efp_pc_pd(energy)
       end if
       !Charge-quadrupole
       if(N_pq > 0) then
          call efp_efp_pc_pq(energy)
       end if
       !Charge-octopole
       if(N_po > 0) then
          call efp_efp_pc_po(energy)
       end if
    end if

    if(N_pd > 0) then
       !Dipole-dipole
       call efp_efp_pd_pd(energy)

       !Dipole-quadrupole
       if(N_pq > 0) then
          call efp_efp_pd_pq(energy)
       end if
       !Dipole-octopole
       if(N_po > 0) then
       end if
    end if

    if(N_pq > 0) then
       !Quadrupole-quadrupole
       call efp_efp_pq_pq(energy)

       !Quadrupole-octopole
       if(N_po > 0) then
       end if
    end if

    if(N_po > 0) then
       !Octopole-octopole
    end if

    if(N_rcf > 0) then
       !Repulsive interaction
       call efp_efp_repulsion(energy)
    end if

    if(operations_solvation_effect) then
       call calc_mp_potential()
       if(N_unique_atoms == 0) call efp_mp_solv_energy(efp_efp_en)
    end if

    if((N_ipd > 0  .and. N_unique_atoms == 0) .or.  &
       (calc_polar .and. N_ipd > 0)) then
       !Polarization
       call efp_efp_polar(energy)
    end if

    if(operations_solvation_effect .and. N_unique_atoms == 0) &
         call deallocate_V_mp()

    if(print_energy) print*,'EFP-EFP energy',efp_efp_en

  end subroutine efp_efp_energy
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_gradients()
    !  Purpose: Calculate gradients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, efp_fixed
    use point_dqo_module, only: N_pd,N_pq,N_po
    use induced_dipoles_module, only: N_ipd
    use efp_rep_module, only: N_rcf
    use unique_atom_module, only: N_unique_atoms
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------

    if(efp_fixed) return

    if(pointcharge_N > 0) then
       !Charge-charge
       call efp_efp_pc_pc(gradient)

       !Charge-dipole
       if(N_pd > 0) then
          call efp_efp_pc_pd(gradient)
       end if

       !Charge-quadrupole
       if(N_pq > 0) then
          call efp_efp_pc_pq(gradient)
       end if

       !Charge-octopole
       if(N_po > 0) then
          call efp_efp_pc_po(gradient)
       end if
    end if

    if(N_pd > 0) then
       !Dipole-dipole
       call efp_efp_pd_pd(gradient)

       !Dipole-quadrupole
       if(N_pq > 0) then
          call efp_efp_pd_pq(gradient)
       end if

       !Dipole-octopole
       if(N_po > 0) then
       end if
    end if

    if(N_pq > 0) then
       !Quadrupole-quadrupole
       call efp_efp_pq_pq(gradient)

       !Quadrupole-octopole
       if(N_po > 0) then
       end if
    end if

    if(N_po > 0) then
       !Octopole-octopole
    end if

    if(N_rcf > 0) then
       !Repulsive interaction
       call efp_efp_repulsion(gradient)
    end if

    if(N_ipd > 0) then
       !Polarization
       call efp_efp_polar(gradient)
    end if

  end subroutine efp_efp_gradients
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
  subroutine efp_efp_repulsion(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used -----------------------------------
    use efp_rep_module, only: N_rcf, rcf_array
    use efp_rep_module, only: repf_grad_info, repf_grad_index
    use efp_rep_module, only: gradient_repf_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i1,i2
    integer(i4_kind) :: n_equal1,n_equal2,ind1,ind2
    integer(i4_kind) :: i_grp1,i_grp2,ityp
    real(r8_kind) :: r1(3),r2(3),r12(3),dr,A,C,E
    real(r8_kind) :: efp_efp_rep_energy,gradient1(3),gradient2(3)
    logical :: do_energy, do_grads
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_rep_energy=zero

    ind1=0
    cu1: do i_unique1=1,N_rcf

       if(do_grads) then
          grad_dim1=repf_grad_index(i_unique1+1)-repf_grad_index(i_unique1)
          index1=repf_grad_index(i_unique1)-1
       end if

       n_equal1=rcf_array(i_unique1)%N_equal_centers
       ce1: do i_equal1=1,n_equal1

          if(do_grads) rotmat1=>repf_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero

          r1=rcf_array(i_unique1)%position(:,i_equal1)

          ind1=ind1+1
          i_grp1=rcf_array(i_unique1)%group(i_equal1)

          ind2=0
          cu2 :do i_unique2=1,N_rcf

             if(do_grads) then
                grad_dim2=repf_grad_index(i_unique2+1)-repf_grad_index(i_unique2)
                index2=repf_grad_index(i_unique2)-1
             end if

             n_equal2=rcf_array(i_unique2)%N_equal_centers
             ce2: do i_equal2=1,n_equal2

                ind2=ind2+1
                if(ind2 <= ind1) cycle ce2
                i_grp2=rcf_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>repf_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero

                r2=rcf_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                ityp=rcf_array(i_unique2)%type(i_equal2)
                C=rcf_array(i_unique1)%C(ityp)
                A=rcf_array(i_unique1)%A(ityp)

                E=C*exp(-A*dr)
                if(do_energy) efp_efp_rep_energy=efp_efp_rep_energy+E

                if(do_grads) then
                   gradient1=gradient1-A*E*r12/dr

                   gradient2=gradient2+A*E*r12/dr

                   do i2=1,grad_dim2
                      gradient_repf_totalsym(index2+i2)=gradient_repf_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_repf_totalsym(index1+i1)=gradient_repf_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
             enddo
          end if

       end do ce1
    end do cu1

    if(do_energy .and. print_energy) print*,'EFP-EFP REP ENERGY:', efp_efp_rep_energy

    efp_efp_en=efp_efp_en+efp_efp_rep_energy

  end subroutine efp_efp_repulsion
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pc_pc(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use pointcharge_module, only: unique_pc_grad_info,unique_pc_index
    use pointcharge_module, only: gradient_pc_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i1,i2
    integer(i4_kind) :: n_equal1,n_equal2,ind1,ind2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: Z1,Z2,A1,C1,A2,C2,E,G(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    logical :: do_energy, do_grads
    real(r8_kind) :: efp_efp_pc_pc_energy,gradient1(3),gradient2(3)
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pc_pc_energy=zero

    ind1=0
    cu1: do i_unique1=1,pointcharge_N

       if(do_grads) then
          grad_dim1=unique_pc_index(i_unique1+1)-unique_pc_index(i_unique1)
          index1=unique_pc_index(i_unique1)-1
       end if

       Z1=pointcharge_array(i_unique1)%Z
       C1=pointcharge_array(i_unique1)%Cf
       A1=pointcharge_array(i_unique1)%Af
       n_equal1=pointcharge_array(i_unique1)%N_equal_charges

       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>unique_pc_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero

          r1=pointcharge_array(i_unique1)%position(:,i_equal1)

          ind1=ind1+1
          i_grp1=pointcharge_array(i_unique1)%group(i_equal1)

          ind2=0
          cu2 :do i_unique2=1,pointcharge_N

             if(do_grads) then
                grad_dim2=unique_pc_index(i_unique2+1)-unique_pc_index(i_unique2)
                index2=unique_pc_index(i_unique2)-1
             end if

             Z2=pointcharge_array(i_unique2)%Z
             C2=pointcharge_array(i_unique2)%Cf
             A2=pointcharge_array(i_unique2)%Af
             n_equal2=pointcharge_array(i_unique2)%N_equal_charges

             ce2: do i_equal2=1,n_equal2

                ind2=ind2+1
                if(ind2 <= ind1) cycle ce2
                i_grp2=pointcharge_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>unique_pc_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero

                r2=pointcharge_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                E=Z1*Z2/dr
                if(C1 /= zero .and. C2 == zero) E=E-Z1*Z2*exp(-A1*dr)/dr
                if(C2 /= zero .and. C1 == zero) E=E-Z1*Z2*exp(-A2*dr)/dr
                if(C2 /= zero .and. C1 /= zero) then 
                   if(A1 /= A2) E=E-Z1*Z2*(A1**2*exp(-A2*dr)-A2**2*exp(-A1*dr))/ &
                     (dr*(A1**2-A2**2))
                   if(A1 == A2) E=E-Z1*Z2*(one+A1*dr/two)*exp(-A1*dr)/dr
                end if
                if(do_energy) efp_efp_pc_pc_energy=efp_efp_pc_pc_energy+E

                if(do_grads) then
                   G=-Z1*Z2*r12/dr**3
                   if(C1 /= zero .and. C2 == zero) G=G+Z1*Z2*exp(-A1*dr)*r12*(one/dr**3+A1/dr**2)
                   if(C2 /= zero .and. C1 == zero) G=G+Z1*Z2*exp(-A2*dr)*r12*(one/dr**3+A2/dr**2)
                   if(C2 /= zero .and. C1 /= zero) then 
                      if(A1 /= A2) G=G+Z1*Z2*(A1**2*exp(-A2*dr)-A2**2*exp(-A1*dr))*r12/ &
                           (dr**3*(A1**2-A2**2))+ &
                           Z1*Z2*(A1**2*A2*exp(-A2*dr)-A2**2*A1*exp(-A1*dr))*r12/(dr**2*(A1**2-A2**2))
                      if(A1 == A2) G=G+Z1*Z2*exp(-A1*dr)*r12*(one/dr**3+A1/dr**2+A1**2/(dr*two))
                   end if

                   gradient1=gradient1+G

                   gradient2=gradient2-G

                   do i2=1,grad_dim2
                      gradient_pc_totalsym(index2+i2)=gradient_pc_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_pc_totalsym(index1+i1)=gradient_pc_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP CHARGE-CHARGE ENERGY:', efp_efp_pc_pc_energy

    efp_efp_en=efp_efp_en+efp_efp_pc_pc_energy

  end subroutine efp_efp_pc_pc
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pc_pd(to_do)
    !  Purpose: Calculate energy and gradients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use pointcharge_module, only: unique_pc_grad_info,unique_pc_index
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: dipoles_grad_index,dipoles_grad_info
    use point_dqo_module, only: gradient_dip_totalsym,torque_dip_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i1,i2
    integer(i4_kind) :: n_equal1,n_equal2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: Z1,D2(3),G(3),E(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    logical :: do_energy, do_grads
    real(r8_kind) :: efp_efp_pc_pd_energy,gradient1(3),gradient2(3),torque2(3)
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pc_pd_energy=zero

    cu1: do i_unique1=1,pointcharge_N

       if(do_grads) then
          grad_dim1=unique_pc_index(i_unique1+1)-unique_pc_index(i_unique1)
          index1=unique_pc_index(i_unique1)-1
       end if

       n_equal1=pointcharge_array(i_unique1)%N_equal_charges
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>unique_pc_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero

          r1=pointcharge_array(i_unique1)%position(:,i_equal1)
          Z1=pointcharge_array(i_unique1)%Z

          i_grp1=pointcharge_array(i_unique1)%group(i_equal1)

          cu2 :do i_unique2=1,N_pd

             if(do_grads) then
                grad_dim2=dipoles_grad_index(i_unique2+1)-dipoles_grad_index(i_unique2)
                index2=dipoles_grad_index(i_unique2)-1
             end if

             n_equal2=pd_array(i_unique2)%N_equal_dipoles
             ce2: do i_equal2=1,n_equal2

                i_grp2=pd_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>dipoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=pd_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                D2=pd_array(i_unique2)%dipole(:,i_equal2)
                
                E=Z1*r12/dr**3

                if(do_energy) efp_efp_pc_pd_energy=efp_efp_pc_pd_energy+dot_product(D2,E)

                if(do_grads) then
                   G=-three*Z1*dot_product(D2,r12)*r12/dr**5+Z1*D2/dr**3

                   gradient1=gradient1+G

                   gradient2=gradient2-G

                   torque2=torque2+vector_product(D2,E)

                   do i2=1,grad_dim2
                      gradient_dip_totalsym(index2+i2)=gradient_dip_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_dip_totalsym(index2+i2)=torque_dip_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_pc_totalsym(index1+i1)=gradient_pc_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP CHARGE-DIPOLE ENERGY:', efp_efp_pc_pd_energy

    efp_efp_en=efp_efp_en+efp_efp_pc_pd_energy

  end subroutine efp_efp_pc_pd
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pc_pq(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use pointcharge_module, only: unique_pc_grad_info,unique_pc_index
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: qpoles_grad_index,qpoles_grad_info
    use point_dqo_module, only: gradient_quad_totalsym,torque_quad_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i,j,k
    integer(i4_kind) :: n_equal1,n_equal2,i1,i2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: Z1,Q2(3,3),Eij,G(3),Qr(3),EE(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    real(r8_kind) :: efp_efp_pc_pq_energy,gradient1(3),gradient2(3),torque2(3)
    logical :: do_energy, do_grads
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pc_pq_energy=zero

    cu1: do i_unique1=1,pointcharge_N

       if(do_grads) then
          grad_dim1=unique_pc_index(i_unique1+1)-unique_pc_index(i_unique1)
          index1=unique_pc_index(i_unique1)-1
       end if

       n_equal1=pointcharge_array(i_unique1)%N_equal_charges
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>unique_pc_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero

          r1=pointcharge_array(i_unique1)%position(:,i_equal1)
          Z1=pointcharge_array(i_unique1)%Z

          i_grp1=pointcharge_array(i_unique1)%group(i_equal1)

          cu2 :do i_unique2=1,N_pq

             if(do_grads) then
                grad_dim2=qpoles_grad_index(i_unique2+1)-qpoles_grad_index(i_unique2)
                index2=qpoles_grad_index(i_unique2)-1
             end if

             n_equal2=pq_array(i_unique2)%N_equal_qpoles
             ce2: do i_equal2=1,n_equal2

                i_grp2=pq_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>qpoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=pq_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                Q2=pq_array(i_unique2)%quadrupole(:,:,i_equal2)
                
                Eij=zero
                G=zero
                do i=1,3
                   do j=1,3
                      Eij=Eij+three*r12(i)*r12(j)*Q2(i,j)/dr**5
                      if(i==j) Eij=Eij-Q2(i,j)/dr**3

                      if(do_grads) then
                         G=G-(fifteen*r12(i)*r12(j)*r12/dr**7)*Q2(i,j)
                         do k=1,3
                            if(k==i) G(k)=G(k)+(three*r12(j)/dr**5)*Q2(i,j)
                            if(k==j) G(k)=G(k)+(three*r12(i)/dr**5)*Q2(i,j)
                         end do
                         if(i==j) G=G+(three*r12/dr**5)*Q2(i,j)
                      end if
                   end do
                end do
                if(do_energy) efp_efp_pc_pq_energy=efp_efp_pc_pq_energy+ &
                     Z1*Eij/three

                if(do_grads) then
                   gradient1=gradient1+Z1*G/three

                   gradient2=gradient2-Z1*G/three

                   Qr=matmul(Q2,r12)
                   EE=two*Z1*r12/dr**5
                   torque2=torque2+vector_product(Qr,EE)

                   do i2=1,grad_dim2
                      gradient_quad_totalsym(index2+i2)=gradient_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_quad_totalsym(index2+i2)=torque_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_pc_totalsym(index1+i1)=gradient_pc_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy)print*,'EFP-EFP CHARGE-QUADRUPOLE ENERGY:', efp_efp_pc_pq_energy

    efp_efp_en=efp_efp_en+efp_efp_pc_pq_energy

  end subroutine efp_efp_pc_pq
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pc_po(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use pointcharge_module, only: pointcharge_N, pointcharge_array
    use pointcharge_module, only: unique_pc_grad_info,unique_pc_index
    use pointcharge_module, only: gradient_pc_totalsym
    use point_dqo_module, only: po_array, N_po
    use point_dqo_module, only: opoles_grad_index,opoles_grad_info
    use point_dqo_module, only: gradient_oct_totalsym,torque_oct_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i,j,k,l
    integer(i4_kind) :: n_equal1,n_equal2,i1,i2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: Z1,O2(3,3,3),Eijk,G(3),EE(3),Orr(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    logical :: do_energy, do_grads
    real(r8_kind) :: efp_efp_pc_po_energy,gradient1(3),gradient2(3),torque2(3)
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pc_po_energy=zero

    cu1: do i_unique1=1,pointcharge_N

       if(do_grads) then
          grad_dim1=unique_pc_index(i_unique1+1)-unique_pc_index(i_unique1)
          index1=unique_pc_index(i_unique1)-1
       end if

       n_equal1=pointcharge_array(i_unique1)%N_equal_charges
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>unique_pc_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero

          r1=pointcharge_array(i_unique1)%position(:,i_equal1)
          Z1=pointcharge_array(i_unique1)%Z

          i_grp1=pointcharge_array(i_unique1)%group(i_equal1)

          cu2 :do i_unique2=1,N_po

             if(do_grads) then
                grad_dim2=opoles_grad_index(i_unique2+1)-opoles_grad_index(i_unique2)
                index2=opoles_grad_index(i_unique2)-1
             end if

             n_equal2=po_array(i_unique2)%N_equal_opoles
             ce2: do i_equal2=1,n_equal2

                i_grp2=po_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>opoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=po_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                O2=po_array(i_unique2)%octopole(:,:,:,i_equal2)
             
                Eijk=zero
                G=zero
                Orr=zero
                do i=1,3
                   do j=1,3
                      do k=1,3
                         Eijk=Eijk+fifteen*r12(i)*r12(j)*r12(k)*O2(i,j,k)/dr**7
                         if(i==j) Eijk=Eijk-three*r12(k)*O2(i,j,k)/dr**5
                         if(i==k) Eijk=Eijk-three*r12(j)*O2(i,j,k)/dr**5
                         if(j==k) Eijk=Eijk-three*r12(i)*O2(i,j,k)/dr**5

                         if(do_grads) then
                            G=G-(105.0_r8_kind*r12(i)*r12(j)*r12(k)/dr**9)*r12*O2(i,j,k)
                            do l=1,3
                               if(l==i) G(l)=G(l)+(fifteen*r12(j)*r12(k)/dr**7)*O2(i,j,k)
                               if(l==j) G(l)=G(l)+(fifteen*r12(i)*r12(k)/dr**7)*O2(i,j,k)
                               if(l==k) G(l)=G(l)+(fifteen*r12(i)*r12(j)/dr**7)*O2(i,j,k)

                               if(l==i) Orr(l)=Orr(l)+O2(i,j,k)*r12(j)*r12(k)
                               if(l==j) Orr(l)=Orr(l)+O2(i,j,k)*r12(i)*r12(k)
                               if(l==k) Orr(l)=Orr(l)+O2(i,j,k)*r12(i)*r12(j)
                            end do
                            if(i==j) G=G+(fifteen*r12(k)*r12/dr**7)*O2(i,j,k)
                            if(i==k) G=G+(fifteen*r12(j)*r12/dr**7)*O2(i,j,k)
                            if(j==k) G=G+(fifteen*r12(i)*r12/dr**7)*O2(i,j,k)
                            do l=1,3
                               if(i==j .and. k==l) G(l)=G(l)-(three*O2(i,j,k)/dr**5)
                               if(i==k .and. j==l) G(l)=G(l)-(three*O2(i,j,k)/dr**5)
                               if(j==k .and. i==l) G(l)=G(l)-(three*O2(i,j,k)/dr**5)
                            end do
                         end if
                      end do
                   end do
                end do
                if(do_energy) efp_efp_pc_po_energy=efp_efp_pc_po_energy+ &
                     Z1*Eijk/fifteen

                if(do_grads) then
                   gradient1=gradient1+Z1*G/fifteen

                   gradient2=gradient2-Z1*G/fifteen

                   EE=Z1*r12/dr**7
                   torque2=torque2+vector_product(Orr,EE)

                   do i2=1,grad_dim2
                      gradient_oct_totalsym(index2+i2)=gradient_oct_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_oct_totalsym(index2+i2)=torque_oct_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_pc_totalsym(index1+i1)=gradient_pc_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP CHARGE-OCTOPOLE ENERGY:', efp_efp_pc_po_energy

    efp_efp_en=efp_efp_en+efp_efp_pc_po_energy

  end subroutine efp_efp_pc_po
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pd_pd(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used -----------------------------------
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: dipoles_grad_index,dipoles_grad_info
    use point_dqo_module, only: gradient_dip_totalsym,torque_dip_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i1,i2
    integer(i4_kind) :: n_equal1,n_equal2,ind1,ind2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: D1(3),D2(3),G(3),E1(3),E2(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    real(r8_kind) :: efp_efp_pd_pd_energy,gradient1(3),gradient2(3),torque1(3),torque2(3)
    logical :: do_energy, do_grads
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pd_pd_energy=zero

    ind1=0
    cu1: do i_unique1=1,N_pd

       if(do_grads) then
          grad_dim1=dipoles_grad_index(i_unique1+1)-dipoles_grad_index(i_unique1)
          index1=dipoles_grad_index(i_unique1)-1
       end if

       n_equal1=pd_array(i_unique1)%N_equal_dipoles
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>dipoles_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero
          torque1=zero

          r1=pd_array(i_unique1)%position(:,i_equal1)
          D1=pd_array(i_unique1)%dipole(:,i_equal1)

          ind1=ind1+1
          i_grp1=pd_array(i_unique1)%group(i_equal1)

          ind2=0
          cu2 :do i_unique2=1,N_pd

             if(do_grads) then
                grad_dim2=dipoles_grad_index(i_unique2+1)-dipoles_grad_index(i_unique2)
                index2=dipoles_grad_index(i_unique2)-1
             end if

             n_equal2=pd_array(i_unique2)%N_equal_dipoles
             ce2: do i_equal2=1,n_equal2

                ind2=ind2+1
                if(ind2 <= ind1) cycle ce2
                i_grp2=pd_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>dipoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=pd_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                D2=pd_array(i_unique2)%dipole(:,i_equal2)

                if(do_energy) efp_efp_pd_pd_energy=efp_efp_pd_pd_energy+dot_product(D1,D2)/dr**3- &
                     three*dot_product(D1,r12)*dot_product(D2,r12)/dr**5

                if(do_grads) then
                   G=-three*dot_product(D1,D2)*r12/dr**5- &
                        three*(D1*dot_product(D2,r12)+D2*dot_product(D1,r12))/dr**5+ &
                        fifteen*dot_product(D1,r12)*dot_product(D2,r12)*r12/dr**7

                   gradient1=gradient1+G

                   gradient2=gradient2-G

                   E2=D2/dr**3-three*r12*dot_product(D2,r12)/dr**5
                   E1=D1/dr**3-three*r12*dot_product(D1,r12)/dr**5

                   torque1=torque1+vector_product(D1,E2)

                   torque2=torque2+vector_product(D2,E1)
                
                   do i2=1,grad_dim2
                      gradient_dip_totalsym(index2+i2)=gradient_dip_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_dip_totalsym(index2+i2)=torque_dip_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_dip_totalsym(index1+i1)=gradient_dip_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
                torque_dip_totalsym(index1+i1)=torque_dip_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*torque1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP DIPOLE-DIPOLE ENERGY:', efp_efp_pd_pd_energy

    efp_efp_en=efp_efp_en+efp_efp_pd_pd_energy

  end subroutine efp_efp_pd_pd
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pd_pq(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use point_dqo_module, only: pd_array, N_pd
    use point_dqo_module, only: dipoles_grad_index,dipoles_grad_info
    use point_dqo_module, only: gradient_dip_totalsym,torque_dip_totalsym
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: qpoles_grad_index,qpoles_grad_info
    use point_dqo_module, only: gradient_quad_totalsym,torque_quad_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i,j,k
    integer(i4_kind) :: n_equal1,n_equal2,i1,i2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: D1(3),Q2(3,3),Eij,G(3),E2(3),Qr(3),EE(3)
    real(r8_kind) :: r1(3),r2(3),r12(3),dr
    real(r8_kind) :: efp_efp_pd_pq_energy,gradient1(3),gradient2(3),torque1(3),torque2(3)
    logical :: do_energy, do_grads
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pd_pq_energy=zero

    cu1: do i_unique1=1,N_pd

       if(do_grads) then
          grad_dim1=dipoles_grad_index(i_unique1+1)-dipoles_grad_index(i_unique1)
          index1=dipoles_grad_index(i_unique1)-1
       end if

       n_equal1=pd_array(i_unique1)%N_equal_dipoles
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>dipoles_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero
          torque1=zero

          r1=pd_array(i_unique1)%position(:,i_equal1)
          D1=pd_array(i_unique1)%dipole(:,i_equal1)

          i_grp1=pd_array(i_unique1)%group(i_equal1)

          cu2 :do i_unique2=1,N_pq

             if(do_grads) then
                grad_dim2=qpoles_grad_index(i_unique2+1)-qpoles_grad_index(i_unique2)
                index2=qpoles_grad_index(i_unique2)-1
             end if

             n_equal2=pq_array(i_unique2)%N_equal_qpoles
             ce2: do i_equal2=1,n_equal2

                i_grp2=pd_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>qpoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=pq_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                Q2=pq_array(i_unique2)%quadrupole(:,:,i_equal2)
                
                Eij=zero
                G=zero
                E2=zero
                EE=zero
                do i=1,3
                   do j=1,3
                      Eij=Eij+(three*(D1(i)*r12(j)+D1(j)*r12(i))/dr**5- &
                           fifteen*r12(i)*r12(j)*dot_product(r12,D1)/dr**7)*Q2(i,j)

                      if(do_grads) then
                         G=G-(fifteen*(D1(i)*r12(j)+D1(j)*r12(i))*r12/dr**7)*Q2(i,j)
                         do k=1,3
                            if(k==i) G(k)=G(k)+(three*D1(j)/dr**5)*Q2(i,j)
                            if(k==j) G(k)=G(k)+(three*D1(i)/dr**5)*Q2(i,j)
                         end do
                         G=G+(105.0_r8_kind*r12(i)*r12(j)*dot_product(r12,D1)*r12/dr**9)*Q2(i,j)
                         do k=1,3
                            if(k==i) G(k)=G(k)-(fifteen*r12(j)*dot_product(r12,D1)/dr**7)*Q2(i,j)
                            if(k==j) G(k)=G(k)-(fifteen*r12(i)*dot_product(r12,D1)/dr**7)*Q2(i,j)
                         end do
                         G=G-(fifteen*r12(i)*r12(j)*D1/dr**7)*Q2(i,j)

                         E2=E2-fifteen*r12(i)*r12(j)*r12*Q2(i,j)/dr**7
                         do k=1,3
                            if(i==k) E2(k)=E2(k)+three*r12(j)*Q2(i,j)/dr**5
                            if(j==k) E2(k)=E2(k)+three*r12(i)*Q2(i,j)/dr**5
                         end do
                      end if
                   end do
                   if(do_grads) then
                      EE(1)=EE(1)+(Q2(2,i)*r12(3)-Q2(3,i)*r12(2))*D1(i)/dr**5
                      EE(2)=EE(2)+(Q2(3,i)*r12(1)-Q2(1,i)*r12(3))*D1(i)/dr**5
                      EE(3)=EE(3)+(Q2(1,i)*r12(2)-Q2(2,i)*r12(1))*D1(i)/dr**5
                   end if
                end do
                if(do_energy) efp_efp_pd_pq_energy=efp_efp_pd_pq_energy+Eij/three

                if(do_grads) then
                   gradient1=gradient1+G/three

                   gradient2=gradient2-G/three

                   E2=E2/three
                   torque1=torque1+vector_product(D1,E2)

                   Qr=matmul(Q2,r12)
                   torque2=torque2+two*EE
                   EE=two*(D1/dr**5-five*dot_product(D1,r12)*r12/dr**7)
                   torque2=torque2+vector_product(Qr,EE)

                   do i2=1,grad_dim2
                      gradient_quad_totalsym(index2+i2)=gradient_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_quad_totalsym(index2+i2)=torque_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if

             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_dip_totalsym(index1+i1)=gradient_dip_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
                torque_dip_totalsym(index1+i1)=torque_dip_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*torque1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP DIPOLE-QUADRUPOLE ENERGY:', efp_efp_pd_pq_energy

    efp_efp_en=efp_efp_en+efp_efp_pd_pq_energy

  end subroutine efp_efp_pd_pq
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_pq_pq(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use point_dqo_module, only: pq_array, N_pq
    use point_dqo_module, only: qpoles_grad_index,qpoles_grad_info
    use point_dqo_module, only: gradient_quad_totalsym,torque_quad_totalsym
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: i_unique1,i_equal1,i_unique2,i_equal2,i,j,k,l,m,n,i1,i2
    integer(i4_kind) :: n_equal1,n_equal2,ind1,ind2
    integer(i4_kind) :: i_grp1,i_grp2
    real(r8_kind) :: Q1(3,3),Q2(3,3),Eij,Eijk,Eijkl,Gij(3),Gijk(3),Gijkl(3),Qr1(3),Qr2(3),Qrr1,Qrr2
    real(r8_kind) :: r1(3),r2(3),r12(3),dr,T1(3),T2(3)
    real(r8_kind) :: efp_efp_pq_pq_energy,gradient1(3),gradient2(3),torque1(3),torque2(3)
    logical :: do_energy, do_grads
    real(r8_kind), pointer :: rotmat1(:,:),rotmat2(:,:)
    integer(i4_kind) :: grad_dim1,index1,grad_dim2,index2
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    efp_efp_pq_pq_energy=zero

    ind1=0
    cu1: do i_unique1=1,N_pq

       if(do_grads) then
          grad_dim1=qpoles_grad_index(i_unique1+1)-qpoles_grad_index(i_unique1)
          index1=qpoles_grad_index(i_unique1)-1
       end if

       n_equal1=pq_array(i_unique1)%N_equal_qpoles
       ce1: do i_equal1=1,n_equal1
          if(do_grads) rotmat1=>qpoles_grad_info(i_unique1)%m(:,:,i_equal1)
          gradient1=zero
          torque1=zero

          r1=pq_array(i_unique1)%position(:,i_equal1)
          Q1=pq_array(i_unique1)%quadrupole(:,:,i_equal1)

          ind1=ind1+1
          i_grp1=pq_array(i_unique1)%group(i_equal1)

          ind2=0
          cu2 :do i_unique2=1,N_pq

             if(do_grads) then
                grad_dim2=qpoles_grad_index(i_unique2+1)-qpoles_grad_index(i_unique2)
                index2=qpoles_grad_index(i_unique2)-1
             end if

             n_equal2=pq_array(i_unique2)%N_equal_qpoles
             ce2: do i_equal2=1,n_equal2

                ind2=ind2+1
                if(ind2 <= ind1) cycle ce2
                i_grp2=pq_array(i_unique2)%group(i_equal2)
                if(i_grp2 == i_grp1) cycle ce2

                if(do_grads) rotmat2=>qpoles_grad_info(i_unique2)%m(:,:,i_equal2)
                gradient2=zero
                torque2=zero

                r2=pq_array(i_unique2)%position(:,i_equal2)
                r12=r1-r2
                dr=sqrt(dot_product(r12,r12))
                if(dr <= 0.00001_r8_kind) dr=0.00001_r8_kind

                Q2=pq_array(i_unique2)%quadrupole(:,:,i_equal2)

                if(do_grads) then
                   Qr1=matmul(Q1,r12)
                   Qrr1=dot_product(r12,Qr1)
                   Qr2=matmul(Q2,r12)
                   Qrr2=dot_product(r12,Qr2)
                end if

                Eij=zero; Eijk=zero; Eijkl=zero
                Gijk=zero; Gijkl=zero
                T1=zero; T2=zero
                do i=1,3
                   do j=1,3
                      Eij=Eij+Q1(i,j)*Q2(i,j)
                      do k=1,3
                         Eijk=Eijk+Q1(i,j)*Q2(i,k)*r12(j)*r12(k)

                         if(do_grads) then
                            do m=1,3
                               if(m==j) Gijk(m)=Gijk(m)+r12(k)*Q1(i,j)*Q2(i,k)
                               if(m==k) Gijk(m)=Gijk(m)+r12(j)*Q1(i,j)*Q2(i,k)
                            end do
                         end if

                         do l=1,3
                            Eijkl=Eijkl+Q1(i,j)*Q2(k,l)*r12(i)*r12(j)*r12(k)*r12(l)

                            if(do_grads) then
                               do n=1,3
                                  if(n==i) Gijkl(n)=Gijkl(n)+r12(j)*r12(k)*r12(l)*Q1(i,j)*Q2(k,l)
                                  if(n==j) Gijkl(n)=Gijkl(n)+r12(i)*r12(k)*r12(l)*Q1(i,j)*Q2(k,l)
                                  if(n==k) Gijkl(n)=Gijkl(n)+r12(i)*r12(j)*r12(l)*Q1(i,j)*Q2(k,l)
                                  if(n==l) Gijkl(n)=Gijkl(n)+r12(i)*r12(j)*r12(k)*Q1(i,j)*Q2(k,l)
                               end do
                            end if
                         end do
                      end do
                   end do

                   if(do_grads) then
                      T1(1)=T1(1)+Q1(2,i)*(-30.0_r8_kind*(r12(3)*Qr2(i)+r12(i)*Qr2(3)))/dr**7- &
                                  Q1(3,i)*(-30.0_r8_kind*(r12(2)*Qr2(i)+r12(i)*Qr2(2)))/dr**7+ &
                                  Q1(2,i)*(105.0_r8_kind*r12(i)*r12(3)*Qrr2)/dr**9- &
                                  Q1(3,i)*(105.0_r8_kind*r12(i)*r12(2)*Qrr2)/dr**9+ &
                                  Q1(2,i)*six*Q2(i,3)/dr**5- &
                                  Q1(3,i)*six*Q2(i,2)/dr**5

                      T1(2)=T1(2)+Q1(3,i)*(-30.0_r8_kind*(r12(1)*Qr2(i)+r12(i)*Qr2(1)))/dr**7- &
                                  Q1(1,i)*(-30.0_r8_kind*(r12(3)*Qr2(i)+r12(i)*Qr2(3)))/dr**7+ &
                                  Q1(3,i)*(105.0_r8_kind*r12(i)*r12(1)*Qrr2)/dr**9- &
                                  Q1(1,i)*(105.0_r8_kind*r12(i)*r12(3)*Qrr2)/dr**9+ &
                                  Q1(3,i)*six*Q2(i,1)/dr**5- &
                                  Q1(1,i)*six*Q2(i,3)/dr**5

                      T1(3)=T1(3)+Q1(1,i)*(-30.0_r8_kind*(r12(2)*Qr2(i)+r12(i)*Qr2(2)))/dr**7- &
                                  Q1(2,i)*(-30.0_r8_kind*(r12(1)*Qr2(i)+r12(i)*Qr2(1)))/dr**7+ &
                                  Q1(1,i)*(105.0_r8_kind*r12(i)*r12(2)*Qrr2)/dr**9- &
                                  Q1(2,i)*(105.0_r8_kind*r12(i)*r12(1)*Qrr2)/dr**9+ &
                                  Q1(1,i)*six*Q2(i,2)/dr**5- &
                                  Q1(2,i)*six*Q2(i,1)/dr**5

                      T2(1)=T2(1)+Q2(2,i)*(-30.0_r8_kind*(r12(3)*Qr1(i)+r12(i)*Qr1(3)))/dr**7- &
                                  Q2(3,i)*(-30.0_r8_kind*(r12(2)*Qr1(i)+r12(i)*Qr1(2)))/dr**7+ &
                                  Q2(2,i)*(105.0_r8_kind*r12(i)*r12(3)*Qrr1)/dr**9- &
                                  Q2(3,i)*(105.0_r8_kind*r12(i)*r12(2)*Qrr1)/dr**9+ &
                                  Q2(2,i)*six*Q1(i,3)/dr**5- &
                                  Q2(3,i)*six*Q1(i,2)/dr**5

                      T2(2)=T2(2)+Q2(3,i)*(-30.0_r8_kind*(r12(1)*Qr1(i)+r12(i)*Qr1(1)))/dr**7- &
                                  Q2(1,i)*(-30.0_r8_kind*(r12(3)*Qr1(i)+r12(i)*Qr1(3)))/dr**7+ &
                                  Q2(3,i)*(105.0_r8_kind*r12(i)*r12(1)*Qrr1)/dr**9- &
                                  Q2(1,i)*(105.0_r8_kind*r12(i)*r12(3)*Qrr1)/dr**9+ &
                                  Q2(3,i)*six*Q1(i,1)/dr**5- &
                                  Q2(1,i)*six*Q1(i,3)/dr**5

                      T2(3)=T2(3)+Q2(1,i)*(-30.0_r8_kind*(r12(2)*Qr1(i)+r12(i)*Qr1(2)))/dr**7- &
                                  Q2(2,i)*(-30.0_r8_kind*(r12(1)*Qr1(i)+r12(i)*Qr1(1)))/dr**7+ &
                                  Q2(1,i)*(105.0_r8_kind*r12(i)*r12(2)*Qrr1)/dr**9- &
                                  Q2(2,i)*(105.0_r8_kind*r12(i)*r12(1)*Qrr1)/dr**9+ &
                                  Q2(1,i)*six*Q1(i,2)/dr**5- &
                                  Q2(2,i)*six*Q1(i,1)/dr**5

                   end if
                end do

                if(do_energy) efp_efp_pq_pq_energy=efp_efp_pq_pq_energy+(six*Eij/dr**5- &
                     60.0_r8_kind*Eijk/dr**7+105.0_r8_kind*Eijkl/dr**9)/nine

                if(do_grads) then
                   Gij=-30.0_r8_kind*Eij*r12/dr**7

                   Gijk=-420.0_r8_kind*Eijk*r12/dr**9+60.0_r8_kind*Gijk/dr**7

                   Gijkl=-945_r8_kind*Eijkl*r12/dr**11+105.0_r8_kind*Gijkl/dr**9

                   gradient1=gradient1+(Gij-Gijk+Gijkl)/nine

                   gradient2=gradient2-(Gij-Gijk+Gijkl)/nine

                   torque1=torque1+two*T1/nine

                   torque2=torque2+two*T2/nine

                   do i2=1,grad_dim2
                      gradient_quad_totalsym(index2+i2)=gradient_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*gradient2(:))
                      torque_quad_totalsym(index2+i2)=torque_quad_totalsym(index2+i2)+ &
                           sum(rotmat2(i2,:)*torque2(:))
                   enddo
                end if
             end do ce2
          end do cu2

          if(do_grads) then
             do i1=1,grad_dim1
                gradient_quad_totalsym(index1+i1)=gradient_quad_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*gradient1(:))
                torque_quad_totalsym(index1+i1)=torque_quad_totalsym(index1+i1)+ &
                     sum(rotmat1(i1,:)*torque1(:))
             enddo
          end if

       end do ce1
    end do cu1


    if(do_energy .and. print_energy) print*,'EFP-EFP QUADRUPOLE-QUADRUPOLE ENERGY:', efp_efp_pq_pq_energy

    efp_efp_en=efp_efp_en+efp_efp_pq_pq_energy

  end subroutine efp_efp_pq_pq
  !*************************************************************

  !*************************************************************
  subroutine efp_efp_polar(to_do)
    !  Purpose: Calculate energy and gardients of EFP-EFP interactions
    !------------ Modules used ------------------- ---------------
    use efp_polar_module, only: calc_E_mp,calc_E_id,allocate_Efield,deallocate_Efield,calc_id1
    use efp_polar_module, only: calc_id,calc_id_FFenergy,calc_id_FFgrads,total_id,calc_polar
    use efp_solv_module, only: efp_id_mp_solv_energy,allocate_V_and_Q_id
    use efp_solv_module, only: allocate_E_cav,deallocate_E_cav,deallocate_V_id
    use induced_dipoles_module, only: free_field_arrays,do_pol_pcm
    use operations_module, only: operations_solvation_effect
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: to_do
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: efp_efp_pol_energy,efp_efp_pol_energy_old,d_En,efp_pol_solv_energy
    logical :: do_energy, do_grads
    integer(i4_kind) :: id_cycle
    real(r8_kind), parameter :: small_en=1.0e-14_r8_kind
    !------------ Executable code ------------------------------------

    do_energy=to_do==energy
    do_grads=to_do==gradient

    if(do_energy) then 
       call allocate_Efield()
       call calc_E_mp()
       if(operations_solvation_effect .and. do_pol_pcm) then
          call allocate_V_and_Q_id()
          call allocate_E_cav()
       end if

       efp_efp_pol_energy_old=zero
       id_cycle=0
       if(print_id_scc) then 
          print*,'===================================='
          print*,'Calculation of induced dipoles'
          print*,'------------------------------------'
          print*,'cycle      energy          D_E/E'
       end if

       calc_idip: do
          id_cycle=id_cycle+1
!          if(id_cycle > 1) call calc_E_id()
!          call calc_id(id_cycle)
          call calc_id1()
          call calc_id_FFenergy(efp_efp_pol_energy)
          efp_pol_solv_energy=zero
          if(operations_solvation_effect) then 
             call efp_id_mp_solv_energy(efp_pol_solv_energy)
             efp_efp_pol_energy=efp_efp_pol_energy+efp_pol_solv_energy
          end if

          d_En=abs((efp_efp_pol_energy-efp_efp_pol_energy_old)/efp_efp_pol_energy)

          if(print_id_scc) print('(i2,3x,f18.15,3x,e11.4)'),id_cycle, &
               efp_efp_pol_energy,d_En

          if(d_En <= small1) exit calc_idip

          efp_efp_pol_energy_old=efp_efp_pol_energy
       end do calc_idip
       if(print_id_scc) then 
          print*,'===================================='
       end if

!       print*,'Total induced dipole:', total_id()
       if(print_energy) print*,'EFP_EFP POLARIZABLE ENERGY:', efp_efp_pol_energy-efp_pol_solv_energy
       efp_efp_en=efp_efp_en+efp_efp_pol_energy
       if(print_energy .and. operations_solvation_effect) then 
          print*,'SOLVATION ENERGY DUE TO EFP INDUCED DIPOLES AND MULTIPOLES -' 
          print*,'0.5*Qid*Vmp:', efp_pol_solv_energy
       end if

       if(.not.calc_polar) then
          call deallocate_Efield()
          call free_field_arrays()
          if(operations_solvation_effect .and. do_pol_pcm) then
             call deallocate_V_id()
             call deallocate_E_cav()
          end if
       end if
    else if(do_grads) then
       call calc_id_FFgrads() 
    end if

  end subroutine efp_efp_polar
  !*************************************************************

  !--------------- End of module -------------------------------------
end module efp_efp_module

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
module efp_solv_grad_module
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
#include <def.h>

  use type_module ! type specification parameters
  use datatype
  use common_data_module

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------
  !------------ public functions and subroutines ------------------
  public init_X_solv_grads, X_solv_grads_shutdown
  public totsym_X_solv_grad_pack, totsym_X_solv_grad_unpack
  public X_solv_grad_cart_write, transform_X_solv_grad_to_cart
  public dV_dr, dV_dtet
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------
  !------------ Declaration of constants and variables ------------
  type(arrmat2),public,allocatable :: gradient_mpole_solv_cartesian(:)
  type(arrmat2),public,allocatable :: torque_mpole_solv_cartesian(:)
  real(r8_kind),public,allocatable,target :: gradient_mpole_solv_totalsym(:)
  real(r8_kind),public,allocatable,target :: torque_mpole_solv_totalsym(:)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains
  !*************************************************************
  subroutine init_X_solv_grads()
    ! Purpose: allocate arrays to calculate solv gradients on X centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use pointcharge_module, only : pointcharge_N,pointcharge_array,totsym_grad_pc_length
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i,n_eq
    !------------ Executable code --------------------------------

    allocate(gradient_mpole_solv_cartesian(pointcharge_N), stat=status)
    ASSERT(status==0)
    allocate(torque_mpole_solv_cartesian(pointcharge_N), stat=status)
    ASSERT(status==0)

    do i=1,pointcharge_N
       n_eq=pointcharge_array(i)%N_equal_charges
       allocate(gradient_mpole_solv_cartesian(i)%m(3,n_eq),stat=status)
       ASSERT(status==0)
       gradient_mpole_solv_cartesian(i)%m=0.0_r8_kind
       allocate(torque_mpole_solv_cartesian(i)%m(3,n_eq),stat=status)
       ASSERT(status==0)
       torque_mpole_solv_cartesian(i)%m=0.0_r8_kind
    end do

    allocate(gradient_mpole_solv_totalsym(totsym_grad_pc_length),stat=status)
    ASSERT(status==0)
    gradient_mpole_solv_totalsym=0.0_r8_kind
    allocate(torque_mpole_solv_totalsym(totsym_grad_pc_length),stat=status)
    ASSERT(status==0)
    torque_mpole_solv_totalsym=0.0_r8_kind

  end subroutine init_X_solv_grads
  !*************************************************************

  !*************************************************************
  subroutine X_solv_grads_shutdown()
    ! Purpose: free arrays for calculating solv gradients on X centers
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use pointcharge_module, only : pointcharge_N
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status,i
    !------------ Executable code --------------------------------

    do i=1,pointcharge_N
       deallocate(gradient_mpole_solv_cartesian(i)%m,stat=status)
       ASSERT(status==0)
       deallocate(torque_mpole_solv_cartesian(i)%m,stat=status)
       ASSERT(status==0)
    end do

    deallocate(gradient_mpole_solv_cartesian, stat=status)
    ASSERT(status==0)
    deallocate(torque_mpole_solv_cartesian, stat=status)
    ASSERT(status==0)

    deallocate(gradient_mpole_solv_totalsym,stat=status)
    ASSERT(status==0)
    deallocate(torque_mpole_solv_totalsym,stat=status)
    ASSERT(status==0)

  end subroutine X_solv_grads_shutdown
  !*************************************************************

  !*************************************************************
  subroutine totsym_X_solv_grad_pack()
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use pointcharge_module, only : totsym_grad_pc_length
    use comm_module, only: commpack
    !------------ Declaration of local variables ----------------
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    call commpack(gradient_mpole_solv_totalsym,totsym_grad_pc_length,1,info)
    ASSERT(info==0)
    call commpack(torque_mpole_solv_totalsym,totsym_grad_pc_length,1,info)
    ASSERT(info==0)

  end subroutine totsym_X_solv_grad_pack
  !*************************************************************

  !*************************************************************
  subroutine totsym_X_solv_grad_unpack()
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use pointcharge_module, only : totsym_grad_pc_length
    use comm_module, only: communpack
    !------------ Declaration of local variables ----------------
    real(r8_kind) :: help_arr(totsym_grad_pc_length)
    integer(i4_kind) :: info
    !------------ Executable code -------------------------------

    call communpack(help_arr,totsym_grad_pc_length,1,info)
    ASSERT(info==0)
    gradient_mpole_solv_totalsym=gradient_mpole_solv_totalsym+help_arr
    call communpack(help_arr,totsym_grad_pc_length,1,info)
    ASSERT(info==0)
    torque_mpole_solv_totalsym=torque_mpole_solv_totalsym+help_arr

  end subroutine totsym_X_solv_grad_unpack
  !*************************************************************

  !*************************************************************
  subroutine transform_X_solv_grad_to_cart()
    ! purpose: transform symmetry adapted gradient components to cartesian
    !          coordinates and add them to an array of cartesian gradients
    !------------ modules used -----------------------------------
    use pointcharge_module, only : pointcharge_array, pointcharge_N, unique_pc_grad_info, &
         unique_pc_index
   !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: i_unique,i_equal,index,i,grad_dim
    real(kind=r8_kind),pointer :: rotmat(:,:)
    !------------ Executable code -------------------------------

    do i_unique=1,pointcharge_N
       index=unique_pc_index(i_unique)
       grad_dim=unique_pc_index(i_unique+1)-index
       index=index-1
       do i_equal=1,pointcharge_array(i_unique)%n_equal_charges
          rotmat=>unique_pc_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             gradient_mpole_solv_cartesian(i_unique)%m(:,i_equal) = &
                  gradient_mpole_solv_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*gradient_mpole_solv_totalsym(index+i)

             torque_mpole_solv_cartesian(i_unique)%m(:,i_equal) = &
                  torque_mpole_solv_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*torque_mpole_solv_totalsym(index+i)
          end do
       end do
    enddo

  end subroutine transform_X_solv_grad_to_cart
  !*************************************************************  

  !*************************************************************
  subroutine X_solv_grad_cart_write()
    !  Purpose: writing the cartesian gradients
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
    use iounitadmin_module, only: output_unit
    use pointcharge_module, only : pointcharge_array, pointcharge_N
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,j
    !------------ Executable code --------------------------------

    if(pointcharge_N > 0) &
         write(output_unit,'(/A)') 'Solvation Cartesian gradients and torques on EFP Centers'
    do i=1,pointcharge_N
       write(output_unit,*) 'Unique Center:',i
       do j=1,pointcharge_array(i)%N_equal_charges
          write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
               gradient_mpole_solv_cartesian(i)%m(:,j),' / ',torque_mpole_solv_cartesian(i)%m(:,j)
       end do
    end do

  end subroutine X_solv_grad_cart_write
  !*************************************************************

  !*************************************************************
  function dV_dr(i_uc,r,dr,QQ_s,Z,D1,Q,D2,Q1_s,Q2_s) result(grad)
    ! Purpose :
    !------------ Modules used -----------------------------------
    use unique_atom_module, only: N_unique_atoms
    use pointcharge_module, only: pointcharge_N
    use point_dqo_module, only: N_pd, N_pq
    use induced_dipoles_module, only: N_ipd
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: grad(3)
    integer(i4_kind), intent(in) :: i_uc
    real(r8_kind), intent(in) :: r(3),dr,QQ_s,Z,D1(3),Q(3,3),D2(3),Q1_s,Q2_s
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: E1(3), E2(3)
    integer(i4_kind) :: i,j,k
    !------------ Executable code --------------------------------

    grad=zero
    if    (i_uc <= N_unique_atoms) then
       grad=-(QQ_s*Z/dr**3)*r

    elseif(i_uc >  N_unique_atoms .and. i_uc <= N_unique_atoms+pointcharge_N) then
       grad=grad-(QQ_s*Z/dr**3)*r

       if(N_pd > 0) grad=grad+QQ_s*(three*dot_product(D1,r)*r/dr**5-D1/dr**3)

       if(N_pq > 0) then
          E1=zero
          do i=1,3
             do j=1,3
                E1=E1-(fifteen*r(i)*r(j)*r/dr**7)*Q(i,j)
                do k=1,3
                   if(k==i) E1(k)=E1(k)+(three*r(j)/dr**5)*Q(i,j)
                   if(k==j) E1(k)=E1(k)+(three*r(i)/dr**5)*Q(i,j)
                end do
                if(i==j) E1=E1+(three*r/dr**5)*Q(i,j)
             end do
          end do
          grad=grad+QQ_s*E1/three       
       end if

    elseif(i_uc >  N_unique_atoms+pointcharge_N) then
       
       E1=three*dot_product(D1,r)*r/dr**5-D1/dr**3
       E2=three*dot_product(D2,r)*r/dr**5-D2/dr**3
       grad=(E1*Q2_s+E2*Q1_s)*half
    end if

  end function dV_dr
  !*************************************************************

  !*************************************************************
  function dV_dtet(i_uc,r,dr,QQ_s,Z,D1,Q,D2,Q1_s,Q2_s) result(torque)
    ! Purpose :
    !------------ Modules used -----------------------------------
    use unique_atom_module, only: N_unique_atoms
    use pointcharge_module, only: pointcharge_N
    use point_dqo_module, only: N_pd, N_pq
    use induced_dipoles_module, only: N_ipd
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: torque(3)
    integer(i4_kind), intent(in) :: i_uc
    real(r8_kind), intent(in) :: r(3),dr,QQ_s,Z,D1(3),Q(3,3),D2(3),Q1_s,Q2_s
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: Qr(3),E1(3),E2(3)
    !------------ Executable code --------------------------------

    torque=zero
    if    (i_uc <= N_unique_atoms) then

    elseif(i_uc >  N_unique_atoms .and. i_uc <= N_unique_atoms+pointcharge_N) then

       if(N_pd > 0) torque=torque-(QQ_s/dr**3)*vector_product(D1,r)

       if(N_pq > 0) then
          Qr=matmul(Q,r)
          E1=two*QQ_s*r/dr**5
          torque=torque+vector_product(Qr,E1)
       end if

    elseif(i_uc >  N_unique_atoms+pointcharge_N) then
       
       E1=(Q2_s/dr**3)*vector_product(D1,r)
       E2=(Q1_s/dr**3)*vector_product(D2,r)
       torque=-(E1+E2)*half
    end if

  end function dV_dtet
  !*************************************************************

  !*************************************************************
  function vector_product(v1,v2)
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: vector_product(3)
    real(r8_kind) :: v1(3),v2(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------
     
    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !*************************************************************

  !*************************************************************
end module efp_solv_grad_module

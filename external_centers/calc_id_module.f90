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
module calc_id_module
  !---------------------------------------------------------------
  !
  !  Purpose: calculate induced point dipoles
  !           (accounting for polarization of eXternal centers)
  !
  !
  !  Module called by:
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AS
  !  Date: 12.02.2008
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
  use unique_atom_module, only: unique_atom_type
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use pointcharge_module, only: rcm
#endif
  use induced_dipoles_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements -----------------------------

  !------------ public functions and subroutines ------------------
  public calc_nuc_id_grads, calc_id_energy, transform_id_grad_to_cart
  public id_grad_cart_write, calc_id_grads, calc_induced_dipmom
  public build_Pol_ham, calc_Pol_ham, dealloc_Pol_ham
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !----------------------------------------------------------------

  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
#ifdef WITH_EFP
  subroutine calc_induced_dipmom(E_tot,E_tot1,E_id)
#else
  subroutine calc_induced_dipmom(E_tot,E_tot1)
#endif
    ! Purpose: Calculate induced dipole moment of external polarizable
    ! centers
    !------------ Declaration of formal parameters ---------------
    type(arrmat2), intent(in) :: E_tot(:),E_tot1(:)   ! cartesian field array
#ifdef WITH_EFP
    type(arrmat2), intent(in) :: E_id(:)   ! cartesian field array
#endif
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i,j,n_eq
    real(r8_kind) :: alpha(3,3),E(3)
    !--- executable code-------------------------------------

    ASSERT(size(E_tot)==N_ipd)
    E=0.0_r8_kind

    do i=1,N_ipd
       n_eq=ipd_array(i)%N_equal_dipoles
       do j=1,n_eq
          alpha=ipd_array(i)%pol_tensor(:,:,j)
          E=E_tot(i)%m(:,j)
          ipd_array(i)%idipole(:,j)=matmul(alpha,E)
          E=E_tot1(i)%m(:,j)
          ipd_array(i)%idipole1(:,j)=matmul(transpose(alpha),E)
       end do
    end do

#ifdef WITH_EFP
    call cart_idip_to_totsym(E_id)
#else
    call cart_idip_to_totsym()
#endif

    call send_receive_id()

  end subroutine calc_induced_dipmom
  !*************************************************************

  !*************************************************************
#ifdef WITH_EFP
  subroutine calc_ind_dip_iter(E_tot,E_tot1,print_iter)
    ! Purpose: Calculate induced dipole moment of external polarizable
    ! centers
    use efp_polar_module, only: calc_id_FFenergy,calc_id1,calc_polar,calc_E_id
    use elec_static_field_module, only : E_id
    !------------ Declaration of formal parameters ---------------
    type(arrmat2), intent(in) :: E_tot(:),E_tot1(:)   ! cartesian field array
    logical, intent(in) :: print_iter
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    real(r8_kind) :: energy,energy_old,d_En
    integer(i4_kind) :: id_cycle
    real(r8_kind), parameter :: small_en=1.0e-16_r8_kind
    !--- executable code-------------------------------------

    calc_polar=.true.
    energy_old=0.0_r8_kind
    id_cycle=0
    if(print_iter) then 
       print*,'===================================='
       print*,'Calculation of induced dipoles'
       print*,'------------------------------------'
       print*,'cycle      energy          D_E/E'
    end if

    calc_idip: do
       id_cycle=id_cycle+1
       call calc_id1()
       call calc_id_FFenergy(energy)

       d_En=abs((energy-energy_old)/energy)

       if(print_iter) print('(i2,3x,f18.15,3x,e11.4)'),id_cycle,energy,d_En

       if(d_En <= small_en) exit calc_idip

       energy_old=energy
    end do calc_idip
    if(print_iter) then 
       print*,'===================================='
    end if
    calc_polar=.false.

    call calc_E_id()
    call cart_idip_to_totsym(E_id)

  end subroutine calc_ind_dip_iter
#endif
  !*************************************************************

  !*************************************************************
#ifdef WITH_EFP
  subroutine cart_idip_to_totsym(E_id)
    use efp_module, only: n_efp
    type(arrmat2), intent(in) :: E_id(:)   ! cartesian field array
#else
  subroutine cart_idip_to_totsym()
#endif
    ! Purpose : Convert cartesian induced dipoles to total symmetric ones
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ux,i_ex,n_ex,i
    real(r8_kind), pointer :: rotmat(:,:)
    integer(i4_kind) :: grad_dim,index
    real(r8_kind),pointer  :: idipole(:,:),idipole1(:,:)
    real(r8_kind) :: dip_sum(3),DD
    !--- executable code-------------------------------------

    totsym_ind_dip=0.0_r8_kind

    do i_ux=1,N_ipd
       n_ex=ipd_array(i_ux)%N_equal_dipoles
       idipole=>ipd_array(i_ux)%idipole
       idipole1=>ipd_array(i_ux)%idipole1
       
       grad_dim=ind_dip_grad_index(i_ux+1)-ind_dip_grad_index(i_ux)
       index=ind_dip_grad_index(i_ux)-1
        
       do i_ex=1,n_ex
          DD=sqrt(dot_product(idipole(:,i_ex),idipole(:,i_ex)))
          if(DD >= 0.2_r8_kind) then
             print*,'WARNING: Induced dipole too large (> 0.2 au): ',DD,' Center: ', i_ux 
          end if
          dip_sum=idipole(:,i_ex)+idipole1(:,i_ex)
#if 0
#ifdef WITH_EFP
! the six lines below look as wrong implentation. But may be not. Future will show
          if(n_efp > 1) then
             E=E_id(i_ux)%m(:,i_ex)
             alpha=ipd_array(i_ux)%pol_tensor(:,:,i_ex)
             idip_id1=matmul(transpose(alpha),E)
             dip_sum=dip_sum-idip_id1
          end if
#endif
#endif

          rotmat=>ind_dip_grad_info(i_ux)%m(:,:,i_ex)

          do i=1,grad_dim
             totsym_ind_dip(index+i)=totsym_ind_dip(index+i)+ &
                  sum(rotmat(i,:)*dip_sum)
          enddo

       end do
    end do

  end subroutine cart_idip_to_totsym
  !*************************************************************

  !*************************************************************
#ifdef WITH_EFP
  subroutine build_Pol_ham(print_micro_iter,update_id)
#else
  subroutine build_Pol_ham(update_id)
#endif
    ! Purpose: Initiate calculation of additional contribution into the total
    ! hamiltonian of the system produced by polarizable external
    ! centers
    use elec_static_field_module, only : start_read_field_e, E_ele, E_nuc
#ifdef WITH_EFP
    use elec_static_field_module, only : E_mp, E_id, E_id1, E_cav, E_cav1
    use efp_module, only: n_efp
    use efp_polar_module, only: calc_E_id
    use operations_module, only: operations_solvation_effect
    use efp_solv_module, only: calc_E_cav, do_solv
#endif
    use casc_logic_module
    use comm_module, only : comm_get_n_processors, comm_init_send, comm_send, &
         comm_parallel, commpack, communpack
    use msgtag_module, only: msgtag_Pol_ham
    !------------ Declaration of formal parameters ---------------
#ifdef WITH_EFP
    logical :: print_micro_iter
#endif
    logical :: update_id
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    type(arrmat2), allocatable :: E_tot(:),E_tot1(:)
    integer(i4_kind) :: i,status,m1,m2,i_pr
    integer(i4_kind), allocatable :: vpn(:)
    !--- executable code-------------------------------------

    if(.not. update_id) goto 100
    call start_read_field_e()
#ifdef WITH_EFP
    do_iter:if(do_iterations .and. n_efp > 1) then

       call calc_ind_dip_iter(E_tot,E_tot1,print_micro_iter)

    else !do_iter
#endif
#ifdef WITH_EFP
       if(n_efp > 1) call calc_E_id()
       if(operations_solvation_effect .and. n_efp > 0 .and. do_solv .and. do_pol_pcm) then
          call calc_E_cav()
       end if
#endif

       allocate(E_tot(N_ipd),E_tot1(N_ipd),stat=status)
       ASSERT(status==0)

       do i=1,N_ipd
          m1=size(E_ele(i)%m,1)
          m2=size(E_ele(i)%m,2)
          allocate(E_tot(i)%m(m1,m2),E_tot1(i)%m(m1,m2),stat=status)
          ASSERT(status==0)
          E_tot(i)%m=E_ele(i)%m+E_nuc(i)%m
          E_tot1(i)%m=E_ele(i)%m+E_nuc(i)%m
#ifdef WITH_EFP
          if(n_efp > 1) then 
             E_tot(i)%m=E_tot(i)%m+E_mp(i)%m+E_id(i)%m
             E_tot1(i)%m=E_tot1(i)%m+E_mp(i)%m+E_id1(i)%m
          end if
          if(operations_solvation_effect .and. n_efp > 0 .and. do_solv .and. do_pol_pcm) then
             E_tot(i)%m=E_tot(i)%m+E_cav(i)%m
             E_tot1(i)%m=E_tot1(i)%m+E_cav1(i)%m
          end if
#endif
       end do

#ifdef WITH_EFP
       call calc_induced_dipmom(E_tot,E_tot1,E_id)
#else
       call calc_induced_dipmom(E_tot,E_tot1)
#endif

       do i=1,N_ipd
          deallocate(E_tot(i)%m,E_tot1(i)%m,stat=status)
          ASSERT(status==0)
       end do
       deallocate(E_tot,E_tot1,stat=status)
       ASSERT(status==0)
#ifdef WITH_EFP
    end if do_iter
#endif

100 continue

    i_pr=comm_get_n_processors()

    allocate(vpn(i_pr),stat=status)
    ASSERT(status==0)

    call comp_vpn(vpn)

    if(comm_parallel()) then
       do i=2,i_pr
          call comm_init_send(i,msgtag_Pol_ham)
          call commpack(vpn(i),status)
          ASSERT(status==0)
          call comm_send()
       enddo
    endif

    call calc_Pol_ham(vpn(1))

    deallocate(vpn,stat=status)
    ASSERT(status==0)

  end subroutine build_Pol_ham
  !*************************************************************

  !*************************************************************
  subroutine alloc_Pol_ham
    use symmetry_data_module, only : symmetry_data_n_irreps, symmetry_data_dimension
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind)       :: n_irrep,status
    integer(i4_kind)       :: n,i
    !--- executable code-------------------------------------

    n_irrep = symmetry_data_n_irreps()

    allocate (ham_Pol(n_irrep),STAT=status)
    ASSERT(status==0)

    do i=1,n_irrep
       n=symmetry_data_dimension(i)
       allocate(ham_Pol(i)%m(n,n),STAT=status)
       ASSERT(status==0)
       ham_Pol(i)%m=0.0_r8_kind
    enddo

  end subroutine alloc_Pol_ham
  !*************************************************************

  !*************************************************************
  subroutine dealloc_Pol_ham
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: status,i
    !--- executable code-------------------------------------

    do i=1,size(ham_Pol)
       deallocate(ham_Pol(i)%m,stat=status)
       ASSERT(status==0)
    end do

    deallocate(ham_Pol,stat=status)
    ASSERT(status==0)

  end subroutine dealloc_Pol_ham
  !*************************************************************

  !*************************************************************
  subroutine calc_Pol_ham(vpnx)
    use options_module, only : options_integrals_on_file
    use elec_static_field_module,  only : bounds, field_integral_open, field_integral_close
    use readwriteblocked_module
    use msgtag_module, only : msgtag_back_Pol_ham
    use symmetry_data_module, only : symmetry_data_n_irreps, symmetry_data_dimension
    use comm_module, only : comm_myindex, comm_get_n_processors, comm_i_am_master, &
    comm_save_recv, comm_init_send, comm_send, commpack, communpack
    use casc_logic_module
    use integralstore_module, only : integralstore_3c_field
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),optional :: vpnx
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    type(readwriteblocked_tapehandle) :: th_field
    real(r8_kind),pointer :: field_int(:)
    real(r8_kind), allocatable :: sum_x(:)
    real(r8_kind), allocatable :: help_arr(:,:,:)
    integer(i4_kind), allocatable :: dim_irrep(:)
    integer(i4_kind), allocatable :: rpn(:)
    integer(i4_kind)       :: item_arr_field,my_ind,n_irrep,i_gamma, &
         first_i,vpn1,vpn2,i_pr,di,nitem1, &
         lower_node, upper_node
    logical ::  integrals_on_file
    integer(i4_kind) :: n,status,info,i_mn,m,n_mn,k,i,i_meta,i_last
    !--- executable code-------------------------------------

    integrals_on_file=options_integrals_on_file()

    n_irrep = symmetry_data_n_irreps()

    allocate(dim_irrep(n_irrep),stat=status)
    ASSERT(status==0)

    do n=1,n_irrep
       dim_irrep(n) = symmetry_data_dimension(n)
    enddo

    my_ind=comm_myindex()
    item_arr_field=bounds%item_arr(my_ind)
    first_i=1
    do i=1,my_ind-1
       first_i=first_i+bounds%item_arr(i)
    enddo

    if(comm_i_am_master()) then
       ASSERT(present(vpnx))
       vpn1=vpnx
    else
       call communpack(vpn1,1,1,info)
       ASSERT(info==0)
    endif

    i_pr=comm_get_n_processors()
    allocate(rpn(2*i_pr-1),stat=status)
    ASSERT(status==0)
    call comp_rpn(rpn)

    call alloc_Pol_ham()

    if ( integrals_on_file ) then
       ! open the integral files ----------------------
       if (item_arr_field /= 0) call field_integral_open(th_field)
    else
       i_meta=1
    endif

    if (integrals_on_file) then 
       allocate( field_int(item_arr_field),STAT=status)
       ASSERT(status==0)
    endif

    i_gamma_lab: do i_gamma = 1,n_irrep
       vpn2=vpn1
       n_mn = ( dim_irrep(i_gamma) * (dim_irrep(i_gamma) + 1) ) / 2

       if (item_arr_field /= 0) then
          allocate( sum_x(n_mn), STAT=status )
          ASSERT(status==0)
          sum_x = 0.0_r8_kind

          i_mn = 1
          do m=1,dim_irrep(i_gamma)
             do n=1,m
                ! read in integrals file
                if (integrals_on_file) then
                   call readwriteblocked_read(field_int,th_field)
                else
                   i_last=i_meta+item_arr_field-1
                   field_int => integralstore_3c_field(i_meta:i_last)
                   i_meta=i_last+1
                endif
                do k=1,item_arr_field
                   sum_x(i_mn) =sum_x(i_mn)+field_int(k)*totsym_ind_dip(k+first_i-1)*0.5_r8_kind
                enddo
                i_mn = i_mn + 1
             enddo
          enddo

          i_mn = 1
          do m=1,dim_irrep(i_gamma)
             do n=1,m-1
                ham_Pol(i_gamma)%m(m,n) = sum_x(i_mn)
                ham_Pol(i_gamma)%m(n,m) = sum_x(i_mn)
                i_mn = i_mn + 1
             enddo
             ham_Pol(i_gamma)%m(m,m) = sum_x(i_mn)
             i_mn = i_mn + 1
          enddo

          deallocate( sum_x, STAT=status )
          ASSERT(status==0)
       end if

       !cascadic procedure sending Hamiltonian data to master
       if(2*(vpn1/2) == vpn1) then
          allocate(help_arr(dim_irrep(i_gamma),dim_irrep(i_gamma),1),stat=status)
          ASSERT(status==0)
       endif

       di=dim_irrep(i_gamma)
       nitem1=di*di

       vpn2_lab: do while(vpn2 > 1)
          if(mod(vpn2,2) == 0) then
             vpn2=vpn2/2
             lower_node=rpn(2*vpn2+1)

             call comm_save_recv(lower_node,msgtag_back_Pol_ham)

             call communpack(help_arr(1,1,1),nitem1,1,info)
             ASSERT(info==0)

             ham_Pol(i_gamma)%m=ham_Pol(i_gamma)%m+help_arr(:,:,1)
          else
             upper_node=rpn((vpn2-1)/2)

             call comm_init_send(upper_node,msgtag_back_Pol_ham)

             call commpack(ham_Pol(i_gamma)%m(1,1),nitem1,1,info)
             ASSERT(info==0)

             call comm_send()
             vpn2=0
          endif
       enddo vpn2_lab
       if(mod(vpn1,2)==0) then
          deallocate(help_arr,stat=status)
          ASSERT(status==0)
       endif
    end do i_gamma_lab

    if (integrals_on_file) then
       deallocate(field_int,STAT=status)
       ASSERT(status==0)
    endif
    if ( integrals_on_file ) then
       ! close the integral files ----------------------
       if (item_arr_field /= 0) call field_integral_close(th_field)
    endif

    if(.not.comm_i_am_master()) call dealloc_Pol_ham()

    deallocate(dim_irrep,stat=status)
    ASSERT(status==0)
    deallocate(rpn,stat=status)
    ASSERT(status==0)

  end subroutine calc_Pol_ham
  !*************************************************************

  !*************************************************************
  subroutine calc_id_energy()
    ! Purpose: Calculate interaction energy between QM atoms
    ! and induced dipoles
    use elec_static_field_module, only : E_nuc,E_ele
#ifdef WITH_EFP
    use elec_static_field_module, only : E_mp
    use efp_module, only: n_efp
#endif
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
             E_f=E_nuc(i_uc)%m(:,i_ec)+E_ele(i_uc)%m(:,i_ec)
#ifdef WITH_EFP
             if(n_efp > 1) E_f=E_f+E_mp(i_uc)%m(:,i_ec)
#endif
             en_ind_dipole=en_ind_dipole- &
                  dot_product(E_f,idipole(:,i_ec))*0.5_r8_kind
          end do
       end do
!print*,en_ind_dipole,'en_ind_dipole'
    end if

  end subroutine calc_id_energy
  !*************************************************************

  !*************************************************************
  subroutine calc_nuc_id_grads(totsym_grad,gradient_index)
    ! Purpose: Calculate nuclear gradients of interaction energy between atomic nuclears
    ! end external induced dipoles
    use unique_atom_module, only : unique_atoms, &
         pseudopot_present,N_moving_unique_atoms,moving_unique_atom_index, &
         unique_atom_grad_info
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), intent(inout) :: totsym_grad(:)
    integer(i4_kind), intent(in) :: gradient_index(N_moving_unique_atoms + 1)
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ua,i_ea,i_uc,i_ec,ma,i
    integer(i4_kind) :: n_ec,n_ea
    real(r8_kind),pointer  :: xa(:,:),xc(:,:),rotmat(:,:)
    real(r8_kind),pointer  :: idipole(:,:),idipole1(:,:)
    real(r8_kind) :: gradient(3)
    real(r8_kind) :: Z, Zc,r_ac(3),d_ac
    integer(i4_kind) :: grad_dim,index
    !--- executable code-------------------------------------

    uniq_at: do ma=1,N_moving_unique_atoms
       i_ua=moving_unique_atom_index(ma)
       n_ea=unique_atoms(i_ua)%n_equal_atoms
       Z=unique_atoms(i_ua)%Z
       Zc=unique_atoms(i_ua)%Zc
       if (.not.pseudopot_present) Zc = 0.0_r8_kind
       Z=Z-Zc
       xa=>unique_atoms(i_ua)%position

       grad_dim=gradient_index(ma+1)-gradient_index(ma)
       index=gradient_index(ma)-1
       
       eq_at: do i_ea=1,n_ea
          rotmat=>unique_atom_grad_info(ma)%m(:,:,i_ea)
          gradient=0.0_r8_kind

          uniq_pd: do i_uc=1,N_ipd
             n_ec=ipd_array(i_uc)%N_equal_dipoles
             xc=>ipd_array(i_uc)%position
             idipole=>ipd_array(i_uc)%idipole
             idipole1=>ipd_array(i_uc)%idipole1

             eq_pd: do i_ec=1,n_ec
                r_ac=xa(:,i_ea)-xc(:,i_ec)
                d_ac=sqrt(dot_product(r_ac,r_ac))
                if(d_ac < 0.001_r8_kind) d_ac=0.001_r8_kind

                gradient=gradient-(Z*idipole(:,i_ec)/d_ac**3- &
                     3.0_r8_kind*Z*dot_product(r_ac,idipole(:,i_ec))*r_ac/d_ac**5)- &
                     (Z*idipole1(:,i_ec)/d_ac**3- &
                     3.0_r8_kind*Z*dot_product(r_ac,idipole1(:,i_ec))*r_ac/d_ac**5)
             end do eq_pd
          end do uniq_pd

!print*,ma,i_ea,gradient
          do i=1,grad_dim
             totsym_grad(index+i) = totsym_grad(index+i)+ &
                  sum(rotmat(i,:)*gradient(:))*0.5_r8_kind
             
          enddo
       end do eq_at
    end do uniq_at

  end subroutine calc_nuc_id_grads
  !*************************************************************


  !*************************************************************
  subroutine calc_id_grads()
    ! Purpose: Calculate gradients on Pol. centers (interaction with atomic nuclears)
    use unique_atom_module, only : unique_atoms, N_unique_atoms, &
         pseudopot_present
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !----------- declaration of local variables -------------
    integer(i4_kind) :: i_ux,i_ex,n_ex,i_ua,i_ea,n_ea,i
    real(r8_kind), pointer :: xx(:,:),xa(:,:),rotmat(:,:)
    integer(i4_kind) :: grad_dim,index
    real(r8_kind) :: gradient(3),torque(3),r_ax(3),d_ax
    real(r8_kind),pointer  :: idipole(:,:),idipole1(:,:)
    real(r8_kind) :: Z_a
    !---executable code--------------------------------------

    ! Point dipoles
    do i_ux=1,N_ipd
       n_ex=ipd_array(i_ux)%N_equal_dipoles
       xx=>ipd_array(i_ux)%position
       idipole=>ipd_array(i_ux)%idipole
       idipole1=>ipd_array(i_ux)%idipole1
       
       grad_dim=ind_dip_grad_index(i_ux+1)-ind_dip_grad_index(i_ux)
       index=ind_dip_grad_index(i_ux)-1
        
       do i_ex=1,n_ex
          rotmat=>ind_dip_grad_info(i_ux)%m(:,:,i_ex)
          gradient=0.0_r8_kind
          torque=0.0_r8_kind

          do i_ua=1,N_unique_atoms
             n_ea=unique_atoms(i_ua)%N_equal_atoms
             Z_a=unique_atoms(i_ua)%Z
             if(pseudopot_present) Z_a=Z_a-unique_atoms(i_ua)%Zc
             xa=>unique_atoms(i_ua)%position
              
             do i_ea=1,n_ea
                r_ax=xa(:,i_ea)-xx(:,i_ex)
                d_ax=sqrt(dot_product(r_ax,r_ax))
                if(d_ax < 0.001_r8_kind) d_ax=0.001_r8_kind

                gradient=gradient-(Z_a*idipole(:,i_ex)/d_ax**3- &
                        3.0_r8_kind*Z_a*dot_product(r_ax,idipole(:,i_ex))*r_ax/d_ax**5)- &
                        (Z_a*idipole1(:,i_ex)/d_ax**3- &
                        3.0_r8_kind*Z_a*dot_product(r_ax,idipole1(:,i_ex))*r_ax/d_ax**5)

#ifdef WITH_EFP
                torque=torque-vector_product(idipole(:,i_ex)+idipole1(:,i_ex),r_ax)*Z_a/d_ax**3
#endif
             end do
          end do

          do i=1,grad_dim
             grad_idip_totalsym(index+i)=grad_idip_totalsym(index+i)- &
                  sum(rotmat(i,:)*gradient(:))*0.5_r8_kind

#ifdef WITH_EFP
             torque_idip_totalsym(index+i)=torque_idip_totalsym(index+i)+ & !-
                  sum(rotmat(i,:)*torque(:))*0.5_r8_kind
#endif
          enddo

       end do
    end do

  end subroutine calc_id_grads
  !*************************************************************

  !*************************************************************
  subroutine transform_id_grad_to_cart()
    ! purpose: transform symmetry adapted gradient components to cartesian
    !          coordinates and add them to an array of cartesian gradients
    !
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: i_unique,i_equal,index,i,grad_dim
    real(kind=r8_kind),pointer :: rotmat(:,:)
#ifdef WITH_EFP
    integer(kind=i4_kind) :: i_group
    real(kind=r8_kind) :: r12(3),r1(3),r2(3),f(3)
#endif
    !------------ Executable code -------------------------------

    do i_unique=1,N_ipd
       index=ind_dip_grad_index(i_unique)
       grad_dim=ind_dip_grad_index(i_unique+1)-index
       index=index-1
       do i_equal=1,ipd_array(i_unique)%n_equal_dipoles
          rotmat=>ind_dip_grad_info(i_unique)%m(:,:,i_equal)
          do i=1,grad_dim
             grad_idip_cartesian(i_unique)%m(:,i_equal) = &
                  grad_idip_cartesian(i_unique)%m(:,i_equal) + &
                  rotmat(i,:)*grad_idip_totalsym(index+i)
#ifdef WITH_EFP
             if(efp) then
                torque_idip_cartesian(i_unique)%m(:,i_equal) = &
                     torque_idip_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*torque_idip_totalsym(index+i)
             end if
#endif
          end do
       end do
    enddo
#ifdef WITH_EFP
    if(efp) then
       do i_unique=1,N_ipd
          do i_equal=1,ipd_array(i_unique)%n_equal_dipoles
             r1=ipd_array(i_unique)%position(:,i_equal)
             i_group=ipd_array(i_unique)%group(i_equal)
             r2=rcm(:,i_group)
             r12=r2-r1
             f=-grad_idip_cartesian(i_unique)%m(:,i_equal)

             torque_idip_cartesian(i_unique)%m(:,i_equal) = &
                  torque_idip_cartesian(i_unique)%m(:,i_equal) + vector_product(r12,f)
          end do
       end do
    end if
#endif

  end subroutine transform_id_grad_to_cart
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
  subroutine id_grad_cart_write()
    !  Purpose: writing the cartesian gradients
    !** End of interface *****************************************
    !------------ modules used -----------------------------------
     use iounitadmin_module, only: output_unit
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,j
    !------------ Executable code --------------------------------

    if(N_ipd > 0) then
#ifdef WITH_EFP
       if(efp) then
          write(output_unit,'(/A)') 'Cartesian gradients and torques on Induced Point Dipoles'
          do i=1,N_ipd
             write(output_unit,*) 'Unique Induced Point Dipole:',i 
             do j=1,ipd_array(i)%n_equal_dipoles
                write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                     grad_idip_cartesian(i)%m(:,j),' / ',torque_idip_cartesian(i)%m(:,j)
             enddo
          end do
       else
#endif
          write(output_unit,'(/A)') 'Cartesian gradients on Induced Point Dipoles'
          do i=1,N_ipd
             write(output_unit,*) 'Unique Induced Point Dipole:',i 
             do j=1,ipd_array(i)%n_equal_dipoles
                write(output_unit,'(A14,3F15.10)') 'Equal Center:',&
                     grad_idip_cartesian(i)%m(:,j)
             enddo
          end do
#ifdef WITH_EFP
       end if
#endif
    end if

  end subroutine id_grad_cart_write
  !*************************************************************
  !--------------- End of module ----------------------------------
end module calc_id_module

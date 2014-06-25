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
!================================================================
! Public interface of module
!================================================================
module xcfit_hamiltonian
  !---------------------------------------------------------------
  !
  !  Purpose: This module creates the XC-part of the hamiltionian by
  !           the fitting technique
  !
  !  References : Old lcgto,P. Knappe thesis, O. Haeberlen thesis,
  !               Kobenka, thesis
  !
  !  Author: MS
  !  Date: 2/96
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification
  ! Author: UB
  ! Date:   6/97
  ! Description: Introduction of xcfit_allocate and xcfit_deallocate
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module  ! type specification parameters
  use orbitalstore_module
  use orbital_module
  use density_data_module
  use comm_module
  use msgtag_module
  use hamiltonian_module
  use occupied_levels_module
  use iounitadmin_module
  use mixing_module, only: mixing_xc
  use input_module
  use becke_perdew_module
! use density_calc_module ! use in actual subs only!
  use linsys_module, only: decompose_linsys, solve_linsys
  implicit none
  private
  save

  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public :: build_xcfit_main, xcfit_setup, build_xcfit, xcfit_close

!================================================================
! End of public interface of module
!================================================================

  real(kind=r8_kind),allocatable :: rho(:,:),&
         dfdrho(:,:),&  ! derivative of f with respect to rho
         rhs(:,:),lhs(:,:), &
         fxc(:),&        ! funktion f according to JPG
         eps(:)

  real(kind=r8_kind),pointer  :: wts(:,:)
  integer(kind=i4_kind) :: ispin,n_irrep,n_xc,n_rhs,n_lhs,lhs_length,&
       rhs_length,vec_length
  integer(kind=i4_kind),allocatable :: dims(:),partners(:)
  real(kind=r8_kind)  :: charge_int,exc_int  ! numerically integrated charge
  ! and xc-energy
  type(orbital_type),pointer :: orbs_ob(:)   ! values of the basis functions
  type(orbital_type) :: orbs_xc ! values of the fitting functions


contains

  subroutine build_xcfit_main(loop)
    ! purpose : wrapper for build_xcfit; it runs only on the master and sends
    !           the message "execute build_xcfit" to the slaves. Subsequently
    !           build_xc is called. build_xcfit_main is called in every scf-
    !           cycle by main_scf
    integer(kind=i4_kind),optional :: loop ! number of the actual scf-cycle
    !** End of interface *****************************************
    if ( comm_parallel() ) then
       call comm_init_send(comm_all_other_hosts,msgtag_build_xcfit)
       call comm_send()
    endif
    call build_xcfit(loop)
  end subroutine build_xcfit_main

  !***********************************************************

  subroutine xcfit_allocate()
    use machineparameters_module, only: machineparameters_veclen
    ! purpose: allocates non-permanent data

    call orbital_setup(machineparameters_veclen)
    call orbital_allocate(orbs_ob=orbs_ob,orb_xc=orbs_xc)
  end subroutine xcfit_allocate
  !*******************************************************

  !*******************************************************
  subroutine xcfit_deallocate()
    ! purpose: deallocates non-permanent data

    call orbital_free(orbs_ob=orbs_ob,orb_xc=orbs_xc)
    call orbital_shutdown()
  end subroutine xcfit_deallocate
  !*******************************************************

  !*******************************************************
  subroutine xcfit_setup()
    ! purpose : routine performs the necessary allocations
    !           additionally some informations from
    !           symmetry_data_module are stored in module private
    !           variables
    !** End of interface **************************************
    use symmetry_data_module ! description of irreps
    use fit_coeff_module
    use xc_hamiltonian
    use machineparameters_module, only: machineparameters_veclen
    use xc_cntrl, only: is_on, xc_gga, xc_vwn, xc_xalpha
    use density_calc_module, only: density_calc_setup
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind) :: alloc_stat

    if(is_on(xc_gga))call error_handler( &
         "GGA functional not possible for fitted V_xc potentials")
    ASSERT(is_on(xc_vwn))
    ASSERT(is_on(xc_xalpha))

    vec_length=machineparameters_veclen

    n_xc=fit_coeff_n_xc()
    ispin=ssym%n_spin
    n_irrep=ssym%n_irrep
    n_rhs=ispin+1
    n_lhs=n_rhs
    rhs_length=n_rhs*n_xc
    lhs_length=n_xc*(n_xc+1)/2
    allocate( &
         dims(n_irrep),&
         partners(n_irrep), &
         rho(vec_length,ispin),&
         fxc(vec_length),&
         eps(vec_length),&
         rhs(n_xc,n_rhs),&
         lhs(lhs_length,n_lhs), &
         dfdrho(vec_length,ispin),&
         wts(vec_length,n_rhs),&
         stat=alloc_stat)
    if(alloc_stat/=0) call error_handler('xcfit_setup: allocation failed')
    dims=ssym%dim
    partners=ssym%partner

    call density_calc_setup()
  end subroutine xcfit_setup

  !***************************************************************

  subroutine xcfit_close()
    ! purpose : perform the deallocation
    use density_calc_module, only: density_calc_close
    implicit none
    !** End of interface *****************************************
    integer(kind=i4_kind) :: alloc_stat
    deallocate(dims,partners,rho,fxc,eps,rhs,lhs,dfdrho,wts,stat=alloc_stat)
    if(alloc_stat/=0) call error_handler('xcfit_close deallocation')
    call density_calc_close()
  end subroutine xcfit_close

  !***************************************************************

  subroutine xcfit_clear()
    ! purpose : set some variables to zero
    !** End of interface *****************************************
    rho=0.0_r8_kind
    dfdrho=0.0_r8_kind
    fxc=0.0_r8_kind
    eps=0.0_r8_kind
  end subroutine xcfit_clear

  !***************************************************************

   subroutine xcfit_send_tomaster()
     ! purpose : send hamiltonians from the slave to the master
     !** End of interface ***************************************
     integer(kind=i4_kind) :: info

     call comm_init_send(comm_master_host,msgtag_xcfitham_send)
     call commpack(charge_int,info)
     if(info/=0) call error_handler&
          ('Error packing the charge in su xcfit_send_tomaster')
     call commpack(reshape(lhs,(/size(lhs)/)),lhs_length*n_lhs,1,info)
     if(info/=0) call error_handler&
          ('Error packing the  left hand side of xc-fit les')
     call commpack(reshape(rhs,(/rhs_length/)),rhs_length,1,info)
     if(info/=0) call error_handler&
          ('Error packing the rigth side of xc-fit')
     call comm_send()
   end subroutine xcfit_send_tomaster

   !***************************************************************

   subroutine xcfit_receive()
   ! purpose : receive the hamiltonians from the slaves
   !** End of interface *****************************************
     real(kind=r8_kind)  :: help_arr(lhs_length*n_lhs),help_arr2(rhs_length),&
          help_real
     integer(kind=i4_kind) :: i,info
     do i=1,comm_get_n_processors()-1
        call comm_save_recv(comm_all_other_hosts,msgtag_xcfitham_send)
        if(comm_msgtag()/=msgtag_xcfitham_send) call error_handler&
             ('Wrong msgtag in su xcfit_receive')
        call communpack(help_real,info) ! unpacking integrated charge
        if(info/=0) call error_handler&
             ('Error unpacking the charge in su xcfit_receive')
        charge_int=charge_int+help_real
        call communpack(help_arr,lhs_length*n_lhs,1,info) ! unpacking lhs
        if(info/=0) call error_handler&
             ('Error unpacking the lhs of xc-fit les')
        lhs=lhs+reshape(help_arr,(/lhs_length,n_lhs/))
        call communpack(help_arr2,rhs_length,1,info) ! unpacking rhs
        rhs=rhs+reshape(help_arr2,(/size(rhs,1),size(rhs,2)/))
        if(info/=0) call error_handler&
             ('Error unpacking the rhs of xc-fit les')
    end do
    write(output_unit,'("The numerically integrated charge is:",F20.16)') &
         charge_int
   end subroutine xcfit_receive

   !***************************************************************

   subroutine build_xcfit(loop)
     ! purpose : main routine for fitting the xc-hamiltonian
     !           XC-Hamiltonian
     !
     use vwnc
     use time_module
     use fit_coeff_module
     use options_module
     use density_calc_module, only: density_calc
     use grid_module, only: more_grid, grid_loop_setup
     implicit none
     integer(kind=i4_kind),optional :: loop ! number of the current scf cycle
     !** End of interface *****************************************
     real(kind=r8_kind),parameter :: deps=1.0e-50_r8_kind
     integer(kind=i4_kind) :: i,j,s,vec_length_act,counter
     logical:: use_model_density
     real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)

     use_model_density = options_xcmode() == xcmode_model_density .or. &
                         options_xcmode() == xcmode_extended_mda
     call xcfit_allocate()
     call xcfit_clear()
     lhs=0.0_r8_kind
     rhs=0.0_r8_kind
     charge_int=0.0_r8_kind
     exc_int=0.0_r8_kind

     call grid_loop_setup()
     ! loop over gridpoints:
     do while( more_grid(vec_length,grdpts,grdwts) ) ! fetching part of the grid
        ! more_grid() will return
        !   grdpts => new batch of coordinates
        !   grdwts => corresponding weights
        ! not more than "vec_length" long

        ! this may be, in general, below vec_length, particularly in last
        ! iteration:
        vec_length_act = size(grdpts,1)

        ! zero rho(:), fxc(:,), etc. FIXME: use local vars!
        call xcfit_clear()

        ! now calculating the orbitals
        call orbital_calculate(grdpts(1:vec_length_act,1:3),&
             vec_length_act,orbs_ob=orbs_ob,orb_xc=orbs_xc)
        ! now calculating the density
        call density_calc(vec_length_act,rho,orbs_ob)

        !
        ! FIXME: use xc_functionals(...) here, until then always do VWN ...
        !
        ! now evaluating the functionals
        if ( use_model_density ) then
            call vwn_calcMDA(rho,dfdrho,ispin,fxc,vec_length_act,eps)
            call xalpha_calcMDA(rho,dfdrho,ispin,fxc,vec_length_act,eps)
        else
            call vwn_calc(rho,dfdrho,ispin,fxc,vec_length_act,eps)
            call xalpha_calc(rho,dfdrho,ispin,fxc,vec_length_act,eps)
        endif

        ! performing integration over the grid
        charge_int=charge_int+sum(rho(1:vec_length_act,:)*spread&
             (grdwts(1:vec_length_act),2,ispin))
        exc_int=exc_int+sum(fxc(1:vec_length_act)*grdwts(1:vec_length_act))
        counter=1
        ! now the weights are scaled
        wts(1:vec_length_act,1:ispin)=spread(grdwts(1:vec_length_act),2,&
             ispin)*abs(rho(1:vec_length_act,1:ispin)/&
             (dfdrho(1:vec_length_act,1:ispin)+deps))
        wts(1:vec_length_act,ispin+1)=grdwts(1:vec_length_act)&
             *abs(sum(rho(1:vec_length_act,:),2)/&
             (eps(1:vec_length_act)+deps))
        !   now building up the left andf right sides of the linear equation
        !   system
        !   We have ispin+1 left and right sides
        !   The ispin+1 eqs is for fitting  e_xc
        do i=1,n_xc
           ! Actually, (-1)*<V_xc|g_k> and (-1)*<eps_xc|g_k> is loaded here
           do s=1,ispin
              rhs(i,s)=rhs(i,s)-sum(dfdrho(1:vec_length_act,s)*wts&
                   (1:vec_length_act,s)*orbs_xc%o(1:vec_length_act,i,1))
           enddo
           rhs(i,ispin+1)=rhs(i,ispin+1)-sum(eps(1:vec_length_act)*wts&
                (1:vec_length_act,&
                ispin+1)*orbs_xc%o(1:vec_length_act,i,1))
           do j=1,i
              do s=1,ispin+1
                 lhs(counter,s)=lhs(counter,s)+sum&
                      (orbs_xc%o(1:vec_length_act,i,1)*&
                      orbs_xc%o(1:vec_length_act,j,1)*wts&
                      (1:vec_length_act,s))
              end do
              counter=counter+1
           enddo
        enddo
        !       call add_timetype(tt,build_xc_time)
     enddo

     if(comm_i_am_master()) then
        ! receiving parts of the les from the slaves
        call xcfit_receive()
        ! now solving linear equation system in order to get
        ! fitting coefficients
        call copy_coeff_xc_old
        do s=1,ispin
           call decompose_linsys(n_xc,lhs(:,s))
           call solve_linsys(n_xc,lhs(:,s),rhs(:,s))
        enddo
        call decompose_linsys(n_xc,lhs(:,ispin+1))
        call solve_linsys(n_xc,lhs(:,ispin+1),rhs(:,ispin+1))
        ! Actually, coeff_xc [c_k] and coeff_xc_en [d_k] are loaded here such
        ! that V_xc = (-1)*Sum(k) c_k g_k and eps_xc = (-1)*Sum(k) d_k g_k
        coeff_xc=rhs(:,1:ispin)
        xc_initialized = .false.
        coeff_en_xc=rhs(:,ispin+1)
        ! now mixing the xc-coefficients
        call mixing_xc(loop-1) ! iter=loop-1 instead of iter=loop : in order
                               ! to stay consistent with the V14alpha version
                               ! of the PARAGAU program, where build_xcfit was
                               ! still located at the  e n d  of the SCF loop.
                               ! This way iter=1 at the first call of mixing_xc
        ! the old exchange coefficients are not needed anymore
        call free_coeff_xc_old
     else
        ! sending parts of the les to the master
        call xcfit_send_tomaster()
     end if
     call xcfit_deallocate()

   end subroutine build_xcfit

   !***************************************************************

 end module xcfit_hamiltonian

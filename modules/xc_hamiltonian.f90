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
module xc_hamiltonian
!---------------------------------------------------------------
!
!  Purpose: This module creates the XC-part of the hamiltionian by
!           numerical integration over the grid
!
!           The calculation is accordiarningng to
!  References: Johnson,Gill,Pople,J. Chem. Phys. 98 (7) 1.4.1993
!
!   Module called by: main_scf,energy_calc_module
!
!           The module is organized as follows:
!
!           There  is  one startup  routine  called xc_setup().   This
!           routine is by every worker.
!
!           There are routines available for reading and writing the input
!           and packing and unpacking this data (xc_input_read,
!           xc_input_write,xc_input_bcast)
!
!           the routine xc_clean is called internally between scf_cycles
!           in order to reset variables to zero.
!
!           The main routine is called build_xc ( the wrapper build_xc_main)
!           In this routine the hamiltonian is stored in one linear array
!           The following tasks are performed:
!            ! calculation of the orbitals and orbitalgradients
!            ! calculation of the density and the density gradient
!            ! evaluation of the xc_functionals
!            ! numerical integration of the charge over the grid
!            ! numerical integration of the xc_energy over the grid
!            ! assembling of the hamiltonian
!            ! sending the hamiltonian to the master (slave)
!            ! receiving the hamiltonians from the slaves (master)
!            ! perform the mixing of the new hamiltonian with the old
!            ! hamiltonian
!            ! resorting the elements of the hamiltonian to
!            ! the normal quadratic structure
!
!  Author: MS
!  Date: 2/96
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: MS
! Date:   10/96
! Description: Performance has been improved by rewriting
!              the integration of matrixelements
!
! Modification (Please copy before editing)
! Author: FN
! Date:   11/96
! Description: The ability to continue an old calculation
!              using the XC-hamiltonian (and the fitcoefficients)
!              has been implemented.
!              The organization is as follows:
!
! Modification (Please copy before editing)
! Author: TG
! Date:   12/96
! Description: The parts of the Hamiltonian are received in
!              a cascadic prodedure.
!
! Modification (Please copy before editing)
! Author: MS
! Date:   8/97
! Description: A bug concerning the treatment of pseudo 2D irreps
!              has been fixed
!
! Modification (Please copy before editing)
! Author: MM
! Date:   10/97
! Description: extension to spin orbit
!
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!
!------------ Modules used --------------------------------------
!---------------------------------------------------------------
!

# include "def.h"
  use type_module  ! type specification parameters
  use orbitalstore_module
  use orbital_module
  use density_data_module
  use comm_module
  use hamiltonian_module
  use occupied_levels_module
  use mixing_module, only: mixing_ham, mixing_discard_init
! use density_calc_module ! use in actual subs only!
  use machineparameters_module
  use options_module, only: options_recover, &
                            recover_scfstate,  &
                            options_xcmode, xcmode_numeric_exch, &
                            options_n_spin,options_spin_orbit,   &
                            options_orbitals_in_memory,          &
                            options_orbitals_on_file
  use unique_atom_module
  use orbitalprojection_module, only: orbitalprojection_ob
  use spin_orbit_module, only: spin_orbit_polarized
  use strings, only: strlen=>stringlength_string
  implicit none
  private
  save

!== Interrupt end of public interface of module =================

!------------ public variables ----------------------------------
  real(kind=r8_kind), allocatable, public :: ham_xc_arr(:)
  ! Elements of the kohn-Sham-Matrix
  real(kind=r8_kind), allocatable, public :: ham_xc_arr_real(:),ham_xc_arr_imag(:)
  ! Elements of the kohn-Sham-Matrix (in case of spin orbit the KS-Matrix is complex)
  logical, public :: mat_initialized, matold_initialized
  real(kind=r8_kind),public  :: s_average

!------------ public functions and subroutines ------------------
  public :: xc_setup, build_xc_main, build_xc, xc_close, &
       & xc_hamiltonian_store, xc_hamiltonian_recover,&
       & xc_get_exc

!================================================================
! End of public interface of module
!================================================================

  integer(i4_kind), parameter, private ::&
       & X=1, Y=2, Z=3,&
       & UP=1, DN=2,&
       & UPUP=1, DNDN=2, UPDN=3,&
       & RE=1, IM=2

  real(r8_kind), parameter ::&
         QUARTER = 0.25_r8_kind ,&
         HALF = 0.5_r8_kind     ,&
         ONE  = 1.0_r8_kind     ,&
         TWO  = 2.0_r8_kind     ,&
         ZERO = 0.0_r8_kind

  real(kind=r8_kind),allocatable :: rho(:,:),&
         dfdrho(:,:),&  ! derivative of f with respect to rho
         ham_xc_arr_old(:), ham_xc_arr_old_real(:),ham_xc_arr_old_imag(:),&
          !!$         & help_arr(:),&
         & help_arr_real(:),help_arr_imag(:),help_arr2(:,:),& !<< GET RID OF THEM, they must be local vars
         & fxc(:),fxc_gga(:),&        ! funktion f according to JPG
         gamma(:,:),& ! norms of the density gradients
         ! 1. column : grarho(:,1)**2
         ! 2. column : grarho(:,2)**2
         ! 3. column : grarho(:,1)*grarho(:,2)
         dfdgrarho(:,:)  ,& ! derivatives of f with respect to gamma
         tau(:,:)        ,& ! kinetic energy density
         dfdtau(:,:)        ! derivatives of f with respect to tau
#ifdef WITH_CORE_DENS
  real(kind=r8_kind),allocatable :: rho_core(:,:), grad_rho_core(:,:,:)
#endif
  ! in case of spin polarized spin orbit we need
  real(kind=r8_kind),allocatable :: rho_ud(:,:),help_arr_uu_real(:),help_arr_uu_imag(:),&
       help_arr_ud_real(:),help_arr_ud_imag(:),&
       help_arr_dd_real(:),help_arr_dd_imag(:),&
       help_arr_du_real(:),help_arr_du_imag(:),&
       s_abs(:),s_dens(:,:),dfdn(:),dfds(:),& !s_exp(:),&
       grarho_udx(:,:),grarho_udy(:,:),grarho_udz(:,:),&
       gras_densx(:),gras_densy(:),gras_densz(:)

! These were global module vars, now made local subroutine vars:
! real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)
  integer(kind=i4_kind) :: ispin,vec_length

  integer(kind=i4_kind),pointer :: dims(:),partners(:)

  integer(kind=i4_kind) ::&
       & xc_length      = -1,&
       & xc_length_vec  = -1,&
       & xc_length_proj = -1
  integer(kind=i4_kind) ::&
       & n_irrep      = -1,&
       & n_vec_irrep  = -1,&
       & n_proj_irrep = -1
  integer(kind=i4_kind),pointer :: vec_dims(:), vec_partners(:)
  integer(kind=i4_kind),pointer :: proj_dims(:),proj_partners(:)

  logical, allocatable :: pseudo(:)
  real(kind=r8_kind)  :: charge_int,exc_int
  real(kind=r8_kind)  :: spin_vector(3)
  type(orbital_type),pointer :: orbs_ob(:)
  type(spinor_type),pointer  :: orbs_spinor_ob(:)
  type(orbital_gradient_type),pointer :: orbs_grads(:)
  type(spinor_gradient_type),pointer :: spinor_grads(:)
  real(kind=r8_kind)  :: exc_gga_int
#ifdef WITH_CORE_DENS
  type(core_orbital_type)          :: orbs_ob_core
  type(core_orbital_gradient_type) :: grads_core
#endif

contains


  subroutine build_xc_main(loop)
    ! purpose : wrapper for build_xc; it runs only on the master and sends
    !           the message "execute build_xc" to the slaves. Subsequently
    !           build_xc is called. build_xc_main is called in every scf-
    !           cycle by main_scf.
    ! ------ Modules -------------------------------------------
    use msgtag_module, only: msgtag_build_xc
    implicit none
    integer(kind=i4_kind),optional :: loop ! number of the actual scf-cycle
                                           ! (within the current SCF run)
    !** End of interface *****************************************

    if ( comm_parallel() ) then
       call comm_init_send(comm_all_other_hosts, msgtag_build_xc)
       call comm_send()
    endif

    call build_xc(loop=loop)
  end subroutine build_xc_main

  !***********************************************************

  subroutine xc_allocate()
    use xc_cntrl, only: is_on, xc_so_spatial, xc_nl_calc
    implicit none
    ! purpose: allocates non-permanent data
    !** End of interface **************************************
    if(.not.xc_nl_calc) then
       ! initialise orbitals_module
       call orbital_setup(machineparameters_veclen)
       if (options_spin_orbit) then
          DPRINT 'xch/xc_allocate: allocate spinors'
          call orbital_allocate(orbs_spinor_ob=orbs_spinor_ob)
          if(is_on(xc_so_spatial))then
             DPRINT 'xch/xc_allocate: ... and orbitals'
             call orbital_allocate(orbs_ob)
          end if
       else
          call orbital_allocate(orbs_ob)
       endif
#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
          call core_orbital_allocate(orbs_ob_core)
       end if
#endif
    else
       ! initialise orbitals_module
       call orbital_setup(machineparameters_veclen,do_gradients=.true.)
       if (options_spin_orbit) then
          call orbital_allocate(orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
          DPRINT 'xch/xc_allocate: allocate spinors and grads'
          if(is_on(xc_so_spatial))then
             DPRINT 'xch/xc_allocate: ... and orbitals and grads'
             call orbital_allocate(orbs_ob,orbs_grads)
          endif
       else
          call orbital_allocate(orbs_ob,orbs_grads)
       endif

#ifdef WITH_CORE_DENS
       if (pseudopot_present.and.core_density_setup) then
          call core_orbital_allocate(orbs_ob_core=orbs_ob_core,&
                                     grads_core=grads_core)
       end if
#endif
    endif
  end subroutine xc_allocate

  !*******************************************************

  subroutine xc_deallocate()
    use xc_cntrl, only: is_on, xc_so_spatial, xc_nl_calc
    ! purpose: deallocates non-permanent data
    !** End of interface **************************************

    if(.not.xc_nl_calc) then
       if (options_spin_orbit) then
          call orbital_free(orbs_spinor_ob=orbs_spinor_ob)
          if(is_on(xc_so_spatial))then
             call orbital_free(orbs_ob)
          endif
       else
          call orbital_free(orbs_ob)
#ifdef WITH_CORE_DENS
          if (pseudopot_present.and.core_density_setup) then
             call core_orbital_free(orbs_ob_core)
          end if
#endif
       endif
       call orbital_shutdown()
    else
       if (options_spin_orbit) then
          call orbital_free(orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
          if(is_on(xc_so_spatial))then
             call orbital_free(orbs_ob,orbs_grads)
          endif
       else
          call orbital_free(orbs_ob,orbs_grads)
       endif
       call orbital_shutdown()

#ifdef WITH_CORE_DENS
      if (pseudopot_present.and.core_density_setup) then
         call core_orbital_free(orbs_ob_core,grads_core)
      end if
#endif
    endif
  end subroutine xc_deallocate

  !*******************************************************

  subroutine xc_setup()
    ! purpose : routine performs the necessary allocations
    !           additionally some informations from
    !           symmetry_data_module are stored in module private
    !           variables
    !------------ Modules used -----------------------------------
    use error_module
    use operations_module, only: operations_geo_opt, operations_qm_mm
    use symmetry_data_module ! description of irreps
    use xc_cntrl, only: xc_nl_calc, xc_is_on=>is_on, xc_so_spatial, xc_gga, xc_mgga
    use spin_orbit_module, only: so_is_on=>is_on, op_polarized
    use density_calc_module, only: density_calc_setup
#ifdef WITH_CORE_DENS
    use density_calc_module, only: fitted_core_dens_calc_setup
#endif
    implicit none
    !** End of interface **************************************

    integer(kind=i4_kind) :: i,alloc_stat
    !------------------ Executable code ---------------------------

    DPRINT 'xch/xc_setup: entered'

    ! Initialize orbitalstore module
    if(options_orbitals_in_memory() .or. options_orbitals_on_file() ) then
       call orbitalstore_setup(machineparameters_veclen)
    end if

    if(options_spin_orbit)then
       if(so_is_on(op_polarized))then
          if(xc_is_on(xc_so_spatial))&
               & call error("xch/xc_setup: conflicting: Polarized SO & Spatial SO")
          if(xc_is_on(xc_gga))&
               & call error("xch/xc_setup: not implemented: Polarized SO & GGA")
       endif
    endif

    ! reset the exchange energy
    exc_int = 0.0_r8_kind
    vec_length=machineparameters_veclen

    ispin=ssym%n_spin
    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       exc_gga_int = 0.0_r8_kind
       n_proj_irrep = ssym%n_proj_irrep
       allocate(proj_dims(n_proj_irrep),proj_partners(n_proj_irrep), &
            stat=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('allocation (1) failed in su xc_setup')
       proj_dims     = ssym%dim_proj
       proj_partners = ssym%partner_proj
    endif
!!$    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_vec_irrep = ssym%n_irrep
       allocate(vec_dims(n_vec_irrep),vec_partners(n_vec_irrep),pseudo(n_vec_irrep), &
            stat=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('allocation (1) failed in su xc_setup')
       vec_dims     = ssym%dim
       vec_partners = ssym%partner
       pseudo=ssym%pseudo
!!$    endif ! options_spin_orbit

    xc_length_vec=0
    do i=1,n_vec_irrep
       xc_length_vec = xc_length_vec + (vec_dims(i)*(vec_dims(i)+1))/2
    end do
    xc_length_vec = xc_length_vec*ispin

    DPRINT 'xch/xc_setup: xc_length_vec =',xc_length_vec
    DPRINT 'xch/xc_setup: vec_dims      =',vec_dims
    DPRINT 'xch/xc_setup: vec_partners  =',vec_partners

    if(options_spin_orbit)then
       xc_length_proj=0
       do i=1,n_proj_irrep
          xc_length_proj = xc_length_proj + (proj_dims(i)*(proj_dims(i)+1))/2
       end do
!!$       xc_length_proj = xc_length_proj * ispin
    endif

    if(options_spin_orbit)then
       n_irrep   =  n_proj_irrep
       xc_length =  xc_length_proj
       dims      => proj_dims
       partners  => proj_partners

    DPRINT 'xch/xc_setup: xc_length_proj =',xc_length_proj
    DPRINT 'xch/xc_setup: proj_dims      =',proj_dims
    DPRINT 'xch/xc_setup: proj_partners  =',proj_partners

    else
       n_irrep   =  n_vec_irrep
       xc_length =  xc_length_vec
       dims      => vec_dims
       partners  => vec_partners
    endif

    DPRINT 'xch/xc_setup: xc_length =',xc_length
    DPRINT 'xch/xc_setup: dims      =',dims
    DPRINT 'xch/xc_setup: partners  =',partners

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if (spin_orbit_polarized) then
          !
          ! OPEN SHELL
          !
          allocate(dfdrho(vec_length,2),rho(vec_length,ispin),fxc(vec_length),&
               dfdn(vec_length),dfds(vec_length),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation failed (2) in su xc_setup')
!!$          allocate(s_dens(vec_length,3),s_abs(vec_length),rho_ud(vec_length,2),&
!!$               s_exp(vec_length),stat=alloc_stat)
          allocate(s_dens(vec_length,3),s_abs(vec_length),rho_ud(vec_length,2),&
               & stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation of s_dens failed in su xc_setup')
          allocate(help_arr_uu_real(vec_length),help_arr_uu_imag(vec_length),&
               help_arr_ud_real(vec_length),help_arr_ud_imag(vec_length),&
               help_arr_du_real(vec_length),help_arr_du_imag(vec_length),&
               help_arr_dd_real(vec_length),help_arr_dd_imag(vec_length),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (3) failed in su xc_setup')
          if (xc_nl_calc) then
             allocate(fxc_gga(vec_length),gamma(vec_length,3),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation failed (4) in su xc_setup')
             allocate(gras_densx(vec_length),gras_densy(vec_length),gras_densz(vec_length),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation of gras_dens failed in su xc_setup')
             allocate(grarho_udx(vec_length,2),grarho_udy(vec_length,2),grarho_udz(vec_length,2),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation of grarho_ud failed in su xc_setup')
             allocate(dfdgrarho(vec_length,3),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation of dfdgrarho failed in su xc_setup')
          endif
       else
          !
          ! CLOSED SHELL
          !
          allocate(dfdrho(vec_length,ispin),rho(vec_length,ispin),fxc(vec_length),&
               stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation failed (2) in su xc_setup')
          allocate(help_arr_real(vec_length),help_arr_imag(vec_length),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (3) failed in su xc_setup')
          if (xc_nl_calc) then
             allocate(gras_densx(vec_length),gras_densy(vec_length),gras_densz(vec_length),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation of gras_dens failed in su xc_setup')
             allocate(fxc_gga(vec_length),gamma(vec_length,1),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation failed (4) in su xc_setup')
             allocate(dfdgrarho(vec_length,1),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ('allocation of dfdgrarho failed in su xc_setup')
          endif
       endif

#ifdef WITH_CORE_DENS
    if (pseudopot_present.and.core_density_setup) then
      allocate(grad_rho_core(vec_length,3,ispin),rho_core(vec_length,ispin),&
           stat=alloc_stat)
      if(alloc_stat/=0) call error_handler&
         ('allocation failed (2core) in su xc_setup')
    end if
#endif
       if(comm_i_am_master()) then
          allocate(ham_xc_arr_real(xc_length),ham_xc_arr_imag(xc_length),&
               stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (3) failed in su xc_setup')
          if (options_recover() /= recover_scfstate .or. &
              operations_geo_opt .or. operations_qm_mm) then
             ham_xc_arr_real = 0.0_r8_kind
             ham_xc_arr_imag = 0.0_r8_kind
             mat_initialized = .true.
          end if
       endif
    else ! options_spin_orbit
       allocate(dfdrho(vec_length,ispin),rho(vec_length,ispin),fxc(vec_length),&
            stat=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('allocation failed (2) in su xc_setup')
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
!!$       allocate(help_arr(vec_length),stat=alloc_stat)
       if(comm_i_am_master()) then
          ! the old XC matrix now is a temporary array of build_xc
          ! and the XC matrix is no longer allocated in main_scf
          allocate(ham_xc_arr(xc_length),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (3) failed in su xc_setup')
          if (options_recover() /= recover_scfstate .or. &
              operations_geo_opt .or. operations_qm_mm) then
             ham_xc_arr = 0.0_r8_kind
             mat_initialized = .true.
          end if
       end if
       if(xc_nl_calc) then
          allocate(gamma(vec_length,1+2*(ispin-1)),STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (4) failed in su xc_setup')
          allocate(dfdgrarho(vec_length,1+2*(ispin-1)),STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('allocation (5) failed in su xc_setup')
             !-----------------------------------------------------------+
             ! Allocation of MGGA variables                              |
             !-----------------------------------------------------------+
             if (xc_is_on(xc_mgga)) then
                allocate(tau(vec_length,ispin),STAT=alloc_stat)
                if(alloc_stat/=0) call error_handler&
                   ('allocation (7) failed in su xc_setup')
                allocate(dfdtau(vec_length,ispin),STAT=alloc_stat)
                if(alloc_stat/=0) call error_handler&
                   ('allocation (8) failed in su xc_setup')
             end if
       endif
    endif !  options_spin_orbit
    ! for a NLDA calculation
    if(xc_nl_calc) then
       !??? gamma, dgdgrarho ???
       allocate(help_arr2(vec_length,3),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ('allocation (6) failed in su xc_setup')
    endif
    call density_calc_setup()

#ifdef WITH_CORE_DENS
    if ( pseudopot_present.and.core_density_setup) then
       call fitted_core_dens_calc_setup()
    end if
#endif
   end subroutine xc_setup

   !***************************************************************

   subroutine xc_close()
     ! purpose : perform the deallocation
     !** End of interface *****************************************
     use comm_module, only: comm_i_am_master
     use xc_cntrl, only: xc_nl_calc, xc_is_on=>is_on, xc_mgga
     use density_calc_module, only: density_calc_close
#ifdef WITH_CORE_DENS
     use density_calc_module, only: fitted_core_dens_calc_close
#endif
     implicit none
     ! *** end of interface ***

     integer(kind=i4_kind):: alloc_stat

     ! Close orbitalstore_module
     call orbitalstore_shutdown()
     if (options_spin_orbit) then
        if (xc_nl_calc) then
           deallocate(fxc_gga,stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ('deallocation (1a) failed in su xc_close')
        endif
        if (spin_orbit_polarized) then
           !
           ! OPEN SHELL
           !
!!$           deallocate(s_dens,s_abs,rho_ud,dfdn,dfds,s_exp,stat=alloc_stat)
           deallocate(s_dens,s_abs,rho_ud,dfdn,dfds,stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ('deallocation of s_dens failed in su xc_close')
           if (xc_nl_calc) then
              deallocate(gras_densx,gras_densy,gras_densz,stat=alloc_stat)
              if(alloc_stat/=0) call error_handler &
                   ('deallocation of gras_dens failed in su xc_close')
              deallocate(grarho_udx,grarho_udy,grarho_udz,stat=alloc_stat)
              if(alloc_stat/=0) call error_handler &
                   ('deallocation of grarho_ud failed in su xc_close')
           endif
           deallocate(help_arr_uu_real,help_arr_uu_imag,&
                help_arr_ud_real,help_arr_ud_imag,&
                help_arr_du_real,help_arr_du_imag,&
                help_arr_dd_real,help_arr_dd_imag,&
                stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ('deallocation (1a spinor) failed in su xc_close')
        else
           !
           ! CLOSED SHELL
           !
           deallocate(help_arr_real,help_arr_imag,&
                stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ('deallocation (1a spinor) failed in su xc_close')
        endif
     else
!!$        deallocate(help_arr,&
!!$             stat=alloc_stat)
!!$        if(alloc_stat/=0) call error_handler &
!!$             ('deallocation (1a) failed in su xc_close')
     endif

!!$     deallocate(rho,dfdrho,fxc,dims,partners,pseudo,&
!!$          stat=alloc_stat)
     deallocate(rho,dfdrho,fxc,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ('deallocation (1) failed in su xc_close')

     deallocate(vec_dims,vec_partners,pseudo,stat=alloc_stat)
     if(alloc_stat/=0) call error_handler &
          ('deallocation (1.1) failed in su xc_close')
     if(options_spin_orbit)then
        deallocate(proj_dims,proj_partners,stat=alloc_stat)
        if(alloc_stat/=0) call error_handler &
             ('deallocation (1.2) failed in su xc_close')
     endif
#ifdef WITH_CORE_DENS
     if (pseudopot_present.and.core_density_setup) then
       deallocate(grad_rho_core,rho_core,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler &
          ('deallocation (1core) failed in su xc_close')
     end if
#endif
     if (xc_nl_calc) then
        deallocate(help_arr2,STAT=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("xc_close : deallocation (2) failed")
        deallocate(dfdgrarho,STAT=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("xc_close : deallocation (3) failed")
        deallocate(gamma,STAT=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("xc_close : deallocation (4) failed")
        !---------------------------------------------------------------+
        ! Deallocation of MGGA variables                                |
        !---------------------------------------------------------------+
        if (xc_is_on(xc_mgga)) then
           deallocate(tau,STAT=alloc_stat)
           if(alloc_stat/=0) call error_handler&
              ('xc_close : deallocation (6,mgga) failed')
           deallocate(dfdtau,STAT=alloc_stat)
           if(alloc_stat/=0) call error_handler&
              ('xc_close : deallocation (7,mgga) failed')
        end if
     endif

     if (comm_i_am_master()) then
        ! the old XC matrix now is a temporary array of build_xc
        ! and the XC matrix is no longer allocated in main_scf
        if (options_spin_orbit) then
           deallocate(ham_xc_arr_real,ham_xc_arr_imag,stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ("xc_close : deallocation (5 spinor) failed ")
        else
           deallocate(ham_xc_arr,stat=alloc_stat)
           if(alloc_stat/=0) call error_handler &
                ("xc_close : deallocation (5) failed ")
        endif
     endif
     call density_calc_close()
#ifdef WITH_CORE_DENS
     if (pseudopot_present.and.core_density_setup) then
        call fitted_core_dens_calc_close()
     end if
#endif
   end subroutine xc_close

   !***************************************************************

   subroutine xc_clear()
     !
     ! FIXME: Initialize vars explicitly!
     !        Make them local vars of respective subs!
     !        Yet better dont rely on their state at all!
     !
     use xc_cntrl, only: xc_nl_calc
     implicit none
     ! purpose : set some variables to zero
     !** End of interface *****************************************
     rho=0.0_r8_kind
     dfdrho=0.0_r8_kind
#ifdef WITH_CORE_DENS
     if (pseudopot_present.and.core_density_setup) then
        rho_core = 0.0_r8_kind
        grad_rho_core = 0.0_r8_kind
     end if
#endif
     fxc=0.0_r8_kind
     if (options_spin_orbit) then
        if (xc_nl_calc) then
           fxc_gga = 0.0_r8_kind
        endif
        if (spin_orbit_polarized) then
           !
           ! OPEN SHELL
           !
           rho_ud = 0.0_r8_kind
           help_arr_uu_real=0.0_r8_kind
           help_arr_uu_imag=0.0_r8_kind
           help_arr_ud_real=0.0_r8_kind
           help_arr_ud_imag=0.0_r8_kind
           help_arr_du_real=0.0_r8_kind
           help_arr_du_imag=0.0_r8_kind
           help_arr_dd_real=0.0_r8_kind
           help_arr_dd_imag=0.0_r8_kind
           dfdn=0.0_r8_kind
           dfds=0.0_r8_kind
           s_dens=0.0_r8_kind
!!$           s_exp = 0.0_r8_kind
        else
           !
           ! CLOSED SHELL
           !
           help_arr_real=0.0_r8_kind
           help_arr_imag=0.0_r8_kind
        endif
     else
!!$        help_arr=0.0_r8_kind
     endif
! These were global module vars, now made local subroutine vars:
!    nullify(grdpts,grdwts)
     if(xc_nl_calc) then
        gamma=0.0_r8_kind
        dfdgrarho=0.0_r8_kind
     endif
   end subroutine xc_clear

   !***************************************************************

   subroutine xc_sndrcv()
     ! purpose : reduces the hamiltonians from all processors
     use xc_cntrl, only: xc_nl_calc
     use comm, only: comm_reduce
     implicit none
     !** End of interface *****************************************

    call comm_reduce(charge_int)
    call comm_reduce(exc_int)

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if (spin_orbit_polarized) then
          call comm_reduce(spin_vector)
          call comm_reduce(s_average)
       endif
       if (xc_nl_calc) then
          call comm_reduce(exc_gga_int)
       endif
       call comm_reduce(ham_xc_arr_real)
       call comm_reduce(ham_xc_arr_imag)
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call comm_reduce(ham_xc_arr)
    endif ! options_spin_orbit
   end subroutine xc_sndrcv

   !***************************************************************

   function xc_get_exc()
     ! purpose : make the module private integrated xc_energy
     !           available to the calling unit
     !** End of interface *****************************************
     real(kind=r8_kind) :: xc_get_exc
     xc_get_exc=exc_int
   end function xc_get_exc

   !***************************************************************

    subroutine calc_xcks(ham_xc_arr, variant)
      !
      ! was a piece of code, now a module subroutine
      ! determine Kohn-Sham-Matrix without spin orbit
      !
      use timer_module
      use time_module, only: start_timer, stop_timer
      use xc_cntrl
      use f77_blas, only: dgemm
      use xc_func, only: xc_functionals
      use density_calc_module, only: density_calc, density_calc_nl
      use grid_module, only: more_grid, grid_loop_setup
      implicit none
      real(kind=r8_kind),intent(inout) :: ham_xc_arr(:) ! ...(xc_length_vec)
      character(len=*), intent(in)     :: variant ! "orbitals" or "spinors"
      ! controls how the density is evaluated, integration is always
      ! over plain orbital functions
      ! *** end of interface ***

      logical :: spin_orbit_density
      integer(i4_kind)      :: alloc_stat
      integer(kind=i4_kind) :: i,s,vla
      integer(kind=i4_kind) :: counter

      real(r8_kind), allocatable :: orbs_help(:,:), ham(:,:)
      ! These were global module vars, now made local subroutine vars:
      real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)
      !
      ! overwrite global variables:
      !
      integer(kind=i4_kind)         :: n_irrep
      integer(kind=i4_kind),pointer :: dims(:),partners(:)

      integer(kind=i4_kind)         :: maxdim ! maxval(dims(:))

      !
      ! Eclipse global module vars:
      !
      real(r8_kind) :: rho(vec_length,ispin)
      real(r8_kind) :: fxc(vec_length)
      real(r8_kind) :: dfdrho(vec_length,ispin)

      FPP_TIMER_DECL(tot)
      FPP_TIMER_DECL(orb)
      FPP_TIMER_DECL(dens)
      FPP_TIMER_DECL(xcfun)
      FPP_TIMER_DECL(intgr)
      FPP_TIMER_DECL(trist)

      DPRINT 'xch/calc_xcks: entered'
      FPP_TIMER_START(tot)

      select case(variant)
      case ("orbitals")
        spin_orbit_density = .false.
      case ("spinors")
        spin_orbit_density = .true.
      case default
        ABORT('no such case')
      end select

      ! use vector irrep dimension:
      n_irrep  =  n_vec_irrep
      dims     => vec_dims
      partners => vec_partners

      ! dimension for intermediate arrays, max of irrep dimensions:
      maxdim = maxval(dims)

      ! allocate space for intermediate arrays
      allocate(orbs_help(vec_length, maxdim), ham(maxdim, maxdim), stat=alloc_stat)
      ASSERT(alloc_stat==0)

      call grid_loop_setup() !static pre-distribution of jobs
      ! loop over gridpoints:
      do while( more_grid(vec_length, grdpts, grdwts) ) ! fetching part of the grid
         ! more_grid() will return
         !   grdpts => new batch of coordinates
         !   grdwts => corresponding weights
         ! not more than "vec_length" long
         ! it will return only false if the proc cannot steal anymore

         ! this may be, in general, below vec_length, particularly in last
         ! iteration:
         vla = size(grdpts,1)

         call start_timer(timer_grid_orbitals)

         FPP_TIMER_START(orb)
         if(.not.spin_orbit_density)then
            ! calculate only spatial orbitals:
            call orbital_calculate(grdpts(1:vla,1:3),vla,orbs_ob)
         else if(is_on(xc_so_orb_to_sporb))then
            ! calculte both at one time:
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_ob=orbs_ob,&
                 & orbs_spinor_ob=orbs_spinor_ob)
         else
            ! calculate both independently:
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_ob)
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_spinor_ob=orbs_spinor_ob)
         endif
         FPP_TIMER_STOP(orb)

         call stop_timer(timer_grid_orbitals)

         call start_timer(timer_grid_density)
         FPP_TIMER_START(dens)
         if(.not.spin_orbit_density)then
            if(options_spin_orbit) then
               call density_calc(vla,rho,orbs_ob)
            else
               ! FIXME: WHY?:
!              call density_calc(vla,rho,orbs_ob)
               call density_calc_nl(vla,rho,orbs_ob=orbs_ob)
            endif
         else
            call density_calc(vla,rho,orbs_spinor_ob=orbs_spinor_ob)
         endif
         FPP_TIMER_STOP(dens)

         ! performing integration over the grid
         charge_int=charge_int+sum(rho(1:vla,:)*spread&
              (grdwts(1:vla),2,ispin))
         call stop_timer(timer_grid_density)

         call start_timer(timer_grid_functionals)

         FPP_TIMER_START(xcfun)
         call xc_functionals(vla,ispin,rho,fxc,dfdrho)
         FPP_TIMER_STOP(xcfun)

         call stop_timer(timer_grid_functionals)

         ! performing integration over the grid
         exc_int = exc_int + sum(fxc(1:vla) * grdwts(1:vla))

         call start_timer(timer_grid_xcbuild)

         do s=1,ispin
            counter = s ! starting offset: 1 for s=1, and 2 for s=2

            dfdrho(1:vla,s) = dfdrho(1:vla,s) * grdwts(1:vla)
            do i=1,n_irrep
               FPP_TIMER_START(intgr)
               call integrate(dims(i), partners(i), dfdrho(:,s), orbs_ob(i)%o, ham, orbs_help)
               FPP_TIMER_STOP(intgr)
               FPP_TIMER_START(trist)
               call tristore(dims(i), counter, ispin, ham_xc_arr, ham)
               FPP_TIMER_STOP(trist)
            enddo! loop over irreps
         enddo! loop over spins

         ! FIXME: Initialize vars explicitly!
         ! NO MORE: call xc_clear()
         call stop_timer(timer_grid_xcbuild)
      enddo

      deallocate(orbs_help, ham, stat=alloc_stat)
      ASSERT(alloc_stat==0)

      FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
    print *,'calc_xcks: TIMING total          =',FPP_TIMER_VALUE(tot)
    print *,'calc_xcks: TIMING |-orbitals     =',FPP_TIMER_VALUE(orb)
    print *,'calc_xcks: TIMING |-density      =',FPP_TIMER_VALUE(dens)
    print *,'calc_xcks: TIMING |-xc func      =',FPP_TIMER_VALUE(xcfun)
    print *,'calc_xcks: TIMING |-integrate    =',FPP_TIMER_VALUE(intgr)
    print *,'calc_xcks: TIMING |-tri-store    =',FPP_TIMER_VALUE(trist)
#endif
      DPRINT 'xch/calc_xcks: exit'
    contains

      subroutine integrate(n, np, dfdrho, orbs, ham, orbs_help)
        integer(i4_kind), intent(in) :: n  ! irrep dimension
        integer(i4_kind), intent(in) :: np ! number of partners
        real(r8_kind), intent(in)    :: dfdrho(:) ! (>=vla)
        real(r8_kind), intent(in)    :: orbs(:,:,:)
        real(r8_kind), intent(out)   :: ham(:,:)  ! used only (:n,:n)
        real(r8_kind)                :: orbs_help(:,:) ! used only (:vla,:n), work array
        ! *** end of interface ***

        integer(i4_kind) :: m,j
        real(r8_kind)    :: alpha
        real(r8_kind)    :: beta

        alpha = 1.0_r8_kind/REAL(np,r8_kind)
        beta  = 0.0_r8_kind
        do m=1,np ! number of partners
           ! prepeare orbs_help
           do j=1,n
              orbs_help(1:vla,j) = orbs(1:vla,j,m) * dfdrho(1:vla)
           end do

           ! now do actual integration with the help of blas
           call dgemm('t','n', n, n, vla, alpha, &
                orbs_help(:,:)  , size(orbs_help,1), &
                orbs(:,:,m)     , size(orbs,1), beta, &
                ham(:,:)        , size(ham, 1))
           ! for other ms add up:
           beta = 1.0_r8_kind
        end do
      end subroutine integrate

      subroutine tristore(n, counter, ispin, ham_xc_arr, ham)
        integer(i4_kind), intent(in)    :: n ! irrep dimension
        integer(i4_kind), intent(inout) :: counter ! offset into ham_xc_arr
        integer(i4_kind), intent(in)    :: ispin
        real(r8_kind),    intent(inout) :: ham_xc_arr(:)
        real(r8_kind),    intent(in)    :: ham(:,:)
        ! *** end of interface ***

        integer(i4_kind) :: j,k

        do j=1,n
           do k=1,j
              ham_xc_arr(counter) = ham_xc_arr(counter) + ham(k,j)
              counter = counter + ispin ! step by 2 if polarized
           end do
        end do
      end subroutine tristore

    end subroutine calc_xcks

    subroutine build_xc(loop)
    ! purpose : main routine for calculation of the numerical
    !           XC-Hamiltonian
    !           The elements of the hamiltonian are stored in the
    !           linear array ham_xc_arr
    ! Input-Parameter:
    !   loop               number of SCF-cycle (within current SCF run)
    !   vpnm (optional)    if called from master, this is the VPN
    !
    use timer_module
    use time_module
    use occupation_module, only: get_n_elec
    use output_module, only: output_grid
    use xc_cntrl, only: is_on, whatis, xc_so_spatial,&
         & xc_nl_calc, xc_gga_version, xc_sop_version
    use iounitadmin_module, only: output_unit
!!$    use xc_ham_trafo, only: show_ham
    integer(kind=i4_kind),optional :: loop
    !** End of interface *****************************************

    integer(kind=i4_kind) :: alloc_stat

    call xc_allocate()
    if (.not.comm_i_am_master()) then
       if (options_spin_orbit) then
          allocate(ham_xc_arr_real(xc_length),ham_xc_arr_imag(xc_length),stat=alloc_stat)
          ASSERT (alloc_stat==0)
       else
          allocate(ham_xc_arr(xc_length),stat=alloc_stat)
          ASSERT (alloc_stat==0)
       endif
    else ! on the master
       ! ham_xc_arr_old is now a temporary array which is deallocate
       ! directly after mixing again
       if (options_spin_orbit) then
          allocate(ham_xc_arr_old_real(xc_length),ham_xc_arr_old_imag(xc_length),stat=alloc_stat)
          ASSERT (alloc_stat==0)
          ham_xc_arr_old_real = ham_xc_arr_real
          ham_xc_arr_old_imag = ham_xc_arr_imag
       else
          allocate(ham_xc_arr_old(xc_length),stat=alloc_stat)
          ASSERT (alloc_stat==0)
          ham_xc_arr_old = ham_xc_arr
       endif
       matold_initialized = mat_initialized
    endif

    if(options_spin_orbit.and.is_on(xc_so_spatial))then
       allocate(ham_xc_arr(xc_length_vec),stat=alloc_stat)
       ASSERT (alloc_stat==0)
       ham_xc_arr = 0.0_r8_kind
    endif


    call init_timer(timer_grid_orbitals)
    call init_timer(timer_grid_density)
    call init_timer(timer_grid_functionals)
    call init_timer(timer_grid_xcbuild)
    ! FIXME: Initialize vars explicitly!
    call xc_clear()
    if (options_spin_orbit) then
       ham_xc_arr_real=0.0_r8_kind
       ham_xc_arr_imag=0.0_r8_kind
       if (spin_orbit_polarized) then
          spin_vector = 0.0
          s_average = 0.0_r8_kind
       endif
       if (xc_nl_calc) then
          exc_gga_int=0.0_r8_kind
       endif
    else
       ham_xc_arr=0.0_r8_kind
    endif
    charge_int=0.0_r8_kind
    exc_int=0.0_r8_kind

    if ( .not. options_spin_orbit ) then
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       if ( .not. xc_nl_calc ) then
          ! LDA calculation
          !
          call calc_xcks(ham_xc_arr, "orbitals")
       else
          ! GGA calculation
          !
          ASSERT(whatis(xc_GGA_version)==2)
          DPRINT "xch/build_xc: calling GGA ..."
          call calc_xcks_nl(ham_xc_arr)
          DPRINT "xch/build_xc: ... done"
       endif
    else
       ! SPIN ORBIT
       if (spin_orbit_polarized) then
          ! OPEN SHELL
          ASSERT(whatis(xc_sop_version)==2)
          DPRINT 'xch/build_xc: call calc_xcks_sop() ...'
          call calc_xcks_sop()
          DPRINT 'xch/build_xc: ... done'
       else
          ! CLOSED SHELL
          if(is_on(xc_so_spatial))then
             call calc_xcks_space_so(ham_xc_arr_real, ham_xc_arr)
          else
             if(xc_nl_calc)then
                call error_handler("xch/build_xc: for SO-GGA use eg. XC='bp,spatial' key")
             endif
             call calc_xcks_so()
          endif
       endif
    endif

    ! now check if we have pseudo 2d irreps and if necesarry perform an additional
    ! symmetrization
    if (.not.options_spin_orbit) then
       call xc_pseudo2d_sym(ham_xc_arr)
    endif

    ! sum up contributions to xc hamiltonian (in global array ham_xc_arr)
    ! from all workers:
    call xc_sndrcv()

    if(comm_i_am_master()) then
       mat_initialized = .false.
       ! now mixing of the hamiltonian
       if (mixing_discard_init()) then
          ! use the automatic mixing strategy
          if (options_spin_orbit) then
             call mixing_ham(ham_xc_arr_real, ham_xc_arr_old_real, matold_initialized)
             call mixing_ham(ham_xc_arr_imag, ham_xc_arr_old_imag, matold_initialized)
          else
             call mixing_ham(ham_xc_arr, ham_xc_arr_old, matold_initialized)
          endif
       else
          ! use the original hard coded mixing strategy
          if ( loop > 3 ) then
             if (options_spin_orbit) then
                call mixing_ham(ham_xc_arr_real, ham_xc_arr_old_real)
                call mixing_ham(ham_xc_arr_imag, ham_xc_arr_old_imag)
             else
                call mixing_ham(ham_xc_arr, ham_xc_arr_old)
             endif
          endif
       endif
       ! the old XC matrix is not needed anymore. Thus:
       if (options_spin_orbit) then
          deallocate(ham_xc_arr_old_real,ham_xc_arr_old_imag,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation failed in su build_xc')
       else
          deallocate(ham_xc_arr_old,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation failed in su build_xc')
       endif
    end if
    call xc_deallocate()
    if(.not.comm_i_am_master()) then
       ! On the master processor, we need ham_xc_arr to
       ! build up ham_tot
       if (options_spin_orbit) then
          deallocate(ham_xc_arr_real,ham_xc_arr_imag,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation failed in su build_xc')
       else
          deallocate(ham_xc_arr,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('deallocation failed in su build_xc')
       endif
    endif

    if(options_spin_orbit.and.is_on(xc_so_spatial))then
       deallocate(ham_xc_arr,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ('xcm/build_xc: ham_xc_arr dealloc failed')
    endif

    if(comm_i_am_master()) then
       if (output_grid) then
          write(output_unit, &
               '("The numerically integrated charge is   :",F21.14)')&
               charge_int
          write(output_unit, &
               '("The numerically integrated xc-energy is:",F21.14)')&
               exc_int
          if (options_spin_orbit.and.xc_nl_calc) then
             write(output_unit, &
                  '("The numerically integrated GGA xc-energy is:",F21.14)')&
                  exc_gga_int
             write(output_unit, &
                  '("The numerically integrated LDA/GGA xc-energy is:",F21.14)')&
                  exc_gga_int+exc_int
          endif
          if (options_spin_orbit.and.spin_orbit_polarized) then
             !
             ! SPIN ORBIT
             !
             write(output_unit, &
                  '("The numerically integrated spin is:",3F21.14)')&
                  spin_vector
             write(output_unit, &
                  '("The numerically averaged spin is:",F21.14)')&
                  s_average
          endif
       endif
    endif
  contains

    subroutine calc_xcks_space_so(ham_xc_p, ham_xc_v)
      ! a piece of code ...
      ! determine Kohn-Sham-Matrix without spin orbit
      ! and then transform it to spin-orbit
      use error_module
      use xc_ham_trafo, only: xc_so_trafo=>so_trafo
      use xc_cntrl,     only: whatis, xc_gga_version, xc_nl_calc
      implicit none
      real(r8_kind), intent(inout) :: ham_xc_p(:), ham_xc_v(:)
      ! *** end of interface ***

      if(.not.xc_nl_calc)then
         call calc_xcks(ham_xc_v, "spinors")
      else
         if(whatis(xc_gga_version).eq.2)then
            DPRINT 'xch/calc_xcks_space_so: call calc_xcks_nl(ham_xc_v,spor_dens=.true.)'
            DPRINT 'xch/calc_xcks_space_so: ( size(ham_xc_v)=',size(ham_xc_v),')'
            call calc_xcks_nl(ham_xc_v, spor_dens=.true.)
         endif
         if ( whatis(xc_gga_version) .eq. 1 ) then
            ABORT("GGAv1 and Spatial do not work together")
         endif
      endif

      !
      ! Transform orbital rep "ham_xc_v" of XC into spinor rep "ham_xc_p".
      ! Both are real, so there is no storage for imaginary parts.
      !
      call xc_so_trafo(ham_xc_p, ham_xc_v)
    end subroutine calc_xcks_space_so

    subroutine calc_xcks_nl(ham_xc_arr,spor_dens)
      ! a piece of code ...
      use error_module
      use f77_blas, only: dgemm
      use xc_cntrl,            only: xc_denscalc_version                       &
                                   , xc_mgga                                   &
                                   , xc_so_orb_to_sporb
      use xc_func, only: xc_functionals
      use density_calc_module, only: density_calc_nl, density_calc_nl_v2
      use grid_module, only: more_grid, grid_loop_setup
      implicit none
      real(r8_kind),intent(inout)   :: ham_xc_arr(:) !(xc_length)
      logical, intent(in), optional :: spor_dens
      ! *** end of interface ***

      !
      ! Calculate matrix elements of Vxc
      !
      ! <i|Vxc|j> = SUM_GRID[ (1/2)F_i*A_SCAL*F_j + F_i*B_VECT*GRAD(F_j) ] +
      !           + SAME(i<->j)
      ! A_SCAL =def= df/drho                 * weights
      ! B_VECT =def= df/dgamma * 2 GRAD(rho) * weights
      !

      logical          :: spin_orbit_density
      integer(I4_kind) :: memstat

      real(r8_kind),dimension(vec_length,    ispin) :: A_SCAL,C_SCAL
      real(r8_kind),dimension(vec_length,X:Z,ispin) :: B_VECT,&
           & grarho

      real(r8_kind), allocatable :: orbs_help(:,:), ham(:,:)

      integer(i4_kind) :: irr, xyz, vla, S, counter, j, k, m
#ifdef FPP_NODGEMM
      integer(i4_kind) :: i
#endif

      ! These were global module vars, now made local subroutine vars:
      real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)

      ! overwrite global variables:
      integer(kind=i4_kind)         :: n_irrep
      integer(kind=i4_kind),pointer :: dims(:),partners(:)

      integer(kind=i4_kind)         :: maxdim ! maxval(dims(:))

      DPRINT 'xch/calc_xcks_nl: entered'

      if(present(spor_dens))then
         spin_orbit_density = spor_dens
      else
         spin_orbit_density = .false.
      endif

      n_irrep  =  n_vec_irrep
      dims     => vec_dims
      partners => vec_partners

      ! dimension for intermediate arrays, max of irrep dimensions:
      maxdim = maxval(dims)

      ! allocate space for intermediate arrays
      allocate(orbs_help(vec_length, maxdim), ham(maxdim, maxdim), stat=memstat)
      ASSERT(memstat==0)

      call grid_loop_setup() !static pre-distribution of jobs
      ! loop over gridpoints:
      do while( more_grid(vec_length, grdpts, grdwts) ) ! fetching part of the grid
         ! more_grid() will return
         !   grdpts => new batch of coordinates
         !   grdwts => corresponding weights
         ! not more than "vec_length" long
         ! it will return only false if the proc cannot steal anymore

         ! this may be, in general, below vec_length, particularly in last
         ! iteration:
         vla = size(grdpts,1)
         DPRINT 'xch/calc_xcks_nl: got vla=',vla

         ! calculate orbitals:
         call start_timer(timer_grid_orbitals)

         if(.not.spin_orbit_density)then
            ! calculate only spatial orbitals:
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_ob, grads=orbs_grads)
         else if(is_on(xc_so_orb_to_sporb))then
            call error("xch/calc_xcks_nl: orb->sporb doesnt work yet")
!!$            ! calculte both at one time:
!!$            call orbital_calculate(grdpts(1:vla,1:3),vla,&
!!$                 & orbs_ob=orbs_ob,&
!!$                 & orbs_spinor_ob=orbs_spinor_ob)
         else
            ! calculate both independently:
            DPRINT 'xch/calc_xcks_nl: call orbital_calculate(orbitals and grads)...'
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_ob, grads=orbs_grads)
            DPRINT 'xch/calc_xcks_nl: call orbital_calculate(spinors and grads)...'
            call orbital_calculate(grdpts(1:vla,1:3),vla,&
                 & orbs_spinor_ob=orbs_spinor_ob, spinor_grads=spinor_grads)
            DPRINT 'xch/calc_xcks_nl: ... done'
         endif

         call stop_timer(timer_grid_orbitals)

         ! calculate density:
         call start_timer(timer_grid_density)

         if(.not.spin_orbit_density)then
           if (.not.is_on(xc_mgga)) then
             call density_calc_nl(vla,rho,gamma,grarho,orbs_ob,orbs_grads)
           else
             call density_calc_nl(vla,rho,gamma,grarho,orbs_ob,orbs_grads,tau=tau)
           endif
         else
            if(whatis(xc_denscalc_version).eq.2)then
               DPRINT 'xch/calc_xcks_nl: call density_calc_nl_v2(spinors and grads)...(dgemm)'
               call density_calc_nl_v2(vla,rho,gamma,grarho,&
                 & orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads)
            else
               DPRINT 'xch/calc_xcks_nl: call density_calc_nl_v1(spinors and grads)...'
               call density_calc_nl(vla,rho,gamma,grarho,&
                 & orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads)
            endif
            DPRINT 'xch/calc_xcks_nl: ... done'
         endif
         ! performing charge integration over the grid
         charge_int=charge_int+sum(rho(1:vla,:)*spread&
              (grdwts(1:vla),2,ispin))
         call stop_timer(timer_grid_density)

         ! values of the functional and its derivatives:
         call start_timer(timer_grid_functionals)

         ! code separated into a sub:
         if (.not.is_on(xc_mgga)) then
           call xc_functionals(vla,ispin,rho,fxc,dfdrho,gamma,dfdgrarho)
         else
           call xc_functionals(vla,ispin,rho,fxc,dfdrho,gamma,dfdgrarho,tau=tau,ft=dfdtau)
         end if

         ! performing exc integration over the grid
         exc_int=exc_int+sum(fxc(1:vla)*grdwts(1:vla))
         call stop_timer(timer_grid_functionals)

         DPRINT 'xch/calc_xcks_nl: INTEGRATE ...'
         call start_timer(timer_grid_xcbuild)

         do S=1,ispin
            A_SCAL(:vla,S) = HALF * dfdrho(:vla,S) * grdwts(:vla)
         end do
         if ( ispin == 1 ) then
            do xyz=X,Z
               B_VECT(:vla,xyz,UP) = grdwts(:vla) *&
                    & TWO * dfdgrarho(:vla,UPUP) * grarho(:vla,xyz,UP)
            enddo
         else
            do xyz=X,Z
               B_VECT(:vla,xyz,UP) = grdwts(:vla) * &
                    & (TWO * dfdgrarho(:vla,UPUP) * grarho(:vla,xyz,UP)&
                    &      + dfdgrarho(:vla,UPDN) * grarho(:vla,xyz,DN))
               B_VECT(:vla,xyz,DN) = grdwts(:vla) * &
                    & (TWO * dfdgrarho(:vla,DNDN) * grarho(:vla,xyz,DN)&
                    &      + dfdgrarho(:vla,UPDN) * grarho(:vla,xyz,UP))
            enddo
         endif
         !--------------------------------------------------------------------------------------------+
         ! calculation of non scalar factor of tau-term of hamiltonian                                |
         !--------------------------------------------------------------------------------------------+
         if (is_on(xc_mgga)) then
           do S=1,ispin
              C_SCAL(:vla,S) = dfdtau(:vla,S) * grdwts(:vla)
           end do
         endif

         do s = 1, ispin
            ! Index into ham_xc_arr(:), alpha/beta at odd/even positions,
            ! irrep blocks one after another:
            counter = s

            do irr = 1, n_irrep
               associate (p => orbs_ob(irr) % o, q => orbs_grads(irr) % o)

                 !
                 ! NOTE: we need to  clear accumulator before sum over
                 ! partner  index m.  This can  be easily  achieved by
                 ! plain simple
                 !
                 !   ham(1:dims(irr), 1:dims(irr)) = 0.0
                 !
                 ! Profiler tools indicate  that such an "inessential"
                 ! clearing   operation   consumes  unnecessary   many
                 ! cycles. So instead, we introduce a branch depending
                 ! on  the partner  index m  ==  1 so  that the  first
                 ! iteration is slightly different:
                 !
               do m = 1, partners(irr)
                  do j = 1, dims(irr)
                     orbs_help(:vla, j) = A_SCAL(:vla, s) * p(:vla, j, m) + &
                          B_VECT(:vla, X, s) * q(:vla, X, j, m) + &
                          B_VECT(:vla, Y, s) * q(:vla, Y, j, m) + &
                          B_VECT(:vla, Z, s) * q(:vla, Z, j, m)
                  enddo
#ifdef FPP_NODGEMM
                do i = 1, dims(irr)
                   do j = 1, dims(irr)
                      acc = dot_product (p(:vla, i, m), orbs_help(:vla, j))
                      if (m == 1) then
                         ham(i, j) = acc
                      else
                         ham(i, j) = ham(i, j) + acc
                      endif
                   enddo
                enddo

                !-------------------------------------------------------------------------------------+
                ! Non-BLAS addition of tau-term to hamiltonian                                        |
                ! dot products calculated first:                                                      |
                !                                                                                     |
                !          H_ij = H_ij + 0.5 * grdwt * dF_dtau * (gradphi^T*gradphi)                  |
                ! STILL NOT SURE IF FACTOR 0.5 IS CORRECT!!!!!!!!!!!!!!!!!!!                          |
                !-------------------------------------------------------------------------------------+
                if (is_on(xc_mgga)) then
                   do i=1,dims(irr)
                      do j=1,dims(irr)
                         ham(i,j) = ham(i,j) + SUM(C_SCAL(:vla,S) * SUM(q(:vla,X:Z,i,m) * q(:vla,X:Z,j,m),DIM=2))*QUARTER
                      enddo
                   enddo
                endif
#else
                block
                   real (r8_kind) :: beta

                   !
                   ! At least  the linux man page says  that with beta
                   ! == 0  the initial value of  accumulator ham(:, :)
                   ! does not matter:
                   !
                   !  "When BETA is supplied as zero then C need not
                   !   be set on input."
                   !
                   if (m == 1) then
                      beta = 0.0
                   else
                      beta = 1.0
                   endif

                   call dgemm ('t', 'n', dims(irr), dims(irr), vla, &
                        ONE, orbs_help(:, :), size (orbs_help, 1), &
                        p(:, :, m), size (p, 1), &
                        beta, ham(:, :), size (ham, 1))
                end block

                if (is_on(xc_mgga)) then
                   !----------------------------------------------------------------------------------+
                   ! BLAS calculation of tau term in the hamiltonian                                  |
                   ! X,Y and Z components added successively:                                         |
                   !                                                                                  |
                   !      H_ij = H_ij + (0.5 * grdwt * dF_dtau * gradphi(X)) * gradphi(X)             |
                   !                  + (0.5 * grdwt * dF_dtau * gradphi(Y)) * gradphi(Y)             |
                   !                  + (0.5 * grdwt * dF_dtau * gradphi(Z)) * gradphi(Z)             |
                   ! STILL NOT SURE IF FACTOR 0.5 IS CORRECT!!!!!!!!!!!!!!!!!!!                       |
                   !----------------------------------------------------------------------------------+
                   do xyz=X,Z
                      !-------------------------------------------------------------------------------+
                      ! scalar factor * gradphi                                                       |
                      !-------------------------------------------------------------------------------+
                      do j=1,dims(irr)
                         orbs_help(:vla,j) = C_SCAL(:vla,S) * q(:vla,xyz,j,m)
                      enddo
                      call dgemm('t', 'n', dims(irr), dims(irr), vla, &
                           QUARTER , orbs_help(:,:)  , size(orbs_help,1)      ,&
                                     q(1:1,xyz,1,m) , size(q,1) * size(q,2) ,&
                           ONE     , ham(:,:)       , size(ham,1))
                   end do
                endif
#endif
               end do! loop over partner

               ! remap stuff to ham_xc_arr
               do j=1,dims(irr)
                  do k=1,j
                     ham_xc_arr(counter)=ham_xc_arr(counter)+(ham(j,k)+ham(k,j))/&
                          partners(irr)
                     counter=counter+ispin
                  end do
               end do
               end associate    ! p => ..., q => ...
            end do              ! irr = 1, n_irrep
         end do                 ! s = 1, ispin
         DPRINT 'xch/calc_xcks_nl: ...done'

         ! FIXME: Initialize vars explicitly!
         call xc_clear()
         call stop_timer(timer_grid_xcbuild)
      end do ! while ( more_grid(...) )

      deallocate(orbs_help, ham, stat=memstat)
      ASSERT(memstat==0)
    end subroutine calc_xcks_nl

    subroutine calc_xcks_so()
      ! a pice of code ...
      use xc_cntrl
      use xc_func, only: xc_functionals
      use density_calc_module, only: density_calc
      use grid_module, only: more_grid, grid_loop_setup
      implicit none
      ! *** end of interface ***

      integer(kind=i4_kind) :: i,j,k,m,vla
      integer(kind=i4_kind) :: i_grid,counter
      real(kind=r8_kind),pointer ::&
           & p_up_real(:,:,:),p_up_imag(:,:,:),&
           & p_dn_real(:,:,:),p_dn_imag(:,:,:) ! intermediate pointers
      real(kind=r8_kind) :: help_xc_real,help_xc_imag
      ! These were global module vars, now made local subroutine vars:
      real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)

!!$      call error_handler("xch/calc_xcks_so: outdated")
!!$
      DPRINT 'xch/calc_xcks_so: entered'
      call grid_loop_setup() !static pre-distribution of jobs
      ! loop over gridpoints:
      do while( more_grid(vec_length, grdpts, grdwts) ) ! fetching part of the grid
         ! more_grid() will return
         !   grdpts => new batch of coordinates
         !   grdwts => corresponding weights
         ! not more than "vec_length" long
         ! it will return only false if the proc cannot steal anymore

         ! this may be, in general, below vec_length, particularly in last
         ! iteration:
         vla = size(grdpts,1)

         call start_timer(timer_grid_orbitals)
         call orbital_calculate(grdpts(1:vla,1:3),&
              vla,orbs_spinor_ob=orbs_spinor_ob)
         call stop_timer(timer_grid_orbitals)
         call start_timer(timer_grid_density)

         call density_calc(vla,rho,orbs_spinor_ob=orbs_spinor_ob)

         call stop_timer(timer_grid_density)
         call start_timer(timer_grid_functionals)

         call xc_functionals(vla,ispin,rho,fxc,dfdrho)

         call stop_timer(timer_grid_functionals)
         ! performing integration over the grid
         charge_int=charge_int+sum(rho(1:vla,1)*grdwts(1:vla))
         exc_int=exc_int+sum(fxc(1:vla)*grdwts(1:vla))
         call start_timer(timer_grid_xcbuild)
         counter=1
         do i=1,n_irrep
            p_up_real=>orbs_spinor_ob(i)%o(:,:,:,1,1) !spinor(1)%o
            p_up_imag=>orbs_spinor_ob(i)%o(:,:,:,2,1) !spinor(1)%o_imag
            p_dn_real=>orbs_spinor_ob(i)%o(:,:,:,1,2) !spinor(2)%o
            p_dn_imag=>orbs_spinor_ob(i)%o(:,:,:,2,2) !spinor(2)%o_imag
            do j=1,dims(i)
               do k=1,j
                  help_xc_real=0.0_r8_kind
                  help_xc_imag=0.0_r8_kind
                  help_arr_real(1:vla)=&
                       ( p_up_real(1:vla,j,1)*p_up_real(1:vla,k,1)+&
                       p_up_imag(1:vla,j,1)*p_up_imag(1:vla,k,1)+&
                       p_dn_real(1:vla,j,1)*p_dn_real(1:vla,k,1)+&
                       p_dn_imag(1:vla,j,1)*p_dn_imag(1:vla,k,1))*&
                       grdwts(1:vla)
                  help_arr_imag(1:vla)=&
                       ( p_up_real(1:vla,j,1)*p_up_imag(1:vla,k,1)-&
                       p_up_imag(1:vla,j,1)*p_up_real(1:vla,k,1)+&
                       p_dn_real(1:vla,j,1)*p_dn_imag(1:vla,k,1)-&
                       p_dn_imag(1:vla,j,1)*p_dn_real(1:vla,k,1))*&
                       grdwts(1:vla)
                  !call buggy("build_xc: now loop over partners")
                  do m=2,partners(i)
                     help_arr_real(1:vla)=help_arr_real(1:vla)+&
                          ( p_up_real(1:vla,j,m)*p_up_real(1:vla,k,m)+&
                          p_up_imag(1:vla,j,m)*p_up_imag(1:vla,k,m)+&
                          p_dn_real(1:vla,j,m)*p_dn_real(1:vla,k,m)+&
                          p_dn_imag(1:vla,j,m)*p_dn_imag(1:vla,k,m))*&
                          grdwts(1:vla)
                     help_arr_imag(1:vla)=help_arr_imag(1:vla)+&
                          ( p_up_real(1:vla,j,m)*p_up_imag(1:vla,k,m)-&
                          p_up_imag(1:vla,j,m)*p_up_real(1:vla,k,m)+&
                          p_dn_real(1:vla,j,m)*p_dn_imag(1:vla,k,m)-&
                          p_dn_imag(1:vla,j,m)*p_dn_real(1:vla,k,m))*&
                          grdwts(1:vla)
                  enddo! loop over partners
                  ! building of xc_hamiltonian
                  do i_grid=1,vla
                     help_xc_real=help_xc_real+&
                          dfdrho(i_grid,1)*help_arr_real&
                          (i_grid)
                     help_xc_imag=help_xc_imag+&
                          dfdrho(i_grid,1)*help_arr_imag&
                          (i_grid)
                  enddo
                  ! building of xc_hamiltonian finished
                  ham_xc_arr_real(counter)=ham_xc_arr_real(counter)+help_xc_real/partners(i)
                  !ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)+help_xc_imag/partners(i)
                  ! ***
                  ! Note: Here we are calculating the matrix elements <k|Vxc|j> not <j|Vxc|k>
                  ! ***
                  ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)-help_xc_imag/partners(i)
                  counter=counter+1
               enddo!loop over nu
            enddo!loop over mu
         enddo! loop over irreps
         ! FIXME: Initialize vars explicitly!
         call xc_clear()
         call stop_timer(timer_grid_xcbuild)
      enddo
      DPRINT 'xch/calc_xcks_so: exit'
    end subroutine calc_xcks_so

    subroutine calc_xcks_sop()
      ! a piece of the code ...
      use error_module
      use xc_cntrl, only:&
           & is_on, xc_xalpha, xc_vwn, cutoff=>xc_sdens_cutoff
!!$      use options_module, only: options_debug_key ! <- backdoor
      use f77_blas, only: dgemm
      use xc_func, only: xc_functionals
      use density_calc_module, only: density_calc
      use grid_module, only: more_grid, grid_loop_setup
      implicit none
      ! *** end of interface ***

      ! These were global module vars, now made local subroutine vars:
      real(kind=r8_kind),pointer  :: grdpts(:,:),grdwts(:)

      integer(kind=i4_kind) :: vl,i,j,m,irr,memstat
      integer(kind=i4_kind) :: counter
      real(kind=r8_kind),pointer ::&
           & p_up_real(:,:,:),p_up_imag(:,:,:),&
           & p_dn_real(:,:,:),p_dn_imag(:,:,:) ! intermediate pointers
      real(kind=r8_kind), allocatable :: ptmp(:,:,:,:) ! (vl,dim(irr),UP:DN,RE:IM)
      real(kind=r8_kind), allocatable :: ham(:,:,:)    ! (dim(irr),dim(irr),RE:IM)

      call grid_loop_setup() !static pre-distribution of jobs
      ! loop over gridpoints:
      do while( more_grid(vec_length, grdpts, grdwts) ) ! fetching part of the grid
         ! more_gridd() will return
         !   grdpts => new batch of coordinates
         !   grdwts => corresponding weights
         ! not more than "vec_length" long
         ! it will return only false if the proc cannot steal anymore

         ! this may be, in general, below vec_length, particularly in last
         ! iteration:
         vl = size(grdpts,1)

         call start_timer(timer_grid_orbitals)
         call orbital_calculate(grdpts(1:vl,1:3),&
              vl,orbs_spinor_ob=orbs_spinor_ob)
         call stop_timer(timer_grid_orbitals)

         call start_timer(timer_grid_density)
         call density_calc(vl,rho,orbs_spinor_ob=orbs_spinor_ob,s_dens=s_dens)
         call stop_timer(timer_grid_density)

         call start_timer(timer_grid_functionals)
         !
         ! calculate absolute value of spin density |s|
         !
         s_abs(1:vl) = sqrt(s_dens(1:vl,1)**2+&
              &             s_dens(1:vl,2)**2+&
              &             s_dens(1:vl,3)**2)
         !
         ! determine "spin up" and "spin down" densities
         ! following the transformations:
         ! rho_up   = (rho + |s|)/2
         ! rho_down = (rho - |s|)/2
         rho_ud(1:vl,1) =  0.5_r8_kind*(rho(1:vl,1) + s_abs(1:vl)) ! spin up
         rho_ud(1:vl,2) =  0.5_r8_kind*(rho(1:vl,1) - s_abs(1:vl)) ! spin down

         call xc_functionals(vl,2,rho_ud,fxc,dfdrho)

         ! now calculate derivatives by n and s
         ! following the formulas:
         ! df/drho = 1/2*(df/d(rho_up) + df/d(rho_down))
         ! df/|s|  = 1/2*(df/d(rho_up) - df/d(rho_down))
         ! (and multiply them with weights)
         dfdn(1:vl) = 0.5_r8_kind*(dfdrho(1:vl,1) + dfdrho(1:vl,2)) * grdwts(:vl)
         dfds(1:vl) = 0.5_r8_kind*(dfdrho(1:vl,1) - dfdrho(1:vl,2)) * grdwts(:vl)

         ! calculate total spin
         spin_vector(1) = spin_vector(1) + sum(abs(s_dens(1:vl,1))*grdwts(1:vl))
         spin_vector(2) = spin_vector(2) + sum(abs(s_dens(1:vl,2))*grdwts(1:vl))
         spin_vector(3) = spin_vector(3) + sum(abs(s_dens(1:vl,3))*grdwts(1:vl))

         ! normalize spin density:
         where( s_abs(1:vl) .lt. cutoff )
            s_dens(1:vl,1) = zero ! direction
            s_dens(1:vl,2) = zero ! is not
            s_dens(1:vl,3) = zero ! defined
         elsewhere
            s_dens(1:vl,1) = s_dens(1:vl,1)/s_abs(1:vl)
            s_dens(1:vl,2) = s_dens(1:vl,2)/s_abs(1:vl)
            s_dens(1:vl,3) = s_dens(1:vl,3)/s_abs(1:vl)
         end where

         ! performing integration over the grid
         charge_int=charge_int+sum(rho(1:vl,1)*grdwts(1:vl))
         exc_int=exc_int+sum(fxc(1:vl)*grdwts(1:vl))
         s_average = s_average + sum(s_abs(1:vl)*grdwts(1:vl))

         ! preevaluate df/ds = df/d|s| * s/|s|:
         s_dens(:vl,1) = dfds(:vl) * s_dens(:vl,1)
         s_dens(:vl,2) = dfds(:vl) * s_dens(:vl,2)
         s_dens(:vl,3) = dfds(:vl) * s_dens(:vl,3)

         call stop_timer(timer_grid_functionals)

         call start_timer(timer_grid_xcbuild)
         counter=1
         do irr=1,n_irrep
            p_up_real=>orbs_spinor_ob(irr)%o(:,:,:,1,1) !spinor(1)%o
            p_up_imag=>orbs_spinor_ob(irr)%o(:,:,:,2,1) !spinor(1)%o_imag
            p_dn_real=>orbs_spinor_ob(irr)%o(:,:,:,1,2) !spinor(2)%o
            p_dn_imag=>orbs_spinor_ob(irr)%o(:,:,:,2,2) !spinor(2)%o_imag

            allocate( ptmp(vl,dims(irr),UP:DN,RE:IM), ham(dims(irr),dims(irr),RE:IM),&
                 & stat=memstat)
            if(memstat.ne.0) call error("xch/calc_xcks_sop: alloc ptmp failed")

            ham = zero
            do m=1,partners(irr)
               ! multiplying with a + b*sigma, a = df/dn, b = df/d|s| * s/|s|
               ! | a  + bz  bx - iby |
               ! | bx + iby a  - bz  |
               do i=1,dims(irr)
                  ptmp(:vl,i,UP,RE) =&
                       & (dfdn(:vl) + s_dens(:vl,Z)) * p_up_real(:vl,i,m)&
                       &            + s_dens(:vl,X)  * p_dn_real(:vl,i,m)&
                       &            + s_dens(:vl,Y)  * p_dn_imag(:vl,i,m)

                  ptmp(:vl,i,UP,IM) =&
                       & (dfdn(:vl) + s_dens(:vl,Z)) * p_up_imag(:vl,i,m)&
                       &            + s_dens(:vl,X)  * p_dn_imag(:vl,i,m)&
                       &            - s_dens(:vl,Y)  * p_dn_real(:vl,i,m)

                  ptmp(:vl,i,DN,RE) =&
                       & (dfdn(:vl) - s_dens(:vl,Z)) * p_dn_real(:vl,i,m)&
                       &            + s_dens(:vl,X)  * p_up_real(:vl,i,m)&
                       &            - s_dens(:vl,Y)  * p_up_imag(:vl,i,m)

                  ptmp(:vl,i,DN,IM) =&
                       & (dfdn(:vl) - s_dens(:vl,Z)) * p_dn_imag(:vl,i,m)&
                       &            + s_dens(:vl,X)  * p_up_imag(:vl,i,m)&
                       &            + s_dens(:vl,Y)  * p_up_real(:vl,i,m)
               enddo
#ifdef FPP_NOBLAS
               ! NO BLAS VERSION:
               do i=1,dims(irr)
                  do j=1,i
                     do k=1,vl
                        ham(j,i,RE) = ham(j,i,RE)&
                             & + p_up_real(k,j,m) * ptmp(k,i,UP,RE)&
                             & + p_up_imag(k,j,m) * ptmp(k,i,UP,IM)&
                             & + p_dn_real(k,j,m) * ptmp(k,i,DN,RE)&
                             & + p_dn_imag(k,j,m) * ptmp(k,i,DN,IM)
                        ham(j,i,IM) = ham(j,i,IM)&
                             & + p_up_real(k,j,m) * ptmp(k,i,UP,IM)&
                             & - p_up_imag(k,j,m) * ptmp(k,i,UP,RE)&
                             & + p_dn_real(k,j,m) * ptmp(k,i,DN,IM)&
                             & - p_dn_imag(k,j,m) * ptmp(k,i,DN,RE)
                     enddo
                  enddo
               enddo
#else
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_up_real(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,UP,RE) , vl,&
                    &  one, ham(:,:,RE)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_up_imag(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,UP,IM) , vl,&
                    &  one, ham(:,:,RE)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_dn_real(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,DN,RE) , vl,&
                    &  one, ham(:,:,RE)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_dn_imag(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,DN,IM) , vl,&
                    &  one, ham(:,:,RE)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_up_real(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,UP,IM) , vl,&
                    &  one, ham(:,:,IM)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & -one, p_up_imag(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,UP,RE) , vl,&
                    &  one, ham(:,:,IM)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & +one, p_dn_real(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,DN,IM) , vl,&
                    &  one, ham(:,:,IM)     , dims(irr))
               call dgemm('t','n',dims(irr),dims(irr),vl,&
                    & -one, p_dn_imag(:,:,m), size(p_up_real,1),&
                    &       ptmp(:,:,DN,RE) , vl,&
                    &  one, ham(:,:,IM)     , dims(irr))
#endif
            enddo! partners

            do i=1,dims(irr)
               do j=1,i
                  ham_xc_arr_real(counter) = ham_xc_arr_real(counter)&
                       & + ham(j,i,RE)/partners(irr)
                  ham_xc_arr_imag(counter) = ham_xc_arr_imag(counter)&
                       & + ham(j,i,IM)/partners(irr)
                  counter = counter + 1
               enddo
            enddo

            deallocate(ptmp,ham,stat=memstat)
            if(memstat.ne.0) call error("xch/calc_xcks_sop: dealloc ptmp failed")
         enddo! irreps

         ! FIXME: Initialize vars explicitly!
         call xc_clear()
         call stop_timer(timer_grid_xcbuild)
      enddo
    end subroutine calc_xcks_sop

    subroutine calc_xcks_nl_so
      call error_handler("xch/calc_xcks_nl_so: need update")
!!$      ! determine Kohn-Sham Matrix including Spin-Orbit and nonlocal potentials
!!$      exc_gga_int = 0.0_r8_kind
!!$      if (spin_orbit_polarized) then
!!$         !
!!$         ! OPEN SHELL
!!$         !
!!$         do
!!$            last=give_grid(vec_length,grdpts,grdwts)
!!$            vla=size(grdpts,1)
!!$            call start_timer(timer_grid_orbitals)
!!$            call orbital_calculate(grdpts(1:vla,1:3),&
!!$                 vla,orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
!!$            call stop_timer(timer_grid_orbitals)
!!$            call start_timer(timer_grid_density)
!!$            call density_calc_nl(vla,rho,gamma,grarhox,grarhoy,&
!!$                 grarhoz,orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads,&
!!$                 s_abs=s_abs,s_dens=s_dens,gras_densx=gras_densx,gras_densy=gras_densy,&
!!$                 gras_densz=gras_densz)
!!$            ! determine "spin up" and "spin down" densities
!!$            ! following the transformations:
!!$            ! rho_up   = (rho + |s|)/2
!!$            ! rho_down = (rho - |s|)/2
!!$
!!$            rho_ud(1:vla,1) =  0.5_r8_kind*(rho(1:vla,1) + s_abs(1:vla)) ! spin up
!!$            rho_ud(1:vla,2) =  0.5_r8_kind*(rho(1:vla,1) - s_abs(1:vla)) ! spin down
!!$            where (abs(rho_ud).lt.1.0e-40)
!!$               rho_ud = 0.0_r8_kind
!!$            end where
!!$
!!$            ! determine gradients of "spin up" and "spin down" densities
!!$            grarho_udx(1:vla,1) =  0.5_r8_kind*(grarhox(1:vla,1) + gras_densx(1:vla)) ! spin up
!!$            grarho_udx(1:vla,2) =  0.5_r8_kind*(grarhox(1:vla,1) - gras_densx(1:vla)) ! spin down
!!$            grarho_udy(1:vla,1) =  0.5_r8_kind*(grarhoy(1:vla,1) + gras_densy(1:vla)) ! spin up
!!$            grarho_udy(1:vla,2) =  0.5_r8_kind*(grarhoy(1:vla,1) - gras_densy(1:vla)) ! spin down
!!$            grarho_udz(1:vla,1) =  0.5_r8_kind*(grarhoz(1:vla,1) + gras_densz(1:vla)) ! spin up
!!$            grarho_udz(1:vla,2) =  0.5_r8_kind*(grarhoz(1:vla,1) - gras_densz(1:vla)) ! spin down
!!$
!!$            gamma(1:vla,1) = grarho_udx(1:vla,1)*grarho_udx(1:vla,1)+&
!!$                 grarho_udy(1:vla,1)*grarho_udy(1:vla,1)+&
!!$                 grarho_udz(1:vla,1)*grarho_udz(1:vla,1)
!!$            gamma(1:vla,2) = grarho_udx(1:vla,2)*grarho_udx(1:vla,2)+&
!!$                 grarho_udy(1:vla,2)*grarho_udy(1:vla,2)+&
!!$                 grarho_udz(1:vla,2)*grarho_udz(1:vla,2)
!!$            gamma(1:vla,3) = grarho_udx(1:vla,1)*grarho_udx(1:vla,2)+&
!!$                 grarho_udy(1:vla,1)*grarho_udy(1:vla,2)+&
!!$                 grarho_udz(1:vla,1)*grarho_udz(1:vla,2)
!!$
!!$            where (abs(gamma).lt.1.0e-40)
!!$               gamma = 0.0_r8_kind
!!$            end where
!!$
!!$            call stop_timer(timer_grid_density)
!!$            call start_timer(timer_grid_functionals)
!!$            if (vwn) then
!!$               call vwn_calc(rho_ud,dfdrho,2,fxc,vla)
!!$            end if
!!$            if (xalpha) call xalpha_calc(rho_ud,dfdrho,2,fxc,vla)
!!$            if(perdew)  call perdew_calc(rho_ud,gamma,dfdrho_dummy,2,fxc_gga,dfdgrarho,&
!!$                 vla)
!!$            if(becke88)  call becke88_calc&
!!$                 (rho_ud,gamma,dfdrho_dummy,2,fxc_gga,dfdgrarho,vla)
!!$            if(perdewwang91c)  call pw91c_calc&
!!$                 (rho_ud,gamma,dfdrho_dummy,2,fxc_gga,dfdgrarho,vla,&
!!$                 fxc_lda,dfdrho_lda)
!!$            if(perdewwang91x)  call pw91x_calc&
!!$                 (rho_ud,gamma,dfdrho_dummy,2,fxc_gga,dfdgrarho,vla)
!!$            if(baerends94)  call baerends94_calc&
!!$                 (rho_ud,gamma,dfdrho_dummy,2,fxc_gga,dfdgrarho,vla)
!!$            ! determine spin expectation value
!!$            s_exp(1:vla) = s_abs(1:vla)/rho(1:vla,1)
!!$
!!$            ! now calculate derivatives by n and s
!!$            ! following the formulas:
!!$            ! df/drho = 1/2*(df/d(rho_up) + df/d(rho_down))
!!$            ! df/|s|  = 1/2*(df/d(rho_up) - df/d(rho_down))
!!$            dfdn(1:vla) = 0.5_r8_kind*(dfdrho(1:vla,1) + dfdrho(1:vla,2))
!!$            dfds(1:vla) = 0.5_r8_kind*(dfdrho(1:vla,1) - dfdrho(1:vla,2))
!!$            ! calculate total spin
!!$            sx_int = sx_int + sum(s_dens(1:vla,1)*grdwts(1:vla))
!!$            sy_int = sy_int + sum(s_dens(1:vla,2)*grdwts(1:vla))
!!$            sz_int = sz_int + sum(s_dens(1:vla,3)*grdwts(1:vla))
!!$
!!$            ! normalize spin density:
!!$            s_dens(1:vla,1) = s_dens(1:vla,1)/s_abs(1:vla)
!!$            s_dens(1:vla,2) = s_dens(1:vla,2)/s_abs(1:vla)
!!$            s_dens(1:vla,3) = s_dens(1:vla,3)/s_abs(1:vla)
!!$
!!$
!!$            call stop_timer(timer_grid_functionals)
!!$            ! performing integration over the grid
!!$            charge_int=charge_int+sum(rho(1:vla,1)*grdwts(1:vla))
!!$            exc_int=exc_int+sum(fxc(1:vla)*grdwts(1:vla))
!!$            exc_gga_int=exc_gga_int+sum(fxc_gga(1:vla)*grdwts(1:vla))
!!$            s_average = s_average + sum(s_abs(1:vla)*grdwts(1:vla))
!!$            call start_timer(timer_grid_xcbuild)
!!$            counter=1
!!$            do i=1,n_irrep
!!$               p_up_real=>orbs_spinor_ob(i)%spinor(1)%o
!!$               p_up_imag=>orbs_spinor_ob(i)%spinor(1)%o_imag
!!$               p_dn_real=>orbs_spinor_ob(i)%spinor(2)%o
!!$               p_dn_imag=>orbs_spinor_ob(i)%spinor(2)%o_imag
!!$               do j=1,dims(i)
!!$                  do k=1,j
!!$                     help_xc_real=0.0_r8_kind
!!$                     help_xc_imag=0.0_r8_kind
!!$                     !call buggy("build_xc: calculate help_arr_real")
!!$                     ! <chij,up|chik,up>
!!$                     help_arr_uu_real(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_up_real(1:vla,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*p_up_imag(1:vla,k,1))*grdwts(1:vla)
!!$                     help_arr_uu_imag(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*p_up_real(1:vla,k,1))*grdwts(1:vla)
!!$                     ! <chij,up|chik,dn>
!!$                     help_arr_ud_real(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*p_dn_imag(1:vla,k,1))*grdwts(1:vla)
!!$                     help_arr_ud_imag(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*p_dn_real(1:vla,k,1))*grdwts(1:vla)
!!$                     ! <chij,dn|chik,up>
!!$                     help_arr_du_real(1:vla)=&
!!$                          ( p_dn_real(1:vla,j,1)*p_up_real(1:vla,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*p_up_imag(1:vla,k,1))*grdwts(1:vla)
!!$                     help_arr_du_imag(1:vla)=&
!!$                          ( p_dn_real(1:vla,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*p_up_real(1:vla,k,1))*grdwts(1:vla)
!!$                     ! <chij,dn|chik,dn>
!!$                     help_arr_dd_real(1:vla)=&
!!$                          (p_dn_real(1:vla,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*p_dn_imag(1:vla,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arr_dd_imag(1:vla)=&
!!$                          (p_dn_real(1:vla,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*p_dn_real(1:vla,k,1))*&
!!$                          grdwts(1:vla)
!!$                     !call buggy("build_xc: now loop over partners")
!!$                     do m=2,partners(i)
!!$                        ! <chij,up|chik,up>
!!$                        help_arr_uu_real(1:vla)=help_arr_uu_real(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_up_real(1:vla,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*p_up_imag(1:vla,k,m))*grdwts(1:vla)
!!$                        help_arr_uu_imag(1:vla)=help_arr_uu_imag(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*p_up_real(1:vla,k,m))*grdwts(1:vla)
!!$                        ! <chij,up|chik,dn>
!!$                        help_arr_ud_real(1:vla)=help_arr_ud_real(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*p_dn_imag(1:vla,k,m))*grdwts(1:vla)
!!$                        help_arr_ud_imag(1:vla)=help_arr_ud_imag(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*p_dn_real(1:vla,k,m))*grdwts(1:vla)
!!$                        ! <chij,dn|chik,up>
!!$                        help_arr_du_real(1:vla)=help_arr_du_real(1:vla)+&
!!$                             ( p_dn_real(1:vla,j,m)*p_up_real(1:vla,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*p_up_imag(1:vla,k,m))*grdwts(1:vla)
!!$                        help_arr_du_imag(1:vla)=help_arr_du_imag(1:vla)+&
!!$                             ( p_dn_real(1:vla,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*p_up_real(1:vla,k,m))*grdwts(1:vla)
!!$                        ! <chij,dn|chik,dn>
!!$                        help_arr_dd_real(1:vla)=help_arr_dd_real(1:vla)+&
!!$                             (p_dn_real(1:vla,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*p_dn_imag(1:vla,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arr_dd_imag(1:vla)=help_arr_dd_imag(1:vla)+&
!!$                             (p_dn_real(1:vla,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*p_dn_real(1:vla,k,m))*&
!!$                             grdwts(1:vla)
!!$                     enddo! loop over partners
!!$                     !call buggy("build_xc: building of xc_hamiltonian")
!!$                     ! building of xc_hamiltonian
!!$                     !
!!$                     ! df/dn
!!$                     !
!!$                     do i_grid=1,vla
!!$                        help_xc_real=help_xc_real+&
!!$                             dfdn(i_grid)*(help_arr_uu_real&
!!$                             (i_grid) + help_arr_dd_real(i_grid))
!!$                        help_xc_imag=help_xc_imag+&
!!$                             dfdn(i_grid)*(help_arr_uu_imag&
!!$                             (i_grid) + help_arr_dd_imag(i_grid))
!!$                     enddo
!!$                     !
!!$                     ! df/ds * s/|s| * sigma
!!$                     !
!!$                     ! sx * sigmax:
!!$                     do i_grid=1,vla
!!$                        help_xc_real=help_xc_real+&
!!$                             dfds(i_grid)*s_dens(i_grid,1)*(help_arr_ud_real&
!!$                             (i_grid) + help_arr_du_real(i_grid))
!!$                        help_xc_imag=help_xc_imag+&
!!$                             dfds(i_grid)*s_dens(i_grid,1)*(help_arr_ud_imag&
!!$                             (i_grid) + help_arr_du_imag(i_grid))
!!$                     enddo
!!$                     ! sy * sigmay:
!!$                     do i_grid=1,vla
!!$                        help_xc_real=help_xc_real+&
!!$                             dfds(i_grid)*s_dens(i_grid,2)*(help_arr_ud_imag&
!!$                             (i_grid) - help_arr_du_imag(i_grid))
!!$                        help_xc_imag=help_xc_imag+&
!!$                             dfds(i_grid)*s_dens(i_grid,2)*(help_arr_du_real&
!!$                             (i_grid) - help_arr_ud_real(i_grid))
!!$                     enddo
!!$                     ! sz * sigmaz:
!!$                     do i_grid=1,vla
!!$                        help_xc_real=help_xc_real+&
!!$                             dfds(i_grid)*s_dens(i_grid,3)*(help_arr_uu_real&
!!$                             (i_grid) - help_arr_dd_real(i_grid))
!!$                        help_xc_imag=help_xc_imag+&
!!$                             dfds(i_grid)*s_dens(i_grid,3)*(help_arr_uu_imag&
!!$                             (i_grid) - help_arr_dd_imag(i_grid))
!!$                     enddo
!!$                     ! building of xc_hamiltonian finished
!!$                     ham_xc_arr_real(counter)=ham_xc_arr_real(counter)+help_xc_real/partners(i)
!!$                     !ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)+help_xc_imag/partners(i)
!!$                     ! ***
!!$                     ! Note: Here we are calculating the matrix elements <k|Vxc|j> not <j|Vxc|k>
!!$                     ! ***
!!$                     ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)-help_xc_imag/partners(i)
!!$                     counter=counter+1
!!$                  enddo!loop over nu
!!$               enddo!loop over mu
!!$            enddo! loop over irreps
!!$            call xc_clear()
!!$            call stop_timer(timer_grid_xcbuild)
!!$            if(last) exit
!!$         enddo
!!$      else ! spin_orbit_polarized
!!$         !
!!$         ! CLOSED SHELL
!!$         !
!!$         do
!!$            last=give_grid(vec_length,grdpts,grdwts)
!!$            vla=size(grdpts,1)
!!$            call start_timer(timer_grid_orbitals)
!!$            call orbital_calculate(grdpts(1:vla,1:3),&
!!$                 vla,orbs_spinor_ob=orbs_spinor_ob,spinor_grads=spinor_grads)
!!$            call stop_timer(timer_grid_orbitals)
!!$            call start_timer(timer_grid_density)
!!$            call density_calc_nl(vla,rho,gamma,grarhox,grarhoy,&
!!$                 grarhoz,orbs_spinor_ob=orbs_spinor_ob,orbs_spinor_grads=spinor_grads)
!!$            call stop_timer(timer_grid_density)
!!$            call start_timer(timer_grid_functionals)
!!$            if (vwn) then
!!$               call vwn_calc(rho,dfdrho,ispin,fxc,vla)
!!$            end if
!!$            if (xalpha) call xalpha_calc(rho,dfdrho,ispin,fxc,vla)
!!$            if(perdew)  call perdew_calc(rho,gamma,dfdrho_dummy,ispin,fxc_gga,dfdgrarho,&
!!$                 vla)
!!$            if(becke88)  call becke88_calc&
!!$                 (rho,gamma,dfdrho_dummy,ispin,fxc_gga,dfdgrarho,vla)
!!$            if(perdewwang91c)  call pw91c_calc&
!!$                 (rho,gamma,dfdrho_dummy,ispin,fxc_gga,dfdgrarho,vla,&
!!$                 fxc_lda,dfdrho_lda)
!!$            if(perdewwang91x)  call pw91x_calc&
!!$                 (rho,gamma,dfdrho_dummy,ispin,fxc_gga,dfdgrarho,vla)
!!$            if(baerends94)  call baerends94_calc&
!!$                 (rho,gamma,dfdrho_dummy,ispin,fxc_gga,dfdgrarho,vla)
!!$            call stop_timer(timer_grid_functionals)
!!$            ! performing charge and exc integration over the grid
!!$            charge_int=charge_int+sum(rho(1:vla,:)*spread&
!!$                 (grdwts(1:vla),2,ispin))
!!$            exc_int=exc_int+sum(fxc(1:vla)*grdwts(1:vla))
!!$            exc_gga_int=exc_gga_int+sum(fxc_gga(1:vla)*grdwts(1:vla))
!!$            call start_timer(timer_grid_xcbuild)
!!$            counter=1
!!$            ! calculate 2 * df/d(grad[rho]) * grad[rho]
!!$            help_x_up(1:vla)=2.0_r8_kind*&
!!$                 dfdgrarho(1:vla,1)*&
!!$                 grarhox(1:vla,1)
!!$            help_y_up(1:vla)=2.0_r8_kind*&
!!$                 dfdgrarho(1:vla,1)*&
!!$                 grarhoy(1:vla,1)
!!$            help_z_up(1:vla)=2.0_r8_kind*&
!!$                 dfdgrarho(1:vla,1)*&
!!$                 grarhoz(1:vla,1)
!!$            do i=1,n_irrep
!!$               p_up_real=>orbs_spinor_ob(i)%spinor(1)%o
!!$               p_up_imag=>orbs_spinor_ob(i)%spinor(1)%o_imag
!!$               p_dn_real=>orbs_spinor_ob(i)%spinor(2)%o
!!$               p_dn_imag=>orbs_spinor_ob(i)%spinor(2)%o_imag
!!$               q_up_real=>spinor_grads(i)%spinor(1)%o
!!$               q_up_imag=>spinor_grads(i)%spinor(1)%o_imag
!!$               q_dn_real=>spinor_grads(i)%spinor(2)%o
!!$               q_dn_imag=>spinor_grads(i)%spinor(2)%o_imag
!!$               ! now calculate all matrix element (j,k) for irrep i
!!$               do j=1,dims(i)
!!$                  do k=1,j
!!$                     help_xc_real=0.0_r8_kind
!!$                     help_xc_imag=0.0_r8_kind
!!$                     !call buggy("build_xc: calculate help_arr_real")
!!$                     !print*,"build_xc: calculate help_arr_real"
!!$                     help_arr_real(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_up_real(1:vla,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*p_up_imag(1:vla,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*p_dn_imag(1:vla,k,1))*&
!!$                          grdwts(1:vla)
!!$                     !call buggy("build_xc: calculate help_arr_imag")
!!$                     help_arr_imag(1:vla)=&
!!$                          ( p_up_real(1:vla,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*p_up_real(1:vla,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*p_dn_real(1:vla,k,1))*&
!!$                          grdwts(1:vla)
!!$                     !
!!$                     ! calculate grad(psi(k)*psi(j))
!!$                     !
!!$                     help_arrx_real(1:vla)=&
!!$                          ( q_up_real(1:vla,1,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_up_imag(1:vla,1,j,1)*p_up_imag(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,1,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          q_dn_imag(1:vla,1,j,1)*p_dn_imag(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_real(1:vla,1,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*q_up_imag(1:vla,1,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_real(1:vla,1,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_imag(1:vla,1,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arrx_imag(1:vla)=&
!!$                          ( q_up_real(1:vla,1,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          q_up_imag(1:vla,1,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,1,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          q_dn_imag(1:vla,1,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_imag(1:vla,1,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*q_up_real(1:vla,1,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_imag(1:vla,1,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_real(1:vla,1,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arry_real(1:vla)=&
!!$                          ( q_up_real(1:vla,2,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_up_imag(1:vla,2,j,1)*p_up_imag(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,2,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          q_dn_imag(1:vla,2,j,1)*p_dn_imag(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_real(1:vla,2,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*q_up_imag(1:vla,2,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_real(1:vla,2,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_imag(1:vla,2,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arry_imag(1:vla)=&
!!$                          ( q_up_real(1:vla,2,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          q_up_imag(1:vla,2,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,2,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          q_dn_imag(1:vla,2,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_imag(1:vla,2,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*q_up_real(1:vla,2,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_imag(1:vla,2,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_real(1:vla,2,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arrz_real(1:vla)=&
!!$                          ( q_up_real(1:vla,3,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_up_imag(1:vla,3,j,1)*p_up_imag(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,3,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          q_dn_imag(1:vla,3,j,1)*p_dn_imag(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_real(1:vla,3,k,1)+&
!!$                          p_up_imag(1:vla,j,1)*q_up_imag(1:vla,3,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_real(1:vla,3,k,1)+&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_imag(1:vla,3,k,1))*&
!!$                          grdwts(1:vla)
!!$                     help_arrz_imag(1:vla)=&
!!$                          ( q_up_real(1:vla,3,j,1)*p_up_imag(1:vla,k,1)-&
!!$                          q_up_imag(1:vla,3,j,1)*p_up_real(1:vla,k,1)+&
!!$                          q_dn_real(1:vla,3,j,1)*p_dn_imag(1:vla,k,1)-&
!!$                          q_dn_imag(1:vla,3,j,1)*p_dn_real(1:vla,k,1)+&
!!$                          p_up_real(1:vla,j,1)*q_up_imag(1:vla,3,k,1)-&
!!$                          p_up_imag(1:vla,j,1)*q_up_real(1:vla,3,k,1)+&
!!$                          p_dn_real(1:vla,j,1)*q_dn_imag(1:vla,3,k,1)-&
!!$                          p_dn_imag(1:vla,j,1)*q_dn_real(1:vla,3,k,1))*&
!!$                          grdwts(1:vla)
!!$
!!$                     !call buggy("build_xc: now loop over partners")
!!$                     do m=2,partners(i)
!!$                        help_arr_real(1:vla)=help_arr_real(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_up_real(1:vla,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*p_up_imag(1:vla,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*p_dn_imag(1:vla,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arr_imag(1:vla)=help_arr_imag(1:vla)+&
!!$                             ( p_up_real(1:vla,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*p_up_real(1:vla,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*p_dn_real(1:vla,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arrx_real(1:vla)=help_arrx_real(1:vla) +&
!!$                             ( q_up_real(1:vla,1,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_up_imag(1:vla,1,j,m)*p_up_imag(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,1,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             q_dn_imag(1:vla,1,j,m)*p_dn_imag(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_real(1:vla,1,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*q_up_imag(1:vla,1,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_real(1:vla,1,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_imag(1:vla,1,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arrx_imag(1:vla)=help_arry_imag(1:vla) +&
!!$                             ( q_up_real(1:vla,1,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             q_up_imag(1:vla,1,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,1,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             q_dn_imag(1:vla,1,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_imag(1:vla,1,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*q_up_real(1:vla,1,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_imag(1:vla,1,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_real(1:vla,1,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arry_real(1:vla)=help_arry_real(1:vla) +&
!!$                             ( q_up_real(1:vla,2,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_up_imag(1:vla,2,j,m)*p_up_imag(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,2,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             q_dn_imag(1:vla,2,j,m)*p_dn_imag(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_real(1:vla,2,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*q_up_imag(1:vla,2,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_real(1:vla,2,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_imag(1:vla,2,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arry_imag(1:vla)=help_arry_imag(1:vla) +&
!!$                             ( q_up_real(1:vla,2,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             q_up_imag(1:vla,2,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,2,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             q_dn_imag(1:vla,2,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_imag(1:vla,2,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*q_up_real(1:vla,2,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_imag(1:vla,2,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_real(1:vla,2,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arrz_real(1:vla)=help_arrz_real(1:vla) +&
!!$                             ( q_up_real(1:vla,3,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_up_imag(1:vla,3,j,m)*p_up_imag(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,3,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             q_dn_imag(1:vla,3,j,m)*p_dn_imag(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_real(1:vla,3,k,m)+&
!!$                             p_up_imag(1:vla,j,m)*q_up_imag(1:vla,3,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_real(1:vla,3,k,m)+&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_imag(1:vla,3,k,m))*&
!!$                             grdwts(1:vla)
!!$                        help_arrz_imag(1:vla)=help_arrz_imag(1:vla) +&
!!$                             ( q_up_real(1:vla,3,j,m)*p_up_imag(1:vla,k,m)-&
!!$                             q_up_imag(1:vla,3,j,m)*p_up_real(1:vla,k,m)+&
!!$                             q_dn_real(1:vla,3,j,m)*p_dn_imag(1:vla,k,m)-&
!!$                             q_dn_imag(1:vla,3,j,m)*p_dn_real(1:vla,k,m)+&
!!$                             p_up_real(1:vla,j,m)*q_up_imag(1:vla,3,k,m)-&
!!$                             p_up_imag(1:vla,j,m)*q_up_real(1:vla,3,k,m)+&
!!$                             p_dn_real(1:vla,j,m)*q_dn_imag(1:vla,3,k,m)-&
!!$                             p_dn_imag(1:vla,j,m)*q_dn_real(1:vla,3,k,m))*&
!!$                             grdwts(1:vla)
!!$                     enddo! loop over partners
!!$                     !call buggy("build_xc: building of xc_hamiltonian")
!!$                     ! building of xc_hamiltonian
!!$                     ! LDA-part of XC-Matrix
!!$                     help_arr_real(1:vla) = help_arr_real(1:vla)*dfdrho(1:vla,1)
!!$                     help_arr_imag(1:vla) = help_arr_imag(1:vla)*dfdrho(1:vla,1)
!!$                     ! GGA-part of XC-Matrix
!!$                     help_arr_real(1:vla) = help_arr_real(1:vla) +&
!!$                          help_x_up(1:vla)*&
!!$                          help_arrx_real(1:vla)
!!$                     help_arr_imag(1:vla) = help_arr_imag(1:vla) +&
!!$                          help_x_up(1:vla)*&
!!$                          help_arrx_imag(1:vla)
!!$                     help_arr_real(1:vla) = help_arr_real(1:vla) +&
!!$                          help_y_up(1:vla)*&
!!$                          help_arry_real(1:vla)
!!$                     help_arr_imag(1:vla) = help_arr_imag(1:vla) +&
!!$                          help_y_up(1:vla)*&
!!$                          help_arry_imag(1:vla)
!!$                     help_arr_real(1:vla) = help_arr_real(1:vla) +&
!!$                          help_z_up(1:vla)*&
!!$                          help_arrz_real(1:vla)
!!$                     help_arr_imag(1:vla) = help_arr_imag(1:vla) +&
!!$                          help_z_up(1:vla)*&
!!$                          help_arrz_imag(1:vla)
!!$                     do i_grid=1,vla
!!$                        help_xc_real=help_xc_real+&
!!$                             help_arr_real(i_grid)*grdwts(i_grid)
!!$                        help_xc_imag=help_xc_imag+&
!!$                             help_arr_imag(i_grid)*grdwts(i_grid)
!!$                     enddo
!!$                     ! building of xc_hamiltonian finished
!!$                     ham_xc_arr_real(counter)=ham_xc_arr_real(counter)+help_xc_real/partners(i)
!!$                     !ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)+help_xc_imag/partners(i)
!!$                     ! ***
!!$                     ! Note: Here we are calculating the matrix elements <k|Vxc|j> not <j|Vxc|k>
!!$                     ! ***
!!$                     ham_xc_arr_imag(counter)=ham_xc_arr_imag(counter)-help_xc_imag/partners(i)
!!$                     counter=counter+1
!!$                  enddo!loop over nu
!!$               enddo!loop over mu
!!$            enddo! loop over irreps
!!$            call xc_clear()
!!$            call stop_timer(timer_grid_xcbuild)
!!$            if(last) exit
!!$         enddo
!!$      endif
    end subroutine calc_xcks_nl_so

  end subroutine build_xc

   !*********************************************************

   subroutine xc_hamiltonian_store(th)
     ! Purpose: store the XC hamiltonian and XC energy to a
     !          readwriteblocked file in case the input switch
     !          "save_scfstate" is set.
     !
     ! subroutine called by: 'main_scf'
     !
     ! UB, 8/97
     !------------ Modules ----------------------------------------
     use iounitadmin_module, only: write_to_output_units, output_unit
     use output_module     , only: output_main_scf, output_data_saved
     use readwriteblocked_module
     !------------ Declaration of formal paramaters ---------------
     type(readwriteblocked_tapehandle), intent(inout) :: th
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind) :: n_spin
     real(kind=r8_kind) :: yes(1) = (/1.0_r8_kind/), no(1) = (/0.0_r8_kind/)
     !------------ Executable code --------------------------------
     n_spin = options_n_spin()

     if (output_main_scf) call write_to_output_units &
          ("XC_HAMILTONIAN_STORE: saving XC hamiltonian")
     if (output_data_saved) then
        write(output_unit,'(/ a     )')'Stored XC(num) matrix :'
     endif

     if (options_xcmode() /= xcmode_numeric_exch) then
        call readwriteblocked_write(no,th)
        if (output_data_saved) then
           write(output_unit,'(  a     )')'numeric exchange used ?'
           write(output_unit,'(4es20.13)')no
        endif
     else
        call readwriteblocked_write(yes,th)
        call readwriteblocked_write((/real(n_spin,r8_kind), &
                                      real(xc_length,r8_kind)/),th)
        if (options_spin_orbit) then
           call readwriteblocked_write(ham_xc_arr_real(1:xc_length:n_spin),th)
           call readwriteblocked_write(ham_xc_arr_imag(1:xc_length:n_spin),th)
        else
           call readwriteblocked_write(ham_xc_arr(1:xc_length:n_spin),th)
        endif
        if (output_data_saved) then
           write(output_unit,'(  a     )')'numeric exchange used ?'
           write(output_unit,'(4es20.13)')yes
           write(output_unit,'(  a     )')'n_spin, xc_length'
           write(output_unit,'(4es20.13)')(/real(n_spin,r8_kind), &
                                            real(xc_length,r8_kind)/)
           write(output_unit,'( a,i1,a )')'ham_xc_arr(1:xc_length:',n_spin,')'
           if (options_spin_orbit) then
              write(output_unit,'(4es20.13)')ham_xc_arr_real(1:xc_length:n_spin)
              write(output_unit,'(4es20.13)')ham_xc_arr_imag(1:xc_length:n_spin)
           else
              write(output_unit,'(4es20.13)')ham_xc_arr(1:xc_length:n_spin)
           endif
        endif
        if (n_spin > 1) then
           call readwriteblocked_write(ham_xc_arr(2:xc_length:2),th)
           if (output_data_saved) then
              write(output_unit,'(  a     )')'ham_xc_arr(2:xc_length:2)'
              write(output_unit,'(4es20.13)')ham_xc_arr(2:xc_length:2)
           endif
        endif
        call readwriteblocked_write((/exc_int/),th)
        if (output_data_saved) then
           write(output_unit,'(  a     )')'exc_int'
           write(output_unit,'(4es20.13)')(/exc_int/)
       endif
     endif

   end subroutine xc_hamiltonian_store

!*************************************************************
! record  1: numeric_exch
! record  2: n_spin, xc_length   [if numeric_exch]
! record  3: ham_xc_arr(spin=1)  [if numeric_exch]
! record  4: ham_xc_arr(spin=2)  [if numeric_exch & if n_spin > 1]
! record  5: exc_int             [if numeric_exch]
!*************************************************************

   subroutine xc_hamiltonian_recover(th)
     ! Purpose: recovers the XC hamiltonian and XC energy from a
     !          readwriteblocked file in case the input switch
     !          "read_scfstate" is set.
     !
     ! subroutine called by: 'main_scf'
     !
     ! UB 8/97
     !------------ Modules ----------------------------------------
     use iounitadmin_module, only: write_to_output_units, output_unit
     use output_module     , only: output_main_scf, output_data_read
     use readwriteblocked_module
     !------------ Declaration of formal_parameters ---------------
     type(readwriteblocked_tapehandle), intent(inout) :: th
     !** End of interface *****************************************
     !------------ Declaration of local variables   ---------------
     allocatable           :: buffer
     integer(kind=i4_kind) :: n_spin, spin_stored, arr_length
     real(kind=r8_kind)    :: dummy(1), arr_dim(2), num_exch(1), buffer(:), &
                              zero = 0.0_r8_kind, half = 0.5_r8_kind
     logical               :: numeric_exch
     integer               :: alloc_stat
     !------------ Executable code --------------------------------
     n_spin = options_n_spin()
     numeric_exch = options_xcmode() == xcmode_numeric_exch

     call readwriteblocked_read(num_exch,th)
     if (output_data_read) then
        write(output_unit,'(/ a     )')'Recovered XC(num) matrix :'
        write(output_unit,'(  a     )')'numeric exchange used ?'
        write(output_unit,'(4es20.13)')num_exch(1)
     endif
     if (num_exch(1) /= zero) then
        call readwriteblocked_read(arr_dim,th)
        spin_stored = int(arr_dim(1),i4_kind)
        arr_length = int(arr_dim(2),i4_kind)
        if (output_data_read) then
           write(output_unit,'(  a     )')'n_spin, xc_length'
           write(output_unit,'(4es20.13)')arr_dim(1:2)
        endif
        if (numeric_exch) then
           if (output_main_scf) call write_to_output_units &
                ("XC_HAMILTONIAN_REC: reading XC hamiltonian")
           mat_initialized = .false.
        else
           if (output_main_scf) call write_to_output_units &
                ("XC_HAMILTONIAN_REC: XC hamiltonian ignored")
           call readwriteblocked_skipread(arr_length+1,th)
           if (output_data_read) then
              write(output_unit,'(  a     )')'ham_xc_arr skipped'
              write(output_unit,'(  a     )')'exc_int    skipped'
           endif
           return
        endif
     else
        if (numeric_exch) then
           if (output_main_scf) call write_to_output_units &
                ("XC_HAMILTONIAN_REC: XC hamiltonian set to zero")
           ham_xc_arr = zero
           mat_initialized = .true.
        endif
        return
     endif

     ! start recovering XC hamiltonian now
     if (spin_stored > n_spin .and. output_main_scf) then
        call write_to_output_units &
             ("XC_HAMILTONIAN_REC: trying to convert spin-polarized data from")
        call write_to_output_units &
             ("                    tape into the spin-restricted data required")
        ! ham_xc_arr(:) := 1/2 * sum(s) ham_xc_arr(:,s)
     endif
     if (spin_stored < n_spin .and. output_main_scf) then
        call write_to_output_units &
             ("XC_HAMILTONIAN_REC: trying to convert spin-restricted data from")
        call write_to_output_units &
             ("                    tape into the spin-polarized data required")
        ! ham_xc_arr(:,s) := ham_xc_arr(:)
     endif

     if (options_spin_orbit) then
        call readwriteblocked_read(ham_xc_arr_real(1:xc_length:1),th)
        call readwriteblocked_read(ham_xc_arr_imag(1:xc_length:1),th)
     else
        call readwriteblocked_read(ham_xc_arr(1:xc_length:n_spin),th)
     endif
     if (output_data_read) then
        write(output_unit,'( a,i1,a )')'ham_xc_arr(1:xc_length:',n_spin,')'
        write(output_unit,'(4es20.13)')ham_xc_arr(1:xc_length:n_spin)
     endif
     if (spin_stored > 1) then
        if (n_spin > 1) then
           call readwriteblocked_read(ham_xc_arr(2:xc_length:2),th)
           if (output_data_read) then
              write(output_unit,'(  a     )')'ham_xc_arr(2:xc_length:2)'
              write(output_unit,'(4es20.13)')ham_xc_arr(2:xc_length:2)
           endif
        else
           allocate(buffer(xc_length),stat=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("XC_HAMILTONIAN_REC: allocation of buffer failed")
           call readwriteblocked_read(buffer,th)
           if (output_data_read) then
              write(output_unit,'(  a     )')'ham_xc_arr(2:xc_length:2)'
              write(output_unit,'(4es20.13)')buffer
           endif
           ham_xc_arr = ( ham_xc_arr + buffer ) * half
           deallocate(buffer,stat=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("XC_HAMILTONIAN_REC: deallocation of buffer failed")
        endif
     elseif (n_spin > 1) then
        ham_xc_arr(2:xc_length:2) = ham_xc_arr(1:xc_length:2)
        mat_initialized = .true. ! because the spin polarization is missing
     endif
     call readwriteblocked_read(dummy(1:1),th)
     exc_int = dummy(1)
     if (output_data_read) then
        write(output_unit,'(  a     )')'exc_int'
        write(output_unit,'(4es20.13)')dummy(1)
     endif

   end subroutine xc_hamiltonian_recover
   !*************************************************************

   !*************************************************************
   subroutine xc_pseudo2d_sym(ham_xc_arr)
     ! Purpose: in the case of pseudo 2d irreps this subroutine performs
     !          an additional symmetrization of the matrix elements, which is
     !          because the integration grid is not total symmetric.
     !          The procedure is analogous to the one for the analytical
     !          matrixelements. ( see integral_calc_quad_2cob3c)
     !
     ! subroutine called by: 'build_xc'
     !
     ! MS,8/97
     !------------ Modules ----------------------------------------
!    use iounitadmin_module
     implicit none
     real(r8_kind), intent(inout) :: ham_xc_arr(:)
     !** End of interface *****************************************
     !------------ Declaration of local variables   ---------------
     integer(kind=i4_kind) :: i_ir,i_ua_1,i_l_1,i_ua_2,i_l_2,n_exp1,n_exp2,&
          n_ind1d_1,n_ind1d_2,counter,counter2,counter3,counter4,offset,stride,&
          i_meta_1,i_meta_2,lmax,i_spin,i_ind_1,i_ind_2,i_exp1,i_exp2,ob_1,ob_2
     !------------ Executable code --------------------------------
     offset=0
     do i_ir=1,n_irrep ! loop over irreps
        if(pseudo(i_ir)) then
           ob_1=orbitalprojection_ob(i_ir,0,1)-1
           ob_2=orbitalprojection_ob(i_ir,0,1)-1
           do i_ua_1=1,n_unique_atoms
              do i_l_1=0,unique_atoms(i_ua_1)%lmax_ob
                 n_ind1d_1=unique_atoms(i_ua_1)%symadapt_partner(i_ir,i_l_1)%&
                      n_independent_fcts/2
                 n_exp1=unique_atoms(i_ua_1)%l_ob(i_l_1)%n_contracted_fcts+&
                      unique_atoms(i_ua_1)%l_ob(i_l_1)%n_uncontracted_fcts
                 do i_ua_2=1,i_ua_1
                    if(i_ua_1==i_ua_2) then
                       lmax=i_l_1
                    else
                       lmax=unique_atoms(i_ua_2)%lmax_ob
                    end if
                    do i_l_2=0,lmax
                       n_exp2=unique_atoms(i_ua_2)%l_ob(i_l_2)%n_contracted_fcts+&
                            unique_atoms(i_ua_2)%l_ob(i_l_2)%n_uncontracted_fcts
                       n_ind1d_2=unique_atoms(i_ua_2)%symadapt_partner&
                            (i_ir,i_l_2)%n_independent_fcts/2
                       if(i_ua_1==i_ua_2.and.i_l_1==i_l_2) then
                          do i_spin=1,ispin
                             ! first treat diagonal independent functions
                             counter=((orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-1-ob_1))/2+&
                                  orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2
                             counter2=n_exp2*n_ind1d_2+&
                                  orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2+&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)+&
                                  n_exp1*n_ind1d_1-ob_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+&
                                  n_exp1*n_ind1d_1-1-ob_1)/2
                             stride=(orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1)*ispin-ispin
                             counter=(counter-1)*ispin+i_spin+offset
                             counter2=(counter2-1)*ispin+i_spin+offset
                             do i_meta_1=1,n_exp1*n_ind1d_1
                                do i_meta_2=1,i_meta_1
                                   ham_xc_arr(counter)=0.5_r8_kind*(ham_xc_arr(counter)+&
                                        ham_xc_arr(counter2))
                                   ham_xc_arr(counter2)=ham_xc_arr(counter)
                                   counter=counter+ispin
                                   counter2=counter2+ispin
                                end do
                                counter=counter+stride
                                counter2=counter2+stride+n_exp1*n_ind1d_1*ispin
                             end do
                             ! now off diagonals
                             do i_ind_1=1,n_ind1d_1*n_exp1
                                counter4=(orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+n_exp1*n_ind1d_1+&
                                     i_ind_1-1)*(orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+&
                                     n_exp1*n_ind1d_1+i_ind_1-2)/2+orbitalprojection_ob(i_ir,i_l_2,i_ua_2)+&
                                     i_ind_1-1-ob_2
                                counter4=(counter4-1)*ispin+i_spin+offset
                                counter3=counter4
                                stride=(orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+&
                                     n_exp1*n_ind1d_1+i_ind_1-1)*ispin
                                do i_ind_2=i_ind_1,n_ind1d_1*n_exp1
                                   ham_xc_arr(counter4)=0.5_r8_kind*(ham_xc_arr(counter4)-&
                                        ham_xc_arr(counter3))
                                   ham_xc_arr(counter3)=-ham_xc_arr(counter4)
                                   counter4=counter4+ispin
                                   counter3=counter3+stride
                                   stride=stride+ispin
                                end do
                             end do
                          end do! loop over spin
                       else
                          do i_spin=1,ispin
                             stride=(orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1-&
                                  n_exp2*n_ind1d_2+1)*ispin-ispin
                             counter=((orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1-1))/2+&
                                  orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2
                             counter2=n_exp2*n_ind1d_2+orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2+&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+n_exp1*n_ind1d_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+n_exp1*n_ind1d_1-1)/2
                             counter3=((orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1-1))/2+&
                                  orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2+n_exp2*n_ind1d_2
                             counter4=orbitalprojection_ob(i_ir,i_l_2,i_ua_2)-ob_2+&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+n_exp1*n_ind1d_1)*&
                                  (orbitalprojection_ob(i_ir,i_l_1,i_ua_1)-ob_1+n_exp1*n_ind1d_1-1)/2
                             counter=(counter-1)*ispin+i_spin+offset
                             counter2=(counter2-1)*ispin+i_spin+offset
                             counter3=(counter3-1)*ispin+i_spin+offset
                             counter4=(counter4-1)*ispin+i_spin+offset
                             do i_ind_1=1,n_ind1d_1
                                do i_exp1=1,n_exp1
                                   do i_ind_2=1,n_ind1d_2
                                      do i_exp2=1,n_exp2
                                         ! first diagonal elements
                                         ham_xc_arr(counter)=0.5_r8_kind*&
                                              (ham_xc_arr(counter)+ham_xc_arr(counter2))
                                         ham_xc_arr(counter2)=ham_xc_arr(counter)
                                         ! now off diagonals
                                         ham_xc_arr(counter3)=0.5_r8_kind*&
                                              (ham_xc_arr(counter3)-ham_xc_arr(counter4))
                                         ham_xc_arr(counter4)=-ham_xc_arr(counter3)
                                          counter=counter+ispin
                                          counter2=counter2+ispin
                                          counter3=counter3+ispin
                                          counter4=counter4+ispin
                                      end do ! loop over i_exp2
                                   end do ! loop over i_ind_2
                                   counter=counter+stride
                                   counter2=counter2+stride+n_exp1*n_ind1d_1*ispin
                                   counter3=counter3+stride
                                   counter4=counter4+stride+n_exp1*n_ind1d_1*ispin
                                   stride=stride+ispin
                                end do ! loop over i_exp1
                             end do ! loop over i_ind1
                          end do ! loop over spin
                       end if
                    end do
                end do
            end do
         end do
       endif
    offset=offset+(dims(i_ir)*(dims(i_ir)+1))/2*ispin
   end do  ! loop over irreps
  end subroutine xc_pseudo2d_sym
   !************************************************************

 end module xc_hamiltonian

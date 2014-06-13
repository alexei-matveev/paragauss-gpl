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
!================================================================
! Public interface of module
!================================================================
module density_calc_module
!---------------------------------------------------------------
!
!  Purpose: Calculation of the density and gradient of gensity
!           for a given set of gridpoints
!
!  Author: MS
!  Date: 4/96
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
! Description: Performance of routines has been increased
!              ordering of indices in gradient variables has
!              been changed
!              Subroutines with nl in the name are necessary for
!              GGA functionals
!
! Modification (Please copy before editing)
! Author: MS
! Date:   3/97
! Description: additional subroutines, which are necessary for calculation
!              of energy gradients (density_calc_ph, density_calc_nl_ph)
!              In these subroutines also gradients with respect to nuclear
!              are calculated displacements
!
! Modification (Please copy before editing)
! Author: UB
! Date:   6/97 (GGA version: 8/98)
! Description: Subroutines to compute the fitted charge density and its
!              various gradients added
!
! Modification (Please copy before editing)
! Author: HH
! Date:   11/97
! Description: Added new subroutine "density_calc_response()" which is an
!              modifed version of "density_calc()":
!              In addition to the density "rho" also occupied and unoccupied
!              orbitals are calculated. These quantities are needed in the
!              response_module.
!              Based on "density_calc_response()" from modified V11alpha
!              version by MS (see view "ttfs_hh_V11alpha").
!
! Modification (Please copy before editing)
! Author: MM
! Date:   12/97
! Description: calculation of the density for spin orbit
!
! Modification
! Author: TS
! Date:   07/09
! Description: extension to mgga
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
! define FPP_TIMERS 2
# include "def.h"
!## define FPP_NOBLAS 0
  use type_module  ! type specification parameters
  use orbitalstore_module
  use orbital_module
  use iounitadmin_module
  use occupied_levels_module
  use datatype, only: arrmat4,arrmat5,arrmat6,arrmat7
  use options_module, only: options_spin_orbit
  use unique_atom_module
  use comm_module
  USE_MEMLOG
  implicit none
  private
  save

#ifdef WITH_EXPERIMENTAL
# undef WITH_SECDER
#endif

     type, public :: density_matrix_type
        ! to contain calculated sym. orbitals of one irrep i
        real(kind=r8_kind), pointer :: d_ur_ur(:,:,:,:) !density matrix - the product of up real and up real wave functions
        real(kind=r8_kind), pointer :: d_ui_ui(:,:,:,:)
        real(kind=r8_kind), pointer :: d_dr_dr(:,:,:,:)
        real(kind=r8_kind), pointer :: d_di_di(:,:,:,:)
        real(kind=r8_kind), pointer :: d_ur_ui(:,:,:,:)
        real(kind=r8_kind), pointer :: d_dr_di(:,:,:,:)
        real(kind=r8_kind), pointer :: d_ur_di(:,:,:,:)
        real(kind=r8_kind), pointer :: d_ui_dr(:,:,:,:)
        real(kind=r8_kind), pointer :: d_ui_di(:,:,:,:)
        real(kind=r8_kind), pointer :: d_ur_dr(:,:,:,:)
        real(kind=r8_kind), pointer :: d(:,:,:,:)
     end type density_matrix_type

  ! FIXME: make them public to be used also in density_calc_cpks:
  integer(i4_kind), public, protected :: n_irrep
  integer(i4_kind), allocatable, public, protected :: partners(:) ! (n_irrep)


!== Interrupt end of public interface of module =================

  integer(kind=i4_kind) :: alloc_stat(34)=1
  integer(kind=i4_kind) :: ispin, max_equal
  integer(kind=i4_kind),allocatable :: dims(:), &
       & dims_core(:),dims_fit(:)
  real(kind=r8_kind), allocatable :: coeff_core(:) !AM???

  integer(i4_kind), parameter, private ::&
       & VALUE=0, X=1, Y=2, Z=3,&
       & UP=1, DN=2,&
       & UPUP=1, DNDN=2, UPDN=3,&
       & RE=1, IM=2

  real(r8_kind), parameter :: ZERO    = 0.0_r8_kind  &
                            , ONEHALF = 0.5_r8_kind  &
                            , ONE     = 1.0_r8_kind  &
                            , TWO     = 2.0_r8_kind
  logical :: fitted_density_is_closed_priv=.true.

  !------------ public functions and subroutines ------------------
  public :: density_calc,&
       density_calc_nl,&
       density_calc_nl_v2,&
       density_calc_setup,&
       density_calc_close,&
       density_calc_ph,&
       density_calc_ph_nl,&
       fitted_density_calc_setup,&
       fitted_density_calc_close,&
       fitted_density_calc,&
       fitted_core_dens_calc_close,&
       fitted_density_is_closed,&
       print_alloc_density_calc

#ifdef WITH_CORE_DENS
  public :: fitted_core_dens_calc_setup
  public :: fitted_core_density_calc
#endif

#ifdef WITH_GTENSOR
  public :: density_calc_hfc
#endif

#ifdef WITH_RESPONSE
  public :: density_calc_response
#endif

!================================================================
! End of public interface of module
!================================================================

  !------------ public types of subroutines formal parameters ---
  public :: orbital_type, orbital_gradient_type

  integer, private :: refcount = 0 ! save, counts module users

contains


  subroutine density_calc_setup()
    !
    ! Performs   the    necessary   allocations;   additionally   some
    ! informations  from  symmetry_data_module  are stored  in  module
    ! private variables.
    !
    ! FIXME: is there a way to get rid of the module state completely?
    !
    use symmetry_data_module, only: ssym ! description of irreps
    implicit none
    !** End of interface **************************************

    ! One more user of the module:
    refcount = refcount + 1

    if (refcount > 1) then
       ! Someone already initialized the module, nothing to be done:
       return
    endif

    ispin = ssym % n_spin
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym % n_proj_irrep

       allocate (dims(n_irrep), partners(n_irrep), stat=alloc_stat(1))
       ASSERT(alloc_stat(1).eq.0)

       dims = ssym % dim_proj
       partners = ssym % partner_proj
    else
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       n_irrep = ssym % n_irrep

       allocate (dims(n_irrep), partners(n_irrep), stat=alloc_stat(1))
       ASSERT(alloc_stat(1).eq.0)

       dims = ssym % dim
       partners = ssym % partner
    endif
    alloc_stat(2) = 0 ! partners

    ASSERT(n_unique_atoms>0)
    max_equal = maxval (unique_atoms(:) % n_equal_atoms)
  end subroutine density_calc_setup


  subroutine density_calc_close()
    !
    ! Deallocation of module private variables.
    !
    implicit none
    !** End of interface **************************************

    ! One user of the module fewer:
    refcount = refcount - 1

    if (refcount > 0) then
       ! Someone still  needs the  module, postpone cleanup  until the
       ! next call:
       return
    endif

    deallocate (dims, partners, stat=alloc_stat(1))
    ASSERT(alloc_stat(1).eq.0)
    alloc_stat(1) = 1 ! dims
    alloc_stat(2) = 1 ! partners
  end subroutine density_calc_close


  subroutine fitted_density_calc_setup()
    ! purpose : routine to perform the necessary allocations;
    !           additionally some private variable as initialized
    !** End of interface **************************************
    use options_module, only: options_n_spin
    use symmetry_data_module, only: get_totalsymmetric_irrep

    type(unique_atom_type)        , pointer :: ua
    type(unique_atom_basis_type)  , pointer :: uab
    type(unique_atom_partner_type), pointer :: uap
    integer(kind=i4_kind) :: i_ir,i_ua,i_ma,i_bas,n_orbs

!..............................................................................
! dims_fit(i_ma)        : the number of local fitting functions
!                     associated with the moving unique atom i_ma
! ua%N_glob_cons_ch : the number of gobally contracted fitting functions
!                     associated with the unique atom ua = unique_atoms(i_ua)
!..............................................................................
!MF: Why outcommended???
    ! ispin=options_n_spin()
    ispin=options_n_spin()
!MF: is dims (or new dims_fit) ever needed?
    allocate( dims_fit(N_moving_unique_atoms), stat = alloc_stat(13) )
    ASSERT(alloc_stat(13).eq.0)
    max_equal=0
    i_ir = get_totalsymmetric_irrep()
    do i_ua=1,n_unique_atoms
       ua => unique_atoms(i_ua)
       max_equal = max(max_equal,ua%n_equal_atoms)
       i_ma = ua%moving_atom
       if (i_ma > 0) then
          n_orbs = 0
          do i_bas=-1,ua%lmax_ch
            select case (i_bas)
            case (-1) ! s-type
              uab => ua%l_ch(0)
              uap => ua%symadapt_partner(i_ir,0)
            case (0) ! r^2-type
              uab => ua%r2_ch
              uap => ua%symadapt_partner(i_ir,0)
            case default ! l>0-type
              uab => ua%l_ch(i_bas)
              uap => ua%symadapt_partner(i_ir,i_bas)
            end select
            n_orbs = n_orbs + uap%N_independent_fcts &
                   * ( uab%N_uncontracted_fcts + uab%N_contracted_fcts )
          end do
          dims_fit(i_ma) = n_orbs
       end if
    end do
    fitted_density_is_closed_priv=.false.

  end subroutine fitted_density_calc_setup
  !***************************************************************

#ifdef WITH_CORE_DENS
  subroutine fitted_core_dens_calc_setup()
    ! purpose : routine to perform the necessary allocations
    !           additionally some private variable as initialized
    !** End of interface **************************************
    use options_module, only: options_n_spin
    use operations_module,only: operations_core_density
    use fit_coeff_module, only: fit_coeff_n_cd
    type(unique_atom_type)        , pointer :: ua
    type(unique_atom_atomic_dens_type),pointer :: uac
    integer(kind=i4_kind) :: i_ua,i_bas,&
                             n_core_orbs, n_coeff_index, i_exp

   !------------- excutable code -------------------------------------
!MF: Why outcommended?
    !ispin=options_n_spin()
    ispin=options_n_spin()
    allocate( dims_core(N_unique_atoms), stat = alloc_stat(3) )
    ASSERT(alloc_stat(3).eq.0)
     allocate( coeff_core(fit_coeff_n_cd()), stat = alloc_stat(4) )
    if( alloc_stat(4) /= 0) call error_handler &
        ('fitted_density_calc_setup: allocation of coeff_core failed')

    n_coeff_index = 0
    do i_ua=1,n_unique_atoms
       ua => unique_atoms(i_ua)
          if ( ua%zc /= 0.0_r8_kind .and. .not.operations_core_density ) then
             n_core_orbs = 0
             do i_bas=-1,0
                select case (i_bas)
                case(-1) ! s-type for core density
                    uac => ua%s_core
                case(0)  ! r^2-type for core density
                    uac => ua%r2_core
                end select
                n_core_orbs = n_core_orbs + uac%N_exponents
                do i_exp = 1, uac%n_exponents
                    n_coeff_index = n_coeff_index + 1
                    coeff_core(n_coeff_index) = uac%contractions(i_exp)
                end do
             end do
            dims_core(i_ua) = n_core_orbs
          end if
    end do
    if ( fit_coeff_n_cd() /= n_coeff_index ) call error_handler(&
            'fitted_density_calc_setup:fit_coeff_n_cd /= n_coeff_index')

  end subroutine fitted_core_dens_calc_setup
#endif

  !***************************************************************

  subroutine fitted_core_dens_calc_close()
    ! purpose : deallocation of module private arrays
    !** End of interface **************************************

    deallocate( dims_core, stat = alloc_stat(3) )
    if( alloc_stat(3) /= 0) call error_handler &
        ('fitted_density_calc_close: deallocation of dims_core failed')
    alloc_stat(3)=1
    deallocate( coeff_core, stat = alloc_stat(4) )
    if( alloc_stat(4) /= 0) call error_handler &
        ('fitted_density_calc_close: deallocation of coeff_core failed')
    alloc_stat(4)=1
  end subroutine fitted_core_dens_calc_close

  !***************************************************************
  subroutine fitted_density_calc_close()
    ! purpose : deallocation of module private arrays
    !** End of interface **************************************

    deallocate( dims_fit, stat = alloc_stat(1)  )
    if( alloc_stat(1) /= 0) call error_handler &
        ('fitted_density_calc_close: deallocation of dims failed')
    alloc_stat(1)=1
    fitted_density_is_closed_priv=.true.
  end subroutine fitted_density_calc_close

  !***************************************************************

  subroutine density_calc (vl, rho, orbs_ob, orbs_spinor_ob, s_dens)
    !
    ! Calculation of the density from orbitlas (or spinors).
    !
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:, :)
    ! spin density s_dens(vl, 3):
    real(r8_kind), intent(out), optional :: s_dens(:, :)
    type(orbital_type), intent(in), optional :: orbs_ob(:)
    type(spinor_type), intent(in), optional :: orbs_spinor_ob(:)
    !** End of interface *****************************************


#ifndef FPP_NOBLAS
           !           Blas Version

           real(kind=r8_kind)  :: occ_real
           real(kind=r8_kind), allocatable :: phi_arr(:,:)
           real(kind=r8_kind), allocatable :: phi_up_arr_real(:,:),phi_up_arr_imag(:,:)
           real(kind=r8_kind), allocatable :: phi_dn_arr_real(:,:),phi_dn_arr_imag(:,:)
           ! phi_arr(vl, number of occupied orbitals)
           real(kind=r8_kind),pointer :: ob(:,:,:)
!!$           real(kind=r8_kind),pointer :: ob_up_real(:,:,:),ob_up_imag(:,:,:),ob_dn_real(:,:,:),ob_dn_imag(:,:,:)
           real(kind=r8_kind), allocatable :: rho_up(:),rho_dn(:)
           integer(i4_kind) :: s, i, j, l, eig_dim, occ_dim
           rho(1:vl,:)=0.0_r8_kind

           ! NOTE:  vl is  allways  smaller or  equal  then the  lower
           ! dimension  of  orbs_ob()%o,  orbs_spinor_ob()%o. Rho  has
           ! dimension vl (FIXME: is it used anywhere?).
           DPRINT 'use Blas option'

           DPRINT 'density_calc', vl, shape(rho)
           if (present(s_dens)) then
              DPRINT 'density_calc',shape(s_dens)
           end if
           if (options_spin_orbit) then
              !
              ! SPIN ORBIT
              !
!!$              if (spin_orbit_polarized) then
              if (present(s_dens)) then
                 s_dens(1:vl,:)=0.0_r8_kind
                 allocate(rho_up(vl),rho_dn(vl),stat=alloc_stat(5))
                 if(alloc_stat(5)/=0) call error_handler&
                      ('density_calc: allocation rho_up failed')
              endif
              do i=1,n_irrep
                 DPRINT ' IRREP = ', i, shape(eigvec_occ_real(i)%m)
!!$                  ob_up_real=> orbs_spinor_ob(i)%spinor(1)%o
!!$                  ob_up_imag=> orbs_spinor_ob(i)%spinor(1)%o_imag
!!$                  ob_dn_real=> orbs_spinor_ob(i)%spinor(2)%o
!!$                  ob_dn_imag=> orbs_spinor_ob(i)%spinor(2)%o_imag
                  associate( occ=>occ_num_occ(i)%m                             &
                           , eigv_real=>eigvec_occ_real(i)%m                   &
                           , eigv_imag=>eigvec_occ_imag(i)%m )
                  eig_dim=size(eigv_real,1)
                  occ_dim=size(eigv_real,2)
                  if (occ_dim > 0) then
                  allocate(phi_up_arr_real(vl,occ_dim),&
                      phi_up_arr_imag(vl,occ_dim),&
                      phi_dn_arr_real(vl,occ_dim),&
                      phi_dn_arr_imag(vl,occ_dim),&
                      stat=alloc_stat(6))
                  if(alloc_stat(6)/=0) call error_handler&
                       ('density_calc: allocation phi_arr failed')
                  phi_up_arr_real = 0.0_r8_kind
                  phi_up_arr_imag = 0.0_r8_kind
                  phi_dn_arr_real = 0.0_r8_kind
                  phi_dn_arr_imag = 0.0_r8_kind
                  do l=1,partners(i)
                      ! phi_up_real
                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, RE, UP), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_real(1, 1), eig_dim, 0.0_r8_kind, &
                           phi_up_arr_real, vl)

                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, -1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, IM, UP), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_imag(1, 1), eig_dim, 1.0_r8_kind, &
                           phi_up_arr_real, vl)

                      ! phi_up_imag
                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, IM, UP), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_real(1, 1), eig_dim, 0.0_r8_kind, &
                           phi_up_arr_imag, vl)

                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, RE, UP), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_imag(1, 1), eig_dim, 1.0_r8_kind, &
                           phi_up_arr_imag, vl)

                      ! phi_dn_real
                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, RE, DN), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_real(1, 1), eig_dim, 0.0_r8_kind, &
                           phi_dn_arr_real, vl)

                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, -1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, IM, DN), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_imag(1, 1), eig_dim, 1.0_r8_kind, &
                           phi_dn_arr_real, vl)

                      ! phi_dn_imag
                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, IM, DN), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_real(1, 1), eig_dim, 0.0_r8_kind, &
                           phi_dn_arr_imag, vl)

                      call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                           orbs_spinor_ob(i)%o(1, 1, l, RE, DN), size(orbs_spinor_ob(i)%o, 1), &
                           eigv_imag(1, 1), eig_dim, 1.0_r8_kind, &
                           phi_dn_arr_imag, vl)

                     ! if (spin_orbit_polarized) then
                         if (present(s_dens)) then
                         !
                         ! OPEN SHELL
                         !
                         ! calculate not only density but also spin density
                         do j=1,occ_dim
                            occ_real=occ(j,1)/partners(i)
                            rho_up(1:vl)=&
                                 (phi_up_arr_real(:,j)*phi_up_arr_real(:,j)+&
                                  phi_up_arr_imag(:,j)*phi_up_arr_imag(:,j))*occ_real
                            rho_dn(1:vl)=&
                                  (phi_dn_arr_real(:,j)*phi_dn_arr_real(:,j)+&
                                  phi_dn_arr_imag(:,j)*phi_dn_arr_imag(:,j))*occ_real
                            rho(1:vl,1)    =  rho(1:vl,1)    + rho_up(1:vl) + rho_dn(1:vl)
                            s_dens(1:vl,3) =  s_dens(1:vl,3) + rho_up(1:vl) - rho_dn(1:vl)
                            s_dens(1:vl,1) =  s_dens(1:vl,1) + (phi_up_arr_real(:,j)* phi_dn_arr_real(:,j) +&
                                  phi_up_arr_imag(:,j)* phi_dn_arr_imag(:,j))*occ_real*2
                            s_dens(1:vl,2) =  s_dens(1:vl,2) + (phi_up_arr_real(:,j)* phi_dn_arr_imag(:,j) -&
                                  phi_up_arr_imag(:,j)* phi_dn_arr_real(:,j))*occ_real*2
                         end do! loop over j
                       else ! spin_orbit_polarized
                         !
                         ! CLOSED SHELL
                         !
                         do j=1,occ_dim
                            occ_real=occ(j,1)/partners(i)
                            rho(1:vl,1)=rho(1:vl,1)+&
                                 (phi_up_arr_real(:,j)*phi_up_arr_real(:,j)+&
                                  phi_up_arr_imag(:,j)*phi_up_arr_imag(:,j)+&
                                  phi_dn_arr_real(:,j)*phi_dn_arr_real(:,j)+&
                                  phi_dn_arr_imag(:,j)*phi_dn_arr_imag(:,j))*occ_real
                         end do! loop over j
                       endif ! spin_orbit_polarized
                     end do! loop over partners
                  deallocate(phi_up_arr_real,phi_up_arr_imag,&
                             phi_dn_arr_real,phi_dn_arr_imag,stat=alloc_stat(6))
                  if(alloc_stat(6)/=0) call error_handler&
                       ('density_calc: deallocation phi_arr_real/imag failed')
                  alloc_stat(6)=1
                  endif
                  end associate
              end do
!!$              if (spin_orbit_polarized) then
              if (present(s_dens)) then
                 deallocate(rho_up,rho_dn,stat=alloc_stat(5))
                 if(alloc_stat(5)/=0) call error_handler&
                      ('density_calc: deallocation rho_up failed')
                 alloc_stat(5)=1
              endif
           else ! options_spin_orbit
              !
              ! STANDARD SCF (NO SPIN ORBIT)
              !
              do i=1,n_irrep
                 ob=> orbs_ob(i)%o
                 associate( occ=>occ_num_occ(i)%m )
                 eig_dim=size(eigvec_occ(i)%m,1)
                 occ_dim=size(eigvec_occ(i)%m,2)

                 if (occ_dim > 0) then
                 allocate(phi_arr(vl,occ_dim),stat=alloc_stat(7))
                 if(alloc_stat(7)/=0) call error_handler&
                      ('density_calc: allocation phi_arr failed')
                 do s=1,ispin
                    do l=1,partners(i)
                       call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                            ob(:, :, l), size(ob, 1), &
                            eigvec_occ(i)%m(:, :, s), eig_dim, 0.0_r8_kind, &
                            phi_arr, vl)
                       do j=1,occ_dim
                          occ_real=occ(j,s)/partners(i)
                          rho(1:vl,s)=rho(1:vl,s)+&
                               phi_arr(:,j)*phi_arr(:,j)*occ_real
                       end do! loop over j
                    end do! loop over partners
                 end do
                 deallocate(phi_arr,stat=alloc_stat(7))
                 if(alloc_stat(7)/=0) call error_handler&
                      ('density_calc: deallocation phi_arr failed')
                 alloc_stat(7)=1
                 endif
                 end associate
              end do
           endif ! options_spin_orbit

#else
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning

           !           f77-like version

           real(kind=r8_kind)  :: phi(vl),occ_real
           real(kind=r8_kind)  :: phi_up_real(vl),phi_up_imag(vl),&
                phi_down_real(vl),phi_down_imag(vl)
           real(kind=r8_kind),pointer :: ob(:,:,:)
           real(kind=r8_kind),pointer :: ob_up_real(:,:,:),ob_up_imag(:,:,:),&
                ob_dn_real(:,:,:),ob_dn_imag(:,:,:)
           integer(i4_kind) :: s,i,j,l,m,eig_dim,occ_dim
           rho(1:vl,:)=0.0_r8_kind

           if (options_spin_orbit) then
              !
              ! SPIN ORBIT
              !
              do i=1,n_irrep ! irreps
                 ob_up_real=> orbs_spinor_ob(i)%o(:,:,:,1,1) !spinor(1)%o
                 ob_up_imag=> orbs_spinor_ob(i)%o(:,:,:,2,1) !spinor(1)%o_imag
                 ob_dn_real=> orbs_spinor_ob(i)%o(:,:,:,1,2) !spinor(2)%o
                 ob_dn_imag=> orbs_spinor_ob(i)%o(:,:,:,2,2) !spinor(2)%o_imag
                 associate( occ=>occ_num_occ(i)%m                              &
                          , eigv_real=>eigvec_occ_real(i)%m                    &
                          , eigv_imag=>eigvec_occ_imag(i)%m )
                 eig_dim=size(eigv_real,1)
                 occ_dim=size(eigv_real,2)
                 if (occ_dim > 0) then
                    do j=1,occ_dim ! orbitals
                       occ_real=occ(j,1)/partners(i)
                       do l=1,partners(i) ! partners
                          phi_up_real=0.0_r8_kind
                          phi_up_imag=0.0_r8_kind
                          phi_down_real=0.0_r8_kind
                          phi_down_imag=0.0_r8_kind
                          do m=1,eig_dim ! basisfunctions
                             phi_up_real = phi_up_real +&
                                  ob_up_real(1:vl,m,l)*eigv_real(m,j)&
                                  - ob_up_imag(1:vl,m,l)*eigv_imag(m,j)
                             phi_up_imag = phi_up_imag +&
                                  ob_up_imag(1:vl,m,l)*eigv_real(m,j)&
                                  + ob_up_real(1:vl,m,l)*eigv_imag(m,j)
                             phi_down_real = phi_down_real +&
                                  ob_dn_real(1:vl,m,l)*eigv_real(m,j)&
                                  - ob_dn_imag(1:vl,m,l)*eigv_imag(m,j)
                             phi_down_imag = phi_down_imag +&
                                  ob_dn_imag(1:vl,m,l)*eigv_real(m,j)&
                                 + ob_dn_real(1:vl,m,l)*eigv_imag(m,j)
                          end do! basisfunctions
                          ! calculation of density
                          rho(1:vl,1)=rho(1:vl,1)+&
                               (phi_up_real*phi_up_real + phi_up_imag*phi_up_imag +&
                               phi_down_real*phi_down_real + phi_down_imag*phi_down_imag)*&
                               occ_real
                       end do! partners
                    end do! orbitals
                 end if
                 end associate
              end do ! irreps
           else ! options_spin_orbit
             !
             ! STANDARD SCF (NO SPIN ORBIT)
             !
             do i=1,n_irrep
                ob=> orbs_ob(i)%o
                associate( occ=>occ_num_occ(i)%m )
                eig_dim=size(eigvec_occ(i)%m,1)
                occ_dim=size(eigvec_occ(i)%m,2)
                if (occ_dim > 0) then
                do s=1,ispin
                   do j=1,occ_dim
                      occ_real=occ(j,s)/partners(i)
                      do l=1,partners(i)
                         phi=0.0_r8_kind
                         do m=1,eig_dim
                            phi=phi+eigvec_occ(i)%m(m,j,s)&
                                 *ob(1:vl,m,l)
                         end do
                         ! calculation of density
                         rho(1:vl,s)=rho(1:vl,s)+phi*phi*occ_real
                      end do
                   end do
                end do
                endif
                end associate
             end do
           endif! options_spin_orbit
#endif


 end subroutine density_calc

#ifdef WITH_GTENSOR
 subroutine density_calc_hfc (vl, rho, dens, s_dens)
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:, :)
    ! spin density s_dens(vl, 3):
    real(r8_kind), intent(out), optional :: s_dens(:,:)
    type(density_matrix_type), pointer, optional :: dens(:)
    !** End of interface *****************************************

               !           f77-like version

           real(kind=r8_kind)  :: phi(vl),occ_real
           real(kind=r8_kind)  :: rho_up(vl), rho_dn(vl)
           real(kind=r8_kind)  :: phi_ui_phi_ui(vl), phi_ur_phi_ur(vl),phi_dr_phi_dr(vl),&
                phi_di_phi_di(vl), phi_ur_phi_dr(vl), phi_ui_phi_di(vl), phi_ui_phi_dr(vl),  phi_ur_phi_di(vl)
           real(kind=r8_kind),pointer :: d_ur_ur(:,:,:,:), d_ui_ui(:,:,:,:), d_dr_dr(:,:,:,:), d_di_di(:,:,:,:), &
                d_ur_ui(:,:,:,:), d_dr_di(:,:,:,:), d_ur_di(:,:,:,:), d_ui_dr(:,:,:,:), d_ui_di(:,:,:,:), d_ur_dr(:,:,:,:)
           integer(i4_kind) :: s,i,j,l,m,n,eig_dim,occ_dim


!init
           rho(1:vl,:)=0.0_r8_kind

              if (options_spin_orbit) then
              !
              ! SPIN ORBIT
              !

              do i=1,n_irrep ! irreps

                 d_ur_ur => dens(i)%d_ur_ur
                 d_ui_ui => dens(i)%d_ui_ui
                 d_dr_dr => dens(i)%d_dr_dr
                 d_di_di => dens(i)%d_di_di
                 d_ur_ui => dens(i)%d_ur_ui
                 d_dr_di => dens(i)%d_dr_di
                 d_ur_di => dens(i)%d_ur_di
                 d_ui_dr => dens(i)%d_ui_dr
                 d_ui_di => dens(i)%d_ui_di
                 d_ur_dr => dens(i)%d_ur_dr
                 associate( occ=>occ_num_occ(i)%m                              &
                          , eigv_real=>eigvec_occ_real(i)%m                    &
                          , eigv_imag=>eigvec_occ_imag(i)%m )
                 eig_dim=size(eigv_real,1)
                 occ_dim=size(eigv_real,2)
                 if (occ_dim > 0) then
                    do j=1,occ_dim ! orbitals
                       occ_real=occ(j,1)/partners(i)
                       do l=1,partners(i) ! partners

                          phi_ur_phi_ur = 0.0_r8_kind
                          phi_ui_phi_ui = 0.0_r8_kind
                          phi_dr_phi_dr = 0.0_r8_kind
                          phi_di_phi_di = 0.0_r8_kind
                          phi_ur_phi_dr = 0.0_r8_kind
                          phi_ui_phi_di = 0.0_r8_kind
                          phi_ui_phi_dr = 0.0_r8_kind
                          phi_ur_phi_di = 0.0_r8_kind

                          do m=1,eig_dim ! basisfunctions
                             do n=1,eig_dim
                             phi_ur_phi_ur = phi_ur_phi_ur +&
                                   eigv_real(m,j)*eigv_real(n,j)*d_ur_ur(1:vl,m,n,l)&
                                  + eigv_imag(n,j)*eigv_imag(n,j)*d_ui_ui(1:vl,m,n,l)&
                                  - eigv_real(m,j)*eigv_imag(n,j)*d_ur_ui(1:vl,m,n,l)&
                                  - eigv_real(n,j)*eigv_imag(m,j)*d_ur_ui(1:vl,n,m,l)

                               phi_ui_phi_ui = phi_ui_phi_ui +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_ui_ui(1:vl,m,n,l)&
                                  +eigv_imag(m,j)*eigv_imag(n,j)*d_ur_ur(1:vl,m,n,l)&
                                   +eigv_real(m,j)*eigv_imag(n,j)*d_ur_ui(1:vl,n,m,l)&
                                  +eigv_real(n,j)*eigv_imag(m,j)*d_ur_ui(1:vl,m,n,l)

                                  phi_dr_phi_dr = phi_dr_phi_dr +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_dr_dr(1:vl,m,n,l)&
                                 +eigv_imag(m,j)*eigv_imag(n,j)*d_di_di(1:vl,m,n,l)&
                                 -eigv_real(m,j)*eigv_imag(n,j)*d_dr_di(1:vl,m,n,l)&
                                 -eigv_real(n,j)*eigv_imag(m,j)*d_dr_di(1:vl,n,m,l)

                             phi_di_phi_di = phi_di_phi_di +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_di_di(1:vl,m,n,l)&
                                 +eigv_imag(m,j)*eigv_imag(n,j)*d_dr_dr(1:vl,m,n,l)&
                                 +eigv_real(m,j)*eigv_imag(n,j)*d_dr_di(1:vl,n,m,l)&
                                 +eigv_real(n,j)*eigv_imag(m,j)*d_dr_di(1:vl,m,n,l)


                             phi_ur_phi_dr = phi_ur_phi_dr +&
                                 eigv_real(m,j)*eigv_real(n,j)*d_ur_dr(1:vl,m,n,l)&
                                -eigv_real(m,j)*eigv_imag(n,j)*d_ur_di(1:vl,m,n,l)&
                                -eigv_imag(m,j)*eigv_real(n,j)*d_ui_dr(1:vl,m,n,l)&
                                +eigv_imag(m,j)*eigv_imag(n,j)*d_ui_di(1:vl,m,n,l)



                             phi_ui_phi_di = phi_ui_phi_di +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_ui_di(1:vl,m,n,l)&
                                 +eigv_real(m,j)*eigv_imag(n,j)*d_ui_dr(1:vl,m,n,l)&
                                 +eigv_imag(m,j)*eigv_real(n,j)*d_ur_di(1:vl,m,n,l)&
                                 +eigv_imag(m,j)*eigv_imag(n,j)*d_ur_dr(1:vl,n,m,l)


                             phi_ui_phi_dr = phi_ui_phi_dr +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_ui_dr(1:vl,m,n,l)&
                                 -eigv_real(m,j)*eigv_imag(n,j)*d_ui_di(1:vl,m,n,l)&
                                 +eigv_imag(m,j)*eigv_real(n,j)*d_ur_dr(1:vl,m,n,l)&
                                 -eigv_imag(m,j)*eigv_imag(n,j)*d_ur_di(1:vl,m,n,l)


                              phi_ur_phi_di = phi_ur_phi_di +&
                                  eigv_real(m,j)*eigv_real(n,j)*d_ur_di(1:vl,m,n,l)&
                                 +eigv_real(m,j)*eigv_imag(n,j)*d_ur_dr(1:vl,m,n,l)&
                                 -eigv_imag(m,j)*eigv_real(n,j)*d_ui_di(1:vl,m,n,l)&
                                 -eigv_imag(m,j)*eigv_imag(n,j)*d_ui_dr(1:vl,m,n,l)

                          end do
                       end do! basisfunctions
                          ! calculation of spin density
                          rho_up(1:vl)=&
                               (phi_ur_phi_ur + &
                               phi_ui_phi_ui ) * occ_real
                          rho_dn(1:vl)=&
                               (phi_dr_phi_dr  + &
                               phi_di_phi_di ) * occ_real
                          rho(1:vl,1)    =  rho(1:vl,1)    + rho_up(1:vl) + rho_dn(1:vl)
                          s_dens(1:vl,3) =  s_dens(1:vl,3) + rho_up(1:vl) - rho_dn(1:vl)
                          s_dens(1:vl,1) =  s_dens(1:vl,1) + (phi_ur_phi_dr  + &
                               phi_ui_phi_di )*occ_real*2
                          s_dens(1:vl,2) =  s_dens(1:vl,2) + (phi_ur_phi_di  -&
                               phi_ui_phi_dr )*occ_real*2
                       end do! partners
                    end do! orbitals
                 end if
                 end associate
              end do ! irreps

           else ! options_spin_orbit

             do i=1,n_irrep

                associate( occ=>occ_num_occ(i)%m )
                eig_dim=size(eigvec_occ(i)%m,1)
                occ_dim=size(eigvec_occ(i)%m,2)
                if (occ_dim > 0) then
                do s=1,ispin
                   do j=1,occ_dim
                      occ_real=occ(j,s)/partners(i)
                      do l=1,partners(i)
                         phi=0.0_r8_kind
                         do m=1,eig_dim
                            do n=1,eig_dim
                               phi=phi+eigvec_occ(i)%m(m,j,s)*eigvec_occ(i)%m(n,j,s)*dens(i)%d(1:vl,n,m,l)
                            end do
                         end do
                         ! calculation of density
                         rho(1:vl,s)=rho(1:vl,s)+phi*&
                              occ_real
                      end do
                   end do
                end do
                endif
                end associate
             end do
           endif! options_spin_orbit

end subroutine density_calc_hfc
#endif

  !***************************************************************
#ifdef WITH_RESPONSE
  subroutine density_calc_response (vec_length, rho, orbs_ob, phi_ob)
    !
    ! Calculation of  the density and  of the occupied  and unoccupied
    ! orbitals needed in the response_module.
    !
    use eigen_data_module, only: eigvec
    implicit none
    integer(i4_kind), intent(in) :: vec_length
    real(r8_kind), intent(out) :: rho(:, :)
    type(orbital_type), intent(in) :: orbs_ob(:)
    type(orbital_spin_type), intent(in) :: phi_ob(:)
    !** End of interface *****************************************

#ifndef FPP_NOBLAS
    !           Blas Version

    real(kind=r8_kind)  :: occ_real
    real(kind=r8_kind), allocatable :: phi_arr(:,:)
    ! phi_arr(vec_length, number of occupied orbitals)
    real(kind=r8_kind), pointer :: ob(:,:,:) ,phi_p(:,:,:)
    integer(i4_kind) :: s, i, j, l, eig_dim, aux1

    rho(1:vec_length,:)=0.0_r8_kind

    ! NOTE:  vec_length  is  allways   smaller  or  equal  then  lower
    ! dimension of orbs_ob()%o. Rho have dimension vec_length.
    do i=1,n_irrep
       ob=> orbs_ob(i)%o
       associate( occ=>occ_num_occ(i)%m )
       eig_dim=size(eigvec(i)%m,1)
       aux1=size(eigvec_occ(i)%m,2) ! not eigv !!
       if (aux1 > 0) then ! occ. orbitals
       allocate(phi_arr(vec_length,aux1),stat=alloc_stat(33))
       ASSERT(alloc_stat(33).eq.0)
       do s=1,ispin
          phi_p=>phi_ob(i)%o(:,:,:,s)     ! unocc. orbitals for this spin
          do l=1,partners(i)
             call dgemm ('n', 'n', vec_length, aux1, eig_dim, 1.0_r8_kind, &
                  ob(1, 1, l), size(ob, 1), &
                  eigvec(i)%m(1, 1, s), eig_dim, 0.0_r8_kind, &
                  phi_arr, vec_length)
             phi_p(1:vec_length,1:aux1,l)=phi_arr(:,:)
             do j=1,aux1
                occ_real=occ(j,s)/partners(i)
                rho(1:vec_length,s)=rho(1:vec_length,s)+&
                     phi_arr(:,j)*phi_arr(:,j)*occ_real
             end do! loop over j
          end do! loop over partners
       end do
       deallocate(phi_arr,stat=alloc_stat(33))
       ASSERT(alloc_stat(33).eq.0)
       alloc_stat(33)=1
       endif ! occ. orbitals
       if (aux1 < eig_dim) then ! vir. orbitals
       allocate(phi_arr(vec_length,eig_dim-aux1),stat=alloc_stat(33))
       ASSERT(alloc_stat(33).eq.0)
       do s=1,ispin
          phi_p=>phi_ob(i)%o(:,:,:,s)     ! unocc. orbitals for this spin
          do l=1,partners(i)
!!$             call dgemm('n','n',vec_length,eig_dim-aux1,&
!!$                  eig_dim,1.0_r8_kind,ob(1,1,l), size(ob, 1),&
!!$                  eigv(1,aux1+1:eig_dim,s),&
!!$                  eig_dim,0.0_r8_kind,phi_arr,vec_length)
             call dgemm ('n', 'n', vec_length, eig_dim - aux1, eig_dim, 1.0_r8_kind, &
                  ob(1, 1, l), size(ob, 1), &
                  eigvec(i)%m(1, aux1 + 1, s), eig_dim, 0.0_r8_kind, &
                  phi_arr, vec_length)
             phi_p(1:vec_length,aux1+1:eig_dim,l)=phi_arr(:,:)
          end do! loop over partners
       end do
       deallocate(phi_arr,stat=alloc_stat(33))
       ASSERT(alloc_stat(33).eq.0)
       alloc_stat(33)=1
       endif ! vir. orbitals
       end associate
    end do
#else
#warning
#warning "You rely on slow F77-code for matrix multiplications!"
#warning


    !           f77-like version

    real(kind=r8_kind)  :: phi(vec_length),occ_real
    real(kind=r8_kind),pointer :: ob(:,:,:), phi_p(:,:,:)
    integer(i4_kind) :: s,i,j,l,m,eig_dim,counter,aux1
    rho(1:vec_length,:)=0.0_r8_kind

    do i=1,n_irrep
       ob=> orbs_ob(i)%o
       associate( occ=>occ_num_occ(i)%m )
       eig_dim=size(eigvec(i)%m,1)
       aux1=size(eigvec_occ(i)%m,2) ! not eigv !!
       do s=1,ispin
          phi_p=>phi_ob(i)%o(:,:,:,s)
          if (aux1 > 0) then ! occ. orbitals
          do j=1,aux1
             occ_real=occ(j,s)/partners(i)
             do l=1,partners(i)
                phi=0.0_r8_kind
                do m=1,eig_dim
                   phi=phi+eigvec(i)%m(m,j,s)&
                        *ob(1:vec_length,m,l)
                end do
                phi_p(1:vec_length,j,l)=phi(1:vec_length)
                ! calculation of density
                rho(1:vec_length,s)=rho(1:vec_length,s)+phi*phi*&
                     occ_real
             end do
          end do
          endif ! occ. orbitals
          counter=1
          if (aux1 < eig_dim) then ! vir. orbitals
          do j=aux1+1,eig_dim
             do l=1,partners(i)
                phi=0.0_r8_kind
                do m=1,eig_dim
                   phi=phi+eigvec(i)%m(m,j,s)&
                        *ob(1:vec_length,m,l)
                end do
                phi_p(1:vec_length,j,l)=phi(1:vec_length)
             end do
          enddo
          endif ! vir. orbitals
       end do
       end associate
    end do
#endif

  end subroutine density_calc_response
#endif

  subroutine density_calc_nl_v2 (vl, rho, gamma, grarho, &
       orbs_spinor_ob, orbs_spinor_grads)
    !
    ! Calculation of the density.  Also the gradient of the density is
    ! calculated.
    !
    use spin_orbit_module, only: spin_orbit_polarized
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:,:), gamma(:,:), grarho(:,:,:) ! (vl,3,ispin)
    type(spinor_type), intent(in), optional :: orbs_spinor_ob(:)
    type(spinor_gradient_type), intent(in), optional :: orbs_spinor_grads(:)
    !** End of interface *****************************************

    real(r8_kind), allocatable :: phi(:,:,:,:,:) !(vl,n_ev,0:4,UP:DN,RE:IM)
    real(r8_kind), pointer :: o(:,:,:,:,:), g(:,:,:,:,:,:)

    real(r8_kind)    :: occ_real

    integer(i4_kind) :: irr, dim_irr, n_ev, p, i, j, xyz, ovl
    integer(i4_kind) :: ldo
    integer(i4_kind) :: memstat

    DPRINT 'dcm/density_calc_nl_v2: entered'
    if(spin_orbit_polarized)then
       call error_handler("dcm/density_calc_nl_v2: polarized version not yet supported")
    endif

    rho    = zero
    grarho = zero

    do irr=1,size(orbs_spinor_ob) !==n_irr
       associate( eigv_real => eigvec_occ_real(irr)%m                & ! (dims(irr),occ_ev_dim)
                , eigv_imag => eigvec_occ_imag(irr)%m                &
                , occ => occ_num_occ(irr)%m )

       dim_irr   =  size(eigv_real,1)
       n_ev      =  size(eigv_real,2)
       if (n_ev.eq.0) cycle ! next irr

       o => orbs_spinor_ob(irr)%o
       g => orbs_spinor_grads(irr)%o

       ldo = size(o, 1)
       ASSERT(ldo==size(g,1))

       allocate(phi(vl,0:3,n_ev,UP:DN,RE:IM),stat=memstat)
       ASSERT(memstat==0)

       do p=1,partners(irr)

          ! phi(up,re) = o(up,re)*e(re) - o(up,im)*e(im):
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,RE,UP)      , ldo,    &
                       &       eigv_real(:,:)      , dim_irr,&
                       & zero, phi(1,VALUE,1,UP,RE), 4*vl)
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & -one, o(:,:,p,IM,UP)      , ldo,    &
                       &       eigv_imag(:,:)      , dim_irr,&
                       &  one, phi(1,VALUE,1,UP,RE), 4*vl)

          ! phi(up,im) = o(up,im)*e(re) + o(up,im)*e(re):
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,IM,UP)      , ldo,    &
                       &       eigv_real(:,:)      , dim_irr,&
                       & zero, phi(1,VALUE,1,UP,IM), 4*vl)
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,RE,UP)      , ldo,    &
                       &       eigv_imag(:,:)      , dim_irr,&
                       &  one, phi(1,VALUE,1,UP,IM), 4*vl)

          ! phi(dn,re) = o(dn,re)*e(re) - o(dn,im)*e(im):
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,RE,DN)      , ldo,    &
                       &       eigv_real(:,:)      , dim_irr,&
                       & zero, phi(1,VALUE,1,DN,RE), 4*vl)
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & -one, o(:,:,p,IM,DN)      , ldo,    &
                       &       eigv_imag(:,:)      , dim_irr,&
                       &  one, phi(1,VALUE,1,DN,RE), 4*vl)

          ! phi(dn,im) = o(dn,im)*e(re) + o(dn,im)*e(re):
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,IM,DN)      , ldo,    &
                       &       eigv_real(:,:)      , dim_irr,&
                       & zero, phi(1,VALUE,1,DN,IM), 4*vl)
          call dgemm('n','n',vl,n_ev,dim_irr,&
                       & +one, o(:,:,p,RE,DN)      , ldo,    &
                       &       eigv_imag(:,:)      , dim_irr,&
                       &  one, phi(1,VALUE,1,DN,IM), 4*vl)

          if( vl .eq. ldo )then
             ovl = vl * 3
          else
             ovl = vl
          endif
          do XYZ=X,Z
             ! dphi(up,re) = o(up,re)*e(re) - o(up,im)*e(im):
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,RE,UP)       , 3*ldo,  &
                  &       eigv_real(:,:)           , dim_irr,&
                  & zero, phi(1,XYZ,1,UP,RE)       , 4*vl)
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & -one, g(1,XYZ,1,p,IM,UP)       , 3*ldo,  &
                  &       eigv_imag(:,:)           , dim_irr,&
                  &  one, phi(1,XYZ,1,UP,RE)       , 4*vl)

             ! dphi(up,im) = o(up,im)*e(re) + o(up,im)*e(re):
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,IM,UP)       , 3*ldo,  &
                  &       eigv_real(:,:)           , dim_irr,&
                  & zero, phi(1,XYZ,1,UP,IM)       , 4*vl)
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,RE,UP)       , 3*ldo,  &
                  &       eigv_imag(:,:)           , dim_irr,&
                  &  one, phi(1,XYZ,1,UP,IM)       , 4*vl)

             ! dphi(dn,re) = o(dn,re)*e(re) - o(dn,im)*e(im):
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,RE,DN)       , 3*ldo,  &
                  &       eigv_real(:,:)           , dim_irr,&
                  & zero, phi(1,XYZ,1,DN,RE)       , 4*vl)
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & -one, g(1,XYZ,1,p,IM,DN)       , 3*ldo,  &
                  &       eigv_imag(:,:)           , dim_irr,&
                  &  one, phi(1,XYZ,1,DN,RE)       , 4*vl)

             ! dphi(dn,im) = o(dn,im)*e(re) + o(dn,im)*e(re):
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,IM,DN)       , 3*ldo,  &
                  &       eigv_real(:,:)           , dim_irr,&
                  & zero, phi(1,XYZ,1,DN,IM)       , 4*vl)
             call dgemm('n','n',ovl,n_ev,dim_irr,&
                  & +one, g(1,XYZ,1,p,RE,DN)       , 3*ldo,  &
                  &       eigv_imag(:,:)           , dim_irr,&
                  &  one, phi(1,XYZ,1,DN,IM)       , 4*vl)
             if (ovl .gt. vl) exit
          enddo
          do j=1,n_ev
             occ_real = occ(j, 1) / partners(irr)
             do i = 1, vl
                rho(i, 1) = rho(i, 1) + occ_real * SUM(phi(i, VALUE, j, UP:DN, RE:IM)**2)
             enddo
             do xyz = X, Z
                do i = 1, vl
                   grarho(i, xyz, 1) = grarho(i, xyz, 1) + &
                        (2 * occ_real) * SUM(phi(i, VALUE, j, UP:DN, RE:IM) * &
                                             phi(i, xyz, j, UP:DN, RE:IM))
                enddo
             enddo
          enddo
       enddo ! partner
       end associate
       deallocate(phi,stat=memstat)
       ASSERT(memstat==0)
    enddo ! irrep

    gamma(:vl,1) = SUM(grarho(:vl,X:Z,1)**2, DIM=2)
    DPRINT 'dcm/density_calc_nl_v2: exit'
  end subroutine density_calc_nl_v2

  subroutine density_calc_nl (vl, rho, gamma, grarho, orbs_ob, &
       orbs_grads, orbs_spinor_ob, orbs_spinor_grads, s_dens, s_abs, &
       gras_dens, tau)
    !
    ! Calculation  of  the  density.   Depending on  the  presense  of
    ! optional  arguments   also  the  gradient  of   the  density  is
    ! calculated.
    !
    use spin_orbit_module, only: spin_orbit_polarized
    use f77_blas, only: dgemm
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:, :)
    real(r8_kind), intent(out), optional :: gamma(:,:), grarho(:,:,:)
    type(orbital_type), intent(in), optional :: orbs_ob(:) ! (n_irrep)
    type(orbital_gradient_type), intent(in), optional :: orbs_grads(:) ! (n_irrep)
    type(spinor_type), intent(in), optional :: orbs_spinor_ob(:) ! (n_proj_irreps?)
    type(spinor_gradient_type), intent(in), optional :: orbs_spinor_grads(:) ! (n_proj_irreps?)
    real(r8_kind), intent(out), optional :: s_dens(:,:) ! (:vl,3) spin density
    real(r8_kind), intent(out), optional :: s_abs(:) ! (:vl)   == |s_dens|
    real(r8_kind), intent(out), optional :: gras_dens(:,:) ! (:vl,3) == d/dR s_abs
    real(r8_kind), intent(out), optional :: tau(:,:) ! (:vl,s) kinetic energy density for MGGA
    !** End of interface *****************************************

    integer(i4_kind) :: g,xyz,C,R

    real(kind=r8_kind)  :: &
         occ_real,occ_real_2,help_real,&
         & help_imag,&
         help_vec_ud(vl,UP:DN,RE:IM)
    real(kind=r8_kind)  :: graphi_xyz_ud(vl,X:Z,UP:DN,RE:IM)

    real(kind=r8_kind)  :: phi_ud(vl,UP:DN,RE:IM)

    real(kind=r8_kind), pointer :: ob(:,:,:), grads(:,:,:,:)
    real(kind=r8_kind),pointer :: ob_up_real(:,:,:),ob_up_imag(:,:,:),&
         ob_dn_real(:,:,:),ob_dn_imag(:,:,:)
    real(kind=r8_kind),pointer :: grads_up_real(:,:,:,:),grads_up_imag(:,:,:,:),&
         grads_dn_real(:,:,:,:),grads_dn_imag(:,:,:,:)

    integer(i4_kind) :: s, i, j, l, m, eig_dim, occ_dim
    real(r8_kind),allocatable :: phi_arr(:,:), graphi_arr(:,:,:), grad_help(:,:,:)
    !
    ! Variables for Spin Polarized Version
    real(kind=r8_kind)  :: rho_ud(vl,UP:DN) !
    real(kind=r8_kind)  :: grs_dens(vl,3,3)

    rho(1:vl,:)=0.0_r8_kind
    if(present(grarho)) grarho = 0.0_r8_kind
    if(present(tau))     tau   = ZERO

    if (options_spin_orbit) then
      !
      ! SPIN ORBIT
      !
!!$       write(1,*) 'Warning: Spin orbit with GGA not yet optimized for Linux !'
       if (spin_orbit_polarized) then
          s_dens(1:vl,:)=0.0_r8_kind
          grs_dens(1:vl,1:3,1:3)=0.0_r8_kind
          gras_dens = 0.0_r8_kind
       endif

      do i=1,n_irrep
       ob_up_real => orbs_spinor_ob(i)%o(:, :, :, 1, 1) ! spinor(1)%o
       ob_up_imag => orbs_spinor_ob(i)%o(:, :, :, 2, 1) ! spinor(1)%o_imag
       ob_dn_real => orbs_spinor_ob(i)%o(:, :, :, 1, 2) ! spinor(2)%o
       ob_dn_imag => orbs_spinor_ob(i)%o(:, :, :, 2, 2) ! spinor(2)%o_imag

       grads_up_real => orbs_spinor_grads(i)%o(:, :, :, :, 1, 1) ! spinor(1)%o
       grads_up_imag => orbs_spinor_grads(i)%o(:, :, :, :, 2, 1) ! spinor(1)%o_imag
       grads_dn_real => orbs_spinor_grads(i)%o(:, :, :, :, 1, 2) ! spinor(2)%o
       grads_dn_imag => orbs_spinor_grads(i)%o(:, :, :, :, 2, 2) ! spinor(2)%o_imag
       associate( occ => occ_num_occ(i)%m                                      &
                , eigv_real => eigvec_occ_real(i)%m                            &
                , eigv_imag => eigvec_occ_imag(i)%m )
       eig_dim = size(eigv_real, 1)
          do j=1,size(eigv_real,2)
             occ_real=occ(j,1)/partners(i)
             occ_real_2=occ_real*2.0_r8_kind
             do l=1,partners(i)
                phi_ud        = 0.0_r8_kind
                graphi_xyz_ud = 0.0_r8_kind
                ! calculation of orbital and its gradient
                do m=1,eig_dim
                   help_real=eigv_real(m,j)
                   help_imag=eigv_imag(m,j)

                   phi_ud(1:vl,UP,RE) = phi_ud(1:vl,UP,RE)&
                        & + ob_up_real(1:vl,m,l)*help_real&
                        & - ob_up_imag(1:vl,m,l)*help_imag

                   phi_ud(1:vl,UP,IM) = phi_ud(1:vl,UP,IM)&
                        & + ob_up_imag(1:vl,m,l)*eigv_real(m,j)&
                        & + ob_up_real(1:vl,m,l)*eigv_imag(m,j)

                   phi_ud(1:vl,DN,RE) =  phi_ud(1:vl,DN,RE)&
                        & + ob_dn_real(1:vl,m,l)*help_real&
                        & - ob_dn_imag(1:vl,m,l)*help_imag

                   phi_ud(1:vl,DN,IM) =  phi_ud(1:vl,DN,IM)&
                        & + ob_dn_imag(1:vl,m,l)*help_real&
                        & + ob_dn_real(1:vl,m,l)*help_imag

                   do g=X,Z
                      graphi_xyz_ud(1:vl,g,UP,RE) = graphi_xyz_ud(1:vl,g,UP,RE)&
                           & + grads_up_real(1:vl,g,m,l)*help_real&
                           & - grads_up_imag(1:vl,g,m,l)*help_imag

                      graphi_xyz_ud(1:vl,g,UP,IM) = graphi_xyz_ud(1:vl,g,UP,IM)&
                           & + grads_up_imag(1:vl,g,m,l)*help_real&
                           & + grads_up_real(1:vl,g,m,l)*help_imag

                      graphi_xyz_ud(1:vl,g,DN,RE) = graphi_xyz_ud(1:vl,g,DN,RE)&
                           & + grads_dn_real(1:vl,g,m,l)*help_real&
                           & - grads_dn_imag(1:vl,g,m,l)*help_imag

                      graphi_xyz_ud(1:vl,g,DN,IM) = graphi_xyz_ud(1:vl,g,DN,IM)&
                           & + grads_dn_imag(1:vl,g,m,l)*help_real&
                           & + grads_dn_real(1:vl,g,m,l)*help_imag
                   enddo

                end do
                if (spin_orbit_polarized) then
                   !
                   ! OPEN SHELL
                   !
                   rho_ud(1:vl,UP:DN) = occ_real * SUM(phi_ud(1:vl,UP:DN,RE:IM)**2, DIM=3)
                   rho(1:vl,1) = rho(1:vl,1) + SUM(rho_ud(1:vl,UP:DN), DIM=2)

                   s_dens(1:vl,Z) = s_dens(1:vl,Z) + (rho_ud(1:vl,UP) - rho_ud(1:vl,DN))
                   s_dens(1:vl,X) = s_dens(1:vl,X) + (&
                        &   phi_ud(1:vl,UP,RE) * phi_ud(1:vl,DN,RE)&
                        & + phi_ud(1:vl,UP,IM) * phi_ud(1:vl,DN,IM)&
                        &   ) * occ_real * 2
                   s_dens(1:vl,Y) =  s_dens(1:vl,Y) + (&
                        &   phi_ud(1:vl,UP,RE) * phi_ud(1:vl,DN,IM)&
                        & - phi_ud(1:vl,UP,IM) * phi_ud(1:vl,DN,RE)&
                        &   ) * occ_real * 2

                   ! gradients of rho
                   do g=X,Z
                      do C=UP,DN
                         do R=RE,IM
                            grarho(1:vl,g,1) = grarho(1:vl,g,1)&
                                 & + ( phi_ud(1:vl,C,R) * graphi_xyz_ud(1:vl,g,C,R) )&
                                 &   * occ_real_2
                         enddo
                      enddo
                   enddo

                   ! gradients of sz_dens
!!$                   call error_handler("dcm/density_calc_nl: need modification here ...")

!!$                   grs_dens(1:vl,3,1)=grs_dens(1:vl,3,1)+&
!!$                        (phi_up_real(:)*graphix_up_real(:)+&
!!$                        phi_up_imag(:)*graphix_up_imag(:)-&
!!$                        phi_down_real(:)*graphix_down_real(:)-&
!!$                        phi_down_imag(:)*graphix_down_imag(:))*occ_real_2
!!$                   grs_dens(1:vl,3,2)=grs_dens(1:vl,3,2)+&
!!$                        (phi_up_real(:)*graphiy_up_real(:)+&
!!$                        phi_up_imag(:)*graphiy_up_imag(:)-&
!!$                        phi_down_real(:)*graphiy_down_real(:)-&
!!$                        phi_down_imag(:)*graphiy_down_imag(:))*occ_real_2
!!$                   grs_dens(1:vl,3,3)=grs_dens(1:vl,3,3)+&
!!$                        (phi_up_real(:)*graphiz_up_real(:)+&
!!$                        phi_up_imag(:)*graphiz_up_imag(:)-&
!!$                        phi_down_real(:)*graphiz_down_real(:)-&
!!$                        phi_down_imag(:)*graphiz_down_imag(:))*occ_real_2

                   do g=X,Z ! grads of SZ
                      grs_dens(:vl,Z,g) = grs_dens(:vl,Z,g) + (&
                           &   phi_ud(:vl,UP,RE) * graphi_xyz_ud(:vl,g,UP,RE)&
                           & + phi_ud(:vl,UP,IM) * graphi_xyz_ud(:vl,g,UP,IM)&
                           & - phi_ud(:vl,DN,RE) * graphi_xyz_ud(:vl,g,DN,RE)&
                           & - phi_ud(:vl,DN,IM) * graphi_xyz_ud(:vl,g,DN,IM)&
                           & ) * occ_real_2
                   enddo

                   ! gradients of sx_dens
!!$                   grs_dens(1:vl,1,1) =  grs_dens(1:vl,1,1) + &
!!$                        (graphix_up_real(:)* phi_down_real(:) +&
!!$                        graphix_up_imag(:)* phi_down_imag(:)+ &
!!$                        phi_up_real(:)* graphix_down_real(:) +&
!!$                        phi_up_imag(:)* graphix_down_imag(:))*occ_real_2
!!$                   grs_dens(1:vl,1,2) =  grs_dens(1:vl,1,2) + &
!!$                        (graphiy_up_real(:)* phi_down_real(:) +&
!!$                        graphiy_up_imag(:)* phi_down_imag(:)+ &
!!$                        phi_up_real(:)* graphiy_down_real(:) +&
!!$                        phi_up_imag(:)* graphiy_down_imag(:))*occ_real_2
!!$                   grs_dens(1:vl,1,3) =  grs_dens(1:vl,1,3) + &
!!$                        (graphiz_up_real(:)* phi_down_real(:) +&
!!$                        graphiz_up_imag(:)* phi_down_imag(:)+ &
!!$                        phi_up_real(:)* graphiz_down_real(:) +&
!!$                        phi_up_imag(:)* graphiz_down_imag(:))*occ_real_2
                   do g=X,Z ! grads of SX
                      grs_dens(:vl,X,g) = grs_dens(:vl,X,g) + (&
                           &   phi_ud(:vl,DN,RE) * graphi_xyz_ud(:vl,g,UP,RE)&
                           & + phi_ud(:vl,DN,IM) * graphi_xyz_ud(:vl,g,UP,IM)&
                           & + phi_ud(:vl,UP,RE) * graphi_xyz_ud(:vl,g,DN,RE)&
                           & + phi_ud(:vl,UP,IM) * graphi_xyz_ud(:vl,g,DN,IM)&
                           & ) * occ_real_2
                   enddo
!!$
!!$                   ! gradients of sy_dens
!!$                   grs_dens(1:vl,2,1) =  grs_dens(1:vl,2,1) + &
!!$                        (graphix_up_real(:)* phi_down_imag(:) -&
!!$                        graphix_up_imag(:)* phi_down_real(:) +&
!!$                        phi_up_real(:)* graphix_down_imag(:) -&
!!$                        phi_up_imag(:)* graphix_down_real(:))*occ_real_2
!!$                   grs_dens(1:vl,2,2) =  grs_dens(1:vl,2,2) + &
!!$                        (graphiy_up_real(:)* phi_down_imag(:) -&
!!$                        graphiy_up_imag(:)* phi_down_real(:) +&
!!$                        phi_up_real(:)* graphiy_down_imag(:) -&
!!$                        phi_up_imag(:)* graphiy_down_real(:))*occ_real_2
!!$                   grs_dens(1:vl,2,3) =  grs_dens(1:vl,2,3) + &
!!$                        (graphiz_up_real(:)* phi_down_imag(:) -&
!!$                        graphiz_up_imag(:)* phi_down_real(:) +&
!!$                        phi_up_real(:)* graphiz_down_imag(:) -&
!!$                        phi_up_imag(:)* graphiz_down_real(:))*occ_real_2
                   do g=X,Z ! grads of SY
                      grs_dens(:vl,Y,g) = grs_dens(:vl,Y,g) + (&
                           &   phi_ud(:vl,DN,IM) * graphi_xyz_ud(:vl,g,UP,RE)&
                           & - phi_ud(:vl,DN,RE) * graphi_xyz_ud(:vl,g,UP,IM)&
                           & + phi_ud(:vl,UP,RE) * graphi_xyz_ud(:vl,g,DN,IM)&
                           & - phi_ud(:vl,UP,IM) * graphi_xyz_ud(:vl,g,DN,RE)&
                           & ) * occ_real_2
                   enddo
                else
                   !
                   ! CLOSED SHELL
                   !
                   rho(1:vl,1) = rho(1:vl,1) + occ_real&
                        & * ( SUM(SUM(phi_ud(1:vl,UP:DN,RE:IM)**2, DIM=3), DIM=2) )

                   ! calculation of the gradient of the density
                   help_vec_ud(1:vl,UP:DN,RE:IM) = phi_ud(1:vl,UP:DN,RE:IM)*occ_real_2

                   do g=X,Z
                      grarho(1:vl,g,1) = grarho(1:vl,g,1) + &
                           & SUM(&
                           &     SUM(&
                           &         graphi_xyz_ud(1:vl,g,UP:DN,RE:IM)*help_vec_ud(1:vl,UP:DN,RE:IM),&
                           &         DIM=3),&
                           &     DIM=2&
                           &    )
                   enddo
                endif
             enddo
          end do
          end associate
      end do
      if (spin_orbit_polarized) then
         !
         ! OPEN SHELL
         !

         ! now calculate absolute value of density and its gradients:
         s_abs(:vl)=sqrt( s_dens(:vl,X)**2 &
                        + s_dens(:vl,Y)**2 &
                        + s_dens(:vl,Z)**2 &
                        )
         do g=X,Z ! gradients
            do xyz=X,Z ! components of S
               gras_dens(1:vl,g) = gras_dens(1:vl,g) &
                 & + grs_dens(1:vl,xyz,g)*s_dens(1:vl,xyz)
            enddo

            ! safe division by length of polarization vector:
            WHERE( s_abs(:vl) /= 0.0_r8_kind )
              gras_dens(:vl,g) = gras_dens(:vl,g) / s_abs(:vl)
            ELSEWHERE
              gras_dens(:vl,g) = 0.0_r8_kind
            ENDWHERE
         enddo
      endif

    else ! options_spin_orbit
      !
      ! STANDARD SCF (NO SPIN ORBIT)
      !
        if(.not.present(gamma)) then
              do i=1,n_irrep
                 ob => orbs_ob(i)%o
                 associate( occ => occ_num_occ(i)%m ) ! TODO: maybe occ_num_occ as intent(in)?
                 eig_dim = size(eigvec_occ(i)%m, 1)
                 occ_dim = size(eigvec_occ(i)%m, 2)

                 if (occ_dim > 0) then
                 allocate(phi_arr(vl,occ_dim),stat=alloc_stat(7))
                 if(alloc_stat(7)/=0) call error_handler&
                      ('density_calc: allocation phi_arr failed')

                 do s=1,ispin
                    do l=1,partners(i)
                       call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                            ob(:, :, l), size(ob, 1), &
                            eigvec_occ(i)%m(:, :, s), eig_dim, 0.0_r8_kind, &
                            phi_arr, vl)
                       do j=1,occ_dim
                          occ_real=occ(j,s)/partners(i)
                          rho(1:vl,s)=rho(1:vl,s)+&
                               phi_arr(:,j)*phi_arr(:,j)*occ_real
                       end do! loop over j
                    end do! loop over partners
                 end do
                 deallocate(phi_arr,stat=alloc_stat(7))
                 ASSERT(alloc_stat(7).eq.0)
                 alloc_stat(7)=1
                 endif
                 end associate
              end do

      else ! present gamma
      do i=1,n_irrep
       ob => orbs_ob(i)%o
       associate( occ => occ_num_occ(i)%m )
       grads => orbs_grads(i)%o

       eig_dim = size(eigvec_occ(i)%m, 1)
       occ_dim = size(eigvec_occ(i)%m, 2)

       if (occ_dim > 0) then
       allocate(phi_arr(vl,occ_dim),&
                graphi_arr(vl,3,occ_dim),stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(34)=0
       MEMLOG(size(phi_arr)+size(graphi_arr))

       do s=1,ispin
             do l=1,partners(i)
                ! calculation of orbital

                call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                     ob(:, :, l), size(ob, 1), &
                     eigvec_occ(i)%m(:, :, s), eig_dim, 0.0_r8_kind, &
                     phi_arr, vl)

                ! and gradient
                if (vl == size(grads, 1)) then

                   call dgemm ('n', 'n', vl*3, occ_dim, eig_dim, 1.0_r8_kind, &
                        grads(:, :, :, l), 3 * size(grads, 1), &
                        eigvec_occ(i)%m(:, :, s), eig_dim, 0.0_r8_kind, &
                        graphi_arr(:,:,:), vl*3)
                else
                   allocate(grad_help(vl,3,eig_dim),stat=alloc_stat(9))
                   ASSERT(alloc_stat(9).eq.0)
                   MEMLOG(size(grad_help))

                   grad_help(1:vl,:,:)=grads(1:vl,:,:,l)

                   call dgemm('n','n',vl*3,occ_dim,&
                        eig_dim,1.0_r8_kind,grad_help(:,:,:),vl*3,eigvec_occ(i)%m(:,:,s),&
                        eig_dim,0.0_r8_kind,graphi_arr(:,:,:),vl*3)

                   MEMLOG(-size(grad_help))
                   deallocate(grad_help,stat=alloc_stat(9))
                   ASSERT(alloc_stat(9).eq.0)
                   alloc_stat(9)=1
                endif

                ! sum up the squares of the orbitals into density,
                ! also compute the density gradient:
                do j=1,occ_dim

                   occ_real = occ(j,s) / partners(i)

                   ! Calculation of the density of the orbital
                   rho(1:vl,s) = rho(1:vl,s) + phi_arr(:,j) * phi_arr(:,j) * occ_real

                   do g=1,3
                      grarho(1:vl,g,s) = grarho(1:vl,g,s) &
                                       + 2 * phi_arr(:,j) * graphi_arr(:,g,j) * occ_real
                   enddo

                   !----------------------------------------------------------------------------------+
                   ! Calculation of the kinetic energy density from the gradients of the orbitals     |
                   ! Note that the normal definition is used here:                                    |
                   !                                                                                  |
                   !                     tau = 0.5*sum_{i}|\nabla\phi_{i}|^2                          |
                   !                                                                                  |
                   ! MAY BE CHANGED LATER ON                                                          |
                   !----------------------------------------------------------------------------------+
                   if (present(tau)) then
                      tau(1:vl,s) = tau(1:vl,s) &
                          + sum(graphi_arr(1:vl,:,j) * graphi_arr(1:vl,:,j),DIM=2) * ONEHALF * occ_real
                   endif
                enddo
          end do ! partner
       end do ! spin
       MEMLOG(-size(phi_arr)-size(graphi_arr))
       deallocate(phi_arr,graphi_arr,stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(8)=1
       alloc_stat(34)=1 !graphi_arr
       endif
       end associate
      end do !irrep
     endif ! gamma
    endif ! options_spin_orbit

    if(present(gamma)) then
       if(ispin==1) then
          if (.not.spin_orbit_polarized) then
             gamma(1:vl,1) = SUM(grarho(1:vl,X:Z,1)**2, DIM=2)
          endif
       else
          gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)
          gamma(1:vl,DNDN) = SUM(grarho(1:vl,X:Z,DN)**2, DIM=2)
          gamma(1:vl,UPDN) = SUM(grarho(1:vl,X:Z,UP)*grarho(1:vl,X:Z,DN), DIM=2)
       endif
    endif

  end subroutine density_calc_nl


  subroutine fitted_density_calc(vec_length,&
                                 rho,grad,sec_der,nuc_grad,sec_nuc_der,&
                                 fcts,grads,sec_ders,nuc_grads,sec_nuc_ders,&
                                 altern_coeff)
    ! Purpose : Calculation of the fitted density (and its derivatives)
    !           rho(r,s) = Sum(k=1,N) a_k,s f_k(r)
    !           Input:  a_k,tot  = coeff_charge(k)
    !                   a_k,spin = coeff_spin(k)
    !                   a_k,up   = [ a_k,tot + a_k,spin ] / 2
    !                   a_k,down = [ a_k,tot - a_k,spin ] / 2
    !           Output: rho(r,tot)  = rho(:,1)   if spin-restricted
    !                   rho(r,up)   = rho(:,1)   if spin-polarized
    !                   rho(r,down) = rho(:,2)   if spin_polarized
    use fit_coeff_module
    use orbitalprojection_module, only: orbitalprojection_ch, &
                                        orbitalprojection_globcontr_ch
    use time_module, only: start_timer, stop_timer
    use timer_module

    integer(kind=i4_kind),intent(in) :: vec_length
    ! output arrays
    real(kind=r8_kind), optional,intent(out) :: rho (:,:)
    real(kind=r8_kind), optional,intent(out) :: grad(:,:,:)
    real(kind=r8_kind), optional,intent(out) :: sec_der(:,:,:)
    type(arrmat4)     , optional, target :: nuc_grad(:)
    type(arrmat5)     , optional, target :: sec_nuc_der(:)
    ! input arrays
    type(orbital_type)                 , optional, target :: fcts
    type(orbital_gradient_type)        , optional, target :: grads
    type(orbital_sec_der_type)         , optional, target :: sec_ders
    type(orbital_nuclear_gradient_type), optional, target :: nuc_grads(:)
    type(orbital_nuclear_sec_der_type) , optional, target :: sec_nuc_ders(:)
    real(kind=r8_kind), optional         :: altern_coeff(:,:)

    !** End of interface *****************************************
    integer(i4_kind)              :: vl, s, i, j, i_ua, i_ma, i_ea
    logical                       :: spin, elec, nuc
    real(kind=r8_kind)            :: signum, coeff
    real(kind=r8_kind),allocatable:: coeffv(:)
    real(kind=r8_kind), parameter :: half = 0.5_r8_kind
    real(kind=r8_kind), pointer   :: gs(:,:), snds(:,:), &
                                            g (:,:), sd (:,:,:)
    if ( ( present(rho        ) .neqv. present(fcts        ) ) .or. &
         ( present(grad       ) .neqv. present(grads       ) ) .or. &
         ( present(sec_der    ) .neqv. present(sec_ders    ) ) .or. &
         ( present(nuc_grad   ) .neqv. present(nuc_grads   ) ) .or. &
         ( present(sec_nuc_der) .neqv. present(sec_nuc_ders) ) ) &
         call error_handler("fitted_density_calc: illegal set of arguments")

        allocate(coeffv(fit_coeff_n_ch()))

    if(present(rho)) rho(:,:)=0.0_r8_kind

!AG #ifdef FPP_FAST_COMPILE
!AG     call error_handler("dcm/fitted_density_calc: recompile without FPP_FAST_COMPILE")
!AG #else

    call start_timer(timer_int_fit_density_calc)
    vl   = vec_length
    if(present(grad) ) then
        DPRINT  sum(grads%o(:,:,:,1)),sum(grad(:vl,:,1)),present(altern_coeff)
    endif
    DPRINT shape(fcts%o(:,:,1)),vl,present(sec_der)
    DPRINT sum(fcts%o(:,:,1))
    spin = ispin > 1
    elec = present(rho) .or. present(grad) .or. present(sec_der)
    nuc  = present(nuc_grad) .or. present(sec_nuc_der)
    do s=1,ispin
       signum = real(3-2*s,r8_kind)
       if (elec) then
             if (spin) then
                coeffv = half*( coeff_charge(:) + signum * coeff_spin(:) )
             else
                coeffv = coeff_charge(:)
             endif

             if(present(altern_coeff)) coeffv=altern_coeff(:,s)

             ! compute the density and its electronic derivatives
             if (present(rho)) then
               do i=1,fit_coeff_n_ch()
                rho(:vl,s) = rho(:vl,s) + coeffv(i) * fcts%o(:vl,i,1)
               enddo
             endif

             if (present(grad)) then
          do i=1,fit_coeff_n_ch()
                 grad(:vl,:,s) = grad(:vl,:,s) + coeffv(i) * grads%o(:vl,:,i,1)
          enddo
             endif
             if (present(sec_der)) then
          do i=1,fit_coeff_n_ch()
                sec_der(:vl,:,s) = sec_der(:vl,:,s) + coeffv(i) * sec_ders%o(:vl,:,i,1)
          enddo
             endif
       endif

       if (nuc) then
          do i_ma=1,N_moving_unique_atoms
             i_ua = moving_unique_atom_index(i_ma)
             i = orbitalprojection_ch(-1,i_ua)
             do j=1,dims_fit(i_ma)+unique_atoms(i_ua)%N_glob_cons_ch
                if (spin) then
                   coeff = half*( coeff_charge(i) + signum * coeff_spin(i) )
                else
                   coeff = coeff_charge(i)
                endif
                do i_ea=1,unique_atoms(i_ua)%N_equal_atoms
                   if (present(nuc_grad)) then
                      gs => nuc_grads(i_ma)%o(:,:,i_ea,j,1)
                      g  => nuc_grad (i_ua)%m(:,:,i_ea,s)
                      g(:vl,:) = g(:vl,:) + coeff * gs(:vl,:)
                   endif
                   if (present(sec_nuc_der)) then
                      snds => sec_nuc_ders(i_ma)%o(:,:,i_ea,j,1) ! 1:vl,1:6
                      sd   => sec_nuc_der (i_ua)%m(:,:,:,i_ea,s) ! 1:vl,1:3,1:3
                      sd(:vl,1,1) = sd(:vl,1,1) + coeff * snds(:vl,1)
                      sd(:vl,1,2) = sd(:vl,1,2) + coeff * snds(:vl,2)
                      sd(:vl,1,3) = sd(:vl,1,3) + coeff * snds(:vl,3)
                      sd(:vl,2,2) = sd(:vl,2,2) + coeff * snds(:vl,4)
                      sd(:vl,2,3) = sd(:vl,2,3) + coeff * snds(:vl,5)
                      sd(:vl,3,3) = sd(:vl,3,3) + coeff * snds(:vl,6)
                   endif
                enddo
                i = i + 1
                if (j == dims_fit(i_ma)) i = orbitalprojection_globcontr_ch(i_ua)
             enddo
             if (present(sec_nuc_der)) then
                do i_ea=1,unique_atoms(i_ua)%N_equal_atoms
                   sd => sec_nuc_der(i_ua)%m(:,:,:,i_ea,s)
                   sd(:vl,2,1) = sd(:vl,1,2)
                   sd(:vl,3,1) = sd(:vl,1,3)
                   sd(:vl,3,2) = sd(:vl,2,3)
                enddo
             endif
          enddo ! i_ma
       endif
    end do ! s
    call stop_timer(timer_int_fit_density_calc)

        deallocate(coeffv)
  end subroutine fitted_density_calc

  subroutine fitted_core_density_calc(vec_length,&
                      rho,grad,sec_der,nuc_grad,sec_nuc_der,&
                      fcts,grads,sec_ders,nuc_grads,sec_nuc_ders)
    ! Purpose : Calculation of the fitted core density (and its derivatives)
    !           rho(r,s) = Sum(k=1,N) a_k,s f_k(r)
    !           Input:  a_k,tot  = coeff_charge(k)
    !                   a_k,spin = coeff_spin(k)
    !                   a_k,up   = [ a_k,tot + a_k,spin ] / 2
    !                   a_k,down = [ a_k,tot - a_k,spin ] / 2
    !           Output: rho(r,tot)  = rho(:,1)   if spin-restricted
    !                   rho(r,up)   = rho(:,1)   if spin-polarized
    !                   rho(r,down) = rho(:,2)   if spin_polarized
    use fit_coeff_module, only: fit_coeff_n_cd
    use orbitalprojection_module, only: orbitalprojection_cd
    integer(kind=i4_kind),intent(in) :: vec_length
    ! output arrays
    real(kind=r8_kind), optional         :: rho (:,:)
    real(kind=r8_kind), optional         :: grad(:,:,:)
    real(kind=r8_kind), optional         :: sec_der(:,:,:)
    type(arrmat4)     , optional, target :: nuc_grad(:)
    type(arrmat5)     , optional, target :: sec_nuc_der(:)
    ! input arrays
    type(core_orbital_type)                 , optional, target :: fcts
    type(core_orbital_gradient_type)        , optional, target :: grads
    type(core_orbital_sec_der_type)         , optional, target :: sec_ders
    type(core_orbital_nuc_gradient_type), optional, target :: nuc_grads(:)
    type(core_orbital_nuc_sec_der_type) , optional, target :: sec_nuc_ders(:)
    !** End of interface *****************************************
    integer(i4_kind)              :: vl, s, i, j, i_ua, i_ma, i_ea
    logical                       :: spin, elec, nuc
    real(kind=r8_kind)            :: coeff_ispin, coeff
    real(kind=r8_kind), parameter :: half = 0.5_r8_kind
    real(kind=r8_kind), pointer   :: fs(:), gs(:,:), sds(:,:), snds(:,:), &
                                            g (:,:), sd (:,:,:)

    if ( ( present(rho        ) .neqv. present(fcts        ) ) .or. &
         ( present(grad       ) .neqv. present(grads       ) ) .or. &
         ( present(sec_der    ) .neqv. present(sec_ders    ) ) .or. &
         ( present(nuc_grad   ) .neqv. present(nuc_grads   ) ) .or. &
         ( present(sec_nuc_der) .neqv. present(sec_nuc_ders) ) ) &
         call error_handler("fitted_density_calc: illegal set of arguments")

    vl   = vec_length
    spin = ispin > 1
    elec = present(rho) .or. present(grad) .or. present(sec_der)
    nuc  = present(nuc_grad) .or. present(sec_nuc_der)

    coeff_ispin = 0.0_r8_kind /real(ispin, r8_kind) ! HHHHHHH

    do s=1,ispin
       if (elec) then
          do i = 1, fit_coeff_n_cd()
             ! compute the density and its electronic derivatives
                coeff = coeff_core(i)*coeff_ispin
                if (present(rho)) then
                   fs => fcts%o(:,i)
                   rho(:vl,s) = rho(:vl,s) + coeff * fs(:vl)
                endif
                if (present(grad)) then
                   gs => grads%o(:,:,i)
                   grad(:vl,:,s) = grad(:vl,:,s) + coeff * gs(:vl,:)
                endif
                if (present(sec_der)) then
                   sds => sec_ders%o(:,:,i)
                   sec_der(:vl,:,s) = sec_der(:vl,:,s) + coeff * sds(:vl,:)
                endif
          end do
       endif
       if (nuc) then
          do i_ma=1,N_moving_unique_atoms
             i_ua = moving_unique_atom_index(i_ma)
             i = orbitalprojection_cd(-1,i_ua)
             do j=1,dims_core(i_ma)
                   coeff = coeff_core(i)*coeff_ispin
                do i_ea=1,unique_atoms(i_ua)%N_equal_atoms
                   if (present(nuc_grad)) then
                      gs => nuc_grads(i_ma)%o(:,:,i_ea,j)
                      g  => nuc_grad (i_ua)%m(:,:,i_ea,s)
                      g(:vl,:) = g(:vl,:) + coeff * gs(:vl,:)
                   endif
                   if (present(sec_nuc_der)) then
                      snds => sec_nuc_ders(i_ma)%o(:,:,i_ea,j)   ! 1:vl,1:6
                      sd   => sec_nuc_der (i_ua)%m(:,:,:,i_ea,s) ! 1:vl,1:3,1:3
                      sd(:vl,1,1) = sd(:vl,1,1) + coeff * snds(:vl,1)
                      sd(:vl,1,2) = sd(:vl,1,2) + coeff * snds(:vl,2)
                      sd(:vl,1,3) = sd(:vl,1,3) + coeff * snds(:vl,3)
                      sd(:vl,2,2) = sd(:vl,2,2) + coeff * snds(:vl,4)
                      sd(:vl,2,3) = sd(:vl,2,3) + coeff * snds(:vl,5)
                      sd(:vl,3,3) = sd(:vl,3,3) + coeff * snds(:vl,6)
                   endif
                enddo
                i = i + 1
             enddo
             if (present(sec_nuc_der)) then
                do i_ea=1,unique_atoms(i_ua)%N_equal_atoms
                   sd => sec_nuc_der(i_ua)%m(:,:,:,i_ea,s)
                   sd(:vl,2,1) = sd(:vl,1,2)
                   sd(:vl,3,1) = sd(:vl,1,3)
                   sd(:vl,3,2) = sd(:vl,2,3)
                enddo
             endif
          enddo ! i_ma
       endif
    end do ! s

  end subroutine fitted_core_density_calc

  subroutine density_calc_ph (vl, rho, gamma, grarho, &
       nuc_grarho, orbs_ob, orbs_grads, nuc_grads)
    !
    ! Calculation of the density. Also the gradient of the density is
    ! calculated.
    !
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:,:), gamma(:,:), grarho(:,:,:)
    type(arrmat4), allocatable, target :: nuc_grarho(:) ! gradient of rho with
    ! respect to nuclear displacements
    type(orbital_type), intent(in) :: orbs_ob(:) ! orbitals
    type(orbital_gradient_type), intent(in) :: orbs_grads(:) ! gradients of orbitals
    type(orbital_nuclear_gradient_type), pointer :: nuc_grads(:,:)
    !** End of interface *****************************************

    real(kind=r8_kind)  :: &
         occ_real,occ_real_2,help_vec(vl)
    real(kind=r8_kind),pointer  :: ob(:,:,:),&
         grads(:,:,:,:),nuc(:,:,:,:,:),nuc_grarho_p(:,:,:,:)
    integer(i4_kind) :: s, i, j, l, eig_dim, counter, i_ua, i_ea, occ_dim, orb_dim
    real(kind=r8_kind), allocatable :: graphi_arr(:,:,:), phi_arr(:,:), &
         grad_help(:,:,:), nuc_help(:,:,:,:), help_nuc_arr(:,:,:,:)

    integer(i4_kind) :: xyz

    rho(1:vl,:)=0.0_r8_kind
    grarho = ZERO
    do i_ua=1,n_unique_atoms
       nuc_grarho(i_ua)%m(1:vl,:,:,:)=0.0_r8_kind
    end do

    irr: do i=1,n_irrep
       ob=> orbs_ob(i)%o
       associate( occ=>occ_num_occ(i)%m )
       grads=>orbs_grads(i)%o
       eig_dim=size(eigvec_occ(i)%m,1)
       occ_dim = size(eigvec_occ(i)%m,2)
       if (occ_dim > 0) then
       allocate(phi_arr(vl,occ_dim),&
             graphi_arr(vl,3,occ_dim),stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(34)=0 ! graphi_arr
       MEMLOG(size(phi_arr)+size(graphi_arr))
       spin: do s=1,ispin
         part: do l=1,partners(i)
             ! calculation of orbital

             call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                  ob(1, 1, l), size(ob, 1), &
                  eigvec_occ(i)%m(1, 1, s), eig_dim, 0.0_r8_kind, &
                  phi_arr, vl)

             ! and gradient
             if (vl == size(grads, 1)) then
                call dgemm ('n', 'n', vl*3, occ_dim, eig_dim, 1.0_r8_kind, &
                     grads(1, 1, 1, l), 3 * size(grads, 1), &
                     eigvec_occ(i)%m(1, 1, s), eig_dim, 0.0_r8_kind, &
                     graphi_arr(1, 1, 1), vl*3)
             else
                allocate(grad_help(vl,3,eig_dim),stat=alloc_stat(9))
                ASSERT(alloc_stat(9).eq.0)
                grad_help(1:vl,:,:)=grads(1:vl,:,:,l)
                call dgemm('n','n',vl*3,occ_dim,&
                     eig_dim,1.0_r8_kind,grad_help(1,1,1),vl*3,eigvec_occ(i)%m(1,1,s),&
                     eig_dim,0.0_r8_kind,graphi_arr(1,1,1),vl*3)
                deallocate(grad_help,stat=alloc_stat(9))
             ASSERT(alloc_stat(9).eq.0)
             alloc_stat(9)=1
             endif


             do j=1,occ_dim
                occ_real=occ(j,s)/partners(i)
                occ_real_2=occ_real*2.0_r8_kind
                ! Calculation of the density of the orbital
                rho(1:vl,s)=rho(1:vl,s)+phi_arr(1:vl,j)*phi_arr(1:vl,j)*occ_real
                ! calculation of the gradient of the density
                help_vec=phi_arr(:,j)*occ_real_2
                phi_arr(:,j)=phi_arr(:,j)*occ_real_2
                do xyz=X,Z
                   grarho(1:vl,xyz,s) = grarho(1:vl,xyz,s)&
                        & + help_vec*graphi_arr(:,xyz,j)
                enddo
             end do

             counter=0
             ! now calculate gradient with respect to nuclear displacements
           ua:  do i_ua=1,n_unique_atoms
                nuc => nuc_grads(i_ua, i)%o
                orb_dim = size(nuc, 4)
                if (orb_dim > 0) then
                allocate(help_nuc_arr(vl,3,unique_atoms(i_ua)%n_equal_atoms,occ_dim),&
                         stat=alloc_stat(10)) ! to hold orbital nuclear grads
                ASSERT(alloc_stat(10).eq.0)

                if (vl == size(nuc, 1)) then
                   ! help_nuc_arr = gradients of KS-orbitals wrt nuclear displacements
                   call dgemm ('n', 'n', vl*3*unique_atoms(i_ua)%n_equal_atoms, occ_dim, orb_dim, &
                        1.0_r8_kind, &
                        nuc(1, 1, 1, 1, l), 3 * unique_atoms(i_ua)%n_equal_atoms * size(nuc, 1),&
                        eigvec_occ(i)%m(1 + counter, 1, s), eig_dim, 0.0_r8_kind, &
                        help_nuc_arr(1, 1, 1, 1), vl*3*unique_atoms(i_ua)%n_equal_atoms)
                else
                   ! copy data first
                   allocate (nuc_help(vl, 3, unique_atoms(i_ua)%n_equal_atoms, orb_dim), &
                        stat=alloc_stat(12))
                   ASSERT(alloc_stat(12).eq.0)

                   nuc_help(1:vl, :, :, :) = nuc(1:vl, :, :, :, l)
                   call dgemm ('n', 'n', vl*3*unique_atoms(i_ua)%n_equal_atoms, occ_dim, orb_dim, &
                        1.0_r8_kind, &
                        nuc_help(1, 1, 1, 1), 3 * unique_atoms(i_ua)%n_equal_atoms * size(nuc_help, 1), &
                        eigvec_occ(i)%m(1 + counter, 1, s), eig_dim, 0.0_r8_kind, &
                        help_nuc_arr(1, 1, 1, 1), vl*3*unique_atoms(i_ua)%n_equal_atoms)
                   deallocate (nuc_help, stat=alloc_stat(12))
                   ASSERT(alloc_stat(12).eq.0)
                   alloc_stat(12) = 1
                end if

                ! help_nuc_arr=0.0_r8_kind
                nuc_grarho_p=>nuc_grarho(i_ua)%m
                do j=1,occ_dim
                   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                      nuc_grarho_p(1:vl,1,i_ea,s)=&
                           nuc_grarho_p(1:vl,1,i_ea,s)+&
                           help_nuc_arr(1:vl,1,i_ea,j)*phi_arr(:,j)
                      nuc_grarho_p(1:vl,2,i_ea,s)=&
                           nuc_grarho_p(1:vl,2,i_ea,s)+&
                           help_nuc_arr(1:vl,2,i_ea,j)*phi_arr(:,j)
                      nuc_grarho_p(1:vl,3,i_ea,s)=&
                           nuc_grarho_p(1:vl,3,i_ea,s)+&
                           help_nuc_arr(1:vl,3,i_ea,j)*phi_arr(:,j)
                   end do
                end do
                counter=counter+orb_dim
                deallocate(help_nuc_arr,stat=alloc_stat(10))
                ASSERT(alloc_stat(10).eq.0)
                alloc_stat(10)=1
                endif ! orb_dim
             enddo ua

          enddo part
       enddo spin
       MEMLOG(-size(phi_arr)-size(graphi_arr))
       deallocate(phi_arr,graphi_arr,stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(8)=1
       alloc_stat(34)=1 ! graphi_arr
       endif
       end associate
    enddo irr



    if(ispin==1) then
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)
    else
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)

       gamma(1:vl,DNDN) = SUM(grarho(1:vl,X:Z,DN)**2, DIM=2)

       gamma(1:vl,UPDN) = SUM(grarho(1:vl,X:Z,UP)*grarho(1:vl,X:Z,DN), DIM=2)
    endif

  end subroutine density_calc_ph


  subroutine density_calc_ph_nl (vl, rho, gamma, grarho, nuc_grarho, &
       orbs_ob, orbs_grads, &
       nuc_grads, nuc_sec_der, &
       nuc_3rd_der, nuc_sec_derrho, &
       nuc_dervsrho,                    & !!! dervs only, complete partial rho dervs
       nuc_dervs_grarho,                &
       graphi,                          & !!! dervs only
       tau,                             & !!! mgga mode only
       nuc_grad_tau)
    !
    ! Calculation of the density (rho). Also the gradients of the
    ! density (grarho,gamma,nuc_grarho) is calculated.
    !
    use eigen_data_module, only: eigvec
    implicit none
    integer(i4_kind), intent(in) :: vl
    real(r8_kind), intent(out) :: rho(:,:), gamma(:,:), grarho(:,:,:)
    real(r8_kind), intent(out), optional :: tau(:,:)
    !-------------------------------------------------------------------+
    ! gradient of rho with respect to nuclear displacements:            |
    !            [-d / dR] rho                                          |
    !-------------------------------------------------------------------+
    type(arrmat4), allocatable, target :: nuc_grarho(:)
    !-------------------------------------------------------------------+
    ! gradients of components of nablarho wrt nuclear displacements     |
    ! is a nonsymmetric 3x3-matrix:                                     |
    !            [-d / dR] (d rho / dr)                                 |
    !-------------------------------------------------------------------+
    ! old designation:                                                  |
    ! 2nd gradient of rho with respect to nuclear displacements         |
    !-------------------------------------------------------------------+
    type(arrmat5), allocatable, target, optional :: nuc_sec_derrho(:)
    !-------------------------------------------------------------------+
    ! gradient of tau with respect to nuclear displacements:            |
    !            [-d / dR] tau                                          |
    ! note: during the calculationa implicite usage of 3x3 matrix       |
    !-------------------------------------------------------------------+
    type(arrmat4), allocatable,target, optional :: nuc_grad_tau(:)
    !-------- to use this route in all density calcs calls
    type(orbital_type), intent(in) :: orbs_ob(:) ! orbitals
    type(orbital_gradient_type), intent(in) :: orbs_grads(:) ! gradients of
    type(orbital_nuclear_gradient_type), pointer :: nuc_grads(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_sec_der(:,:)
    type(orbital_nuclear_sec_der_type), pointer, optional :: nuc_3rd_der(:,:)
    type(arrmat6), allocatable, optional :: nuc_dervsrho(:,:) !! dervs only, partial dervs of rho
    type(arrmat7), allocatable, optional :: nuc_dervs_grarho(:,:) !! dervs only, partial dervs of rho
    !** End of interface *****************************************

    real(kind=r8_kind)  :: occ_real,occ_real_2,help_vec(vl)
    real(kind=r8_kind), allocatable :: eigv(:,:,:)
    real(kind=r8_kind),pointer  :: ob(:,:,:),&
         grads(:,:,:,:),nuc(:,:,:,:,:),sec_der(:,:,:,:,:),nuc_grarho_p(:,:,:,:),&
         nuc_sec_derrho_p(:,:,:,:,:),nuc_grad_tau_p(:,:,:,:)
    integer(i4_kind) :: i_ua2
    integer(i4_kind) :: s, i, j, l, eig_dim, counter, i_ua, i_ea, &
         occ_dim, orb_dim
    integer(i4_kind) :: i_ea2,orb_dim2,counter2,i_irr
    real(kind=r8_kind),parameter:: two=2.0_r8_kind

    real(kind=r8_kind), pointer :: graphi_arr(:,:,:)
    real(kind=r8_kind), allocatable :: phi_arr(:,:), &
         grad_help(:,:,:), nuc_help(:,:,:,:), sec_der_help(:,:,:,:), &
         help_nuc_arr(:,:,:,:), help_sec_der_arr(:,:,:,:), &
         help_3rd_der_arr(:,:,:,:)

   type(arrmat4),optional, target, intent(inout):: graphi(:,:)

   type(arrmat4),allocatable :: help_nuc_arr2(:)

    integer(i4_kind) :: kk,k
    integer(i4_kind) :: xyz
    integer(i4_kind) :: mo_dim

    real(kind=r8_kind):: dervs_grarho_help(vl,3),help_dervsrho(vl)

    rho(:vl,:)=0.0_r8_kind
    if(present(tau)) tau(:vl,:)=ZERO
    grarho = ZERO

        if(present(nuc_dervsrho)) then
         allocate(help_nuc_arr2(n_unique_atoms),stat=alloc_stat(16))
         MEMLOG(size(help_nuc_arr2))
         ASSERT(alloc_stat(16).eq.0)
        endif

    do i_ua=1,n_unique_atoms
       nuc_grarho(i_ua)%m(:vl,:,:,:)=0.0_r8_kind

       if(present(nuc_sec_derrho)) nuc_sec_derrho(i_ua)%m(:vl,:,:,:,:)=0.0_r8_kind
       if(present(nuc_grad_tau))   nuc_grad_tau(i_ua)%m(:vl,:,:,:)=0.0_r8_kind

    if(present(nuc_dervsrho)) then
     do i_ua2=1,n_unique_atoms
       nuc_dervsrho(i_ua,i_ua2)%m(:vl,:,:,:,:,:)=0.0_r8_kind
       if(present(nuc_dervs_grarho)) &
        nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,:,:,:,:,:,:)=0.0_r8_kind
     enddo
    endif

    enddo

    irr: do i=1,n_irrep
       i_irr=i
       ob=> orbs_ob(i_irr)%o
       associate( occ=>occ_num_occ(i_irr)%m )
       ! FIXME ---v
       if ( allocated( eigv ) ) deallocate( eigv )
       if(present(graphi)) then
          allocate( eigv( size( eigvec(i_irr)%m, 1 ),  size( eigvec(i_irr)%m, 2 ),  size( eigvec(i_irr)%m, 3 ) ) )
          eigv = eigvec(i_irr)%m
       else
          allocate( eigv( size( eigvec_occ(i_irr)%m, 1 ),  size( eigvec_occ(i_irr)%m, 2 ),  size( eigvec_occ(i_irr)%m, 3 ) ) )
          eigv = eigvec_occ(i_irr)%m
       endif
!!!                eigv=1.0_r8_kind !!!!!!!!! to fix coeefs

       grads=>orbs_grads(i)%o
       eig_dim=size(eigv,1)
       occ_dim=size(eigvec_occ(i_irr)%m,2)

      if (occ_dim > 0) then

      if(present(graphi)) then
       allocate(phi_arr(vl,occ_dim),stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       MEMLOG(size(phi_arr))
      else
       allocate(phi_arr(vl,occ_dim),&
                graphi_arr(vl,3,occ_dim),stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(34)=0 ! graphi_arr
       MEMLOG(size(phi_arr)+size(graphi_arr))
      endif

       spins: do s=1,ispin
          partns: do l=1,partners(i_irr)
             ! calculation of orbital
             call dgemm ('n', 'n', vl, occ_dim, eig_dim, 1.0_r8_kind, &
                  ob(1, 1, l), size(ob, 1), &
                  eigv(1, 1, s), eig_dim, 0.0_r8_kind, &
                  phi_arr, vl)

       mo_dim=occ_dim

       if(present(graphi)) then
        ! If the calling program wants them back, use the version it provided.
        ! Otherwise it has been allocated by before.
        graphi_arr=>graphi(i_irr,s)%m(:,:,:,l)
        mo_dim=eig_dim
       endif

             ! and gradient
             if (vl == size(grads, 1)) then
                call dgemm ('n', 'n', vl*3, mo_dim, eig_dim, 1.0_r8_kind, &
                     grads(1, 1, 1, l), 3 * size(grads, 1), &
                     eigv(1, 1, s), eig_dim, 0.0_r8_kind, &
                     graphi_arr(1, 1, 1), vl*3)

             else
                allocate(grad_help(vl,3,eig_dim),stat=alloc_stat(9))
                ASSERT(alloc_stat(9).eq.0)
                MEMLOG(size(grad_help))
                grad_help(:vl,:,:)=grads(:vl,:,:,l)
                call dgemm('n','n',vl*3,mo_dim,eig_dim, &
                           1.0_r8_kind,grad_help(1,1,1),vl*3, &
                           eigv(1,1,s),eig_dim, &
                           0.0_r8_kind,graphi_arr(1,1,1),vl*3)
                MEMLOG(-size(grad_help))
                deallocate(grad_help,stat=alloc_stat(9))
                ASSERT(alloc_stat(9).eq.0)
                alloc_stat(9)=1
             endif

             calc_grarho:do j=1,occ_dim
                occ_real=occ(j,s)/partners(i_irr)
                occ_real_2=occ_real*2.0_r8_kind
                ! Calculation of the density of the orbital
                rho(:vl,s)=rho(:vl,s)+phi_arr(:,j)*phi_arr(:,j)*occ_real


                ! calculation of the gradient of the density
                phi_arr(:,j)=phi_arr(:,j)*occ_real_2
                help_vec=phi_arr(:,j)

                do xyz=X,Z
                   grarho(1:vl,xyz,s)=grarho(1:vl,xyz,s)&
                        & + help_vec*graphi_arr(:,xyz,j)
                enddo

                !----------------------------------------------------------------------------------+
                ! Calculation of the kinetic energy density from the gradients of the orbitals     |
                ! Note that the normal definition is used here:                                    |
                !                                                                                  |
                !                     tau = 0.5*sum_{i}|\nabla\phi_{i}|^2                          |
                !                                                                                  |
                ! MAY BE CHANGED LATER ON                                                          |
                !----------------------------------------------------------------------------------+
                if (present(tau)) then
                   tau(1:vl,s) = tau(1:vl,s) &
                       + sum(graphi_arr(1:vl,:,j) * graphi_arr(1:vl,:,j),DIM=2) * ONEHALF * occ_real
                endif
             enddo calc_grarho

             counter=0
             ! now calculate gradient with respect to nuclear displacements
   uniques: do i_ua=1,n_unique_atoms

                nuc=> nuc_grads(i_ua,i)%o
                if(present(nuc_sec_der)) sec_der=> nuc_sec_der(i_ua,i)%o

                orb_dim=size(nuc,4)
               orb_dim1: if (orb_dim > 0) then

                allocate(help_nuc_arr(vl,3,unique_atoms(i_ua)%n_equal_atoms,occ_dim),&
                         stat=alloc_stat(10))
                         ! to hold orbital nuclear diplacement grads
                ASSERT(alloc_stat(10).eq.0)
                MEMLOG(size(help_nuc_arr))

                  if(present(nuc_dervsrho).or.present(nuc_sec_derrho)) then
                   allocate(help_sec_der_arr(vl,6,unique_atoms(i_ua)%n_equal_atoms,occ_dim),&
                            stat=alloc_stat(18))
                   ASSERT(alloc_stat(18).eq.0)
                   MEMLOG(size(help_sec_der_arr))
                  endif

                if(present(nuc_dervs_grarho)) then
                 allocate(help_3rd_der_arr(vl,10,unique_atoms(i_ua)%n_equal_atoms,occ_dim),&
                            stat=alloc_stat(19))
                   ASSERT(alloc_stat(19).eq.0)
                   MEMLOG(size(help_3rd_der_arr))
                   help_3rd_der_arr=0.0_r8_kind !!! to be deleted
                endif

                help_arr: if (vl == size(nuc, 1)) then

                   ! calculate help_nuc_arr to calc partional nuc_grarho

                   call dgemm ('n', 'n', vl*3*unique_atoms(i_ua)%n_equal_atoms, &
                        occ_dim, orb_dim, &
                        1.0_r8_kind, &
                        nuc(1, 1, 1, 1, l), &
                        vl*3*unique_atoms(i_ua)%n_equal_atoms, & ! vl == lower dim, see if
                        eigv(1 + counter, 1, s), eig_dim, 0.0_r8_kind, &
                        help_nuc_arr(1, 1, 1, 1),&
                        vl*3*unique_atoms(i_ua)%n_equal_atoms)

                  if(present(nuc_dervsrho).or.present(nuc_sec_derrho)) then
                   ! now help_sec_der_arr
                    call dgemm ('n', 'n', vl*6*unique_atoms(i_ua)%n_equal_atoms, &
                        occ_dim, orb_dim, &
                        1.0_r8_kind, &
                        sec_der(1, 1, 1, 1, l), &
                        vl*6*unique_atoms(i_ua)%n_equal_atoms,& ! vl == lower dim, see if
                        eigv(1 + counter, 1, s), eig_dim, 0.0_r8_kind, &
                        help_sec_der_arr(1, 1, 1, 1),&
                        vl*6*unique_atoms(i_ua)%n_equal_atoms)
                  endif

                  if(present(nuc_3rd_der)) then
 FPP_TIMER_START(t_nuc3rd)
                   call dgemm ('n', 'n', vl*10*unique_atoms(i_ua)%n_equal_atoms, &
                        occ_dim, orb_dim, &
                        1.0_r8_kind, &
                        nuc_3rd_der(i_ua,i_irr)%o(1, 1, 1, 1, l), &
                        vl*10*unique_atoms(i_ua)%n_equal_atoms, & ! vl ==lower dim, see if
                        eigv(1 + counter, 1, s), eig_dim, 0.0_r8_kind, &
                        help_3rd_der_arr(1, 1, 1, 1), &
                        vl*10*unique_atoms(i_ua)%n_equal_atoms)
 FPP_TIMER_STOP(t_nuc3rd)
                  endif

                else help_arr
                   ! copy data first

                    allocate(nuc_help(vl,3,unique_atoms(i_ua)%n_equal_atoms,orb_dim), &
                             stat=alloc_stat(12))
                    ASSERT(alloc_stat(12).eq.0)
                    MEMLOG(size(nuc_help))

                   nuc_help(1:vl,:,:,:)=nuc(1:vl,:,:,:,l)

                   ! calculate help_nuc_arr
                   call dgemm('n','n',vl*3*unique_atoms(i_ua)%n_equal_atoms,&
                        occ_dim,orb_dim,   1.0_r8_kind,&
                        nuc_help(1,1,1,1),vl*3*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim,0.0_r8_kind, &
                        help_nuc_arr(1,1,1,1),&
                        vl*3*unique_atoms(i_ua)%n_equal_atoms)

                 dervs: if(present(nuc_dervsrho).or.present(nuc_sec_derrho)) then
                   allocate(sec_der_help(vl,6,unique_atoms(i_ua)%n_equal_atoms,orb_dim),&
                            stat=alloc_stat(17))
                    ASSERT(alloc_stat(17).eq.0)
                    MEMLOG(size(sec_der_help))

                   ! now help_sec_der_arr
                   sec_der_help(1:vl,:,:,:)=sec_der(1:vl,:,:,:,l)

                   call dgemm('n','n',vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                        occ_dim,orb_dim,1.0_r8_kind,&
                        sec_der_help(1,1,1,1),vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim,0.0_r8_kind, &
                        help_sec_der_arr(1,1,1,1),&                !to calc nuc_dervsrho
                        vl*6*unique_atoms(i_ua)%n_equal_atoms)


                    MEMLOG(-size(sec_der_help))
                    deallocate(sec_der_help, stat=alloc_stat(17))
                    ASSERT(alloc_stat(17).eq.0)
                    alloc_stat(17)=1
                 endif dervs

                third_dervs: if(present(nuc_3rd_der)) then
 FPP_TIMER_START(t_nuc3rd)

                   allocate(sec_der_help(vl,10,unique_atoms(i_ua)%n_equal_atoms,orb_dim),&
                             stat=alloc_stat(20))
                    ASSERT(alloc_stat(20).eq.0)
                    MEMLOG(size(sec_der_help))

                   ! now help_sec_der_arr
                   sec_der_help(:vl,:,:,:)=nuc_3rd_der(i_ua,i_irr)%o(:vl,:,:,:,l)

                   call dgemm('n','n',vl*10*unique_atoms(i_ua)%n_equal_atoms,&
                        occ_dim,orb_dim,1.0_r8_kind,&
                        sec_der_help(1,1,1,1),vl*10*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim,0.0_r8_kind, &
                        help_3rd_der_arr(1,1,1,1),&                !to calc nuc_dervsrho
                        vl*10*unique_atoms(i_ua)%n_equal_atoms)


                    MEMLOG(-size(sec_der_help))
                    deallocate(sec_der_help, stat=alloc_stat(20))
                    ASSERT(alloc_stat(20).eq.0)
                    alloc_stat(20)=1


 FPP_TIMER_STOP(t_nuc3rd)
                endif third_dervs

                  MEMLOG(-size(nuc_help))
                  deallocate(nuc_help, stat=alloc_stat(12))
                     ASSERT(alloc_stat(12).eq.0)
                  alloc_stat(12)=1
                endif help_arr
!                if(present(nuc_3rd_der)) then
!                 print*,'sec_der_arr 3rd_der_arr',  &
!                        sum(help_sec_der_arr(:vl,1,1,:occ_dim)), &
!                        sum(help_3rd_der_arr(:vl,1,1,:occ_dim)),i_ua
!                endif



                   ! now final building of displacement gradient follows
                   nuc_grarho_p=>nuc_grarho(i_ua)%m

                   if(present(nuc_sec_derrho)) nuc_sec_derrho_p=>nuc_sec_derrho(i_ua)%m
                   ! compare with nuc_dervsrho(i_ua,i_ua)%m(1:vl,1,i_ea,1,i_ea,s) used in nonlocal dervs calcs
                   if(present(nuc_grad_tau)) nuc_grad_tau_p=>nuc_grad_tau(i_ua)%m

                   mo: do j=1,occ_dim

                      occ_real=occ(j,s)/partners(i)
                      occ_real_2=occ_real*2.0_r8_kind

                      equals: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms

                      !-------------------------------------------------------------------------------------------------+
                      ! calculate drho_dR                                                                               |
                      !-------------------------------------------------------------------------------------------------+
                      ! NOTE: phi_arr contains occ_real_2                                                               |
                      !-------------------------------------------------------------------------------------------------+
                      nuc_grarho(i_ua)%m(1:vl,:,i_ea,s)=&
                           nuc_grarho(i_ua)%m(1:vl,:,i_ea,s)+&
                           help_nuc_arr(1:vl,:,i_ea,j)*spread(phi_arr(:,j),2,3)



                         ! second derivatives follow
                         ! Note: Here we have to treat the full matrix, because
                         !       nuc_sec_derrho is not symmetric ( it is not

                         !       exactly the Hessian matrix, but sum^prime( Pmunu*Xmunu)
                         !       see Pople CPL 199 (1992) 557

                      derrho: if(present(nuc_sec_derrho)) then

                         !!! help_nuc_arr * graphi_arr !!! (1)
                         !----------------------------------------------------------------------------------------------+
                         ! calculate first term of "dgamma_dR":                                                         |
                         !----------------------------------------------------------------------------------------------+
                         nuc_sec_derrho_p(:vl,1,:,i_ea,s)=nuc_sec_derrho_p(:vl,1,:,i_ea,s)+&
                              occ_real_2*spread(help_nuc_arr(:vl,1,i_ea,j),2,3)*graphi_arr(:vl,:3,j)

                         nuc_sec_derrho_p(:vl,2,:,i_ea,s)=nuc_sec_derrho_p(:vl,2,:,i_ea,s)+&
                              occ_real_2*spread(help_nuc_arr(:vl,2,i_ea,j),2,3)*graphi_arr(:vl,:,j)

                         nuc_sec_derrho_p(:vl,3,:,i_ea,s)=nuc_sec_derrho_p(:vl,3,:,i_ea,s)+&
                              occ_real_2*spread(help_nuc_arr(:vl,3,i_ea,j),2,3)*graphi_arr(:vl,:,j)

                         !!! help_sec_der_arr * phi_arr(:,j)  !!! (2)
                         !----------------------------------------------------------------------------------------------+
                         ! calculate second term of "dgamma_dR":                                                        |
                         !----------------------------------------------------------------------------------------------+
                         !                                                                                              |
                         !                       /  1  2  3  \                                                          |
                         !    help_sec_der_arr = |  2  4  5  |                                                          |
                         !                       \  3  5  6  /                                                          |
                         !                                                                                              |
                         !                                                                                              |
                         !----------------------------------------------------------------------------------------------+
                         ! NOTE: phi_arr contains occ_real_2                                                            |
                         !----------------------------------------------------------------------------------------------+
                         nuc_sec_derrho_p(:vl,1,1,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,1,1,i_ea,s)+&
                              help_sec_der_arr(:vl,1,i_ea,j)*phi_arr(:,j)
                              !-------------- nuc_gra_graphi
                         nuc_sec_derrho_p(:vl,1,2,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,1,2,i_ea,s)+&
                              help_sec_der_arr(:vl,2,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,2,1,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,2,1,i_ea,s)+&
                              help_sec_der_arr(:vl,2,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,1,3,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,1,3,i_ea,s)+&
                              help_sec_der_arr(:vl,3,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,3,1,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,3,1,i_ea,s)+&
                              help_sec_der_arr(:vl,3,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,2,2,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,2,2,i_ea,s)+&
                              help_sec_der_arr(:vl,4,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,2,3,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,2,3,i_ea,s)+&
                              help_sec_der_arr(:vl,5,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,3,2,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,3,2,i_ea,s)+&
                              help_sec_der_arr(:vl,5,i_ea,j)*phi_arr(:,j)
                         nuc_sec_derrho_p(:vl,3,3,i_ea,s)=&
                              nuc_sec_derrho_p(:vl,3,3,i_ea,s)+&
                              help_sec_der_arr(:vl,6,i_ea,j)*phi_arr(:,j)

                         !----------------------------------------------------------------------------------------------+
                         ! calculate "dtau_dR":                                                                         |
                         !----------------------------------------------------------------------------------------------+
                         if (present(nuc_grad_tau)) then
                            nuc_grad_tau_p(1:vl,1,i_ea,s) = nuc_grad_tau_p(1:vl,1,i_ea,s) + &
                                       ( help_sec_der_arr(:vl,1,i_ea,j)*graphi_arr(:vl,1,j) &
                                       + help_sec_der_arr(:vl,2,i_ea,j)*graphi_arr(:vl,2,j) &
                                       + help_sec_der_arr(:vl,3,i_ea,j)*graphi_arr(:vl,3,j))&
                                       * occ_real
                            nuc_grad_tau_p(1:vl,2,i_ea,s) = nuc_grad_tau_p(1:vl,2,i_ea,s) + &
                                       ( help_sec_der_arr(:vl,2,i_ea,j)*graphi_arr(:vl,1,j) &
                                       + help_sec_der_arr(:vl,4,i_ea,j)*graphi_arr(:vl,2,j) &
                                       + help_sec_der_arr(:vl,5,i_ea,j)*graphi_arr(:vl,3,j))&
                                       * occ_real
                            nuc_grad_tau_p(1:vl,3,i_ea,s) = nuc_grad_tau_p(1:vl,3,i_ea,s) + &
                                       ( help_sec_der_arr(:vl,3,i_ea,j)*graphi_arr(:vl,1,j) &
                                       + help_sec_der_arr(:vl,5,i_ea,j)*graphi_arr(:vl,2,j) &
                                       + help_sec_der_arr(:vl,6,i_ea,j)*graphi_arr(:vl,3,j))&
                                       * occ_real
                         end if
                      endif derrho

FPP_TIMER_START(t_nuc_dervsrho)
                     dervsrho2: if(present(nuc_dervsrho)) then ! first of two contribs

FPP_TIMER_START(t_graphi_sec_der)
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,s)+ &
                              help_sec_der_arr(:vl,1,i_ea,j)*phi_arr(:,j)

                      help_dervsrho=help_sec_der_arr(:vl,2,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,s)+help_dervsrho
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,s)+help_dervsrho

                      help_dervsrho=help_sec_der_arr(:vl,3,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,s)+help_dervsrho
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,s)+help_dervsrho

                      help_dervsrho=help_sec_der_arr(:vl,4,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,s)+help_dervsrho

                      help_dervsrho=help_sec_der_arr(:vl,5,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,s)+help_dervsrho
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,s)+help_dervsrho

                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,s)= &
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,s)+ &
                              help_sec_der_arr(:vl,6,i_ea,j)*phi_arr(:,j)

FPP_TIMER_STOP(t_graphi_sec_der)

                     endif dervsrho2

                  dervs_grarho2: if(present(nuc_dervs_grarho)) then

                     !!! help_nuc_arr * graphi_arr ->
                     !!! help_sec_der_arr * graphi_arr  !!!       (3)

FPP_TIMER_START(t_graphi_sec_der)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,:3,s)+ &
                      occ_real_2*spread(help_sec_der_arr(:vl,1,i_ea,j),2,3)*graphi_arr(:vl,:3,j)

                     dervs_grarho_help(:vl,:3)=occ_real_2*graphi_arr(:vl,:3,j)* &
                              spread(help_sec_der_arr(:vl,2,i_ea,j),2,3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)

                     dervs_grarho_help(:vl,:3)=occ_real_2*graphi_arr(:vl,:3,j)* &
                        spread(help_sec_der_arr(:vl,3,i_ea,j),2,3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,:3,s)+ &
                      occ_real_2*spread(help_sec_der_arr(:vl,4,i_ea,j),2,3)*graphi_arr(:vl,:3,j)

                     dervs_grarho_help(:vl,:3)=occ_real_2*graphi_arr(:vl,:3,j)* &
                        spread(help_sec_der_arr(:vl,5,i_ea,j),2,3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,:3,s)+ &
                        dervs_grarho_help(:vl,:3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,:3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,:3,s)+ &
                      occ_real_2*spread(help_sec_der_arr(:vl,6,i_ea,j),2,3)*graphi_arr(:vl,:3,j)
FPP_TIMER_STOP(t_graphi_sec_der)

                    if(.true.) then  !!! help_3rd_der_arr * phi_arr !!! (2)
FPP_TIMER_START(t_nuc3rd)
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,1,s)= & !!! (1) xxx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,1,s)+ &
                              help_3rd_der_arr(:vl,1,i_ea,j)*phi_arr(:vl,j)

                     dervs_grarho_help(:vl,2)=help_3rd_der_arr(:vl,2,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,1,s)= & !!! (2) yxx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,2,s)= & !!! (3) xxy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,1,s)= & !!! (4) xyx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,2)

                     dervs_grarho_help(:vl,3)=help_3rd_der_arr(:vl,3,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,1,s)= & !!! (5) xzx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,1,s)= & !!! (6) zxx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,3,s)= & !!! (7) xxz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,3)

                     dervs_grarho_help(:vl,1)=help_3rd_der_arr(:vl,4,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,1,s)= & !!! (8) yyx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,1)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,2,s)= & !!! (9) yxy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,1)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,2,s)= & !!! (10) xyy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,1)

                     dervs_grarho_help(:vl,2)=help_3rd_der_arr(:vl,5,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,3,s)= & !!! (11) xyz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,2,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,3,s)= & !!! (12) yxz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,1,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,2,s)= & !!! (13) zxy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,1,s)= & !!! (14) zyx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,2,s)= & !!! (15) xzy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,1,s)= & !!! (16) yzx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,2)

                     dervs_grarho_help(:vl,1)=help_3rd_der_arr(:vl,6,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,1,s)= & !!! (17) zzx
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,1,s)+ &
                              dervs_grarho_help(:vl,1)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,3,s)= & !!! (18) xzz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,1,i_ea,3,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,1)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,3,s)= & !!! (19) zxz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,1,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,1)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,2,s)= & !!! (20) yyy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,2,s)+ &
                              help_3rd_der_arr(:vl,7,i_ea,j)*phi_arr(:,j)

                     dervs_grarho_help(:vl,3)=help_3rd_der_arr(:vl,8,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,2,s)= & !!! (21) yzy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,2,s)= & !!! (22) zyy
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,3)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,3,s)= & !!! (23) yyz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,2,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,3)

                     dervs_grarho_help(:vl,2)=help_3rd_der_arr(:vl,9,i_ea,j)*phi_arr(:,j)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,3,s)= & !!! (24) yzz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,2,i_ea,3,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,3,s)= & !!! (25) ZyZ
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,2,i_ea,3,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,2,s)= & !!! (26) ZyZ
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,2,s)+ &
                              dervs_grarho_help(:vl,2)

                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,3,s)= & !!! (27) zzz
                     nuc_dervs_grarho(i_ua,i_ua)%m(:vl,3,i_ea,3,i_ea,3,s)+ &
                              help_3rd_der_arr(:vl,10,i_ea,j)*phi_arr(:,j)
FPP_TIMER_STOP(t_nuc3rd)
                   endif
FPP_TIMER_STOP(t_nuc_dervsrho)

                  endif dervs_grarho2

                   enddo equals ! (i_ea)
                enddo mo

         dervsrho: if(.true..and.present(nuc_dervsrho)) then

FPP_TIMER_START(t_nuc_dervsrho)

               counter2=0
          ua2: do i_ua2=1,n_unique_atoms
               orb_dim2=size(nuc_grads(i_ua2,i)%o,4)
               orbdim2: if (orb_dim2 > 0) then

FPP_TIMER_START(t_help_nuc_arr2)
               calc_help_nuc_arr2: if(.not.allocated(help_nuc_arr2(i_ua2)%m)) then

               allocate(help_nuc_arr2(i_ua2)%m(vl,3,unique_atoms(i_ua2)%n_equal_atoms,occ_dim),&
                        stat=alloc_stat(14))
                       ! to hold orbital nuc displacement grads for second unique
               ASSERT(alloc_stat(14).eq.0)
               MEMLOG(size(help_nuc_arr2(i_ua2)%m))

               fullorelse: if (vl == size(nuc_grads(i_ua2, i)%o, 1)) then
                  call dgemm ('n', 'n', vl*3*unique_atoms(i_ua2)%n_equal_atoms,occ_dim,orb_dim2, &
                       1.0_r8_kind, &
                       nuc_grads(i_ua2,i)%o(1, 1, 1, 1, l), &
                       vl*3*unique_atoms(i_ua2)%n_equal_atoms, & ! vl == lower dim, see if
                       eigv(1 + counter2, 1, s), eig_dim, 0.0_r8_kind, &
                       help_nuc_arr2(i_ua2)%m(1, 1, 1, 1), &
                       vl*3*unique_atoms(i_ua2)%n_equal_atoms)

               else fullorelse
                  ! copy data first
                  allocate(nuc_help(vl,3,unique_atoms(i_ua2)%n_equal_atoms,&
                                                   orb_dim2), &
                           sec_der_help(vl,6,unique_atoms(i_ua2)%n_equal_atoms,&
                                                       orb_dim2),&
                           stat=alloc_stat(15))
                  ASSERT(alloc_stat(15).eq.0)
                  MEMLOG(size(nuc_help)+size(sec_der_help))

                   nuc_help(1:vl,:,:,:)=      nuc_grads(i_ua2,i)%o(1:vl,:,:,:,l)
                   sec_der_help(1:vl,:,:,:)=nuc_sec_der(i_ua2,i)%o(1:vl,:,:,:,l)

                   call dgemm('n','n',vl*3*unique_atoms(i_ua2)%n_equal_atoms,&
                                             occ_dim,orb_dim2, &
                        1.0_r8_kind,&
                        nuc_help(1,1,1,1),vl*3*unique_atoms(i_ua2)%n_equal_atoms,&
                        eigv(1+counter2,1,s),eig_dim, &
                        0.0_r8_kind,help_nuc_arr2(i_ua2)%m(1,1,1,1),&
                          vl*3*unique_atoms(i_ua2)%n_equal_atoms)

                   MEMLOG(-size(nuc_help)-size(sec_der_help))
                   deallocate(nuc_help,sec_der_help, stat=alloc_stat(15))
                   ASSERT(alloc_stat(15).eq.0)
                   alloc_stat(15)=1
                endif fullorelse
#if 1
               do j=1,occ_dim
                occ_real=occ(j,s)/partners(i)
                occ_real_2=occ_real*2.0_r8_kind
                help_nuc_arr2(i_ua2)%m(:vl,:,:,j)=help_nuc_arr2(i_ua2)%m(:vl,:,:,j)*occ_real_2
               enddo
#endif
               endif calc_help_nuc_arr2
FPP_TIMER_STOP(t_help_nuc_arr2)


               ! cross contrib to nuc_dervsrho, 2nd of two contribs
#if 1
             if(i_ua2.ge.i_ua) then
FPP_TIMER_START(t_cross_nuc_dervsrho)
              do j=1,occ_dim
#if 0
               occ_real=occ(j,s)/partners(i)
               occ_real_2=occ_real*2.0_r8_kind
               occ_real_2=1.0_r8_kind
#endif

              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
              do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,1,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,1,i_ea2,s)+ &
                         help_nuc_arr(:vl,1,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,s)+ &
                         help_nuc_arr(:vl,1,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,1,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,1,i_ea2,s)+ &
                         help_nuc_arr(:vl,2,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,3,i_ea2,s)= &
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,3,i_ea2,s)+ &
                          help_nuc_arr(:vl,1,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,1,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,1,i_ea2,s)+ &
                          help_nuc_arr(:vl,3,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,2,i_ea2,s)= &
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,2,i_ea2,s)+ &
                          help_nuc_arr(:vl,2,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,3,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,3,i_ea2,s)+ &
                         help_nuc_arr(:vl,2,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,2,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,2,i_ea2,s)+ &
                         help_nuc_arr(:vl,3,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,3,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,3,i_ea2,s)+ &
                         help_nuc_arr(:vl,3,i_ea,j)*help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
              enddo
              enddo
              enddo
FPP_TIMER_STOP(t_cross_nuc_dervsrho)
             endif

          dervs_grarho: if(present(nuc_dervs_grarho)) then
FPP_TIMER_START(t_cross_dervs_grarho)
            if(i_ua2.eq.i_ua) then
              do j=1,occ_dim
#if 0
               occ_real=occ(j,s)/partners(i)
               occ_real_2=occ_real*2.0_r8_kind
               occ_real_2=1.0_r8_kind
#endif

              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
              do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms

                   !!! help_sec_der_arr * help_nuc_arr !!! (1)

                   do kk =1,3
                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,1,i_ea,j)

                     !first index - help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)+ &
                       dervs_grarho_help(:vl,kk)
                     !first index - help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)+ &
                       dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,2,i_ea,j) !<-graphi_arr(:vl,:,j)
                     !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,3,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)= help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,4,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,5,i_ea,j) !<-graphi_arr(:vl,:,j)
                      !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,6,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                    enddo
              enddo
              enddo
              enddo
            elseif(i_ua.gt.i_ua2) then
              do j=1,occ_dim
#if 0
               occ_real=occ(j,s)/partners(i)
               occ_real_2=occ_real*2.0_r8_kind
               occ_real_2=1.0_r8_kind
#endif

              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
              do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms
                   do kk =1,3
                     !!!first index - help_nuc_arr
                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,1,i_ea,j)

                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)+ &
                       dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,2,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,3,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)= help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,4,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,5,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,6,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                    enddo
              enddo
              enddo
              enddo
            else
              do j=1,occ_dim
#if 0
               occ_real=occ(j,s)/partners(i)
               occ_real_2=occ_real*2.0_r8_kind
               occ_real_2=1.0_r8_kind
#endif

              do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
              do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms

                   do kk =1,3
                     !first index - help_sec_der_arr
                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,1,i_ea,j)

                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)+ &
                       dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,2,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,3,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)= help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,4,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,5,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,6,i_ea,j) !<-graphi_arr(:vl,:,j)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                    enddo
              enddo
              enddo
              enddo
            endif
FPP_TIMER_STOP(t_cross_dervs_grarho)
       endif dervs_grarho
#endif

#if 0
                mo2: do j=1,occ_dim

                occ_real=occ(j,s)/partners(i)
                occ_real_2=occ_real*2.0_r8_kind

                   equals1: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   equals2: do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms



                  ! in expresions below help_nuc_arr and help_nuc_arr2 do not
                  ! depend on occupation, thus occ_real_2 factor is used
                  ! to account for occupation dependence

FPP_TIMER_START(t_cross_nuc_dervsrho)
                  if(i_ua2.ge.i_ua) then
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,1,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,1,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,1,i_ea,j)* &
                          help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,1,i_ea,j)* &
                          help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,1,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,1,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,2,i_ea,j)* &
                          help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,3,i_ea2,s)= &
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,3,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(:vl,1,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,1,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,1,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(:vl,3,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,2,i_ea2,s)= &
                      nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,2,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(:vl,2,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,3,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,2,i_ea,3,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,2,i_ea,j)* &
                          help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,2,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,2,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,3,i_ea,j)* &
                          help_nuc_arr2(i_ua2)%m(:vl,2,i_ea2,j)

                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,3,i_ea2,s)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,3,i_ea,3,i_ea2,s)+ &
                         occ_real_2*help_nuc_arr(:vl,3,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(:vl,3,i_ea2,j)
                    endif
FPP_TIMER_STOP(t_cross_nuc_dervsrho)
                 dervs_grarho: if(present(nuc_dervs_grarho)) then
FPP_TIMER_START(t_cross_dervs_grarho)


                     !!! help_sec_der_arr * help_nuc_arr !!! (1)

                   do kk =1,3
                     dervs_grarho_help(:vl,kk)=occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,1,i_ea,j)

                     !first index - help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,1,s)+ &
                       dervs_grarho_help(:vl,kk)
                     !first index - help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,1,s)+ &
                       dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,2,i_ea,j) !<-graphi_arr(:vl,:,j)
                     !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,3,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,1,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,1,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,1,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,1,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)= occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,4,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,5,i_ea,j) !<-graphi_arr(:vl,:,j)
                      !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,2,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,2,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,2,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,2,s)+ &
                      dervs_grarho_help(:vl,kk)

                     dervs_grarho_help(:vl,kk)=occ_real_2*help_nuc_arr2(i_ua2)%m(:vl,kk,i_ea2,j)* &
                              help_sec_der_arr(:vl,6,i_ea,j) !<-graphi_arr(:vl,:,j)
                       !!! first index help_nuc_arr
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)= &
                     nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,3,i_ea,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                       !!! first index help_sec_der_arr
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,3,i_ea,kk,i_ea2,3,s)+ &
                      dervs_grarho_help(:vl,kk)
                    enddo

FPP_TIMER_STOP(t_cross_dervs_grarho)

                 endif dervs_grarho

                  enddo equals2
                  enddo equals1
               enddo mo2
#endif

               counter2=counter2+orb_dim2
               endif orbdim2

              enddo ua2

 FPP_TIMER_STOP(t_nuc_dervsrho)
         endif dervsrho

                counter=counter+orb_dim

                 MEMLOG(-size(help_nuc_arr))
                 deallocate(help_nuc_arr, stat=alloc_stat(10))
                 ASSERT(alloc_stat(10).eq.0)
                 alloc_stat(10)=1

                 if(present(nuc_dervsrho).or.present(nuc_sec_derrho)) then
                  MEMLOG(-size(help_sec_der_arr))
                  deallocate(help_sec_der_arr, stat=alloc_stat(18))
                  ASSERT(alloc_stat(18).eq.0)
                  alloc_stat(18)=1
                 endif

                 if(present(nuc_dervs_grarho)) then
                  MEMLOG(-size(help_3rd_der_arr))
                  deallocate(help_3rd_der_arr, stat=alloc_stat(19))
                  ASSERT(alloc_stat(19).eq.0)
                  alloc_stat(19)=1
                 endif

                endif orb_dim1
    enddo uniques ! (i_ua)

             if(present(nuc_dervsrho)) then
              do i_ua2=1,n_unique_atoms
               if(allocated(help_nuc_arr2(i_ua2)%m)) then
                MEMLOG(-size(help_nuc_arr2(i_ua2)%m))
                deallocate(help_nuc_arr2(i_ua2)%m, stat=alloc_stat(14))
                ASSERT(alloc_stat(14).eq.0)
                alloc_stat(14)=1
               endif
               enddo
             endif

          enddo partns
       enddo spins


     if(present(graphi)) then
       MEMLOG(-size(phi_arr))
       deallocate(phi_arr,stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(8)=1
     else
       MEMLOG(-size(phi_arr)-size(graphi_arr))
       deallocate(phi_arr,graphi_arr,stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(8)=1
       alloc_stat(34)=1 ! graphi_arr
     endif

       endif
       deallocate( eigv )
       end associate
    enddo irr

         if(present(nuc_dervsrho)) then

          MEMLOG(-size(help_nuc_arr2))
          deallocate(help_nuc_arr2,stat=alloc_stat(16))
          ASSERT(alloc_stat(16).eq.0)
          alloc_stat(16)=1
#if 1
             do i_ua=1,n_unique_atoms
              do i_ua2=i_ua+1,n_unique_atoms
                   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms
                   do k=1,3
                   do kk=1,3
!                   print*,'nuc_dervsrho', i_ua,i_ea, i_ua2,i_ea2,&
!                   sum(nuc_dervsrho(i_ua2,i_ua)%m(:vl,2,i_ea2,1,i_ea,1)), &
!                   sum(nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,1)), &
!                   sum(abs(nuc_dervsrho(i_ua2,i_ua)%m(:vl,2,i_ea2,1,i_ea,1) &
!                          -nuc_dervsrho(i_ua,i_ua2)%m(:vl,1,i_ea,2,i_ea2,1)))

                    nuc_dervsrho(i_ua2,i_ua)%m(:vl,kk,i_ea2,k,i_ea,:)= &
                     nuc_dervsrho(i_ua,i_ua2)%m(:vl,k,i_ea,kk,i_ea2,:)
                   enddo
                   enddo
                   enddo
                   enddo
              enddo
             enddo
#endif
         endif

#if 1 /* if  nuc_dervs_grarho symmetric */
            if(present(nuc_dervs_grarho)) then
             do i_ua=1,n_unique_atoms
              do i_ua2=i_ua+1,n_unique_atoms
                   do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
                   do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms
                   do k=1,3
                   do kk=1,3
                    nuc_dervs_grarho(i_ua2,i_ua)%m(:vl,kk,i_ea2,k,i_ea,:,:)= &
                     nuc_dervs_grarho(i_ua,i_ua2)%m(:vl,k,i_ea,kk,i_ea2,:,:)
                   enddo
                   enddo
                   enddo
                   enddo
              enddo
             enddo
             endif
#endif


    if(ispin==1) then
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)
    else
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)

       gamma(1:vl,DNDN) = SUM(grarho(1:vl,X:Z,DN)**2, DIM=2)

       gamma(1:vl,UPDN) = SUM(grarho(1:vl,X:Z,UP)*grarho(1:vl,X:Z,DN), DIM=2)
    endif

  end subroutine density_calc_ph_nl

  subroutine print_alloc_density_calc()
    integer(kind=i4_kind):: i_all
    DPRINT 'density module alloc control list'
    do i_all=1,size(alloc_stat)
    if(alloc_stat(i_all).eq.0) print*,i_all,' alloc_stat eq 0'
    enddo
    DPRINT 'density_calc deallocations checked'
  end subroutine print_alloc_density_calc

   function fitted_density_is_closed()
       logical:: fitted_density_is_closed
       fitted_density_is_closed=fitted_density_is_closed_priv
   end function fitted_density_is_closed

end module density_calc_module

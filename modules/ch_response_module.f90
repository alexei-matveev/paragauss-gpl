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
module ch_response_module
  !---------------------------------------------------------------
  !
  !  Purpose: Correct allocated, calculate and deallocated charge
  !    fit basis function with all irreps and partners.
  !    Previous it was only for full symmetrical case, but
  !    for the implementation of symmetry into TDDFT it is
  !    NOT enough. One should have charge fit basis of all
  !    symmetry.
  !
  !
  !  Module called by: response_module.f90
  !
  !
  !  References:
  !   based on
  !           orbital_module.f90   -> fit_fct_allocate, fit_fct_calculate
  !           fit_coeff_module.f90 -> dimensions_calc
  !
  !
  !  Author: SB
  !  Date: 12/2004
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
  use symmetry_data_module ! description of irreps
  use orbitalstore_module
  use machineparameters_module
  use unique_atom_module
  use solid_harmonics_module ! calculates solid harmonics

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  type solid_harmonics_array_type
     ! solid harmonics of all equal atoms of one unique atome
     type(solid_harmonics_type), pointer :: sh(:)
     type(solid_harmonics_grads_type), pointer :: shg(:)
     type(solid_harmonics_sec_der_type), pointer :: shs(:)
  end type solid_harmonics_array_type


!!$  type glob_con_intermediate_l_type
!!$     ! to store intermediate values of global contraction of one unique atom:
!!$     real(kind=r8_kind), pointer :: oo  (:,:,:) ! kept for compatibility
!!$     ! oo  (grid_length,N_exponents,N_independent_fcts)
!!$     real(kind=r8_kind), pointer :: o   (:,:,:)
!!$     ! o   (grid_length,N_exponents,N_independent_fcts)
!!$     real(kind=r8_kind), pointer :: og  (:,:,:,:)
!!$     ! og  (grid_length,3,N_exponents,N_independent_fcts)
!!$     real(kind=r8_kind), pointer :: osd (:,:,:,:)
!!$     ! osd (grid_length,6,N_exponents,N_independent_fcts)
!!$     ! 1: xx  2: xy  3: xz  4: yy  5: yz  6: zz
!!$  end type glob_con_intermediate_l_type

  type(solid_harmonics_array_type), allocatable, target, private :: &
       solid_harmonics(:) ! solid_harmonics(N_unique_atoms)

!!$  type glob_con_intermediate_type
!!$     ! to store intermediate values of global contraction of one unique atom:
!!$     type(glob_con_intermediate_l_type), pointer :: l(:) ! l(-1:lmax) s,r2,p,...
!!$  end type glob_con_intermediate_type

!!$  type(glob_con_intermediate_type), pointer :: uncont_orbs_ch(:)
!!$  ! uncontracted primitive basis functions needed as intermediates
!!$  ! for subroutine calculate_orbital if global contractions are to be done

  real(kind=r8_kind), allocatable, private :: &
       prim_orb_fac(:,:,:),  cont_orb_fac(:,:,:), &
       prim_grad_fac(:,:,:), cont_grad_fac(:,:,:), &
       prim_sd_fac(:,:,:),   cont_sd_fac(:,:,:)

!!$  integer, pointer, public :: orbitalprojection_ch(:,:,:)
  ! orbital_index_??(-1:N_max_l,N_unique_atoms)
  ! first index of given l (-1 means r**2) and ua
  ! only totalsymmetric Irrep allowed
!!$  integer, pointer, public :: orbitalprojection_globcontr_ch(:,:)


  integer(kind=i4_kind)           :: N_max_exponents, N_max_equal_atoms,as, N_ea
!!$  logical                         :: do_glob_con_ch


  integer(i4_kind), parameter ::&
       & XX = 1, &
       & XY = 2, &
       & XZ = 3, &
       & YY = 4, &
       & YZ = 5, &
       & ZZ = 6

  real(r8_kind), parameter ::&
       & ZERO = 0.0_r8_kind,&
       & ONE  = 1.0_r8_kind,&
       & TWO  = 2.0_r8_kind

  logical :: parameters_calculated = .false.

  !------------ public functions and subroutines ------------------
  public fit_fct_allocate_response,fit_fct_free_response,fit_fct_calculate_response, &
       orbital_setup_response, dimension_of_fit_ch, orb_position, fit_position

  !================================================================
  ! End of public interface of module
  !================================================================


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine fit_fct_allocate_response(fcts_ch, grads, sec_ders)
    !  Purpose: allocates the arrays to hold fit functions
    !  and their various derivates on a set of grid points,
    !  only for ch-type, in other fords, charge fit basis
    !
    !  based on: fit_fct_allocate from orbital_module.f90
    !

    implicit none
    !------------ Declaration of formal parameters ---------------
    type(orbital_type)                 , pointer , optional :: fcts_ch(:)
    type(orbital_gradient_type)        , pointer , optional :: grads(:)
    type(orbital_sec_der_type)         , pointer , optional :: sec_ders(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer :: n_dim, i_ir, n_partn
    integer :: stat_alloc, vec_len, n_irreps
    !------------ Declaration of local handles -------------------
!!$type(unique_atom_type), pointer :: ua
    !------------ Executable code --------------------------------

    n_irreps = symmetry_data_n_irreps() !! Number of irreps

    if (present(fcts_ch)) then
       allocate(fcts_ch(n_irreps), stat=stat_alloc)
       if (stat_alloc /= 0) call error_handler( &
            "fit_fct_allocate_repsonse: allocate of fcts_ch failed")
    end if

    if (present(grads)) then
       allocate(grads(n_irreps), stat=stat_alloc)
       if (stat_alloc /= 0) call error_handler( &
            "fit_fct_allocate_response: allocate of grads failed")
    end if

    if (present(sec_ders)) then
       allocate(sec_ders(n_irreps), stat=stat_alloc)
       if (stat_alloc /= 0) call error_handler( &
            "fit_fct_allocate_response: allocate of sec_ders failed")
    end if


    vec_len = machineparameters_veclen
    do i_ir=1,n_irreps
       n_dim   = dimension_of_fit_ch(i_ir)
       n_partn = symmetry_data_n_partners(i_ir)

       if (present(fcts_ch)) then
          allocate( fcts_ch(i_ir)%o(vec_len,n_dim,n_partn), stat=stat_alloc)
          if (stat_alloc /= 0) call error_handler( &
               "fit_fct_allocate: allocate of fcts%o failed")
       endif

       if (present(grads)) then
          allocate( grads(i_ir)%o(vec_len,3,n_dim,n_partn), stat=stat_alloc)
          if (stat_alloc /= 0) call error_handler( &
               "fit_fct_allocate: allocate of grads%o failed")
       endif

       if (present(sec_ders)) then
          allocate( sec_ders(i_ir)%o(vec_len,6,n_dim,n_partn), stat=stat_alloc)
          if (stat_alloc /= 0) call error_handler( &
               "fit_fct_allocate: allocate of sec_ders%o failed")
       endif
    end do

  end subroutine fit_fct_allocate_response
  !*************************************************************

  !*************************************************************
  subroutine fit_fct_free_response(fcts, grads, sec_ders)
    !  Purpose: deallocates the arrays ment to hold fit functions
    !  and their various derivates on a set of grid points
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(orbital_type)                          , optional :: fcts(:)
    type(orbital_gradient_type)                 , optional :: grads(:)
    type(orbital_sec_der_type)                  , optional :: sec_ders(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer :: as, i_ir, n_irreps
    !------------ Executable code --------------------------------

    n_irreps = symmetry_data_n_irreps()

    if (present(fcts)) then
       do i_ir = 1,n_irreps
          deallocate( fcts(i_ir)%o, stat=as); if (as /= 0) call error_handler("fffr: deallocate of fcts%o failed")
       end do
    endif
    if (present(grads)) then
       do i_ir = 1,n_irreps
          deallocate( grads(i_ir)%o, stat=as); if (as /= 0) call error_handler("fffr: deallocate of grads%o failed")
       end do
    endif
    if (present(sec_ders)) then
       do i_ir = 1,n_irreps
          deallocate( sec_ders(i_ir)%o, stat=as); if (as /= 0) call error_handler( "fffr: deallocate of sec_ders%o failed")
       end do
    endif
  end subroutine fit_fct_free_response
  !*************************************************************


  !*************************************************************
  subroutine orbital_setup_response(vl,do_gradients,do_sec_der)
    !  Purpose: initialising and allocating module-privat variables
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: vl
    logical, optional, intent(in) :: do_gradients, do_sec_der
    ! default: false
    !** End of interface ***************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)       :: i_ua,i_ea,N_ea,i_l
    type(unique_atom_type), pointer   :: ua
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------

    !!    vl = machineparameters_veclen

    if ( .not. parameters_calculated ) then
       ! determine N_max_exponents, N_equal_atoms, do_glob_con_xx and N_max_l_xx
       N_max_exponents = 0
       N_max_equal_atoms = 0
!!$       do_glob_con_ch = .false.
       do i_ua = 1, N_unique_atoms
          ua => unique_atoms(i_ua)
          do i_l = 0, ua%lmax_ob
             N_max_exponents = max( N_max_exponents, ua%l_ob(i_l)%N_exponents )
          enddo
          do i_l = 0, ua%lmax_ch
             N_max_exponents = max( N_max_exponents, ua%l_ch(i_l)%N_exponents )
          enddo
          N_max_exponents = max( N_max_exponents, ua%r2_ch%N_exponents )
          N_max_equal_atoms = max( N_max_equal_atoms, ua%N_equal_atoms )
!!$          do_glob_con_ch = do_glob_con_ch .or. ua%N_glob_cons_ch .gt. 0
       enddo
       parameters_calculated = .true.
    endif

    ! initialise solid_harmonics_module
    call solid_harmonics_setup( vl )

    ! allocate solid_harmonics
    allocate( solid_harmonics(N_unique_atoms), stat=as )
    if ( as .ne. 0 ) call error_handler( &
         "orbital_setup: allocate of solid harmonics failed")
    if ( do_sec_der ) then
       do i_ua = 1, N_unique_atoms
          N_ea = unique_atoms(i_ua)%N_equal_atoms
          allocate( solid_harmonics(i_ua)%sh(N_ea), &
               solid_harmonics(i_ua)%shg(N_ea), &
               solid_harmonics(i_ua)%shs(N_ea), &
               stat = as )
          if ( as .ne. 0 ) call error_handler( &
               "orbital_setup: allocate of solid harmonics%sh failed")
          do i_ea = 1,N_ea
             call solid_harmonics_allocate( vl, &
                  unique_atoms(i_ua)%lmax_all, &
                  solid_harmonics(i_ua)%sh(i_ea), &
                  solid_harmonics(i_ua)%shg(i_ea), &
                  solid_harmonics(i_ua)%shs(i_ea) )
          enddo
       enddo
    elseif ( do_gradients ) then
       do i_ua = 1, N_unique_atoms
          N_ea = unique_atoms(i_ua)%N_equal_atoms
          allocate( solid_harmonics(i_ua)%sh(N_ea), &
               solid_harmonics(i_ua)%shg(N_ea), &
               stat=as )
          if ( as .ne. 0 ) call error_handler( &
               "orbital_setup: allocate of solid harmonics%sh failed")
          do i_ea = 1,N_ea
             call solid_harmonics_allocate( vl, &
                  unique_atoms(i_ua)%lmax_all, &
                  solid_harmonics(i_ua)%sh(i_ea), &
                  solid_harmonics(i_ua)%shg(i_ea) )
          enddo
       enddo
    endif

    ! allocate prim_orb_fac(vec_len,N_max_exponents,N_max_equal_atoms)
    allocate( prim_orb_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
    if ( as .ne. 0 ) call error_handler( &
         "orbital__setup: allocate of prim_orb_fac failed")

    ! allocate cont_orb_fac(v_len,N_max_exponents,N_max_equal_atoms)
    allocate( cont_orb_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
    if ( as .ne. 0 ) call error_handler( &
         "orbital_setup: allocate of cont_orb_fac failed")


    if ( do_gradients ) then

       ! allocate prim_grad_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( prim_grad_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
       if ( as .ne. 0 ) call error_handler( &
            "orbital__setup: allocate of prim_grad_fac failed")

       ! allocate cont_grad_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( cont_grad_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
       if ( as .ne. 0 ) call error_handler( &
            "orbital_setup: allocate of cont_grad_fac failed")


    else
       ! dummy allocation to facilitate call of combine_rad_and_y(00/lm)
       allocate( prim_grad_fac(1,1,1), stat=as )
       allocate( cont_grad_fac(1,1,1), stat=as )
       if (as /= 0) call error_handler( &
            "orbital_setup: dummy allocate of (prim/cont)_grad_fac failed")
    endif

    if ( do_sec_der ) then

       ! allocate prim_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( prim_sd_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
       if ( as .ne. 0 ) call error_handler( &
            "orbital__setup: allocate of prim_grad_fac failed")

       ! allocate cont_sd_fac(vec_len,N_max_exponents,N_max_equal_atoms)
       allocate( cont_sd_fac(vl,N_max_exponents,N_max_equal_atoms), stat=as )
       if ( as .ne. 0 ) call error_handler( &
            "orbital_setup: allocate of cont_grad_fac failed")

    else
       ! dummy allocation to facilitate call of combine_rad_and_y(00/lm)
       allocate( prim_sd_fac(1,1,1), stat=as )
       allocate( cont_sd_fac(1,1,1), stat=as )
       if (as /= 0) call error_handler( &
            "orbital_setup: dummy allocate of (prim/cont)_sd_fac failed")
    endif


  end subroutine orbital_setup_response
  !*************************************************************


  !*************************************************************
  subroutine fit_fct_calculate_response(grid_points,N_points, &
       fcts, grads, sec_ders)
    !  Purpose: calculates the values of the fit functions
    !  and their various derivates on a set of grid points
    !  for all symmetry
    !
    !  based on: fit_fct_calculate from orbital_module.f90
    !
    USE DEBUG
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in) :: grid_points(:,:)
    ! give this argument only at first call for given grid_points.
    ! Solid harmonics are evaluated when this argument is present and
    ! remain stored
    integer, intent(in) :: N_points
    ! default: vector_length from orbital_setup
    ! either "ch" for charge fit functions or "xc" for exchange fit functions
    type(orbital_type)                 , target           :: fcts(:)     ! (n_irr)
    type(orbital_gradient_type)        , target, optional :: grads(:)    ! (n_irr)
    type(orbital_sec_der_type)         , target, optional :: sec_ders(:) ! (n_irr)
    ! the various types of fit function arrays to load (all intent(out))
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind), parameter :: zero = 0.0_r8_kind
    integer(kind=i4_kind) :: vl, i_ua, i_ea, i_ir, &
         i_bas, i_l, i_exp, i_c, i_if, i_fu, m, &
         i_pa, i_sa, i_sa0, i_step
    real(kind=r8_kind)    :: coef
    logical               :: do_sec_ders, do_grads, do_fcts, &
         is_r2function, orbital_read
    ! pointers to the output variables
    real(kind=r8_kind), pointer :: so   (:,:,:)     ! => fcts%o
    real(kind=r8_kind), pointer :: sog  (:,:,:,:)   ! => grads%o
    real(kind=r8_kind), pointer :: sosd (:,:,:,:)   ! => sec_ders%o
    ! basis set information
    type(unique_atom_type)              , pointer :: ua
    integer(kind=i4_kind)                         :: lmax
    type(unique_atom_basis_type)        , pointer :: r2_bas
    type(unique_atom_basis_type)        , pointer :: l_bas(:)
    type(unique_atom_basis_type)        , pointer :: uab
    ! symmetrization information
    type(unique_atom_partner_type)      , pointer :: uap
    type(unique_atom_symadapt_type)     , pointer :: uas
    ! solid harmonics information
    type(solid_harmonics_array_type)    , pointer :: sha
    type(solid_harmonics_type)          , pointer :: sh
    type(solid_harmonics_grads_type)    , pointer :: shg
    type(solid_harmonics_sec_der_type)  , pointer :: shsd
    real(kind=r8_kind)                  , pointer :: r2(:)
    real(kind=r8_kind)                  , pointer :: x(:), y(:), z(:)
    real(kind=r8_kind)                  , pointer :: shm(:,:)
    real(kind=r8_kind)                  , pointer :: shm_grad(:,:,:)
    real(kind=r8_kind)                  , pointer :: shm_sd(:,:,:)
    ! working arrays
    real(kind=r8_kind), allocatable :: temp(:)
    real(kind=r8_kind), allocatable :: help(:), help_g(:), help_gy(:), &
         help_x(:), help_y(:), help_z(:), &
         help_xx(:), help_xy(:), help_xz(:), &
         help_yy(:), help_yz(:), help_zz(:)
    integer(i4_kind) :: indx(size(fcts)) ! (n_irr)
    real(kind=r8_kind) :: NORMAB, n_eq
    !------------ Declaration of external functions --------------
    external error_handler
    !------------ Executable code --------------------------------

    do_sec_ders  = present(sec_ders)
    do_grads     = do_sec_ders .or. present(grads)
    do_fcts      = .true.

    if (N_points > machineparameters_veclen) call error_handler( &
            "fit_fct_calculate: N_points larger than vec_len")
    vl = N_points

    ! set pointers on the output arrays (if present) else nullify
    ! first for the fit functions and their electronic derivatives

    ! Return point for calculated orbitals
    if(orbitalstore_initialized .and. orbital_read) then
       return
    end if

    ! allocate working arrays
    allocate( temp(N_max_exponents), stat=as )
    if (as /= 0) call error_handler( &
         "fit_fct_calculate: allocate temp failed")
    if (do_fcts) then
       allocate( help(vl), stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: allocate help failed")
    end if
    if (do_grads) then
       allocate( help_x(vl), help_y(vl), help_z(vl), &
            help_g(vl), help_gy(vl), stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{x,y,z} and help_g{,y} failed")
    end if
    if (do_sec_ders) then
       allocate( help_xx(vl), help_xy(vl), help_xz(vl), &
            help_yy(vl), help_yz(vl), help_zz(vl), stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: allocate help_{xx,xy,xz,yy,yz,zz} failed")
    end if

    indx = 0
    ! loop over unique atoms
    i_ua_: do i_ua=1,N_unique_atoms
       ua  => unique_atoms(i_ua)
       sha => solid_harmonics(i_ua)

       ! load solid harmonics and their derivatives on the grid
       do i_ea=1,ua%N_equal_atoms
          ! load solid_harmonics even if do_moving_fcts = .false.
          call solid_harmonics_calculate( &
               sha%sh(i_ea), grid_points, ua%position(:,i_ea),N_points )
          call solid_harmonics_calculate_grads( &
               sha%sh(i_ea), sha%shg(i_ea), N_points )
          call solid_harmonics_calc_sec_der( &
               sha%shg(i_ea), sha%shs(i_ea), N_points )
       end do

       ! set pointer on the basis set information
       lmax   =  ua%lmax_ch
       r2_bas => ua%r2_ch
       l_bas  => ua%l_ch

       ! start evaluating the fit functions on the grid


       ! loop over all types of basis functions (s,r2,p,d,...)
       i_bas_: do i_bas = -1,lmax
          ! check i_bas to set i_l and is_r2function
          NORMAB = 1.0_r8_kind
          n_eq = unique_atoms(i_ua)%n_equal_atoms
          select case (i_bas)
          case (-1) ! s-type
             i_l = 0
             is_r2function = .false.
             uab => l_bas(0)
             NORMAB = sqrt(n_eq)
          case (0) ! r2-type
             i_l = 0
             is_r2function = .true.
             uab => r2_bas
             NORMAB = sqrt(n_eq)
          case default
             i_l = i_bas
             is_r2function = .false.
             uab => l_bas(i_l)
          end select

          ! first evaluate the radial part

          temp(:uab%N_exponents) = - ( uab%exponents + uab%exponents )

          do i_ea = 1,ua%N_equal_atoms
             sh => sha%sh(i_ea)
             r2 => sh%r2

             ! calculate the radial part of the primitive orbitals
             if (is_r2function) then
                !              R(r) = r^2 exp(-ar^2)
                ! [1/r d/dr]   R(r) = ( 2 - 2ar^2 ) exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = ( 4a^2r^2 - 8a ) exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   help = exp( - uab%exponents(i_exp) * r2(:vl) ) * NORMAB !! * NORMAB added
                   prim_orb_fac(:vl,i_exp,i_ea) = r2(:vl) * help
                   help = help + help
                   prim_grad_fac(:vl,i_exp,i_ea) = help + &
                        temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   prim_sd_fac(:vl,i_exp,i_ea) = temp(i_exp) * ( help + &
                        prim_grad_fac(:vl,i_exp,i_ea) )
                end do
             else
                !              R(r) = N exp(-ar^2)
                ! [1/r d/dr]   R(r) = -2a N exp(-ar^2)
                ! [1/r d/dr]^2 R(r) = 4a^2 N exp(-ar^2)
                do i_exp = 1,uab%N_exponents
                   prim_orb_fac(:vl,i_exp,i_ea) = &
                        uab%norms(i_exp) * exp( - uab%exponents(i_exp) * r2(:vl) ) * NORMAB !! * NORMAB added
                   prim_grad_fac(:vl,i_exp,i_ea) = &
                        temp(i_exp) * prim_orb_fac(:vl,i_exp,i_ea)
                   prim_sd_fac(:vl,i_exp,i_ea) = &
                        temp(i_exp) * prim_grad_fac(:vl,i_exp,i_ea)
                end do
             end if

             ! calculate the radial part of the local contractions
             ! [1/r d/dr]^n C_i(r) = Sum(a) coeff(a,i) [1/r d/dr]^n R_a(r)
             do i_c = 1,uab%N_contracted_fcts ! contractions
                cont_orb_fac(:vl,i_c,i_ea) = zero
                do i_exp = 1,uab%N_exponents ! contributing exponents
                   coef = uab%contractions(i_exp,i_c)
                   cont_orb_fac(:vl,i_c,i_ea) = cont_orb_fac(:vl,i_c,i_ea)+&
                        coef * prim_orb_fac(:vl,i_exp,i_ea)
                end do ! contributing exponents
                cont_grad_fac(:vl,i_c,i_ea) = zero
                do i_exp = 1,uab%N_exponents ! contributing exponents
                   coef = uab%contractions(i_exp,i_c)
                   cont_grad_fac(:vl,i_c,i_ea)=cont_grad_fac(:vl,i_c,i_ea)+&
                        coef * prim_grad_fac(:vl,i_exp,i_ea)
                end do! contributing exponents
                cont_sd_fac(:vl,i_c,i_ea) = zero
                do i_exp = 1,uab%N_exponents ! contributing exponents
                   coef = uab%contractions(i_exp,i_c)
                   cont_sd_fac(:vl,i_c,i_ea) = cont_sd_fac(:vl,i_c,i_ea) + &
                        coef * prim_sd_fac(:vl,i_exp,i_ea)
                end do! contributing exponents
             end do! contractions

          end do! equal atoms

          ! symmetry adaption:
          ! for all independent fcts,
          !         uncontracted exponents and contractions,
          !         gridpoints:
          !   sum over equal atoms and m
          !
          ! do renormalisation of calculated orbitals in the same sweep


          i_ir_: do i_ir = 1, symmetry_data_n_irreps()
             so => fcts(i_ir)%o
             if (present(grads)) then
                sog => grads(i_ir)%o
             else
                nullify(sog)
             end if
             if (present(sec_ders)) then
                sosd => sec_ders(i_ir)%o
             else
                nullify(sosd)
             end if

             i_partner_: do i_pa = 1, symmetry_data_n_partners(i_ir)
                i_sa0 = indx(i_ir)

                uap    => ua%symadapt_partner(i_ir,i_l)
                i_step =  uap%N_independent_fcts
                i_if_: do i_if = 1,uap%N_independent_fcts ! independent fcts
                   uas => uap%symadapt(i_if,i_pa) !! i_pa = 1, for C1

                   i_fu_: do i_fu = 1,uas%N_fcts ! contributing functions
                      i_sa = i_sa0 + i_if

                      i_ea = uas%I_equal_atom(i_fu)
                      coef = uas%c(i_fu)
!!! x,y,z needed for grad and sec_ders
                      x => sha%sh(i_ea)%l(1)%m(:,2)
                      y => sha%sh(i_ea)%l(1)%m(:,3)
                      z => sha%sh(i_ea)%l(1)%m(:,1)

                      if (i_l .eq. 0) then
                         ! accumulate primitives in fcts, grads, ...
                         call combine_rad_and_y00(i_sa, uab%N_uncontracted_fcts, coef, i_step, i_pa, i_ea, & !! 1.0_r8_kind = coeff
                              prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                              so, sog, sosd)
                         ! accumulate contractions in fcts, grads, ...
                         call combine_rad_and_y00(i_sa, uab%N_contracted_fcts,   coef, i_step, i_pa, i_ea, &   !! 1.0_r8_kind = coeff
                              cont_orb_fac, cont_grad_fac, cont_sd_fac, &
                              so, sog, sosd)
                      else ! i_l > 0
                         ! the standard case including the Y_lm contributions
                         sh   => sha%sh(i_ea)
                         shm  => sh%l(i_l)%m
                         m    = uas%m(i_fu)
                         shg  => sha%shg(i_ea)
                         shm_grad => shg%l(i_l)%m
                         shsd => sha%shs(i_ea)
                         shm_sd   => shsd%l(i_l)%m
                         call combine_rad_and_ylm(i_sa, uab%N_uncontracted_fcts, &
                              coef, i_step, i_pa, i_ea, &              !! SB: added i_pa, previous 1
                              prim_orb_fac, prim_grad_fac, prim_sd_fac, &
                              so, sog, sosd)
                         ! accumulate contractions in fcts, grads, ...
                         call combine_rad_and_ylm(i_sa, uab%N_contracted_fcts,   &
                              coef, i_step, i_pa, i_ea, &              !! SB: added i_pa, previous 1
                              cont_orb_fac, cont_grad_fac, cont_sd_fac, &
                              so, sog, sosd)
                      end if !! i_l .eq. 0

                   enddo i_fu_ ! contributing functions
                enddo i_if_  ! independent fcts

             enddo i_partner_ ! nimber of partners
             indx(i_ir) = indx(i_ir) + uap%N_independent_fcts*(uab%N_uncontracted_fcts + uab%N_contracted_fcts)

          enddo i_ir_ ! number of irreps
       enddo i_bas_ ! i_bas or i_l
    enddo i_ua_ ! unique atom

    ! deallocate working arrays
    if (do_sec_ders) then
       deallocate( help_xx, help_xy, help_xz, &
            help_yy, help_yz, help_zz, stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help_{xx,xy,xz,yy,yz,zz} failed")
    end if
    if (do_grads) then
       deallocate( help_x, help_y, help_z, &
            help_g, help_gy, stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help_{x,y,z} and help_g{,y} failed")
    end if
    if (do_fcts) then
       deallocate( help, stat=as )
       if (as /= 0) call error_handler( &
            "fit_fct_calculate: deallocate help failed")
    end if
    deallocate( temp, stat=as )
    if (as /= 0) call error_handler( &
         "fit_fct_calculate: deallocate temp failed")

  contains

    subroutine combine_rad_and_ylm(i_sa, N_fcts, coeff, i_incr, i_pa, i_ea, &
         rad, rad_g, rad_sd, o, og, osd)
      ! Purpose : evaluation of the primitive Gaussians and their various
      !           derivatives on a set of grid points
      ! ----- Declaration of formal parameters -----------------------------
      ! alternating control parameter
      integer(kind=i4_kind), intent(inout) :: i_sa ! index of the S-A function to start with
      integer(kind=i4_kind), intent(in) :: N_fcts
      real   (kind=r8_kind), intent(in) :: coeff
      integer(kind=i4_kind), intent(in) :: i_incr
      integer(kind=i4_kind), intent(in) :: i_pa  !! partner!?
      integer(kind=i4_kind), intent(in) :: i_ea   !! equivalent_atoms
      ! the radial functions to process
      real(kind=r8_kind), intent(in) :: rad   (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_g (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_sd(:,:,:)
      ! the result types (in arbitrary occurance)
      real(kind=r8_kind), optional, pointer :: o   (:,:,:)
      real(kind=r8_kind), optional, pointer :: og  (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: osd (:,:,:,:)
      ! ----- Declaration of local variables -------------------------------
      integer(kind=i4_kind) :: i_fct
      !---------------------------------------------------------------------
      !               phi(r) = C R(r) Y_lm(r)
      !        d/dr_i phi(r) = C [1/r d/dr] R(r) r_i Y_lm(r)
      !                      + C R(r) d/dr_i Y_lm(r)
      ! d/dr_i d/dr_j phi(r) = C [1/r d/dr]^2 R(r) r_j r_i Y_lm(r)
      !                      + C [1/r d/dr] R(r) delta_ij Y_lm(r)
      !                      + C [1/r d/dr] R(r) r_i d/dr_j Y_lm(r)
      !                      + C [1/r d/dr] R(r) r_j d/dr_i Y_lm(r)
      !                      + C R(r) d/dr_j d/dr_i Y_lm(r)
      !------------ Executable code --------------------------------

!!$print *,'ylm: 1 sa=', i_sa
      do i_fct = 1,N_fcts

         if (do_fcts) then
            ! help(r) := C * R(r)
            help = coeff * rad(:vl,i_fct,i_ea)
            if (associated(o)) then
               o(:vl,i_sa,i_pa) = o(:vl,i_sa,i_pa) + &
                    help * shm(:vl,m)
            end if
         end if

         if (do_grads) then
            ! help_g(r)   := C * [1/r d/dr] R(r)
            ! help_gy(r)  := C * [1/r d/dr] R(r) * Y_lm(r)
            ! help_r_i(r) := d/dr_i phi(r)
            help_g  = coeff * rad_g(:vl,i_fct,i_ea)
            help_gy = help_g * shm(:vl,m)
            help_x  = help_gy * x(:vl) + &
                 help * shm_grad(:vl,1,m)
            help_y  = help_gy * y(:vl) + &
                 help * shm_grad(:vl,2,m)
            help_z  = help_gy * z(:vl) + &
                 help * shm_grad(:vl,3,m)
            if (associated(og)) then
               og(:vl,1,i_sa,i_pa) = og(:vl,1,i_sa,i_pa) + help_x
               og(:vl,2,i_sa,i_pa) = og(:vl,2,i_sa,i_pa) + help_y
               og(:vl,3,i_sa,i_pa) = og(:vl,3,i_sa,i_pa) + help_z
            end if
         end if

         if (do_sec_ders) then
            ! help_r_ir_j(r) := d/dr_i d/dr_j phi(r)
            help_xx = help * shm_sd(:vl,XX,m) + help_gy
            help_xy = help * shm_sd(:vl,XY,m)
            help_xz = help * shm_sd(:vl,XZ,m)
            help_yy = help * shm_sd(:vl,YY,m) + help_gy
            help_yz = help * shm_sd(:vl,YZ,m)
            help_zz = help * shm_sd(:vl,ZZ,m) + help_gy
            ! help_r_i(r)    := C * [1/r d/dr] R(r) * r_i
            help_x  = help_g * x(:vl)
            help_y  = help_g * y(:vl)
            help_z  = help_g * z(:vl)
            help_xx = help_xx + help_x * shm_grad(:vl,1,m) * 2.0_r8_kind
            help_xy = help_xy + help_x * shm_grad(:vl,2,m) &
                 + help_y * shm_grad(:vl,1,m)
            help_xz = help_xz + help_x * shm_grad(:vl,3,m) &
                 + help_z * shm_grad(:vl,1,m)
            help_yy = help_yy + help_y * shm_grad(:vl,2,m) * 2.0_r8_kind
            help_yz = help_yz + help_y * shm_grad(:vl,3,m) &
                 + help_z * shm_grad(:vl,2,m)
            help_zz = help_zz + help_z * shm_grad(:vl,3,m) * 2.0_r8_kind
            ! help_gy(r)     := C * [1/r d/dr]^2 R(r) * Y_lm(r)
            ! help_r_i(r)    := C * [1/r d/dr]^2 R(r) * r_i * Y_lm(r)
            help_gy = coeff * rad_sd(:vl,i_fct,i_ea) * shm(:vl,m)
            help_x  = help_gy * x(:vl)
            help_y  = help_gy * y(:vl)
            help_z  = help_gy * z(:vl)
            help_xx = help_xx + help_x * x(:vl)
            help_xy = help_xy + help_x * y(:vl)
            help_xz = help_xz + help_x * z(:vl)
            help_yy = help_yy + help_y * y(:vl)
            help_yz = help_yz + help_y * z(:vl)
            help_zz = help_zz + help_z * z(:vl)
            if (associated(osd)) then
               osd(:vl,XX,i_sa,i_pa) = osd(:vl,XX,i_sa,i_pa) + help_xx(:vl)
               osd(:vl,XY,i_sa,i_pa) = osd(:vl,XY,i_sa,i_pa) + help_xy(:vl)
               osd(:vl,XZ,i_sa,i_pa) = osd(:vl,XZ,i_sa,i_pa) + help_xz(:vl)
               osd(:vl,YY,i_sa,i_pa) = osd(:vl,YY,i_sa,i_pa) + help_yy(:vl)
               osd(:vl,YZ,i_sa,i_pa) = osd(:vl,YZ,i_sa,i_pa) + help_yz(:vl)
               osd(:vl,ZZ,i_sa,i_pa) = osd(:vl,ZZ,i_sa,i_pa) + help_zz(:vl)

!!$               osd(:vl,2,1,i_sa,i_pa) = osd(:vl,2,1,i_sa,i_pa) + help_xy
!!$               osd(:vl,3,1,i_sa,i_pa) = osd(:vl,3,1,i_sa,i_pa) + help_xz
!!$               osd(:vl,3,2,i_sa,i_pa) = osd(:vl,3,2,i_sa,i_pa) + help_yz
            end if
         end if

!!$         i_sa = i_sa + i_incr
!!$         i_sa = i_sa + 1
         i_sa = i_sa + i_incr
      end do
!!$print *,'ylm: 2 sa=', i_sa

    end subroutine combine_rad_and_ylm

    subroutine combine_rad_and_y00(i_sa, N_fcts, coeff, i_incr, i_pa, i_ea, &
         rad, rad_g, rad_sd, o, og, osd)
      ! Purpose : evaluation of the primitive Gaussians and their various
      !           derivatives on a set of grid points
      ! ----- Declaration of formal parameters -----------------------------
      ! alternating control parameter
      integer(kind=i4_kind), intent(inout) :: i_sa ! index of the S-A function to start with
      integer(kind=i4_kind), intent(in) :: N_fcts
      real   (kind=r8_kind), intent(in) :: coeff
      integer(kind=i4_kind), intent(in) :: i_incr !! SB: i_incr = #_ind_fcts
      integer(kind=i4_kind), intent(in) :: i_pa   !! SB: partners
      integer(kind=i4_kind), intent(in) :: i_ea   !! SB: equivalent_atoms
      ! the radial functions to process
      real(kind=r8_kind), intent(in) :: rad   (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_g (:,:,:)
      real(kind=r8_kind), intent(in) :: rad_sd(:,:,:)
      ! the result types (in arbitrary occurance)
      real(kind=r8_kind), optional, pointer :: o   (:,:,:)
      real(kind=r8_kind), optional, pointer :: og  (:,:,:,:)
      real(kind=r8_kind), optional, pointer :: osd (:,:,:,:)
      ! ----- Declaration of local variables -------------------------------
      integer(kind=i4_kind) :: i_fct
      !---------------------------------------------------------------------
      !               phi(r) = C R(r)
      !        d/dr_i phi(r) = C [1/r d/dr] R(r) r_i
      ! d/dr_i d/dr_j phi(r) = C [1/r d/dr]^2 R(r) r_j r_i
      !                      + C [1/r d/dr] R(r) delta_ij
      !------------ Executable code --------------------------------

!!$print *,'y00: 1 sa=', i_sa
      do i_fct = 1,N_fcts

         if (do_fcts) then
            if (associated(o)) then
               o(:vl,i_sa,i_pa) = o(:vl,i_sa,i_pa) + &
                    coeff * rad(:vl,i_fct,i_ea)
            end if
         end if

         if (do_grads) then
            ! help_g(r)   := C * [1/r d/dr] R(r)
            ! help_r_i(r) := d/dr_i phi(r)
            help_g  = coeff * rad_g(:vl,i_fct,i_ea)
            help_x  = help_g * x(:vl)
            help_y  = help_g * y(:vl)
            help_z  = help_g * z(:vl)
            if (associated(og)) then
               og(:vl,1,i_sa,i_pa) = og(:vl,1,i_sa,i_pa) + help_x
               og(:vl,2,i_sa,i_pa) = og(:vl,2,i_sa,i_pa) + help_y
               og(:vl,3,i_sa,i_pa) = og(:vl,3,i_sa,i_pa) + help_z
            end if
         end if

         if (do_sec_ders) then
            ! help(r)        := C * [1/r d/dr]^2 R(r)
            ! help_r_i(r)    := C * [1/r d/dr]^2 R(r) * r_i
            ! help_r_ir_j(r) := d/dr_i d/dr_j phi(r)
            help    = coeff * rad_sd(:vl,i_fct,i_ea)
            help_x  = help * x(:vl)
            help_y  = help * y(:vl)
            help_z  = help * z(:vl)
            help_xx = help_x * x(:vl) + help_g
            help_xy = help_x * y(:vl)
            help_xz = help_x * z(:vl)
            help_yy = help_y * y(:vl) + help_g
            help_yz = help_y * z(:vl)
            help_zz = help_z * z(:vl) + help_g
            if (associated(osd)) then
               osd(:vl,XX,i_sa,i_pa) = osd(:vl,XX,i_sa,i_pa) + help_xx
               osd(:vl,XY,i_sa,i_pa) = osd(:vl,XY,i_sa,i_pa) + help_xy
               osd(:vl,XZ,i_sa,i_pa) = osd(:vl,XZ,i_sa,i_pa) + help_xz
               osd(:vl,YY,i_sa,i_pa) = osd(:vl,YY,i_sa,i_pa) + help_yy
               osd(:vl,YZ,i_sa,i_pa) = osd(:vl,YZ,i_sa,i_pa) + help_yz
               osd(:vl,ZZ,i_sa,i_pa) = osd(:vl,ZZ,i_sa,i_pa) + help_zz
!!$               osd(:vl,2,1,i_sa,1) = osd(:vl,2,1,i_sa,1) + help_xy
!!$               osd(:vl,3,1,i_sa,1) = osd(:vl,3,1,i_sa,1) + help_xz
!!$               osd(:vl,3,2,i_sa,1) = osd(:vl,3,2,i_sa,1) + help_yz
            end if
         end if

!!$         i_so = i_so + 1
         !!i_sa = i_sa + 1
         i_sa = i_sa + i_incr
      end do
!!$print *,'y00: 2 sa=', i_sa

    end subroutine combine_rad_and_y00
  end subroutine fit_fct_calculate_response
  !*************************************************************

  !*************************************************************
  integer(i4_kind) function dimension_of_fit_ch(i_ir)
    ! Purpose: calculate the dimension of the irrep i_ir
    !          for charge fit basis.
    !
    ! based on dimensions_calc from fit_coeff_module.f90
    !
    !** End of interface **************************************
    !------------ Modules used --------------------------------
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind),intent(in) :: i_ir
     !------------ Declaration of local handles -------------------
    integer(kind=i4_kind)  :: i_ua,i_l
    integer(kind=i4_kind)  :: nff
    type(unique_atom_basis_type),pointer :: uab
    type(unique_atom_partner_type),pointer :: uap
    integer(i4_kind) :: counter
    type(unique_atom_type), pointer :: ua
    !------------ Executable code --------------------------------

    ! charge fitfunctions
    counter = 0
    do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       do i_l = -1, ua%lmax_ch ! s, r2, p, ...
          select case(i_l)
          case (-1) ! s-
             uab => ua%l_ch(0)
             uap => ua%symadapt_partner(i_ir,0)
          case (0)  ! r2-
             uab => ua%r2_ch
             uap => ua%symadapt_partner(i_ir,0)
          case (1:) ! p-, d-, ...
             uab => ua%l_ch(i_l)
             uap => ua%symadapt_partner(i_ir,i_l)
          case default
             ABORT('no such case')
          end select
          nff = (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
               * uap%N_independent_fcts

          counter = counter + nff

       enddo
    enddo
    dimension_of_fit_ch = counter

  end function dimension_of_fit_ch
  !*************************************************************

  integer(i4_kind) function orb_position(i_ir, ua, ll)
    use unique_atom_module, only: unique_atoms
    implicit none
    integer(kind = i4_kind), intent(in) :: ua   ! unique atom
    integer(kind = i4_kind), intent(in) :: i_ir ! irrep
    integer(kind = i4_kind), intent(in) :: ll
    !------------ Declaration of local types ---------------------
    integer(kind = i4_kind) :: loop, ua_loop
    integer(kind = i4_kind) :: ncntr, nunct, nifct

    orb_position = 0

    do ua_loop = 1, ua-1

       do loop = 0, unique_atoms(ua_loop)%lmax_ob !!CC
          ncntr = unique_atoms(ua_loop)%l_ob(loop)%N_contracted_fcts
          nunct = unique_atoms(ua_loop)%l_ob(loop)%N_uncontracted_fcts
          nifct = unique_atoms(ua_loop)%symadapt_partner(i_ir,loop)%N_independent_fcts
          orb_position = orb_position + (ncntr + nunct) * nifct
       end do

    end do

    do loop = 0, ll
       ncntr = unique_atoms(ua)%l_ob(loop)%N_contracted_fcts
       nunct = unique_atoms(ua)%l_ob(loop)%N_uncontracted_fcts
       nifct = unique_atoms(ua)%symadapt_partner(i_ir,loop)%N_independent_fcts
       orb_position = orb_position + (ncntr + nunct) * nifct
    end do

  end function orb_position

  integer(i4_kind) function fit_position(i_ir, ua, ll)
    use unique_atom_module, only: unique_atoms
    implicit none
    integer(kind = i4_kind), intent(in) :: ua   ! unique atom
    integer(kind = i4_kind), intent(in) :: i_ir ! irrep
    integer(kind = i4_kind), intent(in) :: ll
    !------------ Declaration of local types ---------------------
    integer(kind = i4_kind) :: i_l, i_ua
    type(unique_atom_basis_type),  pointer :: uab
    type(unique_atom_partner_type),pointer :: uap
    type(unique_atom_type),        pointer :: uam

    fit_position = 0

    ! charge fitfunctions
    do i_ua = 1, ua-1
       uam => unique_atoms(i_ua)
       do i_l = -1, uam%lmax_ch ! s, r2, p, ...
          select case(i_l)
          case (-1) ! s-
             uab => uam%l_ch(0)
             uap => uam%symadapt_partner(i_ir,0)
          case (0)  ! r2-
             uab => uam%r2_ch
             uap => uam%symadapt_partner(i_ir,0)
          case (1:) ! p-, d-, ...
             uab => uam%l_ch(i_l)
             uap => uam%symadapt_partner(i_ir,i_l)
          case default
             ABORT('no such case')
          end select
          fit_position = fit_position &
               &       + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
               &       * uap%N_independent_fcts
       end do
    end do

    uam => unique_atoms(ua)
    do i_l = -1, ll
       select case(i_l)
       case (-1) ! s-
          uab => uam%l_ch(0)
          uap => uam%symadapt_partner(i_ir,0)
       case (0)  ! r2-
          uab => uam%r2_ch
          uap => uam%symadapt_partner(i_ir,0)
       case (1:) ! p-, d-, ...
          uab => uam%l_ch(i_l)
          uap => uam%symadapt_partner(i_ir,i_l)
       case default
          ABORT('no such case')
       end select
       fit_position = fit_position &
            &       + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
            &       * uap%N_independent_fcts
    end do

  end function fit_position


!!  !*************************************************************
!!  subroutine <name>_
!!    !  Purpose: ..
!!    !------------ Modules used ------------------- ---------------
!!    use
!!    implicit none
!!    !------------ Declaration of formal parameters ---------------
!!    integer(kind=i4_kind), intent(     ) ::
!!    real(kind=r8_kind),    intent(     ) ::
!!    logical,               intent(     ) ::
!!    character,             intent(     ) ::
!!    !** End of interface *****************************************
!!    !------------ Declaration of local variables -----------------
!!    integer(kind=i4_kind)                ::
!!    real(kind=r8_kind)                   ::
!!    logical                              ::
!!    character                            ::
!!    !------------ Executable code --------------------------------


!!  end subroutine <name>_
!!  !*************************************************************


  !--------------- End of module ----------------------------------
end module ch_response_module

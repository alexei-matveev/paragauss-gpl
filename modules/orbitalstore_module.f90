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
module orbitalstore_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains the orbitals calculated on grid
  !           &MAIN_OPTIONS :
  !           ORBITALS_IN_MEMORY = .TRUE. (orbitals are kept in memory)
  !           ORBITALS_ON_FILE   = .TRUE. (       ="=        on file  )
  !           The part of the code is duplicated to minimize changes
  !           in the "orbital_module"
  !
  !  Module called by: orbital_calculate
  !                    [xc_hamiltonian, xcmda_hamiltonian,
  !                     xcfit_hamiltonian]
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
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
  use type_module            ! type specification parameters
  use options_module, only : options_orbitals_in_memory, &
                             options_orbitals_on_file
  use symmetry_data_module
  use unique_atom_module
  USE_MEMLOG
  implicit none
  save                       ! save all variables defined in this module
  private                    ! by default, all names are private

  !------------ Declaration of public types -----------------------
  !             (moved from orbital_module)

  ! xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,lin(ua,i),prt)
  type, public :: orbital_type
     ! to contain calculated sym. orbitals of one irrep i
     real(kind=r8_kind), pointer :: o(:,:,:)
     ! o(vec_len, symmetry_data_dimension(i), symmetry_data_n_partners(i))
     ! symertry-adapted molecular orbitals on vec_len gridpoints
!    ! ordering within dimension see orbitalprojection_module
!    real(kind=r8_kind), pointer :: o_imag(:,:,:)
!    ! imaginary part of orbitals
  end type orbital_type

  type, public :: spinor_type
!    ! contains two component orbitals of one irrep i
!    type(orbital_type) :: spinor(2)
     ! actual storage:
     real(kind=r8_kind), pointer :: o(:,:,:,:,:)
     ! o(vec_len,bas_dim,n_partners,2,2)
     ! spinor(1)%o      => o(:,:,:,1,1)
     ! spinor(1)%o_imag => o(:,:,:,2,1)
     ! spinor(2)%o      => o(:,:,:,1,2)
     ! spinor(2)%o_imag => o(:,:,:,2,2)
  end type spinor_type

  type, public :: orbital_spin_type ! for the response part
     ! to contain calculated sym. orbitals of one irrep irreps(i)
     real(kind=r8_kind), pointer :: o(:,:,:,:)
     ! o(vec_len,symmetry_data_dimension(i),symmetry_data_n_partners(i),n_spin)
     ! symertry-adapted molecular orbitals on vec_len gridpoints
     ! ordering within dimension see orbitalprojection_module
  end type orbital_spin_type

  ! d/dr_n xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,lin(ua,i),prt)
  type, public :: orbital_gradient_type
     ! to contain calculated sym. orbital gradients of one irrep i
     real(kind=r8_kind), pointer :: o(:,:,:,:)
     ! o(vec_len, 3, symmetry_data_dimension(i), symmetry_data_n_partners(i))
     ! symertry-adapted molecular orbital gradients  on vec_len gridpoints
     ! ordering within dimension see orbitalprojection_module
!    real(kind=r8_kind), pointer :: o_imag(:,:,:,:)
  end type orbital_gradient_type

  ! d/dr_n xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,lin(ua,i),prt)
  type, public :: spinor_gradient_type
!    ! contains two component orbitals gradients of one irrep i
!    type(orbital_gradient_type) :: spinor(2)
     ! actual storage:
     real(kind=r8_kind), pointer :: o(:,:,:,:,:,:) ! (vl,3,NBAS,NPART,RE:IM,UP:DN)
     ! spinor(1)%o      => o(:,:,:,:,1,1)
     ! spinor(1)%o_imag => o(:,:,:,:,2,1)
     ! spinor(2)%o      => o(:,:,:,:,1,2)
     ! spinor(2)%o_imag => o(:,:,:,:,2,2)
  end type spinor_gradient_type

  ! d/dr_n d/dr_m xi[ua,irr,prt]_i(r_a) = <var>(irr)%o(a,n,m,lin(ua,i),prt)
  type, public :: orbital_sec_der_type
     ! to contain calculated sym. orbital secound derivatives of one irrep i
     real(kind=r8_kind), pointer :: o(:,:,:,:)
     ! o(vec_len, 6, symmetry_data_dimension(i), symmetry_data_n_partners(i))
     ! secound dimension: 1:6 double derivatives -- xx xy xz yy yx zz
     ! symertry-adapted molecular orbital secound derivatives on vec_len gridpoints
     ! ordering within dimension see orbitalprojection_module
  end type orbital_sec_der_type

  ! [-d/dR_n] xi[ua,irr,prt]_i(r_a) = <var>(irr,ua)%o(a,n,R,i,prt)
  type, public :: orbital_nuclear_gradient_type
     ! to contain calculated gradient with respect
     ! to nuclear displacement of sym. orbitals of one irrep irreps(i)
     real(kind=r8_kind), pointer :: o(:,:,:,:,:)
     ! o(vec_len, 3,n_equal_atoms,n_orbitals,
     !       symmetry_data_n_partners(i))
     ! molecular orbital gradients  on vec_len gridpoints with respect to
     ! nuclear displacements
     ! n_orbitals: Number of orbitals  per u_nique_atoms
  end type orbital_nuclear_gradient_type

  ! d/dR_n d/dR_m xi[ua,irr,prt]_i(r_a) = <var>(irr,ua)%o(a,lin(n,m),R,i,prt)
  type, public :: orbital_nuclear_sec_der_type
     ! to contain calculated 2nd derivatives with respect
     ! to nuclear displacement of sym. orbitals of one irrep irreps(i)
     real(kind=r8_kind), pointer :: o(:,:,:,:,:)
     ! o(vec_len, 6,n_equal_atoms,n_orbitals,
     !       symmetry_data_n_partners(i))
     ! n_orbitals: Number of orbitals  per u_nique_atoms
     ! only the upper half of the hessian matrix is stored
  end type orbital_nuclear_sec_der_type

  ! reduced variants of the orbital types used for core density fitting fcts.
  ! xi[ua]_i(r_a) = <var>%o(a,lin(ua,i))
  type, public :: core_orbital_type
     ! to contain the length of grid points
     ! o(vec_len,fit%n_cd)
     real(kind=r8_kind), pointer :: o(:,:)
  end type core_orbital_type

  ! d/dr_n xi[ua]_i(r_a) = <var>%o(a,n,lin(ua,i))
  type, public :: core_orbital_gradient_type
     ! to contain the length of grid points and coordinates
     ! o(vec_len,3,fit%n_cd)
     real(kind=r8_kind), pointer :: o(:,:,:)
  end type core_orbital_gradient_type

  ! d/dr_n d/dr_m xi[ua]_i(r_a) = <var>%o(a,n,m,lin(ua,i))
  type, public :: core_orbital_sec_der_type
     ! to contain the length of grid points and two coordinates
     ! o(vec_len,6,fit%n_cd)
     real(kind=r8_kind), pointer :: o(:,:,:)
     ! secound dimension: 1-6: xx xy xz yy yz zz
  end type core_orbital_sec_der_type

  ! [-d/dR_n] xi[ua]_i(r_a) = <var>(ua)%o(a,n,R,i)
  type, public :: core_orbital_nuc_gradient_type
     ! to contain the length of grid points
     ! o(vec_len, 3,n_equal_atoms,fit%n_cd)
     real(kind=r8_kind), pointer :: o(:,:,:,:)
     ! molecular orbital gradients  on vec_len gridpoints with respect to
     ! nuclear displacements
  end  type core_orbital_nuc_gradient_type

  ! d/dR_n d/dR_m xi[ua]_i(r_a) = <var>(ua)%o(a,lin(n,m),R,i)
  type, public :: core_orbital_nuc_sec_der_type
     ! to contain the length of grid points
     ! o(vec_len,6,n_equal_atoms,fit%n_cd)
     real(kind=r8_kind), pointer :: o(:,:,:,:)
     ! n_orbitals: Number of orbitals  per u_nique_atoms
     ! only the upper half of the hessian matrix is stored
  end type core_orbital_nuc_sec_der_type

  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------
  type orbitalstore_type
     type(orbital_type),                  pointer :: orbs_ob(:)
     type(orbital_type)                           :: fcts
     type(spinor_type),                   pointer :: orbs_spinor_ob(:)
     type(orbital_gradient_type),         pointer :: grads(:)
     type(orbital_gradient_type)                  :: fcts_grads
     type(spinor_gradient_type),          pointer :: spinor_grads(:)
     type(orbital_sec_der_type),          pointer :: sec_der(:)
     type(orbital_nuclear_gradient_type), pointer :: nuc_grads(:,:)
     type(orbital_spin_type),             pointer :: phi_ob(:)
     type(orbital_nuclear_sec_der_type),  pointer :: nuc_sec_der(:,:)
     type(orbital_type),                  pointer :: fcts_ch(:) !! SB for response
  end type orbitalstore_type

!!$  type chanel_type
!!$     integer(kind=i4_kind) :: orbs_ob
!!$     integer(kind=i4_kind) :: grads
!!$     integer(kind=i4_kind) :: spinor_grads
!!$     integer(kind=i4_kind) :: sec_der
!!$     integer(kind=i4_kind) :: nuc_grads
!!$     integer(kind=i4_kind) :: nuc_sec_der
!!$     integer(kind=i4_kind) :: phi_ob
!!$     integer(kind=i4_kind) :: orbs_spinor_ob
!!$  end type chanel_type

  integer(kind=i4_kind), parameter ::&
       & OS_ORBS_OB        =  1, &
       & OS_FCTS           =  2, &
       & OS_ORBS_SPINOR_OB =  3, &
       & OS_GRADS          =  4, &
       & OS_FCTS_GRADS     =  5, &
       & OS_SPINOR_GRADS   =  6, &
       & OS_SEC_DER        =  7, &
       & OS_NUC_GRADS      =  8, &
       & OS_PHI_OB         =  9, &
       & OS_NUC_SEC_DER    = 10, &
       & OS_N_OPTIONS      = 10
  logical :: using(OS_N_OPTIONS) = .false.

  real(r8_kind), parameter :: cutoff=1.0E-10_r8_kind ! see how much of storage is nonzero
  integer(i4_kind)         :: counts(2,OS_N_OPTIONS)

  !------------ Declaration of constants and variables ------------
  integer(i4_kind) :: orbitalstore_mem
  logical, public, protected :: orbitalstore_initialized = .false., &
                      stored_orbs_ob         , &
                      stored_fcts            , &
                      stored_orbs_spinor_ob  , &
                      stored_grads           , &
                      stored_fcts_grads      , &
                      stored_spinor_grads    , &
                      stored_sec_der         , &
                      stored_nuc_grads       , &
                      stored_phi_ob          , &
                      stored_nuc_sec_der

  !------------ public functions and subroutines ------------------
  public :: orbitalstore_setup,    orbitalstore_shutdown, &
            orbitalstore_allocate, orbitalstore_free,  &
            orbitalstore_rw

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------
  type(orbitalstore_type), private, allocatable    :: orbitalstore(:)
!!$  type(chanel_type),                     parameter :: chanel = (1,2,3,4,5,6,7,8)
!!$  character(len=15),       dimension(8), parameter :: ch_name = (/&
!!$                                                      "orbs_ob        ",&
!!$                                                      "grads          ",&
!!$                                                      "spinor_grads   ",&
!!$                                                      "sec_der        ",&
!!$                                                      "nuc_grads      ",&
!!$                                                      "nuc_sec_der    ",&
!!$                                                      "phi_ob         ",&
!!$                                                      "orbs_spinor_ob "/)

  !------------ Declaration of constants and variables ----
  integer(kind=i4_kind) :: vec_len
  integer(kind=i4_kind) :: &
       N_grid_blocks,         &
       grid_block_orbs_ob,    &
       grid_block_fcts,       &
       grid_block_grads,      &
       grid_block_fcts_grads, &
       grid_block_sec_der,    &
       grid_block_nuc_grads,  &
       grid_block_nuc_sec_der,&
       grid_block_phi_ob,     &
       grid_block_orbs_spinor_ob, &
       grid_block_spinor_grads

  !------------ Subroutines ---------------------------------------
contains


  !****************************************************************
  subroutine orbitalstore_setup(vector_length)
    ! Purpose : allocates memory or open files to save
    !           the orbitals on grid
    !           Called by xc_setup : xcfit_setup : xcmda_setup
    use grid_module, only : more_grid, grid_loop_setup
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: vector_length
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind),dimension(:,:),pointer :: grdpts
    real(kind=r8_kind),dimension(:),  pointer :: grdwts
    integer(kind=i4_kind) :: alloc_stat
    !------------ Executable code --------------------------------
    if(.not.orbitalstore_initialized) then
       vec_len = vector_length
       ! simply to avoid writting new function in "grid_module"

       ! nothing is used yet:
       using(:) = .false.

       N_grid_blocks=0

       !FIXME: won't work with dynamic load balancing, no way of knowing bevore
       !       how many are there
       call grid_loop_setup()
       ! count number of grid batches:
       do while( more_grid(vec_len,grdpts,grdwts) ) ! fetching part of the grid
          ASSERT(associated(grdpts))
          ASSERT(size(grdpts,1)>0)

          N_grid_blocks = N_grid_blocks + 1
       end do

       if(options_orbitals_in_memory() ) then
          !FIXME: this won't have the size needed for all the jobs
          !       to do on the grids, as dynamic load balancing might
          !       change the number, should be aborted earlier already,
          !       but just in case:
          call error_handler("orbitalstore_init : orbitals_in_memory broken")
          allocate( orbitalstore(N_grid_blocks), stat=alloc_stat)
          if(alloc_stat /= 0) call error_handler &
               ("orbitalstore_init_orbs_ob : allocation orbitalstore failed")
          orbitalstore_initialized = .true.
          orbitalstore_mem         = 0
          grid_block_orbs_ob       = 0;  stored_orbs_ob     = .false.
          grid_block_fcts          = 0;  stored_fcts        = .false.
          grid_block_grads         = 0;  stored_grads       = .false.
          grid_block_fcts_grads    = 0;  stored_fcts_grads  = .false.
          grid_block_sec_der       = 0;  stored_sec_der     = .false.
          grid_block_nuc_grads     = 0;  stored_nuc_grads   = .false.
          grid_block_nuc_sec_der   = 0;  stored_nuc_sec_der = .false.
          grid_block_phi_ob        = 0;  stored_phi_ob      = .false.
          grid_block_orbs_spinor_ob= 0;  stored_orbs_spinor_ob = .false.
          grid_block_spinor_grads  = 0;  stored_spinor_grads   = .false.
       else if( options_orbitals_on_file() ) then
          call error_handler("orbitalstore_init : ORBITLAS_ON_FILE not yet implemented")
          orbitalstore_initialized = .true.
       end if
    end if

  end subroutine orbitalstore_setup
  !****************************************************************

  !****************************************************************
  subroutine orbitalstore_shutdown()
    ! Purpose : deallocate memory or close files for grid-orbitals
    ! Called by : xc_close : xcfit_close : xcmda_close
    use error_module, only: MyID
    implicit none
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: alloc_stat
    !------------ Executable code --------------------------------

    if( .not. orbitalstore_initialized ) RETURN

       counts = 0
       call orbitalstore_free()
       if( options_orbitals_in_memory() ) then
          deallocate( orbitalstore, stat=alloc_stat)
          if(alloc_stat /= 0) call error_handler &
               ("orbitalstore_free : deallocation of orbitalstore failed")
       else if( options_orbitals_on_file() ) then
          call error_handler("orbitalstore_free : ORBITLAS_ON_FILE not yet implemented")
       end if
       print *,MyID,'orbitalstore:', orbitalstore_mem, ' b used ' &
                   , counts(2,OS_ORBS_OB),' of ', counts(1,OS_ORBS_OB)        &
                   , ' below ', cutoff
       orbitalstore_initialized = .false.
  end subroutine orbitalstore_shutdown
  !****************************************************************

  !****************************************************************
  subroutine orbitalstore_allocate(orbs_ob,grads,sec_der,nuc_grads,  &
                               nuc_sec_der,phi_ob,orbs_spinor_ob,&
                               spinor_grads,fcts,fcts_grads,fcts_dim)
    !  Purpose: allocating memory or open file to save
    !           the orbitals on the grid
    !           Called by "xc_hamiltonian" ...
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(orbital_type),                  pointer, optional :: orbs_ob(:)        !
    type(orbital_gradient_type),         pointer, optional :: grads(:)          !
    type(spinor_gradient_type),          pointer, optional :: spinor_grads(:)
    type(orbital_sec_der_type),          pointer, optional :: sec_der(:)        !
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:,:)    !
    type(orbital_nuclear_sec_der_type),  pointer, optional :: nuc_sec_der(:,:)  !
    type(orbital_spin_type),             pointer, optional :: phi_ob(:)         !
    type(spinor_type),                   pointer, optional :: orbs_spinor_ob(:)
    type(orbital_type),                           optional :: fcts
    type(orbital_gradient_type),                  optional :: fcts_grads
    integer(kind=i4_kind),            intent(in), optional :: fcts_dim
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind),  parameter            :: zero = 0.0_r8_kind
    integer(kind=i4_kind)                     :: i_block, i_ir, i_unique, i_l, &
                                                 n_orbitals, alloc_stat
    logical                                   :: first_block
    integer(kind=i4_kind)                     :: mem
    !------------ Executable code --------------------------------
    DPRINT 'orbstore_alloc: initialized=', orbitalstore_initialized,' n_b=', N_grid_blocks

    if( present(orbs_ob)        ) using(OS_ORBS_OB)        = .true.
    if( present(grads)          ) using(OS_GRADS)          = .true.
    if( present(spinor_grads)   ) using(OS_SPINOR_GRADS)   = .true.
    if( present(sec_der)        ) using(OS_SEC_DER)        = .true.
    if( present(nuc_grads)      ) using(OS_NUC_GRADS)      = .true.
    if( present(nuc_sec_der)    ) using(OS_NUC_SEC_DER)    = .true.
    if( present(phi_ob)         ) using(OS_PHI_OB)         = .true.
    if( present(orbs_spinor_ob) ) using(OS_ORBS_SPINOR_OB) = .true.
    if( present(fcts)           ) using(OS_FCTS)           = .true.
    if( present(fcts_grads)     ) using(OS_FCTS_GRADS)     = .true.

    if(.not.orbitalstore_initialized) call error_handler &
         ("orbitalstore_allocate : Orbitalstore not initialized")
    if( options_orbitals_in_memory() ) then
       DPRINT 'orbstore_alloc: orbitals_in_memory!'
       grid_blocks : do i_block = 1,N_grid_blocks
          first_block = i_block==1

          if( present(orbs_ob) .and. .not.grid_block_orbs_ob==N_grid_blocks) then
             DPRINT 'orbstore_alloc: +orbs_ob'
             allocate( orbitalstore(i_block)%orbs_ob(symmetry_data_n_irreps()), &
                                                                stat=alloc_stat )
             if ( alloc_stat /= 0) call error_handler &
                  ( "orbitalstore_init : allocate of orbs_ob failed")
             do i_ir = 1, symmetry_data_n_irreps()
                allocate( orbitalstore(i_block)%                             &
                          orbs_ob(i_ir)%o( vec_len,                          &
                                           symmetry_data_dimension(i_ir),    &
                                           symmetry_data_n_partners(i_ir) ), &
                                                             stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ( "orbitalstore_init : allocate of orbs_ob%o failed")
                orbitalstore(i_block)% orbs_ob(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%orbs_ob(i_ir)%o)
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)
             end do
             if(first_block) orbs_ob => orbitalstore(i_block)%orbs_ob
          end if

          if( present(fcts) .and. .not.grid_block_fcts==N_grid_blocks )then
             DPRINT 'orbstore_alloc: +fcts'
             if(.NOT.present(fcts_dim)) call error_handler("N_DIM needed for fcts")
             allocate( orbitalstore(i_block)% &
                       fcts%o(vec_len,fcts_dim,1), stat=alloc_stat)
             orbitalstore(i_block)%fcts%o=zero
             mem = size(orbitalstore(i_block)%fcts%o)
             orbitalstore_mem = orbitalstore_mem + mem*8
             MEMLOG(mem)
             if (alloc_stat /= 0) call error_handler &
                  ("orbitalstore_allocate: allocate of fcts%o failed")
          end if

          if( present(grads) .and. .not.grid_block_grads==N_grid_blocks )then
             DPRINT 'orbstore_alloc: +grads'
             allocate(orbitalstore(i_block)%grads(symmetry_data_n_irreps()), stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler( &
                  "orbitalstore_init: allocate of grads failed")
             do i_ir = 1, symmetry_data_n_irreps()
                allocate( orbitalstore(i_block)%                          &
                          grads(i_ir)%o(vec_len,3,                        &
                                        symmetry_data_dimension(i_ir),    &
                                        symmetry_data_n_partners(i_ir) ), &
                                                          stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ("orbitalstore_init: allocate of grads%o failed")
                orbitalstore(i_block)%grads(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%grads(i_ir)%o)
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)
             enddo
             if(first_block) grads => orbitalstore(i_block)%grads
          end if

          if(present(fcts_grads) .and. .not.grid_block_fcts==N_grid_blocks ) then
             DPRINT 'orbstore_alloc: +fcts_grads'
             if(.NOT.present(fcts_dim)) call error_handler("N_dim needed for fcts")
             allocate( orbitalstore(i_block)% &
                       fcts_grads%o(vec_len,3,fcts_dim,1), stat=alloc_stat)
             orbitalstore(i_block)%fcts_grads%o = zero
             mem = size(orbitalstore(i_block)%fcts_grads%o)
             orbitalstore_mem = orbitalstore_mem + mem*8
             MEMLOG(mem)
             if (alloc_stat /= 0) call error_handler &
                  ("orbitalstore_allocate: allocate of grads%o failed")
          end if

          if ( present(sec_der) .and. .not.grid_block_sec_der==N_grid_blocks ) then
             DPRINT 'orbstore_alloc: +sec_der'
             allocate( orbitalstore(i_block)%sec_der(symmetry_data_n_irreps()), stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ( "orbitalstore_init :  allocate of sec_der failed")
             do i_ir = 1, symmetry_data_n_irreps()
                allocate( orbitalstore(i_block)%                             &
                          sec_der(i_ir)%o( vec_len,6,                        &
                                           symmetry_data_dimension(i_ir),    &
                                           symmetry_data_n_partners(i_ir) ), &
                                                             stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler&
                     ( "orbitalstore_init :  allocate of sec_der%o failed")
                orbitalstore(i_block)%sec_der(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%sec_der(i_ir)%o)
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)
             end do
             if(first_block) sec_der = orbitalstore(i_block)%sec_der
          endif

          if ( present(nuc_grads) .and. .not.grid_block_nuc_grads==N_grid_blocks ) then
             DPRINT 'orbstore_alloc: +nuc_grads'
             allocate( orbitalstore(i_block) %                             &
                       nuc_grads(n_unique_atoms,symmetry_data_n_irreps()), &
                       stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_init : allocate of nuc_grads failed")
             do i_ir = 1, symmetry_data_n_irreps()
                do i_unique=1,n_unique_atoms
                   n_orbitals=0
                   do i_l=0,unique_atoms(i_unique)%lmax_ob
                      n_orbitals=n_orbitals+unique_atoms(i_unique)%&
                           symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                           (unique_atoms(i_unique)%l_ob(i_l)%n_contracted_fcts+&
                           unique_atoms(i_unique)%l_ob(i_l)%n_uncontracted_fcts)
                   end do
                   allocate( orbitalstore(i_block)%                                           &
                             nuc_grads(i_unique,i_ir)%o( vec_len, 3 ,                         &
                                                         unique_atoms(i_unique)%n_equal_atoms,&
                                                         n_orbitals,                          &
                                                         symmetry_data_n_partners(i_ir) ),    &
                             stat=alloc_stat )
                   if ( alloc_stat /= 0 ) call error_handler &
                        ("orbitalstore_init: allocate of nuc_grads%o failed")
                   orbitalstore(i_block)%nuc_grads(i_unique,i_ir)%o = zero
                   mem = size(orbitalstore(i_block)%nuc_grads(i_unique,i_ir)%o)
                   orbitalstore_mem = orbitalstore_mem + mem*8
                   MEMLOG(mem)
                enddo
             end do
             if(first_block) nuc_grads => orbitalstore(i_block)%nuc_grads
          endif

          if ( present(nuc_sec_der) .and. .not.grid_block_nuc_sec_der==N_grid_blocks ) then
             DPRINT 'orbstore_alloc: +nuc_sec_der'
             allocate( orbitalstore(i_block)%                               &
                       nuc_sec_der(n_unique_atoms,symmetry_data_n_irreps()), &
                       stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler( &
                  "orbitalstore_init : allocate of nuc_sec_der failed")
             do i_ir = 1, symmetry_data_n_irreps()
                do i_unique=1,n_unique_atoms
                   n_orbitals=0
                   do i_l=0,unique_atoms(i_unique)%lmax_ob
                      n_orbitals=n_orbitals+unique_atoms(i_unique)%&
                           symadapt_partner(i_ir,i_l)%n_independent_fcts*&
                           (unique_atoms(i_unique)%l_ob(i_l)%n_contracted_fcts+&
                           unique_atoms(i_unique)%l_ob(i_l)%n_uncontracted_fcts)
                   end do
                   allocate( orbitalstore(i_block)%                            &
                             nuc_sec_der(i_unique,i_ir)%o( vec_len, 6 ,        &
                                         unique_atoms(i_unique)%n_equal_atoms, &
                                         n_orbitals,                           &
                                         symmetry_data_n_partners(i_ir) ),     &
                             stat=alloc_stat )
                   if ( alloc_stat /= 0 ) call error_handler( &
                        "orbitalstore_init : allocate of nuc_sec_der%o failed")
                   orbitalstore(i_block)%nuc_sec_der(i_unique,i_ir)%o = zero
                   mem = size(orbitalstore(i_block)%nuc_sec_der(i_unique,i_ir)%o)
                   orbitalstore_mem = orbitalstore_mem + mem*8
                   MEMLOG(mem)
                end do
             end do
             if(first_block) nuc_sec_der => orbitalstore(i_block)%nuc_sec_der
          end if

          if ( present(phi_ob) .and. .not.grid_block_phi_ob==N_grid_blocks ) then
             DPRINT 'orbstore_alloc: +phi_ob'
             allocate( orbitalstore(i_block)%phi_ob(symmetry_data_n_irreps()), &
                       stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_init : allocate of phi failed")
             do i_ir = 1, symmetry_data_n_irreps()
                allocate( orbitalstore(i_block)%                           &
                          phi_ob(i_ir)%o( vec_len,                         &
                                          symmetry_data_dimension(i_ir),   &
                                          symmetry_data_n_partners(i_ir) , &
                                          symmetry_data_n_spin() ),        &
                          stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler( &
                     "orbitastore_init : allocate of phi%o failed")
                orbitalstore(i_block)%phi_ob(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%phi_ob(i_ir)%o )
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)
             enddo
             if(first_block) phi_ob => orbitalstore(i_block)%phi_ob
          end if

          ! SPINORS, FOR SO ONLY:
          if( present(orbs_spinor_ob) .and. .not.grid_block_orbs_spinor_ob==N_grid_blocks) then
             DPRINT 'orbstore_alloc: +orbs_spinor_ob'
             allocate( orbitalstore(i_block)%orbs_spinor_ob(symmetry_data_n_proj_irreps()), &
                                                                stat=alloc_stat )
             ASSERT(alloc_stat==0)
             do i_ir = 1, symmetry_data_n_proj_irreps()
                allocate( orbitalstore(i_block)%                                  &
                          orbs_spinor_ob(i_ir)%o( vec_len,                        &
                                           symmetry_data_dimension_proj(i_ir),    &
                                           symmetry_data_n_partners_proj(i_ir),   &
                                           2, 2 ),                                &
                          stat=alloc_stat )
                ASSERT(alloc_stat==0)
                orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o)
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)

!               ! associate alternative addressing of components:
!               orbitalstore(i_block)%orbs_spinor_ob(i_ir)%spinor(1)%o_imag => &
!                    & orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o(:,:,:,2,1)
!               orbitalstore(i_block)%orbs_spinor_ob(i_ir)%spinor(2)%o      => &
!                    & orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o(:,:,:,1,2)
!               orbitalstore(i_block)%orbs_spinor_ob(i_ir)%spinor(2)%o_imag => &
!                    & orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o(:,:,:,2,2)
             end do
             if(first_block) orbs_spinor_ob => orbitalstore(i_block)%orbs_spinor_ob
          end if

          ! SPINOR GRADS, FOR SO ONLY:
          if( present(spinor_grads) .and. .not.grid_block_spinor_grads==N_grid_blocks) then
             DPRINT 'orbstore_alloc: +spinor_grads'
             allocate( orbitalstore(i_block)%spinor_grads(symmetry_data_n_proj_irreps()), &
                                                                stat=alloc_stat )
             ASSERT(alloc_stat==0)
             do i_ir = 1, symmetry_data_n_proj_irreps()
                allocate( orbitalstore(i_block)%                                  &
                          spinor_grads(i_ir)%o( vec_len, 3,                       &
                                           symmetry_data_dimension_proj(i_ir),    &
                                           symmetry_data_n_partners_proj(i_ir),   &
                                           2, 2 ),                                &
                          stat=alloc_stat )
                ASSERT(alloc_stat==0)
                orbitalstore(i_block)%spinor_grads(i_ir)%o = zero
                mem = size(orbitalstore(i_block)%spinor_grads(i_ir)%o)
                orbitalstore_mem = orbitalstore_mem + mem*8
                MEMLOG(mem)

!               ! associate alternative addressing of components:
!               orbitalstore(i_block)%spinor_grads(i_ir)%spinor(1)%o      => &
!                    & orbitalstore(i_block)%spinor_grads(i_ir)%o(:,:,:,:,1,1)
!               orbitalstore(i_block)%spinor_grads(i_ir)%spinor(1)%o_imag => &
!                    & orbitalstore(i_block)%spinor_grads(i_ir)%o(:,:,:,:,2,1)
!               orbitalstore(i_block)%spinor_grads(i_ir)%spinor(2)%o      => &
!                    & orbitalstore(i_block)%spinor_grads(i_ir)%o(:,:,:,:,1,2)
!               orbitalstore(i_block)%spinor_grads(i_ir)%spinor(2)%o_imag => &
!                    & orbitalstore(i_block)%spinor_grads(i_ir)%o(:,:,:,:,2,2)
             end do
             if(first_block) spinor_grads => orbitalstore(i_block)%spinor_grads
          end if

       end do grid_blocks
    end if
    if( options_orbitals_on_file() ) then
       call error_handler("orbitalstore_init : ORBITLAS_ON_FILE not yet implemented")
    end if
    if(.not. orbitalstore_initialized) orbitalstore_initialized = .true.

  end subroutine orbitalstore_allocate
  !****************************************************************

  !****************************************************************
  subroutine orbitalstore_free()
    !  Purpose: deallocating memory or delete file for grid orbitals
    !           Called by orbitalstore_shutdown
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_block, i_ir, i_unique,alloc_stat
    integer(kind=i4_kind) :: mem
    !------------ Executable code --------------------------------

    if( options_orbitals_in_memory() ) then
       grid_blocks : do i_block = 1,N_grid_blocks
          DPRINT 'orbstore_free: block=',i_block

          if ( using(OS_ORBS_OB) ) then
             DPRINT 'orbstore_free: -orbs_ob'
             do i_ir = 1, symmetry_data_n_irreps()
                mem = size(orbitalstore(i_block)%orbs_ob(i_ir)%o)
                MEMLOG(-mem)

                counts(1,OS_ORBS_OB) = counts(1,OS_ORBS_OB) &
                     + size(orbitalstore(i_block)%orbs_ob(i_ir)%o)
                counts(2,OS_ORBS_OB) = counts(2,OS_ORBS_OB) &
                     + count(abs(orbitalstore(i_block)%orbs_ob(i_ir)%o) < cutoff )

                deallocate( orbitalstore(i_block)%orbs_ob(i_ir)%o, stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ("orbitalstore_free : deaallocation of orbs_ob%o failed")
             end do
             deallocate( orbitalstore(i_block)%orbs_ob, stat=alloc_stat )
             if ( alloc_stat /= 0) call error_handler &
                  ("orbitalstore_free : deaallocation of orbs_ob failed")
             stored_orbs_ob = .false.
          end if

          if ( using(OS_FCTS) ) then
             DPRINT 'orbstore_free: -fcts'
             mem = size(orbitalstore(i_block)%fcts%o)
             MEMLOG(-mem)
             deallocate(orbitalstore(i_block)%fcts%o, stat=alloc_stat)
             if(alloc_stat /= 0) call error_handler &
                  ("orbitalstore_free : deallocation of fcts%o failed")
          end if

          if ( using(OS_GRADS) ) then
             DPRINT 'orbstore_free: -grads'
             do i_ir = 1, symmetry_data_n_irreps()
                mem = size(orbitalstore(i_block)%grads(i_ir)%o)
                MEMLOG(-mem)
                deallocate( orbitalstore(i_block)%grads(i_ir)%o, stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ("orbitalstore_free : deallocation of grads%o failed")
             enddo
             deallocate(orbitalstore(i_block)%grads, stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_free: deallocation of grads failed")
             stored_grads = .false.
          end if

          if ( using(OS_FCTS_GRADS) ) then
             DPRINT 'orbstore_free: -fcts_grads'
             mem = size(orbitalstore(i_block)%fcts_grads%o)
             MEMLOG(-mem)
             deallocate(orbitalstore(i_block)%fcts_grads%o, stat=alloc_stat)
             if(alloc_stat /= 0) call error_handler &
                  ("orbitalstore_free : deallocation of fcts_grads%o failed")
          end if

          if ( using(OS_SEC_DER) ) then
             DPRINT 'orbstore_free: -sec_der'
             do i_ir = 1, symmetry_data_n_irreps()
                mem = size(orbitalstore(i_block)%sec_der(i_ir)%o)
                MEMLOG(-mem)
                deallocate( orbitalstore(i_block)%sec_der(i_ir)%o, stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ("orbitalstore_free :  deallocation of sec_der%o failed")
             end do
             deallocate( orbitalstore(i_block)%sec_der, stat=alloc_stat)
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_init :  deallocation of sec_der failed")
             stored_sec_der = .false.
          endif

          if ( using(OS_NUC_GRADS) ) then
             DPRINT 'orbstore_free: -nuc_grads'
             do i_ir = 1, symmetry_data_n_irreps()
                do i_unique=1,n_unique_atoms
                   mem = size(orbitalstore(i_block)%nuc_grads(i_unique,i_ir)%o)
                   MEMLOG(-mem)
                   deallocate( orbitalstore(i_block)% &
                               nuc_grads(i_unique,i_ir)%o ,stat=alloc_stat )
                   if ( alloc_stat /= 0 ) call error_handler &
                        ("orbitalstore_free: deallocation of nuc_grads%o failed")
                end do
             end do
             deallocate( orbitalstore(i_block)%nuc_grads,stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_init : allocate of nuc_grads failed")
             stored_nuc_grads = .false.
          endif

          if ( using(OS_NUC_SEC_DER) ) then
             DPRINT 'orbstore_free: -nuc_sec_der'
             do i_ir = 1, symmetry_data_n_irreps()
                do i_unique=1,n_unique_atoms
                   mem = size(orbitalstore(i_block)%nuc_sec_der(i_unique,i_ir)%o)
                   MEMLOG(-mem)
                   deallocate( orbitalstore(i_block)%                        &
                               nuc_sec_der(i_unique,i_ir)%o, stat=alloc_stat )
                   if ( alloc_stat /= 0 ) call error_handler &
                        ("orbitalstore_free : deallocation of nuc_sec_der%o failed")
                end do
             end do
             deallocate( orbitalstore(i_block)%nuc_sec_der, stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_free : deallocation of nuc_sec_der failed")
             stored_nuc_sec_der = .false.
          end if

          if ( using(OS_PHI_OB) ) then
             DPRINT 'orbstore_free: -phi_ob'
             do i_ir = 1, symmetry_data_n_irreps()
                mem = size(orbitalstore(i_block)%phi_ob(i_ir)%o)
                MEMLOG(-mem)
                deallocate( orbitalstore(i_block)%phi_ob(i_ir)%o, stat=alloc_stat )
                if ( alloc_stat /= 0 ) call error_handler &
                     ("orbitastore_free : deallocation of phi%o failed")
             enddo
             deallocate( orbitalstore(i_block)%phi_ob, stat=alloc_stat )
             if ( alloc_stat /= 0 ) call error_handler &
                  ("orbitalstore_free : deallocation of phi failed")
             stored_phi_ob = .false.
          end if

          if ( using(OS_ORBS_SPINOR_OB) ) then
             DPRINT 'orbstore_free: -orbs_spinor_ob'
             do i_ir = 1, symmetry_data_n_proj_irreps()
                mem = size(orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o)
                MEMLOG(-mem)
                deallocate( orbitalstore(i_block)%orbs_spinor_ob(i_ir)%o, stat=alloc_stat )
                ASSERT(alloc_stat==0)
             end do
             deallocate( orbitalstore(i_block)%orbs_spinor_ob, stat=alloc_stat )
             ASSERT(alloc_stat==0)
             stored_orbs_spinor_ob = .false.
          end if

          if ( using(OS_SPINOR_GRADS) ) then
             DPRINT 'orbstore_free: -spinor_grads'
             do i_ir = 1, symmetry_data_n_proj_irreps()
                mem = size(orbitalstore(i_block)%spinor_grads(i_ir)%o)
                MEMLOG(-mem)
                deallocate( orbitalstore(i_block)%spinor_grads(i_ir)%o, stat=alloc_stat )
                ASSERT(alloc_stat==0)
             end do
             deallocate( orbitalstore(i_block)%spinor_grads, stat=alloc_stat )
             ASSERT(alloc_stat==0)
             stored_spinor_grads = .false.
          end if

      end do grid_blocks
    end if

  end subroutine orbitalstore_free
  !****************************************************************

  !****************************************************************
  subroutine orbitalstore_rw(orbs_ob,grads,sec_der,       &
                             nuc_grads,nuc_sec_der,phi_ob,&
                             orbs_spinor_ob,spinor_grads, &
                             fcts, fcts_grads)
    !  Purpose: Save (rw_type="save") or read (rw_type="read")
    !           grid-orbitals
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(orbital_type),                  pointer, optional :: orbs_ob(:)        !
    type(orbital_type),                           optional :: fcts              !
    type(orbital_gradient_type),         pointer, optional :: grads(:)          !
    type(orbital_gradient_type),                  optional :: fcts_grads        !
    type(spinor_gradient_type),          pointer, optional :: spinor_grads(:)
    type(orbital_sec_der_type),          pointer, optional :: sec_der(:)        !
    type(orbital_nuclear_gradient_type), pointer, optional :: nuc_grads(:,:)    !
    type(orbital_nuclear_sec_der_type),  pointer, optional :: nuc_sec_der(:,:)  !
    type(orbital_spin_type),             pointer, optional :: phi_ob(:)         !
    type(spinor_type),                   pointer, optional :: orbs_spinor_ob(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------
    !        CLEAN ORBITALS BEFORE READING FROM FILE !!!!!!!

    if( present(orbs_ob)        )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_orbs_ob,grid_block_orbs_ob)
          orbs_ob  => orbitalstore(grid_block_orbs_ob)%orbs_ob
       end if
    end if

    if( present(fcts)           )  then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_fcts, grid_block_fcts)
          fcts = orbitalstore(grid_block_fcts)%fcts
       end if
    end if

    if( present(grads)          )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_grads,grid_block_grads)
          grads => orbitalstore(grid_block_grads)%grads
       end if
    end if

    if( present(fcts_grads)     )  then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_fcts_grads, grid_block_fcts_grads)
          fcts_grads = orbitalstore(grid_block_fcts_grads)%fcts_grads
       end if
    end if

    if( present(sec_der)        )   then
       if( options_orbitals_in_memory() ) then
          sec_der => orbitalstore(grid_block_sec_der)%sec_der
          call rw_set_counter(stored_sec_der,grid_block_sec_der)
       end if
    end if

    if( present(nuc_grads)      )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_nuc_grads,grid_block_nuc_grads)
          nuc_grads  => orbitalstore(grid_block_nuc_grads)%nuc_grads
       end if
    end if

    if( present(nuc_sec_der)    )   then
       if(options_orbitals_in_memory() ) then
          call rw_set_counter(stored_nuc_sec_der,grid_block_nuc_sec_der)
          nuc_sec_der => orbitalstore(grid_block_nuc_sec_der)%nuc_sec_der
       end if
    end if

    if( present(phi_ob)         )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_phi_ob,grid_block_phi_ob)
          phi_ob => orbitalstore(grid_block_phi_ob)%phi_ob
       end if
    end if

    if( present(orbs_spinor_ob) )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_orbs_spinor_ob,grid_block_orbs_spinor_ob)
          orbs_spinor_ob  => orbitalstore(grid_block_orbs_spinor_ob)%orbs_spinor_ob
       end if
    end if

    if( present(spinor_grads) )   then
       if( options_orbitals_in_memory() ) then
          call rw_set_counter(stored_spinor_grads,grid_block_spinor_grads)
          spinor_grads  => orbitalstore(grid_block_spinor_grads)%spinor_grads
       end if
    end if

    if(options_orbitals_on_file() ) then
       call error_handler("orbitalstore_free : ORBITLAS_ON_FILE not yet implemented")
    end if

  contains

    !****************************************************************
    subroutine rw_set_counter(stored_orbital,grid_block)
      !------------ Declaration of formal parameters ---------------
      logical,               intent(inout) :: stored_orbital
      integer(kind=i4_kind), intent(inout) :: grid_block
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      !------------ Executable code --------------------------------
      if(grid_block==N_grid_blocks) then
         if(.not.stored_orbital) stored_orbital = .true.
         grid_block = 0
      end if
      grid_block=grid_block+1
    end subroutine rw_set_counter
    !****************************************************************

  end subroutine orbitalstore_rw
  !****************************************************************

  !--------------- End of module ----------------------------------
end module orbitalstore_module

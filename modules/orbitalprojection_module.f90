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
module orbitalprojection_module
!---------------------------------------------------------------
!
!  Purpose: describes the projection of contracted and symmetry- 
!           adapted orbital fcts and charge and exchange
!           fitfunctions to the storage format used in scf-part.
!           This storage format only seperates Irreps and Partners.
!           All other Indizes are mapped to the Meta_Index running
!           over the dimension of the Irrep in the following way
!           (Index of outermost loop given first):
!
!           For orbital functions:
!             1. unique atom 
!             2. l
!             3. fct index of unique_atom_symadapt_type
!                "independent functions"
!             4. uncontracted exponents and contractions
!
!           For charge and exchange fitfunctions:
!             1. unique atom 
!             2. l (l=0, r**2, l=1, l=2, ...)
!             3. uncontracted exponents and contractions
!             4. fct index of unique_atom_symadapt_type
!                "independent functions"
!             Global contractions appended to the end:
!             1. unique atom
!             2. index of global contraction
!
!            The arrays orbitalprojection_?? defined here
!            allow to look up for a given Irrep (which partner
!            makes no difference) the first index of a given
!            unique atom and l
!            The arrays orbitalprojection_globcontr_?? defined 
!            here allow to look up the first index of the global
!            contractions of a given unique atom
!
!  Module called by: orbital_module
!
!  Author: TB (don`t blame me for the arbitrary order,
!              this is simply old lcgto)
!  Date: 11/95
!
!---------------------------------------------------------------------
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
!
! Modification
! Author: MM
! Date:   12/97
! Description: Introduction of orbitalprojection_spor_ob
!
!
! Modification put if(options_kinematic_factors) in some places
! Author: DG
! Date:  21/03/2001 
! Description: I need to pass matrix elements in uncontracted dimensions,
! I implemented orbitalprojection_spor_ob_k for it
!(where uncontracted dimensions are stored)
!ToDo: for the case of uncontructed dimensions put the switch if (options_kinematic_factors)
!into orbitalprojection_module.f90 and rid off the _k index everyware
!---------------------------------------------------------------------

#include "def.h"
use type_module ! type specification parameters
use options_module, only: options_spin_orbit, options_kinematic_factors ! spin orbit flag
implicit none
private         ! by default, all names are private
save
!== Interrupt end of public interface of module =================

!------------ Declaration of constants and variables ------------
integer, pointer, public :: orbitalprojection_ob(:,:,:) => NULL()
      ! orbital_index_ob(N_irreps,0:N_max_l,N_unique_atoms)
      ! first index of given l and ua for the given Irrep
integer, pointer, public :: orbitalprojection_spor_ob_k(:,:,:) => NULL() ! for the case of uncontructed dimension 

integer, pointer, public :: orbitalprojection_spor_ob(:,:,:) => NULL()
      ! orbital_index_ob(N_irreps,0:N_max_l,N_unique_atoms)
      ! first index of given l,spin (=j) and ua for the given Irrep

integer, pointer, public :: orbitalprojection_ch(:,:) => NULL()
integer, pointer, public :: orbitalprojection_cd(:,:) => NULL()
integer, pointer, public :: orbitalprojection_xc(:,:) => NULL()
integer, pointer, public :: orbitalprojection_ch_eperef(:,:) => NULL()
      ! orbital_index_??(-1:N_max_l,N_unique_atoms)
      ! first index of given l (-1 means r**2) and ua
      ! only totalsymmetric Irrep allowed
integer, pointer, public :: orbitalprojection_globcontr_ch(:) => NULL()
integer, pointer, public :: orbitalprojection_globcontr_xc(:) => NULL()
      ! orbital_index_globcontr_??(N_unique_atoms)

!------------ public functions and subroutines ------------------
public :: orbitalprojection_setup
public :: orbitalprojection_close
public :: orbitalprojection_print


  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Declaration of constants and variables ----
integer, private :: N_max_l_ob, N_max_l_ch, N_max_l_xc
      ! maximal number of l for all unique atoms



!---------------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

   subroutine orbitalprojection_setup()
    use comm, only: comm_rank
    use output_module, only: output_orbitalprojections
    use iounitadmin_module, only: output_unit
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call orbitalprojection_calculate()

        if ( output_orbitalprojections ) &
            call orbitalprojection_print(output_unit)
    endif

    call orbitalprojection_bcast()
   end subroutine orbitalprojection_setup


   !*************************************************************
   subroutine orbitalprojection_calculate()
     !  Purpose: calculates orbs_dim_indices and does allocation
     use symmetry_data_module
     use unique_atom_module 
     use operations_module, only:operations_core_density
     ! the contents of symmetry_data_module and unique_atom_module
     ! are used and must be calculated before !!!
     !** End of interface ***************************************
     implicit none
     !------------ Declaration of local variables ---------------
     integer :: i_sa_ob, i_sa_ch, i_sa_xc, i_ua, i_l, i_ir, i_bas,&
                i_sa_cd
#ifdef WITH_CORE_DENS
     type(unique_atom_atomic_dens_type), pointer :: uac
#endif
     logical                                     :: pseudopot
     type(unique_atom_type),             pointer :: ua
     type(unique_atom_basis_type),       pointer :: uab
     type(unique_atom_partner_type),     pointer :: uap
     !------------ Executable code -------------------------------

     call orbitalprojection_calcdims()
     call orbitalprojection_allocate()

     ! atomic orbitals
     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! loop over irreps
        do i_ir = 1, symmetry_data_n_proj_irreps()
           i_sa_ob = 1
           ! loop over unique atoms
           do i_ua = 1, N_unique_atoms
              ua => unique_atoms(i_ua)
              ! loop over l ob
              do i_l = 0, ua%lmax_ob
                 uab => ua%l_ob(i_l)
                 uap => ua%symadapt_spor_partner(i_ir,i_l)
                 orbitalprojection_spor_ob(i_ir,i_l,i_ua) = i_sa_ob
                 ! add space for uncontracted exponents and contractions of all independent fcts
                 i_sa_ob = i_sa_ob + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
                      * uap%N_independent_fcts
           enddo! loop over l ob
           enddo! loop over unique atoms
        enddo! loop over irreps
!!$     else ! options_spin_orbit
     endif
        if (options_kinematic_factors .and. options_spin_orbit) then
           !
        ! SPIN ORBIT
        !
        ! loop over irreps
        do i_ir = 1, symmetry_data_n_proj_irreps()
           i_sa_ob = 1
           ! loop over unique atoms
           do i_ua = 1, N_unique_atoms
              ua => unique_atoms(i_ua)
              ! loop over l ob
              do i_l = 0, ua%lmax_ob
                 uab => ua%l_ob(i_l)
                 uap => ua%symadapt_spor_partner(i_ir,i_l)
                 orbitalprojection_spor_ob_k(i_ir,i_l,i_ua) = i_sa_ob
                 ! add space for uncontracted exponents  of all independent fcts
                 i_sa_ob = i_sa_ob + (uab%N_exponents ) &
                      * uap%N_independent_fcts
                 
           enddo! loop over l ob
           enddo! loop over unique atoms
        enddo! loop over irreps
!!$!
!!$        ! SPIN ORBIT
!!$        !
!!$        ! loop over irreps
!!$        do i_ir = 1, symmetry_data_n_proj_irreps()
!!$           i_sa_ob = 1
!!$           ! loop over unique atoms
!!$           do i_ua = 1, N_unique_atoms
!!$              ua => unique_atoms(i_ua)
!!$              ! loop over l ob
!!$              do i_l = 0, ua%lmax_ob
!!$                 uab => ua%l_ob(i_l)
!!$                 do i_spin = 1,2
!!$                    if ((i_l.eq.0).and.(i_spin.eq.1)) then
!!$                       cycle
!!$                    endif
!!$                    uap => ua%symadapt_spor_partner(i_ir,i_l)
!!$                    orbitalprojection_spor_ob_k(i_ir,i_l,i_ua) = i_sa_ob
!!$                    ! add space for uncontracted exponents and contractions of all independent fcts
!!$                    i_sa_ob = i_sa_ob + (uab%N_exponents) &
!!$                         * uap%N_independent_fcts
!!$                 enddo ! loop over spin
!!$              enddo! loop over l ob
!!$           enddo! loop over unique atoms
!!$        enddo! loop over irreps
        end if
        !
        ! STANDARD SCF (NO SPIN ORBIT)
        !
        ! loop over irreps
        do i_ir = 1, symmetry_data_n_irreps()
           i_sa_ob = 1
           ! loop over unique atoms
           do i_ua = 1, N_unique_atoms
              ua => unique_atoms(i_ua)
              ! loop over l ob
              do i_l = 0, ua%lmax_ob
                 uab => ua%l_ob(i_l)
                 uap => ua%symadapt_partner(i_ir,i_l)
                 orbitalprojection_ob(i_ir,i_l,i_ua) = i_sa_ob
                 ! add space for uncontracted exponents and contractions of all independent fcts
                 i_sa_ob = i_sa_ob + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
                      * uap%N_independent_fcts
              enddo! loop over l ob
           enddo! loop over unique atoms
        enddo! loop over irreps
!!$     endif ! options_spin_orbit

     ! fitfunctions ch and xc and core
     i_ir = get_totalsymmetric_irrep()
     i_sa_ch = 1
     i_sa_xc = 1
     i_sa_cd = 1
     ! loop over unique atoms
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        pseudopot = &
        pseudopot_present.and.ua%zc /= 0.0_r8_kind .and. &
             .not.operations_core_density
        ! loop over bas ch
        do i_bas = -1, ua%lmax_ch
           if (i_bas .eq. -1) then
              uab => ua%l_ch(0)
              uap => ua%symadapt_partner(i_ir,0)
           elseif (i_bas .eq. 0) then
              uab => ua%r2_ch
              uap => ua%symadapt_partner(i_ir,0)
           else
              i_l = i_bas
              uab => ua%l_ch(i_l)
              uap => ua%symadapt_partner(i_ir,i_l)
           endif
           orbitalprojection_ch(i_bas,i_ua) = i_sa_ch
           ! add space for uncontracted exponents and contractions of all independent fcts
           i_sa_ch = i_sa_ch + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
                               * uap%N_independent_fcts
        enddo! loop over l ch
#ifdef WITH_CORE_DENS
        ! loop over basis functions for core density
        if ( pseudopot.and.core_density_setup ) then
          do i_bas = -1, 0
             if (i_bas .eq. -1) then
                uac => ua%s_core
             elseif (i_bas .eq. 0) then
                uac => ua%r2_core
             endif
             orbitalprojection_cd(i_bas,i_ua) = i_sa_cd
             ! add space for uncontracted fitting functions 
             ! uap%N_independent_fcts = 1 for s- and r^2-type fit functions 
             i_sa_cd = i_sa_cd + uac%N_exponents
          enddo ! loop over core functions with s_core and r2_core
        
        endif ! end of core orbitals
#endif
        ! loop over bas xc
        do i_bas = -1, ua%lmax_xc
           if (i_bas .eq. -1) then
              uab => ua%l_xc(0)
              uap => ua%symadapt_partner(i_ir,0)
           elseif (i_bas .eq. 0) then
              uab => ua%r2_xc
              uap => ua%symadapt_partner(i_ir,0)
           else
              i_l = i_bas
              uab => ua%l_xc(i_l)
              uap => ua%symadapt_partner(i_ir,i_l)
           endif
           orbitalprojection_xc(i_bas,i_ua) = i_sa_xc
           ! add space for uncontracted exponents and contractions of all independent fcts
           i_sa_xc = i_sa_xc + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
                               * uap%N_independent_fcts
        enddo! loop over l xc
     enddo! loop over unique atoms
     ! global contractions
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        orbitalprojection_globcontr_ch(i_ua) = i_sa_ch
        orbitalprojection_globcontr_xc(i_ua) = i_sa_xc
        i_sa_ch = i_sa_ch + ua%N_glob_cons_ch
        i_sa_xc = i_sa_xc + ua%N_glob_cons_xc
     enddo! loop over unique atoms
     
   end subroutine orbitalprojection_calculate
   !*************************************************************



   !*************************************************************
   subroutine orbitalprojection_print(io_unit)
     !  Purpose: prints orbs_dim_indices
     use symmetry_data_module
     use unique_atom_module
     ! the contents of symmetry_data_module and unique_atom_module
     ! are used and must be calculated before !!!
     implicit none
     integer, intent(in) :: io_unit
     !** End of interface ***************************************
     !------------ Declaration of local variables ---------------
     integer :: i_ua, i_l, i_ir, i_bas
     type(unique_atom_type),             pointer :: ua
     !------------ Executable code -------------------------------

     write (io_unit,*)
     write (io_unit,*) "Orbitalprojections Orbitals"
     write (io_unit,*)

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! loop over irreps
        do i_ir = 1, symmetry_data_n_proj_irreps()
           write (io_unit,*) "Irrep ", i_ir
           ! loop over unique atoms
           do i_ua = 1, N_unique_atoms
              ua => unique_atoms(i_ua)
              write (io_unit,'("unique atom ",I3,"  l = 1,2, ... L: ", 20(I4,","))') &
                   i_ua, (orbitalprojection_spor_ob(i_ir,i_l,i_ua), i_l = 1, ua%lmax_ob)
           enddo! loop over unique atoms
        enddo! loop over irreps
     else ! options_spin_orbit
        !
        ! STANDARD SCF (NO SPIN ORBIT)
        !
        ! loop over irreps
        do i_ir = 1, symmetry_data_n_irreps()
           write (io_unit,*) "Irrep ", i_ir
           ! loop over unique atoms
           do i_ua = 1, N_unique_atoms
              ua => unique_atoms(i_ua)
              write (io_unit,'("unique atom ",I3,"  l = 0,1,2, ... : ", 20(I4,","))') &
                   i_ua, (orbitalprojection_ob(i_ir,i_l,i_ua), i_l = 0, ua%lmax_ob)
           enddo! loop over unique atoms
        enddo! loop over irreps
     endif ! options_spin_orbit


     write (io_unit,*)
     write (io_unit,*) "Orbitalprojections Fitfunctions Charge "
     write (io_unit,*)

     ! fitfunctions ch and xc
     ! loop over unique atoms
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        write (io_unit,'("unique atom ",I3,"  l=0, r**2, l=1,2, ... : ", 20(I4,","))') &
             i_ua, (orbitalprojection_ch(i_bas,i_ua), i_bas = -1, ua%lmax_ch)
        if (ua%zc/=0.0_r8_kind) then
        write (io_unit,*)
        write (io_unit,*) "Orbitalprojections Fitfunctions Core "
        write (io_unit,*)
        write (io_unit,'("unique atom ",I3,"  l=0, r**2,       ", 20(I4,","))') &
             i_ua, (orbitalprojection_cd(i_bas,i_ua), i_bas = -1, 0)
        endif
     enddo! loop over unique atoms
     ! global contractions
     write (io_unit,*) "Global Projections Start" 
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        write (io_unit,'("unique atom ",I3," : ", I4)') &
             i_ua, orbitalprojection_globcontr_ch(i_ua)
     enddo! loop over unique atoms
     

     write (io_unit,*)
     write (io_unit,*) "Orbitalprojections Fitfunctions Exchange"
     write (io_unit,*)

     ! fitfunctions ch and xc
     ! loop over unique atoms
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        write (io_unit,'("unique atom ",I3,"  l=0, r**2, l=1,2, ... : ", 20(I4,","))') &
             i_ua, (orbitalprojection_xc(i_bas,i_ua), i_bas = -1, ua%lmax_xc)
     enddo! loop over unique atoms
     ! global contractions
     write (io_unit,*) "Global Projections Start" 
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        write (io_unit,'("unique atom ",I3," : ", I4)') &
             i_ua, orbitalprojection_globcontr_xc(i_ua)
     enddo! loop over unique atoms

     write (io_unit,*)

   end subroutine orbitalprojection_print
   !*************************************************************


   !*************************************************************
   subroutine orbitalprojection_allocate
     !  Purpose: allocate orbital_index_??
     !** End of interface ***************************************
     use unique_atom_module 
     use symmetry_data_module
     implicit none
     !------------ Declaration of local variables ---------------
     integer              :: status, ts
     logical:: a_switch
     !------------ Executable code ------------------------------
     ts = get_totalsymmetric_irrep()
     if (options_spin_orbit) then

        allocate( orbitalprojection_spor_ob( symmetry_data_n_proj_irreps(), &
             0:N_max_l_ob, &
             N_unique_atoms             ), &
             stat=status )
        orbitalprojection_spor_ob = -1
!!$        if (options_kinematic_factors) then
!!$           allocate( orbitalprojection_spor_ob_k( symmetry_data_n_proj_irreps(), &
!!$             0:N_max_l_ob,2, &
!!$             N_unique_atoms             ), &
!!$             stat=status )
!!$        end if
        if (options_kinematic_factors .and. options_spin_orbit) then
           allocate( orbitalprojection_spor_ob_k( symmetry_data_n_proj_irreps(), &
             0:N_max_l_ob, &
             N_unique_atoms             ), &
             stat=status )
        end if
        
!!$     else
     endif
!      if(.not.associated(orbitalprojection_ob)) then
        a_switch=.not.associated(orbitalprojection_ob)
        ASSERT(a_switch)
        allocate( orbitalprojection_ob( symmetry_data_n_irreps(), &
             0:N_max_l_ob, &
             N_unique_atoms             ), &
             stat=status )
        ASSERT(status.eq.0)
        orbitalprojection_ob = -1
!       endif
!!$     endif

     ts = get_totalsymmetric_irrep()
!     if(.not.associated(orbitalprojection_ch)) then
      a_switch=.not.associated(orbitalprojection_ch)
      ASSERT(a_switch)
      allocate( orbitalprojection_ch( -1:N_max_l_ch,N_unique_atoms ),stat=status )
      ASSERT(status.eq.0)
      orbitalprojection_ch = -1
!     endif

     if(pseudopot_present.and.core_density_setup)then
     allocate( orbitalprojection_cd(-1:0, N_unique_atoms ),&
               stat = status )
     if ( status .ne. 0 ) call error_handler( &
          "orbitalprojection_allocate: allocate of orbital_core_indices_ch failed")
     endif

!     if(.not.associated(orbitalprojection_xc)) then
      a_switch=.not.associated(orbitalprojection_xc)
      ASSERT(a_switch)
      allocate( orbitalprojection_xc( -1:N_max_l_xc, N_unique_atoms ), stat=status )
      ASSERT(status.eq.0)
      orbitalprojection_xc = -1
!     endif

!    if(.not.associated(orbitalprojection_globcontr_ch)) then
     a_switch=.not.associated(orbitalprojection_globcontr_ch)
     ASSERT(a_switch)
     allocate( orbitalprojection_globcontr_ch(N_unique_atoms), stat=status )
     ASSERT(status.eq.0)
     orbitalprojection_globcontr_ch = -1
!    endif

!    if(.not.associated(orbitalprojection_globcontr_xc)) then
     a_switch=.not.associated(orbitalprojection_globcontr_xc)
     ASSERT(a_switch)
     allocate( orbitalprojection_globcontr_xc(N_unique_atoms), stat=status )
     ASSERT(status.eq.0)
     orbitalprojection_globcontr_xc = -1
!    endif

   end subroutine orbitalprojection_allocate
   !*************************************************************


   !*************************************************************
   subroutine orbitalprojection_close()
     !  Purpose: deallocate orbital_index_?? 
     ! Subroutine called by: main_master/main_slave
     !** End of interface ***************************************
     use unique_atom_module, only: pseudopot_present,&
                                   core_density_setup
     !------------ Declaration of local variables ---------------
     integer              ::  alloc_stat

     if (options_spin_orbit) then
        deallocate(orbitalprojection_spor_ob,stat=alloc_stat)
        if (options_kinematic_factors) deallocate(orbitalprojection_spor_ob_k,stat=alloc_stat)
        if (alloc_stat/=0) call error_handler &
             ("orbitalprojection_close : allocation (1) failed")
!!$     else
     endif

        if(associated(orbitalprojection_ob)) then
         deallocate(orbitalprojection_ob,stat=alloc_stat)
         ASSERT(alloc_stat.eq.0)
        endif

!!$     endif
     if(associated(orbitalprojection_ch) ) then
      deallocate(orbitalprojection_ch,stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
     endif

     if ( pseudopot_present.and.core_density_setup)then
     deallocate(orbitalprojection_cd, stat=alloc_stat)
     if (alloc_stat/=0) call error_handler &
          ("orbitalprojection_close : allocation (2c) failed")
     endif

    if(associated(orbitalprojection_xc)) then
     deallocate(orbitalprojection_xc,stat=alloc_stat)
     ASSERT(alloc_stat.eq.0)
    endif

    
    if(associated(orbitalprojection_globcontr_ch)) then
     deallocate(orbitalprojection_globcontr_ch,stat=alloc_stat)
     ASSERT(alloc_stat.eq.0)
    endif

    if(associated(orbitalprojection_globcontr_xc)) then
     deallocate(orbitalprojection_globcontr_xc,stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
    endif

   end subroutine orbitalprojection_close
   !*************************************************************


   !*************************************************************
   subroutine orbitalprojection_calcdims
     ! Purpose: determine N_max_l_??
     use unique_atom_module 
     ! desription of unique atoms, must be calculated before
     !** End of interface ***************************************
     implicit none
     integer(kind=i4_kind)           :: i_ua
     type(unique_atom_type), pointer :: ua
     !------------ Executable code ------------------------------
     N_max_l_ob = 0
     N_max_l_ch = 0
     N_max_l_xc = 0
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        N_max_l_ob = max( N_max_l_ob, ua%lmax_ob )
        N_max_l_xc = max( N_max_l_xc, ua%lmax_xc )
        N_max_l_ch = max( N_max_l_ch, ua%lmax_ch )
     enddo
   end subroutine orbitalprojection_calcdims
   !*************************************************************


   !****************************************************************************
   subroutine orbitalprojection_bcast
     ! Purpose: bcasts data
     !** End of interface ******************************************************
     use unique_atom_module, only: pseudopot_present                           &
                                 , core_density_setup
     use comm,               only: comm_bcast                                  &
                                 , comm_rank
     !------------ Executable code ------------------------------
     !
     call comm_bcast(N_max_l_ob)
     call comm_bcast(N_max_l_ch)
     call comm_bcast(N_max_l_xc)
     !
     if ( comm_rank() /= 0 ) then
       call orbitalprojection_allocate()
     endif
     !
     if (options_spin_orbit) then
       call comm_bcast(orbitalprojection_spor_ob)
     endif
     !
     call comm_bcast(orbitalprojection_ob)
     !
     call comm_bcast(orbitalprojection_ch)
     call comm_bcast(orbitalprojection_xc)
     call comm_bcast(orbitalprojection_globcontr_ch)
     call comm_bcast(orbitalprojection_globcontr_xc)
     !
     if (pseudopot_present.and.core_density_setup) then
       call comm_bcast(orbitalprojection_cd)
     endif
     !
     if (options_kinematic_factors .and. options_spin_orbit) then
       call comm_bcast(orbitalprojection_spor_ob_k)
     endif
     !
   end subroutine orbitalprojection_bcast
   !****************************************************************************

!--------------- End of module ----------------------------------
end module orbitalprojection_module

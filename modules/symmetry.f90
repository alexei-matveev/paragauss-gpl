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
module symmetry

# include "def.h"
  use type_module, only: i4_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private

public :: main_symm

! and these are arguments it accepts:
integer(i4_kind), parameter, public :: MAIN_SYMM_INIT = 1
integer(i4_kind), parameter, public :: MAIN_SYMM_DONE = 2


contains

subroutine  main_symm (do_what)
!
!  Purpose: main routine of the symmetry part
!
!
!  Subroutine called by: main / main_master
!
!
!  References: publisher document: symdoc
!
!
!  Author: MM
!  Date: 12/96
!
!================================================================
! End of public interface of module
!================================================================
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

  !------------ Modules used ------------------------------------
  use iounitadmin_module, only: output_unit, write_to_output_units
  use output_module, only: output_main_symm, output_symadapt
  use time_module, only: start_timer, stop_timer
  use timer_module, only: timer_symm
  use group_module, only: group_name, group_num_cl, group_num_ir, &
        group_num_re, irrep_can, proj_irrep_can, group_symmel_alloc, &
        group_symmel_read, group_symmel_multab, group_class_gen, &
        group_char_gen, group_proj_char_gen, group_generator_gen, &
        group_irrep_label, group_proj_irrep_pre_label, &
        group_symmel_dealloc, subgroup_chain_gen, group_close
  use efm_module, only: efm_irrep_gen, efm_proj_irrep_gen, &
       efm_proj_irrep_label, efm_cg_gen, efm_proj_cg_gen, efm_done
  use symm_module, only: symm_get_pcoupling, symm_get_cccoupling, &
       symm_init, symm_done, symm_symmadapt, symm_dipole_selrules_gen
  use symm_positions, only: symm_equivalents_gen, &
       symm_positions_init => init, symm_positions_done => done
  use symmetry_data_module, only: symmetry_data_point_group, &
       symmetry_data_set, symmetry_data_set_pcoupling, &
       symmetry_data_set_cccoupling
  use unique_atom_methods, only: unique_atom_symadapt_write, &
       unique_atom_symequiv_write
  use operations_module, only: operations_dipole, &
       operations_properties, operations_gtensor
  use efield_module, only: efield_applied
  use options_module, only: options_spin_orbit
  implicit none
  integer(i4_kind), intent(in), optional :: do_what
  ! *** end of interface ***



  !------------ Declaration of local variables --------------------
  !integer(i4_kind), parameter :: MainSymmInit = 1, MainSymmDone = 2
  integer(i4_kind)            :: do_what_
  integer(kind=i4_kind)       :: i
  ! counter
!!$  logical                               :: group_only
  ! only process group part of program


!----------------------------------------------------------------
!------------ Executable code -----------------------------------

  DPRINT 'main_symm: entered'
#ifdef FPP_DEBUG
  do i=1,size(ptgrps)
     DPRINT i,'=>'//ptgrps(i)//'<'
  enddo
#endif
  do_what_ = MAIN_SYMM_INIT
  if( present(do_what) ) do_what_ = do_what

  if (output_main_symm) call write_to_output_units("main_symm: entered")

  if(do_what_.eq.MAIN_SYMM_INIT)then
     call init_ms()
     if (output_symadapt) then
        call unique_atom_symadapt_write(output_unit)
        call unique_atom_symequiv_write(output_unit)
     endif
  endif

  if(do_what_.eq.MAIN_SYMM_DONE)then
     call done_ms()
  endif

  if (output_main_symm) call write_to_output_units("main_symm: exited")

contains

  subroutine init_ms()
    implicit none
    ! internal sub of main_symm
    ! **** end of interface ***

    DPRINT 'main_symm::init: entered'

  call start_timer(timer_symm)

  !
  ! determine group properties (class structure,
  ! characters,irreps)
  !

  ! set pointgroup

  ! determine symmetry group from input
  group_name = symmetry_data_point_group()

!!$  group_only = .false.
!!$
!!$  if (group_only) then
!!$     call error_handler("main_symm: I DO enter this branch!")
!!$     ! allocate symmetry elements
!!$     call group_symmel_alloc
!!$
!!$     ! load symmetry elements
!!$     call group_symmel_read
!!$
!!$     ! generate group multiplication table
!!$     call group_symmel_multab
!!$
!!$     ! generate classes
!!$     call group_class_gen
!!$
!!$     ! generates characters
!!$     call group_char_gen
!!$
!!$     ! generates characters of projective representations
!!$     call group_proj_char_gen
!!$
!!$     ! generate generator
!!$     call group_generator_gen
!!$
!!$     ! label irreps
!!$     call group_irrep_label
!!$
!!$     ! label projective irreps
!!$     call group_proj_irrep_label
!!$
!!$     ! generate subgroup chain
!!$     call subgroup_chain_gen
!!$
!!$     ! deallocate symmetry elements
!!$     call group_symmel_dealloc
!!$
!!$  else
     ! -------------------------------------------------------------------------------------------
     ! Determine Group Properties and Character Table
     ! -------------------------------------------------------------------------------------------

     ! allocate symmetry elements
     call group_symmel_alloc()

     ! load symmetry elements
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating symmetry elements")
     call group_symmel_read()

     ! generate group multiplication table
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating multiplication table")
     call group_symmel_multab()

     ! generate classes
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating class structure")
     call group_class_gen()

     ! generates characters
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating character table")
     call group_char_gen(group_num_cl)

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! generates characters of projective representations
        call group_proj_char_gen(group_num_re)
     endif

     ! generate generator
     call group_generator_gen()


     ! -------------------------------------------------------------------------------------------
     ! generate rotation matrices
     ! -------------------------------------------------------------------------------------------

     ! generate and allocate rotation matrices
     DPRINT 'main_symm::init: call symm_init(...)'
     call symm_init(spin_orbit=options_spin_orbit)


     ! -------------------------------------------------------------------------------------------
     ! label irreps
     ! -------------------------------------------------------------------------------------------

     ! label irreps
     if (output_main_symm) call write_to_output_units( &
          "main_symm: labeling irreps")
     call group_irrep_label

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! generates characters of projective representations
        call group_proj_irrep_pre_label
     endif


     ! generate subgroup chain
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating subgroup chains")
     call subgroup_chain_gen

     ! -------------------------------------------------------------------------------------------
     ! investigate unique atoms
     ! -------------------------------------------------------------------------------------------


     ! generate symmetry equivalent atoms and distances
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating symmetryequivalent atoms and distances")
     DPRINT 'main_symm::init: call symm_positions_init()'
     call symm_positions_init()
     DPRINT 'main_symm::init: call symm_equivalents_gen()'
     call symm_equivalents_gen()
     ! -------------------------------------------------------------------------------------------
     ! determine irrep representation matrices
     ! -------------------------------------------------------------------------------------------

     ! generate irreps
     if (output_main_symm) call write_to_output_units( &
          "main_symm: generating irreducible representation matrices")
     call efm_irrep_gen

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! generate projective irreps
        call write_to_output_units( "main_symm: generating projective irreducible representation matrices")
        call efm_proj_irrep_gen
        ! now do relabeling of irreps
        call efm_proj_irrep_label
     endif
     ! test rotation matrix generator
     !call efm_test_ylm

     ! -------------------------------------------------------------------------------------------
     ! set symmetry data of cluster code/1
     ! -------------------------------------------------------------------------------------------

     ! allocate space
     call symmetry_data_set(n_irrep=group_num_ir,point_group=group_name)

     ! fill symmetry information
     do i=1,group_num_ir
        call symmetry_data_set(index=i,&
             &partner=irrep_can(i)%time_dimension,name=irrep_can(i)%label,&
             pseudo=irrep_can(i)%pseudo)
     end do

     call symmetry_data_set(SPIN_ORBIT=options_spin_orbit )

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! allocate space for projective irreps
        call symmetry_data_set(n_proj_irrep=group_num_re,proj=.true.)

        ! fill symmetry information for projective irreps
        do i=1,group_num_re
           call symmetry_data_set(index=i,&
                &partner_proj=proj_irrep_can(i)%dimension,name_proj=proj_irrep_can(i)%label,&
                &jz=proj_irrep_can(i)%sign_of_jz)
        end do
     endif

     ! -------------------------------------------------------------------------------------------
     ! calculate cg coefficients
     ! -------------------------------------------------------------------------------------------

     if (output_main_symm) call write_to_output_units( &
          "main_symm: calculating cg coefficients")
     call efm_cg_gen

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! calculate cg coefficients of projective irreps
        call write_to_output_units( "main_symm: calculating projective cg coefficients")
        call efm_proj_cg_gen
     endif

     ! -------------------------------------------------------------------------------------------
     ! determine dipole selection rules
     ! -------------------------------------------------------------------------------------------

     if ( operations_dipole .or. efield_applied() .or. &
          operations_properties ) then
        ! calculate dipole selection rules
        if (output_main_symm) call write_to_output_units( &
             "main_symm: calculating dipole selection rules")
        if(.not.operations_gtensor) call symm_dipole_selrules_gen
     endif

     ! -------------------------------------------------------------------------------------------
     ! do symmetry adaption
     ! -------------------------------------------------------------------------------------------

     if (output_main_symm) call write_to_output_units( &
          "main_symm: doing symmetry adaption " )
     call symm_symmadapt(spor=options_spin_orbit)
     ! -------------------------------------------------------------------------------------------
     ! set symmetry data of cluster code/2
     ! -------------------------------------------------------------------------------------------

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        call symmetry_data_set_pcoupling(symm_get_pcoupling())
        call symmetry_data_set_cccoupling(symm_get_cccoupling())
     endif

     call stop_timer(timer_symm)

   end subroutine init_ms

   subroutine done_ms()
     implicit none
     ! internal of mayn_symm
     ! *** end of interface ***

     DPRINT 'main_symm::done: entered'

     ! -----------------------------------------------
     ! deallocate data
     ! -----------------------------------------------

     ! deallocate symmetry elements
     DPRINT 'main_symm::done: call group_symmel_dealloc()'
     call group_symmel_dealloc()

     ! deallocate rotation matrices
     DPRINT 'main_symm::done: call symm_done(...)'
     call symm_done(spin_orbit=options_spin_orbit)

     DPRINT 'main_symm::done: call symm_positions_done()'
     call symm_positions_done()

     DPRINT 'main_symm::done: call efm_done(...)'
     call efm_done(spin_orbit=options_spin_orbit)

     ! deallocate all remaining pointers
     DPRINT 'main_symm::done: call group_close(...)'
     call group_close()

   end subroutine done_ms

end subroutine main_symm

end module symmetry

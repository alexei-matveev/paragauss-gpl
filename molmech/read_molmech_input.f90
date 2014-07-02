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
subroutine read_molmech_input()

  !===================================================================
  ! End of public interface of subroutine
  !===================================================================
  !------------ Modules used --- -----------------------------------
  use type_module
  use inp_out_module
  use common_data_module
  use tasks_main_options_module
  use slab_module
  use species_module
  use external_field_module
  use potentials_module
  use n_body_lists_module
  use coulomb_module
  use hess_and_opt_module
  use pc_array_module
  use qmmm1_interface_module
  use cavity_module, only: read_solv, write_solv

  implicit none

  logical :: read_nm
  real(r8_kind) :: qdm(4)
  !------------ Executable code -----------------------------------

  call get_input_device(qmmm)
  call get_output_device(qmmm)

  call write_banner()

  read_nm=read_topic()
  if(read_nm) call write_topic_to_output()

  read_nm=read_tasks()
  call write_tasks_to_output()

  if(calc_energy) then
     read_nm=read_main_options()
     call write_options_to_output()

     if(solvent) then
        read_nm=read_solv()
        if(solvent) call write_solv()
     end if

     read_nm=read_nb_options()
     call write_nb_options()
  end if

  read_nm=read_cell()
  if(.not. read_nm) read_nm=read_vectors()
  if(read_nm) then 
     call write_cell_vect_to_output()
     lattice_calc=.true.
     slab_calc=.false.
  else
     if(calc_ppc_array) &
          call error_handler("MolMech: Generation of PC array is possible only for 3-D case")
     read_nm=read_slab_cell()
     if(.not. read_nm) read_nm=read_slab_vectors()
     if(read_nm) then
        call write_slab_cell_vect_to_output()
        lattice_calc=.false.
        slab_calc=.true.
     else
        lattice_calc=.false.
        minimal_image=.true.
        slab_calc=.false.
     end if
  end if

  if((lattice_calc .or. slab_calc) .and. with_optimizer) &
       call error_handler("MolMech: Only cluster can be optimized with PG Optimizer")

  if(lattice_calc .and. minimal_image) then
!!$     if(valid_cutoff(rcut_vdw) /= zero) then 
!!$        write(message,'(f13.8)') valid_cutoff(rcut_vdw)
!!$        call error_handler &
!!$             ("MolMech: Van der Waals cutoff radius exceeds half"// &
!!$             " the shortest current allowed value "//trim(message)// &
!!$             " of used unit cell")
!!$     end if
  else if(slab_calc .and. minimal_image) then
!!$     if(valid_cutoff_slab(rcut_vdw) /= zero) then 
!!$        write(message,'(f13.8)') valid_cutoff(rcut_vdw)
!!$        call error_handler &
!!$             ("MolMech: Van der Waals cutoff radius exceeds half"// &
!!$             " the shortest current allowed value "//trim(message)// &
!!$             " of used slab unit cell")
!!$     end if
  end if

  if(calc_energy) then
     read_nm=read_ext_field_type()
     if(read_nm) then
        if (ext_field) then 
           if (lattice_calc .or. slab_calc) then 
              call error_handler &
                   ("MolMech: External field currently can be apply to nonperiodic systems")
           else
              call write_ext_field_type_to_output()
           end if
        end if
     end if
  end if

  read_nm=read_coord()
  if(.not. read_nm) call error_handler &
       ("MolMech: Where are your coordinates?")

  if(qmmm == 2 .or. qmmm == 3) then
     read_nm=read_qmmm()
  end if
  if(qmmm == -1) then
     read_nm=read_qmmm1()
  end if

  if(read_gx) then
     read_nm=read_coord_gx()
     if(.not. read_nm) call error_handler &
          ("MolMech: Something wrong with GX file that You have tried to read in")
  end if

  call read_species()
  if(with_optimizer) call read_pdc()
  call write_species_to_output()

  if(calc_ppc_array) then
     read_nm=read_pc_array_options()
     read_nm=read_embed_cluster()
  end if

  if(calc_energy) then
     if(ext_field) call read_ext_field()

     qdm=check_total_charge_and_dipmom()
     if(lattice_calc .or. slab_calc) then
        if(abs(qdm(1)) >= small) call error_handler &
             ("MolMech: The total charge of system is not equal ZERO")
     end if

     write(output_device,'(a17,/)') ' Potentials used:'

     call potential_list_init()
     read_nm=read_potential()
!!$  if(.not. read_nm) call error_handler &
!!$       ("MolMech: Where are your potentials?")
     if(trim(energy_unit) /= "KJ/MOL" .or. trim(length_unit) /= "ANGSTROM") then
        call convert_ff_parameters()
     end if
     call write_potential_to_output()

     if(automatic_nb_lists) then
        call autobuilding_nb_lists()
     else
        nb_lists_from_input=read_Nb_lists()
     end if

     if(coulomb) then
        read_nm=read_coulomb_options()
        call write_coulomb_to_output()
        if(trim(coulomb_type) == "DIPOLE-DIPOLE") then
           if(ext_field) call error_handler &
                ("MolMech: Extenal PC: currently no possibility"// &
                " to calculate PC - bond dipoles interaction")
           read_nm=read_dipoles()
           if(read_nm) call write_dipoles_to_output()
        end if
     end if

     if(calc_optimization .or. with_optimizer) then
        read_nm=read_opt_options()
        if(.not.with_optimizer) call write_opt_options_to_output()
     end if
     if(calc_optimization .or. calc_hessian) then
        read_nm=read_hess_options()
        call write_hess_options_to_output()
     end if
  end if
     
  call close_input_device()

contains

  !******************************************************************
  subroutine write_banner()

    write(output_device,'(20x,a40)') '****************************************'
    write(output_device,'(20x,a40)') '*                                      *'
    write(output_device,'(20x,a40)') '*        ParaGauss (vers. 3.0)         *'
    write(output_device,'(20x,a40)') '* Molecular Mechanics Part (vers. 0.95)*'
    write(output_device,'(20x,a40)') '*                                      *'
    write(output_device,'(20x,a40)') '****************************************'
    write(output_device,'(a1)') ' '

  end subroutine write_banner
  !******************************************************************

end subroutine read_molmech_input




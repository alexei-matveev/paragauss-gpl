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
module  efield_module
  !-------------------------------------------------------------------
  !
  !  Purpose: This module calculates integrals of
  !           constant external electric field
  !           from dipol integrals:
  !
  !  Int = sum[i=x,y,z] ( E  * Int       )
  !     E                  i      Dipol,i
  !
  !  Those integrals are either stored integralstore_module
  !  or written to tape efield.dat. They are used in scf
  !  part by ham_calc_module.
  !  The required dipol integrals are calculated as well
  !  but deallocated afterwards.
  !
  !  Routines to read / write the input namelist EFIELD
  !  and to print applied electrical field are also
  !  included (input unit: a.u.).
  !  All variables are private, put there are functions to
  !  inquire if and which electrical field is switched on.
  !
  !  To investigate singlet triplet excitations, the
  !  (non-physical) option to apply the same field with
  !  reverse orientation to both spins in an open shell
  !  calculation is added.
  !
  !  Please note that applied electrical fields may breake the
  !  symmetry of the problem. Th symmetry part should check this.
  !
  !  Module called by: main_master (before scf part)
  !  (all data in this module only exist on master!)
  !
  !  Author: TB
  !  Date:   11/97
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ public functions and subroutines ---------------------
  public efield_read, efield_write, efield_calculate_integrals, &
       efield_applied, efield_field, efield_print, &
       efield_reverse_for_spins, efield_gradient, efield_intensity, &
       efield_change, efield_der, efield_gradient_store
  public :: efield_send_recv
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------
  real(kind=r8_kind), dimension(3) :: field = 0.0_r8_kind
  logical, parameter :: &
       df_apply = .false., &
       df_reverse_for_spins  = .false., &
       df_intensity  = .false.
  logical :: apply, reverse_for_spins, intensity

  namelist /efield/ apply, reverse_for_spins, intensity


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  !*************************************************************
  subroutine efield_read()
    !  Purpose: reads namelist efield from input
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use input_module
    use operations_module, only: operations_dipole
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: status, unit, i
    !------------ Executable code ------------------------------------
    apply = df_apply
    reverse_for_spins = df_reverse_for_spins
    intensity = df_intensity
    if ( input_line_is_namelist("EFIELD") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=efield, iostat=status)
       if (status .gt. 0) call input_error( &
            "EFIELD_READ: namelist EFIELD")
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit,*,iostat=status) (field(i), i=1,3)
       if (status .gt. 0) call input_error( &
            "EFIELD_READ: electrical field vector")
     endif

     if(apply.and..not.operations_dipole)then
       call input_error("efield_read: enable dipole moments to use electric field!")
       ! otherwise in integral_calc_quad_dipole the relevant integral
       ! calculation routines are not called at all. This results in
       ! junk & NaN numbers.
     endif
     if(intensity.and..not.operations_dipole)then
       call input_error("efield_read: enable dipole moments to use electric field!")
     endif
  end subroutine efield_read
  !*************************************************************


  !*************************************************************
  subroutine efield_write(iounit)
    !
    ! Purpose: writes input namelist efield.
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full, echo
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    !** End of interface *****************************************

    integer :: status

    call start("EFIELD","EFIELD_WRITE", &
         iounit,operations_echo_input_level)
    call flag("APPLY             ", apply            , df_apply             )
    call flag("REVERSE_FOR_SPINS ", reverse_for_spins, df_reverse_for_spins )
    call flag("INTENSITY         ", intensity        , df_intensity         )
    call stop(empty_line=.false.)
    if (echo()) then
       write (iounit,fmt='(3F15.6)',iostat=status) field
       if (status /= 0) call error_handler( &
            "EFIELD_WRITE: write of field failed")
       write (iounit,'()')
    endif
  end subroutine efield_write
  !*************************************************************


  !*************************************************************
  subroutine efield_print(iounit)
    !  Purpose: writes information about applied efield to iounit
    !------------ Modules used ----------------------------------
    use options_module, only:  options_spin_restricted
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer, intent(in) :: iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind), parameter :: au2Vm = 5.14224382e+11
    integer :: status
    !------------ Executable code ------------------------------------
    write (iounit,err=999,fmt=*)
    if (apply) then
       write (iounit,iostat=status, &
            fmt='("The following homogenous external electrical field is applied:")')
       write (iounit,err=999,fmt='(16X,"x",14X,"y",14X,"z")')
       write (iounit,err=999,fmt='(" (a.u.)  ",3F15.6)') field
       write (iounit,err=999,fmt='(" (V/m)   ",3E15.6)') field * au2Vm
       if ( reverse_for_spins .and. .not. options_spin_restricted() ) then
          write (iounit,err=999, &
               fmt='("The field will be applied in opposite direction for secound spin")')
       endif
    else
       write (iounit,err=999, &
            fmt='("No external electrical field is applied")')
    endif
    write (iounit,err=999,fmt=*)
    return
999 call error_handler("EFIELD_PRINT: write error")
  end subroutine efield_print
  !*************************************************************


  !*************************************************************
  logical function efield_applied()
    !  Purpose: returns if external electrical field is applied
    !** End of interface ***************************************
    efield_applied = apply
  end function efield_applied
  !*************************************************************


  !*************************************************************
  logical function efield_reverse_for_spins()
    !  Purpose: returns if external electrical field is will be
    !  applied in reverse orientation for spin up and down.
    !  This (non-physical) option allows to investigate
    !  singlet-triplet excitations
    !** End of interface ***************************************
    efield_reverse_for_spins = reverse_for_spins
  end function efield_reverse_for_spins
  !*************************************************************

  !*************************************************************
  logical function efield_intensity()
    !  Purpose: switch for IR intensity calculation by finite difference
    !           method of forces in presence of electric field
    !** End of interface ***************************************
    efield_intensity = intensity
  end function efield_intensity
  !*************************************************************

  !*************************************************************
  subroutine  hessian_update
    !  Purpose: calculation of the mode vector for calculation of intensity 
    !           without presence of  electric field
    !** End of interface ***************************************
    use options_module, only: update_hessian_iteration
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------
    if (efield_intensity()) then
      update_hessian_iteration = 1
    endif
  end subroutine  hessian_update
  !*************************************************************

  !*************************************************************
  subroutine efield_change(a)
    !
    ! Passing the electric field strength in six direction
    !
    !** End of interface *****************************************
    implicit none
    integer(kind=i4_kind), intent(in)      :: a
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind)  :: mag
    !------------ Executable code ------------------------------------
    mag =  sqrt(dot_product(field,field))
    ! magnitude of field strength

    select case (a)
      case (1)
         field(1) =  mag
         field(2) =  0.0
         field(3) =  0.0
      case (2)
         field(1) = - mag
         field(2) =  0.0
         field(3) =  0.0
      case (3)
         field(1) =  0.0
         field(2) =  mag
         field(3) =  0.0
      case (4)
         field(1) =  0.0
         field(2) = - mag
         field(3) =  0.0
      case (5)
         field(1) =  0.0
         field(2) =  0.0
         field(3) =  mag
      case (6)
         field(1) =  0.0
         field(2) =  0.0
         field(3) = - mag
      case (7)
         apply    = .false. ! here field strength is zero as apply is false
         call hessian_update
      case default
        print*, "error in input"
        ABORT('no such loop is there')
    end select
  end subroutine efield_change
  !*************************************************************

  !*************************************************************
  subroutine efield_der(gradder,error)
    !
    ! This subroutine is for making finite difference of gradient in
    ! presence of electric field with respect to field strenth for IR intensity
    !
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use unique_atom_module
    use iounitadmin_module
    implicit none
    logical, intent(out)               :: error
    real(kind=r8_kind), intent(out)    :: gradder(:,:,:)
    ! --- Declaration of local variables ----------------
    real(kind=r8_kind)    :: gradder_dis(3,6,N_unique_atoms), de ! field strength
    integer(kind=i4_kind) :: io_grad,i,j
    logical               :: exist
    ! --- Executable code -------------------------------

    inquire (EXIST = exist,file="gradient.dat")
    if (.not.exist) then
      print*,'efield_module: efield_gradient: gradient_der is not present'
      error = .true.
      RETURN
    endif
    error = .false.
    io_grad = openget_iounit(status='old',form='formatted',&
            file="gradient.dat")
    gradder_dis = 0.0
    ! now we read the gradient for all atoms and calculate the
    ! derivatives with respect to electric field
    do i = 1, 6
      do j = 1, N_unique_atoms
        read(io_grad,*) gradder_dis(:,i,j)
      end do
    end do
    call returnclose_iounit(io_grad)

    ! magnitude of field strength for finite difference
    de = sqrt(dot_product(field,field))

    gradder = 0.0
    ! taking the finite difference of gradient with respect to field
    gradder(:,1,:) = (gradder_dis(:,1,:)-gradder_dis(:,2,:))/(2*de)
    gradder(:,2,:) = (gradder_dis(:,3,:)-gradder_dis(:,4,:))/(2*de)
    gradder(:,3,:) = (gradder_dis(:,5,:)-gradder_dis(:,6,:))/(2*de)
    print*, "the derivative of gradient is"
    write (*,'(3F20.12)')  gradder
  end subroutine efield_der
  !*************************************************************

  !*************************************************************
  function efield_field()
    !  Purpose: returns applied external electrical field
    real(kind=r8_kind), dimension(3) :: efield_field
    !** End of interface ***************************************
    efield_field = field
  end function efield_field
  !*************************************************************


  !*************************************************************
  subroutine efield_calculate_integrals ()
    !
    ! Purpose: calculates integrals of external electrical field.  for
    ! this purpose, dipole integral  part is run, but dipole integrals
    ! are deallocated afterwards again The integrals are either stored
    ! integralstore_module  or written to  tape efield.dat.   They are
    ! used in scf part by ham_calc_module.
    !
    ! Runs on all workers.
    !
    use integralstore_module, only: integralstore_2cob_efield                  &
                                  , integralstore_allocate_efield              &
                                  , integralstore_deallocate_efield
    use dipole_module
    use output_module, only: output_efield_calc
    use options_module, only: options_integrals_on_file
    use operations_module, only: operations_dipole
    use integralpar_module, only: integralpar_set
    use timer_module, only: timer_efield, timer_dipole_integral
    use time_module
    use iounitadmin_module, only: output_unit, write_to_output_units
    use interfaces, only: main_integral
    implicit none
    !** End of interface *****************************************

    if (.not. apply) call error_handler( &
         "efield_calculate_integrals: no external electrival field applied")

    if(.not.operations_dipole)then
      call error_handler("efield_calculate_integrals: enable dipole moments to use electric field!")
    endif

    call start_timer(timer_efield)

    call efield_print(output_unit)

    if ( .not. allocated(dipole_integrals) ) then
       ! allocate storage for dipole integrals
       if (output_efield_calc) call write_to_output_units( &
            "efield_calculate_integrals: allocate storage for dipole integrals")
       call dipole_allocate(offdiagonals=.false.)
       ! run dipol integral part
       if (output_efield_calc) call write_to_output_units( &
            "efield_calculate_integrals: run dipol integral part")
       call start_timer(timer_dipole_integral)
       call integralpar_set('Dipole')

       call main_integral ()
       call stop_timer(timer_dipole_integral)
       call print_timer(timer_dipole_integral,output_unit, &
         "Calculation of dipole integrals", &
         diff=.true.,absolut=.false.)
    else
      WARN('efield_calculate_integrals: not calculating ints')
    endif

    ! allocate efield integrals
    if (output_efield_calc) call write_to_output_units( &
         "efield_calculate_integrals: allocate efield integrals")
    call integralstore_allocate_efield( integralsize() )

    ! calculate efield integrals from dipol integrals
    if (output_efield_calc) call write_to_output_units( &
         "efield_calculate_integrals: calculate efield integrals")
    call calculate()

    ! deallocate dipol integrals
    if (output_efield_calc) call write_to_output_units( &
         "efield_calculate_integrals: deallocate dipol integrals")
    call dipole_free()

    ! write efield integrals to tape and deallocate them
    if ( options_integrals_on_file() ) then
       if (output_efield_calc) call write_to_output_units( &
            "efield_calculate_integrals: write efield integrals to tape")
       call write_integrals()
       if (output_efield_calc) call write_to_output_units( &
            "efield_calculate_integrals: deallocate efield integrals")
       call integralstore_deallocate_efield()
    endif

    call stop_timer(timer_efield)
    call print_timer(timer_efield,output_unit, &
         "Calculation of integrals of external electrical field (including dipole integrals)", &
         diff=.true.,absolut=.false.)

  contains


    integer(kind=i4_kind) function integralsize()
      ! returns total size of efield integrals in triangular storage
      use symmetry_data_module
      implicit none
      !------------ Declaration of local variables -----------------
      integer(kind=i4_kind) :: i_ir, dim_ir
      !------------ Executable code --------------------------------
      integralsize = 0
      do i_ir = 1, symmetry_data_n_irreps()
         dim_ir = symmetry_data_dimension(i_ir)
         integralsize = integralsize + dim_ir*(dim_ir+1)/2
      enddo
    end function integralsize


    subroutine write_integrals()
      ! write efield integrals to tape
      use readwriteblocked_module
      use filename_module, only: tmpfile
      implicit none
      !------------ Declaration of local variables -----------------
      type(readwriteblocked_tapehandle) :: th
      !------------ Executable code --------------------------------
      call readwriteblocked_startwrite(trim(tmpfile("efield.dat")), th)
      call readwriteblocked_write(integralstore_2cob_efield,th)
      call readwriteblocked_stopwrite(th)
    end subroutine write_integrals


    subroutine calculate()
      ! calculate efield integrals from dipol integrals
      use symmetry_data_module
      implicit none
      !------------ Declaration of local variables -----------------
      integer(kind=i4_kind) :: i_ir,dim_ir,i_xyz,i_ip,i_meta,dim_diag
      !------------ Executable code --------------------------------
      integralstore_2cob_efield = 0.0_r8_kind
      i_meta = 1
      i_ip = 1
      irrep: do i_ir = 1, symmetry_data_n_irreps()
         dim_ir = symmetry_data_dimension(i_ir)
         dim_diag = dim_ir*(dim_ir+1)/2
         do i_xyz = 1, 3
            integralstore_2cob_efield(i_meta:i_meta+dim_diag-1) = &
                 integralstore_2cob_efield(i_meta:i_meta+dim_diag-1) + &
                 dipole_integrals(i_xyz,i_ip,i_ip)%diagonal * field(i_xyz)
         enddo
         i_ip = i_ip + symmetry_data_n_partners(i_ir)
         i_meta = i_meta + dim_diag
      enddo irrep  
    end subroutine calculate


  end subroutine efield_calculate_integrals
  !*************************************************************
  subroutine efield_gradient(grad_final)
    !
    ! contribution of nuclear force tothe total force in presence of
    ! electric field
    !
    use unique_atom_module
    use iounitadmin_module
    use datatype, only: arrmat2
    implicit none
    type(arrmat2), intent(inout)      :: grad_final(:)
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)            :: i_unique,i_equal
    real(kind=r8_kind)               :: nuc_efield(3)
    !------------ Executable code ------------------------------------
    do i_unique = 1, N_unique_atoms
      nuc_efield(:) = - efield_field() * unique_atoms (i_unique)%z
    !
    ! Now I have confirm the sign of the nuclear contribution to the total gradient
    ! in presence of electric field.
    !
      write(output_unit,*) 'Unique Center:',i_unique
        do i_equal=1,unique_atoms(i_unique)%n_equal_atoms
          print *, "contribution from nucleus to the total force for unique", & 
                                                 i_unique, i_equal, nuc_efield
          write(output_unit,'(A20,3F20.12)')  &
                  'Electronic Contri.:', grad_final(i_unique)%m(:,i_equal)

          write(output_unit,'(A20,3F20.12)')  &
                  'Nucleus Contri.:', nuc_efield

          grad_final(i_unique)%m(:,i_equal) = grad_final(i_unique)%m(:,i_equal) + nuc_efield

          write(output_unit,'(A20,3F20.12)') 'Total Contri.:',&
                           grad_final(i_unique)%m(:,i_equal)
        enddo
    enddo
  end subroutine efield_gradient
  !*************************************************************

  !*************************************************************
  subroutine efield_gradient_store(grad_final)
    !
    ! storing the gradient in presence of electric field for IR intensities
    !
    use unique_atom_module
    use iounitadmin_module
    use datatype, only: arrmat2
    implicit none
    type(arrmat2), intent(in)      :: grad_final(:)
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)            :: io_grad
    integer(kind=i4_kind)            :: i_unique,i_equal
    !------------ Executable code ------------------------------------
    io_grad=openget_iounit(status='unknown',form='formatted',position='append',&
       file="gradient.dat")
    do i_unique = 1, N_unique_atoms
      do i_equal=1, unique_atoms(i_unique)%n_equal_atoms
        write(io_grad,'(3F20.12)') grad_final(i_unique)%m(:,i_equal)
      enddo
    enddo
    call returnclose_iounit(io_grad)
  end subroutine efield_gradient_store
  !*************************************************************
 
  !*************************************************************
  subroutine efield_send_recv()
    !
    ! Passing some of variable from this module to the slaves for parallel jobs
    !
    use comm_module
    use msgtag_module,     only: msgtag_packed_message
    use xpack,             only: pck, upck
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------
   
    if ( comm_i_am_master() ) then
      call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
      call pck(apply)
      call pck(field)
      call comm_send()
    else
      call comm_save_recv(comm_master_host, msgtag_packed_message)
      call upck(apply)
      call upck(field)
    endif
  end subroutine efield_send_recv
  !*************************************************************

  !--------------- End of module -------------------------------------
end module efield_module

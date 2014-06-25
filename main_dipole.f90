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
subroutine  main_dipole
  !
  !  Main  routine  of  dipole   part.   Runs  dipole  integral  part,
  !  calculates dipole  moments and prints  them.  Runs in  a parallel
  !  context.   Called  by:   main_master().   Important  calls:  call
  !  dipoleg_calculate().
  !
  !  References: Ph.D. Thesis A Goerling
  !
  !  Author: TB
  !  Date: 8/97
  !
  !
  !----------------------------------------------------------------
  !================================================================
  ! End of public interface of module
  !================================================================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Author: HH
  ! Date:   10/97
  ! Description: For response module add following:
  ! If operations_response=T then call subroutine
  ! "dipole_trans_response" from dipol_module.f90.
  ! This subroutine writes all transition dip.mom.
  ! in an temporary, unformatted file in temp dir.
  !
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
# include "def.h"
  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use comm, only: comm_rank
  use dipole_module, only: dipole_integrals, &
       dipole_offdiagonals_calculated, dipole_xes_spectra, dipole_allocate, &
       dipole_calculate, dipole_print, dipole_free, dipole_make_xes_spectra, &
       dipole_trans_response, dipole_transitionmoment_uf, &
       dipole_transitionmoment_f, dipole_write_simol
#ifdef WITH_GTENSOR
  use hfc_module
  use gtensor_module
#endif
  use output_module
  use integralpar_module, only: integralpar_set
  use iounitadmin_module
  use timer_module
  use time_module
  use filename_module
  use operations_module, only: operations_response, operations_dipole, &
       operations_gtensor, operations_hfc
#ifdef WITH_RESPONSE
  use response_module, only: response_do_dipole
#endif
  use options_module, only: options_spin_orbit
  use interfaces, only: main_integral

  implicit none

  !------------ Declaration of local variables --------------------
  logical :: offdiagonals_required, new_integrals_required,  &
       old_integrals_exist
  integer :: iounit

  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------


  hfcso: if (operations_hfc) then
#ifdef WITH_GTENSOR
     call start_timer (timer_dipole)

     DPRINT 'main_dipole : hfc part'

     call start_timer (timer_dipole_integral)
     call say ("new_integrals_required")

     ! allocate storage for integrals
     call say ("allocate storage for integrals")
     call hfc_allocate ()


     ! run integral part
     call say ("integralpar_set_dipole")
     call integralpar_set ('DipoleOff')

     call say ("main_integral")
     call main_integral ()

     call stop_timer (timer_dipole_integral)


     ! calculate gtensors
     call start_timer (timer_dipole_calculate)


     call say ("calculate hfc")
     call hfc_calculate ()

     call stop_timer (timer_dipole_calculate)

     call hfc_free ()

     call stop_timer (timer_dipole)
  elseif (operations_hfc .and.  .not. options_spin_orbit) then
     ABORT ('dead branch')
#else
     ABORT ('recompile -DWITH_GTENSOR')
#endif
  end if hfcso

  opgt: if (operations_gtensor) then
#ifdef WITH_GTENSOR
     call start_timer (timer_dipole)

     DPRINT 'main_dipole : gtensor part'

     call start_timer (timer_dipole_integral)
     call say ("main_gten: new_integrals_required")

     ! allocate storage for integrals
     call say ("main_gten: allocate storage for integrals")
    ! call dipoleg_allocate ()
      call gtensor_allocate ()

     ! run integral part
     call say ("main_gten: integralpar_set_dipole")
     call integralpar_set ('DipoleOff')

     call say ("main_gten: main_integral")
     call main_integral ()

     call stop_timer (timer_dipole_integral)


     ! calculate gtensors
     call start_timer (timer_dipole_calculate)


     call say ("calculate g-tensor")
     call gtensor_calculate ()

     call stop_timer (timer_dipole_calculate)

    ! call dipoleg_free ()
     call gtensor_free ()

     call stop_timer (timer_dipole)
#else
     ABORT ('recompile -DWITH_GTENSOR')
#endif
  end if opgt


  opdip: if (operations_dipole) then

  call start_timer (timer_dipole)

  ! determine what needs to be done
  offdiagonals_required = &
       output_dipole_transitionm_f .or. output_dipole_transitionm_uf &
       & .or. operations_response.or. dipole_xes_spectra ()
  old_integrals_exist = allocated (dipole_integrals)
  if (old_integrals_exist) then
     new_integrals_required = offdiagonals_required .and.  &
          .not. dipole_offdiagonals_calculated
  else
     new_integrals_required = .true.
  endif


  if (new_integrals_required) then
     call start_timer (timer_dipole_integral)
     call say ("new_integrals_required")

     ! deallocate old integrals before new ones are calculated
     if (old_integrals_exist) then
        call say ("deallocate old integrals")
        call dipole_free ()
     endif

     ! allocate storage for integrals
     call say ("allocate storage for integrals")
     call dipole_allocate (offdiagonals_required)

     ! run integral part
     call say ("integralpar_set_dipole")
     if (offdiagonals_required) then
       call integralpar_set ('DipoleOff')
     else
       call integralpar_set ('Dipole')
     endif
     call say ("main_integral")
     call main_integral ()

     call stop_timer (timer_dipole_integral)
  endif


  ! calculate dipole moments
  call start_timer (timer_dipole_calculate)
  call say ("calculate dipole moment")
  call dipole_calculate ()
  call stop_timer (timer_dipole_calculate)


  call start_timer (timer_dipole_print)

  !
  ! Print dipole  moments. The struct with  the data to  be printed is
  ! only valid on master. See dipole_calculate().
  !
  if (comm_rank() == 0) then
     call say ("call dipole_print()")
     call dipole_print (output_unit, output_dipole_detailed)
     call say ("done dipole_print()")
  endif


  ! print total dipol moment to a special file used for simol frequency calculations
  if (output_dipole_simol) then
     call say ("write simol file")
     call dipole_write_simol ()
  endif


  ! calculate transition dipole moments and print them formatted
  if (output_dipole_transitionm_f) then
     call say ("calculate transition dipole moments and print them formatted")
     call dipole_transitionmoment_f (output_unit)
  endif


  ! calculate transition dipole moments and write them unformatted
  if (output_dipole_transitionm_uf) then
     call say ("calculate transition dipole moments and write them unformatted")
     iounit = openget_iounit (trim (outfile ("dipoletransitionmoments")), &
          form="UNFORMATTED", status="REPLACE", action="WRITE")
     call dipole_transitionmoment_uf (iounit)
     call returnclose_iounit (iounit, status="KEEP")
  endif

  !  calculate transition dipol moments and write them unformatted to
  ! tmp-file for response module
  if (operations_response) then
#ifdef WITH_RESPONSE
     if (response_do_dipole ()) then
     call say ("calculate transition dipole moments for response module")
     iounit = openget_iounit (trim (tmpfile ("resp_dipoles_tmp.dat")), &
          form="UNFORMATTED", status="REPLACE", action="WRITE")
     call dipole_trans_response (iounit)
     call returnclose_iounit (iounit, status="KEEP")
     endif
#else
     ABORT ('recompile -DWITH_RESPONSE')
#endif
  endif

  ! calculate xes_spectra
  if (dipole_xes_spectra ()) then
     call say ("dipole_make_xes_spectra")
     iounit = openget_iounit (trim (outfile ("xes.dat")), &
          form="FORMATTED", status="REPLACE", action="WRITE")
     call dipole_make_xes_spectra (iounit)
     call returnclose_iounit (iounit, status="KEEP")
  end if
  call stop_timer (timer_dipole_print)

  ! deallocate integrals and calculated dipole moments
  call say ("deallocate integrals and calculated dipole moments")
  call dipole_free ()

  call stop_timer (timer_dipole)


end if opdip

contains

  subroutine say (phrase)
    use output_module, only: output_main_dipole
    use iounitadmin_module, only: write_to_output_units
    implicit none
    character (len=*), intent (in) :: phrase
    ! *** end of interface ***

    if (output_main_dipole) then
        call write_to_output_units ("main_dipole: " // phrase)
    endif
  end subroutine say
end subroutine main_dipole

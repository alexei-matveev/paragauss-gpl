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
module fermi_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains routines for the level broadening
  !
  !
  !  Module called by: main_scf
  !
  !  References: Old LCGTO and references therein
  !
  !  Author: FN
  !  Date: 4/96
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification
  !
  ! Author: Thomas Seemueller
  ! Date:   11/99
  ! Description: newton and newton_spin: added step restriction into
  !              the newton procedures; redesigned the way of deter-
  !              mining the Newton step delta_e
  !              get_cutoff: case(fermi_broadening): changed n to
  !              (1 - 1e-12) and (1e-12)
  !              added: get_spin_diff_parameters
  !
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined datatypes
  use symmetry_data_module, only: ssym, &
       symmetry_data_n_spin,symmetry_data_n_irr,symmetry_data_dimension,&
       symmetry_data_n_partners
  use occupation_module
  use iounitadmin_module, only: write_to_output_units,output_unit,&
       stdout_unit
  use eigen_data_module, only: eigval,holes_per_irrep,spin_of_hole,&
       n_holes_per_irrep
  use output_module, only: output_fermi,output_fermi_newton
  use options_module, only: options_spin_orbit

  use operations_module, only: operations_echo_input_level, operations_core_density

  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module =================

  !------------ Declaration of public variables -------------------
  public :: fermi_level_broad, &
            fermi_fix_up_and_down ! kept for compatibility

  !------------ public functions and subroutines ------------------
  public fermi_reoccup, fermi_read, fermi_write_input, fermi_get_entropy, &
       fermi_read_scfcontrol, check_occ_fermi


!================================================================
! End of public interface of module
!================================================================

  ! --------- private interface ----------------
  interface e_fermi_start
     module procedure ef_start
     module procedure ef_start_spin
  end interface

  !------------ Declaration of constants and variables ----
  ! these constants are mainly used in the erfc-function
  real(kind=r8_kind), parameter    :: zero = 0.0_r8_kind,       &
       half=0.5_r8_kind,                                        &
       one=1.0_r8_kind, two=2.0_r8_kind, three=3.0_r8_kind,     &
       four=4.0_r8_kind, five=5.0_r8_kind,                      &
       sqpi2 = 2.506628274631001_r8_kind,                       &
       sqrt2 = 1.414213562373095_r8_kind,                       &
       pi    = 3.141592653589793_r8_kind,                       &
       g_fac = 2.563102894212887_r8_kind ! depends on n1 and n2 !

  real(kind=r8_kind), parameter :: sigma_tol = 4.0e-5_r8_kind
  real(kind=r8_kind), parameter :: convert1 = 27.211652_r8_kind

  !---- private parameters for enumeration type LEVEL_BROADENING
  integer(kind=i4_kind), parameter :: gauss_broadening = 1, &
                                      fermi_broadening = 2, &
                                      sinus_broadening = 3
  ! ----------- Default values for input parameters -------------
  logical :: df_fermi_level_broad = .false. , & ! kept for compatibility
             df_fermi_gauss       = .false. , &
             df_fermi_fermi       = .false. , &
             df_fermi_sinus       = .false. , &
             df_fermi_up_down     = .false. , & ! kept for compatibility
             df_fermi_common_ef   =  .true.     ! kept for compatibility
  integer(kind=i4_kind) :: df_fermi_unpaired   =   0, & ! kept for compatibility
                           df_fermi_newton_max = 100
  real(kind=r8_kind) :: df_fermi_e_range  =     0.0_r8_kind, & ! eV
                        df_fermi_sigma    =     0.0_r8_kind, & ! eV
                        df_fermi_cutoff   =     5.0_r8_kind, & ! Ht
                        df_fermi_accuracy = 1.0E-10_r8_kind

  ! ----------- fermi input parameters ---------------------------
  logical :: fermi_level_broad    , & ! kept for compatibility (in input)
             fermi_gauss          , &
             fermi_fermi          , &
             fermi_sinus          , &
             fermi_fix_up_and_down, & ! kept for compatibility (in input)
             fermi_common_ef          ! kept for compatibility (in input)
  integer(kind=i4_kind) :: fermi_unpaired  , & ! kept for compatibility (in input)
                           fermi_newton_max
  real(kind=r8_kind) :: fermi_energy_range, &
                        fermi_sigma       , &
                        fermi_cutoff      , &
                        fermi_accuracy

  namelist /fermi/ fermi_level_broad    , & ! kept for compatibility (in input)
                   fermi_gauss          , &
                   fermi_fermi          , &
                   fermi_sinus          , &
                   fermi_energy_range   , &
                   fermi_sigma          , &
                   fermi_cutoff         , &
                   fermi_fix_up_and_down, & ! kept for compatibility (in input)
                   fermi_common_ef      , & ! kept for compatibility (in input)
                   fermi_unpaired       , & ! kept for compatibility (in input)
                   fermi_newton_max     , &
                   fermi_accuracy

  ! internal enumeration type to control the level_broadening procedure
  integer(kind=i4_kind) :: level_broadening
  ! entropy contribution, /= zero only for 'fermi' and 'sinus'
  real(kind=r8_kind)    :: entropy = 0.0_r8_kind
  ! these are useful in several routines
  integer(i4_kind) :: num_spin  ! FIXME: does it need to be global?
  integer(kind=i4_kind) :: num_orb,num_orb_spin

  ! dimensions of irreps
  ! (in order to account for SPIN ORBIT more easily)
  integer(kind=i4_kind)                :: n_irrep = -1
  integer(kind=i4_kind),allocatable    :: dim(:)
  integer(kind=i4_kind),allocatable    :: partner(:)
  ! n_irrep    : number of irreps
  ! dim_irrep : number of independent functions in irrep


  !----------------------------------------------------------------
contains


  !*************************************************************

  subroutine fermi_reoccup()
    !  Purpose: main routine for the reoccupation including
    !           gaussian level broadening around the
    !           fermi energy.
    !** End of interface *****************************************
    !
    !  Author: FN
    !  Date: 4/96
    !
    !------------ Modules used ----------------------------------
    use iounitadmin_module, only: output_unit
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)             :: efs = zero
    real(kind=r8_kind),allocatable :: efs_spin(:),e_fermi_spin(:)
    integer(kind=i4_kind)          :: alloc_stat,is
    real(kind=r8_kind)             :: e_fermi
    real(kind=r8_kind),parameter   :: convert=27.211652_r8_kind
    logical                        :: fixed_spin

    external error_handler

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim(n_irrep),partner(n_irrep))
       dim = ssym%dim_proj
       partner = ssym%partner_proj
    else
       n_irrep = ssym%n_irrep
       allocate(dim(n_irrep),partner(n_irrep))
       dim  = ssym%dim
       partner = ssym%partner
    endif

    select case (level_broadening)
    case (gauss_broadening)
       call write_to_output_units("FERMI_REOCUP : &
            &GAUSSIAN level broadening with sigma = ",re=fermi_sigma)
       call write_to_output_units("               &
            &                    and energy_range = ",re=fermi_energy_range)
    case(fermi_broadening)
       call write_to_output_units &
            ("FERMI_REOCUP : FERMI function ")
    case(sinus_broadening)
       call write_to_output_units &
            ("FERMI_REOCUP : SINUS function ")
    end select

    ! make the switch fixed_spin_diff (occupation_module) available
    call get_fix_spin_switch(fixed_spin)

    if (fixed_spin) then
       ! if we have a fixed spin configuration

       allocate (efs_spin(ssym%n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fermi_reoccup : allocation 1 failed")
       efs_spin=zero
       allocate(e_fermi_spin(ssym%n_spin),STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("fermi_reoccup: alloction 2 failed")
       e_fermi_spin=zero

       ! get a start value for e_fermi
       ! FIXME: maybe (HOMO+LUMO)/2 is a better start?
       call e_fermi_start(efs_spin)

       ! find the 'real' e_fermi
       call newton_spin(efs_spin,e_fermi_spin)

       if(output_fermi) then
          do is=1,ssym%n_spin
             write(output_unit,1100)is,convert*efs_spin(is),&
                  convert*e_fermi_spin(is)
          enddo
          write(output_unit,*)" "
       endif

       deallocate(efs_spin,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("fermi_reoccup : deallocation (1) failed")
       deallocate(e_fermi_spin,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("fermi_reoccup : deallocation (2) failed")
    else

       ! we let the electron configuration relax to
       ! whatever it wants to
       ! FIXME: maybe (HOMO+LUMO)/2 is a better start?
       call e_fermi_start(efs)

       call newton(efs,e_fermi)

       if(output_fermi) then
          write(output_unit,1000)convert*efs,convert*e_fermi
          write(output_unit,*)" "
       endif

    endif
1000 FORMAT (' highest occ. level=',F13.6,&
          '      fermi energy=',F13.6)
1100 FORMAT (' Spin ',I2,'  highest occ. level=',F13.6,&
          '      fermi energy=',F13.6)

    entropy = entropy_contribution()

    deallocate(dim,partner,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler( &
         "occupation_level_sort: deallocation of dim_irrep failed")
    n_irrep = -1

  end subroutine fermi_reoccup

  !*************************************************************

  subroutine ef_start(efs)
    !  Purpose:  produce the start fermi energy 'efs' for
    !            the Newton-Rhaphson procedure.
    !            This start fermi energy should be
    !            the highest occupied level (counting
    !            from the Aufbauprinciple) in order to reach
    !            convergence in 4-5 steps ( see P.Knappe,
    !            Dissertation, 1990)
    !** End of interface *****************************************
    !
    !------------ Modules used ----------------------------------
!    use occupation_module, only: occupation_level_sort, &
!         get_n_elec,num_list

    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),    intent(out) :: efs
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)               :: n_elec,remmel,degen,delta_hole,&
                                        occ_hole,remmel_core
    integer(kind=i4_kind)            :: n_elec_int,step_dummy
    integer(kind=i4_kind)            :: i,m,is,index
    logical                          :: hole

    external error_handler
    !------------ Executable code --------------------------------

    ASSERT(n_irrep/=-1)

    ! this fills the variable 'num_list', where
    ! num_list(i)%m(m,is) contains the number of the
    ! eigenvalue eigval(i)%m(m,is)
    call occupation_level_sort(eigval)

    !---- now perform the Aufbau-Principle without actually
    !     filling the electrons into the levels

    call get_n_elec(n_elec)
    n_elec_int = nint(n_elec,kind=i4_kind)

    num_orb=0
    do i=1,n_irrep                            ! this is the
       num_orb=num_orb+dim(i)*partner(i) ! total number
    enddo                                          ! of orbitals

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       num_spin = 1             ! number of el.
    else
       num_orb = num_orb * 2 ! 2 electrons fit into each orbital
       num_spin = 3 - ssym%n_spin ! number of el.
    endif

    if(n_elec_int.gt.num_orb) call error_handler&
         ( "ef_start :more electrons than orbitals" )
    step_dummy=1

    remmel_core = 0.0
    if(operations_core_density) call get_n_core_elec(remmel_core)
    remmel = n_elec                              ! remmel =
    orbital : do index=1,num_orb                 ! REMaining
       !                                         ! MELectrons
       do i=1,n_irrep
          do is=1,ssym%n_spin
             do m=1,dim(i)
                if(num_list(i)%m(m,is).eq.index) then
                   degen = num_spin * partner(i)
                   if (operations_core_density .and. remmel_core > 0.0_r8_kind) then
                      if (remmel_core > degen) then
                         occ_num_core(i)%m(m,is) = degen
                         n_occo_core(is,i) = n_occo_core(is,i) + 1
                         remmel_core = remmel_core - degen
                      else
                         occ_num_core(i)%m(m,is) = remmel_core
                         n_occo_core(is,i) = n_occo_core(is,i) + 1
                         remmel_core = 0.0_r8_kind
                      endif
                   endif
                   if (occupation_fixed_hole()) then
                      degen = num_spin * partner(i)
                      delta_hole=0.0_r8_kind
                      hole=.false.
                      call fermi_hole(i,m,is,step_dummy,hole,delta_hole)
                      occ_hole=degen-delta_hole
                      if(remmel<=zero) then
                         if(hole.and.(occ_hole>zero)) call error_handler &
                              ("FERMI/ef_start: an hole is assigned a non-zero occ.number but no electrons are left")
                         exit orbital
                      elseif ( remmel<occ_hole.and.hole ) then
                         call error_handler("FERMI/ef_start: no more electrons left to occupy hole")
                      elseif ( remmel<occ_hole.and.remmel/=zero.and..not.hole) then
                         efs = eigval(i)%m(m,is)
                         exit orbital
                      elseif ( remmel==occ_hole) then
                         efs = eigval(i)%m(m,is)
                         exit orbital
                      else
                         remmel = max( remmel - occ_hole, zero )
                         if (remmel == zero ) then
                            efs = eigval(i)%m(m,is)
                            exit orbital
                         endif
                         cycle orbital
                      endif
                   else
                      degen = num_spin * partner(i)

                      if(remmel.le.zero) then
                         exit orbital
                      elseif(remmel.eq.degen) then
                         efs = eigval(i)%m(m,is)
                         remmel = zero
                         exit orbital
                      else
                         remmel = max(remmel - degen, zero)
                         if(remmel.eq.zero) then
                            efs = eigval(i)%m(m,is)
                            exit orbital
                         endif
                         cycle orbital
                      endif
                   endif
                end if
             enddo
          enddo
       enddo

    enddo orbital

  end subroutine ef_start

  !*************************************************************

  subroutine ef_start_spin(efs)
    !----------------------------------------------------------------
    !  Purpose: produce two start fermi energies for the
    !           Newton-Raphson procedure. This desgined for the
    !           case of two separate fermi energies for
    !           each spin. For more information see header
    !           of 'ef_start'.
    !
    !  Subroutine called by: fermi_reoccup
    !  Author: FN
    !  Date: 4/96
    !------------ Modules used --------------------------------------
!    use occupation_module, only: occupation_level_sort, &
!         get_n_elec, num_list
    !------------ Declaration of formal parameters ------------------
    real(kind=r8_kind),    intent(out ) :: efs(:)
    !** End of interface *****************************************
    !------------ Declaration of subroutines ------------------------
    external error_handler
    !------------ Declaration of local variables --------------------
    real(kind=r8_kind)               :: n_elec,remmel,degen,delta_hole,&
                                        occ_hole,remmel_core,spin_diff
    integer(kind=i4_kind)            :: n_elec_int
    integer(kind=i4_kind)            :: i,m,is,index,step_dummy
    integer(i4_kind) :: n_elec_spin(2) ! FIXME: why integer?
    logical                          :: hole
    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    ! this fills the variable 'num_list', where
    ! num_list(i)%m(m,is) contains the number of the
    ! eigenvalue eigval(i)%m(m,is), here the ordering is done
    ! for each spin separately
    call occupation_level_sort(eigval,separate_spin=.true.)

    call get_n_elec(n_elec)
    n_elec_int = nint(n_elec,kind=i4_kind)

    call get_spin_diff(spin_diff)

    ! calculate the electron numbers for each spin
    n_elec_spin(1) = half * (n_elec - spin_diff) ! Spin1 = minority
    n_elec_spin(2) = n_elec - n_elec_spin(1) ! Spin2 = majority

    !---- now perform the Aufbau-Principle without actually
    !     filling the electrons into the levels
    num_orb=0
    do i=1,n_irrep                            ! this is the
       num_orb=num_orb+dim(i)*partner(i) ! total number
    enddo   ! of orbitals
    num_orb_spin = num_orb

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       num_spin = 1             ! number of el.
    else
       num_orb = num_orb * 2      ! 2 electrons fit into each orbital
       num_spin = 3 - ssym%n_spin ! number of el.
    endif

    if(n_elec_int.gt.num_orb) call error_handler &
         ( "ef_start_spin :more electrons than orbitals" )

    if (operations_core_density) call get_n_core_elec(remmel_core)
    step_dummy=1
    do is=1,ssym%n_spin

       remmel = n_elec_spin(is)
       spin_orbital : do index=1,num_orb_spin

          do i=1,n_irrep
             do m=1,dim(i)

                if (num_list(i)%m(m,is).eq.index) then
                   degen = partner(i)
                   if (operations_core_density .and. remmel_core > 0.0_r8_kind) then
                      if (remmel_core > degen) then
                         occ_num_core(i)%m(m,is) = degen
                         n_occo_core(is,i) = n_occo_core(is,i) + 1
                         remmel_core = remmel_core - degen
                      else
                         occ_num_core(i)%m(m,is) = remmel_core
                         n_occo_core(is,i) = n_occo_core(is,i) + 1
                         remmel_core = 0.0_r8_kind
                      endif
                   endif
                   if(fixed_hole) then
                      delta_hole=0.0_r8_kind
                      call fermi_hole(i,m,is,step_dummy,hole,delta_hole)
                      occ_hole=degen-delta_hole
                      if(remmel<=zero) then
                         if(hole.and.(occ_hole>zero)) call error_handler &
                              ("FERMI/ef_start: an hole is assigned a non-zero occ.number but no electrons are left")
                         exit spin_orbital
                      elseif ( remmel<occ_hole ) then
                         call error_handler("FERMI/ef_start: no more electrons left to occupy orbital")
                      elseif ( remmel==occ_hole) then
                         efs(is) = eigval(i)%m(m,is)
                         exit spin_orbital
                      else
                         remmel = max( remmel - occ_hole, zero )
                         if (remmel == zero ) then
                            efs(is) = eigval(i)%m(m,is)
                            exit spin_orbital
                         endif
                         cycle spin_orbital
                      endif
                   else
                      if (remmel.le.zero) then
                         exit spin_orbital
                      elseif ( remmel.eq.degen) then
                         efs(is) = eigval(i)%m(m,is)
                         remmel = zero
                         exit spin_orbital
                      else
                         remmel = max(remmel - degen, zero)
                         if (remmel.eq.zero) then
                            efs(is) = eigval(i)%m(m,is)
                            exit spin_orbital
                         endif
                         cycle spin_orbital
                      endif
                   endif

                endif

             enddo
          enddo

       enddo spin_orbital
    enddo
   end subroutine ef_start_spin

  !*************************************************************

  subroutine fermi_read(scfcontrol)
    ! Purpose: read in the input for the gaussian level broadening
    !          around the fermi energy.
    !          It consists of the following variables:
    !
    !         level_broadening  Can be one of 'gauss_broadening',
    !                           'fermi_broadening' or 'sinus_broadening'
    !                           This chooses different types of level
    !                           broadeing. Only the last two methods
    !                           produce an energy which is consistent
    !                           with the gradients.
    !             fermi_sigma : width of level broadening.
    !                           'gauss' : width of the gauss-function
    !                           DEFAULT: 0.3 eV
    !                           'fermi' : width of energy-range for
    !                           values of the fermi-function between
    !                           n1=0.1 and n2=0.9.
    !                           DEFAULT:
    !                           'sinus' : width of energy-range for
    !                           values of the sin-function between
    !                           n1=0.1 and n2=0.9.
    !                           DEFAULT:
    !        fermi_newton_max : max number of Newton-Rhaphson steps
    !                           to calculate e_fermi. DEFAULT: 100
    !            fermi_cutoff : cutoff value for the occupation-function:
    !                           'gauss': Only values of eigval between
    !                           +/- sqrt(2)*sigma*cutoff are taken
    !                           into account for the level broadening
    !                           DEFAULT: 5.0
    !                           'fermi': below the cutoff value the
    !                           fermi-function is set to one. If 'cutoff'
    !                           is larger than the calculated cutoff-value
    !                           ef-ln(n/1-n)/(kT) a WARNING is printed
    !                           out and the calculated value is taken
    !                           rather than the 'cutoff'
    !                           DEFAULT: calculated value
    !                           'sinus': meaningless, since the natural
    !                           cutoff is (ef-e)/(kT)=-pi/4.
    !          fermi_accuracy : accuracy to which the Newton-Rhaphson
    !                           is to converge: The criterion is
    !                           delta_efermi <= accuracy*abs(e_fermi)
    !                           DEFAULT: 1.0e-10
    !** End of interface ***************************************
    !
    !  Author: FN
    !  Date: 4/96
    ! ---- Modules ---------------------------------------------
    use input_module
#if 0
    use options_module, only: update_hessian_iteration
#endif
    logical,optional,intent(in) :: scfcontrol
    integer(kind=i4_kind) :: unit,status
    logical check_sigma

    if (present(scfcontrol)) then
       check_sigma = .not.scfcontrol
    else
       check_sigma = .true.
    endif


    !these are the recommended default values:
    ! get a common fermi energy for both up
    ! and down spin
    fermi_level_broad     = df_fermi_level_broad
    fermi_gauss           = df_fermi_gauss
    fermi_fermi           = df_fermi_fermi
    fermi_sinus           = df_fermi_sinus
    fermi_sigma           = df_fermi_sigma
    fermi_energy_range    = df_fermi_e_range
    fermi_newton_max      = df_fermi_newton_max
    fermi_cutoff          = df_fermi_cutoff
    fermi_accuracy        = df_fermi_accuracy

    if ( input_line_is_namelist("fermi") ) then
          call input_read_to_intermediate
          unit= input_intermediate_unit()
          read(unit, nml=fermi, iostat=status)
          if (status .gt. 0) call input_error( &
               "fermi_read: namelist fermi")
#if 0 /* level_broad not  implemented for dervs */
      if (update_hessian_iteration .gt. 0) then
         fermi_level_broad     =.false.
         fermi_gauss=.false.
         fermi_fermi=.false.
         fermi_energy_range    = 0.0_r8_kind
      endif
#endif
    end if
    ! to ensure compatibility
    if (fermi_level_broad) fermi_gauss = .true.

    status = 0
    if (fermi_gauss) then
       level_broadening = gauss_broadening
       status = status + 1
    endif
    if (fermi_fermi) then
       level_broadening = fermi_broadening
       status = status + 1
    endif
    if (fermi_sinus) then
       level_broadening = sinus_broadening
       status = status + 1
    endif
    fermi_level_broad = status > 0
    if (status > 1) then
       call input_error &
            ("fermi_read: please choose only ONE method for level broadening")
    endif


    ! in case of NO level broadening, set some switches to
    ! false and exit

    if (.not.fermi_level_broad) then
!       fermi_fix_up_and_down = .false.
!       fermi_common_ef=.false.
!       fermi_unpaired=0
       return
    endif

    ! ensure that either fermi_sigma OR fermi_energy_range is set
    ! OR that the two are consistent


    if (fermi_sigma == df_fermi_sigma ) then
       if (fermi_energy_range == df_fermi_e_range) then
          call input_error &
               ("fermi_read: please specify either a SIGMA or an ENERGY_RANGE")
       else
          fermi_sigma = fermi_energy_range/g_fac
       endif
    else ! fermi_sigma /= df_fermi_sigma
       if (fermi_energy_range == df_fermi_e_range) then
          fermi_energy_range = fermi_sigma*g_fac
       elseif (abs(fermi_sigma-fermi_energy_range/g_fac)>=sigma_tol) then
          if (check_sigma) call input_error &
               ("fermi_read: ENERGY_RANGE and SIGMA are not consistent.")
       endif
    endif

!    if (fermi_fix_up_and_down) fermi_common_ef = .false.

!    if( fermi_unpaired.lt.0) call input_error &
!         ("fermi_read: please supply a positive value fermi_unpaired")

!    if (fermi_unpaired.ne.0.and.fermi_common_ef) call input_error &
!         ("fermi_read: either unpaired elecs or common_ef ")

    if (fermi_sigma.le.zero) call input_error &
         ("fermi_read: please supply positive sigma ")

    if (fermi_newton_max.le.0) call input_error &
         ("fermi_read: please supply newton_max >= 0")

    if (fermi_cutoff.lt.zero) call input_error &
         ("fermi_read : please supply a positive fermi_cutoff")

    if (fermi_accuracy.lt.zero) call input_error &
         ("fermi_read : please supply a positive fermi_accuracy ")

    return
  end subroutine fermi_read

  !*************************************************************

  subroutine fermi_read_scfcontrol(warning,new_trace_format,new_scfcontrol_data)
    ! purpose: read in the namelist fermi from the
    !          input file  in $TTFSSTART/scfcontrol
    !          at each scf-cycle, giving a notice if
    !          parameters were changed.
    use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters -------------
    logical, intent(in   ), optional :: warning ! if present and .false. warning
                                                ! messages are suppressed
    logical, intent(inout), optional :: new_trace_format,new_scfcontrol_data
    !** End of interface ***************************************
    logical :: &
         old_fermi_level_broad    , &
         old_fermi_gauss          , &
         old_fermi_fermi          , &
         old_fermi_sinus
!         old_fermi_fix_up_and_down, &
!         old_fermi_common_ef
    integer(kind=i4_kind) :: &
!         old_fermi_unpaired  , &
         old_fermi_newton_max
    real(kind=r8_kind) :: &
         old_fermi_energy_range, &
         old_fermi_sigma       , &
         old_fermi_cutoff      , &
         old_fermi_accuracy
    logical :: changed,changed_e,changed_s

    old_fermi_level_broad     = fermi_level_broad
    old_fermi_gauss           = fermi_gauss
    old_fermi_fermi           = fermi_fermi
    old_fermi_sinus           = fermi_sinus
!    old_fermi_fix_up_and_down = fermi_fix_up_and_down
!    old_fermi_common_ef       = fermi_common_ef
!    old_fermi_unpaired        = fermi_unpaired
    old_fermi_sigma           = fermi_sigma
    old_fermi_energy_range    = fermi_energy_range
    old_fermi_newton_max      = fermi_newton_max
    old_fermi_cutoff          = fermi_cutoff
    old_fermi_accuracy        = fermi_accuracy

    call fermi_read(scfcontrol=.true.)

    if (present(new_trace_format)) then
       new_trace_format = new_trace_format .or. &
            ( old_fermi_level_broad .neqv. fermi_level_broad )
    endif

    if (present(warning)) then
       if (.not.warning) return
    endif

    changed = check_logical("fermi_level_broad", &
         fermi_level_broad,old_fermi_level_broad)
    changed = check_logical("fermi_gauss", &
         fermi_gauss,old_fermi_gauss)
    changed = check_logical("fermi_fermi", &
         fermi_fermi,old_fermi_fermi)
    changed = check_logical("fermi_sinus", &
         fermi_sinus,old_fermi_sinus)
!    changed = check_logical("fermi_fix_up_and_down", &
!         fermi_fix_up_and_down,old_fermi_fix_up_and_down)
!    changed = check_logical("fermi_common_ef", &
!         fermi_common_ef,old_fermi_common_ef)
!    changed = check_intg("fermi_unpaired", &
!         fermi_unpaired,old_fermi_unpaired)
    changed = check_intg("fermi_newton_max", &
         fermi_newton_max,old_fermi_newton_max)

    changed_e = check_real("fermi_energy_range", &
         fermi_energy_range,old_fermi_energy_range)

    changed_s = check_real("fermi_sigma", &
         fermi_sigma,old_fermi_sigma)

    changed = check_real("fermi_cutoff", &
         fermi_cutoff,old_fermi_cutoff)
    changed = check_real("fermi_accuracy", &
         fermi_accuracy,old_fermi_accuracy)

    if (changed_e .and. changed_s) then
       if (abs(fermi_sigma-fermi_energy_range/g_fac)>=sigma_tol) then
          call error_handler("fermi_read_scfcontrol: &
               &inconsistent simultaneous change of sigma and energy_range")
       endif
    elseif (changed_e) then
       fermi_sigma = fermi_energy_range/g_fac
    elseif (changed_s) then
       fermi_energy_range = fermi_sigma*g_fac
    endif

    if (present(new_scfcontrol_data)) then
       new_scfcontrol_data = new_scfcontrol_data .or. &
            ( changed_e .neqv. changed_s )
    endif

  contains

    logical function check_intg(name,val,val_old)
      character(len=*)     , intent(in) :: name
      integer(kind=i4_kind), intent(in) :: val, val_old
      check_intg=.false.
      if ( val .ne. val_old ) then
         call write_to_output_units("fermi_read_scfcontrol: "// &
              name//" was altered. New value: ",inte=val)
         call write_to_trace_unit("convergence_read_scfcontrol: "// &
              name//" was altered. New value: ",inte=val)
         check_intg=.true.
       endif
     end function check_intg

     logical function check_logical(name,val,val_old)
       character(len=*)     , intent(in) :: name
       logical, intent(in) :: val, val_old
       check_logical=.false.
       if ( val .neqv. val_old ) then
          if (val) then
             call write_to_output_units("fermi_read_scfcontrol: Flag "// &
                  name//" was switched on.")
             call write_to_trace_unit("fermi_read_scfcontrol: Flag "// &
                  name//" was switched on.")
          else
             call write_to_output_units("fermi_read_scfcontrol: Flag "// &
                  name//" was switched off.")
             call write_to_trace_unit("fermi_read_scfcontrol: Flag "// &
                  name//" was switched off.")
          endif
          check_logical=.true.
       endif
     end function check_logical

     logical function check_real(name,val,val_old)
       character(len=*)  , intent(in) :: name
       real(kind=r8_kind), intent(in) :: val, val_old
       check_real=.false.
       if ( val .ne. val_old ) then
          call write_to_output_units("fermi_read_scfcontrol: "// &
               name//" was altered. New value: ",re=val)
          call write_to_trace_unit("fermi_read_scfcontrol: "// &
               name//" was altered. New value: ",real=val)
          check_real=.true.
       endif
     end function check_real

  end subroutine fermi_read_scfcontrol

  !*************************************************************

  subroutine newton(e_start,e_end)
    ! Purpose: Newton-Rhaphson procedure to calculate the
    !          Fermi energy.
    !          This procedure seeks the zeros of the
    !          function L(ef) where
    !          L(ef) = 1/2 sum_i[ m_i (1+erfc(arg_i)) ] - N
    !
    !          i : Orbitals
    !          m_i: degeneracy of orbital i
    !          N  : total number of electrons
    !          arg_i : (ef - e_i/(sqrt(2)*sigma))
    !
    !------------ Modules used ----------------------------------
!    use occupation_module, only: get_n_elec, occ_num, &
!         alloc_occ_num,n_occo,alloc_n_occo
    use init_module, only: init
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in)  :: e_start
    real(kind=r8_kind), intent(out) :: e_end
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)     :: step,i,m,is,counter
    real(kind=r8_kind)        :: ef,degen,fn,dfn, &
         df, deltan, delta_e, n_elec, delta_hole
    real(kind=r8_kind)           :: summ
    real(kind=r8_kind),parameter :: small = 1.e-8_r8_kind
    real(kind=r8_kind)           :: gamma_plus,gamma_minus,&
         remmel,stepmax,e_diff
    real(kind=r8_kind),parameter :: max_step = 0.0918724_r8_kind ! =2.5 eV
    real(kind=r8_kind),parameter :: atomic_unit = 27.211652_r8_kind
    real(kind=r8_kind)           :: ec_lower,ec_upper
    logical                      :: hole

    external error_handler

    ASSERT(n_irrep/=-1)

    e_end = zero

    call get_n_elec(n_elec)

    if (.not.allocated(occ_num)) then
       call alloc_occ_num()
    endif

    if (.not.allocated(n_occo)) then
       call alloc_n_occo(ssym)
    endif

    ! num_orb and num_spin should have already
    ! been calculated in 'ef_start'

    ef = e_start

    newton_loop : do step = 1,fermi_newton_max
       call init(n_occo)

       call get_cutoff(ef,ec_lower,ec_upper)
       ! calculate L/L`
       deltan = zero
       dfn = zero
       df = zero
       delta_e = zero
       summ = zero
       if (output_fermi_newton) then
          write (output_unit,*)&
               "Newton: ----- Newton-loop nr. ",step
          write (output_unit,*)&
               "Newton:  current E_fermi [eV] = ", ef * atomic_unit
       endif

       counter=1
       irrep : do i=1,n_irrep
          element : do m=1,dim(i)
             spin : do is=1,ssym%n_spin

                degen = num_spin * partner(i)

                delta_hole=0.0_r8_kind
                hole=.false.
                if (occupation_fixed_hole()) then
                   call fermi_hole(i,m,is,step,hole,delta_hole)
                endif

                ! provide an 'Abschneidekriterium'
                !  case 1: eigenvalue is greater than
                !  e_fermi + cutoff*sqpi2*sigma
                !  since eigenvalues are not sorted we
                !  cannot exit here
                if (eigval(i)%m(m,is).gt.ec_upper) then
                   if (hole.and.(delta_hole/=degen)) call error_handler &
                        ("FERMI_REOCCUP:Newton : no holes with nonzero occupation number at above fermi level")
                   occ_num(i)%m(m,is) = zero
                   fn = zero
                   dfn =zero
                   cycle spin
                elseif(eigval(i)%m(m,is).le.ec_upper .and. &
                     eigval(i)%m(m,is).gt.ec_lower) then
                   ! case 2: eigenvalue is between e_fermi
                   ! +cutoff*... and e_fermi - cutoff*...
                   fn=lfunc(ef,eigval(i)%m(m,is))
                   dfn=dlfunc(ef,eigval(i)%m(m,is))
                elseif ( eigval(i)%m(m,is).le.ec_lower) then

                   ! case 3: eigenvalue is below
                   ! e_fermi - cutoff*sqpi2*sigma
                   fn = one
                   dfn = zero
                endif

                if(delta_hole>degen) call error_handler &
                     ("FERMI_REOCCUP/Newton : level is occupied with negative number of electrons")
                fn=fn*(degen-delta_hole)
                dfn=dfn*(degen-delta_hole)

                occ_num(i)%m(m,is) = fn
                n_occo(is,i)=n_occo(is,i) + 1
                deltan = deltan + fn
                df = df + dfn

             enddo spin
          enddo element
       enddo irrep


       ! calculate the energetic distances gamma_plus and
       ! gamma_minus between ef and the nearest
       ! eigenvalue both for higher and for lower lying levels
       gamma_plus  = HUGE(zero)
       gamma_minus = HUGE(zero)
       irrp: do i=1,n_irrep
          elemt: do m=1,dim(i)
             spn: do is=1,ssym%n_spin
                e_diff = eigval(i)%m(m,is) - ef
                if( e_diff .eq. zero ) cycle ! the case when ef==E(HOMO)
                if (e_diff .gt. zero) then
                   gamma_plus  = min( e_diff,gamma_plus)
                else
                   gamma_minus = min(-e_diff,gamma_minus)
                endif
             enddo spn
          enddo elemt
       enddo irrp
       ! test, if gamma_plus and gamma_minus got new values:
       if ( (gamma_plus .eq. HUGE(zero)) .or. &
            (gamma_minus .eq. HUGE(zero)) ) &
            call error_handler &
            ("NEWTON: E_fermi lies outside of the spectrum")

       ! remmel = number of remaining electrons
       remmel = -(deltan - n_elec)
       ! determine proper step restrictions
       if (remmel .gt. zero) then
          ! ef has to be shifted towards higher energies
          stepmax = gamma_plus
       elseif (remmel .lt. zero) then
          ! ef has to be shifted towards lower energies
          stepmax = gamma_minus
       else
          ! remmel=0
          ! Newton procedure should converge after next step
          stepmax = max_step
       endif
       ! check, if stepmax <= max_step
       if (stepmax .gt. max_step) then
          ! restrict stepmax to the value of max_step
          stepmax=max_step
       endif

       ! now calculate delta_e
       if (df .gt. zero) then
          delta_e = remmel/df  ! df is always positive
          if (abs(delta_e).gt. stepmax) then
             ! newton step too long,
             ! set step length equal to stepmax
             delta_e = sign(stepmax,remmel)
          endif
       elseif (df .eq. zero) then
          ! no levels affected by broadening
          if (abs(remmel) .lt. small) then
             ! remmel < small AND no level affected by broadening:
             ! set ef exactly to (LUMO-HOMO)/2
             delta_e = (gamma_plus - gamma_minus)/2.0_r8_kind
          else ! abs(remmel) => small
             delta_e = sign(stepmax,remmel)
          endif
       else ! df < 0 error!
          call error_handler ("Newton: unexpected value: df < 0")
       endif

       if (output_fermi_newton) then
          write(output_unit,*)&
               "Newton: Number of remaining electrons: remmel = ",remmel
          write(output_unit,*)&
               "Newton: delta_e[eV] = ", delta_e * atomic_unit
       endif

       ef = ef + delta_e

       if ( (abs(delta_e).lt.small*abs(ef)) .and. &
            (abs(remmel).lt.small) ) then
          if (output_fermi) call write_to_output_units &
               ("Newton: convergence reached after cycle",step)
          exit newton_loop
       endif

    enddo newton_loop

    if (output_fermi_newton) then
       write(output_unit,*)&
            "Newton: final E_fermi[eV] = ", ef * atomic_unit
    endif

    if (output_fermi) then
       write(output_unit,*)"NEWTON: Indices of highest occupied MOs per Irrep and Spin"
       write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
       do i=1,symmetry_data_n_irr()
          write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo(is,i),is=1,ssym%n_spin)
       enddo
       write(output_unit,*)" "
       if(operations_core_density)then
          write(output_unit,*)"FERMI_START: Indices of highest occupied core orbital per Irrep and Spin"
          write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
          do i=1,symmetry_data_n_irr()
             write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo_core(is,i),is=1,ssym%n_spin)
!:TST[ occ_num_core
!            do m=1,maxval(n_occo_core(:,i))
!               write(output_unit,*)(occ_num_core(i)%m(m,is),is=1,ssym%n_spin)
!            enddo
!:TST]
          enddo
          write(output_unit,*)" "
       endif
    endif

    e_end = ef

    return
  end subroutine newton

  !*************************************************************

  subroutine newton_spin(e_start_spin,e_end_spin)
    ! Purpose: Newton-Rhaphson procedure to calculate the
    !          Fermi energy for each spin separately.
    !          Either the two spin sorts are simply kept
    !          fixed separately (fermi_unpaired = 0) or
    !          additionally a number of unpaired electrons
    !          can be specified (fermi_unpaired > 0).
    !
    !------------ Modules used ----------------------------------
!    use occupation_module, only: get_n_elec, occ_num, &
!         alloc_occ_num,n_occo,alloc_n_occo
    use init_module, only: init
    use machineparameters_module, only: machineparameters_DimCheck
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in)  :: e_start_spin(:)
    real(kind=r8_kind), intent(out) :: e_end_spin(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)     :: step,i,m,is,alloc_stat
    real(kind=r8_kind)        :: degen,fn,dfn, &
         df,deltan,delta_e,n_elec,delta_hole
    real(r8_kind), allocatable :: ef(:), n_elec_spin(:) ! Real here!
    real(kind=r8_kind)           :: ec_upper,ec_lower,spin_diff,&
         gamma_plus,gamma_minus,stepmax,remmel_m(2),e_diff
    real(kind=r8_kind),parameter :: small = 1.e-8_r8_kind
    real(kind=r8_kind),parameter :: max_step = 0.0918724_r8_kind ! =2.5 eV
    real(kind=r8_kind),parameter :: atomic_unit = 27.211652_r8_kind
    logical                      :: hole

    external error_handler

    ASSERT(n_irrep/=-1)

    call get_n_elec(n_elec)

    if (machineparameters_DimCheck) then
       if (ssym%n_spin.le.1.0) call error_handler &
            ("newton_spin: n_spin is not 2 !!!")
       if (size(e_start_spin,1).ne.ssym%n_spin) &
            call error_handler("newton_spin: dim 1 wrong")
       if (size(e_end_spin,1).ne.ssym%n_spin) &
            call error_handler("newton_spin: dim 2 wrong")
    endif
    if (.not.allocated(occ_num)) &
         call alloc_occ_num()
    if (.not.allocated(n_occo)) &
         call alloc_n_occo(ssym)

    allocate(ef(ssym%n_spin),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("newton_spin: allocation 1 failed")
    allocate(n_elec_spin(ssym%n_spin),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("newton_spin: allocation 2 failed")

    ! get the number of unpaired electrons defined by
    ! the user through the Input parameter MAGN_MOMENT
    ! (namelist OCCUPATION)
    call get_spin_diff(spin_diff)

    ! calculate the electron numbers for each spin
    n_elec_spin(1) = half * (n_elec - spin_diff) ! Spin 1 = minority
    n_elec_spin(2) = n_elec - n_elec_spin(1) ! Spin 2 = majority

    spin: do is=1,ssym%n_spin
       ef(is) = e_start_spin(is)
       if (output_fermi_newton) then
             write(output_unit,*)&
                  "Newton_Spin: perform Newton procedure for Spin ",is
          endif

       newton_loop: do step =1,fermi_newton_max

          call init(n_occo(is,:))

          call get_cutoff(ef(is),ec_lower,ec_upper)
          deltan = zero
          dfn = zero
          df = zero
          delta_e = zero
          delta_hole=0.0_r8_kind

          if (output_fermi_newton) then
             write (output_unit,*)&
                  "Newton_Spin: ----- Newton-loop nr. ",step
             write (output_unit,*)&
                  "Newton_Spin:  current E_fermi [eV] = ", ef(is) * atomic_unit
          endif

          hole = .false.  ! added (TS)
          if (occupation_fixed_hole()) then
             call fermi_hole(i,m,is,step,hole,delta_hole)
          endif
          irrep: do i=1,n_irrep
             element: do m=1,dim(i)
                degen = num_spin * partner(i)
                if (eigval(i)%m(m,is).gt.ec_upper) then
                   if (hole.and.(delta_hole/=degen)) call error_handler &
                        ("FERMI_REOCCUP:Newton_spin : no holes with nonzero occupation number at above fermi level")
                   occ_num(i)%m(m,is) = zero
                   fn = zero
                   dfn =zero
                   cycle element
                elseif(eigval(i)%m(m,is).le.ec_upper .and. &
                     eigval(i)%m(m,is).gt.ec_lower) then
                   fn = lfunc(ef(is),eigval(i)%m(m,is))
                   dfn = dlfunc(ef(is),eigval(i)%m(m,is))
                elseif( eigval(i)%m(m,is).le.ec_lower ) then
                   fn = one
                   dfn = zero
                endif
                if(delta_hole>real(degen,kind=r8_kind)) call error_handler &
                     ("FERMI_REOCCUP/Newton_spin: level will be occupied with negative number of electrons")

                fn=fn*(degen-delta_hole)
                dfn=dfn*(degen-delta_hole)

                occ_num(i)%m(m,is) = fn
                n_occo(is,i) = n_occo(is,i) + 1
                deltan = deltan + fn
                df = df + dfn

             enddo element
          enddo irrep

          ! calculate the energetic distances gamma_plus and
          ! gamma_minus between ef(is) and the nearest
          ! eigenvalue both for higher and for lower lying levels
          gamma_plus = zero
          gamma_minus = zero
          remmel_m(is) = -(deltan - n_elec_spin(is))
          irrp: do i=1,n_irrep
             elemt: do m=1,dim(i)
                e_diff = eigval(i)%m(m,is) - ef(is)
                if (e_diff .gt. zero) then
                   if ( (e_diff .lt. gamma_plus) .or. &
                        (gamma_plus .eq. zero) ) gamma_plus = e_diff
                elseif (e_diff .lt. zero) then
                   e_diff = -e_diff
                   if ( (e_diff .lt. gamma_minus).or. &
                        (gamma_minus .eq. zero) ) gamma_minus = e_diff
              ! else-part: e_diff=0. In this case gamma_plus and
              ! gamma_minus are left unchanged
                endif
             enddo elemt
          enddo irrp
          ! test, if gamma_plus and gamma_minus got new values:
          if ( (gamma_plus .eq. zero) .or. &
               (gamma_minus .eq. zero) ) &
               call error_handler &
               ("Newton_spin: E_fermi lies outside of the spectrum")


          ! determine proper step restrictions
          stepmax=max_step
          if (remmel_m(is) .gt. zero) then
             ! ef has to be shifted towards higher energies
             stepmax = gamma_plus
          elseif (remmel_m(is) .lt. zero) then
             ! ef has to be shifted towards lower energies
             stepmax = gamma_minus
          else
             ! remmel_m(is)=0
          endif
          ! check, if stepmax <= max_step
          if (stepmax .gt. max_step) then
             stepmax=max_step
          endif

          ! now calculate delta_e
          if (df .gt. zero) then
             delta_e = remmel_m(is)/df  ! df is always positive
             if (abs(delta_e).gt. stepmax) then
                ! newton step too long,
                ! set step length equal to stepmax
                delta_e = sign(stepmax,remmel_m(is))
             endif
          elseif (df .eq. zero) then
             ! no levels affected by broadening
             if (abs(remmel_m(is)) .lt. small) then
                ! remmel = 0 AND no level affected by broadening:
                ! set ef exactly to (LUMO-HOMO)/2
                delta_e = (gamma_plus - gamma_minus)/2.0_r8_kind
             else ! abs(remmel) > small
                delta_e = sign(stepmax,remmel_m(is))
             endif
          else ! df < 0 error!
             call error_handler ("Newton_spin: unexpected value: df < 0")
          endif

          if (output_fermi_newton) then
             write(output_unit,*)&
                  "Newton_spin: Spin ",is
             write(output_unit,*)&
                  "Newton_spin: Number of remaining electrons: remmel_m = "&
                  ,remmel_m(is)
             write(output_unit,*)&
                  "Newton_spin: delta_e[eV] = ", delta_e * atomic_unit
          endif

          ef(is) = ef(is) + delta_e

          if ( (abs(delta_e).lt.small*abs(ef(is))) .and. &
               (abs(remmel_m(is)).lt.small) ) then
             if (output_fermi) call write_to_output_units &
                  ("Newton_spin: convergence reached after cycle",step)
             exit newton_loop
          endif
       enddo newton_loop

       e_end_spin(is) = ef(is)

       if (output_fermi_newton) then
          write(output_unit,*)&
               "Newton: final E_fermi[eV] = ", ef * atomic_unit
       endif
    enddo spin

    if (output_fermi) then
       write(output_unit,*)"NEWTON: Indices of highest occupied MOs per Irrep and Spin"
       write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
       do i=1,symmetry_data_n_irr()
          write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo(is,i),is=1,ssym%n_spin)
       enddo
       write(output_unit,*)" "
       if(operations_core_density)then
          write(output_unit,*)"FERMI_START: Indices of highest occupied core orbital per Irrep and Spin"
          write(output_unit,'(" Irrep ",2(a8,i2))')("Spin   ",is,is=1,ssym%n_spin)
          do i=1,symmetry_data_n_irr()
             write(output_unit,'(i4,4x,2(3x,i5))')i,(n_occo_core(is,i),is=1,ssym%n_spin)
!:TST[ occ_num_core
!            do m=1,maxval(n_occo_core(:,i))
!               write(output_unit,*)(occ_num_core(i)%m(m,is),is=1,ssym%n_spin)
!            enddo
!:TST]
          enddo
          write(output_unit,*)" "
       endif
    endif


    deallocate(ef,n_elec_spin,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("newton_spin : deallocation failed ")

  end subroutine newton_spin

  !*************************************************************

  function erfc(x)
    ! Purpose : calculate the complementary error-function
    !           using the chebyshev approximation
    !
    ! References: YL Luke 1975, pp 123-4
    !             Num. Rec. 1988 "CHEBEV"
    !
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind), intent(in) :: x
    real(kind=r8_kind)             :: erfc
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind),parameter :: na=25,nc=22
    real(kind=r8_kind)              :: a(0:na),c(0:nc)
    real(kind=r8_kind)              :: d,dd,alpha,z,sv
    integer(kind=i4_kind)           :: j
    real(kind=r8_kind),parameter    :: sqpi2_loc=1.128379167095513e0_r8_kind

    data a / &
         .109547129977762e+1_r8_kind, -.289175401126989e+0_r8_kind, &
         .110456398633795e+0_r8_kind, -.412531882278565e-1_r8_kind, &
         .140828380706516e-1_r8_kind, -.432929544743143e-2_r8_kind, &
         .119827190159228e-2_r8_kind, -.299972962353249e-3_r8_kind, &
         .683258603788747e-4_r8_kind, -.142469884548677e-4_r8_kind, &
         .273540877283989e-5_r8_kind, -.048619128719754e-5_r8_kind, &
         .008038727621172e-5_r8_kind, -.001241841831213e-5_r8_kind, &
         .000179953258879e-5_r8_kind, -.000024547948775e-5_r8_kind, &
         .000003162508603e-5_r8_kind, -.000000385902200e-5_r8_kind, &
         .000000044720291e-5_r8_kind, -.000000004933613e-5_r8_kind, &
         .000000000519303e-5_r8_kind, -.000000000052258e-5_r8_kind, &
         .000000000005037e-5_r8_kind, -.000000000000466e-5_r8_kind, &
         .000000000000041e-5_r8_kind, -.000000000000004e-5_r8_kind /
    DATA c / &
         .975083423708556e+0_r8_kind, -.240493938504146e-1_r8_kind, &
         .820452240880432e-3_r8_kind, -.434293081303427e-4_r8_kind, &
         .301844703403493e-5_r8_kind, -.025447331925082e-5_r8_kind, &
         .002485835302051e-5_r8_kind, -.000273172013238e-5_r8_kind, &
         .000033084722797e-5_r8_kind, -.000004350549080e-5_r8_kind, &
         .000000614121457e-5_r8_kind, -.000000092236928e-5_r8_kind, &
         .000000014635665e-5_r8_kind, -.000000002439278e-5_r8_kind, &
         .000000000424976e-5_r8_kind, -.000000000077084e-5_r8_kind, &
         .000000000014507e-5_r8_kind, -.000000000002824e-5_r8_kind, &
         .000000000000567e-5_r8_kind, -.000000000000117e-5_r8_kind, &
         .000000000000025e-5_r8_kind, -.000000000000005e-5_r8_kind, &
         .000000000000001e-5_r8_kind /

    d = zero
    dd = zero
    erfc = zero
    if (abs(x).lt.three) then
       !         CALCULATE VIA ERF
       z=x/three
       alpha=two-four*z*z
       do  J=NA,0,-1
          SV=D
          D=-ALPHA*D-DD+A(J)
          DD=SV
       enddo

       erfc=ONE-sqpi2_loc*Z*(D-DD)
    else
       !      CALCULATE DIRECTLY
       Z=abs(THREE/X)
       ALPHA=TWO-FOUR*Z*Z
       do J=NC,0,-1
          SV=D
          D=-ALPHA*D-DD+C(J)
          DD=SV
       enddo
       if (X.GT.ZERO) then
          ERFC=HALF*exp(-X*X)/X*(D+HALF*ALPHA*DD)*sqpi2_loc
       else
          ERFC=TWO-HALF*exp(-X*X)/(-X)*(D+HALF*ALPHA*DD)*sqpi2_loc
       endif
    endif

  end function erfc

  !*************************************************************

  subroutine fermi_write_input(iounit,scfcontrol)
    !
    ! Purpose: write the fermi namelist with its default values to the
    ! input.out file.
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full, real_format1, real_format2
    implicit none
    integer, intent(in) :: iounit
    logical,optional,intent(in) :: scfcontrol
    !** End of interface ***************************************

    integer(kind=i4_kind) :: level
    logical :: scfcontrol_

    real_format1 = '("    ",a," = ", f9.5:" # ",a)'
    real_format2 = '("    ",a," = ",es9.3:" # ",a)'

    level = operations_echo_input_level
    scfcontrol_ = .false.
    if (present(scfcontrol)) then
       scfcontrol_ = scfcontrol
       if (scfcontrol) level = echo_level_full
    endif

    call start("FERMI","FERMI_WRITE_INPUT",iounit,level)
    call flag("FERMI_GAUSS          ",fermi_gauss          ,df_fermi_gauss     )
    call flag("FERMI_FERMI          ",fermi_fermi          ,df_fermi_fermi     )
    call flag("FERMI_SINUS          ",fermi_sinus          ,df_fermi_sinus     )
    if (.not.scfcontrol_) then
    call real("FERMI_ENERGY_RANGE   ",fermi_energy_range   ,df_fermi_e_range ,1)
    endif
    call real("FERMI_SIGMA          ",fermi_sigma          ,df_fermi_sigma   ,1)
    call real("FERMI_CUTOFF         ",fermi_cutoff         ,df_fermi_cutoff  ,1)
!    call flag("FERMI_FIX_UP_AND_DOWN",fermi_fix_up_and_down,df_fermi_up_down   )
!    call flag("FERMI_COMMON_EF      ",fermi_common_ef      ,df_fermi_common_ef )
!    call intg("FERMI_UNPAIRED       ",fermi_unpaired       ,df_fermi_unpaired  )
    call intg("FERMI_NEWTON_MAX     ",fermi_newton_max     ,df_fermi_newton_max)
    call real("FERMI_ACCURACY       ",fermi_accuracy       ,df_fermi_accuracy,2)
    if(.not.scfcontrol_) then
       call stop()
    else
       call stop(empty_line=.false.)    !!!!!!!!!!!!!!!!AS
    end if

  end subroutine fermi_write_input

  !*************************************************************

  subroutine read_eigvals()
    ! purpose: just for debugging. This routine read the file
    !          $TTFSDATADIR/eigevalues.dat, which contains the
    !          the eigenvalues of the old program and copies
    !          them to the eigen_data_module.
    !** End of interface *****************************************
    use iounitadmin_module
    use filename_module, only: tmpfile
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: io_u,i,m,is,counter
    type(arrmat2),allocatable :: eigval_help(:)
    !------------ Declaration of subroutines used ----------------

    io_u = get_iounit()
    open(unit=io_u,status='old',form='formatted', &
         file=trim(tmpfile('eigenvalues.dat')))

    allocate(eigval_help(ssym%n_irrep))
    do i=1,ssym%n_irrep
       allocate(eigval_help(i)%m(ssym%dim(i),ssym%n_spin))
    enddo
    counter=1
    do i=1,ssym%n_irrep
       do m=1,ssym%dim(i)
          do is=1,ssym%n_spin

             read(io_u,*) eigval_help(i)%m(m,is)
             !             print*,' Robital No. ',counter,' Eigval ',&
             !                  eigval_help(i)%m(m,is)
             counter=counter+1

          enddo
       enddo
       eigval(i)%m = eigval_help(i)%m
    enddo

    deallocate(eigval_help)
    return
  end subroutine read_eigvals
  !*************************************************************

  subroutine get_cutoff(ef,ec_lower,ec_upper)
    ! Purpose: depending on the chosen 'level_broadening' either return
    !          the user specified cutoff or calculate a cutoff
    !          for the fermi-function with:
    !          e_cutoff = ef - ln(n/(1-n)*kT
    !          where kT = 1/get_sigma(), n=0.99999
    ! -------------------------------------------------------
    real(kind=r8_kind),intent(in)  :: ef
    real(kind=r8_kind),intent(out) :: ec_lower, ec_upper
    real(kind=r8_kind) :: n
    ! -------------------------------------------------------
    select case (level_broadening)
    case (gauss_broadening)
       ec_upper = ef + sqrt2*get_sigma()*fermi_cutoff
       ec_lower = ef - sqrt2*get_sigma()*fermi_cutoff
    case (fermi_broadening)
       n=(1.0_r8_kind - 1e-12_r8_kind)  ! formerly 1 - 1e-6
       ec_lower = ef - log(n/(1.0_r8_kind-n))/get_sigma()
       n=1.0e-12_r8_kind  ! formerly 1e-6
       ec_upper = ef - log(n/(1.0_r8_kind-n))/get_sigma()
    case (sinus_broadening)
       ec_upper = ef + pi/(2.0_r8_kind*get_sigma())
       ec_lower = ef - pi/(2.0_r8_kind*get_sigma())
    case default
       call error_handler ("get_cutoff : no level broadening method is set")
    end select
    return
  end subroutine get_cutoff

  !*************************************************************

  function get_sigma()
    ! Purpose: depending on the chosen method the measure of
    !          the broadening is returned:
    !          'gauss' : fermi_sigma, as supplied in the input
    !                    including unit conversion
    !          'fermi' : 90%-10% value, conversion from energy
    !                    in eV
    !          'sinus' : same as 'fermi'
    ! --------------------------------------------------------
    real(kind=r8_kind)   :: get_sigma
    !** End of interface *****************************************
    real(kind=r8_kind)   :: n1,n2,factor
    real(kind=r8_kind), parameter :: atomic_unit = 27.211652_r8_kind
    select case (level_broadening)
    case (gauss_broadening)
       get_sigma = fermi_sigma/atomic_unit
    case (fermi_broadening)
       n1=0.9_r8_kind
       n2=0.1_r8_kind
       factor = log((1.0_r8_kind-n2)/n2) - log((1.0_r8_kind-n1)/n1)
       get_sigma = factor / (fermi_energy_range/atomic_unit)
    case (sinus_broadening)
       n1=0.9_r8_kind
       n2=0.1_r8_kind
       factor = asin(one-two*n2) - asin(one-two*n1)
       get_sigma = factor / (fermi_energy_range/atomic_unit)
    case default
       call error_handler ("get_sigma: no broadening method specified")
    end select
    return
  end function get_sigma

  !*************************************************************

  function lfunc(ef,ei)
    ! Purpose : depending on the chosen method lfunc returns
    !           the value of the function to be evaluated for
    !           determining e_fermi and occ_num:
    !           'gauss': lfunc = 1.0-0.5*erfc((ef-ei)/(sqrt(2)*sigma)
    !           'fermi': lfunc = (1+exp(-(ef-ei)/sigma))**(-1)
    !           'sinus': lfunc = (...)
    ! ---------------------------------------------------------
    real(kind=r8_kind), intent(in) :: ef,ei
    real(kind=r8_kind)             :: lfunc
    !** End of interface *****************************************
    real(kind=r8_kind)             :: arg

    select case (level_broadening)
    case (gauss_broadening)
       arg = (ef-ei)/(sqrt2*get_sigma())
       lfunc=one-half*erfc(arg)
    case (fermi_broadening)
       arg = (ef-ei)*get_sigma()
       lfunc = one/(one + exp(-arg))
    case (sinus_broadening)
       lfunc = half*(one - sin(get_sigma()*(ei-ef)))
    end select

  end function lfunc

  !*************************************************************

  function dlfunc(ef,ei)
    ! Purpose: depending on the chosen method 'dlfunc' returns
    !          the value of the first derivative of 'lfunc'
    !          to be evaluated for determining e_fermi:
    !          'gauss': 1/(sqpi2)*exp(-[(ef-ei)/(sqrt(2)*sigma)]**2)
    !          'fermi': ...
    !          'sinus': ...
    ! ------------------------------------------------------
    real(kind=r8_kind), intent(in) :: ef,ei
    real(kind=r8_kind)             :: dlfunc
    !** End of interface *****************************************
    real(kind=r8_kind)             :: arg

    select case (level_broadening)
    case (gauss_broadening)
       arg = (ef-ei)/(sqrt2*get_sigma())
       dlfunc = one/(sqpi2*get_sigma())*exp(-arg**2)
    case (fermi_broadening)
       arg = get_sigma()*(ef-ei)
       dlfunc = get_sigma()*exp(-arg)/(one+exp(-arg))**2
    case (sinus_broadening)
       dlfunc=  half*get_sigma()*cos(get_sigma()*(ei-ef))
    end select
  end function dlfunc

  !*************************************************************

  function entropy_contribution()
    ! Purpose : calculates the entropy contribution to the
    !           free energy functional
    !           F[rho] = E_tot[rho] - T*S[occ_num]
    !           The calculated contribution is
    !           +T*S[occ_num] such that it must be substracted
    !           from the total energy to yield the free energy.
    !           Depending on the chosen method the contribution
    !           is:
    !           'fermi': -1/sigma*sum_{orbitals}(ni*ln(ni)-(1-ni)ln(1-ni))
    !-------------------------------------------------------------
    real(kind=r8_kind)  :: entropy_contribution
    !** End of interface *****************************************
    real(kind=r8_kind)  :: arg, summ
    integer(kind=i4_kind) :: is,i,mu,degen


    select case (level_broadening)
    case (gauss_broadening)
       entropy_contribution=zero
    case (fermi_broadening)
       summ=zero
       do is=1,symmetry_data_n_spin()
          do i=1,n_irrep
             degen = partner(i) * num_spin
             do mu=1,dim(i)
                arg = occ_num(i)%m(mu,is)/degen
                summ = summ + (x_lnx(arg) + &
                     x_lnx(one-arg))*degen
             enddo

          enddo
       enddo
       entropy_contribution = -summ/get_sigma()
    case (sinus_broadening)
       summ = zero
       do is=1,symmetry_data_n_spin()
          do i=1,n_irrep
             degen = partner(i) * num_spin
             do mu=1,dim(i)
                arg = one - two*occ_num(i)%m(mu,is)/degen
                if (abs(arg)==-one) then
                   summ=summ+zero
                else
                   summ = summ + ((arg*asin(arg) + sqrt(one-arg**2))*two - pi)*degen/four
                endif
             enddo
          enddo
       enddo
       entropy_contribution = -summ/get_sigma()
    case default
       entropy_contribution=zero
    end select
  contains
    function x_lnx(arg)
      ! Purpose: calculates log(arg). If 'arg' is lower than
      !          1.0E-8, a linear curve is used to fit x*ln(x)
      !----------------------------------------------------
      real(kind=r8_kind) :: arg,x_lnx
      real(kind=r8_kind) :: x0,m

      x0 = 1.0e-8_r8_kind
      if (arg == zero ) then
         x_lnx = zero
      elseif ( arg <= x0 ) then
         m = log(x0)
         x_lnx = m*arg
      else
         x_lnx = arg*log(arg)
      endif

    end function x_lnx

  end function entropy_contribution

  !*************************************************************

  function fermi_get_entropy()
    ! Purpose: return the private variable 'entropy'.
    !------------------------------------------------
    real(kind=r8_kind)   :: fermi_get_entropy
    !** End of interface *****************************************
    fermi_get_entropy=entropy
  end function fermi_get_entropy

  !*************************************************************
  subroutine fermi_hole(irrep,orb,is,step,hole,delta)
    ! Purpose: look up if the current orbital is assigned a hole.
    !          This is done by looping through all holes in the
    !          current Irrep (holes_per_irrep(i_irr)%m and see if
    !          if the orbital with indices 'irrep','orb' and 'is'
    !          is a hole.
    !          'delta' returns the amount by which the degeneracy
    !          in the calling routine has to be decreased to achieve the
    !          user-specified occupation number for the hole.
    !
    ! Subroutine called by : newton,newton_spin
    ! --- Declaration of formal parameters --------------------
    integer(kind=i4_kind),intent(in)  :: irrep,orb,is,step
    logical,intent(out)               :: hole
    real(kind=r8_kind),intent(out)    :: delta
    ! --- declaration of local variables ----------------------
    integer(kind=i4_kind)  :: i_hole,degen
    hole=.false.
    degen = symmetry_data_n_partners(irrep) * num_spin
    holes: do i_hole=1,n_holes_per_irrep(irrep)
       if (holes_per_irrep(irrep)%m(i_hole).eq.orb .and. &
            spin_of_hole(irrep)%m(i_hole).eq.is) then
          if (output_fermi) then
             write(output_unit,'("Newton-Step ",i3,":   Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                  step,orb,irrep,eigval(irrep)%m(orb,is)*convert1
             write(stdout_unit,'("Newton-Step ",i3,":   Orbital ",i4," in Irrep ",i2," with eigenvalue ",ES13.3," is a hole")')&
                  step,orb,irrep,eigval(irrep)%m(orb,is)*convert1
          endif
          hole=.true.
          delta=degen-occnum_hole(irrep)%m(i_hole)
          exit holes
       endif
    enddo holes
  end subroutine fermi_hole
  !*************************************************************

![TS:
  !*************************************************************
  subroutine check_occ_fermi()
    ! Purpose: check, if the input switches fermi_fix_up_and_down
    ! and fermi_unpaired have been set by the user (i. e. changed
    ! their default values). If so, then check if fixed_spin_diff
    ! and magn_moment (both namelist occupation) still have their
    ! default values. In this case copy the values from
    ! fermi_fix_up_and_down to fixed_spin_diff  and
    ! fermi_unpaired to magn_moment.
    ! In every other case ignore the values of the two fermi_...
    ! switches
    !
    ! Function called by: read_input
    ! Author: TS
    ! Date  : 12/99
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    logical   :: fixed_spdiff,df_fixed_spdiff
    real(kind=r8_kind) :: magn_momt,df_magn_momt
    !------------ Executable code --------------------------------

    ! make the input switches fixed_spin_diff and
    ! magn_moment (occuptaion_module) available
    call get_fix_spin_switch(fixed_spdiff)
    call get_df_fix_spsw(df_fixed_spdiff)
    call get_spin_diff(magn_momt)
    call get_df_spin_diff(df_magn_momt)

    if (fermi_fix_up_and_down .neqv. df_fermi_up_down) then
       ! old switch 'fermi_fix_up_and_down' was set by user
       if (fixed_spdiff .eqv. df_fixed_spdiff) then
          ! the new switch 'fixed_spin_diff' was ignored by
          ! an ignorant user. Copy value:
          call put_fix_spin_switch(fermi_fix_up_and_down)
       else
          ! both old and new  switch were changed
          ! print out warning
          write(output_unit,*)'Warning: value of FERMI_FIX_UP_AND_DOWN will be ignored;'
          write(output_unit,*)'         instead, value of FIXED_SPIN_DIFF is used.'
       endif
    else
       ! old switch wasn`t used; evertyhing ok
    endif

    if (fermi_unpaired .ne. df_fermi_unpaired) then
       ! old switch 'fermi_unpaired' was set by user
       if (magn_momt .eq. df_magn_momt) then
          ! the new switch 'magn_moment' was ignored by
          ! an ignorant user. Copy value:
          call put_magn_moment(real(fermi_unpaired,kind=r8_kind))
       else
          ! both old and new  switch were changed
          ! print out warning
          write(output_unit,*)'Warning: value of FERMI_UNPAIRED will be ignored;'
          write(output_unit,*)'         instead, value of MAGN_MOMENT is used.'
       endif
    else
       ! old switch wasn`t used; evertyhing ok
    endif

  end subroutine check_occ_fermi
  !*************************************************************

  !--------------- End of module ----------------------------------
end module fermi_module

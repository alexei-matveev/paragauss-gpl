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
module diis_fock_module
  !---------------------------------------------------------------
  !
  !  Purpose:  direct inversion of  iterative subspace  routine (DIIS)
  !           which  accelerates the convergence  of the  SCF-cycle of
  !           the parameter Fock matrix  there will be some older sets
  !           kept and  a linear combination of them  will replace the
  !           actual  Fock  matrix  the  coeffiecents for  the  linear
  !           combination will be  gotten through an minimalization of
  !           the norm of an error vector with  Bij = ei ej for Bc = 0
  !           under the restriction sum ci  = 1; error vector given by
  !           e = HDS - SDH transformed to orthogonal basis.
  !
  !  SET THE RIGHT PARAMETERS:
  !
  !  The default values are set  for Perturbation Theory = true and no
  !  FERMI: they are gotten mostly  by looking on small systems if you
  !  set Perturbation  Theory false set  everystep = 0  for additional
  !  FERMI set the endthreshold (for example 10E-10)
  !
  !  For big systems:  make threshold bigger if DIIS  starts too late,
  !    make it smaller if it  prevents convergence, it only works near
  !    the solution If  in your system all but  the COEFF_CRITERION is
  !    fullfiled  easily and then  the system  needs long  to fullfill
  !    these, set  the endthreshold.  It will make  sure, that  in the
  !    region where there  is no more advancement, DIIS  stops. So the
  !    other convergence methods get a chance. Charge mixing is better
  !    in getting the COEFF_CRITERION low. Set endthreshold not to big
  !    or it  will change the method  too early and loose  too much of
  !    its good  convergence for the other coefficents,  but make sure
  !    it starts somehow. So far 10E-10 was a good value for it.
  !
  !  For very  small system:  if the expected  loops without  DIIS are
  !    below 20, DIIS will not help much, so you do best without it.
  !
  !  To  learn more  about the  parameter  look below,  when they  are
  !  defined.
  !
  !  Module called by: main_scf, mixing_module
  !
  !
  !  References: P. Pulay, Chem. Phys.Lett. 73 (1979),393
  !              P. Pulay, J. Comp. Chem. 3 (1982), 556
  !
  !
  !  Author: AN
  !  Date: 05.05.09
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification include stop_pert
  ! Author: AN
  ! Date:   7.2009
  ! Description:  DIIS is  now able  to stop  the  perturbation theory
  !              calculation when it starts  to work it sets therefore
  !              the perturbation theory flag  to false when it starts
  !              to set a new hamiltonian
  !
  !              ATTENTION: the flag will  be only reset on the master
  !              therefore  the slaves  now have  to get  their values
  !              everytime they need it from the master
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use iounitadmin_module
  use debug, only: octave
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !------------ public functions and subroutines ------------------
  public :: diis_read_input, diis_write_input

  public :: diis_fock_module_close!(), cleans up the module state

  public :: diis_fock_matrix!(ham_tot, overlap, densmat, loop)
  public :: fixed_fock_matrix!(ham_tot, loop)
  public :: diis_charge_coeff!(coeff_charge, coeff_charge_old, metric, beta) -> coeff_charge
  public :: diis_pertstop

  logical, public, protected :: diis_on       ! set on read input, changed occasionally
  logical, public, protected :: diis_step     ! flipping regularly
  logical, public, protected :: diis_charge_mixing ! for coeff_charge, set once on read input

  !
  ! Self-mixing of Fock matrix in the last DIIS step
  ! (it makes no sense to use it before diis_fock_matrix() is called):
  !
  real(r8_kind), public, protected :: diis_fock_beta = HUGE(1.0D0)

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Default values for input parameters -------------
  integer(i4_kind), parameter       :: df_mmax = 8
  integer(i4_kind), parameter       :: df_start_steps = 1
  integer(i4_kind), parameter       :: df_everystep = 0
  integer(i4_kind), parameter       :: df_loop_start = 0
  integer(i4_kind), parameter       :: df_cc_wait = 0
  real(r8_kind), parameter          :: df_cfix = 0.1_r8_kind
  real(r8_kind), parameter          :: df_threshold = 0.1_r8_kind
  real(r8_kind), parameter          :: df_diagdamp = 0.02_r8_kind
  real(r8_kind), parameter          :: df_endthreshold = 0.0E-10_r8_kind
  logical, parameter                :: df_diis_on = .false.
  logical, parameter                :: df_stopchar = .true.
  logical, parameter                :: df_cc_diismixing = .false.
  logical, parameter                :: df_stop_pert = .true.
  !------------ diis input parameters -------------------------
  integer(i4_kind)                  :: mmax
  integer(i4_kind)                  :: start_steps
  integer(i4_kind)                  :: everystep
  integer(i4_kind)                  :: loop_start
  integer(i4_kind)                  :: cc_wait
  real(r8_kind)                     :: threshold
  real(r8_kind)                     :: cfix
  real(r8_kind)                     :: diagdamp
  real(r8_kind)                     :: endthreshold
  real(r8_kind)                     :: ermaxold
  real(r8_kind)                     :: diagdamp2
  logical                           :: diis_nmix = .false. ! kept for compatibility
  logical                           :: stopchar
  logical                           :: stop_pert
  logical                           :: cc_diismixing ! == diis_charge_mixing

  !
  ! INPUT PARAMETER
  !
  ! diis_on: if true diis is performed, if false nothing is done
  !
  ! mmax: is maximum number of stored old hamiltonians
  !
  ! mfix: is the number of initial fixed mixing steps for hamiltonian
  !
  ! cfix: is the corresponding fixed mixing coefficient
  !
  ! start_steps: when there are more than start_steps old hamiltonians
  !   stored, mixing starts
  !
  ! everystep: if > 0 there will be only the next diis_step after
  !   everystep steps of relaxation
  !
  ! threshold:> maxval(errorvec), then hamiltonians will be stored
  !   after the routine fell back over that value, the storage starts
  !   by 1 again
  !
  ! diagdamp: damping value, compare the end of the routine
  !   diis_update
  !
  ! diis_nmix: ignored, kept for input compatibility. FIXME: remove it
  !
  ! stopchar: if true dynamical mixing (mixing_module) will be stopped
  !   when a diis step is performed
  !
  ! endthreshold: if < maxval(errovec) DIIS will stop if the previous
  !   step does not improve the errorvector
  !

  namelist /diis/ mmax, &
                  cfix, &
                  start_steps, &
                  everystep, &
                  threshold, &
                  diagdamp, &
                  diis_on, &
                  diis_nmix, &
                  loop_start, &
                  endthreshold, &
                  cc_wait, &
                  cc_diismixing, &
                  stop_pert, &
                  stopchar

  integer(kind=i4_kind)           :: diis_loop ! counts only loops, when diis actuals his storage
  real(kind=r8_kind),allocatable       :: B(:,:) ! mmax x mmax
  type(arrmat3),allocatable            :: ermat(:,:), hamold(:,:) ! mmax x n_irrep
  type(arrmat2),allocatable         :: Xortho(:) ! n_irrep
  logical                              :: diis_start, diis_onlyonce
  logical                              :: onlyonce, pert_store
  real(kind=r8_kind)                :: cc_thres
  ! for the diismixing routine
  real(kind=r8_kind),allocatable       :: b_diis(:,:)
  real(kind=r8_kind),allocatable :: coeff_charge_diisdp(:,:), coeff_charge_diis(:,:)
  integer(kind=i4_kind)       :: diis_ndim, cc_estep, cc_control, cc_danf
  integer(kind=i4_kind)           :: cc_diis_loop ! diis_loop for diismixing

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !**************************diis_fock********************************

  subroutine diis_fock_matrix(hamact, overlap, densmatact, loop_scf)
    !
    ! Purpose: replaces the actual fockmatrix by a linear combination
    !          of older Fock matrices via the DIIS algorithm
    !
    ! This sub is supposed to set
    !
    !   diis_fock_beta
    !
    ! to the self-mixing proportion.
    !
    ! FIXME: This global public variable is used elsewhere at another
    !        times, therefore it may appear inconvenient to pass it to
    !        outside via the intent(out) argument.
    !
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: loop_scf
    type(arrmat3),         intent(inout) :: hamact(:)
    type(arrmat3),         intent(in) :: densmatact(:)
    type(arrmat2),         intent(in) :: overlap(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: diis_allocstat, i, diis_loop_m, n_irrep, n_dim, n_spin, actmatsize
    real(kind=r8_kind)                   :: c(mmax)
    type(arrmat3), allocatable           :: ermatact(:) !n_irrep
    !------------ Executable code --------------------------------

    if(loop_scf > loop_start) then
      print*, "DIIS-Fock matrix: entered", loop_scf
      if(.not. allocated(B)) call diis_initalise(overlap, hamact)

      ! allocate the actual error vector
      n_irrep = size(hamact)
      allocate(ermatact(n_irrep), stat = diis_allocstat)
      ASSERT(diis_allocstat .eq. 0)
      do i = 1, n_irrep
        n_dim = size(hamact(i)%m,1)
        n_spin = size(hamact(i)%m,3)
        allocate(ermatact(i)%m(n_dim,n_dim, n_spin), stat = diis_allocstat)
        ASSERT(diis_allocstat .eq. 0)
      enddo

      ! calculates the actual error vector
      call calc_act_errmat(hamact, overlap, densmatact, ermatact)
      ! decides whether and which part of the diis-routine should run
      ! ATTENTION: in the first iteration the convergence measure will be 0, thus DIIS
      !           would start to collect results. Thus skip that iteration.
      if (loop_scf > 1) call diis_control(ermatact, diis_start, diis_step, diis_loop, start_steps)

      ! updates the matrix B and the remembered error matrices
      if(diis_start) call diis_update(ermatact, ermat, B, hamold, hamact, diis_loop, diis_loop_m, diis_step)

      ! makes the diis-step
      if (diis_step) then
        actmatsize = min(diis_loop,mmax)
        c(:mmax) = 0.0_r8_kind
        c(:actmatsize) = diis_calc_c(B(:actmatsize, :actmatsize), diagdamp)
        call write_to_output_units('coefficients for new matrix:')
        write(output_unit,*) c
        print*, c
! call octave("c", c)

        ! mixing of hamiltonian here:
        call diis_update_ham( hamact, hamold, c )

        !
        ! This data is used outside of the module, set public global var ...
        !
        diis_fock_beta = c(diis_loop_m)
      else
        !
        ! No-mixing is equivalent to beta=1.0 (learn to love branches)
        !
        diis_fock_beta = 1.0
      endif

      ! the actual error vector is not needed for this loop any more
      do i = 1, n_irrep
        deallocate(ermatact(i)%m, stat = diis_allocstat)
        ASSERT(diis_allocstat .eq. 0)
      enddo
      deallocate(ermatact, stat = diis_allocstat)
      ASSERT(diis_allocstat .eq. 0)
    else
      diis_step = .false.
      !
      ! No-mixing is equivalent to beta=1.0 (learn to love branches)
      !
      diis_fock_beta = 1.0
    endif

    if(.not. stopchar) diis_step = .false.
  end subroutine diis_fock_matrix

  !*************************************************************

  subroutine fixed_fock_matrix(hamact, loop_scf)
    !
    ! Purpose: replaces the actual Fock matrix by a linear combination
    ! of older Fock matrices via a fixed mixing ratio.
    !
    implicit none
    integer(kind=i4_kind), intent(in)    :: loop_scf
    type(arrmat3),         intent(inout) :: hamact(:)
    ! *** end of interface ***

    !
    ! FIXME:  beware  of the  two  SAVEd  variables preserved  between
    ! invokations.  Saved data  is  initialized (for  integers ---  to
    ! zero) only ONCE at the program startup.
    !
    type(arrmat3), allocatable, save :: hamlast(:)
    integer(i4_kind), save :: was_empty ! 0 when HAMLAST is not in
                                        ! use; 1 when HAMLAST is
                                        ! allocated; > 1 when HAMLAST
                                        ! is filled with data.
    integer(i4_kind) :: i_gamma

    if (.not. diis_step .and. cfix /= 0.0) then
       !
       ! Fixed     mixing    for     initial     steps.     Subroutine
       ! fixed_update_ham() used below  replaces HAMACT and HAMLAST by
       ! a  linear combination  of the  two  with the  mixing ratio  C
       ! depending on the value of "was_empty":
       !
       !    HAMACT := C * HAMACT + (1 - C) * HAMLAST
       !    HAMLAST := HAMACT
       !
       ! Note that "was_empty" is  NOT an "iteration count" but ranges
       ! from 0 to 3.
       !
       if (was_empty > 1) then
          !
          ! Perform all other fixed mixing steps with cfix
          !
          ASSERT(allocated(hamlast))
          call fixed_update_ham(cfix, hamact, hamlast)
       else
          !
          ! Allocate and perform first fixed mixing steps with cfix = 1
          !
          ASSERT(was_empty==0)

          !
          ! Allocate array to hold the last hamiltonian
          !
          allocate(hamlast(size(hamact)))
          do i_gamma = 1, size(hamact)
             allocate(hamlast(i_gamma)%m(size(hamact(i_gamma)%m, 1) &
                                       , size(hamact(i_gamma)%m, 2) &
                                       , size(hamact(i_gamma)%m, 3)))
          enddo
          was_empty = 2

          !
          ! First iteration is different, just save HAMACT, there is no
          ! meaningfull last Fock matrix:
          !
          do i_gamma = 1, size(hamact)
             hamlast(i_gamma)%m = hamact(i_gamma)%m
          enddo
      endif
    endif

    if (diis_step) then
       !
       ! Not needed any more (if everything converges fine). Clean up:
       !
       if (allocated(hamlast)) then
          do i_gamma = 1, size(hamlast)
             deallocate(hamlast(i_gamma)%m)
          enddo
          deallocate(hamlast)

          ! Reset this for consistency
          was_empty = 0
       endif
    endif
  end subroutine fixed_fock_matrix

  !*************************************************************

  subroutine calc_act_errmat(ham, overlap, dens, ermatact)
    ! purpose: new errormatrix according to e = S^-1/2 (FDS - k.c) S^-1/2
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: matmult
    use symmetry_data_module, only: symmetry_data_n_partners
    implicit none
    type(arrmat3),         intent(in) :: ham(:), dens(:)
    type(arrmat2),         intent(in) :: overlap(:)
    type(arrmat3),         intent(inout) :: ermatact(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)          :: i_gamma, sval, n_irrep, n_spin, n_partner
  !------------ Executable code --------------------------------
    ! first: e = (FDS - k.c.)
    n_irrep = size(ham)
    n_spin = size(ham(1)%m,3)
    do i_gamma = 1, n_irrep
      n_partner = symmetry_data_n_partners(i_gamma)
      do sval = 1, n_spin
        ! n_partner is used, because of the density
        ermatact(i_gamma)%m(:,:,sval) = matmult(ham(i_gamma)%m(:,:,sval), dens(i_gamma)%m(:,:,sval))/ n_partner
        ermatact(i_gamma)%m(:,:,sval) = matmult(ermatact(i_gamma)%m(:,:,sval), overlap(i_gamma)%m(:,:))
        ermatact(i_gamma)%m(:,:,sval) = ermatact(i_gamma)%m(:,:,sval) - transpose(ermatact(i_gamma)%m(:,:,sval))
      enddo
    enddo
! call octave("D", dens(1)%m(:,:,1))
! call octave("Error", ermatact(1)%m(:,:,1))

    ! convert to orthonormal basis with Xortho = S^-1/2
    do i_gamma = 1, n_irrep
      do sval = 1, n_spin
        ermatact(i_gamma)%m(:,:,sval) = matmult(ermatact(i_gamma)%m(:,:,sval), Xortho(i_gamma)%m(:,:))
        ermatact(i_gamma)%m(:,:,sval) = matmult(Xortho(i_gamma)%m(:,:), ermatact(i_gamma)%m(:,:,sval))
      enddo
    enddo
! call octave("Errorort", ermatact(1)%m(:,:,1))
  end subroutine calc_act_errmat

  !*************************************************************

  subroutine diis_control(ermatact, diis_start, diis_step, diis_loop, start_stepsin)
    ! purpose: decides what of the diis- cycle will be done:
    !          diis_start true means that the rememberd errormatrices and the
    !          residium matrix will be updated with the new error matrix
    !          diis_step true means that the diis-step will be performed and the
    !          parameter (FOCK-matrix) will be updated
    use options_module, only: options_switch_pertth
    implicit none
    type(arrmat3),         intent(in) :: ermatact(:)
    integer(kind=i4_kind),    intent(in) :: start_stepsin
    logical,               intent(inout) :: diis_start, diis_step
    integer(kind=i4_kind) , intent(inout)            :: diis_loop
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: i_gamma, sval, n_irrep, n_spin
    real(kind=r8_kind)                :: ermax
  !------------ Executable code --------------------------------
    !the maximum value of the errormatrix is needed
    n_irrep = size(ermatact)
    n_spin = size(ermatact(1)%m,3)
    ermax = 0.0_r8_kind
    do i_gamma =1, n_irrep
      do sval = 1, n_spin
         ermax = max(ermax, maxval(abs(ermatact(i_gamma)%m(:,:,sval))))
      enddo
    enddo

      ! stop_pert true means no perturbation_theory calculations
      ! when diis updates the hessian
      ! if diis does not update the hessian this iteration but has done so
      ! the last one it has to reset the perturbation_theory flag to its old value
      ! later it will decide if it stops it for this calculation
      if(stop_pert) then
        call options_switch_pertth( pert_store)
      endif

    ! the maximum value of the errormatrix should be below a threshold
    ! diis_loop only counts the loops when values for diis are stored, else its reset to zero
    if (ermax < threshold) then
      diis_start = .true.
      diis_loop = diis_loop + 1
    else if (diis_start .and. ermax < threshold * 10) then
      ! If DIIS has started do only stop it again after it got much worse.
!     print *, "keep DIIS", ermax
      diis_loop = diis_loop + 1
    else
      diis_start = .false.
      diis_loop = 0
    endif
    ! steps will performed when matrices were updated and when there have been more than start_steps
    ! steps and additional only after everystep cycles of relaxation
    if(diis_start .and. (diis_loop > start_stepsin) .and. mod((diis_loop -start_stepsin -1),&
             (everystep +1)) .eq. 0) then
      diis_step = .true.
      ! stop_pert true means no perturbation_theory calculations
      ! when diis updates the hessian
      ! so set the its flag to false right now
      if(stop_pert) then
        call options_switch_pertth(.false.)
      endif
    else
      diis_step = .false.
    endif
    write(output_unit,*) 'DIIS - loop ', diis_loop, 'has maximum value', ermax
    print *, 'DIIS - loop ', diis_loop, 'has maximum value', ermax

    ! tell output about control result
    if(diis_start) then
      call  write_to_output_units( 'max(error) < threshold, therefore diis runs')
      if(diis_onlyonce) then
        call write_to_trace_unit("DIIS started, mixing after steps: ", start_stepsin)
        diis_onlyonce = .false.
      endif
    endif
    if(diis_start .and. (.not. diis_step)) then
      call write_to_output_units( 'no diis Fock matrix update because:')
      if(diis_loop .le. start_stepsin) then
        call write_to_output_units( 'not enough values for diis-mixing')
      elseif(.not. (mod((diis_loop -start_stepsin -1), (everystep +1)) .eq. 0)) then
        call write_to_output_units( 'relaxation of system')
        call write_to_trace_unit("DIIS: relaxation of system")
      endif
    endif
    if(.not. diis_start .and. .not. diis_onlyonce) then
      diis_onlyonce = .true.
      call write_to_trace_unit("DIIS stopped")
    endif

    if( ermax < endthreshold .and. ermax > ermaxold) then
      diis_on = .false.
      diis_step = .false.
      ! stop_pert true means no perturbation_theory calculations
      ! when diis updates the hessian
      ! if diis stops the perturbation_theory flag can get its old value back
      if(stop_pert) call options_switch_pertth( pert_store)
      call write_to_output_units( 'DIIS converged, choose other convergence method')
      call write_to_trace_unit("DIIS: stopped, endconvergence reached")
      if ( diis_charge_mixing ) then
        call write_to_output_units( 'new method: charge coefff. DIIS')
      call write_to_trace_unit("cc_DIIS: started instead")
      else
        call write_to_output_units( 'return to starting method')
      endif
    endif
    ermaxold = ermax

    if ( diis_charge_mixing ) then
      ! decide for diismix what to do
      if(cc_control < 4) then
        ! this is the case where diis_fock_matrix starts and diismix takes over
        ! when diis_fock_matrix ends
        if(ermax < endthreshold .and. diis_loop > 2) then
          if(diis_on) then
            ! still diis_fock_matrix used but store values for diismix
            cc_control = 2
          else
            ! diismix starts now and replaces diis_fock_matrix
            cc_control = 3
          endif
        else
          cc_control = 1
        endif
      else
      ! in this case the two diis routines are performed alternatly
      ! everystep gives amount of diismix loops after one from the diis_fock_matrix
      cc_control = 4
        if(diis_start) cc_control = cc_control + 1
        if(diis_start .and. (.not. diis_step) .and. &
           (diis_loop > start_stepsin)) cc_control = cc_control + 1
      endif
    endif

  end subroutine diis_control

  !*************************************************************

  subroutine diis_update(ermatact, ermat, B, hamold, hamact, diis_loop, diis_loop_m, diis_step)
    ! purpose: ermatact is put in ermat, ermat is filled in B(residiuum matrix)
    !------------ Modules used ------------------- ---------------
    use symmetry_data_module, only: symmetry_data_n_partners
    implicit none
    type(arrmat3),         intent(inout) :: ermat(:,:), hamold(:,:)
    type(arrmat3),         intent(in) :: ermatact(:)
    real(kind=r8_kind),    intent(inout) :: B(:,:)
    type(arrmat3),         intent(in) :: hamact(:)
    integer(kind=i4_kind), intent(inout) :: diis_loop, diis_loop_m
    logical,               intent(inout) :: diis_step
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: i_gamma, sval, j, n_irrep, n_spin, n_partner
    real(kind=r8_kind)                :: bwork
  !------------ Executable code --------------------------------
    n_irrep = size(ermatact)
    n_spin = size(ermatact(1)%m,3)
    ! diis_loop_m is the working storage number
    ! diis_loop starts with 1, diis_loop = n* mmax -> diis_loop_m = mmax
    diis_loop_m = mod(diis_loop -1,mmax) +1

    ! the actual error matrix is put in the error matrix storage, same for hamiltonian
    do i_gamma =1, n_irrep
      do sval = 1, n_spin
        ermat(diis_loop_m, i_gamma)%m(:,:, sval) = ermatact(i_gamma)%m(:,:, sval)
        hamold(diis_loop_m, i_gamma)%m(:,:, sval) = hamact(i_gamma)%m(:,:, sval)
      enddo
    enddo

! do j =1, mmax
!   call octave("Ham", hamold(j,1)%m(:,:,1))
! enddo
    ! calculate new values for B
    do j = 1, min(diis_loop,mmax)
     ! if B is not yet fully occupied there need only a part of the row to be replaced/set
     ! insertion starts on the left top of B and goes (jmax) steps down
     bwork = 0.0_r8_kind
     ! B(i,j) = trace(e_i, e_j^T)
     do i_gamma =1, n_irrep
       ! some irreps would appear more than once in the full hamiltonian
       n_partner = symmetry_data_n_partners(i_gamma)
       do sval = 1, n_spin
         bwork = bwork + sum( ermat(j, i_gamma)%m(:,:,sval) *ermatact(i_gamma)%m(:,:,sval))*n_partner
       enddo
     enddo
     B(j, diis_loop_m) = bwork
     B(diis_loop_m, j) = bwork
    enddo

    ! B(i,i) = 0 means, error(i) = 0, H(i) best next trial vector
    if(B(diis_loop_m, diis_loop_m) .eq. 0) then
      if(diis_step) call write_to_output_units( 'in this loop no diis mixing, because error = 0')
      diis_step = .false.
      diis_loop = diis_loop -1
    endif
! call octave("B",B)
  end subroutine diis_update

  !*************************************************************

  function diis_calc_c(B, diagdampin) result(c)
    ! purpose: solves for the new coefficients c (Bc = rhs)
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: linsolve
    implicit none
    real(kind=r8_kind),    intent(in   ) :: B(:,:), diagdampin
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)                   :: ccal(size(B,1)+1), rhs(size(B,1)+1)
    real(kind=r8_kind)                   :: Bcal(size(B,1)+1,size(B,2)+1), c(size(B,1))
    integer(kind=i4_kind)                :: sb, i
    !------------ Executable code --------------------------------
    sb = size(B,1)
    ! first the right hand side of the equation is build
    rhs(1:sb) = 0.0_r8_kind
    rhs(sb+1) = -1.0_r8_kind

    ! the matrix is extended
    Bcal(:sb, :sb) = B
    Bcal(:sb, sb +1) = -1.0_r8_kind
    Bcal( sb +1, :sb) = -1.0_r8_kind
    Bcal(sb +1, sb +1 ) = 0.0_r8_kind

    ! something like a damping factor should improve convergence for nearly linear behaviour
    ! of the error vector
    ! reference:   T. P.  Hamilton, P.Pulay J. Chem. Phys. 84, 5728 (1986); DOI:10.1063/1.449880
    do i = 1, sb
      Bcal(i,i) = (1.0_r8_kind + diagdampin) *Bcal(i,i)
    enddo

! call octave("Bcal", Bcal)
    ccal = linsolve(Bcal, rhs)
    c = ccal(1:sb)
! call octave("ccal", ccal)
  end function diis_calc_c

  !*************************************************************

  subroutine fixed_update_ham(c, hamact, hamold)
    !
    ! Purpose: gives new hamiltonian as H = sum(H_i c_i).
    !
    implicit none
    real(kind=r8_kind), intent(in)    :: c
    type(arrmat3),      intent(inout) :: hamact(:)
    type(arrmat3),      intent(inout) :: hamold(:)
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i_gamma

    do i_gamma = 1, size(hamact)
      !
      ! Combine actual with last hamiltonian:
      !
      hamact(i_gamma)%m   =                c  * hamact(i_gamma)%m              &
                          + (1.0_r8_kind - c) * hamold(i_gamma)%m
      !
      ! Update last hamiltonian:
      !
      hamold(i_gamma)%m = hamact(i_gamma)%m
    enddo
  end subroutine fixed_update_ham

  !*************************************************************

  subroutine diis_update_ham( hamact, hamold, c )
    ! purpose: gives new hamiltonian as H = sum(H_i c_i)
    implicit none
    type(arrmat3),         intent(inout) :: hamact(:)
    type(arrmat3),         intent(in) :: hamold(:,:)
    real(kind=r8_kind),        intent(in) :: c(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: i_gamma, j, n_irrep, n_spin, n_dim
    !------------ Executable code --------------------------------
! call octave("Hold",hamact(1)%m(:,:,1))
    n_irrep = size(hamact)
    n_spin = size(hamact(1)%m,3)
    do i_gamma = 1, n_irrep
      n_dim = size(hamact(i_gamma)%m,1)
      hamact(i_gamma)%m(1:n_dim, 1:n_dim, 1:n_spin) = 0
      do j = 1, min(diis_loop,mmax)
        ! add all hamolds
        hamact(i_gamma)%m(1:n_dim, 1:n_dim, 1:n_spin) = hamact(i_gamma)%m(1:n_dim, 1:n_dim, 1:n_spin) &
        +  c(j)*hamold(j,i_gamma)%m(1:n_dim, 1:n_dim, 1:n_spin)
      enddo
    enddo
! call octave("Hnew",hamact(1)%m(:,:,1))
    print*, "DIIS-Fock matrix: Fock-matrix updated"
  end subroutine diis_update_ham

  !****************************start read and end***********************************

  subroutine diis_initalise(overlap, hamact)
    ! purpose: initialisation of diis, allocates needed stuff, starts loop, makes S^-1/2 for changing basis
    !          needed; overlap S for making transfer matrix, hamact for getting size of system
    !------------ Modules used ------------------- ---------------
    use matrix_functions, only: funm
    use options_module, only: options_perturbation_theory
    implicit none
    type(arrmat2),         intent(in) :: overlap(:)
    type(arrmat3),         intent(in) :: hamact(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)              :: allocate_state, i_gamma, j, n, n_irrep, n_spin
  !------------ Executable code --------------------------------
    print*, "DIIS-Fock matrix: start allocation"
    n_irrep = size(hamact)
    n_spin = size(hamact(1)%m,3)
    diis_loop = 0
    cc_diis_loop = 0
    ! allocates ermat (storage for error matrix), B (residiuum matrix of error matrices)
    ! hamold (storage of old hamiltonians), Xortho for transforming basis
    allocate(B(mmax, mmax), &
        ermat( mmax,n_irrep), hamold( max(mmax,1),n_irrep),&
        Xortho(n_irrep), stat = allocate_state)
      ASSERT(allocate_state .eq. 0)
      B(:,:) = 0.0_r8_kind
    do i_gamma = 1, n_irrep
      n = size(hamact(i_gamma)%m,1)
      do j = 1, mmax
        allocate( ermat(j,i_gamma)%m(n, n, n_spin), stat = allocate_state )
        ASSERT(allocate_state .eq. 0)
        ermat(j,i_gamma)%m(:,:,:) = 0.0_r8_kind
      enddo
      do j = 1, mmax
        allocate( hamold(j,i_gamma)%m(n, n, n_spin), stat = allocate_state)
        ASSERT(allocate_state .eq. 0)
        hamold(j,i_gamma)%m(:,:,:) = 0.0_r8_kind
      enddo
      allocate(Xortho(i_gamma)%m(n, n), stat = allocate_state)
      ASSERT(allocate_state .eq. 0)
    enddo

    ! Xortho = S^-1/2 is build of overlap matrix with lapack routine
    print*, "DIIS-Fock matrix: build transfer matrix"
    do i_gamma = 1, n_irrep
      n = size(hamact(i_gamma)%m,1)
      Xortho(i_gamma)%m(:, :) = funm(overlap(i_gamma)%m(:, :), invsqrt)
    enddo
! call octave("S", overlap(1)%m)
! call octave("X", Xortho(1)%m)
    print*, "DIIS-Fock matrix: initalisation ended"
    diis_onlyonce = .true.
    onlyonce = .true.
    ! diis_onlyonce for diis_fock, onlyonce for diismix
    ! for control when to write in trace-output

    ! stop_pert true means no perturbation_theory calculations
    ! when diis updates the hessian
    ! if diis will stop again it needs to know which value the
    ! perturbation_theory flag had had
    ! take care: only master has the most recent value of the flag
    if(stop_pert) pert_store = options_perturbation_theory()

  end subroutine diis_initalise

  !*************************************************************

  subroutine diis_read_input(scfcontrol)
    ! reads in the input for the diis_fock_module
    !------------ Modules used ----------------------------------
    use input_module, only : input_line_is_namelist                            &
                           , input_read_to_intermediate                        &
                           , input_intermediate_unit                           &
                           , input_error
    implicit none
    logical, intent(in), optional :: scfcontrol
    !** End of interface ***************************************
    integer(kind=i4_kind) :: unit,status
    !------------ Executable code --------------------------------
    mmax = df_mmax
    start_steps = df_start_steps
    everystep = df_everystep
    threshold = df_threshold
    diagdamp = df_diagdamp
    diis_on = df_diis_on
    stopchar = df_stopchar
    loop_start = df_loop_start
    endthreshold = df_endthreshold
    cc_wait = df_cc_wait
    ! this is the name used in the input:
    cc_diismixing = df_cc_diismixing
    stop_pert = df_stop_pert

    ! needed for mixing module, there diis_step = true should prevent charge mixing
    ! when the Fock matrix is mixed by diis; here it is ensured that charge mixing can start
    ! and run properly, without care what diis does or it if is ever started
    diis_step = .false.
    ermaxold = 100.0

    if ( input_line_is_namelist("diis") ) then
          call input_read_to_intermediate
          unit= input_intermediate_unit()
          read(unit, nml=diis, iostat=status)
          if (status .gt. 0) call input_error( &
                "diis: namelist diis")
    endif

    ! this is the name used in the code:
    diis_charge_mixing = cc_diismixing

    if (mmax .ge. 100) then
      call input_error &
            ("DIIS_READ_INPUT: mmax too large")
    elseif (mmax .lt. 1) then
       call input_error &
            ("DIIS_READ_INPUT: mmax too small")
    endif
    if (start_steps .ge. mmax) then
      call input_error &
            ("DIIS_READ_INPUT: start_steps too large")
    elseif (start_steps .lt. 1) then
       call input_error &
            ("DIIS_READ_INPUT: start_steps too small")
    endif
    if (everystep .ge. mmax) then
      call input_error &
            ("DIIS_READ_INPUT: everystep too large")
    elseif (everystep .lt. 0) then
       call input_error &
            ("DIIS_READ_INPUT: everystep too small")
    endif
   if (threshold .ge. 1.0E200_r8_kind) then
      call input_error &
            ("DIIS_READ_INPUT: threshold too large")
    elseif (threshold .lt. 0.0_r8_kind) then
       call input_error &
            ("DIIS_READ_INPUT: threshold too small")
    endif
   if (endthreshold .ge. threshold) then
      call input_error &
            ("DIIS_READ_INPUT: endthreshold too large")
    elseif (endthreshold .lt. 0.0_r8_kind) then
       call input_error &
            ("DIIS_READ_INPUT: endthreshold too small")
    endif
    if (diagdamp .lt. 0.0_r8_kind) then
       call input_error &
            ("DIIS_READ_INPUT: diagdamp too small")
    endif

    if (loop_start .lt. 0) then
       call input_error &
            ("DIIS_READ_INPUT: loop_start negativ")
    endif

    ! for the diismix, if no diis_fock mixing (diis_on) the parameter
    ! belong to diismix, else it will be controled by diis-fock
    if ( diis_charge_mixing ) then
      if(diis_on) then
        cc_estep = cc_wait
        cc_control = 1
        if( cc_wait < 0) cc_control = 4
      else
        cc_control = 0
        cc_thres = threshold
        cc_estep = everystep
      endif
      diagdamp2 = diagdamp
      diis_ndim = mmax
    endif

  end subroutine diis_read_input

  !*************************************************************

  subroutine diis_write_input(iounit,full,scfcontrol)
    ! purpose: write the namelist fock_diis
    !------------ Modules used -----------------------------------
    use echo_input_module, only: start, flag, real, intg, stop                 &
                               , echo_level_full
    use operations_module, only: operations_echo_input_level
    implicit none
    !------------ Declaration of formal parameters -------------
    integer, intent(in)           :: iounit
    logical, intent(in), optional :: full
    logical, intent(in), optional :: scfcontrol
    !** End of interface ***************************************
    !------------ Declaration of subroutines used --------------
    external error_handler
    !
    character(len=35)             :: real_format
    character(len=35)             :: intg_format
    character(len=35)             :: flag_format
    !------------ Executable code ------------------------------
    real_format = '("    ",a," = ",es9.3:" # ",a)'
    intg_format = '("    ",a," = ", i9  :" # ",a)'
    flag_format = '("    ",a," = ",4x,a5:" # ",a)'

    if (present(full)) then
       call start("DIIS","DIIS_WRITE_INPUT", &
            iounit,echo_level_full)
    else
       call start("DIIS","DIIS_WRITE_INPUT", &
            iounit,operations_echo_input_level)

    endif
    call flag("DIIS_ON              ", &
               diis_on, df_diis_on)
    call flag("CC_DIISMIXING        ", &
               cc_diismixing, df_cc_diismixing)
    call flag("STOPCHAR             ", &
               stopchar, df_stopchar)
    call flag("STOP_PERT            ", &
               stop_pert, df_stop_pert)
    call intg("MMAX                 ", &
               mmax        , df_mmax )
    call intg("LOOP_START           ", &
               loop_start             , df_loop_start )
    call intg("START_STEPS          ", &
               start_steps           , df_start_steps )
    call intg("EVERYSTEP            ", &
               everystep             , df_everystep )
    call intg("CC_WAIT              ", &
               cc_wait             , df_cc_wait )
    call real("THRESHOLD            ", &
               threshold    , df_threshold )
    call real("ENDTHRESHOLD         ", &
               endthreshold    , df_endthreshold )
    call real("DIAGDAMP             ", &
               diagdamp       , df_diagdamp )
    call real("CFIX                 ", &
               cfix           , df_cfix )

    if (present(scfcontrol)) then
       if (scfcontrol) then
          call stop(empty_line=.false.)
       else
          call stop()
       end if
    else
       call stop()
    endif

  end subroutine diis_write_input

  !*************************************************************

  subroutine diis_fock_module_close ()
    !
    ! Purpose: deallocates diis matrices, cleans up module state.
    !
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: allocate_state, n_irrep, i_gamma, j

    if(allocated(B)) then
      n_irrep=size(Xortho)
      do i_gamma = 1, n_irrep
        do j = 1, mmax
          deallocate(ermat(j,i_gamma)%m, hamold(j,i_gamma)%m, stat = allocate_state)
          ASSERT(allocate_state .eq. 0)
        enddo
        deallocate(Xortho(i_gamma)%m, stat = allocate_state)
        ASSERT(allocate_state .eq. 0)
      enddo
      deallocate(B, ermat, hamold, Xortho, stat = allocate_state)
      ASSERT(allocate_state .eq. 0)
      ! in the case that endthreshold was set, diis_on was put false during the
      ! scf part, here it is reset for the next one
      diis_on = .true.
    endif
    if(allocated(b_diis)) then
      deallocate(b_diis,coeff_charge_diis, coeff_charge_diisdp, stat = allocate_state)
      ASSERT(allocate_state .eq. 0)
    endif
  end subroutine diis_fock_module_close

  !*************************************************************

  function invsqrt(x) result (f)
    ! function = x^-1/2
    implicit none
    real(kind=r8_kind), intent(in) :: x
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)             :: f
    !------------ Executable code --------------------------------
    f = 1.0_r8_kind / sqrt( x )
  end function invsqrt

  !**************************************diismixing***************

  function diis_charge_coeff( coeff_charge_now, coeff_charge_old &
                                          , metric, beta )  result(coeff_charge)
    !  Purpose: replaces the actual coeff_charge by a linear combination of
    !           older ones via the DIIS algorithm
    use fit_coeff_module, only: ff_map => fit_coeff_ff_map
    implicit none
    real(kind=r8_kind), intent(in)  :: coeff_charge_now(:)
    real(kind=r8_kind), intent(in)  :: coeff_charge_old(:)
    real(kind=r8_kind), intent(in)  :: metric(:)
    real(r8_kind),      intent(out) :: beta
    ! Separate output variable to make the real coeff_charge protected
    real(kind=r8_kind)              :: coeff_charge(size(coeff_charge_now))! out
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: alloc_stat, c_dim, actsize, cc_diis_loop_m, i
    real(kind=r8_kind)                   :: c(diis_ndim)
    real(kind=r8_kind)                   :: coeff_chargdis(size(coeff_charge_now))
    logical                              :: cc_start, cc_step
    !------------ Executable code --------------------------------
    c_dim = size(coeff_charge_now)

    ! in the first loop there need some parameters to be allocated, they will be
    ! disallocated together with the ones from the Fock_diis mixing
    if( .not.allocated(coeff_charge_diis) ) then
      allocate(coeff_charge_diis(c_dim,diis_ndim), stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      allocate(coeff_charge_diisdp(c_dim,diis_ndim), stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      allocate(b_diis(diis_ndim, diis_ndim), stat=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      b_diis(:,:)=0.0_r8_kind
    endif
    ! the new error vector
    coeff_chargdis = coeff_charge_now - coeff_charge_old

    ! here the error vector is needed to decide if diis should start, if it has to work alone
    ! the coul_norm term makes the scale better, so that values with few important will not dominate
    ! the decicion when to start, the decicion how diismix runs could also come from outside the
    ! routine via cc_control > 0
    call cc_diis_control(coeff_chargdis * sqrt(ff_map(:)%COUL_NORM),&
                         cc_start, cc_step, cc_diis_loop, cc_control, cc_thres, cc_estep)

    if(cc_start) then
      cc_diis_loop_m = mod(cc_diis_loop - 1, diis_ndim) + 1
      actsize = min(cc_diis_loop, diis_ndim)
      call write_to_output_units('cc_DIIS: update data')
      ! the b_matrix and the storage vectors are updated
      ! metric is needed for the norm of the error vectors to make b
      call cc_diis_update(b_diis, coeff_charge_diisdp, coeff_chargdis, coeff_charge_diis,&
                           coeff_charge_old, cc_diis_loop_m, actsize, metric, cc_diis_loop, cc_step)
    endif

    if(cc_step) then
      c(:diis_ndim) = 0.0_r8_kind
      ! get the new diis_coeff but only for the used matrix size
      c(:actsize) = diis_calc_c(b_diis(:actsize, :actsize), diagdamp2)
      call write_to_output_units('coefficients for new matrix:')
      write(output_unit,*) c

      ! the new coeff_charge were build, coeff_charge_diis holds the input values
      ! of the coeff_charge, meaning the coeff_charge_old
      coeff_charge = 0.0_r8_kind
      do i = 1, actsize
        coeff_charge = coeff_charge + c(i) * coeff_charge_diis(:, i)
      enddo

      ! if diismix should work in turn with charge mixing
      ! there is the need of beta as old value for charge mixing in the next loop
      beta = c(cc_diis_loop_m)
    else
      !
      ! In this branch coeff_charge remains unmodified, effective beta = 1:
      !
      coeff_charge = coeff_charge_now
      !
      beta = 1.0
    endif
  end function diis_charge_coeff

  !*************************************************************

  subroutine cc_diis_update(b, val_diisdp, valdif, val_diis, val, loop_m, asize, metric, loop, step)
    !Purpose: stores the new value and the new error vector and updates the matrix b (b_diis)
    real(kind=r8_kind), intent(inout) :: val_diisdp(:,:), val_diis(:,:) !size(coeff_charge) x n_diis
    real(kind=r8_kind), intent(inout) :: b(:,:)
    real(kind=r8_kind), intent(in) :: val(:), valdif(:)
    integer(kind=i4_kind), intent(in)       :: loop_m, asize
    integer(kind=i4_kind), intent(inout)       :: loop
    real(kind=r8_kind),intent(in) :: metric(:)
    logical,               intent(inout) :: step
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: j, i, k, ik, c_dim
    real(kind=r8_kind)                :: bwork
    !------------ Executable code --------------------------------
    c_dim = size(valdif)
    val_diisdp(:, loop_m) = valdif
    val_diis(:, loop_m) = val

    do j = 1, asize
      bwork = 0.0_r8_kind

      ! bwork = < valdif | metric | valdif* >
      ik = 0
      do i =1, c_dim
        do k =1, i-1
          ik = ik + 1
          bwork = bwork + metric(ik) * &
                  (val_diisdp(i, loop_m) * val_diisdp(k, j) &
                 + val_diisdp(k, loop_m) * val_diisdp(i, j))
        enddo
        ik = ik + 1
        bwork = bwork + metric(ik) * (val_diisdp(i, loop_m) *val_diisdp(i,j))
      enddo

      b(j, loop_m) = bwork
      b(loop_m, j) = bwork
    enddo

    ! b(i,i) = 0 means, error(i) = 0, here the best is no update and no storage
    ! the next loop will overwrite the values of b and the storage which were performed
    ! before
    if(b(loop_m, loop_m) .eq. 0) then
      if(step) call write_to_output_units( 'in this loop no diis charge coeff mixing, because error = 0')
      step = .false.
      loop = loop - 1
    endif

  end subroutine cc_diis_update

  !*************************************************************

 subroutine cc_diis_control(dif, start, step, d_loop, control, thresh, estep)
    ! purpose: decides what of the diis- cycle will be done:
    !          start true means that the rememberd errormatrices and the
    !          residium matrix will be updated with the new error matrix
    !          step true means that the diis-step will be performed and the
    !          parameter (charge coeff. vector) will be updated
    !          if diismix should decide on its own what to do, cc_control = 0 leads to
    !          another control routine, which does exactly this
    !          otherwise cc_control sets the decicions, cc_control is set by the diis_fock
    !          routine, which now also decides what diismix should do
    implicit none
    real(kind=r8_kind),         intent(in) :: dif(:), thresh
    integer(kind=i4_kind),    intent(in) :: control, estep
    integer(kind=i4_kind),    intent(inout) :: d_loop
    logical,               intent(inout) :: start, step
    !------------ Executable code --------------------------------
    select case (control)

      ! only charge coeff. diis mixing, so it has to deside on its own when to work
      case(0)
        call cc_own_diis_control(dif, start, step, d_loop, thresh, estep)

      ! the next three cases refer to the situation where there is first diis_fock mixing
      ! performed and after it "converged" diismixing should try to converge the charge fit coeffs. also
      case(1)
        ! diismixing not yet needed
        start = .false.
        step = .false.
        d_loop = 0
      case(2)
        ! start filling the storage
        start = .true.
        d_loop = d_loop + 1
        cc_danf = d_loop + 1
        step = .false.
      case(3)
        ! let diismixing work
        start = .true.
        d_loop = d_loop + 1
        step = .false.
        if( mod((d_loop - cc_danf),(estep + 1)) .eq. 0) step = .true.
        if(step) then
          diis_step = .true.
        else
          diis_step = .false.
        endif

      ! the next three cases refer to the situation where diis_fock mixing and
      ! diismixing are performed alternately
      ! decision whats performed when is done in the diis_control
      case(4)
        ! no diis mixing
        start = .false.
        step = .false.
        d_loop = 0
      case(5)
        ! diis_fock mixings term
        start = .true.
        d_loop = d_loop + 1
        step = .false.
      case(6)
        ! diismixings term
        start = .true.
        d_loop = d_loop + 1
        step = .true.
        if(.not. diis_step) diis_step = .true.
      case default
        call error_handler ("diismixing: no valid working case")
    end select

 end subroutine cc_diis_control

 !*************************************************************

 subroutine cc_own_diis_control(dif, start, step, d_loop, thresh, estep)
    ! purpose: decides what of the diis- cycle will be done:
    !          diis_start true means that the rememberd errormatrices and the
    !          residium matrix will be updated with the new error matrix
    !          diis_step true means that the diis-step will be performed and the
    !          parameter (charge coeff. vector) will be updated
    !          this routine decides the values on its own, the same way the diis_fock
    !          mixing routine does its control
    use options_module, only: options_switch_pertth
    implicit none
    real(kind=r8_kind),         intent(in) :: dif(:), thresh
    integer(kind=i4_kind),    intent(in) :: estep
    integer(kind=i4_kind),    intent(inout) :: d_loop
    logical,               intent(inout) :: start, step
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)                :: errmax
  !------------ Executable code --------------------------------

      ! stop_pert true means no perturbation_theory calculations
      ! when diis updates the charge coeff.
      ! this option is in this context not recommended, at least for the systems tested so far
      ! for futher details see above the case for the fock matrix
    if(stop_pert) call options_switch_pertth( pert_store)

    errmax = maxval(abs(dif))
    ! the maximum value of the errorvector should be below a threshold
    ! the loops are only counted when values for diis are stored, else it is reset to zero
    if (errmax < thresh) then
      start = .true.
      ! stop_pert true means no perturbation_theory calculations
      ! when diis updates the charge coeff
      ! for futher details see above the case for the fock matrix
      call options_switch_pertth( .false.)
      d_loop = d_loop + 1
    else
      start = .false.
      d_loop = 0
    endif
    ! steps will performed when matrices were updated and when there have been more than one
    ! step and additional only after everystep cycles of relaxation
    if(start .and. (d_loop > 1) .and. mod((d_loop - 2),&
             (estep + 1)) .eq. 0) then
      step = .true.
    else
      step = .false.
    endif
    write(output_unit,*) 'charge coeff. DIIS - loop ', d_loop, 'has maximum value', errmax

    ! tell output about control result
    if(start) then
      call  write_to_output_units( 'max(error) < threshold, therefore charge coeff. diis runs')
      if(onlyonce) then
        call write_to_trace_unit("charge coeff. DIIS started")
        onlyonce = .false.
      endif
    endif
    if(start .and. (.not. step)) then
      call write_to_output_units( 'no charge coeff. update with diis because:')
      if(d_loop .le. 2) then
        call write_to_output_units( 'not enough values for diis-mixing')
      elseif(.not. (mod((d_loop - 2), (estep + 1)) .eq. 0)) then
        call write_to_output_units( 'relaxation of system')
        call write_to_trace_unit("charge coeff. DIIS: relaxation of system")
      endif
    endif
    if(.not. start .and. .not. onlyonce) then
      onlyonce = .true.
      call write_to_trace_unit("charge coeff. DIIS stopped")
    endif

  end subroutine cc_own_diis_control

 !*************************************************************

  function diis_pertstop () result (ps_log)
  ! Purpose: gives to other modules the value of the variable stop_pert
  !          so they can find out if diis want to stop perturbation_theory
  !          calculations later
    logical           :: ps_log
  !------------ Executable code --------------------------------
    ps_log = stop_pert

  end function diis_pertstop

  !--------------- End of module ----------------------------------
end module diis_fock_module

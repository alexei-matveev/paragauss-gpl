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
module mixing_module
!-------------- Module specification ---------------------------
!
!  Purpose: mixing of the charge and exchane fitting functions and
!           of the numerical XC-Hamiltonian
!
!  Module called by: main_scf,xc_hamiltonian
!
!  Reference:   Zerner and Hehenberger CPL 62, 550 (1979)
!
!  Description:
!
!  AIN(N+1)=AOUT(N)*(1-ALPHA)+ALPHA*AIN(N); A=CHarge or eXChange
!  BETA=1-ALPHA
!
!  FIXED MIXING: BETA=CHMIX or BETA=XCMIX
!
!  DYNAMIC MIXING ACCORDING ZERNER AND HEHENBERGER CPL 62, 550 (1979)
!
!  ALPHA=MMEAN/(MMEAN-1)
!
!  MMEAN=(AOUT(N)-AOUT(N-1))/(AIN(N)-AIN(N-1))
!  AVERAGED OVER ALL COEFFICIENTS
!
!    MS 4.9.1995
!
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
! Modification (Please copy before editing)
! Author: UB
! Date:   9/7/97
! Description: Option "DISCARD_INITialized fit coefficients after
!              the first scf cycle" introduced.
!              Option "mix SPIN_SEPERATEly" introduced.
!
! Modification (Please copy before editing)
! Author: MM
! Date:   6/98
! Description: Extension to Spin Orbit
!
! Modification (Please copy before editing)
! Author: AN
! Date:   ..5/09.
! Description: exchange old DIIS routine through a new one,
!              situated in the diis_fock_module
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!------------ Modules used --------------------------------------
# include"def.h"
  use type_module
  use input_module
  use fit_coeff_module,   only : coeff_charge                                  &
                               , coeff_charge_old                              &
                               , coeff_spin                                    &
                               , coeff_spin_old                                &
                               , coeff_xc                                      &
                               , coeff_xc_old                                  &
                               , coeff_proj                                    &
                               , chold_initialized                             &
                               , spold_initialized                             &
                               , xcold_initialized                             &
                               , pv_initialized                                &
                               , fit_coeff_set_ch                              &
                               , fit_coeff_n_ch                                &
                               , fit_coeff_n_xc
  use iounitadmin_module
  use operations_module,  only : operations_fitbasis_opt
  implicit none
  private
  save
!== Interrupt end of public interface of module =================


!------------ public functions and subroutines ------------------
  public :: mixing_ch
  public :: mixing_xc
  public :: mixing_ham
  public :: mixing_proj
  public :: mixing_read
  public :: mixing_write_input
  public :: mixing_close
  public :: mixing_discard_init
  public :: mixing_state_store
  public :: mixing_state_recover
  public :: mixing_ch_ignore_coeff_old
#ifdef WITH_SCFCONTROL
  public :: mixing_read_scfcontrol
#endif
  public :: mixing_reset_buffers
  public :: mixing_beta_ch
  public :: mixing_beta_sp
  public :: mixing_beta_xc


!================================================================
! End of public interface of module
!================================================================

  !------------ Default values for input parameters -------------
  real(kind=r8_kind) :: df_chmix = 0.03_r8_kind, &
                        df_spmix = 0.03_r8_kind, &
                        df_xcmix = 0.03_r8_kind
  real(kind=r8_kind) :: df_lvshift_energy = 0.0_r8_kind
  integer(kind=i4_kind) :: df_start_after_cycle = 10
  logical :: df_discard_init  = .false., &
             df_spin_separate = .false., &
             df_global_mixing = .false.

  !------------ mixing input parameters -------------------------
  real(kind=r8_kind) :: chmix, &
                        spmix, &
                        xcmix
  real(kind=r8_kind),public :: lvshift_energy,level_shift
  integer(kind=i4_kind)       ::    idystr            ! kept for compatibility
  integer(kind=i4_kind),public::    start_after_cycle
  logical :: discard_init , &
             spin_separate, &
             global_mixing

  real (r8_kind), allocatable :: CHOLDIN(:), CHOLDOUT(:), &
       XCOLDIN(:,:), XCOLDOUT(:,:), &
       SPOLDIN(:), SPOLDOUT(:)
  real(kind=r8_kind)    :: BETACHOLD,BETAXCOLD,BETASPOLD
  real(kind=r8_kind)    :: zero = 0.0_r8_kind, one = 1.0_r8_kind

  namelist /mixing/ chmix,xcmix,spmix,idystr,discard_init,spin_separate, &
                    global_mixing,start_after_cycle, &
                    lvshift_energy

  character(len=7) :: method ! 'Dynamic', 'Global', or 'Fixed'
                             ! 'Diismix'
  ! variables to hold the new state of the spin_separate and idystr variable
  ! until a call of mixing_reset_buffers
  integer(kind=i4_kind) :: idystr_new
  logical               :: spin_separate_new

contains

  subroutine mixing_read (init_beta_values)
    ! Purpose: read the input for the mixing. It consits of three variables:
    !          chmix: mixing factor for total charge density coefficients
    !          spmix: mixing factor for spin charge density coefficients
    !                 (only used in case of a spin-polarized MDA calculation)
    !          xcmix: mixing factor for xc-coefficients; obsolete in the case
    !                 of a numerical calculation of xc-Hamiltonian
    !          idystr: number of the scf-cycle after which dynamical mixing
    !                  should start(>=3)
    use input_module
    use options_module, only: options_xcmode, &
                              xcmode_model_density, xcmode_extended_mda, &
                              options_n_spin, options_recover, &
                              recover_nothing, &
                              recover_fragment,lvshift_mixing
    logical, optional :: init_beta_values ! [default = .true.]
    !** End of interface ***************************************

    integer(kind=i4_kind) :: unit,status
    logical               :: init_beta

    if (present(init_beta_values)) then
       init_beta = init_beta_values
    else
       init_beta = .true.
    endif

    ! ------ Re-Define input dependent default values ----------
    if ( options_xcmode() == xcmode_model_density .or. &
         options_xcmode() == xcmode_extended_mda ) then
       df_discard_init = .true.
       df_spin_separate = .true.
       df_global_mixing = .true.
    elseif (options_recover() /= recover_nothing) then
       df_discard_init = .true.
    endif

    chmix         = df_chmix
    spmix         = df_spmix
    xcmix         = df_xcmix
    discard_init  = df_discard_init
    spin_separate = df_spin_separate
    global_mixing = df_global_mixing
    idystr            = -1 ! kept for compatibility
    start_after_cycle = -1 ! temporary default (see below)
    ! -----------------------------------

     if ( input_line_is_namelist("mixing") ) then
        call input_read_to_intermediate
        unit= input_intermediate_unit()
        read(unit, nml=mixing, iostat=status)
        if (status .gt. 0) call input_error( &
             "mixing_read: namelist mixing")
     endif

    lvshift_mixing=lvshift_energy.gt.0.01_r8_kind
    if(lvshift_mixing) level_shift=lvshift_energy/27.211658_r8_kind

    spin_separate = spin_separate .and. options_n_spin() > 1 .and. &
                    ( options_xcmode() == xcmode_model_density .or. &
                      options_xcmode() == xcmode_extended_mda )
    if (.not.spin_separate) spmix = chmix

    if (start_after_cycle == -1) then
       if (idystr == -1) then
          start_after_cycle = df_start_after_cycle
          idystr            = df_start_after_cycle
       else ! idystr has been specified in the input
          start_after_cycle = idystr
       endif
    else ! start_after_cycle has been specified in the input
       if (idystr == -1) then
          idystr = start_after_cycle
       else ! idystr has also been specified in the input
          if (idystr /= start_after_cycle) call input_error("mixing_read: &
               &Incompatible setting of IDYSTR and START_AFTER_CYCLE")
       endif
    endif

   if (options_recover() == recover_fragment) discard_init = .true.


    if (idystr.le.1) then
       call input_error("mixing_read: please choose START_AFTER_CYCLE > 1 ")
    endif
    if (init_beta) then
       betachold=chmix
       betaspold=spmix
       betaxcold=xcmix
    endif
    idystr_new = idystr
  end subroutine mixing_read


  subroutine mixing_write_input(iounit,scfcontrol)
    ! purpose: write the mixing namelist with its
    !          default values to the input.out
    !          file.
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    logical,optional,intent(in) :: scfcontrol
    !** End of interface ***************************************

    if (present(scfcontrol)) then
       if (scfcontrol) goto 1000
    endif

    call start("MIXING","MIXING_WRITE_INPUT", &
         iounit,operations_echo_input_level)
    call real("CHMIX            ",chmix        ,df_chmix            )
    call real("SPMIX            ",spmix        ,df_spmix            )
    call real("XCMIX            ",xcmix        ,df_xcmix            )
    call real("LVSHIFT_ENERGY   ",lvshift_energy,df_lvshift_energy  )
    call intg("START_AFTER_CYCLE",idystr       ,df_start_after_cycle)
    call flag("DISCARD_INIT     ",discard_init ,df_discard_init     )
    call flag("SPIN_SEPARATE    ",spin_separate,df_spin_separate    )
    call flag("GLOBAL_MIXING    ",global_mixing,df_global_mixing    )
    call stop(empty_line=.false.) !!!!!!!!!!!!!AS
    write(iounit,'()')

    return

    1000 continue ! entry point for the "scfcontrol" mode

    call start("MIXING","MIXING_WRITE_INPUT",iounit,echo_level_full)
    call real("CHMIX            ",chmix        ,df_chmix            )
    call real("SPMIX            ",spmix        ,df_spmix            )
    call real("XCMIX            ",xcmix        ,df_xcmix            )
    call intg("START_AFTER_CYCLE",idystr       ,df_start_after_cycle)
    call flag("SPIN_SEPARATE    ",spin_separate,df_spin_separate    )
    call flag("GLOBAL_MIXING    ",global_mixing,df_global_mixing    )
    call stop()

  end subroutine mixing_write_input


  subroutine mixing_close (loop, exch_fit)
    ! Purpose : Deallocate some module private variables
    ! -------- declaration of formail parameters --------------
    integer (i4_kind), intent (in) :: loop
    logical, intent (in), optional :: exch_fit
    !** End of interface ***************************************

    integer (i4_kind) :: iter_ch, iter_xc, alloc_stat
    logical :: xcfit
    external error_handler

    iter_ch = loop
    iter_xc = loop - 1 ! see call of mixing_xc in build_xcfit
    if (present (exch_fit)) then
      xcfit = exch_fit
    else
      xcfit = .false.
    endif

    ! Note, that mixing_ch has not yet been called in the last SCF cycle
    if (iter_ch > idystr - 1) then
       if (allocated (choldin)) then
          deallocate (choldin, stat=alloc_stat)
          ASSERT (alloc_stat==0)
       endif
       if (allocated (choldout)) then
          deallocate (choldout, stat=alloc_stat)
          ASSERT (alloc_stat==0)
       endif

       if (spin_separate) then
          if (allocated (spoldin)) then
             deallocate (spoldin, stat=alloc_stat)
             ASSERT (alloc_stat==0)
          endif
          if (allocated (spoldout)) then
             deallocate (spoldout, stat=alloc_stat)
             ASSERT (alloc_stat==0)
          endif
       endif
    endif

    ! Note, that mixing_xc has already been called in the last SCF cycle
    if (xcfit .and. (iter_xc >= idystr - 1) )then
       if (allocated (xcoldin)) then
          deallocate (xcoldin, stat=alloc_stat)
          ASSERT (alloc_stat==0)
       endif
       if (allocated (xcoldout)) then
          deallocate (xcoldout, stat=alloc_stat)
          ASSERT (alloc_stat==0)
       endif
    endif
  end subroutine mixing_close

  logical function mixing_ch_ignore_coeff_old(iter)
    ! Purpose: .TRUE. if the old charge density fit coefficients are
    !                 discarded during mixing in the current SCF cycle
    implicit none
    integer(kind=i4_kind),intent(in) :: iter ! number of the current scf-cycle
    !** End of interface ***************************************
    mixing_ch_ignore_coeff_old = iter <= idystr-1 .and. &
         chold_initialized .and. discard_init
  end function mixing_ch_ignore_coeff_old

  subroutine mixing_ch(iter, b, metric)
    ! Purpose: Mixing of charge coefficients
    use options_module, only: options_n_spin, options_xcmode, &
                              xcmode_model_density, xcmode_extended_mda
    use diis_fock_module, only: diis_charge_coeff, diis_step, diis_charge_mixing &
                              , diis_fock_beta
    implicit none
    integer(i4_kind), intent(in)  :: iter ! number of the current scf-cycle
    real(r8_kind),    intent(out) :: b ! returns the actual value of beta
    real(r8_kind),    intent(in)  :: metric(:) ! the metric tensor
    !** End of interface ***************************************

    integer   :: ALLOC_ERROR
    logical   :: spin_coeff_also
    external error_handler

    !
    ! Do DIIS  charge mixing  (only if requested)  and set  the mixing
    ! ratio b:
    !
    if (diis_charge_mixing) then
      !
      ! Function  diis_charge_coeff() return  coeff_charge  array that
      ! was generated  by mixing with  a diis routine,  implemented in
      ! diis_fock module.  It uses  the difference between old and new
      ! coeff_charge as  error vector. There is  no additional mixing,
      ! because diis_charge_coeff() uses  for the following mixing the
      ! case with true diis_step, like in the case of the Fock matrix,
      ! see below
      !
      ! FIXME:  this sub  either does  mixing, or  just  saves results
      ! internally to  perform mixing later.  There is no way  to tell
      ! which choice  was made, at  the moment except  to (unreliably)
      ! guess from the return value of mixing ratio b:
      !
      call fit_coeff_set_ch (diis_charge_coeff (coeff_charge, &
           coeff_charge_old, metric, b))

      !
      ! FIXME:   Is    it   so?    We   are   getting    "beta"   from
      ! diis_charge_coeff() ...   How do "events" diis_step  = T/F and
      ! beta =  1.0/not 1.0 correlate?  If  they do, why  isnt it made
      ! explicit?
      !
      if (.not. diis_step) then
        ASSERT(b==1.0)
      endif

      if (spin_separate) then
        WARN('WARNING: DIIS not adapted for the spin_separate')
      endif
    else
      !
      ! FIXME: communication via global module state!
      !
      b = diis_fock_beta
    endif

    ! FIXME:  empirical  magic:  function  DYNDAMP  shows  that  b  is
    ! restricted to range: 0.1 0.9  do not produce smaller values, for
    ! normal charge mixing it is  not possible to have negative mixing
    ! values, thus the maximal value here is 1
    if (b < 0.1_r8_kind) b = 0.1_r8_kind
    if (b > 1.0_r8_kind) b = 1.0_r8_kind

    !
    ! In case DIIS  took over the job of  mixing (either charge coeffs
    ! or the Fock matrix) do only the necessary minimum:
    !
    if ( diis_step ) then
      !
      ! * diis_step   is    flipped   to   true    either   when   (A)
      !   diis_charge_coeff() did  non-trivial mixing (correlates with
      !   beta returned by it being  not identically equal to 1.0), or
      !   (B)  when diis_fock_matrix() that  is called  from elsewhere
      !   did  non-trivial mixing  of the  hamiltonian. In  the latter
      !   case  we have  no  access to  the corresponding  self-mixing
      !   coefficient. We need to request  it from the global state of
      !   the diis_fock_module (see above)
      !
      ! * When DIIS  mixes the Fock matrix, there  shouldn't be charge
      !   mixing.  In this case b is set to diis_fock_beta, a variable
      !   provided by the DIIS module.  This helps in the case that in
      !   the next loop there should be charge mixing, when beta would
      !   be  used as  the  old mixing  coefficient  in the  dynamical
      !   mixing.  Don't bother about it with DIIS_ON = false (no DIIS
      !   running) then diis_step is always false.
      !
      BETACHOLD = b

      !
      ! ... do nothing else:
      !
      return
    endif

    !
    ! FRACTION MIXING (several ways to define the fraction)
    !

    spin_coeff_also = options_n_spin() > 1 .and. &
         (options_xcmode() == xcmode_model_density .or. &
         options_xcmode() == xcmode_extended_mda)

    !
    ! FIXME: find a better blace for initialization:
    !
    if (.not. allocated (CHOLDIN)) then
      allocate(CHOLDIN(fit_coeff_n_ch()), CHOLDOUT(fit_coeff_n_ch()), STAT=ALLOC_ERROR)
      if(ALLOC_ERROR/=0) call error_handler&
           ('Allocation failed in subroutine mixing')

      if (spin_separate) then
         allocate(SPOLDIN(fit_coeff_n_ch()), SPOLDOUT(fit_coeff_n_ch()), STAT=ALLOC_ERROR)
         if(ALLOC_ERROR/=0) call error_handler&
              ('Allocation failed in subroutine mixing')
      endif
    endif

    !
    ! Determine the fraction:
    !
    if (ITER < IDYSTR) then
      if (ITER <= 1) then
        !
        ! Initial coefficients may not even be normalized, discard:
        !
        BETACHOLD = ONE
      else
        BETACHOLD = CHMIX
      endif
    else ! if ITER >= IDYSTR
      !
      ! Update BETACHOLD in place:
      !
      BETACHOLD = legacy_beta(coeff_charge_old, coeff_charge, CHOLDIN, CHOLDOUT, metric, CHMIX, BETACHOLD)
    endif

    !
    ! Mix old and new into new:
    !
    call fit_coeff_set_ch(mix(BETACHOLD, coeff_charge_old, coeff_charge))

    if (spin_coeff_also) then
      if (spin_separate) then
        !
        ! Analogously to BETACHOLD
        !
        if (ITER < IDYSTR) then
          if (ITER <= 1) then
            BETASPOLD = ONE
          else
            BETASPOLD = SPMIX
          endif
        else ! if ITER >= IDYSTR
          !
          ! Update BETASPOLD in place:
          !
          BETASPOLD = legacy_beta(coeff_spin_old, coeff_spin, SPOLDIN, SPOLDOUT, metric, SPMIX, BETASPOLD)
        endif
      else
        !
        ! use same as for charge:
        !
        BETASPOLD = BETACHOLD
      endif

      coeff_spin = mix(BETASPOLD, coeff_spin_old, coeff_spin)
    endif

    !
    !  PREPARATION FOR DYNAMIC MIXING IN NEXT ITERATION
    !
    ! PREPARE FOR DYNAMICAL MIXING IN NEXT CYCLE
    CHOLDIN = coeff_charge_old
    CHOLDOUT = coeff_charge

    if (spin_separate) then
       SPOLDIN = coeff_spin_old
       SPOLDOUT = coeff_spin
    endif

    ! used for orbital mixing in PT:
    b = BETACHOLD
  end subroutine mixing_ch

  subroutine legacy(ain, aout, aoldin, aoldout, metric, beta0, beta)
    !
    ! DYNAMIC MIXING
    !
    implicit none
    real(r8_kind), intent(in) :: ain(:)
    real(r8_kind), intent(inout) :: aout(:)
    real(r8_kind), intent(inout) :: aoldin(:), aoldout(:)
    real(r8_kind), intent(in) :: metric(:)
    real(r8_kind), intent(in) :: beta0
    real(r8_kind), intent(inout) :: beta
    ! *** end of interface ***

    ! update beta in place:
    beta = legacy_beta(ain, aout, aoldin, aoldout, metric, beta0, beta)

    aoldin = ain
    aoldout = aout

    aout = mix(beta, ain, aout)
  end subroutine legacy

  function legacy_beta(ain, aout, aoldin, aoldout, metric, beta0, betaold) result(beta)
    !
    ! DYNAMIC MIXING
    !
    implicit none
    real(r8_kind), intent(in) :: ain(:), aout(:)
    real(r8_kind), intent(in) :: aoldin(:), aoldout(:)
    real(r8_kind), intent(in) :: metric(:)
    real(r8_kind), intent(in) :: beta0, betaold
    real(r8_kind)             :: beta ! result
    ! *** end of interface ***

    if ( global_mixing ) then
      beta = GLOBAL(ain, aout, aoldin, aoldout, metric, beta0, "coefficients:")
    else
      beta = DYNDAMP(ain, aout, aoldin, aoldout)

      if (beta == -one) beta = beta0

      beta = (beta + betaold) / 2
    endif
  end function legacy_beta

  function mix(beta, aold, acur) result(aout)
    implicit none
    real(r8_kind), intent(in) :: beta
    real(r8_kind), intent(in) :: aold(:)
    real(r8_kind), intent(in) :: acur(:)
    real(r8_kind) :: aout(size(acur))
    ! *** end of interface ***

    aout = acur * beta + (one - beta) * aold
  end function mix

  subroutine mixing_xc(ITER,metric)
    ! Purpose: mixing of the xc-coefficients
    use output_module, only: output_scfloops
    use diis_fock_module, only: diis_on
    implicit none
    integer(kind=i4_kind),intent(IN)   :: ITER
    real(kind=r8_kind),intent(in),optional :: metric(:) ! the metric tensor
    !** End of interface **************************************
    real(kind=r8_kind)              ::  BETA,BETAXC,b,BETAUP,BETADN
    integer         ::  ALLOC_ERROR
    external error_handler

    ! FIXME: needs some work:
    ASSERT(.not.diis_on)

    if (ITER<=IDYSTR-1) then
       !mixing_xc
       !  FIXED MIXING
       !
       if (ITER==IDYSTR-1)  then
          !
          !  PREPARATION FOR DYNAMIC MIXING IN NEXT ITERATION
          !
          allocate( XCOLDIN(fit_coeff_n_xc(),size(coeff_xc,2))                 &
                  , XCOLDOUT(fit_coeff_n_xc(),size(coeff_xc,2))                &
                  , STAT=ALLOC_ERROR                                           )
          if(ALLOC_ERROR/=0) call error_handler&
               ('Allocation failed in subroutine mixing')
          XCOLDIN=coeff_xc_old
          XCOLDOUT=coeff_xc
       endif
       if (xcold_initialized .and. discard_init) then
          b=one
          if (iter == idystr-1) betaxcold=one
       else
          coeff_xc=coeff_xc*XCMIX+(one-XCMIX)*coeff_xc_old
          b=xcmix
       endif
       BETAXCOLD=b
       method = 'Fixed'
    else
       !
       ! DYNAMIC MIXING
       !
       if (global_mixing .and. present(metric)) then
          if (size(coeff_xc,1) == 1) then
             BETA=GLOBAL(coeff_xc_old(:,1),coeff_xc(:,1),XCOLDIN(:,1), &
                         XCOLDOUT(:,1),metric,xcmix,'Exchange coefficients:')
          else
             BETAUP=GLOBAL(coeff_xc_old(:,1),coeff_xc(:,1),XCOLDIN(:,1), &
                           XCOLDOUT(:,1),metric,xcmix,'XC(up)   coefficients:')
             BETADN=GLOBAL(coeff_xc_old(:,2),coeff_xc(:,2),XCOLDIN(:,2), &
                           XCOLDOUT(:,2),metric,xcmix,'XC(down) coefficients:')
             BETA=(BETAUP+BETADN)*0.5_r8_kind
          endif
          method = 'Global'
       else
       if (output_scfloops) &
            write(output_unit,'(A)')'Dynamic damping for exchange coefficients:'
       BETAXC=DYNDAMP(reshape(coeff_xc_old,(/size(coeff_xc)/)),&
            reshape(coeff_xc,(/size(coeff_xc)/)),&
            reshape(XCOLDIN,(/size(coeff_xc)/))&
            ,reshape(XCOLDOUT,(/size(coeff_xc)/)))
       if (BETAXC==-one) then
          BETAXC=XCMIX
          BETA=(BETAXC+BETAXCOLD)/2
          if (output_scfloops) &
               write(output_unit,'(1X,A25,F8.3)') 'RESULTING BETA:',BETA
       else
          BETA=(BETAXC+BETAXCOLD)/2
       endif
          method = 'Dynamic'
       endif
       XCOLDIN=coeff_xc_old
       XCOLDOUT=coeff_xc
       coeff_xc=coeff_xc*BETA+(one-BETA)*coeff_xc_old
       b=BETA
       BETAXCOLD=BETA
    end if
    if (output_scfloops) &
         write(output_unit,'(1X,A,1X,A,F8.3)') &
         trim(method),'mixing for exchange coefficients, Beta:',b

  end subroutine mixing_xc

  !***********************************************************

  subroutine mixing_ham(ham_xc_arr, ham_xc_arr_old, matold_initialized)
    !    PURPOSE:
    !    Mixing of the numerical xc-Hamiltonion. Beta is the same as
    !    for the chargefit case.
    !------------ Modules ----------------------------------------
    use output_module, only: output_scfloops
    use diis_fock_module, only: diis_step
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(inout) :: ham_xc_arr(:)
    real(kind=r8_kind),intent(in)    :: ham_xc_arr_old(:)
    logical, optional, intent(in)    :: matold_initialized
    !** End of interface ***************************************

    integer, save :: loop = 0
    logical :: discard

    loop = loop + 1

    if (present(matold_initialized)) then
       discard = matold_initialized .and. discard_init
    else
       discard = .false.
    endif

    ! no additional mixing if DIIS mixes the Fock_matrix:
    discard = discard .or. diis_step

    if ( loop <= 1 ) then
        WARN('FIXME: skip xc mixing in first SCF')
        discard = .true.
    endif

    if (discard) then
       ! nothing
       if (output_scfloops) write(output_unit,'(A40)') &
            'Mixing for numerical Hamiltonian, Skip'
    else
       ham_xc_arr = ham_xc_arr * betachold + (one - betachold) * ham_xc_arr_old

       if (output_scfloops) write(output_unit,'(A40, F8.3)') &
            'Mixing for numerical Hamiltonian, Beta:', betachold
    endif
  end subroutine mixing_ham

 !**************************************************************
!:UB[ added

  subroutine mixing_proj(coeff_proj_new)
    !    PURPOSE:
    !    Mixing of the charge projection for the ext. MDA.
    !    Beta is the same as for the chargefit case.
    !------------ Modules ----------------------------------------
    use output_module, only: output_scfloops
    use diis_fock_module, only: diis_step
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in) :: coeff_proj_new(:,:)
    !** End of interface ***************************************
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind) :: b
    !------------ Executable code --------------------------------
    if (pv_initialized .and. discard_init) then
!:UB[ added
       coeff_proj = coeff_proj_new
!:UB]
       b = one
    elseif (diis_step) then
       b = one
    else
       coeff_proj = coeff_proj_new*betachold + (one-betachold)*coeff_proj
       b = betachold
    endif
    if (output_scfloops) &
         write(output_unit,'(A40,F8.3)') &
         'Mixing for charge projection, Beta:', b

  end subroutine mixing_proj

 !**************************************************************
!:UB]

 function DYNDAMP(AIN,AOUT,AOLDIN,AOLDOUT)
   ! Purpose: Calculation of the optimal mixing parameter
   use output_module, only: output_scfloops
   implicit none
   real(kind=r8_kind) :: DYNDAMP
   real(kind=r8_kind),intent(IN),dimension(:)  :: AIN,AOUT,AOLDIN,AOLDOUT
   !** End of interface **************************************

   real(r8_kind), dimension(size(AIN)) :: NOM, DNOM, M
   real(kind=r8_kind):: MMEAN,MSQUARE,SIGMA
   integer(i4_kind) :: H(size(AIN)), NR

   H = 1
   DNOM=AIN-AOLDIN

   NOM=AOUT-AOLDOUT
   where(abs(DNOM)<1E-10_r8_kind)
      DNOM=one
      NOM=zero
      H = 0
   ENDWHERE

   M=NOM/DNOM

   where((M>zero).OR.(M<-19.0_r8_kind))
      M=zero
      H = 0
   ENDWHERE

   NR = sum(H)

   if(NR>1) then
      MMEAN=sum(M)/NR
      MSQUARE=sum(M*M)/NR
!!$   SIGMA=sqrt(MSQUARE-MMEAN*MMEAN)/NR
      SIGMA=sqrt(MSQUARE-MMEAN*MMEAN)
      DYNDAMP=one-MMEAN/(MMEAN-one)
      if(DYNDAMP>0.9_r8_kind) DYNDAMP=0.9_r8_kind
      if(DYNDAMP<0.1_r8_kind) DYNDAMP=0.1_r8_kind
      if (output_scfloops) &
         write(output_unit,1000)NR,DYNDAMP,MMEAN,SIGMA
   else
      DYNDAMP=-one
      if (output_scfloops) &
         write(output_unit,1000)NR,DYNDAMP
   endif

   1000 Format('AK_USED = ',I5,'  DYNBETA = ',F8.3: &
               '  AVERAGE SLOPE = ',F8.3,'  SIGMA = ',F8.3)
 end function DYNDAMP

  !***********************************************************

 function GLOBAL(AIN,AOUT,AOLDIN,AOLDOUT,METRIC,DEFAULT,TEXT)
   ! Purpose: Calculation of an optimal global mixing parameter
   !          using the following algorithm: (Uwe Birkenheuer, 14/10/97)
   !
   !          da_in  := aold_in  - a_in
   !          da_out := aold_out - a_out
   !          da     := a_out - a_in
   !
   !          a_new(a) := a_out + da_out * < da_in | metric | a - a_in > / N
   !          with
   !          N = < da_in | metric | da_in >
   !
   !          Remark: a_new(a_in   ) = a_out
   !                  a_new(aold_in) = aold_out
   !
   !          a_next := beta*a_out + (1-beta)*a_in = a_in + beta*da
   !
   !          a_calc := a_new(a_next)
   !                  = a_out + beta * da_out * < da_in | metric | da > / N
   !
   !          da_new := a_calc - a_next
   !                  = da + beta*{ da_out * < da_in | metric | da > / N - da }
   !                  = da + beta * da_eff
   !          with
   !          da_eff = da_out * < da_in | metric | da > / N - da
   !
   !          target := < da_new | metric | da_new >
   !
   !          target =! min  ==>  d/d beta target(beta) =! 0 :
   !
   !          beta = - < da | metric | da_eff > / < da_eff | metric | da_eff >
   !
   use output_module, only: output_scfloops
   implicit none
   real(kind=r8_kind)                          :: GLOBAL
   real(kind=r8_kind),intent(IN),dimension(:)  :: AIN,AOUT,AOLDIN,AOLDOUT
   real(kind=r8_kind),intent(IN),dimension(:)  :: METRIC ! packed storage mode
   real(kind=r8_kind),intent(IN)               :: DEFAULT
   character(len=*)  ,intent(IN)               :: TEXT
   !** End of interface **************************************
   real(kind=r8_kind), dimension(size(ain)) :: da, da_eff ! automatic arrays
   real(kind=r8_kind)                       :: norm, proj

   real(kind=r8_kind), parameter            :: tiny = 1.0E-16_r8_kind

   ! first compute da_in and N = < da_in | metric | da_in >
   ! using da_eff as working array
   da_eff = aoldin - ain
   norm = dot (da_eff, da_eff)
   if (output_scfloops) write(output_unit,1000)trim(text),sqrt(norm)
   1000 format(1X,A,' ||a_in - a_in(old)|| = ',ES9.2)
   if (norm < tiny*tiny) then ! a_new(a) := a_out  =>  beta = 1.0
      global = 1.0_r8_kind
      return
   endif

   ! now compute da and da_eff
   da = aout - ain
   proj = dot (da_eff, da) / norm
   da_eff = ( aoldout - aout ) * proj - da

   ! finally compute beta
   norm = dot (da_eff, da_eff)
   if (norm < tiny*tiny ) then
      if (output_scfloops) write(output_unit,'(1X,A,1X,A)')trim(text), &
           'No a_k progression along a_out-a_in; default value taken'
      global = default
   else
      global = - dot (da_eff, da) / norm
   endif

   contains

   function dot (x, y)
     !
     ! Scalar product with a metric <x|Gy>.
     !
     real (r8_kind), intent (in) :: x(:), y(:)
     real (r8_kind) :: dot
     ! *** end of interface ***

     integer (i4_kind) :: i, j, lin

     dot = 0.0
     lin = 0
     do i = 1, size (x)
        do j = 1, i - 1
           lin = lin + 1
           dot = dot + (x(i) * y(j) + x(j) * y(i)) * metric(lin)
        enddo
        lin = lin + 1
        dot = dot + x(i) * y(i) * metric(lin)
     enddo
   end function dot

 end function GLOBAL

  !***********************************************************

  logical function mixing_discard_init()
    !    PURPOSE: gives access to the control parameter discard_init
    !    TRUE  -- whenever the old data to be mixed has only been
    !             initialized (arbitrarily) use a mixing factor of 1.
    !    FALSE -- always use the current mixing factor
    !** End of interface *****************************************
    mixing_discard_init = discard_init
  end function mixing_discard_init

 !**************************************************************

  subroutine mixing_state_store(loop,th,mode)
    ! Purpose: stores the current state of the internal buffer of
    !          the mixing_module to a readwriteblocked file in case
    !          any of the "save_..." options is used.
    !
    !
    ! << th >>     << mode >>   action
    ! PRESENT      NOT PRESENT  store present state immediatelly
    ! NOT PRESENT  PRESENT      keep present state for later storage
    ! PRESENT      PRESENT      store previously saved state
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB, 8/97
    !------------ Modules ----------------------------------------
    use options_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_saved
    use readwriteblocked_module
    !------------ Declaration of formal paramaters ---------------
    integer(kind=i4_kind), intent(in) :: loop
    type(readwriteblocked_tapehandle), optional, intent(inout) :: th
    integer(kind=i4_kind), optional, intent(in) :: mode
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: iter_ch, iter_xc, n_spin
    real(kind=r8_kind) :: yes(1) = (/1.0_r8_kind/), no(1) = (/0.0_r8_kind/)
    logical :: spin_coeff_also
    real(r8_kind), save, pointer :: chold_kept(:,:), spold_kept(:,:), &
                                    xcold_kept(:,:,:)
    real(r8_kind), save :: betach_kept, betasp_kept, betaxc_kept
    logical, save :: charge_and_spin_kept = .false., &
                     xc_and_en_xc_kept    = .false.
    logical  :: reset
    integer  :: status
    external error_handler
    !------------ Executable code --------------------------------
    n_spin = options_n_spin()
    iter_ch = loop
    iter_xc = loop - 1 ! see call of mixing_xc in build_xcfit
    spin_coeff_also = n_spin > 1 .and. &
                      ( options_xcmode() == xcmode_model_density .or. &
                        options_xcmode() == xcmode_extended_mda )

    if (.not.present(th)) then
       ! keep modus
       if (mode == recover_eigenvec) then
          if (iter_ch > idystr - 1) then
             ! save beta(ch/sp), chold(in/out) and spold(in/out)
             betach_kept = betachold
             if (.not.charge_and_spin_kept) then
                allocate(chold_kept(fit_coeff_n_ch(),2),stat=status)
                if (status /= 0) call error_handler &
                     ("MIXING_STATE_STORE: allocation of chold_kept failed")
             endif
             chold_kept(:,1) = choldin
             chold_kept(:,2) = choldout
             if (spin_separate) then
                betasp_kept = betaspold
                if (.not.charge_and_spin_kept) then
                   allocate(spold_kept(fit_coeff_n_ch(),2),stat=status)
                   if (status /= 0) call error_handler &
                        ("MIXING_STATE_STORE: allocation of spold_kept failed")
                endif
                spold_kept(:,1) = spoldin
                spold_kept(:,2) = spoldout
             endif
             charge_and_spin_kept = .true.
          endif
          if (options_xcmode() == xcmode_exchange_fit .and. &
               iter_xc >= idystr - 1) then
             ! save betaxc and xcold(in/out)
             betaxc_kept = betaxcold
             if (.not.xc_and_en_xc_kept) then
                allocate(xcold_kept(fit_coeff_n_xc(),n_spin,2),stat=status)
                if (status /= 0) call error_handler &
                     ("MIXING_STATE_STORE: allocation of xcold_kept failed")
             endif
             xcold_kept(:,:,1) = xcoldin
             xcold_kept(:,:,2) = xcoldout
             xc_and_en_xc_kept = .true.
          endif
       endif
       return
    endif

    if (output_main_scf) call write_to_output_units &
         ("MIXING_STATE_STORE: saving current mixing state")

    if (present(mode)) then
       reset = mode == recover_eigenvec
    else
       reset = .false.
    endif

    ! store current SCF loop
    call readwriteblocked_write((/real(loop,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored mixing state :'
       write(output_unit,'(  a     )')'loop'
       write(output_unit,'(4es20.13)')(/real(loop,r8_kind)/)
    endif

    ! store charge fit state:
    ! note, that mixing_ch has not yet been called in the current SCF loop
    if (iter_ch <= idystr - 1) then
       call readwriteblocked_write(no,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'dynamic charge mixing ?'
          write(output_unit,'(4es20.13)')no
       endif
    else
       call readwriteblocked_write(yes,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'dynamic charge mixing ?'
          write(output_unit,'(4es20.13)')yes
       endif
       if (reset) then
          charge_and_spin_kept = .false.
          call readwriteblocked_write((/betach_kept/),th)
          call readwriteblocked_write(chold_kept(:,1),th)
          call readwriteblocked_write(chold_kept(:,2),th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'betachold'
             write(output_unit,'(4es20.13)')(/betach_kept/)
             write(output_unit,'(  a     )')'choldin'
             write(output_unit,'(4es20.13)')chold_kept(:,1)
             write(output_unit,'(  a     )')'choldout'
             write(output_unit,'(4es20.13)')chold_kept(:,2)
          endif
          deallocate(chold_kept,stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_STORE: deallocation of chold_kept failed")
       else
          call readwriteblocked_write((/betachold/),th)
          call readwriteblocked_write(choldin,th)
          call readwriteblocked_write(choldout,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'betachold'
             write(output_unit,'(4es20.13)')(/betachold/)
             write(output_unit,'(  a     )')'choldin'
             write(output_unit,'(4es20.13)')choldin
             write(output_unit,'(  a     )')'choldout'
             write(output_unit,'(4es20.13)')choldout
          endif
       endif
       if (.not.spin_separate) then
          call readwriteblocked_write(no,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'mix coeff_spin separately ?'
             write(output_unit,'(4es20.13)')no
          endif
       else
          call readwriteblocked_write(yes,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'mix coeff_spin separately ?'
             write(output_unit,'(4es20.13)')yes
          endif
          if (reset) then
             call readwriteblocked_write((/betasp_kept/),th)
             call readwriteblocked_write(spold_kept(:,1),th)
             call readwriteblocked_write(spold_kept(:,2),th)
             if (output_data_saved) then
                write(output_unit,'(  a     )')'betaspold'
                write(output_unit,'(4es20.13)')(/betasp_kept/)
                write(output_unit,'(  a     )')'spoldin'
                write(output_unit,'(4es20.13)')spold_kept(:,1)
                write(output_unit,'(  a     )')'spoldout'
                write(output_unit,'(4es20.13)')spold_kept(:,2)
             endif
             deallocate(spold_kept,stat=status)
             if (status /= 0) call error_handler &
                  ("MIXING_STATE_STORE: deallocation of spold_kept failed")
          else
             call readwriteblocked_write((/betaspold/),th)
             call readwriteblocked_write(spoldin,th)
             call readwriteblocked_write(spoldout,th)
             if (output_data_saved) then
                write(output_unit,'(  a     )')'betaspold'
                write(output_unit,'(4es20.13)')(/betaspold/)
                write(output_unit,'(  a     )')'spoldin'
                write(output_unit,'(4es20.13)')spoldin
                write(output_unit,'(  a     )')'spoldout'
                write(output_unit,'(4es20.13)')spoldout
             endif
          endif
       endif
    endif

    ! store exchange fit state:
    ! note, that mixing_xc has already been called in the current SCF loop
    if (options_xcmode() /= xcmode_exchange_fit .or. iter_xc < idystr - 1) then
       call readwriteblocked_write(no,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'dynamic exchange mixing ?'
          write(output_unit,'(4es20.13)')no
       endif
    else
       call readwriteblocked_write(yes,th)
       call readwriteblocked_write((/real(fit_coeff_n_xc(),r8_kind)/),th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'dynamic exchange mixing ?'
          write(output_unit,'(4es20.13)')yes
          write(output_unit,'(  a     )')'n_xc'
          write(output_unit,'(4es20.13)')(/real(fit_coeff_n_xc(),r8_kind)/)
       endif
       if (reset) then
          call readwriteblocked_write((/betaxc_kept/),th)
          call readwriteblocked_write(xcold_kept(:,1,1),th)
          call readwriteblocked_write(xcold_kept(:,1,2),th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'betaxcold'
             write(output_unit,'(4es20.13)')(/betaxc_kept/)
             write(output_unit,'(  a     )')'xcoldin(:,1)'
             write(output_unit,'(4es20.13)')xcold_kept(:,1,1)
             write(output_unit,'(  a     )')'xcoldout(:,1)'
             write(output_unit,'(4es20.13)')xcold_kept(:,1,2)
          endif
       else
          call readwriteblocked_write((/betaxcold/),th)
          call readwriteblocked_write(xcoldin(:,1),th)
          call readwriteblocked_write(xcoldout(:,1),th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'betaxcold'
             write(output_unit,'(4es20.13)')(/betaxcold/)
             write(output_unit,'(  a     )')'xcoldin(:,1)'
             write(output_unit,'(4es20.13)')xcoldin(:,1)
             write(output_unit,'(  a     )')'xcoldout(:,1)'
             write(output_unit,'(4es20.13)')xcoldout(:,1)
          endif
       endif
       if (n_spin == 1) then
          call readwriteblocked_write(no,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'coeff_xc(:,2) required ?'
             write(output_unit,'(4es20.13)')no
          endif
       else
          call readwriteblocked_write(yes,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'coeff_xc(:,2) required ?'
             write(output_unit,'(4es20.13)')yes
          endif
          if (reset) then
             call readwriteblocked_write(xcold_kept(:,2,1),th)
             call readwriteblocked_write(xcold_kept(:,2,2),th)
             if (output_data_saved) then
                write(output_unit,'(  a     )')'xcoldin(:,2)'
                write(output_unit,'(4es20.13)')xcold_kept(:,2,1)
                write(output_unit,'(  a     )')'xcoldout(:,2)'
                write(output_unit,'(4es20.13)')xcold_kept(:,2,2)
             endif
          else
             call readwriteblocked_write(xcoldin(:,2),th)
             call readwriteblocked_write(xcoldout(:,2),th)
             if (output_data_saved) then
                write(output_unit,'(  a     )')'xcoldin(:,2)'
                write(output_unit,'(4es20.13)')xcoldin(:,2)
                write(output_unit,'(  a     )')'xcoldout(:,2)'
                write(output_unit,'(4es20.13)')xcoldout(:,2)
             endif
          endif
       endif
       if (reset) then
          xc_and_en_xc_kept = .false.
          deallocate(xcold_kept,stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_STORE: deallocation of xcold_kept failed")
       endif
    endif

  end subroutine mixing_state_store

!*************************************************************
! record  1: loop
! record  2: dynamic_mixing
! record  3: betachold      [if dynamic_mixing]
! record  4: choldin        [if dynamic_mixing]
! record  5: choldout       [if dynamic_mixing]
! record  6: spin_separate  [if dynamic_mixing]
! record  7: betaspold      [if dynamic_mixing & spin_separate]
! record  8: spoldin        [if dynamic_mixing & spin_separate]
! record  9: spoldout       [if dynamic_mixing & spin_separate]
! record 10: exchange_fit & dynamic_mixing
! record 11: n_xc           [if exchange_fit & dynamic_mixing]
! record 12: betaxcold      [if exchange_fit & dynamic_mixing]
! record 13: xcoldin(:,1)   [if exchange_fit & dynamic_mixing]
! record 14: xcoldout(:,1)  [if exchange_fit & dynamic_mixing]
! record 15: spin_polarized [if exchange_fit & dynamic_mixing]
! record 16: xcoldin(:,2)   [if exchange_fit & dynamic_mixing & spin_polarized]
! record 17: xcoldout(:,2)  [if exchange_fit & dynamic_mixing & spin_polarized]
!*************************************************************

  subroutine mixing_state_recover(loop,th)
    ! Purpose: recovers the current state of the internal buffer of
    !          the mixing_module from a readwriteblocked file in case
    !          one of the "read_..." options is used.
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB, 8/97
    !------------ Modules ----------------------------------------
    use options_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module, only: output_main_scf, output_data_read
    use readwriteblocked_module
    !------------ Declaration of formal paramaters ---------------
    integer(kind=i4_kind), intent(out) :: loop
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: iter_ch, iter_xc, n_spin, n_xc_stored
    real(kind=r8_kind)    :: dummy(1), nxc(1), zero = 0.0_r8_kind, &
                             dyn_mix(1), spin_sep(1), exch_dyn(1), spin_pol(1)
    logical               :: spin_coeff_also, exchange_fit, dynamic_xc_mixing
    integer               :: status
    external error_handler
    !------------ Executable code --------------------------------
    n_spin = options_n_spin()
    exchange_fit = options_xcmode() == xcmode_exchange_fit
    spin_coeff_also = n_spin > 1 .and. &
                      ( options_xcmode() == xcmode_model_density .or. &
                        options_xcmode() == xcmode_extended_mda )
    if (output_main_scf) call write_to_output_units &
         ("MIXING_STATE_RECOVER: reading current mixing state")

    ! first recover current SCF loop
    call readwriteblocked_read(dummy,th)
    loop = int(dummy(1),i4_kind)
    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered mixing state :'
       write(output_unit,'(  a     )')'loop'
       write(output_unit,'(4es20.13)')dummy(1)
    endif
    if (options_reset_scfcycle()) loop = 1
    iter_ch = loop
    iter_xc = loop - 1 ! see call of mixing_xc in build_xcfit

    ! recover charge fit state:
    ! note, that mixing_ch will still be called in the current SCF cycle
    call readwriteblocked_read(dyn_mix,th)
    if (output_data_read) then
       write(output_unit,'(  a     )')'dynamic charge mixing ?'
       write(output_unit,'(4es20.13)')dyn_mix(1)
    endif
    if (dyn_mix(1) /= zero) then
       if (iter_ch > idystr - 1) then
          call readwriteblocked_read(dummy,th)
          betachold = dummy(1)
          allocate(choldin(fit_coeff_n_ch()),choldout(fit_coeff_n_ch()),stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_RECOVER: allocation (1) failed")
          call readwriteblocked_read(choldin,th)
          call readwriteblocked_read(choldout,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betachold'
             write(output_unit,'(4es20.13)')dummy(1)
             write(output_unit,'(  a     )')'choldin'
             write(output_unit,'(4es20.13)')choldin
             write(output_unit,'(  a     )')'choldout'
             write(output_unit,'(4es20.13)')choldout
          endif
       else
          call readwriteblocked_skipread(2*fit_coeff_n_ch()+1,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betachold skipped'
             write(output_unit,'(  a     )')'choldin   skipped'
             write(output_unit,'(  a     )')'choldout  skipped'
          endif
       endif
       call readwriteblocked_read(spin_sep,th)
       if (output_data_read) then
          write(output_unit,'(  a     )')'mix coeff_spin separately ?'
          write(output_unit,'(4es20.13)')spin_sep(1)
       endif
    else
       if (iter_ch > idystr - 1) then
          if (output_main_scf) call write_to_output_units &
               ("MIXING_STATE_RECOVER: starting of dynamic mixing reset")
          idystr = iter_ch + 1
       endif
       spin_sep(1) = zero
    endif

    if (dyn_mix(1) /= zero .and. spin_sep(1) /= zero) then
       if (iter_ch > idystr - 1 .and. spin_separate) then
          call readwriteblocked_read(dummy,th)
          betaspold = dummy(1)
          allocate(spoldin(fit_coeff_n_ch()),spoldout(fit_coeff_n_ch()),stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_RECOVER: allocation (1') failed")
          call readwriteblocked_read(spoldin,th)
          call readwriteblocked_read(spoldout,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betaspold'
             write(output_unit,'(4es20.13)')dummy(1)
             write(output_unit,'(  a     )')'spoldin'
             write(output_unit,'(4es20.13)')spoldin
             write(output_unit,'(  a     )')'spoldout'
             write(output_unit,'(4es20.13)')spoldout
          endif
       else
          call readwriteblocked_skipread(2*fit_coeff_n_ch()+1,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betaspold skipped'
             write(output_unit,'(  a     )')'spoldin   skipped'
             write(output_unit,'(  a     )')'spoldout  skipped'
          endif
       endif
    elseif (iter_ch > idystr - 1 .and. spin_separate) then
       if (output_main_scf) call write_to_output_units &
            ("MIXING_STATE_RECOVER: starting of dynamic mixing reset")
       deallocate(choldin,choldout,stat=status)
       if (status /= 0) call error_handler &
            ("MIXING_STATE_RECOVER: deallocation (1a) failed")
       idystr = iter_ch + 1
    endif

    ! recover exchange fit state:
    ! note, that mixing_xc has already been called in current SCF cycle
    dynamic_xc_mixing = exchange_fit .and. iter_xc >= idystr - 1
    call readwriteblocked_read(exch_dyn,th)
    if (output_data_read) then
       write(output_unit,'(  a     )')'dynamic exchange mixing ?'
       write(output_unit,'(4es20.13)')exch_dyn(1)
    endif
    if (exch_dyn(1) /= zero) then
       call readwriteblocked_read(nxc,th)
       n_xc_stored = int(nxc(1),i4_kind)
       if (output_data_read) then
          write(output_unit,'(  a     )')'n_xc'
          write(output_unit,'(4es20.13)')nxc(1)
       endif
       if (dynamic_xc_mixing) then
          call readwriteblocked_read(dummy,th)
          betaxcold = dummy(1)
          allocate(xcoldin(fit_coeff_n_ch(),n_spin),xcoldout(fit_coeff_n_ch(),n_spin), &
               stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_RECOVER: allocation (2) failed")
          call readwriteblocked_read(xcoldin(:,1),th)
          call readwriteblocked_read(xcoldout(:,1),th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betaxcold'
             write(output_unit,'(4es20.13)')dummy(1)
             write(output_unit,'(  a     )')'xcoldin(:,1)'
             write(output_unit,'(4es20.13)')xcoldin(:,1)
             write(output_unit,'(  a     )')'xcoldout(:,1)'
             write(output_unit,'(4es20.13)')xcoldout(:,1)
          endif
       else
          call readwriteblocked_skipread(2*n_xc_stored+1,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'betaxcold     skipped'
             write(output_unit,'(  a     )')'xcoldin(:,1)  skipped'
             write(output_unit,'(  a     )')'xcoldout(:,1) skipped'
          endif
       endif
       call readwriteblocked_read(spin_pol,th)
       if (output_data_read) then
          write(output_unit,'(  a     )')'coeff_xc(:,2) required ?'
          write(output_unit,'(4es20.13)')spin_pol
       endif
    else
       if (dynamic_xc_mixing) then
          if (output_main_scf) call write_to_output_units &
               ("MIXING_STATE_RECOVER: starting of dynamic mixing reset")
          if (iter_ch > idystr - 1) then
             deallocate(choldin,choldout,stat=status)
             if (status /= 0) call error_handler &
                  ("MIXING_STATE_RECOVER: deallocation (1b) failed")
             if (spin_separate) then
                deallocate(spoldin,spoldout,stat=status)
                if (status /= 0) call error_handler &
                     ("MIXING_STATE_RECOVER: deallocation (1'b) failed")
             endif
          endif
          idystr = max( iter_ch + 1 , iter_xc + 2 )
          dynamic_xc_mixing = .false.
       endif
       spin_pol(1) = zero
    endif

    if (exch_dyn(1) /= zero .and. spin_pol(1) /= zero) then
       if (dynamic_xc_mixing .and. n_spin > 1) then
          call readwriteblocked_read(xcoldin(:,2),th)
          call readwriteblocked_read(xcoldout(:,2),th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'xcoldin(:,2)'
             write(output_unit,'(4es20.13)')xcoldin(:,2)
             write(output_unit,'(  a     )')'xcoldout(:,2)'
             write(output_unit,'(4es20.13)')xcoldout(:,2)
          endif
       else
          call readwriteblocked_skipread(2*n_xc_stored,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'xcoldin(:,2)  skipped'
             write(output_unit,'(  a     )')'xcoldout(:,2) skipped'
          endif
       endif
    elseif (dynamic_xc_mixing .and. n_spin > 1) then
       if (output_main_scf) call write_to_output_units &
            ("MIXING_STATE_RECOVER: starting of dynamic mixing reset")
       if (iter_ch > idystr - 1) then
          deallocate(choldin,choldout,stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_STATE_RECOVER: deallocation (1c) failed")
          if (spin_separate) then
             deallocate(spoldin,spoldout,stat=status)
             if (status /= 0) call error_handler &
                  ("MIXING_STATE_RECOVER: deallocation (1'c) failed")
          endif
       endif
       deallocate(xcoldin,xcoldout,stat=status)
       if (status /= 0) call error_handler &
            ("MIXING_STATE_RECOVER: deallocation (2c) failed")
       idystr = max( iter_ch + 1 , iter_xc + 2 )
       dynamic_xc_mixing = .false.
    endif

  end subroutine mixing_state_recover

#ifdef WITH_SCFCONTROL
  subroutine mixing_read_scfcontrol(warning)
    ! purpose: read in the namelist 'mixing' from the file
    ! $TTFSSTART/scfcontrol at each scf-cycle giving a
    ! notice if parameters have changed.
    !------------ Declaration of formal parameters -------------
    logical, intent(in), optional :: warning ! if present and .false. warning
                                             ! messages are suppressed
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind)    :: chmix_old
    real(kind=r8_kind)    :: xcmix_old
    real(kind=r8_kind)    :: spmix_old
    integer(kind=i4_kind) :: idystr_old
    logical               :: discard_init_old
    logical               :: spin_separate_old
    logical               :: global_mixing_old
    !------------ Executable code ------------------------------

    ! save input variables
    chmix_old=chmix
    xcmix_old=xcmix
    spmix_old=spmix
    idystr_old=idystr
    discard_init_old=discard_init
    spin_separate_old=spin_separate
    global_mixing_old=global_mixing

    call mixing_read (init_beta_values=.false.)
    discard_init = discard_init_old

    if (present(warning)) then
       if (.not.warning) return
    endif

    call check_real("CHMIX",chmix,chmix_old)
    if (spin_separate) then
       call check_real("SPMIX",spmix,spmix_old)
    endif
    call check_real("XCMIX",xcmix,xcmix_old)
    call check_intg("IDYSTR",idystr,idystr_old)
    call check_flag("SPIN_SEPARATE",spin_separate,spin_separate_old)
    call check_flag("GLOBAL_MIXING",global_mixing,global_mixing_old)

    ! keep state of idystr and spin_separate until call of mixing_reset buffer
    idystr_new = idystr
    spin_separate_new = spin_separate
    idystr = idystr_old
    spin_separate = spin_separate_old

  contains

    subroutine check_intg(name,val,val_old)
      character(len=*)     , intent(in) :: name
      integer(kind=i4_kind), intent(in) :: val, val_old
      if ( val .ne. val_old ) then
         call write_to_output_units("mixing_read_scfcontrol: "// &
              name//" was altered. New value: ",inte=val)
         call write_to_trace_unit("mixing_read_scfcontrol: "// &
              name//" was altered. New value: ",inte=val)
      endif
    end subroutine check_intg

    subroutine check_real(name,val,val_old)
      character(len=*)     , intent(in) :: name
      real(kind=r8_kind)   , intent(in) :: val, val_old
      if ( val .ne. val_old ) then
         call write_to_output_units("mixing_read_scfcontrol: "// &
              name//" was altered. New value: ",re=val)
         call write_to_trace_unit("mixing_read_scfcontrol: "// &
              name//" was altered. New value: ",real=val)
      endif
    end subroutine check_real

    subroutine check_flag(name,val,val_old)
      character(len=*)     , intent(in) :: name
      logical              , intent(in) :: val, val_old
      character(len=5) :: val_str
      if ( val .neqv. val_old ) then
         if (val) then
            val_str = ' TRUE'
         else
            val_str = 'FALSE'
         endif
         call write_to_output_units("mixing_read_scfcontrol: "// &
              name//" was altered. New value: "//val_str)
         call write_to_trace_unit("mixing_read_scfcontrol: "// &
              name//" was altered. New value: "//val_str)
      endif
    end subroutine check_flag

  end subroutine mixing_read_scfcontrol
#endif

  subroutine mixing_reset_buffers(loop,update_scfcontrol)
    ! purpose: resets the mixing buffers according to the new values
    !          of spin_separate and idystr
    !------------ Modules ----------------------------------------
    use options_module, only : options_xcmode, xcmode_exchange_fit
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind), intent(in) :: loop
    logical, optional, intent(inout) :: update_scfcontrol
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    logical               :: xcfit, spin_separate_old
    integer(kind=i4_kind) :: iter_ch, iter_xc, status, idystr_old
    !------------ Executable code ------------------------------

    ! now set idystr and spin_separate according to the scfcontrol file
    spin_separate_old = spin_separate
    idystr_old = idystr
    spin_separate = spin_separate_new
    idystr = idystr_new

    ! note, that mixing_xc has already been called in the current SCF loop
    ! whereas mixing_ch has not yet been called in the current SCF loop
    iter_ch = loop
    iter_xc = loop - 1 ! see call of mixing_xc in build_xcfit

    xcfit = options_xcmode() == xcmode_exchange_fit

    if (spin_separate .and. .not. spin_separate_old) then
       ! spin fit buffer currently not activated
       idystr = max(idystr,iter_ch+1) ! => iter_ch <= idystr-1
    endif
    if (iter_ch <= idystr_old-1) then
       ! charge fit buffer not yet allocated
       idystr = max(idystr,iter_ch+1) ! => iter_ch <= idystr-1
    endif
    if (xcfit .and. iter_xc < idystr_old-1) then
       ! exchange fit buffer not yet allocated
       idystr = max(idystr,iter_xc+2) ! => iter_xc < idystr-1
    endif
    if (idystr /= idystr_new) then
       call write_to_output_units("mixing_reset_buffers: "// &
            "IDYSTR had to be altered. New value: ",inte=idystr)
       call write_to_trace_unit("mixing_reset_buffers: "// &
            "IDYSTR had to be altered. New value: ",inte=idystr)
       if (present(update_scfcontrol)) update_scfcontrol = .true.
    endif

    if (iter_ch <= idystr-1 .and. iter_ch > idystr_old-1) then
       deallocate(choldin,choldout,stat=status)
       if (status /= 0) call error_handler &
            ("MIXING_RESET_BUFFERS: deallocation (1) failed")
       if (spin_separate_old) then
          deallocate(spoldin,spoldout,stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_RESET_BUFFERS: deallocation (1') failed")
       endif
    endif
    if (xcfit .and. iter_xc < idystr_old-1 .and. iter_xc >= idystr_old-1) then
       deallocate(xcoldin,xcoldout,stat=status)
       if (status /= 0) call error_handler &
            ("MIXING_RESET_BUFFERS: deallocation (2) failed")
    endif
    if (.not.spin_separate .and. spin_separate_old) then
       ! if iter_ch <= idystr-1 spin buffers have already been deallocated
       if (iter_ch > idystr_old-1 .and. iter_ch > idystr-1) then
          deallocate(spoldin,spoldout,stat=status)
          if (status /= 0) call error_handler &
               ("MIXING_RESET_BUFFERS: deallocation (1'') failed")
       endif
    endif

  end subroutine mixing_reset_buffers

  function mixing_beta_ch()
    real(kind=r8_kind)  :: mixing_beta_ch
    !** End of interface ***************************************
    !------- executable code -------------------
    mixing_beta_ch = betachold
  end function mixing_beta_ch

  function mixing_beta_sp()
    real(kind=r8_kind)  :: mixing_beta_sp
    !** End of interface ***************************************
    !------- executable code -------------------
    mixing_beta_sp = betaspold
  end function mixing_beta_sp

  function mixing_beta_xc()
    real(kind=r8_kind)  :: mixing_beta_xc
    !** End of interface ***************************************
    !------- executable code -------------------
    mixing_beta_xc = betaxcold
  end function mixing_beta_xc

end module mixing_module

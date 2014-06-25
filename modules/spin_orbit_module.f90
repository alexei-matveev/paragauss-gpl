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
module spin_orbit_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains steering parameters of spin orbit variant
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: MM
  !  Date: 7/98
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
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public :: spin_orbit_read_input
  public :: spin_orbit_write_input
  public :: spin_orbit_input_bcast
  public :: WhatIs
  public :: is_on
  public :: max_alloc_reals

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of constants and variables ----

  ! ----------- Spin Orbit control parameters ---------------------------

  logical,public :: spin_orbit_polarized !!$,spin_orbit_no_pvxp

  integer(IK),parameter,public ::&
       op_SpinOrbit   = 1,&
       op_Polarized   = 2,&
       op_BackTrafo   = 3,&
       op_FitTrafo    = 4,&
       op_NoPVxP      = 5,&
       op_Eigensolver = 6,&
       op_BSOA        = 7, &
       op_RelFit      = 8, &
       op_FinNuc      = 9, &
       op_NoPFyP      =10
  integer(IK),parameter :: &
       NOp            =10 ! NumberOfOptions, the last above
  integer(IK)           :: ops(NOp) = -1

  integer(IK),parameter, public :: &
       SO_RELFIT    = 1, & ! default, compute all integrals
       SO_RELFIT_1C = 2    ! compute only 1C (same atom) integrals

  integer(IK), parameter          :: strlen=256
  character(len=strlen)           :: LEVEL    = "undef"
  character(len=strlen),parameter :: df_LEVEL = "none"
  character(len=strlen)           :: RunLevel = "undef"

  real(RK), parameter             :: const_speed_of_light=  137.03604_rk
  real(RK), public                :: speed_of_light =  const_speed_of_light
  real(RK), private               :: VC             =  1.0_rk ! (V/C) ratio
  real(RK), private               :: df_VC          =  1.0_rk
  real(RK), public                :: relfit_minexp  =       0.0_rk
  real(RK), public                :: relfit_maxexp  =  HUGE(1.0_rk)
  real(RK), public                :: relfit_minZ    =  -1.0_rk
  integer(IK),private             :: max_alloc_reals_ = -1
  ! default in fit_trafo_module

  ! ----------- Default values for input parameters -------------



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine fingerprint()
    use error_module
    use iounitadmin_module, only: write_to_output_units
    implicit none
    ! *** end of interface ***
    character(len=strlen) :: append
    character(len=strlen) :: buf
    integer(IK)           :: op

    RunLevel = " "

    if(is_on(op_SpinOrbit))then
       RunLevel = ccat(RunLevel,"SpinOrbit")
    else
       RunLevel = ccat(RunLevel,"none")
    endif

    if(is_on(op_NoPVxP))then
       RunLevel = ccat(RunLevel,",NoPVxP")
    endif

    select case(whatis(op_NoPFyP))
    case (1)
       RunLevel = ccat(RunLevel,",NoPFsP")
    case (2)
       RunLevel = ccat(RunLevel,",NoPFxP")
    case (3)
       RunLevel = ccat(RunLevel,",NoPFyP")
    case (0)
       ! default
    case default
       ABORT('no such case')
    end select

    if(is_on(op_BSOA))then
       RunLevel = ccat(RunLevel,",BSOA")
    endif

    if(is_on(op_Polarized))then
       RunLevel = ccat(RunLevel,",Polarized")
    endif

    if(is_on(op_BackTrafo))then
       select case(WhatIs(op_BackTrafo))
       case (2)
          RunLevel = ccat(RunLevel,",BackTrafo(VH)")
       case (0)
          ! do nothing
       case default
          ABORT('no such use of BackTrafo')
       end select
    endif

    if(is_on(op_FitTrafo))then
       RunLevel = ccat(RunLevel,",FitTrafo")

       select case(WhatIs(op_FitTrafo))
       case(1)
          append = ":DKH(1)"
       case(2)
          append = ":DKH(2)[=oDKH(2)]"
       case(3)
          append = ":hDKH(2)"
       case default
          call error("som/fingerprint: no such case")
       end select

       RunLevel = ccat(RunLevel,append)
    endif

    if(is_on(op_RelFit))then
       select case(whatis(op_RelFit))
       case (SO_RELFIT)
          RunLevel = ccat(RunLevel,",RelFit")
       case (SO_RELFIT_1C)
          RunLevel = ccat(RunLevel,",RelFit(1c)")
       case default
          ABORT('whatis(op_RelFit)')
       end select
    endif

    ! ====== Finite Nuclei ===========
    if(is_on(op_FinNuc))then
       op = whatis(op_FinNuc)
       if( op == 7 )then
          RunLevel = ccat(RunLevel,",FinNuc")
       else if( op == 6 )then
          RunLevel = ccat(RunLevel,",FinRel")
       else
          if( IAND(op,1) == 1 )then
             RunLevel = ccat(RunLevel,",FinNucOnly")
          endif
          if( IAND(op,2) == 2 )then
             RunLevel = ccat(RunLevel,",FinPVSP")
          endif
          if( IAND(op,4) == 4 )then
             RunLevel = ccat(RunLevel,",FinPVXP")
          endif
       endif
       ASSERT(op<=7)
    endif
    ! ================================

    if(.not.is_on(op_Eigensolver))then
       RunLevel = ccat(RunLevel,",EIS=EISPACK")
    endif

    call say ("")
    if (is_on (op_SpinOrbit)) then
       call say ("SpinOrbit Run")
    else
       call say ("Standard  Run")
    endif
    call say ("")
    call say (trim (RunLevel))
    if (VC /= df_VC) then
       call say ("SPEED OF LIGHT CHANGED!")
    endif
    write (buf, '("V/C = ",1PE16.8)') VC
    buf = adjustl (buf)
    call say (trim (buf))
    call say ("")

  contains

    subroutine say (phrase)
      use comm, only: comm_rank
      use iounitadmin_module, only: write_to_output_units
      implicit none
      character (len=*), intent (in) :: phrase
      ! *** end of interface ***

      if (comm_rank() == 0) then
         call write_to_output_units ("SO: " // phrase)
      endif
    end subroutine say

    function ccat(s1,s2) result(s)
      implicit none
      character(len=*),intent(in) :: s1,s2
      character(LEN=strlen)       :: s !<<<result
      ! *** end of interface ***

      s = " "
      if(len_trim(s1)>0) s = trim(s1)
      if(len_trim(s2)>0) s = trim(s) // trim(s2)
      s = trim(s)
    end function ccat
  end subroutine fingerprint

  subroutine set_defaults()
    implicit none
    ! *** end of interface ***

    LEVEL = df_LEVEL
    VC    = df_VC
  end subroutine set_defaults

  subroutine set_options()
    use error_module
    use strings
    use options_module, only: options_relativistic
    implicit none
    ! *** end of interface ***

    integer(IK) :: op

    if(spresent(LEVEL,"EIS=EISPACK"))then
       ops(op_Eigensolver) = 0
    else
       ops(op_Eigensolver) = 1
    endif

    if(spresent(LEVEL,"EIS=LAPACK"))then
       WARN('EIS=LAPACK is default anyway')
    endif

    if(wpresent(LEVEL,"SpinOrbit"))then
       if(spresent(LEVEL,"SpinOrbit:DKH(3)"))then
          ops(op_SpinOrbit) = 3
       elseif(spresent(LEVEL,"SpinOrbit:DKH(1)"))then
          ops(op_SpinOrbit) = 1
          print *,'DKHHHHHHHHHHHHHHH1'
       else
          ops(op_SpinOrbit) = 2
       endif
    else
       ops(op_SpinOrbit) = 0
    endif
    DPRINT 'ops(op_SpinOrbit)=',ops(op_SpinOrbit),whatis(op_SpinOrbit)

    if(wpresent(LEVEL,"NoPVxP"))then
       ops(op_NoPVxP) = 1
    else
       ops(op_NoPVxP) = 0
    endif

    ops(op_NoPFyP) = 0
    if(wpresent(LEVEL,"NoPFsP"))then
       ops(op_NoPFyP) = 1
    endif
    if(wpresent(LEVEL,"NoPFxP"))then
       ops(op_NoPFyP) = 2
    endif
    if(wpresent(LEVEL,"NoPFyP"))then
       ops(op_NoPFyP) = 3
    endif

    if(wpresent(LEVEL,"BSOA"))then
       ops(op_BSOA) = 1
    else
       ops(op_BSOA) = 0
    endif

    if(wpresent(LEVEL,"Polarized"))then
       ops(op_Polarized) = 1
    else
       ops(op_Polarized) = 0
    endif
    spin_orbit_polarized = is_on(op_Polarized)

    if(wpresent(LEVEL,"BackTrafo"))then
       if(spresent(LEVEL,"BackTrafo(VH)"))then
          ops(op_BackTrafo) = 2
          WARN('BackTrafo set to 2')
       else
          WARN('use FitTrafo instead of BackTrafo')
          ops(op_BackTrafo) = 0
       endif
    else
       ops(op_BackTrafo) = 0
    endif

    if(wpresent(LEVEL,"FitTrafo"))then
       ops(op_FitTrafo)     = 1

       if(spresent(LEVEL,"FitTrafo:DKH(2)"))then
          ops(op_FitTrafo)  = 2
       endif

       if(spresent(LEVEL,"FitTrafo:oDKH(2)"))then
          ops(op_FitTrafo)  = 2
       endif

       if(spresent(LEVEL,"FitTrafo:hDKH(2)"))then
          ops(op_FitTrafo)  = 3
       endif
    else
       ops(op_FitTrafo)     = 0
    end if

    if(wpresent(LEVEL,"RelFit"))then
       if(spresent(LEVEL,"RelFit(1c)"))then
          ops(op_RelFit)  = SO_RELFIT_1C
       else
          ops(op_RelFit)  = SO_RELFIT
       endif
    else
       ops(op_RelFit)     = 0
    endif

    if(.not.is_on(op_SpinOrbit))then
       if(    is_on(op_Polarized).or. is_on(op_BackTrafo).or.&
            & is_on(op_FitTrafo) )then
          ops(op_SpinOrbit) = 2
       else
          ops(op_SpinOrbit) = 0
       endif
    endif

    ! ====== Finite nuclei ===========
    op = 0 ! default overwritten
    if(wpresent(LEVEL,"FinNuc"))then
       op = 7         ! 1 + 2 + 4
    endif
    if(wpresent(LEVEL,"FinRel"))then
       op = 6         ! 2 + 4 = PVSP+PVXP
    endif
    if(wpresent(LEVEL,"FinNucOnly"))then
       op = IOR(op,1) ! 1 = V_{fin}(ua)
    endif
    if(wpresent(LEVEL,"FinPVSP"))then
       op = IOR(op,2) ! 2 = PV_{fin}SP(ua)
    endif
    if(wpresent(LEVEL,"FinPVXP"))then
       op = IOR(op,4) ! 4 = PV_{fin}XP(ua)
    endif
    if(.not.is_on(op_SpinOrbit))then
       ! unset (eventually) PVXP
       if( IAND(op,4) == 4 )then
          WARN('FinPVXP w/o SO?')
       endif
       op = IAND(op,3)
    endif
    if(.not.options_relativistic)then
       ! unset (eventually) PVSP
       if( IAND(op,2) == 2 )then
          WARN('FinPVSP w/o SR?')
       endif
       op = IAND(op,1)
    endif
    ops(op_FinNuc) = op
    ! ================================

!!$    write(*,nml=spin_orbit)
    speed_of_light = const_speed_of_light * VC

    call check_options()
  end subroutine set_options

  subroutine check_options()
    use error_module
    use options_module, only: options_spin_orbit
    implicit none
    ! *** end of interface ***

    if(LEVEL.ne.df_LEVEL)then
       if(is_on(op_SpinOrbit).neqv.options_spin_orbit)&
            & call error("som/check_options: options_spin_orbit & SpinOrbit conflict")
    endif

    ! reload op_SpinOrbit:
    if(options_spin_orbit)then
       if(.not.is_on(op_SpinOrbit))then
          ops(op_SpinOrbit) = 2
       endif
    else
       ops(op_SpinOrbit) = 0
    endif

    if(is_on(op_FitTrafo).and.is_on(op_BSOA))&
         & call error("som/check_options: FitTrafo & BSOA do not make sence together")

    if(is_on(op_NoPVxP).and.is_on(op_BSOA))&
         & call error("som/check_options: NoPVxP & BSOA do not make sence together")

    if( VC < 0.0 )&
         & call error("som/check_options: negative V/C ? abort...")
  end subroutine check_options

  function is_on(which) result(op)
    implicit none
    integer(IK),intent(in) :: which
    logical                :: op !<<< result
    ! *** end of intrface ***

    op = (WhatIs(which) /= 0)
  end function is_on

  function WhatIs(which) result(op)
    use error_module, only: MyID
    implicit none
    integer(IK),intent(in) :: which
    integer(IK)            :: op !<<< result
    ! *** end of intrface ***

!   print *, MyID, "whatis(", which, ")=", ops(which)

    if(ops(which).eq.-1)then
       print *, MyID, 'som/whatis: requested option ', which, ' was not set'
       print *, MyID, "ops=", ops
       ABORT('som/whatis: op is not yet set, see tty')
    endif

    op = ops(which)
  end function WhatIs

  !*************************************************************
  subroutine spin_orbit_read_input
    !  Purpose: read the input concerning spin_orbit
    !** End of interface *****************************************
    use input_module
    implicit none
    integer(IK) :: unit,status, max_alloc
    !
    ! only to read options in >>>
    ! inquire always through  WhatIs() or is_on()!
    !
    namelist /spin_orbit/&
         LEVEL,&
         VC, &
         relfit_minexp, &
         relfit_maxexp, &
         relfit_minZ  , &
         max_alloc ! let it be undocumented
    !------------ Executable code --------------------------------

    ! default values >>>
    call set_defaults()
    max_alloc = -1

    if ( input_line_is_namelist("spin_orbit") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=spin_orbit, iostat=status)
       if (status .gt. 0) call input_error( &
            "spin_orbit_read_input: namelist spin_orbit")
    endif
!   write (*,nml=spin_orbit)

    call set_options()
    if( max_alloc .ne. -1 )then
       max_alloc_reals_ = max_alloc
    endif

    call fingerprint()
  end subroutine spin_orbit_read_input

  subroutine spin_orbit_write_input(iounit)
    !
    ! Purpose: write the input concerning spin_orbit to output
    !
    use echo_input_module, only: start, real, strng, stop
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(IK), intent(in) :: iounit
    !** End of interface *****************************************

    call start("SPIN_ORBIT","SPIN_ORBIT_WRITE_INPUT", &
         iounit,operations_echo_input_level)
    call strng("LEVEL", LEVEL, df_LEVEL)
    call real( "VC   ",    VC, df_VC   )
    call real( "relfit_minexp ",    relfit_minexp,      0.0_rk  , format=4)
    call real( "relfit_maxexp ",    relfit_maxexp, HUGE(1.0_rk) , format=4)
    call real( "relfit_minZ   ",    relfit_minZ  ,     -1.0_rk  , format=4)
    call stop()
  end subroutine spin_orbit_write_input

  subroutine spin_orbit_input_bcast()
    ! Purpose : broadcasting the spin_orbit_input
    use comm, only: comm_bcast, comm_rank
    implicit none
    !** End of interface *****************************************

    call comm_bcast( LEVEL            )
    call comm_bcast( VC               )
    call comm_bcast( max_alloc_reals_ )
    call comm_bcast( relfit_minexp    )
    call comm_bcast( relfit_maxexp    )
    call comm_bcast( relfit_minZ      )

    if ( comm_rank() /= 0 ) call set_options()

  end subroutine spin_orbit_input_bcast

  function max_alloc_reals() result(n)
    implicit none
    integer(IK) :: n !<<< result
    ! *** end of interface ***

    n = max_alloc_reals_
  end function max_alloc_reals

  !*************************************************************
  !--------------- End of module ----------------------------------
end module spin_orbit_module

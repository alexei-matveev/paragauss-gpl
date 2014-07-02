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
module solv_charge_mixing_module
  !---------------------------------------------------------------
  !
  !  Purpose: Perform mixing of solvation surface charges
  !           during SCF procedure to improve convergence (?)
  !
  !
  !  Module called by: main_scf
  !
  !
  !  References: ...
  !
  !
  !  Author: Aleksey Shor
  !  Date: 03.04.2007
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
#include "def.h"
  use type_module ! type specification parameters
  use datatype
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  real(r8_kind),public :: Qs_mix=0.0_r8_kind
  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public read_solv_mix,write_solv_mix,solv_charge_mix,mix_charges,Q_dealloc

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----
  real(r8_kind),allocatable :: Q_in(:),Q_out(:)
  real(r8_kind)             :: beta_old

  logical,       parameter :: df_do_mixing=.false.
  real(r8_kind), parameter :: df_mixing=0.5_r8_kind
  logical,       parameter :: df_dynamic_mix=.false.

  logical       :: do_mixing
  real(r8_kind) :: mixing
  logical       :: dynamic_mix

  namelist /solv_charge_mixing/ do_mixing,mixing,dynamic_mix

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine read_solv_mix
    !  Purpose: read in mixing parameters from input
    !------------ Modules used ------------------- ---------------
    use input_module
    implicit none
    !------------ Declaration of formal parameters ---------------

    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: unit,status
    !------------ Executable code --------------------------------

    do_mixing=df_do_mixing
    mixing=df_mixing
    dynamic_mix=df_dynamic_mix

    if ( input_line_is_namelist("solv_charge_mixing") ) then
       unit = input_intermediate_unit()
       call input_read_to_intermediate
       read(unit, nml=solv_charge_mixing, iostat=status)
       if (status .ne. 0) call input_error( &
            "read_solv_mix: namelist solv_charge_mixing.")

       if(mixing <= 0.0_r8_kind) mixing=df_mixing
       if(mixing >= 1.0_r8_kind) mixing=1.0_r8_kind
    end if

  end subroutine read_solv_mix
  !*************************************************************

  !*************************************************************
  subroutine write_solv_mix(unit)
    !
    ! Purpose: write mixing parameters into input.out.
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full, real_format1
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(i4_kind), intent(in) :: unit
    !** End of interface *****************************************

    real_format1 = '("    ",a," = ",f7.3:" # ",a)'

    call start("SOLV_CHARGE_MIXING","WRITE_SOLV_MIX",unit,operations_echo_input_level)
    call flag("DO_MIXING  ",do_mixing  ,df_do_mixing  )
    call real("MIXING     ",mixing     ,df_mixing     ,1)
    call flag("DINAMIC_MIX",dynamic_mix,df_dynamic_mix  )
    call stop()

  end subroutine write_solv_mix
  !*************************************************************

  !*************************************************************
  subroutine solv_charge_mix(charges,charges_old,n_ch,iter,iter_start)
    !  Purpose: perform static (fixed) and dynamical mixing
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in)    :: n_ch,iter,iter_start
    real(r8_kind),    intent(inout) :: charges(n_ch),charges_old(n_ch)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind)              :: charges_buf(n_ch)
    real(r8_kind)              :: beta,d_mixing
    integer(i4_kind),parameter :: ndyn=1
    integer(i4_kind)           :: status
    !------------ Executable code --------------------------------


    if(.not.dynamic_mix) then
       if(iter==iter_start) then
          charges_old=charges
       else
          charges_buf=charges
          charges=charges*mixing+(1.0_r8_kind-mixing)*charges_old
          charges_old=charges_buf
       end if
       Qs_mix=mixing
    else
       if(iter==iter_start) then
          charges_old=charges
          Qs_mix=mixing
       else if(iter > iter_start .and. iter <= iter_start+ndyn) then
          if(iter==iter_start+ndyn) then
             allocate(Q_in(n_ch),Q_out(n_ch),stat=status)
             ASSERT(status==0)
             Q_in=charges_old
             Q_out=charges
          end if
          charges_buf=charges
          charges=charges*mixing+(1.0_r8_kind-mixing)*charges_old
          charges_old=charges_buf

          Qs_mix=mixing

          beta_old=mixing
       else if(iter > iter_start+ndyn) then
          beta=dyndamp(charges_old,charges,Q_in,Q_out)
          if(beta==-1.0_r8_kind) then
             beta=mixing
          end if
          d_mixing=(beta+beta_old)/2.0_r8_kind

          Qs_mix=d_mixing

          Q_in=charges_old
          Q_out=charges

          charges_buf=charges
          charges=charges*d_mixing+(1.0_r8_kind-d_mixing)*charges_old
          charges_old=charges_buf

          beta_old=d_mixing
       end if
    end if

  end subroutine solv_charge_mix
  !*************************************************************

  !*************************************************************
  subroutine Q_dealloc()
    !
    ! Does no communication, idempotent.
    !
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: status

    if (allocated (Q_in)) then
       deallocate (Q_in, Q_out, stat=status)
       ASSERT(status==0)
    end if
  end subroutine Q_dealloc
  !*************************************************************

  !*************************************************************
  function mix_charges()
    !  Purpose: pass information on mixing out of module
    !------------ Modules used -----------------------------------
    implicit none
    logical :: mix_charges
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    mix_charges=do_mixing

  end function mix_charges
  !*************************************************************

  !*************************************************************
  function DYNDAMP(AIN,AOUT,AOLDIN,AOLDOUT)
    ! Purpose: Calculation of the optimal mixing parameter.
    ! Taken from mixing_module
    implicit none
    real(r8_kind)                          :: DYNDAMP
    real(r8_kind),intent(IN),dimension(:)  :: AIN,AOUT,AOLDIN,AOLDOUT
    !** End of interface **************************************
    real (r8_kind), dimension (size (AIN)) :: NOM, DNOM, M
    integer (i4_kind), dimension (size (AIN)) :: H
    real(r8_kind)                      :: MMEAN,MSQUARE,SIGMA
    integer(i4_kind)                   :: NR
    real(r8_kind),parameter            :: zero=0.0_r8_kind
    real(r8_kind),parameter            :: one=1.0_r8_kind

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

    NR = sum (H)

    if(NR>1) then
       MMEAN=sum(M)/NR
       MSQUARE=sum(M*M)/NR
!!$   SIGMA=sqrt(MSQUARE-MMEAN*MMEAN)/NR
       SIGMA=sqrt(MSQUARE-MMEAN*MMEAN)
       DYNDAMP=one-MMEAN/(MMEAN-one)
       if(DYNDAMP>0.9_r8_kind) DYNDAMP=0.9_r8_kind
       if(DYNDAMP<0.1_r8_kind) DYNDAMP=0.1_r8_kind
    else
       DYNDAMP=-one
    endif

  end function DYNDAMP
  !--------------- End of module ----------------------------------
end module solv_charge_mixing_module

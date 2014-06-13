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
module interfaces
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
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

  use type_module, only:&
       & IK => i4_kind,&
       & RK => r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  public          ! by default, all names are public
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------


  interface
    subroutine chargefit(loop, coeff_dev, coulomb_dev)
      use type_module
      implicit none
      integer(i4_kind), intent(in) :: loop
      real(r8_kind), intent(out) :: coeff_dev, coulomb_dev
    end subroutine chargefit
  end interface

  interface
     subroutine main_integral(context)
       use type_module, only: IK=>i4_kind
       implicit none
       integer, intent(in) :: context
     end subroutine main_integral
  end interface

  ! values for context in main_integral(context):
  integer(IK), parameter :: &
       IMAST = 1,           & ! calling in master context
       ISLAV = 2,           & ! calling in slave  context
       IPARA = IMAST+ISLAV    ! calling in parallel context, bith on mast. and slave

  interface
     subroutine integral_trafo(mode)
       use type_module, only: IK=>i4_kind
       implicit none
       integer(IK),intent(in) :: mode
     end subroutine integral_trafo
  end interface
  ! and these are arguments it accepts:
  integer(IK),parameter ::&
       RELENRG = 1, &
       RELGRAD = 2, &
       RELSDER = 4

  interface
     subroutine calc_3center(na,la,nb,lb,equalb,equala,IMODE)
       use type_module ! type specification parameters
       integer(kind=i4_kind), optional,intent(in) :: equalb,equala
       integer(kind=i4_kind), intent(in)     :: na ! number of unique atom a
       integer(kind=i4_kind), intent(in)     :: nb ! number of unique atom b
       integer(kind=i4_kind), intent(in)     :: la ! angular momentum of unique atom a
       integer(kind=i4_kind), intent(in)     :: lb ! angular momentum of unique atom b
       integer(kind=i8_kind), optional,intent(in) :: IMODE ! bitmask
     end subroutine calc_3center
  end interface

  interface
     subroutine ll_calculate_grads(na,nb,equalb,la,lb,imode)
       use type_module ! type specification parameters
       integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
       integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
       integer(kind=i4_kind),intent(in) :: equalb ! number of equal atom b
       integer(kind=i4_kind),intent(in) :: la ! angular momentum of unique atom a
       integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
       integer(kind=i8_kind),intent(in) :: imode ! for control
     end subroutine ll_calculate_grads

     subroutine ls_calculate_grads(na,nb,equalb_,la_in,lb,imode)
       use type_module ! type specification parameters
       integer(kind=i4_kind),intent(in) :: na ! number of unique atom a
       integer(kind=i4_kind),intent(in) :: nb ! number of unique atom b
       integer(kind=i4_kind),intent(in) :: la_in ! angular momentum of unique atom a
       integer(kind=i4_kind),intent(in) :: equalb_ ! number of equal atom b
       integer(kind=i4_kind),intent(in) :: lb ! angular momentum of unique atom b
       integer(kind=i8_kind),intent(in) :: imode ! for control
     end subroutine ls_calculate_grads

     subroutine ss_calculate_grads(i_ea1,i_ea2,imode)
       use type_module ! type specification parameters
       integer(kind=i4_kind),intent(in)  :: i_ea1,i_ea2
       integer(kind=i8_kind),intent(in)  :: imode ! for control
     end subroutine ss_calculate_grads
  end interface

  interface
     subroutine potential_calculate(task)
       character(len=*), intent(in) :: task
     end subroutine potential_calculate

     subroutine grad_solv_calculate(task)
       character(len=*), intent(in) :: task
     end subroutine grad_solv_calculate
  end interface

  interface
     subroutine main_molmech(job,iwork,geo_steps)
       use type_module
       implicit none
       integer(i4_kind)           :: job
       integer(i4_kind), optional :: iwork,geo_steps
     end subroutine main_molmech
  end interface

  interface
     subroutine calc_cpks_gvec(ilower, iupper)
       use type_module
       implicit none
       integer(i4_kind), intent(in) :: ilower, iupper
     end subroutine calc_cpks_gvec
  end interface

  interface
     subroutine calc_cpks_h1imp(ilower, iupper)
       use type_module
       implicit none
       integer(i4_kind), intent(in) :: ilower, iupper
     end subroutine calc_cpks_h1imp
  end interface


  !------------ public functions and subroutines ------------------

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------

  !--------------- End of module ----------------------------------
end module interfaces

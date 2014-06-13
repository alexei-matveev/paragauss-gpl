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
module mat_charge_module
  !---------------------------------------------------------------
  !-------------- Module specification ---------------------------
  !---------------------------------------------------------------
  !
  !  Purpose: Provide the variable MAT_CHARGE ( [F_k|F_l])
  !           (Charge overlap matrix)
  !           and a routine to read it in from an
  !           external file
  !
  !
  !  Module called by: energy_calc_module/energ_coul_2z
  !                    chargefit and ...
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   7/7/97
  ! Description: Packed storage mode for mat_charge introduced
  !              ij=0; DO i=1,n_ch; DO j=1,i; ij=ij+1; ij=>(i,j); ENDDO; ENDDO
  !              Precalculation of mat_charge decomposition introduced
  !              Precalculation of dual_charge_norm introduced
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: Extended to the requirements of the potential extended
  !              model density approch:
  !              mat_exchange (analytical) introduced
  !              mat_xc_ch introduced
  !              matinv_exchange (numerical) introduced
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !------------ public functions and subroutines ------------------

  !------------ Modules used --------------------------------------
#include "def.h"
  use type_module ! type specification parameters
  use fit_coeff_module, only: fit,get_fit,charge_norm
  use options_module, only: options_perturbation_theory, options_xcmode, &
                            xcmode_model_density, xcmode_extended_mda
  implicit none
  private
  save

  !------------ Declaration of public constants and variables -----
  real(kind=r8_kind), allocatable, public, protected :: mat_charge(:)
  real(kind=r8_kind), allocatable, public, protected :: matinv_charge(:)
  real(kind=r8_kind), allocatable, public, protected :: dual_charge_norm(:)
  real(kind=r8_kind), allocatable, public, protected :: mat_exchange(:)
  real(kind=r8_kind), allocatable, public            :: matinv_exchange(:)
  real(kind=r8_kind), allocatable, public, protected :: mat_xc_ch(:, :)

  !------------ Declaration of subroutines used -------------------
  public :: mat_charge_set
  public :: read_mat_charge,free_mat_charge


!================================================================
! End of public interface of module
!================================================================


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  subroutine mat_charge_set(overlap)
    !
    ! Called from  int_send_2cff_module when charge  overlap matrix is
    ! ready.  Executed on master  only, thus, mat_charge is only valid
    ! there. Symmetric  charge overlap matrix is packed  into an array
    ! of length n*(n+1)/2 with n being the number of fit functions.
    !
    implicit none
    real(r8_kind), intent(in) :: overlap(:) ! n*(n+1)/2
    ! *** end of interface ***

    integer :: memstat

    allocate(mat_charge(size(overlap)), stat=memstat)
    ASSERT(memstat==0)

    mat_charge = overlap

    ! FIXME: provide charge_norm and do factorizaton right here.
  end subroutine mat_charge_set

  subroutine read_mat_charge()
    !
    ! FIXME: description outdated. Read  in MAT_CHARGE from file.  The
    ! file  is  named  :  $TTFSSTART/mat_charge.dat and  is  currently
    ! produced by the old lcgto.  The  routine is now only used by the
    ! master and since MAT_CHARGE  is used for the energy calculation.
    ! This  in turn  is only  done by  the master.   Nevertheless this
    ! routine  could  as  well  be  called by  the  slaves  under  the
    ! assumption, that  the slave processor possesses  the above file.
    ! FN, 10/95
    !
    ! If perturbation theory is not  turned on or if the model density
    ! approach is used, then:
    !
    ! -- Precalculation of mat_charge decomposition
    ! -- Precalculation of dual_charge_norm
    !
    !  UB, 7/97
    !
    use iounitadmin_module, only: openget_iounit, returnclose_iounit
    use filename_module, only: tmpfile
    use linsys_module, only: decompose_linsys
    use diis_fock_module, only: diis_pertstop
    implicit none
    !** End of interface *****************************************

    type(fit) :: n_fit
    real(r8_kind), parameter :: half = 0.5_r8_kind
    integer(i4_kind), parameter :: IO_len = 1024 ! rec_len = IO_len + 1
    integer(i4_kind) :: i, nmat, imax, io_in, u_g, status
    real(r8_kind), allocatable  :: mat_mixed(:)
    integer(i4_kind) :: alloc_stat

    call get_fit(n_fit)

    ! If   by   some   error   this  routine   was   called   although
    ! BOUNDS_CH%ITEM_ARR was 0 for this processor, print out a wraning
    ! and proceed ...
    if (n_fit%n_ch == 0) then
        WARN('charge basis empty')
    endif

    !
    ! The   last    two   lines   (options_perturbation_theory()   and
    ! diis_pertstop() belong  together).  It is now  possible that the
    ! perturbation_theory flag  will be changed during the  SCF run by
    ! changing it  from true to false  by the DIIS  routine.  There is
    ! then need  of the  matrices allocated here  as this  function is
    ! only called ones.  It is  important that it knows then that DIIS
    ! may do  this, therefore diis_pertstop  which gives the  value of
    ! the DIIS variable pert_stop (perturbation theory stop).
    !
    if (options_xcmode() == xcmode_model_density .or. &
        options_xcmode() == xcmode_extended_mda .or. &
       .not.options_perturbation_theory() .or. &
       diis_pertstop()) then

       !
       ! The  CHARGE_NORM must  be available  to decompose  the charge
       ! overlap matrix in advance:
       !
       ASSERT(allocated(charge_norm))
       ASSERT(allocated(mat_charge))

       !
       ! I guess decompose_linsys does an in-place factorization:
       !
       allocate(matinv_charge(size(mat_charge)), STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       matinv_charge = mat_charge

       ! note that size of dual_charge_norm is n_ch + 1 !
       allocate(dual_charge_norm(n_fit%n_ch + 1), STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       call decompose_linsys(n_fit%n_ch, matinv_charge, charge_norm, &
                             dual_charge_norm)
    endif

    ! ---------------------------------------------------------

    ! if MAT_CHARGE (and MATINV_CHARGE) is longer than the maximum
    ! block size, read it in in pieces of size rec_len = IO_LEN+1
    ! Obviously the old LCGTO has a different max. Record length
    ! and thus we have to set a special IO_len (grrmph) --------

    if (options_xcmode() == xcmode_extended_mda) then

       nmat = (n_fit%n_xc*(n_fit%n_xc+1))/2

       allocate(mat_exchange(nmat),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       allocate(matinv_exchange(nmat),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       allocate(mat_xc_ch(n_fit%n_xc,n_fit%n_ch),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

!      << first read in the analytical exchange overlap matrix >>
       io_in = openget_iounit(form='unformatted',status='old',&
            file=trim(tmpfile('mat_exchange.dat')))

       imax = nmat/(IO_len+1)
       u_g=1
       do i=1,imax
          read(io_in,err=98) mat_exchange(u_g:u_g+IO_len)
          u_g=u_g+IO_len+1
       enddo
       read(io_in,err=98) mat_exchange(u_g:nmat)

       call returnclose_iounit(io_in)

!      << then read in the mixed overlap matrix >>
       nmat = n_fit%n_xc * n_fit%n_ch

       allocate( mat_mixed(nmat), stat=status)
       if (status /= 0) call error_handler( &
            "read_mat_charge : allocation of mat_mixed failed")

       io_in = openget_iounit(form='unformatted',status='old',&
            file=trim(tmpfile('mat_mixed.dat')))

       imax = nmat/(IO_len+1)
       u_g=1
       do i=1,imax
          read(io_in,err=97) mat_mixed(u_g:u_g+IO_len)
          u_g=u_g+IO_len+1
       enddo
       read(io_in,err=97) mat_mixed(u_g:nmat)

       call returnclose_iounit(io_in)

       u_g=1
       do i=1,n_fit%n_ch
          mat_xc_ch(:,i) = mat_mixed(u_g:u_g+n_fit%n_xc-1)
          u_g=u_g+n_fit%n_xc
       end do

       deallocate( mat_mixed, stat=status)
       if (status /= 0) call error_handler( &
            "read_mat_charge : deallocation of mat_mixed failed")

    endif

    return

98  call error_handler("read_mat_charge: reading mat_exchange.dat failed")
97  call error_handler("read_mat_charge: reading mat_mixed.dat failed")
  end subroutine read_mat_charge

  subroutine free_mat_charge()
    !  Purpose: deallocates MAT_CHARGE
    !           This is written as a separate routine for
    !           reasons of 'better' programming style.
    !           FN 10/95
    !
    ! Is called from a parallel context, but does something
    ! meaningfull only on master.
    !
    ! Idempotent.
    !
    implicit none
    !** End of interface *****************************************

    if(allocated(mat_charge))then
      deallocate(mat_charge)
    endif

    if(allocated(matinv_charge))then
      deallocate(matinv_charge, dual_charge_norm)
    endif

    if(allocated(mat_exchange))then
      deallocate(mat_exchange, matinv_exchange, mat_xc_ch)
    endif
  end subroutine free_mat_charge
  !*************************************************************

!--------------- End of module ----------------------------------
end module mat_charge_module

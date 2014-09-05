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
module se_timefunction_module
  !-------------------------------------------------------------------
  !
  !  Purpose: provides functions to read the coefficients of a polynomial
  !           cost function from a standardized ASCII file.
  !
  !
  !  Module called by: se_scheduling_module
  !
  !
  !  Author: Martin Roderus
  !  Date: 04/2010
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

  use type_module ! type specification parameters
  !use iounitadmin_module, only: openget_iounit, returnclose_iounit
  implicit none
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------
  type, public :: se_timefunction_type
    integer(kind=i4_kind)           :: procCount
    real(kind=r8_kind), allocatable :: coefficients(:)
    real(kind=r8_kind)              :: scaleFactor
    integer(kind=i4_kind)           :: polyDegr
    integer(kind=i4_kind)           :: minSize
  end type se_timefunction_type



  !------------ public functions and subroutines ---------------------
  public :: se_timefunction_init

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------
  integer(kind=i4_kind), parameter :: fileUnit = 9898
  integer(kind=i4_kind), parameter :: dummyLen = 200



  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine se_timefunction_init(timeFunctions)
    !
    !  Purpose: initialize 'timeFunctions' of type se_timefunction_type with
    !  the coefficients of a cost function. In case a filename (with absolute
    !  path) is given, the coefficients are read from an external ascii file
    !  with defined syntax. Otherwise 'timeFunctions' is initialized with pre-
    !  defined values.
    !
    implicit none
    type(se_timefunction_type), allocatable, intent(out) :: timeFunctions(:)
    ! *** end of interface ***

    character(len=0) :: filename = '' ! FIXME: use hardwired coeefficients
    integer :: iostatus

    integer :: nFunctions, polyDegr, ios

    if( len(filename) .gt. 0 ) then
      call getfunctioninfo( filename, nFunctions, polyDegr, ios )
      if( ios .eq. 0 ) then
        call timefunction_allocate( timeFunctions, nFunctions, polyDegr )
        call readcoefficients( filename, timeFunctions, ios )
      end if
      iostatus = ios
    end if

!... <<START>> this section is automatically generated
    if( len(filename) .eq. 0   .or.   ios .ne. 0 ) then
      call timefunction_allocate( timeFunctions, 12, 3 )

      timeFunctions(1)%procCount = 1
      timeFunctions(1)%coefficients(1) = -1.064553E-01
      timeFunctions(1)%coefficients(2) = +1.729892E-03
      timeFunctions(1)%coefficients(3) = -2.868018E-06
      timeFunctions(1)%coefficients(4) = +3.477505E-09
      timeFunctions(1)%scaleFactor = 1.0
      timeFunctions(1)%polyDegr = 3
      timeFunctions(1)%minSize = getMinSize( 1 )

      timeFunctions(2)%procCount = 2
      timeFunctions(2)%coefficients(1) = +7.729889E-01
      timeFunctions(2)%coefficients(2) = -2.245244E-03
      timeFunctions(2)%coefficients(3) = +1.893119E-06
      timeFunctions(2)%coefficients(4) = +9.138625E-10
      timeFunctions(2)%scaleFactor = 1.0
      timeFunctions(2)%polyDegr = 3
      timeFunctions(2)%minSize = getMinSize( 2 )

      timeFunctions(3)%procCount = 3
      timeFunctions(3)%coefficients(1) = +8.430134E-03
      timeFunctions(3)%coefficients(2) = +4.603614E-04
      timeFunctions(3)%coefficients(3) = -6.308973E-07
      timeFunctions(3)%coefficients(4) = +1.195275E-09
      timeFunctions(3)%scaleFactor = 1.0
      timeFunctions(3)%polyDegr = 3
      timeFunctions(3)%minSize = getMinSize( 3 )

      timeFunctions(4)%procCount = 4
      timeFunctions(4)%coefficients(1) = -3.154133E-01
      timeFunctions(4)%coefficients(2) = +1.368635E-03
      timeFunctions(4)%coefficients(3) = -9.459603E-07
      timeFunctions(4)%coefficients(4) = +7.612922E-10
      timeFunctions(4)%scaleFactor = 1.0
      timeFunctions(4)%polyDegr = 3
      timeFunctions(4)%minSize = getMinSize( 4 )

      timeFunctions(5)%procCount = 5
      timeFunctions(5)%coefficients(1) = -2.816272E-01
      timeFunctions(5)%coefficients(2) = +1.106293E-03
      timeFunctions(5)%coefficients(3) = -5.797010E-07
      timeFunctions(5)%coefficients(4) = +5.520187E-10
      timeFunctions(5)%scaleFactor = 1.0
      timeFunctions(5)%polyDegr = 3
      timeFunctions(5)%minSize = getMinSize( 5 )

      timeFunctions(6)%procCount = 6
      timeFunctions(6)%coefficients(1) = -1.205034E+00
      timeFunctions(6)%coefficients(2) = +3.945849E-03
      timeFunctions(6)%coefficients(3) = -3.123970E-06
      timeFunctions(6)%coefficients(4) = +1.179315E-09
      timeFunctions(6)%scaleFactor = 1.0
      timeFunctions(6)%polyDegr = 3
      timeFunctions(6)%minSize = getMinSize( 6 )

      timeFunctions(7)%procCount = 7
      timeFunctions(7)%coefficients(1) = +1.171662E-01
      timeFunctions(7)%coefficients(2) = +1.624827E-04
      timeFunctions(7)%coefficients(3) = +1.525649E-07
      timeFunctions(7)%coefficients(4) = +2.780782E-10
      timeFunctions(7)%scaleFactor = 1.0
      timeFunctions(7)%polyDegr = 3
      timeFunctions(7)%minSize = getMinSize( 7 )

      timeFunctions(8)%procCount = 8
      timeFunctions(8)%coefficients(1) = -8.191168E-02
      timeFunctions(8)%coefficients(2) = +3.340531E-04
      timeFunctions(8)%coefficients(3) = +1.797315E-07
      timeFunctions(8)%coefficients(4) = +1.811571E-10
      timeFunctions(8)%scaleFactor = 1.0
      timeFunctions(8)%polyDegr = 3
      timeFunctions(8)%minSize = getMinSize( 8 )

      timeFunctions(9)%procCount = 9
      timeFunctions(9)%coefficients(1) = -1.111392E-02
      timeFunctions(9)%coefficients(2) = +1.591703E-04
      timeFunctions(9)%coefficients(3) = +2.306133E-07
      timeFunctions(9)%coefficients(4) = +1.395544E-10
      timeFunctions(9)%scaleFactor = 1.0
      timeFunctions(9)%polyDegr = 3
      timeFunctions(9)%minSize = getMinSize( 9 )

      timeFunctions(10)%procCount = 16
      timeFunctions(10)%coefficients(1) = +1.438538E-02
      timeFunctions(10)%coefficients(2) = +6.021485E-05
      timeFunctions(10)%coefficients(3) = +2.514895E-07
      timeFunctions(10)%coefficients(4) = +5.872729E-11
      timeFunctions(10)%scaleFactor = 1.0
      timeFunctions(10)%polyDegr = 3
      timeFunctions(10)%minSize = getMinSize( 16 )

      timeFunctions(11)%procCount = 25
      timeFunctions(11)%coefficients(1) = +1.348232E-02
      timeFunctions(11)%coefficients(2) = +1.017546E-04
      timeFunctions(11)%coefficients(3) = +1.831714E-07
      timeFunctions(11)%coefficients(4) = +4.214535E-11
      timeFunctions(11)%scaleFactor = 1.0
      timeFunctions(11)%polyDegr = 3
      timeFunctions(11)%minSize = getMinSize( 25 )

      timeFunctions(12)%procCount = 36
      timeFunctions(12)%coefficients(1) = +2.629501E-02
      timeFunctions(12)%coefficients(2) = +8.394168E-05
      timeFunctions(12)%coefficients(3) = +1.715824E-07
      timeFunctions(12)%coefficients(4) = +2.408075E-11
      timeFunctions(12)%scaleFactor = 1.0
      timeFunctions(12)%polyDegr = 3
      timeFunctions(12)%minSize = getMinSize( 36 )

    end if
!... <<END>> this section is automatically generated

  end subroutine se_timefunction_init


  subroutine timefunction_allocate( timeFunctions, nFunctions, polyDegr )
    !  Purpose: allocate timeFunctions
    implicit none
    type(se_timefunction_type), allocatable, intent(out)  :: timeFunctions(:)
    integer, intent(in) :: nFunctions
    integer, intent(in) :: polyDegr
    !** End of interface *****************************************

    integer :: i

    allocate( timeFunctions(nFunctions) )

    do i=1, nFunctions
      allocate( timeFunctions(i)%coefficients(polyDegr+1) )
    end do
  end subroutine timefunction_allocate

  subroutine getfunctioninfo( filename, nFunctions, polyDegr, iostatus )
    !  Purpose: read header from file and retrieve the following information:
    !           the total number of functions and the polynomial degree of a function.
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*),      intent(in)  :: filename
    integer(kind=i4_kind), intent(out) :: nFunctions
    integer(kind=i4_kind), intent(out) :: polyDegr
    integer(kind=i4_kind), intent(out) :: iostatus
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)   :: ios
    character(len=dummyLen) :: stringDummy
    !------------ Executable code ------------------------------------

    !fileUnit = openget_iounit( trim(filename), status='old' ) !!FIXME: in paragauss ersetzen
    open(unit=fileUnit, file=trim(filename), status='old', iostat=iostatus)!!!FIXME: zeile loeschen
    if( iostatus .ne. 0 ) then    !# I/O error handling FIXME: Stderr von Paragauss???
      write(*,*) '[SchedEig] Error opening file: ', trim(filename)
      return
    end if

    ios = 0
    do while (ios .eq. 0)
      read( unit=fileUnit, fmt='(A)', iostat=ios ) stringDummy
      if( stringDummy(:2)=='N:' ) exit
    end do

    read( stringDummy(3:), fmt='(I5)' ) nFunctions

    rewind( fileUnit )

    ios = 0
    do while (ios .eq. 0)
      read( unit=fileUnit, fmt='(A)', iostat=ios ) stringDummy
      if( stringDummy(:4)=='DEG:' ) exit
    end do

    read( stringDummy(5:), fmt='(I5)' ) polyDegr

    rewind( fileUnit )
    close( fileUnit )
    !call returnclose_iounit( fileUnit )

  end subroutine getfunctioninfo
  !*************************************************************


  !*************************************************************
  subroutine readcoefficients( filename, timeFunctions, iostatus )
    !  Purpose: read the coefficients from the file
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*),           intent(in)    :: filename
    type(se_timefunction_type), intent(inout) :: timeFunctions(:)
    integer(kind=i4_kind),      intent(out)   :: iostatus
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)   :: i, j, ios
    character(len=dummyLen) :: stringDummy
    !------------ Executable code ------------------------------------

    open(unit=fileUnit, file=trim(filename), status='old', iostat=iostatus)
    if( iostatus .ne. 0 ) then    ! I/O error handling FIXME: Stderr von Paragauss???
      write(*,*) '[SchedEig] Error opening file: ', trim(filename)
      return
    end if

    ios = 0; i=1
    do  while (ios .eq. 0)
      read( unit=fileUnit, fmt='(A)', iostat=ios ) stringDummy
      if( stringDummy(:2)=='P:' ) then
        j=1
        read( stringDummy(3:), fmt='(I5)' ) timeFunctions(i)%procCount
        read( unit=fileUnit, fmt='(A3)', iostat=ios, advance='no' ) stringDummy
        do j=1, timeFunctions(i)%polyDegr + 1
          read( unit=fileUnit, fmt='(D13.6)', advance='no' ) timeFunctions(i)%coefficients(j)
          if(   .not.   j .eq. timeFunctions(i)%polyDegr + 1 ) then
            read( unit=fileUnit, fmt='(A4)', advance='no' ) stringDummy
          end if
        end do
        timeFunctions(i)%scaleFactor = 1.
        timeFunctions(i)%minSize = getMinSize( timeFunctions(i)%procCount )

        i = i+1
      end if
    end do

    close( fileUnit )

  end subroutine readcoefficients
  !*************************************************************


  !*************************************************************
  integer function getMinSize( procCount )
    !  Purpose: return a defined minimum size of a matrix which a
    !           given number of processors can compute
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(in) :: procCount
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind), parameter :: blocksize=64
    !------------ Executable code ------------------------------------


    if( .not. procCount .eq. 1 ) then
      if( procCount .le. 9 ) then
        getMinSize = procCount*blocksize
      else
        getMinSize = int (sqrt (real (procCount)) * blocksize)
      end if
    else
      getMinSize = 20
    end if

  end function getMinSize
  !*************************************************************


  !--------------- End of module -------------------------------------
end module se_timefunction_module

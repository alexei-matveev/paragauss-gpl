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
!=========================================================================================
! Public interface of module
!=========================================================================================
module thermodyn_prop_module
  !---------------------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------------------
  !== Interrupt of public interface of module ============================================
  !---------------------------------------------------------------------------------------
  ! Modifications
  !---------------------------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !---------------------------------------------------------------------------------------
  use type_module                                          ! type specification parameters

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ========================================


  !------------ public functions and subroutines ------------------
  public                        :: thermodynamic_properties

  !------------ Declaration of constants and variables -----------------------------------
  !*** For the declaration of mathematical and fundamental constants as well as for the **
  !*** declaration of the conversion factors see "/modules/constants.f90" ****************

  !=======================================================================================
  ! End of public interface of module
  !=======================================================================================

  !---------------------------------------------------------------------------------------
  !------------ Subroutines --------------------------------------------------------------
  contains


  !***************************************************************************************
  !FIXME: add later on electronic components?
  subroutine thermodynamic_properties( output_unit, fre, Mom_inert, total_mass           &
                                     , tot_energy, sigma, T, P, deltaT, NTpoints)

  !---------------------------------------------------------------------------------------
  !                                                                                      !
  !  Purpose: Calculation of translational, electronic, vibrational and rotational       !
  !           contributions to internal thermal energy "U" and entropy "S"               !
  !                                                                                      !
  !                                                                                      !
  !  Module called by: "frequency_main subroutine"  in file  "frequency_module.f90"      !
  !                                                                                      !
  !                                                                                      !
  !  References: Florian Schlosser: PhD-Thesis.                                          !
  !                                                                                      !
  !              D.A. McQuarrie, J.D. Simon: Molecular Thermodynamics,                   !
  !              University Science Books, Sausalito, California                         !
  !                                                                                      !
  !                                                                                      !
  !  Author: TMS                                                                         !
  !  Date: February 2009                                                                 !
  !                                                                                      !
  !                                                                                      !
  !---------------------------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: TMS
  ! Date:   June 2009
  ! Description: Extension to temperature intervalls
  !
  ! Modification (Please copy before editing)
  ! Author:
  ! Date:
  ! Description:
  !
  !----------------------------------------------------------------
    use constants , only             : ZERO, HALF ,ONE , THREE_HALF, TWO, PI, FIVE_HALF &
                                     , THREE, TWO_PI, THOUSAND &
                                     , R_gas &                        ! Ideal gas constant
                                     , k_B &                          ! Boltzmann constant
                                     , c_speedoflight &            ! vacuum speed of light
                                     , h_Planck &                        ! Planck constant
                                     , hbar_Planck &               ! Planck constant / 2pi
                                     , u_Atom_Mass &                ! atomic mass constant
                                     , E_h2Jmol                   ! Conversion Eh -> J/mol
    implicit none
    !------------ Declaration of formal parameters ---------------------------------------
    integer(kind=i4_kind),intent(in)        :: output_unit &
                                             , NTpoints &    ! Number of Temperature steps
                                             , sigma ! symmetry index: number of identical
                                                      ! positions accessible via rotations
    real(kind=r8_kind),intent(in)           :: fre(:)  ! vector containing all frequencies
                                                             ! of the normal modes in [Hz]
    real(r8_kind),dimension(3),intent(in)   :: Mom_inert          ! Vector with Moments of
                                                                     ! inertia in [kg*m^2]

    real(kind=r8_kind),intent(in)           :: total_mass  &         ! molecules mass [kg]
                                             , tot_energy  &            ! tot_energy [E_h]
                                             , T           &             ! temperature [K]
                                             , P                           ! pressure [Pa]

    real(kind=r8_kind),intent(in)           :: deltaT               ! Temperature step [K]

    real(kind=r8_kind), parameter           :: liC = 5.0E-51_r8_kind ! Linearity Criterion

    !** End of interface *****************************************************************
    !------------ Declaration of local variables -----------------------------------------
    integer(kind=i4_kind)                   :: K = 0 &                   ! summation index
                                             , J = 0 &                   ! auxiliary index
                                             , num_pos_freq = 0 ! number of positive modes
    real(r8_kind),dimension(NTpoints)       :: tmpT   ! Temperature points
    real(kind=r8_kind),dimension(NTpoints)  :: U            &    ! internal thermal energy
                                             , S            &                    ! entropy
                                             , entH         &                   ! enthalpy
                                             , Gfe          &          ! Gibbs free energy
                                             , Gfe_tot                           ! + e_sum
    real(kind=r8_kind),dimension(NTpoints)  :: q_t          &    ! Translational variables
                                             , E_t          &
                                             , S_t
    real(kind=r8_kind),dimension(NTpoints)  :: E_e          &    ! electr motion variables
                                             , S_e
    real(kind=r8_kind),dimension(NTpoints)  :: E_v          &      ! Vibrational variables
                                             , S_v
    real(kind=r8_kind)                      :: ZPE_v
    real(kind=r8_kind),allocatable,dimension(:) :: T_v &        ! Vibrational Temperatures
                                                 , freq             ! positive frequencies
    real(kind=r8_kind) ,dimension(NTpoints) :: q_r          &       ! Rotational variables
                                             , E_r          &
                                             , S_r
    real(kind=r8_kind), dimension(3)        :: mom_I        & ! ordered moments of inertia
                                             , T_r                ! Rotational Temperature

    logical                                 :: lin_mol = .FALSE.

    character(9)                            :: rot_eps_shape
    !------------ Executable code --------------------------------------------------------
    !*** initializing variables **********************************************************
    tmpT      = ZERO
    U         = ZERO
    S         = ZERO
    entH      = ZERO
    Gfe       = ZERO
    Gfe_tot   = ZERO
    q_t       = ZERO
    E_t       = ZERO
    S_t       = ZERO
    E_e       = ZERO
    S_e       = ZERO
    E_v       = ZERO
    S_v       = ZERO
    ZPE_v     = ZERO
    q_r       = ZERO
    E_r       = ZERO
    S_r       = ZERO
    mom_I     = ZERO
    T_r       = ZERO
    !*** Define temperature points *******************************************************
    DO K = 1,NTpoints
      tmpT(K) = T + real(K-1,r8_kind)*deltaT
    END DO
    !*** count positive frequencies ******************************************************
    num_pos_freq = 0
    DO K = 1,size(fre)
      IF (fre(K) > 0) THEN
        num_pos_freq = num_pos_freq + 1
      END IF
    END DO
    !*** preparations ********************************************************************
    allocate(freq(num_pos_freq),T_v(num_pos_freq))
    !*** complete initializing variables *************************************************
    T_v       = ZERO
    freq      = ZERO
    !*** write positive frequencies into "freq"-vector ***********************************
    J = 0
    DO K = 1,size(fre)
      IF (fre(K) > 0) THEN
        J = J + 1
        freq(J) = fre(K)
      END IF
    END DO
    !*** Sort eigenvalues ****************************************************************
    mom_I = sort3(Mom_inert)

    !*** Write initial output ************************************************************
    write(output_unit,*)
    write(output_unit,01)
    01 format('  +---------------------------------------------------------------------------------------+')
    write(output_unit,02)
    02 format('  |',87X,'|')
    write(output_unit,03)
    03 format('  |     Thermodynamic properties section:',49X,'|')
    write(output_unit,02)
    write(output_unit,04) T, P
    04 format('  |',5X,'T   = ',F8.3,' K',21X,'p = ',F9.1,' Pa',29X,'|')
    do K = 2,NTpoints
      write(output_unit,0401) tmpT(K)
      0401 format('  |',11X,F8.3,' K',66X,'|')
    enddo
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*************************************************************************************
    !    calculation of translational contributions using the expressions                !
    !                                                                                    !
    !         q_t = (2*pi*m*k_B*T/h^2)^(3/2)*(k_B*T/P)                                   !
    !                                                                                    !
    !         E_t = (3/2)*R*T                                                            !
    !                                                                                    !
    !         S_t = R*(ln(q_t)+5/2)                                                      !
    !                                                                                    !
    !*************************************************************************************
    !*** Calculations of translational contributions !************************************
    q_t = (sqrt(TWO_PI*total_mass)/h_Planck)**THREE*(k_B*tmpT)**FIVE_HALF/P

    E_t = THREE_HALF*R_gas*tmpT

    S_t = R_gas*(log(q_t)+FIVE_HALF)

    !*** Output of translational contributions *******************************************
    write(output_unit,05) total_mass/u_Atom_Mass
    05 format('  |     Translational Contributions',10X,'total mass: ',F9.3,'  au',20X'|')
    write(output_unit,02)

    write(output_unit,06) E_t(1)/THOUSAND, S_t(1)/THOUSAND
    06 format('  |',5X,'E_t  =  ',E12.6,'  kJ/mol',9X,'S_t  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_t(1)/E_h2Jmol, S_t(1)/E_h2Jmol
    07 format('  |',12X,'(',E12.6,'      au)',15X,'(',E12.6,'      au/K)',14X,'|')
    write(output_unit,02)
    !*** Loop over rest of temperature points ********************************************
    do K = 2,NTpoints
      write(output_unit,0601) E_t(K)/THOUSAND, S_t(K)/THOUSAND
      0601 format('  |',13X,E12.6,'  kJ/mol',17X,E12.6,'  kJ/K*mol',15X,'|')
      write(output_unit,07) E_t(K)/E_h2Jmol, S_t(K)/E_h2Jmol
      write(output_unit,02)
    enddo
    write(output_unit,01)
    write(output_unit,02)

    !*************************************************************************************
    !    calculation of the electronic contributions are not carried out so far          !
    !    the relevant expressions would be:                                              !
    !                                                                                    !
    !         q_e = g                                                                    !
    !                                                                                    !
    !         E_e := 0                                                                   !
    !                                                                                    !
    !         S_e = R*ln(q_e)                                                            !
    !                                                                                    !
    !*************************************************************************************
    !q_e = real(g,r8_kind)

    E_e = ZERO

    S_e = R_gas*log(ONE)

    !*** Output of electronic contributions **********************************************
    write(output_unit,08)
    08 format('  |     Electronic Contributions',13X,'(correct manually if necessary)',14X,'|')
    write(output_unit,02)
    !write(output_unit,09) E_e/THOUSAND, S_e/THOUSAND
    !09 format('  |',5X,'E_e  =  ',E12.6,'  kJ/mol',9X,'S_e  =  ',E12.6,'  kJ/K*mol',15X,'|')
    !write(output_unit,07) E_e/E_h2Jmol, S_e/E_h2Jmol
    !write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*************************************************************************************
    !    calculation of vibrational contributions using the expressions                  !
    !                                                                                    !
    !         T_v(K) = h*v_K/k_B                     (vibrational temperature)           !
    !                                                                                    !
    !         E_v = R*\sum_{K}T_v(K)*(1/2 + 1/(exp(T_v(K)/T)-1))                         !
    !                                                                                    !
    !         S_v = R*\sum_{K}((T_v(K)/T)/(exp(T_v(K)/T)-1) + ln(1-1/(exp(T_v(K)/T))))   !
    !                                                                                    !
    !*************************************************************************************
    T_v = freq*h_Planck/k_B
    ZPE_v = R_gas*sum(T_v)/TWO

    !*** Loop over temperature points ****************************************************
    do J = 1,NTpoints
      DO K = 1,size(freq)
        E_v(J) = E_v(J) + T_v(K)*(HALF + ONE/(exp(T_v(K)/tmpT(J))-ONE))
        S_v(J) = S_v(J) + ((T_v(K)/tmpT(J))/(exp(T_v(K)/tmpT(J))-ONE)) - log(ONE - exp(-T_v(K)/tmpT(J)))
      END DO
    enddo
    E_v = E_v*R_gas
    S_v = S_v*R_gas

    !*** Output of vibrational contributions *********************************************
    write(output_unit,10) num_pos_freq
    10 format('  |     Vibrational Contributions',10X,I4,' positive frequencies',22X,'|')
    write(output_unit,11) size(fre) - num_pos_freq
    11 format('  |',40X,I4,' negative frequencies (omitted)',12X,'|')
    write(output_unit,02)
    write(output_unit,12) ZPE_v/THOUSAND
    12 format('  |',5X,'ZPE  =  ',E12.6,'  kJ/mol',54X,'|')
    write(output_unit,13) ZPE_v/E_h2Jmol
    13 format('  |',12X,'(',E12.6,'      au)',53X'|')
    write(output_unit,02)
    write(output_unit,14) E_v(1)/THOUSAND, S_v(1)/THOUSAND
    14 format('  |',5X,'E_v  =  ',E12.6,'  kJ/mol',9X,'S_v  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_v(1)/E_h2Jmol, S_v(1)/E_h2Jmol
    write(output_unit,02)
    !*** Loop over temperature points ****************************************************
    do J = 2,NTpoints
      write(output_unit,0601) E_v(J)/THOUSAND, S_v(J)/THOUSAND
      write(output_unit,07) E_v(J)/E_h2Jmol, S_v(J)/E_h2Jmol
      write(output_unit,02)
    enddo
    write(output_unit,01)
    write(output_unit,02)

    !*************************************************************************************
    !    calculation of rotational contributions using the expressions                   !
    !                                                                                    !
    !         T_r(x,y,z) = h^2/(8*pi^2*k_B*I(x,y,z))                                     !
    !                                                                                    !
    !         q_r = (pi^0.5/sigma)*(T^1.5/(T_r(x)*T_r(y)*T_r(z))^0.5)                    !
    !                                                                                    !
    !         E_r = 3/2*R*T                                                              !
    !                                                                                    !
    !         S_r = R*(ln(q_r)+3/2)                                                      !
    !                                                                                    !
    !*************************************************************************************

    !*** determine linearity of molecule and check if Mom_Inert positive *****************
    lin_mol = .FALSE.
    IF (minval(mom_I) < liC) THEN
      lin_mol = .TRUE.
      rot_eps_shape = 'disk'
    ELSE IF (minval(mom_I) < ZERO) THEN
      write(06,*) "ERROR: negative moments of inertia!"
    END IF

    !*** determine shape of molecule *****************************************************
    !*** Spherical molecule? *************************************************************
    IF (mom_I(1)-mom_I(2) < liC .AND. mom_I(2)-mom_I(3) < liC) THEN
      rot_eps_shape = 'spherical'
    !*** Two axis equal ? ****************************************************************
    ELSE IF ((mom_I(1)-mom_I(2)) < liC) THEN
        rot_eps_shape = 'oblate'
    ELSE IF ((mom_I(2)-mom_I(3)) < liC) THEN
        rot_eps_shape = 'prolate'
    !*** Triaxial? ***********************************************************************
    ELSE
      rot_eps_shape = 'triaxial'
    END IF

    !*** Calculate rotational temperatures ***********************************************
    DO K = 1,3
      T_r(K) = hbar_Planck*hbar_Planck/(TWO*k_B*mom_I(K))
    END DO

    !*** calculate rotational contributions **********************************************
    IF (lin_mol) THEN
      q_r = tmpT/(sigma*T_r(1))
      E_r = R_gas*tmpT
      S_r = R_gas*(log(q_r)+ONE)
    ELSE
      q_r = sqrt(PI)*(tmpT**THREE_HALF)/(real(sigma,r8_kind)*sqrt(product(T_r)))
      E_r = THREE_HALF*R_gas*tmpT
      S_r = R_gas*(log(q_r)+THREE_HALF)
    END IF

    !*** Output of rotational contributions and shape of rotational ellipsoid ************
    write(output_unit,15) rot_eps_shape
    15 format('  |     Rotational Contributions',13X,A9,' rotational ellipsoid',15X,'|')
    write(output_unit,16) sigma
    16 format('  |',42X,'symmetry index : ',I3,25X,'|')
    write(output_unit,02)

    write(output_unit,17) E_r(1)/THOUSAND, S_r(1)/THOUSAND
    17 format('  |',5X,'E_r  =  ',E12.6,'  kJ/mol',9X,'S_r  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_r(1)/E_h2Jmol, S_r(1)/E_h2Jmol
    write(output_unit,02)
    !*** Loop over temperature points ****************************************************
    do K = 2,NTpoints
      write(output_unit,0601) E_r(K)/THOUSAND, S_r(K)/THOUSAND
      write(output_unit,07) E_r(K)/E_h2Jmol, S_r(K)/E_h2Jmol
      write(output_unit,02)
    enddo
    write(output_unit,01)
    write(output_unit,02)

    !*** Sum calculated contributions to obtain Thermodynamic_Properties *****************
    U = E_t + E_r
    U = U + E_v
    S = S_v + S_r
    S = S + S_t
    entH = U + R_gas*tmpT
    Gfe = entH - tmpT*S

    !*** Output of total thermodynamic quantities ****************************************
    write(output_unit,18)
    18 format('  |     Summary: THERMODYNAMICAL PROPERTIES',47X,'|')
    write(output_unit,02)
    write(output_unit,19) U(1)/THOUSAND, S(1)/THOUSAND
    19 format('  |',5X,'U    =  ',E12.6,'  kJ/mol',9X,'S    =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) U(1)/E_h2Jmol, S(1)/E_h2Jmol
    write(output_unit,02)
    !*** Loop over temperature points ****************************************************
    do K = 2,NTpoints
      write(output_unit,0601) U(K)/THOUSAND, S(K)/THOUSAND
      write(output_unit,07)   U(K)/E_h2Jmol, S(K)/E_h2Jmol
      write(output_unit,02)
    enddo

    write(output_unit,02)

    write(output_unit,20) entH(1)/THOUSAND, Gfe(1)/THOUSAND
    20 format('  |',5X,'H    =  ',E12.6,'  kJ/mol',9X,'G    =  ',E12.6,'  kJ/mol',17X,'|')
    write(output_unit,07) entH(1)/E_h2Jmol, Gfe(1)/E_h2Jmol
    write(output_unit,02)
    !*** Loop over temperature points ****************************************************
    do K = 2,NTpoints
      write(output_unit,0601) entH(K)/THOUSAND, Gfe(K)/THOUSAND
      write(output_unit,07)   entH(K)/E_h2Jmol, Gfe(K)/E_h2Jmol
      write(output_unit,02)
    enddo
    write(output_unit,01)
    write(output_unit,02)

    !*** Adding electronic contributions *************************************************
    Gfe_tot = Gfe + tot_energy

    !*** Output of corrected e_sum value *************************************************
    write(output_unit,21)
    21 format('  |     End of thermodynamic properties section',43X,'|')
    write(output_unit,02)
    write(output_unit,22)  Gfe_tot(1)/THOUSAND
    22 format('  |     G_tot = ',E17.10,' kJ/mol',13X,'(including "e_sum")',18X,'|')
    write(output_unit,23)  Gfe_tot(1)/E_h2Jmol
    23 format('  |',12X,'(',E17.10,'     au)',49X,'|')
    write(output_unit,02)
    2201 format('  |',13X,E17.10,' kJ/mol',50X,'|')
    do K = 2,NTpoints
      write(output_unit,2201) Gfe_tot(K)/THOUSAND
      write(output_unit,23)   Gfe_tot(K)/E_h2Jmol
      write(output_unit,02)
    enddo
    write(output_unit,01)
    write(output_unit,*)

    deallocate(freq,T_v)
    return
  end subroutine thermodynamic_properties
  !***************************************************************************************

#if 0
  subroutine thermodynamic_properties_print(output_unit,T,P,tot_mass,E_t,S_t,E_e,S_e&
                                           ,num_freq,num_pos_freq,E_v,S_v,ZPE_v&
                                           ,rot_eps_shape,sigma,E_r,S_r,U,S,H,Gfe,Gfe_tot)
  !---------------------------------------------------------------------------------------
  !                                                                                      !
  !  Purpose: Print results of thermodynamic properties module into "output_unit"        !
  !                                                                                      !
  !                                                                                      !
  !  Module called by: "thermodynamic_properties"  in this file                          !
  !                                                                                      !
  !                                                                                      !
  !  References: none                                                                    !
  !                                                                                      !
  !                                                                                      !
  !  Author: TMS                                                                         !
  !  Date: March 2009                                                                    !
  !                                                                                      !
  !                                                                                      !
  !---------------------------------------------------------------------------------------
    use constants , only             : THOUSAND, E_h2Jmol,u_Atom_Mass
    !** End of interface *****************************************************************
    integer(kind=i4_kind),intent(in)        :: output_unit        ! Output file designator

    real(kind=r8_kind),intent(in)           :: T &                       ! temperature [K]
                                             , P                           ! pressure [Pa]

    real(kind=r8_kind),intent(in)           :: tot_mass, E_t, S_t

    real(kind=r8_kind),intent(in)           :: E_e, S_e

    integer(kind=i4_kind),intent(in)        :: num_freq, num_pos_freq

    real(kind=r8_kind),intent(in)           :: E_v,S_v,ZPE_v

    integer(kind=i4_kind),intent(in)        :: sigma

    real(kind=r8_kind),intent(in)           :: E_r,S_r

    character(9),intent(in)                 :: rot_eps_shape

    real(kind=r8_kind),intent(in)           :: U,S,H,Gfe,Gfe_tot

    !*** Write initial output ************************************************************
    write(output_unit,*)
    write(output_unit,01)
    01 format('  +---------------------------------------------------------------------------------------+')
    write(output_unit,02)
    02 format('  |',87X,'|')
    write(output_unit,03)
    03 format('  |     Thermodynamic properties section:',49X,'|')
    write(output_unit,02)
    write(output_unit,04) T, P
    04 format('  |',5X,'T   = ',F8.3,' K',21X,'p = ',F9.1,' Pa',29X,'|')
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Output of translational contributions *******************************************
    write(output_unit,05) tot_mass/u_Atom_Mass
    05 format('  |     Translational Contributions',10X,F4.1,20X'|')
    write(output_unit,02)
    write(output_unit,06) E_t/THOUSAND, S_t/THOUSAND
    06 format('  |',5X,'E_t  =  ',E12.6,'  kJ/mol',9X,'S_t  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_t/E_h2Jmol, S_t/E_h2Jmol
    07 format('  |',12X,'(',E12.6,'      au)',15X,'(',E12.6,'      au/K)',14X,'|')
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Output of electronic contributions **********************************************
    write(output_unit,08)
    08 format('  |     Electronic Contributions',13X,'(correct manually if necessary)',14X,'|')
    write(output_unit,02)
    write(output_unit,09) E_e/THOUSAND, S_e/THOUSAND
    09 format('  |',5X,'E_e  =  ',E12.6,'  kJ/mol',9X,'S_e  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_e/E_h2Jmol, S_e/E_h2Jmol
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Output of vibrational contributions *********************************************
    write(output_unit,10) num_pos_freq
    10 format('  |     Vibrational Contributions',10X,I4,' positive frequencies',22X,'|')
    write(output_unit,11) num_freq - num_pos_freq
    11 format('  |'40X,I4,' negative frequencies (omitted)',12X,'|')
    write(output_unit,02)
    write(output_unit,12) ZPE_v/THOUSAND
    12 format('  |',5X,'ZPE  =  ',E12.6,'  kJ/mol',54X,'|')
    write(output_unit,13) ZPE_v/E_h2Jmol
    13 format('  |',12X,'(',E12.6,'      au)',53X'|')
    write(output_unit,02)
    write(output_unit,14) E_v/THOUSAND, S_v/THOUSAND
    14 format('  |',5X,'E_v  =  ',E12.6,'  kJ/mol',9X,'S_v  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_v/E_h2Jmol, S_v/E_h2Jmol
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Output of shape contributions ***************************************************
    write(output_unit,15) rot_eps_shape
    15 format('  |     Rotational Contributions',13X,A9,' rotational ellipsoid',15X,'|')
    write(output_unit,16) sigma
    16 format('  |',42X,'symmetry index : ',I3,25X,'|')
    write(output_unit,02)
    write(output_unit,17) E_r/THOUSAND, S_r/THOUSAND
    17 format('  |',5X,'E_r  =  ',E12.6,'  kJ/mol',9X,'S_r  =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) E_r/E_h2Jmol, S_r/E_h2Jmol
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Output of total thermodynamic quantities ****************************************
    write(output_unit,18)
    18 format('  |     Summary: THERMODYNAMICAL PROPERTIES',47X,'|')
    write(output_unit,02)
    write(output_unit,19) U/THOUSAND, S/THOUSAND
    19 format('  |',5X,'U    =  ',E12.6,'  kJ/mol',9X,'S    =  ',E12.6,'  kJ/K*mol',15X,'|')
    write(output_unit,07) U/E_h2Jmol, S/E_h2Jmol
    write(output_unit,20) H/THOUSAND, Gfe/THOUSAND
    20 format('  |',5X,'H    =  ',E12.6,'  kJ/mol',9X,'G    =  ',E12.6,'  kJ/mol',17X,'|')
    write(output_unit,07) H/E_h2Jmol, Gfe/E_h2Jmol
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,02)

    !*** Adding electronic contributions *************************************************
    write(output_unit,21)
    21 format('  |     End of thermodynamic properties section',43X,'|')
    write(output_unit,02)
    write(output_unit,22)  Gfe_tot/THOUSAND
    22 format('  |     G_tot = ',E17.10,' kJ/mol',13X,'(including "e_sum")',18X,'|')
    write(output_unit,23)  Gfe_tot/E_h2Jmol
    23 format('  |',12X,'(',E17.10,'     au)',49X,'|')
    write(output_unit,02)
    write(output_unit,01)
    write(output_unit,*)
  end subroutine thermodynamic_properties_print
  !***************************************************************************************
#endif

  function sort3(tempvec) result(sortvec)
    implicit none
    real(r8_kind), intent(in) :: tempvec(3)
    real(r8_kind)             :: sortvec(3)
    !** End of interface *****************************************************************
    sortvec(1) = maxval(tempvec)
    sortvec(3) = minval(tempvec)
    IF (tempvec(1) < sortvec(1) .AND. tempvec(1) > sortvec(3)) THEN
      sortvec(2) = tempvec(1)
    ELSE IF (tempvec(2) < sortvec(1) .AND. tempvec(2) > sortvec(3)) THEN
      sortvec(2) = tempvec(2)
    ELSE
      sortvec(2) = tempvec(3)
    END IF
  END function sort3

  !--------------- End of module ---------------------------------------------------------
end module thermodyn_prop_module

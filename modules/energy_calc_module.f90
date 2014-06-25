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
module energy_calc_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: Contains all routines necessary to calculate
  !           the two- and three-center energies
  !           - kinetic energy           E_KIN
  !           - nuclear attr. energy     E_NUC
  !           - coulomb energy           E_COUL
  !           - exch.-corr. energy       E_XC
  !           - e-e part of the coulomb
  !             energy                   E_COUL_2Z
  !           - core-core interaction
  !             energy                   E_NUC_CORE
  !
  !  Contents:
  !     Variables:
  !          (i) PRIVATE
  !                      e_kin,e_nuc,e_coul,e_xc,e_coul_2z,e_nuc_core
  !
  !     Routines:
  !          (i) PUBLIC
  !                      energ_calc_2z  -> calculates e_kin,e_nuc
  !                      energ_coul_2z  -> calcuates 2z-part of e_coul
  !                      init_energy    -> sets all parts of energy to zero
  !                      get_energy     -> get total energy
  !                      read_core_core -> as long as there is no integral
  !                                        part, we have to read this in
  !                                        from the old LCGTO
  !                      write_energies -> prints all energies
  !                      write_energy   -> prints a single energy value
  !
  !
  !  Module called by: build_hamiltonian
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   5/97
  ! Description: added point charges to core_interaction_calc
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   7/7/97
  ! Description: - packed storage mode for mat_charge introduced
  !                ij=0; DO i=1,n_ch; DO j=1,i; ij=ij+1; ij=>(i,j); ENDDO; ENDDO
  !              - function add_energy_contributions introduced
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   10/97
  ! Description: extension to spin orbit
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   12/97
  ! Description: output of electrostatic contribution to free solvation
  !              energy was added
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined datatypes
  use symmetry_data_module, only :sym
  use fit_coeff_module, only: fit

  implicit none
  private
  save
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public :: energ_coul_2z                                                      &
          , init_energy                                                        &
          , write_energies                                                     &
          , write_energy                                                       &
          , core_interaction_calc                                              &
          , direct_2c_energy_calc_and_add                                      &
          , direct_mat_energy_calc_and_add

  public :: get_energy!(key=val)
  public :: set_energy!(key=val)

  public :: energy_calc_reduce!()

#ifdef WITH_EPE
  public :: addional_core_ewpc_interaction
#endif

  public :: sym,fit,arrmat3,arrmat2

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----

  !
  ! Components of the energy, either computed here
  ! or set with set_energy(name=value)
  !
  real(r8_kind) :: e_kin = 0.0
  real(r8_kind) :: e_nuc = 0.0
  real(r8_kind) :: e_coul = 0.0
  real(r8_kind) :: e_exex = 0.0
  real(r8_kind) :: e_xc = 0.0
  real(r8_kind) :: e_dft_plus_u = 0.0 ! DFT+U energy correction, computed in dft_plus_u_module
  real(r8_kind) :: e_efield = 0.0
  real(r8_kind) :: e_coul_2z = 0.0  ! 2-center part of E_COUL
  real(r8_kind) :: e_nuc_core = 0.0 ! core-core interaction energy
  real(r8_kind) :: e_rism = 0.0 ! RISM solvation energy, not PCM!
  real(r8_kind) :: e_tot = 0.0

  !
  ! Added-up energies, these are derived quantities, for legacy reasons
  ! computed in update_energies() but returned by get_energy()
  !
  real(kind=r8_kind)     :: a_kin,a_nuc,a_coul,a_xc,a_tot,a_energy_solv_el,a_energy_solv_nonel
  real(kind=r8_kind)     :: a_efield = 0.0_r8_kind,e_nuc_ewpc_add = 0.0_r8_kind

  ! array of energies for the different post_scf functionals
  real(kind=r8_kind), allocatable :: a_xc_arr(:), a_tot_arr(:)
  ! array of names for the different post_scf functionals
  character(len=32), allocatable :: name_xc_arr(:)
  ! (pointers because they will be alocated by external xc_func_get_exc())
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine energ_coul_2z(coeff_charge, e)
    !  Purpose: calculate the 2-center part of the coulomb energy
    !           e_coul_2z = -1/2 sum_{k,l} (ak [Fk|Fl] al).
    !           Please note the sign!
    !           The necessary two-center integrals MAT_CHARGE
    !           are provided by the module mat_charge_module
    !  Subroutine called by: build_hamiltonian
    !
    use mat_charge_module, only: mat_charge
    implicit none
    real(r8_kind), intent(in)  :: coeff_charge(:)
    real(r8_kind), intent(out) :: e
    ! fit coefficients for charge density
    !** End of interface *****************************************

    real(kind=r8_kind), parameter :: half = 0.5_r8_kind
    integer(i4_kind) :: n, i, j, ij

    n = size(coeff_charge)
    ASSERT(n*(n+1)/2==size(mat_charge))

    ! now calculate e_coul_2z
    ! using a very elegant algorithm to do
    ! sum_[k,l](a_k [f_k|f_l] a_l)

    e = 0.0

    ij = 0
    do i = 1, n
       ! off-diagonal contributions:
       do j = 1, i-1
          ij = ij + 1

          e = e - coeff_charge(i) * coeff_charge(j) * mat_charge(ij) ! (j,i)

       enddo
       ! diagonal contribution:
       ij = ij + 1
       e = e - half * coeff_charge(i) * coeff_charge(i) * mat_charge(ij) ! (i,i)
    enddo

  end subroutine energ_coul_2z

  !*************************************************************

  subroutine update_energies(post_scf)
    ! Purpose: update the interesting energies
    !          E_kin, E_nuc, E_coul, E_xc and E_tot
    use xc_hamiltonian, only: xc_get_exc
    use xcmda_hamiltonian, only: xcmda_get_exc
    use options_module, only: options_xcmode, &
                              xcmode_numeric_exch, xcmode_exchange_fit, &
                              xcmode_model_density, xcmode_extended_mda
    use operations_module, only: operations_solvation_effect
    use xc_func, only: xc_func_get_exc
    use solv_cavity_module, only: disp_rep_energy, cavitation_energy, energy_cav
    use solv_electrostat_module,  only: energy_solv_el
    use disp_rep_module, only: E_disp, E_rep
    use induced_dipoles_module, only: en_ind_dipole
#ifdef WITH_EFP
    use efp_efp_module, only: efp_efp_en
#endif
    implicit none
    logical, intent(in), optional :: post_scf
    !** End of interface *****************************************

    integer(kind=i4_kind) :: alloc_stat

    !
    ! FIXME: make sure the additions are zero if contributions are off
    ! instead of introducing many  nested ifs. For example, energy_cav
    ! is   only   set   in   energy_and_grad_of_cavity   called   from
    ! main_gradient  only if  solvation is  on.  Otherwise  it remains
    ! uninitialized. It  would be nice  if each term  just contributed
    ! zero, when turned off.
    !
    a_energy_solv_nonel = 0.0
    a_energy_solv_el = 0.0
    if (operations_solvation_effect) then
       if (cavitation_energy) then
          a_energy_solv_nonel = energy_cav
       endif

       if (disp_rep_energy) then
          a_energy_solv_nonel = a_energy_solv_nonel + E_disp + E_rep
       endif

       a_energy_solv_el = energy_solv_el
    endif

    a_kin    = e_kin
    a_nuc    = e_nuc
    a_coul   = e_coul + e_coul_2z
    a_efield = e_efield
    !
    ! To make it compatible with the legacy expression
    ! e_tot = sum_ij { P_ij ( T_ij + Z_ij + J_ij )}
    e_tot  = e_kin                                                             &
           + e_nuc                                                             &
           + e_coul

    a_tot  = e_tot                                                             &
           + e_exex                                                            &
           + e_coul_2z + e_nuc_core                                            &
           + a_efield &
           + a_energy_solv_el + a_energy_solv_nonel &
           + e_nuc_ewpc_add

    a_tot = a_tot + en_ind_dipole
#ifdef WITH_EFP
    a_tot = a_tot +efp_efp_en
#endif

    !
    ! This contribution is nonzero only if DFT+U is in effect
    ! (e_dft_plus_u is set by set_energy() from ham_calc_module)
    !
    a_tot = a_tot + e_dft_plus_u

    ! now load SCF XC-Energies
    select case (options_xcmode())
    case (xcmode_model_density,xcmode_extended_mda)
       a_xc = xcmda_get_exc()
    case (xcmode_numeric_exch)
       a_xc = xc_get_exc()
    case (xcmode_exchange_fit)
       a_xc = e_xc
    end select

    !
    ! The term should be zero in a regular case or #ifndef WITH_BGY3D.
    ! See main_gradient() where this addition is computed and saved in
    ! the module private variable:
    !
    a_tot = a_tot + e_rism

    !
    ! FIXME: a singularity here! In SCF the e_xc energy is returned by
    ! xc_get_exc(), see above.  However after SCF it is  assumed to be
    ! the   last   entry  in   the   array   of   many  different   XC
    ! funcitonals.  This  may  be  an atavism  of  so-called  post-scf
    ! calculations  where a  functional employed  in post-scf  was not
    ! necessarily  the same  as  in SCF.  Say  a GGA  post-scf on  LDA
    ! density.
    !
    if (present(post_scf)) then
       if (post_scf) then
          call xc_func_get_exc(a_xc_arr,name_xc_arr)
          allocate(a_tot_arr(size(a_xc_arr)),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               "update_energies: allocation failed")
          a_tot_arr(:) = a_tot + a_xc_arr(:)
          a_xc=a_xc_arr(size(a_xc_arr))
       endif
    endif
    a_tot = a_tot + a_xc

  end subroutine update_energies

  !*************************************************************

  subroutine get_energy(tot)
    !
    ! Returns the total energy making sure all terms have been added.
    !
    implicit none
    real(r8_kind), intent(out) :: tot
    !** End of interface *****************************************

    ! This sets module private vars a_*, including a_tot,
    ! Please make sure this sub is idempotent (calling more than
    ! one time should have no net effect):
    call update_energies()

    tot = a_tot
  end subroutine get_energy

  !*************************************************************

  subroutine core_interaction_calc()
    !  Purpose: This routine actually calculates the
    !           electrostatic core-core interaction.
    !           The energy is stored in the (private)
    !           variable e_nuc_core
    !           The electrostatic core-pointcharge interaction
    !           is also included in this value.
    !           Note that the pointcharge-pointcharge interaction
    !           is omitted.
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use unique_atom_module, only : unique_atoms, N_unique_atoms, &
                                   pseudopot_present
    use pointcharge_module, only : pointcharge_array, pointcharge_N, &
         n_timps, unique_timps
    use point_dqo_module, only : calc_nuc_dqo_energy
#ifdef WITH_EFP
    use pointcharge_module, only : print_energy
    use qmmm_interface_module, only : efp
#endif
#ifdef WITH_EPE
    use ewaldpc_module
#endif
    use iounitadmin_module
    implicit none
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)          :: z1,z2,dist,dist2
#ifdef WITH_EPE
    real(kind=r8_kind)          :: e_nuc_ewpc
#endif
    real(kind=r8_kind)          :: zc1,zc2
    integer(kind=i4_kind)       :: na, nb, eq_a, eq_b, n_equal_charges
    integer(kind=i4_kind)       :: N_equal_atoms_a, N_equal_atoms_b
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:)
    real(kind=r8_kind)          :: e_nuc_pc_efp,C,A
    !------------ Executable code --------------------------------

    e_nuc_core = 0.0_r8_kind
    e_nuc_ewpc_add = 0.0_r8_kind !!!!!!!!!!!! AS

    ! atom - atom interaction
    unique1: do na=1,N_unique_atoms+n_timps
       if(na.gt.N_unique_atoms) then
          if(unique_timps(na-N_unique_atoms)%moving_atom.eq.0) cycle unique1
          z1 = unique_timps(na-N_unique_atoms)%Z
          zc1= unique_timps(na-N_unique_atoms)%ZC
          xa => unique_timps(na-N_unique_atoms)%position
          N_equal_atoms_a=unique_timps(na-N_unique_atoms)%N_equal_atoms
       else
          z1 = unique_atoms(na)%Z
          zc1= unique_atoms(na)%ZC
          xa => unique_atoms(na)%position
          N_equal_atoms_a=unique_atoms(na)%N_equal_atoms
       end if

       unique2: do nb=1,N_unique_atoms+n_timps
          if(nb.gt.N_unique_atoms) then
             if(unique_timps(nb-N_unique_atoms)%moving_atom.eq.0) cycle unique2
             z2 = unique_timps(nb-N_unique_atoms)%Z
             zc2= unique_timps(nb-N_unique_atoms)%ZC
             xb => unique_timps(nb-N_unique_atoms)%position
             N_equal_atoms_b=unique_timps(nb-N_unique_atoms)%N_equal_atoms
          else
             z2 = unique_atoms(nb)%Z
             zc2= unique_atoms(nb)%ZC
             xb => unique_atoms(nb)%position
             N_equal_atoms_b=unique_atoms(nb)%N_equal_atoms
          end if


          equal1: do eq_a=1,N_equal_atoms_a
             equal2: do eq_b=1,N_equal_atoms_b

                if((na.eq.nb).and.(eq_a.eq.eq_b)) then
                   cycle equal2
                endif

                dist = sqrt(sum((xa(:,eq_a)-xb(:,eq_b))**2))
                if (pseudopot_present) then
                   e_nuc_core = e_nuc_core + (z1-zc1)*(z2-zc2)/dist
                else
                   e_nuc_core = e_nuc_core + z1*z2/dist
                end if
             enddo equal2
          enddo equal1

       enddo unique2
    enddo unique1
    e_nuc_core=0.5_r8_kind*e_nuc_core
!!$        print*,'e_nuc_core',e_nuc_core
!!$        write(output_unit,*)'e_nuc_core', e_nuc_core
    ! add atom - pointcharge interaction

    e_nuc_pc_efp=0.0_r8_kind
    unique_pc: do nb=1,pointcharge_N + n_timps
       if(nb<=n_timps) then
          if(unique_timps(nb)%moving_atom.ne.0) cycle unique_pc
          z2 =  unique_timps(nb)%Z - unique_timps(nb)%ZC
          xb => unique_timps(nb)%position
          n_equal_charges = unique_timps(nb)%N_equal_atoms
       else
          z2 = pointcharge_array(nb-n_timps)%Z
          xb => pointcharge_array(nb-n_timps)%position
          n_equal_charges = pointcharge_array(nb-n_timps)%N_equal_charges
          C=pointcharge_array(nb-n_timps)%C
          A=pointcharge_array(nb-n_timps)%A
       end if
       unique_a: do na=1,N_unique_atoms+n_timps
          if(na.gt.N_unique_atoms) then
             if(unique_timps(na-N_unique_atoms)%moving_atom.eq.0) cycle unique_a
             z1 = unique_timps(na-N_unique_atoms)%Z - &
                  unique_timps(na-N_unique_atoms)%ZC
             xa => unique_timps(na-N_unique_atoms)%position
             N_equal_atoms_a=unique_timps(na-N_unique_atoms)%N_equal_atoms
          else
             z1 = unique_atoms(na)%Z-unique_atoms(na)%ZC
             xa => unique_atoms(na)%position
             N_equal_atoms_a=unique_atoms(na)%N_equal_atoms
          end if


          equal_a: do eq_a=1,N_equal_atoms_a
             equal_pc: do eq_b=1,n_equal_charges
                dist2 = sum((xa(:,eq_a)-xb(:,eq_b))**2)
                dist = sqrt(dist2)
                if(dist <= 1.0e-10_r8_kind) cycle equal_pc
                e_nuc_core = e_nuc_core + z1*z2/dist
#ifdef WITH_EFP
                e_nuc_pc_efp = e_nuc_pc_efp + z1*z2/dist
                if(C /= 0.0_r8_kind) then
                   e_nuc_core = e_nuc_core - C*exp(-A*dist2)*z1*z2/dist
                   e_nuc_pc_efp = e_nuc_pc_efp - C*exp(-A*dist2)*z1*z2/dist
                end if
#endif
             enddo equal_pc
          enddo equal_a

       enddo unique_a
    enddo unique_pc
#ifdef WITH_EFP
    if((pointcharge_N > 0) .and. efp .and. print_energy) &
         print*,'NUC-PC  EFP ',e_nuc_pc_efp
#endif

#ifdef WITH_EPE
    ! add atom - ewald file  pointcharge interaction
    if(EWPC_n.gt.0) then
       e_nuc_ewpc = 0.0_r8_kind
       unique_ewpc: do nb=1,EWPC_N
          z2 = ewpc_array(nb)%Z
          xb => ewpc_array(nb)%position
          unique_at: do na=1,N_unique_atoms+n_timps
             if(na.gt.N_unique_atoms) then
                if(unique_timps(na-N_unique_atoms)%moving_atom.eq.0) cycle unique_at
                z1 = unique_timps(na-N_unique_atoms)%Z- &
                     unique_timps(na-N_unique_atoms)%ZC
                xa => unique_timps(na-N_unique_atoms)%position
                N_equal_atoms_a=unique_timps(na-N_unique_atoms)%N_equal_atoms
             else
                z1 = unique_atoms(na)%Z - unique_atoms(na)%ZC
                xa => unique_atoms(na)%position
                N_equal_atoms_a=unique_atoms(na)%N_equal_atoms
             end if

             equal_at: do eq_a=1,N_equal_atoms_a
                equal_ewpc: do eq_b=1,ewpc_array(nb)%N_equal_charges

                   dist = sqrt(sum((xa(:,eq_a)-xb(:,eq_b))**2))
                   e_nuc_ewpc = e_nuc_ewpc + z1*z2/dist

                enddo equal_ewpc
             enddo equal_at

          enddo unique_at
       enddo unique_ewpc
!!$    print* ,'energy of interaction of atomic nuclei with ewald file PC ', e_nuc_ewpc
       e_nuc_core=e_nuc_core+e_nuc_ewpc
    endif !ewpc_n
#endif

    e_nuc_core=e_nuc_core+calc_nuc_dqo_energy()

  end subroutine core_interaction_calc

#ifdef WITH_EPE
  !*************************************************************
  subroutine addional_core_ewpc_interaction(m_charge)

    use unique_atom_module, only: unique_atoms, N_unique_atoms
    use ewaldpc_module
    implicit none
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind),intent(in) :: m_charge(N_unique_atoms)

    real(kind=r8_kind)          :: z1,z2,dist
    integer(kind=i4_kind)       :: na, nb, eq_a, eq_b
    integer(kind=i4_kind)       :: N_equal_atoms_a
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:)
    !------------ Executable code --------------------------------

       e_nuc_ewpc_add = 0.0_r8_kind
    if(EWPC_n.gt.0) then
       unique_ewpc: do nb=1,EWPC_N
          z2 = ewpc_array(nb)%Z
          xb => ewpc_array(nb)%position
          unique_at: do na=1,N_unique_atoms
             z1 =m_charge(na)
             xa => unique_atoms(na)%position
             N_equal_atoms_a=unique_atoms(na)%N_equal_atoms

             equal_at: do eq_a=1,N_equal_atoms_a
                equal_ewpc: do eq_b=1,ewpc_array(nb)%N_equal_charges

                   dist = sqrt(sum((xa(:,eq_a)-xb(:,eq_b))**2))
                   e_nuc_ewpc_add = e_nuc_ewpc_add + z1*z2/dist

             enddo equal_ewpc
          enddo equal_at

          enddo unique_at
       enddo unique_ewpc
    endif !ewpc_n
  end subroutine addional_core_ewpc_interaction
#endif

  !*************************************************************
  subroutine write_energies (unit, post_scf)
    !
    ! Write the private variables E_KIN, E_NUC, E_COUL, E_XC to unit.
    !
    ! FIXME: not only printing, but also a side effect of deallocating
    ! name_xc_arr, a_xc_arr, and a_tot_arr in the post-scf call.
    !
    use constants, only: kcal ! == 1 / 627.49
    use solv_cavity_module, only: energy_cav, cavitation_energy, disp_rep_energy
    ! use solv_electrostat_module, only: energy_solv_el
    use disp_rep_module, only: E_disp, E_rep
    use operations_module, only: operations_solvation_effect
    implicit none
    integer (i4_kind), intent(in) :: unit
    logical, intent (in), optional :: post_scf
    !** End of interface *****************************************

    logical :: post_scf_energy
    integer (i4_kind) :: i, alloc_stat
    real (r8_kind) :: e_sum

    !
    ! FIXME: Using save'ed variables to print the energy differences
    ! between SCF iterations:
    !
    real (r8_kind) :: pre_xc  = 0.0_r8_kind ! exchange-correlation
    real (r8_kind) :: pre_el  = 0.0_r8_kind ! electronic
    real (r8_kind) :: pre_tot = 0.0_r8_kind ! total
    real (r8_kind) :: pre_sel = 0.0_r8_kind ! solvation-electronic
    real (r8_kind) :: pre_dft_plus_u = 0.0_r8_kind ! DFT+U

    ! Status of optional arguments propagates, think of passing
    ! null/non-null pointers:
    call update_energies(post_scf)

    if (present(post_scf)) then
       post_scf_energy = post_scf
    else
       post_scf_energy = .false.
    endif

    1000 format(2X,A6,' =  ',F25.12:'  (',A,')')
    1001 format(2X,A6,' =  ',F25.12:'  [',F20.12,']')

    !
    ! Here "a_tot"  is a  module global variable  that was  updated in
    ! update_energies()  to  contain   total  energy.   However,  some
    ! contributions (such  as DFT-D) are handled here  (FIXME!).  I am
    ! not  sure  if  DFT-D  energy  has  a  meaningfull  value  during
    ! SCF. "e_sum" is just a  local variable. Read-only "a_tot" is not
    ! referenced below this line:
    !
    e_sum = a_tot

    if (post_scf_energy) then
       ASSERT(allocated(a_xc_arr))
       ASSERT(allocated(a_tot_arr))

       if (unit > 0) then
          write(unit,*) " "
          write(unit,*) "Post SCF Energies:"

          do i = 1, size(a_xc_arr)
             write(unit, '()')
             write(unit, 1000) "e_xc  ", a_xc_arr(i), trim(name_xc_arr(i))
             write(unit, 1000) "e_sum ", a_tot_arr(i), trim(name_xc_arr(i))
          end do
       endif

       ! FIXME: get rid of this,  this sub should only be printing. No
       ! other side effects:
       deallocate (name_xc_arr, a_xc_arr, a_tot_arr, stat=alloc_stat)
       ASSERT(alloc_stat==0)

       if (unit > 0) then
          write(unit,'()')

          if (operations_solvation_effect) then
             write(unit, 1000) "efr_el", a_energy_solv_el, 'FINAL'
             write(unit, 1000) "efrcdr", a_energy_solv_nonel, 'FINAL'
          endif

          ! only if DFT+U is in effect:
          if (e_dft_plus_u /= 0) then
             write(unit, 1000) "Edft+u", e_dft_plus_u
          endif

          !
          ! Here e_sum = a_tot + E(DFT-D), the latter is zero if DFT-D is
          ! off:
          !
          write(unit, 1000) "e_xc  ", a_xc, 'FINAL'
          write(unit, 1000) "e_sum ", e_sum, 'FINAL'

          if (operations_solvation_effect) then
             write (unit, *)
             write (unit, *) '===================================================================='
             write (unit, *) '        THE SOLVATION EFFECT - energy contributions'
             write (unit, *) '--------------------------------------------------------------------'

2000         format (a19, F25.12, " a.u. ", F10.2, " kcal/mol")

             write (unit, 2000) "Coulomb energy", a_energy_solv_el, a_energy_solv_el / kcal

             if (cavitation_energy) then
                write (unit, 2000) "Cavitation energy", energy_cav, energy_cav / kcal
             else
                ASSERT (energy_cav==0.0)
             endif

             if (disp_rep_energy) then
                write (unit, 2000) "Dispersion energy", E_disp, E_disp / kcal
                write (unit, 2000) "Repulsion energy", E_rep, E_rep / kcal
             else
                ASSERT (E_disp==0.0)
                ASSERT (E_rep==0.0)
             endif
             write (unit, 2000) "Sum", &
                  (a_energy_solv_el + energy_cav + E_disp + E_rep), &
                  (a_energy_solv_el + energy_cav + E_disp + E_rep) / kcal
             write (unit, *)
             write (unit, *) "The sum is missing electron reorganization energy."
             write (unit, *) "Use total energy differences instead."
          endif
       endif
    else
       if (unit > 0 ) then
          if (operations_solvation_effect) then
             write (unit, 1001) "efr_el", a_energy_solv_el, a_energy_solv_el - pre_sel
             if (cavitation_energy .or. disp_rep_energy) then
                write (unit, 1000) "efrcdr", a_energy_solv_nonel
             endif
          endif

          ! only if DFT+U is in effect:
          if (e_dft_plus_u /= 0)then
             write (unit, 1001) "Edft+u", e_dft_plus_u, e_dft_plus_u - pre_dft_plus_u
          endif

          write(unit,1001) "e_xc  ", a_xc, a_xc - pre_xc
          write(unit,1001) "e_el  ", e_sum - e_nuc_core, e_sum - e_nuc_core - pre_el
          write(unit,1001) "e_sum ", e_sum, e_sum - pre_tot
          write(unit,'()')
       endif

       if (unit == 6) then
          ! FIXME: this is  ugly. We print the same  energies first to
          ! output, then  on tty. The saved variables  for values from
          ! the previous iteration need to be updated once:
          pre_xc  = a_xc
          pre_el  = e_sum - e_nuc_core
          pre_tot = e_sum
          pre_sel = a_energy_solv_el
          pre_dft_plus_u = e_dft_plus_u
       endif
    end if
  end subroutine write_energies
  !*************************************************************

  !*************************************************************
  subroutine write_energy( io_unit                                             &
                         , print_kin                                           &
                         , print_nuc                                           &
                         , print_cou                                           &
                         , print_exx                                           )
    !
    ! Purpose: write the private variables to io_unit separately.
    !
    implicit none
    integer(kind=i4_kind), intent(in) :: io_unit
    logical, optional,     intent(in) :: print_kin
    logical, optional,     intent(in) :: print_nuc
    logical, optional,     intent(in) :: print_cou
    logical, optional,     intent(in) :: print_exx
    !** End of interface *****************************************

    1000 format(2X,A6,' =  ',F25.12:'  (',A,')')

    if (present(print_kin)) then
       if ( print_kin ) then
          !
          write(io_unit,1000)"e_kin ", e_kin
          !
       endif
    endif
    !
    if (present(print_nuc)) then
       if ( print_nuc ) then
          !
          write(io_unit,1000)"e_nuc ", e_nuc
          !
       endif
    endif
    !
    if (present(print_cou)) then
       if ( print_cou ) then
          !
          write(io_unit,1000)"e_coul", e_coul + e_coul_2z
          !
       endif
    endif
    !
    if (present(print_exx)) then
       if ( print_exx ) then
          !
          write(io_unit,1000)"e_exex ", e_exex
          !
       endif
    endif
    !
  end subroutine write_energy
  !*************************************************************

  !*************************************************************
  subroutine init_energy()
    !  Purpose: set the energies to zero
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)            :: zero = 0.0_r8_kind
    !------------ Executable code --------------------------------
    e_kin               = zero
    e_nuc               = zero
    e_efield            = zero
    e_coul              = zero
    e_coul_2z           = zero
    e_xc                = zero
    e_tot               = zero
    a_energy_solv_el    = zero
    a_energy_solv_nonel = zero
    e_rism = 0.0
  end subroutine init_energy

  subroutine set_energy (kin, nuc, coul, c_2z, exex, xc, &
       efield, dft_plus_u, rism)
    !
    ! Set energy contributions stored as the private variables in this
    ! module (for the  'direct_energy_calc' option). The variables set
    ! here  are   added  up  in  update_energies()   and  returned  by
    ! get_energy().
    !
    real(r8_kind), optional, intent(in) :: kin
    real(r8_kind), optional, intent(in) :: nuc
    real(r8_kind), optional, intent(in) :: xc
    real(r8_kind), optional, intent(in) :: coul
    real(r8_kind), optional, intent(in) :: c_2z
    real(r8_kind), optional, intent(in) :: exex
    real(r8_kind), optional, intent(in) :: efield
    real(r8_kind), optional, intent(in) :: dft_plus_u
    real(r8_kind), optional, intent(in) :: rism
    !** End of interface *****************************************


    if (present(kin) )   e_kin            = kin
    if (present(nuc) )   e_nuc            = nuc
    if (present(coul))   e_coul           = coul
    if (present(c_2z))   e_coul_2z        = c_2z
    if (present(exex))   e_exex           = exex
    if (present(xc)  )   e_xc             = xc
    if (present(efield)) e_efield         = efield
    if (present(dft_plus_u)) e_dft_plus_u = dft_plus_u
    if (present(rism)) e_rism = rism
  end subroutine set_energy

  subroutine energy_calc_reduce()
    !
    ! Sup up partial contributions to energies over workers. Note that
    ! some energies are only computed  by master, still we rely on the
    ! fact that in such case the slave contributions are proper zeros.
    !
    ! FIXME: used once in ham_calc_module.f90, usefulness
    ! questioned. Should it die?
    !
    use comm, only: comm_reduce, comm_rank, comm_same
    implicit none
    ! *** end of interface ***

    ! Some of the  energy contributions are computed by  master, so no
    ! need to reduce. To make sure  we get it right, verify this every
    ! time:
    if( comm_rank() /= 0 ) then
        ASSERT(e_kin==0.0)
        ASSERT(e_nuc==0.0)
        ASSERT(e_coul_2z==0.0)
        ASSERT(e_dft_plus_u==0.0)
        ASSERT(e_nuc_core==0.0)
        ASSERT(e_efield==0.0)
    endif

    ! Other energy contributions, such as e_rism, returned by external
    ! (serial or parallel) library with its own conventions are set on
    ! all  workers  to the  same  final value  and  thus  do not  need
    ! reduction either:
    ASSERT(comm_same(e_rism))

    ! There are still  some that need to be  summed up over processors
    ! and collected on master (this is not idempotent, of course):
    call comm_reduce(e_coul)
    call comm_reduce(e_xc)
    call comm_reduce(e_tot)
  end subroutine energy_calc_reduce

  subroutine direct_2c_energy_calc_and_add( f_t_arr, f_trace, t_arr            &
                                          , d_mat_real, d_mat_imag             &
                                          , e_2c, h_mat_real, h_mat_imag       )
    !--------------------------------------------------------------------------+
    !  Purpose: add 2-center terms (e.g. kin or nuc) to matrix ham
    !           and calculate expectation value.
    !           for spin orbit calculations the convention is to
    !           calculate
    !                          __
    !                         \  '
    !             tr{ P T } =  >   Re{ P  } Re{ T  } +  Im{ P  } Im{ T  }
    !                         /__,      ij       ij          ij       ij
    !                         k,l
    !
    !           whereas t_arr is supposed to be the lower triangular part
    !           of the matrix T
    !
    !------------ Declaration of interface variables --------------------------+
    ! prefactor for addition to hamiltonian matrix
    real(kind=r8_kind),           intent(in)    :: f_t_arr
    real(kind=r8_kind),           intent(in)    :: f_trace
    ! subsection of term of hamiltonian matrix
    real(kind=r8_kind),           intent(in)    :: t_arr(:)
    ! subsection of density matrix
    real(kind=r8_kind),           intent(in)    :: d_mat_real(:,:)
    real(kind=r8_kind), optional, intent(in)    :: d_mat_imag(:,:)
    ! subsection of hamiltonian matrix
    real(kind=r8_kind),           intent(inout) :: h_mat_real(:,:)
    real(kind=r8_kind), optional, intent(inout) :: h_mat_imag(:,:)
    !
    real(kind=r8_kind),           intent(inout) :: e_2c
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                       :: i_meta
    integer(kind=i4_kind)                       :: m
    integer(kind=i4_kind)                       :: n
    !------------ Executable code --------------------------------
    !
    if ( .not. present(h_mat_imag)) then
      !
      ! STANDARD SCF
      !
      i_meta = 0
      do m = 1, size(d_mat_real, 1)
        do n = 1, m - 1
          i_meta = i_meta + 1
          !
          ! OFF DIAGONAL ELEMENTS
          !
          ! calculate expectation value. Factor 2 due to treat upper triangle
          e_2c = e_2c + f_trace * d_mat_real(n,m) * t_arr(i_meta) * 2.0_r8_kind
          !
          ! add to hamiltonian matrix
          h_mat_real(n,m) = h_mat_real(n,m) + f_t_arr * t_arr(i_meta)
          h_mat_real(m,n) = h_mat_real(n,m)
          !
        enddo
        i_meta = i_meta + 1
        !
        ! DIAGONAL ELEMENTS
        !
        ! calculate expectation value
        e_2c = e_2c + f_trace * d_mat_real(m,m) * t_arr(i_meta)
        !
        ! add to hamiltonian matrix
        h_mat_real(m,m) = h_mat_real(m,m) + f_t_arr * t_arr(i_meta)
        !
      enddo
      !
    else
      !
      ! SPIN ORBIT
      !
      ASSERT(present(d_mat_imag))
      !
      i_meta = -1
      do m = 1, size(d_mat_real, 1)
        do n = 1, m - 1
          i_meta = i_meta + 2
          !
          ! OFF DIAGONAL ELEMENTS
          !
          ! calculate expectation value. Factor 2 due to treat upper triangle
          e_2c = e_2c + f_trace * ( d_mat_real(n,m) * t_arr(i_meta)            &
                                  + d_mat_imag(n,m) * t_arr(i_meta+1))         &
                                * 2.0_r8_kind
          !
          ! add to hamiltonian matrix
          h_mat_real(n,m) =   h_mat_real(n,m) + f_t_arr * t_arr(i_meta)
          h_mat_real(m,n) =   h_mat_real(n,m)
          h_mat_imag(n,m) =   h_mat_imag(n,m) - f_t_arr * t_arr(i_meta+1)
          h_mat_imag(m,n) = - h_mat_imag(n,m)
          !
        enddo
        i_meta = i_meta + 2
        !
        ! DIAGONAL ELEMENTS
        !
        ! calculate expectation value
        e_2c = e_2c + f_trace * ( d_mat_real(m,m) * t_arr(i_meta)              &
                                + d_mat_imag(m,m) * t_arr(i_meta+1)            )
        !
        ! add to hamiltonian matrix
        h_mat_real(m,m) =   h_mat_real(m,m) + f_t_arr * t_arr(i_meta)
        h_mat_imag(m,m) =   h_mat_imag(m,m) + f_t_arr * t_arr(i_meta+1)
        !
      enddo
      !
    endif
    !
  end subroutine direct_2c_energy_calc_and_add
  !*************************************************************

  !*************************************************************
  subroutine direct_mat_energy_calc_and_add( f_trace, f_t_mat                  &
                                          , t_mat_real, t_mat_imag             &
                                          , d_mat_real, d_mat_imag             &
                                          , e_tr, h_mat_real, h_mat_imag       )
    !  Purpose: add matrix term t_mat to matrix ham
    !           and calculate expectation value
    !------------ Declaration of interface variables -------------
    ! prefactor for addition to hamiltonian matrix
    real(kind=r8_kind),           intent(in)    :: f_trace
    real(kind=r8_kind),           intent(in)    :: f_t_mat
    ! subsection of term of hamiltonian matrix
    real(kind=r8_kind),           intent(in)    :: t_mat_real(:,:)
    real(kind=r8_kind), optional, intent(in)    :: t_mat_imag(:,:)
    ! subsection of density matrix
    real(kind=r8_kind),           intent(in)    :: d_mat_real(:,:)
    real(kind=r8_kind), optional, intent(in)    :: d_mat_imag(:,:)
    ! subsection of hamiltonian matrix
    real(kind=r8_kind),           intent(inout) :: h_mat_real(:,:)
    real(kind=r8_kind), optional, intent(inout) :: h_mat_imag(:,:)
    !
    real(kind=r8_kind),           intent(inout) :: e_tr
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                       :: m
    !------------ Executable code --------------------------------
    !
    if ( .not. present(h_mat_imag)) then
      !
      ! STANDARD SCF
      !
      do m = 1, size(d_mat_real, 1)
        !
        ! calculate expectation value
        e_tr = e_tr + f_trace * sum(d_mat_real(:,m) * t_mat_real(:,m))
        !
        ! add to hamiltonian matrix
        h_mat_real(:,m) = h_mat_real(:,m) + f_t_mat * t_mat_real(:,m)
        !
      enddo
      !
    else
      !
      ! SPIN ORBIT
      !
      ASSERT(present(d_mat_imag))
      ASSERT(present(t_mat_imag))
      !
      do m = 1, size(d_mat_real, 1)
        !
        ! OFF DIAGONAL ELEMENTS
        !
        ! calculate expectation value
        e_tr = e_tr + f_trace * sum( d_mat_real(:,m) * t_mat_real(:,m)         &
                                   - d_mat_imag(:,m) * t_mat_imag(:,m)         )
        !
        ! add to hamiltonian matrix
        h_mat_real(:,m) =   h_mat_real(:,m) + f_t_mat * t_mat_real(:,m)
        h_mat_imag(:,m) =   h_mat_imag(:,m) + f_t_mat * t_mat_imag(:,m)
        !
      enddo
      !
    endif
    !
  end subroutine direct_mat_energy_calc_and_add
  !*************************************************************

  !--------------- End of module ----------------------------------
end module energy_calc_module

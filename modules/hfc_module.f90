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
Module hfc_module
  !-------------------------------------------------------------------
  !
  !  Purpose: database to hold hfc matrix elements
  !           as well as calculated hfc matrix elements.
  !
  !
  !  The data in this module only exist on the master !!!
  !
  !  Module called by: ...
  !
  !  References: ...
  !  Author: DG
  !  Date: 1/2003
  !
  !
  ! Modification
  ! Author:
  ! Date:
  ! Description:
  !
  !-------------------------------------------------------------------
#include "def.h"
#define DEBUG
      Use type_module ! type specification parameters
      Use options_module, Only: options_spin_orbit, &
     & options_kinematic_factors, options_finite_nucleus, options_relativistic
      Use dipole_module
      Use unique_atom_module
      Use symmetry_data_module
      Use eigen_data_module
      Use occupation_module
      Use filename_module
      Use iounitadmin_module
      Use dimensions, Only: n_irr => number_of_irreps, IrrUDim => &
     & IrrUBasDimSpor, IrrCDim => IrrBasDimSpor ! 1 means one-component

       Use contraction_module, Only: contract_matrix_spor => contract


      Implicit None
      Save ! save all variables defined in this module
      Private ! by default, all names are private

       Type, public:: g_matrix_element
       !DG type g_matrix_element
       !Contain data consernied to calculated g-tensors, matrix elements <PHI(i)|(dH^Z/dB^0(k))|PHI(j)>
       ! Sigma_Real, Sigma_Im is matrix elements <|Sigma|>
       ! L_Real, L_Im is matrix elements <|L|>
!
         Integer (Kind=i4_kind) :: i_irr, i_ip, i_orb, n_orb, i_pa
         Complex(c16_kind)      :: vecSigma11(3), vecSigma12(3)
         ! matrix elements of VecSigmaXY = <X| vecSigma |Y>
         ! the Hamiltonian is built of scalar product:
         ! | vecSigma11,   vecSigma12 |
         ! | vecSigma12*, -vecSigma11 |
         ! with manetic field vecB
         Complex(c16_kind)      :: vecL11(3), vecL12(3)
         ! see vecSigmaXY
         Real (Kind=r8_kind), Dimension (3) :: Sigma_Real_diag, L_Real_diag, Sigma_Im_diag, L_Im_diag
         Real (Kind=r8_kind), Dimension (3) :: Sigma_Real_offdiag, &
              & Sigma_Im_offdiag, L_Real_offdiag, L_Im_offdiag

         Real (Kind=r8_kind) ::  D_real_diag, D_Real_offdiag, D_Im_offdiag
         Real (Kind=r8_kind) :: eigenvalue
         Real (Kind=r8_kind) :: norm, norm_offdiag
      End Type g_matrix_element

      Type, public :: properties_integrals_type
         Type (dipole_integral_type),pointer :: integrals(:,:,:) ! one less
      end Type properties_integrals_type
       Type (properties_integrals_type), Public :: hfc_integrals_diag,  hfc_integrals_offdiag
       Type (properties_integrals_type), Private:: hfc_integrals_rem,  hfc_integrals_cont_rem,hfc_integrals_cont_offdiag
       Type (dipole_integral_type), Allocatable, Target, Public :: &
     & hfc_integrals (:, :, :, :)!,hfc_integrals_diag(:,:,:,:), hfc_integrals_offdiag(:,:,:,:)
      Type (dipole_integral_type), Allocatable, Target ::  hfc_integrals_cont (:, :, :, :) !contracted
       ! Default valuse of the input variables
       Integer (Kind=i4_kind), Pointer, Dimension (:) :: dim_irr_arr, &
            & dim_irr_c_arr, dimensions_arr
       Type (g_matrix_element) a_mx_el
       Type (g_matrix_element) iso_mx_el
       Real (Kind=r8_kind), Parameter :: PI = 3.1415926535897932368
       Real (Kind=r8_kind), Parameter :: Mp = 1836.12
       Real (Kind=r8_kind), Parameter :: ge = 2.00231929 !value from Jone E. Harriman
       Real (Kind=r8_kind), Parameter :: speed_of_light =  137.03604
       Real (Kind=r8_kind), Parameter :: coeff_au = 2 * PI * ge / (Mp*speed_of_light*speed_of_light*3)
       Real (Kind=r8_kind), Parameter :: K_convert_to_MHz = &
            & 6.57968374E+09
       Real (Kind=r8_kind), Parameter :: K_convert_to_Gauss = &
            & 2.347828566E+09
  !  "Theoretical Foundations of Electron Spin Resonance"
        Real (Kind=r8_kind), Parameter :: very_small = 0.000001 ! precision parameter
       Logical, Parameter, Public :: CONTRACTED = .false.
       Logical, Parameter, Public :: UNCONTRACTED = .true.
        Integer (Kind=i4_kind), Public :: n_exclude_atom_orbitals,n_orbitals, dkh_level
  Real (Kind=r8_kind), Public :: mix_coeff
  Logical, Public :: read_from_list
  Logical, Public :: exclude_atom_orbitals, include_key
  Namelist / magnetic_properties /&
       & read_from_list,&
       & n_orbitals, &
       & exclude_atom_orbitals,&
       & include_key,&
       & n_exclude_atom_orbitals, &
       & mix_coeff,&
       & dkh_level ! DKH level for magnetic operator
  Integer (Kind=i4_kind), Allocatable,Public :: input_irreps (:), &
       & input_partners (:), input_orbitals (:)
  Integer (Kind=i4_kind), Allocatable, Public :: &
       & input_exclude_unique_atoms_1 (:), input_exclude_atom_l_1 (:), &
       & input_exclude_atom_n_1 (:), input_exclude_unique_atoms_2 (:), &
       & input_exclude_atom_l_2 (:), input_exclude_atom_n_2 (:)
  Public   magnetic_read, magnetic_write
  Logical :: df_read_from_list = .False., df_exclude_atom_orbitals &
       & = .False., df_include_key = .False.
  Integer (Kind=i4_kind) :: df_n_orbitals = 0, &
       & df_n_exclude_atom_orbitals = 0
  Integer (Kind=i4_kind) :: df_dkh_level = 1
  Real (Kind=r8_kind) :: df_mix_coeff = 0.5_r8_kind
!== Interrupt end of public interface of module =================

 interface load_array
    module procedure load_matrix
    module procedure load_vector
 end interface

interface new
   module procedure init_hfc_integrals
end interface

interface delete
   module procedure null_hfc_integrals
end interface

Public :: &
     hfc_init, & ! FIXME: called from hfc_allocate AND gten_allocate
     hfc_done, & ! FIXME: called from hfc_free     AND gten_free
     hfc_calculate, &
     hfc_allocate, &
     hfc_free, &
     load_array, &
     new, &
     delete, &
     remap_to_full_store, &
     remap_to_compact_store, &
     magnetic_properties_shutdown

Contains

  subroutine hfc_init()
    implicit none
    ! *** end of interface ***

    DPRINT 'hfc_init: call dimensions_calculate()'
    call dimensions_calculate()
  end subroutine hfc_init

  subroutine hfc_done()
    implicit none
    ! *** end of interface ***

    Integer (Kind=i4_kind) :: status

    DPRINT 'hfc_init: deallocate dimensions'
    deallocate (dim_irr_arr, dim_irr_c_arr, Stat = status)
    ASSERT(status==0)
  end subroutine hfc_done

  subroutine dimensions_calculate()
    implicit none
    ! *** end of interface ***

    !------------ Declaration of local variables ---------------------
    Integer (Kind=i4_kind) :: i_irrep, i_unique, i_l, n_irreps
    Integer (Kind=i4_kind) :: status

    !------------ Executable code ------------------------------------

    DPRINT 'hfc::dimensions_calculate: entered'

    if(options_spin_orbit) then
       n_irreps = symmetry_data_n_proj_irreps()
    else
       n_irreps = symmetry_data_n_irreps ()
    endif
    allocate( &
         dim_irr_arr(n_irreps), &
         dim_irr_c_arr(n_irreps), &
         Stat = status)
    ASSERT(status==0)

    if (.not. options_spin_orbit) then
       ! calculate uncontructed dimensions
       do i_irrep = 1, n_irreps
          dim_irr_arr(i_irrep) = 0
          do i_unique=1,n_unique_atoms
             do i_l=0,unique_atoms(i_unique)%lmax_ob
                dim_irr_arr(i_irrep)=dim_irr_arr(i_irrep)+unique_atoms(i_unique)%&
                     symadapt_partner(i_irrep,i_l)%n_independent_fcts*&
                     unique_atoms(i_unique)%l_ob(i_l)%n_exponents
             end do
          end do
       end do
       ! calculate_contracted_dimensions
       do i_irrep = 1, n_irreps
          dim_irr_c_arr(i_irrep) = symmetry_data_dimension(i_irrep)
       end do
    else !options_spin_orbit
       DPRINT 'hfc::dimensions_calculate: IrrUDim=',IrrUDim
       DPRINT 'hfc::dimensions_calculate: IrrCDim=',IrrCDim
       dim_irr_arr   = IrrUDim(:)
       dim_irr_c_arr = IrrCDim(:)
    end if
  end subroutine dimensions_calculate

  Subroutine magnetic_read
    ! purpose: read namelist for gtensors
    !
    ! routine called by: read_input
    !** End of interface *****************************************
    Use input_module
    !------------ Declaration of local variables ---------------------
    Integer (Kind=i4_kind) :: unit, status, counter
    External error_handler
    !------------ Executable code ------------------------------------
    read_from_list = df_read_from_list
    n_orbitals = df_n_orbitals
    exclude_atom_orbitals = df_exclude_atom_orbitals
    include_key = df_include_key
    n_exclude_atom_orbitals = df_n_exclude_atom_orbitals
    mix_coeff = df_mix_coeff
    dkh_level = df_dkh_level
    If (input_line_is_namelist("magnetic_properties")) Then
       Call input_read_to_intermediate ()
       unit = input_intermediate_unit ()
       Read (Unit, Nml=magnetic_properties, IoStat=Status)
       If (status .Gt. 0) Call input_error ("occupation_read: name&
            &list magnetic_properties.")
       !
       If (read_from_list) Then
          Allocate (input_irreps(n_orbitals), Stat=status)
          If (status .Ne. 0) Call error_handler ("hfc_module:&
               & allocate of input_irreps")
          !
          Allocate (input_partners(n_orbitals), Stat=status)
          If (status .Ne. 0) Call error_handler ("hfc_module:&
               & allocate of input_partners")
          !
          Allocate (input_orbitals(n_orbitals), Stat=status)
          If (status .Ne. 0) Call error_handler ("hfc_module:&
               & allocate of input_orbitals failed")
          Do counter = 1, n_orbitals
             Call input_read_to_intermediate
             unit = input_intermediate_unit ()
             Read (Unit, Fmt=*, IoStat=Status) input_irreps &
                  & (counter), input_partners (counter), input_orbitals &
                  & (counter)
             If (status > 0) Call input_error ("magnetic: n&
                  &oc_start")
          End Do
       End If
       !
       If (exclude_atom_orbitals) Then
          Allocate &
               & (input_exclude_unique_atoms_1(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          !
          Allocate &
               & (input_exclude_atom_l_1(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          Allocate &
               & (input_exclude_atom_n_1(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          !
          Allocate &
               & (input_exclude_unique_atoms_2(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          !
          Allocate &
               & (input_exclude_atom_l_2(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          Allocate &
               & (input_exclude_atom_n_2(n_exclude_atom_orbitals), &
               & Stat=status)
          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
               & allocate of dipole_integrals failed")
          !
          !
          Do counter = 1, n_exclude_atom_orbitals
             Call input_read_to_intermediate
             unit = input_intermediate_unit ()
             Read (Unit, Fmt=*, IoStat=Status) &
                  & input_exclude_unique_atoms_1 (counter), &
                  & input_exclude_atom_l_1 (counter), &
                  & input_exclude_atom_n_1 (counter), &
                  & input_exclude_unique_atoms_2 (counter), &
                  & input_exclude_atom_l_2 (counter), &
                  & input_exclude_atom_n_2 (counter)
             If (status > 0) Call input_error ("occupation_read: n&
                  &oc_start")
          End Do
       End If
       !
    End If
  End Subroutine magnetic_read

  subroutine magnetic_properties_shutdown
    Integer (Kind=i4_kind) ::status

    Deallocate (input_irreps, Stat=status)
    If (status /= 0) Call error_handler ('gtensor_calculate: dellocation of input_irreps failed')
    Deallocate (input_partners, Stat=status)
    If (status /= 0) Call error_handler ('gtensor_calculate: dellocation of input_irreps failed')
    Deallocate (input_orbitals, Stat=status)
    If (status /= 0) Call error_handler ('gtensor_calculate: dellocation of input_irreps failed')
  end subroutine magnetic_properties_shutdown
  !*********************************************
  Subroutine magnetic_write (iounit)
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full
    Use operations_module, Only: operations_echo_input_level
    implicit none
    Integer, Intent (In) :: iounit
    ! *** end of interface ***

    Call start ("MAGNETIC_PROPERTIES", "MAGNETIC_WRITE", iounit, &
         & operations_echo_input_level)
    Call flag ("READ_FROM_LIST ", read_from_list, &
         & df_read_from_list)
    Call intg ("N_ORBITALS     ", n_orbitals, df_n_orbitals)
    Call flag ("EXCLUDE_ATOM_ORBITALS   ", exclude_atom_orbitals, &
         & df_exclude_atom_orbitals)
    Call flag ("INCLUDE_KEY   ", include_key, df_include_key)
    Call intg ("N_EXCLUDE_ATOM_ORBITALS ", &
         & n_exclude_atom_orbitals, df_n_exclude_atom_orbitals)
    Call real ("MIX_COEFF ", mix_coeff, df_mix_coeff)
    !
    Call stop ()
  End Subroutine magnetic_write


Subroutine reltrafo_hfc
  !  Purpose: makes relativistic transformations of magnetic matrix elements
  !
  !
  !
  !
  !
  !
  !  Author: DG
  !  Date:   6/2003
  !
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  Use symmetry_data_module
  use matrix_module, only: &
        pgfree, tr, sim, eigs, alloc, cmatrix, rdmatrix, rmatrix, &
        operator(*), operator(**), assignment(=), operator(-), &
        operator(+), mult
  use reltrafo, only: p2_diag

   Implicit None

  Integer (i4_kind) :: dim_irrep, n_irreps, n_partners, i_ir,i_ip, i_xyz
  Integer (i4_kind) :: i_partn, n_orb
   Real (Kind=r8_kind) ::c2, c4, deltaE
   Real (Kind=r8_kind), parameter :: half = 0.5_r8_kind
   Integer (i4_kind) :: i,j,i_ua,n,m, diagoffdiag
   Type (dipole_integral_type), Pointer :: hfc_p, hfc_cont_p, hfc_tmp
   type(cmatrix)  :: UF,UB, Nuc, PVYP, AVA, ARVRA,ARVRAt, DKH1, hpp, GTI, DKH2, WW
   type(rdmatrix) :: t_diag,Tp, Ep,Ap,Kp,ApKp,K2p2, k2p2_rev

!-----------executable code---------------------
   nullify (hfc_p,hfc_cont_p)

  c2 = speed_of_light ** 2
  c4 = speed_of_light ** 4
  !dim_irr_arr => IrrUDim(:)
  !dim_irr_c_arr => IrrCDim(:)
  kf: If (options_kinematic_factors) Then
!
     n_irreps = symmetry_data_n_proj_irreps ()
     i_ip=1
               cycle_kf: Do i_ir = 1, n_irreps

                  dim_irrep = dim_irr_arr(i_ir)
                  n_orb = dim_irr_c_arr (i_ir)
!
!

                  call alloc(dim_irrep,UF,UB)

                  call load_array(dim_irrep,UF%re, i_ir,"forward_real" )

                  call load_array(dim_irrep,UF%im, i_ir,"forward_imag" )

                  call load_array(dim_irrep,UB%re, i_ir,"U_backward_real" )

                  call load_array(dim_irrep,UB%im, i_ir,"U_backward_imag" )

                  call alloc (dim_irrep, t_diag)

                  call load_array(dim_irrep,t_diag%d, i_ir,"gt_t_diag" )

                  call p2_diag(t_diag, Tp, Ep, Ap, Kp, ApKp, K2p2)

                  call alloc(dim_irrep, Nuc, PVYP, AVA, ARVRA,ARVRAt,DKH1, hpp, GTI, WW,DKH2)

                  ! call alloc(dim_irrep, Tpp,Upp)

                  call load_array(dim_irrep,Nuc%re,i_ir,"gt_Nuc_real")
                  call load_array(dim_irrep,Nuc%im,i_ir,"gt_Nuc_imag")
                  call load_array(dim_irrep,PVYP%re,i_ir,"gt_PVYP_real")
                  call load_array(dim_irrep,PVYP%im,i_ir,"gt_PVYP_imag")

                  AVA   =  mult(Ap,Nuc,Ap)
                  ARVRA =  mult(ApKp,PVYP,ApKp)
                  K2p2_rev = k2p2**(-1)
                  do j=1,dim_irrep
                     do i=1,dim_irrep
                        DeltaE = Ep%d(i)+Ep%d(j)
                        AVA%re(i,j)   = AVA%re(i,j)/DeltaE
                        AVA%im(i,j)   = AVA%im(i,j)/DeltaE

                        ARVRA%re(i,j) = ARVRA%re(i,j)/DeltaE
                        ARVRA%im(i,j) = ARVRA%im(i,j)/DeltaE
                     enddo
                  enddo

                  ARVRAt%re = ARVRA%re
                  ARVRAt%im = - ARVRA%im

     !             do j=1,dim_irrep
     !                do i=1,dim_irrep
     !             Tpp%re(i,j) = K2p2%d(i) * AVA%re(i,j) - ARVRA%re(i,j)
     !             Tpp%im(i,j) = K2p2%d(i) * AVA%im(i,j) - ARVRA%im(i,j)
     !
     !             Upp%re(i,j) = AVA%re(i,j) - ARVRA%re(i,j) * (1_r8_kind/(K2p2%d(j)))
     !             Upp%im(i,j) = AVA%im(i,j) - ARVRA%im(i,j) * (1_r8_kind/(K2p2%d(j)))
     !          end do
     !       end do

                  ! Efficient way
                  !  Tpp = K2p2 * AVA - ARVRA
                  !
                  !  Upp = AVA - ARVRA * K2p2_rev

                  n_partners = symmetry_data_n_partners_proj (i_ir)
                npart: Do i_partn = 1, n_partners
                  xyz: Do i_xyz = 1, 9
                     nua:   Do i_ua = 1, N_unique_atoms
                        diag: Do diagoffdiag = 1, 2
                           if (diagoffdiag == 1 )then
                              hfc_tmp => hfc_integrals_diag%integrals (i_xyz, i_ip, i_ua)
                              hfc_p => hfc_integrals_rem%integrals (i_xyz, i_ip, i_ua)
                              hfc_cont_p => hfc_integrals_cont_rem%integrals(i_xyz, i_ip, i_ua)
                              !  Diagonal part

                              Call remap_to_full_store (hfc_tmp,hfc_p)

                              GTI%re(:,:) = hfc_p%offdiagonal(:,:)
                              GTI%im(:,:) = hfc_p%offdiagonal_imag(:,:)
                           else
                              hfc_p => hfc_integrals_offdiag%integrals (i_xyz, i_ip, i_ua)
                              hfc_cont_p => hfc_integrals_cont_offdiag%integrals(i_xyz, i_ip, i_ua)

                              GTI%re(:,:) = hfc_p%offdiagonal(:,:)
                              GTI%im(:,:) = hfc_p%offdiagonal_imag(:,:)
                           end if

                           !transormations into momentum space
                           GTI = tr(UF) * GTI * UF

                           !DKH1 transformation
                           DKH1 = ApKp * GTI * Ap + Ap * GTI * ApKp

                           DKH1 = DKH1 * speed_of_light

                           if (dkh_level == 1) then

                              GTI = DKH1

                           elseif (dkh_level == 2) then
                              !DKH2 transformation
                              !init


                              Do n = 1, dim_irrep
                                 Do m = 1, dim_irrep
                                    hpp%re(m,n)= ApKp%d(m)*GTI%re(m,n)*Ap%d(n)/(Ep%d(m)+Ep%d(n))
                                    hpp%im(m,n)= ApKp%d(m)*GTI%im(m,n)*Ap%d(n)/(Ep%d(m)+Ep%d(n))
                                 end Do
                              end Do




                  !  hpp =( - half * mult((mult(Upp,Ep)+mult(Ep,Upp)), hpp) &
                  !         - half * mult(Upp , mult(hpp,Ep)+mult(Ep,hpp)) &
                  !         - half * mult(tr(hpp) , mult(tr(Upp),Ep) + mult(Ep,tr(Upp))) &
                  !         - half * mult(mult(tr(hpp),Ep)+mult(Ep,tr(hpp)),tr(Upp)) &
                  !         + half * mult(hpp , mult(Tpp, Ep) + mult(Ep,Tpp)) &
                  !        + half * mult(mult(hpp,Ep)+mult(Ep,hpp),Tpp) &
                  !         + half * mult(mult(tr(Tpp),Ep) + mult(Ep,tr(Tpp)), tr(hpp)) &
                  !         + half * mult(tr(Tpp),mult(tr(hpp),Ep) + mult(Ep,tr(hpp)))) * speed_of_light

                              WW = - AVA * hpp - hpp * AVA + AVA * k2p2 * hpp + hpp * k2p2 * AVA
                              WW = WW * Ep + Ep * WW
                              DKH2 = WW + 2.0_r8_kind *&
                                   ( - AVA * Ep * hpp - hpp * Ep * AVA + AVA * Ep * k2p2 * hpp + hpp * Ep * k2p2 * AVA)


                              if ( diagoffdiag == 1 ) then
                                 WW = ARVRA * k2p2_rev * hpp + hpp * k2p2_rev * ARVRA - ARVRA * hpp - hpp * ARVRA
                                 WW = WW * Ep + Ep * WW
                                 WW = WW + 2.0_r8_kind * (ARVRA * k2p2_rev * Ep * hpp + hpp * k2p2_rev * Ep * ARVRA &
                                      - ARVRA * Ep * hpp - hpp * Ep * ARVRA)

                              else
                                 WW = ARVRAt * k2p2_rev * hpp + hpp * k2p2_rev * ARVRA - ARVRAt * hpp - hpp * ARVRA
                                 WW = WW * Ep + Ep * WW
                                 WW = WW + 2.0_r8_kind * (ARVRAt * k2p2_rev * Ep * hpp + hpp * k2p2_rev * Ep * ARVRA &
                                      - ARVRAt * Ep * hpp - hpp * Ep * ARVRA)
                              end if

                               GTI = - (DKH2 + WW) * speed_of_light * half + DKH1


                           end if

                   GTI = tr(UB) * GTI * UB


                   Call contract_matrix_spor(i_ir, GTI%re, hfc_cont_p%offdiagonal)


                   Call contract_matrix_spor(i_ir, GTI%im,hfc_cont_p%offdiagonal_imag)

                     if(diagoffdiag==1)   Call remap_to_compact_store (hfc_cont_p, hfc_tmp, n_orb)

                  end Do diag
               end Do nua
            end Do xyz
!
            i_ip=i_ip+1
         end Do npart

!

              call pgfree(UF,UB)
              call pgfree(Nuc,PVYP, AVA, ARVRA, ARVRAt,WW, DKH1, hpp, DKH2)
              call pgfree (t_diag, Ep, Tp,Ap, Kp, ApKp,K2p2, k2p2_rev)
           !   call pgfree (Upp,Tpp)

               End Do cycle_kf
            End If kf

end Subroutine reltrafo_hfc

Subroutine reltrafo_hfc_one_comp
  !  Purpose: makes relativistic transformations of magnetic matrix elements
  !
  !
  !
  !
  !
  !
  !  Author: DG
  !  Date:   8/2004
  !  A one component version
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  Use symmetry_data_module
  use matrix_module, only: &
        pgfree, tr, sim, eigs, alloc, cmatrix, rmatrix, rdmatrix, &
        operator(*), operator(**), assignment(=), operator(-), &
        operator(+), mult
  use reltrafo, only: p2_diag
  use contraction_module, only: contract_matrix=>contract

   Implicit None

  Integer (i4_kind) :: dim_irrep, n_irreps, n_partners, i_ir,i_ip, i_xyz
  Integer (i4_kind) :: i_partn, dim_c_irrep
   Real (Kind=r8_kind) ::c2, c4, deltaE
   Real (Kind=r8_kind), parameter :: half = 0.5_r8_kind
   Integer (i4_kind) :: i,j,i_ua,n,m
   Type (dipole_integral_type), Pointer :: hfc_p, hfc_cont_p, hfc_tmp
   type(rmatrix)  :: UF,UB, Nuc, PVYP, AVA, ARVRA,ARVRAt, DKH1, hpp, GTI, DKH2, WW
   type(rdmatrix) :: t_diag,Tp, Ep,Ap,Kp,ApKp,K2p2, k2p2_rev

!-----------executable code---------------------
   nullify (hfc_p,hfc_cont_p)

  c2 = speed_of_light ** 2
  c4 = speed_of_light ** 4
!!$  dim_irr_arr => IrrUDim(:)
!!$  dim_irr_c_arr => IrrCDim(:)
  kf: If (options_kinematic_factors) Then
!
     n_irreps = symmetry_data_n_irreps ()
     i_ip=1
               cycle_kf: Do i_ir = 1, n_irreps

                !  dim_irrep = dim_irr_arr(i_ir)
                !  n_orb = dim_irr_c_arr (i_ir)
                   dim_irrep = dim_irr_arr(i_ir)
                   dim_c_irrep = dim_irr_c_arr(i_ir) ! DG to be modified

                  call alloc(dim_irrep,UF)
                  call alloc(dim_irrep,UB)

                  call load_array(dim_irrep,UF%m, i_ir,"u_forward" )

                  call load_array(dim_irrep,UB%m, i_ir,"u_backward" )


                  call alloc (dim_irrep, t_diag)

                  call load_array(dim_irrep,t_diag%d, i_ir,"hfc_t_diag" )

                  call p2_diag(t_diag, Tp, Ep, Ap, Kp, ApKp, K2p2)

                  call alloc(dim_irrep, Nuc)
                   call alloc(dim_irrep, PVYP)
                    call alloc(dim_irrep, AVA)
                     call alloc(dim_irrep, ARVRA)
                      call alloc(dim_irrep, ARVRAt)
                       call alloc(dim_irrep, DKH1)
                        call alloc(dim_irrep, hpp)
                         call alloc(dim_irrep, GTI)
                          call alloc(dim_irrep, WW )
                           call alloc(dim_irrep, DKH2 )


                  call load_array(dim_irrep,Nuc%m,i_ir,"hfc_nuc")

                  call load_array(dim_irrep,PVYP%m,i_ir,"hfc_pvsp")

                  !test
                !  Nuc%m  = 0.0_r8_kind
                !  PVYP%m = 0.0_r8_kind
!!$
!!$
                  AVA  = mult(Nuc, Ap) !AVA   =  mult(Ap,Nuc,Ap)
                  AVA  = mult(Ap, AVA)

                  ARVRA = mult(PVYP, ApKp) ! ARVRA =  mult(ApKp,PVYP,ApKp)
                  ARVRA = mult(ApKp, ARVRA)

                  K2p2_rev = k2p2**(-1)
                  do j=1,dim_irrep
                     do i=1,dim_irrep
                        DeltaE = Ep%d(i)+Ep%d(j)
                        AVA%m(i,j)   = AVA%m(i,j)/DeltaE


                        ARVRA%m(i,j) = ARVRA%m(i,j)/DeltaE

                     enddo
                  enddo

                  ARVRAt%m = ARVRA%m



                  ! Efficient way
                  !  Tpp = K2p2 * AVA - ARVRA
                  !
                  !  Upp = AVA - ARVRA * K2p2_rev

                  n_partners = symmetry_data_n_partners (i_ir)
                npart: Do i_partn = 1, n_partners
                  xyz: Do i_xyz = 1, 7 !
                     nua:   Do i_ua = 1, N_unique_atoms


                              ! triangular uncontracted matrix:
                              hfc_tmp => hfc_integrals_diag%integrals (i_xyz, i_ip, i_ua)

                              ! remapped from triangular form to square form matrix:
                              hfc_p => hfc_integrals_rem%integrals (i_xyz, i_ip, i_ua)

                              ! contacted square matrix:
                              hfc_cont_p => hfc_integrals_cont_rem%integrals(i_xyz, i_ip, i_ua)
                              !  Diagonal part

                              Call remap_to_full_store (hfc_tmp,hfc_p)

                              GTI%m(:,:) = hfc_p%offdiagonal(:,:)



                           !transormations into momentum space
                           !GTI = tr(UF) * GTI * UF
                           GTI = mult(GTI, UF)
                           GTI = mult(tr(UF), GTI)

                           !DKH1 transformation
                           !DKH1 = ApKp * GTI * Ap + Ap * GTI * ApKp
                            DKH1 = mult(ApKp, mult(GTI, Ap)) + mult(Ap, mult(GTI, ApKp))

                           DKH1 = DKH1 * speed_of_light

                           if (dkh_level == 1) then

                              GTI = DKH1

                           elseif (dkh_level == 2) then
                              !DKH2 transformation
                              !init


                              Do n = 1, dim_irrep
                                 Do m = 1, dim_irrep
                                    hpp%m(m,n)= ApKp%d(m)*GTI%m(m,n)*Ap%d(n)/(Ep%d(m)+Ep%d(n))

                                 end Do
                              end Do





                              WW = - mult(AVA, hpp) - mult(hpp, AVA) + mult(AVA, mult(k2p2, hpp)) + mult(hpp, mult(k2p2, AVA))
                              WW = mult(WW, Ep) + mult(Ep, WW)
                              DKH2 = WW + 2.0_r8_kind *&
                                   ( - mult(AVA, mult(Ep, hpp)) - mult(hpp, mult(Ep, AVA)) +&
                                   mult(AVA, mult(Ep, mult(k2p2, hpp)))+&
                                   mult(hpp, mult(Ep, mult(k2p2, AVA))))



                                 WW = mult(ARVRA, mult(k2p2_rev, hpp)) + mult(hpp, mult(k2p2_rev, ARVRA)) - mult(ARVRA, hpp)&
                                      - mult(hpp, ARVRA)
                                 WW = mult(WW, Ep) + mult(Ep, WW)
                                 WW = WW + 2.0_r8_kind * (mult(ARVRA, mult( k2p2_rev, mult(Ep, hpp)))&
                                      + mult(hpp, mult(k2p2_rev, mult(Ep, ARVRA))) &
                                      - mult(ARVRA, mult(Ep, hpp)) - mult(hpp, mult(Ep, ARVRA)))



                            !   GTI = - (DKH2 + WW) * speed_of_light * half + DKH1
                                 GTI = (DKH2 + WW) * speed_of_light * half + DKH1


                            end if

                   GTI = tr(UB) * GTI * UB


                   call  contract_matrix(i_ir, GTI%m, hfc_cont_p%offdiagonal)



                   Call remap_to_compact_store (hfc_cont_p, hfc_tmp, dim_c_irrep)

                end Do nua

             end Do xyz
!
             i_ip=i_ip+1
          end Do npart

!

              call pgfree(UF)
              call pgfree(UB)
              call pgfree(Nuc)
              call pgfree(PVYP)
              call pgfree(AVA)
              call pgfree(ARVRA)
              call pgfree(ARVRAt)
              call pgfree(WW)
              call pgfree(DKH1)
              call pgfree(hpp)
              call pgfree(DKH2)
              call pgfree (t_diag, Ep, Tp,Ap, Kp, ApKp,K2p2, k2p2_rev)


               End Do cycle_kf
            End If kf

          end Subroutine reltrafo_hfc_one_comp

Subroutine init_hfc_integrals(this,n_xyz,n_ua,usedoff,usedim,components)
Type (properties_integrals_type) :: this
Integer (Kind=i4_kind), intent(in) :: n_xyz, usedoff, n_ua
Logical, intent(in) :: usedim
Integer (Kind=i4_kind) :: Status, n_ip, n_ir, i_ip,i_ir,&
     i_pa, i_ua, dim_ir, i_xyz, components
!------------------------------------

      If (usedim) Then
         dimensions_arr => dim_irr_arr(:)
      Else
         dimensions_arr => dim_irr_c_arr(:)
      End If

select case (components)

   case (2)
      n_ip = symmetry_data_n_ip_proj ()
      n_ir = symmetry_data_n_proj_irreps ()
      Allocate (this%integrals(n_xyz, n_ip,n_ua), Stat=status)
      If (status .Ne. 0) Call error_handler ("dipoleg_allocate: allo&
           &cate of dipoleg_integrals failed")

      ual:Do i_ua = 1, n_ua
         Do i_xyz = 1,n_xyz
            i_ip = 1
            ir1l:Do i_ir = 1, n_ir
               dim_ir = dimensions_arr (i_ir)

               pa1l: Do i_pa = 1, symmetry_data_n_partners_proj  (i_ir)

                  call init_dipole_integral_type(this%integrals(i_xyz, i_ip,i_ua),dim_ir,usedoff )

                  i_ip = i_ip + 1
               end Do pa1l
            end Do ir1l
         end Do
      end Do ual
case (1)
 n_ip = symmetry_data_n_ip ()
 n_ir = symmetry_data_n_irreps ()

Allocate (this%integrals(n_xyz, n_ip,n_ua), Stat=status)
If (status .Ne. 0) Call error_handler ("dipoleg_allocate: allo&
     &cate of dipoleg_integrals failed")

Do i_ua = 1, n_ua
   Do i_xyz = 1,n_xyz
      i_ip = 1
      Do i_ir = 1, n_ir
         dim_ir = dimensions_arr (i_ir)
         !dim_ir =  symmetry_data_dimension (i_ir) ! contracted or uncontracted?

          Do i_pa = 1, symmetry_data_n_partners (i_ir)

                  call init_dipole_integral_type(this%integrals(i_xyz, i_ip,i_ua),dim_ir,usedoff )

                  i_ip = i_ip + 1
               end Do
            end Do
         end Do
      end Do
case default
Call error_handler ("init_hfc_integrals: &
     & forbiden case")
end select
end Subroutine init_hfc_integrals

Subroutine null_hfc_integrals(this)
Type (properties_integrals_type):: this
Integer (Kind=i4_kind) :: Status, i_ip, i_ua, i_xyz
Integer (Kind=i4_kind) :: dimensions(3)
!------------------------------------

dimensions = shape(this%integrals)

ual:Do i_ua = 1,dimensions(3)
   pa2l:Do i_ip = 1, dimensions(2)
         Do i_xyz = 1, dimensions(1)
            call null_dipole_integral_type(this%integrals(i_xyz,i_ip,i_ua))
         end Do
   end Do pa2l
end Do ual
if (associated(this%integrals)) then
deallocate (this%integrals, Stat=status)
If (status .Ne. 0) Call error_handler ("dipoleg_allocate: allo&
     &cate of dipoleg_integrals failed")
end if
end Subroutine null_hfc_integrals


Subroutine init_dipole_integral_type(this,n,use )
type(dipole_integral_type), intent(inout) :: this
Integer (Kind=i4_kind), intent(in) ::n, use
!------------end of interface----------------
Integer (Kind=i4_kind) :: status
select case (use)
case (1) !diagonal two component triangular matrix
   Allocate(this%diagonal(n*(n+1)/2), this%diagonal_imag(n*(n+1)/2),Stat=status )
   If (status .Ne. 0) Call error_handler ("hfc_alloc&
        &ate: allocate dipole_integrals%diagonal failed")
   !init
   this%use = DIPOLE_DIAGONAL
   this%diagonal = 0.0_r8_kind
   this%diagonal_imag = 0.0_r8_kind
   !nullify(this%offdiagonal, this%offdiagonal_imag)
case(2)! offdiagonal square matrix
   Allocate(this%offdiagonal(n,n), this%offdiagonal_imag(n,n),Stat=status )
   If (status .Ne. 0) Call error_handler ("hfc_alloc&
        &ate: allocate dipole_integrals%diagonal failed")
   !init
   this%use = DIPOLE_OFFDIAGONAL
   this%offdiagonal = 0.0_r8_kind
   this%offdiagonal_imag = 0.0_r8_kind
   !nullify(this%diagonal, this%diagonal_imag)
case (3) !one component triangular matrix
Allocate(this%diagonal(n*(n+1)/2), Stat=status )
   If (status .Ne. 0) Call error_handler ("hfc_alloc&
        &ate: allocate dipole_integrals%diagonal failed")
   !init
   this%use = DIPOLE_DIAGONAL
   this%diagonal = 0.0_r8_kind
   !nullify(this%offdiagonal, this%offdiagonal_imag)
case (4) !one component square matrix
   Allocate(this%offdiagonal(n,n), Stat=status )
   If (status .Ne. 0) Call error_handler ("hfc_alloc&
        &ate: allocate dipole_integrals%diagonal failed")
   !init
   this%use = DIPOLE_OFFDIAGONAL
   this%offdiagonal = 0.0_r8_kind
   !nullify(this%diagonal, this%diagonal_imag)
case default

Call error_handler ("hfc_alloc&
     &ate: forbidden case")
end select
end Subroutine init_dipole_integral_type

Subroutine null_dipole_integral_type(this )
type(dipole_integral_type), intent(inout) :: this
!------------end of interface----------------
Integer (Kind=i4_kind) :: status
if (allocated(this%diagonal)) then
   Deallocate (this%diagonal, Stat=status)
   If (status .Ne. 0) Call error_handler ("dipoleg_free:&
        & deallocate of dipoleg_integrals%diagonal failed")
   !nullify(this%diagonal)
end if
if (allocated(this%diagonal_imag)) then
   Deallocate (this%diagonal_imag, Stat=status)
   If (status .Ne. 0) Call error_handler ("dipoleg_free:&
        & deallocate of dipoleg_integrals%diagonal failed")
   !nullify(this%diagonal_imag)
end if
if (allocated(this%offdiagonal)) then
   Deallocate (this%offdiagonal, Stat=status)
   If (status .Ne. 0) Call error_handler ("dipoleg_free:&
        & deallocate of dipoleg_integrals%diagonal failed")
   !nullify(this%offdiagonal)
end if
if (allocated(this%offdiagonal_imag)) then
   Deallocate (this%offdiagonal_imag, Stat=status)
   If (status .Ne. 0) Call error_handler ("dipoleg_free:&
        & deallocate of dipoleg_integrals%diagonal failed")
   !nullify(this%offdiagonal_imag)
end if
this%use = dipole_unused
end Subroutine null_dipole_integral_type

Subroutine load_matrix (n, A, i_ir, filename)
  !  Purpose: reads matrix from file to 2D array
  !
  !
  !
  !
  !
  !
  !  Author: DG
  !  Date:   6/2003
  !
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  Use quadrupel_fname, Only: qfilename
  Use io, Only: read_buffer
  use iounitadmin_module
  Implicit none
  integer(IK), intent(in)         ::n    !dimension
  integer(IK), intent(in)         ::i_ir !dimension
  real(RK), intent(out) :: A(:,:)             !Matrix
  character(len=*) :: filename

Call read_buffer(qfilename(filename, i_ir, "dat") , A)
end Subroutine load_matrix

Subroutine load_vector(n, A, i_ir, filename)
  !  Purpose: reads vector from file to  array
  !
  !
  !
  !
  !
  !
  !  Author: DG
  !  Date:   6/2003
  !
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  Use quadrupel_fname, Only: qfilename
  Use io, Only: read_buffer
  Implicit none
  integer(IK), intent(in)         ::n    !dimension
  integer(IK), intent(in)         ::i_ir !dimension
  real(RK), intent(out) :: A(:)             !Matrix
  character(len=*) :: filename
  ! *** end of interface ***

  call read_buffer(qfilename(filename, i_ir, "dat"), A(:n))
End Subroutine load_vector

Subroutine hfc_calculate
  use spin_orbit_module, only: spin_orbit_polarized
  Implicit None
  Integer (Kind=i4_kind) :: n_ip, n_ir, i_ip1, n_orb, i_ir, n_pa, i_pa,i_orb
  Real (kind=r8_kind):: occ_num_current
  Logical  :: spin_polarized
  !-----------------------------------------

  ! CC (Kramers) Irrep Coupling
  !cc_coupling => symmetry_data_get_cccoupling()
  ! Is it non-trivial?:
  !spin_polarized = any(cc_coupling.NE.(/(i,i=1,size(cc_coupling))/))
  spin_polarized = spin_orbit_polarized
  DPRINT 'gtensor_calculate: spin_polarized=',spin_polarized

  if (options_spin_orbit) then
     n_ip = symmetry_data_n_ip_proj ()
     n_ir = symmetry_data_n_proj_irreps ()
     !
     i_ip1 = 0
     n_orb = symmetry_data_dimension_proj (1)
     !Searching for the first partner
     a_mx_el%eigenvalue = eigval(1)%m(1, 1)
     irrep: Do i_ir = 1, n_ir
        n_orb = symmetry_data_dimension_proj (i_ir)
        n_pa = symmetry_data_n_partners_proj (i_ir)
        partner1: Do i_pa = 1, n_pa
           i_ip1 = i_ip1 + 1
           orbitals_sog: Do i_orb = n_orb, 1, - 1 !all occ and empty orb`s
              occ_num_current = occ_num(i_ir)%m(i_orb, 1)
              If ((occ_num_current < (1+very_small) .And. &
                   & occ_num_current > (1-very_small)) .Or. (occ_num_current &
                   & < (2+very_small) .And. occ_num_current > (2-very_small))) Then
                 If (eigval(i_ir)%m(i_orb, 1) > a_mx_el%eigenvalue .or.&
                      abs(eigval(i_ir)%m(i_orb, 1)-a_mx_el%eigenvalue )< very_small) Then
                    a_mx_el%eigenvalue = &
                         & eigval(i_ir)%m(i_orb, 1)
                    a_mx_el%i_ip = i_ip1
                    a_mx_el%i_irr = i_ir
                    a_mx_el%i_orb = i_orb
                    a_mx_el%n_orb = n_orb
                    a_mx_el%i_pa = i_pa
                 End If
                 Cycle orbitals_sog
              End If
           End Do orbitals_sog
        End Do partner1
     End Do irrep
  else ! no spin orbit
     n_ip = symmetry_data_n_ip ()
     n_ir = symmetry_data_n_irreps ()
     !
     i_ip1 = 0
     n_orb = symmetry_data_dimension (1)
     !Searching for the first partner
     a_mx_el%eigenvalue = eigval(1)%m(1, 1)
     irr: Do i_ir = 1, n_ir
        n_orb = symmetry_data_dimension (i_ir)
        n_pa = symmetry_data_n_partners (i_ir)
        part: Do i_pa = 1, n_pa
           i_ip1 = i_ip1 + 1
           orbitals: Do i_orb = n_orb, 1, - 1 !all occ and empty orb`s
              occ_num_current = occ_num(i_ir)%m(i_orb, 1)
              If ((occ_num_current < (1+very_small) .And. &
                   & occ_num_current > (1-very_small)) .Or. (occ_num_current &
                   & < (2+very_small) .And. occ_num_current > (2-very_small))) Then
                 If (eigval(i_ir)%m(i_orb, 1) > a_mx_el%eigenvalue .or.&
                      abs(eigval(i_ir)%m(i_orb, 1)-a_mx_el%eigenvalue )< very_small) Then
                    a_mx_el%eigenvalue = &
                         & eigval(i_ir)%m(i_orb, 1)
                    a_mx_el%i_ip = i_ip1
                    a_mx_el%i_irr = i_ir
                    a_mx_el%i_orb = i_orb
                    a_mx_el%n_orb = n_orb
                    a_mx_el%i_pa = i_pa
                 End If
                 Cycle orbitals
              End If
           End Do orbitals
        End Do part
     End Do irr
  end if

        if (options_spin_orbit .and. options_kinematic_factors) then

              call reltrafo_hfc()

        end if

         if (options_relativistic .and. options_kinematic_factors &
              .and. (.not. options_spin_orbit) ) then
           call reltrafo_hfc_one_comp()
         end if

         if (options_spin_orbit) then
            SP:   If(spin_polarized) then
               call  calculate_hfc_dipole_so_ncsdft()
            else
               call calculate_hfc_dipole_so()
            end If SP
            !call  calculate_hfc_dipole_so_ncsdft()
            !call calculate_hfc_dipole_so()
        else !non-relat/scalar relat
            call calculate_hfc_dipole_one_comp()
        end if
      end Subroutine hfc_calculate


       Subroutine remap_to_full_store (in,out)
            Type (dipole_integral_type), Pointer:: in, out
            ! *** end of interface ***


            Integer (Kind=i4_kind) :: n
            Integer (Kind=i4_kind) :: i_bas12, i_bas1, i_bas2
#if _DPRINT
            real(r8_kind) :: diff
#endif

            DPRINT 'remap_to_full_store entered'
            !hfc => hfc_integrals (i_xyz, i_ip2, i_ip1, i_ua)

            integralstorage: Select Case (in%use)
            Case (DIPOLE_DIAGONAL) integralstorage

               n = size(out%offdiagonal,1)
               DPRINT 'remap_to_full_store: n=',n
ASSERT(n==size(out%offdiagonal,2))
ASSERT(n*(n+1)==2*size(in%diagonal))

               !triangular storage
               i_bas12 = 0
!
               bas1_diagonal_so: Do i_bas2 = 1, n
                  bas2_diagonal_so: Do i_bas1 = 1, i_bas2
                     i_bas12 = i_bas12 + 1
                     out%offdiagonal (i_bas1, i_bas2) = &
                    & in%diagonal(i_bas12)
                     out%offdiagonal (i_bas2, i_bas1) = &
                    & in%diagonal(i_bas12)
!
                     if (allocated(in%diagonal_imag)) then

                     out%offdiagonal_imag (i_bas1, i_bas2) = &
                    & in%diagonal_imag(i_bas12)
                     out%offdiagonal_imag (i_bas2, i_bas1) = - & !matrix is Hermitian
                    & in%diagonal_imag(i_bas12)

                  end if
                  End Do bas2_diagonal_so
               End Do bas1_diagonal_so
#if _DPRINT
               diff = maxval(abs(out%offdiagonal-transpose(out%offdiagonal)))
               DPRINT 'remap_to_full_store: hermitean.re?=',diff
               diff = maxval(abs(out%offdiagonal_imag+transpose(out%offdiagonal_imag)))
               DPRINT 'remap_to_full_store: hermitean.im?=',diff
               diff = maxval(abs(out%offdiagonal))
               DPRINT 'remap_to_full_store: maxval.im=',diff
               diff = maxval(abs(out%offdiagonal_imag))
               DPRINT 'remap_to_full_store: maxval.re=',diff
#endif
!
            Case (DIPOLE_OFFDIAGONAL) integralstorage
               DPRINT 'remap_to_full_store: case DIPOLE_OFFDIAGONAL'
            Case (dipole_unused) integralstorage
               DPRINT 'remap_to_full_store: case dipole_unused'
            Case Default integralstorage
               Call error_handler&
                    ("from remap_to_full_store: dipole_transitionmoment_f: forbidden case")
            End Select integralstorage
            DPRINT 'remap_to_full_store done'
         End Subroutine remap_to_full_store

         Subroutine remap_to_compact_store (hfc,hfc_out, n_orb)
      ! pice of code to remap gti_contracted to  the compact store
            Integer (Kind=i4_kind), Intent (In) :: n_orb
            Type (dipole_integral_type), Pointer:: hfc, hfc_out
            Integer (Kind=i4_kind) :: i_bas12, i_bas1, i_bas2
            DPRINT'remap_to_compact_store entered'

            integralstorage: Select Case (hfc%use)
            Case (DIPOLE_OFFDIAGONAL) integralstorage
         !triangular storage
               i_bas12 = 0
!
               bas1_diagonal_so: Do i_bas2 = 1, n_orb
                  bas2_diagonal_so: Do i_bas1 = 1, i_bas2
                     i_bas12 = i_bas12 + 1
                     hfc_out%diagonal (i_bas12) = hfc%offdiagonal(i_bas1, &
                    & i_bas2)
!
                     if (allocated(hfc%offdiagonal_imag)) then

                    hfc_out%diagonal_imag (i_bas12) = &
                    & hfc%offdiagonal_imag(i_bas1, i_bas2)

                  end if
                  End Do bas2_diagonal_so
               End Do bas1_diagonal_so
!
            Case (DIPOLE_DIAGONAL) integralstorage
         ! diagonal storage
            Case (dipole_unused) integralstorage
            Case Default integralstorage
               Call error_handler ("remap_to_compact_store: forbidden c&
              &ase")
            End Select integralstorage
            DPRINT'remap_to_compact_store done'
         End Subroutine remap_to_compact_store

!*************************************************************
Subroutine hfc_allocate()
Implicit None
!-----------------------------------

call hfc_init()

if(options_spin_orbit) then

if(options_kinematic_factors) then
! uncontracted dimensions:
call new(hfc_integrals_diag,        9, n_unique_atoms, DIPOLE_DIAGONAL,   UNCONTRACTED, 2)

! uncontracted dimensions:
call new(hfc_integrals_offdiag,     9, n_unique_atoms, DIPOLE_OFFDIAGONAL,UNCONTRACTED, 2)

! remapped from diagonal storage to offdiagonal:
call new(hfc_integrals_rem,         9, n_unique_atoms, DIPOLE_OFFDIAGONAL,UNCONTRACTED, 2)

! for kinematic factors multiplication etc.
call new(hfc_integrals_cont_rem,    9, n_unique_atoms, DIPOLE_OFFDIAGONAL,CONTRACTED,   2)
call new(hfc_integrals_cont_offdiag,9, n_unique_atoms, DIPOLE_OFFDIAGONAL,CONTRACTED,   2)
else
! contracted dimension:
call new(hfc_integrals_diag,        9, n_unique_atoms, DIPOLE_DIAGONAL,   CONTRACTED,   2)

! contracted dimension:
call new(hfc_integrals_offdiag,     9, n_unique_atoms, DIPOLE_OFFDIAGONAL,CONTRACTED,   2)
end if
else !not spin_orbit

if(options_kinematic_factors) then
! The initial triangular matrix is stored here:
call new(hfc_integrals_diag,        7, n_unique_atoms, 3,                 UNCONTRACTED, 1)

! real square matrix:
call new(hfc_integrals_rem,         7, n_unique_atoms, 4,                 UNCONTRACTED, 1)

! 4 is an index for the square one-component matrix:
call new(hfc_integrals_cont_rem,    7, n_unique_atoms, 4,                 CONTRACTED,   1)
else
! 3 is the index for  dipole diagonal one component:
call new(hfc_integrals_diag,        7, n_unique_atoms, 3,                 CONTRACTED,   1)
end if
end if

end Subroutine hfc_allocate
!*************************************************************

Subroutine hfc_free()
Implicit None
!-----------------------------------
if (options_spin_orbit) then
call delete(hfc_integrals_diag)
call delete(hfc_integrals_offdiag)
if(options_kinematic_factors) then
call delete(hfc_integrals_rem)
call delete(hfc_integrals_cont_rem)
call delete(hfc_integrals_cont_offdiag)
end if
else! not SO

if(options_kinematic_factors) then
call delete (hfc_integrals_rem)
call delete (hfc_integrals_cont_rem)
end if
call delete(hfc_integrals_diag)
end if

call hfc_done()
end Subroutine hfc_free


  !*************************************************************
     Subroutine calc_write_tmgt (hfc, n_orb1,ev_real,ev_imag, tm_real,tm_imag )
      ! arguments:
      ! i_ip1 - irrep-partner
      ! i_ip2=1 for diagonal matrix elements
      ! i_ip2=2 for offdiagonal
      ! n_orb1 - number of orbitals (the size of matrix)
      ! i_irb - ...
      ! Porpose: calculation  of
      ! transition moments (matrix elements) between  orbital and the same orbital (diagonal)
      ! and Kramers cogugation (matrix offdiagonal)
          !  Use orbitalprojection_module, Only: &
          ! & orbitalprojection_spor_ob, orbitalprojection_spor_ob_k
          !  Use unique_atom_module
          !  Use dimensions, Only: BasDimSpor, SymDimSpor
            Integer (Kind=i4_kind), Intent (In) :: n_orb1
            Real (Kind=r8_kind), Intent(Out) :: tm_real, tm_imag
            Real (Kind=r8_kind), Dimension (:), Pointer :: ev_real, ev_imag

            Type (dipole_integral_type), Pointer:: hfc

            Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12
            Real (Kind=r8_kind) :: coeff_real_a, &
           & coeff_imag_a, coef1, coef2

            Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind

!
!
!
            tm_real = 0.0_r8_kind
            tm_imag = 0.0_r8_kind
!
            ! mapping
            ! 1 = L1/r^3 ;  4 = (3Sx * x^2 - Sx * r^2)/r^5 + 3Sy *yz/r^5 + 3Sx *xz/r^5
            ! 2 = L2/r^3 ;  5 = (3Sy * r^2 - Sy * r^2)/r^5 + 3Sx *yx/r^5 + 3Sy *yz/r^5
            ! 3 = L3/r^3 ;  6 = (3Sz * z^2 - Sz * r^2)/r^5 + 3Sx *xz/r^5 + 3Sy *yz/r^5
            ! 7,8,9 = iso term




                integralstorage_s: Select Case ( hfc%use)

                Case (1) integralstorage_s
!
            !triangular storage
                  i_bas12 = 0


                  bas1_diagonal_so_d: Do i_bas2 = 1, n_orb1
!
                     bas2_diagonal_so_d: Do i_bas1 = 1, i_bas2 - 1
!
                        i_bas12 = i_bas12 + 1
!
                        coef1 = (ev_real(i_bas2)*&
                       & ev_real(i_bas1)+&
                       & ev_imag(i_bas2)*&
                       & ev_imag(i_bas1)+&
                       & ev_real(i_bas1)*&
                       & ev_real(i_bas2)+&
                       & ev_imag(i_bas1)*&
                       & ev_imag(i_bas2))
!
                        coef2 = (ev_real(i_bas2)*&
                       & ev_imag(i_bas1)-&
                       & ev_imag(i_bas2)*&
                       & ev_real(i_bas1)-&
                       & ev_real(i_bas1)*&
                       & ev_imag(i_bas2)+&
                       & ev_imag(i_bas1)*&
                       & ev_real(i_bas2))
!
!
                        tm_real = tm_real &
                       &  +  hfc%diagonal(i_bas12) * coef1 + &
                       &  hfc%diagonal_imag(i_bas12) * coef2
!
!
                     End Do bas2_diagonal_so_d
!
                     i_bas12 = i_bas12 + 1
!
!
                     coef1 = (ev_real(i_bas2)*&
                    & ev_real(i_bas1)+&
                    & ev_imag(i_bas2)*&
                    & ev_imag(i_bas1))
!
                     coef2 = (ev_real(i_bas2)*&
                    & ev_imag(i_bas1)-&
                    & ev_imag(i_bas2)*&
                    & ev_real(i_bas1))
!
!
                     tm_real = tm_real &
                    &  +  hfc%diagonal(i_bas12) * coef1 + &
                    &  hfc%diagonal_imag(i_bas12) * coef2
!
                  End Do bas1_diagonal_so_d
!
               Case (2) integralstorage_s
            ! full storage

!
                  Do i_bas2 = 1, n_orb1
!
                     Do i_bas1 = 1, n_orb1

                           coeff_real_a = (ev_real(i_bas2)*&
                          & ev_real(i_bas1)-&
                          & ev_imag(i_bas2)*&
                          & ev_imag(i_bas1))
                           coeff_imag_a = (ev_real(i_bas2)*&
                          & ev_imag(i_bas1)+&
                          & ev_imag(i_bas2)*&
                          & ev_real(i_bas1))
!
!
                           tm_real = &
                          & tm_real + &
                          &  hfc%offdiagonal(i_bas1, i_bas2) * &
                          & coeff_real_a
                           tm_imag  = &
                          & tm_imag  + &
                          &  hfc%offdiagonal(i_bas1, i_bas2) * &
                          & coeff_imag_a
!
                           tm_real  = &
                          & tm_real - &
                          & hfc%offdiagonal_imag(i_bas1, i_bas2) * &
                          & coeff_imag_a
!
                           tm_imag  = &
                          & tm_imag  + &
                          &  hfc%offdiagonal_imag(i_bas1, i_bas2) * &
                          & coeff_real_a

!
                     End Do
                  End Do
!
               Case (0) integralstorage_s
                  tm_real = 0.0_r8_kind
                  tm_imag = 0.0_r8_kind
               Case Default integralstorage_s
                  Call error_handler&
                       ("calc_write_tmgt: sxyz:  dipole_transitionmoment_f: forbidden case")
               End Select integralstorage_s
           ! end Do sxyz
          end Subroutine calc_write_tmgt

       subroutine  calculate_hfc_dipole_so
         !  Purpose: calculates a-tensor components in SO case
         !
         !** End of interface *****************************************


           Implicit None
         Integer (Kind=i4_kind) :: i_xyz
         Integer (Kind=i4_kind) :: i_ua
          Real (Kind=r8_kind), Pointer :: eigenvector_real1 (:, :), &
        & eigenvector_imag1 (:, :)
         Real (Kind=r8_kind), Dimension (:), Pointer :: &
        & eigenvector1_real, eigenvector1_imag

         Real (Kind=r8_kind) ::tm_r, tm_i
         Type (dipole_integral_type), Pointer:: hfcint

               eigenvector_real1 => eigvec_real(a_mx_el%i_irr)%m(:, :)
               eigenvector_imag1 => eigvec_imag(a_mx_el%i_irr)%m(:, :)
               !
               eigenvector1_real => eigenvector_real1 (:, a_mx_el%i_orb)
               eigenvector1_imag => eigenvector_imag1 (:, a_mx_el%i_orb)

               ualoop: Do i_ua = 1 ,n_unique_atoms
                 sxyz: Do i_xyz = 1, 9
                  okf: If (  options_kinematic_factors) Then

                     hfcint => hfc_integrals_diag%integrals (i_xyz, a_mx_el%i_ip , i_ua)

                  else
                     hfcint => hfc_integrals_diag%integrals (i_xyz, a_mx_el%i_ip , i_ua)

                  end If okf
               Call calc_write_tmgt (hfcint, a_mx_el%n_orb, eigenvector1_real, eigenvector1_imag, tm_r, tm_i )
!
               select case(i_xyz)
                  case (4:6)
               a_mx_el%Sigma_Real_diag(i_xyz-3) = tm_r
               a_mx_el%Sigma_Im_diag(i_xyz-3) = tm_i
               case (1:3)
               a_mx_el%L_Real_diag(i_xyz) = tm_r
               a_mx_el%L_Im_diag(i_xyz) = tm_i
               case (7:9)
               iso_mx_el%sigma_real_diag(i_xyz-6) = tm_r

               !
!
            end select

            If (  options_kinematic_factors) Then

               hfcint => hfc_integrals_cont_offdiag%integrals (i_xyz, a_mx_el%i_ip , i_ua)

                  else

                     hfcint => hfc_integrals_offdiag%integrals (i_xyz, a_mx_el%i_ip , i_ua)
                  end If
               Call calc_write_tmgt (hfcint, a_mx_el%n_orb, eigenvector1_real, eigenvector1_imag, tm_r, tm_i )
               select case(i_xyz)
                  case (4:6)
               a_mx_el%Sigma_Real_offdiag(i_xyz-3) = tm_r
               a_mx_el%Sigma_Im_offdiag(i_xyz-3) = tm_i
               case (1:3)
               a_mx_el%L_Real_offdiag(i_xyz) = tm_r
               a_mx_el%L_Im_offdiag(i_xyz) = tm_i

            end select
            end Do sxyz
               Call calc_A_tensor ()
            end Do ualoop


           Contains


          Subroutine calc_A_tensor ()
!
            Real (Kind=r8_kind), Dimension (3, 3) :: g_Sigma, g_L, g
            Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      !--------------------------------------------------------------------------
            Write (output_unit,*) 'Kinematic factors is ', &
           & options_kinematic_factors
            Write (output_unit,*) 'Unique atom : ', &
                 & i_ua
            Write (output_unit, Fmt='(4(A,I3,2x),A,F14.7,A)') "irrep: ",&
           &  a_mx_el%i_irr, "partner: ", a_mx_el%i_pa, "irrep-partner:&
           & ", a_mx_el%i_ip, "orb: ", a_mx_el%i_orb, " E  =", &
           & a_mx_el%eigenvalue * convert1, " eV"


            Write (output_unit,&
                 Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
                 & a_mx_el%Sigma_Real_diag , &
           & a_mx_el%Sigma_Im_diag ,&
           & a_mx_el%L_Real_diag , &
           & a_mx_el%L_Im_diag
            Write (output_unit,&
           & Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
           & a_mx_el%Sigma_Real_offdiag , &
           & a_mx_el%Sigma_Im_offdiag ,&
           & a_mx_el%L_Real_offdiag , &
           & a_mx_el%L_Im_offdiag

            g_Sigma (:, 1) = a_mx_el%Sigma_Real_offdiag(:)
            g_Sigma (:, 2) = -a_mx_el%Sigma_Im_offdiag(:)
            g_Sigma (:, 3) = a_mx_el%Sigma_Real_diag(:)
           ! g_Sigma (:, 3) = - g_mx_el%Sigma_Real_diag(:)
!
            g_L (:, 1) = a_mx_el%L_Real_offdiag(:)
            g_L (:, 2) = -a_mx_el%L_Im_offdiag(:)
            g_L (:, 3) = a_mx_el%L_Real_diag(:)

           ! g_L (:, 3) = - g_mx_el%L_Real_diag(:)
            call print_iso()
            g = 0.0_r8_kind !DG
            g (:, :) = (g_Sigma(:, :)+g_L(:, :)) !coeff!
            Write (output_unit,*) "A-tensor" !
            Call print_A_tensors (g)
            g = 0.0_r8_kind
            g (:, :) = g_Sigma (:, :)!coeff
            Write (output_unit,*) "dipole contribution"
            Call print_A_tensors (g)
            g = 0.0_r8_kind
!
            g (:, :) = g_L (:, :) !coeff
            Write (output_unit,*) "dipole SO contribution"
            Call print_A_tensors (g)

          End Subroutine calc_A_tensor
          Subroutine print_iso
          Real (Kind=r8_kind) :: coeff2, coeff_au2, tmp
          !----------------------------------------
          if (options_finite_nucleus) then
            coeff_au2 = coeff_au
          else
            coeff_au2 = -coeff_au * 1/(4*PI)
          end if
          coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
               & unique_atoms(i_ua)%nuclear_spin
          Write (output_unit,*)&
                 & 'Aiso_x   a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2,&
                    iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2 * K_convert_to_Gauss
               Write (output_unit,*)&
                 & 'Aiso_y   a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') iso_mx_el%Sigma_real_diag(2) * coeff_au2 * coeff2,&
                    iso_mx_el%Sigma_real_diag(2) * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    iso_mx_el%Sigma_real_diag(2) * coeff_au2 * coeff2 * K_convert_to_Gauss
               Write (output_unit,*)&
                 & 'Aiso_z   a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') iso_mx_el%Sigma_real_diag(3) * coeff_au2 * coeff2,&
                    iso_mx_el%Sigma_real_diag(3) * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    iso_mx_el%Sigma_real_diag(3) * coeff_au2 * coeff2 * K_convert_to_Gauss

               tmp = sqrt(iso_mx_el%Sigma_real_diag(1)**2 +iso_mx_el%Sigma_real_diag(2)**2 + iso_mx_el%Sigma_real_diag(3)**2 )

               Write (output_unit,*)&
                 & 'Aiso_module a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') tmp * coeff_au2 * coeff2,&
                    tmp * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    tmp * coeff_au2 * coeff2 * K_convert_to_Gauss

          end Subroutine print_iso
         Subroutine print_A_tensors (g)
            Real (Kind=r8_kind), Dimension (3, 3), Intent (In) :: g

      !-------------------------------------------------
            Real (Kind=r8_kind), Dimension (3, 3) :: G_tensor
            Real (Kind=r8_kind), Parameter :: &
           & diagonal_request_precision = 0.0000001
            Real (Kind=r8_kind), Dimension (3) :: diag
            Real (Kind=r8_kind) :: W (3), WORK (9), coeff2, coeff_au1, coeff_au2
            Logical diagonalisation_request
            Integer i, j, m, info
            Real (Kind=r8_kind) sm, buf
            diagonalisation_request = .False.
            coeff_au1 = coeff_au * 3/(8*PI)
            coeff_au2 = -coeff_au * 1/(4*PI)
            coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
                 & unique_atoms(i_ua)%nuclear_spin


            Write (output_unit,*) "pre a"
            Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') g (1, :), &
           & g (2, :), g (3, :)
!
            Do i = 1, 3
               Do j = 1, 3
                  sm = 0.0
                  Do m = 1, 3
                     buf = g (i, m) * g (j, m)
                     sm = sm + buf
                  End Do
                  G_tensor (i, j) = sm
                  If (i /= j .And. Abs(sm) > &
                 & diagonal_request_precision) diagonalisation_request &
                 & = .True.
               End Do
            End Do
            Write (output_unit,*) "Big A"
            Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') G_tensor &
           & (1, :), G_tensor (2, :), G_tensor (3, :)
!
            If (diagonalisation_request) Then
               Write (output_unit,*) "A tensor is not diagonal"
               Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') &
              & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
               Call DSYEV ('V', 'L', 3, G_tensor, 3, W, WORK, 9, info)
               If (info /= 0) Then
                  Write (output_unit,*) "Could not diagonalize G-tensor&
                 &"
                  Return
               End If
               Write (output_unit,*) "Diagonal components of A"
               Write (output_unit, Fmt='(3F14.7)') W
               Write (output_unit,*) "Directions of the main axes of g &
              &tensor (eigen vectors)"
               Write (output_unit, Fmt='(3F14.7/3F14.7/3F14.7)') &
              & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
               Write (output_unit,*) "Angles in degrees"
               Write (output_unit, Fmt='(3F14.7)') (Acos(G_tensor(1, &
              & 1))) * 180 / PI, (Acos(G_tensor(2, 2))) * 180 / PI, &
              & (Acos(G_tensor(3, 3))) * 180 / PI
               Write (output_unit,*) "A-tensor"
              ! Write (output_unit, Fmt='(3F14.7)') Sqrt (Abs(W(:)))
               diag(:)= Sqrt (Abs(W(:)))

             !  Return
            Else
               Do i = 1, 3
                  diag (i) = Sqrt (Abs(G_tensor(i, i)))
               End Do
            end If
               diag(:) = diag (:)-sum(diag(:)/3) ! Traceless, sign of A should be taken from g-tensor
               Write (output_unit,*) "au         xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2
               Write (output_unit,*) "MHz        xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2 *  K_convert_to_MHz
               Write (output_unit,*) "Gauss      xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2 *  K_convert_to_Gauss

          !  End If
          End Subroutine print_A_tensors

        end subroutine calculate_hfc_dipole_so

       subroutine  calculate_hfc_dipole_so_ncsdft
         !  Purpose: calculates a-tensor components in SO NCSDFT case, linear molecules
         !
         !** End of interface *****************************************


         Implicit None
         Integer (Kind=i4_kind) :: i_ir, n_ir, &
              & i_orb, n_orb, i_xyz
         Integer (Kind=i4_kind) :: i_ua
         Real (Kind=r8_kind), Pointer :: eigenvector_real1 (:, :), &
              & eigenvector_imag1 (:, :)
         Real (Kind=r8_kind), Dimension (:), Pointer :: &
              & eigenvector1_real, eigenvector1_imag

         Real (Kind=r8_kind) ::tm_r, tm_i
         Type (dipole_integral_type), Pointer:: hfcint

         n_ir = symmetry_data_n_proj_irreps ()

         ualoop_sp: Do i_ua = 1 ,n_unique_atoms


            a_mx_el%Sigma_Real_diag = 0.0_r8_kind
            a_mx_el%Sigma_Im_diag = 0.0_r8_kind
            iso_mx_el%sigma_real_diag = 0.0_r8_kind
            iso_mx_el%sigma_Im_diag = 0.0_r8_kind

         irr_sp: do i_ir = 1, n_ir !only for one-dimensional irrep, every irrep has only one parner
       !     irr_sp: do i_ir = a_mx_el%i_irr, a_mx_el%i_irr
            n_orb = symmetry_data_dimension_proj(i_ir)
            orb_sp: do i_orb = 1, n_orb
            ! orb_sp: do i_orb = a_mx_el%i_orb, a_mx_el%i_orb
               if ( occ_num(i_ir)%m(i_orb, 1) < very_small) cycle

               eigenvector_real1 => eigvec_real(i_ir)%m(:, :)
               eigenvector_imag1 => eigvec_imag(i_ir)%m(:, :)
               !
               eigenvector1_real => eigenvector_real1 (:, i_orb)
               eigenvector1_imag => eigenvector_imag1 (:, i_orb)



               sxyz: Do i_xyz = 1, 9


                     okf: If (  options_kinematic_factors) Then

                        hfcint => hfc_integrals_diag%integrals (i_xyz, i_ir , i_ua) !One-dimensional irreps only
                     else
                        hfcint => hfc_integrals_diag%integrals (i_xyz, i_ir , i_ua)

                     end If okf
                     Call calc_write_tmgt (hfcint, n_orb, eigenvector1_real, eigenvector1_imag, tm_r, tm_i )
                     !

                     select case(i_xyz)
                     case (4:6)
                        a_mx_el%Sigma_Real_diag(i_xyz-3) = a_mx_el%Sigma_Real_diag(i_xyz-3) + tm_r
                        a_mx_el%Sigma_Im_diag(i_xyz-3) = a_mx_el%Sigma_Im_diag(i_xyz-3) + tm_i
                     case (1:3)
                        a_mx_el%L_Real_diag(i_xyz) = a_mx_el%L_Real_diag(i_xyz)+ tm_r
                        a_mx_el%L_Im_diag(i_xyz) = a_mx_el%L_Im_diag(i_xyz) + tm_i
                     case(7:9)
                        iso_mx_el%sigma_real_diag(i_xyz - 6) = iso_mx_el%sigma_real_diag(i_xyz - 6) + tm_r
                        iso_mx_el%Sigma_Im_diag(i_xyz - 6) = iso_mx_el%Sigma_Im_diag(i_xyz - 6) + tm_i

                     end select

                  end Do sxyz


               end do orb_sp
            end do irr_sp
              Call calc_A_tensor ()
         end Do ualoop_sp

       Contains


          Subroutine calc_A_tensor ()
            ! FIXME: who needs it here???
            Real (Kind=r8_kind), Dimension (3, 3) :: g
            Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      !--------------------------------------------------------------------------
            Write (output_unit,*) 'Kinematic factors is ', &
           & options_kinematic_factors
            Write (output_unit,*) 'Unique atom : ', &
                 & i_ua
            Write (output_unit, Fmt='(4(A,I3,2x),A,F14.7,A)') "irrep: ",&
           &  a_mx_el%i_irr, "orb: ", a_mx_el%i_orb


            Write (output_unit,&
                 Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
                 & a_mx_el%Sigma_Real_diag , &
                 & a_mx_el%Sigma_Im_diag ,&
                 & a_mx_el%L_Real_diag , &
                 & a_mx_el%L_Im_diag



            call print_iso()

            Write (output_unit,*) "dipole SO contribution"
            ! FIXME: who sets this g???
            Call print_A_tensors (g)

          End Subroutine calc_A_tensor

          Subroutine print_iso
          Real (Kind=r8_kind) :: coeff2, coeff_au2, tmp
          !----------------------------------------
          if (options_finite_nucleus) then
            coeff_au2 = coeff_au
          else
            coeff_au2 = -coeff_au * 1/(4*PI)
          end if
          coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
               & unique_atoms(i_ua)%nuclear_spin
          Write (output_unit,*)&
                 & 'Aiso_x   '
               Write (output_unit, Fmt='(F18.7)') iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2

                Write (output_unit,*)&
                 & 'Aiso_y '
               Write (output_unit, Fmt='(F18.7)') iso_mx_el%Sigma_real_diag(2) * coeff_au2 * coeff2

               Write (output_unit,*)&
                 & 'Aiso_z '
               Write (output_unit, Fmt='(F18.7)') iso_mx_el%Sigma_real_diag(3) * coeff_au2 * coeff2


               tmp = sqrt(iso_mx_el%Sigma_real_diag(1)**2 +iso_mx_el%Sigma_real_diag(2)**2 + iso_mx_el%Sigma_real_diag(3)**2 )

               Write (output_unit,*)&
                 & 'Aiso_module a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') tmp * coeff_au2 * coeff2,&
                    tmp * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    tmp * coeff_au2 * coeff2 * K_convert_to_Gauss
          end Subroutine print_iso
         Subroutine print_A_tensors (g)
            Real (Kind=r8_kind), Dimension (3, 3), Intent (In) :: g

      !-------------------------------------------------
            Real (Kind=r8_kind), Parameter :: &
           & diagonal_request_precision = 0.0000001
            Real (Kind=r8_kind), Dimension (3) :: diag
            Real (Kind=r8_kind) :: coeff2, coeff_au1, coeff_au2
            Logical diagonalisation_request

            diagonalisation_request = .False.
            coeff_au1 = coeff_au * 3/(8*PI)
            coeff_au2 = -coeff_au * 1/(4*PI)
            coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
                 & unique_atoms(i_ua)%nuclear_spin

            diag(:) = a_mx_el%Sigma_real_diag(:) + a_mx_el%L_real_diag(:)

               Write (output_unit,*) "diagonal component au "
               Write (output_unit, Fmt='(3F18.7)') diag(:) * coeff_au1 * coeff2
               Write (output_unit,*)  "diagonal component au MHz"
               Write (output_unit, Fmt='(3F18.7)') diag(:) * coeff_au1 * coeff2 *  K_convert_to_MHz
               Write (output_unit,*) "diagonal component Gauss"
               Write (output_unit, Fmt='(3F18.7)') diag(:) * coeff_au1 * coeff2 *  K_convert_to_Gauss

          !  End If
          End Subroutine print_A_tensors

        end subroutine calculate_hfc_dipole_so_ncsdft


        subroutine  calculate_hfc_dipole_one_comp
         !  Purpose: calculates a-tensor components for one-component wave function
         !
         !** End of interface *****************************************
           Implicit None
         Integer (Kind=i4_kind) :: i_ir, n_ir, i_pa, n_pa, i_ip, n_ip, &
        & i_orb, n_orb, i_bas1, i_bas2, i_bas12, i_xyz
         Integer (Kind=i4_kind) :: i_ua, n_spin, i_spin
          Real (Kind=r8_kind), Pointer :: eigenvector (:,:)
         Real (Kind=r8_kind) ::tm_r, coeff
         Type (dipole_integral_type), Pointer:: hfcint
         Real (Kind=r8_kind), Dimension (3, 3) :: g
         Real (Kind=r8_kind), Dimension (3, 3,2) :: g_total
         Real (Kind=r8_kind), Dimension (2) :: iso_total

         !---------------Executable code--------------------

         iso_mx_el%sigma_real_diag = 0.0_r8_kind

          n_spin = symmetry_data_n_spin ()
          n_ip = symmetry_data_n_ip ()    !one-component case only
          n_ir = symmetry_data_n_irreps ()

          ualoop: Do i_ua = 1 ,n_unique_atoms
             iso_mx_el%sigma_real_diag = 0.0_r8_kind
             g = 0.0_r8_kind
             g_total = 0.0_r8_kind
             iso_total = 0.0_r8_kind
             sxyz: Do i_xyz = 1, 7

            i_ip = 1
            irnr:    Do i_ir = 1, n_ir


             n_orb = symmetry_data_dimension (i_ir)
             n_pa = symmetry_data_n_partners (i_ir)

             partner: Do i_pa = 1, n_pa

                spin: Do i_spin = 1, n_spin
                eigenvector => eigvec(i_ir)%m(:, :, i_spin)  !!Correct here DG

                      okf: If (  options_kinematic_factors) Then

                         !  hfcint => hfc_integrals_cont_rem%integrals (i_xyz, 1, a_mx_el%i_ip , i_ua)
                       !  hfcint => hfc_integrals_diag%integrals (i_xyz, a_mx_el%i_ip , i_ua)
                         hfcint => hfc_integrals_diag%integrals (i_xyz, i_ip, i_ua)
                      else

                         hfcint => hfc_integrals_diag%integrals (i_xyz, i_ip, i_ua)

                      end If okf
                   !   Call calc_write_tmgt

                      orbitals: Do i_orb = 1, n_orb

                           !Triangular storage

                    needs_sum : if ( abs( occ_num(i_ir)%m(i_orb, i_spin)-1.0_r8_kind) < very_small) then
                        ! search for the single occupied orbitals:

                           tm_r = 0.0_r8_kind
                           i_bas12 = 0

                        Do i_bas1 = 1, n_orb

                           Do i_bas2 = 1, i_bas1 - 1

                              i_bas12 = i_bas12 + 1

                              coeff = eigenvector (i_bas1, i_orb) * &
                                   & eigenvector (i_bas2, i_orb) *  2.0_r8_kind

                              tm_r = tm_r + hfcint%diagonal(i_bas12) * coeff


                           End Do

                           i_bas12 = i_bas12 + 1

                           coeff = eigenvector (i_bas1, i_orb) * &
                                & eigenvector (i_bas2, i_orb)

                           tm_r = tm_r +  hfcint%diagonal(i_bas12)* coeff

                        End Do

                      !
                      select case(i_xyz)

                      case (4)
                         g_total(1, 2,i_spin) = g_total(1, 2, i_spin) +  tm_r
                         g_total(2, 1, i_spin) = g_total(2, 1, i_spin) +  tm_r

                      case (5)
                         g_total(1, 3, i_spin) = g_total(1, 3, i_spin) +  tm_r
                         g_total(3, 1, i_spin) = g_total(3, 1, i_spin) +  tm_r

                      case (6)
                         g_total(2, 3, i_spin) = g_total(2, 3, i_spin) + tm_r
                         g_total(3, 2, i_spin) = g_total(3, 2, i_spin) + tm_r

                      case (1:3)
                         g_total(i_xyz, i_xyz, i_spin) = g_total(i_xyz, i_xyz, i_spin) + tm_r

                      case (7)
                         iso_total(i_spin) =  iso_total(i_spin)  +  tm_r

                      end select

                   end if needs_sum

                end Do orbitals
             end Do spin
                   i_ip = i_ip + 1
                end Do partner
                end Do irnr
             end Do sxyz
             g(:,:) = g_total(:,:,1) - g_total(:,:,2)
             iso_mx_el%sigma_real_diag(1) = iso_total(1)-iso_total(2)
             call calc_A_tensor
          end Do ualoop

contains

          Subroutine calc_A_tensor ()
!

            Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      !--------------------------------------------------------------------------
            Write (output_unit,*) 'Kinematic factors is ', &
           & options_kinematic_factors
            Write (output_unit,*) 'Unique atom : ', &
                 & i_ua
            Write (output_unit, Fmt='(4(A,I3,2x),A,F14.7,A)') "irrep: ",&
           &  a_mx_el%i_irr, "partner: ", a_mx_el%i_pa, "irrep-partner:&
           & ", a_mx_el%i_ip, "orb: ", a_mx_el%i_orb, " E  =", &
           & a_mx_el%eigenvalue * convert1, " eV"


!!$            Write (output_unit,&
!!$                 Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
!!$                 & a_mx_el%Sigma_Real_diag , &
!!$           & a_mx_el%Sigma_Im_diag ,&
!!$           & a_mx_el%L_Real_diag , &
!!$           & a_mx_el%L_Im_diag
!!$            Write (output_unit,&
!!$           & Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
!!$           & a_mx_el%Sigma_Real_offdiag , &
!!$           & a_mx_el%Sigma_Im_offdiag ,&
!!$           & a_mx_el%L_Real_offdiag , &
!!$           & a_mx_el%L_Im_offdiag
!!$
!!$            g (:, 1) = a_mx_el%Sigma_Real_offdiag(:)
!!$            g (:, 2) = -a_mx_el%Sigma_Im_offdiag(:)
!!$            g (:, 3) = a_mx_el%Sigma_Real_diag(:)

            call print_iso()
            !g = 0.0_r8_kind !DG
            !g (:, :) = g_Sigma(:, :) !coeff!
            Write (output_unit,*) "A-tensor" !
            Call print_A_tensors (g)

          End Subroutine calc_A_tensor
          Subroutine print_iso
          Real (Kind=r8_kind) :: coeff2, coeff_au2
          !----------------------------------------
          if (options_finite_nucleus) then
            coeff_au2 = coeff_au
          else
            coeff_au2 = -coeff_au * 1/(4*PI)
          end if
          coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
               & unique_atoms(i_ua)%nuclear_spin
          Write (output_unit,*)&
                 & 'Aiso   a.u.                MHz               Gauss '
               Write (output_unit, Fmt='(3F18.7)') iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2,&
                    iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2 * K_convert_to_MHz ,&
                    iso_mx_el%Sigma_real_diag(1) * coeff_au2 * coeff2 * K_convert_to_Gauss
             end Subroutine print_iso

             Subroutine print_A_tensors (g)
            Real (Kind=r8_kind), Dimension (3, 3), Intent (In) :: g

      !-------------------------------------------------
            Real (Kind=r8_kind), Dimension (3, 3) :: G_tensor
            Real (Kind=r8_kind), Parameter :: &
           & diagonal_request_precision = 0.0000001
            Real (Kind=r8_kind), Dimension (3) :: diag
            Real (Kind=r8_kind) :: W (3), WORK (9), coeff2, coeff_au1, coeff_au2
            Logical diagonalisation_request
            Integer i, j, m, info
            Real (Kind=r8_kind) sm, buf
            diagonalisation_request = .False.
            coeff_au1 = coeff_au * 3/(8*PI)
            coeff_au2 = -coeff_au * 1/(4*PI)
            coeff2 = unique_atoms(i_ua)%nuclear_magnetic_moment * 2 / &
                 & unique_atoms(i_ua)%nuclear_spin


            Write (output_unit,*) "pre a"
            Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') g (1, :), &
           & g (2, :), g (3, :)
!
            Do i = 1, 3
               Do j = 1, 3
                  sm = 0.0
                  Do m = 1, 3
                     buf = g (i, m) * g (j, m)
                     sm = sm + buf
                  End Do
                  G_tensor (i, j) = sm
                  If (i /= j .And. Abs(sm) > &
                 & diagonal_request_precision) diagonalisation_request &
                 & = .True.
               End Do
            End Do
            Write (output_unit,*) "Big A"
            Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') G_tensor &
           & (1, :), G_tensor (2, :), G_tensor (3, :)
!
            If (diagonalisation_request) Then
               Write (output_unit,*) "A tensor is not diagonal"
               Write (output_unit, Fmt='(3F18.7/3F18.7/3F18.7)') &
              & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
               Call DSYEV ('V', 'L', 3, G_tensor, 3, W, WORK, 9, info)
               If (info /= 0) Then
                  Write (output_unit,*) "Could not diagonalize G-tensor&
                 &"
                  Return
               End If
               Write (output_unit,*) "Diagonal components of A"
               Write (output_unit, Fmt='(3F14.7)') W
               Write (output_unit,*) "Directions of the main axes of g &
              &tensor (eigen vectors)"
               Write (output_unit, Fmt='(3F14.7/3F14.7/3F14.7)') &
              & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
               Write (output_unit,*) "Angles in degrees"
               Write (output_unit, Fmt='(3F14.7)') (Acos(G_tensor(1, &
              & 1))) * 180 / PI, (Acos(G_tensor(2, 2))) * 180 / PI, &
              & (Acos(G_tensor(3, 3))) * 180 / PI
               Write (output_unit,*) "A-tensor"
              ! Write (output_unit, Fmt='(3F14.7)') Sqrt (Abs(W(:)))
               diag(:)= Sqrt (Abs(W(:)))

             !  Return
            Else
               Do i = 1, 3
                  diag (i) = Sqrt (Abs(G_tensor(i, i)))
               End Do
            end If
               diag(:) = diag (:)-sum(diag(:)/3) ! Traceless, sign of A should be taken from g-tensor
               Write (output_unit,*) "au         xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2
               Write (output_unit,*) "MHz        xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2 *  K_convert_to_MHz
               Write (output_unit,*) "Gauss      xx                 yy               zz"
               Write (output_unit, Fmt='(3F18.7)') diag (:) * coeff_au1 * coeff2 *  K_convert_to_Gauss

          !  End If
             End Subroutine print_A_tensors
           end subroutine calculate_hfc_dipole_one_comp


         End Module hfc_module

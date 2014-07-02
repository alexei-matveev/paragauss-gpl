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
Module gtensor_module
  !
  !
  !  Purpose: calculates g-tensor components
  !
  !
  !
  !
  !
  !
  !  Author: DG
  !  Date:   1/2003
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
# include "def.h"
 !
  Use type_module
  Use symmetry_data_module
  Use occupation_module
  Use eigen_data_module
  Use options_module, Only: options_spin_orbit, &
       & options_kinematic_factors
  Use iounitadmin_module
  Use dipole_module
  use hfc_module
  Implicit None
  Save ! save all variables defined in this module
  Private ! by default, all names are private

  type, private :: spectrum
     integer(i4_kind)          :: n_occ
     real(r8_kind), pointer    :: energy(:) ! (n_occ)
     integer(i4_kind), pointer :: i_orb(:)  ! (n_occ), real index of level
  end type spectrum

  !------------ Declaration of variables --------------------------
!!$      Type (dipole_integral_type), Allocatable, Target, Public :: &
!!$     & gti_contracted (:, :, :)
  !    (7,2,2) - Kramers type blocks
  Type (dipole_integral_type), Allocatable, Target, Public :: &
       & dipoleg_integrals (:, :, :)
  ! (3, n_ip, n_ip)
  Type (properties_integrals_type), Public :: gten_integrals_diag,  gten_integrals_offdiag,&
       gten_integrals_rem,  gten_integrals_cont_rem,gten_integrals_cont_offdiag

!!$  Integer (Kind=i4_kind), Public :: n_exclude_atom_orbitals,n_orbitals, dkh_level
!!$  Real (Kind=r8_kind), Public :: mix_coeff
!!$  Logical, Public :: read_from_list
!!$  Logical, Public :: exclude_atom_orbitals, include_key
!!$  Namelist / gtensor_orbitals /&
!!$       & read_from_list,&
!!$       & n_orbitals, &
!!$       & exclude_atom_orbitals,&
!!$       & include_key,&
!!$       & n_exclude_atom_orbitals, &
!!$       & mix_coeff,&
!!$       & dkh_level ! DKH level for magnetic operator
!!$  Integer (Kind=i4_kind), Allocatable :: input_irreps (:), &
!!$       & input_partners (:), input_orbitals (:)
!!$  Integer (Kind=i4_kind), Allocatable, Public :: &
!!$       & input_exclude_unique_atoms_1 (:), input_exclude_atom_l_1 (:), &
!!$       & input_exclude_atom_n_1 (:), input_exclude_unique_atoms_2 (:), &
!!$       & input_exclude_atom_l_2 (:), input_exclude_atom_n_2 (:)
  Real (Kind=r8_kind), Parameter :: ge = 2.00231929
  !value from Jone E. Harriman "Theoretical Foundations of Electron Spin Resonance"
  Real (Kind=r8_kind), Parameter :: PI = 3.1415926535897932368
  Real (Kind=r8_kind), Parameter :: Mp = 1836.12
  Real (Kind=r8_kind), Parameter :: speed_of_light =  137.03604
  Public gtensor_allocate, gtensor_free, gtensor_calculate

!!$  Logical :: df_read_from_list = .False., df_exclude_atom_orbitals &
!!$       & = .False., df_include_key = .False.
!!$  Integer (Kind=i4_kind) :: df_n_orbitals = 0, &
!!$       & df_n_exclude_atom_orbitals = 0
!!$  Integer (Kind=i4_kind) :: df_dkh_level = 1
!!$  Real (Kind=r8_kind) :: df_mix_coeff = 0.5_r8_kind
  !
  Logical  :: spin_polarized =  .true.
Contains

  function calculate_matrix_element(ev1,gti, ev2) result(c)
    use matrix_module
    type(cmatrix), intent(in)::ev1,ev2,gti
    type(cmatrix):: temp
    complex (Kind=c16_kind):: c
    !code
    call alloc(1,temp)
    temp = mult(tr(ev1),mult(gti,ev2))

    c = cmplx(temp%re(1,1), temp%im(1,1),c16_kind)
    call pgfree(temp)
  end function calculate_matrix_element


  Subroutine reltrafo_gten
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
    !== Interrupt of public interface of module =====================
    !----------------------------------------------------------------
    ! Modifications
    !----------------------------------------------------------------
    !
    ! Modification (Please copy before editing)
    ! Author: ...
    ! Date:   ...
    ! Description: ...
    !----------------------------------------------------------------
    Use symmetry_data_module
    Use quadrupel_fname, Only: qfilename
    use reltrafo, Only: p2_diag
    use matrix_module,dummy=>diag
    Use contraction_module, Only: contract_matrix_spor => contract
    Use dimensions, Only:  dim_irr_arr => IrrUBasDimSpor, dim_irr_c_arr => IrrBasDimSpor

    Implicit None

    Integer (i4_kind) :: dim_irrep, n_irreps, n_partners, i_ir,i_ip, i_xyz
    Integer (i4_kind) :: i_partn, n_orb
    Real (Kind=r8_kind) ::c2, c4, deltaE
    Real (Kind=r8_kind), parameter :: half = 0.5_r8_kind
    Integer (i4_kind) :: i,j,n,m, diagoffdiag
    Type (dipole_integral_type), Pointer :: hfc_p, hfc_cont_p, hfc_tmp
    type(cmatrix)  :: UF,UB, Nuc, PVYP, AVA, ARVRA, DKH1,DKH2, hpp,&
         & GTI, WW, ARVRAt
    type(rdmatrix) :: t_diag,Tp, Ep,Ap,Kp,ApKp,K2p2, K2p2_rev

    !-----------executable code---------------------
    nullify (hfc_p,hfc_cont_p)

    c2 = speed_of_light ** 2
    c4 = speed_of_light ** 4
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

          call alloc (dim_irrep,t_diag, Ep, Tp,Ap, Kp, ApKp,K2p2, k2p2_rev)

          call load_array(dim_irrep,t_diag%d, i_ir,"gt_t_diag" )

          call p2_diag(t_diag,Tp,Ep,Ap,Kp,ApKp,K2p2)

          call alloc(dim_irrep, Nuc, PVYP,AVA, ARVRA, DKH1,DKH2, hpp, GTI)
          call alloc(dim_irrep, WW, ARVRAt)
        !  call alloc(dim_irrep, Upp, Tpp)

          call load_array(dim_irrep,Nuc%re,i_ir,"gt_Nuc_real")
          call load_array(dim_irrep,Nuc%im,i_ir,"gt_Nuc_imag")
          call load_array(dim_irrep,PVYP%re,i_ir,"gt_PVYP_real")
          call load_array(dim_irrep,PVYP%im,i_ir,"gt_PVYP_imag")

          AVA   = mult(Ap,Nuc,Ap)
          ARVRA = mult(ApKp,PVYP,ApKp)
          K2p2_rev%d = 1.0_r8_kind/k2p2%d
          do j=1,dim_irrep
             do i=1,dim_irrep
                DeltaE = Ep%d(i)+Ep%d(j)
                AVA%re(i,j)   = AVA%re(i,j)/DeltaE
                AVA%im(i,j)   = AVA%im(i,j)/DeltaE

                ARVRA%re(i,j) = ARVRA%re(i,j)/DeltaE
                ARVRA%im(i,j) = ARVRA%im(i,j)/DeltaE !FIXME:dg
             enddo
          enddo

          ARVRAt%re = ARVRA%re
          ARVRAt%im = - ARVRA%im


!          do j=1,dim_irrep
!             do i=1,dim_irrep
!                Tpp%re(i,j) = K2p2%d(i) * AVA%re(i,j) - ARVRA%re(i,j)
!                Tpp%im(i,j) = K2p2%d(i) * AVA%im(i,j) - ARVRA%im(i,j)
!
!
!               Upp%re(i,j) = AVA%re(i,j) - ARVRA%re(i,j) * (1_r8_kind/(K2p2%d(j)))
!               Upp%im(i,j) = AVA%im(i,j) - ARVRA%im(i,j) * (1_r8_kind/(K2p2%d(j)))
!
!            end do
!         end do

! Efficient way
 !  Tpp = K2p2 * AVA - ARVRA
 !
 !  Upp = AVA -ARVRA *  K2p2_rev




          n_partners = symmetry_data_n_partners_proj (i_ir)
          npart: Do i_partn = 1, n_partners
             xyz: Do i_xyz = 1, 7 !
                if (i_xyz == 4) cycle

                diag: Do diagoffdiag = 1, 2
                   if (diagoffdiag == 1 )then
                      hfc_tmp => gten_integrals_diag%integrals (i_xyz, i_ip, 1)
                      hfc_p => gten_integrals_rem%integrals (i_xyz, i_ip, 1)
                      hfc_cont_p => gten_integrals_cont_rem%integrals(i_xyz, i_ip, 1)
                      !  Diagonal part

                      Call remap_to_full_store (hfc_tmp,hfc_p)

                      GTI%re(:,:) = hfc_p%offdiagonal(:,:)
                      GTI%im(:,:) = hfc_p%offdiagonal_imag(:,:)
                   else
                      hfc_p => gten_integrals_offdiag%integrals (i_xyz, i_ip, 1)
                      hfc_cont_p => gten_integrals_cont_offdiag%integrals(i_xyz, i_ip, 1)

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
                      !
                      !DKH2 transformation
                      !init

                      Do n = 1, dim_irrep
                         Do m = 1, dim_irrep
                            hpp%re(m,n)= ApKp%d(m)*GTI%re(m,n)*Ap%d(n)/(Ep%d(m)+Ep%d(n))
                            hpp%im(m,n)= ApKp%d(m)*GTI%im(m,n)*Ap%d(n)/(Ep%d(m)+Ep%d(n))
                         end Do
                      end Do

                   !   DKH2 =( - half * mult((mult(Upp,Ep)+mult(Ep,Upp)), hpp) &
                   !        - half * mult(Upp , mult(hpp,Ep)+mult(Ep,hpp)) &
                   !        - half * mult(tr(hpp) , mult(tr(Upp),Ep) + mult(Ep,tr(Upp))) &
                   !        - half * mult(mult(tr(hpp),Ep)+mult(Ep,tr(hpp)),tr(Upp)) &
                   !        + half * mult(hpp , mult(Tpp, Ep) + mult(Ep,Tpp)) &
                   !        + half * mult(mult(hpp,Ep)+mult(Ep,hpp),Tpp) &
                   !        + half * mult(mult(tr(Tpp),Ep) + mult(Ep,tr(Tpp)), tr(hpp)) &
                   !        + half * mult(tr(Tpp),mult(tr(hpp),Ep) + mult(Ep,tr(hpp)))) * speed_of_light

                         WW = - AVA * hpp - hpp * AVA + AVA * k2p2 * hpp + hpp * k2p2 * AVA
                         WW = WW * Ep + Ep * WW
                         DKH2 = WW + 2.0_r8_kind *&
                              ( - AVA * Ep * hpp - hpp * Ep * AVA + AVA * Ep * k2p2 * hpp + hpp * Ep * k2p2 * AVA)


                         if (diagoffdiag ==1 ) then
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


                   Call contract_matrix_spor(i_ir, GTI%im, hfc_cont_p%offdiagonal_imag)


                   if(diagoffdiag==1) then
                      hfc_tmp%diagonal = 0.0_r8_kind
                      Call remap_to_compact_store (hfc_cont_p,hfc_tmp, n_orb)

                   end if

                end Do diag


             end Do xyz
             !
             i_ip=i_ip+1
          end Do npart

          !

          call pgfree(UF,UB, GTI)
          call pgfree(Nuc,PVYP, AVA, ARVRA, DKH1,DKH2, hpp)
          call pgfree(WW, ARVRAt)
         ! call pgfree(Upp, Tpp)
          call pgfree (t_diag, Ep, Tp,Ap, Kp, ApKp,K2p2, k2p2_rev)

       End Do cycle_kf
    End If kf

  end Subroutine reltrafo_gten

  Subroutine gtensor_allocate()
    use hfc_module, only: hfc_init
    Implicit None
    !-----------------------------------

    DPRINT 'gtensor_allocate: call hfc_init()'
    call hfc_init() ! FIXME: initialize dimensions, why to keep them there?

    if(options_kinematic_factors) then
       ! uncontracted dimensions:
       call new(gten_integrals_diag,         10, 1, DIPOLE_DIAGONAL,    UNCONTRACTED, 2)

       ! uncontracted dimensions:
       call new(gten_integrals_offdiag,      10, 1, DIPOLE_OFFDIAGONAL, UNCONTRACTED, 2)

       ! remapped from diagonal storage to offdiagonal:
       call new(gten_integrals_rem,           7, 1, DIPOLE_OFFDIAGONAL, UNCONTRACTED, 2)

       ! for kinematic factors multiplication etc:
       call new(gten_integrals_cont_rem,     7,  1, DIPOLE_OFFDIAGONAL, CONTRACTED,   2)
       call new(gten_integrals_cont_offdiag, 7,  1, DIPOLE_OFFDIAGONAL, CONTRACTED,   2)
    else
       ! contracted dimension:
       call new(gten_integrals_diag,        10,  1, DIPOLE_DIAGONAL,    CONTRACTED,   2)

       ! contracted dimension:
       call new(gten_integrals_offdiag,     10,  1, DIPOLE_OFFDIAGONAL, CONTRACTED,   2)

       ! remapped for matrix multiplication:
       call new(gten_integrals_cont_rem,     7,  1, DIPOLE_OFFDIAGONAL, CONTRACTED,   2)
    end if
  end Subroutine gtensor_allocate

  Subroutine gtensor_free()
    use hfc_module, only: hfc_done
    Implicit None
    !-----------------------------------

    call delete(gten_integrals_diag)
    call delete(gten_integrals_offdiag)
    call delete(gten_integrals_cont_rem)
    if(options_kinematic_factors) then
       call delete(gten_integrals_rem)
       call delete(gten_integrals_cont_offdiag)
    end if

    DPRINT 'gtensor_allocate: call hfc_done()'
    call hfc_done()
  end Subroutine gtensor_free

!!$  Subroutine gtensor_read
!!$    ! purpose: read namelist for gtensors
!!$    !
!!$    ! routine called by: read_input
!!$    !** End of interface *****************************************
!!$    Use input_module
!!$    !------------ Declaration of local variables -----------------
!!$    Integer (Kind=i4_kind) :: unit, status, counter
!!$    Real (Kind=r8_kind) :: norm
!!$    External error_handler
!!$    !------------ Executable code --------------------------------
!!$    read_from_list = df_read_from_list
!!$    n_orbitals = df_n_orbitals
!!$    exclude_atom_orbitals = df_exclude_atom_orbitals
!!$    include_key = df_include_key
!!$    n_exclude_atom_orbitals = df_n_exclude_atom_orbitals
!!$    mix_coeff = df_mix_coeff
!!$    dkh_level = df_dkh_level
!!$    If (input_line_is_namelist("gtensor_orbitals")) Then
!!$       Call input_read_to_intermediate ()
!!$       unit = input_intermediate_unit ()
!!$       Read (Unit, Nml=gtensor_orbitals, IoStat=Status)
!!$       If (status .Gt. 0) Call input_error ("occupation_read: name&
!!$            &list gtensor_orbitals.")
!!$       !
!!$       If (read_from_list) Then
!!$          Allocate (input_irreps(n_orbitals), Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          Allocate (input_partners(n_orbitals), Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          Allocate (input_orbitals(n_orbitals), Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          Do counter = 1, n_orbitals
!!$             Call input_read_to_intermediate
!!$             unit = input_intermediate_unit ()
!!$             Read (Unit, Fmt=*, IoStat=Status) input_irreps &
!!$                  & (counter), input_partners (counter), input_orbitals &
!!$                  & (counter)
!!$             If (status > 0) Call input_error ("occupation_read: n&
!!$                  &oc_start")
!!$          End Do
!!$       End If
!!$       !
!!$       If (exclude_atom_orbitals) Then
!!$          Allocate &
!!$               & (input_exclude_unique_atoms_1(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          Allocate &
!!$               & (input_exclude_atom_l_1(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          Allocate &
!!$               & (input_exclude_atom_n_1(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          Allocate &
!!$               & (input_exclude_unique_atoms_2(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          Allocate &
!!$               & (input_exclude_atom_l_2(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          Allocate &
!!$               & (input_exclude_atom_n_2(n_exclude_atom_orbitals), &
!!$               & Stat=status)
!!$          If (status .Ne. 0) Call error_handler ("dipole_allocate:&
!!$               & allocate of dipole_integrals failed")
!!$          !
!!$          !
!!$          Do counter = 1, n_exclude_atom_orbitals
!!$             Call input_read_to_intermediate
!!$             unit = input_intermediate_unit ()
!!$             Read (Unit, Fmt=*, IoStat=Status) &
!!$                  & input_exclude_unique_atoms_1 (counter), &
!!$                  & input_exclude_atom_l_1 (counter), &
!!$                  & input_exclude_atom_n_1 (counter), &
!!$                  & input_exclude_unique_atoms_2 (counter), &
!!$                  & input_exclude_atom_l_2 (counter), &
!!$                  & input_exclude_atom_n_2 (counter)
!!$             If (status > 0) Call input_error ("occupation_read: n&
!!$                  &oc_start")
!!$          End Do
!!$       End If
!!$       !
!!$    End If
!!$  End Subroutine gtensor_read
!!$
!!$  !*********************************************
!!$  Subroutine gtensor_write (iounit)
!!$    Use echo_input_module
!!$    Use operations_module, Only: operations_echo_input_level
!!$    Integer (Kind=i4_kind) :: status
!!$    Integer, Intent (In) :: iounit
!!$    !
!!$    flag_format = '("    ",a," = ",7x,a5:" # ",a)'
!!$    intg_format = '("    ",a," = ",i12  :" # ",a)'
!!$    real_format = '("    ",a," = ",f12.5:" # ",a)'
!!$    Call start ("GTENSOR_ORBITALS", "GTENSOR_WRITE", iounit, &
!!$         & operations_echo_input_level)
!!$    Call flag ("READ_FROM_LIST ", read_from_list, &
!!$         & df_read_from_list)
!!$    Call intg ("N_ORBITALS     ", n_orbitals, df_n_orbitals)
!!$    Call flag ("EXCLUDE_ATOM_ORBITALS   ", exclude_atom_orbitals, &
!!$         & df_exclude_atom_orbitals)
!!$    Call flag ("INCLUDE_KEY   ", include_key, df_include_key)
!!$    Call intg ("N_EXCLUDE_ATOM_ORBITALS ", &
!!$         & n_exclude_atom_orbitals, df_n_exclude_atom_orbitals)
!!$    Call real ("MIX_COEFF ", mix_coeff, df_mix_coeff)
!!$    !
!!$    Call stop ()
!!$  End Subroutine gtensor_write

  !**********************************************
  Subroutine gtensor_calculate
    !  Purpose: calculates g-tensor components
    !
    !** End of interface *****************************************
    !
    Use filename_module
    Use input_module
    Use unique_atom_module
    Use dimensions, Only: BasDimSpor, SymDimSpor
    Use symmetry_data_module, only: symmetry_data_get_cccoupling
!!$         use matrix_module, only: pgfree,tr,sim,inv,eigs,alloc,cmatrix,rdmatrix, &
!!$              operator(*),assignment(=),operator(-)
    Implicit None

    ! transformations to momentum space and back
    !------------ Declaration of local variables ---------------------
    Integer (Kind=i4_kind) :: i_ir, n_ir, n_pa, n_ip, &
         & i_orb, n_orb, i_xyz !, i_pa
    Integer (Kind=i4_kind) :: i_ip1
    Integer (Kind=i4_kind) :: i, j

    Real (Kind=r8_kind), Parameter :: very_small = 0.000001_r8_kind, zero = 0.0_r8_kind ! precision parameter
    Real (Kind=r8_kind) :: occ_num_current
    Type (g_matrix_element) g_mx_el

    Type (unique_atom_type), Pointer :: ua1, ua2
    Integer (Kind=i4_kind) :: i_ua1, i_ua2, i_l1, i_l2, i_f1, &
         & i_f2, n_f1, n_f2
    Integer (Kind=i4_kind) :: exclude_unique_atoms_1, &
         & exclude_unique_atoms_2, exclude_atom_l_1, exclude_atom_l_2, &
         & exclude_atom_n_1, exclude_atom_n_2
    Integer(i4_kind), pointer :: cc_coupling(:)

    !------------ Executable code ------------------------------------
    !
    DPRINT 'gtensor_calculate: entered'
    spin_orbit: If (.not.options_spin_orbit) then
       call error_handler("gtm::gtensor_calculate: need SO")
    endif spin_orbit
    ! CC (Kramers) Irrep Coupling
    cc_coupling => symmetry_data_get_cccoupling()
    ! Is it non-trivial?:
    spin_polarized = any(cc_coupling.NE.(/(i,i=1,size(cc_coupling))/))
    DPRINT 'gtensor_calculate: spin_polarized=',spin_polarized


       n_ip = symmetry_data_n_ip_proj ()
       n_ir = symmetry_data_n_proj_irreps ()
       !
       i_ip1 = 0
       n_orb = symmetry_data_dimension_proj (1)
       !Searching for the first partner
       g_mx_el%eigenvalue = eigval(1)%m(1, 1)
       !



       ReadList:   If (read_from_list) Then
          Do i = 1, n_orbitals
             g_mx_el%i_irr = input_irreps (i)
             g_mx_el%i_pa = input_partners (i)
             g_mx_el%i_ip = symmetry_data_i_ip_proj &
                  & (g_mx_el%i_irr, g_mx_el%i_pa)
             g_mx_el%i_orb = input_orbitals (i)
             g_mx_el%eigenvalue = &
                  & eigval(g_mx_el%i_irr)%m(g_mx_el%i_orb, 1)
             g_mx_el%n_orb = symmetry_data_dimension_proj &
                  & (g_mx_el%i_irr)

          End Do
       Else ReadList
          ! Find the HOMO == topmost orbital with 1 or 2 electrons,
          ! always takes the first partner within the irrep (trace i_ip1):
          i_ip1 = 1
          irrep: Do i_ir = 1, n_ir
             n_orb = symmetry_data_dimension_proj (i_ir)
             n_pa = symmetry_data_n_partners_proj (i_ir)
!!$             partner1: Do i_pa = 1, n_pa
!!$                i_ip1 = i_ip1 + 1
                orbitals_sog: Do i_orb = n_orb, 1, - 1 !all occ and empty orb`s
                   occ_num_current = occ_num(i_ir)%m(i_orb, 1)
                   If ((occ_num_current < (1+very_small) .And. &
                        & occ_num_current > (1-very_small)) .Or. (occ_num_current &
                        & < (2+very_small) .And. occ_num_current > (2-very_small))) Then
                      ! occ of the ?topmost? is either 1 or 2

                      If (eigval(i_ir)%m(i_orb, 1) > g_mx_el%eigenvalue .or. &
                           abs(eigval(i_ir)%m(i_orb, 1)-g_mx_el%eigenvalue) < very_small) Then
                         ! eigval of the ?topmost? is higher then in other irreps

                         g_mx_el%eigenvalue = eigval(i_ir)%m(i_orb, 1)
                         g_mx_el%i_irr = i_ir
                         g_mx_el%i_pa  = 1
                         g_mx_el%i_ip  = i_ip1 ! linearized (i_irr, i_pa)
                         g_mx_el%i_orb = i_orb
                         g_mx_el%n_orb = n_orb
!!$                         g_mx_el%i_pa = i_pa
                         DPRINT 'gtensor_calculate: Irr =',i_ir
                         DPRINT 'gtensor_calculate: IrPr=',i_ip1
                         DPRINT 'gtensor_calculate: Orb =',i_orb
                         DPRINT 'gtensor_calculate: NOrb=',n_orb
                         DPRINT 'gtensor_calculate: E  =',g_mx_el%eigenvalue
                         exit orbitals_sog ! found the ?topmost?
                      End If
                      !NOOP: Cycle orbitals_sog
                   End If
                End Do orbitals_sog
!!$             End Do partner1
                i_ip1 = i_ip1 + n_pa
          End Do irrep
       end If ReadList
       !
       !
       if (options_kinematic_factors)then
          DPRINT 'gtensor_calculate: call reltrafo_gten()'
          call reltrafo_gten()
          DPRINT 'gtensor_calculate: .'
       endif


       SP:   If(spin_polarized) then
        !  DPRINT 'gtensor_calculate: call spin_polarized_calc_v1()'
        !  call spin_polarized_calc_v1()
           Call spin_restricted_calc()
!!$          DPRINT 'gtensor_calculate: call spin_polarized_calc_v2()'
!!$          call spin_polarized_calc_v2()
       else SP
          DPRINT 'gtensor_calculate: Call spin_restricted_calc()'
          Call spin_restricted_calc()
       end If SP
      ! DPRINT 'gtensor_calculate: call spin_polarized_calc_v2()'
      ! call spin_polarized_calc_v2()

       !
       list: If (read_from_list) Then

          call magnetic_properties_shutdown
       End If list

  Contains


    subroutine gti_sum_expect(spectr,vJ,vLSum11,vLSum12)
      use matrix_module
      implicit none
      type(spectrum), intent(in)     :: spectr(:)
      type(cmatrix), intent(in)      :: vJ(:,:,:)
      real(r8_kind),intent(out)      :: vLSum11(6)
      complex(c16_kind),intent(out)  :: vLSum12(6)
      ! *** end of interface ***

      integer(i4_kind)  :: i_ir,n_ir,i_orb
      real(r8_kind)     :: vecLSigma11(6)
      complex(c16_kind) :: vecLSigma12(6)

      DPRINT 'gti_sum_expect: entered'
      n_ir = size(vJ,3)
      ASSERT(6==size(vJ,1))
      ASSERT(2==size(vJ,2))
      ASSERT(n_ir==size(spectr))

      ! "expectation"  values of g-tensor over eigenfunctions:
      vLSum11 = 0.0
      vLSum12 = 0.0
      irr_: do i_ir=1,n_ir
         do i_orb=1,spectr(i_ir)%n_occ

            call gti_expect(i_orb,vJ(:,:,i_ir),vecLSigma11,vecLSigma12)

            vLSum11 = vLSum11 + vecLSigma11
            vLSum12 = vlSum12 + vecLSigma12
#if _DPRINT

            if(g_mx_el%i_irr == i_ir&
                 & .and.&
                 & g_mx_el%i_orb == spectr(i_ir)%i_orb(i_orb))then
               print *,'****** THIS SEEMS TO BE THE ONE ******'
            endif
            print *,'Irrep=',i_ir,&
                 & ' Orb=', i_orb,&
                 & ' E=', spectr(i_ir)%energy(i_orb)

            call calc_G_tensor_v2(&
                 & vecLSigma11(4:6),vecLSigma12(4:6),&
                 & vecLSigma11(1:3),vecLSigma12(1:3),&
                 & i_ir,&
                 & i_orb,&
                 & spectr(i_ir)%energy(i_orb))
#endif
         enddo
      end do irr_
#if _DPRINT
            DPRINT 'TOTAL:'
            call calc_G_tensor_v2(&
                 & vLSum11(4:6),vLSum12(4:6),&
                 & vLSum11(1:3),vLSum12(1:3),&
                 & 0,&
                 & 0,&
                 & 0.0_r8_kind)
#endif
    end subroutine gti_sum_expect

    subroutine gti_expect(i_orb,vJ,vLS11,vLS12)
      use matrix_module
      implicit none
      integer(i4_kind), intent(in)   :: i_orb
      type(cmatrix), intent(in)      :: vJ(:,:)
      real(r8_kind),intent(out)      :: vLS11(6)
      complex(c16_kind),intent(out)  :: vLS12(6)
      ! *** end of interface ***

      integer(i4_kind)  :: xyz
      real(r8_kind)     :: im

      ASSERT(6==size(vJ,1))
      ASSERT(2==size(vJ,2))

      ! "expectation"  values of g-tensor over eigenfunctions:
      do xyz=1,6
         ASSERT(square(vJ(xyz,1)))
         ASSERT(square(vJ(xyz,2)))
         ASSERT(i_orb<=size(vJ(xyz,1)%re,1))
         ASSERT(i_orb<=size(vJ(xyz,2)%re,1))

         ! check image:
         im = vJ(xyz,1)%im(i_orb,i_orb)
         ASSERT(abs(im)<1e-9)

         ! assign real:
         vLS11(xyz) = vJ(xyz,1)%re(i_orb,i_orb)
         ! this is complex:
         vLS12(xyz) = CMPLX( vJ(xyz,2)%re(i_orb,i_orb), vJ(xyz,2)%im(i_orb,i_orb), r8_kind )
      enddo
    end subroutine gti_expect

    subroutine gti_trafo(out, in, U)
      ! move me out of scope!
      use matrix_module
      implicit none
      type(cmatrix), intent(out)     :: out(:,:,:) ! (6,2,n_irr)
      type(cmatrix), intent(in)      :: in(:,:,:) ! (6,2,n_irr)
      type(cmatrix), intent(in)      :: U(:) ! (n_irr)
      ! *** end of interface ***

      integer(i4_kind) :: n_irr,irr,xyz

      n_irr = size(U)
      ASSERT(all(shape(in)==shape(out)))
      ASSERT(n_irr==size(in,3))
      ASSERT(2==size(in,2))
      ASSERT(6==size(in,1))

      do irr=1,n_irr
         do xyz=1,6
            out(xyz,1,irr) = tr(U(irr),'h') * in(xyz,1,irr) * U(irr)
            out(xyz,2,irr) = tr(U(irr),'t') * in(xyz,2,irr) * U(irr)
         enddo
      enddo
    end subroutine gti_trafo


    Subroutine spin_polarized_calc_v2()
      use matrix_module, dummy=>diag
      implicit none
      !-----------------------
      type(cmatrix)  :: D12
      Integer (Kind=i4_kind)::&
           & i_ir1, i_ir2, dim1, dim2,&
           & counter
      real(r8_kind)      :: vLS11a(6) ! L+Sigma
      complex(c16_kind)  :: vLS12a(6) ! L+Sigma
      real(r8_kind)      :: vLS11b(6) ! L+Sigma
      complex(c16_kind)  :: vLS12b(6) ! L+Sigma
      real(r8_kind)      :: vLS11c(6) ! L+Sigma
      complex(c16_kind)  :: vLS12c(6) ! L+Sigma
      ! matrix elements of VecSigmaXY = <X| vecSigma |Y>
      ! the Hamiltonian is built of scalar product:
      !
      ! | vecSigma11,   vecSigma12 |
      ! | vecSigma12*, -vecSigma11 |
      !
      real(r8_kind), parameter :: zero = 0.0_r8_kind, one=1.0_r8_kind
      integer(i4_kind)         :: n_occ,off,dim0,irr
      integer(i4_kind)         :: memstat
      real(r8_kind)            :: prod
      !
      type(spectrum), allocatable :: spectr(:) ! (n_irr)
      type(cmatrix), allocatable  ::&
           & GTI(:,:,:), & ! (6,2,n_irr), each (dim(irr),dim(irr))
           & vJ(:,:,:), &  ! (6,2,n_irr), each (occ(irr),occ(irr))
           & eigVec(:) ! (n_irr), each (dim(irr),occ(irr))
      type(cmatrix) :: &
           & OVI ! (dim(irr1),dim(irr2)), actualy square
      type(cmatrix),allocatable  :: U(:) ! (n_irr)
      type(rdmatrix),allocatable :: diag(:) ! (n_irr), used every second
      !---------------------------------------------

      n_ir = symmetry_data_n_proj_irreps ()
      allocate(&
           & spectr(n_ir),&
           & GTI(6,2,n_ir),&
           & vJ(6,2,n_ir),&
           & eigVec(n_ir),&
           & STAT=memstat)
      ASSERT(memstat.eq.0)

      prep: do i_ir = 1, n_ir
         n_orb = symmetry_data_dimension_proj(i_ir)
         n_occ = count(occ_num(i_ir)%m(:, 1).GE.very_small)

         spectr(i_ir)%n_occ = n_occ
         allocate(&
              & spectr(i_ir)%energy(n_occ),&
              & spectr(i_ir)%i_orb(n_occ),&
              & STAT=memstat)
         ASSERT(memstat.eq.0)

         call alloc(n_orb,n_occ,eigVec(i_ir))

         counter = 0
         do i_orb = 1, n_orb
            if ( occ_num(i_ir)%m(i_orb, 1) < very_small) cycle
            counter = counter + 1

            spectr(i_ir)%i_orb(counter) = i_orb
            spectr(i_ir)%energy(counter) = eigval(i_ir)%m(i_orb,1)

            eigVec(i_ir)%re(:,counter) = eigvec_real(i_ir)%m(:, i_orb)
            eigVec(i_ir)%im(:,counter) = eigvec_imag(i_ir)%m(:, i_orb)
         enddo
         ASSERT(counter.eq.n_occ)

         ! FIXME: lots of copying here:
         do off = 1, 2
            call gti_map(i_ir, 1, off, GTI(1, off, i_ir)) ! Lx
            call gti_map(i_ir, 2, off, GTI(2, off, i_ir)) ! Ly
            call gti_map(i_ir, 3, off, GTI(3, off, i_ir)) ! Lz
            ! 4 is an overlap?
            call gti_map(i_ir, 5, off, GTI(4, off, i_ir)) ! Sx
            call gti_map(i_ir, 6, off, GTI(5, off, i_ir)) ! Sy
            call gti_map(i_ir, 7, off, GTI(6, off, i_ir)) ! Sz
         enddo
      enddo prep

      ! storage for GTI integrals on eigenfunctions:
      do i_ir=1,n_ir
         dim1 = spectr(i_ir)%n_occ
         do i_xyz=1,6
            call alloc(dim1, vJ(i_xyz,1,i_ir))
            call alloc(dim1, vJ(i_xyz,2,i_ir))
         enddo
      enddo

      ! G-tensor matrix elements over eigenfunctions:
      call gti_trafo(vJ, GTI, eigVec)

#if _DPRINT
      do i_ir=1,n_ir
         i_ir1 = i_ir
         i_ir2 = cc_coupling(i_ir1)
         if(i_ir2<i_ir1) cycle ! already processed

         DPRINT 'EIG: vJ(xyz,off,irr):'
         do i_xyz=1,6
            DPRINT 'xyz=',i_xyz," (off,irr)=...:"
            call print_cmatrix("(1,1)",vJ(i_xyz,1,i_ir1))
            call print_cmatrix("(2,1)",vJ(i_xyz,2,i_ir1))
            call print_cmatrix("(1,2)",vJ(i_xyz,1,i_ir2))
            call print_cmatrix("(2,2)",vJ(i_xyz,2,i_ir2))
         enddo
      enddo
#endif

      ! "expectation"  values of g-tensor over eigenfunctions:
      DPRINT 'G-TENSOR ON EIGENFUNCTIONS:'
      call gti_sum_expect(spectr,vJ,vLS11a,vLS12a)
      ! vLS12 is summed over all, wrong

      ! check if it is a polarized calculation:
      if(.not.spin_polarized) goto 9999 ! cleanup and exit

      ! only every second of diag(:) will be used
      allocate(U(n_ir),diag(n_ir),STAT=memstat)
      ASSERT(memstat==0)

      ! compute overlap of <irr1|K*Irr2>
      svd_: do i_ir=1,n_ir
         i_ir1 = i_ir
         i_ir2 = cc_coupling(i_ir1)
         if(i_ir2<i_ir1) cycle ! already processed

         dim1 = spectr(i_ir1)%n_occ
         dim2 = spectr(i_ir2)%n_occ

         DPRINT 'cross(1): processing irreps ', i_ir1, i_ir2
         DPRINT 'cross(1):               dim=', dim1, dim2

         ! FIXME: copying here:
         call gti_map(i_ir, 4, 2, OVI)

         call alloc(dim1,dim2,D12)

!!$         D12 = mult(tr(eigVec(i_ir1)),mult(OVI,eigVec(i_ir2),"nc"))
!!$         call print_cmatrix("?D12",D12)
!!$         ! FIXME: why not CC?:
!!$         D12 = tr(eigVec(i_ir1)) * OVI * eigVec(i_ir2)
         ! FIXME: Why first CC?:
         D12 = tr(eigVec(i_ir1),"t") * OVI * eigVec(i_ir2)
#if _DPRINT
         call print_cmatrix("v2:D12",D12)
#endif

         call alloc(dim1,U(i_ir1))
         call alloc(dim2,U(i_ir2))
         call alloc(min(dim1,dim2),diag(i_ir))

         call svd(D12,U(i_ir1),diag(i_ir),U(i_ir2))
         ! SVD(D,U,d,VH) : d = UH*D*V;  D = U*d*VH
         ! that is why:
         U(i_ir2) = tr(U(i_ir2))

#if _DPRINT
         DPRINT 'cross: number of zeros=',count(abs(diag(i_ir)%d)<1e-9)
         ! test:
         D12 = tr(U(i_ir1)) * D12 * U(i_ir2)
         call print_cmatrix("U^H*D12*V (diag?)", D12)
#endif
         call pgfree(D12)

      enddo svd_

      !
      ! Transfrom to svd basis:
      !
      call gti_trafo(vJ, (vJ), U)
      !
      ! FIXME: it is no quite kosher to alias intent(in) and intent(out)
      !        parameters!

#if _DPRINT
      do i_ir=1,n_ir
         i_ir1 = i_ir
         i_ir2 = cc_coupling(i_ir1)
         if(i_ir2<i_ir1) cycle ! already processed

         DPRINT 'SVD: vJ(xyz,off,irr):'
         do i_xyz=1,6
            DPRINT 'xyz=',i_xyz," (off,irr)=...:"
            call print_cmatrix("(1,1)",vJ(i_xyz,1,i_ir1))
            call print_cmatrix("(2,1)",vJ(i_xyz,2,i_ir1))
            call print_cmatrix("(1,2)",vJ(i_xyz,1,i_ir2))
            call print_cmatrix("(2,2)",vJ(i_xyz,2,i_ir2))
         enddo
      enddo
#endif

      DPRINT 'G-TENSOR OVER SVD FUNCITONS (dismiss energies!)'
      call gti_sum_expect(spectr,vJ,vLS11b,vLS12b)
      ! vecLSigma12 is summed over all, wrong
      DPRINT 'UNTRANS/TRANS DIFF LS11=',maxval(abs(vLS11a-vLS11b))
      DPRINT 'UNTRANS/TRANS DIFF LS12=',maxval(abs(vLS12a-vLS12b))

      ! clean
      do i_ir=1,size(U)
         call pgfree(U(i_ir))
      enddo

      prod = zero
      corrorb_: do i_ir=1,n_ir
         i_ir1 = i_ir
         i_ir2 = cc_coupling(i_ir1)
         if( i_ir2<i_ir1 ) cycle ! already processed

         dim1 = spectr(i_ir1)%n_occ
         dim2 = spectr(i_ir2)%n_occ

         dim0 = count(abs(diag(i_ir)%d)>=1e-9) ! FIXME: small number here
         DPRINT 'cross(2): number of non-zeros=',dim0
         ASSERT(dim0==min(dim1,dim2))
         if(dim0/=max(dim1,dim2))then
            ! here the unpaired electron
            ASSERT(abs(dim2-dim1)==1)
            ! make sure there is, only one unpaired:
            ASSERT(prod==zero)
            if(dim2>dim1)then
               irr   = i_ir2
               i_orb = dim2
            else
               irr   = i_ir1
               i_orb = dim1
            endif

            DPRINT 'G-TENSOR OVER UNPAIRED EIGENFUNCITON'
            call gti_expect(i_orb,vJ(:,:,irr),vLS11c,vLS12c)

            call calc_G_tensor_v2(&
                 & vLS11c(4:6),vLS12c(4:6),&
                 & vLS11c(1:3),vLS12c(1:3),&
                 & 0,&
                 & 1,&
                 & -1.0_r8_kind)

            prod = product(diag(i_ir)%d)
            DPRINT 'd_ii=',diag(i_ir)%d
            DPRINT 'prod=',prod,' prod^2=',prod*prod
         endif
      enddo corrorb_

      DPRINT 'prod=',prod,' prod^2=',prod*prod
      prod = prod*prod
      DPRINT 'G-TENSOR AS SUGGESTED BY DG'
      call calc_G_tensor_v2(&
           & vLS11a(4:6), vLS12c(4:6),&
           & vLS11a(1:3), vLS12c(1:3)*prod,&
           & 0,&
           & 2,&
           & -2.0_r8_kind)
      DPRINT 'G-TENSOR SHOULD BE AS BEFORE'
      call calc_G_tensor_v2(&
           & vLS11b(4:6), vLS12c(4:6),&
           & vLS11b(1:3), vLS12c(1:3)*prod,&
           & 0,&
           & 3,&
           & -3.0_r8_kind)
      DPRINT 'G-TENSOR FORMALLY FROM DETERMINANTS'
      call calc_G_tensor_v2(&
           & vLS11a(4:6), vLS12c(4:6)*prod,&
           & vLS11a(1:3), vLS12c(1:3)*prod,&
           & 0,&
           & 4,&
           & -4.0_r8_kind)

      ! clean:
      do i_ir=1,n_ir
         i_ir1 = i_ir
         i_ir2 = cc_coupling(i_ir1)
         if( i_ir2<i_ir1 ) cycle ! already processed
         call pgfree(diag(i_ir))
      enddo

      ! clean up and exit
9999 continue
      do i_ir=1,n_ir
         deallocate(spectr(i_ir)%energy,spectr(i_ir)%i_orb,STAT=memstat)
         ASSERT(memstat.eq.0)
         do off=1,2
            do i_xyz=1,6
               call pgfree(vJ(i_xyz,off,i_ir))
            enddo
         enddo
         call pgfree(eigVec(i_ir))
      enddo
      deallocate(spectr,GTI,eigVec,STAT=memstat)
      ASSERT(memstat.eq.0)
    end Subroutine spin_polarized_calc_v2

#if _DPRINT
    subroutine print_vecHam(name,v11,v12)
      implicit none
      character(len=*), intent(in)           :: name
      real(r8_kind), intent(in)              :: v11(3)
      complex(c16_kind), intent(in),optional :: v12(3)
      ! *** end of interfacce ***

!!$      write(*,'(A)')  name//':'
         write(*,'(A,T15,3F15.10)') trim(name)//" 11", v11
      if(present(v12))then
         write(*,'(A,T15,3F15.10)') trim(name)//" 12.re", real(v12)
         write(*,'(A,T15,3F15.10)') trim(name)//" 12.im", aimag(v12)
      endif
    end subroutine print_vecHam
#endif

Subroutine spin_polarized_calc_v1()
  use matrix_module
  !-----------------------
  type(cmatrix)  :: eigvectors1,eigvectors2,eigvectors_trans1,eigvectors_trans2,&
       D12, D21 , DD1,DD2, U1,U2, V1,V2,  D_diag12, eigvector_trans1,eigvector_trans2,S, U_vector1,&
       U_vector2, scalar_product, U_vector_sum, lambda
  type(rdmatrix) :: U_diag1,V_diag2, V_diag1, U_diag2
  ! Works only for one dementional irreps
  Real (Kind=r8_kind)::Sigma_diag(3), L_diag(3), sigma_real_offdiag(3),&
       Sigma_im_offdiag(3),L_real_offdiag(3), L_im_offdiag(3), product, &
       Sigma_real_offdiag_sum(3),Sigma_im_offdiag_sum(3),L_real_offdiag_sum(3), L_im_offdiag_sum(3)
  Integer (Kind=i4_kind):: i_orb1, i_orb2, n_occ_arr(8), counter, &
       i_ir1, i_ir2, dim1, dim2, dim_sum
  Type (dipole_integral_type), Pointer:: hfc
  Real (Kind=r8_kind)::mx_re, mx_im
  Type (dipole_integral_type), Pointer:: hfcint, hfc_tmp

  complex(kind = c16_kind)::c

  real(r8_kind)     :: vecLSigma11(6),vLSum11(6)
  complex(c16_kind) :: vLSum12(6)
  !---------------------------------------------
  ! first diagonal part
  Sigma_diag = 0.0_r8_kind
  L_diag = 0.0_r8_kind
  sigma_real_offdiag_sum = 0.0_r8_kind
  Sigma_im_offdiag_sum = 0.0_r8_kind
  L_real_offdiag_sum = 0.0_r8_kind
  L_im_offdiag_sum = 0.0_r8_kind

  vLSum11 = 0.0_r8_kind
  vlSum12 = 0.0_r8_kind

  n_ir = symmetry_data_n_proj_irreps ()

  irr_diag:  do i_ir = 1, n_ir
     n_orb = symmetry_data_dimension_proj(i_ir)
     call alloc(n_orb,1, eigvector_trans1)
     call alloc(n_orb, S)
     counter = 0
     do i_orb = 1, n_orb
        if ( occ_num(i_ir)%m(i_orb, 1) < very_small) cycle

        eigvector_trans1%re(:,1)= eigvec_real(i_ir)%m(:, i_orb)
        eigvector_trans1%im(:,1)= eigvec_imag(i_ir)%m(:, i_orb)
        do i_xyz = 1,7
           if(i_xyz==4) cycle


           If (  options_kinematic_factors) Then

              hfc_tmp=> gten_integrals_cont_rem%integrals(i_xyz,i_ir, 1)
              S%re(:,:) = hfc_tmp%offdiagonal(:,:)
              S%im(:,:) = hfc_tmp%offdiagonal_imag(:,:)
           else

              hfcint => gten_integrals_diag%integrals (i_xyz, i_ir , 1)
              hfc_tmp=> gten_integrals_cont_rem%integrals(i_xyz, i_ir , 1)
              Call remap_to_full_store (hfcint,hfc_tmp)
              S%re(:,:) = hfc_tmp%offdiagonal(:,:)
              S%im(:,:) = hfc_tmp%offdiagonal_imag(:,:)
           end If

!!$           ! check Hermitean?
!!$           mx_im = maxval(abs(S-tr(S)))
!!$           DPRINT 'hermitean? =',mx_im


           c = calculate_matrix_element(eigvector_trans1,S,eigvector_trans1)
           mx_re = Real(c)
           mx_im = Aimag(c)
!!$           DPRINT 'mx_im=',mx_im
           ASSERT(abs(mx_im)<1e-9)

           select case(i_xyz)
           case (1)
              L_diag(1) = L_diag(1) + mx_re
              vecLSigma11(1) = mx_re
           case(2)
              L_diag(2) = L_diag(2) + mx_re
              vecLSigma11(2) = mx_re
           case(3)
              L_diag(3) = L_diag(3) + mx_re
              vecLSigma11(3) = mx_re

           case(5)
              Sigma_diag(1) =  Sigma_diag(1) + mx_re
              vecLSigma11(4) = mx_re
           case(6)
              Sigma_diag(2) =  Sigma_diag(2) + mx_re
              vecLSigma11(5) = mx_re
           case(7)
              Sigma_diag(3) =  Sigma_diag(3) + mx_re
              vecLSigma11(6) = mx_re
           end select

        end do

        vLSum11 = vLSum11 + vecLSigma11
#if _DPRINT
        DPRINT 'v1: irr=',i_ir,' orb=',i_orb
        call print_vecHam("Sigma",vecLSigma11(4:6))
        call print_vecHam("L    ",vecLSigma11(1:3))
#endif

        counter = counter + 1
     end do
     n_occ_arr(i_ir)= counter
     call pgfree(eigvector_trans1)
     call pgfree(S)
  end do irr_diag

#if _DPRINT
  DPRINT 'v1: TOTAL'
  call print_vecHam("Sigma",vLSum11(4:6))
  call print_vecHam("L    ",vLSum11(1:3))
#endif


  ! Calculate offdiagonal part
  !load matrix before
  do i_ir = 1, n_ir

     n_orb = symmetry_data_dimension_proj(i_ir)

     If (  options_kinematic_factors) Then

        hfc => gten_integrals_cont_offdiag%integrals (4, i_ir , 1)
     else

        hfc => gten_integrals_offdiag%integrals (4, i_ir , 1)

     end If


  end do
  !For 2-irrep symmetries only!
  product = 1.0_r8_kind
  Sigma_real_offdiag = 0.0_r8_kind
  Sigma_im_offdiag = 0.0_r8_kind
  L_real_offdiag = 0.0_r8_kind
  L_im_offdiag = 0.0_r8_kind
  irr_loop:   do i_ir = 1, n_ir,2 !pairs of irreps
     i_ir1 = i_ir
     i_ir2 = i_ir +1
     if(n_occ_arr(i_ir1) == 0 .and. n_occ_arr(i_ir2) == 0 ) cycle
     n_orb = symmetry_data_dimension_proj(i_ir)

     dim_sum = n_occ_arr(i_ir1) + n_occ_arr(i_ir2)


     dim1 = n_occ_arr(i_ir1)
     dim2 = n_occ_arr(i_ir2)
     if(dim1 == dim2) cycle
     call alloc(n_orb, S)
     call alloc(dim1,dim2, D12)
     call alloc(dim2, dim1, D21)
     call alloc(n_orb,dim1, eigvectors1)
     call alloc(n_orb,dim2, eigvectors2)
     call alloc(n_orb,dim1, eigvectors_trans1)
     call alloc(n_orb,dim2, eigvectors_trans2)
     call alloc(n_orb,1, eigvector_trans1)
     call alloc(n_orb,1, eigvector_trans2)
     call alloc(dim1, DD1)
     call alloc(dim2, DD2)
     call alloc(dim1, U1)
     call alloc(dim2, V1)
     call alloc(dim1, U2)
     call alloc(dim2, V2)
     call alloc(dim1, U_diag2)
     call alloc(dim2, V_diag2)
     call alloc(dim2, V_diag1)
     call alloc(dim1, U_diag1)
     call alloc(dim1,dim2, D_diag12)
     call alloc(dim_sum,1, U_vector1)
     call alloc(dim_sum,1, U_vector2)
     call alloc(dim_sum,1, U_vector_sum)
     call alloc(1,1,scalar_product)
     call alloc(dim2,dim1,lambda)




     DPRINT "v1:processing irreps  ", i_ir1, i_ir2
     DPRINT "v1:dimension = ", dim1, dim2


     do i_orb1 = 1, dim1 !eigenvector1 - Kramers conj.

        eigvectors1%re(:,i_orb1) =  eigvec_real(i_ir1)%m(:, i_orb1)
        eigvectors1%im(:,i_orb1) =  - eigvec_imag(i_ir1)%m(:, i_orb1)
     end do
     do i_orb2 = 1, dim2

        eigvectors2%re(:, i_orb2) =  eigvec_real(i_ir2)%m(:, i_orb2)
        eigvectors2%im(:, i_orb2) =  eigvec_imag(i_ir2)%m(:, i_orb2)
     end do


     If (  options_kinematic_factors) Then


        hfcint => gten_integrals_cont_offdiag%integrals (4, i_ir1 , 1)
     else

        hfcint => gten_integrals_offdiag%integrals (4, i_ir1 , 1)

     end If

     S%re(:,:) = hfcint%offdiagonal(:,:)
     S%im(:,:) = hfcint%offdiagonal_imag(:,:)


     D12 = mult(tr(eigvectors1),mult(S,eigvectors2))


!!$     DPRINT "Real part D12"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D12%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part D12"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D12%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:D12",D12)
#endif

     D21 = tr(D12)      !Constructed transformed matrix
     D21%im(:,:) = - D21%im(:,:)

!!$     DPRINT "Real part D21"
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') D21%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part D21"
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') D21%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:D21",D21)
#endif

     DD2 = mult(tr(D12),D12)
     DD1 = mult(tr(D21),D21)

!!$     DPRINT "Real part D12 * D12"
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') DD2%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part D12 * D12"
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') DD2%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:D12^H*D12",DD2)
#endif

!!$     DPRINT "Real part D21 * D21"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') DD1%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part D21 * D21"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') DD1%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:D21^H*D21",DD1)
#endif

     call eigs(DD2,V_diag2,V2)
     call eigs(DD1,U_diag2,U2)


     do i_orb2 = 1,dim2
        V1%re(:,dim2 - i_orb2 + 1) = V2%re(:,i_orb2)
        V1%im(:,dim2 - i_orb2 + 1) = V2%im(:,i_orb2)
        V_diag1%d(dim2 - i_orb2 + 1) = V_diag2%d(i_orb2)
     end do

     do i_orb2 = 1,dim1
        U1%re(:,dim1 - i_orb2 + 1) = U2%re(:,i_orb2)
        U1%im(:,dim1 - i_orb2 + 1) = U2%im(:,i_orb2)
        U_diag1%d(dim1 - i_orb2 + 1) = U_diag2%d(i_orb2)
     end do

     DPRINT " eigenvalues of D12 * D12"
     write(*,'(15F10.6)') V_diag1%d(:)

!!$     DPRINT "eigenvectors of D12 * D12 Real part "
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') V1%re(i_orb2,:)
!!$     end do
!!$     DPRINT "eigenvalues of D12 * D12 Imag part "
!!$     do i_orb2 = 1, dim2
!!$        write(*,'(15F10.6)') V1%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:eigenvectors of D12^H*D12",V1)
#endif

     DPRINT " eigenvectors of D21 * D21"
     write(*,'(15F10.6)') U_diag1%d(:)

!!$     DPRINT "Eigenvectors of  D21 * D21  Real part "
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') U1%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Eigenvectors of  D21 * D21 Imag part"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') U1%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:eigenvectors of D21^H*D21",U1)
#endif


     do i_orb1 = 1, dim1
        lambda%re(i_orb1, i_orb1) = 1/Sqrt(V_diag1%d(i_orb1))
     end do
     DPRINT "Real part lambda"
     do i_orb2 = 1, dim2
        write(*,'(15F10.6)') lambda%re(i_orb2,:)
     end do
     U1 = mult(mult(D12,V1), lambda)

!!$     DPRINT "Real part U1"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') U1%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part U1"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') U1%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:U1",U1)
#endif

     D_diag12 = mult(tr(U1), mult(D12, V1))


!!$     DPRINT "Real part D_diag12"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D_diag12%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part D_diag12"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D_diag12%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:D_diag12",D_diag12)
#endif

     !check the othogonality
     DD1 = mult(tr(U1),U1)
!!$     DPRINT "Real part trU1U1"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') DD1%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part trU1U1"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') DD1%im(i_orb2,:)
!!$     end do
#if _DPRINT
     call print_cmatrix("v1:U1^H*U1",DD1)
#endif

     !   Transformed functions
     eigvectors_trans1 = mult(eigvectors1, U1)
     eigvectors_trans2 = mult(eigvectors2, V1)

     !Check the overlap of the transformed functions

     D12 = mult(tr(eigvectors_trans1),mult(S,eigvectors_trans2))


!!$     DPRINT "Real part of the overlap of the transformed functions"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D_diag12%re(i_orb2,:)
!!$     end do
!!$     DPRINT "Imag part of the overlap of the transformed functions"
!!$     do i_orb2 = 1, dim1
!!$        write(*,'(15F10.6)') D_diag12%im(i_orb2,:)
!!$     end do
#if _DPRINT
    call print_cmatrix("v1:transf. overlap",D_diag12)
#endif


     if(dim1 /= dim2 ) then
        do i_orb2 = 1 , dim1
           product = product * D_diag12%re(i_orb2 ,i_orb2)
        end do
     else
        do i_orb2 = 1, dim1
           product = product * D_diag12%re(i_orb2 ,i_orb2)
        end do
     end if

     DPRINT "Product = ", product, " Product^2 = ", product*product


     gtencalc:if (dim2 /= dim1 ) then !the SOMO is in the irrep 2

        eigvector_trans1%re(:,1) = eigvectors_trans2%re(:,dim2)
        eigvector_trans1%im(:,1) = -eigvectors_trans2%im(:,dim2)
        ! eigvector_trans1%im(:,1) = eigvectors_trans2%im(:,dim2)
        eigvector_trans2%re(:,1) = eigvectors_trans2%re(:,dim2)
        eigvector_trans2%im(:,1) = eigvectors_trans2%im(:,dim2)



        xyz:   do i_xyz = 1,7
           if(i_xyz == 4) cycle
           if (dim1 > dim2 ) then
              If (  options_kinematic_factors) Then


                 hfcint => gten_integrals_cont_offdiag%integrals (i_xyz, i_ir1 , 1)
              else

                 hfcint => gten_integrals_offdiag%integrals (i_xyz, i_ir1 , 1)

              end If
           else
              If (  options_kinematic_factors) Then


                 hfcint => gten_integrals_cont_offdiag%integrals (i_xyz, i_ir2 , 1)
              else

                 hfcint => gten_integrals_offdiag%integrals (i_xyz, i_ir2 , 1)

              end If
           end if
           S%re(:,:) = hfcint%offdiagonal(:,:)
           S%im(:,:) = hfcint%offdiagonal_imag(:,:)

           c = calculate_matrix_element(eigvector_trans1,S,eigvector_trans2)

           mx_re = Real(c)
           mx_im = Aimag(c)

           select case(i_xyz)
           case (1)
              L_real_offdiag(1) =  mx_re !* product * product
              L_im_offdiag(1) = mx_im !* product * product

           case(2)

              L_real_offdiag(2) = mx_re !* product * product
              L_im_offdiag(2) = mx_im !* product * product

           case(3)
              L_real_offdiag(3) = mx_re !* product * product
              L_im_offdiag(3) = mx_im !* product * product

           case(5)

              Sigma_real_offdiag(1) =  mx_re
              Sigma_im_offdiag(1) =  mx_im


           case(6)


              Sigma_real_offdiag(2) =  mx_re
              Sigma_im_offdiag(2) =  mx_im

           case(7)

              Sigma_real_offdiag(3) =  mx_re
              Sigma_im_offdiag(3) =  mx_im

           end select

        end do xyz

     end if gtencalc

     call pgfree(S,D12,D21, U1, V1,U2,V2,DD1, DD2)
     call pgfree(eigvectors1)
     call pgfree(eigvectors2)
     call pgfree(eigvectors_trans1)
     call pgfree(eigvectors_trans2)
     call pgfree(eigvector_trans1)
     call pgfree(eigvector_trans2)
     call pgfree(D_diag12)
     call pgfree(V_diag2)
     call pgfree(V_diag1)
     call pgfree(U_diag1)
     call pgfree(U_diag2)
     call pgfree(U_vector1)
     call pgfree(U_vector2)
     call pgfree(U_vector_sum)
     call pgfree(scalar_product)
     call pgfree(lambda)
  end do irr_loop


#if _DPRINT
  g_mx_el%Sigma_Real_diag = Sigma_diag
  g_mx_el%L_Real_diag = L_diag
  g_mx_el%Sigma_Im_diag = 0.0_r8_kind
  g_mx_el%L_Im_diag = 0.0_r8_kind

  g_mx_el%Sigma_real_offdiag = Sigma_real_offdiag
  g_mx_el%Sigma_im_offdiag = Sigma_im_offdiag
  g_mx_el%L_real_offdiag = L_real_offdiag
  g_mx_el%L_im_offdiag = L_im_offdiag
  Call calc_G_tensor(g_mx_el)
#endif

  g_mx_el%Sigma_Real_diag = Sigma_diag
  g_mx_el%L_Real_diag = L_diag
  g_mx_el%Sigma_Im_diag = 0.0_r8_kind
  g_mx_el%L_Im_diag = 0.0_r8_kind

  g_mx_el%Sigma_real_offdiag = Sigma_real_offdiag
  g_mx_el%Sigma_im_offdiag = Sigma_im_offdiag
  g_mx_el%L_real_offdiag = L_real_offdiag * product * product
  g_mx_el%L_im_offdiag = L_im_offdiag     * product * product
  Call calc_G_tensor(g_mx_el)

end Subroutine spin_polarized_calc_v1


    Subroutine calc_matrix_elements
      Implicit None
      Real (Kind=r8_kind), Pointer :: eigenvector_real1 (:, :), &
           & eigenvector_imag1 (:, :)
      Real (Kind=r8_kind), Dimension (:), Pointer :: &
           & eigenvector1_real, eigenvector1_imag

      Real (Kind=r8_kind) ::tm_r, tm_i
      Type (dipole_integral_type), Pointer:: hfcint

      real(r8_kind)     :: vecSigma11(3)
      complex(c16_kind) :: vecSigma12(3)
      real(r8_kind)     :: vecL11(3)
      complex(c16_kind) :: vecL12(3)
      !---------------------------------------------------------

      DPRINT 'calc_matrix_elements: entered'

      !

      eigenvector_real1 => eigvec_real(g_mx_el%i_irr)%m(:, :)
      eigenvector_imag1 => eigvec_imag(g_mx_el%i_irr)%m(:, :)
      !
      eigenvector1_real => eigenvector_real1 (:, g_mx_el%i_orb)
      eigenvector1_imag => eigenvector_imag1 (:, g_mx_el%i_orb)


      !
      sxyz: Do i_xyz = 1, 7
         if (i_xyz==4) cycle
         okf: If (  options_kinematic_factors) Then

            ! hfcint => gten_integrals_cont_rem%integrals (i_xyz, 1, g_mx_el%i_ip , 1)
            hfcint => gten_integrals_diag%integrals (i_xyz, g_mx_el%i_ip , 1)
         else
            hfcint => gten_integrals_diag%integrals (i_xyz, g_mx_el%i_ip , 1)

         end If okf
         Call calc_write_tmgt (hfcint, g_mx_el%n_orb,g_mx_el%i_irr,i_xyz, eigenvector1_real, eigenvector1_imag, tm_r, tm_i )

         ! diagonal is real:
         ASSERT(tm_i<1e-7)

         select case(i_xyz)
         case (5:7)
            g_mx_el%Sigma_Real_diag(i_xyz-4) = tm_r
            g_mx_el%Sigma_Im_diag(i_xyz-4) = tm_i
            vecSigma11(i_xyz-4) = tm_r
         case (1:3)
            g_mx_el%L_Real_diag(i_xyz) = tm_r
            g_mx_el%L_Im_diag(i_xyz) = tm_i
            vecL11(i_xyz) = tm_r
         case default
            ABORT("no such case")
            !
            !
         end select

         If (  options_kinematic_factors) Then

            hfcint => gten_integrals_cont_offdiag%integrals (i_xyz, g_mx_el%i_ip , 1)

         else

            hfcint => gten_integrals_offdiag%integrals (i_xyz, g_mx_el%i_ip , 1)
         end If

         Call calc_write_tmgt(&
              & hfcint,&
              & g_mx_el%n_orb,&
              & g_mx_el%i_irr,&
              & i_xyz,&
              & eigenvector1_real,&
              & eigenvector1_imag,&
              & tm_r, tm_i )
         select case(i_xyz)
         case (5:7)
            g_mx_el%Sigma_Real_offdiag(i_xyz-4) = tm_r
            g_mx_el%Sigma_Im_offdiag(i_xyz-4) = tm_i
            vecSigma12(i_xyz-4) = cmplx(tm_r,tm_i,r8_kind)
         case (1:3)
            g_mx_el%L_Real_offdiag(i_xyz) = tm_r
            g_mx_el%L_Im_offdiag(i_xyz) = tm_i
            vecL12(i_xyz) = cmplx(tm_r,tm_i,r8_kind)
         case default
            ABORT("no such case")
         end select

      end Do sxyz

#if _DPRINT
            print *,'Irrep=',g_mx_el%i_irr,&
                 & ' Orb=',g_mx_el%i_orb ,' E=',&
                 & eigval(g_mx_el%i_irr)%m(g_mx_el%i_orb, 1)
            call print_vecHam("Sigma",vecSigma11,vecSigma12)
            call print_vecHam("L",vecL11,vecL12)
!!$            print *,' vecSigma11    = ',vecSigma11
!!$            print *,' vecL11        = ',vecL11
!!$            print *,' vecSigma12.re = ',Real(vecSigma12)
!!$            print *,' vecSigma12.im = ',aImag(vecSigma12)
!!$            print *,' vecL12.re     = ',Real(vecL12)
!!$            print *,' vecL12.im     = ',aImag(vecL12)
#endif

    End Subroutine calc_matrix_elements

    Subroutine spin_restricted_calc()
      Real (Kind=r8_kind), Parameter :: cut_off_parameter = 10e-5
      Real (Kind=r8_kind) :: max_elem, min_elem
      !  If (options_kinematic_factors) Call &
           ! & multiply_kinematic_factors (g_mx_el%i_irr, g_mx_el%i_ip, &
      ! & g_mx_el%i_ip, g_mx_el%n_orb)
      !

      !
      If (exclude_atom_orbitals) Then
         If (n_exclude_atom_orbitals == 1 .And. &
              & (input_exclude_unique_atoms_1(1) == 0 .And. &
              & input_exclude_atom_l_1(1) == 0 .And. &
              & input_exclude_atom_n_1(1) == 0)) Then
            Write (output_unit,*) 'Atom shell contributions to L:'
            !
            Do i_ua1 = 1, N_unique_atoms
               Do i_ua2 = 1, N_unique_atoms
                  ua1 => unique_atoms (i_ua1)
                  ua2 => unique_atoms (i_ua2)
                  Do i_l1 = 0, ua1%lmax_ob
                     !
                     Do i_l2 = 0, ua2%lmax_ob
                        n_f1 = BasDimSpor(i_ua1)%LM(i_l1, &
                             & g_mx_el%i_irr) / &
                             & SymDimSpor(i_ua1)%LM(i_l1, &
                             & g_mx_el%i_irr)!number of independent functions
                        n_f2 = BasDimSpor(i_ua2)%LM(i_l2, &
                             & g_mx_el%i_irr) / &
                             & SymDimSpor(i_ua2)%LM(i_l2, &
                             & g_mx_el%i_irr)
                        !
                        Do i_f1 = 1, n_f1
                           Do i_f2 = 1, n_f2
                              exclude_unique_atoms_1 = i_ua1
                              exclude_unique_atoms_2 = i_ua2
                              exclude_atom_l_1 = i_l1
                              exclude_atom_l_2 = i_l2
                              exclude_atom_n_1 = i_f1
                              exclude_atom_n_2 = i_f2
                              !
                              Call calc_matrix_elements
                              !
                              !
                              Write (output_unit, Fmt='(6(A,I3,2x))') "ua1 = ", i_ua1, "n1 = ", &
                                   & i_f1, "l1 = ", i_l1, "ua2 = ", &
                                   & i_ua2, "n2 = ", i_f2, "l2 = ", &
                                   & i_l2
                              !
                              !Cut off very small contributions
                              max_elem = Max (MAXVAL(g_mx_el%Sigma_Real_offdiag(:)), MAXVAL(g_mx_el%Sigma_Im_offdiag(:)), &
                                   & MAXVAL(g_mx_el%L_Real_offdiag(:)), MAXVAL(g_mx_el%L_Im_offdiag(:)), &
                                   & MAXVAL(g_mx_el%Sigma_Real_diag(:)), MAXVAL(g_mx_el%Sigma_Im_diag(:)), &
                                   & MAXVAL(g_mx_el%L_Real_diag(:)), MAXVAL(g_mx_el%L_Im_diag(:)))
                              !
                              min_elem = Min (MINVAL(g_mx_el%Sigma_Real_offdiag(:)), MINVAL(g_mx_el%Sigma_Im_offdiag(:)), &
                                   & MINVAL(g_mx_el%L_Real_offdiag(:)), MINVAL(g_mx_el%L_Im_offdiag(:)), &
                                   & MINVAL(g_mx_el%Sigma_Real_diag(:)), MINVAL(g_mx_el%Sigma_Im_diag(:)), &
                                   & MINVAL(g_mx_el%L_Real_diag(:)), MINVAL(g_mx_el%L_Im_diag(:)))
                              !
                              If ((Abs(max_elem) < &
                                   & cut_off_parameter) .Or. &
                                   & (Abs(min_elem) < &
                                   & cut_off_parameter)) Then
                                 Write (output_unit,*) 'No contribution '
                                 Cycle
                              End If
                              !
                              Call calc_G_tensor(g_mx_el)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         Else If ((input_exclude_unique_atoms_1(1) .Ne. 0 .Or. &
              & input_exclude_atom_l_1(1) .Ne. 0 .Or. &
              & input_exclude_atom_n_1(1) .Ne. 0) .And. &
              & (input_exclude_unique_atoms_2(1) == 0 .And. &
              & input_exclude_atom_l_2(1) == 0 .And. &
              & input_exclude_atom_n_2(1) == 0)) Then
            !
            !
            Do j = 1, n_exclude_atom_orbitals
               Do i = 1, j
                  !
                  exclude_unique_atoms_1 = &
                       & input_exclude_unique_atoms_1 (i)
                  exclude_unique_atoms_2 = &
                       & input_exclude_unique_atoms_1 (j)
                  exclude_atom_l_1 = input_exclude_atom_l_1 (i)
                  exclude_atom_l_2 = input_exclude_atom_l_1 (j)
                  exclude_atom_n_1 = input_exclude_atom_n_1 (i)
                  exclude_atom_n_2 = input_exclude_atom_n_1 (j)
                  Write (output_unit,*) 'Exctracted contribution '
                  Write (output_unit, Fmt='(6(A,I3,2x))') "ua1 = &
                       &", exclude_unique_atoms_1, "n1 = ", &
                       & exclude_atom_n_1, "l1 = ", exclude_atom_l_1, "&
                       &ua2 = ", exclude_unique_atoms_2, "n2 = ", &
                       & exclude_atom_n_2, "l2 = ", exclude_atom_l_2
                  !
                  !
                  Call calc_matrix_elements ()
                  Call calc_G_tensor(g_mx_el)
                  !
               End Do
            End Do
            Write (output_unit,*) 'All contribution '
            exclude_atom_orbitals = .False.
            include_key = .True.
            Call calc_matrix_elements ()
            Call calc_G_tensor(g_mx_el)
         Else
            Do i = 1, n_exclude_atom_orbitals
               exclude_unique_atoms_1 = &
                    & input_exclude_unique_atoms_1 (i)
               exclude_unique_atoms_2 = &
                    & input_exclude_unique_atoms_2 (i)
               exclude_atom_l_1 = input_exclude_atom_l_1 (i)
               exclude_atom_l_2 = input_exclude_atom_l_2 (i)
               exclude_atom_n_1 = input_exclude_atom_n_1 (i)
               exclude_atom_n_2 = input_exclude_atom_n_2 (i)
               Write (output_unit,*) 'Exctracted contribution '
               Write (output_unit, Fmt='(6(A,I3,2x))') "ua1 = ", &
                    & exclude_unique_atoms_1, "n1 = ", &
                    & exclude_atom_n_1, "l1 = ", exclude_atom_l_1, "ua2&
                    & = ", exclude_unique_atoms_2, "n2 = ", &
                    & exclude_atom_n_2, "l2 = ", exclude_atom_l_2
               !
               !
               Call calc_matrix_elements ()
               Call calc_G_tensor(g_mx_el)
            End Do
         End If
         !
         !
      Else
         Call calc_matrix_elements ()
         Call calc_G_tensor(g_mx_el)
      End If
    End Subroutine spin_restricted_calc

    Subroutine calc_G_tensor(g_mx_el)
      !
      Type (g_matrix_element), intent(in) :: g_mx_el
      ! *** end of interface ***

      Real (Kind=r8_kind), Dimension (3, 3) :: g_Sigma, g_L, g
      Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      integer :: unit
      !--------------------------------------------------------------------------
#if _DPRINT
      unit = 6
#else
      unit = output_unit
#endif
!!$      Write (output_unit,*) 'Kinematic factors is ', &
!!$           & options_kinematic_factors
      Write (unit, Fmt='(4(A,I3,2x),A,F14.7,A)') "irrep: ",&
           &  g_mx_el%i_irr, "partner: ", g_mx_el%i_pa, "irrep-partner:&
           & ", g_mx_el%i_ip, "orb: ", g_mx_el%i_orb, " E  =", &
           & g_mx_el%eigenvalue * convert1, " eV"


      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S11    ", g_mx_el%Sigma_Real_diag
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S12.re ", g_mx_el%Sigma_Real_offdiag
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S12.im ", g_mx_el%Sigma_Im_offdiag

      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L11    ", g_mx_el%L_Real_diag
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L12.re ", g_mx_el%L_Real_offdiag
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L12.im ", g_mx_el%L_Im_offdiag

!!$      Write (output_unit,&
!!$           Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
!!$           & g_mx_el%Sigma_Real_diag, &
!!$           & g_mx_el%Sigma_Im_diag, g_mx_el%L_Real_diag, &
!!$           & g_mx_el%L_Im_diag
!!$      Write (output_unit,&
!!$           & Fmt='(3F14.7,3x,3F14.7/3F14.7,3x,3F14.7)')&
!!$           & g_mx_el%Sigma_Real_offdiag, &
!!$           & g_mx_el%Sigma_Im_offdiag, g_mx_el%L_Real_offdiag, &
!!$           & g_mx_el%L_Im_offdiag

      g_Sigma (:, 1) = g_mx_el%Sigma_Real_offdiag(:)
      !g_Sigma (:, 2) = g_mx_el%Sigma_Im_offdiag(:)
      g_Sigma (:, 2) = -g_mx_el%Sigma_Im_offdiag(:)
      g_Sigma (:, 3) = g_mx_el%Sigma_Real_diag(:)
      ! g_Sigma (:, 3) = - g_mx_el%Sigma_Real_diag(:)
      !
      g_L (:, 1) = g_mx_el%L_Real_offdiag(:)
      g_L (:, 2) = -g_mx_el%L_Im_offdiag(:)
      !g_L (:, 2) = g_mx_el%L_Im_offdiag(:)
      g_L (:, 3) = g_mx_el%L_Real_diag(:)
      ! g_L (:, 3) = - g_mx_el%L_Real_diag(:)

      !
      g = 0.0_r8_kind !DG
      g (:, :) = ge * (g_Sigma(:, :)+g_L(:, :))
!!$      Write (output_unit,*) "g-tensor"
      Call print_G_tensors (g, 'S',"full")
      g = 0.0_r8_kind
      g (:, :) = ge * g_Sigma (:, :)
!!$      Write (output_unit,*) "g Sigma"
      Call print_G_tensors (g, 'S',"sigma")
      g = 0.0_r8_kind
      !
      g (:, :) = ge * g_L (:, :)
!!$      Write (output_unit,*) "Absolute values of g L"
      Call print_G_tensors (g, 'L',"L")
      !
    End Subroutine calc_G_tensor

    Subroutine calc_G_tensor_v2(vS11,vS12,vL11,vL12,&
         & irr,i_orb,eny)
      !
      real(r8_kind), intent(in)     :: vS11(3), vL11(3)
      complex(c16_kind), intent(in) :: vS12(3), vL12(3)
      integer(i4_kind), intent(in)  :: i_orb, irr
      real(r8_kind), intent(in)     :: eny
      ! *** end of interface ***

      Real (Kind=r8_kind), Dimension (3, 3) :: g_Sigma, g_L, g
      Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      integer :: unit
      !--------------------------------------------------------------------------
#if _DPRINT
      unit = 6
#else
      unit = output_unit
#endif
!!$      Write (unit,*) 'Kinematic factors is ', &
!!$           & options_kinematic_factors
      Write (unit, Fmt='(2(A,I3,2x),A,F14.7,A)')&
           & "irrep: ", irr,&
           & "orb: ",i_orb, &
           & " E  =", eny * convert1, " eV"

      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S11    ", vS11
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S12.re ", real(vS12)
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "S12.im ", aimag(vS12)

      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L11    ", vL11
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L12.re ", real(vL12)
      Write (unit, Fmt='(A,1X,T15,3F15.10)') "L12.im ", aimag(vL12)

      g_Sigma (:, 1) = real(vS12) ! g_mx_el%Sigma_Real_offdiag(:)
      ! FIXME: why I had to change sign?
!!$      g_Sigma (:, 2) = -aimag(vS12) ! -g_mx_el%Sigma_Im_offdiag(:)
      g_Sigma (:, 2) = aimag(vS12)
      g_Sigma (:, 3) = vS11 ! g_mx_el%Sigma_Real_diag(:)
      !
      g_L (:, 1) = real(vL12)  ! g_mx_el%L_Real_offdiag(:)
      ! FIXME: why I had to change sign?
!!$      g_L (:, 2) = -aimag(vL12) ! -g_mx_el%L_Im_offdiag(:)
      g_L (:, 2) = aimag(vL12)
      g_L (:, 3) = vL11 ! g_mx_el%L_Real_diag(:)
      !
      ! FIXME: why minus?:
      g = ge * (g_Sigma + g_L)
      Call print_G_tensors (g, 'S',"full")
      !
      g = ge * g_Sigma
      Call print_G_tensors (g, 'S',"sigma")
      !
      g = ge * g_L
      Call print_G_tensors (g, 'L',"L")
      !
    End Subroutine calc_G_tensor_v2

    Subroutine print_G_tensors (g,m_type,title,verbose)
      implicit none
      Real (Kind=r8_kind), Dimension (3, 3), Intent (In) :: g
      Character (Len=1), Intent (In) :: m_type
      character(len=*), intent(in),optional :: title
      integer(i4_kind), intent(in),optional :: verbose
      !-------------------------------------------------
      Real (Kind=r8_kind), Dimension (3, 3) :: G_tensor
      Real (Kind=r8_kind), Parameter :: &
           & diagonal_request_precision = 0.0000001
      Real (Kind=r8_kind), Dimension (3) :: diag
      Real (Kind=r8_kind) :: W (3), WORK (9)
      Logical  :: diagonalisation_request
      Integer  :: i, j, m, info
      Real(Kind=r8_kind) :: sm, buf
      integer :: unit,verbose_
      character(len=15) :: title_

#if _DPRINT
      unit = 6
#else
      unit = output_unit
#endif

      verbose_ = 2
      if(present(verbose)) verbose_ = verbose
      title_ = ""
      if(present(title)) title_ = title

      diagonalisation_request = .False.
      if(verbose_.ge.2)then
!!$         Write (unit,*) "pre g"
         Write (unit, Fmt='(A,1X,(T15,3F15.10))')&
              & trim(title_)//" pre g",&
              & g (1, :), g (2, :), g (3, :)
      endif
      !
      Do i = 1, 3
         Do j = 1, 3
            sm = 0.0
            Do m = 1, 3
               buf = g (i, m) * g (j, m)
               sm = sm + buf
            End Do
            G_tensor (i, j) = sm
            if (i /= j .And. Abs(sm) > diagonal_request_precision)then
               diagonalisation_request = .True.
            endif
         End Do
      End Do
      if(verbose_.ge.2)then
         Write (unit, Fmt='(A,1X,(T15,3F15.10))')&
              & trim(title_)//" big G",&
              & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
      endif
      !
      If (diagonalisation_request) Then
         if(verbose_.ge.2)then
            Write (unit,*) "G tensor is not diagonal"
            Write (unit, Fmt='(3F18.7/3F18.7/3F18.7)') &
                 & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
         endif
         Call DSYEV ('V', 'L', 3, G_tensor, 3, W, WORK, 9, info)
         If (info /= 0) Then
            Write (unit,*) "Could not diagonalize G-tensor&
                 &"
            Return
         End If
         if(verbose_.ge.2)then
            Write (unit,*) "Diagonal components of G"
            Write (unit, Fmt='(3F14.7)') W
            Write (unit,*) "Directions of the main axes of g &
                 &tensor (eigen vectors)"
            Write (unit, Fmt='(3F14.7/3F14.7/3F14.7)') &
                 & G_tensor (1, :), G_tensor (2, :), G_tensor (3, :)
            Write (unit,*) "Angles in degrees"
            Write (unit, Fmt='(3F14.7)') (Acos(G_tensor(1, &
              & 1))) * 180 / PI, (Acos(G_tensor(2, 2))) * 180 / PI, &
              & (Acos(G_tensor(3, 3))) * 180 / PI
            Write (unit,*) "g-tensor"
            Write (unit, Fmt='(3F14.7)') Sqrt (Abs(W(:)))
         endif
!!$         Write (unit,*) "Shift"
         Select Case (m_type)
         Case ('S')
            Write (unit, Fmt='(A,1X,T15,3F15.10)')&
                 & trim(title_)//" shift", Sqrt (Abs(W(:))) - ge
         Case ('L')
            Write (unit, Fmt='(A,1X,T15,3F15.10)')&
                 & trim(title_)//" shift", Sqrt (Abs(W(:)))
         End Select
         !
         Return
      Else
         Do i = 1, 3
            diag (i) = Sqrt (Abs(G_tensor(i, i)))
         End Do
         if(verbose_.ge.2)then
            Write (unit,*) "           xx                 yy &
                 &              zz"
            Write (unit, Fmt='(3F18.7)') diag (:)
            Write (unit,*) "Shift      xx                 yy &
                 &              zz"
         endif
         !
         Select Case (m_type)
         Case ('S')
            Write (unit, Fmt='(A,1X,T15,3F15.10)')&
                 & trim(title_)//" shift", diag (:) - ge
         Case ('L')
            Write (unit, Fmt='(A,1X,T15,3F15.10)')&
                 & trim(title_)//" shift", diag (:)
         End Select
         !
      End If
    End Subroutine print_G_tensors

    Logical Function exclude_contribution (bas1, bas2, i, &
         & pattern_1, pattern_2)
      Integer (Kind=i4_kind), Intent (In) :: bas1, bas2, i
      Logical, Dimension (:), Intent (In) :: pattern_1, pattern_2
      !
      !
      !!  logical :: include_contribution
      !--------------------------------
      If ( .Not. exclude_atom_orbitals) Then
         exclude_contribution = .False.
      Else
         If (i > 4) Then
            exclude_contribution = .False.
         Else
            !Only for L component
            exclude_contribution = (exclude_atom_orbitals .And. &
                 & ((pattern_1(bas1) .And. pattern_2(bas2)) .Or. &
                 & (pattern_2(bas1) .And. pattern_1(bas2)))) .Eqv. &
                 & include_key
         End If
      End If
    End Function exclude_contribution

    !
    Subroutine calc_write_tmgt (gti, n_orb1,i_irr,i_xyz, ev_real,ev_imag, tm_real,tm_imag )
      ! arguments:
      ! i_ip1 - irrep-partner
      ! i_ip2=1 for diagonal matrix elements
      ! i_ip2=2 for offdiagonal
      ! n_orb1 - number of orbitals (the size of matrix)
      ! i_irb - ...
      ! Porpose: calculation  of
      ! transition moments (matrix elements) between  orbital and the same orbital (diagonal)
      ! and Kramers cogugation (matrix offdiagonal)
      Use unique_atom_module
      Use dimensions, Only: BasDimSpor, SymDimSpor
      Integer (Kind=i4_kind), Intent (In) :: n_orb1, i_irr,i_xyz
      Real (Kind=r8_kind), Intent(Out) :: tm_real, tm_imag
      Real (Kind=r8_kind), Dimension (:), Pointer :: ev_real, ev_imag

      Type (dipole_integral_type), Pointer:: gti

      Integer (Kind=i4_kind) :: i_bas1, i_bas2, i_bas12,  &
           & i_lower, i_upper, i_ua, i_l, i_n, counter, n_f
      Real (Kind=r8_kind) :: coeff_real_a, &
           & coeff_imag_a, coef1, coef2

      Real (Kind=r8_kind), Parameter :: convert1 = &
           & 27.211652_r8_kind
      Logical :: pattern_1 (1:n_orb1), pattern_2 (1:n_orb1)
      Type (unique_atom_type), Pointer :: ua
      !
      !
      !
      tm_real = 0.0_r8_kind
      tm_imag = 0.0_r8_kind
      !
      pattern_1 (:) = .False.
      pattern_2 (:) = .False.
      If (exclude_atom_orbitals) Then
         i_lower = 1_i4_kind
         i_upper = 0
         !
         Do i_ua = 1, N_unique_atoms
            ua => unique_atoms (i_ua)
            Do i_l = 0, ua%lmax_ob
               n_f = BasDimSpor(i_ua)%LM(i_l, i_irr) / &
                    & SymDimSpor(i_ua)%LM(i_l, i_irr)! number of independent functions
               i_upper = i_upper + BasDimSpor(i_ua)%LM(i_l, &
                    & i_irr)
               !
               Do i_n = 1, n_f
                  !
                  If (i_ua == exclude_unique_atoms_1 .And. i_l == &
                       & exclude_atom_l_1 .And. i_n == &
                       & exclude_atom_n_1) Then
                     !
                     Do counter = 0, SymDimSpor(i_ua)%LM(i_l, &
                          & i_irr) - 1
                        pattern_1 (i_lower+(i_n-1)+counter*n_f) = &
                             & .True.
                     End Do
                  End If
               End Do
               i_lower = i_upper + 1
            End Do
         End Do
         i_lower = 1_i4_kind
         i_upper = 0
         !
         Do i_ua = 1, N_unique_atoms
            ua => unique_atoms (i_ua)
            Do i_l = 0, ua%lmax_ob
               n_f = BasDimSpor(i_ua)%LM(i_l, i_irr) / &
                    & SymDimSpor(i_ua)%LM(i_l, i_irr)! number of independent functions
               i_upper = i_upper + BasDimSpor(i_ua)%LM(i_l, &
                    & i_irr)
               !
               Do i_n = 1, n_f
                  !
                  If (i_ua == exclude_unique_atoms_2 .And. i_l == &
                       & exclude_atom_l_2 .And. i_n == &
                       & exclude_atom_n_2) Then
                     !
                     Do counter = 0, SymDimSpor(i_ua)%LM(i_l, &
                          & i_irr) - 1
                        pattern_2 (i_lower+(i_n-1)+counter*n_f) = &
                             & .True.
                     End Do
                  End If
               End Do
               i_lower = i_upper + 1
            End Do
         End Do
         !
      End if

      integralstorage_s: Select Case (gti%use)
      Case (DIPOLE_DIAGONAL) integralstorage_s
         !
         !triangular storage
         i_bas12 = 0

         bas1_diagonal_so_d: Do i_bas2 = 1, n_orb1
            !
            bas2_diagonal_so_d: Do i_bas1 = 1, i_bas2 - 1
               !
               i_bas12 = i_bas12 + 1
               !
               If (exclude_contribution(i_bas1, i_bas2, i_xyz, &
                    & pattern_1, pattern_2)) Cycle
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
                    & + gti%diagonal(i_bas12) * coef1 + &
                    & gti%diagonal_imag(i_bas12) * coef2
               !
               !
            End Do bas2_diagonal_so_d
            !
            i_bas12 = i_bas12 + 1
            !
            !
            If (exclude_contribution(i_bas1, i_bas2, i_xyz, &
                 & pattern_1, pattern_2)) Cycle
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
                 & + gti%diagonal(i_bas12) * coef1 + &
                 & gti%diagonal_imag(i_bas12) * coef2
            !
         End Do bas1_diagonal_so_d
         !
      Case (DIPOLE_OFFDIAGONAL) integralstorage_s
         ! full storage


         !
         Do i_bas2 = 1, n_orb1
            !
            Do i_bas1 = 1, n_orb1
               !
               If (exclude_contribution(i_bas1, i_bas2, i_xyz, &
                    & pattern_1, pattern_2)) Cycle
               !
               If (i_xyz == 4) Then
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
                  tm_real = &
                       & tm_real + &
                       & gti%offdiagonal(i_bas1, i_bas2) * coef1 + &
                       & gti%offdiagonal_imag(i_bas1, i_bas2) * &
                       & coef2
               Else
                  !
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
                       & gti%offdiagonal(i_bas1, i_bas2) * &
                       & coeff_real_a
                  tm_imag  = &
                       & tm_imag + &
                       & gti%offdiagonal(i_bas1, i_bas2) * &
                       & coeff_imag_a
                  !
                  tm_real = &
                       & tm_real - &
                       & gti%offdiagonal_imag(i_bas1, i_bas2) * &
                       & coeff_imag_a
                  !
                  tm_imag  = &
                       & tm_imag  + &
                       & gti%offdiagonal_imag(i_bas1, i_bas2) * &
                       & coeff_real_a
               End If
               !
            End Do
         End Do
         !
      Case (dipole_unused) integralstorage_s
         tm_real = 0.0_r8_kind
         tm_imag = 0.0_r8_kind
      Case Default integralstorage_s
         Call error_handler ("calc_write_tmgt: sxyz:  dipole_t&
              &ransitionmoment_f: forbidden case")
      End Select integralstorage_s

         End Subroutine calc_write_tmgt
!
!
      End Subroutine gtensor_calculate

    subroutine print_cmatrix(name,A)
      use matrix_module, only : cmatrix
      implicit none
      character(len=*), intent(in) :: name
      type(cmatrix), intent(in)    :: A
      ! *** end of interface ***

      integer(i4_kind) :: i, n,m, unit

      unit = output_unit

      write(unit,*) "Real part ",name,' shape=',shape(A%re)
      n = size(A%re, 1)
      m = size(A%re, 2)
      if(n>10) n = 10
      if(m>10) m = 10
      do i = 1, n
         write(unit,'(15F15.10)') A%re(i,:m)
      end do
      write(unit,*) "Imag part ",name,' shape=',shape(A%im)
      do i = 1, n
         write(unit,'(15F15.10)') A%im(i,:m)
      end do
    end subroutine print_cmatrix

    subroutine gti_map(i_ir, xyz, off, GTI, uc)
      !
      ! Fill GTI matrix with integrals from the integral storage
      !
      ! FIXME: I assume one was not going to modify the integrals in
      !        gten_integrals_* when one linked to the storage
      !        with pointers as in original code.
      !
      use matrix_module, only: cmatrix
      implicit none
      integer(i4_kind), intent(in) :: i_ir, xyz, off
      type(cmatrix), intent(out) :: GTI
      logical, optional, intent(in) :: uc
      ! *** end of interface ***

      Type (dipole_integral_type), Pointer :: hfc,hfc1
      logical :: uc_

      ABORT('add contr/uncontr')

      uc_ = .false.
      if(present(uc)) uc_ = uc

      ASSERT(xyz>0)
      ASSERT(xyz<=7)
      ASSERT(off>0)
      ASSERT(off<=2)

      DPRINT 'gti_map: irr=',i_ir,' xyz=',xyz,' kin. fac.=',options_kinematic_factors,' off=',off
      !if(options_kinematic_factors)then
      if(.not.options_kinematic_factors)then
         if(off.eq.2)then
            DPRINT 'gti_map: cont_offdiag'
               hfc => gten_integrals_cont_offdiag%integrals(xyz,i_ir,1)
         else
            DPRINT 'gti_map: cont_rem(1)'

            hfc => gten_integrals_cont_rem%integrals(xyz,i_ir,1)
         endif
      else
         if(off.eq.2)then
            DPRINT 'gti_map: offdiag'
            if(uc_)then
               hfc => gten_integrals_offdiag%integrals(xyz,i_ir,1)
            else
               hfc => gten_integrals_offdiag%integrals(xyz,i_ir,1)
            endif
         else
            hfc1 => gten_integrals_diag%integrals(xyz,i_ir ,1)
            hfc => gten_integrals_cont_rem%integrals(xyz,i_ir ,1)
            Call remap_to_full_store(hfc1,hfc)
            DPRINT 'gti_map: cont_rem(2)'
         endif
      endif

      !
      ! This is a derived type constructor, compare with the type declaration:
      !
      GTI = cmatrix(hfc%offdiagonal, hfc%offdiagonal_imag)

      ! Original code did not involve copying:
!     GTI%re => hfc%offdiagonal
!     GTI%im => hfc%offdiagonal_imag
    end subroutine gti_map

End Module gtensor_module

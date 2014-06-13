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
!=====================================================================================================#
! public interface of module
!=====================================================================================================#
module dft_d_module
!-----------------------------------------------------------------------
!
!  Purpose: This module performs the calculations necessary for the
!           evaluation of the empirical dispersion term.
!  Reference: S.Grimme;
!             J.Comput.Chem. 27 (15), 1787, 2006.
!
!           Input variables:
!           atom      : Nuclear charges / element identifier
!           pos       : Nuclear positions
!
!           Output variables:
!           Edisp     : energy value of the dispersion term
!           N_grads   : Nuclear gradients
!
!  Module called by:
!
!  Author: TS
!  Date: August, 2009
!
!== Interrupt of public interface of module ============================
!-----------------------------------------------------------------------
! Modifications
!-----------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author:      ...
! Date:        ...
! Description: ...
!
!------------ Modules used ---------------------------------------------

  use type_module

  implicit none
  private
  save

!------------ public functions and subroutines -------------------------
  public :: dft_d, E_dft_d, empirical_methods_read
  public :: empirical_methods_write
  public :: dft_d_applied

!=====================================================================================================#
! End of public interface of module
!=====================================================================================================#


  integer, parameter       :: RK=r8_kind
  integer, parameter       :: IK=i4_kind

  real(RK), parameter      :: ZERO     = 0.00_RK    ,&
                              Tenth    = 0.10_RK    ,&
                              ONE      = 1.00_RK    ,&
                              SIX      = 6.00_RK    ,&
                              THOUSAND = 1000.00_RK

  ! FIXME: Assign parameter automatically to used functional!
  real(RK), parameter      :: s6           = 0.75_RK             ! for PBE
                                                                 !
  real(RK), parameter      :: alpha        = 20.0_RK           & ! steepness of damping function
                            , lower_cutoff = 0.30_RK             ! Cutoff factor

  real(RK), parameter      :: Jmol_E_h     = 2625499.53918948_RK ! conversion factor E_h --> J/mol

  type                     :: gx4c
    integer(IK)            :: Z
    real(RK),dimension(3)  :: R
  end type gx4c

  real(RK)                 :: E_dft_d = ZERO

  ! variables and parameters for PG I/O
  logical, parameter       :: df_grimme06            = .FALSE. & !
                            , df_interaction_matrix  = .FALSE.

  logical                  :: grimme06                         & ! control variable for DFT-D
                            , interaction_matrix                 ! control variable for definition of interaction matrix

  logical, allocatable, dimension(:)     :: group_int_mat        ! interaction matrix for the different groups

  integer(IK), parameter   :: df_N_interaction_groups = 0        ! default number of interaction groups

  integer(IK)              :: N_interaction_groups               ! number of interaction groups

  integer(IK), allocatable, dimension(:) :: group_assignment     ! Number of unique atoms in group

  namelist /empirical_methods/ grimme06 &
                             , N_interaction_groups &
                             , interaction_matrix
  !===================================================================================================#
  contains
  !===================================================================================================#
  subroutine empirical_methods_read()
    !-------------------------------------------------------------------------------------------------+
    !                                                                                                 |
    ! Purpose: reads namelist empirical_methods from input                                            |
    !                                                                                                 |
    !-------------------------------------------------------------------------------------------------+
    use input_module
    use unique_atom_module , only : N_unique_atoms
    implicit none
    integer(IK)              :: stat, unit
    integer(IK)              :: i1, s_point, e_point
    !-------------------------------------------------------------------------------------------------+
    ! end of interface                                                                                |
    !-------------------------------------------------------------------------------------------------+
    !
    ! set defaults
    grimme06              = df_grimme06
    N_interaction_groups  = df_N_interaction_groups
    interaction_matrix    = df_interaction_matrix
    !
    ! read namelist
    if ( input_line_is_namelist("EMPIRICAL_METHODS") ) then
       if (N_unique_atoms .le. 0) call input_error("N_unique_atoms .LE. 0! Check if UNIQUE_ATOM_NUMBER &
                                                  & is set correctly and before EMPIRICAL_METHODS")
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=empirical_methods, iostat=stat)
       if (stat .gt. 0) call input_error("EMPIRICAL_METHODS_READ: namelist EMPIRICAL_METHODS")
       ! read interaction grouping when N_interaction_groups is set
       s_point = 1
       if (N_interaction_groups .GT. 1) then
         ! Allocate group interaction variables if input read for the first time
         if(.NOT. allocated(group_assignment)) allocate(group_assignment(N_unique_atoms))
         if(.NOT. allocated(group_int_mat))    allocate(group_int_mat(N_interaction_groups*(N_interaction_groups+1)/2))
         ! initialize group interaction matrix
         group_int_mat = .TRUE.
         call input_read_to_intermediate
         read(unit,fmt=*,iostat=stat) group_assignment
         ! check group_assignment
         if (minval(group_assignment) .LE. 0) then
           call error_handler("DFT-D interaction group assignment less than 1")
         else if (maxval(group_assignment) .GT. N_interaction_groups) then
           call error_handler("DFT-D interaction group assignment exceeds number of interaction groups")
         end if
         ! start with building/reading group interaction matrix
         if (stat .gt. 0) call input_error("EMPIRICAL_METHODS_READ: interaction group assignment")
         do i1=1,N_interaction_groups
           e_point = s_point+i1-1
           ! finish initialization of group interaction matrix as (1 - delta_ij)
           group_int_mat(e_point) = .FALSE.
           if (interaction_matrix) then
             ! read line of group interaction matrix
             call input_read_to_intermediate
             read(unit,fmt=*,iostat=stat) group_int_mat(s_point:e_point)
             if (stat .gt. 0) call input_error("EMPIRICAL_METHODS_READ: read group interaction matrix")
           end if
           s_point = s_point + i1
         end do
       end if
    endif
  end subroutine empirical_methods_read
  !===================================================================================================#
  !
  !===================================================================================================#
  subroutine empirical_methods_write(iounit)
    !-------------------------------------------------------------------------------------------------+
    !                                                                                                 |
    ! Purpose: writes input namelist empirical_methods                                                |
    !                                                                                                 |
    !-------------------------------------------------------------------------------------------------+
    use echo_input_module
    use unique_atom_module , only : N_unique_atoms
    use operations_module  , only : operations_echo_input_level
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer, intent(in)      :: iounit
    integer(IK)              :: stat
    integer(IK)              :: i1, i2
    !
    ! end of interface
    !-------------------------------------------------------------------------------------------------+
    ! Write Empirical Methods Namelist                                                                |
    !-------------------------------------------------------------------------------------------------+
    flag_format = '("    ",a," = ",a5:" # ",a)'
    call start("EMPIRICAL_METHODS","EMPIRICAL_METHODS_WRITE",iounit,operations_echo_input_level)
    call flag("GRIMME06          ", grimme06         , df_grimme06           )
    call intg("N_INTERACTION_GROUPS", N_interaction_groups , df_N_interaction_groups )
    call flag("INTERACTION_MATRIX", interaction_matrix, df_interaction_matrix )
    call stop(empty_line=.false.)
    !
    if (N_interaction_groups .GT. 1) then
      !-----------------------------------------------------------------------------------------------+
      ! Write Assignment of the Unique Atoms to the Interaction Groups                                |
      !-----------------------------------------------------------------------------------------------+
      write(iounit,fmt='("# Group assignment")',advance='YES',iostat=stat)
      write(iounit,fmt='("#")',advance='NO',iostat=stat)
      do i1 = 1,N_unique_atoms
        write(iounit,fmt='(I4)',advance='NO',iostat=stat) i1
        if (stat .gt. 0) call error_handler("EMPIRICAL_METHODS_WRITE: error before group assignment.")
      end do
      write(iounit,fmt='("   # Unique atom")',advance='YES',iostat=stat)
      write(iounit,fmt='(" ")',advance='NO',iostat=stat)
      do i1 = 1,N_unique_atoms
        write(iounit,fmt='(I4)',advance='NO',iostat=stat) group_assignment(i1)
        if (stat .gt. 0) call error_handler("EMPIRICAL_METHODS_WRITE: error at group assignment.")
      end do
      write(iounit,fmt='("   # Interaction Group")',advance='YES',iostat=stat)
      write(iounit,fmt='("#")',advance='YES',iostat=stat)
      !
      !-----------------------------------------------------------------------------------------------+
      ! Write Assignment of the Unique Atoms to the Interaction Groups                                |
      !-----------------------------------------------------------------------------------------------+
      write(iounit,fmt='("# Group Interaction Matrix")',advance='YES',iostat=stat)
      write(iounit,fmt='("#")',advance='NO',iostat=stat)
      do i1 = 1,N_interaction_groups
        write(iounit,fmt='(I4)',advance='NO',iostat=stat) i1
        if (stat .gt. 0) call error_handler("EMPIRICAL_METHODS_WRITE: error before interaction matrix.")
      end do
      write(iounit,fmt='("   # Interaction Group")',advance='YES',iostat=stat)
      do i1=1,N_interaction_groups
        write(iounit,fmt='(" ")',advance='NO',iostat=stat)
        do i2=1,i1
          write(iounit,fmt='(L4)',advance='NO',iostat=stat) group_int_mat(i1 * (i1-1) / 2 + i2)
        end do
        write(iounit,fmt='()',advance='YES',iostat=stat)
      end do
    endif
    write(iounit,fmt='()')
  end subroutine empirical_methods_write
  !===================================================================================================#
  !
  !===================================================================================================#
  subroutine dft_d(energy,gradients)
    use iounitadmin_module ,   only : openget_iounit,write_to_output_units,write_to_trace_unit
    use unique_atom_module ,   only : N_unique_atoms, unique_atoms
    use datatype ,             only : arrmat2
    !
    implicit none
    !-------------------------------------------------------------------------------------------------+
    !                                                                                                 |
    !  Purpose: wrapper for dft_d_grimme06. reads in gxfile and adds dispersion energy and gradients  |
    !           to e_sum and grads                                                                    |
    !                                                                                                 |
    !  Author: TS                                                                                     |
    !  Date: August, 2009                                                                             |
    !                                                                                                 |
    !-------------------------------------------------------------------------------------------------+
    !-------------------------------------------------------------------------------------------------+
    ! I/O variables                                                                                   |
    !------------------------------------------------------------------+                              |
    real(RK),optional,intent(inout)                   :: energy        !                              |
                                                                       !                              |
    type(arrmat2),optional,intent(inout),dimension(:) :: gradients     !                              |
                                                                       !                              |
    integer(IK)                                       :: atom_number  &!                              |
                                                       , i_atom       &!                              |
                                                       , i_uni        &!                              |
                                                       , i_equ        &!                              |
                                                       , i_grp        &!                              |
                                                       , i_sta        &!                              |
                                                       , i_end        &!                              |
                                                       , neqat         !                              |
                                                                       !                              |
    integer(IK), allocatable, dimension(:)            :: int_assign    !                              |
                                                                       !                              |
    logical, allocatable,dimension(:,:)               :: sq_int_mat    !                              |
                                                                       !                              |
    real(RK),allocatable,dimension(:,:),target        :: grds          !                              |
                                                                       !                              |
    type(gx4c),allocatable,dimension(:)               :: nuc           !                              |
                                                                       !                              |
    real(RK),dimension(:,:),pointer                   :: unigrds       !                              |
    !------------------------------------------------------------------+                              |
    ! Initialize variables                                                                            |
    !-------------------------------------------------------------------------------------------------+
    atom_number = 0
    i_atom      = 0
    !
    !-------------------------------------------------------------------------------------------------+
    ! Count total number of atoms and allocate variables                                              |
    !-------------------------------------------------------------------------------------------------+
    do i_uni = 1, N_unique_atoms
      atom_number = atom_number + unique_atoms(i_uni)%n_equal_atoms
    end do ! uniques
    !
    allocate(nuc(atom_number),grds(3,atom_number))
    !
    if (N_interaction_groups .GT. 1) then
      allocate(int_assign(atom_number),sq_int_mat(N_interaction_groups,N_interaction_groups))
    end if
    !
    !-------------------------------------------------------------------------------------------------+
    ! Unpack unique atom data - Z, coordinates, interaction group                                     |
    !-------------------------------------------------------------------------------------------------+
    do i_uni = 1, N_unique_atoms
      do i_equ = 1, unique_atoms(i_uni)%n_equal_atoms
        i_atom = i_atom + 1
        nuc(i_atom)%Z = unique_atoms(i_uni)%z
        nuc(i_atom)%R = unique_atoms(i_uni)%position(1:3,i_equ)
        if (N_interaction_groups .GT. 1) then
          !
          if (group_assignment(i_uni) .GE. 1) then
            int_assign(i_atom) = group_assignment(i_uni)
          end if
        end if
      end do ! equals
    end do ! uniques
    !
    !-------------------------------------------------------------------------------------------------+
    ! Reconstruct group interaction matrix                                                            |
    !-------------------------------------------------------------------------------------------------+
    if (N_interaction_groups .GT. 1) then
      do i_grp = 1, N_interaction_groups
        i_sta = i_grp * (i_grp - 1) / 2 + 1
        i_end = i_grp * (i_grp + 1) / 2
        sq_int_mat(1:i_grp,i_grp) = group_int_mat(i_sta:i_end)
        sq_int_mat(i_grp,1:i_grp-1) = group_int_mat(i_sta:i_end-1)
      end do
    end if
    !
    !-------------------------------------------------------------------------------------------------+
    ! Call dft_d module                                                                               |
    !-------------------------------------------------------------------------------------------------+
    if (N_interaction_groups .GT. 1) then
      call dft_d_grimme06(nuc,E_dft_d,grds,int_assign,sq_int_mat)
    else
      call dft_d_grimme06(nuc,E_dft_d,grds)
    end if
    !
    !-------------------------------------------------------------------------------------------------+
    ! Repack gradient correction and add to gradients                                                 |
    !-------------------------------------------------------------------------------------------------+
    if(present(gradients)) then
      i_atom = 0
      write(*,*) " DFT-D Gradient Contribution:"
      do i_uni = 1, N_unique_atoms
        neqat = unique_atoms(i_uni)%n_equal_atoms
        unigrds => grds(1:3,i_atom+1:i_atom+neqat)
        i_atom = i_atom + neqat
        gradients(i_uni)%m = gradients(i_uni)%m + unigrds
        write(*,*) 'Unique Center:',i_uni
        do i_equ = 1,neqat
          write(*,'(A20,3F20.12)') 'Equal Center:',unigrds(1:3,i_equ)
        end do
      end do ! uniques
      call write_to_trace_unit("DFT_D: gradients updated")
      call write_to_output_units("dftd_gradient: end")
    end if
    !
    !-------------------------------------------------------------------------------------------------+
    ! Add dispersion correction to total energy and print correction and updated energy               |
    !-------------------------------------------------------------------------------------------------+
    if(present(energy)) then
      energy = energy + E_dft_d
      call write_to_trace_unit("DFT_D: dispersion correction: ",real=E_dft_d)
      call write_to_trace_unit("DFT_D: corrected total energy:",real=energy)
      call write_to_output_units("DFT-D Correction")
      call write_to_output_units(" e_disp =           ",re=E_dft_d)
      call write_to_output_units(" e_sum  =           ",re=energy)
      call write_to_output_units(" ")
      call write_to_output_units(" ")
    end if
    !
    deallocate(nuc,grds)
    !
    if (N_interaction_groups .GT. 1) then
      deallocate(int_assign)
    end if
    !
  end subroutine dft_d
  !===================================================================================================#
  !
  !===================================================================================================#
  subroutine dft_d_grimme06(nuc,Edisp,N_grads,int_assign,sq_int_mat)
    use iounitadmin_module ,   only : write_to_output_units
    implicit none
    !-------------------------------------------------------------------------------------------------+
    !                                                                                                 |
    !  Purpose: Main module of DFT-D : Performs evaluation of the empirical dispersion terms          |
    !                                                                                                 |
    !  Reference : S. Grimme; J.Comput.Chem. 27 (2006); pages 1787-1799                               |
    !                                                                                                 |
    !  Schematics of procedure:                                                                       |
    !                                                                                                 |
    !  Author: TS                                                                                     |
    !  Date: August, 2009                                                                             |
    !                                                                                                 |
    !-------------------------------------------------------------------------------------------------+
    !-------------------------------------------------------------------------------------------------+
    ! I/O variables                                                                                   |
    !-------------------------------------------------------------+                                   |
    type(gx4c), dimension(:), intent(in)            :: nuc        ! nuclear data: charge, coords      |
                                                                  !                                   |
    integer(IK), dimension(:), optional, intent(in) :: int_assign !                                   |
                                                                  !                                   |
    logical, dimension(:,:), optional, intent(in)   :: sq_int_mat !                                   |
                                                                  !                                   |
    real(RK), intent(out)                           :: Edisp      ! Value of dispersion term          |
                                                                  !                                   |
    real(RK), dimension(:,:), intent(out)           :: N_grads    ! Nuclear gradients                 |
    !-------------------------------------------------------------+                                   |
    ! end of interface                                                                                |
    !-------------------------------------------------------------------------------------------------+
    !
    !-------------------------------------------------------------------------------------------------+
    ! auxiliary variables                                                                             |
    !-------------------------------------------------------+                                         |
    integer(IK)                             :: i1          &! indices                                 |
                                             , i2           !                                         |
                                                            !                                         |
    integer(IK)                             :: igroup       ! indices                                 |
                                                            !                                         |
    real(RK)                                :: Rij         &! distance ij                             |
                                             , Rijvdw      &! sum of vdW radii                        |
                                             , eij         &! exponential term ij                     |
                                             , Cij         &! coefficient ij                          |
                                             , Dij         &! dispersion contribution ij              |
                                             , D           &! unscaled total dispersion               |
                                             , deij_dRij   &! derivative of eij wrt Rij               |
                                             , dDij_dRij    !                                         |
                                                            !                                         |
    real(RK), dimension(3)                  :: dirij        ! direction i --> j (i1 --> i2)           |
                                                            !                                         |
    real(RK), dimension(size(nuc))          :: C6          &! vectors containing C6 and R0 parameters |
                                             , R0           ! for all atoms                           |
                                                            !                                         |
    real(RK), dimension(3,size(nuc))        :: grads        !                                         |
    !-------------------------------------------------------+                                         |
    ! executable lines                                                                                |
    !-------------------------------------------------------------------------------------------------+
    !
    !-------------------------------------------------------------------------------------------------+
    ! initializing variables                                                                          |
    !-------------------------------------------------------------------------------------------------+
    dirij   = ZERO
    Rij     = ZERO
    Rijvdw  = ZERO
    eij     = ZERO
    Cij     = ZERO
    Dij     = ZERO
    Edisp   = ZERO
    grads   = ZERO
    N_grads = ZERO
    D       = ZERO
    !
    !-------------------------------------------------------------------------------------------------+
    ! assigning parameters                                                                            |
    !-------------------------------------------------------------------------------------------------+
    if (present(int_assign) .AND. present(sq_int_mat)) then
      do i1 = 1,size(nuc)
        igroup = int_assign(i1)
        !
        !---------------------------------------------------------------------------------------------+
        ! Check if any interaction between actual group and groups                                    |
        !---------------------------------------------------------------------------------------------+
        if (ANY(sq_int_mat(:,igroup))) then
          call get_parameters(int(nuc(i1)%Z),C6(i1),R0(i1))
        else
          C6(i1) = ZERO
          R0(i1) = ZERO
        end if
        !
      end do
    else if (.NOT. (present(int_assign) .AND. present(sq_int_mat))) then
      do i1 = 1,size(nuc)
        call get_parameters(int(nuc(i1)%Z),C6(i1),R0(i1))
      end do
    end if
    !
    !-------------------------------------------------------------------------------------------------+
    ! loop over pairs of atoms                                                                        |
    !-------------------------------------------------------------------------------------------------+
    do i1 = 1,size(nuc)-1
      do i2 = i1+1,size(nuc)
        if (present(int_assign) .AND. present(sq_int_mat)) then
          if(.NOT. sq_int_mat(int_assign(i1),int_assign(i2))) then
            cycle ! to next interaction pair
          end if ! interaction considered
        end if ! present int_assign
        !---------------------------------------------------------------------------------------------+
        ! direction from atom i1 to atom i2                                                           |
        !---------------------------------------------------------------------------------------------+
        dirij = nuc(i2)%R - nuc(i1)%R
        !---------------------------------------------------------------------------------------------+
        ! distance between atom i1 and atom i2                                                        |
        !---------------------------------------------------------------------------------------------+
        Rij = sqrt(sum(dirij**2))
        !---------------------------------------------------------------------------------------------+
        ! sum of vdW radii of atom i1 and atom i2                                                     |
        !---------------------------------------------------------------------------------------------+
        Rijvdw = R0(i1) + R0(i2)
        !---------------------------------------------------------------------------------------------+
        ! Lower cutoff: function assumes very large values at small distances.                        |
        !               Dij and dDij_dRij set to zero if                                              |
        !                                                                                             |
        !                           Rij < Rij0 * 0.15                                                 |
        !                                                                                             |
        !               A warning is given out in such cases an nothing is incremented                |
        !---------------------------------------------------------------------------------------------+
        if (Rij .le. lower_cutoff * Rijvdw) then
          call write_to_output_units("DFT-D: WARNING!!! Rij under treshold value. Check gxfile.")
          cycle
        end if
        !---------------------------------------------------------------------------------------------+
        ! normalize direction                                                                         |
        !---------------------------------------------------------------------------------------------+
        dirij = dirij / Rij
        !---------------------------------------------------------------------------------------------+
        ! calculate exponential term               - d (Rij/Rij_0 - 1)                                |
        !                                   eij = e                                                   |
        !---------------------------------------------------------------------------------------------+
        eij = exp(- alpha * (Rij / Rijvdw - ONE))
        !---------------------------------------------------------------------------------------------+
        ! calculate parameter Cij                                                                     |
        !---------------------------------------------------------------------------------------------+
        Cij = sqrt(C6(i1)*C6(i2))
        !---------------------------------------------------------------------------------------------+
        ! calculate contribution Dij of pair i1,i2 to total dispersion                                |
        !                                                                                             |
        !                                           Cij       1                                       |
        !                                   Dij = ------- ---------                                   |
        !                                          Rij^6   1 + eij                                    |
        !---------------------------------------------------------------------------------------------+
        Dij = Cij / ((ONE + eij) * Rij**6)
        !---------------------------------------------------------------------------------------------+
        ! sum up to total dispersion D. Note the (-) sign                                             |
        !---------------------------------------------------------------------------------------------+
        D = D - Dij
        !---------------------------------------------------------------------------------------------+
        ! calculate derivative of eij                                - d (Rij/Rij_0 - 1)              |
        !                                 deij_dRij = - (d / Rij_0) e                                 |
        !                                                                                             |
        !---------------------------------------------------------------------------------------------+
        deij_dRij = alpha * eij / Rijvdw
        !---------------------------------------------------------------------------------------------+
        ! calculate derivative of Dij                                                                 |
        !                                                                                             |
        !                                                Cij        /   6      deij_dRij  \           |
        !                             dDij_dRij = ----------------- | ----- - ----------- |           |
        !                                          (1 + eij) Rij^6  \  Rij      1 + eij   /           |
        !                                                                                             |
        !---------------------------------------------------------------------------------------------+
        dDij_dRij = Dij * (SIX / Rij - deij_dRij / (ONE + eij))
        !---------------------------------------------------------------------------------------------+
        ! sum up gradients of atoms i1 and i2                                                         |
        !---------------------------------------------------------------------------------------------+
        grads(:,i1) = grads(:,i1) - dirij * dDij_dRij ! Points away from 2nd atom if dDij_dRij positive
        grads(:,i2) = grads(:,i2) + dirij * dDij_dRij ! Points away from 1st atom if dDij_dRij positive
      end do ! inner loop over pairs
    end do ! outer loop over pairs
    !
    !-------------------------------------------------------------------------------------------------+
    ! obtain total dispersion energy and gradients by scaling with s6 factor                          |
    !-------------------------------------------------------------------------------------------------+
    ! NOTE: changement needed, if other functional than PBE is used
    Edisp = s6 * D
    N_Grads = s6 * grads
    !
  end subroutine dft_d_grimme06
  !===================================================================================================#
  !
  !===================================================================================================#
  subroutine get_parameters(element,C6el,R0el)
    implicit none
    !-------------------------------------------------------------------------------------------------+
    !                                                                                                 |
    !  Purpose: Provides C6 and R0 parameters for the elements up to Xe                               |
    !                                                                                                 |
    !  Author: TS                                                                                     |
    !  Date: August, 2009                                                                             |
    !                                                                                                 |
    !-------------------------------------------------------------------------------------------------+
    !-------------------------------------------------------------------------------------------------+
    ! I/O variables                                                                                   |
    !-------------------------------------------------------+                                         |
    integer(IK), intent(in)                 :: element      !                                         |
    real(RK), intent(out)                   :: R0el       ,&! R0-parameter for inquired element       |
                                               C6el         ! C6-parameter for inquired element       |
    !-------------------------------------------------------+                                         |
    ! end of interface                                                                                |
    !-------------------------------------------------------------------------------------------------+

    !-------------------------------------------------------------------------------------------------+
    ! parameter vector: C6 and R0 in a.u.                                                             |
    !-------------------------------------------------------+                                         |
    real(RK), dimension(2,54), parameter :: CR = reshape((/&!  Element  C6    R0                      |
                                       0.14_RK , 1.001_RK ,&!H     2.42833898337275  1.89161585090411 |
                                       0.08_RK , 1.012_RK ,&!He    1.38762227621300  1.91240283827668 |
                                       1.61_RK , 0.825_RK ,&!Li   27.92589830878664  1.55902405294294 |
                                       1.61_RK , 1.408_RK ,&!Be   27.92589830878664  2.66073438368929 |
                                       3.13_RK , 1.485_RK ,&!B    54.29072155683364  2.80624329529730 |
                                       1.75_RK , 1.452_RK ,&!C    30.35423729215938  2.74388233317958 |
                                       1.23_RK , 1.397_RK ,&!N    21.33469249677488  2.63994739631672 |
                                       0.70_RK , 1.342_RK ,&!O    12.14169491686375  2.53601245945386 |
                                       0.75_RK , 1.287_RK ,&!F    13.00895883949688  2.43207752259099 |
                                       0.63_RK , 1.243_RK ,&!Ne   10.92752542517738  2.34892957310070 |
                                       5.71_RK , 1.144_RK ,&!Na   99.04153996470291  2.16184668674755 |
                                       5.71_RK , 1.364_RK ,&!Mg   99.04153996470291  2.57758643419900 |
                                      10.79_RK , 1.639_RK ,&!Al  187.15555450422841  3.09726111851332 |
                                       9.23_RK , 1.716_RK ,&!Si  160.09692011807493  3.24277003012132 |
                                       7.84_RK , 1.705_RK ,&!P   135.98698306887403  3.22198304274875 |
                                       5.57_RK , 1.683_RK ,&!S    96.61320098133015  3.18040906800361 |
                                       5.07_RK , 1.639_RK ,&!Cl   87.94056175499891  3.09726111851332 |
                                       4.61_RK , 1.595_RK ,&!Ar   79.96173366677415  3.01411316902302 |
                                      10.80_RK , 1.485_RK ,&!K   187.32900728875507  2.80624329529730 |
                                      10.80_RK , 1.474_RK ,&!Ca  187.32900728875507  2.78545630792473 |
                                      10.80_RK , 1.562_RK ,&!Sc  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Ti  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!V   187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Cr  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Mn  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Fe  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Co  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Ni  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Cu  187.32900728875507  2.95175220690531 |
                                      10.80_RK , 1.562_RK ,&!Zn  187.32900728875507  2.95175220690531 |
                                      16.99_RK , 1.650_RK ,&!Ga  294.69628091073594  3.11804810588589 |
                                      17.10_RK , 1.727_RK ,&!Ge  296.60426154052885  3.26355701749390 |
                                      16.37_RK , 1.760_RK ,&!As  283.94220827008525  3.32591797961161 |
                                      12.64_RK , 1.771_RK ,&!Se  219.24431964165407  3.34670496698419 |
                                      12.47_RK , 1.749_RK ,&!Br  216.29562230470145  3.30513099223904 |
                                      12.01_RK , 1.727_RK ,&!Kr  208.31679421647669  3.26355701749390 |
                                      24.67_RK , 1.628_RK ,&!Rb  427.90801942718406  3.07647413114074 |
                                      24.67_RK , 1.606_RK ,&!Sr  427.90801942718406  3.03490015639560 |
                                      24.67_RK , 1.639_RK ,&!Y   427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Zr  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Nb  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Mo  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Tc  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Ru  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Rh  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Pd  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Ag  427.90801942718406  3.09726111851332 |
                                      24.67_RK , 1.639_RK ,&!Cd  427.90801942718406  3.09726111851332 |
                                      37.32_RK , 1.672_RK ,&!In  647.32579185336465  3.15962208063103 |
                                      38.71_RK , 1.804_RK ,&!Sn  671.43572890256564  3.40906592910190 |
                                      38.44_RK , 1.881_RK ,&!Sb  666.75250372034668  3.55457484070991 |
                                      31.74_RK , 1.892_RK ,&!Te  550.53913808750792  3.57536182808248 |
                                      31.50_RK , 1.892_RK ,&! I  546.37627125886888  3.57536182808248 |
                                      29.99_RK , 1.881_RK  &!Xe  520.18490079534854  3.55457484070991 |
                                               /),(/2,54/)) !                                         |
    !-------------------------------------------------------+                                         |
    ! Conversion factors                                                                              |
    !-------------------------------------------------------+                                         |
    real(RK), parameter    :: C_con = 17.3452779073911_RK ,&! J nm^6 / mol --> E_h * a0^6             |
                              R_con = 1.88972612477933_RK   ! Angstroem -- > a0                       |
    !-------------------------------------------------------+                                         |
    ! executable lines                                                                                |
    !-------------------------------------------------------------------------------------------------+

    !-------------------------------------------------------------------------------------------------+
    ! assign parameters to elements with Z <= 54 otherwise stop and call error handler                |
    !-------------------------------------------------------------------------------------------------+
    if(element .gt. 0 .and. element .le. 54) then
      C6el = CR(1,element) * C_con
      R0el = CR(2,element) * R_con
    else
      call error_handler("input contains elements not supported by DFT-D")
    end if
  end subroutine get_parameters
  !===================================================================================================#
  !
  !===================================================================================================#
  logical function dft_d_applied()
    !-------------------------------------------------------------------------------------------------+
    ! Purpose: returns if empirical DFT-D dispersion correction is applied                            |
    !-------------------------------------------------------------------------------------------------+
    dft_d_applied = grimme06
  end function dft_d_applied
  !===================================================================================================#
end module DFT_D_module


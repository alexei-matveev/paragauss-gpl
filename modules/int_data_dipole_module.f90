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
module  int_data_dipole_module
  !---------------------------------------------------------------
  !
  !  Purpose: This is the module containing the data used when
  !           calculating the 2 center orbital and three center
  !           integrals to one specific quadrupel
  !           ( unique atom 1, l 1, unique atom 2, l 2 ).
  !           Which qudrupel is actually processed and the
  !           main loop variables and pointers to relevant data
  !           structures are to be contained here as well.
  !           Variables to hold primitive, contracted
  !           and symmetry adapted integarals are contained.
  !           Variables of primitive and contracted integrals
  !           only hold the values for one specific non-symmetry-
  !           equivalent connection vector.
  !
  !           Methods for allocating and freeing all necessary
  !           integals (int_data_dipole_setup and 
  !           int_data_dipole_shutdown) are included.
  !           int_data_dipole_setup should be called after
  !           quadrupel is set to present value.
  !           There is a reset method (int_data_dipole_reset)
  !           to initialize all required primitive and
  !           contracted integrals with zero.
  !
  !           To shorten list of use statements of calling
  !           routines and modules, this module is public.
  !
  !
  !  Module called by: ...
  !
  !
  !  References: Publisher Document: Concepts of Integral Part
  ! 
  !
  !  Author: TB
  !  Date: 5/96
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
  !== Interrupt end of public interface of module =================
  !------------ all modules used are public !! ------------------
#include <def.h>

  use type_module
  use datatype
  use quadrupel_module
  use unique_atom_module
  use symmetry_data_module
  use integralpar_module
  use fit_coeff_module
  use options_module, only: options_spin_orbit, options_kinematic_factors
  use iounitadmin_module

  !== Interrupt of public interface of module =====================
  implicit none
  save ! save all variables defined in this module
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ----------------------------

  type, public :: symadapt_nottotsym_2c_int_type
     ! to hold symmetry-adapted not total symmetric
     ! two center integrals of one IRREP
     real(kind=r8_kind), pointer :: int(:,:,:,:) => NULL()
     ! int(n_c2, n_c1, n_if2, n_if1)
     ! for relativistic integrals, n_c = n_exp
     real(kind=r8_kind), pointer :: int_imag(:,:,:,:) => NULL()
     ! int(n_c2, n_c1, n_if2, n_if1)
     ! imaginary part of integrals in case of spin orbit
  end type symadapt_nottotsym_2c_int_type

  !------------ Declaration of constants and variables ----------
  type(quadrupel_type),          public  :: quadrupel
  ! quadrupel ( unique atom 1, l 1, unique atom 2, l 2 )
  ! that is processed right now
  type(unique_atom_type), pointer,   public  :: ua1, ua2
  ! unique atoms processed right now
  type(unique_atom_basis_type), pointer, public :: ua1_basis, ua2_basis
  ! basis of unique atoms processed right now
  real(kind=r8_kind), pointer :: contractions1(:,:), contractions2(:,:)
  ! contractions(N_exponents:N_contracted_fcts)
  ! contraction matrices of unique atoms basis 1 and 2
  real(kind=r8_kind), dimension(3) :: center1, center2
  ! coordinates of centers
  integer(kind=i4_kind),             public  :: n_m1, n_m2
  ! number of mangnetical quantum numbers for l in quadrupel
  integer(kind=i4_kind),             public  :: n_exp1, n_exp2
  ! number of primitive exponents in ua1_basis and ua2_basis
  integer(kind=i4_kind),             public  :: n_contr1, n_contr2
  ! number of contractions in ua1_basis and ua2_basis
  integer(kind=i4_kind),             public  :: n_uncontr1, n_uncontr2
  ! number of uncontracted exponents used as trivial contractions
  ! in ua1_basis and ua2_basis
  integer(kind=i4_kind),             public  :: n_c1, n_c2
  ! number of contracted basis functions in ua1_basis and ua2_basis
  ! n_cx = n_contrx + n_uncontrx
  ! attention: this is not a dimension for symmetry_adapted 
  ! relativistic integrals !!
  logical,                           public :: diagonal
  ! true if (ua1 == ua2) .and. (l1 == l2) for quadrupel
  !---- Primitive integrals
  ! dimensions: ( n_exp2, n_exp1, n_m2, n_m1, 3 )
  ! last index: x, y, z component
  real(kind=r8_kind), pointer, public   :: prim_int_2cob_dipole(:,:,:,:,:)
  ! 2-center orbital integral of dipoles
  ! dimensions: ( n_exp2, n_exp1, n_m2, n_m1, 3 ) !:gten
  ! last index: x, y, z component
  real(kind=r8_kind), pointer, public   :: prim_int_2cob_dipoleg(:,:,:,:,:)
  ! 2-center orbital integral of L opetator
  real(kind=r8_kind), pointer, public   :: prim_int_2cob_hfc(:,:,:,:,:,:)
  ! 2-center orbital integral of hfc opetators

  !---- Contracted integrals
  ! dimensions:( n_c2, n_c1, n_m2, n_m1, 3 )
  ! last index: x, y, z component
  real(kind=r8_kind), pointer, public   :: cont_int_2cob_dipole(:,:,:,:,:)
  ! 2-center orbital integral of dipoles
  ! dimensions:( n_c2, n_c1, n_m2, n_m1, 3 )
  ! last index: x, y, z component
  real(kind=r8_kind), pointer, public   :: cont_int_2cob_dipoleg(:,:,:,:,:)
  ! 2-center orbital integral of the L operator
  real(kind=r8_kind), pointer, public   :: cont_int_2cob_hfc(:,:,:,:,:,:)
  ! 2-center orbital integral of the hfc operators
  !---- Symmetry-adapted integrals
  ! dimension: (n_ip,n_ip,3) combined index for irreps and partner used.
  ! See symmetry_data_module.
  ! Only part of the value range of inner index is really used:
  ! inner index runs from 1 to outer index if off-diagonal blocks required
  ! and is the same as outer otherwise.
  ! See symmetry_data_dipoles_exist(i_ir,i_ir,i_xyz).
  ! attention: not everything inside will be associated !!
  ! last index: x, y, z component
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_dipole(:,:,:)
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_dipole_p(:,:,:)
  ! 2-center orbital integral of kinetic energy
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_dipoleg(:,:,:)
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_dipoleg_p(:,:,:)
  ! 2-center orbital integral of  L operator
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_hfc(:,:,:,:)
  type(symadapt_nottotsym_2c_int_type), allocatable, target :: symadapt_int_2cob_hfc_p(:,:,:,:)
  ! 2-center orbital integral of  L operator

  !------------ public functions and subroutines ----------------
  public int_data_dipole_setup, int_data_dipole_shutdown

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Subroutines -------------------------------------
contains


  !**************************************************************
  subroutine int_data_dipole_setup()
    !  Purpose: setting of public variables and
    !           allocation of storare required all over the module
    !           (in particular symmetry-adapted integrals)
    !** End of interface ****************************************
    use operations_module, only: operations_gtensor, operations_dipole,operations_hfc
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status, i_xyz, n_ir, n_ip, i_ir1, i_ir2,diag_offdiag, &
         i_ip1, i_ip2, n_if1, n_if2, i_pa1, i_pa2, dim_c1, dim_c2, dim1, dim2, i_ua
    !------------ Executable code -------------------------------
    ! init public variables and error checks
    if ( quadrupel%ua1 .lt. 1 .or. quadrupel%ua1 .gt. N_unique_atoms) &
         call error_handler("int_data_dipole_setup: wrong number of unique atom 1")
    if ( quadrupel%ua2 .lt. 1 .or. quadrupel%ua2 .gt. N_unique_atoms) &
         call error_handler("int_data_dipole_setup: wrong number of unique atom 2")
    ua1 => unique_atoms(quadrupel%ua1)
    ua2 => unique_atoms(quadrupel%ua2)

    if ( quadrupel%l1 .lt. 0 .or. quadrupel%l1 .gt. ua1%lmax_ob) &
         call error_handler("int_data_dipole_setup: wrong l for unique atom 1")
    if ( quadrupel%l2 .lt. 0 .or. quadrupel%l2 .gt. ua2%lmax_ob) &
         call error_handler("int_data_dipole_setup: wrong l for unique atom 2")
    ua1_basis => ua1%l_ob(quadrupel%l1)
    ua2_basis => ua2%l_ob(quadrupel%l2)
    contractions1 => ua1_basis%contractions
    contractions2 => ua2_basis%contractions

    n_m1 = ( 2 * quadrupel%l1 ) + 1
    n_m2 = ( 2 * quadrupel%l2 ) + 1
    n_exp1 = ua1_basis%N_exponents
    n_exp2 = ua2_basis%N_exponents
    n_contr1 = ua1_basis%N_contracted_fcts
    n_contr2 = ua2_basis%N_contracted_fcts
    n_uncontr1 = ua1_basis%N_uncontracted_fcts
    n_uncontr2 = ua2_basis%N_uncontracted_fcts
    n_c1 = n_contr1 + n_uncontr1
    n_c2 = n_contr2 + n_uncontr2
    diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
         (quadrupel%l1 == quadrupel%l2)

    ! local variables required for allocation
    if ( integralpar_relativistic) then
       dim_c1 = n_exp1
       dim_c2 = n_exp2
    else
       dim_c1 = n_c1
       dim_c2 = n_c2
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_ip = symmetry_data_n_ip_proj()
       n_ir = symmetry_data_n_proj_irreps() 
    else
       n_ip = symmetry_data_n_ip()
       n_ir = symmetry_data_n_irreps()
    endif


    ! allocate symmetry adapted integrals
    if ( integralpar_2cob_dipole ) then
       notso:  if(.not.options_spin_orbit) then
       opdip:  if(operations_dipole) then
          allocate(symadapt_int_2cob_dipole(n_ip,n_ip,3), stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_data_dipole_setup: allocate of symadapt_int_2cob_dipole failed" )
          xyz: do i_xyz = 1, 3
             i_ip1 = 1
             irrep1: do i_ir1 = 1, n_ir
                n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                   i_ip2 = 1
                   irrep2: do i_ir2 = 1, n_ir
                      n_if2 = ua2%symadapt_partner&
                           (i_ir2,quadrupel%l2)%N_independent_fcts
                      partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                         if ( (i_ip1 .eq. i_ip2) .or. &
                              (integralpar_offdiag_dipoles &
                              .and. (i_ip2 .lt. i_ip1) .and. &
                              symmetry_data_dipoles_exist&
                              (i_ir2,i_ir1,i_xyz)) )then
                            allocate( symadapt_int_2cob_dipole&
                                 (i_ip2,i_ip1,i_xyz)%int&
                                 (dim_c2,dim_c1,n_if2,n_if1), stat=status )
                            if ( status .ne. 0 ) call error_handler( & 
                                 "int_data_dipole_setup: allocate of symadapt_int_2cob_dipole%int failed" )
                            symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int = 0.0_r8_kind
                            
                         else
                            nullify(symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int)
                         endif
                         i_ip2 = i_ip2 + 1
                      enddo partner2
                   enddo irrep2
                   i_ip1 = i_ip1 + 1
                enddo partner1
             enddo irrep1
          enddo xyz
       end if opdip

       if(operations_hfc) then
          !Not spin-orbit case
             allocate( symadapt_int_2cob_hfc_p(1,n_ip,10,n_unique_atoms), stat=status ) 
             ! d2(1/r)/dxdx, xydy,dzdz,dxy, dxz, dyz  = 6 + Lk/r^3 = 9; 10 -for isotropic part
             ! 1 - diagonal
             if ( status .ne. 0 ) call error_handler( &   
                  "int_data_dipole_setup: allocate of symadapt_int_2cob_dipoleg failed" )
             do i_ua = 1, n_unique_atoms
             do i_xyz = 1, 10
                i_ip1 = 1
                do i_ir1 = 1, n_ir
                   do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                      do diag_offdiag = 1,1
                         nullify( symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz,i_ua)%int)
                         
                      end do
                      i_ip1 = i_ip1 + 1
                   enddo
                enddo
             enddo
          end do
             if(options_kinematic_factors) then
                dim1 = n_exp1
                dim2 = n_exp2
             else
                dim1 = dim_c1
                dim2 = dim_c2
             end if
              do  i_ua = 1, n_unique_atoms
             xyz_2sohfc: do i_xyz = 1, 10
                do diag_offdiag = 1,1
                   i_ip1 = 1
                   do i_ir1 = 1, n_ir

                      n_if1 = ua1%symadapt_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                       do i_pa1 = 1, symmetry_data_n_partners(i_ir1)

                         i_ir2 = i_ir1
                         i_pa2 = i_pa1

                         n_if2 = ua2%symadapt_partner(i_ir2,quadrupel%l2)%N_independent_fcts

                         allocate( symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz, i_ua)%&
                              &int(dim2,dim1,n_if2,n_if1), &
                              stat=status ) 

                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_hfc_setup: allocate of symadapt_int_2cob_hfc%int failed" )
                         symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int &
                              = 0.0_r8_kind

                         i_ip1 = i_ip1 + 1
                      enddo 
                   enddo
                end do 
             enddo xyz_2sohfc
          end do 

       end if
       
    end if notso

       so: if (options_spin_orbit) then
         
          n_ir = symmetry_data_n_proj_irreps()
          n_ip = symmetry_data_n_ip_proj()
          op_dipole:if (operations_dipole) then
          allocate( symadapt_int_2cob_dipole_p(n_ip,n_ip,3), stat=status )
          if(status.ne.0) call error_handler('symadapt_int_2cob_dipole_p allocation failed')

          xyz_2: do i_xyz = 1, 3
             i_ip1 = 1
             irrep1_2: do i_ir1 = 1, n_ir
                n_if1 = ua1%symadapt_spor_partner&
                     (i_ir1,quadrupel%l1)%N_independent_fcts
                partner1_2: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                   i_ip2 = 1
                   irrep2_2: do i_ir2 = 1, n_ir
                      n_if2 = ua2%symadapt_spor_partner&
                           (i_ir2,quadrupel%l2)%N_independent_fcts
                      partner2_2: do i_pa2 = 1, symmetry_data_n_partners_proj(i_ir2)


                         if ( (i_ip1 .eq. i_ip2) .or. &
                              (integralpar_offdiag_dipoles .and. (i_ip2 .lt. i_ip1) .and. &
                              symmetry_data_dipoles_exist(i_ir2,i_ir1,i_xyz)) )then

                            allocate&
                                 ( symadapt_int_2cob_dipole_p&
                                 (i_ip2,i_ip1,i_xyz)%&
                                 &int(dim_c2,dim_c1,n_if2,n_if1), &
                                 stat=status )
                            if ( status .ne. 0 ) call error_handler( & 
                                 "int_data_dipole_setup: allocate of symadapt_int_2cob_dipole%int failed" )
                            symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int = 0.0_r8_kind

                            allocate( symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)&
                                 &%int_imag(dim_c2,dim_c1,n_if2,n_if1), &
                                 stat=status )
                            if ( status .ne. 0 ) call error_handler( & 
                                 "int_data_dipole_setup: allocate of symadapt_int_2cob_dipole%int failed" )
                            symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int_imag = 0.0_r8_kind
                        
                         endif
                         i_ip2 = i_ip2 + 1
                      enddo partner2_2
                   enddo irrep2_2
                   i_ip1 = i_ip1 + 1
                enddo partner1_2
             enddo irrep1_2
          enddo xyz_2
       end if op_dipole

       op_gten: if(operations_gtensor) then
             allocate( symadapt_int_2cob_dipoleg_p(2,n_ip,10), stat=status ) !DG 7+3 for gauge correction term 
             if ( status .ne. 0 ) call error_handler( &                      ! 1 - diagonal, 2 -for constracted Kramers pair
                  "int_data_dipole_setup: allocate of symadapt_int_2cob_dipoleg failed" )

             do i_xyz = 1, 10!DG  7+3 !3->4
                i_ip1 = 1
                do i_ir1 = 1, n_ir
                   do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                      do diag_offdiag = 1,2
                         nullify( symadapt_int_2cob_dipoleg_p&
                              (diag_offdiag,i_ip1,i_xyz)%int)
                         nullify(symadapt_int_2cob_dipoleg_p&
                              (diag_offdiag,i_ip1,i_xyz)%int_imag)  !for diagonal and offdiagonal
                      end do
                      i_ip1 = i_ip1 + 1
                   enddo
                enddo
             enddo

             if(options_kinematic_factors) then
                dim1 = n_exp1
                dim2 = n_exp2
             else
                dim1 = dim_c1
                dim2 = dim_c2
             end if
             xyz_2sog: do i_xyz = 1, 10 !DG 10=7+3 !3->4
                diag_offdiag_alloc:   do diag_offdiag = 1,2
                   i_ip1 = 1
                   irrep1_2g: do i_ir1 = 1, n_ir

                      n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                      partner1_2g: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)

                         i_ir2 = i_ir1
                         i_pa2 = i_pa1

                         n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts

                         allocate( symadapt_int_2cob_dipoleg_p&
                              (diag_offdiag,i_ip1,i_xyz)%&
                              &int(dim2,dim1,n_if2,n_if1), &
                              stat=status ) 

                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_dipole_setup: allocate of symadapt_int_2cob_dipoleg%int failed" )
                         symadapt_int_2cob_dipoleg_p(diag_offdiag,i_ip1,i_xyz)%int &
                              = 0.0_r8_kind

                         allocate( symadapt_int_2cob_dipoleg_p&
                              (diag_offdiag,i_ip1,i_xyz)&
                              &%int_imag(dim2,dim1,n_if2,n_if1), &
                              stat=status )

                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_dipoleg_setup: allocate of symadapt_int_2cob_dipoleg%int_imag failed" )
                         symadapt_int_2cob_dipoleg_p&
                              (diag_offdiag,i_ip1,i_xyz)%int_imag &
                              = 0.0_r8_kind

                         i_ip1 = i_ip1 + 1
                      enddo partner1_2g
                   enddo irrep1_2g
                end do diag_offdiag_alloc
             enddo xyz_2sog
          endif op_gten

       op_hfc: if(operations_hfc) then
             allocate( symadapt_int_2cob_hfc_p(2,n_ip,10,n_unique_atoms), stat=status ) 
             ! hyperfine_dipole_xyz + L_xyz/r^3 + f(r) * sigma_xyz(isotropic part) = 9
             ! 1 - diagonal, 2 -for constracted Kramers pair
             if ( status .ne. 0 ) call error_handler( &   
                  "int_data_dipole_setup: allocate of symadapt_int_2cob_dipoleg failed" )
             do i_ua = 1, n_unique_atoms
             do i_xyz = 1, 9
                i_ip1 = 1
                do i_ir1 = 1, n_ir
                   do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                      do diag_offdiag = 1,2
                         nullify( symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz,i_ua)%int)
                         nullify(symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz,i_ua)%int_imag)  !for diagonal and offdiagonal
                      end do
                      i_ip1 = i_ip1 + 1
                   enddo
                enddo
             enddo
          end do
             if(options_kinematic_factors) then
                dim1 = n_exp1
                dim2 = n_exp2
             else
                dim1 = dim_c1
                dim2 = dim_c2
             end if
             ua: do  i_ua = 1, n_unique_atoms
                 do i_xyz = 1, 10
                diag_offdiag_alloc_hfc:   do diag_offdiag = 1,2
                   i_ip1 = 1
                   irrep1_2_hfc: do i_ir1 = 1, n_ir

                      n_if1 = ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)%N_independent_fcts
                      partner1_2_hfc: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)

                         i_ir2 = i_ir1
                         i_pa2 = i_pa1

                         n_if2 = ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)%N_independent_fcts

                         allocate( symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz, i_ua)%&
                              &int(dim2,dim1,n_if2,n_if1), &
                              stat=status ) 

                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_hfc_setup: allocate of symadapt_int_2cob_hfc%int failed" )
                         symadapt_int_2cob_hfc_p(diag_offdiag,i_ip1,i_xyz,i_ua)%int &
                              = 0.0_r8_kind

                         allocate( symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz, i_ua)&
                              &%int_imag(dim2,dim1,n_if2,n_if1), &
                              stat=status )

                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_dipoleg_setup: allocate of symadapt_int_2cob_dipoleg%int_imag failed" )
                         symadapt_int_2cob_hfc_p&
                              (diag_offdiag,i_ip1,i_xyz, i_ua)%int_imag &
                              = 0.0_r8_kind

                         i_ip1 = i_ip1 + 1
                      enddo partner1_2_hfc
                   enddo irrep1_2_hfc
                end do diag_offdiag_alloc_hfc
             enddo
          end do ua
       endif op_hfc
 
       endif so! options_spin_orbit
    endif 
  end subroutine int_data_dipole_setup

  subroutine int_data_dipole_shutdown()
    use operations_module, only: operations_hfc
    !  Purpose: shutdown work + deallocating
    !           (in particular symmetry-adapted integrals)
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status, i_ip1, i_ip2, i_xyz, n_ip, i_ua
    !------------ Executable code -------------------------------
    ! deallocate symmetry adapted integrals

    n_ip = symmetry_data_n_ip()

    if ( integralpar_2cob_dipole ) then

       if(.not.options_spin_orbit) then
          if ( allocated(symadapt_int_2cob_dipole) ) then
          n_ip = symmetry_data_n_ip()
          do i_xyz = 1, 3
             do i_ip1 = 1, n_ip
                do i_ip2 = 1, i_ip1
                   if ( associated(symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int) ) then
                      deallocate( symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int, stat=status )
                      if ( status .ne. 0 ) call error_handler( & 
                           "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipole%int single failed" )
                   end if

                enddo
             enddo
          enddo
          deallocate( symadapt_int_2cob_dipole, stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipole failed" )
       end if
       
       if (operations_hfc) then 
            if ( allocated(symadapt_int_2cob_hfc_p) ) then
               do i_ua = 1, n_unique_atoms
                  do i_xyz = 1, 9
                     do i_ip1 = 1, n_ip
                        do i_ip2 = 1,1
                           if( associated(&
                                symadapt_int_2cob_hfc_p(i_ip2,i_ip1,i_xyz,i_ua)%int) ) then
                              deallocate( &
                                   symadapt_int_2cob_hfc_p(i_ip2,i_ip1,i_xyz,i_ua)%int, stat=status )
                              if ( status .ne. 0 ) call error_handler( & 
                                   "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p%int failed" )
                           end if
                        end do
                     enddo
                  enddo
               end do
               deallocate( symadapt_int_2cob_hfc_p, stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p failed" )
            end if
            
         end if
      end if

       SO:if (options_spin_orbit) then
          n_ip = symmetry_data_n_ip_proj()
      opdip:if( allocated(symadapt_int_2cob_dipole_p) ) then
             do i_xyz = 1, 3
                do i_ip1 = 1, n_ip
                   do i_ip2 = 1, i_ip1
                      if ( associated(symadapt_int_2cob_dipole_p&
                           (i_ip2,i_ip1,i_xyz)%int) ) then
                         deallocate( symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int, stat=status )
                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipole%int failed" )
                         deallocate( symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int_imag, stat=status )
                         if ( status .ne. 0 ) call error_handler( & 
                              "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipole%int_imag failed" )
                      endif
                   enddo
                enddo
             enddo
             deallocate( symadapt_int_2cob_dipole_p, stat=status )
         if ( status .ne. 0 ) call error_handler( &
            "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipole_p failed" )
      end if opdip

     opgten:if ( allocated(symadapt_int_2cob_dipoleg_p) ) then
          xyzs:do i_xyz = 1, 10!DG 7+3
             do i_ip1 = 1, n_ip
                do i_ip2 = 1,2
                   if( associated(&
                        symadapt_int_2cob_dipoleg_p(i_ip2,i_ip1,i_xyz)%int) ) then
                      deallocate( &
                           symadapt_int_2cob_dipoleg_p(i_ip2,i_ip1,i_xyz)%int, stat=status )
                      if ( status .ne. 0 ) call error_handler( & 
                           "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p%int failed" )
                      deallocate( &
                           symadapt_int_2cob_dipoleg_p(i_ip2,i_ip1,i_xyz)%int_imag, stat=status )
                      if ( status .ne. 0 ) call error_handler( & 
                           "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p%int_imag failed" )
                   endif
                end do
             enddo
          enddo xyzs

          deallocate( symadapt_int_2cob_dipoleg_p, stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p failed" )
       end if opgten
       
  ophfc:if ( allocated(symadapt_int_2cob_hfc_p) ) then
   uahfc:  do i_ua = 1, n_unique_atoms
          xyzshfc:do i_xyz = 1, 10
             do i_ip1 = 1, n_ip
                do i_ip2 = 1,2
                   if( associated(&
                        symadapt_int_2cob_hfc_p(i_ip2,i_ip1,i_xyz,i_ua)%int) ) then
                      deallocate( &
                           symadapt_int_2cob_hfc_p(i_ip2,i_ip1,i_xyz,i_ua)%int, stat=status )
                      if ( status .ne. 0 ) call error_handler( & 
                           "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p%int failed" )
                      deallocate( &
                           symadapt_int_2cob_hfc_p(i_ip2,i_ip1,i_xyz,i_ua)%int_imag, stat=status )
                      if ( status .ne. 0 ) call error_handler( & 
                           "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p%int_imag failed" )
                   endif
                end do
             enddo
          enddo xyzshfc
       end do uahfc
          deallocate( symadapt_int_2cob_hfc_p, stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "int_data_dipole_shutdown: deallocate of symadapt_int_2cob_dipoleg_p failed" )
       end if ophfc

    end if SO
 endif

end subroutine int_data_dipole_shutdown
!**************************************************************
end module int_data_dipole_module

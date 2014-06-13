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
module  pseudo_potential_module
!---------------------------------------------------------------
!
!   Purpose: information about the basis functions of
!            pseudo potentials, which are used for
!            integral part
!
!  modules called by: data used in integral
!
!  Author: HU
!  Data:   10/98
!
!
!---------------------------------------------------------------
!== Interrupt of public interface of module ====================
!---------------------------------------------------------------
! Modification
!---------------------------------------------------------------
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

use type_module  ! type specification parameters
use comm_module   ! pvm related variables and routines
use unique_atom_module, only : unique_atoms
implicit none
save             ! save all variables defined in this module
private          ! by defaut, all names are private
!== Interrupt end of public interface of module ==================


!--------- Declaration of types ----------------------------------
type, public :: unique_pseudo_atoms_basis_type
   ! angular quantum number dependent information for either local
   ! potential ( type one ) and semilocal potential ( type two )
   integer(kind=i4_kind)                 :: n_terms
           ! number of the terms for Gaussians
   integer(kind=i4_kind), pointer        :: n_or_sym(:)
           ! symmetry of corresponding Gaussians
   real(kind=r8_kind), pointer           :: cexponents(:)
           ! Gaussian exponentials
   real(kind=r8_kind), pointer           :: dcoefficients(:)
           ! combination coefficients of Gaussians
end type unique_pseudo_atoms_basis_type

type, public :: unique_pseudo_atoms_type
   ! describes all properties of primitive, contracted basis functions
   ! of an unique atom with pseudopotentials
   character(len=12)                     :: name
            ! name of unique pseudo_atom
   real(kind=r8_kind)               :: zc
            ! nuclear charge for the core
   integer(kind=i4_kind)            :: n_equal_pseudo_atoms
            ! number of partners of unique pseudo atoms
   real(kind=r8_kind), pointer       :: position_pseudo(:,:)
            ! coordinates of all centers
   real(kind=r8_kind)                :: position_first_pseudo(3)
            ! positions of first partner of unique pseudo atom
   integer(kind=i4_kind)             :: l_max1 ! l_max + 1
       ! angular momentum plus one
   type(unique_pseudo_atoms_basis_type), pointer :: l_ob_pseudo(:) ! l_or_pseudo(0:l_max1-1)
       ! information of pseudopotential basis set for semilocal case
 end type unique_pseudo_atoms_type


   !------------------ Declaration of public variables ------------------------------
   integer(kind=i4_kind), public                   :: N_unique_pseudo_atoms
          ! Number of unique atoms
   type(unique_pseudo_atoms_type), pointer, public :: unique_pseudo_atoms(:)

   !----------------- public functions and subroutine -------------------------------
   public unique_pseudo_atom_read_basis,   unique_pseudo_atom_basis_read,    &
          unique_pseudo_atom_pack,         unique_pseudo_atom_unpack,        &
          unique_pseudo_atom_read,         unique_pseudo_atom_write,         &
          unique_pseudo_atom_alloc,        unique_pseudo_atom_basis_alloc,   &
          unique_pseudo_atom_basis_pack,   unique_pseudo_unique_atom_write,  &
          unique_pseudo_unique_atom_read,  unique_pseudo_atom_basis_unpack

   !=================================================================================
   ! End of public interface of module
   !=================================================================================

   !---------------- declaration of constants and variables -------------------------
   ! namelists for reading and writing

   namelist /unique_pseudo_atom_number/ N_unique_pseudo_atoms
   character(len=12), private                      :: name
   real(kind=r8_kind), private                :: zc
   integer(kind=i4_kind), private             :: N_equal_pseudo_atoms
   namelist /unique_pseudo_atom/zc, N_equal_pseudo_atoms, name

   integer(kind=i4_kind), private             :: l_max1

   namelist /unique_pseudo_atom_basisset/l_max1

   integer(kind=i4_kind), private             :: n_terms

   namelist /unique_pseudo_basis/n_terms

   ! defaults for namelist input
   integer(kind=i4_kind), private :: df_N_unique_pseudo_atoms  =  0
   character(len=12)         , private :: df_name                   =  "            "
   real(kind=r8_kind)   , private :: df_zc                     =  0.0_r8_kind
   integer(kind=i4_kind), private :: df_N_equal_pseudo_atoms   =  0
   integer(kind=i4_kind), private :: df_l_max1                 =  0
   integer(kind=i4_kind), private :: df_n_terms                =  0

   !------------- Initialising Variables -----------------------------------------
   data N_unique_pseudo_atoms,N_equal_pseudo_atoms,l_max1,n_terms,&
        name,zc/4*0," ",0.0E0_r8_kind/

   !------------------------------------------------------------------------------

   !---------------- Subroutines -------------------------------------------------

   contains

   !*******************************************************************************
   subroutine unique_pseudo_atom_alloc(ua_pseudo,l_max1,n_equal_pseudo_atoms)
   !  Purpose: allocates arrays contained in ua and initialises ua
   !           with the exception of symadapt and renorm
   !           ua%lmax_all must be defined before first call
   !------------ Declaration of formal parameters ----------------------------------
   type(unique_pseudo_atoms_type), intent(inout) :: ua_pseudo
   integer, optional,      intent(in)    :: l_max1
   integer, optional,      intent(in)    :: n_equal_pseudo_atoms
   !** End of interface ***********************************************************
   !------------ Declaration of local variables -----------------------------------
   integer(kind=i4_kind)                 :: status
   !------------ Declaration of subroutines used ----------------------------------
   external error_handler
   !------------ Executable code --------------------------------------------------
   if ( present(l_max1) ) then
        ua_pseudo%l_max1 = l_max1
      allocate( ua_pseudo%l_ob_pseudo(0:l_max1), stat=status )
   if(status.ne.0)call error_handler(&
     "unique_pseudo_atom_alloc: allocate of l_ob_pseudo failed")
   endif

   if( present(n_equal_pseudo_atoms) ) then
      ua_pseudo%n_equal_pseudo_atoms = n_equal_pseudo_atoms
      allocate(ua_pseudo%position_pseudo(3,n_equal_pseudo_atoms),stat=&
               status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_alloc: allocate of position_pseudo failed")
   end if
   end subroutine unique_pseudo_atom_alloc
   !=============================================================================

   subroutine unique_pseudo_atom_read()
   ! purpose: reads in names and zc of unique pseudo atoms
   !          and coordinates of all centers ( " equal pseudo atoms " )
   ! called by read_input()
   !---- end of interface ---------------------------------------------------------

   !---- used modules -------------------------------------------------------------
   use input_module
   implicit none
   !---- declaration of local variables -------------------------------------------
   type(unique_pseudo_atoms_type), pointer   :: ua_pseudo
   integer(kind=i4_kind)                     :: i,j,k,status,unit
   !---- declaration of subroutines used ------------------------------------------
   external error_handler
   !---- Excutable codes ----------------------------------------------------------

   if ( .not. input_line_is_namelist("unique_pseudo_atom_number") )&
           call input_error( &
    "unique_pseudo_atom_read: namelist unique_pseudo_atom_number expected")

   unit = input_intermediate_unit()

   n_unique_pseudo_atoms = df_n_unique_pseudo_atoms
   if ( input_line_is_namelist("unique_pseudo_atom_number") ) then
      call input_read_to_intermediate()
      read(unit, nml=unique_pseudo_atom_number, iostat=status)
      if (status .ne. 0) call input_error( &
           "unique_pseudo_atom_read: namelist unique_pseudo_atom_number.")
   endif
   if(n_unique_pseudo_atoms.lt.0)call input_error( &
   "unique_pseudo_atom_read: namelist unique_pseudo_atom_number.")

   allocate( unique_pseudo_atoms(n_unique_pseudo_atoms), stat=status )
   if (status .ne. 0) call error_handler( &
        "unique_pseudo_atom_read: allocate failed.")

   ! loop over unique pseudo atoms to read nml=unique_pseudo_atom

   do i = 1, n_unique_pseudo_atoms

   if(.not.input_line_is_namelist("unique_pseudo_atom")) call input_error( &
       "unique_pseudo_atom_read: &
       &namelist unique__pseudo_atom expected")

   ua_pseudo => unique_pseudo_atoms(i)

   name                 = df_name
   zc                   = df_zc
   n_equal_pseudo_atoms = df_n_equal_pseudo_atoms

   call input_read_to_intermediate()
   read(unit, nml=unique_pseudo_atom, iostat=status)
   if (status .gt. 0) call input_error( &
      "unique__pseudo_atom_read: namelist unique_pseudo_atom.")

   ua_pseudo%name                 = name
   ua_pseudo%zc                   = zc
   ua_pseudo%n_equal_pseudo_atoms = n_equal_pseudo_atoms

      allocate( ua_pseudo%position_pseudo(3,n_equal_pseudo_atoms), stat=status )
      if ( status .ne. 0 ) call error_handler( &
           "unique_pseudo_atom_read: allocate of position failed")

   ! read coordinates of the first equal pseudo atom

   do j = 1, n_equal_pseudo_atoms
   call input_read_to_intermediate()
   read(unit, fmt=*, iostat=status) (ua_pseudo%position_pseudo(k,j), k=1,3)
   if(status.ne.0)call input_error("unique_pseudo_atom_read: position.")
   enddo

      ua_pseudo%position_first_pseudo = ua_pseudo%position_pseudo(:,1)

   end do

   end subroutine unique_pseudo_atom_read
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine unique_pseudo_atom_read_basis()
   ! purpose: reads in basis set for pseudo potentials
   !
   ! called by read_input()
   !---- end of interface ---------------------------------------------------------

   !---- used modules -------------------------------------------------------------
   use input_module
   implicit none
   !---- declaration of local variables -------------------------------------------
   type(unique_pseudo_atoms_type), pointer   :: ua_pseudo
   integer(kind=i4_kind)                     :: i,j,unit,status
   !---- declaration of subroutines used ------------------------------------------
   external error_handler
   !---- Excutable codes ----------------------------------------------------------

   ! loop over unique pseudo atoms to read nml=unique_pseudo_atom_basisset
   ! and basissets for pseudo potentials

   do i = 1, N_equal_pseudo_atoms

     if(.not. input_line_is_namelist("unique_pseudo_atom_basisset") ) &
       call  input_error("unique_pseudo_atom_read_basis: &
             &namelist unique_pseudo_atom_basisset expected")

      ua_pseudo => unique_pseudo_atoms(i)

      l_max1 = df_l_max1

     call input_read_to_intermediate()

     unit = input_intermediate_unit()

     read(unit, nml=unique_pseudo_atom_basisset, iostat=status)
     if(status .gt. 0) call input_error( &
          "unique_pseudo_atom_read_basis:&
          & namelist unique_pseudo_atom_basisset.")
     if(l_max1.lt.0) call input_error( &
           "unique_pseudo_atom_read_basis:&
           & namelist unique_pseudo_atom_basisset: l_max1")
      call unique_pseudo_atom_alloc(ua_pseudo,l_max1)

      ! read l dependent basis information for pseudo potentials
      do j=0,l_max1
         call unique_pseudo_atom_basis_read(&
                     ua_pseudo%l_ob_pseudo(j))
      enddo
   end do

   end subroutine unique_pseudo_atom_read_basis
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine unique_pseudo_atom_basis_read(l)
   !  Purpose: read in all information in l and does allocation
   use input_module
   !---------------------- Declaration of formal parameters --------------------
   type(unique_pseudo_atoms_basis_type), intent(inout) :: l
   !--- end of interface -------------------------------------------------------
   !--- declaration of local variables -----------------------------------------
   integer(kind=i4_kind)                              :: status,unit
   !--- declaration of subroutines used ----------------------------------------
   external error_handler
   !--- executable codes -------------------------------------------------------

   if(.not.input_line_is_namelist("unique_pseudo_basis") ) call &
    input_error("unique_pseudo_basis_read:&
    & namelist  unique_pseudo_basis expected")
    unit = input_intermediate_unit()

   n_terms  = df_n_terms
   call input_read_to_intermediate()
   read(unit, nml=unique_pseudo_basis, iostat=status)
   if (status.gt.0) call input_error( &
        "unique_pseudo_atom_basis_read: namelist unique_pseudo_basis.")
   if (n_terms.le.0) call input_error( &
        "unique_pseudo_atom_basis_read: namelist unique_pseudo_basis:&
         & n_terms")
   call unique_pseudo_atom_basis_alloc(l,n_terms)
   ! read symmetry of Gaussian, exponentials and coefficients
   read(unit, fmt=*, iostat=status) l%n_or_sym, l%cexponents, l%dcoefficients
         if(status.gt.0) call input_error( &
           "unique_pseudo_basis_read: n_or_sym,exponents,dcoefficients.")

   end subroutine unique_pseudo_atom_basis_read
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine unique_pseudo_atom_basis_alloc(l,n_terms)
   !  Purpose: allocates arrays contained in l and initialises l
   ! ----- Declaration of formal parameters ------------------------------------
   type(unique_pseudo_atoms_basis_type),  intent(inout)   :: l
   integer(kind=i4_kind),                intent(in)      :: n_terms
   ! ==== end of interface =====================================================
   ! ---- Declaration of local variables ---------------------------------------
   integer(kind=i4_kind)                                 :: status
   ! ---- Declaration of subroutines used --------------------------------------
   external error_handler
   ! ---- Executable codes -----------------------------------------------------

   l%n_terms   = n_terms

   allocate(l%n_or_sym(n_terms),l%cexponents(n_terms),&
            l%dcoefficients(n_terms), stat=status )
   if(status.ne.0) call error_handler( &
   "unique_pseudo_atom_basis_alloc: allocate of exponents failed")

   end subroutine unique_pseudo_atom_basis_alloc
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine unique_pseudo_atom_pack()
   !  purpose : packs all information contained in
   !            unique_pseudo_atoms(:)
   !---------- end of interface -------------------------------------------------
   !--- modules used  -----------------------------------------------------------
   use comm_module, only: commpack
   !--- Declaration of local variables ------------------------------------------
   type(unique_pseudo_atoms_type), pointer    :: ua_pseudo
   integer(kind=i4_kind)                     :: i, j, status, length
   !--- Declaration of subroutines used -----------------------------------------
   external error_handler
   !---------- Executable codes -------------------------------------------------
   call commpack(n_unique_pseudo_atoms,status)
   if(status.gt.0) call error_handler( &
     "unique_pseudo_atom_pack: error at N_unique_pseudo_atoms.")

   ! loop over unique_pseudo atoms
   do i = 1, n_unique_pseudo_atoms

      ua_pseudo => unique_pseudo_atoms(i)

      call commpack(ua_pseudo%zc, status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_pack: error at ZC.")
      call commpack(ua_pseudo%N_equal_pseudo_atoms,status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_pack: error at N_equal_pseudo_atoms.")
      call commpack(ua_pseudo%l_max1, status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_pack: error at l_max1.")
      call commpack(ua_pseudo%name,status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_pack: error at name.")

      ! pack atom coordinates of all equal atoms
        length = ua_pseudo%N_equal_pseudo_atoms * 3
      call commpack(ua_pseudo%position_pseudo(1,1), length, 1, status)
      if(status.ne.0) call error_handler(&
       "unique_pseudo_atom_pack: error at position_pseudo.")

      ! pack l dependent basis set for pseudo potentials

      do j = 0, ua_pseudo%l_max1
         call unique_pseudo_atom_basis_pack(ua_pseudo%l_ob_pseudo(j))
      end do

    end do

  end subroutine unique_pseudo_atom_pack
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unique_pseudo_atom_unpack()
  ! Purpose: unpacks all information contained in
  !          unique_pseudo_atoms(:) and does allocation
  !-----------------end of interface --------------------------------------------
  !----- modules used -----------------------------------------------------------
  use comm_module, only: communpack
  !----- Declaration of local variables -----------------------------------------
  type(unique_pseudo_atoms_type), pointer   :: ua_pseudo
  integer(kind=i4_kind)                    :: i, j, status, length
  !--- Declaration of subroutines used -----------------------------------------
   external error_handler
  !---------- Executable codes -------------------------------------------------

  call communpack(n_unique_pseudo_atoms,status)
   if(status.gt.0) call error_handler( &
     "unique_pseudo_atom_unpack: error at N_unique_pseudo_atoms.")

   allocate( unique_pseudo_atoms(N_unique_pseudo_atoms), stat=status )
   if (status .ne. 0) call error_handler( &
        "unique_pseudo_atom_unpack: allocate failed.")


   ! loop over unique_pseudo atoms
   do i = 1, n_unique_pseudo_atoms

      ua_pseudo => unique_pseudo_atoms(i)

      call communpack(ua_pseudo%zc, status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_unpack: error at ZC.")
      call communpack(ua_pseudo%N_equal_pseudo_atoms,status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_unpack: error at N_equal_pseudo_atoms.")
      call communpack(ua_pseudo%l_max1, status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_unpack: error at l_max1.")
      call communpack(ua_pseudo%name,status)
      if(status.ne.0) call error_handler(&
        "unique_pseudo_atom_unpack: error at name.")

      call unique_pseudo_atom_alloc(ua_pseudo,l_max1,n_equal_pseudo_atoms)

      ! read atom coordinates of all equal atoms
        length = ua_pseudo%N_equal_pseudo_atoms * 3
      call communpack(ua_pseudo%position_pseudo(1,1), length, 1, status)
      if(status.ne.0) call error_handler(&
       "unique_pseudo_atom_unpack: error at position_pseudo.")

      ! read l dependent basis set for pseudo potentials
      do j = 0, ua_pseudo%l_max1
         call unique_pseudo_atom_basis_unpack(ua_pseudo%l_ob_pseudo(j))
      end do

    end do

  end subroutine unique_pseudo_atom_unpack
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unique_pseudo_atom_basis_pack(l)
  !  Purpose: pack all information in l
  !--- end of interface --------------------------------------------------------
  !--- modules used ------------------------------------------------------------
  use comm_module, only: commpack
  !--- Declaration of formal parameters ----------------------------------------
  type(unique_pseudo_atoms_basis_type), intent(inout) :: l
  !--- Declaration of local variables ------------------------------------------
  integer(kind=i4_kind)                         :: status
  !--- Declaration of subroutines used -----------------------------------------
  external error_handler
  !--- Executable codes --------------------------------------------------------

  call commpack(l%n_terms,status)
    if(status.ne.0) call error_handler(&
    "unique_pseudo_atom_basis_pack: error at n_terms.")

  ! pack symmetry index of Gaussians

  call commpack(l%n_or_sym, l%n_terms, 1, status)
    if(status.ne.0) call error_handler(&
   "unique_pseudo_atom_basis_pack: error at n_or_sym.")

  ! pack exponents for Gaussians

  call commpack(l%cexponents, l%n_terms, 1, status)
   if(status.ne.0) call error_handler(&
  "unique_pseudo_atom_basis_pack: error at cexponents.")

  ! pack coefficients for Gaussians

  call commpack(l%dcoefficients, l%n_terms, 1, status)
   if(status.ne.0) call error_handler(&
  "unique_pseudo_atom_basis_pack: error at coefficients.")

  end subroutine unique_pseudo_atom_basis_pack
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unique_pseudo_atom_basis_unpack(l)
  !  purpose : read in all information in l and dose allocation
  !--- end of interface --------------------------------------------------------
  !---------------------- modules used -----------------------------------------
  use comm_module, only: communpack
  !----- Declaration of formal parameters --------------------------------------
  type(unique_pseudo_atoms_basis_type),  intent(inout) :: l
  !----- Declaration of local variables ----------------------------------------
  integer(kind=i4_kind)                               :: status
  !----- Declaration of subroutines used ---------------------------------------
  external error_handler
  !----- Executable codes ------------------------------------------------------

  call communpack(n_terms,status)
    if(status.ne.0) call error_handler(&
    "unique_pseudo_atom_basis_unpack: error at n_terms.")

  call unique_pseudo_atom_basis_alloc(l,n_terms)

  ! unpack symmetry index of Gaussians

  call communpack(l%n_or_sym, l%n_terms, 1, status)
    if(status.ne.0) call error_handler(&
   "unique_pseudo_atom_basis_unpack: error at n_or_sym.")

  ! unpack exponents for Gaussians

  call communpack(l%cexponents, l%n_terms, 1, status)
   if(status.ne.0) call error_handler(&
  "unique_pseudo_atom_basis_unpack: error at cexponents.")

  ! unpack coefficients for Gaussians

  call communpack(l%dcoefficients, l%n_terms, 1, status)
   if(status.ne.0) call error_handler(&
  "unique_pseudo_atom_basis_unpack: error at coefficients.")

  end subroutine unique_pseudo_atom_basis_unpack
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unique_pseudo_atom_write(iounit)
  !
  ! Purpose: write name and zc and coordinates of all centers for each
  ! unique pseudo atoms.
  !
  use echo_input_module
  use operations_module, only: operations_echo_input_level
  implicit none
  integer(kind=i4_kind)             :: iounit
  !---------------------- end of interface -------------------------------------

  type(unique_pseudo_atoms_type),     pointer :: ua_pseudo
  integer(kind=i4_kind)                      :: i, j, k, status
  character(len=29)                          :: header
  !------------------------------------------------------------------------------
  external error_handler
  !----------- Executable codes -------------------------------------------------

  call start("UNIQUE_PSEUDO_ATOM_NUMBER", "UNIQUE_PSEUDO_ATOM_WRITE", &
       iounit,operations_echo_input_level)
  call intg("N_UNIQUE_PSEUDO_ATOMS", n_unique_pseudo_atoms,&
             df_n_unique_pseudo_atoms)
  call stop

  write(iounit,fmt='(A/)', iostat=status) " # << Geometry >>"
  if (status /= 0) call error_handler &
        ("UNIQUE_PSEUDO_ATOM_WRITE: writing of geometry banner failed")

   word_format = '("    ",a," = ",a14  :" # ",a)' ! including quotes

  ! loop over unique_pseudo atoms

  do i = 1, n_unique_pseudo_atoms

     ua_pseudo => unique_pseudo_atoms(i)

      write(header,'("UNIQUE_PSEUDO_ATOM # unique atom",i4)') i
      call start(header,"UNIQUE_PSEUDO_ATOM_WRITE", &
           iounit,operations_echo_input_level)
      call word("NAME          ", ua_pseudo%name                ,&
           &df_name                )
      call real("Z             ", ua_pseudo%zc                  ,&
           &df_Zc                  )
      call intg("N_pseudo_equal", ua_pseudo%N_equal_pseudo_atoms,&
           &df_N_equal_pseudo_atoms)
      call stop(empty_line=.false.)

      ! write pseudo atom coordinates of all pseudo equal atoms
      do j=1,ua_pseudo%N_equal_pseudo_atoms
         write(iounit, fmt='(3ES25.15)', iostat=status)&
              (ua_pseudo%position_pseudo(k,j), k=1,3)
         if (status .gt. 0) call error_handler(&
            "unique_pseudo_atom_write: error at positions.")
      enddo
      write(iounit, fmt='()', iostat=status)
      if (status .gt. 0) call error_handler(&
         "unique_pseudo_atom_write: error after positions.")

   enddo

   end subroutine unique_pseudo_atom_write
   !===============================================================================


   subroutine unique_pseudo_unique_atom_read()
   ! purpose: reads in names and zc of unique pseudo atoms
   !          and coordinates of first centers ( " equal pseudo atoms " )
   ! called by read_input()
   !---- end of interface ---------------------------------------------------------

   !---- used modules -------------------------------------------------------------
   use input_module
   use iounitadmin_module
   use filename_module
   !-------------------------------------------------------------------------------
   implicit none
   !---- declaration of local variables -------------------------------------------
   type(unique_pseudo_atoms_type), pointer     :: ua_pseudo
   integer(kind=i4_kind)                     :: i,k,status,unit
   !---- declaration of subroutines used ------------------------------------------
   external error_handler
   !---- Excutable codes ----------------------------------------------------------

   if ( .not. input_line_is_namelist("unique_pseudo_unique_atom_number") )&
           call input_error( &
    "unique_pseudo_uniqe_atom_read: namelist unique_pseudo_atom_number expected")

   unit = input_intermediate_unit()

   n_unique_pseudo_atoms = df_n_unique_pseudo_atoms
      call input_read_to_intermediate()
      read(unit, nml=unique_pseudo_atom_number, iostat=status)
      if (status .ne. 0) call input_error( &
           "unique_pseudo_unique_atom_read: namelist unique_pseudo_atom_number.")
   if(n_unique_pseudo_atoms.lt.0)call input_error( &
   "unique_pseudo_atom_read: namelist unique_pseudo_atom_number.")

   allocate( unique_pseudo_atoms(n_unique_pseudo_atoms), stat=status )
   if (status .ne. 0) call error_handler( &
        "unique_pseudo_unique_atom_read: allocate failed.")

   ! loop over unique pseudo atoms to read nml=unique_pseudo_atom

      do i = 1, n_unique_pseudo_atoms

        if(.not.input_line_is_namelist("unique_pseudo_atom")) call input_error( &
        "unique_pseudo_unique_atom_read: namelist unique__pseudo_atom expected")

        ua_pseudo => unique_pseudo_atoms(i)

        name                 = df_name
        zc                   = df_zc
        n_equal_pseudo_atoms = df_n_equal_pseudo_atoms

        call input_read_to_intermediate()
        read(unit, nml=unique_pseudo_atom, iostat=status)
        if(status.gt.0) call input_error( &
        "unique__pseudo_atom_read: namelist unique_pseudo_atom.")

        ua_pseudo%name                 = name
        ua_pseudo%zc                   = zc
        ua_pseudo%n_equal_pseudo_atoms = n_equal_pseudo_atoms

        ! read coordinates of the first equal pseudo atom

        call input_read_to_intermediate()
        read(unit, fmt=*, iostat=status)&
        (ua_pseudo%position_first_pseudo(k), k=1,3)
        if(status.gt.0) call input_error(&
        "unique_pseudo_unique_atom_read: positions.")

     enddo

   end subroutine unique_pseudo_unique_atom_read
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine unique_pseudo_unique_atom_write(iounit)
  !
  ! Purpose: write name and zc and coordinates of all centers for each
  ! unique pseudo atoms.
  !
  use echo_input_module
  use operations_module, only: operations_echo_input_level
  implicit none
  integer(kind=i4_kind)             :: iounit
  !---------------------- end of interface -------------------------------------

  type(unique_pseudo_atoms_type),     pointer :: ua_pseudo
  integer(kind=i4_kind)                      :: i, k, status
  character(len=29)                          :: header
  !------------------------------------------------------------------------------
  external error_handler
  !----------- Executable codes -------------------------------------------------

  call start("UNIQUE_PSEUDO_ATOM_NUMBER", "UNIQUE_PSEUDO_ATOM_WRITE", &
       iounit,operations_echo_input_level)
  call intg("N_UNIQUE_PSEUDO_ATOMS", n_unique_pseudo_atoms,&
             df_n_unique_pseudo_atoms)
  call stop

  write(iounit,fmt='(A/)', iostat=status) " # << Geometry >>"
  if (status /= 0) call error_handler &
        ("UNIQUE_PSEUDO_UNIQUE_ATOM_WRITE: writing of geometry banner failed")

   word_format = '("    ",a," = ",a14  :" # ",a)' ! including quotes

  ! loop over unique_pseudo atoms

  do i = 1, n_unique_pseudo_atoms

     ua_pseudo => unique_pseudo_atoms(i)

      write(header,'("UNIQUE_PSEUDO_ATOM # unique atom",i4)') i
      call start(header,"UNIQUE_PSEUDO_ATOM_WRITE", &
           iounit,operations_echo_input_level)
      call word("NAME          ", ua_pseudo%name                ,&
           &df_name                )
      call real("Z             ", ua_pseudo%zc                  ,&
           &df_Zc                  )
      call intg("N_pseudo_equal", ua_pseudo%N_equal_pseudo_atoms,&
           &df_N_equal_pseudo_atoms)
      call stop(empty_line=.false.)

      ! write pseudo atom coordinates of all pseudo equal atoms
         write(iounit, fmt='(3ES25.15)', iostat=status)&
              (ua_pseudo%position_first_pseudo(k), k=1,3)
         if (status.gt.0) call error_handler(&
            "unique_pseudo_atom_write: error at positions.")
      write(iounit, fmt='()', iostat=status)
      if (status.gt.0) call error_handler(&
         "unique_pseudo_atom_write: error after positions.")

   enddo

   end subroutine unique_pseudo_unique_atom_write
   !==============================================================================
   end module  pseudo_potential_module


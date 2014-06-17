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
!==============================================================
! Public interface of module
!===============================================================
module unique_atom_methods
  !
  ! SEE UNIQUE_ATOM_MODULE ...
  !
  !  Database   for  Coordinates   and   Coefficients  of   primitive,
  !  contracted and symmetry adapted basis functions.  Subroutines for
  !  reading  and  writing this  database  are  also defined.   Memory
  !  allocation   is   automatically   done  by   these   subroutines.
  !  Subroutines doing memory allocation are made public for data that
  !  will be calculated externally.
  !
  !  Due to possible needs to replace parts of reading by internal and
  !  possibly parallel calculations, the  reading, writing of data was
  !  divided  over several  routines that  have  to be  called in  the
  !  following sequence (replace xxx by read or write):
  !
  !    unique_atom_xxx
  !
  !  or for new symetry part:
  !
  !    unique_atom_unique_xxx
  !    unique_atom_xxx_basis
  !    unique_atom_symadapt_xxx
  !    unique_atom_symequiv_xxx
  !
  !  unique_atom_setup() must be called after symmetry part or reading
  !  in  input to calculate  some quantities  needed in  integlral and
  !  scf-part!
  !
  !  References: publisher document: orbital_calculation
  !
  !
  !  Author: TB
  !  Date: 18.07.95
  !
  !== Interrupt of public interface of module =====================
  !
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   29.4.96
  ! Description: implementation of contractions changed
  !
  ! Modification (Please copy before editing)
  ! Author: FN
  ! Date:   8/96
  ! Description: Shift the calculation of the fitbasis
  !              dimensions to the fit_coeff_module due
  !              to cross-dependencies (use statements)
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   11/96
  ! Description:
  !   moving the reading of global contractions in input file:
  !   they should be read in together with the rest of basis set
  !   for a given unique atom.
  !   Deleted unique_atom_glob_cons_read/write.
  !   Reading/Writing of global contractions now in unique_atom_read/write.
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: Uwe Birkenheuer
  ! Date:   June 1998
  ! Description: Moving_Unique_Atom concept introduced
  !
  ! Modification
  ! Author: AM
  ! Date:   05.99
  ! Description: modified subs:
  !                 unique_atom_... _setup,
  !                                 _calc_irrep_dims
  !                                 _elim_irreps
  !              last splitted into two:
  !                                 _elim_empty_syadapt
  !
  !----
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

!------------ Modules used --------------------------------------
#include "def.h"
use type_module  ! type specification parameters
use unique_atom_module,  only : unique_atom_type                               &
                              , unique_atom_basis_type                         &
                              , unique_atom_glob_con_type                      &
                              , unique_atom_partner_type                       &
                              , unique_atom_symadapt_type                      &
                              , unique_atom_pseudopot_type                     &
                              , unique_atom_atomic_dens_type                   &
                              , unique_atoms                                   &
                              , unique_atom_symequiv                           &
                              , unique_atom_grad_info                          &
                              , moving_unique_atom_index                       &
                              , N_unique_atoms                                 &
                              , N_moving_unique_atoms                          &
                              , core_density_setup                             &
                              , pseudopot_present                              &
                              , unique_atom_lmax_all                           &
                              , unique_atom_lmax_ob                            &
                              , unique_atom_lmax_pseudo                        &
                              , unique_atom_lmax_ch                            &
                              , unique_atom_lmax_xc
use symmetry_data_module ! database for irreps
use datatype, only : arrmat3
#ifdef FPP_DEBUG
     use error_module, only: MyID
#endif
implicit none
save ! save all variables defined in this module
private ! by default, everything is private
!== Interrupt end of public interface of module =================

!------------ Declaration of types ------------------------------
!
! ALL TYPES ARE DECLARED IN UNIQUE_ATOM_MODULE ...
!

!------------ public functions and subroutines ------------------


public :: unique_atom_core_charge ! type (unique_atom_type) -> real

public :: unique_atom_read
public :: unique_atom_write

public :: unique_atom_read_basis!()
public :: unique_atom_write_basis!()

#ifdef WITH_GUILE
public :: unique_atom_find_basis!(ua)
#endif

public :: unique_atom_setup
public :: unique_atom_close

public :: unique_atom_grad_information

public :: unique_atom_pseudopot_read
public :: unique_atom_pseudopot_bcast

public :: unique_atom_symequiv_write
public :: unique_atom_symadapt_write

public :: unique_atom_calc_lmax
public :: unique_atom_make_gx

public :: unique_atom_alloc
public :: unique_atom_assign_symm

!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of private variables ------------------
! defaults for namelist input
integer           , private, parameter :: df_N_unique_atoms      =  0
character(len=12)      , private, parameter :: df_name                =  "            "
real(kind=r8_kind), private, parameter :: df_Z                   =  0.0_r8_kind
real(kind=r8_kind), private, parameter :: df_ZC                  =  0.0_r8_kind
real(kind=r8_kind), private, parameter :: df_nuclear_magnetic_moment =  0.0_r8_kind
real(kind=r8_kind), private, parameter :: df_nuclear_radius      =  0.0_r8_kind
integer           , private, parameter :: df_nuclear_spin        =  0
integer           , private, parameter :: df_N_equal_atoms       =  0
logical           , private, parameter :: df_fixed               = .false.
integer           , private, parameter :: df_lmax_ob             = -1
integer           , private, parameter :: df_lmax_pseudo         = -1
integer           , private, parameter :: df_lmax_ch             = -2
integer           , private, parameter :: df_lmax_xc             = -2
integer           , private, parameter :: df_N_glob_cons_ch      =  0
integer           , private, parameter :: df_N_glob_cons_xc      =  0
integer           , private, parameter :: df_N_exponents         =  0
logical           , private, parameter :: df_automatic           =  .false.
integer           , private, parameter :: df_N_uncontracted_fcts =  0
integer           , private, parameter :: df_N_contracted_fcts   =  0
integer           , private, parameter :: df_N_fcts              =  0
integer           , private, parameter :: df_N_contributing_fcts =  0
integer           , private, parameter :: df_N_exponents_core    =  0
logical           , private, parameter :: df_core_density_setup  = .false.
logical           , private, parameter :: df_core_density_use    = .false.

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  elemental function unique_atom_core_charge (ua) result (q)
    implicit none
    type (unique_atom_type), intent (in) :: ua
    real (r8_kind) :: q
    ! *** end of interface ***

    ! q is the nuclear charge Z screened by ZC core electrons:
    q = ua % z - ua % zc
  end function unique_atom_core_charge

   !*************************************************************
   subroutine unique_atom_alloc(ua, &
        lmax_ob, lmax_ch, lmax_xc, N_equal_atoms, lmax_pseudo )
   !  Purpose: allocates arrays contained in ua and initialises ua
   !           with the exception of symadapt and renorm
   !           ua%lmax_all must be defined before first call
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_type), intent(inout) :: ua
   integer, optional,      intent(in)    :: lmax_ob
   integer, optional,      intent(in)    :: lmax_ch
   integer, optional,      intent(in)    :: lmax_xc
   integer, optional,      intent(in)    :: N_equal_atoms
   integer, optional,      intent(in)    :: lmax_pseudo
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer                                :: status
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if ( present(lmax_ob) ) then
      ua%lmax_ob = lmax_ob
      allocate( ua%l_ob(0:lmax_ob), stat=status )
      if ( status .ne. 0 ) call error_handler( &
           "unique_atom_alloc: allocate of l_ob failed")
      ua%lmax_all = max(ua%lmax_ob,ua%lmax_all)
   endif

   if( present(lmax_pseudo) )then
      ua%lmax_pseudo = lmax_pseudo
      allocate(ua%l_pseudopot(0:ua%lmax_pseudo),stat=status )
      if ( status .ne. 0 ) call error_handler( &
         "unique_atom_alloc: allocate of l_pseudopot failed")
   endif

   if ( present(lmax_ch) ) then
      ua%lmax_ch = lmax_ch
      allocate( ua%l_ch(0:lmax_ch), stat=status )
      if ( status .ne. 0 ) call error_handler( &
           "unique_atom_alloc: allocate of l_ch failed")
      ua%lmax_all = max(ua%lmax_ch,ua%lmax_all)
   endif

   if ( present(lmax_xc) ) then
      ua%lmax_xc = lmax_xc
      allocate( ua%l_xc(0:lmax_xc), stat=status )
      if ( status .ne. 0 ) call error_handler( &
           "unique_atom_alloc: allocate of l_xc failed")
      ua%lmax_all = max(ua%lmax_xc,ua%lmax_all)
   endif

   if ( present(N_equal_atoms) ) then
      ua%N_equal_atoms = N_equal_atoms
      allocate( ua%position(3,N_equal_atoms), stat=status )
      if ( status .ne. 0 ) call error_handler( &
           "unique_atom_alloc: allocate of position failed")
   endif
   end subroutine unique_atom_alloc
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_basis_alloc( l, N_exponents, &
                    N_uncontracted_fcts, N_contracted_fcts, N_fcts )
   !  Purpose: allocates arrays contained in l and initialises l
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_basis_type), intent(inout) :: l
   integer(kind=i4_kind),  intent(in)          :: N_exponents
   integer(kind=i4_kind),  intent(in)          :: N_uncontracted_fcts
   integer(kind=i4_kind),  intent(in)          :: N_contracted_fcts
   integer(kind=i4_kind),  intent(in)          :: N_fcts
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer                                  :: status
   integer :: N_bas
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if ( N_contracted_fcts .gt. N_exponents ) call error_handler( &
        "unique_atom_basis_alloc: invalid N_contracted_fcts" )
   l%N_exponents = N_exponents
   l%N_contracted_fcts = N_contracted_fcts
   l%N_uncontracted_fcts = N_uncontracted_fcts
   allocate( l%exponents(N_exponents), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_basis_alloc: allocate of exponents failed")
   allocate( l%norms(N_exponents), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_basis_alloc: allocate of norms failed")
   allocate(l%contractions(N_exponents,N_contracted_fcts),stat=status)
   if ( status .ne. 0 ) call error_handler(&
        & "unique_atom_basis_alloc: allocate of contractions failed")

   ! basis dimension:
   N_bas = N_uncontracted_fcts + N_contracted_fcts

   allocate(l%functions(N_bas, N_fcts), stat=status)
   ASSERT(status==0)
   end subroutine unique_atom_basis_alloc
   !*************************************************************

   !*************************************************************
   subroutine unique_atom_glob_con_alloc(c, N_contributing_fcts)
   !  Purpose: allocates arrays contained in c and initialises c
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_glob_con_type), intent(inout)  :: c
   integer,                 intent(in)     :: N_contributing_fcts
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer                                  :: status
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   c%N_contributing_fcts = N_contributing_fcts
   allocate( c%l(N_contributing_fcts), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_glob_con_alloc: allocate of l failed")
   allocate( c%index_ind_fct(N_contributing_fcts), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_glob_con_alloc: allocate of index_ind_fct failed")
   allocate( c%index_exp(N_contributing_fcts), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_glob_con_alloc: allocate of index_exp failed")
   allocate( c%coefs(N_contributing_fcts), stat=status )
   if ( status .ne. 0 ) call error_handler( &
        "unique_atom_glob_con_alloc: allocate of coefs failed")
    end subroutine unique_atom_glob_con_alloc
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_glob_cons_alloc( ua, &
        N_glob_cons_ch, N_glob_cons_xc )
   !  Purpose: allocates arrays contained in ua and initialises ua
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_type), intent(inout)  :: ua
   integer, optional, intent(in) :: N_glob_cons_ch, N_glob_cons_xc
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer                                :: status
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if ( present(N_glob_cons_ch) ) then
      ua%N_glob_cons_ch = N_glob_cons_ch
      if ( N_glob_cons_ch .gt. 0 ) then
         allocate( ua%glob_con_ch(N_glob_cons_ch), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "unique_atom_alloc: allocate of glob_cons_ch failed")
      else
         ! leave it deallocated (ua%glob_con_ch)
      endif
   endif
   if ( present(N_glob_cons_xc) ) then
      ua%N_glob_cons_xc = N_glob_cons_xc
      if ( N_glob_cons_xc .gt. 0 ) then
         allocate( ua%glob_con_xc(N_glob_cons_xc), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "unique_atom_alloc: allocate of glob_cons_xc failed")
      else
         ! leave it deallocated (ua%glob_con_xc)
      endif
   endif
   end subroutine unique_atom_glob_cons_alloc
   !*************************************************************

   !*************************************************************
   subroutine  unique_atom_read (loop) ! was unique_atom_unique_read()
     !
     ! Reads  in  names, fixed  option,  and  Z  of unique  atoms  and
     ! coordinates  of   first  centers  ("equal   atoms")  Also  load
     ! N_moving_unique_atoms.   Also allocates  and  reads gxepe_array
     ! and gxepe_impu for epe treatment
     !
   use input_module, only: input_line_is_namelist, &
        input_intermediate_unit, input_error, input_read_to_intermediate
   use operations_module, only: operations_geo_opt, operations_core_density, &
        operations_read_gx, operations_qm_mm
#ifdef WITH_EFP
   use operations_module, only: operations_integral, operations_scf, &
        operations_post_scf
#endif
   use options_module, only: options_relativistic
   use unique_atom_module, only: unique_atom_iwork, &
        GLOB_n_unique_atoms => n_unique_atoms
   !^^^^ make it visible with another name, see comment below
   use atom_data_module, only: nuc_radius
   implicit none
   integer(kind=i4_kind), intent(in):: loop
   !** End of interface *****************************************

   type(unique_atom_type), pointer :: ua
   integer ::  i, k, status, unit
   !
   ! FIXME: Does this still affect us?  Somehow
   ! UNIQUE_ATOM_MODULE::N_UNIQUE_ATOMS is not changed by
   ! read(nml=...) below trying a workaround (Linux/MPI,Absoft6.0)
   !
   integer(i4_kind)   :: n_unique_atoms ! local variable, hides the module variable

   namelist /unique_atom_number/ N_unique_atoms

   character(len=12)                            :: name = ""
   real(kind=r8_kind)                      :: Z    = 0.0_r8_kind
   integer                                 :: N_equal_atoms = 0
   logical                                 :: fixed = .false.
   integer                                 :: nuclear_spin = 0
   real(kind=r8_kind)                      :: nuclear_magnetic_moment = 0.0_r8_kind
   real(kind=r8_kind)                      :: nuclear_radius = 0.0_r8_kind
   namelist /unique_atom/ Z, N_equal_atoms, name, fixed, &
                          nuclear_spin, nuclear_magnetic_moment, nuclear_radius
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if ( .not. input_line_is_namelist("unique_atom_number") ) call input_error( &
        "unique_atom_read: namelist unique_atom_number expected")

   unit = input_intermediate_unit()

   N_unique_atoms = df_N_unique_atoms
   call input_read_to_intermediate()

   read(unit, nml=unique_atom_number, iostat=status)

   if (status .gt. 0) call input_error( &
        "unique_atom_read: namelist unique_atom_number.")

#ifdef WITH_EFP
   if ( N_unique_atoms .le. 0 .and. .not. efp) call input_error( &
        "unique_atom_read: namelist unique_atom_number: N_unique_atoms .LE. 0")
#else
   if ( N_unique_atoms .le. 0 ) call input_error( &
        "unique_atom_read: namelist unique_atom_number: N_unique_atoms .LE. 0")
#endif
   if (operations_core_density .and. N_unique_atoms > 1) call input_error( &
        "unique_atom_read: core densities creation only possible in single atom calculations")
   ! set the global varibale:
   GLOB_n_unique_atoms = n_unique_atoms

#ifdef WITH_EFP
   if(N_unique_atoms == 0) then
      operations_integral = .false.
      operations_scf = .false.
      operations_post_scf = .false.
   end if
#endif

   allocate( unique_atoms(N_unique_atoms), stat=status )
   if (status .ne. 0) call error_handler( &
        "unique_atom_read: allocate failed.")

   ! loop over unique atoms to read nml=unique_atom and
   ! positions of equal atoms
   N_moving_unique_atoms = 0
   do i=1,N_unique_atoms

      if ( .not. input_line_is_namelist("unique_atom") ) call input_error( &
           "unique_atom_read: namelist unique_atom expected")

      ua => unique_atoms(i)

      name          = df_name
      Z             = df_Z
      N_equal_atoms = df_N_equal_atoms
      fixed         = df_fixed
      nuclear_spin  = df_nuclear_spin
      nuclear_magnetic_moment = df_nuclear_magnetic_moment
      nuclear_radius = df_nuclear_radius
      call input_read_to_intermediate
      read(unit, nml=unique_atom, iostat=status)
      if (status .gt. 0) call input_error( &
           "unique_atom_read: namelist unique_atom.")
      if (N_equal_atoms .le. 0) call input_error( &
           "unique_atom_read: namelist unique_atom: N_equal_atoms")
      if (operations_core_density .and. N_equal_atoms > 1) call input_error( &
         "unique_atom_read: core densities creation only possible in single atom calculations")
      ua%name = name
      ua%Z = Z
      ua%N_equal_atoms = N_equal_atoms
      ua%nuclear_spin = nuclear_spin
      ua%nuclear_magnetic_moment = nuclear_magnetic_moment
      if(nuclear_radius==df_nuclear_radius)then
         nuclear_radius = - nuc_radius(NINT(Z))
      endif
      ua%nuclear_radius = nuclear_radius
      if (fixed) then
         ua%moving_atom = 0
      else
         N_moving_unique_atoms = N_moving_unique_atoms + 1
         ua%moving_atom = N_moving_unique_atoms
      endif
      ua%N_glob_cons_ch = 0
      ua%N_glob_cons_xc = 0
      ! read atom coordinates of first equal atoms
      call input_read_to_intermediate
      read(unit, fmt=*, iostat=status) (ua%position_first_ea(k), k=1,3)
      if (status .gt. 0) call input_error("unique_atom_read: positions.")

      ! intialize UA %components:
      ua%Zc = 0.0_r8_kind

   enddo
   if (N_moving_unique_atoms < N_unique_atoms .and. options_relativistic) &
        call input_error("unique_atom_read: Sorry! &
        &Fixed atoms not yet completely implemented for relativistic forces")

   ! For geometry  optimization read coordinates  from gxfile.  FIXME:
   ! Happens    to   also    set   the    foreign    module   variable
   ! unique_atom_iwork!
   if (operations_geo_opt .or. operations_qm_mm .or. operations_read_gx) then
      call read_geometry (loop, unique_atoms, unique_atom_iwork)
   endif
   end subroutine unique_atom_read
   !*************************************************************

   !*************************************************************
   subroutine read_geometry (loop, unique_atoms, iwork)
     !
     ! Read  positions of uniqe  atoms from  gxfile or  elsewhere. Was
     ! part of unique_atom_read() before.
     !
     ! FIXME:  Happens  to  also   set  the  foreign  module  variable
     ! unique_atom_iwork!
     !
     use iounitadmin_module, only: openget_iounit, returnclose_iounit, &
          write_to_output_units, write_to_trace_unit, output_unit
     use operations_module, only: operations_transit, operations_gx_epeformat
#ifdef WITH_MOLMECH
     use operations_module, only: operations_qm_mm_new
#endif
     use filename_module, only: inpfile
#ifdef WITH_MOLMECH
     use qmmm_interface_module, only: gx, gx_qm, imomm, qm_mm, read_qmmm_input
#endif
#ifdef WITH_EPE
     use ewaldpc_module, only: gxepe_array, EX_GXEPE, gxepe_impu, n_epe_r
#endif
     implicit none
     integer, intent (in) :: loop
     type (unique_atom_type), intent (inout) :: unique_atoms(:)
     integer, intent (out) :: iwork
     ! *** end of interface ***

     real (r8_kind) :: x_coord, y_coord, z_coord, z_dummy, rwork
     character (len=1) :: cha_loop
     character (len=5) ::  gx_buff
     integer :: ieq_dummy, indexes(7), impu
     integer :: i, i_ua, io_u, status, counter_equal, gx_count
     logical :: ex_gx
#ifdef WITH_EPE
     integer :: io_gxepe
     real (r8_kind), dimension(3) :: r_gxepe
#endif

#ifdef WITH_MOLMECH
      nqm_mm_new: if (.not. operations_qm_mm_new) then
#endif
#ifdef WITH_EPE
      inquire(file= trim(inpfile('epe.r')), exist=ex_gxepe)
      if(ex_gxepe) then
         call gxepe_allocate
         io_gxepe = openget_iounit (file=trim (inpfile ('epe.r')), &
              status='old', form='formatted')
      endif ! ex_gxepe
#endif

      if (operations_transit) then
         write(cha_loop,'(i1)') loop
         inquire(file= trim(inpfile('gx.'//cha_loop)), exist=ex_gx)
      else
         inquire(file= trim(inpfile('gxfile')), exist=ex_gx)
      endif

      read_in_gxfile: if (ex_gx) then

         if (operations_transit) then
            io_u = openget_iounit(file=trim(inpfile('gx.'//cha_loop)), status='old', &
                 form='formatted')
         else
            io_u = openget_iounit(file=trim(inpfile('gxfile')), status='old', &
                 form='formatted')
         endif

         ! This is used for error reporting only:
         gx_count = 0
#ifdef WITH_EPE
         n_epe_r = 0
#endif
         do i_ua = 1, N_unique_atoms
            counter_equal = 1
            do while (counter_equal <= unique_atoms(i_ua) % n_equal_atoms)

               gx_count = gx_count + 1
               if (operations_gx_epeformat) then
                  read(io_u, *, ERR=100, END=101) z_dummy, x_coord, y_coord, z_coord, ieq_dummy, indexes, impu
#ifdef WITH_EPE
                  if (ex_gxepe .and. impu .ne. 0) read(io_gxepe, *) r_gxepe
#endif
               else
                  read(io_u, *, ERR=100, END=101) z_dummy, x_coord, y_coord, z_coord, ieq_dummy
               endif

               ! Ignore dummy  atoms which  are identified not  by the
               ! value of Z, but rather  by zero index of the group of
               ! equivalent atoms. Note how  the value of Z is ignored
               ! altogether, not only for dummy atoms.
               if (ieq_dummy == 0) cycle

               ! Note  how the xyz  coordinates of  all but  the first
               ! atom  in the  group of  symmetry equivalent  ones are
               ! ignored:
               if (counter_equal == 1) then
                  unique_atoms(i_ua) % position_first_ea(1) = x_coord
                  unique_atoms(i_ua) % position_first_ea(2) = y_coord
                  unique_atoms(i_ua) % position_first_ea(3) = z_coord
               endif
#ifdef WITH_EPE
               if (ex_gxepe) then
                  gxepe_impu(i_ua) = impu
                  if (impu .ne. 0) then
                     gxepe_array(i_ua) % position(:, counter_equal) = r_gxepe
                     n_epe_r = n_epe_r + 1
                  endif
               endif
#endif
               counter_equal = counter_equal + 1
            enddo
         enddo ! i_ua=1,N_unique_atoms


         gx_count = gx_count + 1
         read(io_u, *, iostat=status) rwork
         if (status .gt. 0) call error_handler &
              ("unique_atom_read: reading iwork from gxfile")
         iwork = int (-rwork, i4_kind)
         call returnclose_iounit (io_u)
      else !read_in_gxfile
         if(operations_gx_epeformat) call error_handler("unique_atom_read: GXFILE is absent")
         write(*,*)'**********************************************'
         write(*,*)'WARNING: !!!!!!!!!!!!!!!!!!'
         write(*,*)'GXFILE has to be read in, but not exists.'
         write(*,*)'Atomic coordinates from PG input will be used.'
         write(*,*)'Be carefull if you use internal coordinates,'
         write(*,*)'optimization cannot be done.'
         write(*,*)'**********************************************'
         call write_to_output_units('**********************************************')
         call write_to_output_units('WARNING: !!!!!!!!!!!!!!!!!!')
         call write_to_output_units('GXFILE has to be read in, but not exists.')
         call write_to_output_units('Atomic coordinates from PG input will be used.')
         call write_to_output_units('Be carefull if you use internal coordinates,')
         call write_to_output_units('optimization cannot be done.')
         call write_to_output_units('**********************************************')
         call write_to_trace_unit('**********************************************')
         call write_to_trace_unit('WARNING: !!!!!!!!!!!!!!!!!!')
         call write_to_trace_unit('GXFILE has to be read in, but not exists.')
         call write_to_trace_unit('Atomic coordinates from PG input will be used.')
         call write_to_trace_unit('Be carefull if you use internal coordinates,')
         call write_to_trace_unit('optimization cannot be done.')
         call write_to_trace_unit('**********************************************')
         iwork = 1
      end if read_in_gxfile

#ifdef WITH_EPE
      if (ex_gxepe) then
         call returnclose_iounit (io_gxepe)

         ! This only writes  to output. FIXME: not the  right place to
         ! do this.
         write(output_unit,*) &
              ' gxepe array to relate  QM and impurity cluster centers '
         do i_ua = 1, N_unique_atoms
            do counter_equal = 1, unique_atoms(i_ua) % n_equal_atoms
               write(output_unit,*) gxepe_array(i_ua) % position(:, counter_equal)
            enddo ! counter_equal=1,n_equal_atoms
         enddo ! i_ua=1,N_unique_atoms
      endif ! ex_gxepe
#endif

#ifdef WITH_MOLMECH
       else nqm_mm_new
          if (imomm) then
             call read_qmmm_input (rwork)
             iwork = int(-rwork, i4_kind)
             i = 0
             do i_ua = 1, N_unique_atoms
                counter_equal = 1
                do while (counter_equal <= unique_atoms(i_ua) % n_equal_atoms)
                   i = i + 1
                   if (gx(i) % i_unique .ne. 0) then
                      if (counter_equal == 1) then
                         unique_atoms(i_ua) % position_first_ea(1) = gx_qm(i) % x
                         unique_atoms(i_ua) % position_first_ea(2) = gx_qm(i) % y
                         unique_atoms(i_ua) % position_first_ea(3) = gx_qm(i) % z
                      endif
                      counter_equal = counter_equal + 1
                   endif
                enddo
             enddo! i_ua=1,N_unique_atoms
          elseif(qm_mm) then
          end if
       endif nqm_mm_new
#endif
   return

100     write(gx_buff,'(i5)') gx_count
        call error_handler("unique_atom_read: wrong format of GX file. Line - "//trim(gx_buff))
101     call error_handler("unique_atom_read:  attempt to read after end of GX file")
   end subroutine read_geometry
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_read_basis()
   !  Purpose: reads in basis sets for all unique atoms
   use unique_atom_module, only: unique_atom_make_empty_basis
   use input_module
   use operations_module, only: operations_gradients, operations_core_density,&
        operations_pseudobonds
#ifdef WITH_ERI4C
   use eri4c_options,     only : J_exact
#endif
   !** End of interface *****************************************

   !------------ Declaration of local variables -----------------
   type(unique_atom_type), pointer :: ua
   integer             ::  i,j,status,unit,i_gc
   logical             ::  pseudopot, core_dens
   logical             :: nofit_ch, nofit_xc
   logical             :: skipfit = .FALSE.

   integer            :: lmax_ob          = 0
   integer            :: lmax_pseudo      = 0
   integer            :: lmax_ch          = 0
   integer            :: lmax_xc          = 0
   integer            :: N_glob_cons_ch   = 0
   integer            :: N_glob_cons_xc   = 0
   real(kind=r8_kind) :: ZC               = 0.0_r8_kind
   logical            :: core_density_use = .false.

   namelist /unique_atom_basisset/        &
                         lmax_ob          &
                       , lmax_pseudo      &
                       , lmax_ch          &
                       , lmax_xc          &
                       , N_glob_cons_ch   &
                       , N_glob_cons_xc   &
                       , ZC               &
                       , core_density_use

   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------

   unit = input_intermediate_unit()
   core_density_setup= df_core_density_setup
   ! loop over unique atoms to read nml=unique_atom_basiset
   ! and basissets for orbitals, charge and exchange fitfunctions
   do i=1,N_unique_atoms

      if ( .not. input_line_is_namelist("unique_atom_basisset") ) &
           call input_error("unique_atom_read_basis: &
           &namelist unique_atom_basisset expected")

      ua => unique_atoms(i)

      lmax_ob        = df_lmax_ob
      lmax_pseudo    = df_lmax_pseudo
      lmax_ch        = df_lmax_ch
      lmax_xc        = df_lmax_xc
      N_glob_cons_ch = df_N_glob_cons_ch
      N_glob_cons_xc = df_N_glob_cons_xc
      core_density_use = df_core_density_use
      zc=df_zc
      call input_read_to_intermediate()
      read(unit, nml=unique_atom_basisset, iostat=status)

      if (status .gt. 0) call input_error( &
           "unique_atom_read_basis: namelist unique_atom_basisset.")
      if (lmax_ob .lt. 0) call input_error( &
           "unique_atom_read_basis: namelist unique_atom_basisset: lmax_ob")
      if (N_glob_cons_ch .lt. 0) call input_error( &
           "unique_atom_read_basis: namelist unique_atom_basisset: N_glob_cons_ch")
      if (N_glob_cons_xc .lt. 0) call input_error( &
           "unique_atom_read_basis: namelist unique_atom_basisset: N_glob_cons_xc")
      ua%lmax_all = 0
      if (operations_gradients .and. ua%moving_atom > 0) ua%lmax_all = 1

      if(.not.operations_pseudobonds) then
         ! pseudo stuff
         if(zc>ua%z) call input_error &
              ("unique_atom_read_basis namelist unique_atom_basisset: &
              & ZC larger than Z")
         if(zc<0) call input_error &
              ("unique_atom_read_basis namelist unique_atom_basisset: &
              & ZC below zero")
      end if

      ! record the charge of pp-core:
      ua%zc     = zc
      pseudopot = zc /= 0.0_r8_kind .and. .not.operations_core_density

      ! (fake?) core density is included with every pp-parametrization, expect it to come:
      core_dens = zc /= 0.0_r8_kind ! and not only when operations_core_density

      if (pseudopot .and. lmax_pseudo .lt. 0) call input_error( &
           "unique_atom_read_basis: namelist unique_atom_basisset: lmax_pseudo")

      ! set the global flag if ANY of the atoms matches the condition:
      pseudopot_present = pseudopot .and. (lmax_pseudo/=0) & ! FIXME: lmax_pseudo /= -1 ???
                          .or. pseudopot_present             ! if the flag was already set
      ! pseudopot_present is declared in unique_atom_module

!     FIXME: Logic behind is not clear:
!     if ((.not.pseudopot_present).and.(.not.options_spin_orbit)) then
!        if(lmax_pseudo/=0.and.pseudopot) &
!             pseudopot_present = .true.
!       ! results in problems for spin orbit run
!       ! because of different value on master and slaves
!       ! but required for regular run with pseudopotentials
!     end if
!
!     if (.not.core_density_setup) then
!        if(pseudopot_present.and.core_density_use) &
!             core_density_setup = .true.
!     endif

      ! expect core-density to come if input contained core_density_use=T:
      core_dens = core_dens .or. core_density_use

      nofit_ch = lmax_ch == df_lmax_ch
      if ( nofit_ch ) lmax_ch = 0
      nofit_xc = lmax_xc == df_lmax_xc
      if ( nofit_xc ) lmax_xc = 0
#ifdef WITH_ERI4C   
      skipfit = J_exact
#endif

      call unique_atom_alloc(ua, lmax_ob, lmax_ch, lmax_xc)

      if (pseudopot) then
         call unique_atom_alloc(ua, lmax_pseudo=lmax_pseudo )
      else
         ua%lmax_pseudo = -1
      endif

      call unique_atom_glob_cons_alloc(ua,N_glob_cons_ch, N_glob_cons_xc)

      ! read l dependent basis information for orbitals
      do j = 0, lmax_ob
         call read_shell("ob", j, ua%l_ob(j))
      enddo

      ! read pseudo potential data (if required)
      if (pseudopot) then
         ! read the local pseudo potential first
         call unique_atom_pseudopot_read(ua%l_pseudopot(lmax_pseudo))
         ! then read the semi-local contributions
         do j=0,lmax_pseudo-1
            call unique_atom_pseudopot_read(ua%l_pseudopot(j))
         enddo
      endif

#ifdef WITH_CORE_DENS
      ! read core density data (if required)
      if (core_dens) then
         ! first read r^2 -type atomic core fitting basis functions
         call unique_atom_core_density_read(ua%r2_core)
         ! then read s-type atomic core fititng basis functions
         call unique_atom_core_density_read(ua%s_core)
      else
         ! PP-bases contain UNIQUE_ATOM_CORE_DENSITY namelists,
         ! one for r2 ...
         call compat_core_density_skip(ua%r2_core) ! if any
         ! ... and one for s-fitfunctions:
         call compat_core_density_skip(ua%s_core) ! if any
      endif
#else
      ! PP-bases contain UNIQUE_ATOM_CORE_DENSITY namelists,
      ! one for r2 ...
      call compat_core_density_skip(ua%r2_core) ! if any
      ! ... and one for s-fitfunctions:
      call compat_core_density_skip(ua%s_core) ! if any
#endif

      if ( nofit_ch ) then
         !
         ! FIXME: some code  in the wild assumes r2-  and s- exponents
         ! are  always   present.   That  code  should   be  fixed  by
         ! respecting  %lmax_ch. This  is  a workaround  to make  that
         ! broken code happy:
         !
         ua%r2_ch   = unique_atom_make_empty_basis()
         ASSERT(size(ua%l_ch)>0)
         ua%l_ch(0) = unique_atom_make_empty_basis()
      else
         !
         ! Read  one-or more shells  from the  input, depending  on what
         ! lmax_ch says:
         !
         ASSERT(lmax_ch>-2)

         !
         ! Read  r**2 basis  information for  charge  fitfcts.  FIXME:
         ! consider treating r2  as yet another member of  an array at
         ! position -1.
         !
         call read_shell("ch", -1, ua%r2_ch, ua=ua)
         if (skipfit) then
            ua%r2_ch = unique_atom_make_empty_basis()
         endif

         do j = 0, lmax_ch
            ! ua is used to eventually autogenerate exponents:
            call read_shell("ch", j, ua%l_ch(j), ua=ua)
            if ( skipfit ) then
               ua%l_ch(j) = unique_atom_make_empty_basis()
            endif
         enddo
      endif

      ! read global contractions for charge fitfcts
      do i_gc = 1, N_glob_cons_ch
         call unique_atom_glob_con_read(ua%glob_con_ch(i_gc))
      enddo

#ifdef WITH_CORE_DENS
      ! link the core fit functions to the charge fit functions
      if (core_dens .and. operations_core_density) then
         call unique_atom_link_core_fitfct(ua%r2_core,ua%r2_ch)
         call unique_atom_link_core_fitfct(ua%s_core,ua%l_ch(0))
      endif
#endif

      if ( .not. nofit_xc ) then
         ! read r**2 basis information for exchange fitfcts
         call read_shell("xc", -1, ua%r2_xc)

         ! read l dependent basis information for exchange fitfcts
         do j = 0, lmax_xc
            call read_shell("xc", j, ua%l_xc(j))
         enddo
      else
         ua%r2_xc = unique_atom_make_empty_basis()
         do j=0,lmax_xc
            ua%l_xc(j) = unique_atom_make_empty_basis()
         enddo
      endif

      ! read global contractions for exchange fitfcts
      do i_gc = 1, N_glob_cons_xc
         call unique_atom_glob_con_read(ua%glob_con_xc(i_gc))
      enddo
   enddo
   end subroutine unique_atom_read_basis
   !*************************************************************

#ifdef WITH_GUILE
   subroutine unique_atom_find_basis(ua)
     !
     ! Purpose: find basis sets for a unique atom.
     !
     use unique_atom_module, only: unique_atom_type, &
          unique_atom_make_empty_basis
     use scheme, only: scheme_find_basis
     implicit none
     type(unique_atom_type), intent(inout) :: ua
     ! *** end of interface ***

     call scheme_find_basis (ua%name, "def-SVP", ua%l_ob)
     call scheme_find_basis (ua%name, "fake", ua%l_ch)
     call scheme_find_basis (ua%name, "fake", ua%l_xc)

     ua%lmax_ob = size(ua%l_ob) - 1
     ua%lmax_ch = size(ua%l_ch) - 1
     ua%lmax_xc = size(ua%l_xc) - 1
     ua%lmax_all = max(ua%lmax_ob, ua%lmax_ch, ua%lmax_xc)

     ! record the charge of pp-core:
     ua%zc = 0.0

     ! FIXME: not supported:
     call unique_atom_atomic_dens_alloc(ua%r2_core, 0)
     call unique_atom_atomic_dens_alloc(ua%s_core, 0)
     call unique_atom_glob_cons_alloc(ua, 0, 0)

     ! FIXME: get rid of special treatment of r2-bases:
     ua%r2_ch = unique_atom_make_empty_basis()
     ua%r2_xc = unique_atom_make_empty_basis()
   end subroutine unique_atom_find_basis
#endif

   !*************************************************************
   subroutine unique_atom_symequiv_free
     ! Purpose : deallocates sym.equivalent distances
     !** End of interface *****************************************
     implicit none
     !------------ Declaration of variables ---------------
     integer(kind=i4_kind)              :: status
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------

     if ( allocated(unique_atom_symequiv) ) then
        !
        ! Recursive deallocation of nested structures with allocatable
        ! components:
        !
        deallocate(unique_atom_symequiv, stat=status)
        ASSERT(status.eq.0)
     endif
   end subroutine unique_atom_symequiv_free
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_symequiv_write(io_unit)
     ! Purpose : writes out the sym.equivalent distance
     !           information in 'input-format' in close
     !           analogy to those routines designed by TB.
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     !------------ Declaration of formal parameters ---------------
     integer(kind=i4_kind),intent(in)    :: io_unit
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind)   :: ia,ib,k,status
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------


     write(io_unit, fmt=*, iostat=status)""
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status) &
          "# symmetry equivalent distances"
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status)""
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status) &
          "# Format: for each pair of unique atoms ua1 >= ua2:"
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status) &
          "#   first line: Number n not symmetry equivalent connection vectors"
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status) &
          "#   n following lines: index of equal atom of ua2 and wight"
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")
     write(io_unit, fmt=*, iostat=status)""
     if(status.ne.0) call error_handler &
          ("unique_atom_symequiv_write: write failed")

     do ia=1,N_unique_atoms
        do ib=1,ia
           write(io_unit, fmt=*, iostat=status) &
                "# between unique atom ",ia," and unique atom ",ib
           if(status.ne.0) call error_handler &
                ("unique_atom_symequiv_write: write failed")
           write(io_unit, fmt=*, iostat=status) &
                unique_atom_symequiv(ia,ib)%n_equiv_dist
           if(status.ne.0) call error_handler &
                ("unique_atom_symequiv_write: write failed")

           do k=1,unique_atom_symequiv(ia,ib)%n_equiv_dist
              write(io_unit, fmt=*, iostat=status) &
                   unique_atom_symequiv(ia,ib)%index_eq_atom2(k), &
                   unique_atom_symequiv(ia,ib)%weight(k)
              if(status.ne.0) call error_handler &
                   ("unique_atom_symequiv_write: write failed")
           enddo

        enddo
     enddo

     write(io_unit, fmt=*, iostat=status)""

   end subroutine unique_atom_symequiv_write
   !*************************************************************


  subroutine unique_atom_bcast()
    implicit none
    ! *** end of interface ***

  end subroutine unique_atom_bcast


   subroutine unique_atom_symadapt_free
     !
     !  Purpose: deallocates information about group structure
     !           and symmetry adaption
     !
     use symmetry_data_module, only: symmetry_data_close
     implicit none
     !** End of interface *****************************************

     !
     ! unique_atoms hold copies of components
     !
     !  %symadapt_partner
     !  %symadapt_spor_partner
     !
     ! It as copied from uatom_symmadapt module,
     ! so they also may need to free the storage!
     !

     !
     ! We rely here on recursive deallocation, all components
     ! of unique_atom_type are allocatable, so the statment
     !
     !          deallocate(unique_atoms)
     !
     ! in unique_atom_close() will recursively deallocate everything.
     !

     !
     ! Deallocate information about group structure
     !
     call symmetry_data_close()
   end subroutine unique_atom_symadapt_free
   !*************************************************************

   !*************************************************************
   subroutine read_shell(typ, L, uab, ua)
   !
   ! Purpose:   read   in   all   information  in   "uab"   and   does
   ! allocation. FIXME:  is overloaded with the  logic to autogenerate
   ! fit exponents. Should be moved elsewhere.
   !
   use input_module
   use unique_atom_module, only: unique_atom_basis_norm_ob, &
        unique_atom_basis_norm_ch
   implicit none
   character(len=2), intent(in) :: typ ! "ob", "ch" or "xc"
   integer(i4_kind), intent(in) :: L ! -1 for r2-basis, otherwise >= 0
   type(unique_atom_basis_type), intent(inout) :: uab
   type(unique_atom_type), intent(in), optional :: ua ! only used to autogenerate exponents
   !** End of interface *****************************************

   integer             ::  i,j,status,unit
   ! function have to be read in
   integer             :: N_exponents = 0
   integer             :: N_uncontracted_fcts = 0
   integer             :: N_contracted_fcts = 0
   integer             :: N_fcts = 0
   logical             :: automatic ! deduce fit exponents from AO basis
                                    ! obtained before. for r2 and s only
   namelist /unique_atom_basis/ N_exponents, N_uncontracted_fcts, &
   N_contracted_fcts, N_fcts, automatic

   integer             :: N_bas
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------

   if ( .not. input_line_is_namelist("unique_atom_basis") ) call input_error( &
        "read_shell: namelist  unique_atom_basis expected")

   unit = input_intermediate_unit()

   N_exponents         = df_N_exponents
   N_uncontracted_fcts = df_N_uncontracted_fcts
   N_contracted_fcts   = df_N_contracted_fcts
   N_fcts = df_N_fcts
   automatic           = df_automatic
   call input_read_to_intermediate()
   read(unit, nml=unique_atom_basis, iostat=status)
   if (status .gt. 0) call input_error( &
        "read_shell: namelist unique_atom_basis.")

   if (N_exponents < 0 ) call input_error( &
        "read_shell: namelist unique_atom_basis: N_exponents")

   if (N_uncontracted_fcts < 0 .or. N_uncontracted_fcts > N_exponents) &
        call input_error( &
        "read_shell: namelist unique_atom_basis: N_uncontracted_fcts")

   !
   ! Legacy automatic generation of fit basis is only supported for s-
   ! and r2 fit functions:
   !
   if (automatic .and. (typ == "ob" .or. L > 0)) then
      call input_error( &
           "read_shell: automatic makes only sense for s and r2 fitting functions")
   endif

   if (automatic) then
      ! this is used to get orbital exponents:
      ASSERT(present(ua))

      if (L == 0) then
         call unique_atom_basis_alloc(uab, ua%l_ob(0)%N_exponents, ua%l_ob(0)%N_exponents, 0, N_fcts)
         ! FIXME: for auto XC fit factor 2 is questionable:
         uab%exponents=2.0_r8_kind*ua%l_ob(0)%exponents
      else
         ASSERT(L==-1)
         call unique_atom_basis_alloc(uab, ua%l_ob(1)%N_exponents, ua%l_ob(1)%N_exponents, 0, N_fcts)
         ! FIXME: for auto XC fit factor 2 is questionable:
         uab%exponents=2.0_r8_kind*ua%l_ob(1)%exponents
      end if
   else
      call unique_atom_basis_alloc(uab, N_exponents, N_uncontracted_fcts, N_contracted_fcts, N_fcts)
      ! read exponents
      if ( N_exponents /= 0 ) then
         call input_read_to_intermediate
         read(unit, fmt=*, iostat=status) uab%exponents
         if (status .gt. 0) call input_error( &
              "read_shell: exponents.")

         ! loop over contractions
         do i=1,N_contracted_fcts
            ! read coefs
            call input_read_to_intermediate
            read(unit, fmt=*, iostat=status) ( uab%contractions(j,i), j=1,N_exponents )
            if (status .gt. 0) call input_error( &
                 "read_shell: contraction coefs.")
         enddo

         ! loop over linear combinations:
         N_bas = N_uncontracted_fcts + N_contracted_fcts
         do i = 1, N_fcts
            ! read coefs
            call input_read_to_intermediate
            read(unit, fmt=*, iostat=status) ( uab%functions(j, i), j = 1, N_bas )
            if (status .gt. 0) call input_error( &
                 "read_shell: linear combinaiton coefs.")
         enddo
      endif
   end if

   !
   ! Compute  legacy norms  which are  supposed to  correspond  to the
   ! local    normalization    convention    used   for    contraction
   ! coefficients. Norms  are used explicitly in grid  based code when
   ! computing  values  of  the   functions  at  a  space  point  and,
   ! implicitly, in the code for analytic integrals.
   !
   if (typ == "ob") then
      uab%norms = unique_atom_basis_norm_ob(L, uab%exponents)
   else if (typ == "ch" .or. typ == "xc") then
      !
      ! Here  norm_ch  must  also  handle  the case  of  L=-1  for  r2
      ! exponents:
      !
      uab%norms = unique_atom_basis_norm_ch(L, uab%exponents)
   else
      ABORT("no such typ")
   endif
   end subroutine read_shell
   !*************************************************************

   !*************************************************************
   subroutine unique_atom_write(iounit) ! was unique_atom_unique_write(...)
   !
   !  Purpose: writes name and Z and coordinates of first
   !           centers for each unique atoms
   use echo_input_module, only: start, real, flag, intg, word, stop, &
        echo_level_full, word_format
   use operations_module, only: operations_echo_input_level
   use atom_data_module, only: nuc_radius
   implicit none
   integer, intent(in) :: iounit
   !** End of interface *****************************************

   type(unique_atom_type), pointer :: ua
   integer             ::  i,k, status
   character(len=29)   ::  header

   external error_handler

   call start("UNIQUE_ATOM_NUMBER","UNIQUE_ATOM_UNIQUE_WRITE", &
        iounit,operations_echo_input_level)
   call intg("N_UNIQUE_ATOMS", N_unique_atoms, df_N_unique_atoms)
   call stop

   write(iounit, fmt='(A/)', iostat=status) " # << Geometry >>"
   if (status /= 0) call error_handler &
        ("UNIQUE_ATOM_UNIQUE_WRITE: writing of geometry banner failed")

   word_format = '("    ",a," = ",a14  :" # ",a)' ! including quotes

   ! loop over unique atoms
   do i=1,N_unique_atoms
      ua => unique_atoms(i)

      write(header,'("UNIQUE_ATOM # unique atom",i4)') i
      call start(header,"UNIQUE_ATOM_UNIQUE_WRITE", &
           iounit,operations_echo_input_level)
      call word("NAME         ", ua%name            , df_name         )
      call real("Z            ", ua%Z               , df_Z            )
      call intg("nuclear_spin ", ua%nuclear_spin    , df_nuclear_spin )
      call real("nuclear_magnetic_moment ", ua%nuclear_magnetic_moment , df_nuclear_magnetic_moment)
      call real("nuclear_radius ", ua%nuclear_radius , - nuc_radius(NINT(ua%Z)), &
           fmt='("    ",a," = ",es16.8:" # ",a)' )
      call intg("N_equal_atoms", ua%N_equal_atoms   , df_N_equal_atoms)
      call flag("FIXED        ", ua%moving_atom == 0, df_fixed        )
      call stop(empty_line=.false.)

      ! write atom coordinates of all equal atoms
      write(iounit, fmt='(3ES25.15)', iostat=status) (ua%position_first_ea(k), k=1,3)
      if (status .gt. 0) call error_handler("unique_atom_unique_write: error at positions.")
      write(iounit, fmt='()', iostat=status)
      if (status .gt. 0) call error_handler("unique_atom_unique_write: error after positions.")

   enddo

 end subroutine unique_atom_write
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_write_basis(iounit)
   !  Purpose: writes basis sets for all unique atoms
   !------------ Modules used -----------------------------------
   use echo_input_module
   use operations_module, only: operations_echo_input_level, operations_core_density
   !------------ Declaration of formal parameters ---------------
   integer, intent(in) :: iounit
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   type(unique_atom_type), pointer :: ua
   integer             ::  i,j,i_gc, status
   character(len=77)   ::  header
   logical             ::  pseudopot, core_dens
   integer             :: lmax_ob = 0
   integer             :: lmax_pseudo = 0
   integer             :: lmax_ch = 0
   integer             :: lmax_xc = 0
   !------------ Executable code --------------------------------
   external error_handler
   !------------ Executable code --------------------------------

   write(iounit, fmt='(A/)', iostat=status) " # << Basis Sets >>"
   if (status /= 0) call error_handler &
        ("UNIQUE_ATOM_WRITE_BASIS: writing of basis set banner failed")

   ! loop over unique atoms
   do i=1,N_unique_atoms
      ua => unique_atoms(i)
      pseudopot = ua%zc /= 0.0_r8_kind .and. .not.operations_core_density
      core_dens = ua%zc /= 0.0_r8_kind .and.      operations_core_density

      write(header,'("UNIQUE_ATOM_BASISSET # unique atom",i4)') i
      call start(header,"UNIQUE_ATOM_WRITE_BASIS", &
           iounit,operations_echo_input_level)
      call intg("LMAX_OB       ", ua%lmax_ob       , df_lmax_ob       )
      call intg("LMAX_PSEUDO   ", ua%lmax_pseudo   , df_lmax_pseudo   )
      call intg("LMAX_CH       ", ua%lmax_ch       , df_lmax_ch       )
      call intg("LMAX_XC       ", ua%lmax_xc       , df_lmax_xc       )
      call intg("N_GLOB_CONS_CH", ua%N_glob_cons_ch, df_N_glob_cons_ch)
      call intg("N_GLOB_CONS_XC", ua%N_glob_cons_xc, df_N_glob_cons_xc)
      call real("ZC            ", ua%ZC            , df_ZC, 2         )
      call stop()

      lmax_ob = ua%lmax_ob
      lmax_ch = ua%lmax_ch
      lmax_xc = ua%lmax_xc
      lmax_pseudo = ua%lmax_pseudo

      ! write l dependent basis information for orbitals
      do j=0,lmax_ob
         write(header,'( "UNIQUE_ATOM_BASIS # unique atom",i4,", orbital basis set : l = ",i1)') i,j
         call unique_atom_basis_write(iounit, ua%l_ob(j), header)
      enddo

      ! write pseudo potential data (if required)
      if (pseudopot) then
         ! first write the local pseudo potential data
         write(header, '("UNIQUE_ATOM_PSEUDOPOT # unique atom",i4, ",  pseudopot. representation: local")') i
         call unique_atom_pseudopot_write(iounit, ua%l_pseudopot(lmax_pseudo), header)
         ! then write the semi-local contributions
         do j=0,lmax_pseudo-1
            write(header, '("UNIQUE_ATOM_PSEUDOPOT # unique atom",i4, ",  pseudopot. representation: l = ", i1)') i, j
            call unique_atom_pseudopot_write(iounit, ua%l_pseudopot(j), header)
         enddo
      endif

#ifdef WITH_CORE_DENS
      ! write the core density basis set (if required)
      if (pseudopot .or. core_dens) then
         ! first write the r^2-type date
         write(header, '("UNIQUE_ATOM_CORE_DENSITY # unique atom",i4, ",  r**2-type core density")') i
         call unique_atom_core_density_write(iounit, ua%r2_core, header, basis_only=core_dens )
         ! then write the s-type date
         write(header, '("UNIQUE_ATOM_CORE_DENSITY # unique atom",i4, ",  s-type core density")') i
         call unique_atom_core_density_write(iounit, ua%s_core , header, basis_only=core_dens )
      endif
#endif

      ! write r**2 basis information for charge fitfcts
      write(header,'("UNIQUE_ATOM_BASIS # unique atom",i4, ", charge fit basis set : r**2")') i
      call unique_atom_basis_write(iounit, ua%r2_ch, header)

      ! write l dependent basis information for charge fitfcts
      do j=0,lmax_ch
         write(header,'("UNIQUE_ATOM_BASIS # unique atom",i4, ", charge fit basis set : l = ",i1)') i,j
         call unique_atom_basis_write(iounit, ua%l_ch(j), header)
      enddo

      ! write global contractions for charge fitfcts
      do i_gc = 1, ua%N_glob_cons_ch
         write(header,'("UNIQUE_ATOM_GLOB_CON # unique atom",i4, ", global charge fit contraction :",i4)') i,i_gc
         call unique_atom_glob_con_write(ua%glob_con_ch(i_gc),iounit,header)
      enddo

      ! write r**2 basis information for exchange fitfcts
      write(header,'("UNIQUE_ATOM_BASIS # unique atom",i4, ", exchange fit basis set : r**2")') i
      call unique_atom_basis_write(iounit, ua%r2_xc, header)

      ! write l dependent basis information for exchange fitfcts
      do j=0,lmax_xc
         write(header,'("UNIQUE_ATOM_BASIS # unique atom",i4, ", exchange fit basis set : l = ",i1)') i,j
         call unique_atom_basis_write(iounit, ua%l_xc(j), header)
      enddo

      ! write global contractions for exchange fitfcts
      do i_gc = 1, ua%N_glob_cons_xc
         write(header,'("UNIQUE_ATOM_GLOB_CON # unique atom",i4, ", global exchange fit contraction :",i4)') i,i_gc
         call unique_atom_glob_con_write(ua%glob_con_xc(i_gc),iounit,header)
      enddo

   enddo

   end subroutine unique_atom_write_basis
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_symadapt_write(output_unit)
   !  Purpose: writes information about symmetry adaption for all
   !           unique atoms in unique_atoms(:)
   !------------ Declaration of formal parameters ---------------
   integer, intent(in) :: output_unit
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   type(unique_atom_type), pointer :: ua
   type(unique_atom_symadapt_type), pointer :: sa
   type(unique_atom_partner_type), pointer :: sap
   integer             :: i,j,k,l,m,n,lmax,status
   integer             :: N_fcts
   namelist /unique_atom_indp_fct/  N_fcts
   namelist /unique_atom_symadapt/  N_fcts
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   DPRINT 'unique_atom_symadapt_write: entered'
   write(output_unit, fmt=*, iostat=status)
   if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
   write(output_unit, fmt=*, iostat=status)
   if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
   write(output_unit, fmt=*, iostat=status) "# symmetry adaption"
   if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
   write(output_unit, fmt=*, iostat=status)
   if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
    ! loop over unique atoms
   do i=1,N_unique_atoms
      DPRINT 'unique_atom_symadapt_write: i=',i
      write(output_unit, fmt=*, iostat=status)
      if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
      ua => unique_atoms(i)
      DPRINT 'Shape(ua%symadapt_partner)=',Shape(ua%symadapt_partner)
      ! loop over angular momentum
      lmax = ua%lmax_all
      do l = 0, lmax
         DPRINT 'unique_atom_symadapt_write: l=',l
         write(output_unit, fmt=*, iostat=status)
         if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
         ! loop over irreps
         do j = 1, symmetry_data_n_irreps()
            DPRINT 'unique_atom_symadapt_write: j=',j
            sap => ua%symadapt_partner(j,l)
            write(output_unit, fmt=*, iostat=status) "# symmetry adaption  ua ", i, &
                 " l ", l, " Irrep ", j
            if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
            N_fcts = sap%N_independent_fcts
            write(output_unit, nml=unique_atom_indp_fct, iostat=status)
            if ( status .ne. 0 ) call error_handler( &
                 "unique_atom_symadapt_write: namelist unique_atom_indp_fct")
            ! loop over partners
            do k = 1, symmetry_data_n_partners(j)
               DPRINT 'unique_atom_symadapt_write: k=',k
               ! loop over independent fcts
               do m = 1, sap%N_independent_fcts
                  DPRINT 'unique_atom_symadapt_write: m=',m
                  sa => sap%symadapt(m,k)
                  N_fcts = sa%N_fcts
                  write(output_unit, fmt=*, iostat=status) "# partner ", k, " ind. fct ", m
                  if (status .gt. 0) call error_handler("unique_atom_symadapt_write: comment")
                  write(output_unit, nml=unique_atom_symadapt, iostat=status)
                  if ( status .ne. 0 ) call error_handler( &
                       "unique_atom_symadapt_write: namelist unique_atom_symadapt")
                  do n = 1, N_fcts
                     DPRINT 'unique_atom_symadapt_write: n=',n
                     write(output_unit, fmt='(2I4,ES25.15)', iostat=status) &
                          sa%I_equal_atom(n), sa%m(n),  sa%c(n)
                     if ( status .ne. 0 ) call error_handler( &
                          "unique_atom_symadapt_write: write failed")
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   end subroutine unique_atom_symadapt_write
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_basis_write(iounit, l, header)
   !
   ! Purpose: write all information in l
   !
   use echo_input_module, only: start, real, flag, intg, strng, stop, &
        echo_level_full
   use operations_module, only: operations_echo_input_level
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_basis_type), intent(inout)  :: l
   integer, intent(in) :: iounit
   character(len=*) :: header
   !** End of interface *****************************************

   integer             ::  i,status
   integer             :: N_exponents = 0
   integer             :: N_uncontracted_fcts = 0
   integer             :: N_contracted_fcts = 0
   integer             :: N_fcts
   logical             :: automatic = .false.
   !------------ Declaration of subroutines used ----------------
   external error_handler

   N_exponents = l%N_exponents
   N_uncontracted_fcts = l%N_uncontracted_fcts
   N_contracted_fcts = l%N_contracted_fcts
   N_fcts = size(l%functions, 2)

   call start(header,"UNIQUE_ATOM_BASIS_WRITE", &
        iounit,operations_echo_input_level)
   call intg("N_EXPONENTS        ",N_exponents        ,df_N_exponents        )
   call intg("N_UNCONTRACTED_FCTS",N_uncontracted_fcts,df_N_uncontracted_fcts)
   call intg("N_CONTRACTED_FCTS  ",N_contracted_fcts  ,df_N_contracted_fcts  )
   call intg("N_FCTS             ",N_fcts             ,df_N_fcts )

   ! clear this issue: where should it come from?
   ASSERT(.not.automatic)
   call flag("AUTOMATIC          ",automatic          ,df_automatic          )
   call stop(empty_line=.false.)

   ! write exponents
   write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%exponents
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at exponents.")

   ! loop over contractions
   do i=1,N_contracted_fcts
      write(iounit, fmt='("  # Contraction",i4)', iostat=status) i
      if (status .gt. 0) call error_handler( &
           "unique_atom_basis_write: write error at contraction headers.")
      write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%contractions(:,i)
      if (status .gt. 0) call error_handler( &
           "unique_atom_basis_write: write error at contractions.")
   enddo

   ! loop over combinations
   do i = 1, size(l%functions, 2)
      write(iounit, fmt='("  # Combination",i4)', iostat=status) i
      if (status .gt. 0) call error_handler( &
           "unique_atom_basis_write: write error at combination headers.")
      write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%functions(:, i)
      if (status .gt. 0) call error_handler( &
           "unique_atom_basis_write: write error at combination.")
   enddo

   ! terminate writing
   write(iounit, fmt='()', iostat=status)
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at empty line.")

   end subroutine unique_atom_basis_write
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_pseudopot_write(iounit, l, header)
   !
   ! Purpose: write all information in l
   !
   use echo_input_module, only: start, real, flag, intg, strng, stop, &
        echo_level_full
   use operations_module, only: operations_echo_input_level
   implicit none
   type(unique_atom_pseudopot_type), intent(inout)  :: l
   integer, intent(in) :: iounit
   character(len=*) :: header
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer             ::  status
   integer             :: N_exponents = 0
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   N_exponents = l%N_exponents

   call start(header,"UNIQUE_ATOM_PSEUDOPOT_WRITE", &
        iounit,operations_echo_input_level)
   call intg("N_EXPONENTS",N_exponents,df_N_exponents)
   call stop(empty_line=.false.)

   ! write powers of prefactors
   write(iounit, fmt='("  # powers of the prefactors")')
   write(iounit, fmt='(15I5:" %")', iostat=status) l%powers
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at powers.")
   ! write exponents
   write(iounit, fmt='("  # exponents")')
   write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%exponents
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at exponents.")
   ! write coefficients
   write(iounit, fmt='("  # coefficients")')
   write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%coefficients
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at coefficients.")

   ! terminate writing
   write(iounit, fmt='()', iostat=status)
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at empty line.")

   end subroutine unique_atom_pseudopot_write
   !*************************************************************


   !*************************************************************
#ifdef WITH_CORE_DENS
   subroutine unique_atom_core_density_write(iounit, l, header, basis_only)
   !
   ! Purpose: write all information in l
   !
   use echo_input_module, only: start, real, flag, intg, strng, stop, &
        echo_level_full
   use operations_module, only: operations_echo_input_level
   implicit none
   type(unique_atom_atomic_dens_type), intent(inout)  :: l
   integer, intent(in) :: iounit
   character(len=*) :: header
   logical, optional :: basis_only
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer             ::  status
   logical             ::  write_contractions
   integer             :: N_exponents = 0
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if (present(basis_only)) then
      write_contractions = .not. basis_only
   else
      write_contractions = .true.
   endif

   N_exponents = l%N_exponents

   call start(header,"UNIQUE_ATOM_CORE_DENSITY_WRITE", &
        iounit,operations_echo_input_level)
   call intg("N_EXPONENTS",N_exponents,df_N_exponents)
   call stop(empty_line=.false.)

   ! write exponents
   write(iounit, fmt='("  # exponents")')
   write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%exponents
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at exponents.")
   if (write_contractions) then
      ! write contractions
      write(iounit, fmt='("  # contraction")')
      write(iounit, fmt='(3ES25.15:"  %")', iostat=status) l%contractions
      if (status .gt. 0) call error_handler( &
           "unique_atom_basis_write: write error at contractions.")
   endif
   ! terminate writing
   write(iounit, fmt='()', iostat=status)
   if (status .gt. 0) call error_handler( &
        "unique_atom_basis_write: write error at empty line.")

   end subroutine unique_atom_core_density_write
#endif
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_setup()
     !
     !  Purpose: Does several setup task necessary between
     !           symmetry part and integral part:
     !  + calculates lmax of all unique atoms for ob and fitfcts
     !  + calculates norms of all basis functions
     !  + calculates dimensions of all irreps and dimension
     !    of fitfunction basis of ch and xc
     !  + elims empty irreps
     !
     ! Runs in parallel context on all workers.
     !
     use type_module, only: IK=>i4_kind
     use options_module, only: options_spin_orbit
     use symmetry_data_module, only:&
          & symmetry_data_n_irreps,symmetry_data_n_proj_irreps,&
          & symmetry_data_set,symmetry_data_elim_irreps,&
          & symmetry_data_elim_pirreps
     use uatom_symmadapt, only: uaSymm, uaSymmSpor, uas_Large, uas_Small
     use clebsch_gordan, only: clebsch_gordan_bcast, clebsch_gordan_eliminate
     use spin_orbit_module, only: is_on, op_FitTrafo
     implicit none
     !** End of interface ***************************************

     integer :: memstat
     integer(IK) :: nirr, i
     integer(IK), allocatable :: dims_vec(:), dims_proj(:)

     DPRINT MyID, "unique_atom_setup: entered"

     !
     ! We seem to use CG module below, so let it do its
     ! parallel setup first:
     !
     call clebsch_gordan_bcast(options_spin_orbit)

     !
     ! Broadacst the necessary minimum, let the slaves derive
     ! the rest:
     !
     call unique_atom_bcast()

     !
     ! Now compute the missing parts, MUSTDIE!
     !
     call unique_atom_calc_lmax()

     !
     ! moved from unique_atom_calc_irrep_dims:
     !
     call find_totalsymmetric_irrep()

     DPRINT MyID,'unique_atom_setup: shape(unique_atoms)=',shape(unique_atoms)
     DPRINT MyID,'unique_atom_setup: shape(uaSymm)=',shape(uaSymm)

     !
     ! FIXME: Zero as an irrep dimension is a perfectly valid number,
     !        why do we go through the pain of eliminating them?
     !

     !
     ! Compute dimensions of vector irreps (some may be empty)
     !
     nirr = symmetry_data_n_irreps()

     allocate(dims_vec(nirr), STAT=memstat)
     ASSERT(memstat==0)

     !
     ! Take the basis info from unique_atoms and symmetry info from uaSymm:
     !
     call unique_atom_calc_irrep_dims(unique_atoms, uaSymm, dims_vec, spor=.false.)

     DPRINT MyID,'unique_atom_setup: vec_dims=',dims_vec

     !
     ! This continues initialization of symmetry_data_module.
     !
     do i = 1, nirr
        call symmetry_data_set(index=i, dim=dims_vec(i))
     enddo

     !
     ! Eliminate empty irreps from uaSymm:
     !
     call unique_atom_elim_empty_syadapt(uaSymm, mask=(dims_vec.eq.0))

     DPRINT MyID,'unique_atom_setup: 1. call unique_atom_assign_symm(unique_atoms ...)'
     !
     ! Set %symadapt_partner component of unique_atoms:
     !
     do i = 1, n_unique_atoms
        call unique_atom_assign_symm(unique_atoms(i), uaSymm(i))
     enddo

     call symmetry_data_elim_irreps()

     DPRINT MyID,'unique_atom_setup: call clebsch_gordan_eliminate(SR)'
     call clebsch_gordan_eliminate(which_vec=(dims_vec.eq.0))
     DPRINT MyID,'unique_atom_setup: done clebsch_gordan_eliminate(...)'

     !---------------------------------

     !---------------------------------
     ! spin-orbit:
     if(options_spin_orbit)then

        nirr = symmetry_data_n_proj_irreps()

        allocate(dims_proj(nirr), STAT=memstat)
        ASSERT(memstat==0)

        call unique_atom_calc_irrep_dims(unique_atoms, uas_Large, dims_proj, spor=.true.)

        do i = 1, nirr
           call symmetry_data_set(index=i, dim_proj=dims_proj(i))
        enddo

        DPRINT MyID,'unique_atom_setup: eliminate Large...'
        call unique_atom_elim_empty_syadapt(uas_Large, mask=(dims_proj.eq.0))
        if(is_on(op_FitTrafo))then
           DPRINT MyID,'unique_atom_setup: eliminate Small...'
           call unique_atom_elim_empty_syadapt(uas_Small, mask=(dims_proj.eq.0))
        endif
        DPRINT MyID,'unique_atom_setup: 2. call unique_atom_assign_symm(unique_atoms ...)'
        !
        ! set %symadapt_spor_partner component of unique atoms:
        !
        ! uaSymmSpor => uas_Large
        ! load to unique_atoms:
        do i = 1, n_unique_atoms
           call unique_atom_assign_symm(unique_atoms(i), uaSymmSpor(i))
        enddo

        call symmetry_data_elim_pirreps()

        DPRINT MyID,'unique_atom_setup: call clebsch_gordan_eliminate(SO)'
        call clebsch_gordan_eliminate(which_vec=(dims_vec.eq.0), which_proj=(dims_proj.eq.0))
        DPRINT MyID,'unique_atom_setup: done clebsch_gordan_eliminate(...)'
     endif

     call unique_atom_grad_information()

     DPRINT MyID, "unique_atom_setup: exit"
    end subroutine unique_atom_setup
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_calc_lmax
     !  Purpose: calculates lmax of all unique atoms for ob and fitfcts
     !** End of interface ***************************************
     use operations_module, only: operations_gradients
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind)           :: i_ua
     type(unique_atom_type), pointer :: ua
     !------------ Executable code --------------------------------

     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        unique_atom_lmax_ob = max( unique_atom_lmax_ob, ua%lmax_ob )
        unique_atom_lmax_ch = max( unique_atom_lmax_ch, ua%lmax_ch )
        unique_atom_lmax_xc = max( unique_atom_lmax_xc, ua%lmax_xc )
        unique_atom_lmax_pseudo = &
        max(unique_atom_lmax_pseudo,ua%lmax_pseudo)
     enddo
     if (operations_gradients) unique_atom_lmax_all = 1
     unique_atom_lmax_all = max( unique_atom_lmax_all,&
                                 unique_atom_lmax_ob, &
                                 unique_atom_lmax_ch, &
                                 unique_atom_lmax_xc, &
                                 unique_atom_lmax_pseudo )
   end subroutine unique_atom_calc_lmax
   !*************************************************************

   !*************************************************************
   subroutine unique_atom_calc_irrep_dims(unique_atoms, uatoms, dims, spor)
     !
     !  Purpose: calculates dimensions of all irreps and dimension
     !
     ! Gets the basis size info from "unique_atoms" and symmetry info from
     ! "uatoms"
     !
     use uatom_symmadapt, only: uatom
     type(unique_atom_type), intent(in) :: unique_atoms(:) ! (n_unique_atoms)
     type(uatom), intent(in) :: uatoms(:) ! (n_unique_atoms)
     integer(i4_kind), intent(out) :: dims(:) ! (n_irrep)
     logical, intent(in) :: spor
     ! *** end of interface ***

     integer :: irrep, ua, L

     do irrep = 1, size(dims)
         dims(irrep) = 0
         do ua = 1, size(unique_atoms)
             do L = 0, unique_atoms(ua)%lmax_ob
                 dims(irrep) = dims(irrep) &
                     + nrad(ua, L) * nindep(ua, L, irrep, spor)
             enddo
         enddo
     enddo

     contains

      integer function nrad(ua, L)
        !
        ! Number of radial funcitons in a shell.
        !
        implicit none
        integer, intent(in) :: ua, L
        ! *** end of interface ***

        nrad = unique_atoms(ua)%l_ob(L)%N_uncontracted_fcts &
             + unique_atoms(ua)%l_ob(L)%N_contracted_fcts
      end function nrad

      integer function nindep(ua, L, irrep, spor)
        !
        ! Number of symmetry adapted combinations for
        ! an orbital of this shell.
        !
        implicit none
        integer, intent(in) :: ua, L, irrep
        logical, intent(in) :: spor
        ! *** end of interface ***

        if ( spor ) then
            nindep = uatoms(ua)%symadapt_spor_partner(irrep, L)%N_independent_fcts
        else
            nindep = uatoms(ua)%symadapt_partner(irrep, L)%N_independent_fcts
        endif
      end function nindep

   end subroutine unique_atom_calc_irrep_dims
   !*************************************************************

   subroutine unique_atom_elim_empty_syadapt(uas, mask)
     !
     !  Purpose: elliminates all irreps with index i such that
     !
     !                  mask(irr) == true
     !
     !           Is called after orbital and fitfunction basis
     !           and symadapt is read in.
     !
     use uatom_symmadapt, only: uatom, uatom_prune
     type(uatom), intent(inout) :: uas(:) ! (n_unique_atoms)
     logical, intent(in) :: mask(:) ! (n_irreps)
     ! *** end of interface ***

     integer(i4_kind) :: i_ua

     DPRINT   "unique_atom_elim_empty_syadapt: mask=", mask
     do i_ua = 1, size(uas)
        uas(i_ua) = uatom_prune(uas(i_ua), mask)
     enddo
     DPRINT   "unique_atom_elim_empty_syadapt: done"
   end subroutine unique_atom_elim_empty_syadapt

   !*************************************************************
   subroutine unique_atom_glob_con_read(c)
   !  Purpose: read in all information in c and does allocation
   use input_module
   !------------ Declaration of formal parameters ---------------
   type(unique_atom_glob_con_type), intent(inout)  :: c
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer             ::  i,status,unit
   integer             :: N_contributing_fcts
   namelist /unique_atom_glob_con/ N_contributing_fcts
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   if ( .not. input_line_is_namelist("unique_atom_glob_con") ) &
        call input_error("unique_atom_glob_con_read: &
        &namelist unique_atom_glob_con expected")

   unit = input_intermediate_unit()

   N_contributing_fcts = 0
   call input_read_to_intermediate()
   read(unit, nml=unique_atom_glob_con, iostat=status)
   if (status .gt. 0) call input_error( &
        "unique_atom_glob_con_read: namelist unique_atom_glob_con.")
   if (N_contributing_fcts .le. 0) call input_error( &
        "unique_atom_glob_con_read: namelist unique_atom_glob_con: N_contributing_fcts")
   call unique_atom_glob_con_alloc(c, N_contributing_fcts)
   ! read indices and coef
   do i=1, N_contributing_fcts
      call input_read_to_intermediate
      read(unit, fmt=*, iostat=status) &
           c%l(i), c%index_ind_fct(i), c%index_exp(i), c%coefs(i)
      if (status .gt. 0) call input_error( &
           "unique_atom_glob_con_read: contraction indices and coefs.")
   enddo

   end subroutine unique_atom_glob_con_read
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_glob_con_write(c,io_unit,header)
   !
   !  Purpose: writes all information in c to io_unit
   !
   use echo_input_module, only: start, real, flag, intg, strng, stop, &
        echo_level_full
   use operations_module, only: operations_echo_input_level
   implicit none
   type(unique_atom_glob_con_type), intent(in) :: c
   integer,                         intent(in) :: io_unit
   character(len=*),                intent(in) :: header
   !** End of interface *****************************************
   !------------ Declaration of local variables -----------------
   integer             ::  i,status
   integer             :: N_contributing_fcts
   !------------ Declaration of subroutines used ----------------
   external error_handler
   !------------ Executable code --------------------------------
   N_contributing_fcts = c%N_contributing_fcts

   call start(header,"UNIQUE_ATOM_GLOB_CON_WRITE", &
        io_unit,operations_echo_input_level)
   call intg("N_CONTRIBUTING_FCTS",N_contributing_fcts,df_N_contributing_fcts)
   call stop(empty_line=.false.)

   ! write indices and coef
   write(io_unit, fmt='(A)',iostat=status) &
        ' # l ind exp     coefficient          [ l = s(-1),r2(0),p(1),... ]'
   if (status .gt. 0) call error_handler( &
        "unique_atom_glob_con_write: contraction header.")
   do i=1, N_contributing_fcts
      write(io_unit, fmt='(3I4,ES25.15)', iostat=status) &
           c%l(i), c%index_ind_fct(i), c%index_exp(i), c%coefs(i)
      if (status .gt. 0) call error_handler( &
           "unique_atom_glob_con_write: contraction indices and coefs.")
   enddo
   write(io_unit, fmt='()', iostat=status)
   if (status .gt. 0) call error_handler( &
        "unique_atom_glob_con_write: empty line.")

   end subroutine unique_atom_glob_con_write
   !*************************************************************

   !*************************************************************
   subroutine unique_atom_gradinfo_dealloc()
     ! Purpose: deallocate the variable 'unique_atom_grad_info'
     !          and the pointer 'moving_unique_atom_index'
     ! subroutine called by: main_gradient
     implicit none
     !** End of interface *****************************************

     integer(i4_kind) :: i, alloc_stat

     !
     ! FIXME:   arrmat3  has  pointer   components.  Change   that  to
     ! allocatable  and  rely  on  recursive  deallocation  of  nested
     ! structures with allocatable components:
     !
     if (allocated (unique_atom_grad_info)) then
        do i = 1, size (unique_atom_grad_info)
           deallocate (unique_atom_grad_info(i) % m, stat=alloc_stat)
           ASSERT(alloc_stat==0)
        enddo
        deallocate (unique_atom_grad_info, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
     end if

     if (allocated (moving_unique_atom_index)) then
        deallocate(moving_unique_atom_index, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
     end if
   end subroutine unique_atom_gradinfo_dealloc
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_grad_information()
     !
     ! Purpose: change the variable
     !
     !   unique_atom(i) % symapdapt_partner(irrep, l) %
     !   symadapt(n_indep, n_partner) % {N_fcts, I_equal_atom, c, ...}
     !
     ! to a variable better suited for the gradient part where we loop
     ! over  equal  atoms.  (->   better  performance)  and  only  the
     ! non-fixed unique_atoms are considered.
     !
     ! Also loads the moving_unique_atom_index.
     !
     ! Subroutine called by: subroutine post_scf (if
     ! operations_gradients)
     !
     implicit none
     !** End of interface *****************************************

     type(unique_atom_type)         , pointer :: ua
     type(unique_atom_partner_type) , pointer :: sap
     type(unique_atom_symadapt_type), pointer :: sa
     type(arrmat3)                  , pointer :: gi
     integer(kind=i4_kind) :: i_ua,i_ma,i_if,i_cd,i_ea,ts, &
                              n_indep,n_equals,alloc_stat
     external error_handler
     !
     ! ----------- executable code -----------------------------
     !

     DPRINT   "unique_atom_grad_information: entered"

     allocate (moving_unique_atom_index(N_moving_unique_atoms), STAT=alloc_stat)
     ASSERT(alloc_stat==0)

     allocate (unique_atom_grad_info(N_moving_unique_atoms), STAT=alloc_stat)
     ASSERT(alloc_stat==0)

     ts = get_totalsymmetric_irrep()
     do i_ua = 1, N_unique_atoms
        ua => unique_atoms(i_ua)
        sap => ua % symadapt_partner(ts, 1) ! l = 1
        i_ma = ua % moving_atom

        ! FIXME:  there  is no  reason  why  symmetry  info should  be
        ! available only for  selected (moving) unique atoms. Somewhat
        ! confusingly   the   structure  unique_atom_grad_info(:)   is
        ! indexed by moving atom index. See how often this fires:
        ASSERT(i_ua==i_ma)

        if (i_ma <= 0) cycle

        moving_unique_atom_index(i_ma) = i_ua
        gi => unique_atom_grad_info(i_ma)

        n_indep = sap % N_independent_fcts
        n_equals = ua % N_equal_atoms

        allocate (gi % m(n_indep, 3, n_equals), STAT=alloc_stat)
        ASSERT(alloc_stat==0)

        gi % m(:, :, :) = 0.0

        do i_if = 1, n_indep ! loop over all sym. adapted gradients of ua
           sa => sap % symadapt(i_if, 1) ! first partner
           do i_cd = 1, sa % N_fcts ! loop over all contributing derivatives
              i_ea = sa % I_equal_atom(i_cd)

              select case (sa % m(i_cd))
              case (2); gi % m(i_if, 1, i_ea) = sa % c(i_cd)
              case (3); gi % m(i_if, 2, i_ea) = sa % c(i_cd)
              case (1); gi % m(i_if, 3, i_ea) = sa % c(i_cd)
              case default
                 call error_handler("unique_atom_grad_information: sth. wrong")
              end select

           enddo
        enddo
     enddo
   end subroutine unique_atom_grad_information
   !*************************************************************


   !*************************************************************
   subroutine unique_atom_make_gx (iloop)
     !
     ! Purpose: write a template for gx  file the gx file is used as a
     ! input  for  the  geometry  optimizer.   Subroutine  called  by:
     ! main_gradient().
     !
     use iounitadmin_module, only: openget_iounit, returnclose_iounit
     use filename_module, only: inpfile
     use operations_module, only: operations_gx_highprec
     use operations_module, only: OPERATIONS_GX_EPEFORMAT
     implicit none
     integer (i4_kind), intent (in) :: iloop
     !** End of interface *****************************************

     integer (i4_kind) :: i_ua, i_ea, io_u, counter
     integer (i4_kind), parameter :: deck(6) = 0
     external error_handler

     ! Creates a file in the input directory:
     io_u = openget_iounit (file= trim (inpfile('gxfile')), &
          status='new', form='formatted')

     counter = 0
     do i_ua = 1, n_unique_atoms
        do i_ea = 1, unique_atoms(i_ua)%n_equal_atoms
           counter = counter + 1
           if (operations_gx_highprec) then
              if (operations_gx_epeformat) then
                 write (io_u, '(F5.2,3(2x,F21.12),2I3,2X,3I3,2X,3I3,i5)')  &
                      unique_atoms(i_ua)%z, unique_atoms(i_ua)%position(:, i_ea),&
                      i_ua, counter, deck, -1
              else
                 write (io_u, '(F5.2,3(2x,F21.12),2I3,2X,3I3,2X,3I3)') &
                      unique_atoms(i_ua)%z, unique_atoms(i_ua)%position(:, i_ea),&
                      i_ua, counter, deck
              end if

           else
              if (OPERATIONS_GX_EPEFORMAT) then
                 write (io_u, '(F5.2,3(F15.7),2I4,i5,2I3,i5,2I3,i5)') &
                      unique_atoms(i_ua)%z, unique_atoms(i_ua)%position(:, i_ea), &
                      i_ua, counter, deck, -1
              else
                 write (io_u, '(F5.2,3(2x,F13.7),2I3,2X,3I3,2X,3I3)') &
                      unique_atoms(i_ua)%z, unique_atoms(i_ua)%position(:, i_ea), &
                      i_ua, counter, deck
              endif
           endif
        end do
     end do
     write (io_u, '(F5.1)') real (-iloop, kind= r8_kind)

     call returnclose_iounit (io_u)
   end subroutine unique_atom_make_gx
   !*************************************************************

   subroutine unique_atom_pseudopot_alloc(l,N_exponents)
   !  Purpose: allocates arrays contained in l and initialises l
   ! ----- Declaration of formal parameters ------------------------------------
   type(unique_atom_pseudopot_type),  intent(inout)   :: l
   integer(kind=i4_kind),             intent(in)      :: N_exponents
   ! ==== end of interface =====================================================
   ! ---- Declaration of local variables ---------------------------------------
   integer(kind=i4_kind)                              :: status
   ! ---- Declaration of subroutines used --------------------------------------
   external error_handler
   ! ---- Executable codes -----------------------------------------------------

   l%N_exponents = N_exponents

   allocate(l%powers(N_exponents),l%exponents(N_exponents),&
            l%coefficients(N_exponents), stat=status )
   if(status.ne.0) call error_handler( &
   "unique_atom_pseudopot_alloc: allocate of exponents,coefficients,powers failed")

   end subroutine unique_atom_pseudopot_alloc
   !*************************************************************

#ifdef WITH_CORE_DENS
   subroutine unique_atom_atomic_dens_alloc(l,N_exponents)
   !  Purpose: allocates arrays contained in l and initialises l
   ! ----- Declaration of formal parameters ------------------------------------
   type(unique_atom_atomic_dens_type), intent(inout)   :: l
   integer(kind=i4_kind)             , intent(in   )   :: N_exponents
   ! ==== end of interface =====================================================
   ! ---- Declaration of local variables ---------------------------------------
   integer(kind=i4_kind)                                 :: status
   ! ---- Declaration of subroutines used --------------------------------------
   external error_handler
   ! ---- Executable codes -----------------------------------------------------

   l%N_exponents   = N_exponents

   allocate(l%exponents(N_exponents), stat=status )
   if(status.ne.0) call error_handler( &
        "unique_atom_atomic_dens_alloc: allocate of exponents failed")

   ! single contraction so far:
   allocate(l%contractions(N_exponents), stat=status )
   if(status.ne.0) call error_handler( &
        "unique_atom_atomic_dens_alloc: allocate of contractions failed")
   end subroutine unique_atom_atomic_dens_alloc
#endif
   !*************************************************************

   subroutine unique_atom_pseudopot_read(l)
   !  Purpose: read in all information in l and does allocation
   use input_module
   !---------------------- Declaration of formal parameters --------------------
   type(unique_atom_pseudopot_type), intent(inout) :: l
   !--- end of interface -------------------------------------------------------
   !--- declaration of local variables -----------------------------------------
   integer(kind=i4_kind)                              :: status,unit
   integer                                 :: N_exponents = 0
   namelist /unique_atom_pseudopot/N_exponents
   !--- declaration of subroutines used ----------------------------------------
   external error_handler
   !--- executable codes -------------------------------------------------------

   if(.not.input_line_is_namelist("unique_atom_pseudopot") ) call &
    input_error("unique_atom_pseudopot_read:&
    & namelist unique_atom_pseudopot expected")
    unit = input_intermediate_unit()

   N_exponents = df_N_exponents
   call input_read_to_intermediate()
   read(unit, nml=unique_atom_pseudopot, iostat=status)
   if (status.gt.0) call input_error( &
        "unique_atom_pseudopot_read: namelist unique_atom_pseudopot.")
   if (N_exponents < 0) call input_error( &
        "unique_atom_pseudopot_read: namelist unique_atom_pseudopot:&
         & N_exponents")
   call unique_atom_pseudopot_alloc(l,N_exponents)

   ! FIXME: why doesnt it work for N_exponents == 0?
   if (N_exponents > 0) then
      !
      ! Read power prefactors, Gaussian exponents, and coefficients:
      !
      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) l%powers
      ASSERT(status==0)

      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) l%exponents
      ASSERT(status==0)

      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) l%coefficients
      ASSERT(status==0)
   endif

   end subroutine unique_atom_pseudopot_read
   !*************************************************************

#ifdef WITH_CORE_DENS
   subroutine unique_atom_core_density_read(l)
   !  Purpose: read in all information in l and does allocation
   use input_module
   !---------------------- Declaration of formal parameters --------------------
   type(unique_atom_atomic_dens_type), intent(inout) :: l
   !--- end of interface -------------------------------------------------------
   !--- declaration of local variables -----------------------------------------
   integer(kind=i4_kind)                              :: status,unit
   integer                                            :: N_exponents = 0
   namelist /unique_atom_core_density/N_exponents
   !--- declaration of subroutines used ----------------------------------------
   external error_handler
   !--- executable codes -------------------------------------------------------

   if(.not.input_line_is_namelist("unique_atom_core_density") ) call &
    input_error("unique_atom_core_density_read:&
    & namelist  unique_atom_core_density expected")
    unit = input_intermediate_unit()

   N_exponents = df_N_exponents_core
   call input_read_to_intermediate()
   read(unit, nml=unique_atom_core_density, iostat=status)
   if (status.gt.0) call input_error( &
        "unique_atom_core_density_read: namelist unique_atom_core_density.")

   if (N_exponents<0) call input_error( &
        "unique_atom_core_density_read: namelist unique_atom_core_density: &
        &N_exponents<0")

   call unique_atom_atomic_dens_alloc(l,N_exponents)

   ! read Gaussian exponents and coefficients
   call input_read_to_intermediate()
   read(unit, fmt=*, iostat=status) l%exponents
   if(status.gt.0) call input_error( &
        "unique_atom_core_density_read: exponents.")

   call input_read_to_intermediate()
   read(unit, fmt=*, iostat=status) l%contractions
   if(status.gt.0) call input_error( &
        "unique_atom_core_density_read: contractions.")

   end subroutine unique_atom_core_density_read
#endif
   subroutine compat_core_density_skip(l)
      ! skips /UNIQUE_ATOM_CORE_DENSITY/ namelist and data
      use input_module
      implicit none
      type(unique_atom_atomic_dens_type), intent(inout) :: l
      ! *** end of interface ***

      integer(kind=i4_kind) :: status,unit
      integer               :: N_exponents = 0
      namelist /unique_atom_core_density/ N_exponents

      integer(i4_kind) :: i
      real(r8_kind)    :: dumb

      if(.not.input_line_is_namelist("unique_atom_core_density") )then
         ! nothing to be done, silently ignore ...
         l%n_exponents=N_exponents
         RETURN
      endif

      unit = input_intermediate_unit()

      ! read namelist namelist itself:
      call input_read_to_intermediate()
      read(unit, nml=unique_atom_core_density, iostat=status)
      if (status.gt.0) call input_error( &
           "compat_core_density_skip: namelist unique_atom_core_density.")

      ! read exponents and coefficients
      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) (dumb,i=1,N_exponents)
      if(status.gt.0) call input_error( &
           "compat_core_density_skip: exponents.")

      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) (dumb,i=1,N_exponents)
      if(status.gt.0) call input_error( &
           "compat_core_density_skip: contractions.")
   end subroutine compat_core_density_skip
   !***********************************************************************

  !*****************************************************************************
  subroutine unique_atom_pseudopot_bcast(l)
    !  purpose : broadcasting all information in l and allocation for slaves
    !---------------------- modules used ---------------------------------------
    use comm, only: comm_bcast                                                 &
                  , comm_rank
    !----- Declaration of formal parameters ------------------------------------
    type(unique_atom_pseudopot_type),  intent(inout) :: l
    !----- Declaration of local variables --------------------------------------
    integer                                          :: N_exponents = 0
    !----- Executable code -----------------------------------------------------
    !
    call comm_bcast( n_exponents    )
    !
    if( comm_rank() /= 0 ) then
      call unique_atom_pseudopot_alloc(l,n_exponents)
    endif
    !
    ! broadcast symmetry index of Gaussians
    call comm_bcast( l%powers       )
    !
    ! broadcast exponents for Gaussians
    call comm_bcast( l%exponents    )
    !
    ! broadcast coefficients for Gaussians
    call comm_bcast( l%coefficients )
    !
  end subroutine unique_atom_pseudopot_bcast
  !*****************************************************************************

   !*************************************************************
   subroutine unique_atom_close()
     ! Purpose: reset the (private) variable 'int_done'
     !          at the end of the integralpart, such that
     !          unique_atom_calc_symadapt can be called again
     !          at the next run of the integral part
     !
     ! Runs in parallel context on all workers.
     !
     use uatom_symmadapt, only: uatom_symmadapt_done
     use clebsch_gordan, only: clebsch_gordan_close
     implicit none
     !** End of interface *****************************************

     integer(kind=i4_kind) :: alloc_stat

     !
     ! calls call symmetry_data_close()
     !
     call unique_atom_symadapt_free()

     !
     ! Deallocate
     !
     !          unique_atom_symequiv(:,:)
     !
     call unique_atom_symequiv_free()

     !
     ! Recursive deallocation of nested structures with allocatable
     ! components:
     !
     if(allocated(unique_atoms)) then
        deallocate(unique_atoms, stat=alloc_stat)
        ASSERT(alloc_stat.eq.0)
     end if

     !
     ! Deallocate
     !
     !          unique_atom_grad_info(:)
     !          moving_unique_atom_index(:)
     !
     call unique_atom_gradinfo_dealloc()

     !
     ! FIXME: what about resetting these from unique_atom_module?
     !
     ! N_unique_atoms = 0
     ! N_moving_unique_atoms
     ! unique_atoms_eperef(:)
     ! unique_atom_lmax_all = -1, &
     ! unique_atom_lmax_ob = -1, &
     ! unique_atom_lmax_ch = -1, &
     ! unique_atom_lmax_xc = -1, &
     ! unique_atom_lmax_pseudo=-1
     ! unique_atom_iwork
     ! pseudopot_present  = .false.
     ! core_density_setup = .false.

     !
     ! Intiate finalization in modules that we were using,
     ! see matching code in unique_atom_setup() ...
     !

     !
     ! Symmetry info is stored also in uatoms:
     !
     call  uatom_symmadapt_done()

     !
     ! Clebsch-Gordan coefficients are of no use:
     !
     call clebsch_gordan_close()
   end subroutine unique_atom_close
   !*************************************************************


#ifdef WITH_CORE_DENS
   !*************************************************************
   subroutine unique_atom_link_core_fitfct(core,fitfct)
   !  Purpose: establishes a link between core density fitting functions and
   !           the corresponding charge fitting functions
   !---------------------- Declaration of formal parameters --------------------
   type(unique_atom_atomic_dens_type), intent(inout) :: core
   type(unique_atom_basis_type)      , intent(in   ) :: fitfct
   !--- end of interface -------------------------------------------------------
   !--- declaration of local variables -----------------------------------------
   integer(kind=i4_kind) :: i_cexp, i_fexp, N_exponents, indfct
   real   (kind=r8_kind) :: target, exponent, exptol, eps = 1.0E-5_r8_kind
   !--- declaration of subroutines used ----------------------------------------
   external error_handler
   !--- executable codes -------------------------------------------------------
   N_exponents = fitfct%N_exponents
   do i_cexp=1,core%N_exponents
      target = core%exponents(i_cexp)
      exptol = target * eps
      indfct = 0
      fit: do i_fexp=1,fitfct%N_uncontracted_fcts
         exponent = fitfct%exponents(i_fexp)
         if ( abs(exponent-target) < exptol ) then
            indfct = i_fexp
            exit fit
         endif
      end do fit
      if (indfct == 0) call error_handler( &
           "unique_atom_read_basis: core fit basis must be a subset of uncontr. &
           &charge fit basis")
      core%linked_fitfcts(i_cexp) = indfct
   end do

   end subroutine unique_atom_link_core_fitfct
   !***********************************************************************
#endif

   subroutine unique_atom_assign_symm(lhs, rhs)
     !
     ! Sets the symmetry information from specialized
     ! uatom struct to a all-inclusive unique_atom_type
     !
     ! Affects only %symadapt_partner or %symadapt_spor_partner
     !
     use uatom_symmadapt, only: uatom
     use unique_atom_module, only: unique_atom_type
     implicit none
     type(unique_atom_type), intent(inout) :: lhs ! %symadapt*partner is intent(out)
     type(uatom), intent(in) :: rhs
     ! *** end of interface **

     integer(i4_kind) :: shp(2), memstat

     DPRINT 'unique_atom_assign_symm: entered'

     ASSERT(rhs%lmax==lhs%lmax_all)

     if (allocated(rhs%symadapt_partner)) then
        !
        ! STANDARD SCF:
        !

        !
        ! Array shape is (n_irrep, 0:lmax):
        !
        shp = shape(rhs%symadapt_partner)
        DPRINT   "allocated(lhs%symadapt_partner)=", allocated(lhs%symadapt_partner)
        allocate(lhs%symadapt_partner(shp(1), 0:shp(2)-1), stat=memstat)
        ASSERT(memstat.eq.0)

        ! actual copy:
        lhs%symadapt_partner = rhs%symadapt_partner
     endif

     if ( allocated(rhs%symadapt_spor_partner) ) then
        !
        ! SPIN-ORBIT:
        !

        !
        ! Array shape is (n_irrep, 0:lmax):
        !
        shp = shape(rhs%symadapt_spor_partner)
        allocate(lhs%symadapt_spor_partner(shp(1), 0:shp(2)-1), stat=memstat)
        ASSERT(memstat.eq.0)

        ! actual copy:
        lhs%symadapt_spor_partner = rhs%symadapt_spor_partner
     endif
   end subroutine unique_atom_assign_symm

!--------------- End of module ----------------------------------
end module unique_atom_methods

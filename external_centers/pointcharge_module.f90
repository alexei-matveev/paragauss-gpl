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
module  pointcharge_module
!---------------------------------------------------------------
!
!  Purpose: implements point charges, i. e. centers that
!           hold a charge but no basis functions.
!
!  Like atoms, point charges are grouped in collections of
!  symmetry equivalent point charges.
!
!  Module called by: data used in integral and symmetry part
!
!  Author: TB
!  Date:   5/97
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

#include "def.h"
use type_module ! type specification parameters
use comm_module
use datatype
use unique_atom_module,    only: unique_atom_type                              &
                               , unique_atom_partner_type                      &
                               , unique_atom_symadapt_type                     &
                               , unique_atom_sa_int_type                       &
                               , unique_atoms                                  &
                               , pseudopot_present                             &
                               , N_unique_atoms                                &
                               , unique_atom_iwork
#ifdef WITH_EFP
use qmmm_interface_module, only: efp
#endif
implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================


!!------------ Declaration of types ------------------------------
!type, public ::  pointcharge_type
!   character(len=12)                      :: name
!       ! name of unique point charge
!   real(kind=r8_kind)                :: Z
!       ! charge of point charge
!   integer                           :: N_equal_charges
!       ! Number of partners of unique point charge
!   real(kind=r8_kind), pointer       :: position(:,:)
!       !position(3,N_equal_charges)
!       ! positions of partners ("equal charges") of unique charge
!   real(kind=r8_kind)                :: position_first_ec(3)
!       ! positions of first partner ("equal charge") of unique charge
!end type pointcharge_type
!!!VN in request to numerous requestes of user this difinition above is moved
!!! to datatype module
!AM: make it visible for later use:
!!!public pointcharge_type

!------------ Declaration of constants and variables ------------
integer(i4_kind), public :: pointcharge_N ! modified from point_dqo_module.f90
integer(i4_kind), public, protected :: n_timps
type(pointcharge_type), pointer, public  :: pointcharge_array(:)
                ! pointcharge_array(pointcharge_number)
type(unique_atom_type), pointer, public  :: unique_timps(:)
type(unique_atom_type), pointer, public  :: unique_pc(:)
                ! unique_timps(n_unique_timps)
integer(kind=i4_kind), public           :: N_moving_unique_timps = -1
   ! Number of non-fixed unique atoms

type(arrmat3),pointer,public :: unique_timp_grad_info(:)
   ! unique_timp_grad_info(N_moving_unique_timps)%m(n_indep,3,n_equals)
   !  secound dimension for magnetical quantum numbers:  1: m=2, 2: m=3, 3: m=1
   ! Equivalent to variable
   !  unique_timp(i)%symapdapt_partner(irrep,l)%symadapt(n_indep,n_partner)%
   !    N_fcts
   !    I_equal_atom
   !    c ...
   ! with i_irrep = 1, l = 1
   ! but better suited for the gradient part where
   ! we loop over equal timps. (-> better performance)
type(arrmat3),pointer,public :: unique_pc_grad_info(:)

integer(kind=i4_kind),allocatable,public :: moving_unique_timp_index(:)
   ! moving_unique_timp_index(N_moving_unique_timps)
integer(i4_kind),allocatable,public :: unique_pc_index(:)

logical, public :: moving_pc=.false.
integer(i4_kind), public :: totsym_grad_pc_length=0
real(r8_kind),public,allocatable,target :: gradient_pc_totalsym(:)
real(r8_kind),public,allocatable,target :: torque_pc_totalsym(:)
type(arrmat2),public,allocatable :: gradient_pc_cartesian(:)
type(arrmat2),public,allocatable :: torque_pc_cartesian(:)

logical, public :: print_pc_grad=.false.
#ifdef WITH_EFP
real(r8_kind), public, allocatable :: rcm(:,:) ! mass center of efragment (better if
                                               ! define in efp_module
logical, public :: print_energy ! to print out QM-EFP and EFP-EFP energies (again not
                                ! very suitable place)
logical, public :: no_fixed, qm_fixed, efp_fixed
#endif
integer(i4_kind),public :: present_X_centers

!------------ public functions and subroutines ------------------
public :: pointcharge_read
public :: pointcharge_write
public :: init_pc_array
public :: pointcharge_bcast
public :: pointcharge_close
public :: init_pointcharges_grads
public :: pc_grads_shutdown
public :: pc_grad_cart_write
public :: totsym_PC_grad_unpack
public :: totsym_PC_grad_pack
public :: transform_PC_grad_to_cart
public :: calc_PC_grads
public :: calc_nuc_pc_energy
public :: unique_timp_grad_information
public :: unique_timp_symadapt_bcast

!================================================================
! End of public interface of module
!================================================================


!------------ Declaration of constants and variables ----
! namelists for reading and writing

namelist /pointcharge_number/ pointcharge_N, n_timps

character(len=12), private                   :: name
real(kind=r8_kind), private             :: Z, ZC
integer, private                        :: N_equal_charges
integer, private                        :: N_equal_timps=0
integer, private                        :: lmax_pseudo
logical, private                        :: fixed=.false.
real(kind=r8_kind), private             :: c = 0.0_r8_kind
real(kind=r8_kind), private             :: a = 0.0_r8_kind
namelist /pointcharge/ Z, N_equal_charges, name
namelist /timp/ lmax_pseudo, n_equal_timps, name, z, zc, fixed

!------------ Definition of Default Input Values -------
integer(kind=i4_kind), private :: df_pointcharge_N   = 0, &
                                  df_N_equal_charges = 0, &
                                  df_n_timps = 0, &
                                  df_n_equal_timps =0 , &
                                  df_lmax_pseudo = -1
real(kind=r8_kind), private    :: df_Z = 1.0_r8_kind
real(kind=r8_kind), private    :: df_ZC = 1.0_r8_kind
character(len=12), private          :: df_name = "            "
logical           , private :: df_fixed               = .false.
real(kind=r8_kind), private    :: df_c = 0.0_r8_kind
real(kind=r8_kind), private    :: df_a = 0.0_r8_kind
!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


   !*************************************************************
   subroutine pointcharge_read()
   !  Purpose: reads in names and Z of unique charges
   !           and coordinates of first centers ("equal charges")
   ! called by read_input()
   !** End of interface *****************************************
   !------------ Modules used ------------------- ---------------
   use input_module
   use operations_module
   use iounitadmin_module
   use filename_module
#ifdef WITH_EPE
   use ewaldpc_module
#endif
   use unique_atom_methods
   implicit none
   !------------ Declaration of local variables -----------------
   type(pointcharge_type), pointer :: pc
   type(unique_atom_type), pointer :: ut
   integer(kind=i4_kind)           :: i_pc, k, status, i_timp, l_pseudo, unit
   integer (i4_kind) :: io_u
#ifdef WITH_EPE
   integer (i4_kind) :: counter_equal, i_ua, io_gxepe
#endif
   real(kind=r8_kind)              :: iwork
   !------------ Executable code --------------------------------
   unit = input_intermediate_unit()

   pointcharge_N = df_pointcharge_N
   n_timps= df_n_timps
   if ( input_line_is_namelist("pointcharge_number") ) then
      call input_read_to_intermediate()
      read(unit, nml=pointcharge_number, iostat=status)
      if (status .ne. 0) call input_error( &
           "pointcharge_read: namelist pointcharge_number.")
   endif
   if ( pointcharge_N .lt. 0 ) call input_error( &
        "pointcharge_read: namelist pointcharge_number: pointcharge_N")
   if ( pointcharge_N .gt. 0 .and. .not. operations_symm ) call input_error( &
        "pointcharge_read: point charges require run of new symmetry part.")
   if ( n_timps .lt. 0 ) call input_error( &
        "pointcharge_read: namelist pointcharge_number: pointcharge_N")
   if ( n_timps .gt. 0 .and. .not. operations_symm ) call input_error( &
        "pointcharge_read: point charges require run of new symmetry part.")

   ! Zero-sized array are legal:
   ASSERT(pointcharge_N>=0)
   allocate (pointcharge_array(pointcharge_N), stat=status)
   ASSERT(status==0)

!!$        if((.not.options_spin_orbit).and.pointcharge_N>0) pseudopot_present=.true.
        ! results in problems for spin orbit run
        ! becouse of different value on master and slaves
        ! but required for regular run with pseudopotentials
   if(n_timps > 0) then
      allocate( unique_timps(n_timps), stat=status )
      if (status .ne. 0) call error_handler( &
           "pointcharge_read: allocate timps failed.")
   end if

   if(n_timps>0) pseudopot_present=.true.
   ! loop over unique timps
   N_moving_unique_timps = 0
   do i_timp=1,n_timps
      if ( .not. input_line_is_namelist("timp") ) call input_error( &
        "pointcharge_read: namelist pointcharge exspected")
      ut => unique_timps(i_timp)
      n_equal_timps = df_n_equal_timps
      name = df_name
      fixed         = df_fixed
      z = df_z
      zc = df_zc
      lmax_pseudo =df_lmax_pseudo
      call input_read_to_intermediate()
      read(unit, nml=timp, iostat=status)
      if (status .ne. 0) call input_error( &
           "pointcharge_read: namelist timp")
      if(z < 0.0)  call input_error( &
           "pointcharge_read: z lower than 0")
      if(zc < 0.0)  call input_error( &
           "pointcharge_read: zc lower than 0")
      if(zc > z) call input_error( &
           "pointcharge_read: zc greater than z")
      if (N_equal_timps .le. 0) call input_error( &
           "pointcharge_read: namelist timp: N_equal_timps")
      if (lmax_pseudo .le. 0) call input_error( &
           "pointcharge_read: namelist timp: lmax_pseudo")
      ut%name = name
      ut%lmax_pseudo = lmax_pseudo
      ut%n_equal_atoms = n_equal_timps
      ut%z = z
      ut%zc = zc
      ! in epe calc grads on timps are required
      if (fixed) then
         ut%moving_atom = 0
      else
         N_moving_unique_timps = N_moving_unique_timps + 1
         ut%moving_atom = N_moving_unique_timps!!!+N_moving_unique_atoms
      endif
      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) (ut%position_first_ea(k), k=1,3)
      if (status .ne. 0) call input_error("pointcharge_read: position timp.")
   end do

   ! now we read pseudo potentials associated with timps
   do i_timp=1,n_timps
      ut => unique_timps(i_timp)
      call unique_atom_alloc(ut,lmax_pseudo=ut%lmax_pseudo)
      call unique_atom_pseudopot_read(ut%l_pseudopot(ut%lmax_pseudo))
      do l_pseudo=0,ut%lmax_pseudo-1
         call unique_atom_pseudopot_read(ut%l_pseudopot(l_pseudo))
      end do
   end do

   ! read coordinats of timps from gxfile

   if(n_moving_unique_timps.gt.0) then
      if (operations_geo_opt.or.operations_read_gx) then
#ifdef WITH_EPE
         inquire(file= trim(inpfile('epe.r')), exist=ex_gxepe)
         if(ex_gxepe) then
            call gxepe_allocate
            io_gxepe=get_iounit()
            open (io_gxepe,status='old',form='formatted',file=trim(inpfile('epe.r')))
         endif ! ex_gxepe

         io_u=get_iounit()
         open (io_u,status='old',form='formatted',file=&
              trim(inpfile('gxfile')))
         call read_gxtimps(io_u,io_gxepe) !coordinats of reg positions
#else
         ABORT('recompile -DWITH_EPE')
#endif
      end if

      read(io_u,*,iostat=status) iwork
      if (status .gt. 0) call error_handler("unique_atom_unique_read: reading iwork from gxfile")
      if(iwork>=0.0_r8_kind) &
           call error_handler("unique_atom_unique_read: non negative value for iwork")
      unique_atom_iwork=int(-iwork,i4_kind)
      close(io_u)
      call return_iounit(io_u)
#ifdef WITH_EPE
      if(ex_gxepe) then
         write(output_unit,*) 'gxepe array centers'
         do i_ua=1,N_unique_atoms
            do counter_equal=1,unique_atoms(i_ua)%n_equal_atoms
               write(output_unit,*) gxepe_array(i_ua)%position(:,counter_equal)
            enddo ! counter_equal=1,n_equal_atoms
         enddo ! i_ua=1,N_unique_atoms
         close(io_gxepe)
         call return_iounit(io_gxepe)
      endif ! ex_gxepe
#endif
  end if

   ! loop over unique charges to read nml=pointcharge and
   ! positions of first equal charge
   do i_pc=1,pointcharge_N

      if ( .not. input_line_is_namelist("pointcharge") ) call input_error( &
        "pointcharge_read: namelist pointcharge exspected")

      pc => pointcharge_array(i_pc)

      N_equal_charges = df_N_equal_charges
      Z               = df_Z
      name            = df_name
      c               = df_c
      a               = df_a
      call input_read_to_intermediate()
      read(unit, nml=pointcharge, iostat=status)
      if (status .ne. 0) call input_error( &
           "pointcharge_read: namelist pointcharge.")
      if (N_equal_charges .le. 0) call input_error( &
           "pointcharge_read: namelist pointcharge: N_equal_charges")
      pc%name = name
      pc%Z = Z
      pc%N_equal_charges = N_equal_charges
      pc%c = c
      pc%a = a
      pc%cf = c
      pc%af = a
      ! read coordinates of first equal charge
      call input_read_to_intermediate()
      read(unit, fmt=*, iostat=status) (pc%position_first_ec(k), k=1,3)
      if (status .ne. 0) call input_error("pointcharge_read: position.")

   enddo

   end subroutine pointcharge_read
   !*************************************************************


   !*************************************************************
   subroutine pointcharge_write(unit)
   !  Purpose: writes names and Z of unique charges
   !           and coordinates of first centers ("equal charges")
   !           in format suitable for input to unit.
   ! called by write_input()
   !
   use echo_input_module, only: start, real, flag, intg, word, stop, &
        echo_level_full, word_format
   use operations_module, only: operations_echo_input_level
   implicit none
   integer(kind=i4_kind), intent(in) :: unit
   !** End of interface *****************************************

   !------------ Declaration of local variables -----------------
   type(pointcharge_type), pointer :: pc
   integer(kind=i4_kind)           :: i_pc, status, i_timp, lmax_pseudo
   character(len=32)               :: header
   type(unique_atom_type), pointer :: ut
   !------------ Executable code --------------------------------

   call start("POINTCHARGE_NUMBER","POINTCHARGE_WRITE",&
        unit,operations_echo_input_level)
   call intg("N_TIMPS",n_timps,df_n_timps)
   call intg("POINTCHARGE_N",pointcharge_N,df_pointcharge_N)
   call stop()

   word_format = '("    ",a," = ",a14  :" # ",a)' ! including the quotes

   ! loop over unique timps
   do i_timp=1,n_timps
      ut=> unique_timps(i_timp)
      name=ut%name
      z = ut%z
      zc = ut%zc
      n_equal_timps = ut%n_equal_atoms
      lmax_pseudo = ut%lmax_pseudo
      write(header,'("TIMP # timp ",i5)') i_timp
      call start(header,"POINTCHARGE_WRITE",unit,operations_echo_input_level)
      call word("NAME           ",name           ,df_name         )
      call real("Z              ",Z              ,df_Z            )
      call real("ZC             ",ZC             ,df_ZC           )
      call intg("N_EQUAL_TIMPS  ",N_equal_timps  ,df_N_equal_timps)
      call intg("LMAX_PSEUDO    ",lmax_pseudo    ,df_lmax_pseudo)
      call flag("FIXED        ", ut%moving_atom == 0, df_fixed        )
      call stop(empty_line=.false.)

      ! write coordinates of first equal charge
      write(unit, fmt='(3E25.15E2/)',iostat=status) ut%position_first_ea
      if (status .gt. 0) call error_handler("pointcharge_write: position.")

   end do
   ! loop over unique charges to read nml=pointcharge and
   ! positions of first equal charge
   do i_pc=1,pointcharge_N

      pc => pointcharge_array(i_pc)

      name = pc%name
      Z = pc%Z
      N_equal_charges = pc%N_equal_charges

      write(header,'("POINTCHARGE # point charge ",i5)')i_pc
      call start(header,"POINTCHARGE_WRITE",unit,operations_echo_input_level)
      call word("NAME           ",name           ,df_name           )
      call real("Z              ",Z              ,df_Z              )
      call intg("N_EQUAL_CHARGES",N_equal_charges,df_N_equal_charges)
      call stop(empty_line=.false.)

      ! write coordinates of first equal charge
      write(unit, fmt='(3E25.15E2/)',iostat=status) pc%position_first_ec
      if (status .gt. 0) call error_handler("pointcharge_write: position.")

   enddo

   end subroutine pointcharge_write
   !*************************************************************


   !*************************************************************
   subroutine init_pc_array(N_unique_pc,Uniq_pos,N_eq,Q,name,C,A,Cf,Af)
     ! Purpose: initialize pointcharge_array not from input
     ! Has to be called before MAIN_SYMM
     !** End of interface *****************************************
     !------------ Modules used -----------------------------------
     !------------ Declaration of formal parameters ---------------
     integer(i4_kind), intent(in) :: N_unique_pc
     real(r8_kind),    intent(in) :: Uniq_pos(3,N_unique_pc)
     integer(i4_kind), intent(in) :: N_eq(N_unique_pc)
     real(r8_kind),    intent(in) :: Q(N_unique_pc)
     character(len=12),     intent(in) :: name(N_unique_pc)
     real(r8_kind),    intent(in) :: C(N_unique_pc)
     real(r8_kind),    intent(in) :: A(N_unique_pc)
     real(r8_kind),    intent(in) :: Cf(N_unique_pc)
     real(r8_kind),    intent(in) :: Af(N_unique_pc)
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: status,i
     !------------ Executable code --------------------------------

     pointcharge_N=N_unique_pc
     allocate(pointcharge_array(pointcharge_N), stat=status)
     ASSERT(status==0)

     do i=1,pointcharge_N
        pointcharge_array(i)%Z=Q(i)
        pointcharge_array(i)%name=trim(name(i))
        pointcharge_array(i)%position_first_ec=Uniq_pos(:,i)
        pointcharge_array(i)%n_equal_charges=N_eq(i)
        pointcharge_array(i)%c=C(i)
        pointcharge_array(i)%a=A(i)
        pointcharge_array(i)%cf=Cf(i)
        pointcharge_array(i)%af=Af(i)
     end do

   end subroutine init_pc_array
   !*************************************************************

  !*****************************************************************************
  subroutine pointcharge_bcast()
    !  Purpose:  broadcasts all information contained in
    !            pointcharge_array(:)
    ! called by send_recv_init_options
    !------------ Modules used -------------------------------------------------
#ifdef WITH_EFP
    use operations_module, only: operations_solvation_effect
#endif
    use unique_atom_methods, only: unique_atom_pseudopot_bcast                 &
                                 , unique_atom_alloc
    use comm,                only: comm_bcast                                  &
                                 , comm_rank
    implicit none
    !------------ Declaration of local variables -----------------
    type(pointcharge_type), pointer :: pc
    type(unique_atom_type), pointer :: ut
    integer(i4_kind)                :: i_pc, i_timp, l_pseudo, status
    logical                         :: do_alloc
    !------------ Executable code --------------------------------
    !
    call comm_bcast( n_timps)
    call comm_bcast( N_moving_unique_timps)
    if( n_timps > 0 ) pseudopot_present=.true.
    call comm_bcast( pointcharge_N)
    call comm_bcast( moving_pc)
    !
    if( comm_rank() /= 0 .and. .not. associated(unique_timps)      ) then
      allocate( unique_timps(n_timps)           , stat=status )
      ASSERT(status==0)
    end if
    !
    do_alloc=.false.
    !
    if( comm_rank() /= 0 .and. .not. associated(pointcharge_array) ) then
      allocate( pointcharge_array(pointcharge_N), stat=status )
      ASSERT(status==0)
      do_alloc=.true.
    end if
    !
    ! loop over unique timps
    do i_timp = 1, n_timps
      ut => unique_timps(i_timp)
      !
      call comm_bcast( ut%z             )
      call comm_bcast( ut%zc            )
      call comm_bcast( ut%moving_atom   )
      call comm_bcast( ut%lmax_pseudo   )
      call comm_bcast( ut%N_equal_atoms )
      !
      if( comm_rank() /= 0 ) then
        allocate( ut%position(3,ut%N_equal_atoms), stat=status )
      endif
      !
      ! unpack atom coordinates of all equal atoms
      call comm_bcast( ut%position      )
      ut%position_first_ea = ut%position(:,1)
      !
      call comm_bcast( ut%name          )
      !
      if( comm_rank() == 0 ) then
        call unique_atom_alloc( ut, lmax_pseudo=ut%lmax_pseudo )
      endif
      do l_pseudo = 0, ut%lmax_pseudo
        call unique_atom_pseudopot_bcast( ut%l_pseudopot(l_pseudo) )
      end do
      !
    end do
    !
    ! loop over unique point charges
    if( pointcharge_N > 0 ) then
      do i_pc = 1, pointcharge_N
        !
        pc => pointcharge_array(i_pc)
        !
        call comm_bcast( pc%Z               )
        call comm_bcast( pc%C               )
        call comm_bcast( pc%A               )
        call comm_bcast( pc%N_equal_charges )
        call comm_bcast( pc%name            )
        !
        if( comm_rank() /= 0 .and. do_alloc ) then
          allocate( pc%position(3,pc%N_equal_charges), stat=status )
          if (status .ne. 0) call error_handler( &
              "pointcharge_unpack: allocate of position failed.")
        end if
        !
        ! broadcast atom coordinates of all equal atoms
        call comm_bcast( pc%position        )
        !
        pc%position_first_ec = pc%position(:,1)
        !
#ifdef WITH_EFP
        ! pack efp groups of all equal pc
        if( efp .and. operations_solvation_effect ) then
          if( comm_rank() /= 0 .and. do_alloc ) then
            allocate( pc%group(pc%N_equal_charges), stat=status )
            if (status .ne. 0) call error_handler( &
                "pointcharge_unpack: allocate of pc group failed.")
          end if
          call comm_bcast( pc%group         )
        end if
#endif
     enddo
     !
     if(moving_pc) call pc_points_info_bcast()
     !
   endif

   contains

     !*************************************************************
     subroutine pc_points_info_bcast()
       !------------ Declaration of local variables -----------------
       integer(i4_kind) :: i, status
       integer(i4_kind) :: n1, n2, n3
       logical          :: do_alloc
       !------------ Executable code --------------------------------
       !
       do_alloc=.false.
       !
       if( comm_rank() /= 0 .and. .not. associated(unique_pc_grad_info) ) then
         allocate( unique_pc_grad_info(pointcharge_N), stat=status )
         ASSERT(status==0)
         do_alloc=.true.
       end if
       !
       do i = 1, pointcharge_N
         if ( comm_rank() == 0 ) then
           n1=size(unique_pc_grad_info(i)%m,1)
           n2=size(unique_pc_grad_info(i)%m,2)
           n3=size(unique_pc_grad_info(i)%m,3)
         endif
         call comm_bcast( n1 )
         call comm_bcast( n2 )
         call comm_bcast( n3 )
         !
         if( comm_rank() /= 0 .and. do_alloc ) then
           allocate( unique_pc_grad_info(i)%m(n1,n2,n3), stat=status )
           ASSERT(status==0)
         end if
         !
         call comm_bcast( unique_pc_grad_info(i)%m )
         !
       end do
       !
       if( comm_rank() /= 0 .and. do_alloc ) then
         allocate( unique_pc_index(pointcharge_N+1), stat=status )
         ASSERT(status==0)
       end if
       !
       call comm_bcast( unique_pc_index       )
       call comm_bcast( totsym_grad_pc_length )
       !
     end subroutine pc_points_info_bcast
     !
   end subroutine pointcharge_bcast
   !****************************************************************************

   !*************************************************************
   subroutine pointcharge_close()
     !  Purpose:  packs all information contained in
     !            pointcharge_array(:)
     ! called by send_recv_init_options
     !** End of interface *****************************************
     !------------ Modules used ------------------- ---------------
     use uatom_symmadapt, only: sa_free
     implicit none
     !------------ Declaration of local variables -----------------
     type(pointcharge_type), pointer :: pc
     type(unique_atom_type), pointer :: ut
     integer(kind=i4_kind)           :: i_pc, i_timp, l_pseudo, status
     integer(kind=i4_kind)           :: i,j,k
     !------------ Executable code --------------------------------
     do i_timp=1,n_timps
        ut=> unique_timps(i_timp)
        deallocate(ut%position,stat=status)
        if(status/=0) call error_handler("pointcharge_close: deallocating unique_timps &
             & position failed")
        do l_pseudo=0,ut%lmax_pseudo
           deallocate(ut%l_pseudopot(l_pseudo)%exponents,ut%l_pseudopot(l_pseudo)%coefficients,&
                ut%l_pseudopot(l_pseudo)%powers, stat=status )
           if(status/=0) call error_handler("pointcharge_close: deallocating exponents &
                & failed")
        end do
     end do

     do i_pc=1,pointcharge_N
        pc  => pointcharge_array(i_pc)
        deallocate(pc%position,stat=status)
        if(status/=0) call error_handler("pointcharge_close: deallocating pc position &
             & failed")
        if(comm_i_am_master()) then
           if(associated(pc%group)) then
              deallocate(pc%group,stat=status)
              ASSERT(status==0)
           endif
        end if
     end do
     if(associated(unique_timps)) then
        deallocate(unique_timps,stat=status)
        ASSERT(status==0)
     end if
     if (associated (pointcharge_array)) then
        deallocate (pointcharge_array, stat=status)
        ASSERT(status==0)
     end if
     if(associated(unique_pc)) then
        do i=1,pointcharge_N
           ut=> unique_pc(i)
           deallocate(ut%position,stat=status)
           ASSERT(status==0)
           do j=1,size(ut%symadapt_partner,1)
              do k=1,ut%lmax_all
                 call sa_free(ut%symadapt_partner(j,k))
              end do
           end do
           deallocate(ut%symadapt_partner,STAT=status)
           ASSERT(status==0)
        end do
        deallocate(unique_pc,stat=status)
        ASSERT(status==0)
     end if

     call unique_timp_gradinfo_dealloc()
   end subroutine pointcharge_close

   !*************************************************************
   subroutine unique_timp_gradinfo_dealloc()
     ! Purpose: deallocate the variable 'unique_timp_grad_info'
     !          and the pointer 'moving_unique_timp_index'
     ! subroutine called by: main_gradient
     !** End of interface *****************************************
     ! ----------- modules ------------------------------------
     ! ----------- declaration of local variables -------------
     integer(kind=i4_kind) :: i, alloc_stat
     external error_handler
      ! ----------- executable code -----------------------------

     if(N_moving_unique_timps > 0) then
        do i=1,N_moving_unique_timps
           deallocate(unique_timp_grad_info(i)%m,STAT=alloc_stat)
           if (alloc_stat /= 0 ) call error_handler &
                ("unique_timp_gradinfo_dealloc: deallocation (1) failed")
        enddo
        deallocate(unique_timp_grad_info,STAT=alloc_stat)
        if (alloc_stat /= 0 ) call error_handler &
             ("unique_timp_gradinfo_dealloc: deallocation (2) failed")
        deallocate(moving_unique_timp_index,STAT=alloc_stat)
        if (alloc_stat /= 0 ) call error_handler &
             ("unique_timp_gradinfo_dealloc: deallocation (3) failed")
     end if

     if(moving_pc) then
        do i=1,pointcharge_N
           deallocate(unique_pc_grad_info(i)%m,STAT=alloc_stat)
           ASSERT(alloc_stat==0)
        enddo
        deallocate(unique_pc_grad_info,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        deallocate(unique_pc_index,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
     end if
   end subroutine unique_timp_gradinfo_dealloc
   !*************************************************************

   subroutine unique_timp_grad_information()
     use symmetry_data_module, only : get_totalsymmetric_irrep
     ! Purpose: change the variable
     !  unique_timp(i)%symapdapt_partner(irrep,l)%symadapt(n_indep,n_partner)%
     !  N_fcts
     !  I_equal_atom
     !  c ...
     ! to a variable better suited for the gradient part where
     ! we loop over equal atoms. (-> better performance)
     ! and only the non-fixed unique_timps are considered.
     !
     ! Also loads the moving_unique_timp_index.
     !
     ! Subroutine called by: subroutine post_scf (if operations_gradients)
     !** End of interface *****************************************
 !*************************************************************
     ! ----------- declaration of local variables -------------
     type(unique_atom_type)         , pointer :: ua
     type(unique_atom_partner_type) , pointer :: sap
     type(unique_atom_symadapt_type), pointer :: sa
     type(arrmat3)                  , pointer :: gi
     integer(kind=i4_kind) :: i_ua,i_ma,i_if,i_cd,i_ea,ts, &
                              n_indep,n_equals,alloc_stat

     external error_handler
     ! ----------- executable code -----------------------------
     !

     if(N_moving_unique_timps > 0) then
     allocate( moving_unique_timp_index(N_moving_unique_timps),&
          STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("unique_timp_grad_information: allocation moving_unique_timp_index failed")
     moving_unique_timp_index=0_i4_kind

     allocate( unique_timp_grad_info(N_moving_unique_timps),STAT=alloc_stat)
     if (alloc_stat.ne.0) call error_handler &
          ("unique_timp_grad_information: allocation (1) failed")
     ts = get_totalsymmetric_irrep()

     do i_ua=1,n_timps
        ua => unique_timps(i_ua)
        sap => ua%symadapt_partner(ts,1) ! l = 1
        i_ma = ua%moving_atom
        if (i_ma > 0) then
           moving_unique_timp_index(i_ma) = i_ua
           gi => unique_timp_grad_info(i_ma)

           n_indep = sap%N_independent_fcts
           n_equals = ua%N_equal_atoms
           allocate ( gi%m(n_indep,3,n_equals),STAT=alloc_stat)
           if (alloc_stat.ne.0) call error_handler &
                ("unique_timp_grad_information: allocation (2) failed")
           gi%m = 0.0_r8_kind

           do i_if=1,n_indep ! loop over all sym. adapted gradients of ua
              sa => sap%symadapt(i_if,1) ! first partner
              do i_cd=1,sa%N_fcts ! loop over all contributing derivatives
                 i_ea = sa%I_equal_atom(i_cd)

                 select case (sa%m(i_cd))
                 case (2); gi%m(i_if,1,i_ea) = sa%c(i_cd)
                 case (3); gi%m(i_if,2,i_ea) = sa%c(i_cd)
                 case (1); gi%m(i_if,3,i_ea) = sa%c(i_cd)
                 case default
                    call error_handler("unique_timp_grad_information: sth. wrong")
                 end select

              enddo
           enddo
        end if
     end do
     endif

   end subroutine unique_timp_grad_information

   !****************************************************************************
   subroutine unique_timp_symadapt_bcast()
   !  Purpose: broadcasts in information about symmetry adaption for all
   !           unique atoms in unique_atoms(:) and does allocation on slaves
   !           Informations in irrep_module are used and must be existant
   !------------ modules used --------------------------------------------------
   use comm,                 only: comm_bcast                                  &
                                 , comm_rank
   use symmetry_data_module, only: get_totalsymmetric_irrep
   use uatom_symmadapt,      only: sa_cbcast
   implicit none
   !------------ Declaration of local variables --------------------------------
   type(unique_atom_type),          pointer :: ua
   type(unique_atom_partner_type),  pointer :: sap
   integer                                  :: i,j,l,status
   integer(kind=i4_kind) :: ts
   !------------ Declaration of subroutines used -------------------------------
   external error_handler
   !------------ Executable code -----------------------------------------------
   !
   ! loop over unique atoms
   ts = get_totalsymmetric_irrep()
   do i = 1, N_timps
     ua => unique_timps(i)
     !
     if ( comm_rank() /= 0 ) then
       allocate( ua%symadapt_partner(ts:ts,1:1), stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "unique_timp_symadapt_unpack: ua%symadapt_partner allocation failed")
       ! is left unallocated: ua%symadapt_partner(ts,1)%symadapt
     endif
     !
     ! loop over angular momentum
     do l = 1, 1 !??????
       ! loop over irreps
       do j = ts, ts !??????
         !
         sap => ua%symadapt_partner(j,l)
         call sa_cbcast( sap )
         !
       enddo
     enddo
   enddo
   !
   end subroutine unique_timp_symadapt_bcast
   !****************************************************************************

   !*************************************************************
   subroutine init_pointcharges_grads()
     !  Purpose: allocate arrays for calculating gradients on PC
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: status,i,n_eq
     !------------ Executable code --------------------------------

     allocate(gradient_pc_cartesian(pointcharge_N),stat=status)
     ASSERT(status==0)
#ifdef WITH_EFP
     if(efp) then
        allocate(torque_pc_cartesian(pointcharge_N),stat=status)
        ASSERT(status==0)
     end if
#endif

     do i=1,pointcharge_N
        n_eq=pointcharge_array(i)%N_equal_charges
        allocate(gradient_pc_cartesian(i)%m(3,n_eq),stat=status)
        ASSERT(status==0)
        gradient_pc_cartesian(i)%m=0.0_r8_kind
#ifdef WITH_EFP
        if(efp) then
           allocate(torque_pc_cartesian(i)%m(3,n_eq),stat=status)
           ASSERT(status==0)
           torque_pc_cartesian(i)%m=0.0_r8_kind
        end if
#endif
     end do

     allocate(gradient_pc_totalsym(totsym_grad_pc_length),stat=status)
     ASSERT(status==0)
     gradient_pc_totalsym=0.0_r8_kind
#ifdef WITH_EFP
     if(efp) then
        allocate(torque_pc_totalsym(totsym_grad_pc_length),stat=status)
        ASSERT(status==0)
        torque_pc_totalsym=0.0_r8_kind
     end if
#endif

   end subroutine init_pointcharges_grads

   !*************************************************************

   subroutine pc_grads_shutdown()
     !  Purpose: free arrays for calculating gradients on PC
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: status,i
     !------------ Executable code --------------------------------

     do i=1,pointcharge_N
        deallocate(gradient_pc_cartesian(i)%m,stat=status)
        ASSERT(status==0)
#ifdef WITH_EFP
        if(efp) then
           deallocate(torque_pc_cartesian(i)%m,stat=status)
           ASSERT(status==0)
        end if
#endif
     end do
     deallocate(gradient_pc_cartesian,stat=status)
     ASSERT(status==0)
#ifdef WITH_EFP
     if(efp) then
        deallocate(torque_pc_cartesian,stat=status)
        ASSERT(status==0)
     end if
#endif

     deallocate(gradient_pc_totalsym,stat=status)
     ASSERT(status==0)
#ifdef WITH_EFP
     if(efp) then
        deallocate(torque_pc_totalsym,stat=status)
        ASSERT(status==0)
     end if
#endif

   end subroutine pc_grads_shutdown

   !*************************************************************

   subroutine totsym_PC_grad_unpack()
     !  Purpose: unpack 2 center contribution to the PC gradient from the
     !           slave
     !------------ Declaration of local variables ----------------
     real(r8_kind) :: help_arr(totsym_grad_pc_length)
     integer(i4_kind) :: info
     !------------ Executable code -------------------------------

     call communpack(help_arr,totsym_grad_pc_length,1,info)
     ASSERT(info==0)
     gradient_pc_totalsym=gradient_pc_totalsym+help_arr

   end subroutine totsym_PC_grad_unpack

   !*************************************************************

   subroutine totsym_PC_grad_pack()
     !  Purpose: pack 2 center contribution to the PC gradient on the
     !           slave
     !------------ Declaration of local variables ----------------
     integer(i4_kind) :: info
     !------------ Executable code -------------------------------

     call commpack(gradient_pc_totalsym,totsym_grad_pc_length,1,info)
     ASSERT(info==0)

   end subroutine totsym_PC_grad_pack

   !*************************************************************

   subroutine transform_PC_grad_to_cart()
     ! purpose: transform symmetry adapted gradient components to cartesian
     !          coordinates and add them to an array of cartesian gradients
     !
     !------------ Declaration of local variables ----------------
     integer(kind=i4_kind) :: i_unique,i_equal,index,i,grad_dim
     real(kind=r8_kind),pointer :: rotmat(:,:)
#ifdef WITH_EFP
     integer(kind=i4_kind) :: i_group
     real(kind=r8_kind) :: r12(3),r1(3),r2(3),f(3)
#endif
     !------------ Executable code -------------------------------

     do i_unique=1,pointcharge_N
        index=unique_pc_index(i_unique)
        grad_dim=unique_pc_index(i_unique+1)-index
        index=index-1
        do i_equal=1,pointcharge_array(i_unique)%n_equal_charges
           rotmat=>unique_pc_grad_info(i_unique)%m(:,:,i_equal)
           do i=1,grad_dim
              gradient_pc_cartesian(i_unique)%m(:,i_equal) = &
                   gradient_pc_cartesian(i_unique)%m(:,i_equal) + &
                   rotmat(i,:)*gradient_pc_totalsym(index+i)
#ifdef WITH_EFP
             if(efp) then
                torque_pc_cartesian(i_unique)%m(:,i_equal) = &
                     torque_pc_cartesian(i_unique)%m(:,i_equal) + &
                     rotmat(i,:)*torque_pc_totalsym(index+i)
             end if
#endif
           end do
        end do
     enddo

#ifdef WITH_EFP
     if(efp) then
        do i_unique=1,pointcharge_N
           do i_equal=1,pointcharge_array(i_unique)%n_equal_charges
              r1=pointcharge_array(i_unique)%position(:,i_equal)
              i_group=pointcharge_array(i_unique)%group(i_equal)
              r2=rcm(:,i_group)
              r12=r2-r1
              f=-gradient_pc_cartesian(i_unique)%m(:,i_equal)

              torque_pc_cartesian(i_unique)%m(:,i_equal) = &
                   torque_pc_cartesian(i_unique)%m(:,i_equal)+vector_product(r12,f)
           end do
        end do
     end if
#endif

   end subroutine transform_PC_grad_to_cart

   !*************************************************************

   function vector_product(v1,v2)
     !------------ Modules used ------------------- ---------------
     !------------ Declaration of formal parameters ---------------
     real(r8_kind) :: vector_product(3)
     real(r8_kind) :: v1(3),v2(3)
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     !------------ Executable code --------------------------------

     vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
     vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
     vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

   end function vector_product

   !*************************************************************

   subroutine pc_grad_cart_write()
     !  Purpose: writing the cartesian gradients
     !** End of interface *****************************************
     !------------ modules used -----------------------------------
     use iounitadmin_module, only: output_unit
     !------------ Declaration of local variables -----------------
     integer(i4_kind) :: i,j
     !------------ Executable code --------------------------------

#ifdef WITH_EFP
     if(efp) then
        write(output_unit,'(/A)') 'Cartesian gradients and torques on Point Charges'
        do i=1,pointcharge_N
           write(output_unit,*) 'Unique Point Charge:',i
           do j=1,pointcharge_array(i)%n_equal_charges
              write(output_unit,'(A14,3F15.10,a3,3F15.10)') 'Equal Center: ',&
                   gradient_pc_cartesian(i)%m(:,j),' / ',torque_pc_cartesian(i)%m(:,j)
           enddo
        end do
     else
#endif
        write(output_unit,'(/A)') 'Cartesian gradients on Point Charges'
        do i=1,pointcharge_N
           write(output_unit,*) 'Unique Point Charge:',i
           do j=1,pointcharge_array(i)%n_equal_charges
              write(output_unit,'(A14,3F15.10)') 'Equal Center: ',&
                   gradient_pc_cartesian(i)%m(:,j)
           enddo
        end do
#ifdef WITH_EFP
     end if
#endif

   end subroutine pc_grad_cart_write

   !*************************************************************
  function calc_nuc_pc_energy() result (e_nuc_pc)
    ! Purpose: Calculate interaction energy between atomic nuclears
    ! end external centers (dipoles, quadrupoles, octopoles)
    use unique_atom_module, only : unique_atoms, N_unique_atoms
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: e_nuc_pc
    !** End of interface *****************************************
    real(kind=r8_kind)          :: z1,z2,dist,dist2
    integer(kind=i4_kind)       :: na, nb, eq_a, eq_b, n_equal_charges
    integer(kind=i4_kind)       :: N_equal_atoms_a
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:)
    real(kind=r8_kind)          :: C,A
    !----------- declaration of local variables -------------

    e_nuc_pc=0.0_r8_kind
    unique_pc: do nb=1,pointcharge_N
       z2 = pointcharge_array(nb-n_timps)%Z
       xb => pointcharge_array(nb-n_timps)%position
       n_equal_charges = pointcharge_array(nb)%N_equal_charges
       C=pointcharge_array(nb-n_timps)%C
       A=pointcharge_array(nb-n_timps)%A

       unique_a: do na=1,N_unique_atoms
          z1 = unique_atoms(na)%Z-unique_atoms(na)%ZC
          xa => unique_atoms(na)%position
          N_equal_atoms_a=unique_atoms(na)%N_equal_atoms

          equal_a: do eq_a=1,N_equal_atoms_a
             equal_pc: do eq_b=1,n_equal_charges
                dist2 = sum((xa(:,eq_a)-xb(:,eq_b))**2)
                dist = sqrt(dist2)
                if(dist <= 1.0e-10_r8_kind) cycle equal_pc
                e_nuc_pc = e_nuc_pc + z1*z2/dist
#ifdef WITH_EFP
                if(C /= 0.0_r8_kind) then
                   e_nuc_pc = e_nuc_pc - C*exp(-A*dist2)*z1*z2/dist
                end if
#endif
             enddo equal_pc
          enddo equal_a

       enddo unique_a
    enddo unique_pc
#ifdef WITH_EFP
    if((pointcharge_N > 0) .and. efp .and. print_energy) &
         print*,'NUC-PC  EFP ',e_nuc_pc
#endif

  end function calc_nuc_pc_energy

   !*************************************************************
   subroutine calc_PC_grads()
     ! Purpose: Calculate gradients on PC (interaction with atomic nuclears)
     use unique_atom_module, only : unique_atoms, N_unique_atoms, &
          pseudopot_present
     !------------ Declaration of formal parameters ---------------
     !** End of interface *****************************************
     !----------- declaration of local variables -------------
     integer(i4_kind) :: i_up,i_ep,n_ep,i_ua,i_ea,n_ea,i
     real(r8_kind) :: Z_pc,Z_a,C,A
     real(r8_kind), pointer :: xp(:,:),xa(:,:),rotmat(:,:)
     integer(i4_kind) :: grad_dim,index
     real(r8_kind) :: gradient(3),r_ap(3),d_ap,d_ap2
     !--- executable code-------------------------------------

     do i_up=1,pointcharge_N
        n_ep=pointcharge_array(i_up)%N_equal_charges
        Z_pc=pointcharge_array(i_up)%Z
        xp=>pointcharge_array(i_up)%position
        C=pointcharge_array(i_up)%C
        A=pointcharge_array(i_up)%A

        grad_dim=unique_pc_index(i_up+1)-unique_pc_index(i_up)
        index=unique_pc_index(i_up)-1

        do i_ep=1,n_ep
           rotmat=>unique_pc_grad_info(i_up)%m(:,:,i_ep)
           gradient=0.0_r8_kind

           do i_ua=1,N_unique_atoms
              n_ea=unique_atoms(i_ua)%N_equal_atoms
              Z_a=unique_atoms(i_ua)%Z
              if(pseudopot_present) Z_a=Z_a-unique_atoms(i_ua)%Zc
              xa=>unique_atoms(i_ua)%position

              do i_ea=1,n_ea
                 r_ap=xa(:,i_ea)-xp(:,i_ep)
                 d_ap2=dot_product(r_ap,r_ap)
                 d_ap=sqrt(d_ap2)
                 if(d_ap < 0.001_r8_kind) d_ap=0.001_r8_kind
                 gradient=gradient-Z_a*Z_pc*r_ap/d_ap**3
                 if(C /= 0.0_r8_kind) then
                    gradient=gradient+C*exp(-A*d_ap2)*Z_a*Z_pc*(1.0+2.0_r8_kind*A*d_ap2)*r_ap/d_ap**3
                 end if
              end do
           end do

           do i=1,grad_dim
              gradient_pc_totalsym(index+i)=gradient_pc_totalsym(index+i)- &
                   sum(rotmat(i,:)*gradient(:))
           enddo

        end do
     end do

   end subroutine calc_PC_grads

end module pointcharge_module

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
module grid_module
  !
  !  This modules generates the 3D numerical integration grid which is
  !  used  for integration  of the  xc-contributions.  The  grid  is a
  !  superposition  of  atom centered  grids.   The  atomic grids  are
  !  stored in  the variable  agrid.  Every atomic  grid consits  of a
  !  radial  and  a angulargrid.   The  radial  part  is calculted  in
  !  routine radgrid, the angular part in inigrid.
  !
  !  The orientation  of the  grid is given  by a rotation  matrix ro.
  !  This is  calculated in  routine get_symmetry. If  input parameter
  !  sym_reduce is  set to false, the  grid is exactly the  same as in
  !  the  old lcgto.  The z-axis  is  pointing towards  the center  of
  !  nuclear charge. If sym_reduce is  true, the grid is orientated in
  !  such a way, that the  symmetry can be exploided best. This option
  !  should be used with care, because not all possible cases could be
  !  tested.
  !
  !  The integration weight is a product of geometric weight and becke
  !  weight.   The  geometric  weight  is determined  in  radgrid  and
  !  inigrid,  the  becke  weight,  which  is  a  consequence  of  the
  !  superposition  of the  atomic  grids, is  calculated in  function
  !  atomicweight.
  !
  !  The atomic grids are  stored in variable agrid. After calculation
  !  of  the  becke weights  the  gridpoints  are  copied to  variable
  !  grdpts.
  !
  !  All parameters  describing the grid  are stored in  the structure
  !  grid_par.
  !
  !
  !
  !
  !  Module called by: main_scf,post_scf_main,main_slave
  !
  !
  !  References: Diplom thesis P. Ulbricht and references therein
  !
  !
  !  Author: MS
  !  Date: 11/95
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification: - Parallel algorithm for generation of the grid
  !               - post_scf grid implemented
  ! Author: MS
  ! Date:   9/96
  ! Description: radialgrid and angulargrid are still calculated on
  !              the master. Before calculation of the becke weights
  !              the grid is spread to the slaves. The calculation of
  !              the becke weights is done on the slaves.
  !
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   3/97
  ! Description: Subroutine for calculation of derivatives of the grid
  !              added Also steering parameter weight_grads is true.
  !              If this parameter is true, gradients of weight will
  !              be calculated and no Becke scaling and rotation of
  !              the grid will be performed
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        17/7/1998
  ! Description: separate variables sym_reduce_scf and sym_reduce_ph
  !              introduced check routine for equal SCF and POST SCF
  !              grids introduced
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !---------------------------------------------------------------

# include "def.h"
  use type_module
  use unique_atom_module
  use datatype
  use iounitadmin_module
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
#endif
  USE_MEMLOG
  implicit none
  private
  save
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  public :: grid_main
  public :: grid_read
  public :: grid_write
  public :: grid_read_ph
  public :: grid_write_ph
  public :: atomicweight
  public :: grid_close
  public :: grid_bcast
  public :: atomicweight_and_grad
  public :: grid_copy_to_ph
  public :: equal_scf_and_ph_grids
  public :: atomicweight_and_dervs

  public :: grid_main1
  public :: grid_close1 ! for call from master-only context

  !
  ! For SCF:
  !
  public :: grid_loop_setup!()
  public :: more_grid!(n, pts, wts) -> logical

  !
  ! For Post-SCF:
  !
  public :: grid_loop_setup_atom!()
  public :: more_grid_atom!(n, ua, pts, wts) -> logical

#ifdef WITH_GUILE
  public :: guile_pople_radius  ! (SCM int) -> SCM double
  public :: guile_slater_radius ! (SCM int) -> SCM double
  public :: guile_ionic_radius  ! (SCM int) -> SCM double
#endif

  public print_alloc_grid

  logical, public :: weight_grads ! decide if for gradient calculation also
  ! the gradients of integration weights have to be calculated

!================================================================
! End of public interface of module
!================================================================

  logical, private :: adjust_cell_size = .false.
  ! Base grid weights on Voronoj cells that are built by
  ! planes cutting the bonds not just in the middle.
  ! Instead, account for different atomic sizes.
  ! In older builds
  !            adjust_cell_size == .not. weight_grads
  ! was hardwired, as the cell scaling was not implemented
  ! for gradients of weights.

  integer(kind=i4_kind):: alloc_grid(4)
  ! 3 - agrid
  type GRID_ATOM  ! structure contains datas for grid_construction
     integer            :: NRAD    ! measure for the number of radial shells
     integer            :: NANG    ! angular momentum of the Lebedev grid
     real(kind=R8_KIND) :: RADIUS  ! OF THE ATOM
     integer            :: NPART   ! NUMBER OF DOMAINS FOR RADIAL INTEGRATION
     integer,dimension(5) :: NANGA ! L VALUES FOR DIFFERENT DOMAINS
     real(kind=R8_KIND),dimension(4) :: ALPHA
     integer                           :: NRADPTS
  end type GRID_ATOM

  type(grid_atom),dimension(:),allocatable,private,target :: &
       grid_par_scf,grid_par_ph

  type(grid_atom),dimension(:),pointer :: grid_par

  type(arrmat2), allocatable, target, private :: agrid(:) ! (n_ua), complete atomic grid
  real(kind=r8_kind), allocatable, target, private :: grdpts(:,:)

  !
  ! FIXME: target is here only to make things like
  !
  !     p_points  => agrid(unique_atom)%m(jobs(1)+1:jobs(2), 1:3)
  !
  ! compile also when arrmat2 is declared with
  ! allocatable component instead of pointer.
  !

  real(kind=R8_KIND),dimension(3)  :: CHARGE_CENTER

  logical,private            :: sym_reduce_scf, sym_reduce_ph
  logical,private            :: weight_grads_scf, weight_grads_ph
  ! controll switch which deceides
  ! if grid reduction due to the local symmetry is used
  logical,private            :: gridatom_found, gridatom_ph_found
  ! control switch for grid_write and grid_write_ph

  ! Upper Limits of a few input parameter
  integer(kind=i4_kind), parameter, private :: max_nrad     = 200, &
                                               max_nang     = 291, &
                                               max_npart    = 5  , &
                                               max_nrad_ph  = 200, &
                                               max_nang_ph  = 291, &
                                               max_npart_ph = 5
  real(kind=r8_kind), parameter, private :: max_radius    = 10.0_r8_kind, &
                                            max_radius_ph = 10.0_r8_kind

  ! The Default input parameters
  integer(kind=i4_kind), private :: df_nrad     =  30, &
                                    df_nang     = 131, &
                                    df_npart    =   1, &
                                    df_nrad_ph  =  50, &
                                    df_nang_ph  = 291, &
                                    df_npart_ph =   1
  real(kind=r8_kind), private :: df_radius    = 0.0_r8_kind, &
                                 df_radius_ph = 0.0_r8_kind
  logical, private :: df_rslater       = .true. , &
                      df_rpople        = .false., &
                      df_rionic        = .false., &
                      df_sym_reduce    = .false., &
                      df_rslater_ph    = .TRUE. , &
                      df_rpople_ph     = .false., &
                      df_rionic_ph     = .false., &
                      df_weight_grads  = .true., & ! vn: many parts
                                                   ! just require this
                                                   ! for better more
                                                   ! exact work
                      df_sym_reduce_ph = .false., & ! = df_sym_reduce!
                      df_weight_grads_ph = .true. ! = df_weight_grads!

  ! The local input paramters
  integer(kind=i4_kind), private :: nrad , &
                                    nang , &
                                    npart
  real(kind=r8_kind), private :: radius
  logical, private :: rslater, &
                      rpople , &
                      rionic , &
                      sym_reduce

  ! The input namelists
  namelist /grid/ sym_reduce, weight_grads
  namelist /gridatom/ nrad   , &
                      nang   , &
                      npart  , &
                      radius , &
                      rslater, &
                      rpople , &
                      rionic
  namelist /grid_ph/ sym_reduce  , &
                     weight_grads
  namelist /gridatom_ph/ nrad   , &
                         nang   , &
                         npart  , &
                         radius , &
                         rslater, &
                         rpople , &
                         rionic

  real(kind=r8_kind),parameter ::&
       & one = 1.0_r8_kind,&
       & two = 2.0_r8_kind

  logical, private :: output_grid = .false. ! will be reset in grid_main()

  FPP_TIMER_DECL(gmain)
  FPP_TIMER_DECL(mkgrd)
  FPP_TIMER_DECL(awt)
  FPP_TIMER_DECL(awt0)
  FPP_TIMER_DECL(awt1)
  FPP_TIMER_DECL(awt2)

contains

  subroutine grid_main1 (post_scf)
    use comm_module, only: comm_init_send, comm_all_other_hosts, comm_send, &
         comm_i_am_master
    use msgtag_module, only: msgtag_gr_send
    implicit none
    logical, intent( in  ) :: post_scf
    ! *** end of interface ***

    if (comm_i_am_master()) then
      call comm_init_send (comm_all_other_hosts, msgtag_gr_send)
      call comm_send()
    end if

    call grid_main (post_scf)
  end subroutine grid_main1

  subroutine grid_main(post_scf)
    use output_module, only: output_module_output_grid => output_grid
    implicit none
    logical, optional, intent(in) :: post_scf
    !** End of interface *****************************************

    type(grid_atom),dimension(:),pointer :: grid_par1
    logical :: sym_reduce1, weight_grads1

    FPP_TIMER_START(gmain)

    !
    ! Only  master  is  allowed  to  use output_unit,  on  slaves  the
    ! output_unit  may not  even be  opened.  By  convention  the unit
    ! number will  be negative in this  case. This is  a global module
    ! varibale,  named for  historical  reasons the  same  way as  the
    ! public variable in the output_module:
    !
    output_grid = output_module_output_grid .and. output_unit > 0

    if(present(post_scf)) then
       if(post_scf) then
          if (output_grid) write(output_unit,*)  'Entering the grid section for post scf'
          grid_par1 => grid_par_ph
          sym_reduce1 = sym_reduce_ph
          weight_grads1 = weight_grads_ph
       else
          if (output_grid) write(output_unit,*) &
               'Entering the grid section'
          grid_par1 => grid_par_scf
          sym_reduce1 = sym_reduce_scf
          ! for compatibility with older versions
          weight_grads1 = weight_grads_ph ! not weight_grads_scf
       end if
    else
       if (output_grid) write(output_unit,*) &
            'Entering the grid section'
       grid_par1 => grid_par_scf
       sym_reduce1 = sym_reduce_scf
       ! for compatibility with older versions
       weight_grads1 = weight_grads_ph ! not weight_grads_scf
    end if

    ! constructs the grid --- points and weights:
    call grid_main_calc(grid_par1, sym_reduce1, weight_grads1)

    !
    ! Next   mkgrid()  will  scale   quadrature  weights   by  (Becke)
    ! partitioning  function and  collect  all grid  points  are in  a
    ! continuous array.  The jobs will be  distributed dynamically, so
    ! there is no chance to know beforehand which one will stay
    !
    if ( .not. post_scf ) then
        !
        ! Note: in  Post SCF the  atomic grid structure AGRID  is used
        ! directly.
        !
        call mkgrid(agrid, grdpts)
    endif

    FPP_TIMER_STOP(gmain)
  end subroutine grid_main

  subroutine grid_main_calc(grid_par1, sym_reduce1, weight_grads1 )
    ! main subroutine for creation of integration gird executed only
    ! on the master
    use TIME_MODULE
    use symmetry_element, only: set_sym_element
    use symmetry_data_module, only: symmetry_data_point_group
    implicit none
    type(grid_atom),dimension(:),target, intent(in) :: grid_par1
    logical, intent(in) :: sym_reduce1, weight_grads1
    !** End of interface *****************************************

    real(kind=R8_KIND),dimension(:,:),pointer :: RADGRD,ANGGRD,TOTGRD, helpgrd
    integer :: I, J, ALLOCSTAT, L, GR_SIZE, START, angsize
    integer,dimension(5)   :: GR_SIZE_INC
    integer,dimension(2,11) :: nang_par=reshape((/31,6,51,14,71,26,91,38,111,&
         50,131,78,151,86,171,110,191,146,231,194,291,302/),(/2,11/))
    grid_par => grid_par1
    sym_reduce = sym_reduce1
    weight_grads = weight_grads1

    if (output_grid) write(output_unit,*) &
         '================================================'
    if (output_grid) then
       write(output_unit,*) '================================================'
       write(output_unit,'("sym_reduce    = ",L1)') sym_reduce
       write(output_unit,'("weight_grads  = ",L1)') weight_grads
    endif

    ! Use atomic radii to adjust the Voronoj cells:
    adjust_cell_size = .not. weight_grads
    !WARN('adjust_cell_size==true')

    call calc_charge()
    if (output_grid) then
       write(output_unit,'("charge_center = ",3F16.10)') charge_center
    endif
    call set_grid_par()

    ! initializes symmetry_element module:
    call set_sym_element(symmetry_data_point_group())

    if (.not. allocated(AGRID)) then
      allocate(AGRID(N_UNIQUE_ATOMS),STAT=ALLOCSTAT)
      ASSERT (ALLOCSTAT==0)
    endif

    do J=1,N_UNIQUE_ATOMS
       if (output_grid) write(output_unit,*) '============================================'
       if (output_grid) write(output_unit,*) 'Creating the gridpoints for atom no.',j
       if (output_grid) write(output_unit,*) 'Input Parameters are:'
       if (output_grid) write(output_unit,'(3(a7,i4,2x),a7,f7.4)') 'NRAD:',grid_par(j)%nrad,&
            'NANG:',grid_par(j)%nang,'NPART:',grid_par(j)%npart,'Radius:',&
            grid_par(j)%radius

       ! Calculation of the radialgrid
       RADGRD=>RADIALGRID(0,GRID_PAR(J)%NRADPTS,GRID_PAR(J)%RADIUS)
       !
       if (output_grid) write(output_unit,*) 'The grid contains ',size(radgrd,1),'shells'

       if(GRID_PAR(J)%NPART==1) then    ! NO GRID PRUNING

          ! calculate angular grid
          ANGGRD => angular_grid(J, GRID_PAR(J)%NANG)

          allocate(AGRID(J)%M(size(RADGRD,1)*size(ANGGRD,1),4),STAT=ALLOCSTAT)
          ASSERT (ALLOCSTAT==0)
          TOTGRD=>AGRID(J)%M        ! JUST FOR SIMPLIFICATION

          ! build product of angular and radial grid
          do I=1,size(RADGRD,1)
             ! positions totgrd(:,1:3)
             TOTGRD((I-1)*size(ANGGRD,1)+1:I*size(ANGGRD,1),1:3)=ANGGRD(:,1:3)* &
                  RADGRD(I,1)
             ! one-cente integration weights weights totgrd(:,4)
             TOTGRD((I-1)*size(ANGGRD,1)+1:I*size(ANGGRD,1),4)=ANGGRD(:,4)* &
                  RADGRD(I,2)*16.0_R8_KIND*atan(1.0_R8_KIND)
          enddo
       else                     !GRID PRUNING
          do L=1,5
             do I=1,11
                if(GRID_PAR(J)%NANGA(L).LE.NANG_PAR(1,I)) then
                   GRID_PAR(J)%NANGA(L)=NANG_PAR(1,I)
                   GR_SIZE_INC(L)=NANG_PAR(2,I)
                   exit
                endif
             enddo

          enddo
          GR_SIZE=0
          L=1
          do I=1,size(RADGRD,1)
             if(L==5) then
                GR_SIZE=GR_SIZE+GR_SIZE_INC(L)
             else
                if(RADGRD(I,1)<=GRID_PAR(J)%ALPHA(L)) then
                   GR_SIZE=GR_SIZE+GR_SIZE_INC(L)
                else
                   L=L+1
                   GR_SIZE=GR_SIZE+GR_SIZE_INC(L)
                endif
             endif
          enddo
          allocate(helpgrd(GR_SIZE,4),STAT=ALLOCSTAT)
          ASSERT (ALLOCSTAT==0)
          START=1
          L=1
          ANGGRD => angular_grid(J, GRID_PAR(J)%NANGA(L))
          do I=1,size(RADGRD,1)
             if(L/=5) then
                if(RADGRD(I,1)>GRID_PAR(J)%ALPHA(L)*GRID_PAR(J)%RADIUS) then
                   L=L+1
                   deallocate(ANGGRD)
                   ANGGRD => angular_grid(J, GRID_PAR(J)%NANGA(L))
                end if
             endif
             angsize=size(ANGGRD,1)
             helpgrd(START:START+angsize-1,1:3)=ANGGRD(:,1:3)*RADGRD(I,1)
             helpgrd(START:START+angsize-1,4)=ANGGRD(:,4)* &
                  RADGRD(I,2)*16.0_R8_KIND*atan(1.0_R8_KIND)
             START=START+angsize
          enddo
          allocate(AGRID(J)%M(start-1,4),STAT=ALLOCSTAT)
          ASSERT (ALLOCSTAT==0)
          TOTGRD=>AGRID(J)%M        ! JUST FOR SIMPLIFICATION

          totgrd(1:start-1,:)=helpgrd(1:start-1,:)
          deallocate(helpgrd,stat=allocstat)
          ASSERT (allocstat==0)
       endif

       if (output_grid) then
          write(output_unit,'("The grid contains ",I6," points")') size(totgrd,1)
       endif

       totgrd(:,1:3)=totgrd(:,1:3)+spread(unique_atoms(J)%position(:,1),1,&
            size(totgrd,1))

       deallocate(RADGRD)
       deallocate(ANGGRD)
    enddo

    if (output_grid) write(output_unit,*) &
         '=================================================='
    if (output_grid) write(output_unit,*) &
         'Construction of the grid finished'
    if (output_grid) write(output_unit,*) &
         '=================================================='
  end subroutine grid_main_calc
  !--------------------------------------------------------------------------


  subroutine grid_write(unit)
    !
    ! Writes input read by  grid_read() concerning the construction of
    ! the grid. The input consits of one general namelist (grid). This
    ! namelist contains the logical  sym_reduce, which selects, if the
    ! grid  is reduced  due to  local  symmetry.  For  every unique  a
    ! namelist  follows,  which contains  nrad,  nang,  npart and  the
    ! atomic radius. This radius can  be specified directly in a.u. or
    ! the kind of radius can be selected.
    !
    use input_module
    use echo_input_module, only: start, real, flag, intg, strng, &
         stop, echo_level_full
    use operations_module, only: operations_scf, operations_echo_input_level
    implicit none
    integer(kind=i4_kind), intent(in) :: unit
    !** End of interface **************************************

    integer(kind=i4_kind) :: status,i
    character(len=26) :: header

    call start("GRID","GRID_WRITE",unit,operations_echo_input_level)
    call flag("SYM_REDUCE      ",sym_reduce_scf  ,df_sym_reduce      )
    call flag("WEIGHT_GRADS    ",weight_grads_scf,df_weight_grads    )
    call stop()

    do i=1,n_unique_atoms
       nrad=grid_par_scf(i)%nrad
       nang=grid_par_scf(i)%nang
       npart=grid_par_scf(i)%npart
       radius=grid_par_scf(i)%radius
       rslater = radius == rslater_f(int(abs(unique_atoms(i)%z)))
       if (rslater) radius = df_radius
       rpople = radius == rpople_f(int(abs(unique_atoms(i)%z)))
       if (rpople) radius = df_radius
       rionic = radius == rionic_f(int(abs(unique_atoms(i)%z)))
       if (rionic) radius = df_radius
       if (rslater) rslater=df_rslater

       write(header,'("GRIDATOM # unique atom",i4)') i
       call start(header,"GRID_WRITE",unit,operations_echo_input_level)
       call intg("NRAD   ",nrad   ,df_nrad   )
       call intg("NANG   ",nang   ,df_nang   )
       call intg("NPART  ",npart  ,df_npart  )
       call real("RADIUS ",radius ,df_radius )
       call flag("RSLATER",rslater,df_rslater)
       call flag("RPOPLE ",rpople ,df_rpople )
       call flag("RIONIC ",rionic ,df_rionic )
       call stop(keep_namelist=gridatom_found)
    enddo


    if (.not.operations_scf) then
       call write_to_output_units("grid_write: deallocating")
       if (allocated(grid_par_scf)) then
          deallocate(grid_par_scf,STAT=status)
          if (status.ne.0) call error_handler &
               ("grid_write: deallocation (1) failed ")
       end if
    endif

  end subroutine grid_write
  !--------------------------------------------------------------------------



  subroutine grid_write_ph(unit)
    !
    ! Writes input read by  grid_read() concerning the construction of
    ! the grid. The input consits of one general namelist (grid). This
    ! namelist contains the logical  sym_reduce, which selects, if the
    ! grid  is reduced  due to  local  symmetry.  For  every unique  a
    ! namelist  follows,  which contains  nrad,  nang,  npart and  the
    ! atomic radius. This radius can  be specified directly in a.u. or
    ! the kind of radius can be selected.
    !
    use input_module
    use echo_input_module, only: start, real, flag, intg, strng, &
         stop, echo_level_full
    use operations_module, only: operations_post_scf, &
                                 operations_echo_input_level, &
                                 operations_response
    use comm_module, only: comm_i_am_master
    implicit none
    integer(kind=i4_kind), intent(in) :: unit
    !** End of interface *****************************************

    integer(i4_kind) :: i
    character(len=29) :: header
    logical :: equal_gridatoms

    call start("GRID_PH","GRID_WRITE_PH",unit,operations_echo_input_level)
    call flag("SYM_REDUCE  ",sym_reduce_ph  ,df_sym_reduce_ph  )
    call flag("WEIGHT_GRADS",weight_grads_ph,df_weight_grads_ph)
    call stop()

    if ( .not.comm_i_am_master() .and. &
         & .not.(operations_post_scf .or. operations_response)) then
       write(unit,'(a/)')" # Namelists GRIDATOM_PH : not passed to the slaves"
       goto 1000
    endif

    if (.not.(operations_post_scf)) then
       equal_gridatoms=.true.
    else
       equal_gridatoms = equal_scf_and_ph_grids(gridatom_only=.true.)
    end if

    do i=1,n_unique_atoms
       nrad=grid_par_ph(i)%nrad
       nang=grid_par_ph(i)%nang
       npart=grid_par_ph(i)%npart
       radius=grid_par_ph(i)%radius
       rslater = radius == rslater_f(int(abs(unique_atoms(i)%z)))
       if (rslater) radius = df_radius_ph
       rpople = radius == rpople_f(int(abs(unique_atoms(i)%z)))
       if (rpople) radius = df_radius_ph
       rionic = radius == rionic_f(int(abs(unique_atoms(i)%z)))
       if (rionic) radius = df_radius_ph
       if (rslater) rslater=df_rslater_ph

       if (equal_gridatoms) then
          ! the  entire set  of GRIDATOM_PH  namelists may  be skipped
          ! because  grid_par_scf(*)%xxx equals  grid_par_ph(*)%xxx in
          ! all cases
          write(header,'("GRIDATOM_PH # unique atom",i4)') i
          call start(header,"GRID_WRITE_PH",unit,operations_echo_input_level)
          call intg("NRAD   ",nrad   ,nrad   )
          call intg("NANG   ",nang   ,nang   )
          call intg("NPART  ",npart  ,npart  )
          call real("RADIUS ",radius ,radius )
          call flag("RSLATER",rslater,rslater)
          call flag("RPOPLE ",rpople ,rpople )
          call flag("RIONIC ",rionic ,rionic )
          call stop(drop_namelist=.true.)
       else
          write(header,'("GRIDATOM_PH # unique atom",i4)') i
          call start(header,"GRID_WRITE_PH",unit,operations_echo_input_level)
          call intg("NRAD   ",nrad   ,df_nrad_ph   )
          call intg("NANG   ",nang   ,df_nang_ph   )
          call intg("NPART  ",npart  ,df_npart_ph  )
          call real("RADIUS ",radius ,df_radius_ph )
          call flag("RSLATER",rslater,df_rslater_ph)
          call flag("RPOPLE ",rpople ,df_rpople_ph )
          call flag("RIONIC ",rionic ,df_rionic_ph )
          call stop(keep_namelist=gridatom_ph_found)
       endif
    enddo

    1000 continue
  end subroutine grid_write_ph
  !--------------------------------------------------------------------------



  subroutine grid_read()
    !
    ! Reads  input concerning  the construction  of the  scf-grid. The
    ! input  consits of  one  general namelist  (grid). This  namelist
    ! contains the  logical sym_reduce, which selects, if  the grid is
    ! reduced  due to  local symmetry.   For every  unique  a namelist
    ! follows,  which  contains  nrad,  nang,  npart  and  the  atomic
    ! radius. This  radius can  be specified directly  in a.u.  or the
    ! kind of radius can be selected.
    !
    use input_module
    use options_module,only: update_hessian_iteration
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: allocstat,unit,status,i

    if(n_unique_atoms > 0)  then
       allocate(grid_par_scf(n_unique_atoms),STAT=allocstat)
       if(allocstat/=0) then
          call error_handler('Allocation failed in routine: grid_read')
       end if
    end if

    unit = input_intermediate_unit()

    sym_reduce = df_sym_reduce
    weight_grads = df_weight_grads
    if ( input_line_is_namelist("grid") ) then
       call input_read_to_intermediate
       read(unit, nml=grid, iostat=status)
       if (status .gt. 0) call input_error( &
            "grid_read: namelist grid")
    if(update_hessian_iteration.gt.0) weight_grads=.true.
    endif
    sym_reduce_scf = sym_reduce
    df_sym_reduce_ph = sym_reduce_scf
    weight_grads_scf = weight_grads
    df_weight_grads_ph = weight_grads_scf

    gridatom_found = input_line_is_namelist("gridatom")

    do i=1,n_unique_atoms
       nrad    = df_nrad
       nang    = df_nang
       npart   = df_npart
       radius  = df_radius
       rslater = df_rslater
       rpople  = df_rpople
       rionic  = df_rionic
       if ( gridatom_found ) then
          if ( .not. input_line_is_namelist("gridatom") ) call input_error( &
               "grid_read: namelist gridatom expected")
          call input_read_to_intermediate
          read(unit, nml=gridatom, iostat=status)
          if (status .gt. 0) call input_error( &
               "grid_read: namelist gridatom")
       endif
       if(npart/=1.and.npart/=5) call input_error &
            ('grid_read: Sorry, until now only npart=1 or npart=5 is implemented')
       if(nrad<0.or.nrad>max_nrad) call input_error &
            ('wrong nrad value specified')
       if(nang<0.or.nang>max_nang) call input_error &
            ('wrong nang value specified')
       if(npart<0.or.npart>max_npart) call input_error &
            ('wrong npart value specified')
       if(radius==df_radius.and..not.(rslater.or.rpople.or.rionic)) &
            rslater=.true.
       if(radius/=df_radius) then
          if(radius<=0.0_r8_kind.or.radius>=max_radius) &
               call input_error('Wrong value for radius')
          if(rslater.or.rpople.or.rionic) call input_error &
               ('Dont know what to do, multiple radius definition')
       endif
       if(rslater) then
          if(rpople.or.rionic) call input_error &
               ('Dont know what to do, multiple radius definition')
          radius=rslater_f(int(unique_atoms(i)%z))
       endif
       if(rpople) then
          if(rionic.or.rslater) call input_error &
               ('Dont know what to do, multiple radius definition')
          radius=rpople_f(int(unique_atoms(i)%z))
       endif
       if(rionic)then
          if(rpople.or.rslater) call input_error &
               ('Dont know what to do, multiple radius definition')
          radius=rionic_f(int(unique_atoms(i)%z))
       endif
       if(.not.radius.gt.0.0_r8_kind)then
          WARN('force radius = 1.0')
          radius = 1.0
       endif
       grid_par_scf(i)%nrad=nrad
       grid_par_scf(i)%nang=nang
       grid_par_scf(i)%npart=npart
       grid_par_scf(i)%radius=radius

    enddo

  end subroutine grid_read
  !---------------------------------------------------------------------------

  subroutine grid_read_ph(grid_ph_only)
    !
    ! Reads  input concerning  the  construction of  the ph-grid.  The
    ! input consits  of one general namelist  (grid_ph). This namelist
    ! contains the  logical sym_reduce, which selects, if  the grid is
    ! reduced  due to  local symmetry.   For every  unique  a namelist
    ! follows,  which  contains  nrad,  nang,  npart  and  the  atomic
    ! radius. This  radius can  be specified directly  in a.u.  or the
    ! kind of radius can be selected.
    !
    use input_module
    implicit none
    logical, intent(out) :: grid_ph_only
    !** End of interface *****************************************

    integer(kind=i4_kind) :: unit,status,i

    if(n_unique_atoms > 0) then
       allocate(grid_par_ph(n_unique_atoms),STAT=alloc_grid(2))
       ASSERT(alloc_grid(2).eq.0)
       alloc_grid(2)=1
    end if

    unit = input_intermediate_unit()

    sym_reduce = df_sym_reduce_ph
    weight_grads = df_weight_grads_ph
    if ( input_line_is_namelist("grid_ph") ) then
       call input_read_to_intermediate
       read(unit, nml=grid_ph, iostat=status)
       if (status .gt. 0) call input_error( &
            "grid_read: namelist grid_ph")
    endif
    weight_grads_ph = weight_grads
    sym_reduce_ph = sym_reduce

    gridatom_ph_found = input_line_is_namelist("gridatom_ph")

    if ( gridatom_ph_found ) then
    do i=1,n_unique_atoms
       nrad    = df_nrad_ph
       nang    = df_nang_ph
       npart   = df_npart_ph
       radius  = df_radius_ph
       rslater = df_rslater_ph
       rpople  = df_rpople_ph
       rionic  = df_rionic_ph
       !if ( gridatom_ph_found ) then ! moved outside the UA loop
          if ( .not. input_line_is_namelist("gridatom_ph") ) call input_error( &
               "grid_read_ph: namelist gridatom_ph expected")
          call input_read_to_intermediate
          read(unit, nml=gridatom_ph, iostat=status)
          if (status .gt. 0) call input_error( &
               "grid_read_ph: namelist gridatom")
       !endif ! moved down outside UA loop
       if(npart/=1.and.npart/=5) call input_error &
            ('grid_read: Sorry, until now only npart=1 or npart=5 is implemented')
       if(nrad<0.or.nrad>max_nrad_ph) call input_error &
            ('wrong nrad value specified')
       if(nang<0.or.nang>max_nang_ph) call input_error &
            ('wrong nang value specified')
       if(npart<0.or.npart>max_npart_ph) call input_error &
            ('wrong npart value specified')

       ! if not given in input, take assume they want rslater:
       if(radius==df_radius_ph.and..not.(rslater.or.rpople.or.rionic)) &
            rslater=.true.

       if(radius/=df_radius_ph) then
          if(radius<=0.0_r8_kind.or.radius>=max_radius_ph) &
               call input_error('Wrong value for radius')
          if(rslater.or.rpople.or.rionic) call input_error &
               ('Dont know what to do, multiple radius definition')
       endif

       if(rslater) then
          if(rpople.or.rionic) call input_error &
               ('Dont know what to do, multiple radius definition')
          radius=rslater_f(int(unique_atoms(i)%z))
       endif

       if(rpople) then
          if(rionic) call input_error &
               ('Dont know what to do, multiple radius definition')
          radius=rpople_f(int(unique_atoms(i)%z))
       endif

       if(rionic) radius=rionic_f(int(unique_atoms(i)%z))

       ASSERT(radius>0.0_r8_kind)
       grid_par_ph(i)%nrad=nrad
       grid_par_ph(i)%nang=nang
       grid_par_ph(i)%npart=npart
       grid_par_ph(i)%radius=radius

    enddo
    endif ! if ( gridatom_ph_found ) then

    grid_ph_only = .not. gridatom_ph_found
  end subroutine grid_read_ph
  !---------------------------------------------------------------------------


  subroutine grid_copy_to_ph(gridatom_only)
    !--------------------------------------------------------------------------
    ! Purpose : Copies setting for ordinary grid to post hov grid
    implicit none
    logical, optional :: gridatom_only ! default = .false.
    !** End of interface *****************************************

    logical :: full_copy
    integer(kind=i4_kind) :: i

    if (present(gridatom_only)) then
       full_copy = .not.gridatom_only
    else
       full_copy = .true.
    endif
    if (full_copy) then
       if(n_unique_atoms > 0) then
          allocate(grid_par_ph(n_unique_atoms),STAT=alloc_grid(2))
          ASSERT(alloc_grid(2).eq.0)
          alloc_grid(2)=1
!        print*,'alloc_grid(2) on in grid_copy_to_ph' !
       end if
       df_sym_reduce_ph = sym_reduce_scf
       sym_reduce_ph = df_sym_reduce_ph
       df_weight_grads_ph = weight_grads_scf
       weight_grads_ph = df_weight_grads_ph
       weight_grads = weight_grads_ph ! weight_grads is public
    endif
    do i=1,n_unique_atoms
       grid_par_ph(i)%nrad   = grid_par_scf(i)%nrad
       grid_par_ph(i)%nang   = grid_par_scf(i)%nang
       grid_par_ph(i)%npart  = grid_par_scf(i)%npart
       grid_par_ph(i)%radius = grid_par_scf(i)%radius
    enddo

  end subroutine grid_copy_to_ph
  !---------------------------------------------------------------------------


  function equal_scf_and_ph_grids(gridatom_only) result(equal)
    !
    ! Purpose : Checks if identical grid parameters are used in the
    !           SCF and in the post SCF part
    implicit none
    logical, optional :: gridatom_only ! default = .false.
    logical           :: equal
    !** End of interface *****************************************
    logical               :: check_all
    integer(kind=i4_kind) :: i

    if (present(gridatom_only)) then
       check_all = .not.gridatom_only
    else
       check_all = .true.
    endif

    equal = .true.
    if (check_all) then
       equal = sym_reduce_ph   .eqv. sym_reduce_scf ! .and. &
               ! for compatibility with older versions
               ! weight_grads_ph .eqv. weight_grads_ph ! not weight_grads_scf
    endif
    do i=1,n_unique_atoms
       equal = equal .and. grid_par_ph(i)%nrad   == grid_par_scf(i)%nrad   &
                     .and. grid_par_ph(i)%nang   == grid_par_scf(i)%nang   &
                     .and. grid_par_ph(i)%npart  == grid_par_scf(i)%npart  &
                     .and. grid_par_ph(i)%radius == grid_par_scf(i)%radius
    enddo

  end function equal_scf_and_ph_grids
  !---------------------------------------------------------------------------

  function rpople_f (z)
    !
    ! Atomic radii from pople in a.u.
    !
    implicit none
    integer (i4_kind), intent (in) :: z
    real (r8_kind) :: rpople_f
    !** End of interface *****************************************

    integer, parameter  :: p=R8_KIND  ! TO SIMPLIFY WRITING
    integer (i4_kind) :: i, zz

    ! atomic radii from pople in a.u.
    ! Slater Radii from Phys. Rev. 36 (1930),  p.57]
    integer, parameter :: NA = 92
    real (r8_kind), dimension (0:NA) :: RPOPLE_arr = [ &
         1.0_p, &               ! fake entry for Z < 1
         1.0000_p, 0.5882_p, &  ! H, He
         3.0769_p, 2.0513_p, 1.5383_p, 1.2308_p, 1.0256_p, 0.8791_p, 0.7692_p, &
         0.6838_p, 4.0909_p, 3.1579_p, 2.5714_p, 2.1687_p, 1.8750_p, 1.6514_p, 1.4754_p, &
         1.3333_p, 6.2227_p, 4.8035_p, 4.5633_p, 4.3460_p, 4.1484_p, 4.6407_p, 3.8028_p, &
         3.6507_p, 3.5103_p, 3.3803_p, 3.7000_p, 3.1471_p, 2.7380_p, 2.4230_p, 2.1730_p, &
         1.9698_p, 1.8013_p, 1.6594_p, 7.2727_p, 5.6140_p, 5.3333_p, 5.0794_p, 5.7143_p, &
         5.4237_p, 5.1613_p, 4.9231_p, 4.7059_p, 4.5070_p, 4.3243_p, 3.6782_p, 3.2000_p, &
         2.8319_p, 2.5397_p, 2.3022_p, 2.1053_p, 1.9394_p, 8.0182_p, 6.1895_p, 5.8800_p, &
         (6.1895_p, I=1, 6), 5.8800_p, (6.1895_p, I=1, 6), 5.8800_p, 5.6000_p, 5.3455_p, &
         5.1130_p, 4.9000_p, 4.7040_p, 4.5231_p, 4.9690_p, 4.7676_p, 4.0552_p, 3.5280_p, &
         3.1221_p, 2.8000_p, 2.5381_p, 2.3211_p, 2.1382_p, 8.4046_p, 6.4877_p, 6.1633_p, &
         5.8698_p, 6.1633_p, 6.1633_p] ! ..., U

    if (1 <= z .and. z <= NA) then
       zz = z
    else if (z > NA) then
       WARN("rpople_f: Z>92, use Z=92") ! FIXME: update msg with NA!
       zz = NA
    else
       WARN("rpople_f: Z out of range, force to 1.0")
       zz = 0
    end if

    rpople_f = rpople_arr(zz)
  end function rpople_f

  function rslater_f (z)
    !
    ! Slater radii in a. u.
    !
    use constants, only: angstrom
    implicit none
    integer (i4_kind), intent (in) :: z
    real (r8_kind) :: rslater_f
    !** End of interface *****************************************

    integer, parameter  :: p=R8_KIND  ! TO SIMPLIFY WRITING
    integer (i4_kind) :: i, zz

    !
    ! Slater radii in Angstroem.
    !
    ! Slater, J. C. (1964).  "Atomic Radii in Crystals".  Journal of
    ! Chemical Physics 41 (10): 3199–3205.
    ! http://dx.doi.org/10.1063/1.1725697
    !
    integer, parameter :: NA = 95
    real (r8_kind), dimension (0:NA) :: RSLATER_arr = [ &
         1 / angstrom, &        ! fake entry for Z < 1
         0.53_p, 0.25_p, &      ! H, He
         1.45_p, 1.05_p, 0.85_p, 0.70_p, 0.65_p, 0.60_p, 0.50_p, 0.50_p, &
         1.80_p, 1.50_p, 1.25_p, 1.10_p, 1.00_p, 1.00_p, 1.00_p, 1.00_p, 2.20_p, 1.80_p, &
         1.60_p, 1.40_p, 1.35_p, 1.40_p, 1.40_p, 1.40_p, 1.35_p, 1.35_p, 1.35_p, 1.35_p, &
         1.30_p, 1.25_p, 1.15_p, 1.15_p, 1.15_p, 1.15_p, 2.35_p, 2.00_p, 1.80_p, 1.55_p, &
         1.45_p, 1.45_p, 1.35_p, 1.30_p, 1.35_p, 1.40_p, 1.60_p, 1.55_p, 1.55_p, 1.45_p, &
         1.45_p, 1.40_p, 1.40_p, 1.40_p, 2.60_p, 2.15_p, 1.95_p, (1.85_p, i=1, 6), 1.80_p, &
         (1.75_p, i=1, 7), 1.55_p, 1.45_p, 1.35_p, 1.35_p, 1.30_p, 1.35_p, 1.35_p, 1.35_p, &
         1.50_p, 1.90_p, 1.80_p, 1.60_p, 1.90_p, 1.90_p, 1.90_p, 2.60_p, 2.15_p, &
         1.95_p, 1.80_p, 1.80_p, (1.75_p, i = 1, 4)] ! Ac, Th, Pa, U, Np, Pu, Am

    if (1 <= z .and. z <= NA) then
       zz = z
    else if (z > NA) then
       WARN("rslater_f: Z>95, use Z=95") ! FIXME: update msg with NA!
       zz = NA
    else
       WARN("rslater_f: Z out of range, force to 1.0")
       zz = 0
    end if

    rslater_f = rslater_arr(zz) * angstrom
  end function rslater_f

  function rionic_f (z)
    !
    ! Ionic radii in a. u.
    !
    implicit none
    integer (i4_kind), intent (in) :: z
    real (r8_kind) :: rionic_f
    !** End of interface *****************************************

    integer, parameter  :: p=R8_KIND  ! TO SIMPLIFY WRITING
    integer (i4_kind) :: zz

    ! Ionic radii in a. u.
    ! RD Shannon: Acta Cristal. A 32, Part5, 751-767 (1976)
    ! LH Ahrens: Geochim. Cosmochim. Acta 2, 155-169 (1952)
    ! WA Zachariasen: Structure Reports 13, 435 (1950)
    ! NN Greenwood, A Earnshaw: Chemistry of the elements,
    ! Pergamon, Oxford (1984)
    integer, parameter :: NA = 92
    real (r8_kind), dimension (0:NA) :: rionic_arr= [ &
         1.0_p, &               ! fake entry for Z < 1
         0.370_p, 1.400_p, &    ! H, He
         0.760_p, 0.450_p, 0.270_p, 0.160_p, 0.160_p, 1.400_p, 1.330_p, 1.540_p, &
         1.020_p, 0.720_p, 0.535_p, 0.400_p, 0.440_p, 1.840_p, 1.810_p, 1.890_p, &
         1.380_p, 1.000_p, 0.745_p, 0.860_p, 0.790_p, 0.800_p, 0.830_p, 0.780_p, 0.745_p, &
         0.690_p, 0.770_p, 0.740_p, 0.620_p, 0.730_p, 0.580_p, 1.980_p, 1.900_p, 2.020_p, &
         1.520_p, 1.180_p, 0.900_p, 0.720_p, 0.720_p, 0.690_p, 0.645_p, 0.680_p, 0.665_p, &
         0.860_p, 1.150_p, 0.950_p, 0.800_p, 0.930_p, 0.760_p, 2.210_p, 2.200_p, 2.160_p, &
         1.670_p, 1.350_p, 1.032_p, 1.010_p, 0.990_p, 0.983_p, 0.970_p, 0.958_p, &
         1.170_p, 0.938_p, 0.923_p, 1.070_p, 0.901_p, 0.890_p, 1.030_p, 1.020_p, &
         0.861_p, 0.710_p, 0.720_p, 0.660_p, 0.630_p, 0.630_p, 0.680_p, 0.800_p, &
         1.370_p, 1.190_p, 1.500_p, 1.190_p, 1.030_p, 2.300_p, 2.270_p, 2.400_p, &
         1.800_p, 1.430_p, 1.120_p, 1.080_p, 1.040_p, 1.025_p] ! ..., U

    if (1 <= z .and. z <= NA) then
       zz = z
    else if (z > NA) then
       WARN("rionic_f: Z>92, use Z=92") ! FIXME: update msg with NA!
       zz = NA
    else
       WARN("rionic_f: Z out of range, force to 1.0")
       zz = 0
    end if

    rionic_f = rionic_arr(zz)
  end function rionic_f

  logical function more_grid(n, p_points, p_weights)
    !
    ! Gives access  to the gridpoints,  which are kept in  the private
    ! array grdpts (1:N,  1:4).  The first argument "n"  is the number
    ! of requested  gridpoints.  This  procedure returns "true"  if it
    ! actually returns non-zero number of grid points.
    !
    ! When there are  no points left it returns  "false" once and then
    ! upon next invocation it will  proceed from the very beginning as
    ! usual.
    !
    ! This differs  from the behaviour of  earlier implementation that
    ! returned  "true"  on  the  last  batch  of  points  and  "false"
    ! otherwise.
    !
    ! It uses  DLB algorithm for dynamically work  stealing.  Thus the
    ! task indices which should be done by this worker are handed over
    ! from the DLB  routine dlb_give_more(). The only thing  to do for
    ! this routine  is to  direct the pointer  to the needed  data. As
    ! here the atoms  are not needed, they are  taken from the general
    ! grdpts. Note  that grdpts(:,  :) now keeps  the whole  grid, not
    ! only the part statically assigned to the current processor as in
    ! earlier implementation.
    !
    !-------------------------------------------------------------------------
    use dlb, only: dlb_give_more, idlb_kind
    implicit none
    integer(i4_kind), intent(in) :: n ! number of requested gridpoints
    real(r8_kind), dimension(:, :), pointer :: p_points
    real(r8_kind), dimension(:), pointer :: p_weights
    !** End of interface *****************************************

    integer(idlb_kind) :: jobs(2)
    integer(idlb_kind) :: n_dlb

    !
    ! True if the job interval is not empty:
    !
    n_dlb = n
    more_grid = dlb_give_more(n_dlb, jobs)

    if (more_grid) then
      p_points  => grdpts(jobs(1)+1:jobs(2), 1:3)
      p_weights => grdpts(jobs(1)+1:jobs(2), 4)
    else
      p_points  => NULL()
      p_weights => NULL()
    endif
  end function more_grid


  !-----------------------------------------------------------------------
  logical function more_grid_atom(n, unique_atom, p_points, p_weights)
    !
    ! Gives access  to the gridpoints,  which are kept in  the private
    ! array grdpts. n is the number of requested gridpoints. more_grid
    ! returns  a  logical, which  is  .true.  if  it actually  returns
    ! non-zero number of grid points.
    !
    ! When there are  no points left it returns  "false" once and then
    ! upon next invocation it will  proceed from the very beginning as
    ! usual.
    !
    ! This differs  from the behaviour of  earlier implementation that
    ! returned  "true"  on  the  last  batch  of  points  and  "false"
    ! otherwise.
    !
    ! It uses  dlb algorithm for  dynamically work stealing.  Thus the
    ! job numbers  which should be done  by this proc  are handed over
    ! from the  DLB routine  dlb_give_more. The only  thing to  do for
    ! this routine  is to  direct the pointer  to the needed  data. As
    ! here the atoms are needed, they are taken from the agrid.
    !
    !------------------------------------------------------------------
    use dlb, only: dlb_give_more_color, long => idlb_kind
    implicit none
    integer (i4_kind), intent (in) :: n ! number of requested gridpoints
    integer (i4_kind), intent (out) :: unique_atom ! number of atom
    real (r8_kind), dimension (:, :), pointer :: p_points
    real (r8_kind), dimension (:), pointer   :: p_weights
    !** End of interface *****************************************

    integer (long) :: n_dlb     ! number of requested gridpoints
    integer (long) :: color     ! number of atom in other format
    integer (long) :: jobs(2)

    n_dlb = n
    more_grid_atom = dlb_give_more_color(n_dlb, color, jobs)
    unique_atom = int (color, kind = i4_kind) ! interface obligation

    if (more_grid_atom) then
      p_points  => agrid(unique_atom)%m(jobs(1)+1:jobs(2), 1:3)
      p_weights => agrid(unique_atom)%m(jobs(1)+1:jobs(2), 4)
    else
      p_points  => NULL()
      p_weights => NULL()
    endif
  end function more_grid_atom

  !-----------------------------------------------------------------------

  subroutine grid_bcast (post_scf)
    !
    ! Purpose: broadcast structure grid_par to the slaves
    !
    implicit none
    logical, optional, intent (in) :: post_scf
    !** End of interface *******************************************************

  end subroutine grid_bcast

  !-----------------------------------------------------------------------

  function RADIALGRID(NRADTYPE,NRADPTS,RCUR) RESULT(P)
    !
    ! This  routine  calculates the  radial  part  of the  integration
    ! grid. In the moment only the  radialgrid of Pople ( Gill et al.,
    ! Chem.  Phys. Lett. 209,506(1993)  ) is  available. By  using the
    ! switch NRADTYPE other radial grids may be implemented
    !
    use TYPE_MODULE
    implicit none
    real(R8_KIND), dimension(:,:), pointer :: P ! on exit p points on
    ! the radialgrid, p(:,1) are the radii, p(:,2) the weights
    integer, intent(IN) :: NRADTYPE, NRADPTS ! NRADTYPE is a switch to
    ! select different radialgrids, NRADTYPE contains the number of
    ! radial gridpoints. It is calculated in the routine set_grid_par
    ! out of NRAD and the radius of the unique atom. Note: On exit the
    ! number of gridpoints does not have to be equal to nradpts due to
    ! a cut off. (cmp.; RAD_LIMIT)
    real(R8_KIND), intent(IN) :: RCUR ! length which scales the shape of
    ! the radialgriid. Normally it is the atomic radius.
    !** End of interface *****************************************

    real(kind=R8_KIND),parameter    :: RAD_LIMIT= 40.0  ! cut off border

    real(kind=R8_KIND)     :: R,RHN,RHNI
    integer(kind=I4_KIND)  :: ALLOC_STAT,I,II

    ASSERT(RCUR>0.0_r8_kind)

    ! implemented only for Pople radial grid:
    ASSERT(NRADTYPE.EQ.0)

    if (output_grid) write (output_unit,*) 'Pople radial grid selected'
    if (output_grid) write (output_unit,*) 'Radius:', RCUR

    RHN = NRADPTS + 1
    do I = NRADPTS, 1, -1
       RHNI = RHN - I
       R = RCUR * (I / RHNI)**2

       if ( R < RAD_LIMIT ) exit
    enddo

    allocate(P(I, 2), STAT=ALLOC_STAT)
    ASSERT (ALLOC_STAT==0)

    if (output_grid) write(output_unit,*) 'The radialgrid contains', size(P,1), 'points'

    P(I, 1) = R
    P(I, 2) = 2 * RCUR**3 * RHN * (I / RHNI)**5 / RHNI**2
    II=I
    do I = II - 1, 1, -1
       RHNI = RHN - I
       P(I, 1) = RCUR * (I / RHNI)**2
       P(I, 2) = 2 * RCUR**3 * RHN * (I / RHNI)**5 / RHNI**2
    enddo
  end function RADIALGRID

  !--------------------------------------------------------------------------

  function angular_grid(J, ICGRD) result(pts)
    !
    ! This  routine calculates  the  angular part  of the  integration
    ! grid,  using   the  formulas  from  Lebedev   and  McLaren.  The
    ! orientation of  the grid is also  performed. In order  to do so,
    ! the rotation  matrix ro is  used. If the variable  sym_reduce is
    ! .true. points  which are  equivalent to ann  other point  due to
    ! symmetry reasons are eliminated
    !
    ! Warning: returns allocated pointer, beware of memory leaks!
    !
    use symmetry_element, only: sym_element_local ! for prune()
    implicit none
    real(r8_kind), pointer       :: pts(:, :) ! on exit
    ! pts points on the angular grid; pts(:, 1:3) contains  the
    ! coordinates of the gridpoints, pts(:, 4) the weight
    integer(i4_kind), intent(in) :: J      ! unique atom index
    integer(i4_kind), intent(in) :: ICGRD  ! order of the grid
    ! *** end of interface ***

    real(r8_kind) :: ro(3, 3) ! rotation matrix
    real(r8_kind), pointer :: p(:, :)

    ! Lebedev-McLaren grid:
    p => INIGRD(ICGRD)

    if (output_grid) write(output_unit,*) 'CHEKSUM: (should be close to one)', sum(p)

    if (sym_reduce) then
      ! sets sym_element_local and returns suitable rotation matrix:
      call get_symmetry(J, ro)
    else
      ro = rotmat(J)
    endif

    ! ro is not used if weight_grads = T, complicates derivatives:
    if(.not.weight_grads)then
      p(:, 1:3) = rotate(ro, p(:, 1:3))
    endif

    ! eliminating the gridpoints which are equivalent due to symmetry:
    if (sym_reduce) then
      ! uses sym_element_local that was set by get_symmetry():
      pts => prune(sym_element_local, p)

      deallocate(p)
    else
      pts => p
    endif

    if (output_grid) write(output_unit,*) 'Number of points per shell is',size(pts, 1)
    if (output_grid) write(output_unit,*) '=============================================='
  end function angular_grid

  function INIGRD(ICGRD) result(pts)
    !
    ! Purpose: This routine calculates the angular part of the
    !          integration grid, using the formulas from Lebedev and
    !          McLaren.
    !
    ! Warning: returns allocated pointer, beware of memory leaks!
    !
    implicit none
    real(r8_kind), pointer         :: pts(:, :) ! on exit pts points
    ! on the angular grid; pts(:, 1:3) contains the coordinates of the
    ! gridpoints, pts(:, 4) the weight
    integer(i4_kind), intent(IN)   :: icgrd  ! switch, that determinates
      ! the angular order of the grid
    !** End of interface *****************************************

    ! help variables
    real(kind=r8_kind), pointer :: p(:, :)

    real(kind=R8_KIND) :: A1(6,3),A2(12,3),A3(8,3),B1(24,3),B2(24,3),B3(24,3),&
         B4(24,3),B5(24,3),B6(24,3),C1(24,3),C2(24,3),D1(48,3),D2(48,3)
    real(kind=R8_KIND) :: WA1,WA2,WA3,WB1,WB2,WB3,WB4,WB5,WB6,WC1,WC2,WD1,WD2,&
         TAU,V1,V2,M1,M2,M3,M4,M5,M6,U1,U2,W1,W2
    integer  :: NPTS,ALLOC_STAT

    ! end of help variables

    NPTS=0

    if (ICGRD==31) then

       !  A.D.McLaren, Math. Comp. 17, 361 (1963).

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING MCLAREN GRID WITH L = 3, N = 6'
       if (output_grid) write (output_unit,*)

       allocate(p(6,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            (' MEMORY ALLOCATION FAILED IN SU INIGRD')


       WA1 = 1.0_r8_kind/6.0_r8_kind

       call GETA1(A1)

       call ADDPTS(NPTS, 6,A1,WA1,P)

       elseif (ICGRD.EQ.32) then

          !       A.D.McLaren, Math. Comp. 17, 361 (1963).

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING MCLAREN GRID WITH L = 3, N = 8'
       if (output_grid) write (output_unit,*)

       allocate(p(8,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA2 = 1.0_r8_kind/8.0_r8_kind

       call GETA3(A3)

       call ADDPTS(NPTS, 8,A3,WA2,P)


       elseif (ICGRD.EQ.33) then

          !       A.D.McLaren, Math. Comp. 17, 361 (1963).

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING MCLAREN GRID WITH L = 3, N = 12'
       if (output_grid) write (output_unit,*)

       allocate(p(12,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA3 = 1.0_r8_kind/12.0_r8_kind

       call GETA2(A2)

       call ADDPTS(NPTS,12,A2,WA3,P)


       elseif (ICGRD.EQ.51) then

          !       A.D.McLaren, Math. Comp. 17, 361 (1963).

       allocate(p(14,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING MCLAREN GRID WITH L = 5, N = 14'
       if (output_grid) write (output_unit,*)

       WA1 = 8.0_r8_kind/120.0_r8_kind
       WA3 = 9.0_r8_kind/120.0_r8_kind

       call GETA1(A1)
       call GETA3(A3)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)


       elseif (ICGRD.EQ.71) then

          !       A.D.McLaren, Math. Comp. 17, 361 (1963).

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING MCLAREN GRID WITH L = 7, N = 26'
       if (output_grid) write (output_unit,*)

       allocate(p(26,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA1 = 40.0_r8_kind/840.0_r8_kind
       WA2 = 32.0_r8_kind/840.0_r8_kind
       WA3 = 27.0_r8_kind/840.0_r8_kind

       call GETA1(A1)
       call GETA2(A2)
       call GETA3(A3)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS,12,A2,WA2,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)


       elseif (ICGRD.EQ.91) then

       allocate(p(38,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 9, N = 38'
       if (output_grid) write (output_unit,*)

       WA1 = 1.0_r8_kind/105.0_r8_kind
       WA3 = 9.0_r8_kind/280.0_r8_kind
       WC1 = 1.0_r8_kind/35.0_r8_kind

       V1 = 1.0_r8_kind/6.0_r8_kind

       call GETA1(A1)
       call GETA3(A3)
       call GETC1(V1,C1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,C1,WC1,P)

       elseif (ICGRD.EQ.111) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 11, N = 50'
       if (output_grid) write (output_unit,*)

       allocate(p(50,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA1 = 4.0_r8_kind/315.0_r8_kind
       WA2 = 64.0_r8_kind/2835.0_r8_kind
       WA3 = 27.0_r8_kind/1280.0_r8_kind
       WB1 = 14641.0_r8_kind/725760.0_r8_kind

       M1 = 3.0_r8_kind/sqrt(11.0_r8_kind)

       call GETA1(A1)
       call GETA2(A2)
       call GETA3(A3)
       call GETBK(M1,B1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS,12,A2,WA2,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)

       elseif (ICGRD.EQ.131) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) 'USING THE LEBEDEV GRID WITH L = 13, N = 78'
       if (output_grid) write (output_unit,*)

       allocate(p(78,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0)  call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA1 = 0.0138665921047_r8_kind
       WB1 = 0.0130509318626_r8_kind
       WB2 = 0.0132064232231_r8_kind
       WC1 = 0.0119426635549_r8_kind

       M1 = 0.914152532416_r8_kind
       M2 = 0.359236381200_r8_kind
       V1 = 0.206339716383_r8_kind

       call GETA1(A1)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETC1(V1,C1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,C1,WC1,P)

       elseif (ICGRD.EQ.151) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 15, N = 86'
       if (output_grid) write (output_unit,*)

       allocate(p(86,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       TAU = 2.0_R8_KIND*(20.0_r8_kind-sqrt(13.0_r8_kind))/43.0_r8_kind

       !        WA1 = 16.0D+00*(TAU*TAU-20.0D+00*TAU+22.0D+00)/
       !     &           (3465.0D+00*(3.0_R8_KIND*TAU-4.0_R8_KIND)**2)
       !        WA3 = -243.0D+00*(27.0D+00*TAU*TAU-20.0D+00*TAU-4.0_R8_KIND)/
       !     &           (24640.0D+00*(3.0_R8_KIND*TAU-4.0_R8_KIND)**2)
       WA1 = 8.0D+00/693.0D+00
       WA3 = 0.0119439090859_r8_kind
       WB1 = 0.0111105557106_r8_kind
       WB2 = 0.0118765012945_r8_kind
       WC1 = 0.0118123037469_r8_kind

       M1 = 0.852518311701_r8_kind
       M2 = 0.189063552885_r8_kind
       V1 = 0.120441650315_r8_kind

       call GETA1(A1)
       call GETA3(A3)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETC1(V1,C1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,C1,WC1,P)


       elseif (ICGRD.EQ.171) then

       allocate(p(110,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 17, N = 110'
       if (output_grid) write (output_unit,*)

       WA1 = 0.00382827049494_r8_kind
       WA3 = 0.00979373751249_r8_kind
       WB1 = 0.00821173728319_r8_kind
       WB2 = 0.00959547133607_r8_kind
       WB3 = 0.00994281489118_r8_kind
       WC1 = 19652.0_r8_kind/2027025.0_r8_kind

       M1 = 0.965124035087_r8_kind
       M2 = 0.828769981253_r8_kind
       M3 = 0.215957291846_r8_kind
       V1 = 3.0_R8_KIND/17.0_r8_kind

       call GETA1(A1)
       call GETA3(A3)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETBK(M3,B3)
       call GETC1(V1,C1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,B3,WB3,P)
       call ADDPTS(NPTS,24,C1,WC1,P)

       elseif (ICGRD.EQ.191) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 19, N = 146'
       if (output_grid) write (output_unit,*)

       allocate(p(146,4),STAT=ALLOC_STAT)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA1 = 1856.0_r8_kind/3095235.0_r8_kind
       WA2 = 606208.0_r8_kind/82219995.0_r8_kind
       WA3 = 6490935.0_r8_kind/900204032.0_r8_kind
       WB1 = 7.57439415905E-03_r8_kind
       WB2 = 6.75382948631E-03_r8_kind
       WB3 = 7.11635549312E-03_r8_kind
       WD1 = 1773593.0_r8_kind/253693440.0_r8_kind

       M1 = 0.974888643677_r8_kind
       M2 = 0.807089818360_r8_kind
       M3 = 0.291298882210_r8_kind
       U1 = 0.140355381171_r8_kind
       W1 = 0.449332832327_r8_kind

       call GETA1(A1)
       call GETA2(A2)
       call GETA3(A3)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETBK(M3,B3)
       call GETDK(U1,W1,D1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS,12,A2,WA2,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,B3,WB3,P)
       call ADDPTS(NPTS,48,D1,WD1,P)

       elseif (ICGRD.EQ.231) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 23, N = 194'
       if (output_grid) write (output_unit,*)

       allocate(p(194,4),STAT=alloc_stat)

       if(alloc_stat/=0) call error_handler&
            ('MEMORY ALLOCATION FAILED IN SU INIGRD')

       !        WA1 = 9344.0D+00/5242545.0D+00
       !        WA2 = 27246592.0D+00/4765968207.0D+00
       !        WA3 = 59049.0D+00*1599797.0D+00/(15323648.0D+00*1106105.0D+00)
       WA1 = 1.78234044724E-03_r8_kind
       WA2 = 5.71690594998E-03_r8_kind
       WA3 = 5.57338317884E-03_r8_kind
       WB1 = 5.51877146727E-03_r8_kind
       WB2 = 5.15823771181E-03_r8_kind
       WB3 = 5.60870408259E-03_r8_kind
       WB4 = 4.10677702817E-03_r8_kind
       WC1 = 5.05184606462E-03_r8_kind
       WD1 = 5.53024891623E-03_r8_kind
       !       WC1 = 2085136.0_r8_kind/412747335.0_r8_kind
       !       WD1 = 231173.0_r8_kind*231173.0_r8_kind/&
       !(1580544.0_r8_kind*6113965.0_r8_kind)

       M1 = 0.777493219315_r8_kind
       M2 = 0.912509096867_r8_kind
       M3 = 0.314196994183_r8_kind
       M4 = 0.982972302707_r8_kind
       V1 = 2.00_r8_kind/19.0_r8_kind
       U1 = 0.159041710538_r8_kind
       W1 = 0.525118572443_r8_kind

       call GETA1(A1)
       call GETA2(A2)
       call GETA3(A3)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETBK(M3,B3)
       call GETBK(M4,B4)
       call GETC1(V1,C1)
       call GETDK(U1,W1,D1)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS,12,A2,WA2,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,B3,WB3,P)
       call ADDPTS(NPTS,24,B4,WB4,P)
       call ADDPTS(NPTS,24,C1,WC1,P)
       call ADDPTS(NPTS,48,D1,WD1,P)

       elseif (ICGRD.EQ.291) then

       if (output_grid) write (output_unit,*)
       if (output_grid) write (output_unit,*) ' USING THE LEBEDEV GRID WITH L = 29, N = 302'
       if (output_grid) write (output_unit,*)

       allocate(p(302,4),STAT=alloc_stat)
       if(alloc_stat/=0) &
            call error_handler('MEMORY ALLOCATION FAILED IN SU INIGRD')

       WA1 = 8.54591172878E-04_r8_kind
       WA3 = 3.59911928502E-03_r8_kind
       WB1 = 3.65004580768E-03_r8_kind
       WB2 = 3.60482260142E-03_r8_kind
       WB3 = 3.57672966173E-03_r8_kind
       WB4 = 3.44978842429E-03_r8_kind
       WB5 = 3.10895312238E-03_r8_kind
       WB6 = 2.35210141366E-03_r8_kind
       WC1 = 3.60082093222E-03_r8_kind
       WC2 = 2.98234496317E-03_r8_kind
       WD1 = 3.57154055427E-03_r8_kind
       WD2 = 3.39231220501E-03_r8_kind

       M1 = 0.129238672710_r8_kind
       M2 = 0.371034178385_r8_kind
       M3 = 0.743452042987_r8_kind
       M4 = 0.867643624544_r8_kind
       M5 = 0.949454317226_r8_kind
       M6 = 0.990705621379_r8_kind

       V1 = 0.220093335298_r8_kind
       V2 = 0.650272754659E-01_r8_kind
       U1 = 0.800072749407_r8_kind
       W1 = 0.544867737258_r8_kind
       U2 = 0.902442529533_r8_kind
       W2 = 0.412772408317_r8_kind

       call GETA1(A1)
       call GETA3(A3)
       call GETBK(M1,B1)
       call GETBK(M2,B2)
       call GETBK(M3,B3)
       call GETBK(M4,B4)
       call GETBK(M5,B5)
       call GETBK(M6,B6)
       call GETC1(V1,C1)
       call GETC1(V2,C2)
       call GETDK(U1,W1,D1)
       call GETDK(U2,W2,D2)

       call ADDPTS(NPTS, 6,A1,WA1,P)
       call ADDPTS(NPTS, 8,A3,WA3,P)
       call ADDPTS(NPTS,24,B1,WB1,P)
       call ADDPTS(NPTS,24,B2,WB2,P)
       call ADDPTS(NPTS,24,B3,WB3,P)
       call ADDPTS(NPTS,24,B4,WB4,P)
       call ADDPTS(NPTS,24,B5,WB5,P)
       call ADDPTS(NPTS,24,B6,WB6,P)
       call ADDPTS(NPTS,24,C1,WC1,P)
       call ADDPTS(NPTS,24,C2,WC2,P)
       call ADDPTS(NPTS,48,D1,WD1,P)
       call ADDPTS(NPTS,48,D2,WD2,P)

    endif

    ! output pointer:
    pts => p
  contains

    subroutine ADDPTS(NPTS,N,A,WA,FF)
      implicit none
      integer,intent(IN) ::N
      integer,intent(INOUT) ::NPTS
      real(kind=r8_kind),intent(INOUT),dimension(:,:) ::FF
      real(kind=r8_kind),intent(IN),dimension(:,:)    ::A
      real(kind=r8_kind),intent(IN)     ::  WA


      FF(NPTS+1:NPTS+N,1:3)=A
      FF(NPTS+1:NPTS+N,4)=WA

      NPTS=NPTS+N

    end subroutine ADDPTS

    subroutine GETA1(A1)
      implicit none

      real(kind=R8_KIND),intent(OUT)    ::A1(6,3)
      real(kind=R8_KIND),parameter      ::ONE=1.0_R8_KIND,NULL=0.0_R8_KIND


      A1 = NULL

      A1(1,3) =  ONE
      A1(2,3) = -ONE
      A1(3,2) =  ONE
      A1(4,2) = -ONE
      A1(5,1) =  ONE
      A1(6,1) = -ONE

      return

    end subroutine GETA1

    subroutine GETA2(A2)
      implicit none

      real(kind=R8_KIND),intent(OUT)    :: A2(12,3)
      real(kind=R8_KIND),parameter      :: ONE=1.0_R8_KIND,NULL=0.0_R8_KIND

      real(kind=R8_KIND)                :: RSQ2


      RSQ2=ONE/sqrt(2.0_R8_KIND)

      A2=NULL

      A2(1,1)  =  RSQ2
      A2(2,1)  =  RSQ2
      A2(1,2)  =  RSQ2
      A2(2,2)  = -RSQ2
      A2(3,1)  = -RSQ2
      A2(3,2)  =  RSQ2
      A2(4,1)  = -RSQ2
      A2(4,2)  = -RSQ2

      A2(5,1)  =  RSQ2
      A2(5,3)  =  RSQ2
      A2(6,1)  =  RSQ2
      A2(6,3)  = -RSQ2
      A2(7,1)  = -RSQ2
      A2(7,3)  =  RSQ2
      A2(8,1)  = -RSQ2
      A2(8,3)  = -RSQ2

      A2(9,2) =  RSQ2
      A2(9,3) =  RSQ2
      A2(10,2) =  RSQ2
      A2(10,3) = -RSQ2
      A2(11,2) = -RSQ2
      A2(11,3) =  RSQ2
      A2(12,2) = -RSQ2
      A2(12,3) = -RSQ2

      return

    end subroutine GETA2

    subroutine GETA3(A3)
      implicit none

      real(kind=R8_KIND),intent(OUT)    :: A3(8,3)
      real(kind=R8_KIND),parameter      :: ONE=1.0_R8_KIND,NULL=0.0_R8_KIND

      real(kind=R8_KIND)                :: RSQ3

      integer                           :: N,I,J,K


      RSQ3=ONE/sqrt(3.0_R8_KIND)

      N = 0

      do I=-1,1,2
         do J=-1,1,2
            do K=-1,1,2
               N = N + 1
               A3(N,1) = I*RSQ3
               A3(N,2) = J*RSQ3
               A3(N,3) = K*RSQ3
            enddo
         enddo
      enddo

      return

    end subroutine GETA3

    subroutine GETBK(MK,BK)
      implicit none

      real(kind=R8_KIND),intent(OUT)    :: BK(24,3)
      real(kind=R8_KIND),intent(IN)     :: MK

      real(kind=R8_KIND)                :: LK
      real(kind=R8_KIND),parameter      :: ONE=1.0_R8_KIND,NULL=0.0_R8_KIND

      integer                           :: I,J,K,N1,N2,N3


      LK=sqrt((ONE-MK*MK)/2.0_R8_KIND)

      N1 = 0
      N2 = 8
      N3 = 16

      do I=-1,1,2
         do J=-1,1,2
            do K=-1,1,2

               N1 = N1 + 1
               N2 = N2 + 1
               N3 = N3 + 1

               BK(N1,1) = I*LK
               BK(N1,2) = J*LK
               BK(N1,3) = K*MK

               BK(N2,1) = I*LK
               BK(N2,2) = J*MK
               BK(N2,3) = K*LK

               BK(N3,1) = I*MK
               BK(N3,2) = J*LK
               BK(N3,3) = K*LK

            enddo
         enddo
      enddo

    end subroutine GETBK

    subroutine GETC1(V,C1)
      implicit none

      real(kind=R8_KIND),intent(OUT)    :: C1(24,3)
      real(kind=R8_KIND),intent(IN)     :: V

      real(kind=R8_KIND)                :: D,P1,Q1
      real(kind=R8_KIND),parameter      :: ONE=1.0_R8_KIND


      D=ONE-4.0_R8_KIND*V

      D=sqrt(D)

      P1=sqrt((ONE+D)/2.0_R8_KIND)
      Q1=sqrt((ONE-D)/2.0_R8_KIND)

      C1=0.0_R8_KIND
      C1(1,1) =  P1
      C1(1,2) =  Q1
      C1(2,1) =  P1
      C1(2,2) = -Q1
      C1(3,1) = -P1
      C1(3,2) =  Q1
      C1(4,1) = -P1
      C1(4,2) = -Q1

      C1(5,1) =  P1
      C1(5,3) =  Q1
      C1(6,1) =  P1
      C1(6,3) = -Q1
      C1(7,1) = -P1
      C1(7,3) =  Q1
      C1(8,1) = -P1
      C1(8,3) = -Q1

      C1(9,2) =  P1
      C1(9,3) =  Q1
      C1(10,2) =  P1
      C1(10,3) = -Q1
      C1(11,2) = -P1
      C1(11,3) =  Q1
      C1(12,2) = -P1
      C1(12,3) = -Q1

      C1(13,1) =  Q1
      C1(13,2) =  P1
      C1(14,1) =  Q1
      C1(14,2) = -P1
      C1(15,1) = -Q1
      C1(15,2) =  P1
      C1(16,1) = -Q1
      C1(16,2) = -P1

      C1(17,1) =  Q1
      C1(17,3) =  P1
      C1(18,1) =  Q1
      C1(18,3) = -P1
      C1(19,1) = -Q1
      C1(19,3) =  P1
      C1(20,1) = -Q1
      C1(20,3) = -P1

      C1(21,2) =  Q1
      C1(21,3) =  P1
      C1(22,2) =  Q1
      C1(22,3) = -P1
      C1(23,2) = -Q1
      C1(23,3) =  P1
      C1(24,2) = -Q1
      C1(24,3) = -P1

      return

    end subroutine GETC1

    subroutine GETDK(UK,WK,DK)
      implicit none

      real(kind=R8_KIND),intent(OUT)    :: DK(48,3)
      real(kind=R8_KIND),intent(IN)     :: UK,WK

      real(kind=R8_KIND)                :: RK

      real(kind=R8_KIND),parameter      :: ONE=1.0_R8_KIND
      integer                           :: I,J,K,N1,N2,N3


      RK=sqrt(ONE-UK*UK-WK*WK)

      N1 = 0
      N2 = 8
      N3 = 16

      do I=-1,1,2
         do J=-1,1,2
            do K=-1,1,2

               N1 = N1 + 1
               N2 = N2 + 1
               N3 = N3 + 1

               DK(N1,1) = I*UK
               DK(N1,2) = J*WK
               DK(N1,3) = K*RK

               DK(N2,1) = I*WK
               DK(N2,2) = J*RK
               DK(N2,3) = K*UK

               DK(N3,1) = I*RK
               DK(N3,2) = J*UK
               DK(N3,3) = K*WK

            enddo
         enddo
      enddo
      DK(25:48,1)=DK(1:24,1)
      DK(25:48,2)=DK(1:24,3)
      DK(25:48,3)=DK(1:24,2)
      return

    end subroutine GETDK

  end function INIGRD

  !---------------------------------------------------------------------

  function rotate(ro, pts) result(rpts)
    !
    ! Rotate atomic grid with 3x3 ro-matrix, leave weights intact.
    !
    implicit none
    real(r8_kind), intent(in) :: ro(3, 3)
    real(r8_kind), intent(in) :: pts(:, :) ! (npts, 3)
    real(r8_kind)             :: rpts(size(pts,1), size(pts, 2))
    ! *** end of interface ***

    integer(i4_kind) :: i, k

    do i = 1, size(pts, 1)
      do k = 1, 3
        rpts(i, k) = ro(k, 1) * pts(i, 1) &
                   + ro(k, 2) * pts(i, 2) &
                   + ro(k, 3) * pts(i, 3)
      enddo
      ! copy the weights in case they were provided:
      rpts(i, 4:) = pts(i, 4:)
    enddo
  end function rotate

  !---------------------------------------------------------------------

  function prune(sym_element_local, p) result(p1)
    !
    ! Eliminate the gridpoints which are equivalent due to symmetry.
    !
    ! WARINING: returns an allocated pointer, beware of memory leaks!
    !
    use symmetry_element, only: s_elements ! type
    implicit none
    type(s_elements), intent(in) :: sym_element_local
    real(r8_kind),    intent(in) :: p(:, :) ! (npts, 4)
    real(r8_kind), pointer       :: p1(:, :) ! (npts1, 4), npts1 <= npts
    ! *** end of interface ***

    integer(i4_kind)   :: alloc_stat
    integer(i4_kind)   :: i, j, k, l
    real(kind=r8_kind) :: grid_temp(size(p,1), size(p,2))
    real(kind=r8_kind) :: x_trans(3), x_trans2(3)
    real(kind=r8_kind) :: ro2(3, 3)
    real(kind=r8_kind) :: c, s
    real(kind=r8_kind) :: zweipi

    zweipi = 8.0_r8_kind * atan(1.0_r8_kind)

    if (output_grid) write(output_unit,*) '=========================================='
    if (output_grid) write(output_unit,*) &
         'Starting reduction of the angular grid due to symmetry'

    grid_temp(:, 1:4) = p(:, 1:4)

    ! mirrorplanes
    do i = 1, sym_element_local%num_sigma
       do j = 1, size(p,1)
          if ( grid_temp(j,4) /= 0.0_r8_kind ) then
             x_trans = p(j, 1:3) - two * sym_element_local%sigma(:, i)&
                  & * sum(p(j, 1:3) * sym_element_local%sigma(:, i))
             do k = j + 1, size(p,1)
                if ( sum((x_trans - p(k, 1:3))**2) <= 1.0D-14 ) then
                   grid_temp(j, 4) = grid_temp(j, 4) + grid_temp(k, 4)
                   grid_temp(k, 4) = 0.0_r8_kind
                endif
             end do
          end if
       end do
    end do

    ! rotation axes
    do i = 1, sym_element_local%num_axis
       ro2 = rotmat2(sym_element_local%rotaxis(:, i))
       do j = 1, size(p,1)
          if ( grid_temp(j,4) /= 0.0_r8_kind ) then
             x_trans = MATMUL(p(j, 1:3), ro2)
             s = sin(zweipi / sym_element_local%n_axis(i))
             c = cos(zweipi / sym_element_local%n_axis(i))
             do k = 1, sym_element_local%n_axis(i) - 1
                x_trans2(1) = c * x_trans(1) + s * x_trans(2)
                x_trans2(2) = c * x_trans(2) - s * x_trans(1)
                x_trans2(3) = x_trans(3)
                x_trans = MATMUL(ro2, x_trans2)
                do l = j + 1, size(p,1)
                   if ( sum((x_trans - p(l, 1:3))**2 ) <= 1.0D-12 ) then
                      grid_temp(j,4) = grid_temp(j,4) + grid_temp(l,4)
                      grid_temp(l,4) = 0.0_r8_kind
                   endif
                end do
                x_trans = x_trans2
             end do
          end if
       end do
    end do

    ! count grid points with non-zero weights:
    j=0
    do i = 1, size(p, 1)
      if ( grid_temp(i, 4) == 0.0_r8_kind ) cycle
      j = j + 1
      grid_temp(j, 1:3) = p(i, 1:3)
      grid_temp(j, 4) = grid_temp(i, 4)
    end do

    ! we return an allocated pointer from here:
    allocate(p1(j, 4), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    p1 = grid_temp(1:j, 1:4)
  end function prune

  !---------------------------------------------------------------------

  function atomicweight(atnum, x) result(w)
    !-------------------------------------------------------------------
    !  Purpose: Calculates the weight of every gridpoint according to
    !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988))
    !-------------------------------------------------------------------
    implicit none
    integer(i4_kind), intent(in)    :: atnum ! index of atom to which current grid points belong
    real(r8_kind),    intent(in)    :: x(:, :)
    real(r8_kind)                   :: w(size(x, 1))
    !** End of interface *****************************************

    call atomicweight2(atnum, x, grid_par(:)%radius, w)
  end function atomicweight

  function atomicweight_and_grad(atnum, x, graw) result(w)
    !-------------------------------------------------------------------
    !  Purpose: Calculates the weight of every gridpoint according to
    !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988)) and
    !           the gradient of the weight with respect to
    !-------------------------------------------------------------------
    implicit none
    integer(i4_kind), intent(in)    :: atnum ! index of atom to which current grid points belong
    real(r8_kind),    intent(in)    :: x(:, :)
    type(arrmat3),    intent(inout) :: graw(:)
    real(r8_kind)                   :: w(size(x, 1))
    !** End of interface *****************************************

    call atomicweight2(atnum, x, grid_par(:)%radius, w, graw)
  end function atomicweight_and_grad

  function atomicweight_and_dervs(atnum, x, graw, dervsw) result(w)
    !-------------------------------------------------------------------
    !  Purpose: Calculates the weight of every gridpoint according to
    !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988)) and
    !           the gradient of the weight with respect to
    !-------------------------------------------------------------------
    implicit none
    integer(i4_kind), intent(in)    :: atnum ! index of atom to which current grid points belong
    real(r8_kind),    intent(in)    :: x(:, :)
    type(arrmat3),    intent(inout) :: graw(:)
    type(arrmat5),    intent(inout) :: dervsw(:,:)
    real(r8_kind)                   :: w(size(x, 1))
    !** End of interface *****************************************

    call atomicweight2(atnum, x, grid_par(:)%radius, w, graw, dervsw)
  end function atomicweight_and_dervs

  subroutine atomicweight2(ua, x, radii, awt, graw, derw)
    !  Purpose: Calculates the weight of every gridpoint according to
    !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988))
    implicit none
    integer(i4_kind), intent(in)    :: ua      ! grid hosting atom
    real(r8_kind)   , intent(in)    :: x(:,:)  ! (npts,1:3) ! FIXME: reverse!
    real(r8_kind)   , intent(in)    :: radii(:) ! (nua)
    real(r8_kind)   , intent(out)   :: awt(:)   ! (npts), weights of (ua,1) at x
    type(arrmat3)   , intent(inout) :: graw(:) ! (nua)%m(npts,3,nea), d awt / d R
    type(arrmat5)   , intent(inout) :: derw(:,:) ! (nua,nua)%m(npts,3,nea,3,neb), d2 awt / d Ra d Rb
    optional :: graw, derw
    ! *** end of interface ***

    integer(i4_kind) :: npts, natm
    integer(i4_kind) :: i,j,k,ia
    integer(i4_kind) :: u1,e1,u2,e2
    real(r8_kind), allocatable :: atm(:,:) ! (3,natm)    ! atomic coords
    real(r8_kind), allocatable :: rad(:)   ! (natm)      ! atomic radii
    real(r8_kind), allocatable :: dst(:,:) ! (natm,natm) ! interatomic distances
    real(r8_kind), allocatable :: art(:,:) ! (natm,natm) ! atomic radii ``ratios''

    real(r8_kind), allocatable :: gwt(:,:,:)! (npts,3,natm) ! temp for grads of weight
    real(r8_kind), allocatable :: dwt(:,:,:,:,:)! (npts,3,3,natm,natm) ! temp for dervs of weight

    integer :: memstat

    FPP_TIMER_START(awt)

    ! number of (grid) points:
    npts = size(x,1)

    ! number of atoms:
    natm = sum(UNIQUE_ATOMS(:)%N_EQUAL_ATOMS)

    allocate(atm(3,natm),rad(natm),STAT=memstat)
    ASSERT(memstat==0)

    if(present(graw))then
      allocate(gwt(npts,3,natm),STAT=memstat)
      ASSERT(memstat==0)
    endif

    if(present(derw))then
      ASSERT(present(graw))
      allocate(dwt(npts,3,3,natm,natm),STAT=memstat)
      ASSERT(memstat==0)
    endif

    ! collect atomic coordiantes and radii into linear arrays:
    ia = 0
    k  = 0
    do u1=1,size(UNIQUE_ATOMS)
      do e1=1,UNIQUE_ATOMS(u1)%N_EQUAL_ATOMS
        k = k + 1
        atm(:,k) = UNIQUE_ATOMS(u1)%POSITION(:,e1)
        rad(k)   = radii(u1)
        ! note the index of the grid host:
        if( u1==ua .and. e1==1 ) ia=k
      enddo
    enddo
    ASSERT(ia>0)

    ! precompute interatomic distances:
    allocate(dst(natm,natm))
    do j=1,natm
      do i=1,j
        dst(i,j) = sqrt(sum( (atm(:,i)-atm(:,j))**2 ))
        dst(j,i) = dst(i,j)
      enddo
    enddo

    ! precompute atomic radii ratios:
    allocate(art(natm, natm))
    if(adjust_cell_size)then
      do j=1,natm
        do i=1,j
          art(i,j) = ( rad(i) - rad(j) ) / ( rad(i) + rad(j) )
          art(j,i) = - art(i,j)
        enddo
      enddo
    else
      ! as if all grid radii were equal:
      art(:, :) = 0.0
    endif

    !
    ! the atm variable seems to be "optimized" away
    ! by intel compiler at -O2, even though it is
    ! used within "contained" functions.
    !
    ! From now on we pass everything explicitly
    ! to the contained functions.
    !

    if(present(derw))then
      FPP_TIMER_START(awt2)

      ! calculate becke weights/grads/dervs for atom ``ia'' at each grid point:
      do k=1,npts
        call becke2(ia, x(k,1:3), awt(k), gwt(k,:,:), dwt(k,:,:,:,:) & ! dwt(npts,3,3,natm,natm)
                   , natm, atm, dst, art)
                   ! pass explicitly, instead of host association
      enddo

      FPP_TIMER_STOP(awt2)
    else if(present(graw))then
      FPP_TIMER_START(awt1)

      ! calculate becke weights/grads for atom ``ia'' at each grid point:
      do k=1,npts
        call becke1(ia, x(k,1:3), awt(k), gwt(k,:,:) & ! gwt(npts,3,natm)
                   , natm, atm, dst, art)
                   ! pass explicitly, instead of host association
      enddo

      FPP_TIMER_STOP(awt1)
    else
      FPP_TIMER_START(awt0)

      ! calculate becke weights for atom ``ia'' at each grid point:
      do k=1,npts
        call becke0(ia, x(k,1:3), awt(k) &
                   , natm, atm, dst, art)
                   ! pass explicitly, instead of host association
      enddo

      FPP_TIMER_STOP(awt0)
    endif

    ! copy gradients into legacy datastructure:
    if(present(graw))then
      FPP_TIMER_START(awt1)
      k  = 0
      do u1=1,size(UNIQUE_ATOMS)
        do e1=1,UNIQUE_ATOMS(u1)%N_EQUAL_ATOMS
          k = k + 1
          graw(u1)%m(:npts,:,e1) = gwt(:,:,k) ! gwt(npts,3,natm)
        enddo
      enddo
      FPP_TIMER_STOP(awt1)
    endif

    ! copy derivatives into legacy datastructure:
    if(present(derw))then
      FPP_TIMER_START(awt2)
      i  = 0
      do u1=1,size(UNIQUE_ATOMS)
        do e1=1,UNIQUE_ATOMS(u1)%N_EQUAL_ATOMS
          i = i + 1
          k = 0
          do u2=1,size(UNIQUE_ATOMS)
            do e2=1,UNIQUE_ATOMS(u2)%N_EQUAL_ATOMS
              k = k + 1
              derw(u1,u2)%m(:npts,:,e1,:,e2) = dwt(:,:,:,i,k) ! dwt(npts,3,3,natm,natm)
            enddo
          enddo
        enddo
      enddo
      FPP_TIMER_STOP(awt2)
    endif

    deallocate(atm,dst,rad,STAT=memstat)
    ASSERT(memstat==0)
    if(allocated(gwt))then
      deallocate(gwt,STAT=memstat)
      ASSERT(memstat==0)
    endif
    if(allocated(dwt))then
      deallocate(dwt,STAT=memstat)
      ASSERT(memstat==0)
    endif

    deallocate(art)

    FPP_TIMER_STOP(awt)

  contains

    subroutine becke0(a, x, w, natm, atm, dst, art)
      !  Purpose: Calculates the weight of a gridpoint according to
      !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988))
      use becke_step_func, only: mu0
      implicit none
      integer(i4_kind), intent(in)  :: a        ! grid hosting atom == Voronoj cell index
      real(r8_kind)   , intent(in)  :: x(:)     ! (1:3), grid point
      real(r8_kind)   , intent(out) :: w        ! weight of Voronoj cell "a" at "x"
      ! were host associated from above, now passed explicitly:
      integer(i4_kind), intent(in)  :: natm     ! number of atoms
      real(r8_kind)   , intent(in)  :: atm(:,:) ! atm(:,:) ! (1:3,natm), atomic coords
      real(r8_kind)   , intent(in)  :: dst(:,:) ! dst(:,:) ! (natm,natm), precomputed interatomic distances
      real(r8_kind)   , intent(in)  :: art(:,:) ! art(:,:) ! (natm,natm), precomputed atomic radii ratios
      ! *** end of interface ***

      integer(i4_kind) :: i,k
      real(r8_kind)    :: v(3)
      real(r8_kind)    :: r(natm)
      real(r8_kind)    :: sik,ski
      real(r8_kind)    :: p(natm)
      real(r8_kind)    :: mu

      ! compute distances to atomic centers:
      do i=1,natm
        v    = atm(:,i) - x
        r(i) = sqrt(sum(v**2))

        ! a neutral value:
        p(i) = 1
      enddo

      ! loop over atomic pairs:
      do i=1,natm
        do k=1,i-1

          ! compute elliptic coordinates of a grid point:
          mu = ( r(i) - r(k) ) / dst(k,i)
          ! note that mu(i,k) = - mu(k,i):

          ! account for different atomic sizes,
          ! increase the cell (region with mu<0) for big (art>0) atoms:
          mu = mu0(mu, art(i, k))

          ! compute the pair cutoff function:
          sik = (1 - mu) / 2
          ski = (1 + mu) / 2

          ! accumulate the products of the pair cutoff functions:
          p(i) = p(i) * sik
          p(k) = p(k) * ski
        enddo
      enddo

      ! compute the weight for atom ``a'':
      w = p(a) / sum(p)
    end subroutine becke0

    subroutine becke1(a, x, w, w1, natm, atm, dst, art)
      !  Purpose: Calculates the weight of a gridpoint according to
      !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988))
      use becke_step_func, only: mu1
      implicit none
      integer(i4_kind), intent(in)  :: a       ! grid hosting atom
      real(r8_kind)   , intent(in)  :: x(:)    ! (1:3), grid point
      real(r8_kind)   , intent(out) :: w       ! resulting weight
      real(r8_kind)   , intent(out) :: w1(:,:) ! (3,natm) weight gradients
      ! were host associated from above, now passed explicitly:
      integer(i4_kind), intent(in)  :: natm     ! number of atoms
      real(r8_kind)   , intent(in)  :: atm(:,:) ! atm(:,:) ! (1:3,natm), atomic coords
      real(r8_kind)   , intent(in)  :: dst(:,:) ! dst(:,:) ! (natm,natm), precomputed interatomic distances
      real(r8_kind)   , intent(in)  :: art(:,:) ! art(:,:) ! (natm,natm), precomputed atomic radii ratios
      ! *** end of interface ***

      real(r8_kind), parameter :: eps=1.0e-50_r8_kind ! offset denominators

      integer(i4_kind) :: i,j,k
      real(r8_kind)    :: v(3,natm)
      real(r8_kind)    :: r(natm)
      real(r8_kind)    :: s(natm,natm,0:1)
      real(r8_kind)    :: p(natm), sump
      real(r8_kind)    :: p1(3,natm)
      real(r8_kind)    :: mu(natm,natm)
      real(r8_kind)    :: m(0:1)
      real(r8_kind)    :: m1(3)

      ! compute distances to atomic centers:
      do i=1,natm
        v(:,i) = atm(:,i) - x
        r(i)   = sqrt(sum(v(:,i)**2))
        ASSERT(r(i)>0)
      enddo

      ! loop over atomic pairs:
      do j=1,natm
          mu(j,j) = 0.0_r8_kind
        do i=1,j-1

          ! compute elliptic coordinate of a grid point:
          mu(i,j) = ( r(i) - r(j) ) / dst(i,j)
          mu(j,i) = - mu(i,j)

          ! compute the partition fanction for the pair (i,j):
          ! account for different atomic sizes:
          m(0:1) = mu1(mu(i, j), art(i, j))

          ! compute the pair cutoff function:
          s(i,j,0) = ( 1 - m(0) ) / 2
          s(j,i,0) = ( 1 + m(0) ) / 2
          ! NOTE: s(-m) = 1 - s(m)

          ! logarithmic derivative s(1)/s(0), s(m) is monotonousely decreasing:
          s(i,j,1) = - m(1) / ( 2 * s(i,j,0) + eps )
          s(j,i,1) = - m(1) / ( 2 * s(j,i,0) + eps )
          ! NOTE: d/dm ( 1 - s(-m) ) = ds/dm
        enddo
        ! some neutral values:
        s(i,i,0) = 1.0_r8_kind
        s(i,i,1) = 0.0_r8_kind
      enddo

      ! compute the products of the pair functions:
      do i=1,natm
        p(i) = product(s(i,:,0)) ! there is an extra factor s(i,i)
      enddo

      ! normalization factor:
      sump = sum(p(:)) ! all Ps have an extra factor s(i,i)

      ! scale partition functions so that they sum to one:
      do i=1,natm
        p(i) = p(i) / sump
      enddo

      ! return only the weight for atom ``a'':
      w = p(a) ! / sump

      ! compute derivatives dW(a)/dR(ib) of a-partition function wrt all atoms k:
      w1(:,a) = 0.0 ! gradient wrt grid host a is computed using translational invariance
      do k=1,natm
        ! working with one k-gradient at a time...
        if( k == a ) cycle ! computed using d/da = - sum(k/=a) d/dk

        ! we will need derivatives dP(i)/dk of all partition functions P(i):
        p1(:,k) = 0 ! will accumulate sum over i/=k
        do i=1,natm
          if( i == k  ) cycle ! dP(k)/dk is special

          ! dMU(i,k)/dRk:
          m1(:)  = - ( v(:,k)/r(k) - mu(i,k) * ( v(:,i) - v(:,k) ) / dst(i,k) ) / dst(i,k)

          ! dP(i)/dRk, i/=k :
          p1(:,i) = p(i) * s(i,k,1) * m1(:)

          ! dP(k)/dRk has contributions from all i/=k:
          p1(:,k) = p1(:,k) &
                  - p(k) * s(k,i,1) * m1(:) ! FIXME: s(k,i,1) ?
        enddo
        ! derivative of all P(i) wrt Rk are ready in p1(:,i) ...

        ! derivatives dW(a)/dk for k/=a:
        w1(:,k) = p1(:,a) ! / sump
        do i=1,natm
          w1(:,k) = w1(:,k) - w * p1(:,i) ! / sump
        enddo
!       w1(:,k) = 0.0
!       do i=1,natm
!         if( i == a ) cycle
!         w1(:,k) = w1(:,k) + ( p(i) * p1(:,a) &
!                             - p(a) * p1(:,i) )
!       enddo

        ! derivative dW(a)/da:
        w1(:,a) = w1(:,a) - w1(:,k) ! d/da = - sum(k/=a) d/dk
      enddo
    end subroutine becke1

    subroutine becke2(a, x, w, w1, w2, natm, atm, dst, art)
      !  Purpose: Calculates the weight of a gridpoint according to
      !           the scheme of Becke ( J. Chem. Phys. 88,2547 (1988))
      use becke_step_func, only: mu2, muab
      implicit none
      integer(i4_kind), intent(in)  :: a       ! grid hosting atom
      real(r8_kind)   , intent(in)  :: x(:)    ! (1:3), grid point
      real(r8_kind)   , intent(out) :: w       ! resulting weight
      real(r8_kind)   , intent(out) :: w1(:,:) ! (3,natm) weight gradients
      real(r8_kind)   , intent(out) :: w2(:,:,:,:) ! (3,3,natm,natm) weight 2-derivatives
      ! were host associated from above, now passed explicitly:
      integer(i4_kind), intent(in)  :: natm     ! number of atoms
      real(r8_kind)   , intent(in)  :: atm(:,:) ! atm(:,:) ! (1:3,natm), atomic coords
      real(r8_kind)   , intent(in)  :: dst(:,:) ! dst(:,:) ! (natm,natm), precomputed interatomic distances
      real(r8_kind)   , intent(in)  :: art(:,:) ! art(:,:) ! (natm,natm), precomputed atomic radii ratios
      ! *** end of interface ***

      real(r8_kind), parameter :: eps=1.0e-50_r8_kind ! offset denominators

      integer(i4_kind) :: i,k,l
      integer(i4_kind) :: s,t
      real(r8_kind)    :: v(3,natm)
      real(r8_kind)    :: r(natm)
      real(r8_kind)    :: sik(0:2), ski(0:2)
      real(r8_kind)    :: p(natm), sump
      real(r8_kind)    :: p1(3,natm,natm), sump1(3,natm)
      real(r8_kind)    :: p2(3,3,natm,natm,natm), sump2
      real(r8_kind)    :: mu!(natm,natm)
      real(r8_kind)    :: m(0:2)
!     real(r8_kind)    :: m1(3)
      real(r8_kind)    :: mi(3),mk(3),mii(3,3),mik(3,3),mkk(3,3)
!     real(r8_kind)    :: na(3),nb(3),naa(3,3),nab(3,3),nbb(3,3)

      ! compute distances to atomic centers:
      do i=1,natm
        v(:,i) = atm(:,i) - x
        r(i)   = sqrt(sum(v(:,i)**2))
        ASSERT(r(i)>0)
      enddo

      ! initialize arrays with neutral values:
      p(:) = 1
      p1(:,:,:) = 0
      p2(:,:,:,:,:) = 0

      ! loop over atomic pairs:
      do k=1,natm
        do i=1,k-1

          ! compute elliptic coordinate of a grid point:
          mu = ( r(i) - r(k) ) / dst(i,k)

          ! compute the partition fanction for the pair (i,k):
          ! account for different atomic sizes:
          m(0:2) = mu2(mu, art(i, k))

          ! compute the pair cutoff function:
          sik(0) = ( 1 - m(0) ) / 2
          ski(0) = ( 1 + m(0) ) / 2
          ! NOTE: s(m) + s(-m) = 1

          sik(1) = - m(1) / 2
          ski(1) =   m(1) / 2
          ! NOTE: in fact, s(m) is monotonousely decreasing both at m=mu and m=-mu.
          !       However, we technically treat both as functions of single mu!

          sik(2) = - m(2) / 2
          ski(2) =   m(2) / 2
          ! NOTE: We technically treat both as functions of single mu!
          !       Although, s(m) + s(-m) = 1 --- where one is convex another is concave.

          call muab(v(:,i),v(:,k),mi,mk,mii,mik,mkk)

          do t=1,3 ! (s,t) are Cartesian components
          do s=1,3
            ! scale second derivatives:
            p2(s,t,:,:,i) = p2(s,t,:,:,i)             * sik(0)
            p2(s,t,:,:,k) = p2(s,t,:,:,k)             * ski(0)

            ! more terms in the rows/columns i and k (li,il,lk,kl, for all l):
            do l=1,natm
            p2(s,t,l,i,i) = p2(s,t,l,i,i) + p1(s,l,i) * sik(1) * mi(t)
            p2(s,t,l,k,i) = p2(s,t,l,k,i) + p1(s,l,i) * sik(1) * mk(t)

            p2(s,t,i,l,i) = p2(s,t,i,l,i) + p1(t,l,i) * sik(1) * mi(s)
            p2(s,t,k,l,i) = p2(s,t,k,l,i) + p1(t,l,i) * sik(1) * mk(s)

            p2(s,t,l,i,k) = p2(s,t,l,i,k) + p1(s,l,k) * ski(1) * mi(t)
            p2(s,t,l,k,k) = p2(s,t,l,k,k) + p1(s,l,k) * ski(1) * mk(t)

            p2(s,t,i,l,k) = p2(s,t,i,l,k) + p1(t,l,k) * ski(1) * mi(s)
            p2(s,t,k,l,k) = p2(s,t,k,l,k) + p1(t,l,k) * ski(1) * mk(s)
            enddo

            ! more terms in the derivatives of p(i) wrt  ik:
            p2(s,t,i,i,i) = p2(s,t,i,i,i) +    p(i) * ( sik(2) * mi(s) * mi(t) &
                                                      + sik(1) * mii(s,t) )
            p2(s,t,i,k,i) = p2(s,t,i,k,i) +    p(i) * ( sik(2) * mi(s) * mk(t) &
                                                      + sik(1) * mik(s,t) )
            p2(s,t,k,i,i) = p2(s,t,k,i,i) +    p(i) * ( sik(2) * mk(s) * mi(t) &
                                                      + sik(1) * mik(t,s) )
            p2(s,t,k,k,i) = p2(s,t,k,k,i) +    p(i) * ( sik(2) * mk(s) * mk(t) &
                                                      + sik(1) * mkk(s,t) )

            ! scale second derivatives:

            ! more terms in the derivatives of p(k) wrt  ki:
            p2(s,t,i,i,k) = p2(s,t,i,i,k) +    p(k) * ( ski(2) * mi(s) * mi(t) &
                                                      + ski(1) * mii(s,t) )
            p2(s,t,i,k,k) = p2(s,t,i,k,k) +    p(k) * ( ski(2) * mi(s) * mk(t) &
                                                      + ski(1) * mik(s,t) )
            p2(s,t,k,i,k) = p2(s,t,k,i,k) +    p(k) * ( ski(2) * mk(s) * mi(t) &
                                                      + ski(1) * mik(t,s) )
            p2(s,t,k,k,k) = p2(s,t,k,k,k) +    p(k) * ( ski(2) * mk(s) * mk(t) &
                                                      + ski(1) * mkk(s,t) )
          enddo
          enddo

          ! scale first derivatives of p(i), older values have to be used above for p2:
          p1(:,:,i) = p1(:,:,i) * sik(0)

          ! one more term in first derivatives of p(i/k) wrt i and k:
          p1(:,i,i) = p1(:,i,i)          + p(i) * sik(1) * mi(:)
          p1(:,k,i) = p1(:,k,i)          + p(i) * sik(1) * mk(:)

          ! scale first derivatives of p(k):
          p1(:,:,k) = p1(:,:,k) * ski(0)

          ! one more term in first derivatives of p(k) wrt i and k:
          p1(:,i,k) = p1(:,i,k)          + p(k) * ski(1) * mi(:)
          p1(:,k,k) = p1(:,k,k)          + p(k) * ski(1) * mk(:)

          ! update values of the product, older values to be used above for p1:
          p(i) = p(i) * sik(0)
          p(k) = p(k) * ski(0)
        enddo
      enddo

      ! normalization factor:
      sump  = sum(p(:))

      ! return only the weight for atom ``a'':
      w = p(a) / sump

      ! ... and its gradients:
      do k=1,natm
        sump1(:,k) = sum(p1(:,k,:),2)
        w1(:,k)    = p1(:,k,a) / sump - w * sump1(:,k) / sump
      enddo

      ! ... and the second derivatives:
      do k=1,natm
      do i=1,natm
        do t=1,3
        do s=1,3
          sump2   =   sum(p2(s,t,i,k,:))
          w2(s,t,i,k) = ( p2(s,t,i,k,a)      &
                      - w * sump2            &
                      - w1(s,i) * sump1(t,k) &
                      - w1(t,k) * sump1(s,i) &
                      ) / sump

        enddo
        enddo
      enddo
      enddo

      ! The grid-point x is attached to the grid host a:
      ! so that we need to substitute:
      !
      ! d/da -> d/da + d/dx = d/da - sum d/dk = - sum d/dk
      !                                 k          k/=a
      !

      ! Weight w(a) is translationally invariant, use this for d/da:
      w1(:,a)    = 0 ! contributes to the sum in next stmt
      w1(:,a)    = - sum(w1(:,:),2)    ! d/da = - sum(k/=a) d/dk

      w2(:,:,:,a) = 0
      w2(:,:,:,a) = - sum(w2(:,:,:,:),4)
      w2(:,:,a,:) = 0
      w2(:,:,a,:) = - sum(w2(:,:,:,:),3)
    end subroutine becke2
  end subroutine atomicweight2

  subroutine SET_GRID_PAR()
    !---------------------------------------------------------------------
    !   SETS GRID_PARAMETER WHICH ARE NEEDED FOR BUILDING THE GRID, NAMELY
    !   THE ATOMIC RADII, AND POPLES ALPHAPARAMETER IF NEEDED
    !-----------------------------------------------------------------------
    !** End of interface *****************************************
    implicit none

    integer :: I,J,K
    real(kind=r8_kind),dimension(4)  :: ALPHADAT


    do J=1,N_UNIQUE_ATOMS

       if(GRID_PAR(J)%NPART==0)  GRID_PAR(J)%NPART=4
       K=int(UNIQUE_ATOMS(J)%Z+0.5)
       GRID_PAR(J) % NRADPTS = int (GRID_PAR(J) % NRAD * GRID_PAR(J) % RADIUS)

       if(GRID_PAR(J)%NPART/=1) then
          write(output_unit,*) 'NPART:',GRID_PAR(J)%NPART,&
               ' ,PRUNING FOR UNIQUE ',J&
               ,'ACTIVATED'
          if(K<=2) then          ! DETERMINATION OF POPLE ALPHA VALUES
             ALPHADAT(1) = 0.25_r8_kind
             ALPHADAT(2) = 0.50_r8_kind
             ALPHADAT(3) = 1.00_r8_kind
             ALPHADAT(4) = 4.50_r8_kind
          elseif (K<=10) then
             ALPHADAT(1) = 0.16666666666666666666667_r8_kind
             ALPHADAT(2) = 0.500_r8_kind
             ALPHADAT(3) = 0.900_r8_kind
             ALPHADAT(4) = 3.500_r8_kind
          elseif (K<=18) then
             ALPHADAT(1) = 0.100_r8_kind
             ALPHADAT(2) = 0.400_r8_kind
             ALPHADAT(3) = 0.800_r8_kind
             ALPHADAT(4) = 2.500_r8_kind
          elseif (K<=36) then
             ALPHADAT(1) = 0.025_r8_kind
             ALPHADAT(2) = 0.200_r8_kind
             ALPHADAT(3) = 0.500_r8_kind
             ALPHADAT(4) = 2.500_r8_kind
          elseif (K<=54) then
             ALPHADAT(1) = 0.015_r8_kind
             ALPHADAT(2) = 0.040_r8_kind
             ALPHADAT(3) = 0.200_r8_kind
             ALPHADAT(4) = 2.500_r8_kind
          elseif (K<=86) then
             ALPHADAT(1) = 0.005_r8_kind
             ALPHADAT(2) = 0.015_r8_kind
             ALPHADAT(3) = 0.040_r8_kind
             ALPHADAT(4) = 2.500_r8_kind
          else
             ALPHADAT(1) = 0.003_r8_kind
             ALPHADAT(2) = 0.008_r8_kind
             ALPHADAT(3) = 0.020_r8_kind
             ALPHADAT(4) = 2.000_r8_kind
          endif
          GRID_PAR(J)%ALPHA=ALPHADAT
          do I = 1, 4
             GRID_PAR(J) % NANGA(I) = 31 + int (((I - 1.0d0) / 3) * (GRID_PAR(J) % NANG - 31))
          enddo
          GRID_PAR(J) % NANGA(5) = GRID_PAR(J) % NANGA(3)
       endif
    enddo
    !---------------------------------------------------------------------
  end subroutine SET_GRID_PAR
  !---------------------------------------------------------------------


  function rotmat(I)
    !---------------------------------------------------------------------
    !  PURPOSE: BUILDING AB THE ROTATIONMATRIX WHICH IS NEEDED FOR
    !           THE CORRECT ORIENTATION OF THE GRID
    !           THE Z-AXIS OF THE ATOMIC GRID IS THE CONNECTION LINE OF
    !           THE CENTER OF ATOMIC CHARGE AND THE ATOMIC POSITION
    !           THIS PROVIDES TRANSLATION AND ROTATIONAL INVARIANCE OF THE GRID
    !
    !  INPUT PARAMETER:  I, NUMBER OF THE UNIQUE ATOM
    !
    !  SUBROUTINE USES CHARGE_CENTER AND UNIQUE_ATOMS(I)%POSITION
    !---------------------------------------------------------------------
    implicit none
    real(kind=R8_KIND),dimension(3,3) :: ROTMAT
    integer,intent(IN)              :: I  ! UNIQUE ATOM
    !** End of interface *****************************************
    real(kind=R8_KIND)              :: ZDIST,ODIST,ZSIGN
    real(kind=R8_KIND),dimension(3) :: DIF !COORDINATES RELATIVE TO&
         !  CHARGE_CENTER

    ROTMAT=reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))

    DIF=UNIQUE_ATOMS(I)%POSITION(:,1)-CHARGE_CENTER



    ZDIST= DIF(1)*DIF(1)+DIF(2)*DIF(2)

    if(ZDIST<1.0D-12.AND.DIF(3)<0.0_R8_KIND) ROTMAT(3,3)=-1.0_R8_KIND

    if(ZDIST>=1.0D-12) then

       ZSIGN = 1.0_R8_KIND
       if(DIF(3).LT.0) ZSIGN = -1.0_R8_KIND
       ODIST = ZDIST + DIF(3)*DIF(3)
       ODIST = sqrt(ODIST)
       ZDIST = sqrt(ZDIST)
       ROTMAT(1,1) = ZSIGN*DIF(1)*DIF(3)/ZDIST/ODIST
       ROTMAT(2,1) = ZSIGN*DIF(2)*DIF(3)/ZDIST/ODIST
       ROTMAT(3,1) = -ZSIGN*ZDIST/ODIST
       ROTMAT(1,2) = -DIF(2)/ZDIST
       ROTMAT(2,2) = DIF(1)/ZDIST
       ROTMAT(3,2) = 0.0_R8_KIND
       ROTMAT(1,3) = DIF(1)/ODIST
       ROTMAT(2,3) = DIF(2)/ODIST
       ROTMAT(3,3) = DIF(3)/ODIST
    endif
    !---------------------------------------------------------------------
  end function rotmat



  function ROTMAT2(z)
    !---------------------------------------------------------------------
    !  PURPOSE: BUILDING AB THE ROTATIONMATRIX WHICH IS NEEDED FOR
    !           THE CORRECT ORIENTATION OF THE GRID
    !           THE Z-AXIS OF THE ATOMIC GRID IS THE CONNECTION LINE OF
    !           THE CENTER OF ATOMIC CHARGE AND THE ATOMIC POSITION
    !           THIS PROVIDES TRANSLATION AND ROTATIONAL INVARIANCE OF THE GRID
    !
    !  INPUT PARAMETER:  z, coordinates of the z-axis
    !
    !  SUBROUTINE USES CHARGE_CENTER AND UNIQUE_ATOMS(I)%POSITION
    !---------------------------------------------------------------------
    implicit none
    real(kind=R8_KIND),dimension(3,3) :: rotmat2
    real(kind=r8_kind),dimension(3) :: z  ! UNIQUE ATOM
    !** End of interface *****************************************
    real(kind=R8_KIND)              :: ZDIST,ODIST,ZSIGN
    real(kind=R8_KIND),dimension(3) :: DIF !COORDINATES RELATIVE &
         !  TO CHARGE_CENTER

    rotmat2=reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))

    DIF=z


    ZDIST= DIF(1)*DIF(1)+DIF(2)*DIF(2)

    if(ZDIST<1.0D-12.AND.DIF(3)<0.0_R8_KIND) ROTMAT2(3,3)=-1.0_R8_KIND

    if(ZDIST>=1.0D-12) then

       ZSIGN = 1.0_R8_KIND
       if(DIF(3).LT.0) ZSIGN = -1.0_R8_KIND
       ODIST = ZDIST + DIF(3)*DIF(3)
       ODIST = sqrt(ODIST)
       ZDIST = sqrt(ZDIST)
       ROTMAT2(1,1) = ZSIGN*DIF(1)*DIF(3)/ZDIST/ODIST
       ROTMAT2(2,1) = ZSIGN*DIF(2)*DIF(3)/ZDIST/ODIST
       ROTMAT2(3,1) = -ZSIGN*ZDIST/ODIST
       ROTMAT2(1,2) = -DIF(2)/ZDIST
       ROTMAT2(2,2) = DIF(1)/ZDIST
       ROTMAT2(3,2) = 0.0_R8_KIND
       ROTMAT2(1,3) = DIF(1)/ODIST
       ROTMAT2(2,3) = DIF(2)/ODIST
       ROTMAT2(3,3) = DIF(3)/ODIST
    endif

  end function ROTMAT2
  !---------------------------------------------------------------------



  subroutine CALC_CHARGE()
    !---------------------------------------------------------------------
    ! Purpose:  SUBROUTINE CALCULATES THE CENTER OF NUCLEAR CHARGE
    !  RESULT IS STORED IN CHARGE_CENTER
    !---------------------------------------------------------------------
    !** End of interface *****************************************
    implicit none

    real(kind=R8_KIND)    :: CHARGE

    integer  :: I,J


    CHARGE_CENTER=0.0_R8_KIND
    CHARGE=0.0_R8_KIND

    do I=1,N_UNIQUE_ATOMS
       do J=1,UNIQUE_ATOMS(I)%N_EQUAL_ATOMS
          CHARGE_CENTER=CHARGE_CENTER+UNIQUE_ATOMS(I)%POSITION(:,J)*&
               UNIQUE_ATOMS(I)%Z
          CHARGE=CHARGE+UNIQUE_ATOMS(I)%Z
       enddo
    enddo

    CHARGE_CENTER=CHARGE_CENTER/CHARGE
    !---------------------------------------------------------------------
  end subroutine calc_charge

  !*************************************************************
  subroutine grid_loop_setup()
    !
    ! Setup of  the dlb- routines,  which will be used  to dynamically
    ! (re)-distribute the jobs, the dlb needs for a start only to know
    ! the job distribution  as here all jobs are  equal, there is only
    ! need of the total number of them
    !
    use dlb, only: dlb_setup, idlb_kind
    implicit none
    ! *** end of interface ***

    integer(idlb_kind) :: number_of_batches

    ! transform it in the integer kind used by dlb
    number_of_batches = size(grdpts, 1)
    call dlb_setup(number_of_batches)
  end subroutine grid_loop_setup

  subroutine grid_loop_setup_atom()
    !
    ! Setup of  the dlb- routines,  which will be used  to dynamically
    ! (re)-distribute the jobs, the dlb needs for a start only to know
    ! the job distribution
    !
    use dlb, only: dlb_setup_color, idlb_kind
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: i
    integer(idlb_kind) :: atomic_grid_sizes(size(agrid)) ! (NUA)

    do i = 1, size(agrid)
        atomic_grid_sizes(i) = size(agrid(i)%m, 1)
    enddo

    call dlb_setup_color(atomic_grid_sizes)
  end subroutine grid_loop_setup_atom


  subroutine mkgrid(agrid, grdpts)
    !
    ! Compute Becke  weights. In  SCF there is  no point  to recompute
    ! them every iteration as the geometry does not change. Also fills
    ! linear array of grid point positions.
    !
    ! NOTE: requires allocatable arguments capability
    !
    use comm, only: comm_allreduce
    use dlb, only: dlb_setup_color, dlb_give_more_color, long => idlb_kind
    use machineparameters_module, only: maxpts => machineparameters_veclen
    implicit none
    type(arrmat2), intent(inout) :: agrid(:)     ! (n_unique_atoms)
    real(r8_kind), allocatable   :: grdpts(:, :) ! (npts, 4) on output
    ! *** end of interface ***

    integer (long) :: gcounts(size(agrid)) ! (n_unique_atoms)
    integer (long) :: offsets(size(agrid)) ! (n_unique_atoms)
    integer (long) :: npts
    integer (long) :: slice(2), maxpts_dlb_kind

    ! "atomic coordinates" of grid slice
    integer (long) :: ua_dlb
    integer (long) :: lo, hi
    integer (i4_kind) :: ua     ! == ua_dlb

    ! "molecular coordinates" of grid slice
    integer (long) :: glo, ghi

    FPP_TIMER_START(mkgrd)

    ASSERT(.not.allocated(grdpts))

    !
    ! Number of grid points for each atom and offsets into
    ! a linear array:
    !
    do ua = 1, size(agrid)
        gcounts(ua) = size(agrid(ua)%m, 1)
        offsets(ua) = sum(gcounts(1:ua-1))
    end do

    ! total number of points:
    npts = sum(gcounts)

    !
    ! Allocate output array:
    !
    allocate(grdpts(npts, 4))

    !
    ! We will abuse comm_allreduce(grdwts) to collect computed
    ! weights, zero is idepmotent wrt addition:
    !
    grdpts(:, 4) = 0.0

    call dlb_setup_color(gcounts)
    maxpts_dlb_kind = maxpts
    ! FIXME: one might need to adjust maxpts here:
    do while ( dlb_give_more_color(maxpts_dlb_kind, ua_dlb, slice) )
        !
        ! DLB returns half open intervals
        !
        !       (slice(1), slice(2)]
        !
        ! thus granting all "k" such that
        !
        !       slice(1) < k <= slice(2)
        !
        ua = int (ua_dlb, kind = i4_kind) ! passed to atomicweight()
        lo = slice(1) + 1
        hi = slice(2)
        ASSERT(hi-lo+1<=maxpts)

        !
        ! I was granted this section of AGRID to process:
        !
        ! x => agrid(ua)%m(lo:hi, 1:3)
        ! w => agrid(ua)%m(lo:hi, 4)

        !
        ! This corresponds to this section in a linear array grdpts:
        !
        glo = offsets(ua) + lo
        ghi = offsets(ua) + hi

        !
        ! Atomic quadrature weights multiplied by
        ! (Becke) space partitioning function:
        !
        grdpts(glo:ghi, 4) =  agrid(ua)%m(lo:hi, 4) &
            * atomicweight(ua, agrid(ua)%m(lo:hi, 1:3)) &
            * unique_atoms(ua)%n_equal_atoms
        ! In case it was internal updated
        maxpts_dlb_kind = maxpts
    enddo

    !
    ! Collect all computed weights and bcast result:
    !
    call comm_allreduce(grdpts(:, 4))

    !
    ! Update local structure AGRID with weights scaled by
    ! (Becke) partitioning function, copy positions into
    ! linear array:
    !
    do ua = 1, size(agrid)
        lo = offsets(ua) + 1
        hi = offsets(ua) + gcounts(ua)

        ! update weights:
        agrid(ua)%m(:, 4) = grdpts(lo:hi, 4)

        ! copy positions:
        grdpts(lo:hi, 1:3) = agrid(ua)%m(:, 1:3)
    end do

    FPP_TIMER_STOP(mkgrd)
  end subroutine mkgrid

  !-------------------------------------------------------------------------
  subroutine grid_close1(keep_grid_par_ph)
    ! Purpose: wrapper around grid_close for using it in a master-only context
    use comm_module, only: comm_init_send, comm_all_other_hosts, comm_send, &
         comm_i_am_master
    use msgtag_module, only: msgtag_gridph_close, msgtag_grid_close
    implicit none
    logical, intent(in) :: keep_grid_par_ph
    !** End of interface *****************************************

    if (comm_i_am_master()) then
      if (keep_grid_par_ph) then
        call comm_init_send(comm_all_other_hosts, msgtag_gridph_close)
      else
        call comm_init_send(comm_all_other_hosts, msgtag_grid_close)
      endif
      call comm_send()
    endif

    call grid_close(keep_grid_par_ph)
  end subroutine grid_close1

  subroutine grid_close(keep_grid_par_ph)
    !
    ! Purpose:  deallocate  variables,  cleans   up  the  module  state.  This
    ! subroutine is called from
    !
    !     main_scf
    !     main_slave
    !     post_scf_module
    !     response_module
    !     modules/initialization.f90
    !
    ! Idempotent.
    !
    ! The difference between the calls is only, to know if the grid_par_ph has
    ! to  be  kept  for  further  use  this  is  decided  by  keep_grid_par_ph
    ! everything else which is there  is deallocated. FIXME: all other modules
    ! clean everything unconditionally.
    !
    implicit none
    logical, intent(in) :: keep_grid_par_ph
    !** End of interface *****************************************

    integer(kind=i4_kind) :: alloc_stat, i

    !
    ! FIXME: use of alloc_grid(4) not quite transparent.
    !
    if (allocated(agrid)) then
      do i = 1, size(agrid)
        deallocate(agrid(i)%m, stat=alloc_grid(4))
        ASSERT(alloc_grid(4).eq.0)
      end do
      deallocate(agrid, stat=alloc_grid(3))
      ASSERT(alloc_grid(3).eq.0)
    endif

    if (allocated(grdpts)) then
      deallocate(grdpts, stat=alloc_stat)
      ASSERT (alloc_stat==0)
    endif

    if (allocated(grid_par_scf)) then
      deallocate(grid_par_scf, stat=alloc_stat)
      ASSERT (alloc_stat==0)
    endif

    if (allocated(grid_par_ph) .and. .not. keep_grid_par_ph) then
      deallocate(grid_par_ph, stat=alloc_grid(2))
      ASSERT(alloc_grid(2).eq.0)
    endif

#ifdef FPP_TIMERS
    print *,'Timings:'
    print *,'grid_main = ',FPP_TIMER_VALUE(gmain)
    print *,'mkgrid    = ',FPP_TIMER_VALUE(mkgrd)
    print *,'AWT2      = ',FPP_TIMER_VALUE(awt)
    print *,' |-0      = ',FPP_TIMER_VALUE(awt0)
    print *,' |-1      = ',FPP_TIMER_VALUE(awt1)
    print *,' `-2      = ',FPP_TIMER_VALUE(awt2)
#endif
  end subroutine grid_close

  subroutine find_circular_set(uatom,RO,roz)
    !------------------------------------------------------------------------
    implicit none
    integer,intent(IN)  :: uatom
    real(kind=R8_KIND),intent(IN),dimension(3,3) :: RO
    real(kind=R8_KIND),intent(out),dimension(3,3) :: ROZ
    !** End of interface *****************************************
    integer      :: I,K,ALLOCSTAT,z_max_nunique,z_max_nequal

    type HELP_COORD
       real(kind=R8_KIND), dimension(:,:), allocatable :: POSITION
    end type HELP_COORD

    type(HELP_COORD),dimension(N_UNIQUE_ATOMS) :: H_COORD
    logical,save   :: LINEAR=.false.
    logical :: first
    real(kind=r8_kind) :: cospi,sinpi,z_dist,z_dist_max
    real(kind=r8_kind),parameter :: one=1.0_r8_kind,zero=0.0_r8_kind

    if (output_grid) write(output_unit,*) 'Find circular set entered'

    roz=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))

    if (linear) then
       return
    endif
    linear=.true.
    first=.true.
    z_dist_max=0.0_r8_kind


    do I=1,N_UNIQUE_ATOMS
       allocate(H_COORD(I)%POSITION(3,UNIQUE_ATOMS(I)%N_EQUAL_ATOMS), &
            STAT=ALLOCSTAT)
       if(ALLOCSTAT/=0) call error_handler&
            ('ALLOCATION FAILED IN SU GET_SYMMETRY')
       do K=1,UNIQUE_ATOMS(I)%N_EQUAL_ATOMS
          H_COORD(I)%POSITION(:,K)=sum(RO*spread(UNIQUE_ATOMS(I)%POSITION &
               (:,K)-UNIQUE_ATOMS(uatom)%POSITION(:,1),2,3),1)
          z_dist=H_COORD(I)%POSITION(1,K)**2+H_COORD(I)%&
               POSITION(2,K)**2    ! distance from z-axis
          linear=linear.and.(z_dist<=1.0D-16)
          if(.not.(i==uatom.and.k==1)) then
             if((z_dist-z_dist_max)>=1.0D-16.or.first) then
                z_dist_max=z_dist
                z_max_nunique=i
                z_max_nequal=k
                first=.false.
             else
                if(abs(z_dist-z_dist_max)<=1.0D-16) then
                   if(i/=z_max_nunique) then
                      if((abs(h_coord(i)%position(3,k))-abs(h_coord(&
                           z_max_nunique)%position(3,z_max_nequal)))>=1.0D-10)&
                           then
                         z_dist_max=z_dist
                         z_max_nunique=i
                         z_max_nequal=k
                      else
                         if(h_coord(i)%position(3,k)>0.and.h_coord&
                              (z_max_nunique)%position(3,z_max_nequal)<0) then
                            z_dist_max=z_dist
                            z_max_nunique=i
                            z_max_nequal=k
                         else
                            if(unique_atoms(i)%z<unique_atoms&
                                 (z_max_nunique)%z) then
                               z_dist_max=z_dist
                               z_max_nunique=i
                               z_max_nequal=k
                            else
                               if (output_grid) write(output_unit,*) &
                                    'CAUTION: ORIENTATION OF THE GRID',&
                                    ' DEPENDS OF THE ORDER OF THE INPUT'
                            endif
                         endif
                      endif
                   endif
                endif
             endif
          endif
       enddo
    enddo

    if (LINEAR) then
       if (output_grid) write(output_unit,*) 'Linear Molecule'
       return
    endif

    cospi=h_coord(z_max_nunique)%position(1,z_max_nequal)/sqrt(z_dist_max)
    sinpi=h_coord(z_max_nunique)%position(2,z_max_nequal)/sqrt(z_dist_max)

    roz(1,1)=cospi
    roz(2,2)=cospi
    roz(1,2)=-sinpi
    roz(2,1)=sinpi

    if (output_grid) write(output_unit,*) 'ROZ:',cospi,sinpi

    do I = 1, size(H_COORD)
       deallocate(H_COORD(I)%POSITION, stat=ALLOCSTAT)
       ASSERT(ALLOCSTAT==0)
    enddo
  end subroutine find_circular_set

  !-------------------------------------------------------------------------

  subroutine get_symmetry(uatom, ro)
    !
    ! This  routine detects the  local symmetry  of the  unique uatom.
    ! The local symmetry is stored in the structure sym_element_local.
    ! This routine is called by  inigrd, which uses the information to
    ! reduce the size of the  integration grid. It also determines the
    ! orientation  of the  grid, which  is  passed to  inigrd via  the
    ! rotation matrix ro. The orientation is as follows. The z-axis is
    ! pointing from  the unique  to the center  of charge.  The second
    ! direction is  perpendicular to  the first and  lying in  a local
    ! symmetry-plane. If no local  symmetry-plane is present, then the
    ! second direction  is towards an  other atom, which is  chosen in
    ! the routine get_circular_set.
    !
    use symmetry_element, only: sym_element, sym_element_local
    ! modifies sym_element_local
    implicit none
    integer(kind=i4_kind)   :: uatom
    real(kind=r8_kind),intent(out),dimension(3,3) :: ro
    !** End of interface *****************************************

    real(kind=r8_kind),dimension(3,3)  :: roz
    real(kind=r8_kind) :: xnorm,skalar,sqrt3m1,sqrt2,xnorm_rev
    real(kind=r8_kind),dimension(3) :: x,x_axis,x_axis_trans,z_axis
    integer(kind=i4_kind)  ::  n_axis_max,i_n_axis,i,counter,n_help

    real(kind=r8_kind),parameter :: one=1.0_r8_kind,zero=0.0_r8_kind
    real(kind=r8_kind),parameter :: small=1.0e-14_r8_kind
#ifndef FPP_AIX_XLF
    intrinsic matmul
#endif
    real(kind=r8_kind),dimension(3,3) :: mat

    ASSERT(sym_reduce)

    counter=0
    if (output_grid) write(output_unit,*) 'Get Symmetry entered'
    x=unique_atoms(uatom)%position(:,1)
    xnorm=sqrt(sum(x*x))
    n_axis_max=0

    skalar=zero
    n_help=0
    if(associated(sym_element%rotaxis)) then
      n_help= size(sym_element%rotaxis,2)
    endif

    !needed because of agressive optimization on VPP
    if(xnorm>small) xnorm_rev=one/xnorm

    if ( output_grid) write(output_unit,*) 'nhelp',n_help

!   do i=1,size(sym_element%rotaxis,2)
    do i=1,n_help
!      if(xnorm/=0.0_r8_kind) then
       if(xnorm>small) then
!         skalar=sum(sym_element%rotaxis(:,i)*x)/xnorm
          skalar=sum(sym_element%rotaxis(:,i)*x)*xnorm_rev
       endif
!      if(abs(abs(skalar)-1.0_r8_kind)<=1.0D-12.or.xnorm==0.0_r8_kind) then
       if(abs(abs(skalar)-1.0_r8_kind)<=small.or.xnorm<=small) then
          counter=counter+1
          sym_element_local%rotaxis(:,counter)=sym_element%rotaxis(:,i)
          sym_element_local%n_axis(counter)=sym_element%n_axis(i)

          if(sym_element_local%n_axis(counter)>n_axis_max) then
             i_n_axis = counter
             n_axis_max=sym_element_local%n_axis(counter)
          endif
       endif
    enddo
    sym_element_local%num_axis=counter

    if(n_axis_max>0) then
       z_axis=sym_element_local%rotaxis(:,i_n_axis)
    endif

    if (output_grid) then
       write(output_unit,*) 'Number of local rotaxes',sym_element_local%num_axis
    endif

    counter=0
    if( associated(sym_element%sigma) )then
      ! FIXME: how does it happen that %sigma is not associated?
      do i=1,size(sym_element%sigma,2)
         skalar=sum(sym_element%sigma(:,i)*x)
         if(abs(skalar)<=1.0D-12) then
            counter=counter+1
            sym_element_local%sigma(:,counter)=sym_element%sigma(:,i)
         endif
      enddo
    endif ! associated(sym_element%sigma)
    sym_element_local%num_sigma=counter
    if (output_grid)then
       write(output_unit,*) 'Number of local mirror planes', sym_element_local%num_sigma
    endif

    !
    ! Now calculate rotation matrix compatible with this site symmetry:
    !

    !
    ! rotate oktaeder if rotationaxis is threefold
    !
    if(n_axis_max==3) then
       ro=rotmat2(z_axis)
       if(n_unique_atoms/=1) then
          call find_circular_set(uatom,ro,roz)
       else
          roz=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))
       end if
       sqrt3m1=1.0_r8_kind/sqrt(3.0_r8_kind)
       sqrt2=sqrt(2.0_r8_kind)
       z_axis=(/one,one,one/)
       mat = reshape((/sqrt2*sqrt3m1,0.0_r8_kind,sqrt3m1,-sqrt3m1/sqrt2,&
            1/sqrt2,sqrt3m1,-sqrt3m1/sqrt2,-1/sqrt2,sqrt3m1/),(/3,3/))
       ro=MATMUL(ro,mat)
       ro=MATMUL(ro,roz)
       return
    endif

    if(n_axis_max/=0.and.counter/=0) then
       skalar=sum(sym_element_local%rotaxis(:,i_n_axis)*sym_element_local%&
            sigma(:,1))
       if(abs(skalar-1)>1.0E-12_r8_kind) then
          x_axis=sym_element_local%sigma(:,1)
          ro=rotmat2(z_axis)
          x_axis_trans=MATMUL(x_axis,ro)
          roz=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))
          roz(1,1)=x_axis_trans(1)
          roz(2,2)=x_axis_trans(1)
          roz(1,2)=x_axis_trans(2)
          roz(2,1)=-x_axis_trans(2)
       else
          ro=rotmat(uatom)
          if(n_unique_atoms/=1) then
             call find_circular_set(uatom,ro,roz)
          else
             roz=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))
          end if
       endif
    else
       ro=rotmat(uatom)
       if(n_unique_atoms/=1) then
          call find_circular_set(uatom,ro,roz)
       else
          roz=reshape((/one,zero,zero,zero,one,zero,zero,zero,one/),(/3,3/))
       end if
    endif
    ro=MATMUL(ro,roz)
  end subroutine get_symmetry

  subroutine print_alloc_grid()
  integer(kind=i4_kind):: i
  do i=1,size(alloc_grid)
  if(alloc_grid(i).ne.0) print*,'alloc_grid allocated',i
  enddo
  end subroutine print_alloc_grid


#ifdef WITH_GUILE
  function guile_pople_radius (z) result (r) bind (c)
    !
    ! Export rpople_f() to Scheme.
    !
    use scm, only: scm_t, scm_to_int, assignment(=)
    implicit none
    type(scm_t), intent(in), value :: z ! SCM int
    type(scm_t) :: r                    ! SCM double
    ! *** end of interface ***

    r = rpople_f (scm_to_int (z))
  end function guile_pople_radius

  function guile_slater_radius (z) result (r) bind (c)
    !
    ! Export rslater_f() to Scheme.
    !
    use scm, only: scm_t, scm_to_int, assignment(=)
    implicit none
    type(scm_t), intent(in), value :: z ! SCM int
    type(scm_t) :: r                    ! SCM double
    ! *** end of interface ***

    r = rslater_f (scm_to_int (z))
  end function guile_slater_radius

  function guile_ionic_radius (z) result (r) bind (c)
    !
    ! Export rionic_f() to Scheme.
    !
    use scm, only: scm_t, scm_to_int, assignment(=)
    implicit none
    type(scm_t), intent(in), value :: z ! SCM int
    type(scm_t) :: r                    ! SCM double
    ! *** end of interface ***

    r = rionic_f (scm_to_int (z))
  end function guile_ionic_radius
#endif

end module grid_module


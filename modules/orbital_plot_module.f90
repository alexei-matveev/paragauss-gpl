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
!=====================================================================
! Public interface of module
!=====================================================================
module orbital_plot_module
  !-------------------------------------------------------------------
  !
  !  Contains subroutines for plotting orbitals on a special grid. The
  !  orbitals, specified in the input,  evaluated on the plot grid and
  !  then written to a file which subsequently has to be passed to the
  !  final  plotting  programm.   With  the option  split_mode  it  is
  !  possible  to print the  contributions of  every atom  and angular
  !  momentum seperately.  It is  also possible to calculate (valence)
  !  densities and spin density differences.
  !
  !  Module called by: properties_module
  !
  !  References:
  !
  !
  !  Author: MS
  !  Date: 12/97
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: MH
  ! Date:   08/2007
  ! Description: Added NTO plots. The modifications can be seen
  !              searching for "MODIFICATION FOR NTO"
  !-------------------------------------------------------------------
# include "def.h"
  use type_module, only: i4_kind, r8_kind
  use symmetry_data_module
  use orbitalstore_module, only: orbital_type
  use unique_atom_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ public functions and subroutines ---------------------
  public :: orbital_plot_main
  public :: orbital_plot_read
  public :: orbital_plot_write
  public :: orbital_plot_bcast
  public :: orbital_plot_check


  !===================================================================
  ! End of public interface of module
  !===================================================================


  type :: inp ! USED FOR INPUT, UGLY!
     ! type to store which orbitals have to be plotted
     integer(i4_kind) :: n_orbs ! number of orbitals
     integer(i4_kind), allocatable :: list(:) ! indices of orbitals
     integer(i4_kind) :: irrep
     character(len=3) :: irrep_name
     integer(i4_kind) :: partner
     integer(i4_kind) :: spin
     integer(i4_kind) :: trans ! MODIFICATION FOR NTO
     logical :: range          ! if true not a complete list, but only
                               ! lowest and highet orbital have to be
                               ! given in the input
     logical :: split_mode     ! if true, the wavefunctions will be
                               ! split into angular and atomic
                               ! contributions for the atoms given in
                               ! split_list
     integer(i4_kind) :: n_split_atoms ! dimension of split_list
     integer(i4_kind), allocatable :: split_list(:)
     integer(i4_kind) :: n_orbs_split ! number of single contributions
                                      ! to one orbital
  end type inp

  type orb ! LINEAR IN FINAL GRID-ORBITALS
     integer(i4_kind) :: principal
     integer(i4_kind) :: irrep
     integer(i4_kind) :: partner
     integer(i4_kind) :: spin
     character(len=32) :: name
  end type orb

  type(orbital_type) :: fcts_p
  integer(i4_kind) :: th = -1
  integer(i4_kind) :: num_nto
  integer(i4_kind), allocatable :: size_nto(:)
  type(inp), allocatable :: orb_list(:)
  type(orb), allocatable :: list_orb(:) ! "reverse" mapping

  ! defaults for input
  real(r8_kind), parameter :: DD = 3.0
  real(r8_kind), parameter :: df_x0(3) = [-DD, -DD, -DD]
  real(r8_kind), parameter :: df_x1(3) = [+DD, -DD, -DD]
  real(r8_kind), parameter :: df_x2(3) = [-DD, +DD, -DD]
  real(r8_kind), parameter :: df_x3(3) = [-DD, -DD, +DD]

  integer(i4_kind) :: RES(3) ! res in X, Y, and Z
  integer(i4_kind), parameter :: df_RES(3) = [20, 20, 20]

  integer(i4_kind), parameter :: &
       df_irrep = 0, &
       df_partner = 1, &
       df_spin = 0,  &
       df_n_orb = 0, &
       df_n_input_lines = 0

  integer(i4_kind), parameter :: df_n_split_atoms = -1 ! FIXME: was
                                                       ! not given a
                                                       ! value
  logical, parameter :: &
       df_range           = .false., &
       df_split_mode      = .false., &
       df_density_plot    = .false., &
       df_density_tot     = .false., &
       df_spin_difference = .false.

  logical, parameter :: &
       df_alternative_options = .false., &
       df_calc_rho            = .false., &
       df_calc_rho_fit        = .false., &
       df_calc_xc_pot         = .false., &
       df_calc_xc_pot_mda     = .false., &
       df_nto                 = .false. ! MODIFICATION FOR NTO

  character(len=3), parameter :: df_irrep_name = '   '
  character(len=2), parameter :: df_trans = 'SS'

  real(r8_kind) :: x0(3), &     ! origin
                   x1(3), &     ! corner 1
                   x2(3), &     ! corner 2
                   x3(3)        ! corner 3

  integer(i4_kind) :: n_input_lines = 0
  integer(i4_kind) :: n_orb     ! number of orbitals

  logical :: density_plot, &  ! if true not single orbitals are
                              ! plotted, but the density defined by
                              ! all (!) orbitals defined in the input
                              ! is plotted
             spin_difference,&! if true not the sum over the two spin
                              ! componets, but the difference is
                              ! evaluated
             density_tot      ! if true the whole density is plotted,
                              ! makes only sense if density_plot is
                              ! true

  logical :: alternative_options,                    &
             calc_rho,        &  ! altern. option to draw density
             calc_rho_fit,    &  ! altern. option to draw fitted density
             calc_xc_pot,     &  ! altern. option to draw V_xc[rho_exact]
             calc_xc_pot_mda, &  ! altern. option to draw V_xc[fitted_rho]
             nto                 ! MODIFICATION FOR NTO

  !
  ! These are not  part of the input, but  historically kept as global
  ! vars. FIXME: avoid them:
  !
  integer(i4_kind) :: grid_length ! length of the plotting grid

  type(orbital_type), pointer :: orbs_ob(:) ! see orbital_allocate()
  real(r8_kind), allocatable :: phi(:,:)
  ! phi(grid_length,n_orbs) value of orbitals on grid data structure
  ! to keep prim orbitals
  real(r8_kind), allocatable :: rho(:) ! rho(grid_length): density
  real(r8_kind), allocatable, target :: rho_s(:,:), rho_fit(:,:), &
       xc_pot(:,:), xc_pot_mda(:,:), fxc(:)
  real(r8_kind), allocatable, target :: gamm(:,:), xcpotgr(:,:)

contains

  !*************************************************************
  subroutine orbital_plot_read()
    ! Purpose: read input connected with plotting
    !          orbitals
    use input_module
    implicit none
    !** End of interface *****************************************

    integer :: unit, status
    integer(i4_kind) :: i_input, i_low, i_high, alloc_stat, i, &
         i_split

    integer(i4_kind) :: &
         irrep, &   ! irrep
         partner, & ! partner
         spin,  &   ! spin
         n_split_atoms ! number of atoms to be included in splitted
                       ! orbital plot

    logical :: range ! if true not a complete list of orbital number
                     ! must be given, but only lower and upper limit

    logical :: split_mode ! if true the orbitals are split in to
                          ! atomic and angular moment contributions

    character(len=3) :: irrep_name
    character(len=2) :: trans

    namelist /orbital_plot/ n_input_lines, &
         density_plot, &
         spin_difference, &
         density_tot, &
         nto,      &
         calc_rho, &
         calc_rho_fit, &
         calc_xc_pot,  &
         calc_xc_pot_mda, &
         X0, X1, X2, X3, &
         RES

    namelist /orbital_list/ irrep, irrep_name, partner, spin, n_orb, &
         trans, range, split_mode, n_split_atoms


    ! read namelist orbital_plot
    n_input_lines=df_n_input_lines
    density_plot=df_density_plot
    spin_difference=df_spin_difference
    alternative_options = df_alternative_options
    calc_rho            = df_calc_rho
    calc_rho_fit        = df_calc_rho_fit
    calc_xc_pot         = df_calc_xc_pot
    calc_xc_pot_mda     = df_calc_xc_pot_mda
    nto                 = df_nto               !MODIFICATION FOR NTO
    trans               = df_trans
    x0=df_x0
    x1=df_x1
    x2=df_x2
    x3=df_x3
    RES = df_RES

    if ( input_line_is_namelist("orbital_plot") ) then
       call input_read_to_intermediate()
       unit = input_intermediate_unit()
       read(unit, nml=orbital_plot, iostat=status)
       if( status /= 0 )then
          call error_handler("orbital_plot_read : error reading namelist ORBITAL_PLOT")
       endif
       ASSERT(status==0)

       ! write (*, nml=orbital_plot)

       alternative_options = calc_rho    .or. calc_rho_fit .or. &
                             calc_xc_pot .or. calc_xc_pot_mda
       if(alternative_options) then
          if (density_tot) &
               call error_handler &
               ("orbital_plot_read : DENSITY_TOT and CALC_(...) are the same")
          density_plot = alternative_options
          density_tot  = alternative_options
          if(n_input_lines /= 0) call error_handler(&
               "orbital_plot_read : Inconsistent input - N_INPUT_LINES must be zero")
          if(calc_xc_pot)     calc_rho     = alternative_options
       end if
    end if
    if(n_input_lines>0) then
       allocate(orb_list(n_input_lines), stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if
    if(spin_difference.and..not.density_plot) &
         call input_error(&
         'orbital_plot_read: sorry, you must set density_plot to true, if you want to calculate spin_difference')
    if (density_tot .and. .not. density_plot) &
        call input_error &
        ('orbital_plot_read: sorry, you must set density_plot to true, if you set density_tot to true')
    if (density_tot .and. n_input_lines /= 0) &
         call input_error &
         ('orbital_plot_read: sorry, n_input_lines together with density_tot does not make sense')
    do i_input=1,n_input_lines
       irrep=df_irrep
       irrep_name = df_irrep_name
       partner=df_partner
       spin=df_spin
       range=df_range
       n_orb=df_n_orb
       split_mode=df_split_mode
       n_split_atoms=df_n_split_atoms
       if ( input_line_is_namelist("orbital_list") ) then
          call input_read_to_intermediate()
          unit = input_intermediate_unit()
          read(unit, nml=orbital_list, iostat=status)
          if (status .gt. 0) call input_error( &
               "orbital_plot_read: namelist orbital_plot_read")
          orb_list(i_input)%irrep_name=irrep_name
          !----------------------!
          ! MODIFICATION FOR NTO !
          !----------------------!
          orb_list(i_input)%trans=1
          if(trans.eq.'ST') orb_list(i_input)%trans=2
          !--------------------------!
          ! END MODIFICATION FOR NTO !
          !--------------------------!
          if(range) then
             ! number of  orbitals is given by a  range but internally
             ! stored as a usuall list
             orb_list(i_input)%range=.true.
             call input_read_to_intermediate()
             read(unit,*) i_low,i_high
             if(i_low<1.or.i_high<i_low) call input_error(&
                  'Wrong values for range')
             orb_list(i_input)%n_orbs=1+i_high-i_low
             if(irrep<1) call input_error( &
                  'orbital_plot_read:  wrong value for irrep')
             orb_list(i_input)%irrep=irrep
             if(partner<1) call input_error( &
                  'orbital_plot_read:  wrong value for partner')
             if(partner/=1.and.density_plot) call input_error( &
                  'orbital_plot_read: sorry, density_plot and specifiying a partner does not make sense')
             orb_list(i_input)%partner=partner
             if(spin<1.or.spin>2) call input_error( &
                  'orbital_plot_read:  wrong value for spin')
             orb_list(i_input)%spin=spin
             allocate(orb_list(i_input)%list(1+i_high-i_low), stat=&
                  alloc_stat)
             ASSERT(alloc_stat==0)
             do i=1,i_high-i_low+1
                orb_list(i_input)%list(i)=i-1+i_low
             end do
          else
             if(n_orb<=0) call input_error(&
                  'orbital_plot_read: non positive value for n_orb')
             orb_list(i_input)%range=.false.
             orb_list(i_input)%irrep=irrep
             orb_list(i_input)%partner=partner
             orb_list(i_input)%spin=spin
             orb_list(i_input)%n_orbs=n_orb
             allocate(orb_list(i_input)%list(n_orb), stat=alloc_stat)
             ASSERT(alloc_stat==0)
             call input_read_to_intermediate()
             read(unit,*)  orb_list(i_input)%list
             if(minval(orb_list(i_input)%list)<1) call error_handler(&
                  'orbital_plot_read: sorry your list contains negative  values')
          end if
          orb_list(i_input)%split_mode=.false.
          orb_list(i_input)%n_split_atoms=0
          if(split_mode) then
             if(density_plot) call input_error(&
                  'Sorry, you can not use split_mode for density_plot')
             if(n_split_atoms<=0) call input_error(&
                  'orbital_plot_read: n_split_atoms must be positive')
             orb_list(i_input)%n_split_atoms=n_split_atoms
             orb_list(i_input)%split_mode=split_mode
             allocate(orb_list(i_input)%split_list(n_split_atoms), stat=alloc_stat)
             ASSERT(alloc_stat==0)
             call input_read_to_intermediate()
             read(unit,*) orb_list(i_input)%split_list
             do i_split=1,n_split_atoms
                if(orb_list(i_input)%split_list(i_split)<=0.or.&
                     orb_list(i_input)%split_list(i_split)> n_unique_atoms) &
                     call input_error('split_list has wrong value')
             end do
          endif
       else
          call input_error('Sorry, at least one additional namelist orbital_list is required')
       end if
    end do ! i_input
    !----------------------!
    ! MODIFICATION FOR NTO !
    !----------------------!
#ifdef WITH_RESPONSE
    if(nto) call parse_nto_to_mo()
#endif
    !---------------------------!
    ! END MODIFICATION FOR NTO  !
    !---------------------------!
  end subroutine orbital_plot_read

  !----------------------!
  ! MODIFICATION FOR NTO !
  !----------------------!
#ifdef WITH_RESPONSE
  subroutine parse_nto_to_mo()
      ! Purpose: Expand the NTOs into MOs
      use nto_plot_module
      implicit none
      !** End of interface *****************************************

      integer(i4_kind)                :: i_input

      !** Calculate the nto expansion **
      num_nto=size(orb_list)

      !** Determine which MO expand all the NTOs selected
      do i_input=1,num_nto
         call nto_plot_search( orb_list(i_input)%irrep , &  ! Irrep number
                               orb_list(i_input)%spin  , &  ! Spin number
                               orb_list(i_input)%list  , &  ! NTO list
                               i_input                 , &  ! NTO index
                               num_nto                 , &  ! Total number of NTOs
                               orb_list(i_input)%trans    ) ! SS[1] or ST[2]?
      end do
!
      call nto_plot_order()      ! Order the MOs according to 1) Spin,
                                 ! 2) Symmetry, 3) Index

!
      call nto_plot_merge()      ! Put the MO values into orb_list,
                                 ! matrix that the rest of this module
                                 ! use to calculate the grid points
!

      !****************************************************************
      CONTAINS
      SUBROUTINE nto_plot_merge()
          !  Purpose: Expand the NTOs into MOs
          IMPLICIT NONE
          !** End of interface *****************************************

          INTEGER(i4_kind)                  :: idx,jdx
          INTEGER(i4_kind)                  :: up_lim,low_lim
          INTEGER(i4_kind)                  :: current_irrep,break_value,idx_break
          INTEGER(i4_kind)                  :: control_value, alloc_stat
          !------------ Executable code --------------------------------

          ! Recursive  deallocation   of  structure  with  allocatable
          ! components:
          deallocate(orb_list, stat=alloc_stat)
          ASSERT(alloc_stat==0)

          ! Set the new number of input lines and where the first spin
          ! and the second are located

          break_value=1
          idx_break=1
          control_value=nto_orb(1,2)
          do i_input=2,size(nto_orb,1)
              if(nto_orb(i_input,3).ne.nto_orb(1,3)) exit
              break_value=break_value+1
              if(nto_orb(i_input,2).ne.control_value) then
                  n_input_lines=n_input_lines+1
                  control_value=nto_orb(i_input,2)
                  if(nto_orb(i_input,3).eq.nto_orb(1,3)) idx_break=idx_break+1
              end if
          end do

          if(any(nto_orb(:,3).eq.1).and.any(nto_orb(:,3).eq.2)) then
              control_value=nto_orb(break_value+1,2)
              do i_input=break_value+2,size(nto_orb,1)
                  if(nto_orb(i_input,2).ne.control_value) then
                      n_input_lines=n_input_lines+1
                      control_value=nto_orb(i_input,2)
                      if(nto_orb(i_input,3).eq.nto_orb(1,3)) idx_break=idx_break+1
                  end if
              end do
          end if
!
          !** Allocate the new orb_list
          allocate(orb_list(n_input_lines))
!
          low_lim=1
          up_lim=1
          current_irrep=nto_orb(1,2)
          do idx=1,idx_break
             orb_list(idx)%irrep=current_irrep
             orb_list(idx)%spin = 1
             orb_list(idx)%partner = 1
             orb_list(idx)%range=.false.
             orb_list(idx)%split_mode=.false.
             do jdx=1,break_value
                 if(nto_orb(jdx,2).lt.current_irrep) cycle
                 if(current_irrep.ne.nto_orb(jdx,2)) then
                     current_irrep = nto_orb(jdx,2)
                     exit
                 end if
                 up_lim=up_lim+1
             end do
             allocate(orb_list(idx)%list(up_lim-low_lim))
             orb_list(idx)%n_orbs=up_lim-low_lim
             orb_list(idx)%list(:)=nto_orb(low_lim:up_lim-1,1)
             low_lim=up_lim
          end do
!
          if(any(nto_orb(:,3).eq.1).and.any(nto_orb(:,3).eq.2)) then
              current_irrep=nto_orb(low_lim,2)
              do idx=idx_break+1,size(orb_list)
                  orb_list(idx)%irrep=current_irrep
                  orb_list(idx)%spin = 2
                  orb_list(idx)%partner = 1
                  orb_list(idx)%range=.false.
                  orb_list(idx)%split_mode=.false.
                  do jdx=break_value+1,size(nto_orb,1)
                      if(nto_orb(jdx,2).lt.current_irrep) cycle
                      if(current_irrep.ne.nto_orb(jdx,2)) then
                           current_irrep = nto_orb(jdx,2)
                           exit
                      end if
                      up_lim=up_lim+1
                  end do
                  allocate(orb_list(idx)%list(up_lim-low_lim))
                  orb_list(idx)%n_orbs=up_lim-low_lim
                  orb_list(idx)%list(:)=nto_orb(low_lim:up_lim-1,1)
                  low_lim=up_lim
              end do
          end if
!
    END SUBROUTINE nto_plot_merge
  END SUBROUTINE parse_nto_to_mo
#endif
  !----------------------!
  ! END MODIFICATION NTO !
  !----------------------!
  !*************************************************************
  subroutine orbital_plot_check()
    !
    ! Check  the input  this is  done in  this  seperatate subroutine,
    ! because it can take place only after the symmetry part.
    !
    use datatype
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: i_input, alloc_stat, i_irrep, i_orb, i
    type(arrmat3log), allocatable :: help_log(:)

    if(n_input_lines==0) return

    if(density_plot.and.spin_difference.and.ssym%n_spin==1) &
         call error_handler('orbital_plot_check: spin_restricted must be false to calculate spin difference')
    allocate(help_log(ssym%n_irrep),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    do i_irrep=1,ssym%n_irrep
       allocate(help_log(i_irrep)%m(ssym%dim(i_irrep),ssym%partner(i_irrep),ssym%n_spin),&
            stat=alloc_stat)
       ASSERT(alloc_stat==0)

       help_log(i_irrep)%m=.false.
    end do
    do i_input=1,n_input_lines
       if(orb_list(i_input)%irrep_name/=df_irrep_name) then ! find out number of irrep
          do i = 1 , ssym%n_irrep
             if(adjustl(orb_list(i_input)%irrep_name)==adjustl(ssym%name(i))) then
                if(orb_list(i_input)%irrep == df_irrep) then
                   orb_list(i_input)%irrep = i
                else
                   if(orb_list(i_input)%irrep /= i) &
                        call error_handler('orbital_plot_check: incosistency between irrep and irrep name')
                end if
             end if
          end do
       end if
       if(orb_list(i_input)%irrep<1.or.&
           (orb_list(i_input)%irrep> ssym%n_irrep)) &
           call error_handler(&
           'orbital_plot_check: wrong value for irrep')
       if(orb_list(i_input)%partner<1.or.&
            (orb_list(i_input)%partner> ssym%partner(orb_list(i_input)%irrep))) &
           call error_handler(&
           'orbital_plot_check: wrong value for partner')
       if(orb_list(i_input)%spin<1.or.orb_list(i_input)%spin>ssym%n_spin) &
            call error_handler&
            ('orbital_plot_check: wrong value for spin')
       do i_orb=1,orb_list(i_input)%n_orbs
          if(help_log(orb_list(i_input)%irrep)%m&
               (orb_list(i_input)%list(i_orb),orb_list(i_input)%partner,orb_list(i_input)%spin)) &
               call error_handler( &
               'orbital_plot_check: sorry you try to plot the same orbital two times')
          help_log(orb_list(i_input)%irrep)%m&
               (orb_list(i_input)%list(i_orb),orb_list(i_input)%partner,orb_list(i_input)%spin)=.true.
       end do
    end do

    do i_irrep=1,ssym%n_irrep
       deallocate(help_log(i_irrep)%m, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end do
    deallocate(help_log, stat=alloc_stat)
    ASSERT(alloc_stat==0)
  end subroutine orbital_plot_check


  subroutine local_grid (n, off, len)
    !
    ! For grid size  n compute the local portion  to be processed by
    ! this  worker. With  base-0 indices  the local  portion  is the
    ! interval [off, off + len).
    !
    use comm, only: comm_size, comm_rank
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: off, len
    ! *** end of interface ***

    integer :: rank, np
    integer :: F!(n, np, rank), statment function

    ! Statement function definition, note that F(n, np, np) == n:
    F(n, np, rank) = (n * rank) / np

    ! Number of workers and the local rank: 0 <= rank < np:
    np = comm_size()
    rank = comm_rank()

    ! The grid for this worker starts at "off" with an extent given by
    ! "len":
    off = F(n, np, rank)
    len = F(n, np, rank + 1) - F(n, np, rank)
  end subroutine local_grid


  pure function point (ijk, N, mesh)
    !
    ! Returns i * mesh(:, 1) + j * mesh(:, 2) + k * mesh(:, 3)
    !
    implicit none
    integer, intent(in) :: ijk, N(3)
    real(r8_kind), intent(in) :: mesh(3, 3)
    real(r8_kind) :: point(3)
    ! *** end of interface ***

    integer :: i, j, k

    call map3 (ijk, N, i, j, k)

    point = i * mesh(:, 1) + j * mesh(:, 2) + k * mesh(:, 3)
  end function point


  pure subroutine map2 (ij, N, i, j)
    !
    ! Convert linear index ij -> [i, j] such that ij = i + N * j, 0 <=
    ! i < N.
    !
    implicit none
    integer, intent(in) :: ij, N
    integer, intent(out) :: i, j
    ! *** end of interface ***

    j = ij / N
    i = ij - j * N
  end subroutine map2


  pure subroutine map3 (ijk, N, i, j, k)
    !
    ! Convert linear index ijk -> [i, j, k] such that ijk = i + N(1) *
    ! (j + N(2) * k ), 0 <= i < N(1), 0 <= j < N(2).
    !
    implicit none
    integer, intent(in) :: ijk, N(3)
    integer, intent(out) :: i, j, k
    ! *** end of interface ***

    ! NOTE: N(3) is not used.
    integer :: ij

    !  ijk == ij + k * N(1) * N(2):
    call map2 (ijk, N(1) * N(2), ij, k)

    ! ij == i + j * N(1)
    call map2 (ij, N(1), i, j)
  end subroutine map3


  subroutine orbital_plot_bcast()
    ! Purpose: transfer options from master to slaves
    use comm, only: comm_bcast, comm_rank
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: alloc_stat, i_input

    call comm_bcast(density_plot)
    call comm_bcast(density_tot)
    call comm_bcast(spin_difference)
    call comm_bcast(n_input_lines)

    if( comm_rank() /= 0 .and. n_input_lines > 0 ) then
      allocate(orb_list(n_input_lines), stat=alloc_stat)
      ASSERT(alloc_stat==0)
    end if

    call comm_bcast( x0 )
    call comm_bcast( x1 )
    call comm_bcast( x2 )
    call comm_bcast( x3 )
    call comm_bcast( RES )
    ! Unpacking alternative options
    call comm_bcast( alternative_options )
    call comm_bcast( calc_rho )
    call comm_bcast( calc_rho_fit )
    call comm_bcast( calc_xc_pot )
    call comm_bcast( calc_xc_pot_mda )

    do i_input = 1, n_input_lines

      call comm_bcast (orb_list(i_input)%n_orbs)
      if (comm_rank() /= 0) then
        allocate (orb_list(i_input)%list(orb_list(i_input)%n_orbs), &
             stat=alloc_stat )
        ASSERT(alloc_stat==0)
      endif

      call comm_bcast (orb_list(i_input)%irrep)
      call comm_bcast (orb_list(i_input)%partner)
      call comm_bcast (orb_list(i_input)%spin)
      call comm_bcast (orb_list(i_input)%list)
      call comm_bcast (orb_list(i_input)%split_mode)

      if (orb_list(i_input)%split_mode) then

        call comm_bcast (orb_list(i_input)%n_split_atoms)

        if (comm_rank() /= 0) then
          allocate (orb_list(i_input)%split_list(orb_list(i_input)%n_split_atoms), &
               stat=alloc_stat)
          ASSERT(alloc_stat==0)
        endif

        call comm_bcast (orb_list(i_input)%split_list)
      end if
    end do
  end subroutine orbital_plot_bcast
  !*************************************************************

  !*************************************************************
  subroutine orbital_plot_write(iounit)
    !
    ! Purpose: echo input.
    !
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    !** End of interface *****************************************

    integer(i4_kind) :: i_input

    word_format = '("    ",a," = ",a6:" # ",a)' ! including quotes

    call start("ORBITAL_PLOT" ,"ORBITAL_PLOT_WRITE", &
         iounit,operations_echo_input_level)
    call intg("N_INPUT_LINES  ", n_input_lines,  df_n_input_lines   )
    call flag("DENSITY_PLOT   ", density_plot,   df_density_plot    )
    call flag("DENSITY_TOT    ", density_tot,    df_density_tot     )
    call flag("SPIN_DIFFERENCE", spin_difference,df_spin_difference )
    call flag("CALC_RHO       ", calc_rho,       df_calc_rho        )
    call flag("CALC_RHO_FIT   ", calc_rho_fit,   df_calc_rho_fit    )
    call flag("CALC_XC_POT    ", calc_xc_pot,    df_calc_xc_pot     )
    call flag("CALC_XC_POT_MDA", calc_xc_pot_mda,df_calc_xc_pot_mda )
    if (echo()) then
       write(iounit,'(4X,A                 )') '# the origin:'
       write(iounit,'(4X,A,"=",3(F13.8,","))') 'X0', X0
       write(iounit,'(4X,A                 )') '# 1. corner:'
       write(iounit,'(4X,A,"=",3(F13.8,","))') 'X1', X1
       write(iounit,'(4X,A                 )') '# 2. corner:'
       write(iounit,'(4X,A,"=",3(F13.8,","))') 'X2', X2
       write(iounit,'(4X,A                 )') '# 3. corner:'
       write(iounit,'(4X,A,"=",3(F13.8,","))') 'X3', X3
       write(iounit,'(4X,A,"=",3(I4   ,","))') 'RES', RES
    end if
    call stop(empty_line=.false.)
    write(iounit,'()')

    do i_input=1,n_input_lines
       call start("ORBITAL_LIST","ORBITAL_PLOT_WRITE", &
            iounit,operations_echo_input_level)
       call flag("RANGE        ",orb_list(i_input)%range,df_range)
       call intg("IRREP        ",orb_list(i_input)%irrep,df_irrep)
       call word("IRREP_NAME   ",orb_list(i_input)%irrep_name,df_irrep_name)
       call intg("PARTNER      ",orb_list(i_input)%partner,df_partner)
       call intg("SPIN         ",orb_list(i_input)%spin,df_spin)
       if (orb_list(i_input)%range) then
       ! df_n_orb is automatically set to orb_list(i_input)%n_orbs
       call intg("N_ORB        ",orb_list(i_input)%n_orbs, &
                                 orb_list(i_input)%n_orbs)
       else
       call intg("N_ORB        ",orb_list(i_input)%n_orbs,df_n_orb)
       endif
       call flag("SPLIT_MODE   ",orb_list(i_input)%split_mode,df_split_mode)
       call intg("N_SPLIT_ATOMS",orb_list(i_input)%n_split_atoms,df_n_split_atoms)
       call stop(empty_line=.false.)
       1000 format(2I4," # index of the first and the last orbital")
       2000 format(" # List of orbitals"/(20I4))
       3000 format(" # List of atoms"/(20I4))
       if (echo()) then
          if(orb_list(i_input)%range) then
             write(iounit,1000) orb_list(i_input)%list(1),&
                  orb_list(i_input)%list(size(orb_list(i_input)%list))
          else
             write(iounit,2000) orb_list(i_input)%list
          end if
          if(orb_list(i_input)%split_mode) &
               write(iounit,3000) orb_list(i_input)%split_list
          write(iounit,'()')
       end if
    end do
  end subroutine orbital_plot_write
  !*************************************************************


  !*************************************************************
  subroutine orbital_plot_setup(vec_length)
    !
    ! Purpose: do some preparations, do allocation work.
    !
    use xc_cntrl, only: xc_nl_calc
    use orbital_module, only: orbital_setup, orbital_allocate
    implicit none
    integer(i4_kind), intent(in) :: vec_length
    !** End of interface *****************************************

    integer(i4_kind) :: i_input, i_unique, counter, &
         i_atom, i_l, alloc_stat, n_orbs_split, i_list, i_ir, i_spin
    integer :: grid_offset ! set but not used
    integer :: n_orb_tot ! Total number of orbitals to plot INCLUDING
                         ! the splitted orbitals

    ! now we decide how many orbitals  we have to plot in total in the
    ! case of a density plot  we have set and allocate the appropriate
    ! datastructures
    if (density_tot) then
       n_input_lines = symmetry_data_n_irreps() * symmetry_data_n_spin()
       allocate(orb_list(n_input_lines),stat=alloc_stat)
       ASSERT(alloc_stat==0)

       counter=1
       do i_ir=1,symmetry_data_n_irreps()
          do i_spin=1,symmetry_data_n_spin()
             orb_list(counter)%split_mode=.false.
             orb_list(counter)%irrep=i_ir
             orb_list(counter)%partner=1
             orb_list(counter)%spin=i_spin
             orb_list(counter)%n_orbs=symmetry_data_dimension(i_ir)
             allocate(orb_list(counter)%list(symmetry_data_dimension(i_ir)), &
                  stat=alloc_stat)
             ASSERT(alloc_stat==0)

             do i_list=1,symmetry_data_dimension(i_ir)
                orb_list(counter)%list(i_list)=i_list
             end do
             counter=counter+1
          end do
       end do
    end if
    n_orb_tot=0
    do i_input=1,n_input_lines
       if(orb_list(i_input)%split_mode) then
          n_orbs_split=0
          do i_atom=1,orb_list(i_input)%n_split_atoms
             i_unique=orb_list(i_input)%split_list(i_atom)
             do i_l=0,unique_atoms(i_unique)%lmax_ob
                if(unique_atoms(i_unique)%symadapt_partner&
                     (orb_list(i_input)%irrep,i_l)%N_independent_fcts>0) &
                     n_orbs_split=n_orbs_split+1
             end do
          end do
          orb_list(i_input)%n_orbs_split=n_orbs_split
          n_orb_tot=n_orb_tot+n_orbs_split*orb_list(i_input)%n_orbs
       else
          n_orb_tot=n_orb_tot+orb_list(i_input)%n_orbs
       end if
    end do

    ! allocate and populate reverse map:
    allocate(list_orb(n_orb_tot)) ! FIXME: deallocate
    call make_list_orb(orb_list,list_orb)

    ! This sets the global var  "grid_length" to the size of the local
    ! portion of the grid:
    call local_grid (product (RES), grid_offset, grid_length)

    ! allocate space for all plotted orbitals and all gridpoints
    call orbital_setup(vec_length)
    call orbital_allocate(orbs_ob=orbs_ob)

    if (density_plot) then
       if(.not. alternative_options) then
          allocate(rho(grid_length), stat=alloc_stat)
          ASSERT(alloc_stat==0)

          rho=0.0_r8_kind
       else ! ALLOCATIONS FOR AN ALTERNATIVE PLOTTING
          if(calc_rho) then
             allocate(rho_s(grid_length,symmetry_data_n_spin()), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             rho_s = 0.0_r8_kind
          end if
          if(calc_rho_fit) then
             allocate(rho_fit(grid_length,symmetry_data_n_spin()), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             rho_fit = 0.0_r8_kind
          end if
          if(calc_xc_pot) then
             if(xc_nl_calc)then
                ABORT('GGA not yet')
             endif
             allocate(xc_pot(grid_length,symmetry_data_n_spin()), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             allocate(fxc(grid_length), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             allocate(gamm(grid_length,2*symmetry_data_n_spin()-1), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             if(xc_nl_calc) then
                allocate(xcpotgr(grid_length,2*symmetry_data_n_spin()-1), stat=alloc_stat)
                ASSERT(alloc_stat==0)

                xcpotgr(:,:) = 0.0_r8_kind
             end if
             xc_pot       = 0.0_r8_kind
          end if
          if(calc_xc_pot_mda) then
             allocate(xc_pot_mda(grid_length,symmetry_data_n_spin()), stat=alloc_stat)
             ASSERT(alloc_stat==0)

             xc_pot_mda = 0.0_r8_kind
          end if
       end if
    else
       allocate(phi(grid_length,n_orb_tot),stat=alloc_stat)
       ASSERT(alloc_stat==0)

       phi=0.0_r8_kind
    end if
    ! now all preperations are finished
  end subroutine orbital_plot_setup

  subroutine make_list_orb(orb_list, list_orb)
    implicit none
    type(inp), intent(in)  :: orb_list(:)
    type(orb), intent(out) :: list_orb(:)
    ! *** end of interface ***

    integer(i4_kind) :: i_orb
    integer(i4_kind) :: i_input,o

    n_orb = size(list_orb)

    i_orb = 0
    do i_input = 1, size(orb_list)
       if (orb_list(i_input)%split_mode) then
          ABORT('need update')
       else
          do o = 1, orb_list(i_input)%n_orbs
             i_orb = i_orb + 1
             ASSERT(i_orb <= n_orb)
             list_orb(i_orb)%principal = orb_list(i_input)%list(o)
             ! all "orbitals" or "contributions to orbitals" have same
             list_orb(i_orb)%irrep     = orb_list(i_input)%irrep
             list_orb(i_orb)%partner   = orb_list(i_input)%partner
             list_orb(i_orb)%spin      = orb_list(i_input)%spin
          enddo
       end if
    enddo

    call make_orb_names(list_orb)
  end subroutine make_list_orb

  subroutine make_orb_names(list)
    implicit none
    type(orb), intent(inout) :: list(:)
    ! *** end of interface ***

    integer(i4_kind) :: i_orb, n, irr, p, s
    character(len=2), parameter :: ud(2) = (/'up','dn'/)
    character(len=32)           :: name
    character(len=4)            :: symm, part, spin

    do i_orb=1,size(list)
       n   = list(i_orb)%principal
       irr = list(i_orb)%irrep
       p   = list(i_orb)%partner
       s  = list(i_orb)%spin
       symm = symmetry_data_irrepname(irr)

       if(symmetry_data_n_partners(irr) > 1 )then
          part = '(' // trim(adjustl(itos(p))) // ')'
       else
          part = ''
       endif

       if(symmetry_data_n_spin() > 1 )then
          spin = ', ' // ud(s)
       else
          spin = ''
       endif

       name = &
            trim(adjustl(itos(n))) // trim(adjustl(symm))    // &
            trim(part)             // trim(spin)

       list(i_orb)%name = name
    enddo
    ! === end of code ===
  contains
    function itos(i)
      integer(i4_kind), intent(in) :: i
      character(len=4)                  :: itos
      ! *** end of interface ***

      character(len=4) :: buf
      write(buf,'(I4)') i
      itos = buf
    end function itos
  end subroutine make_orb_names
  !*************************************************************


  !*************************************************************
  subroutine orbital_plot_close()
    ! Purpose: do the shutdown work
    use orbital_module, only: orbital_free, orbital_shutdown
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: alloc_stat

    call orbital_free()
    call orbital_shutdown()

    if (density_plot) then
       if(.not. alternative_options) then
          deallocate(rho, stat=alloc_stat)
          ASSERT(alloc_stat==0)

       else ! DEALLOCATE ALTERNATIVE ARRAYS
          if(calc_rho) then
             deallocate(rho_s, stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          if(calc_rho_fit) then
             deallocate(rho_fit, stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          if(calc_xc_pot) then
             deallocate(xc_pot, stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
          if(calc_xc_pot_mda) then
             deallocate(xc_pot_mda, stat=alloc_stat)
             ASSERT(alloc_stat==0)
          end if
       end if
    else
          deallocate(phi, stat=alloc_stat)
          ASSERT(alloc_stat==0)
    end if

    ! Recursive   deallocations   of   structures   with   allocatable
    ! components:
    if (allocated(orb_list)) then
       deallocate(orb_list,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if

    if (allocated(list_orb)) then
       deallocate(list_orb,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if
  end subroutine orbital_plot_close
  !*************************************************************

  !*************************************************************
  subroutine orbital_plot_fit_bcast()
    ! Purpose: send the fitting coefficients to slaves
    !
    use comm, only: comm_size, comm_rank
    implicit none
    !** End of interface *****************************************

    if (comm_size() <= 1) return

    ! FIXME: convert these to MPI_BCAST/comm_bcast:
    if (comm_rank() == 0) then
       call orbital_plot_fit_send()
    else
       call orbital_plot_fit_receive()
    endif
  end subroutine orbital_plot_fit_bcast
  !*************************************************************

  !*************************************************************
  subroutine orbital_plot_fit_send()
    !
    ! Purpose: send the fitting coefficients to slaves
    !
    use msgtag_module
    use comm_module
    use fit_coeff_module
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: info , ch_len, s

    ch_len = fit_coeff_n_ch()
    call comm_init_send(comm_all_other_hosts,msgtag_packed_message)
    call commpack(ch_len,info)
    ASSERT(allocated(coeff_charge))
    call commpack(coeff_charge(1),ch_len,1,info)
    if(info /= 0) &
         call error_handler("orbital_plot_fit_send : packing coeff_charge")
    if( calc_xc_pot_mda ) then
       do s = 1, symmetry_data_n_spin()
          call commpack(coeff_xcmda(1,s), ch_len,1,info)
          if(info /= 0) &
               call error_handler("orbital_plot_fit_sevd : packing coeff_xcmda(s)")
       end do
    end if
    if( symmetry_data_n_spin() > 1 ) then
       ASSERT(allocated(coeff_spin))
       call commpack(coeff_spin(1),  ch_len,1,info)
       if(info /= 0) &
            call error_handler("orbital_plot_fit_send : packing coeff_spin")
    end if
    call comm_send()
 end subroutine orbital_plot_fit_send
  !*************************************************************

 !*************************************************************
  subroutine orbital_plot_fit_receive()
    !
    ! Purpose: Receives fitting coefficients from master
    !
    use msgtag_module
    use comm_module
    use fit_coeff_module
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: info, ch_len, s, alloc_stat

    call comm_save_recv(comm_master_host,msgtag_packed_message)
    call communpack(ch_len,info)
    if(info /= 0) &
         call error_handler("Orbital_plot_fit_receive: unpacking ch_len")
    if(.not.allocated(coeff_charge) ) then
       ! Just to allow compilation
       call fit_coeff_initialize()
!      allocate(coeff_charge(ch_len), stat=alloc_stat)
!      if(alloc_stat /= 0) &
!           call error_handler("Orbital_plot_fit_receive: coeff_charge allocation failed")
    end if
    call communpack(coeff_charge(1),ch_len,1,info)
    if(info /= 0) &
         call error_handler("orbital_plot_fit_receive : unpacking coeff_charge")
    if( calc_xc_pot_mda ) then
       if(.not.allocated(coeff_xcmda)) then
          allocate(coeff_xcmda(ch_len,symmetry_data_n_spin()), stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end if
       do s = 1, symmetry_data_n_spin()
          call communpack(coeff_xcmda(1,s), ch_len,1,info)
          if(info /= 0) &
               call error_handler("orbital_plot_fit_receive : unpacking coeff_xcmda(s)")
       end do
    end if
    if(symmetry_data_n_spin() > 1 ) then
       if(.not.allocated(coeff_spin) ) then
          allocate(coeff_spin(ch_len), stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end if
       call communpack(coeff_spin(1), ch_len,1,info)
       if(info /= 0) &
            call error_handler("orbital_plot_fit_receive : unpacking coeff_spin")
    end if
  end subroutine orbital_plot_fit_receive
  !*************************************************************

  !*************************************************************
  subroutine oplot_collect (rho, th)
    !
    ! Receive plooted  orbitals from the  slaves and write the  to the
    ! punchfile in proper order.
    !
    use comm, only: comm_size, comm_rank, comm_send, comm_recv
    implicit none
    real(r8_kind), intent(in) :: rho(:)
    integer(i4_kind), intent(in) :: th ! tapehandle for writing
    !** End of interface *****************************************

    integer :: rank, vl

    ! Grid shares may differ by one  unit. So in case master has less
    ! than any of the slaves reserve one more:
    real(r8_kind) :: buffer(size(rho) + 1)

    ASSERT(size(rho)==grid_length)

    ! NOTE: These sends and receives use the default MPI tag:
    if (comm_rank() == 0) then
       ! first we write master contribution
       call oplot_write (rho, th)

       ! then get points from slaves, if any:
       do rank = 1, comm_size() - 1
          call comm_recv (vl, rank)
          ASSERT(vl <= size(buffer))
          call comm_recv (buffer(:vl), rank)

          call oplot_write (buffer(:vl), th)
       end do
    else ! and for slaves:
       ! send my portion to master, all of it:
       call comm_send (size(rho), 0)
       call comm_send (rho, 0)
    endif
  end subroutine oplot_collect

  ! contents ot tape
  ! 1. title: not implemented yet
  ! 2. grid_parameters: x0,x1,x2 origin of grid and two orthogonal
  !                     direction vectors
  !                     step_x,step_y,length_x,length_y
  ! 3. density_plot
  !    if true
  !       4. density on the grid
  !    else
  !       4. n_orb_tot :number of orbitals
  !       5. for every orbital: name of the irrep ! not yet implemented
  !                             spin
  !                             number of orbital
  !                             split_mode
  !                             n_split_orbs
  !       6. for every orbital and every splitted contribution:
  !                            value of orbital on the grid

  !*************************************************************
  subroutine oplot_startwrite (filename, io)
    !
    ! Purpose: open the file for writing grid data.
    !
    use comm, only: comm_rank
    use iounitadmin_module, only: openget_iounit
    use filename_module, only: outfile
    use unique_atom_module , only: unique_atoms
    use symm_positions, only: symm_nabor
    use punchfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in)  :: filename
    integer(i4_kind), intent(out) :: io ! tapehandle for writing
    !** End of interface *****************************************

    if (comm_rank() /= 0) then
       io = -1
       RETURN
    endif

    ! Dont  specify status,  overwrite the  file, if  it  does already
    ! exist:
    io = openget_iounit (FILE=trim(outfile(filename)), FORM='formatted')

    call pun_title (io, 'This file was generated by ParaGauss')
    call pun_block_header (io, 'fragment', 0)
    call pun_title (io, 'molecule')
    call symm_nabor (io, unique_atoms, 1)
  end subroutine oplot_startwrite

  subroutine oplot_startblock (title, io)
    ! Purpose: so far just close the file after e.g. all orbitals have
    !          been written
    use punchfile
    implicit none
    character(len=*), intent(in) :: title
    integer(i4_kind), intent(in) :: io ! tapehandle for writing
    ! *** end of interface ***

    real(r8_kind)    :: D(3)
    integer(i4_kind) :: NPT

    if (io < 0) RETURN

    call pun_block_header (io, 'data', 0)
    call pun_grid_title (io, title)
    D(1) = SQRT (SUM ((X1 - X0)**2)) ! length_x
    D(2) = SQRT (SUM ((X2 - X0)**2)) ! length_y
    D(3) = SQRT (SUM ((X3 - X0)**2)) ! length_z
    call pun_grid_axes (io, RES, D)
    call pun_grid_mapping (io, X0, X1, X2, X3)
    NPT = PRODUCT (RES)
    call pun_block_header (io, 'grid_data', NPT, 'elements', 1)
  end subroutine oplot_startblock

  subroutine oplot_write (buf, io)
    ! Purpose: so far just close the file after e.g. all orbitals have
    !          been written
    use punchfile
    implicit none
    real(r8_kind), intent(in) :: buf(:)
    integer(i4_kind), intent(in) :: io ! tapehandle for writing
    ! *** end of interface ***

    if (io < 0) RETURN

    call pun_grid_data (io, buf)
  end subroutine oplot_write

  subroutine oplot_stopwrite (io)
    ! Purpose: so far just close the file after e.g. all orbitals have
    !          been written
    use iounitadmin_module, only: returnclose_iounit
    implicit none
    integer(i4_kind), intent(in) :: io ! tapehandle for writing
    ! *** end of interface ***

    if (io < 0) RETURN

    ! so far just close the file
    call returnclose_iounit (io)
  end subroutine oplot_stopwrite
  !*************************************************************

  !*************************************************************
  subroutine orbital_plot_main()
    !
    ! Purpose: main routine for plotting orbitals.
    !
    use orbitalprojection_module, only: orbitalprojection_ob
    use gamma_module, only: gamma_is_closed, gamma_setup
    use machineparameters_module, only: vec_length => machineparameters_veclen
    use orbital_module, only: orbital_calculate, &
         fit_fct_allocate, fit_fct_calculate, fit_fct_free
    use eigen_data_module, only: eigvec
    use occupied_levels_module, only: occ_num_occ
    use density_calc_module, only: fitted_density_is_closed, &
         fitted_density_calc, fitted_density_calc_setup, &
         fitted_density_calc_close
    use fit_coeff_module, only: coeff_xcmda
    use xc_cntrl, only: xc_nl_calc
    use xc_func, only: xc_functionals
#ifdef WITH_RESPONSE
    use nto_plot_module, only: nto_plot_grid, nto_coeff
#endif
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: orb_index, grid_start, grid_end, &
         i_spin, i_partner, i_orb, i_ir, i_input, vla, m, &
         index, i_unique, m_min, m_max, n_ind, i_atom, i_l
    integer(i4_kind) :: spin
    real(r8_kind), pointer :: orb(:,:,:), eigv(:,:,:)
    real(r8_kind) :: occ_real
    character(len=80) :: title
#ifdef WITH_RESPONSE
    real(r8_kind), allocatable :: nphi(:)       ! ADDED FOR NTO
    character(len=10) :: label, label1          ! ADDED FOR NTO
#endif
    real(r8_kind), pointer :: p_rho(:,:), p_fxc(:), p_vxc(:,:), &
         p_xcpotgr(:,:), help_array(:,:)

    ! relatively small automatic arrays:
    real(r8_kind) :: phi_help(vec_length)
    real(r8_kind) :: grdpts(vec_length, 3) ! gridpoints
    real(r8_kind) :: mesh(3, 3)      ! mesh (cell) params
    integer :: i, off, len, ijk

    ! This also calls orbital_setup() and allocates global orbs_ob:
    call orbital_plot_setup (vec_length)

    ! in case of fitted density  and related potentials are plotted we
    ! need to have fitting coefficients on slaves (again).
    if(calc_rho_fit .or. calc_xc_pot) then
       call orbital_plot_fit_bcast()
    end if

    ! Get the  range [off,  off +  len) for the  local portion  of the
    ! grid:
    call local_grid (product (RES), off, len)

    if (RES(1) > 1) then
       mesh(:, 1) = (X1 - X0) / (RES(1) - 1)
    else
       mesh(:, 1) = 0.0
    endif

    if (RES(2) > 1) then
       mesh(:, 2) = (X2 - X0) / (RES(2) - 1)
    else
       mesh(:, 2) = 0.0
    endif

    if (RES(3) > 1) then
       mesh(:, 3) = (X3 - X0) / (RES(3) - 1)
    else
       mesh(:, 3) = 0.0
    endif

    ! print *, "grid size=", product (RES), "local len=", len, "off=", off
    grid_start = 1

    ! loop over gridpoints:
    do while (grid_start <= len) ! fetching part of the grid

       ! Process a batch  of grid points, no more  than vla.  This may
       ! be,  in  general,  below  vec_length,  particularly  in  last
       ! iteration:
       vla = min (1 + len - grid_start, vec_length)
       grid_end = grid_start + vla - 1

       ! print *, "processing [", off + grid_start, off + grid_start + vla, ")"
       do i = 1, vla
          ! First  argument  to  point(),  the grid  point  index,  is
          ! base-0:
          ijk = off + (grid_start - 1) + (i - 1)
          grdpts(i, :) = X0 + point (ijk, RES, mesh)
       enddo

       ! calculate primitive, contracted, symmetryadapted orbitals
       call orbital_calculate (grdpts(1:vla, 1:3), vla, orbs_ob)
       ! now build actual orbitals
       orb_index=1
       do i_input = 1, n_input_lines
          i_ir=orb_list(i_input)%irrep
          i_partner=orb_list(i_input)%partner
          i_spin=orb_list(i_input)%spin

          ! check user input:
          ASSERT(i_ir>0)
          ASSERT(i_ir<=size(eigvec))

          orb  => orbs_ob(i_ir)%o
          eigv => eigvec(i_ir)%m
          associate( occ => occ_num_occ(i_ir)%m )

          ASSERT(i_partner>0)
          ASSERT(i_partner<=size(orb,3))
          ASSERT(i_spin>0)
          ASSERT(i_spin<=size(eigv,3))

          if(orb_list(i_input)%split_mode) then
             do i_orb=1,orb_list(i_input)%n_orbs
                index=orb_list(i_input)%list(i_orb)
                ! check user input:
                ASSERT(index>0)
                !----------------------!
                ! MODIFICATION FOR NTO !
                !----------------------!
                if(nto) then
                else
                    ASSERT(index<=size(eigv,1))
                end if
                !---------------------------!
                !  END MODIFICATION FOR NTO !
                !---------------------------!
                do i_atom=1,orb_list(i_input)%n_split_atoms
                   i_unique=orb_list(i_input)%split_list(i_atom)
                   do i_l=0,unique_atoms(i_atom)%lmax_ob
                      n_ind=unique_atoms(i_unique)%symadapt_partner&
                           (orb_list(i_input)%irrep,i_l)%N_independent_fcts
                      if(n_ind>0) then
                         m_min=orbitalprojection_ob(i_ir,i_l,i_unique)
                         m_max=m_min-1+n_ind*&
                              (unique_atoms(i_unique)%l_ob(i_l)%N_contracted_fcts+&
                              unique_atoms(i_unique)%l_ob(i_l)%N_uncontracted_fcts)
                         do m=m_min,m_max
                            phi(grid_start:grid_end,orb_index)=&
                                 phi(grid_start:grid_end,orb_index)+&
                                 orb(1:vla,m,i_partner)*eigv(m,index,i_spin)
                         end do
                         orb_index=orb_index+1
                      end if
                   end do
                end do
             end do
          else  ! non split-mode
             do i_orb=1,orb_list(i_input)%n_orbs
                index=orb_list(i_input)%list(i_orb)
                ! check user input:
                ASSERT(index>0)
                !----------------------!
                ! MODIFICATION FOR NTO !
                !----------------------!
                if(nto) then
                else
                    ASSERT(index<=size(eigv,1))
                end if
                !----------------------------!
                !  END MODIFICATION FOR NTO  !
                !----------------------------!
                if (density_plot) then
                   if (index<=size(occ,1)) then
                     ! do i_spin=1,ssym%n_spin
                      occ_real = occ(index,i_spin) / ssym%partner(i_ir)
                      do i_partner=1,ssym%partner(i_ir)
                         phi_help=0.0_r8_kind
                         do m=1,ssym%dim(i_ir)
                            phi_help(1:vla)=phi_help(1:vla)+&
                                 orb(1:vla,m,i_partner)*eigv(m,index,i_spin)
                         end do
                         phi_help(1:vla)=phi_help(1:vla)&
                              *phi_help(1:vla)*occ_real
                         if(.not. alternative_options) then
                            if(i_spin==2.and.spin_difference) then
                               rho(grid_start:grid_end)=rho(grid_start:grid_end)-&
                                    phi_help(1:vla)
                            else
                               rho(grid_start:grid_end)=rho(grid_start:grid_end)+&
                                    phi_help(1:vla)
                            end if
                         else
                            if(calc_rho) then
                               rho_s(grid_start:grid_end,i_spin) = rho_s(grid_start:grid_end,i_spin) + &
                                    phi_help(1:vla)
                            end if
                         end if
                      end do! loop over partners
                   end if
                else
                   do m=1,ssym%dim(i_ir)
                      phi(grid_start:grid_end,orb_index)=&
                           phi(grid_start:grid_end,orb_index)+&
                           orb(1:vla,m,i_partner)*eigv(m,index,i_spin)
                   end do
                end if
                orb_index=orb_index+1
             end do! loop over orbitals
          end if! split mode
          end associate
       end do! loop over inputlines
       if(vla == 0) exit

       ! calculate fitted density , xc-potentials
       spin = symmetry_data_n_spin()
       if( calc_xc_pot) then
          p_rho => rho_s (grid_start:grid_end, 1:spin )
          p_fxc => fxc   (grid_start:grid_end)
          p_vxc => xc_pot(grid_start:grid_end, 1:spin )
          p_fxc = 0.0_r8_kind
          p_vxc = 0.0_r8_kind
          if(xc_nl_calc) then
             p_xcpotgr => xcpotgr(grid_start:grid_end,:)
             call xc_functionals(vla,spin,p_rho,p_fxc,p_vxc,gamm,p_xcpotgr)
          else
             call xc_functionals(vla,spin,p_rho,p_fxc,p_vxc)
          endif
       end if
       if( calc_rho_fit ) then
          if(fitted_density_is_closed()) call fitted_density_calc_setup()
          call fit_fct_allocate("ch",fcts=fcts_p)
          call fit_fct_calculate(grdpts,vla,"ch",fcts=fcts_p)
          help_array => rho_fit(grid_start:grid_end,1:symmetry_data_n_spin())
          call fitted_density_calc(vla,rho=help_array,fcts=fcts_p)
          call fit_fct_free("ch",fcts=fcts_p)
          if(.not.fitted_density_is_closed()) call fitted_density_calc_close()
       end if
       if( calc_xc_pot_mda ) then
          if(fitted_density_is_closed()) call fitted_density_calc_setup()

          if(gamma_is_closed()) call gamma_setup(16)
          ! cp. 16 with default numj=17 in old vers. of gamma_module

          call fit_fct_allocate("ch",fcts=fcts_p)
          call fit_fct_calculate(grdpts,vla,"cl",fcts=fcts_p)  ! "ch" ----> "cl"
          help_array => xc_pot_mda(grid_start:grid_end,1:symmetry_data_n_spin())
          call fitted_density_calc(vla,rho=help_array,fcts=fcts_p,altern_coeff=coeff_xcmda)
          call fit_fct_free("ch",fcts=fcts_p)
          if(.not.fitted_density_is_closed()) call fitted_density_calc_close()
       end if
       ! end of fitted density or xc-potentials calculations
       grid_start = grid_end + 1
    end do! loop over gridpoints

    ! now we can receive slave contributions and save the stuff to tape
    ! this does nothing on slave except returns th < 0 :
    call oplot_startwrite ('plot.pun', th)

    if (density_plot) then
       if (.not. alternative_options) then
          call oplot_startblock ('density', th)
          call oplot_collect (rho, th)
       else
          do spin = 1, symmetry_data_n_spin()
             if (calc_rho) then
                call oplot_startblock ('rho_s', th)    !_rho(spin) )
                call oplot_collect (rho_s(:,spin), th) !_rho(spin) )
             end if
             if (calc_rho_fit) then
                call oplot_startblock ('rho_fit', th) !_rho_fit(spin) )
                call oplot_collect (rho_fit(:,spin), th) !_rho_fit(spin) )
             end if
             if (calc_xc_pot) then
                call oplot_startblock ('xc_pot', th) !_xc_pot(spin) )
                call oplot_collect (xc_pot(:,spin), th) !_xc_pot(spin) )
             end if
             if (calc_xc_pot_mda) then
                call oplot_startblock ('xc_pot_mda', th) !_xc_pot_mda(spin) )
                call oplot_collect (xc_pot_mda(:,spin), th) !_xc_pot_mda(spin) )
             end if
          enddo
       end if
    else
       if (.not. nto) then
          do i_orb = 1, size(list_orb)
             title = 'orbital ' // trim(list_orb(i_orb)%name)
             call oplot_startblock (trim(title), th)
             call oplot_collect (phi(:,i_orb), th)
          enddo
       else
       !---------------------!
       !MODIFICATION FOR NTO !
       !---------------------!
#ifndef WITH_RESPONSE
         ABORT('recompile -DWITH_RESPONSE')
#else
         allocate(nphi(size(phi,1)))
         do i_orb=1,size(nto_coeff)
           do i_spin=1,size(nto_coeff(i_orb)%nto_list)
!
               !!---- Grid for occupied NTOs  ----!!
!
               call nto_plot_grid(i_orb    , &           !Calculation of the grid points
                                  i_spin   , &
                                  .true.   , &
                                  phi      , &
                                  nphi        )
               write(label,'(a)')ssym%name(nto_coeff(i_orb)%irrep) ! Final Irrep Name in which the NTO belongs to
               write(label1,'(i4)')nto_coeff(i_orb)%nto_list(i_spin)%index ! Index of the NTO
               if(nto_coeff(i_orb)%spin.eq.1) then
                  if(nto_coeff(i_orb)%trans.eq.1) then
                      title = trim(label)//'(up,SS)/'//trim(label1)//'NTO(OCC)' ! Title for NTO with spin up, SS
                  else
                      title = trim(label)//'(up,ST)/'//trim(label1)//'NTO(OCC)' ! Title for NTO with spin up, ST
                  end if
               else
                  if(nto_coeff(i_orb)%trans.eq.1) then
                      title = trim(label)//'(down,SS)/'//trim(label1)//'NTO(OCC)' ! Title for NTO with spin down, SS
                  else
                      title = trim(label)//'(down,ST)/'//trim(label1)//'NTO(OCC)' ! Title for NTO with spin down, ST
                  end if
               end if
               call oplot_startblock(trim(title), th ) ! Create the title and write it into the pun file
               call oplot_collect( nphi, th ) ! Put the points in the grid file
!
               !!---- Grid for unoccupied NTOs  ----!!
!
               call nto_plot_grid(i_orb    , & ! Calculation of the grid points
                                  i_spin   , &
                                  .false.  , &
                                  phi      , &
                                  nphi        )
               write(label,'(a)')ssym%name(nto_coeff(i_orb)%irrep) ! Final Irrep Name in which the NTO belongs to
               write(label1,'(i4)')nto_coeff(i_orb)%nto_list(i_spin)%index ! Index of the NTO
               if(nto_coeff(i_orb)%spin.eq.1) then
                  if(nto_coeff(i_orb)%trans.eq.1) then
                      title = trim(label)//'(up,SS)/'//trim(label1)//'NTO(UNOCC)' ! Title for NTO with spin up, SS
                  else
                      title = trim(label)//'(up,ST)/'//trim(label1)//'NTO(UNOCC)' ! Title for NTO with spin up, ST
                  end if
               else
                  if(nto_coeff(i_orb)%trans.eq.1) then
                      title = trim(label)//'(down,SS)/'//trim(label1)//'NTO(UNOCC)' ! Title for NTO with spin down, SS
                  else
                      title = trim(label)//'(down,ST)/'//trim(label1)//'NTO(UNOCC)' ! Title for NTO with spin down, ST
                  end if
               end if
               call oplot_startblock(trim(title), th ) ! Create the title and write it into the pun file
               call oplot_collect( nphi, th ) ! Put the points in the grid file
!
            end do
         end do
         deallocate(nphi)
#endif
       !---------------------------!
       ! END MODIFICATION FOR NTO  !
       !---------------------------!
       end if
    end if

    call oplot_stopwrite(th)

    call orbital_plot_close()

#if 0
    if(xsf_output)then
       call xsf_rewrite(th)
    endif
#endif
  end subroutine orbital_plot_main

#if 0
  subroutine xsf_rewrite(th)
    !
    ! Purpose: rewriting plot-files to XSF input
    !
    !** End of interface *****************************************
    integer(i4_kind)           :: th
    type(readwriteblocked_tapehandle) :: th !
    !------------ Declaration of local variables ---------------------
    integer(i4_kind)           :: i_input, &
                                       alloc_stat, &
                                       i_list, &
                                       n_orb_plot, n_points, &
                                       n_points_x,n_points_y,n_points_z,&
                                       xsf_unit,ua,ea
    real(r8_kind)              :: help_name(20),hv(4)
    real(r8_kind), allocatable :: rho_tot(:)
    real(r8_kind), parameter   :: au_A=0.52918_r8_kind
    character(len=200)              :: xsf_file
    !------------ Executable code ------------------------------------

    call readwriteblocked_startread(readwriteblocked_filename(th),th)
    ! Read headers written by oplot_startwrite
    call readwriteblocked_read(help_name,th)
    call readwriteblocked_read(x0,th)
    call readwriteblocked_read(x1,th)
    call readwriteblocked_read(x2,th)
    call readwriteblocked_read(hv,th)
    step_x   = hv(1); step_y   = hv(2)
    length_x = hv(3); length_y = hv(4)
    if(density_plot .or. alternative_options) then
       call readwriteblocked_read(hv(1:1),th)
    else
       call readwriteblocked_read(hv(1:1),th)
       ! Total Number of orbitales contained in the tape
       n_orb_plot=sum(orb_list(:)%n_orbs)
       call readwriteblocked_read(hv(1:1),th)
       do i_input=1, n_input_lines
          do i_list=1,orb_list(i_input)%n_orbs
             call readwriteblocked_read(hv(1:1),th)
             call readwriteblocked_read(hv(1:3),th)
             call readwriteblocked_read(hv(1:1),th)
             call readwriteblocked_read(hv(1:1),th)
          end do
       end do
    end if
    n_points_x = int(length_x/step_x)+1; write(*,*) "N_x ", n_points_x
    n_points_y = int(length_y/step_y)+1; write(*,*) "N_y ", n_points_y
    n_points_z = int(length_z/step_z)+1; write(*,*) "N_z ", n_points_z
!!$    n_points = n_points_x * n_points_y
!!$    if(XSF_3D) n_points = n_points  * n_points_z
    n_points = n_points_x * n_points_y * n_points_z
    write(*,*) "N_  ", n_points
    if(density_plot) then
       allocate(rho_tot(N_points), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       call readwriteblocked_read(rho_tot,th)
    end if
    xsf_file = readwriteblocked_filename(th)
    xsf_file = xsf_file(1:len_trim(xsf_file)-4)//".xsf"
    write(*,*) "XSF_FILE  : ", trim(xsf_file)
    xsf_unit = openget_iounit(trim(xsf_file))

    ! Start writing XSF-file

    write(xsf_unit,*) " ATOMS"
    do ua=1,N_unique_atoms
       do ea=1,unique_atoms(ua)%N_equal_atoms
          write(xsf_unit,"(i3,3x,3(f14.10,3x))") int(unique_atoms(ua)%Z), &
                                                 unique_atoms(ua)%position(:,ea)*au_A
       end do
    end do
    if(xsf_3D) then
       write(xsf_unit,*) " BEGIN_BLOCK_DATAGRID_3D"
       write(xsf_unit,*) " Just a comment"
       write(xsf_unit,*) " DATAGRID_3D_UNKNOWN"
       write(xsf_unit,"(3i6)") n_points_x, n_points_y,n_points_z
       write(xsf_unit,"(3(f16.10))") x0*au_A
       write(xsf_unit,"(3(f16.10))") x1*length_x*au_A
       write(xsf_unit,"(3(f16.10))") x2*length_y*au_A
       write(xsf_unit,"(3(f16.10))") x3*length_z*au_A
    else if(xsf_2D) then
       write(xsf_unit,*) " BEGIN_BLOCK_DATAGRID_2D"
       write(xsf_unit,*) " Just a comment"
       write(xsf_unit,*) " DATAGRID_2D_UNKNOWN"
       write(xsf_unit,"(2i6)") n_points_x, n_points_y
       write(xsf_unit,"(3(f16.10))") x0*au_A
       write(xsf_unit,"(3(f16.10))") x1*length_x*au_A
       write(xsf_unit,"(3(f16.10))") x2*length_y*au_A
    else
       call error_handler("Orbital_plot xsf_rewrite : XSF_2D or XSF_3D not set")
    end if
    write(xsf_unit,"(5(f16.10,3x))") rho_tot
    if(xsf_2D) then
       write(xsf_unit,*) " END_DATAGRID_2D"
       write(xsf_unit,*) " END_BLOCK_DATAGRID_2D"
    end if
    if(xsf_3D) then
       write(xsf_unit,*) " END_DATAGRID_3D"
       write(xsf_unit,*) " END_BLOCK_DATAGRID_3D"
    end if

    if(density_plot) then
       deallocate(rho_tot, stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if
  end subroutine xsf_rewrite
#endif

end module orbital_plot_module


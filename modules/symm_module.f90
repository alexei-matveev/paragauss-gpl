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
module  symm_module
  !-------------------------------------------------------------------
  !
  !  Purpose: runs the symmetry part 
  !           on top of the modules group and efm, where
  !           group provides essential symmetry data and
  !           efm is the central symmetrizer
  !           operates on the nuclear positions
  !           generates symmetry equivalent atoms
  !           and symmetry equivalent distances
  !
  !
  !  Author: MM
  !  Date: 11/96
  !
  !
  !-------------------------------------------------------------------
!== Interrupt of public interface of module =====================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification 
  ! Author: AM
  ! Date:  04-05.99
  ! Description: symm_proj_symmadapt has been divided into several
  !              subs. Generator large_to_small has been added.
  !
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module, only: &
        IK=>i4_kind, &
        RK=>r8_kind, &
        CK=>c16_kind ! type specification parameters
  use error_module
  use group_module
  use efm_module, only: efm_cscoII_diag_invspace ! more in subs !!
  use clebsch_gordan, only: sym_prod ! more in subs !!
  use options_module, only: options_spin_orbit

  implicit none
  save            ! keep contents of module
  private         ! by default, all names are private
!== Interrupt end of public interface of module =================



  !*********************************
  type,private :: vsalcs_ua
     !
     integer(IK)                            :: lmax,n_ea,n_irr
     !<<< it`s always good to have them together
     type(sym_prod),pointer :: mos(:,:)
     !
     ! mos(n_irr,0:lmax), M_olecular O_rbitals  
  end type vsalcs_ua

  type,private :: psalcs_ua
     !
     integer(IK)                            :: lmax,n_ea,n_irr
     !
     type(sym_prod),pointer :: mos(:,:,:)
     !
     ! mos(n_irr,0:lmax,2)
     ! M_olecular O_rbitals, spread over equivalent atoms
  end type psalcs_ua
  !*********************************


  !------------ Declaration of constants and variables ---------------

  !------------ public functions and subroutines ---------------------
  public :: &
       & symm_trafo_alloc,&
       & symm_proj_trafo_alloc,&
       & symm_trafo_dealloc,&
       & symm_symmadapt,&
       & symm_dipole_selrules_gen,&
       & symm_get_pcoupling,&
       & symm_get_cccoupling,&
       & symm_init, symm_done,&
       & symm_adapt_centers


  !===================================================================
  ! End of public interface of module
  !===================================================================

  interface load_salcs
     module procedure load_vsalcs_to_uatom
     module procedure load_psalcs_to_uatom
  end interface

  interface show_fine
     module procedure show_fine_r
     module procedure show_fine_c
  end interface

  !------------ Declaration of constants and variables ---------------

  logical                          :: initialized = .false.

  integer(IK),dimension(:),pointer :: pcoupling
  integer(IK),dimension(:),pointer :: cccoupling

  integer(IK)                      :: LMAX_ALL=0
  
  !-------------------------------------------------------------------
  real(RK),parameter :: epsln   = 1.0E-10_rk
  logical, parameter :: verbose = .false.
  !------------ Subroutines ------------------------------------------
contains

  subroutine symm_init(spin_orbit)
    use unique_atom_module
    use unique_atom_methods, only: unique_atom_calc_lmax
    use group_module, only: group_num_re
    use uatom_symmadapt, only: uatom_symmadapt_init
    implicit none
    logical,optional,intent(in)  :: spin_orbit
    ! *** end of interface ***

    logical :: SO=.false.
    integer :: memstat

    if(present(spin_orbit))SO=spin_orbit

    call unique_atom_calc_lmax()     !<<< from UA module

    LMAX_ALL = unique_atom_lmax_all !<<< from UA module

    call symm_trafo_alloc()

    if(SO)then
       call symm_proj_trafo_alloc()

       allocate(pcoupling(group_num_re),STAT=memstat)
       if(memstat/=0)call error_handler("sm/symm_init: alloc pcoupling failed")
       pcoupling = -1

       allocate(cccoupling(group_num_re),STAT=memstat)
       if(memstat/=0)call error_handler("sm/symm_init: alloc cccoupling failed")
       cccoupling = -1
    endif

    call uatom_symmadapt_init(spor=SO)
    call uatom_symmadapt_init(n_ua=n_unique_atoms)

    initialized = .true.
  end subroutine symm_init

  subroutine symm_done(spin_orbit)
    implicit none
    logical,optional,intent(in)  :: spin_orbit
    ! *** end of interface ***

    logical :: SO=.false.
    integer :: memstat

    if(present(spin_orbit))SO=spin_orbit

    call symm_trafo_dealloc

    if(SO)then
       ! call symm_proj_trafo_dealloc !??? doesnt exist ???

       deallocate(pcoupling,STAT=memstat)
       if(memstat/=0)call error_handler("sm/symm_done: alloc pcoupling failed")

       deallocate(cccoupling,STAT=memstat)
       if(memstat/=0)call error_handler("sm/symm_done: alloc cccoupling failed")
    endif

    initialized = .false.
  end subroutine symm_done

  function symm_get_pcoupling() result(pc)
    use error_module
    implicit none
    integer(IK),pointer :: pc(:) !<<< result
    ! *** end of interface ***

    call error(.not.initialized,&
         & "sm/symm_get_pcoupling: not initialized")

    if(any(pcoupling.eq.-1))&
         & call warn("sm/symm_get_pcoupling: pcoupling is undefined")

    pc => pcoupling
  end function symm_get_pcoupling

  function symm_get_cccoupling() result(ccc)
    use error_module
    implicit none
    integer(IK),pointer :: ccc(:) !<<< result
    ! *** end of interface ***

    call error(.not.initialized,&
         & "sm/symm_get_cccoupling: not initialized")

    if(any(cccoupling.eq.-1))&
         & call warn("sm/symm_get_cccoupling: cccoupling is undefined")

    ccc => cccoupling
  end function symm_get_cccoupling

  !*************************************************************
  subroutine symm_trafo_alloc
    !  Purpose: generate and allocate transformation matrices for all
    !           angular momenta l=0,unique_atom_lmax_all
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    use efm_module
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: l
    ! counter
    integer(IK)                      :: l_max
    ! angular momentum
    integer(IK)                      :: dim_trafo
    ! dimension of transformation matrix
    integer(IK)                      :: status
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    l_max = max(9,LMAX_ALL)

    ! allocate transformation matrices
    allocate(ylm_trafos(0:l_max),stat = status)
    if (status.ne.0) then
       call error_handler("symm_trafo_alloc: allocation failed")
    end if

    ! loop over all angular momenta l needed
    do l=0,l_max
       dim_trafo = 2*l+1
       ! allocate transformation matrices
       allocate(ylm_trafos(l)%matrix(dim_trafo,dim_trafo,group_num_el))
       ! generate transformation matrices
       call efm_trafo_gen(l,ylm_trafos(l)%matrix)
    end do

  end subroutine symm_trafo_alloc
  !*************************************************************

  !*************************************************************
  subroutine symm_proj_trafo_alloc
    !  Purpose: generate and allocate transformation matrices for all
    !           angular momenta l=0,unique_atom_lmax_all
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    use efm_module
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: l,s
    ! counter
    integer(IK)                      :: l_max
    ! angular momentum
    integer(IK)                      :: dim_trafo
    ! dimension of transformation matrix
    integer(IK)                      :: status
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    ! we need at least an angular momentum of 9 because we employ
    ! the angular momenta also for irrep labeling
    l_max = max(9,LMAX_ALL) !<<<mdf

    ! allocate transformation matrices
    allocate(ylm_proj_trafos(0:l_max,2),comp_to_real(0:l_max,2,2),stat = status)
    if (status.ne.0) then
       call error_handler("symm_proj_trafo_alloc: allocation failed")
    end if

    ! loop over all angular momenta l needed
    do l=0,l_max
       do s=1,2
          if ((l.eq.0).and.(s.eq.1)) then
             cycle
          end if
          dim_trafo = 2*l + 2*s - 2
          ! allocate transformation matrices
          allocate(ylm_proj_trafos(l,s)%matrix(dim_trafo,dim_trafo,group_num_el))
          allocate(comp_to_real(l, s, 1)%matrix(2*l + 1, dim_trafo))
          allocate(comp_to_real(l, s, 2)%matrix(2*l + 1, dim_trafo))
          ! generate transformation matrices
          call efm_proj_trafo_gen(l,s,ylm_proj_trafos(l,s)%matrix,comp_to_real(l,s,2)%matrix,comp_to_real(l,s,1)%matrix)
       end do
    end do

  end subroutine symm_proj_trafo_alloc
  !*************************************************************

  !*************************************************************
  subroutine symm_trafo_dealloc
    !  Purpose: deallocate transformation matrices
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    ! use unique_atom_module !<<<mdf
    use group_module
    use efm_decl
    use efm_module
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: l,s
    ! angular momentum counter

    ! deallocate transformation matrices
    do l=0,ubound(ylm_trafos,1)
       deallocate(ylm_trafos(l)%matrix)
    end do
    deallocate(ylm_trafos)

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do l=0,ubound(ylm_proj_trafos,1)
          do s=1,2
             if (s.eq.1) then
                cycle ! FIXME: what is this?
             endif
             deallocate(ylm_proj_trafos(l,s)%matrix)
             deallocate(comp_to_real(l, s, 1)%matrix)
             deallocate(comp_to_real(l, s, 2)%matrix)
          end do
       end do
       deallocate(ylm_proj_trafos,comp_to_real)
    endif

    ! deallocation of point_trafos is
    ! done by symm_positions::done(), called from main_symm()
  end subroutine symm_trafo_dealloc
  !*************************************************************

  subroutine symm_symmadapt(spor)
    use symm_positions, only: point_trafos_timps,point_trafos_pc1
    use pointcharge_module, only:&
         & n_timps, n_moving_unique_timps, unique_timps
    use pointcharge_module, only: unique_pc, moving_pc, pointcharge_n, pointcharge_array
    implicit none
    logical, intent(in) :: spor
    integer(IK) :: i,status
    ! *** end of interface ***

    call symm_adapt_uatoms(spor)

    if(n_timps.ne.0.and.n_moving_unique_timps.ne.0) then
       call symm_adapt_centers(unique_timps, point_trafos_timps, L_MAX=1)
    endif

    if(pointcharge_n > 0 .and. moving_pc) then
       allocate(unique_pc(pointcharge_n),stat=status)
       ASSERT(status==0)
       do i=1,pointcharge_n
          unique_pc(i)%n_equal_atoms=pointcharge_array(i)%n_equal_charges
          unique_pc(i)%lmax_all=1
          allocate(unique_pc(i)%position(3,pointcharge_array(i)%N_equal_charges), &
               stat=status)
          ASSERT(status==0)
          unique_pc(i)%position=pointcharge_array(i)%position
       end do
       call symm_adapt_centers(unique_pc, point_trafos_pc1, L_MAX=1)
    end if

  end subroutine symm_symmadapt

  subroutine symm_adapt_uatoms(spor)
    use error_module, only: error
    use uatom_symmadapt, only: uaSymm, sa_calc_sa_int
    use unique_atom_module, only: unique_atom_lmax_all, unique_atoms, n_unique_atoms
    use symm_positions, only: point_trafos
    implicit none
    logical,optional,intent(in) :: spor
    ! *** end of interface ***

    logical                 :: spin_orbit
    integer(IK)             :: memstat
    integer(IK)             :: ua
    type(vsalcs_ua),pointer :: vsalcs(:) ! vsalcs(n_ua)

    
    DPRINT 'sm/symm_adapt_uatoms: entered'

    spin_orbit = .false.
    if(present(spor)) spin_orbit = spor

    allocate(vsalcs(n_unique_atoms), STAT=memstat)
    ASSERT(memstat==0)

    DPRINT 'sm/symm_adapt_uatoms: call symm_symmadapt_alt(...)'
    call symm_symmadapt_alt(&
         & unique_atom_lmax_all,&
         & unique_atoms%lmax_all,&
         & ylm_trafos,&
         & n_unique_atoms,&
         & unique_atoms%n_equal_atoms,&
         & point_trafos,&
         & vsalcs)

    if(spin_orbit)then
       DPRINT 'sm/symm_adapt_uatoms: call symm_proj_symmadapt(...)'
       call symm_proj_symmadapt(vsalcs)
    endif

    DPRINT 'sm/symm_adapt_uatoms: call pseudo2d(...)'
    call pseudo2d(&
         & unique_atom_lmax_all,&
         & unique_atoms%lmax_all,&
         & n_unique_atoms,&
         & unique_atoms%n_equal_atoms,&
         & vsalcs)

    DPRINT 'sm/symm_adapt_uatoms: call load_vsalcs_to_uatom(...)'
    call load_salcs(&
         & n_unique_atoms,&
         & unique_atoms%n_equal_atoms,&
         & unique_atoms%lmax_all,&
         & vsalcs,&
         & uaSymm)

    DPRINT 'sm/symm_adapt_uatoms: call sa_calc_sa_int(uaSymm(ua)):'
    do ua = 1, n_unique_atoms
       ! translate symm-adapt coeff to sa_int form:
       call sa_calc_sa_int(uaSymm(ua))
    enddo

    DPRINT 'sm/symm_adapt_uatoms: call free_vsalcs_ua()'
    do ua = 1, size(vsalcs)
       call free_vsalcs_ua(vsalcs(ua))
    enddo
    deallocate(vsalcs,STAT=memstat)
    ASSERT(memstat==0)
    DPRINT 'sm/symm_adapt_uatoms: exit'
  end subroutine symm_adapt_uatoms

  subroutine symm_adapt_centers(ucenters, pnt_trafos, l_max)
    use error_module, only: error
    use uatom_symmadapt
    use unique_atom_module,  only: unique_atom_type
    use unique_atom_methods, only: unique_atom_assign_symm
    implicit none
    type(unique_atom_type), intent(inout)     :: ucenters(:)
    type(symm_transformation_int), intent(in) :: pnt_trafos(:)
    integer(IK), intent(in), optional         :: l_max
    ! *** end of interface ***

    integer(IK)             :: memstat
    integer(IK)             :: ua,lmax,n_ua
    integer(IK),allocatable :: angs(:)
    type(vsalcs_ua),pointer :: vsalcs(:) ! vsalcs(n_ua)
    type(uatom),allocatable :: ucSymm(:)

    DPRINT 'sm/symm_symmadapt_centers: entered'

    n_ua = size(ucenters)
    allocate(vsalcs(n_ua),angs(n_ua),ucSymm(n_ua),STAT=memstat)
    if(memstat/=0) call error("sm/symm_adapt_centers: alloc failed")

    if(present(l_max))then
       angs(:) = l_max
    else
       do ua=1,n_ua
          angs(ua) = ucenters(ua)%lmax_ob
       enddo
    endif
    lmax = maxval(angs(:))

    DPRINT 'sm/symm_adapt_centers: call symm_symmadapt_alt(...)'
    call symm_symmadapt_alt(&
         & lmax,&
         & angs,&
         & ylm_trafos,&
         & n_ua,&
         & ucenters%n_equal_atoms,&
         & pnt_trafos,&
         & vsalcs)

    DPRINT 'sm/symm_adapt_centers: call pseudo2d(...)'
    call pseudo2d(&
         & lmax,&
         & angs,&
         & n_ua,&
         & ucenters%n_equal_atoms,&
         & vsalcs)

    DPRINT 'sm/symm_adapt_centers: call load_vsalcs_to_uatom(...)'
    call load_salcs(&
         & n_ua,&
         & ucenters%n_equal_atoms,&
         & angs,&
         & vsalcs,&
         & ucSymm)

    DPRINT 'sm/symm_adapt_centers: call sa_calc_sa_int(uaSymm(ua)):'
    do ua = 1, n_ua
       ! translate symm-adapt coeff to sa_int form:
       call sa_calc_sa_int(ucSymm(ua))

       ! copy uatom components to unique_atom_components
       call unique_atom_assign_symm(ucenters(ua), ucSymm(ua))
    enddo

    deallocate(vsalcs,angs,ucSymm,STAT=memstat)
    if(memstat/=0) call error("symm_symmadapt: dealloc failed")

    DPRINT 'sm/symm_adapt_centers: exit'
  end subroutine symm_adapt_centers


  !*************************************************************
  subroutine symm_symmadapt_alt(lmax,angs,ylm_trafos,n_ua,nea,pnt_trafos,salcs)
    !  Purpose: symmetryadaption of basisfunctions on
    !           all symmetry equivalent atoms and for all
    !           angular momenta
    !           THUS: this is the central routine of the 
    !                 Symmetry part!
    !------------ Modules used -----------------------------------
    use iounitadmin_module
    use group_module, only: n_irr => group_num_ir, n_el => group_num_el
    use efm_decl
    use efm_module
    use clebsch_gordan
    use unique_atom_module
    implicit none
    integer(IK),intent(in)                   :: lmax,angs(:),n_ua,nea(:) ! (n_ua)
    type(symm_transformation),intent(in)     :: ylm_trafos(0:) ! ylm_trafos(0:lmax)
    type(symm_transformation_int),intent(in) :: pnt_trafos(:)  ! pnt_trafos(n_ua)
    type(vsalcs_ua),intent(inout)            :: salcs(:)       ! salcs(n_ua)
    !** End of interface *****************************************
    !------------ Declaration of local constants  ----------------
    real(RK),parameter :: small = 1.e-10 
    ! very small value
    !------------ Declaration of local variables ---------------------
    integer(IK)                            ::&
         & i,j,k,&
         & i_m,j_m,k_m,&
         & i_p,j_p,k_p,&
         & i_dim,j_dim,k_dim,&
         & n,m
    ! counters
    integer(IK)                            :: ua,n_ea,ang,n_lm
    ! unique atom and angular momentum
    integer(IK),allocatable                :: multipl(:),multipl_total(:)
    ! multiplicity of irreps
    real(RK),pointer                       :: vector_point(:),vector_ylm(:)
    ! work array
    real(RK),allocatable                   :: trafo_matrices(:,:,:), dproduct(:,:)
    ! trafo_matrices
    integer(IK)                            :: mult_point,mult_ylm,mult_irrep
    ! multiplicties of irreps 

    type(efm_cscoII_diag_invspace) :: ylm_adapt(0:lmax)
    ! symmetryadapted spherical harmonics
    type(efm_cscoII_diag_invspace) :: point_adapt(n_ua)
    ! symmetryadapted points

    type(sym_prod),pointer  :: salc(:)
    type(sym_prod),pointer  :: eigs
    type(prod_bas), pointer :: pb

    integer(IK) :: alloc_stat

    logical :: io
    !------------ Declaration of subroutines used ----------------
    external error_handler

    !------------ Executable code ------------------------------------

    DPRINT 'sm/symm_symmadapt_alt: entered'

    ! Slaves do not always have output_unit open:
    io = output_unit > 0 .and. .not. no_output_unit_output

    ! allocate work arrays
    allocate(multipl(n_irr),multipl_total(n_irr),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symm_symmadapt_alt: allocation of multiplicity failed")

    DPRINT 'sm/symm_symmadapt_alt: symmetrize spherical harmonics'
    !
    ! symmetrize spherical harmonics
    !
100 FORMAT('*', 2X,99A)
101 FORMAT('*', 20(A6))
102 FORMAT('*', 20(I6))

    if (io) then
       write(output_unit,100) ('*',i=1,6*(n_irr+2)-2)
       write(output_unit,100) "1. SYMMETRIZING SPHERICAL HARMONICS,"
       write(output_unit,100) "   IRREP MULTIPLICITIES:"
       write(output_unit,101) "L","2L+1",(trim(irrep_can(i)%label),i=1,n_irr)
    endif

    ! loop over angular momenta
    do ang=0,lmax

       n_lm = 2*ang + 1       
       call efm_alloc(n_irr,ylm_adapt(ang),dim=n_lm)
       ! symmetrize
       call efm_cscoII_solve(ylm_trafos(ang)%matrix,n_lm,ylm_adapt(ang))

       if (io) then
          write(output_unit,102) ang,n_lm,multiplicities(ylm_adapt(ang)%csco_eigenspaces)
       endif
    end do

    DPRINT 'sm/symm_symmadapt_alt: symmetrize nuclear positions'
    ! 
    ! symmetrize nuclear positions
    !

    if (io) then
       write(output_unit,100) ('*',i=1,6*(n_irr+2)-2)
       write(output_unit,100) "2. SYMMETRIZING ATOMIC POSITIONS,"
       write(output_unit,100) "   IRREP MULTIPLICITIES:"
       write(output_unit,101) "UA","NEA",(trim(irrep_can(i)%label),i=1,n_irr)
    endif

    ! loop over all unique atoms
    do ua=1,n_ua
       n_ea = nea(ua) !unique_atoms(unique)%N_equal_atoms

       ! allocate real transformation matrices
       allocate(trafo_matrices(n_ea,n_ea,n_el),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("symm_symmadapt_alt: allocation of trafo_matrices failed")

       call efm_alloc(n_irr,point_adapt(ua),dim=n_ea)

       trafo_matrices = REAL(pnt_trafos(ua)%matrix,KIND=RK)
       ! symmetrize
       call efm_cscoII_solve(trafo_matrices, n_ea, point_adapt(ua))

       if (io) then
          write(output_unit,102) ua,n_ea,multiplicities(point_adapt(ua)%csco_eigenspaces)
       endif

       ! deallocate real transformation matrices
       deallocate(trafo_matrices)
    end do

    DPRINT 'sm/symm_symmadapt_alt: now alloc SALCS'
    !
    ! now construct SALCS
    !
    do ua=1,n_ua
       ang  = angs(ua)
       n_ea = nea(ua)
       call alloc_vsalcs_ua(ang, n_ea, n_irr, salcs(ua))
    enddo
    
    DPRINT 'sm/symm_symmadapt_alt: now construct SALCS'

    if (io) then
       write(output_unit,100) ('*',i=1,6*(n_irr+2)-2)
       write(output_unit,100) "3. SYMMETRIZING MOLECULAR SHELLS,"
       write(output_unit,100) "   IRREP MULTIPLICITIES, AND THE SUM:"
       write(output_unit,101) "UA","L",(trim(irrep_can(i)%label),i=1,n_irr)
    endif

    multipl_total = 0
    ! loop over all unique atoms
    ua_: do ua=1,n_ua
       ! loop over angular momenta
       ang_: do ang=0,angs(ua)

!!$          DPRINT 'sm/symm_symmadapt_alt: ua=',ua,' L=',ang

          salc => salcs(ua)%mos(:,ang)

          ! dimension of the direct product space
          n_lm = 2*ang + 1
          n_ea = nea(ua)

          !
          ! determine multiplicity of the irreps
          !
!!$          DPRINT 'sm/symm_symmadapt_alt: determine multiplicity:'
          multipl = 0

          ! loop over all combinations of irreps of point trafos and 
          ! spherical harmonics
          do k=1,n_irr
             k_dim      = irrep_can(k)%dimension

             do i=1,n_irr
                do j=1,n_irr

                   mult_point = multiplicity(point_adapt(ua)%csco_eigenspaces(i))
                   mult_ylm   = multiplicity(ylm_adapt(ang)%csco_eigenspaces(j))

                   mult_irrep = cg(k,i,j)%mult

                   if(mult_irrep/=multiplicity(efm_cg(i,j)%csco_eigenspaces(k)))&
                        & call error("sm/symm_symmadapt_alt: multiplicity?")

                   mult_irrep = mult_irrep*mult_point*mult_ylm

                   multipl(k) = multipl(k) + mult_irrep
                enddo
             enddo

             call cg_alloc(multipl(k),salc(k))
             do k_m=1,multipl(k)
                call cg_alloc(&
                     & k_dim,&
                     & n_ea,&
                     & n_lm,&
                     & salc(k)%sub(k_m)&
                     & ) 
             enddo
          enddo

          ! increment total multiplicity
          multipl_total = multipl_total + multipl

!!$          DPRINT 'sm/symm_symmadapt_alt: Multiplicities=',multipl_total

          ! show multiplicities
          if (io) then
             write(output_unit,102) ua,ang,(multipl(k),k=1,n_irr)
          endif

          ! allocate work arrays
          allocate(dproduct(n_ea,n_lm),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("symm_symmadapt_alt: allocation of work arrays failed")

          !
          ! calculate SALCS for current unique atom and angular momentum
          !

          ! init:
          do k=1,n_irr
             do k_m=1,salc(k)%mult
                salc(k)%sub(k_m)%c = 0.0_rk
             enddo
          enddo
          multipl(:) = 0 ! offset now

          ! loop over all combinations of irreps of point trafos and 
          ! spherical harmonics
          
          i_:do i=1,n_irr
             i_dim      = irrep_can(i)%dimension
             mult_point = multiplicity(point_adapt(ua)%csco_eigenspaces(i))

             j_:do j=1,n_irr
                j_dim    = irrep_can(j)%dimension
                mult_ylm = multiplicity(ylm_adapt(ang)%csco_eigenspaces(j))

                
                ! loop over all combinations of multiplicities for the
                ! current irrep combination
                i_m_:do i_m=1,mult_point
                   j_m_:do j_m=1,mult_ylm
                      ! --- begin of loops over combinations of multiplicities 

                      i_p_:do i_p=1,i_dim
                         j_p_:do j_p=1,j_dim

                            ! we consider the product space of:
                            ! point: (irrep i,i_m) x ylm: (irrep j,j_m)
                            vector_point => point_adapt(ua)&
                                 & %csco_eigenspaces(i)&
                                 & %cscoII_subspaces(i_p)&
                                 & %eigenfunction(:,i_m)
                            vector_ylm => ylm_adapt(ang)&
                                 & %csco_eigenspaces(j)&
                                 & %cscoII_subspaces(j_p)&
                                 & %eigenfunction(:,j_m)
                            do n=1,n_ea
                               do m=1,n_lm
                                  dproduct(n,m) = vector_point(n) * vector_ylm(m)
                               enddo
                            enddo

                            ! now symmetry adaption of the product space is performed
                            ! employing cg-coefficients
                            !

                            ! loop over product irreps:
                            k_:do k=1,n_irr
                               k_dim    = irrep_can(k)%dimension
                               eigs     => salc(k)
                               ! loop over all partners of the irrep
                               k_p_:do k_p=1,k_dim
                                  ! loop over multiplicity of product:
                                  k_m_:do k_m=1,cg(k,i,j)%mult
                                     pb  => eigs%sub(multipl(k)+k_m) ! Product Basis
                                     !
                                     ! this is where the SALCS are assembled!
                                     !
                                     pb%c(k_p,1:n_ea,1:n_lm) = pb%c(k_p,1:n_ea,1:n_lm)&
                                          & + cg(k,i,j)%sub(k_m)%c(k_p,i_p,j_p)&
                                          & * dproduct
                                  enddo k_m_
                               enddo k_p_
                            enddo k_

                         end do j_p_
                      end do i_p_

                      ! shift to account for calculated subspaces:
                      do k=1,n_irr
                         multipl(k) = multipl(k) + cg(k,i,j)%mult
                      enddo
                      ! --- end of loops over combinations of multiplicities 
                   enddo j_m_
                end do i_m_
             end do j_
             ! end of loops over combinations of irreps
          end do i_

       ! deallocate work array
       deallocate(dproduct,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("symm_symmadapt_alt: deallocation of work arrays failed")

       end do ang_
    end do ua_

    ! free efm_ structs:
    do ua=1,n_ua
       call efm_free(point_adapt(ua))
    enddo
    do ang=0,lmax
       call efm_free(ylm_adapt(ang))
    enddo

    ! show multiplicities
    if (io) then
       write(output_unit,102) -1,-1,(multipl_total(k),k=1,n_irr)
       write(output_unit,100) ('*',i=1,6*(n_irr+2)-2)
    endif

    ! deallocate work arrays
    deallocate(multipl,multipl_total)

    DPRINT 'sm/symm_symmadapt_alt: exit'
  end subroutine symm_symmadapt_alt

  subroutine pseudo2d(lmax, angs, n_ua, nea, salcs)
    use group_module, only: n_irr => group_num_ir
    use clebsch_gordan
    implicit none
    integer(IK), intent(in) :: lmax, angs(:), n_ua, nea(:) ! (n_ua)
    type(vsalcs_ua), intent(inout) :: salcs(:) ! salcs(n_ua)
    !** End of interface *****************************************

    integer(IK) :: irr, k_m, k_p
    integer(IK) :: k_m1D ! counters
    integer(IK) :: ua, ang ! unique atom and angular momentum
    integer(IK) :: mult, na, nb, nc

    type(sym_prod) :: bas1D

    !
    ! special treatment of pseudo-2D irreps
    !

    ! pseudo-2D irreps were considered as two dimensional irreps.
    ! But since the grand orthogonality theorem doesn`t apply on them
    ! we actually must consider them as one dimensional, i.e. in the
    ! Hamiltonian matrix the rows can interact. Thus for use in the
    ! Clustercode the datastructure must be changed.
    DPRINT 'sm/pseudo2d: entered'

    ! loop over all irreps
    do irr = 1, n_irr
          if ( .not. irrep_can(irr)%pseudo ) cycle ! irrep
          DPRINT 'sm/pseudo2d: irrep=', irr ,' is a pseudo irrep'
          do ua = 1, n_ua
             do ang = 0, angs(ua)
                DPRINT 'sm/pseudo2d: ua=', ua, ' L=', ang

                !
                ! Extract dimensions from the irreducible eigenspace
                ! of a (pseudo) irrep "irr" originating from the shell
                ! "ang" of unique atom group "ua".  Dont trust them
                ! yet if multiplicity is zero:
                !
                call clebsch_gordan_dimensions(salcs(ua)%mos(irr, ang), mult, na, nb, nc)

                DPRINT 'sm/pseudo2d: mult=', mult * na, "=", mult, "*", na
                if ( mult /= 0 ) then
                   ASSERT(na==2)
                endif

                !
                ! Two partners of the pseudo-2D will become
                ! independent instances:
                !
                bas1D = clebsch_gordan_new(mult * na, 1, nb, nc)

                !
                ! Now copy the data into the fresh allocated sotrage,
                ! adjusting for the new shape of the structure:
                !
                k_m1D = 0
                do k_p = 1, 2 ! na, actually
                   do k_m = 1, mult
                      k_m1D = k_m1D + 1
                      DPRINT 'sm/pseudo2d: k_m=',k_m,' k_p=',k_p,' // k_m1D=',k_m1D

                      !
                      ! Not the other way around:
                      !
                      !   do k_m = 1, mult
                      !     do k_p = 1, 2
                      !
                      ! AM, VN 07.10.2001.
                      !

                      bas1D%sub(k_m1D)%c(1, :, :) = salcs(ua)%mos(irr, ang)%sub(k_m)%c(k_p, :, :)
                   enddo
                enddo

                !
                ! Replace the item in the inten(input) argument array:
                !
                DPRINT 'sm/pseudo2d: replace'

!!$                ! FIXME: is deallocation necessary if we do not use
!!$                !        pointers but allocatable components?
!!$                call cg_free(salcs(ua)%mos(irr, ang))

                !
                ! This is the irreducible eigenspace of a (pseudo)
                ! irrep "irr" originating from the shell "ang" of
                ! unique atom group "ua". Now not "degenerate"
                ! anymore.
                !
                ! Gfortran 4.3 (Debian Lenny) has troubles
                ! transferring allocation status of zero-sized
                ! allocatable components. This should be a user
                ! defined assignemnt(=) in such case:
                !
                salcs(ua)%mos(irr, ang) = bas1D
             enddo ! angular momenta
          enddo ! unique atoms
       enddo ! irreps
    DPRINT 'sm/pseudo2d: exit'
  end subroutine pseudo2d

  subroutine load_vsalcs_to_uatom(n_ua,nea,angs,salcs,ua)
    ! 
    ! write SALCS to datastructures of unique_atom_module
    !
    use group_module, only:&
         & n_irr=>group_num_ir,irrep_can
    use efm_decl
    use clebsch_gordan
    use uatom_symmadapt
    use matrix_sparse,&
         & sparse_alloc => alloc,&
         & sparse_free  => free,&
         & sparse_pack  => spack
    implicit none
    integer(IK),intent(in)                     :: n_ua
    integer(IK),intent(in)                     :: angs(:),nea(:)
    type(vsalcs_ua),intent(in)                 :: salcs(:)
    type(uatom),intent(inout)                  :: ua(:)
    target :: ua
    ! *** end of intarface ***

    integer(IK) :: k,unique,ang,n_equiv,irrep_dim,irrep_mult

    integer(IK) :: k_m,k_p
    integer(IK) :: lmax,nm,m_ylm,m_clm

    type(sym_prod),pointer  :: salc(:)
    type(partner_type),pointer  :: sapartn
    type(SparseRMatrix)         :: c_pck

    DPRINT 'sm/load_vsalcs_to_uatom: entered'

      ! loop over all unique atoms
    do unique=1,n_ua
       lmax = angs(unique)
       n_equiv = nea(unique)

!!$       DPRINT 'sm/load_vsalcs_to_uatom: unique=',unique,' lmax=',lmax

       call sa_alloc(n_irr, lmax, ua(unique))
       
       ! loop over angular momenta
       do ang=0,lmax
!!$          DPRINT 'sm/load_vsalcs_to_uatom: ang=',ang
          salc => salcs(unique)%mos(:,ang)
             
          ! loop over all irreps
          k_:do k=1,n_irr
             irrep_dim  = irrep_can(k)%time_dimension
             irrep_mult = salc(k)%mult

             sapartn => ua(unique)%symadapt_partner(k,ang)

             call sa_alloc(n_equiv,irrep_mult,irrep_dim,sapartn)

             k_m_:do k_m=1,irrep_mult

                ! over partner within the irrep:
                k_p_:do k_p=1,irrep_dim

                   call sparse_pack(&
                        & salc(k)%sub(k_m)%c(k_p,:,:),&
                        & c_pck&
                        & )

                   ! change convention:
                   ! m_Ylm - convention of symmetry part
                   ! m-Clm - convention of the rest of the code
                   do nm=1,c_pck%n
                      m_ylm = c_pck%j(nm)
                      m_clm = m_ylm - ang - 1
                      if(m_clm<0)then
                         m_clm = - 2 * m_clm
                      else
                         m_clm =   2 * m_clm + 1
                      endif
                      c_pck%j(nm) = m_clm
                   enddo

                   call sa_alloc(c_pck%n,sapartn%symadapt(k_m,k_p))

                   ! fill in new datastructure
                   sapartn%symadapt(k_m,k_p)%I_equal_atom = c_pck%i
                   sapartn%symadapt(k_m,k_p)%m            = c_pck%j
                   sapartn%symadapt(k_m,k_p)%c            = c_pck%c

                   call sparse_free(c_pck)
                enddo k_p_
             enddo k_m_
          enddo k_
  
       end do! loop over angular momenta
    end do! loop over all unique atoms
  end subroutine load_vsalcs_to_uatom

  subroutine symm_proj_symmadapt(vsalcs)
    use uatom_symmadapt, only: uatom, uatom_symmadapt_init
    use unique_atom_module
    use group_module, only:&
         & irrep_can,&
         & proj_irrep_can,&
         & n_preps=>group_num_re
    use efm_module, only:&
         & efm_proj_cg
    use error_module, only:&
         & error
    use spin_orbit_module, only:&
         & is_on,&
         & op_FitTrafo
    implicit none
    type(vsalcs_ua),intent(in) :: vsalcs(:)
    ! *** end of interface ***

    real(RK),parameter             :: zero = 0.0_rk
    type(uatom),pointer            :: uaL(:),uaS(:)
    type(psalcs_ua),pointer        :: salcsL(:),salcsS(:)

    type(symm_basis_transformation_c),pointer :: us(:)
    ! U-matriceS(n-proj-irreps)
    integer  :: i,memstat

    !---------------------------------
    ! preparations:
    !

    allocate(&
         & salcsL(n_unique_atoms), salcsS(n_unique_atoms),&
         & uaL(n_unique_atoms),   uaS(n_unique_atoms),&
         & us(n_preps),&
         & STAT=memstat)
    call error(memstat,"sm/symm_proj_symmadapt: (wrapper) alloc failed")

    !---------------------------------
    ! now actual symmetry adaption:
    !

    call proj_symmadapt_su2(vsalcs,salcsL)

    DPRINT 'sm/symm_proj_symmadapt: LOADING uas_... Large'

    call load_salcs(&
         & n_unique_atoms, unique_atoms%n_equal_atoms, unique_atoms%lmax_all,&
         & salcsL,         uaL )

    call uatom_symmadapt_init(Large=uaL)

!!$    if(is_on(op_FitTrafo))then

    ! PSEUDO-SCALAR COUPLING:
    call alloc_arr_trafo_c(n_preps, proj_irrep_can%dimension, us)

    call pseudo_scalar_coupling(irrep_can,proj_irrep_can,efm_proj_cg,pcoupling,us)

    DPRINT 'sm/symm_proj_symmadapt: PCOUPLING=',(i,i=1,size(pcoupling))
    DPRINT '----------------------------------',pcoupling

    call rephase(pcoupling,us)

    if(is_on(op_FitTrafo))then
       !
       ! Small Component Symmetryzation has to be done
       !
       call large_to_small_salcs(pcoupling, us, salcsL, salcsS)

!!$       call free(us)

       DPRINT 'sm/symm_proj_symmadapt: LOADING uas_... Small'
       call load_salcs(&
            & n_unique_atoms, unique_atoms%n_equal_atoms, unique_atoms%lmax_all,&
            & salcsS,         uaS )

       call uatom_symmadapt_init(Small=uaS)
       
       do i=1,size(salcsS)
          call free_psalcs_ua(salcsS(i))
       enddo
    endif

    ! Complex-Conjugate (Kramers) Coupling:
    call cc_coupling(proj_irrep_can,cccoupling,us)

    DPRINT 'sm/symm_proj_symmadapt: CCCOUPLING=',(i,i=1,size(cccoupling))
    DPRINT '-----------------------------------',cccoupling


    !---------------------------------
    ! clean up:
    !
    do i=1,size(salcsL)
       call free_psalcs_ua(salcsL(i))
    enddo
    call free_arr_trafo_c(us)
    deallocate(salcsL,salcsS,uaL,uaS,us,STAT=memstat)
    call error(memstat,"sm/symm_proj_symmadapt: (wrapper) dealloc failed")
  end subroutine symm_proj_symmadapt

  subroutine proj_symmadapt_su2(vsalcs,psalcs)
    use group_module, only:&
         & n_pirr=>group_num_re, n_virr=>group_num_ir &
         , proj_irrep_can
    use clebsch_gordan
    use error_module
    use iounitadmin_module, only: output_unit
    implicit none
    type(vsalcs_ua),intent(in)    ::  vsalcs(:)
    type(psalcs_ua),intent(inout) ::  psalcs(:)
    ! *** end of interface ***

    integer(IK) :: lmax,n_ua,n_ea,l,j
    integer(IK) :: ua

    if(size(vsalcs).ne.size(psalcs))&
         & call error("sm/proj_symmadapt_su2: sizes?")

    n_ua = size(vsalcs)

    do ua = 1, n_ua
       lmax = vsalcs(ua)%lmax
       n_ea = vsalcs(ua)%n_ea
       call alloc_psalcs_ua(lmax, n_ea, n_pirr, psalcs(ua))
    enddo

100 FORMAT('*', 2X,99A)
101 FORMAT('*', 20(A6))
102 FORMAT('*', 20(I6))

    if (output_unit > 0) then
       write(output_unit,100) ('*',j=1,6*(n_pirr+2)-2)
       write(output_unit,100) "1. SYMMETRIZING MOLECULAR SPINORS,"
       write(output_unit,100) "   IRREP MULTIPLICITIES:"
       write(output_unit,101) "UA","L",(trim(proj_irrep_can(j)%label),j=1,n_pirr)
    endif

    do ua = 1, n_ua
       lmax = vsalcs(ua)%lmax
       do l = 0, lmax
          do j = 1, n_pirr
             DPRINT 'sm/proj_symmadapt_su2: ua=', ua, ' L=', l, ' pirr=', j 
             call xsu2(vsalcs(ua)%mos(:, l), psalcs(ua)%mos(j, l, :), vsu2cg(j, :))
             DPRINT 'sm/proj_symmadapt_su2: MULT=', psalcs(ua)%mos(j, l, 1)%mult
          enddo
          if (output_unit > 0) then
             write(output_unit,102) ua,l, (psalcs(ua)%mos(j,l,1)%mult,j=1,n_pirr)
          endif
       enddo
    enddo
  contains 
    subroutine xsu2(vs, ps, cg)
      implicit none
      type(sym_prod), intent(in) ::  vs(:) ! vs(n_virr)
      type(sym_prod), intent(inout) ::  ps(:) ! ps(2)
      type(sym_prod), intent(in) ::  cg(:) ! cg(n_virr)
      ! *** end of interface ***

      integer(IK) :: m_cg,m_vs,i,c,&
           & k_ps,k_vs,k_cg,&
           & n_vp,n_ea,n_lm,n_pp,&
           & p_p,p_v
      integer(IK) :: mult

      ASSERT(n_virr==size(vs))
      ASSERT(2==size(ps))

      ! count subspaces:
      mult = 0
      do i=1,n_virr
         DPRINT 'sm/xsu2: virr=',i
         m_vs = vs(i)%mult
         m_cg = cg(i)%mult
         mult = mult +  m_vs * m_cg
         DPRINT 'sm/xsu2: +(',m_vs,' * ',m_cg,')=',m_vs * m_cg
!!$         if(m_cg>1) call error("sm/xsu2: evrika !")
         if(m_cg>1)then
            WARN("VIrrep*SU(2) -> PIrrep x 2")
         endif
      enddo
      DPRINT 'sm/xsu2: mult=',mult

      ! pre allocate projective salcs:
      do c = 1, 2
         call cg_alloc(mult, ps(c))
      enddo

      ! fill in projective salcs
      k_ps = 0
      irrep: do i = 1, n_virr
         m_vs = vs(i)%mult
         m_cg = cg(i)%mult

         vec_mult: do k_vs = 1, m_vs
            cg_mult: do k_cg = 1, m_cg
               ! proj_mult: do k_ps=...
               k_ps = k_ps + 1
               component: do c = 1, 2

                  n_vp = size(vs(i)%sub(k_vs)%c, 1) ! number of vector partners
                  n_ea = size(vs(i)%sub(k_vs)%c, 2) ! number of equivalent atoms
                  n_lm = size(vs(i)%sub(k_vs)%c, 3) ! number of LM harmonics
                  n_pp = size(cg(i)%sub(k_cg)%z, 1) ! number of projective partners

                  call cg_alloc(n_pp, n_ea, n_lm, ps(c)%sub(k_ps), cmplx=.true.)

                  !
                  !    mu           mu
                  ! PSI   = SUM   CG       * PHI
                  !    pp      pv   pp, pv      pv
                  !
                  ps(c)%sub(k_ps)%z(:, :, :) = 0.0_rk
                  proj_partner: do p_p = 1, n_pp
                     vec_partner: do p_v = 1, n_vp
                        ps(c)%sub(k_ps)%z(p_p, :, :) = &
                            ps(c)%sub(k_ps)%z(p_p, :, :) + &
                            cg(i)%sub(k_cg)%z(p_p, p_v, c) * vs(i)%sub(k_vs)%c(p_v, :, :)
                     enddo vec_partner
                  enddo proj_partner
                  ps(c)%sub(k_ps)%re =  REAL(ps(c)%sub(k_ps)%z, kind=RK)
                  ps(c)%sub(k_ps)%im = AIMAG(ps(c)%sub(k_ps)%z)

               enddo component
            enddo cg_mult
         enddo vec_mult
      enddo irrep
    end subroutine xsu2
  end subroutine proj_symmadapt_su2

  subroutine load_psalcs_to_uatom(n_ua,nea,angs,salcs,uat)
    ! 
    ! write SALCS to datastructures of unique_atom_module
    !
    use group_module, only: n_irr=>group_num_re, proj_irrep_can
    use clebsch_gordan
    use uatom_symmadapt
    use matrix_sparse,&
         & sparse_alloc => alloc,&
         & sparse_free  => free,&
         & sparse_pack  => spack
    implicit none
    integer(IK),intent(in)                     :: n_ua
    integer(IK),intent(in)                     :: angs(:),nea(:)
    type(psalcs_ua),intent(in)                 :: salcs(:)
    type(uatom),intent(inout)                  :: uat(:)
    target :: uat
    ! *** end of intarface ***

    integer(IK) :: k,ua,ang,n_ea,irrep_dim,irrep_mult,n_lm_fcts

    integer(IK) :: k_m,k_p,c,ea
    integer(IK) :: lmax,nm,m_ylm,m_clm

    type(sym_prod),pointer      :: salc(:)
    type(partner_type),pointer  :: sapartn
    type(sa_int_type),pointer   :: sat(:)
    type(SparseCMatrix)         :: c_pck

    DPRINT 'sm/load_psalcs_to_uatom: entered'

      ! loop over all unique atoms
    do ua=1,n_ua
       lmax = angs(ua)
       n_ea = nea(ua)

!!$       DPRINT 'sm/load_psalcs_to_uatom: unique=',ua,' lmax=',lmax

       call sa_alloc(n_irr,lmax,uat(ua),SPOR=.true.)
       
       ! loop over angular momenta
       do ang=0,lmax
!!$          DPRINT 'sm/load_psalcs_to_uatom: ang=',ang
 
          ! loop over all irreps
          k_:do k=1,n_irr
             irrep_dim  = proj_irrep_can(k)%dimension
             salc => salcs(ua)%mos(k,ang,:)
             irrep_mult = salc(1)%mult ! == ...(2)%mult

             sapartn => uat(ua)%symadapt_spor_partner(k,ang)

!!$             DPRINT 'sm/load_psalcs_to_uatom: irr=',k,' (mult ',irrep_mult,')'
             call sa_alloc(n_ea,irrep_mult,irrep_dim,sapartn,SPOR=.true.)

             k_m_:do k_m=1,irrep_mult

                ! over partner within the irrep:
                k_p_:do k_p=1,irrep_dim

                   component: do c=1,2
                      sat  => sapartn%sa_spor_int(:,c,k_m,k_p) !..(n_ea)

!!$                      DPRINT 'sm/load_psalcs_to_uatom: k_m=',k_m,' k_p=',k_p,' c=',c,' (mult ',salc(c)%mult,')'
!!$
!!$                      DPRINT 'sm/load_psalcs_to_uatom: shape(salc%sub(k_m)%z(:,:,:))=',&
!!$                           & shape(salc(c)%sub(k_m)%z(:,:,:))

                      call sparse_pack(&
                           & salc(c)%sub(k_m)%z(k_p,:,:),&
                           & c_pck&
                           & )

                      ! change convention:
                      ! m_Ylm - convention of symmetry part
                      ! m-Clm - convention of the rest of the code
                      do nm=1,c_pck%n
                         m_ylm = c_pck%j(nm)
                         m_clm = m_ylm - ang - 1
                         if(m_clm<0)then
                            m_clm = - 2 * m_clm
                         else
                            m_clm =   2 * m_clm + 1
                         endif
                         c_pck%j(nm) = m_clm
                      enddo
                      ! c_pck%I -- EA-index
                      ! c_pck%J -- LM-index
                      ! c_pck%Z -- coeff

                      ! fill in uatom data structure:
                      do ea=1,n_ea
                         n_lm_fcts = count(c_pck%I.eq.ea)
                         call sa_alloc(n_lm_fcts,sat(ea),cmplx=.true.)

                         if(n_lm_fcts.eq.0) cycle
                         sat(ea)%M  = pack(c_pck%J, (c_pck%I.eq.ea))
                         sat(ea)%RE = pack(c_pck%RE,(c_pck%I.eq.ea))
                         sat(ea)%IM = pack(c_pck%IM,(c_pck%I.eq.ea))
!!$                         j = 0
!!$                         do i=1,c_pck%N
!!$                            if(c_pck%I(i).eq.ea)then
!!$                               j = j + 1
!!$                               sat(ea)%M(j)  = c_pck%J(i)
!!$                               sat(ea)%RE(j) = c_pck%RE(i)
!!$                               sat(ea)%IM(j) = c_pck%IM(i)
!!$                            endif
!!$                         enddo
                      enddo

                      call sparse_free(c_pck)
                   enddo component
                enddo k_p_
             enddo k_m_
          enddo k_
  
       end do! loop over angular momenta
    end do! loop over all unique atoms
  end subroutine load_psalcs_to_uatom

  subroutine pseudo_scalar_coupling(virreps,pirreps,vp_cg,coupl,us)
    use group_module, only: nvirr => group_num_ir,&
         &                  npirr => group_num_re,&
         &                  nel   => group_num_el
    use efm_decl
    use efm_module
    implicit none
    type(group_irrep),                intent(in)  :: virreps(:) !(nvirr)
    type(group_proj_irrep),           intent(in)  :: pirreps(:) !(npirr)
    type(efm_cscoII_diag_invspace),   intent(in)  :: vp_cg(:,:) !(nvirr,npirr)
    integer(IK),                      intent(out) :: coupl(:)   !(npirr)
    type(symm_basis_transformation_c),intent(inout) :: us(:)      !(npirr)
    ! *** end of interface ***

    real(RK),parameter                     :: one=1.0_rk
    integer(IK)                            :: i,j,jp,dim
    real(RK),dimension(:,:,:),allocatable  :: trafos
    type(efm_cscoII_diag_invspace)         :: sa

    integer(IK) :: memstat

    if(size(coupl)   /= npirr)call error_handler("sm/psclr...: coupl   wrong")
    if(size(virreps) /= nvirr)call error_handler('sm/psclr...: virreps wrong')
    if(size(pirreps) /= npirr)call error_handler('sm/psclr...: pirreps wrong')
    if(    size(vp_cg,1) /= nvirr .OR.&
         & size(vp_cg,2) /= npirr)call error_handler('sm/psclr...: vp_cg wrong')

    allocate(trafos(1,1,nel),STAT=memstat)
    if(memstat/=0)call error_handler("sm/psclr...: alloc 0 failed")

    call efm_alloc(nvirr,sa,dim=1)

    trafos = one                ! trivial representation
    call efm_gauge(trafos,-one) ! representation of pscalar
    call efm_cscoII_solve(trafos,1,sa)

    jp = efm_find_one(sa)
    
    DPRINT 'PSEUDOSCALAR TRANSFORMS AS ',virreps(jp)%label

    deallocate(trafos, STAT=memstat)
    if(memstat/=0)call error_handler("sm/psclr...: dealloc 0 failed")
    call efm_free(sa)
    

    DWRITE(*,'(A)') '---  THUS  IT  COUPLES   ---'
    do i=1, npirr

       dim = pirreps(i)%dimension

       j = efm_find_one(vp_cg(jp,i))

       coupl(i)=j

       DWRITE(*,'(A8," and ",A8)') pirreps(i)%label,pirreps(j)%label

       if(pirreps(i)%dimension /= pirreps(j)%dimension)&
            & call error_handler("sm/psclr...: dims of i and j conflict")

       call extract_basis(us(i)%matrix, vp_cg(jp, i)%csco_eigenspaces(j))

    enddo

    if(any( coupl(coupl) /= (/ (j,j=1,npirr) /) ))&
         & call error_handler("sm/psclr...: strange type of coupling (2)")

    DWRITE(*,'(A)') '----------------------------'
  end subroutine pseudo_scalar_coupling

  subroutine cc_coupling(pirreps,coupl,us)
    use group_module, only: npirr => group_num_re
    use efm_decl
    use efm_module
    implicit none
    type(group_proj_irrep),           intent(in)  :: pirreps(:) !(npirr)
    integer(IK),                      intent(out) :: coupl(:)   !(npirr)
    type(symm_basis_transformation_c),intent(inout) :: us(:)      !(npirr)
    ! *** end of interface ***

    integer(IK)                            :: i,j,dim
    type(efm_cscoII_diag_invspace)         :: sa

    ASSERT(size(coupl)==npirr)
    ASSERT(size(pirreps)==npirr)


    DWRITE(*,'(A)') '---  KRAMERS CONJUGATION COUPLES  ---'

    do i=1, npirr
       DPRINT 'ccc: irr1=',i
       dim = pirreps(i)%dimension

       call efm_alloc(npirr,sa,proj=.true.,dim=dim)

       DPRINT 'ccc: call efm_cscoII_solve()'
       call efm_proj_cscoII_solve&
            & (&
            & conjg(proj_irrep_can(i)%irrep_matrix),&
            & dim,&
            & sa &
            & )

       j = efm_find_one(sa)
       DPRINT 'ccc: irr2=',j
       coupl(i)=j
       DWRITE(*,'(A8," and ",A8)') pirreps(i)%label,pirreps(j)%label

       if(pirreps(i)%dimension /= pirreps(j)%dimension)&
            & call error_handler("sm/ccc...: dims of i and j conflict")

       call extract_basis(us(i)%matrix, sa%csco_eigenspaces(j))

       DPRINT '--- KRAMERS COUPLING ---'
       DCALL show_fine(us(i)%matrix)

       call efm_free(sa)
    enddo

    if(any( coupl(coupl) /= (/ (j,j=1,npirr) /) ))&
         & call error_handler("sm/ccc...: strange type of coupling (2)")

    DWRITE(*,'(A)') '----------------------------'
  end subroutine cc_coupling

  subroutine large_to_small_salcs(pcp, us, saL, saS)
    use group_module, only: n_p_irr=>group_num_re
    use clebsch_gordan, only: cg_alloc,prod_bas
    implicit none
    integer(IK),intent(in)                        :: pcp(:)
    !^ pseudo-coupling: pcp(n-proj-irreps)
    type(symm_basis_transformation_c),intent(in)  :: us(:) ! us(n-proj-irreps)
    !^ transformation of pseudo-scaled basis to a standard form
    type(psalcs_ua),intent(in)                     :: saL(:)
    type(psalcs_ua),intent(inout)                  :: saS(:)
    !^ SALC`s for large and small components
    ! *** end of interface ***

    type(sym_prod), pointer :: i_sp, j_sp

    integer(IK)         :: i
    integer(IK)         :: lmax,n_ua,ua,L,i_irr,j_irr,i_c,mult,m
    integer(IK)         :: n_ea

    DPRINT 'sm/L2S_salcs: entered'

    if(n_p_irr /= size(pcp) .or. n_p_irr /= size(us))&
         & call error_handler("sm/L2S_salcs: shapes wrong")
    if(any(shape(saL) /= shape(saS)))&
         & call error_handler("sm/L2S_salcs: saL & saS shapes differ")
    if(any(pcp(pcp) /= (/(i,i=1,n_p_irr)/)))&
         & call error_handler("sm/L2S_salcs: pcoupling is strange")

    n_ua  = size(saL)

    do ua=1,n_ua
       lmax = saL(ua)%lmax
       n_ea = saL(ua)%n_ea

       call alloc_psalcs_ua(lmax, n_ea, n_p_irr, saS(ua))

       irr_:do i_irr=1,n_p_irr
          j_irr = pcp(i_irr)    ! coupled irrep
          ! U => us(i_irr)%matrix ! unitary adjustment matrix

          do i_c=1,2 ! spinor component
             do L=0,saL(ua)%lmax
                i_sp => saL(ua)%mos(i_irr, L, i_c)
                j_sp => saS(ua)%mos(j_irr, L, i_c)

                mult = i_sp%mult ! same for all i_c,L
                call cg_alloc(mult, j_sp)

!!$                DPRINT 'sm/L2S_salcs: ua=',ua,' i_irr=',i_irr,' j_irr=',j_irr,&
!!$                     & ' i_c=',i_c,' L=',L,' mult=',mult

                do m = 1, mult ! multiplicity of the subspace
                   j_sp%sub(m) = rotate(i_sp%sub(m), us(i_irr)%matrix)
!!$                   DPRINT 'sm/L2S_salcs: mult=',m,':'
!!$                   DPRINT 'sm/L2S_salcs: j_sp%sub(m)%re=',j_sp%sub(m)%re
!!$                   DPRINT 'sm/L2S_salcs: j_sp%sub(m)%im=',j_sp%sub(m)%im
                enddo
             enddo
          enddo

       enddo irr_
    enddo
  contains

    function rotate(pbL, U) result(pbS)
      implicit none
      type(prod_bas), intent(in) :: pbL
      complex(CK), intent(in) :: U(:,:)
      type(prod_bas) :: pbS
      ! *** end of interface ***

      integer(IK) :: M, N, ns(3)

      ! ns(1) == na == n_partn ! same
      ! ns(2) == nb == n_ea    ! for all
      ! ns(3) == nc == n_lm    ! subspaces
      ns = shape(pbL%z)

      !
      ! See intent(out) for the type with allocatable components:
      !
      call cg_alloc(ns(1), ns(2), ns(3), pbS, cmplx=.true.)

      M = ns(1) ! na
      N = ns(2) * ns(3) ! nb *  nc

      !
      ! Custom matrix multiplication:
      !
      call mm(M, N, U, pbL%Z, pbS%Z)

      pbS%RE =  REAL(pbS%Z)
      pbS%IM = AIMAG(pbS%Z)
    end function rotate

    subroutine mm(M,N,A,B,C)
      ! very specific matrix multiplication:
      ! C(M,N) = tr(A(M,M)) * B(M,N)
      implicit none
      integer(IK),intent(in)  :: M,N
      complex(CK),intent(in)  :: A(M,*), B(M,*) ! A(M,M),B(M,N)
      complex(CK),intent(out) :: C(M,*)         !        C(M,N)
      ! *** end of interface ***

      complex(CK),parameter :: zero = (0.0_rk,0.0_rk)
      integer(IK)           :: i,j,k

      do i=1,M
         do k=1,N
            C(i,k) = zero
            do j=1,M
               C(i,k) = C(i,k) + A(j,i) * B(j,k)
            enddo
         enddo
      enddo

!!$      do i=1,M
!!$         do k=1,N
!!$            C(i,k) = zero
!!$            do j=1,N
!!$               C(i,k) = C(i,k) + A(i,j)*B(j,k)
!!$            enddo
!!$         enddo
!!$      enddo
    end subroutine mm

  end subroutine large_to_small_salcs

  !*************************************************************
  subroutine symm_dipole_selrules_gen
    !  Purpose: generates the dipole selection rules
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    use efm_decl
    use efm_module
    use symmetry_data_module
    implicit none
    !------------ Declaration of local constants  ----------------
    real(RK),parameter :: small = 1.e-10_rk
    ! very small value
    type(efm_cscoII_diag_invspace)  :: op_adapt
    ! symmetry adaption of dipole operators
    integer(IK)                :: i,alloc_stat,i_irrep,&
         j_irrep,i_partner,dim_irrep,component,vector
    ! counters
    integer(IK)  :: dipoleop_irrep(3)
    ! irrep of dipole operators

    integer(IK)  :: mult,isym,xyz
    real(RK)     :: coeff

    logical :: io
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code ------------------------------------

    ! Slaves do not always have output_unit open:
    io = output_unit > 0 .and. .not. no_output_unit_output

    DPRINT  'symm_dipole_selrules_gen: entered'

    ! allocLate eigenspaces of symmetrized dipole operators
    allocate(op_adapt%csco_eigenspaces(group_num_ir),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symm_symmadapt: allocation of op_adapt(ang_mom)%csco_eigenspaces failed")

    ! allocate data structure of symmetry data module
    allocate(symmetry_data_dipoles_exist(group_num_ir,group_num_ir,3),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symm_symmadapt: allocation of symmetry_data_dipole_exist failed")
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       ! allocate data structure of symmetry data module
       allocate(symmetry_data_pdipoles_exist(group_num_pir,group_num_pir,3),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("symm_symmadapt: allocation of symmetry_data_dipole_exist failed")
       ! initialize symmetry_data_dipole_exist
       symmetry_data_pdipoles_exist = .false.
    endif

    ! initialize symmetry_data_dipole_exist
    symmetry_data_dipoles_exist = .false.
    
    ! deallocate dipole operators
    call efm_cscoII_solve(ylm_trafos(1)%matrix,3,op_adapt)


    symmetry_data_dip_irrep_mult = 0
    symmetry_data_dip_components = 0.0_RK
    isym = 0

    ! determine irreps of dipole operators
    ! loop over irreps
    do i_irrep =1,group_num_ir
       if (.not.op_adapt%csco_eigenspaces(i_irrep)%exists) cycle

       dim_irrep = irrep_can(i_irrep)%dimension

       ! loop over partners
       do i_partner =1,dim_irrep

          mult = op_adapt%csco_eigenspaces(i_irrep)%cscoII_subspaces(i_partner)%multiplicity

          symmetry_data_dip_irrep_mult(i_irrep) = mult
          DPRINT  'symm_dipole_selrules_gen: mult=',mult,' of (',i_irrep,i_partner,')'

          ! loop over multiplicity
          do vector = 1,mult
             isym = isym + 1

             ! loop over components x,y,z
             do component = 1,3
                select case(component)
                case (1)
                   xyz = 1
                case (2)
                   xyz = 3
                case (3)
                   xyz = 2
                end select

                coeff = op_adapt%csco_eigenspaces(i_irrep)%cscoII_subspaces(i_partner)%&
                     &eigenfunction(component,vector)

                symmetry_data_dip_components(isym,xyz) = coeff

                if (abs(coeff).gt.small) then
                   dipoleop_irrep(xyz) = i_irrep
!!$                      exit
                end if
             enddo
             DPRINT  'symm_dipole_selrules_gen: isym=',isym,' (',i_irrep,i_partner,vector,')'
             DPRINT  'symm_dipole_selrules_gen: coef=',symmetry_data_dip_components(isym,:)
          enddo
       enddo
    enddo
    ASSERT(isym==3)

    ! determine irreps of dipole operators
    !do i_irrep =1,group_num_ir
    !   if (op_adapt%csco_eigenspaces(i_irrep)%exists) then
    !      dim_irrep = irrep_can(i_irrep)%dimension
    !      do i_partner =1,dim_irrep
    !         if (abs(op_adapt%csco_eigenspaces(i_irrep)%cscoII_subspaces(i_partner)%&
    !              &eigenfunction(1,1)).gt.small) then
    !            dipoleop_irrep(1) = i_irrep
    !         end if
    !         if (abs(op_adapt%csco_eigenspaces(i_irrep)%cscoII_subspaces(i_partner)%&
    !              &eigenfunction(3,1)).gt.small) then
    !            dipoleop_irrep(2) = i_irrep
    !         end if
    !         if (abs(op_adapt%csco_eigenspaces(i_irrep)%cscoII_subspaces(i_partner)%&
    !              &eigenfunction(2,1)).gt.small) then
    !            dipoleop_irrep(3) = i_irrep
    !         end if
    !      end do
    !   end if
    !end do

    ! now determine dipole selection rules
    if (io) then
        write(output_unit,*) "Dipole Selection Rules"
        write(output_unit,*)
    endif
    ! loop over dipole component
    do i=1,3
       if (io) then
           if (i.eq.1) then
              write(output_unit,*) "X couples:"
           elseif (i.eq.2) then
              write(output_unit,*) "Y couples:"
           elseif (i.eq.3) then
              write(output_unit,*) "Z couples:"
           endif
       endif
       ! loop over irreps
       do i_irrep =1,group_num_ir
          do j_irrep =1,group_num_ir
             if ( efm_cg( i_irrep, dipoleop_irrep(i) )%csco_eigenspaces(j_irrep)%exists ) then
                if (io) then
                    write(output_unit,*) irrep_can(i_irrep)%label," and ",irrep_can(j_irrep)%label
                endif
                symmetry_data_dipoles_exist(i_irrep,j_irrep,i) = .true.
             end if
          end do
       end do
    end do

    ! now determine dipole selection rules for projective irreps
    if (io) then
        write(output_unit,*) "Dipole Selection Rules for projective Irreps"
        write(output_unit,*)
    endif
    ! loop over dipole component
    do i=1,3
       if (io) then
           if (i.eq.1) then
              write(output_unit,*) "X couples:"
           elseif (i.eq.2) then
              write(output_unit,*) "Y couples:"
           elseif (i.eq.3) then
              write(output_unit,*) "Z couples:"
           endif
       endif
       ! loop over irreps
       do i_irrep =1,group_num_pir
          do j_irrep =1,group_num_pir
             if ( efm_proj_cg( dipoleop_irrep(i), i_irrep  )%csco_eigenspaces(j_irrep)%exists ) then
                if (io) then
                    write(output_unit,*) proj_irrep_can(i_irrep)%label," and ",proj_irrep_can(j_irrep)%label
                endif
                symmetry_data_pdipoles_exist(i_irrep,j_irrep,i) = .true.
             end if
          end do
       end do
    end do

    if (io) then
        write(output_unit,*) "Symmetrized Dipole Operators"
    endif
    ! show symmetrized dipole operators
    call efm_show_invspace(op_adapt)

    ! deallocate symmetrized dipole operators
    call efm_free(op_adapt)
    
  end subroutine symm_dipole_selrules_gen
  !*************************************************************

  !------------------------------------------------------------------
  !
  ! MEMORY MANAGEMENT (SPOR ONLY!!!)
  !
  !

  subroutine alloc_vsalcs_ua(lmax,n_ea,n_irr,salcs)
    use error_module
    implicit none
    integer(IK),intent(in)        :: lmax,n_ea,n_irr
    type(vsalcs_ua),intent(inout) :: salcs
    ! *** end of interface ***

  
    integer(IK) :: memstat

    salcs%lmax  = lmax
    salcs%n_irr = n_irr
    salcs%n_ea  = n_ea

    allocate(salcs%mos(n_irr,0:lmax),STAT=memstat)
    call error(memstat,"sm/alloc_vsalcs_ua: alloc failed")
  end subroutine alloc_vsalcs_ua

  subroutine free_vsalcs_ua(salcs)
    use clebsch_gordan
    implicit none
    type(vsalcs_ua),intent(inout) :: salcs
    ! *** end of interface ***

    integer(IK) :: irr,l
    integer(IK) :: memstat

    if(salcs%lmax/=ubound(salcs%mos,2))&
         & call error_handler("sm/free_vsalcs_ua: wrong mos shape")

    do irr=1,salcs%n_irr
       do l=0,salcs%lmax
          call cg_free(salcs%mos(irr,l))
       enddo
    enddo
    deallocate(salcs%mos,STAT=memstat)
    ASSERT(memstat==0)
    salcs%lmax  = -1
    salcs%n_irr = -1
    salcs%n_ea  = -1
  end subroutine free_vsalcs_ua

  subroutine alloc_psalcs_ua(lmax,n_ea,n_irr,salcs)
    use clebsch_gordan
    use error_module
    implicit none
    integer(IK),intent(in)       :: lmax,n_ea,n_irr
    type(psalcs_ua),intent(inout) :: salcs
    ! *** end of interface ***

    integer :: memstat

    salcs%lmax  = lmax
    salcs%n_irr = n_irr
    salcs%n_ea  = n_ea

    allocate(salcs%mos(n_irr,0:lmax,2),STAT=memstat)
    call error(memstat,"sm/alloc_psalcs_ua: alloc failed")
  end subroutine alloc_psalcs_ua

  subroutine free_psalcs_ua(salcs)
    use clebsch_gordan
    implicit none
    type(psalcs_ua),intent(inout) :: salcs
    ! *** end of interface ***

    integer(IK) :: i,l,c
    integer(IK) :: memstat

    if(salcs%lmax/=ubound(salcs%mos,2))&
         & call error_handler("sm/free_psalcs_ua: wrong mos shape")
    
    do i=1,salcs%n_irr
       do l=0,salcs%lmax
          do c=1,2
             call cg_free(salcs%mos(i,l,c))
          enddo
       enddo
    enddo
    deallocate(salcs%mos,STAT=memstat)
    ASSERT(memstat==0)
    salcs%lmax  = -1
    salcs%n_irr = -1
    salcs%n_ea  = -1
  end subroutine free_psalcs_ua

  !-------------------------------------------------------
  !
  ! the pair for type(symm_basis_transformation_c) arrays
  !
  subroutine alloc_arr_trafo_c(n,dims,ts)
    use group_module, only: symm_basis_transformation_c
    implicit none
    integer(IK),intent(in)                          :: n,dims(:)
    type(symm_basis_transformation_c),intent(inout) :: ts(:)
    ! *** end of interface ***

    integer(IK) :: i
    integer :: memstat

    if(n /= size(ts) .or. n/= size(dims))&
         & call error_handler("sm/alloc_atrafo_c: wrong shapes")

    do i=1,n
       allocate(ts(i)%matrix(dims(i),dims(i)),STAT=memstat)
       if(memstat/=0)call error_handler("sm/alloc_atrafo_c: alloc failed")
    enddo
  end subroutine alloc_arr_trafo_c

  subroutine free_arr_trafo_c(ts)
    use group_module, only: symm_basis_transformation_c
    implicit none
    type(symm_basis_transformation_c),intent(inout) :: ts(:)
    ! *** end of interface ***

    integer(IK) :: i,n
    integer :: memstat

    n=size(ts)

    do i=1,n
       deallocate(ts(i)%matrix,STAT=memstat)
       if(memstat/=0)call error_handler("sm/free_atrafo_c: dealloc failed")
    enddo
  end subroutine free_arr_trafo_c
  
!!$  subroutine show_invspace(inv)
!!$    use group_module, only: irrep_can
!!$    use efm_module, only: efm_cscoII_diag_invspace
!!$    implicit none
!!$    type(efm_cscoII_diag_invspace),intent(in) :: inv
!!$    ! *** end of interface ***
!!$
!!$    integer(IK) :: irr,n,m,dim_irrep,mult_irrep,dim_trafo
!!$
!!$    DPRINT 'shis: projective=',inv%projective
!!$    do irr=1,size(inv%csco_eigenspaces)
!!$       if(inv%csco_eigenspaces(irr)%exists)then
!!$          DPRINT 'shis: irrep ',irrep_can(irr)%label,' exists'
!!$       else
!!$          cycle
!!$       endif
!!$
!!$       if(irrep_can(irr)%pseudo)then
!!$          DPRINT 'shis: (pseudo)'
!!$       endif
!!$
!!$       dim_irrep  = size(inv%csco_eigenspaces(irr)%cscoII_subspaces)
!!$       dim_trafo  = size(inv%csco_eigenspaces(irr)%cscoII_subspaces(1)%eigenfunction,1)
!!$       mult_irrep = size(inv%csco_eigenspaces(irr)%cscoII_subspaces(1)%eigenfunction,2)
!!$
!!$
!!$       ! loop over the multiples
!!$       do n=1,mult_irrep
!!$          ! loop over all partners of the irrep
!!$          do m=1,dim_irrep
!!$             write(*,'(" shis:",2I4,(11F7.3),/:(14X,11F7.3))') n,m,&
!!$                  & inv%csco_eigenspaces(irr)%cscoII_subspaces(m)&
!!$                  & %eigenfunction(:,n)
!!$          end do
!!$       end do
!!$    enddo
!!$  end subroutine show_invspace
!!$
  subroutine show_fine_c(U)
    implicit none
    complex(CK),intent(in) :: U(:,:)
    ! *** end of interface ***

    integer(IK) :: n,m,i,j

    n = size(U,1)
    m = size(U,2)

    do j=1,m
#ifndef _COMPAC_FORTRAN
       write(*,'(20(SSF6.3,SPF6.3,"i  "))')&
            & ( U(i,j), i=1,n )
#else
       write(*,'(20(F6.3,F6.3,"i  "))')&
            & ( U(i,j), i=1,n )
#endif
    enddo
  end subroutine show_fine_c

  subroutine show_fine_r(U)
    implicit none
    real(RK),intent(in) :: U(:,:)
    ! *** end of interface ***

    integer(IK) :: n,m,i,j

    n = size(U,1)
    m = size(U,2)

    do j=1,m
#ifndef _COMPAC_FORTRAN
       write(*,'(20(SSF6.3,6X,"   "))')&
            & ( U(i,j), i=1,n )
#else
       write(*,'(20(F6.3,6X,"   "))')&
            & ( U(i,j), i=1,n )
#endif
    enddo
  end subroutine show_fine_r
!!$
!!$  function examine_cmplx(M_real,M_imag) result(res)
!!$    implicit none
!!$    real(RK),dimension(:,:),intent(in) :: M_real,M_imag
!!$    real(RK),dimension(4) :: res ! <<<result
!!$    ! *** end of interface ***
!!$
!!$    real(RK) ::&
!!$         & diag_norm_real, diag_norm_imag,&
!!$         & offdiag_norm_real, offdiag_norm_imag
!!$    integer(IK) :: n
!!$
!!$    n = size(M_real,1)
!!$
!!$    diag_norm_real =   sum(pack(M_real,diag_mask(n))**2)
!!$    diag_norm_imag =   sum(pack(M_imag,diag_mask(n))**2)
!!$    offdiag_norm_real= sum(pack(M_real,offdiag_mask(n))**2)
!!$    offdiag_norm_imag= sum(pack(M_imag,offdiag_mask(n))**2)
!!$
!!$    res=(/ diag_norm_real,offdiag_norm_real,&
!!$         & diag_norm_imag,offdiag_norm_imag /)
!!$    return
!!$  end function examine_cmplx
!!$
!!$  function examine_real(M) result(res)
!!$    implicit none
!!$    real(RK),dimension(:,:),intent(in) :: M
!!$    real(RK),dimension(2) :: res ! <<<result
!!$    ! *** end of interface ***
!!$
!!$    real(RK) ::&
!!$         & diag_norm_real,offdiag_norm_real
!!$    integer(IK) :: n
!!$
!!$    n = size(M,1)
!!$
!!$    diag_norm_real =   sum(pack(M,diag_mask(n))**2)
!!$    offdiag_norm_real= sum(pack(M,offdiag_mask(n))**2)
!!$
!!$    res=(/ diag_norm_real,offdiag_norm_real /)
!!$    return
!!$  end function examine_real
!!$
  function diag_mask(n) result(m)
    integer(IK),intent(in) :: n
    logical,dimension(n,n) :: m !<<<result
    ! *** end of interface ***

    integer(IK) :: i

    m = .false.
    do i=1,n
       m(i,i) = .true.
    enddo
  end function diag_mask

  function offdiag_mask(n) result(m)
    integer(IK),intent(in) :: n
    logical,dimension(n,n) :: m !<<<result
    ! *** end of interface ***

    integer(IK) :: i

    m = .true.
    do i=1,n
       m(i,i) = .false.
    enddo
  end function offdiag_mask

  subroutine rephase(pcp,us)
    use efm_module, only: phase
    implicit none
    integer(IK),intent(in)                          :: pcp(:)
    type(symm_basis_transformation_c),intent(inout) :: us(:)
    ! *** end of interface ***
    
    logical     :: processed(size(pcp)),trivial
    complex(CK), allocatable :: uw(:, :)
    complex(CK) :: ph,phi,phj
    integer(IK) :: n,m,i,j
    real(RK)    :: eps
    integer     :: memstat

    n = size(pcp)

    if(n /= size(us))&
         & call error_handler("sm/rephase: sizes differ")
    if(any(pcp(pcp) /= (/(i,i=1,n)/)))&
         & call error_handler("sm/rephase: strange coupling")

    eps = 10000.0*epsilon(1.0_rk)

    processed = .false.
    do i = 1,n
       j = pcp(i)
       if(processed(i).or.processed(j))cycle

       processed(i) = .true.
       processed(j) = .true.

       DPRINT 'sm/rephase: processing irreps ',i,' <-> ',j

       trivial = (i==j)

       if(trivial)then

          DPRINT '(trivial coupling)'

          ! ui => us(i)%matrix

          DPRINT '---------- Ui before ---'
          DCALL show_fine(us(i)%matrix)

          if(.not.unitary(us(i)%matrix))&
               & call error_handler("sm/rephase: ui is not unitary")

          ph = phase(us(i)%matrix)

          us(i)%matrix = ph * us(i)%matrix

          if(.not.unity(us(i)%matrix))then
             DPRINT 'sm/rephase, WARNING: coupling is trivial BUT basis changes'
          endif

          DPRINT 'sm/rephase: ui rephased with ',ph

          DPRINT '---------- Ui after ----'
          DCALL show_fine(us(i)%matrix)
       else

          ! ui => us(i)%matrix
          ! uj => us(j)%matrix
          DPRINT '---------- Ui before ---'
          DCALL show_fine(us(i)%matrix)
          DPRINT '---------- Uj before ---'
          DCALL show_fine(us(j)%matrix)

          m = size(us(i)%matrix, 1)

          allocate(uw(m,m),STAT=memstat)
          if(memstat/=0)call error_handler("sm/rephase: alloc failed")

          if(.not.(unitary(us(i)%matrix) .and. unitary(us(j)%matrix)))&
               & call error_handler("sm/rephase: ui .or uj is not unitary")

          phi = phase(us(i)%matrix)
          phj = phase(us(j)%matrix)

          us(i)%matrix  = phi * us(i)%matrix
          us(j)%matrix  = phj * us(j)%matrix

          DPRINT 'sm/rephase: ui and uj rephased with '
          DPRINT phi,' and ',phj

          uw = matmul(us(i)%matrix, us(j)%matrix)

          ph = phase(uw)

          uw = ph * uw

          if(.not.unity(uw))then
             call error_handler("sm/rephase: ui * uj is not a pure phase factor")
          else
             DPRINT 'sm/rephase: ui * uj = ph*1, where ph=',ph
          endif

          ph = sqrt(ph)

          us(i)%matrix  = ph * us(i)%matrix
          us(j)%matrix  = ph * us(j)%matrix

          DPRINT 'sm/rephase: both ui and uj rephased again with ',ph

          DPRINT '---------- Ui after ----'
          DCALL show_fine(us(i)%matrix)
          DPRINT '---------- Uj after ----'
          DCALL show_fine(us(j)%matrix)

          uw  = us(i)%matrix - conjg(transpose(us(j)%matrix))

          if(any(abs(uw) > eps))&
               & call error_handler("sm/rephase: why ui /= uj+ ???")

          deallocate(uw,STAT=memstat)
          if(memstat/=0)call error_handler("sm/rephase: dealloc failed")
       endif
    enddo
  end subroutine rephase

  function unitary(u) result(yes)
    implicit none
    complex(CK),intent(in) :: u(:,:)
    logical                :: yes !<<< result
    ! *** end of interface ***
    
    complex(CK) :: w(size(u,1),size(u,2))

    yes = (size(u,1) == size(u,2))

    if(.not.yes)return

    w = matmul(u,conjg(transpose(u)))
    yes = unity(w)

    if(.not.yes)return

    w = matmul(conjg(transpose(u)),u)
    yes = unity(w)
  end function unitary

  function unity(u) result(yes)
    implicit none
    complex(CK),intent(in) :: u(:,:)
    logical                :: yes !<<< result
    ! *** end of interface ***
    
    complex(CK),parameter :: one  = (1.0_ck,0.0_ck)
    complex(CK),parameter :: zero = (0.0_ck,0.0_ck)
    real(RK)              :: eps
    integer(IK)           :: n

    eps = 10000.0*epsilon(1.0_rk)

    yes = (size(u,1) == size(u,2))

    if(.not.yes)return

    n = size(u,1)

    yes = all(abs( u - merge(one,zero,diag_mask(n)) ) < eps)
  end function unity

end module symm_module

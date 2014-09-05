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
module density_data_module
  !-------------------------------------------------------------------
  !-------------- Module specification ---------------------------
  !-------------------------------------------------------------------
  !
  !  Purpose: contains the density matrix DENSMAT and all
  !           routines for
  !                - reading it from file
  !                - setting up an arbitrary initial matrix
  !                - writing it (to somewhere ...)
  !                - allocating and freeing storage for it
  !                - generating the density matrix
  !                  (originally done in density_module and chargefit)
  !
  !
  !  Contents:Variables(PUBLIC)  densmat
  !                               first index : IRREP
  !                               m:
  !                               first and second index : Orbitals
  !                                         third index  : Spin
  !
  !           Routines (PUBLIC)  density_data_alloc
  !                              density_data_free
  !                              print_densmat -> print densat to file
  !                                               $TTFSOUT/densmat.new
  !                             gendensmat_occ -> generate the density matrix
  !                                               based on the data available
  !                                               in the occupied_levels_module
  !
  !  Module called by: all routines that need DENSMAT
  !
  !
  !  References: none
  !
  !
  !  Author: Folke Noertemann
  !  Date: 11/95
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   2/7/97
  ! Description: GENDENSMAT and PGENDENSMAT from DENSITY_MODULE
  !              included here. The parts of CHARGEFIT that setup
  !              the density matrix on each slave (and the master)
  !              are collected in a new subroutine: CALC_DENSMAT
  !
  !              Generation of the density matrix is now based on the
  !              occupied orbitals only (as send by send_eigvec_occ)
  !              Date: 23/7/97
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification
  ! Author: MM
  ! Date:   10/97
  ! Description: Extension to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   3/02
  ! Description: subroutines save_densmat and open_densmat were included
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !------------ Modules used -----------------------------------------
# include "def.h"
! define FPP_SERIAL /* if parallel version of gendensmat() breaks */
  use type_module ! type specification parameters
  use datatype, only: arrmat2, arrmat3
  use symmetry_data_module, only: ssym
  ! use convergence_module leads to cyclic dependencies, thus
  ! communication with the convergence_module is set up by means
  ! of the formal parameter density_dev passed to gendensmat_occ
  implicit none
  save
  private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of public constants and variables -----
  type(arrmat3), allocatable, public, protected :: densmat(:)

  ! in case of spin orbit the density matrix is complex and there is no spin component
  type(arrmat2), allocatable, public, protected :: densmat_real(:)
  type(arrmat2), allocatable, public, protected :: densmat_imag(:)

  ! in case of creating a core density
  type(arrmat3), allocatable, public, protected :: core_densmat(:)

  !------------ public functions and subroutines ---------------------
  public :: print_densmat, &
       save_densmat, open_densmat

  public :: density_data_alloc!(), does no communication
  public :: density_data_free!(), does no communication
  public :: gendensmat_occ!(density_deviation), to be executed in parallel context
  public :: gendensmat_tot!( irrep, dmat ), for one irepp only

  !===================================================================
  ! End of public interface of module
  !===================================================================
  public :: arrmat3 !,sym

  integer(i4_kind), parameter :: &
       DM_NRDM = 1, &
       DM_SODM = 2, &
       DM_CORE = 3


  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************

  subroutine do_alloc_densmat(what, dims, n_spin)
    !
    !  Purpose: allocate the apporpriate space for the
    !           density matrix densmat and initialize it
    !           just in case ...
    !
    use init_module, only: init
    implicit none
    integer(i4_kind), intent(in) :: what
    integer(i4_kind), intent(in) :: dims(:) ! (n_irr)
    integer(i4_kind), intent(in) :: n_spin
    ! *** end of interface ***

    select case(what)
    case(DM_SODM)
       !
       ! SPIN ORBIT
       !
       call alloc2(dims, densmat_real)
       call alloc2(dims, densmat_imag)

       ! FIXME: set to zero?
       call init(densmat_real)
       call init(densmat_imag)

    case(DM_NRDM)
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call alloc3(dims, n_spin, densmat)

       ! FIXME: set to zero?
       call init(densmat)

    case(DM_CORE)
       call alloc3(dims, n_spin, core_densmat)

       ! FIXME: set to zero?
       call init(core_densmat)

    case default
       ABORT('no such case')
    end select
  end subroutine do_alloc_densmat

  subroutine density_data_alloc()
    !
    ! Purpose: wrapper for the do_alloc_densmat()
    !
    use options_module, only: options_spin_orbit
    use operations_module, only: operations_core_density
    use symmetry_data_module, only: ssym
    implicit none
    !** End of interface *****************************************

    ! Set appropriate  dimensions of irreps (use  projective irreps in
    ! case of spin orbit)
    if (options_spin_orbit) then
       call do_alloc_densmat(DM_SODM, ssym%dim_proj, 1)
    else
       call do_alloc_densmat(DM_NRDM, ssym%dim     , ssym%n_spin)
    endif

    if(operations_core_density)then
       ASSERT(.not.options_spin_orbit)
       call do_alloc_densmat(DM_CORE, ssym%dim     , ssym%n_spin)
    endif
  end subroutine density_data_alloc

  !*************************************************************

  subroutine do_free_densmat(what)
    !
    ! purpose: free the space allocated in  density_data_alloc
    !
    implicit none
    integer(i4_kind), intent(in) :: what
    ! *** end of interface ***

    select case(what)
    case(DM_SODM)
       !
       ! SPIN ORBIT
       !
       call free2(densmat_real)
       call free2(densmat_imag)

    case(DM_NRDM)
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call free3(densmat)

    case(DM_CORE)
       call free3(core_densmat)

    case default
       ABORT('no such case')
    end select
  end subroutine do_free_densmat

  subroutine density_data_free()
    !
    ! Purpose: wrapper for the do_free_densmat()
    !
    use options_module, only: options_spin_orbit
    use operations_module, ONLY : operations_core_density
    implicit none
    ! *** end of interface ***

    if (options_spin_orbit) then
       call do_free_densmat(DM_SODM)
    else
       call do_free_densmat(DM_NRDM)
    endif

    if(operations_core_density)then
       ASSERT(.not.options_spin_orbit)
       call do_free_densmat(DM_CORE)
    endif
  end subroutine density_data_free


  subroutine print_densmat(loop)
    ! Purpose: print out the density matrix to file
    !          $TTFSOUT/densmat.new.
    ! Modules used -----------------------------------------------
    use filename_module, only: outfile
    use iounitadmin_module  ! get_iounit,return_iounit
    use symmetry_data_module, only: symmetry_data_n_irreps, &
         symmetry_data_dimension, &
         symmetry_data_n_spin
    use operations_module, ONLY : operations_core_density
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in),optional :: loop
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: io_u,count,i,is,m,n,i_dim
    data count / 0 /

    if(.not.present(loop)) count=count+1

    io_u=get_iounit()
    open(io_u,status='unknown',form='formatted', &
         position='append', &
         file=trim(outfile('densmat.new')))

    write(io_u,*)' '
    if(present(loop)) then
       write(io_u,*)' +++++++++++++ Loop ',loop,' +++++++++++++++'
    else
       write(io_u,*)' +++++++++++++ Loop ',count,' +++++++++++++++'
    endif

    do i=1,symmetry_data_n_irreps()
       write(io_u,*)' --- Irrep ',i,' ---'

       do is=1,symmetry_data_n_spin()

          i_dim=symmetry_data_dimension(i)
          do m=1,i_dim
             write(io_u,'(4(A,F13.10))')(' ',densmat(i)%m(m,n,is), n=1,i_dim)
          enddo
          if (operations_core_density) then
             write(io_u,*) 'core density:'
             do m=1,i_dim
                write(io_u,'(4(A,F13.10))')(' ',core_densmat(i)%m(m,n,is), n=1,i_dim)
             enddo
          endif
          write(io_u,*) ' '

       enddo!is

    enddo!i


    close(io_u)
    call return_iounit(io_u)

    return
  end subroutine print_densmat


  subroutine gendensmat_tot(irr, P)
    !
    ! FIXME: Duplicates functionality of gendensmat_occ()
    !        needs to be parallelized.
    !
    use occupied_levels_module, only: eigvec_occ, occ_num_occ
    implicit none
    integer(i4_kind), intent(in)  :: irr
    real(r8_kind)   , intent(out) :: P(:,:)
    ! *** end of interface ***

    integer(i4_kind) :: spin, n, m

    n = size(eigvec_occ(irr)%m, 1)
    ASSERT(size(P,1)==n)
    ASSERT(size(P,2)==n)
    m = size(eigvec_occ(irr)%m,2)
    n = size(occ_num_occ(irr)%m,1)
    ASSERT(m<=n)

    P = 0.0
    do spin = 1, size(eigvec_occ(irr)%m, 3) ! UP and DOWN
       ! use default occupation numbers:
       P = P + abat(eigvec_occ(irr)%m(:, :, spin), occ_num_occ(irr)%m(:m, spin))
    enddo
  end subroutine gendensmat_tot

  !*************************************************************

  subroutine gendensmat_occ (density_dev)
    !
    ! Purpose:  Assemble density matrix  from occupied  eigenstates as
    ! stored in  eigvec_occ and occ_num_occ.   Allocates denity matrix
    ! if it is not already allocated.
    !
    ! Input parameter (not modified on output):
    !
    !   Name:        Description/Range:
    !   ssym         ssymmetry stuff
    !   n_occo       the index of the heighest non-empty orbital
    !   occ_num_occ  occupation number of the occupied states
    !   eigvec_occ   the orbital coefficient of the occupied states
    !
    ! Subroutine called by: main_scf(), main_gradient(),
    ! properties_main().
    !
    ! FIXME: the presence of the argument toggles more allocations for
    ! space to  store the old  density matrix in  order to be  able to
    ! compute  the difference.   This argument  was present  on master
    ! (call   from  main_scf())   but   not  on   slaves  (call   from
    ! main_slave()).  Now that it  is called  from a  paralle context,
    ! this argument is present on every worker.
    !
    !  Author: Uwe Birkenheuer (based on gendensmat of TG)
    !  Date: 7/97
    !
    !---------------------------------------------------------------
    ! Modifications
    !---------------------------------------------------------------
    !
    ! Modification (Please copy before editing)
    ! Author: ...
    ! Date:   ...
    ! Description: ...
    !--------------------------------------------------------------
    ! MODULES
    !--------------------------------------------------------------
    use iounitadmin_module, only: output_unit
    use options_module, only: options_spin_orbit
    use output_module, only: output_n_density_dev
    implicit none
    real (r8_kind), optional, intent (out) :: density_dev
    ! *** end of interface ***

    integer(kind=i4_kind)       :: i_gamma
    logical                     :: present_density_dev

    !
    ! Temp storage for computing the differences:
    !
    type(arrmat3), allocatable :: save_densmat(:)
    type(arrmat2), allocatable :: save_densmat_real(:)
    type(arrmat2), allocatable :: save_densmat_imag(:)

    real(r8_kind) :: dev(max(output_n_density_dev, 1))

    ! make sure that densmat is really allocated
    ! this should be the case only before the first cycle
    if (options_spin_orbit) then
       if( .not.allocated(densmat_real) )then
          call density_data_alloc()
       endif
    else
       if( .not.allocated(densmat) ) then
          call density_data_alloc()
       endif
    endif

    !
    ! This is the case only on master who is also interested in
    ! how different the density matrix from the previous scf
    ! iteration:
    !
    present_density_dev = present(density_dev)

    !
    ! Save old density matrix:
    !
    if (present_density_dev) then
       if (options_spin_orbit) then
          call copy2(densmat_real, save_densmat_real)
          call copy2(densmat_imag, save_densmat_imag)
       else
          call copy3(densmat, save_densmat)
       endif
    endif

    !
    ! Now set up  the new density matrix, this sub  does the real work
    ! and updates global densmat/densmat_real/densmat_imag:
    !
    call gendensmat()

    !
    ! Compute differences to the old density matrix:
    !
    if (present_density_dev) then
       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          do i_gamma = 1, size(densmat_real)
            save_densmat_real(i_gamma)%m = save_densmat_real(i_gamma)%m - densmat_real(i_gamma)%m
            save_densmat_imag(i_gamma)%m = save_densmat_imag(i_gamma)%m - densmat_imag(i_gamma)%m

            save_densmat_real(i_gamma)%m = sqrt(save_densmat_real(i_gamma)%m**2 + save_densmat_imag(i_gamma)%m**2)
          enddo
          call maxvals2(save_densmat_real, dev)
       else
          !
          ! STANDARD SCF (NO SPIN ORBIT)
          !
          do i_gamma = 1, size(densmat)
            save_densmat(i_gamma)%m = abs(save_densmat(i_gamma)%m - densmat(i_gamma)%m)
          enddo
          call maxvals3(save_densmat, dev)
       endif

       if (options_spin_orbit) then
          call free2(save_densmat_real)
          call free2(save_densmat_imag)
       else
          call free3(save_densmat)
       endif

       !
       ! Return max deviation:
       !
       density_dev = dev(1)

       !
       ! FIXME:  Most of  the complexity  is  for these  lines in  the
       ! output which I never looked at.
       !
       if (output_n_density_dev > 0) then
          if (output_unit > 0) then
             write (output_unit, *) " "
             write (output_unit, *) "MAX. DIFF. IN DENSITY MATRIX      "
             write (output_unit, 1000) dev
             write (output_unit, *) " "
1000         format (2x, 6ES13.3)
          endif
       endif
    endif
  end subroutine gendensmat_occ

  subroutine gendensmat()
    !
    ! Assemble density  matrix from occupied eigenstates  as stored in
    ! eigvec_occ and occ_num_occ.  UB 7/97
    !
    !  Input parameter (not modified on output):
    !
    !  Name:              Description/Range:
    !  ssym               ssymmetry stuff
    !  n_occo             the index of the heighest non-empty orbital
    !  occ_num_occ        occupation number of the occupied states
    !  eigvec_occ         the orbital coefficient of the occupied states
    !
    !  Subroutine called by: pre_dens_master,chargefit
    !
    !  Author: Uwe Birkenheuer (based on gendensmat of TG)
    !  Date: 7/97
    !
    !---------------------------------------------------------------
    ! Modifications
    !---------------------------------------------------------------
    !
    ! Modification (Please copy before editing)
    ! Author: ...
    ! Date:   ...
    ! Description: ...
    !--------------------------------------------------------------
    ! MODULES
    !--------------------------------------------------------------
    use occupation_module, only: n_occo, n_occo_core
    use operations_module, only: operations_core_density
    use occupied_levels_module, only: eigvec_occ, occ_num_occ,eigvec_occ_real,eigvec_occ_imag, &
                                      occ_num_core_occ
    use options_module, only: options_spin_orbit
    implicit none
    ! *** end of interface ***

    integer(i4_kind)       :: i_gamma, i_spin
    integer(i4_kind)       :: nocco

    ! make sure that densmat is really allocated
    if (options_spin_orbit) then
      ASSERT(allocated(densmat_real))
      ASSERT(allocated(densmat_imag))
    else
      ASSERT(allocated(densmat))
    endif

    ! now set up the density matrix
    do i_gamma = 1, size(occ_num_occ) ! n_irrep

       do i_spin = 1, size(occ_num_occ(i_gamma)%m, 2) ! n_spin
             !
             ! Number of occupied orbitals for this (irrep, spin) pair:
             !
             nocco = n_occo(i_spin, i_gamma)

             !
             ! Occupation numbers (of occupied levels only):
             !
             ! occn => occ_num_occ(i_gamma)%m(:nocco, i_spin)

             if (options_spin_orbit) then
                !
                ! SPIN ORBIT
                !

                !
                ! Original code computed density matrix by the following
                ! procedure (loop over a <= b and sum over k is assumed)
                !
                ! dens_real(a, b) = eigv_real(a, k) * n(k) * eigv_real(b, k)
                !                 + eigv_imag(a, k) * n(k) * eigv_imag(b, k)
                !
                ! dens_imag(a, b) = eigv_real(a, k) * n(k) * eigv_imag(b, k)
                !                 - eigv_imag(a, k) * n(k) * eigv_real(a, k)
                !
                ! dens_real(b, a) =   dens_real(a, b)
                ! dens_imag(b, a) = - dens_imag(a, b)
                !
                ! FIXME: shouldnt it rather be "imag = imag * real - real * imag"
                !        as in imaginary part of D = V * n * V^H ?
                !        Do we flip the sign of I here?
                !        See also the comment about suspicious sign for Coulomb
                !        matrix in ham_calc_module!
                !

                !
                ! Compute contribution to the densmat due to these
                ! occupied levels, use BLAS if necessary:
                !
#ifndef FPP_SERIAL
                densmat_real(i_gamma)%m = pabat(eigvec_occ_real(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin)) &
                                        + pabat(eigvec_occ_imag(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin))
                densmat_imag(i_gamma)%m = pabct(eigvec_occ_real(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin), &
                                                eigvec_occ_imag(i_gamma)%m)
#else
                densmat_real(i_gamma)%m = abat(eigvec_occ_real(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin)) &
                                        + abat(eigvec_occ_imag(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin))
                densmat_imag(i_gamma)%m = abct(eigvec_occ_real(i_gamma)%m, occ_num_occ(i_gamma)%m(:nocco, i_spin), &
                                               eigvec_occ_imag(i_gamma)%m)
#endif
                densmat_imag(i_gamma)%m = densmat_imag(i_gamma)%m - transpose(densmat_imag(i_gamma)%m)
             else ! options_spin_orbit
                !
                ! STANDARD SCF (NO SPIN ORBIT)
                !
                !
                ! Compute contribution to the densmat due to these
                ! occupied levels, use BLAS if necessary:
                !
#ifndef FPP_SERIAL
                densmat(i_gamma)%m(:,:,i_spin) = pabat( eigvec_occ(i_gamma)%m(:,:,i_spin),&
                                                 occ_num_occ(i_gamma)%m(:nocco, i_spin))
#else
                densmat(i_gamma)%m(:,:,i_spin) = abat(eigvec_occ(i_gamma)%m(:,:,i_spin),&
                                                 occ_num_occ(i_gamma)%m(:nocco, i_spin))
#endif
             endif ! options_spin_orbit

             ! load the core density matrix (if required)
             if (operations_core_density) then
  ASSERT(.not. options_spin_orbit)
                nocco = n_occo_core(i_spin, i_gamma)
                ! occn => occ_num_core_occ(i_gamma)%m(:nocco, i_spin)
                WARN('untested')
#ifndef FPP_SERIAL
                core_densmat(i_gamma)%m(:,:,i_spin) = pabat(eigvec_occ(i_gamma)%m(:,:,i_spin),&
                                                    occ_num_core_occ(i_gamma)%m(:nocco, i_spin))
#else
                core_densmat(i_gamma)%m(:,:,i_spin) = abat(eigvec_occ(i_gamma)%m(:,:,i_spin),&
                                              occ_num_core_occ(i_gamma)%m(:nocco, i_spin))
#endif
             endif

          enddo
    enddo
  end subroutine gendensmat

#ifndef FPP_SERIAL
  function pabat(a, b) result(c)
    !
    ! Compute in parallel
    !                 T
    !    C = A * B * A
    !
    ! with diagonal B. Missing elements of B assumed
    ! to be zeros.
    !
    ! NOTE: it is assumed that input A and B is avaialable
    !       to every worker in full.
    !
    use comm
    implicit none
    real(r8_kind), intent(in)  :: a(:, :) ! (m, >=n)
    real(r8_kind), intent(in)  :: b(:)    ! (n)
    real(r8_kind)              :: c(size(a, 1), size(a, 1)) ! (m, m)
    ! *** end of interface ***

    integer(i4_kind) :: n, np, rank, lo, hi
    integer(i4_kind), allocatable :: counts(:) ! (0:np-1)

    np = comm_size()
    rank = comm_rank()

    ! Use base-0 bounds to index with base-0 ranks:
    allocate(counts(0:np-1))

    !
    ! This will be the parallelization axis:
    !
    n = size(b)
    ASSERT(n<=size(a,2))

    !
    ! FIXME: the length of this axis is the number of occupied levels
    ! (correlates with the number of electrons) in the system.
    ! For small systems this number isnt large.
    !

    ! (almost) equal shares:
    counts = work_shares(n, np)

    !
    ! Compute offsets into b(:) using cumulative sums:
    !
    lo = sum(counts(0:rank-1)) + 1
    hi = sum(counts(0:rank-1)) + counts(rank)

    !
    ! The range b(lo:hi) is to be processed by "rank"
    !
    c = abat(a(:, lo:hi), b(lo:hi))

    !
    ! Sum partial contributions over processors, bcast result:
    !
    call comm_allreduce(c)

    deallocate(counts)
  end function pabat

  function pabct(a, b, c) result(d)
    !
    ! Compute in parallel
    !                 T
    !    D = A * B * C
    !
    ! with diagonal B. Missing elements of B assumed
    ! to be zeros.
    !
    ! NOTE: it is assumed that input A, B, and C is avaialable
    !       to every worker in full.
    !
    use comm
    implicit none
    real(r8_kind), intent(in)  :: a(:, :) ! (m, >=n)
    real(r8_kind), intent(in)  :: b(:)    ! (n)
    real(r8_kind), intent(in)  :: c(:, :) ! (m, >=n)
    real(r8_kind)              :: d(size(a, 1), size(c, 1)) ! (m, m)
    ! *** end of interface ***

    integer(i4_kind) :: n, np, rank, lo, hi
    integer(i4_kind), allocatable :: counts(:) ! (0:np-1)

    np = comm_size()
    rank = comm_rank()

    ! Use base-0 bounds to index with base-0 ranks:
    allocate(counts(0:np-1))

    !
    ! This will be the parallelization axis:
    !
    n = size(b)
    ASSERT(n<=size(a,2))
    ASSERT(n<=size(c,2))

    !
    ! FIXME: the length of this axis is the number of occupied levels
    ! (correlates with the number of electrons) in the system.
    ! For small systems this number isnt large.
    !

    ! (almost) equal shares:
    counts = work_shares(n, np)

    !
    ! Compute offsets into b(:) using cumulative sums:
    !
    lo = sum(counts(0:rank-1)) + 1
    hi = sum(counts(0:rank-1)) + counts(rank)

    !
    ! The range b(lo:hi) is to be processed by "rank"
    !
    d = abct(a(:, lo:hi), b(lo:hi), c(:, lo:hi))

    !
    ! Sum partial contributions over processors, bcast result:
    !
    call comm_allreduce(d)

    deallocate(counts)
  end function pabct
#endif

  function abat(a, b) result(c)
    !
    ! Compute
    !                 T
    !    C = A * B * A
    !
    ! with diagonal B. Missing elements of B assumed
    ! to be zeros.
    !
    use f77_blas, only: dgemm
    implicit none
    real(r8_kind), intent(in)  :: a(:, :) ! (m, >=n)
    real(r8_kind), intent(in)  :: b(:)    ! (n)
    real(r8_kind)              :: c(size(a, 1), size(a, 1)) ! (m, m)
    ! *** end of interface ***

    integer(i4_kind) :: k, n, m

#ifdef FPP_NODGEMM
    integer(i4_kind) :: i, j

    n = size(b)
    ASSERT(n<=size(a, 2))
    m = size(a, 1)

    do j = 1, m
      do i = 1, j
        c(i, j) = 0.0
        do k = 1, n
          c(i, j) = c(i, j) + b(k) * a(i, k) * a(j, k)
        enddo
        c(j, i) = c(i, j)
      enddo
    enddo
#else
    real(r8_kind) :: ab(size(a, 1), size(b))

    n = size(b)
    ASSERT(n<=size(a, 2))
    m = size(a, 1)

    ! dont call BLAS with 0-sized arrays
    if ( n == 0 ) then
       ! c is not necessarily 0-sized:
       c = 0.0
       RETURN
    endif
    if ( m == 0 ) RETURN ! c is 0-sized

    do k = 1, n
      ab(:, k) = a(:, k) * b(k)
    enddo

    call dgemm( 'n', 't', m, m, n, 1.0D0    &
              , ab(:, :), size(ab, 1)       &
              , a(:, :),  size(a, 1), 0.0D0 &
              , c(:, :),  size(c, 1)        )
#endif
  end function abat

  function abct(a, b, c) result(d)
    !
    ! Compute
    !                 T
    !    D = A * B * C
    !
    ! with diagonal B. Missing elements of B assumed
    ! to be zeros.
    !
    use f77_blas, only: dgemm
    implicit none
    real(r8_kind), intent(in)  :: a(:, :) ! (m, >=n)
    real(r8_kind), intent(in)  :: b(:)    ! (n)
    real(r8_kind), intent(in)  :: c(:, :) ! (m, >=n)
    real(r8_kind)              :: d(size(a, 1), size(c, 1)) ! (m, m)
    ! *** end of interface ***

    integer(i4_kind) :: k, n, ma, mc

#ifdef FPP_NODGEMM
    integer(i4_kind) :: i, j

    n = size(b)
    ASSERT(n<=size(a, 2))
    ASSERT(n<=size(c, 2))
    ma = size(a, 1)
    mc = size(c, 1)
    ASSERT(ma==mc)

    do j = 1, mc
      do i = 1, ma
        d(i, j) = 0.0
        do k = 1, n
          d(i, j) = d(i, j) + b(k) * a(i, k) * c(j, k)
        enddo
      enddo
    enddo
#else
    real(r8_kind) :: ab(size(a, 1), size(b))

    n = size(b)
    ASSERT(n<=size(a, 2))
    ASSERT(n<=size(c, 2))
    ma = size(a, 1)
    mc = size(c, 1)
    ASSERT(ma==mc)

    ! dont call BLAS with 0-sized arrays
    if ( n == 0 ) then
       ! d is not necessarily 0-sized:
       d = 0.0
       RETURN
    endif
    if ( ma == 0 ) RETURN ! d is 0-sized

    do k = 1, n
      ab(:, k) = a(:, k) * b(k)
    enddo

    call dgemm( 'n', 't', ma, mc, n, 1.0D0  &
              , ab(:, :), size(ab, 1)       &
              , c(:, :),  size(c, 1), 0.0D0 &
              , d(:, :),  size(d, 1)        )
#endif
  end function abct

  !*************************************************************
  subroutine save_densmat
    ! Purpose: save the density matrix in INPUT directory
    ! so that it can be used later if no nessesity to perform
    ! SCF calculations. At the moment it was realized for not
    ! spin_orbit
    ! Modules used -----------------------------------------------
    use filename_module, only: inpfile
    use iounitadmin_module, only : get_iounit,return_iounit
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i,j,k,io_u
    !------------ Executable code ------------------------------------

    io_u = get_iounit()
    open(unit= io_u, form='unformatted',status='replace', &
         file=trim(inpfile('densmat.save')))

    do i = 1, size(densmat)
       do j = 1, size(densmat(i)%m, 1)
          do k = 1, size(densmat(i)%m, 2)
            write(io_u) densmat(i)%m(j, k, :)
          end do
       end do
    end do

    close(io_u)
    call return_iounit(io_u)
  end subroutine save_densmat


  subroutine open_densmat()
    !
    ! Read the density matrix earlier  saved in INPUT directory At the
    ! moment it was realized for not spin_orbit
    !
    use filename_module, only: inpfile
    use iounitadmin_module, only : openget_iounit, returnclose_iounit
    ! *** end of interface ***

    integer (i4_kind) :: i, j, k, io_u
    logical :: yes

    inquire (file= trim (inpfile ('densmat.save')), exist= yes)
    if (.not. yes) call error_handler &
         ("open_densmat: file densmat.save cannot be read. It is absent")

    io_u = openget_iounit (file= trim (inpfile ('densmat.save')), &
         form= 'unformatted', status= 'old')

    ASSERT(allocated(densmat))
    do i = 1, ssym % n_irrep
       ASSERT(allocated(densmat(i)%m))
       do j = 1, ssym % dim(i)
          do k = 1, ssym % dim(i)
             read (io_u) densmat(i) % m(j, k, 1: ssym % n_spin)
          enddo
       enddo
    enddo

    call returnclose_iounit (io_u)
  end subroutine open_densmat


  subroutine copy2(a, b)
    implicit none
    type(arrmat2), intent(in)  :: a(:)
    type(arrmat2), allocatable :: b(:)
    ! *** end of interface ***

    integer(i4_kind) :: dims(size(a))
    integer(i4_kind) :: i

    do i = 1, size(a)
      dims(i) = size(a(i)%m, 1)
    enddo

    call alloc2(dims, b)

    do i = 1, size(a)
      b(i)%m = a(i)%m
    enddo
  end subroutine copy2

  subroutine copy3(a, b)
    implicit none
    type(arrmat3), intent(in)  :: a(:)
    type(arrmat3), allocatable :: b(:)
    ! *** end of interface ***

    integer(i4_kind) :: dims(size(a))
    integer(i4_kind) :: i, n_spin

    do i = 1, size(a)
      dims(i) = size(a(i)%m, 1)
      n_spin = size(a(i)%m, 3)
    enddo

    call alloc3(dims, n_spin, b)

    do i = 1, size(a)
      b(i)%m = a(i)%m
    enddo
  end subroutine copy3

  subroutine alloc2(dims, dmat)
    implicit none
    integer(i4_kind), intent(in) :: dims(:) ! (n_irrep)
    type(arrmat2), allocatable   :: dmat(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    ASSERT(.not.allocated(dmat))

    allocate(dmat(size(dims)))

    do i = 1, size(dmat)
      allocate(dmat(i)%m(dims(i), dims(i)))
    enddo
  end subroutine alloc2

  subroutine alloc3(dims, n_spin, dmat)
    implicit none
    integer(i4_kind), intent(in) :: dims(:) ! (n_irrep)
    integer(i4_kind), intent(in) :: n_spin
    type(arrmat3), allocatable   :: dmat(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    ASSERT(.not.allocated(dmat))

    allocate(dmat(size(dims)))

    do i = 1, size(dmat)
      allocate(dmat(i)%m(dims(i), dims(i), n_spin))
    enddo
  end subroutine alloc3

  subroutine free2(dmat)
    implicit none
    type(arrmat2), allocatable :: dmat(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    ASSERT(allocated(dmat))

    do i = 1, size(dmat)
      deallocate(dmat(i)%m)
    enddo

    deallocate(dmat)
  end subroutine free2

  subroutine free3(dmat)
    implicit none
    type(arrmat3), allocatable :: dmat(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    ASSERT(allocated(dmat))

    do i = 1, size(dmat)
      deallocate(dmat(i)%m)
    enddo

    deallocate(dmat)
  end subroutine free3

  subroutine maxvals2(a, vals)
    implicit none
    type(arrmat2), intent(in)  :: a(:)
    real(r8_kind), intent(out) :: vals(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    vals = - huge(1.0_r8_kind)

    do i = 1, size(a)
      call maxvals_update_buf(a(i)%m, size(a(i)%m), vals)
    enddo
  end subroutine maxvals2

  subroutine maxvals3(a, vals)
    implicit none
    type(arrmat3), intent(in)  :: a(:)
    real(r8_kind), intent(out) :: vals(:)
    ! *** end of interface ***

    integer(i4_kind) :: i

    vals = - huge(1.0_r8_kind)

    do i = 1, size(a)
      call maxvals_update_buf(a(i)%m, size(a(i)%m), vals)
    enddo
  end subroutine maxvals3

  subroutine maxvals_update_buf(a, n, vals)
    !
    ! Compute "size(vals)" largest values in array "a"
    ! of size "n".
    !
    implicit none
    real(r8_kind), intent(in)    :: a(*) ! (n)
    integer(i4_kind), intent(in) :: n
    real(r8_kind), intent(inout) :: vals(:)
    ! *** end of interface ***

    integer(i4_kind) :: m, i, k

    m = size(vals)
    ASSERT(m>0)

    do i = 1, n
      if ( a(i) <= vals(m) ) cycle

      ! find position k of a(i) in vals(1) >= ... >= vals(m)
      k = 1
      do while ( a(i) < vals(k) ) ! terminates at k==m at the latest
        k = k + 1
      enddo

      ! shift values towards the end, pops vals(m):
      vals(k+1:m) = vals(k:m-1)
      vals(k) = a(i)
    enddo
  end subroutine maxvals_update_buf

  function work_shares(n, np) result(work)
    !
    ! Share N jobs among NP workers
    !
    implicit none
    integer(i4_kind), intent(in) :: n, np
    integer(i4_kind)             :: work(np) ! sum(work) == n
    ! *** end of interface ***

    integer(i4_kind) :: rem

    ! assign a fair fraction:
    work(1:np) = n / np

    ! remainder:
    rem = n - sum(work)
    ASSERT(rem<np)

    ! distribute remainder:
    work(1:rem) = work(1:rem) + 1

    ASSERT(sum(work)==n)
  end function work_shares

!--------------- End of module ----------------------------------
end module density_data_module

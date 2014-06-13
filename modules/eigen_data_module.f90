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
module eigen_data_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose:
  !           Contains the following PUBLIC variables:
  !
  !           (i)  Eigenvalues     EIGVAL
  !           (ii) Eigenvectors    EIGVEC
  !
  !           Contains the following routines:
  !
  !           (i) PUBLIC
  !               - eigen_data_solve: decides which of the
  !                                   solving strategies is used
  !               - eigen_data_alloc: if necessary, called by
  !                                   the above routines
  !
  !          - print_eigendata     -> prints out eigenvectors
  !                                   to the file $TTFSOUT/eigvec.dat
  !
  !        - read_eigenvectors     -> reads in the eigenvectors from
  !                                   a file written by the old lcgto
  !                                   This is needed for comparison to
  !                                   the old lcgto and lateron for the
  !                                   restart mechanism.
  !
  !  Module called by: build_hamiltonian
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification
  ! Author: TG
  ! Date:   9/96
  ! Description:
  ! 1. Parallel solution of the eigen problem if the product of the
  !    number of spins with the number of IRREPS is at least two.
  ! 2. Rewriting of solve_para.
  ! 3. Integrating the routine eigen_solver from eigen_solver.f90
  !    as routine eigensolve_slave.
  !
  ! Modification
  ! Author: HH
  ! Date:   11/97
  ! Description:
  ! 1. New subroutine "send_eigvec_all" for sending the full "eigvec"
  !    information to slaves -> needed by the response module.
  ! 2. New subroutine "receive_eigvec_all" for receiving the "eigvec"
  !    information on the slaves.
  !
  ! Modification
  ! Author: MM
  ! Date:   10/97
  ! Description: extension to spin orbit
   !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined types
  use comm_module  ! comm information
  use iounitadmin_module
  use output_module, only: output_main_scf,output_eigen_strategy
  use hamiltonian_module, only: ham_tot,ham_tot_real,ham_tot_imag, &
                                ham_lsft, lsft_allocstat
  use options_module, only: options_spin_orbit, lvshift_mixing
  implicit none
  private
  save

  !== Interrupt end of public interface of module =================
  integer(kind=i4_kind) :: allocate_stat(5)
  !------------ Declaration of public constants and variables -----
  public :: arrmat2,arrmat3
  type(arrmat2), allocatable, target, public :: eigval(:)
  ! eigval(i_ir)%m(i_orbital,i_spin)
  type(arrmat3),allocatable,target, public   :: eigvec(:)
  ! eigvec(i_ir)%m(i_basis,i_orbital,i_spin)
  ! in case of spin orbit orbital coefficients are complex
  type(arrmat2), allocatable, target, public :: eigvec_real(:)
  type(arrmat2), allocatable, target, public :: eigvec_imag(:)

  type(arrmat2),allocatable,private          :: eigvec_hole_real(:) ! n_irr
  type(arrmat2),allocatable,private          :: eigvec_hole_imag(:) ! n_irr

  !
  ! FIXME: target is here only to make things like
  !
  !     v => eigvec_occ_real(irrep)%m
  !     n => occ_num_occ(irrep)%m
  !
  ! compile also when arrmat2 is declared with
  ! allocatable component instead of pointer.
  !

  !------------ public functions and subroutines ------------------

  public :: eigen_data_solve!(), to be called from a parallel context
  public :: eigen_data_solve1!(), trampoline for master/main_slave, must die
  public :: eigen_data_alloc!(), does no communication
  public :: eigen_data_free!(), does no communication

  public :: &
       & build_lvsft_hamiltonian,&
       & print_eigendata, &
       & eigen_hole_setup,eigen_hole_shutdown,&
       & eigvec_write, eigvec_read

  public :: eigen_data_bcast

  public :: eigen_data_dump!(file), to dump hamiltonian, eigenvectors and eigenvalues to file


!================================================================
! End of public interface of module
!================================================================

! interface rsg
!    SUBROUTINE RSG(NM,N,A,B,W,MATZ,Z,IERR)
!      use type_module
!      INTEGER(kind=i4_kind) N,NM,IERR,MATZ
!      REAL(kind=r8_kind) A(:,:),B(:,:),W(:),Z(:,:)
!    end SUBROUTINE RSG
! end interface


  ! variables needed for HOLES -----------------------------
  logical :: fixed_hole = .false.
  !This variable is a 'mirror'-variable of the fixed_hole-switch
  !in the occupation_module. It is set to .true. in the routine
  !'eigen_hole_setup', which is called only if the fixed_hole-variable
  ! in the occupation_modules is .true. This trick is needed, because
  ! we cannot use the occupation_module here.
  logical :: forced_hole = .false., update_projector=.false.
  ! same as 'fixed_hole'
  integer(kind=i4_kind),allocatable,public :: n_holes_per_irrep(:)
  ! 'n_holes_per_irrep(n_irreps)' contains the number of hole per irrep i
  type(arrmat1int),allocatable,public :: holes_per_irrep(:),&
       spin_of_hole(:)
  ! 'holes_per_irrep(n_irreps)' contains the list of hole indices in irrep i
  ! 'spin_of_hole(n_irreps)' contains the spin of ho,e with indices i_irrep,m_orb
  type(arrmat3),allocatable,public         :: eigvec_hole(:)
  ! eigvec_hole(irrep)%m(i_hole,m_irrep) contains eigenvector of MO
  ! with hole
  integer(kind=i4_kind),parameter  :: N_LIST=5
  ! determines the size of the array where the N_LIST largest dot_products
  ! for the hole localization are stored
  real(kind=r8_kind),parameter     :: hole_tolerance=1.0e-3_r8_kind

  logical, private :: EIGDATA_ALLOCATED = .false.
  logical, private :: EIGDATA_VALID     = .false.

  character(len=*), parameter, private :: EIGDATA_FILENAME  = 'saved_eigdata.dat'

  logical, public :: make_level_shift
  integer(kind=i4_kind), private :: ls_gamma,lsis



!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine eigen_data_solve1()
    !
    ! Separate the legacy parallel context from eigen_data_solve()
    !
    use comm_module, only: comm_init_send, comm_send, comm_all_other_hosts, &
         comm_i_am_master
    use msgtag_module, only: msgtag_eigen_data_solve
    implicit none
    !** End of interface *****************************************

    if( comm_i_am_master() )then
       !
       ! First tell the slaves to enter this subroutine:
       !
       call comm_init_send(comm_all_other_hosts, msgtag_eigen_data_solve)
       call comm_send()
    endif

    call eigen_data_solve()
  end subroutine eigen_data_solve1
  !*************************************************************

  !*************************************************************
  subroutine eigen_data_solve()
    !
    !  Purpose: decide which of the eigensolving strategies to use and
    !           calling aproriate routine.  Runs on all workes.
    !
    !  author: Folke Noertemann
    !  date  : 10/95
    !
    !----------------Modification-------------------------------
    ! Strategy was changed:
    ! Eigenproblems are solved in parallel if
    ! n_irreps*n_spin .gt. 1
    ! TG, 9/96
    !
    use symmetry_data_module, only: ssym  ! symmetry information
    use comm, only: comm_bcast, comm_rank
    use overlap_module, only: overlap
    use time_module, only: start_timer, stop_timer
    use timer_module, only: timer_scf_eigen
#ifdef WITH_SCHEDEIG
    use comm, only: comm_world
    use se_eigen_module, only: se_eigen_compeigs
#endif
#ifdef FPP_DEBUG
    use error_module, only: MyID
#endif
    implicit none
    !** End of interface *****************************************

    TRACE("eigen_data_solve/entered")
    DPRINT MyID, "eigen_data_solve: entered"

    call start_timer(timer_scf_eigen)

    if ( lvshift_mixing ) then
       !
       ! FIXME: this global var is changed ocassionally, possibly on
       !        master only, do this only when changing its value.
       !        Moreover it looks like there is some confusion about
       !        which one of the two variables, lvshift_mixing or
       !        make_level_shift, is to be used in which case.
       !
       call comm_bcast(make_level_shift)
    else
       make_level_shift = .false.
    endif

!   ASSERT(allocated(ham_tot)) ! not with SO
!   ASSERT(allocated(overlap)) ! not on slaves
!   ASSERT(allocated(eigval))
!   ASSERT(allocated(eigvec))

#if defined(WITH_MATRIX_PARALLEL) || defined(WITH_SCHEDEIG)
    !
    ! Parallel eigensolver v2 and v3.  See modules/matrix_parallel.f90
    ! and ./schedeig/
    !

    !
    ! FIXME: does not yet handle spin-orbit case:
    !
    if ( options_spin_orbit ) then
       !
       ! Call the "regular" legacy implementation:
       !
       call solve_legacy(ssym)
    else
#ifdef WITH_SCHEDEIG
       DPRINT "eigen_data_solve: call se_eigen_compeigs(...)"
       call se_eigen_compeigs(ham_tot, overlap, eigval, eigvec, comm_world)
#else
       DPRINT "eigen_data_solve: call geigs_parallel(...)"
       call geigs_parallel(ham_tot, overlap, eigval, eigvec)
#endif
       DPRINT MyID, "eigen_data_solve: done!"
    endif
#else
    !
    ! Traditional solver: wrapper for solve_para_master() and
    ! solve_para_slave().  Runs in a parallel context.  Performs
    ! block-wise parallelization with serial diagonalization of each
    ! block:
    !
    call solve_legacy(ssym)
#endif

    call stop_timer(timer_scf_eigen)
    DPRINT  MyID, "eigen_data_solve: exit"
    TRACE("eigen_data_solve/return")
  end subroutine eigen_data_solve
  !*************************************************************

  subroutine build_lvsft_hamiltonian(n_occo, level_shift, set_start_lvshift)
    !
    ! input:  ham_lsft matrixes to be modifed
    !         level_shift - energy shifts of vacant levels
    !
    ! result: modified ham_lsft for all representations
    !
    use symmetry_data_module, only: ssym
    implicit none
    integer(kind=i4_kind), intent(in), dimension(:,:) :: n_occo
    real(kind=r8_kind), intent(in) :: level_shift
    logical, intent(in) :: set_start_lvshift
    ! *** end of interface ***

    integer(kind=i4_kind) :: i_gamma,is,nm,vac
    real(kind=r8_kind), allocatable, dimension(:,:):: shifts,temp

    make_level_shift = set_start_lvshift

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
     do i_gamma = 1,ssym%n_irrep
       if ( ssym%dim(i_gamma) == 0 ) cycle

       nm = ssym%dim(i_gamma)
       do is = 1, ssym%n_spin
          !             call rsg(nm,nn,ham_tot(i_gamma)%m(1:nn,1:nn,is), &
          !                  overlap(i_gamma)%m,eigval(i_gamma)%m(1:nn,is),matz, &
          !                  eigvec(i_gamma)%m(1:nn,1:nn,is),ierr)
          allocate(shifts(nm, nm), temp(nm, nm), stat=allocate_stat(1))
          ASSERT(allocate_stat(1).eq.0)
          allocate_stat(1) = size(shifts) * 2

          shifts = 0.0_r8_kind

          do vac = n_occo(is, i_gamma) + 1, nm
             shifts(vac, vac) = level_shift
          enddo

          temp = matmul(shifts, transpose(ham_lsft(i_gamma)%m(:, :, is)))

          ham_lsft(i_gamma)%m(1:nm, 1:nm, is) = &
               matmul(ham_lsft(i_gamma)%m(:, :, is), temp)

          deallocate(shifts, temp, stat=allocate_stat(1))
          ASSERT(allocate_stat(1) .eq. 0)
       enddo
    enddo
  end subroutine build_lvsft_hamiltonian

  subroutine sort_irreps(dim_irrep, sorted)
    !  Takes information in ssym and builds the vector
    !  'sorted' from it which contains the IRREP numbers
    !  descendingly ordered by associated matrix dimension.
    !  Sorting method:  Straight injection
    !
    !  Input parameter:
    !  ssym        symmetry information
    !  Output parameter:
    !  sorted      sorted IRREP numbers
    !
    !  Subroutine called by: solve_para_master
    !
    !  TG 9/96
    !
    !------------ Declaration of formal parameters ---------------
    implicit none
    ! dim_irrep : number of independent functions in irrep
    integer(kind=i4_kind), intent(in)  :: dim_irrep(:) ! (n_irreps)
    integer(kind=i4_kind), intent(out) :: sorted(:)    ! (n_irreps)
    ! *** end of interface ***

    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: help,i,j,sj,si
    integer(kind=i4_kind) :: n_irrep ! number of irreps
    !------------ Executable code --------------------------------

    n_irrep = size(dim_irrep)

    do i = 1,n_irrep
       sorted(i) = i
    enddo

    do i=1,n_irrep
       do j=i+1,n_irrep
          sj=sorted(j)
          si=sorted(i)
          if (dim_irrep(sj).gt.dim_irrep(si)) then
             help=sorted(i)
             sorted(i)=sorted(j)
             sorted(j)=help
          endif
       enddo
    enddo
  end subroutine sort_irreps

  subroutine solve_legacy(ssym)
    !  Purpose: provide logic for the solving one spin component
    !           of an IRREP serially on  a slave processor.
    !           This routine is called in the case:
    !           (Number of IRREPS)*(Number of spins) >= 2
    !  Input parameter (not modified on output):
    !  Name:              Description/Range:
    !  ssym               symmetry information
    !  ham_tot            total hamiltonian (type arrmat3)
    !                     including all IRREPs
    !  Subroutine called by: solve_eigen
    !
    !  Modifications:
    !  Totally rewritten
    !  9/96
    !
    ! Context: master only
    !
    !------------ Modules ---------------------------------------
    use symmetry_data_module, only: sym
    use comm, only: comm_rank
    use msgtag_module, only: msgtag_eigensolve, msgtag_SendBackEig, &
         msgtag_stop_daemon
    use spin_orbit_module, only: is_on,op_Eigensolver
    use matrix_module,     only: geigs
    use overlap_module, only: overlap, overlap_real, overlap_imag
    implicit none
    type(sym), intent(in) :: ssym
    !** End of interface *****************************************

    integer(kind=i4_kind) :: n_tasks, n_to_start, n_ready, &
         i_spin, i_gamma, ready_tid, alloc_stat, &
         i_g, i_s, i, non_null, n

    integer(kind=i4_kind), allocatable :: sorted(:) ! (n_irrep)

    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind) :: n_irrep
    integer(kind=i4_kind), allocatable :: dim_irrep(:) ! (n_irrep)
    ! n_irrep    : number of irreps
    ! dim_irrep  : number of independent functions in irrep

    ! MOVED TO BEFORE SCF LOOP: call read_overlap()

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif

    !
    ! Slaves will listen to what master will tell them and exit when
    ! told so.  FIXME: convert it to SPMD style.
    !
    if ( comm_rank() /= 0 ) then
       call slave()
       return ! *** RETURN POINT (SLAVES ONLY) ***
    endif

    if (options_spin_orbit) then
       ASSERT(allocated(eigvec_real))
       ASSERT(allocated(eigvec_imag))
       ASSERT(allocated(overlap_real))
       ASSERT(allocated(overlap_imag))
    else
       ASSERT(allocated(eigvec))
       ASSERT(allocated(overlap))
    endif
       ASSERT(allocated(eigval))

    ! allocate index list of irreps sorted by size
    allocate(sorted(n_irrep),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("SOLVE_LEGACY: allocation of sorted failed")

    ! in "sorted" returns the irrep indices sorted by irrep size:
    call sort_irreps(dim_irrep, sorted)

    ! non-zero sized irreps:
    non_null = count( dim_irrep /= 0 )

    n_tasks = non_null*ssym%n_spin

    !
    ! The tuple (i_gamma, i_spin) identifies the task to be solved.
    ! Note, however, that the actual irrep is taken from the sorted
    ! sequence as sorted(i_gamma):
    !
    i_spin = 1
    i_gamma = 1
    n_to_start = n_tasks

    ! first send starting values to as much slaves as possible
    do i = 2, min0(n_tasks + 1, comm_get_n_processors())
       !
       ! Post the task (i_gamma, i_spin) for solving by worker "i":
       !
       call post(i, sorted(i_gamma), i_spin)

       !
       ! Advance job descriptor: (i_gamma, i_spin):
       !
       n_to_start = n_to_start - 1
       i_spin = i_spin + 1
       if ( i_spin .gt. ssym%n_spin ) then
          i_spin = 1
          i_gamma = i_gamma + 1
       endif
    enddo

    n_ready = 0
    do while ( n_ready .lt. n_tasks )
       !
       ! Check feedback from other workers, if there is none, do one
       ! block myself:
       !
       if ( comm_save_recv_nonblocking(comm_all_other_hosts, msgtag_SendBackEig) ) then ! solve_legacy
          !
          ! Obtain results for intent(out) but unused job descriptor
          ! (i_g, i_s) from intent(out) worker "ready_tid".
          !
          call wait(ready_tid, i_g, i_s)

          if ( n_to_start .gt. 0 ) then ! there are tasks left
             !
             ! Post the task (i_gamma, i_spin) for solving by worker
             ! that just finished:
             !
             call post(ready_tid, sorted(i_gamma), i_spin)
          endif

          !
          ! One more job is ready now:
          !
          n_ready = n_ready + 1

       else ! comm_save_recv_nonblocking

          if ( n_to_start .gt. 0 ) then ! there are tasks left
             !
             ! Solve generalized hermitean-definite eigenvalue
             ! problem H * x = S * x * eigenvalue:
             !
             call solve(sorted(i_gamma), i_spin)

             !
             ! One more job is ready now:
             !
             n_ready = n_ready + 1
          endif
       endif

       !
       ! Advance job descriptor: (i_gamma, i_spin):
       !
       if ( n_to_start .gt. 0 ) then ! there are tasks left
          n_to_start = n_to_start - 1
          i_spin = i_spin + 1
          if ( i_spin .gt. ssym%n_spin ) then
             i_spin = 1
             i_gamma = i_gamma + 1
          endif
       endif
    enddo
    DPRINT 'all tasks done', n_tasks

    ! FIXME: find a better place:
    if(fixed_hole.and..not.forced_hole) then
       ! localize holes before overlap is deallocated:
       call eigen_localize_hole()
    endif

    ! MOVED TO AFTER SCF LOOP: call dealloc_overlap()

    ! finally tell the slaves to exit slave() daemon:
    call comm_init_send(comm_all_other_hosts, msgtag_stop_daemon)
    call comm_send()

    contains

      subroutine solve(irrep, spin)
        !
        ! Diagonalze matrix of irrep "irrep" and spin "spin".
        !
        implicit none
        integer(i4_kind), intent(in) :: irrep, spin
        ! *** end of interface ***

        integer(i4_kind) :: nn, nm, matz, ierr

        if ( options_spin_orbit ) then
           !
           ! SPIN ORBIT
           !
           ASSERT(is_on(op_Eigensolver))
           call geigs(ham_tot_real(irrep)%m, ham_tot_imag(irrep)%m, &
                overlap_real(irrep)%m, overlap_imag(irrep)%m, &
                eigval(irrep)%m(:, 1), &
                eigvec_real(irrep)%m, eigvec_imag(irrep)%m)
        else ! options_spin_orbit
           !
           ! STANDARD SCF (NO SPIN ORBIT)
           !
           nn = dim_irrep(irrep)
           nm = nn
           matz = 1
           if ( lvshift_mixing ) then
              call rsg_lsft(nm, nn, ham_tot(irrep)%m(:, :, spin), & !=>rsg_lsft
                   overlap(irrep)%m, eigval(irrep)%m(:, spin), matz, &
                   eigvec(irrep)%m(:, :, spin), &
                   ham_lsft(irrep)%m(:, :, spin), ierr)
           else
              call geigs(ham_tot(irrep)%m(:, :, spin), overlap(irrep)%m, &
                   eigval(irrep)%m(:, spin), eigvec(irrep)%m(:, :, spin))
           endif
        endif ! options_spin_orbit
      end subroutine solve

      subroutine post(worker, irrep, spin)
        !
        ! Tell worker with base-1 index "worker" to diagonalze matrix
        ! of irrep "irrep" and spin "spin".  See receiver side in
        ! recv().
        !
        implicit none
        integer(i4_kind), intent(in) :: worker, irrep, spin
        ! *** end of interface ***

        integer(i4_kind) :: info

        call comm_init_send(worker, msgtag_eigensolve)

        if ( options_spin_orbit ) then
           !
           ! SPIN ORBIT
           !
           call commpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call commpack(ham_tot_real(irrep)%m(1, 1), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)

           call commpack(ham_tot_imag(irrep)%m(1, 1), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)
        else ! options_spin_orbit
           !
           ! STANDARD SCF (NO SPIN ORBIT)
           !
           call commpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call commpack(spin, 1, 1, info)
           ASSERT(info==0)

           call commpack(ham_tot(irrep)%m(1, 1, spin), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)

           if ( make_level_shift ) then
              call commpack(ham_lsft(irrep)%m(1, 1, spin), dim_irrep(irrep)**2, 1, info)
              ASSERT(info==0)
           endif
        endif ! options_spin_orbit

        call comm_send()
      end subroutine post

      subroutine recv(irrep, spin, ham_real, ham_imag)
        !
        ! Return input sent by the master.  We declare array
        ! allocatable as their shape is only know after irrep index is
        ! received.  See the sender side in post().
        !
        use msgtag_module, only: msgtag_SendBackEig
        implicit none
        integer(i4_kind), intent(out) :: irrep, spin
        real(r8_kind), intent(out), allocatable :: ham_real(:, :)
        real(r8_kind), intent(out), allocatable :: ham_imag(:, :)
        ! *** end of interface ***

        integer(i4_kind) :: n, info

        if ( options_spin_orbit ) then
           call communpack(irrep, 1, 1, info)
           ASSERT(info==0)

           ! this is not sent in SO case:
           spin = 1

           n = dim_irrep(irrep)
           allocate(ham_real(n, n), ham_imag(n, n))

           call communpack(ham_real(1, 1), size(ham_real), 1, info)
           ASSERT(info==0)

           call communpack(ham_imag(1, 1), size(ham_imag), 1, info)
           ASSERT(info==0)
        else
           call communpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call communpack(spin, 1, 1, info)
           ASSERT(info==0)

           n = dim_irrep(irrep)
           allocate(ham_real(n, n))

           call communpack(ham_real(1, 1), size(ham_real), 1, info)
           ASSERT(info==0)

           if ( make_level_shift ) then
              ! the name of the variable is misleading:
              allocate(ham_imag(n, n))
              call communpack(ham_imag(1, 1), size(ham_imag), 1, info)
              ASSERT(info==0)
           endif
        endif
      end subroutine recv

      subroutine wait(worker, irrep, spin)
        !
        ! Receive results and return base-1 index "worker", "irrep"
        ! and "spin" that identify them.  See the sender part in
        ! reply().
        !
        implicit none
        integer(i4_kind), intent(out) :: worker, irrep, spin
        ! *** end of interface ***

        integer(i4_kind) :: info

        if (options_spin_orbit) then
           !
           ! SPIN ORBIT
           !
           call communpack(irrep, 1, 1, info)
           ASSERT(info==0)

           ! not used in SO:
           spin = -1

           call communpack(eigvec_real(irrep)%m(1, 1), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)

           call communpack(eigvec_imag(irrep)%m(1, 1), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)

           call communpack(eigval(irrep)%m(1, 1), dim_irrep(irrep), 1, info)
           ASSERT(info==0)
        else ! options_spin_orbit
           !
           ! STANDARD SCF (NO SPIN ORBIT)
           !
           call communpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call communpack(spin, 1, 1, info)
           ASSERT(info==0)

           call communpack(eigvec(irrep)%m(1, 1, spin), dim_irrep(irrep)**2, 1, info)
           ASSERT(info==0)

           call communpack(eigval(irrep)%m(1, spin), dim_irrep(irrep), 1, info)
           ASSERT(info==0)

           if ( lvshift_mixing ) then
              call communpack(ham_lsft(irrep)%m(1, 1, spin), dim_irrep(irrep)**2, 1, info)
              ASSERT(info==0)
           endif
        endif ! options_spin_orbit

        worker = comm_sendinghost()
      end subroutine wait

      subroutine reply(irrep, spin, eigval, eigvec_real, eigvec_imag)
        !
        ! Send results back to master.  We declare array allocatable
        ! as not all of the are always allocated. See the receiver
        ! side in wait().
        !
        implicit none
        integer(i4_kind), intent(in) :: irrep, spin
        real(r8_kind), intent(in), allocatable :: eigval(:)
        real(r8_kind), intent(in), allocatable :: eigvec_real(:, :)
        real(r8_kind), intent(in), allocatable :: eigvec_imag(:, :)
        ! *** end of interface ***

        integer(i4_kind) :: info

        call comm_init_send(comm_master_host, msgtag_SendBackEig)

        if ( options_spin_orbit ) then
           call commpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call commpack(eigvec_real(1, 1), size(eigvec_real), 1, info)
           ASSERT(info==0)

           call commpack(eigvec_imag(1, 1), size(eigvec_imag), 1, info)
           ASSERT(info==0)

           call commpack(eigval(1), size(eigval), 1, info)
           ASSERT(info==0)
        else
           call commpack(irrep, 1, 1, info)
           ASSERT(info==0)

           call commpack(spin, 1, 1, info)
           ASSERT(info==0)

           call commpack(eigvec_real(1, 1), size(eigvec_real), 1, info)
           ASSERT(info==0)

           call commpack(eigval(1), size(eigval), 1, info)
           ASSERT(info==0)

           if ( lvshift_mixing ) then
              ! the name of the variable is misleading:
              call commpack(eigvec_imag(1, 1), size(eigvec_imag), 1, info)
              ASSERT(info==0)
           endif
        endif

        call comm_send()
      end subroutine reply

      subroutine slave()
        !
        ! Slave daemon that proceses the orders of master,
        ! was part of main_slave()
        !
        ! Context: slaves only
        !
        use error_module, only: MyID
        use comm_module
        use msgtag_module, only: msgtag_eigensolve, msgtag_stop_daemon, msgtag_name
        implicit none
        ! *** end of interface ***

        integer(i4_kind) :: msgtag

        ! MOVED TO prescf_init(): call read_overlap()

        do ! receive and handle messages until told to exit:

           call comm_save_recv(comm_master_host, comm_any_message)

           msgtag = comm_msgtag()
           DPRINT  MyID//'solve_para_slave: received msgtag=',msgtag_name(msgtag)

           select case ( msgtag )

           case (msgtag_eigensolve)

              DPRINT  MyID//'solve_para_slave: call solve_slave()'
              ! diagonalize one block and send results back:
              call solve_slave()

           case (msgtag_stop_daemon)
              DPRINT  MyID//'solve_para_slave: terminating'
              exit ! do loop

           case default
              print *,MyID//'solve_para_slave:: received msgtag=',msgtag,': ',msgtag_name(msgtag)
              call error_handler('solve_para_slave:: wrong message tag: '//msgtag_name(msgtag))
           end select
        enddo

        ! MOVED TO prescf_finalize(): call dealloc_overlap()
      end subroutine slave

      subroutine solve_slave()
        !-------------- Module specification -----------------------
        !
        !  Purpose: put one spin component of an IRREP of the total
        !           hamiltonian into EISPACKS eigensolver.
        !           IRREP is received via comm and Eigenvalues
        !           and Eigenvectors are sent back.
        !           Due to practical considerations, the overlap matrix
        !           is simply read in again. Therefore the complete
        !           symmetry information has to be present.
        !           Additionally the number of the actual IRREP
        !           has to be present. This piece of information also has to be
        !           sent back to enable the master to receive the eigenstuff
        !           directly in eigvec(i_gamma)%m etc.
        !
        !  Context: slaves only
        !
        !  Author: Folke Noertemann
        !  Date: 8/95
        !
        !---------------------------------------------------------------
        ! Modifications
        ! Modification: This subroutine was in an own file in older
        !               versions.
        !               Solving the eigen problem only for one
        !               spin component of an IRREP.
        !               Minor modifications were necessary.
        ! Author: TG
        ! Date: 9/96
        !** End of interface *****************************************
        !--------------------------------------------------------------
        ! MODULES
        !--------------------------------------------------------------
        use overlap_module, only: overlap, overlap_real, overlap_imag, &
            read_overlap, dealloc_overlap
        use matrix_module, only: geigs
        implicit none

        !------------ Declaration of local variables --------------------
        integer(kind=i4_kind) :: i_gamma, i_spin, n

        !---- help variables for eigensolver
        real(r8_kind), allocatable :: a(:,:), z(:,:), a_imag(:,:), z_imag(:,:)
        real(r8_kind), allocatable :: w(:)
        integer(i4_kind) :: matz, ierr

        !
        ! Receive data from master:
        !
        call recv(i_gamma, i_spin, a, a_imag)

        ! dimension of output:
        n = dim_irrep(i_gamma)

#ifdef WITH_BGY3D
        if (.not.allocated(overlap)) then
           call read_overlap()
        endif
#endif

        if (options_spin_orbit) then
           !
           ! SPIN ORBIT
           !

           ! allocate and read overlap matrix
           ! MOVED TO BEFORE SCF LOOP: call read_overlap()
           ASSERT(allocated(overlap_real))
           ASSERT(allocated(overlap_imag))

           allocate(w(n), z(n, n), z_imag(n, n))

           !
           ! Solve generalized hermitean-definite eigenvalue problem:
           !
           call geigs(a, a_imag, overlap_real(i_gamma)%m, overlap_imag(i_gamma)%m, w, z, z_imag)
        else ! options_spin_orbit
           !
           ! STANDARD SCF (NO SPIN ORBIT)
           !

           ! allocate and read overlap matrix
           ! MOVED TO BEFORE SCF LOOP: call read_overlap()
           ASSERT(allocated(overlap))

           allocate(w(n), z(n, n))

           ! here the eigensolver is called
           matz = 1 ! Yes, we also want the eigenvectors
           if ( lvshift_mixing ) then
              ! FIXME: make sure rsg_lsft does not modify overlap:
              call rsg_lsft(n, n, a, overlap(i_gamma)%m, w, matz, z, z_imag, ierr)
           else
              call geigs(a, overlap(i_gamma)%m, w, z)
           endif
        endif ! options_spin_orbit

        !
        ! Send back i_gamma, i_spin, eigenvectors and eigenvalues.
        ! Spin is ignored in SO case. The last argument (ham_lsftM) is
        ! ignored if .not. make_level_shift in the regular (no SO)
        ! case:
        !
        call reply(i_gamma, i_spin, w, z, z_imag)

        ! MOVED TO AFTER SCF LOOP: call dealloc_overlap()
#ifdef WITH_BGY3D
        call dealloc_overlap()
#endif
      end subroutine solve_slave

  end subroutine solve_legacy


  subroutine eigen_data_alloc()
    !
    ! Purpose: allocate the appropriate space for the eigenstuff
    !
    ! FN 10/95
    !
    use symmetry_data_module
    use init_module, only : init
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind)         :: i_gamma,alloc_stat,n

    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep : number of independent functions in irrep
     !------------ Executable code ------------------------------

    DPRINT 'eigen_data_alloc: entered'

    ! Set  to  true  at  the  end  of this  sub,  reset  to  false  by
    ! eigen_data_free():
    ASSERT(.not.EIGDATA_ALLOCATED)

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       allocate (eigvec_real(n_irrep),eigvec_imag(n_irrep),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       allocate (eigval(n_irrep),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       do i_gamma=1,n_irrep
          allocate( eigvec_real(i_gamma)%m(dim_irrep(i_gamma), &
               dim_irrep(i_gamma)),eigvec_imag(i_gamma)%m(dim_irrep(i_gamma),&
               dim_irrep(i_gamma)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)

          allocate( eigval(i_gamma)%m(dim_irrep(i_gamma), &
               1),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo

       call init(eigvec_real)
       call init(eigvec_imag)

    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       allocate (eigvec(ssym%n_irrep),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       allocate (eigval(ssym%n_irrep),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       do i_gamma=1,ssym%n_irrep
          allocate( eigvec(i_gamma)%m(ssym%dim(i_gamma), &
               ssym%dim(i_gamma),ssym%n_spin),STAT=alloc_stat)
          ASSERT(alloc_stat==0)

          allocate( eigval(i_gamma)%m(ssym%dim(i_gamma), &
               ssym%n_spin),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
       call init(eigvec)
    endif! options_spin_orbit

    call init(eigval)

    ! deallocate appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    deallocate(dim_irrep)

    EIGDATA_ALLOCATED = .true.
  end subroutine eigen_data_alloc


  subroutine eigen_data_free()
    !
    ! Purpose:  complimentary to eigen_data_alloc().   Deallocates the
    ! eigenvectors where necessary.
    !
    ! Subroutine called by: finalize_geometry() at least.
    !
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind) :: i, alloc_stat

    DPRINT 'eigen_data_free: entered'

    ASSERT(EIGDATA_ALLOCATED)

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i = 1, size(eigval)
          deallocate(eigvec_real(i)%m, eigvec_imag(i)%m, STAT=alloc_stat)
          ASSERT(alloc_stat==0)

          deallocate(eigval(i)%m ,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo

       deallocate(eigvec_real, eigvec_imag, STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       deallocate(eigval, STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    else ! options_spin_orbit
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       do i = 1, size(eigval)
          deallocate(eigvec(i)%m, STAT=alloc_stat)
          ASSERT(alloc_stat==0)

          deallocate(eigval(i)%m, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo

       deallocate(eigvec, STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       deallocate(eigval, STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    endif! options_spin_orbit

    EIGDATA_ALLOCATED = .false.
  end subroutine eigen_data_free

  subroutine eigen_data_dump(filename)
    !
    ! Dump Hamiltonian, eigenvectors and eigenvalues into a file
    ! for debugging only.
    !
    use iounitadmin_module
    use filename_module, only: outfile
    use overlap_module, only: overlap
    implicit none
    character(len=*), intent(in) :: filename
    ! *** end of interface ***

    integer(kind=i4_kind)  :: io_u
    character(len=6), save :: position = "asis"
    !------------ Executable code ------------------------------

    io_u = openget_iounit( file=trim(outfile(filename)) &
                         , form='formatted', status='unknown', position=position )

    ! on subsequend calls append:
    position = "append"

    ! for SO need to dump real/imag:
    ASSERT(.not.options_spin_orbit)

    ASSERT(allocated(ham_tot))
    ASSERT(allocated(overlap))
    ASSERT(allocated(eigvec))
    ASSERT(allocated(eigval))

    call data_dump(io_u, ham_tot, overlap, eigvec, eigval)

    call returnclose_iounit(io_u)
  end subroutine eigen_data_dump

  subroutine data_dump(iou, h, s, v, e)
    !
    ! Dump data to unit "iou"
    !
    use debug, only: dump
    implicit none
    integer(i4_kind), intent(in) :: iou
    type(arrmat3), intent(in)    :: h(:) ! (n_irreps)%m(:,:,:), hamiltonian
    type(arrmat2), intent(in)    :: s(:) ! (n_irreps)%m(:,:), overlap
    type(arrmat3), intent(in)    :: v(:) ! (n_irreps)%m(:,:,:), eigenvectors
    type(arrmat2), intent(in)    :: e(:) ! (n_irreps)%m(:,:), eigenvalues
    ! *** end of interface ***

    integer(I4_kind) :: irr

    call dump(iou,size(h))
    do irr=1,size(h)
      call dump(iou,h(irr)%m)
    enddo

    call dump(iou,size(s))
    do irr=1,size(s)
      call dump(iou,s(irr)%m)
    enddo

    call dump(iou,size(v))
    do irr=1,size(v)
      call dump(iou,v(irr)%m)
    enddo

    call dump(iou,size(e))
    do irr=1,size(e)
      call dump(iou,e(irr)%m)
    enddo
  end subroutine data_dump

  !*************************************************************
  subroutine print_eigendata()
    !  Purpose: print out the eigenvectors in a well-suited
    !           format to the file $TTFSOUT/eigemist.
    !           The format is supposed to be comparable to
    !           the eigenvectors of the old lcgto
    !           The eigenvectors are printed out as LINES.
    !  subroutine called by: 'main_scf'
    !** End of interface ***************************************
    !  Author: Folke Noertemann
    !  Date: 8/95
    !
    !------------ Modules --------------------------------------
    use symmetry_data_module
    use iounitadmin_module
    use filename_module, only: outfile
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: i_gamma,io_u,i_dim,is,m,mm,count
    data count / 0 /
    ! number and dimensions of irreps (spin orbit and vector case)
    integer(kind=i4_kind)               :: n_irrep
    integer(kind=i4_kind),allocatable   :: dim_irrep(:)
    !------------ Executable code ------------------------------

    io_u=get_iounit()
    open(io_u,form='formatted',status='unknown',&
         position='append',file=&
         trim(outfile('eigvec.dat')))

    if (options_spin_orbit) then
       n_irrep = symmetry_data_n_proj_irreps()
       allocate(dim_irrep(n_irrep))
       do m=1,n_irrep
          dim_irrep(m) = symmetry_data_dimension_proj(m)
       enddo
    else
       n_irrep = symmetry_data_n_irreps()
       allocate(dim_irrep(n_irrep))
       do m=1,n_irrep
          dim_irrep(m) = symmetry_data_dimension(m)
       enddo
    endif


    count=count+1
    write(io_u,*)' '
    write(io_u,*)' '
    write(io_u,*)' '
!    write(io_u,*)'+++++++++++++ Loop ',count,' ++++++++++++++++'
    do i_gamma=1,n_irrep
       i_dim=dim_irrep(i_gamma)
       write(io_u,*)' '
       write(io_u,*)' '
       write(io_u,*)'---------- Irrep ',i_gamma,' -------------'

       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          write(io_u,*)' '
          do m=1,i_dim
             write(io_u,*)' Eigenvalue (',m,') :',eigval(i_gamma)%m(m,1)
             write(io_u,*)' real part:'
             write(io_u,'(4(A,E19.12))')&
                  ( ' ',eigvec_real(i_gamma)%m(mm,m),mm=1,i_dim)
             write(io_u,*)' imaginary part:'
             write(io_u,'(4(A,E19.12))')&
                  ( ' ',eigvec_imag(i_gamma)%m(mm,m),mm=1,i_dim)
          enddo
       else ! options_spin_orbit
          !
          ! STANDARD SCF (NO SPIN ORBIT)
          !
          do is=1,symmetry_data_n_spin()
             write(io_u,*)' Spin : ',is
             write(io_u,*)' '


             do m=1,i_dim
                write(io_u,*)' '
                write(io_u,*)' Eigenvalue :',eigval(i_gamma)%m(m,is)
                write(io_u,*)' '
                write(io_u,'(4(A,E19.12))')&
                     ( ' ',eigvec(i_gamma)%m(mm,m,is),mm=1,i_dim)
             enddo

          enddo
       endif! options_spin_orbit
    enddo
    write(io_u,*)' '
    !write(io_u,*)'+++++++++++ End of loop ',count,' +++++++++++++'
    close(io_u)
    call return_iounit(io_u)
    deallocate(dim_irrep)
  end subroutine print_eigendata

  subroutine eigvec_write()
    use iounitadmin_module, only: openget_iounit, returnclose_iounit
    use filename_module, only: outfile
    implicit none
    ! *** end of interface ***
    integer(i4_kind) :: irr,n_irr,unit,sh(2)

    DPRINT 'eid::eigvec_write: entered'

    ASSERT(options_spin_orbit)

    n_irr = size(eigvec_real)
    ASSERT(n_irr==size(eigvec_imag))

    unit = openget_iounit(&
         & trim(outfile(EIGDATA_FILENAME)), &
         & form='unformatted')

    write (unit) n_irr
    do irr=1,n_irr
       sh = shape(eigvec_real(irr)%m)
       write (unit) sh
       write(unit) eigvec_real(irr)%m
       write(unit) eigvec_imag(irr)%m
    enddo

    call returnclose_iounit(unit)
  end subroutine eigvec_write

  subroutine eigvec_read()
    use iounitadmin_module, only: openget_iounit, returnclose_iounit
    use filename_module, only: recfile
    implicit none
    ! *** end of interface ***
    integer(i4_kind) :: irr,n_irr,unit
    integer(i4_kind) :: sh1(2),sh2(2)

    DPRINT 'eid::eigvec_read: entered'

    ASSERT(options_spin_orbit)
    ASSERT(EIGDATA_ALLOCATED)
    ASSERT(allocated(eigvec_real))
    ASSERT(allocated(eigvec_imag))

    unit = openget_iounit (recfile (EIGDATA_FILENAME), form='unformatted')

    read (unit) n_irr
    ASSERT(n_irr==size(eigvec_real))
    ASSERT(n_irr==size(eigvec_imag))

    do irr=1,n_irr
       sh1 = shape(eigvec_real(irr)%m)
       read (unit) sh2
       ASSERT(all(sh1==sh2))

       read(unit) eigvec_real(irr)%m
       read(unit) eigvec_imag(irr)%m
    enddo

    EIGDATA_VALID = .true.

    call returnclose_iounit(unit)
  end subroutine eigvec_read

  !*************************************************************
  subroutine eigen_hole_setup (force_hole, n_holes, hole_list, &
       hole_irrep, force_update, hole_spin)
    !
    ! Purpose: fill the variables 'n_holes_per_irrep',
    ! 'holes_per_irrep' and eigvec_hole using the variables
    ! fixed_hole, hole_list and hole_irrep from the occupation
    ! module.
    !
    ! Subroutine called by: 'occupation_get_holes', which in turn
    ! is called by 'do_recover' in 'main_scf' to ensure that
    ! a) fixed_hole modus works only with recovered eigenvectors
    ! b) eigenvectors are allocated and initialized
    !-------------------------------------------------------
    use symmetry_data_module
    implicit none
    ! --- Declaration of formal parameters -----------------
    logical,intent(in)               :: force_hole
    integer(kind=i4_kind),intent(in) :: n_holes
    integer(kind=i4_kind),intent(in) :: hole_list(:)
    integer(kind=i4_kind),intent(in) :: hole_irrep(:)
    integer(i4_kind), intent(in), optional :: hole_spin(:)
    logical,intent(in)               :: force_update
    !** End of interface ***************************************
    ! --- Declaration of local variables -------------------
    integer(kind=i4_kind) :: counter,alloc_stat,i,j,num_hole&
         ,num_spin,i_hole
    logical :: spin_restricted
    integer(i4_kind) :: irr,n_irr,n_func,n_hole,hole
    !---- Executable code ----------------------------------

    if (present(hole_spin)) then
       spin_restricted=.false.
    else
       spin_restricted=.true.
    endif

    if(force_hole) then
       forced_hole=.true.
    else
       forced_hole=.false.
    endif

    fixed_hole=.true.
    ! This variable determines if the projector, i.e. the
    ! hole eigenvector originating from the reference
    ! calculation is to be updated during the SCF-Cycles
    ! even if the hole does not *jump* from one orbital
    ! to another. This can be set in order to take care of
    ! cases where the holes moves slowly and then suddenly
    ! cannot be idetnified anymore.
    update_projector=force_update

    if(options_spin_orbit)then
       n_irr = symmetry_data_n_proj_irreps()
    else
       n_irr = symmetry_data_n_irreps()
    endif

    if(options_spin_orbit.and..not.EIGDATA_VALID)then
       DPRINT 'eid::eigen_hole_setup: call eigvec_read()'
       call eigvec_read()
       DPRINT 'eid::eigen_hole_setup: .'
    endif

    ! allocate the array 'n_holes_per_irrep' for later use
    allocate(n_holes_per_irrep(n_irr),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("eigen_hole_setup: allocation of n_holes_per_irrep failed")
    n_holes_per_irrep = 0
    allocate(holes_per_irrep(n_irr),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("eigen_hole_setup: allocation of holes_per_irrep failed")
    allocate(spin_of_hole(n_irr),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("eigen_hole_setup: allocation of spin_of_hole failed")

    n_holes_per_irrep = 0
    do i=1,n_holes
       n_holes_per_irrep(hole_irrep(i)) = &
            n_holes_per_irrep(hole_irrep(i)) + 1
    enddo

    irreps: do i=1,n_irr !symmetry_data_n_irreps()
!!$       ! zero-length arrays are allowed:
!!$       if (n_holes_per_irrep(i) == 0) cycle irreps

       allocate(holes_per_irrep(i)%m(n_holes_per_irrep(i)),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("eigen_hole_setup: allocation (1a) failed")

       allocate(spin_of_hole(i)%m(n_holes_per_irrep(i)),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("eigen_hole_setup: allocation (2a) failed")

       counter=0
       do j=1,n_holes
          if (hole_irrep(j)==i) then
             counter=counter+1
             holes_per_irrep(i)%m(counter)=hole_list(j)
             if (present(hole_spin)) then
                spin_of_hole(i)%m(counter) = hole_spin(j)
             else
                spin_of_hole(i)%m(counter) = 1
             endif
          endif
       enddo
    enddo irreps

    if(.not.forced_hole) then
       if(.not.options_spin_orbit)then
          allocate(eigvec_hole(n_irr),STAT=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ("eigen_hole_setup: allocation (1) failed")!
          irreps_2:  do i=1,n_irr !symmetry_data_n_irreps()
             allocate(eigvec_hole(i)%m(symmetry_data_dimension(i),&
                  n_holes_per_irrep(i),symmetry_data_n_spin()),STAT=alloc_stat)
             if(alloc_stat/=0) call error_handler &
                  ("eigen_hole_setup: allocation (2) failed")

!!$             ! zero-length loops are allowed:
!!$             if (n_holes_per_irrep(i) == 0) cycle irreps_2

             do i_hole=1,n_holes_per_irrep(i)
                num_hole=holes_per_irrep(i)%m(i_hole)
                if(spin_restricted) then
                   eigvec_hole(i)%m(:,i_hole,1) = eigvec(i)%m(:,num_hole,1)
                else
                   num_spin=spin_of_hole(i)%m(i_hole)
                   eigvec_hole(i)%m(:,i_hole,num_spin) = &
                        eigvec(i)%m(:,num_hole,num_spin)
                endif
             enddo
             ! now eigvec_hole(irrep)%m(basis,orbital_with_hole,spin)
             ! contains the eigenvector for the orbitals specified
             ! to have a hole.
          enddo irreps_2
       else
          DPRINT 'eid::eigen_hole_setup: spin-orbit branch...'
          ASSERT(EIGDATA_VALID)
!!$          DPRINT 'eid::eigen_hole_setup: n_irr=',n_irr,&
!!$               & 'size(n_holes_per_irrep)=',size(n_holes_per_irrep)
!!$          DPRINT 'eid::eigen_hole_setup: n_holes_per_irrep=',n_holes_per_irrep
          allocate(eigvec_hole_real(n_irr),eigvec_hole_imag(n_irr),&
               & STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          do irr=1,n_irr
             n_func = symmetry_data_dimension_proj(irr)
             n_hole = n_holes_per_irrep(irr)
!!$             DPRINT 'eid::eigen_hole_setup: n_func=',n_func,'n_hole=',n_hole,&
!!$                  & 'size(holes_per_irrep(irr)%m=',size(holes_per_irrep(irr)%m)
             allocate(&
                  & eigvec_hole_real(irr)%m(n_func,n_hole),&
                  & eigvec_hole_imag(irr)%m(n_func,n_hole),&
                  & STAT= alloc_stat)
             ASSERT(alloc_stat==0)
             do i_hole=1,n_hole
                hole = holes_per_irrep(irr)%m(i_hole)
                DPRINT 'eid::eigen_hole_setup: store hole(irr,i,hole)=',irr,i_hole,hole
                eigvec_hole_real(irr)%m(:,i_hole) = eigvec_real(irr)%m(:,hole)
                eigvec_hole_imag(irr)%m(:,i_hole) = eigvec_imag(irr)%m(:,hole)
             enddo
          enddo
       endif
    endif

  end subroutine eigen_hole_setup
  !*************************************************************


  !*************************************************************
  subroutine eigen_hole_shutdown()
    ! Purpose: deallocate all variables that are necessary for
    ! the fixed_hole modus, i.e. all book-keeping variables as
    ! well es the data-structures to keep those eigenvectors
    ! belonging to holes.
    !
    ! Subroutine called by: 'main_scf'.
    !** End of interface ***************************************
    !-------------------------------------------------------
    integer(kind=i4_kind) :: i,alloc_stat,n_irreps
    ! --- Executable code ----------------------------------

!!$    n_irreps=symmetry_data_n_irreps()
    n_irreps = size(n_holes_per_irrep)

    do i=1,n_irreps
!!$       ! zero-length arrays
!!$       if (n_holes_per_irrep(i) == 0) cycle
       deallocate(holes_per_irrep(i)%m,STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("eigen_hole_shutdown: deallocation (2) failed")

       deallocate(spin_of_hole(i)%m,STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("eigen_hole_shutdown: deallocation (3) failed")
    enddo

    deallocate(n_holes_per_irrep,STAT=alloc_stat)
    if(alloc_stat/=0) call error_handler&
         ("eigen_hole_shutdown: deallocation (1) failed")


    deallocate(holes_per_irrep,STAT=alloc_stat)
    if(alloc_stat/=0) call error_handler&
         ("eigen_hole_shutdown: deallocation (4) failed")
    deallocate(spin_of_hole,STAT=alloc_stat)
    if(alloc_stat/=0) call error_handler&
         ("eigen_hole_shutdown: deallocation (5) failed")

    if (.not.forced_hole) then

       if(allocated(eigvec_hole))then
          !.not.options_spin_orbit:
          do i=1,n_irreps
             deallocate(eigvec_hole(i)%m,STAT=alloc_stat)
             if(alloc_stat/=0) call error_handler&
                  ("eigen_hole_shutdown: deallocation (6) failed")
          enddo
          deallocate(eigvec_hole,STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler&
               ("eigen_hole_shutdown: deallocation (7) failed")
       endif

       if(allocated(eigvec_hole_real))then
          ! options_spin_orbit:
          do i=1,n_irreps
             deallocate(eigvec_hole_real(i)%m,STAT=alloc_stat)
             ASSERT(alloc_stat==0)
             deallocate(eigvec_hole_imag(i)%m,STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(eigvec_hole_real,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(eigvec_hole_imag,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif
    endif

  end subroutine eigen_hole_shutdown
  !*************************************************************


  !*************************************************************
  subroutine eigen_localize_hole()
    ! Purpose: try to localize the user specified holes by biulding
    ! the dot_product:
    ! c = sum(i,j){ c_hole(i)<Chi_i|Chi_j>c_act(j)}
    ! where
    ! c_hole         is the eigenvector of the user-specified hole
    ! c_act          is the actual eigenvector
    ! <Chi_i|Chi_j>  is the overlap matrix
    ! The actual eigenvector with the largest dot_product c is
    ! identified as the actual hole-Orbital. If this is a different
    ! orbital than the user specified, a warning is printed out.
    !
    ! Subroutine called by: solve_serial, solve_para
    use symmetry_data_module
    use overlap_module, only: overlap, overlap_real, overlap_imag
    implicit none
    !** End of interface *****************************************

    integer(i4_kind) :: i, n_holes, alloc_stat, i_hole, &
         index_largest, j, is, i_list, i1, i2, help_index
    real(r8_kind) :: c, c_max, diff, help
    real(r8_kind), allocatable :: help_hole(:)
    logical :: spin_restricted
    real(r8_kind), dimension(N_LIST) :: c_list, help_list
    integer(i4_kind), dimension(N_LIST) :: index_list, help_i_list
    integer(i4_kind), allocatable :: dims(:) ! n_irr
    real(r8_kind), allocatable :: help_hole_real(:), help_hole_imag(:)
    integer(i4_kind) :: irr, n_irr

    DPRINT 'eid::eigen_localize_hole: entered'

    n_holes = sum(n_holes_per_irrep)
    if(symmetry_data_n_spin() == 2) then
       spin_restricted=.false.
    else
       spin_restricted=.true.
    endif

    if(options_spin_orbit)then
       n_irr = symmetry_data_n_proj_irreps()
    else
       n_irr = symmetry_data_n_irreps()
    endif
    allocate(dims(n_irr),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    do irr=1,n_irr
       if(options_spin_orbit)then
          dims(irr) = symmetry_data_dimension_proj(irr)
       else
          dims(irr) = symmetry_data_dimension(irr)
       endif
    enddo

    irreps: do i=1,n_irr !symmetry_data_n_irreps()
       if (n_holes_per_irrep(i) == 0) cycle irreps

       if(options_spin_orbit)then
          allocate(help_hole_real(dims(i)),help_hole_imag(dims(i)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       else
          allocate(help_hole(dims(i)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif

       holes: do i_hole=1,n_holes_per_irrep(i)
          if(spin_restricted) then
             is=1
          else
             is=spin_of_hole(i)%m(i_hole)
          endif

!!$          do j=1,symmetry_data_dimension(i)
!!$             help_hole(j) = sum(eigvec_hole(i)%m(:,i_hole,is)*overlap(i)%m(j,:))
!!$          enddo
          call pre_compute(i,i_hole,is) ! see contains

          c_max=0.0_r8_kind
          c_list=-huge(1.0_r8_kind)
          index_largest=0
          do j=1,dims(i) !symmetry_data_dimension(i) ! Loop over MO`s in Irrep i

             !c=sum(help_hole(:)*eigvec(i)%m(:,j,is))
             c = scalar_product(i,j,is) ! with current hole, reterns abs
!!$             DPRINT 'eid::eigen_localize_hole: c(',i_hole,',',j,')=',c

!!$             ! We want to compare only the absolute value of the dot_product
!!$             if (c<0.0_r8_kind) then
!!$                c=-c
!!$             endif
             if (minval(c_list)==-huge(1.0_r8_kind)) then
                fill: do i_list=1,N_LIST ! first fill 'c_list' completely and if it is filled
                   if(c_list(i_list)==-huge(1.0_r8_kind)) then ! sort its entries
                      c_list(i_list) = c
                      index_list(i_list) = j
                      if(i_list==N_LIST) then
                         ! sort c_list
                         do i1=1,N_LIST
                            do i2=i1+1,N_LIST
                               if (c_list(i2)>c_list(i1)) then
                                  help=c_list(i1)
                                  help_index=index_list(i1)
                                  c_list(i1)=c_list(i2)
                                  index_list(i1)=index_list(i2)
                                  c_list(i2)=help
                                  index_list(i2)=help_index
                               endif
                            enddo
                         enddo! end of sort of c_list
                      endif
                      exit fill
                   endif
                enddo fill
             endif
             if ( minval(c_list)/=-huge(1.0_r8_kind)) then
                refill :do i_list=1,N_LIST
                   if (c>c_list(i_list)) then
                      help_list = c_list ! help_list is needed for technical reasons
                      help_i_list = index_list
                      do i1=i_list+1,N_LIST ! shift all entries back by one index
                         c_list(i1)=help_list(i1-1) ! thereby keeping the N_LIST largest elements
                         index_list(i1) = help_i_list(i1-1)
                      enddo
                      c_list(i_list)=c
                      index_list(i_list)=j
                      exit refill
                   endif
                enddo refill
             endif
             if (c>c_max) then
                c_max=c
                index_largest=j
             endif
          enddo!                            ! Loop over MO`s in Irrep i
          write(output_unit,*)"EIGEN_LOCALIZE_HOLE: Largest dot-products for hole identification:"
          write(output_unit,'(10(2x,6ES13.3))')(c_list(i_list),i_list=1,N_LIST)
          do i_list=1,N_LIST-1
             diff=c_list(i_list)-c_list(i_list+1)
             if ((diff < hole_tolerance) .and. (min(c_list(i_list),-c_list(i_list+1))>hole_tolerance) ) then
                write(output_unit,*)" Warning: impossible to decide if hole is orbital",index_list(i_list),&
                     " or ",index_list(i_list+1)
             endif
          enddo

          write(output_unit,*)" Hole No. ",i_hole," in Irrep ",i," was identified as Orbital No. ",index_largest
          if ( index_largest/=holes_per_irrep(i)%m(i_hole) ) then
             call write_to_output_units&
                  ("eigen_localize_hole: Warning - hole was identified with a different orbital then before")
             ! that is why the whole story -- walking orbitals ..., that is normal! See below:
             write(output_unit,*)" "
             write(output_unit,*)" eigen_localize_hole: "
             if(.not.options_spin_orbit)then
             write(output_unit,*)" Index of Hole specified for Irrep ",symmetry_data_irrepname(i)
             else
             write(output_unit,*)" Index of Hole specified for Irrep ",symmetry_data_irrepname_proj(i)
             endif
             write(output_unit,*)" ",holes_per_irrep(i)%m(i_hole),"   Spin: ",spin_of_hole(i)%m(i_hole)
             write(output_unit,*)" Index of Hole found for Irrep     "
             write(output_unit,*)" ",index_largest,"   Spin: ",is,", will be used"
!!$             ! WHAT IS THIS ALL FOR?, WHAT DOES UPDATE_HOLE THEN MEAN???:
!!$             write(output_unit,*)" changing indices and eigenvector to those corresponding to current hole"
!!$             write(output_unit,*)" "
!!$             eigvec_hole(i)%m(:,i_hole,is) = &
!!$                  eigvec(i)%m(:,index_largest,is)
!!$             !             eigvec_hole(i)%m(:,holes_per_irrep(i)%m(i_hole),is) = &  ! bug
!!$             !                  eigvec(i)%m(:,index_largest,is)! bug ?
             ! the corresponding eigenvector (read at start) we dont want to change:
             ! will be used by occupation_module:
             holes_per_irrep(i)%m(i_hole)=index_largest
!!$             spin_of_hole(i)%m(i_hole) = is  ! actually this should not change anyway.
          else
             if (update_projector) then
                write(output_unit,*)" eigen_localize_hole: updating projector of hole"
                write(output_unit,*)"                      with indices   Irrep ",i,"   Index ",index_largest,"  Spin ",is
                write(output_unit,*)" "
                eigvec_hole(i)%m(:,i_hole,is) = eigvec(i)%m(:,index_largest,is)
             endif
          endif

       enddo holes

       if(allocated(help_hole))then
          deallocate(help_hole,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif
       if(allocated(help_hole_real))then
          deallocate(help_hole_real,help_hole_imag,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif
    enddo irreps

    deallocate(dims,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
  contains

    function scalar_product(irr,j,is) result(p)
      implicit none
      integer(i4_kind),intent(in) :: irr,j,is
      real(r8_kind)               :: p ! result
      ! *** end of interface ***

      real(r8_kind) :: re,im

!!$      DPRINT 'eid::scalar_product: irr,j,is=',irr,j,is

      if(options_spin_orbit)then
         ASSERT(is==1)
         re =   + sum(eigvec_real(irr)%m(:,j)*help_hole_real(:))&
              & + sum(eigvec_imag(irr)%m(:,j)*help_hole_imag(:))
         im =   + sum(eigvec_real(irr)%m(:,j)*help_hole_imag(:))&
              & - sum(eigvec_imag(irr)%m(:,j)*help_hole_real(:))
         p = sqrt(re**2+im**2)
      else
         p = sum(eigvec(irr)%m(:,j,is)*help_hole(:))
         p = abs(p)
      endif
    end function scalar_product

    subroutine pre_compute(irr,ih,is)
      implicit none
      integer(i4_kind),intent(in) :: irr,ih,is
      ! *** end of interface ***

      integer(i4_kind) :: j

!!$      DPRINT 'eid::pre_compute: irr,ih,is=',irr,ih,is

      if(options_spin_orbit)then
         ASSERT(is==1)
         do j=1,symmetry_data_dimension_proj(irr)
            help_hole_real(j) =&
                 & + sum(overlap_real(irr)%m(j,:)*eigvec_hole_real(irr)%m(:,ih))&
                 & - sum(overlap_imag(irr)%m(j,:)*eigvec_hole_imag(irr)%m(:,ih))
            help_hole_imag(j) =&
                 & + sum(overlap_real(irr)%m(j,:)*eigvec_hole_imag(irr)%m(:,ih))&
                 & + sum(overlap_imag(irr)%m(j,:)*eigvec_hole_real(irr)%m(:,ih))
         enddo
!!$         DPRINT 'eid::pre_compute: sum(help_hole_real)=',sum(help_hole_real)
!!$         DPRINT 'eid::pre_compute: sum(help_hole_imag)=',sum(help_hole_imag)
      else
         do j=1,symmetry_data_dimension(irr)
            help_hole(j) = sum(eigvec_hole(irr)%m(:,ih,is)*overlap(irr)%m(j,:))
         enddo
      endif
    end subroutine pre_compute
  end subroutine eigen_localize_hole


  subroutine eigen_data_bcast()
    ! call it from parallel context
    use comm, only: comm_parallel, comm_rank, comm_bcast
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: irr

    if( .not. comm_parallel() ) RETURN

    ! allocation on slaves, reallocate them if they allready exist,
    ! as the size may have changed
    if (comm_rank() /= 0) call realloc_eigvec_all()

    ! test if eigvec exists on master, then send to slaves
    ASSERT(allocated(eigvec))
    do irr=1,size(eigvec)
       ASSERT(allocated(eigvec(irr)%m))
       call comm_bcast(eigvec(irr)%m)
    enddo

    ! also send/receive eigenvalues, why not?
    ASSERT(allocated(eigval))
    do irr = 1, size(eigval)
      ASSERT(allocated(eigval(irr)%m))
      call comm_bcast(eigval(irr)%m)
    enddo
  end subroutine eigen_data_bcast

  !*************************************************************
  subroutine realloc_eigvec_all()
    ! Purpose: ensures that eigvec exists (designed for running on the
    !          slaves) and that they have the actual number, thus
    !          allocate them no matter if they already exists or not
    !------------ Modules used ------------------- ---------------
    use symmetry_data_module
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i,alloc_stat

    external error_handler
    !------------ Executable code --------------------------------

    ! if variables were already allocated then reallocate them to
    ! ensure proper size
    if(allocated(eigvec)) then
       do i=1,ssym%n_irrep
          deallocate(eigvec(i)%m,stat=alloc_stat)
          ASSERT (alloc_stat==0)
       end do
       deallocate(eigvec,stat=alloc_stat)
       ASSERT (alloc_stat==0)
    end if

    allocate(eigvec(ssym%n_irrep),stat=alloc_stat)
    ASSERT (alloc_stat==0)
    do i=1,ssym%n_irrep
       ! first allocate the memory
       allocate(eigvec(i)%m(ssym%dim(i),ssym%dim(i),ssym%n_spin), stat=alloc_stat)
       ASSERT (alloc_stat==0)
       ! then initialize it to zero
       eigvec(i)%m=0.0_r8_kind
    enddo
  end subroutine realloc_eigvec_all

#ifdef WITH_MATRIX_PARALLEL
  subroutine geigs_parallel(ham_tot, overlap, eigval, eigvec)
    use comm, only: comm_bcast
    use matrix_parallel, only: rmatrix, rdmatrix, matrix, array, geigs
    implicit none
    type(arrmat3), intent(inout) :: ham_tot(:) ! (n_irr), FIXME: intent, see bcast
    type(arrmat2), intent(in) :: overlap(:) ! (n_irr)
    type(arrmat2), intent(inout) :: eigval(:) ! (n_irr), needs to be allocated
    type(arrmat3), intent(inout) :: eigvec(:) ! (n_irr), needs to be allocated
    ! *** end of interface ***

    integer :: irr, spin
    type(rmatrix) :: H, S, V
    type(rdmatrix) :: e

    do irr = 1, size(overlap) ! n_irr

       S = matrix(overlap(irr)%m)

       !
       ! FIXME: Fock matrix is only valid on master so far, so
       ! broadcast:
       !
       call comm_bcast(ham_tot(irr)%m)

       do spin = 1, size(ham_tot(irr)%m, 3) ! n_spin

          H = matrix(ham_tot(irr)%m(:, :, spin))

          call geigs(H, S, e, V)

          eigval(irr)%m(:, spin) = array(e)
          eigvec(irr)%m(:, :, spin) = array(V)
       enddo
    enddo
  end subroutine geigs_parallel
#endif

      SUBROUTINE rsg_lsft(NM, N, A, B, W, MATZ, Z, ham_lsftM, IERR)
      use type_module
      implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) N,NM,IERR,MATZ!,ind
      REAL(kind=r8_kind) A(:,:),B(:,:),W(:),Z(:,:)
      REAL(kind=r8_kind) ham_lsftM(:, :)
      REAL(kind=r8_kind),allocatable :: W_loc(:),Z_loc(:,:),FV1(:),FV2(:)
!     REAL(kind=r8_kind),allocatable :: temp(:,:),temp1(:,:)

!
!    THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
!    SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
!    TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
!    FOR THE REAL SYMMETRI!GENERALIZED EIGENPROBLEM  AX = (LAMBDA)BX.
!
!    ON INPUT
!
!       NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
!       ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!       DIMENSION STATEMENT.
!
!       N  IS THE ORDER OF THE MATRICES  A  AND  B.
!
!       A  CONTAINS A REAL SYMMETRI!MATRIX.
!
!       B  CONTAINS A POSITIVE DEFINITE REAL SYMMETRIC MATRIX.
!
!       MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
!       ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
!       ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
!
!    ON OUTPUT
!
!       W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
!
!       Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
!
!       IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
!          COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
!          AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
!
!    QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!    MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!   XS THIS VERSION DATED AUGUST 1983.
!
!   adapted to f90 and parameter list modified TB 5/97
!   requires interface block in calling code !!!
!
!     ------------------------------------------------------------------

      WARN('RSG modifies arguments!')
      allocate (W_loc(N),Z_loc(NM,N),FV1(N),FV2(N))
!
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10   continue

        CALL  REDUC(NM,N,A,B,FV2,IERR)
        if(make_level_shift) &
                     A=A+ham_lsftM
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W_loc,FV1,FV2)
      CALL  TQLRAT(N,W_loc,FV2,IERR)
      GO TO 50

   20 continue
      CALL  TRED2(NM,N,A,W_loc,FV1,Z_loc)
      CALL  TQL2(NM,N,W_loc,FV1,Z_loc,IERR)
            if(lvshift_mixing) &
                     ham_lsftM=Z_loc
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAK(NM,N,B,FV2,N,Z_loc)

   50 continue
      Z = Z_loc(1:NM,1:N)
      W = W_loc(1:N)
      deallocate(W_loc,Z_loc,FV1,FV2)
      RETURN

      END subroutine rsg_lsft


!--------------- End of module ----------------------------------
end module eigen_data_module

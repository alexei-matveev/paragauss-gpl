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
module matrix_parallel
  !-------------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !  Names are  mangled by suffixes  indicating the types  of accepted
  !  arguments:
  !
  !     0       double precision scalar
  !     1       double precision array
  !             interpreted as diagonal matrix
  !     2       double precision array with two axes
  !             interpreted as a rectangular matrix
  !     r       type(rmatrix)
  !     c       type(cmatrix)
  !     d       type(rdmatrix)
  !     h       type(chmatrix)
  !     i       integer
  !
  !  Return value, when indicated, follows the double underscore.
  !
  ! Copyright (c) 2011-2012 Martin Roderus
  ! Copyright (c) 2011-2012 Alexei Matveev
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  !-------------------------------------------------------------------
# include "def.h"
  use type_module, only: RK => r8_kind  ! type specification parameters
  USE_MPI, only: MPI_COMM_NULL
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------

  type, public :: pm_ctxt
     private
     integer :: blacs_ctxt = -1
     integer :: mpi_comm = MPI_COMM_NULL
  end type pm_ctxt

  type, public :: rmatrix
     private
     integer :: mpi_comm = MPI_COMM_NULL
     integer :: desc(9) = -1
     real(RK), allocatable :: m(:,:)
  end type rmatrix

  type, public :: rdmatrix
     private
     real(RK), allocatable :: d(:)
  end type rdmatrix

  interface matrix
     !
     ! Convert array into dense distributed or diagonal matrix:
     !
     module procedure matrix_2
     module procedure matrix_1
  end interface

  public :: matrix

  interface array
     !
     ! Convert an eventually distributed object to a plain array:
     !
     module procedure array_r
     module procedure array_d
  end interface

  public :: array

  interface operator(*)
     module procedure mult_r_r  ! RK matrix * matrix

     module procedure mult_d_d  ! RK diagonal * diagonal

     module procedure mult_r_d  ! RK matrix * diagonal
     module procedure mult_d_r  ! RK diagonal * matrix

     module procedure mult_r_0  ! RK matrix * scalar
     module procedure mult_0_r  ! RK scalar * matrix

     module procedure mult_d_0  ! RK diagonal * scalar
     module procedure mult_0_d  ! RK scalar * diagonal
  end interface operator(*)

  public :: operator(*)

  interface mult
     module procedure mult_d_r_d
  end interface

  public :: mult

  interface operator(+)
     module procedure plus_r_r  ! RK matrix + matrix

     module procedure plus_d_d  ! RK diagonal + diagonal

     module procedure plus_r_d  ! RK matrix + diagonal
!    module procedure plus_d_r  ! RK diagonal + matrix, FIXME: missing
  end interface operator(+)

  public :: operator(+)

  interface operator(-)
     module procedure minus_r_r
     module procedure minus_r
     module procedure minus_d_d
     module procedure minus_d
  end interface

  public :: operator(-)

  interface operator(**)
     module procedure pow_d_i
  end interface

  public :: operator(**)

  interface size
     module procedure size_d!(rdmatrix)
     module procedure size_r_i!(rmatrix, axis)
  end interface

  public :: size

  public :: geigs               ! (H, S, e, V)
  public :: rpt                 ! (e, V) -> rmatrix
  public :: tr                  ! rmatrix -> rmatrix
  public :: sim                 ! (U, V) -> rmatrix
  public :: maxabs              ! rmatrix -> double
  public :: show_m              ! rmatrix -> IO

  public :: pm_create_ctxt
  public :: pm_get_ctxt

  !===================================================================
  ! End of public interface of module
  !===================================================================

  integer, parameter :: BLOCKSIZE = 64
  integer, parameter :: ctxt_=2, m_=3, n_=4

  !
  ! Deafult context, set automatically to some meanigfull value on the
  ! first use of matrix_2() without context:
  !
  type(pm_ctxt) :: default_context = pm_ctxt(-1, MPI_COMM_NULL)

contains

  function pm_create_ctxt(mpi_comm) result(pmctxt)
    !
    ! Generate  BLACS  context from  MPI  communicator.  Note that  in
    ! general not all MPI process will be part of the BLACS grid.
    !
    USE_MPI, only: MPI_SUCCESS
    use f77_scalapack, only: blacs_gridmap
    implicit none
    integer, intent(in) :: mpi_comm
    type(pm_ctxt) :: pmctxt
    ! *** end of interface ***

    integer :: comm_size, ierr ! MPI
    integer :: nprow, npcol ! BLACS
    integer, allocatable :: proc_grid(:,:)
    integer :: i, j

    pmctxt%mpi_comm = mpi_comm

    call MPI_Comm_size(mpi_comm, comm_size, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    call blacs_proc_grid(comm_size, nprow, npcol)
    allocate(proc_grid(0:nprow-1, 0:npcol-1))

    !
    ! rank_of ()  encodes either column-major,  row-major or otherwise
    ! grid:
    !
    forall (i = 0:nprow-1, j = 0:npcol-1)
       proc_grid(i, j) = rank_of (nprow, npcol, i, j)
    end forall

    !
    ! First argument of blacs_gridmap is inout: reads the system context
    ! (MPI communicator) and returns a BLACS context
    !
    pmctxt%blacs_ctxt = mpi_comm
    call blacs_gridmap (pmctxt%blacs_ctxt, proc_grid, size(proc_grid, 1), nprow, npcol)
  end function pm_create_ctxt


  pure function rank_of (nprow, npcol, row, col) result (rank)
    !
    ! Encodes the grid map.
    !
    ! Use  row-major order  for  cartesian grids,  which  is also  the
    ! default for MPI. Original code used column-major order.
    !
    ! MPI_TYPE_CREATE_DARRAY assumes row-major cartesian grid, if this
    ! does  not hold,  the custom  data  type used  for one-bcast  per
    ! worker implementation of collective array_r () is wrong.
    !
    implicit none
    integer, intent(in) :: nprow, npcol
    integer, intent(in) :: row, col ! base-0
    integer :: rank
    ! *** end of interface ***

    !
    ! Rank of a (row, col) processor in the BLACS grid:
    !
    rank = row * npcol + col ! row-major
    !ank = row + col * nprow ! column-major
  end function rank_of


  pure subroutine blacs_proc_grid(nproc, nprow, npcol)
    !
    ! Purpose: Input:  total number of processes p.  Output: number of
    ! row- and  column-processes of the  BLACS process grid.   If p<9,
    ! nprow=1 and npcol=p. Else nprow=npcol=floor(sqrt(p))
    !
    implicit none
    integer, intent(in)  :: nproc
    integer, intent(out) :: nprow, npcol
    ! *** end of interface ***

    if (nproc < 9) then
       ! Less  than  9  processors   alloted  to  the  task,  generate
       ! 1-dimensional process grid:
       nprow = 1
       npcol = nproc
    else
       ! 9 processors or more alloted, generate square process grid:
       nprow = sqrt(real(nproc))
       npcol = nprow
    end if
  end subroutine blacs_proc_grid


  function pm_in_procpool(pmctxt) result(answer)
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(pm_ctxt), intent(in) :: pmctxt
    logical :: answer
    ! *** end of interface ***

    integer :: nprow, npcol, myrow, mycol ! BLACS

    call blacs_gridinfo (pmctxt%blacs_ctxt, nprow, npcol, myrow, mycol)

    answer = (myrow >= 0)
  end function pm_in_procpool


  function matrix_2(array, context) result(matrix)
    USE_MPI, only: MPI_COMM_NULL
    use comm, only: comm_world
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    real(RK), intent(in) :: array(:,:)
    type(pm_ctxt), intent(in), optional :: context
    type(rmatrix) :: matrix
    ! *** end of interface ***

    type(pm_ctxt) :: ctxt
    integer :: nprow, npcol, myrow, mycol

    if (present(context)) then
       ctxt = context
    else
       !
       ! FIXME:  this is  an ugly  hack to  allow use  of constructors
       ! without explicitly managing contexts:
       !
       if (default_context%mpi_comm == MPI_COMM_NULL) then
          default_context = pm_create_ctxt(comm_world)
       endif
       ctxt = default_context
    endif

    call allocate_rmatrix(ctxt, size(array, 1), size(array, 2), matrix)

    !
    ! Now each  worker takes its  section from the  array (if it  is a
    ! part of the BLACS grid, of course:
    !
    call blacs_gridinfo (ctxt % blacs_ctxt, nprow, npcol, myrow, mycol)

    if (myrow < 0 .or. mycol < 0) return

    call extract (array, matrix % m, nprow, npcol, myrow, mycol)
  end function matrix_2


  pure function matrix_1(array) result(matrix)
    implicit none
    real(RK), intent(in) :: array(:)
    type(rdmatrix) :: matrix
    ! *** end of interface ***

    matrix = rdmatrix(array)
  end function matrix_1


  function block_is_in_local_pmat(pmctxt, i, j) result(answer)
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(pm_ctxt), intent(in)  :: pmctxt
    integer, intent(in) :: i, j
    ! *** end of interface ***

    logical :: answer
    integer :: nprow, npcol, myrow, mycol


    call Blacs_gridinfo (pmctxt%blacs_ctxt, nprow, npcol, myrow, mycol)

    answer = myrow >= 0  .and.  mod(i-1, nprow) == myrow  .and.  mod(j-1, npcol) == mycol
  end function block_is_in_local_pmat


  function pm_get_ctxt(pmat) result(ctxt)
    implicit none
    type(rmatrix), intent(in) :: pmat
    type(pm_ctxt) :: ctxt
    ! *** end of interface ***

    ! type constructor here:
    ctxt = pm_ctxt(pmat%desc(ctxt_), pmat%mpi_comm)
  end function pm_get_ctxt

#if 0
  subroutine cb_collect_array(pmctxt, mat, pmat, i, j, m, n)
    !
    ! Collective routine to collect a parallel matrix (pmat) and store
    ! it  in a  local array.   After the  routine call,  the  array is
    ! duplicated on each process.
    !
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    type(pm_ctxt), intent(in)  :: pmctxt
    real(RK), intent(out)  :: mat(:,:)
    real(RK), intent(in) :: pmat(:,:)
    integer, intent(in) :: i, j, m, n
    ! *** end of interface ***

    integer :: nprow, npcol
    integer :: rootrank, rootrow, rootcol, mpi_size, ierr
    integer :: i_start_mat, i_end_mat, j_start_mat, j_end_mat
    integer :: i_start_pmat, i_end_pmat, j_start_pmat, j_end_pmat
    integer :: nr_of_block_elements

    !integer :: myrank !MR_DBG

    call MPI_Comm_size(pmctxt%mpi_comm, mpi_size, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    call blacs_proc_grid(mpi_size, nprow, npcol)


    ! Set start- and end-indices of the block to be processed
    i_start_mat = get_coordinate_mat(i)
    i_end_mat   = i_start_mat + get_tmp_bs(i, m) - 1
    j_start_mat = get_coordinate_mat(j)
    j_end_mat   = j_start_mat + get_tmp_bs(j, n) - 1

    i_start_pmat = get_coordinate_pmat(i, nprow)
    i_end_pmat   = i_start_pmat + get_tmp_bs(i, m) - 1
    j_start_pmat = get_coordinate_pmat(j, npcol)
    j_end_pmat   = j_start_pmat + get_tmp_bs(j, n) - 1


    if(block_is_in_local_pmat(pmctxt, i, j)) then
       mat( i_start_mat : i_end_mat , j_start_mat : j_end_mat ) = &
            pmat( i_start_pmat : i_end_pmat , j_start_pmat : j_end_pmat)
    end if

    !
    ! Here the column-major, row-major  or otherwise grid order of the
    ! BLACS grid is assumed by use of root_of ():
    !
    rootrow = mod(i - 1, nprow)
    rootcol = mod(j - 1, npcol)
    rootrank = rank_of (nprow, npcol, rootrow, rootcol)

    nr_of_block_elements = get_tmp_bs(i, m)*get_tmp_bs(j, n)

    call MPI_Bcast( mat( i_start_mat : i_end_mat , j_start_mat : j_end_mat ), &
            & nr_of_block_elements, MPI_DOUBLE_PRECISION, rootrank,&
            & pmctxt%mpi_comm, ierr)
    ASSERT(ierr==MPI_SUCCESS)

  contains
    pure function get_tmp_bs(index, mat_size) result(tmp_bs)
      implicit none
      integer, intent(in) :: index, mat_size
      integer :: tmp_bs
      ! *** end of interface ***

      tmp_bs = min(BLOCKSIZE, mat_size - get_coordinate_mat(index) + 1)
    end function get_tmp_bs

    pure function get_coordinate_mat(index) result(coordinate)
      implicit none
      integer, intent(in) :: index
      integer :: coordinate
      ! *** end of interface ***

      coordinate = (index - 1) * BLOCKSIZE + 1
    end function get_coordinate_mat

    pure function get_coordinate_pmat(index, nprowcol) result(coordinate)
      implicit none
      integer, intent(in) :: index
      integer, intent(in) :: nprowcol
      integer :: coordinate
      ! *** end of interface ***

      coordinate = (index - 1) / nprowcol * blocksize + 1
    end function get_coordinate_pmat
  end subroutine cb_collect_array


  function array_r(matrix) result(array)
    implicit none
    type(rmatrix), intent(in) :: matrix
    real(RK) :: array(size(matrix, 1), size(matrix, 2)) ! FIXME: use of pure func here
    ! *** end of itnerface ***

    integer :: i, j
    integer :: m, n, m_blocks, n_blocks
    type(pm_ctxt) :: ctxt

    ctxt = pm_get_ctxt(matrix)

    m = size(matrix, 1)
    n = size(matrix, 2)
    m_blocks = (m + BLOCKSIZE - 1) / BLOCKSIZE
    n_blocks = (n + BLOCKSIZE - 1) / BLOCKSIZE

    do j = 1, n_blocks
       do i = 1, m_blocks
          call cb_collect_array(ctxt, array, matrix%m, i, j, m, n)
       enddo
    enddo
  end function array_r
#else
  function array_r (matrix) result (array)
    !
    ! Convert a  distributed matrix to a  replicated array. Collective
    ! operation.
    !
    ! FIXME: The block cyclic distribution  of the array over the grid
    ! of processes  assumes the row-major  topology of the  BLACS grid
    ! (this is  hardwired into MPI_TYPE_CREATE_DARRAY).  So this order
    ! should be consistent with the one implemented by rank_of ().
    !
    USE_MPI, only: MPI_SUCCESS, MPI_DOUBLE_PRECISION
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(rmatrix), intent(in) :: matrix
    real(RK) :: array(size(matrix, 1), size(matrix, 2)) ! user-defined pure func here
    ! *** end of interface ***

    integer, parameter :: ndims = 2
    integer :: psizes(ndims)
    integer :: rank, row, col, ierr
    integer :: myrow, mycol, myrank
    integer :: dtype
    type(pm_ctxt) :: ctxt

    ctxt = pm_get_ctxt (matrix)

    call MPI_COMM_RANK (ctxt % mpi_comm, myrank, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    !
    ! Fill  grid   shape  into  psizes,  should   be  consistent  with
    ! conventions in pm_create_ctxt():
    !
    call blacs_gridinfo (ctxt % blacs_ctxt, psizes(1), psizes(2), myrow, mycol)

    ! FIXME: is this info valid on those that are not part of the grid?
    ASSERT(all(psizes>0))

    !
    ! Each member of the process grid holds a part of data, bcast that
    ! to everyone  in the MPI  world (not only  to those in  the BLACS
    ! grid):
    !
    rank = -1
    do row = 0, psizes (1) - 1
       do col = 0, psizes (2) - 1
          rank = rank + 1

          !
          ! MPI_TYPE_CREATE_DARRAY  assumes row-major  cartesian grid,
          ! if this does not hold the custom dtype is wrong:
          !
          if (rank /= rank_of (psizes (1), psizes (2), row, col)) then
             ABORT('need row-major')
          endif

          !
          ! This type should represent the "section" of array that has to
          ! be received from rank:
          !
          dtype = make_type (product (psizes), rank, shape (array), psizes)

          if (rank == myrank) then
             ASSERT(row==myrow)
             ASSERT(col==mycol)
             !
             ! FIXME: check if this works reliably with your MPI. This
             ! did not  work quite reliably with a  6-digit tag. Since
             ! then the tag was changed to a smaller one, but no tests
             ! have been done.
             !
             ! call copy_type (matrix % m, array, dtype)

             !
             ! This is a manual copy:
             !
             call insert (matrix % m, array, psizes(1), psizes(2), row, col)
          endif

          !
          ! Now that  the data on  rank was already copied  into array
          ! section,  bcast that data  to other  workers in  MPI world
          ! (not only to those in the BLACS grid):
          !
          call MPI_BCAST (array, 1, dtype, rank, ctxt % mpi_comm, ierr)
          ASSERT(ierr==MPI_SUCCESS)

          call free_type (dtype)
       enddo
    enddo

  contains

    function make_type (size, rank, gsizes, psizes) result (newtype)
      USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_DISTRIBUTE_CYCLIC, &
           MPI_ORDER_FORTRAN, MPI_SUCCESS
      implicit none
      integer, intent(in) :: size, rank
      integer, intent(in) :: gsizes(:) !(ndims)
      integer, intent(in) :: psizes(:) !(ndims)
      integer :: newtype
      ! *** end of interface ***

      integer, parameter :: distribs(ndims) = [MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC]
      integer, parameter :: dargs(ndims) = [BLOCKSIZE, BLOCKSIZE]
      integer :: ierr

      call MPI_TYPE_CREATE_DARRAY (size, rank, ndims, gsizes, distribs, dargs, psizes, &
           MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  newtype, ierr)
      ASSERT(ierr==MPI_SUCCESS)

      call MPI_TYPE_COMMIT (newtype, ierr)
      ASSERT(ierr==MPI_SUCCESS)
    end function make_type

    subroutine free_type (newtype)
      USE_MPI, only: MPI_SUCCESS
      implicit none
      integer, intent(in) :: newtype
      ! *** end of interface ***

      integer :: ierr

      call MPI_TYPE_FREE (newtype, ierr)
      ASSERT(ierr==MPI_SUCCESS)
    end subroutine free_type

    subroutine copy_type (section, array, dtype)
      !
      ! Copy the data into proper section of the array by sending to
      ! myself.
      !
      USE_MPI, only: MPI_SUCCESS, MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE, MPI_COMM_SELF
      implicit none
      real(RK), intent(in) :: section (:, :)
      real(RK), intent(inout) :: array (:, :)
      integer, intent(in) :: dtype
      ! *** end of interface ***

      integer, parameter :: tag = 21042 !reduce below the guaranteed value of 32767
                                        !+ 012 just a date
      integer :: req, status(MPI_STATUS_SIZE)

      !
      ! FIXME:  both  variants  fail  occasionally on  Debian  Lenny
      ! (OpenMPI 1.2.7 RC2) for unknown reason:
      !
#if 1
      call MPI_IRECV (array, 1, dtype, 0, tag, MPI_COMM_SELF, req, ierr)
      ASSERT(ierr==MPI_SUCCESS)

      ! sending a continous array:
      call MPI_SEND (section, size (section), MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_SELF, ierr)
      ASSERT(ierr==MPI_SUCCESS)
#else
      ! sending a continous array:
      call MPI_ISEND (section, size(section), MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_SELF, req, ierr)
      ASSERT(ierr==MPI_SUCCESS)

      call MPI_RECV (array, 1, dtype, 0, tag, MPI_COMM_SELF, ierr)
      ASSERT(ierr==MPI_SUCCESS)
#endif

      call MPI_WAIT (req, status, ierr)
      ASSERT(ierr==MPI_SUCCESS)
    end subroutine copy_type
  end function array_r
#endif

  pure function array_d(a) result(array)
    !
    ! We  represent  diagonal matrices  by  arrays  replicated on  all
    ! workers in the MPI world.   NOTE: no operation should produce an
    ! rdmatrix that  violates this  constrain. See eigensolver  for an
    ! example of how to deal with that.
    !
    USE_MPI, only: MPI_DOUBLE_PRECISION
    implicit none
    type(rdmatrix), intent(in) :: a
    real(RK) :: array(size(a))
    ! *** end of interface ***

    !
    ! Due to  implicit constrain (to  be) fullfilled by  all functions
    ! that return a  valid rdmatrix, we assume here  the data is valid
    ! on all workers in MPI world:
    !
    array = a%d
  end function array_d


  function pm_index_mat_to_index_pmat(pmctxt, i_mat, j_mat, i_pmat, j_pmat) result(is_in_pmat)
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(pm_ctxt), intent(in) :: pmctxt
    integer, intent(in)       :: i_mat, j_mat
    integer, intent(out)      :: i_pmat, j_pmat
    logical                   :: is_in_pmat

    integer :: i_off, j_off, i_block, j_block
    integer :: nprow, npcol, myrow, mycol


    i_block = index_mat_to_blockindex(i_mat)
    j_block = index_mat_to_blockindex(j_mat)

    is_in_pmat = block_is_in_local_pmat(pmctxt, i_block, j_block)

    if(is_in_pmat) then
       call Blacs_gridinfo (pmctxt%blacs_ctxt, nprow, npcol, myrow, mycol)
       i_off = (i_block-1)/nprow*blocksize
       i_pmat = i_off + i_mat-(i_block-1)*blocksize

       j_off = (j_block-1)/npcol*blocksize
       j_pmat = j_off + j_mat-(j_block-1)*blocksize
    else
       i_pmat = -1
       j_pmat = -1
    end if

  end function pm_index_mat_to_index_pmat


  function index_mat_to_blockindex(i_mat) result(i_block)
    implicit none
    integer, intent(in) :: i_mat
    integer :: i_block

    i_block = (i_mat-1)/blocksize+1

  end function index_mat_to_blockindex


  subroutine allocate_rmatrix(pmctxt, m, n, pmat)
    !
    ! On  the  workers  not in  the  BLACS  grid  %m is  allocated  to
    ! zero-sized array and %desc is only partially initialized.
    !
    use f77_scalapack, only: descinit, blacs_gridinfo, numroc
    implicit none
    type(pm_ctxt), intent(in) :: pmctxt
    integer, intent(in) :: m, n ! matrix dimensions
    type(rmatrix), intent(inout) :: pmat ! FIXME: intent(out)?
    ! *** end of interface ***

    ! local variables
    integer :: nprow, npcol, myrow, mycol ! BLACS
    integer :: loc_rows, loc_cols, info

    pmat%mpi_comm = pmctxt%mpi_comm

    call blacs_gridinfo (pmctxt%blacs_ctxt, nprow, npcol, myrow, mycol)

    loc_rows = max (numroc (m, blocksize, myrow, 0, nprow), 0)
    loc_cols = max (numroc (n, blocksize, mycol, 0, npcol), 0)

    !loc_rows = max(loc_rows, 0)
    !loc_cols = max(loc_cols, 0)

    !print *, myrow, mycol, ': allocate_rmatrix loc_rows: ', loc_rows !MR_DBG
    !print *, myrow, mycol, ': allocate_rmatrix loc_cols: ', loc_cols !MR_DBG

    ! FIXME: check for allocation status of p* data structure
    allocate(pmat%m(loc_rows, loc_cols))

    if(pm_in_procpool(pmctxt)) then
       call descinit(pmat%desc, m, n, blocksize, blocksize, 0, 0, pmctxt%blacs_ctxt, max(1, loc_rows), info)
    else
       pmat%desc(ctxt_) = pmctxt%blacs_ctxt
       pmat%desc(m_) = m
       pmat%desc(n_) = n
    end if
  end subroutine allocate_rmatrix


  function pm_get_bs() result(bs)
    implicit none
    integer :: bs

    bs = blocksize

  end function pm_get_bs


  subroutine show_m(a)
    implicit none
    type(rmatrix), intent(in) :: a
    ! *** end of interface ***

    real(RK) :: a_(size(a, 1), size(a, 2))
    integer :: j

    a_ = array(a)
    do j = 1, size(a_, 2)
       write(*, "(20F20.10)") a_(:, j)
    enddo
  end subroutine show_m


  function maxabs(a) result(b)
    !
    ! Same as maxval(abs(array(a)))
    !
    USE_MPI, only: MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX
    implicit none
    type(rmatrix), intent(in) :: a
    real(RK) :: b
    ! *** end of interface ***

    integer :: ierr

    b = maxval(abs(a%m))

    call MPI_ALLREDUCE(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, a%mpi_comm, ierr)
    if ( ierr /= 0 ) stop "maxabs: error in MPI_ALLREDUCE"
    ! print *, "diff=", b - maxval(abs(array(a)))
  end function maxabs


  subroutine geigs(A, B, e, V)
    !
    ! Returns  a freshly allocated  rdmatrix e  and rmatrix  V holding
    ! eigenvalues and eigenvectors, respectively.
    !
    USE_MPI, only: MPI_DOUBLE_PRECISION
    use f77_scalapack, only: pdsygvx, blacs_gridinfo
    implicit none
    type(rmatrix), intent(in) :: A, B
    type(rdmatrix), intent(out) :: e
    type(rmatrix), intent(out) :: V
    ! *** end of interface ***

    ! BLACS variables
    integer :: myrow, mycol
    integer :: nprow, npcol

    ! PDSYGVX temporary storage:
    integer :: lwork = -1, liwork = -1
    integer, allocatable :: iwork(:), ifail(:), iclustr(:)
    integer :: iworkdummy(10)
    double precision, allocatable :: work(:)
    double precision              :: workdummy(10)

    ! PDSYGVX parameters
    integer,     parameter :: itype = 1
    character, parameter :: jobz = 'V'
    character, parameter :: range = 'A'
    character, parameter :: uplo = 'L'
    integer, parameter     :: ia=1, ja=1, ib=1, jb=1, iz=1
    integer, parameter     :: jz=1, il=1, iu=1
    double precision, parameter :: abstol = 0.0

    !  ORFAC  (global input)  Specifies which  eigenvectors  should be
    !          reorthogonalized.   Eigenvectors   that  correspond  to
    !          eigenvalues which are  within tol=ORFAC*norm(A) of each
    !          other  are  to be  reorthogonalized.   However, if  the
    !          workspace  is  insufficient  (see  LWORK), tol  may  be
    !          decreased until all eigenvectors to be reorthogonalized
    !          can be  stored in one  process.  No reorthogonalization
    !          will be done if ORFAC  equals zero.  A default value of
    !          10^-3 is  used if ORFAC  is negative.  ORFAC  should be
    !          identical on all processes.
    !
    ! FIXME: literal constants for eigenvalue spacing and "accidental"
    ! degeneracies here:
    !
    double precision, parameter :: orfac = 1.0d-12
    integer, parameter :: cluster_size = 7
    double precision, parameter :: vl = 0.0, vu = 0.0

    ! PDSYGV variables
    integer :: n
    integer :: m  ! Number of eigenvalues found
    integer :: nz  ! Number of eigenvectors found
    double precision, allocatable :: gap(:)
    integer :: info

    !
    ! Note  the  difference between  size(B%m,  1),  size(B%m, 2)  and
    ! size(B, 1) == size(B, 2):
    !
    real(RK) :: A_(size(A%m, 1), size(A%m, 2))
    real(RK) :: B_(size(B%m, 1), size(B%m, 2))
    real(RK) :: e_(size(A, 1)) ! matrices should be square

    TRACE("geigs_pdsygvx/entered")
    !
    ! Egiensolvers used  to "destroy" the input matrices,  make a copy
    ! of the actual data. Note that Gfortran 4.3 seems to erroneousely
    ! optimize away  a copy if you do  the for the whole  rmatrix as a
    ! struct:
    !
    A_ = A%m
    B_ = B%m

    call blacs_gridinfo (A%desc(ctxt_), nprow, npcol, myrow, mycol)

    n = size(A, 1)

    !ASSERT(size(A, 2) == n)
    !ASSERT(size(B, 1) == n)
    !ASSERT(size(B, 2) == n)

    !
    ! FIXME: comment on what this does or supposed to do on workers
    ! that are not involved:
    !
    call allocate_rmatrix(pm_get_ctxt(A), n, n, V)

    !
    ! Implicit constrains  on the representation  of diagonal matrices
    ! require us to return valid eigenvalues on all workers in the MPI
    ! world. So even if some workers do not participate in solving the
    ! eigenvalue problem, they still have to return a valid rdmatrix:
    !
    if (myrow < 0) goto 100 ! bcast eigenvalues, make rdmatrix, return

    ! First run for workspace querry
    lwork = -1
    liwork = -1

    allocate(iclustr(2*nprow*npcol))
    allocate(gap(nprow*npcol))
    allocate(ifail(n))

    call pdsygvx(itype, jobz, range, uplo, n, A_, ia, ja, A%desc, &
         & B_, ib, jb, B%desc, vl, vu, il, iu, abstol, m, nz, e_,&
         & orfac, V%m, iz, jz, V%desc, workdummy, lwork, iworkdummy,&
         & liwork, ifail, iclustr, gap, info)

    if (info /= 0) then
       print *, "geigs: ERROR: pdsygvx() returned info=", info, "in inquiry run"
       stop "error in pdsygvx, see tty"
    endif

    !
    ! Quote    from   http://www.netlib.org/scalapack/double/pdsyevx.f
    ! (reference  implementation). On  output "if  JOBZ='V'  WORK(1) =
    ! optimal amount of workspace  required to compute eigenvalues and
    ! eigenvectors efficiently with  *no guarantee on orthogonality*."
    ! Emphasis mine. So the value  returned in workdummy(1) may not be
    ! what one needs.
    !
    ! Another quote: "The computed  eigenvectors may not be orthogonal
    ! if the minimal workspace is supplied and ORFAC is too small.  If
    ! you want to guarantee  orthogonality (at the cost of potentially
    ! poor  performance)  you  should  add  the  following  to  LWORK:
    ! (CLUSTERSIZE-1)*N"
    !
    lwork = workdummy(1) + (cluster_size - 1) * n
    liwork = iworkdummy(1)

    ! Allocate workspace
    allocate(work(lwork))
    allocate(iwork(liwork))

    call pdsygvx(itype, jobz, range, uplo, n, A_, ia, ja, A%desc, &
         & B_, ib, jb, B%desc, vl, vu, il, iu, abstol, m, nz, e_,&
         & orfac, V%m, iz, jz, V%desc, work, lwork, iwork,&
         & liwork, ifail, iclustr, gap, info)

    if (info /= 0) then
       if (mod(info/2, 2) /= 0) then
          print *, "iclustr=", iclustr
          print *, "gap=", gap
       endif
       print *, "geigs: ERROR: pdsygvx() returned info=", info
       stop "error in pdsygvx, see tty"
    endif

100 CONTINUE
    !
    ! Process  0 is  always part  of  the process  grid see  allotment
    ! algorithm  in pm_create_ctxt.   FIXME: bcast  is only  needed if
    ! BLACS grid /= MPI world:
    !
    call MPI_Bcast(e_, size(e_), MPI_DOUBLE_PRECISION, 0, A%mpi_comm, info)
    if (info /= 0) stop "MPI_Bcast"

    e = matrix(e_)
    TRACE("geigs_pdsygvx/return")
  end subroutine geigs


  function plus_r_d(A, B) result(C)
    implicit none
    type(rmatrix), intent(in) :: A
    type(rdmatrix), intent(in) :: B
    type(rmatrix) :: C
    ! *** end of interface ***

    integer :: i
    integer :: i_pmat, j_pmat
    integer :: m_a, n_a, m_b

    m_a = size(A%m, 1)
    n_a = size(A%m, 2)
    m_b = size(B%d)

    !ASSERT(m_b == m_a  .or.  m_b == n_a)

    C=A

    do i=1, size(B%d)
       if(pm_index_mat_to_index_pmat(pm_get_ctxt(A), i, i, i_pmat, j_pmat)) then
          C%m(i_pmat, j_pmat) = C%m(i_pmat, j_pmat) + B%d(i)
       end if
    end do
  end function plus_r_d


  pure function plus_d_d(A, B) result(C)
    implicit none
    type(rdmatrix), intent(in) :: A, B
    type(rdmatrix) :: C
    ! *** end of interface ***

    C = matrix(array(A) + array(B))
  end function plus_d_d


  pure function minus_d_d(A, B) result(C)
    implicit none
    type(rdmatrix), intent(in) :: A, B
    type(rdmatrix) :: C
    ! *** end of interface ***

    C = matrix(array(A) - array(B))
  end function minus_d_d


  pure function minus_d(A) result(C)
    implicit none
    type(rdmatrix), intent(in) :: A
    type(rdmatrix) :: C
    ! *** end of interface ***

    C = matrix(-array(A))
  end function minus_d


  pure function mult_d_0(A, alpha) result(B)
    ! FIXME: not tested yet
    implicit none
    type(rdmatrix), intent(in) :: A
    real(RK), intent(in) :: alpha
    type(rdmatrix) :: B
    ! *** end of interface ***

    B = matrix(array(A) * alpha)
  end function mult_d_0


  pure function mult_0_d(alpha, A) result(B)
    ! FIXME: not tested yet
    implicit none
    real(RK), intent(in) :: alpha
    type(rdmatrix), intent(in) :: A
    type(rdmatrix) :: B
    ! *** end of interface ***

    B = matrix(alpha * array(A))
  end function mult_0_d


  pure function mult_d_d(A, B) result(C)
    implicit none
    type(rdmatrix), intent(in) :: A, B
    type(rdmatrix) :: C
    ! *** end of interface ***

    C = matrix(array(A) * array(B))
  end function mult_d_d


  function mult_r_0(A, alpha) result(B)
    implicit none
    type(rmatrix), intent(in) :: A
    real(RK), intent(in) :: alpha
    type(rmatrix) :: B
    ! *** end of interface ***

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(A, 2), B)

    B%m = A%m * alpha
  end function mult_r_0


  function mult_0_r(alpha, A) result(B)
    implicit none
    real(RK), intent(in) :: alpha
    type(rmatrix), intent(in) :: A
    type(rmatrix) :: B
    ! *** end of interface ***

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(A, 2), B)

    B%m = alpha * A%m
  end function mult_0_r


  function plus_r_r(A, B) result(C)
    implicit none
    type(rmatrix), intent(in) :: A, B
    type(rmatrix) :: C
    ! *** end of interface ***

    !ASSERT(size(A,1)==size(B,1))
    !ASSERT(size(A,2)==size(B,2))
    !ASSERT(A%mpi_comm) == B%mpi_comm)
    !ASSERT(A%desc(ctxt_)) == B%desc(ctxt_))

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(A, 2), C)

    C%m = A%m + B%m
  end function plus_r_r


  function minus_r_r(A, B) result(C)
    implicit none
    type(rmatrix), intent(in) :: A, B
    type(rmatrix) :: C
    ! *** end of interface ***

    !ASSERT(size(A,1)==size(B,1))
    !ASSERT(size(A,2)==size(B,2))
    !ASSERT(A%mpi_comm) == B%mpi_comm)
    !ASSERT(A%desc(ctxt_)) == B%desc(ctxt_))

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(A, 2), C)

    C%m = A%m - B%m
  end function minus_r_r


  function minus_r(A) result(C)
    implicit none
    type(rmatrix), intent(in) :: A
    type(rmatrix) :: C
    ! *** end of interface ***

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(A, 2), C)

    C%m = - A%m
  end function minus_r


  function mult_r_r(A, B) result(C)
    use f77_scalapack, only: pdgemm, blacs_gridinfo
    implicit none
    type(rmatrix), intent(in) :: A, B
    type(rmatrix) :: C
    ! *** end of interface ***

    integer :: nprow, npcol, myrow, mycol ! BLACS

    double precision, parameter :: ONE = 1.0D+0
    double precision, parameter :: zero = 0.0D+0

    !ASSERT(size(A,2)==size(B,1))
    !ASSERT(A%mpi_comm) == B%mpi_comm)

    if ( A%desc(ctxt_) /= B%desc(ctxt_) ) then
       stop "Matrix contexts differ!"
    endif

    !
    ! FIXME: comment on what this does or supposed to do on workers
    ! that are not involved:
    !
    call allocate_rmatrix(pm_get_ctxt(A), size(A, 1), size(B, 2), C)

    !
    ! FIXME: Is this necessary, even after the call to PBLAS was fixed
    ! with proper beta, or just "for the case"?
    !
    C%m = 0.0_RK

    !
    ! FIXME: when two contexts are the same, why double checking?
    !
    call blacs_gridinfo (A%desc(ctxt_), nprow, npcol, myrow, mycol)
    if(myrow < 0) return ! MR_DBG: vielleicht durch Routine von rmatrix_types ersetzen
    call blacs_gridinfo (B%desc(ctxt_), nprow, npcol, myrow, mycol)
    if(myrow < 0) return ! MR_DBG: vielleicht durch Routine von rmatrix_types ersetzen

    ! PDGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
    call pdgemm('N', 'N', A%desc(m_), B%desc(n_), A%desc(n_), &
         ONE, A%m, 1, 1, A%desc, B%m, 1, 1, B%desc, zero, C%m, 1, 1, C%desc)
  end function mult_r_r


  function mult_r_d(A, B) result(C)
    use f77_scalapack, only: blacs_gridinfo, indxl2g
    implicit none
    type(rmatrix), intent(in) :: A
    type(rdmatrix), intent(in) :: B
    type(rmatrix) :: C
    ! *** end of interface ***

    ! BLACS variables
    integer :: myrow, mycol
    integer :: nprow, npcol

    integer :: m, n, m_loc, n_loc, j_loc, bs

    m = size(A, 1)
    n = size(A, 2)

    !ASSERT(n==size(B))

    call allocate_rmatrix(pm_get_ctxt(A), m, n, C)

    call blacs_gridinfo (A%desc(ctxt_), nprow, npcol, myrow, mycol)
    bs = pm_get_bs()

    m_loc = size(C%m, 1)
    n_loc = size(C%m, 2)

    forall(j_loc = 1:n_loc)
       C%m(1:m_loc, j_loc) = A%m(1:m_loc, j_loc) *&
            & B%d(indxl2g (j_loc, bs, mycol, 0, npcol))
    end forall
  end function mult_r_d


  function mult_d_r(A, B) result(C)
    use f77_scalapack, only: blacs_gridinfo, indxl2g
    implicit none
    type(rdmatrix), intent(in) :: A
    type(rmatrix), intent(in) :: B
    type(rmatrix) :: C
    ! *** end of interface ***

    ! BLACS variables
    integer :: myrow, mycol
    integer :: nprow, npcol

    integer :: m, n, m_loc, n_loc, i_loc, j_loc, bs

    m = size(B, 1)
    n = size(B, 2)

    !ASSERT(m==size(A))

    call allocate_rmatrix(pm_get_ctxt(B), m, n, C)

    call blacs_gridinfo (B%desc(ctxt_), nprow, npcol, myrow, mycol)
    bs = pm_get_bs()

    m_loc = size(C%m, 1)
    n_loc = size(C%m, 2)

    forall(i_loc = 1:m_loc, j_loc = 1:n_loc)
          C%m(i_loc, j_loc) =&
               & A%d(indxl2g (i_loc, bs, myrow, 0, nprow)) *&
               & B%m(i_loc, j_loc)
    end forall
  end function mult_d_r


  function mult_d_r_d(a, B, c) result(D)
    implicit none
    type(rdmatrix), intent(in) :: a, c
    type(rmatrix), intent(in) :: B
    type(rmatrix) :: D
    ! *** end of interface ***

    D = a * B * c
  end function mult_d_r_d


  function rpt(X, e) result(Y)
    use f77_scalapack, only: blacs_gridinfo, indxl2g
    implicit none
    type(rmatrix), intent(in) :: X
    type(rdmatrix), intent(in) :: e
    type(rmatrix) :: Y
    ! *** end of interface ***

    ! BLACS variables
    integer :: myrow, mycol
    integer :: nprow, npcol

    integer :: n, bs
    integer :: m_loc, n_loc, i_loc, j_loc

    n = size(e)

    !ASSERT(n==size(X,1))
    !ASSERT(n==size(X,2))

    call allocate_rmatrix(pm_get_ctxt(X), n, n, Y)

    call blacs_gridinfo (X%desc(ctxt_), nprow, npcol, myrow, mycol)
    bs = pm_get_bs()

    m_loc = size(X%m, 1)
    n_loc = size(X%m, 2)

    forall(i_loc = 1:m_loc, j_loc = 1:n_loc)
       Y%m(i_loc, j_loc) = X%m(i_loc, j_loc) / &
            &(E%d(indxl2g (i_loc, bs, myrow, 0, nprow)) + &
            &E%d(indxl2g (j_loc, bs, mycol, 0, npcol)))
    end forall
  end function rpt


  function tr(A) result(C)
    !
    ! Transpose a rmatrix.
    !
    use f77_scalapack, only: pdtran, blacs_gridinfo
    implicit none
    type(rmatrix), intent(in) :: A
    type(rmatrix) :: C
    ! *** end of interface ***

    ! local variables
    integer :: nprow, npcol, myrow, mycol ! BLACS
    integer :: m_C, n_C
    double precision, parameter :: ONE = 1.0D+0
    double precision, parameter :: ZERO = 0.0D+0

    n_C = size(A, 1)
    m_C = size(A, 2)

    call allocate_rmatrix(pm_get_ctxt(A), size(A, 2), size(A, 1), C)

    call blacs_gridinfo (A%desc(ctxt_), nprow, npcol, myrow, mycol)
    if(myrow < 0) return ! MR_DBG: vielleicht durch Routine von rmatrix_types ersetzen

    call pdtran (m_C, n_C, ONE, A%m, 1, 1, A%desc, ZERO, C%m, 1, 1, C%desc)
  end function tr


  pure function pow_d_i(a, p) result(ap)
    !
    !           p
    ! Computes a  for diagonal a.
    !
    implicit none
    type(rdmatrix), intent(in) :: a
    integer, intent(in) :: p
    type(rdmatrix) :: ap
    ! *** end of interface ***

    ap = matrix(array(a)**p)
  end function pow_d_i


  pure function size_d(a) result(n)
    !
    ! Size of a diagonal matrix, here length of a diagonal.
    !
    implicit none
    type(rdmatrix), intent(in) :: a
    integer :: n
    ! *** end of interface ***

    n = size(a%d)
  end function size_d


  pure function size_r_i(a, axis) result(n)
    !
    ! Dimension of an rmatrix along the axis.
    !
    implicit none
    type(rmatrix), intent(in) :: a
    integer, intent(in) :: axis
    integer :: n
    ! *** end of interface ***

    select case (axis)
    case (1)
       n = a%desc(m_)
    case (2)
       n = a%desc(n_)
    case default
       !
       ! This is a pure function, we cannot abort here, return invalid
       ! size instead:
       !
       n = -1
    end select
  end function size_r_i

  function sim(V, U) result(VP)
    !
    ! A more efficient implementation is welcome.
    !
    implicit none
    type(rmatrix), intent(in) :: V, U
    type(rmatrix) :: VP
    ! *** end of interface ***

    VP = tr(U) * V * U
  end function sim

  subroutine extract (array, section, nprow, npcol, row, col)
    !
    ! Extract a section from the array.
    !
    implicit none
    real(RK), intent(in) :: array (:, :)
    real(RK), intent(out) :: section (:, :)
    integer, intent(in) :: nprow, npcol
    integer, intent(in) :: row, col ! base-0
    ! *** end of interface ***

    integer :: i, j
    integer :: m, n, m_blocks, n_blocks

    integer :: iarr, jarr
    integer :: isec, jsec
    integer :: isiz, jsiz

    m = size(array, 1)
    n = size(array, 2)
    m_blocks = (m + BLOCKSIZE - 1) / BLOCKSIZE
    n_blocks = (n + BLOCKSIZE - 1) / BLOCKSIZE

    do j = col, n_blocks - 1, npcol
       do i = row, m_blocks - 1, nprow
          !
          ! Compute  base-0  offsets  into array/section  and  the
          ! block sizes:
          !
          iarr = i * BLOCKSIZE
          jarr = j * BLOCKSIZE

          isec = i / nprow * BLOCKSIZE
          jsec = j / npcol * BLOCKSIZE

          isiz = min(BLOCKSIZE, m - iarr)
          jsiz = min(BLOCKSIZE, n - jarr)

          section (isec + 1:isec + isiz, jsec + 1:jsec + jsiz) = &
               array (iarr + 1:iarr + isiz, jarr + 1:jarr + jsiz)
       enddo
    enddo
  end subroutine extract


  subroutine insert (section, array, nprow, npcol, row, col)
    !
    ! Copy the data into proper section of the array.
    !
    implicit none
    real(RK), intent(in) :: section (:, :)
    real(RK), intent(inout) :: array (:, :)
    integer, intent(in) :: nprow, npcol
    integer, intent(in) :: row, col ! base-0
    ! *** end of interface ***

    integer :: i, j
    integer :: m, n, m_blocks, n_blocks

    integer :: iarr, jarr
    integer :: isec, jsec
    integer :: isiz, jsiz

    m = size(array, 1)
    n = size(array, 2)
    m_blocks = (m + BLOCKSIZE - 1) / BLOCKSIZE
    n_blocks = (n + BLOCKSIZE - 1) / BLOCKSIZE

    do j = col, n_blocks - 1, npcol
       do i = row, m_blocks - 1, nprow
          !
          ! Compute  base-0  offsets  into array/section  and  the
          ! block sizes:
          !
          iarr = i * BLOCKSIZE
          jarr = j * BLOCKSIZE

          isec = i / nprow * BLOCKSIZE
          jsec = j / npcol * BLOCKSIZE

          isiz = min(BLOCKSIZE, m - iarr)
          jsiz = min(BLOCKSIZE, n - jarr)

          array (iarr + 1:iarr + isiz, jarr + 1:jarr + jsiz) = &
               section (isec + 1:isec + isiz, jsec + 1:jsec + jsiz)
       enddo
    enddo
  end subroutine insert


end module matrix_parallel

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
module runrec_module

	use options_module
	use time_module
	use mpi

	implicit none

!----------------------------------------------------------------------------------
!# Module data

	!# Array with desired processor configurations (automatically generated)
	integer, private, allocatable :: procConf(:)


	!# MPI variables
	integer, private :: np, myrank


!----------------------------------------------------------------------------------
!# Public Module subroutines

contains

	subroutine rt_setup
		!# Data
		!--------------------------------------------------------
		integer :: i, offset, add, ierr
		integer :: procConfSize, minProcConf, maxProcConf		! Number of processor configurations
		character(len=50) :: inputArgument

		!# Code
		!--------------------------------------------------------


		!# Initialize module parameters
		call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
		call MPI_Comm_size( MPI_COMM_WORLD, np, ierr )


		!--------- Reading command line arguments ----------
		if( iargc() .ne. 6 ) then
			if( myrank .eq. 0 ) write( unit=0, fmt=* ) 'ERROR: Wrong number of command line arguments. Aborting...'

			call MPI_Finalize( ierr )
			call exit(-1)
		end if

		do i=1, iargc()
			call getarg(i, inputArgument)
			if( index( inputArgument, 'minprocconf' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%minProcConf )
			else if( index( inputArgument, 'maxprocconf' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%maxProcConf )
			else if( index( inputArgument, 'blocksize' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%blocksize )
			else if( index( inputArgument, 'minmatsize' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%minMatSize )
			else if( index( inputArgument, 'maxmatsize' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%maxMatSize )
			else if( index( inputArgument, 'stepsize' ) .gt. 0 ) then
				call readInputArgument( inputArgument, rtOpt_Obj%stepsize )
			else
				if( myrank .eq. 0 ) write( unit=0, fmt=* ) 'ERROR: Wrong command line argument: ', trim(inputArgument), '. Aborting...'

				call MPI_Finalize( ierr )
				call exit(-1)
			end if
		end do


		maxProcConf = min(np, rtOpt_Obj%maxProcConf)
		minProcConf = max(1, rtOpt_Obj%minProcConf)
		if( np .lt. minProcConf ) then
			write(*,*) 'ERROR: Minimum processor count greater than available MPI processes. Aborting...'

			call MPI_Finalize( ierr )
			call exit(-1)
		end if
		if( minProcConf .gt. maxProcConf ) then
			write(*,*) 'ERROR: Minimum processor count greater than maximum processor count. Aborting...'

			call MPI_Finalize( ierr )
			call exit(-1)
		end if



		!# Compute how many processor configurations there are and allocate corresponding array
		procConfSize = min( maxProcConf, 9 - minProcConf )
		if( procConfSize .lt. 0 ) procConfSize = 0
		if( maxProcConf > 9 ) then
                   procConfSize = procConfSize + floor( sqrt( real(maxProcConf) ) )&
                        & - max( floor( sqrt( real(minProcConf) ) ), 2 )
                end if
		if( sqrt( real(minProcConf) )  .gt.  2 ) then
			if( floor( sqrt( real(minProcConf) ) ) .eq. ceiling( sqrt( real(minProcConf) ) ) ) procConfSize = procConfSize + 1
		end if


		allocate( procConf(procConfSize) )

		!# Initialize array with processor configurations
		if( min(8, maxProcConf)-minProcConf+1 .gt. 0 ) then
			do i=1, min(8, maxProcConf)-minProcConf+1
				procConf(i) = i+minProcConf-1
			end do
			offset = i-1
		else
			offset=0
		end if


		if( maxProcConf .gt. 8 ) then
			if( offset .ne. 0 ) then
				add = 2
			else
				add = ceiling( sqrt( real(minProcConf) ) ) - 1
			end if


			do i=1, procConfSize-offset
				procConf(i+offset) = (i + add )**2
			end do
		end if

	end subroutine rt_setup


	subroutine rt_deallocate

		deallocate( procConf )

	end subroutine rt_deallocate



!----------------------------------------------------------------------------------
!# Private Module Main Subroutines

	subroutine rt_runRecord
		!# Data
		!--------------------------------------------------------
		integer :: i, ii
		integer :: thisProcConf				! Processor configuration in each iteration
		integer :: blocksize
		integer :: stepsize
		!integer :: matSize							! Matrix size in each iteration
		integer :: minMatSize, maxMatSize			! Smallest and greatest possible matrix size
		integer :: npRow, npCol						! Number of rows and columns in the BLACS process grid
		integer :: matSizesSize						! Size of array with matrix sizes
		integer :: numMatIter						! Number of matrices which can be processed at once
		integer :: numMatThisIter					! Temporary number of matrices in a iteration step
		integer :: thisMatIndex						! Index in global array matSizes for this process in an iteration step
		integer :: thisMatSize						! Matrix size for one iteration (from array matSizes)
		!integer :: locProcRank, locBlacsRankRow, locBlacsRankCol	! Blacs rank for this process in an iteration step
		integer :: blacsContext
		integer :: myRankRow, myRankCol
		integer :: nLocRow, nLocCol

		!# Array with desired matrix sizes (automatically generated)
		integer, dimension(:), allocatable :: matSizes

		!# Data Matrices for the generalized eigenproblem  A * z = lambda * B * z
		real(kind=r8_kind), allocatable :: A(:,:), B(:,:)

		!# MPI variables
		integer :: ierr


		!# Time evaluation
		type(timeObj) :: t1, t2
		character(len=8) :: start_date
		character(len=10) :: start_time

		character(len=8) :: end_date
		character(len=10) :: end_time

		!# External functions
		integer :: numroc


		!# Code
		!--------------------------------------------------------
		blocksize = rtOpt_Obj%blocksize
		stepsize = rtOpt_Obj%stepsize


		!# Iterate over each "procConf" element. Perform runtime measurements to each processor configuration
		do i=1, size(procConf)

			!# Processes are synchronized
			call MPI_Barrier( MPI_COMM_WORLD, ierr )

			!# Processor configuration for this iteration step
			thisProcConf = procConf(i)
			!write(*,*) 'Prozess ', myrank, ' thisProcConf1: ', thisProcConf

			!# Get BLACS process grid dimensions
			call blacsProcGridDim( thisProcConf, npRow, npCol )


			!# Determine smallest and greatest possible matrix size
			minMatSize = max(npCol*blocksize, rtOpt_Obj%minMatSize)
			maxMatSize = rtOpt_Obj%maxMatSize

			matSizesSize = (maxMatSize - minMatSize) / stepsize + 1

			allocate( matSizes(matSizesSize) )

			do ii=1, matSizesSize
				matSizes(ii) = minMatSize + (ii-1)*stepsize
			end do


			!# Maximum number of matrices which can be processed in one iteration step
			numMatIter = np / thisProcConf

			!#
			do ii=1, matSizesSize, numMatIter

				!# Processes are synchronized
				call MPI_Barrier( MPI_COMM_WORLD, ierr )


				!# Number of matrices which are processed in THIS iteration step
				numMatThisIter = min(matSizesSize - ii + 1, numMatIter)


				!# Index for array "matSizes" for this process in this iteration step
				thisMatIndex = ii + myrank / thisProcConf


				!# Matrix size for this process in this iteration step
				if( thisMatIndex .le. matSizesSize ) then
					thisMatSize = matSizes( thisMatIndex )
				else
					thisMatSize = -1
				end if


				!write(*,*) 'Prozess ', myrank, ' thisProcConf2: ', thisProcConf
				if( thisProcConf == 1 ) then
					!# Sequential processing, use LAPACK
					!------------------------------------------------------------------

					!# Only continue if this process is involved in matrix processing
					if( myrank < numMatThisIter * thisProcConf ) then


						allocate( A(thisMatSize, thisMatSize) )
						allocate( B(thisMatSize, thisMatSize) )

						myRankRow = 0
						myRankCol = 0

						call fillMatrices( A, B, thisProcConf, myRankRow, myRankCol, blocksize )

						t1 = timeconvert( start_date, start_time )

						call execLapack( A, B, thisMatSize )

						t2 = timeconvert( end_date, end_time )


						write(*,*) "-------------------------------------------------------------------------"
						write(*,*) "Matrix dimension: ", thisMatSize, "; Process count: ", thisProcConf
						write(*,*) "Processes involved: ", myrank
						write(*, fmt='(X, A, F6.3)' ) "Processing time: ", t2%t - t1%t
						write(*,*) "Start date and time: ", start_date, "; ", start_time
						write(*,*) "End date and time: ", end_date, "; ", end_time
						!write(*, fmt='(X, A, I3, X, I6, F7.3)') "Code: ", thisProcConf, thisMatSize, t2%t - t1%t
						write(*,*) "-------------------------------------------------------------------------"


						deallocate( A )
						deallocate( B )

					end if		! Only continue if this process is involved in matrix processing


				else
					!# Parallel processing, use ScaLAPACK
					!------------------------------------------------------------------

					!# Get a BLACS context for this run
					call generateBlacsContext( myrank, thisProcConf, numMatThisIter, blacsContext )


					!# Only continue if this process is involved in matrix processing
					if( thisMatSize .gt. 0  .and.  myrank < numMatThisIter * thisProcConf ) then

						call blacs_gridinfo( blacsContext, npRow, npCol, myRankRow, myRankCol )

						!call getLocalMatSize( thisMatSize, blocksize, thisProcConf, myRankRow, myRankCol, nLocRow, nLocCol )
						nLocRow = numroc( thisMatSize, blocksize, myRankRow, 0, npRow )
						nLocCol = numroc( thisMatSize, blocksize, myRankCol, 0, npCol )

						allocate( A(nLocRow, nLocCol) )
						allocate( B(nLocRow, nLocCol) )

						call fillMatrices( A, B, thisProcConf, myRankRow, myRankCol, blocksize )

						t1 = timeconvert( start_date, start_time )

						call execScaLapack( A, B, thisMatSize, blocksize, blacsContext )

						t2 = timeconvert( end_date, end_time )


						if( myRankRow .eq. 0  .and.  myRankCol .eq. 0 ) then
							write(*,*) "-------------------------------------------------------------------------"
							write(*,*) "Matrix dimension: ", thisMatSize, "; Process count: ", thisProcConf
							write(*,*) "Processes involved: ", myrank
							write(*, fmt='(X, A, F7.3)' ) "Processing time: ", t2%t - t1%t
							write(*,*) "Start date and time: ", start_date, "; ", start_time
							write(*,*) "End date and time: ", end_date, "; ", end_time
							write(*, fmt='(X, A, I3, X, I6, F7.3)') "Code: ", thisProcConf, thisMatSize, t2%t - t1%t
							write(*,*) "-------------------------------------------------------------------------"
						end if

						deallocate( A )
						deallocate( B )

						!# The generated BLACS handle is freed but MPI continues afterwards (non-neg. argument)
						call blacs_gridexit( blacsContext )

					end if		! Only continue if this process is involved into matrix processing

				end if		! ( thisProcConf == 1 )

			end do

			deallocate( matSizes )

		end do
	end subroutine rt_runRecord



!----------------------------------------------------------------------------------
!# Private Module Helper Subroutines

	subroutine blacsProcGridDim( nproc, nprow, npcol )
		!# In-Output parameter
		integer, intent(in) :: nproc
		integer, intent(out) :: nprow, npcol

		if( nproc .lt. 9 ) then
			!# Less than 9 processors alloted to the task, generate 1-dimensiona process grid
			nprow = 1; npcol = nproc
		else
			!# 9 processors or more alloted, generate square process grid
			nprow = sqrt( real( nproc ) )
			npcol = nprow
		end if

	end subroutine blacsProcGridDim


	subroutine generateBlacsContext( myrank, procConf, nMatrices, blacsContext )
		!# In-Output parameter
		integer, intent(in) :: myrank, procConf, nMatrices
		integer, intent(inout) :: blacsContext   !!! nur out?????

		!# Local variables
		integer :: i, j
		integer :: tmpBlacsContext, myLocRank, npRow, npCol, myRankRow, myRankCol
		integer :: lowerProcRanks	! Number of processes with higher and lower rank
		integer, dimension(:,:), allocatable ::  blacsProcGrid



		!# Extract BLACS process grid grid dimensions
		call blacsProcGridDim( procConf, npRow, npCol )

		!# Allocate BLACS process grid
		allocate( blacsProcGrid(npRow, npCol) )


		myLocRank = mod( myrank, procConf )
		!lowerProcRanks = myLocRank
		!higherProcRanks = procConf - lowerProcRanks - 1

		do i=1, nMatrices
			!blacsProcGrid = reshape((/j,j=(i-1)*procConf,procConf*i-1/),(/npRow,npCol/))
			blacsProcGrid = reshape( (/(j,j=(i-1)*procConf,procConf*i-1)/),(/npRow,npCol/) )

			call blacs_get( -1, 0, tmpBlacsContext )

			call blacs_gridmap( tmpBlacsContext, blacsProcGrid, npRow, npRow, npCol )

			if( myrank .ge. (i-1)*procConf  .and.  myrank .le. procConf*i-1 ) blacsContext = tmpBlacsContext
		end do

		deallocate(blacsProcGrid)

	end subroutine generateBlacsContext


	!!!!!!!!!!! Parameterliste überpruefen !!!!!!
	subroutine getLocalMatSize( globMatSize, blocksize, procConf, myRankRow, myRankCol, nLocRow, nLocCol )
		!# In-Output parameter
		integer, intent(in) :: globMatSize, blocksize, procConf, myRankRow, myRankCol
		integer, intent(out) :: nLocRow, nLocCol
		!integer, intent(out), optional :: nBlocksRow_out, nLocBlocksRow_out

		!# Local variables
		integer :: i, npRow, npCol
		integer :: nBlocksRow, nLocBlocksRow
		integer :: nBlocksCol, nLocBlocksCol


		call blacsProcGridDim( procConf, npRow, npCol )

		nBlocksRow = int( ceiling( real(globMatSize) / real(blocksize) ) )
		nBlocksCol = nBlocksRow

		nLocRow = 0
		nLocCol = 0

		do i=myRankRow+1, nBlocksRow, npRow
			if( i*blocksize .gt. globMatSize ) then
				nLocRow = nLocRow + globMatSize - (i-1)*blocksize
			else
				nLocRow = nLocRow + blocksize
			end if
		end do

		do i=myRankCol+1, nBlocksCol, npCol
			if( i*blocksize .gt. globMatSize ) then
				nLocCol = nLocCol + globMatSize - (i-1)*blocksize
			else
				nLocCol = nLocCol + blocksize
			end if
		end do

	end subroutine getLocalMatSize


	subroutine fillMatrices( A, B, procConf, myRankRow, myRankCol, blocksize )
		!# In-Output parameter
		real(kind=r8_kind), dimension(:,:), allocatable, intent(inout) :: A, B
		integer, intent(in) :: procConf, myRankRow, myRankCol, blocksize
		!integer, dimension(9), intent(out), optional :: desca, descb


		!# Local variables
		integer :: i, j, rows, cols, npRow, npCol
		integer :: i_block, j_block, thisBlocksizeRow, thisBlocksizeCol
		integer :: globBlockIndexRow, globBlockIndexCol, nLocBlocksRow, nLocBlocksCol

		!# RNG variables
		integer, dimension(:), allocatable :: seedArr
		integer :: n, seedSize, seedSca
		real(kind=r8_kind) :: dummy


		!# Seed the random number generator
		call random_seed( size=seedSize )
		allocate( seedArr(seedSize) )
		call system_clock( count=seedSca )
		seedSca = seedSca + myrank
		!seedArr = seedSca + 37 * (/ (i - 1, i = 1, n) /)

		do i=1, seedSize
			seedArr(i) = seedSca * myrank + i
		end do

		!write(*,*) seedArr

		call random_seed( put=seedArr )


		!# Fill matrix A
		call random_number( A )


		!# Fill matrix B
		rows = size( B, 1 )
		cols = size( B, 2 )

		call blacsProcGridDim( procConf, npRow, npCol )

		nLocBlocksRow = int( ceiling( real(rows) / real(blocksize) ) )
		nLocBlocksCol = int( ceiling( real(cols) / real(blocksize) ) )

		do j_block=1, nLocBlocksCol
			do i_block=1, nLocBlocksRow
				globBlockIndexRow = myRankRow + npRow*(i_block-1)
				globBlockIndexCol = myRankCol + npCol*(j_block-1)

				if( i_block*blocksize .gt. rows ) then
					thisBlocksizeRow = rows - blocksize*(i_block-1)
				else
					thisBlocksizeRow = blocksize
				end if

				if( j_block*blocksize .gt. cols ) then
					thisBlocksizeCol = cols - blocksize*(j_block-1)
				else
					thisBlocksizeCol = blocksize
				end if


				do j=1, thisBlocksizeCol
					do i=1, thisBlocksizeRow

						dummy = 0_r8_kind

						if( globBlockIndexCol .eq. globBlockIndexRow  .and.  i .eq. j ) call random_number(dummy)

						B( (i_block-1)*blocksize + i, (j_block-1)*blocksize + j ) = dummy

					end do
				end do
			end do
		end do

		deallocate( seedArr )

	end subroutine fillMatrices



	subroutine readInputArgument( characterArgument, numericalArgument )
		!# In-Output parameter
		character(len=*), intent(in)  :: characterArgument
		integer,          intent(out) :: numericalArgument

		!# Local variables
		integer :: equalIndex


		equalIndex = index( characterArgument, '=', .true. )
		read( unit=characterArgument(equalIndex+1:), fmt='(I10.1)' ) numericalArgument


	end subroutine readInputArgument



	subroutine execLapack( A, B, n )
		!# In-Output parameter
		real(kind=r8_kind), allocatable, intent(inout) :: A(:,:), B(:,:)
		integer, intent(in) :: n


		!# DSYGV arguments
		real(kind=r8_kind), dimension(:), allocatable :: W, work
		real(kind=r8_kind), dimension(1) :: workdummy

		!# DSYGV input parameters
		character, parameter :: jobz = 'V', uplo = 'L'
		integer, parameter :: itype = 1
		integer :: lda, ldb, info, lwork = -1

		!# Time evaluation
		type(timeObj) :: t1, t2


		lda = n; ldb = n
		allocate( W( n ) )
		lwork = -1
		call dsygv( itype, jobz, uplo, n, A, lda, B, ldb, W, workdummy, lwork, info )


		lwork = workdummy(1)
		allocate( work( lwork ) )

		t1 = timeconvert()

		call dsygv( itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, info )

		t2 = timeconvert()

		write(*, fmt='(X, A, I3, X, I6, F7.3)') "Code: ", 1, n, t2%t - t1%t

		deallocate(W); deallocate( work )

	end subroutine execLapack


	subroutine execScaLapack( A, B, n, blocksize, blacsContext )
		!# In-Output parameter
		real(kind=r8_kind), dimension(:,:), allocatable, intent(in) :: A, B
		integer, intent(in) :: n, blocksize, blacsContext


		!# Routine variables
		integer :: i

		!# Matrix variables
		integer :: lld, nLocRows, nLocCols

		!# BLACS variables
		integer :: npRow, npCol, myRankRow, myRankCol

		!# ScaLAPACK parameters
		integer, parameter :: ibtype=1, ia=1, ja=1, ib=1, jb=1, iz=1, jz=1, vl=1, vu=1, il=1, iu=1
		character(len=1), parameter :: eigrange='A', jobz='V', uplo='L'
		real(kind=r8_kind), parameter :: abstol=0._r8_kind, orfac=0._r8_kind

		!# ScaLAPACK variables
		real(kind=r8_kind), dimension(:,:), allocatable :: lZ
		real(kind=r8_kind), dimension(:), allocatable :: W, work, gap
		real(kind=r8_kind) :: workdummy(10)

		integer, dimension(:), allocatable :: iwork, ifail, iclustr
		integer :: m, nz, lwork=-1, liwork=-1, iworkdummy(10), info
		integer, dimension(9) :: desca, descb, descz

		!# External functions
		integer :: numroc



		call blacs_gridinfo( blacsContext, npRow, npCol, myRankRow, myRankCol )

		nLocRows = numroc( n, blocksize, myRankRow, 0, npRow )
		nLocCols = numroc( n, blocksize, myRankCol, 0, npCol )
		lld = nLocRows

		!# local Z
		allocate( lZ( nLocRows,nLocCols ), stat=info )
		if( info .ne. 0 ) then
			write (*,*) 'Allocation error. Status ', info
			call blacs_abort( blacsContext, -1 )
			call exit( -1 )
		end if

		!write(*,*) "Myrank = ", myrank, "; n = ", n, "; blocksize = ", blocksize, "; blacsContext = ", blacsContext,&
		!	& "; lld = ", lld, "; vor"
 		call descinit( desca, n, n, blocksize, blocksize, 0, 0, blacsContext, lld, info )
		call descinit( descb, n, n, blocksize, blocksize, 0, 0, blacsContext, lld, info )
		call descinit( descz, n, n, blocksize, blocksize, 0, 0, blacsContext, lld, info )
		!write(*,*) "Myrank = ", myrank, "; nach; ", desca


		!# Workspace allocation
		allocate( iclustr( 2*nprow*npcol ) ); allocate( gap( nprow*npcol ) )
		allocate( W(n), stat=info ); allocate( ifail(n) )
		if( info .ne. 0 ) then
   			write (*,*) 'Allocation error. Status ', info
   			call blacs_abort( blacsContext, -1 )
   			call exit( -1 )
		end if

		do i=1, 10
			workdummy(i)=0._r8_kind; iworkdummy(i)=0
		end do


		!# First run for workspace querry
		lwork=-1; liwork=-1
		call pdsygvx( ibtype, jobz, eigrange, uplo, n, A, ia, ja, desca, B, ib, jb, descb,&
	     	& vl, vu, il, iu, abstol, m, nz, W, orfac, lZ, iz, jz, descz, workdummy, lwork,&
     		& iworkdummy, liwork, ifail, iclustr, gap, info )
		lwork = int(workdummy(1)); liwork = iworkdummy(1)

		allocate( iwork( liwork ), stat=info )
		if( info .ne. 0 ) then
   			write (*,*) 'Allocation error. Status ', info
   			call blacs_abort( blacsContext, -1 )
			call exit( -1 )
		end if
		allocate( work( lwork ), stat=info )
		if( info .ne. 0 ) then
			write (*,*) 'Allocation error. Status ', info
			call blacs_abort( blacsContext, -1 )
			call exit( -1 )
		end if

		!write(*,*)

		call pdsygvx( ibtype, jobz, eigrange, uplo, n, A, ia, ja, desca, B, ib, jb, descb,&
			& vl, vu, il, iu, abstol, m, nz, W, orfac, lZ, iz, jz, descz, work, lwork,&
			& iwork, liwork, ifail, iclustr, gap, info )


		deallocate ( lZ )
		deallocate( work ); deallocate( iwork )
		deallocate( W ); deallocate( iclustr ); deallocate( ifail ); deallocate( gap )


	end subroutine execScaLapack

end module runrec_module

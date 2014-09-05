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
module peigval
  !-------------------------------------------------------------------
  !
  !  Purpose: Distribute a matrices between processors and find
  !  an eigenvalue in optimal time
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: SB, MG
  !  Date: 2.3.2004
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

#include "def.h"
  use type_module, only: IK=>i4_kind, RK=>r8_kind, RK_=>r4_kind ! type specification parameters
  use datatype
  use hamiltonian_module, only: ham_tot,ham_tot_real,ham_tot_imag
  use overlap_module

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

!  public :: arrmat2, arrmat3
!  type(arrmat2),allocatable,public           :: eigval(:)
!  type(arrmat3),allocatable,public           :: eigvec(:)


  !------------ public functions and subroutines ---------------------
  public parallelization_subroutine, eigs_entry 

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine eigs_entry(eigvec, eigval) 
    use symmetry_data_module
    use msgtag_module
    use comm_module
    implicit none
    type(arrmat3), intent(inout) :: eigvec(:) ! (n_irreps)
    type(arrmat2), intent(inout) :: eigval(:) ! (n_irreps)
    ! *** end of interface ***

    type(arrmat2), allocatable :: MATRIXS(:)

    integer(IK) :: n_irrep, alloc_stat, n, i, n_tasks, j, l
    integer(IK),allocatable    :: sorted(:),dim_irrep(:)
    integer(IK) :: info, i_spin, ig

    logical :: i_am_master

    i_spin = 1 ! only that case, must be reconsidered

!!! NON SPIN ORBIT CASE

    n_irrep = ssym%n_spin*ssym%n_irrep
    if (ssym%n_spin == 1) then
      allocate(dim_irrep(n_irrep))
      do n=1,n_irrep
        dim_irrep(n)  = ssym%dim(n)
      enddo
    else
      print *,"ssym%n_spin = ",ssym%n_spin
      allocate(dim_irrep(n_irrep))
      i = 1
      do n=1,ssym%n_irrep
        dim_irrep(i)  = ssym%dim(n)
        dim_irrep(i+1) = dim_irrep(i)
        i = i + 2
      enddo
    end if

    allocate(sorted(ssym%n_irrep),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("EIGS_ENTRY: allocation of sorted failed")

    call sort_irreps(ssym%dim(:), sorted)
   
    i_am_master = comm_i_am_master()

    allocate(MATRIXS(n_irrep))

    i_spin = 1

    if (i_am_master) then
     !! send ham_tot to all processors
       if (comm_get_n_processors() /= 1) then
         do j=2,comm_get_n_processors()
            WARN('FIXME: needs verification')
            ! send data to processor "j"
            call comm_init_send(j, msgtag_packed_message)
            i = 1
            do  l = 1, n_irrep 
              call commpack(ham_tot(i)%m(1,1,i_spin),ssym%dim(i)*ssym%dim(i),1,info)
              allocate(MATRIXS(l)%m(dim_irrep(l),dim_irrep(l)))
              MATRIXS(l)%m(1:dim_irrep(l),1:dim_irrep(l)) = ham_tot(i)%m(1:ssym%dim(i),1:ssym%dim(i),i_spin)
              if((ssym%n_spin /= 1 .and. i_spin == 2) .or. ssym%n_spin == 1) then 
                i = i + 1
                i_spin = 1
              else
                i_spin = i_spin + 1
              end if
              end do
            call comm_send()
         end do
       else
         i = 1
         do  l = 1, n_irrep
            allocate(MATRIXS(l)%m(dim_irrep(l),dim_irrep(l)))
            MATRIXS(l)%m(1:dim_irrep(l),1:dim_irrep(l)) = ham_tot(i)%m(1:ssym%dim(i),1:ssym%dim(i),i_spin)
            if((ssym%n_spin /= 1 .and. i_spin == 2) .or. ssym%n_spin == 1) then 
              i = i + 1
              i_spin = 1
            else
              i_spin = i_spin + 1
            end if
         end do
       end if
    else
          WARN('FIXME: needs verification')
          ! recv data from master:
          call comm_save_recv(comm_master_host, msgtag_packed_message)
          do  l = 1, n_irrep
            print *,"RECEIVE: dim_irrep(l)=",dim_irrep(l)
            allocate(MATRIXS(l)%m(dim_irrep(l),dim_irrep(l)))
            call communpack(MATRIXS(l)%m(1,1),dim_irrep(l)*dim_irrep(l),1,info)
!            print *,"l=",l
          end do
    endif

    print *,"EIGS_ENTRY: # of IRREPS",n_irrep
    print *,"EIGS_ENTRY: SORTED IRREPS",sorted
    print *,"EIGS_ENTRY: DIMENSIONAL OF IRREPS",dim_irrep

    call read_overlap()

    call parallelization_subroutine(eigvec, eigval, ssym%n_spin, MATRIXS, overlap, dim_irrep, n_irrep)

!    call dealloc_overlap()

    deallocate(sorted,STAT=alloc_stat)
    deallocate(dim_irrep,STAT=alloc_stat)    
    deallocate(MATRIXS,STAT=alloc_stat)

    if (alloc_stat.ne.0) call error_handler &
         ("EIGS_ENTRY: deallocation of sorted failed")
  end subroutine eigs_entry

  !*************************************************************
  subroutine parallelization_subroutine(eigvec, eigval, n_spin, MATRIXS, overlap, dim_irrep, N )
    !  Purpose: 
    !    This subroutine obtain optimal map of distribution
    !    matrices by processors from algorithm_module, reorganized it
    !    and solved eigenvalues via mirceas_solver.
    !    
    !  On enter: 
    !       MATRIXS is a set of matrices to be solved
    !       DIM is a set of dimensions of matrices
    !       N is a number of matrices to be solved
    !
    !  On exit:
    !       EIGENVALUES is a set of eigenvalues
    !
    !  Prediction time and real time will be obtained and showed too
    !
    !------------ Modules used ------------------- ---------------
    use algorithm_main
    use linux_part
    use datatype, only: arrmat2, arrmat3
    use comm_module
    use msgtag_module
    implicit none
    type(arrmat3), intent(inout) :: eigvec(:) ! (n_irreps)
    type(arrmat2), intent(inout) :: eigval(:) ! (n_irreps)
    ! *** end of interface ***

      integer(IK), intent(in) :: n,  dim_irrep(:)
      type(arrmat2) :: overlap(:)
      type(arrmat2), intent(in) :: MATRIXS(:)
      integer(IK), intent(in) ::  n_spin

      integer(IK) :: step, info, p
      integer(kind=IK), allocatable, dimension(:) :: input_vector
      integer(kind=IK), allocatable, dimension(:,:) :: distribution, final_result
      real(RK) :: final_time
      real(RK), allocatable :: prediction_time_vector(:)
      logical :: not_that_case,loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for eigensolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(kind=RK), allocatable, dimension(:,:) :: A, B, Z
      real(kind=RK), allocatable, dimension(:) :: W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for grid_mapping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(IK) :: alloc_stat,i , k, l, f_loop, procID, j, s, nn
      integer(IK) :: IAM
      integer(IK) :: icontext0
      integer(IK),allocatable :: map(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Variable - queue_matrix (#_of_proc,#_of queue) for the queue of the task
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     integer(IK), allocatable :: queue_matrix(:,:), final_result_correct(:,:)
     real(RK), allocatable :: time_vector1(:),time_vector2(:)
     real(RK) :: vspom
     logical :: flag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Variables for timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(IK) :: time_array_0(8), time_array_1(8)
      real(RK) ::  start_time, finish_time

      integer(IK) :: i_spin, ig, non_zero, slave_finished
      
!      ASSERT(allocated(eigvec))
!      ASSERT(allocated(eigval))

      allocate(input_vector(n),stat=alloc_stat)

      allocate(distribution(n,n),final_result(n,n), stat=alloc_stat)

!      print *,iam,'call mpi_init(ierr)'
!      call MPI_Init(ierr)
!      print *,iam,'1) ierr=',ierr

      print *,iam,'call BLACS_PINFO()'
      call BLACS_PINFO( IAM, P )
      print *,iam,'iam=',iam,' nprocs=',P

      final_result = 0.0
      loop = .true.
      step = 2

      do while(loop)
!      print *,"Processing...."

       not_that_case = .false.
!        print *,iam,"||","||",n,"||",p
        call algorithm1(n, p, dim_irrep, distribution, step)
!        print *,distribution

        do i=1,n
          if (maxval(distribution(1:n,i)) .eq. 0) then
            not_that_case = .true.
          end if
        end do
        if (.not. not_that_case .and. maxval(final_result) .eq. 0) then
           final_result = distribution
           loop = .false.
        end if
        step = step+1
      end do

!      print *,"MIRCEA'S ALG SAID:"

!      do i = 1,n
!        print *,final_result(:,i)
!      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Make a queue_matrix, don't change i
!  (i = # of column in queue_matrix & final_result_correct)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      i = 0
      do j = 1,n
       if ( maxval( final_result(j,:) ) /= 0 )  i = i + 1
      end do

!        print *,"Start allocating ", p, i, n
        allocate ( queue_matrix(p,i), final_result_correct(n,i), stat = alloc_stat )
!        if (alloc_stat == 0) print *,"Allocation correct"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final result matrix without columns of 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      s = 1
      do j = 1,n
        if ( maxval(final_result(j,:)) /= 0 ) then
          final_result_correct(:,s) = final_result(j,:)
          s = s + 1
        end if
      end do

      do j = 1,n
        print *,final_result_correct(j,:)
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      queue_matrix = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All table for queue_matrix, works not properly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do f_loop = 1,i
        s = 1
        do j = 1,n
          k = final_result_correct(j,f_loop)
          if ( k /= 0 ) then
            do l = s,( s + k ) - 1
               queue_matrix(l,f_loop) = j
            end do
            s = s + k
          end if
        end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start ordering
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( i /= 1 ) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form time_vector1 and 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do s=1,i-1
          allocate(time_vector1(p),time_vector2(p), stat = alloc_stat)
          time_vector1 = 0
          time_vector2 = 0
          do f_loop = 1,p
            if (queue_matrix(f_loop,s) /= 0 ) time_vector1(f_loop) = &
              & time_vector1(f_loop) + time_fit(final_result_correct(queue_matrix(f_loop,1),s), &
              & dim_irrep(queue_matrix(f_loop,s)))
            if (queue_matrix(f_loop,s+1) /= 0 ) time_vector2(f_loop) = &
              & time_fit(final_result_correct(queue_matrix(f_loop,2),s+1),&
                              & dim_irrep(queue_matrix(f_loop,s+1)))
          end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find minimal value in time_vector1 and maximal_value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          f_loop = 1
          do while ( maxval(time_vector2) /= 0 )
               if ( time_vector1(f_loop) == minval(time_vector1) ) then
!                 print *,"f=",f
                 flag = .true.
                 do l = 1,p
                   if ( time_vector2(l) == maxval(time_vector2) .and. flag .and. maxval(time_vector2) /= 0) then
                     if ( l /= f_loop ) then
!                       print *,"l=",l
                       time_vector2(l) = time_vector2(f_loop)
                       vspom = queue_matrix(f_loop,s+1)
                       queue_matrix(f_loop,s+1) = queue_matrix(l,s+1)
                       queue_matrix(l,s+1) = vspom
                       flag = .false.
                     end if
                     time_vector2(f_loop) = 0

                   end if
                 end do
                 time_vector1(f_loop) = 9999.9
               end if
          f_loop = f_loop + 1
          if (f_loop == p + 1) f_loop = 1
          end do
        end do
      end if

      print *,"Debug for queue matrix : "
      do j = 1,p
        print *,queue_matrix(j,:)
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prediction time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(prediction_time_vector(p), stat = alloc_stat)

      prediction_time_vector = 0.

      do l = 1,p
        do j = 1,i
          if (queue_matrix(l,j) /= 0 ) prediction_time_vector(l) = &
          prediction_time_vector(l) + time_fit(final_result_correct(queue_matrix(l,j),j), &
          dim_irrep(queue_matrix(l,j)))
        end do
      end do

      final_time = maxval(prediction_time_vector)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Main part: We make a grid
!  decomposition, setting up map grid on all machines
!  N_of_working_matrices is the order number of the matrix with
!      which we are working on the #iam processor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      print *,"Prediction time is ",final_time

      call date_and_time(values=time_array_0)
        start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
          + time_array_0 (7) + 0.001 * time_array_0 (8)

      do s = 1,i ! case should run via all queue_matrix row

      print *,iam,"DEBUG: Step # ",s, "i = ",i

      print *,iam,'call BLACS_GET(1.1)'
      call BLACS_GET( 0, 0, icontext0 )
      print *,iam,'icontext0=',icontext0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Need to choose appropriate input_vecter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        input_vector = final_result_correct(:,s)

         if ( queue_matrix(iam+1,s) /= 0 ) then
           k = input_vector(queue_matrix(iam+1,s))

           if ( k /= 1 ) then
             allocate(map(k),stat = alloc_stat)
             f_loop = 1
             do j = 1,k
!              map(j) = ((queue_matrix(iam+1,s)*k)-k)+(j-1)
               loop = .true.
               l = f_loop
               do while (loop)
                 if (queue_matrix(iam+1,s) == queue_matrix(l,s)) then
                   map(j) = l-1
                   f_loop = l + 1
                   loop = .false.
                 end if
                 l = l + 1
               end do
             end do
             procID = k
           else
             allocate(map(k),stat = alloc_stat)
             map(1) = iam
             procID = 1
           end if

           print *,iam,' || ',procID,' || ',map,' || ',NROW(procID),'||',NCOL(procID)
           print *,iam,'call blacs_gridmap(icontext0)'
           call blacs_gridmap(icontext0,map,NROW(procID),NROW(procID),NCOL(procID))
           print *,iam,'icontext0=',icontext0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting up matrices A,B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         nn = dim_irrep(queue_matrix(iam+1,s))
         print *,"nn=",nn
         allocate(A(1:nn,1:nn),B(1:nn,1:nn),Z(1:nn,1:nn),W(1:nn))
         print *,"Set up A and B"
         A = MATRIXS(queue_matrix(iam+1,s))%m(1:nn,1:nn)
         if (n_spin == 1) then
           B = overlap(queue_matrix(iam+1,s))%m(1:nn,1:nn)
         else
           if (mod(queue_matrix(iam+1,s),2) /= 0) then
             l = (queue_matrix(iam+1,s)+1)/2
           else
             l = queue_matrix(iam+1,s)/2
           end if
            B = overlap(l)%m(1:nn,1:nn)
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now we must start our calculations for eigensolver
!  call eigensolve (icontext0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           print *,iam, 'lets start our calculations', dim_irrep(queue_matrix(iam+1,s)), procID

           print *,iam,' || ', dim_irrep(queue_matrix(iam+1,s)),' || ', procID, ' || ',icontext0
           call mirceas_solver(dim_irrep(queue_matrix(iam+1,s)), A, B, Z, W, procID, icontext0, info)
           print *,iam, 'JOB WELL DONE!!!'

         if (n_spin == 1) then 
           i_spin = 1
           l = queue_matrix(iam+1,s)
           eigval(queue_matrix(iam+1,s))%m(1:nn,i_spin) = W
           eigvec(queue_matrix(iam+1,s))%m(1:nn,1:nn,i_spin) = Z
         else
           if (mod(queue_matrix(iam+1,s),2) /= 0) then
             i_spin = 1
             l = (queue_matrix(iam+1,s)+1)/2
           else
             i_spin = 2
             l = queue_matrix(iam+1,s)/2
           end if
            print *,"l=",l
           eigval(l)%m(1:nn,i_spin) = W
           eigvec(l)%m(1:nn,1:nn,i_spin) = Z
         end if
 
         if ( iam/=0 ) then
           print *,"SEND: ",nn
           call comm_init_send(comm_master_host,msgtag_SendBackEig)
           call commpack(l,1,1,info)
           call commpack(nn,1,1,info)
           call commpack(i_spin,1,1,info)
           call commpack(W(1),nn,1,info)
           call commpack(Z(1,1),nn*nn,1,info)
           call comm_send()
         end if
 
        else

          allocate(map(1),stat = alloc_stat)
          map(1) = iam
          procID = 1

!          print *,iam,'call blacs_gridmap(icontext0)'
          call blacs_gridmap(icontext0,map,1,1,procID)
!          print *,iam,'icontext0=',icontext0

          print *,iam," IGNORED CASE"
     
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finalizing all processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(map,stat = alloc_stat)

      deallocate(A,B,Z,W,stat = alloc_stat)

      print *,iam,'call blacs_gridexit(icontext0)'
      call blacs_gridexit(icontext0)
      print *,iam,'icontext0=',icontext0
      end do

      call   date_and_time(values=time_array_1)
      finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
        + time_array_1 (7) + 0.001 * time_array_1 (8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! All matrices received by master
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      non_zero = 0
      do s = 1,i
        do l = 2,p
          if (queue_matrix(l,s)/=0) non_zero = non_zero + 1
        end do
      end do

         if ( iam == 0 ) then
           do while (non_zero/=0)
             call comm_save_recv(comm_all_other_hosts,msgtag_SendBackEig)
             print *,"non_zero = ",non_zero
             call communpack(l,1,1,info)
             print *,"RECV: ",l
             call communpack(nn,1,1,info)
             print *,"RECV: ",nn
             call communpack(i_spin,1,1,info)
             print *,"RECV: ",i_spin
             call communpack(eigval(l)%m(1,i_spin),nn,1,info)
             call communpack(eigvec(l)%m(1,1,i_spin),nn*nn,1,info)
             non_zero = non_zero - 1
           end do
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!      print *,iam,' Real time is:',finish_time - start_time

!      call blacs_exit(1)

      deallocate(queue_matrix, final_result_correct, stat = alloc_stat)

      deallocate(input_vector,stat=alloc_stat)

      deallocate(distribution,final_result, stat=alloc_stat)

      deallocate(time_vector1,time_vector2, stat = alloc_stat)

      deallocate(prediction_time_vector, stat = alloc_stat)

      print *,iam,' Real time is:',finish_time - start_time
      

  end subroutine parallelization_subroutine
  !*************************************************************
  


  !*************************************************************
  subroutine mirceas_solver( N, A, B, Z, W, P, CONTEXT, INFO )
    !  Purpose: 
    !     =======
    !     This subroutine calls ScaLAPACK or LAPACK to solve
    !      a generalized eigenproblem of the type AZ = (lambda)BZ,
    !      where A is a symmetric matrix and B is a symmetric positive definite
    !      matrix (for real problem) or A is Hermitian and B is Hermitian and
    !      positive definite.
    !
    !     The subroutine should be part of the module "eigen_data_module".
    !
    !
    !      N                Dimension of the matrices
    !      For ScaLAPACK:
    !      NB               Blocking factor to distribute the raws of the array;
    !                        the same value will be used for the blocking factor
    !                        to distribute the culumns of the array. The default
    !                        here is 64.
    !      NPROW            Number of process rows in the process grid to be
    !                        created; optional parameter: present only if
    !                        ScaLAPACK is called.
    !      NPCOL            Number of process culumns in the process grid to be
    !                        created; optional parameter: present only if
    !                        ScaLAPACK is called.
    !                       The total number of processes involved is:
    !                        P = NPROW * NPCOL.
    !
    !     If package = 'ScaLAPACK' either PDSYGVX (real problem) will be called;
    !     if package = 'LAPACK_X' LAPACK driver expert routine will be called; !NOT USED BECAUSE SLOW!
    !     if package = 'LAPACK_D' LAPACK devide & conquer routine will be called.
    !
    !     Routine called by "master_solver".

    !------------ Modules used ------------------- ---------------
    implicit none

!     .. ARGUMENTS ..
      integer(kind=IK), intent(in) :: N, P
      real(kind=RK), dimension(:,:), intent(in) :: A, B
      real(kind=RK), dimension(:,:), intent(out) :: Z
      real(kind=RK), dimension(:), intent(out) :: W
      integer(kind=IK), intent(out) :: info


      integer(kind=IK), parameter :: NB = 64


!     .. Other parameters ..
      integer(kind=IK) :: l, i
      character(3) :: package
      integer(kind=IK), parameter :: DUMMY = 1
      real(kind=RK) :: ZERO = 0.0D+00, &
                       ORFAC = -1.0D+00, &
                       ABSTOL = -1.0D+00
      integer(kind=IK) :: allocate_status
      integer(kind=IK) :: LWORK = 100000, LIWORK = 50000, &
                          LDA = 100, ANB, NPS, SQNPC

!     .. Local Scalars ..
      integer(kind=IK) :: CONTEXT, IAM, IBTYPE = 1, &
                          M, MYPCOL, MYPROW, &
                          NZ, NPROW, NPCOL

!     .. Local Arrays ..
      real(kind=RK), allocatable, dimension(:,:) :: AA, BB
      integer(kind=IK), allocatable, dimension(:) :: IWORK, DESCA, DESCB, &
                                                     DESCZ, ICLUSTR, IFAIL
      real(kind=RK), allocatable, dimension(:) :: GAP, WORK

!     .. ScaLAPACK and LAPACK tool functions ..
      real(kind=RK), external :: PDLAMCH, DLAMCH, PDLANSY
      integer(kind=IK), external :: NUMROC, ICEIL, PJLAENV

      intrinsic          MAX, CEILING, REAL, INT, SQRT

!     .. External Subroutines ..
      external           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO, &
                         DESCINIT, PDSYGVX, DSYGVX, DSYGVD, PDSYTTRD, &
                         BLACS_SETUP

    !------------ Executable code ------------------------------------


      if( p == 1 ) then
         package = 'L_D'
      else
         package = 'S__'
      end if

      if( package == 'S__' ) then

!        Get information about the context
         call BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYPROW, MYPCOL )

!         print *,"DEBUG: ",CONTEXT,"||", NPROW,"||", NPCOL,"||", MYPROW,"||", MYPCOL

!         print *,"start from allocating"
         allocate( WORK(LWORK), IWORK(LIWORK), &
                   ICLUSTR(NPROW*NPCOL*2), &
                   DESCA(9), DESCB(9), DESCZ(9), &
                   IFAIL(N), GAP(NPROW*NPCOL), AA(N,N), BB(N,N), &
                   STAT=allocate_status )

         if( allocate_status /= 0 ) then
!            call error_handler( "solver, PDSYGVX: first allocation &
!                               &of arrays WORK, IWORK, ICLUSTR, &
!                               &DESCA, DESCB, DESCZ, AA, BB, IFAIL and/or &
!                               &GAP failed" )
         end if

         LDA = MAX( 1, NUMROC( N, NB, MYPROW, 0, NPROW ) )
         if( LDA < CEILING( REAL( CEILING( REAL( N )/NB ) )/NPROW ) * NB ) then
             LDA = CEILING( REAL( CEILING( REAL( N )/NB ) )/NPROW ) * NB
         end if
!         print *,"DEBUG: LDA is ",LDA


!         print *,"DESCINIT( DESCA, ...  )"
!        These are basic array descriptors for dense matrices.
         call DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
!         print *,"DESCINIT( DESCA, ...  )",INFO
         call DESCINIT( DESCB, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
!         print *,"DESCINIT( DESCB, ...  )",INFO
         call DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
!         print *,"DESCINIT( DESCZ, ...  )",INFO

         if( INFO /= 0 ) then

!            call error_handler( "solver, PDSYGVX: descriptor arrays could &
!                               &not be initialized" )
            deallocate( WORK, IWORK, ICLUSTR, DESCA, DESCB, DESCZ, IFAIL, &
                        GAP, AA, BB, STAT=allocate_status )
            if( allocate_status /= 0 ) then
!               call error_handler( "solver, PDSYGVX: deallocation &
!                                  &of arrays WORK, IWORK, ICLUSTR, &
!                                  &DESCA, DESCB, DESCZ, AA, BB, IFAIL &
!                                  &and/or GAP failed" )
            end if
         end if
!        Distribute the matrices A and B which will be stored in
!        the local ScaLAPACK arrays AA and BB, respectivly
         do l = 1, N
             do i = 1, N
                     call PDELSET( AA, l, i, DESCA, A(l,i) )
                     call PDELSET( BB, l, i, DESCB, B(l,i) )
             end do
         end do

!        Get a good value for ABSTOL
!        Comment the next line if you want a default value.
         ABSTOL = 2 * PDLAMCH( CONTEXT, 'S' )

!        Get a good value for ORFAC to relax the reorthogonalization
!        Comment the next two lines for a default value for ORFAC
         ORFAC = PDLAMCH( CONTEXT, 'E' ) * PDLANSY( 'O', 'U', N, A, 1, 1, &
                                                    DESCA, WORK )
!        Unomment the next line if reorthogonalization should be skipped
!         ORFAC = ZERO


!          Get optimal values for LWORK and LIWORK.
           LWORK = -1

!           print *,"Finaly it started to calculate"

           call PDSYGVX( IBTYPE, 'V', 'A', 'U', N, AA, 1, 1, DESCA, BB, 1, &
                         1, DESCB, ZERO, ZERO, DUMMY, DUMMY, ABSTOL, M, NZ, &
                         W, ORFAC, Z, 1, 1, DESCZ, WORK, LWORK, IWORK, &
                         LIWORK, IFAIL, ICLUSTR, GAP, INFO )

!           print *,"SCALAPACK info = ",info
           if( INFO /= 0 ) then
!              call error_handler( "solver, PDSYGVX: failed to find optimal &
!                                 &values for LWORK and LIWORK" )
              deallocate( WORK, IWORK, ICLUSTR, DESCA, DESCB, DESCZ, IFAIL, &
                          GAP, AA, BB, STAT=allocate_status )
              if( allocate_status /= 0 ) then
!                 call error_handler( "solver, PDSYGVX: deallocation &
!                                    &of arrays WORK, IWORK, ICLUSTR, &
!                                    &DESCA, DESCB, DESCZ, AA, BB, IFAIL &
!                                    &and/or GAP failed" )
              end if
           end if

           LWORK = WORK(1)
           LIWORK = IWORK(1)

           deallocate( WORK, IWORK, STAT=allocate_status )
           if( allocate_status /= 0 ) then
!              call error_handler( "solver, PDSYGVX: first allocation &
!                                 &of array(s) WORK and/or IWORK failed" )
              deallocate( ICLUSTR, DESCA, DESCB, DESCZ, IFAIL, &
                          GAP, AA, BB, STAT=allocate_status )
              if( allocate_status /= 0 ) then
!                 call error_handler( "solver, PDSYGVX: deallocation &
!                                    &of arrays WORK, IWORK, ICLUSTR, &
!                                    &DESCA, DESCB, DESCZ, AA, BB, IFAIL &
!                                    &and/or GAP failed" )
              end if
           end if


           if( LIWORK < 6 * MAX( N, NPROW * NPCOL + 1, 4 ) ) then
               LIWORK = 6 * MAX( N, NPROW * NPCOL + 1, 4 )
           end if

           allocate( WORK(LWORK), IWORK(LIWORK), STAT=allocate_status )

           if( allocate_status /= 0 ) then
!              call error_handler( "solver, PDSYGVX: second allocation &
!                                 &of array(s) WORK and/or IWORK failed" )
              deallocate( ICLUSTR, DESCA, DESCB, DESCZ, IFAIL, &
                          GAP, AA, BB, STAT=allocate_status )
              if( allocate_status /= 0 ) then
!                 call error_handler( "solver, PDSYGVX: deallocation &
!                                    &of arrays WORK, IWORK, ICLUSTR, &
!                                    &DESCA, DESCB, DESCZ, AA, BB, IFAIL &
!                                    &and/or GAP failed" )
              end if
           end if

!          Ask PDSYGVX to solve the eigenproblem.
!           print *,"DEBUG: Starting SCALAPACK calc"
           print *,"IBTYPE =",IBTYPE
           call PDSYGVX( IBTYPE, 'V', 'A', 'U', N, AA, 1, 1, DESCA, BB, 1, &
                         1, DESCB, ZERO, ZERO, DUMMY, DUMMY, ABSTOL, M, NZ, &
                         W, ORFAC, Z, 1, 1, DESCZ, WORK, LWORK, IWORK, &
                         LIWORK, IFAIL, ICLUSTR, GAP, INFO )
           print *,"SCALAPACK info = ",info

         if( MYPROW == 0 .and. MYPCOL == 0 ) then
            ! communication here
!            write(*,*) ' ScaLAPACK has been used'
!            write(*,*) ' Number of processes:', NPROW * NPCOL
!            write(*,*) ' Grid:', '(', NPROW, 'x', NPCOL, ')'
!            write(*,*) ' Dimension of the problem: N =', N
!        write(*,*) ' Number of times the eigensolver has been used: k =',&
!                       k
!            write(*,*) ' Blocking factor to distribute the &
!                        &matrices: NB =', NB
!            write(*,*) ' Number of eigenvalues found: M =', M
!            write(*,*) ' Number of eigenvectors computed: NZ =', NZ
!            write(*,*) ' Integer work space needed:', IWORK(1)
            if( INFO /= 0 ) then
               write(*,*) ' IFAIL =', IFAIL
               write(*,*) ' ICLUSTR =', ICLUSTR
            end if
!            write(*,*) ' Optimal amount of workspace &
!                         &(minimum value for LWORK):', LWORK, ';'
!            write(*,*) 'if orthogonality of the eigenvectors is affected &
!                        &then a larger value must be provided.'
!            write(*,*) ' INFO =', INFO

!           Print out the random number used to construct the matrix A.
!           TESTING ONLY.
!            write(*,*) ' Random number:', r
         end if

         if( info /= 0 ) then
!            write(*,*) 'solver: stopped due to error in &
!                                &eigensolver INFO =', INFO
!            call error_handler( "solver, PDSYGVX: error PDSYGVX" )
         end if

         deallocate( WORK, IWORK, ICLUSTR, DESCA, DESCB, DESCZ, &
                     IFAIL, GAP, AA, BB, STAT=allocate_status )

         if( allocate_status /= 0 ) then

!            call error_handler( "solver, PDSYGVX: final deallocation &
!                               &or array(s) WORK, IWORK, ICLUSTR, &
!                               &DESCA, DESCB, DESCZ, AA, BB, IFAIL and/or &
!                               &GAP failed" )

         end if

      else if( package == 'L_D' ) then

!         print *," Starting Lapack_D"

         allocate( WORK(LWORK), IWORK(LIWORK), STAT=allocate_status )
         if( allocate_status /= 0 ) then

!            call error_handler( "solver, DSYGVD: first allocation &
!                               &of array(s) WORK and/or IWORK failed" )
             stop
         end if

!          Get optimal values for LWORK and LIWORK.
           LWORK = -1
!           print *,'Lapack D: DSYGVD( ... )'
           call DSYGVD( IBTYPE, 'V', 'U', N, A, N, B, N, W, WORK, LWORK, &
                        IWORK, LIWORK, INFO )
           if( INFO /= 0 ) then
!              call error_handler( "solver, DSYGVD: failed to find optimal &
!                                 &values for LWORK and LIWORK" )
              deallocate( WORK, IWORK, STAT=allocate_status )
              if( allocate_status /= 0 ) then
!                 call error_handler( "solver, DSYGVD: deallocation &
!                                    &of array(s) WORK and/or IWORK failed" )

              end if

              stop
           end if

           LWORK = WORK(1)
           LIWORK = IWORK(1)
           deallocate( WORK, IWORK, STAT=allocate_status )

           if( allocate_status /= 0 ) then
!              call error_handler( "solver, DSYGVD: first allocation &
!                                 &of array(s) WORK and/or IWORK failed" )
               stop
           end if

           allocate( WORK(LWORK), IWORK(LIWORK), STAT=allocate_status )

           if( allocate_status /= 0 ) then

!              call error_handler( "solver, DSYGVD: second allocation &
!                                 &of array(s) WORK and/or IWORK failed" )
               stop

           end if
!          Ask DSYGVD to solve the eigenproblem.
!           print *,'Lapack D 2 : DSYGVD( ... )'
           call DSYGVD( IBTYPE, 'V', 'U', N, A, N, B, N, W, WORK, LWORK, &
                        IWORK, LIWORK, INFO )
           Z = A


!         if( info /= 0 ) then
!             print *,'solver: stopped due to error in &
!                                &eigensolver INFO =', INFO
!            call error_handler( "solver, DSYGVD: error DSYGVD" )
!         end if

!         print *,"Last deallocate, may be bug?"
         deallocate( WORK, IWORK, STAT=allocate_status )
 

         if( allocate_status /= 0 ) then
!            call error_handler( "solver, DSYGVD: deallocation of array(s) &
!                              &WORK and/or IWORK failed" )

         end if


      else

!         call error_handler( "solver: variable package is not in range" )

      end if

  end subroutine mirceas_solver
  !*************************************************************

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
    !  Subroutine called by: solve_para
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
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: help,i,j,sj,si
    integer(kind=i4_kind) :: n_irrep ! number of irreps
    !------------ Executable code ------------------------------------

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

  !--------------- End of module -------------------------------------
end module peigval


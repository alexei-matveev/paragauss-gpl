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
MODULE  result_module
  !-------------------------------------------------------------------
  !
  !  Purpose: 
  !  The routines in this module are used to create the readable
  !  textfiles which contain the relevant results from the 
  !  calculation.
  !
  !  The module contains routines for 
  !   - calculating S->S oscillator strengths from eigenvectors
  !   - assigning single particle transitions to excitation 
  !     energies
  !   - calculating Cauchy expansion coefficients from oscillator
  !     strengths
  !
  !  This module is called by "main_master()" only.
  !  The input data, i.e.
  !   - eigenvectors 
  !   - eigenvalues
  !   - IPA single particle MO transition
  !  is contained in the global data module.
  !
  !  Names of output files are defined within this module. 
  !  Subroutines/functions which do not begin with a module
  !  name - such as "get_main_SP_contributions()" - are 
  !  private to this module.  
  !
  !
  !  Module called by: main_master()
  !
  !
  !  References: Notes HH
  ! 
  !
  !  Author: HH
  !  Date:   2/1999
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: 
  ! Date:   ...
  ! Description: rebuilded by SB
  !
  !-------------------------------------------------------------------
#include <def.h>
  USE type_module          ! type specification parameters
  USE iounitadmin_module   ! routines for I/O
  USE output_module        ! contains output options

  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ public functions and subroutines ---------------------
  PUBLIC result_main
  PUBLIC result_osc_strengths
  PUBLIC results_exc
  PUBLIC results_SPc

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of constants and variables ---------------


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
CONTAINS


  !*************************************************************
  SUBROUTINE result_main(i_ir,N,eps,eta,MO_ALL,IRR)
    !  Purpose: 
    !------------ Modules used ------------------- ---------------
    USE global_module, ONLY: gl_NTO, gl_Nlow, gl_oscstr, gl_max_sp_trans
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),INTENT(IN) :: i_ir, N
    REAL   (KIND=r8_kind),INTENT(IN) :: eps(:), eta(:)
    INTEGER(KIND=i4_kind),INTENT(IN) :: MO_ALL(:,:), IRR(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)             :: status, N1
    INTEGER(KIND=i4_kind),ALLOCATABLE :: as_index(:,:,:)
    INTEGER(KIND=i4_kind),ALLOCATABLE :: as_irrep(:,:,:)
    INTEGER(KIND=i4_kind),ALLOCATABLE :: as_spin(:,:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: as_delta_eps(:,:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: as_val(:,:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: f_osc(:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: f_comp(:,:)
    INTEGER(KIND=i4_kind)             :: NSPTR
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    ! get the eigenvector
    N1 = min(gl_NLow,N)

    NSPTR = min(gl_max_sp_trans,N1)

    if(gl_NTO) NSPTR=N1

    ALLOCATE(&
         as_index    (N1,2,NSPTR),&
         as_irrep    (N1,2,NSPTR),&
         as_spin     (N1,  NSPTR),&
         as_delta_eps(N1,  NSPTR),&
         as_val      (N1,  NSPTR),&
         stat=status)
    ASSERT(status==0)

    CALL results_SPc(i_ir,N,eps,eta,MO_ALL,IRR,&
         as_index, as_irrep, as_spin, as_delta_eps, as_val)

    if (gl_oscstr) then 
       ALLOCATE(f_osc(N), f_comp(N,3),stat=status)
       ASSERT(status==0)
       f_osc     = 0.0_r8_kind
       f_comp    = 0.0_r8_kind

       CALL result_osc_strengths(i_ir,N,eps,eta,f_osc,f_comp)
       CALL results_exc(i_ir,N,as_index, as_irrep, as_spin, as_delta_eps, as_val, f_osc, f_comp)

       DEALLOCATE(f_osc, f_comp, stat=status)
       ASSERT(status==0)
    else
       CALL results_exc(i_ir,N,as_index, as_irrep, as_spin, as_delta_eps, as_val)
    end if

    DEALLOCATE(as_index, as_irrep, as_spin, as_delta_eps, as_val, stat=status)
    ASSERT(status==0)

  END SUBROUTINE result_main
  !*************************************************************

  !*************************************************************
  SUBROUTINE results_SPc(i_ir,N,eps,eta,MO,IRR,as_index, as_irrep, as_spin, as_delta_eps, as_val)
    ! Purpose: 
    ! Go through the absolute value of each eigenvector and find 
    ! the 1,...,NUM_INDICES largest components. Here
    !     NUM_INDICES = size(as_index,3)=size(as_val,2)
    ! For the j-th of these components of the i-th eigenvector 
    ! do save
    ! a) the original eigvec-value in absval(i,j)
    ! b) the a-index in as_index(i,1,j)
    ! c) the s-index in as_index(i,2,j)
    !------------ Modules used ------------------- ---------------
    USE global_module, ONLY: gl_calcall, gl_vec, gl_Nlow, gl_N_as, gl_N_spin
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind), INTENT( in) :: i_ir,N
    REAL(KIND=r8_kind),    INTENT( in) :: eps(:), eta(:)
    INTEGER(KIND=i4_kind), INTENT( in) :: MO(:,:),IRR(:,:)
    INTEGER(KIND=i4_kind), INTENT(out) :: as_index(:,:,:)
    INTEGER(KIND=i4_kind), INTENT(out) :: as_irrep(:,:,:)
    INTEGER(KIND=i4_kind), INTENT(out) :: as_spin(:,:)
    REAL(KIND=r8_kind),    INTENT(out) :: as_delta_eps(:,:)
    REAL(KIND=r8_kind),    INTENT(out) :: as_val(:,:)    
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)                :: status,n_as,i_as,as_max(1), j
    REAL   (KIND=r8_kind),ALLOCATABLE    :: eigvec(:,:), local_vec(:)
    LOGICAL              ,ALLOCATABLE    :: mask(:)
    INTEGER(KIND=i4_kind)                :: num_indices
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    IF(gl_calcall) THEN
       ALLOCATE(eigvec(N,N),STAT = status)
    else
       ALLOCATE(eigvec(N,min(gl_Nlow,N)),STAT=status) 
    END IF
    ASSERT(status==0)

    eigvec = gl_vec

    n_as = SIZE(eigvec,1)

    num_indices = SIZE(as_index,3)

    ALLOCATE(local_vec(n_as),&
         &   mask     (n_as),stat = status)
    ASSERT(status == 0)

    DO i_as=1,SIZE(eigvec,2) 

       local_vec = ABS(eigvec(:,i_as))
       mask = .true.
       do j = 1, num_indices
          as_max = MAXLOC(local_vec, mask)
          if ( (as_max(1) <= gl_N_as(i_ir,1)) .or. (gl_N_spin == 1) ) then
             as_spin(i_as,j) = 1_i4_kind
          else
             as_spin(i_as,j) = 2_i4_kind
          end if
          as_index    (i_as,1,j) = MO (as_max(1),1) 
          as_irrep    (i_as,1,j) = IRR(as_max(1),1)
          as_index    (i_as,2,j) = MO (as_max(1),2)
          as_irrep    (i_as,2,j) = IRR(as_max(1),2)
          as_delta_eps(i_as,j)   = eps(as_max(1))
          as_val      (i_as,j)   = eigvec(as_max(1),i_as)
          mask(as_max(1))        = .false.
       end do
    END DO

    DEALLOCATE(mask, eigvec, local_vec, stat=status)
    ASSERT(status==0)

  END SUBROUTINE results_SPc
  !*************************************************************

  !*************************************************************
  SUBROUTINE results_exc(i_ir,N,as_index, as_irrep, as_spin, as_delta_eps, as_val, f_osc, f_comp)
    !  Purpose: 
    !  Process the real results of the density matrix based
    !  methods which are
    !  (1) the squares of excitation energies (possibly
    !      real and imaginary part)
    !  (2) the corresponding oscillator strengths
    !  (3) the x,y,z contributions to the oscil. strengths
    !  and write the results into outputfiles.
    !  Note: 
    !  - the results are sorted in the ascending order of the
    !    array "real_eig" and of these only the first
    !    "num_eigen" values will be printed
    !  - if "print_all"=TRUE then ALL eigenvalues will be printed,
    !    NOT ONLY the first few "num_eigen"...
    !  - Works for both ss/st eigenvalues
    !  
    !------------ Modules used ------------------- ---------------
    USE global_module,     ONLY: gl_val, gl_Nlow,gl_Ir,gl_SS,gl_ST, gl_NTO, gl_max_sp_trans
    USE phys_param_module, ONLY: hartree2ev
    USE filename_module,   ONLY: outfile
    USE iounitadmin_module,ONLY: output_unit
    USE symmetry_data_module
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),INTENT(in)             :: i_ir, N
    INTEGER(KIND=i4_kind),INTENT(in)             :: as_index(:,:,:)
    INTEGER(KIND=i4_kind),INTENT(in)             :: as_irrep(:,:,:)
    INTEGER(KIND=i4_kind),INTENT(in)             :: as_spin (:,:)
    REAL(KIND=r8_kind),   INTENT(in)             :: as_delta_eps(:,:)
    REAL(KIND=r8_kind),   INTENT(in)             :: as_val      (:,:)
    REAL(KIND=r8_kind),   INTENT(in)   ,OPTIONAL :: f_osc(:)
    REAL(KIND=r8_kind),   INTENT(in)   ,OPTIONAL :: f_comp(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)                :: alloc_stat, io_unit, wirr, nto_unit
    INTEGER(KIND=i4_kind)                :: iter, n_iter, j, il
    REAL(KIND=r8_kind)                   :: omega_value
    INTEGER(KIND=i4_kind),ALLOCATABLE    :: index_list(:)
    REAL(KIND=r8_kind),ALLOCATABLE       :: eigval(:)
    CHARACTER(LEN=100)                   :: header_str,FMT4OUT
    LOGICAL                              :: exists

    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    ALLOCATE(eigval(min(gl_NLow,N)),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    eigval = gl_val
    wirr = gl_Ir

    ALLOCATE(index_list(SIZE(eigval,1)),stat=alloc_stat)
    ASSERT(alloc_stat == 0)

    DO iter=1, SIZE(index_list,1)
       index_list(iter) = iter
    END DO

    ! first sort the eigenvalues by the size of the real part
    ! prior to sorting:
    ! real_eig(i),f_osc(i) belong to the same index i
    ! after sorting:
    ! real_eig(i), f_osc(index_list(i)), ... belong to the same index i
    if (size(index_list,1) /= 1) CALL sort_by_first(eigval,index_list)

    io_unit=openget_iounit(file=TRIM(outfile(filename(i_ir))), &
                            status='unknown',form='formatted')

    if(gl_NTO) then
       inquire(file=TRIM(outfile("nto.tmp")), exist=exists)
       nto_unit=openget_iounit(file=TRIM(outfile("nto.tmp")), &
                               status='unknown',form='formatted',position='append')
       if(wirr==1) then
          if(gl_ST.and..not.gl_SS.and..not.exists) then
             write(nto_unit,*)2,ssym%n_irrep               !If only ST transitions are present
          else if(.not.gl_ST.and.gl_SS) then
             write(nto_unit,*)1,ssym%n_irrep           !If only SS transitions are present
          else if(gl_ST.and.gl_SS) then
             write(nto_unit,*)3,ssym%n_irrep
          end if
       end if
    end if


    if (gl_SS) then
       header_str = " ========================== EXCITATIONS AND OSCILLATOR STRENGTHS S->S ============================="
    elseif (gl_ST) then
       header_str = " ========================== EXCITATIONS AND OSCILLATOR STRENGTHS S->T ============================="
    else
       header_str = " =============================== EXCITATIONS AND OSCILLATOR STRENGTHS ============================="
    end if
    WRITE (output_unit,*) TRIM(header_str)

    header_str = "  SP   N  SYM    ENERGY    INIT    FINA   AMPLT      OSCSTR     OSCX      OSCY      OSCZ      DELTA"
    WRITE (io_unit,*)     TRIM(header_str)
    WRITE (output_unit,*) TRIM(header_str)

    header_str = "                  [eV]                                [au]      [au]      [au]      [au]      [eV]"
    WRITE (io_unit,*)     TRIM(header_str)
    WRITE (output_unit,*) TRIM(header_str)
    

    ! finally print the table
    n_iter = SIZE(eigval)

    if(gl_NTO) write(nto_unit,*)n_iter
    i_as_loop_: DO iter=1, n_iter
       if(gl_NTO) write(nto_unit,*)size(as_index,3)
       j_loop_: DO j = 1, size(as_index,3)
          if(gl_NTO) then
                if(j==1) then
                     il = index_list(iter)
                     if(eigval(iter)>=0.0_r8_kind) THEN
                        omega_value=SQRT(eigval(iter))
                     else
                        omega_value=eigval(iter)
                     end if
                     if(present(f_osc)) then
                         write(nto_unit,*)   omega_value*hartree2ev ,&
                                            & f_osc (il)
                     else 
                         write(nto_unit,*)   omega_value*hartree2ev
                     end if    
                end if    
                write(nto_unit,*)      as_spin(il,j),&
                                  & as_index(il,1,j),&
                                  & as_irrep(il,1,j),&
                                  & as_index(il,2,j),&
                                  & as_irrep(il,2,j),&
                                  & as_val  (il,  j)
                if(j.gt.min(size(as_index,3),gl_max_sp_trans)) cycle                  
          end if
          j_if_: if (j==1) then 

             IF(eigval(iter)>=0.0_r8_kind) THEN
                omega_value=SQRT(eigval(iter))
             else
                omega_value=eigval(iter)
             end IF

             il = index_list(iter)

             ! write frequency in atomic units and electron volts
             f_osc_if_: IF(PRESENT(f_osc)) THEN

                if (as_spin(il,1)==1) then
                    FMT4OUT = "(3X,'A',1X,I4,2X,A3,2X,F10.4,2(1X,I3,A3),2X,F6.2,2X,4(2X,F8.4),2X,F10.4)"
                else
                    FMT4OUT = "(3X,'B',1X,I4,2X,A3,2X,F10.4,2(1X,I3,A3),2X,F6.2,2X,4(2X,F8.4),2X,F10.4)"
                end if
                !! IN exc_SS_
                WRITE (io_unit,FMT4OUT) &
                        & iter,&
                        & ssym%name(wirr),&
                        & omega_value*hartree2ev,&
                        & as_index(il,1,1),&
                        & ssym%name(as_irrep(il,1,1)),&
                        & as_index(il,2,1),&
                        & ssym%name(as_irrep(il,2,1)),&
                        & as_val  (il,  1),&
                        & f_osc (il), &
                        & f_comp(il,1),&
                        & f_comp(il,2),&
                        & f_comp(il,3),&
                        & as_delta_eps(il,1)*hartree2ev
                !! IN output
                WRITE (output_unit,FMT4OUT) &
                        & iter,&
                        & ssym%name(wirr),&
                        & omega_value*hartree2ev,&
                        & as_index(il,1,1),&
                        & ssym%name(as_irrep(il,1,1)),&
                        & as_index(il,2,1),&
                        & ssym%name(as_irrep(il,2,1)),&
                        & as_val  (il,  1),&
                        & f_osc (il), &
                        & f_comp(il,1),&
                        & f_comp(il,2),&
                        & f_comp(il,3),&
                        & as_delta_eps(il,1)*hartree2ev
             else
                if (as_spin(index_list(iter),1)==1) then
                   FMT4OUT = "(3X,'A',1X,I4,2X,A3,2X,F10.4,2(1X,I3,A3),2X,F6.2,2X)"
                else
                   FMT4OUT = "(3X,'B',1X,I4,2X,A3,2X,F10.4,2(1X,I3,A3),2X,F6.2,2X)"
                end if
                !! IN exc_SS_
                WRITE (io_unit,FMT4OUT) &
                     & iter,&
                     & ssym%name(wirr),&
                     & omega_value*hartree2ev,&
                     & as_index(il,1,1),&
                     & ssym%name(as_irrep(il,1,1)),&
                     & as_index(il,2,1),&
                     & ssym%name(as_irrep(il,2,1)),&
                     & as_val  (il,  1)
                !! IN output
                WRITE (output_unit,FMT4OUT) &
                     & iter,&
                     & ssym%name(wirr),&
                     & omega_value*hartree2ev,&
                     & as_index(il,1,1),&
                     & ssym%name(as_irrep(il,1,1)),&
                     & as_index(il,2,1),&
                     & ssym%name(as_irrep(il,2,1)),&
                     & as_val  (il,  1)
             END IF f_osc_if_
          else !!$j_if_
             if (as_spin(il,j)==1) then
                FMT4OUT = "(3X,'A',1X,I4,2X,A3,12X,2(1X,I3,A3),2X,F6.2,44X,F10.4)"
             else
                FMT4OUT = "(3X,'B',1X,I4,2X,A3,12X,2(1X,I3,A3),2X,F6.2,44X,F10.4)"
             end if
             !! IN exc_SS_
             WRITE (io_unit,FMT4OUT) &
                        & iter,&
                        & ssym%name(wirr),&
                        & as_index(il,1,j),&
                        & ssym%name(as_irrep(il,1,j)),&
                        & as_index(il,2,j),&
                        & ssym%name(as_irrep(il,2,j)),&
                        & as_val  (il,  j),&
                        & as_delta_eps(il,j)*hartree2ev
             !! IN output
             WRITE (output_unit,FMT4OUT) &
                        & iter,&
                        & ssym%name(wirr),&
                        & as_index(il,1,j),&
                        & ssym%name(as_irrep(il,1,j)),&
                        & as_index(il,2,j),&
                        & ssym%name(as_irrep(il,2,j)),&
                        & as_val  (il,  j),&
                        & as_delta_eps(il,j)*hartree2ev
          end if j_if_
       end do j_loop_
       WRITE (io_unit,*) " === "
       WRITE (output_unit,*) " === "


    END DO i_as_loop_

    header_str = " =====END====="
    WRITE (output_unit,*) TRIM(header_str)

    ! close data file
    CALL returnclose_iounit(io_unit)
    if(gl_NTO) CALL returnclose_iounit(nto_unit)
    ! deallocate the index list
    DEALLOCATE(eigval, stat=alloc_stat)
    ASSERT(alloc_stat==0)
    DEALLOCATE(index_list, stat=alloc_stat)
    ASSERT(alloc_stat==0)

  END SUBROUTINE results_exc

  !*************************************************************
  SUBROUTINE result_osc_strengths(i_ir,N,eps,eta,f_osc,f_comp)
    ! Purpose: 
    ! Calculate oscillator strengths
    !
    !                 3
    !    f_i = (2/3) Sum  [transpose(d_m)*n^(1/2)*epsilon^(1/2)*x_i]^2    
    !                m=1
    ! where 
    !     x_i                 = ith eigenvector       
    !     d_m                 = transition dipole moment in direction m
    !     epsilon^(1/2)_as,bt = delta(as,bt)/(epsilon_s - epsilon_a)
    !     delta(i,j)          = Kronecker symbol with indices i,j
    !     n_as,bt             = delta(as,bt)*(n_occ(a) - n_occ(b))
    !------------ Modules used ------------------- ---------------
    USE global_module
    USE filename_module, ONLY: data_dir
    USE readwriteblocked_module, ONLY: readwriteblocked_startread, &
         readwriteblocked_read, readwriteblocked_stopread,  &
         readwriteblocked_tapehandle
    USE symmetry_data_module, ONLY: symmetry_data_n_partners
    USE resp_util_module, ONLY: resp_util_fname
    USE constants,ONLY: ONE, TWO, THREE, FOUR
    USE io
    USE debug
    USE comm_module
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind), INTENT(IN) :: N, i_ir
    REAL(KIND=r8_kind), INTENT(   IN) :: eps(:), eta(:)
    REAL(KIND=r8_kind), INTENT(inout) :: f_osc(:)
    REAL(KIND=r8_kind), INTENT(INOUT) :: f_comp(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)                :: alloc_stat, i_f
    INTEGER(KIND=i4_kind)                :: NS, NF, i_spin
    REAL   (KIND=r8_kind),ALLOCATABLE    :: h1_vec(:)
    REAL(KIND=r8_kind), ALLOCATABLE      :: trans_dip_mom(:,:),tdm_tmp(:,:)
    REAL   (KIND=r8_kind),ALLOCATABLE    :: eigvec(:,:),vec_i(:)
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    if ((.not. gl_SS) .and. gl_ST) return 

    ASSERT(comm_i_am_master())

    ALLOCATE(eigvec(N,min(N,gl_Nlow)),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    ALLOCATE(trans_dip_mom(N,3),vec_i(N),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    eigvec = gl_vec

    ALLOCATE(h1_vec(N), stat=alloc_stat)
    ASSERT(alloc_stat==0)
    
    h1_vec = SQRT(REAL(eps * eta,r8_kind)) !! eta

    NS = 1
    NF = gl_N_as(i_ir,1)
    do i_spin = 1, gl_N_spin
       ALLOCATE(tdm_tmp(gl_N_as(i_ir,i_spin),3),STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       call read_buffer( trim(data_dir)//'/'//resp_util_fname('dipol',i_ir,i_spin) &
            , tdm_tmp &
            )
       
       trans_dip_mom(NS:NF,1:3) = tdm_tmp
       NS = NF+1
       NF = N
       DEALLOCATE(tdm_tmp,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    end do

    DO i_f = 1, size(eigvec,2) ! loop over all osc. strengths

       vec_i = eigvec(:,i_f)
       
       vec_i = vec_i * h1_vec
       f_comp(i_f,1) = DOT_PRODUCT(trans_dip_mom(:,1),vec_i)       
       f_comp(i_f,2) = DOT_PRODUCT(trans_dip_mom(:,2),vec_i)
       f_comp(i_f,3) = DOT_PRODUCT(trans_dip_mom(:,3),vec_i)

       f_comp(i_f,1:3) = TWO * (f_comp(i_f,1:3)**2) / THREE
       f_osc(i_f)      = f_comp(i_f,1) + f_comp(i_f,2) + f_comp(i_f,3)

!!$       auxsum =  f_comp(i_f,1)**2 + f_comp(i_f,2)**2 + f_comp(i_f,3)**2
!!$       f_osc(i_f) = CC * auxsum / THREE

    END DO

    deallocate(eigvec,vec_i,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    DEALLOCATE(h1_vec,stat=alloc_stat)
    ASSERT(alloc_stat==0)

  END SUBROUTINE result_osc_strengths

  !*************************************************************
  SUBROUTINE sort_by_first(real_arr,int_arr)
    !  Purpose: 
    !  Sort the real array "real_arr" and the integer array 
    !  "int_arr" in the ascending order determined by "real_arr".
    !  
    !  Note: HEAPSORT is used
    ! (I hate to admit, it but this code is from NRecipies)
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    REAL(KIND=r8_kind),    INTENT(inout) :: real_arr(:)
    INTEGER(KIND=i4_kind), INTENT(inout) :: int_arr(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind), PARAMETER ::  max_iterations=10000_i4_kind
    INTEGER(KIND=i4_kind) :: N,I,IR,J,L,iter,iter_2, help_int
    REAL(KIND=r8_kind)    :: help_real

    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------


    ! first get the number of eigenvalues
    N = SIZE(real_arr)

    IF(N/=SIZE(int_arr,1)) CALL error_handler("out_tape_module:&
         &sort_by_first: real and integer parameter arrays are of &
         &different size !")

    ! calculate some auxiliary quantities
    L=N/2+1
    IR=N

    DO iter=1,max_iterations

       IF (L > 1) THEN
          L=L-1
          help_real=real_arr(L)
          help_int =int_arr (L)
       ELSE
          help_real=real_arr(IR)
          help_int =int_arr (IR)
          real_arr(IR)=real_arr(1)
          int_arr (IR)=int_arr (1)
          IR=IR-1
          IF (IR == 1) THEN
             real_arr(1)=help_real
             int_arr (1)=help_int
             EXIT
          END IF
       END IF

       I=L
       J=L+L
       DO iter_2=1,max_iterations
          IF(J>IR) EXIT
          IF (J < IR) THEN
             IF (real_arr(J) < real_arr(J+1) ) J=J+1
          END IF
          IF (help_real< real_arr(J)) THEN
             real_arr(I)=real_arr(J)
             int_arr (I)=int_arr (J)
             I=J
             J=J+J
          ELSE
             J=IR+1
          END IF
       END DO
       ! did we exit the loop correctly ?
       IF(iter_2>max_iterations+1) CALL error_handler(&
            "    out_tape_module: sort_by_first: max. number&
            & of iterations reached in inner loop")

       real_arr(I)=help_real
       int_arr (I)=help_int
    END DO

    ! did we exit the loop correctly ?
    IF(iter>max_iterations+1) CALL error_handler(&
         "    out_tape_module: sort_by_first: max. number&
         & of iterations reached in outer loop")

  END SUBROUTINE sort_by_first
  !*************************************************************

  !*************************************************************
  character(len=15) function filename(irr_c)
    !  Purpose: 
    !  Create a filename  
    !------------ Modules used ------------------- ---------------
    use global_module, only: gl_S_App, gl_S_APP_XC, gl_SS, gl_ST
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER (KIND=i4_kind),INTENT(in )   :: irr_c
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    character(len=4) :: irc_char
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------
    
    filename="exc_"

    if(gl_SS) then 
       filename=TRIM(filename)//"SS_"
    elseif (gl_ST) then
       filename=TRIM(filename)//"ST_"
    end if

    if (gl_S_APP) then
       if (gl_S_APP_XC) then
          filename=TRIM(filename)//"XC_"
       else
          filename=TRIM(filename)//"Q_"
       end if
    end if

    write (irc_char, '(i4)') irr_c
    irc_char = adjustl(irc_char)

    filename=TRIM(filename)//trim(irc_char)//".dat"

  END FUNCTION filename
  !*************************************************************

  !--------------- End of module -------------------------------------
END MODULE result_module

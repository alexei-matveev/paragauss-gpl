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
MODULE nto_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Calculate the NTO expansion and output it into files 
  !           
  !  References: Martin, Richard L., J. Chem. Phys., Volume 118, Issue 11, pp. 4775-4777 (2003) 
  !
  !  Author: Huix-Rotllant, Miquel
  !  Date: 08/2007
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
#include <def.h>
  USE type_module          ! type specification parameters
  USE iounitadmin_module   ! routines for I/O
  USE output_module        ! contains output options
  implicit none            ! by default, all the variables have to be specifically declared
  save                     ! save all variables defined in this module
  private                  ! by default, all names are private
  !------------ Declaration of types ---------------------------------
  !------------ Declaration of constants and variables ---------------
  !------------ Interface statements ---------------------------------
  !------------ public functions and subroutines ---------------------
  public nto_module_main                      !Subroutine called by tddft_diag
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !------------ Declaration of types ---------------------------------
   type ntoinfo
            real(r8_kind),allocatable      ::   subT(:,:),subU(:,:),subVT(:,:),subTDIAG(:)
            real(r8_kind),allocatable      ::   as_val(:)
            integer(i4_kind),allocatable   ::   occ(:),unocc(:)
            integer(i4_kind),allocatable   ::   pairs(:)
   end type ntoinfo
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
CONTAINS
  !*************************************************************
  SUBROUTINE nto_module_main()
    ! Purpouse:
    !     1) Read the information coming from the 
    !        response module contained in the temporal file nto.tmp
    !     2) Organize the results from the previous step 
    !        ==> Done in subroutine "reorder"
    !     3) Construct the density Matrix
    !     4) Call the subroutine dgesvd90 (lapack90 implementation 
    !        of Singular Value Decompositions for real matrices)
    !     5) Write the results in two ways:
    !         - Formatted  ==> readable for the user (basic output). 
    !                          There're so many files as irreps 
    !                          and spins we are calculating
    !                          The notation is: nto_SPIN_IRREP.dat where
    !
    !                                           SPIN = SS for Singlet-Singlet
    !                                                  ST for Singlet-Triplet
    !
    !                                           IRREP = Number of irrep which 
    !                                                   the results belong to
    !
    !         - Unformatted ==> to be used in the nto_plot_module to make 3D 
    !                           plots of NTOs. The filename is: nto.dat 
    !                           It conatins all the results of the previous files.
    !
    !------------ Modules used ------------------- ---------------
    USE filename_module     ,   ONLY: outfile ! Final output

    USE symmetry_data_module                             !Symmetry information

    USE matrix_eigenval     ,   ONLY: dgesvd90           !SVD routine

    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)                :: io_f, &             !Formatted file with nto output
                                            io_unf, &           !Unformatted file with nto information for plotting
                                            nto_info            !Temporal file with the information necessary for nto

    INTEGER(KIND=i4_kind)                :: n_spins  , &        !number of spin
                                            n_irreps , &        !number of irreps
                                            n_values , &        !number of excitation energies(per irrep/spin)
                                            n_nto    , &        !number of ntos that survive from the values
                                            n_idx    , &        !index associated with the spin
                                            n_pairs

    INTEGER(KIND=i4_kind)                :: i_irrep , &         !irrep loop variable
                                            i_spin  , &         !spin loop variable
                                            i_val               !excitation energy loop variable

    REAL(KIND=r8_kind)   , ALLOCATABLE   :: as_val(:),totU(:,:),totVT(:,:),TDIAG(:)               !Density vector

    INTEGER(KIND=i4_kind),ALLOCATABLE    :: pairs(:,:),tmp(:),fin_occ(:,:),fin_unocc(:,:)           

    type(ntoinfo),allocatable            :: nto_inf(:)

    LOGICAL                              :: exists      !Determine whether a file exist or not

    INTEGER(KIND=i4_kind)                :: alloc_stat               !Allocation status
    INTEGER(KIND=i4_kind)                :: counter                  !Counter used in the reorder process
    REAL(KIND=r8_kind)                   :: exc_energy,osc_value     !Excitation energy / oscillator strength
    REAL(KIND=r8_kind)                   :: max_value
    INTEGER(KIND=i4_kind)                :: a,s,ap,sp,idx,jdx,kdx,ldx,spin_ini !Loop variables

    INTEGER(KIND=i4_kind)                :: dimocc,dimunocc,mindim   !Occupied, unoccupied dimensions and minimum value of the both
    CHARACTER(100)                       :: label                    !Name of the file
    CHARACTER(3), ALLOCATABLE            :: symmetries(:)
    !------------ Declaration of external procedures ----------------
    !------------ Executable code ------------------------------------
!
    !Get the eigenvalues and the irrep we are working on
    allocate(symmetries(ssym%n_irrep))
    do idx=1,ssym%n_irrep
           symmetries(idx)=ssym%name(idx)
    end do
!
    !Clear the nto.dat file if there exists a previous one 
    inquire(file=TRIM(outfile("nto.dat")), exist=exists)
    if(exists) CALL system("rm "//TRIM(outfile("nto.dat")))
!
    io_unf=openget_iounit(file=TRIM(outfile("nto.dat")), &    !!Unformatted file with the information for making plots 
         status='new',form='formatted')                               
!
    nto_info=openget_iounit(file=TRIM(outfile("nto.tmp")), &  !!Unformatted file with information of the excitations 
         status='old',form='formatted') 
!
    read(nto_info,*)n_spins,n_irreps   !!Number of spins (SS,ST) of the calculation and number of irreps
    write(io_unf,*)n_spins,n_irreps    !!Written in the unformatted file
    if(n_spins.eq.3) then             !!If the file contains SS and ST, run from 1 to 2
          n_spins=2
          spin_ini=1
    else if(n_spins.eq.2) then        !!If the file contains only ST, stay in 2
          n_spins=2
          spin_ini=2
    else if(n_spins.eq.1) then        !!If the file contains only SS, stay in 1
          n_spins=1
          spin_ini=1
    end if
!    
    !!Reading the information of the excitations in the nto.tmp file, which is created by the result_module.f90 
    !!after a response calculation. It contains all the necessary information
    do  i_irrep=1,n_irreps                       !!Run over each irrep present in the file
        do  i_spin=spin_ini,n_spins              !!Run over SS and ST
             read(nto_info,*)n_values  !!Read how many values are there
         !    
           ! n_nto=get_num_ntos(nto_info,n_values) !!Calculate how many NTOs are going to be calculated
             write(io_unf,*)n_values                  !Write the number of ntos to an irrep contain  
                                                   !(written for all irreps and spins)
           ! if(n_nto.eq.0) exit                   !If no NTO survive, don't start the calculation
         !    
             !Create the file where the information of the NTO is going to be written and open it 
             !The name is 'nto_' + spin(SS/ST) + irrep + '.dat'
             write(label,'(i4)')i_irrep
             if(i_spin.eq.1) then  !!For SS
                inquire(file=TRIM(outfile('nto_SS_'//TRIM(label)//'.dat')), exist=exists)
                if(exists) call system('rm '//TRIM(outfile('nto_SS_'//TRIM(label)//'.dat')))
                io_f=openget_iounit(file=TRIM(outfile('nto_SS_'//TRIM(label)//'.dat')), &
                    status='new', form='formatted')
             else  !!For ST
                inquire(file=TRIM(outfile('nto_ST_'//TRIM(label)//'.dat')), exist=exists)
                if(exists) call system('rm '//TRIM(outfile('nto_ST_'//TRIM(label)//'.dat')))
                io_f=openget_iounit(file=TRIM(outfile('nto_ST_'//TRIM(label)//'.dat')), &
                    status='new', form='formatted')
             end if
         !   
             n_nto=0      
             do  i_val=1,n_values            !!Run over all the excitation energies
                   read(nto_info,*)n_idx     !!Determine how many values are there
                   allocate(pairs(n_idx,6),as_val(n_idx),stat=alloc_stat)
                   ASSERT(alloc_stat == 0)
                   read(nto_info,*)exc_energy,osc_value   !!Read excitation energy and oscillator strenght
                   do idx=1,n_idx
                        read(nto_info,*)pairs(idx,6),(pairs(idx,jdx),jdx=1,4),as_val(idx)
                   end do
!
                   pairs(:,5)=0_i4_kind               
                   n_pairs=0_i4_kind
                   do idx=1,n_idx
                       if(all(pairs(:,5).ne.0_i4_kind)) exit
                       if(pairs(idx,5).ne.0_i4_kind) cycle
                       n_pairs=n_pairs+1
                       do jdx=idx,n_idx
                           if(pairs(jdx,2).eq.pairs(idx,2).and. & 
                              pairs(jdx,4).eq.pairs(idx,4).and. &
                              pairs(jdx,6).eq.pairs(idx,6)      ) then
                              pairs(jdx,5)=n_pairs
                           end if
                       end do
                       pairs(idx,5)=n_pairs
                   end do
!
                   n_nto=n_nto+1
!
                   allocate(nto_inf(n_pairs),stat=alloc_stat)       
                   ASSERT(alloc_stat == 0)
!                   
                   do idx=1,n_pairs
                       counter=0
                       do jdx=1,n_idx
                           if(idx.ne.pairs(jdx,5)) cycle
                           counter=counter+1
                       end do
                       allocate(nto_inf(idx)%occ(counter),nto_inf(idx)%unocc(counter),stat=alloc_stat)
                       ASSERT(alloc_stat == 0)
                   !
                       counter=0
                       do jdx=1,n_idx
                            if(idx.ne.pairs(jdx,5)) cycle
                            counter=counter+1
                            nto_inf(idx)%occ(counter)=pairs(jdx,1)
                            nto_inf(idx)%unocc(counter)=pairs(jdx,3)
                       end do

                       do jdx=1,size(nto_inf(idx)%occ)
                            do kdx=jdx,size(nto_inf(idx)%occ)
                               if(nto_inf(idx)%occ(kdx).le.nto_inf(idx)%occ(jdx)) then
                                   counter = nto_inf(idx)%occ(kdx)
                                   nto_inf(idx)%occ(kdx)=nto_inf(idx)%occ(jdx)
                                   nto_inf(idx)%occ(jdx)=counter
                               end if
                            end do
                       end do

                       do jdx=1,size(nto_inf(idx)%unocc)
                            do kdx=jdx,size(nto_inf(idx)%unocc)
                                 if(nto_inf(idx)%unocc(kdx).le.nto_inf(idx)%unocc(jdx)) then
                                        counter = nto_inf(idx)%unocc(kdx)
                                        nto_inf(idx)%unocc(kdx)=nto_inf(idx)%unocc(jdx)
                                        nto_inf(idx)%unocc(jdx)=counter
                                 end if
                            end do
                       end do     

                       dimocc=1 
                       if(size(nto_inf(idx)%occ).gt.1) then
                           allocate(tmp(size(nto_inf(idx)%occ)))
                           tmp(1)=nto_inf(idx)%occ(1)
                           do jdx=2,size(nto_inf(idx)%occ)
                              if(nto_inf(idx)%occ(jdx).ne.nto_inf(idx)%occ(jdx-1)) then
                                  dimocc=dimocc+1
                                  tmp(dimocc)=nto_inf(idx)%occ(jdx)
                              end if
                           end do
                           deallocate(nto_inf(idx)%occ)
                           allocate(nto_inf(idx)%occ(dimocc))
                           nto_inf(idx)%occ(:)=tmp(1:dimocc)
                           deallocate(tmp)
                       end if

                       dimunocc=1
                       if(size(nto_inf(idx)%unocc).gt.1) then
                           allocate(tmp(size(nto_inf(idx)%unocc)))
                           tmp(1)=nto_inf(idx)%unocc(1) 
                           do jdx=2,size(nto_inf(idx)%unocc)
                              if(nto_inf(idx)%unocc(jdx).ne.nto_inf(idx)%unocc(jdx-1)) then
                                    dimunocc=dimunocc+1
                                    tmp(dimunocc)=nto_inf(idx)%unocc(jdx)
                              end if
                           end do
                           deallocate(nto_inf(idx)%unocc)
                           allocate(nto_inf(idx)%unocc(dimunocc))
                           nto_inf(idx)%unocc(:)=tmp(1:dimunocc)
                           deallocate(tmp)
                       end if

                       allocate( nto_inf(idx)%subT(dimocc,dimunocc)          , &
                                 nto_inf(idx)%subU(dimocc,dimocc)            , &
                                 nto_inf(idx)%subVT(dimunocc,dimunocc)       , &
                                 nto_inf(idx)%subTDIAG(min(dimocc,dimunocc)) , & 
                                 stat=alloc_stat                             )
                       nto_inf(idx)%subT=0_r8_kind
                       do a=1,dimocc
                         do s=1,dimunocc
                              nto_inf(idx)%subT(a,s)=exc_ampl(nto_inf(idx)%occ(a)    , &
                                                              nto_inf(idx)%unocc(s)  , &
                                                              idx                    , &
                                                              pairs                  , &
                                                              as_val                    )
                         end do
                       end do
                   !-----------------------------------------------------------------!
                   !SVD for real matrices. Implemented in modules/matrix_eigenval.f90!
                   !								     !
                   !Input: T(dimocc,dimunocc) (see Construc Transition Matrix above) !
                   !         dimocc = dimension of occupied states		     !
                   !         dimunocc = dimension of unoccupied states		     !
                   !								     !
                   !Output:    TDIAG: Natural Transition Matrix (vector of dimension !
                   !                  min(dimocc,dimunocc) )			     !
                   !               totU: Occupied Natural Transition Orbitals (matrix   !
                   !                  of dimension dimocc)			     !
                   !              totVT: Tranposed of the Unoccupied Natural Transition !
                   !                  Orbitals (matrix of dimension dimunocc)	     !
                   !								     !
                   !-----------------------------------------------------------------!
                !
                       CALL  dgesvd90(nto_inf(idx)%subT,nto_inf(idx)%subU,nto_inf(idx) &
                                      %subTDIAG,nto_inf(idx)%subVT)
                !       
                   end do
!
                   dimocc=0
                   dimunocc=0
                   do idx=1,size(nto_inf)
                       dimocc=dimocc+size(nto_inf(idx)%subU,1)
                       dimunocc=dimunocc+size(nto_inf(idx)%subVT,1)
                   end do 

                   mindim=0
                   do idx=1,size(nto_inf)
                      do jdx=1,size(nto_inf(idx)%subTDIAG)
                         if(nto_inf(idx)%subTDIAG(jdx).eq.0_r8_kind) cycle
                         mindim=mindim+1
                      end do
                   end do

                   allocate(totU(dimocc,mindim),totVT(mindim,dimunocc),TDIAG(mindim), &
                            fin_occ(dimocc,3),fin_unocc(dimunocc,3))
                   
                   counter=0
                   do idx=1,size(nto_inf)
                       do jdx=1,size(nto_inf(idx)%subTDIAG)
                           if(nto_inf(idx)%subTDIAG(jdx).eq.0_r8_kind) cycle
                           counter=counter+1
                           TDIAG(counter)=nto_inf(idx)%subTDIAG(jdx)
                       end do
                   end do

                   do idx=1,mindim
                       do jdx=idx,mindim
                           if(TDIAG(jdx).ge.TDIAG(idx)) then
                                   max_value=TDIAG(idx)
                                   TDIAG(idx)=TDIAG(jdx)
                                   TDIAG(jdx)=max_value
                           end if
                       end do
                   end do

                   totU=0_r8_kind
                   dimocc=0
                   do idx=1,size(nto_inf)
                       fin_occ(dimocc+1:dimocc+size(nto_inf(idx)%subU,1),1)=nto_inf(idx)%occ(:)
                       fin_occ(dimocc+1:dimocc+size(nto_inf(idx)%subU,1),2)=idx
                       dimocc=dimocc+size(nto_inf(idx)%subU,1)
                   end do

                   do idx=1,mindim
                      do jdx=1,n_pairs
                          do kdx=1,size(nto_inf(jdx)%subTDIAG)
                         !
                             if(nto_inf(jdx)%subTDIAG(kdx).eq.0_r8_kind) cycle
                             if(nto_inf(jdx)%subTDIAG(kdx).eq.TDIAG(idx)) then 
                                 do ldx = 1,dimocc
                                      if(fin_occ(ldx,2).ne.jdx) cycle
                                      counter= ldx+size(nto_inf(jdx)%subU,1)-1
                                      totU(ldx:counter,idx)=nto_inf(jdx)%subU(:,kdx)
                                      exit
                                 end do
                            end if    
                         !     
                          end do 
                      end do
                   end do

                   totVT=0_r8_kind
                   dimunocc=0
                   do idx=1,size(nto_inf)
                       fin_unocc(dimunocc+1:dimunocc+size(nto_inf(idx)%subVT,1),1)=nto_inf(idx)%unocc(:)
                       fin_unocc(dimunocc+1:dimunocc+size(nto_inf(idx)%subVT,1),2)=idx
                       dimunocc=dimunocc+size(nto_inf(idx)%subVT,1)
                   end do

                   do idx=1,min(dimocc,dimunocc)
                      do jdx=1,size(nto_inf) 
                          do kdx=1,size(nto_inf(jdx)%subTDIAG)
                         !
                             if(nto_inf(jdx)%subTDIAG(kdx).eq.0_r8_kind) cycle
                             if(nto_inf(jdx)%subTDIAG(kdx).eq.TDIAG(idx)) then 
                                 do ldx = 1,dimunocc
                                      if(fin_unocc(ldx,2).ne.jdx) cycle
                                         counter= ldx+size(nto_inf(jdx)%subVT,2)-1
                                         totVT(idx,ldx:counter)=nto_inf(jdx)%subVT(kdx,:)
                                         exit
                                 end do
                             end if    
                         !     
                          end do
                      end do
                   end do

                   do idx=1,size(fin_occ,1)
                       do jdx=1,size(pairs,1)
                             if(fin_occ(idx,2).ne.pairs(jdx,5)) cycle 
                             fin_occ(idx,2)=pairs(jdx,2)
                             fin_occ(idx,3)=pairs(jdx,6)
                             exit
                       end do 
                   end do

                   do idx=1,size(fin_unocc,1)
                       do jdx=1,size(pairs,1)
                             if(fin_unocc(idx,2).ne.pairs(jdx,5)) cycle
                             fin_unocc(idx,2)=pairs(jdx,4)
                             fin_unocc(idx,3)=pairs(jdx,6)
                             exit
                       end do
                   end do

                   !-------------------------------
                   !End SVD calculation           !
                   !-------------------------------
                !
                   !--------------------------------------------------------------------------------------!
                   !Write the information in a file							  !
                   !											  !
                   !io_unf: Unformatted file, for doing a NTO plot					  !
                   !io_f:   Formatted file,  output for the user					  !
                   !											  !
                   ! Output Structure in formatted files(for each irrep):				  !
                   !											  !
                   !      1) Number of NTO								  !
                   !      2) Transition Energy in eV							  !
                   !      3) Transition Matrix in vector form						  !
                   !         (the original one, obtained after a TDDFT calculation)			  !
                   !      4) Natural Transition Matrix (TDIAG)						  !
                   !      5) Occupied NTOs (totU)								  !
                   !      6) Unoccupied NTOs (totVT). The transpose of this matrix is printed out		  !
                   !											  !
                   !      3) in Matrix form is the input of the SVD subroutine 				  !
                   !         (dgesvd90, modules/matrix_eigenval.f90)					  !
                   !      4),5),6) are the outputs of the SVD subroutine 				  !
                   !               (dgesvd90, modules/matrix_eigenval.f90)				  !
                   !											  !
                   !											  !
                   ! Output Structure in unformatted files :						  !
                   !											  !
                   !      1) irrep									  !
                   !      2) number of NTO, occupied dimension,unoccupied dimension			  !
                   !      4) occupied   NTO expansion labels (number MO, number symm), 			  !
                   !         expansion coefficients (dimension = minimum value of 2) )			  !
                   !      5) unoccupied NTO expansion labels (number MO, number symm), 			  !
                   !         expansion coefficients (dimension = minimum value of 2) )			  !
                   ! 											  !
                   !   /****  NOTE  ********************************************************************\ !
                   !   *										* !
                   !   *    If SS and ST  the results are alternated                     	        * !
                   !   *        Irrep 1 , SS								* !
                   !   *        Irrep 1 , ST								* !
                   !   *        Irrep 2 , SS								* !
                   !   *        Irrep 2 , ST                                                      	* !
                   !   \********************************************************************************/ !
                   ! 											  !
                   !--------------------------------------------------------------------------------------!

                   do idx=1,size(nto_inf)
                       deallocate(nto_inf(idx)%occ      , &
                                  nto_inf(idx)%unocc    , &
                                  nto_inf(idx)%subT     , &
                                  nto_inf(idx)%subU     , & 
                                  nto_inf(idx)%subVT    , &
                                  nto_inf(idx)%subTDIAG , &
                                  stat=alloc_stat) 
                       ASSERT(alloc_stat == 0) 
                   end do
                   deallocate(nto_inf,stat = alloc_stat)
                   ASSERT(alloc_stat == 0)
                  
                   mindim=0
                   do idx=1,size(TDIAG)
                       if(TDIAG(idx).lt.0.004) exit
                       mindim=mindim+1
                   end do 
                   if(mindim.gt.10) mindim = 10
                   
                   WRITE(io_f,*)
                   WRITE(io_f,*)'-----------------------------------------------------------'
                   write(io_unf,*)dimocc,dimunocc,mindim  !! Write 2) into unformatted file 
                   
                   !!Number of NTO (1)
                   WRITE(io_f,*) n_nto,' NTO'                 

                   !!Transition Energy(2)
                   WRITE(io_f,'(A,F8.4)')'Transition Energy [eV]:  ', exc_energy 
                   WRITE(io_f,*)
                   WRITE(io_f,*) 

                 
                   !!Original Transition Matrix(3)
              !    WRITE(io_f,*)'Transition Matrix (vector form)'
              !    WRITE(label,*)'(A,2(1X,I3,A3),f9.6)'
              !    DO a=1,size(pairs,1)
              !       IF(pairs(a,6).eq.1) then
              !           WRITE(io_f,label)'A',pairs(a,1),TRIM(symmetries(pairs(a,2))),&
              !                             & pairs(a,3),TRIM(symmetries(pairs(a,4))), &
              !                             & as_val(a)
              !       ELSE
              !           WRITE(io_f,label)'B',pairs(a,1),TRIM(symmetries(pairs(a,2))),&
              !                             & pairs(a,3),TRIM(symmetries(pairs(a,4))), &
              !                             & as_val(a)
              !       END IF
              !    END DO
                   WRITE(io_f,*)
                   WRITE(io_f,*)      
                  
                   !!Natural Transition Matrix (4)
                   WRITE(io_f,*)'Natural Transition Matrix (vector form)'
                   WRITE(label,*)'(8X',mindim,'(f9.6))'
                   WRITE(io_f,label)(TDIAG(a),a=1,mindim)
                   write(io_f,*)sum(TDIAG**2)
                   WRITE(io_f,*)
                   WRITE(io_f,*) 

                   !!Occupied Natural Transition Orbitals (5)
                   WRITE(io_f,*)'Occupied Natural Transition Orbitals'
                   WRITE(label,*)'(A,I3,A3,',mindim,'(F9.6))'
                   DO a=1,dimocc
                       if(all(abs(totU(a,:mindim)).lt.5d-6)) cycle 
                       if(fin_occ(a,3).eq.1) then
                        WRITE(io_f,label)'A',fin_occ(a,1),TRIM(symmetries(fin_occ(a,2))), &
                                       (totU(a,ap),ap=1,mindim)
                       else                
                        WRITE(io_f,label)'B',fin_occ(a,1),TRIM(symmetries(fin_occ(a,2))), &
                                       (totU(a,ap),ap=1,mindim)
                       end if
                        write(io_unf,*)fin_occ(a,3),fin_occ(a,1),fin_occ(a,2),(totU(a,ap),ap=1,mindim)
                   END DO
                   WRITE(io_f,*)
                   WRITE(io_f,*)

                   !!Unoccupied Natural Transition Orbitals (6)
                   WRITE(io_f,*)'Unoccupied Natural Transition Orbitals'
                   WRITE(label,*)'(A,I3,A3,',mindim,'(F9.6))'
                   DO s=1,dimunocc
                        if(all(abs(totVT(:mindim,s)).lt.5d-6)) cycle 
                        if(fin_unocc(s,3).eq.1) then
                        WRITE(io_f,label)'A',fin_unocc(s,1),TRIM(symmetries(fin_unocc(s,2))), &
                                         (totVT(sp,s),sp=1,mindim)
                        else
                        WRITE(io_f,label)'B',fin_unocc(s,1),TRIM(symmetries(fin_unocc(s,2))), &
                                         (totVT(sp,s),sp=1,mindim)
                        end if
                        WRITE(io_unf,*)fin_unocc(s,1),fin_unocc(s,2), &
                                       (totVT(sp,s),sp=1,mindim)
                   END DO
                   WRITE(io_f,*)
                   WRITE(io_f,*)
                   !-------------------------------!
                   !End write information in a file!
                   !-------------------------------!

                   deallocate(pairs,as_val,stat=alloc_stat)
                   ASSERT(alloc_stat == 0)
                   deallocate(totU,totVT,fin_occ,fin_unocc,TDIAG,stat=alloc_stat)
                   ASSERT(alloc_stat == 0)
             END DO
             !Close formated output file
             CALL returnclose_iounit(io_f)
        END DO 
    END DO 

    !Close output unformatted file
    CALL returnclose_iounit(io_unf)
    !Close temporal file
    CALL returnclose_iounit(nto_info)

  END SUBROUTINE nto_module_main 
  !********************************************************************************************
  REAL(kind=r8_kind) FUNCTION exc_ampl(occ,unocc,irrep,pairs,as_val)
    implicit none
    INTEGER(KIND=i4_kind),INTENT(in)  :: occ,unocc,irrep
    INTEGER(KIND=i4_kind),INTENT(in)  :: pairs(:,:)
    REAL(KIND=r8_kind),INTENT(in)     :: as_val(:)
    INTEGER(KIND=i4_kind)             :: indx
      exc_ampl=0_r8_kind
      do indx=1,size(pairs,1)
             if(pairs(indx,5).ne.irrep) cycle
             if(pairs(indx,1).eq.occ.and.pairs(indx,3).eq.unocc) then
                  exc_ampl=exc_ampl+as_val(indx)
             end if
      end do
  end function exc_ampl

  !********************************************************************************************
  INTEGER(kind=i4_kind) FUNCTION get_num_ntos(nto_info,n_values)
    !  Purpose:  Get the numer of NTOs that will survive after some conditions are met                                                     
    !------------ Modules used ------------------- ---------------
    USE global_module, only: gl_oscstr
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),INTENT(in)  :: nto_info      !File unit where the information is
    INTEGER(KIND=i4_kind),INTENT(in)  :: n_values      !Number of values that are present in the file
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind)                        :: idx,jdx    !Variables for the loops
    INTEGER(KIND=i4_kind)                        :: go_back    !Number of positions to rewind the file
    INTEGER(KIND=i4_kind)                        :: irr1,idx1,irr2,idx2,n_idx
    INTEGER(KIND=i4_kind)                        :: num_nto
    REAL(KIND=r8_kind)                           :: exc_energy   , &
                                                    osc_value 
    REAL(KIND=r8_kind), ALLOCATABLE              :: exc_amplitude(:)
    !------------ Declaration of subroutines used ----------------
    !------------ Executable code ------------------------------------
    num_nto=0
    go_back=0
    do idx=1,n_values
!  
       read(nto_info,*)n_idx
       if(allocated(exc_amplitude)) deallocate(exc_amplitude)
       allocate(exc_amplitude(n_idx))
       if(gl_oscstr) then
              read(nto_info,*)exc_energy,osc_value
       else
              read(nto_info,*)exc_energy
       end if
!
       do jdx=1,n_idx
               read(nto_info,*)irr1,idx1,irr2,idx2,exc_amplitude(jdx)
       end do
!
       go_back=go_back+n_idx+2
!
       if(any(exc_amplitude.eq.0_r8_kind)) cycle
!
       num_nto=num_nto+1
    end do
!
    if(num_nto.ne.0) then
        do idx=1,go_back
           backspace nto_info
        end do
    end if
!
    get_num_ntos=num_nto
!
  END FUNCTION get_num_ntos
!************************************************************************
END MODULE nto_module

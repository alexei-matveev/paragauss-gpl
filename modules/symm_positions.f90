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
module symm_positions
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

# include "def.h"
  use type_module, only:&
       & IK => i4_kind,&
       & RK => r8_kind ! type specification parameters
  use group_module ,only: symm_transformation_int
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  integer(IK),dimension(:),allocatable,public:: pcr_no,psb_ind
  type(symm_transformation_int),pointer,public      :: point_trafos(:)
  ! transformations of the symmetry equivalent atoms among each other
  ! by all symmetry operations
  ! point_trafos(N_unique_atoms)
  type(symm_transformation_int),pointer,public      :: point_trafos_timps(:)
  ! transformations of the symmetry equivalent point charges
  ! among each other by all symmetry operations
  ! point_trafos(n_timps)
  type(symm_transformation_int),pointer,public      :: point_trafos_pc1(:)

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public init,done,symm_equivalents_gen
  public symm_mapping_gen
  public :: symm_mapping_shells

  public :: symm_nabor

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init()
    use unique_atom_module, only: n_unique_atoms
    use pointcharge_module, only: n_timps, n_moving_unique_timps
    use pointcharge_module, only: pointcharge_n
    implicit none
    ! *** end of interface ***
    integer(IK) :: memstat

    DPRINT 'sp/init: entered'

    ! allocate transformation matrix of the equal atoms
    allocate(point_trafos(N_unique_atoms),STAT=memstat)
    if(memstat/=0) call error_handler("sp/init: alloc 1 failed")

!!$    if((n_timps.ne.0.and.N_moving_unique_timps.ne.0) .or. pointcharge_N.ne.0)then
!!$       allocate(point_trafos_timps(n_timps+pointcharge_N),STAT=memstat)
    if(n_timps.ne.0.and.N_moving_unique_timps.ne.0)then
       allocate(point_trafos_timps(n_timps),STAT=memstat)
       if(memstat/=0) call error_handler("sp/init: alloc 2 failed")
    endif

    if(pointcharge_N > 0) then
       allocate(point_trafos_pc1(pointcharge_N),STAT=memstat)
       ASSERT(memstat==0)
    end if

  end subroutine init

  subroutine done()
    implicit none
    ! *** end of interface ***
    integer(IK) :: memstat,i

    DPRINT 'sp/done: entered'

    do i=1,size(point_trafos)
       deallocate(point_trafos(i)%matrix,STAT=memstat)
       if(memstat/=0) call error_handler("sp/done: dealloc 1 failed")
    enddo
    deallocate(point_trafos,STAT=memstat)
    if(memstat/=0) call error_handler("sp/done: dealloc 2 failed")

    if(associated(point_trafos_timps))then
       do i=1,size(point_trafos_timps)
          deallocate(point_trafos_timps(i)%matrix,STAT=memstat)
          if(memstat/=0) call error_handler("sp/done: dealloc 3 failed")
       enddo
       deallocate(point_trafos_timps,STAT=memstat)
       if(memstat/=0) call error_handler("sp/done: dealloc 4 failed")
    endif

    if(associated(point_trafos_pc1))then
       do i=1,size(point_trafos_pc1)
          deallocate(point_trafos_pc1(i)%matrix,STAT=memstat)
          ASSERT(memstat==0)
       enddo
       deallocate(point_trafos_pc1,STAT=memstat)
       ASSERT(memstat==0)
    endif

  end subroutine

  !*************************************************************
  subroutine symm_equivalents_gen
    !  Purpose: generates all equal atoms from one
    !           equal atom given
    !           generates symmetry equivalent distances
    !           treats EPE related PC arrays
    !           based on data in gxepe_array leaves in
    !           EPE related arrays only PC centers of
    !           environmen
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module      , only: sub_group                                    &
                                , group_coset                                  &
                                , group_coset_decomp                           &
                                , group_num_el                                 &
                                , ylm_trafos
    use efm_module
    use unique_atom_module
    use pointcharge_module
    use datatype
    use operations_module, only : operations_qm_epe
    implicit none
    !------------ Declaration of local constants  ----------------
    real(RK),parameter :: small = 1.e-10_rk
    ! very small value
    !------------ Declaration of local variables -----------------
    integer(IK)                :: i,j,k, n,l,alloc_stat,error
    ! counters
    integer(IK)                :: n_equal,n_equiv,n_centers
    ! number of equal atoms
    real(RK)                   :: position(3),position2(3),position3(3)
    ! coordinates of unique atom in the order (x,z,y)
    type(sub_group),pointer              :: local_groups(:)
    ! local symmetry group of the unique atoms
    ! local_groups(N_unique_atoms)
    type(group_coset),pointer            :: cosets(:)
    ! cosets of the local symmetry groups
    integer(IK),allocatable    :: sum_matrix(:,:)
    ! sum_matrix(n_equal,n_equal)
    integer(IK)                :: sum_int
    ! auxiliary variable
    type(sub_group)                      :: local_groups_pc
    ! local symmetry group of the unique point charges
    type(group_coset)                    :: cosets_pc
    ! coset of the local symmetry group of the unique point charges

    type(pointcharge_type), pointer      :: pc
    type(unique_atom_type), pointer      :: ut ! pointer to point on timps
    logical:: a_switch
    !------------ Declaration of subroutines used ----------------
    external error_handler
    !------------ Executable code --------------------------------

    DPRINT 'sp/symm_equivalents_gen: entered'

    ! allocate local arrays
    allocate(local_groups(N_unique_atoms))
    allocate(cosets(N_unique_atoms))

    !
    ! determine local symmetry group of uniques
    !


    ! loop over all uniques
    n_centers=0
    do i=1,N_unique_atoms
       ! reorder coordinates of unique atom as
       ! (x,y,z) --> (x,z,y) in order to comply with the
       ! convention for angular momentum l=1
       position(1) = unique_atoms(i)%position_first_ea(1)
       position(2) = unique_atoms(i)%position_first_ea(3)
       position(3) = unique_atoms(i)%position_first_ea(2)

       !
       ! determine local symmetry groups
       !

       ! now apply all symmetry operations to the position of the
       ! unique atom
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
          end if
       end do

       ! allocate group elements
       local_groups(i)%num_el = n
       allocate(local_groups(i)%elements(n))

       ! fill up group elements
       n = 0
       position3 = 0.0_rk
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups(i)%elements(n) = j
             position3 = position3 + position2
          end if
       end do
       ! adjust position of atom to lie exactly on a symmetry plane/axis
       position = position3/local_groups(i)%num_el
       unique_atoms(i)%position_first_ea(1) = position(1)
       unique_atoms(i)%position_first_ea(2) = position(3)
       unique_atoms(i)%position_first_ea(3) = position(2)

       !
       ! now determine symmetry equivalent atoms
       !

       call group_coset_decomp(n_equal,local_groups(i),&
            cosets(i),point_trafos(i)%matrix)

       ! allocate positions of equal atoms
       allocate(unique_atoms(i)%position(3,n_equal),stat=alloc_stat)
       if (alloc_stat .ne. 0 ) call error_handler( &
            "symm_equivalents_gen: allocation of unique_atom(unique)%%position failed")
       !n_equal = unique_atoms(i)%N_equal_atoms
       n_centers=n_centers+n_equal
       ! determine positions of equal atoms
       if (output_unit > 0) then
          write(output_unit,*) " Equal atoms of Type ",i
       endif
       do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets(i)%elements(1,j)),position)
          unique_atoms(i)%position(1,j) = position2(1)
          unique_atoms(i)%position(2,j) = position2(3)
          unique_atoms(i)%position(3,j) = position2(2)
          if (output_unit > 0) then
             write(output_unit,*) "x: ",position2(1)," y: ",position2(3)," z: ",position2(2)
          endif
       end do

       if (n_equal .ne. unique_atoms(i)%N_equal_atoms) call error_handler( &
            "symm_equivalents_gen: calculate number of equal atoms and input number differ" )

    end do

    ! output section -----------------------------------------------------
    if (output_unit > 0) then
       call output_geometry (output_unit, unique_atoms)
       call symm_nabor (output_unit, unique_atoms, 0)
       write (output_unit, *) " "
    endif
    ! --------------------------------------------------------------------

    !
    ! determine symmetry equivalent distances
    !

    ! allocate symmequiv table
    a_switch=.not.allocated(unique_atom_symequiv)
    ASSERT(a_switch)
    allocate(unique_atom_symequiv(N_unique_atoms,N_unique_atoms),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! loop over all unique atoms
    do i=1,N_unique_atoms
       ! loop over all unique atoms
       do j=1,N_unique_atoms
          n_equal = unique_atoms(j)%N_equal_atoms

          ! allocate work array
          allocate(sum_matrix(n_equal,n_equal))
          sum_matrix = 0

          ! add all point transformation matrices of the subgroup of i
          do k=1,local_groups(i)%num_el
             sum_matrix = sum_matrix + point_trafos(j)%matrix(:,:,&
                  local_groups(i)%elements(k))
          end do

          n_equiv  = 0
          do k=1,n_equal
             sum_int  = 0
             do l=1,n_equal
                if (sum_matrix(l,k).ne.0) then
                   if(k.ne.l) then
                      sum_matrix(:,l) = 0
                   end if
                   sum_int = sum_int + 1
                end if
             end do
             sum_matrix(1,k) = sum_int
             if (sum_int.gt.0) then
                n_equiv = n_equiv + 1
             end if
          end do

          ! allocate symmequiv distances
          unique_atom_symequiv(i,j)%n_equiv_dist = n_equiv

          allocate(unique_atom_symequiv(i,j)%index_eq_atom2(n_equiv),STAT=alloc_stat)
          ASSERT(alloc_stat.eq.0)

          a_switch=.not.allocated(unique_atom_symequiv(i,j)%weight)
          ASSERT(a_switch)
          allocate(unique_atom_symequiv(i,j)%weight(n_equiv), STAT=alloc_stat) ! (1)
          ASSERT(alloc_stat.eq.0)

          ! now fill in symmequiv distances
          if (output_unit > 0) then
             write(output_unit,*) "symmetry equivalent distances from type ",i," to type ",j
             write(output_unit,*) "there are ",n_equiv," different distances"
          endif
          sum_int = 0
          do k=1,n_equal
             if (sum_matrix(1,k).gt.0) then
                sum_int = sum_int + 1
                unique_atom_symequiv(i,j)%index_eq_atom2(sum_int) = k
                unique_atom_symequiv(i,j)%weight(sum_int) = unique_atoms(i)%N_equal_atoms &
                     * sum_matrix(1,k)
                if (output_unit > 0) then
                   write(output_unit,*) "# ",sum_int," representant: ",k," weight: ", &
                        unique_atom_symequiv(i,j)%weight(sum_int)
                endif
             end if
          end do
          ! deallocate temporary array
          deallocate(sum_matrix)
       end do
    end do

    ! deallocate local arrays
    do i=1,ubound(local_groups,1)
       deallocate(local_groups(i)%elements,stat=error)
       ASSERT(error.eq.0)
    enddo
    deallocate(local_groups,stat=error)
    if (error/=0) call error_handler &
         ("symm_equivalents_gen: deallocation (2) failed")
    do i=1,ubound(cosets,1)
       deallocate(cosets(i)%elements,stat=error)
       if (error/=0) call error_handler &
            ("symm_equivalents_gen: deallocation (3) failed")
    enddo
    deallocate(cosets,stat=error)
    if (error/=0) call error_handler &
         ("symm_equivalents_gen: deallocation (4) failed")

    !
    ! determine local symmetry group of unique point charges
    !

    ! loop over all unique charges
    do i=1,pointcharge_N+n_timps
       if(i<=n_timps) then
          ut => unique_timps(i)
          position(1) = ut%position_first_ea(1)
          position(2) = ut%position_first_ea(3)
          position(3) = ut%position_first_ea(2)

       else

          pc => pointcharge_array(i-n_timps)

          ! reorder coordinates of unique atom as
          ! (x,y,z) --> (x,z,y) in order to comply with the
          ! convention for angular momentum l=1
          position(1) = pc%position_first_ec(1)
          position(2) = pc%position_first_ec(3)
          position(3) = pc%position_first_ec(2)
       end if
       !
       ! determine local symmetry groups
       !

       ! now apply all symmetry operations to the position of the
       ! unique atom
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
          end if
       end do

       ! allocate group elements
       local_groups_pc%num_el = n
       allocate(local_groups_pc%elements(n))

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups_pc%elements(n) = j
          end if
       end do

       !
       ! now determine symmetry equivalent atoms
       !

       if(i<=n_timps) then
          call group_coset_decomp(n_equal,local_groups_pc,&
               cosets_pc,point_trafos_timps(i)%matrix)
       else
          call group_coset_decomp(n_equal,local_groups_pc,&
               cosets_pc,point_trafos_pc1(i-n_timps)%matrix)
       end if


       ! allocate positions of equal charges
       if(i<=n_timps) then
          allocate(ut%position(3,ut%n_equal_atoms),stat=alloc_stat)
          if (alloc_stat .ne. 0 ) call error_handler( &
               "symm_equivalents_gen: allocation of ut%position failed")
       else
          allocate(pc%position(3,pc%N_equal_charges),stat=alloc_stat)
          if (alloc_stat .ne. 0 ) call error_handler( &
               "symm_equivalents_gen: allocation of pc%position failed")
       end if

       ! determine positions of equal atoms
       if(i<=n_timps) then
          if (output_unit > 0) then
             write(output_unit,*) " Equal timp ",i
          endif
       else
          if (output_unit > 0) then
             write(output_unit,*) " Equal charges of Type ",i-n_timps
          endif
       end if
       do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
          if(i<=n_timps) then
             ut%position(1,j) = position2(1)
             ut%position(2,j) = position2(3)
             ut%position(3,j) = position2(2)
          else
             pc%position(1,j) = position2(1)
             pc%position(2,j) = position2(3)
             pc%position(3,j) = position2(2)
          end if
          !write(output_unit,*) "x: ",position2(1)," y: ",position2(3)," z: ",position2(2)
       end do

!!$       deallocate(local_groups_pc%elements, &
!!$            point_trafos_pc%matrix,stat=alloc_stat) ! to be used in symm_symmadapt
       deallocate(local_groups_pc%elements, stat=alloc_stat)

       if (alloc_stat .ne. 0 ) call error_handler( &
            "symm_equivalents_gen: deallocation of point charge helpers failed")

       if(i<=n_timps) then
          if (n_equal .ne. ut%n_equal_atoms) call error_handler( &
               "symm_equivalents_gen: calculate number of equal timps and input number differ" )
       else
          if (n_equal .ne. pc%N_equal_charges) call error_handler( &
               "symm_equivalents_gen: calculate number of equal charges and input number differ" )
       end if
    end do

    if (n_timps.gt. 0) call output_geometry_timp()
    if (pointcharge_N.gt. 0) call output_geometry_pc()

    if(operations_qm_epe) then
#ifdef WITH_EPE
       call symm_epe()
       call symm_ewpc_gen()
#else
       ABORT('recompile -DWITH_EPE')
#endif
    endif

 contains

#if SYMM_EPE
  subroutine symm_epe
! treats EPE centers
  use ewaldpc_module,ex_pcr_dummy=>ex_pcr
  use filename_module
  type(unique_atom_type), pointer :: ua
        logical:: ex_pcr
        real (RK):: e_nuc_pcr,dist
        integer(IK):: na,nb,eq_a,pcr_unit,n_gxat,kl
  integer(IK)::  it_pcr,ii,i_ua,i_eq
  real(RK), allocatable, dimension(:,:) ::pcs_temp
  real(RK), allocatable, dimension(:,:) ::pcc_temp
  real(RK),parameter:: small_pcr=5.d-1

  inquire (file= trim(input_dir)//'/epe.pcr',exist=ex_pcr)
  write(*,*)  ' epe.pcr has been found'

  pcr_ex: if(ex_pcr.and.epe_embedding) then
  inquire (file=trim(input_dir)//'/epe.pcs',exist=ex_pcs)
  write(*,*)  ' epe.pcs has been found'
  inquire (file=trim(input_dir)//'/epe.pcc',exist=ex_pcc)
  write(*,*)  ' epe.pcc has been found'

  if(ex_pcs) then
   pcr_unit=openget_iounit(file=trim(input_dir)//'/epe.pcs',&
                               form='FORMATTED',status='old')
   read(pcr_unit,*)  pcs_N      !module variable
   write(output_unit,*)  ' Number of EPE shell pcs_N ',pcs_N
   allocate(pcs_temp(4,pcs_N),stat=ewa_allocstat(10))
   ASSERT(ewa_allocstat(10).eq.0)
   ewa_allocstat(10)=1

   do i=1,pcs_N
    read(pcr_unit,*) pcs_temp(:,i)
   enddo ! i=1,pcs_N
    call returnclose_iounit(pcr_unit)
  endif ! ex_pcs

  if(ex_pcc) then
   pcr_unit=openget_iounit(trim(input_dir)//'/epe.pcc',&
                          form='FORMATTED',status='old')
   read(pcr_unit,*)  pcc_N      !module variable
   write(output_unit,*) ' Number of EPE cores pcc_N ',pcc_N
   allocate(pcc_temp(4,pcs_N),stat=ewa_allocstat(11))
   ASSERT(ewa_allocstat(11).eq.0)
          ewa_allocstat(11)=1
        do i=1,pcc_N
         read(pcr_unit,*) pcc_temp(:,i)
        enddo ! i=1,pcc_N
    call returnclose_iounit(pcr_unit)
    if(pcc_N.ne.pcs_N) stop 'pcs_n ne pcc_n'
  endif ! ex_pcc

   nullify(pc)

   pcr_unit=openget_iounit(trim(input_dir)//'/epe.pcr',&
                          form='FORMATTED',status='old')
   read(pcr_unit,*)  pcr_N  !module variable
   allocate(pc,pcr_temp(4,pcr_N),pcr_no(pcr_N),psb_ind(pcr_N),stat=ewa_allocstat(9)) !4
   if(ewa_allocstat(9).ne.0) call error_handler("allocate pcr_temp failed")
      ewa_allocstat(9)=1 ! pcr_temp pcr_no psb_ind
      ewa_allocstat(16)=1 ! pc


    ! loop over all centers in file and sort out PC which positions coincide
    ! with the positions of regular atoms
     write(output_unit,*) 'Number of EPE reference points pcr_n ' ,pcr_n

    do i=1,pcr_n
        read(pcr_unit,*)  pcr_temp(:,i),index
        psb_ind(i)=index
    enddo ! i=1,pcr_n
    kl=0
!       n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)
    if(.not.allocated(gxepe_impu)) &
         call error_handler('gxepe_impu does not allocated, check data')
        n_gxat=0
        do i=1,N_unique_atoms
         if(gxepe_impu(i).ne.0) n_gxat=n_gxat+unique_atoms(i)%N_equal_atoms
        enddo
        write(output_unit,*) 'number of atoms in gx-file', n_gxat
        do i=1,n_timps
         if(gxepe_impu(i).ne.0) n_gxat=n_gxat+unique_timps(i)%n_equal_atoms
        enddo
        write(output_unit,*)'number of atoms and n_timps in gx file ',n_gxat

        if(.not.ex_gxepe) call error_handler(" file epe.r not found")
        md:    do i=1,pcr_n
           do i_ua=1,N_unique_atoms+n_timps
              if(gxepe_impu(i_ua).ne.0) then
                 if(i_ua.gt.N_unique_atoms) then
                    ua=>unique_timps(i_ua-N_unique_atoms)
                 else
                    ua=>unique_atoms(i_ua)
                 end if

                 do i_eq=1,ua%N_equal_atoms
                    if (dot_product(gxepe_array(i_ua)%position(:,i_eq)&
                         -pcr_temp(1:3,i), &
                         gxepe_array(i_ua)%position(:,i_eq)-pcr_temp(1:3,i)) &
                         .lt.small_pcr) then
                       kl=kl+1
                       ! swap pcr positions
                       pc%position_first_ec(:)=pcr_temp(1:3,kl)
                       pc%z=pcr_temp(4,kl)
                       pc%c=0.0_r8_kind
                       pcr_temp(1:3,kl)=pcr_temp(1:3,i)
                       pcr_temp(4,kl)=pcr_temp(4,i)
                       pcr_temp(1:3,i)=pc%position_first_ec
                       pcr_temp(4,i)=pc%z
                       index=psb_ind(kl)
                       psb_ind(kl)=psb_ind(i)
                       psb_ind(i)=index

                       pc%position_first_ec(:)=pcs_temp(1:3,kl)
                       pc%z=pcs_temp(4,kl)
                       pc%c=0.0_r8_kind
                       pcs_temp(1:3,kl)=pcs_temp(1:3,i)
                       pcs_temp(4,kl)=pcs_temp(4,i)
                       pcs_temp(1:3,i)=pc%position_first_ec
                       pcs_temp(4,i)=pc%z

                       pc%position_first_ec(:)=pcc_temp(1:3,kl)
                       pc%z=pcc_temp(4,kl)
                       pc%c=0.0_r8_kind
                       pcc_temp(1:3,kl)=pcc_temp(1:3,i)
                       pcc_temp(4,kl)=pcc_temp(4,i)
                       pcc_temp(1:3,i)=pc%position_first_ec
                       pcc_temp(4,i)=pc%z

                       if(kl.eq.n_gxat) exit md
                       cycle md
                    endif ! dot_product
                 enddo ! i_eq=1,ua%N_equal_atoms

              endif
           enddo ! i_ua=1,N_unique_atoms
        enddo md
        n_gxat_pcr=kl

        call returnclose_iounit(pcr_unit)

        write(output_unit,*) &
             'pcr_n after sorting out PC in atomic positions' ,pcr_n-kl

        ! check potential of PC
        pc%z=0.0_rk
        pc%c=0.0_rk
        e_nuc_pcr = 0.0_rk
        do nb=kl+1,pcr_n
           pc%z=pc%z+pcs_temp(4,nb)+pcc_temp(4,nb)
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 e_nuc_pcr = e_nuc_pcr+pcs_temp(4,nb)*unique_atoms(na)%Z/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                      -pcs_temp(1:3,nb))**2))
                 e_nuc_pcr = e_nuc_pcr+pcc_temp(4,nb)*unique_atoms(na)%Z/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                      -pcc_temp(1:3,nb))**2))
              enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo ! na=1,N_unique_atoms
        enddo ! 1,pcr_n

        write(output_unit,*) ' Energy of interaction of the EPE centers and '
        write(output_unit,*) ' atoms of cluster e_nuc_pcr calculsted  with use '
        write(output_unit,*) ' use of pcr_temp coordinates and  total charge'
        write(output_unit,*) ' of the EPE centers Z', e_nuc_pcr,pc%z

        i=1+kl
        it_pcr=0
    pcr_wh:    do while(i.le.pcr_n)
           pc%position_first_ec(1:3)=pcr_temp(1:3,i)
           pc%name='EPE PC '
           it_pcr=it_pcr+1

           position(1) = pc%position_first_ec(1)
           position(2) = pc%position_first_ec(3)
           position(3) = pc%position_first_ec(2)


           ! now apply all symmetry operations to the position of the
           ! unique atom
           n = 0
           do j=1,group_num_el
              position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
              if (dot_product(position2-position,position2-position) < small_pcr) &
                   then
                 n = n+1
              end if
           enddo


           ! allocate group elements
           local_groups_pc%num_el = n
           allocate(local_groups_pc%elements(n),stat=ewa_allocstat(12))
           ASSERT(ewa_allocstat(12).eq.0)
                  ewa_allocstat(12)=1

           ! fill up group elements
           n = 0
           do j=1,group_num_el
              position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
              if (dot_product(position2-position,position2-position) < small_pcr) then
                 n = n+1
                 local_groups_pc%elements(n) = j
              end if
           end do
           !
           ! now determine symmetry equivalent atoms
           !

           call group_coset_decomp(n_equal,local_groups_pc,&
                cosets_pc,point_trafos_pc%matrix)

                ewa_allocstat(13)=1 ! point_trafos_pc%matrix
                ewa_allocstat(14)=1 ! cosets_pc%elements

           pc%N_equal_charges=n_equal

           ! allocate positions of equal charges
           allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
           if (ewa_allocstat(15) .ne. 0 ) call error_handler( &
                "symm_equivalents_gen: allocation of pc%position failed")
               ewa_allocstat(15)=1

           ! determine positions of equal atoms
           do j=1,n_equal
              position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
              pc%position(1,j) = position2(1)
              pc%position(2,j) = position2(3)
              pc%position(3,j) = position2(2)
              do ii=i,pcr_n
                 if(dot_product(pcr_temp(1:3,ii)-pc%position(:,j) &
                      ,pcr_temp(1:3,ii)-pc%position(:,j)) < small_pcr) then

                    pcr_temp(1:3,ii)=pcr_temp(1:3,i)
                    pcr_temp(1:3,i)=pc%position(1:3,j)
                    pc%z=pcr_temp(4,ii)
                    pc%c=0.0_r8_kind
                    pcr_temp(4,ii)=pcr_temp(4,i)
                    pcr_temp(4,i)=pc%z
                    pcr_no(i)=it_pcr
                    index=psb_ind(ii)
                    psb_ind(ii)=psb_ind(i)
                    psb_ind(i)=index

                    pc%position(1:3,j)=pcs_temp(1:3,ii)
                    pc%z=pcs_temp(4,ii)
                    pc%c=0.0_r8_kind
                    pcs_temp(1:3,ii)=pcs_temp(1:3,i)
                    pcs_temp(4,ii)=pcs_temp(4,i)
                    pcs_temp(1:3,i)=pc%position(1:3,j)
                    pcs_temp(4,i)=pc%z

                    pc%position(1:3,j)=pcc_temp(1:3,ii)
                    pc%z=pcc_temp(4,ii)
                    pc%c=0.0_r8_kind
                    pcc_temp(1:3,ii)=pcc_temp(1:3,i)
                    pcc_temp(4,ii)=pcc_temp(4,i)
                    pcc_temp(1:3,i)=pc%position(1:3,j)
                    pcc_temp(4,i)=pc%z
                    exit
                 endif !
              enddo ! ii=i+1,pcr_n
              if(i>pcr_n) stop ' i>pcr_n, check PCR array'
              i=i+1
           enddo ! j=1,n_equal

           deallocate(local_groups_pc%elements, &
                      point_trafos_pc%matrix, &
                      cosets_pc%elements, &
                      pc%position, &
                      stat=ewa_allocstat(12))      ! local_groups_pc%elements
                ASSERT(ewa_allocstat(12).eq.0)
                ewa_allocstat(13)=0 ! point_trafos_pc%matrix
                ewa_allocstat(14)=0 ! cosets_pc%elements
                ewa_allocstat(15)=0 ! pc%position

           if (n_equal .ne. pc%N_equal_charges) call error_handler( &
                "symm_equivalents_gen: calculate number of equal charges and input number differ" )
    enddo  pcr_wh! while

    i=1+kl
    it_pcr=0  !types of EWPC
    pcs_wh: do while(i.le.pcr_n)
       pc%position_first_ec(1:3)=pcs_temp(1:3,i)
       pc%name='EPE PCS '
       it_pcr=it_pcr+1

       position(1) = pc%position_first_ec(1)
       position(2) = pc%position_first_ec(3)
       position(3) = pc%position_first_ec(2)


       ! now apply all symmetry operations to the position of the
       ! unique atom
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small_pcr) then
             n = n+1
          end if
       enddo


       ! allocate group elements
       local_groups_pc%num_el = n
       allocate(local_groups_pc%elements(n),stat=ewa_allocstat(12))
       ASSERT(ewa_allocstat(12).eq.0)
              ewa_allocstat(12)=1

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small_pcr) then
             n = n+1
             local_groups_pc%elements(n) = j
          end if
       end do
       !
       ! now determine symmetry equivalent atoms
       !

       call group_coset_decomp(n_equal,local_groups_pc,&
            cosets_pc,point_trafos_pc%matrix)

            ewa_allocstat(13)=1 ! point_trafos_pc%matrix
            ewa_allocstat(14)=1 ! cosets_pc%elements

       pc%N_equal_charges=n_equal

       ! allocate positions of equal charges
       allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))

       ASSERT(ewa_allocstat(15).eq.0)
              ewa_allocstat(15)=1

       ! determine positions of equal atoms
       do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
          pc%position(1,j) = position2(1)
          pc%position(2,j) = position2(3)
          pc%position(3,j) = position2(2)
        do ii=i,pcr_n
        if(dot_product(pcs_temp(1:3,ii)-pc%position(:,j) &
                      ,pcs_temp(1:3,ii)-pc%position(:,j)) < small_pcr) then

        pcs_temp(1:3,ii)=pcs_temp(1:3,i)
        pcs_temp(1:3,i)=pc%position(1:3,j)
        pc%z=pcs_temp(4,ii)
        pc%c=0.0_r8_kind
        pcs_temp(4,ii)=pcs_temp(4,i)
        pcs_temp(4,i)=pc%z
        exit
        endif !
        enddo ! ii=i+1,pcr_n
        if(i>pcr_n) stop ' i>pcr_n, check PCR array'
        i=i+1
        enddo ! j=1,n_equal

       deallocate(local_groups_pc%elements, &
            cosets_pc%elements, &
            pc%position, &
            point_trafos_pc%matrix,stat=ewa_allocstat(12))
       if (ewa_allocstat(12) .ne. 0 ) call error_handler( &
            "symm_equivalents_gen: deallocation of point charge helpers failed")
           ewa_allocstat(13)=0 ! point_trafos_pc%matrix
           ewa_allocstat(14)=0 ! cosets_pc%elements
           ewa_allocstat(15)=0 ! pc%position

       if (n_equal .ne. pc%N_equal_charges) call error_handler( &
            "symm_equivalents_gen: calculate number of equal charges and input number differ" )
    enddo pcs_wh ! while

     i=1+kl
     it_pcr=0  !types of EWPC
    pcc_wh: do while(i.le.pcr_n)
       pc%position_first_ec(1:3)=pcc_temp(1:3,i)
       pc%name='EPE PCC '
       it_pcr=it_pcr+1

       position(1) = pc%position_first_ec(1)
       position(2) = pc%position_first_ec(3)
       position(3) = pc%position_first_ec(2)


       ! now apply all symmetry operations to the position of the
       ! unique atom
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
       if (dot_product(position2-position,position2-position) < small_pcr) then
             n = n+1
          end if
       enddo


       ! allocate group elements
       local_groups_pc%num_el = n
       allocate(local_groups_pc%elements(n),stat=ewa_allocstat(12))
       ASSERT(ewa_allocstat(12).eq.0)
       ewa_allocstat(12)=1

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
      if (dot_product(position2-position,position2-position) < small_pcr) then
             n = n+1
             local_groups_pc%elements(n) = j
          end if
       end do
       !
       ! now determine symmetry equivalent atoms
       !

       call group_coset_decomp(n_equal,local_groups_pc,&
            cosets_pc,point_trafos_pc%matrix)

            ewa_allocstat(13)=1 ! point_trafos_pc%matrix
            ewa_allocstat(14)=1 ! cosets_pc%elements

       pc%N_equal_charges=n_equal

       ! allocate positions of equal charges
       allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
       if (alloc_stat .ne. 0 ) call error_handler( &
            "symm_equivalents_gen: allocation of pc%position failed")

       ASSERT(ewa_allocstat(15).eq.0)
              ewa_allocstat(15)=1

       ! determine positions of equal atoms
       do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
          pc%position(1,j) = position2(1)
          pc%position(2,j) = position2(3)
          pc%position(3,j) = position2(2)
        do ii=i,pcr_n
        if(dot_product(pcc_temp(1:3,ii)-pc%position(:,j) &
                      ,pcc_temp(1:3,ii)-pc%position(:,j)) < small_pcr) then

        pcc_temp(1:3,ii)=pcc_temp(1:3,i)
        pcc_temp(1:3,i)=pc%position(1:3,j)
        pc%z=pcc_temp(4,ii)
        pc%c=0.0_r8_kind
        pcc_temp(4,ii)=pcc_temp(4,i)
        pcc_temp(4,i)=pc%z
        exit
        endif !
        enddo ! ii=i+1,pcr_n
        if(i>pcr_n) stop ' i>pcr_n, check PCR array'
        i=i+1
        enddo ! j=1,n_equal

       deallocate(local_groups_pc%elements, &
            pc%position, &
            cosets_pc%elements, &
            point_trafos_pc%matrix, &
                      stat=ewa_allocstat(12)) ! local_groups_pc%elements
       if (ewa_allocstat(12) .ne. 0 ) call error_handler( &
            "symm_equivalents_gen: deallocation of point charge helpers failed")
           ewa_allocstat(13)=0 ! point_trafos_pc%matrix
           ewa_allocstat(14)=0 ! cosets_pc%elements
           ewa_allocstat(15)=0 ! pc%position

       if (n_equal .ne. pc%N_equal_charges) call error_handler( &
            "symm_equivalents_gen: calculate number of equal charges and input number differ" )
    enddo pcc_wh ! while

      ! check potential of PC
        e_nuc_pcr = 0.0_rk
        deallocate(pc,stat=ewa_allocstat(16))
        nullify(pc)

        allocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)
        ewa_allocstat(16)=1

        pc%z=0.0_rk
        do nb=1+kl,pcr_n

        if(pcs_n.ne.0) then
           pc%z=pc%z+pcs_temp(4,nb)+pcc_temp(4,nb)
             do na=1,N_unique_atoms
                   do eq_a=1,unique_atoms(na)%N_equal_atoms
                      e_nuc_pcr = e_nuc_pcr+pcs_temp(4,nb)*unique_atoms(na)%Z/ &
                           sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                           -pcs_temp(1:3,nb))**2))
                      e_nuc_pcr = e_nuc_pcr+pcc_temp(4,nb)*unique_atoms(na)%Z/ &
                           sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                           -pcc_temp(1:3,nb))**2))
                   enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
             enddo ! na=1,N_unique_atoms

         else
           pc%z=pc%z+pcr_temp(4,nb)
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 dist = sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                      -pcr_temp(1:3,nb))**2))
                 e_nuc_pcr = e_nuc_pcr+pcr_temp(4,nb)*unique_atoms(na)%Z/dist
              enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo ! na=1,N_unique_atoms
      endif ! pcs_n.ne.0/else
     enddo ! 1,pcr_n
     write(output_unit,*) 'e_nuc_pcr  and Z after enforced symmetrization ', &
          e_nuc_pcr,pc%z

     if(pcs_n.eq.0.and.pcc_n.eq.0) then
        allocate(pcr_array(it_pcr),stat=ewa_allocstat(17))
        ASSERT(ewa_allocstat(17).eq.0)
               ewa_allocstat(17)=1

     else
      if(pcs_n.ne.0) then
         allocate(pcs_array(it_pcr),stat=ewa_allocstat(18))
         ASSERT(ewa_allocstat(18).eq.0)
         ewa_allocstat(18)=1
      endif

      if(pcc_n.ne.0) then
        allocate(pcc_array(it_pcr),stat=ewa_allocstat(19))
         ASSERT(ewa_allocstat(19).eq.0)
         ewa_allocstat(19)=1
      endif

     endif ! pcs_n.eq.0.and.pcc_n.eq.0/else

     z_pcs_pcc: if(pcs_n.eq.0.or.pcc_n.eq.0) then
        i=1+kl
        do while(i.le.pcr_n)
           pc=>pcr_array(pcr_no(i))
           n_equal=0
           pc%z=0.0_rk

           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcr_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n

           pc%z=pc%z/n_equal
           pc%c=0.0_rk

           allocate(pc%position(3,n_equal),stat=ewa_allocstat(24))
           ASSERT(ewa_allocstat(24).eq.0)
                  ewa_allocstat(24)=1 ! pcr_array(:)%position

           pc%name='pcr    '
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcr_temp(1:3,i+ii-1)
!!!     print all generated charges
!!!     print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$              write(output_unit, '(4f15.8,i3,4i2)') pc%position(:,ii), &
!!$                   pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
           dist=sqrt(sum((unique_atoms(1)%position(:,1) &
                -pc%position(:,1))**2))
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)' ) &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0,dist
           i=i+n_equal

        enddo !while

     else z_pcs_pcc! basic option

!!$        write(output_unit,*) '** treat regular positions'
!!$        i=1+kl
!!$        do while(i.le.pcr_n)
!!$           pc=>pcr_array(pcr_no(i))
!!$           n_equal=0
!!$           pc%z=0.0_rk
!!$
!!$           do ii=i,pcr_n
!!$              n_equal=n_equal+1
!!$              pc%z=pc%z+pcr_temp(4,ii)
!!$              if(ii.eq.pcr_n) exit
!!$              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
!!$           enddo !ii=i,pcr_n
!!$
!!$           pc%z=pc%z/n_equal
!!$
!!$           allocate(pc%position(3,n_equal),stat=alloc_stat)
!!$           pc%name='pcr    '
!!$           pc%N_equal_charges= n_equal
!!$           do ii=1,n_equal
!!$              pc%position(:,ii)=pcr_temp(1:3,i+ii-1)
!!$!!!  print all generated charges
!!$!!!          print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$              write(output_unit, '(4f15.8,i3,4i2)') pc%position(:,ii), &
!!$                   pc%z,n_equal,0,0,1,0
!!$           enddo ! ii=1,equal
!!$           dist=sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
!!$           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)' ) &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0,dist
!!$           i=i+n_equal
!!$
!!$        enddo !while

        i=1+kl
        nullify(pc)
        fill_pcs:do while(i.le.pcr_n)
           pc=>pcs_array(pcr_no(i))
           n_equal=0

           pc%z=0.0_rk
           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcs_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n

           pc%z=pc%z/n_equal
           pc%c=0.0_r8_kind

           allocate(pc%position(3,n_equal),stat=ewa_allocstat(25))
           ASSERT(ewa_allocstat(25).eq.0)
                  ewa_allocstat(25)=1 ! pcs_array(:)%position

           pc%name='pcs    '
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcs_temp(1:3,i+ii-1)
!!!     print all generated charges
!!!     print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$              write(output_unit, '(4f15.8,i3,4i2)')&
!!$                   pc%position(:,ii),pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)'  &
!!$                ) pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$                  sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
           i=i+n_equal
        enddo fill_pcs
        !** done

        nullify(pc)
        i=1+kl
        fill_pcc: do while(i.le.pcr_n)
           pc=>pcc_array(pcr_no(i))
           n_equal=0
           pc%z=0.0_rk

           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcc_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n

           pc%z=pc%z/n_equal
           pc%c=0.0_r8_kind

           allocate(pc%position(3,n_equal),stat=ewa_allocstat(26))
                    ASSERT(ewa_allocstat(26).eq.0)
                           ewa_allocstat(26)=1
           pc%name='pcc    '
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcc_temp(1:3,i+ii-1)
!!!           print all generated charges
!!!           print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$           write(output_unit, '(4f15.8,i3,4i2)') pc%position(:,ii), &
!!$                   pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)' ) &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$                sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
           i=i+n_equal
        enddo fill_pcc !while


    endif z_pcs_pcc!    (pcs_n.eq.0.or.pcc_n.eq.0) then

    write(output_unit,*) 'number of groups of symmetry equavalent PC', &
             it_pcr
    if(pcs_n.ne.0) pcs_n=it_pcr
    if(pcc_n.ne.0) pcc_n=it_pcr
    n_unique_pcr=it_pcr

        if(pcs_n.ne.0) then
           deallocate(pcs_temp,stat=ewa_allocstat(10))
           ASSERT(ewa_allocstat(10).eq.0)
        endif
        if(pcc_n.ne.0) then
           deallocate(pcc_temp,stat=ewa_allocstat(11))
           ASSERT(ewa_allocstat(11).eq.0)
        endif
     endif pcr_ex ! ex_pcr.and.epe_embedding
 end subroutine symm_epe
#endif

#if SYMMEWPCGEN
  subroutine symm_ewpc_gen
    use filename_module
    use pointcharge_module
    use operations_module
! treats PC of ewald program generated file
  use ewaldpc_module, ex_pcew_dumm=>ex_pcew
  type(unique_atom_type), pointer  :: ua
  logical:: ex_pcew
  real (RK):: e_nuc_ewpc,dist
  integer(IK):: na,nb,eq_a,kl_ew,kl1,kl_pcr,n_gxat,kl_ew_eq
  integer(IK)::  it_ewpc,it_pcr,ii,i_ua,i_eq,eq_pcc,eq_pcs
  integer(IK),dimension(:),allocatable:: pcew_no
  integer(IK):: status
  real(RK), allocatable, dimension(:,:) ::pcew_temp
  real(RK),parameter:: small_ew=1.d-2

  inquire (file= trim(input_dir)//'/ewald.pcr',exist=ex_pcew)

  write(*,*)  ' ewald.pcr has been found'
  if(ex_pcew.and.ewaldpc) then

   !called after symm_epe and thus ponter was  on  pcc_array pcs_array

   ewpc_unit=openget_iounit(trim(input_dir)//'/ewald.pcr'&
                           ,form='FORMATTED',status='old')
   read(ewpc_unit,*)  EWPC_N
   nullify(pc)
   allocate(pc,pcew_temp(4,EWPC_N),pcew_no(EWPC_N),stat=ewa_allocstat(21))
           ASSERT(ewa_allocstat(21).eq.0)
                  ewa_allocstat(21)=1 ! pcew_temp
                  ewa_allocstat(27)=1 ! pcew_no

                  ewa_allocstat(16)=1 ! PC as a temp


   !** loop over all centers in file and sort out PC which positions coincide
   !   with the positions of regular atoms of QM cluster
   write(output_unit,*) &
        ' Number of centers generated with program ewald,  EWPC_N ' ,EWPC_N
   if(ex_gxepe) then
    write(output_unit,*) &
        'as file epe.r with regular atomic positions is located'
    write(output_unit,*) &
        ' EWPC centers coinciding with atoms will be sorted out'
!!$      write(output_unit,*)  gxepe_array(1)%position(:,1)
   end if

   ii=1
   it_ewpc=EWPC_N ! it_ewpc here is help variable to define number of PC
   mc:    do i=1,EWPC_N
      read(ewpc_unit,*)  pcew_temp(:,ii)
!!$      if(ii.eq.1) write(output_unit,*) pcew_temp(:,ii)
      do i_ua=1,N_unique_atoms+n_timps
         if(ex_gxepe) then
            if(i_ua.gt.N_unique_atoms) then
               do i_eq=1,unique_timps(i_ua-N_unique_atoms)%N_equal_atoms
                  if (dot_product(gxepe_array(i_ua)%position(:,i_eq) &
                       -pcew_temp(1:3,ii), &
                       gxepe_array(i_ua)%position(:,i_eq)-pcew_temp(1:3,ii)) &
                       .lt.small_ew) then
                     it_ewpc=it_ewpc-1
                     cycle  mc
                  endif ! dot_product
               end do
            else
               do i_eq=1,unique_atoms(i_ua)%N_equal_atoms
                  if (dot_product(gxepe_array(i_ua)%position(:,i_eq)- &
                       pcew_temp(1:3,ii), &
                       gxepe_array(i_ua)%position(:,i_eq)-pcew_temp(1:3,ii)) &
                       .lt.small_ew) then
                     it_ewpc=it_ewpc-1
                     cycle  mc
                  endif ! dot_product
               enddo ! i_eq=1,unique_atoms(i_ua)%N_equal_atoms
            end if


           else ! .not.ex_gxepe
              do i_eq=1,unique_atoms(i_ua)%N_equal_atoms
                 if (dot_product(unique_atoms(i_ua)%position(:,i_eq)- &
                      pcew_temp(1:3,ii), &
                      unique_atoms(i_ua)%position(:,i_eq)-pcew_temp(1:3,ii)) &
                      .lt.small_ew) then
                    it_ewpc=it_ewpc-1
                    cycle  mc
                 endif !
                 enddo
              endif ! ex_gxepe/else
           enddo
           ii=ii+1
        enddo mc
        call returnclose_iounit(ewpc_unit)
        !** now pcew_temp contains only PC which do not related to centers of QM cluster

        EWPC_N=it_ewpc
        write(output_unit,*) 'EWPC_N after sorting out PC in atomic positions' ,&
        EWPC_N
        if(ex_gxepe) then
         do i_ua=1,N_unique_atoms
           deallocate(gxepe_array(i_ua)%position,stat=ewa_allocstat(23))
           ASSERT(ewa_allocstat(23).eq.0)
         end do
           deallocate(gxepe_array,gxepe_impu,stat=ewa_allocstat(22))
           ASSERT(ewa_allocstat(23).eq.0)
         endif ! ex_gxepe
!** done

        kl_ew=0
        kl_ew_eq=0
        if(pcr_n.ne.0) then  ! if epe charges exist
           if(n_timps.gt.0) then
              n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)+&
                   sum(unique_timps(1:n_timps)%N_equal_atoms)
           else
              n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)
           end if

           print*,'n_gxat pcr_n n_gxat_pcr ',n_gxat, pcr_n, n_gxat_pcr
           kl_pcr=n_gxat_pcr
           write(output_unit,*) 'initial value  of kl_pcr',kl_pcr
           do i=1,EWPC_N
              kl1=kl_pcr+1
              do k=kl1,pcr_n
                 if(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i)).lt.small_ew &
                      .and. ((abs(pcr_temp(4,k)-pcew_temp(4,i)).lt.0.01_rk) &
                      .or.(OPERATIONS_PSEUDOBONDS.and.psb_ind(k).ne.0)) ) then
                    !**** centers coincide and charges are equal
                    ! such centers generated with ewald can be neglected
                    kl_pcr=kl_pcr+1
                    kl_ew=kl_ew+1
                    pc%position_first_ec=pcr_temp(1:3,kl_pcr)
                    pcr_temp(1:3,kl_pcr)=pcr_temp(1:3,k)
                    pcr_temp(1:3,k)=pc%position_first_ec
                    pc%z=pcr_temp(4,kl_pcr)
                    pc%c=0.0_r8_kind
                    pcr_temp(4,kl_pcr)=pcr_temp(4,k)
                    pcr_temp(4,k)=pc%z
                    index=psb_ind(kl_pcr)
                    psb_ind(kl_pcr)=psb_ind(k)
                    psb_ind(k)=index

                    pc%position_first_ec=pcew_temp(1:3,kl_ew)
                    pcew_temp(1:3,kl_ew)=pcew_temp(1:3,i)
                    pcew_temp(1:3,i)=pc%position_first_ec

                    pc%z=pcew_temp(4,kl_ew)
                    pcew_temp(4,kl_ew)=pcew_temp(4,i)
                    pcew_temp(4,i)=pc%z
                    exit
                 endif ! dot_product
              enddo ! kl1,pcr_n
           enddo ! i=1,EWPC_N
           write(output_unit,*) 'number of found coinciding in charge centers in pcr and ewpc', &
                kl_pcr,kl_ew
           write(output_unit,*) 'indexes  for  last coinciding elements ' ,kl_ew,kl_pcr
           kl_ew_eq=kl_ew ! centers starting from this index will not be neglected
           ! but treated with modified charge now this centers are sorted
           ! and the charges are modified

           !locate centers which coinside but have different charges
           do i=kl_ew+1,EWPC_N
              kl1=kl_pcr+1
              do k=kl1,pcr_n
                 if(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i)).lt.small_ew )  then
                    kl_pcr=kl_pcr+1
                    kl_ew=kl_ew+1
                    pc%position_first_ec=pcr_temp(1:3,kl_pcr)
                    pcr_temp(1:3,kl_pcr)=pcr_temp(1:3,k)
                    pcr_temp(1:3,k)=pc%position_first_ec
                    pc%z=pcr_temp(4,kl_pcr)
                    pc%c=0.0_r8_kind
                    pcr_temp(4,kl_pcr)=pcr_temp(4,k)
                    pcr_temp(4,k)=pc%z

                    pc%position_first_ec=pcew_temp(1:3,kl_ew)
                    pcew_temp(1:3,kl_ew)=pcew_temp(1:3,i)
                    pcew_temp(1:3,i)=pc%position_first_ec

                    pc%z=pcew_temp(4,kl_ew)
                    pc%c=0.0_rk
                    pcew_temp(4,kl_ew)=pcew_temp(4,i)
                    pcew_temp(4,i)=pc%z
                    pcew_temp(4,kl_ew)=pcew_temp(4,kl_ew)-pcr_temp(4,kl_pcr)
                    pcr_temp(4,kl_pcr)=0.0_rk
                    exit
                 endif! dot_product
              enddo! kl1,pcr_n
           enddo! i=1,EWPC_N
           write(output_unit,*) &
                'number of coinciding centers in ewpc and pcr arrays & charge ',kl_ew &
                ,sum(pcew_temp(4,1:kl_ew)),sum(pcr_temp(4,1:n_gxat:kl_pcr))

        endif! pcr_n.ne.0

        ! **   now first kl centers in pcew_temp and pcr_temp coincide

        ! **  check potential of PC
        e_nuc_ewpc = 0.0_rk
        do nb=1+kl_ew,EWPC_N
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 e_nuc_ewpc = e_nuc_ewpc+pcew_temp(4,nb)*unique_atoms(na)%Z/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                      -pcew_temp(1:3,nb))**2))
              enddo!  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo! na=1,N_unique_atoms
        enddo! 1,EWPC_N
        write(output_unit,*) 'e_nuc_ewpc with use of pcew_temp and Z', e_nuc_ewpc, &
             sum(pcew_temp(4,1+kl_ew:EWPC_N))
        !** done

   if(pcr_n.ne.0) then
      write(output_unit,*) '** take additional centers in pcew_temp',1+kl_ew
      else
         print*,'take PC from pcew_temp starting from ',1+kl_ew
      end if

   i=1+kl_ew_eq
   it_ewpc=0    !types of EWPC
   ewpc_wh: do while(i.le.EWPC_N)
      pc%position_first_ec(1:3)=pcew_temp(1:3,i)
      pc%name='ewald PC '
      it_ewpc=it_ewpc+1
      position(1) = pc%position_first_ec(1)
      position(2) = pc%position_first_ec(3)
      position(3) = pc%position_first_ec(2)
      ! now apply all symmetry operations to the position of the
      ! unique atom
      n = 0
      do j=1,group_num_el
         position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
         if (dot_product(position2-position,position2-position) < small_ew) then
            n = n+1
         end if
      enddo
      ! allocate group elements
      local_groups_pc%num_el = n
      allocate(local_groups_pc%elements(n),stat=ewa_allocstat(12))
      ASSERT(ewa_allocstat(12).eq.0)
      ewa_allocstat(12)=1

      ! fill up group elements
      n = 0
      do j=1,group_num_el
         position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
         if (dot_product(position2-position,position2-position) < small_ew) then
            n = n+1
            local_groups_pc%elements(n) = j
         end if
      end do
      !
      ! determine symmetry equivalent atoms
      !

      call group_coset_decomp(n_equal,local_groups_pc,&
           cosets_pc,point_trafos_pc%matrix)
      pc%N_equal_charges=n_equal

      ewa_allocstat(13)=1 ! point_trafos_pc%matrix
      ewa_allocstat(14)=1 ! cosets_pc%elements


      ! allocate positions of equal charges
      allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
      ASSERT(ewa_allocstat(15).eq.0)
      ewa_allocstat(15)=1 ! pc%position in temp

      ! determine positions of equal atoms
      !       write(output_unit,*) " Equal charges of Type ",i
      do j=1,n_equal
         position2 = &
              MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
         pc%position(1,j) = position2(1)
         pc%position(2,j) = position2(3)
         pc%position(3,j) = position2(2)
         do ii=i,EWPC_N
            if(dot_product(pcew_temp(1:3,ii)-pc%position(:,j) &
                 ,pcew_temp(1:3,ii)-pc%position(:,j)) < small_ew) then

               pcew_temp(1:3,ii)=pcew_temp(1:3,i)
               pcew_temp(1:3,i)=pc%position(1:3,j)
               pc%z=pcew_temp(4,ii)
               pc%c=0.0_r8_kind
               pcew_temp(4,ii)=pcew_temp(4,i)
               pcew_temp(4,i)=pc%z
               pcew_no(i)=it_ewpc

               exit
            endif !
         enddo ! ii=i+1,EWPC_N
         if(i>EWPC_N) stop ' i>EWPC_N, check EWPC array'
         i=i+1
      enddo ! j=1,n_equal

      deallocate(local_groups_pc%elements, point_trafos_pc%matrix, &
                 cosets_pc%elements, &
                 pc%position, stat=ewa_allocstat(12))
      ASSERT(ewa_allocstat(12).eq.0) ! local_groups_pc%elements
      ewa_allocstat(13)=0            ! point_trafos_pc%matrix
      ewa_allocstat(14)=0            ! coset%elements
      ewa_allocstat(15)=0            ! pc%position in temp

      if (n_equal .ne. pc%N_equal_charges) call error_handler( &
           "symm_equivalents_gen: calculate number of equal charges and input number differ" )
    enddo ewpc_wh ! while
    !** done

    !**  check potential of PC
    e_nuc_ewpc = 0.0_rk
    do nb=1+kl_ew,EWPC_N
       do na=1,N_unique_atoms
          do eq_a=1,unique_atoms(na)%N_equal_atoms
             e_nuc_ewpc = e_nuc_ewpc+pcew_temp(4,nb)*unique_atoms(na)%Z/ &
                  sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                  -pcew_temp(1:3,nb))**2))
          enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
       enddo ! na=1,N_unique_atoms
    enddo ! 1,EWPC_N
    write(output_unit,*) 'e_nuc_ewpc with use of pcew_temp & Z', e_nuc_ewpc, &
         sum(pcew_temp(4,1+kl_ew:EWPC_N))
    !** done
    it_pcr=0    !types of pcr
    epe_ex:if(pcr_n.gt.kl_pcr) then
       write(output_unit,*) &
            ' not coinciding with ewpc environment epe centers  exist in  No=', &
            pcr_n-kl_pcr
       write(output_unit,*) 'charge of these PC ', sum(pcr_temp(4,1+kl_pcr:pcr_n))
       !** take additional centers in pcr_temp
       i=1+kl_pcr
      pcr_wh: do while(i.le.pcr_n)
          pc%position_first_ec(1:3)=pcr_temp(1:3,i)
          pc%name='pcr '
          it_pcr=it_pcr+1

          position(1) = pc%position_first_ec(1)
          position(2) = pc%position_first_ec(3)
          position(3) = pc%position_first_ec(2)
          ! now apply all symmetry operations to the position of the
          ! unique atom
          n = 0
          do j=1,group_num_el
             position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
             if (dot_product(position2-position,position2-position) < small_ew) then
                n = n+1
             end if
          enddo
          ! allocate group elements
          local_groups_pc%num_el = n
          allocate(local_groups_pc%elements(n),stat=ewa_allocstat(12))
          ASSERT(ewa_allocstat(12).eq.0)
          ewa_allocstat(12)=1
          ! fill up group elements
          n = 0
          do j=1,group_num_el
             position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
             if (dot_product(position2-position,position2-position) < small_ew) then
                n = n+1
                local_groups_pc%elements(n) = j
             end if
          end do
          !
          ! determine symmetry equivalent atoms
          !
          call group_coset_decomp(n_equal,local_groups_pc,&
               cosets_pc,point_trafos_pc%matrix)
          ewa_allocstat(13)=1 ! point_trafos_pc%matrix
          ewa_allocstat(14)=1 ! cosets_pc

          pc%N_equal_charges=n_equal

          ! allocate positions of equal charges
          allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
          ASSERT(ewa_allocstat(15).eq.0)
          ewa_allocstat(15)=1

          ! determine positions of equal atoms
          do j=1,n_equal
             position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
             pc%position(1,j) = position2(1)
             pc%position(2,j) = position2(3)
             pc%position(3,j) = position2(2)
             do ii=i,pcr_n
                if(dot_product(pcr_temp(1:3,ii)-pc%position(:,j) &
                     ,pcr_temp(1:3,ii)-pc%position(:,j)) < small_ew) then

                   pcr_temp(1:3,ii)=pcr_temp(1:3,i)
                   pcr_temp(1:3,i)=pc%position(1:3,j)
                   pc%z=pcr_temp(4,ii)
                   pc%c=0.0_r8_kind
                   pcr_temp(4,ii)=pcr_temp(4,i)
                   pcr_temp(4,i)=pc%z
                   pcr_no(i)=it_pcr
                   exit
                endif !
             enddo ! ii=i+1,pcr_n
             if(i>pcr_n) stop ' i>EWPC_N, check EWPC array'
             i=i+1
          enddo ! j=1,n_equal

          deallocate(local_groups_pc%elements, point_trafos_pc%matrix, &
                     cosets_pc%elements,&
                     pc%position,                stat=ewa_allocstat(13))
          ASSERT(ewa_allocstat(13).eq.0) ! point_trafos_pc%matrix
          ewa_allocstat(12)=0 ! local_groups_pc%elements
          ewa_allocstat(14)=0 ! cosets_pc%elements
          ewa_allocstat(15)=0 ! pc%position

          if (n_equal .ne. pc%N_equal_charges) call error_handler( &
               "symm_equivalents_gen: calculate number of equal charges and input number differ" )
    enddo pcr_wh! while
    !** done
    endif epe_ex

    n_unique_pcr=pcs_n+pcc_n
    if(pcr_n.gt.kl_pcr) n_unique_pcr=n_unique_pcr+it_pcr
    allocate(ewpc_array(it_ewpc+n_unique_pcr),ewpc_arrel_used(it_ewpc+n_unique_pcr), &
                                                                stat=ewa_allocstat(3))
        if(ewa_allocstat(3).ne.0) call error_handler("allocate of ewpc_array failed")
        ewa_allocstat(3)=1
        ewpc_arrel_used=.false.
    DPRINT 'size of ewpc_array',size(ewpc_array)

        deallocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)

!!$        write(output_unit,  *)  '** store not coinciding pcew centers'
        i=1+kl_ew_eq
       ewpcarr_wh: do while(i.le.EWPC_N)
           pc=>ewpc_array(pcew_no(i))
           n_equal=0
           pc%z=0.0_rk
           do ii=i,ewpc_N
              n_equal=n_equal+1
              pc%z=pc%z+pcew_temp(4,ii)
              if(ii.eq.ewpc_N) exit
              if(pcew_no(ii+1).ne.pcew_no(ii)) exit
           enddo !ii=i,ewpc_N
           pc%z=pc%z/n_equal
           pc%c=0.0_r8_kind
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(4))
           ASSERT(ewa_allocstat(4).eq.0)
           ewa_allocstat(4)=1
           ewpc_arrel_used(pcew_no(i))=.true.
           pc%name='ewpc    '
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcew_temp(1:3,i+ii-1)
!!!     print all generated charges
!!!     print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$              write(output_unit, '(4f15.8,i3,4i2)') pc%position(:,ii), &
!!$                   pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)')  &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$             sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8,i4)') &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$                sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2)),(pcew_no(i))
        i=i+n_equal
        enddo ewpcarr_wh

        write(output_unit,*) &
             'number of groups of symmetry equavalent ewpc PC',it_ewpc
       ewpc_n=it_ewpc

        !** done
     pcr_n_ex:   if(pcr_n.ne.0) then
        write(output_unit,*) &
             '** store not coinciding pcr centers'
        print*,'** store not coinciding pcr centers'
        i=1+kl_pcr
       nullify(pc)
       ewpcarr_pcr: do while(i.le.pcr_n)
           pc=>ewpc_array(pcr_no(i)+ewpc_n)
           n_equal=0
           pc%z=0.0_rk
           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z-pcr_temp(4,ii) ! negative charge for virtual PC
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo
           pc%z=pc%z/n_equal
           pc%c=0.0_r8_kind
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(4))
           ASSERT(ewa_allocstat(4).eq.0)
           ewa_allocstat(4)=1
           ewpc_arrel_used(pcr_no(i)+ewpc_n)=.true.
           pc%name='pcr    '
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcr_temp(1:3,i+ii-1)
!!!     print all generated charges
!!!     print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!!     write(output_unit, '(4f15.8,i3,4i2)'),pc%position(:,ii),&
!!!                                             pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)')  &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$             sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
!!$           print '(4f15.8,i3,4i2,f15.8)',pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!!$                sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
           i=i+n_equal

        enddo  ewpcarr_pcr
         write(output_unit,*) &
              'number of groups of symmetry equavalent pcr PC',it_pcr
         ewpc_n=ewpc_n+it_pcr
        !** done
     endif pcr_n_ex


    start_regI_epevar=ewpc_n+1
    cluster_nuc_epe_en=0.0_rk
    if(pcs_n.ne.0) then !** store pcs centers and calculate
                        ! cluster_nuc_epe_en for epe shells
       nullify(pc)
       do na=1,pcs_n
          pc=>ewpc_array(na+ewpc_n)
          pc%n_equal_charges=pcs_array(na)%n_equal_charges
          allocate(pc%position(3,pc%n_equal_charges),stat=ewa_allocstat(4))
          ASSERT(ewa_allocstat(4).eq.0)
          ewa_allocstat(4)=1
          pc=pcs_array(na)
          ewpc_arrel_used(na+ewpc_n)=.true.
          DPRINT 'pcs ewpc_arrel_used', na+ewpc_n
          if(na+ewpc_n.eq.413) then
           DPRINT 'sum position',sum(pc%position)
          endif

          do eq_pcs=1,pc%n_equal_charges
             do i=1,N_unique_atoms+n_timps
                if(i.gt.N_unique_atoms) then
                   ua=>unique_timps(i-N_unique_atoms)
                   do eq_a=1,ua%N_equal_atoms
                      dist=sqrt(sum((ua%position(:,eq_a) &
                           -pc%position(:,eq_pcs))**2))
                      cluster_nuc_epe_en=cluster_nuc_epe_en+ &
                           pc%z*(ua%Z-ua%zc)/dist
                   enddo

                else
                   do eq_a=1,unique_atoms(i)%N_equal_atoms
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                           -pc%position(:,eq_pcs))**2))
                      cluster_nuc_epe_en=cluster_nuc_epe_en+ &
                           pc%z*(unique_atoms(i)%Z-unique_atoms(i)%zc)/dist
                   enddo
             end if

             enddo
          enddo

       enddo ! N_unique_pcr
       ewpc_n=ewpc_n+pcs_n
    endif ! pcs_n.ne.0

    if(pcc_n.ne.0) then !** store pcc centers and calculate
                        ! cluster_nuc_epe_en for epe cores
       nullify(pc)
       do na=1,pcc_n
          pc=>ewpc_array(na+ewpc_n)
          pc%n_equal_charges=pcc_array(na)%n_equal_charges
          allocate(pc%position(3,pc%n_equal_charges),stat=ewa_allocstat(4))
          ASSERT(ewa_allocstat(4).eq.0)
          ewa_allocstat(4)=1
          pc=pcc_array(na)
          ewpc_arrel_used(na+ewpc_n)=.true.
          DPRINT 'pcc ewpc_arrel_used', na+ewpc_n

          do eq_pcc=1,pc%n_equal_charges
             do i=1,N_unique_atoms+n_timps
                if(i.gt.N_unique_atoms) then
                   ua=>unique_timps(i-N_unique_atoms)
                   do eq_a=1,ua%N_equal_atoms
                      dist=sqrt(sum((ua%position(:,eq_a) &
                           -pc%position(:,eq_pcc))**2))
                      cluster_nuc_epe_en=cluster_nuc_epe_en &
                           +pc%z*(ua%Z-ua%ZC)/dist
                   enddo
                else
                   do eq_a=1,unique_atoms(i)%N_equal_atoms
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                           -pc%position(:,eq_pcc))**2))
                      cluster_nuc_epe_en=cluster_nuc_epe_en &
                           +pc%z*(unique_atoms(i)%Z-unique_atoms(i)%ZC)/dist
                   enddo
                end if

             enddo
          enddo
       enddo ! N_unique_pcr
       ewpc_n=ewpc_n+pcc_n
    endif ! pcc_n.ne.0
        write(output_unit,*) 'final  number of ewald and epe centers '
        write(output_unit,*) 'to model external field actin on claster '
        write(output_unit,*) ewpc_n
        write(output_unit,*)
        write(output_unit,*) 'energy of interaction of the claster nuclei'
        write(output_unit,*) 'with epe only centers, cluster_nuc_epe_en'
        write(output_unit,*)  cluster_nuc_epe_en

       deallocate(pcew_temp,stat=ewa_allocstat(21))
        ASSERT(ewa_allocstat(21).eq.0)
       if(pcr_n.ne.0) then
           deallocate(pcr_temp,pcr_no,psb_ind, &
                          stat=ewa_allocstat(9))
           if(ewa_allocstat(9).ne.0)  &
                call error_handler("deallocate pcr_temp failed")
        end if


!       if(pcr_n.ne.0.and.pcs_n.eq.0.and.pcc_n.eq.0) then
!          do na=1,size(pcr_array)
!           deallocate(pcr_array(na)%position,stat=ewa_allocstat(24))
!           ASSERT(ewa_allocstat(24).eq.0)
!          enddo
!          deallocate(pcr_array,stat=ewa_allocstat(17))
!          ASSERT(ewa_allocstat(17).eq.0)
!        endif

!        if(pcs_n.ne.0) then
!          do na=1,size(pcs_array)
!           deallocate(pcs_array(na)%position,stat=ewa_allocstat(25))
!           ASSERT(ewa_allocstat(25).eq.0)
!          enddo
!          deallocate(pcs_array,stat=ewa_allocstat(18))
!          ASSERT(ewa_allocstat(18).eq.0)
!        endif

!        if(pcc_n.ne.0)  then
!         do na=1,size(pcc_array)
!          deallocate(pcc_array(na)%position,stat=ewa_allocstat(26))
!          ASSERT(ewa_allocstat(26).eq.0)
!         enddo
!         deallocate(pcc_array,stat=ewa_allocstat(19))
!         ASSERT(ewa_allocstat(19).eq.0)
!        endif

        if(allocated(pcew_no)) then
         deallocate(pcew_no,stat=ewa_allocstat(27))
         ASSERT(ewa_allocstat(27).eq.0)
        endif

  endif ! ex_pcew.and.ewaldpc
 end subroutine symm_ewpc_gen
#endif


    subroutine output_geometry (output_unit, unique_atoms)
      ! Purpose: write the geomtry to output in a user-friendly,
      !          script-friendly format. Special thanks to K.Albert!
      ! --------------------------------------------------
      use constants, only: angstrom
      use unique_atom_module, only: unique_atom_type
      implicit none
      integer, intent (in) :: output_unit
      type (unique_atom_type), intent (in) :: unique_atoms(:)
      ! *** end of interface ***

      integer(IK) :: i_eq, i_ua

      write (output_unit, *) "-- Geometry in au -----------------------------------------------------"

      do i_ua = 1, size (unique_atoms)
         do i_eq=1,unique_atoms(i_ua)%N_equal_atoms
            write(output_unit,1000)unique_atoms(i_ua)%Z,&
                 trim(unique_atoms(i_ua)%name),&
                 unique_atoms(i_ua)%position(1,i_eq),&
                 unique_atoms(i_ua)%position(2,i_eq),&
                 unique_atoms(i_ua)%position(3,i_eq)
         enddo
      enddo

      write (output_unit, *) " "
      write (output_unit, *) "-- Geometry in Angstroms ----------------------------------------------"
      do i_ua = 1, size (unique_atoms)
         do i_eq=1,unique_atoms(i_ua)%N_equal_atoms
            write(output_unit,1100)unique_atoms(i_ua)%Z,&
                 trim(unique_atoms(i_ua)%name),&
                 unique_atoms(i_ua)%position(1,i_eq) / angstrom,&
                 unique_atoms(i_ua)%position(2,i_eq) / angstrom,&
                 unique_atoms(i_ua)%position(3,i_eq) / angstrom
         enddo
      enddo

      write (output_unit, *) " "
      write (output_unit, *) "-- Geometry in xyz format ---------------------------------------------"
      write (output_unit, *) sum (unique_atoms % N_equal_atoms)
      write (output_unit, *) "# comment line"
      do i_ua = 1, size (unique_atoms)
         do i_eq = 1, unique_atoms(i_ua) % N_equal_atoms
            ! FIXME: the  names are supplied  in the input and  may be
            ! arbitrary.    They  may   not  be   understood   by  the
            ! viewer. Should we rather derive them from Z?
            write (output_unit, 1200) trim (unique_atoms(i_ua) % name), &
                 unique_atoms(i_ua) % position(1, i_eq) / angstrom, &
                 unique_atoms(i_ua) % position(2, i_eq) / angstrom, &
                 unique_atoms(i_ua) % position(3, i_eq) / angstrom
         enddo
      enddo
      write(output_unit,*)"-----------------------------------------------------------------------"

1000  format('au         ',F5.1,2X,A4,2X,3(F11.6,7X))
1100  format('angstrom   ',F5.1,2X,A4,2X,3(F11.6,7X))
1200  format(A4, 2X, 3(F11.6, 7X))
    end subroutine output_geometry

    subroutine output_geometry_pc()
      ! Purpose: write the geomtry to output in a user-friendly,
      !          script-friendly format. Special thanks to K.Albert!
      ! --------------------------------------------------
      use constants, only: angstrom
      implicit none
      ! *** end of interface ***

      integer(IK) :: i_eq, i_ua

      write(output_unit,*)"-- Geometry of Point Charges in au ------------------------------------"

      do i_ua=1,pointcharge_N
         do i_eq=1,pointcharge_array(i_ua)%N_equal_charges
            write(output_unit,1000)pointcharge_array(i_ua)%Z,&
                 trim(pointcharge_array(i_ua)%name),&
                 pointcharge_array(i_ua)%position(1,i_eq),&
                 pointcharge_array(i_ua)%position(2,i_eq),&
                 pointcharge_array(i_ua)%position(3,i_eq)
         enddo
      enddo

      write(output_unit,*)" "
      write(output_unit,*)"-- Geometry of Point Charges in Angstroms -----------------------------"
      do i_ua=1,pointcharge_N
         do i_eq=1,pointcharge_array(i_ua)%N_equal_charges
            write(output_unit,1100)pointcharge_array(i_ua)%Z,&
                 trim(pointcharge_array(i_ua)%name),&
                 pointcharge_array(i_ua)%position(1,i_eq) / angstrom,&
                 pointcharge_array(i_ua)%position(2,i_eq) / angstrom,&
                 pointcharge_array(i_ua)%position(3,i_eq) / angstrom
         enddo
      enddo
      write(output_unit,*)" "
      write(output_unit,*)"-----------------------------------------------------------------------"

1000  format('au         ',F5.1,2X,A4,2X,3(F11.6,7X))
1100  format('angstrom   ',F5.1,2X,A4,2X,3(F11.6,7X))
    end subroutine output_geometry_pc

    subroutine output_geometry_timp()
      ! Purpose: write the geomtry to output in a user-friendly,
      !          script-friendly format. Special thanks to K.Albert!
      ! --------------------------------------------------
      use constants, only: angstrom
      implicit none
      ! *** end of interface ***

      integer(IK) :: i_eq, i_ua

      write(output_unit,*)"-- Geometry of timps in au ---------------------------------------------"

      do i_ua=1,n_timps
         do i_eq=1,unique_timps(i_ua)%n_equal_atoms
            write(output_unit,1000) unique_timps(i_ua)%Z,&
                 trim(unique_timps(i_ua)%name),&
                 unique_timps(i_ua)%position(1,i_eq),&
                 unique_timps(i_ua)%position(2,i_eq),&
                 unique_timps(i_ua)%position(3,i_eq)
         enddo
      enddo

      write(output_unit,*)" "
      write(output_unit,*)"-- Geometry of timps Charges in Angstroms ------------------------------"
      do i_ua=1,n_timps
         do i_eq=1,unique_timps(i_ua)%N_equal_atoms
            write(output_unit,1100) unique_timps(i_ua)%Z,&
                 trim(unique_timps(i_ua)%name),&
                 unique_timps(i_ua)%position(1,i_eq) / angstrom,&
                 unique_timps(i_ua)%position(2,i_eq) / angstrom,&
                 unique_timps(i_ua)%position(3,i_eq) / angstrom
         enddo
      enddo
      write(output_unit,*)" "
      write(output_unit,*)"-----------------------------------------------------------------------"

1000  format('au         ',F5.1,2X,A4,2X,3(F11.6,7X))
1100  format('angstrom   ',F5.1,2X,A4,2X,3(F11.6,7X))
    end subroutine output_geometry_timp

  end subroutine symm_equivalents_gen

  !*************************************************************

  subroutine symm_mapping_gen( uas,                         symm_positions_map )
    !  Purpose: generates index mappings within all equal atom groups
    !------------ Modules used ------------------- ---------------
    use type_module
    use datatype          , only : arrmat2int
    use group_module      , only : group_num_el                                &
                                 , ylm_trafos
    use unique_atom_module, only : unique_atom_type
    implicit none
    type(unique_atom_type)       , intent(in)  :: uas(:)
    type(arrmat2int), allocatable, intent(out) :: symm_positions_map(:)
    !------------ Declaration of local constants  ----------------
    real(RK),parameter :: small = 1.e-10_rk
    ! very small value
    !------------ Declaration of local variables -----------------
    integer(IK)                   :: i, j, k, l
    integer(IK)                   :: N_equal_atoms
    real(RK)                      :: position1(3)                              &
                                   , position2(3)                              &
                                   , position3(3)
    !------------ Executable code --------------------------------
    allocate( symm_positions_map(size( uas )) )
    ! loop over all uniques
    do i = 1, size( uas )
      N_equal_atoms = uas(i)%N_equal_atoms
      allocate( symm_positions_map(i)%m(group_num_el, N_equal_atoms) )
      do j = 1, N_equal_atoms
        ! reorder coordinates of unique atom as
        ! (x,y,z) --> (x,z,y) in order to comply with the
        ! convention for angular momentum l=1
        position1(1) = uas(i)%position(1, j)
        position1(2) = uas(i)%position(3, j)
        position1(3) = uas(i)%position(2, j)
        !
        do k = 1, group_num_el
          ! do transformation by applying transformation matrix
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,k),position1)
          !
          do l = 1, N_equal_atoms
            ! reorder coordinates of unique atom as
            position3(1) = uas(i)%position(1, l) - position2(1)
            position3(2) = uas(i)%position(3, l) - position2(2)
            position3(3) = uas(i)%position(2, l) - position2(3)
            ! look if transformation target matches
            if ( dot_product( position3, position3 ) < small ) then
              ! mapping found. set index and leave
              symm_positions_map(i)%m(k, j) = l
              exit
            end if
          end do
        end do
      enddo
    enddo
  end subroutine symm_mapping_gen

  !*************************************************************

  subroutine symm_mapping_shells( uas, n_atm,                         symm_map )
    !  Purpose: generates index mappings within all equal atom groups
    !------------ Modules used -----------------------------------
    use type_module
    use datatype          , only : arrmat2int
    use group_module      , only : group_num_el                                &
                                 , ylm_trafos
    use unique_atom_module, only : unique_atom_type
    implicit none
    type(unique_atom_type)       , intent(in)  :: uas(:)
    integer(IK)                  , intent(in)  :: n_atm
    integer(IK)     , allocatable, intent(out) :: symm_map(:,:)
    !------------ Declaration of local constants  ----------------
    real(RK), parameter           :: small = 1.0e-11_rk ! very small value
    !------------ Declaration of local variables -----------------
    integer(IK)                   :: i, j, k, l, m
    integer(IK)                   :: n_equal_atoms
    real(RK)                      :: position1(3)                              &
                                   , position2(3)                              &
                                   , position3(3)
    logical                       :: found
    !------------ Executable code --------------------------------
    !
    ALLOCATE( symm_map(group_num_el, n_atm) )
    !
    ! loop over all uniques
    m = 0
    DO i = 1, size( uas )
      n_equal_atoms = uas(i)%N_equal_atoms
      !
      ! the l = 0 case (assumed to be always present)
      DO j = 1, n_equal_atoms
        m = m + 1
        ! reorder coordinates of unique atom as
        ! (x,y,z) --> (x,z,y) in order to comply with the
        ! convention for angular momentum l=1
        position1(1) = uas(i)%position(1, j)
        position1(2) = uas(i)%position(3, j)
        position1(3) = uas(i)%position(2, j)
        !
        DO k = 1, group_num_el
          !
          found = .FALSE.
          !
          ! do transformation by applying transformation matrix
          position2 = MATMUL( ylm_trafos(1)%matrix(:,:,k), position1 )
          !
          DO l = 1, n_equal_atoms
            !
            ! reorder coordinates of equal atom
            position3(1) = ABS( uas(i)%position(1, l) - position2(1) )
            position3(2) = ABS( uas(i)%position(3, l) - position2(2) )
            position3(3) = ABS( uas(i)%position(2, l) - position2(3) )
            !
            ! look where transformation target matches
            found = all( position3 < small )
            IF ( found ) THEN
              ! mapping found. set index and leave
              symm_map(k, m) = l - j
              EXIT
            END IF
          END DO
          IF ( .not. found ) stop 'unable to map u, e, sym'
        END DO
      END DO
    ENDDO
  end subroutine symm_mapping_shells

  !*************************************************************

  subroutine symm_nabor(io, unique_atoms, IMOD)
    use unique_atom_module, only: unique_atom_type
    use punchfile, only: pun_coordinates, pun_connectivity
    implicit none
    integer(IK)           , intent(in)  :: io
    type(unique_atom_type), intent(in)  :: unique_atoms(:)
    integer(IK)           , intent(in)  :: IMOD
    optional :: IMOD
    ! *** end of interface ***

    integer(IK)           :: IMODE
    integer(IK)           :: NA, ic,i,j
    integer(IK)           :: alloc_stat
    real(RK), allocatable :: coords(:,:)
    real(RK), allocatable :: atnum(:)
    integer(IK) , allocatable :: bonds(:,:)
    character(len=12), allocatable :: names(:)

    IMODE = 0
    if( present(IMOD) ) IMODE = IMOD

    NA = SUM( unique_atoms(:)%N_EQUAL_ATOMS )

    allocate(coords(3,NA),atnum(NA),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    ic=0
    do i=1,size(unique_atoms)
       do j=1,unique_atoms(i)%N_equal_atoms
          ic=ic+1
          coords(:,ic) = unique_atoms(i)%position(:,j)
          atnum(ic) = unique_atoms(i)%Z
       enddo
    enddo

    select case(IMODE)
    case (0)
       call NABOR(io,NA,coords,atnum)
    case (1)
       allocate(bonds(NA,NA),names(NA))

       ic = 0
       do i=1,size(unique_atoms)
          do j=1,unique_atoms(i)%N_equal_atoms
             ic = ic + 1
             names(ic) = unique_atoms(i)%name
          enddo
       enddo

       call NABOR(io, NA, coords, atnum, IMODE, bonds)

       call pun_coordinates(io,names,coords)
       call pun_connectivity(io,bonds)

       deallocate(bonds,names)
    end select

    deallocate(coords,atnum,STAT=alloc_stat)
    ASSERT(alloc_stat==0)
  end subroutine symm_nabor

  !*************************************************************
  subroutine NABOR(io,NA,CC,Z,IMOD,BMAT)
    !-----------------------------------------------------------------------
    !
    !     T.A.HALGREN, THE CITY COLLEGE OF THE CITY UNIVERSITY OF NEW YORK
    !     CALCULATE INTERATOMIC DISTANCES
    !     CALCULATE AND PRINT BOND ANGLES AND DIHEDRAL ANGLES FOR
    !     NEAREST NEIGHBOR ATOMS (DEFINED AS BEING CLOSER THAN ABOUT
    !     1.2 TIMES THE SUM OF THE NORMAL COVALENT RADII).
    !-----------------------------------------------------------------------
    use constants, only: angstrom
    implicit none
    integer(IK), intent(in)  :: io
    integer(IK), intent(in)  :: na
    real(RK)   , intent(in)  :: cc(:,:), z(:) ! (na,na), (na)
    integer(IK), intent(in)   :: IMOD
    integer(IK), intent(out)  :: BMAT(:,:) ! bond matrix (na,na)
    optional :: IMOD,BMAT
    ! *** end of interface ***

    ! --- declaration of local variables -----------------------------------
!   real(RK)    :: angle(4),covrad(100)
!   integer(IK) :: ma(4),mb(4),mc(4)
!   real(RK)    :: nbor(na,na),r(na,na)
    character(len=2)      :: ATOM(100)
    integer(IK) :: i,j,ikan,jkan,nam1,nijk,jp1,k, n
    real(RK)    :: raij,cov,rij,rik,rjk,theta,costh
    real(RK),parameter  :: scale_fac=1.2_rk
    integer(IK) :: NBOR(NA,NA),MA(4),MB(4),MC(4) !,MD(4)
    real(RK)    :: R(NA,NA),ANGLE(4)
    !         CC(3,NA),Z(NA)
    real(RK) :: COVRAD(100)
    !PK--
    data ATOM/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na', &
         'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr', &
         'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
         'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',&
         'Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',&
         'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re',&
         'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',&
         'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/
    data COVRAD /0.3200_rk,0.9300_rk,1.2300_rk,0.9000_rk,0.8200_rk,0.7700_rk,0.7500_rk, &
         .7300_rk,0.7200_rk,0.7100_rk,1.5400_rk,1.3600_rk,1.1800_rk,1.1100_rk,1.0600_rk,&
        1.0200_rk,0.9900_rk,0.9800_rk,2.0300_rk,1.7400_rk,1.4400_rk,1.3200_rk,1.2200_rk,&
        1.1800_rk,1.1700_rk,1.1700_rk,1.1600_rk,1.1500_rk,1.1700_rk,1.2500_rk,1.2600_rk,&
        1.2200_rk,1.2000_rk,1.1600_rk,1.1400_rk,1.1200_rk,2.1600_rk,1.9100_rk,1.6200_rk,&
        1.4500_rk,1.3400_rk,1.3000_rk,1.2700_rk,1.2500_rk,1.2500_rk,1.2800_rk,1.3400_rk,&
        1.4800_rk,1.4400_rk,1.4100_rk,1.4000_rk,1.3600_rk,1.3300_rk,1.3100_rk,2.3500_rk,&
        1.9800_rk,1.6900_rk,1.6500_rk,1.6500_rk,1.6400_rk,1.6300_rk,1.6200_rk,1.8500_rk,&
        1.6100_rk,1.5900_rk,1.5900_rk,1.5800_rk,1.5700_rk,1.5600_rk,1.5600_rk,1.5600_rk,&
        1.4400_rk,1.3400_rk,1.3000_rk,1.2800_rk,1.2600_rk,1.2700_rk,1.3000_rk,1.3400_rk,&
        1.4900_rk,1.4800_rk,1.4700_rk,1.4600_rk,1.4600_rk,1.4500_rk,1.0000_rk,1.0000_rk,&
        1.0000_rk,1.0000_rk,1.6500_rk,1.0000_rk,1.4200_rk,1.0000_rk,1.0000_rk,1.0000_rk,&
        1.0000_rk,1.0000_rk,1.0000_rk,1.0000_rk,0.8000_rk/

    INTEGER(IK) :: IMODE

    IMODE = 0
    IF( PRESENT(IMOD) ) IMODE = IMOD
    !PK--
    !PK-- ---------------------------------------------------------------
    !PK   FIRST BOND DISTANCES
    !PK-- ---------------------------------------------------------------

    SELECT CASE(IMODE)
    CASE (0)
       ! CASE 0: output distances and angles
       if (NA .LE. 2) RETURN
       write (io,1000)
    CASE (1)
       ! CASE 1: output connectivitiy matrix
       ASSERT(PRESENT(BMAT))
    END SELECT

    do 120 I = 1,NA
       IKAN = int(abs(Z(I)) + 0.5_rk)
       do 110 J = 1,I
          JKAN = int(abs(Z(J)) + 0.5_rk)
          NBOR(I,J) = 0
          RAIJ = sqrt((CC(1,I)-CC(1,J))**2+(CC(2,I)-CC(2,J))**2 &
               + (CC(3,I)-CC(3,J))**2 )
          R(I,J) = RAIJ
          R(J,I) = RAIJ
          COV = (COVRAD(IKAN) + COVRAD(JKAN)) * 1.2_rk * angstrom
          if (RAIJ .GT. scale_fac*COV)   GO TO 100
          NBOR(I,J) = 1
          if (I .EQ. J) GO TO 100
          RIJ = RAIJ / angstrom

          SELECT CASE(IMODE)
          CASE (0)
             write (io,2000) I,J,ATOM(IKAN),ATOM(JKAN),RAIJ,RIJ
             if ( RIJ < 0.3 ) then
                print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                print *, "XX                                                            XX"
                print *, "XX                    W A R N I N G !                         XX"
                print *, "XX                                                            XX"
                print *, "XX               Atoms too close! See below.                  XX"
                print *, "XX                                                            XX"
                print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                print *, ""
                print *, I, ATOM(IKAN), CC(1:3, I)
                print *, J, ATOM(JKAN), CC(1:3, J)
                print *, ""
                ABORT("Atoms to close, see tty!")
             endif
!         CASE (1)
!            BMAT(I,J) = 1
!            BMAT(J,I) = 1
          END SELECT

100       NBOR(J,I) = NBOR(I,J)
110    continue
120 continue

    ! output connectivity matrix if required:
    if(present(BMAT)) BMAT = NBOR

    IF ( IMODE > 0 ) RETURN
    !PK-- ---------------------------------------------------------------
    !PK   NOW BOND ANGLES
    !PK-- ---------------------------------------------------------------
    write (io,3000)
    NAM1 = NA-1
    NIJK = 0
    do 160 I = 1,NA
       do 150 J = 1,NAM1
          if (J .EQ. I) GO TO 150
          if (NBOR(I,J) .NE. 1) GO TO 150
          JP1 = J+1
          RIJ = R(I,J)
          do 140 K = JP1,NA
             if (K .EQ. I) GO TO 140
             if (NBOR(I,K) .NE. 1) GO TO 140
             RIK = R(I,K)
             RJK = R(J,K)
             !-----------------------------------------------------------------------
             !     COMPUTE J-I-K BOND ANGLE
             !-----------------------------------------------------------------------
             if (RIK*RIJ .EQ. 0.0_rk) GO TO 130
             COSTH = (RIK**2+RIJ**2-RJK**2)/(2.0_rk*RIK*RIJ)
             if (COSTH .LT. -1.0_rk) COSTH = -1.0_rk
             if (COSTH .GT. 1.0_rk) COSTH = 1.0_rk
             THETA =  acos(COSTH)*180.0_rk/3.141592653589793_rk
             NIJK = NIJK+1
             MA(NIJK) = I
             MB(NIJK) = J
             MC(NIJK) = K
             ANGLE(NIJK) = THETA
             if (NIJK .NE. 4) GO TO 130
             write (io,4000) (MB(N),MA(N),MC(N),ANGLE(N),N = 1,4)
             NIJK = 0
130       continue
140    continue
150 continue
160 continue
 if (NIJK .NE. 0) write (io,4000) (MB(N),MA(N),MC(N),ANGLE(N),N = 1, &
      NIJK)


 !PK-- ------------------------------------------------------------------
 !PK   NOW PROCEED WITH DIHEDRAL ANGLES
 !PK-- ------------------------------------------------------------------
 !PK--
 !     GO TO 297
 !PK--
 !     IF (NA.LT.4) GO TO 297
 !     write(io,190)
 ! 190 FORMAT(///48X,'         DIHEDRAL ANGLES (DEGREES)'//,
 !    .    7X,'A - B - C - D    ANGLE',8X,'A - B - C - D    ANGLE',
 !    .    8X,'A - B - C - D    ANGLE',8X,'A - B - C - D    ANGLE'/)
 !     NIJK=0
 !     DO 200 I=2,NA
 !     IM1=I-1
 !     DO 210 L=1,IM1
 !     DO 220 J=1,NA
 !     IF(J.EQ.I.OR.J.EQ.L) GO TO 220
 !     IF(NBOR(I,J).NE.1) GO TO 220
 !     DO 230 K=1,NA
 !     IF(K.EQ.I.OR.K.EQ.J.OR.K.EQ.L) GO TO 230
 !     IF(NBOR(J,K)*NBOR(K,L).NE.1) GO TO 230
 !     AX=((CC(2,K)-CC(2,J))*(CC(3,I)-CC(3,J))-(CC(2,I)-CC(2,J))*
 !    1 (CC(3,K)-CC(3,J)))
 !     BX=((CC(2,J)-CC(2,K))*(CC(3,L)-CC(3,K))-(CC(2,L)-CC(2,K))*
 !    1(CC(3,J)-CC(3,K)))
 !     AY=((CC(1,I)-CC(1,J))*(CC(3,K)-CC(3,J))-(CC(1,K)-CC(1,J))*
 !    1 (CC(3,I)-CC(3,J)))
 !     BY=((CC(1,L)-CC(1,K))*(CC(3,J)-CC(3,K))-(CC(1,J)-CC(1,K))*
 !    1 (CC(3,L)-CC(3,K)))
 !     AZ=((CC(1,K)-CC(1,J))*(CC(2,I)-CC(2,J))-(CC(1,I)-CC(1,J))*
 !    1 (CC(2,K)-CC(2,J)))
 !     BZ=((CC(1,J)-CC(1,K))*(CC(2,L)-CC(2,K))-(CC(1,L)-CC(1,K))*
 !    1 (CC(2,J)-CC(2,K)))
 !     ABSA= SQRT(AX**2+AY**2+AZ**2)
 !     ABSB= SQRT(BX**2+BY**2+BZ**2)
 !     IF(ABSA.LT.1.0D-6.OR.ABSB.LT.1.0D-6) GO TO 230
 !     COSTH=(AX*BX+AY*BY+AZ*BZ)/ABSA/ABSB
 !     IF(COSTH.LT.-1.) COSTH=-1.
 !     IF(COSTH.GT.1.) COSTH=1.
 !     THETA= ACOS(COSTH)*180.D0/3.14159 26535 89793 E0
 !     THETA = 180. - THETA
 !     NIJK=NIJK+1
 !     MA(NIJK)=I
 !     MB(NIJK)=J
 !     MC(NIJK)=K
 !     MD(NIJK)=L
 !     ANGLE(NIJK)=THETA
 !     IF(NIJK.NE.4) GO TO 250
 !     write(io,255) (MA(N),MB(N),MC(N),MD(N),ANGLE(N),N=1,4)
 ! 255 FORMAT(4(7X,I2,1H-,I2,I5,1H-,I2,F10.4))
 !     NIJK=0
 ! 250 CONTINUE
 ! 230 CONTINUE
 ! 220 CONTINUE
 ! 210 CONTINUE
 ! 200 CONTINUE
 !     IF(NIJK.NE.0) write(io,255) (MA(N),MB(N),MC(N),MD(N),ANGLE(N),
 !    1 N=1,NIJK)

1000 format(31X,'NEAREST NEIGHBOR INTERATOMIC DISTANCES',//,35X,'(ATOMIC UNITS)   *  (ANGSTROMS) ')
2000 format(18X,I3,' - ',I3,2X,A2,',',A2,2X,F10.6,8X,F10.6)
3000 format(47X,'NEAREST NEIGHBOR BOND ANGLES',//,12X,'ANGLE',7X,'DEGREES',12X,&
          'ANGLE',7X,'DEGREES',12X,'ANGLE',7X,'DEGREES',12X,'ANGLE',7X,'DEGREES')
!3000 format(10X,'ANGLE',7X,'DEGREES',11X,'ANGLE',7X,'DEGREES',11X,'ANGLE',7X,'DEGREES',11X,'ANGLE',7X,'DEGREES')
4000 format(4(7X,I3,' -',I3,' -',I3,F11.4))
end subroutine

  !*************************************************************


  !--------------- End of module ----------------------------------
end module symm_positions

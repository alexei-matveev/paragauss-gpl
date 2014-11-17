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
  subroutine symm_epe
! treats EPE centers
#include <def.h>
  use datatype
  use calc3c_switches
  use iounitadmin_module
  use ewaldpc_module,ex_pcr_dummy=>ex_pcr
  use filename_module
  use group_module
  use symm_positions
  use unique_atom_module
  use pointcharge_module
  use symm_positions
  use type_module, only:IK => i4_kind, RK => r8_kind
  use comm, only: comm_rank

implicit none

    type(symm_transformation_int)        :: point_trafos_pc
    ! transformations of the symmetry equivalent point charges

    type(sub_group)                      :: local_groups_pc
    ! local symmetry group of the unique point charges
    type(group_coset)                    :: cosets_pc
    ! coset of the local symmetry group of the unique point charges

  real(RK)                   :: position(3),position2(3)
  type(pointcharge_type), pointer      :: pc,pcr
        logical:: ex_pcr
        real (RK):: e_nuc_pcr,dist
        integer(IK):: na,nb,eq_a,pcr_unit,n_gxat,kl
  integer(IK)::  it_pcr,ii
  integer(IK):: i,n,j,n_equal
  real(RK), allocatable, dimension(:,:) ::pcs_temp
  real(RK), allocatable, dimension(:,:) ::pcc_temp
  real(RK),parameter:: small_pcr=5.d-1
  real(RK):: min_distepe


!  this does not work on Hitachi
!  inquire (file= trim(inpfile('epe.pcr')), exist=ex_pcr)
!  print*,  ' epe.pcr has been found',ex_pcr

!  pcr_ex: if(ex_pcr.and.epe_embedding) then
  pcr_ex: if(epe_embedding) then
   ex_pcr=.true.
   ex_pcs=.true.
   ex_pcc=.true.

  if(ex_pcs) then
   pcr_unit=openget_iounit(file=trim(inpfile('epe.pcs')), &
                               form='FORMATTED',status='old')
   read(pcr_unit,*)  pcs_N      !module variable
   if(output_unit >  0) write(output_unit,*)  ' Number of EPE shell pcs_N ',pcs_N
   allocate(pcs_temp(4,pcs_N),stat=ewa_allocstat(10)) 
   ASSERT(ewa_allocstat(10).eq.0)
   ewa_allocstat(10)=1

   do i=1,pcs_N
    read(pcr_unit,*) pcs_temp(:,i)
   enddo ! i=1,pcs_N
    call returnclose_iounit(pcr_unit)
  endif ! ex_pcs

  if(ex_pcc) then
   pcr_unit=openget_iounit(trim(inpfile('epe.pcc')), &
                          form='FORMATTED',status='old')
   read(pcr_unit,*)  pcc_N      !module variable
   if(output_unit >  0) write(output_unit,*) ' Number of EPE cores pcc_N ',pcc_N
   allocate(pcc_temp(4,pcc_N),stat=ewa_allocstat(11))
   ASSERT(ewa_allocstat(11).eq.0)
          ewa_allocstat(11)=1
        do i=1,pcc_N
         read(pcr_unit,*) pcc_temp(:,i)
        enddo ! i=1,pcc_N
    call returnclose_iounit(pcr_unit)
    if(pcc_N.ne.pcs_N) stop 'pcs_n ne pcc_n'
  endif ! ex_pcc

   nullify(pc)
   !first it will be used as ia temp and then as element of array
   
   pcr_unit=openget_iounit(trim(inpfile('epe.pcr')), &
                          form='FORMATTED',status='old')
   read(pcr_unit,*)  pcr_N  !module variable 
   DPRINT pcr_N ,'pcr_N read from epe.pcr'

!   allocate(pc,pcr_temp(4,pcr_N),pcr_no(pcr_N),psb_ind(pcr_N),stat=ewa_allocstat(9)) !4
   allocate(pc,pcr_temp(4,pcr_N),pcr_no(pcr_N),stat=ewa_allocstat(9)) !4
   if(ewa_allocstat(9).ne.0) call error_handler("allocate pcr_temp failed")
      ewa_allocstat(9)=1 ! pcr_temp pcr_no psb_ind
      ewa_allocstat(16)=1 ! pc

  
    ! loop over all centers in file and sort out PC which positions coincide
    ! with the positions of regular atoms
    if(output_unit >  0)  write(output_unit,*) 'Number of EPE reference points pcr_n ' ,pcr_n

    do i=1,pcr_n
!        read(pcr_unit,*)  pcr_temp(:,i),psb_ind(i)
        read(pcr_unit,*)  pcr_temp(:,i)
    enddo ! i=1,pcr_n
!       n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)
    ASSERT(allocated(gxepe_impu))
     n_gxat=0
     do i=1,N_unique_atoms
      if(gxepe_impu(i).ne.0) n_gxat=n_gxat+unique_atoms(i)%N_equal_atoms
     enddo
     if(output_unit >  0)    write(output_unit,*) 'number of atoms in gx-file', n_gxat
        do i=1,n_timps
         if(gxepe_impu(i).ne.0) n_gxat=n_gxat+unique_timps(i)%n_equal_atoms
        enddo
     if(output_unit >  0)  write(output_unit,*)'number of atoms and n_timps in gx file ',n_gxat

        if(.not.ex_gxepe) call error_handler(" file epe.r not found")

           kl=0
#if 0 /* no epe centers which would coinside with atoms of cluster */
        md:do i=1,pcr_n
           do i_ua=1,N_unique_atoms+n_timps
              if(gxepe_impu(i_ua).ne.0) then
                 if(i_ua.gt.N_unique_atoms) then
                    ua=>unique_timps(i_ua-N_unique_atoms)
                 else
                    ua=>unique_atoms(i_ua)
                 end if
           
                 do i_eq=1,ua%N_equal_atoms
                    if (dot_product( &
                         gxepe_array(i_ua)%position(:,i_eq)-pcr_temp(1:3,i), &
                         gxepe_array(i_ua)%position(:,i_eq)-pcr_temp(1:3,i)) &
                                      .lt.small_pcr) then
                       kl=kl+1
                       ! swap pcr positions
                       pc%position_first_ec(:)=pcr_temp(1:3,kl)
                       pc%z=pcr_temp(4,kl)
                       pc%c=0.0_rk
                       pcr_temp(1:3,kl)=pcr_temp(1:3,i)
                       pcr_temp(4,kl)=pcr_temp(4,i)
                       pcr_temp(1:3,i)=pc%position_first_ec
                       pcr_temp(4,i)=pc%z
#if 0
                       index=psb_ind(kl)
                       psb_ind(kl)=psb_ind(i)
                       psb_ind(i)=index
#endif
                       
                       pc%position_first_ec(:)=pcs_temp(1:3,kl)
                       pc%z=pcs_temp(4,kl)
                       pc%c=0.0_rk
                       pcs_temp(1:3,kl)=pcs_temp(1:3,i)
                       pcs_temp(4,kl)=pcs_temp(4,i)
                       pcs_temp(1:3,i)=pc%position_first_ec
                       pcs_temp(4,i)=pc%z
                       
                       pc%position_first_ec(:)=pcc_temp(1:3,kl)
                       pc%z=pcc_temp(4,kl)
                       pc%c=0.0_rk
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
#endif

        n_gxat_pcr=kl

        call returnclose_iounit(pcr_unit)

        if(output_unit >  0) write(output_unit,*) &
             'pcr_n after sorting out PC in atomic positions' ,pcr_n-kl

#if 1 /*  check potential of PC */
        pc%z=0.0_rk
        e_nuc_pcr = 0.0_rk
        do nb=kl+1,pcr_n
           pc%z=pc%z+pcs_temp(4,nb)+pcc_temp(4,nb)
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 e_nuc_pcr = e_nuc_pcr+pcs_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a)-pcs_temp(1:3,nb))**2))
                 e_nuc_pcr = e_nuc_pcr+pcc_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a)-pcc_temp(1:3,nb))**2))
              enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo ! na=1,N_unique_atoms
        enddo ! 1,pcr_n

        if(output_unit >  0) then
           write(output_unit,*) ' Energy of interaction of the EPE centers and '
           write(output_unit,*) ' atoms of cluster e_nuc_pcr calculsted  with use '
           write(output_unit,*) ' use of pcr_temp coordinates and  total charge'
           write(output_unit,*) ' of the EPE centers Z', e_nuc_pcr,pc%z
           DPRINT 'e_nuc_pcr and sum Z', e_nuc_pcr, pc%z
           DPRINT 'to be compared with corresp part of ewpc_array'
        endif
#endif

#if 1 /*  now oder pcr_temp on groups of symm equivalent centers */
        i=1+kl
        it_pcr=0
    pcr_wh:    do while(i.le.pcr_n)
           pc%position_first_ec(1:3)=pcr_temp(1:3,i)
           pc%name='EPE PC '
           it_pcr=it_pcr+1
           
           call setup_reorder()
           
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
                    pc%c=0.0_rk
                    pcr_temp(4,ii)=pcr_temp(4,i)
                    pcr_temp(4,i)=pc%z
                    pcr_no(i)=it_pcr
#if 0
                    index=psb_ind(ii)
                    psb_ind(ii)=psb_ind(i)
                    psb_ind(i)=index
#endif
                    
                    pc%position(1:3,j)=pcs_temp(1:3,ii)
                    if(j.eq.1) pc%position_first_ec=pc%position(1:3,j)
                    pc%z=pcs_temp(4,ii)
                    pc%c=0.0_rk
                    pcs_temp(1:3,ii)=pcs_temp(1:3,i)
                    pcs_temp(4,ii)=pcs_temp(4,i)
                    pcs_temp(1:3,i)=pc%position(1:3,j)
                    pcs_temp(4,i)=pc%z
                    
                    pc%position(1:3,j)=pcc_temp(1:3,ii)
                    if(j.eq.1) pc%position_first_ec=pc%position(1:3,j)
                    pc%z=pcc_temp(4,ii)
                    pc%c=0.0_rk
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

   call close_reorder()

    enddo  pcr_wh
#if 0
    call write_pcr_temp()
#endif
#endif

#if 1 /* check e_nuc_pcr */
    e_nuc_pcr=0.0 
    min_distepe=100.0
    do nb=1,pcr_n
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 dist=sqrt(sum((unique_atoms(na)%position(:,eq_a) -pcs_temp(1:3,nb))**2))
                 min_distepe=min(min_distepe,dist)
                 if(min_distepe.lt.1.7) then
                    if(comm_rank() == 0) then
                       print*,'min_distepe too small unique_atom:equal_atom',min_distepe,na,eq_a
                       print*,'pcs position',pcs_temp(1:3,nb)
                    endif
                    stop 'error in epe.pcs position'
                 endif
                 e_nuc_pcr = e_nuc_pcr+ &
                             pcs_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/dist 
                 dist=sqrt(sum((unique_atoms(na)%position(:,eq_a)-pcc_temp(1:3,nb))**2))
                 min_distepe=min(min_distepe,dist)
                 if(min_distepe.lt.1.7) then
                    if(comm_rank() == 0) then
                       print*,'min_distepe too small unique_atom:equal_atom',min_distepe,na,eq_a
                       print*,'pcc position',pcc_temp(1:3,nb)
                    endif
                    stop 'error in epe.pcc position'
                 endif
                 e_nuc_pcr = e_nuc_pcr+ &
                             pcc_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/dist
              enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo ! na=1,N_unique_atoms
    enddo
    DPRINT 'e_nuc_pcr min_distepe after symmetry reordering',e_nuc_pcr,min_distepe
#endif

#if 1 /* reorder shells */
    i=1+kl
    it_pcr=0  !types of EWPC
    pcs_wh: do while(i.le.pcr_n)
       pc%position_first_ec(1:3)=pcs_temp(1:3,i)
       pc%name='EPE PCS '
       it_pcr=it_pcr+1
       
       call setup_reorder()

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
  pc%c=0.0_rk
  pcs_temp(4,ii)=pcs_temp(4,i)
  pcs_temp(4,i)=pc%z
  exit
  endif ! 
 enddo ! ii=i+1,pcr_n
  if(i>pcr_n) stop ' i>pcr_n, check PCR array'
  i=i+1
enddo ! j=1,n_equal
 
  call close_reorder()

enddo pcs_wh ! while 
#endif

#if 1 /* reorder cores */
     i=1+kl
     it_pcr=0  !types of EWPC
    pcc_wh: do while(i.le.pcr_n)
       pc%position_first_ec(1:3)=pcc_temp(1:3,i)
       pc%name='EPE PCC '
       it_pcr=it_pcr+1

       call setup_reorder()

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
 pc%c=0.0_rk
 pcc_temp(4,ii)=pcc_temp(4,i)
 pcc_temp(4,i)=pc%z
 exit
 endif ! 
enddo ! ii=i+1,pcr_n
if(i>pcr_n) stop ' i>pcr_n, check PCR array'
i=i+1
enddo ! j=1,n_equal

  call close_reorder()

 enddo pcc_wh ! while 
#endif

        deallocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)
        nullify(pc)

#if 1 /* check e_nuc_pcr */
      ! check potential of PC
        allocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)
        ewa_allocstat(16)=1

        e_nuc_pcr = 0.0_rk
        pc%z=0.0_rk
        do nb=1+kl,pcr_n

        if(pcs_n.ne.0) then
           pc%z=pc%z+pcs_temp(4,nb)+pcc_temp(4,nb)
             do na=1,N_unique_atoms
                   do eq_a=1,unique_atoms(na)%N_equal_atoms
                      e_nuc_pcr = e_nuc_pcr+pcs_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                           sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                                                 -pcs_temp(1:3,nb))**2))
                      e_nuc_pcr = e_nuc_pcr+pcc_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
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
     if(output_unit >  0) then 
        write(output_unit,*) 'e_nuc_pcr and Z after enforced symm',e_nuc_pcr,pc%z
        DPRINT 'e_nuc_pcr and Z after enforced symm',e_nuc_pcr,pc%z
     endif

        ! use pc as a tem finished here
        deallocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)
        nullify(pc)
#endif

     if(pcs_n.eq.0.and.pcc_n.eq.0) then
        allocate(pcr_array(it_pcr),stat=ewa_allocstat(17))
        ASSERT(ewa_allocstat(17).eq.0) 
               ewa_allocstat(17)=1
        
     else
#if 1 /* no pcr_array */
        allocate(pcr_array(it_pcr),stat=ewa_allocstat(17))
        ASSERT(ewa_allocstat(17).eq.0) 
               ewa_allocstat(17)=1
#endif

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

#if 1 /* no z_pcs_pcc */
     z_pcs_pcc: if(pcs_n.eq.0.or.pcc_n.eq.0) then
#if 0 /* not basic option */
        i=1+kl
        do while(i.le.pcr_n)
           pc=>pcr_array(pcr_no(i))
           n_equal=0
           pc%z=0.0_rk
           pc%c=0.0_rk

           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcr_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n

           pc%z=pc%z/n_equal

           allocate(pcr_array(pcr_no(i))%position(3,n_equal),stat=ewa_allocstat(24))
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
#endif

     else z_pcs_pcc! i.e. basic option
        
#if 1 /* no pcr_array */
        if(comm_rank() == 0) print*, '** treat regular positions'
        i=1+kl
        nullify(pc)
        do while(i.le.pcr_n)
           pc=>pcr_array(pcr_no(i))
           n_equal=0
           pc%z=0.0_rk
           pc%c=0.0_rk

           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcr_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n

           pc%z=pc%z/n_equal

           allocate(pc%position(3,n_equal),stat=ewa_allocstat(24))
           ASSERT(ewa_allocstat(24).eq.0)
           ewa_allocstat(24)=1
           pc%name='pcr    ' 
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcr_temp(1:3,i+ii-1)
!!!     print all generated charges
!        print '(4f15.8,i3,4i2)',pc%position(:,ii),pc%z,n_equal,0,0,1,0
!!$              write(output_unit, '(4f15.8,i3,4i2)') pc%position(:,ii), &
!!$                   pc%z,n_equal,0,0,1,0
           enddo ! ii=1,equal
!!!           dist=sqrt(sum((unique_atoms(1)%position(:,1) &
!!$                -pc%position(:,1))**2))
!!$           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)' ) &
!!$                pc%position(:,1),pc%z,n_equal,0,0,1,0,dist
           i=i+n_equal
        enddo !while
#endif

        i=1+kl
        nullify(pc)
        if(comm_rank() == 0) then
           if(print_epe) print*,'atoms of pcs_array'
        endif
        fill_pcs:do while(i.le.pcr_n)
           pc=>pcs_array(pcr_no(i))
#if 1 /* no pcr_array */
           pcr=>pcr_array(pcr_no(i))
#endif

           n_equal=0
           pc%z=0.0_rk
           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcs_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n
           
           pc%z=pc%z/n_equal
           pc%c=0.0_rk

           allocate(pcs_array(pcr_no(i))%position(3,n_equal),stat=ewa_allocstat(25))
           ASSERT(ewa_allocstat(25).eq.0) 
                  ewa_allocstat(25)=1 ! pcs_array(:)%position
           
           pc%name='pcs    ' 
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcs_temp(1:3,i+ii-1)
!!!     print all generated charges
#if 1
        dist=sqrt(sum((pcr%position(:,ii)-pc%position(:,ii))**2))
#endif
#if 0 /* print all pcs */
        if(print_epe) print '(4f15.8,2i4,f15.8)',pc%position(:,ii),pc%z,n_equal,i+ii-1,dist
#endif
!!$              write(output_unit, '(4f15.8,i3,4i2)')&
!!$                   pc%position(:,ii),pc%z,n_equal,0,0,1,0

           enddo ! ii=1,equal
           pc%position_first_ec=pc%position(:,1)
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)'  &
#if 0
        if(print_epe) print '(4f15.8,i4,4i4,f15.8)', &
                 pc%position(:,1),pc%z,n_equal,0,0,1,0, &
                  sqrt(sum((unique_atoms(1)%position(:,1) &
                                    -pc%position(:,1))**2))
#endif
           i=i+n_equal
        enddo fill_pcs
        !** done

        nullify(pc)
        i=1+kl
        if(comm_rank() == 0) then
           if(print_epe) print*,'atoms of pcc_array'
        endif
        fill_pcc: do while(i.le.pcr_n)
           pc=>pcc_array(pcr_no(i))
#if 1
           pcr=>pcs_array(pcr_no(i))
#endif
           n_equal=0
           pc%z=0.0_rk
           
           do ii=i,pcr_n
              n_equal=n_equal+1
              pc%z=pc%z+pcc_temp(4,ii)
              if(ii.eq.pcr_n) exit
              if(pcr_no(ii+1).ne.pcr_no(ii)) exit
           enddo !ii=i,pcr_n
           
           pc%z=pc%z/n_equal
           pc%c=0.0_rk
           
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(26))
                    ASSERT(ewa_allocstat(26).eq.0)
                           ewa_allocstat(26)=1
           pc%name='pcc    ' 
           pc%N_equal_charges= n_equal
           do ii=1,n_equal
              pc%position(:,ii)=pcc_temp(1:3,i+ii-1)
!!!           print all generated charges
#if 1
        dist=sqrt(sum((pcr%position(:,ii)-pc%position(:,ii))**2))
#endif
#if 0 /* print all pcc */
           if(print_epe) &
             print '(4f15.8,2i5,f15.8)',pc%position(:,ii),pc%z,n_equal,i+ii-1,dist
#endif
           enddo ! ii=1,equal
           pc%position_first_ec=pc%position(:,1)
           !       print only first atoms in the group
!           write(output_unit, '(4f15.8,i3,4i2,f15.8)' ) &
!     if(print_epe) print '(4f15.8,i4,4i4,f15.8)', &
!                 pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!                 sqrt(sum((unique_atoms(1)%position(:,1) &
!                                   -pc%position(:,1))**2))
           i=i+n_equal
        enddo fill_pcc !while

#if 1 /* e_nuc_pcr after storing */
    e_nuc_pcr=0.0
    do nb=1,it_pcr
       pc=>pcc_array(nb)
       do ii=1,pc%N_equal_charges
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 e_nuc_pcr = e_nuc_pcr+pc%z*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                                            -pc%position(:,ii))**2))
                 e_nuc_pcr = e_nuc_pcr+pcs_array(nb)%z*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                      sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                                            -pcs_array(nb)%position(:,ii))**2))
              enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo ! na=1,N_unique_atoms
     enddo
!     DPRINT nb,e_nuc_pcr,pc%N_equal_charges,pcs_array(nb)%z,sum(pcs_array(nb)%position(:,pc%N_equal_charges))
     
    enddo
    if(comm_rank() == 0) then
       if(print_epe) print*, 'e_nuc_pcr after storing  ',e_nuc_pcr
    endif
#endif
        
        
    endif z_pcs_pcc!    (pcs_n.eq.0.or.pcc_n.eq.0) then
#endif

    if(output_unit >  0) then 
       write(output_unit,*) 'number of groups of symmetry equavalent PC', it_pcr
       print*,'number of groups of symmetry equavalent PC', it_pcr
    endif
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
#if 0
     call write_pc_array(pcc_array,'pcc_array')
     call write_pc_array(pcs_array,'pcs_array')
#endif
     endif pcr_ex ! ex_pcr.and.epe_embedding
 contains

#if 0
subroutine write_pcr_temp()
integer(i4_kind):: pc_unit,i
pc_unit=openget_iounit(file=trim(inpfile('pcr_temp')), form='FORMATTED',status='unknown')
write(pc_unit,*) size(pcr_temp,2)
do i=1,size(pcr_temp,2)
write(pc_unit,*) pcr_temp(:,i),pcr_no(i)
enddo
call returnclose_iounit(pc_unit,'keep')
end subroutine write_pcr_temp

subroutine write_pc_array(pc_array,name)
type(pointcharge_type):: pc_array(:)
character(len=9), intent(in):: name
integer(i4_kind):: pc_unit,i,eq
  pc_unit=openget_iounit(file=trim(inpfile(name)), form='FORMATTED',status='unknown')
write(pc_unit,*) size(pc_array)
do i=1,size(pc_array)
write(pc_unit,*) pc_array(i)%z,pc_array(i)%n_equal_charges,pc_array(i)%position_first_ec,pc_array(i)%name
do eq=1,pc_array(i)%n_equal_charges
write(pc_unit,*) pc_array(i)%position(:,eq)
enddo
enddo
call returnclose_iounit(pc_unit,'keep')
print *,name//' written '//trim(inpfile(name))

end subroutine write_pc_array
#endif

subroutine setup_reorder()
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

           call group_coset_decomp(n_equal,local_groups_pc, cosets_pc,point_trafos_pc%matrix)

                ewa_allocstat(13)=1 ! point_trafos_pc%matrix
                ewa_allocstat(14)=1 ! cosets_pc%elements

           pc%N_equal_charges=n_equal

           ! allocate positions of equal charges
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(15))
           ASSERT(ewa_allocstat(15).eq.0)
                  ewa_allocstat(15)=1
 end subroutine setup_reorder

 subroutine close_reorder()
           deallocate(local_groups_pc%elements, point_trafos_pc%matrix, &
                      cosets_pc%elements, pc%position, &
                      stat=ewa_allocstat(12))      ! local_groups_pc%elements
                ASSERT(ewa_allocstat(12).eq.0) 
                ewa_allocstat(13)=0 ! point_trafos_pc%matrix
                ewa_allocstat(14)=0 ! cosets_pc%elements
                ewa_allocstat(15)=0 ! pc%position
 end subroutine close_reorder

 end subroutine symm_epe

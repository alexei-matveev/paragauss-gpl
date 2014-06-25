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
subroutine symm_ewpc_gen ()
  !
  ! Purpose:  ewpc_array containing  all information  about external
  ! field is filled in
  !
# include "def.h"
  use datatype
  use calc3c_switches
  use filename_module
  use symm_positions
  use pointcharge_module
  use operations_module
  use iounitadmin_module
  use group_module, only: group_num_el, symm_transformation_int, &
       sub_group, group_coset, ylm_trafos, group_coset_decomp
  use unique_atom_module
  use ewaldpc_module, ex_pcew_dumm=>ex_pcew
  use type_module, only: IK => i4_kind, RK => r8_kind
  implicit none
  ! *** end of interface ***

    type(symm_transformation_int)        :: point_trafos_pc
    ! transformations of the symmetry equivalent point charges

    type(sub_group)                      :: local_groups_pc
    ! local symmetry group of the unique point charges
    type(group_coset)                    :: cosets_pc
    ! coset of the local symmetry group of the unique point charges

  real (RK) :: position(3), position2(3)
  type(pointcharge_type), pointer      :: pc

! treats PC of ewald program generated file
  type(unique_atom_type), pointer  :: ua
  logical:: ex_pcew
  real (RK):: e_nuc_ewpc,dist
  integer(IK):: na, nb, eq_a, kl_ew, kl1, kl_pcr, n_gxat, kl_ew_eq
  integer(IK)::  it_ewpc,it_pcr,ii,i_ua,i_eq,eq_pcc,eq_pcs
  integer(IK),dimension(:),allocatable:: pcew_no
  integer (IK):: n_ewpc_arr, i, j, n_equal, n, j_epe_r, epe_r(300)
  real(RK), allocatable, dimension(:,:) ::pcew_temp
  real(RK),parameter:: small_ew=1.d-2
  real(RK)::ewpc_array_charge,e_nuc_pcr,min_pccdist,min_pcsdist,min_pcrdist
  real (RK) :: min_ewpcdist, gv_nuc_ewpc(3)
  real(RK)::max_ewpc_pcr_dev=0.0_rk
  character(len=256) :: trace_message
  integer (IK) :: k
#if 0
  real(RK) :: r_screep, rc_screep(3), distance_to_screep
#endif

  print*, 'symm_ewpc_gen ______________________________________________'

  kl_pcr=0

  if(ewaldpc) then
  ex_pcew=.true.

   ewpc_unit=openget_iounit(trim(inpfile('ewald.pcr')) &
                           ,form='FORMATTED',status='old')
#if 0 /* control that cluster is inside screep */
   read(ewpc_unit,*)  EWPC_N,r_screep,rc_screep
       ii=0
       distance_to_screep=100.0_rk
       do i=1,N_unique_atoms+n_timps
       do eq_a=1,unique_atoms(i)%N_equal_atoms
       ii=ii+1
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                                           -rc_screep)**2))
       distance_to_screep=min(r_screep-dist,distance_to_screep)
                   enddo
             enddo
print*, 'minimum distance_to_screep' ,distance_to_screep
ASSERT(distance_to_screep.gt.2.0_rk)
#else
read(ewpc_unit,*)  EWPC_N
#endif

   nullify(pc) !called after symm_epe and thus ponter was  on  pcc_array pcs_array
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
    write(output_unit,*) 'as file epe.r with regular atomic positions is located'
    write(output_unit,*) ' EWPC centers coinciding with atoms will be sorted out'
!!$      write(output_unit,*)  gxepe_array(1)%position(:,1)
   end if

   epe_r=0
   ii=1
   it_ewpc=EWPC_N ! it_ewpc here is help variable to define number of PC
   mc:    do i=1,EWPC_N
      read(ewpc_unit,*)  pcew_temp(:,ii)
!!$      if(ii.eq.1) write(output_unit,*) pcew_temp(:,ii)
      j_epe_r=0
      do i_ua=1,N_unique_atoms+n_timps
         if(gxepe_impu(i_ua)==0) cycle  !!!!!!!!!!!!!!AS important
         if(ex_gxepe) then
            if(i_ua.gt.N_unique_atoms) then
               do i_eq=1,unique_timps(i_ua-N_unique_atoms)%N_equal_atoms
                  j_epe_r=j_epe_r+1
                  if (dot_product(gxepe_array(i_ua)%position(:,i_eq) &
                                                 -pcew_temp(1:3,ii), &
                       gxepe_array(i_ua)%position(:,i_eq)-pcew_temp(1:3,ii)) &
                       .lt.small_ew) then
                     epe_r(j_epe_r)=1
                     it_ewpc=it_ewpc-1
                     cycle  mc
                  endif ! dot_product
               end do
            else
               do i_eq=1,unique_atoms(i_ua)%N_equal_atoms
                  j_epe_r=j_epe_r+1
                  if (dot_product(gxepe_array(i_ua)%position(:,i_eq)- &
                       pcew_temp(1:3,ii), &
                       gxepe_array(i_ua)%position(:,i_eq)-pcew_temp(1:3,ii)) &
                       .lt.small_ew) then
                     epe_r(j_epe_r)=1
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
        write(output_unit,*) 'EWPC_N after sorting out PC in atomic positions' ,EWPC_N
        do j=1,n_epe_r
           if(epe_r(j)==0) print*,'WARNING : ',j,' EPE.R center has no EWALD center coincided !!!!!!!!!!!!'
        end do

       gx_epe_ex: if(ex_gxepe) then
         do i_ua=1,N_unique_atoms
           deallocate(gxepe_array(i_ua)%position,stat=ewa_allocstat(23))
           ASSERT(ewa_allocstat(23).eq.0)
         end do
           deallocate(gxepe_array,gxepe_impu,stat=ewa_allocstat(22))
           ASSERT(ewa_allocstat(23).eq.0)
         endif gx_epe_ex
!** done
#if 1
        e_nuc_ewpc=0.0_rk
        do na=1,N_unique_atoms
        do eq_a=1,unique_atoms(na)%N_equal_atoms
        gv_nuc_ewpc=0.0_rk
        do nb=1,EWPC_N
        dist=sqrt(sum((unique_atoms(na)%position(:,eq_a)  &
                                            -pcew_temp(1:3,nb))**2))
        e_nuc_ewpc=e_nuc_ewpc+pcew_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/dist
        gv_nuc_ewpc=gv_nuc_ewpc+ pcew_temp(4,nb)*&
                       (unique_atoms(na)%position(:,eq_a)-pcew_temp(1:3,nb))/dist**3
          enddo
          print*,'gv_nuc_ewpc: ',gv_nuc_ewpc,sum(gv_nuc_ewpc**2)
         enddo
        enddo
    print*,'e_nuc_ewpc with initial pcew_temp',e_nuc_ewpc
#endif

!**   reorder atoms in  pcr_temp and pcew_temp
!**  first go coinsiding atoms with the same charges
#if 0
call read_pcr_temp(pcr_n)
#endif


        kl_ew=0
        kl_ew_eq=0
        ex_pcr_n: if(pcr_n.ne.0) then  ! if epe charges exist
           if(n_timps.gt.0) then
              n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)+&
                     sum(unique_timps(1:n_timps)%N_equal_atoms)
           else
              n_gxat=sum(unique_atoms(1:N_unique_atoms)%N_equal_atoms)
           end if

           if(print_epe) print*,'n_gxat pcr_n n_gxat_pcr ',n_gxat, pcr_n, n_gxat_pcr
           kl_pcr=n_gxat_pcr
           write(output_unit,*) 'initial values  of kl_pcr and kl_ew',kl_pcr,kl_ew
           print*, 'initial values  of kl_pcr and kl_ew',kl_pcr,kl_ew
           do i=1,EWPC_N
              kl1=kl_pcr+1
              do k=kl1,pcr_n
                 if(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i)).lt.small_ew &
                      .and. ((abs(pcr_temp(4,k)-pcew_temp(4,i)).lt.0.01_rk) ) ) then
                    !**** centers coincide and charges are equal
                    ! such centers generated with ewald can be neglected
                    max_ewpc_pcr_dev=max(sqrt(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i))),max_ewpc_pcr_dev)
                    kl_pcr=kl_pcr+1
                    kl_ew=kl_ew+1
                    pc%position_first_ec=pcr_temp(1:3,kl_pcr)
                    pcr_temp(1:3,kl_pcr)=pcr_temp(1:3,k)
                    pcr_temp(1:3,k)=pc%position_first_ec
                    pc%z=pcr_temp(4,kl_pcr)
                    pcr_temp(4,kl_pcr)=pcr_temp(4,k)
                    pcr_temp(4,k)=pc%z
#if 0
                    index=psb_ind(kl_pcr)
                    psb_ind(kl_pcr)=psb_ind(k)
                    psb_ind(k)=index
#endif

                    pc%position_first_ec=pcew_temp(1:3,kl_ew)
                    pcew_temp(1:3,kl_ew)=pcew_temp(1:3,i)
                    pcew_temp(1:3,i)=pc%position_first_ec

                    pc%z=pcew_temp(4,kl_ew)
                    pcew_temp(4,kl_ew)=pcew_temp(4,i)
                    pcew_temp(4,i)=pc%z
                    pc%c=0.0_rk
                    exit
                 endif ! dot_product
              enddo ! kl1,pcr_n
           enddo ! i=1,EWPC_N
           print*,'max_ewpc_pcr_dev for equivalent centers', max_ewpc_pcr_dev

           if(print_epe) &
           print*, 'found coinciding in charge centers in pcr and ewpc (kl_pcr,kl_ew)', kl_pcr,kl_ew

           write(output_unit,*) 'number of found coinciding in charge centers in pcr and ewpc',kl_pcr,kl_ew
           write(output_unit,*) 'indexes  for  last coinciding elements ' ,kl_ew,kl_pcr
!!** done equivalent atoms go first

           kl_ew_eq=kl_ew ! centers starting from this index will not be neglected
                          ! but treated with modified charge now this centers are sorted
                          ! and the charges are modified

           print*, 'locate centers which coinside but have different charges'
           do i=kl_ew+1,EWPC_N
              kl1=kl_pcr+1
              do k=kl1,pcr_n
                 if(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i)).lt.small_ew )  then
                    max_ewpc_pcr_dev=max(sqrt(dot_product(pcr_temp(1:3,k)-pcew_temp(1:3,i)  &
                      ,pcr_temp(1:3,k)-pcew_temp(1:3,i))),max_ewpc_pcr_dev)
                    kl_pcr=kl_pcr+1
                    kl_ew=kl_ew+1
                    pc%position_first_ec=pcr_temp(1:3,kl_pcr)
                    pcr_temp(1:3,kl_pcr)=pcr_temp(1:3,k)
                    pcr_temp(1:3,k)=pc%position_first_ec
                    pc%z=pcr_temp(4,kl_pcr)
                    pcr_temp(4,kl_pcr)=pcr_temp(4,k)
                    pcr_temp(4,k)=pc%z

                    pc%position_first_ec=pcew_temp(1:3,kl_ew)
                    pcew_temp(1:3,kl_ew)=pcew_temp(1:3,i)
                    pcew_temp(1:3,i)=pc%position_first_ec

                    pc%z=pcew_temp(4,kl_ew)
                    pcew_temp(4,kl_ew)=pcew_temp(4,i)
                    pcew_temp(4,i)=pc%z
                    pcew_temp(4,kl_ew)=pcew_temp(4,kl_ew)-pcr_temp(4,kl_pcr)
                    pcr_temp(4,kl_pcr)=0.0_rk
                    pc%c=0.0_rk
                    exit
                 endif! dot_product
              enddo! kl1,pcr_n
           enddo! i=1,EWPC_N
           if(print_epe) then
            print*,'max_ewpc_pcr_dev for other coinciding centers', max_ewpc_pcr_dev
            print*, 'coinciding centers in ewpc and pcr arrays  kl_ew kl_pcr',kl_ew,kl_pcr
            print*, 'charge of collected PCs, pcew_temp vs pcr_temp', &
                     sum(pcew_temp(4,1:kl_ew)),sum(pcr_temp(4,1+n_gxat_pcr:kl_pcr))
            DPRINT  'if not changed from previous value then no such centers'
           endif

           write(output_unit,*) &
                'number of coinciding centers in ewpc and pcr arrays & charges ',kl_ew &
                ,sum(pcew_temp(4,1:kl_ew)),sum(pcr_temp(4,1+n_gxat_pcr:kl_pcr))
!!** done second go centers with modified charges

        endif ex_pcr_n

        ! **   now first kl centers in pcew_temp and pcr_temp coincide
        !      calculate its interact with cluster to be compared with cluster_epe interact
        e_nuc_pcr=0.0_rk
        min_pcrdist=100.0
        do nb=1,kl_ew
         do na=1,N_unique_atoms
          do eq_a=1,unique_atoms(na)%N_equal_atoms
           dist=sqrt(sum((unique_atoms(na)%position(:,eq_a)  &
                                            -pcew_temp(1:3,nb))**2))
           min_pcrdist=min(min_pcrdist,dist)
           e_nuc_pcr=e_nuc_pcr+ &
                     pcew_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/dist
          enddo
         enddo
          DPRINT pcew_temp(1:3,nb),pcew_temp(4,nb),nb
        enddo
        print*,'min_pcrdist', min_pcrdist
        if(print_epe) &
        print*,'e_nuc_pcr to be compared with e_nuc_epe',e_nuc_pcr

        ! **  check potential of PC
        e_nuc_ewpc = 0.0_rk
        min_ewpcdist=100.0_rk
        do nb=1+kl_ew,EWPC_N
           j=0
           do na=1,N_unique_atoms
              do eq_a=1,unique_atoms(na)%N_equal_atoms
                 j=j+1
                 dist=sqrt(sum((unique_atoms(na)%position(:,eq_a) &
                      -pcew_temp(1:3,nb))**2))
                 if(dist < min_ewpcdist) then
                    min_ewpcdist=dist
                    ii=j
                 end if
                 e_nuc_ewpc = e_nuc_ewpc+pcew_temp(4,nb)*unique_atoms(na)%Z/dist
              enddo!  eq_a=1,unique_atoms(na)%N_equal_atoms
           enddo! na=1,N_unique_atoms
        enddo! 1,EWPC_N
        min_ewpcdist= min_ewpcdist*0.529177_rk
        if(min_ewpcdist < 1.0_rk) then
           write(trace_message,*) "WARNING: minimal distance between EWALD.PCR and QM cluster is ",min_ewpcdist, &
                "Angstrom for QM center ",ii
           call write_to_trace_unit(trim(trace_message))
        end if

#if 0
        print*,'min_ewpcdist',min_ewpcdist
        ASSERT(min_ewpcdist.gt.1.0_rk)
#endif
        DPRINT 'now first kl centers in pcew_temp and pcr_temp coincide'

        write(output_unit,*) 'e_nuc_ewpc with use of pcew_temp and Z', e_nuc_ewpc, &
                                                    sum(pcew_temp(4,1+kl_ew:EWPC_N))
        if(print_epe) print*,'e_nuc_ewpc with use of pcew_temp and Z', e_nuc_ewpc, &
                                                    sum(pcew_temp(4,1+kl_ew:EWPC_N))
        !** done

      if(pcr_n.ne.0) then
         if(print_epe) write(output_unit,*) '** take additional centers in pcew_temp',1+kl_ew
      else
         if(print_epe) write(output_unit,*) '** take PC from pcew_temp starting from ',1+kl_ew
      end if

!!** now reorder all atoms in ewpc_temp  starting from first non equivalent one
   i=1+kl_ew_eq
   it_ewpc=0   !types of EWPC
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

      call group_coset_decomp(n_equal,local_groups_pc,cosets_pc,point_trafos_pc%matrix)
      pc%N_equal_charges=n_equal

      ewa_allocstat(13)=1 ! point_trafos_pc%matrix
      ewa_allocstat(14)=1 ! cosets_pc%elements


      allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
      ASSERT(ewa_allocstat(15).eq.0)
      ewa_allocstat(15)=1 ! pc%position in temp

      ! determine positions of equal atoms
      !       write(output_unit,*) " Equal charges of Type ",i
      do j=1,n_equal
         position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets_pc%elements(1,j)),position)
         pc%position(1,j) = position2(1)
         pc%position(2,j) = position2(3)
         pc%position(3,j) = position2(2)
         do ii=i,EWPC_N
            if(dot_product(pcew_temp(1:3,ii)-pc%position(:,j) &
                 ,pcew_temp(1:3,ii)-pc%position(:,j)) < small_ew) then

               pcew_temp(1:3,ii)=pcew_temp(1:3,i)
               pcew_temp(1:3,i)=pc%position(1:3,j)

               pc%z=pcew_temp(4,ii)
               pc%c=0.0_rk
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
                 cosets_pc%elements, pc%position, stat=ewa_allocstat(12))
      ASSERT(ewa_allocstat(12).eq.0) ! local_groups_pc%elements
      ewa_allocstat(13)=0            ! point_trafos_pc%matrix
      ewa_allocstat(14)=0            ! coset%elements
      ewa_allocstat(15)=0            ! pc%position in temp

      ASSERT(n_equal.eq.pc%N_equal_charges)
    enddo ewpc_wh ! while
    !** done

    !**  check potential of PC
    e_nuc_ewpc = 0.0_rk
    do nb=1+kl_ew,EWPC_N
       do na=1,N_unique_atoms
          do eq_a=1,unique_atoms(na)%N_equal_atoms
           e_nuc_ewpc = e_nuc_ewpc+pcew_temp(4,nb)*(unique_atoms(na)%Z-unique_atoms(na)%ZC)/ &
                  sqrt(sum((unique_atoms(na)%position(:,eq_a)-pcew_temp(1:3,nb))**2))
          enddo !  eq_a=1,unique_atoms(na)%N_equal_atoms
       enddo ! na=1,N_unique_atoms
    enddo ! 1,EWPC_N
    write(output_unit,*) 'e_nuc_ewpc with use of pcew_temp & Z', e_nuc_ewpc, &
                                              sum(pcew_temp(4,1+kl_ew:EWPC_N))
    if(print_epe) then
    print*,'e_nuc_ewpc and z', e_nuc_ewpc,sum(pcew_temp(4,1+kl_ew:EWPC_N))
    print*,'epe coinciding centers are not taken into account'
    endif
    !** done

!!** now reorder still left  atoms in  pcr_temp
    it_pcr=0   !types of pcr
    epe_ex:if(pcr_n.gt.kl_pcr) then
       if(print_epe) print*, pcr_n-kl_pcr, &
                             ' not coinciding with ewpc environment epe centers  are found of total  charge', &
                             sum(pcr_temp(4,1+kl_pcr:pcr_n))
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
          call group_coset_decomp(n_equal,local_groups_pc,cosets_pc,point_trafos_pc%matrix)
          ewa_allocstat(13)=1 ! point_trafos_pc%matrix
          ewa_allocstat(14)=1 ! cosets_pc

          pc%N_equal_charges=n_equal

          allocate(pc%position(3,pc%N_equal_charges),stat=ewa_allocstat(15))
          ASSERT(ewa_allocstat(15).eq.0)
          ewa_allocstat(15)=1

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
                   pc%c=0.0_rk
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

    enddo pcr_wh! while
    !** done
    else epe_ex
       DPRINT 'no not coinciding with ewpc environment epe centers  exist'
    endif epe_ex

        deallocate(pc,stat=ewa_allocstat(16))
        ASSERT(ewa_allocstat(16).eq.0)
        nullify(pc)

#if 0
print*,'read_pc_array'
call read_pc_array(pcc_array,'pcc_array',pcc_n)
call read_pc_array(pcs_array,'pcs_array',pcs_n)
print*,'read_pc_array done'
#endif

    n_unique_pcr=pcs_n+pcc_n
    if(pcr_n.gt.kl_pcr) n_unique_pcr=n_unique_pcr+it_pcr
    allocate(ewpc_array(it_ewpc+n_unique_pcr),  stat=ewa_allocstat(3))
    ASSERT(ewa_allocstat(3).eq.0)
           ewa_allocstat(3)=1




!!$        write(output_unit,  *)  '** store not coinciding pcew centers'
       n_ewpc_arr=0
        i=1+kl_ew_eq
       if(print_epe) print*,' ewpc_array (ewald.pcr part) :'

!      first go ewald_pcr unique centers
!      second pcs and pcc centers


       ewpcarr_wh: do while(i.le.EWPC_N)
           ASSERT(pcew_no(i).le.size(ewpc_array))
           pc=>ewpc_array(pcew_no(i))
           n_equal=0
           pc%z=0.0_rk
           pc%c=0.0_rk
           do ii=i,ewpc_N
              n_equal=n_equal+1
              pc%z=pc%z+pcew_temp(4,ii)
              if(ii.eq.ewpc_N) exit
              if(pcew_no(ii+1).ne.pcew_no(ii)) exit
           enddo !ii=i,ewpc_N
           pc%z=pc%z/n_equal
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(4))
           ASSERT(ewa_allocstat(4).eq.0)
           ewa_allocstat(4)=1
           pc%name='ewpc    '
           pc%N_equal_charges= n_equal
           pc%position_first_ec=pcew_temp(1:3,i)
           do ii=1,n_equal
              pc%position(:,ii)=pcew_temp(1:3,i+ii-1)
              n_ewpc_arr=n_ewpc_arr+1
!!!     print all generated charges
#if 0 /* print all ewpc */
              if(print_epe) &
               print '(4f15.8,i4,4i4)',pc%position(:,ii),pc%z,n_equal,ii,pcew_no(i),n_ewpc_arr
               write(output_unit, '(4f15.8,i4,4i4)') pc%position(:,ii), &
                  pc%z,n_equal,0,0,1,0
#endif
           enddo ! ii=1,equal
           !       print only first atoms in the group
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8)')  &
!       print '(4f15.8,i3,4i2,f15.8)' &
!                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!                                  sqrt(sum((unique_atoms(1)%position(:,1) &
!                                                    -pc%position(:,1))**2))
!!$           write(output_unit, '(4f15.8,i3,4i2,f15.8,i4)') &
!       if(print_epe) print '(4f15.8,i4,4i4,f15.8,i4)', &
!                pc%position(:,1),pc%z,n_equal,0,0,1,0, &
!                sqrt(sum((unique_atoms(1)%position(:,1) &
!                     -pc%position(:,1))**2)),(pcew_no(i))
        i=i+n_equal
        enddo ewpcarr_wh

        write(output_unit,*) &
             'number of groups of symmetry equavalent ewpc PC',it_ewpc
       ewpc_n=it_ewpc
       DPRINT 'it_ewpc',it_ewpc

        !** done

     pcr_n_ex:   if(pcr_n.ne.0) then
        write(output_unit,*) &
             '** store not coinciding pcr centers'
        if(print_epe) print*,'** store not coinciding pcr centers_____________', &
                               size(pcc_array),size(pcs_array)
        i=1+kl_pcr
       ewpcarr_pcr: do while(i.le.pcr_n) ! total number of atoms in epe.pcr
ASSERT(pcr_no(i)+ewpc_n.le.size(ewpc_array))
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
           pc%c=0.0_rk
           allocate(pc%position(3,n_equal),stat=ewa_allocstat(4))
           ASSERT(ewa_allocstat(4).eq.0)
           ewa_allocstat(4)=1
           pc%name='pcr    '
           pc%N_equal_charges= n_equal
           pc%position_first_ec=pcr_temp(1:3,i)
           do ii=1,n_equal
              pc%position(:,ii)=pcr_temp(1:3,i+ii-1)
              n_ewpc_arr=n_ewpc_arr+1
!!!     print all generated charges
#if 0 /* print all pcr */
        if(print_epe) print '(4f15.8,i4,4i4)',pc%position(:,ii),pc%z,n_equal,ii,pcr_no(i)+ewpc_n,n_ewpc_arr
#endif
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
         write(output_unit,*) 'number of groups of symmetry equavalent pcr PC',it_pcr
         ewpc_n=ewpc_n+it_pcr
         DPRINT 'it_ewpc + it_pcr', ewpc_n
        !** done
     endif pcr_n_ex



    start_regI_epevar=ewpc_n+1
    cluster_nuc_epe_en=0.0_rk

    store_pcs: if(pcs_n.ne.0) then !** store pcs centers and calculate
                                   !   cluster_nuc_epe_en for epe shells
       if(print_epe) print*,' ewpc_array (epe.pcs part) :'
       min_pcsdist=100.0
       do na=1,size(pcs_array)
          pc=>ewpc_array(na+ewpc_n)
          allocate(pc%position(3,pcs_array(na)%n_equal_charges))
!          pc=pcs_array(na)
          pc%z=pcs_array(na)%z
          pc%c=0.0_rk
          pc%name=pcs_array(na)%name
          pc%n_equal_charges=pcs_array(na)%n_equal_charges
          pc%position=pcs_array(na)%position
          pc%position_first_ec=pcs_array(na)%position_first_ec

#if 0 /* delete deference with ewpc_temp */
         pc%position=pcr_array(na)%position
#endif

#if 1 /* control cluster_nuc_epe_en */
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

                else ! reg atom
                   do eq_a=1,unique_atoms(i)%N_equal_atoms
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                                           -pc%position(:,eq_pcs))**2))
                      min_pcsdist=min(min_pcsdist,dist)
                      cluster_nuc_epe_en=cluster_nuc_epe_en+ &
                       pc%z*(unique_atoms(i)%Z-unique_atoms(i)%zc)/dist
                   enddo
             end if

             enddo
            n_ewpc_arr=n_ewpc_arr+1
#if 0 /* print all pcc */
            if(print_epe) print '(4f15.8,i4,4i4,i5)', &
                    pc%position(:,eq_pcs),pc%z,pc%n_equal_charges,eq_pcs,na+ewpc_n,n_ewpc_arr
#endif
          enddo
!           print '(4f15.8,i3,4i2,i5)',pc%position(:,1),pc%z,pc%n_equal_charges,0,0,1,0,na+ewpc_n

        DPRINT na, cluster_nuc_epe_en,pc%N_equal_charges,pc%z,sum(pc%position(:,pc%N_equal_charges))
#endif
       enddo ! N_unique_pcr
       DPRINT 'min_pcsdist', min_pcsdist
       DPRINT 'cluster_nuc_epe_en pcs', cluster_nuc_epe_en
       ewpc_n=ewpc_n+size(pcs_array)
    endif store_pcs
#if 0
call read_pc_array(ewpc_array,'ewpc_arra',ewpc_n)
         if(associated(pcc_array)) then
          do na=1,size(pcc_array)
           deallocate(pcc_array(na)%position,stat=ewa_allocstat(26))
           ASSERT(ewa_allocstat(26).eq.0)
          enddo
         deallocate(pcc_array,stat=ewa_allocstat(19))
         ASSERT(ewa_allocstat(19).eq.0)
        endif

!        if(pcs_n.ne.0) then
        if(associated(pcs_array)) then
          do na=1,size(pcs_array)
           deallocate(pcs_array(na)%position,stat=ewa_allocstat(25))
           ASSERT(ewa_allocstat(25).eq.0)
          enddo
           deallocate(pcs_array,stat=ewa_allocstat(18))
          ASSERT(ewa_allocstat(18).eq.0)
        endif
call print_e_nuc_ewpc()
return
#endif

    store_pcc: if(pcc_n.ne.0) then !** store pcc centers and calculate
                                   !   cluster_nuc_epe_en for epe cores
       if(print_epe) print*,' ewpc_array (epe.pcc part) :'
       nullify(pc)
       min_pccdist=100.0

       do na=1,size(pcc_array)
ASSERT(na+ewpc_n.le.size(ewpc_array))
          pc=>ewpc_array(na+ewpc_n)
!          pc=pcc_array(na)
          allocate(pc%position(3,pcc_array(na)%n_equal_charges))
!          pc=pcs_array(na)
          pc%z=pcc_array(na)%z
          pc%c=0.0_rk
          pc%name=pcc_array(na)%name
          pc%n_equal_charges=pcc_array(na)%n_equal_charges
          pc%position=pcc_array(na)%position
          pc%position_first_ec=pcc_array(na)%position_first_ec
#if 0 /* delete deference with ewpc_temp */
         pc%position=pcr_array(na)%position
#endif

#if 1 /* control cluster_nuc_epe_en */
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
                      min_pccdist=min(dist,min_pccdist)
                      cluster_nuc_epe_en=cluster_nuc_epe_en &
                      +pc%z*(unique_atoms(i)%Z-unique_atoms(i)%ZC)/dist
                   enddo
                end if
             enddo
           n_ewpc_arr=n_ewpc_arr+1
#if 0 /* print all pcc */
           if(print_epe) print '(4f15.8,i4,4i4,i5)', &
                         pc%position(:,eq_pcc),pc%z,pc%n_equal_charges,eq_pcc,na+ewpc_n,n_ewpc_arr
#endif
          enddo
!       print '(4f15.8,i3,4i2,i5)',pc%position(:,1),pc%z,pc%n_equal_charges,0,0,1,0,na+ewpc_n
#endif
       enddo ! N_unique_pcr
       DPRINT 'min_pccdist',min_pccdist
       ewpc_n=ewpc_n+size(pcc_array)
    endif store_pcc

       deallocate(pcew_temp,stat=ewa_allocstat(21))
        ASSERT(ewa_allocstat(21).eq.0)

    ewpc_array_charge=0.0_rk
    do  na=1,ewpc_n
    ewpc_array_charge = ewpc_array_charge+ewpc_array(na)%z*ewpc_array(na)%n_equal_charges
    enddo

   if(print_epe) print*, 'ewpc_array_charge ', ewpc_array_charge

        write(output_unit,*) 'final  number of ewald and epe centers '
        write(output_unit,*) 'to model external field actin on claster '
        write(output_unit,*) ewpc_n
        write(output_unit,*)
        write(output_unit,*) 'energy of interaction of the claster nuclei'
        write(output_unit,*) 'with epe only centers, cluster_nuc_epe_en'
        write(output_unit,*)  cluster_nuc_epe_en

        if(print_epe) &
        print*, 'cluster_nuc_epe_en calculated with ewpc_array coinciding and the same charges ', cluster_nuc_epe_en


       if(print_epe) print*, 'PCR_N----', pcr_n
       if(pcr_n.ne.0) then
!           deallocate(pcr_temp,pcr_no,psb_ind, stat=ewa_allocstat(9))
           deallocate(pcr_temp,pcr_no, stat=ewa_allocstat(9))
           ASSERT(ewa_allocstat(9).eq.0)
        end if


!       if(pcr_n.ne.0) then
       if(associated(pcr_array)) then
          do na=1,size(pcr_array)
           deallocate(pcr_array(na)%position,stat=ewa_allocstat(24))
           ASSERT(ewa_allocstat(24).eq.0)
          enddo
          deallocate(pcr_array,stat=ewa_allocstat(17))
          ASSERT(ewa_allocstat(17).eq.0)
          print*,'pcr_array deallocated'
        endif

!        if(pcc_n.ne.0)  then
         if(associated(pcc_array)) then
          do na=1,size(pcc_array)
           deallocate(pcc_array(na)%position,stat=ewa_allocstat(26))
           ASSERT(ewa_allocstat(26).eq.0)
          enddo
         deallocate(pcc_array,stat=ewa_allocstat(19))
         ASSERT(ewa_allocstat(19).eq.0)
        endif

!        if(pcs_n.ne.0) then
        if(associated(pcs_array)) then
          do na=1,size(pcs_array)
           deallocate(pcs_array(na)%position,stat=ewa_allocstat(25))
           ASSERT(ewa_allocstat(25).eq.0)
          enddo
           deallocate(pcs_array,stat=ewa_allocstat(18))
          ASSERT(ewa_allocstat(18).eq.0)
        endif



        if(allocated(pcew_no)) then
         deallocate(pcew_no,stat=ewa_allocstat(27))
         ASSERT(ewa_allocstat(27).eq.0)
        endif

       nullify(pc)
       e_nuc_ewpc=0.0_rk
!       print*, n_ewpc_arr, 'n_ewpc_arr xx'
       print*

       min_ewpcdist=100.0_rk
       do na=1,ewpc_n
          pc=>ewpc_array(na)


          do eq_pcc=1,pc%n_equal_charges
!         print*,' F',pc%position(:,eq_pcc)*0.529,pc%z,na
             j=0
             do i=1,N_unique_atoms+n_timps
                if(i.gt.N_unique_atoms) then
                   ua=>unique_timps(i-N_unique_atoms)
                   do eq_a=1,ua%N_equal_atoms
                      dist=sqrt(sum((ua%position(:,eq_a) &
                           -pc%position(:,eq_pcc))**2))
                      e_nuc_ewpc=e_nuc_ewpc+pc%z*(ua%Z-ua%ZC)/dist
                   enddo
                else
                   do eq_a=1,unique_atoms(i)%N_equal_atoms
                      j=j+1
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                                           -pc%position(:,eq_pcc))**2))
                      e_nuc_ewpc=e_nuc_ewpc &
                      +pc%z*(unique_atoms(i)%Z-unique_atoms(i)%ZC)/dist

                      if(dist < min_ewpcdist) then
                         min_ewpcdist=dist
                         ii=j
                      end if
                   enddo
                end if
             enddo
          enddo
       enddo ! N_unique_pcr
       min_ewpcdist= min_ewpcdist*0.529177_rk
       if(min_ewpcdist < 1.0_rk) then
          write(trace_message,*) "WARNING: minimal distance between EPE point charges and QM cluster is ",min_ewpcdist, &
               "Angstrom for QM center ",ii
          call write_to_trace_unit(trim(trace_message))
       end if
   print*,'e_nuc_ewpc with final ewpc_array',e_nuc_ewpc,size(ewpc_array),ewpc_n

#if 1 /* check complite external field */
   print*, 'complete potentials on atoms'
       ii=0
       do i=1,N_unique_atoms+n_timps
       do eq_a=1,unique_atoms(i)%N_equal_atoms
       ii=ii+1
       e_nuc_ewpc=0.0_rk
       gv_nuc_ewpc=0.0_rk
       do na=1,ewpc_n
          pc=>ewpc_array(na)
          do eq_pcc=1,pc%n_equal_charges
                      dist=sqrt(sum((unique_atoms(i)%position(:,eq_a) &
                                           -pc%position(:,eq_pcc))**2))
                      e_nuc_ewpc=e_nuc_ewpc+pc%z/dist
                      gv_nuc_ewpc=gv_nuc_ewpc+ pc%z* &
                       (unique_atoms(i)%position(:,eq_a)-pc%position(:,eq_pcc))/dist**3
                   enddo
             enddo
        print *, ii,e_nuc_ewpc,gv_nuc_ewpc,sum(gv_nuc_ewpc(:)**2)
          enddo
       enddo ! N_unique_pcr
#endif

!call write_pc_array('ewpc_arra')

  endif ! ewaldpc
#if 0
contains

subroutine print_e_nuc_ewpc()
use  energy_calc_module, only: e_nuc_ewpc
print*, e_nuc_ewpc()
end subroutine print_e_nuc_ewpc

subroutine write_pc_array(name)
!type(pointcharge_type),intent(in):: pc_array(:)
character(len=9), intent(in):: name
integer(i4_kind):: pc_unit,i,eq
  pc_unit=openget_iounit(file=trim(inpfile(name)), form='FORMATTED',status='unknown')
write(pc_unit,*) size(ewpc_array)
do i=1,size(ewpc_array)
ewpc_array(i)%position_first_ec=ewpc_array(i)%position(:,1)
write(pc_unit,*) ewpc_array(i)%z,ewpc_array(i)%n_equal_charges,ewpc_array(i)%position_first_ec,ewpc_array(i)%name,i
do eq=1,ewpc_array(i)%n_equal_charges
write(pc_unit,*) ewpc_array(i)%position(:,eq)
enddo
enddo
call returnclose_iounit(pc_unit,'keep')
print *,name//' written '//trim(inpfile(name))

end subroutine write_pc_array


subroutine read_pcr_temp(pcr_n)
integer(i4_kind):: pc_unit,i,pcr_n
pc_unit=openget_iounit(file=trim(inpfile('pcr_temp')), form='FORMATTED',status='unknown')
read(pc_unit,*) pcr_n
if(.not.allocated(pcr_temp)) allocate(pcr_temp(4,pcr_n),pcr_no(pcr_n))
do i=1,size(pcr_temp,2)
read(pc_unit,*) pcr_temp(:,i),pcr_no(i)
enddo
call returnclose_iounit(pc_unit,'keep')
end subroutine read_pcr_temp

subroutine read_pc_array(pc_array,name,pc_n)
integer(i4_kind),intent(out):: pc_n
type(pointcharge_type),pointer:: pc_array(:)
character(len=9), intent(in):: name
integer(i4_kind):: pc_unit,i,eq
  pc_unit=openget_iounit(file=trim(inpfile(name)), form='FORMATTED', status='unknown')
read(pc_unit,*) pc_n
print*,'pc_n',pc_n
if(.not.associated(pc_array)) allocate(pc_array(pc_n))
do i=1,pc_n
read(pc_unit,*) pc_array(i)%z,pc_array(i)%n_equal_charges,pc_array(i)%position_first_ec,pc_array(i)%name
if(.not.associated(pc_array(i)%position))&
allocate(pc_array(i)%position(3,pc_array(i)%n_equal_charges))
do eq=1,pc_array(i)%n_equal_charges
read(pc_unit,*) pc_array(i)%position(:,eq)
enddo
enddo
call returnclose_iounit(pc_unit,'keep')
print *,name//' read '//trim(inpfile(name))

end subroutine read_pc_array
#endif

 end subroutine symm_ewpc_gen

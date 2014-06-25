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
module epe_pg_module

# include <def.h>
use type_module
use iounitadmin_module
  use   epecom_module
  use   culon_module
  use   mol_module, only : epe_generater
  use   epepar_module, only: ec_name_of_type=>name_of_type, &
        ec_max_type_ions=>max_type_ions

  implicit none
   save
  private

public :: pggener, get_pgcent, epe2pg, an_epe_type, nhds, lato,epe2pg_epe

  integer(kind=i4_kind) :: ncb,ncr,i,j
  real(kind=r8_kind) :: ckk,skl
  real(kind=r8_kind) :: gma(3,3)
  real (r8_kind), dimension(ndt)::an_epe_type =0.0_r8_kind
  integer (kind=i4_kind):: gx_unit,r_unit,opt_unit &
                          ,pcs_unit,pcr_unit,pcc_unit
  integer (kind=i4_kind),allocatable,dimension(:)::nhds
  logical, allocatable, dimension(:)::lato
  logical, public :: operations_make_reg_reference=.false.
contains

subroutine pggener

!-------------------------------------------------------------------------
!       *** Parameters  and discription ***
!       cluster atoms are specified by use of parameter in a vector
!       impu. If impu(i).lt.0 no new hades model parameters are 
!       inputed
!       at imput atoms of cluster are rotated and shifted to the hades
!       cristall system
!       at output all atoms of the set n1 and N2A are shifted and rotated
!       back
!
!       interface block connecting hades with quantum chemical prog.
!
!       interfaced_mode=1
!       interfaced_mode=2
!       interfaced_mode.lt.0
!
!       gx file contain atoms to be considered as impurities !
!------------------------------------------------------------------------


  integer(kind=i4_kind) :: L(3),M(3),nct
  integer(kind=i4_kind) :: istat
  logical, allocatable, dimension(:)::limpu
  logical :: f_exist
  character(len=4) bl(10) !/10*'    '/
  bl='    '

  DPRINT 'pggener: '//trim(epe_input_dir)//'/gxfile'
  inquire(file=trim(epe_input_dir)//'/gxfile',exist=f_exist)
  if(.not.f_exist) then
     call error_handler("pggener: gxfile not existed")
  end if
  gx_unit=openget_iounit(trim(epe_input_dir)//'/gxfile', &
                           status='old',form='formatted')

  natoms=0      ! number of regular atoms
  i=0
  do
    i=i+1
    read(gx_unit,*,err=100) an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
    write(output_epe,203) an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
    DPRINT '(f7.2,3(2X,f13.7),2i4,2x,3I3,2X,3I3,2x,i3)', an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
203 FORMAT((f7.2,3(2X,f13.7),2i4,2x,3I4,2X,3I4,2x,i3))
    if(an(i).le.0.5) exit
  enddo 

        natoms=i-1
        allocate(nhds(natoms),limpu(natoms),stat=epealloc_stat(12))
        ASSERT(epealloc_stat(12).eq.0)
               epealloc_stat(12)=1

        limpu=impu(1:natoms).lt.0
        call returnclose_iounit(gx_unit)

  inquire(file=trim(epe_input_dir)//'/epe.r',exist=f_exist)
  if(.not.f_exist) then
#ifdef NEW_EPE
     if(operations_make_reg_reference) then
        r_unit=openget_iounit(trim(epe_input_dir)//'/epe.r')

        do i=1,natoms
           if(impu(i) == -1) write(r_unit,'(3f25.12)') gx_centers(:,i)
        end do

        call returnclose_iounit(r_unit)
        print*, 'epe.r file has been created'
     else
#endif
        call error_handler("pggener: epe.r not existed")
#ifdef NEW_EPE
     end if
#endif
  end if

  r_unit=openget_iounit(trim(epe_input_dir)//'/epe.r', &
                          status='old', form='formatted')
  write(output_epe,*) ' atoms of epe.r file'

       allocate(gx_imp_r(3,natoms),stat=epealloc_stat(5))
        if(epealloc_stat(5).ne.0) call error_handler("allocate gx_imp_r failed")
        epealloc_stat(5)=1

  gx_imp_r=0.0_r8_kind
  do i=1,natoms
       if(limpu(i)) then
            read(r_unit,*,err=101) gx_imp_r(:,i)
            write(output_epe,*) gx_imp_r(:,i)
        endif 
  enddo ! i=1,natoms
  call returnclose_iounit(r_unit)

! **searching of atomic centers relations
! **shift and rotate cluster  to the epe lattice system
  call search_centre_and_rotation &
       (natoms,gx_imp_r,rot_gto_to_epe,shft_gto_to_epe,1)       
  GMA=rot_gto_to_epe

  CALL MINV(GMA,3,VC,L,M)
  ncb=n_centres_of_generation   
                ! trough the regular EPE input (0 is regular value)
  ncr=n_centres_of_generation

  ncb=ncb+count(limpu)

 allocate(lato(natoms),stat=epealloc_stat(11))
 ASSERT(epealloc_stat(11).eq.0)
        epealloc_stat(11)=1
 lato=indgx(1,1:natoms).ne.0

  do i=1,natoms
    if(limpu(i)) then    
      if(impu(i).gt.0) then
        ncb=ncb+1
        nct=ncb
      else      ! impu(i).lt.0
        n_centres_of_generation=n_centres_of_generation+1
        nct=n_centres_of_generation
      endif     !impu(i).gt.0/else

      R_CENT_GENER(nct,:)= gx_imp_r(:,i)*auangs
    endif       ! limpu(i)
  enddo !  i=1,natoms
 deallocate(limpu,stat=istat)
 if(istat.ne.0) call error_handler("deallocate limpu failed")

  write(output_epe,*)'coordinates and numbers of the generation centers '
        do nct=1,ncb
         write(output_epe, '((3F12.6,I5))')R_CENT_GENER(nct,:),nct
        enddo

! **initial generation of lattice to define positions   !
! **of generation centers in the regular structure      !
! **no generation center for impu.ge.0                  !

  write(output_epe,*) 'gener klgen', klgen
 nhds=0

 ! print*,'call epe_generater(klgen,.false.)'

  call epe_generater(klgen,.false.)
 ! print*,'done epe_generater'
  n_centres_of_generation=ncr
  do j=1,natoms
    if(impu(j).ge.0) cycle
    do i=1,klgen
      if(dot_product(gx_imp_r(:,j)*auangs-epe(i)%r, &
                        gx_imp_r(:,j)*auangs-epe(i)%r).gt.0.1_r8_kind) cycle
      n_centres_of_generation=n_centres_of_generation+1
      R_CENT_GENER(n_centres_of_generation,:)=epe(i)%r
    enddo !i=1,klgen
  enddo !j=1,natoms
  n_centres_of_generation=ncb
  write(output_epe,*) 'number of PG related centers of generation', n_centres_of_generation 
  do i=1,n_centres_of_generation
    write(output_epe,*) R_CENT_GENER(i,:)
  enddo ! i=,n_centres_of_generation

  return

100 call error_handler("epe_pg_module: pggener: error of reading in GXFILE ")
101 call error_handler("epe_pg_module: pggener: error of reading in EPE.R ")

end subroutine pggener
!--------------------------------------------------------------

subroutine   get_pgcent
! ** treat specification of impurities  **
  use   epecom_module
  use   culon_module
  use comm_module

  real(kind=r8_kind) :: q_imp,q_sh_imp,pk_imp
  real(kind=r8_kind) :: bim,roim,cim,dim
  integer(kind=i4_kind),dimension(:),allocatable :: n_im2epe
  integer(kind=i4_kind) :: istat
  integer(kind=i4_kind) :: j,n_pgimp,NVC,i,k,type_imp,mti, &
                           reg_ref_unit,epe_ref_unit

  namelist /dvmhds_impurities/ type_imp,q_imp,q_sh_imp,pk_imp
  namelist /dvmhds_imp_param_of_potential/ bim,roim,cim,dim

  mti=max_type_ions

  write(output_epe,*) ' -----------------------------------------------'
  write(output_epe,*) ' Specification of impurities '
  write(output_epe,*) (impu(j),j=1,natoms)

  n_pgimp=0
  nvc=0

  allocate(n_im2epe(natoms),stat=epealloc_stat(20))
  ASSERT(epealloc_stat(20).eq.0)
  epealloc_stat(20)=1

  if(make_epe_reference) then
     epe_ref_unit=openget_iounit(trim(epe_input_dir)//'/epe_reference',&
          status='unknown', form='formatted')
  endif

!!$  print*,'define impurities types'
   if(explicit_coupling) then
     do i=1,reg_I_n_ions
                 do k=1,ec_max_type_ions
                    if(abs(an_epe_type(epe(i)%k)-ec_name_of_type(k)).lt.0.0001) then
                     epe(i)%ec=k
                       exit
                    endif
                    if(k.eq.ec_max_type_ions) then
                         print*,'an in gxfile ',an_epe_type(epe(i)%k)
                         call error_handler('epe_pg_module: type not found')
                    endif
                 enddo

     enddo
    endif

  do j=1,natoms
     do i=1,reg_I_n_ions
        if(impu(j).le.0.and.lato(j)) then 
           if(dot_product(gx_imp_r(:,j)*auangs-epe(i)%r, &
                gx_imp_r(:,j)*auangs-epe(i)%r).gt.0.1) cycle
           nhds(j)=i
           nvc=nvc+1
           if(impu(j).eq.0) cycle
           n_pgimp=n_pgimp+1      ! number of PG related vacancies
           n_im2epe(n_pgimp)=i  
           n_impurities=n_impurities+1
           ASSERT(n_impurities.le.ndpt)

           nhds(j)=n_impurities
           !      **set positions of impurites to those of gx file **

            r_imp(n_impurities,:)=epe(i)%r

           if(make_epe_reference)  then
            print*,'write epe_reference',epe_ref_unit
                write(epe_ref_unit,'(3f25.14)') gx_centers(:,j)*auangs ! r_imp(n_impurities,:)
            print*,'write epe_reference done'
           endif

           R_NUC_IMP (n_impurities,:)=gx_centers(:,j)*auangs ! R_IMP(n_impurities,:)
           R_NUC_IMPO(n_impurities,:)=gx_centers(:,j)*auangs ! R_IMP(n_impurities,:)
           R_SH_IMP(n_impurities,:)=gx_centers(:,j)*auangs ! R_IMP(n_impurities,:)
           R_SH_IMPO(n_impurities,:)=gx_centers(:,j)*auangs ! R_IMP(n_impurities,:)

           if(qm_interfaced_mode.and.impu(j).lt.0) then 
              !take parameters of corresponding  parent epe
              q_impurity(n_impurities)=q_ion(epe(i)%k)  !ions
              q_sh_impurity(n_impurities)=q_shell(epe(i)%k)
              q_nuc_impurity(n_impurities)=q_nuclear(epe(i)%k)
              pk_impurity(n_impurities)=pk(epe(i)%k)
              TYPE_IMPURITY(n_impurities)=epe(i)%k
              imp2center(n_impurities)=i ! temporarly commented
              if(explicit_coupling) then
                 do k=1,ec_max_type_ions
                    if(an(j)-aint(an(j)).lt.0.0001) cycle
                    if(abs(an(j)-ec_name_of_type(k)).lt.0.0001) then
                       explicit_coupling_type(n_impurities)=k
                       epe(i)%ec=k
                       exit
                    endif
                    if(k.eq.ec_max_type_ions) then
                         print*,'an in gxfile ',an(j)
                         call error_handler('epe_pg_module: type not found')
                    endif
                 enddo

              endif
           endif
        endif
     enddo !i
  enddo !j=1,natoms

!  deallocate(lato,stat=epealloc_stat(11))
!  ASSERT(epealloc_stat(11).eq.0)
!  print*,'lato deallocated'

        DPRINT ' types of impurities'
        if(explicit_coupling) then
        do i=1,n_impurities
           DPRINT  i, explicit_coupling_type(i)
        end do

        else
        do i=1,n_impurities
           DPRINT  i, type_impurity(i)
        end do
     endif

     if(make_epe_reference) then
        call returnclose_iounit(epe_ref_unit)
        print*,'______________________________________epe_reference created'
        call write_to_trace_unit('epe_reference created')
     end if


        if(operations_make_reg_reference) then
           reg_ref_unit=openget_iounit(trim(epe_input_dir)//'/reg_reference', &
                status='unknown', form='formatted')
           do j=1,n_pgimp
              write(reg_ref_unit,'(3f25.14)') r_sh_ion(n_im2epe(j),:)
              write(reg_ref_unit,'(3f25.14)') r_nuc_ion(n_im2epe(j),:)
              R_SH_IMP(j,:)=r_sh_ion(n_im2epe(j),:)
              R_NUC_IMP(j,:)=r_nuc_ion(n_im2epe(j),:)
           enddo
           call returnclose_iounit(reg_ref_unit)
           print*,' reg_reference file with equilibrium positions of'
           print*,' EPE centers of the impurity cluster is created'
        endif

        if(use_epe_reference) then
           allocate(reg_reference(n_pgimp),epe_reference(n_pgimp),stat=istat)
           if(istat.ne.0) call error_handler(" allocate reg_reference failed")
!!$     print*,'get_pgcent: reg_reference: ', &
!!$     trim(epe_input_dir)//'/reg_reference'
           reg_ref_unit=openget_iounit(trim(epe_input_dir)//'/reg_reference', &
                status='old', form='formatted')
           do j=1,n_pgimp
              read(reg_ref_unit,*,err=110) reg_reference(j)%rs
              read(reg_ref_unit,*,err=110) reg_reference(j)%rc
           enddo
           call returnclose_iounit(reg_ref_unit)
!!$           print *,'reg_reference is read'

           epe_ref_unit=openget_iounit(trim(epe_input_dir)//'/epe_reference', &
                status='unknown', form='formatted')
           do j=1,n_pgimp
              read(epe_ref_unit,*,err=111) epe_reference(j)%rs
              epe_reference(j)%rc=epe_reference(j)%rs
           enddo
           ! For epe reference coordinatas of cores and shells coincide and
           ! they are actually read from gx file.

           call returnclose_iounit(epe_ref_unit)
!!$           print *,'epe_reference is read'
        endif ! use_epe_reference
        ! regular reference is defined with cores and shells

        deallocate(n_im2epe,stat=istat)
        if(istat.ne.0) call error_handler("deallocate n_im2epe failed")

        write(output_epe,*)'coordinates of EPE vacancies'
        do i=1,n_pgimp
        write(output_epe,*) epe(i)%r
        enddo 
  n_vacancies=n_vacancies+n_pgimp
  write (output_epe,102)n_pgimp
102 format(' number of escapac related substituens ',i3)
  write (output_epe,101) nvc,nhds(1:natoms)
101 format(' number of escapac related centers ', i3  &
         /' their serial numbers'/(1x,35i3))
  write (output_epe,*) '-------------------------------------------------'

  return

110 call error_handler("epe_pg_module: get_pgcent: error of reading in REG_REFERENCE")
111 call error_handler("epe_pg_module: get_pgcent: error of reading in EPE_REFERENCE")
end subroutine get_pgcent
!----------------------------------------------------

subroutine epe2pg_epe()
! **define corrections to the energy and the potentials for all
! **atoms of cluster
! **write coordinates of the Hades centers to file for use in
! **PG code
  use   epecom_module
  use   culon_module

  real(kind=r8_kind),allocatable,dimension(:,:):: r_in_pg_system
  integer(kind=i4_kind) :: i,up_limit
  integer(kind=i4_kind) :: istat,nn_tmp
  real(kind=r8_kind),allocatable :: q_tmp(:)
  logical :: stat
  logical exist_pcr

   type psb_centers
      real(kind=r8_kind)::rb(3),qb,x
      integer(kind=i4_kind):: linkto,linkfrom(3)
   end type psb_centers
   type(psb_centers),allocatable ::psb_pc(:)
   integer(kind=i4_kind)::link_count
   real(kind=r8_kind) ::qs_psbond,rc_psbond(3),qc_psbond

  !print*,' **make backward transformation for atoms of cluster ' 
  call search_centre_and_rotation &
        (natoms,gx_centers,gma,shft_gto_to_epe,-1)
  gx_unit=openget_iounit(trim(epe_input_dir)//'/gxfile', &
                          status='old',form='formatted')
  DPRINT 'epe2pg: '//trim(epe_input_dir)//'/gxfile unit=',gx_unit
  do i=1,natoms+1
203 FORMAT((f7.2,3(2X,f13.7),2i4,2x,3I4,2X,3I4,2x,i3))
        if(gx_highprec) then
        write(gx_unit,"(f7.2,3(2X,f21.12),2i4,2x,3I4,2X,3I4,2x,i3)") &
             an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
        else
        write(gx_unit,203) an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
        endif
  enddo !i=1,natoms+1 

  if(qm_interfaced_mode) then
    write(gx_unit,'(2f20.7,3i5)') etot_epe, etot_epe,81,1,0      
  endif !abs(interfaced_mode).eq.1) 
  call returnclose_iounit(gx_unit)

  !create/treat files to impose PC field on the 
  !EPE embedded cluster

  if(extended_epe) then
!!$  write(output_epe,*) 'extended_epe---------------',extended_epe
     up_limit=reg_2a_n_ions
  else
     up_limit=reg_I_n_ions   
  end if

  if(reg_2a_treated) then
!!$     print*,'reg_2a_treated------------------',reg_2a_treated
     end_treated_region= reg_2a_n_ions
  else
     end_treated_region=reg_I_n_ions
  end if

  inquire(file=trim(epe_input_dir)//'/epe.pcr', exist=exist_pcr)

  if(exist_pcr.and.qm_interfaced_mode.and.operations_make_reg_reference) &
       then
     call system('rm '//trim(epe_input_dir)//'/epe.pcr')
     exist_pcr=.false.
   endif

  pcr_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcr', &
                                status='unknown',form='formatted')
  !print*,'pseudobond_approach',  pseudobond_approach
  if(pseudobond_approach) then
     allocate(psb_pc(reg_I_n_ions),stat=istat)
     if(istat.ne.0) then
        call error_handler("allocate psb_pc failed")
     end if
     psb_pc(:)%x=9.0_r8_kind
     psb_pc(:)%linkto=0_i4_kind
     psb_pc(:)%linkfrom(1)=0_i4_kind
     psb_pc(:)%linkfrom(2)=0_i4_kind
     psb_pc(:)%linkfrom(3)=0_i4_kind

     if(exist_pcr) then
        read(pcr_unit,*) i
        if(i.ne.up_limit-n_vacancies) then
           print*,' inconsistent number in epe.pcr file'
        end if
        do i=n_vacancies+1,reg_I_n_ions
           read(pcr_unit,*)  &
                psb_pc(i)%rb ,psb_pc(i)%qb ,psb_pc(i)%linkto,psb_pc(i)%linkfrom
        end do
     end if

  end if
  
        epe(:reg_2a_n_ions)%s(1)=r_sh_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%s(2)=r_sh_ion(1:reg_2a_n_ions,2)
        epe(:reg_2a_n_ions)%s(3)=r_sh_ion(1:reg_2a_n_ions,3)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)

  allocate(r_in_pg_system(3,reg_2a_n_ions),stat=istat)
  if(istat.ne.0) call error_handler("allocate r_in_pg_system failed")
  !print*,' !**   start pcs treatment'
  r_in_pg_system(1,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,1)/auangs
  r_in_pg_system(2,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,2)/auangs
  r_in_pg_system(3,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,3)/auangs
  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)

  reg_I_pg(1:end_treated_region)%rs(1)=r_in_pg_system(1,1:end_treated_region)
  reg_I_pg(1:end_treated_region)%rs(2)=r_in_pg_system(2,1:end_treated_region)
  reg_I_pg(1:end_treated_region)%rs(3)=r_in_pg_system(3,1:end_treated_region)

  if(pc_aswritemode) then
     inquire(file=trim(epe_input_dir)//'/epe.pcs',exist=stat)

     pcs_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcs', &
          status='unknown',form='formatted')

     if(stat .and. .not.operations_make_reg_reference) then
        read(pcs_unit,'(i5)') nn_tmp
        allocate(q_tmp(nn_tmp))
        do i=1,nn_tmp
           read(pcs_unit,'(51x,f17.8)') q_tmp(i)
        enddo
        rewind pcs_unit
     endif

     ! their positions under optimization
     write(pcs_unit,'(i5)') up_limit-n_vacancies
     do i=1+n_vacancies,up_limit
       if(stat .and. .not.operations_make_reg_reference) then
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_tmp(i-n_vacancies) &
                ,1,0,0,1,0, an_epe_type(epe(i)%k)
        else
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,out_charges(epe(i)%k)-q_nuclear(epe(i)%k)/qau_qepe &
                ,1,0,0,1,0, an_epe_type(epe(i)%k)
        endif
     enddo
     DPRINT 'epe2pg pc_aswritemode',an_epe_type(epe(1+n_vacancies)%k)
     call returnclose_iounit(pcs_unit)
     if(allocated(q_tmp)) deallocate(q_tmp)

  else !pcs_vnwritemode
  pcs_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcs', &
                                status='unknown',form='formatted')
                                ! their positions under optimization
 
  write(pcs_unit,'(i5)') up_limit-n_vacancies   ! number of ions which change
  do i=1+n_vacancies,up_limit
     ind_charges: if(independent_epe_charges) then

        if(pseudobond_approach) then
           if(psb_pc(i)%linkto.gt.0) then
              psb_pc(i)%rb=r_in_pg_system(1:3,i)
              write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                   0.0_r8_kind,1,0,0,1,epe(i)%k , an_epe_type(epe(i)%k)
           elseif(psb_pc(i)%linkto.lt.0) then
              qs_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_shell(epe(i)%k)
              do link_count=1,-psb_pc(i)%linkto
                qs_psbond= qs_psbond-&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x* &
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb/&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto
              end do
              
            write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                qs_psbond, 1,0,0,1,epe(i)%k , an_epe_type(epe(i)%k)
           else
              write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                   ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_shell(epe(i)%k) &
                   ,1,0,0,1,epe(i)%k, an_epe_type(epe(i)%k)

           end if

        else ! not pseudobond_approach
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                *q_shell(epe(i)%k) &
                ,0,0,0,0,epe(i)%k , an_epe_type(epe(i)%k)
        end if

     else ind_charges
        write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(:,i) &
             ,q_shell(epe(i)%k)/qau_qepe,1,0,0,1,epe(i)%k, an_epe_type(epe(i)%k)
     endif ind_charges

  enddo

  call returnclose_iounit(pcs_unit)
  endif

  !print*,' **make file with coordinates of cores positioned to the PG '
  ! **cluster coordinate system


  r_in_pg_system(1,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,1)/auangs
  r_in_pg_system(2,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,2)/auangs
  r_in_pg_system(3,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,3)/auangs

  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)
  reg_I_pg(:end_treated_region)%rc(1)=r_in_pg_system(1,:end_treated_region)
  reg_I_pg(:end_treated_region)%rc(2)=r_in_pg_system(2,:end_treated_region)
  reg_I_pg(:end_treated_region)%rc(3)=r_in_pg_system(3,:end_treated_region)

        
  ! **now write PG PC related to the epe  cores
  
  pcc_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcc', &
       status='unknown',form='formatted')
  write(pcc_unit,'(i5)') up_limit-n_vacancies 
  do i=1+n_vacancies,up_limit

     if(independent_epe_charges) then
        if(pseudobond_approach) then
           if(psb_pc(i)%linkto.gt.0) then
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                   0.0_r8_kind,0,0,0,0,0 
           elseif(psb_pc(i)%linkto.lt.0) then
              qc_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k)
              rc_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k) *r_in_pg_system(1:3,i)
              do link_count=1,-psb_pc(i)%linkto
                qc_psbond= qc_psbond+&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb*(1.0_r8_kind &
                   +psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto
                rc_psbond=rc_psbond+&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb*(1.0_r8_kind &
                   +psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto &
                   *(psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%rb+&
                   r_in_pg_system(1:3,i)* &
                   psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /(1.0_r8_kind+psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)
              end do
              rc_psbond=rc_psbond/qc_psbond
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') rc_psbond, qc_psbond ,1,0,0,1,0  

           else
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                   ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k) ,1,0,0,1,0   

           end if


        else
           write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                *q_nuclear(epe(i)%k) ,1,0,0,1,0           
        end if
     else
        write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
             ,q_nuclear(epe(i)%k)/qau_qepe,1,0,0,1,0
     end if

  enddo
 
  call returnclose_iounit(pcc_unit)


      r_in_pg_system(1,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(1)/auangs
      r_in_pg_system(2,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(2)/auangs
      r_in_pg_system(3,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(3)/auangs
  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)

        if(pc_aswritemode) then
           write(pcr_unit,'(i5)') up_limit-n_vacancies
           do i=1+n_vacancies,up_limit
              write(pcr_unit,'(4f17.8,I3,4I2)') r_in_pg_system(:,i),  &
                   out_charges(epe(i)%k),1,0,0,1,0
           enddo
        else
           if(.not.exist_pcr) then
              print* , '**** file epe.pcr not found and will be created'



              write(pcr_unit,'(i5)') up_limit-n_vacancies
              do i=1+n_vacancies,up_limit
                 if(independent_epe_charges) then
                    write(pcr_unit,'(4f17.8,I5,4I5)') r_in_pg_system(:,i),  &
                         q_epecl_coupl(epe(i)%k) ,0,0,0,i-n_vacancies,i

                 else
                    write(pcr_unit,'(4f17.8,I5,4I5)') r_in_pg_system(:,i),  &
                         (q_nuclear(epe(i)%k)+q_shell(epe(i)%k))/qau_qepe  & 
                         ,0,0,0,i-n_vacancies,i
                 end if
              enddo


           end if
        endif

  call returnclose_iounit(pcr_unit)

 deallocate(r_in_pg_system,stat=istat)
 if(istat.ne.0) call error_handler(" deallocate r_in_pg_system failed")
 if(pseudobond_approach) then
    deallocate(psb_pc,stat=istat)
    if(istat.ne.0) call error_handler("deallocate psb_pc failed")
 end if

end subroutine epe2pg_epe

subroutine epe2pg(imp_core_grad,imp_shell_grad)
! **define corrections to the energy and the potentials for all
! **atoms of cluster
! **write coordinates of the Hades centers to file for use in
! **PG code
  use   epecom_module
  use   culon_module

  real(kind=r8_kind),allocatable,dimension(:,:):: r_in_pg_system
  real(kind=r8_kind), intent(in), dimension(N_IMPURITIES,3):: imp_core_grad, &
       imp_shell_grad
  integer(kind=i4_kind) :: i,up_limit
  integer(kind=i4_kind) :: istat,nn_tmp
  real(kind=r8_kind),allocatable :: q_tmp(:)
  logical :: stat
  logical exist_pcr

   type psb_centers
      real(kind=r8_kind)::rb(3),qb,x
      integer(kind=i4_kind):: linkto,linkfrom(3)
   end type psb_centers
   type(psb_centers),allocatable ::psb_pc(:)
   integer(kind=i4_kind)::link_count
   real(kind=r8_kind) ::qs_psbond,rc_psbond(3),qc_psbond

 !print*,' **make backward transformation for atoms of cluster ' 
  call search_centre_and_rotation &
        (natoms,gx_centers,gma,shft_gto_to_epe,-1)
  gx_unit=openget_iounit(trim(epe_input_dir)//'/gxfile', &
                          status='old',form='formatted')
  DPRINT 'epe2pg: '//trim(epe_input_dir)//'/gxfile unit=',gx_unit
  do i=1,natoms+1
203 FORMAT((f7.2,3(2X,f13.7),2i4,2x,3I4,2X,3I4,2x,i3))
        if(gx_highprec) then
        write(gx_unit,"(f7.2,3(2X,f21.12),2i4,2x,3I4,2X,3I4,2x,i3)") &
             an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
        else
        write(gx_unit,203) an(i),gx_centers(:,i),indgx(1:8,i),impu(i)
        endif
  enddo !i=1,natoms+1 

  if(qm_interfaced_mode) then
    write(gx_unit,'(2f20.7,3i5)') etot_epe, etot_epe,81,1,0      
    do i=1,n_impurities
      write(gx_unit,'(i5,3f15.8)') I,(imp_core_grad(i,:)+imp_shell_grad(i,:))*auangs
    enddo !1,n_impurities 
  endif !abs(interfaced_mode).eq.1) 
  call returnclose_iounit(gx_unit)

  !create/treat files to impose PC field on the 
  !EPE embedded cluster

  if(extended_epe) then
!!$  write(output_epe,*) 'extended_epe---------------',extended_epe
     up_limit=reg_2a_n_ions
  else
     up_limit=reg_I_n_ions   
  end if

  if(reg_2a_treated) then
!!$     print*,'reg_2a_treated------------------',reg_2a_treated
     end_treated_region= reg_2a_n_ions
  else
     end_treated_region=reg_I_n_ions
  end if

  inquire(file=trim(epe_input_dir)//'/epe.pcr', exist=exist_pcr)

  if(exist_pcr.and.qm_interfaced_mode.and.operations_make_reg_reference) &
       then
     call system('rm '//trim(epe_input_dir)//'/epe.pcr')
     exist_pcr=.false.
   endif

  pcr_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcr', &
                                status='unknown',form='formatted')
 !print*,'pseudobond_approach',  pseudobond_approach
  if(pseudobond_approach) then
     allocate(psb_pc(reg_I_n_ions),stat=istat)
     if(istat.ne.0) then
        call error_handler("allocate psb_pc failed")
     end if
     psb_pc(:)%x=9.0_r8_kind
     psb_pc(:)%linkto=0_i4_kind
     psb_pc(:)%linkfrom(1)=0_i4_kind
     psb_pc(:)%linkfrom(2)=0_i4_kind
     psb_pc(:)%linkfrom(3)=0_i4_kind

     if(exist_pcr) then
        read(pcr_unit,*) i
        if(i.ne.up_limit-n_vacancies) then
           print*,' inconsistent number in epe.pcr file'
        end if
        do i=n_vacancies+1,reg_I_n_ions
           read(pcr_unit,*)  &
                psb_pc(i)%rb ,psb_pc(i)%qb ,psb_pc(i)%linkto,psb_pc(i)%linkfrom
        end do
     end if

  end if
  
        epe(:reg_2a_n_ions)%s(1)=r_sh_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%s(2)=r_sh_ion(1:reg_2a_n_ions,2)
        epe(:reg_2a_n_ions)%s(3)=r_sh_ion(1:reg_2a_n_ions,3)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)
        epe(:reg_2a_n_ions)%c(1)=r_nuc_ion(1:reg_2a_n_ions,1)

  allocate(r_in_pg_system(3,reg_2a_n_ions),stat=istat)
  if(istat.ne.0) call error_handler("allocate r_in_pg_system failed")
 !print*,' !**   start pcs treatment'
  r_in_pg_system(1,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,1)/auangs
  r_in_pg_system(2,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,2)/auangs
  r_in_pg_system(3,1:reg_2a_n_ions)=r_sh_ion(1:reg_2a_n_ions,3)/auangs
  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)

  reg_I_pg(1:end_treated_region)%rs(1)=r_in_pg_system(1,1:end_treated_region)
  reg_I_pg(1:end_treated_region)%rs(2)=r_in_pg_system(2,1:end_treated_region)
  reg_I_pg(1:end_treated_region)%rs(3)=r_in_pg_system(3,1:end_treated_region)

  if(pc_aswritemode) then
     inquire(file=trim(epe_input_dir)//'/epe.pcs',exist=stat)

     pcs_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcs', &
          status='unknown',form='formatted')

     if(stat .and. .not.operations_make_reg_reference) then
        read(pcs_unit,'(i5)') nn_tmp
        allocate(q_tmp(nn_tmp))
        do i=1,nn_tmp
           read(pcs_unit,'(51x,f17.8)') q_tmp(i)
        enddo
        rewind pcs_unit
     endif

     ! their positions under optimization
     write(pcs_unit,'(i5)') up_limit-n_vacancies
     do i=1+n_vacancies,up_limit
       if(stat .and. .not.operations_make_reg_reference) then
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_tmp(i-n_vacancies) &
                ,1,0,0,1,0, an_epe_type(epe(i)%k)
        else
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,out_charges(epe(i)%k)-q_nuclear(epe(i)%k)/qau_qepe &
                ,1,0,0,1,0, an_epe_type(epe(i)%k)
        endif
     enddo
     DPRINT 'epe2pg pc_aswritemode',an_epe_type(epe(1+n_vacancies)%k)
     call returnclose_iounit(pcs_unit)
     if(allocated(q_tmp)) deallocate(q_tmp)

  else !pcs_vnwritemode
  pcs_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcs', &
                                status='unknown',form='formatted')
                                ! their positions under optimization
 
  write(pcs_unit,'(i5)') up_limit-n_vacancies   ! number of ions which change
  do i=1+n_vacancies,up_limit
     ind_charges: if(independent_epe_charges) then

        if(pseudobond_approach) then
           if(psb_pc(i)%linkto.gt.0) then
              psb_pc(i)%rb=r_in_pg_system(1:3,i)
              write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                   0.0_r8_kind,1,0,0,1,epe(i)%k , an_epe_type(epe(i)%k)
           elseif(psb_pc(i)%linkto.lt.0) then
              qs_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_shell(epe(i)%k)
              do link_count=1,-psb_pc(i)%linkto
                qs_psbond= qs_psbond-&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x* &
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb/&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto
              end do
              
            write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                qs_psbond, 1,0,0,1,epe(i)%k , an_epe_type(epe(i)%k)
           else
              write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                   ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_shell(epe(i)%k) &
                   ,1,0,0,1,epe(i)%k, an_epe_type(epe(i)%k)

           end if

        else ! not pseudobond_approach
           write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                *q_shell(epe(i)%k) &
                ,0,0,0,0,epe(i)%k , an_epe_type(epe(i)%k)
        end if

     else ind_charges
        write(pcs_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(:,i) &
             ,q_shell(epe(i)%k)/qau_qepe,1,0,0,1,epe(i)%k, an_epe_type(epe(i)%k)
     endif ind_charges

  enddo

  call returnclose_iounit(pcs_unit)
  endif

 !print*,' **make file with coordinates of cores positioned to the PG '
! **cluster coordinate system


  r_in_pg_system(1,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,1)/auangs
  r_in_pg_system(2,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,2)/auangs
  r_in_pg_system(3,1:reg_2a_n_ions)=r_nuc_ion(1:reg_2a_n_ions,3)/auangs

  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)
  reg_I_pg(:end_treated_region)%rc(1)=r_in_pg_system(1,:end_treated_region)
  reg_I_pg(:end_treated_region)%rc(2)=r_in_pg_system(2,:end_treated_region)
  reg_I_pg(:end_treated_region)%rc(3)=r_in_pg_system(3,:end_treated_region)

        
  ! **now write PG PC related to the epe  cores
  
  pcc_unit=openget_iounit(trim(epe_input_dir)//'/epe.pcc', &
       status='unknown',form='formatted')
  write(pcc_unit,'(i5)') up_limit-n_vacancies 
  do i=1+n_vacancies,up_limit

     if(independent_epe_charges) then
        if(pseudobond_approach) then
           if(psb_pc(i)%linkto.gt.0) then
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i),  &
                   0.0_r8_kind,0,0,0,0,0 
           elseif(psb_pc(i)%linkto.lt.0) then
              qc_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k)
              rc_psbond=&
                   q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k) *r_in_pg_system(1:3,i)
              do link_count=1,-psb_pc(i)%linkto
                qc_psbond= qc_psbond+&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb*(1.0_r8_kind &
                   +psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto
                rc_psbond=rc_psbond+&
                     psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%qb*(1.0_r8_kind &
                   +psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%linkto &
                   *(psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%rb+&
                   r_in_pg_system(1:3,i)* &
                   psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)&
                   /(1.0_r8_kind+psb_pc(psb_pc(i)%linkfrom(link_count)+n_vacancies)%x)
              end do
              rc_psbond=rc_psbond/qc_psbond
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') rc_psbond, qc_psbond ,1,0,0,1,0  

           else
              write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                   ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                   *q_nuclear(epe(i)%k) ,1,0,0,1,0   

           end if


        else
           write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
                ,q_epecl_coupl(epe(i)%k)/(q_nuclear(epe(i)%k)+q_shell(epe(i)%k))&
                *q_nuclear(epe(i)%k) ,1,0,0,1,0           
        end if
     else
        write(pcc_unit,'(4f17.8,I3,4I2,f6.2)') r_in_pg_system(1:3,i)  &
             ,q_nuclear(epe(i)%k)/qau_qepe,1,0,0,1,0
     end if

  enddo
 
  call returnclose_iounit(pcc_unit)


      r_in_pg_system(1,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(1)/auangs
      r_in_pg_system(2,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(2)/auangs
      r_in_pg_system(3,1:reg_2a_n_ions)=epe(1:reg_2a_n_ions)%r(3)/auangs
  call search_centre_and_rotation &
       (reg_2a_n_ions,r_in_pg_system,gma,shft_gto_to_epe,-1)

        if(pc_aswritemode) then
           write(pcr_unit,'(i5)') up_limit-n_vacancies
           do i=1+n_vacancies,up_limit
              write(pcr_unit,'(4f17.8,I3,4I2)') r_in_pg_system(:,i),  &
                   out_charges(epe(i)%k),1,0,0,1,0
           enddo
        else
           if(.not.exist_pcr) then
              print* , '**** file epe.pcr not found and will be created'



              write(pcr_unit,'(i5)') up_limit-n_vacancies
              do i=1+n_vacancies,up_limit
                 if(independent_epe_charges) then
                    write(pcr_unit,'(4f17.8,I5,4I5)') r_in_pg_system(:,i),  &
                         q_epecl_coupl(epe(i)%k) ,0,0,0,i-n_vacancies,i

                 else
                    write(pcr_unit,'(4f17.8,I5,4I5)') r_in_pg_system(:,i),  &
                         (q_nuclear(epe(i)%k)+q_shell(epe(i)%k))/qau_qepe  & 
                         ,0,0,0,i-n_vacancies,i
                 end if
              enddo


           end if
        endif

  call returnclose_iounit(pcr_unit)

 deallocate(r_in_pg_system,stat=istat)
 if(istat.ne.0) call error_handler(" deallocate r_in_pg_system failed")
 if(pseudobond_approach) then
    deallocate(psb_pc,stat=istat)
    if(istat.ne.0) call error_handler("deallocate psb_pc failed")
 end if

end subroutine epe2pg
!---------------------------------------------------------

SUBROUTINE search_centre_and_rotation(NN,ABC,rot,shift,is)

  real(kind=r8_kind), intent(inout), DIMENSION(3,*)::ABC 
  real(kind=r8_kind), intent(in),dimension(3,3):: rot
  real(kind=r8_kind), intent(in),dimension(3):: shift
  real(kind=r8_kind) :: RW(3)
  integer(kind=i4_kind) :: NN,is,j,n

  DO N=1,NN
    if(is.lt.0) then
      abc(1:3,n)=abc(1:3,n)-shift(1:3)
      rw(1:3)=0.d0
    else
      rw(1:3)=shift(1:3)
    endif
    DO J=1,3
        RW(J)=RW(J)+sum(ABC(:,N)*rot(J,:))
    enddo       !       J=1,3                   
      ABC(1:3,N)=RW(1:3)
  enddo ! N=1,NN

END SUBROUTINE search_centre_and_rotation
end module epe_pg_module

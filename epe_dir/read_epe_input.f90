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
subroutine read_epe_input
  ! **this procedure reads input file
# include <def.h>   
  use type_module
  use epecom_module
  use mol_module
  use str_module
  use atoms_parameters_module
  use epe_pg_module,only:an_epe_type,operations_make_reg_reference
  use iounitadmin_module
  use epepar_module, only: list_epe_par, treat_epepar_namelist &
       ,par, ec_name_of_type=>name_of_type, get_slsp, &
       ec_max_type_ions=>max_type_ions

  implicit none

  integer(kind=i4_kind) :: nr,r0im,r1im,kim,k1im
  real(kind=r8_kind):: mn,mx
  real(kind=r8_kind), dimension(3):: rl,xyz_center_gener,xyz_step_of_shift
  real(kind=r8_kind), dimension(3,ngxat) :: r_ion_in_cell_gx
  real(kind=r8_kind) :: Q_imp,Q_shell_imp,PK_imp,Rimp(3)
  real(kind=r8_kind) :: bim,roim,cim,dim
  logical :: next,lgxcell
  integer(kind=i4_kind) :: type_imp
  character(len=3) :: name_of_imp,imp_name(2)
  integer(kind=i4_kind) :: counter
  character(len=4), dimension(18):: comments
  integer(kind=i4_kind), dimension(ndt):: ion_types,position
  integer(kind=i4_kind) :: i,ii,jj,j,k,ix,li,lj,na,l,m, &
       STEP_OF_DISPLACEMENT=200,mti,output_level=3
  integer(kind=i4_kind) :: status
  logical :: epein,options_read_parameters
  logical :: lcell_a=.false. ! translation cell in units of A
  real(kind=r8_kind), dimension(3):: xx, yy, zz
  real(kind=r8_kind):: q_ml,r_ml(3),r_first_sphere=5.0,r_2A_sphere=15.0 &
       ,r_first_short_interact=10.0,r_second_short_interact=10.0
  character(len=120) :: q_name
  character(len=20) :: q_access, q_position
  integer(kind=i4_kind) :: n_3_body,i1,i2,i3
  integer(kind=i4_kind), allocatable :: n_3b(:),index(:,:,:)
  character(len=3) :: atm_nm(6,3)
  real(kind=r8_kind) :: k_i(6),r_3b,theta_i(6),ec_k_i(6),ec_theta_i(6)
  integer(kind=i4_kind):: interfaced_mode = 0
  character(len=80) :: message

  namelist /rotation_shift/ xx, yy, zz, shft_gto_to_epe, shft_gto_to_hds
  namelist /n_gener_ions/ n_gen_ions
  namelist /n_ions_in_primitive_cell/ n_ions_cell
  namelist /n_translation_of_primitive_cell/ n_trans_primitive_cell
  namelist /parameters_of_generation/ n_centres_of_generation
  namelist /centres_of_generation/ xyz_center_gener,xyz_step_of_shift
  namelist /main_parameters/ radius_first_sphere,radius_2A_sphere, &
           radius_long_interact,error_function_parameter
  namelist /vacancy_impurity/ n_vacancies,n_impurities,n_types_central_atoms_3body_im
  namelist /impurities/ name_of_imp,type_imp,Q_imp,Q_shell_imp,PK_imp,rimp

  namelist /work_options/ periodic_optimization,use_pgdata, &
                          lcell_a,interfaced_mode, qm_interfaced_mode, &
                          basic_action, &
                          output_level, &
                          read_configuration,write_configuration, &
                          write_config_to_lcgto_format,print_deform, &
                          print_ions_from_first_sphere, &
                          n_print_additional_ions,step_of_displacement, &
	                  an_epe_type,options_read_parameters,scaling_factor, &
			  make_epe_reference, operations_make_reg_reference,&
			  use_epe_reference,n_pgepe_iterations, &
                          extended_epe, &
			  option_c3_symm, &
			  ml_displacements, &
			  fixed_dielectric_const, fx_dielectric_const, &
                          ml_cluster_simulated, n_ml_cluster_simulators, &
			  ml_tensors, &
			  GX_HIGHPREC , &
                          error_function_parameter, &
                          r_first_sphere, r_2A_sphere &
                          ,r_first_short_interact,r_second_short_interact , &
                          embed_convergence_check, &
                          embed_convergence_limit, & ! energy gain in epe optimization
                          n_types_central_atoms_3body, &
                          abs_g, &
			  weight_hess, &
			  n_hess_update_cycles, &
                          independent_epe_charges, &
                          pseudobond_approach, &
			  explicit_coupling

  namelist /epe_parameters/ parv,options_read_parameters
  namelist /imp_param_of_potential/ imp_name,bim,roim,cim,dim,next
  namelist /cluster_simulators/ q_ml,r_ml
  namelist /tree_body_interaction/ n_3_body,atm_nm,k_i,r_3b,theta_i,ec_k_i,ec_theta_i
  namelist /optimization_parameters/ n_iterations,abs_g,n_hess_update_cycles, &
            print_gradients,weight_hess, n_pgepe_iterations
  namelist /xyz_output/ unit_cell_output,region_I_output,output_step


  print*, '****** READ_EPE_INPUT ******', trim(epe_input_dir)
  core_shell=.false.
  input_epe=openget_iounit(trim(epe_input_dir)//'/epe.input', &
				 form='formatted', status='old')

  call init_epe_par()
	
  output_epe=openget_iounit(trim(epe_input_dir)//'/epe.out', form='formatted',&
                                            status='unknown', position='append')
 
  options_read_parameters=.false.

  inquire(output_epe, name=q_name, access=q_access, position=q_position)

  scaling_factor=1.889726877774_r8_kind*auangs
  radius_long_interact= 100.0_r8_kind
  r_2A_sphere=18.0_r8_kind
  r_first_sphere=5.0_r8_kind
  r_first_short_interact=10.0_r8_kind
  r_second_short_interact=10.0_r8_kind
  ERROR_FUNCTION_PARAMETER=df_ERROR_FUNCTION_PARAMETER

  read (input_epe,nml=work_options)
	write(output_epe,nml=work_options)
	
	if(interfaced_mode.ne.0) qm_interfaced_mode=.true.
  reg_2a_treated=ml_displacements.and..not.ml_cluster_simulated
  if(operations_make_reg_reference) then
     make_epe_reference=.false.
     use_epe_reference=.false.   
     use_pgdata = .false.
  end if
  if(make_epe_reference) embed_convergence_check=.false.
  
  DPRINT 'write work_options to output_epe',output_epe
  write(output_epe,nml=work_options)
  DPRINT 'done'
  
  DPRINT 'read nml=optimization_parameters'
  read (input_epe,nml=optimization_parameters)
  DPRINT 'n_iterations,abs_g,n_hess_update_cycles,weight_hess,print_gradients,n_pgepe_iterations'
  DPRINT n_iterations,abs_g,n_hess_update_cycles,weight_hess,print_gradients,n_pgepe_iterations
  write(output_epe,nml=optimization_parameters)
  DPRINT 'done optimization_parameters'

  read (input_epe,nml=xyz_output)
  
  write(output_epe,nml=xyz_output)
  DPRINT 'write nml=xyz_output done'

  radius_first_sphere=r_first_sphere 
  radius_2A_sphere=r_2A_sphere

  if(ml_cluster_simulated.and.n_ml_cluster_simulators.ne.0) then
     allocate (ml_cluster_simulators(n_ml_cluster_simulators),stat=status)
     if(status.ne.0) call error_handler("allocate ml_cluster_simulators failed")
     do counter=1,n_ml_cluster_simulators
        read (input_epe,nml=cluster_simulators)
        ml_cluster_simulators(counter)%q=q_ml*qau_qepe
        ml_cluster_simulators(counter)%r=r_ml*auangs
     end do
  end if

  call set_default_parv
  DPRINT 'set_default_parv done'

  do while(options_read_parameters)
     options_read_parameters=.false.
     read (input_epe,nml=epe_parameters)
     allocate(p)
     p=epe_par(parv,top)
     top=>p
  enddo ! while


  read (input_epe,nml=rotation_shift)
  write (output_epe,nml=rotation_shift)

  rot_gto_to_epe(:,1)=xx
  rot_gto_to_epe(:,2)=yy
  rot_gto_to_epe(:,3)=zz

  read(input_epe,nml=n_gener_ions)
  write(output_epe, nml=n_gener_ions)
  allocate(r_nuc_ion(n_gen_ions,3),r_sh_ion(n_gen_ions,3),epe(n_gen_ions), &
                                     q_zl(n_gen_ions),stat=epealloc_stat(9))
           ASSERT(epealloc_stat(9).eq.0)
                  epealloc_stat(9)=1
                 epealloc_stat(10)=1 ! q_zl

101 format(18a4) 
  read(input_epe, 101)COMMENTS
  write(output_epe, 101)COMMENTS

  DO I=1,3
    read(input_epe,*)(VECTORS_TRANS(I,J),J=1,3)
    write(output_epe,*)(VECTORS_TRANS(I,J),J=1,3)
  enddo	! I=1,3

  DPRINT 'read nml=n_ions_in_primitive_cell'
  read(input_epe,nml=n_ions_in_primitive_cell)
   write(output_epe, nml=n_ions_in_primitive_cell)
   allocate(which_epe_ion(n_ions_cell),stat=epealloc_stat(7))
   ASSERT(epealloc_stat(7).eq.0)
   epealloc_stat(7)=1

  read(input_epe, 101)COMMENTS
  write(output_epe, 101)COMMENTS

  fixed(1:N_IONS_CELL)=.false.

  max_type_ions=0

  do I=1,N_IONS_CELL
     read(input_epe,'(a80)')message
     READ (message,*,end=32) NA,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:),fixed(i)
     goto 30
32   READ (message,*) NA,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
30   if(fixed(i))  write(output_epe,1016) NA,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:),fixed(i)
     if(.not.fixed(i)) write(output_epe,1006) NA,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
     ion_types(TYPE_OF_ION(I))=type_of_ion(i)
     name_of_type(TYPE_OF_ION(I))=NAME_OF_IONS(i)
     if(max_type_ions.lt.TYPE_OF_ION(I)) max_type_ions=TYPE_OF_ION(I)
  enddo! I=1,N_IONS_CELL
1006 FORMAT (I3,I2,1X,A2,3F10.5)
1016 FORMAT (I3,I2,1X,A3,3F10.5,1x,l1)

  allocate(ml_dprime_factors(max_type_ions),ml_fac(max_type_ions),stat=epealloc_stat(8))
  ASSERT(epealloc_stat(8).eq.0)
         epealloc_stat(8)=1
	write(output_epe,*) 'explicit_coupling',explicit_coupling
	if(explicit_coupling) then
	write(output_epe,*) 'call treat_epepar_namelist'
	 call treat_epepar_namelist()
	 call list_epe_par
	 print*,'done'
	 call list_epe_par
	endif

  if(lcell_a) then
     write(output_epe,*) ' translation cell in units of A'
     do i=1,n_ions_cell
        rl(:)=r_ion_in_cell(i,:)
        r_ion_in_cell(i,1)=dot_product(rl,vectors_trans(1,:))
        r_ion_in_cell(i,2)=dot_product(rl,vectors_trans(2,:))
        r_ion_in_cell(i,3)=dot_product(rl,vectors_trans(3,:))
        write (output_epe,1006) i,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
     enddo	! i=1,n_ions_cell
  endif ! lcell_a

  VECTORS_TRANS=VECTORS_TRANS*scaling_factor
  DPRINT 'inquire gxcell file'
  inquire (file=trim(epe_input_dir)//'/gxcell',exist=lgxcell)
  
  if(lgxcell) then
     write(output_epe,*) 'read coordinates of cell from gx.c file'
     print*,'old inp: read coordinates of cell from gx.c file'
     gxcell_unit=openget_iounit(trim(epe_input_dir)//'/gxcell',  &
          form='formatted', status='unknown')
     ii=0
2000 ii=ii+1
     if(ii.gt.ngxat) then
        write(output_epe,*) 'read_input: ngxat is too small'
     else	! ii.le.ngxat
        read(gxcell_unit,*)an(ii),r_ion_in_cell_gx(:,ii),indgx(:,ii),impu(ii)
        write(output_epe,1006) ii,type_of_ion(ii),NAME_OF_IONS(ii)&
             ,r_ion_in_cell_gx(:,ii)
        if(an(ii).gt.0.5d0) go to 2000
        write(output_epe,*) ' read_input: No of centers in gx.c file',ii-1
        jj=0
        do i=1,ii-1
           if(indgx(1,i).ne.0) then
              jj=jj+1
              r_ion_in_cell(jj,:)=r_ion_in_cell_gx(:,i)	
           endif	! indgx(1,i).ne.0
        enddo	!i=1,ii-1 
     endif	! ii.gt.ngxat
    
     print*,'gx.c file read'
     rewind gxcell_unit
     do i=1,ii
        write(gxcell_unit,203) an(i),r_ion_in_cell_gx(:,i),indgx(:,i),impu(i)
203     FORMAT((f5.2,3(2X,f13.7),2i4,2x,3I3,2X,3I3,2x,2i3))
     enddo !i=1,ii 
     print*,'gx.c file is updated'
	
     if(.not.periodic_optimization) call returnclose_iounit(gxcell_unit) !nothing to read or write
  else	! not.lgxcell fill in indgx matrix
     if(periodic_optimization) then	
        gxcell_unit=openget_iounit(trim(epe_input_dir)//'/gxcell', &
             form='formatted', status='unknown')
        do i=1,n_ions_cell
           indgx(1,i)=i
        enddo	! i=1,n_ions_cell
        do i=1,n_ions_cell
           an(i)=TYPE_OF_ION(i)
           write(gxcell_unit,203) an(i),(r_ion_in_cell(i,j),j=1,3)  &
                ,(indgx(k,i),k=1,8),impu(i)
        enddo	! i=1,n_ions_cell
        write(gxcell_unit,203) zero,zero,zero,zero,0,0,0,0,0,0,0,0,0
     endif	! periodic_optimization
  endif	! lgxcell

  R_ION_IN_CELL(1:n_ions_cell,:)=R_ION_IN_CELL(1:n_ions_cell,:)*scaling_factor
  write(output_epe,*)	' transcell in atomic units '
  do i=1,n_ions_cell
     write(output_epe,203) an(i),(r_ion_in_cell(i,j)/auangs,j=1,3)  &
          ,i,indgx(2:8,i)
  enddo	! i=1,n_ions_cell 

  if (n_types_central_atoms_3body /= 0) then
     atm_nm='  '
     allocate(ki(max_type_ions,max_type_ions,max_type_ions), &
              theta_0(max_type_ions,max_type_ions,max_type_ions),stat=epealloc_stat(13))
              ASSERT(epealloc_stat(13).eq.0)
                     epealloc_stat(13)=1
     if(explicit_coupling) then
      allocate(ec%ki(max_type_ions,max_type_ions,max_type_ions), &
              ec%theta_0(max_type_ions,max_type_ions,max_type_ions), &
                                               stat=epealloc_stat(16))
              ASSERT(epealloc_stat(16).eq.0)
                     epealloc_stat(16)=1
     endif
     ki=0.0_r8_kind
     allocate(types(n_types_central_atoms_3body,5),r3b(n_types_central_atoms_3body), &
          n_3b(n_types_central_atoms_3body),index(n_types_central_atoms_3body,6,3), &
          stat=epealloc_stat(14))
          ASSERT(epealloc_stat(14).eq.0)
                 epealloc_stat(14)=1
                 epealloc_stat(18)=1 ! n_3b index

     types=0
     do i=1,n_types_central_atoms_3body
        read(input_epe,nml=tree_body_interaction) ! 2
!!$        write (output_epe,nml=tree_body_interaction)
        n_3b(i)=n_3_body
        do j=1,n_3_body
           do k=1,max_type_ions
              if(atm_nm(j,1) == name_of_type(k)) i1=k
              if(atm_nm(j,2) == name_of_type(k)) i2=k
              if(atm_nm(j,3) == name_of_type(k)) i3=k
           enddo
           ki(i1,i2,i3)=k_i(j)
           ki(i3,i2,i1)=k_i(j)
           theta_0(i1,i2,i3)=theta_i(j)
           theta_0(i3,i2,i1)=theta_i(j)
         if(explicit_coupling) then
           ec%ki(i1,i2,i3)=k_i(j)
           ec%ki(i3,i2,i1)=k_i(j)
           ec%theta_0(i1,i2,i3)=theta_i(j)
           ec%theta_0(i3,i2,i1)=theta_i(j)
         endif
           index(i,j,1)=i1
           index(i,j,2)=i2
           index(i,j,3)=i3
           if(j==1) then
              types(i,1)=i2
           endif
           do l=2,5
              if(types(i,l)==i1) exit
              if(types(i,l)==0) then
                 types(i,l)=i1
                 exit
              endif
           enddo
           do l=2,5
              if(types(i,l)==i3) exit
              if(types(i,l)==0) then
                 types(i,l)=i3
                 exit
              endif
           enddo
        enddo
        r3b(i)=r_3b
     enddo
  endif

	do i=1,n_ions_cell
	do j=i+1,n_ions_cell
	 if(dot_product(R_ION_IN_CELL(i,:)-R_ION_IN_CELL(j,:), &
			R_ION_IN_CELL(i,:)-R_ION_IN_CELL(j,:)) &
		.lt.0.001_r8_kind) then
	print*,'atoms in unit cell coincide',i,j
	stop
	endif
	enddo
	enddo
	
  DPRINT 'read nml=n_translation_of_primitive_cell'
  read (input_epe,nml=n_translation_of_primitive_cell)
  write (output_epe,nml=n_translation_of_primitive_cell)

  read (input_epe,nml=parameters_of_generation)
  write (output_epe,nml=parameters_of_generation)

  if(n_centres_of_generation.ne.0) then
     DO I=1,N_CENTRES_OF_GENERATION
        read (input_epe,nml=centres_of_generation)
        write (output_epe,nml=centres_of_generation)
        R_CENT_GENER(I,:)=xyz_center_gener(:)*scaling_factor
        shift_cent_gener(i,:)=xyz_step_of_shift(:)*scaling_factor
     enddo ! I=1,N_CENTRES_OF_GENERATION
  endif

  read (input_epe,nml=main_parameters)
  write (output_epe,nml=main_parameters)

  RADIUS_FIRST_SPHERE=RADIUS_FIRST_SPHERE*scaling_factor
  RADIUS_2A_SPHERE=RADIUS_2A_SPHERE*scaling_factor
  RADIUS_LONG_INTERACT=(RADIUS_LONG_INTERACT*scaling_factor)**2
  
  if(explicit_coupling) then
	i=max(max_type_ions,ec_max_type_ions)
     else
	i=max_type_ions
   endif

  allocate(host%sr1(i,i,i+1), host%b(i,i,i+1),host%ro(i,i,i+1),&
       host%c(i,i,i+1),host%d(i,i,i+1),host%sr2(i,i,i+1),&
       host%k(i,i,i+1),host%r0(i,i,i+1),host%k1(i,i,i+1),host%r1(i,i,i+1), & 
       stat=epealloc_stat(1))
  if(epealloc_stat(1).ne.0) call error_handler(" 1 allocate host failed")
  epealloc_stat(1)=1 ! host%sr1 host%sr2
  epealloc_stat(2)=1 ! host% b ro c d 
  host%k=0.0_r8_kind 
  host%r0=0.0_r8_kind 
  host%k1=0.0_r8_kind 
  host%r1=0.0_r8_kind 
  host%b=0.0_r8_kind 
  host%ro=0.0_r8_kind; host%ro(:,:,1)=1.0_r8_kind     !!!!!!!!!!!!!AS
  host%c=0.0_r8_kind; host%d=0.0_r8_kind
  host%sr1=10.0_r8_kind
  host%sr2=10.0_r8_kind
   
  if(explicit_coupling) then
    do i=1,ec_max_type_ions
     do j=i,ec_max_type_ions
	mn=min(ec_name_of_type(i),ec_name_of_type(j))
	mx=max(ec_name_of_type(i),ec_name_of_type(j))
	call get_slsp(mn,mx,0,nr)
	if(nr.eq.0) then
	print*,'read_epe_input: no pair potential is found'
	print*,mn,mx
	else
	 host%k(i,j,1)=par%item%k
	 host%k(j,i,1)=par%item%k
	 host%k1(i,j,1)=par%item%k1
	 host%k1(j,i,1)=par%item%k1
	 host%r1(i,j,1)=par%item%r1
	 host%r1(j,i,1)=par%item%r1
	 host%r0(i,j,1)=par%item%r0
	 host%r0(j,i,1)=par%item%r0

	 host%b(i,j,1)=par%item%b
	 host%b(j,i,1)=par%item%b
	 host%ro(i,j,1)=par%item%r
	 host%ro(j,i,1)=par%item%r
	 host%c(j,i,1)=par%item%c
	 host%c(i,j,1)=par%item%c
	 host%d(j,i,1)=par%item%d
	 host%d(i,j,1)=par%item%d
	 host%sr1(j,i,1)=par%item%cutoff
	 host%sr1(i,j,1)=par%item%cutoff
	 host%sr2(j,i,1)=par%item%cutoff
	 host%sr2(i,j,1)=par%item%cutoff
	 
	endif
     enddo
    enddo
  endif
  do i=1,max_type_ions
     do j=i,max_type_ions
           position(i)=0
           position(j)=0
	p=>top
	do while(associated(p)) 
           position(i)=0
           position(j)=0
           do k=1,p%item%nlth
              if(name_of_type(i)==p%item%atom_name(k))  then
                 position(i)=k
                 exit
              endif ! 
           enddo ! k=1,p%item%nlth
           if(position(i).ne.0) then
              do k=1,p%item%nlth
                 if(name_of_type(j)==p%item%atom_name(k))  then
                    position(j)=k
                    exit
                 endif !
              enddo ! k=1,p%nlth
           endif ! position(i).ne.0
           if(position(i).ne.0.and.position(j).ne.0) exit
           p=>p%next
	enddo ! while
        if(position(i).eq.0.or.position(j).eq.0)  then
           print*,name_of_type(j)
           stop 'no epe parameters found'
           else
            print*,'search is finished' , position(i),  position(j)
	endif

	
        li=position(i)
        lj=position(j)
	write(output_epe,*)li,lj,p%item%bi(li,lj)	
!       no harmonic pptential input via epe.input
        print*,li,lj,p%item%bi(li,lj)
        host%b(i,j,1)=p%item%bi(li,lj)
        host%b(j,i,1)=p%item%bi(li,lj)
        host%ro(I,J,1)=p%item%roi(li,lj)
        host%ro(j,i,1)=p%item%roi(li,lj)
        host%c(i,j,1)=p%item%ci(li,lj)
        host%c(j,i,1)=p%item%ci(li,lj)
        host%d(i,j,1)=p%item%di(li,lj)
        host%d(j,i,1)=p%item%di(li,lj)
	 host%sr1(j,i,1)=p%item%cutoff(li,lj)
	 host%sr1(i,j,1)=p%item%cutoff(li,lj)
	 host%sr2(j,i,1)=p%item%cutoff(li,lj)
	 host%sr2(i,j,1)=p%item%cutoff(li,lj)
     enddo ! j=i,max_type_ions

     q_shell(i)=p%item%q_si(li)*qau_qepe
     q_nuclear(i)=(p%item%q_ni(li)-p%item%q_si(li))*qau_qepe
     pk(i)=p%item%pki(li)
     q_ion(i)=q_shell(i)+q_nuclear(i)
     q_epecl_coupl(i)=p%item%q_zi(i)
  enddo ! i=1,max_type_ions

  write(output_epe,*) '=========================================================='
  write(output_epe,*) '        The Potential Parameters of the Lattice'

  write(output_epe,*) '----------------------------------------------------------' 
  write(output_epe,*) '                 Two-body parameters'
  write(output_epe,*) '----------------------------------------------------------' 
  write(output_epe,*) '              b           ro          c           d' 

  if(explicit_coupling) then
	print*,'explicit_coupling  ec_max_type_ions ', ec_max_type_ions
  do i=1,ec_max_type_ions
     do j=i,ec_max_type_ions
        write(output_epe,'(2f6.2,4f12.5)') &
	ec_name_of_type(i),ec_name_of_type(j), &
	host%b(i,j,1),host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1)
       write(output_epe,'(12x,4f12.5)') &
       host%k(i,j,1),host%r0(i,j,1),host%k1(i,j,1),host%r1(i,j,1)
     enddo 
  enddo 
  endif

  do i=1,max_type_ions
     do j=i,max_type_ions
        write(output_epe,'(2a4,4f12.5)') name_of_type(i),name_of_type(j), &
	host%b(i,j,1),host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1)
       write(output_epe,'(12x,4f12.5)') &
       host%k(i,j,1),host%r0(i,j,1),host%k1(i,j,1),host%r1(i,j,1)
     enddo ! j=i,max_type_ions
  enddo ! i=1,max_type_ions
  
  write(output_epe,*) '-------------------------------------------'
  write(output_epe,*) '        Shell-model parameters'
  write(output_epe,*) '-------------------------------------------'
  write(output_epe,*) '      q_nuclear     q_shell       pk'
  do i=1,max_type_ions
     write(output_epe,'(a4,4f12.5)') name_of_type(i),&
          q_nuclear(i)/qau_qepe,q_shell(i)/qau_qepe,pk(i)
  enddo ! i=1,max_type_ions

  if (n_types_central_atoms_3body /= 0) then
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '              Three-body parameters'
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '                  ki        theta       r3b'
     do i=1,n_types_central_atoms_3body
        do j=1,n_3b(i)
           write(output_epe,'(3a4,3f12.5)') name_of_type(index(i,j,1)), &
                name_of_type(index(i,j,2)),name_of_type(index(i,j,3)), &
                ki(index(i,j,1),index(i,j,2),index(i,j,3)), &
                theta_0(index(i,j,1),index(i,j,2),index(i,j,3)), r3b(i)
        enddo
     enddo
     deallocate(n_3b,index)
  endif
  write(output_epe,*) '=========================================================='

  do i=1,n_ions_cell
     do j=1,max_type_ions
        if(name_of_ions(i)==name_of_type(j)) then
           q_z(i)=q_ion(j)/qau_qepe
           exit
        endif
     enddo
  enddo
  if(independent_epe_charges) then
print*,' independent_epe_charges'
     write(output_epe,*) & 
          '    q_nuclear     q_shell       pk   q_epecl'
  else
     write(output_epe,*) '    q_nuclear     q_shell       pk'
  end if

  do i=1,max_type_ions
     if(independent_epe_charges) then 
        write(output_epe,'(i2,5f12.5)') i,&
             q_nuclear(i)/qau_qepe,q_shell(i)/qau_qepe,pk(i),q_epecl_coupl(i)  
     else
        write(output_epe,'(i2,4f12.5)') i,&
             q_nuclear(i)/qau_qepe,q_shell(i)/qau_qepe,pk(i)
     end if
 
  enddo ! i=1,max_type_ions
  write(output_epe,*) '-------------------------------'

92  FORMAT(I5,6F9.3,3I3) !!!!!!!!!!!!!!!!!!!!


      do i=1,n_ions_cell
      do j=1,max_type_ions
       if(name_of_ions(i)==name_of_type(j)) then
       q_z(i)=q_ion(j)/qau_qepe
       exit
       endif
      enddo
      enddo

! **parameters for impurities   

  read (input_epe,nml=vacancy_impurity)
  write (output_epe,nml=vacancy_impurity)

  IF(N_IMPURITIES.EQ.0) then
     write(output_epe, *)" there are no basic impurities "

  else
     read(input_epe, 101)COMMENTS
     read(input_epe, 101)COMMENTS
     mti=max_type_ions
     DO I=1,N_IMPURITIES
        read (input_epe,nml=impurities)

        TYPE_IMPURITY(I)=type_imp
        if(TYPE_IMPURITY(I).le.mti) then
           write(output_epe,*) 'read_epe_input:type of impurity cannot be smaller then'
           write(output_epe,*) 'largest type of ions'
           call error_handler("read_epe_input:input error ")
        endif
        if(max_type_ions < TYPE_IMPURITY(I)) then
           max_type_ions=TYPE_IMPURITY(I)
           name_of_type(max_type_ions)=name_of_imp
           do j=1,max_type_ions-1
              if(name_of_type(max_type_ions) == name_of_type(j)) &
              call error_handler &
	("read_epe_input: Different types of ions must have different names")
           enddo
        endif
        Q_IMPURITY(I)=Qau_qepe*Q_imp
        Q_SH_IMPURITY(I)=Qau_qepe*q_shell_imp
        Q_NUC_IMPURITY(I)=Q_IMPURITY(I)-Q_SH_IMPURITY(I)
        pk_impurity(i)=PK_imp
        R_imp(I,:)=rimp(:)*scaling_factor
        R_SH_IMP(I,:)=R_imp(I,:)
        R_SH_IMPO(I,:)=R_imp(I,:)
        R_NUC_IMP(I,:)=R_imp(I,:)
        R_NUC_IMPO(I,:)=R_imp(I,:)
     enddo	!  I=1,N_IMPURITIES

     do i=1,max_type_ions
        k=mti+1
        if(k < i) k=i
        do j=k,max_type_ions
           host%k(J,I,1)=0.0_r8_kind
           host%k(I,J,1)=0.0_r8_kind
           host%k1(J,I,1)=0.0_r8_kind
           host%k1(I,J,1)=0.0_r8_kind
           host%r1(J,I,1)=0.0_r8_kind
           host%r1(I,J,1)=0.0_r8_kind
           host%r0(I,J,1)=0.0_r8_kind
           host%r0(J,I,1)=0.0_r8_kind

           host%B(J,I,1)=0.0_r8_kind
           host%B(I,J,1)=0.0_r8_kind
           host%Ro(J,I,1)=1.0_r8_kind
           host%RO(I,J,1)=1.0_r8_kind
           host%C(J,I,1)=0.0_r8_kind
           host%C(I,J,1)=0.0_r8_kind
           host%D(J,I,1)=0.0_r8_kind
           host%D(I,J,1)=0.0_r8_kind
           host%sr1(i,J,1)=3.5_r8_kind
           host%sr1(j,i,1)=3.5_r8_kind
           host%sr2(i,J,1)=3.5_r8_kind
           host%sr2(j,i,1)=3.5_r8_kind
        enddo
     enddo

     m=max_type_ions*(max_type_ions-1)/2+max_type_ions
     m=m-(mti*(mti-1)/2+mti)
     do i=1,m
        next=.true.

        kim=0.0_r8_kind
        k1im=0.0_r8_kind
        r1im=0.0_r8_kind
        r0im=0.0_r8_kind

        bim=0.0_r8_kind
        roim=1.0_r8_kind
        cim=0.0_r8_kind
        dim=0.0_r8_kind

        read (input_epe,nml=imp_param_of_potential)
        do j=1,max_type_ions
           if(imp_name(1) == name_of_type(j)) exit
        enddo
        do k=1,max_type_ions
           if(imp_name(2) == name_of_type(k)) exit
        enddo
        host%k(J,K,1)=kim
        host%k(K,J,1)=kim
        host%k1(J,K,1)=k1im
        host%k1(K,J,1)=k1im
        host%r1(J,K,1)=r1im
        host%r1(K,J,1)=r1im
        host%r0(J,K,1)=r0im
        host%r0(K,J,1)=r0im

        host%B(J,K,1)=bim
        host%B(K,J,1)=bim
        host%Ro(J,K,1)=roim
        host%RO(K,J,1)=roim
        host%C(J,K,1)=cim
        host%C(K,J,1)=cim
        host%D(J,K,1)=dim
        host%D(K,J,1)=dim
        if(.not. next) exit
     enddo

     write(output_epe,*) '=========================================================='
     write(output_epe,*) '    The Potential Parameters of the Impurities'
     write(output_epe,*) '----------------------------------------------------------'
     write(output_epe,*) '                 Two-body parameters'
     write(output_epe,*) '----------------------------------------------------------'
     write(output_epe,*) '            b_im        ro_im       c_im        d_im' 
     
     do i=1,max_type_ions
        k=mti+1
        if(k < i) k=i
        do j=k,max_type_ions
           write(output_epe,'(2a4,4f12.5)')name_of_type(i),name_of_type(j),host%b(i,j,1) &
	,host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1)
           write(output_epe,'(8x ,4f12.5)') host%k(i,j,1) &
	,host%r0(i,j,1),host%k1(i,j,1),host%r1(i,j,1)
        enddo
     enddo
     write(output_epe,*) '----------------------------------------------------------'

     write(output_epe,*) '        Shell-model parameters'
     write(output_epe,*) '-------------------------------------------'
     write(output_epe,*) '      q_nuclear     q_shell       pk'
     do i=mti,max_type_ions
        do j=1,n_impurities
           if(i == type_impurity(j)) then
              write(output_epe,'(a4,4f12.5)') name_of_type(i),&
                   q_nuc_impurity(j)/qau_qepe,q_sh_impurity(j)/qau_qepe,pk_impurity(j)
              exit
           endif
        enddo
     enddo
     
     if(n_types_central_atoms_3body_im > 0) then
        atm_nm='  '
        allocate(types_im(n_types_central_atoms_3body_im,5),r3b_im(n_types_central_atoms_3body_im), &
             n_3b(n_types_central_atoms_3body_im),index(n_types_central_atoms_3body_im,6,3), &
             stat=status)
        if(status.ne.0) call error_handler("allocate types and r3b failed")
        types_im=0
        do i=1,n_types_central_atoms_3body_im
           read (input_epe,nml=tree_body_interaction)
           n_3b(i)=n_3_body
           do j=1,n_3_body
              do k=1,max_type_ions
                 if(atm_nm(j,1) == name_of_type(k)) i1=k
                 if(atm_nm(j,2) == name_of_type(k)) i2=k
                 if(atm_nm(j,3) == name_of_type(k)) i3=k
              enddo
              ki(i1,i2,i3)=k_i(j)
              ki(i3,i2,i1)=k_i(j)
              theta_0(i1,i2,i3)=theta_i(j)
              theta_0(i3,i2,i1)=theta_i(j)
              index(i,j,1)=i1
              index(i,j,2)=i2
              index(i,j,3)=i3
              if(j==1) then
                 types_im(i,1)=i2
              endif
              do l=2,5
                 if(types_im(i,l)==i1) exit
                 if(types_im(i,l)==0) then
                    types_im(i,l)=i1
                    exit
                 endif
              enddo
              do l=2,5
                 if(types_im(i,l)==i3) exit
                 if(types_im(i,l)==0) then
                    types_im(i,l)=i3
                    exit
                 endif
              enddo
           enddo
           r3b_im(i)=r_3b
        enddo
     endif
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '              Three-body parameters'
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '                  ki        alpha       r3b'
     do i=1,n_types_central_atoms_3body_im
        do j=1,n_3b(i)
           write(output_epe,'(3a4,3f12.5)') name_of_type(index(i,j,1)), &
                name_of_type(index(i,j,2)),name_of_type(index(i,j,3)), &
                ki(index(i,j,1),index(i,j,2),index(i,j,3)), &
                theta_0(index(i,j,1),index(i,j,2),index(i,j,3)), r3b_im(i)
        enddo
     enddo
     deallocate(n_3b,index)

     write(output_epe,*) '=========================================================='
  endif

  if(n_centres_of_generation.eq.0.and.pg_interfaced_mode.eq.0) pg_interfaced_mode=1



  E2BIND=zero
  E2BDIS=zero
  E2BDEF=zero
  RECA=1.0_r8_kind/scaling_factor
  DMAX=scaling_factor/4.0_r8_kind
  epeit_count=0
  DEL=scaling_factor/STEP_OF_DISPLACEMENT
  ERROR_FUNCTION_PARAMETER=ERROR_FUNCTION_PARAMETER/4
  PIS=SQRT(PI)
  ET2=ERROR_FUNCTION_PARAMETER**2
  ERFO=-2.0_r8_kind*ERROR_FUNCTION_PARAMETER/PIS

  call returnclose_iounit(input_epe)

end subroutine read_epe_input
!******************************************************************************

!******************************************************************************
subroutine read_epe_input_new
  ! **this procedure read input file

  use type_module
  use epecom_module
  use mol_module
  use str_module
  use epe_pg_module,only:an_epe_type,operations_make_reg_reference
  use iounitadmin_module
  use epepar_module, only: treat_epepar_namelist &
        ,par, ec_name_of_type=>name_of_type, get_slsp, &
         ec_max_type_ions=>max_type_ions

  implicit none

  integer(kind=i4_kind) :: nr
  real(kind=r8_kind):: mn,mx
  real(kind=r8_kind), dimension(3):: rl,xyz_center_gener,xyz_step_of_shift
  real(kind=r8_kind), dimension(3,ngxat) :: r_ion_in_cell_gx
  real(kind=r8_kind) :: Q_imp,Q_shell_imp,PK_imp,Rimp(3)
  real(kind=r8_kind) :: Q_core_epe,Q_shell_epe,PK_epe,q_coupl_qmepe
  real(kind=r8_kind) :: bim,roim,cim,dim
  real(kind=r8_kind) :: b0,ro0,c0,d0
  integer(kind=i4_kind) :: common_atom !!!!!!!!!!!!!!AS
  real(kind=r8_kind) :: ann,r_sh(3)
  logical :: type_number(ndt,ndcell)
  real(kind=r8_kind) :: cutoff
  logical :: next
  logical :: lgxcell,lcellvec,lgulpv,lgulpp,lcell
  integer(kind=i4_kind) :: type_imp
  character(len=3) :: name_of_imp,imp_name(2),atm_name(2)
  integer(kind=i4_kind) :: counter
  character(len=80) :: comments
  integer(kind=i4_kind), dimension(ndt):: ion_types,position
  integer(kind=i4_kind) :: i,ii,jj,j,k,ix,li,lj,na,na1,l,m, &
       STEP_OF_DISPLACEMENT=200,mti,output_level=3
  integer(kind=i4_kind) :: status
  logical :: epein,options_read_parameters
  logical :: lcell_a=.false. ! translation cell in units of A
  real(kind=r8_kind), dimension(3):: xx, yy, zz
  real(kind=r8_kind):: q_ml,r_ml(3),r_first_sphere=5.0,r_2A_sphere=15.0 &
       ,r_first_short_interact=10.0,r_second_short_interact=10.0
  character(len=120) :: q_name; character(len=20) :: q_access, q_position
  integer(kind=i4_kind) :: n_3_body=0,i1,i2,i3
  integer(kind=i4_kind), allocatable :: n_3b(:),indexx(:,:,:)
  character(len=3) :: atm_nm(6,3)
  real(kind=r8_kind) :: k_i(6)
  real(kind=r8_kind) :: theta_i(6)
  real(kind=r8_kind) :: r_3b
  character(len=80) :: task 
  integer(kind=i4_kind) :: interfaced_mode=0_i4_kind
  character(len=80) :: message
  integer(i4_kind) :: istat

  namelist /tasks/ task
  namelist /rotation_shift/ xx, yy, zz, shft_gto_to_epe, shft_gto_to_hds
  namelist /n_gener_ions/ n_gen_ions
  namelist /n_ions_in_primitive_cell/ n_ions_cell
  namelist /n_translation_of_primitive_cell/ n_trans_primitive_cell
  namelist /parameters_of_generation/ n_centres_of_generation
  namelist /centres_of_generation/ xyz_center_gener,xyz_step_of_shift
  namelist /main_parameters/ radius_first_sphere,radius_2A_sphere, &
           n_types_central_atoms_3body,radius_long_interact,error_function_parameter, &
           option_c3_symm
  namelist /epe_relax_work_options/ use_pgdata, &
                          lcell_a,interfaced_mode, qm_interfaced_mode, &
                          basic_action, &
                          moving_epecenter, &
                          coordinates_moving_epecenter, &
                          output_level, &
                          read_configuration,write_configuration, &
                          print_ions_from_first_sphere, &
                          step_of_displacement, &
                          make_epe_reference, use_epe_reference, &
                          operations_make_reg_reference,&
			  ml_displacements, &
			  fixed_dielectric_const, fx_dielectric_const, &
                          ml_cluster_simulated, n_ml_cluster_simulators, &
                          ml_tensors, &
                          GX_HIGHPREC, extended_epe
  namelist /epe_paragauss/ make_epe_reference, use_epe_reference, &
                           use_epe_reference, &
                           use_pgdata, &
                           extended_epe, &
                           embed_convergence_check, &
                           embed_convergence_limit, &  ! energy gain in epe optimization
                           explicit_coupling
  namelist /epe_output_charges_and_types/ out_charges,an_epe_type
  namelist /epe_center/ Q_core_epe,Q_shell_epe,PK_epe,q_coupl_qmepe
  namelist /epe_param_of_potential/ atm_name,b0,ro0,c0,d0,next,cutoff,common_atom !!!!!!!!!!!!!!!!AS
  namelist /vacancy_impurity/ n_vacancies,n_impurities,n_types_central_atoms_3body_im
  namelist /impurities/ name_of_imp,type_imp,Q_imp,Q_shell_imp,PK_imp,rimp
  namelist /imp_param_of_potential/ imp_name,bim,roim,cim,dim,next,cutoff,common_atom !!!!!!!!!!!!!!!AS
  namelist /three_body_interaction/ atm_nm,k_i,theta_i,r_3b,n_3_body
  namelist /cluster_simulators/ q_ml,r_ml
  namelist /optimization_parameters/ n_iterations,abs_g,n_hess_update_cycles,weight_hess, &
       print_gradients,n_pgepe_iterations
  namelist /xyz_output/ unit_cell_output,region_I_output,output_step
  namelist /calc_pot/ point00,point01,point02,dxaxis,dyaxis,nxgrid,nygrid,output_units, &
       v_abs_limit,output_format
  !========================================================================
  call convert_epe_input()
  input_epe=openget_iounit(trim(epe_input_dir)//'/epe.inp.interm', &
				 form='formatted', status='old')
	
  output_epe=get_iounit()
  open(output_epe,file=trim(epe_input_dir)//'/epe.out', form='formatted',&
 				            status='unknown', position='append')

  pc_aswritemode=.true.

  scaling_factor=1.889726877774_r8_kind*auangs
  radius_long_interact= 100.0_r8_kind
  r_2A_sphere=18.0_r8_kind
  r_first_sphere=5.0_r8_kind
  ERROR_FUNCTION_PARAMETER=df_ERROR_FUNCTION_PARAMETER

  radius_first_sphere=r_first_sphere 
  radius_2A_sphere=r_2A_sphere

  ! work parameters !!!!!!!!!!!!!!!!!!!!!!!!
  read (input_epe,nml=tasks)

  if (index(task,'periodic_optimization') /= 0 .or. &
       index(task,'PERIODIC_OPTIMIZATION') /= 0) then
     periodic_optimization=.true.
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '       The Optimization of The Unit Cell         '
     write(output_epe,*) '-------------------------------------------------'
  else if (index(task,'epe_relaxation') /= 0 .or. &
       index(task,'EPE_RELAXATION') /= 0  .or. &
       index(task,'EPE_relaxation') /= 0) then
     periodic_optimization=.false.
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '       The Relaxation of EPE Enviroment          '
     write(output_epe,*) '-------------------------------------------------'
  else
     call error_handler("read_epe_input_new:NAMELIST TASK: wrong TASK parameter")
  endif

  lgulpp=.false.; lgulpv=.false.; lcell=.false.
  lgulpp=(index(task,'gulpp') /= 0 .or. &
       index(task,'GULPP') /= 0) 
  lgulpv=(index(task,'gulpv') /= 0 .or. &
       index(task,'GULPV') /= 0)
  lcell=(index(task,'gxcell') /= 0 .or. &
       index(task,'GXCELL') /= 0)
  if((lgulpp .and. lgulpv).or. & 
     (lgulpp .and. lcell).or.  & 
     (lgulpv .and. lcell)) call error_handler( &
     "either gulpp or gulpv or lcell can be present at the same time")

  lpotcalc=(index(task,'calc_pot') /= 0 .or. &
       index(task,'CALC_POT') /= 0)

  if(periodic_optimization) then
     basic_action = 1
     interfaced_mode = 0
  endif

  read (input_epe,nml=optimization_parameters)
  write(output_epe,*) ' '
  write(output_epe,nml=optimization_parameters)

  read (input_epe,nml=xyz_output)
  write(output_epe,*) ' '
  write(output_epe,nml=xyz_output)

  ! unit cell  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lgxcell=.false.; lcellvec=.false.
  if(lgulpp) then
     inquire (file=trim(epe_input_dir)//'/cellvec',exist=lcellvec)
     if(lcellvec) gxcell_unit=openget_iounit(trim(epe_input_dir)//'/cellvec',  &
          form='formatted', status='unknown')
     print*, 'gxcell_unit', gxcell_unit, trim(epe_input_dir)//'/cellvec',lcellvec
  else if(lgulpv .or. lcell) then
     inquire (file=trim(epe_input_dir)//'/gxcell',exist=lgxcell)
     if(lgxcell) gxcell_unit=openget_iounit(trim(epe_input_dir)//'/gxcell',  &
          form='formatted', status='unknown')
  endif

  read(input_epe,'(a80)')COMMENTS
  if (trim(adjustl(COMMENTS)) == "TRANSLATION VECTORS" .or. &
       trim(adjustl(COMMENTS)) == "translation vectors") then
     write(output_epe,*)COMMENTS
     DO j=1,3
        read(input_epe,*)(VECTORS_TRANS(I,J),i=1,3)
     enddo	! I=1,3
     if(lgulpp) then
        write(output_epe,*) '  vectors are taken from gx.cv file'
        DO j=1,3
           read(gxcell_unit,*)(VECTORS_TRANS(I,J),i=1,3)
        enddo	! I=1,3
     endif
     do j=1,3
        write(output_epe,'(3(3x,F15.9))')(VECTORS_TRANS(I,J),i=1,3)
     enddo
     VECTORS_TRANS=VECTORS_TRANS*scaling_factor
  else
     call error_handler("read_epe_input_new: Not TRANSLATION VECTORS")
  endif

  read(input_epe,nml=n_ions_in_primitive_cell)
  write(output_epe,*) ' '
  write(output_epe, nml=n_ions_in_primitive_cell)
  allocate(which_epe_ion(n_ions_cell),stat=epealloc_stat(7))
  ASSERT(epealloc_stat(7).eq.0)
  epealloc_stat(7)=1

101 format (18a4)

  read(input_epe,'(a80)')COMMENTS
  if (trim(adjustl(COMMENTS)) == "number, type, name and coordinates of ions in the primitive cell" .or. &
       trim(adjustl(COMMENTS)) == "NUMBER, TYPE, NAME AND COORDINATES OF IONS IN THE PRIMITIVE CELL") then
     write(output_epe, *)trim(adjustl(COMMENTS))
     DPRINT trim(adjustl(COMMENTS))
     na1=1
  elseif(trim(adjustl(COMMENTS)) == "type, name and coordinates of ions in the primitive cell" .or. &
       trim(adjustl(COMMENTS)) == "TYPE, NAME AND COORDINATES OF IONS IN THE PRIMITIVE CELL") then
     write(output_epe, *)trim(adjustl(COMMENTS))
     DPRINT trim(adjustl(COMMENTS))
     na1=0
  else
     call error_handler( &
          "read_epe_input_new: Not number, type, name and coordinates of ions in the primitive cell")
  endif

!!$  inquire (file=trim(epedata_dir)//'/gxcell',exist=lgxcell)
  if(lgxcell) then
     write(output_epe,*) 'new input: read coordinates of cell from gx.c file'
     print*,'new input: read coordinates of cell from gx.c file'
  endif
  if(lcellvec) then
     write(output_epe,*) 'read coordinates of cell from gx.cv file'
     print*,'read coordinates of cell from gx.cv file',N_IONS_CELL,na1
  endif

  fixed(1:N_IONS_CELL)=.false.

  max_type_ions=0
  do I=1,N_IONS_CELL
     read(input_epe,'(a80)')message
     if (na1==1) then 
        READ (message,*,end=32) na,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:),fixed(i)
        goto 30
32      READ (message,*) na,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
     end if
     if (na1==0) then 
        READ (message,*,end=33) type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:),fixed(i)
        goto 30
33      READ (message,*) type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
     end if
30   if(.not.lgxcell .and. .not.lcellvec) then
        if(fixed(i)) write(output_epe,1016) i,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:),fixed(i)
        if(.not.fixed(i)) write(output_epe,1006) i,type_of_ion(i),NAME_OF_IONS(i),R_ION_IN_CELL(I,:)
     endif
     ion_types(TYPE_OF_ION(I))=type_of_ion(i)
     name_of_type(TYPE_OF_ION(I))=NAME_OF_IONS(i)
     if(max_type_ions.lt.TYPE_OF_ION(I)) max_type_ions=TYPE_OF_ION(I)
  enddo! I=1,N_IONS_CELL
1006 FORMAT (1x,I3,6x,I2,4X,A3,5x,3F10.5)
1016 FORMAT (1x,I3,6x,I2,4X,A3,5x,3F10.5,1x,l1)

  if(lgxcell .or. lcellvec) then
     type_number=.false.
     R_ION_IN_CELL=zero; R_SHELL_IN_CELL=zero; R_CORE_IN_CELL=zero
     do I=1,N_IONS_CELL
        read(gxcell_unit,'(f5.2,3f15.7)') an(i),r_core_in_cell(i,:)
        DPRINT i, an(i),r_core_in_cell(i,:)
        TYPE_OF_ION(i)=int(an(i))
        type_number(TYPE_OF_ION(i),i)=.true.
     enddo
     if(.not.lcell) then
        l1: do
           read(gxcell_unit,*,end=22)ann,r_sh
           DPRINT  ann,r_sh
           if(ann == zero) exit l1
           i=ann
           l2: do j=1,ndcell
              if(type_number(i,j))then
                 R_SHELL_IN_CELL(j,:)=r_sh
                 type_number(i,j)=.false.
                 exit l2
              endif
           enddo l2
        enddo l1
        DPRINT 'l1 loop done'
22      backspace(gxcell_unit)
        DPRINT 'backspace(gxcell_unit) done'
     end if
     if(.not. periodic_optimization) call returnclose_iounit(gxcell_unit)
     DPRINT 'gxcell_unit returnclose_iounit done'
     do i=1,N_IONS_CELL
        if(fixed(i)) then
           write(output_epe,1017) i,type_of_ion(i),NAME_OF_IONS(i),' c ',&
                r_core_in_cell(i,:),fixed(i)
           if(dot_product(R_SHELL_IN_CELL(i,:),R_SHELL_IN_CELL(i,:)) /= zero) &
                write(output_epe,1017) i,type_of_ion(i),NAME_OF_IONS(i),' s ',&
                r_shell_in_cell(i,:),fixed(i)
        else
           write(output_epe,1007) i,type_of_ion(i),NAME_OF_IONS(i),' c ',&
                r_core_in_cell(i,:)
           if(dot_product(R_SHELL_IN_CELL(i,:),R_SHELL_IN_CELL(i,:)) /= zero) &
                write(output_epe,1007) i,type_of_ion(i),NAME_OF_IONS(i),' s ',&
                r_shell_in_cell(i,:)
        end if
     enddo
1007 FORMAT (1x,I3,6x,I2,4X,2A3,5x,3F10.5)
1017 FORMAT (1x,I3,6x,I2,4X,2A3,5x,3F10.5,1x,l1)
     core_shell=.false.
     do i=1,n_ions_cell
        if(dot_product(R_SHELL_IN_CELL(i,:),R_SHELL_IN_CELL(i,:)) /= zero) then 
           core_shell=.true.
           r_ion_in_cell(i,:)=(r_core_in_cell(i,:)+r_shell_in_cell(i,:))/2.0_r8_kind
        else
           r_ion_in_cell(i,:)=r_core_in_cell(i,:)
           r_shell_in_cell(i,:)=r_core_in_cell(i,:)
        endif
     enddo
  else

     if(periodic_optimization.and..false.) then
        core_shell=.false.
        if(.not. lgxcell) gxcell_unit=openget_iounit(trim(epe_input_dir)//'/gxcell', &
             form='formatted', status='unknown')
        indgx=0
        impu=0
        do i=1,n_ions_cell
           indgx(1,i)=i
        enddo	! i=1,n_ions_cell
        do i=1,n_ions_cell
           an(i)=TYPE_OF_ION(i)
           write(gxcell_unit,203) an(i),(r_ion_in_cell(i,j),j=1,3)  &
                ,(indgx(k,i),k=1,8),impu(i)
        enddo	! i=1,n_ions_cell
        write(gxcell_unit,203) zero,zero,zero,zero,0,0,0,0,0,0,0,0,0
203     FORMAT((f5.2,3(2X,f13.7),2i4,2x,3I3,2X,3I3,2x,2i3))
     endif
  endif

! epe potential parameters !!!!!!!!!!!!!!!!!!!
  print*,'read nml=main_parameters'
  read (input_epe,nml=main_parameters)
  write (output_epe,nml=main_parameters)

  RADIUS_FIRST_SPHERE=RADIUS_FIRST_SPHERE*scaling_factor
  RADIUS_2A_SPHERE=RADIUS_2A_SPHERE*scaling_factor
  RADIUS_LONG_INTERACT=(RADIUS_LONG_INTERACT*scaling_factor)**2

  ALLOCATE(host_tmp%sr1(30,30,31), &
       host_tmp%k(30,30,31),&
       host_tmp%k1(30,30,31),&
       host_tmp%r1(30,30,31),&
       host_tmp%r0(30,30,31),&
       host_tmp%b(30,30,31),&
       host_tmp%ro(30,30,31),&
       host_tmp%c(30,30,31),&
       host_tmp%d(30,30,31),&
       host_tmp%sr2(30,30,31),&
       stat=epealloc_stat(3))
  if(epealloc_stat(3).ne.0) call error_handler("allocate host_tmp%sr1 failed")
  epealloc_stat(3)=1
  
  host_tmp%k= 0.0_r8_kind 
  host_tmp%k1= 0.0_r8_kind 
  host_tmp%r1= 0.0_r8_kind 
  host_tmp%r0= 0.0_r8_kind 

  host_tmp%b= 0.0_r8_kind 
  host_tmp%ro= 0.0_r8_kind 
  host_tmp%ro(:,:,1) = 1.0_r8_kind !!!!!!!!!!!!AS
  host_tmp%c= 0.0_r8_kind 
  host_tmp%d=0.0_r8_kind
  host_tmp%sr1= 4.5_r8_kind
  host_tmp%sr2= 4.5_r8_kind
  q_coupl_qmepe = 0.0_r8_kind

  do i=1,max_type_ions
     read (input_epe,nml=epe_center)
     pk(i)=pk_epe
     q_shell(i)=q_shell_epe*qau_qepe
     q_nuclear(i)=q_core_epe*qau_qepe
     q_ion(i)=q_shell(i)+q_nuclear(i)
     q_epecl_coupl(i)=q_coupl_qmepe
  enddo

  m=max_type_ions*(max_type_ions-1)/2+max_type_ions
  i=0
  do 
     next=.true.
     b0=0.0_r8_kind
     ro0=1.0_r8_kind
     c0=0.0_r8_kind
     d0=0.0_r8_kind
     cutoff=10.0_r8_kind
     common_atom=0  !!!!!!!!!!!!!!!!!AS
     print*,'read nml=epe_param_of_potential'
     read (input_epe,nml=epe_param_of_potential)
     i=i+1
     do j=1,max_type_ions
        if(trim(atm_name(1)) == trim(name_of_type(j))) exit
     enddo
     do k=1,max_type_ions
        if(trim(atm_name(2)) == trim(name_of_type(k))) exit
     enddo
     if(j > max_type_ions .or. k > max_type_ions) then
        write(message,'(i2)') i
        call error_handler("Your "//trim(message)//"-th namelist EPE_PARAM_OF_POTENTIAL"// &
             " contains wrong atom name")
     end if
     jj=common_atom+1                !!!!!!!!!!!!!!!AS

     host_tmp%k(J,K,jj)=0.0
     host_tmp%k(K,J,jj)=0.0
     host_tmp%k1(J,K,jj)=0.0
     host_tmp%k1(K,J,jj)=0.0
     host_tmp%r1(J,K,jj)=0.0
     host_tmp%r1(K,J,jj)=0.0
     host_tmp%r0(J,K,jj)=0.0
     host_tmp%r0(K,J,jj)=0.0

     host_tmp%B(J,K,jj)=b0           !!!!!!!!!!!!!!!AS
     host_tmp%B(K,J,jj)=b0           !!!!!!!!!!!!!!!AS
     host_tmp%Ro(J,K,jj)=ro0         !!!!!!!!!!!!!!!AS
     host_tmp%RO(K,J,jj)=ro0         !!!!!!!!!!!!!!!AS
     host_tmp%C(J,K,jj)=c0           !!!!!!!!!!!!!!!AS
     host_tmp%C(K,J,jj)=c0           !!!!!!!!!!!!!!!AS
     host_tmp%D(J,K,jj)=d0           !!!!!!!!!!!!!!!AS
     host_tmp%D(K,J,jj)=d0           !!!!!!!!!!!!!!!AS
     host_tmp%sr1(K,J,jj)=cutoff     !!!!!!!!!!!!!!!AS
     host_tmp%sr1(j,k,jj)=cutoff     !!!!!!!!!!!!!!!AS
     host_tmp%sr2(K,J,jj)=cutoff     !!!!!!!!!!!!!!!AS
     host_tmp%sr2(j,k,jj)=cutoff     !!!!!!!!!!!!!!!AS
     if(.not. next) exit
  enddo

  if (n_types_central_atoms_3body /= 0) then
     allocate(ki(max_type_ions,max_type_ions,max_type_ions), &
          theta_0(max_type_ions,max_type_ions,max_type_ions),stat=status)
     if(status.ne.0) call error_handler("allocate ki failed")
     ki=0.0_r8_kind
     allocate(types(n_types_central_atoms_3body,5), r3b(n_types_central_atoms_3body), &
          n_3b(n_types_central_atoms_3body),indexx(n_types_central_atoms_3body,6,3), &
          stat=status)
     if(status.ne.0) call error_handler("allocate types and r3b failed")

     types=0
     do i=1,n_types_central_atoms_3body
        read (input_epe,nml=three_body_interaction)
        n_3b(i)=n_3_body
        do j=1,n_3_body
           do k=1,max_type_ions
              if(atm_nm(j,1) == name_of_type(k)) i1=k
              if(atm_nm(j,2) == name_of_type(k)) i2=k
              if(atm_nm(j,3) == name_of_type(k)) i3=k
           enddo
           if(i1 > max_type_ions .or. i2 > max_type_ions .or. i3 > max_type_ions) then
              write(message,'(i2)') i
              call error_handler("Your "//trim(message)//"-th namelist THREE_BODY_INTRACTION"// &
                   " contains wrong atom name")
           end if
           ki(i1,i2,i3)=k_i(j)
           ki(i3,i2,i1)=k_i(j)
           theta_0(i1,i2,i3)=theta_i(j)
           theta_0(i3,i2,i1)=theta_i(j)
           indexx(i,j,1)=i1
           indexx(i,j,2)=i2
           indexx(i,j,3)=i3
           if(j==1) then
              types(i,1)=i2
           endif
           do l=2,5
              if(types(i,l)==i1) exit
              if(types(i,l)==0) then
                 types(i,l)=i1
                 exit
              endif
           enddo
           do l=2,5
              if(types(i,l)==i3) exit
              if(types(i,l)==0) then
                 types(i,l)=i3
                 exit
              endif
           enddo
        enddo
        r3b(i)=r_3b
     enddo
  endif

  write(output_epe,*) '==============================================================================='
  write(output_epe,*) '              The Potential Parameters of the Lattice'
  write(output_epe,*) '-------------------------------------------------------------------------------' 
  write(output_epe,*) '                        Two-body parameters'
  write(output_epe,*) '-------------------------------------------------------------------------------' 
  write(output_epe,*) '              b           ro          c           d        cutoff   common atom' 
  do i=1,max_type_ions
     do j=i,max_type_ions
        k=0
        do                        !<<<<<<<<<<<<<<<<AS
           k=k+1
           if(k > 30 ) exit
           if(k==1) then
              write(output_epe,'(2a4,5f12.5)') name_of_type(i),name_of_type(j) &
                   ,host_tmp%b(i,j,k),host_tmp%ro(i,j,k),host_tmp%c(i,j,k),host_tmp%d(i,j,k), &
                   host_tmp%sr1(i,j,k)
           else                                  
              if(host_tmp%b(i,j,k) == 0.0_r8_kind) cycle
              write(output_epe,'(2a4,5f12.5,4x,i4)') name_of_type(i),name_of_type(j) &
                   ,host_tmp%b(i,j,k),host_tmp%ro(i,j,k),host_tmp%c(i,j,k),host_tmp%d(i,j,k), &
                   host_tmp%sr1(i,j,k),k-1
           end if
        enddo                     !>>>>>>>>>>>>>>>>AS
     end do
  enddo
  
  write(output_epe,*) '-------------------------------------------'
  write(output_epe,*) '          Shell-model parameters        '
  write(output_epe,*) '-------------------------------------------'
  write(output_epe,*) '      q_nuclear     q_shell       pki   '
  do i=1,max_type_ions
     write(output_epe,'(a4,3f12.5)') name_of_type(i),&
          q_nuclear(i)/qau_qepe,q_shell(i)/qau_qepe,pk(i)
  enddo


  if (n_types_central_atoms_3body /= 0) then
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '              Three-body parameters'
     write(output_epe,*) '-------------------------------------------------'
     write(output_epe,*) '                  ki        theta       r3b'
     do i=1,n_types_central_atoms_3body
        do j=1,n_3b(i)
           write(output_epe,'(3a4,3f12.5)') name_of_type(indexx(i,j,1)), &
                name_of_type(indexx(i,j,2)),name_of_type(indexx(i,j,3)), &
                ki(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), &
                theta_0(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), r3b(i)
        enddo
     enddo
     deallocate(n_3b,indexx)
  endif
  write(output_epe,*) '====================================================================='

  do i=1,n_ions_cell
     do j=1,max_type_ions
        if(name_of_ions(i)==name_of_type(j)) then
           q_z(i)=q_ion(j)/qau_qepe
           exit
        endif
     enddo
  enddo

  if(lpotcalc) then
     point00=zero
     point01=zero
     point02=zero
     dyaxis=zero
     dxaxis=zero
     nxgrid=50
     nygrid=50
     v_abs_limit=zero
     output_units="eV-a"
     output_format="gnuplot"
     read (input_epe,nml=calc_pot)
     write (output_epe,nml=calc_pot)
     if(dot_product(point00-point01,point00-point01) == zero .or. &
          dot_product(point02-point00,point02-point00) == zero .or. &
          dot_product(point01-point02,point01-point02) == zero) call error_handler( &
          "Potential input: points defining the grid plane must be different")
     if(dyaxis(1) == dyaxis(2) .or. &
          dxaxis(1) == dxaxis(2)) call error_handler("Potential input: check dxaxis or dyaxis ")
     if(nxgrid <= 1 .or. nygrid <= 1) call error_handler("Potential input: grid is equal ZERO")
     if(v_abs_limit < zero) call error_handler("Potential_input: V_abs_limit < 0")
     if(index(output_units,"a.u.")==0 .and. index(output_units,"eV-a")==0) &
          call error_handler("Potential_input: Only a.u. or eV-a units can be applied")
     if(index(output_format,"scilab")==0 .and. index(output_format,"SCILAB")==0 .and. &
          index(output_format,"worksheet")==0 .and. index(output_format,"WORKSHEET")==0 .and. &
          index(output_format,"gnuplot")==0 .and. index(output_format,"GNUPLOT")==0) &
          call error_handler( &
          "Potential_input: Permited keywords are scilab(SCILAB), worksheet(WORKSHEET), gnulpot(GNUPLOT)")
     n_iterations=0
  endif

  ! epe options !!!!!!!!!!!!!!!!!!!!!
  if (periodic_optimization) then
     ALLOCATE(host%sr1(max_type_ions,max_type_ions,max_type_ions+1), &        !!!!!!!!!!!!!!!!AS
          host%b(max_type_ions,max_type_ions,max_type_ions+1),&
          host%ro(max_type_ions,max_type_ions,max_type_ions+1),&
          host%c(max_type_ions,max_type_ions,max_type_ions+1),&
          host%d(max_type_ions,max_type_ions,max_type_ions+1),&
          host%sr2(max_type_ions,max_type_ions,max_type_ions+1),&
          host%k(max_type_ions,max_type_ions,max_type_ions+1),&
          host%r0(max_type_ions,max_type_ions,max_type_ions+1),&
          host%k1(max_type_ions,max_type_ions,max_type_ions+1),&
          host%r1(max_type_ions,max_type_ions,max_type_ions+1), & 
          stat=epealloc_stat(1))
     if(epealloc_stat(1).ne.0) &
        call error_handler("deallocate host%sr1 failed(b)")
     epealloc_stat(1)=1 !sr1 sr2
     epealloc_stat(2)=1 ! b ro c d

     host%k=0.0_r8_kind
     host%k1=0.0_r8_kind
     host%r1=0.0_r8_kind
     host%r0=0.0_r8_kind

     host%b=0.0_r8_kind; host%ro=0.0_r8_kind; host%ro(:,:,1)=1.0_r8_kind        !!!!!!!!!!!!!!!!!!AS
     host%c=0.0_r8_kind; host%d=0.0_r8_kind
     host%sr1=10.0_r8_kind
     host%sr2=10.0_r8_kind

     host%k=host_tmp%k(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)      
     host%k1=host_tmp%k1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)   
     host%r1=host_tmp%r1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)  
     host%r0=host_tmp%r0(1:max_type_ions,1:max_type_ions,1:max_type_ions+1) 

     host%sr1=host_tmp%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)       !!!!!!!!!!!!!!!!AS
     host%b=host_tmp%b(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)           !!!!!!!!!!!!!!!!AS
     host%ro=host_tmp%ro(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)         !!!!!!!!!!!!!!!!AS
     host%c=host_tmp%c(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)           !!!!!!!!!!!!!!!!AS
     host%d=host_tmp%d(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)           !!!!!!!!!!!!!!!!AS
     host%sr2=host_tmp%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)       !!!!!!!!!!!!!!!!AS

     n_trans_primitive_cell=6
     n_gen_ions=n_ions_cell*30
     if(n_gen_ions < 1500) n_gen_ions=1500
     if(n_gen_ions > 6000) n_gen_ions=6000
     if(n_ions_cell <= 10) then
        n_centres_of_generation=n_ions_cell
        DO I=1,N_CENTRES_OF_GENERATION
           R_CENT_GENER(I,:)=r_ion_in_cell(i,:)*scaling_factor
           shift_cent_gener(i,:)=zero
        enddo
     else
        n_centres_of_generation=10
        DO I=1,N_CENTRES_OF_GENERATION
           j=int((n_ions_cell*1.0_r8_kind/10.0_r8_kind)*i)
           R_CENT_GENER(I,:)=r_ion_in_cell(j,:)*scaling_factor
           shift_cent_gener(i,:)=zero
        enddo
     endif
     allocate(r_nuc_ion(n_gen_ions,3), &
          r_sh_ion(n_gen_ions,3), &
          epe(n_gen_ions),q_zl(n_gen_ions),stat=epealloc_stat(9))
      ASSERT(epealloc_stat(9).eq.0)
             epealloc_stat(9)=1 ! r_nuc_ion r_sh_ion epe 
             epealloc_stat(10)=1 ! q_zl
  else
     print*,' read epe_relax_work_options'
     read (input_epe,nml=epe_relax_work_options)
     if(interfaced_mode.ne.0) qm_interfaced_mode=.true.
     reg_2a_treated=ml_displacements.and..not.ml_cluster_simulated

     if(operations_make_reg_reference) then
        make_epe_reference=.false.
        use_epe_reference=.false.   
        use_pgdata = .false.
     end if
     write(output_epe,nml=epe_relax_work_options,IOSTAT=istat)
     if( istat /= 0 )then
       ABORT('error writing nml=epe_relax_work_options')
     endif

     DPRINT 'epe_output_charges_and_types use_epe_reference',use_epe_reference
     if(operations_make_reg_reference .or. use_pgdata.or.use_epe_reference) then
        do i=1,max_type_ions
           out_charges(i)=q_ion(i)/qau_qepe
        enddo
        print*,' read epe_output_charges_and_types'
        read (input_epe,nml=epe_output_charges_and_types)
        write(output_epe,*) '================================================'
        write(output_epe,*) 'Charges that will be used in QM-EPE calculations'
        write(output_epe,*) '------------------------------------------------'
        write(output_epe,*) '               Q_out       Qs_out      Qc_out'
        write(output_epe,*) '------------------------------------------------'
        do i=1,max_type_ions
           write(output_epe,'(4x,a4,4x,6f12.5)') name_of_type(i),out_charges(i), &
                out_charges(i)-q_nuclear(i)/qau_qepe,q_nuclear(i)/qau_qepe
        enddo
        write(output_epe,*) '------------------------------------------'
        write(output_epe,*) 'Types of EPE atoms in  QM-EPE calculations'
        write(output_epe,*) '------------------------------------------'
        do i=1,max_type_ions
           write(output_epe,'(10x,a4,10x,f5.2)') name_of_type(i),an_epe_type(i)
        enddo
        write(output_epe,*) '=========================================='
        scale_factor=out_charges(1)/(q_ion(1)/qau_qepe)
     endif

     if(ml_cluster_simulated.and.n_ml_cluster_simulators.ne.0) then
        allocate (ml_cluster_simulators(n_ml_cluster_simulators),stat=status)
        if(status.ne.0) call error_handler("allocate ml_cluster_simulators failed")
        do counter=1,n_ml_cluster_simulators
           read (input_epe,nml=cluster_simulators)
           ml_cluster_simulators(counter)%q=q_ml*qau_qepe
           ml_cluster_simulators(counter)%r=r_ml*auangs
        end do
     end if
     allocate(ml_dprime_factors(max_type_ions),ml_fac(max_type_ions),stat=epealloc_stat(8))
     ASSERT(epealloc_stat(8).eq.0)
            epealloc_stat(8)=1

     explicit_coupling=.false.
     ! note that use_epe_reference may be redefined here if use_pgdata eq true
     ! but use_epe_reference can not be else
     if(use_pgdata.or.use_epe_reference) then
       print*,'read nml=epe_paragauss'
        read (input_epe,nml=epe_paragauss)
        write(output_epe,nml=epe_paragauss)
        
        if(explicit_coupling) then
            call treat_epepar_namelist()

           ALLOCATE(host%sr1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1), & !!!!!!!!!!!!!!!!AS
                host%b(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%ro(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%c(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%d(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%k(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%r0(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%k1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%r1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                host%sr2(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& !!!!!!!!!!!!!!!!AS
                stat=epealloc_stat(1))
           if(epealloc_stat(1).ne.0) &
              call error_handler("1b allocate host failed(a)")
           epealloc_stat(1)=1 !sr1 sr2
           epealloc_stat(2)=1 ! b ro c d

           host%k=0.0_r8_kind 
           host%k1=0.0_r8_kind 
           host%r1=0.0_r8_kind 
           host%r0=0.0_r8_kind 

           host%b=0.0_r8_kind 
           host%ro=0.0_r8_kind; host%ro(:,:,1)=1.0_r8_kind !!!!!!!!!!!!!AS
           host%c=0.0_r8_kind; host%d=0.0_r8_kind
           host%sr1=10.0_r8_kind
           host%sr2=10.0_r8_kind
           host%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &               !<<<<<<<<<<<<<<AS
                host_tmp%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)

           host%k(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%k(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%k1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%k1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%r1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%r1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%r0(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%r0(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)

           host%b(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%b(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%ro(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%ro(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%c(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%c(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%d(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%d(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)         !>>>>>>>>>>>AS

           do i=1,ec_max_type_ions
              do j=i,ec_max_type_ions
                 if(i<=max_type_ions .and. j<=max_type_ions) cycle
                 mn=min(ec_name_of_type(i),ec_name_of_type(j))
                 mx=max(ec_name_of_type(i),ec_name_of_type(j))
                 call get_slsp(mn,mx,0,nr)
                 if(nr == 0) cycle
                 host%k(i,j,:)=par%item%k
                 host%k(j,i,:)=par%item%k
                 host%k1(i,j,:)=par%item%k1
                 host%k1(j,i,:)=par%item%k1
                 host%r1(i,j,:)=par%item%r1
                 host%r1(j,i,:)=par%item%r1
                 host%r0(i,j,:)=par%item%r0
                 host%r0(j,i,:)=par%item%r0

                 host%b(i,j,:)=par%item%b
                 host%b(j,i,:)=par%item%b
                 host%ro(i,j,:)=par%item%r
                 host%ro(j,i,:)=par%item%r
                 host%c(j,i,:)=par%item%c
                 host%c(i,j,:)=par%item%c
                 host%d(j,i,:)=par%item%d
                 host%d(i,j,:)=par%item%d
                 host%sr1(j,i,:)=par%item%cutoff
                 host%sr1(i,j,:)=par%item%cutoff
                 host%sr2(j,i,:)=par%item%cutoff
                 host%sr2(i,j,:)=par%item%cutoff
              enddo
           enddo

           write(output_epe,*) '========================================================================='
           write(output_epe,*) '              The Potential Parameters of the Explicit_coupling'
           write(output_epe,*) '-------------------------------------------------------------------------' 
           write(output_epe,*) '                          Two-body parameters'
           write(output_epe,*) '-------------------------------------------------------------------------' 
           write(output_epe,*) '                  b           ro          c           d         cutoff' 
           do i=1,ec_max_type_ions
              do j=i,ec_max_type_ions
                 if(i<=max_type_ions .and. j<=max_type_ions) cycle
                 write(output_epe,'(2f6.2,5f12.5)') ec_name_of_type(i),ec_name_of_type(j) &
                      ,host%b(i,j,1),host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1),host%sr1(i,j,1)
                 write(output_epe,'(12x,5f12.5)') &
                      host%k(i,j,1),host%r0(i,j,1),host%k1(i,j,1),host%r1(i,j,1)
              enddo
           enddo
           write(output_epe,*) '========================================================================='
        endif
     endif

     print*,'read nml=n_gener_ions'
     read(input_epe,nml=n_gener_ions)
     write(output_epe, nml=n_gener_ions)
     allocate(r_nuc_ion(n_gen_ions,3), &
          r_sh_ion(n_gen_ions,3), &
          epe(n_gen_ions),q_zl(n_gen_ions),stat=epealloc_stat(9))
          ASSERT(epealloc_stat(9).eq.0)
                 epealloc_stat(9)=1
                 epealloc_stat(10)=1 ! q_zl

     print*,'read nml=n_translation_of_primitive_cell'
     read (input_epe,nml=n_translation_of_primitive_cell)
     write (output_epe,nml=n_translation_of_primitive_cell)

     print*,'read nml=parameters_of_generation'
     read (input_epe,nml=parameters_of_generation)
     write (output_epe,nml=parameters_of_generation)

     if(n_centres_of_generation.ne.0) then
        DO I=1,N_CENTRES_OF_GENERATION
           read (input_epe,nml=centres_of_generation)
           write (output_epe,nml=centres_of_generation)
           R_CENT_GENER(I,:)=xyz_center_gener(:)*scaling_factor
           shift_cent_gener(i,:)=xyz_step_of_shift(:)*scaling_factor
        enddo! I=1,N_CENTRES_OF_GENERATION
     endif

     print*,'read nml=rotation_shift'
     read (input_epe,nml=rotation_shift)
     write (output_epe,nml=rotation_shift)

     rot_gto_to_epe(:,1)=xx
     rot_gto_to_epe(:,2)=yy
     rot_gto_to_epe(:,3)=zz

     ! parameters for impurities !!!!!!!!!!!!!!!!  

     print*,'read nml=vacancy_impurity'
     read (input_epe,nml=vacancy_impurity)
     write (output_epe,nml=vacancy_impurity)

     IF(N_IMPURITIES.EQ.0) then
        write(output_epe, 63)
63      format(' there are no basic impurities ')
        if(.not.explicit_coupling) then 
           ALLOCATE(host%sr1(max_type_ions,max_type_ions,max_type_ions+1), &      !<<<AS
                host%b(max_type_ions,max_type_ions,max_type_ions+1),&
                host%k(max_type_ions,max_type_ions,max_type_ions+1),&
                host%k1(max_type_ions,max_type_ions,max_type_ions+1),&
                host%r1(max_type_ions,max_type_ions,max_type_ions+1),&
                host%r0(max_type_ions,max_type_ions,max_type_ions+1),&
                host%ro(max_type_ions,max_type_ions,max_type_ions+1),&
                host%c(max_type_ions,max_type_ions,max_type_ions+1),&
                host%d(max_type_ions,max_type_ions,max_type_ions+1),&
                host%sr2(max_type_ions,max_type_ions,max_type_ions+1),&          !>>>>>>>>AS
                stat=status)
           if(status.ne.0) then
              call error_handler("deallocate host%sr1 failed(c)")
           end if
           host%k=0.0_r8_kind 
           host%k1=0.0_r8_kind 
           host%r1=0.0_r8_kind 
           host%r0=0.0_r8_kind 

           host%b=0.0_r8_kind 
           host%ro=0.0_r8_kind 
           host%ro(:,:,1)=1.0_r8_kind     !!!!!!!!!!!!!AS
           host%c=0.0_r8_kind; host%d=0.0_r8_kind
           host%sr1=10.0_r8_kind
           host%sr2=10.0_r8_kind

           host%k(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
                host_tmp%k(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%k1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
                host_tmp%k1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%r0(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
                host_tmp%r0(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%r1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
                host_tmp%r1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)

           host%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &              !<<<<<<<<<<<AS
                host_tmp%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%b(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
                host_tmp%b(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%ro(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%ro(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%c(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%c(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%d(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%d(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
           host%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &
                host_tmp%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)         !>>>>>>>>>>>>>AS
        endif
     else
        read(input_epe, 101)COMMENTS
        read(input_epe, 101)COMMENTS
        mti=max_type_ions
        DO I=1,N_IMPURITIES
           read (input_epe,nml=impurities)

           TYPE_IMPURITY(I)=type_imp
           if(TYPE_IMPURITY(I).le.mti) then
              write(output_epe,*) 'read_epe_input:type of impurity cannot be smaller then'
              write(output_epe,*) 'largest type of ions'
              call error_handler("read_epe_input:input error ")
           endif
           if(max_type_ions < TYPE_IMPURITY(I)) then
              max_type_ions=TYPE_IMPURITY(I)
              name_of_type(max_type_ions)=name_of_imp
              do j=1,max_type_ions-1
                 if(name_of_type(max_type_ions) == name_of_type(j)) &
                      call error_handler("read_epe_input: Different types of ions must have different names")
              enddo
           endif
           Q_IMPURITY(I)=Qau_qepe*Q_imp
           Q_SH_IMPURITY(I)=Qau_qepe*q_shell_imp
           Q_NUC_IMPURITY(I)=Q_IMPURITY(I)-Q_SH_IMPURITY(I)
           pk_impurity(i)=PK_imp
           R_imp(I,:)=rimp(:)*scaling_factor
           R_SH_IMP(I,:)=R_imp(I,:)
           R_SH_IMPO(I,:)=R_imp(I,:)
           R_NUC_IMP(I,:)=R_imp(I,:)
           R_NUC_IMPO(I,:)=R_imp(I,:)
        enddo        !  I=1,N_IMPURITIES

        if(.not.explicit_coupling) then 
           ALLOCATE(host%sr1(max_type_ions,max_type_ions,max_type_ions+1), &      !<<<<<<<<<<<AS
                host%k(max_type_ions,max_type_ions,max_type_ions+1),&
                host%k1(max_type_ions,max_type_ions,max_type_ions+1),&
                host%r1(max_type_ions,max_type_ions,max_type_ions+1),&
                host%r0(max_type_ions,max_type_ions,max_type_ions+1),&
                host%b(max_type_ions,max_type_ions,max_type_ions+1),&
                host%ro(max_type_ions,max_type_ions,max_type_ions+1),&
                host%c(max_type_ions,max_type_ions,max_type_ions+1),&
                host%d(max_type_ions,max_type_ions,max_type_ions+1),&
                host%sr2(max_type_ions,max_type_ions,max_type_ions+1),&           !>>>>>>>>>>>>AS
                stat=status)
           if(status.ne.0) then
              call error_handler("deallocate host%sr1 failed(d)")
           end if
           host%k=0.0_r8_kind 
           host%k1=0.0_r8_kind 
           host%r1=0.0_r8_kind 
           host%r0=0.0_r8_kind 

           host%b=0.0_r8_kind 
           host%ro=0.0_r8_kind; host%ro(:,:,1)=1.0_r8_kind
           host%c=0.0_r8_kind; host%d=0.0_r8_kind
           host%sr1=10.0_r8_kind
           host%sr2=10.0_r8_kind
           host%sr1(1:mti,1:mti,1:mti+1)=host_tmp%sr1(1:mti,1:mti,1:mti+1)   !<<<<<<<<<<<AS

           host%k(1:mti,1:mti,1:mti+1)=host_tmp%k(1:mti,1:mti,1:mti+1)
           host%k1(1:mti,1:mti,1:mti+1)=host_tmp%k1(1:mti,1:mti,1:mti+1)
           host%r1(1:mti,1:mti,1:mti+1)=host_tmp%r1(1:mti,1:mti,1:mti+1)
           host%r0(1:mti,1:mti,1:mti+1)=host_tmp%r0(1:mti,1:mti,1:mti+1)

           host%b(1:mti,1:mti,1:mti+1)=host_tmp%b(1:mti,1:mti,1:mti+1)
           host%ro(1:mti,1:mti,1:mti+1)=host_tmp%ro(1:mti,1:mti,1:mti+1)
           host%c(1:mti,1:mti,1:mti+1)=host_tmp%c(1:mti,1:mti,1:mti+1)
           host%d(1:mti,1:mti,1:mti+1)=host_tmp%d(1:mti,1:mti,1:mti+1)
           host%sr2(1:mti,1:mti,1:mti+1)=host_tmp%sr2(1:mti,1:mti,1:mti+1)    !>>>>>>>>>AS
        endif

        m=max_type_ions*(max_type_ions-1)/2+max_type_ions
        m=m-(mti*(mti-1)/2+mti)
        do i=1,m
           next=.true.
           bim=0.0_r8_kind
           roim=1.0_r8_kind
           cim=0.0_r8_kind
           dim=0.0_r8_kind
           cutoff=3.5_r8_kind
           common_atom=0         !!!!!!!!!!!!!AS
           read (input_epe,nml=imp_param_of_potential)
           do j=1,max_type_ions
              if(imp_name(1) == name_of_type(j)) exit
           enddo
           do k=1,max_type_ions
              if(imp_name(2) == name_of_type(k)) exit
           enddo
           jj=common_atom+1      !!!!!!!!!!!!!!AS

           host%k(J,K,jj)=0.0
           host%k(K,J,jj)=0.0
           host%k1(J,K,jj)=0.0
           host%k1(K,J,jj)=0.0
           host%r1(J,K,jj)=0.0
           host%r1(K,J,jj)=0.0
           host%r0(J,K,jj)=0.0
           host%r0(K,J,jj)=0.0

           host%B(J,K,jj)=bim      !<<<<<<<<<<<AS
           host%B(K,J,jj)=bim
           host%Ro(J,K,jj)=roim
           host%RO(K,J,jj)=roim
           host%C(J,K,jj)=cim
           host%C(K,J,jj)=cim
           host%D(J,K,jj)=dim
           host%D(K,J,jj)=dim
           host%sr1(K,J,jj)=cutoff
           host%sr1(j,k,jj)=cutoff
           host%sr2(K,J,jj)=cutoff
           host%sr2(j,k,jj)=cutoff  !>>>>>>>>>>>AS
           if(.not. next) exit
        enddo

        write(output_epe,*) '==================================================================='
        write(output_epe,*) '    The Potential Parameters of the Impurities'
        write(output_epe,*) '-------------------------------------------------------------------'
        write(output_epe,*) '                 Two-body parameters'
        write(output_epe,*) '-------------------------------------------------------------------'
        write(output_epe,*) '            b_im        ro_im       c_im        d_im    common_atom' 
        
        do i=1,max_type_ions
           k=mti+1
           if(k < i) k=i
           do j=k,max_type_ions
              do jj=1,max_type_ions             !!!!!!!!!!!!!!!AS
                 if(jj==1) then
                    write(output_epe,'(2a4,4f12.5)')name_of_type(i),name_of_type(j),&
                         host%b(i,j,1),host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1)
                    write(output_epe,'(8x,4f12.5)')&
                         host%k(i,j,1),host%r0(i,j,1),host%k1(i,j,1),host%r1(i,j,1)
                 else
                    if(host%b(i,j,jj)==0.0_r8_kind) cycle
                    write(output_epe,'(2a4,4f12.5,4x,i4)')name_of_type(i),name_of_type(j),&
                         host%b(i,j,jj),host%ro(i,j,jj),host%c(i,j,jj),host%d(i,j,jj),jj-1
                    write(output_epe,'(8x,4f12.5,4x,i4)')&
                         host%k(i,j,jj),host%r0(i,j,jj),host%k1(i,j,jj),host%r1(i,j,jj)
                 end if
              end do
           enddo
        enddo
        write(output_epe,*) '-------------------------------------------------------------------'

        write(output_epe,*) '        Shell-model parameters 1'
        write(output_epe,*) '-------------------------------------------'
        write(output_epe,*) '      q_nuclear     q_shell       pk'
        do i=mti,max_type_ions
           do j=1,n_impurities
              if(i == type_impurity(j)) then
                 write(output_epe,'(a4,4f12.5)') name_of_type(i),&
                      q_nuc_impurity(j)/qau_qepe,q_sh_impurity(j)/qau_qepe,pk_impurity(j)
                 exit
              endif
           enddo
        enddo
     
        if(n_types_central_atoms_3body_im > 0) then
           atm_nm='  '
           allocate(types_im(n_types_central_atoms_3body_im,5),r3b_im(n_types_central_atoms_3body_im), &
                n_3b(n_types_central_atoms_3body_im),indexx(n_types_central_atoms_3body_im,6,3), &
                stat=status)
           if(status.ne.0) call error_handler("allocate types and r3b failed")
           types_im=0
           do i=1,n_types_central_atoms_3body_im
              print*,'read nml=three_body_interaction'
              read (input_epe,nml=three_body_interaction)
              n_3b(i)=n_3_body
              do j=1,n_3_body
                 do k=1,max_type_ions
                    if(atm_nm(j,1) == name_of_type(k)) i1=k
                    if(atm_nm(j,2) == name_of_type(k)) i2=k
                    if(atm_nm(j,3) == name_of_type(k)) i3=k
                 enddo
                 ki(i1,i2,i3)=k_i(j)
                 ki(i3,i2,i1)=k_i(j)
                 theta_0(i1,i2,i3)=theta_i(j)
                 theta_0(i3,i2,i1)=theta_i(j)
                 indexx(i,j,1)=i1
                 indexx(i,j,2)=i2
                 indexx(i,j,3)=i3
                 if(j==1) then
                    types_im(i,1)=i2
                 endif
                 do l=2,5
                    if(types_im(i,l)==i1) exit
                    if(types_im(i,l)==0) then
                       types_im(i,l)=i1
                       exit
                    endif
                 enddo
                 do l=2,5
                    if(types_im(i,l)==i3) exit
                    if(types_im(i,l)==0) then
                       types_im(i,l)=i3
                       exit
                    endif
                 enddo
              enddo
              r3b_im(i)=r_3b
           enddo

           write(output_epe,*) '-------------------------------------------------'
           write(output_epe,*) '              Three-body parameters'
           write(output_epe,*) '-------------------------------------------------'
           write(output_epe,*) '                  ki        alpha       r3b'
           do i=1,n_types_central_atoms_3body_im
              do j=1,n_3b(i)
                 write(output_epe,'(3a4,3f12.5)') name_of_type(indexx(i,j,1)), &
                      name_of_type(indexx(i,j,2)),name_of_type(indexx(i,j,3)), &
                      ki(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), &
                      theta_0(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)),r3b_im(i)
              enddo
           enddo
           deallocate(n_3b,indexx)

        endif

        write(output_epe,*) '=========================================================='
     endif
  endif

  if(n_centres_of_generation.eq.0.and.pg_interfaced_mode.eq.0) pg_interfaced_mode=1

  deallocate(host_tmp%sr1, &
       host_tmp%k,host_tmp%k1,host_tmp%r1,host_tmp%r0,&
       host_tmp%b,&
       host_tmp%ro,&
       host_tmp%c,&
       host_tmp%d,&
       host_tmp%sr2,&
       stat=epealloc_stat(3))
  if(epealloc_stat(3).ne.0) call error_handler("deallocate host_tmp%sr1 failed")

  E2BIND=zero
  E2BDIS=zero
  E2BDEF=zero
  RECA=1./scaling_factor
  DMAX=scaling_factor/4.
  epeit_count=0
  DEL=scaling_factor/STEP_OF_DISPLACEMENT
  ERROR_FUNCTION_PARAMETER=ERROR_FUNCTION_PARAMETER/4
  PIS=SQRT(PI)
  ET2=ERROR_FUNCTION_PARAMETER**2
  ERFO=-2.*ERROR_FUNCTION_PARAMETER/PIS

  call returnclose_iounit(input_epe)
  print*, 'returnclose input_epe'

end subroutine read_epe_input_new

subroutine convert_epe_input()
  use type_module
  use iounitadmin_module
  use epecom_module
  
  integer(kind=i4_kind) :: input_interm,input_epe1,ind
  character(len=200) :: input_line
  character(len=3) :: buf
  character(len=6) :: form

  input_epe1=openget_iounit(trim(epe_input_dir)//'/epe.input', &
                                 form='formatted', status='old')
  input_interm=openget_iounit(trim(epe_input_dir)//'/epe.inp.interm', &
                                 form='formatted')

  do 
     read(input_epe1,'(a200)',end=100) input_line
     ind=index(input_line,'#')
     if (ind == 1) cycle
     if (ind == 0) ind=len_trim(input_line)+1
     if (ind == 1) cycle
     write(buf,'(i3)') ind-1
     form='(a'//trim(buf)//')'
     write(input_interm,form) trim(input_line)
  enddo

100  call returnclose_iounit(input_epe1)
  call returnclose_iounit(input_interm)

end subroutine convert_epe_input


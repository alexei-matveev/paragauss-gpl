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
subroutine read_epeinput
!the procedure that reads in new EPE input file
# include <def.h>

  use type_module
  use epecom_module
  use ewaldpc_module, only: get_epe_references,get_qm_references,no_pg,use_epe_qm_references,no_relax
  use inp_out_module, only: upcase,check_string !from MM
  use qm_epe_interface_module
  use mol_module
  use str_module
  use epe_pg_module,only:an_epe_type,operations_make_reg_reference
  use iounitadmin_module
  use epepar_module, only: treat_epepar_namelist, &
       ec_name_of_type => name_of_type, get_slsp, &
         ec_max_type_ions=>max_type_ions
  use filename_module, only: env => filename_env, filename_namelengthmax

  implicit none

  integer(i4_kind) :: ion_types(ndt)
  logical :: type_number(ndt,ndcell)
  logical :: new_reg

  integer(i4_kind) :: cell_unit
  integer(i4_kind) :: status
  logical :: lcellvec
  logical :: read_at

  character(len=filename_namelengthmax) :: libdir

  integer(i4_kind) :: m
  integer(i4_kind) ::  STEP_OF_DISPLACEMENT
!------------------------------------------------------------------------

  lpotcalc=.false.

  scaling_factor=1.889726877774_r8_kind*auangs
  radius_long_interact= 100.0_r8_kind
  ERROR_FUNCTION_PARAMETER=2.917719_r8_kind !df_ERROR_FUNCTION_PARAMETER
  n_centres_of_generation=0
  n_vacancies=0
  n_impurities=0
  n_types_central_atoms_3body_im=0
  STEP_OF_DISPLACEMENT=200
  pc_aswritemode=.true.

  input_epe=openget_iounit(trim(epe_input_dir)//'/epe.input', &
                                 form='formatted', status='old')

  output_epe=openget_iounit(trim(epe_input_dir)//'/epe.out', &
                                 form='formatted', status='unknown')

  inquire (file=trim(epe_input_dir)//'/cellvec',exist=lcellvec)
  if(lcellvec) cell_unit=openget_iounit(trim(epe_input_dir)//'/cellvec',  &
       form='formatted', status='unknown')

  if(get_epe_references) then
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) '  Calculation of EPE environment reference files  '
     write(output_epe,*) '        epe.r, epe.pcr(cs), reg_reference         '
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) ''
  else if(get_qm_references) then
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) '    Calculation of QM cluster reference files     '
     write(output_epe,*) '          epe_reference, pgepe_reference          '
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) ''
  else if(use_epe_qm_references) then
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) '          Relaxation of EPE environment           '
     write(output_epe,*) '                around QM cluster                 '
     write(output_epe,*) '--------------------------------------------------'
     write(output_epe,*) ''
  end if
  write(output_epe,*) '=================================================================='
  write(output_epe,*) 'PERIODIC SYSTEM'

  if(lcellvec) call read_atom_names()

  call read_cell_vectors()

  call read_coordinates()

  if(lcellvec) call returnclose_iounit(cell_unit)

  call read_force_field()

  call read_epe_param()

  if(lpotcalc) then
  end if

  if(get_epe_references) then
     operations_make_reg_reference=.true.
     basic_action=1
     qm_interfaced_mode=.false.
     read_configuration=.true.
     write_configuration=.true.
     make_epe_reference=.false.
     use_epe_reference=.false.   
     use_pgdata=.false.

     ml_displacements=.false.
     ml_cluster_simulated=.false.
     explicit_coupling=.false.
  else if(get_qm_references) then
     operations_make_reg_reference=.false.
     basic_action=0
     qm_interfaced_mode=.true.
     read_configuration=.true.
     write_configuration=.true.
     make_epe_reference=.true.
     use_epe_reference=.false.   
     use_pgdata=.true.
     explicit_coupling=.true.

     ml_displacements=.false.
     ml_cluster_simulated=.false.
  else if(use_epe_qm_references.and.no_pg) then
     operations_make_reg_reference=.false.
     basic_action=0
     qm_interfaced_mode=.true.
     read_configuration=.true.
     write_configuration=.false.
     make_epe_reference=.false.
     use_epe_reference=.true.   
     use_pgdata=.false.
     explicit_coupling=.true.
     n_iterations=1

     ml_displacements=.false.
     ml_cluster_simulated=.false.
  else if(use_epe_qm_references) then
     operations_make_reg_reference=.false.
     basic_action=1
     if(no_relax) basic_action=0
     qm_interfaced_mode=.true.
     read_configuration=.true.
     write_configuration=.true.
     make_epe_reference=.false.
     use_epe_reference=.true.   
     use_pgdata=.true.
     explicit_coupling=.true.

     ml_displacements=.false.
     ml_cluster_simulated=.false.
  end if
  reg_2a_treated=ml_displacements.and..not.ml_cluster_simulated

  if(operations_make_reg_reference .or. use_pgdata.or.use_epe_reference) then
     call define_ref_charges()
  endif

  if(explicit_coupling) call read_ec_ff()

  if(get_epe_references .or. get_qm_references .or. use_epe_qm_references) then
     call read_optimization_param()
  end if
  if(get_qm_references .and. no_pg) n_iterations=0

  call read_xyz_output()

  if(ml_cluster_simulated.and.n_ml_cluster_simulators.ne.0) then
  end if
  allocate(ml_dprime_factors(max_type_ions),ml_fac(max_type_ions),stat=status)
  ASSERT(status.eq.0)
  
  if(use_pgdata.or.use_epe_reference) then
  end if

  if(n_impurities==0) then
     if(.not.explicit_coupling) then 
        ALLOCATE(host%sr1(max_type_ions,max_type_ions,max_type_ions+1), & 
             host%b(max_type_ions,max_type_ions,max_type_ions+1),&
             host%k(max_type_ions,max_type_ions,max_type_ions+1),&
             host%k1(max_type_ions,max_type_ions,max_type_ions+1),&
             host%r1(max_type_ions,max_type_ions,max_type_ions+1),&
             host%r0(max_type_ions,max_type_ions,max_type_ions+1),&
             host%ro(max_type_ions,max_type_ions,max_type_ions+1),&
             host%c(max_type_ions,max_type_ions,max_type_ions+1),&
             host%d(max_type_ions,max_type_ions,max_type_ions+1),&
             host%sr2(max_type_ions,max_type_ions,max_type_ions+1),& 
             stat=status)
        ASSERT(status == 0) 

        host%k=0.0_r8_kind 
        host%k1=0.0_r8_kind 
        host%r1=0.0_r8_kind 
        host%r0=0.0_r8_kind 

        host%b=0.0_r8_kind 
        host%ro=0.0_r8_kind 
        host%ro(:,:,1)=1.0_r8_kind  
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

        host%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= &       
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
             host_tmp%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)
     end if
  end if

  if(n_centres_of_generation.eq.0.and.pg_interfaced_mode.eq.0) pg_interfaced_mode=1

  deallocate(host_tmp%sr1, &
       host_tmp%k,host_tmp%k1,host_tmp%r1,host_tmp%r0,&
       host_tmp%b,&
       host_tmp%ro,&
       host_tmp%c,&
       host_tmp%d,&
       host_tmp%sr2,&
       stat=status)
  ASSERT(status == 0)

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

contains

!##############################################################################

!##############################################################################
  subroutine read_atom_names()

    character(len=80) :: buffer
    character(len=10) :: buffer1
    integer(i4_kind) :: i,j
    logical, parameter :: df_new_reg=.true.

    new_reg=df_new_reg

    read_block=find_block(input_epe,"&ATOM_NAMES",epeinp)
    if(.not.read_block) then
       read_at=.false.
    else
       read_at=.true.
       i=0; j=0
       do
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          i=i+1
          if(i==1) then
             !reading in number of atoms in the unit cell
             call upcase(buffer)
             if(check_string(buffer,"N_ATOMS")) then
                read(buffer,*,err=100)  buffer1, n_ions_cell
                if(check_string(buffer,"OLD")) then
                   new_reg=.false.
                else if(check_string(buffer,"NEW")) then
                   new_reg=.true.
                end if
             else
                call error_handler("Read_epeinput: The first line in &ATOM_NAMES"// &
                     " block has to be number of atoms (N_ATOMS  NNN)")
             end if
          else if(i==2) then
             call upcase(buffer)
             if(check_string(buffer,"N_ATOM_TYPES")) then
                read(buffer,*,err=100)  buffer1, max_type_ions
             else
                call error_handler("Read_epeinput: The second line in &ATOM_NAMES"// &
                     " block has to be number of atomic types (N_ATOM_TYPES  NNN)")
             end if
          else
             j=j+1
             read(buffer,*,err=100)  name_of_type(j)
             ion_types(j)=j
          end if
       end do
       if(j /= max_type_ions) call error_handler("Read_epeinput: Check number of atomic types!!!")
    end if

    return

100 call wrong_line(epeinp)

  end subroutine read_atom_names
!====================================================================================

!====================================================================================
  subroutine read_cell_vectors()

    character(len=80) :: buffer
    integer(i4_kind) :: i,j

    if(read_at) goto 10

    read_block=find_block(input_epe,"&CELL_VECTORS",epeinp)
    if(.not.read_block) then
       call error_handler("Read_epeinput: No cell vectors found")
    else
       j=0
       do 
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          j=j+1
          if(j > 3) exit
          read(buffer,*,err=100) (vectors_trans(i,j),i=1,3)
       end do
       if(j < 3) call error_handler("Read_epeinput: Number of cell vectors < 3")
    end if

10  if(lcellvec) then
       do j=1,3
          read(cell_unit,*,err=200)(VECTORS_TRANS(I,J),i=1,3)
       enddo
    end if

    write(output_epe,*) '--------------------------------------------------------------------'
    write(output_epe,*) 'Cell Vectors'
    if(lcellvec) write(output_epe,*) 'Vectors are taken from cellvec file'
    do j=1,3
       write(output_epe,'(3(3x,F15.9))')(VECTORS_TRANS(I,J),i=1,3)
    enddo
    write(output_epe,*) '--------------------------------------------------------------------'
    VECTORS_TRANS=VECTORS_TRANS*scaling_factor

    return

100 call wrong_line(epeinp)
200 call error_handler("Read_epeinput: Wrong cellvec file")

  end subroutine read_cell_vectors
!====================================================================================

!====================================================================================
  subroutine read_coordinates()

    character(len=80) :: buffer
    character(len=10) :: buffer1
    real(r8_kind) :: r_sh(3),ann
    integer(i4_kind) :: i,j,l,maxt
    character(len=5) :: message

    if(read_at) goto 20

    read_block=find_block(input_epe,"&COORDINATES",epeinp)
    if(.not.read_block) then
       call error_handler("Read_epeinput: No cell coordinates found")
    else
       i=0; j=0
       do
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          i=i+1
          if(i==1) then
             !reading in number of atoms in the unit cell
             call upcase(buffer)
             if(check_string(buffer,"N_ATOMS")) then
                read(buffer,*,err=100)  buffer1, n_ions_cell
                if(check_string(buffer,"OLD")) then
                   new_reg=.false.
                else if(check_string(buffer,"NEW")) then
                   new_reg=.true.
                end if
             else
                call error_handler("Read_epeinput: The first line in &COORINATES"// &
                     " block has to be number of atoms (N_ATOMS  NNN)")
             end if
          else
             !reading atomic coordinates
             j=j+1
             fixed(j)=.false.
             read(buffer,*,err=100,end=200)  name_of_ions(j),r_ion_in_cell(j,:),fixed(j)
             goto 10
200          read(buffer,*,err=100)  name_of_ions(j),r_ion_in_cell(j,:)
10           r_core_in_cell(j,:)=r_ion_in_cell(j,:)
             r_shell_in_cell(j,:)=r_ion_in_cell(j,:)
          end if
       end do
       if(j /= n_ions_cell) call error_handler("Read_epeinput: Check number of atoms!!!")

       max_type_ions=0
       do i=1,n_ions_cell
          call get_type(i)
       end do
    end if

20  if(lcellvec) then
       maxt=0
       type_number=.false.
       R_ION_IN_CELL=zero; R_SHELL_IN_CELL=zero; R_CORE_IN_CELL=zero
       do I=1,N_IONS_CELL
          read(cell_unit,'(f5.2,3f15.7)',err=300) an(i),r_core_in_cell(i,:)

          l=int(an(i))
          if(read_at) then
             maxt=max(l,maxt)
             if(maxt > max_type_ions) call error_handler( &
                  "Number of atomic types in CELLVEC > defined in EPE.INPUT")
             name_of_ions(i)=name_of_type(l)
             TYPE_OF_ION(i)=l
             fixed(i)=.false.
          else
             if(TYPE_OF_ION(i) /= l) then
                write(message,'(i5)') i
                call error_handler( &
                     "Non coincided types of atom "//message//" in EPE.INPUT and CELLVEC files")
             end if
          end if

          type_number(TYPE_OF_ION(i),i)=.true.
       enddo
       if(read_at) then
          if(maxt < max_type_ions) max_type_ions=maxt
       end if
       l1: do
          read(cell_unit,'(f5.2,3f15.7)',end=22,err=300)ann,r_sh
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
22     core_shell=.false.
       do i=1,n_ions_cell
          if(dot_product(R_SHELL_IN_CELL(i,:),R_SHELL_IN_CELL(i,:)) /= zero) then
             core_shell=.true.
             if(new_reg) then 
                r_ion_in_cell(i,:)=r_core_in_cell(i,:)
             else
                r_ion_in_cell(i,:)=(r_core_in_cell(i,:)+r_shell_in_cell(i,:))/2.0_r8_kind
             end if
          else
             r_ion_in_cell(i,:)=r_core_in_cell(i,:)
             r_shell_in_cell(i,:)=r_core_in_cell(i,:)
          endif
       enddo
    end if

    if(lcellvec) write(output_epe,*) 'Read coordinates of cell from cellvec file'
    write(output_epe,*) 'Number  Type   Name                  Coordinates (angstrom)'
    do i=1,n_ions_cell
       if(fixed(i)) then
          write(output_epe,'(1x,I3,6x,I2,4X,2A3,5x,3F14.9,1x,l1)') &
               i,type_of_ion(i),NAME_OF_IONS(i),' c ',&
               r_core_in_cell(i,:),fixed(i)
          write(output_epe,'(1x,I3,6x,I2,4X,2A3,5x,3F14.9,1x,l1)') &
               i,type_of_ion(i),NAME_OF_IONS(i),' s ',&
               r_shell_in_cell(i,:),fixed(i)
       else
          write(output_epe,'(1x,I3,6x,I2,4X,2A3,5x,3F14.9)') &
               i,type_of_ion(i),NAME_OF_IONS(i),' c ',&
               r_core_in_cell(i,:)
          write(output_epe,'(1x,I3,6x,I2,4X,2A3,5x,3F14.9)') &
               i,type_of_ion(i),NAME_OF_IONS(i),' s ',&
               r_shell_in_cell(i,:)
       end if
    end do
    write(output_epe,*) '--------------------------------------------------------------------'

    return

100 call wrong_line(epeinp)
300 call error_handler("Read_epeinput: Wrong cellvec file")

  end subroutine read_coordinates
!====================================================================================

!====================================================================================
  subroutine get_type(at_num)

    integer(i4_kind) :: at_num
    integer(i4_kind) :: ig
    logical :: previous_type

    if(max_type_ions == 0) then
       type_of_ion(at_num)=1
       max_type_ions=1
    else
       previous_type=.false.
       do ig=1,at_num-1
          if(trim(name_of_ions(at_num))==trim(name_of_ions(ig))) then
             previous_type=.true.
             type_of_ion(at_num)=type_of_ion(ig)
             exit
          end if
       end do
       if(.not.previous_type) then
          max_type_ions=max_type_ions+1
          type_of_ion(at_num)=max_type_ions
       end if
    end if
    ion_types(TYPE_OF_ION(at_num))=type_of_ion(at_num)
    name_of_type(TYPE_OF_ION(at_num))=NAME_OF_IONS(at_num)

  end subroutine get_type
!====================================================================================

!====================================================================================
  subroutine read_force_field()

    character(len=80) :: buffer
    integer(i4_kind) :: i

    if(.not.env("TTFSLIBS",libdir)) call error_handler( &
         "Read_epeinput: environment variable TTFSLIBS is not set")

    read_block=find_block(input_epe,"&FORCE_FIELD",epeinp)
    if(.not.read_block) then
       call error_handler("Read_epeinput: No force field found")
    else
       i=0
       do
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          i=i+1
          if(i > 1) cycle
          read(buffer,*,err=100) ff_name
       end do
    end if

    ff_name=adjustl(ff_name)
    inquire (file=trim(libdir)//"/"//trim(ff_name),exist=l_ff)
    if(l_ff) then 
       ff_unit=openget_iounit(trim(libdir)//"/"//trim(ff_name),  &
            form='formatted', status='unknown',action='read')
!    print*,'ff_lib ', trim(libdir)//"/"//trim(ff_name)
    else
       call error_handler("Read_epeinput: file "//trim(libdir)//"/"//trim(ff_name)//" not found")
    end if
    
    call FF_read_atom_param()

    call FF_read_2body_param()

    call FF_read_3body_param()

    call returnclose_iounit(ff_unit)

    return

100 call wrong_line(epeinp)

  end subroutine read_force_field
!====================================================================================

!====================================================================================
  subroutine FF_read_atom_param()

    character(len=80) :: buffer
    character(len=6) :: atname
    real(r8_kind) :: q_core_epe,q_shell_epe,pk_epe,q_coupl_qmepe,epe_type
    integer(i4_kind) :: i,j
    logical :: found_type(100)

    q_coupl_qmepe = 0.0_r8_kind

    read_block=find_block(ff_unit,"&ATOM_PROPERTIES",fflib)
    if(.not.read_block) then
       call error_handler("Read_fflib: No atomic properties found")
    else
       found_type=.false.
       i=0
       do
          call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          i=i+1
          read(buffer,*,err=100) atname, q_core_epe, q_shell_epe, pk_epe, epe_type

          do j=1,max_type_ions
             if(trim(atname)==trim(name_of_type(j))) then
                pk(j)=pk_epe
                q_shell(j)=q_shell_epe*qau_qepe
                q_nuclear(j)=q_core_epe*qau_qepe
                q_ion(j)=q_shell(j)+q_nuclear(j)
                q_epecl_coupl(j)=q_coupl_qmepe
                an_epe_type(j)=epe_type
                found_type(j)=.true.
                exit
             end if
          end do
       end do
    end if

    do i=1,max_type_ions
       if(.not.found_type(i)) then
          print*,trim(name_of_type(i))//" atom is not defined in Force Field used"
          call error_handler("Read_fflib: Perhaps unit cell used contains wrong atom name")
       end if
    end do

    write(output_epe,*) 'Force Field data taken from '//trim(ff_name)
    write(output_epe,*) '-------------------------------------------'
    write(output_epe,*) '          Shell-model parameters        '
    write(output_epe,*) '-------------------------------------------'
    write(output_epe,*) '        q_core      q_shell      pki    '
    do i=1,max_type_ions
       write(output_epe,'(a4,3f12.5)') name_of_type(i),&
            q_nuclear(i)/qau_qepe,q_shell(i)/qau_qepe,pk(i)
    enddo

    do i=1,n_ions_cell
       do j=1,max_type_ions
          if(name_of_ions(i)==name_of_type(j)) then
             q_z(i)=q_ion(j)/qau_qepe
             exit
          endif
       enddo
    enddo

    return

100 call wrong_line(fflib)

  end subroutine FF_read_atom_param
!====================================================================================

!====================================================================================
  subroutine define_ref_charges

    integer(i4_kind) :: i

    do i=1,max_type_ions
       out_charges(i)=q_ion(i)/qau_qepe
    enddo

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

  end subroutine define_ref_charges
!====================================================================================

!====================================================================================
  subroutine FF_read_2body_param()

    character(len=80) :: buffer
    character(len=6) :: atm_name(2),atm_common
    real(r8_kind) :: b0,ro0,c0,cutoff
    integer(i4_kind) :: common_atom
    integer(i4_kind) :: i,j,k,jj

    allocate(host_tmp%sr1(30,30,31), &
         host_tmp%k(30,30,31),&
         host_tmp%k1(30,30,31),&
         host_tmp%r1(30,30,31),&
         host_tmp%r0(30,30,31),&
         host_tmp%b(30,30,31),&
         host_tmp%ro(30,30,31),&
         host_tmp%c(30,30,31),&
         host_tmp%d(30,30,31),&
         host_tmp%sr2(30,30,31),&
         stat=status)
    ASSERT(status.eq.0)

    host_tmp%k= 0.0_r8_kind
    host_tmp%k1= 0.0_r8_kind
    host_tmp%r1= 0.0_r8_kind
    host_tmp%r0= 0.0_r8_kind
    host_tmp%b= 0.0_r8_kind
    host_tmp%ro= 0.0_r8_kind
    host_tmp%ro(:,:,1) = 1.0_r8_kind
    host_tmp%c= 0.0_r8_kind
    host_tmp%d=0.0_r8_kind
    host_tmp%sr1= 4.5_r8_kind
    host_tmp%sr2= 4.5_r8_kind

    m=max_type_ions*(max_type_ions-1)/2+max_type_ions

    read_block=find_block(ff_unit,"&TWO-BODY",fflib)
    if(.not.read_block) then
       call error_handler("Read_fflib: No atomic two body parameters found")
    else
       do
          call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          read(buffer,*,err=100) atm_name,b0,ro0,c0,cutoff,atm_common
          do j=1,max_type_ions
             if(trim(atm_name(1)) == trim(name_of_type(j))) exit
          end do
          do k=1,max_type_ions
             if(trim(atm_name(2)) == trim(name_of_type(k))) exit
          end do
          if(j > max_type_ions .or. k > max_type_ions) cycle

          if(trim(atm_common)=="NO") then
             common_atom=0
          else
             do i=1,max_type_ions
                if(trim(atm_common) == trim(name_of_type(i))) exit
             end do
             if(i > max_type_ions) then 
                common_atom=0
             else
                common_atom=i
             end if
          end if

          jj=common_atom+1
          
          host_tmp%k(J,K,jj)=0.0
          host_tmp%k(K,J,jj)=0.0
          host_tmp%k1(J,K,jj)=0.0
          host_tmp%k1(K,J,jj)=0.0
          host_tmp%r1(J,K,jj)=0.0
          host_tmp%r1(K,J,jj)=0.0
          host_tmp%r0(J,K,jj)=0.0
          host_tmp%r0(K,J,jj)=0.0

          host_tmp%B(J,K,jj)=b0
          host_tmp%B(K,J,jj)=b0 
          host_tmp%Ro(J,K,jj)=ro0 
          host_tmp%RO(K,J,jj)=ro0 
          host_tmp%C(J,K,jj)=c0  
          host_tmp%C(K,J,jj)=c0 
          host_tmp%sr1(K,J,jj)=cutoff 
          host_tmp%sr1(j,k,jj)=cutoff 
          host_tmp%sr2(K,J,jj)=cutoff 
          host_tmp%sr2(j,k,jj)=cutoff 
       end do
    end if

    write(output_epe,*) '-------------------------------------------------------------------------------'
    write(output_epe,*) '              The Potential Parameters of the Lattice'
    write(output_epe,*) '-------------------------------------------------------------------------------'
    write(output_epe,*) '                        Two-body parameters'
    write(output_epe,*) '-------------------------------------------------------------------------------'
    write(output_epe,*) '              b           ro          c           d        cutoff   common atom'
    do i=1,max_type_ions
       do j=i,max_type_ions
          k=0
          do 
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
          enddo
       end do
    enddo
    write(output_epe,*) '-------------------------------------------------------------------------------'

    return

100 call wrong_line(fflib)

  end subroutine FF_read_2body_param
!====================================================================================

!====================================================================================
  subroutine FF_read_3body_param

    character(len=80) :: buffer
    integer(i4_kind) :: N_central_at
    character(len=12) :: aaa
    character(len=6) :: central_nm(20)
    integer(kind=i4_kind), allocatable :: n_3b(:),indexx(:,:,:)
    integer(kind=i4_kind) :: t_3b(20)
    character(len=3) :: atm_nm(3)
    real(kind=r8_kind) :: k_i,theta_i,r_3b

    integer(i4_kind) :: i,j,k,l,i1,i2,i3

    n_types_central_atoms_3body=0
!    print*,'FF_read_3body_param'

    read_block=find_block(ff_unit,"&THREE-BODY",fflib)
    if(.not.read_block) then
       return
    else
       i=0
       do
          call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          i=i+1
          if(i==1) then
             read(buffer,*,err=100) aaa, N_central_at

          else if(i==2) then
             read(buffer,*,err=100) central_nm(1:N_central_at)
             l=0
             do j=1,N_central_at
                do k=1,max_type_ions
                   if(trim(central_nm(j)) == trim(name_of_type(k))) then 
                      n_types_central_atoms_3body=n_types_central_atoms_3body+1
                      l=l+1
                      t_3b(l)=k
                      exit
                   end if
                end do
             end do

          else
             if (n_types_central_atoms_3body == 0) cycle
             if(i==3) then
                allocate(ki(max_type_ions,max_type_ions,max_type_ions), &
                     theta_0(max_type_ions,max_type_ions,max_type_ions),stat=status)
                ASSERT(status.eq.0)
                allocate(types(n_types_central_atoms_3body,5), r3b(n_types_central_atoms_3body), &
                     n_3b(n_types_central_atoms_3body),indexx(n_types_central_atoms_3body,6,3), &
                     stat=status)
                ASSERT(status.eq.0)
                ki=0.0_r8_kind
                n_3b=0
                types=0
             end if

             read(buffer,*,err=100) atm_nm,k_i,theta_i,r_3b
             i1=0; i2=0; i3=0
             do k=1,max_type_ions
                if(atm_nm(1) == name_of_type(k)) i1=k
                if(atm_nm(2) == name_of_type(k)) i2=k
                if(atm_nm(3) == name_of_type(k)) i3=k
             enddo
             if(i1 ==0 .or. i2 == 0 .or. i3 == 0) cycle

             do j=1,n_types_central_atoms_3body
                if(i2 == t_3b(j)) then 
                   n_3b(j)=n_3b(j)+1
                   exit
                end if
             end do
             l=n_3b(j)
             ki(i1,i2,i3)=k_i
             ki(i3,i2,i1)=k_i
             theta_0(i1,i2,i3)=theta_i
             theta_0(i3,i2,i1)=theta_i
             r3b(j)=r_3b
             indexx(j,l,1)=i1
             indexx(j,l,2)=i2
             indexx(j,l,3)=i3

             if(l==1) then
                types(j,1)=i2
             endif
             do k=2,5
                if(types(j,k)==i1) exit
                if(types(j,k)==0) then
                   types(j,k)=i1
                   exit
                endif
             enddo
             do k=2,5
                if(types(j,k)==i3) exit
                if(types(j,k)==0) then
                   types(j,k)=i3
                   exit
                endif
             enddo
          end if
       end do
    end if

    if (n_types_central_atoms_3body /= 0) then
       write(output_epe,*) '              Three-body parameters(R)'
       write(output_epe,*) '-------------------------------------------------'
       write(output_epe,*) '                  ki        theta       r3b'
!       print*,name_of_type(:)
       do i=1,n_types_central_atoms_3body
!          print*,i,n_3b(i),' type_central_atoms_3body'
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

    return

100 call wrong_line(fflib)

  end subroutine FF_read_3body_param
!====================================================================================

!====================================================================================
  subroutine read_ec_ff()

    ff_name=adjustl(ff_name)
    inquire (file=trim(libdir)//"/"//trim(ff_name),exist=l_ff)
    if(l_ff) then 
       ff_unit=openget_iounit(trim(libdir)//"/"//trim(ff_name),  &
            form='formatted', status='unknown',action='read')
    else
       call error_handler("Read_epeinput:read_ec_ff: file "//trim(libdir)//"/"//trim(ff_name)//" not found")
    end if

    call get_ec_type()

    call read_ec_2body_param()

    call read_ec_3body_param()

    call returnclose_iounit(ff_unit)

  end subroutine read_ec_ff
!====================================================================================

!====================================================================================
  subroutine get_ec_type()

    logical :: f_exist,previous_type
    integer(i4_kind) :: gx_unit,status,i,j,nn
    real(r8_kind),allocatable :: ec_name(:)
    real(r8_kind) :: an,x,y,z
    integer(i4_kind) :: num(9)
    character(len=256) :: buffer
    real(r8_kind), parameter :: small=1.0e-3_r8_kind
    
    inquire(file=trim(epe_input_dir)//'/gxfile',exist=f_exist)
    if(.not.f_exist) then
       call error_handler("Read_epeinput:read_ec_ff:get_ec_type: gxfile not existed")
    end if
    gx_unit=openget_iounit(trim(epe_input_dir)//'/gxfile')

    allocate(ec_name(500),stat=status)
    ASSERT(status==0)

    do i=1,max_type_ions
       ec_name(i)=an_epe_type(i)
    end do

    i=max_type_ions
    do
       read(gx_unit,'(a256)',err=100) buffer
       read(buffer,*,err=100) an
       if(an <= 0.5_r8_kind) exit

       read(buffer,*,err=100) an,x,y,z,num

       do j=1,max_type_ions
          if(abs(an-ec_name(j)) <= small) call error_handler( &
               "Read_epeinput:read_ec_ff:get_ec_type: Coincided atom names in GXFILE and EPE.PCS files")
       end do
       if(num(9) == -1) then
          i=i+1
          ec_name(i)=an
       end if
    end do

    nn=i

    call returnclose_iounit(gx_unit)

    ec_max_type_ions=0
    do i=1,nn
       previous_type=.false.
       if(ec_max_type_ions == 0) then
          ec_max_type_ions=1
       else
          do j=1,i-1
             if(ec_name(i)==ec_name(j)) then
                previous_type=.true.
                exit
             end if
          end do
          if(.not.previous_type) then
             ec_max_type_ions=ec_max_type_ions+1
          end if
       end if
       if(.not.previous_type) ec_name_of_type(ec_max_type_ions)=ec_name(i)
    end do
    ASSERT(ec_max_type_ions>max_type_ions)

    deallocate(ec_name,stat=status)
    ASSERT(status==0)

    return

100 call error_handler("Read_epeinput:read_ec_ff:get_ec_type: error of reading in gxfile")

  end subroutine get_ec_type
!====================================================================================

!====================================================================================
  subroutine read_ec_2body_param()

    character(len=80) :: buffer
    character(len=6) :: atm_name(2)
    real(r8_kind) :: b0,r0,c0,cutoff,r_at_name(2)
    real(r8_kind), parameter :: small=1.0e-3_r8_kind
    integer(i4_kind) :: status,i,j,n_2b

    allocate(host%sr1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1), & 
         host%b(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%ro(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%c(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& 
         host%d(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),& 
         host%k(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%r0(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%k1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%r1(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         host%sr2(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions+1),&
         stat=status)
    ASSERT(status==0)

    host%k=0.0_r8_kind 
    host%k1=0.0_r8_kind 
    host%r1=0.0_r8_kind 
    host%r0=0.0_r8_kind 

    host%b=0.0_r8_kind 
    host%ro=0.0_r8_kind; host%ro(:,:,1)=1.0_r8_kind
    host%c=0.0_r8_kind; host%d=0.0_r8_kind
    host%sr1=10.0_r8_kind
    host%sr2=10.0_r8_kind

    host%sr1(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)= & 
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
         host_tmp%sr2(1:max_type_ions,1:max_type_ions,1:max_type_ions+1)

    read_block=find_block(ff_unit,"&TWO-BODY",fflib)
    if(.not.read_block) then
       call error_handler("Read_ec_2body_param: No atomic two body parameters found")
    else
       n_2b=0
       do
          call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit

          read(buffer,*,err=100) atm_name,b0,r0,c0,cutoff
          read(atm_name(1),*,iostat=status) r_at_name(1)
          if(status /= 0) cycle
          read(atm_name(2),*,iostat=status) r_at_name(2)
          if(status /= 0) cycle

          do i=1,ec_max_type_ions
             if(abs(r_at_name(1)-ec_name_of_type(i)) <= small) exit
          end do
          do j=1,ec_max_type_ions
             if(abs(r_at_name(2)-ec_name_of_type(j)) <= small) exit
          end do
          if(i > ec_max_type_ions .or. j > ec_max_type_ions) cycle
          if(i<=max_type_ions .and. j<=max_type_ions) cycle
          n_2b=n_2b+1

          host%b(i,j,:)=b0
          host%b(j,i,:)=b0
          host%ro(i,j,:)=r0
          host%ro(j,i,:)=r0
          host%c(j,i,:)=c0
          host%c(i,j,:)=c0
          host%d(j,i,:)=0.0_r8_kind
          host%d(i,j,:)=0.0_r8_kind
          host%sr1(j,i,:)=cutoff
          host%sr1(i,j,:)=cutoff
          host%sr2(j,i,:)=cutoff
          host%sr2(i,j,:)=cutoff
       end do
    end if

    if(n_2b==0) then
       call error_handler( &
            "Read_ec_2body_param: No 2body interaction between QM cluster and EPE")
    else
       write(output_epe,*) '========================================================================='
       write(output_epe,*) '              The Potential Parameters of the Explicit_coupling'
       write(output_epe,*) '-------------------------------------------------------------------------' 
       write(output_epe,*) '                          Two-body parameters'
       write(output_epe,*) '-------------------------------------------------------------------------' 
       write(output_epe,*) '                  b           ro          c           d         cutoff' 
!       print*,'  ec_name_of_type'
       do i=1,ec_max_type_ions
!        print*, i, ec_name_of_type(i)
       enddo
       do i=1,ec_max_type_ions
          do j=i,ec_max_type_ions
             if(i<=max_type_ions .and. j<=max_type_ions) cycle
             if(host%b(i,j,1) == 0.0_r8_kind) cycle
             write(output_epe,'(2f6.2,5f12.5)') ec_name_of_type(i),ec_name_of_type(j) &
                  ,host%b(i,j,1),host%ro(i,j,1),host%c(i,j,1),host%d(i,j,1),host%sr1(i,j,1)
          enddo
       enddo
       write(output_epe,*) '========================================================================='
    end if

    return

100 call wrong_line(fflib)

  end subroutine read_ec_2body_param

    subroutine read_ec_3body_param()

      character(len=80) :: buffer
      integer(i4_kind) :: N_central_at
      character(len=12) :: aaa
      character(len=6) :: central_nm(20)
      real(kind=r8_kind) :: r_c_nm
      integer(kind=i4_kind), allocatable :: n_3b(:),indexx(:,:,:)
      integer(kind=i4_kind) :: t_3b(20)
      character(len=6) :: atm_nm(3)
      real(kind=r8_kind) :: at_name(3)
      real(kind=r8_kind) :: k_i,theta_i,r_3b

      real(r8_kind), parameter :: small=1.0e-3_r8_kind
      integer(i4_kind) :: i,j,k,l,i1,i2,i3,status

      n_ec_types_central_atoms_3body=0
      read_block=find_block(ff_unit,"&THREE-BODY",fflib)
      if(.not.read_block) then
         return
      else
         i=0
         cycle1: do
            call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
            if(empty_line) cycle cycle1
            if(comment_line) cycle cycle1
            if(end_block) exit cycle1
            i=i+1
            if(i==1) then
               read(buffer,*,err=100) aaa, N_central_at

            else if(i==2) then
               read(buffer,*,err=100) central_nm(1:N_central_at)
               l=0
               do j=1,N_central_at
                  read(central_nm(j),*,iostat=status) r_c_nm
                  if(status /= 0) cycle
                  do k=1,ec_max_type_ions
                     if(abs(r_c_nm-ec_name_of_type(k)) <= small) then
                        n_ec_types_central_atoms_3body=n_ec_types_central_atoms_3body+1
                        l=l+1
                        t_3b(l)=k
                        exit
                     end if
                  end do
               end do
            else
               if (n_ec_types_central_atoms_3body == 0) cycle
               if(i==3) then
                  allocate(ec%ki(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions), &
                       ec%theta_0(ec_max_type_ions,ec_max_type_ions,ec_max_type_ions),stat=status)
                  ASSERT(status.eq.0)
                  allocate(ec%types(n_ec_types_central_atoms_3body,5), ec%r3b(n_ec_types_central_atoms_3body), &
                       n_3b(n_ec_types_central_atoms_3body),indexx(n_ec_types_central_atoms_3body,20,3), &
                       stat=status)
                  ASSERT(status.eq.0)
                  ec%types=0
                  ec%ki=0.0_r8_kind
                  n_3b = 0
                  ec%types = 0
               end if

               read(buffer,*,err=100) atm_nm,k_i,theta_i,r_3b
               read(atm_nm(1),*,iostat=status) at_name(1)
               if(status /= 0) cycle cycle1
               read(atm_nm(2),*,iostat=status) at_name(2)
               if(status /= 0) cycle cycle1
               read(atm_nm(3),*,iostat=status) at_name(3)
               if(status /= 0) cycle cycle1
               i1=0; i2=0; i3=0
               do k=1,ec_max_type_ions
                  if(abs(at_name(1) - ec_name_of_type(k)) <= small) i1=k
                  if(abs(at_name(2) - ec_name_of_type(k)) <= small) i2=k
                  if(abs(at_name(3) - ec_name_of_type(k)) <= small) i3=k
               enddo
               if(i1 ==0 .or. i2 == 0 .or. i3 == 0) cycle

               do j=1,n_ec_types_central_atoms_3body
                  if(i2 == t_3b(j)) then
                     n_3b(j)=n_3b(j)+1
                     exit
                  end if
               end do
               l=n_3b(j)
               ec%ki(i1,i2,i3)=k_i
               ec%ki(i3,i2,i1)=k_i
               ec%theta_0(i1,i2,i3)=theta_i
               ec%theta_0(i3,i2,i1)=theta_i
               ec%r3b(j)=r_3b
               indexx(j,l,1)=i1
               indexx(j,l,2)=i2
               indexx(j,l,3)=i3

               if(l==1) then
                  ec%types(j,1)=i2
               endif
               do k=2,5
                  if(ec%types(j,k)==i1) exit
                  if(ec%types(j,k)==0) then
                     ec%types(j,k)=i1
                     exit
                  endif
               enddo
               do k=2,5
                  if(ec%types(j,k)==i3) exit
                  if(ec%types(j,k)==0) then
                     ec%types(j,k)=i3
                     exit
                  endif
               enddo
            end if
         end do cycle1
      end if

      if (n_ec_types_central_atoms_3body /= 0) then
         write(*,*) '              Three-body parameters'
         write(*,*) '------------------------------------------------------'
         write(*,*) '                       ki        theta         r3b'
         do i=1,n_ec_types_central_atoms_3body
            do j=1,n_3b(i)
               write(*,'(3f6.2,3f12.5)') ec_name_of_type(indexx(i,j,1)), &
                    ec_name_of_type(indexx(i,j,2)),ec_name_of_type(indexx(i,j,3)), &
                    ec%ki(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), &
                    ec%theta_0(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), ec%r3b(i)
            enddo
         enddo
         deallocate(n_3b,indexx,stat=status)
         ASSERT(status == 0)
      endif
      write(*,*) '------------------------------------------------------'

      return

100   call wrong_line(fflib)

    end subroutine read_ec_3body_param

!====================================================================================

!====================================================================================
  subroutine read_epe_param

    character(len=80) :: buffer,buffer1
    integer(i4_kind) :: j

    real(r8_kind),parameter :: df_r_2A_sphere=22.0_r8_kind
    real(r8_kind),parameter :: df_r_first_sphere=10.0_r8_kind
    integer(i4_kind),parameter :: df_n_gen_ions=6000
    integer(i4_kind),parameter :: df_n_trans_primitive_cell(3)=(/5,5,5/)
    logical, parameter :: df_option_c3_symm=.false.
    
    real(kind=r8_kind), dimension(3):: xx, yy, zz
    integer(i4_kind) :: status

    radius_2A_sphere=df_r_2A_sphere
    radius_first_sphere=df_r_first_sphere
    option_c3_symm=df_option_c3_symm
    n_gen_ions=df_n_gen_ions
    n_trans_primitive_cell=df_n_trans_primitive_cell

    read_block=find_block(input_epe,"&EPE_PARAMETERS",epeinp)
    if(.not.read_block) then
       goto 10
    else
       j=0
       do 
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          
          call upcase(buffer)
          if(check_string(buffer,"RADIUS_2A_REGION")) then
             read(buffer,*,err=100) buffer1, radius_first_sphere
          else if(check_string(buffer,"RADIUS_2B_REGION")) then
             read(buffer,*,err=100) buffer1, radius_2A_sphere
          else if(check_string(buffer,"N_GEN_IONS")) then
             read(buffer,*,err=100) buffer1, n_gen_ions
          else if(check_string(buffer,"N_PRIMITIVE_CELL_TRANS")) then
             read(buffer,*,err=100) buffer1, n_trans_primitive_cell
          else if(check_string(buffer,"NO_C3_SYMM")) then
             option_c3_symm=.false.
          else if(check_string(buffer,"C3_SYMM")) then
             option_c3_symm=.true.
          end if
       end do
    end if

10  if(radius_first_sphere <= 0.0_r8_kind) &
         call error_handler("Read_epeinput: Radius of 2A region <= 0.0_r8_kind")
    if(radius_2A_sphere <= 0.0_r8_kind) &
         call error_handler("Read_epeinput: Radius of 2B region <= 0.0_r8_kind")
    if(radius_2A_sphere <= radius_first_sphere) &
         call error_handler("Read_epeinput: Radius of 2B region <= Radius of 2A region")

    write(output_epe,*) "Parameters of the EPE environment"
    write(output_epe,*) "----------------------------------------------"
    write(output_epe,*) "Radius of 2A region: ",radius_first_sphere," angstroms"
    write(output_epe,*) "Radius of 2B region: ",radius_2A_sphere," angstroms"
    if(option_c3_symm) then
       write(output_epe,*) "C3 symmetrization"
    end if
    write(output_epe,*) "Number of generated EPE centers: ",n_gen_ions
    write(output_epe,*) "Number of unit cell tanslations: ",n_trans_primitive_cell
    write(output_epe,*) "====================================================================="

    RADIUS_FIRST_SPHERE=RADIUS_FIRST_SPHERE*scaling_factor
    RADIUS_2A_SPHERE=RADIUS_2A_SPHERE*scaling_factor
    RADIUS_LONG_INTERACT=(RADIUS_LONG_INTERACT*scaling_factor)**2

    allocate(r_nuc_ion(n_gen_ions,3), &
         r_sh_ion(n_gen_ions,3), &
         epe(n_gen_ions),q_zl(n_gen_ions),stat=status)
    ASSERT(status.eq.0)

    xx(1)=1.0_r8_kind; xx(2)=0.0_r8_kind;xx(3)=0.0_r8_kind
    yy(1)=0.0_r8_kind; yy(2)=1.0_r8_kind;yy(3)=0.0_r8_kind
    zz(1)=0.0_r8_kind; zz(2)=0.0_r8_kind;zz(3)=1.0_r8_kind

    rot_gto_to_epe(:,1)=xx
    rot_gto_to_epe(:,2)=yy
    rot_gto_to_epe(:,3)=zz

    shft_gto_to_epe=0.0_r8_kind
    shft_gto_to_hds=0.0_r8_kind

    return

100 call wrong_line(epeinp)

  end subroutine read_epe_param
!====================================================================================

!====================================================================================
  subroutine read_optimization_param()

    character(len=80) :: buffer,buffer1

    integer(i4_kind), parameter :: df_n_iterations=100
    integer(i4_kind), parameter :: df_n_hess_update_cycles=10
    real(r8_kind), parameter :: df_weight_hess=1.0_r8_kind
    real(r8_kind), parameter :: df_abs_g=0.0001_r8_kind
    logical, parameter :: df_print_gradients=.false.

    
    n_iterations=df_n_iterations
    n_hess_update_cycles=df_n_hess_update_cycles
    weight_hess=df_weight_hess
    abs_g=df_abs_g
    print_gradients=df_print_gradients

    read_block=find_block(input_epe,"&OPTIMIZATION",epeinp)
    if(.not.read_block) then
       goto 10
    else
       do 
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit
          
          call upcase(buffer)
          if(check_string(buffer,"ITERATIONS")) then
             read(buffer,*,err=100) buffer1, n_iterations
          else if(check_string(buffer,"MAX_GRAD")) then
             read(buffer,*,err=100) buffer1, abs_g
          else if(check_string(buffer,"HESS_UPDATE")) then
             read(buffer,*,err=100) buffer1, n_hess_update_cycles
          else if(check_string(buffer,"HESS_WEIGHT")) then
             read(buffer,*,err=100) buffer1, weight_hess
          else if(check_string(buffer,"NO_PRINT_GRADS")) then
             print_gradients=.false.
          else if(check_string(buffer,"PRINT_GRADS")) then
             print_gradients=.true.
          end if
       end do
    end if


10  if(n_iterations < 0) n_iterations=df_n_iterations
    n_pgepe_iterations=n_iterations
    if(abs_g < 0.0_r8_kind) abs_g=df_abs_g
    if(n_hess_update_cycles < 1) n_hess_update_cycles=df_n_hess_update_cycles
    if(weight_hess <=  0.0_r8_kind) weight_hess=df_weight_hess

    write(output_epe,*) "Optimization parameters of EPE environment"
    write(output_epe,*) "------------------------------------------"
    write(output_epe,*) "Number of iterations: ",n_iterations
    write(output_epe,*) "Maximal gradient to stop: ", abs_g
    write(output_epe,*) "Update Hessian every ",n_hess_update_cycles," cycles"
    write(output_epe,*) "Hessian weight: ",weight_hess
    write(output_epe,*) "====================================================================="

    return

100 call wrong_line(epeinp)

  end subroutine read_optimization_param
!====================================================================================

!====================================================================================
  subroutine read_xyz_output()

    character(len=80) :: buffer

    logical, parameter :: df_region_I_output=.false.

    region_I_output=df_region_I_output

    read_block=find_block(input_epe,"&XYZ_OUTPUT",epeinp)
    if(.not.read_block) then
       goto 10
    else
       do 
          call read_block_line(input_epe,buffer,end_block,comment_line,empty_line,epeinp)
          if(empty_line) cycle
          if(comment_line) cycle
          if(end_block) exit

          call upcase(buffer)
          if(check_string(buffer,"NO_OUTPUT_IIA")) then
             region_I_output=.false.
          else if(check_string(buffer,"OUTPUT_IIA")) then
             region_I_output=.true.
          end if
       end do
    end if

10  if(region_I_output) then
       write(output_epe,*) "Saving IIa region in XYZ format"
       write(output_epe,*) "====================================================================="
    end if

  return

100 call wrong_line(epeinp)

  end subroutine read_xyz_output
!====================================================================================

!====================================================================================

end subroutine read_epeinput
